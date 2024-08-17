from rdkit import Chem
from rdkit.Chem import AllChem
import vina
from meeko import MoleculePreparation, PDBQTWriterLegacy, PDBQTMolecule, RDKitMolCreate

import numpy as np
from glob import glob
import os
import tempfile
from sklearn.decomposition import PCA

from .reward_utils import is_exo

# Define class for docking
class DockLigand():
    """ Class for carrying out docking operations
        Requires guestmol, hostmol and conf
        The required contents of conf depends on the docking method
    """

    def __init__(self, hostmol, conf):
        self.conf = conf
        # Non optional variables assigned as class attributes
        self.hostmol = hostmol

    """ VINA SCORING """
    def score_map_vina(self, mol):
        """ Screens across a set of rotations and translations
        """
        hostmol = self.hostmol
        num_rot = self.conf["vina_num_rotations"]
        num_tra = self.conf["vina_num_translations"]
        hostpdbqtfile = self.conf["host_pdbqt"]

        aligned_coords = self.align_mol(mol)
        
        # Create 2d grid of rotations and translations
        rotations = np.linspace(0,90,num_rot)
        translations = np.linspace(0,4,num_tra)
        scores = []
        optmols = []

        # Initialize and setup vina optimizer
        vinaobj = vina.Vina(verbosity=0)
        vinaobj.set_receptor(hostpdbqtfile)
        vinaobj.compute_vina_maps(center=[0.0,0.0,0.0], box_size=[24.0, 24.0, 24.0])

        # Score original conformer
        score, optmol = self.vina_opt(mol, vinaobj)
        scores.append(score)
        optmols.append(optmol)

        # Conformer 2: flip the guest 90 degrees about the x axis (means it probably wont fit, so will not screen in this direction)
        rotation_matrix = np.array([[np.cos(90),-np.sin(90),0], [np.sin(90), np.cos(90), 0], [0,0,1]])
        rotation_coords = np.matmul(rotation_matrix,aligned_coords.T).T

        rotconf = Chem.Conformer(mol.GetConformer(0))
        for index, atm in enumerate(mol.GetAtoms()):
            rotconf.SetAtomPosition(atm.GetIdx(), rotation_coords[index])
        rotmol = Chem.Mol(mol)
        rotmol.RemoveAllConformers()
        rotmol.AddConformer(rotconf, assignId=True)
        # Score flipped conformer
        score, optmol = self.vina_opt(rotmol, vinaobj)
        scores.append(score)
        optmols.append(optmol)

        # Screen across rotations and translations
        for rotate_by in rotations:
            for translate_by in translations:
                new_coords = self.transform_coords(aligned_coords, rotate_by=rotate_by, translate_by=translate_by)

                # Assign conformer
                conf = Chem.Conformer(mol.GetConformer(0))
                for index, atm in enumerate(mol.GetAtoms()):
                    conf.SetAtomPosition(atm.GetIdx(), new_coords[index])
                tramol = Chem.Mol(mol)
                tramol.RemoveAllConformers()
                tramol.AddConformer(conf, assignId=True)

                # Score conformer
                score, optmol = self.vina_opt(tramol, vinaobj)
                scores.append(score)
                optmols.append(optmol)

        # Combine each of the poses with the host and calculate scores
        onehostmol = Chem.Mol(hostmol)
        onehostmol.RemoveAllConformers()
        onehostmol.AddConformer(hostmol.GetConformer(0), assignId=True)
        
        scores = np.array(scores)
        optmols = [Chem.CombineMols(onehostmol, i) for i in optmols]
        # Sort the conformers by score and assign IDs in ascending order
        score_confids = np.argsort(scores)
        scores = [float(i) for i in np.sort(scores)]
        
        finalcomplex = Chem.Mol(optmols[0])
        finalcomplex.RemoveAllConformers()
        for i in score_confids:
            finalcomplex.AddConformer(optmols[i].GetConformer(0), assignId=True)

        return finalcomplex, scores

    """ COMBINED MMFF94 + vina SCORING METHODS """

    def score_map_comb(self, mol):
        """ Uses the mmff94 force field to optimise just a few poses, instead of using vina score.
        """
        hostmol = self.hostmol
        num_rot = self.conf["vina_num_rotations"]
        num_tra = self.conf["vina_num_translations"]
        hostpdbqtfile = self.conf["host_pdbqt"]

        # Conformer 1: principal axis of guest aligned with axis through the cavity of the CB (z-axis)
        aligned_coords = self.align_mol(mol)
        
        # Conformer 2: flip the guest 90 degrees about the x axis (means it probably wont fit, so will not screen in this direction)
        rotation_matrix = np.array([[np.cos(90),-np.sin(90),0], [np.sin(90), np.cos(90), 0], [0,0,1]])
        rotation_coords = np.matmul(rotation_matrix,aligned_coords.T).T

        rotconf = Chem.Conformer(mol.GetConformer(0))
        for index, atm in enumerate(mol.GetAtoms()):
            rotconf.SetAtomPosition(atm.GetIdx(), rotation_coords[index])
        mol.AddConformer(rotconf, assignId=True)
        
        # Create 2d grid of rotations and translations
        rotations = np.linspace(0,90, num_rot)
        translations = np.linspace(0,4, num_tra)

        # Screen across rotations and translations
        for rotate_by in rotations:
            for translate_by in translations:
                # Add the new conformers as conformers to the mol object
                new_coords = self.transform_coords(aligned_coords, rotate_by=rotate_by, translate_by=translate_by)
                # Get a copy of the original conformer to edit with the transformed coordinates
                conf = Chem.Conformer(mol.GetConformer(0))
                for index, atm in enumerate(mol.GetAtoms()):
                    conf.SetAtomPosition(atm.GetIdx(), new_coords[index])
                mol.AddConformer(conf, assignId=True)
        
        # Initialize and setup vina scorer, set no refine as this is a scoring only job
        vinaobj = vina.Vina(verbosity=0, no_refine=True)
        vinaobj.set_receptor(hostpdbqtfile)
        vinaobj.compute_vina_maps(center=[0.0,0.0,0.0], box_size=[24.0, 24.0, 24.0])
        
        # Combine each of the poses with the host and calculate scores
        complexconfs = Chem.CombineMols(hostmol, mol)

        # Optimise the confomers
        Chem.GetSSSR(complexconfs)
        mmff_res = AllChem.MMFFOptimizeMoleculeConfs(complexconfs, numThreads=0, ignoreInterfragInteractions=False)
        scores = np.array([i[1] for i in mmff_res])

        # Use vina scoring to confirm the best MMFF optimised structure
        vina_scores = []
        for pose in range(mol.GetNumConformers()):
            emptymol = Chem.GetMolFrags(complexconfs, asMols=True)[1]
            emptymol.RemoveAllConformers()
            emptymol.AddConformer(Chem.GetMolFrags(complexconfs, asMols=True)[1].GetConformer(pose), assignId=True)
            vina_score = self.vina_scoring(emptymol, vinaobj)
            vina_scores.append(vina_score)

        scores = scores + vina_scores
        # Sort the conformers by score and assign IDs in ascending order
        score_confids = [int(i) for i in np.argsort(scores)]
        scores = [float(i) for i in np.sort(scores)]
        
        finalcomplex = Chem.Mol(complexconfs)
        finalcomplex.RemoveAllConformers()
        for i in score_confids:
            finalcomplex.AddConformer(complexconfs.GetConformer(i), assignId=True)
        
        return finalcomplex, scores
    
    """VINA DOCKING VIA GENERATIVE METHOD"""

    def vina_dock(self, mol):
        """ Docking only using vina
            Intended for use with too_large molecules
        """
        hostmol = self.hostmol
        exhaustiveness = self.conf["exhaustiveness"]
        n_poses = self.conf["n_poses"]
        min_rmsd = self.conf["min_rmsd"]
        hostpdbqtfile = self.conf["host_pdbqt"]

        preparator = MoleculePreparation(merge_these_atom_types=[])
        mol_setups = preparator.prepare(mol)

        # Initialize and setup vina optimizer
        vinaobj = vina.Vina(verbosity=0)
        vinaobj.set_receptor(hostpdbqtfile)
        vinaobj.compute_vina_maps(center=[0.0,0.0,0.0], box_size=[24.0, 24.0, 24.0])

        for setup in mol_setups:
            pdbqt_string, is_ok, err_msg = PDBQTWriterLegacy.write_string(setup, bad_charge_ok=True)

        if not is_ok:
            print(err_msg)

        vinaobj.set_ligand_from_string(pdbqt_string)

        _ = vinaobj.optimize()
        vinaobj.dock(exhaustiveness=exhaustiveness, n_poses=n_poses, min_rmsd=min_rmsd)

        vina_output_str = vinaobj.poses()
        if vina_output_str == "":
            raise ValueError("No poses found")

        pdbqt_mol = PDBQTMolecule(vina_output_str, is_dlg=False, skip_typing=True)
        rdkitmol = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)

        hostmol_confs = Chem.Mol(hostmol)
        hostmol_confs.RemoveAllConformers()
        for i in range(rdkitmol[0].GetNumConformers()):
            hostmol_confs.AddConformer(hostmol.GetConformer(), assignId=True)
            
        complex = Chem.CombineMols(hostmol_confs, rdkitmol[0])

        return complex, vinaobj.energies()

    """HELPER FUNCTIONS FOR DOCKING METHODS"""
    
    def align_mol(self, mol):
        """ Align principal axis of a molecule along the z axis 
        """
        guest_atoms = [atm.GetSymbol() for atm in mol.GetAtoms()]
        guest_coords = np.array([mol.GetConformer().GetAtomPosition(atm.GetIdx()) for atm in mol.GetAtoms()])

        #add clouds of points around the atoms of the guest molecule, 3 sections in polar (theta) and 6 in azimuthal (phi)
        #This prepares the coordinates for the PCA
        #   x = r*sin(theta)*cos(phi)
        #	y = r*sin(theta)*sin(phi)
        #	z = r*cos(theta)
        atomic_radii = {'H' :1.2, 'C':1.7, 'N':1.55, 'O':1.52, 'F':1.47, 'Cl':1.75, 'Br':1.85, 'I':1.98, 'P':1.8, 'S':1.8, 'As':1.85, 'B':2.13, 'Si':2.1, 'Se':1.9, 'Te':2.06}
        
        cloud_points = [] 
        for index,atm in enumerate(guest_atoms):
            cloud_radius = atomic_radii[atm]
            # Add a cloud point above and below the atom
            cloud_points.append(guest_coords[index] + np.array([0,0,cloud_radius]))
            cloud_points.append(guest_coords[index] + np.array([0,0,-cloud_radius]))
            # Add a cloud of points around the atom, 3 in theta and 6 in phi
            for theta in np.linspace(0,np.pi,3):
                for phi in np.linspace(0,2*np.pi,6):
                    cloud_points.append(guest_coords[index] + np.array([cloud_radius*np.sin(theta)*np.cos(phi),cloud_radius*np.sin(theta)*np.sin(phi),cloud_radius*np.cos(theta)]))

        # Add cloud points to the guest coordinates
        cloud_points = np.array(cloud_points)
        guest_coords_with_clouds = np.concatenate((guest_coords,cloud_points),axis=0)

        # Initiliase PCA
        pca = PCA(n_components=3)
        # Fit PCA to guest coordinates and transform the coordinates
        pca.fit_transform(guest_coords_with_clouds)
        transform_coord = pca.transform(guest_coords)

        # Direct the principal axis of the guest molecule towards the z-axis (the axis pointing through the cavity of the host)
        theta = np.arctan2(transform_coord[0,0],transform_coord[0,2])
        rotation_matrix = np.array([[np.cos(theta),0,np.sin(theta)],[0,1,0],[-np.sin(theta),0,np.cos(theta)]])
        transform_coord = np.matmul(rotation_matrix,transform_coord.T).T

        # Centre the transformed coordinates on the host centroid (0,0,0)
        transform_coord_centered = transform_coord.copy()
        transform_coord_centered[:,0] = transform_coord[:,0] - np.mean(transform_coord[:,0])
        transform_coord_centered[:,1] = transform_coord[:,1] - np.mean(transform_coord[:,1])
        transform_coord_centered[:,2] = transform_coord[:,2] - np.mean(transform_coord[:,2])

        return transform_coord_centered

    def transform_coords(self, coords, rotate_by, translate_by):
        """ Transforms the coordinates
            Rotates up to 90 degrees in one direction (symmetry of CBs)
            and translates along z-axis up to +4 (cavity has height ~9 A)
        """
        # Rotate coordinates about y axis
        rotation_matrix = np.array([[1,0,0],[0,np.cos(rotate_by),-np.sin(rotate_by)],[0,np.sin(rotate_by),np.cos(rotate_by)]])
        coords = np.matmul(rotation_matrix,coords.T).T

        # Translate coordinates
        coords[:,1] = coords[:,1] + translate_by

        return coords

    def modify_pdbqt_str(self, pdbqt_str, new_coords):
        """ Changes the coordinates of a pdbqt string to new coordinates
        """
        pdbqt_lines = pdbqt_str.split('\n')
        new_pdbqt_lines = []
        i = 0
        for line in pdbqt_lines:
            if line[:4]=='ATOM': 
                line = line[:30] + '%8.3f%8.3f%8.3f' % (new_coords[i,0], new_coords[i,1], new_coords[i,2]) + line[54:]
                i += 1
            new_pdbqt_lines.append(line)
        new_pdbqt_str = '\n'.join(new_pdbqt_lines)
        return new_pdbqt_str

    def vina_opt(self, guestmol, vinaobj):
        preparator = MoleculePreparation(merge_these_atom_types=[])
        mol_setups = preparator.prepare(guestmol)
        for setup in mol_setups:
            pdbqt_string, is_ok, err_msg = PDBQTWriterLegacy.write_string(setup, bad_charge_ok=True)

        if not is_ok:
            print(err_msg)
        
        vinaobj.set_ligand_from_string(pdbqt_string)
        opt_ens = vinaobj.optimize()

        with tempfile.NamedTemporaryFile(dir=self.conf["output_dir"]) as tf:
            vinaobj.write_pose(tf.name, overwrite=True, remarks="")
            pdbqt_mol = PDBQTMolecule.from_file(tf.name, skip_typing=True)

        # Returns a list of rdkit mols
        rdkitmols = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
        
        return opt_ens[0], rdkitmols[0]

    def vina_scoring(self, guestmol, vinaobj):
            preparator = MoleculePreparation()
            mol_setups = preparator.prepare(guestmol)
            for setup in mol_setups:
                pdbqt_string, is_ok, err_msg = PDBQTWriterLegacy.write_string(setup, bad_charge_ok=True)

            if not is_ok:
                print(err_msg)
        
            vinaobj.set_ligand_from_string(pdbqt_string)
            vina_ens = vinaobj.score()

            return vina_ens[0]
    
    def get_best_pose(self, complexmols, scores):
        """ From list of poses returns the best pose
            Rejects poses that are exo (unless its the final pose, which is caught by chem_sim)
        """
         # Complexes are returned as conformers
        for i in range(complexmols.GetNumConformers()):
            exo = is_exo(complexmols, self.hostmol, self.conf, confId=i)

            if not exo or i == complexmols.GetNumConformers()-1:
                complexmol = Chem.Mol(complexmols)
                complexmol.RemoveAllConformers()
                complexmol.AddConformer(complexmols.GetConformer(i), assignId=True)
                break

        return complexmol

if __name__ == "__main__":
    # working on this on branch docking_module_edit
    # access it with git checkout docking_module_edit
    # push it with git -u origin docking_module_edit
    # merge later on
    from smi2sdf import process_smi
    from reward_utils import is_exo
    from glob import glob
    import time
    import pandas as pd

    # Grab benchmark guests
    mols = [Chem.MolFromMolFile(i, removeHs=False, sanitize=False) for i in glob("./../data/top_10_spartan/*guest*")]
    hostmol = Chem.MolFromMolFile("data/host_aligned.sdf", removeHs=False, sanitize=False)
    hostpdbqtfile = "data/host_aligned.pdbqt"
    print(mols)

    smis = [Chem.MolToSmiles(i) for i in mols]
    smis = ["c12=CC=CC=c1cccc2", "CC"]
    guests = [process_smi(Chem.MolFromSmiles(i), 1, 0.35) for i in smis]
    print("confs generated")

    # Code for testing docking
    conf = {"vina_num_rotations":4,
            "vina_num_translations":4, 
            "centroid_diff_threshold": 4,
            "cavity_atoms_threshold": 6,
            "exhaustiveness": 32,
            "n_poses": 10,
            "min_rmsd": 3.0,
            "host_pdbqt":"data/host_aligned.pdbqt"}

    # Add a conformer to the host of each generated pose
    n_poses = (conf["vina_num_rotations"])*(conf["vina_num_translations"])
    for i in range(n_poses+1):
        newhost = Chem.Conformer(hostmol.GetConformer(0))
        hostmol.AddConformer(newhost, assignId=True)

    dock = DockLigand(hostmol, conf)
    
    # Comparing vina opt and vina dock
    print("vina docking")
    vina_dock = pd.DataFrame(columns=["complex", "score", "is_exo"])
    start_time_dock = time.time()
    for count, guestmol in enumerate(guests):
        complexmols, scores = dock.vina_dock(guestmol)
        pose_scores = []
        pose_exo = []
        for i in range(complexmols.GetNumConformers()):
            complexmol = Chem.Mol(complexmols)
            complexmol.RemoveAllConformers()
            complexmol.AddConformer(complexmols.GetConformer(i), assignId=True)
            exo = is_exo(complexmol, hostmol, conf, confId=0)

            pose_scores.append(scores[i][0])
            pose_exo.append(exo)

        vina_dock = vina_dock.append({"complex":complexmol, "smiles":smis[count], "score":pose_scores, "is_exo":pose_exo}, ignore_index=True)
        print(count, "done")
    end_time_dock = time.time()

     # Comparing vina opt and vina dock
    print("vina scoring")
    vina_score = pd.DataFrame(columns=["complex", "score", "is_exo"])
    start_time_score = time.time()
    for count, guestmol in enumerate(guests):
        complexmols, scores = dock.score_map_vina(guestmol)
        pose_scores = []
        pose_exo = []
        for i in range(complexmols.GetNumConformers()):
            complexmol = Chem.Mol(complexmols)
            complexmol.RemoveAllConformers()
            complexmol.AddConformer(complexmols.GetConformer(i), assignId=True)

            exo = is_exo(complexmol, hostmol, conf, confId=0)

            pose_scores.append(scores[i])
            pose_exo.append(exo)

        vina_score = vina_score.append({"complex":complexmol, "smiles":smis[count], "score":pose_scores, "is_exo":pose_exo}, ignore_index=True)
        print(count, "done")
    end_time_score = time.time()

    # Calculate the elapsed time for each loop
    elapsed_time_dock = end_time_dock - start_time_dock
    elapsed_time_score = end_time_score - start_time_score

    print("Elapsed time for vina_dock loop:", elapsed_time_dock, "seconds")
    print("Elapsed time for vina_score loop:", elapsed_time_score, "seconds")
    print(vina_dock)
    print(vina_score)

    pd.to_pickle(vina_dock, "vina_dock_results.pkl")
    pd.to_pickle(vina_score, "vina_score_results.pkl")
