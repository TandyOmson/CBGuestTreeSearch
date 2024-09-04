""" Handles docking via a three step process
    1. Pose generation of initial guesses
    2. Pre optimisation via cheap forcefield
    3. Scoring and selection of best poses
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms
import vina
from meeko import MoleculePreparation, PDBQTWriterLegacy, PDBQTMolecule, RDKitMolCreate
from scipy.spatial import Delaunay
from sklearn.decomposition import PCA
import math
import numpy as np

from glob import glob
import os
import tempfile
import subprocess as sp

"""HELPER FUNCTIONS FOR DOCKING METHODS"""

def is_small_cylinder(guestmol, cutoff_r=3.65, cutoff_h=4.55):
    """ Determines if the guest is small
    """
    # Find length of guest molecule by finding the distance between the centroid and the furthest atom
    centroid = rdMolTransforms.ComputeCentroid(guestmol.GetConformer())
    max_dist = 0
    vector = None
    
    for atom in guestmol.GetAtoms():
        if atom.GetAtomicNum() != 1:
            # Get the distance between the centroid and the atom
            dist = np.linalg.norm(guestmol.GetConformer().GetAtomPosition(atom.GetIdx()) - centroid)
            if dist > max_dist:
                max_dist = dist
                # Get the vector between the centroid and the furthest atom
                vector = guestmol.GetConformer().GetAtomPosition(atom.GetIdx()) - centroid
    
    # Find the radius of the molecule as the furthest atom from the centroid in the plane perpendicular to the vector
    max_rad = 0
    for atom in guestmol.GetAtoms():
        if atom.GetAtomicNum() != 1:
            # Vector projection of the atom onto the plane perpendicular to the vector
            v = guestmol.GetConformer().GetAtomPosition(atom.GetIdx()) - centroid
            proj = (np.dot(v, vector) / np.dot(vector, vector)) * vector
            x = centroid + v - proj
            # Get the distance between the projected point and the centroid
            proj_rad = np.linalg.norm(x-centroid)
            if proj_rad > max_rad:
                max_rad = proj_rad
    
    if max_rad < cutoff_r and max_dist < cutoff_h:
        return True
    else:
        return False


def is_exo(mol, confId=-1):
    """ Checks for exo complex
    """
    centroiddiffthreshold=4
    cavityatomsthreshold=6
    
    hostmol = Chem.GetMolFrags(mol, asMols=True)[0]

    Chem.RemoveHs(mol)
    Chem.RemoveHs(hostmol)

    # Separate host and guest, get their coordinates
    guest_coords = np.array([mol.GetConformer().GetAtomPosition(atm.GetIdx()) for count, atm in enumerate(mol.GetAtoms()) if count >= hostmol.GetNumAtoms()])
    host_coords = np.array([mol.GetConformer().GetAtomPosition(atm.GetIdx()) for count, atm in enumerate(mol.GetAtoms()) if not count >= hostmol.GetNumAtoms()])

    # Get host and guest centroid
    guest_centroid = np.array(guest_coords).mean(axis=0)
    host_centroid = np.array(host_coords).mean(axis=0)

    # Calculate distance between guest and host centroid
    centroid_diff = np.linalg.norm(guest_centroid - host_centroid)

    # Delauny defines the convex hull of the host atoms
    # Delauny is a triangulation such that none of the host atoms are inside the circumsphere of any tetrahedron in the traingulation
    hull = Delaunay(host_coords)
    
    # Calculate number of atoms in cavity
    cavity_atoms = 0
    for atm in guest_coords:
        # Points outside the triangulation return -1
        if hull.find_simplex(atm) >= 0:
            cavity_atoms += 1

    # Check if exclusion complex
    isExo = False
    if centroid_diff > centroiddiffthreshold or cavity_atoms < cavityatomsthreshold:
        isExo = True

    return isExo

def align_mol(mol):
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

def transform_coords(coords, rotate_by, translate_by):
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

def modify_pdbqt_str(pdbqt_str, new_coords):
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

class DockMol():
    """ Class that stores all the information involved in docking on rdkit mol
        In addition to object required during pose generation and pose scoring
    """
    
    def __init__(self, mol, **kwargs):
        self.mol = mol

        # complex are keys, values are scores
        self.complexes = {}

        # Leaving kwargs in for redundancy
        self.__dict__.update(kwargs)

class DockLigand():
    """ Class for carrying out docking operations
        Requires guestmol, hostmol and conf
        The required contents of conf depends on the docking method
    """

    def __init__(self, conf):
        self.conf = conf
        # Non optional variables assigned as class attributes
        self.hostmol = Chem.MolFromMolFile(conf["host_sdf"], removeHs=False)
        self.hostpdbqtfile = conf["host_pdbqt"]

    """ Calculator functions """

    def get_docked_complex(self, mol):
        """ Docks a guest molecule mol
        """
        try:
            dockmol = self.get_poses(mol)
            complexmol = self.get_best_pose(dockmol)
        except:
            raise ValueError("Error in docking")

        return complexmol

    def get_poses(self, mol):
        """ Gets a list of poses for a molecule
            Stores them in the "complexes" attribute of DockMol object
        """
        # Decide whether to use alignment docking
        try:
            is_small = is_small_cylinder(mol)
        except:
            is_small = False
        
        dockmol = DockMol(mol)
        
        # If the molecule is small, apply alignment docking
        if is_small:
            self.alignment_poses(dockmol, self.conf["vina_num_rotations"], self.conf["vina_num_translations"])
            self.MMFF94_vina_opt(dockmol)
        else:
            self.vina_dock(dockmol, n_poses=10, min_rmsd=0.5)

        return dockmol

    def get_best_pose(self, dockmol):
        """ Chooses the best pose that is an "endo" complex
        """
        for key, complex_pose in sorted(dockmol.complexes.items()):
            complexmol = None
            if not is_exo(complex_pose):
                complexmol = complex_pose
                break

            if not complexmol:
                raise ValueError(f"No endo complexes found ({len(dockmol.complexes)} exo)")
        
        return complexmol

    """ Internal docking functions """

    def alignment_poses(self, dockmol, num_rot, num_tra):
        """ Generates binding poses via alignment docking
            Adds them to conformers of dockmol
        """
        hostmol = self.hostmol
        hostpdbqtfile = self.hostpdbqtfile
        
        aligned_coords = align_mol(dockmol.mol)
        
        # Create 2d grid of rotations and translations
        rotations = np.linspace(0,90,num_rot)
        translations = np.linspace(0,4,num_tra)

        # Conformer 2: flip the guest 90 degrees about the x axis (means it probably wont fit, so will not screen in this direction)
        rotation_matrix = np.array([[np.cos(90),-np.sin(90),0], [np.sin(90), np.cos(90), 0], [0,0,1]])
        rotation_coords = np.matmul(rotation_matrix,aligned_coords.T).T

        rotconf = Chem.Conformer(dockmol.mol.GetConformer(0))
        for index, atm in enumerate(dockmol.mol.GetAtoms()):
            rotconf.SetAtomPosition(atm.GetIdx(), rotation_coords[index])
        dockmol.mol.AddConformer(rotconf, assignId=True)

        # Screen across rotations and translations
        for rotate_by in rotations:
            for translate_by in translations:
                new_coords = transform_coords(aligned_coords, rotate_by=rotate_by, translate_by=translate_by)

                # Assign conformer
                traconf = Chem.Conformer(dockmol.mol.GetConformer(0))
                for index, atm in enumerate(dockmol.mol.GetAtoms()):
                    traconf.SetAtomPosition(atm.GetIdx(), new_coords[index])
                dockmol.mol.AddConformer(traconf, assignId=True)

        # Combine each of the poses with the host and calculate scores
        hostmols = Chem.Mol(hostmol)
        hostmols.RemoveAllConformers()
        for i in range(dockmol.mol.GetNumConformers()):
            hostmols.AddConformer(hostmol.GetConformer(0), assignId=True)

        binding_poses = Chem.CombineMols(hostmols, dockmol.mol)
        Chem.SanitizeMol(binding_poses)

        # Add new class attribute if one doesn't exist already
        if hasattr(dockmol, "binding_poses"):
            for conf in binding_poses.GetConformers():
                dockmol.binding_poses.AddConformer(conf, assignId=True)
        else:
            setattr(dockmol, "binding_poses", binding_poses)
            

    def MMFF94_vina_opt(self, dockmol):
        """ Runs MMFF94, then optimisation on a binding pose
            Modifies the pose in the conformers and adds a vina score
        """
        hostpdbqtfile = self.hostpdbqtfile

        mmff_res = AllChem.MMFFOptimizeMoleculeConfs(dockmol.binding_poses, numThreads=0, ignoreInterfragInteractions=False)

        # In order to do vina operations, I need to use separate mol objects for each guest pose
        mols = []
        for conf in dockmol.mol.GetConformers():
            mol = Chem.Mol(dockmol.mol)
            mol.RemoveAllConformers()
            mol.AddConformer(conf)
            mols.append(mol)

        for mol in mols:
            preparator = MoleculePreparation()
            mol_setups = preparator.prepare(mol)
            
            for setup in mol_setups:
                pdbqt_string, is_ok, err_msg = PDBQTWriterLegacy.write_string(setup, bad_charge_ok=True)

            if not is_ok:
                print(err_msg)

            vinaobj = vina.Vina(verbosity=0)
            vinaobj.set_receptor(hostpdbqtfile)
            vinaobj.compute_vina_maps(center=[0.0,0.0,0.0], box_size=[24.0, 24.0, 24.0])
        
            vinaobj.set_ligand_from_string(pdbqt_string)
            vina_decomp = vinaobj.optimize()
            vina_en = vina_decomp[0]
            
            with tempfile.NamedTemporaryFile(dir=self.conf["output_dir"]) as tf:
                vinaobj.write_pose(tf.name, overwrite=True, remarks="")
                pdbqt_mol = PDBQTMolecule.from_file(tf.name, skip_typing=True)

            complexmol = Chem.CombineMols(self.hostmol, RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)[0])

            # Returns an rdkit mol
            dockmol.complexes[vina_en] = complexmol

    """VINA DOCKING VIA GENERATIVE METHOD"""

    def vina_dock(self, dockmol, n_poses, min_rmsd):
        """ Docking only using vina
            Intended for use with too_large molecules
        """
        hostmol = self.hostmol
        exhaustiveness = self.conf["exhaustiveness"]
        n_poses = self.conf["n_poses"]
        min_rmsd = self.conf["min_rmsd"]
        hostpdbqtfile = self.conf["host_pdbqt"]

        preparator = MoleculePreparation(merge_these_atom_types=[])
        mol_setups = preparator.prepare(dockmol.mol)

        # Initialize and setup vina optimizer
        vinaobj = vina.Vina(verbosity=0)
        vinaobj.set_receptor(hostpdbqtfile)
        vinaobj.compute_vina_maps(center=[0.0,0.0,0.0], box_size=[24.0, 24.0, 24.0])

        for setup in mol_setups:
            pdbqt_string, is_ok, err_msg = PDBQTWriterLegacy.write_string(setup, bad_charge_ok=True)

        if not is_ok:
            print(err_msg)

        vinaobj.set_ligand_from_string(pdbqt_string)
        
        vinaobj.dock(exhaustiveness=exhaustiveness, n_poses=n_poses, min_rmsd=min_rmsd)
        
        vina_output_str = vinaobj.poses(n_poses=n_poses)
                
        if vina_output_str == "":
            raise ValueError("No poses found")
        pdbqt_mol = PDBQTMolecule(vina_output_str, is_dlg=False, skip_typing=True)
        rdkitmol = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)[0]
        vina_ens = vinaobj.energies(n_poses=n_poses)[:,0]

        for conf, vina_en in zip(rdkitmol.GetConformers(), vina_ens):
            mol = Chem.Mol(rdkitmol)
            mol.RemoveAllConformers()
            mol.AddConformer(conf, assignId=True)
            
            complexmol = Chem.CombineMols(self.hostmol, mol)
            # Returns an rdkit mol
            dockmol.complexes[vina_en] = complexmol
    
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

    smis = [""]
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
            "host_pdbqt":"data/host_aligned.pdbqt",
            "host_sdf":"data/host_aligned.sdf",
            "output_dir":"./temp",
    }

    dock = DockLigand(conf)
    
    # Comparing vina opt and vina dock
    print("vina docking")
    vina_dock = pd.DataFrame(columns=["complex", "score", "is_exo"])
    start_time_dock = time.time()
    for count, guestmol in enumerate(guests):
        dockmol = dock.get_poses(guestmol)
        complexmol = dock.get_best_pose(dockmol)
        Chem.MolToMolFile(complexmol, "./temp/tempo2.sdf")
        print(count, "done")
        
    end_time_dock = time.time()

    """
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

            exo = is_exo(complexmol, confId=0)

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
`    """
