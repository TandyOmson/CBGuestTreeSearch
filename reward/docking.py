from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

import numpy as np
from sklearn.decomposition import PCA

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

""" MMFF94 METHOD """

def score_map_mmff94(mol, hostmol, num_rot, num_tra):
    """ Uses the mmff94 force field to optimise just a few poses, instead of using vina score.
    """

    #print("DOCKING - aligning guest")
    # Conformer 1: principal axis of guest aligned with axis through the cavity of the CB (z-axis)
    aligned_coords = align_mol(mol)
    
    # Conformer 2: flip the guest 90 degrees about the x axis (means it probably wont fit, so will not screen in this direction)
    rotation_matrix = np.array([[np.cos(90),-np.sin(90),0], [np.sin(90), np.cos(90), 0], [0,0,1]])
    rotation_coords = np.matmul(rotation_matrix,aligned_coords.T).T

    rotconf = Chem.Conformer(mol.GetConformer(0))
    for index, atm in enumerate(mol.GetAtoms()):
        rotconf.SetAtomPosition(atm.GetIdx(), rotation_coords[index])
    mol.AddConformer(rotconf, assignId=True)
    
    # Create 2d grid of rotations and translations
    rotations = np.linspace(0,90,num_rot)
    translations = np.linspace(0,4,num_tra)

    #print("DOCKING - performing transformations")

    # Screen across rotations and translations
    for rotate_by in rotations:
        for translate_by in translations:
            # Add the new conformers as conformers to the mol object
            new_coords = transform_coords(aligned_coords, rotate_by=rotate_by, translate_by=translate_by)
            # Get a copy of the original conformer to edit with the transformed coordinates
            conf = Chem.Conformer(mol.GetConformer(0))
            for index, atm in enumerate(mol.GetAtoms()):
                conf.SetAtomPosition(atm.GetIdx(), new_coords[index])
            mol.AddConformer(conf, assignId=True)

    # Combine each of the poses with the host and calculate scores
    complexconfs = Chem.CombineMols(hostmol, mol)

    # Add relevant info for optimisation
    Chem.GetSSSR(complexconfs)

    #print("DOCKING - performing MMFF94 optimsiation on poses")
    
    # Optimise the conformers    
    mmff_res = AllChem.MMFFOptimizeMoleculeConfs(complexconfs, numThreads=0)
    scores = np.array([i[1] for i in mmff_res])
    
    #min_score = np.min(scores)
    #min_conf = complexconfs.GetConformer(int(np.where(scores==min_score)[0][0]))

    #print("DOCKING - scoring poses")

    # Sort the conformers by score and assign IDs in ascending order
    score_confids = [int(i) for i in np.argsort(scores)]
    scores = [float(i) for i in np.sort(scores)]
    
    finalcomplex = Chem.Mol(complexconfs)
    finalcomplex.RemoveAllConformers()
    for i in score_confids:
        finalcomplex.AddConformer(complexconfs.GetConformer(i), assignId=True)
    
    return finalcomplex, scores

""" VINA METHODS """

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

def score_map_vina(mol, hostmol, num_rot, num_tra, vinaobj, pdbqt_str):
    """ Screens across a set of rotations and translations
    """
    aligned_coords = align_mol(mol)
    
    meeko_prep = MoleculePreparation(merge_these_atom_types=[])
    mol_setups = meeko_prep.prepare(mol)
    for setup in mol_setups:
        pdbqt_setup, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
        if is_ok:
            pdbqt_setup_final = pdbqt_setup
    pdbqt_str = modify_pdbqt_str(pdbqt_setup_final, align_mol(mol))
    aligned_coords = align_mol(mol)
    
    # Create 2d grid of rotations and translations
    rotations = np.linspace(0,90,num_rot)
    translations = np.linspace(0,4,num_tra)
    scores = np.zeros((num_rot,num_tra))

    # Screen across rotations and translations
    for i,rotate_by in enumerate(rotations):
        for j,translate_by in enumerate(translations):
            new_coords = transform_coords(aligned_coords, rotate_by=rotate_by, translate_by=translate_by)
            vinaobj.set_ligand_from_string(modify_pdbqt_str(pdbqt_str, new_coords))
            vinaobj.compute_vina_maps(center=[0.0,0.0,0.0], box_size=[24.0, 24.0, 24.0])
            scores[i,j] = vinaobj.score()[0]

    # Find minimum score, and the rotation and translation that gave it
    min_score = np.min(scores)
    min_score_index = np.where(scores==min_score)

    min_transformation = (rotations[min_score_index[0][0]], translations[min_score_index[1][0]])
    transformed_coords = transform_coords(aligned_coords, min_transformation[0], min_transformation[1])
    
    newmol = Chem.RWMol(mol)
    for index,atm in enumerate(mol.GetAtoms()):
            newmol.GetConformer().SetAtomPosition(atm.GetIdx(),transformed_coords[index])
            min_conf = newmol.GetConformer()
    newmol.RemoveAllConformers()
    newmol.AddConformer(min_conf)

    min_complex = Chem.CombineMols(hostmol, newmol)

    return min_complex, min_score