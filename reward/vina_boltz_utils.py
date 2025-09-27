from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import rdMolAlign
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.rdchem import PeriodicTable

import vina
from meeko import MoleculePreparation, PDBQTWriterLegacy, PDBQTMolecule, RDKitMolCreate

import random
from collections import defaultdict
from scipy.spatial import Delaunay
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA
import numpy as np
import math

import subprocess as sp
import tempfile
import re
import json

# thermodynamic temp. in kcal/mol (at 298.15 K)
beta = 1.68787

# Disable RDKit logging
RDLogger.DisableLog('rdApp.*')

def add_nitrogen_charges(m):
    m.UpdatePropertyCache(strict=False)
    ps = Chem.DetectChemistryProblems(m)
    if not ps:
        Chem.SanitizeMol(m)
        return m
    for p in ps:
        if p.GetType()=='AtomValenceException':
            at = m.GetAtomWithIdx(p.GetAtomIdx())
            if at.GetAtomicNum()==7 and at.GetFormalCharge()==0 and at.GetExplicitValence()==4:
                at.SetFormalCharge(1)
            if at.GetAtomicNum()==7 and at.GetFormalCharge()==0:
                bondcount = 0
                for b in at.GetBonds():
                    bondcount += b.GetBondTypeAsDouble()
                if int(boundcount) > 3:
                    at.SetFormalCharge(1)
                
    Chem.SanitizeMol(m)
    return m    

#########
# CREST #
#########
def stratified_weighted_sample(lst, total_samples=5, bias_strength=1.0):
    """
    Stratified sampling of list indices, slightly biased toward lower indices, always includes index 0.
    
    Parameters:
    - lst: List of values
    - total_samples: Total number of indices to sample (including index 0)
    - bias_strength: Float > 0. Higher values bias more toward the front
    
    Returns:
    - Sorted list of sampled indices
    """
    if total_samples < 1:
        return []

    # Stratify indices by value
    strata = defaultdict(list)
    for i, val in enumerate(lst):
        strata[val].append(i)

    sampled_indices = set([0])  # Always include the first index
    remaining_samples = total_samples - 1

    # Flatten all indices except 0
    all_other_indices = [i for i in range(1, len(lst))]
    
    # Create weights: lower index → higher weight
    weights = [1 / ((i + 1) ** bias_strength) for i in all_other_indices]
    total_weight = sum(weights)
    probabilities = [w / total_weight for w in weights]

    # Draw samples from remaining indices
    if remaining_samples > 0 and all_other_indices:
        sampled = random.choices(all_other_indices, weights=probabilities, k=remaining_samples)
        sampled_indices.update(sampled)

    return sorted(sampled_indices)

########
# VINA #
########
def write_out_confs(confs, writer, n_poses):
    n, i = 0, 0
    while n < n_poses:
        try:
            writer.write(confs[i])
            n += 1
            i += 1
        except:
            break
    return n

def parse_vina_output(filename, exo_ignore):
    affinities = []
    rmsd_lbs = []
    rmsd_ubs = []

    with open(filename, 'r') as file:
        lines = file.readlines()

    parsing = False
    count = 0
    for line in lines:
        # Start parsing after the table header
        if re.match(r'^\s*-+\+-+', line):
            parsing = True
            continue
        if parsing:
            if line.strip().split()[0] == "Writing":
                break  # Stop at first empty line after data
            parts = line.strip().split()
            if len(parts) == 4:
                if count not in exo_ignore:
                    mode = int(parts[0])
                    affinity = float(parts[1])
                    rmsd_lb = float(parts[2])
                    rmsd_ub = float(parts[3])

                    affinities.append(affinity)
                    rmsd_lbs.append(rmsd_lb)
                    rmsd_ubs.append(rmsd_ub)
                count += 1

    return affinities, rmsd_lbs, rmsd_ubs

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
    if mol.GetNumAtoms() > 141:
        if centroid_diff > centroiddiffthreshold or cavity_atoms < cavityatomsthreshold:
            isExo = True
    else:
        if centroid_diff > centroiddiffthreshold:
            isExo = True

    return isExo

# ALIGNMENT DOCKING FUNCTIONS
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

def PCA_align_pose(hostmol, guest):
    """ Generates binding poses via alignment docking
    """
    aligned_coords = align_mol(guest)
    
    conf = Chem.Conformer(guest.GetConformer(0))
    for index, atm in enumerate(guest.GetAtoms()):
        conf.SetAtomPosition(atm.GetIdx(), aligned_coords[index])
    guest.RemoveAllConformers()
    guest.AddConformer(conf, assignId=True)

    binding_pose = Chem.CombineMols(hostmol, guest)
    binding_pose = add_nitrogen_charges(binding_pose)

    return binding_pose

def filter_by_rmsd(complexmols, vina_ens, n_poses):
    """ Filter complexmols by rmsd 
    """
    guests = [Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)[-1] for mol in complexmols]
    guests = [add_nitrogen_charges(mol) for mol in guests]

    match = guests[0].GetSubstructMatch(guests[0])
    atom_map = list(enumerate(match))

    N = len(guests)
    rmsd_matrix = np.zeros((N, N))
        
    # Compute RMSD between all pairs
    for i in range(N):
        for j in range(i+1, N):
            rmsd = rdMolAlign.GetBestRMS(guests[i], guests[j], map=[atom_map])
            rmsd_matrix[i, j] = rmsd
            rmsd_matrix[j, i] = rmsd

    # Convert to condensed form (1D) for linkage
    condensed = squareform(rmsd_matrix)
    Z = linkage(condensed, method='average')
    
    # Assign cluster labels (e.g., 5 clusters or a distance threshold)
    labels = fcluster(Z, t=n_poses, criterion='maxclust')

    # select from each cluster successively, keeping the structure that has an affinity that is not seen already
    seen_affins = set()
    selected_indices = []
    selected_affins = []
    for clust in range(1, len(np.unique(labels))+1):
        # Get indices in this cluster
        cluster_indices = np.where(labels == clust)[0]
        
        # Sort these indices by affinity
        sorted_cluster_indices = sorted(cluster_indices, key=lambda i: vina_ens[i])

        # Pick first unseen index
        for i in sorted_cluster_indices:
            aff = vina_ens[i]

            if aff not in seen_affins:
                selected_indices.append(int(i))
                selected_affins.append(float(aff))
                seen_affins.add(aff)
                break  # move to next cluster
            else:
                # if affins is already seen, just add the lowest
                selected_indices.append(int(sorted_cluster_indices[0]))
                selected_affins.append(float(vina_ens[0]))
                break

    complexmols = [complexmols[i] for i in selected_indices]
    vina_ens = [vina_ens[i] for i in selected_indices]
            
    return complexmols, vina_ens

def get_symmetry_number(sym_prog, mol, util_dir):
    try:
        atomic_numbers = [atm.GetAtomicNum() for atm in mol.GetAtoms()]
        positions = [[p[0], p[1], p[2]] for p in mol.GetConformer().GetPositions()]
    except:
        return 1

    with tempfile.NamedTemporaryFile(mode="w", delete=False) as fw:
        fw.write(f"{str(len(atomic_numbers))}\n\n")
        for atomic_num, p in zip(atomic_numbers, positions):
            atom_line = str(atomic_num)
            if p[0] > 0.0:
                atom_line += 10*" "
            else:
                atom_line += 9*" "
                
            atom_line += f"{p[0]:.5f}"

            if p[1] > 0.0:
                atom_line += 8*" "
            else:
                atom_line += 7*" "
                
            atom_line += f"{p[1]:.5f}"

            if p[2] > 0.0:
                atom_line += 8*" "
            else:
                atom_line += 7*" "
                
            atom_line += f"{p[2]:.5f}"

            atom_line += "\n"
            
            fw.write(atom_line)

        fw.close()

        # Still in tempfile context manager
        result = sp.run(sym_prog.split(" ") + [fw.name], stdout=sp.PIPE)

    # search for point group in output string
    m = re.search(r'symmetry elements:\s*(.*)', result.stdout.decode('utf-8'))
    if m:
        symmetry_elements = m.group(1)
    else:
        symmetry_elements = " "

    with open(f"{util_dir}/point_group_lookup.json", "r") as fr:
        point_group_table = json.load(fr)

    try:
        point_group = point_group_table[symmetry_elements]
    except:
        try:
            point_group = point_group_table[symmetry_elements.split(" ")[0] + " "]
        except:
            point_group = "C1"

    with open(f"{util_dir}/rot_num_lookup.json", "r") as fr:
        rot_num_table = json.load(fr)

    try:
        rot_num = rot_num_table[point_group]
    except:
        rot_num = 1
        
    return rot_num

def xb(en):
    return math.exp(-beta*en)

def get_energy_change(states, degens):
    Zr = 1
    for e, n in zip(states, [1]*len(states)):
        # "relative" partition function
        Zr += n*xb(e)
    
    dE = 0
    S = 1.98624e-3*math.log(Zr)
    for e, n in zip(states, degens):
        # "relative" partition function
        dE += n*e*xb(e)/Zr
        S += 1.98624e-3*beta*((n*e*xb(e))/Zr)
        
    return Zr, dE, S*298.15
