""" METHODS FOR REWARD CALCULATION
fix_charge: looks for atoms with mismatching explicit valences and charges them (for example N)

get_opt: optimises the xTB energy and retrieves energy, returning the optimised molecule and the energy
xtb_opt: caller for xtb optimisation
get_binding: reads binding energy (can add other properties here)

is_small_cylinder: determines if the guest is small
isExo: determines if the guest is exohedrally bound
"""

from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import PropertyMol
from rdkit.Chem import rdDetermineBonds
from scipy.spatial import Delaunay
import math

import subprocess as sp
import os
import numpy as np
from shutil import rmtree
from glob import glob

""" MOLECULE PREP METHODS """

def fix_charge(mol):
    """ looks for atoms with mismatching explicit valences and charges them (for example N) 
    """
    mol.UpdatePropertyCache(strict=False)
    for at in mol.GetAtoms():
        if at.GetAtomicNum() == 7 and at.GetExplicitValence()==4 and at.GetFormalCharge()==0:
            at.SetFormalCharge(1)
    Chem.SanitizeMol(mol)

""" CHECKS ON MOLECULES """

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
    
def is_exo(complexmol, hostmol, conf, confId=-1):
    """ Checks for exo complex
    """
    if not conf["centroid_diff_threshold"]:
        centroiddiffthreshold=4
        cavityatomsthreshold=6
    else: 
        centroiddiffthreshold = conf["centroid_diff_threshold"]
        cavityatomsthreshold = conf["cavity_atoms_threshold"]

    Chem.RemoveHs(complexmol)
    Chem.RemoveHs(hostmol)

    # Separate host and guest, get their coordinates
    guest_coords = np.array([complexmol.GetConformer(confId).GetAtomPosition(atm.GetIdx()) for count, atm in enumerate(complexmol.GetAtoms()) if count >= hostmol.GetNumAtoms()])
    host_coords = np.array([complexmol.GetConformer(confId).GetAtomPosition(atm.GetIdx()) for count, atm in enumerate(complexmol.GetAtoms()) if not count >= hostmol.GetNumAtoms()])

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

def check_smiles_change(beforemol, aftermol):
    """ Post docking and optimisation, checks if the optimised guest is the same as the original guest
        i.e. no bond changes. If this is the case, the molecule should be flagged and get a bad score.
    """
    # Get the smiles of the original guest
    before_smiles = Chem.MolToSmiles(beforemol, canonical=True)

    # Get the smiles of the optimised guest
    conn_mol = Chem.Mol(aftermol)
    # Simple determination of connectivity based on bond lengths and chemistry rules
    # This should pick up chemical changes in the molecule
    rdDetermineBonds.DetermineConnectivity(conn_mol)
    try:
        rdDetermineBonds.DetermineBondOrders(conn_mol)
        after_smiles = Chem.MolToSmiles(conn_mol, canonical=True)
    except:
        is_changed = True
        return is_changed

    # Check if the smiles are the same
    is_changed = False
    if before_smiles != after_smiles:
        is_changed = True

    return is_changed

def covalent_CB(mol):
    covalent_CB = False

    frags = Chem.GetMolFrags(mol, asMols=True)

    # Check if any atom in the guest is within 1.4 angstroms of any atom in the host
    for count1, positions in enumerate(frags[1].GetConformer().GetPositions()):
        for count2, host_positions in enumerate(frags[0].GetConformer().GetPositions()):
            if np.linalg.norm(positions - host_positions) < 1.5:
                covalent_CB = True
                return covalent_CB

    return covalent_CB

def get_incorrect_bond_angle(mol):
    """ Takes a mol and returns whether there are any tight angles between any two bonds
    """
    mol = Chem.GetMolFrags(mol, asMols=True)[-1]

    incorrect_angle = False
    mol = Chem.RemoveHs(mol)
    for atom in mol.GetAtoms():
        if atom.GetDegree() > 1:
            # Run through all the permuatation of the other atom indices of the bonds associated with this atom
            atom_bonds = []
            for bond in atom.GetBonds():
                atom_bonds.append(bond.GetOtherAtomIdx(atom.GetIdx()))
            for i in atom_bonds:
                for j in atom_bonds:
                        if i != j:
                            if rdMolTransforms.GetAngleDeg(mol.GetConformer(), i, atom.GetIdx(), j) < 30.0:
                                incorrect_angle = True
                                return incorrect_angle

    return incorrect_angle

def get_incorrect_bond_length(mol):
    mol = Chem.GetMolFrags(mol, asMols=True)[-1]

    wrong_bond_length = False
    for bond in mol.GetBonds():
        bond_length = np.linalg.norm(mol.GetConformer().GetAtomPosition(bond.GetBeginAtom().GetIdx()) -
                                mol.GetConformer().GetAtomPosition(bond.GetEndAtom().GetIdx()))
        if bond_length > 1.6:
            wrong_bond_length = True
            return wrong_bond_length

    return wrong_bond_length

def get_property_mol(mol):
    """ Helper function: for regular RDMol objects, properties are deleted when the molecule is modified or pickled.
    This function returns a copy of the molecule as a PropertyMol, to be used when returning a Mol.
    """
    props = mol.GetPropsAsDict()
    propmol = PropertyMol.PropertyMol(mol)
    for key,val in props.items():
        propmol.SetProp(key, val)

    return propmol
