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
from scipy.spatial import Delaunay

import subprocess as sp
import os
import numpy as np

""" MOLECULE PREP METHODS """

def fix_charge(mol):
    """ looks for atoms with mismatching explicit valences and charges them (for example N) 
    """
    mol.UpdatePropertyCache(strict=False)
    for at in mol.GetAtoms():
        if at.GetAtomicNum() == 7 and at.GetExplicitValence()==4 and at.GetFormalCharge()==0:
            at.SetFormalCharge(1)
    Chem.SanitizeMol(mol)

""" XTB METHODS """

def get_opt(mol, outfile, conf):
    """ Calls methods to optimise a mol and retrieve energy
    """
    orgdir = os.getcwd()
    os.chdir(conf["xtb_tempdir"])

    Chem.MolToMolFile(mol,"mol.sdf",kekulize=False)
    xtb_opt("mol.sdf", outfile)
    try:
        finalmol = Chem.MolFromMolFile("xtbopt.sdf",removeHs=False,sanitize=False)
    except:
        print("didnt optimise")

    # Put the outfile somewhere? (Write a function to get the relevant information and deposit it in a filewriter)
    en = get_binding(outfile)
    os.chdir(orgdir)
    
    return finalmol, en

def xtb_opt(filename, outfile):
    """ Optimises from an sdf file """

    sp.run(["obabel","-isdf","-osdf",f"{filename}","-O",f"{filename}"],stderr=sp.DEVNULL)
    sp.run(["xtb",f"{filename}","--ohess","normal","--alpb","water"],stdout=open(outfile,"w"),stderr=sp.DEVNULL)
    
    return outfile

def get_binding(filename):
    """ Gets binding energy from an sdf file """

    gen = (i.split() for i in reversed(open(filename,"r").readlines()))

    binding = 0
    for i in gen:
        if i:
            if " ".join(i[:4]) == ":: total free energy":
                binding = float(i[4])*627.5095
                break

    return binding

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
    
def is_exo(complexmol, hostmol, conf, confId=0):
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
