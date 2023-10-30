# note: Can I store host information when I start the program? This would require the creation of module level (ie here) or global variables.
# Module level variables would be calculated and stored when this class is first referenced (ie when getattr() returns the reward class below)
import vina
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
import numpy as np
from scipy.spatial import Delaunay
import os
import subprocess as sp

# Methods for reward calculations
from reward.smi2sdf import process_smi
from reward.docking import score_map, score_map_mmff94
from reward.xtb_bind import xtb_opt, get_binding
from reward.sascorer import calculateScore

# Can inherit from the batchReward abstract class
from reward.reward import Reward

class CBDock_reward(Reward):
    """ Reward Class in config. A reference to this class from config is passed through the MCTS class to the evaluate_node method in utils 
        The following methods are static (no instantiation of class required) and abstract (must be included in the class) methods.
    """
    def get_objective_functions(conf):
        """ Must return a list of functions """

        def binding_en(mol):
            """ Calculate values contributing to the reward arising from binding
            """
            _host_initialised = False
            if not _host_initialised:
                initialise_host(conf)

            # 1. Generate conformers
            guestmol = smi2sdf(mol, conf)
            fix_charge(guestmol)

            # If molecule is small, set the "max score" (code 1000)
            try:
                is_small = is_small_cylinder(guestmol)
            except:
                print("bad conformer")
                return 20.2
                
            if not is_small:
                print("Guest too large")
                return 20.0
            
            # 2. Dock the best conformer
            # note that currently this takes the best result from docking. 
            # I may want to do a quick optimisation on the top few complexes if they are very close together, but have large difference in pose
            complexmol, min_score = score_map_mmff94(
                guestmol,
                hostmol,
                conf["vina_num_rotations"],
                conf["vina_num_translations"]
                )
            
            # 3. Calculate binding energy
            try:
                print("optimising guest")
                optguestmol, guest_en = get_opt(guestmol, "guestopt.out", conf)
                print("optimising complex")
                optcomplexmol, complex_en = get_opt(complexmol, "complexopt.out", conf)
            except:
                print("couldn't optimise")
                return 20.1

            try:
                optcomplexmol.SetDoubleProp("binding", complex_en)
                final_complex_writer.write(optcomplexmol)
            except:
                print("could not write complex")

            # If the resulting complex is endo, set the "max score" (code 1001)
            is_exo = isExo(optcomplexmol, hostmol, conf)
            if is_exo:
                print("Exo complex")
                return 20.2
            
            return complex_en - guest_en - host_en

        def sa_scorer(mol):
            """ Calculate synthetic accessibility score
            """
            sa_score = calculateScore(mol)
            
            return sa_score
            
        return [binding_en, sa_scorer]
    
    def calc_reward_from_objective_values(values, conf):
        """ Must return a float based on results of objective functions (values) 
        """
        binding_en = values[0]
        #sa_score = values[1]
        print("score: ", binding_en)

        # Use base score as rough binding energy of adamantane 
        base_score = -10.0

        score_diff = binding_en - base_score

        return - score_diff * 0.1 / (1 + abs(score_diff) * 0.1)

def initialise_host(conf):
    """ A "pseudo constructor" for the class containing only static methods, will only work these out at the start
    """
    global _host_initialised, v, hostmol, host_en, final_complex_writer

    v = vina.Vina(verbosity=0)
    # Note host must me aligned with z-axis
    v.set_receptor(conf["host_pdbqt"])
    hostmol = Chem.MolFromMolFile(conf["host_sdf"],removeHs=False,sanitize=False)
    # If using multiple properties, turn this into a dictionary
    host_en = conf["host_en"]
    final_complex_writer = Chem.SDWriter(conf["complexsdfout"])
    
    _host_initialised = True

# Methods for reward calculations
def smi2sdf(mol, conf):
    confmol = process_smi(mol, conf["molgen_n_confs"],conf["molgen_rmsd_threshold"])
    return confmol

def fix_charge(mol):
    """ looks for atoms with mismatching explicit valences and charges them (for example N) """
    mol.UpdatePropertyCache(strict=False)
    for at in mol.GetAtoms():
        if at.GetAtomicNum() == 7 and at.GetExplicitValence()==4 and at.GetFormalCharge()==0:
            at.SetFormalCharge(1)
    Chem.SanitizeMol(mol)
    
def get_opt(mol, outfile, conf):
    """ Optimises the xTB energy and retrieves energy
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
    
def isExo(complexmol, hostmol, conf):
    """ POST DOCKING AND OPTIMISATION 
        Checks for exo complex, if exo sets a bad score
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
    guest_coords = np.array([complexmol.GetConformer().GetAtomPosition(atm.GetIdx()) for count, atm in enumerate(complexmol.GetAtoms()) if count >= hostmol.GetNumAtoms()])
    host_coords = np.array([complexmol.GetConformer().GetAtomPosition(atm.GetIdx()) for count, atm in enumerate(complexmol.GetAtoms()) if not count >= hostmol.GetNumAtoms()])

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

# Test
if __name__ == "__main__":
    # Manual config
    conf = {"host_pdbqt":"data/host_aligned.pdbqt",
            "host_sdf":"data/host_aligned.sdf",
            "host_en":-157089.5126460131,
            "molgen_n_confs":1, 
            "molgen_rmsd_threshold":0.35, 
            "vina_num_rotations":4, 
            "vina_num_translations":4, 
            "xtb_tempdir":"/home/spine/DProjects/DMCTS/DChemSim/xtbtmp",
            "complexsdfout":"data/complexes.sdf",
            "centroid_diff_threshold": 4,
            "cavity_atoms_threshold": 6
            }
    
    # Test molecule "Tropylium cation"
    testsmi = "C1=CC=C[CH+]C=C1"
    testmol = Chem.MolFromSmiles(testsmi)

    # Funcntion (almost) as it appears in code
    import time

    start = time.time()
    def _get_objective_values(mol, conf):
        return [f(mol) for f in CBDock_reward.get_objective_functions(conf)]
    end = time.time()
    print(CBDock_reward.calc_reward_from_objective_values(_get_objective_values(testmol, conf), conf))
    print("time:",end-start)
