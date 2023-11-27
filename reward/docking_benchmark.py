""" This file takes an example dataset 
    (This would usually be some of the data from the molecules used to train the RNN),
    docks the guest in the host, and runs an xTB optimisation calculating the binding energy.
    This effectively runs the chemistry simulation reward function, but without the RNN and MCTS.
"""
from rdkit import Chem
import pandas as pd

# Methods for reward calculations
from smi2sdf import process_smi
from docking import score_map_mmff94, score_map_vina, score_map_comb
from reward_utils import fix_charge, get_opt, is_small_cylinder, is_exo, get_property_mol
from sascorer import calculateScore

def binding_en(mol, dockfunc):
    """ Calculate values contributing to the reward arising from binding
    """

    # Create a null molecule out (a molecule with no atoms) to keep indexes in output structure files
    nullmol = Chem.MolFromSmiles("C")
    nullmol.SetDoubleProp("en", 20.0)

    # 1. Generate conformers
    print("STATUS - generating conformers")
    guestmol = process_smi(mol, conf["molgen_n_confs"],conf["molgen_rmsd_threshold"])
    fix_charge(guestmol)

    # If the guest is too large, set bad score
    try:
        is_small = is_small_cylinder(guestmol)
    except:
        print("END - bad conformer")
        nullmol.SetDoubleProp("en", 20.1)
        return get_property_mol(nullmol)
        
    if not is_small:
        nullmol.SetDoubleProp("en", 20.2)
        print("END = - guest too large")
        return get_property_mol(nullmol)
    
    # 2. Dock the best conformer
    # note that currently this takes the best result from docking. 
    # I may want to do a quick optimisation on the top few complexes if they are very close together, but have large difference in pose
    print("STATUS - docking")
    complexmols, scores = dockfunc(
        guestmol,
        hostmol,
        conf["vina_num_rotations"],
        conf["vina_num_translations"],
        host_pdbqt
        )
    
    # Complexmols are ordered by score, so check through until an exo complex (pre xtb optimisation) is found
    for i in range(complexmols.GetNumConformers()):
        exo = is_exo(complexmols, hostmol, conf, confId=i)
        if not exo:
            complexmol = Chem.Mol(complexmols)
            complexmol.RemoveAllConformers()
            complexmol.AddConformer(complexmols.GetConformer(i), assignId=True)
            break
        elif i == complexmols.GetNumConformers()-1:
            print("END - Exo complex")
            nullmol.SetDoubleProp("en", 20.3)
            return get_property_mol(nullmol)
    
    # 3. Calculate binding energy
    try:
        print("STATUS - optimising guest")
        optguestmol, guest_en = get_opt(guestmol, "guestopt.out", conf)
        print("STATUS - optimising complex")
        optcomplexmol, complex_en = get_opt(complexmol, "complexopt.out", conf)
    except:
        print("END - couldn't optimise")
        nullmol.SetDoubleProp("en", 20.4)
        return get_property_mol(nullmol)
    
    # If the result of xTB optimisation is exo, set bad score
    exo = is_exo(optcomplexmol, hostmol, conf)
    if exo:
        print("END - Exo complex")
        nullmol.SetDoubleProp("en", 20.5)
        return get_property_mol(nullmol)
    
    bind_en = complex_en - guest_en - host_en
    optcomplexmol.SetDoubleProp("en", bind_en)

    # For testing, add scores from binding poses
    try:
        for rank, i in enumerate(scores):
            optcomplexmol.SetDoubleProp("en_" + str(rank), i)
    except:
        print("NOTE - couldn't add pose scores to mol object")

    return get_property_mol(optcomplexmol)

def sa_scorer(mol):
    """ Calculate synthetic accessibility score
    """
    sa_score = calculateScore(mol)
    
    return sa_score

if __name__ == "__main__":

    conf = {"host_pdbqt":"data/host_aligned.pdbqt",
            "host_en":-157089.5126460131,
            "molgen_n_confs":1, 
            "molgen_rmsd_threshold":0.35, 
            "vina_num_rotations":4, 
            "vina_num_translations":4, 
            "centroid_diff_threshold": 4,
            "cavity_atoms_threshold": 6,
            "molecules_to_dock":"data/TestData/fewmols.smi"
            }

    global _host_initialised, hostmol, host_pdbqt, host_en

    host_pdbqt = conf["host_pdbqt"]
    hostmol = Chem.MolFromMolFile(conf["host_sdf"],removeHs=False,sanitize=False)

    confs_per_guest = conf["vina_num_rotations"] * conf["vina_num_translations"] + 1
    for i in range(confs_per_guest):
        newhost = Chem.Conformer(hostmol.GetConformer(0))
        hostmol.AddConformer(newhost, assignId=True)

    # If using multiple properties, turn this into a dictionary
    host_en = conf["host_en"]
    print("globvars set")

    df = pd.read_csv(conf["molecules_to_dock"], sep="\t", names=["name","smiles"], index_col=False)
    print(df)
    
    for count,smi in enumerate(df["smiles"]):
        mol = Chem.MolFromSmiles(smi)
        guestmol = process_smi(mol, 1, 0.35)

        print("STATUS - docking ", count)
        try:
            complexmolvina = binding_en(mol, score_map_vina)
            complexmolcomb = binding_en(mol, score_map_comb)
            sa_score = sa_scorer(mol)
        except:
            continue

        df.loc[count,"vina"] = complexmolvina
        df.loc[count,"v_en"] = complexmolvina.GetDoubleProp("en")
        df.loc[count,"comb"] = complexmolcomb
        df.loc[count,"c_en"] = complexmolcomb.GetDoubleProp("en")
        df.loc[count,"sa_score"] = sa_score

        df.to_pickle("data/TestData/fewmols.pkl")

