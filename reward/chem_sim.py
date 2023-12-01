""" Module for storing the modular chemistry simulator either applied in chem_sim_reward.py
    or used standalone for data collection and benchmarking
    Three key phases:
    - Molgen
    - Docking
    - Binding energy calculation

    The chemistry simulator logs data separately to the tree search (which only needs outputs relevant to the reward function)
    However config for standalone running, and ChemSim logging are found in the same config file as the MCTS
    Usage:
    (chemts) $ python3 chem_sim.py -c <mcts_config>
"""
from rdkit import Chem

# Methods for binding calculations
from reward.smi2sdf import process_smi
from reward.docking import score_map_vina, score_map_comb
from reward.reward_utils import fix_charge, is_small_cylinder, is_exo, get_property_mol
from reward.xtb_opt import xtbEnergy

class ChemSim():
    """ Contains the chemitry simulator
    """
    def __init__(self, conf, hostmol):
        self.conf = conf
        # Non optional variables set as class attributes
        self.hostmol = hostmol
        self.host_en = conf["host_en"]

    # ... output setup functions

    def run(self, mol):
        """ Runs the chemistry simulator based on the setup for one mol
            The relevant output is the return property mol containing binding energy and complex geoemtry
            All other ouput is handled by the ChemSim class methods
        """

        # Create a null molecule out (a molecule with no atoms) to keep indexes in output structure files
        nullmol = Chem.MolFromSmiles("C")
        nullmol.SetDoubleProp("en", 20.0)

        # 1. Generate conformers
        print("STATUS - generating conformers")
        guestmol = process_smi(mol, self.conf["molgen_n_confs"], self.conf["molgen_rmsd_threshold"])
        fix_charge(guestmol)

        # If the guest is too large, set bad score
        try:
            is_small = is_small_cylinder(guestmol)
        except:
            print("END - bad conformer")
            guestmol.SetDoubleProp("en", 20.1)
            return get_property_mol(nullmol)
            
        if not is_small:
            guestmol.SetDoubleProp("en", 20.2)
            print("END = - guest too large")
            return get_property_mol(guestmol)
        
        # 2. Dock the best conformer
        # note that currently this takes the best result from docking. 
        # I may want to do a quick optimisation on the top few complexes if they are very close together, but have large difference in pose
        print("STATUS - docking")
        try:
            complexmols, scores = score_map_vina(
                guestmol,
                self.hostmol,
                self.conf["vina_num_rotations"],
                self.conf["vina_num_translations"],
                self.conf["host_pdbqt"]
                )
        except:
            complexmols, scores = score_map_comb(
                guestmol,
                self.hostmol,
                self.conf["vina_num_rotations"],
                self.conf["vina_num_translations"],
                self.conf["host_pdbqt"]
                ) 
        
        # Complexmols are ordered by score, so check through until an exo complex (pre xtb optimisation) is found
        for i in range(complexmols.GetNumConformers()):
            exo = is_exo(complexmols, self.hostmol, self.conf, confId=i)
            if not exo:
                complexmol = Chem.Mol(complexmols)
                complexmol.RemoveAllConformers()
                complexmol.AddConformer(complexmols.GetConformer(i), assignId=True)
                break
            elif i == complexmols.GetNumConformers()-1:
                print("END - Exo complex")
                guestmol.SetDoubleProp("en", 20.3)
                return get_property_mol(guestmol)
        
        # 3. Calculate binding energy
        # Definitely set this up to reuse this object otherwise memory usage will get out of hand
        calc = xtbEnergy(self.conf)
        # try:
        print("STATUS - optimising guest")
        optguestmol, guest_en = calc.get_opt(guestmol)
        print("STATUS - optimising complex")
        optcomplexmol, complex_en = calc.get_opt(complexmol)
        # except:
        #     print("END - couldn't optimise")
        #     nullmol.SetDoubleProp("en", 20.4)
        #     return get_property_mol(nullmol)
        
        # If the result of xTB optimisation is exo, set bad score
        exo = is_exo(optcomplexmol, self.hostmol, self.conf)
        if exo:
            print("END - Exo complex")
            guestmol.SetDoubleProp("en", 20.5)
            return get_property_mol(guestmol)
        
        bind_en = complex_en - guest_en - self.host_en
        optcomplexmol.SetDoubleProp("en", bind_en)

        # For docking testing, add scores from binding poses
        try:
            for rank, i in enumerate(scores):
                optcomplexmol.SetDoubleProp("en_" + str(rank), i)
        except:
            print("NOTE - couldn't add pose scores to mol object")

        # Additional optional outputs from binding energy calculation

        # Additional outputs

        return get_property_mol(optcomplexmol)