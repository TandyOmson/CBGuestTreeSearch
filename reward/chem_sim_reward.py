from rdkit import Chem
import traceback
import numpy as np

# Methods for reward calculations
from reward.chem_sim import ChemSim
from reward.sascorer import calculateScore

# Can inherit from the batchReward abstract class
from reward.reward import Reward

# On module call, _host_initialised is set to False, then becomes true on first call to binding_mol
_host_initialised = False

class CBDock_reward(Reward):
    """ Reward Class in config. A reference to this class from config is passed through the MCTS class to the evaluate_node method in utils 
        The following methods are static (no instantiation of class required) and abstract (must be included in the class) methods.
    """
    def get_objective_functions(conf):
        """ Must return a list of functions """

        def binding_mol(mol):
            """ Calculate values contributing to the reward arising from binding
            """
            if not _host_initialised:
                _initialise_host(conf)

            simulator = ChemSim(conf, hostmol)
            simulator.setup()
            smi = mol.GetProp("_Smiles")

            try:
                confs = simulator.get_confs(mol)
                
                # Try top conformers (number determined by conf["molgen_n_confs"])
                # Set optlevel
                if conf["molgen_n_confs"] > 1:
                    org_optlevel = conf["optlevel"]
                    org_thermo = conf["thermo"]
                    conf["optlevel"] = "loose"
                    conf["thermo"] = False

                    molsout = []
                    guestmolsout = []

                    molsoutdock = []
                    guestmolsoutdock = []
                    for i in confs:
                        confmoldock, confguestmoldock = simulator.run_dock(i)
                        confmolout, confguestmolout = simulator.run_opt(confmoldock, confguestmoldock)

                        molsoutdock.append(confmoldock)
                        guestmolsoutdock.append(confguestmoldock)

                        molsout.append(confmolout)
                        guestmolsout.append(confguestmolout)

                    # Identify the best binding energy from crude optimisation
                    bind_ens = [float(i.GetDoubleProp("en")) for i in molsout]
                    best_idx = np.argmin(bind_ens)
                    bestconf, bestguestconf = molsoutdock[best_idx], guestmolsoutdock[best_idx]

                    conf["optlevel"] = org_optlevel
                    conf["thermo"] = org_thermo
                
                molout, guestmolout = simulator.run_opt(bestconf, bestguestconf)
                molout.SetProp("smiles", smi)
                guestmolout.SetProp("smiles", smi)

                # Record additional properties from ChemSim (only complex included in reward functions)
                simulator.flush_csv_no_mol(molout)
                simulator.flush_csv_no_mol(guestmolout, guest=True)

                return molout

            except Exception as e:
                print("Problem with smiles: ", smi)
                print(e)
                print(traceback.format_exc())

                # Still have to flush to keep SMILES in MCTS results corresponding to ChemSim results
                nullmol = Chem.MolFromSmiles("C")
                nullmol.SetDoubleProp("en", 20.0)
                nullmol.SetProp("smiles", smi)
                simulator.flush_csv_no_mol(nullmol)
                simulator.flush_csv_no_mol(nullmol, guest=True)
                
                return None

        def sa_scorer(mol):
            """ Calculate synthetic accessibility score
            """
            sa_score = calculateScore(mol)
            
            return sa_score
            
        return [binding_mol, sa_scorer]
    
    def calc_reward_from_objective_values(values, conf):
        """ Must return a float based on results of objective functions (values) 
        """
        if values[0] is None:
            return -1.0

        binding_en = float(values[0].GetProp("en"))
        #sa_score = values[1]
        print("binding: ", binding_en)

        # Use base score as rough binding energy of adamantane 
        base_score = -20.0
        score_diff = binding_en - base_score
        
        # + sa_score/5.0
        return - score_diff * 0.1 / (1 + abs(score_diff) * 0.1)

def _initialise_host(conf):
    """ A "pseudo constructor" for the class containing only static methods, will only work these out at the start
    """
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
    
    _host_initialised = True

# Test
if __name__ == "__main__":
    import os
    os.environ["OPENBLAS_NUM_THREADS"] = "1"

    # Manual config
    conf = {"host_pdbqt":"data/host_aligned.pdbqt",
            "host_sdf":"data/host_aligned.sdf",
            "host_en":-157089.5126460131,
            "molgen_n_confs":1, 
            "molgen_rmsd_threshold":0.35, 
            "vina_num_rotations":4, 
            "vina_num_translations":4, 
            "centroid_diff_threshold": 4,
            "cavity_atoms_threshold": 6,
            "thermo":True,
            "solvent":True,
            "optlevel":"normal",
            "additional_ens":False,
            "partial_charges":False
            }
    
    # Test molecule "Tropylium cation"
    testsmi = "CCC"
    testmol = Chem.MolFromSmiles(testsmi)

    testsmi2 = "C(CCCCC#C)(=O)NCCC=1C=CC(=C(C=1)O)O"
    testmol2 = Chem.MolFromSmiles(testsmi2)

    # _get_objective_values (almost) as it appears in code
    def _get_objective_values(mol, conf):
        return [f(mol) for f in CBDock_reward.get_objective_functions(conf)]
    print(CBDock_reward.calc_reward_from_objective_values(_get_objective_values(testmol, conf), conf))
    #print(CBDock_reward.calc_reward_from_objective_values(_get_objective_values(testmol2, conf), conf))
