from rdkit import Chem
from rdkit.Chem import PropertyMol
import traceback
import numpy as np
import os
import pickle

# Methods for reward calculations
from reward.chem_sim import ChemSim
from reward.sascorer import calculateScore

# Can inherit from the batchReward abstract class
from reward.reward import Reward

class CBDock_reward(Reward):
    """ Reward Class in config. A reference to this class from config is passed through the MCTS class to the evaluate_node method in utils 
        The following methods are static (no instantiation of class required) and abstract (must be included in the class) methods.
    """
    def get_objective_functions(conf, simulator):
        """ Must return a list of functions """

        def binding_mol(mol):
            """ Calculate values contributing to the reward arising from binding
            """
                
            smi = mol.GetProp("_Smiles")
            print("running", smi)

            try:
                confs = simulator.get_confs(mol)
                
                bestconf, bestguestconf = simulator.run_dock(confs[0])
                molout, guestmolout = simulator.run_opt(bestconf, bestguestconf)
                
                molout.SetProp("smiles", smi)
                guestmolout.SetProp("smiles", smi)

                # If post filters were triggered, set a bad score
                """
                if molout.HasProp("isExo") == 1:
                    molout.SetDoubleProp("en", 25.0)

                if molout.HasProp("covalent_CB") == 1:
                    molout.SetDoubleProp("en", 25.0)

                if molout.HasProp("bad_angle") == 1:
                    molout.SetDoubleProp("en", 25.0)

                if molout.HasProp("bad_length") == 1:
                    molout.SetDoubleProp("en", 25.0)
                """

                # Required as flush wipes mol properties
                #en = float(molout.GetProp("en"))
                binding_mol = PropertyMol.PropertyMol(molout)
                #binding_mol.SetProp("en", en)

                guestmolout = PropertyMol.PropertyMol(guestmolout)

                # Record additional properties from ChemSim (only complex included in reward functions)
                #simulator.flush(molout)
                #simulator.flush(guestmolout, guest=True)

                # Return pickled object so that properties persist
                return pickle.dumps(binding_mol), pickle.dumps(guestmolout)

            except Exception as e:
                print("Problem with smiles: ", smi)
                print(e)
                print(traceback.format_exc())

                with open(f"{os.path.abspath(conf['output_dir'])}/failed.smi", "a") as fa:
                    fa.write(f"{smi}\n")
                
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
        if not values[0]:
            return -1.0
        else:
            try:
                binding_en = float(pickle.loads(values[0][0]).GetProp("en"))
            except:
                return -1.0
            
        #sa_score = values[1]
        print("binding: ", binding_en)

        # Use base score as rough binding energy of adamantane 
        base_score = -20.0
        score_diff = binding_en - base_score
        
        # + sa_score/5.0
        reward_score = - score_diff * 0.1 / (1 + abs(score_diff) * 0.1)
        print("reward score:", reward_score)

        if np.isnan(reward_score):
            print("reward returned nan")
            return -1.0
        
        return reward_score

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
