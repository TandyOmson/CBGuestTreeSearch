import vina
from rdkit import Chem

# Methods for reward calculations
from reward.smi2sdf import process_smi
from reward.docking import score_map_mmff94, score_map_vina
from reward.reward_utils import fix_charge, get_opt, is_small_cylinder, is_exo, get_property_mol
from reward.sascorer import calculateScore

# Can inherit from the batchReward abstract class
from reward.reward import Reward

# On module call, _host_initialised is set to False, then becomes true on first call to binding_en
_host_initialised = False

class CBDock_reward(Reward):
    """ Reward Class in config. A reference to this class from config is passed through the MCTS class to the evaluate_node method in utils 
        The following methods are static (no instantiation of class required) and abstract (must be included in the class) methods.
    """
    def get_objective_functions(conf):
        """ Must return a list of functions """

        def binding_en(mol):
            """ Calculate values contributing to the reward arising from binding
            """

            if not _host_initialised:
                _initialise_host(conf)

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
            complexmols, scores = score_map_mmff94(
                guestmol,
                hostmol,
                conf["vina_num_rotations"],
                conf["vina_num_translations"]
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
        
        def binding_en_VINA(mol):
            """ Calculate values contributing to the reward arising from binding
            """

            if not _host_initialised:
                _initialise_host(conf)

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
                return get_property_mol(nullmol)
                
            if not is_small:
                print("END = - guest too large")
                return get_property_mol(nullmol)
            
            # 2. Dock the best conformer
            # note that currently this takes the best result from docking. 
            # I may want to do a quick optimisation on the top few complexes if they are very close together, but have large difference in pose
            print("STATUS - docking")
            complexmols, scores = score_map_vina(
                guestmol,
                hostmol,
                conf["vina_num_rotations"],
                conf["vina_num_translations"],
                v
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
                    return get_property_mol(nullmol)
            
            # 3. Calculate binding energy
            try:
                print("STATUS - optimising guest")
                optguestmol, guest_en = get_opt(guestmol, "guestopt.out", conf)
                print("STATUS - optimising complex")
                optcomplexmol, complex_en = get_opt(complexmol, "complexopt.out", conf)
            except:
                print("END - couldn't optimise")
                return get_property_mol(nullmol)
            
            # If the result of xTB optimisation is exo, set bad score
            exo = is_exo(optcomplexmol, hostmol, conf)
            if exo:
                print("END - Exo complex")
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
            
        return [binding_en, sa_scorer, binding_en_VINA]
    
    def calc_reward_from_objective_values(values, conf):
        """ Must return a float based on results of objective functions (values) 
        """
        binding_en = float(values[0].GetProp("en"))
        #sa_score = values[1]
        print("score: ", binding_en)

        # Use base score as rough binding energy of adamantane 
        base_score = -10.0

        score_diff = binding_en - base_score

        return - score_diff * 0.1 / (1 + abs(score_diff) * 0.1)

def _initialise_host(conf):
    """ A "pseudo constructor" for the class containing only static methods, will only work these out at the start
    """
    global _host_initialised, v, hostmol, host_en
    #global guestwriter, complexwriter, posewriter

    v = vina.Vina(verbosity=0)
    # Note host must me aligned with z-axis
    v.set_receptor(conf["host_pdbqt"])
    hostmol = Chem.MolFromMolFile(conf["host_sdf"],removeHs=False,sanitize=False)

    confs_per_guest = conf["vina_num_rotations"] * conf["vina_num_translations"] + 1
    for i in range(confs_per_guest):
        newhost = Chem.Conformer(hostmol.GetConformer(0))
        hostmol.AddConformer(newhost, assignId=True)
    
    # If using multiple properties, turn this into a dictionary
    host_en = conf["host_en"]

    #guestwriter = Chem.SDWriter(conf["molsoutdir"] + "/guests.sdf")
    #complexwriter = Chem.SDWriter(conf["molsoutdir"] + "/complexes.sdf")
    #posewriter = Chem.SDWriter(conf["molsoutdir"] + "/poses.sdf")

    print("globvars set")
    
    _host_initialised = True

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
            "centroid_diff_threshold": 4,
            "cavity_atoms_threshold": 6,
            }
    
    # Test molecule "Tropylium cation"
    testsmi = "C1=CC=C[CH+]C=C1"
    testmol = Chem.MolFromSmiles(testsmi)

    testsmi2 = "C(CCCCC#C)(=O)NCCC=1C=CC(=C(C=1)O)O"
    testmol2 = Chem.MolFromSmiles(testsmi2)

    # _get_objective_values (almost) as it appears in code
    import time

    start = time.time()
    def _get_objective_values(mol, conf):
        return [f(mol) for f in CBDock_reward.get_objective_functions(conf)]
    end = time.time()
    print(CBDock_reward.calc_reward_from_objective_values(_get_objective_values(testmol2, conf), conf))
    print(CBDock_reward.calc_reward_from_objective_values(_get_objective_values(testmol, conf), conf))
    print("time:",end-start)