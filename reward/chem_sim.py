""" Module for storing the modular chemistry simulator either applied in chem_sim_reward.py
    or used standalone for data collection and benchmarking
    Three key phases:
    - Molgen
    - Docking
    - Binding energy calculation

    The chemistry simulator logs data separately to the tree search (which only needs outputs relevant to the reward function)
    However config for standalone running, and ChemSim logging are found in the same config file as the MCTS
    Usage:
    (chemts) $ python3 chem_sim.py -i <.smi file> -c <mcts_config>
"""
from rdkit import Chem
import pandas as pd
import os
from joblib import Parallel, delayed
import traceback
import numpy as np

import traceback

import time

# Methods for binding calculations
if __name__ != "__main__":
    from reward.smi2sdf import process_smi
    from reward.reward_utils import fix_charge, is_small_cylinder, is_exo, check_smiles_change, get_property_mol
    from reward.docking import DockLigand
    from reward.xtb_opt import xtbEnergy

# Class for diagnostic error handling
class ChemSimError(Exception):
    def __init__(self, message, original_exception=None):
        self.message = message
        self.original_exception = original_exception
        super().__init__(self.message)

    def __str__(self):
        if self.original_exception:
            return f"{self.message}. Original exception: {type(self.original_exception).__name__}: {str(self.original_exception)}"
        else:
            return self.message

class ChemSim():
    """ Contains the chemitry simulator
    """
    def __init__(self, conf, hostmol, **kwargs):
        self.conf = conf
        # Non optional variables set as class attributes
        self.hostmol = hostmol
        self.host_en = conf["host_en"]
        self.outdir = conf["output_dir"]
        # SET AS STANDALONE CHEMISTRY SIMULATOR
        self.is_standalone = kwargs.get("standalone", False)
        self.proplist = []

    def setup(self):
        """ Sets up output directory for chemistry simulation data
            Maybe try the global design pattern for this later if ChemSim is only initialized once
        """
        if self.is_standalone:
            try:
                os.mkdir(self.outdir)
            except:
                print("WARNING - output directory already exists")
        
        # if standalone, make it a pickle (to include structures, else just a csv)
        if self.is_standalone:
            self.outfile = os.path.join(self.outdir, "chem_sim.pkl")
            self.guestfile = os.path.join(self.outdir, "guest_sim.pkl")
            self.df = pd.DataFrame(columns=self.proplist.append("binding_mol"))
            self.guestdf = pd.DataFrame(columns=self.proplist.append("binding_mol"))
        else:
            self.outfile = os.path.join(self.outdir, "chem_sim.csv")
            self.guestfile = os.path.join(self.outdir, "guest_sim.csv")
            self.df = pd.DataFrame(columns=self.proplist)
            self.guestdf = pd.DataFrame(columns=self.proplist)

    def flush(self, propertymol, guest=False):
        """ Writes the output of a molecule to the output file
            Properties are written to the propertymol based on the config
        """
        moldict = {}
        for i in propertymol.GetPropNames():
            try:
                moldict[i] = float(propertymol.GetProp(i))
            except:
                moldict[i] = propertymol.GetProp(i)
            
        # Clear properties in molecule object to save space
        for i in propertymol.GetPropNames():
            propertymol.ClearProp(i)
        # Retain molecule structure
        moldict["mol"] = propertymol

        if guest:
            idx = len(self.guestdf) + 1
            df_mol = pd.DataFrame(moldict, index=[idx])
            self.guestdf = pd.concat([self.guestdf, df_mol], axis=0)
            self.guestdf.to_pickle(self.guestfile)
        else:
            idx = len(self.df) + 1
            df_mol = pd.DataFrame(moldict, index=[idx])
            self.df = pd.concat([self.df, df_mol], axis=0)
            self.df.to_pickle(self.outfile)

    def flush_csv_no_mol(self, propertymol, guest=False):
        """ Writes the output of a molecule to the output file - for additional energies
            Mol structure is NOT saved, as this is in the MCTS .pkl
        """
        moldict = {}
        for i in propertymol.GetPropNames():
            try:
                moldict[i] = float(propertymol.GetProp(i))
            except:
                moldict[i] = propertymol.GetProp(i)
            
        # Clear properties in molecule object to save space
        for i in propertymol.GetPropNames():
            propertymol.ClearProp(i)
        # Retain molecule structure
        moldict["mol"] = propertymol
        
        if guest:
            num = len(self.guestdf) + 1
            df_mol = pd.DataFrame(moldict, index=[num])
            self.guestdf = pd.concat([self.guestdf, df_mol], axis=0)
            self.guestdf.to_csv(self.guestfile)
        else:
            num = len(self.df) + 1
            df_mol = pd.DataFrame(moldict, index=[num])
            self.df = pd.concat([self.df, df_mol], axis=0)
            self.df.to_csv(self.outfile)

    def get_confs(self, mol):
        """ Generates conformers of a molecule
            The molgen_n_confs parameter determines how many per molecule to generated
        """
        try:
            guestmol = process_smi(mol, self.conf["molgen_n_confs"], self.conf["molgen_rmsd_threshold"])
        except Exception as e:
            raise ChemSimError("Error in conformer generation") from e
        
        fix_charge(guestmol)

        confs_as_mols = []
        for i in guestmol.GetConformers():
            conf = Chem.Mol(guestmol)
            conf.RemoveAllConformers()
            conf.AddConformer(i, assignId=True)
            confs_as_mols.append(conf)

        return confs_as_mols

    def run_dock(self, guestmol):
        """ Runs the chemistry simulator based on the setup for one mol
            The relevant output is the return property mol containing binding energy and complex geoemtry
            All other ouput is handled by the ChemSim class methods
        """
        # If the guest is too large, set bad score
        try:
            is_small = is_small_cylinder(guestmol)
        except:
            # If there is an error in the check, assume it is small
            is_small = True

        # 2. Dock the best conformer
        #print("STATUS - docking")
        dock = DockLigand(self.hostmol, self.conf)
        
        if not is_small:
            if self.conf["vina_large_guest"]:
                try:
                    complexmols, scores = dock.vina_dock(guestmol)
                except Exception as e:
                    raise ChemSimError("Error in vina docking of large guest") from e

            else:
                raise ChemSimError("Guest too large, vina_large_guest parameter set to False")
        
        # 3 attempts at docking for small molecules
        # if is_small:
        #     try:
        #         complexmols, scores = dock.score_map_comb(guestmol)
        #     except:
        #         try:
        #             complexmols, scores = dock.score_map_vina(guestmol)
        #         except:
        #             try:
        #                 complexmols, scores = dock.vina_dock(guestmol)
        #             except Exception as e:
        #                 raise ChemSimError("Error in docking of small guest") from e
        complexmols = []
        scores = []
        docking_failures = 0
        if is_small:
            try:
                start = time.time()
                complexmol_dock, score = dock.score_map_comb(guestmol)
                complexmols.append(complexmol_dock)
                scores.extend(score)
                end = time.time()
                print("comb", end-start)
            except:
                docking_failures += 1
            try:
                start = time.time()
                complexmol_dock, score = dock.vina_dock(guestmol)
                complexmols.append(complexmol_dock)
                scores.extend(score)
                end = time.time()
                print("vina", end-start)
            except:
                docking_failures += 1
            try:
                start = time.time()
                complexmol_dock, score = dock.score_map_vina(guestmol)
                complexmols.append(complexmol_dock)
                scores.extend(score)
                end = time.time()
                print("vina_score", end-start)
            except:
                docking_failures += 1

            if docking_failures == 3:
                raise ValueError("All docking methods failed")

        complexmol_confs = Chem.Mol(complexmols[0])
        complexmol_confs.RemoveAllConformers()       
        for i in complexmols:
            for c in i.GetConformers():
                complexmol_confs.AddConformer(c, assignId=True)

        try:
            complexmol = dock.get_best_pose(complexmol_confs, scores)
        except ValueError as e:
            raise ChemSimError("Error in getting best pose") from e

        return complexmol, guestmol

    def run_opt(self, complexmol, guestmol):
        # 3. Calculate binding energy
        # Definitely set this up to reuse this object otherwise memory usage will get out of hand
        calc = xtbEnergy(self.conf)
        try:
            #print("STATUS - optimising guest")
            optguestmol, guest_en = calc.get_opt(guestmol)
            #print("STATUS - optimising complex")
            optcomplexmol, complex_en = calc.get_opt(complexmol)
        except Exception as e:
            raise ChemSimError("Error in xTB optimisation") from e
        
        # If the result of xTB optimisation is exo, set bad score
        exo = is_exo(optcomplexmol, self.hostmol, self.conf)
        if exo:
            optcomplexmol.SetProp("is_exo", "True")
        
        # # Check if the smiles of the guest have changed
        # changed = check_smiles_change(optguestmol, Chem.GetMolFrags(optcomplexmol, asMols=True)[1])
        # if changed:
        #     optcomplexmol.SetProp("is_changed", "True")
        
        bind_en = complex_en - guest_en - self.host_en
        optcomplexmol.SetDoubleProp("en", bind_en)
        optguestmol.SetDoubleProp("en", guest_en)

        return get_property_mol(optcomplexmol), get_property_mol(optguestmol)
    
if __name__ == "__main__":
    import os
    import argparse
    import yaml

    from smi2sdf import process_smi
    from docking import DockLigand
    from reward_utils import fix_charge, is_small_cylinder, is_exo, check_smiles_change, get_property_mol
    from xtb_opt import xtbEnergy
    
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["OMP_NUM_THREADS"] = "16"

    parser = argparse.ArgumentParser(
        description="",
        usage=f"python {os.path.basename(__file__)} -c <config_file> -i <input_file>"
    )
    parser.add_argument(
        "-c", "--config", type=str, required=True,
        help="path to a .yaml file"
    )
    parser.add_argument(
        "-i", "--input", type=str, required=True,
        help="path to a .smi file"
    )
    args = parser.parse_args()
    with open(args.config, 'r') as f:
        conf = yaml.load(f, Loader=yaml.SafeLoader)

    global _host_initialised, hostmol, host_pdbqt, host_en

    host_pdbqt = conf["host_pdbqt"]
    hostmol = Chem.MolFromMolFile(conf["host_sdf"],removeHs=False,sanitize=False)

    confs_per_guest = conf["vina_num_rotations"] * conf["vina_num_translations"] + 1
    for i in range(confs_per_guest):
        newhost = Chem.Conformer(hostmol.GetConformer(0))
        hostmol.AddConformer(newhost, assignId=True)
    
    host_en = conf["host_en"]
    print("globvars set")

    df = pd.read_csv(args.input, sep=" ", names=["smiles"], index_col=False)
    smi_gen = (i for i in df["smiles"])

    simulator = ChemSim(conf, hostmol, standalone=True)
    simulator.setup()

    # Wrapper method with callback for parallel processing
    def process_molecule_wrapper(simulator, smi):
        try:
            mol = Chem.MolFromSmiles(smi)
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

            return molout, guestmolout

        except Exception as e:
            print("Problem with smiles: ", smi)
            print(e)
            print(traceback.format_exc())
            return None
        
    with Parallel(n_jobs=8, prefer="processes", return_as="generator", verbose=51) as parallel:
        for result in parallel(delayed(process_molecule_wrapper)(simulator, smi) for smi in smi_gen):
            if result:
                simulator.flush(result[0])
                simulator.flush(result[1], guest=True)
