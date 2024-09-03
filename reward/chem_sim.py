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
    from reward.reward_utils import fix_charge, is_small_cylinder, is_exo, get_incorrect_bond_angle, get_incorrect_bond_length, covalent_CB, get_property_mol
    from reward.docking import DockLigand
    from reward.xtb_opt import xtbEnergy
    from reward.gaff_opt import AmberCalculator

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
        self.outdir = conf["output_dir"]
        
        # SET AS STANDALONE CHEMISTRY SIMULATOR
        self.is_standalone = kwargs.get("standalone", False)
        self.proplist = []
        self.calculator = AmberCalculator(conf)
        self.docking = DockLigand(conf)

    def setup(self):
        """ Sets up output directory for chemistry simulation data
            Maybe try the global design pattern for this later if ChemSim is only initialized once
        """
        if self.is_standalone:
            try:
                os.mkdir(self.outdir)
            except:
                print("WARNING - output directory already exists")

        self.bindingfile = os.path.join(self.outdir, "binding_sim.csv")
        self.complexfile = os.path.join(self.outdir, "complex_sim.pkl")
        self.guestfile = os.path.join(self.outdir, "guest_sim.pkl")
        self.complexdf = pd.DataFrame(columns=self.proplist)
        self.guestdf = pd.DataFrame(columns=self.proplist)
        self.bindingdf = pd.DataFrame(columns=self.proplist)

        # Set the paths of the host_pdbqt and host_sdf to absolute paths inside the scope of the object
        self.conf["host_pdbqt"] = os.path.realpath(self.conf["host_pdbqt"])
        self.conf["host_sdf"] = os.path.realpath(self.conf["host_sdf"])
        self.conf["output_dir"] = os.path.realpath(self.conf["output_dir"])

    def flush(self, propertymol, guest=False):
        """ Writes the output of a molecule to the output file
            Properties are written to the propertymol based on the config
        """
        moldict = {}
        bindingdict = {}
        for i in propertymol.GetPropNames():
            if i.split("_")[0] == "binding":
                try:
                    bindingdict[i] = float(propertymol.GetProp(i))
                except:
                    bindingdict[i] = propertymol.GetProp(i)
            else:
                try:
                    moldict[i] = float(propertymol.GetProp(i))
                except:
                    moldict[i] = propertymol.GetProp(i)

        # Add smiles into bindingdict
        bindingdict["smiles"] = propertymol.GetProp("smiles")
            
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
            idx = len(self.complexdf) + 1
            df_mol = pd.DataFrame(moldict, index=[idx])
            self.complexdf = pd.concat([self.complexdf, df_mol], axis=0)
            self.complexdf.to_pickle(self.complexfile)

            idx = len(self.bindingdf) + 1
            df_mol = pd.DataFrame(bindingdict, index=[idx])
            self.bindingdf = pd.concat([self.bindingdf, df_mol], axis=0)
            self.bindingdf.to_csv(self.bindingfile)
            

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
        try:
            complexmol = self.docking.get_docked_complex(guestmol)
        except Exception as e:
            raise ChemSimError("Error in docking") from e
                    
        return complexmol, guestmol

    def run_opt(self, complexmol, guestmol):
        # 3. Calculate binding energy
        # Definitely set this up to reuse this object otherwise memory usage will get out of hand
        try:
            optcomplexmol, optguestmol = self.calculator.get_guest_complex_opt(complexmol, guestmol)
        except Exception as e:
            raise ChemSimError("Error in optimisation") from e
        
        # POST FILTERS (sets bad score in MCTS, sets flags in ChemSim)
        # CURRENTLY ONLY APPLY TO HCs DIFFERENT POST FITLERS NEEDED FOR OTHER TYPES OF MOLECULES
        
        # If the result of xTB optimisation is exo, set bad score
        """
        exo = is_exo(optcomplexmol, self.hostmol, self.conf)
        if exo:
            optcomplexmol.SetProp("is_exo", "True")

        covalent = covalent_CB(optcomplexmol)
        if covalent:
            optcomplexmol.SetProp("covalent_CB", "True")

        bad_angle = get_incorrect_bond_angle(optcomplexmol)
        if bad_angle:
            optcomplexmol.SetProp("bad_angle", "True")

        bad_length = get_incorrect_bond_length(optcomplexmol)
        if bad_length:
            optcomplexmol.SetProp("bad_length", "True")
        """

        # Get decomposition of binding energies:
        complex_en_dict = {propname : optcomplexmol.GetProp(propname) for propname in optcomplexmol.GetPropNames()}
        guest_en_dict = {propname : optguestmol.GetProp(propname) for propname in optguestmol.GetPropNames()}

        binding_en_dict = {f"binding_{key_comp}" : float(val_comp) - float(val_host) - float(val_guest) for (key_comp, val_comp), (key_host, val_host), (key_guest, val_guest) in zip(complex_en_dict.items(), self.calculator.host_en_dict.items(), guest_en_dict.items())}

        for key, val in binding_en_dict.items():
            optcomplexmol.SetDoubleProp(key, val)

        # Set en for MCTS
        if self.conf["thermo"]:
            en = [val for key, val in binding_en_dict.items() if key.split("_")[-1] == "Ecorr" ][0]
        else:
            en = [val for key, val in binding_en_dict.items() if key.split("_")[-1] == "en" ][0]
            
        optcomplexmol.SetDoubleProp("en", en)

        return get_property_mol(optcomplexmol), get_property_mol(optguestmol)
    
if __name__ == "__main__":
    import os
    import argparse
    import yaml

    from smi2sdf import process_smi
    from docking import DockLigand
    from reward_utils import fix_charge, is_small_cylinder, is_exo, get_incorrect_bond_angle, get_incorrect_bond_length, covalent_CB, get_property_mol
    from xtb_opt import xtbEnergy
    from gaff_opt import AmberCalculator
    
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["OMP_NUM_THREADS"] = "8"

    print("Starting ChemSim binding energy evaluation")

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

    hostmol = Chem.MolFromMolFile(conf["host_sdf"], removeHs=False)

    """
    confs_per_guest = conf["vina_num_rotations"] * conf["vina_num_translations"] + 1
    for i in range(confs_per_guest):
        newhost = Chem.Conformer(hostmol.GetConformer(0))
        hostmol.AddConformer(newhost, assignId=True)
    """

    print("host initialised")

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
            if conf["docking_n_confs"] > 1:
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

            else:
                bestconf, bestguestconf = simulator.run_dock(confs[0])
                
            molout, guestmolout = simulator.run_opt(bestconf, bestguestconf)
            molout.SetProp("smiles", smi)
            guestmolout.SetProp("smiles", smi)

            return molout, guestmolout

        except Exception as e:
            print("Problem with smiles: ", smi)
            print(e)
            print(traceback.format_exc())
            return None
        
    with Parallel(n_jobs=conf["n_jobs"], prefer="processes", return_as="generator", verbose=51) as parallel:
        for result in parallel(delayed(process_molecule_wrapper)(simulator, smi) for smi in smi_gen):
            if result:
                simulator.flush(result[0])
                simulator.flush(result[1], guest=True)
