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

# Methods for binding calculations
if __name__ != "__main__":
    from reward.smi2sdf import process_smi
    from reward.reward_utils import fix_charge, is_small_cylinder, is_exo, check_smiles_change, get_property_mol
    from reward.docking import DockLigand
    from reward.xtb_opt import xtbEnergy

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
        if self.is_standalone:
            # Convert the propertymol to a dictionary
            moldict = {}
            for i in propertymol.GetPropNames():
                try:
                    moldict[i] = float(propertymol.GetProp(i))
                except:
                    moldict[i] = propertymol.GetProp(i)
            # Remove the props from the mol to save space
            for i in propertymol.GetPropNames():
                propertymol.ClearProp(i)
            # Add the mol back to get the complexed structure
            moldict["mol"] = propertymol
            num = len(self.df) + 1
            df_mol = pd.DataFrame(moldict, index=[num])
            if guest:
                self.guestdf = pd.concat([self.guestdf, df_mol], axis=0)
                self.guestdf.to_pickle(self.guestfile)
            else:
                self.df = pd.concat([self.df, df_mol], axis=0)
                self.df.to_pickle(self.outfile)
        else:
            # Convert the propertymol to a dictionary
            moldict = {}
            for i in propertymol.GetPropNames():
                try:
                    moldict[i] = float(propertymol.GetProp(i))
                except:
                    moldict[i] = propertymol.GetProp(i)
            num = len(self.df) + 1
            df_mol = pd.DataFrame(moldict, index=[num])
            if guest:
                self.guestdf = pd.concat([self.guestdf, df_mol], axis=0)
                self.guestdf.to_csv(self.guestfile)
            else:
                self.df = pd.concat([self.df, df_mol], axis=0)
                self.df.to_csv(self.outfile)

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
        try:
            guestmol = process_smi(mol, self.conf["molgen_n_confs"], self.conf["molgen_rmsd_threshold"])
        except:
            print("END - bad conformer")
            nullmol.SetDoubleProp("en", 20.1)
            return get_property_mol(nullmol), get_property_mol(nullmol)
        fix_charge(guestmol)

        # If the guest is too large, set bad score
        try:
            is_small = is_small_cylinder(guestmol)
        except:
            print("END - bad conformer")
            guestmol.SetDoubleProp("en", 20.1)
            return get_property_mol(nullmol), get_property_mol(nullmol)

        # 2. Dock the best conformer
        print("STATUS - docking")
        dock = DockLigand(self.hostmol, self.conf)
        
        if not is_small:
            if self.conf["vina_large_guest"]:
                try:
                    complexmols, scores = dock.vina_dock(guestmol)
                except:
                    guestmol.SetDoubleProp("en", 20.2)
                    print("END = - guest too large, couldnt dock by vina")
                    return get_property_mol(guestmol), get_property_mol(nullmol)

            else:
                guestmol.SetDoubleProp("en", 20.2)
                print("END = - guest too large")
                return get_property_mol(guestmol), get_property_mol(nullmol)
        
        if is_small:
            try:
                complexmols, scores = dock.score_map_vina(guestmol)
            except:
                complexmols, scores = dock.score_map_comb(guestmol)

        # Complexes are returned as conformers
        for i in range(complexmols.GetNumConformers()):
            exo = is_exo(complexmols, self.hostmol, self.conf, confId=i)
            if not exo:
                complexmol = Chem.Mol(complexmols)
                complexmol.RemoveAllConformers()
                complexmol.AddConformer(complexmols.GetConformer(i), assignId=True)
                break
            # If all conformers are exo, end
            elif i == complexmols.GetNumConformers()-1:
                print("END - Exo complex")
                guestmol.SetDoubleProp("en", 20.3)
                return get_property_mol(guestmol), get_property_mol(guestmol)
        
        # 3. Calculate binding energy
        # Definitely set this up to reuse this object otherwise memory usage will get out of hand
        calc = xtbEnergy(self.conf)
        try:
            print("STATUS - optimising guest")
            optguestmol, guest_en = calc.get_opt(guestmol)
            print("STATUS - optimising complex")
            optcomplexmol, complex_en = calc.get_opt(complexmol)
        except:
            print("END - couldn't optimise")
            nullmol.SetDoubleProp("en", 20.4)
            return get_property_mol(nullmol), get_property_mol(nullmol)
        
        # If the result of xTB optimisation is exo, set bad score
        exo = is_exo(optcomplexmol, self.hostmol, self.conf)
        if exo:
            print("END - Exo complex")
            guestmol.SetDoubleProp("en", 20.5)
            return get_property_mol(guestmol), get_property_mol(guestmol)
        
        # Check if the smiles of the guest have changed
        changed = check_smiles_change(optguestmol, Chem.GetMolFrags(optcomplexmol, asMols=True)[1])
        if changed:
            print("END - Smiles changed")
            optcomplexmol.SetDoubleProp("en", 20.6)
            return get_property_mol(optcomplexmol), get_property_mol(optguestmol)
        
        bind_en = complex_en - guest_en - self.host_en
        optcomplexmol.SetDoubleProp("en", bind_en)
        optguestmol.SetDoubleProp("en", guest_en)

        # For docking testing, add scores from binding poses
        # try:
        #     for rank, i in enumerate(scores):
        #         optcomplexmol.SetDoubleProp("en_" + str(rank), i)
        # except:
        #     print("NOTE - couldn't add pose scores to mol object")

        return get_property_mol(optcomplexmol), get_property_mol(optguestmol)
    
if __name__ == "__main__":
    import os
    import argparse
    import yaml

    from smi2sdf import process_smi
    from docking import score_map_vina, score_map_comb, vina_dock
    from reward_utils import fix_charge, is_small_cylinder, is_exo, check_smiles_change, get_property_mol
    from xtb_opt import xtbEnergy
    
    os.environ["OPENBLAS_NUM_THREADS"] = "1"

    parser = argparse.ArgumentParser(
        description="",
        usage=f"python {os.path.basename(__file__)} -c <config_file> i <input_file>"
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
    
    # If using multiple properties, turn this into a dictionary
    host_en = conf["host_en"]
    print("globvars set")

    df = pd.read_csv(args.input, sep=" ", names=["smiles"], index_col=False)

    simulator = ChemSim(conf, hostmol, standalone=True)
    simulator.setup()
    for count, i in enumerate(df["smiles"],1):
        mol = Chem.MolFromSmiles(i)
        molout, guestmolout = simulator.run(mol)
        simulator.flush(molout)
        simulator.flush(guestmolout, guest=True)
        print("done ", count)