import os
import tempfile
import vina
import subprocess as sp
import traceback
import pickle

import vina
from meeko import MoleculePreparation, PDBQTWriterLegacy, PDBQTMolecule, RDKitMolCreate
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFingerprintGenerator

from chemtsv2.abc import Reward
from reward.vina_boltz_utils import *
from reward.descriptor_calcs import *

# TEMPORARY FIX FOR INCONSISTENT VERSION WARNING WHEN WORKING WITH PIPE AND KDE
# I need to just create pipe and kde in sklearn version 1.7.0
import warnings
from sklearn.exceptions import InconsistentVersionWarning

# Suppress only InconsistentVersionWarning
warnings.filterwarnings('ignore', category=InconsistentVersionWarning)

def init_pipe_kde(conf):
    with open(conf["mol_embed_pipeline"], "rb") as fr:
        pipe = pickle.load(fr)

    with open(conf["kde"], "rb") as fr:
        kde = pickle.load(fr)

    return pipe, kde
            
# Class for functions
class CrestVinaCalc():
    def __init__(self, rundir):
        # rundir is the directory to place all files and run subprocesses
        self.rundir = rundir

    def run_rdkit_confgen(self, guest, n_confs, nprocs):
        """ Embeds mutliple conformers and optimizes them with MMFF94
            Returns conformers and energies
        """
        AllChem.EmbedMultipleConfs(guest, n_confs)
        ens = [i[1] for i in AllChem.MMFFOptimizeMoleculeConfs(guest, numThreads=nprocs)]
        confs = [conf for conf in guest.GetConformers()]
        
        mols = []
        for i, conf in enumerate(confs):
            mol = Chem.Mol(guest)
            mol.RemoveAllConformers()
            mol.AddConformer(conf)
            mols.append(mol)

        return mols, ens

    def run_crest(self, guest, n_confs, nprocs):
        """ Runs CREST on a guest
            Returns conformers and energies
        """
        Chem.MolToMolFile(guest, f"{self.rundir}/guest_for_crest.sdf")
        sp.run(f"crest guest_for_crest.sdf --mdlen x0.5 gfn2 --gbsa h2o -T {nprocs}".split(" "), cwd=self.rundir, stdout=open(f"{self.rundir}/crest.out", "w"))
        try:
            mols = Chem.SDMolSupplier(f"{self.rundir}/crest_conformers.sdf", removeHs=False, sanitize=False)
            mols = [add_nitrogen_charges(mol) for mol in mols]
            ens = [float(i.split()[1]) for i in open(f"{self.rundir}/crest.energies", "r").readlines()]
        except:
            print("Error in CREST output")
            raise Exception
        
        sample_indices = stratified_weighted_sample(ens, total_samples=n_confs, bias_strength=0.05)
        mols = [mols[i] for i in sample_indices]
        ens = [ens[i] for i in  sample_indices]
                
        return mols, ens

    def run_vina(self, repsdffile, ligsdffile, n_comps, inp):
        """ Dock a molecule into a host
        hostfile is a .pdbqt file, use command "obabel host.pdb -O host.pdbqt -xr -xh" to convert
        Sets the vina interaction energy as a MolDoubleProp
        """        
        sp.run(["obabel", repsdffile, "-O", "host.pdbqt", "-xr", "-xh"], cwd=self.rundir, stderr=open(f"{self.rundir}/vina.log", "w"))
        sp.run(["obabel", ligsdffile, "-O", "guest.pdbqt", "-xh"], cwd=self.rundir, stderr=open(f"{self.rundir}/vina.log", "a"))
        
        # Run Vina
        sp.run(["vina", f"--config={inp}", "--num_modes=100", "--energy_range=25", "--receptor=host.pdbqt", f"--ligand=guest.pdbqt", f"--out=poses.pdbqt"], cwd=self.rundir, stdout=open(f"{self.rundir}/vina.log","a"))
        
        # Revert output .pdbqt to .sdf
        sp.run(["obabel", "poses.pdbqt", "-O", "poses.sdf"], cwd=self.rundir, stderr=open(f"{self.rundir}/vina.log","a"))
        
        # Read in all binding poses
        guestmols = Chem.SDMolSupplier(f"{self.rundir}/poses.sdf",removeHs=False,sanitize=False)
        hostmol = Chem.MolFromMolFile(f"{repsdffile}",removeHs=False,sanitize=False)
        
        complexmols = []
        exos = []
        for count, mol in enumerate(guestmols):
            comp = Chem.CombineMols(hostmol, mol)
            comp = add_nitrogen_charges(comp) 
            if is_exo(comp) is True:
                exos.append(count)
            else:
                complexmols.append(comp)
        
        if len(complexmols) != 0:
            vina_ens, rmsd_lbs, rmsd_ubs = parse_vina_output(f"{self.rundir}/vina.log", exo_ignore=exos)        
        
        # use aligment docking as as last resort
        else:
            print("docking by alignment as a last resort...")
            one_guestmol = Chem.MolFromMolFile(ligsdffile, removeHs=False)
            one_guestmol = add_nitrogen_charges(one_guestmol)
            align_pose = PCA_align_pose(hostmol, one_guestmol)
            complexmols, vina_ens = self.MMFF94_vina_opt(align_pose, f"{self.rundir}/host.pdbqt")
                
        return complexmols, vina_ens

    def MMFF94_vina_opt(self, mol, hostpdbqtfile):
        """ Runs MMFF94, then optimisation on a binding pose
            Modifies the pose in the conformers and adds a vina score
        """
        hostmol, guestmol = Chem.GetMolFrags(mol, asMols=True)
    
        mmff_res = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0, ignoreInterfragInteractions=False)
        
        preparator = MoleculePreparation()
        mol_setups = preparator.prepare(guestmol)
        
        for setup in mol_setups:
            pdbqt_string, is_ok, err_msg = PDBQTWriterLegacy.write_string(setup, bad_charge_ok=True)
    
        if not is_ok:
            print(err_msg)
    
        vinaobj = vina.Vina(verbosity=0)
        vinaobj.set_receptor(hostpdbqtfile)
        vinaobj.compute_vina_maps(center=[0.0,0.0,0.0], box_size=[20.0, 20.0, 20.0])
        
        vinaobj.set_ligand_from_string(pdbqt_string)
        vina_decomp = vinaobj.optimize()
        vina_en = -vina_decomp[0]
    
        vinaobj.write_pose(f"{self.rundir}/temp_align_pose.pdbqt", overwrite=True, remarks="")
        pdbqt_mol = PDBQTMolecule.from_file(f"{self.rundir}/temp_align_pose.pdbqt", skip_typing=True)
        complexmol = Chem.CombineMols(hostmol, RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)[0])
    
        # Returns an rdkit mol
        return [complexmol], [vina_en]

class Vina_reward(Reward):
    def get_objective_functions(conf):
        def VinaScore(mol):
            conf["hostfile"] = f"{conf['vina_boltz_utils']}/host.sdf"
            conf["vina_input"] = f"{conf['vina_boltz_utils']}/vinadock.inp"
            # Command for calculating symmetry number
            conf["sym_cmd"] = f"{conf['vina_boltz_utils']}/get_symmetry -maxaxisorder 10 -maxoptcycles 100 -same 0.001 -primary 0.5 -final 0.1 -maxoptstep 0.5 -minoptstep 1.0e-07 -gradstep 1.0e-7 -minchange 1.0e-8 -minchgcycles 5 --"
            
            mol = Chem.AddHs(mol)
            with tempfile.TemporaryDirectory() as tmpdir:
                if conf["debug"]:
                    print(f"[DEBUG] Using temporary directory: {tmpdir}")
                    
                calc = CrestVinaCalc(tmpdir)
                try:
                    res = AllChem.EmbedMolecule(mol)
                    if res == -1:
                        res = AllChem.EmbedMolecule(mol, useBasicKnowledge=False)
                        if res == -1:
                            raise Exception
                        else:
                            AllChem.MMFFOptimizeMolecule(mol)
                            
                    mol = add_nitrogen_charges(mol)
                
                    if conf["debug"]:
                        print("Running CREST...")

                    if conf["confgen"] == "crest":
                        crest_confs, crest_ens = calc.run_crest(mol, conf["num_crest_confs"], conf["num_crest_procs"])
                    elif conf["confgen"] == "rdkit":
                        crest_confs, crest_ens = calc.run_rdkit_confgen(mol, conf["num_crest_confs"], conf["num_crest_procs"])
                    else:
                        print("conformer generation (confgen) option not specified, closing...")
                        raise Exception

                    complexmols = []
                    complexens = []
                
                    if conf["debug"]:
                        print(f"{len(crest_confs)} crest_ens: {crest_ens}") 
                    
                    for m in crest_confs:
                        repsdffile, ligsdffile = conf["hostfile"], f"{tmpdir}/vina_temp.sdf"
                        Chem.MolToMolFile(m, ligsdffile)
                        vina_confs, vina_ens = calc.run_vina(repsdffile, ligsdffile, conf["vina_comps"], conf["vina_input"])
                        complexmols.extend(vina_confs)
                        complexens.extend(vina_ens)
                        
                    if conf["debug"]:
                        print(f"Vina Docking energies: {vina_ens}")
                
                    boltzmols, boltzens = filter_by_rmsd(complexmols, complexens, conf["boltz_comps"])
                    if conf["debug"]:
                        print(f"Sampled boltzmann energies: {boltzens}")
                
                    sym_num = get_symmetry_number(conf["sym_cmd"], crest_confs[0], conf["vina_boltz_utils"])
                    degen = sym_num
                    if conf["debug"]:
                        print(f"Symmetry number: {sym_num}")
                
                    boltzdiff = [i-min(boltzens) for i in boltzens]
                    Zr, dE, TS = get_energy_change(boltzdiff, [degen]*len(boltzdiff))
                    if conf["debug"]:
                        print("min", min(boltzens))
                        print("boltzmann", Zr, dE, TS)
                
                    final_en = min(boltzens) + dE - TS

                    # LFER calculations on the best docked vina pose
                    docked_guest = Chem.GetMolFrags(boltzmols[0], asMols=True, sanitizeFrags=False)[-1]
                    docked_guest = add_nitrogen_charges(docked_guest)

                    Mor29p, Mor24u, H046, Km = get_descs_hydrocarbons(docked_guest)

                    final_en = final_en - Mor29p*3.52995297 + Mor24u*1.3107841 - H046*0.27674704 + Km*1.65390875 - 1.72209064
                        
                    # save best pose
                    pose_dir = f"{conf['output_dir']}/3D_pose"
                    if not os.path.exists(pose_dir):
                        os.mkdir(pose_dir)
                    pose_file_name = f"{pose_dir}/mol_{conf['gid']}.sdf"
                    Chem.MolToMolFile(boltzmols[0], pose_file_name)
#                    writer = Chem.SDWriter(pose_file_name)
#                    for i in boltzmols:
#                        writer.write(i)
#                    writer.close()
                
                    if conf["debug"]:
                        print(f"final vina score: {final_en}")
                
                    return final_en
                
                except Exception as e:
                    print(f"Error SMILES: {Chem.MolToSmiles(mol)}")
                    pose_dir = f"{conf['output_dir']}/3D_pose"
                    pose_file_name = f"{pose_dir}/mol_{conf['gid']}.sdf"
                    nullmol = Chem.MolFromSmiles("C")
                    Chem.MolToMolFile(nullmol, pose_file_name)
                    traceback.print_exc()
                    return None

        def uniq_score(mol):
            """ An objective function calculating nearby density and decomposed space position
                based on training set kde fitting
            """
            
            pipe, kde = init_pipe_kde(conf)
            # molecule fingerprint
            fp = [rdFingerprintGenerator.GetAtomPairGenerator().GetFingerprint(mol)]
            # position in pca space
            X_pca = pipe.transform(fp)

            nearby_density = np.exp(kde.score_samples(X_pca)[0])
            
            return nearby_density, conf["max_density"], X_pca[0]
        
        return [VinaScore, uniq_score]

    def calc_reward_from_objective_values(values, conf):
        min_inter_score = values[0]

        # Options of how to use kde density data
        max_density = values[1][1]
        
        # maximise distance from a chosen point in decomposed space
        if conf["max_density_point"] != -1:
            X_pca = values[1][2]
            distance = 1 - float(np.linalg.norm(np.array(conf["max_density_point"]) - X_pca))
            u_base = 0.6
        # minimise distnace from a chosen point in decomposed space
        elif conf["min_density_point"] != -1:
            X_pca= values[1][2]
            distance = float(np.linalg.norm(np.array(conf["min_density_point"]) - X_pca))
            u_base = 0.3
        # minimise density
        elif conf["max_density"] != -1:
            u_score = values[1][0]
            u_base = 0.3
            
        if min_inter_score is None:
            return -1
        
        score_diff = min_inter_score - conf["vina_base_score"]
        # Determines how quickly reward changes around the base score in either direction (very sensitive)
        steepness = 0.3
        vina_reward =  -score_diff * steepness / (1 + abs(score_diff) * steepness)

        u_reward = 1 - ((u_score*(1+((1/u_base)*max_density)))/(u_score + max_density))

        return vina_reward + u_reward
