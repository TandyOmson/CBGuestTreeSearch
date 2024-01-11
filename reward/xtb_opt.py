""" XTB METHODS """
from rdkit import Chem
import os
import subprocess as sp
from glob import glob
from shutil import rmtree

class xtbEnergy():
    """ Class for handling xTB energy calculation and output
    """
    def __init__(self, conf):
        self.conf = conf
        # Non optional variables assigned as class attributes
        self.is_solvent = conf["solvent"]
        self.optlevel = conf["optlevel"]
        self.is_thermo = conf["thermo"]
        self.outdir = conf["output_dir"]

    def get_opt(self, mol):
        """ Calls methods to optimise a mol and retrieve energy
        """
        orgdir = os.getcwd()
        # Get current temp dirs
        curr = glob(f"{self.outdir}/xtbtmp_*")
        # Make new temp dir
        dirId = len(curr)+1
        os.mkdir(f"{self.outdir}/xtbtmp_{dirId}")
        os.chdir(f"{self.outdir}/xtbtmp_{dirId}")

        Chem.MolToMolFile(mol,"mol.sdf",kekulize=False)
        self.xtb_opt("mol.sdf")
        try:
            finalmol = Chem.MolFromMolFile("xtbopt.sdf",removeHs=False,sanitize=False)
        except:
            print("didnt optimise")

        en = self.get_en()

        # Additional outputs are assigned as class attributes
        # if self.conf["additional_ens"]:

        # if self.conf["partial_charges"]:

        os.chdir(orgdir)
        rmtree(f"{self.outdir}/xtbtmp_{dirId}")
        
        return finalmol, en

    def xtb_opt(self, filename):
        """ Optimises from an sdf file """

        sp.run(["obabel","-isdf","-osdf",f"{filename}","-O",f"{filename}"],stderr=sp.DEVNULL)
        cmd = ["xtb", f"{filename}"]

        if self.is_thermo:
            cmd.append("--ohess")
        else:
            cmd.append("--opt")

        cmd.append(self.optlevel)

        if self.is_solvent:
            cmd.append("--alpb")
            cmd.append("water")

        sp.run(cmd, stdout=open("opt.out","w"), stderr=sp.DEVNULL)
        
        return

    def get_en(self):
        """ Gets total energy from an xtb output file """
        gen = (i.split() for i in reversed(open("opt.out","r").readlines()))

        binding = 0
        for i in gen:
            if i:
                if self.is_thermo:
                    if " ".join(i[:4]) == ":: total free energy":
                        binding = float(i[4])*627.5095
                        break
                else:
                    if " ".join(i[:3]) == ":: total energy":
                        binding = float(i[3])*627.5095
                        break
        return binding
        
        # def get_additional_ens

        # def get_partial_charges        
