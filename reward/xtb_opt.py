""" XTB METHODS """
from rdkit import Chem
import os
import subprocess as sp
from glob import glob
from shutil import rmtree
import tempfile

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
        with tempfile.TemporaryDirectory(dir=os.abspath(f"{self.outdir}")) as d:
            orgdir = os.getcwd()
            os.chdir(d)
            Chem.MolToMolFile(mol, "mol.sdf", kekulize=False)

            # Check whether to add .CHRG file
            chrg = Chem.GetFormalCharge(mol)
            if chrg != 0:
                with open(".CHRG", "w") as fw:
                    fw.write(str(chrg))
                      
            self.xtb_opt("mol.sdf")

            try:
                finalmol = Chem.MolFromMolFile("xtbopt.sdf",removeHs=False,sanitize=False)
            except:
                os.chdir(orgdir)
                raise ValueError

            en = self.get_en()

            # Additional outputs are assigned as class attributes
            if self.conf["additional_ens"]:
                self.get_additional_ens(finalmol)
            os.chdir(orgdir)

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
            os.chdir(orgdir)
            rmtree(f"{self.outdir}/xtbtmp_{dirId}")
            raise ValueError

        en = self.get_en()

        # Additional outputs are assigned as class attributes
        if self.conf["additional_ens"]:
            self.get_additional_ens(finalmol)

        # if self.conf["partial_charges"]:

        os.chdir(orgdir)
        rmtree(f"{self.outdir}/xtbtmp_{dirId}")
        """
        return finalmol, en

    def xtb_opt(self, filename):
        """ Optimises from an sdf file """

        sp.run(["obabel","-isdf","-osdf",f"{filename}","-O",f"{filename}"],stderr=sp.DEVNULL)
        if self.conf["additional_ens"]:
            cmd = ["xtb", "--input", "./../../../reward/xtb_additional.inp", f"{filename}"]
        else:
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
                    if " ".join(i[:4]) == "| TOTAL FREE ENERGY":
                        binding = float(i[4])*627.5095
                        break
                else:
                    if " ".join(i[:3]) == "| TOTAL ENERGY":
                        binding = float(i[3])*627.5095
                        break
        return binding
        
    def get_additional_ens(self, mol):
        """ Gets additional energies from an xtb output file """
        gen = (i.split() for i in reversed(open("opt.out","r").readlines()))
        prop_gen = (i.split() for i in open("properties.out","r").readlines())

        # Get molecule name
        en_dict = {}

        for i in gen:
            if " ".join(i[:3]) == ":: G(RRHO) contrib.":
                en_dict["thermo"] = float(i[3])

            if " ".join(i[:3]) == ":: zero point":
                en_dict["zero_point"] = float(i[4])

            if " ".join(i[:3]) == ":: repulsion energy":
                en_dict["repulsion"]  = float(i[3])

            if " ".join(i[:3]) == ":: -> Gshift":
                en_dict["Gshift"] = float(i[3])

            if " ".join(i[:3]) == ":: -> Ghb":
                en_dict["Ghb"] = float(i[3])

            if " ".join(i[:3]) == ":: -> Gsasa":
                en_dict["Gsasa"] = float(i[3])

            if " ".join(i[:3]) == ":: -> Gelec":
                en_dict["Gelec"] = float(i[3])

            if " ".join(i[:3]) == ":: -> Gsolv":
                en_dict["Gsolv"] = float(i[3])

            if " ".join(i[:3]) == ":: -> dispersion":
                en_dict["disp"] = float(i[3])

            if " ".join(i[:4]) == ":: -> anisotropic XC":
                en_dict["aniso_XC"] = float(i[4])

            if " ".join(i[:4]) == ":: -> anisotropic ES":
                en_dict["aniso_ES"] = float(i[4])

            if " ".join(i[:4]) == ":: -> isotropic ES":
                en_dict["iso_ES"] = float(i[4])
        
        for i in prop_gen:
            count = 0
            try:
                if i[0] == "full:":
                    if count == 0:
                        en_dict["dipole_x"] = float(i[1])
                        en_dict["dipole_y"] = float(i[2])
                        en_dict["dipole_z"] = float(i[3])
                        count += 1
                    if count == 1:
                        en_dict["quadrupole_xx"] = float(i[1])
                        en_dict["quadrupole_xy"] = float(i[2])
                        en_dict["quadrupole_yy"] = float(i[3])
                        en_dict["quadrupole_xz"] = float(i[4])
                        en_dict["quadrupole_yz"] = float(i[5])
                        en_dict["quadrupole_zz"] = float(i[6])
            except:
                continue

        for key, val in en_dict.items():
            mol.SetProp(f"xtb_{key}", str(val))    

        return

    # def get_partial_charges        
