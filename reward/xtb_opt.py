""" XTB METHODS """
from rdkit import Chem
import os
import subprocess as sp
from glob import glob
from shutil import rmtree
import tempfile
import contextlib

""" HELPER FUNCTIONS FOR XTB OPTIMISATION METHODS """

@contextlib.contextmanager
def my_temp_dir(parentdir, delete=False):
    # Later python versions have a delete option on tempdir for debugging, creating one myself
    temp_dir = tempfile.mkdtemp(dir=parentdir)
    try:
        yield temp_dir
    finally:
        if delete:
            shutil.rmtree(temp_dir)


class xtbMol():
    """ Stores the rdkit mol
        In addition to objects required through xTB calculations
    """

    def __init__(self, mol, name, **kwargs):
        self.mol = mol
        self.name = name

        self.files = {}

        # Leaving kwargs in for redundancy
        self.__dict__.update(kwargs)

class xtbCalculator():
    """ Class for handling xTB energy calculation and output
    """
    def __init__(self, conf):
        self.conf = conf

        self.outdir = os.path.abspath(conf["output_dir"])
        self.hostdir = os.path.abspath(conf["host_dir"])
        self.host_sdf = os.path.abspath(conf["host_sdf"])
        self.xtb_inp = os.path.abspath(conf["xtb_inp"])

        # Get host energy
        hostmol = Chem.MolFromMolFile(conf["host_sdf"], removeHs=False)
        self.host = xtbMol(hostmol, "cb7")

        print("optimising host")
        orgdir = os.getcwd()
        os.chdir(self.hostdir)
        hostmol, host_en = self.get_opt(self.host)
        os.chdir(orgdir)
        
        self.hostmol = hostmol
        self.host_en = host_en
        self.host_en_dict = {propname : hostmol.GetProp(propname) for propname in hostmol.GetPropNames()}
        
        print("host optimised")

    def get_guest_complex_opt(self, comp, guest):
        """ Optimises the guest and complex
        """
        with my_temp_dir(parentdir=f"{self.outdir}", delete=False) as d:
            try:
                orgdir = os.getcwd()
                os.chdir(d)
                guestxtbmol = xtbMol(guest, "guest")
                complexxtbmol = xtbMol(comp, "complex")

                finalguestmol, guest_en = self.get_opt(guestxtbmol)
                finalcomplex, complex_en = self.get_opt(complexambermol)
                os.chdir(orgdir)

            except:
                os.chdir(orgdir)
                raise ValueError("Error in optimisation")

        return finalcomplexmol, finalguestmol

    def get_opt(self, xtbmol):
        """ Calls methods to optimise a mol and retrieve energy
        """
        with my_temp_dir(parentdir="./", delete=False) as d:
            orgdir = os.getcwd()
            os.chdir(d)
            Chem.MolToMolFile(xtbmol.mol, "mol.sdf", kekulize=False)

            # Check whether to add .CHRG file
            chrg = Chem.GetFormalCharge(xtbmol.mol)
            if chrg != 0:
                with open(".CHRG", "w") as fw:
                    fw.write(str(chrg))
                      
            self.xtb_opt("mol.sdf", self.xtb_inp, self.conf["optlevel"], self.conf["thermo"], self.conf["solvent"], self.conf["additional_ens"])
            finalmol = Chem.MolFromMolFile("xtbopt.sdf",removeHs=False,sanitize=False)
            
            en = self.get_en("opt.out", self.conf["thermo"])

            # Additional outputs are assigned as class attributes
            if self.conf["additional_ens"]:
                self.get_additional_ens(xtbmol)
                for key, val in xtbmol.en_dict.items():
                    finalmol.SetProp(f"xtb_{key}", str(val))
                
            os.chdir(orgdir)

        return finalmol, en

    @staticmethod
    def xtb_opt(filename, xtb_inp, optlevel, is_thermo, is_solvent, additional_ens):
        """ Optimises from an sdf file """

        sp.run(["obabel","-isdf","-osdf",f"{filename}","-O",f"{filename}"],stderr=sp.DEVNULL)
        if additional_ens:
            cmd = ["xtb", "--input", xtb_inp, f"{filename}"]
        else:
            cmd = ["xtb", f"{filename}"]

        if is_thermo:
            cmd.append("--ohess")
        else:
            cmd.append("--opt")

        cmd.append(optlevel)

        if is_solvent:
            cmd.append("--alpb")
            cmd.append("water")

        sp.run(cmd, stdout=open("opt.out","w"), stderr=sp.DEVNULL)
        
        return

    @staticmethod
    def get_en(outfile, is_thermo):
        """ Gets total energy from an xtb output file """
        gen = (i.split() for i in reversed(open(outfile,"r").readlines()))

        binding = 0
        for i in gen:
            if i:
                if is_thermo:
                    if " ".join(i[:4]) == "| TOTAL FREE ENERGY":
                        binding = float(i[4])*627.5095
                        break
                else:
                    if " ".join(i[:3]) == "| TOTAL ENERGY":
                        binding = float(i[3])*627.5095
                        break
        return binding

    @staticmethod
    def get_additional_ens(xtbmol):
        """ Gets additional energies from an xtb output file """
        gen = (i.split() for i in reversed(open("opt.out","r").readlines()))
        prop_gen = (i.split() for i in open("properties.out","r").readlines())

        # Get molecule name
        en_dict = {}

        for i in gen:
            if " ".join(i[:3]) == "| TOTAL ENERGY":
                en_dict["scf_en"] = float(i[3])

            if " ".join(i[:3]) == "| TOTAL ENTHALPY":
                en_dict["enthal"] = float(i[3])

            if " ".join(i[:3]) == ":: SCC energy":
                en_dict["scc_en"] = float(i[3])
                
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

                if " ".join(i[:3]) == ":: G(RRHO) contrib.":
                    en_dict["thermo"] = float(i[3])

                if " ".join(i[:3]) == ":: zero point":
                    en_dict["zero_point"] = float(i[4])

            except:
                continue

        setattr(xtbmol, "en_dict", en_dict)

        return

    # def get_partial_charges        
