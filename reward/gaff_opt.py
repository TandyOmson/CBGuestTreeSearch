""" Collection of python functions for processing files through AMBER with GAFF2
    - inputs are always files and a molecule name (converted to a 3 letter residue name in the program)
    - output are always files, which can be read by supplementary functions
"""
import subprocess as sp
import tempfile
import traceback
from collections import defaultdict
import os
import math

from rdkit import Chem

# TEMP FOR TESTING
from reward_utils import get_property_mol

def nonblank_lines(f, reverse=False):
    if not reverse:
        for l in f:
            line = l.rstrip()
            if line:
                yield line
    else:
        for l in reversed(f):
            line = l.rstrip()
            if line:
                yield line

def edit_resname(pdbfile, *constituents):
    with open(pdbfile, 'r') as f:
        lines = f.readlines()
    # Create a dictioanry with default value 0
    atom_counts = defaultdict(int)
    
    # Get the number of digits in the number of atoms to assign whitespace after atom name
    for line in lines:
        if line.startswith('HETATM'):
            atom_counts[line.split()[2]] += 1

    # Get max atom count    
    max_atom_count = max(atom_counts.values())
    whitespace = len(str(max_atom_count))

    # Get the column number of the residue number
    res_column = 24 + whitespace

    # Reset atom counts
    atom_counts = defaultdict(int)

    frag_num = 0
    j = 0
    for i, line in enumerate(lines):
        if line.startswith('HETATM'):
            # Remove species at the end
            line = line[:-4] + line[-3:]

            species = line.split()[2]

            # If the .pdb is a single molcule
            if len(constituents)==1:
                if atom_counts[species] == 0:
                    lines[i] = line[:13] + species + " "*(whitespace+1) + constituents[0].residue_name + line[27:]
                    
                else:
                    lines[i] =  line[:13] + species + str(atom_counts[species]) + " "*(whitespace+1-len(str(atom_counts[species]))) + constituents[0].residue_name + line[27:]
                atom_counts[species] += 1

                # Add residue numner at column 24
                lines[i] = lines[i][:res_column-1] + "1" + 6*" " + lines[i][res_column-1:]

            # If the .pdb is a complex
            elif len(constituents)>1:
                mol_frag = constituents[frag_num] 
                if j < mol_frag.mol.GetNumAtoms():
                    if atom_counts[species] == 0:
                        lines[i] = line[:13] + species + " "*(whitespace+1) + mol_frag.residue_name + line[27:]

                    else:
                        lines[i] = line[:13] + species + str(atom_counts[species]) + " "*(whitespace+1-len(str(atom_counts[species]))) + mol_frag.residue_name + line[27:]

                    atom_counts[species] += 1
                    # Add residue numner at relvant column
                    lines[i] = lines[i][:res_column-1] + str(frag_num+1) + 6*" " + lines[i][res_column-1:]
                    j += 1
                
                # Reset atom counts
                if i == mol_frag.mol.GetNumAtoms() - 1:
                    # Insert "TER" after the host molecule
                    lines.insert(i+1, "TER\n")
                    atom_counts = defaultdict(int)
                    frag_num += 1
                    j = 0

    #lines.insert(len(lines)-1, "TER\n")
                    
    with open(pdbfile, 'w') as f:
        f.writelines(lines)
    return

class AmberMol():
    """ Stores the rdkit mol (which stores properties)
        In addition to objects required through GAFF calculations
    """

    def __init__(self, mol, name, chgfile=None, **kwargs):
        # Currently fileinput must be an sdf
        self.mol = mol
        self.name = name
        self.residue_name = name.upper()[:3]
        
        self.files = {}
        # Add a partial charge file if available
        if chgfile:
            self.files["chgfile"] = chgfile
        
        # Leaving kwargs in for redundancy
        self.__dict__.update(kwargs)
        
class AmberCalculator():
    """ Class for running structure preparation and calculation in Amber24
        Inidividual components are static methods
        Combined calculations use the config in the object state
    """

    def __init__(self, conf):
        # Set options from configfile
        self.conf = conf

        # Non optional variables assigned as class attributes
        #self.is_solvent = conf["solvent"]
        #self.is_thermo = conf["thermo"]
        
        self.outdir = os.path.abspath(conf["output_dir"])
        self.hostdir = os.path.abspath(conf["host_dir"])
        self.host_sdf = os.path.abspath(conf["host_sdf"])
        self.min_file = os.path.abspath(conf["min_file"])

        # Get host energy
        hostmol = Chem.MolFromMolFile(conf["host_sdf"], removeHs=False)
        self.host = AmberMol(hostmol, "cb7")

        print("optimising host")
        orgdir = os.getcwd()
        os.chdir(self.hostdir)
        hostmol, host_en = self.get_opt(self.host)
        os.chdir(orgdir)
        
        self.hostmol = hostmol
        self.host_en = host_en
        self.host_en_dict = {propname : hostmol.GetProp(propname) for propname in hostmol.GetPropNames()}
        
        self.host.files["prmtop"] = self.hostdir + "/" + self.host.files["prmtop"]
        self.host.files["lib"] = self.hostdir + "/" + self.host.files["lib"]
        self.host.files["frcmod"] = self.hostdir + "/" + self.host.files["frcmod"]
        self.host.files["mol2"] = self.hostdir + "/" + self.host.files["mol2"]
        
        print("host optimised")
        
    # For gaff, the routine needs to keep the tempfiles of the guest to have the prmtop for the complex
    def get_guest_complex_opt(self, comp, guest):
        """ optimises the guest, keeping the prmtop for the complex, then optimises the complex
        """
        with tempfile.TemporaryDirectory(dir=f"{self.outdir}") as d:
            try:
                orgdir = os.getcwd()
                os.chdir(d)
                guestambermol = AmberMol(guest, "guest")
                complexambermol = AmberMol(comp, "complex")
                setattr(complexambermol, "constituents", [self.host, guestambermol])

                finalguestmol, guest_en = self.get_opt(guestambermol)
                finalcomplexmol, complex_en = self.get_opt(complexambermol)
                os.chdir(orgdir)
                
            except:
                os.chdir(orgdir)
                raise ValueError
            
        return finalcomplexmol, finalguestmol
            
    def get_opt(self, ambermol):
        """ Calls methods to optimise mol and retrieve energy
            This is all done in the same directory for a specific ambermol
        """
        Chem.MolToMolFile(ambermol.mol, f"{ambermol.name}.sdf", kekulize=False)

        if Chem.GetFormalCharge(ambermol.mol) != 0:
            self.conf["chg_method"] = "abcg2"

            self.ambermol_preprocess(ambermol)
            self.minimize_sander(ambermol, sander_file=self.min_file)
            self.get_opt_pdb(ambermol)
            self.get_en_sander(ambermol, outfile=ambermol.files["min_out_sander"])
            self.nmode_nab(ambermol)
            #breakpoint()
            self.get_nmode(ambermol, outfile=ambermol.files['nmode_out'])

            self.conf["chg_method"] = "gas"

        else:
            self.ambermol_preprocess(ambermol)
            self.minimize_sander(ambermol, sander_file=self.min_file)
            self.get_opt_pdb(ambermol)
            self.get_en_sander(ambermol, outfile=ambermol.files["min_out_sander"])
            self.nmode_nab(ambermol)
            #breakpoint()
            self.get_nmode(ambermol, outfile=ambermol.files['nmode_out'])

        en = ambermol.en_dict["en"]

        # Get final mol to return to program - ambermol is discarded
        sp.run(["obabel", "-ipdb", ambermol.files["min_pdb"], "-osdf", "-O", f"{ambermol.name}_min.sdf"], stdout=sp.DEVNULL, stderr=sp.DEVNULL)
        finalmol = Chem.MolFromMolFile(f"{ambermol.name}_min.sdf", removeHs=False, sanitize=False)
        
        for key, val in ambermol.en_dict.items():
            finalmol.SetProp(f"gaff_{key}", str(val))

        return finalmol, en

    @staticmethod
    def get_pdb(ambermol):
        """ in: rdkit mol
            out: amber formatted pdb file
        """
        ambermol.files["sdf"] = f"{ambermol.name}.sdf"
        ambermol.files["pdb"] = f"{ambermol.name}.pdb"
        ambermol.files["pdb4amber"] = f"{ambermol.name}_4amber.pdb"
        
        Chem.MolToMolFile(ambermol.mol, ambermol.files["sdf"])
        sp.run(["obabel", "-isdf", ambermol.files["sdf"], "-opdb", "-O", ambermol.files["pdb"]], stdout=sp.DEVNULL, stderr=sp.DEVNULL)
        sp.run(["pdb4amber", "-i", ambermol.files["pdb"], "-o", ambermol.files["pdb4amber"]], stdout=sp.DEVNULL, stderr=sp.DEVNULL)

        # Need to edit residue names if the molecule is a complex
        if hasattr(ambermol, "constituents"):
            edit_resname(ambermol.files["pdb4amber"], *ambermol.constituents)
        else:
            edit_resname(ambermol.files["pdb4amber"], ambermol)

    @staticmethod
    def sovlate_packmol_memgen(ambermol):
        pass

    @staticmethod
    def antechamber(ambermol, conf):
        """ in: amber formated pdb
            out: antechamber made mol2

        """
        ambermol.files["mol2"] = f"{ambermol.name}.mol2"
               
        ante_cmd = ["antechamber", "-i", ambermol.files["pdb4amber"], "-fi", "pdb", "-fo", "mol2", "-o", ambermol.files["mol2"]]
            
        if conf["chg_method"] == "read":
            ante_cmd.extend(["-c", "rc", "-cf", ambermol.files["chgfile"]]) 
            
        else:
            chg_method = conf["chg_method"]
            ante_cmd.extend(["-c", chg_method])

        ante_cmd.extend(["-nc", str(Chem.GetFormalCharge(ambermol.mol)), "-s", "2", "-pf", "yes", "-rn", ambermol.residue_name, "-at", "gaff2"])

        sp.run(ante_cmd, stdout=sp.DEVNULL)

    @staticmethod
    def parmchk(ambermol):
        """ in: antechamber formatted mol2 file
            out: frcmod (modifications to force field parameters)
        """
        ambermol.files["frcmod"] = f"{ambermol.name}.frcmod"
    
        sp.run(["parmchk2", "-i", ambermol.files["mol2"], "-f", "mol2", "-o", ambermol.files["frcmod"]], stdout=sp.DEVNULL)

    @staticmethod
    def tleap_single(ambermol):
        """ in: frcmod, mol2
            out: prmtop and rst7 files for input into sander
        """
        ambermol.files["lib"] = f"{ambermol.name}.lib"
        ambermol.files["prmtop"] = f"{ambermol.name}.prmtop"
        ambermol.files["rst7"] = f"{ambermol.name}.rst7"
               
        tleap_script = ["source leaprc.gaff2",
                        f"loadamberparams {ambermol.files['frcmod']}",
                        f"{ambermol.residue_name} = loadmol2 {ambermol.files['mol2']}",
                        f"check {ambermol.residue_name}",
                        f"saveoff {ambermol.residue_name} {ambermol.files['lib']}",
                        f"saveamberparm {ambermol.residue_name} {ambermol.files['prmtop']} {ambermol.files['rst7']}",
                        "quit"
                        ]

        with open(f"leap.in", "w") as fw:
            for i in tleap_script:
               fw.write(f"{i}\n")

        sp.run(["tleap", "-f", f"leap.in"], stdout=sp.DEVNULL)

    @staticmethod
    def tleap_complex(complexambermol, *constituents):
        # Instead of creating a lib file, load in the lib files from the individual constituents
        # These should be in the ambermol.files dictionary of the complex ambermol object
        complexambermol.files["mol2"] = f"{complexambermol.name}.mol2"
        complexambermol.files["prmtop"] = f"{complexambermol.name}.prmtop"
        complexambermol.files["rst7"] = f"{complexambermol.name}.rst7"
        complexambermol.residue_name = "COM"

        tleap_script = ["source leaprc.gaff2"]

        tleap_script.extend([f"loadoff {ambermol.files['lib']}" for ambermol in constituents])
        tleap_script.extend([f"loadamberparams {ambermol.files['frcmod']}" for ambermol in constituents])

        tleap_script.extend([f"{complexambermol.residue_name} = loadPDB {complexambermol.files['pdb4amber']}",
                             f"savemol2 {complexambermol.residue_name} {complexambermol.files['mol2']} 1",
                             f"saveamberparm {complexambermol.residue_name} {complexambermol.files['prmtop']} {complexambermol.files['rst7']}",
                             "quit"
        ])

        with open(f"leap.in", "w") as fw:
            for i in tleap_script:
               fw.write(f"{i}\n")

        sp.run(["tleap", "-f", f"leap.in"], stdout=sp.DEVNULL)

    @staticmethod
    def minimize_sander(ambermol, sander_file):
        ambermol.files["ncrst"] = f"{ambermol.name}.ncrst"
        ambermol.files["min_traj"] = f"{ambermol.name}.traj"
        ambermol.files["min_out_sander"] = f"{ambermol.name}_sander_min.log"

        sp.run(["sander", "-O",
                "-i", sander_file,
                "-o", ambermol.files['min_traj'],
                "-p", ambermol.files['prmtop'],
                "-c", ambermol.files['rst7'],
                "-r", ambermol.files['ncrst']],
               stdout=open(f"min.err", "w"),
               stderr=open(f"min.err", "w"),
        )
        sp.run(["mv", "mdinfo", ambermol.files["min_out_sander"]])

    @staticmethod
    def minimize_nab(ambermol, script):
        ambermol.files["ncrst"] = f"{ambermol.name}.ncrst"
        ambermol.files["min_out"] = f"{ambermol.name}_min.log"

        # NOTE the outfiles with have nothing in them if there is a segmentation fault during nab script execution
        sp.run([os.path.abspath(script),
                ambermol.files["prmtop"],
                ambermol.files["rst7"],
                ambermol.files["ncrst"]],
               stdout=open(ambermol.files["min_out"], "w"),
               stderr=open(ambermol.files["min_out"], "w")
        )

    @staticmethod
    def nmode_nab_legacy(ambermol, script):
        # Calculates normal modes from the nmode nab script
        ambermol.files["nmode_out"] = f"{ambermol.name}_nmode.log"

        sp.run(["./nab_scripts/nmode",
                ambermol.files["prmtop"],
                ambermol.files["ncrst"]],
               stdout=open(ambermol.files["nmode_out"], "w"),
               stderr=open(ambermol.files["nmode_out"], "w")
        )

        # Cleanup of undeeded file
        sp.run(["rm", "vecs"])

    @staticmethod
    def nmode_nab(ambermol):
        # Runs conjugate gradient descent minimisation, followed by newton raphson minimisation, followed by nmode
        # Uses a nab script, requires original nab compiler. Am planning to get a nabc script working if possible
        ambermol.files["nmode_out"] = f"{ambermol.name}_nmode.log"
        ambermol.files["min_pdb"] = f"{ambermol.name}_min.pdb"
        
        nab_script = [
            'molecule m;',
            'float x[dynamic], f[dynamic], fret;',
        
            'm = getpdb( "{0}" );',
            'readparm( m, "{1}");',
            'allocate x[ 3*m.natoms ];  allocate f[ 3*m.natoms ];'
            'setxyz_from_mol( m, NULL, x );',

            'mm_options( "ntpr=5000, nsnb=99999, gb=1, gbsa=1, cut=99.0, rgbmax=99.0, diel=C" );',
            'mme_init( m, NULL, ":::", x, NULL);',
        
            '// conjugate gradient minimization',
            # 'dgrad = 3*natoms*0.001;',
            '//conjgrad(x, 3*m.natoms, fret, mme, 0.00001, 0.00001, 100000 );',
        
            '// Newton-Raphson minimization\fP',
            'mm_options( "ntpr=100" );',
            '//newton( x, 3*m.natoms, fret, mme, mme2, 0.0000001, 0.0, 2000 );',

            '// Write out the minimized pdb',
            '//putpdb( "{2}", m );',
        
            '// get the normal modes:',
            'nmode( x, 3*m.natoms, mme2, 0, 0, 0.0, 0.0, 0);'
        ]
        
        nab_script = [x.format(ambermol.files["min_pdb"], ambermol.files["prmtop"], ambermol.files["min_pdb"])+'\n' for x in nab_script]
        
        # Generate nab script (for some reason have to hard code argv, another thing to try and fix
        with open(f"min_prog.nab", 'w') as w:
            w.writelines(nab_script)
        
        # Compile nab script (ensure the original nab compiler is in path)
        sp.run(["nab", f"min_prog.nab", "-o", f"min_prog"])

        # Run nab script
        sp.run([f"./min_prog"],
               stdout=open(ambermol.files["nmode_out"], "w"),
               stderr=open(ambermol.files["nmode_out"], "w")
        )
        
        return

    # CLASS METHODS FOR AMBER CALCULATIONS
    def ambermol_preprocess(self, ambermol):

        # Minimsise single ambermol
        if len(Chem.GetMolFrags(ambermol.mol, asMols=True)) == 1:
            AmberCalculator.get_pdb(ambermol)
            AmberCalculator.antechamber(ambermol, self.conf)
            AmberCalculator.parmchk(ambermol)
            AmberCalculator.tleap_single(ambermol)

        # Minmise a complex of moecules that are not bound covalently
        if len(Chem.GetMolFrags(ambermol.mol, asMols=True)) > 1:
            AmberCalculator.get_pdb(ambermol)
            AmberCalculator.tleap_complex(ambermol, *ambermol.constituents)

    # CLASS METHODS FOR RETRIEVAL OF OUTPUTS
    def get_opt_pdb(self, ambermol):
        ambermol.files["min_pdb"] = f"{ambermol.name}_min.pdb"
        
        sp.run(["ambpdb", "-p", ambermol.files['prmtop'], "-c", ambermol.files["ncrst"]], stdout=open(ambermol.files["min_pdb"], "w"))
        
    def get_en_sander(self, ambermol, outfile):
        """ Retrieves energy from the "mdinfo" file output by sander (this should be renamed by now)
        """

        min_out = [i.split() for i in open(outfile, "r")]
        en_dict = {}

        en_dict["en"] = float(min_out[3][1])
        en_dict["Ebond"] = float(min_out[5][2]) # bond stretching
        en_dict["Eangle"] = float(min_out[5][5]) # angle bending
        en_dict["Edihed"] = float(min_out[5][8]) # intrinsic dihedral
        en_dict["Evdw"] = float(min_out[6][2]) # van der waals
        en_dict["Ecoul"] = float(min_out[6][5]) # coloumb energy
        en_dict["Egb"] = float(min_out[6][8]) # Generalised born energy
        en_dict["Evdw_14"] = float(min_out[7][3]) # topological 1-4 van der waals (atoms separated by 3 bonds)
        en_dict["Ecoul_14"] = float(min_out[7][7]) # topological 1-4 coulomb (atoms separated by 3 bonds)
        en_dict["Enp"] = float(min_out[8][2]) # Non polar solvent energy, associate with surface area in contact with solvent

        setattr(ambermol, "en_dict", en_dict)

    def get_en_nab(self, ambermol, outfile):
        """ Energy from nab xmin script
        """

        min_out = [i.split() for i in open(outfile, "r")]
        en_dict = {}
        
        for i in min_out:
            if i[0] == "bond:":
                en_dict["Ebond"] = float(i[1]) # bond stretching
            if i[0] == "angle:":
                en_dict["Eangle"] = float(i[1]) # angle bending
            if i[0] == "dihedral:":
                en_dict["Edihed"] = float(i[1]) # intrinsic dihedral
            if i[0] == "enb_14:":
                en_dict["Evdw_14"] = float(i[1]) # topological 1-4 van der waals (atoms separated by 3 bonds)
            if i[0] == "eel14:":
                en_dict["Ecoul_14"] = float(i[1]) # topological 1-4 coulomb (atoms separated by 3 bonds)
            if i[0] == "enb:":
                en_dict["Evdw"] = float(i[1]) # Van der Waals (non-bonded)
            if i[0] == "eel:":
                en_dict["Ecoul"] = float(i[1]) # coloumb energy
            if i[0] == "egb:":
                en_dict["Egb"] = float(i[1]) # Generalised born energy
            if i[0] == "esurf:":
                en_dict["Enp"] = float(i[1]) # Non polar solvent energy, associate with surface area in contact with solvent
            if i[0] == "Total:":
                en_dict["en"] = float(i[1])
                
        setattr(ambermol, "en_dict", en_dict)

        return

    def get_en_conj_newt(self, ambermol, outfile):
        """ Gets energy AND NORMAL MODES from the conjugate gradient + newton raphson routine
        """
        min_out = reversed([i.split() for i in open(outfile, "r") if i])

        # Minimised energies
        en_dict = {}
        for i in min_out:
            if i == []:
                continue
            if i[0] == "ff:":
                en_dict["en"] = float(i[2])
                en_dict["Ebad"] = float(i[3]) # bond angle dihedral
                en_dict["Evdw"] = float(i[4]) # van der waals
                en_dict["Ecoul"] = float(i[5])
                en_dict["Enp"] = float(i[6])
                en_dict["Egb"] = float(i[7])
                break

        # Get vibrational mode frequencies
        frequencies = []

        nmode_out = (i.split() for i in open(outfile, "r"))
        i = next(nmode_out)

        while i[0] != "ff":
            i = next(nmode_out)
            if not bool(i):
                i = next(nmode_out)
        done = False
        while done == False:
            try:
                i = next(nmode_out)
                frequencies.append(i[1])
            except StopIteration:
                done = True
            
        # Valid frequencies will depend on number of degrees of freedom lost whether it is a complex, 6 for guest
        valid_freqs = 6
        try:
            extra_frags = GetMolFrags(ambermol.mol, asMols=True) - 1
        except:
            extra_frags = 0
        
        valid_freqs += extra_frags*6

        # Convert frequencies to appropraite units, account for loss of degrees of freedom
        frequencies = [float(x)*1e2*2.99e8 for x in frequencies[valid_freqs:]]
        
        k = 1.38E-23
        R = 8.3145
        T = 298.15
        beta = 1/(k*T)
        h = 6.626E-34
        Srrho = 0
        for f in frequencies:
            inc = R*(beta*h*f*(math.exp(beta*h*f)-1)**-1 - math.log(1-math.exp(-beta*h*f)))
            Srrho = inc+Srrho
        
        #  kcal/mol
        Srrho = Srrho/4.184*T/1000
        en_dict["Esrrho"] = Srrho
        
        en_dict["Ecorr"]  = en_dict["en"] - en_dict["Esrrho"]

        setattr(ambermol, "en_dict", en_dict)
            
        return

    def get_nmode(self, ambermol, outfile):
        """ Retrieves only normal modes from output of a nab nmode calculation
        """
        min_out = [i.split() for i in open(outfile, "r")]
        
        # Get vibrational mode frequencies
        frequencies = []

        nmode_out = (i.split() for i in open(outfile, "r"))
        i = next(nmode_out)

        while i[0] != "ff":
            i = next(nmode_out)
            if not bool(i):
                i = next(nmode_out)
        done = False
        while done == False:
            try:
                i = next(nmode_out)
                frequencies.append(i[1])
            except StopIteration:
                done = True
            
        # Valid frequencies will depend on number of degrees of freedom lost whether it is a complex, 6 for guest
        valid_freqs = 6
        try:
            extra_frags = GetMolFrags(ambermol.mol, asMols=True) - 1
        except:
            extra_frags = 0
        
        valid_freqs += extra_frags*6

        # Convert frequencies to appropraite units, account for loss of degrees of freedom
        frequencies = [float(x)*1e2*2.99e8 for x in frequencies[valid_freqs:]]
        
        k = 1.38E-23
        R = 8.3145
        T = 298.15
        beta = 1/(k*T)
        h = 6.626E-34
        Srrho = 0
        for f in frequencies:
            inc = R*(beta*h*f*(math.exp(beta*h*f)-1)**-1 - math.log(1-math.exp(-beta*h*f)))
            Srrho = inc+Srrho
        
        #  kcal/mol
        Srrho = Srrho/4.184*T/1000
        ambermol.en_dict["Esrrho"] = Srrho

        ambermol.en_dict["Ecorr"]  = ambermol.en_dict["en"] + ambermol.en_dict["Esrrho"]

        return

if __name__ == "__main__":

    # May wish to change antechamber conf to file input as with sander
    gaff_conf = {"chg_method" : "gas",
                 "min_file" : os.path.abspath("min.in"),
                 "output_dir" : "temp",
                 "host_dir" : os.path.abspath("data/gaff_host"),
                 "host_sdf" : os.path.abspath("data/host_aligned.sdf"),
    }

    rundir = os.getcwd()
    # Initialising the caclulator optimises the host
    gaff_calc = AmberCalculator(gaff_conf)
    propertymol = get_property_mol(gaff_calc.hostmol)
    hostmoldict = {}
    for i in propertymol.GetPropNames():
        try:
            hostmoldict[i] = float(propertymol.GetProp(i))
        except:
            hostmoldict[i] = propertymol.GetProp(i)

    # Test some guests
    import pandas as pd

    example_complexmols = Chem.SDMolSupplier("./temp/10_complexes.sdf", removeHs=False, sanitize=False)
    example_guestmols = Chem.SDMolSupplier("./temp/10_guests.sdf", removeHs=False, sanitize=False)

    bindings = {}

    for count, (complexmol, guestmol) in enumerate(zip(example_complexmols, example_guestmols)):
        try:
            # Gasteiger charges assume total charge is 0, need to switch to a slower method
            if Chem.GetFormalCharge(complexmol) != 0:
                gaff_conf["chg_method"] = "abcg2"
                
            finalcomplex, finalguests = gaff_calc.get_guest_complex_opt(complexmol, guestmol)

            propertymol = get_property_mol(finalcomplex)
            complexmoldict = {}
            for i in propertymol.GetPropNames():
                try:
                    complexmoldict[i] = float(propertymol.GetProp(i))
                except:
                    complexmoldict[i] = propertymol.GetProp(i)

            print("host:", hostmoldict)
            print("complex:", complexmoldict)

            guestmoldicts = []
            propertymols = [get_property_mol(guest) for guest in finalguests]
            for propmol in propertymols:
                guestmoldict = {}
                for i in propmol.GetPropNames():
                    try:
                        guestmoldict[i] = float(propmol.GetProp(i))
                    except:
                        guestmoldict[i] = propmol.GetProp(i)

                print("guest:", guestmoldict)
                guestmoldicts.append(guestmoldict)

            binding_dict = {key_c: c - h - g for ((key_c, c), (key_h, h), (key_g, g)) in zip(complexmoldict.items(), hostmoldict.items(), guestmoldicts[0].items())}
            print(binding_dict)
            bindings[count] = binding_dict["gaff_en"]

            with open("./temp/gaff_ens.csv", "w") as fw:
                for key, val in bindings.items():
                    fw.write(f"{key}\t{val}\n")
                    
        except Exception as e:
            print(e)
            print(traceback.format_exc())
            os.chdir(rundir)
            continue
