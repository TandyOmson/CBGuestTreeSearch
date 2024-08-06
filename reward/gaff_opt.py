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

class Molecule():
    """ Stores the rdkit mol (which stores properties) as all files relating to the molecule
    """

    def __init__(self, mol, name, chgfile=None, **kwargs):
        # Currently fileinput must be an sdf
        self.mol = mol
        self.name = name
        self.residue_name = name.upper()[:3]
        # Dictionary of file format, file location
        self.files = {}

        # Add a partial charge file if available
        if chgfile:
            self.files["chgfile"] = chgfile
        
        # Remaining kwargs define output directories for calculations
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
        self.outdir = conf["output_dir"]
        self.hostdir = conf["host_dir"]
        self.host_sdf = conf["host_sdf"]

        # Get host energy
        hostmol = Chem.MolFromMolFile(conf["host_sdf"], removeHs=False)
        self.host = Molecule(hostmol, "cb7")
        print("optimising host")
        hostmol, host_en = self.get_opt(self.host, conf["host_dir"])
        self.hostmol = hostmol
        self.host_en = host_en
        
    # For gaff, the routine needs to keep the tempfiles of the guest to have the prmtop for the complex
    def get_guest_complex_opt(self, comp, *guests):
        """ optimises the guest, keeping the prmtop for the complex, then optimises the complex
        """
        with tempfile.TemporaryDirectory(dir=os.path.abspath(f"{self.outdir}")) as d:
            orgdir = os.getcwd()
            os.chdir(d)

            finalguestmols = []
            finalguestens = []
            guestmolecules = []
            # Get opt for the guests
            for count, i in enumerate(guests):
                guest = Molecule(i, f"guest_{count}")
                guestmolecules.append(guest)
                finalguestmol, guest_en = self.get_opt(guest, d)
                finalguestmols.append(finalguestmol)
                finalguestens.append(finalguestens)

            complexmolecule = Molecule(comp, "complex")
            setattr(complexmolecule, "constituents", [self.host, *guestmolecules])

            finalcomplexmol, complex_en = self.get_opt(complexmolecule, d)
            os.chdir(orgdir)
            # A breakpoint can be used to check temporary directory before it is deleted for debugging
            breakpoint()

        return finalcomplexmol, finalguestmols
            
    def get_opt(self, molecule, min_outdir):
        """ Calls methods to optimise mol and retrieve energy
        """
        orgdir = os.getcwd()
        os.chdir(min_outdir)
        Chem.MolToMolFile(molecule.mol, f"{molecule.name}.sdf", kekulize=False)

        self.molecule_preprocess(molecule, preprocessing_outdir=os.path.abspath(min_outdir))
        self.conj_newt_nab(molecule, min_outdir=os.path.abspath(min_outdir))
        self.get_en_conj_newt(molecule, outfile=os.path.abspath(f"{molecule.files['min_out']}"))

        en = molecule.en_dict["en"]
        
        # NEED TO HANDLE CHARGE HERE
        try:
            sp.run(["obabel", "-ipdb", molecule.files["min_pdb"], "-osdf", "-O", f"{molecule.name}_min.sdf"], stdout=sp.DEVNULL, stderr=sp.DEVNULL)
            finalmol = Chem.MolFromMolFile(f"{molecule.name}_min.sdf", removeHs=False, sanitize=False)
        except:
            os.chdir(orgdir)
            raise ValueError

        for key, val in molecule.en_dict.items():
            finalmol.SetProp(f"gaff_{key}", str(val))

        os.chdir(orgdir)

        return finalmol, en

    @staticmethod
    def get_pdb(molecule, preprocessing_outdir):
        #inp = f"{molecule.fileinput}"
        molecule.files["pdb"] = f"{preprocessing_outdir}/{molecule.name}.pdb"
        molecule.files["pdb4amber"] = f"{preprocessing_outdir}/{molecule.name}_4amber.pdb"

        Chem.MolToMolFile(molecule.mol, f"{molecule.name}.sdf")
        sp.run(["obabel", "-isdf", f"{molecule.name}.sdf", "-opdb", "-O", molecule.files["pdb"]], stdout=sp.DEVNULL, stderr=sp.DEVNULL)
        sp.run(["pdb4amber", "-i", molecule.files["pdb"], "-o", molecule.files["pdb4amber"]], stdout=sp.DEVNULL, stderr=sp.DEVNULL)

        # Need to edit residue names if the molecule is a complex
        if hasattr(molecule, "constituents"):
            edit_resname(molecule.files["pdb4amber"], *molecule.constituents)
        else:
            edit_resname(molecule.files["pdb4amber"], molecule)

    @staticmethod
    def sovlate_packmol_memgen(molecule, preprocessing_outdir):
        pass

    @staticmethod
    def antechamber(molecule, preprocessing_outdir, antechamber_conf):
        molecule.files["mol2"] = f"{preprocessing_outdir}/{molecule.name}.mol2"
               
        ante_cmd = ["antechamber", "-i", molecule.files["pdb4amber"], "-fi", "pdb", "-fo", "mol2", "-o", molecule.files["mol2"]]
            
        if antechamber_conf["chg_method"] == "read":
            ante_cmd.extend(["-c", "rc", "-cf", molecule.files["chgfile"]]) 
            
        else:
            chg_method = antechamber_conf["chg_method"]
            ante_cmd.extend(["-c", chg_method])

        ante_cmd.extend(["-s", "2", "-pf", "yes", "-rn", molecule.residue_name, "-at", "gaff2"])

        sp.run(ante_cmd, stdout=sp.DEVNULL)

    @staticmethod
    def parmchk(molecule, preprocessing_outdir):
        molecule.files["frcmod"] = f"{preprocessing_outdir}/{molecule.name}.frcmod"
    
        sp.run(["parmchk2", "-i", molecule.files["mol2"], "-f", "mol2", "-o", molecule.files["frcmod"]], stdout=sp.DEVNULL)

    @staticmethod
    def tleap_single(molecule, preprocessing_outdir):
        molecule.files["lib"] = f"{preprocessing_outdir}/{molecule.name}.lib"
        molecule.files["prmtop"] = f"{preprocessing_outdir}/{molecule.name}.prmtop"
        molecule.files["rst7"] = f"{preprocessing_outdir}/{molecule.name}.rst7"
               
        tleap_script = ["source leaprc.gaff2",
                        f"loadamberparams {molecule.files['frcmod']}",
                        f"{molecule.residue_name} = loadmol2 {molecule.files['mol2']}",
                        f"check {molecule.residue_name}",
                        f"saveoff {molecule.residue_name} {molecule.files['lib']}",
                        f"saveamberparm {molecule.residue_name} {molecule.files['prmtop']} {molecule.files['rst7']}",
                        "quit"
                        ]

        with open(f"{preprocessing_outdir}/leap.in", "w") as fw:
            for i in tleap_script:
               fw.write(f"{i}\n")

        sp.run(["tleap", "-f", f"{preprocessing_outdir}/leap.in"], stdout=sp.DEVNULL)
        sp.run(["mv", "leap.log", f"{preprocessing_outdir}"])

    @staticmethod
    def tleap_complex(complexmolecule, preprocessing_outdir, *constituents):
        # Instead of creating a lib file, load in the lib files from the individual constituents
        # These should be in the molecule.files dictionary of the complex molecule object
        complexmolecule.files["mol2"] = f"{preprocessing_outdir}/{complexmolecule.name}.mol2"
        complexmolecule.files["prmtop"] = f"{preprocessing_outdir}/{complexmolecule.name}.prmtop"
        complexmolecule.files["rst7"] = f"{preprocessing_outdir}/{complexmolecule.name}.rst7"
        complexmolecule.residue_name = "COM"

        tleap_script = ["source leaprc.gaff2"]

        tleap_script.extend([f"loadoff {molecule.files['lib']}" for molecule in constituents])
        tleap_script.extend([f"loadamberparams {molecule.files['frcmod']}" for molecule in constituents])

        tleap_script.extend([f"{complexmolecule.residue_name} = loadPDB {complexmolecule.files['pdb4amber']}",
                             f"savemol2 {complexmolecule.residue_name} {complexmolecule.files['mol2']} 1",
                             f"saveamberparm {complexmolecule.residue_name} {complexmolecule.files['prmtop']} {complexmolecule.files['rst7']}",
                             "quit"
        ])

        with open(f"{preprocessing_outdir}/leap.in", "w") as fw:
            for i in tleap_script:
               fw.write(f"{i}\n")

        sp.run(["tleap", "-f", f"{preprocessing_outdir}/leap.in"], stdout=sp.DEVNULL)
        sp.run(["mv", "leap.log", f"{preprocessing_outdir}"])

    @staticmethod
    def minimize_sander(molecule, min_outdir, sander_file):
        molecule.files["ncrst"] = f"{min_outdir}/{molecule.name}.ncrst"
        molecule.files["min_out_sander"] = f"{min_outdir}/mdinfo"
        molecule.files["min_log"] = f"{min_outdir}/{molecule.name}_sander_min.log"

        sp.run(["sander", "-O",
                "-i", sander_file,
                "-o", molecule.files['min_log'],
                "-p", molecule.files['prmtop'],
                "-c", molecule.files['rst7'],
                "-r", molecule.files['ncrst']],
               stdout=open(f"{min_outdir}/min.err", "w")
        )
        sp.run(["mv", "mdinfo", f"{min_outdir}"])

    @staticmethod
    def minimize_nab(molecule, min_outdir):
        molecule.files["ncrst"] = f"{min_outdir}/{molecule.name}.ncrst"
        molecule.files["min_out"] = f"{min_outdir}/{molecule.name}_min.log"

        # NOTE the outfiles with have nothing in them if there is a segmentation fault during nab script execution
        sp.run(["./nab_scripts/xmin",
                molecule.files["prmtop"],
                molecule.files["rst7"],
                molecule.files["ncrst"]],
               stdout=open(molecule.files["min_out"], "w"),
               stderr=open(molecule.files["min_out"], "w")
        )

    @staticmethod
    def nmode_nab(molecule, nmode_outdir):
        # Calculates normal modes from the nmode nab script
        molecule.files["nmode_out"] = f"{nmode_outdir}/{molecule.name}_nmode.log"

        sp.run(["./nab_scripts/nmode",
                molecule.files["prmtop"],
                molecule.files["ncrst"]],
               stdout=open(molecule.files["nmode_out"], "w"),
               stderr=open(molecule.files["nmode_out"], "w")
        )

        # Cleanup of undeeded file
        sp.run(["rm", "vecs"])

    @staticmethod
    def conj_newt_nab(molecule, min_outdir):
        # Runs conjugate gradient descent minimisation, followed by newton raphson minimisation, followed by nmode
        # Uses a nab script, requires original nab compiler. Am planning to get a nabc script working if possible
        molecule.files["min_out"] = f"{min_outdir}/{molecule.name}_min.log"
        molecule.files["min_pdb"] = f"{min_outdir}/{molecule.name}_min.pdb"
        
        nab_script = [
            'molecule m;',
            'float x[dynamic], f[dynamic], fret;',
        
            'm = getpdb( "{0}" );',
            'readparm( m, "{1}");',
            'allocate x[ 3*m.natoms ];  allocate f[ 3*m.natoms ];'
            'setxyz_from_mol( m, NULL, x );',

            'mm_options( "ntpr=50, nsnb=99999, gb=1, gbsa=1, cut=999.0, rgbmax=999.0, diel=C, dielc=1.0" );',
            'mme_init( m, NULL, "::Z", x, NULL);',
        
            '// conjugate gradient minimization',
            # 'dgrad = 3*natoms*0.001;',
            'conjgrad(x, 3*m.natoms, fret, mme, 0.00001, 0.01, 100000 );',
        
            '// Newton-Raphson minimization\fP',
            'mm_options( "ntpr=1" );',
            'newton( x, 3*m.natoms, fret, mme, mme2, 0.0000001, 0.0, 100 );',

            '// Write out the minimized pdb',
            'putpdb( "{2}", m );',
        
            '// get the normal modes:',
            'nmode( x, 3*m.natoms, mme2, 0, 0, 0.0, 0.0, 0);'
        ]
        
        nab_script = [x.format(molecule.files["pdb4amber"], molecule.files["prmtop"], molecule.files["min_pdb"])+'\n' for x in nab_script]
        
        # Generate nab script (for some reason have to hard code argv, another thing to try and fix
        with open(f"{min_outdir}/min_prog.nab", 'w') as w:
            w.writelines(nab_script)
        
        # Compile nab script (ensure the original nab compiler is in path)
        sp.run(["nab", f"{min_outdir}/min_prog.nab", "-o", f"{min_outdir}/min_prog"])

        # Run nab script
        sp.run([f"{min_outdir}/min_prog"],
               stdout=open(molecule.files["min_out"], "w"),
               stderr=open(molecule.files["min_out"], "w")
        )
        
        return

    # CLASS METHODS FOR AMBER CALCULATIONS
    def molecule_preprocess(self, molecule, preprocessing_outdir):

        # Minimsise single molecule
        if len(Chem.GetMolFrags(molecule.mol, asMols=True)) == 1:
            AmberCalculator.get_pdb(molecule, preprocessing_outdir)
            AmberCalculator.antechamber(molecule, preprocessing_outdir, self.conf["antechamber_conf"])
            AmberCalculator.parmchk(molecule, preprocessing_outdir)
            AmberCalculator.tleap_single(molecule, preprocessing_outdir)

        # Minmise a complex of moecules that are not bound covalently
        if len(Chem.GetMolFrags(molecule.mol, asMols=True)) > 1:
            AmberCalculator.get_pdb(molecule, preprocessing_outdir)
            AmberCalculator.tleap_complex(molecule, preprocessing_outdir, *molecule.constituents)

    # CLASS METHODS FOR RETRIEVAL OF OUTPUTS
    def get_opt_pdb(self, molecule, min_outdir):
        molecule.files["min_pdb"] = f"{min_outdir}/{molecule.name}_min.pdb"
        
        try:
            sp.run(["ambpdb", "-p", molecule.files['prmtop'], "-c", molecule.files["ncrst"]], stdout=open(molecule.files["min_pdb"], "w"))
        except FileNotFoundError:
            print("no minimisation output file")

    def get_en_sander(self, molecule, outfile):

        min_out = [i.split() for i in open(outfile, "r")]
        en_dict = {}

        en_dict["Etot"] = float(min_out[3][1])
        en_dict["Ebond"] = float(min_out[5][2]) # bond stretching
        en_dict["Eangle"] = float(min_out[5][5]) # angle bending
        en_dict["Edihed"] = float(min_out[5][8]) # intrinsic dihedral
        en_dict["Evdw"] = float(min_out[6][2]) # van der waals
        en_dict["Ecoul"] = float(min_out[6][5]) # coloumb energy
        en_dict["Egb"] = float(min_out[6][8]) # Generalised born energy
        en_dict["Evdw_14"] = float(min_out[7][3]) # topological 1-4 van der waals (atoms separated by 3 bonds)
        en_dict["Ecoul_14"] = float(min_out[7][7]) # topological 1-4 coulomb (atoms separated by 3 bonds)
        en_dict["Enp"] = float(min_out[8][2]) # Non polar solvent energy, associate with surface area in contact with solvent

        setattr(molecule, "en_dict", en_dict)

    def get_en_nab(self, molecule, outfile):

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
                en_dict["Etot"] = float(i[1])
                
        setattr(molecule, "en_dict", en_dict)

        return

    def get_en_conj_newt(self, molecule, outfile):
        # Gets energy from the conjugate gradient + newton raphson routine
        # Also retrieves normal modes
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
            extra_frags = GetMolFrags(molecule.mol, asMols=True) - 1
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

        setattr(molecule, "en_dict", en_dict)
            
        return

    def get_nmode(self, molecule, outfile):
        # Retrieves normal modes from the nmode calculation and calculates configurational entropy contribution to energy
        min_out = [i.split() for i in open(outfile, "r")]
        
        min_out.remove([])

        E_zpve = 0
        E_entropic = 0
        for i in min_out:
            try:
                if i[0] == "Zero-point":
                    E_zpve = float(i[3])
                if i[0] == "Total:":
                    E_entropic = float(i[1])
            except:
                continue

        molecule.en_dict["Etot_corr"] = molecule.en_dict["Etot"] - (E_entropic - E_zpve)

        return

if __name__ == "__main__":

    # May wish to change antechamber conf to file input as with sander
    gaff_conf = {"antechamber_conf" : {"chg_method" : "gas"},
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
                gaff_conf["antechamber_conf"] = {"chg_method" : "bcc"}
                
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
