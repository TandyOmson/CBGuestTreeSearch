""" Gets binding for the dataset 
1. Move complex sdf to new folder, create guest sdf in the folder
2. Optimise complex and guest
3. Calculate binding
4. Add binding to dataframe

Usage:
python xtb_bind.py <complex_sdf_folder> <host_xtb_opthess.out>
"""

from rdkit import Chem
import subprocess as sp
from glob import glob
import os
import sys

def xtb_opt(filename, outfile):
    """ Optimises from an sdf file """

    sp.run(["obabel","-isdf","-osdf",f"{filename}","-O",f"{filename}"])
    sp.run(["xtb",f"{filename}","--ohess","normal","--alpb","water"],stdout=open(outfile,"w"),stderr=sp.DEVNULL)
    
    return outfile

def get_binding(filename):
    """ Gets binding energy from an sdf file """

    gen = (i.split() for i in reversed(open(filename,"r").readlines()))

    binding = 0
    for i in gen:
        if i:
            if " ".join(i[:4]) == ":: total free energy":
                binding = float(i[4])*627.5095
                break

    return binding

