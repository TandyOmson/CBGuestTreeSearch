""" Duplication checker
"""

import yaml
import argparse

from logging import getLogger
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolops 
from rdkit import RDLogger
import selfies as sf
import numpy as np
from numpy.random import default_rng
import pickle

from chemtsv2.preprocessing import smi_tokenizer, selfies_tokenizer_from_smiles

if __name__ == "__main__":

    gen_smis = [i.rstrip() for i in open("results/03_12_25_gen.smi", "r").readlines()]
    #check for duplicates
    duplicate_count = 0
    canon_gen_smis = []
    for count, smi in enumerate(gen_smis):
        print(f"Checking {count}", end="\r")
        mol = Chem.MolFromSmiles(smi)
        canon_smi = Chem.MolToSmiles(mol, isomericSmiles=True, kekuleSmiles=True, canonical=True)
        if canon_smi in canon_gen_smis:
            duplicate_count += 1

    print(f"{duplicate_count} duplicate molecules")
