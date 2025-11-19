from rdkit import Chem
from rdkit.Chem import AllChem
from joblib import Parallel, delayed
import pandas as pd
import sys
import os
import subprocess as sp

def robust_embed(mol):
    if mol == None:
        raise Exception
    
    mol = Chem.AddHs(mol)

    res = AllChem.EmbedMolecule(mol)
    if res == -1:
        res = AllChem.EmbedMolecule(mol, useBasicKnowledge=False)
        if res == -1:
            raise Exception
        else:
            AllChem.MMFFOptimizeMolecule(mol)

    return mol

def openbabel_embed(smi):
    obabel_cmd = 'echo "{}" | obabel -ismi -osdf --gen3d best --minimize'
    res = sp.run(obabel_cmd.format(smi),
                 shell=True,
                 capture_output=True,
                 text=True,
                 )
    
    mol = Chem.MolFromMolBlock(res.stdout, removeHs=False)
    if mol == None:
        print(smi)
        print(res.stderr)
        raise Exception

    return mol

def process_smi(smi, idx, nullmol, outdir, write_individual=False):
    mol = Chem.MolFromSmiles(smi)
    
    if mol is None:
        mol = nullmol
    try:
        mol = robust_embed(mol)
    except:
        try:
            mol = openbabel_embed(smi)
        except:
            mol = nullmol
        
    if write_individual:
        Chem.MolToMolFile(mol, os.path.join(outdir, f"mol_{idx}.sdf"))
    return mol

if __name__ == "__main__":
    smis = [i.rstrip() for i in open(sys.argv[1], "r").readlines()]
    outdir = sys.argv[2]

    nullmol = Chem.MolFromSmiles("C")
    
# Parallel execution
    results = Parallel(n_jobs=-1, verbose=5)(
        delayed(process_smi)(smi, idx, nullmol, outdir)
        for idx, smi in enumerate(smis)
    )

    # Collect into single SDF
    writer = Chem.SDWriter(os.path.join(outdir, "all.sdf"))
    for count, mol in enumerate(results):
        writer.write(mol)           
    writer.close()
