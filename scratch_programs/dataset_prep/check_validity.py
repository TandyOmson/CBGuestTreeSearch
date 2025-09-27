from rdkit import Chem
from rdkit import RDLogger
import sys

RDLogger.DisableLog('rdApp.*')

smis = [i.rstrip() for i in open(sys.argv[1], "r").readlines()]

invalid_mol_count = 0
too_big_count = 0
unique_mol_count = 0
duplicate_mol_count = 0

dupes = []
canon_smiles_map = []
inchi_map = []
for count, smi in enumerate(smis):
    print("SMILES string {} of {}, ({} duplicate, {} invalid, {} too big)".format(count, len(smis), duplicate_mol_count, invalid_mol_count, too_big_count), end="\r")

    mol = Chem.MolFromSmiles(smi)
    if mol == None:
        invalid_mol_count += 1
        continue
    elif mol.GetNumHeavyAtoms() > 20:
        too_big_mol_count += 1    
        continue
    canon_smi = Chem.MolToSmiles(mol, isomericSmiles=True, kekuleSmiles=True, canonical=True)

    # Inchi is much better for finding duplicates than SMILES
    inchi = Chem.MolToInchi(mol)
    if inchi in inchi_map:
        duplicate_mol_count += 1
        dupes.append(mol)
    else:
        unique_mol_count += 1
        canon_smiles_map.append(canon_smi)
        inchi_map.append(inchi)

print("retained {}. ({} duplicate, {} invalid, {} too big)".format(len(canon_smiles_map), duplicate_mol_count, invalid_mol_count, too_big_count))

with open(sys.argv[2], "w") as fw:
    for i in canon_smiles_map:
        fw.write(f"{i}\n")

#from rdkit.Chem import AllChem        
#writer = Chem.SDWriter("dupes.sdf")
#for m in dupes:
#    m_H = Chem.AddHs(m)
#    AllChem.EmbedMolecule(m_H)
#    writer.write(m_H)
