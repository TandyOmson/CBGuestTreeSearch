from rdkit import Chem
import sys

mol = Chem.MolFromMolFile(sys.argv[1], removeHs=False)
host, guest = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)

Chem.MolToMolFile(host, "host.sdf")
Chem.MolToMolFile(guest, "guest.sdf")
