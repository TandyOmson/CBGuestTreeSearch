""" Filter out molecules containing triple bonds
    This has come about when designing hydrocarbons,
    often triple bonds or strange double + 2 bonds are generated
    This results in strong, but unrealistic interactions with the host
"""

from filter.filter import Filter
from rdkit import Chem

class TripleBondFilter(Filter):
    def check(mol, config):
        num_triple_bonds = 0
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
                num_triple_bonds += 1
        return num_triple_bonds == 0