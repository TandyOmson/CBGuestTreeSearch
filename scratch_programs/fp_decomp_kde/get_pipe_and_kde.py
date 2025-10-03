""" Minimal code for getting mol embedding pipe and kde from a SMILES set
"""
import sys
import pickle

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.manifold import Isomap
from sklearn.preprocessing import MinMaxScaler
from sklearn.pipeline import Pipeline
from sklearn.neighbors import KernelDensity

smis = [i.rstrip() for i in open(sys.argv[1], "r").readlines()]
mols = [Chem.MolFromSmiles(s) for s in smis]
fpgen = rdFingerprintGenerator.GetAtomPairGenerator()
fps = [fpgen.GetFingerprint(m) for m in mols]

reducer = Isomap(n_components=2, n_neighbors=10, n_jobs=8)
scaler = MinMaxScaler((-1,1))
pipe = Pipeline([("reducer", reducer), ("scaler", scaler)])
X = pipe.fit_transform(fps)

kde = KernelDensity(kernel="gaussian", bandwidth=0.5)
kde.fit(X)

outdir = sys.argv[2]

with open(outdir + "pipe.pkl", "wb") as fw:
    pickle.dump(pipe, fw)
    
with open(outdir + "kde.pkl", "wb") as fw:
    pickle.dump(kde, fw)
