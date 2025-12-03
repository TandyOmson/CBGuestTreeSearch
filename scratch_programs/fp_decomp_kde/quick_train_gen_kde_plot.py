""" Program for quickly plotting a training set and generated molecules in kde space, given a decomposition pipeline and kde object
"""

import pickle
import argparse

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
import numpy as np

import matplotlib.pyplot as plt

def decomp_smis(smis, pipe):
    mols = [Chem.MolFromSmiles(s) for s in smis]
    
    fpgen = rdFingerprintGenerator.GetAtomPairGenerator()
    fps = [fpgen.GetFingerprint(m) for m in mols]
    
    X = pipe.transform(fps)

    return X

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="quick plot for training and generated smiles in kde space",
        usage="python {os.path.basename(__file__)} -p PIPE_KDE_DIRECTORY -t TRAINING_SMIS -g GEN_SMIS"
    )
    parser.add_argument(
        "-p",
        "--pipedir",
        type=str,
        required=True,
        help="path to directory with pipe.pkl and kde.pkl in it",
    )
    parser.add_argument(
        "-t",
        "--trainfile",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-g",
        "--genfile",
        type=str,
        required=True,
    )
    args = parser.parse_args()
    
    pipefile = args.pipedir + "/" + "pipe.pkl"
    kdefile = args.pipedir + "/" + "kde.pkl"
    
    with open(pipefile, 'rb') as fr:
        pipe = pickle.load(fr)

    with open(kdefile, 'rb') as fr:
        kde = pickle.load(fr)

    train_smis = [i.rstrip() for i in open(args.trainfile, 'r').readlines()]
    print("Processing training set...")
    X_train = decomp_smis(train_smis, pipe)

    gen_smis = [i.rstrip() for i in open(args.genfile, 'r').readlines()]
    print("Processing generated set...")
    X_gen = decomp_smis(gen_smis, pipe)

    fig, ax = plt.subplots()
    
    xx, yy = np.meshgrid(np.linspace(-1.1, 1.1, 100), np.linspace(-1.1, 1.1, 100))
    Z = np.vstack([yy.ravel(), xx.ravel()]).T
    Z = np.exp(kde.score_samples(Z)).reshape(100, 100)

    levels = np.linspace(0, Z.max(), 25)
    m = ax.contour(xx, yy, Z, levels=levels, cmap=plt.cm.Reds)

    fig = ax.figure
    cbar = fig.colorbar(m, ax=ax)
    cbar.set_label('Gaussian Kernel Density')

    ax.scatter(X_train[:, 0], X_train[:, 1], s=2, c='blue', label='training')
    ax.scatter(X_gen[:, 0], X_gen[:, 1], s=2, c='orange', label='generated')
    ax.legend()
    ax.grid(True)

    plt.show()
