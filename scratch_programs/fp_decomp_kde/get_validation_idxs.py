""" Uses kde contours to cluster molceules, which are then used to split a validation set
    Prints validation set indices out
    Plots the clusters and decomposed space points
"""

import pickle
import argparse
import random

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.path import Path

def decomp_smis(smis, pipe):
    mols = [Chem.MolFromSmiles(s) for s in smis]
    
    fpgen = rdFingerprintGenerator.GetAtomPairGenerator()
    fps = [fpgen.GetFingerprint(m) for m in mols]
    
    X = pipe.transform(fps)

    return X

def cluster_by_kde_contours(X_decomp, kde, n_groups, min_cluster_len=100):
    xx, yy = np.meshgrid(
        np.linspace(-1.01, 1.01, 100), 
        np.linspace(-1.01, 1.01, 100)
    )
    Z = np.vstack([yy.ravel(), xx.ravel()]).T
    # density at gridpoints
    Z = np.exp(kde.score_samples(Z)).reshape(xx.shape)

    levels = np.linspace(Z.min(), Z.max(), n_groups+2)

    fig, ax = plt.subplots()
    contours = ax.contourf(xx, yy, Z, levels=levels)
    plt.close(fig)

    # contour polygon containment clustering
    cluster_labels = np.full(X_decomp.shape[0], fill_value=-1)
    for cluster_id, path in enumerate(contours.get_paths()):
            poly = Path(path.vertices)
            inside = poly.contains_points(X_decomp)
            cluster_labels[inside] = cluster_id
            
    # if cluster size is smaller than minimum, combine it with the next cluster along
    smallest_cluster_len = -np.inf
    while smallest_cluster_len < min_cluster_len:
        cluster_lens = []
        for i in np.unique(cluster_labels):
            c_idxs = np.where(cluster_labels==i)[0]
            c_len = len(c_idxs)
            cluster_lens.append(c_len)
            if c_len < min_cluster_len:
                for j in c_idxs:
                    cluster_labels[j] = i+1

        smallest_cluster_len = min(cluster_lens)
                    
    # Alternative fast binning by quantiles (doesn't match contours exactly, but is cheap)
    #thresholds = np.quantile(density, np.linspace(0, 1, n_groups + 1))
    #cluster_labels = np.digitize(density, thresholds[1:], right=True)
    
    return Z, contours, cluster_labels, cluster_lens

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="",
        usage="python {os.path.basename(__file__)} -p PIPE_KDE_DIRECTORY -t TRAINING_SMIS -o VAL_SET_INDICES_OUT"
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
        "-o",
        "--outfile",
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

    val_split = 0.2
    n_clusters = 10

    X = decomp_smis(train_smis, pipe)
    Z, contours, cluster_labels, cluster_lens = cluster_by_kde_contours(X, kde, n_clusters, min_cluster_len=100)

    val_set_idxs = []
    empty_clusters = 0
    for i in range(1, n_clusters+1):
        if len(np.where(cluster_labels==i)[0]) < int((len(cluster_labels)/n_clusters)*val_split):
            empty_clusters += 1

    samples_per_cluster = int((len(cluster_labels)/(n_clusters - empty_clusters))*val_split)

    for i in range(1, n_clusters+1):
        if len(np.where(cluster_labels==i)[0]) > samples_per_cluster:
            val_set_idxs.extend(random.sample(list(np.where(cluster_labels==i)[0]), samples_per_cluster))

    print(f"Training set is size {len(train_smis)}")
    print(f"Validation set is size {len(val_set_idxs)}")

    with open(args.outfile, "w") as fw:
        for i in val_set_idxs:
            fw.write(f"{i}\n")

    fig, ax = plt.subplots(constrained_layout=True)
    
    xx, yy = np.meshgrid(np.linspace(-1.1, 1.1, 100), np.linspace(-1.1, 1.1, 100))
    Z = np.vstack([yy.ravel(), xx.ravel()]).T
    Z = np.exp(kde.score_samples(Z)).reshape(100, 100)

    levels = np.linspace(0, Z.max(), 25)
    m = ax.contour(xx, yy, Z, levels=levels, cmap=plt.cm.Reds)

    fig = ax.figure
    cbar = fig.colorbar(m, ax=ax)
    cbar.set_label('Gaussian Kernel Density')

    ax.scatter(X[:, 0], X[:, 1], s=2, c='blue')
    ax.grid(True)

    inset_ax = ax.inset_axes([0.1, 1.05, 0.8, 0.2])
    inset_ax.bar(np.unique(cluster_labels), cluster_lens, color="blue", alpha=0.7)
    inset_ax.set_ylabel("Cluster Size")

    plt.show()
