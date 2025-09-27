from typing import Callable, Dict
import pickle
import pandas as pd
import numpy as np
import random

from rdkit import Chem
from rdkit import DataStructs
from sklearn.pipeline import Pipeline
from sklearn.neighbors import KernelDensity
from rdkit.ML.Cluster import Butina

from sklearn.decomposition import PCA, KernelPCA
from sklearn.manifold import LocallyLinearEmbedding, Isomap, TSNE
import umap
from sklearn.preprocessing import MinMaxScaler

from fingerprints.fingerprints import FingerprintGenerators
from utils import cluster_by_kde_contours, sample_by_cluster, retrain_decomper, butina_thresh_search
from evaluators.gaussianity import skew, kurtosis, kl_divergence
from evaluators.projection import fp_decomp_distance_correlation, stability_under_noise_perturbation, split_mean_shift, outlier_mean_shift
from evaluators.clusters import cluster_overlap, get_cluster_sizes
from evaluators.intrinsic_dim import twonn, mle
from evaluators.plots import plot_decomp_kde_2D, plot_decomp_kde_2D_val, plot_decomp_kde_2D_clustered

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec

debug = False
def dprint(s):
    if debug == True:
        print(s)
    else:
        pass

fpgen = FingerprintGenerators()

FINGERPRINTS: Dict[str, Callable] = {
    'morgan_ECFP4_2048' : lambda mol: fpgen.get_morgan_ecfp4(mol),
    'morgan_ECFP6_2048' : lambda mol: fpgen.get_morgan_ecfp6(mol),
    'atom_pair_2048' : lambda mol: fpgen.get_atom_pair_fp(mol),
    'topo_tors_2048' : lambda mol: fpgen.get_torsion_fp(mol),
    'rdkit_2048' : lambda mol: fpgen.get_rdkit_fp(mol),
    'avalon_512' : lambda mol: fpgen.get_avalon_fp(mol),
    'MACCS_166' : lambda mol: fpgen.get_maccs_fp(mol),
}

DECOMPOSERS: Dict[str, Pipeline] = {
    "pca" : Pipeline([("decomp", PCA(n_components=2)), ("scaler", MinMaxScaler((-1,1)))]),
    "pca_rbf" : Pipeline([("decomp", KernelPCA(n_components=2, kernel="rbf")),("scaler", MinMaxScaler((-1,1)))]),
    "isomap" : Pipeline([("decomp", Isomap(n_components=2)),("scaler", MinMaxScaler((-1,1)))]),
    "umap" : Pipeline([("decomp", umap.UMAP(n_components=2)),("scaler", MinMaxScaler((-1,1)))]),
    "lle" : Pipeline([("decomp", LocallyLinearEmbedding(n_components=2)),("scaler", MinMaxScaler((-1,1)))]),
    "t-sne" : Pipeline([("decomp", TSNE(n_components=2)),("scaler", MinMaxScaler((-1,1)))]),
}

# Load SMILES data
smis = [i.rstrip() for i in open("data/hydros8992.smi", "r").readlines()][:2000]
mols = [Chem.MolFromSmiles(s) for s in smis]

results = []
# Run combinations
run_count = 1
for fp_name, fp_func in FINGERPRINTS.items():
    fps_all_vect = [fp_func(m) for m in mols]
    fps_all = np.array(fps_all_vect)

    # avoid running more than once per fingerprint
    dprint("TwoNN")
    two_nearest_neighbour_estimate = twonn(fps_all)
    dprint("MLE")
    maximum_likelihood_estimate = mle(fps_all)
    
    for decomp_name, decomp_func in DECOMPOSERS.items():
        print("Running {} / {}, ({} of {})".format(fp_name, decomp_name, run_count, len(FINGERPRINTS)*len(DECOMPOSERS)))

        # Generate data for tests
        try:
            X_decomp_all = decomp_func.fit_transform(fps_all)
        except:
            print("skipped {} / {}".format(fp_name, decomp_name))
            run_count += 1
            continue
        print("decomposed training set to embedding: {}".format(X_decomp_all.shape))
        
        # Gaussian kernel denisty estimation
        kde = KernelDensity(kernel="gaussian", bandwidth=0.5)
        kde.fit(X_decomp_all)

        # Cluster by kde
        dprint("kde clustering")
        n_clusters = 10
        Z_all, contours, kde_cluster_labels, kde_cluster_lens = cluster_by_kde_contours(X_decomp_all, kde, 10, min_cluster_len=1)

        # Cluster by butina
        dprint("butina clustering")
        dist_mat = []
        for i in range(X_decomp_all.shape[0]):
            for j in range(i):
                sim = DataStructs.TanimotoSimilarity(fps_all_vect[i], fps_all_vect[j])
                dist_mat.append(1-sim)

        butina_clusters = Butina.ClusterData(data=dist_mat, nPts=len(fps_all), distThresh=.65, isDistData=True)
        #butina_clusters = butina_thresh_search(dist_mat, len(fps_all), 12, max_iter=5, tol=4, max_cluster_frac=0.5, max_singleton_frac=0.2)
        
        butina_clusters = sorted(butina_clusters, key=len, reverse=True)
        butina_cluster_labels = np.full((X_decomp_all.shape[0]), -1)
        for cluster_idx in range(len(butina_clusters)):
            for idx in butina_clusters[cluster_idx]:
                butina_cluster_labels[idx] = cluster_idx + 1

        butina_cluster_lens = get_cluster_sizes(butina_cluster_labels)

        # Get samples:
        # random sample
        sample_idxs_random = random.sample(range(X_decomp_all.shape[0]), int(X_decomp_all.shape[0]*0.1))
        fps_sample_random = fps_all[sample_idxs_random, :]
        X_decomp_sample_random = X_decomp_all[sample_idxs_random, :]

        # stratified by kde clusters sample
        sample_idxs_kde = sample_by_cluster(kde_cluster_labels, val_split=0.2)

        # stratified by butina clusters sample
        sample_idxs_butina = sample_by_cluster(butina_cluster_labels, val_split=0.2)
        # TODO: TEMP function
        if len(sample_idxs_kde) == 0:
            print("skipping kde for {} / {}".format(fp_name, decomp_name))
            sample_idxs_kde = sample_idxs_random
        
        if len(sample_idxs_butina) == 0:
            print("skipping butina for {} / {}".format(fp_name, decomp_name))
            sample_idxs_butina = sample_idxs_random

        # Evaluations
        metrics = {}

        # gaussianity
        dprint("skew")
        metrics["skew_1"], metrics["skew_2"] = skew(X_decomp_all)
        dprint("kurtosis")
        metrics["kurtosis_1"], metrics["kurtosis_2"] = kurtosis(X_decomp_all)
        dprint("kl div")
        metrics["kl_div"] = kl_divergence(X_decomp_all)

        # projection
        dprint("corr (sample)")
        metrics["spearmanr"], metrics["spearmanp"], metrics["pearsonr"], metrics["pearsonp"] = fp_decomp_distance_correlation(fps_sample_random, X_decomp_sample_random)
        dprint("noise pert")
        metrics["noise_pert"] = stability_under_noise_perturbation(fps_all, X_decomp_all, decomp_func)
        dprint("outlier mean shift")
        metrics["outlier_mean_shift"] = outlier_mean_shift(X_decomp_all, kde)

        # shifts in training sets (removing validation sets)
        metrics["kde_mean_shift"] = split_mean_shift(X_decomp_all, sample_idxs_kde)
        metrics["butina_mean_shift"] = split_mean_shift(X_decomp_all, sample_idxs_butina)

        # clusters
        metrics["avg_cluster_size"] = np.mean(kde_cluster_lens)
        metrics["max_cluster_size"] = max(kde_cluster_lens)
        metrics["min_cluster_size"] = min(kde_cluster_lens)

        metrics["cluster_overlap"] = cluster_overlap(kde_cluster_labels, butina_cluster_labels)

        # intrinsic dimension estimation (calculated earlier for the fingerprint)
        metrics["twoNN"] = two_nearest_neighbour_estimate
        metrics["MLE"] = maximum_likelihood_estimate

        results.append({'fingerprint' : fp_name,
                        'decomp' : decomp_name,
                        **metrics,
                        })
        
        # Retrain the pipeline omitting a random sample, kde cluster sample and butina cluster sample
        # Transform training and val set using fitted pipeline, then plot
        X_train_random, X_val_random, kde_random = retrain_decomper(fps_all, sample_idxs_random, decomp_func)
        X_train_kde, X_val_kde, kde_kde = retrain_decomper(fps_all, sample_idxs_kde, decomp_func)
        X_train_butina, X_val_butina, kde_butina = retrain_decomper(fps_all, sample_idxs_butina, decomp_func)
        
        # Visualisations
        # 3x2 grid
        # no clustering , kde clustering, butina clustering
        # random val set, kde val set   , butina val set

        # Make overarching figure
        fig = plt.figure(figsize=(12, 8), constrained_layout=True)
        fig.suptitle(fp_name + " / " + decomp_name, size="xx-large")
        gs = GridSpec(2, 3, figure=fig,
                      height_ratios=[1,1],
                      width_ratios=[1,1,1])

        ax_1 = fig.add_subplot(gs[0, 0])
        plot_decomp_kde_2D(X_decomp_all, Z_all, ax_1)

        ax_2 = fig.add_subplot(gs[0, 1])
        plot_decomp_kde_2D_clustered(X_decomp_all, Z_all, ax_2, kde_cluster_labels, kde_cluster_lens)

        ax_3 = fig.add_subplot(gs[0, 2])
        plot_decomp_kde_2D_clustered(X_decomp_all, Z_all, ax_3, butina_cluster_labels, butina_cluster_lens)
        ax_4 = fig.add_subplot(gs[1,0])
        plot_decomp_kde_2D_val(X_train_random, kde_random, ax_4, X_val_random)
        ax_4.set_title("Random sample")

        ax_5 = fig.add_subplot(gs[1,1])
        plot_decomp_kde_2D_val(X_train_kde, kde_kde, ax_5, X_val_kde)
        ax_5.set_title("KDE Clusters")

        ax_6 = fig.add_subplot(gs[1,2])
        plot_decomp_kde_2D_val(X_train_butina, kde_butina, ax_6, X_val_butina)
        ax_6.set_title("Butina Clusters")

        fig.savefig("results/" + fp_name + "_" + decomp_name + ".png", dpi=300)
        plt.close(fig)

        run_count += 1
        df_results = pd.DataFrame(results)
        df_results.to_csv("temp_results.csv")
        
df_results = pd.DataFrame(results)
df_results.to_csv("results/metrics.csv")
