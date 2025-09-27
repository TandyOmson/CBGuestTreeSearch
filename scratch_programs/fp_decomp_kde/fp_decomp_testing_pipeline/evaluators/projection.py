import numpy as np
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import pairwise_distances, adjusted_mutual_info_score, jaccard_score

def fp_decomp_distance_correlation(fps, X):
    """ Gets correlations between pairwise tanimoto distances of fps
        and euclidian distances in decomposition space
    """
    # jaccard is approx equal to tanimoto for bit vectors
    fp_tanimoto = 1 - pairwise_distances(fps, metric='jaccard')
    decomp = pairwise_distances(X, metric="euclidean")
    spearman = spearmanr(fp_tanimoto.flatten(), decomp.flatten())
    pearson = pearsonr(fp_tanimoto.flatten(), decomp.flatten())
    return float(spearman[0]),float(spearman[1]), float(pearson[0]), float(pearson[1])

def stability_under_noise_perturbation(fps, X, pipe):
    """ Adds noise, checks mean projection shift
    """
    noisy_fps = fps.copy()
    noise_mask = np.random.binomial(1, 0.9, fps.shape)
    noisy_fps *= noise_mask
    noisy_X = pipe.fit_transform(noisy_fps)

    stability = np.mean(np.linalg.norm(X - noisy_X, axis=1))
    return stability

def outlier_mean_shift(X, kde):
    """ Removes outliers, checks mean projection shift
    """
    density = kde.score_samples(X)
    non_outliers_idxs = density > np.percentile(density, 5)
    X_non_outliers = X[non_outliers_idxs]

    mean_shift = np.linalg.norm(np.mean(X, axis=0) - np.mean(X_non_outliers, axis=0))
    return mean_shift

def split_mean_shift(X, X_val_idxs):
    """ Checks difference in mean projection between X and X without X_val_idxs
    """
    mean_shift = np.linalg.norm(np.mean(X, axis=0) - np.mean(X[~np.isin(np.arange(X.shape[0]), X_val_idxs)], axis=0))
    return mean_shift
