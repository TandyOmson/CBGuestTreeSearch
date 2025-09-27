import matplotlib.pyplot as plt
from matplotlib.path import Path
import numpy as np
import random
from sklearn.neighbors import KernelDensity
from rdkit.ML.Cluster import Butina

def cluster_by_kde_contours(X_decomp, kde, n_groups, min_cluster_len=100):

    X, Y = np.meshgrid(
        np.linspace(-1.01, 1.01, 100), 
        np.linspace(-1.01, 1.01, 100)
    )
    Z = np.vstack([Y.ravel(), X.ravel()]).T
    # density at gridpoints
    Z = np.exp(kde.score_samples(Z)).reshape(X.shape)

    levels = np.linspace(Z.min(), Z.max(), n_groups+2)

    fig, ax = plt.subplots()
    contours = ax.contourf(X, Y, Z, levels=levels)
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

def sample_by_cluster(cluster_labels, val_split=0.1):
    n_clusters = len(np.unique(cluster_labels))
    val_set_idxs = []
    empty_clusters = 0
    for i in range(1, n_clusters+1):
        if len(np.where(cluster_labels==i)[0]) < int((len(cluster_labels)/n_clusters)*val_split):
            empty_clusters += 1

    try:
        samples_per_cluster = int((len(cluster_labels)/(n_clusters-empty_clusters))*val_split)
    except:
        samples_per_cluster = len(cluster_labels)*0.01
    
    
    for i in range(1, n_clusters+1):
        # If not in a cluster, leave out of the validation set
        if len(np.where(cluster_labels==i)[0]) > samples_per_cluster:
            val_set_idxs.extend(random.sample(list(np.where(cluster_labels==i)[0]), samples_per_cluster))
    return val_set_idxs

def retrain_decomper(fps, val_idxs, decomp_func):
    """ Removes X_val_idxs from fps, then retrains decomper
        Returns X_train, X_val and a KernelDensity object trained on X_train
    """
    fps_train = fps[~np.isin(np.arange(fps.shape[0]), val_idxs)]
    X_train = decomp_func.fit_transform(fps_train)
    X_val = decomp_func.transform(fps[val_idxs,:])
        
    # Gaussian kernel denisty estimation
    kde = KernelDensity(kernel="gaussian", bandwidth=0.5)
    kde.fit(X_train)
    
    return X_train, X_val, kde

def butina_thresh_search(distmat, n_pts, target_n_clusters, max_iter=20, tol=0, max_cluster_frac=0.5, max_singleton_frac=0.5):
    """ Best butina threshold
    """
    low, high = 0.0, 1.0
    best_clusters = None
    best_diff = float('inf')
    best_balanced_clusters = None
    best_balanced_diff = float('inf')

    count = 0
    while count < max_iter:
        for thresh in list(reversed(np.arange(0.3, 0.9, 0.1))):
            clusters = Butina.ClusterData(distmat, n_pts, thresh, isDistData=True)
            n = len(clusters)
            diff = abs(n - target_n_clusters)
            
            # Evaluate cluster balance
            sizes = [len(c) for c in clusters]
            largest = max(sizes)
            singleton_count = sum(1 for c in clusters if len(c) == 1)
            singleton_frac = singleton_count / n if n > 0 else 0
            
            is_balanced = (
                largest <= max_cluster_frac * n_pts and
                singleton_frac <= max_singleton_frac
            )
            
            # Track best balanced result
            if is_balanced and diff < best_balanced_diff:
                best_balanced_clusters = clusters
                best_balanced_diff = diff
            
            # Always track best match regardless of balance
            if diff < best_diff:
                best_clusters = clusters
                best_diff = diff
            
            # Adjust search window
            if n > target_n_clusters:
                continue
            elif n < target_n_clusters:
                continue
            else:
                return best_balanced_clusters  # exact match
            
            if best_balanced_diff <= tol:
                return best_balanced_clusters

    # Prefer balanced clusters, fall back to best overall
    return best_balanced_clusters if best_balanced_clusters is not None else best_clusters

