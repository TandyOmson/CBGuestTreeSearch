from sklearn.metrics import adjusted_mutual_info_score
import numpy as np

def cluster_overlap(kde_labels, clust_labels):
    ami = adjusted_mutual_info_score(kde_labels, clust_labels)
    return ami

def get_cluster_sizes(cluster_labels):
    cluster_lens = []
    for i in np.unique(cluster_labels):
        cl = len(np.where(cluster_labels==i)[0])
        cluster_lens.append(cl)

    return cluster_lens
