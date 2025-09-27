import numpy as np
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from scipy.stats import norm
from skdim.id import MLE, KNN

def twonn(X):
    knn = KNN(k=2)
    return knn.fit_transform(X)

def mle(X):
    """MLE intrinsic dimensionality estimator from scikit-dimension
    """
    mle = MLE()
    return mle.fit_transform(X)
