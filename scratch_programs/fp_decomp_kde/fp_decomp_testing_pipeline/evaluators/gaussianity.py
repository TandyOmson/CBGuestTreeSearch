import numpy as np
from sklearn.neighbors import KernelDensity
from scipy.stats import skew as scipy_skew
from scipy.stats import kurtosis as scipy_kurtosis
from scipy.stats import  multivariate_normal

def skew(X):
    return scipy_skew(X)

def kurtosis(X):
    return scipy_kurtosis(X)

def kl_divergence(X, bandwidth=0.1):
    """ Measures the difference between a the kernel density of X and a gaussian fitted to X
    """
    
    # Scaling data allows for consistent bandwidth
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(X)
    
    mean = np.mean(X, axis=0)
    cov = np.cov(X, rowvar=False)
    mvn = multivariate_normal(mean, cov)
    
    log_p = kde.score_samples(X) # estimated density at points X
    log_q = mvn.logpdf(X) # log probability density at points X on fitted multivariate normal distribution
    
    kl_div = np.mean(log_p - log_q)
    
    return kl_div
