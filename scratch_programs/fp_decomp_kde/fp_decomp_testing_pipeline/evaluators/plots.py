""" Plotting functions
"""
import numpy as np
import matplotlib.pyplot as plt

def plot_decomp_kde_2D(X, Z, ax):
    """ Plot decomposition with KDE contours
    """
    xx, yy = np.meshgrid(
        np.linspace(-1.01, 1.01, 100), 
        np.linspace(-1.01, 1.01, 100)
    )

    levels = np.linspace(Z.min(), Z.max(), 12)
    ax.contourf(xx, yy, Z, levels=levels, cmap=plt.cm.Reds)
    ax.scatter(X[:,0], X[:,1], s=5)

    ax.grid(True)
    
    return


def plot_decomp_kde_2D_clustered(X, Z, ax, cluster_labels, cluster_lens):
    """ Plot decomposition with KDE contours
        color X by cluster label
        plot a mini bar graph on top of the axes for cluster lens
    """
    xx, yy = np.meshgrid(
        np.linspace(-1.01, 1.01, 100), 
        np.linspace(-1.01, 1.01, 100)
    )

    levels = np.linspace(Z.min(), Z.max(), 12)
    ax.contourf(xx, yy, Z, levels=levels, cmap=plt.cm.Reds)
    ax.scatter(X[:,0], X[:,1], s=5, c=cluster_labels, cmap="viridis")
    ax.grid(True)
    
    inset_ax = ax.inset_axes([0.1, 1.05, 0.8, 0.2])
    inset_ax.bar(np.unique(cluster_labels), cluster_lens, color="blue", alpha=0.7)    
    inset_ax.set_ylabel("Cluster Size")
    
    return

def plot_decomp_kde_2D_val(X, kde, ax, X_val):
    """ Calculate KDE contours on training set
        plot transparent X with X_val on top
    """
    xx, yy = np.meshgrid(
        np.linspace(-1.01, 1.01, 100), 
        np.linspace(-1.01, 1.01, 100)
    )

    Z = np.vstack([yy.ravel(), xx.ravel()]).T
    # density at gridpoints
    Z = np.exp(kde.score_samples(Z)).reshape(xx.shape)

    levels = np.linspace(Z.min(), Z.max(), 12)
    ax.contourf(xx, yy, Z, levels=levels, cmap=plt.cm.Reds)
    ax.scatter(X[:,0], X[:,1], s=5, alpha=0.5, color="gray")
    ax.scatter(X_val[:,0], X_val[:,1], s=5, alpha=0.75, color="black")
    
    ax.grid(True)
    
    return
