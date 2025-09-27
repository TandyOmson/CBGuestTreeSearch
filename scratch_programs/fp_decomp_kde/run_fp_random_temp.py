""" Script for testing fingerprint and decomposition candidates
    NOTE: CDK python wrapper (its in java) can't pipe the output
    in a shell without breaking the program. So don't do that.
"""
import os
import pickle

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem import MACCSkeys
from rdkit.Chem import AllChem
from rdkit import Avalon
from rdkit.Avalon import pyAvalonTools

from CDK_pywrapper import CDK, FPType

from sklearn.decomposition import PCA, KernelPCA
from sklearn.manifold import LocallyLinearEmbedding, Isomap, TSNE
import umap
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler
from sklearn.neighbors import KernelDensity
from sklearn.pipeline import Pipeline
from scipy.stats import skew, kurtosis, multivariate_normal
from sklearn.base import BaseEstimator, TransformerMixin

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec

# --------------
# Util Functions
# --------------
def mardia_tests(X):
    """ Returns two numpy arrays with skew along each axis and kurtosis along each axis
    """
    return kurtosis(X), skew(X)

def kl_divergence(X):
    """ Measures the difference between a the kernel density of X and a gaussian fitted to X
    """
    
    # Scaling data allows for consistent bandwidth
    kde = KernelDensity(kernel="gaussian", bandwidth=0.1).fit(X)
    
    mean = np.mean(X, axis=0)
    cov = np.cov(X, rowvar=False)
    mvn = multivariate_normal(mean, cov)
    
    log_p = kde.score_samples(X) # estimated density at points X
    log_q = mvn.logpdf(X) # log probability density at points X on fitted multivariate normal distribution
    
    kl_div = np.mean(log_p - log_q)
    
    return kl_div

# ---------
# Load data
# ---------

# PCCP data
df_pccp = pd.read_csv("/home/andyt/DProjects/DMCTS/VINA_ChemTSv2/data/dataset_benchmarks/HCs/PCCP_all_benchmark/out.csv", index_col=0, names=["smiles", "reward","VinaScore"])
df_pccp = df_pccp.sort_index()
df_pccp["mol"] = Chem.SDMolSupplier("/home/andyt/DProjects/DMCTS/VINA_ChemTSv2/scratch_programs/rdkit_embed/PCCP_top_100_first/all.sdf", removeHs=False)

mols_top = df_pccp.sort_values(by=["VinaScore"]).iloc[:100]

# Data generated randomly by the RNN
random_data_smiles = "/home/andyt/DProjects/Dcurated_datasets/Drandom_small_hydrocarbons_for_cb7/validity_tests/cleaned_data.smi"
print("Comapring training set with SMILES from {}".format(random_data_smiles))
df_mcts = pd.read_csv(random_data_smiles, names=["smiles"])
df_mcts["reward"] = len(df_mcts)*[0.01]
df_mcts["VinaScore"] = len(df_mcts)*[-10]

# Remove duplicate SMILES strings
invalid_mol_count = 0
unique_mol_count = 0
duplicate_mol_count = 0

df_mols = []
canon_smiles_map = []
for count, (smi, reward) in enumerate(zip(df_mcts["smiles"], df_mcts["reward"])):
    if reward == -1.0:
        invalid_mol_count += 1
        df_mcts = df_mcts.drop(count)
        continue
    
    mol = Chem.MolFromSmiles(smi)
    canon_smi = Chem.MolToSmiles(mol, isomericSmiles=True, kekuleSmiles=True, canonical=True)
    if canon_smi in canon_smiles_map:
        duplicate_mol_count += 1
        df_mcts = df_mcts.drop(count)
        
    else:
        unique_mol_count += 1
        canon_smiles_map.append(canon_smi)
        df_mols.append(mol)

df_mcts["mol"] = df_mols

class FingerprintDecomposition():

    def __init__(self, pipe, decomp_name, fp_func, fp_name):
        self.pipe = pipe
        self.kde = None
        self.decomp_name = decomp_name
        self.fp_func = fp_func
        self.fp_name = fp_name        
    
    def fit_decompose(self, df):
        """ Full function for checking out a single fingerprint with single (fitted) decomposition pipe
            - Fingerprint calculation and decomposition
            - Calculate gaussianity stats
            - Calculate kde
    
            # NOTE: scaling MUST be done before decomposition for vectors that have varying scales
                    obvs for BitVectors this is not required
        """

        try:
            fps = [self.fp_func(i) for i in df["mol"]]
        except:
            fps = self.fp_func(list(df["mol"])).values.tolist()

        X_decomp = self.pipe.fit_transform(fps)
    
        # Gaussianity tests
        # -----------------
        gauss_res = {}
        kurt, sk = mardia_tests(X_decomp)
        try:
            kl_div = kl_divergence(X_decomp)
        except:
            kl_div = 1.0
        gauss_res[self.fp_name + self.decomp_name] = [kurt, sk, kl_div]
        
        gauss_df = pd.DataFrame(gauss_res, index=["kurtosis", "skewness", "KL_divergence"]).T
    
        # KDE
        # ---
        kde = KernelDensity(kernel="gaussian", bandwidth=0.5).fit(X_decomp)
        
        return X_decomp, gauss_df, kde

    def decompose(self, df):
        """ Decomposes without fitting or all the other stuff
        """
        try:
            fps = [self.fp_func(i) for i in df["mol"]]
        except:
            fps = self.fp_func(list(df["mol"])).values.tolist()

        X_decomp = self.pipe.transform(fps)
        
        return X_decomp
    
    def plot_decomp(self, df, X_decomp, kde, ax):
        """
        df has VinaScore
        X_decomp is the decomposed mol fingerprints
        kde is a fitted kde object
        Plots:
            - Plot decomposition
            - plots of: affinity map, top 100
        """
        
        # 2D decomp plot
        #---------------
        sns.scatterplot(x=X_decomp[:, 0], y=X_decomp[:, 1], s=20, ax=ax)
        ax.set_xlabel(self.decomp_name + "1")
        ax.set_ylabel(self.decomp_name + "2")
        ax.grid(True)
        
        # mini histograms above 2D plot
        inset_ax = ax.inset_axes([0.1, 1.05, 0.8, 0.2])  # position and size
        inset_ax.hist(X_decomp[:, 0], bins=20, color='blue', alpha=0.7)
        inset_ax.set_ylabel("1D")
        inset_ax.set_xticklabels([])
        inset_ax.set_yticklabels([])
        inset_ax.grid(True)

        return

    def plot_kde(self, df, X_decomp, kde, ax):
        # KDE plot
        # --------
        # evaluate density across latent space
        X, Y = np.meshgrid(np.linspace(-1.1, 1.1, 100), np.linspace(-1.1, 1.1, 100))
        Z = np.vstack([Y.ravel(), X.ravel()]).T
        # density at gridpoints
        Z = np.exp(kde.score_samples(Z)).reshape(100,100)
    
        # kernel density countour plot
        levels = np.linspace(0, Z.max(), 25)
        m = ax.contourf(X, Y, Z, levels=levels, cmap=plt.cm.Reds)

#        ax.set_xlabel(self.decomp_name + "1")
#        ax.set_ylabel(self.decomp_name + "2")
        ax.grid(True)
#        plt.colorbar(m, label="Gaussian Kernel Density", ax=ax)

        # affinity shading
        sns.scatterplot(x=X_decomp[:, 0], y=X_decomp[:, 1], hue=df["VinaScore"], palette='viridis', s=20, ax=ax)
        norm = plt.Normalize(df['VinaScore'].min(), df['VinaScore'].max())
        sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
#        plt.colorbar(sm, label='Affinity', ax=ax)

        ax.legend([])

        return m

    def plot_top_100(self, df, X_decomp, kde, ax):
        # Top 100 on KDE plot
        # -------------------
        X, Y = np.meshgrid(np.linspace(-1.1, 1.1, 100), np.linspace(-1.1, 1.1, 100))
        Z = np.vstack([Y.ravel(), X.ravel()]).T
        # density at gridpoints
        Z = np.exp(kde.score_samples(Z)).reshape(100,100)

        # kernel density countour plot
        levels = np.linspace(0, Z.max(), 25)
        m = ax.contourf(X, Y, Z, levels=levels, cmap=plt.cm.Reds)

#        ax.set_xlabel(self.decomp_name + "1")
#        ax.set_ylabel(self.decomp_name + "2")
        ax.grid(True)
#        plt.colorbar(m, label="Gaussian Kernel Density", ax=ax)

        # get top 100 from df
        topmols = df.sort_values(by=["VinaScore"]).reset_index().iloc[:100]

        # top 100 affinity shading
        sns.scatterplot(x=X_decomp[:, 0], y=X_decomp[:, 1], s=10, ax=ax, alpha=0.4, color="gray")
        top_X = [X_decomp[i, 0] for i in topmols.index]
        top_y = [X_decomp[i, 1] for i in topmols.index]
        sns.scatterplot(x=top_X, y=top_y, hue=topmols["VinaScore"], palette='viridis', s=20, ax=ax, marker="s")
        
        norm = plt.Normalize(topmols['VinaScore'].min(), topmols['VinaScore'].max())
        sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
        plt.colorbar(sm, label='Affinity', ax=ax)
        ax.legend([])

        return

    def plot_master_figure(self, df_train, X_decomp_train, kde, df_test, X_decomp_test, df_gauss, outdir=""):
        fig = plt.figure(figsize=(12, 8), constrained_layout=True)
        fig.suptitle(self.fp_name + "_" + self.decomp_name, size="xx-large")
        gs = GridSpec(4, 4, figure=fig,
                      height_ratios=[1,1,1,1],  # bottom row smaller for table
                      width_ratios=[1,1,2,2]) # right side wider
        
        # Create axes
        ax_big   = fig.add_subplot(gs[0:2, 0:2])
        ax_table = fig.add_subplot(gs[2, 0:2])
        ax_cmaps = fig.add_subplot(gs[3, 0:2])
        ax_2     = fig.add_subplot(gs[0:2, 2])
        ax_3     = fig.add_subplot(gs[0:2, 3])
        ax_4     = fig.add_subplot(gs[2:4, 2])
        ax_5     = fig.add_subplot(gs[2:4, 3])

        self.plot_decomp(df_train, X_decomp_train, kde, ax=ax_big)

        # returns scalar mappable for colorbar
        sm = self.plot_kde(df_train, X_decomp_train, kde, ax=ax_2)
        ax_2.set_title("PCCP")
        self.plot_top_100(df_train, X_decomp_train, kde, ax=ax_3)
        ax_3.set_title("PCCP Top 100")

        # Plot the coloarmap for kde on a separate axis
        cbar = fig.colorbar(sm, cax=ax_cmaps, orientation="horizontal")
        cbar.set_label("Kernel Density Estimation")
        ax_cmaps.xaxis.set_label_position("top")
        ax_cmaps.tick_params(axis='x', rotation=90)

        self.plot_kde(df_test, X_decomp_test, kde, ax=ax_4)
        ax_4.set_title("MCTS")
        self.plot_top_100(df_test, X_decomp_test, kde, ax=ax_5)
        ax_5.set_title("MCTS Top 100")
        
        # Place gaussianity metrics underneath left hand side
        ax_table.axis("off")
#        tbl = pd.plotting.table(ax_table, df_gauss, loc="center")
#        tbl.auto_set_font_size(False)
#        tbl.set_fontsize(8)
        
        for i, row in enumerate(df_gauss.columns):
            if i != 2:
                row_text = "{}: {}".format(row, [f"{j:.3f}" for j in list(df_gauss[row].values[0])])
            else:
                row_text = "{}: {:.3f}".format(row, df_gauss[row].values[0])
            ax_table.text(0, 0.7 - i*0.2, row_text, fontsize=12, va='top', ha='left')
        ax_table.text(0, 0.9, "Normality Metrics", fontsize=12, va='top', ha='left')
                
#        plt.tight_layout()
        fig.savefig(outdir + "/" + self.fp_name + "_" + self.decomp_name + ".png", dpi=300)
        plt.close(fig)
        return

    def save_pipeline_and_kde(self, pipeline, kde, outdir):
        pipeout = outdir + "/" + self.fp_name + "_" + self.decomp_name + "_" + "pipe.pkl"
        kdeout = outdir + "/" + self.fp_name + "_" + self.decomp_name + "_" + "kde.pkl"

        with open(pipeout, "wb") as fw:
            pickle.dump(pipeline, fw)

        with open(kdeout, "wb") as fw:
            pickle.dump(kde, fw)
        
        return

    def save_pipeline_and_kde_w_fp_func(self, pipeline, kde, fp_func, outdir):
        """ Currently doesnt work
        """
        fp_func = make_pickleable_fp_func(fp_func)
        fp_transformer = FingerprintTransformer(fp_func)
        
        pipeline = Pipeline([(self.fp_name, fp_transformer)] + pipeline.steps)
        
        pipeout = outdir + "/" + self.fp_name + "_" + self.decomp_name + "_" + "pipe.pkl"
        kdeout = outdir + "/" + self.fp_name + "_" + self.decomp_name + "_" + "kde.pkl"

        with open(pipeout, "wb") as fw:
            pickle.dump(pipeline, fw)

        with open(kdeout, "wb") as fw:
            pickle.dump(kde, fw)
        
        return

    def get_max_density_point(self, X_decomp, kde):
        """ Get the point of max density and density value from a set of X decomposed points and a trained kde object
        """
        kde_at_decomp = np.exp(kde.score_samples(X_decomp))
        max_idx = np.argmax(kde_at_decomp)

        max_point = X_decomp[max_idx]
        max_density = kde_at_decomp.max()

        return max_point, max_density

# Fingerprint decomposition consists of a decomposition pipeline and a fingerprinting function that takes a mol object

# Dictonary of decomposition alogrithms
alg_dict = { 'PCA' : Pipeline([('decomp', PCA(n_components=2)), ('scaler', MinMaxScaler((-1,1)))]),
             'PCA_RBF' : Pipeline([('decomp', KernelPCA(n_components=2, kernel="rbf", gamma=0.2)), ('scaler', MinMaxScaler((-1,1)))]),
             'isomap' : Pipeline([('decomp', Isomap(n_components=2, n_neighbors=200)), ('scaler', MinMaxScaler((-1,1)))]),
             'umap' : Pipeline([('decomp', umap.UMAP(n_components=2)), ('scaler', MinMaxScaler((-1,1)))]),
             'lle' : Pipeline([('decomp', LocallyLinearEmbedding(n_components=2, n_neighbors=500)), ('scaler', MinMaxScaler((-1,1)))]),
             'tnse' : Pipeline([('decomp', TSNE(n_components=2)), ('scaler', MinMaxScaler((-1,1)))]),
           }

# Need to initialiase CDK Fingerprints first so the library won't initialise every time
cdk_pubchem = CDK(fingerprint=FPType.PubchemFP)
cdk_cdk = CDK(fingerprint=FPType.FP)
cdk_estate = CDK(fingerprint=FPType.EStateFP)
cdk_kr = CDK(fingerprint=FPType.KRFP)
cdk_lingo = CDK(fingerprint=FPType.LingoFP)

# Dictionary of molecular fingerprinting functions
fp_dict = { 'morgan_ECFP4' : rdFingerprintGenerator.GetMorganGenerator(radius=2).GetFingerprint,
            'morgan_FCFP4' : rdFingerprintGenerator.GetMorganGenerator(radius=2, atomInvariantsGenerator=AllChem.GetMorganFeatureAtomInvGen()).GetFingerprint,
            'MACCS' : MACCSkeys.GenMACCSKeys,
            'rdkit' : AllChem.GetRDKitFPGenerator().GetFingerprint,
            'avalon' : pyAvalonTools.GetAvalonFP,
            'topo_torsion' : rdFingerprintGenerator.GetTopologicalTorsionGenerator().GetFingerprint,
            'atom_pair' : rdFingerprintGenerator.GetAtomPairGenerator().GetFingerprint,
            'pubchem' : cdk_pubchem.calculate,
            'cdk' : cdk_cdk.calculate,
            'EState' : cdk_estate.calculate,
            'klekota_roth' : cdk_kr.calculate,
            'lingo' : cdk_lingo.calculate,
           }

decompers = {}

# All combinations of decompers
#for alg_name, alg in alg_dict.items():
#    for fp_name, fp_func in fp_dict.items():
#        d = FingerprintDecomposition(alg, alg_name, fp_func, fp_name)
#        decompers[fp_name + "_" + alg_name] = d
# outdir = decomper.decomp_name

# Specfic decompers to use
#decompers["MACCS" + "_" + "PCA"] = FingerprintDecomposition(alg_dict["PCA"], "PCA", fp_dict["MACCS"], "MACCS")
decompers["morgan_ECFP4" + "_" + "PCA"] = FingerprintDecomposition(alg_dict["PCA"], "PCA", fp_dict["morgan_ECFP4"], "morgan_ECFP4")
#decompers["atom_pair" + "_" + "PCA"] = FingerprintDecomposition(alg_dict["PCA"], "PCA", fp_dict["atom_pair"], "atom_pair")
#decompers["MACCS" + "_" + "PCA_RBF"] = FingerprintDecomposition(alg_dict["PCA_RBF"], "PCA_RBF", fp_dict["MACCS"], "MACCS")
#decompers["klekota_roth" + "_" + "PCA_RBF"] = FingerprintDecomposition(alg_dict["PCA_RBF"], "PCA_RBF", fp_dict["klekota_roth"], "klekota_roth")
#decompers["MACCS" + "_" + "isomap"] = FingerprintDecomposition(alg_dict["isomap"], "isomap", fp_dict["MACCS"], "MACCS")
#decompers["avalon" + "_" + "isomap"] = FingerprintDecomposition(alg_dict["isomap"], "isomap", fp_dict["avalon"], "avalon")
#decompers["atom_pair" + "_" + "isomap"] = FingerprintDecomposition(alg_dict["isomap"], "isomap", fp_dict["atom_pair"], "atom_pair")
outdir = "selected_candidates_RNN_random_generated"

for count, (key, decomper) in enumerate(decompers.items()):
    if os.path.exists(outdir + "/" + decomper.fp_name + "_" + decomper.decomp_name + ".png" ):
        print(key, "exists, not running...")
        continue
    else:
        try:
            os.makedirs(outdir)
        except:
            pass
    
    print("running {}, {} of {}".format(key, count, len(decompers)))
    # fit on PCCP
    X_decomp_pccp, gauss_df, kde = decomper.fit_decompose(df_pccp)
    # decompose MCTS
    X_decomp_mcts = decomper.decompose(df_mcts)

    decomper.plot_master_figure(df_pccp, X_decomp_pccp, kde, df_mcts, X_decomp_mcts, gauss_df, outdir=outdir)
    # Save fitted pipeline and kde objects
    decomper.save_pipeline_and_kde(decomper.pipe, kde, outdir)
    
    max_point, max_density = decomper.get_max_density_point(X_decomp_pccp, kde)
        
    print(key, gauss_df)
    print("max density at {}, value: {}".format(max_point, max_density))
