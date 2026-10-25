import anndata as ad
ad.settings.allow_write_nullable_strings = True
import numpy as np
import scanpy as sc
import pandas as pd
from scipy.sparse import issparse
import warnings
warnings.filterwarnings("ignore", message=".*encountered in (vecdot|matmul)", category=RuntimeWarning)
from pathlib import Path
import matplotlib.pyplot as plt

# function load data
def loadData(mtxPath, genePath, barcodePath):

    # mtx is genes x cells --> transposing  so cells are rows
    adata = sc.read_mtx(mtxPath).T

    # tsv files have no header and only one column
    genes = pd.read_csv(genePath, header=None, sep="\t")
    barcodes = pd.read_csv(barcodePath, header=None, sep="\t")

    adata.obs_names = barcodes[0].values
    adata.var_names = genes[0].values
    adata.var_names_make_unique()

    # getting sample data for cells -- in format Ex. N7.EpiA.barcode --> dropping barcode to keep N7.EpiA
    adata.obs["samples"] = adata.obs_names.str.rsplit(".", n=1).str[0]

    return adata


# quality check
def qualityCheck(adata):

    # checking if data is already filtered/normalized
    X = adata.X.data if issparse(adata.X) else adata.X.ravel()
    print("Min:", X.min(), "Max:", X.max())
    print("All integers:", np.all(np.mod(X, 1) == 0))

    # checking number of samples and stats
    counts = adata.obs["samples"].value_counts()
    print(f"Number of samples: {len(counts)}")
    print(counts.describe())

    # mitochondrial rna
    adata.var["mt"] = adata.var_names.str.startswith("MT-")

    # getting qc metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    qc = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    print(adata.obs[qc].describe(percentiles = [0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]))

    # violin plot of qc metrics
    v = sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], stripplot=False,
                     multi_panel=True, show=False)
    titles = ["Genes per Cell", "Total Counts", "Mitochondrial Reads"]
    ylabels = ["Number of Genes", "Counts", "Percent (%)"]

    for ax, title, label in zip(v.axes.flat, titles, ylabels):
        ax.set_title(title)
        ax.set_ylabel(label)
        ax.set_xlabel("")

    v.fig.suptitle("scRNA-Seq QC Metrics", fontsize=14, y=1.05)
    v.fig.savefig(str(script_dir / "figures/qc_violin.png"), dpi=300, bbox_inches="tight")
    plt.close()

    # checking how low quality cells there are
    print(f"Total Cells: {adata.n_obs}")
    print(f"Cells with <200 genes: {(adata.obs.n_genes_by_counts < 200).sum()}")
    print(f"Cells with >5000 genes: {(adata.obs.n_genes_by_counts > 5000).sum()}")
    print(f"Cells with >20% mitochondrial RNA: {(adata.obs.pct_counts_mt > 20).sum()}")


# filtering cells based on cutoffs
def filterCells(adata, minGeneCutoff = 200, minCellCutoff = 3, mitoCutoff = 0.2):

    # filtering low expression cells and cells with too much mitochondrial genes
    sc.pp.filter_cells(adata, min_genes = minGeneCutoff)
    sc.pp.filter_genes(adata, min_cells = minCellCutoff)
    adata = adata[adata.obs.pct_counts_mt < mitoCutoff].copy()

    return adata


# normalizing and transforming gene counts
def normalize(adata, normalizingTargetNum, highlyVarGenesNum):

    # saving raw counts
    adata.layers["counts"] = adata.X.copy()

    # normalizing gene counts for cells
    sc.pp.normalize_total(adata, target_sum = normalizingTargetNum)

    # log 1 transforming
    sc.pp.log1p(adata)

    # saving "raw" adata
    adata.raw = adata

    # finding highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=highlyVarGenesNum)

    # filtering to only highly variable genes
    adata = adata[:, adata.var.highly_variable].copy()

    return adata

def runPCA(adata, nComps = 50, scale = True):

    # scaling data
    if scale:
        sc.pp.scale(adata, max_value = 10)

    sc.pp.pca(adata, n_comps = nComps)

    varExplained = adata.uns["pca"]["variance_ratio"]

    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False)
    plt.ylabel("Principal Components")
    plt.xlabel("Fraction of Total Variance Explained")
    plt.title("PCA Scree Plot")

    plt.savefig(str(script_dir / "figures/pca_expl_var.png"), dpi=300, bbox_inches='tight')
    plt.close()

    return adata

### RUNNING ON EPITHELIAL TISSUE

script_dir = Path(__file__).resolve().parent
project_dir = script_dir.parent

mtxPath = project_dir / "data/gene_sorted-Epi.matrix.mtx"
genePath = project_dir / "data/Epi.genes.tsv"
barcodePath = project_dir / "data/Epi.barcodes2.tsv"

print("Reading in 10X files...")

adata = loadData(mtxPath, genePath, barcodePath)

print("Assessing data quality...")

qualityCheck(adata)

adata = filterCells(adata)

adata = normalize(adata, 1000, 2000)

print("Data filtered and normalized, saving data object...")

adata.write_h5ad(str(project_dir / "processed_data/uc_Epi_preprocessed_noPCA.h5ad"), compression = "gzip")

print("Computing PCs, saving second data object...")

adata_pca = runPCA(adata)

adata_pca.write_h5ad(str(project_dir / "processed_data/uc_Epi_preprocessed_PCA.h5ad"), compression = "gzip")

print("Done!")