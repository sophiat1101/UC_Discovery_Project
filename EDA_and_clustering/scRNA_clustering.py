import anndata as ad
ad.settings.allow_write_nullable_strings = True
import scanpy as sc
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import sys
from sklearn.metrics import adjusted_rand_score
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "preprocess"))
from cell_type_grouping import groupCellTypes

script_dir = Path(__file__).resolve().parent
project_dir = script_dir.parent
processed_dir = project_dir / "processed_data"
figures_dir = script_dir / "figures"
figures_dir.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(processed_dir / "uc_Epi_preprocessed_PCA.h5ad")

# grouping cell type annotations and adding as new obs to adata as cell_type_group
adata = groupCellTypes(adata)

# EDA -- TSNE

# colored by grouped cell type
sc.tl.tsne(adata, n_pcs = 30, use_rep = "X_pca_harmony")
sc.pl.tsne(adata, color = "cell_type_group", title = "tSNE Colored by Broad Cell Type", show = False)
plt.savefig(figures_dir / "tsne_grouped_cell_type.png", dpi = 300)
plt.tight_layout()
plt.show()

# colored by granular cell type
fig, ax = plt.subplots(figsize = (10, 7))
sc.pl.tsne(adata, color = "Cluster", title = "t-SNE Colored by Specific Cell Type", ax = ax, show = False)
# locking plot size
ax.set_box_aspect(1)
plt.savefig(figures_dir / "tsne_specific_cell_type.png", dpi = 300, bbox_inches = "tight")
plt.show()

# colored by disease state
sc.pl.tsne(adata, color = "Health", title = "t-SNE Colored by Disease State", show = False)
plt.savefig(figures_dir / "tsne_disease_state.png", dpi = 300)
plt.tight_layout()
plt.show()


# LEIDEN CLUSTERING -- SWEEP OVER N_NEIGHBORS, RESOLUTION AND PCS

n_neighbors = [10, 15, 20, 30]
pcValues = [10, 15, 20, 30, 40, 50]
resValues = [0.25, 0.5, 0.75, 1, 1.5, 2]
results = []

for k in n_neighbors:
    for npcs in pcValues:
        nnKey = f"nn{k}_pc{npcs}"
        sc.pp.neighbors(adata, n_neighbors = k, n_pcs = npcs, use_rep = "X_pca_harmony", key_added= nnKey)

        for r in resValues:
            leidenKey = f"leiden_nn{k}_pc{npcs}_res{r}"
            sc.tl.leiden(adata, resolution = r, key_added= leidenKey, neighbors_key = nnKey, flavor = "igraph",
                         n_iterations = 2, directed = False)

            results.append({"n_neighbors": k,
                            "n_pcs": npcs,
                            "resolution": r,
                            "n_clusters": adata.obs[leidenKey].nunique(),
                            "ARI_specificCells": adjusted_rand_score(adata.obs["Cluster"], adata.obs[leidenKey]),
                            "ARI_groupedCells": adjusted_rand_score(adata.obs["cell_type_group"], adata.obs[leidenKey])})

        # deleting sweep graph from adata
        del adata.obsp[f"{nnKey}_distances"], adata.obsp[f"{nnKey}_connectivities"], adata.uns[nnKey]

sweepResults = pd.DataFrame(results)
print(sweepResults.sort_values("ARI_specificCells", ascending = False).head(10))

metrics = ["ARI_specificCells", "ARI_groupedCells", "n_clusters"]
titles = ["ARI vs Specific Cell Type", "ARI vs Broad Cell Type", "Number of Leiden Clusters"]

fig, axes = plt.subplots(len(n_neighbors), 3, figsize = (18, 5*len(n_neighbors)))
for row, k in zip(axes, n_neighbors):
    subset = sweepResults[sweepResults["n_neighbors"] == k]
    for ax, col, title in zip(row, metrics, titles):
        table = subset.pivot(index = "n_pcs", columns = "resolution", values = col)
        sns.heatmap(table, annot = True, fmt = ".2f" if col != "n_clusters" else "d", cmap = "viridis", ax = ax)
        ax.set_title(f"{title} ({k} neighbors)")
        ax.set_xlabel("Resolution")
        ax.set_ylabel("Number of PCs")
plt.savefig(figures_dir / "tsne_sweep_heatmap.png", dpi = 300, bbox_inches = "tight")
plt.show()


# TSNE COLORED BY LEIDEN CLUSTER ASSIGNMENT AFTER SWEEP -- maximize ARI score of leiden clusters with known
# labels for both specific cell types and broad cell types

neighbors = 20
pcs = 30
res = 0.25

# parameters chosen from sweep: 20 neighbors, 30 pcs, 0.25 resolution
sc.pp.neighbors(adata, n_neighbors = neighbors, n_pcs = pcs, use_rep = "X_pca_harmony")

# n_iterations = -1 so it keeps clustered until no improvement
# directed = False so connections between cells is two-way
sc.tl.leiden(adata, resolution = res, key_added = "leiden", flavor = "igraph", n_iterations = -1, directed = False)

# plotting t-SNE colored by leiden cluster assignment
sc.pl.tsne(adata, color = "leiden", title = f"t-SNE Colored by Leiden Cluster\nn_neighbors={neighbors}, resolution={res}, n_pcs={pcs}", show = False)
plt.savefig(figures_dir / "tsne_leiden.png", dpi = 300, bbox_inches = "tight")
plt.show()

# SAVING FINAL LEIDEN CLUSTER ASSIGNMENTS

# dropping cluster cols
adata.obs = adata.obs.loc[:, ~adata.obs.columns.str.startswith("leiden_nn")]

# recording leiden parameters
adata.uns["leiden_params"] = {"n_neighbors": neighbors, "n_pcs": pcs, "resolution": res, "use_rep": "X_pca_harmony"}

# writing object with clusters, neighbor group and tSNE -- contains leiden column, tSNE coordinates, neighbor graph,
# and chosen parameters for leiden
adata.write_h5ad(processed_dir / "uc_Epi_clustered.h5ad", compression = "gzip")

# making csv with cell barcode mapped to leiden cluster table with original and grouped cell type, patient id and disease state
(adata.obs[["leiden", "Cluster", "cell_type_group", "Subject", "Health"]].rename_axis("barcode").to_csv(processed_dir / "uc_Epi_clustered.csv", index = False))

