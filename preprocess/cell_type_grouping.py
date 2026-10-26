# TO USE:
# add this to top of script when need
# import sys
# sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "preprocess") OR change for wherever location ur using
# from cell_type_grouping import groupCellTypes

import pandas as pd

# fine-grained cell types from metaData.txt "Cluster" column --> broader cell type groups
cellTypeGroups = {
    # epithelial
    "Stem": "Stem",
    "TA 1": "Transit Amplifying",
    "TA 2": "Transit Amplifying",
    "Cycling TA": "Transit Amplifying",
    "Secretory TA": "Transit Amplifying",
    "Enterocyte Progenitors": "Enterocytes",
    "Immature Enterocytes 1": "Enterocytes",
    "Immature Enterocytes 2": "Enterocytes",
    "Enterocytes": "Enterocytes",
    "Best4+ Enterocytes": "Enterocytes",
    "Immature Goblet": "Goblet",
    "Goblet": "Goblet",
    "Enteroendocrine": "Enteroendocrine",
    "Tuft": "Tuft",
    "M cells": "M cells",

    # immune
    "CD4+ Activated Fos-hi": "T cells",
    "CD4+ Activated Fos-lo": "T cells",
    "CD4+ Memory": "T cells",
    "CD4+ PD1+": "T cells",
    "Tregs": "T cells",
    "CD8+ IELs": "T cells",
    "CD8+ IL17+": "T cells",
    "CD8+ LP": "T cells",
    "Cycling T": "T cells",
    "NKs": "NK/ILCs",
    "ILCs": "NK/ILCs",
    "Follicular": "B cells",
    "GC": "B cells",
    "Cycling B": "B cells",
    "Plasma": "Plasma",
    "Macrophages": "Myeloid",
    "Inflammatory Monocytes": "Myeloid",
    "Cycling Monocytes": "Myeloid",
    "DC1": "Myeloid",
    "DC2": "Myeloid",
    "CD69+ Mast": "Mast",
    "CD69- Mast": "Mast",
    "MT-hi": "MT-hi",

    # stromal
    "WNT2B+ Fos-hi": "Fibroblasts",
    "WNT2B+ Fos-lo 1": "Fibroblasts",
    "WNT2B+ Fos-lo 2": "Fibroblasts",
    "WNT5B+ 1": "Fibroblasts",
    "WNT5B+ 2": "Fibroblasts",
    "RSPO3+": "Fibroblasts",
    "Inflammatory Fibroblasts": "Fibroblasts",
    "Myofibroblasts": "Fibroblasts",
    "Endothelial": "Endothelial",
    "Microvascular": "Endothelial",
    "Post-capillary Venules": "Endothelial",
    "Pericytes": "Pericytes",
    "Glia": "Glia",
}


# adding grouped cell type column to adata.obs --> cell types not in the dict keep their original label
def groupCellTypes(adata, cellTypeCol = "Cluster", groupCol = "cell_type_group", groups = cellTypeGroups):

    cellTypes = adata.obs[cellTypeCol].astype(str)
    adata.obs[groupCol] = pd.Categorical(cellTypes.map(groups).fillna(cellTypes))

    return adata
