# UC_Discovery_Project

### Project Overview
#### Investigation Goals
This project aims to investigate the genetic drivers of ulcerative colitis, an inflammatory bowel disease, through the use of both traditional differential expression analysis (DEA) and unsupervised ML techniques (classifier models). Comparison will be made between the disease relevant genes and cell types derived from both methodologies.

### Dataset 
123006 cells from epithelial tissue,, 21784 genes, 3 disease states (healthy, inflammatory UC tissue, non-inflammatory UC tissue)
across 18 diseased patients and 12 healthy patients

Source data from the Broad institute: https://singlecell.broadinstitute.org/single_cell/study/SCP259/intra-and-inter-cellular-rewiring-of-the-human-colon-during-ulcerative-colitis?label=Best4%2B%20Enterocytes&cluster=Epithelial%20tSNE&spatialGroups=--&annotation=Cluster--group--study&subsample=100000

#### Analysis Pipeline
1. Cell quality control and count normalization/transformation
2. Cell type re-annotation into broad categories 
3. PCA reduction and tSNE embedding
4. Leiden clustering with hyperparameter sweep
5. Sample pseudobulking (cell type level and cluster assignment level)
6. DEA (cluster level pseudobulked samples)
7. Classifier training/testing (cell type level psuedobulked samples, 1 model per cell type)
8. Comparison of results between methods and literature review of top genes

#### Repository structure notes
All AnnData objects can be downloaded from git release
