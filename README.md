## Code for Ahmad et. al. 2023
This repository contains the code to reproduce all analyses used for *Ahmad et. al. 2023*. 

The analysis started with the raw output of [CellRanger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) run on the FASTQs available at (GEO link).

The main script is `all_clustering/all_clustering.R`.  The scripts `macrophage_clustering/macrophage_clustering.R` and `t_cell_clustering/t_cell_clustering.R` were used for the macrophage/cDC and T-cell/NK cell clustering.

The `frozen_clustering` folder has an identical structure to the main folder and was used for the frozen-only analysis.

The `figure` folders and `supplemental_table` folder contain the code used to generate the figures and supplemental tables.

All R packages used in the analysis are available on CRAN or Bioconductor, and provided as packages on conda-forge or bioconda.
