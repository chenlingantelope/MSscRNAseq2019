# MSscRNAseq2019
Analysis for 2019 submission "Integrated single cell analysis of blood and cerebrospinal fluid leukocytes in multiple sclerosis" Schafflick, Xu, Hartlehnert et. al

* DE 
  * EdgeR: R scripts and bash script for running EdgeR to obtain DE genes
  * Wilcoxon: Python script for running EdgeR to obtain DE genes
* Vision
  * Create Vision Session for CD4(Vision_CD4.R) and all cells (Vision.R)
* prop_test
  * BetaBinom.R: R script for testing proportional differences between cell type abundances between tissue and disease states
* signatures
  * unmodified: original .gmt files 
  * allsigs.txt: all signatures sets used for the CSEA analysis
  * randsigs.txt: 1-to-1 matched signature set to allsigs.txt 
  * TFH.matched.txt: Mean matched random signature sets to TFH gene sets
  * Th1.matched.txt: Mean matched random signature sets to Th1 gene sets
* Notebooks
  * All Ipython Notebooks used to generate relevant analysis and figures, supplementary figures in the paper. Note that these will not be able to cloned and run locally because not all raw data files are available online pre-publication. They are meant to show the analysis pipeline and details. Also some notebooks are named by the Figure and Table name when their role is simply generating figure and table output, while others are named by the name of the analysis (DE, CSEA, PropTest etc.) 
  * meta/: clustering results, cell disease state and batch id for our dataset
  * utils.py: Functions for generating some of the custom plots such as dotplots, gene expression heatmaps etc.

## External Data Hosting
Integrated dataset with cell level metadata and UMAP 
https://figshare.com/articles/dataset/MS_CSF_h5ad/14356661
Raw single-cell RNA-seq data from the present study have been deposited in GEO repository with the accession code GSE138266.
