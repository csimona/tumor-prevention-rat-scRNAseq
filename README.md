Breast cancer prevention by short-term inhibition of TGFB signaling
--------

This repository contains code and plots to support genomic analyses and reproduce results from the paper [Breast cancer prevention by short-term inhibition of TGFB signaling](https://).


Abstract
--------
Cancer prevention has a profound impact on cancer-associated mortality and morbidity. We previously identified TGFB signaling as a candidate regulator of mammary epithelial cells associated with breast cancer risk. Here, we show that short-term TGFBR inhibitor (TGFBRi) treatment of peripubertal ACI inbred and Sprague Dawley outbred rats induces lasting changes and prevents estrogen- and carcinogen-induced mammary tumors, respectively. We identify TGFBRi-responsive cell populations by single cell RNA-sequencing, including a unique epithelial subpopulation designated secretory basal cells (SBCs) with progenitor features. We detect SBCs in normal human breast tissues and find them to be associated with breast cancer risk. Interactome analysis identifies SBCs as the most interactive cell population and the main source of insulin-IGF signaling. Accordingly, inhibition of TGFBR and IGF1R decrease proliferation of organoid cultures. Our results reveal a critical role for TGFB in regulating mammary epithelial cells relevant to breast cancer and serve as a proof-of-principle cancer prevention strategy.


<p align="center"> 
<img src="overview_github.png" width="900" height="400"/>
</p>


Content
-------
* `/01_cellranger/`: cellranger quality control plots on the three datasets: sorted ACI epithelial cells, sorted SD epithelial cells, whole mammary gland SD cells. 
* `/02_R_single_cell/`: R scripts to reproduce the single cell analyses.
* `/03_R_bulk/`: R scripts to reproduce the bulk RNASeq analyses.
* `/data/`: preprocessed bulk outputs, cellranger count outputs and htseq count outputs.
* `/plots/`: output plots from the single cell and bulk analyses, organized by folder.
* `/tables/`: number of cells per cluster, relevant gene lists used throughout the analysis, differentially expressed genes from bulk analyses, and sample metadata.


Data
-------
The raw scRNAseq and bulk RNASeq data is accessible at [GEO GSEA]().

Preprocessed RDS objects can be downloaded from the following [zenodo link](https://zenodo.org/record/7293642#.Y4krEi-B19f):

HTSeq counts inputs from the bulk RNASeq analysis can be downloaded from `data\`.


## Contact
For any questions, comments or suggestions, please feel free out to me via GitHub, Twitter or email (scristea@jimmy.harvard.edu).