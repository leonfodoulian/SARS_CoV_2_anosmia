# SARS_CoV_2_anosmia
This repository contains all codes used to analyse the data reported in Fodoulian et al. (2020) iScience. Helper R functions called with *source()* in the scripts will be uploaded in a separate repository as these are generic functions.

The scripts, their description and corresponding figure are:
* Scripts used to analyse the bulk RNA-seq data of human olfactory and respiratory biopsies and produce the plots from Figure 1.
	* human_moe_and_re.star_mapping.sh

	STAR mapping shell script.

	* human_moe_and_re.featurecounts.sh

	featureCounts quantification shell script.

	* human_moe_and_re.featureCounts_output_to_TPM_matrix.R

	R script used to reformat featureCounts output files for DESeq2 analysis and plotting.

	* human_moe_and_re.plot_TPMs.R

	R script used to produce the plots from Figure 1 D-G.

	* markers_list.txt

	List of genes used to produce the plots from Figure 1 D-G.

	* Human_MOE_RE_RNAseq_DESeq2_v1.22.2_apeglm_v1.4.2.R

	R script used to run the DESeq2 analysis and produce the plot from Figure 1 C.

* Script used to analyse the 10X Genomics single-cell RNA-seq data of human olfactory epithelial cells reported in Durante et al. (2020) Nat Neurosci and produce the plots from Figure 2.

	* DuranteNatNeurosci2020_seurat_clustering_as_reported_in_paper_seurat_v3.1.4.R

* Script used to analyse the SMART-Seq v4 single-nucleus RNA-seq data of human brain cells from the Allen Brain Map cell types database and produce the plots from Figure 6.
	
	* Human_Multiple_Cortical_Areas_SMART-seq.R

* Script used to analyse the Droplet-based single-nucleus RNA-seq data of human brain cells reported in Lake et al. (2018) Nat Biotechnol and produce the plots from Supplementary Figure 6.

	* LakeNatBiotechnol2018_pagoda2_v0.1.1.R

Note that the scripts are not "cleaned" and contain many lines that are commented (i.e. which were not used for the analysis). Please don't hesitate to contact me in case some sections of the scripts are not clear. 
