
# CroCoNet analyses

This repository contains the code to reproduce all analyses in the
manuscript

### **CroCoNet: a quantitative approach for comparing gene regulatory networks across species**

by Anita Térmeg, Johanna Geuder, Vladyslav Storozhuk, Zane Kliesmete,
Fiona C. Edenhofer, Beate Vieth, Philipp Janssen, Ines Hellmann

# <img src="pipeline.png" align="center" width="1000" />

The data necessary to reproduce this analysis can be found on
ArrayExpress and GEO:

| Accession                   | Dataset          |
|-----------------------------|------------------|
| E-MTAB-XXXXX                | scRNA-seq data   |
| E-MTAB-13373 & E-MTAB-15654 | ATAC-seq data    |
| GSE298717                   | Perturb-seq data |

To be able to smoothly run all analyses, please follow these steps:

1.  Start a new R project in a new directory

2.  Clone this Github repository to a subdirectory “scripts”

3.  Download the linked Zenodo repository to a subdirectory “data”

### 1. Neural differentiation dataset

The main example data throughout the paper is a scRNA-seq dataset
obtained during the early neural differentiation of human, gorilla and
cynomolgus macaque induced pluripotent (iPS) cell lines. The relevant
scripts to analyze this dataset are the following:  

- [Mapping & QC](1.neural_differentiation_dataset/1.1.mapping_and_QC/)  
  - [Download the reference
    genomes](1.neural_differentiation_dataset/1.1.mapping_and_QC/1.1.1.download_ref_genomes.sh)  
  - [Create non-human primate annotations using
    Liftoff](1.neural_differentiation_dataset/1.1.mapping_and_QC/1.1.2.run_liftoff.sh)  
  - [Filter Liftoff
    annotations](1.neural_differentiation_dataset/1.1.mapping_and_QC/1.1.3.filter_liftoff_gtfs.R)  
  - [Remove small contigs from the gorGor6 genome for efficient
    computation](1.neural_differentiation_dataset/1.1.mapping_and_QC/1.1.4.remove_small_contigs.R)  
  - [Create STAR
    indices](1.neural_differentiation_dataset/1.1.mapping_and_QC/1.1.5.create_STAR_indices.sh)  
  - [Download FASTQ data of the scRNA-seq
    experiment](1.neural_differentiation_dataset/1.1.mapping_and_QC/1.1.6.download_FASTQ.sh)  
  - [Trim poly-A
    tails](1.neural_differentiation_dataset/1.1.mapping_and_QC/1.1.7.trimming.sh)  
  - [Map reads using the zUMIs
    pipeline](1.neural_differentiation_dataset/1.1.mapping_and_QC/1.1.8.mapping.sh)  
  - [YAML for mapping to the hg38
    genome](1.neural_differentiation_dataset/1.1.mapping_and_QC/hg38.yaml)  
  - [YAML for mapping to the gorGor6
    genome](1.neural_differentiation_dataset/1.1.mapping_and_QC/gorGor6.yaml)  
  - [YAML for mapping to the macFas6
    genome](1.neural_differentiation_dataset/1.1.mapping_and_QC/macFas6.yaml)  
- [Network
  inference](1.neural_differentiation_dataset/1.2.network_inference/)  
- [CroCoNet
  analysis](1.neural_differentiation_dataset/1.3.CroCoNet_analysis/)  
- [Additional
  analyses](1.neural_differentiation_dataset/1.4.additional_analysis/)  

### 2. Validations

To validate the results of the CroCoNet analysis on the neural
differentiation dataset, we used various additional datasets and
databases. These lines of analysis can be reproduced using the following
scripts:  

- [Regulator interactions and module
  overlaps](2.validations/2.1.regulator_interactions_and_module_overlaps/)  
- [Pathway enrichment](.validations/2.2.pathway_enrichment/)  
- [Binding site enrichment and
  divergence](2.validations/2.3.binding_site_enrichment_and_divergence/)  
- [Sequence divergence](2.validations/2.4.sequence_divergence/)  
- [Expression pattern
  divergence](2.validations/2.5.expression_pattern_divergence/)  
- [Analysis of POU5F1 ChIP-seq
  data](2.validations/2.6.POU5F1_ChIP_seq/)  
- [Enrichment of LTR7 elements near POU5F1 module
  members](2.validations/2.7.POU5F1_LTR7_enrichment/)  
- [Aanalysis of POU5F1 Perturb-seq
  data](2.validations/2.8.POU5F1_Perturb_seq/)  
