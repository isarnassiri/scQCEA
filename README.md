scQCEA
==========
* [Introduction](#introduction)
* [Installation](#Installation)
* [Usage](#Usage)
* [Citation](#citation)
<a name="introduction"/>

### Introduction

scQCEA (acronym of the single-cell RNA sequencing Quality Control and Enrichment Analysis) is an R package for annotation and quality control report of scRNA-Seq profiles, which performs a probabilistic assignment of the reference cell types to identify clusters, before downstream analysis such as gene network inference. scQCEA provides automated cell type annotation on scRNA-seq data and identifies differential patterns in gene expression. scQCEA generates an interactive report of quality control metrics which allows visual evaluation of QC metrics, objective selection of insightful optimal cluster numbers and discrimination between true variation and background noise. 

Please see the [`manual`](https://isarnassiri.github.io/scQCEA/) for the usage of scQCEA including the explanations of the HTML report and how to prepare data input files.

<a name="installation"/>

### Installation
1. Install the R [(LINK)](https://cran.r-project.org/)
2. Install the free version of rStudio [(LINK)](https://www.rstudio.com/products/rstudio/download/)
3. Run the following command in rStudio to install scQCEA as an R package:

```{r,eval=FALSE}
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
    
library(devtools)
devtools::install_github("walkerke/bsselectR")
install_github("isarnassiri/scQCEA")
```

### Usage

It is easy to create an interactive QC report for those who possess little or no programming language skills. To run and generate an interactive QC report on your computer please install and call the scQCEA using rStudio, select all scripts incluidng `GenerateInteractiveQCReport()` function, and click on the "Run" button at the top right of the Source tab. An interactive QC report automatically will be generated in one HTML file, including four sections: experimental workflow, data processing workflow, sample information and QC metrics, data analysis and quality control.

```{r,eval=FALSE}

#########################################################################
# Please execute the code in the RStudio IDE (https://www.rstudio.com/) #
#########################################################################

library("scQCEA")
InputDir=system.file("extdata", package = "scQCEA")
GenerateInteractiveQCReport(InputDir)

############################################################ 
#  Find the "Interactive QC Report" in the Outputs/ folder #
############################################################

```

By default, the HTML report will be written in /Outputs directory named `CLICK_ME.html`. You can open `CLICK_ME.html` without using rStudio/R. In addition, you can find a zip file in the /Outputs directory which is particularly useful to share or store the QC reports. The content of the "Data processing Workflow" section is automatically adjusted based on the type of application (s) and the "Library Type" column in "samples.metadata" file.

Please see the [`manual`](https://isarnassiri.github.io/scQCEA/) for the usage of scQCEA to run the pipeline on your own data.

### Cell Type Enrichment Analysis
Cell type annotation on scRNA-Seq data is a pre-step for generating an interactive QC report with scQCEA. This step requires some bioinformatics efforts, but scQCEA provides `CellTypeEnrichment()` functions, for cell-type enrichment analysis at the single-cell level that comprises all the intermediate steps including visualization:

```{r,eval=FALSE}

##### Cell Type Enrichment Analysis #####

library("scQCEA")

csQCEAdir <- system.file("extdata", package = "scQCEA")
# A directory path incluidng input files/folders

DataTyep <- '10X-gex'
# Name of a folder including input files

SampleName <- '481207_03' 
# Name of an indicated sample

SamplesMetadata = paste(csQCEAdir, 'Inputs/samples.metadata', sep = '/' )
# Metadata of samples including the following headers: Project Number,	LIMS ID,	Sample Name,	Index	Library Type,	Genome,	Flowcell ID,	Lane Number,	Sequencing ID

ReadCount = paste(csQCEAdir, 'Inputs', DataTyep, SampleName, 'outs', sep = '/')
# Gene-cell count matrix from 10X CellRanger count

GTF = paste(csQCEAdir, 'ensembl_human.txt', sep = '/')
# We convert Ensembl ids to gene names/symbols by parsing this GTF (General Transfer Format) file

BackendDataDir = paste(csQCEAdir, 'ReferenceGeneSets/', sep = '/')
# We used Human Protein Atlas database (version 21.0) to generate a repository of reference gene sets that are exclusively expressed in each cell type

tSNECellranger = paste(csQCEAdir, 'Inputs', DataTyep, SampleName, '/outs/analysis/tsne/gene_expression_2_components', sep = '/')
# tSNE projections from 10X CellRanger count

UMAPCellranger =  paste(csQCEAdir, 'Inputs', DataTyep, SampleName, '/outs/analysis/umap/gene_expression_2_components', sep = '/')
# UMAP projections from 10X CellRanger count

RawFeatureDir = paste(csQCEAdir, 'Inputs', DataTyep, SampleName, 'outs/raw_feature_bc_matrix', sep = '/')
# A folder including raw feature-barcode matrices from 10X CellRanger count (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)

FilteredFeatureBarcodes = paste(csQCEAdir, 'Inputs', DataTyep, SampleName, 'outs/filtered_feature_bc_matrix', sep = '/')
# A folder including raw feature-barcode matrices from 10X CellRanger count (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)

CellTypeEnrichment(SampleName, SamplesMetadata, ReadCount, GTF, BackendDataDir, tSNECellranger, UMAPCellranger, RawFeatureDir, FilteredFeatureBarcodes ) 
``` 

`GenerateInteractiveQCReport()` function uses these output files and generates an interactive QC report for multiple samples to compare and examine biases and outliers over biological and technical measures.

### Citation
Isar Nassiri, Benjamin Fairfax, Andrew J Kwok, Aneesha Bhandari, Katherine Bull, Angela Lee, Yanxia Wu, Julian Knight, David Buck, Paolo Piazza. Demultiplexing of Single Cell RNA Sequencing Data using Interindividual Variation in Gene Expression. 
