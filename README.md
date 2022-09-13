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
3. Download scQCEA from GitHub [(LINK)](https://github.com/isarnassiri/scQCEA/), and unzip the folder
4. To install scQCEA, run the `RUN_ME.R` script from the RStudio. All dependency packages automatically will be downloaded, installed and loaded from CRAN-like repositories.

### Usage

It is easy to create an interactive QC report for those who possess little or no programming language skills. To run and generate an interactive QC report on your computer please open the `RUN_ME.R` file using rStudio, select all scripts incluidng `GenerateInteractiveQCReport()` function, and click on the "Run" button at the top right of the Source tab. An interactive QC report automatically will be generated in one HTML file, including four sections: experimental workflow, data processing workflow, sample information and QC metrics, data analysis and quality control.

```{r,eval=FALSE}

#########################################################################
# Please execute the code in the RStudio IDE (https://www.rstudio.com/) #
#########################################################################

##### Install and load R packages #####
if(!("rstudioapi" %in% installed.packages()[,"Package"])) install.packages("rstudioapi", repos = "http://cran.us.r-project.org", update = FALSE); library("rstudioapi")

##### Generate an "Interactive QC Report" #####
setwd("~/"); setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), '/Scripts/')); 
source("GenerateInteractiveQCReport.R")
GenerateInteractiveQCReport()

############################################################ 
#  Find the "Interactive QC Report" in the Outputs/ folder #
############################################################

```

By default, the HTML report will be written in /Outputs directory named `CLICK_ME.html`. You can open `CLICK_ME.html` without using rStudio/R. In addition, you can find a zip file in the /Outputs directory which is particularly useful to share or store the QC reports. The content of the "Data processing Workflow" section is automatically adjusted based on the type of application (s) and the "Library Type" column in "samples.metadata" file.

### Cell Type Enrichment Analysis
Cell type annotation on scRNA-Seq data is a pre-step for generating an interactive QC report with scQCEA. This step requires some bioinformatics efforts, but scQCEA provides `CellTypeEnrichment()` functions, for cell-type enrichment analysis at the single-cell level that comprises all the intermediate steps including visualization (you can find the code in `RUN_ME.R` file):

```{r,eval=FALSE}

##### Cell Type Enrichment Analysis #####
setwd("~/"); setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), '/Scripts/')); 
source("CellTypeEnrichment.R")
CellTypeEnrichment()

``` 

`GenerateInteractiveQCReport()` function uses these output files and generates an interactive QC report for multiple samples to compare and examine biases and outliers over biological and technical measures.

### Citation

Isar Nassiri, Benjamin Fairfax, Angela Lee, Yanxia Wu, David Buck, Paolo Piazza. scQCEA: A Framework for Annotation and Quality Control Report of Single-Cell RNA-Sequencing Data. 
