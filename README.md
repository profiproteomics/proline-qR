
# Proline qR

## Overview

Proline-qR is a set of R scripts and Rmarkdown reports for generating reports of controlled quantitation datasets from [Proline](http://www.profiproteomics.fr/proline/). Some data are included but you can provide your own data (exported as xlsx document from Proline) and tell Proline-qR to generate the appropriate report depending on the reference quantitation datasets.

The generated reports and figures are broadly inspired by the figures shown in [Proline paper](https://doi.org/10.1093/bioinformatics/btaa118) (Bouyssié,D. et al. (2020) Proline: an efficient and user-friendly software suite for large-scale proteomics. Bioinformatics, 36, 3148–3155). 

The supported datasets:

+ UPS: is a proteomic standard composed of an equimolarmixture of 48 human proteins (Sigma UPS1) spiked at different concentrations into a background of yeast cell lysate. It has been described in Ramus,C. et al. (2016) Benchmarking quantitative label-free LC–MS data processing workflows using a complex spiked proteomic standard dataset. Journal of Proteomics, 132, 51–62 and in Bouyssié,D. et al. (2020) Proline: an efficient and user-friendly software suite for large-scale proteomics. Bioinformatics, 36, 3148–3155.

  + UPS_10Concentrations: 10 different concentrations

  + UPS_4Concentrations: a subset of 4 concentrations 

+ ABRF: is the dataset from the ABRF 2015 Study described in Choi,M. et al. (2017) ABRF Proteome Informatics Research Group (iPRG) 2015 Study: Detection of Differentially Abundant Proteins in Label-Free Quantitative LC–MS/MS Experiments. J. Proteome Res., 16, 945–957.

+ PME12: is a EuPA Standardization Initiative. Proteomics Multicentric Experiment 12 (PME12) was launched last April 2019 with the aim of comparing the performance of different label free quantification methods/pipelines. HeLa, Saccharomyces cerevisiae and Escherichia coli tryptic digests were purchased, resuspended in milli Q water with 0.1% FA and sonicated for optimal recovery. Peptide extracts were then mixed in different proportions as described by Navarro et al. (Nature Biotech 2016).


## Installation

Proline-qR has been tested with R version 3.6.0. The following packages/version are required: 

```r
rmdformats_0.3.5 matrixcalc_1.0-3 reshape2_1.4.3   ggplot2_3.2.0    dplyr_1.0.0      openxlsx_4.1.0.1 rmarkdown_1.13   knitr_1.23
```

In the current version the code is delivered as an archive of the source project. To install Proline-qR uncompress the archive on your disk and open `Proline-qR.Rproj` with RStudio. Once opened load the required library and generate a report from a folder containing an xlsx export from Proline:      

```r
# source script files 
source('R/file-utils.R')
source('R/commons_functions.R')
source('R/generate-reports.R')

```

## Usage

To use Proline-qR start by loading the required dependencies

```r
# load required libraries

library(knitr)
library(rmarkdown)
library(dplyr)
library(ggplot2)
library(openxlsx)

```

The main functions are located in the file `R/generate-reports.R`. The following line generate in the `output`folder a report for the PME12 dataset data from the Proline xlsx that are in the supplied `data` folder. The search rule for relative forlder path is to search in the `data/<dataset_name>/<relative folder path>` folder. 

```r
generate_dataset_report("PME12", "20191216/spec_bestion_sum")

```

An absolute filepath can also be provided. In this case the report is generate into the same folder. 

```r
generate_dataset_report("PME12", "C:/Local/tmp/20191216/spec_bestion_sum")

```


## Getting help

To get help send a email to the authors.