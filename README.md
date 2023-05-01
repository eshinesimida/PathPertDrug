# PathPertDrug:a software package for drug repurposing in cancer based on pathway perturbation

PathPertDrug contains the basic R functions and sample data for running the PathPertDrug algorithm. After installing and loading the package, users will be able to explore the framework of PathPertDrug


# More about PathPertDrug
The PathPertDrug is a novel tool used to identify the potentail candidate drugs for cancers based on pathway perturbation.

# Getting started

## Step 1. Install the package
The dependency `EnrichmentBrowser`, `KEGGdzPathwaysGEO`, `KEGGandMetacoreDzPathwaysGEO` and `SPIA` are unavailable on the CRAN but available on [BioConductor](https://www.bioconductor.org/). So we need to install the BiocManager manually. 

``` r
if (!"BiocManager" %in% as.data.frame(installed.packages())$Package)
  install.packages("BiocManager")
BiocManager::install(c("EnrichmentBrowser", "KEGGdzPathwaysGEO","KEGGandMetacoreDzPathwaysGEO","SPIA"))
```
Then you can install the development version of PathPerDrug from [GitHub](https://github.com/) with:

``` r
if (!"devtools" %in% as.data.frame(installed.packages())$Package)
  install.packages("devtools")
devtools::install_github("eshinesimida/PathPertDrug")

### Examples

Below is a basic example that shows how to solve a common problem:

``` r
library(PathPertDrug)

# Load the data
data(example_disease_list)
data(example_drug_target_list)
data(example_ppi)

# Perform a simple PathPertDrug analysis

result <- DTSEA(
    network = example_ppi,
    disease = example_disease_list,
    drugs = example_drug_target_list
)
```

You can enable the multicore feature to utilize the multicore advantages. Here is the benchmark. 

``` r
