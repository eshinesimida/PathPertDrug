# PathPertDrug:a software package for drug repurposing in cancer based on pathway perturbation

PathPertDrug contains the basic R functions and sample data for running the PathPertDrug algorithm. After installing and loading the package, users will be able to explore the framework of PathPertDrug


library(devtools) 
devtools::install_github('eshinesimida/PathPertDrug')

# More about PathPertDrug
The PathPertDrug is a novel tool used to identify the potentail candidate drugs for cancers based on pathway perturbation.

# Getting started

Step 1. Install the package
The dependency EnrichmentBrowser, KEGGdzPathwaysGEO,KEGGandMetacoreDzPathwaysGEO, and SPIA are unavailable on the CRAN but avaailable on BioConductor. So we need to install the BioManager manually.
