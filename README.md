# PathPertDrug:a software package for drug repurposing in cancer based on pathway perturbation

PathPertDrug contains the basic R functions and sample data for running the PathPertDrug algorithm. After installing and loading the package, users will be able to explore the framework of PathPertDrug


# More about PathPertDrug
The PathPertDrug is a novel tool used to identify the potentail candidate drugs for cancers based on pathway perturbation.

# Getting started

## Step 1. Pre run the method for installation

You should ensure that you have the necessary system dependencies configured.

For Windows (8.1 / 10 / 11): Rtools should be installed to the system path.

The latest base R is recommended. The compatibility of the earlier version (v4.0.x) is under evaluation.
We use R version is [64-bit] d:\Program Files\R\R-4.1.2

## Step 2. Install the package
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

```
## Examples

Below is a basic example that shows how to obtain candidate drugs of colorecatal cancer:

``` r
#Load require package
library(PathPertDrug)
library(EnrichmentBrowser)
library(KEGGdzPathwaysGEO)
library(KEGGandMetacoreDzPathwaysGEO)
library(SPIA)


#Load the data
#--------1.GSE8671 differnential analysis--
data("GSE8671")
exprs_all <- exprs(GSE8671)
# Add the gene symbol
all.eset <- probe2gene(GSE8671)
before.norm <- assay(all.eset)
# Gene normalization
all.eset <- normalize(all.eset, norm.method="quantile")
after.norm <- assay(all.eset)

exprs_all1 <- data.frame(after.norm)
table(colData(all.eset)$Group)
colData(all.eset)$GROUP <- ifelse(colData(all.eset)$Group == "d", 1, 0)
normal <- length(which(colData(all.eset)$GROUP == '0'))
tumor <- length(which(colData(all.eset)$GROUP == '1'))
# Get the differential expression genes in limmar package
all.eset <- deAna(all.eset, padj.method="BH")
all_de <- rowData(all.eset, use.names=TRUE)
all_de <- data.frame(all_de)
tg <- all_de[order(all_de$ADJ.PVAL, decreasing = FALSE),]

tg3 <- all_de[all_de$ADJ.PVAL<0.05,]


DE_Colorectal = tg3$FC
names(DE_Colorectal)<-rownames(tg3)
ALL_Colorectal <- rownames(tg)


# pathway analysis based on combined evidence; # use nB=2000 or more for more accurate results
# obtain disease-related significant pathways
res=spia(de=DE_Colorectal,all=ALL_Colorectal,organism="hsa",nB=2000,plots=FALSE,beta=NULL,
         combine="fisher",verbose=TRUE)


res <- res[res$pG <0.05,]
res$importance <- -log(res$pG)

#Load drug-induced differentially expressed genes

drug_genes <- read.csv('drug_gene_all1.csv', header = T, stringsAsFactors = F)
drug_genes_MCF7 <- drug_genes[which(drug_genes$cell_line == 'MCF7'),]
drug_genes_MCF7_1 <- drug_genes_MCF7[drug_genes_MCF7$fdr<0.05,]
ALL= unique(drug_genes$ENTREZID)
drugs <- unique(drug_genes_MCF7_1$drug)

pb <- txtProgressBar(min = 0, max = length(drugs), style = 3)

#This process need several hours to get the results of all drug-disease 
scores <- c()
for(ii in 1:length(drugs)){
  setTxtProgressBar(pb, ii)
  i1 <- drugs[ii]
  cat('i=',i1,'\n')
  drug_gene <- drug_genes_MCF7_1[which(drug_genes_MCF7_1$drug == i1),]
  DE <- drug_gene$fc
  names(DE) <- drug_gene$ENTREZID
  res1=spia(de=DE,all=ALL,organism="hsa",nB=2000,plots=FALSE,beta=NULL,
            combine="fisher",verbose=TRUE)
  res_disease <- res[,c(2,11,13)]
  res_disease$status_disease <- ifelse(res_disease$Status=='Activated',1,-1)
  res_disease$status_disease1 <- res_disease$status_disease*res_disease$importance
  res_drug <- res1[,c(2,11)]
  res_drug$status_drug <- ifelse(res_drug$Status=='Activated',1,-1)
  B <- merge(res_disease, res_drug, by='ID')
  scores <- c(scores,sum(B$status_disease1*B$status_drug))
}

# the reverse scores of drug-disease
dataframe_col3 <- data.frame(drugs=drugs, scores = scores)
```

