
#' Calculate drug-disease reverse scores based on pathway perturbation
#'
#' This function uses the previous SPIA method and integrate the change of
#' of genes Pearson coefficient(PCC) from two groups. We proposed a set of three
#' pathway analysis methods based on the change of PCC. We applied these approaches
#' to colorectal cancer, lung cancer and Alzheimer's disease datasets and so on.
#'
#'
#' We used a compendium of 22 GEO datasets obtained from the KEGGdzPathwaysGEO
#' and the KEGGandMetacoreDzPathwaysGEO benchmark sets(Tarca et al., 2012; Tarca et al., 2013)
#' in this study. These datasets have been specifically chosen as they study a
#' certain human disease in a corresponding KEGG pathway (e.g. Alzheimer's disease).
#' These pathways are regarded as the target pathways in the following. We investigate
#' first how well the individual set- and network-based methods detect the target
#'
#' The expression data sets involved 11 conditions and 17 tissues.
#'
#'
#' @param de   The number of differential genes
#' @param all  All genes in human
#' @param normal the number of normal samples
#' @param tumor the number of tumor samples
#' @param flag flag = 1,0,-1 , if flag = 1 from normal to tumor, flag = -1 from tumor to normal, flag = 0 stand for absolute value between tow groups
#' @export
#' @examples
#' #import EnrichmentBrowser, KEGGandMetacoreDzPathwaysGEO, KEGGdzPathwaysGEO and SPIA package
#' library(EnrichmentBrowser)
#' library(KEGGandMetacoreDzPathwaysGEO)
#' library(SPIA)
#' data("GSE8671")
#' # Get expression profile of GSE8671
#' exprs_all <- exprs(GSE8671)
#' # Add the gene symbol
#' all.eset <- probe.2.gene.eset(GSE8671)
#' head(featureNames(all.eset))
#' before.norm <- exprs(all.eset)
#' # Gene normalization
#' all.eset <- normalize(all.eset, norm.method="quantile")
#' after.norm <- exprs(all.eset)
#' exprs_all1 <- data.frame(after.norm)
#' table(pData(all.eset)$Group)
#' pData(all.eset)$GROUP <- ifelse(pData(all.eset)$Group == "d", 1, 0)
#' normal <- length(which(pData(all.eset)$GROUP == '0'))
#' tumor <- length(which(pData(all.eset)$GROUP == '1'))
#' # Get the differential expression genes in limmar package
#' all.eset <- de.ana(all.eset)
#' head(fData(all.eset), n=4)
#' all_de <- fData(all.eset)
#' tg <- all_de[all_de$ADJ.PVAL < 0.1,]
#' DE_colorectal = tg$FC
#' names(DE_colorectal)<-as.vector(rownames(tg))
#' ALL_colorectal = rownames(all_de)
#' #The result of spia method
#' res_spia = spia(de = DE_colorectal, all=ALL_colorectal, organism="hsa",nB=2000,plots=FALSE,beta=NULL,combine="fisher",verbose=TRUE)
#' gse_madat2 <- exprs_all1
#' # The results of spia_nt method
#' res_nt = spiap(de=DE_colorectal, all=ALL_colorectal, gse_madat2 = gse_madat2,normal = normal, tumor = tumor,norganism="hsa",nB=2000,plots=FALSE,
#'                   beta=NULL,combine="fisher",verbose=T, flag = 1)
#' #The results of spia_tn method
#' res_tn = spiap(de=DE_colorectal, all=ALL_colorectal, gse_madat2 = gse_madat2, normal = normal, tumor = tumor,organism="hsa",nB=2000,plots=FALSE,
#'                   beta=NULL,combine="fisher",verbose=T, flag = -1)
#' #The results of spia_abs method
#' res_abs = spiap(de = DE_colorectal, all=ALL_colorectal , gse_madat2 = gse_madat2,normal = normal, tumor = tumor,organism="hsa",nB=2000,plots=FALSE,
#'                   beta=NULL,combine="fisher",verbose=T, flag = 0)

#spiap <- function(de = NULL, all = NULL,gse_madat2 = gse_madat2, normal = NULL, tumor = NULL, organism = "hsa", data.dir = NULL,
#                  pathids = NULL, nB = 2000, plots = FALSE, verbose = TRUE,
#                  beta = NULL, combine = "fisher"){
# UseMethod("spiap")
#}

#
# This is an example function named 'colorectal cancer'
# which prints 'candidate drugs of CRC'.

#disease-related pathways based on SPIA method
library(SPIA)
data(colorectalcancer)
options(digits=3)
head(top)

#1.obtain differentially expressed genes of colorectal cancer
library(hgu133plus2.db)
x <- hgu133plus2ENTREZID
top$ENTREZ<-unlist(as.list(x[top$ID]))
top<-top[!is.na(top$ENTREZ),]
top<-top[!duplicated(top$ENTREZ),]
tg1<-top[top$adj.P.Val<0.1,]
DE_Colorectal=tg1$logFC
names(DE_Colorectal)<-as.vector(tg1$ENTREZ)
ALL_Colorectal=top$ENTREZ

#2.enrichment analysis of CRC
# pathway analysis based on combined evidence; # use nB=2000 or more for more accurate results
res=spia(de=DE_Colorectal,all=ALL_Colorectal,organism="hsa",nB=2000,plots=FALSE,beta=NULL,
         combine="fisher",verbose=TRUE)


res <- res[res$pG <0.05,]
res$importance <- -log(res$pG)

#----3.drug-related genes in three cell lines
drug_genes <- read.csv('drug_gene_all1.csv', header = T, stringsAsFactors = F)

#---------------------MCF,regard threhold as fdr<0.1 ---
drug_genes_MCF7 <- drug_genes[which(drug_genes$cell_line == 'MCF7'),]

ALL= unique(drug_genes$ENTREZID)
#res <- res3
drugs <- unique(drug_genes_MCF7$drug)
pb <- txtProgressBar(min = 0, max = length(drugs), style = 3)
scores <- c()
for(ii in 1:length(drugs)){
  setTxtProgressBar(pb, ii)
  i1 <- drugs[ii]
  cat('i=',i1,'\n')
  drug_gene <- drug_genes_MCF7[which(drug_genes_MCF7$drug == i1),]
  DE <- drug_gene$fc
  names(DE) <- drug_gene$ENTREZID
  #reverse score between drug and disease
  score <- RS(res = res, DE_Drug=DE, ALL_Drug=ALL)
  scores <- c(scores,score)
}
dataframe_col3 <- data.frame(drugs=drugs, scores = scores)
write.csv(dataframe_col3,'drug_col_sig.csv')




