
#' Get disease-related data
#'
#' @return Fold change of DEGs
#' @export
#'
#' @examples
#' DE_Colorectal <- get_data1()
get_data1 <- function(){
  library(EnrichmentBrowser)
  library(KEGGdzPathwaysGEO)
  library(KEGGandMetacoreDzPathwaysGEO)
  library(SPIA)
  #----------------------------------------
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
  DE_Colorectal

}


#' Get disease-related data
#'
#' @return All gene
#' @export
#'
#' @examples
#' ALL_Colorectal <- get_data2()
get_data2 <- function(){
  library(EnrichmentBrowser)
  library(KEGGdzPathwaysGEO)
  library(KEGGandMetacoreDzPathwaysGEO)
  library(SPIA)
  #----------------------------------------
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
  ALL_Colorectal

}


