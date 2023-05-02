#------------.drug-related genes in three cell lines

#' Get DEGs of drug
#'
#' @return DEGs of drug
#' @export
#'
#' @examples
#' drug_exp <- get_drug_data1()
get_drug_data1 <- function(){
  drug_genes <- read.csv('drug_gene_all1.csv', header = T, stringsAsFactors = F)

  #---------------------MCF, fdr<0.01, fc>1---
  drug_genes_MCF7 <- drug_genes[which(drug_genes$cell_line == 'MCF7'),]

  drug_genes_MCF7_1 <- drug_genes_MCF7[drug_genes_MCF7$fdr<0.05,]
  ALL= unique(drug_genes$ENTREZID)
  drug_exp <- drug_genes_MCF7_1
  drug_exp

}

#' Get all gene of drug
#'
#' @return all gene of drug
#' @export
#'
#' @examples
#' ALL <- get_drug_data2()
get_drug_data2 <- function(){
  drug_genes <- read.csv('drug_gene_all1.csv', header = T, stringsAsFactors = F)
  ALL <- unique(drug_genes$ENTREZID)

  ALL

}



