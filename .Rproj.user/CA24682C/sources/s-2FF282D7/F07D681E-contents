
#' enrichment analysis of CRC
#'
#' @param de,FC of differentially expressed genes
#' @param all, all gene ID
#'
#' @return
#' @export
#'
#' @examples
get_pathway_disease <- function(de, all){
  library(SPIA)
  # pathway analysis based on combined evidence; # use nB=2000 or more for more accurate results
  res=spia(de=DE_Colorectal,all=ALL_Colorectal,organism="hsa",nB=2000,plots=FALSE,beta=NULL,
           combine="fisher",verbose=TRUE)
  res <- res[res$pG <0.05,]
  res$importance <- -log(res$pG)
  res

}

