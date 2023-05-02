#' Title
#'
#' @param res disease-related pathways
#' @param drug_name drug name
#' @param all all gene of drug
#'
#' @return reverse score of drug-disease
#' @export
#'
#' @examples
#' score  <- RS(res, de, all = ALL)
#'
RS <- function(res, drug_name, all){
  #drug-perturbed genes enrichment analysis
  drug_gene <- drug_exp[which(drug_exp$drug == drug_name),]
  DE <- drug_gene$fc
  names(DE) <- drug_gene$ENTREZID
  res1<-spia(de=DE,all=ALL,organism="hsa",nB=2000,plots=FALSE,
             beta=NULL,combine="fisher",verbose=TRUE)
  res_disease <- res[,c(2,11,13)]
  res_disease$status_disease <- ifelse(res_disease$Status=='Activated',1,-1)
  res_disease$status_disease1 <- res_disease$status_disease*res_disease$importance

  res_drug <- res1[,c(2,11)]
  res_drug$status_drug <- ifelse(res_drug$Status=='Activated',1,-1)
  B <- merge(res_disease, res_drug, by='ID')
  #calculate reverse score
  score <- sum(B$status_disease1*B$status_drug)
  score
}
