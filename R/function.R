#function of reverse function
RS <- function(res, DE_Drug, ALL_Drug){
  #drug-perturbed genes enrichment analysis
  res1=spia(de=DE_Drug,all=ALL_Drug,organism="hsa",nB=2000,plots=FALSE,beta=NULL,combine="fisher",verbose=TRUE)
  res_disease <- res[,c(2,11)]
  res_disease$status_disease <- ifelse(res_disease$Status=='Activated',1,-1)
  res_drug <- res1[,c(2,11)]
  res_drug$status_drug <- ifelse(res_drug$Status=='Activated',1,-1)
  B <- merge(res_disease, res_drug, by='ID')
  #calculate reverse score
  score <- -sum(B$status_disease*B$status_drug)
  score
}
