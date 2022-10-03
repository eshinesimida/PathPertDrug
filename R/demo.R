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

#----3.drug-related genes in three cell lines
drug_genes <- read.csv('drug_gene_all1.csv', header = T, stringsAsFactors = F)

#---------------------MCF,regard threhold as fdr<0.1 ---
drug_genes_MCF7 <- drug_genes[which(drug_genes$cell_line == 'MCF7'),]
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
  #drug-perturbed genes enrichment analysis
  res1=spia(de=DE,all=ALL,organism="hsa",nB=2000,plots=FALSE,beta=NULL,
            combine="fisher",verbose=TRUE)
  res_disease <- res[,c(2,11)]
  res_disease$status_disease <- ifelse(res_disease$Status=='Activated',1,-1)
  res_drug <- res1[,c(2,11)]
  res_drug$status_drug <- ifelse(res_drug$Status=='Activated',1,-1)
  B <- merge(res_disease, res_drug, by='ID')
  scores <- c(scores,-sum(B$status_disease*B$status_drug))
}
dataframe_col3 <- data.frame(drugs=drugs, scores = scores)
write.csv(dataframe_col3,'drug_col_sig.csv')


#reverse score between drug and disease

#disease-related genes
library(hgu133plus2.db)
x <- hgu133plus2ENTREZID
top$ENTREZ<-unlist(as.list(x[top$ID]))
top<-top[!is.na(top$ENTREZ),]
top<-top[!duplicated(top$ENTREZ),]
tg1<-top[top$adj.P.Val<0.1,]
DE_Colorectal=tg1$logFC
names(DE_Colorectal)<-as.vector(tg1$ENTREZ)
ALL_Colorectal=top$ENTREZ

#disease-related genes enrichment analysis
res=spia(de=DE_Disease,all=ALL_Disease,organism="hsa",nB=2000,plots=FALSE,beta=NULL,combine="fisher",verbose=TRUE)

#drug-perturbed genes
drug_genes <- read.csv('drug_gene_all1.csv', header = T, stringsAsFactors = F)
drug_genes_MCF7 <- drug_genes[which(drug_genes$cell_line == 'MCF7'),]
drug_gene <- drug_genes_MCF7[which(drug_genes_MCF7$drug == 'metformin'),]
DE <- drug_gene$fc
names(DE) <- drug_gene$ENTREZID
ALL= unique(drug_genes$ENTREZID)


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

score <- RS(res = res, DE_Drug=DE, ALL_Drug=ALL)


