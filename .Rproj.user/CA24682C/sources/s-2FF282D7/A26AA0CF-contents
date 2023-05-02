
library(KEGGdzPathwaysGEO)
library(KEGGandMetacoreDzPathwaysGEO)
data('GSE8671')
plot_pathway_disease <- function(ExpData = GSE8671,pathwayID = 1){

  library(EnrichmentBrowser)
  library(SPIA)
  library(pathview)
  data(demo.paths)
  exprs_all <- exprs(ExpData)
  # Add the gene symbol
  all.eset <- probe2gene(ExpData)
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
  tg1 <- tg3
  tg1 <- as.matrix(tg1)
  pv.out <- pathview(gene.data = tg1[, 1], pathway.id = demo.paths$sel.paths[pathwayID],
                     species = "hsa", out.suffix = "gse8671", kegg.native = F,
                     sign.pos = demo.paths$spos[pathwayID], same.layer = T)

}


#plot significant pathway by drug
plot_pathway_drug <- function(drug = "metformin",pathwayID = 1){

  drug_gene <- drug_genes_MCF7_1[which(drug_genes_MCF7_1$drug == drug),]
  rownames(drug_gene) <- drug_gene$ENTREZID
  drug_gene$fc
  drug_gene1 <- drug_gene[,6:7]
  drug_gene1 <- as.matrix(drug_gene1)
  pv.out <- pathview(gene.data = drug_gene1[, 1], pathway.id = demo.paths$sel.paths[pathwayID],
                     species = "hsa", out.suffix = "gse8671_drug1", kegg.native = F,
                     sign.pos = demo.paths$spos[pathwayID], same.layer = T)

}

#drug structure

getDrug_structure <- function (drugname = "", main = "", sub = "")
{
  haveChemmineR <- PackageLoaded("ChemmineR")
  havervest <- PackageLoaded("rvest")
  if (haveChemmineR == FALSE) {
    stop("The 'ChemmineR' library, should be loaded first")
  }
  if (havervest == FALSE) {
    stop("The 'rvest' library, should be loaded first")
  }
  drugname <- unlist(strsplit(drugname, "\\("))[1]
  drugname1 <- tolower(drugname)
  Drugs_CID <- GetExample("Drugs_CID")
  drugCid <- Drugs_CID[which(Drugs_CID[, 1] == drugname1),
                       2]
  drug_url <- paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/",
                    drugCid, "/record/SDF/?record_type=2d&response_type=display",
                    sep = "")
  cw <- try(read_html(drug_url))
  if ("try-error" %in% class(cw)) {
    stop("Please ensure smooth network connection")
  }
  drugnr <- html_text(cw)
  drugnr <- strsplit(drugnr, "\n")
  drugnr <- unlist(drugnr)
  sdfset <- read.SDFset(drugnr)
  if (main == "") {
    sdfset@ID <- drugname
  }
  else {
    sdfset@ID <- main
  }
  return(sdfset)
}
