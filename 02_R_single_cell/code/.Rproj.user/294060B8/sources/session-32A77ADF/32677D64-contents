## Creating the SBC signature by merging DEGs among SBC cells vs. remainder basal common to both ACI&SD
## Generating ranked list files for GSEA analysis of the enrichment of SBC signatures in high-risk signatures
## Generating ranked list files for GSEA analysis of the enrichment of PROCR+ Wang2015 signature in the SBC signature, as well as in SBC subclusters
## Author: Simona Cristea
## Date: October 2022
## Publication: Breast cancer prevention by short-term inhibition of TGFβ signaling, Nat Comm 2022

## load libraries 
library(tidyverse)
library(readxl)

###########################################################################
###########################################################################
###                                                                     ###
###                            SBC SIGNATURE                            ###
###                                                                     ###
###########################################################################
###########################################################################
secretory.sd4 <- read.table(file = "../../tables/secretory-sd4.txt", sep = "\t")
secretory.sd4$gene <- rownames(secretory.sd4)
secretory.sd4$rank <- c(1:nrow(secretory.sd4))

secretory.aci <- read.table(file = "../../tables/secretory-aci.txt", sep = "\t")
secretory.aci$gene <- rownames(secretory.aci)
secretory.aci$rank <- c(1:nrow(secretory.aci))


combined <- full_join(secretory.sd4, secretory.aci, by = "gene", suffix = c(".sd", ".aci"))
combined$avg_strains_log2FC <- (combined$avg_log2FC.aci + combined$avg_log2FC.sd)/2

combined <- as_tibble(combined)



############################################################################
############################################################################
###                                                                      ###
###                      GSEA RISK ASSESSMENT FILES                      ###
###                                                                      ###
############################################################################
############################################################################
####### read in risk assessment lists
risk.assessment <- list()
risk.assessment$PvsNP <- read_excel("../../tables/gene-lists/Risk assessment_Masa.xlsx", sheet = 1)
risk.assessment$BRCA1 <- read_excel("../../tables/gene-lists/Risk assessment_Masa.xlsx", sheet = 2)
risk.assessment$BRCA2 <- read_excel("../../tables/gene-lists/Risk assessment_Masa.xlsx", sheet = 3)
risk.assessment <- lapply(risk.assessment, function(x){
  colnames(x) <- x[3,]
  x <- x[-c(1:3),]
  return(x)
})

#p.val.thresh <- 0.05
q.val.thresh <- 0.1
risk.assessment$BRCA1 <- risk.assessment$BRCA1[which(as.numeric(risk.assessment$BRCA1$`q-value`) < q.val.thresh),]
risk.assessment$BRCA2 <- risk.assessment$BRCA2[which(as.numeric(risk.assessment$BRCA2$`q-value`) < q.val.thresh),]

risk.assessment$BRCA1pos <- risk.assessment$BRCA1[which(as.numeric(risk.assessment$BRCA1$PseudoFC) > 0),]
risk.assessment$BRCA1pos <- risk.assessment$BRCA1pos[order(as.numeric(risk.assessment$BRCA1pos$PseudoFC), decreasing = TRUE),]

risk.assessment$BRCA1neg <- risk.assessment$BRCA1[which(as.numeric(risk.assessment$BRCA1$PseudoFC) <= 0),]
risk.assessment$BRCA1neg <- risk.assessment$BRCA1neg[order(as.numeric(risk.assessment$BRCA1neg$PseudoFC)),]

risk.assessment$BRCA2pos <- risk.assessment$BRCA2[which(as.numeric(risk.assessment$BRCA2$PseudoFC) > 0),]
risk.assessment$BRCA2pos <- risk.assessment$BRCA2pos[order(as.numeric(risk.assessment$BRCA2pos$PseudoFC), decreasing = TRUE),]

risk.assessment$BRCA2neg <- risk.assessment$BRCA2[which(as.numeric(risk.assessment$BRCA2$PseudoFC) <= 0),]
risk.assessment$BRCA2neg <- risk.assessment$BRCA2neg[order(as.numeric(risk.assessment$BRCA2neg$PseudoFC)),]


risk.assessment$PvsNPpos <- risk.assessment$PvsNP[which(as.numeric(risk.assessment$PvsNP$PseudoFC) > 0),]
risk.assessment$PvsNPpos <- risk.assessment$PvsNPpos[order(as.numeric(risk.assessment$PvsNPpos$PseudoFC), decreasing = TRUE),]

risk.assessment$PvsNPneg <- risk.assessment$PvsNP[which(as.numeric(risk.assessment$PvsNP$PseudoFC) <= 0),]
risk.assessment$PvsNPneg <- risk.assessment$PvsNPneg[order(as.numeric(risk.assessment$PvsNPneg$PseudoFC)),]


write.table(file = "../../tables/gene-lists/BRCA1pos.txt", risk.assessment$BRCA1pos, quote = FALSE, sep = "\t")
write.table(file = "../../tables/gene-lists/BRCA1neg.txt", risk.assessment$BRCA1neg, quote = FALSE, sep = "\t")
write.table(file = "../../tables/gene-lists/BRCA2pos.txt", risk.assessment$BRCA2pos, quote = FALSE, sep = "\t")
write.table(file = "../../tables/gene-lists/BRCA2neg.txt", risk.assessment$BRCA2neg, quote = FALSE, sep = "\t")
write.table(file = "../../tables/gene-lists/PvsNPpos.txt", risk.assessment$PvsNPpos, quote = FALSE, sep = "\t")
write.table(file = "../../tables/gene-lists/PvsNPneg.txt", risk.assessment$PvsNPneg, quote = FALSE, sep = "\t")



###########################################################################
###########################################################################
###                                                                     ###
###                           SBC SUBCLUSTERS                           ###
###                                                                     ###
###########################################################################
###########################################################################
## creating the SBC subclusters
rat <- readRDS(file = "../../data/rdata/rat-initial-sd.RDS")
rat.fractions.detailed <- SplitObject(rat, split.by = "celltype.detailed") 
rat <- rat.fractions.detailed$secretory


## renormalizing
rat <- SCTransform(rat, verbose = TRUE)
rat <- RunPCA(rat, features = VariableFeatures(object = rat, assay = "SCT"), assay = "SCT")
rat <- FindNeighbors(rat, dims = 1:30, verbose = TRUE, assay = "SCT", reduction = "pca", graph.name = "SNN.SCT")
rat <- FindClusters(rat, verbose = TRUE, graph.name = "SNN.SCT")
rat <- RunUMAP(rat, dims = 1:30, verbose = TRUE, reduction.name = "UMAP.SCT", assay = "SCT")


## clustree to understand impact of resolution on clustering
for (res in seq(from = 0.1, to = 1, by = 0.1)){
  rat <- FindClusters(rat, verbose = TRUE, graph.name = "SNN.SCT", resolution = res)
  rat <- RunUMAP(rat, dims = 1:30)
}
pdf("../../plots/sd-sorted/secretory-clustree-umap.pdf", width = 12, height = 10)
clustree::clustree(rat, prefix = "SNN.SCT_res.")
dev.off()


## integrate
options(future.globals.maxSize = 4000 * 1024^2)
rat.list <- SplitObject(rat, split.by = "animal")
rat.list <- lapply(rat.list, FUN = function(x) {
  x <- SCTransform(x, verbose = TRUE)
})
rat.features <- SelectIntegrationFeatures(object.list = rat.list, nfeatures = 3000)
rat.list <- PrepSCTIntegration(object.list = rat.list, anchor.features = rat.features, verbose = TRUE)
rat.anchors <- FindIntegrationAnchors(object.list = rat.list, normalization.method = "SCT",
                                      anchor.features = rat.features, verbose = TRUE, dims = 1:5,
                                      k.anchor = 5, k.filter = 5, k.score = 5)
integrated <- IntegrateData(anchorset = rat.anchors, normalization.method = "SCT", verbose = TRUE, dims = 1:10, k.weight = 5)
integrated <- RunPCA(integrated, verbose = TRUE)
integrated <- FindNeighbors(integrated, dims = 1:10, verbose = FALSE, reduction = "pca", graph.name = "SNN.INTEGRATE")
integrated <- FindClusters(integrated, verbose = TRUE, graph.name = "SNN.INTEGRATE", resolution = 0.6)
integrated <- RunUMAP(integrated, dims = 1:10)


## clustree on the integrated data
for (res in seq(from = 0.1, to = 1, by = 0.1)){
  integrated <- FindClusters(integrated, verbose = TRUE, graph.name = "SNN.INTEGRATE", resolution = res)
  integrated <- RunUMAP(integrated, dims = 1:30)
  pdf(paste0("../../plots/sd-sorted/secretory-integrated-umap-res", res, ".pdf"))
  print(DimPlot(integrated))
  dev.off()
}
pdf(paste0("../../plots/sd-sorted/secretory-integrated-clustree-umap.pdf"), width = 12, height = 10)
clustree::clustree(integrated, prefix = "SNN.INTEGRATE_res.")
dev.off()


## DEG tables for resolution 0.2 for integrated
DefaultAssay(integrated) <- "RNA"
rat <- NormalizeData(integrated, assay = "RNA")
Idents(integrated) <- integrated$SNN.INTEGRATE_res.0.2
markers.clust0 <- FindMarkers(integrated, ident.1 = 0, ident.2 = NULL, assay = "RNA", slot = "data")
write.table(markers.clust0, file = "../../tables/gene-lists/SBCs-SD-cluster0.txt", quote = FALSE, sep = "\t")
markers.clust1 <- FindMarkers(integrated, ident.1 = 1, ident.2 = NULL, assay = "RNA", slot = "data")
write.table(markers.clust1, file = "../../tables/gene-lists/SBCs-SD-cluster1.txt", quote = FALSE, sep = "\t")
markers.clust2 <- FindMarkers(integrated, ident.1 = 2, ident.2 = NULL, assay = "RNA", slot = "data")
write.table(markers.clust2, file = "../../tables/gene-lists/SBCs-SD-cluster2.txt", quote = FALSE, sep = "\t")
markers.clust3 <- FindMarkers(integrated, ident.1 = 3, ident.2 = NULL, assay = "RNA", slot = "data")
write.table(markers.clust3, file = "../../tables/gene-lists/SBCs-SD-cluster3.txt", quote = FALSE, sep = "\t")


## markers to plot
additional.markers <- c("SPARCL1", "ACTA2", "LALBA", "ESR1", "KRT5", "KRT7", "PGR", "KRT15", "KRT17", "COL1A1", 
                        "MKI67", "CDH1", "LUM","ID3", "PROCR", "CDH5", "TWIST1")
for (m in additional.markers){
  pdf(paste0("../../plots/sd-sorted/secretory-", m, ".pdf"))
  print(FeaturePlot(integrated, features = m))
  dev.off()
}

rat <- integrated
saveRDS(rat, file = "../../data/rdata/sd4-secretory.RDS")
