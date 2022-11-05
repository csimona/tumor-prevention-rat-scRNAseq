## Single cell analysis on whole mammary gland SD data, separately per cell type
## Author: Simona Cristea
## Date: October 2022
## Publication: Breast cancer prevention by short-term inhibition of TGFβ signaling, Nat Comm 2022


## load libraries
source(here::here("libraries.R"))

## read in data
rat <- readRDS(file = "../../data/rdata/rat-whole-gland.RDS")



###########################################################################
###########################################################################
###                                                                     ###
###                         SPLIT PER CELL TYPE                         ###
###                                                                     ###
###########################################################################
###########################################################################
rat.split <- SplitObject(rat, split.by = "celltype")
rat.split$`nerve-myelin` <- NULL
rat.split$vein <- NULL

rat.split <- lapply(rat.split, function(x){
  x <- SCTransform(x, verbose = TRUE)
  x <- RunPCA(x, features = VariableFeatures(object = x, assay = "SCT"), assay = "SCT")
  x <- FindNeighbors(x, dims = 1:30, verbose = FALSE, assay = "SCT", reduction = "pca", graph.name = "SNN.SCT")
  x <- FindClusters(x, verbose = TRUE, graph.name = "SNN.SCT")
  x <- RunUMAP(x, dims = 1:30, verbose = TRUE, reduction.name = "UMAP.SCT", assay = "SCT")
  return(x)
})

canonical.markers <- read.csv(file ="../../tables/gene-lists/canonical_markers.csv")
canonical.markers <- canonical.markers$marker



############################################################################
############################################################################
###                                                                      ###
###                                STROMA                                ###
###                                                                      ###
############################################################################
############################################################################
pdf("../../plots/whole-gland/stroma/umap-clusters.pdf")
DimPlot(rat.split$stroma, label = TRUE)
dev.off()
pdf("../../plots/whole-gland/stroma/umap-drug.pdf")
DimPlot(rat.split$stroma, group.by = "drug", cols = c("control" = "#C0C0C0", "drug" = "#FED242"))
dev.off()
pdf("../../plots/whole-gland/stroma/umap-animal.pdf")
DimPlot(rat.split$stroma, group.by = "sample")
dev.off()
pdf("../../plots/whole-gland/stroma/umap-animal.pdf")
DimPlot(rat.split$stroma, group.by = "sample")
dev.off()
pdf("../../plots/whole-gland/stroma/umap-celltype-ex-classification.pdf")
DimPlot(rat.split$stroma, group.by = "celltype.detailed")
dev.off()
pdf("../../plots/whole-gland/stroma/ridge.qc.pdf")
RidgePlot(rat.split$stroma, features = c("percent.mt", "nFeature_SCT"))
dev.off()



############################################################################
############################################################################
###                                                                      ###
###                         STROMA CELL SUBTYPES                         ###
###                                                                      ###
############################################################################
############################################################################
diff.expr.stroma <- FindAllMarkers(rat.split$stroma)
write.table(diff.expr.stroma, quote = FALSE, row.names = FALSE, sep = "\t",
            file = "../../tables/gene-lists/whole.gland.diff.expr.stroma.txt")

DefaultAssay(rat.split$stroma) <- "SCT"
pdf("../../plots/whole-gland/stroma/heatmap.pdf", height = 10, width = 20)
DoHeatmap(rat.split$stroma, features = canonical.markers)
dev.off()
pdf("../../plots/whole-gland/stroma/heatmap-zoom.pdf", height = 10, width = 20)
DoHeatmap(rat.split$stroma[,which(as.numeric(rat.split$stroma@meta.data$seurat_clusters) > 15)], features = canonical.markers, assay = "SCT", slot = "scale.data")
dev.off()

DefaultAssay(rat.split$stroma) <- "RNA"
rat.split$stroma <- NormalizeData(rat.split$stroma, assay = "RNA")
pdf("../../plots/whole-gland/stroma/markers1.pdf")
FeaturePlot(rat.split$stroma, features = c("ACTA2", "TAGLN", "VWF", "MKI67", "SPARCL1", "COL1A1"), label = TRUE)
dev.off()

rat.split$stroma$celltype.stroma <- dplyr::recode(rat.split$stroma$seurat_clusters, "0" = "fibroblasts", "1" = "fibroblasts", 
                                                  "2" = "endothelium", "3" = "endothelium",
                                                  "4" = "fibroblasts", "5" = "endothelium",
                                                  "6" = "fibroblasts", "7" = "fibroblasts", 
                                                  "8" = "endothelium", "9" = "endothelium", 
                                                  "10" = "muscle", "11" = "endothelium",
                                                  "12" = "fibroblasts", "13" = "immune", 
                                                  "14" = "epithelium", "15" = "stroma prolifrative", 
                                                  "16" = "fibroblasts", "17" = "endothelium",
                                                  "18" = "fibroblasts", "19" = "immune", 
                                                  "20" = "fibroblasts", "21" = "fibroblasts", 
                                                  "22" = "fibroblasts")
rat.split$stroma$celltype.stroma2 <- dplyr::recode(rat.split$stroma$seurat_clusters, "0" = "right-fibroblasts", "1" = "right-fibroblasts", 
                                                   "2" = "endothelium", "3" = "endothelium",
                                                   "4" = "right-fibroblasts", "5" = "endothelium",
                                                   "6" = "left-fibroblasts", "7" = "left-fibroblasts", 
                                                   "8" = "endothelium", "9" = "endothelium", 
                                                   "10" = "muscle", "11" = "endothelium",
                                                   "12" = "left-fibroblasts", "13" = "immune", 
                                                   "14" = "epithelium", "15" = "proliferative", 
                                                   "16" = "right-fibroblasts", "17" = "endothelium",
                                                   "18" = "right-fibroblasts", "19" = "immune", 
                                                   "20" = "right-fibroblasts", "21" = "left-fibroblasts", 
                                                   "22" = "left-fibroblasts")
pdf("../../plots/whole-gland/stroma/umap.celltype.stroma.pdf")
DimPlot(rat.split$stroma, group.by = "celltype.stroma")
dev.off()
pdf("../../plots/whole-gland/stroma/umap.celltype.stroma2.pdf")
DimPlot(rat.split$stroma, group.by = "celltype.stroma2")
dev.off()


DefaultAssay(rat.split$stroma) <- "RNA"
rat.split$stroma <- NormalizeData(rat.split$stroma, assay = "RNA")
markers.stroma <- c("SELE", "VWF", "LYVE1", "PDPN", "CDH5", "CD44", "VIM", "PDGFRA", "MYL1", "MKI67")
markers.stroma <- c("COL1A1", "COL3A1", "COL5A1", "LUM", "COL1A2", "COL6A1", "COL6A2", "FAP")
for (m in markers.stroma){
  pdf(paste0("../../plots/whole-gland/stroma/markers-", m, ".pdf"))
  print(FeaturePlot(rat.split$stroma, features = m))
  dev.off()
  pdf(paste0("../../plots/whole-gland/stroma/markers-", m, ".label.pdf"))
  print(FeaturePlot(rat.split$stroma, features = m, label = TRUE))
  dev.off()
}



############################################################################
############################################################################
###                                                                      ###
###                          STROMA DRUG EFFECT                          ###
###                                                                      ###
############################################################################
############################################################################
DefaultAssay(rat.split$stroma) <- "RNA"
rat.split$stroma <- NormalizeData(rat.split$stroma, assay = "RNA")
Idents(rat.split$stroma) <- rat.split$stroma$celltype.stroma2
n <- 0
p.val.thresh <- 0.1
stimulation.markers <- list()
for (c in c("left-fibroblasts", "right-fibroblasts", "endothelium")){
  Idents(rat.split$stroma) <- rat.split$stroma$celltype.stroma2
  rat.split$stroma@meta.data$clusterdrug <- paste(Idents(rat.split$stroma), rat.split$stroma$drug, sep = "_")
  Idents(rat.split$stroma) <- "clusterdrug"
  stimulation.markers[[n+1]] <- FindMarkers(rat.split$stroma, ident.1 = paste(c, "drug", sep = "_"), 
                                            ident.2 = paste(c, "control", sep = "_"), assay = "RNA")
  stimulation.markers[[n+1]] <- cbind("gene" = rownames(stimulation.markers[[n+1]]), stimulation.markers[[n+1]])
  stimulation.markers[[n+1]] <- as.data.frame(stimulation.markers[[n+1]])
  stimulation.markers[[n+1]] <- stimulation.markers[[n+1]][which(stimulation.markers[[n+1]]$p_val_adj < p.val.thresh),]
  stimulation.markers[[n+1]] <- stimulation.markers[[n+1]][order(stimulation.markers[[n+1]]$avg_log2FC, decreasing = TRUE),]
  names(stimulation.markers)[n+1] <- c
}




############################################################################
############################################################################
###                                                                      ###
###                                IMMUNE                                ###
###                                                                      ###
############################################################################
############################################################################
pdf("../../plots/whole-gland/immune/umap-clusters.pdf")
DimPlot(rat.split$immune, label = TRUE)
dev.off()
pdf("../../plots/whole-gland/immune/umap-drug.pdf")
DimPlot(rat.split$immune, group.by = "drug", cols = c("control" = "#C0C0C0", "drug" = "#FED242"))
dev.off()
pdf("../../plots/whole-gland/immune/umap-animal.pdf")
DimPlot(rat.split$immune, group.by = "sample")
dev.off()
pdf("../../plots/whole-gland/immune/umap-animal.pdf")
DimPlot(rat.split$immune, group.by = "sample")
dev.off()
pdf("../../plots/whole-gland/immune/umap-celltype-ex-classification.pdf")
DimPlot(rat.split$immune, group.by = "celltype.detailed")
dev.off()
pdf("../../plots/whole-gland/immune/ridge.qc.pdf")
RidgePlot(rat.split$immune, features = c("percent.mt", "nFeature_SCT"))
dev.off()



############################################################################
############################################################################
###                                                                      ###
###                         IMMUNE CELL SUBTYPES                         ###
###                                                                      ###
############################################################################
############################################################################
diff.expr.immune <- FindAllMarkers(rat.split$immune)
write.table(diff.expr.immune, quote = FALSE, row.names = FALSE, sep = "\t",
            file = "../../tables/gene-lists/whole.gland.diff.expr.immune.txt")

DefaultAssay(rat.split$immune) <- "SCT"
pdf("../../plots/whole-gland/immune/heatmap.pdf", height = 10, width = 20)
DoHeatmap(rat.split$immune, features = canonical.markers)
dev.off()
pdf("../../plots/whole-gland/immune/heatmap-zoom.pdf", height = 10, width = 20)
DoHeatmap(rat.split$immune[,which(as.numeric(rat.split$immune@meta.data$seurat_clusters) > 15)], features = canonical.markers, assay = "SCT", slot = "scale.data")
dev.off()

DefaultAssay(rat.split$immune) <- "RNA"
rat.split$immune <- NormalizeData(rat.split$immune, assay = "RNA")
pdf("../../plots/whole-gland/immune/markers1.pdf")
FeaturePlot(rat.split$immune, features = c("PTPRC", "CD8A", "MKI67", "MS4A1"), label = TRUE)
dev.off()

rat.split$immune$celltype.immune <- dplyr::recode(rat.split$immune$seurat_clusters, "0" = "T cell", "1" = "B cell", 
                                                  "2" = "B cell", "3" = "B cell",
                                                  "4" = "T cell", "5" = "B cell",
                                                  "6" = "B cell", "7" = "macrophage", 
                                                  "8" = "NK", "9" = "epithelium", 
                                                  "10" = "B cell", "11" = "T cell",
                                                  "12" = "macrophage", "13" = "T cell", 
                                                  "14" = "macrophage", "15" = "macrophage", 
                                                  "16" = "B cell", "17" = "macrophage",
                                                  "18" = "proliferative", "19" = "B cell", 
                                                  "20" = "B cell", "21" = "macrophage", 
                                                  "22" = "dendritic cell", "23" = "T cell", "24" = "T cell")
rat.split$immune$celltype.immune2 <- dplyr::recode(rat.split$immune$seurat_clusters, "0" = "T cell", "1" = "B cell", 
                                                   "2" = "B cell", "3" = "B cell",
                                                   "4" = "T cell", "5" = "B cell",
                                                   "6" = "B cell", "7" = "macrophage-affected", 
                                                   "8" = "NK", "9" = "epithelium", 
                                                   "10" = "B cell", "11" = "T cell",
                                                   "12" = "macrophage-affected", "13" = "T cell", 
                                                   "14" = "macrophage", "15" = "macrophage", 
                                                   "16" = "B cell", "17" = "macrophage",
                                                   "18" = "proliferative", "19" = "B cell", 
                                                   "20" = "B cell", "21" = "macrophage", 
                                                   "22" = "dendritic cell", "23" = "T cell", "24" = "T cell")

pdf("../../plots/whole-gland/immune/umap.celltype.immune.pdf")
DimPlot(rat.split$immune, group.by = "celltype.immune")
dev.off()
pdf("../../plots/whole-gland/immune/umap.celltype.immune2.pdf")
DimPlot(rat.split$immune, group.by = "celltype.immune2")
dev.off()


DefaultAssay(rat.split$immune) <- "RNA"
rat.split$immune <- NormalizeData(rat.split$immune, assay = "RNA")
markers.immune <- c("PTPRC", "CD3D", "CD4", "CD8A", "CD19", "MS4A1", "CSF1R", "NKG7", "CLEC9A", "MKI67")
for (m in markers.immune){
  pdf(paste0("../../plots/whole-gland/immune/markers-", m, ".pdf"))
  print(FeaturePlot(rat.split$immune, features = m))
  dev.off()
  pdf(paste0("../../plots/whole-gland/immune/markers-", m, ".label.pdf"))
  print(FeaturePlot(rat.split$immune, features = m, label = TRUE))
  dev.off()
}



############################################################################
############################################################################
###                                                                      ###
###                          IMMUNE DRUG EFFECT                          ###
###                                                                      ###
############################################################################
############################################################################
DefaultAssay(rat.split$immune) <- "RNA"
rat.split$immune <- NormalizeData(rat.split$immune, assay = "RNA")
Idents(rat.split$immune) <- rat.split$stroma$celltype.immune2
n <- 0
p.val.thresh <- 0.1
stimulation.markers <- list()
for (c in c("macrophage-affected")){
  Idents(rat.split$immune) <- rat.split$immune$celltype.immune2
  rat.split$immune@meta.data$clusterdrug <- paste(Idents(rat.split$immune), rat.split$immune$drug, sep = "_")
  Idents(rat.split$immune) <- "clusterdrug"
  stimulation.markers[[n+1]] <- FindMarkers(rat.split$immune, ident.1 = paste(c, "drug", sep = "_"), 
                                            ident.2 = paste(c, "control", sep = "_"), assay = "RNA")
  stimulation.markers[[n+1]] <- cbind("gene" = rownames(stimulation.markers[[n+1]]), stimulation.markers[[n+1]])
  stimulation.markers[[n+1]] <- as.data.frame(stimulation.markers[[n+1]])
  stimulation.markers[[n+1]] <- stimulation.markers[[n+1]][which(stimulation.markers[[n+1]]$p_val_adj < p.val.thresh),]
  stimulation.markers[[n+1]] <- stimulation.markers[[n+1]][order(stimulation.markers[[n+1]]$avg_log2FC, decreasing = TRUE),]
  names(stimulation.markers)[n+1] <- c
}



############################################################################
############################################################################
###                                                                      ###
###                              EPITHELIUM                              ###
###                                                                      ###
############################################################################
############################################################################
pdf("../../plots/whole-gland/epithelium/umap-clusters.pdf")
DimPlot(rat.split$epithelial, label = TRUE)
dev.off()
pdf("../../plots/whole-gland/epithelium/umap-drug.pdf")
DimPlot(rat.split$epithelial, group.by = "drug", cols = c("control" = "#C0C0C0", "drug" = "#FED242"))
dev.off()
pdf("../../plots/whole-gland/epithelium/umap-animal.pdf")
DimPlot(rat.split$epithelial, group.by = "sample")
dev.off()
pdf("../../plots/whole-gland/epithelium/umap-animal.pdf")
DimPlot(rat.split$epithelial, group.by = "sample")
dev.off()
pdf("../../plots/whole-gland/epithelium/umap-celltype-ex-classification.pdf")
DimPlot(rat.split$epithelial, group.by = "celltype.detailed")
dev.off()
pdf("../../plots/whole-gland/epithelium/ridge.qc.pdf")
RidgePlot(rat.split$epithelial, features = c("percent.mt", "nFeature_SCT"))
dev.off()



############################################################################
############################################################################
###                                                                      ###
###                       EPITHELIUM CELL SUBTYPES                       ###
###                                                                      ###
############################################################################
############################################################################
diff.expr.epithelial <- FindAllMarkers(rat.split$epithelial)
write.table(diff.expr.epithelial, quote = FALSE, row.names = FALSE, sep = "\t",
            file = "../../tables/gene-lists/whole.gland.diff.expr.epithelium.txt")

DefaultAssay(rat.split$epithelial) <- "SCT"
pdf("../../plots/whole-gland/epithelium/heatmap.pdf", height = 10, width = 20)
DoHeatmap(rat.split$epithelial, features = canonical.markers)
dev.off()
pdf("../../plots/whole-gland/epithelium/heatmap-zoom.pdf", height = 10, width = 20)
DoHeatmap(rat.split$epithelial[,which(as.numeric(rat.split$epithelial@meta.data$seurat_clusters) > 15)], features = canonical.markers, assay = "SCT", slot = "scale.data")
dev.off()

DefaultAssay(rat.split$epithelial) <- "RNA"
rat.split$epithelial <- NormalizeData(rat.split$epithelial, assay = "RNA")
pdf("../../plots/whole-gland/epithelium/markers1.pdf")
FeaturePlot(rat.split$epithelial, features = c("LALBA", "PGR", "KIT", "MKI67", "ACTA2", "EPCAM"), label = TRUE)
dev.off()


rat.split$epithelial$celltype.epithelial <- dplyr::recode(rat.split$epithelial$seurat_clusters, "0" = "basal", 
                                                          "1" = "basal", 
                                                          "2" = "luminal progenitor", "3" = "mature luminal",
                                                          "4" = "luminal progenitor", "5" = "luminal progenitor",
                                                          "6" = "luminal progenitor", "7" = "mature luminal", 
                                                          "8" = "luminal progenitor", "9" = "immune", 
                                                          "10" = "luminal progenitor", "11" = "mature luminal",
                                                          "12" = "luminal progenitor", "13" = "luminal progenitor", 
                                                          "14" = "basal", "15" = "basal", 
                                                          "16" = "neural cells", "17" = "immune",
                                                          "18" = "basal", "19" = "basal")
pdf("../../plots/whole-gland/epithelium/umap.celltype.epithelial.pdf")
DimPlot(rat.split$epithelial, group.by = "celltype.epithelial")
dev.off()

rat.split$epithelial$celltype.epithelial2 <- dplyr::recode(rat.split$epithelial$seurat_clusters, "0" = "basal", 
                                                           "1" = "basal", 
                                                           "2" = "luminal progenitor", "3" = "mature luminal",
                                                           "4" = "luminal progenitor", "5" = "luminal progenitor",
                                                           "6" = "luminal progenitor", "7" = "mature luminal", 
                                                           "8" = "luminal progenitor", "9" = "immune", 
                                                           "10" = "luminal progenitor", "11" = "mature luminal",
                                                           "12" = "SBCs", "13" = "luminal progenitor", 
                                                           "14" = "basal", "15" = "basal", 
                                                           "16" = "neural cells", "17" = "immune",
                                                           "18" = "basal", "19" = "SBCs")

pdf("../../plots/whole-gland/epithelium/umap.celltype.epithelial2.pdf")
DimPlot(rat.split$epithelial, group.by = "celltype.epithelial2")
dev.off()

DefaultAssay(rat.split$epithelial) <- "RNA"
rat.split$epithelial <- NormalizeData(rat.split$epithelial, assay = "RNA")
markers.epithelial <- c("KRT5", "PRLR", "KRT18", "MKI67", "ESR1", "TP63", "LALBA", "KIT", "ACTA2", "SPARCL1", "SORBS2", 
                        "ZEB2", "LUM", "LRG1", "NRG1", "CDH5", "PROCR")
for (m in markers.epithelial){
  pdf(paste0("../../plots/whole-gland/epithelium/markers-", m, ".pdf"))
  print(FeaturePlot(rat.split$epithelial, features = m))
  dev.off()
  pdf(paste0("../../plots/whole-gland/epithelium/markers-", m, ".label.pdf"))
  print(FeaturePlot(rat.split$epithelial, features = m, label = TRUE))
  dev.off()
}



############################################################################
############################################################################
###                                                                      ###
###                        EPITHELIUM DRUG EFFECT                        ###
###                                                                      ###
############################################################################
############################################################################
DefaultAssay(rat.split$epithelial) <- "RNA"
rat.split$epithelial <- NormalizeData(rat.split$epithelial, assay = "RNA")
Idents(rat.split$epithelial) <- rat.split$epithelial$celltype.epithelial
n <- 0
p.val.thresh <- 0.1
stimulation.markers <- list()
for (c in c("basal", "luminal progenitor", "mature luminal")){
  Idents(rat.split$epithelial) <- rat.split$epithelial$celltype.epithelial
  rat.split$epithelial@meta.data$clusterdrug <- paste(Idents(rat.split$epithelial), rat.split$epithelial$drug, sep = "_")
  Idents(rat.split$epithelial) <- "clusterdrug"
  stimulation.markers[[n+1]] <- FindMarkers(rat.split$epithelial, ident.1 = paste(c, "drug", sep = "_"), 
                                            ident.2 = paste(c, "control", sep = "_"), assay = "RNA")
  stimulation.markers[[n+1]] <- cbind("gene" = rownames(stimulation.markers[[n+1]]), stimulation.markers[[n+1]])
  stimulation.markers[[n+1]] <- as.data.frame(stimulation.markers[[n+1]])
  stimulation.markers[[n+1]] <- stimulation.markers[[n+1]][which(stimulation.markers[[n+1]]$p_val_adj < p.val.thresh),]
  stimulation.markers[[n+1]] <- stimulation.markers[[n+1]][order(stimulation.markers[[n+1]]$avg_log2FC, decreasing = TRUE),]
  names(stimulation.markers)[n+1] <- c
}
