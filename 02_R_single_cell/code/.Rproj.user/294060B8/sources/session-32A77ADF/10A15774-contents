## Single cell analysis on whole mammary gland SD data, including all cell types
## Author: Simona Cristea
## Date: October 2022
## Publication: Breast cancer prevention by short-term inhibition of TGFβ signaling, Nat Comm 2022


## load libraries
source(here::here("libraries.R"))
sample <- "new"



############################################################################
############################################################################
###                                                                      ###
###                 SECTION 1: READ IN SEQUENCING OUTPUT                 ###
###                                                                      ###
############################################################################
############################################################################
dataset_loc <- "../../data/cellranger-count/"
samples <- c("veh1", "veh2", "veh4", "ly1", "ly2", "ly4", "ly6", "ly11", "ly14")
rat.data <- sapply(samples, function(i){
  d10x <- Read10X(file.path(dataset_loc,i,"outs/filtered_feature_bc_matrix/"), )
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})
rat.data <- do.call("cbind", rat.data)
rownames(rat.data) <- toupper(rownames(rat.data))
rat <- CreateSeuratObject(
  rat.data,
  project = "all_sd4_samples", 
  min.cells = 10,
  min.features = 100,
  names.field = 2,
  names.delim = "\\-")
rat <- AddMetaData(object = rat, metadata = sapply(rownames(rat@meta.data), function(x){paste(strsplit(x, "-")[[1]][2], collapse= "-")}), col.name = "sample")
table(rat$sample)
rat <- AddMetaData(object = rat, metadata = recode(sapply(rat@meta.data$sample, function(y){sub("^([[:alpha:]]*).*", "\\1", y)}),
                                                   "ly" = "drug", "veh" = "control"), col.name = "drug")
Idents(rat) <- rat@meta.data$sample



############################################################################
############################################################################
###                                                                      ###
###                              CLEAN DATA                              ###
###                                                                      ###
############################################################################
############################################################################
mito.genes <- c("ND1", "ND2", "ND3", "ND4L", "ND5", "ND4", "ND6", "CYTB", "COX1", "COX2", "COX3", "ATP6", "ATP8")
rat[["percent.mt"]] <- PercentageFeatureSet(rat, features =  mito.genes)


## remove all cells that have mitochondria percentage higher than a given threshold
percent.mito.thresh <- 15
rat <- rat[,-which(rat$percent.mt > percent.mito.thresh)]



###########################################################################
###########################################################################
###                                                                     ###
###                           PREPROCESS DATA                           ###
###                                                                     ###
###########################################################################
###########################################################################
rat <- SCTransform(rat, verbose = TRUE)
rat <- RunPCA(rat, features = VariableFeatures(object = rat, assay = "SCT"), assay = "SCT")
rat <- FindNeighbors(rat, dims = 1:30, verbose = FALSE, assay = "SCT", reduction = "pca", graph.name = "SNN.SCT")
rat <- FindClusters(rat, verbose = TRUE, graph.name = "SNN.SCT")
rat <- RunUMAP(rat, dims = 1:30, verbose = TRUE, reduction.name = "UMAP.SCT", assay = "SCT")


pdf("../../plots/whole-gland/umap.clusters.pdf")
DimPlot(rat, label = TRUE)
dev.off()

pdf("../../plots/whole-gland/umap.drug.pdf")
DimPlot(rat, group.by = "drug")
dev.off()

pdf("../../plots/whole-gland/ridge.qc.pdf")
RidgePlot(rat, features = c("percent.mt", "nFeature_SCT"))
dev.off()

canonical.markers <- read.csv(file = "../../tables/gene-lists/canonical_markers.csv")
canonical.markers <- canonical.markers$marker
DefaultAssay(rat) <- "SCT"
pdf("../../plots/whole-gland/heatmap.pdf", height = 10, width = 20)
DoHeatmap(rat, features = canonical.markers)
dev.off()
pdf("../../plots/whole-gland/heatmap-zoom.pdf", height = 10, width = 20)
DoHeatmap(rat[,which(as.numeric(rat@meta.data$seurat_clusters) > 15)], features = canonical.markers, assay = "SCT", slot = "scale.data")
dev.off()

DefaultAssay(rat) <- "RNA"
rat <- NormalizeData(rat, assay = "RNA")
pdf("../../plots/whole-gland/markers.pdf", width = 15)
FeaturePlot(rat, features = c("PTPRC", "EPCAM", "LALBA", "ACTA2", "PGR", "KIT", "SPARCL1",
                              "KRT7", "KRT14", "ZEB2"))
dev.off()

pdf("../../plots/whole-gland/umap.sample.pdf")
DimPlot(rat, group.by = "sample")
dev.off()

pdf("../../plots/whole-gland/umap.sample.fig.pdf")
DimPlot(rat, group.by = "sample",
        cols = c("ly1" = "darkorange1", "ly11" = "maroon2", "ly14" = "firebrick1", "ly2" = "firebrick4",
                 "ly4" = "tomato1", "ly6" = "sienna3",
                 "veh1" = "steelblue3", "veh2" = "turquoise2", "veh4" = "mediumpurple3"))
dev.off()

### differential expression for all
diff.expr <- FindAllMarkers(rat)



###########################################################################
###########################################################################
###                                                                     ###
###                   ASSIGN CELL TYPES USING MARKERS                   ###
###                                                                     ###
###########################################################################
###########################################################################
rat$celltype <- dplyr::recode(rat$seurat_clusters, "0" = "epithelial", "1" = "epithelial", "2" = "immune", "3" = "immune", 
       "4" = "stroma", "5" = "immune",
       "6" = "immune", "7" = "immune", "8" = "stroma", "9" = "stroma", 
       "10" = "epithelial", "11" = "immune",
       "12" = "stroma", "13" = "epithelial", 
       "14" = "stroma", "15" = "stroma", "16" = "immune", "17" = "immune",
       "18" = "epithelial", "19" = "immune", "20" = "immune", 
       "21" = "nerve-myelin", "22" = "immune",
       "23" = "stroma", "24" = "immune", "25" = "vein")
pdf("../../plots/whole-gland/umap.celltype.pdf")
DimPlot(rat, group.by = "celltype")
dev.off()

rat$celltype.detailed <- dplyr::recode(rat$seurat_clusters, "0" = "basal", "1" = "LP", "2" = "T", "3" = "B", 
       "4" = "vascular smooth muscle cells", "5" = "B",
       "6" = "B", "7" = "B", "8" = "fibroblasts", "9" = "endothelium TIE1+", 
       "10" = "alveolar", "11" = "NK or T",
       "12" = "endothelium", "13" = "alveolar", 
       "14" = "stroma-progenitor", "15" = "FAP-fibroblasts", "16" = "macrophage", "17" = "macrophage",
       "18" = "LP", "19" = "NK", "20" = "monocyte-macrophage", "21" = "nerve-myelin", "22" = "monocyte-macrophage",
       "23" = "endothelium TIE1+", "24" = "T", "25" = "endothelium TIE1+")
pdf("../../plots/whole-gland/umap.celltype.detailed.pdf")
DimPlot(rat, group.by = "celltype.detailed")
dev.off()



###########################################################################
###########################################################################
###                                                                     ###
###                             DRUG EFFECT                             ###
###                                                                     ###
###########################################################################
###########################################################################
DefaultAssay(rat) <- "RNA"
rat <- NormalizeData(rat, assay = "RNA")
id <- Idents(rat)
no.clusters <- length(unique(Idents(rat)))
p.val.thresh <- 0.1
stimulation.markers <- list()
for (c in 0:(no.clusters-1)){
  Idents(rat) <- id
  rat@meta.data$clusterdrug <- paste(Idents(rat), rat$drug, sep = "_")
  Idents(rat) <- "clusterdrug"
  stimulation.markers[[c+1]] <- FindMarkers(rat, ident.1 = paste(c, "drug", sep = "_"), 
                                            ident.2 = paste(c, "control", sep = "_"), assay = "RNA")
  stimulation.markers[[c+1]] <- cbind("gene" = rownames(stimulation.markers[[c+1]]), stimulation.markers[[c+1]])
  stimulation.markers[[c+1]] <- as.data.frame(stimulation.markers[[c+1]])
  stimulation.markers[[c+1]] <- stimulation.markers[[c+1]][which(stimulation.markers[[c+1]]$p_val_adj < p.val.thresh),]
  stimulation.markers[[c+1]] <- stimulation.markers[[c+1]][order(stimulation.markers[[c+1]]$avg_log2FC, decreasing = TRUE),]
  names(stimulation.markers)[c+1] <- paste0("cluster", c)
}

pdf("../../plots/whole-gland/umap.drug.pdf")
DimPlot(rat, group.by = "drug", cols = c("control" = "#C0C0C0", "drug" = "#FED242"))
dev.off()

## only keep main cell types
rat$maincelltype <- 1
rat$maincelltype[which(Idents(rat) %in% c("nerve-myelin", "vein"))] <- 0
rat.main <- SplitObject(rat, split.by = "maincelltype")
rat.main <- rat.main$`1`



############################################################################
############################################################################
###                                                                      ###
###                         PREVENTION SIGNATURE                         ###
###                                                                      ###
############################################################################
############################################################################
genes <- read.table("../../tables/gene-lists/Positive genes in both strains.txt", sep = "\t", header = TRUE)
genes <- genes$gene
rat$SBCsign <- colSums(rat@assays$RNA@data[na.omit(match(genes, rownames(rat@assays$RNA@data))),])
pdf("../../plots/whole-gland/prevention-sign.pdf")
FeaturePlot(rat, features = c("SBCsign"))
dev.off()


rat <- saveRDS(file = "../../data/rdata/rat-whole-gland.RDS")
