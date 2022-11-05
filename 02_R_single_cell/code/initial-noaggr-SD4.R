## Single cell analysis on sorted epithelium SD data
## Author: Simona Cristea
## Date: October 2022
## Publication: Breast cancer prevention by short-term inhibition of TGFβ signaling, Nat Comm 2022


## load libraries
source(here::here("libraries.R"))
sample <- "sd"



############################################################################
############################################################################
###                                                                      ###
###                 SECTION 1: READ IN SEQUENCING OUTPUT                 ###
###                                                                      ###
############################################################################
############################################################################

# replace with your own path where cellranger count outputs are
dataset_loc <- "../../data/cellranger-count/"
samples <- c("SD4-Bas10", "SD4-Bas11", "SD4-Bas12", "SD4-Bas13", "SD4-Bas14", "SD4-Bas1", "SD4-Bas2", 
             "SD4-Bas3", "SD4-Bas9", "SD4-Lum10", "SD4-Lum11", "SD4-Lum12", "SD4-Lum13", "SD4-Lum14", 
             "SD4-Lum1", "SD4-Lum2", "SD4-Lum3", "SD4-Lum9")

rat.data <- sapply(samples, function(i){
  d10x <- Read10X(file.path(dataset_loc,i,"outs/"), )
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
rat <- AddMetaData(object = rat, metadata = sapply(rownames(rat@meta.data), function(x){paste(strsplit(x, "-")[[1]][2:3], collapse= "-")}), col.name = "sample")
rat <- AddMetaData(object = rat, metadata = recode(tolower(sapply(sapply(rat@meta.data$sample, function(x){paste(strsplit(x, "-")[[1]][2], collapse= "-")}), 
                                                   function(y){return(gsub("[[:digit:]]","",y))})), "bas" = "basal", "lum" = "luminal"),
                                                   col.name = "celltype")
rat <- AddMetaData(object = rat, metadata = paste0("animal", sapply(sapply(rat@meta.data$sample, function(x){paste(strsplit(x, "-")[[1]][2], collapse= "-")}), 
                                                   function(y){return(gsub("[^[:digit:]]","",y))})), col.name = "animal")
rat@meta.data$drug <- sapply(rat@meta.data$sample, function(x){ 
  if (x == "SD4-Bas10")
    return("drug")
  if (x == "SD4-Bas11")
    return("drug")
  if (x == "SD4-Bas12")
    return("drug")
  if (x == "SD4-Bas13")
    return("drug")
  if (x == "SD4-Bas14")
    return("drug")
  if (x == "SD4-Bas1")
    return("control")
  if (x == "SD4-Bas2")
    return("control")
  if (x == "SD4-Bas3")
    return("control")
  if (x == "SD4-Bas9")
    return("drug")
  if (x == "SD4-Lum10")
    return("drug")
  if (x == "SD4-Lum11")
    return("drug")
  if (x == "SD4-Lum12")
    return("drug")
  if (x == "SD4-Lum13")
    return("drug")
  if (x == "SD4-Lum14")
    return("drug")
  if (x == "SD4-Lum1")
    return("control")
  if (x == "SD4-Lum2")
    return("control")
  if (x == "SD4-Lum3")
    return("control")
  if (x == "SD4-Lum9")
    return("drug")
})
Idents(rat) <- rat@meta.data$sample



############################################################################
############################################################################
###                                                                      ###
###                              CLEAN DATA                              ###
###                                                                      ###
############################################################################
############################################################################

## calculate proportion of mitochondrial genes with Seurat
mito.genes <- c("ND1", "ND2", "ND3", "ND4L", "ND5", "ND4", "ND6", "CYTB", "COX1", "COX2", "COX3", "ATP6", "ATP8")
rat[["percent.mt"]] <- PercentageFeatureSet(rat, features =  mito.genes)


## remove all cells that have mitochondria percentage higher than a given threshold
percent.mito.thresh <- 15
rat <- rat[,-which(rat$percent.mt > percent.mito.thresh)]


## save data with only these QC steps applied
saveRDS(rat, file = "../../data/rdata/rat-initial-sd-only-mito.RDS")


## clean data of contamination before normalization (PTPRC+ and KRT- cells)
immune <- which(as.vector(rat@assays$RNA["PTPRC",]!=0))
length(immune)
keratins <- grep("KRT", rownames(rat), value = TRUE)
keratins <- setdiff(keratins, c(grep("KRTAP", keratins, value = TRUE), grep("KRTCAP", keratins, value = TRUE)))
keratins.vals <- colSums(as.matrix(rat@assays$RNA[keratins,]))
no.krt <- which(keratins.vals < 2)
length(no.krt)
to.exclude <- union(immune, no.krt)
length(to.exclude)
rat <- rat[,-to.exclude]
dim(rat)[2]



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



###########################################################################
###########################################################################
###                                                                     ###
###                  DISTRIBUTION OF CELLS PER CLUSTER                  ###
###                                                                     ###
###########################################################################
###########################################################################
distr.abs <- list()
distr.rel <- list()
no.clusters.here <- length(unique(rat@meta.data$seurat_clusters))
for (c in 0:(no.clusters.here-1)){
  cells.idx <- which(rat@meta.data$seurat_clusters == c)
  abs.celltype <- table(rat@meta.data$celltype[cells.idx])
  if (length(abs.celltype) == 1){
    abs.celltype <- c(abs.celltype, 0)
    names(abs.celltype)[2] <- setdiff(c("luminal", "basal"), names(abs.celltype))
    abs.celltype <- abs.celltype[order(names(abs.celltype))]
  }
  abs.celltype <- abs.celltype[order(names(abs.celltype))]
  prop.celltype <- abs.celltype/length(cells.idx)
  abs.animals <- table(rat@meta.data$animal[cells.idx])
  if (length(abs.animals) != length(unique(rat@meta.data$animal))){
    n <- names(abs.animals)
    missing.animals <- setdiff(unique(rat@meta.data$animal), names(abs.animals))
    abs.animals <- c(rep(0, length(missing.animals)), abs.animals)
    names(abs.animals) <- c(missing.animals, n)
  }
  abs.animals <- abs.animals[order(names(abs.animals))]
  prop.animals <- abs.animals/length(cells.idx)
  distr.abs[[c+1]] <- c("total" = length(cells.idx), abs.celltype, abs.animals)
  names(distr.abs)[c+1] <- paste0("cluster", c)
  distr.rel[[c+1]] <- c(prop.celltype, prop.animals)
  names(distr.rel)[c+1] <- paste0("cluster", c)
}
distr.abs <- do.call(rbind, distr.abs)
distr.rel <- do.call(rbind, distr.rel)
distr.rel <- round(distr.rel, 3)
distr.abs <- cbind("cluster" = rownames(distr.abs), distr.abs)
distr.abs <- as.data.frame(distr.abs)
distr.rel <- cbind("cluster" = rownames(distr.rel), distr.rel)
distr.rel <- as.data.frame(distr.rel)
write.table(distr.abs, sep= "\t", quote = FALSE, row.names = FALSE, 
            file = paste0(here::here("../../tables/cells-per-cluster/"), sample, ".cells.per.cluster.abs.v1.txt"))
write.table(distr.rel, sep= "\t", quote = FALSE, row.names = FALSE, 
            file = paste0(here::here("../../tables/cells-per-cluster/"), sample, ".cells.per.cluster.rel.v1.txt"))



############################################################################
############################################################################
###                                                                      ###
###                   REMOVE CELLS LESS THAN THRESHOLD                   ###
###                                                                      ###
############################################################################
############################################################################
thresh.cont <- 0.05
cells.remove <- list()
for (i in 1:nrow(distr.rel)){
  bas <- as.numeric(as.character(distr.rel$basal[i]))
  lum <- as.numeric(as.character(distr.rel$luminal[i]))
  if (bas <= thresh.cont){
    # remove the basal cells from a luminal cluster
    cells.remove[[i]] <- intersect(colnames(rat)[which(rat@meta.data$seurat_clusters == (i-1))], 
                                   colnames(rat)[which(rat@meta.data$celltype == "basal")])
  } else if (lum <= thresh.cont){
    # remove the luminal cells from a basal cluster
    cells.remove[[i]] <- intersect(colnames(rat)[which(rat@meta.data$seurat_clusters == (i-1))], 
                                   colnames(rat)[which(rat@meta.data$celltype == "luminal")])
  } else
    cells.remove[[i]] <- "no cell"
}
names(cells.remove) <- distr.rel$cluster
cells.remove <- setdiff(unlist(cells.remove), "no cell")
table(rat@meta.data$celltype[match(cells.remove, colnames(rat))])

rat.clean <- rat[,-match(cells.remove, colnames(rat))]
rat.full <- rat
rat <- rat.clean
rm(rat.clean)
dim(rat)[2]
table(rat$celltype)


###########################################################################
###########################################################################
###                                                                     ###
###                  REDO PREPROCESSING AFTER CLEANING                  ###
###                                                                     ###
###########################################################################
###########################################################################
rat <- SCTransform(rat, verbose = TRUE)
pdf(paste0(here::here("../../plots/preprocessing/"), sample, "-vargenes-sctransform.pdf"))
LabelPoints(VariableFeaturePlot(rat, assay = "SCT"), points = head(VariableFeatures(rat, assay = "SCT"), 10), repel = TRUE)
dev.off()

rat <- RunPCA(rat, features = VariableFeatures(object = rat, assay = "SCT"), assay = "SCT")
pdf(paste0(here::here("../../plots/preprocessing/"), sample, "-pca-sample-SCT.pdf"))
DimPlot(rat, reduction = "pca", group.by = "sample", label = FALSE)
dev.off()
pdf(paste0(here::here("../../plots/preprocessing/"), sample, "-pca-animal-SCT.pdf"))
DimPlot(rat, reduction = "pca", group.by = "animal", label = FALSE)
dev.off()
pdf(paste0(here::here("../../plots/preprocessing/"), sample, "-pca-drug-SCT.pdf"))
DimPlot(rat, reduction = "pca", group.by = "drug", label = FALSE)
dev.off()
pdf(paste0(here::here("../../plots/preprocessing/"), sample, "-pca-celltype-SCT.pdf"))
DimPlot(rat, reduction = "pca", group.by = "celltype", label = FALSE)
dev.off()

rat <- FindNeighbors(rat, dims = 1:30, verbose = FALSE, assay = "SCT", reduction = "pca", graph.name = "SNN.SCT")
rat <- FindClusters(rat, verbose = TRUE, graph.name = "SNN.SCT")
rat <- RunUMAP(rat, dims = 1:30, verbose = TRUE, reduction.name = "UMAP.SCT", assay = "SCT")

pdf(paste0(here::here("../../plots/aci-sorted/"), sample, "-umap.clusters-SCT.pdf"), width = 8)
DimPlot(rat, label = TRUE)
dev.off()
pdf(paste0(here::here("../../plots/aci-sorted/"), sample, "-umap.sample-SCT.pdf"), width = 8)
DimPlot(rat, label = TRUE, group.by = "sample")
dev.off()
pdf(paste0(here::here("../../plots/aci-sorted/"), sample, "-umap.animal-SCT.pdf"), width = 8)
DimPlot(rat, label = TRUE, group.by = "animal")
dev.off()
pdf(paste0(here::here("../../plots/aci-sorted/"), sample, "-umap.drug-SCT.pdf"), width = 8)
DimPlot(rat, label = TRUE, group.by = "drug")
dev.off()
pdf(paste0(here::here("../../plots/aci-sorted/"), sample, "-umap.celltype-SCT.pdf"), width = 8)
DimPlot(rat, label = TRUE, group.by = "celltype")
dev.off()
pdf(paste0(here::here("../../plots/aci-sorted/"), sample, "-percent.mt.ridge.pdf"), width = 8)
RidgePlot(rat, features = "percent.mt")
dev.off()



###########################################################################
###########################################################################
###                                                                     ###
###              DISTRIBUTION OF CELLS PER CLUSTER: REPEAT              ###
###                                                                     ###
###########################################################################
###########################################################################
distr.abs <- list()
distr.rel <- list()
no.clusters.here <- length(unique(rat@meta.data$seurat_clusters))
for (c in 0:(no.clusters.here-1)){
  cells.idx <- which(rat@meta.data$seurat_clusters == c)
  abs.celltype <- table(rat@meta.data$celltype[cells.idx])
  if (length(abs.celltype) == 1){
    abs.celltype <- c(abs.celltype, 0)
    names(abs.celltype)[2] <- setdiff(c("luminal", "basal"), names(abs.celltype))
    abs.celltype <- abs.celltype[order(names(abs.celltype))]
  }
  abs.celltype <- abs.celltype[order(names(abs.celltype))]
  prop.celltype <- abs.celltype/length(cells.idx)
  abs.animals <- table(rat@meta.data$animal[cells.idx])
  if (length(abs.animals) != length(unique(rat@meta.data$animal))){
    n <- names(abs.animals)
    missing.animals <- setdiff(unique(rat@meta.data$animal), names(abs.animals))
    abs.animals <- c(rep(0, length(missing.animals)), abs.animals)
    names(abs.animals) <- c(missing.animals, n)
  }
  abs.animals <- abs.animals[order(names(abs.animals))]
  prop.animals <- abs.animals/length(cells.idx)
  distr.abs[[c+1]] <- c("total" = length(cells.idx), abs.celltype, abs.animals)
  names(distr.abs)[c+1] <- paste0("cluster", c)
  distr.rel[[c+1]] <- c(prop.celltype, prop.animals)
  names(distr.rel)[c+1] <- paste0("cluster", c)
}
distr.abs <- do.call(rbind, distr.abs)
distr.rel <- do.call(rbind, distr.rel)
distr.rel <- round(distr.rel, 3)
distr.abs <- cbind("cluster" = rownames(distr.abs), distr.abs)
distr.abs <- as.data.frame(distr.abs)
distr.rel <- cbind("cluster" = rownames(distr.rel), distr.rel)
distr.rel <- as.data.frame(distr.rel)
write.table(distr.abs, sep= "\t", quote = FALSE, row.names = FALSE, 
            file = paste0(here::here("../../tables/cells-per-cluster/"), sample, ".cells.per.cluster.abs-v2.txt"))
write.table(distr.rel, sep= "\t", quote = FALSE, row.names = FALSE, 
            file = paste0(here::here("../../tables/cells-per-cluster/"), sample, ".cells.per.cluster.rel-v2.txt"))



############################################################################
############################################################################
###                                                                      ###
###               REMOVE CELLS LESS THAN THRESHOLD: REPEAT               ###
###                                                                      ###
############################################################################
############################################################################
thresh.cont <- 0.05
cells.remove <- list()
for (i in 1:nrow(distr.rel)){
  bas <- as.numeric(as.character(distr.rel$basal[i]))
  lum <- as.numeric(as.character(distr.rel$luminal[i]))
  if (bas <= thresh.cont){
    # remove the basal cells from a luminal cluster
    cells.remove[[i]] <- intersect(colnames(rat)[which(rat@meta.data$seurat_clusters == (i-1))], 
                                   colnames(rat)[which(rat@meta.data$celltype == "basal")])
  } else if (lum <= thresh.cont){
    # remove the luminal cells from a basal cluster
    cells.remove[[i]] <- intersect(colnames(rat)[which(rat@meta.data$seurat_clusters == (i-1))], 
                                   colnames(rat)[which(rat@meta.data$celltype == "luminal")])
  } else
    cells.remove[[i]] <- "no cell"
}
names(cells.remove) <- distr.rel$cluster
cells.remove <- setdiff(unlist(cells.remove), "no cell")
table(rat@meta.data$celltype[match(cells.remove, colnames(rat))])

## remove if anything left to remove
# rat <- rat[,-match(cells.remove, colnames(rat))]
# dim(rat)
# table(rat$celltype)



###########################################################################
###########################################################################
###                                                                     ###
###                          PLOTS FOR FIGURES                          ###
###                                                                     ###
###########################################################################
###########################################################################
markers <- list()
markers$basal <- c("TP63", "EGFR", "KRT14", "KRT5", "KIT", "KRT17", "NGFR", "VIM")
markers$luminal <- c("KRT8", "KRT18", "ESR1", "GATA3", "KRT19", "MUC1", "PGR", "FOXA1", "SPDEF", "XBP1")
markers$gate <- c("ITGB1", "CD24")

pdf(paste0(here::here("../../plots/sd-sorted/"), sample, "-markers.umap.gate-SCT.pdf"), width = 14)
FeaturePlot(rat, features = markers$gate, pt.size = 0.2)
dev.off()

pdf(paste0(here::here("../../plots/sd-sorted/"), sample, "-markers.umap.basal-SCT.pdf"), width = 14)
FeaturePlot(rat, features = markers$basal, pt.size = 0.2, ncol = 4)
dev.off()

pdf(paste0(here::here("../../plots/sd-sorted/"), sample, "-markers.umap.luminal-SCT.pdf"), width = 14)
FeaturePlot(rat, features = markers$luminal, pt.size = 0.2, ncol = 4)
dev.off()


## markers from the Khaled (Bach et al) mouse 2017 paper, defining clusters
markers.khaled <- read.table("../../tables/gene-lists/markers-khaled.txt", fill = TRUE, stringsAsFactors = FALSE, header = TRUE)
for (m in 1:ncol(markers.khaled)){
  markers.here <- setdiff(markers.khaled[,m], "")
  pdf(paste0("../../plots/sd-sorted/khaled_", colnames(markers.khaled)[m], ".pdf"), width = 10)
  print(FeaturePlot(rat, features = markers.here))
  dev.off()
}


## marker plots for paper
for (m in c("LALBA", "ESR1", "PGR", "KRT5", "KRT18", "KRT23", "ACTA2", "KIT", "PRLR", "MKI67", "TP63")){
  pdf(paste0("../../plots/", sample, "-sorted/markers_", m, "-paper.pdf"), width = 8)
  print(FeaturePlot(rat, features = m))
  dev.off()
}

markers.paper <- c("SMAD1", "SMAD2", "SMAD3", "SMAD4", "SMAD5", "SMAD7", "SMAD9", "CDKN1B", "ACVR1B", "TGFBR1", "MINK1", 
                   "TGFBR2", "TGFBR3", "RIPK2", "CSNK1A1", "MAP4K4", "GAK", "CSNK1E", "BMPR1B", "BRAF", "TNIK", "ACVR2B", 
                   "RPS6KA6", "ABL1", "MAP3K20", "NLK", "TP63")
setdiff(markers.paper, rownames(rat))
for (m in markers.paper){
  pdf(paste0("../../plots/sd-sorted/markers_", m, "-paper.pdf"), width = 8)
  print(FeaturePlot(rat, features = m, reduction = "UMAP.SCT"))
  dev.off()
}


## split plots for markers
markers <- c("TGFBR1", "TGFBR2", "TGFBR3", "SMAD2", "SMAD3", "SMAD4", "CRIP1", "ID3", "S100A4")
for (m in markers){
  pdf(paste0("../../plots/sd-sorted/markers_", m, ".split.pdf"), width = 12)
  print(FeaturePlot(rat, features = m, split.by = "drug"))
  dev.off()
}


## for secretory basal cells markers
markers <- toupper(c("Col3a1", "Col1a1", "IGFBP7", "Fabp4", "Fgl2", "Dcn", "Sparcl1", "Sparc", "CRIP1", "ID3", "S100A4"))
for (m in markers){
  pdf(paste0("../../plots/sd-sorted/markers_", m, "-paper.pdf"), width = 8)
  print(FeaturePlot(rat, features = m))
  dev.off()
}


## plots for Carlos' stains
pdf(paste0("../../plots/sd-sorted/markers.staining.pdf"), width = 12, height = 12)
FeaturePlot(rat, features = c("S100A4", "S100A6", "CRIP1", "ZEB2", "PROCR", "ID3", "PMP22"))
dev.off()


## plots for TGFb targets
pdf(paste0("../../plots/sd-sorted/markers.tgfb.targets.pdf"), width = 12, height = 12)
FeaturePlot(rat, features = c("TGFBR1", "TGFBR2", "TGFBR3", "SMAD2", "SMAD3", "SMAD4", "SMAD7"))
dev.off()


## plots for TGFb inhibitors
pdf(paste0("../../plots/sd-sorted/markers.tgfb.inhibitors.pdf"), width = 12, height = 12)
FeaturePlot(rat, features = c("LTBP1", "THBS1", "DCN"))
dev.off()


## plot SD animals with changing color
Idents(rat) <- rat@meta.data$animal
pdf(paste0("../../plots/sd-sorted/umap.animal-colors.pdf"), width = 8)
DimPlot(rat, cols = c("animal7" = "#0BE4F4", "animal8" = "#0C4EF4", "animal9" = "#0BF727", "animal12" = "#E90BF4", 
                      "animal15" = "#F98903", "animal16" = "#FC0303"))
dev.off()


## plots for control and drug
Idents(rat) <- rat$drug
pdf(paste0("../../plots/sd-sorted/umap.drug-figure.pdf"))
DimPlot(rat, cols = c("control" = "#C0C0C0", "drug" = "#FED242"))
dev.off()



############################################################################
############################################################################
###                                                                      ###
###                              CELL TYPES                              ###
###                                                                      ###
############################################################################
############################################################################
# check coordinates of this split if redoing
DimPlot(rat, group.by = "celltype")
rat$celltype.detailed <- rat$celltype
rat$celltype.detailed[intersect(which(rat@reductions$UMAP.SCT@cell.embeddings[,1] > 0), 
                                which(rat@reductions$UMAP.SCT@cell.embeddings[,2] < -5))] <- "secretory"
rat$celltype.detailed[intersect(which(rat$celltype.detailed == "luminal"), which(rat@reductions$UMAP.SCT@cell.embeddings[,1] >0))] <- "luminal-secretory"
rat$celltype.detailed[intersect(which(rat@reductions$UMAP.SCT@cell.embeddings[,1] < 0), 
                                which(rat@reductions$UMAP.SCT@cell.embeddings[,2] < 2))] <- "LP"
rat$celltype.detailed[which(rat$celltype.detailed == "luminal")] <- "ML"
pdf("../../plots/sd-sorted/umap-celltype.detailed.pdf", width = 10)
DimPlot(rat, group.by = "celltype.detailed")
dev.off()


## secretory basal signature
sign <- read.table("../../tables/gene-lists/Positive genes in both strains.txt",sep = "\t", header = TRUE)

rat$sign <- colSums(as.matrix(rat@assays$SCT[match(intersect(sign$gene, rownames(rat@assays$SCT)), rownames(rat@assays$SCT)),]))
pdf(paste0("../../plots/sd-sorted/umap-SBC-signature.pdf"))
FeaturePlot(rat, features = c("sign"))
dev.off()

DefaultAssay(rat) <- "RNA"
rat <- NormalizeData(rat, assay = "RNA")


## differential expression between secretory basal and the rest of basal
Idents(rat) <- rat$celltype.detailed
DefaultAssay(rat) <- "RNA"
rat <- NormalizeData(rat, assay = "RNA")
secretory <- FindMarkers(rat, ident.1 = "secretory", ident.2 = "basal", assay = "RNA")
secretory.pos <- secretory %>% filter(avg_log2FC > 0, p_val_adj < 0.1) %>% arrange(desc(abs(avg_log2FC)))
secretory.neg <- secretory %>% filter(avg_log2FC < 0, p_val_adj < 0.1) %>% arrange(desc(abs(avg_log2FC)))


## keratin levels
keratins <- grep("KRT", rownames(rat), value = TRUE)
keratins <- setdiff(keratins, c(grep("KRTAP", keratins, value = TRUE), grep("KRTCAP", keratins, value = TRUE)))
keratins.counts <- colSums(as.matrix(rat@assays$RNA@counts[keratins,]))
rat$keratins.counts <- keratins.counts
rat$log.keratins.counts <- log(rat$keratins.counts + 1)

pdf("../../plots/sd-sorted/log.keratins.counts.density.all.pdf")
ggplot(rat@meta.data, aes(x=log.keratins.counts, color=drug)) +
  geom_density()
dev.off()


## TGFB-related genes dot plot
pdf("../../plots/sd-sorted/dot-plot-TGFB.pdf")
DotPlot(rat, features = c("TGFBR1", "TGFBR2", "TGFBR3", "SMAD2", "SMAD3", "SMAD4", "TGFB1", "TGFB2", "TGFB3"), split.by = "drug", cols = "Spectral", scale = FALSE) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
dev.off()

saveRDS(rat, file = "../../data/rdata/rat-initial-sd.RDS")



###########################################################################
###########################################################################
###                                                                     ###
###                           SBC SUBCLUSTERS                           ###
###                                                                     ###
###########################################################################
###########################################################################
## differential expression between SBC subclusters and the rest
md <- read.table("../../tables/metadata/secretory.meta.data.txt", sep = "\t")
rat$secretory.subclusters <- "basal"
rat$secretory.subclusters[match(rownames(md %>% filter(SNN.INTEGRATE_res.0.2 == 0)), colnames(rat))] <- "sbc0"
rat$secretory.subclusters[match(rownames(md %>% filter(SNN.INTEGRATE_res.0.2 == 1)), colnames(rat))] <- "sbc1"
rat$secretory.subclusters[match(rownames(md %>% filter(SNN.INTEGRATE_res.0.2 == 2)), colnames(rat))] <- "sbc2"
rat$secretory.subclusters[match(rownames(md %>% filter(SNN.INTEGRATE_res.0.2 == 3)), colnames(rat))] <- "sbc3"

Idents(rat) <- rat$secretory.subclusters
sbc.c0 <- FindMarkers(rat, ident.1 = "sbc0", ident.2 = "basal", assay = "RNA")
sbc.c0.pos <- sbc.c0 %>% filter(avg_log2FC > 0, p_val_adj < 0.1) %>% arrange(desc(abs(avg_log2FC)))
sbc.c0.neg <- sbc.c0 %>% filter(avg_log2FC < 0, p_val_adj < 0.1) %>% arrange(desc(abs(avg_log2FC)))
write.table(cbind(rownames(sbc.c0.pos), abs(sbc.c0.pos$avg_log2FC)), file = "~/Downloads/sbc.c0.pos.rnk", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(cbind(rownames(sbc.c0.neg), abs(sbc.c0.neg$avg_log2FC)), file = "~/Downloads/sbc.c0.neg.abs.rnk", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

sbc.c1 <- FindMarkers(rat, ident.1 = "sbc1", ident.2 = "basal", assay = "RNA")
sbc.c1.pos <- sbc.c1 %>% filter(avg_log2FC > 0, p_val_adj < 0.1) %>% arrange(desc(abs(avg_log2FC)))
sbc.c1.neg <- sbc.c1 %>% filter(avg_log2FC < 0, p_val_adj < 0.1) %>% arrange(desc(abs(avg_log2FC)))
write.table(cbind(rownames(sbc.c1.pos), abs(sbc.c1.pos$avg_log2FC)), file = "~/Downloads/sbc.c1.pos.rnk", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(cbind(rownames(sbc.c1.neg), abs(sbc.c1.neg$avg_log2FC)), file = "~/Downloads/sbc.c1.neg.abs.rnk", sep = "\t", row.names = FALSE, quote = FALSE)

sbc.c2 <- FindMarkers(rat, ident.1 = "sbc2", ident.2 = "basal", assay = "RNA")
sbc.c2.pos <- sbc.c2 %>% filter(avg_log2FC > 0, p_val_adj < 0.1) %>% arrange(desc(abs(avg_log2FC)))
sbc.c2.neg <- sbc.c2 %>% filter(avg_log2FC < 0, p_val_adj < 0.1) %>% arrange(desc(abs(avg_log2FC)))
write.table(cbind(rownames(sbc.c2.pos), abs(sbc.c2.pos$avg_log2FC)), file = "~/Downloads/sbc.c2.pos.rnk", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(cbind(rownames(sbc.c2.neg), abs(sbc.c2.neg$avg_log2FC)), file = "~/Downloads/sbc.c2.neg.abs.rnk", sep = "\t", row.names = FALSE, quote = FALSE)

sbc.c3 <- FindMarkers(rat, ident.1 = "sbc3", ident.2 = "basal", assay = "RNA")
sbc.c3.pos <- sbc.c3 %>% filter(avg_log2FC > 0, p_val_adj < 0.1) %>% arrange(desc(abs(avg_log2FC)))
sbc.c3.neg <- sbc.c3 %>% filter(avg_log2FC < 0, p_val_adj < 0.1) %>% arrange(desc(abs(avg_log2FC)))
write.table(cbind(rownames(sbc.c3.pos), abs(sbc.c3.pos$avg_log2FC)), file = "~/Downloads/sbc.c3.pos.rnk", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(cbind(rownames(sbc.c3.neg), abs(sbc.c3.neg$avg_log2FC)), file = "~/Downloads/sbc.c3.neg.abs.rnk", sep = "\t", row.names = FALSE, quote = FALSE)







