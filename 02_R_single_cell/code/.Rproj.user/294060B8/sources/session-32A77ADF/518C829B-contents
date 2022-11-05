## Single cell analysis on sorted epithelium ACI data
## Author: Simona Cristea
## Date: October 2022
## Publication: Breast cancer prevention by short-term inhibition of TGFβ signaling, Nat Comm 2022

## load libraries
source(here::here("libraries.R"))



############################################################################
############################################################################
###                                                                      ###
###                       READ IN SEQUENCING INPUT                       ###
###                                                                      ###
############################################################################
############################################################################

# replace with your own path where cellranger count outputs 
dataset_loc <- here::here("../../data/cellranger-count/")
samples <- c("ACI2-Bas12", "ACI2-Bas15", "ACI2-Bas16", "ACI2-Bas7", "ACI2-Bas8", "ACI2-Bas9", 
             "ACI2-Lum12", "ACI2-Lum15", "ACI2-Lum16", "ACI2-Lum7", "ACI2-Lum8", "ACI2-Lum9")
rat.data <- sapply(samples, function(i){
  d10x <- Read10X(file.path(dataset_loc,i,"outs/"), )
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})
rat.data <- do.call("cbind", rat.data)
rownames(rat.data) <- toupper(rownames(rat.data))
dim(rat.data)

## create Seurat object
rat <- CreateSeuratObject(
  rat.data,
  project = "all_aci_samples", 
  min.cells = 10,
  min.features = 100,
  names.field = 2,
  names.delim = "\\-")


## extract metadata information from sample name
rat <- AddMetaData(object = rat, metadata = sapply(rownames(rat@meta.data), function(x){paste(strsplit(x, "-")[[1]][2:3], collapse= "-")}), col.name = "sample")
rat <- AddMetaData(object = rat, metadata = recode(tolower(sapply(sapply(rat@meta.data$sample, function(x){paste(strsplit(x, "-")[[1]][2], collapse= "-")}), 
                                                                  function(y){return(gsub("[[:digit:]]","",y))})), "bas" = "basal", "lum" = "luminal"),
                   col.name = "celltype")
rat <- AddMetaData(object = rat, metadata = paste0("animal", sapply(sapply(rat@meta.data$sample, function(x){paste(strsplit(x, "-")[[1]][2], collapse= "-")}), 
                                                                    function(y){return(gsub("[^[:digit:]]","",y))})), col.name = "animal")
rat@meta.data$drug <- sapply(rat@meta.data$sample, function(x){ 
  if (x == "ACI2-Bas12")
    return("drug")
  if (x == "ACI2-Bas15")
    return("drug")
  if (x == "ACI2-Bas16")
    return("drug")
  if (x == "ACI2-Bas7")
    return("control")
  if (x == "ACI2-Bas8")
    return("control")
  if (x == "ACI2-Bas9")
    return("control")
  if (x == "ACI2-Lum12")
    return("drug")
  if (x == "ACI2-Lum15")
    return("drug")
  if (x == "ACI2-Lum16")
    return("drug")
  if (x == "ACI2-Lum7")
    return("control")
  if (x == "ACI2-Lum8")
    return("control")
  if (x == "ACI2-Lum9")
    return("control")
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
saveRDS(rat, file = here::here("../../data/rdata/rat-initial-aci-only-mito.RDS"))


## clean data of contamination before normalization (PTPRC+ and KRT- cells)
immune <- which(as.vector(rat@assays$RNA@counts["PTPRC",]!=0))
length(immune)
keratins <- grep("KRT", rownames(rat), value = TRUE)
keratins <- setdiff(keratins, c(grep("KRTAP", keratins, value = TRUE), grep("KRTCAP", keratins, value = TRUE)))
keratins.vals <- colSums(as.matrix(rat@assays$RNA@counts[keratins,]))
no.krt <- which(keratins.vals < 2)
length(no.krt)
to.exclude <- union(immune, no.krt)
length(to.exclude)
rat <- rat[,-to.exclude]
dim(rat)



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

rat <- rat[,-match(cells.remove, colnames(rat))]
dim(rat)
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

pdf(paste0(here::here("../../plots/aci-sorted/"), sample, "-markers.umap.gate-SCT.pdf"), width = 14)
FeaturePlot(rat, features = markers$gate, pt.size = 0.2)
dev.off()

pdf(paste0(here::here("../../plots/aci-sorted/"), sample, "-markers.umap.basal-SCT.pdf"), width = 14)
FeaturePlot(rat, features = markers$basal, pt.size = 0.2, ncol = 4)
dev.off()

pdf(paste0(here::here("../../plots/aci-sorted/"), sample, "-markers.umap.luminal-SCT.pdf"), width = 14)
FeaturePlot(rat, features = markers$luminal, pt.size = 0.2, ncol = 4)
dev.off()


## markers from the Khaled (Bach et al) mouse 2017 paper, defining clusters
markers.khaled <- read.table("../../tables/gene-lists/markers-khaled.txt", fill = TRUE, stringsAsFactors = FALSE, header = TRUE)
for (m in 1:ncol(markers.khaled)){
  markers.here <- setdiff(markers.khaled[,m], "")
  pdf(paste0("../../plots/aci-sorted/khaled_", colnames(markers.khaled)[m], ".pdf"), width = 10)
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
  pdf(paste0("../../plots/aci-sorted/markers_", m, "-paper.pdf"), width = 8)
  print(FeaturePlot(rat, features = m, reduction = "UMAP.SCT"))
  dev.off()
}


## split plots for markers
markers <- c("TGFBR1", "TGFBR2", "TGFBR3", "SMAD2", "SMAD3", "SMAD4", "CRIP1", "ID3", "S100A4")
for (m in markers){
  pdf(paste0("../../plots/aci-sorted/markers_", m, ".split.pdf"), width = 12)
  print(FeaturePlot(rat, features = m, split.by = "drug"))
  dev.off()
}


## for secretory basal cells markers
markers <- toupper(c("Col3a1", "Col1a1", "IGFBP7", "Fabp4", "Fgl2", "Dcn", "Sparcl1", "Sparc", "CRIP1", "ID3", "S100A4"))
for (m in markers){
  pdf(paste0("../../plots/aci-sorted/markers_", m, "-paper.pdf"), width = 8)
  print(FeaturePlot(rat, features = m))
  dev.off()
}


## plots for Carlos' stains
pdf(paste0("../../plots/aci-sorted/markers.staining.pdf"), width = 12, height = 12)
FeaturePlot(rat, features = c("S100A4", "S100A6", "CRIP1", "ZEB2", "PROCR", "ID3", "PMP22"))
dev.off()


## plots for TGFb targets
pdf(paste0("../../plots/aci-sorted/markers.tgfb.targets.pdf"), width = 12, height = 12)
FeaturePlot(rat, features = c("TGFBR1", "TGFBR2", "TGFBR3", "SMAD2", "SMAD3", "SMAD4", "SMAD7"))
dev.off()


## plots for TGFb inhibitors
pdf(paste0("../../plots/aci-sorted/markers.tgfb.inhibitors.pdf"), width = 12, height = 12)
FeaturePlot(rat, features = c("LTBP1", "THBS1", "DCN"))
dev.off()


## plot ACI animals with changing color
Idents(rat) <- rat@meta.data$animal
pdf(paste0("../../plots/aci-sorted/umap.animal-colors.pdf"), width = 8)
DimPlot(rat, cols = c("animal7" = "#0BE4F4", "animal8" = "#0C4EF4", "animal9" = "#0BF727", "animal12" = "#E90BF4", 
                      "animal15" = "#F98903", "animal16" = "#FC0303"))
dev.off()


## plots for control and drug
Idents(rat) <- rat$drug
pdf(paste0("../../plots/aci-sorted/umap.drug-figure.pdf"))
DimPlot(rat, cols = c("control" = "#C0C0C0", "drug" = "#FED242"))
dev.off()



############################################################################
############################################################################
###                                                                      ###
###                              CELL TYPES                              ###
###                                                                      ###
############################################################################
############################################################################

rat$celltype.detailed <- rat$celltype
DimPlot(rat, group.by = "celltype.detailed")
rat$celltype.detailed[intersect(which(rat@reductions$UMAP.SCT@cell.embeddings[,2]>0), which(rat$celltype == "luminal"))] <- "ML"
rat$celltype.detailed[intersect(which(rat@reductions$UMAP.SCT@cell.embeddings[,2]<0), which(rat$celltype == "luminal"))] <- "LP"
pdf(paste0("../../plots/aci-sorted/umap.celltype-detailed.pdf"))
DimPlot(rat, group.by = "celltype.detailed")
dev.off()


## secretory basal signature
sign <- read.table("../../tables/gene-lists/Positive genes in both strains.txt",sep = "\t", header = TRUE)

rat$sign <- colSums(as.matrix(rat@assays$SCT[match(intersect(sign$gene, rownames(rat@assays$SCT)), rownames(rat@assays$SCT)),]))
pdf(paste0("../../plots/aci-sorted/umap-SBC-signature.pdf"))
FeaturePlot(rat, features = c("sign"))
dev.off()


## keratin levels
keratins <- grep("KRT", rownames(rat), value = TRUE)
keratins <- setdiff(keratins, c(grep("KRTAP", keratins, value = TRUE), grep("KRTCAP", keratins, value = TRUE)))
keratins.counts <- colSums(as.matrix(rat@assays$RNA@counts[keratins,]))
rat$keratins.counts <- keratins.counts
rat$log.keratins.counts <- log(rat$keratins.counts + 1)


pdf("../../plots/aci-sorted/log.keratins.counts.density.all.pdf")
ggplot(rat@meta.data, aes(x=log.keratins.counts, color=drug)) +
  geom_density()
dev.off()


## TGFB-related genes dot plot
pdf("../../plots/aci-sorted/dot-plot-TGFB.pdf")
DotPlot(rat, features = c("TGFBR1", "TGFBR2", "TGFBR3", "SMAD2", "SMAD3", "SMAD4", "TGFB1", "TGFB2", "TGFB3"), split.by = "drug", cols = "Spectral", scale = FALSE) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
dev.off()


saveRDS(rat, file = "../../data/rdata/rat-initial-aci.RDS")



