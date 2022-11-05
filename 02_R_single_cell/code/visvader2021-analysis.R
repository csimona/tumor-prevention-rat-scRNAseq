## Single cell analysis on the Pal et. al 2021 (last author: Visvader Jane) human breast dataset; input data from the authors of the paper
## Author: Simona Cristea
## Date: October 2022
## Publication: Breast cancer prevention by short-term inhibition of TGFβ signaling, Nat Comm 2022

## load libraries
library(Seurat)
library(biomaRt)


###########################################################################
###########################################################################
###                                                                     ###
###                      READ IN DATA FROM AUTHORS                      ###
###                                                                     ###
###########################################################################
###########################################################################

fig4b <- readRDS(file = "../../data/rdata/visvader-figure4b-seurat.rds" ) # all cells from 4BRCA1 and 8 normal premenopausal

## recode sample info
fig4b$sample <- dplyr::recode(fig4b$orig.ident, 
                              "KCF0894" = "B1-0894", "MH0023" = "B1-0023", "MH0023Total" = "B1-0023", 
                              "MH0033Total" = "B1-0033", "MH0064Total" = "N-0064","MH0090" = "B1-0090",
                              "MH0169" = "N-0169", "PM0019" = "N-0019", "PM0092" = "N-0092", 
                              "PM0095Total" = "N-0093", "PM0230" = "N-0230.17", "PM0233" = "N-0123")

fig4b$condition <- dplyr::recode(fig4b$orig.ident, 
                                 "KCF0894" = "BRCA1", "MH0023" = "BRCA1", "MH0023Total" = "BRCA1", "MH0033Total" = "BRCA1", 
                                 "MH0064Total" = "N","MH0090" = "BRCA1", "MH0169" = "N", "PM0019" = "N", "PM0092" = "N", 
                                 "PM0095Total" = "N", "PM0230" = "N", "PM0233" = "N")

## recode celltype
fig4b$celltype <- dplyr::recode(fig4b$seurat_clusters, "0" = "non-ep", "1" = "ep", "2" = "non-ep", "3" = "ep", "4" = "ep", "5" = "non-ep", "6" = "non-ep", "7" = "ep", "8" = "non-ep", "9" = "non-ep")
pdf("../../plots/visvader2021/fig4b-assignment.pdf")
DimPlot(fig4b, group.by = "celltype", reduction = "tsne")
dev.off()

pdf("../../plots/visvader2021/fig4b-clusters.pdf")
DimPlot(fig4b, group.by = "seurat_clusters", reduction = "tsne", label = TRUE)
dev.off()


###########################################################################
###########################################################################
###                                                                     ###
###                        READ IN SBC SIGNATURE                        ###
###                                                                     ###
###########################################################################
###########################################################################
genes <- read.table("../../tables/gene-lists/Positive genes in both strains.txt", sep = "\t", header = TRUE)
hmart <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
rmart <- useMart("ensembl",dataset="rnorvegicus_gene_ensembl")
orth <- getBM(c("external_gene_name","hsapiens_homolog_associated_gene_name"), filters="external_gene_name", 
              values = genes$gene, mart = rmart)
hgenes <- orth$hsapiens_homolog_associated_gene_name
hgenes <- setdiff(hgenes, "")
hgenes <- toupper(hgenes)
genes.4b <- intersect(hgenes, rownames(fig4b@assays$RNA))

fig4b <- NormalizeData(fig4b, assay = "RNA")
fig4b$sign <- colSums(fig4b@assays$RNA@data[match(genes.4b, rownames(fig4b@assays$RNA@data)),])



############################################################################
############################################################################
###                                                                      ###
###            SEPARATE THE EPITHELIUM AND REDO NORMALIZATION            ###
###                                                                      ###
############################################################################
############################################################################
fig4b.split <- SplitObject(fig4b, split.by = "celltype")

fig4b.split$ep <- SCTransform(fig4b.split$ep, verbose = TRUE)
fig4b.split$ep <- RunPCA(fig4b.split$ep, features = VariableFeatures(object = fig4b.split$ep, assay = "SCT"), assay = "SCT")
fig4b.split$ep <- FindNeighbors(fig4b.split$ep, dims = 1:30, verbose = TRUE, assay = "SCT", reduction = "pca", graph.name = "SNN.SCT")
fig4b.split$ep <- FindClusters(fig4b.split$ep, verbose = TRUE, graph.name = "SNN.SCT")
fig4b.split$ep <- RunUMAP(fig4b.split$ep, dims = 1:30, verbose = TRUE, reduction.name = "UMAP.EP", assay = "SCT")
fig4b.split$ep <- RunTSNE(fig4b.split$ep, dims = 1:30, verbose = TRUE, reduction.name = "TSNE.EP", assay = "SCT")
DefaultAssay(fig4b.split$ep) <- "RNA"
fig4b.split$ep <- NormalizeData(fig4b.split$ep, assay = "RNA")
genes.fig4b.ep <- intersect(hgenes, rownames(fig4b.split$ep@assays$RNA))
fig4b.split$ep$sign <- colSums(fig4b.split$ep@assays$RNA@data[match(genes.fig4b.ep, rownames(fig4b.split$ep@assays$RNA@data)),])


fig4b.split$ep$sample <- dplyr::recode(fig4b.split$ep$orig.ident,
                                       "KCF0894" = "BRCA1", "MH0023" = "BRCA1", "MH0023Total" = "BRCA1", "MH0033Total" = "BRCA1", 
                                       "MH0064Total" = "N","MH0090" = "BRCA1", "MH0169" = "N", "PM0019" = "N", "PM0092" = "N", 
                                       "PM0095Total" = "N", "PM0230" = "N", "PM0233" = "N")

fig4b.split$ep$condition <- dplyr::recode(fig4b.split$ep$orig.ident,
                                          "KCF0894" = "BRCA1", "MH0023" = "BRCA1", "MH0023Total" = "BRCA1", "MH0033Total" = "BRCA1", 
                                          "MH0064Total" = "N","MH0090" = "BRCA1", "MH0169" = "N", "PM0019" = "N", "PM0092" = "N", 
                                          "PM0095Total" = "N", "PM0230" = "N", "PM0233" = "N")

fig4b.split$ep$age <- dplyr::recode(fig4b.split$ep$orig.ident,
                                    "KCF0894" = "30", "MH0023" = "42", "MH0023Total" = "42", "MH0033Total" = "31", 
                                    "MH0064Total" = "44","MH0090" = "43", "MH0169" = "35", "PM0019" = "21", "PM0092" = "19", 
                                    "PM0095Total" = "22", "PM0230" = "30", "PM0233" = "34")


### plotting
pdf("../../plots/visvader2021/epithelium-clusters.pdf")
DimPlot(fig4b.split$ep, reduction = "UMAP.EP", label = TRUE)
dev.off()

pdf("../../plots/visvader2021/epithelium-sign.pdf")
FeaturePlot(fig4b.split$ep, reduction = "UMAP.EP", features = c("sign"))
dev.off()

pdf("../../plots/visvader2021/epithelium-condition.pdf")
DimPlot(fig4b.split$ep, group.by = c("condition"), reduction = "UMAP.EP")
dev.off()

pdf("../../plots/visvader2021/epithelium-age.pdf")
DimPlot(fig4b.split$ep, group.by = c("age"), reduction = "UMAP.EP")
dev.off()



###########################################################################
###########################################################################
###                                                                     ###
###                          REMOVE KRT- CELLS                          ###
###                                                                     ###
###########################################################################
###########################################################################
keratins <- grep("KRT", rownames(fig4b.split$ep), value = TRUE)
keratins <- setdiff(keratins, c(grep("KRTAP", keratins, value = TRUE), grep("KRTCAP", keratins, value = TRUE)))
keratins.vals <- colSums(as.matrix(fig4b.split$ep@assays$RNA@counts[keratins,]))
no.krt <- which(keratins.vals < 2)
length(no.krt)
fig4b.split$ep$keratins <- 1
fig4b.split$ep$keratins[no.krt] <- 0
fig4b.split$ep$keratins.vals <- keratins.vals

fig4b.ep.clean <- SplitObject(fig4b.split$ep, split.by = "keratins")
fig4b.ep.clean <- fig4b.ep.clean$`1`
fig4b.ep.clean$sample <- dplyr::recode(fig4b.ep.clean$orig.ident,
                                       "KCF0894" = "KCF0894", "MH0023" = "MH0023", "MH0023Total" = "MH0023", "MH0033Total" = "MH0033Total", 
                                       "MH0064Total" = "MH0064Total","MH0090" = "MH0090", "MH0169" = "MH0169", "PM0019" = "PM0019", "PM0092" = "PM0092", 
                                       "PM0095Total" = "PM0095Total", "PM0230" = "PM0230", "PM0233" = "PM0233")

fig4b.ep.clean$condition <- dplyr::recode(fig4b.ep.clean$orig.ident,
                                          "KCF0894" = "BRCA1", "MH0023" = "BRCA1", "MH0023Total" = "BRCA1", "MH0033Total" = "BRCA1", 
                                          "MH0064Total" = "N","MH0090" = "BRCA1", "MH0169" = "N", "PM0019" = "N", "PM0092" = "N", 
                                          "PM0095Total" = "N", "PM0230" = "N", "PM0233" = "N")

fig4b.ep.clean$age <- dplyr::recode(fig4b.ep.clean$orig.ident,
                                    "KCF0894" = "30", "MH0023" = "42", "MH0023Total" = "42", "MH0033Total" = "31", 
                                    "MH0064Total" = "44","MH0090" = "43", "MH0169" = "35", "PM0019" = "21", "PM0092" = "19", 
                                    "PM0095Total" = "22", "PM0230" = "30", "PM0233" = "34")


fig4b.ep.clean$sample <- dplyr::recode(fig4b.ep.clean$orig.ident,
                                       "KCF0894" = "KCF0894", "MH0023" = "MH0023", "MH0023Total" = "MH0023", "MH0033Total" = "MH0033Total", 
                                       "MH0064Total" = "MH0064Total","MH0090" = "MH0090", "MH0169" = "MH0169", "PM0019" = "PM0019", "PM0092" = "PM0092", 
                                       "PM0095Total" = "PM0095Total", "PM0230" = "PM0230", "PM0233" = "PM0233")

fig4b.ep.clean$condition <- dplyr::recode(fig4b.ep.clean$orig.ident,
                                          "KCF0894" = "BRCA1", "MH0023" = "BRCA1", "MH0023Total" = "BRCA1", "MH0033Total" = "BRCA1", 
                                          "MH0064Total" = "N","MH0090" = "BRCA1", "MH0169" = "N", "PM0019" = "N", "PM0092" = "N", 
                                          "PM0095Total" = "N", "PM0230" = "N", "PM0233" = "N")

fig4b.ep.clean$age <- dplyr::recode(fig4b.ep.clean$orig.ident,
                                    "KCF0894" = "30", "MH0023" = "42", "MH0023Total" = "42", "MH0033Total" = "31", 
                                    "MH0064Total" = "44","MH0090" = "43", "MH0169" = "35", "PM0019" = "21", "PM0092" = "19", 
                                    "PM0095Total" = "22", "PM0230" = "30", "PM0233" = "34")




############################################################################
############################################################################
###                                                                      ###
###                          REDO NORMALIZATION                          ###
###                                                                      ###
############################################################################
############################################################################
fig4b.ep.clean <- SCTransform(fig4b.ep.clean, verbose = TRUE)
fig4b.ep.clean <- RunPCA(fig4b.ep.clean, features = VariableFeatures(object = fig4b.ep.clean, assay = "SCT"), assay = "SCT")
fig4b.ep.clean <- FindNeighbors(fig4b.ep.clean, dims = 1:30, verbose = TRUE, assay = "SCT", reduction = "pca", graph.name = "SNN.SCT")
fig4b.ep.clean <- FindClusters(fig4b.ep.clean, verbose = TRUE, graph.name = "SNN.SCT")
fig4b.ep.clean <- RunUMAP(fig4b.ep.clean, dims = 1:30, verbose = TRUE, reduction.name = "UMAP.EP.CLEAN", assay = "SCT")
fig4b.ep.clean <- RunTSNE(fig4b.ep.clean, dims = 1:30, verbose = TRUE, reduction.name = "TSNE.EP.CLEAN", assay = "SCT")
fig4b.ep.clean <- RunUMAP(fig4b.ep.clean, dims = 1:30, verbose = TRUE, reduction.name = "UMAP.EP.CLEAN.RNA", assay = "RNA")
fig4b.ep.clean <- RunTSNE(fig4b.ep.clean, dims = 1:30, verbose = TRUE, reduction.name = "TSNE.EP.CLEAN.RNA", assay = "RNA")
DefaultAssay(fig4b.ep.clean) <- "RNA"
fig4b.ep.clean <- NormalizeData(fig4b.ep.clean, assay = "RNA")
genes.fig4b.ep.clean <- intersect(hgenes, rownames(fig4b.ep.clean@assays$RNA))
fig4b.ep.clean$sign <- colSums(fig4b.ep.clean@assays$RNA@data[match(genes.fig4b.ep, rownames(fig4b.ep.clean@assays$RNA@data)),])

FeaturePlot(fig4b.ep.clean, features = c("sign"), reduction = "UMAP.EP.CLEAN.RNA")


### plotting
pdf("../../plots/visvader2021/epithelium-clean-clusters.pdf")
DimPlot(fig4b.ep.clean, reduction = "UMAP.EP.CLEAN.RNA", label = TRUE)
dev.off()

pdf("../../plots/visvader2021/epithelium-clean-sign.pdf")
FeaturePlot(fig4b.ep.clean, reduction = "UMAP.EP.CLEAN.RNA", features = c("sign"))
dev.off()

pdf("../../plots/visvader2021/epithelium-clean-SBC-markers.pdf")
FeaturePlot(fig4b.ep.clean, reduction = "UMAP.EP.CLEAN.RNA", features = c("TWIST1", "ZEB1", "ZEB2", "SPARCL1", "NRG1", "LRG1"))
dev.off()

pdf("../../plots/visvader2021/epithelium-clean-markers.pdf")
FeaturePlot(fig4b.ep.clean, reduction = "UMAP.EP.CLEAN.RNA", features = c("LUM", "IGFBP7", "TWIST2", "DCN"))
dev.off()

pdf("../../plots/visvader2021/epithelium-clean-condition.pdf")
DimPlot(fig4b.ep.clean, reduction = "UMAP.EP.CLEAN.RNA", group.by = "condition")
dev.off()

pdf("../../plots/visvader2021/epithelium-clean-sample.pdf")
DimPlot(fig4b.ep.clean, reduction = "UMAP.EP.CLEAN.RNA", group.by = "sample")
dev.off()

pdf("../../plots/visvader2021/epithelium-clean-sample.pdf")
DimPlot(fig4b.ep.clean, reduction = "UMAP.EP.CLEAN.RNA", group.by = "orig.ident")
dev.off()

pdf("../../plots/visvader2021/epithelium-clean-age.pdf")
DimPlot(fig4b.ep.clean, reduction = "UMAP.EP.CLEAN.RNA", group.by = "age")
dev.off()

pdf("../../plots/visvader2021/epithelium-clean-ridge-sign.pdf")
RidgePlot(fig4b.ep.clean, features = c("sign"), group.by = "seurat_clusters")
dev.off()

fig4b.ep.clean$celltype <- dplyr::recode(fig4b.ep.clean$seurat_clusters, 
                                         "0" = "ML", "1" = "basal", "2" = "LP", "3" = "LP", 
                                         "4" = "basal", "5" = "ML",
                                         "6" = "LP", "7" = "basal", "8" = "LP", 
                                         "9" = "LP", "10" = "LP", "11" = "LP", "12" = "LP",
                                         "13" = "LP", "14" = "secretory", "15" = "LP", "16" = "ML", 
                                         "17" = "basal", "18" = "ML")


pdf("../../plots/visvader2021/epithelium-clean-celltype.pdf")
DimPlot(fig4b.ep.clean, group.by = c("celltype"), reduction = "UMAP.EP.CLEAN.RNA")
dev.off()

pdf("../../plots/visvader2021/epithelium-clean-ML-markers.pdf")
FeaturePlot(fig4b.ep.clean, features = c("PGR", "PRLR", "ESR1", "FOXA1"), reduction = "UMAP.EP.CLEAN.RNA")
dev.off()

pdf("../../plots/visvader2021/epithelium-clean-LP-markers.pdf")
FeaturePlot(fig4b.ep.clean, features = c("TNFRSF11A", "KIT", "SOX10"), reduction = "UMAP.EP.CLEAN.RNA")
dev.off()

pdf("../../plots/visvader2021/epithelium-clean-basal-markers.pdf")
FeaturePlot(fig4b.ep.clean, features = c("KRT5", "ACTA2", "MYLK", "SNAI2", "NOTCH4", "DKK3"), reduction = "UMAP.EP.CLEAN.RNA")
dev.off()


## save final Robject
saveRDS(fig4b.ep.clean, file = "../../data/rdata/visvader-epithelium-clean.RDS")
