## Single cell analysis on the Bach et. al 2017 (last author: Khaled Walid) mouse breast dataset; input data from the authors of the paper
## Author: Simona Cristea
## Date: October 2022
## Publication: Breast cancer prevention by short-term inhibition of TGFβ signaling, Nat Comm 2022


###########################################################################
###########################################################################
###                                                                     ###
###                      READ IN DATA FROM AUTHORS                      ###
###                                                                     ###
###########################################################################
###########################################################################
khaled_data <- readRDS(file = "../../data/rdata/khaled_ExpressionList_QC_norm_clustered_clean.rds")



############################################################################
############################################################################
###                                                                      ###
###                     CREATE SEURAT AND PREPROCESS                     ###
###                                                                      ###
############################################################################
############################################################################
s <- CreateSeuratObject(counts = khaled_data$counts)
khaled_data$phenoData$Cluster <- as.numeric(khaled_data$phenoData$Cluster)
all.equal(khaled_data$phenoData$barcode, colnames(s))
s$IsNonEpithelial <- khaled_data$phenoData$IsNonEpithelial
s$Cluster <- khaled_data$phenoData$Cluster
s$SubCluster <- khaled_data$phenoData$SubCluster
s <- SplitObject(s, split.by = "IsNonEpithelial")
s <- s$`FALSE`
s <- SCTransform(s)
s <- RunPCA(s)
s <- FindNeighbors(s)
s <- FindClusters(s)
s <- RunUMAP(s, dims = 1:30)



###########################################################################
###########################################################################
###                                                                     ###
###                            SBC SIGNATURE                            ###
###                                                                     ###
###########################################################################
###########################################################################
sgn <- read.table("../../tables/gene-lists/Positive genes in both strains.txt", header = TRUE, sep = "\t")
sgn <- sgn$gene
id.sgn <- khaled_data$featureData$id[na.omit(match(sgn, toupper(khaled_data$featureData$symbol)))]

# compute signature in data
s$sign <- colSums(s@assays$RNA@data[match(id.sgn,rownames(s@assays$RNA@data)),])

# remove NAs Cluster
s <- s[,-which(is.na(s$Cluster))]

pdf("../../plots/khaled2017/ridge.subcluster.pdf")
RidgePlot(s, features = c("sign"), group.by = "SubCluster")
dev.off()

pdf("../../plots/khaled2017/umap.pdf")
DimPlot(s, group.by = "SubCluster")
dev.off()

pdf("../../plots/khaled2017/sign.pdf")
FeaturePlot(s, features = c("sign"))
dev.off()


saveRDS(file = "../../data/rdata/khaled_sign.RDS", list("s" = s, "khaled_data" = khaled_data))

