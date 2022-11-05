## diffusion map analysis on sorted SD epithlial data (input: data with pseudotime trajectories inferred, on Seurat
## separated by treatment & integrated by animal: Seurat --> monocle --> Seurat; input here is the output of monocle.R script)
## Author: Simona Cristea
## Date: October 2022
## Publication: Breast cancer prevention by short-term inhibition of TGFβ signaling, Nat Comm 2022


library(Seurat)
library(destiny)
library(tidyverse)
tt <- "control"
# tt <- "drug
rat <- readRDS(paste0("../../data/rdata/monocle-SD-", tt, "-seurat.RDS"))



############################################################################
############################################################################
###                                                                      ###
###                          CREATE ESET OBJECT                          ###
###                                                                      ###
############################################################################
############################################################################
ann.data <- rat@meta.data
ann.data$sample <- as.factor(ann.data$sample)
ann.data$celltype <- as.factor(ann.data$celltype)
ann.data$animal <- as.factor(ann.data$animal)
ann.data$percent.mt <- as.factor(ann.data$percent.mt)
ann.data$nFeature_SCT <- as.factor(ann.data$nFeature_SCT)
ann.data$celltype.detailed <- as.factor(ann.data$celltype.detailed)
ann.data$celltype.secretory.cluster <- as.factor(ann.data$celltype.secretory.cluster)
ann.data$monocle3_pseudotime <- as.factor(ann.data$monocle3_pseudotime)
ann.data$SNN.INTEGRATE_res.0.8 <- as.factor(ann.data$SNN.INTEGRATE_res.0.8)
m <- as.matrix(rat@assays$integrated@data)
m.short <- m
m.ann <- cbind(t(m.short), ann.data)
m.eset <- as.ExpressionSet(as.data.frame(m.ann), annotation_cols = colnames(ann.data))



###########################################################################
###########################################################################
###                                                                     ###
###                          RUN DIFFUSION MAP                          ###
###                                                                     ###
###########################################################################
###########################################################################
dm <- DiffusionMap(m.eset, n_pcs = 50)
pdf(paste0("../../plots/trajectory/diff-map-",tt, "-animal.pdf"))
plot(dm, col_by = "animal")
dev.off()
pdf(paste0("../../plots/trajectory/diff-map-",tt, "-celltype.detailed.pdf"))
plot(dm, col_by = "celltype.detailed")
dev.off()
pdf(paste0("../../plots/trajectory/diff-map-",tt, "-SNN.INTEGRATE_res.0.8.pdf"))
plot(dm, col_by = "SNN.INTEGRATE_res.0.8")
dev.off()
pdf(paste0("../../plots/trajectory/diff-map-",tt, "-monocle3_pseudotime.pdf"))
plot(dm, col_by = "monocle3_pseudotime")
dev.off()
pdf(paste0("../../plots/trajectory/diff-map-",tt, "-monocle3_pseudotime.pdf"))
plot(dm, col_by = "monocle3_pseudotime")
dev.off()
pdf(paste0("../../plots/trajectory/diff-map-",tt, "-celltype.secretory.cluster.pdf"))
plot(dm, col_by = "celltype.secretory.cluster")
dev.off()


saveRDS(list(dpt, dm), file = paste0("../../data/rdata/destiny-sd4-", tt, ".RDS"))