## monocle3 pseudotime analysis on sorted SD epithlial data (input: data separated by treatment & integrated by animal)
## Author: Simona Cristea
## Date: October 2022
## Publication: Breast cancer prevention by short-term inhibition of TGFβ signaling, Nat Comm 2022


## load libraries
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)



###########################################################################
###########################################################################
###                                                                     ###
###                      READ IN PREPROCESSED DATA                      ###
###                                                                     ###
###########################################################################
###########################################################################
rat <- readRDS(file = "../../data/rdata/rat-sd-integrated-treatment.RDS")
rat <- rat$drug
tt <- "drug"
# rat <- rat$control
# tt <- "control"



###########################################################################
###########################################################################
###                                                                     ###
###                               MONOCLE                               ###
###                                                                     ###
###########################################################################
###########################################################################
cds <- as.cell_data_set(rat)
cds <- estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(rat)
cds <- preprocess_cds(cds, method = "PCA", norm_method = "none")
cds <- reduce_dimension(cds, reduction_method = "UMAP", max_components = 2)
cds <- cluster_cells(cds)
cds <- learn_graph(cds, use_partition = FALSE)
plot_cells(cds)
pdf(paste0("../../plots/trajectory/monocle-", tt, "-markers.pdf"))
plot_cells(cds, gene = c("ACTA2", "SPARCL1", "LALBA", "PGR"), show_trajectory_graph = FALSE, label_cell_groups = FALSE)
dev.off()
pdf(paste0("../../plots/trajectory/monocle-", tt, "-secretory.pdf"))
plot_cells(cds, gene = c("PROCR", "ZEB2", "ZEB1", "SPARCL1"), show_trajectory_graph = FALSE)
dev.off()
pdf(paste0("../../plots/trajectory/monocle-", tt, "-animal.pdf"))
plot_cells(cds, color_cells_by = "animal", show_trajectory_graph = FALSE, label_cell_groups = FALSE)
dev.off()
pdf(paste0("../../plots/trajectory/monocle-", tt, "-percent.mt.pdf"))
plot_cells(cds, color_cells_by = "percent.mt")
dev.off()
pdf(paste0("../../plots/trajectory/monocle-", tt, "-MME.pdf"))
plot_cells(cds, gene = "MME", show_trajectory_graph = FALSE)
dev.off()

cds <- order_cells(cds, reduction_method = "UMAP")
pdf(paste0("../../plots/trajectory/monocle-", tt, "-pseudotime-SBCs-starting.pdf"))
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
dev.off()


cds.seurat <- as.Seurat(cds, assay = "integrated", counts = "counts")

saveRDS(cds.seurat, file = paste0("../../data/rdata/monocle-SD-", tt, "-seurat.RDS"))
