## Single cell analysis on unfiltered (no PTPRC+, nor KRT- cells removed) sorted epithelium ACI/SD data
## Author: Simona Cristea
## Date: October 2022
## Publication: Breast cancer prevention by short-term inhibition of TGFβ signaling, Nat Comm 2022

source(here::here("libraries.R"))
rat <- readRDS(file="../../data/rdata/rat-initial-aci-only-mito.RDS")
#rat <- readRDS(file="../../data/rdata/rat-initial-sd-only-mito.RDS")



###########################################################################
###########################################################################
###                                                                     ###
###                              INTEGRATE                              ###
###                                                                     ###
###########################################################################
###########################################################################
options(future.globals.maxSize = 4000 * 1024^2)
rat <- SplitObject(rat, split.by = "animal")
rat <- lapply(rat, FUN = function(x) {
  x <- SCTransform(x, verbose = TRUE)
})

integrated.features <- SelectIntegrationFeatures(object.list = rat, nfeatures = 3000)
rat <- PrepSCTIntegration(object.list = rat, anchor.features = integrated.features, verbose = TRUE)
rat <- FindIntegrationAnchors(object.list = rat, normalization.method = "SCT", 
                              anchor.features = integrated.features, verbose = TRUE)
rat <- IntegrateData(anchorset = rat, normalization.method = "SCT", verbose = TRUE)
DefaultAssay(rat) <- "integrated"
rat <- RunPCA(rat, verbose = TRUE)
rat <- FindNeighbors(rat, dims = 1:30, verbose = FALSE, reduction = "pca", graph.name = "SNN.INTEGRATE")
rat <- FindClusters(rat, verbose = TRUE, graph.name = "SNN.INTEGRATE")
rat <- RunUMAP(rat, dims = 1:30)


## keratin levels
keratins <- grep("KRT", rownames(rat), value = TRUE)
keratins <- setdiff(keratins, c(grep("KRTAP", keratins, value = TRUE), grep("KRTCAP", keratins, value = TRUE)))
keratins.counts <- colSums(as.matrix(rat@assays$RNA@counts[keratins,]))
rat$keratins.counts <- keratins.counts
rat$less.than.two.keratins <- 0
rat$less.than.two.keratins[which(keratins.counts < 2)] <- 1


saveRDS(rat, file = "../../data/rdata/rat-aci-integrated-unfiltered.RDS")
#saveRDS(rat, file = "../../data/rdata/rat-sd-integrated-unfiltered.RDS")
