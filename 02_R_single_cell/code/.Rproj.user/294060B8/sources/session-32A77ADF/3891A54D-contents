## Single cell analysis on sorted epithelium data (either ACI or SD), analyzed separately per treatment condition (control/drug) 
## and integrated per animal
## Author: Simona Cristea
## Date: October 2022
## Publication: Breast cancer prevention by short-term inhibition of TGFβ signaling, Nat Comm 2022


## load libraries
source(here::here("libraries.R"))

## read in data
rat <- readRDS(file = "../../data/rdata/rat-initial-sd.RDS")
#rat <- readRDS(file = "../../data/rdata/rat-initial-aci.RDS")



############################################################################
############################################################################
###                                                                      ###
###                      SPLIT OBJECT AND INTEGRATE                      ###
###                                                                      ###
############################################################################
############################################################################
rat <- NormalizeData(rat, assay = "RNA")
rat.drug <- SplitObject(rat, split.by = "drug")
rat.drug <- lapply(rat.drug, FUN = function(x) {
  x <- SplitObject(x, split.by = "animal")
  x <- lapply(x, FUN = function(y) {
    y <- SCTransform(y, verbose = TRUE)
  })
  features <- SelectIntegrationFeatures(object.list = x, nfeatures = 3000)
  x <- PrepSCTIntegration(object.list = x, anchor.features = features, verbose = TRUE)
  anchors <- FindIntegrationAnchors(object.list = x, normalization.method = "SCT", anchor.features = features, verbose = TRUE)
  x <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = TRUE)
  x <- RunPCA(x, verbose = TRUE)
  x <- FindNeighbors(x, dims = 1:30, verbose = FALSE, reduction = "pca", graph.name = "SNN.INTEGRATE")
  x <- FindClusters(x, verbose = TRUE, graph.name = "SNN.INTEGRATE", resolution = 0.8)
  x <- RunUMAP(x, dims = 1:30)
  return(x)
})


rat.drug <- lapply(rat.drug, function(x){
  x <- FindClusters(x, verbose = TRUE, graph.name = "SNN.INTEGRATE", resolution = 0.5)
  x <- FindClusters(x, verbose = TRUE, graph.name = "SNN.INTEGRATE", resolution = 0.6)
  x <- FindClusters(x, verbose = TRUE, graph.name = "SNN.INTEGRATE", resolution = 0.7)
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x, assay = "RNA")
  return(x)
})


## resolutions for the SD data
Idents(rat.drug$drug) <- rat.drug$drug$SNN.INTEGRATE_res.0.7
Idents(rat.drug$control) <- rat.drug$control$SNN.INTEGRATE_res.0.5
rat.drug$control <- subset(rat.drug$control, idents = '10', invert = TRUE)

pdf("../../plots/sd-sorted-integrated/clusters-drug.pdf")
#pdf("../../plots/aci-sorted-integrated/clusters-drug.pdf")
DimPlot(rat.drug$drug, label = TRUE)
dev.off()

pdf("../../plots/sd-sorted-integrated/clusters-control.pdf")
#pdf("../../plots/aci-sorted-integrated/clusters-control.pdf")
DimPlot(rat.drug$control, label = TRUE)
dev.off()

pdf("../../plots/sd-sorted-integrated/markers-drug.pdf")
#pdf("../../plots/aci-sorted-integrated/markers-drug.pdf")
FeaturePlot(rat.drug$drug, features = c("LALBA", "SPARCL1", "ACTA2", "PGR"))
dev.off()


pdf("../../plots/sd-sorted-integrated/markers-control.pdf")
#pdf("../../plots/aci-sorted-integrated/markers-control.pdf")
FeaturePlot(rat.drug$control, features = c("LALBA", "SPARCL1", "ACTA2", "PGR"))
dev.off()


saveRDS(rat.drug, file = "../../data/rdata/rat-sd-integrated-treatment.RDS")
saveRDS(rat.drug, file = "../../data/rdata/rat-aci-integrated-treatment.RDS")