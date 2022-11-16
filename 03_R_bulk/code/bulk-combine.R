## RNASeq analysis of epithelial bulk samples under different treatment/conditions
## Author: Simona Cristea
## Date: October 2022
## Publication: Breast cancer prevention by short-term inhibition of TGFβ signaling, Nat Comm 2022


library(tidyverse)
library(DESeq2)
library(lattice)



###########################################################################
###########################################################################
###                                                                     ###
###                     PREPARE DATA AND PHENOTYPES                     ###
###                                                                     ###
###########################################################################
###########################################################################
# read in data
filenames <- list.files(path = "../../data/htseq-count/", pattern = "\\.txt$", full.names = TRUE)
sample.names <- sapply(filenames, function(x){q <- strsplit(x, ".counts.txt"); return(q[[1]][1])})
sample.names <- sapply(sample.names, function(x){q <- strsplit(x, "//"); return(q[[1]][2])})
names(filenames) <- sample.names
data <- lapply(filenames, function(x){read.table(x, header = FALSE, sep = "\t")})

# check that gene names are all the same
gene.names <- lapply(data, function(x){x[,1]})

# transform to matrix
data <- sapply(data, function(x){x[,2]})
rownames(data) <- toupper(gene.names[[1]])

# read in phenotypes
phenos <- read.table("../../tables/metadata/samples.bulk.txt", sep = "\t", header = TRUE)
phenos <- phenos %>% filter(note %in% c("","normal")) %>% select(-note)



############################################################################
############################################################################
###                                                                      ###
###         SEPARATE PER EXPERIMENT AND PER CELL TYPE AND FILTER         ###
###                                                                      ###
############################################################################
############################################################################
## exp3
phenos.exp3 <- list()
data.exp3 <- list()
i <- 0
for (c in c("Basal", "Luminal", "CD45")){
  i <- i+1
  phenos.exp3[[i]] <- phenos %>% filter(experiment == "exp3", celltype == c)
  phenos.exp3[[i]]$merged <- paste(phenos.exp3[[i]]$condition, phenos.exp3[[i]]$treatment, sep = ".")
  phenos.exp3[[i]]$merged <- as.factor(phenos.exp3[[i]]$merged)
  
  data.exp3[[i]] <- data[, match(phenos.exp3[[i]]$Sample, colnames(data))]
  data.exp3[[i]] <- data.exp3[[i]][which(rowSums(data.exp3[[i]]) > 10),]
}
names(phenos.exp3) <- c("Basal", "Luminal", "CD45")
names(data.exp3) <- c("Basal", "Luminal", "CD45")


## exp7
phenos.exp7 <- list()
data.exp7 <- list()
i <- 0
for (c in c("Luminal", "CD45")){
  i <- i+1
  phenos.exp7[[i]] <- phenos %>% filter(experiment == "exp7", celltype == c)
  phenos.exp7[[i]]$merged <- paste(phenos.exp7[[i]]$condition, phenos.exp7[[i]]$treatment, sep = ".")
  phenos.exp7[[i]]$merged <- as.factor(phenos.exp7[[i]]$merged)
  
  data.exp7[[i]] <- data[, match(phenos.exp7[[i]]$Sample, colnames(data))]
  data.exp7[[i]] <- data.exp7[[i]][which(rowSums(data.exp7[[i]]) > 10),]
}
names(phenos.exp7) <- c("Luminal", "CD45")
names(data.exp7) <- c("Luminal", "CD45")



###########################################################################
###########################################################################
###                                                                     ###
###                  NORMALIZE ALL SAMPLES WITH DESEQ2                  ###
###                                                                     ###
###########################################################################
###########################################################################
## exp3 (basal and luminal)
phenos.all.exp3 <- phenos %>% filter(experiment == "exp3", celltype %in% c("Basal", "Luminal"))
phenos.all.exp3$merged <- paste(phenos.all.exp3$celltype, phenos.all.exp3$condition, phenos.all.exp3$treatment, sep = ".")
phenos.all.exp3$merged <- as.factor(phenos.all.exp3$merged)
phenos.all.exp3$colors <- recode(phenos.all.exp3$merged, 
                                 "Basal.parous.C" = "#F708F7", "Basal.parous.D" = "#FB84FB", 
                                 "Basal.virgin.C" = "#F70606", "Basal.virgin.D" = "#FDC1C1",
                                 "Luminal.parous.C" = "#08D5F7", "Luminal.parous.D" = "#84EAFB",
                                 "Luminal.virgin.C" = "#0E08F9", "Luminal.virgin.D" = "#C3C2FE")
phenos.all.exp3$colors <- as.factor(phenos.all.exp3$colors)
data.all.exp3 <- data[, match(phenos.all.exp3$Sample, colnames(data))]
data.all.exp3 <- data.all.exp3[which(rowSums(data.all.exp3) > 10),]
dds.all.exp3 <- DESeqDataSetFromMatrix(countData = data.all.exp3, colData = phenos.all.exp3, design = ~ merged)
dds.all.exp3 <- estimateSizeFactors(dds.all.exp3)
vsd.all.exp3 <- varianceStabilizingTransformation(dds.all.exp3)
normalizedTableVSD.all.exp3 <- assay(vsd.all.exp3)
normalizedTableVSD.all.exp3 <- normalizedTableVSD.all.exp3[,sort(colnames(normalizedTableVSD.all.exp3))]
normalizedCountsTable.all.exp3 <- as.data.frame(counts(dds.all.exp3, normalized = TRUE))


## exp3 (basal)
phenos.basal.exp3 <- phenos %>% filter(experiment == "exp3", celltype == "Basal")
phenos.basal.exp3$merged <- paste(phenos.basal.exp3$celltype, phenos.basal.exp3$condition, phenos.basal.exp3$treatment, sep = ".")
phenos.basal.exp3$merged <- as.factor(phenos.basal.exp3$merged)
phenos.basal.exp3$colors <- recode(phenos.basal.exp3$merged, 
                                   "Basal.parous.C" = "#F708F7", "Basal.parous.D" = "#FB84FB", 
                                   "Basal.virgin.C" = "#F70606", "Basal.virgin.D" = "#FDC1C1")
phenos.basal.exp3$colors <- as.factor(phenos.basal.exp3$colors)
data.basal.exp3 <- data[, match(phenos.basal.exp3$Sample, colnames(data))]
data.basal.exp3 <- data.basal.exp3[which(rowSums(data.basal.exp3) > 10),]
dds.basal.exp3 <- DESeqDataSetFromMatrix(countData = data.basal.exp3, colData = phenos.basal.exp3, design = ~ merged)
dds.basal.exp3 <- estimateSizeFactors(dds.basal.exp3)
vsd.basal.exp3 <- varianceStabilizingTransformation(dds.basal.exp3)
normalizedTableVSD.basal.exp3 <- assay(vsd.basal.exp3)
normalizedTableVSD.basal.exp3 <- normalizedTableVSD.basal.exp3[,sort(colnames(normalizedTableVSD.basal.exp3))]
normalizedCountsTable.basal.exp3 <- as.data.frame(counts(dds.basal.exp3, normalized = TRUE))


## exp3 (luminal)
phenos.luminal.exp3 <- phenos %>% filter(experiment == "exp3", celltype == "Luminal")
phenos.luminal.exp3$merged <- paste(phenos.luminal.exp3$celltype, phenos.luminal.exp3$condition, phenos.luminal.exp3$treatment, sep = ".")
phenos.luminal.exp3$merged <- as.factor(phenos.luminal.exp3$merged)
phenos.luminal.exp3$colors <- recode(phenos.luminal.exp3$merged, 
                                     "Luminal.parous.C" = "#08D5F7", "Luminal.parous.D" = "#84EAFB",
                                     "Luminal.virgin.C" = "#0E08F9", "Luminal.virgin.D" = "#C3C2FE")
phenos.luminal.exp3$colors <- as.factor(phenos.luminal.exp3$colors)
data.luminal.exp3 <- data[, match(phenos.luminal.exp3$Sample, colnames(data))]
data.luminal.exp3 <- data.luminal.exp3[which(rowSums(data.luminal.exp3) > 10),]
dds.luminal.exp3 <- DESeqDataSetFromMatrix(countData = data.luminal.exp3, colData = phenos.luminal.exp3, design = ~ merged)
dds.luminal.exp3 <- estimateSizeFactors(dds.luminal.exp3)
vsd.luminal.exp3 <- varianceStabilizingTransformation(dds.luminal.exp3)
normalizedTableVSD.luminal.exp3 <- assay(vsd.luminal.exp3)
normalizedTableVSD.luminal.exp3 <- normalizedTableVSD.luminal.exp3[,sort(colnames(normalizedTableVSD.luminal.exp3))]
normalizedCountsTable.luminal.exp3 <- as.data.frame(counts(dds.luminal.exp3, normalized = TRUE))




## exp7 (luminal)
phenos.all.exp7 <- phenos %>% filter(experiment == "exp7", celltype %in% c("Luminal"))
phenos.all.exp7$merged <- paste(phenos.all.exp7$celltype, phenos.all.exp7$condition, phenos.all.exp7$treatment, sep = ".")
phenos.all.exp7$merged <- as.factor(phenos.all.exp7$merged)
phenos.all.exp7$colors <- recode(phenos.all.exp7$merged, 
                                 "Luminal.mE2.C.TF" = "#0E08F9", "Luminal.mE2.D.TF" = "#C3C2FE", 
                                 "Luminal.pE2.C.TB" = "#635FFB", "Luminal.pE2.D.TB" = "#84EAFB",
                                 "Luminal.pE2.D.TF" = "#08D5F7")
data.all.exp7 <- data[, match(phenos.all.exp7$Sample, colnames(data))]
data.all.exp7 <- data.all.exp7[which(rowSums(data.all.exp7) > 10),]
dds.all.exp7 <- DESeqDataSetFromMatrix(countData = data.all.exp7, colData = phenos.all.exp7, design = ~ merged)
dds.all.exp7 <- estimateSizeFactors(dds.all.exp7)
vsd.all.exp7 <- varianceStabilizingTransformation(dds.all.exp7)
normalizedTableVSD.all.exp7 <- assay(vsd.all.exp7)
normalizedTableVSD.all.exp7 <- normalizedTableVSD.all.exp7[,sort(colnames(normalizedTableVSD.all.exp7))]
normalizedCountsTable.all.exp7 <- as.data.frame(counts(dds.all.exp7, normalized = TRUE))



###########################################################################
###########################################################################
###                                                                     ###
###                       PCA PLOT OF ALL SAMPLES                       ###
###                                                                     ###
###########################################################################
###########################################################################
## exp3
pca.all.exp3 <- prcomp(t(assay(vsd.all.exp3)))
importanceIndividual.all.exp3 <- summary(pca.all.exp3)$importance[2,]
importanceCum.all.exp3 <- summary(pca.all.exp3)$importance[3,]
pcaDF.all.exp3 <- as.data.frame(pca.all.exp3$x)
pcaDF.all.exp3 <- pcaDF.all.exp3[,1:3]
pcaDF.all.exp3$Sample <- colnames(vsd.all.exp3)
pcaDF.all.exp3$colors <- vsd.all.exp3$colors
pcaDF.all.exp3 <- pcaDF.all.exp3[sort(rownames(pcaDF.all.exp3)),]
pcaDF.all.exp3$importance1 <- importanceIndividual.all.exp3
pcaDF.all.exp3$importance2 <- importanceCum.all.exp3

x <- c(1:length(pcaDF.all.exp3$importance1))

colorList <- list(col = as.character(pcaDF.all.exp3$colors))
textList <- list(as.character(pcaDF.all.exp3$Sample))
plotPC1PC2 <- xyplot(PC1 ~ PC2, data = pcaDF.all.exp3, pch = 16, cex = 2,  col = as.vector(pcaDF.all.exp3$colors), scales=list(tck=c(1,1), x=list(cex=2), y=list(cex=2)), xlab = list(cex = 2), ylab = list(cex = 2), key = list(rect = colorList, text = textList, columns=2, space = 'top'))
plotPC1PC2 <- xyplot(PC1 ~ PC2, data = pcaDF.all.exp3, pch = 16, cex = 2,  
                     col = as.vector(pcaDF.all.exp3$colors), scales=list(tck=c(1,1), x=list(cex=2), y=list(cex=2)), 
                     xlab = list(cex = 2), ylab = list(cex = 2),
                     key = list(rect = list("col" = c("#F708F7", "#FB84FB", "#F70606", "#FDC1C1", "#08D5F7", "#84EAFB", "#0E08F9", "#C3C2FE")), 
                                text = list("samples" = c("Basal.parous.C", "Basal.parous.D", "Basal.virgin.C", "Basal.virgin.D",
                                                          "Luminal.parous.C", "Luminal.parous.D","Luminal.virgin.C", "Luminal.virgin.D")), 
                                columns=2, space = 'top'))
plotPC1PC3 <- xyplot(PC1 ~ PC3, data = pcaDF.all.exp3, pch = 16, cex = 2,  col = as.vector(pcaDF.all.exp3$colors), scales=list(tck=c(1,1), x=list(cex=2), y=list(cex=2)), xlab = list(cex = 2), ylab = list(cex = 2),
                     key = list(rect = list("col" = c("#F708F7", "#FB84FB", "#F70606", "#FDC1C1", "#08D5F7", "#84EAFB", "#0E08F9", "#C3C2FE")), 
                                text = list("samples" = c("Basal.parous.C", "Basal.parous.D", "Basal.virgin.C", "Basal.virgin.D",
                                                          "Luminal.parous.C", "Luminal.parous.D","Luminal.virgin.C", "Luminal.virgin.D")), 
                                columns=2, space = 'top'))
plotPC2PC3 <- xyplot(PC2 ~ PC3, data = pcaDF.all.exp3, pch = 16, cex = 2,  col = as.vector(pcaDF.all.exp3$colors), scales=list(tck=c(1,1), x=list(cex=2), y=list(cex=2)), xlab = list(cex = 2), ylab = list(cex = 2), 
                     key = list(rect = list("col" = c("#F708F7", "#FB84FB", "#F70606", "#FDC1C1", "#08D5F7", "#84EAFB", "#0E08F9", "#C3C2FE")), 
                                text = list("samples" = c("Basal.parous.C", "Basal.parous.D", "Basal.virgin.C", "Basal.virgin.D",
                                                          "Luminal.parous.C", "Luminal.parous.D","Luminal.virgin.C", "Luminal.virgin.D")), 
                                columns=2, space = 'top'))
pdf("../../plots/bulk/exp3/plotPC1PC2.pdf", width = 10)
try(plotPC1PC2, silent = TRUE)
dev.off()
pdf("../../plots/bulk/exp3/plotPC1PC3.pdf", width = 10)
try(plotPC1PC3, silent = TRUE)
dev.off()
pdf("../../plots/bulk/exp3/plotPC2PC3.pdf", width = 10)
try(plotPC2PC3, silent = TRUE)
dev.off()
pdf("../../plots/bulk/exp3/PCcomponents.pdf", width = 10)
plot(x, pcaDF.all.exp3$importance1, type="b", pch=18, col="red", xlab="Principal Component", ylab="", ylim=c(0,1))
lines(x, pcaDF.all.exp3$importance2, type="b", pch=18, col="blue", lty=2)
legend((length(x)-10), 0.5, legend=c("individual","cumulative"), col=c("red","blue"), lty=1:2, cex=0.8)
dev.off()

plot3d(pcaDF.all.exp3$`PC1`, pcaDF.all.exp3$`PC2`, pcaDF.all.exp3$`PC3`, col = pcaDF.all.exp3$colors, size = 2, type = 's', xlab = 'PC1', ylab = 'PC2', zlab = 'PC3')
q <- rglwidget()
htmlwidgets::saveWidget(q, file = '../../plots/bulk/exp3/all_samples_pca_3D.html')

text3d(pcaDF.all.exp3$`PC1`, pcaDF.all.exp3$`PC2`, pcaDF.all.exp3$`PC3`, pcaDF.all.exp3$Sample, adj = c(-0.1,-0.8))
q <- rglwidget()
htmlwidgets::saveWidget(q, file = '../../plots/bulk/exp3/all_samples_pca_3D_with_label.html')


## exp3 basal
pca.basal.exp3 <- prcomp(t(assay(vsd.basal.exp3)))
importanceIndividual.basal.exp3 <- summary(pca.basal.exp3)$importance[2,]
importanceCum.basal.exp3 <- summary(pca.basal.exp3)$importance[3,]
pcaDF.basal.exp3 <- as.data.frame(pca.basal.exp3$x)
pcaDF.basal.exp3 <- pcaDF.basal.exp3[,1:3]
pcaDF.basal.exp3$Sample <- colnames(vsd.basal.exp3)
pcaDF.basal.exp3$colors <- vsd.basal.exp3$colors
pcaDF.basal.exp3 <- pcaDF.basal.exp3[sort(rownames(pcaDF.basal.exp3)),]
pcaDF.basal.exp3$importance1 <- importanceIndividual.basal.exp3
pcaDF.basal.exp3$importance2 <- importanceCum.basal.exp3

x <- c(1:length(pcaDF.basal.exp3$importance1))

colorList <- list(col = as.character(pcaDF.basal.exp3$colors))
textList <- list(as.character(pcaDF.basal.exp3$Sample))
plotPC1PC2 <- xyplot(PC1 ~ PC2, data = pcaDF.basal.exp3, pch = 16, cex = 2,  col = as.vector(pcaDF.basal.exp3$colors), scales=list(tck=c(1,1), x=list(cex=2), y=list(cex=2)), xlab = list(cex = 2), ylab = list(cex = 2), key = list(rect = colorList, text = textList, columns=2, space = 'top'))
plotPC1PC2 <- xyplot(PC1 ~ PC2, data = pcaDF.basal.exp3, pch = 16, cex = 2,  
                     col = as.vector(pcaDF.basal.exp3$colors), scales=list(tck=c(1,1), x=list(cex=2), y=list(cex=2)), 
                     xlab = list(cex = 2), ylab = list(cex = 2),
                     key = list(rect = list("col" = c("#F708F7", "#FB84FB", "#F70606", "#FDC1C1")), 
                                text = list("samples" = c("Basal.parous.C", "Basal.parous.D", "Basal.virgin.C", "Basal.virgin.D")), 
                                columns=2, space = 'top'))
plotPC1PC3 <- xyplot(PC1 ~ PC3, data = pcaDF.basal.exp3, pch = 16, cex = 2,  col = as.vector(pcaDF.basal.exp3$colors), scales=list(tck=c(1,1), x=list(cex=2), y=list(cex=2)), xlab = list(cex = 2), ylab = list(cex = 2),
                     key = list(rect = list("col" = c("#F708F7", "#FB84FB", "#F70606", "#FDC1C1")), 
                                text = list("samples" = c("Basal.parous.C", "Basal.parous.D", "Basal.virgin.C", "Basal.virgin.D")), 
                                columns=2, space = 'top'))
plotPC2PC3 <- xyplot(PC2 ~ PC3, data = pcaDF.basal.exp3, pch = 16, cex = 2,  col = as.vector(pcaDF.basal.exp3$colors), scales=list(tck=c(1,1), x=list(cex=2), y=list(cex=2)), xlab = list(cex = 2), ylab = list(cex = 2), 
                     key = list(rect = list("col" = c("#F708F7", "#FB84FB", "#F70606", "#FDC1C1")), 
                                text = list("samples" = c("Basal.parous.C", "Basal.parous.D", "Basal.virgin.C", "Basal.virgin.D")), 
                                columns=2, space = 'top'))
pdf("../../plots/bulk/exp3/basal-plotPC1PC2.pdf", width = 10)
try(plotPC1PC2, silent = TRUE)
dev.off()
pdf("../../plots/bulk/exp3/basal-plotPC1PC3.pdf", width = 10)
try(plotPC1PC3, silent = TRUE)
dev.off()
pdf("../../plots/bulk/exp3/basal-plotPC2PC3.pdf", width = 10)
try(plotPC2PC3, silent = TRUE)
dev.off()
pdf("../../plots/bulk/exp3/basal-PCcomponents.pdf", width = 10)
plot(x, pcaDF.basal.exp3$importance1, type="b", pch=18, col="red", xlab="Principal Component", ylab="", ylim=c(0,1))
lines(x, pcaDF.basal.exp3$importance2, type="b", pch=18, col="blue", lty=2)
legend((length(x)-10), 0.5, legend=c("individual","cumulative"), col=c("red","blue"), lty=1:2, cex=0.8)
dev.off()

plot3d(pcaDF.basal.exp3$`PC1`, pcaDF.basal.exp3$`PC2`, pcaDF.basal.exp3$`PC3`, col = pcaDF.basal.exp3$colors, size = 2, type = 's', xlab = 'PC1', ylab = 'PC2', zlab = 'PC3')
q <- rglwidget()
htmlwidgets::saveWidget(q, file = '../../plots/bulk/exp3/basal-all_samples_pca_3D.html')

text3d(pcaDF.basal.exp3$`PC1`, pcaDF.basal.exp3$`PC2`, pcaDF.basal.exp3$`PC3`, pcaDF.basal.exp3$Sample, adj = c(-0.1,-0.8))
q <- rglwidget()
htmlwidgets::saveWidget(q, file = '../../plots/bulk/exp3/basal-all_samples_pca_3D_with_label.html')



## exp3 luminal
pca.luminal.exp3 <- prcomp(t(assay(vsd.luminal.exp3)))
importanceIndividual.luminal.exp3 <- summary(pca.luminal.exp3)$importance[2,]
importanceCum.luminal.exp3 <- summary(pca.luminal.exp3)$importance[3,]
pcaDF.luminal.exp3 <- as.data.frame(pca.luminal.exp3$x)
pcaDF.luminal.exp3 <- pcaDF.luminal.exp3[,1:3]
pcaDF.luminal.exp3$Sample <- colnames(vsd.luminal.exp3)
pcaDF.luminal.exp3$colors <- vsd.luminal.exp3$colors
pcaDF.luminal.exp3 <- pcaDF.luminal.exp3[sort(rownames(pcaDF.luminal.exp3)),]
pcaDF.luminal.exp3$importance1 <- importanceIndividual.luminal.exp3
pcaDF.luminal.exp3$importance2 <- importanceCum.luminal.exp3

x <- c(1:length(pcaDF.luminal.exp3$importance1))

colorList <- list(col = as.character(pcaDF.luminal.exp3$colors))
textList <- list(as.character(pcaDF.luminal.exp3$Sample))
plotPC1PC2 <- xyplot(PC1 ~ PC2, data = pcaDF.luminal.exp3, pch = 16, cex = 2,  col = as.vector(pcaDF.luminal.exp3$colors), scales=list(tck=c(1,1), x=list(cex=2), y=list(cex=2)), xlab = list(cex = 2), ylab = list(cex = 2), key = list(rect = colorList, text = textList, columns=2, space = 'top'))
plotPC1PC2 <- xyplot(PC1 ~ PC2, data = pcaDF.luminal.exp3, pch = 16, cex = 2,  
                     col = as.vector(pcaDF.luminal.exp3$colors), scales=list(tck=c(1,1), x=list(cex=2), y=list(cex=2)), 
                     xlab = list(cex = 2), ylab = list(cex = 2),
                     key = list(rect = list("col" = c("#08D5F7", "#84EAFB", "#0E08F9", "#C3C2FE")), 
                                text = list("samples" = c("Luminal.parous.C", "Luminal.parous.D","Luminal.virgin.C", "Luminal.virgin.D")), 
                                columns=2, space = 'top'))
plotPC1PC3 <- xyplot(PC1 ~ PC3, data = pcaDF.luminal.exp3, pch = 16, cex = 2,  col = as.vector(pcaDF.luminal.exp3$colors), scales=list(tck=c(1,1), x=list(cex=2), y=list(cex=2)), xlab = list(cex = 2), ylab = list(cex = 2),
                     key = list(rect = list("col" = c("#08D5F7", "#84EAFB", "#0E08F9", "#C3C2FE")), 
                                text = list("samples" = c("Luminal.parous.C", "Luminal.parous.D","Luminal.virgin.C", "Luminal.virgin.D")), 
                                columns=2, space = 'top'))
plotPC2PC3 <- xyplot(PC2 ~ PC3, data = pcaDF.luminal.exp3, pch = 16, cex = 2,  col = as.vector(pcaDF.luminal.exp3$colors), scales=list(tck=c(1,1), x=list(cex=2), y=list(cex=2)), xlab = list(cex = 2), ylab = list(cex = 2), 
                     key = list(rect = list("col" = c("#08D5F7", "#84EAFB", "#0E08F9", "#C3C2FE")), 
                                text = list("samples" = c("Luminal.parous.C", "Luminal.parous.D","Luminal.virgin.C", "Luminal.virgin.D")), 
                                columns=2, space = 'top'))

pdf("../../plots/bulk/exp3/luminal-plotPC1PC2.pdf", width = 10)
try(plotPC1PC2, silent = TRUE)
dev.off()
pdf("../../plots/bulk/exp3/luminal-plotPC1PC3.pdf", width = 10)
try(plotPC1PC3, silent = TRUE)
dev.off()
pdf("../../plots/bulk/exp3/luminal-plotPC2PC3.pdf", width = 10)
try(plotPC2PC3, silent = TRUE)
dev.off()
pdf("../../plots/bulk/exp3/luminal-PCcomponents.pdf", width = 10)
plot(x, pcaDF.luminal.exp3$importance1, type="b", pch=18, col="red", xlab="Principal Component", ylab="", ylim=c(0,1))
lines(x, pcaDF.luminal.exp3$importance2, type="b", pch=18, col="blue", lty=2)
legend((length(x)-10), 0.5, legend=c("individual","cumulative"), col=c("red","blue"), lty=1:2, cex=0.8)
dev.off()

plot3d(pcaDF.luminal.exp3$`PC1`, pcaDF.luminal.exp3$`PC2`, pcaDF.luminal.exp3$`PC3`, col = pcaDF.luminal.exp3$colors, size = 2, type = 's', xlab = 'PC1', ylab = 'PC2', zlab = 'PC3')
q <- rglwidget()
htmlwidgets::saveWidget(q, file = '../../plots/bulk/exp3/luminal-all_samples_pca_3D.html')

text3d(pcaDF.luminal.exp3$`PC1`, pcaDF.luminal.exp3$`PC2`, pcaDF.luminal.exp3$`PC3`, pcaDF.luminal.exp3$Sample, adj = c(-0.1,-0.8))
q <- rglwidget()
htmlwidgets::saveWidget(q, file = '../../plots/bulk/exp3//luminal-all_samples_pca_3D_with_label.html')



############################################################################
############################################################################
###                                                                      ###
###                                 EXP7                                 ###
###                                                                      ###
############################################################################
############################################################################
pca.all.exp7 <- prcomp(t(assay(vsd.all.exp7)))
importanceIndividual.all.exp7 <- summary(pca.all.exp7)$importance[2,]
importanceCum.all.exp7 <- summary(pca.all.exp7)$importance[3,]
pcaDF.all.exp7 <- as.data.frame(pca.all.exp7$x)
pcaDF.all.exp7 <- pcaDF.all.exp7[,1:3]
pcaDF.all.exp7$Sample <- colnames(vsd.all.exp7)
pcaDF.all.exp7$colors <- vsd.all.exp7$colors
pcaDF.all.exp7 <- pcaDF.all.exp7[sort(rownames(pcaDF.all.exp7)),]
pcaDF.all.exp7$importance1 <- importanceIndividual.all.exp7
pcaDF.all.exp7$importance2 <- importanceCum.all.exp7

x <- c(1:length(pcaDF.all.exp7$importance1))

colorList <- list(col = as.character(pcaDF.all.exp7$colors))
textList <- list(as.character(pcaDF.all.exp7$Sample))
plotPC1PC2 <- xyplot(PC1 ~ PC2, data = pcaDF.all.exp7, pch = 16, cex = 2,  
                     col = as.vector(pcaDF.all.exp7$colors), scales=list(tck=c(1,1), x=list(cex=2), y=list(cex=2)), 
                     xlab = list(cex = 2), ylab = list(cex = 2),
                     key = list(rect = list("col" = c("#0E08F9", "#C3C2FE", "#635FFB", "#84EAFB", "#08D5F7")), 
                                text = list("samples" = c("Luminal.mE2.C.TF", "Luminal.mE2.D.TF", "Luminal.pE2.C.TB",
                                                          "Luminal.pE2.D.TB", "Luminal.pE2.D.TF")), 
                                columns=2, space = 'top'))
plotPC1PC3 <- xyplot(PC1 ~ PC3, data = pcaDF.all.exp7, pch = 16, cex = 2,  col = as.vector(pcaDF.all.exp7$colors), scales=list(tck=c(1,1), x=list(cex=2), y=list(cex=2)), xlab = list(cex = 2), ylab = list(cex = 2),
                     key = list(rect = list("col" = c("#0E08F9", "#C3C2FE", "#635FFB", "#84EAFB", "#08D5F7")), 
                                text = list("samples" = c("Luminal.mE2.C.TF", "Luminal.mE2.D.TF", "Luminal.pE2.C.TB",
                                                          "Luminal.pE2.D.TB", "Luminal.pE2.D.TF")), 
                                columns=2, space = 'top'))
plotPC2PC3 <- xyplot(PC2 ~ PC3, data = pcaDF.all.exp7, pch = 16, cex = 2,  col = as.vector(pcaDF.all.exp7$colors), scales=list(tck=c(1,1), x=list(cex=2), y=list(cex=2)), xlab = list(cex = 2), ylab = list(cex = 2), 
                     key = list(rect = list("col" = c("#0E08F9", "#C3C2FE", "#635FFB", "#84EAFB", "#08D5F7")), 
                                text = list("samples" = c("Luminal.mE2.C.TF", "Luminal.mE2.D.TF", "Luminal.pE2.C.TB",
                                                          "Luminal.pE2.D.TB", "Luminal.pE2.D.TF")), 
                                columns=2, space = 'top'))
pdf("../../plots/bulk/exp7/plotPC1PC2.pdf", width = 10)
try(plotPC1PC2, silent = TRUE)
dev.off()
pdf("../../plots/bulk/exp7/plotPC1PC3.pdf", width = 10)
try(plotPC1PC3, silent = TRUE)
dev.off()
pdf("../../plots/bulk/exp7/plotPC2PC3.pdf", width = 10)
try(plotPC2PC3, silent = TRUE)
dev.off()
pdf("../../plots/bulk/plots/exp7/PCcomponents.pdf", width = 10)
plot(x, pcaDF.all.exp7$importance1, type="b", pch=18, col="red", xlab="Principal Component", ylab="", ylim=c(0,1))
lines(x, pcaDF.all.exp7$importance2, type="b", pch=18, col="blue", lty=2)
legend((length(x)-10), 0.5, legend=c("individual","cumulative"), col=c("red","blue"), lty=1:2, cex=0.8)
dev.off()


plot3d(pcaDF.all.exp7$`PC1`, pcaDF.all.exp7$`PC2`, pcaDF.all.exp7$`PC3`, col = pcaDF.all.exp7$colors, size = 2, type = 's', xlab = 'PC1', ylab = 'PC2', zlab = 'PC3')
q <- rglwidget()
htmlwidgets::saveWidget(q, file = '../../plots/bulk/exp7/all_samples_pca_3D.html')

text3d(pcaDF.all.exp7$`PC1`, pcaDF.all.exp7$`PC2`, pcaDF.all.exp7$`PC3`, pcaDF.all.exp7$Sample, adj = c(-0.1,-0.8))
q <- rglwidget()
htmlwidgets::saveWidget(q, file = '../../plots/bulk/exp7/all_samples_pca_3D_with_label.html')





############################################################################
############################################################################
###                                                                      ###
###                            DO COMPARISONS                            ###
###                                                                      ###
############################################################################
############################################################################
dds.exp3 <- list()
res.virginD.vs.virginC <- list()
res.parousD.vs.parousC <- list()
res.shrink.virginD.vs.virginC <- list()
res.shrink.parousD.vs.parousC <- list()
i <- 0
for (c in c("Basal", "Luminal", "CD45")){
  i <- i+1
  dds.exp3[[i]] <- DESeqDataSetFromMatrix(countData = data.exp3[[i]], colData = phenos.exp3[[i]], design=~merged)
  dds.exp3[[i]] <- DESeq(dds.exp3[[i]])
  
  res.virginD.vs.virginC[[i]] <- results(dds.exp3[[i]], contrast = c("merged", "virgin.D", "virgin.C"))
  res.parousD.vs.parousC[[i]] <- results(dds.exp3[[i]], contrast = c("merged", "parous.D", "parous.C"))
  
  res.shrink.virginD.vs.virginC[[i]] <- lfcShrink(dds.exp3[[i]], contrast = c("merged", "virgin.D", "virgin.C"), type="ashr")
  res.shrink.parousD.vs.parousC[[i]] <- lfcShrink(dds.exp3[[i]], contrast = c("merged", "parous.D", "parous.C"), type="ashr")
}
names(dds.exp3) <- c("Basal", "Luminal", "CD45")
names(res.virginD.vs.virginC) <- c("Basal", "Luminal", "CD45")
names(res.parousD.vs.parousC) <- c("Basal", "Luminal", "CD45")
names(res.shrink.virginD.vs.virginC) <- c("Basal", "Luminal", "CD45")
names(res.shrink.parousD.vs.parousC) <- c("Basal", "Luminal", "CD45")


dds.exp7 <- list()
res.mE2DTF.vs.mE2CTF <- list()
res.pE2DTB.vs.pE2CTB <- list()
res.pE2DTF.vs.pE2CTB <- list()
res.pE2DTB.vs.pE2DTF <- list()

i <- 0
for (c in c("Luminal", "CD45")){
  i <- i+1
  dds.exp7[[i]] <- DESeqDataSetFromMatrix(countData = data.exp7[[i]], colData = phenos.exp7[[i]], design=~merged)
  dds.exp7[[i]] <- DESeq(dds.exp7[[i]])
  
  res.mE2DTF.vs.mE2CTF[[i]] <- results(dds.exp7[[i]], contrast = c("merged", "mE2.D.TF", "mE2.C.TF"))
  res.pE2DTB.vs.pE2CTB[[i]] <- results(dds.exp7[[i]], contrast = c("merged", "pE2.D.TB", "pE2.C.TB"))
  res.pE2DTF.vs.pE2CTB[[i]] <- results(dds.exp7[[i]], contrast = c("merged", "pE2.D.TF", "pE2.C.TB"))
  res.pE2DTB.vs.pE2DTF[[i]] <- results(dds.exp7[[i]], contrast = c("merged", "pE2.D.TB", "pE2.D.TF"))
}
names(dds.exp7) <- c("Luminal", "CD45")
names(res.mE2DTF.vs.mE2CTF) <- c("Luminal", "CD45")
names(res.pE2DTB.vs.pE2CTB) <- c("Luminal", "CD45")
names(res.pE2DTF.vs.pE2CTB) <- c("Luminal", "CD45")
names(res.pE2DTB.vs.pE2DTF) <- c("Luminal", "CD45")


####### filter and write tables
p.adj.thresh <- 0.05
lfc.thresh <- 1


## exp3
# virginD.vs.virginC
filter.res.virginD.vs.virginC <- lapply(res.virginD.vs.virginC, function(x){
  x <- x[order(x$padj),]
  x <- x[intersect(which(x$padj <= p.adj.thresh), which(abs(x$log2FoldChange) >= lfc.thresh)) ,]
})
lapply(filter.res.virginD.vs.virginC, dim)
lapply(filter.res.virginD.vs.virginC, function(x){
  return(c("pos" = length(which(x$log2FoldChange >= 0)), "neg" = length(which(x$log2FoldChange < 0))))})

list.to.plot <- filter.res.virginD.vs.virginC
name.file <- paste0(strsplit(deparse(substitute(filter.res.virginD.vs.virginC)),split = ".", fixed = TRUE)[[1]][3:5], collapse = ".")
for (i in 1:length(list.to.plot)){
  n <- names(list.to.plot)[i]
  write.table(list.to.plot[[i]], sep = "\t", quote = FALSE, file = paste0("~/Downloads/exp3/",n,".", name.file, ".txt"))
}


# parousD.vs.parousC
filter.res.parousD.vs.parousC <- lapply(res.parousD.vs.parousC, function(x){
  x <- x[order(x$padj),]
  x <- x[intersect(which(x$padj <= p.adj.thresh), which(abs(x$log2FoldChange) >= lfc.thresh)) ,]
})
lapply(filter.res.parousD.vs.parousC, dim)
lapply(filter.res.parousD.vs.parousC, function(x){
  return(c("pos" = length(which(x$log2FoldChange >= 0)), "neg" = length(which(x$log2FoldChange < 0))))})

list.to.plot <- filter.res.parousD.vs.parousC
name.file <- paste0(strsplit(deparse(substitute(filter.res.parousD.vs.parousC)),split = ".", fixed = TRUE)[[1]][3:5], collapse = ".")
for (i in 1:length(list.to.plot)){
  n <- names(list.to.plot)[i]
  write.table(list.to.plot[[i]], sep = "\t", quote = FALSE, file = paste0("../../tables/bulk-DEGs/exp3/",n,".", name.file, ".txt"))
}


# shrink virginD.vs.virginC
filter.res.shrink.virginD.vs.virginC <- lapply(res.shrink.virginD.vs.virginC, function(x){
  x <- x[order(x$padj),]
  x <- x[intersect(which(x$padj <= p.adj.thresh), which(abs(x$log2FoldChange) >= lfc.thresh)) ,]
})
lapply(filter.res.shrink.virginD.vs.virginC, dim)
lapply(filter.res.shrink.virginD.vs.virginC, function(x){
  return(c("pos" = length(which(x$log2FoldChange >= 0)), "neg" = length(which(x$log2FoldChange < 0))))})

list.to.plot <- filter.res.shrink.virginD.vs.virginC
name.file <- paste0(strsplit(deparse(substitute(filter.res.shrink.virginD.vs.virginC)),split = ".", fixed = TRUE)[[1]][3:6], collapse = ".")
for (i in 1:length(list.to.plot)){
  n <- names(list.to.plot)[i]
  write.table(list.to.plot[[i]], sep = "\t", quote = FALSE, file = paste0("../../tables/bulk-DEGs/exp3/",n,".", name.file, ".txt"))
}


# shrink parousD.vs.parousC
filter.res.shrink.parousD.vs.parousC <- lapply(res.shrink.parousD.vs.parousC, function(x){
  x <- x[order(x$padj),]
  x <- x[intersect(which(x$padj <= p.adj.thresh), which(abs(x$log2FoldChange) >= lfc.thresh)) ,]
})
lapply(filter.res.shrink.parousD.vs.parousC, dim)
lapply(filter.res.shrink.parousD.vs.parousC, function(x){
  return(c("pos" = length(which(x$log2FoldChange >= 0)), "neg" = length(which(x$log2FoldChange < 0))))})


list.to.plot <- filter.res.shrink.parousD.vs.parousC
name.file <- paste0(strsplit(deparse(substitute(filter.res.shrink.parousD.vs.parousC)),split = ".", fixed = TRUE)[[1]][3:6], collapse = ".")
for (i in 1:length(list.to.plot)){
  n <- names(list.to.plot)[i]
  write.table(list.to.plot[[i]], sep = "\t", quote = FALSE, file = paste0("../../tables/bulk-DEGs/exp3/",n,".", name.file, ".txt"))
}



## exp7
# mE2DTF.vs.mE2CTF
filter.res.mE2DTF.vs.mE2CTF <- lapply(res.mE2DTF.vs.mE2CTF, function(x){
  x <- x[order(x$padj),]
  x <- x[intersect(which(x$padj <= p.adj.thresh), which(abs(x$log2FoldChange) >= lfc.thresh)) ,]
})
lapply(filter.res.mE2DTF.vs.mE2CTF, dim)
lapply(filter.res.mE2DTF.vs.mE2CTF, function(x){
  return(c("pos" = length(which(x$log2FoldChange >= 0)), "neg" = length(which(x$log2FoldChange < 0))))})

list.to.plot <- filter.res.mE2DTF.vs.mE2CTF
name.file <- paste0(strsplit(deparse(substitute(filter.res.mE2DTF.vs.mE2CTF)),split = ".", fixed = TRUE)[[1]][3:5], collapse = ".")
for (i in 1:length(list.to.plot)){
  n <- names(list.to.plot)[i]
  write.table(list.to.plot[[i]], sep = "\t", quote = FALSE, file = paste0("../../tables/bulk-DEGs/exp7/", n,".", name.file, ".txt"))
}


# pE2DTB.vs.pE2CTB
filter.res.pE2DTB.vs.pE2CTB <- lapply(res.pE2DTB.vs.pE2CTB, function(x){
  x <- x[order(x$padj),]
  x <- x[intersect(which(x$padj <= p.adj.thresh), which(abs(x$log2FoldChange) >= lfc.thresh)) ,]
})
lapply(filter.res.pE2DTB.vs.pE2CTB, dim)
lapply(filter.res.pE2DTB.vs.pE2CTB, function(x){
  return(c("pos" = length(which(x$log2FoldChange >= 0)), "neg" = length(which(x$log2FoldChange < 0))))})

list.to.plot <- filter.res.pE2DTB.vs.pE2CTB
name.file <- paste0(strsplit(deparse(substitute(filter.res.pE2DTB.vs.pE2CTB)),split = ".", fixed = TRUE)[[1]][3:5], collapse = ".")
for (i in 1:length(list.to.plot)){
  n <- names(list.to.plot)[i]
  write.table(list.to.plot[[i]], sep = "\t", quote = FALSE, file = paste0("../../tables/bulk-DEGs/exp7/",n,".", name.file, ".txt"))
}


# pE2DTF.vs.pE2CTB
filter.res.pE2DTF.vs.pE2CTB <- lapply(res.pE2DTF.vs.pE2CTB, function(x){
  x <- x[order(x$padj),]
  x <- x[intersect(which(x$padj <= p.adj.thresh), which(abs(x$log2FoldChange) >= lfc.thresh)) ,]
})
lapply(filter.res.pE2DTF.vs.pE2CTB, dim)
lapply(filter.res.pE2DTF.vs.pE2CTB, function(x){
  return(c("pos" = length(which(x$log2FoldChange >= 0)), "neg" = length(which(x$log2FoldChange < 0))))})

list.to.plot <- filter.res.pE2DTF.vs.pE2CTB
name.file <- paste0(strsplit(deparse(substitute(filter.res.pE2DTF.vs.pE2CTB)),split = ".", fixed = TRUE)[[1]][3:5], collapse = ".")
for (i in 1:length(list.to.plot)){
  n <- names(list.to.plot)[i]
  write.table(list.to.plot[[i]], sep = "\t", quote = FALSE, file = paste0("../../tables/bulk-DEGs/exp7/",n,".", name.file, ".txt"))
}


# pE2DTB.vs.pE2DTF
filter.res.pE2DTB.vs.pE2DTF <- lapply(res.pE2DTB.vs.pE2DTF, function(x){
  x <- x[order(x$padj),]
  x <- x[intersect(which(x$padj <= p.adj.thresh), which(abs(x$log2FoldChange) >= lfc.thresh)) ,]
})
lapply(filter.res.pE2DTB.vs.pE2DTF, dim)
lapply(filter.res.pE2DTB.vs.pE2DTF, function(x){
  return(c("pos" = length(which(x$log2FoldChange >= 0)), "neg" = length(which(x$log2FoldChange < 0))))})

list.to.plot <- filter.res.pE2DTB.vs.pE2DTF
name.file <- paste0(strsplit(deparse(substitute(filter.res.pE2DTB.vs.pE2DTF)),split = ".", fixed = TRUE)[[1]][3:5], collapse = ".")
for (i in 1:length(list.to.plot)){
  n <- names(list.to.plot)[i]
  write.table(list.to.plot[[i]], sep = "\t", quote = FALSE, file = paste0("../../tables/bulk-DEGs/exp7/",n,".", name.file, ".txt"))
}
