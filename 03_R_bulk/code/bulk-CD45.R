## RNASeq analysis of the CD45-sorted (immune) bulk samples
## Author: Simona Cristea
## Date: October 2022
## Publication: Breast cancer prevention by short-term inhibition of TGFβ signaling, Nat Comm 2022


library(biomaRt)
library(pheatmap)
library(DESeq2)
library(lattice)
library(tidyverse)



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
phenos <- phenos %>% dplyr::filter(note %in% c("","normal")) %>% dplyr::select(-note) %>% dplyr::filter(celltype == "CD45")



############################################################################
############################################################################
###                                                                      ###
###                         READ IN IMMUNE GENES                         ###
###                                                                      ###
############################################################################
############################################################################
orth <- read.table("../../tables/gene-lists/LM22_to_rnorvegicus.txt",sep = '\t', header = T)
orth <- orth$Gene.symbol
orth <- toupper(orth)



###########################################################################
###########################################################################
###                                                                     ###
###               SEPARATE PHENOS AND DATA PER EXPERIMENT               ###
###                                                                     ###
###########################################################################
###########################################################################
phenos.list <- list()
data.list <- list()
dds.list <- list()
vsd.list <- list()
normalizedTableVSD.list <- list()
normalizedCountsTable.list <- list()
normalizedTableVSD_2.list <- list()
anno.list <- list()
ph <- list()

i <- 0
for (c in c("exp3", "exp7")){
  i <- i+1
  phenos.list[[i]] <- phenos %>% filter(experiment == c)
  phenos.list[[i]]$merged <- paste(phenos.list[[i]]$condition, phenos.list[[i]]$treatment, sep = ".")
  phenos.list[[i]]$merged <- as.factor(phenos.list[[i]]$merged)

  data.list[[i]] <- data[, match(phenos.list[[i]]$Sample, colnames(data))]
  data.list[[i]] <- data.list[[i]][which(rowSums(data.list[[i]]) > 10),]
  
  dds.list[[i]] <- DESeqDataSetFromMatrix(countData = data.list[[i]], colData = phenos.list[[i]], design = ~merged)
  dds.list[[i]] <- estimateSizeFactors(dds.list[[i]])
  vsd.list[[i]] <- varianceStabilizingTransformation(dds.list[[i]])
  normalizedTableVSD.list[[i]] <- assay(vsd.list[[i]])
  normalizedTableVSD.list[[i]] <- normalizedTableVSD.list[[i]][,sort(colnames(normalizedTableVSD.list[[i]]))]
  normalizedCountsTable.list[[i]] <- as.data.frame(counts(dds.list[[i]], normalized = TRUE))
  
  write.table(normalizedCountsTable.list[[i]], 
              file = paste0("../../data/bulk-outputs/CD45_", c, "_size_normalized.tsv"), sep = '\t', row.names = T, col.names = T, quote = F)
  
  normalizedCountsTable.list[[i]] <- normalizedCountsTable.list[[i]][rownames(normalizedCountsTable.list[[i]]) %in% orth,]
  write.table(normalizedCountsTable.list[[i]], 
              file = paste0("../../data/bulk-outputs/CD45_", c, "_immune_expression_size_normalized.tsv"), sep = '\t', row.names = T, col.names = T, quote = F)
  
  normalizedTableVSD.list[[i]] <- normalizedTableVSD.list[[i]][rownames(normalizedTableVSD.list[[i]]) %in% orth,]
  normalizedTableVSD_2.list[[i]]  <- normalizedTableVSD.list[[i]] - rowMeans(normalizedTableVSD.list[[i]])
  anno.list[[i]] <- data.frame(Merged_factors = phenos.list[[i]]$merged)
  #rownames(anno.list[[i]]) = c(rownames(phenos.list[[i]]))
  rownames(anno.list[[i]]) = phenos.list[[i]]$Sample
  if (c == "exp3"){
    Merged_factors <- c("#F708F7", "#FB84FB", "#F70606", "#FDC1C1")
    names(Merged_factors) <- c("parous.C", "parous.D", "virgin.C", "virgin.D")
    anno_colors <- list(Merged_factors = Merged_factors)
  }
  
  if (c == "exp7"){
    Merged_factors <- c("#0E08F9", "#C3C2FE", "#635FFB", "#FDC1C1", "#84EAFB")
    names(Merged_factors) <- c("mE2.C.TF", "mE2.D.TF", "pE2.C.TB", "pE2.D.TB", "pE2.D.TF")
    anno_colors <- list(Merged_factors = Merged_factors)
  }

  ph[[i]] <- pheatmap(normalizedTableVSD_2.list[[i]], width = 9, height = 11, show_rownames = FALSE, 
                      color = colorRampPalette(c("blue", "white", "red"))( 100 ), 
                      annotation_col = anno.list[[i]], annotation_colors = anno_colors, 
                      main = "relative VST counts of CD45 rat samples clutered with immune related genes", 
                      filename = paste0("../../plots/bulk/", c, "/CD45_immune_expression.pdf"))
}
names(phenos.list) <- c("exp3", "exp7")
names(data.list) <- c("exp3", "exp7")

