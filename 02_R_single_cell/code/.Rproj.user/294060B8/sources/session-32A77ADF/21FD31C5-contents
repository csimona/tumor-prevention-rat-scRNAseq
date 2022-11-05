## interactome analysis on sorted SD epithlial data
## Author: Simona Cristea
## Date: October 2022
## Publication: Breast cancer prevention by short-term inhibition of TGFβ signaling, Nat Comm 2022



############################################################################
############################################################################
###                                                                      ###
###                      READ IN DATA AND RUN LIANA                      ###
###                                                                      ###
############################################################################
############################################################################
library(Seurat)
library(Connectome)
require(liana)
require(tibble)
require(purrr)

rat.drug <- readRDS("../../data/rdata/rat-sd4-integrated-treatment.RDS")
rat.drug.connectome <- liana_wrap(rat.drug$drug, method = c('connectome'), resource = c('all'))
rat.control.connectome <- liana_wrap(rat.drug$control, method = c('connectome'), resource = c('all'))


#saveRDS(file = "../../data/rdata/rat.drug.connectome.RDS")
#saveRDS(file = "../../data/rdata/rat.control.connectome.RDS")


rat.seurat <- rat.drug$drug
rat.pathway.table <- rat.drug.connectome$connectome$OmniPath
type.rat <- "drug"


############################################################################
############################################################################
###                                                                      ###
###                               PATHWAYS                               ###
###                                                                      ###
############################################################################
############################################################################
pathways.path <- "../../tables/gene-lists/pathways/"
pathway.files <- list.files(path = pathways.path, pattern = "*.xlsx")
pathway.files
pathway.genes <- lapply(pathway.files, function(x){
  p <- read_excel(paste0(pathways.path, x))
  p <- p$To
  p <- toupper(p)
  p <- intersect(p, rownames(rat.seurat))
  return(p)
})
names(pathway.genes) <- unlist(sapply(pathway.files, function(x) {strsplit(x, split = "-pathway-genes.xlsx")[1]}))



############################################################################
############################################################################
###                                                                      ###
###                      PROCESS CONNECTOME OUTPUTS                      ###
###                                                                      ###
############################################################################
############################################################################
edges.tables <- lapply(pathway.genes, function(x, type = type.rat){
  outgoing <- rat.pathway.table %>% filter(ligand %in% x, receptor %in% x, weight_norm > 0, weight_sc > 0) %>%
    group_by(source, target) %>%
    tally() %>%
    group_by(source) %>%
    tally(n, name = "outgoing") %>%
    arrange(as.numeric(source))
  
  ingoing <- rat.pathway.table %>% filter(ligand %in% x, receptor %in% x, weight_norm > 0, weight_sc > 0) %>%
    group_by(source, target) %>% 
    tally() %>%
    group_by(target) %>%
    tally(n, name = "ingoing") %>%
    arrange(as.numeric(target))
  
  conns <- full_join(ingoing, outgoing, by = c("target" = "source"))
  conns <- conns %>% dplyr::rename(cluster = target)
  conns$celltype <- conns$cluster
  if (type.rat == "drug")
    conns$celltype <- dplyr::recode(conns$cluster,
                                    "0" = "LP",
                                    "1" = "basal",
                                    "2" = "basal",
                                    "3" = "basal",
                                    "4" = "LP",
                                    "5" = "LP",
                                    "6" = "ML",
                                    "7" = "ML",
                                    "8" = "LP",
                                    "9" = "basal",
                                    "10" = "SBC",
                                    "11" = "basal",
                                    "12" = "LP",
                                    "13" = "LP",
                                    "14" = "ML")
  if (type.rat == "control")
    conns$celltype <- dplyr::recode(conns$cluster,
                                    "0" = "LP",
                                    "1" = "basal",
                                    "2" = "ML",
                                    "3" = "LP",
                                    "4" = "basal",
                                    "5" = "basal",
                                    "6" = "LP",
                                    "7" = "LP",
                                    "8" = "basal",
                                    "9" = "SBC",
                                    "10" = "LP")
  
  conns$cluster <- as.factor(as.numeric(conns$cluster))
  
  conns$fraction.size.ingoing <- conns$ingoing/length(x)
  conns$fraction.size.outgoing <- conns$outgoing/length(x)
  
  conns$fraction.norm.ingoing <- conns$ingoing/max(conns$ingoing, na.rm = TRUE)
  conns$fraction.norm.outgoing <- conns$outgoing/max(conns$outgoing, na.rm = TRUE)
  
  return(conns)
})
names(edges.tables) <- names(pathway.genes)
edges.tables <- bind_rows(edges.tables, .id = "column_label")

tab.seurat.clusters <- table(rat.seurat$seurat_clusters)
edges.tables$norm.factor.per.cell <- tab.seurat.clusters[match(edges.tables$cluster, names(tab.seurat.clusters))]
edges.tables$norm.factor.per.cell <- as.integer(edges.tables$norm.factor.per.cell)
edges.tables$norm.ingoing.per.cell <- edges.tables$ingoing/edges.tables$norm.factor.per.cell
edges.tables$norm.outgoing.per.cell <- edges.tables$outgoing/edges.tables$norm.factor.per.cell

t <- table(rat.seurat$celltype.detailed)
names(t)[which(names(t) == "secretory")] <- "SBC"


pdf(paste0("../../plots/connectome/", type.rat, "-relative.pathways.pdf"), width = 18, height = 12)
edges.tables %>% pivot_longer(c(outgoing, ingoing), names_to = "connection.type", values_to = "size") %>% 
  group_by(column_label, celltype, connection.type) %>% summarise (sum.connections = sum(size, na.rm = TRUE)) %>%
  rowwise() %>% mutate("total.cells.per.celltype" = t[match(celltype, names(t))]) %>% 
  mutate("rel.connections.per.celltype" = sum.connections/total.cells.per.celltype) %>%
  ggplot(aes(x = connection.type, y = rel.connections.per.celltype, fill = celltype)) + geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~column_label) +
  theme_bw() +
  scale_size(range = c(1,8)) +
  theme(text = element_text(size = 16)) 
dev.off()


pdf(paste0("../../plots/connectome/", type.rat, "-relative.pdf"))
edges.tables %>% pivot_longer(c(outgoing, ingoing), names_to = "connection.type", values_to = "size") %>% 
  group_by(celltype, connection.type) %>% summarise (sum.connections = sum(size, na.rm = TRUE)) %>%
  rowwise() %>% mutate("total.cells.per.celltype" = t[match(celltype, names(t))]) %>% mutate("rel.connections.per.celltype" = sum.connections/total.cells.per.celltype) %>%
  ggplot(aes(x = connection.type, y = rel.connections.per.celltype, fill = celltype)) + geom_bar(stat = "identity", position = position_dodge()) +
  theme_bw() +
  scale_size(range = c(1,8)) +
  theme(text = element_text(size = 16))
dev.off()


pdf(paste0("../../plots/connectome/", type.rat, "-pathways-connections.pdf"), width = 18, height = 12)
edges.tables %>% pivot_longer(c(outgoing, ingoing), names_to = "connection.type", values_to = "size") %>% 
  ggplot + geom_point(aes(x = connection.type, y = cluster, size = size, color = celltype)) +
  facet_wrap(~column_label) +
  theme_bw() +
  scale_size(range = c(1,8)) +
  theme(text = element_text(size = 16))   
dev.off()


pdf(paste0("../../plots/connectome/", type.rat, "-pathways-normalized-01.pdf"), width = 18, height = 12)
edges.tables %>% pivot_longer(c(fraction.norm.ingoing, fraction.norm.outgoing), names_to = "connection.type", values_to = "fraction.to.max") %>% 
  ggplot + geom_point(aes(x = connection.type, y = cluster, size = fraction.to.max, color = celltype)) +
  facet_wrap(~column_label) +
  theme_bw() +
  scale_size(range = c(1,8)) +
  theme(text = element_text(size = 16))   
dev.off()


for (p in unique(edges.tables$column_label)){
  p.plot <- edges.tables %>% filter(column_label == p) %>% 
    pivot_longer(c(fraction.norm.ingoing, fraction.norm.outgoing), names_to = "connection.type", values_to = "fraction.to.max") %>%
    ggplot + geom_point(aes(x = connection.type, y = cluster, size = fraction.to.max, color = celltype)) +
    theme_bw() +
    scale_size(range = c(1,8)) +
    theme(text = element_text(size = 16)) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  pdf(paste0("../../plots/connectome/", type.rat, "-single-normalized/", p, ".pdf"))
  print(p.plot)
  dev.off()
}


for (p in unique(edges.tables$column_label)){
  p.plot <- edges.tables %>% filter(column_label == p) %>% 
    pivot_longer(c(outgoing, ingoing), names_to = "connection.type", values_to = "size") %>%
    ggplot + geom_point(aes(x = connection.type, y = cluster, size = size, color = celltype)) +
    theme_bw() +
    scale_size(range = c(1,8)) +
    theme(text = element_text(size = 16)) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  pdf(paste0("../../plots/connectome/", type.rat, "-single-absolute/", p, ".pdf"))
  print(p.plot)
  dev.off()
}

