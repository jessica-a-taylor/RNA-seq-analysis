source("Functions/loadLibraries.R")
loadLibraries()

DEGs <- read.csv("Results/Col-0/DEGs.csv")[,1] # 9084 DEGs
BjornsonData <- as.data.frame(read.csv("Bjornson et al. data analysis/High_conf_DEGs.csv")) # 6140 high-confidence DEGs

# 5062 DEGs in common

# Identify DEGs only detected in the analysis of Bjornson et al. data
Bjornson_only_DEGs <- as.data.frame(gff_genes[which(gff_genes$gene_id %in% BjornsonData$Gene[which(BjornsonData$Gene %!in% DEGs)]),]) %>%
  distinct(gene_id, .keep_all = TRUE)

bedFile <- GRanges(seqnames = Bjornson_only_DEGs$seqnames,
                   IRanges(start = Bjornson_only_DEGs$start, end = Bjornson_only_DEGs$end, width = Bjornson_only_DEGs$width),
                   strand = Bjornson_only_DEGs$strand)

rtracklayer::export.bed(bedFile, "Results/Bjornson-only-genes.bed")

# Identify DEGs only detected in the analysis of my RNA-seq data
Taylor_only_DEGs <- as.data.frame(gff_genes[which(gff_genes$gene_id %in% DEGs[which(DEGs %!in% BjornsonData$Gene)]),]) %>%
  distinct(gene_id, .keep_all = TRUE)

bedFile <- GRanges(seqnames = Taylor_only_DEGs$seqnames,
                   IRanges(start = Taylor_only_DEGs$start, end = Taylor_only_DEGs$end, width = Taylor_only_DEGs$width),
                   strand = Taylor_only_DEGs$strand)

rtracklayer::export.bed(bedFile, "Results/Taylor-only-genes.bed")