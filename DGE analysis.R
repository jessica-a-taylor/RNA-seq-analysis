source("Functions/loadLibraries.R")
loadLibraries()

gff_genes <- rtracklayer::import("ColCEN_GENES.gff3")[,c(2,6)]
gff_genes <- gff_genes[gff_genes$type=="gene",]
gff_genes <- gff_genes[-c(which(is.na(gff_genes$gene_id))),]

allGenotypes <- c("Col", "sdg2", "atx1", "jmj14")

# Import raw read counts.
rawCounts <- data.frame(Gene = read.table("Counts_data/Col/Col_F0_1_readCounts.txt")[,1])

for (genotype in allGenotypes) {
  for (file in list.files(paste0("Counts_data/", genotype))) {
    rawCounts <- cbind(rawCounts, read.table(paste0("Counts_data/", genotype, "/",file))[,2])
  }
}

# Set row names to gene IDs
rownames(rawCounts) <- rawCounts$Gene
rawCounts <- rawCounts[,-1]

# Set column names to sample IDs
colnames(rawCounts) <- c(str_match(list.files("Counts_data/Col"), "^(.*)_readCounts.txt")[,2],
                         str_match(list.files("Counts_data/sdg2"), "^(.*)_readCounts.txt")[,2],
                         str_match(list.files("Counts_data/atx1"), "^(.*)_readCounts.txt")[,2],
                         str_match(list.files("Counts_data/jmj14"), "^(.*)_readCounts.txt")[,2])

# Create metadata file.
colData <- data.frame(Sample = factor(colnames(rawCounts)),
                      Genotype = factor(str_match(colnames(rawCounts),"^([A-Za-z0-9]+)_F[0-9]+_[0-9]+$")[,2],
                                        levels = c("Col", "sdg2", "atx1", "jmj14")),
                      Time = factor(str_match(colnames(rawCounts),"^[A-Za-z0-9]+_(F[0-9]+)_[0-9]+$")[,2],
                                    levels = c("F0", "F30", "F90", "F180")))

rownames(colData) <- colnames(rawCounts)

###############################################################################
# Identify DEGs relative to 0 min for each genotype.
# Construct DESeq dataset.
for (genotype in allGenotypes) {
  dds <- DESeqDataSetFromMatrix(countData = rawCounts[,which(str_detect(colnames(rawCounts),genotype)==TRUE)], 
                                colData = colData[which(colData$Genotype==genotype),], 
                                design = ~ Time)
  
  # Remove genes with < 10 reads.
  dds <- dds[rowSums(counts(dds)) >= 10,]
  
  # Run the DESeq.
  DDS <- DESeq(dds)

  for (time in c("F30", "F90", "F180")) {
  res <- results(DDS, name = paste0("Time_",time,"_vs_F0"), independentFiltering = FALSE)
  res <- res[,c(2,6)]
  colnames(res) <- c(paste0(time, "_FC"), paste0(time, "_padj"))
  
  write.csv(res, paste0("Results/", genotype, "/", time, "_vs_F0.csv"))
  }
}

###############################################################################
# Identify DEGs relative to WT at each time point.
# Construct DESeq dataset.
dds <- DESeqDataSetFromMatrix(countData = rawCounts, colData = colData, design = ~ Genotype+Time+Genotype:Time)

# Remove genes with < 10 reads.
dds <- dds[rowSums(counts(dds)) >= 10,]

# Run the DESeq.
DDS <- DESeq(dds)
normDDS <- counts(DDS, normalized = TRUE) # normalisation with respect to the sequencing depth

for (genotype in allGenotypes) {
  normDDS <- calculate_Zscore(normDDS, genotype)
}
write.csv(normDDS, "Results/Normalised_counts.csv")

# DEGs in Col-0.
DEGs <- c()

for (time in c("F30", "F90", "F180")) {
  res <- results(DDS, name = paste0("Time_",time,"_vs_F0"), independentFiltering = FALSE)
  res <- res[,c(2,6)]
  colnames(res) <- c(paste0(time, "_FC"), paste0(time, "_padj"))
  
  write.csv(res, paste0("Results/Col/", time, "_vs_F0.csv"))
  
  # Get list of DEGs with logFC < -1 or > 1 and padj <= .01
  DEGs <- append(DEGs, rownames(res[which(res[,1]>=1 & res[,2]<=.01 | res[,1]<=-1 & res[,2]<=.01),]))
}

DEGs <- unique(DEGs)
DEGs_normCounts <- normDDS[which(rownames(normDDS) %in% DEGs),]
write.csv(DEGs_normCounts[,which(str_detect(colnames(DEGs_normCounts), "Col")==TRUE)], 
          "Results/Col/DEGs.csv")

bedFile <- as.data.frame(gff_genes[which(gff_genes$gene_id %in% rownames(DEGs_normCounts)),]) %>%
  distinct(gene_id, .keep_all = TRUE)

bedFile <- GRanges(seqnames = bedFile$seqnames,
                   IRanges(start = bedFile$start, end = bedFile$end, width = bedFile$width),
                   strand = bedFile$strand)

rtracklayer::export.bed(bedFile, "Results/Col/DEGs.bed")
rm(bedFile)

# DEGs whose temporal response differs in the mutants.
for (genotype in allGenotypes[-1]) {
  DEGs <- c()
  
  for (time in c("F30", "F90", "F180")) {
    res <- results(DDS, name = paste0("Genotype",genotype,".Time",time), independentFiltering = FALSE)
    res <- res[,c(2,6)]
    colnames(res) <- c(paste0(time, "_FC"), paste0(time, "_padj"))
    
    write.csv(res, paste0("Results/", genotype, "/", time, "_vs_WT.csv"))
    
    # Get list of DEGs with logFC < -1 or > 1 and padj <= .05
    DEGs <- append(DEGs, rownames(res[which(res[,1]>=1 & res[,2]<=.05 | res[,1]<=-1 & res[,2]<=.05),]))
  }
  DEGs <- unique(DEGs)
  DEGs_vs_WT <- normDDS[which(rownames(normDDS) %in% DEGs),]
  write.csv(DEGs_vs_WT[,which(str_detect(colnames(DEGs_vs_WT), genotype)==TRUE)], 
            paste0("Results/", genotype,"/","DEGs_vs_WT.csv"))
  
  bedFile <- as.data.frame(gff_genes[which(gff_genes$gene_id %in% rownames(DEGs_vs_WT)),]) %>%
    distinct(gene_id, .keep_all = TRUE)
  
  bedFile <- GRanges(seqnames = bedFile$seqnames,
                     IRanges(start = bedFile$start, end = bedFile$end, width = bedFile$width),
                     strand = bedFile$strand)
  
  rtracklayer::export.bed(bedFile, paste0("Results/", genotype,"/","DEGs_vs_WT.bed"))
  rm(bedFile)
}

###############################################################################
# PCA
vsd <- vst(dds, blind=TRUE) # transformation of counts data
pcaData <- plotPCA(vsd, intgroup = c("Genotype", "Time"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

plot <- ggplot(pcaData, aes(PC1, PC2, color = Genotype, shape = Time)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw()

plot <- ggplot(pcaData, aes(PC1, PC2, color = Genotype)) +
  geom_point(size = 3) +
  facet_wrap(~ Time) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 16, colour = "black"),
        legend.text = element_text(size = 14, colour = "black"),
        legend.title = element_text(size = 16, colour = "black"),
        strip.text = element_text(size = 16, colour = "black"))

png("Figures/QC/PCA plot.png", width = 750, height = 700)
print(plot)
dev.off()

# Hierarchical clustering
vsd_mat <- assay(vsd)
vsd_cor <- cor(vsd_mat) 
plot <- pheatmap(vsd_cor, fontsize = 14)

png("Figures/QC/Hierarchical clustering heatmap.png", width = 1000, height = 1000)
print(plot)
dev.off()

