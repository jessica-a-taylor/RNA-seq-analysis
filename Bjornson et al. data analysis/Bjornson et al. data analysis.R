source("Functions/loadLibraries.R")
loadLibraries()

gff_genes <- rtracklayer::import("ColCEN_GENES.gff3")[,c(2,6)]
gff_genes <- gff_genes[gff_genes$type=="gene",]
gff_genes <- gff_genes[-c(which(is.na(gff_genes$gene_id))),]

# Import raw read counts.
rawCounts <- read.csv("Bjornson et al. data analysis/PRJEB25079.csv")

# Set row names to gene IDs
rownames(rawCounts) <- rawCounts[,1]
rawCounts <- rawCounts[,-1]

# Filter for mock & flg22 treatment
rawCounts <- rawCounts[,which(str_detect(colnames(rawCounts), "flg")==TRUE |str_detect(colnames(rawCounts), "mock"))]

# Filter for 0, 30, 90 & 180 min time points
rawCounts <- rawCounts[,which(str_match(colnames(rawCounts),"^[A-Za-z]+([0-9]+)_[0-9]+$")[,2] %in% c("0", "30", "90", "180"))]

# Create metadata file.
colData <- data.frame(Sample = factor(colnames(rawCounts)),
                      Treatment = factor(str_match(colnames(rawCounts),"^([A-Za-z]+)[0-9]+_[0-9]+$")[,2],
                                        levels = c("mock", "flg")),
                      Time = factor(str_match(colnames(rawCounts),"^[A-Za-z]+([0-9]+)_[0-9]+$")[,2],
                                    levels = c("0", "30", "90", "180")))

rownames(colData) <- colnames(rawCounts)

################################################
# Construct DESeq dataset for comparison to mock.
dds <- DESeqDataSetFromMatrix(countData = rawCounts, colData = colData, design = ~ Treatment+Time+Treatment:Time)

# Remove genes with < 10 reads.
dds <- dds[rowSums(counts(dds)) >= 10,]

# Run the DESeq.
DDS <- DESeq(dds)
normDDS <- counts(DDS, normalized = TRUE) # normalisation with respect to the sequencing depth
write.csv(normDDS, "Bjornson et al. data analysis/Normalised_counts_including_mock.csv")

# Identify DEGs relative to mock.
DEGs <- c()

for (time in c("30", "90", "180")) {
  res <- results(DDS, name = paste0("Treatmentflg.Time",time), independentFiltering = FALSE)
  res <- res[,c(2,6)]
  colnames(res) <- c(paste0("F",time, "_FC"), paste0(time, "_padj"))
  
  write.csv(res, paste0("Bjornson et al. data analysis/F", time, "_vs_M", time, ".csv"))
  
  # Get list of DEGs with logFC < -1 or > 1 and padj <= .01
  DEGs <- append(DEGs, rownames(res[which(res[,1]>=1 & res[,2]<=.01 | res[,1]<=-1 & res[,2]<=.01),]))
}

DEGs <- unique(DEGs)
DEGs_normCounts <- normDDS[which(rownames(normDDS) %in% DEGs),]
write.csv(DEGs_normCounts, "Bjornson et al. data analysis/DEGs_vs_mock.csv")

################################################
# Construct DESeq dataset for comparison to flg0
dds <- DESeqDataSetFromMatrix(countData = rawCounts[,which(str_detect(colnames(rawCounts), "flg")==TRUE)], 
                              colData = colData[which(str_detect(rownames(colData), "flg")==TRUE),], 
                              design = ~ Time)

# Remove genes with < 10 reads.
dds <- dds[rowSums(counts(dds)) >= 10,]

# Run the DESeq.
DDS <- DESeq(dds)
normDDS <- counts(DDS, normalized = TRUE) # normalisation with respect to the sequencing depth
write.csv(normDDS, "Bjornson et al. data analysis/Normalised_counts_flg_only.csv")

# Identify DEGs relative to flg0
DEGs <- c()

for (time in c("30", "90", "180")) {
  res <- results(DDS, name = paste0("Time_",time, "_vs_0"), independentFiltering = FALSE)
  res <- res[,c(2,6)]
  colnames(res) <- c(paste0("F",time, "_FC"), paste0(time, "_padj"))
  
  write.csv(res, paste0("Bjornson et al. data analysis/Results/F", time, "_vs_F0.csv"))
  
  # Get list of DEGs with logFC < -1 or > 1 and padj <= .01
  DEGs <- append(DEGs, rownames(res[which(res[,1]>=1 & res[,2]<=.01 | res[,1]<=-1 & res[,2]<=.01),]))
}

DEGs <- unique(DEGs)
DEGs_normCounts <- normDDS[which(rownames(normDDS) %in% DEGs),]
write.csv(DEGs_normCounts, "Bjornson et al. data analysis/Results/DEGs_vs_F0.csv")

################################################
# Compare DEG lists
DEGs_vs_mock <- as.data.frame(read.csv("Bjornson et al. data analysis/Results/DEGs_vs_mock.csv")) # 6743 genes
DEGs_vs_F0 <- as.data.frame(read.csv("Bjornson et al. data analysis/Results/DEGs_vs_F0.csv")) # 7941 genes

# 6140 DEGs in common (high-confidence DEGs)
HC_DEGs <- DEGs_vs_F0[which(DEGs_vs_F0$X %in% DEGs_vs_mock$X),]
write.csv(HC_DEGs, "Bjornson et al. data analysis/High_conf_DEGs.csv")

bedFile <- as.data.frame(gff_genes[which(gff_genes$gene_id %in% HC_DEGs$X),]) %>%
  distinct(gene_id, .keep_all = TRUE)

bedFile <- GRanges(seqnames = bedFile$seqnames,
                   IRanges(start = bedFile$start, end = bedFile$end, width = bedFile$width),
                   strand = bedFile$strand)

rtracklayer::export.bed(bedFile, "Bjornson et al. data analysis/High_conf_DEGs.bed")
rm(bedFile)

################################################
normCounts <- as.data.frame(read.csv("Bjornson et al. data analysis/High_conf_DEGs.csv"))[,-1]
colnames(normCounts) <- c("gene_id", paste0("Col_", 
                                            str_match(colnames(normCounts)[-1], "^flg(.*)$")[,2]))
rownames(normCounts) <- normCounts$gene_id

# Remove mitochondrial and chloroplast genes.
normCounts <- normCounts[which(str_detect(normCounts$gene_id, "^(AT[1-5]G).*")==TRUE),]

# Filter for upregulated genes.
upreg_DEGs <- c()

for (time in c("F30", "F90", "F180")) {
  df <- read.csv(paste0("Bjornson et al. data analysis/Results/", time, "_vs_F0.csv"))
  df <- df[which(df[,2]>=1 & df[,3]<=.05),1]
  
  upreg_DEGs <- append(upreg_DEGs, df)
}

normCounts <- normCounts[which(normCounts$gene_id %in% upreg_DEGs),]

# Calculate Z-scores on high-confidence upregulated DEGs.
source("Functions/Calculate_Zscores.R")
normCounts <- calculate_Zscore(normCounts, "Col")

# Fuzzy k-means analysis on all DEGs.
Zscores <- as.matrix(normCounts[,which(str_detect(colnames(normCounts), "Zscore")==TRUE)])

set.seed(20)

fuzzyClusteringData <- ExpressionSet(assayData = as.matrix(Zscores))
m <- mestimate(fuzzyClusteringData)
fuzzClust <- mfuzz(fuzzyClusteringData, centers = 6, m = m)

unfilteredResults <- as.data.frame(fuzzyClusteringData@assayData$exprs)
unfilteredResults$cluster <- fuzzClust$cluster
write.csv(unfilteredResults, "Bjornson et al. data analysis/Results/Clustering_results_unfiltered.csv")

# Remove genes that cannot be assigned to a cluster with a membership value > 0.7.
goodGenes <- acore(fuzzyClusteringData, cl=fuzzClust, min.acore = .7)
genesToCluster <- c()
for (c in 1:length(goodGenes)) {
  genesToCluster <- append(genesToCluster, which(names(fuzzClust$cluster) %in% goodGenes[[c]]$NAME))
}

# Filter for high-confidence genes.
fuzzClust <- data.frame(Gene = names(fuzzClust$cluster[genesToCluster]),
                        Cluster = fuzzClust$cluster[genesToCluster])

fuzzClust <- fuzzClust[order(fuzzClust$Gene),]

fuzzyClusteringData <- as.data.frame(fuzzyClusteringData@assayData$exprs)
fuzzyClusteringData <- fuzzyClusteringData[which(rownames(fuzzyClusteringData) %in% fuzzClust$Gene),]
fuzzyClusteringData$cluster <- fuzzClust$Cluster

write.csv(fuzzyClusteringData, "Bjornson et al. data analysis/Results/Clustering_results_filtered.csv")

# Plot a boxplot of Z-scores at each time point in each cluster.
fuzzyClusteringData <- read.csv("Bjornson et al. data analysis/Results/Clustering_results_filtered.csv")[,-1]

clusterCentres <- data.frame()
for (clust in 1:6) {
  isolateCluster <- fuzzyClusteringData[fuzzyClusteringData$cluster==clust,]
  
  coreGenes <- c()
  for (col in 1:(ncol(isolateCluster)-1)) {
    coreGenes <- append(coreGenes, rownames(isolateCluster[which(isolateCluster[,col] >= quantile(isolateCluster[,col], probs = .25) &
                                                                   isolateCluster[,col] <= quantile(isolateCluster[,col], probs = .75)),]))
  }
  coreGenes <- unique(coreGenes)
  
  isolateCluster <- isolateCluster[which(rownames(isolateCluster) %in% coreGenes),]
  
  df <- data.frame()
  
  for (col in 1:(ncol(isolateCluster)-1)) {
    df <- rbind(df, data.frame(Cluster = rep(clust, times = nrow(isolateCluster)),
                               Time = rep(colnames(isolateCluster)[col], times = nrow(isolateCluster)),
                               Zscore = isolateCluster[,col]))
  } 
  clusterCentres <- rbind(clusterCentres, df)
}
clusterCentres$Time <- str_match(clusterCentres$Time, "^Col_([0-9]+)_Zscore")[,2]
clusterCentres$Cluster <- paste("Cluster", clusterCentres$Cluster)

plot <- ggplot(clusterCentres, aes(x = factor(Time, levels = c("0", "30", "90", "180")), y = Zscore)) + 
  geom_boxplot() + xlab("Time (min)") + ylab("Z-score") + facet_wrap(~Cluster) + theme_bw() +
  theme(axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        title = element_text(size = 14, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        strip.text = element_text(size = 12, colour = "black"))

png("Bjornson et al. data analysis/Figures/Clustering_boxplot.png", width = 800, height = 500)
print(plot)
dev.off()

# Plot a heatmap of Z-scores at each time point in each cluster.
fuzzyClusteringData <- read.csv("Bjornson et al. data analysis/Results/Clustering_results_unfiltered.csv")[,-1]
fuzzyClusteringData$cluster <- factor(fuzzyClusteringData$cluster, 
                                      levels = c(2,3,5,1,4,6))

colnames(fuzzyClusteringData) <- c("0", "30", "90", "180", "cluster")

plot <- pheatmap(fuzzyClusteringData[order(fuzzyClusteringData$cluster),c(1:4)],
                 color = rev(paletteer_d("colorBlindness::ModifiedSpectralScheme11Steps")),
                 cluster_rows = FALSE,
                 cluster_cols = FALSE,
                 show_rownames = FALSE,
                 annotation_names_row = FALSE,
                 annotation_legend = FALSE,
                 scale = "row", fontsize = 14,
                 angle_col = 0)


png("Bjornson et al. data analysis/Figures/Clustering_heatmap.png", width = 600, height = 400)
print(plot)
dev.off()
