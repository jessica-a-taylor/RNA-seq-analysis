source("Functions/loadLibraries.R")
loadLibraries()

# Fuzzy k-means analysis on all upregulated DEGs.
normCounts <- as.data.frame(read.csv("Results/Col/DEGs.csv"))[,1:13]
colnames(normCounts)[1] <- "gene_id"
rownames(normCounts) <- normCounts$gene_id

# Filter for upregulated genes.
upreg_DEGs <- c()

for (time in c("F30", "F90", "F180")) {
  df <- read.csv(paste0("Results/Col/", time, "_vs_F0.csv"))
  df <- df[which(df[,2]>=1 & df[,3]<=.05),1]
  
  upreg_DEGs <- append(upreg_DEGs, df)
}

normCounts <- normCounts[which(normCounts$gene_id %in% upreg_DEGs),]

# Calculate Z-scores.
source("Functions/Calculate_Zscores.R")
normCounts <- calculate_Zscore(normCounts, "Col")

# Remove mitochondrial and chloroplast genes.
normCounts <- normCounts[which(str_detect(normCounts$gene_id, "^(AT[1-5]G).*")==TRUE),]

set.seed(20)

Zscores <- as.matrix(normCounts[,which(str_detect(colnames(normCounts), "Zscore")==TRUE)])

fuzzyClusteringData <- ExpressionSet(assayData = as.matrix(Zscores))
m <- mestimate(fuzzyClusteringData)
fuzzClust <- mfuzz(fuzzyClusteringData, centers = 6, m = m)

unfilteredResults <- as.data.frame(fuzzyClusteringData@assayData$exprs)
unfilteredResults$cluster <- fuzzClust$cluster
write.csv(unfilteredResults, "Results/Col/Clustering_results_unfiltered.csv")

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

write.csv(fuzzyClusteringData, "Results/Col/Clustering_results_filtered.csv")

# Plot a boxplot of Z-scores at each time point in each cluster.
fuzzyClusteringData <- read.csv("Results/Col/Clustering_results_filtered.csv")[,-1]

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

png("Figures/Clustering_boxplot.png", width = 800, height = 500)
print(plot)
dev.off()

# Plot a heatmap of Z-scores at each time point in each cluster.
fuzzyClusteringData <- read.csv("Results/Col/Clustering_results_unfiltered.csv")[,-1]
fuzzyClusteringData$cluster <- factor(fuzzyClusteringData$cluster, 
                                      levels = c(6,3,5,4,2,1))

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


png("Figures/Clustering_heatmap.png", width = 600, height = 400)
print(plot)
dev.off()
