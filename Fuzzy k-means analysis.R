source("Functions/loadLibraries.R")
loadLibraries()

DEGs_normCounts <- as.data.frame(read.csv("Results/Col-0/DEGs.csv"))
colnames(DEGs_normCounts)[1] <- "gene_id"

# Fuzzy k-means analysis.
Zscores <- as.matrix(DEGs_normCounts[,which(str_detect(colnames(DEGs_normCounts), "Zscore")==TRUE)])

set.seed(20)

fuzzyClusteringData <- ExpressionSet(assayData = as.matrix(Zscores))
m <- mestimate(fuzzyClusteringData)
fuzzClust <- mfuzz(fuzzyClusteringData, centers = 6, m = m)

# Remove genes that cannot be assigned to a cluster with a membership value > 0.7.
goodGenes <- acore(fuzzyClusteringData, cl=fuzzClust, min.acore = .5)
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

# Plot a boxplot of Z-scores at each time point in each cluster.
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
clusterCentres$Time <- str_match(clusterCentres$Time, "^[a-zA-Z]+([0-9]+)_Zscore")[,2]
clusterCentres$Cluster <- paste("Cluster", clusterCentres$Cluster)

clusterBoxplot <- ggplot(clusterCentres, aes(x = factor(Time, levels = c("0", "5", "10", "30", "90", "180")), y = Zscore)) + 
  geom_boxplot() + xlab("Time (min)") + ylab("Z-score") + facet_wrap(~Cluster) + theme_bw() 

fuzzyClusteringData$cluster <- factor(fuzzyClusteringData$cluster, 
                                      levels = c(4,5,3,6,2,1))

pheatmap(fuzzyClusteringData[order(fuzzyClusteringData$cluster),c(1:4)],
         color = rev(paletteer_d("colorBlindness::ModifiedSpectralScheme11Steps")),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         annotation_names_row = FALSE,
         annotation_legend = FALSE,
         scale = "row", fontsize = 16)
