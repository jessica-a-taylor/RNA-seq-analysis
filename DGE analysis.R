library(stringr)
library(DESeq2)
library(hash)
library(ggplot2)
library(pheatmap)
library(paletteer)

# Import raw read counts.
rawCounts <- data.frame(Gene = read.table("Counts_data/Col-0/Col_F0_1_readCounts.txt")[,1])
for (file in list.files("Counts_data/Col-0")) {
  rawCounts <- cbind(rawCounts, read.table(paste0("Counts_data/Col-0/",file))[,2])
}

# Set row names to gene IDs
rownames(rawCounts) <- rawCounts$Gene
rawCounts <- rawCounts[,-1]

# Set column names to sample IDs
colnames(rawCounts) <- str_match(list.files("Counts_data/Col-0"), "^(.*)_readCounts.txt")[,2]

# Create metadata file.
colData <- data.frame(Sample = as.factor(colnames(rawCounts)),
                      Time = as.factor(str_match(colnames(rawCounts),"^Col_(F[0-9]+)_[0-9]+$")[,2]))

rownames(colData) <- colnames(rawCounts)

# Construct DESeq dataset.
dds <- DESeqDataSetFromMatrix(countData=rawCounts, colData=colData, design=~Time)
dds$Time <- relevel(dds$Time, ref = "F0")

# Remove genes with < 10 reads.
dds <- dds[rowSums(counts(dds)) >= 10,]

# Run the DESeq.
DDS <- DESeq(dds)
normDDS <- counts(DDS, normalized = TRUE) # normalization with respect to the sequencing depth
write.csv(normDDS, "Results/Col-0/Normalised_counts.csv")

# Save the results.
DEGs <- c()

for (time in c("F30", "F90", "F180")) {
  res <- results(DDS, contrast = c("Time", time, "F0"), independentFiltering = FALSE)
  res <- res[,c(2,6)]
  colnames(res) <- c(paste0(time, "_FC"), paste0(time, "_padj"))
  
  write.csv(res, paste0("Results/Col-0/", time, "_vs_F0.csv"))
  
  # Get list of DEGs with logFC < -1 or > 1 and padj <= .01
  DEGs <- append(DEGs, rownames(res[which(res[,1]>=1 | res[,1]<=-1 & res[,2]<=.01),]))
}

DEGs <- unique(DEGs)
DEGs_normCounts <- normDDS[which(rownames(normDDS) %in% DEGs),]

# Calculate Z-scores
source("Functions/Calculate_Zscores.R")
DEGs_normCounts <- calculate_Zscore(DEGs_normCounts)

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

fuzzClust$cluster <- fuzzClust$cluster[genesToCluster]
fuzzyClusteringData <- ExpressionSet(assayData = fuzzyClusteringData@assayData$exprs[which(rownames(fuzzyClusteringData@assayData$exprs) 
                                                                                           %in% names(fuzzClust$cluster)),])
# Repeat clustering on high-confidence genes.
m <- mestimate(fuzzyClusteringData)
fuzzClust <- mfuzz(fuzzyClusteringData, centers = 6, m = m)

fuzzyClusteringData <- as.data.frame(fuzzyClusteringData@assayData$exprs)
fuzzyClusteringData$cluster <- as.data.frame(fuzzClust$cluster)[,1]

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