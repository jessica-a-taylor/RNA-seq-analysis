loadLibraries <- function() {
  library(stringr)
  library(openxlsx)
  library(DESeq2)
  library(hash)
  library(cluster)
  library(ggplot2)
  library(gplots)
  library(Mfuzz)
}

run_DESeq <- function(dds, treatment) {
  dds_output <- DESeq(dds)
  
  normalisedCounts <- data.frame(rownames(dds))
  
  # For each time point
  for (time in str_match(times, "^([0-9]+).*$")[,2]) {
    
    # Get the results of the DESeq analysis between mock and treatment.
    res <- as.data.frame(results(dds_output, contrast = c("conditions",
                                                          levels(dds$conditions)[which(grepl(paste0(treatment, time), levels(dds$conditions))==TRUE)],
                                                          levels(dds$conditions)[which(grepl(paste0("mock", time), levels(dds$conditions))==TRUE)])))
    
    
    # Save the log2 fold change and p-values.
    res <- res[,c(2,6)]
    colnames(res) <- c(paste0(treatment, time, "_FC"), paste0(treatment, time, "_padj"))
    
    # Save the normalised read counts for the treatment condition.
    res <- cbind(res, as.data.frame(counts(dds_output, normalized=TRUE)[,which(grepl(paste0(treatment, time), colnames(counts(dds_output, normalized=TRUE))) == TRUE)]))
    
    normalisedCounts <- cbind(normalisedCounts, res)
  }
  # Calculate the average normalised read counts for the mock condition.  
  # This will be considered the basal expression level for the enrichment analysis.  
  mockNormCounts <- as.data.frame(counts(dds_output, normalized=TRUE)[,which(grepl("mock", colnames(counts(dds_output, normalized=TRUE))) == TRUE)])
  mockNormCounts <- data.frame(Gene = rownames(mockNormCounts),
                               Counts = as.data.frame(rowMeans(mockNormCounts))[,1])
  write.xlsx(mockNormCounts, paste0("Data/Mock_vs_", treatment, "_normCounts.xlsx"))
  
  return(normalisedCounts)
}

filterDEGs <- function(normalisedCounts, logFC, pvalue, DEG_output, treatment) {
  signif_genes <- c()
  
  for (row in 1:nrow(normalisedCounts)) {
    p <- 0
    q <- 0
    
    for (col in which(grepl("_FC",colnames(normalisedCounts)))) {
      if (is.na(normalisedCounts[row,col])==TRUE) {
        p <- p
      }
      else if (normalisedCounts[row,col] < -logFC | normalisedCounts[row,col] > logFC) {
        p <- p + 1
      }
    }
    
    for (col in which(grepl("_padj",colnames(normalisedCounts)))) {
      if (is.na(normalisedCounts[row,col])==TRUE) {
        q <- q
      }
      else if (normalisedCounts[row,col] < pvalue) {
        q <- q + 1
      }
    }
    
    if (p >= 1 & q >= 1) {
      signif_genes <- append(signif_genes, rownames(normalisedCounts[row,]))
    }
  }
  rm(p, q, row)
  
  normalisedCounts <- normalisedCounts[which(rownames(normalisedCounts) %in% signif_genes),]
  normalisedCounts <- cbind(Gene = rownames(normalisedCounts), normalisedCounts)
  writeData(DEG_output, sheet = paste(treatment, "logFC", sep = "_"), normalisedCounts[,which(grepl("Gene", colnames(normalisedCounts)) |
                                                                                                grepl("FC", colnames(normalisedCounts)))])
  
  normalisedCounts <- normalisedCounts[,-which(grepl("Gene", colnames(normalisedCounts)) |
                                                 grepl("FC", colnames(normalisedCounts)) |
                                                 grepl("padj", colnames(normalisedCounts)))]
  return(normalisedCounts)
}

calculate_Zscores <- function(normalisedCounts, treatment) {
  # Calculate the average counts per time point.
  avgColNames <- c()
  for (time in str_match(times, "^([0-9]+).*$")[,2]) {
    normalisedCounts <- cbind(normalisedCounts, 
                              rowMeans(normalisedCounts[,which(colnames(normalisedCounts) %in% str_match(colnames(normalisedCounts), paste0("^(",treatment,time,").*$"))[,1])]))
    
    avgColNames <- append(avgColNames, paste0(treatment,time,"_avg"))
  }
  colnames(normalisedCounts) <- c(colnames(normalisedCounts)[1:(ncol(normalisedCounts)-length(str_match(names(colData[[treatment]]), "^([0-9]+).*$")[,2]))],
                                  avgColNames)
  
  # Calculate mean and sd across conditions.
  normalisedCounts$pop_avg <- rowMeans(normalisedCounts[,avgColNames])
  normalisedCounts$pop_sd <- rowSds(as.matrix(normalisedCounts[,avgColNames]))
  
  # Calculate Z-scores.
  ZcolNames <- c()
  for (col in avgColNames) {
    normalisedCounts <- cbind(normalisedCounts, (normalisedCounts[,col]-normalisedCounts$pop_avg)/normalisedCounts$pop_sd)
    
    ZcolNames <- append(ZcolNames, paste0(str_match(col, "^([a-zA-Z]+[0-9]+_)avg")[,2],"Zscore"))
  }
  colnames(normalisedCounts) <- c(colnames(normalisedCounts)[1:(ncol(normalisedCounts)-length(str_match(names(colData[[treatment]]), "^([0-9]+).*$")[,2]))],
                                  ZcolNames)
  
  return(normalisedCounts)
}

myheatmap <- function (x,bar) {
  mycolors <- c("darkgoldenrod2","darkslategray3","darkkhaki", "firebrick2", "darkolivegreen3", "darkorange","darkorchid","darksalmon", "darkgray", "maroon3")
  hmcols <- colorRampPalette(c("firebrick", "black", "darkgoldenrod2"))
  
  ngenes = as.character(table(bar)) # number of genes per cluster
  
  # sample 2000 random genes to plot
  # Sample genes best representing each cluster?
  if (length(bar) < 2000) {
    x <- x
  } else {
    ix = sort(sample(1:length(bar),2000)); bar = bar[ix]; x = x[ix,]
  }
  
  # this will cutoff very large values, which could skew the color 
  x = as.matrix(x)-apply(x,1,mean) # does this need to happen twice?
  cutoff = median(unlist(x)) + 3*sd (unlist(x)) 
  x[x>cutoff] <- cutoff
  cutoff = median(unlist(x)) - 3*sd (unlist(x)) 
  x[x< cutoff] <- cutoff
  
  # Change colnames for x axis labels.
  colnames(x) <- str_match(colnames(x), "^[a-zA-Z]+([0-9]+)$")[,2]
  
  heatmap.2(x, Rowv =F,Colv=F, dendrogram ="none",
        col=hmcols, density.info="none", key = F, keysize=.1,
        trace="none", scale="none", labRow = F,
        RowSideColors = mycolors[bar], srtCol = 0,
        xlab = "Time (min)", cexCol = 1.25, adjCol = c(.5,.3),
        margins = c(3,15))

  legend.text = paste("Cluster ", toupper(letters)[unique(bar)], " (n = ", ngenes,")", sep="") 
  
  par(lend = 1)           # square line ends for the color legend
  legend("right",      # location of the legend on the heatmap plot
         legend = legend.text, # category labels
         yjust = 0.5,
         col = mycolors,  # color key
         lty= 1,            # line style
         lwd = 10,
         bty = "n",
         cex = 0.75)
}

howManyClusters <- function(data) {
  sil <- rep(0, 20)
  
  # repeat k-means for 1:20 and extract silhouette:
  for(i in 2:20){
    k1to20 <- kmeans(data, centers = i, nstart = 25, iter.max = 20)
    ss <- silhouette(k1to20$cluster, dist(data))
    sil[i] <- mean(ss[, 3])
  }
  
  # Plot the  average silhouette width
  plot(1:20, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
  abline(v = which.max(sil), lty = 2)
}

kmeansClustering <- function(data, k, treatment) {
  set.seed(20)
  n <- ncol(data)
  
  kClust <- kmeans(data, centers=k, nstart = 1000, iter.max = 50)
  data$cluster <- kClust$cluster
  
  hc <- hclust(as.dist(1-cor(t(kClust$centers-apply(kClust$centers,1,mean)), method="pearson"))) # perform cluster for the reordering of samples
  tem = match(kClust$cluster,hc$order) #  new order 
  data = data[order(tem),] ; 	bar = sort(tem)
  
  clusterCentres <- data.frame()
  for (clust in 1:k) {
    isolateCluster <- data[data$cluster==clust,]
    
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
                                 Treatment = rep(colnames(isolateCluster)[col], times = nrow(isolateCluster)),
                                 Zscore = isolateCluster[,col]))
    } 
    clusterCentres <- rbind(clusterCentres, df)
  }
  
  plot <- ggplot(clusterCentres, aes(x = factor(Treatment, levels = c("flg0", "flg5", "flg10", "flg30", "flg90", "flg180")), y = Zscore)) + 
    geom_boxplot() + xlab("Time") + ylab("Z-score") + facet_wrap(~Cluster)
  
  pdf(paste("Graphs/", treatment, "_Zscore_boxplot", ".pdf", sep = ""), width = 12, height = 5)
  print(plot)
  dev.off() 
  
  heatmap <- myheatmap(data[1:n]-apply(data[1:n],1,mean), bar)
  
  pdf(paste("Graphs/", treatment, "_heatmap", ".pdf", sep = ""), width = 12, height = 5)
  print(heatmap)
  dev.off()
  
  return(data)
}

fuzzykmeansClustering <- function(data, k, treatment) {
  set.seed(20)
  n <- ncol(data)
  
  # Convert data to ExpressionSet for mfuzz.
  fuzzyClusteringData <- ExpressionSet(assayData = as.matrix(data))
  m <- mestimate(fuzzyClusteringData)
  fuzzClust <- mfuzz(fuzzyClusteringData, centers = k, m = m)
  
  # Remove genes that cannot be assigned to a cluster with a membership value > 0.7.
  goodGenes <- acore(fuzzyClusteringData, cl=fuzzClust, min.acore = .7)
  
  genesToCluster <- c()
  for (c in 1:length(goodGenes)) {
    genesToCluster <- append(genesToCluster, which(names(fuzzClust$cluster) %in% goodGenes[[c]]$NAME))
  }
  
  fuzzClust$cluster <- fuzzClust$cluster[genesToCluster]
  fuzzyClusteringData <- ExpressionSet(assayData = fuzzyClusteringData@assayData$exprs[which(rownames(fuzzyClusteringData@assayData$exprs) 
                                                                                             %in% names(fuzzClust$cluster)),])
  
  # Repeat clustering on high-confidence genes.
  m <- mestimate(fuzzyClusteringData)
  fuzzClust <- mfuzz(fuzzyClusteringData, centers = k, m = m)
  
  fuzzyClusteringData <- as.data.frame(fuzzyClusteringData@assayData$exprs)
  fuzzyClusteringData$cluster <- as.data.frame(fuzzClust$cluster)[,1]
  
  # Plot a boxplot of Z-scores at each time point in each cluster.
  clusterCentres <- data.frame()
  for (clust in 1:10) {
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
                                 Treatment = rep(colnames(isolateCluster)[col], times = nrow(isolateCluster)),
                                 Zscore = isolateCluster[,col]))
    } 
    clusterCentres <- rbind(clusterCentres, df)
  }
  
  clusterCentres$Treatment <- str_match(clusterCentres$Treatment, "^[a-zA-Z]+([0-9]+)$")[,2]
  clusterCentres$Cluster <- paste("Cluster", clusterCentres$Cluster)
  
  clusterBoxplot <- ggplot(clusterCentres, aes(x = factor(Treatment, levels = c("0", "5", "10", "30", "90", "180")), y = Zscore)) + 
    geom_boxplot() + xlab("Time (min)") + ylab("Z-score") + facet_wrap(~Cluster) + theme_bw() 
  
  pdf(paste("Graphs/", treatment, "_Zscore_boxplot", ".pdf", sep = ""), width = 8, height = 5)
  print(clusterBoxplot)
  dev.off()
  
  # Plot heatmap.
  hc <- hclust(as.dist(1-cor(t(fuzzClust$centers-apply(fuzzClust$centers,1,mean)), method="pearson"))) # perform cluster for the reordering of samples
  tem = match(fuzzClust$cluster,hc$order) #  new order 
  fuzzyClusteringData = fuzzyClusteringData[order(tem),] ; 	bar = sort(tem)
  
  myheatmap(fuzzyClusteringData[1:n]-apply(fuzzyClusteringData[1:n],1,mean), bar)
  
  return(fuzzyClusteringData)
}
