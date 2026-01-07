library(stringr)
library(DESeq2)
library(hash)
library(ggplot2)

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