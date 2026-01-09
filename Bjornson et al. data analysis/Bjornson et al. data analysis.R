library(stringr)
library(DESeq2)
library(hash)
library(ggplot2)
library(pheatmap)
library(paletteer)
library(GenomicRanges)
library(operators)
library(dplyr)
library(readxl)

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

# Calculate Z-scores
source("Functions/Calculate_Zscores.R")
DEGs_normCounts <- calculate_Zscore(DEGs_normCounts)
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
  
  write.csv(res, paste0("Bjornson et al. data analysis/F", time, "_vs_F0.csv"))
  
  # Get list of DEGs with logFC < -1 or > 1 and padj <= .01
  DEGs <- append(DEGs, rownames(res[which(res[,1]>=1 & res[,2]<=.01 | res[,1]<=-1 & res[,2]<=.01),]))
}

DEGs <- unique(DEGs)
DEGs_normCounts <- normDDS[which(rownames(normDDS) %in% DEGs),]

# Calculate Z-scores
source("Functions/Calculate_Zscores.R")
DEGs_normCounts <- calculate_Zscore(DEGs_normCounts)
write.csv(DEGs_normCounts, "Bjornson et al. data analysis/DEGs_vs_F0.csv")

################################################
# Compare DEG lists
DEGs_vs_mock <- as.data.frame(read.csv("Bjornson et al. data analysis/DEGs_vs_mock.csv")) # 6743 genes
DEGs_vs_F0 <- as.data.frame(read.csv("Bjornson et al. data analysis/DEGs_vs_F0.csv")) # 7941 genes

# 6140 DEGs in common (high-confidence DEGs)
write.csv(DEGs_vs_F0[which(DEGs_vs_F0$X %in% DEGs_vs_mock$X),],
          "Bjornson et al. data analysis/High_conf_DEGs.csv")