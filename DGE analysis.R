setwd("~/nextflow-peakcaller/Differential-expression-analysis")

# Load required libraries.
source("Functions for DGE analysis.R")
loadLibraries()
rm(loadLibraries)

# Import raw read counts.
allCounts <- read.csv("Data/PRJEB25079.csv")
genes <- allCounts[,1]
allCounts <- allCounts[,-1]
row.names(allCounts) <- genes

rm(genes)

# Generate a hash mimicking an annotation file that will define the comparisons across conditions.
colData <- hash()

for (treatment in unique(str_match(colnames(allCounts), "^([a-zA-Z]+).*$")[,2])[-1]) {
  colData[[treatment]] <- hash(`0min` = data.frame(sample = colnames(allCounts)[which(grepl("mock0", colnames(allCounts))==TRUE | grepl(paste0(treatment, "0"), colnames(allCounts))==TRUE)],
                                                condition = str_match(colnames(allCounts)[which(grepl("mock0", colnames(allCounts))==TRUE | grepl(paste0(treatment, "0"), colnames(allCounts))==TRUE)], "^([A-Za-z]+0)_[0-9]+$")[,2]),
                            `5min` = data.frame(sample = colnames(allCounts)[which(grepl("mock5", colnames(allCounts))==TRUE | grepl(paste0(treatment, "5"), colnames(allCounts))==TRUE)],
                                                condition = str_match(colnames(allCounts)[which(grepl("mock5", colnames(allCounts))==TRUE | grepl(paste0(treatment, "5"), colnames(allCounts))==TRUE)], "^([A-Za-z]+5)_[0-9]+$")[,2]),
                            `10min` = data.frame(sample = colnames(allCounts)[which(grepl("mock10", colnames(allCounts))==TRUE | grepl(paste0(treatment, "10"), colnames(allCounts))==TRUE)],
                                                 condition = str_match(colnames(allCounts)[which(grepl("mock10", colnames(allCounts))==TRUE | grepl(paste0(treatment, "10"), colnames(allCounts))==TRUE)], "^([A-Za-z]+10)_[0-9]+$")[,2]),
                            `30min` = data.frame(sample = colnames(allCounts)[which(grepl("mock30", colnames(allCounts))==TRUE | grepl(paste0(treatment, "30"), colnames(allCounts))==TRUE)],
                                                 condition = str_match(colnames(allCounts)[which(grepl("mock30", colnames(allCounts))==TRUE | grepl(paste0(treatment, "30"), colnames(allCounts))==TRUE)], "^([A-Za-z]+30)_[0-9]+$")[,2]),
                            `90min` = data.frame(sample = colnames(allCounts)[which(grepl("mock90", colnames(allCounts))==TRUE | grepl(paste0(treatment, "90"), colnames(allCounts))==TRUE)],
                                                 condition = str_match(colnames(allCounts)[which(grepl("mock90", colnames(allCounts))==TRUE | grepl(paste0(treatment, "90"), colnames(allCounts))==TRUE)], "^([A-Za-z]+90)_[0-9]+$")[,2]),
                            `180min` = data.frame(sample = colnames(allCounts)[which(grepl("mock180", colnames(allCounts))==TRUE | grepl(paste0(treatment, "180"), colnames(allCounts))==TRUE)],
                                                  condition = str_match(colnames(allCounts)[which(grepl("mock180", colnames(allCounts))==TRUE | grepl(paste0(treatment, "180"), colnames(allCounts))==TRUE)], "^([A-Za-z]+180)_[0-9]+$")[,2]))
}


# Create workbook into which the output data will be saved.
DGE_output <- createWorkbook()
DGE_output_filtered <- createWorkbook()

times = c("0min", "5min", "10min", "30min", "90min", "180min")

for (treatment in names(colData)) {
  addWorksheet(DGE_output, sheetName = treatment)
  addWorksheet(DGE_output_filtered, sheetName = paste(treatment, "normCounts", sep = "_"))
  addWorksheet(DGE_output_filtered, sheetName = paste(treatment, "logFC", sep = "_"))
  
  # Construct DESeq dataset.
  comparison <- data.frame()
  for (time in times) {
    comparison <- rbind(comparison, data.frame(samples = colData[[treatment]][[time]]$sample,
                                               conditions = colData[[treatment]][[time]]$condition))
  }
  counts <- allCounts[,comparison$samples]
  dds <- DESeqDataSetFromMatrix(countData=counts, colData=comparison, design=~conditions)
  rm(comparison)
  
  # Remove genes with < 10 reads.
  dds <- dds[rowSums(counts(dds)) >= 10,]
  
  # Run the DESeq pipeline.
  normalisedCounts <- run_DESeq(dds, treatment)
  rm(dds, counts)
  
  allNormalisedCounts <- data.frame(Genes = rownames(normalisedCounts))
  for (time in str_match(times, "^([0-9]+).*$")[,2]) {
    allNormalisedCounts <- cbind(allNormalisedCounts, data.frame(rowMeans(normalisedCounts[,which(!is.na(str_extract(colnames(normalisedCounts), paste0("^(",treatment,time,"_[1-9]+).*$"))))])))
    
  }
  colnames(allNormalisedCounts) <- c("Genes",paste0(treatment, str_match(times, "^([0-9]+).*$")[,2],"_avg"))
  writeData(DGE_output, sheet = treatment, allNormalisedCounts)
  
  # Filter for genes with logFC < -1 or > 1 and padj <= .01 in at least one time point
  normalisedCounts <- filterDEGs(normalisedCounts, 1, .01, DGE_output_filtered, treatment)
  
  # Calculate Z-scores.
  normalisedCounts <- calculate_Zscores(normalisedCounts, treatment)
  
  # Save output.
  normalisedCounts <- cbind(Gene = rownames(normalisedCounts), normalisedCounts)
  writeData(DGE_output_filtered, sheet = paste(treatment, "normCounts", sep = "_"), normalisedCounts)
}
saveWorkbook(DGE_output, "DGE analysis output.xlsx", overwrite = TRUE) 
saveWorkbook(DGE_output_filtered, "DGE analysis output (filtered).xlsx", overwrite = TRUE) 
rm(run_DESeq, filterDEGs, calculate_Zscores)