calculate_Zscore <- function(normCounts) {
  
  # Calculate the average counts per time point.
  avgColNames <- c()
  for (time in c("0","30","90","180")) {
    normCounts <- cbind(normCounts, rowMeans(as.matrix(normCounts[,which(str_detect(colnames(normCounts), time)==TRUE)])))
  
    avgColNames <- append(avgColNames, paste0(time, "_avg"))
  }
  colnames(normCounts) <- c(colnames(normCounts)[1:(ncol(normCounts)-length(avgColNames))], avgColNames)
 
  # Calculate mean and sd across conditions.
  normCounts <- cbind(normCounts, rowMeans(as.matrix(normCounts[,avgColNames])))
  normCounts <- cbind(normCounts, rowSds(as.matrix(normCounts[,avgColNames])))
  
  colnames(normCounts) <- c(colnames(normCounts)[1:(ncol(normCounts)-2)], "pop_avg", "pop_sd")
  
  # Calculate Z-scores.
  ZcolNames <- c()
  for (col in which(str_detect(colnames(normCounts), ".*[0-9]+_avg") == TRUE)) {
    normCounts <- cbind(normCounts, (normCounts[,col]-normCounts[,"pop_avg"])/normCounts[,"pop_sd"])
    
    ZcolNames <- append(ZcolNames, paste0(str_match(colnames(normCounts)[col], "^(.*)avg")[,2],"Zscore"))
  }
  colnames(normCounts) <- c(colnames(normCounts)[1:(ncol(normCounts)-4)], ZcolNames)
  return(normCounts) 
}