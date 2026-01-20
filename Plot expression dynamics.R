source("Functions/loadLibraries.R")
loadLibraries()

# Plot expression dynamics for genes of interest
allGenotypes <- c("Col", "sdg2", "atx1", "jmj14")

genesOfInterest <- c("AT3G25250","AT3G48090","AT1G66600","AT1G18570",
                     "AT4G23550","AT5G13080","AT2G45220"," AT2G30770",
                     "AT3G26830","AT4G39950")

plotData <- data.frame()
for (genotype in allGenotypes) {
  for (time in c("F30", "F90", "F180")) {
    df <- read.csv(paste0("Results/", genotype, "/", time, "_vs_F0.csv"))
    df <- df[which(df$X %in% genesOfInterest),]
    
    plotData <- rbind(plotData, data.frame(Gene = df[,1],
                                           Genotype = genotype,
                                           Time = as.numeric(str_match(time, "F([0-9]+)")[,2]),
                                           logFC = df[,2]))
  }
}

plot <- ggplot(plotData, aes(x = Time, y = logFC,
                             colour = factor(Genotype, levels = c("Col", "sdg2", "atx1", "jmj14")))) +
  geom_point() + geom_line() +
  scale_x_discrete(limits=c(30,90,180)) +
  facet_wrap(~Gene) +
  theme_bw() +
  xlab("Time (min)") +
  labs(colour = "Genotype")

# Plot expression dynamics for memory genes
allGenotypes <- c("Col", "sdg2", "atx1", "jmj14")

genesOfInterest <- read_xlsx("Memory genes/Type II Memory Gene list - Jake Harris.xlsx")

plotData <- data.frame()
for (genotype in allGenotypes) {
  for (time in c("F30", "F90", "F180")) {
    df <- read.csv(paste0("Results/", genotype, "/", time, "_vs_F0.csv"))
    df <- df[which(df$X %in% genesOfInterest$Gene),]
    
    plotData <- rbind(plotData, data.frame(Gene = df[,1],
                                           Genotype = genotype,
                                           Time = as.numeric(str_match(time, "F([0-9]+)")[,2]),
                                           logFC = df[,2]))
  }
}
for (gene in unique(plotData$Gene)) {
  plot <- ggplot(plotData[which(plotData$Gene==gene),],
                 aes(x = Time, y = logFC,
                     colour = factor(Genotype, levels = c("Col", "sdg2", "atx1", "jmj14")))) +
    geom_point() + geom_line() +
    scale_x_discrete(limits=c(30,90,180)) +
    theme_bw() +
    xlab("Time (min)") +
    labs(title = gene, colour = "Genotype") +
    theme(axis.title = element_text(size = 14, colour = "black"),
          axis.text = element_text(size = 12, colour = "black"),
          title = element_text(size = 14, colour = "black"),
          legend.text = element_text(size = 12, colour = "black"),
          legend.title = element_text(size = 14, colour = "black"))
  
  png(paste0("Memory genes/Figures/",gene, "_logFC.png"), width = 600, height = 400)
  print(plot)
  dev.off()
}
