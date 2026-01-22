source("Functions/loadLibraries.R")
loadLibraries()

# Analyse the expression dynamics for memory genes.
allGenotypes <- c("Col", "sdg2", "atx1", "jmj14")
allTimes <- c("F0", "F30", "F90", "F180")

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
# Plot expression change for each gene.
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
# Plot average expression change.
plot <- ggplot(plotData, aes(x = factor(as.character(Time), levels = c("30","90","180")), 
                             y = logFC,
                             fill = factor(Genotype, levels = c("Col", "sdg2", "atx1", "jmj14")))) +
  geom_boxplot() +
  theme_bw() +
  xlab("Time (min)") +
  labs(fill = "Genotype") +
  theme(axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        title = element_text(size = 14, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"))

png("Memory genes/Figures/Average_logFC.png", width = 600, height = 400)
print(plot)
dev.off()

# How many (and which) memory genes are differentially expressed in WT and mutants?
DEG_total <- data.frame()
DEG_counts <- data.frame()

for (genotype in allGenotypes) {
  df <- c()
  for (time in allTimes[-1]) {
    df_temp <- read.csv(paste0("Results/", genotype, "/", time, "_vs_F0.csv"))
    df_temp <- df_temp[which(df_temp[,2]>=1 & df_temp[,3]<=.01 & df_temp[,1] %in% genesOfInterest$Gene |
                               df_temp[,2]<=-1 & df_temp[,3]<=.01 & df_temp[,1] %in% genesOfInterest$Gene),1]

    df <- append(df, df_temp)
    DEG_counts <- rbind(DEG_counts, data.frame(Genotype = genotype,
                                               Time = time,
                                               DEGs = length(df_temp)))
  }
  DEG_total <- rbind(DEG_total, data.frame(Genotype = genotype,
                                           DEGs = length(unique(df))))
}

plot <- ggplot(DEG_counts, aes(x = factor(Time, levels = c("F30", "F90", "F180")), 
                               y = DEGs, 
                               fill = factor(Genotype, levels = c("Col", "sdg2", "atx1", "jmj14")))) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() + xlab("Time (min)") + ylab("Number of DEGs") +
  labs(fill = "Genotype") +
  theme(axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        title = element_text(size = 14, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 14, colour = "black")) +
  scale_fill_manual(values = c("#808080", "#FF7C80", "#6699FF", "#99CC00"))

png("Memory genes/Figures/DEG_counts.png", width = 600, height = 400)
print(plot)
dev.off()

# How many (and which) memory genes have altered expression in mutants?
DEG_total <- data.frame()
DEG_counts <- data.frame()

for (genotype in allGenotypes[-1]) {
  df <- c()
  for (time in allTimes[-1]) {
    df_temp <- read.csv(paste0("Results/", genotype, "/", time, "_vs_WT.csv"))
    df_temp <- df_temp[which(df_temp[,2]>=1 & df_temp[,3]<=.05 & df_temp[,1] %in% genesOfInterest$Gene | 
                               df_temp[,2]<=-1 & df_temp[,3]<=.05 & df_temp[,1] %in% genesOfInterest$Gene),1]
    
    df <- append(df, df_temp)
    DEG_counts <- rbind(DEG_counts, data.frame(Genotype = genotype,
                                               Time = time,
                                               DEGs = length(df_temp)))
  }
  DEG_total <- rbind(DEG_total, data.frame(Genotype = genotype,
                                           DEGs = length(unique(df))))
}

plot <- ggplot(DEG_counts, aes(x = factor(Time, levels = c("F30", "F90", "F180")), 
                               y = DEGs, 
                               fill = factor(Genotype, levels = c("sdg2", "atx1", "jmj14")))) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() + xlab("Time (min)") + ylab("Number of DEGs") +
  labs(fill = "Genotype") +
  theme(axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        title = element_text(size = 14, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 14, colour = "black")) +
  scale_fill_manual(values = c("#FF7C80", "#6699FF", "#99CC00"))


png("Memory genes/Figures/Mutant_vs_WT_DEG_counts.png", width = 600, height = 400)
print(plot)
dev.off()

