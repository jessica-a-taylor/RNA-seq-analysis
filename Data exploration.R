source("Functions/loadLibraries.R")
loadLibraries()

allGenotypes <- c("Col", "sdg2", "atx1", "jmj14")
allTimes <- c("F0", "F30", "F90", "F180")

# 1. How many (and which) genes are differentially expressed in WT and mutants?
DEG_total <- data.frame()
DEG_counts <- data.frame()

for (genotype in allGenotypes) {
  df <- c()
  for (time in allTimes[-1]) {
    df_temp <- read.csv(paste0("Results/", genotype, "/", time, "_vs_F0.csv"))
    df_temp <- df_temp[which(df_temp[,2]>=1 & df_temp[,3]<=.01 | df_temp[,2]<=-1 & df_temp[,3]<=.01),1]
    
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

png("Figures/DEG_counts.png", width = 600, height = 400)
print(plot)
dev.off()

# 2. How many (and which) genes have altered expression in mutants relative to WT?
DEG_total <- data.frame()
DEG_counts <- data.frame()

for (genotype in allGenotypes[-1]) {
  df <- c()
  for (time in allTimes[-1]) {
    df_temp <- read.csv(paste0("Results/", genotype, "/", time, "_vs_WT.csv"))
    df_temp <- df_temp[which(df_temp[,2]>=1 & df_temp[,3]<=.05 | df_temp[,2]<=-1 & df_temp[,3]<=.05),1]
    
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


png("Figures/Mutant_vs_WT_DEG_counts.png", width = 600, height = 400)
print(plot)
dev.off()

