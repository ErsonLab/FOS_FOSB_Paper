#Software Website: https://bioinformatics.mdanderson.org/public-software/estimate/

# Install and load the required packages
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("reshape2", quietly = TRUE)) {
  install.packages("reshape2")
}
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  install.packages("gridExtra")
}
if (!requireNamespace("rstatix", quietly = TRUE)) {
  install.packages("rstatix")
}
if (!requireNamespace("estimate", quietly = TRUE)) {
  rforge <- "http://r-forge.r-project.org"
  install.packages("estimate", repos=rforge, dependencies=TRUE)
}
if (!requireNamespace("utils", quietly = TRUE)) {
  install.packages("utils")
}
if (!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) {
  install.packages("EnsDb.Hsapiens.v86")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}
if (!requireNamespace("ggsignif", quietly = TRUE)) {
  install.packages("ggsignif")
}

library(ggplot2)
library(reshape2)
library(gridExtra)
library(rstatix)
library(estimate)
library(utils)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(RColorBrewer)
library(ggsignif)

#INPUT: COMBAT VST Count Matrix
COMBAT_VST_COUNT_MATRIX <- read.table("combat_counts_vst.txt", sep = "\t", header = TRUE)

#Add the gene symbols and make them row names 
rownames(COMBAT_VST_COUNT_MATRIX) <- gsub("\\.[0-9]*", "", rownames(COMBAT_VST_COUNT_MATRIX))
COMBAT_VST_COUNT_MATRIX$gene_symbol <- mapIds(EnsDb.Hsapiens.v86,
                                              keys = row.names(COMBAT_VST_COUNT_MATRIX),
                                              column = "SYMBOL",
                                              keytype = "GENEID",
                                              multiVals = "first")

aggregated_data <- COMBAT_VST_COUNT_MATRIX %>%
  group_by(gene_symbol) %>%
  summarise(across(everything(), sum, na.rm = TRUE))

aggregated_data <- aggregated_data[-nrow(aggregated_data),]
rownames(aggregated_data) <- aggregated_data$gene_symbol
aggregated_data <- as.data.frame(aggregated_data)
rownames(aggregated_data) <- aggregated_data$gene_symbol
aggregated_data <- aggregated_data[,-1]
COMBAT_VST_COUNT_MATRIX_ESTIMATE <- aggregated_data

write.table(COMBAT_VST_COUNT_MATRIX_ESTIMATE, "ESTIMATE.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
filterCommonGenes("ESTIMATE.txt", output.f="COMBAT_VST_COUNT_MATRIX_ESTIMATE_GENES.gct", id="GeneSymbol")
estimateScore("COMBAT_VST_COUNT_MATRIX_ESTIMATE_GENES.gct", "COMBAT_VST_COUNT_MATRIX_ESTIMATE_SCORE.gct", platform="affymetrix")
plotPurity(scores="COMBAT_VST_COUNT_MATRIX_ESTIMATE_SCORE.gct", samples="n=742", platform="affymetrix")

gct_data <- read.table("COMBAT_VST_COUNT_MATRIX_ESTIMATE_SCORE.gct", skip = 2, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
umap_clusters <- read.table("09172024_umap_df_3d_742_cases_8_clusters.csv", sep = ",", header = TRUE)

# Convert to long format using melt
gct_data_long <- melt(gct_data, id.vars = c("NAME", "Description"), variable.name = "Sample", value.name = "Score")

# Match the values by merging based on Sample names
merged_data <- merge(gct_data_long, umap_clusters, by.x = "Sample", by.y = "Sample_Name")
merged_data <- merged_data %>% filter(Cluster != 0)

all_clusters_ESTIMATE_input <- merged_data

merged_data <- all_clusters_ESTIMATE_input
# Modify the 'Cluster' column, setting Cluster 2 to remain "2" and others to "all"
merged_data$Cluster <- ifelse(merged_data$Cluster == 2, "Cluster 2", "The Rest")

# Check the updated dataframe
head(merged_data)

# Generate separate plots for each score
unique_scores <- unique(merged_data$NAME)

# Provided color codes
color_palette <- c("#B79FEC", "orange")
#color_palette <- c("#FC8D62", "#B79FEC", "#E78AC3", "#C49A6C", 
#                   "#B3DE69", "#FFD92F", "#E5C494", "#B3B3B3")

# Loop through each score to create individual plots and store them in the list
plot_list <- list()
for (score in unique_scores) {
  p <- ggplot(subset(merged_data, NAME == score), aes(x = as.factor(Cluster), y = Score, fill = as.factor(Cluster))) +
    geom_boxplot() +
    labs(title = paste(score, "via ESTIMATE"),
         x = "Cluster", y = score , fill = "Cluster") +
    scale_fill_manual(values = color_palette) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size = 24),
          axis.title = element_text(size = 24, face = "bold"),
          axis.text = element_text(size = 22),
          legend.title = element_text(size = 24, face = "bold"),
          legend.text = element_text(size = 20),
          plot.title = element_text(face = "bold", size = 22))
  #geom_signif(comparisons = list(c("Cluster 2","The Rest")), map_signif_level = FALSE, test = "wilcox.test", step_increase = 0.05)
  #list(c(2,1),c(2,3),c(2,4),c(2,5),c(2,6),c(2,7),c(2,8)), map_signif_level = FALSE, test = "wilcox.test", step_increase = 0.05)
  
  plot_list[[score]] <- p
}

combined_plot <- grid.arrange(grobs = plot_list, ncol = 2)

# Save the combined plot using ggsave
ggsave("ESTIMATE_results_by_cluster.pdf", 
       plot = combined_plot,, width = 12, height = 9, dpi = 300, bg = "white")



