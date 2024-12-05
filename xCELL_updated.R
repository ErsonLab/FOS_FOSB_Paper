# ===================================================================================================================================
#
# Paper Name: Novel Subgroup of Meningiomas Involving FOS and FOSB Gene Fusions
#
# ===================================================================================================================================
# 
# Description: This code performs XCell analysis to deconvolute bulk RNA-seq data and estimate the relative abundance of various cell types within the samples. 
# By leveraging gene expression signatures, it provides insights into the tumor microenvironment and its potential role in shaping the distinct 
# biological characteristics of the analyzed samples.
#
# Author: Kanat Yalcin, MD; Hasan Alanya, MSc
#
# Contact PI: Zeynep Erson-Omay, PhD; Murat Gunel, MD
#
# Date: 12/6/2024 
#
# Group Page: https://ersonlab.org/ and https://medicine.yale.edu/lab/gunel/
#
# GitHub: https://github.com/ErsonLab/FOS_FOSB_Paper
#
# ===================================================================================================================================

# Load required libraries
library(reshape2)
library(dplyr)
library(ggplot2)  # Ensure ggplot2 is loaded

xcell_df <- read.csv("/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/HTseq_700/Data/xCell_res.csv")
umap_df <- read.csv("/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/HTseq_700/Results/09172024_umap_df_3d_742_cases_8_clusters.csv")

# Remove Cluster 0 - this is noise
umap_df <- umap_df[umap_df$Cluster != 0, ]

# Rename columns for consistency
colnames(xcell_df)[1] <- "Cell_Type"
colnames(umap_df)[1] <- "Sample"

# Transpose xCell data so samples are columns
xcell_t <- t(xcell_df[-1])
colnames(xcell_t) <- xcell_df$Cell_Type
xcell_t <- as.data.frame(xcell_t)
xcell_t$Sample <- rownames(xcell_t)

# Merge with UMAP cluster data
merged_df <- merge(xcell_t, umap_df[, c("Sample", "Cluster")], by = "Sample")

# Convert the data from wide to long format for easier comparison
melted_df <- melt(merged_df, id.vars = c("Sample", "Cluster"), 
                  variable.name = "Cell_Type", value.name = "Enrichment")

# Perform statistical comparisons: Cluster 2 vs all others
comparison_results <- data.frame(Cell_Type = character(),
                                 p_value = numeric(),
                                 stringsAsFactors = FALSE)

# Loop through each cell type and perform Wilcoxon test for Cluster 2 vs others
cell_types <- unique(melted_df$Cell_Type)

for (cell_type in cell_types) {
  # Subset the data for cell type
  cell_data <- melted_df[melted_df$Cell_Type == cell_type, ]
  
  # Separate Cluster 2 and all other clusters - c2 vs all
  cluster_2 <- cell_data[cell_data$Cluster == 2, "Enrichment"]
  all_other_clusters <- cell_data[cell_data$Cluster != 2, "Enrichment"]
  
  # Perform Wilcoxon rank-sum test
  test_result <- wilcox.test(cluster_2, all_other_clusters)
  
  # Store the results in the comparison_results data frame
  comparison_results <- rbind(comparison_results, 
                              data.frame(Cell_Type = cell_type, 
                                         p_value = test_result$p.value))
}

# Sort the results by p-value to find the most significant comparisons
comparison_results <- comparison_results %>% arrange(p_value)

write.csv(comparison_results, "/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/HTseq_700/xCELL_cluster2_vs_all_comparison_results.csv", row.names = FALSE)

# Top 15
top_15_cell_types <- head(comparison_results$Cell_Type, 15)
print(top_15_cell_types)

write.csv(top_15_cell_types, "/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/HTseq_700/xCELL_top15_cluster2_vs_all_comparison_results.csv", row.names = FALSE)

# Color palette that Kanat gave 
color_palette <- c("#FC8D62", "#B79FEC", "#E78AC3", "#C49A6C", 
                   "#B3DE69", "#FFD92F", "#E5C494", "#B3B3B3")

# Adjusted plotting loop with matching aesthetics
for (cell_type in top_15_cell_types) {
  # Subset data for the specific cell type
  data_subset <- melted_df[melted_df$Cell_Type == cell_type, ]
  
  p <- ggplot(data_subset, aes(x = factor(Cluster), y = Enrichment, fill = factor(Cluster))) +
    geom_boxplot() +
    scale_fill_manual(values = color_palette) +
    labs(
      x = "Cluster",
      y = "Proportion",
      title = gsub("\\.", " ", cell_type),  # Remove dots from title
      fill = "Cluster"
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 40, face = "bold"),      # Increase tick labels font size
      axis.title = element_text(size = 48, face = "bold"),     # Increase axis labels font size
      legend.title = element_blank(),                          # Remove legend title
      legend.text = element_text(size = 30),                   # Increase legend text size
      legend.position = "bottom",                              # Move legend to bottom
      panel.grid = element_blank(),                            # Remove grid lines
      plot.background = element_rect(fill = "white", colour = "black", linewidth = 2),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
      plot.title = element_text(size = 48, face = "bold", hjust = 0)  # Increase plot title font size, left-aligned
    )
  
  file_name <- paste0("/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/HTseq_700/Figures/", gsub("\\.", "_", cell_type), "_xCELL_boxplot.pdf")
  
  # Adjust the width and height to accommodate the increased font sizes
  ggsave(filename = file_name, plot = p, width = 18, height = 14, dpi = 300, bg = "white")
  
  cat("Saved:", file_name, "\n")
}

######################################################################

# Initialize a list to store the plots
combined_plot_list <- list()

# Flag to extract the legend only once
legend_extracted <- FALSE
legend <- NULL  # Initialize legend variable

# Loop through the top 15 cell types to regenerate and store the plots
for (cell_type in top_15_cell_types) {
  # Subset data for the specific cell type
  data_subset <- melted_df[melted_df$Cell_Type == cell_type, ]
  
  # Create the plot with adjusted aesthetics for multi-plot layout
  p <- ggplot(data_subset, aes(x = factor(Cluster), y = Enrichment, fill = factor(Cluster))) +
    geom_boxplot() +
    scale_fill_manual(values = color_palette) +
    labs(
      x = "Cluster",
      y = "Enrichment",
      title = gsub("\\.", " ", cell_type),  # Remove dots from title
      fill = "Cluster"
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 12, face = "bold"),        # Slightly increased font sizes
      axis.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0),
      legend.position = "none"  # Suppress legend in individual plots
    )
  
  # Extract legend from the first plot
  if (!legend_extracted) {
    legend <- get_legend(
      p + theme(legend.position = "bottom",
                legend.text = element_text(size = 10),
                legend.title = element_blank())
    )
    legend_extracted <- TRUE
  }
  
  # Add the plot to the list
  combined_plot_list[[cell_type]] <- p
}

# Arrange the plots into a grid
combined_plot <- plot_grid(plotlist = combined_plot_list, ncol = 3, align = 'hv')

# Add the legend below the combined plots
final_combined_plot <- plot_grid(
  combined_plot,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.1)  # Adjust the relative heights as needed
)

# Save the combined plot as a PDF
ggsave(
  filename = "/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/HTseq_700/Figures/Combined_Top15_Boxplots.pdf",
  plot = final_combined_plot,
  width = 24,    # Adjust width as needed
  height = 18,   # Adjust height as needed
  dpi = 300,
  bg = "white"
)

# Inform the user
cat("Combined plot saved as Combined_Top15_Boxplots.pdf\n")
