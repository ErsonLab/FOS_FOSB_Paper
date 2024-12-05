# ===================================================================================================================================
#
# Paper Name: Novel Subgroup of Meningiomas Involving FOS and FOSB Gene Fusions
#
# ===================================================================================================================================
# 
# Description: This code utilizes CIBERSORTx to estimate the relative proportions of immune and stromal cell types from bulk RNA-seq data, 
# providing insights into the cellular composition of tumor samples. It includes preprocessing steps to integrate cell-type proportions with 
# clustering results and visualizes the distribution of various cell types across clusters using box plots. Statistical analysis is performed
# using the Wilcoxon test to identify significant differences between clusters. The results are visualized with customized plots, highlighting 
# cell type distributions and significant pairwise comparisons, offering a detailed understanding of the tumor microenvironment's composition 
# and variability.
#
# Author: Hasan Alanya, MSc ; Kanat Yalcin, MD
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

# Load the required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)

# Load the datasets 
cibersortx_results <- read.csv("/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/HTseq_700/Results/CIBERSORTx/CIBERSORTx_Job4_Results.csv")
umap_clusters <- read.csv("/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/HTseq_700/Results/09172024_umap_df_3d_742_cases_8_clusters.csv")

# Rename the column to match
colnames(umap_clusters)[1] <- "Mixture" 

# Remove Cluster 0 - this is noise
umap_clusters <- umap_clusters[umap_clusters$Cluster != 0, ]

# Merge the data on "Mixture"
merged_data <- merge(cibersortx_results, umap_clusters, by = "Mixture")

# Reshape data for plotting
long_data <- merged_data %>%
  dplyr::select(Mixture, Cluster, `B.cells.naive`, `T.cells.CD8`, `T.cells.CD4.memory.resting`,
                `T.cells.follicular.helper`, `NK.cells.resting`, `NK.cells.activated`, `Monocytes`,
                `Macrophages.M2`, `Mast.cells.resting`, `Neutrophils`) %>%
  pivot_longer(cols = -c(Mixture, Cluster), names_to = "Cell_Type", values_to = "Value")


color_palette <- c("#FC8D62", "#B79FEC", "#E78AC3", "#C49A6C", 
                   "#B3DE69", "#FFD92F", "#E5C494", "#B3B3B3")

                    
box_plot <- ggplot(long_data, aes(x = factor(Cluster), y = Value, fill = factor(Cluster))) +
  geom_boxplot() +
  facet_wrap(~Cell_Type, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = color_palette) + 
  labs(
    x = "Cluster", 
    y = "Proportion", 
    title = "Distribution of Cell Types across Clusters", 
    fill = "Cluster"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 30, face = "bold"),      # Match tick labels font size
    axis.title = element_text(size = 40, face = "bold"),     # Match axis labels font size
    legend.title = element_blank(),                          # Match legend title
    legend.text = element_text(size = 30),                   # Match legend text size
    legend.position = "right",                              # Match legend position
    panel.grid = element_blank(),                            # Remove grid lines
    plot.background = element_rect(fill = "white", colour = "black", linewidth = 2),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    plot.title = element_text(size = 48, face = "bold", hjust = 0.5),  # Match plot title font size
    strip.text = element_text(size = 30, face = "bold", hjust = 0),
    panel.spacing = unit(2, "lines"),           # Increase space between facets
    plot.margin = unit(c(1, 1, 1, 1), "cm")     # Match plot title font size
  )

# Save the updated plot
ggsave(
  "/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/HTseq_700/Figures/CIBERSORTx_box_plot_updated.png", 
  plot = box_plot, width = 30, height = 28, dpi = 300, bg = "white"
)


# Wilcoxon test with BH correction
comparison_results <- long_data %>%
  group_by(Cell_Type) %>%
  pairwise_wilcox_test(Value ~ Cluster, p.adjust.method = "BH") %>%
  filter(p.adj.signif != "ns" & p.adj.signif != "*") %>%
  mutate(label = paste0("p.adj = ", format.pval(p.adj, digits = 3)))

# Split comparison results by cell type
list_of_dfs <- split(comparison_results, comparison_results$Cell_Type)

# Create a list of comparison pairs for significance annotation
list_of_lists <- lapply(list_of_dfs, function(df) {
  lapply(seq_len(nrow(df)), function(i) {
    c(as.character(df[i, "group1"]), as.character(df[i, "group2"]))
  })
})

# Create a list to store plots
plot_list <- list()

# Loop through each cell type to create individual plots
for(cell_type in unique(long_data$Cell_Type)) {
  # Filter data for the current cell type
  data_subset <- long_data %>% filter(Cell_Type == cell_type)
  
  # Define the plot for the current cell type
  p <- ggplot(data_subset, aes(x = factor(Cluster), y = Value, fill = factor(Cluster))) +
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
      axis.text = element_text(size = 30, face = "bold"),      # Increase tick labels font size
      axis.title = element_text(size = 35, face = "bold"),     # Increase axis labels font size
      legend.title = element_blank(),                          # Remove legend title
      legend.text = element_text(size = 30),                   # Increase legend text size
      legend.position = "right",                              # Move legend to bottom
      panel.grid = element_blank(),                            # Remove grid lines
      plot.background = element_rect(fill = "white", colour = "white", linewidth = 2),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
      plot.title = element_text(size = 40, face = "bold", hjust = 0)  # Increase plot title font size
    )
  
  # Add significance markers if comparisons exist for the current cell type
  if (!is.null(list_of_lists[[cell_type]])) {
    p <- p + geom_signif(
      comparisons = list_of_lists[[cell_type]],
      map_signif_level = FALSE, 
      step_increase = 0.05
    )
  }
  
  # Add the plot to the list
  plot_list[[cell_type]] <- p
}

# Save each plot individually
ggsave("/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/HTseq_700/Figures/plot_mast_cell_resting.pdf", 
       plot = plot_list[["Mast.cells.resting"]], width = 8, height = 8, dpi = 300, bg = "white")

ggsave("/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/HTseq_700/Figures/plot_t_cells_cd8.pdf", 
       plot = plot_list[["T.cells.CD8"]], width = 8, height = 8, dpi = 300, bg = "white")

ggsave("/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/HTseq_700/Figures/plot_t_cells_follicular_helper.pdf", 
       plot = plot_list[["T.cells.follicular.helper"]], width = 8, height = 8, dpi = 300, bg = "white")

ggsave("/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/HTseq_700/Figures/plot_nk_cells_activated.pdf", 
       plot = plot_list[["NK.cells.activated"]], width = 8, height = 8, dpi = 300, bg = "white")

ggsave("/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/HTseq_700/Figures/plot_monocytes.pdf", 
       plot = plot_list[["Monocytes"]], width = 8, height = 8, dpi = 300, bg = "white")

ggsave("/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/HTseq_700/Figures/plot_macrophages_m2.pdf", 
       plot = plot_list[["Macrophages.M2"]], width = 8, height = 8, dpi = 300, bg = "white")
