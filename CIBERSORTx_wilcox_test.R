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
# =====================================================================
# CIBERSORTx: Cluster 11 vs All-Other Comparison (with p-value legend)
# =====================================================================

# 1) Libraries
library(ggplot2)
library(rstatix)
library(dplyr)
library(tidyr)
library(ggsignif)

# 2) Load data
cibersortx_results <- read.csv("/Users/kanatyalcin/Downloads/CIBERSORTx_Job31_Results.csv")
umap_clusters      <- read.csv("/Users/kanatyalcin/umap_cluster_info_1232.csv")

# 3) Preprocess
colnames(umap_clusters)[1] <- "Mixture"

merged_data <- merge(cibersortx_results, umap_clusters, by = "Mixture")

# 4) Reshape and recode to binary groups
long_data <- merged_data %>%
  select(Mixture, Cluster,
         `B.cells.naive`, `T.cells.CD8`, `T.cells.CD4.memory.resting`,
         `T.cells.follicular.helper`, `NK.cells.resting`, `NK.cells.activated`,
         `Monocytes`, `Macrophages.M2`, `Mast.cells.resting`, `Neutrophils`
  ) %>%
  pivot_longer(cols = -c(Mixture, Cluster),
               names_to = "Cell_Type", values_to = "Value") %>%
  mutate(
    ClusterGroup = factor(
      ifelse(Cluster == 11, "Cluster 11", "Rest"),
      levels = c("Rest", "Cluster 11")
    )
  )

# 5) Color palette for the two groups
binary_palette <- c("Rest" = "#FC8D62", "Cluster 11" = "#B3B3B3")

# 6) Combined facet plot with caption legend
combined_plot <- ggplot(long_data, aes(x = ClusterGroup, y = Value, fill = ClusterGroup)) +
  geom_boxplot() +
  facet_wrap(~Cell_Type, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = binary_palette) +
  labs(
    x       = "Cluster",
    y       = "Proportion",
    title   = "Cluster 11 vs All Rest",
    fill    = "Cluster",
    caption = "*: p < 0.05; **: p < 0.01; ***: p < 0.001"
  ) +
  theme_minimal() +
  theme(
    axis.text        = element_text(size = 30, face = "bold"),
    axis.title       = element_text(size = 40, face = "bold"),
    legend.title     = element_text(size = 35, face = "bold"),
    legend.text      = element_text(size = 30),
    legend.position  = "right",
    panel.grid       = element_blank(),
    plot.background  = element_rect(fill = "white", colour = "black", linewidth = 2),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 2),
    plot.title       = element_text(size = 48, face = "bold", hjust = 0.5),
    strip.text       = element_text(size = 30, face = "bold", hjust = 0),
    panel.spacing    = unit(2, "lines"),
    plot.margin      = unit(c(1,1,1,1), "cm"),
    plot.caption     = element_text(size = 20, hjust = 0, face = "italic")
  ) +
  geom_signif(
    comparisons      = list(c("Rest", "Cluster 11")),
    map_signif_level = TRUE,
    textsize         = 6,
    step_increase    = 0.05
  )

ggsave(
  "CIBERSORTx_cluster11_vs_rest_facet.png",
  combined_plot, width = 36, height = 26, dpi = 300
)

# 7) Individual plots per cell type (with same caption)
plot_list <- list()
for (ct in unique(long_data$Cell_Type)) {
  df_ct <- filter(long_data, Cell_Type == ct)
  p <- ggplot(df_ct, aes(x = ClusterGroup, y = Value, fill = ClusterGroup)) +
    geom_boxplot() +
    scale_fill_manual(values = binary_palette) +
    labs(
      x       = "Cluster",
      y       = "Proportion",
      title   = gsub("\\.", " ", ct),
      fill    = "Cluster",
      caption = "*: p < 0.05; **: p < 0.01; ***: p < 0.001"
    ) +
    theme_minimal() +
    theme(
      axis.text        = element_text(size = 30, face = "bold"),
      axis.title       = element_text(size = 35, face = "bold"),
      legend.title     = element_text(size = 35, face = "bold"),
      legend.text      = element_text(size = 30),
      legend.position  = "right",
      panel.grid       = element_blank(),
      plot.background  = element_rect(fill = "white", colour = "white", linewidth = 2),
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 2),
      plot.title       = element_text(size = 40, face = "bold", hjust = 0),
      plot.caption     = element_text(size = 16, hjust = 0, face = "italic")
    ) +
    geom_signif(
      comparisons      = list(c("Rest", "Cluster 11")),
      map_signif_level = TRUE,
      textsize         = 6,
      step_increase    = 0.05
    )
  plot_list[[ct]] <- p
}

# 8) Save a few example individual plots
ggsave("plot_mast_cells_resting_simplified_c8_vs_all.pdf",  plot_list[["Mast.cells.resting"]],  width=14, height=10, dpi=300)
ggsave("plot_t_cells_cd8_simplified_c8_vs_all.pdf",        plot_list[["T.cells.CD8"]],        width=14, height=10, dpi=300)
ggsave("plot_t_cells_follicular_helper_simplified_c8_vs_all.pdf", plot_list[["T.cells.follicular.helper"]], width=14, height=10, dpi=300)
ggsave("plot_nk_cells_activated_simplified_c8_vs_all.pdf", plot_list[["NK.cells.activated"]], width=14, height=10, dpi=300)
ggsave("plot_monocytes_simplified_c8_vs_all.pdf",          plot_list[["Monocytes"]],          width=14, height=10, dpi=300)
ggsave("plot_macrophages_m2_simplified_c8_vs_all.pdf",     plot_list[["Macrophages.M2"]],     width=14, height=10, dpi=300)
