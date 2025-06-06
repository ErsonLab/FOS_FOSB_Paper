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

library(xCell)
xcell_counts_tpm <- read.table("combat_counts_tpm.tsv", sep = "\t")
rownames(xcell_counts_tpm) <- gsub("\\.[0-9]*","",rownames(xcell_counts_tpm))
head(xcell_counts_tpm)
xcell_counts_tpm <- as.data.frame(xcell_counts_tpm)
xcell_counts_tpm$gene_symbol <- mapIds(EnsDb.Hsapiens.v86,
                                   keys=row.names(xcell_counts_tpm),
                                   column="SYMBOL",
                                   keytype="GENEID",
                                   multiVals="first")
xcell_counts_tpm <- xcell_counts_tpm %>%
  group_by(gene_symbol) %>%
  summarise(across(everything(), sum, na.rm = TRUE))
#xcell_counts <- xcell_counts_backup
xcell_counts_tpm <- xcell_counts_tpm[-nrow(xcell_counts_tpm), ]
rownames(xcell_counts_tpm) <- xcell_counts_tpm$gene_symbol
write.csv(xcell_counts_tpm, "xCell_input.csv")
xcell_counts_tpm <- read.csv("xCell_input.csv", row.names = 1)
xCell_res <- xCellAnalysis(xcell_counts_tpm)
xCell_res <- read.table("./Downloads/xcell_results_1232_combat_TPM_MTremoved.txt", sep ="\t", header = TRUE, row.names = 1)
# Load required libraries
library(reshape2)
library(dplyr)
library(ggplot2)  # Ensure ggplot2 is loaded

xcell_df <- xCell_res
xcell_df <- as.data.frame(xcell_df)
library(tibble)
# moves rownames into a column called "Cell_Type"
#xcell_df <- xcell_df %>%
#  rownames_to_column(var = "Cell_Type")
umap_df <- read.csv("umap_df_2d_12_clusters.csv")

# Remove Cluster 0 - this is noise
umap_df <- umap_df[umap_df$Cluster != 0, ]

# Rename columns for consistency
xcell_df["Cell_Type"] = row.names(xcell_df)
#colnames(xcell_df)[1] <- "Cell_Type"
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

comparison_results <- data.frame(Cell_Type = character(),
                                 p_value = numeric(),
                                 stringsAsFactors = FALSE)

cell_types <- unique(melted_df$Cell_Type)

cell_types

melted_df$Enrichment <- as.numeric(as.character(melted_df$Enrichment))

for (cell_type in cell_types) {
  # Subset the data for cell type
  cell_data <- melted_df[melted_df$Cell_Type == cell_type, ]
  
  cluster_2 <- cell_data[cell_data$Cluster == 11, "Enrichment"]
  all_other_clusters <- cell_data[cell_data$Cluster != 11, "Enrichment"]
  
  # Perform Wilcoxon rank-sum test
  test_result <- wilcox.test(cluster_2, all_other_clusters)
  
  # Store the results in the comparison_results data frame
  comparison_results <- rbind(comparison_results, 
                              data.frame(Cell_Type = cell_type, 
                                         p_value = test_result$p.value))
}

# Sort the results by p-value to find the most significant comparisons
comparison_results <- comparison_results %>% arrange(p_value)

write.csv(comparison_results, "xCELL_cluster11_vs_all_comparison_results.csv", row.names = FALSE)

# Top 15
top_15_cell_types <- head(comparison_results$Cell_Type, 15)
print(top_15_cell_types)

write.csv(top_15_cell_types, "xCELL_top15_cluster11_vs_all_comparison_results.csv", row.names = FALSE)

# Color palette that Kanat gave 
color_palette <- c(
  "0"  = "black",
  "1"  = "#CCA77E",
  "2"  = "#EE9174",
  "3"  = "#B29CB1",
  "4"  = "#AF9AC8",
  "5"  = "#E08CC4",
  "6"  = "#CDB490",
  "7"  = "#AED851",
  "8"  = "#E3D93E",
  "9"  = "#FAD450",
  "10" = "#EAC786",
  "11" = "#D1BDA1",
  "12" = "#B3B3B3"
)

library(ggsignif)

for (cell_type in top_15_cell_types) {
  # 1) subset and make sure Cluster is a factor
  df_sub <- subset(melted_df, Cell_Type == cell_type)
  df_sub$Cluster <- factor(df_sub$Cluster)
  
  # 2) count how many samples in each cluster
  counts <- table(df_sub$Cluster)
  
  # 3) build comparisons of "11" vs every other cluster,
  #    but only if both groups have >= 2 samples
  comp_list <- lapply(names(counts), function(cl) {
    if (cl != "11" && counts[["11"]] >= 2 && counts[[cl]] >= 2) {
      c("11", cl)
    } else {
      NULL
    }
  })
  comp_list <- Filter(Negate(is.null), comp_list)
  
  # 4) base ggplot
  p <- ggplot(df_sub,
              aes(x = Cluster, y = Enrichment, fill = Cluster)) +
    geom_boxplot() +
    scale_fill_manual(values = color_palette) +
    labs(
      x     = "Cluster",
      y     = "Proportion",
      title = gsub("\\.", " ", cell_type),
      fill  = "Cluster"
    ) +
    theme_minimal() +
    theme(
      axis.text       = element_text(size = 40, face = "bold"),
      axis.title      = element_text(size = 48, face = "bold"),
      legend.title    = element_blank(),
      legend.text     = element_text(size = 30),
      legend.position = "bottom",
      panel.grid      = element_blank(),
      plot.background = element_rect(fill = "white", colour = "black", linewidth = 2),
      panel.border    = element_rect(colour = "black", fill = NA, linewidth = 2),
      plot.title      = element_text(size = 48, face = "bold", hjust = 0)
    )
  
  # 5) only add the significance layer if we actually have comparisons
  if (length(comp_list) > 0) {
    p <- p + geom_signif(
      comparisons      = comp_list,
      map_signif_level = FALSE,
      test             = "wilcox.test",
      step_increase    = 0.05
    )
  }
  
  # 6) save to PDF
  file_name <- paste0(gsub("\\.", "_", cell_type), "_xCELL_boxplot.pdf")
  ggsave(
    filename = file_name,
    plot     = p,
    width    = 18,
    height   = 14,
    dpi      = 300,
    bg       = "white"
  )
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
  data_subset <- subset(melted_df, Cell_Type == cell_type)
  
  p <- ggplot(data_subset,
              aes(x = factor(Cluster), y = Enrichment, fill = factor(Cluster))) +
    geom_boxplot() +
    scale_fill_manual(values = color_palette) +
    labs(
      x     = "Cluster",
      y     = "Enrichment",
      title = gsub("\\.", " ", cell_type),
      fill  = "Cluster"
    ) +
    theme_minimal() +
    theme(
      axis.text       = element_text(size = 12, face = "bold"),
      axis.title      = element_text(size = 14, face = "bold"),
      panel.grid      = element_blank(),
      plot.title      = element_text(size = 16, face = "bold", hjust = 0),
      legend.position = "none"
    ) +
    geom_signif(
      comparisons     = list(c(11,1),c(11,2),c(11,3),c(11,4),c(11,5),c(11,6),c(11,7), c(11,8), c(11,10), c(12,11)),
      map_signif_level = FALSE,
      test             = "wilcox.test",
      step_increase    = 0.05
    )
  
  if (!legend_extracted) {
    legend <- get_legend(
      p + theme(
        legend.position = "bottom",
        legend.text     = element_text(size = 10),
        legend.title    = element_blank()
      )
    )
    legend_extracted <- TRUE
  }
  
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
  filename = "Combined_Top15_Boxplots.pdf",
  dpi = 300,
  width    = 30,
  height   = 40,
  bg = "white"
)

# Inform the user
cat("Combined plot saved as Combined_Top15_Boxplots.pdf\n")
# 1) Merge & melt exactly as before
merged_df <- merge(xcell_t, umap_df[, c("Sample", "Cluster")], by = "Sample")
melted_df <- melt(merged_df,
                  id.vars    = c("Sample", "Cluster"),
  plot = final_combined_plot,
  width = 24,    # Adjust width as needed
  height = 18,   # Adjust height as needed
                  variable.name = "Cell_Type",
                  value.name    = "Enrichment")

# 2) Compute your top-15 cell types
comparison_results <- data.frame(Cell_Type = character(),
                                 p_value   = numeric(),
                                 stringsAsFactors = FALSE)

for(ct in unique(melted_df$Cell_Type)) {
  cd   <- filter(melted_df, Cell_Type == ct)
  w   <- wilcox.test(
           cd$Enrichment[cd$Cluster == 11],
           cd$Enrichment[cd$Cluster != 11]
         )
  comparison_results <- bind_rows(
    comparison_results,
    tibble(Cell_Type = ct, p_value = w$p.value)
  )
}

top_15_cell_types <- comparison_results %>%
                      arrange(p_value) %>%
                      slice_head(n = 15) %>%
                      pull(Cell_Type)

# 3) Define a two-colour palette for “Rest” vs “Cluster 11”
group_palette <- c(
  Rest      = "grey80",
  `Cluster 11` = color_palette[["11"]]
)

# 4) Loop to save one PDF per cell type
for(ct in top_15_cell_types) {
  df_sub <- melted_df %>%
              filter(Cell_Type == ct) %>%
              mutate(
                Group = factor(
                  ifelse(Cluster == 11, "Cluster 11", "Rest"),
                  levels = c("Rest", "Cluster 11")
                )
              )

  p <- ggplot(df_sub, aes(x = Group, y = Enrichment, fill = Group)) +
    geom_boxplot() +
    scale_fill_manual(values = group_palette) +
    labs(
      x     = "Group",
      y     = "Proportion",
      title = gsub("\\.", " ", ct)
    ) +
    theme_minimal() +
    theme(
      axis.text       = element_text(size = 40, face = "bold"),
      axis.title      = element_text(size = 48, face = "bold"),
      legend.position = "none",
      panel.grid      = element_blank(),
      plot.background = element_rect(fill = "white", colour = "black", linewidth = 2),
      panel.border    = element_rect(colour = "black", fill = NA, linewidth = 2),
      plot.title      = element_text(size = 48, face = "bold", hjust = 0)
    ) +
    geom_signif(
      comparisons      = list(c("Rest", "Cluster 11")),
      map_signif_level = FALSE,
      test             = "wilcox.test",
      step_increase    = 0.05
    )

  fname <- paste0(gsub("\\.", "_", ct), "_xCELL_11_vs_Rest.pdf")
  ggsave(fname, p, width = 18, height = 14, dpi = 300, bg = "white")
  message("Saved: ", fname)
}

# 5) Now make the combined 3×5 grid with a shared legend
combined <- list()
for(ct in top_15_cell_types) {
  df_sub <- melted_df %>%
              filter(Cell_Type == ct) %>%
              mutate(
                Group = factor(
                  ifelse(Cluster == 11, "Cluster 11", "Rest"),
                  levels = c("Rest", "Cluster 11")
                )
              )

  p <- ggplot(df_sub, aes(x = Group, y = Enrichment, fill = Group)) +
    geom_boxplot() +
    scale_fill_manual(values = group_palette) +
    labs(title = gsub("\\.", " ", ct)) +
    theme_minimal() +
    theme(
      axis.text  = element_text(size = 12, face = "bold"),
      axis.title = element_blank(),
      plot.title = element_text(size = 16, face = "bold", hjust = 0),
      legend.position = "none"
    ) +
    geom_signif(
      comparisons      = list(c("Rest", "Cluster 11")),
      map_signif_level = FALSE,
      test             = "wilcox.test"
    )

  combined[[ct]] <- p
}

# extract a single legend
legend_plot <- combined[[1]] +
  theme(
    legend.position = "bottom",
    legend.text     = element_text(size = 10),
    legend.title    = element_blank()
  )
shared_legend <- get_legend(legend_plot)

# assemble grid + legend
grid <- plot_grid(plotlist = combined, ncol = 3, align = "hv")
final <- plot_grid(grid, shared_legend, ncol = 1, rel_heights = c(1, .1))

ggsave("Combined_11_vs_Rest_top15.pdf", final,
       width = 24, height = 30, dpi = 300, bg = "white")
message("Combined plot saved as Combined_11_vs_Rest_top15.pdf")





library(reshape2)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(cowplot)

# 1) Merge & melt exactly as before
merged_df <- merge(xcell_t, umap_df[, c("Sample", "Cluster")], by = "Sample")
melted_df <- melt(merged_df,
                  id.vars    = c("Sample", "Cluster"),
                  variable.name = "Cell_Type",
                  value.name    = "Enrichment")

# 2) Compute your top-15 cell types
comparison_results <- data.frame(Cell_Type = character(),
                                 p_value   = numeric(),
                                 stringsAsFactors = FALSE)
melted_df$Enrichment <- as.numeric(as.character(melted_df$Enrichment))

for(ct in unique(melted_df$Cell_Type)) {
  cd   <- filter(melted_df, Cell_Type == ct)
  w   <- wilcox.test(
    cd$Enrichment[cd$Cluster == 11],
    cd$Enrichment[cd$Cluster != 11]
  )
  comparison_results <- bind_rows(
    comparison_results,
    tibble(Cell_Type = ct, p_value = w$p.value)
  )
}

top_15_cell_types <- comparison_results %>%
  arrange(p_value) %>%
  slice_head(n = 15) %>%
  pull(Cell_Type)

# 3) Define a two-colour palette for “Rest” vs “Cluster 11”
group_palette <- c(
  Rest      = "grey80",
  `Cluster 11` = color_palette[["11"]]
)

# 4) Loop to save one PDF per cell type
for(ct in top_15_cell_types) {
  df_sub <- melted_df %>%
    filter(Cell_Type == ct) %>%
    mutate(
      Group = factor(
        ifelse(Cluster == 11, "Cluster 11", "Rest"),
        levels = c("Rest", "Cluster 11")
      )
    )
  
  p <- ggplot(df_sub, aes(x = Group, y = Enrichment, fill = Group)) +
    geom_boxplot() +
    scale_fill_manual(values = group_palette) +
    labs(
      x     = "Group",
      y     = "Proportion",
      title = gsub("\\.", " ", ct)
    ) +
    theme_minimal() +
    theme(
      axis.text       = element_text(size = 40, face = "bold"),
      axis.title      = element_text(size = 48, face = "bold"),
      legend.position = "none",
      panel.grid      = element_blank(),
      plot.background = element_rect(fill = "white", colour = "black", linewidth = 2),
      panel.border    = element_rect(colour = "black", fill = NA, linewidth = 2),
      plot.title      = element_text(size = 48, face = "bold", hjust = 0)
    ) +
    geom_signif(
      comparisons      = list(c("Rest", "Cluster 11")),
      map_signif_level = FALSE,
      test             = "wilcox.test",
      step_increase    = 0.05
    )
  
  fname <- paste0(gsub("\\.", "_", ct), "_xCELL_11_vs_Rest.pdf")
  ggsave(fname, p, width = 18, height = 14, dpi = 300, bg = "white")
  message("Saved: ", fname)
}

# 5) Now make the combined 3×5 grid with a shared legend
combined <- list()
for(ct in top_15_cell_types) {
  df_sub <- melted_df %>%
    filter(Cell_Type == ct) %>%
    mutate(
      Group = factor(
        ifelse(Cluster == 11, "Cluster 11", "Rest"),
        levels = c("Rest", "Cluster 11")
      )
    )
  
  p <- ggplot(df_sub, aes(x = Group, y = Enrichment, fill = Group)) +
    geom_boxplot() +
    scale_fill_manual(values = group_palette) +
    labs(title = gsub("\\.", " ", ct)) +
    theme_minimal() +
    theme(
      axis.text  = element_text(size = 12, face = "bold"),
      axis.title = element_blank(),
      panel.grid      = element_blank(),
      plot.title = element_text(size = 16, face = "bold", hjust = 0),
      legend.position = "none"
    ) +
    geom_signif(
      comparisons      = list(c("Rest", "Cluster 11")),
      map_signif_level = TRUE,
      test             = "wilcox.test"
    )
  
  combined[[ct]] <- p
}

# extract a single legend
legend_plot <- combined[[1]] +
  theme(
    legend.position = "bottom",
    legend.text     = element_text(size = 10),
    legend.title    = element_blank()
  )
shared_legend <- get_legend(legend_plot)

# assemble grid + legend
grid <- plot_grid(plotlist = combined, ncol = 3, align = "hv")
final <- plot_grid(grid, shared_legend, ncol = 1, rel_heights = c(1, .1))

ggsave("Combined_11_vs_Rest_top15.pdf", final,
       width = 24, height = 18, dpi = 300, bg = "white")
message("Combined plot saved as Combined_11_vs_Rest_top15.pdf")
