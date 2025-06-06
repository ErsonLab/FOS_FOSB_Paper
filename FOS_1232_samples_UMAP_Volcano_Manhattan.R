# ===================================================================================================================================
#
# Paper Name: Novel Subgroup of Meningiomas Involving FOS and FOSB Gene Fusions
#
# ===================================================================================================================================
# 
# Description: This code facilitates the bulk RNA-seq analysis performed on meningioma samples to uncover novel genomic and transcriptomic 
# features in uncharacterized tumor subgroups. The workflow begins with preprocessing, including quality control, adapter trimming, and 
# alignment of RNA-seq reads to a reference genome, followed by gene expression quantification. Differential expression analysis identifies 
# key genes and pathways, focusing on AP-1 complex targets and immune cell activity in tumors with FOS/FOSB fusion Functional enrichment 
# analysis and clustering highlight unique transcriptomic profiles associated with the novel FOS/FOSB fusion subgroup. 
# Additionally, the code generates high-quality visualizations of differential expression, enriched pathways etc., providing insights into 
# the biology of this meningioma subtype and potential therapeutic targets.
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

#########################################################
### Step 0 : Import Libraries
#########################################################
library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library("org.Hs.eg.db")
library(EnhancedVolcano)
library(textshaping)
library(RColorBrewer)
library(circlize)
library(clusterProfiler)
library(AnnotationDbi)
library(enrichR)
library(ggsci)
library(readxl)
library(Matrix)
library(SummarizedExperiment)
library(factoextra)
library(cowplot)
library(tidyr)
library(biomaRt)
library(dplyr)
library(ggcorrplot)
library(GGally)
library(PCAtools)
library(FactoMineR)
library(gplots)
library(vsn)
library(gridExtra)
library(edgeR)
library(limma)
library(Glimma)
library(sva)
library(tximport)
library(Homo.sapiens)
library(EnsDb.Hsapiens.v86)
library(flashClust)
library(WGCNA)
library(ggalt)
library(ggnewscale)
library(grDevices)
library(viridis)
library(pheatmap)
library(plotly)
library(htmlwidgets)
library(dbscan)
library(stringr)
library(umap)
library(ggrepel)
library(ggsignif)
library(aod)

#########################################################
### Step 1 : Data Import and Pre-processing
#########################################################
set.seed(123)

counts <- read.csv("merged_counts_1232.csv",
                   header = TRUE, 
                   row.names = 1, 
                   sep = ",",
                   check.names = FALSE)


sampleinfo <- read.delim("s_info_1232.txt",
                         header = TRUE, sep = "\t",check.names = FALSE)



# Display the updated sampleinfo DataFrame
print(sampleinfo)

nf2_fusion_cases <- c(
  "SRR10042744",
  "SRR10042745",
  "SRR10042753",
  "SRR21416196",
  "SRR21416257",
  "SRR21416297",
  "SRR21416302",
  "SRR21416387",
  "SRR21416408",
  "SRR21416444",
  "SRR27388530",
  "SRR27388603",
  "SRR27388613",
  "CAM20.0050",
  "CAM20.0080",
  "CAM20.0085",
  "CAM20.0106",
  "CAM20.0155"
)

yap1_maml2_cases <- c(
  "SRR10042672"
)

q235loss <- c("ERR12048884",
              "ERR12048896",
              "ERR12048897",
              "ERR12048891",
              "SRR10042736",
              "SRR21416219",
              "SRR21416235",
              "SRR21416268",
              "SRR21416425",
              "SRR21416461",
              "SRR27388475",
              "SRR27388559",
              "SRR27388612",
              "SRR27388693")

q736loss <- c("SRR10042651",
              "SRR10042713",
              "SRR10042733",
              "SRR21416216",
              "SRR21416225")

target_batches <- c("YALE", "UCSF", "PALACKY", "BAYLOR", "FHCC", "HKU/UCSF")

target_batches <- c("YALE", "UCSF (2022)", "UCSF (2020)", "UCSF (2018)", "PALACKY", "BAYLOR", "FHCC", "UCSD", "HKU/UCSF", "UoT")

# Define the target batches

# Use subset() to filter sampleinfo for rows with desired Batch values
subset_sampleinfo <- subset(sampleinfo, Batch %in% target_batches)

# Subset counts using the Sample column from the filtered sampleinfo
subset_counts <- counts[, as.character(subset_sampleinfo$Sample)]


counts <- subset_counts
sampleinfo <- subset_sampleinfo

determine_fusion_status <- function(srr, fos_status) {
  # Check if FOS_Fusion_Status is "Yes"
  if (fos_status == "Yes") {
    return("FOS")
  }
  
  if (srr %in% nf2_fusion_cases) {
    return("NF2_FUSION")
  }
  
  if (srr %in% q235loss) {
    return("2Q35_FUSION")
  }
  if (srr %in% q736loss) {
    return("7Q36.3_FUSION")
  }
  
  
  if (srr %in% yap1_maml2_cases) {
    return("YAP1_MAML2")
  }
  return("NA")
}

sampleinfo$FUSION_STATUS <- mapply(determine_fusion_status, sampleinfo$Sample, sampleinfo$FOS_Fusion_Status)

sampleinfo$Sample <- as.factor(sampleinfo$Sample)
sampleinfo$FOS_Fusion_Status <- as.factor(sampleinfo$FOS_Fusion_Status)
sampleinfo$Batch <- as.factor(sampleinfo$Batch)
sampleinfo$FUSION_STATUS <- as.factor(sampleinfo$FUSION_STATUS)

sampleinfo$GRADE <- as.factor(sampleinfo$Grade)
sampleinfo$Age <- as.numeric(sampleinfo$Age)
sampleinfo$Mutation <- as.factor(sampleinfo$Mutation)
sampleinfo$Gender <- as.factor(sampleinfo$Gender)

sampleinfo$Baylor_New <- as.factor(sampleinfo$Baylor_New)
sampleinfo$UCSF_Methylation <- as.factor(sampleinfo$UCSF_Methylation)
sampleinfo$Toronto_Methylation <- as.factor(sampleinfo$Toronto_Methylation)
sampleinfo$Merged_Methylation <- as.factor(sampleinfo$Merged_Methylation)

#########################################################
### Step 2 : Batch Correction and Normalization
#########################################################

# Perform Combat-Seq batch correction on raw counts
combat_counts <- ComBat_seq(as.matrix(counts), batch=sampleinfo$Batch, group=NULL)


combat_counts_vst <- vst(round(combat_counts), blind = FALSE)

#########################################################
### Step 3 : UMAP ANALYSIS AND PLOT GENERATION
#########################################################

# Run PCA on the Combat-Seq corrected VST data
pca_combat_first_then_VST <- prcomp(t(combat_counts_vst))


# Extract the PCA components (scores)
pca_data <- pca_combat_first_then_VST$x

custom.config <- umap.defaults
custom.config$random_state <- 123
umap_result_2d <- umap(pca_data, n_neighbors = 13, min_dist = 0.01, n_components = 2, config = custom.config)


umap_coordinates_2d <- umap_result_2d$layout
umap_coordinates_2d <- na.omit(umap_coordinates_2d)  # Remove NA rows
umap_coordinates_2d <- umap_coordinates_2d[apply(umap_coordinates_2d, 1, function(row) all(is.finite(row))), ]  # Remove rows with Inf


dbscan_result_2d <- dbscan(umap_coordinates_2d, eps = 0.28, minPts = 8) ## min pts = 8 

# Access the cluster labels
cluster_labels_2d <- dbscan_result_2d$cluster

umap_df_2d <- as.data.frame(umap_coordinates_2d)
colnames(umap_df_2d) <- c("UMAP_1", "UMAP_2")

# Add cluster labels to the data frame
umap_df_2d$Cluster <- as.factor(cluster_labels_2d)

umap_df_2d$FOS_Fusion_Status <- as.factor(sampleinfo$FOS_Fusion_Status)

umap_df_2d$Sample_Name <- sampleinfo$Sample  # Add sample names

umap_df_2d$Batch <- sampleinfo$Batch

umap_df_2d$FUSION_STATUS <- sampleinfo$FUSION_STATUS

umap_df_2d$GRADE <- sampleinfo$GRADE

umap_df_2d$AGE <- sampleinfo$Age

umap_df_2d$GENDER <- sampleinfo$Gender

umap_df_2d$MUTATION <- sampleinfo$Mutation

umap_df_2d$BAYLOR_METHYLATION <- sampleinfo$Baylor_New

umap_df_2d$UCSF_METHYLATION <- sampleinfo$UCSF_Methylation

umap_df_2d$TORONTO_METHYLATION <- sampleinfo$Toronto_Methylation

umap_df_2d$MERGED_METHYLATION <- sampleinfo$Merged_Methylation

umap_df_2d_noise_removed <- subset(umap_df_2d, Cluster != 0)

umap_df_2d_noise_removed <- umap_df_2d_noise_removed %>%
  mutate(
    MUTATION = as.character(MUTATION),   # ensure it's not a factor
    MUTATION = recode(
      MUTATION,
      "2Q35_FUSION"    = "N/A",
      "7Q36.3_FUSION"  = "N/A",
      "NF2 + TRAF7"    = "NF2",
      "NF2_FUSION"     = "N/A",
      "YaleUnknown"    = "N/A",
      "YAP1_MAML2"     = "N/A",
      .default = MUTATION                   # keep all other values
    ),
    MUTATION = factor(MUTATION)            # optional: back to factor
  )


umap_df_2d_noise_removed <- umap_df_2d_noise_removed %>%
  mutate(
    GENDER = as.character(GENDER),   # ensure it's not a factor
    GENDER = recode(
      GENDER,
      "Female"    = "F",
      "Male"  = "M",
      "Unknown"     = "N/A",
      .default = GENDER                   # keep all other values
    ),
    GENDER = factor(GENDER)            # optional: back to factor
  )

batch_labels <- umap_df_2d_noise_removed$Cluster

clusters <- sort(unique(batch_labels))

# Number of clusters
n_clusters <- length(clusters)
colors <- brewer.pal(n_clusters, "Set2")

two_dplot_dbscan <- plot_ly(umap_df_2d_noise_removed,
                            x = ~UMAP_1, y = ~UMAP_2, 
                            color = ~Cluster, colors = colors, 
                            #symbol = ~Cluster, symbols = symbol_set, # Adjust symbols as needed
                            type = "scatter", mode = "markers",
                            marker = list(size = 10),
                            text = ~paste("Sample: ", Sample_Name, "<br>Mutation Status: ", MUTATION , "<br>Cluster: ", Cluster, "<br>Sequencing Center: ", Batch),  # Concatenate sample info for hover
                            hoverinfo = "text") %>%  # Ensure hover text is enabled
  layout(
    title   = "2D UMAP Plot with DBSCAN Clusters, FOS Fusion Status, and Sample Information",
    xaxis   = list(
      title    = "UMAP 1",
      showgrid = FALSE,
      zeroline = FALSE
    ),
    yaxis   = list(
      title    = "UMAP 2",
      showgrid = FALSE,
      zeroline = FALSE
    ),
    showlegend = TRUE
  )

two_dplot_dbscan
save_image(two_dplot_dbscan, "2D_Umap_cluster.pdf", width = 1250, height = 1028)

library(dplyr)
library(plotly)

# 1) Define the exact draw‐order of statuses:
statuses <- c(
  "NA", 
  "2Q35_FUSION", 
  "7Q36.3_FUSION", 
  "FOS", 
  "NF2_FUSION", 
  "YAP1_MAML2"
)

# 2) Make a named color map (choose any vibrant hexes you prefer):
status_colors <- c(
  "NA"          = "lightgray",
  "2Q35_FUSION" = "#17BECF",
  "7Q36.3_FUSION"= "#9467BD",
  "FOS"         = "#E41A1C",
  "NF2_FUSION"  = "orange",
  "YAP1_MAML2"  = "#4DAF4A"
)

# 3) Build the plot, one trace per status in that exact order:
p <- plot_ly()

for (st in statuses) {
  df_sub <- umap_df_2d_noise_removed %>%
    filter(FUSION_STATUS == st)
  
  # optional: make NA semi‐transparent
  op <- if (st == "NA") 0.4 else 1
  
  p <- p %>%
    add_trace(
      data      = df_sub,
      x         = ~UMAP_1, y = ~UMAP_2,
      type      = "scatter",
      mode      = "markers",
      marker    = list(
        size    = 10,
        color   = status_colors[st],
        opacity = op
      ),
      name      = st,
      hoverinfo = "text",
      text      = ~paste("Sample:", Sample_Name, "<br>Status:", MUTATION)
    )
}

# 4) Final layout tweaks
p <- p %>% 
  layout(
    title   = "2D UMAP Plot with DBSCAN Clusters, FOS Fusion Status, and Sample Information",
    xaxis   = list(
      title    = "UMAP 1",
      showgrid = FALSE,
      zeroline = FALSE
    ),
    yaxis   = list(
      title    = "UMAP 2",
      showgrid = FALSE,
      zeroline = FALSE
    ),
    showlegend = TRUE
  )

p

two_dplot_dbscan <- p
save_image(two_dplot_dbscan, "2D_Umap_fusion_status.pdf", width = 1250, height = 1028)

#########################################################################
### Step 4 : DIFFERENTIAL GENE EXPRESSION AND ENRICHMENT ANALYSIS
#########################################################################

#### Pairwise DESeq ### 

other_clusters <- levels(umap_df_2d$Cluster)[levels(umap_df_2d$Cluster) != "11"]

# Create a list to store results for each comparison
pairwise_results <- list()

# Function to convert gene IDs to gene symbols and filter by protein-coding
convert_and_filter_protein_coding <- function(res) {
  # Remove version numbers from gene IDs
  rownames(res) <- gsub("\\.[0-9]*", "", rownames(res))
  
  # Map gene symbols and biotypes
  res$gene_symbol <- mapIds(EnsDb.Hsapiens.v86,
                            keys = row.names(res),
                            column = "SYMBOL",
                            keytype = "GENEID",
                            multiVals = "first")
  
  res$biotype <- mapIds(EnsDb.Hsapiens.v86,
                        keys = row.names(res),
                        column = "GENEBIOTYPE",
                        keytype = "GENEID",
                        multiVals = "first")
  
  # Filter for protein-coding genes and remove NAs
  res_proteincoding <- dplyr::filter(as.data.frame(res), biotype == "protein_coding")
  res_proteincoding <- na.omit(res_proteincoding)
  
  return(res_proteincoding)
}

convert_but_dont_filter_protein_coding <- function(res) {
  # Remove version numbers from gene IDs
  rownames(res) <- gsub("\\.[0-9]*", "", rownames(res))
  
  # Map gene symbols and biotypes
  res$gene_symbol <- mapIds(EnsDb.Hsapiens.v86,
                            keys = row.names(res),
                            column = "SYMBOL",
                            keytype = "GENEID",
                            multiVals = "first")
  
  res <- na.omit(res)
  
  return(res)
}

# Prepare metadata (sample info), labeling "Cluster_2" vs the current comparison cluster
metadata <- data.frame(
  row.names = colnames(combat_counts),  # Ensure that row names match sample IDs in the count matrix
  cluster = umap_df_2d$Cluster  # Use the cluster info from your UMAP data
)

# Subset data to only include samples from Cluster 2 and the comparison cluster
dds_subset <- DESeqDataSetFromMatrix(countData = combat_counts, 
                                     colData = metadata, 
                                     design = ~ cluster)

# Run DESeq2
dds_subset <- DESeq(dds_subset)

# Create a loop to compare Cluster 2 with each other cluster one by one
for (cluster in other_clusters) {
  

  # Extract results comparing Cluster 2 vs the comparison cluster
  res_c2_vs_others <- results(dds_subset, contrast = c("cluster", "11", cluster))
  
  # Convert gene IDs to gene symbols and filter for protein-coding genes
  res_proteincoding <- convert_and_filter_protein_coding(res_c2_vs_others)
  
  res_raw <- convert_but_dont_filter_protein_coding(res_c2_vs_others)
  
  # Save raw protein-coding filtered results (no padj or log2fold filtering)
  pairwise_results[[paste0("Cluster_11_vs_Cluster_", cluster, "_raw_protein_coding")]] <- res_proteincoding
  
  pairwise_results[[paste0("Cluster_11_vs_Cluster_", cluster, "_raw")]] <- res_raw
  
  # Filter for significant genes based on the different thresholds
  
  # Filtering for padj < 0.05 and abs(log2FoldChange) > 5
  res_filtered_abs_5 <- res_proteincoding[which(res_proteincoding$padj < 0.05 & abs(res_proteincoding$log2FoldChange) > 5), ]
  
  # Filtering for log2FoldChange > 5
  res_filtered_log2fc_5 <- res_proteincoding[which(res_proteincoding$padj < 0.05 & res_proteincoding$log2FoldChange > 5), ]
  
  # Filtering for log2FoldChange > 2
  res_filtered_log2fc_2 <- res_proteincoding[which(res_proteincoding$padj < 0.05 & res_proteincoding$log2FoldChange > 2), ]
  
  # Filtering for padj < 0.05 and abs(log2FoldChange) > 2
  res_filtered_abs_2 <- res_proteincoding[which(res_proteincoding$padj < 0.05 & abs(res_proteincoding$log2FoldChange) > 2), ]
  
  # Save filtered results to the list with a unique name for each filter condition
  pairwise_results[[paste0("Cluster_11_vs_Cluster_", cluster, "_padj_abs_log2fc_5")]] <- res_filtered_abs_5
  pairwise_results[[paste0("Cluster_11_vs_Cluster_", cluster, "_padj_log2fc_5")]] <- res_filtered_log2fc_5
  pairwise_results[[paste0("Cluster_11_vs_Cluster_", cluster, "_padj_log2fc_2")]] <- res_filtered_log2fc_2
  pairwise_results[[paste0("Cluster_11_vs_Cluster_", cluster, "_padj_abs_log2fc_2")]] <- res_filtered_abs_2
  
  # Optionally print a message to track progress
  print(paste("Stored results for Cluster 11 vs Cluster", cluster))
}

# Additionally, compare Cluster 2 vs all other clusters combined
metadata_all <- data.frame(
  row.names = colnames(combat_counts), 
  cluster = umap_df_2d$Cluster
)

# Create a new column marking Cluster 2 as "Cluster_2" and all others as "Other"
metadata_all$group <- ifelse(metadata_all$cluster == "11", "Cluster_11", "Other")

# Subset data to only include samples from Cluster 2 and all other clusters
dds_all <- DESeqDataSetFromMatrix(countData = combat_counts[, !is.na(metadata_all$group)], 
                                  colData = metadata_all[!is.na(metadata_all$group), ], 
                                  design = ~ group)

# Run DESeq2
dds_all <- DESeq(dds_all)

# Extract results comparing Cluster 2 vs all others
res_c2_vs_all <- results(dds_all, contrast = c("group", "Cluster_11", "Other"))

# Convert gene IDs to gene symbols and filter for protein-coding genes
res_proteincoding_all <- convert_and_filter_protein_coding(res_c2_vs_all)

res_raw_all <- convert_but_dont_filter_protein_coding(res_c2_vs_all)

# Save raw protein-coding filtered results (no padj or log2fold filtering) for Cluster 2 vs All
pairwise_results[["Cluster_11_vs_All_raw_protein_coding"]] <- res_proteincoding_all

pairwise_results[["Cluster_11_vs_All_raw"]] <- res_raw_all

# Apply the same filtering as above
res_filtered_all_abs_5 <- res_proteincoding_all[which(res_proteincoding_all$padj < 0.05 & abs(res_proteincoding_all$log2FoldChange) > 5), ]
res_filtered_all_log2fc_5 <- res_proteincoding_all[which(res_proteincoding_all$padj < 0.05 & res_proteincoding_all$log2FoldChange > 5), ]
res_filtered_all_log2fc_2 <- res_proteincoding_all[which(res_proteincoding_all$padj < 0.05 & res_proteincoding_all$log2FoldChange > 2), ]
res_filtered_all_abs_2 <- res_proteincoding_all[which(res_proteincoding_all$padj < 0.05 & abs(res_proteincoding_all$log2FoldChange) > 2), ]

# Save filtered results to the list with unique names for Cluster 2 vs All
pairwise_results[["Cluster_11_vs_All_padj_abs_log2fc_5"]] <- res_filtered_all_abs_5
pairwise_results[["Cluster_11_vs_All_padj_log2fc_5"]] <- res_filtered_all_log2fc_5
pairwise_results[["Cluster_11_vs_All_padj_log2fc_2"]] <- res_filtered_all_log2fc_2
pairwise_results[["Cluster_11_vs_All_padj_abs_log2fc_2"]] <- res_filtered_all_abs_2

print("Stored results for Cluster 11 vs All")

#############################################################################################################################
### Pairwise Enrichment ### 

# List of databases for biological processes and KEGG pathways
selected_dbs_bp_kegg <- c("GO_Biological_Process_2021", "KEGG_2021_Human", "Reactome_2016")

# List of databases for transcription factors
selected_dbs_tf <- c("TRRUST_Transcription_Factors_2019", "ChEA_2016", "TRANSFAC_and_JASPAR_PWMs")

# Create a folder to save plots if it doesn't exist
#output_dir <- "/Users/kanatyalcin/Downloads/Enrichment_Pairwise"
#if (!dir.exists(output_dir)) {
#  dir.create(output_dir)
#}

# Loop through each comparison in pairwise_results
for (comparison_name in names(pairwise_results)) {
  if (grepl("_raw_protein_coding", comparison_name)) {
    next
  }
  if (grepl("_raw", comparison_name)) {
    next
  }
  # Get the filtered genes (gene symbols) for the current comparison
  filtered_genes <- pairwise_results[[comparison_name]]$gene_symbol
  
  # If no genes, skip this iteration
  if (length(filtered_genes) == 0) {
    next
  }
  
  # Run enrichment analysis on Biological Processes and KEGG
  enrich_results <- enrichr(filtered_genes, selected_dbs_bp_kegg)
  
  # Extract the GO and KEGG results
  go_results <- enrich_results[["GO_Biological_Process_2021"]]
  kegg_results <- enrich_results[["KEGG_2021_Human"]]
  reactome_results <- enrich_results[["Reactome_2016"]]
  
  # Select top 10 results for GO and KEGG
  top_go <- go_results %>% top_n(-10, P.value)
  top_kegg <- kegg_results %>% top_n(-10, P.value)
  top_reactome <- reactome_results %>% top_n(-10, P.value)
  
  # Combine GO and KEGG results
  top_enrich <- bind_rows(
    top_go %>% mutate(Database = "GO"),
    top_kegg %>% mutate(Database = "KEGG"),
    top_reactome %>% mutate(Database = "Reactome")
  )
  
  # Create the bar plot for GO and KEGG
  enrichR_top <- ggplot(top_enrich, aes(x = reorder(Term, -log10(P.value)), y = -log10(P.value), fill = Database)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(
      title = paste0("Top Enriched GO Terms, KEGG and Reactome Pathways: ", comparison_name),
      x = "Terms",
      y = "-log10(p-value)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    )
  
  # Run enrichment analysis for transcription factor databases
  tf_enrich_results <- enrichr(filtered_genes, selected_dbs_tf)
  
  # Extract the TRRUST, ChEA, and TRANSFAC results
  trrust_results <- tf_enrich_results[["TRRUST_Transcription_Factors_2019"]]
  chea_results <- tf_enrich_results[["ChEA_2016"]]
  transfac_results <- tf_enrich_results[["TRANSFAC_and_JASPAR_PWMs"]]
  
  # Filter for human-specific transcription factors
  trrust_human <- trrust_results %>% dplyr::filter(str_detect(Term, regex("human", ignore_case = TRUE)))
  chea_human <- chea_results %>% dplyr::filter(str_detect(Term, regex("human", ignore_case = TRUE)))
  transfac_human <- transfac_results %>% dplyr::filter(str_detect(Term, regex("human", ignore_case = TRUE)))
  
  # Select top 10 results for TRRUST, ChEA, and TRANSFAC
  top_trrust <- trrust_human %>% top_n(-10, P.value)
  top_chea <- chea_human %>% top_n(-10, P.value)
  top_transfac <- transfac_human %>% top_n(-10, P.value)
  
  # Combine results from TRRUST, ChEA, and TRANSFAC
  top_tf_enrich <- bind_rows(
    top_trrust %>% mutate(Database = "TRRUST"),
    top_chea %>% mutate(Database = "ChEA"),
    top_transfac %>% mutate(Database = "TRANSFAC_and_JASPAR_PWMs")
  )
  
  # Create the bar plot for transcription factor enrichment
  enrichR_top_TF <- ggplot(top_tf_enrich, aes(x = reorder(Term, -log10(P.value)), y = -log10(P.value), fill = Database)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(
      title = paste0("Top Enriched Human Transcription Factors: ", comparison_name),
      x = "Transcription Factors",
      y = "-log10(p-value)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    )
  
  # Save the plot for transcription factor enrichment
  ggsave(filename = paste0(comparison_name, "_TF.pdf"), plot = enrichR_top_TF, width = 20, height = 8)
  
  print(paste("Enrichment analysis completed for", comparison_name))
}

print("All enrichment analyses and plots have been saved.")


###############################################################################
### Step 5 : VOLCANO PLOT GENERATION BASED ON DIFFERENTIAL EXPRESSION DATA
###############################################################################


# Convert your data to a data frame
data_df <- as.data.frame(pairwise_results[["Cluster_11_vs_All_raw_protein_coding"]])

# Identify the top 10 most upregulated significant genes
top_upregulated <- data_df %>%
  dplyr::filter(padj < 0.05 & log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)

# Identify the top 10 most downregulated significant genes
top_downregulated <- data_df %>%
  dplyr::filter(padj < 0.05 & log2FoldChange < 0) %>%
  arrange(log2FoldChange) %>%
  head(10)

# Combine the top upregulated and downregulated genes
top_genes <- bind_rows(top_upregulated, top_downregulated)

# Create the volcano plot
volcano_plot_top20 <- ggplot(data_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05 & log2FoldChange > 2, "Upregulated",
                                ifelse(padj < 0.05 & log2FoldChange < -2, "Downregulated", "Not significant"))),
             alpha = 0.7, size = 3) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "grey")) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 40, face = "bold"),      # Increased tick labels font size
    axis.title = element_text(size = 48, face = "bold"),     # Increased axis labels font size
    legend.title = element_blank(),
    legend.text = element_text(size = 30),
    legend.position = "bottom",
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", colour = "black", linewidth = 2),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    plot.title = element_text(size = 48, face = "bold", hjust = 0.5)  # Increased plot title font size
  ) +
  labs(x = "log2(Fold Change)", y = "-log10(p-value)",
       title = "Volcano Plot of Differential Expression Analysis (Top 10 Up and Down Regulated Genes)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = c("blue", "red")) +
  ggrepel::geom_label_repel(
    data = top_genes,
    aes(label = gene_symbol),
    size = 12,               # Increased gene labels font size
    fontface = "bold",       # Gene labels bold
    box.padding = 0.5,
    fill = "white",
    color = "black",
    segment.color = 'grey50',
    max.overlaps = Inf
  ) +
  scale_y_continuous(breaks = seq(0, 35, 5), limits = c(0, 35)) +
  scale_x_continuous(breaks = seq(-25, 25, 5), limits = c(-25, 25))

# Display the plot
ggsave("Volcano_Plot_padj005_log2FC2_proteincoding.png", volcano_plot_top20, width = 29, height = 16, dpi = 300)

########################################################################
### Step 5 : MANHATTAN PLOT GENERATION BASED ON ENRICHMENT ANALYSIS
########################################################################

create_manhattan_plot_v12 <- function(enrich_results, comparison_name, output_dir) {
  
  # For Nature color palette
  
  # Extract and combine results from GO, KEGG, TRRUST, ChEA, and TRANSFAC
  go_results <- enrich_results[["GO_Biological_Process_2021"]]
  kegg_results <- enrich_results[["KEGG_2021_Human"]]
  trrust_results <- enrich_results[["TRRUST_Transcription_Factors_2019"]]
  chea_results <- enrich_results[["ChEA_2016"]]
  transfac_results <- enrich_results[["TRANSFAC_and_JASPAR_PWMs"]]
  
  # Apply filters for transcription factors to keep only "human" related terms
  trrust_human <- trrust_results %>% filter(str_detect(Term, regex("human", ignore_case = TRUE)))
  chea_human <- chea_results %>% filter(str_detect(Term, regex("human", ignore_case = TRUE)))
  transfac_human <- transfac_results %>% filter(str_detect(Term, regex("human", ignore_case = TRUE)))
  
  # Select top 10 results for TRRUST, ChEA, and TRANSFAC
  top_trrust <- trrust_human %>% top_n(-10, P.value)
  top_chea <- chea_human %>% top_n(-10, P.value)
  top_transfac <- transfac_human %>% top_n(-10, P.value)
  
  # Combine filtered results from transcription factors
  top_tf_enrich <- bind_rows(
    top_trrust %>% mutate(Database = "TF:TRRUST"),
    top_chea %>% mutate(Database = "TF:ChEA"),
    top_transfac %>% mutate(Database = "TF:TRANSFAC")
  )
  
  # Combine results from GO and KEGG along with the filtered transcription factors
  combined_results <- bind_rows(
    go_results %>% mutate(Database = "GO:BP"),
    kegg_results %>% mutate(Database = "KEGG"),
    top_tf_enrich  # Use the filtered TF results
  )
  
  # Add a -log10(p-adj) column for plotting (if p-adj is not available, use P.value)
  combined_results <- combined_results %>%
    mutate(`-log10(P.value)` = -log10(P.value))
  
  # Remove NA or infinite values
  combined_results <- combined_results %>%
    filter(!is.na(`-log10(P.value)`) & is.finite(`-log10(P.value)`))
  
  # Define the colors for each database using the Nature color palette
  # Install and load ggsci if necessary
  if (!requireNamespace("ggsci", quietly = TRUE)) {
    install.packages("ggsci")
  }
  
  # Get the unique database names
  database_names <- unique(combined_results$Database)
  
  # Generate colors using the Nature palette
  n_colors <- pal_npg()(length(database_names))
  
  # Assign colors to databases
  database_colors <- setNames(n_colors, database_names)
  
  # Convert Database to factor with specified levels and numeric x positions
  combined_results <- combined_results %>% 
    mutate(Database = factor(Database, levels = database_names)) %>% 
    mutate(x_pos = as.numeric(Database))
  
  # Add jittered x positions
  set.seed(123)  # For reproducibility
  combined_results <- combined_results %>% 
    mutate(jittered_x = x_pos + runif(n(), min = -0.3, max = 0.3))
  
  # Create a mapping of x positions to Database names
  x_labels <- combined_results %>% 
    select(x_pos, Database) %>% 
    distinct() %>% 
    arrange(x_pos)
  
  # Labeling logic:
  # 1. Pathways (GO, KEGG): Top 2 most significant
  pathways_top2 <- combined_results %>%
    filter(Database %in% c("GO:BP", "KEGG")) %>%
    group_by(Database) %>%
    top_n(-2, P.value)
  
  pathways_exact_match <- combined_results %>%
    filter(Database %in% c("GO:BP", "KEGG")) %>%
    filter(str_to_lower(Term) == "regulation of response to wounding")
  
  pathways_to_label <- bind_rows(pathways_top2, pathways_exact_match) %>%
    distinct()  # Remove duplicates if any
  
  # 2. Transcription Factors: Label anything with -log10(p-adj) > 3.2 (p-adj < 0.001)
  label_threshold <- 3 # Corresponds to p-adj < 0.001
  tf_to_label <- combined_results %>%
    filter(str_detect(Database, "TF") & `-log10(P.value)` > label_threshold) %>%  # Focus on transcription factors
    mutate(Term = gsub(" .*", "", Term))  # Remove everything after the first space for cleaner labels
  
  # Create the Manhattan-like plot
  manhattan_plot <- ggplot(combined_results, aes(x = jittered_x, y = `-log10(P.value)`, color = Database)) +
    # Add colored tiles (bars) for each database along the x-axis using geom_tile
    geom_tile(data = x_labels, aes(x = x_pos, y = 0, fill = Database), height = 0.05, width = 1, show.legend = FALSE) +
    scale_fill_manual(values = database_colors) +  # Assign colors for each database (bars)
    scale_color_manual(values = database_colors) +  # Ensure dots use the same colors as the bars
    
    geom_point(aes(size = `-log10(P.value)`), alpha = 0.7) +  # Plot the points using jittered x positions
    geom_hline(yintercept = label_threshold, color = "red", linetype = "dashed") +  # Add significance line
    
    # Add labels for the top pathways with reversed arrows
    geom_label_repel(
      data = pathways_to_label, 
      aes(x = jittered_x, label = Term), 
      size = 3.5,
      max.overlaps = Inf, 
      arrow = arrow(length = unit(0.01, "npc"), ends = "first", type = "closed"),
      box.padding = 0.5, 
      point.padding = 0.5, 
      segment.color = "grey50"
    ) +
    
    # Add labels for transcription factors that meet the threshold with reversed arrows
    geom_label_repel(
      data = tf_to_label, 
      aes(x = jittered_x, label = Term), 
      size = 3.5,
      max.overlaps = Inf, 
      arrow = arrow(length = unit(0.01, "npc"), ends = "first", type = "closed"),
      box.padding = 0.5, 
      point.padding = 0.5, 
      segment.color = "grey50"
    ) +
    
    labs(
      title = paste0("Manhattan Plot for ", comparison_name),
      x = "Databases",
      y = "-log10(p-adj)"
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",  # Remove legend
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(size = 16, face = "bold"),  # Increase size and make bold
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),  # Adjust x-axis labels
      axis.text.y = element_text(size = 14, face = "bold"),  # Adjust y-axis labels
      panel.background = element_rect(fill = "white"),     # Set background to white
      plot.background = element_rect(fill = "white"),      # Set background to white
      panel.grid = element_blank()                         # Remove grid lines
    ) +
    scale_x_continuous(breaks = x_labels$x_pos, labels = x_labels$Database)
  # Save the Manhattan plot
  #ggsave(filename = paste0(comparison_name, "_manhattan_plot_v12.png"), plot = manhattan_plot, width = 10, height = 8)
  ggsave(filename = paste0(output_dir, "/", comparison_name, "_manhattan_plot_v12.pdf"), plot = manhattan_plot, width = 10, height = 8, dpi = 300, bg = "white")
  print(paste("Manhattan plot with properly aligned labels saved for", comparison_name))
}


output_dir <- "/Users/kanatyalcin/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Get the filtered genes (gene symbols) for the current comparison
filtered_genes <- pairwise_results[["Cluster_11_vs_All_padj_log2fc_2"]]$gene_symbol

# Run enrichment analysis on Biological Processes, KEGG, and Transcription Factors
enrich_results_bp_kegg <- enrichr(filtered_genes, selected_dbs_bp_kegg)
enrich_results_tf <- enrichr(filtered_genes, selected_dbs_tf)


# Merge the pathway and transcription factor results
enrich_results <- c(enrich_results_bp_kegg, enrich_results_tf)

# Create and save the Manhattan plot for this comparison
create_manhattan_plot_v12(enrich_results, "Cluster_11_vs_All_padj_log2fc_2", output_dir)
