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

#########################################################
### Step 1 : Data Import and Pre-processing
#########################################################
set.seed(123)

counts <- read.csv("/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/FOS_ManuscriptSubmission_2024/Data/merged_counts_clean.csv",
                   header = TRUE, 
                   row.names = 1, 
                   sep = ",",
                   check.names = FALSE)


sampleinfo <- read.delim("/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/FOS_ManuscriptSubmission_2024/Data/s_info.csv",
                         header = TRUE, sep = ",",check.names = FALSE)

sampleinfo$Sample
colnames(counts)


fhcc_nf2_cases <- c(
  "SRR27388428", "SRR27388430", "SRR27388436", "SRR27388493", "SRR27388496",
  "SRR27388537", "SRR27388538", "SRR27388539", "SRR27388543", "SRR27388548",
  "SRR27388588", "SRR27388589", "SRR27388590", "SRR27388591", "SRR27388604",
  "SRR27388610", "SRR27388613", "SRR27388615", "SRR27388641", "SRR27388670",
  "SRR27388679", "SRR27388685", "SRR27388694", "SRR27388700", "SRR27388701"
)

mark_nf2_cases <- c(
  "SRR3995999",
  "SRR3995996",
  "SRR3995988",
  "SRR3995990",
  "SRR3995991",
  "SRR3995992",
  "SRR3995993",
  "ERR12048882"
)

mark_hedgehog_cases <- c(
  "ERR12048892",
  "ERR12048883",
  "ERR12048894",
  "ERR12048887",
  "ERR12048886",
  "ERR12048898",
  "ERR12048925",
  "ERR12048899",
  "ERR12048900",
  "ERR12048885"
)


mark_TRAF7_KLF4_cases <- c(
  "SRR3996000",
  "SRR3996004",
  "SRR3995995",
  "SRR3995989"
)

mark_TRAF7_PI3K_cases <- c(
  "SRR3996003",
  "SRR3995987",
  "SRR3995994"
)

mark_POLR2A_cases <- c(
  "SRR3995998",
  "SRR3996001",
  "SRR3996002",
  "SRR3996005"
)

mark_3PLOSS_cases <- c(
  "ERR12048890",
  "ERR12048893",
  "ERR12048901"
)


# Function to determine NF2_STATUS
determine_nf2_status <- function(srr, fos_status) {
  # Check if FOS_Fusion_Status is "Yes"
  if (fos_status == "Yes") {
    return("FOS")
  }
  
  # Check if SRR is in the list of specific NF2 IDs
  if (srr %in% fhcc_nf2_cases) {
    return("NF2")
  }
  if (srr %in% mark_nf2_cases) {
    return("NF2")
  }
  if (srr %in% mark_hedgehog_cases) {
    return("SMO/HH")
  }
  if (srr %in% mark_TRAF7_KLF4_cases) {
    return("TRAF7/KLF4")
  }
  if (srr %in% mark_TRAF7_PI3K_cases) {
    return("TRAF7/PI3K")
  }
  if (srr %in% mark_POLR2A_cases) {
    return("POLR2A")
  }
  
  
  return("NA")
}

# Applying the function to determine NF2_STATUS for each Sample in sampleinfo
sampleinfo$NF2_STATUS <- mapply(determine_nf2_status, sampleinfo$Sample, sampleinfo$FOS_Fusion_Status)


# Display the updated sampleinfo DataFrame
print(sampleinfo)



sampleinfo$Sample <- as.factor(sampleinfo$Sample)
sampleinfo$FOS_Fusion_Status <- as.factor(sampleinfo$FOS_Fusion_Status)
sampleinfo$Batch <- as.factor(sampleinfo$Batch)
sampleinfo$NF2_STATUS <- as.factor(sampleinfo$NF2_STATUS)

class(sampleinfo)
class(sampleinfo$Sample)
class(sampleinfo$FOS_Fusion_Status)
class(sampleinfo$Batch)

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
umap_result_3d <- umap(pca_data, n_neighbors = 7, min_dist = 0.1, n_components = 3, config = custom.config)
# Extract the layout component, which contains the 3D coordinates
umap_coordinates_3d <- umap_result_3d$layout
# Remove rows with NA or Inf values
umap_coordinates_3d <- na.omit(umap_coordinates_3d)  # Remove NA rows
umap_coordinates_3d <- umap_coordinates_3d[apply(umap_coordinates_3d, 1, function(row) all(is.finite(row))), ]  # Remove rows with Inf

batch_labels <- sampleinfo$Batch
colors <- rainbow(length(unique(batch_labels)))


# Create a data frame for plotting
umap_df <- data.frame(
  UMAP1 = umap_coordinates_3d[, 1],
  UMAP2 = umap_coordinates_3d[, 2],
  UMAP3 = umap_coordinates_3d[, 3],
  FOS_Fusion_Status = factor(batch_labels),
  NF2_Status =sampleinfo$NF2_STATUS,
  Batch = sampleinfo$Batch,
  ID= sampleinfo$Sample
)

dbscan_result_3d <- dbscan(umap_coordinates_3d, eps = 0.65, minPts = 8) ## min pts = 8 

# Access the cluster labels
cluster_labels_3d <- dbscan_result_3d$cluster

umap_df_3d <- as.data.frame(umap_coordinates_3d)
colnames(umap_df_3d) <- c("UMAP_1", "UMAP_2", "UMAP_3")

# Add cluster labels to the data frame
umap_df_3d$Cluster <- as.factor(cluster_labels_3d)

umap_df_3d$FOS_Fusion_Status <- as.factor(sampleinfo$FOS_Fusion_Status)

umap_df_3d$Sample_Name <- sampleinfo$Sample  # Add sample names

umap_df_3d$Batch <- sampleinfo$Batch

umap_df_3d$NF2_STATUS <- sampleinfo$NF2_STATUS


batch_labels <- umap_df_3d$Cluster

# Get unique clusters and ensure consistent ordering
clusters <- sort(unique(batch_labels))

# Number of clusters
n_clusters <- length(clusters)


colors <- brewer.pal(n_clusters, "Set2")


# Create a mapping from clusters to colors
umap_df_3d_noise_removed <- subset(umap_df_3d, Cluster != 0)


batch_labels <- umap_df_3d_noise_removed$Cluster

clusters <- sort(unique(batch_labels))

# Number of clusters
n_clusters <- length(clusters)
colors <- c("red", "gray", "blue", "orange", "cyan", "purple", "green")
# Create a 3D scatter plot with plotly, with colors for clusters, shapes for FOS_Fusion_Status, and sample names on hover
three_dplot_dbscan <- plot_ly(umap_df_3d_noise_removed, 
                              x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
                              color = ~NF2_STATUS, colors = colors, 
                              #symbol = ~Cluster, symbols = symbol_set, # Adjust symbols as needed
                              type = "scatter3d", mode = "markers",
                              marker = list(size = 6),
                              text = ~paste("Sample: ", Sample_Name, "<br>Mutation Status: ", NF2_STATUS , "<br>Cluster: ", Cluster, "<br>Sequencing Center: ", Batch),  # Concatenate sample info for hover
                              hoverinfo = "text") %>%  # Ensure hover text is enabled
  layout(title = "3D UMAP Plot with DBSCAN Clusters, FOS Fusion Status, and Sample Information",
         scene = list(xaxis = list(title = 'UMAP 1'),
                      yaxis = list(title = 'UMAP 2'),
                      zaxis = list(title = 'UMAP 3')),
         showlegend = TRUE)  # Disable legend


# Define the camera settings
camera <- list(
  eye = list(x = -2.0268997808002043, y = 1.3233692539082116, z = 0.8787230490557443)  # Adjust these values to change the angle
)

# Update the layout with the camera settings and customize fonts
three_dplot_dbscan <- three_dplot_dbscan %>%
  layout(
    showlegend = FALSE,
    title = NULL,
    margin = list(l = 80, r = 80, b = 80, t = 80),
    scene = list(
      camera = camera,
      aspectmode = 'cube',
      domain = list(
        x = c(0, 1),
        y = c(0, 1)
      ),
      xaxis = list(
        title = list(
          text = 'UMAP2',
          font = list(
            family = 'Arial Bold, sans-serif',
            size = 40,
            color = 'black'
          ),
          standoff = 30
        ),
        showgrid = TRUE,
        tickfont = list(
          family = 'Arial Bold, sans-serif',
          size = 20,
          color = 'black'
        ),
        showline = TRUE,
        zeroline = FALSE,
        linewidth = 8,
        linecolor = 'black',
        automargin = TRUE
      ),
      yaxis = list(
        title = list(
          text = 'UMAP3',
          font = list(
            family = 'Arial Bold, sans-serif',
            size = 40,
            color = 'black'
          ),
          standoff = 30
        ),
        showgrid = TRUE,
        tickfont = list(
          family = 'Arial Bold, sans-serif',
          size = 20,
          color = 'black'
        ),
        showline = TRUE,
        zeroline = FALSE,
        linewidth = 8,
        linecolor = 'black',
        automargin = TRUE
      ),
      zaxis = list(
        title = list(
          text = 'UMAP1',
          font = list(
            family = 'Arial Bold, sans-serif',
            size = 40,
            color = 'black'
          ),
          standoff = 30
        ),
        showgrid = TRUE,
        tickfont = list(
          family = 'Arial Bold, sans-serif',
          size = 20,
          color = 'black'
        ),
        showline = TRUE,
        zeroline = FALSE,
        linewidth = 8,
        linecolor = 'black',
        automargin = TRUE
      ),
      # Remove the background and grid
      bgcolor = 'rgba(0,0,0,0)'
    )
  )


save_image(three_dplot_dbscan, "/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/FOS_ManuscriptSubmission_2024/Figures/3D_UMAP_plot_mutation_status.pdf", width = 1920, height = 1920)
htmlwidgets::saveWidget(three_dplot_dbscan, "/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/FOS_ManuscriptSubmission_2024/Figures/3D_UMAP_plot_mutation_status.html")

batch_labels <- umap_df_3d_noise_removed$Cluster

clusters <- sort(unique(batch_labels))

# Number of clusters
n_clusters <- length(clusters)
colors <- brewer.pal(n_clusters, "Set2")

# Create a 3D scatter plot with plotly, with colors for clusters, shapes for FOS_Fusion_Status, and sample names on hover
three_dplot_dbscan <- plot_ly(umap_df_3d_noise_removed, 
                              x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
                              color = ~Cluster, colors = colors, 
                              #symbol = ~Cluster, symbols = symbol_set, # Adjust symbols as needed
                              type = "scatter3d", mode = "markers",
                              marker = list(size = 6),
                              text = ~paste("Sample: ", Sample_Name, "<br>Mutation Status: ", NF2_STATUS , "<br>Cluster: ", Cluster, "<br>Sequencing Center: ", Batch),  # Concatenate sample info for hover
                              hoverinfo = "text") %>%  # Ensure hover text is enabled
  layout(title = "3D UMAP Plot with DBSCAN Clusters, FOS Fusion Status, and Sample Information",
         scene = list(xaxis = list(title = 'UMAP 1'),
                      yaxis = list(title = 'UMAP 2'),
                      zaxis = list(title = 'UMAP 3')),
         showlegend = TRUE)  # Disable legend


# Define the camera settings
camera <- list(
  eye = list(x = -2.0268997808002043, y = 1.3233692539082116, z = 0.8787230490557443)  # Adjust these values to change the angle
)

# Update the layout with the camera settings and customize fonts
three_dplot_dbscan <- three_dplot_dbscan %>%
  layout(
    showlegend = FALSE,
    title = NULL,
    margin = list(l = 80, r = 80, b = 80, t = 80),
    scene = list(
      camera = camera,
      aspectmode = 'cube',
      domain = list(
        x = c(0, 1),
        y = c(0, 1)
      ),
      xaxis = list(
        title = list(
          text = 'UMAP2',
          font = list(
            family = 'Arial Bold, sans-serif',
            size = 40,
            color = 'black'
          ),
          standoff = 30
        ),
        showgrid = TRUE,
        tickfont = list(
          family = 'Arial Bold, sans-serif',
          size = 20,
          color = 'black'
        ),
        showline = TRUE,
        zeroline = FALSE,
        linewidth = 8,
        linecolor = 'black',
        automargin = TRUE
      ),
      yaxis = list(
        title = list(
          text = 'UMAP3',
          font = list(
            family = 'Arial Bold, sans-serif',
            size = 40,
            color = 'black'
          ),
          standoff = 30
        ),
        showgrid = TRUE,
        tickfont = list(
          family = 'Arial Bold, sans-serif',
          size = 20,
          color = 'black'
        ),
        showline = TRUE,
        zeroline = FALSE,
        linewidth = 8,
        linecolor = 'black',
        automargin = TRUE
      ),
      zaxis = list(
        title = list(
          text = 'UMAP1',
          font = list(
            family = 'Arial Bold, sans-serif',
            size = 40,
            color = 'black'
          ),
          standoff = 30
        ),
        showgrid = TRUE,
        tickfont = list(
          family = 'Arial Bold, sans-serif',
          size = 20,
          color = 'black'
        ),
        showline = TRUE,
        zeroline = FALSE,
        linewidth = 8,
        linecolor = 'black',
        automargin = TRUE
      ),
      # Remove the background and grid
      bgcolor = 'rgba(0,0,0,0)'
    )
  )

three_dplot_dbscan

save_image(three_dplot_dbscan, "/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/FOS_ManuscriptSubmission_2024/Figures/3D_UMAP_plot_cluster.pdf", width = 1920, height = 1920)
htmlwidgets::saveWidget(three_dplot_dbscan, "/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/FOS_ManuscriptSubmission_2024/Figures/3D_UMAP_plot_cluster.html")

batch_labels <- umap_df_3d_noise_removed$Batch

clusters <- sort(unique(batch_labels))

# Number of clusters
n_clusters <- length(clusters)
colors <- brewer.pal(n_clusters, "Set2")

# Create a 3D scatter plot with plotly, with colors for clusters, shapes for FOS_Fusion_Status, and sample names on hover
three_dplot_dbscan <- plot_ly(umap_df_3d_noise_removed, 
                              x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
                              color = ~Batch, colors = colors, 
                              #symbol = ~Cluster, symbols = symbol_set, # Adjust symbols as needed
                              type = "scatter3d", mode = "markers",
                              marker = list(size = 6),
                              text = ~paste("Sample: ", Sample_Name, "<br>Mutation Status: ", NF2_STATUS , "<br>Cluster: ", Cluster, "<br>Sequencing Center: ", Batch),  # Concatenate sample info for hover
                              hoverinfo = "text") %>%  # Ensure hover text is enabled
  layout(title = "3D UMAP Plot with DBSCAN Clusters, FOS Fusion Status, and Sample Information",
         scene = list(xaxis = list(title = 'UMAP 1'),
                      yaxis = list(title = 'UMAP 2'),
                      zaxis = list(title = 'UMAP 3')),
         showlegend = TRUE)  # Disable legend


# Define the camera settings
camera <- list(
  eye = list(x = -2.0268997808002043, y = 1.3233692539082116, z = 0.8787230490557443)  # Adjust these values to change the angle
)

# Update the layout with the camera settings and customize fonts
three_dplot_dbscan <- three_dplot_dbscan %>%
  layout(
    showlegend = FALSE,
    title = NULL,
    margin = list(l = 80, r = 80, b = 80, t = 80),
    scene = list(
      camera = camera,
      aspectmode = 'cube',
      domain = list(
        x = c(0, 1),
        y = c(0, 1)
      ),
      xaxis = list(
        title = list(
          text = 'UMAP2',
          font = list(
            family = 'Arial Bold, sans-serif',
            size = 40,
            color = 'black'
          ),
          standoff = 30
        ),
        showgrid = TRUE,
        tickfont = list(
          family = 'Arial Bold, sans-serif',
          size = 20,
          color = 'black'
        ),
        showline = TRUE,
        zeroline = FALSE,
        linewidth = 8,
        linecolor = 'black',
        automargin = TRUE
      ),
      yaxis = list(
        title = list(
          text = 'UMAP3',
          font = list(
            family = 'Arial Bold, sans-serif',
            size = 40,
            color = 'black'
          ),
          standoff = 30
        ),
        showgrid = TRUE,
        tickfont = list(
          family = 'Arial Bold, sans-serif',
          size = 20,
          color = 'black'
        ),
        showline = TRUE,
        zeroline = FALSE,
        linewidth = 8,
        linecolor = 'black',
        automargin = TRUE
      ),
      zaxis = list(
        title = list(
          text = 'UMAP1',
          font = list(
            family = 'Arial Bold, sans-serif',
            size = 40,
            color = 'black'
          ),
          standoff = 30
        ),
        showgrid = TRUE,
        tickfont = list(
          family = 'Arial Bold, sans-serif',
          size = 20,
          color = 'black'
        ),
        showline = TRUE,
        zeroline = FALSE,
        linewidth = 8,
        linecolor = 'black',
        automargin = TRUE
      ),
      # Remove the background and grid
      bgcolor = 'rgba(0,0,0,0)'
    )
  )

save_image(three_dplot_dbscan, "/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/FOS_ManuscriptSubmission_2024/Figures/3D_UMAP_plot_batch.pdf", width = 1920, height = 1920)
htmlwidgets::saveWidget(three_dplot_dbscan, "/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/FOS_ManuscriptSubmission_2024/Figures/3D_UMAP_plot_batch.html")

#########################################################################
### Step 4 : DIFFERENTIAL GENE EXPRESSION AND ENRICHMENT ANALYSIS
#########################################################################
metadata <- data.frame(
  row.names = colnames(combat_counts),  # Ensure that row names match sample IDs in count matrix
  cluster = umap_df_3d$Cluster  # Use the cluster info from your UMAP data
)

metadata$group <- ifelse(metadata$cluster == "2", "Cluster_2", "Rest")

metadata$group <- factor(metadata$group, levels = c("Rest", "Cluster_2"))

metadata$group <- relevel(metadata$group, ref = "Rest")

# Create the DESeq2 dataset object
dds_dbscan <- DESeqDataSetFromMatrix(countData = combat_counts, 
                                     colData = metadata, 
                                     design = ~ group)

# Run the DESeq2 differential expression analysis
deseq_dbscan <- DESeq(dds_dbscan)

#res <- results(deseq_dbscan, contrast = c("group", "Cluster_2", "Rest"))
res <- results(deseq_dbscan, contrast = c("group", "Rest", "Cluster_2"))

# Summarize the results
summary(res)

rownames(res) <- gsub("\\.[0-9]*","",rownames(res))

head(res)

res$gene_symbol <- mapIds(EnsDb.Hsapiens.v86,
                          keys=row.names(res),
                          column="SYMBOL",
                          keytype="GENEID",
                          multiVals="first")
n_occur <- data.frame(table(res$gene_symbol))

res$biotype <- mapIds(EnsDb.Hsapiens.v86,
                      keys=row.names(res),
                      column="GENEBIOTYPE",
                      keytype="GENEID",
                      multiVals="first")

res_proteincoding <- dplyr::filter(as.data.frame(res), biotype == "protein_coding")
dim(res_proteincoding)

res_proteincoding <- na.omit(res_proteincoding)

# Apply a stricter p-value threshold and minimum log2 fold change
res_combat_proteincoding_padj_005_log2fold_5_proteincoding <- res_proteincoding %>%
  dplyr::filter(padj < 0.05, abs(log2FoldChange) > 5)

#write.csv(res_combat_proteincoding_padj_005_log2fold_5_proteincoding, file = "/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/HTseq_500/Results/res_combat_proteincoding_padj_005_log2fold_5_proteincoding_FOS_vs_All.csv") 

filtered_genes <- res_combat_proteincoding_padj_005_log2fold_5_proteincoding %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 5) %>%
  pull(gene_symbol)  # Extract the gene column


# List databases available in enrichR
dbs <- listEnrichrDbs()

# Databases
selected_dbs <- c("GO_Biological_Process_2021", "KEGG_2021_Human")
selected_dbs

# Run enrichment analysis on your filtered genes
enrich_results <- enrichr(filtered_genes, selected_dbs)
enrich_results

# Extract the GO results
go_results <- enrich_results[["GO_Biological_Process_2021"]]

# Extract the KEGG results
kegg_results <- enrich_results[["KEGG_2021_Human"]]

# Select top 10 results for GO and KEGG
top_go <- go_results %>% top_n(-10, P.value)
top_kegg <- kegg_results %>% top_n(-10, P.value)

# Combine GO and KEGG results
top_enrich <- bind_rows(
  top_go %>% mutate(Database = "GO"),
  top_kegg %>% mutate(Database = "KEGG")
)

# Create a bar plot
ggplot(top_enrich, aes(x = reorder(Term, -log10(P.value)), y = -log10(P.value), fill = Database)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Top Enriched GO Terms and KEGG Pathways",
    x = "Terms",
    y = "-log10(p-value)"
  ) +
  theme_minimal()

enrichR_top <- ggplot(top_enrich, aes(x = reorder(Term, -log10(P.value)), y = -log10(P.value), fill = Database)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Top Enriched GO Terms and KEGG Pathways",
    x = "Terms",
    y = "-log10(p-value)"
  ) +
  theme_minimal(base_size = 14) +  # IText size - ask zeyep if its okay !!!
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )
#ggsave("/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/HTseq_500/Figures/enrichR_top_FOS_vs_All.png", plot = enrichR_top, width = 12, height = 8, dpi = 300, units = "in")

selected_dbs <- c("TRRUST_Transcription_Factors_2019", "ChEA_2016", "TRANSFAC_and_JASPAR_PWMs")

# Run enrichment analysis again but focusing on transcription factor databases
tf_enrich_results <- enrichr(filtered_genes, selected_dbs)

# Extract the TRRUST and ChEA results
trrust_results <- tf_enrich_results[["TRRUST_Transcription_Factors_2019"]]
chea_results <- tf_enrich_results[["ChEA_2016"]]
trasnfac_results <- tf_enrich_results[["TRANSFAC_and_JASPAR_PWMs"]]

# Filter results where "human" is found in the 'Term' column using grep or str_detect (ignoring case)
trrust_human <- trrust_results %>% dplyr::filter(str_detect(Term, regex("human", ignore_case = TRUE)))
chea_human <- chea_results %>% dplyr::filter(str_detect(Term, regex("human", ignore_case = TRUE)))
transfac_human <- trasnfac_results %>% dplyr::filter(str_detect(Term, regex("human", ignore_case = TRUE)))

# Select top 10 results for TRRUST and ChEA
top_trrust <- trrust_human %>% top_n(-10, P.value)
top_chea <- chea_human %>% top_n(-10, P.value)
top_transfac <- transfac_human %>% top_n(-10, P.value)

# Combine TRRUST and ChEA results
top_tf_enrich <- bind_rows(
  top_trrust %>% mutate(Database = "TRRUST"),
  top_chea %>% mutate(Database = "ChEA"),
  top_transfac %>% mutate(Database = "TRANSFAC_and_JASPAR_PWMs"),
)

# Bar plot for transcription factor enrichment
enrichR_top_TF <- ggplot(top_tf_enrich, aes(x = reorder(Term, -log10(P.value)), y = -log10(P.value), fill = Database)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Top Enriched Human Transcription Factors (TRRUST and ChEA and TRANSFAC)",
    x = "Transcription Factors",
    y = "-log10(p-value)"
  ) +
  theme_minimal()

#### Pairwise DESeq ### 

other_clusters <- levels(umap_df_3d$Cluster)[levels(umap_df_3d$Cluster) != "2"]

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

# Function to convert gene IDs to gene symbols and filter by protein-coding
convert_but_dont_filter_protein_coding <- function(res) {
  # Remove version numbers from gene IDs
  rownames(res) <- gsub("\\.[0-9]*", "", rownames(res))
  
  # Map gene symbols and biotypes
  res$gene_symbol <- mapIds(EnsDb.Hsapiens.v86,
                            keys = row.names(res),
                            column = "SYMBOL",
                            keytype = "GENEID",
                            multiVals = "first")
  
  #res$biotype <- mapIds(EnsDb.Hsapiens.v86,
  #                      keys = row.names(res),
  #                      column = "GENEBIOTYPE",
  #                      keytype = "GENEID",
  #                      multiVals = "first")
  
  # Filter for protein-coding genes and remove NAs
  #res_proteincoding <- dplyr::filter(as.data.frame(res), biotype == "protein_coding")
  res <- na.omit(res)
  
  return(res)
}

# Create a loop to compare Cluster 2 with each other cluster one by one
for (cluster in other_clusters) {
  
  # Prepare metadata (sample info), labeling "Cluster_2" vs the current comparison cluster
  metadata <- data.frame(
    row.names = colnames(combat_counts),  # Ensure that row names match sample IDs in the count matrix
    cluster = umap_df_3d$Cluster  # Use the cluster info from your UMAP data
  )
  
  # Create a new column marking Cluster 2 as "Cluster_2" and the current comparison cluster as "Other"
  metadata$group <- ifelse(metadata$cluster == "2", "Cluster_2", 
                           ifelse(metadata$cluster == cluster, paste0("Cluster_", cluster), NA))
  
  # Subset data to only include samples from Cluster 2 and the comparison cluster
  dds_subset <- DESeqDataSetFromMatrix(countData = combat_counts[, !is.na(metadata$group)], 
                                       colData = metadata[!is.na(metadata$group), ], 
                                       design = ~ group)
  
  # Run DESeq2
  dds_subset <- DESeq(dds_subset)
  
  # Extract results comparing Cluster 2 vs the comparison cluster
  res_c2_vs_others <- results(dds_subset, contrast = c("group", "Cluster_2", paste0("Cluster_", cluster)))
  
  # Convert gene IDs to gene symbols and filter for protein-coding genes
  res_proteincoding <- convert_and_filter_protein_coding(res_c2_vs_others)
  
  res_raw <- convert_but_dont_filter_protein_coding(res_c2_vs_others)
  
  # Save raw protein-coding filtered results (no padj or log2fold filtering)
  pairwise_results[[paste0("Cluster_2_vs_Cluster_", cluster, "_raw_protein_coding")]] <- res_proteincoding
  
  pairwise_results[[paste0("Cluster_2_vs_Cluster_", cluster, "_raw")]] <- res_raw
  
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
  pairwise_results[[paste0("Cluster_2_vs_Cluster_", cluster, "_padj_abs_log2fc_5")]] <- res_filtered_abs_5
  pairwise_results[[paste0("Cluster_2_vs_Cluster_", cluster, "_padj_log2fc_5")]] <- res_filtered_log2fc_5
  pairwise_results[[paste0("Cluster_2_vs_Cluster_", cluster, "_padj_log2fc_2")]] <- res_filtered_log2fc_2
  pairwise_results[[paste0("Cluster_2_vs_Cluster_", cluster, "_padj_abs_log2fc_2")]] <- res_filtered_abs_2
  
  # Optionally print a message to track progress
  print(paste("Stored results for Cluster 2 vs Cluster", cluster))
}

# Additionally, compare Cluster 2 vs all other clusters combined
metadata_all <- data.frame(
  row.names = colnames(combat_counts), 
  cluster = umap_df_3d$Cluster
)

# Create a new column marking Cluster 2 as "Cluster_2" and all others as "Other"
metadata_all$group <- ifelse(metadata_all$cluster == "2", "Cluster_2", "Other")

# Subset data to only include samples from Cluster 2 and all other clusters
dds_all <- DESeqDataSetFromMatrix(countData = combat_counts[, !is.na(metadata_all$group)], 
                                  colData = metadata_all[!is.na(metadata_all$group), ], 
                                  design = ~ group)

# Run DESeq2
dds_all <- DESeq(dds_all)

# Extract results comparing Cluster 2 vs all others
res_c2_vs_all <- results(dds_all, contrast = c("group", "Cluster_2", "Other"))

# Convert gene IDs to gene symbols and filter for protein-coding genes
res_proteincoding_all <- convert_and_filter_protein_coding(res_c2_vs_all)

res_raw_all <- convert_but_dont_filter_protein_coding(res_c2_vs_all)

# Save raw protein-coding filtered results (no padj or log2fold filtering) for Cluster 2 vs All
pairwise_results[["Cluster_2_vs_All_raw_protein_coding"]] <- res_proteincoding_all

pairwise_results[["Cluster_2_vs_All_raw"]] <- res_raw_all

# Apply the same filtering as above
res_filtered_all_abs_5 <- res_proteincoding_all[which(res_proteincoding_all$padj < 0.05 & abs(res_proteincoding_all$log2FoldChange) > 5), ]
res_filtered_all_log2fc_5 <- res_proteincoding_all[which(res_proteincoding_all$padj < 0.05 & res_proteincoding_all$log2FoldChange > 5), ]
res_filtered_all_log2fc_2 <- res_proteincoding_all[which(res_proteincoding_all$padj < 0.05 & res_proteincoding_all$log2FoldChange > 2), ]
res_filtered_all_abs_2 <- res_proteincoding_all[which(res_proteincoding_all$padj < 0.05 & abs(res_proteincoding_all$log2FoldChange) > 2), ]

# Save filtered results to the list with unique names for Cluster 2 vs All
pairwise_results[["Cluster_2_vs_All_padj_abs_log2fc_5"]] <- res_filtered_all_abs_5
pairwise_results[["Cluster_2_vs_All_padj_log2fc_5"]] <- res_filtered_all_log2fc_5
pairwise_results[["Cluster_2_vs_All_padj_log2fc_2"]] <- res_filtered_all_log2fc_2
pairwise_results[["Cluster_2_vs_All_padj_abs_log2fc_2"]] <- res_filtered_all_abs_2

print("Stored results for Cluster 2 vs All")

## <---- BURDAYIZ
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
data_df <- as.data.frame(pairwise_results[["Cluster_2_vs_All_raw_protein_coding"]])

# Identify the top 10 most upregulated significant genes
top_upregulated <- data_df %>%
  filter(padj < 0.05 & log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)

# Identify the top 10 most downregulated significant genes
top_downregulated <- data_df %>%
  filter(padj < 0.05 & log2FoldChange < 0) %>%
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
ggsave("/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/FOS_ManuscriptSubmission_2024/Figures/Volcano_Plot_padj005_log2FC2_proteincoding.png", volcano_plot_top20, width = 29, height = 16, dpi = 300)

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


output_dir <- "/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/FOS_ManuscriptSubmission_2024/Figures"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


gobp_z <- read.csv("/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/FOS_ManuscriptSubmission_2024/Data/Cluster2vsAll_UpregulatedGenes_logFC2_GO_BP_enrich.txt",
                   header = TRUE, 
                   row.names = 1, 
                   sep = "\t",
                   check.names = FALSE)




# Get the filtered genes (gene symbols) for the current comparison
filtered_genes <- pairwise_results[["Cluster_2_vs_All_padj_log2fc_2"]]$gene_symbol

# If no genes, skip this iteration
#if (length(filtered_genes) == 0) {pairwise_results[["Cluster_2_vs_All_padj_log2fc_2"]]
#  next
#}

# Run enrichment analysis on Biological Processes, KEGG, and Transcription Factors
enrich_results_bp_kegg <- enrichr(filtered_genes, selected_dbs_bp_kegg)
enrich_results_tf <- enrichr(filtered_genes, selected_dbs_tf)

enrich_results_bp_kegg$GO_Biological_Process_2021 <- gobp_z

# Merge the pathway and transcription factor results
enrich_results <- c(enrich_results_bp_kegg, enrich_results_tf)

# Create and save the Manhattan plot for this comparison
create_manhattan_plot_v12(enrich_results, "Cluster_2_vs_All_padj_log2fc_2", output_dir)

save.image("/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/FOS_ManuscriptSubmission_2024/Code/FOS_742samples_UMAP_Volcano_Manhattan.RData")

#########################################################################################
### Step 6 : GENE EXPRESSION COMPARISON BETWEEN CLUSTERS FOR SPECIFIC GENES OF INTEREST
#########################################################################################
                                                 
combat_counts_vst_geneSymbol <- as.data.frame(counts(deseq_dds,normalized=TRUE))

rownames(combat_counts_vst_geneSymbol) <- gsub("\\.[0-9]*","",rownames(combat_counts_vst_geneSymbol))
head(combat_counts_vst_geneSymbol)

combat_counts_vst_geneSymbol$gene_symbol <- mapIds(EnsDb.Hsapiens.v86,
                          keys=row.names(combat_counts_vst_geneSymbol),
                          column="SYMBOL",
                          keytype="GENEID",
                          multiVals="first")
combat_counts_vst_geneSymbol$biotype <- mapIds(EnsDb.Hsapiens.v86,
                             keys=row.names(combat_counts_vst_geneSymbol),
                             column="GENEBIOTYPE",
                             keytype="GENEID",
                             multiVals="first")

combat_counts_vst_geneSymbol <- na.omit(combat_counts_vst_geneSymbol)
# Define the genes of interest

genes_of_interest <- c("CRLF2", "ARHGAP36", "PPARG", "DLK1", "CEBPA", "NOTCH3", "JAK3", "BRD4", "FCER1A", "KIT", "IL6", "TNFRSF9", "CFB", "IHH", "DHH", "SHH", "TPSAB1", "TPSB2")

# Identify the sample columns (excluding 'gene_symbol' and 'biotype')
sample_columns <- setdiff(colnames(combat_counts_vst_geneSymbol), c("gene_symbol", "biotype"))

# Extract expression data for the genes of interest
expr_df <- combat_counts_vst_geneSymbol[combat_counts_vst_geneSymbol$gene_symbol %in% genes_of_interest, sample_columns]

# Set the row names to gene symbols for clarity
rownames(expr_df) <- combat_counts_vst_geneSymbol$gene_symbol[combat_counts_vst_geneSymbol$gene_symbol %in% genes_of_interest]

# Transpose the expression dataframe to align samples with `umap_df_3d` row indices
expr_df_t <- t(expr_df)

# Convert the transposed data to a data frame
expr_df_t <- as.data.frame(expr_df_t)

# Ensure that sample names match between the two dataframes
common_samples <- intersect(rownames(umap_df_3d), rownames(expr_df_t))

# Subset both dataframes to include only the common samples
umap_df_3d_sub <- umap_df_3d[common_samples, ]
expr_df_t_sub <- expr_df_t[common_samples, ]

# Combine the dataframes by adding the expression data as new columns
umap_df_3d_with_genes <- cbind(umap_df_3d_sub, expr_df_t_sub)

color_palette <- c("#FC8D62", "#B79FEC", "#E78AC3", "#C49A6C", "#B3DE69", "#FFD92F", "#E5C494", "#B3B3B3")

umap_df_3d_noise_removed_with_genes <- subset(umap_df_3d_with_genes, Cluster != 0)
p <- ggplot(umap_df_3d_noise_removed_with_genes, aes(x = factor(Cluster), y = TPSB2, fill = factor(Cluster))) +
  geom_boxplot() +
  scale_fill_manual(values = color_palette) +
  labs(x = "Cluster", y = "Proportion", title = "TPSB2", fill = "Cluster") +  # Remove dots from title
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(face = "bold", size = 16)) +
  geom_signif(comparisons = list(c(2,1), c(2,3), c(2,4), c(2,5), c(2,6), c(2,7), c(2,8)),
              map_signif_level = FALSE, 
              test = "wilcox.test", 
              step_increase = 0.05)

ggsave("plot_TPSB2.pdf", plot = p, width = 10, height = 9, dpi = 300, bg = "white")

### Log10 ###

combat_counts_deseq_normalized_geneSymbol <- as.data.frame(counts(deseq_dbscan, normalized = TRUE))
rownames(combat_counts_deseq_normalized_geneSymbol) <- gsub("\\.[0-9]*", "", rownames(combat_counts_deseq_normalized_geneSymbol))
head(combat_counts_deseq_normalized_geneSymbol)

combat_counts_deseq_normalized_geneSymbol$gene_symbol <- mapIds(
  EnsDb.Hsapiens.v86,
  keys = row.names(combat_counts_deseq_normalized_geneSymbol),
  column = "SYMBOL",
  keytype = "GENEID",
  multiVals = "first"
)

combat_counts_deseq_normalized_geneSymbol$biotype <- mapIds(
  EnsDb.Hsapiens.v86,
  keys = row.names(combat_counts_deseq_normalized_geneSymbol),
  column = "GENEBIOTYPE",
  keytype = "GENEID",
  multiVals = "first"
)

combat_counts_deseq_normalized_geneSymbol <- na.omit(combat_counts_deseq_normalized_geneSymbol)

# Define the genes of interest
genes_of_interest <- c("CRLF2", "ARHGAP36", "PPARG", "DLK1", "CEBPA", "NOTCH3", "JAK3", "BRD4",
                       "FCER1A", "KIT", "IL6", "TNFRSF9", "CFB", "IHH", "DHH", "SHH", "TPSAB1", "TPSB2")

# Identify the sample columns (excluding 'gene_symbol' and 'biotype')
sample_columns <- setdiff(colnames(combat_counts_deseq_normalized_geneSymbol), c("gene_symbol", "biotype"))

# Extract expression data for the genes of interest
expr_df <- combat_counts_deseq_normalized_geneSymbol[
  combat_counts_deseq_normalized_geneSymbol$gene_symbol %in% genes_of_interest,
  sample_columns
]

# Set the row names to gene symbols for clarity
rownames(expr_df) <- combat_counts_deseq_normalized_geneSymbol$gene_symbol[
  combat_counts_deseq_normalized_geneSymbol$gene_symbol %in% genes_of_interest
]

# Transpose the expression dataframe to align samples with `umap_df_3d` row indices
expr_df_t <- t(expr_df)
expr_df_t <- as.data.frame(expr_df_t)

# Ensure that sample names match between the two dataframes
common_samples <- intersect(rownames(umap_df_3d), rownames(expr_df_t))

# Subset both dataframes to include only the common samples
umap_df_3d_sub <- umap_df_3d[common_samples, ]
expr_df_t_sub <- expr_df_t[common_samples, ]

# Combine the dataframes by adding the expression data as new columns
umap_df_3d_with_genes <- cbind(umap_df_3d_sub, expr_df_t_sub)

color_palette <- c("#FC8D62", "#B79FEC", "#E78AC3", "#C49A6C", "#B3DE69", "#FFD92F", "#E5C494", "#B3B3B3")

# Remove noise cluster
umap_df_3d_noise_removed_with_genes <- subset(umap_df_3d_with_genes, Cluster != 0)

# Add a pseudocount to TPSB2 before taking log10
umap_df_3d_noise_removed_with_genes$TPSB2_log <- umap_df_3d_noise_removed_with_genes$TPSB2 + 1

# First plot with Wilcoxon test
p <- ggplot(umap_df_3d_noise_removed_with_genes, aes(x = factor(Cluster), y = TPSB2_log, fill = factor(Cluster))) +
  geom_boxplot() +
  scale_fill_manual(values = color_palette) +
  labs(x = "Cluster", y = "Expressin Level (log10 scale)", title = "TPSB2", fill = "Cluster") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 16)
  ) +
  scale_y_log10() +
  geom_signif(
    comparisons = list(c("2","1"), c("2","3"), c("2","4"), c("2","5"), c("2","6"), c("2","7"), c("2","8")),
    test = "wilcox.test",
    map_signif_level = FALSE,
    step_increase = 0.05,
    test.args = list(exact = FALSE) # to handle ties gracefully
  )

ggsave(
  "/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/FOS_ManuscriptSubmission_2024/Figures/TPSB2/plot_TPSB2_log10_wilcox_deseq_normalized.pdf",
  plot = p, width = 10, height = 9, dpi = 300, bg = "white"
)

# library(aod)  # for wald.test, must be load if it is not loaded

# Custom Wald test function
wald_test_func <- function(x, y, ...) {
  val <- c(x, y)
  group <- factor(c(rep("Group1", length(x)), rep("Group2", length(y))))
  # Simple linear model
  model <- lm(val ~ group)
  coefs <- coef(model)
  vc <- vcov(model)
  wald_res <- wald.test(Sigma = vc, b = coefs, Terms = 2)
  p_value <- wald_res$result$chi2["P"]
  list(p.value = p_value)
}

my_comparisons <- list(c("2","1"), c("2","3"), c("2","4"), c("2","5"), c("2","6"), c("2","7"), c("2","8"))

# Plot with Wald test
p <- ggplot(umap_df_3d_noise_removed_with_genes, aes(x = factor(Cluster), y = TPSB2_log, fill = factor(Cluster))) +
  geom_boxplot() +
  scale_fill_manual(values = color_palette) +
  labs(x = "Cluster", y = "Expressin Level (log10 scale)", title = "TPSB2", fill = "Cluster") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 16)
  ) +
  scale_y_log10() +
  geom_signif(
    comparisons = my_comparisons,
    test = wald_test_func,
    map_signif_level = FALSE,
    step_increase = 0.05
  )

ggsave(
  "/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/FOS_ManuscriptSubmission_2024/Figures/TPSB2/plot_TPSB2_log10_wald_deseq_normalized.pdf",
  plot = p, width = 10, height = 9, dpi = 300, bg = "white"
)


# ARHGAP36 Wilcoxon Test
# Add pseudocount to avoid log10(0)
umap_df_3d_noise_removed_with_genes$ARHGAP36_log <- umap_df_3d_noise_removed_with_genes$ARHGAP36 + 1

p <- ggplot(umap_df_3d_noise_removed_with_genes, aes(x = factor(Cluster), y = ARHGAP36_log, fill = factor(Cluster))) +
  geom_boxplot() +
  scale_fill_manual(values = color_palette) +
  labs(x = "Cluster", y = "Expressin Level (log10 scale)", title = "ARHGAP36", fill = "Cluster") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 16)
  ) +
  scale_y_log10() +
  geom_signif(
    comparisons = list(c("2","1"), c("2","3"), c("2","4"), c("2","5"), c("2","6"), c("2","7"), c("2","8")),
    test = "wilcox.test",
    map_signif_level = FALSE,
    step_increase = 0.05,
    test.args = list(exact = FALSE)  # Handle ties gracefully
  )

ggsave("/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/FOS_ManuscriptSubmission_2024/Figures/ARHGAP36/plot_ARHGAP36_log10_wilcox_deseq_normalized.pdf", plot = p, width = 10, height = 9, dpi = 300, bg = "white")

p <- ggplot(umap_df_3d_noise_removed_with_genes, aes(x = factor(Cluster), y = ARHGAP36_log, fill = factor(Cluster))) +
  geom_boxplot() +
  scale_fill_manual(values = color_palette) +
  labs(x = "Cluster", y = "Expressin Level (log10 scale)", title = "ARHGAP36", fill = "Cluster") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 16)
  ) +
  scale_y_log10() +
  geom_signif(
    comparisons = my_comparisons,
    test = wald_test_func,  # Already defined above
    map_signif_level = FALSE,
    step_increase = 0.05
  )

ggsave("/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/FOS_ManuscriptSubmission_2024/Figures/ARHGAP36/plot_ARHGAP36_log10_wald_deseq_normalized.pdf", plot = p, width = 10, height = 9, dpi = 300, bg = "white")


# TPSAB1 log10 - wilcox 
# Add pseudocount to avoid log10(0)
umap_df_3d_noise_removed_with_genes$TPSAB1_log <- umap_df_3d_noise_removed_with_genes$TPSAB1 + 1

p <- ggplot(umap_df_3d_noise_removed_with_genes, aes(x = factor(Cluster), y = TPSAB1_log, fill = factor(Cluster))) +
  geom_boxplot() +
  scale_fill_manual(values = color_palette) +
  labs(x = "Cluster", y = "Expressin Level (log10 scale)", title = "TPSAB1", fill = "Cluster") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 16)
  ) +
  scale_y_log10() +
  geom_signif(
    comparisons = my_comparisons,
    test = "wilcox.test",
    map_signif_level = FALSE,
    step_increase = 0.05,
    test.args = list(exact = FALSE)
  )

ggsave("/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/FOS_ManuscriptSubmission_2024/Figures/TPSAB1/plot_TPSAB1_log10_wilcox_deseq_normalized.pdf", plot = p, width = 10, height = 9, dpi = 300, bg = "white")


# TPSAB1 - log10 - wald
p <- ggplot(umap_df_3d_noise_removed_with_genes, aes(x = factor(Cluster), y = TPSAB1_log, fill = factor(Cluster))) +
  geom_boxplot() +
  scale_fill_manual(values = color_palette) +
  labs(x = "Cluster", y = "Expressin Level (log10 scale)", title = "TPSAB1", fill = "Cluster") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 16)
  ) +
  scale_y_log10() +
  geom_signif(
    comparisons = my_comparisons,
    test = wald_test_func,
    map_signif_level = FALSE,
    step_increase = 0.05
  )

ggsave("/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/FOS_ManuscriptSubmission_2024/Figures/TPSAB1/plot_TPSAB1_log10_wald_deseq_normalized.pdf", plot = p, width = 10, height = 9, dpi = 300, bg = "white")

### lfcShrink ### 

resultsNames(deseq_dbscan)

res_2_vs_all <- results(deseq_dbscan, name = "group_Cluster_2_vs_Rest")
library(apeglm) # ensure apeglm is installed and loaded
res_2_vs_all_shrunk <- lfcShrink(deseq_dbscan, 
                                 coef="group_Cluster_2_vs_Rest", 
                                 type="apeglm",
                                lfcThreshold=0.1)

head(res_2_vs_all_shrunk)

# Remove version numbers from row names
rownames(res_2_vs_all_shrunk) <- sub("\\..*", "", rownames(res_2_vs_all_shrunk))

# map the cleaned Ensembl IDs to gene symbols
res_2_vs_all_shrunk$gene_symbol <- mapIds(
  EnsDb.Hsapiens.v86,
  keys = rownames(res_2_vs_all_shrunk),
  column = "SYMBOL",
  keytype = "GENEID",
  multiVals = "first"
)

write.csv(as.data.frame(res_2_vs_all_shrunk), "/Users/hasanalanya/Desktop/ErsonLab-Yale/FOS-Meningioma/BulkRNA/FOS_ManuscriptSubmission_2024/Tables/res_2_vs_all_shrunk.csv")

