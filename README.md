# **Novel Subgroup of Meningiomas Involving FOS and FOSB Gene Fusions**  
### **Bulk RNA-seq Analysis and Tumor Microenvironment Profiling**

---

## **1. Summary**  
This repository provides **R scripts** and **data** used to identify and characterize a **novel meningioma subgroup** harboring **FOS** and **FOSB** gene fusions. The analyses leverage **bulk RNA-seq** data with a focus on **preprocessing**, **differential gene expression**, and **microenvironment** deconvolution.

---

## **2. Pipeline Outline**  
1. **Data Acquisition**  
   - **Public Datasets**: Downloaded from the Sequence Read Archive (SRA) via **SRA Toolkit**  
     - Baylor College of Medicine (160 cases; GEO accession number: **GSE136661**)
     - University of California San Francisco (185 cases; GEO accession number: **GSE183656**)
     - Palacky University and University Hospital (70 cases; ENA accession number: **PRJNA705586**)
     - Fred Hutchinson Cancer Center (279 cases; GEO accession number: **GSE252291**)
     - Yale University School of Medicine (19 cases; GEO accession number: **GSE85133**)
     - Yale University School of Medicine (23 cases; ENA accession number: **PRJEB55424**)
       
   - **Trimming & Alignment**:  
     - **Trim Galore** for adapter trimming  
     - **STAR** for read alignment  
   - **Gene Quantification**: **HTSeq** (GRCh37)

2. **Batch Correction & Normalization**  
   - **ComBat-Seq** for batch adjustment  
   - **VST** (Variance Stabilizing Transformation)  

3. **Dimensionality Reduction & Clustering**  
   - **PCA** and **UMAP** for data visualization  
   - **dbscan** for unsupervised sample clustering  

4. **Differential Expression & Enrichment**  
   - **DESeq2** to identify key differentially expressed genes  
   - **enrichR** for GO, KEGG, and TF enrichment analyses  

5. **Visualization**  
   - **3D UMAP** plots (mutation status, batch, clusters)  
   - **Volcano/Manhattan** plots for DEGs and enriched pathways  
   - **Boxplots** for selected genes of interest  

6. **Immune & Microenvironment Analysis**  
   - **CIBERSORTx** and **xCell** for cell-type deconvolution (immune/stromal)  
   - **ESTIMATE** for tumor purity, stromal, and immune scores  

---

## **3. Script Descriptions**  
- **`FOS_FOSB_RNAseq_Workflow.R`**  
  - Main pipeline: data import, normalization, UMAP clustering, DE analysis, and result visualization.  

- **`CIBERSORTx_wilcox_test.R`**  
  - Integrates **CIBERSORTx** output to compare immune-cell proportions across clusters.  
  - Performs Wilcoxon tests and creates boxplots showing cell-type distribution.  

- **`xCELL_updated.R`**  
  - Runs **xCell** to estimate relative cell abundances.  
  - Merges these estimates with cluster info for statistical comparisons (Wilcoxon).  

- **`ESTIMATE-tumor-purity-stromal-immune-scoring-tool.R`**  
  - Uses **ESTIMATE** algorithm to infer tumor purity and stromal/immune infiltration.  
  - Offers additional microenvironment insights across subgroups.

---

## **4. Usage**  
1. **Clone/Download** this repo.  
2. **Install** required R packages:  
   - Core (e.g., **DESeq2**, **sva**, **dbscan**, **clusterProfiler**, **ggplot2**, **plotly**, etc.)  
   - Deconvolution tools (e.g., **CIBERSORTx**, **xCell**, **ESTIMATE**)  
3. **Edit paths** in each script to match your file system.  
4. **Run** the scripts in order or selectively, depending on your analysis goals.

---

## **5. References & Tools**  
- **Public Datasets (SRA)**: <https://www.ncbi.nlm.nih.gov/sra>  
- **SRA Toolkit**: <https://github.com/ncbi/sra-tools>  
- **Trim Galore**: <https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/>  
- **STAR**: <https://github.com/alexdobin/STAR>  
- **HTSeq**: <https://htseq.readthedocs.io/en/release_0.11.1/>  
- **CIBERSORTx**: <https://cibersortx.stanford.edu/>  
- **xCell**: <https://xcell.ucsf.edu/>  
- **ESTIMATE**: <https://bioinformatics.mdanderson.org/public-software/estimate/>

---

## **6. Authors & Contact**  
- **Lead Authors**:  
  - 
- **Principal Investigators**:  
  - Zeynep Erson-Omay, PhD  
  - Murat Gunel, MD  

Visit our labs for more details:  
- [Erson Lab](https://ersonlab.org/)  
- [Gunel Lab](https://medicine.yale.edu/lab/gunel/)

For questions or feedback, please contact: **zeynep.erson@yale.edu** or **murat.gunel@yale.edu**
