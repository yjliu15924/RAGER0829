# **RAGER table of contents**
1. [Quick start](https://github.com/yjliu15924/RAGER/blob/main/RAGER_github/README.md#quick-start)
2. [Preprocess RNAseq data](https://github.com/yjliu15924/RAGER/blob/main/RAGER_github/Scripts/Preprocess_RNAseq/RNAseq_analysis.md)
3. [Preprocess ATACseq data](https://github.com/yjliu15924/RAGER/blob/main/RAGER_github/Scripts/Preprocess_ATACseq/ATACseq_analysis.md) 
4. [Joint analysis](https://github.com/yjliu15924/RAGER/blob/main/RAGER_github/Scripts/Joint_analysis/Joint_analysis.md)
5. [Custom analysis](https://github.com/yjliu15924/RAGER/blob/main/RAGER_github/Scripts/Custom_analysis/Custom_analysis.md)

## **Quick Start**

Welcome to RAGER! This section will guide through the essential steps to configure and run custom genes analysis pipeline.

Note that the custom analysis is based on the `preprocess RNAseq data`, `preprocess ATACseq data` and `Joint analysis`.
### **Step 1: Open the Configuration File**

Open the `config.yaml` file in the species-specific directory:

```bash
cd ~/PROJECT/RAGER  #Change the directory to PROJECT/RAGER as the current directory
#mouse
vim ./mouse/scripts/snakemake/Custom_genes_analysis/config.yaml  
#human
vim ./human/scripts/snakemake/Custom_genes_analysis/config.yaml 
# you can use any text editor 
```

### **Step 2: Modifying the Configuration File**

The `config.yaml` file contains all parameters needed to customize RAGER for custom genes analysis. Below is a detailed explanation of each section:

#### **Input Files Section**
Specify the input files for custom genes analysis:

```yaml
input_files:
  custom_genes: "datasets/Custom_genes_analysis/custom_genes.txt"    # Custom gene list file
  deg_csv: "datasets/RNAseq/stringtiefile/A_vs_N_DEG.csv"          # Differential expressed genes file
```

**How to modify**: 
- Update `custom_genes` path to point to your gene list file (one gene symbol per line)
- Update `deg_csv` path to match your RNA-seq differential expression analysis output
- Ensure both files are accessible and properly formatted

#### **Paths Section**
This section defines the directory structure for custom genes analysis:

```yaml
paths:
  dir: "datasets/Custom_genes_analysis"                    # Main analysis directory
  scripts_dir: "scripts/snakemake/Custom_genes_analysis/"  # Scripts directory
  
  jaspar_dir: "datasets/JASPAR2024_CORE"                  # JASPAR database directory
  tf_mapping_file: "datasets/Mus_musculus_TF.txt"         # Transcription factor mapping file
```

**How to modify**: 
- All directories will be created by the pipeline if they don't exist
- Update `tf_mapping_file` for different species (e.g.`datasets/Homo_sapiens_TF.txt` for human)

#### **Parameters Section**

**Differential Expression Thresholds:**
```yaml
params:
  log2fc_up: 1                    # Log2 fold change threshold for upregulation
  log2fc_down: -1                 # Log2 fold change threshold for downregulation
  padj: 0.05                      # Adjusted p-value threshold
```

**How to modify**: 
- Adjust `log2fc_up` and `log2fc_down` for more or less stringent fold change requirements
- Modify `padj` threshold for statistical significance (lower = more stringent)
- These parameters will be used to filter the differential expression results for enrichment TFs.

### **Step 3: Prepare Custom Gene List**

Create your custom gene list file:

```bash
# Create the custom genes file
mkdir -p datasets/Custom_genes_analysis
vim datasets/Custom_genes_analysis/custom_genes.txt
```

The custom genes file should contain one gene symbol per line:
```
AC016739
ACTR2
AHSP
ALAS2
APOE
ATP6V1E1
```

**File format requirements**:
- One gene symbol per line
- No headers or additional columns
- UTF-8 text encoding

### **Step 4: Validate Configuration**

Before running the pipeline, ensure:

1. Custom gene list file exists and is properly formatted
2. DEG CSV file from RNA-seq analysis is accessible
3. JASPAR database and TF mapping files are available
4. Differential expression thresholds are appropriate for your analysis
5. All file paths are correct and accessible
6. Output directory permissions allow write access

### **Step 5: Run the Pipeline**

Once configuration is ready, execute:

```bash
#mouse
cd ~/PROJECT/RAGER/mouse  #Change the directory to PROJECT/RAGER/mouse as the current directory
snakemake --snakefile ./scripts/snakemake/Custom_genes_analysis/Custom_gene_analysis.py --configfile ./scripts/snakemake/Custom_genes_analysis/config.yaml -j 10

#human
cd ~/PROJECT/RAGER/human  #Change the directory to PROJECT/RAGER/human as the current directory
snakemake --snakefile ./scripts/snakemake/Custom_genes_analysis/Custom_gene_analysis.py --configfile ./scripts/snakemake/Custom_genes_analysis/config.yaml -j 10
```

### **Expected Outputs**

The pipeline will generate analysis results in the `datasets/Custom_genes_analysis/` directory, including:

- Filtered differential expression results for your custom genes
- Transcription factor binding analysis
- Regulatory network analysis
- Visualization plots and summary statistics

**Note**: This analysis pipeline focuses on your specific genes of interest and provides targeted regulatory insights based on the custom gene list you provide.

## **List of processes**
- [Analysis_of_the_promoter_region](#)
    - [Get_shared_up_gene_InPromoter](#)
    - [Get_shared_down_gene_InPromoter](#)
    - [Expression_shared_peak_annotation](#)
    - [Calculate_shared_peak_TPM_signal](#)
    - [Plot_RNAseq_ATACseq_correlation](#)
    - [Enrichment_analysis_in_promoter_region](#)
        - [Enrichment_analysis_for_shared_gene](#)
    - [Regulatory_Analysis_In_promoter_region](#)
        - [Format_TF_motifs](#)
        - [Extract_shared_Gene_promoterSeq](#)
        - [Motif_enrich](#)
        - [Constructe_TF_Gene_network](#)
        - [Extract_enriched_motifs](#)
        - [Plot_TF_heatmap]()
- [Analysis_of_the_enhancer_region]()
    - [Get_RNAseq_DE_gene_annotation](#)
    - [Calculate_enhancer_DEGene_distance](#)
    - [Get_enhancer_annotation](#)
    - [Calculate_enhancer_TPM_signal](#)
    - [Get_shared_Up_Enhancer_and_shared_Up_gene](#)
    - [Get_shared_Down_Enhancer_and_shared_Down_gene](#)
    - [Enrichment_analysis_in_enhancer_region](#)
        - [Enrichment_analysis_for_shared_gene](#)
    - [Regulatory_Analysis_in_enhancer_region](#)
        - [Extract_enhancer_sequences](#)
        - [Motif_enrich](#)
        - [Constructe_TF_Gene_network](#)
        - [Extract_enriched_motifs](#)
        - [Plot_TF_heatmap](#)

## **Preprocessing**
### **convert_customGene_to_ensemlID**

**Description**

Converts custom gene names/symbols to Ensembl gene IDs for standardized downstream analysis.

**Inputs**
- `custom_genes.txt`: List of custom gene names/symbols

**Outputs**
- `custom_genes_ensemblID.txt`: Custom genes converted to Ensembl IDs

**Output directory**
- `./datasets/Custom_genes_analysis/`

---

### **GSEA_enrichment**

**Description**

Performs Gene Set Enrichment Analysis (GSEA) for custom genes using differential expression data, identifying enriched GO biological processes and KEGG pathways.

**Inputs**
- `custom_genes_ensemblID.txt`: Custom genes with Ensembl IDs
- `RNA_GDC_vs_RNA_ctr_DEG.csv`: Differential expression gene data

**Outputs**
- `custom_gene_gsea_GOenrich.csv`: GO enrichment analysis results
- `custom_gene_gsea_KEGGenrich.csv`: KEGG pathway enrichment analysis results

**Output directory**
- `./datasets/Custom_genes_analysis/`

---

### **extr_promoterSeq**

**Description**

Extracts promoter sequences for custom genes in FASTA format for downstream motif analysis.

**Inputs**
- `custom_genes_ensemblID.txt`: Custom genes with Ensembl IDs


**Outputs**
- `promoter_seqs.fa`: Extracted promoter sequences in FASTA format

**Output directory**
- `./datasets/Custom_genes_analysis/`

---

### **motif_enrich**

**Description**

Performs transcription factor binding motif enrichment analysis on promoter sequences using JASPAR motif database.

**Inputs**
- `promoter_seqs.fa`: Promoter sequences in FASTA format
- `JASPAR2024_CORE/`: JASPAR motif database directory

**Outputs**
- `motif_enrich_result`:Transcription factor binding motif enrichment result directory

**Output directory**
- `./datasets/Custom_genes_analysis/motif_enrich_result/`

---

### **TF_Gene_network**

**Description**

Constructs transcription factor-gene regulatory network based on motif enrichment analysis results.

**Inputs**
- `motif_enrich_result`: Transcription factor binding motif enrichment result directory
**Outputs**
- `TF_gene_network.txt`: Transcription factor-gene regulatory network

**Output directory**
- `./datasets/Custom_genes_analysis/`

---

### **extract_enriched_motifs**

**Description**

Extracts significantly enriched transcription factor motifs and maps them to corresponding transcription factor names.

**Inputs**
- `motif_enrich_result`: Transcription factor binding motif enrichment result directory
- `Homo_sapiens_TF.txt`: Human transcription factor mapping file

**Parameters**
- `FDR < 0.05`:FDR cutoffs for significance

**Outputs**
- `enriched_tfs_list.txt`: List of significantly enriched transcription factors

**Output directory**
- `./datasets/Custom_genes_analysis/`

---

## **TF_heatmap**

**Description**

Generates expression heatmap for enriched transcription factors using differential expression data with specified fold change and significance thresholds.

**Inputs**
- `enriched_tfs_list.txt`: List of significantly enriched transcription factors
- `RNA_GDC_vs_RNA_ctr_DEG.csv`: Differential expression gene data

**Parameters**
- `log2fc_up`:  Log2 fold change threshold for upregulated genes
- `log2fc_down`: Log2 fold change threshold for downregulated genes
- `padj`: 0.05 (Adjusted p-value threshold for statistical significance)

**Outputs**
- `TF_expr_heatmap.pdf`: Expression heatmap of enriched transcription factors

**Output directory**
- `./datasets/Custom_genes_analysis/`