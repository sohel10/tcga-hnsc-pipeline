# TCGA-HNSC RNA-seq Differential Expression Pipeline

Automated, reproducible pipeline for differential gene expression analysis of head and neck squamous cell carcinoma (HNSC) using TCGA RNA-seq data. Built with Nextflow and Bioconductor.

## Overview

Downloads TCGA-HNSC RNA-seq data from the GDC Data Portal, performs differential expression analysis between primary tumors (n=520) and matched normal tissues (n=44), runs Gene Ontology enrichment, and generates publication-quality figures and an interactive Shiny dashboard.

## Quick Start

### Run the analysis directly

```bash
cd ~/tcga-hnsc-pipeline
Rscript run_analysis.R
```

### Run with Nextflow

```bash
curl -s https://get.nextflow.io | bash
./nextflow run main.nf
```

### Run with Docker

```bash
docker build -t sohelahmed/tcga-hnsc-pipeline:latest .
./nextflow run main.nf -profile docker
```

## Output Structure
results/
data/tcga_hnsc_data.rds # Cached TCGA download
deseq2/deseq2_results.csv # Full DE results with gene symbols
deseq2/vst_matrix.rds # Variance-stabilized counts
figures/volcano_plot.png # Volcano plot
figures/heatmap_top50.png # Heatmap of top 50 DE genes
figures/pca_plot.png # PCA sample separation
enrichment/go_enrichment_up.csv # GO terms (upregulated)
enrichment/go_enrichment_down.csv # GO terms (downregulated)

## Interactive Shiny Dashboard

```bash
cd ~/tcga-hnsc-pipeline/shiny
Rscript -e "shiny::runApp('.', port = 3838)"
```

Open http://localhost:3838 in your browser.

## Key Findings

- **27,556 genes** tested after filtering
- **11,531 upregulated** and **8,127 downregulated** genes (padj < 0.1)
- Top upregulated: CA9 (hypoxia marker), MMP11 (matrix metalloproteinase), IL11
- Top downregulated: CST4 (cystatin S), CCL28 (chemokine), PLIN1

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| --padj_cutoff | 0.05 | Adjusted p-value threshold |
| --lfc_cutoff | 1.0 | Log2 fold change threshold |
| --top_n_genes | 50 | Number of top genes for heatmap |

## Author

Sohel Ahmed, Ph.D.

## License

MIT
