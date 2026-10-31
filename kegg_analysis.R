#!/usr/bin/env Rscript
# KEGG Pathway Enrichment Analysis
# Usage: Rscript ~/tcga-hnsc-pipeline/kegg_analysis.R
# Author: Sohel Ahmed, Ph.D.

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

RESULTS_DIR <- path.expand("~/tcga-hnsc-pipeline/results")
ENRICH_DIR  <- file.path(RESULTS_DIR, "enrichment")
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")

cat("=== KEGG Pathway Enrichment Analysis ===\n\n")

res <- read.csv(file.path(RESULTS_DIR, "deseq2/deseq2_results.csv"), row.names = 1)

# Convert Ensembl IDs to Entrez IDs (KEGG requires Entrez)
gene_ids <- sub("\\..*", "", rownames(res))
entrez_map <- mapIds(org.Hs.eg.db, keys = gene_ids,
                     keytype = "ENSEMBL", column = "ENTREZID",
                     multiVals = "first")
res$entrezid <- entrez_map[gene_ids]

# Upregulated genes
sig_up <- res[which(res$padj < 0.05 & res$log2FoldChange > 1), ]
up_entrez <- na.omit(sig_up$entrezid)
cat("Upregulated genes with Entrez IDs:", length(up_entrez), "\n")

kegg_up <- enrichKEGG(
    gene     = up_entrez,
    organism = "hsa",
    pvalueCutoff = 0.05
)
write.csv(as.data.frame(kegg_up), file.path(ENRICH_DIR, "kegg_enrichment_up.csv"))

png(file.path(FIGURES_DIR, "kegg_dotplot_up.png"), width = 12, height = 8, units = "in", res = 300)
print(dotplot(kegg_up, showCategory = 20,
              title = "KEGG Pathways (Upregulated in Tumor)"))
dev.off()
cat("Saved KEGG upregulated dotplot\n")

# Downregulated genes
sig_down <- res[which(res$padj < 0.05 & res$log2FoldChange < -1), ]
down_entrez <- na.omit(sig_down$entrezid)
cat("Downregulated genes with Entrez IDs:", length(down_entrez), "\n")

kegg_down <- enrichKEGG(
    gene     = down_entrez,
    organism = "hsa",
    pvalueCutoff = 0.05
)
write.csv(as.data.frame(kegg_down), file.path(ENRICH_DIR, "kegg_enrichment_down.csv"))

png(file.path(FIGURES_DIR, "kegg_dotplot_down.png"), width = 12, height = 8, units = "in", res = 300)
print(dotplot(kegg_down, showCategory = 20,
              title = "KEGG Pathways (Downregulated in Tumor)"))
dev.off()
cat("Saved KEGG downregulated dotplot\n")

# Barplot comparison
png(file.path(FIGURES_DIR, "kegg_barplot_up.png"), width = 12, height = 8, units = "in", res = 300)
print(barplot(kegg_up, showCategory = 15,
              title = "KEGG Pathways (Upregulated in HNSC Tumors)"))
dev.off()

# Print top pathways
cat("\n--- Top 10 Upregulated KEGG Pathways ---\n")
if (nrow(as.data.frame(kegg_up)) > 0) {
    top_up <- head(as.data.frame(kegg_up)[, c("Description", "Count", "p.adjust")], 10)
    print(top_up)
}

cat("\n--- Top 10 Downregulated KEGG Pathways ---\n")
if (nrow(as.data.frame(kegg_down)) > 0) {
    top_down <- head(as.data.frame(kegg_down)[, c("Description", "Count", "p.adjust")], 10)
    print(top_down)
}

cat("\n=== KEGG Analysis Complete ===\n")
