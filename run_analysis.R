#!/usr/bin/env Rscript

# TCGA-HNSC Differential Gene Expression Analysis
# Usage: Rscript ~/tcga-hnsc-pipeline/run_analysis.R
# (Do NOT use sudo)
# Author: Sohel Ahmed, Ph.D.

PROJECT_DIR <- path.expand("~/tcga-hnsc-pipeline")
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
DATA_DIR    <- file.path(RESULTS_DIR, "data")
DESEQ_DIR   <- file.path(RESULTS_DIR, "deseq2")
ENRICH_DIR  <- file.path(RESULTS_DIR, "enrichment")
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")

PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1.0
TOP_N_GENES <- 50

for (d in c(DATA_DIR, DESEQ_DIR, ENRICH_DIR, FIGURES_DIR)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("All output will be saved under:", RESULTS_DIR, "\n\n")

# Step 1: Load TCGA Data
cat("=== Step 1: Loading TCGA-HNSC data ===\n")
cached_data <- file.path(DATA_DIR, "tcga_hnsc_data.rds")

if (file.exists(cached_data)) {
    cat("Using cached data:", cached_data, "\n")
    data <- readRDS(cached_data)
} else {
    library(TCGAbiolinks)
    query <- GDCquery(
        project = "TCGA-HNSC",
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = "STAR - Counts"
    )
    GDCdownload(query, directory = file.path(DATA_DIR, "GDCdata"))
    data <- GDCprepare(query, directory = file.path(DATA_DIR, "GDCdata"))
    saveRDS(data, cached_data)
    cat("Data saved to:", cached_data, "\n")
}
cat("Samples:", ncol(data), " Genes:", nrow(data), "\n\n")

# Step 2: DESeq2
cat("=== Step 2: Running DESeq2 ===\n")
library(DESeq2)
library(SummarizedExperiment)
library(org.Hs.eg.db)

data_filt <- data[, data$definition %in% c("Primary solid Tumor", "Solid Tissue Normal")]
condition <- factor(
    ifelse(data_filt$definition == "Primary solid Tumor", "Tumor", "Normal"),
    levels = c("Normal", "Tumor")
)

counts_mat <- assay(data_filt, "unstranded")
keep <- rowSums(counts_mat >= 10) >= 20
counts_mat <- counts_mat[keep, ]
cat("Retained", nrow(counts_mat), "genes after filtering\n")

dds <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData   = data.frame(condition = condition),
    design    = ~ condition
)
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "Tumor", "Normal"))
gene_ids <- sub("\\..*", "", rownames(res))
gene_map <- mapIds(org.Hs.eg.db, keys = gene_ids,
                   keytype = "ENSEMBL", column = "SYMBOL",
                   multiVals = "first")
res$symbol <- gene_map[gene_ids]

res_ordered <- res[order(res$padj), ]
results_csv <- file.path(DESEQ_DIR, "deseq2_results.csv")
write.csv(as.data.frame(res_ordered), results_csv)
cat("DESeq2 results saved to:", results_csv, "\n")

saveRDS(dds, file.path(DESEQ_DIR, "deseq2_object.rds"))
vsd <- vst(dds, blind = FALSE)
vsd$condition <- condition
saveRDS(vsd, file.path(DESEQ_DIR, "vst_matrix.rds"))

n_up  <- sum(res$padj < PADJ_CUTOFF & res$log2FoldChange > LFC_CUTOFF, na.rm = TRUE)
n_down <- sum(res$padj < PADJ_CUTOFF & res$log2FoldChange < -LFC_CUTOFF, na.rm = TRUE)
cat("Up:", n_up, "Down:", n_down, "\n\n")

# Step 3: Volcano Plot
cat("=== Step 3: Volcano Plot ===\n")
library(EnhancedVolcano)

volcano_path <- file.path(FIGURES_DIR, "volcano_plot.png")
png(volcano_path, width = 10, height = 8, units = "in", res = 300)
EnhancedVolcano(res,
    lab = res$symbol, x = 'log2FoldChange', y = 'padj',
    title = 'TCGA-HNSC: Tumor vs Normal',
    subtitle = 'Differential Gene Expression (DESeq2)',
    pCutoff = PADJ_CUTOFF, FCcutoff = LFC_CUTOFF,
    pointSize = 1.5, labSize = 3.0)
dev.off()
cat("Saved:", volcano_path, "\n")

# Step 4: Heatmap
cat("=== Step 4: Heatmap ===\n")
library(pheatmap)

top_genes <- rownames(res_ordered)[1:TOP_N_GENES]
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)
top_ids <- sub("\\..*", "", rownames(mat))
top_symbols <- mapIds(org.Hs.eg.db, keys = top_ids,
                      keytype = "ENSEMBL", column = "SYMBOL",
                      multiVals = "first")
rownames(mat) <- ifelse(is.na(top_symbols), rownames(mat), top_symbols)
anno_col <- data.frame(Condition = condition)
rownames(anno_col) <- colnames(mat)

heatmap_path <- file.path(FIGURES_DIR, "heatmap_top50.png")
png(heatmap_path, width = 14, height = 10, units = "in", res = 300)
pheatmap(mat, annotation_col = anno_col, show_colnames = FALSE,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    main = paste("Top", TOP_N_GENES, "DE Genes: TCGA-HNSC"),
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100))
dev.off()
cat("Saved:", heatmap_path, "\n")

# Step 5: PCA Plot
cat("=== Step 5: PCA Plot ===\n")
library(ggplot2)

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

pca_path <- file.path(FIGURES_DIR, "pca_plot.png")
png(pca_path, width = 8, height = 6, units = "in", res = 300)
print(
    ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
        geom_point(size = 3, alpha = 0.7) +
        xlab(paste0("PC1: ", pct_var[1], "% variance")) +
        ylab(paste0("PC2: ", pct_var[2], "% variance")) +
        ggtitle("PCA: TCGA-HNSC Tumor vs Normal") +
        theme_minimal(base_size = 14) +
        scale_color_manual(values = c("Normal" = "#2166AC", "Tumor" = "#B2182B")))
dev.off()
cat("Saved:", pca_path, "\n")

# Step 6: GO Enrichment
cat("=== Step 6: GO Enrichment ===\n")
library(clusterProfiler)

sig_up <- res[which(res$padj < PADJ_CUTOFF & res$log2FoldChange > LFC_CUTOFF), ]
up_genes <- sub("\\..*", "", rownames(sig_up))
ego_up <- enrichGO(gene = up_genes, OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
write.csv(as.data.frame(ego_up), file.path(ENRICH_DIR, "go_enrichment_up.csv"))

png(file.path(ENRICH_DIR, "go_dotplot_up.png"), width = 12, height = 8, units = "in", res = 300)
print(dotplot(ego_up, showCategory = 20, title = "GO BP (Upregulated in Tumor)"))
dev.off()

sig_down <- res[which(res$padj < PADJ_CUTOFF & res$log2FoldChange < -LFC_CUTOFF), ]
down_genes <- sub("\\..*", "", rownames(sig_down))
ego_down <- enrichGO(gene = down_genes, OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
write.csv(as.data.frame(ego_down), file.path(ENRICH_DIR, "go_enrichment_down.csv"))

png(file.path(ENRICH_DIR, "go_dotplot_down.png"), width = 12, height = 8, units = "in", res = 300)
print(dotplot(ego_down, showCategory = 20, title = "GO BP (Downregulated in Tumor)"))
dev.off()

cat("GO enrichment saved to:", ENRICH_DIR, "\n\n")

cat("============================================\n")
cat("DONE! All outputs in:", RESULTS_DIR, "\n")
cat("============================================\n")
