#!/usr/bin/env Rscript
# Survival Analysis — linking gene expression to patient outcomes
# Usage: Rscript ~/tcga-hnsc-pipeline/survival_analysis.R
# Author: Sohel Ahmed, Ph.D.

cat("=== TCGA-HNSC Survival Analysis ===\n\n")

# Install survival packages if needed
for (pkg in c("survival", "survminer")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org")
    }
}

library(survival)
library(survminer)
library(DESeq2)
library(SummarizedExperiment)
library(org.Hs.eg.db)
library(ggplot2)

RESULTS_DIR <- path.expand("~/tcga-hnsc-pipeline/results")
DATA_DIR    <- file.path(RESULTS_DIR, "data")
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
SURV_DIR    <- file.path(RESULTS_DIR, "survival")
dir.create(SURV_DIR, recursive = TRUE, showWarnings = FALSE)

# Load TCGA data (has clinical info in colData)
cat("Loading TCGA data...\n")
data <- readRDS(file.path(DATA_DIR, "tcga_hnsc_data.rds"))

# Load DESeq2 results and VST matrix
res <- read.csv(file.path(RESULTS_DIR, "deseq2/deseq2_results.csv"), row.names = 1)
vsd <- readRDS(file.path(RESULTS_DIR, "deseq2/vst_matrix.rds"))

# Extract clinical data
clin <- as.data.frame(colData(data))

# Check available survival columns
cat("Available survival fields:\n")
surv_cols <- grep("days_to_death|days_to_last_follow|vital_status|OS.time|OS",
                  colnames(clin), value = TRUE)
cat(paste(surv_cols, collapse = ", "), "\n\n")

# Build survival data — use days_to_death and days_to_last_follow_up
# Keep only primary tumors
tumor_idx <- which(clin$definition == "Primary solid Tumor")
clin_tumor <- clin[tumor_idx, ]

# Construct overall survival time and event
clin_tumor$os_time <- as.numeric(clin_tumor$days_to_death)
clin_tumor$os_event <- ifelse(clin_tumor$vital_status == "Dead", 1, 0)

# For alive patients, use days_to_last_follow_up
alive_idx <- which(is.na(clin_tumor$os_time) | clin_tumor$os_event == 0)
clin_tumor$os_time[alive_idx] <- as.numeric(clin_tumor$days_to_last_follow_up[alive_idx])

# Remove patients with missing survival data
valid <- which(!is.na(clin_tumor$os_time) & clin_tumor$os_time > 0)
clin_surv <- clin_tumor[valid, ]
cat("Patients with valid survival data:", nrow(clin_surv), "\n")

# Convert to months
clin_surv$os_months <- clin_surv$os_time / 30.44

# --- Analysis 1: Overall Kaplan-Meier ---
cat("\n--- Overall Kaplan-Meier Curve ---\n")

fit_overall <- survfit(Surv(os_months, os_event) ~ 1, data = clin_surv)
cat("Median OS:", round(summary(fit_overall)$table["median"], 1), "months\n")

png(file.path(FIGURES_DIR, "km_overall.png"), width = 8, height = 6, units = "in", res = 300)
print(ggsurvplot(fit_overall, data = clin_surv,
    title = "TCGA-HNSC: Overall Survival",
    xlab = "Time (months)", ylab = "Survival Probability",
    risk.table = TRUE, risk.table.col = "strata",
    conf.int = TRUE, pval = FALSE,
    ggtheme = theme_minimal(base_size = 14)))
dev.off()
cat("Saved: km_overall.png\n")

# --- Analysis 2: Survival by tumor stage ---
cat("\n--- Survival by Tumor Stage ---\n")

clin_surv$stage_group <- NA
clin_surv$stage_group[grep("Stage I$|Stage IA|Stage IB", clin_surv$ajcc_pathologic_stage)] <- "Stage I-II"
clin_surv$stage_group[grep("Stage II$|Stage IIA|Stage IIB", clin_surv$ajcc_pathologic_stage)] <- "Stage I-II"
clin_surv$stage_group[grep("Stage III$|Stage IIIA|Stage IIIB|Stage IIIC", clin_surv$ajcc_pathologic_stage)] <- "Stage III-IV"
clin_surv$stage_group[grep("Stage IV$|Stage IVA|Stage IVB|Stage IVC", clin_surv$ajcc_pathologic_stage)] <- "Stage III-IV"

stage_valid <- which(!is.na(clin_surv$stage_group))
clin_stage <- clin_surv[stage_valid, ]
cat("Patients with stage info:", nrow(clin_stage), "\n")
cat("Stage I-II:", sum(clin_stage$stage_group == "Stage I-II"),
    " Stage III-IV:", sum(clin_stage$stage_group == "Stage III-IV"), "\n")

fit_stage <- survfit(Surv(os_months, os_event) ~ stage_group, data = clin_stage)

png(file.path(FIGURES_DIR, "km_by_stage.png"), width = 10, height = 7, units = "in", res = 300)
print(ggsurvplot(fit_stage, data = clin_stage,
    title = "TCGA-HNSC: Survival by Tumor Stage",
    xlab = "Time (months)", ylab = "Survival Probability",
    risk.table = TRUE, pval = TRUE,
    palette = c("#2166AC", "#B2182B"),
    legend.labs = c("Stage I-II", "Stage III-IV"),
    ggtheme = theme_minimal(base_size = 14)))
dev.off()
cat("Saved: km_by_stage.png\n")

# --- Analysis 3: Survival by top DE gene expression (high vs low) ---
cat("\n--- Survival by Gene Expression ---\n")

# Pick top 5 DE genes for survival analysis
top5_genes <- head(res[order(res$padj), "symbol"], 5)
top5_genes <- top5_genes[!is.na(top5_genes)]
cat("Analyzing genes:", paste(top5_genes, collapse = ", "), "\n")

# Map gene symbols to Ensembl IDs in VST matrix
gene_ids_vst <- sub("\\..*", "", rownames(vsd))
symbol_map <- mapIds(org.Hs.eg.db, keys = gene_ids_vst,
                     keytype = "ENSEMBL", column = "SYMBOL",
                     multiVals = "first")

# Match VST samples to survival patients
common_samples <- intersect(colnames(vsd), rownames(clin_surv))
cat("Samples with both expression and survival:", length(common_samples), "\n")

surv_results <- data.frame(
    gene = character(),
    hr = numeric(),
    hr_lower = numeric(),
    hr_upper = numeric(),
    pvalue = numeric(),
    stringsAsFactors = FALSE
)

for (gene in top5_genes) {
    gene_idx <- which(symbol_map == gene)[1]
    if (is.na(gene_idx)) next

    expr <- assay(vsd)[gene_idx, common_samples]
    median_expr <- median(expr)

    df <- data.frame(
        os_months = clin_surv[common_samples, "os_months"],
        os_event  = clin_surv[common_samples, "os_event"],
        group     = ifelse(expr > median_expr, "High", "Low")
    )
    df <- df[complete.cases(df), ]

    if (nrow(df) < 20) next

    # Kaplan-Meier plot
    fit_gene <- survfit(Surv(os_months, os_event) ~ group, data = df)

    png(file.path(SURV_DIR, paste0("km_", gene, ".png")),
        width = 10, height = 7, units = "in", res = 300)
    print(ggsurvplot(fit_gene, data = df,
        title = paste("TCGA-HNSC: Survival by", gene, "Expression"),
        xlab = "Time (months)", ylab = "Survival Probability",
        risk.table = TRUE, pval = TRUE,
        palette = c("#B2182B", "#2166AC"),
        legend.labs = c(paste(gene, "High"), paste(gene, "Low")),
        ggtheme = theme_minimal(base_size = 14)))
    dev.off()

    # Cox regression
    cox_fit <- coxph(Surv(os_months, os_event) ~ group, data = df)
    cox_sum <- summary(cox_fit)

    surv_results <- rbind(surv_results, data.frame(
        gene     = gene,
        hr       = round(cox_sum$conf.int[1, 1], 3),
        hr_lower = round(cox_sum$conf.int[1, 3], 3),
        hr_upper = round(cox_sum$conf.int[1, 4], 3),
        pvalue   = signif(cox_sum$coefficients[1, 5], 3)
    ))

    cat("  ", gene, "- HR:", round(cox_sum$conf.int[1, 1], 3),
        "p =", signif(cox_sum$coefficients[1, 5], 3), "\n")
}

# Save Cox results
write.csv(surv_results, file.path(SURV_DIR, "cox_regression_results.csv"), row.names = FALSE)

# --- Analysis 4: Forest plot of hazard ratios ---
cat("\n--- Forest Plot ---\n")

if (nrow(surv_results) >= 2) {
    surv_results$gene <- factor(surv_results$gene, levels = rev(surv_results$gene))

    png(file.path(FIGURES_DIR, "forest_plot_hr.png"), width = 10, height = 6, units = "in", res = 300)
    print(
        ggplot(surv_results, aes(x = hr, y = gene)) +
            geom_point(size = 3) +
            geom_errorbarh(aes(xmin = hr_lower, xmax = hr_upper), height = 0.2) +
            geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
            labs(title = "TCGA-HNSC: Hazard Ratios (High vs Low Expression)",
                 x = "Hazard Ratio (95% CI)", y = "") +
            theme_minimal(base_size = 14) +
            annotate("text", x = max(surv_results$hr_upper) * 1.1,
                     y = surv_results$gene,
                     label = paste("p =", surv_results$pvalue),
                     hjust = 0, size = 3.5)
    )
    dev.off()
    cat("Saved: forest_plot_hr.png\n")
}

cat("\n============================================\n")
cat("Survival analysis complete!\n")
cat("Outputs in:\n")
cat("  ", FIGURES_DIR, " (KM curves, forest plot)\n")
cat("  ", SURV_DIR, " (per-gene KM curves, Cox results)\n")
cat("============================================\n")
