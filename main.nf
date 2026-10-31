#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.project     = "TCGA-HNSC"
params.outdir      = "${launchDir}/results"
params.padj_cutoff = 0.05
params.lfc_cutoff  = 1.0
params.top_n_genes = 50

log.info """
    =============================================
    TCGA RNA-seq DE Analysis Pipeline
    =============================================
    Project       : ${params.project}
    Output dir    : ${params.outdir}
    padj cutoff   : ${params.padj_cutoff}
    LFC cutoff    : ${params.lfc_cutoff}
    Top N genes   : ${params.top_n_genes}
    =============================================
"""

process DOWNLOAD_TCGA {
    publishDir "${params.outdir}/data", mode: 'copy'
    output: path "tcga_data.rds", emit: data
    script:
    """
    #!/usr/bin/env Rscript
    library(TCGAbiolinks)
    query <- GDCquery(project="${params.project}", data.category="Transcriptome Profiling",
        data.type="Gene Expression Quantification", workflow.type="STAR - Counts")
    GDCdownload(query)
    data <- GDCprepare(query)
    saveRDS(data, "tcga_data.rds")
    """
}

process RUN_DESEQ2 {
    publishDir "${params.outdir}/deseq2", mode: 'copy'
    input: path data_rds
    output:
        path "deseq2_results.csv", emit: results
        path "deseq2_object.rds",  emit: dds
        path "vst_matrix.rds",     emit: vst
        path "deseq2_summary.txt", emit: summary
    script:
    """
    #!/usr/bin/env Rscript
    library(DESeq2); library(SummarizedExperiment); library(org.Hs.eg.db)
    data <- readRDS("${data_rds}")
    data_filt <- data[, data\$definition %in% c("Primary solid Tumor","Solid Tissue Normal")]
    condition <- factor(ifelse(data_filt\$definition=="Primary solid Tumor","Tumor","Normal"), levels=c("Normal","Tumor"))
    counts_mat <- assay(data_filt, "unstranded")
    keep <- rowSums(counts_mat >= 10) >= 20; counts_mat <- counts_mat[keep, ]
    dds <- DESeqDataSetFromMatrix(countData=counts_mat, colData=data.frame(condition=condition), design=~condition)
    dds <- DESeq(dds)
    res <- results(dds, contrast=c("condition","Tumor","Normal"))
    gene_ids <- sub("\\\\..*","",rownames(res))
    res\$symbol <- mapIds(org.Hs.eg.db, keys=gene_ids, keytype="ENSEMBL", column="SYMBOL", multiVals="first")
    res_ordered <- res[order(res\$padj),]
    write.csv(as.data.frame(res_ordered), "deseq2_results.csv")
    saveRDS(dds, "deseq2_object.rds")
    vsd <- vst(dds, blind=FALSE); vsd\$condition <- condition; saveRDS(vsd, "vst_matrix.rds")
    sink("deseq2_summary.txt")
    cat("DESeq2 Summary\\nGenes tested:", nrow(res), "\\n")
    cat("Significant (padj<${params.padj_cutoff}):", sum(res\$padj<${params.padj_cutoff}, na.rm=TRUE), "\\n")
    sink()
    """
}

process GO_ENRICHMENT {
    publishDir "${params.outdir}/enrichment", mode: 'copy'
    input: path results_csv
    output:
        path "go_enrichment_up.csv",   emit: go_up
        path "go_enrichment_down.csv", emit: go_down
        path "go_dotplot_up.png",      emit: plot_up
        path "go_dotplot_down.png",    emit: plot_down
    script:
    """
    #!/usr/bin/env Rscript
    library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2)
    res <- read.csv("${results_csv}", row.names=1)
    sig_up <- res[which(res\$padj < ${params.padj_cutoff} & res\$log2FoldChange > ${params.lfc_cutoff}),]
    ego_up <- enrichGO(gene=sub("\\\\..*","",rownames(sig_up)), OrgDb=org.Hs.eg.db, keyType="ENSEMBL", ont="BP", pvalueCutoff=0.05)
    write.csv(as.data.frame(ego_up), "go_enrichment_up.csv")
    png("go_dotplot_up.png", width=12, height=8, units="in", res=300)
    print(dotplot(ego_up, showCategory=20, title="GO BP (Upregulated)")); dev.off()
    sig_down <- res[which(res\$padj < ${params.padj_cutoff} & res\$log2FoldChange < -${params.lfc_cutoff}),]
    ego_down <- enrichGO(gene=sub("\\\\..*","",rownames(sig_down)), OrgDb=org.Hs.eg.db, keyType="ENSEMBL", ont="BP", pvalueCutoff=0.05)
    write.csv(as.data.frame(ego_down), "go_enrichment_down.csv")
    png("go_dotplot_down.png", width=12, height=8, units="in", res=300)
    print(dotplot(ego_down, showCategory=20, title="GO BP (Downregulated)")); dev.off()
    """
}

process GENERATE_FIGURES {
    publishDir "${params.outdir}/figures", mode: 'copy'
    input: path results_csv; path vst_rds
    output: path "volcano_plot.png", emit: volcano; path "heatmap_top50.png", emit: heatmap; path "pca_plot.png", emit: pca
    script:
    """
    #!/usr/bin/env Rscript
    library(DESeq2); library(EnhancedVolcano); library(pheatmap); library(org.Hs.eg.db); library(ggplot2)
    res <- read.csv("${results_csv}", row.names=1); vsd <- readRDS("${vst_rds}")
    png("volcano_plot.png", width=10, height=8, units="in", res=300)
    EnhancedVolcano(res, lab=res\$symbol, x='log2FoldChange', y='padj',
        title='TCGA-HNSC: Tumor vs Normal', pCutoff=${params.padj_cutoff}, FCcutoff=${params.lfc_cutoff}); dev.off()
    top_genes <- rownames(res)[order(res\$padj)][1:${params.top_n_genes}]
    mat <- assay(vsd)[top_genes,]; mat <- mat - rowMeans(mat)
    ids <- sub("\\\\..*","",rownames(mat))
    gm <- mapIds(org.Hs.eg.db, keys=ids, keytype="ENSEMBL", column="SYMBOL", multiVals="first")
    rownames(mat) <- ifelse(is.na(gm[ids]), rownames(mat), gm[ids])
    anno_col <- data.frame(Condition=factor(ifelse(vsd\$condition=="Tumor","Tumor","Normal")))
    rownames(anno_col) <- colnames(mat)
    png("heatmap_top50.png", width=14, height=10, units="in", res=300)
    pheatmap(mat, annotation_col=anno_col, show_colnames=FALSE, main="Top ${params.top_n_genes} DE Genes"); dev.off()
    pca_data <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
    pv <- round(100*attr(pca_data,"percentVar"))
    png("pca_plot.png", width=8, height=6, units="in", res=300)
    print(ggplot(pca_data, aes(PC1,PC2,color=condition))+geom_point(size=3,alpha=0.7)+
        xlab(paste0("PC1: ",pv[1],"%"))+ylab(paste0("PC2: ",pv[2],"%"))+theme_minimal()); dev.off()
    """
}

process GENERATE_REPORT {
    publishDir "${params.outdir}", mode: 'copy'
    input: path summary; path go_up; path go_down
    output: path "pipeline_report.txt", emit: report
    script:
    """
    #!/usr/bin/env Rscript
    cat("Pipeline Report\\n===============\\n\\n", file="pipeline_report.txt")
    cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n\\n", file="pipeline_report.txt", append=TRUE)
    cat(readLines("${summary}"), sep="\\n", file="pipeline_report.txt", append=TRUE)
    """
}

workflow {
    DOWNLOAD_TCGA()
    RUN_DESEQ2(DOWNLOAD_TCGA.out.data)
    GO_ENRICHMENT(RUN_DESEQ2.out.results)
    GENERATE_FIGURES(RUN_DESEQ2.out.results, RUN_DESEQ2.out.vst)
    GENERATE_REPORT(RUN_DESEQ2.out.summary, GO_ENRICHMENT.out.go_up, GO_ENRICHMENT.out.go_down)
}
