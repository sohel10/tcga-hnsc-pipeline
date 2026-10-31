#!/usr/bin/env Rscript
library(shiny)
library(DT)
library(ggplot2)
library(plotly)
library(DESeq2)
library(pheatmap)
library(org.Hs.eg.db)

ui <- fluidPage(
    titlePanel("TCGA-HNSC: Differential Gene Expression Explorer"),
    sidebarLayout(
        sidebarPanel(width = 3,
            h4("Filter Parameters"),
            sliderInput("padj", "Adjusted p-value cutoff:", min = 0.001, max = 0.1, value = 0.05, step = 0.001),
            sliderInput("lfc", "Log2 Fold Change cutoff:", min = 0.5, max = 5, value = 1.0, step = 0.25),
            numericInput("top_n", "Top N genes for heatmap:", value = 50, min = 10, max = 200),
            hr(), h4("Search"),
            textInput("gene_search", "Search gene symbol:", placeholder = "e.g., MMP11, CA9"),
            actionButton("search_btn", "Search", class = "btn-primary"),
            hr(), h4("Dataset Summary"), verbatimTextOutput("summary_text"),
            hr(), downloadButton("download_results", "Download Filtered Results")
        ),
        mainPanel(width = 9,
            tabsetPanel(
                tabPanel("Volcano Plot", plotlyOutput("volcano", height = "600px")),
                tabPanel("Heatmap", plotOutput("heatmap", height = "700px")),
                tabPanel("PCA", plotlyOutput("pca", height = "500px")),
                tabPanel("Gene Table", DTOutput("gene_table")),
                tabPanel("Gene Detail",
                    fluidRow(
                        column(6, plotOutput("gene_boxplot", height = "400px")),
                        column(6, plotOutput("gene_counts", height = "400px"))
                    ), verbatimTextOutput("gene_info")
                ),
                tabPanel("GO Enrichment",
                    fluidRow(
                        column(6, h4("Upregulated in Tumor"), DTOutput("go_up_table")),
                        column(6, h4("Downregulated in Tumor"), DTOutput("go_down_table"))
                    )
                )
            )
        )
    )
)

server <- function(input, output, session) {
    data_path <- reactive({
        base <- path.expand("~/tcga-hnsc-pipeline/results")
        if (file.exists(file.path(base, "deseq2/deseq2_results.csv"))) {
            list(results = file.path(base, "deseq2/deseq2_results.csv"),
                 vst = file.path(base, "deseq2/vst_matrix.rds"),
                 go_up = file.path(base, "enrichment/go_enrichment_up.csv"),
                 go_down = file.path(base, "enrichment/go_enrichment_down.csv"))
        } else { stop("No results found. Run run_analysis.R first.") }
    })
    res_data <- reactive({ read.csv(data_path()$results, row.names = 1) })
    vst_data <- reactive({
        p <- data_path()$vst
        if (!is.null(p) && file.exists(p)) readRDS(p) else NULL
    })
    filtered_res <- reactive({
        df <- res_data()
        sig <- df[which(df$padj < input$padj & abs(df$log2FoldChange) > input$lfc), ]
        sig[order(sig$padj), ]
    })
    output$summary_text <- renderText({
        df <- res_data(); sig <- filtered_res()
        n_up <- sum(sig$log2FoldChange > 0, na.rm = TRUE)
        n_down <- sum(sig$log2FoldChange < 0, na.rm = TRUE)
        paste0("Total genes: ", nrow(df), "\nSignificant: ", nrow(sig),
               "\n  Upregulated: ", n_up, "\n  Downregulated: ", n_down)
    })
    output$volcano <- renderPlotly({
        df <- res_data()
        df$neg_log10_padj <- -log10(df$padj)
        df$significance <- "Not Significant"
        df$significance[which(df$padj < input$padj & df$log2FoldChange > input$lfc)] <- "Upregulated"
        df$significance[which(df$padj < input$padj & df$log2FoldChange < -input$lfc)] <- "Downregulated"
        df$label <- ifelse(is.na(df$symbol), rownames(df), df$symbol)
        p <- ggplot(df, aes(x = log2FoldChange, y = neg_log10_padj, color = significance, text = label)) +
            geom_point(alpha = 0.5, size = 1) +
            scale_color_manual(values = c("Upregulated"="#B2182B","Downregulated"="#2166AC","Not Significant"="grey60")) +
            geom_hline(yintercept = -log10(input$padj), linetype = "dashed", color = "grey40") +
            geom_vline(xintercept = c(-input$lfc, input$lfc), linetype = "dashed", color = "grey40") +
            labs(title = "TCGA-HNSC: Tumor vs Normal", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
            theme_minimal(base_size = 14)
        ggplotly(p, tooltip = "text")
    })
    output$heatmap <- renderPlot({
        vsd <- vst_data()
        if (is.null(vsd)) { plot.new(); text(0.5, 0.5, "Run full pipeline first.", cex = 1.5); return() }
        df <- res_data()
        top_genes <- rownames(df)[order(df$padj)][1:input$top_n]
        mat <- assay(vsd)[top_genes, ]; mat <- mat - rowMeans(mat)
        gene_ids <- sub("\\..*", "", rownames(mat))
        gm <- mapIds(org.Hs.eg.db, keys=gene_ids, keytype="ENSEMBL", column="SYMBOL", multiVals="first")
        rownames(mat) <- ifelse(is.na(gm[gene_ids]), rownames(mat), gm[gene_ids])
        anno_col <- data.frame(Condition = vsd$condition); rownames(anno_col) <- colnames(mat)
        pheatmap(mat, annotation_col=anno_col, show_colnames=FALSE,
            clustering_distance_rows="correlation", clustering_distance_cols="correlation",
            main=paste("Top", input$top_n, "DE Genes"),
            color=colorRampPalette(c("#2166AC","white","#B2182B"))(100))
    })
    output$pca <- renderPlotly({
        vsd <- vst_data(); if (is.null(vsd)) return(plotly_empty())
        pca_data <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
        pct_var <- round(100 * attr(pca_data, "percentVar"))
        p <- ggplot(pca_data, aes(x=PC1, y=PC2, color=condition, text=name)) +
            geom_point(size=3, alpha=0.7) +
            xlab(paste0("PC1: ", pct_var[1], "% variance")) +
            ylab(paste0("PC2: ", pct_var[2], "% variance")) +
            scale_color_manual(values=c("Normal"="#2166AC","Tumor"="#B2182B")) +
            theme_minimal(base_size=14) + ggtitle("PCA: Tumor vs Normal")
        ggplotly(p, tooltip="text")
    })
    output$gene_table <- renderDT({
        df <- filtered_res(); df$padj <- signif(df$padj, 3)
        df$log2FoldChange <- round(df$log2FoldChange, 3); df$baseMean <- round(df$baseMean, 1)
        datatable(df[, c("symbol","baseMean","log2FoldChange","padj")],
            selection="single", options=list(pageLength=25, order=list(list(4,'asc'))))
    })
    selected_gene <- reactiveVal(NULL)
    observeEvent(input$search_btn, { selected_gene(toupper(input$gene_search)) })
    observeEvent(input$gene_table_rows_selected, {
        df <- filtered_res(); gene <- df$symbol[input$gene_table_rows_selected]
        if (!is.na(gene)) selected_gene(gene)
    })
    output$gene_boxplot <- renderPlot({
        gene <- selected_gene()
        if (is.null(gene)) { plot.new(); text(0.5, 0.5, "Select a gene", cex=1.2); return() }
        vsd <- vst_data(); if (is.null(vsd)) return()
        gene_ids <- sub("\\..*", "", rownames(vsd))
        gm <- mapIds(org.Hs.eg.db, keys=gene_ids, keytype="ENSEMBL", column="SYMBOL", multiVals="first")
        idx <- which(gm == gene)[1]; if (is.na(idx)) { plot.new(); text(0.5,0.5,paste(gene,"not found"),cex=1.2); return() }
        df_plot <- data.frame(expression=assay(vsd)[idx,], condition=vsd$condition)
        ggplot(df_plot, aes(x=condition, y=expression, fill=condition)) +
            geom_boxplot(alpha=0.7, outlier.shape=NA) + geom_jitter(width=0.2, alpha=0.3, size=1) +
            scale_fill_manual(values=c("Normal"="#2166AC","Tumor"="#B2182B")) +
            labs(title=paste(gene, "Expression"), y="VST Normalized Expression") +
            theme_minimal(base_size=14) + theme(legend.position="none")
    })
    output$gene_counts <- renderPlot({
        gene <- selected_gene(); if (is.null(gene)) return()
        vsd <- vst_data(); if (is.null(vsd)) return()
        gene_ids <- sub("\\..*", "", rownames(vsd))
        gm <- mapIds(org.Hs.eg.db, keys=gene_ids, keytype="ENSEMBL", column="SYMBOL", multiVals="first")
        idx <- which(gm == gene)[1]; if (is.na(idx)) return()
        df_plot <- data.frame(expression=assay(vsd)[idx,], condition=vsd$condition)
        ggplot(df_plot, aes(x=expression, fill=condition)) + geom_density(alpha=0.5) +
            scale_fill_manual(values=c("Normal"="#2166AC","Tumor"="#B2182B")) +
            labs(title=paste(gene, "Distribution"), x="VST Normalized Expression") +
            theme_minimal(base_size=14)
    })
    output$gene_info <- renderText({
        gene <- selected_gene(); if (is.null(gene)) return("")
        df <- res_data(); row <- df[which(df$symbol == gene), ]
        if (nrow(row) == 0) return(paste("Gene", gene, "not found"))
        row <- row[which.min(row$padj), ]
        paste0("Gene: ", gene, "\nEnsembl ID: ", rownames(row), "\nBase Mean: ", round(row$baseMean,1),
            "\nLog2 Fold Change: ", round(row$log2FoldChange,3), "\nAdj P-value: ", signif(row$padj,3),
            "\nDirection: ", ifelse(row$log2FoldChange > 0, "Upregulated in Tumor", "Downregulated in Tumor"))
    })
    output$go_up_table <- renderDT({
        p <- data_path()$go_up; if (is.null(p) || !file.exists(p)) return(data.frame())
        df <- read.csv(p); df$p.adjust <- signif(df$p.adjust, 3)
        datatable(df[, c("Description","Count","p.adjust")], options=list(pageLength=15))
    })
    output$go_down_table <- renderDT({
        p <- data_path()$go_down; if (is.null(p) || !file.exists(p)) return(data.frame())
        df <- read.csv(p); df$p.adjust <- signif(df$p.adjust, 3)
        datatable(df[, c("Description","Count","p.adjust")], options=list(pageLength=15))
    })
    output$download_results <- downloadHandler(
        filename = function() { paste0("hnsc_de_padj", input$padj, "_lfc", input$lfc, ".csv") },
        content = function(file) { write.csv(filtered_res(), file, row.names=TRUE) }
    )
}
shinyApp(ui = ui, server = server)
