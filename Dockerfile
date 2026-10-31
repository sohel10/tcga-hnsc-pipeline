FROM rocker/r-ver:4.6.1
LABEL maintainer="Sohel Ahmed <sohelcu06@gmail.com>"
LABEL description="TCGA RNA-seq DE analysis pipeline with Bioconductor"
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev \
    libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev \
    libtiff5-dev libjpeg-dev libcairo2-dev procps \
    && rm -rf /var/lib/apt/lists/*
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" && \
    R -e "BiocManager::install(c('DESeq2','edgeR','GenomicRanges','TCGAbiolinks','clusterProfiler','org.Hs.eg.db','EnhancedVolcano','SummarizedExperiment'), ask=FALSE, update=TRUE)"
RUN R -e "install.packages(c('pheatmap','ggplot2','plotly','shiny','DT','rmarkdown'), repos='https://cloud.r-project.org')"
COPY main.nf /pipeline/main.nf
COPY nextflow.config /pipeline/nextflow.config
COPY shiny/ /pipeline/shiny/
WORKDIR /pipeline
CMD ["Rscript", "--version"]
