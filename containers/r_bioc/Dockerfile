# R/Bioconductor container for quantification & statistics (Phase B)
# Based on Bioconductor 3.20 (R 4.4.x)
#
# Packages:
#   - tximport, DESeq2, apeglm, ashr (core DE stack)
#   - GenomicRanges, rtracklayer, GenomicFeatures (cis-association)
#   - clusterProfiler, org.Hs.eg.db, enrichplot, fgsea (functional enrichment)
#   - ggplot2, dplyr, readr, tibble (data wrangling & plotting)

FROM bioconductor/bioconductor_docker:RELEASE_3_20

LABEL org.opencontainers.image.description="R/Bioconductor quantification & statistics for lncRNA-nextflow"
LABEL org.opencontainers.image.source="https://github.com/damouzo/lncRNA-nextflow"
LABEL org.opencontainers.image.licenses="MIT"

RUN R -e " \
    BiocManager::install(c( \
        'tximport', \
        'DESeq2', \
        'apeglm', \
        'ashr', \
        'GenomicRanges', \
        'rtracklayer', \
        'GenomicFeatures', \
        'clusterProfiler', \
        'org.Hs.eg.db', \
        'org.Mm.eg.db', \
        'enrichplot', \
        'fgsea' \
    ), ask=FALSE, update=FALSE) \
"

RUN R -e " \
    install.packages(c( \
        'ggplot2', 'dplyr', 'readr', 'tibble', 'ggrepel' \
    ), repos='https://cloud.r-project.org', Ncpus=4) \
"

# Verify core packages load
RUN R -e " \
    library(tximport); library(DESeq2); library(apeglm); \
    library(GenomicRanges); library(rtracklayer); library(clusterProfiler); \
    library(ggplot2); library(dplyr); library(readr); \
    sessionInfo() \
"

WORKDIR /data
