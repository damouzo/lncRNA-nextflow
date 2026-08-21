# Quarto + R container for the per-contrast reporting engine (Phase B)
# Renders self-contained HTML reports (Quarto) from DESeq2 / cis / enrichment outputs.
#
# Packages:
#   - Quarto CLI 1.6 (HTML renderer)
#   - rmarkdown, knitr (Quarto execution engine)
#   - plotly, DT, kableExtra (self-contained interactive outputs)
#   - readr, dplyr, tibble, ggplot2 (data wrangling & plotting)

FROM rocker/verse:4.4.2

LABEL org.opencontainers.image.description="Quarto reporting for lncRNA-nextflow"
LABEL org.opencontainers.image.source="https://github.com/damouzo/lncRNA-nextflow"
LABEL org.opencontainers.image.licenses="MIT"

# Install Quarto CLI (pinned; overridable via build arg)
ARG QUARTO_VERSION=1.6.39

RUN apt-get update && apt-get install -y --no-install-recommends \
        wget ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN wget -q https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb \
    && dpkg -i quarto-*.deb \
    && rm quarto-*.deb \
    && quarto --version

# Reporting R packages
RUN R -e " \
    install.packages(c( \
        'rmarkdown', 'knitr', 'kableExtra', 'DT', 'plotly', \
        'readr', 'dplyr', 'tibble', 'ggplot2', 'glue' \
    ), repos='https://cloud.r-project.org', Ncpus=4) \
"

# Verify Quarto + R stack
RUN quarto check \
    && R -e "library(rmarkdown); library(knitr); library(DT); library(plotly); \
             library(ggplot2); library(dplyr); library(readr); sessionInfo()"

WORKDIR /data
