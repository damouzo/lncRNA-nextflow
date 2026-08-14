# Discovery container
# Tools for Phase A (Discovery, Curation & Freezing) + Salmon index building

FROM ubuntu:22.04

LABEL org.opencontainers.image.description="Discovery tools for lncRNA-nextflow"
LABEL org.opencontainers.image.source="https://github.com/damouzo/lncRNA-nextflow"
LABEL org.opencontainers.image.licenses="MIT"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    curl \
    ca-certificates \
    build-essential \
    cmake \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libgsl-dev \
    libgsl27 \
    samtools \
    gffread \
    bedtools \
    subread \
    python3 \
    python3-pip \
    python3-numpy \
    && rm -rf /var/lib/apt/lists/*

# ---- StringTie (official binary; apt package is pinned to buggy 2.2.1) ----
# v3.0.0 includes "fixes and improvements for the -e expression quantification mode"
RUN curl -fsSL https://github.com/gpertea/stringtie/releases/download/v3.0.0/stringtie-3.0.0.Linux_x86_64.tar.gz \
    | tar xz -C /opt \
    && ln -s /opt/stringtie-3.0.0.Linux_x86_64/stringtie /usr/local/bin/stringtie

# ---- gffcompare (compiled from source) ----
RUN curl -fsSL https://github.com/gpertea/gffcompare/releases/download/v0.12.6/gffcompare-0.12.6.Linux_x86_64.tar.gz \
    | tar xz -C /opt \
    && ln -s /opt/gffcompare-0.12.6.Linux_x86_64/gffcompare /usr/local/bin/gffcompare

# ---- Salmon 1.10 ----
RUN curl -fsSL https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz \
    | tar xz -C /opt \
    && ln -s /opt/salmon-latest_linux_x86_64/bin/salmon /usr/local/bin/salmon

# ---- RSeQC ----
RUN pip3 install --no-cache-dir RSeQC

# ---- Verify installations ----
RUN samtools --version \
    && stringtie --version \
    && gffcompare --version 2>&1 | head -1 \
    && gffread --version 2>&1 | head -1 \
    && bedtools --version \
    && featureCounts -v 2>&1 | head -1 \
    && salmon --version \
    && infer_experiment.py --help 2>&1 | head -1 \
    && echo "All discovery tools verified successfully"

WORKDIR /data