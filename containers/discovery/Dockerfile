# Discovery container
# Tools for Phase A (Discovery, Curation & Freezing) + Salmon index building
#
# Includes:
#   - samtools, RSeQC (infer_experiment.py), StringTie2, gffcompare,
#     gffread, bedtools, featureCounts, Salmon

FROM ubuntu:22.04

LABEL org.opencontainers.image.description="Discovery tools for lncRNA-nextflow: samtools, RSeQC, StringTie2, gffcompare, gffread, bedtools, featureCounts, Salmon"
LABEL org.opencontainers.image.source="https://github.com/damouzo/lncRNA-nextflow"
LABEL org.opencontainers.image.licenses="MIT"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    curl \
    ca-certificates \
    libcurl4-openssl-dev \
    libssl-dev \
    libgsl-dev \
    libgsl27 \
    && rm -rf /var/lib/apt/lists/*

# ---- Salmon 1.10 ----
RUN curl -fsSL https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz \
    | tar xz -C /opt \
    && ln -s /opt/salmon-latest_linux_x86_64/bin/salmon /usr/local/bin/salmon

# ---- samtools ----
RUN apt-get update && apt-get install -y --no-install-recommends \
    samtools \
    && rm -rf /var/lib/apt/lists/*

# ---- StringTie2 ----
RUN apt-get update && apt-get install -y --no-install-recommends \
    stringtie \
    && rm -rf /var/lib/apt/lists/*

# ---- gffcompare ----
RUN apt-get update && apt-get install -y --no-install-recommends \
    gffcompare \
    && rm -rf /var/lib/apt/lists/*

# ---- gffread ----
RUN apt-get update && apt-get install -y --no-install-recommends \
    gffread \
    && rm -rf /var/lib/apt/lists/*

# ---- bedtools ----
RUN apt-get update && apt-get install -y --no-install-recommends \
    bedtools \
    && rm -rf /var/lib/apt/lists/*

# ---- featureCounts (subread) ----
RUN apt-get update && apt-get install -y --no-install-recommends \
    subread \
    && rm -rf /var/lib/apt/lists/*

# ---- RSeQC (infer_experiment.py) ----
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-numpy \
    && rm -rf /var/lib/apt/lists/* \
    && pip3 install --no-cache-dir RSeQC

# ---- Verify installations ----
RUN samtools --version \
    && stringtie --version \
    && gffcompare --version 2>&1 | head -1 \
    && gffread --version 2>&1 | head -1 \
    && bedtools --version \
    && featureCounts -v 2>&1 | head -1 \
    && salmon --version \
    && infer_experiment.py --help 2>&1 | head -1 \
    && echo "All tools verified"

WORKDIR /data
