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

# ---- pyBigWig (bigWig chrom/length checks + per-interval conservation scoring) ----
RUN pip3 install --no-cache-dir pyBigWig

# ---- UCSC binaries: liftOver (synteny) ----
# bigWig handling is done via pyBigWig above, so only liftOver is fetched here.
# Caveat: UCSC prebuilt binaries (hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64)
# are linked against an older glibc. On glibc >= 2.34 (ubuntu 22.04) liftOver
# may need a rebuild from source; keep the base pinned if it starts failing.
RUN wget -q https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver \
    && chmod +x liftOver \
    && mv liftOver /usr/local/bin/liftOver \
    && liftOver 2>&1 | head -1 || true

# Shared chr-prefix normalization helper used by CONSERVATION_SCORE,
# SYNTENY_CHECK, and CATALOG_OVERLAP — all three call it from PATH.
COPY bin/normalize_chrom.py /usr/local/bin/normalize_chrom.py

# ---- Verify installations ----
RUN samtools --version \
    && stringtie --version \
    && gffcompare --version 2>&1 | head -1 \
    && gffread --version 2>&1 | head -1 \
    && bedtools --version \
    && featureCounts -v 2>&1 | head -1 \
    && salmon --version \
    && infer_experiment.py --help 2>&1 | head -1 \
    && python3 -c "import pyBigWig; print('pyBigWig', pyBigWig.__version__)" \
    && liftOver 2>&1 | head -1 || true \
    && echo "All discovery tools verified successfully"

WORKDIR /data