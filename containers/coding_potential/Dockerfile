# Coding-potential assessment container
# Tools: CPAT + CPC2 for lncRNA vs protein-coding classification

FROM ubuntu:22.04

LABEL org.opencontainers.image.description="Coding-potential tools for lncRNA-nextflow: CPAT + CPC2"
LABEL org.opencontainers.image.source="https://github.com/damouzo/lncRNA-nextflow"
LABEL org.opencontainers.image.licenses="MIT"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-dev \
    gcc \
    g++ \
    make \
    zlib1g-dev \
    python3-numpy \
    python3-scipy \
    python3-sklearn \
    python3-biopython \
    gffread \
    git \
    wget \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# ---- CPAT 3.0.4 ----
RUN pip3 install --no-cache-dir CPAT==3.0.4

# ---- CPC2 (cloned from official GitHub repository) ----
RUN git clone https://github.com/gao-lab/CPC2_standalone.git /opt/CPC2 \
    && chmod +x /opt/CPC2/bin/CPC2.py

ENV PATH="/opt/CPC2/bin:${PATH}"

# ---- Verify ----
RUN python3 -c "import numpy; print('numpy', numpy.__version__)" \
    && python3 -c "import scipy; print('scipy', scipy.__version__)" \
    && python3 -c "import Bio; print('Biopython', Bio.__version__)" \
    && cpat.py --version \
    && CPC2.py --help 2>&1 | head -2 \
    && gffread --version 2>&1 | head -1 \
    && echo "All coding tools verified successfully"

WORKDIR /data