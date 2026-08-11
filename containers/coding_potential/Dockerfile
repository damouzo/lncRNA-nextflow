# Coding-potential assessment container
# Tools: CPAT + CPC2 for lncRNA vs protein-coding classification
#
# Dependencies:
#   - CPAT: Python3 + numpy + scipy + scikit-learn
#   - CPC2: Python3 + BioPython
#   - Both require a FASTA/GTF extraction tool (gffread) to prepare input sequences

FROM ubuntu:22.04

LABEL org.opencontainers.image.description="Coding-potential tools for lncRNA-nextflow: CPAT + CPC2"
LABEL org.opencontainers.image.source="https://github.com/damouzo/lncRNA-nextflow"
LABEL org.opencontainers.image.licenses="MIT"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-numpy \
    python3-scipy \
    python3-sklearn \
    python3-biopython \
    git \
    wget \
    && rm -rf /var/lib/apt/lists/*

# ---- CPAT 3.0.4 ----
RUN wget -qO- https://github.com/sun-a/cpat/archive/refs/tags/v3.0.4.tar.gz \
    | tar xz -C /opt \
    && mv /opt/cpat-3.0.4 /opt/cpat \
    && chmod +x /opt/cpat/bin/*.py

ENV PATH="/opt/cpat/bin:${PATH}"

# ---- CPC2 ----
RUN wget -qO- http://cpc2.gao-lab.org/data/CPC2-beta.tar.gz \
    | tar xz -C /opt \
    && chmod +x /opt/CPC2-beta/bin/CPC2.py

ENV PATH="/opt/CPC2-beta/bin:${PATH}"

# ---- Verify ----
RUN python3 -c "import numpy; print('numpy', numpy.__version__)" \
    && python3 -c "import scipy; print('scipy', scipy.__version__)" \
    && python3 -c "import Bio; print('Biopython', Bio.__version__)" \
    && cpat.py --help 2>&1 | head -2 \
    && CPC2.py --help 2>&1 | head -2 \
    && echo "All tools verified"

WORKDIR /data
