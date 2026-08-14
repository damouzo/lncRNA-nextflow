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

# ---- R 4.4.3 (pinned to match r_bioc container / Bioconductor 3.20) ----
# Required by CPAT's make_logitModel.py, which fits the glm model via Rscript.
# The CRAN jammy repo keeps historical builds, so 4.4.3 is pinned to stay in
# sync with the rest of the pipeline instead of pulling the latest R.
RUN apt-get update && apt-get install -y --no-install-recommends \
        software-properties-common \
        dirmngr \
    && wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
        | tee /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
    && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" \
    && printf 'Package: r-base r-base-core r-base-dev r-base-html r-mathlib r-recommended r-doc-html r-doc-pdf r-doc-info\nPin: version 4.4.3*\nPin-Priority: 1001\n' \
        > /etc/apt/preferences.d/cran-r \
    && apt-get update \
    && apt-get install -y --no-install-recommends r-base r-base-dev \
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
    && Rscript --version \
    && R -e 'glm(mpg ~ wt, data=mtcars)' >/dev/null 2>&1 \
    && echo "All coding tools verified successfully"

WORKDIR /data