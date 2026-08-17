# Coding-potential assessment container
# Tools: CPAT + CPC2 for lncRNA vs protein-coding classification

FROM ubuntu:22.04

LABEL org.opencontainers.image.description="Coding-potential tools for lncRNA-nextflow: CPAT + CPC2"
LABEL org.opencontainers.image.source="https://github.com/damouzo/lncRNA-nextflow"
LABEL org.opencontainers.image.licenses="MIT"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    ca-certificates \
    bzip2 \
    gffread \
    && rm -rf /var/lib/apt/lists/*

# ---- Micromamba: installing CPAT y CPC2 from Bioconda  ----
ENV MAMBA_ROOT_PREFIX=/opt/conda
RUN wget -qO- https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xj -C /usr/local/bin --strip-components=1 bin/micromamba

RUN micromamba create -y -n coding_potential -c conda-forge -c bioconda \
        cpat=3.0.4 \
        cpc2=1.0.1 \
    && micromamba clean -a -y

ENV PATH="/opt/conda/envs/coding_potential/bin:${PATH}"

# ---- Functional verification ----
RUN printf ">test_transcript\nATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGCTAGCTAGCTAGGGATCGATCGATCGA\n" > /tmp/test.fa \
    && cpat.py --version \
    && CPC2.py -i /tmp/test.fa -o /tmp/test_out \
    && test -s /tmp/test_out.txt \
    && gffread --version 2>&1 | head -1 \
    && Rscript --version \
    && Rscript -e 'glm(mpg ~ wt, data=mtcars)' >/dev/null \
    && echo "All coding-potential tools verified successfully"

WORKDIR /data