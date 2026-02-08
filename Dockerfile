FROM ubuntu:22.04

LABEL maintainer="Maryam Lotfi <mariam.lotfi@57357.org>, Deena Jalal <deena.jalal@57357.org>"
LABEL description="Containerized plsMD pipeline for plasmid analysis"
LABEL version="v1.0"

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    PATH="/opt/conda/envs/plsMD/bin:/opt/conda/bin:$PATH" \
    INSTALL_DIR=/opt/plsMD \
    DATA_DIR=/opt/plsMD/data \
    BLAST_DB_DIR=/opt/plsMD/data/blastdb \
    SCRIPT_DIR=/opt/plsMD/scripts \
    ABRICATE_DB_DIR=/opt/conda/envs/plsMD/db

# Install system dependencies
RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends \
    ca-certificates \
    wget \
    git \
    build-essential \
    bzip2 \
    tar \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    /bin/bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh

# Configure conda and create environment
RUN conda config --set channel_priority strict && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda update -n base -y conda && \
    conda create -n plsMD -y python=3.9

# Install bioinformatics tools
RUN conda install -n plsMD -y \
    pandas \
    numpy \
    biopython \
    abricate \
    blast=2.13.0 \
    seqtk \
    mafft \
    iqtree \
    ncbi-amrfinderplus \
    && conda clean -afy

# Create necessary directories
RUN mkdir -p ${INSTALL_DIR} ${DATA_DIR} ${BLAST_DB_DIR} ${SCRIPT_DIR} ${ABRICATE_DB_DIR}

# Download and setup custom ABRicate database
ARG GITHUB_REPO=https://github.com/Genomics-and-Metagenomics-Unit-57357/plsMD
RUN wget ${GITHUB_REPO}/raw/main/rep.mob.typer.tar.gz -O /tmp/rep.mob.typer.tar.gz && \
    tar -xzf /tmp/rep.mob.typer.tar.gz -C ${ABRICATE_DB_DIR} && \
    rm /tmp/rep.mob.typer.tar.gz && \
    /opt/conda/envs/plsMD/bin/abricate --setupdb

# Optional: Download PLSDB database
ARG DOWNLOAD_DB=false
ARG PLSDB_URL=https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download_fasta
RUN if [ "${DOWNLOAD_DB}" = "true" ]; then \
    wget --progress=dot:giga ${PLSDB_URL} -O ${DATA_DIR}/sequences.fasta && \
    makeblastdb -in ${DATA_DIR}/sequences.fasta \
        -dbtype nucl \
        -out ${BLAST_DB_DIR}/plsdb \
        -title "PLSDB" \
        -parse_seqids; \
fi

# Copy scripts
COPY Code/ ${SCRIPT_DIR}/
RUN chmod -R +x ${SCRIPT_DIR}/*.sh && \
    chmod +x ${SCRIPT_DIR}/*.py

# Create plsMD wrapper that properly handles paths
RUN printf '#!/bin/bash\n\
set -e\n\
\n\
# Activate conda environment\n\
source /opt/conda/etc/profile.d/conda.sh\n\
conda activate plsMD\n\
\n\
case "$1" in\n\
  --preprocessing)\n\
    shift\n\
    exec /opt/plsMD/scripts/plsMD_preprocessing.sh "$@"\n\
    ;;\n\
  --processing)\n\
    shift\n\
    exec python /opt/plsMD/scripts/plsMD_processing.py "$@"\n\
    ;;\n\
  --annotation)\n\
    shift\n\
    exec /opt/plsMD/scripts/plsMD_annotation.sh "$@"\n\
    ;;\n\
  --phylogenetics)\n\
    shift\n\
    exec /opt/plsMD/scripts/plsMD_phylogenetics.sh "$@"\n\
    ;;\n\
  --db-path)\n\
    echo "BLAST DB: /opt/plsMD/data/blastdb/plsdb"\n\
    echo "ABRicate DB: /opt/conda/envs/plsMD/db"\n\
    /opt/conda/envs/plsMD/bin/abricate --list\n\
    ;;\n\
  --version)\n\
    echo "plsMD v1.0"\n\
    ;;\n\
  --help|-h)\n\
    cat <<EOF\n\
plsMD v1.0 - Plasmid reconstruction tool from short-read assemblies\n\
\n\
Usage: plsMD [command] [options]\n\
\n\
Commands:\n\
  --preprocessing    Processing & Annotation\n\
  --processing       Processing alignments and reconstruction\n\
  --annotation       Single sample modality\n\
  --phylogenetics    Batch modality\n\
  --db-path          Show database locations\n\
  --version          Show plsMD version\n\
  --help             Show this help message\n\
\n\
Examples:\n\
  plsMD --preprocessing --dir /data/input\n\
  plsMD --annotation --input /data/sample.fasta --output /data/results\n\
\n\
Note: All paths should be relative to /data or use absolute paths from mounted volumes\n\
EOF\n\
    ;;\n\
  *)\n\
    echo "Error: Unknown command '\''$1'\''. Use --help for usage." >&2\n\
    exit 1\n\
    ;;\n\
esac\n' > /usr/local/bin/plsMD && \
    chmod +x /usr/local/bin/plsMD

# Verify installation
RUN /opt/conda/envs/plsMD/bin/python -c "import pandas, numpy, Bio; print('Dependencies verified')" && \
    /opt/conda/envs/plsMD/bin/abricate --version && \
    /opt/conda/envs/plsMD/bin/amrfinder --version && \
    /opt/conda/envs/plsMD/bin/abricate --list

WORKDIR /data
VOLUME ["/data"]

ENTRYPOINT ["plsMD"]
CMD ["--help"]
