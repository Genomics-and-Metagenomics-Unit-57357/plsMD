FROM ubuntu:22.04

LABEL maintainer="Maryam Lotfi <mariam.lotfi@57357.org> , Deena Jalal <deena.jalal@57357.org>"
LABEL description="Containerized plsMD pipeline for plasmid analysis"
LABEL version="v1.0"

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    PATH="/opt/conda/envs/plsMD/bin:/opt/conda/bin:$PATH" \
    INSTALL_DIR=/opt/plsMD \
    DATA_DIR=/opt/plsMD/data \
    BLAST_DB_DIR=/opt/plsMD/data/blastdb \
    SCRIPT_DIR=/opt/plsMD/scripts

RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends \
    ca-certificates \
    wget \
    git \
    build-essential \
    bzip2 \
    && rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    /bin/bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh

RUN conda config --set channel_priority strict && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda update -n base -y conda && \
    conda create -n plsMD -y python=3.9

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

RUN mkdir -p ${INSTALL_DIR} ${DATA_DIR} ${BLAST_DB_DIR} ${SCRIPT_DIR}

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

COPY Code/ ${SCRIPT_DIR}/
RUN chmod -R +x ${SCRIPT_DIR}/*.sh && \
    chmod +x ${SCRIPT_DIR}/*.py

RUN printf '#!/bin/bash\n\
case "$1" in\n\
  --preprocessing) shift; exec /opt/plsMD/scripts/plsMD_preprocessing.sh "$@" ;;\n\
  --processing) shift; exec python /opt/plsMD/scripts/plsMD_processing.py "$@" ;;\n\
  --annotation) shift; exec /opt/plsMD/scripts/plsMD_annotation.sh "$@" ;;\n\
  --phylogenetics) shift; exec /opt/plsMD/scripts/plsMD_phylogenetics.sh "$@" ;;\n\
  --db-path) echo "BLAST DB: /opt/plsMD/data/blastdb/plsdb" ;;\n\
  --version) echo "plsMD v1.0" ;;\n\
  --help|-h) cat <<EOF\n\
plsMD v1.0 - Plasmid reconstruction tool from short-read assemblies\n\
Usage: plsMD [command] [options]\n\
Commands:\n\
  --preprocessing    Processing & Annotation\n\
  --processing       Processing alignments and reconstruction\n\
  --annotation       Single sample modality\n\
  --phylogenetics    Batch modality\n\
  --db-path          database locations\n\
  --version          plsMD version\n\
EOF\n\
;;\n\
  *) echo "Error: Unknown command. Use --help for usage."; exit 1 ;;\n\
esac\n' > /usr/local/bin/plsMD && \
    chmod +x /usr/local/bin/plsMD

RUN /opt/conda/envs/plsMD/bin/python -c "import pandas, numpy, Bio; print('Dependencies verified')" && \
    /opt/conda/envs/plsMD/bin/abricate --version && \
    /opt/conda/envs/plsMD/bin/amrfinder --version

WORKDIR /data
VOLUME ["/data"]
ENTRYPOINT ["plsMD"]
CMD ["--help"]
