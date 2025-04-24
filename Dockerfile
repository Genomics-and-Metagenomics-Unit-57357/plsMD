FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    wget git build-essential ca-certificates curl \
    && rm -rf /var/lib/apt/lists/*

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh

ENV PATH="/opt/conda/bin:$PATH"

RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict

RUN conda create -n plsmd python=3.9 && \
    echo "source activate plsmd" > ~/.bashrc
ENV PATH /opt/conda/envs/plsmd/bin:$PATH

RUN conda install -n plsmd -y \
    abricate blast=2.10 seqtk mafft iqtree biopython wget

ENV INSTALL_DIR=/opt/plsMD
ENV DATA_DIR=$INSTALL_DIR/data
ENV BLAST_DB_DIR=$DATA_DIR/blastdb
RUN mkdir -p $INSTALL_DIR $DATA_DIR $BLAST_DB_DIR

RUN wget https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download_fasta -O $DATA_DIR/sequences.fasta && \
    makeblastdb -in $DATA_DIR/sequences.fasta -dbtype nucl -out $BLAST_DB_DIR/plsdb -title "PLSDB" -parse_seqids

COPY . $INSTALL_DIR/repo

RUN printf '#!/bin/bash\nSCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/repo" && pwd)"\nPLSDB_FASTA="$SCRIPT_DIR/../data/sequences.fasta"\nBLAST_DB="$SCRIPT_DIR/../data/blastdb/plsdb"\n\ncase "$1" in\n    --preprocessing) shift; "$SCRIPT_DIR/plsMD_preprocessing.sh" "$@" ;;\n    --processing) shift; python "$SCRIPT_DIR/plsMD_processing.py" "$@" ;;\n    --annotation) shift; "$SCRIPT_DIR/plsMD_annotation.sh" --db "$BLAST_DB" "$@" ;;\n    --phylogenetics) shift; "$SCRIPT_DIR/plsMD_phylogenetics.sh" "$@" ;;\n    --version|-v) echo "plsMD v1.9" ;;\n    --db-path) echo -e "PLSDB FASTA: $PLSDB_FASTA\\nBLAST DB: $BLAST_DB" ;;\n    --help|-h) echo "Usage: plsMD [command] [options]";;\n    *) echo "Invalid option. Use plsMD --help";;\nesac\n' > $INSTALL_DIR/plsMD && \
    chmod +x $INSTALL_DIR/plsMD

ENV PATH="$INSTALL_DIR:$PATH" \
    PLSDB_PATH="$DATA_DIR/sequences.fasta" \
    BLASTDB="$BLAST_DB_DIR"

WORKDIR /data
CMD ["bash"]
