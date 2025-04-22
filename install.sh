#!/bin/bash

# plsMD Installation Script v1.0

VERSION="v1.0"
REPO_URL="https://github.com/Genomics-and-Metagenomics-Unit-57357/plsMD.git"
PLSDB_URL="https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download_fasta"  # Example PLSDB URL
INSTALL_DIR="/opt/plsMD"
DATA_DIR="$INSTALL_DIR/data"
BLAST_DB_DIR="$DATA_DIR/blastdb"

echo "=== plsMD Installation ($VERSION) ==="
echo "Including PLSDB download and BLAST database setup"


if [ "$(id -u)" -ne 0 ]; then
    echo "Please run as root or with sudo for system-wide installation"
    exit 1
fi


if command -v plsMD &> /dev/null; then
    echo "plsMD is already installed. Checking for updates..."
    plsMD --version
    echo "To reinstall, first remove the existing installation."
    exit 0
fi


if ! command -v conda &> /dev/null; then
    echo "Error: Conda is required but not installed."
    echo "Please install Miniconda first: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi


echo "Creating installation directories..."
mkdir -p "$INSTALL_DIR" "$DATA_DIR" "$BLAST_DB_DIR"
cd "$INSTALL_DIR"


echo "Setting up conda environment..."
conda create -n plsMD python=3.7 -y
conda activate plsMD

echo "Installing dependencies..."
conda install -c bioconda abricate blast=2.10 seqtk mafft iqtree biopython wget -y


echo "Downloading plsMD..."
git clone $REPO_URL "$INSTALL_DIR/repo"


echo "Downloading PLSDB sequences..."
cd "$DATA_DIR"
wget $PLSDB_URL -O plsdb.fna.gz
gunzip plsdb.fna.gz

echo "Building BLAST database..."
makeblastdb -in plsdb.fna -dbtype nucl -out "$BLAST_DB_DIR/plsdb" -title "PLSDB"


echo "Setting up environment variables..."
cat << EOF >> ~/.bashrc
# plsMD configuration
export PLSDB_PATH="$DATA_DIR/plsdb.fna"
export BLASTDB="$BLAST_DB_DIR"
conda activate plsMD
EOF


echo "Configuring scripts..."
cd "$INSTALL_DIR/repo"


sed -i "s|blastdb/plsdb|$BLAST_DB_DIR/plsdb|g" *.sh
sed -i "s|path/to/plsdb.fna|$DATA_DIR/plsdb.fna|g" *.sh *.py

chmod +x *.sh


cat << EOF > /usr/local/bin/plsMD
#!/bin/bash
SCRIPT_DIR="$INSTALL_DIR/repo"
PLSDB_PATH="$DATA_DIR/plsdb.fna"
BLAST_DB="$BLAST_DB_DIR/plsdb"

# Activate conda environment
source \$(conda info --base)/etc/profile.d/conda.sh
conda activate plsMD

case "\$1" in
    --preprocessing) shift; "\$SCRIPT_DIR/plsMD_preprocessing.sh" "\$@" ;;
    --processing) shift; python "\$SCRIPT_DIR/plsMD_processing.py" "\$@" ;;
    --annotation) shift; "\$SCRIPT_DIR/plsMD_annotation.sh" --db "\$BLAST_DB" "\$@" ;;
    --phylogenetics) shift; "\$SCRIPT_DIR/plsMD_phylogenetics.sh" "\$@" ;;
    --version|-v) echo "plsMD $VERSION" ;;
    --db-path) echo -e "PLSDB sequences: \$PLSDB_PATH\\nBLAST database: \$BLAST_DB" ;;
    *) echo "Usage: plsMD [--preprocessing|--processing|--annotation|--phylogenetics|--version|--db-path]" ;;
esac
EOF

chmod +x /usr/local/bin/plsMD


cat << EOF > "$INSTALL_DIR/uninstall_plsMD.sh"
#!/bin/bash
echo "=== plsMD Uninstallation ==="
echo "Removing installation..."
rm -rf "$INSTALL_DIR"
rm -f /usr/local/bin/plsMD
conda env remove -n plsMD
sed -i '/# plsMD configuration/,/conda activate plsMD/d' ~/.bashrc
echo "plsMD has been completely removed"
EOF
chmod +x "$INSTALL_DIR/uninstall_plsMD.sh"

echo ""
echo "Installation complete!"
echo "plsMD $VERSION with PLSDB database is now ready to use."
echo ""
echo "Database locations:"
echo "  - PLSDB sequences: $DATA_DIR/plsdb.fna"
echo "  - BLAST database: $BLAST_DB_DIR/plsdb"
echo ""
echo "Usage:"
echo "    conda activate plsMD"
echo "    plsMD --annotation [options]  # Uses built BLAST database"
echo "    plsMD --db-path              # Show database locations"
echo ""
echo "To uninstall:"
echo "    $INSTALL_DIR/uninstall_plsMD.sh"
