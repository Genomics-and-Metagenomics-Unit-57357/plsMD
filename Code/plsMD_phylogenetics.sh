#!/bin/bash

# Main CLI handler
usage() {
    echo "Usage: $0 --genes_dir <genes_directory> --min_length <min_length> --max_length <max_length> --threads <threads>"
    echo "Flags:"
    echo "  --genes_dir      Root directory containing plasmid folders"
    echo "  --min_length     Minimum sequence length for common sequence search (default: 10)"
    echo "  --max_length     Maximum sequence length for common sequence search (default: 30)"
    echo "  --threads        Number of threads for MAFFT and IQ-TREE (default: 1)"
    echo "  --help           Show this help"
    exit 1
}

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --genes_dir)
            GENES_DIR="$2"
            shift 2
            ;;
        --min_length)
            MIN_LENGTH="$2"
            shift 2
            ;;
        --max_length)
            MAX_LENGTH="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Validate arguments
if [[ -z "$GENES_DIR" || ! -d "$GENES_DIR" ]]; then
    echo "Error: --genes_dir is required and must be a valid directory."
    usage
fi
MIN_LENGTH=${MIN_LENGTH:-10}
MAX_LENGTH=${MAX_LENGTH:-30}
THREADS=${THREADS:-1}

# Python part: Rotate sequences
echo "Running Python script to rotate sequences..."

python3 - <<EOF_PYTHON "$GENES_DIR" "$MIN_LENGTH" "$MAX_LENGTH"
import os
import sys
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord

genes_dir = sys.argv[1]
min_length = int(sys.argv[2])
max_length = int(sys.argv[3])

def find_common_sequence(sequences, gene_directory, min_length=10, max_length=30):
    """Find the longest common sequence among given sequences."""
    for length in range(max_length, min_length - 1, -1):
        common_substrings = None

        for seq, rev_comp in sequences:
            substrings = set()
            for i in range(len(seq) - length + 1):
                substrings.add(seq[i:i+length])
            for i in range(len(rev_comp) - length + 1):
                substrings.add(rev_comp[i:i+length])

            if common_substrings is None:
                common_substrings = substrings
            else:
                common_substrings.intersection_update(substrings)
                if not common_substrings:
                    break

        if not common_substrings:
            continue

        valid_substrings = set()
        for substring in common_substrings:
            valid = True
            for s, rc in sequences:
                total = s.count(substring) + rc.count(substring)
                if total != 1:
                    valid = False
                    break
            if valid:
                valid_substrings.add(substring)

        if valid_substrings:
            longest_common_sequence = next(iter(valid_substrings))
            output_dir = gene_directory
            os.makedirs(output_dir, exist_ok=True)
            output_file = os.path.join(output_dir, f"{os.path.basename(gene_directory)}_common_seq.txt")
            with open(output_file, 'w') as f:
                f.write(longest_common_sequence)
            print(f"Longest common sequence found and saved to: {output_file}")
            return longest_common_sequence

    print("No common sequence found in the specified length range.")
    return None

def rotate_sequences(input_file, output_file, common_seq):
    """Rotate sequences to start with the common sequence."""
    rotated_records = []
    for record in SeqIO.parse(input_file, "fasta"):
        original = str(record.seq)
        rev_comp = str(record.seq.reverse_complement())

        if common_seq in original:
            idx = original.index(common_seq)
            rotated = original[idx:] + original[:idx]
        elif common_seq in rev_comp:
            idx = rev_comp.index(common_seq)
            rotated = rev_comp[idx:] + rev_comp[:idx]
        else:
            print(f"Warning: Common sequence not found in {record.id}. Skipping rotation.")
            rotated = original

        rotated_record = SeqRecord(
            Seq.Seq(rotated),
            id=record.id,
            description=""
        )
        rotated_records.append(rotated_record)

    with open(output_file, 'w') as f:
        SeqIO.write(rotated_records, f, "fasta")

def concatenate_fasta_files_in_plasmid_directory(plasmid_directory, min_length, max_length):
    """Process plasmid directory by concatenating all FASTA files directly in the directory."""
    plasmid_name = os.path.basename(plasmid_directory)
    output_fasta = os.path.join(plasmid_directory, f"{plasmid_name}_concatenated.fasta")
    
    # Find all FASTA files directly in the plasmid directory
    fasta_files = [f for f in os.listdir(plasmid_directory) 
                   if (f.endswith(('.fasta', '.fa', '.fna')) and 
                       not f.endswith('_extracted.fasta') and
                       not f.endswith('_rotated.fasta') and
                       not f.endswith('_concatenated.fasta') and
                       not f.endswith('_common_seq.txt'))]
    
    if not fasta_files:
        print(f"No FASTA files found in {plasmid_directory}")
        return
    
    print(f"Found {len(fasta_files)} FASTA files in {plasmid_directory}")
    
    # Concatenate all FASTA files
    with open(output_fasta, 'w') as out_fh:
        for fasta_file in fasta_files:
            file_path = os.path.join(plasmid_directory, fasta_file)
            try:
                # Read and write each sequence individually to preserve headers
                for record in SeqIO.parse(file_path, "fasta"):
                    # Use original filename + sequence ID to avoid duplicate IDs
                    record.id = f"{fasta_file}_{record.id}"
                    record.description = ""
                    SeqIO.write(record, out_fh, "fasta")
            except Exception as e:
                print(f"Error reading {file_path}: {e}")

    # Check if we actually wrote any sequences
    if os.path.getsize(output_fasta) == 0:
        print(f"No valid sequences found in {plasmid_directory}")
        os.remove(output_fasta)
        return

    # Find common sequence
    try:
        sequences = [(str(r.seq), str(r.seq.reverse_complement()))
                     for r in SeqIO.parse(output_fasta, "fasta")]
        
        if not sequences:
            print(f"No sequences to process in {plasmid_directory}")
            return
            
        common_seq = find_common_sequence(sequences, plasmid_directory, min_length, max_length)

        if common_seq:
            # Rotate sequences
            rotated_fasta = os.path.join(plasmid_directory, f"{plasmid_name}_rotated.fasta")
            rotate_sequences(output_fasta, rotated_fasta, common_seq)
            print(f"Created rotated sequences at: {rotated_fasta}")
        else:
            print(f"No common sequence found for {plasmid_directory}. Skipping rotation.")
    except Exception as e:
        print(f"Error processing {plasmid_directory}: {e}")

def process_all_plasmid_directories(root_dir, min_len, max_len):
    """Process only the direct subdirectories of root directory (plasmid directories)."""
    for item in os.listdir(root_dir):
        plasmid_dir = os.path.join(root_dir, item)
        if os.path.isdir(plasmid_dir):
            print(f"Processing plasmid directory: {plasmid_dir}")
            concatenate_fasta_files_in_plasmid_directory(plasmid_dir, min_len, max_len)
        else:
            print(f"Skipping non-directory: {item}")

# Run the script
process_all_plasmid_directories(genes_dir, min_length, max_length)
EOF_PYTHON

# Bash part: Run MAFFT and IQ-TREE
# Create phylogenetic_tree directory
PHYLO_DIR="${GENES_DIR}/phylogenetic_tree"
mkdir -p "$PHYLO_DIR"
echo "Created phylogenetic_tree directory at: $PHYLO_DIR"
echo "Running MAFFT and IQ-TREE..."

for plasmid_dir in "$GENES_DIR"/*; do
    if [[ -d "$plasmid_dir" ]]; then
        plasmid_name=$(basename "$plasmid_dir")
        rotated_fasta="${plasmid_dir}/${plasmid_name}_rotated.fasta"
        aligned_fasta="${PHYLO_DIR}/${plasmid_name}_aligned.fasta"
        tree_output="${PHYLO_DIR}/${plasmid_name}_tree"

        if [[ -f "$rotated_fasta" ]]; then
            # Run MAFFT
            echo "Running MAFFT on $rotated_fasta..."
            mafft --auto --thread "$THREADS" "$rotated_fasta" > "$aligned_fasta"

            # Run IQ-TREE
            echo "Running IQ-TREE on $aligned_fasta..."
            iqtree -s "$aligned_fasta" -m TEST -bb 1000 -T "$THREADS" -pre "$tree_output"
        else
            echo "Rotated FASTA file not found for plasmid: $plasmid_name"
        fi
    fi
done

echo "Phylogenetic analysis completed. Results saved in: $PHYLO_DIR"
