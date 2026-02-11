#!/bin/bash

IS_DB_DEFAULT="${IS_DB_DIR:-/opt/plsMD/data/blastdb/IS}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --dir)
      ROOT_DIR="$2"
      shift 2
      ;;

    --IS_db)
      IS_DB_OVERRIDE="$2"
      shift 2
      ;;
    *)
      echo "Unknown parameter: $1"
      exit 1
      ;;
  esac
done

if [[ -z "$THREADS" || -z "$ROOT_DIR" ]]; then
  echo "Usage: $0 --threads <threads> --dir <root_directory> [--IS_db <IS_database_path>]"
  exit 1
fi


IS_DB="${IS_DB_OVERRIDE:-$IS_DB_DEFAULT}"


if [[ ! -d "$IS_DB" ]]; then
  echo "Error: IS database directory not found at: $IS_DB"
  echo "Ensure the database was downloaded during Docker build, or pass --IS_db <path>"
  exit 1
fi


IS_DB_PREFIX=$(find "$IS_DB" -name "*.nhr" | head -1 | sed 's/\.nhr$//')
if [[ -z "$IS_DB_PREFIX" ]]; then
  echo "Error: No BLAST index (.nhr) files found in IS database directory: $IS_DB"
  exit 1
fi

echo "Using IS database: $IS_DB_PREFIX"

process_directory() {
  local dir="$1"
  local output_dir="$2"
  
  mkdir -p "$output_dir/annotation/AMR"
  mkdir -p "$output_dir/annotation/VF"
  mkdir -p "$output_dir/annotation/PL"
  mkdir -p "$output_dir/annotation/IS"

  # AMR Genes
  for f in "$dir"/*.fasta; do
    if [[ -f "$f" ]]; then 
        base_name="${f##*/}" 
        base_name="${base_name%.fasta}"
        amrfinder --nucleotide "$f" --threads "$THREADS" --output "$output_dir/annotation/AMR/${base_name}_AMR.txt"
    fi
  done

  # VF
  for f in "$dir"/*.fasta; do
    if [[ -f "$f" ]]; then
        base_name="${f##*/}"
        base_name="${base_name%.fasta}"
        abricate --db vfdb --threads "$THREADS" "$f" > "$output_dir/annotation/VF/${base_name}_VF.txt"
    fi
  done

  # Plasmids
  for f in "$dir"/*.fasta; do
    if [[ -f "$f" ]]; then
        base_name="${f##*/}"
        base_name="${base_name%.fasta}"
        abricate --db plasmidfinder --threads "$THREADS" "$f" > "$output_dir/annotation/PL/${base_name}_PL.txt"
    fi
  done

  # Insertion Sequences
  for f in "$dir"/*.fasta; do
    if [[ -f "$f" ]]; then
        base_name="${f##*/}"
        base_name="${base_name%.fasta}"
        blastn -db "$IS_DB_PREFIX" -query "$f" -num_threads "$THREADS" -outfmt 6 > "$output_dir/annotation/IS/${base_name}_IS.txt"
    fi
  done
}

# Process plasmid_files directory
if [[ -d "$ROOT_DIR/plasmid_files" ]]; then
  process_directory "$ROOT_DIR/plasmid_files" "$ROOT_DIR/plasmid_files"
fi

# Process nonplasmid_files directory
if [[ -d "$ROOT_DIR/nonplasmid_files" ]]; then
  process_directory "$ROOT_DIR/nonplasmid_files" "$ROOT_DIR/nonplasmid_files"
fi

echo "Annotation process completed."
