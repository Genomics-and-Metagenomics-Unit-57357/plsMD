#!/bin/bash
IS_DB_DEFAULT="${IS_DB_DIR:-/opt/plsMD/data/blastdb/IS}"

usage() {
    echo "Usage: $0 --threads <threads> --dir <input_directory> --output <output_directory> [--IS_db <IS_database_path>]"
    echo "  --threads     Number of threads"
    echo "  --dir         Root directory containing plasmid_files/ and nonplasmid_files/ subdirectories"
    echo "  --output      Output directory for annotation results"
    echo "  --IS_db       Path to IS database directory (optional if built into Docker image)"
    exit 1
}

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
    --output)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --IS_db)
      IS_DB_OVERRIDE="$2"
      shift 2
      ;;
    --help|-h)
      usage
      ;;
    *)
      echo "Unknown parameter: $1"
      usage
      ;;
  esac
done

if [[ -z "$THREADS" || -z "$ROOT_DIR" || -z "$OUTPUT_DIR" ]]; then
  echo "Usage: $0 --threads <threads> --dir <root_directory> --output <output_directory> [--IS_db <IS_database_path>]"
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

ROOT_DIR="$(realpath "$ROOT_DIR")"
OUTPUT_DIR="$(realpath "$OUTPUT_DIR")"
mkdir -p "$OUTPUT_DIR" || { echo "Error: Cannot create output directory $OUTPUT_DIR"; exit 1; }

echo "Using IS database: $IS_DB_PREFIX"
echo "Input directory:   $ROOT_DIR"
echo "Output directory:  $OUTPUT_DIR"

process_directory() {
  local input_dir="$1"
  local output_dir="$2"

  if [[ ! -d "$input_dir" ]]; then
    echo "Skipping non-existent directory: $input_dir"
    return
  fi

  mkdir -p "$output_dir/annotation/AMR"
  mkdir -p "$output_dir/annotation/VF"
  mkdir -p "$output_dir/annotation/PL"
  mkdir -p "$output_dir/annotation/IS"

  # AMR Genes
  for f in "$input_dir"/*.fasta; do
    [[ -f "$f" ]] || continue
    base_name="${f##*/}"
    base_name="${base_name%.fasta}"
    amrfinder --nucleotide "$f" --threads "$THREADS" --output "$output_dir/annotation/AMR/${base_name}_AMR.txt"
  done

  # VF
  for f in "$input_dir"/*.fasta; do
    [[ -f "$f" ]] || continue
    base_name="${f##*/}"
    base_name="${base_name%.fasta}"
    abricate --db vfdb --threads "$THREADS" "$f" > "$output_dir/annotation/VF/${base_name}_VF.txt"
  done

  # Plasmids
  for f in "$input_dir"/*.fasta; do
    [[ -f "$f" ]] || continue
    base_name="${f##*/}"
    base_name="${base_name%.fasta}"
    abricate --db plasmidfinder --threads "$THREADS" "$f" > "$output_dir/annotation/PL/${base_name}_PL.txt"
  done

  # Insertion Sequences
  for f in "$input_dir"/*.fasta; do
    [[ -f "$f" ]] || continue
    base_name="${f##*/}"
    base_name="${base_name%.fasta}"
    blastn -db "$IS_DB_PREFIX" -query "$f" -num_threads "$THREADS" -outfmt 6 > "$output_dir/annotation/IS/${base_name}_IS.txt"
  done
}

if [[ -d "$ROOT_DIR/plasmid_files" ]]; then
  echo "Annotating plasmid files..."
  process_directory "$ROOT_DIR/plasmid_files" "$OUTPUT_DIR/plasmid_files"
else
  echo "Warning: $ROOT_DIR/plasmid_files not found, skipping."
fi

if [[ -d "$ROOT_DIR/nonplasmid_files" ]]; then
  echo "Annotating non-plasmid files..."
  process_directory "$ROOT_DIR/nonplasmid_files" "$OUTPUT_DIR/nonplasmid_files"
else
  echo "Warning: $ROOT_DIR/nonplasmid_files not found, skipping."
fi

echo "Annotation process completed. Results saved in: $OUTPUT_DIR"
