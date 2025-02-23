#!/bin/bash

usage() {
    echo "Usage: $0 --dir <input_directory> --threads <num_threads> --db <plsdb_path>"
    echo "  --dir: Specify the input directory containing fasta files."
    echo "  --threads: Number of threads for abricate and blastn."
    echo "  --db: Path to the PLSDB database for blastn."
    echo "  -h, --help: Display this help message."
    exit 1
}

# Check for dependencies
check_dependency() {
    if ! command -v "$1" &> /dev/null; then
        echo "Error: $1 is not installed. Please install it."
        exit 1
    fi
}

check_dependency seqtk
check_dependency abricate
check_dependency blastn

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --dir)
            input_dir="$2"
            shift 2
            ;;
        --threads)
            num_threads="$2"
            shift 2
            ;;
        --db)
            plsdb_path="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

if [[ -z "$input_dir" || ! -d "$input_dir" ]]; then
    echo "Error: Input directory not provided or does not exist."
    usage
fi
if [[ -z "$num_threads" || ! "$num_threads" =~ ^[0-9]+$ ]]; then
    echo "Error: Invalid or missing number of threads."
    usage
fi
if [[ -z "$plsdb_path" || ! -f "$plsdb_path" ]]; then
    echo "Error: PLSDB database path not provided or does not exist."
    usage
fi

cd "$input_dir" || { echo "Error: Unable to change to directory $input_dir"; exit 1; }

if [[ -z "$(find . -maxdepth 1 -name '*.fasta' -print -quit)" ]]; then
    echo "Error: No .fasta files found in $input_dir"
    exit 1
fi


echo "Running abricate on fasta files..."
for f in *.fasta; do
    base_name=$(basename "$f" .fasta)
    abricate --db plasmidfinder "$f" --threads "$num_threads" > "${base_name}_plasmid.txt" || { echo "Error: abricate failed for $f"; exit 1; }
done


echo "Reversing and merging fasta files..."
for file in *.fasta; do
    base_name=$(basename "$file" .fasta)
    reverse_file=$(mktemp) || { echo "Error: Failed to create temporary file"; exit 1; }
    merged_file=$(mktemp) || { echo "Error: Failed to create temporary file"; rm -f "$reverse_file"; exit 1; }

    seqtk seq -r "$file" | awk '/^>/ {gsub(/^>([^ ]+)/, ">"substr($1,2)"_R"); print} !/^>/ {print}' > "$reverse_file" || { echo "Error: seqtk failed for $file"; rm -f "$reverse_file" "$merged_file"; exit 1; }
    cat "$file" "$reverse_file" > "$merged_file" || { echo "Error: cat failed for $file"; rm -f "$reverse_file" "$merged_file"; exit 1; }

    mv "$merged_file" "${base_name}_merged.fasta" || { echo "Error: mv failed for $file"; rm -f "$reverse_file" "$merged_file"; exit 1; }
    rm -f "$reverse_file"

    echo "Processed $file: Created ${base_name}_merged.fasta"
done

echo "Running blastn on merged fasta files..."
for f in *_merged.fasta; do
    base_name=$(basename "$f" _merged.fasta)
    blastn -query "$f" -task blastn -db "$plsdb_path" -out "${base_name}_PLSDB.txt" -num_threads "$num_threads" -perc_identity 90 -qcov_hsp_perc 30 -outfmt '6 qseqid sseqid qstart qend sstart send evalue bitscore pident qcovs slen' || { echo "Error: blastn failed for $f"; exit 1; }
done

echo "Adding headers to BLASTn output files..."
header="qseqid\tsseqid\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tpident\tqcovs\tslen"
for file in *_PLSDB.txt; do
    if [[ -f "$file" ]]; then
        tmp_file=$(mktemp) || { echo "Error creating temp file"; exit 1; }
        echo -e "$header" > "$tmp_file"
        cat "$file" >> "$tmp_file"
        mv "$tmp_file" "$file" || { echo "Error: mv failed"; exit 1; }
    fi
done


echo "Modifying plasmid names in abricate output files..."
declare -A plasmid_map=()

for plasmid_file in *_plasmid.txt; do
    declare -A plasmid_count=()
    temp_file=$(mktemp) || { echo "Error creating temp file"; exit 1; }

    header=$(head -n 1 "$plasmid_file")
    echo "$header" > "$temp_file"


    while IFS= read -r line; do
        if [[ "$line" == "$header" ]]; then
            continue
        fi

        plasmid=$(echo "$line" | awk '{print $6}')
        if [[ -n "$plasmid" ]]; then
            modified_plasmid=$(echo "$plasmid" | sed 's/[()]/_/g')
            if [[ -n "${plasmid_count[$modified_plasmid]}" ]]; then
                plasmid_count[$modified_plasmid]=$((plasmid_count[$modified_plasmid] + 1))
                modified_plasmid="${modified_plasmid}_pld${plasmid_count[$modified_plasmid]}"
            else
                plasmid_count[$modified_plasmid]=1
                modified_plasmid="${modified_plasmid}_pld1"
            fi
            plasmid_map["$modified_plasmid"]=1
            escaped_original_plasmid=$(printf '%s' "$plasmid" | sed -e 's/[]\/$*.^|[]/\\&/g')
            echo "$line" | sed "s|$escaped_original_plasmid|$modified_plasmid|g" >> "$temp_file" || { echo "Error in sed substitution"; exit 1; }
        else
            echo "$line" >> "$temp_file"
        fi
    done < "$plasmid_file"

    mv "$temp_file" "$plasmid_file" || { echo "Error: mv failed"; exit 1; }
done


output_plasmid_list="plasmid_list.txt"
> "$output_plasmid_list"
for plasmid in "${!plasmid_map[@]}"; do
    echo "$plasmid" >> "$output_plasmid_list"
done

echo "Plasmid list extracted to $output_plasmid_list."

mkdir -p unicycler_fasta || { echo "Error creating directory"; exit 1; }
for file in *.fasta; do
    if [[ "$file" != *_merged.fasta ]]; then
        mv "$file" unicycler_fasta/ || { echo "Error moving file"; exit 1; }
    fi
done

echo "plsMD-preprocessing is complete"
