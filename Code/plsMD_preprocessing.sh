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

# Create directories for abricate outputs
mkdir -p plasmidfinder_results rep_typer_results merged_plasmid_results || { echo "Error creating directories"; exit 1; }

echo "Running abricate with plasmidfinder database..."
for f in *.fasta; do
    base_name=$(basename "$f" .fasta)
    abricate --db plasmidfinder "$f" --threads "$num_threads" > "plasmidfinder_results/${base_name}.txt" || { echo "Error: abricate plasmidfinder failed for $f"; exit 1; }
done

echo "Running abricate with rep.mob.typer database..."
for f in *.fasta; do
    base_name=$(basename "$f" .fasta)
    abricate --db rep.mob.typer "$f" --threads "$num_threads" > "rep_typer_results/${base_name}.txt" || { echo "Error: abricate rep.mob.typer failed for $f"; exit 1; }
done


process_sample() {
    local sample_name="$1"
    local plasmidfinder_file="plasmidfinder_results/${sample_name}.txt"
    local reptyper_file="rep_typer_results/${sample_name}.txt"
    local output_file="merged_plasmid_results/${sample_name}_plasmid.txt"
    
    if [[ ! -f "$plasmidfinder_file" ]]; then
        echo "Warning: PlasmidFinder file not found: $plasmidfinder_file"
        return 1
    fi
    
    if [[ ! -f "$reptyper_file" ]]; then
        echo "Warning: RepTyper file not found: $reptyper_file"
        return 1
    fi
    
    local pf_temp=$(mktemp)
    local rt_temp=$(mktemp)
    local output_temp=$(mktemp)
    
    local pf_header=$(head -n 1 "$plasmidfinder_file")
    local rt_header=$(head -n 1 "$reptyper_file")
    
    echo "$pf_header" > "$output_file"
    
    tail -n +2 "$plasmidfinder_file" > "$pf_temp"
    tail -n +2 "$reptyper_file" > "$rt_temp"
    
    #Filter plasmidfinder entries - remove Col plasmids when Inc plasmids exist on same contig
    declare -A contig_has_inc=()
    declare -a filtered_pf_entries=()
    
    while IFS=$'\t' read -r line; do
        [[ -z "$line" ]] && continue
        IFS=$'\t' read -r -a fields <<< "$line"
        if [[ ${#fields[@]} -lt 6 ]]; then
            filtered_pf_entries+=("$line")
            continue
        fi
        
        local sequence="${fields[1]}"
        local plasmid="${fields[5]}"
        
        if [[ "$plasmid" == Inc* ]]; then
            contig_has_inc["$sequence"]=1
        fi
    done < "$pf_temp"
  
    while IFS=$'\t' read -r line; do
        [[ -z "$line" ]] && continue
        IFS=$'\t' read -r -a fields <<< "$line"
        if [[ ${#fields[@]} -lt 6 ]]; then
            filtered_pf_entries+=("$line")
            continue
        fi
        
        local sequence="${fields[1]}"
        local plasmid="${fields[5]}"
        
        if [[ "$plasmid" == Col* ]] && [[ -n "${contig_has_inc[$sequence]}" ]]; then
            echo "Filtering out Col plasmid '$plasmid' from contig '$sequence' (Inc plasmid present)"
            continue
        fi
        
        filtered_pf_entries+=("$line")
    done < "$pf_temp"
    
    for entry in "${filtered_pf_entries[@]}"; do
        echo "$entry" >> "$output_temp"
    done
    
    declare -A pf_contigs=()
    for entry in "${filtered_pf_entries[@]}"; do
        IFS=$'\t' read -r -a fields <<< "$entry"
        if [[ ${#fields[@]} -ge 2 ]]; then
            pf_contigs["${fields[1]}"]=1
        fi
    done
    
    while IFS=$'\t' read -r line; do
        [[ -z "$line" ]] && continue
        IFS=$'\t' read -r -a fields <<< "$line"
        if [[ ${#fields[@]} -ge 2 ]]; then
            contig="${fields[1]}"
            # Only add if contig is NOT in plasmidfinder
            if [[ -z "${pf_contigs[$contig]}" ]]; then
                echo "$line" >> "$output_temp"
            else
                echo "Filtering out rep.typer entry for contig '$contig' (present in plasmidfinder)"
            fi
        fi
    done < "$rt_temp"
    
    if [[ -s "$output_temp" ]]; then
        sort -k1,1 -k3,3n "$output_temp" >> "$output_file"
    fi
    
    rm -f "$pf_temp" "$rt_temp" "$output_temp"
    
    local pf_count=${#filtered_pf_entries[@]}
    local rt_total=$(tail -n +2 "$reptyper_file" | wc -l)
    local rt_filtered=0
    while IFS=$'\t' read -r line; do
        [[ -z "$line" ]] && continue
        IFS=$'\t' read -r -a fields <<< "$line"
        if [[ ${#fields[@]} -ge 2 ]] && [[ -n "${pf_contigs[${fields[1]}]}" ]]; then
            ((rt_filtered++))
        fi
    done < <(tail -n +2 "$reptyper_file")
    local rt_kept=$((rt_total - rt_filtered))
    local total_plasmids=$((pf_count + rt_kept))
    echo "Processed $sample_name: $pf_count from plasmidfinder + $rt_kept from rep.typer = $total_plasmids total"
}

echo "Merging plasmidfinder and rep.mob.typer results (simplified approach)..."
# Main processing loop for merging
for pf_file in plasmidfinder_results/*.txt; do
    if [[ -f "$pf_file" ]]; then
        sample_name=$(basename "$pf_file" .txt)
        echo "Processing sample: $sample_name"
        process_sample "$sample_name" || echo "Error processing $sample_name, continuing..."
    fi
done

for rt_file in rep_typer_results/*.txt; do
    if [[ -f "$rt_file" ]]; then
        sample_name=$(basename "$rt_file" .txt)
        pf_file="plasmidfinder_results/${sample_name}.txt"
        
        if [[ ! -f "$pf_file" ]]; then
            echo "Processing rep.typer-only sample: $sample_name"
            head -n 1 "$rt_file" > "merged_plasmid_results/${sample_name}.txt"
            tail -n +2 "$rt_file" | sort -k1,1 -k3,3n >> "merged_plasmid_results/${sample_name}.txt"
            local entry_count=$(tail -n +2 "$rt_file" | wc -l)
            echo "  -> Copied rep.typer only file: $entry_count entries"
        fi
    fi
done

echo "Copying merged plasmid results to root directory..."
cp merged_plasmid_results/*.txt . || { echo "Error copying merged results"; exit 1; }

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

echo "Modifying plasmid names in merged plasmid results..."
declare -A plasmid_map=()

for plasmid_file in *.txt; do
    # Skip BLAST output files
    if [[ "$plasmid_file" == *"_PLSDB.txt" ]]; then
        continue
    fi
    

    if [[ "$plasmid_file" == *"/"* ]]; then
        continue
    fi
    
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
            modified_plasmid=$(echo "$plasmid" | sed 's/[()\/-]/_/g')
            if [[ -n "${plasmid_count[$modified_plasmid]}" ]]; then
                plasmid_count[$modified_plasmid]=$((plasmid_count[$modified_plasmid] + 1))
                modified_plasmid="${modified_plasmid}_pld${plasmid_count[$modified_plasmid]}"
            else
                plasmid_count[$modified_plasmid]=1
                modified_plasmid="${modified_plasmid}_pld1"
            fi
            plasmid_map["$modified_plasmid"]=1
            escaped_original_plasmid=$(printf '%s' "$plasmid" | sed -e 's/[]\/$*.^|[]/\\&/g')
            echo "$line" | sed "s~\t$plasmid\t~\t$modified_plasmid\t~g" >> "$temp_file" || { echo "Error in sed substitution"; exit 1; }
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
