#!/bin/bash

PLSDB_DEFAULT="${BLAST_DB_DIR:-/opt/plsMD/data/blastdb}/plsdb"

usage() {
    echo ""
    echo "Usage: $0 --input <fasta|fastq> --dir <input_directory> --output <output_directory> --threads <num_threads> [--db <plsdb_prefix>]"
    echo ""
    echo "  --input      : Input type: 'fasta' (default) for unicycler assemblies, or 'fastq' for raw paired-end reads."
    echo "  --dir        : Input directory."
    echo "                   fasta mode: directory containing *.fasta assembly files."
    echo "                   fastq mode: directory containing paired *_1.fastq.gz / *_2.fastq.gz files."
    echo "  --output     : Output directory for all results."
    echo "  --threads    : Number of threads for fastp, unicycler, abricate, and blastn."
    echo "  --db         : Path to the PLSDB BLAST index prefix."
    echo "                 Optional if PLSDB was downloaded during Docker build (--build-arg DOWNLOAD_DB=true)."
    echo "                 Required if using an external database on the host machine."
    echo "  -h, --help   : Display this help message."
    echo ""
    echo "Examples:"
    echo "  # FASTA mode (assemblies already available)"
    echo "  plsMD --preprocessing --input fasta --dir /data/assemblies --output /data/results --threads 8"
    echo ""
    echo "  # FASTQ mode (trim + assemble + analyse in one run)"
    echo "  plsMD --preprocessing --input fastq --dir /data/reads --output /data/results --threads 16"
    echo ""
    exit 1
}

check_dependency() {
    if ! command -v "$1" &> /dev/null; then
        echo "Error: $1 is not installed or not on PATH."
        exit 1
    fi
}

# â”€â”€ Parameter parsing â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

input_mode="fasta"   # default

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --input)
            input_mode="$2"
            shift 2
            ;;
        --dir)
            input_dir="$2"
            shift 2
            ;;
        --output)
            output_dir="$2"
            shift 2
            ;;
        --threads)
            num_threads="$2"
            shift 2
            ;;
        --db)
            plsdb_override="$2"
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

# â”€â”€ Validate common parameters â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

if [[ "$input_mode" != "fasta" && "$input_mode" != "fastq" ]]; then
    echo "Error: --input-mode must be 'fasta' or 'fastq' (got: '$input_mode')."
    usage
fi

if [[ -z "$input_dir" || ! -d "$input_dir" ]]; then
    echo "Error: Input directory not provided or does not exist."
    usage
fi
if [[ -z "$output_dir" ]]; then
    echo "Error: Output directory not provided."
    usage
fi
if [[ -z "$num_threads" || ! "$num_threads" =~ ^[0-9]+$ ]]; then
    echo "Error: Invalid or missing number of threads."
    usage
fi

mkdir -p "$output_dir" || { echo "Error: Unable to create output directory $output_dir"; exit 1; }
output_dir="$(realpath "$output_dir")"
input_dir="$(realpath "$input_dir")"

# â”€â”€ Dependency checks â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

check_dependency seqtk
check_dependency abricate
check_dependency blastn

if [[ "$input_mode" == "fastq" ]]; then
    check_dependency fastp
    check_dependency unicycler
fi

# â”€â”€ PLSDB resolution (shared by both modes) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

if [[ -n "$plsdb_override" ]]; then
    plsdb_path="$plsdb_override"
    if [[ ! -f "${plsdb_path}.nhr" ]]; then
        echo "Error: PLSDB BLAST index not found at provided path: ${plsdb_path}"
        echo "Expected file: ${plsdb_path}.nhr"
        echo "Make sure the path points to the BLAST index prefix, not a directory or fasta file."
        echo "Example: --db /data/plsdb/plsdb  (where plsdb.nhr, plsdb.nin etc. exist)"
        exit 1
    fi
    echo "Using user-provided PLSDB: $plsdb_path"
elif [[ -f "${PLSDB_DEFAULT}.nhr" ]]; then
    plsdb_path="$PLSDB_DEFAULT"
    echo "Using built-in PLSDB: $plsdb_path"
else
    echo "Error: PLSDB database not found."
    echo ""
    echo "Options:"
    echo "  1. Provide your existing database with --db:"
    echo "       plsMD --preprocessing --dir <dir> --output <out> --threads <n> --db /path/to/plsdb"
    echo "       Note: path must be the BLAST index prefix (e.g. /data/plsdb/plsdb)"
    echo ""
    echo "  2. Rebuild the Docker image with the database included:"
    echo "       docker build --build-arg DOWNLOAD_DB=true -t plsmd ."
    exit 1
fi

# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# FASTQ MODE  â€”  trim with fastp, assemble with unicycler, then hand off to
#                the shared fasta pipeline below.
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

if [[ "$input_mode" == "fastq" ]]; then

    echo "Input mode: FASTQ  â€”  trimming and assembly will run before plasmid analysis."

    trimmed_dir="$output_dir/trimmed"
    unicycler_dir="$output_dir/unicycler_assemblies"
    assembled_fasta_dir="$output_dir/assembled_fastas"

    mkdir -p "$trimmed_dir" "$unicycler_dir" "$assembled_fasta_dir"

    # â”€â”€ Auto-detect fastq naming and pair R1/R2 files â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    # Supported patterns (all auto-detected, no flag needed):
    #   <sample>_1.fastq.gz        /  <sample>_2.fastq.gz
    #   <sample>_R1.fastq.gz       /  <sample>_R2.fastq.gz
    #   <sample>_R1_NNN.fastq.gz   /  <sample>_R2_NNN.fastq.gz
    #     (covers _L001_R1_001.fastq.gz and similar Illumina conventions)

    pair_fastq_files() {
        local dir="$1"
        # Prints one tab-separated line per pair:  sample_name <TAB> r1_path <TAB> r2_path
        local paired=0
        declare -A _seen=()

        while IFS= read -r f; do
            local base sample r2 num
            base=$(basename "$f")

            [[ -n "${_seen[$base]}" ]] && continue   # already emitted as part of a pair

            if [[ "$base" =~ ^(.+)_1\.fastq\.gz$ ]]; then
                sample="${BASH_REMATCH[1]}"
                r2="$dir/${sample}_2.fastq.gz"

            elif [[ "$base" =~ ^(.+)_R1_([0-9]+)\.fastq\.gz$ ]]; then
                sample="${BASH_REMATCH[1]}"
                num="${BASH_REMATCH[2]}"
                r2="$dir/${sample}_R2_${num}.fastq.gz"

            elif [[ "$base" =~ ^(.+)_R1\.fastq\.gz$ ]]; then
                sample="${BASH_REMATCH[1]}"
                r2="$dir/${sample}_R2.fastq.gz"

            else
                continue   # not an R1 file â€” skip
            fi

            if [[ -f "$r2" ]]; then
                _seen["$base"]=1
                _seen["$(basename "$r2")"]=1
                printf '%s\t%s\t%s\n' "$sample" "$f" "$r2"
                (( paired++ ))
            else
                echo "Warning: R1 detected but no matching R2 found â€” skipping: $f" >&2
            fi

        done < <(find "$dir" -maxdepth 1 -name "*.fastq.gz" | sort)

        if [[ $paired -eq 0 ]]; then
            return 1
        fi
    }

    # Validate that at least one pair can be found
    if ! pair_fastq_files "$input_dir" > /dev/null 2>&1; then
        echo "Error: No paired fastq.gz files detected in $input_dir"
        echo "Supported naming conventions:"
        echo "  <sample>_1.fastq.gz / <sample>_2.fastq.gz"
        echo "  <sample>_R1.fastq.gz / <sample>_R2.fastq.gz"
        echo "  <sample>_R1_001.fastq.gz / <sample>_R2_001.fastq.gz  (and any _LXXX_ prefix)"
        exit 1
    fi

    # â”€â”€ Step 1: fastp trimming â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

    echo ""
    echo "Step 1/2 â€” Trimming reads with fastp..."

    while IFS=$'\t' read -r sample r1 r2; do
        echo "  Trimming: $sample  ($(basename "$r1") + $(basename "$r2"))"
        fastp \
            -w "$num_threads" \
            -i "$r1" \
            -I "$r2" \
            -o "$trimmed_dir/${sample}_1_trimmed.fastq.gz" \
            -O "$trimmed_dir/${sample}_2_trimmed.fastq.gz" \
            --cut_front \
            --cut_tail \
            --cut_window_size 4 \
            --cut_mean_quality 20 \
            --length_required 50 \
            --html "$trimmed_dir/${sample}_report.html" \
            --json "$trimmed_dir/${sample}_report.json" \
            2>&1 | sed "s/^/    [fastp] /" \
        || { echo "Error: fastp failed for $sample"; exit 1; }
    done < <(pair_fastq_files "$input_dir")

    # â”€â”€ Step 2: unicycler assembly â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

    echo ""
    echo "Step 2/2 â€” Assembling with Unicycler..."

    # Trimmed files always use the normalised _1/_2 convention after fastp
    for r1_trimmed in "$trimmed_dir"/*_1_trimmed.fastq.gz; do
        sample=$(basename "$r1_trimmed" _1_trimmed.fastq.gz)
        r2_trimmed="$trimmed_dir/${sample}_2_trimmed.fastq.gz"

        if [[ ! -f "$r2_trimmed" ]]; then
            echo "Warning: Trimmed R2 not found for $sample, skipping assembly."
            continue
        fi

        echo "  Assembling: $sample"
        unicycler \
            -1 "$r1_trimmed" \
            -2 "$r2_trimmed" \
            --threads "$num_threads" \
            -o "$unicycler_dir/$sample" \
            2>&1 | sed "s/^/    [unicycler] /" \
        || { echo "Error: unicycler failed for $sample"; exit 1; }

        assembly="$unicycler_dir/$sample/assembly.fasta"
        if [[ ! -f "$assembly" ]]; then
            echo "Error: Unicycler finished but assembly.fasta not found for $sample."
            exit 1
        fi

        cp "$assembly" "$assembled_fasta_dir/${sample}.fasta" \
            || { echo "Error: Failed to copy assembly for $sample"; exit 1; }

        echo "  Assembly complete: $sample  â†’  $assembled_fasta_dir/${sample}.fasta"
    done

    # Validate at least one assembly was produced
    if [[ -z "$(find "$assembled_fasta_dir" -maxdepth 1 -name '*.fasta' -print -quit)" ]]; then
        echo "Error: No assembly fasta files were produced. Check fastp/unicycler logs above."
        exit 1
    fi

    # Re-point input_dir so the shared fasta pipeline below uses the assemblies
    input_dir="$assembled_fasta_dir"
    echo ""
    echo "Trimming and assembly complete. Proceeding with plasmid analysis using assembled fastas..."
    echo ""

fi  # end fastq mode

# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# SHARED FASTA PIPELINE
# (reached directly in fasta mode, or after trimming+assembly in fastq mode)
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

echo "Input mode: FASTA pipeline starting..."

# Verify input directory has fasta files
if [[ -z "$(find "$input_dir" -maxdepth 1 -name '*.fasta' -print -quit)" ]]; then
    echo "Error: No .fasta files found in $input_dir"
    exit 1
fi

# Create output subdirectories
mkdir -p "$output_dir/plasmidfinder_results" \
         "$output_dir/rep_typer_results" \
         "$output_dir/merged_plasmid_results" \
         "$output_dir/unicycler_fasta" \
    || { echo "Error creating output subdirectories"; exit 1; }

echo "Running abricate with plasmidfinder database..."
for f in "$input_dir"/*.fasta; do
    base_name=$(basename "$f" .fasta)
    abricate --db plasmidfinder "$f" --threads "$num_threads" \
        > "$output_dir/plasmidfinder_results/${base_name}.txt" \
        || { echo "Error: abricate plasmidfinder failed for $f"; exit 1; }
done

echo "Running abricate with rep.mob.typer database..."
for f in "$input_dir"/*.fasta; do
    base_name=$(basename "$f" .fasta)
    abricate --db rep.mob.typer "$f" --threads "$num_threads" \
        > "$output_dir/rep_typer_results/${base_name}.txt" \
        || { echo "Error: abricate rep.mob.typer failed for $f"; exit 1; }
done


process_sample() {
    local sample_name="$1"
    local plasmidfinder_file="$output_dir/plasmidfinder_results/${sample_name}.txt"
    local reptyper_file="$output_dir/rep_typer_results/${sample_name}.txt"
    local output_file="$output_dir/merged_plasmid_results/${sample_name}_plasmid.txt"

    if [[ ! -f "$plasmidfinder_file" ]]; then
        echo "Warning: PlasmidFinder file not found: $plasmidfinder_file"
        return 1
    fi

    if [[ ! -f "$reptyper_file" ]]; then
        echo "Warning: RepTyper file not found: $reptyper_file"
        return 1
    fi

    local pf_temp rt_temp output_temp
    pf_temp=$(mktemp)
    rt_temp=$(mktemp)
    output_temp=$(mktemp)

    local pf_header
    pf_header=$(head -n 1 "$plasmidfinder_file")

    echo "$pf_header" > "$output_file"

    tail -n +2 "$plasmidfinder_file" > "$pf_temp"
    tail -n +2 "$reptyper_file"      > "$rt_temp"

    # Filter plasmidfinder entries â€” remove Col plasmids when Inc plasmids exist on the same contig
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
            local contig="${fields[1]}"
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
    local rt_total rt_filtered=0
    rt_total=$(tail -n +2 "$reptyper_file" | wc -l)

    while IFS=$'\t' read -r line; do
        [[ -z "$line" ]] && continue
        IFS=$'\t' read -r -a fields <<< "$line"
        if [[ ${#fields[@]} -ge 2 ]] && [[ -n "${pf_contigs[${fields[1]}]}" ]]; then
            ((rt_filtered++))
        fi
    done < <(tail -n +2 "$reptyper_file")

    local rt_kept=$(( rt_total - rt_filtered ))
    local total_plasmids=$(( pf_count + rt_kept ))
    echo "Processed $sample_name: $pf_count from plasmidfinder + $rt_kept from rep.typer = $total_plasmids total"
}

echo "Merging plasmidfinder and rep.mob.typer results..."
for pf_file in "$output_dir/plasmidfinder_results"/*.txt; do
    if [[ -f "$pf_file" ]]; then
        sample_name=$(basename "$pf_file" .txt)
        echo "Processing sample: $sample_name"
        process_sample "$sample_name" || echo "Error processing $sample_name, continuing..."
    fi
done

for rt_file in "$output_dir/rep_typer_results"/*.txt; do
    if [[ -f "$rt_file" ]]; then
        sample_name=$(basename "$rt_file" .txt)
        pf_file="$output_dir/plasmidfinder_results/${sample_name}.txt"

        if [[ ! -f "$pf_file" ]]; then
            echo "Processing rep.typer-only sample: $sample_name"
            head -n 1 "$rt_file" > "$output_dir/merged_plasmid_results/${sample_name}.txt"
            tail -n +2 "$rt_file" | sort -k1,1 -k3,3n \
                >> "$output_dir/merged_plasmid_results/${sample_name}.txt"
            local entry_count
            entry_count=$(tail -n +2 "$rt_file" | wc -l)
            echo "  -> Copied rep.typer only file: $entry_count entries"
        fi
    fi
done

echo "Copying merged plasmid results to output directory..."
cp "$output_dir/merged_plasmid_results/"*.txt "$output_dir/" \
    || { echo "Error copying merged results"; exit 1; }

echo "Reversing and merging fasta files..."
for file in "$input_dir"/*.fasta; do
    base_name=$(basename "$file" .fasta)
    reverse_file=$(mktemp) || { echo "Error: Failed to create temporary file"; exit 1; }
    merged_file=$(mktemp)  || { echo "Error: Failed to create temporary file"; rm -f "$reverse_file"; exit 1; }

    seqtk seq -r "$file" \
        | awk '/^>/ {gsub(/^>([^ ]+)/, ">"substr($1,2)"_R"); print} !/^>/ {print}' \
        > "$reverse_file" \
        || { echo "Error: seqtk failed for $file"; rm -f "$reverse_file" "$merged_file"; exit 1; }

    cat "$file" "$reverse_file" > "$merged_file" \
        || { echo "Error: cat failed for $file"; rm -f "$reverse_file" "$merged_file"; exit 1; }

    mv "$merged_file" "$output_dir/${base_name}_merged.fasta" \
        || { echo "Error: mv failed for $file"; rm -f "$reverse_file" "$merged_file"; exit 1; }
    rm -f "$reverse_file"

    echo "Processed $file: Created ${base_name}_merged.fasta"
done

echo "Running blastn on merged fasta files..."
for f in "$output_dir"/*_merged.fasta; do
    base_name=$(basename "$f" _merged.fasta)
    blastn \
        -query "$f" \
        -task blastn \
        -db "$plsdb_path" \
        -out "$output_dir/${base_name}_PLSDB.txt" \
        -num_threads "$num_threads" \
        -perc_identity 90 \
        -qcov_hsp_perc 30 \
        -outfmt '6 qseqid sseqid qstart qend sstart send evalue bitscore pident qcovs slen' \
        || { echo "Error: blastn failed for $f"; exit 1; }
done

echo "Adding headers to BLASTn output files..."
header="qseqid\tsseqid\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tpident\tqcovs\tslen"
for file in "$output_dir"/*_PLSDB.txt; do
    if [[ -f "$file" ]]; then
        tmp_file=$(mktemp) || { echo "Error creating temp file"; exit 1; }
        echo -e "$header" > "$tmp_file"
        cat "$file" >> "$tmp_file"
        mv "$tmp_file" "$file" || { echo "Error: mv failed"; exit 1; }
    fi
done

echo "Modifying plasmid names in merged plasmid results..."
declare -A plasmid_map=()

for plasmid_file in "$output_dir"/*.txt; do
    if [[ "$plasmid_file" == *"_PLSDB.txt" ]]; then
        continue
    fi
    if [[ "$(basename "$plasmid_file")" == *"/"* ]]; then
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

        plasmid=$(echo "$line" | awk -F'\t' '{print $6}')
        if [[ -n "$plasmid" ]]; then
            modified_plasmid=$(echo "$plasmid" | sed 's/[()\/:;|, -]/_/g')
            modified_plasmid=$(echo "$modified_plasmid" | sed 's/^_*//;s/_*$//')
            if [[ -n "${plasmid_count[$modified_plasmid]}" ]]; then
                plasmid_count[$modified_plasmid]=$(( plasmid_count[$modified_plasmid] + 1 ))
                modified_plasmid="${modified_plasmid}_pld${plasmid_count[$modified_plasmid]}"
            else
                plasmid_count[$modified_plasmid]=1
                modified_plasmid="${modified_plasmid}_pld1"
            fi
            plasmid_map["$modified_plasmid"]=1
            echo "$line" | awk -F'\t' -v old="$plasmid" -v new="$modified_plasmid" \
                'BEGIN{OFS="\t"} {if ($6 == old) $6 = new; print}' >> "$temp_file" \
                || { echo "Error in awk substitution"; exit 1; }
        else
            echo "$line" >> "$temp_file"
        fi
    done < "$plasmid_file"

    mv "$temp_file" "$plasmid_file" || { echo "Error: mv failed"; exit 1; }
done

output_plasmid_list="$output_dir/plasmid_list.txt"
> "$output_plasmid_list"
for plasmid in "${!plasmid_map[@]}"; do
    echo "$plasmid" >> "$output_plasmid_list"
done

echo "Plasmid list extracted to $output_plasmid_list."

# Copy fasta files (assemblies) to unicycler_fasta in output dir
for file in "$input_dir"/*.fasta; do
    cp "$file" "$output_dir/unicycler_fasta/" || { echo "Error copying fasta file $file"; exit 1; }
done

echo "plsMD-preprocessing is complete. All outputs written to: $output_dir"
