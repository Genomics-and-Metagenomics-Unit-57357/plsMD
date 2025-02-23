import os
import re
import shutil
import pandas as pd
from Bio import SeqIO, Entrez
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict
import argparse

def filter_invalid_rows(file_path):
    """Filters out rows in a file where sstart > send."""
    with open(file_path, 'r') as file:
        lines = file.readlines()

    header = lines[0]
    data_lines = lines[1:]

    filtered_lines = [line for line in data_lines if int(line.split('\t')[4]) <= int(line.split('\t')[5])]

    with open(file_path, 'w') as file:
        file.write(header)
        file.writelines(filtered_lines)

def extract_contig_lengths_from_fasta(fasta_file):
    """Extracts the lengths of contigs from a FASTA file."""
    contig_lengths = {}
    with open(fasta_file, 'r') as fasta_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            contig_lengths[record.id] = len(record.seq)
    return contig_lengths

def calculate_query_percentage_aligned(row):
    """Calculates the query percentage aligned based on qstart, qend, and qlen."""
    return (int(row['qend']) - int(row['qstart']) + 1) / int(row['qlen'])

def process_plsdb_file(file_path, fasta_file):
    """Processes a single _PLSDB.txt file."""
    print(f"Processing {file_path}")
    filter_invalid_rows(file_path)

    if not os.path.exists(fasta_file):
        print(f"FASTA file {fasta_file} not found. Skipping {file_path}.")
        return

    print(f"Extracting contig lengths from {fasta_file}")
    contig_lengths = extract_contig_lengths_from_fasta(fasta_file)

    plsdb_df = pd.read_csv(file_path, sep='\t', dtype=str)
    plsdb_df['qlen'] = plsdb_df['qseqid'].apply(lambda x: contig_lengths.get(x.strip(), 'N/A'))
    plsdb_df['q_perc_aligned'] = plsdb_df.apply(calculate_query_percentage_aligned, axis=1)
    plsdb_df.to_csv(file_path, sep='\t', index=False)
    print(f"Updated {file_path} with qlen and q_perc_aligned columns.")

def process_directory(directory_path):
    """Processes all _PLSDB.txt files in a directory."""
    for file in os.listdir(directory_path):
        if file.endswith("_PLSDB.txt"):
            file_path = os.path.join(directory_path, file)
            base_name = file.replace('_PLSDB.txt', '')
            fasta_file = os.path.join(directory_path, f"{base_name}_merged.fasta")
            process_plsdb_file(file_path, fasta_file)

def process_contigs(input_file, output_file, nested_contigs_log, unique_nested_contigs):
    """Process the filtered gene contigs to handle nested and overlapping contigs efficiently."""
    df = pd.read_csv(input_file, sep='\t')
    sample_name = os.path.basename(input_file).replace('_PLSDB.txt', '')
    df.sort_values(by=['sseqid', 'sstart', 'send'], inplace=True)
    plasmid_groups = df.groupby('sseqid')
    final_contigs = []
    nested_contigs = []

    for plasmid, group in plasmid_groups:
        group = group.reset_index(drop=True)
        group['max_send'] = group['send'].cummax().shift(fill_value=group['send'].iloc[0])
        nested_flags = (group['sstart'] >= group['sstart'].shift()) & (group['send'] <= group['max_send'])
        nested_ids = group.loc[nested_flags, 'qseqid'].tolist()
        for contig_id in nested_ids:
            if contig_id not in unique_nested_contigs:
                unique_nested_contigs.add(contig_id)
                nested_contigs.append((sample_name, contig_id))
        filtered_group = group[~nested_flags].copy()
        filtered_group['next_sstart'] = filtered_group['sstart'].shift(-1, fill_value=filtered_group['sstart'].max())
        filtered_group['next_qstart'] = filtered_group['qstart'].shift(-1, fill_value=0)
        filtered_group['overlapped_bases'] = 0
        overlap_condition = (filtered_group['send'] > filtered_group['next_sstart']) & \
                            (filtered_group['next_qstart'].between(1, 10)) & \
                            (filtered_group['qlen'] > filtered_group['qend'])
        filtered_group.loc[overlap_condition, 'overlapped_bases'] = (filtered_group['send'] - filtered_group['next_sstart'] + 1)
        filtered_group.loc[overlap_condition, 'send'] = filtered_group['next_sstart'] - 1
        filtered_group.drop(columns=['next_sstart', 'next_qstart', 'max_send'], inplace=True)
        final_contigs.append(filtered_group)

    final_df = pd.concat(final_contigs)
    if 'overlapped_bases' not in final_df.columns:
        final_df['overlapped_bases'] = 0
    final_df.to_csv(output_file, sep='\t', index=False)

    if nested_contigs:
        with open(nested_contigs_log, 'a') as f:
            for sample, contig_id in nested_contigs:
                f.write(f"{sample}\t{contig_id}\n")

    print(f"Updated filtered contigs saved to {output_file}")
    print(f"Nested contigs logged in {nested_contigs_log}")

def process_files_in_directory(directory_path):
    """Iterate over files in the directory and process _PLSDB.txt files."""
    nested_contigs_log = os.path.join(directory_path, 'nested_contigs.txt')
    with open(nested_contigs_log, 'w') as f:
        f.write("Sample\tContig_ID\n")

    unique_nested_contigs = set()
    for file in os.listdir(directory_path):
        if file.endswith('_PLSDB.txt'):
            input_file = os.path.join(directory_path, file)
            output_file_name = f"{file[:-len('_PLSDB.txt')]}_overlap.txt"
            output_file = os.path.join(directory_path, output_file_name)
            process_contigs(input_file, output_file, nested_contigs_log, unique_nested_contigs)

def add_order_column(data):
    """Add an 'order' column for indexing within each sseqid group."""
    data['order'] = data.groupby('sseqid').cumcount()
    return data

def handle_non_overlapping_contigs(group):
    """Handle non-overlapping contigs within a group by updating q_perc_aligned and marking rows for removal."""
    to_remove = set()
    for i in range(len(group) - 1):
        rep_contig1 = group.iloc[i]
        for j in range(i + 1, len(group)):
            rep_contig2 = group.iloc[j]
            if rep_contig1['qseqid'] != rep_contig2['qseqid']:
                break
            if (rep_contig1['qend'] + rep_contig1['overlapped_bases']) < rep_contig2['qstart']:
                total_aligned = rep_contig1['q_perc_aligned'] + rep_contig2['q_perc_aligned']
                group.at[rep_contig1.name, 'q_perc_aligned'] = total_aligned
                group.at[rep_contig2.name, 'q_perc_aligned'] = total_aligned
                if rep_contig1['order'] < rep_contig2['order']:
                    to_remove.update(group.loc[(group['order'] > rep_contig1['order']) & 
                                               (group['order'] < rep_contig2['order'])].index)
                else:
                    to_remove.update(group.loc[(group['order'] < rep_contig2['order']) | 
                                               (group['order'] > rep_contig1['order'])].index)
                break
    return group, to_remove

def handle_overlapping_contigs(group, to_remove):
    """Handle overlapping contigs within a group by removing lower q_perc_aligned rows below a threshold."""
    for i in range(len(group) - 1):
        rep_contig1 = group.iloc[i]
        for j in range(i + 1, len(group)):
            rep_contig2 = group.iloc[j]
            if rep_contig1['qseqid'].replace('_R', '') != rep_contig2['qseqid'].replace('_R', ''):
                break
            if rep_contig1['qend'] > rep_contig2['qstart']:
                if rep_contig1['q_perc_aligned'] == rep_contig2['q_perc_aligned']:
                    continue
                if (rep_contig1['q_perc_aligned'] < rep_contig2['q_perc_aligned'] and
                        rep_contig1['q_perc_aligned'] < 0.8):
                    to_remove.add(rep_contig1.name)
                elif (rep_contig2['q_perc_aligned'] < rep_contig1['q_perc_aligned'] and
                      rep_contig2['q_perc_aligned'] < 0.8):
                    to_remove.add(rep_contig2.name)
    return group, to_remove

def process_group(group):
    """Process a group of data by handling both non-overlapping and overlapping contigs."""
    group = group.sort_values(by=['qseqid', 'qstart'])
    group, to_remove = handle_non_overlapping_contigs(group)
    group, to_remove = handle_overlapping_contigs(group, to_remove)
    return group.drop(to_remove)

def process_file(input_file, output_file):
    """Process a single input file and save the results to an output file."""
    data = pd.read_csv(input_file, sep='\t')
    data = add_order_column(data)
    processed_data = data.groupby('sseqid').apply(process_group).reset_index(drop=True)
    final_data = processed_data.sort_values(by=['sseqid', 'sstart']).drop(columns=['order'])
    final_data.to_csv(output_file, sep='\t', index=False)
    print(f"Processed file saved to {output_file}")

def process_directory_filtered(input_directory):
    """Process all '_overlap.txt' files in the given directory."""
    for filename in os.listdir(input_directory):
        if filename.endswith('_overlap.txt'):
            input_file = os.path.join(input_directory, filename)
            output_file = os.path.join(input_directory, filename.replace('_overlap.txt', '_overlap_filtered.txt'))
            process_file(input_file, output_file)

def get_circular_contigs(fasta_file):
    """Find circular contigs in the FASTA file."""
    circular_contigs = set()
    for record in SeqIO.parse(fasta_file, "fasta"):
        if "circular=true" in record.description:
            circular_contigs.add(record.id)
    return circular_contigs

def standardize_qseqid(qseqid):
    """Standardize qseqid by stripping '_R' suffix."""
    return qseqid.rstrip('_R')

def process_circular_contigs(input_dir, output_dir):
    """Process circular contigs in the given directory."""
    os.makedirs(output_dir, exist_ok=True)
    for file in os.listdir(input_dir):
        if file.endswith("_overlap_filtered.txt"):
            sample_name = file.replace("_overlap_filtered.txt", "")
            fasta_file = os.path.join(input_dir, f"{sample_name}_merged.fasta")
            if not os.path.exists(fasta_file):
                print(f"FASTA file not found for {sample_name}, skipping...")
                continue
            filtered_file_path = os.path.join(input_dir, file)
            df = pd.read_csv(filtered_file_path, sep="\t")
            circular_contigs = get_circular_contigs(fasta_file)
            if not circular_contigs:
                print(f"No circular contigs found for {sample_name}. Saving unfiltered file.")
                output_file = os.path.join(output_dir, file)
                df.to_csv(output_file, sep="\t", index=False)
                continue
            print(f"Circular contigs found for {sample_name}: {circular_contigs}")
            plasmids_to_exclude = set()
            sseqid_groups = df.groupby('sseqid')
            for sseqid, group in sseqid_groups:
                standardized_qseqids = group['qseqid'].apply(standardize_qseqid).unique()
                if len(standardized_qseqids) > 1:
                    if any(qseqid in circular_contigs for qseqid in standardized_qseqids):
                        plasmids_to_exclude.add(sseqid)
                elif len(standardized_qseqids) == 1:
                    single_qseqid = standardized_qseqids[0]
                    if single_qseqid in circular_contigs:
                        continue
                    elif any(qseqid in circular_contigs for qseqid in group['qseqid']):
                        plasmids_to_exclude.add(sseqid)
            print(f"Plasmids to exclude for {sample_name}: {plasmids_to_exclude}")
            filtered_df = df[~df['sseqid'].isin(plasmids_to_exclude)]
            output_file = os.path.join(output_dir, file)
            filtered_df.to_csv(output_file, sep="\t", index=False)
            print(f"Processed {file}: Excluded {len(plasmids_to_exclude)} plasmids. Saved to {output_file}.")

def merge_intervals(intervals):
    """Merge overlapping intervals."""
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    merged_intervals = []
    for interval in sorted_intervals:
        if not merged_intervals or merged_intervals[-1][1] < interval[0] - 1:
            merged_intervals.append(list(interval))
        else:
            merged_intervals[-1][1] = max(merged_intervals[-1][1], interval[1])
    return merged_intervals

def calculate_coverage_percentage(merged_intervals, subject_length):
    """Calculate the percentage of coverage."""
    if subject_length == 0:
        return 0
    covered_length = sum(end - start + 1 for start, end in merged_intervals)
    return (covered_length / subject_length) * 100

def process_files_coverage(input_directory, coverage_cutoff=0.0):
    """Process files to calculate coverage percentage."""
    for filename in sorted(os.listdir(input_directory)):
        if filename.endswith('_overlap_filtered.txt'):
            input_file = os.path.join(input_directory, filename)
            print(f"Processing file: {filename}")
            with open(input_file, 'r') as infile:
                lines = infile.readlines()
            header = "qseqid\tsseqid\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tpident\tqcovs\tslen\tqlen\tq_perc_aligned\toverlapped_bases"
            if lines[0].strip() != header:
                lines.insert(0, header + '\n')
                with open(input_file, 'w') as infile:
                    infile.writelines(lines)
            plasmid_alignments = defaultdict(list)
            plasmid_lengths = {}
            with open(input_file, 'r') as infile:
                next(infile)
                for line in infile:
                    columns = line.strip().split('\t')
                    try:
                        qseqid, sseqid, qstart, qend, sstart, send, evalue, bitscore, pident, qcovs, slen = columns[:11]
                        sstart, send, slen = int(sstart), int(send), int(slen)
                        plasmid_alignments[sseqid].append((sstart, send))
                        plasmid_lengths[sseqid] = slen
                    except ValueError as ve:
                        print(f"Skipping line due to ValueError: {ve}")
                        continue
            plasmid_coverage = {}
            for sseqid, intervals in plasmid_alignments.items():
                merged_intervals = merge_intervals(intervals)
                coverage_percentage = calculate_coverage_percentage(merged_intervals, plasmid_lengths[sseqid])
                plasmid_coverage[sseqid] = coverage_percentage
            all_lines = []
            with open(input_file, 'r') as infile:
                next(infile)
                for line in infile:
                    columns = line.strip().split('\t')
                    sseqid = columns[1]
                    coverage_percentage = plasmid_coverage.get(sseqid, 0)
                    slen = int(columns[10])
                    bases_covered = slen * coverage_percentage / 100
                    if coverage_percentage >= coverage_cutoff:
                        all_lines.append((sseqid, line.strip(), coverage_percentage, bases_covered))
            all_lines.sort(key=lambda x: x[0])
            with open(input_file, 'w') as outfile:
                header_line = f"{header}\tcoverage_percentage\tbases_covered\n"
                outfile.write(header_line)
                for _, line, coverage_percentage, bases_covered in all_lines:
                    outfile.write(f"{line}\t{coverage_percentage:.2f}%\t{bases_covered:.2f}\n")
            print(f"Processed and updated {input_file}.")

def process_nonreplicon_contigs(input_dir, output_dir):
    """Process non-replicon contigs."""
    os.makedirs(output_dir, exist_ok=True)
    for file in os.listdir(input_dir):
        if file.endswith("_overlap_filtered.txt"):
            sample_name = file.split("_overlap_filtered")[0]
            overlap_file_path = os.path.join(input_dir, file)
            plasmid_file_path = os.path.join(input_dir, f"{sample_name}_plasmid.txt")
            if not os.path.exists(plasmid_file_path):
                print(f"Plasmid file for {sample_name} not found, skipping...")
                continue
            overlap_df = pd.read_csv(overlap_file_path, sep="\t")
            overlap_df['qseqid'] = overlap_df['qseqid'].astype(str)
            overlap_df['qseqid_normalized'] = overlap_df['qseqid'].str.replace("_R", "", regex=False)
            overlap_df['coverage_percentage'] = overlap_df['coverage_percentage'].str.replace('%', '').astype(float)
            overlap_df_filtered = overlap_df[overlap_df['coverage_percentage'] >= 80.0]
            plasmid_df = pd.read_csv(plasmid_file_path, sep="\t")
            plasmid_df['SEQUENCE'] = plasmid_df['SEQUENCE'].astype(str)
            plasmid_df['SEQUENCE_normalized'] = plasmid_df['SEQUENCE'].str.replace("_R", "", regex=False)
            non_replicon_contigs = overlap_df_filtered[~overlap_df_filtered['qseqid_normalized'].isin(plasmid_df['SEQUENCE_normalized'])]
            unique_sseqids = non_replicon_contigs['sseqid'].unique()
            expanded_non_replicon_contigs = overlap_df_filtered[overlap_df_filtered['sseqid'].isin(unique_sseqids)]
            expanded_non_replicon_contigs = expanded_non_replicon_contigs.drop(columns=['qseqid_normalized'])
            output_file_path = os.path.join(output_dir, f"{sample_name}_nonreplicon_contigs.txt")
            expanded_non_replicon_contigs.to_csv(output_file_path, sep="\t", index=False)
            print(f"Processed {sample_name}: Saved to {output_file_path}")

def extract_gene_sample_mapping(input_dir, gene_directories):
    """Extract gene-sample mapping and create directories."""
    gene_samples = {}
    for filename in os.listdir(input_dir):
        if filename.endswith("_plasmid.txt"):
            sample = filename.split("_")[0]
            with open(os.path.join(input_dir, filename), 'r') as file:
                next(file)
                for line in file:
                    gene = line.split()[5]
                    sanitized_gene = ''.join(char if char.isalnum() else '_' for char in gene)
                    if sanitized_gene not in gene_samples:
                        gene_samples[sanitized_gene] = []
                    if sample not in gene_samples[sanitized_gene]:
                        gene_samples[sanitized_gene].append(sample)
    os.makedirs(gene_directories, exist_ok=True)
    for gene, samples in gene_samples.items():
        gene_dir = os.path.join(gene_directories, gene)
        os.makedirs(gene_dir, exist_ok=True)
        print(f"Created directory: {gene_dir}")
        for sample in samples:
            sample_dir = os.path.join(gene_dir, sample)
            os.makedirs(sample_dir, exist_ok=True)
            print(f"Created sample directory: {sample_dir}")
            for ext in ["overlap_filtered.txt", "_plasmid.txt"]:
                for root, _, files in os.walk(input_dir):
                    for file in files:
                        if file.startswith(sample) and file.endswith(ext):
                            src_file = os.path.join(root, file)
                            dest_file = os.path.join(sample_dir, file)
                            if src_file != dest_file:
                                shutil.copy2(src_file, sample_dir)
                                print(f"Copied file: {src_file} to {sample_dir}")
def load_gene_list(file_path):
    """Load gene list from a text file."""
    with open(file_path, 'r') as file:
        gene_list = [line.strip() for line in file if line.strip()]
    return gene_list

def match_contig_id(contig_id, qseqid):
    """Match contig IDs with or without '_R'."""
    return re.fullmatch(f"{re.escape(contig_id)}(_R)?", qseqid) is not None

def process_gene_directories(top_level_directory, gene_list):
    """Process gene directories and extract replicon contigs."""
    for gene_name in gene_list:
        gene_dir_path = os.path.join(top_level_directory, gene_name)
        if os.path.exists(gene_dir_path) and os.path.isdir(gene_dir_path):
            print(f"Processing gene directory: {gene_name}")
            for root, dirs, files in os.walk(gene_dir_path):
                directory_name = os.path.basename(root)
                print(f"Processing subdirectory: {directory_name}")
                abricate_file = next((f for f in files if f.endswith('_plasmid.txt')), None)
                overlap_file = next((f for f in files if f.endswith('_overlap_filtered.txt')), None)
                if abricate_file and overlap_file:
                    print(f"Abricate file: {abricate_file}, Overlap file: {overlap_file}")
                    abricate_path = os.path.join(root, abricate_file)
                    overlap_path = os.path.join(root, overlap_file)
                    gene_output_path = os.path.join(root, f"{directory_name}_replicon.txt")
                    all_plasmids_output_path = os.path.join(root, f"{directory_name}_replicon_contigs.txt")
                    try:
                        abricate_df = pd.read_csv(abricate_path, sep='\t', dtype=str)
                        contig_info = abricate_df[abricate_df['GENE'].str.strip() == gene_name]
                        if not contig_info.empty:
                            for _, row in contig_info.iterrows():
                                contig_id = str(row['SEQUENCE']).strip()
                                overlap_df = pd.read_csv(overlap_path, sep='\t', dtype=str)
                                overlap_df['qseqid'] = overlap_df['qseqid'].str.strip()
                                plasmids_on_contig = overlap_df[overlap_df['qseqid'].apply(lambda x: match_contig_id(contig_id, x))]
                                if not plasmids_on_contig.empty:
                                    if os.path.exists(gene_output_path):
                                        plasmids_on_contig.to_csv(gene_output_path, sep='\t', index=False, mode='a', header=False)
                                    else:
                                        plasmids_on_contig.to_csv(gene_output_path, sep='\t', index=False)
                        if os.path.exists(gene_output_path):
                            gene_df = pd.read_csv(gene_output_path, sep='\t', dtype=str)
                            plasmid_ids = gene_df['sseqid'].unique()
                            all_plasmids_df = overlap_df[overlap_df['sseqid'].isin(plasmid_ids)]
                            all_plasmids_df['gene_name'] = gene_name
                            all_plasmids_df.to_csv(all_plasmids_output_path, sep='\t', index=False)
                            print(f"All occurrences of the identified plasmids have been saved to {all_plasmids_output_path} with gene name.")
                        else:
                            print(f"No plasmids were identified on the specific contig in {directory_name}.")
                        if os.path.exists(gene_output_path):
                            os.remove(gene_output_path)
                            print(f"Deleted the gene output file: {gene_output_path}")
                    except Exception as e:
                        print(f"An error occurred while processing {directory_name}: {e}")
                else:
                    print(f"Required files not found in subdirectory {directory_name}.")
        else:
            print(f"Gene directory not found: {gene_name}")

def process_replicon_filtered_files(top_level_directory):
    """Process files ending with '_replicon_contigs.txt'."""
    for root, dirs, files in os.walk(top_level_directory):
        directory_name = os.path.basename(root)
        print(f"Processing directory: {directory_name}")
        for file_name in files:
            if file_name.endswith('_replicon_contigs.txt'):
                file_path = os.path.join(root, file_name)
                print(f"Opening file: {file_path}")
                try:
                    df = pd.read_csv(file_path, delimiter='\t')
                    required_columns = ['gene_name', 'coverage_percentage', 'pident', 'bases_covered', 'sseqid', 'sstart']
                    if not all(col in df.columns for col in required_columns):
                        print(f"Missing required columns in file {file_path}. Skipping.")
                        continue
                    df['coverage_percentage'] = pd.to_numeric(df['coverage_percentage'].str.replace('%', ''), errors='coerce')
                    df['pident'] = pd.to_numeric(df['pident'], errors='coerce')
                    df['bases_covered'] = pd.to_numeric(df['bases_covered'], errors='coerce')
                    df = df.dropna(subset=['coverage_percentage', 'pident', 'bases_covered'])
                    if df.empty:
                        print(f"File {file_path} is empty after filtering NaN values. Skipping.")
                        continue
                    gene_name = df['gene_name'].values[0]
                    if gene_name.startswith('Col'):
                        print(f"Gene {gene_name} starts with 'Col', selecting plasmid based on highest coverage_percentage.")
                        max_coverage = df['coverage_percentage'].max()
                        highest_coverage_df = df[df['coverage_percentage'] == max_coverage]
                        if len(highest_coverage_df) > 1:
                            max_pident = highest_coverage_df['pident'].max()
                            highest_coverage_df = highest_coverage_df[highest_coverage_df['pident'] == max_pident]
                        if len(highest_coverage_df) > 1:
                            best_plasmid = highest_coverage_df.loc[highest_coverage_df['bases_covered'].idxmax()]
                        else:
                            best_plasmid = highest_coverage_df.iloc[0]
                        best_plasmid_id = best_plasmid['sseqid']
                        filtered_df = df[df['sseqid'] == best_plasmid_id]
                    else:
                        print(f"Gene {gene_name} does not start with 'Col'. Applying progressive coverage filter.")
                        thresholds = [80, 70, 60, 50, 40, 30, 20, 10]
                        for threshold in thresholds:
                            filtered_df = df[df['coverage_percentage'] >= threshold]
                            if not filtered_df.empty:
                                print(f"Found plasmids with coverage_percentage >= {threshold}.")
                                break
                        else:
                            print(f"No plasmids found with coverage_percentage >= 20. Skipping file {file_path}.")
                            continue
                        aggregated_df = filtered_df.groupby('sseqid').agg(
                            total_bases_covered=('bases_covered', 'first'),
                            avg_coverage_percentage=('coverage_percentage', 'first')
                        ).reset_index()
                        max_bases_covered = aggregated_df['total_bases_covered'].max()
                        max_coverage_percentage = aggregated_df['avg_coverage_percentage'].max()
                        aggregated_df['score_bases_covered'] = (aggregated_df['total_bases_covered'] / max_bases_covered) * 100
                        aggregated_df['score_coverage_percentage'] = (aggregated_df['avg_coverage_percentage'] / max_coverage_percentage) * 100
                        aggregated_df['combined_score'] = (aggregated_df['score_bases_covered'] + aggregated_df['score_coverage_percentage']) / 2
                        best_plasmid = aggregated_df.sort_values(by='combined_score', ascending=False).iloc[0]
                        best_plasmid_id = best_plasmid['sseqid']
                        filtered_df = df[df['sseqid'] == best_plasmid_id]
                    print(f"Number of rows after filtering for plasmid '{best_plasmid_id}': {len(filtered_df)}")
                    sorted_df = filtered_df.sort_values(by='sstart')
                    sorted_file_name = file_name.replace('_replicon_contigs.txt', '_replicon_filtered.txt')
                    sorted_file_path = os.path.join(root, sorted_file_name)
                    sorted_df.to_csv(sorted_file_path, sep='\t', index=False)
                    print(f"Processed and saved: {sorted_file_path}")
                except Exception as e:
                    print(f"Error processing file {file_path}: {e}")

def concatenate_sample_files(top_level_directory, output_directory):
    """Concatenate the '_replicon_filtered.txt' files for the same sample across all gene directories."""
    # Ensure the output directory exists
    os.makedirs(output_directory, exist_ok=True)

    sample_files = {}
    for gene_dir in os.listdir(top_level_directory):
        gene_path = os.path.join(top_level_directory, gene_dir)
        if os.path.isdir(gene_path):
            for root, dirs, files in os.walk(gene_path):
                for file in files:
                    if file.endswith('_replicon_filtered.txt'):
                        sample_name = os.path.basename(root)
                        file_path = os.path.join(root, file)
                        df = pd.read_csv(file_path, sep='\t')
                        if sample_name not in sample_files:
                            sample_files[sample_name] = df
                        else:
                            sample_files[sample_name] = pd.concat([sample_files[sample_name], df], ignore_index=True)
    for sample_name, combined_df in sample_files.items():
        output_file = os.path.join(output_directory, f"{sample_name}_plasmid_replicon_filtered.txt")
        combined_df.to_csv(output_file, sep='\t', index=False)
        print(f"Concatenated file saved: {output_file}")

def trim_contig(sequence, overlap):
    """Trim the sequence based on the overlap."""
    return sequence[:-overlap] if overlap > 0 else sequence

def extract_and_update_contigs(input_txt, fasta_file, output_fasta):
    """Extract and update contigs from FASTA file based on the filtered updated file."""
    filtered_df = pd.read_csv(input_txt, sep='\t')
    filtered_df['overlapped_bases'] = filtered_df['overlapped_bases'].fillna(0).astype(int)
    filtered_df['qseqid'] = filtered_df['qseqid'].astype(str).str.strip()
    filtered_df.columns = [filtered_df.columns[0], filtered_df.columns[1], 'qstart', 'qend', *filtered_df.columns[4:]]
    contig_dict = {}
    with open(fasta_file, "r") as fasta_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            contig_id = record.id.strip()
            sequence = str(record.seq)
            contig_dict[contig_id] = sequence
    with open(output_fasta, "w") as output_handle:
        previous_contig_id = None
        previous_qend = None
        contigs_to_write = []
        for _, row in filtered_df.iterrows():
            contig_id = row['qseqid']
            qstart = row['qstart']
            qend = row['qend']
            overlap = int(row['overlapped_bases'])
            if contig_id in contig_dict:
                if contig_id == previous_contig_id:
                    if qstart > previous_qend:
                        continue
                sequence = contig_dict[contig_id]
                trimmed_sequence = trim_contig(sequence, overlap)
                contigs_to_write.append((contig_id, trimmed_sequence))
                previous_contig_id = contig_id
                previous_qend = qend
            else:
                print(f"Contig {contig_id} not found in FASTA file.")
        if len(contigs_to_write) == 1:
            contig_id, trimmed_sequence = contigs_to_write[0]
            output_handle.write(f">{contig_id}\n{trimmed_sequence}\n")
        elif len(contigs_to_write) > 1:
            first_contig_id = contigs_to_write[0][0]
            first_contig_qstart = filtered_df.loc[filtered_df['qseqid'] == first_contig_id, 'qstart'].values[0]
            last_contig = contigs_to_write[-1]
            last_contig_id = last_contig[0]
            last_contig_qend = filtered_df.loc[filtered_df['qseqid'] == last_contig_id, 'qend'].values[-1]
            print(f"First contig qstart: {first_contig_qstart}, Last contig qend: {last_contig_qend}")
            if last_contig_qend < first_contig_qstart:
                print(f"Removing last contig due to overlap: {last_contig_id}")
                contigs_to_write.pop()
            for contig_id, trimmed_sequence in contigs_to_write:
                output_handle.write(f">{contig_id}\n{trimmed_sequence}\n")
    print(f"Extracted contigs saved to {output_fasta}")

def process_directories_for_extraction(top_level_directory, fasta_directory):
    """Iterate over directories and process files in them."""
    for root, dirs, files in os.walk(top_level_directory):
        for file in files:
            if file.endswith("_replicon_filtered.txt"):
                input_txt = os.path.join(root, file)
                basename = os.path.basename(root)

                fasta_file = None
                for f in os.listdir(fasta_directory):
                    if f.startswith(basename) and f.endswith(".fasta") and basename in f: 
                        fasta_file = os.path.join(fasta_directory, f)
                        break

                if fasta_file and os.path.exists(fasta_file): 
                    output_fasta = os.path.join(root, f"{basename}_extracted.fasta")
                    extract_and_update_contigs(input_txt, fasta_file, output_fasta)
                else:
                    print(f"FASTA file not found for {basename} in {fasta_directory}")

def merge_contigs(input_fasta, output_fasta, sample_name, plasmid_name):
    """Merge all contigs into a single sequence and write to a new FASTA file."""
    merged_sequence = ""
    with open(input_fasta, "r") as fasta_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            merged_sequence += str(record.seq)
    new_header = f">{sample_name}_{plasmid_name}"
    with open(output_fasta, "w") as output_handle:
        output_handle.write(f"{new_header}\n{merged_sequence}\n")
    print(f"Merged sequence saved to {output_fasta}")

def process_directories_for_merging(top_level_directory):
    """Process directories to merge contigs in extracted FASTA files."""
    for root, dirs, files in os.walk(top_level_directory):
        directory_name = os.path.basename(root)
        for file in files:
            if file.endswith("_extracted.fasta"):
                input_fasta = os.path.join(root, file)
                overlapped_file = None
                for txt_file in files:
                    if txt_file.endswith("_replicon_filtered.txt"):
                        overlapped_file = os.path.join(root, txt_file)
                        break
                if overlapped_file:
                    df = pd.read_csv(overlapped_file, sep='\t')
                    if 'sseqid' in df.columns and not df.empty:
                        plasmid_name = df['sseqid'].iloc[0]
                    else:
                        print(f"'sseqid' column not found or file is empty: {overlapped_file}")
                        continue
                    sample_name = directory_name
                    output_fasta = os.path.join(root, f"{sample_name}_{plasmid_name}.fasta")
                    merge_contigs(input_fasta, output_fasta, sample_name, plasmid_name)
                else:
                    print(f"No overlapped file found in {root}")

def load_gene_list_from_file(gene_list_file):
    """Load gene names from a text file into a list."""
    with open(gene_list_file, 'r') as f:
        gene_list = [line.strip() for line in f if line.strip()]
    return gene_list

def copy_merged_fasta_files(top_level_directory, destination_directory, gene_list):
    """Copy all specified fasta files into directories specified in gene_list."""
    for root, dirs, files in os.walk(top_level_directory):
        parent_directory = os.path.basename(os.path.dirname(root))
        for file in files:
            if not file.endswith("_extracted.fasta") and file.endswith(".fasta"):
                if parent_directory in gene_list:
                    destination_subdir = os.path.join(destination_directory, parent_directory)
                    if not os.path.exists(destination_subdir):
                        os.makedirs(destination_subdir)
                    source_file = os.path.join(root, file)
                    destination_file = os.path.join(destination_subdir, file)
                    shutil.copy(source_file, destination_file)
                else:
                    print(f"{parent_directory} is not in the gene list, skipping.")

def concatenate_fasta_files(top_level_directory, output_directory):
    """Concatenate the '_extracted.fasta' files for the same sample across all gene directories."""
    sample_fasta_content = {}
    for gene_dir in os.listdir(top_level_directory):
        gene_path = os.path.join(top_level_directory, gene_dir)
        if os.path.isdir(gene_path):
            for root, dirs, files in os.walk(gene_path):
                for file in files:
                    if not file.endswith("_extracted.fasta") and file.endswith(".fasta"):
                        sample_name = os.path.basename(root)
                        file_path = os.path.join(root, file)
                        with open(file_path, 'r') as fasta_file:
                            fasta_records = fasta_file.read().split('>')[1:]
                        if sample_name not in sample_fasta_content:
                            sample_fasta_content[sample_name] = {}
                        for record in fasta_records:
                            if record.strip():
                                header, sequence = record.split('\n', 1)
                                header = header.strip()
                                sequence = sequence.replace('\n', '').strip()
                                if header not in sample_fasta_content[sample_name]:
                                    sample_fasta_content[sample_name][header] = sequence
    for sample_name, records in sample_fasta_content.items():
        output_file = os.path.join(output_directory, f"{sample_name}_plasmid_contigs.fasta")
        with open(output_file, 'w') as output_fasta_file:
            for header, sequence in records.items():
                output_fasta_file.write(f">{header}\n{sequence}\n")
        print(f"Concatenated FASTA file saved: {output_file}")
#
def load_exclusion_list(exclusion_file):
    exclusion_dict = {}
    with open(exclusion_file, "r") as f:
        for line in f:
            sample, contig_id = line.strip().split('\t')
            if sample not in exclusion_dict:
                exclusion_dict[sample] = set()
            exclusion_dict[sample].add(contig_id)
    return exclusion_dict


def exclude_plasmid_contigs(plasmid_files, fasta_directory, output_fasta, exclusion_headers):
    """Exclude plasmid contigs and additional contigs from the exclusion list, and save remaining contigs. """
    excluded_headers = exclusion_headers  
    non_plasmid_headers = set()
    
    with open(output_fasta, "w") as output_handle:
        with open(fasta_directory, "r") as fasta_handle:
            for record in SeqIO.parse(fasta_handle, "fasta"):
                header = record.id.split()[0].replace('_R', '') 
                
                if header not in excluded_headers:
                    SeqIO.write(record, output_handle, "fasta")
                    non_plasmid_headers.add(header)
                else:
                    print(f"Excluded contig: {header}")
    
    print(f"Filtered contigs saved to {output_fasta}")
    return excluded_headers, non_plasmid_headers


def process_nonplasmid_contigs(plasmid_files, fasta_directory, output_directory, directory_path):

    exclusion_file = os.path.join(directory_path, "nested_contigs.txt")
    if not os.path.exists(exclusion_file):
        print(f"Exclusion file not found: {exclusion_file}. Skipping exclusion list processing.")
        exclusion_dict = {}
    else:
        exclusion_dict = load_exclusion_list(exclusion_file)

    os.makedirs(output_directory, exist_ok=True)

    for plasmid_file in os.listdir(plasmid_files):
        if plasmid_file.endswith('_replicon_filtered.txt'):
            sample_name = plasmid_file.replace('_replicon_filtered.txt', '').replace('_plasmid', '')

           
            plasmid_fasta = os.path.join(plasmid_files, plasmid_file)

            
            original_fasta = None
            for merged_file in os.listdir(fasta_directory):
                if merged_file.startswith(sample_name) and merged_file.endswith(".fasta"):
                    original_fasta = os.path.join(fasta_directory, merged_file)
                    break
            
           
            output_fasta = os.path.join(output_directory, f"{sample_name}_nonplasmid_contigs.fasta")
        

            if original_fasta and os.path.exists(original_fasta):
                exclusion_headers = exclusion_dict.get(sample_name, set())
                if os.path.exists(plasmid_fasta):
                    with open(plasmid_fasta, "r") as f:
                        for line in f:
                            if not line.startswith('qseqid'):
                                parts = line.strip().split('\t')
                                header = parts[0].replace('_R', '')
                                exclusion_headers.add(header)

                plasmid_headers, non_plasmid_headers = exclude_plasmid_contigs(
                    plasmid_fasta, original_fasta, output_fasta, exclusion_headers
                )

            else:
                print(f"Original FASTA file not found for sample: {sample_name}")


def save_headers_to_file(headers, output_file):
    """Save headers to a text file."""
    with open(output_file, "w") as f:
        for header in headers:
            f.write(f"{header}\n")
    print(f"Headers saved to {output_file}")


def process_files_for_col_contigs(file1, file2, output_rows):
    """Process files to identify Col and Inc plasmids and append new_sseqid to the replicon file."""
    col_contigs = {}
    with open(file1, 'r') as f1:
        for line in f1:
            parts = line.strip().split('\t')
            if len(parts) > 5:
                gene_name = parts[5].strip()
                contig_number = parts[1].strip()
                if gene_name.startswith("Col"):
                    col_contigs[gene_name] = contig_number

    base_name = os.path.basename(file2)
    sample_name = base_name.split('_plasmid_replicon_filtered.txt')[0]

    replicon_data = defaultdict(list)
    with open(file2, 'r') as f2:
        file2_rows = [line.strip().split('\t') for line in f2]
        if not file2_rows:
            return {} 
        header = file2_rows[0]
        for row in file2_rows[1:]:
            row += [''] * (len(header) - len(row))
            if len(row) < 2:
                continue 
            contig_number = row[0].strip()
            sseqid = row[1].strip()
            last_column = row[-1].strip() if row else ''
            cleaned_gene_name = re.sub(r'_pld\d+', '', last_column)
            cleaned_gene_name = re.sub(r'_[_]*\d+$', '', cleaned_gene_name)
            pld_match = re.search(r'_pld(\d+)', last_column)
            pld_number = int(pld_match.group(1)) if pld_match else 1
            replicon_data[(sseqid, last_column)].append({
                "row": row,
                "gene_name": cleaned_gene_name,
                "pld_number": pld_number,
                "last_column": last_column
            })

    new_sseqid_counter = {}
    x_sseqid_mapping = {}
    for (sseqid, last_column), data_list in replicon_data.items():
        gene_name = re.sub(r'_pld\d+', '', last_column)
        new_sseqid_parts = [re.sub(r'_[_]*1$', '', gene_name)]
        new_sseqid_base = f"{sample_name}_{'_'.join(new_sseqid_parts)}"
        if new_sseqid_base in new_sseqid_counter:
            new_sseqid_counter[new_sseqid_base] += 1
            new_sseqid_group = f"{new_sseqid_base}_{new_sseqid_counter[new_sseqid_base]}"
        else:
            new_sseqid_counter[new_sseqid_base] = 1
            new_sseqid_group = new_sseqid_base
        x_sseqid = f"{sample_name}_{sseqid}"
        x_sseqid_mapping[x_sseqid] = new_sseqid_group
        pld_match = re.search(r'_pld(\d+)', last_column)
        pld_number = int(pld_match.group(1)) if pld_match else 1
        new_sseqid_with_pld = new_sseqid_group
        if pld_number > 1:
            parts = new_sseqid_with_pld.split("_")
            if len(parts) > 1:
                try:
                    int(parts[-1])
                    new_sseqid_with_pld = "_".join(parts[:-1])
                except ValueError:
                    pass
            new_sseqid_with_pld = f"{new_sseqid_with_pld}_{pld_number}"
        elif pld_number == 1:
            parts = new_sseqid_with_pld.split("_")
            if len(parts) > 1:
                try:
                    int(parts[-1])
                    new_sseqid_with_pld = "_".join(parts[:-1])
                except ValueError:
                    pass
        for item in data_list:
            row = item["row"]
            if len(row) >= len(header) + 2:
                row[-2] = new_sseqid_with_pld
                row[-1] = x_sseqid
            else:
                row.append(new_sseqid_with_pld)
                row.append(x_sseqid)

    if "new_sseqid" not in header:
        header.append("new_sseqid")
    if "x_sseqid" not in header:
        header.append("x_sseqid")

    updated_rows = [header]
    for row in file2_rows[1:]:
        processed = False
        for data_list in replicon_data.values():
            for item in data_list:
                if row == item["row"]:
                    updated_rows.append(row)
                    processed = True
                    break
            if processed:
                break
        if not processed:
            row += ['', '']  
            updated_rows.append(row)


    with open(file2, 'w') as f2:
        for row in updated_rows:
            f2.write('\t'.join(row) + '\n')
    return x_sseqid_mapping

def rename_files_in_gene_directories(extracted_fasta_dir, replicon_dir):
    """Rename files in gene directories based on x_sseqid and new_sseqid mapping and update FASTA headers."""
    for gene_dir in os.listdir(extracted_fasta_dir):
        gene_dir_path = os.path.join(extracted_fasta_dir, gene_dir)
        if not os.path.isdir(gene_dir_path):
            continue
        files = os.listdir(gene_dir_path)
        for file in files:
            file_path = os.path.join(gene_dir_path, file)
            if not os.path.isfile(file_path):
                continue
            try:
                file_basename = os.path.basename(file_path)
                file_basename_without_ext, ext = os.path.splitext(file_basename)
                parts = file_basename_without_ext.split('_')
                if len(parts) < 2:
                    continue
                sample_name = parts[0]
                sseqid = '_'.join(parts[1:])
                x_sseqid = f"{sample_name}_{sseqid}"
                replicon_file = os.path.join(replicon_dir, f"{sample_name}_plasmid_replicon_filtered.txt")
                if not os.path.exists(replicon_file):
                    print(f"Error: Replicon file not found: {replicon_file}")
                    continue
                x_sseqid_mapping = {}
                with open(replicon_file, 'r') as f:
                    header = f.readline().strip().split('\t')
                    if "x_sseqid" not in header or "new_sseqid" not in header:
                        continue
                    x_idx = header.index("x_sseqid")
                    new_idx = header.index("new_sseqid")
                    for line in f:
                        line_parts = line.strip().split('\t')
                        if len(line_parts) <= max(x_idx, new_idx):
                            continue
                        x_sseqid_mapping[line_parts[x_idx]] = line_parts[new_idx]
                if x_sseqid not in x_sseqid_mapping:
                    continue
                new_name = x_sseqid_mapping[x_sseqid]
                new_file_path = os.path.join(gene_dir_path, f"{new_name}{ext}")
                shutil.move(file_path, new_file_path)
                with open(new_file_path, 'r') as infile:
                    lines = infile.readlines()
                with open(new_file_path, 'w') as outfile:
                    for line in lines:
                        if line.startswith('>'):
                            outfile.write(f">{new_name}\n")
                        else:
                            outfile.write(line)
            except Exception as e:
                continue

def rename_fasta_headers_in_replicon_dir(replicon_dir):
    """Rename headers in concatenated FASTA files using x_sseqid to new_sseqid mapping."""
    for fasta_file in os.listdir(replicon_dir):
        if not fasta_file.endswith('_plasmid_contigs.fasta'):
            continue
        sample_name = fasta_file.split('_')[0]
        replicon_file = os.path.join(replicon_dir, f"{sample_name}_plasmid_replicon_filtered.txt")
        if not os.path.exists(replicon_file):
            continue
        x_sseqid_mapping = {}
        with open(replicon_file, 'r') as f:
            header = f.readline().strip().split('\t')
            if "x_sseqid" not in header or "new_sseqid" not in header:
                continue
            x_idx = header.index("x_sseqid")
            new_idx = header.index("new_sseqid")
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < max(x_idx, new_idx) + 1:
                    continue
                x_sseqid = parts[x_idx]
                new_sseqid = parts[new_idx]
                x_sseqid_mapping[x_sseqid] = new_sseqid
        fasta_path = os.path.join(replicon_dir, fasta_file)
        temp_path = os.path.join(replicon_dir, f"TEMP_{fasta_file}")
        with open(fasta_path, 'r') as infile, open(temp_path, 'w') as outfile:
            for line in infile:
                if line.startswith('>'):
                    x_sseqid = line[1:].strip()
                    new_sseqid = x_sseqid_mapping.get(x_sseqid, x_sseqid)
                    outfile.write(f">{new_sseqid}\n")
                else:
                    outfile.write(line)
        os.replace(temp_path, fasta_path)

def process_directories_for_col_contigs(plasmid_dir, replicon_dir, combined_output_file, extracted_fasta_dir):
    """Process directories to identify Col and Inc plasmids and update replicon files."""
    output_rows = []
    plasmid_files = {os.path.basename(f): os.path.join(plasmid_dir, f) for f in os.listdir(plasmid_dir) if f.endswith('_plasmid.txt')}
    replicon_files = {os.path.basename(f): os.path.join(replicon_dir, f) for f in os.listdir(replicon_dir) if f.endswith('_plasmid_replicon_filtered.txt')}
    
    print(f"Found {len(plasmid_files)} plasmid files and {len(replicon_files)} replicon files.")  # Debug statement
    
    
    for file_name, plasmid_path in plasmid_files.items():
        base_sample_name = file_name.split('_')[0]
        expected_replicon_file = f"{base_sample_name}_plasmid_replicon_filtered.txt"
        replicon_path = replicon_files.get(expected_replicon_file)
        
        if replicon_path:
            print(f"Processing {plasmid_path} and {replicon_path}") 
            x_sseqid_mapping = process_files_for_col_contigs(plasmid_path, replicon_path, output_rows)
        else:
            print(f"No replicon file found for {file_name}. Expected: {expected_replicon_file}") 

    with open(combined_output_file, 'w') as out:
        for row in output_rows:
            out.write('\t'.join(row) + '\n')
    
    """Rename files in gene directories and update FASTA headers"""
    rename_files_in_gene_directories(extracted_fasta_dir, replicon_dir)
    rename_fasta_headers_in_replicon_dir(replicon_dir)   
 
def process_nonreplicon_contigs_with_circular_check(directory_path, fasta_directory):
    """Process non-replicon contigs and filter for circular contigs."""
    nonreplicon_dir = os.path.join(directory_path, "nonreplicon_contigs") 
    nonplasmid_files = os.path.join(directory_path, "nonplasmid_files")  
    output_dir = os.path.join(directory_path, "nonreplicon_contigs", "nonreplicon_processed")  
    plasmid_files = os.path.join(directory_path, "plasmid_files") 

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(plasmid_files, exist_ok=True)


    for file in os.listdir(nonreplicon_dir):
        if file.endswith(".txt"):

            basename = file.split("_")[0]
            fasta_file = None
            for f in os.listdir(fasta_directory):
                if f.startswith(basename) and f.endswith(".fasta") and basename in f:  
                    fasta_file = os.path.join(fasta_directory, f)
                    break

            if not fasta_file or not os.path.exists(fasta_file): 
                print(f"FASTA file for {basename} not found in {fasta_directory}, skipping...")
                continue
            nonplasmid_fasta_file = os.path.join(nonplasmid_files, f"{basename}_nonplasmid_contigs.fasta")
            combined_fasta_path = os.path.join(plasmid_files, f"{basename}_plasmid_contigs.fasta")

            if not os.path.exists(fasta_file):
                print(f"FASTA file for {basename} not found, skipping...")
                continue

            data_path = os.path.join(nonreplicon_dir, file)
            df = pd.read_csv(data_path, sep="\t")

            df["cleaned_qseqid"] = df["qseqid"].astype(str).str.replace(r"_R$", "", regex=True)

            """Filter conditions:"""
            df = df[(df["coverage_percentage"] == 100) & (df["qcovs"] == 100)]

            sseqid_groups = df.groupby("sseqid")["cleaned_qseqid"].nunique()
            single_contig_sseqids = sseqid_groups[sseqid_groups == 1].index
            df = df[df["sseqid"].isin(single_contig_sseqids)]

            """Parse the FASTA file and filter for circular contigs"""
            selected_contigs = []
            circular_qseqids = set()
            for record in SeqIO.parse(fasta_file, "fasta"):
                header = record.description
                contig_id = header.split()[0].lstrip(">")
                cleaned_contig_id = contig_id.split("_R")[0]  
                if cleaned_contig_id in df["cleaned_qseqid"].values and "circular=true" in header.lower():
                    selected_contigs.append(record)
                    circular_qseqids.add(cleaned_contig_id)  
                    print(f"Selected circular contig: {header}")  

            df = df[df["cleaned_qseqid"].isin(circular_qseqids)]

            """""Exclude chosen `cleaned_qseqid` from the nonplasmid FASTA file"""
            if os.path.exists(nonplasmid_fasta_file):
                nonplasmid_records = []
                for record in SeqIO.parse(nonplasmid_fasta_file, "fasta"):
                    contig_id = record.description.split()[0].lstrip(">")
                    cleaned_contig_id = contig_id.split("_R")[0]
                    if cleaned_contig_id not in circular_qseqids:
                        nonplasmid_records.append(record)

                updated_nonplasmid_path = os.path.join(nonplasmid_files, f"{basename}_nonplasmid_contigs.fasta")
                SeqIO.write(nonplasmid_records, updated_nonplasmid_path, "fasta")
                print(f"Updated nonplasmid FASTA saved for {basename}.")

            if df.empty or not selected_contigs:
                print(f"No matching circular contigs for {basename}, skipping...")
                continue

            output_data_path = os.path.join(output_dir, f"{basename}_filtered.txt")
            df.to_csv(output_data_path, sep="\t", index=False)

            for contig in selected_contigs:
                contig_id = contig.description.split()[0].lstrip(">")
                sseqid = df.loc[df["cleaned_qseqid"] == contig_id, "sseqid"].values[0]
                output_fasta_name = f"{basename}_{sseqid}.fasta"
                output_fasta_path = os.path.join(output_dir, output_fasta_name)

                contig.id = f"{basename}_{sseqid}"
                contig.description = ""  
                SeqIO.write([contig], output_fasta_path, "fasta")

                """Append the contig to the combined FASTA file in the plasmid_files directory"""
                with open(combined_fasta_path, "a") as combined_fasta:
                    SeqIO.write([contig], combined_fasta, "fasta")

            print(f"Processed {basename}: Contigs and filtered data saved in {output_dir} and combined FASTA in {plasmid_files}.")


def generate_sample_report(sample_name, replicon_contigs_file, fasta_file, output_report_file):
    """Generate a report file for a sample."""
    try:
        replicon_df = pd.read_csv(replicon_contigs_file, sep='\t')

        # Check if the 'new_sseqid' column 
        if 'new_sseqid' not in replicon_df.columns:
            print(f"Warning: 'new_sseqid' column not found in {replicon_contigs_file}. Using 'qseqid' instead.")
            if 'qseqid' not in replicon_df.columns:
                print(f"Error: Neither 'new_sseqid' nor 'qseqid' found in {replicon_contigs_file}. Cannot generate report.")
                return  # Exit if neither column exists
            replicon_df['new_sseqid'] = replicon_df['qseqid']  # Use qseqid as a temporary fix

        unique_new_sseqid_df = replicon_df[['new_sseqid', 'sseqid', 'coverage_percentage']].drop_duplicates()
        contig_lengths = {}
        with open(fasta_file, "r") as fasta_handle:
            for record in SeqIO.parse(fasta_handle, "fasta"):
                contig_lengths[record.id] = len(record.seq)

        # Handle cases where some contig IDs might not be in the FASTA file
        unique_new_sseqid_df['bases'] = unique_new_sseqid_df['new_sseqid'].map(contig_lengths)
        unique_new_sseqid_df = unique_new_sseqid_df.dropna(subset=['bases'])

        unique_new_sseqid_df.to_csv(output_report_file, sep='\t', index=False)
        print(f"Report generated successfully for {sample_name}") 
    except Exception as e:
        print(f"Error generating report for {sample_name}: {e}")

def process_samples_for_reports(directory_path, fasta_directory):
    """Process all samples in the directory to generate report files."""
    if not os.path.exists(directory_path) or not os.path.exists(fasta_directory):
        return
    for file in os.listdir(directory_path):
        if file.endswith("_replicon_filtered.txt"):
            sample_name = file.replace("_plasmid_replicon_filtered.txt", "")
            replicon_contigs_file = os.path.join(directory_path, file)
            fasta_file = os.path.join(fasta_directory, f"{sample_name}_plasmid_contigs.fasta")
            output_report_file = os.path.join(directory_path, f"{sample_name}_report.txt")
            if os.path.exists(fasta_file):
                generate_sample_report(sample_name, replicon_contigs_file, fasta_file, output_report_file)

def main():
    parser = argparse.ArgumentParser(description="Process plasmid and replicon data.")
    parser.add_argument('--dir', required=True, help="Directory containing the files and output directories.")

    args = parser.parse_args()

    directory_path = args.directory
    fasta_directory = os.path.join(directory_path, "unicycler_fasta")
    gene_list_file = os.path.join(directory_path, "plasmid_list.txt")

    #Process _PLSDB.txt files
    process_directory(directory_path)

    #Process files for overlapping contigs
    process_files_in_directory(directory_path)

    #Process filtered files
    process_directory_filtered(directory_path)

    #Process circular contigs
    process_circular_contigs(directory_path, directory_path)

    #Process files for coverage
    process_files_coverage(directory_path)

    #Process non-replicon contigs
    process_nonreplicon_contigs(directory_path, os.path.join(directory_path, "nonreplicon_contigs"))

    #Extract gene-sample mapping and create directories
    gene_directories = os.path.join(directory_path, "gene_directories")
    extract_gene_sample_mapping(directory_path, gene_directories)

    # Process gene directories
    gene_list = load_gene_list(gene_list_file)
    process_gene_directories(gene_directories, gene_list)

    # Process replicon filtered files
    process_replicon_filtered_files(gene_directories)

    #Concatenate sample files
    output_directory = os.path.join(directory_path, "plasmid_files")
    concatenate_sample_files(gene_directories, output_directory)

    #Extract and update contigs
    process_directories_for_extraction(gene_directories, directory_path)

    #Merge contigs
    process_directories_for_merging(gene_directories)

    #Copy merged fasta files
    destination_directory = os.path.join(directory_path, "extracted_fasta")
    copy_merged_fasta_files(gene_directories, destination_directory, gene_list)

    #Concatenate fasta files
    concatenate_fasta_files(gene_directories, output_directory)

    #Process non-plasmid contigs
    process_nonplasmid_contigs(output_directory, fasta_directory, os.path.join(directory_path, "nonplasmid_files"), directory_path)

    #Process Col and Inc plasmids
    combined_output_file = os.path.join(directory_path, "Processed_ColContigs.txt")
    process_directories_for_col_contigs(directory_path, output_directory, combined_output_file, destination_directory)

    #Process non-replicon contigs
    process_nonreplicon_contigs_with_circular_check(directory_path, fasta_directory)

    #Generate sample reports
    process_samples_for_reports(output_directory, output_directory)

if __name__ == "__main__":
    main()

