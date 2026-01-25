import os
import re
import shutil
import pandas as pd
import logging
from Bio import SeqIO, Entrez
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict
import argparse
import sys

def setup_logging(log_file='log.txt'):
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )

def filter_invalid_rows(file_path):
    logging.info(f"Filtering invalid rows in {file_path}")
    with open(file_path, 'r') as file:
        lines = file.readlines()

    header = lines[0]
    data_lines = lines[1:]

    filtered_lines = [line for line in data_lines if int(line.split('\t')[4]) <= int(line.split('\t')[5])]

    with open(file_path, 'w') as file:
        file.write(header)
        file.writelines(filtered_lines)

def extract_contig_lengths_from_fasta(fasta_file):
    logging.info(f"Extracting contig lengths from {fasta_file}")
    contig_lengths = {}
    with open(fasta_file, 'r') as fasta_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            contig_lengths[record.id] = len(record.seq)
    return contig_lengths

def calculate_query_percentage_aligned(row):
    return (int(row['qend']) - int(row['qstart']) + 1) / int(row['qlen'])

def process_plsdb_file(file_path, fasta_file):
    logging.info(f"Processing {file_path}")
    filter_invalid_rows(file_path)

    if not os.path.exists(fasta_file):
        logging.error(f"FASTA file {fasta_file} not found. Skipping {file_path}.")
        return

    logging.info(f"Extracting contig lengths from {fasta_file}")
    contig_lengths = extract_contig_lengths_from_fasta(fasta_file)

    plsdb_df = pd.read_csv(file_path, sep='\t', dtype=str)
    plsdb_df['qlen'] = plsdb_df['qseqid'].apply(lambda x: contig_lengths.get(x.strip(), 'N/A'))
    plsdb_df['q_perc_aligned'] = plsdb_df.apply(calculate_query_percentage_aligned, axis=1)
    plsdb_df.to_csv(file_path, sep='\t', index=False)
    logging.info(f"Updated {file_path} with qlen and q_perc_aligned columns.")

def process_directory(directory_path):
    logging.info(f"Processing _PLSDB.txt files in {directory_path}")
    for file in os.listdir(directory_path):
        if file.endswith("_PLSDB.txt"):
            file_path = os.path.join(directory_path, file)
            base_name = file.replace('_PLSDB.txt', '')
            fasta_file = os.path.join(directory_path, f"{base_name}_merged.fasta")
            process_plsdb_file(file_path, fasta_file)

def process_contigs(input_file, output_file, nested_contigs_log, unique_nested_contigs):
    logging.info(f"Processing contigs in {input_file}")
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

    logging.info(f"Updated filtered contigs saved to {output_file}")
    logging.info(f"Nested contigs logged in {nested_contigs_log}")

def process_files_in_directory(directory_path):
    logging.info(f"Processing files for overlapping contigs in {directory_path}")
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
    data['order'] = data.groupby('sseqid').cumcount()
    return data

def handle_non_overlapping_contigs(group):
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
    group = group.sort_values(by=['qseqid', 'qstart'])
    group, to_remove = handle_non_overlapping_contigs(group)
    group, to_remove = handle_overlapping_contigs(group, to_remove)
    return group.drop(to_remove)

def process_file(input_file, output_file):
    logging.info(f"Processing file {input_file}")
    data = pd.read_csv(input_file, sep='\t')
    data = add_order_column(data)
    processed_data = data.groupby('sseqid').apply(process_group).reset_index(drop=True)
    final_data = processed_data.sort_values(by=['sseqid', 'sstart']).drop(columns=['order'])
    final_data.to_csv(output_file, sep='\t', index=False)
    logging.info(f"Processed file saved to {output_file}")

def process_directory_filtered(input_directory):
    logging.info(f"Processing filtered files in {input_directory}")
    for filename in os.listdir(input_directory):
        if filename.endswith('_overlap.txt'):
            input_file = os.path.join(input_directory, filename)
            output_file = os.path.join(input_directory, filename.replace('_overlap.txt', '_overlap_filtered.txt'))
            process_file(input_file, output_file)

def get_circular_contigs(fasta_file):
    circular_contigs = set()
    for record in SeqIO.parse(fasta_file, "fasta"):
        if "circular=true" in record.description:
            circular_contigs.add(record.id)
    return circular_contigs

def standardize_qseqid(qseqid):
    return qseqid.rstrip('_R')

def process_circular_contigs(input_dir, output_dir):
    logging.info(f"Processing circular contigs from {input_dir} to {output_dir}")
    os.makedirs(output_dir, exist_ok=True)
    for file in os.listdir(input_dir):
        if file.endswith("_overlap_filtered.txt"):
            sample_name = file.replace("_overlap_filtered.txt", "")
            fasta_file = os.path.join(input_dir, f"{sample_name}_merged.fasta")
            if not os.path.exists(fasta_file):
                logging.warning(f"FASTA file not found for {sample_name}, skipping...")
                continue
            filtered_file_path = os.path.join(input_dir, file)
            df = pd.read_csv(filtered_file_path, sep="\t")
            circular_contigs = get_circular_contigs(fasta_file)
            if not circular_contigs:
                logging.info(f"No circular contigs found for {sample_name}. Saving unfiltered file.")
                output_file = os.path.join(output_dir, file)
                df.to_csv(output_file, sep="\t", index=False)
                continue
            logging.info(f"Circular contigs found for {sample_name}: {circular_contigs}")
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
            logging.info(f"Plasmids to exclude for {sample_name}: {plasmids_to_exclude}")
            filtered_df = df[~df['sseqid'].isin(plasmids_to_exclude)]
            output_file = os.path.join(output_dir, file)
            filtered_df.to_csv(output_file, sep="\t", index=False)
            logging.info(f"Processed {file}: Excluded {len(plasmids_to_exclude)} plasmids. Saved to {output_file}.")

def merge_intervals(intervals):
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    merged_intervals = []
    for interval in sorted_intervals:
        if not merged_intervals or merged_intervals[-1][1] < interval[0] - 1:
            merged_intervals.append(list(interval))
        else:
            merged_intervals[-1][1] = max(merged_intervals[-1][1], interval[1])
    return merged_intervals

def calculate_coverage_percentage(merged_intervals, subject_length):
    if subject_length == 0:
        return 0
    covered_length = sum(end - start + 1 for start, end in merged_intervals)
    return (covered_length / subject_length) * 100

def process_files_coverage(input_directory, coverage_cutoff=0.0):
    logging.info(f"Processing files for coverage in {input_directory}")
    for filename in sorted(os.listdir(input_directory)):
        if filename.endswith('_overlap_filtered.txt'):
            input_file = os.path.join(input_directory, filename)
            logging.info(f"Processing file: {filename}")
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
                        logging.warning(f"Skipping line due to ValueError: {ve}")
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
            logging.info(f"Processed and updated {input_file}.")

def extract_gene_sample_mapping(input_dir, gene_directories):
    logging.info(f"Extracting gene-sample mapping from {input_dir}")
    gene_samples = {}
    
    for filename in os.listdir(input_dir):
        if filename.endswith("_plasmid.txt"):
            sample_base = filename.replace('_plasmid.txt', '')
            sample_parts = sample_base.split('_')
            sample = sample_parts[0] if len(sample_parts) > 1 else sample_base
            
            file_path = os.path.join(input_dir, filename)
            try:
                with open(file_path, 'r') as file:
                    next(file)
                    for line in file:
                        parts = line.strip().split('\t')
                        if len(parts) > 5:
                            gene = parts[5].strip()
                            sanitized_gene = ''.join(char if char.isalnum() else '_' for char in gene)
                            if sanitized_gene not in gene_samples:
                                gene_samples[sanitized_gene] = set()
                            gene_samples[sanitized_gene].add(sample)
            except Exception as e:
                logging.error(f"Error processing file {filename}: {e}")
    
    os.makedirs(gene_directories, exist_ok=True)
    for gene, samples in gene_samples.items():
        gene_dir = os.path.join(gene_directories, gene)
        os.makedirs(gene_dir, exist_ok=True)
        logging.info(f"Created directory for gene: {gene}")
        
        for sample in samples:
            sample_dir = os.path.join(gene_dir, sample)
            os.makedirs(sample_dir, exist_ok=True)
            logging.info(f"Created sample directory: {sample_dir}")
            
            pattern_start = f"{sample}_"
            found_files = False
            
            for filename in os.listdir(input_dir):
                if filename.endswith("_plasmid.txt") and filename.startswith(pattern_start):
                    source_file = os.path.join(input_dir, filename)
                    dest_file = os.path.join(sample_dir, f"{sample}_plasmid.txt")
                    shutil.copy2(source_file, dest_file)
                    logging.info(f"Copied plasmid file: {source_file} to {dest_file}")
                    found_files = True
                
                elif filename.endswith("_overlap_filtered.txt"):
                    base_name = filename.replace('_overlap_filtered.txt', '')
                    if base_name.startswith(pattern_start):
                        source_file = os.path.join(input_dir, filename)
                        dest_file = os.path.join(sample_dir, f"{sample}_overlap_filtered.txt")
                        shutil.copy2(source_file, dest_file)
                        logging.info(f"Copied overlap file: {source_file} to {dest_file}")
                        found_files = True
            
            if not found_files:
                logging.warning(f"No matching files found for sample {sample} in {input_dir}")


def load_gene_list(file_path):
    logging.info(f"Loading gene list from {file_path}")
    with open(file_path, 'r') as file:
        gene_list = [line.strip() for line in file if line.strip()]
    return gene_list

def match_contig_id(contig_id, qseqid):
    return re.fullmatch(f"{re.escape(contig_id)}(_R)?", qseqid) is not None

def process_gene_directories(top_level_directory, gene_list):
    logging.info(f"Processing gene directories in {top_level_directory}")
    for gene_name in gene_list:
        gene_dir_path = os.path.join(top_level_directory, gene_name)
        if os.path.exists(gene_dir_path) and os.path.isdir(gene_dir_path):
            logging.info(f"Processing gene directory: {gene_name}")
            for root, dirs, files in os.walk(gene_dir_path):
                directory_name = os.path.basename(root)
                # Skip the gene directory itself, only process sample subdirectories
                if directory_name == gene_name:
                    continue
                    
                logging.info(f"Processing subdirectory: {directory_name}")
                abricate_file = next((f for f in files if f.endswith('_plasmid.txt')), None)
                overlap_file = next((f for f in files if f.endswith('_overlap_filtered.txt')), None)
                
                if not abricate_file or not overlap_file:
                    logging.warning(f"Required files not found in subdirectory {directory_name}. Looking for: {directory_name}_plasmid.txt and {directory_name}_overlap_filtered.txt")
                    # List what files are actually present for debugging
                    actual_files = [f for f in files if f.endswith(('.txt', '.fasta'))]
                    if actual_files:
                        logging.info(f"Files found in {directory_name}: {actual_files}")
                    continue

                logging.info(f"Abricate file: {abricate_file}, Overlap file: {overlap_file}")
                abricate_path = os.path.join(root, abricate_file)
                overlap_path = os.path.join(root, overlap_file)
                gene_output_path = os.path.join(root, f"{directory_name}_replicon.txt")
                all_plasmids_output_path = os.path.join(root, f"{directory_name}_replicon_contigs.txt")
                
                try:
                    abricate_df = pd.read_csv(abricate_path, sep='\t', dtype=str, comment=None)
                    # Clean up the GENE column
                    abricate_df['GENE'] = abricate_df['GENE'].str.strip()
                    
                    # Check if the target gene exists in this sample
                    contig_info = abricate_df[abricate_df['GENE'] == gene_name]
                    
                    if contig_info.empty:
                        logging.info(f"Gene {gene_name} not found in sample {directory_name}")
                        continue
                        
                    overlap_df = pd.read_csv(overlap_path, sep='\t', dtype=str, comment=None)
                    overlap_df['qseqid'] = overlap_df['qseqid'].str.strip()
                    
                    plasmids_found = False
                    for _, row in contig_info.iterrows():
                        contig_id = str(row['SEQUENCE']).strip()
                        logging.info(f"Looking for contig {contig_id} in overlap file for gene {gene_name}")
                        
                        plasmids_on_contig = overlap_df[overlap_df['qseqid'].apply(lambda x: match_contig_id(contig_id, x))]
                        
                        if not plasmids_on_contig.empty:
                            plasmids_found = True
                            logging.info(f"Found {len(plasmids_on_contig)} plasmids on contig {contig_id} for gene {gene_name}")
                            
                            if os.path.exists(gene_output_path):
                                plasmids_on_contig.to_csv(gene_output_path, sep='\t', index=False, mode='a', header=False)
                            else:
                                plasmids_on_contig.to_csv(gene_output_path, sep='\t', index=False)
                        else:
                            logging.info(f"No plasmids found on contig {contig_id} for gene {gene_name} in sample {directory_name}")

                    if plasmids_found:
                        if os.path.exists(gene_output_path):
                            gene_df = pd.read_csv(gene_output_path, sep='\t', dtype=str)
                            plasmid_ids = gene_df['sseqid'].unique()
                            all_plasmids_df = overlap_df[overlap_df['sseqid'].isin(plasmid_ids)]
                            all_plasmids_df['gene_name'] = gene_name
                            all_plasmids_df.to_csv(all_plasmids_output_path, sep='\t', index=False)
                            logging.info(f"All occurrences of the identified plasmids have been saved to {all_plasmids_output_path} with gene name.")
                        
                        if os.path.exists(gene_output_path):
                            os.remove(gene_output_path)
                            logging.info(f"Deleted the gene output file: {gene_output_path}")
                    else:
                        logging.info(f"No plasmids were identified for gene {gene_name} in sample {directory_name}")
                        
                except Exception as e:
                    logging.error(f"An error occurred while processing {directory_name}: {e}")
                    import traceback
                    logging.error(traceback.format_exc())
        else:
            logging.warning(f"Gene directory not found: {gene_name}")

def process_replicon_filtered_files(top_level_directory):
    logging.info(f"Processing replicon filtered files in {top_level_directory}")
    for root, dirs, files in os.walk(top_level_directory):
        directory_name = os.path.basename(root)
        logging.info(f"Processing directory: {directory_name}")
        for file_name in files:
            if file_name.endswith('_replicon_contigs.txt'):
                file_path = os.path.join(root, file_name)
                logging.info(f"Opening file: {file_path}")
                try:
                    df = pd.read_csv(file_path, delimiter='\t')
                    required_columns = ['gene_name', 'coverage_percentage', 'pident', 'bases_covered', 'sseqid', 'sstart']
                    if not all(col in df.columns for col in required_columns):
                        logging.warning(f"Missing required columns in file {file_path}. Skipping.")
                        continue
                    df['coverage_percentage'] = pd.to_numeric(df['coverage_percentage'].str.replace('%', ''), errors='coerce')
                    df['pident'] = pd.to_numeric(df['pident'], errors='coerce')
                    df['bases_covered'] = pd.to_numeric(df['bases_covered'], errors='coerce')
                    df = df.dropna(subset=['coverage_percentage', 'pident', 'bases_covered'])
                    if df.empty:
                        logging.warning(f"File {file_path} is empty after filtering NaN values. Skipping.")
                        continue
                    gene_name = df['gene_name'].values[0]
                    if 'Col' in gene_name or 'rep_cluster' in gene_name:
                        logging.info(f"Gene {gene_name} starts with 'Col', selecting plasmid based on highest coverage_percentage.")
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
                        logging.info(f"Gene {gene_name} does not start with 'Col'. Applying progressive coverage filter.")
                        thresholds = [80, 70, 60, 50, 40, 30, 20, 10]
                        for threshold in thresholds:
                            filtered_df = df[df['coverage_percentage'] >= threshold]
                            if not filtered_df.empty:
                                logging.info(f"Found plasmids with coverage_percentage >= {threshold}.")
                                break
                        else:
                            logging.warning(f"No plasmids found with coverage_percentage >= 20. Skipping file {file_path}.")
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
                    logging.info(f"Number of rows after filtering for plasmid '{best_plasmid_id}': {len(filtered_df)}")
                    sorted_df = filtered_df.sort_values(by='sstart')
                    sorted_file_name = file_name.replace('_replicon_contigs.txt', '_replicon_filtered.txt')
                    sorted_file_path = os.path.join(root, sorted_file_name)
                    sorted_df.to_csv(sorted_file_path, sep='\t', index=False)
                    logging.info(f"Processed and saved: {sorted_file_path}")
                except Exception as e:
                    logging.error(f"Error processing file {file_path}: {e}")

def concatenate_sample_files(top_level_directory, output_directory):
    logging.info(f"Concatenating sample files from {top_level_directory} to {output_directory}")
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
        logging.info(f"Concatenated file saved: {output_file}")

def trim_contig(sequence, overlap):
    return sequence[:-overlap] if overlap > 0 else sequence

def extract_and_update_contigs(input_txt, fasta_file, output_fasta):
    logging.info(f"Extracting and updating contigs from {input_txt}")
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
                logging.warning(f"Contig {contig_id} not found in FASTA file.")
        if len(contigs_to_write) == 1:
            contig_id, trimmed_sequence = contigs_to_write[0]
            output_handle.write(f">{contig_id}\n{trimmed_sequence}\n")
        elif len(contigs_to_write) > 1:
            first_contig_id = contigs_to_write[0][0]
            first_contig_qstart = filtered_df.loc[filtered_df['qseqid'] == first_contig_id, 'qstart'].values[0]
            last_contig = contigs_to_write[-1]
            last_contig_id = last_contig[0]
            last_contig_qend = filtered_df.loc[filtered_df['qseqid'] == last_contig_id, 'qend'].values[-1]
            logging.info(f"First contig qstart: {first_contig_qstart}, Last contig qend: {last_contig_qend}")
            if last_contig_qend < first_contig_qstart:
                logging.info(f"Removing last contig due to overlap: {last_contig_id}")
                contigs_to_write.pop()
            for contig_id, trimmed_sequence in contigs_to_write:
                output_handle.write(f">{contig_id}\n{trimmed_sequence}\n")
    logging.info(f"Extracted contigs saved to {output_fasta}")

def process_directories_for_extraction(top_level_directory, fasta_directory):
    logging.info(f"Processing directories for extraction from {top_level_directory}")
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
                    logging.warning(f"FASTA file not found for {basename} in {fasta_directory}")

def merge_contigs(input_fasta, output_fasta, sample_name, plasmid_name):
    logging.info(f"Merging contigs for {sample_name}, plasmid {plasmid_name}")
    merged_sequence = ""
    with open(input_fasta, "r") as fasta_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            merged_sequence += str(record.seq)
    new_header = f">{sample_name}_{plasmid_name}"
    with open(output_fasta, "w") as output_handle:
        output_handle.write(f"{new_header}\n{merged_sequence}\n")
    logging.info(f"Merged sequence saved to {output_fasta}")

def process_directories_for_merging(top_level_directory):
    logging.info(f"Processing directories for merging in {top_level_directory}")
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
                        logging.warning(f"'sseqid' column not found or file is empty: {overlapped_file}")
                        continue
                    sample_name = directory_name
                    output_fasta = os.path.join(root, f"{sample_name}_{plasmid_name}.fasta")
                    merge_contigs(input_fasta, output_fasta, sample_name, plasmid_name)
                else:
                    logging.warning(f"No overlapped file found in {root}")

def load_gene_list_from_file(gene_list_file):
    logging.info(f"Loading gene list from {gene_list_file}")
    with open(gene_list_file, 'r') as f:
        gene_list = [line.strip() for line in f if line.strip()]
    return gene_list

def copy_merged_fasta_files(top_level_directory, destination_directory, gene_list):
    logging.info(f"Copying merged fasta files to {destination_directory}")
    # Ensure destination directory exists
    os.makedirs(destination_directory, exist_ok=True)
    
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
                    logging.info(f"{parent_directory} is not in the gene list, skipping.")

def concatenate_fasta_files(top_level_directory, output_directory):
    logging.info(f"Concatenating fasta files from {top_level_directory} to {output_directory}")
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
        logging.info(f"Concatenated FASTA file saved: {output_file}")

def load_exclusion_list(exclusion_file):
    logging.info(f"Loading exclusion list from {exclusion_file}")
    exclusion_dict = {}
    with open(exclusion_file, "r") as f:
        for line in f:
            sample, contig_id = line.strip().split('\t')
            if sample not in exclusion_dict:
                exclusion_dict[sample] = set()
            exclusion_dict[sample].add(contig_id)
    return exclusion_dict

def exclude_plasmid_contigs(plasmid_files, fasta_file, output_fasta, exclusion_headers):
    logging.info(f"Excluding plasmid contigs from {fasta_file}")
    excluded_headers = exclusion_headers  
    non_plasmid_headers = set()
    
    with open(output_fasta, "w") as output_handle:
        with open(fasta_file, "r") as fasta_handle:
            for record in SeqIO.parse(fasta_handle, "fasta"):
                header = record.id.split()[0].replace('_R', '') 
                
                if header not in excluded_headers:
                    SeqIO.write(record, output_handle, "fasta")
                    non_plasmid_headers.add(header)
                else:
                    logging.info(f"Excluded contig: {header}")
    
    logging.info(f"Filtered contigs saved to {output_fasta}")
    return excluded_headers, non_plasmid_headers

def process_nonplasmid_contigs(plasmid_files, fasta_directory, output_directory, directory_path):
    logging.info(f"Processing non-plasmid contigs")
    
    exclusion_file = os.path.join(directory_path, "nested_contigs.txt")
    exclusion_dict = load_exclusion_list(exclusion_file) if os.path.exists(exclusion_file) else {}
    
    os.makedirs(output_directory, exist_ok=True)

    for plasmid_file in os.listdir(plasmid_files):
        if not plasmid_file.endswith('_replicon_filtered.txt'):
            continue
            
        sample_name = plasmid_file.replace('_plasmid_replicon_filtered.txt', '')
        plasmid_path = os.path.join(plasmid_files, plasmid_file)
        output_fasta = os.path.join(output_directory, f"{sample_name}_nonplasmid_contigs.fasta")
        
        original_fasta = None
        for fasta_file in os.listdir(fasta_directory):
            if not fasta_file.endswith(".fasta"):
                continue
                
            base_name = os.path.splitext(fasta_file)[0]
            
            if base_name == sample_name:
                original_fasta = os.path.join(fasta_directory, fasta_file)
                break
            
            if '_' in base_name:
                prefix = base_name.split('_')[0]
                if prefix == sample_name:
                    original_fasta = os.path.join(fasta_directory, fasta_file)
                    break
        
        if not original_fasta:
            logging.warning(f"FASTA file not found for {sample_name}")
            continue
        
        exclusion_headers = exclusion_dict.get(sample_name, set())
        
        if os.path.exists(plasmid_path):
            with open(plasmid_path, "r") as f:
                next(f)
                for line in f:
                    parts = line.strip().split('\t')
                    if parts:
                        header = parts[0].replace('_R', '')
                        exclusion_headers.add(header)
        
        exclude_plasmid_contigs(plasmid_path, original_fasta, output_fasta, exclusion_headers)

def save_headers_to_file(headers, output_file):
    with open(output_file, "w") as f:
        for header in headers:
            f.write(f"{header}\n")
    logging.info(f"Headers saved to {output_file}")

def process_files_for_col_contigs(file1, file2, output_rows):
    logging.info(f"Processing files for Col contigs: {file1} and {file2}")
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
    logging.info(f"Renaming files in gene directories from {extracted_fasta_dir}")
    
    # Check if directory exists, if not create it
    if not os.path.exists(extracted_fasta_dir):
        logging.warning(f"Extracted FASTA directory {extracted_fasta_dir} does not exist. Creating it.")
        os.makedirs(extracted_fasta_dir, exist_ok=True)
        return  # No files to process if directory was just created
    
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
                    logging.error(f"Error: Replicon file not found: {replicon_file}")
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
                logging.error(f"Error renaming file {file_path}: {e}")

def rename_fasta_headers_in_replicon_dir(replicon_dir):
    logging.info(f"Renaming FASTA headers in {replicon_dir}")
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
    logging.info(f"Processing directories for Col contigs from {plasmid_dir}")
    
    # Ensure extracted_fasta_dir exists
    os.makedirs(extracted_fasta_dir, exist_ok=True)
    
    output_rows = []
    plasmid_files = {os.path.basename(f): os.path.join(plasmid_dir, f) for f in os.listdir(plasmid_dir) if f.endswith('_plasmid.txt')}
    replicon_files = {os.path.basename(f): os.path.join(replicon_dir, f) for f in os.listdir(replicon_dir) if f.endswith('_plasmid_replicon_filtered.txt')}
    
    logging.info(f"Found {len(plasmid_files)} plasmid files and {len(replicon_files)} replicon files.")
    
    for file_name, plasmid_path in plasmid_files.items():
        base_sample_name = file_name.split('_')[0]
        expected_replicon_file = f"{base_sample_name}_plasmid_replicon_filtered.txt"
        replicon_path = replicon_files.get(expected_replicon_file)
        
        if replicon_path:
            logging.info(f"Processing {plasmid_path} and {replicon_path}") 
            x_sseqid_mapping = process_files_for_col_contigs(plasmid_path, replicon_path, output_rows)
        else:
            logging.warning(f"No replicon file found for {file_name}. Expected: {expected_replicon_file}") 

    with open(combined_output_file, 'w') as out:
        for row in output_rows:
            out.write('\t'.join(row) + '\n')
    
    logging.info("Renaming files in gene directories and updating FASTA headers")
    rename_files_in_gene_directories(extracted_fasta_dir, replicon_dir)
    rename_fasta_headers_in_replicon_dir(replicon_dir)   
 
def process_all_circular_contigs(directory_path, fasta_directory):
    logging.info(f"Processing all circular contigs from {directory_path}")
    output_dir = os.path.join(directory_path, "plasmid_files")
    os.makedirs(output_dir, exist_ok=True)

    for file in os.listdir(directory_path):
        if file.endswith("_overlap_filtered.txt"):
            sample_name = file.replace("_overlap_filtered.txt", "")
            
            # Find the corresponding Abricate plasmid file
            plasmid_file = None
            for f in os.listdir(directory_path):
                if f.startswith(sample_name) and f.endswith("_plasmid.txt"):
                    plasmid_file = os.path.join(directory_path, f)
                    break
            
            # Load replicon information from Abricate file
            contigs_with_replicons = set()
            if plasmid_file and os.path.exists(plasmid_file):
                try:
                    plasmid_df = pd.read_csv(plasmid_file, sep='\t')
                    if 'SEQUENCE' in plasmid_df.columns:
                        contigs_with_replicons = set(plasmid_df['SEQUENCE'].astype(str).str.strip().unique())
                        logging.info(f"Found {len(contigs_with_replicons)} contigs with replicons in {sample_name}: {contigs_with_replicons}")
                except Exception as e:
                    logging.warning(f"Error reading plasmid file {plasmid_file}: {e}")
            
            fasta_file = None
            for f in os.listdir(fasta_directory):
                if f.startswith(sample_name) and f.endswith(".fasta"):
                    fasta_file = os.path.join(fasta_directory, f)
                    break
            
            if not fasta_file or not os.path.exists(fasta_file):
                logging.warning(f"FASTA file for {sample_name} not found in {fasta_directory}, skipping...")
                continue
            
            circular_contigs = {}
            for record in SeqIO.parse(fasta_file, "fasta"):
                if "circular=true" in record.description:
                    contig_id = record.id.split()[0]
                    # Check if this contig has replicons in the Abricate file
                    if contig_id not in contigs_with_replicons:
                        circular_contigs[record.id] = {
                            'number': contig_id,
                            'length': len(record.seq)
                        }
                    else:
                        logging.info(f"Skipping circular contig {contig_id} as it has replicons identified by Abricate")
            
            if not circular_contigs:
                logging.info(f"No circular contigs found for {sample_name} (after filtering contigs with replicons)")
                continue
            
            logging.info(f"Found {len(circular_contigs)} circular contigs without replicons for {sample_name}: {list(circular_contigs.keys())}")
            
            # Create entries for circular contigs without alignment data
            circular_entries = []
            for contig_id, contig_info in circular_contigs.items():
                new_sseqid = f"{sample_name}_circular_{contig_info['number']}"
                # Create a minimal entry for the circular contig
                entry = {
                    'qseqid': contig_id,
                    'sseqid': '-',  # No specific plasmid reference
                    'qstart': '1',
                    'qend': str(contig_info['length']),
                    'sstart': '-',
                    'send': '-',
                    'evalue': '-',
                    'bitscore': '-',
                    'pident': '-',
                    'qcovs': '100',  # 100% coverage since it's the whole contig
                    'slen': '-',
                    'qlen': str(contig_info['length']),
                    'q_perc_aligned': '1.0',  # 100% aligned to itself
                    'overlapped_bases': '0',
                    'coverage_percentage': '100.00%',  # 100% coverage
                    'bases_covered': str(contig_info['length']),
                    'new_sseqid': new_sseqid,
                    'x_sseqid': contig_id,
                    'gene_name': 'circular'
                }
                circular_entries.append(entry)
            
            # Convert to DataFrame
            circular_df = pd.DataFrame(circular_entries)
            
            replicon_file = os.path.join(output_dir, f"{sample_name}_plasmid_replicon_filtered.txt")
            
            if os.path.exists(replicon_file):
                existing_df = pd.read_csv(replicon_file, sep="\t")
                # Remove any existing entries for the same circular contigs to avoid duplicates
                existing_circular_mask = existing_df['x_sseqid'].isin(circular_contigs.keys())
                existing_df = existing_df[~existing_circular_mask]
                combined_df = pd.concat([existing_df, circular_df], ignore_index=True)
                combined_df.to_csv(replicon_file, sep="\t", index=False)
            else:
                circular_df.to_csv(replicon_file, sep="\t", index=False)
            
            logging.info(f"Added {len(circular_df)} circular contigs (without replicons) to plasmid file for {sample_name}")
            
            # Also add to the combined FASTA file
            combined_fasta_path = os.path.join(output_dir, f"{sample_name}_plasmid_contigs.fasta")
            
            circular_records = []
            for record in SeqIO.parse(fasta_file, "fasta"):
                if record.id in circular_contigs:
                    new_id = f"{sample_name}_circular_{circular_contigs[record.id]['number']}"
                    new_record = record
                    new_record.id = new_id
                    new_record.description = new_id
                    circular_records.append(new_record)
            
            if os.path.exists(combined_fasta_path):
                # Read existing records and remove any with the same new_sseqid to avoid duplicates
                existing_records = list(SeqIO.parse(combined_fasta_path, "fasta"))
                existing_ids = {rec.id for rec in existing_records}
                
                # Filter out circular records that already exist
                new_circular_records = [rec for rec in circular_records if rec.id not in existing_ids]
                
                if new_circular_records:
                    with open(combined_fasta_path, "a") as combined_fasta:
                        SeqIO.write(new_circular_records, combined_fasta, "fasta")
                    logging.info(f"Added {len(new_circular_records)} new circular contigs (without replicons) to FASTA file for {sample_name}")
                else:
                    logging.info(f"No new circular contigs to add to FASTA file for {sample_name}")
            else:
                SeqIO.write(circular_records, combined_fasta_path, "fasta")
                logging.info(f"Created new FASTA file with {len(circular_records)} circular contigs (without replicons) for {sample_name}")

def generate_sample_report(output_directory, sample_name_or_directory):
    if os.path.isdir(sample_name_or_directory):
        directory = sample_name_or_directory
        logging.info(f"Generating sample reports for all samples in {directory}")
        for filename in os.listdir(directory):
            if filename.endswith('_plasmid_replicon_filtered.txt'):
                sample_name = filename.replace('_plasmid_replicon_filtered.txt', '')
                generate_sample_report(directory, sample_name)
        return
    
    sample_name = sample_name_or_directory
    try:
        logging.info(f"Generating sample report for {sample_name}")
        replicon_contigs_file = os.path.join(output_directory, f"{sample_name}_plasmid_replicon_filtered.txt")
        fasta_file = os.path.join(output_directory, f"{sample_name}_plasmid_contigs.fasta")
        output_report_file = os.path.join(output_directory, f"{sample_name}_report.tsv")
        
        if not os.path.exists(replicon_contigs_file):
            logging.warning(f"Replicon contigs file not found: {replicon_contigs_file}")
            return
        
        replicon_df = pd.read_csv(replicon_contigs_file, sep='\t')

        # Handle circular contigs specially
        if 'new_sseqid' not in replicon_df.columns:
            if 'qseqid' not in replicon_df.columns:
                return
            replicon_df['new_sseqid'] = replicon_df['qseqid']

        # For circular contigs, we want one entry per circular contig, not per alignment
        circular_mask = replicon_df['gene_name'] == 'circular'
        
        # For circular contigs, take the first occurrence (they should all be the same for a given contig)
        circular_df = replicon_df[circular_mask].drop_duplicates('new_sseqid', keep='first')
        
        # For non-circular contigs, use the existing logic
        non_circular_df = replicon_df[~circular_mask]
        
        if not non_circular_df.empty:
            unique_non_circular_df = non_circular_df[['new_sseqid', 'sseqid', 'coverage_percentage']].drop_duplicates()
        else:
            unique_non_circular_df = pd.DataFrame(columns=['new_sseqid', 'sseqid', 'coverage_percentage'])
        
        # Combine circular and non-circular entries
        if not circular_df.empty:
            circular_report_df = circular_df[['new_sseqid', 'sseqid', 'coverage_percentage']].copy()
            combined_df = pd.concat([unique_non_circular_df, circular_report_df], ignore_index=True)
        else:
            combined_df = unique_non_circular_df
        
        # Get contig lengths from FASTA file
        contig_lengths = {}
        if os.path.exists(fasta_file):
            with open(fasta_file, "r") as fasta_handle:
                for record in SeqIO.parse(fasta_handle, "fasta"):
                    contig_lengths[record.id] = len(record.seq)
        else:
            logging.warning(f"FASTA file not found: {fasta_file}")

        combined_df['bases'] = combined_df['new_sseqid'].map(contig_lengths)
        combined_df = combined_df.dropna(subset=['bases'])

        combined_df['type'] = combined_df['new_sseqid'].apply(
            lambda x: 'circular' if 'circular' in x else 'replicon'
        )

        combined_df.to_csv(output_report_file, sep='\t', index=False)
        logging.info(f"Report generated successfully for {sample_name} with {len(combined_df)} entries") 
    except Exception as e:
        logging.error(f"Error generating report for {sample_name}: {e}")

def main():
    # Set up logging in the input directory
    directory_path = None
    if len(sys.argv) > 1 and '--dir' in sys.argv:
        dir_index = sys.argv.index('--dir') + 1
        if dir_index < len(sys.argv):
            directory_path = sys.argv[dir_index]
    
    if directory_path:
        log_file = os.path.join(directory_path, 'log.txt')
    else:
        log_file = 'log.txt'
    
    setup_logging(log_file)
    logging.info("Starting plasmid processing pipeline")
    
    parser = argparse.ArgumentParser(description="Process plasmid and replicon data.")
    parser.add_argument('--dir', required=True, help="Directory containing the files and output directories.")

    args = parser.parse_args()

    directory_path = args.dir
    fasta_directory = os.path.join(directory_path, "unicycler_fasta")
    gene_list_file = os.path.join(directory_path, "plasmid_list.txt")

    logging.info(f"Working directory: {directory_path}")
    logging.info(f"FASTA directory: {fasta_directory}")
    logging.info(f"Gene list file: {gene_list_file}")

    process_directory(directory_path)
    process_files_in_directory(directory_path)
    process_directory_filtered(directory_path)
    process_circular_contigs(directory_path, directory_path)
    process_files_coverage(directory_path)

    gene_directories = os.path.join(directory_path, "gene_directories")
    extract_gene_sample_mapping(directory_path, gene_directories)

    gene_list = load_gene_list(gene_list_file)
    process_gene_directories(gene_directories, gene_list)
    process_replicon_filtered_files(gene_directories)

    output_directory = os.path.join(directory_path, "plasmid_files")
    concatenate_sample_files(gene_directories, output_directory)
    process_directories_for_extraction(gene_directories, directory_path)
    process_directories_for_merging(gene_directories)

    destination_directory = os.path.join(directory_path, "extracted_fasta")
    copy_merged_fasta_files(gene_directories, destination_directory, gene_list)
    concatenate_fasta_files(gene_directories, output_directory)

    process_nonplasmid_contigs(output_directory, fasta_directory, os.path.join(directory_path, "nonplasmid_files"), directory_path)

    combined_output_file = os.path.join(directory_path, "Processed_ColContigs.txt")
    process_directories_for_col_contigs(directory_path, output_directory, combined_output_file, destination_directory)

    process_all_circular_contigs(directory_path, fasta_directory)

    generate_sample_report(output_directory, output_directory)
    
    logging.info("Plasmid processing pipeline completed successfully")

if __name__ == "__main__":
    main()
