#!/usr/bin/env python3

import pandas as pd
import os
import re
import sys
import subprocess
import logging
from collections import defaultdict

def normalize_gene_name(gene_name: str) -> str:

    normalized = re.sub(r'_\d+$', '', gene_name)
    return normalized

def extract_gene_names(gff_file: str) -> list:

    gene_names = []
    
    try:
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 8 and fields[2] == 'CDS':
                    attributes = fields[8]
                    
                    # ONLY process if gene= pattern exists
                    gene_match = re.search(r'gene=([^;]+)', attributes)
                    if gene_match:
                        gene_name = gene_match.group(1)
                        normalized_name = normalize_gene_name(gene_name)
                        gene_names.append(normalized_name)
                    # Ignore lines without gene= attribute
    except Exception as e:
        print(f"Error reading GFF file {gff_file}: {e}")
        return []
    
    print(f"Extracted {len(gene_names)} genes with 'gene=' attribute from {gff_file}")
    return gene_names

def get_shared_and_missing_genes(genes_ref: list, genes_recon: list) -> tuple:

    set_ref = set(genes_ref)
    set_recon = set(genes_recon)
    
    shared_genes = set_ref & set_recon
    missing_in_recon = set_ref - set_recon
    missing_in_ref = set_recon - set_ref
    
    return shared_genes, missing_in_recon, missing_in_ref

def create_circular_pairs(genes: list) -> set:

    pairs = set()
    
    if len(genes) < 2:
        return pairs
    
    # Create linear pairs
    for i in range(len(genes) - 1):
        gene1, gene2 = genes[i], genes[i + 1]
        # Add both directions for bidirectional matching
        pairs.add((gene1, gene2))
        pairs.add((gene2, gene1))
    
    # Add circular pair (last, first) in both directions
    last_gene, first_gene = genes[-1], genes[0]
    pairs.add((last_gene, first_gene))
    pairs.add((first_gene, last_gene))
    
    return pairs

def compute_ngoc_with_circularity(genes_ref: list, genes_recon: list) -> tuple:

    # Get shared genes
    shared_genes, missing_in_recon, missing_in_ref = get_shared_and_missing_genes(genes_ref, genes_recon)
    total_shared = len(shared_genes)
    
    print(f"Shared genes between reference and reconstruction: {total_shared}")
    
    # Create pairs for reference (include ALL genes, not just shared ones)
    ref_pairs = create_circular_pairs(genes_ref)
    
    # Create pairs for reconstruction (include ALL genes)
    recon_pairs = create_circular_pairs(genes_recon)
    
    # Find common pairs (present in both)
    common_pairs = ref_pairs & recon_pairs
    
    # Calculate NGOC - ratio of common pairs to total possible pairs in reference
    total_possible_pairs_ref = len(ref_pairs)
    ngoc = len(common_pairs) / total_possible_pairs_ref if total_possible_pairs_ref > 0 else 0.0
    
    # Create detailed pair analysis - include ALL reference pairs
    pair_details = []
    for pair in sorted(ref_pairs):
        # Check if both genes in the pair are present in reconstruction
        gene1_missing = pair[0] in missing_in_recon
        gene2_missing = pair[1] in missing_in_recon
        
        if gene1_missing or gene2_missing:
            # Mark as BROKEN if either gene is missing
            missing_genes = []
            if gene1_missing:
                missing_genes.append(pair[0])
            if gene2_missing:
                missing_genes.append(pair[1])
            status = f"BROKEN (missing: {', '.join(missing_genes)})"
        elif pair in common_pairs:
            status = "CONCORDANT"
        else:
            status = "BROKEN (order disrupted)"
        
        pair_details.append((pair[0], pair[1], status))
    
    return round(ngoc, 3), shared_genes, pair_details, missing_in_recon, missing_in_ref

def run_synteny_analysis(ref_gff: str, query_gff: str, tool_name: str) -> dict:

    print(f"=== PROPER Synteny Analysis (Reference vs {tool_name}) ===\n")
    
    # Extract gene names (only those with gene= attribute)
    print("Extracting gene names from GFF files...")
    ref_genes = extract_gene_names(ref_gff)
    query_genes = extract_gene_names(query_gff)
    
    if not ref_genes:
        print(f"ERROR: No genes found in reference GFF: {ref_gff}")
        return {'NGOC': 0.0, 'Status': 'Error: No genes in reference'}
    
    if not query_genes:
        print(f"ERROR: No genes found in query GFF: {query_gff}")
        return {'NGOC': 0.0, 'Status': 'Error: No genes in query'}
    
    print(f"Reference genes: {len(ref_genes)}")
    print(f"{tool_name} genes: {len(query_genes)}")
    
    # Compute NGOC with circularity
    ngoc, shared_genes, pair_details, missing_in_recon, missing_in_ref = compute_ngoc_with_circularity(ref_genes, query_genes)
    
    # Display results
    print(f"\n=== RESULTS ===")
    print(f"NGOC Score: {ngoc:.3f}")
    print(f"Total Reference Genes: {len(ref_genes)}")
    print(f"Total {tool_name} Genes: {len(query_genes)}")
    print(f"Shared Genes: {len(shared_genes)}")
    
    # Show missing genes
    if missing_in_recon:
        print(f"\nGenes in Reference but MISSING in Reconstruction: {len(missing_in_recon)}")
        print(f"  - {list(missing_in_recon)}")
    
    if missing_in_ref:
        print(f"Genes in Reconstruction but MISSING in Reference: {len(missing_in_ref)}")
        print(f"  - {list(missing_in_ref)}")
    
    # Show pair analysis summary
    if pair_details:
        concordant_count = sum(1 for _, _, status in pair_details if status == "CONCORDANT")
        broken_count = sum(1 for _, _, status in pair_details if "BROKEN" in status)
        
        print(f"\nPair Analysis Summary:")
        print(f"Concordant pairs: {concordant_count}")
        print(f"Broken pairs: {broken_count}")
        print(f"Total pairs considered: {len(pair_details)}")
    
    # Interpretation
    print(f"\n=== INTERPRETATION ===")
    if ngoc >= 0.9:
        interpretation = "EXCELLENT gene order conservation"
        print(f"✅ {interpretation}")
    elif ngoc >= 0.7:
        interpretation = "GOOD gene order conservation"
        print(f"✅ {interpretation}")
    elif ngoc >= 0.5:
        interpretation = "MODERATE gene order conservation"
        print(f"⚠️ {interpretation}")
    else:
        interpretation = "POOR gene order conservation"
        print(f"❌ {interpretation}")
    
    # Prepare results
    results = {
        'NGOC': ngoc,
        'Shared_Genes': len(shared_genes),
        'Total_Ref_Genes': len(ref_genes),
        f'Total_{tool_name}_Genes': len(query_genes),
        'Missing_in_Reconstruction': len(missing_in_recon),
        'Missing_in_Reference': len(missing_in_ref),
        'Gene_Recovery_Rate': f"{(len(shared_genes)/len(ref_genes))*100:.1f}%" if len(ref_genes) > 0 else "0%",
        'Interpretation': interpretation,
        'Status': 'Success'
    }
    
    return results

def find_tool_gff_file(tool_prokka_dir, sample, plasmid_name, tool_name):
    all_dirs = os.listdir(tool_prokka_dir)

    search_patterns = []
    
    if tool_name == 'mobsuite':
        search_patterns = [
            f"{sample}_{plasmid_name}",  
            plasmid_name, 
            f"{sample}_plasmid_{plasmid_name}", 
        ]
        if plasmid_name.startswith('plasmid_'):
            clean_name = plasmid_name.replace('plasmid_', '')
            search_patterns.extend([
                f"{sample}_{clean_name}",
                clean_name,
            ])
    
    elif tool_name == 'gplas':
        if 'gplas_plasmid_bin_' in plasmid_name:
            bin_number = plasmid_name.replace('gplas_plasmid_bin_', '')
            search_patterns = [
                f"{sample}_bin_{bin_number}",  
                f"{sample}_{plasmid_name}",  
                plasmid_name,  
            ]
        else:
            search_patterns = [
                f"{sample}_{plasmid_name}",
                plasmid_name,
            ]
    
    else:  
        search_patterns = [
            f"{sample}_{plasmid_name}",
            plasmid_name,
            f"{sample}_plasmid_{plasmid_name}",
        ]
    
    for pattern in search_patterns:
        potential_gff_path = os.path.join(tool_prokka_dir, pattern, f"{pattern}.gff")
        if os.path.exists(potential_gff_path):
            logging.info(f"Found GFF file: {potential_gff_path}")
            return potential_gff_path

    logging.info(f"No exact match found for {plasmid_name}, trying partial matches...")
    for dir_name in all_dirs:
        if plasmid_name in dir_name:
            potential_gff_path = os.path.join(tool_prokka_dir, dir_name, f"{dir_name}.gff")
            if os.path.exists(potential_gff_path):
                logging.info(f"Found partial match GFF file: {potential_gff_path}")
                return potential_gff_path

    for dir_name in all_dirs:
        if sample in dir_name and any(part in dir_name for part in plasmid_name.split('_')):
            potential_gff_path = os.path.join(tool_prokka_dir, dir_name, f"{dir_name}.gff")
            if os.path.exists(potential_gff_path):
                logging.info(f"Found sample-based GFF file: {potential_gff_path}")
                return potential_gff_path
    
    return None

def run_synteny_only():
    """
    Run only synteny analysis on existing Prokka outputs using proper pairings
    """
    output_base_dir = "/genomics/Research/pfizer/validation/nov_val/Synteny/final_output"
    pairs_file = "/genomics/Research/pfizer/validation/nov_val/Synteny/sseqid_pairs.txt"

    log_file = os.path.join(output_base_dir, "synteny_only_log.txt")
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='w'),
            logging.StreamHandler()
        ]
    )
    
    logging.info("=== RUNNING PROPER SYNTENY ANALYSIS WITH CIRCULARITY AND BIDIRECTIONAL PAIRS ===")

    pairs_df = pd.read_csv(pairs_file, sep='\t')
    logging.info(f"Loaded {len(pairs_df)} pairs from {pairs_file}")

    prokka_dirs = {
        'ref': os.path.join(output_base_dir, "prokka_ref"),
        'plsmd': os.path.join(output_base_dir, "prokka_plsmd"),
        'mobsuite': os.path.join(output_base_dir, "prokka_mobsuite"),
        'gplas': os.path.join(output_base_dir, "prokka_gplas")
    }

    for tool, dir_path in prokka_dirs.items():
        if os.path.exists(dir_path):
            logging.info(f"✓ Found Prokka output for {tool}: {dir_path}")
            dirs = os.listdir(dir_path)[:5]
            logging.info(f"  First few directories: {dirs}")
        else:
            logging.error(f"✗ Missing Prokka output for {tool}: {dir_path}")
            return
    
    all_results = []
    
    for tool in ['plsmd', 'mobsuite', 'gplas']:
        ref_prokka_dir = prokka_dirs['ref']
        tool_prokka_dir = prokka_dirs[tool]
        output_csv = os.path.join(output_base_dir, f"PROPER_synteny_results_{tool}.csv")
        
        logging.info(f"\nRunning PROPER synteny analysis for {tool}...")
        logging.info(f"Reference Prokka dir: {ref_prokka_dir}")
        logging.info(f"Tool Prokka dir: {tool_prokka_dir}")
        
        results = run_synteny_analysis_with_pairs(ref_prokka_dir, tool_prokka_dir, tool, pairs_df, output_csv)
        all_results.extend(results)

    if all_results:
        combined_df = pd.DataFrame(all_results)
        combined_output = os.path.join(output_base_dir, "PROPER_synteny_results_combined.csv")
        combined_df.to_csv(combined_output, index=False)
        logging.info(f"\nAll PROPER synteny analyses completed! Combined results: {combined_output}")

    success_count = len([r for r in all_results if r['Status'] == 'Success'])
    total_count = len(all_results)
    logging.info(f"\n=== FINAL SUMMARY ===")
    logging.info(f"Successful synteny analyses: {success_count}/{total_count}")
    logging.info(f"Failed synteny analyses: {total_count - success_count}/{total_count}")

def run_synteny_analysis_with_pairs(ref_prokka_dir, tool_prokka_dir, tool_name, pairs_df, output_csv):
    """
    Run PROPER synteny analysis only for the specific pairs defined in the pairs file
    """
    results = []
    tool_col = tool_name.lower() 

    tool_pairs = pairs_df[pairs_df[tool_col].notna() & (pairs_df[tool_col] != '0') & (pairs_df[tool_col] != 0)]
    logging.info(f"Found {len(tool_pairs)} pairs for {tool_name}")

    available_dirs = os.listdir(tool_prokka_dir)
    logging.info(f"Available directories in {tool_prokka_dir}: {len(available_dirs)}")
    logging.info(f"First 10 directories: {available_dirs[:10]}")
    
    for idx, row in tool_pairs.iterrows():
        sample = row['Sample']
        sseqid = row['sseqid']
        plasmid_name = str(row[tool_col]) 

        ref_gff_path = os.path.join(ref_prokka_dir, sseqid, f"{sseqid}.gff")

        tool_gff_path = find_tool_gff_file(tool_prokka_dir, sample, plasmid_name, tool_name)
        
        if not ref_gff_path or not os.path.exists(ref_gff_path):
            logging.warning(f"Reference GFF not found: {ref_gff_path}")
            results.append({
                'Sample': sample,
                'Reference_Plasmid': sseqid,
                'Tool_Plasmid': plasmid_name,
                'Tool': tool_name,
                'NGOC': 0.0,
                'Status': 'Reference GFF not found'
            })
            continue
        
        if not tool_gff_path:
            logging.warning(f"Tool GFF not found for sample {sample}, plasmid {plasmid_name}")
            logging.warning(f"  Looked in: {tool_prokka_dir}")
            logging.warning(f"  Available dirs containing '{plasmid_name}': {[d for d in available_dirs if plasmid_name in d]}")
            results.append({
                'Sample': sample,
                'Reference_Plasmid': sseqid,
                'Tool_Plasmid': plasmid_name,
                'Tool': tool_name,
                'NGOC': 0.0,
                'Status': 'Tool GFF not found'
            })
            continue

        logging.info(f"Running PROPER synteny: {sample} - {sseqid} vs {plasmid_name} ({tool_name})")
        
        try:
            synteny_results = run_synteny_analysis(ref_gff_path, tool_gff_path, tool_name)

            synteny_results.update({
                'Sample': sample,
                'Reference_Plasmid': sseqid,
                'Tool_Plasmid': plasmid_name,
                'Tool': tool_name
            })
            
            results.append(synteny_results)
            logging.info(f"Synteny completed: {sseqid} vs {plasmid_name} - NGOC: {synteny_results['NGOC']}")
            
        except Exception as e:
            logging.error(f"Synteny failed for {sseqid} vs {plasmid_name}: {e}")
            import traceback
            traceback.print_exc()
            results.append({
                'Sample': sample,
                'Reference_Plasmid': sseqid,
                'Tool_Plasmid': plasmid_name,
                'Tool': tool_name,
                'NGOC': 0.0,
                'Status': f'Failed: {e}'
            })

    if results:
        df_results = pd.DataFrame(results)
        os.makedirs(os.path.dirname(output_csv), exist_ok=True)
        df_results.to_csv(output_csv, index=False)
        logging.info(f"Results saved to: {output_csv}")
        success_count = len([r for r in results if r['Status'] == 'Success'])
        logging.info(f"Success rate: {success_count}/{len(results)}")
    
    return results

def main():
    """
    Main function that can be used for individual synteny analysis
    """
    if len(sys.argv) == 4:
        ref_gff = sys.argv[1]
        query_gff = sys.argv[2]
        tool_name = sys.argv[3]
        
        if not os.path.exists(ref_gff):
            print(f"Error: Reference GFF file not found: {ref_gff}")
            sys.exit(1)
        
        if not os.path.exists(query_gff):
            print(f"Error: Query GFF file not found: {query_gff}")
            sys.exit(1)
        
        results = run_synteny_analysis(ref_gff, query_gff, tool_name)
        print(f"\nFinal NGOC Score: {results['NGOC']}")
        
    elif len(sys.argv) == 1:
        run_synteny_only()
    else:
        print("Usage:")
        print("  For individual analysis: python script.py <reference_gff> <query_gff> <tool_name>")
        print("  For batch analysis: python script.py")
        sys.exit(1)

if __name__ == "__main__":
    main()
