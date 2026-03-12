#!/usr/bin/env python3
"""
Usage:
    python blast_pipeline.py \
        --query <query.fasta> \
        --subject <subject.fasta> \
        --output_dir <output_directory> \
        [--sample_name <name>] \
        [--perc_identity 90] \
        [--num_threads 70] \
        [--coverage_cutoff 0.0]

Pipeline steps:
    1. Run blastn
    2. Calculate subject coverage percentage and bases covered  (new_cov)
    3. Calculate query-side aligned bases per contig             (q_prec_cov_alt)
    4. Add qlen, sum_qlen, q_prec_cov, contig_cov, circular     (multi_qlen_cov)
"""

import os
import sys
import argparse
import subprocess
import pandas as pd
import numpy as np
from collections import defaultdict



def merge_intervals(intervals):
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    merged = []
    for iv in sorted_intervals:
        if not merged or merged[-1][1] < iv[0] - 1:
            merged.append(list(iv))
        else:
            merged[-1][1] = max(merged[-1][1], iv[1])
    return merged


def calculate_coverage_percentage(merged_intervals, subject_length):
    if subject_length == 0:
        return 0.0
    covered = sum(end - start + 1 for start, end in merged_intervals)
    return (covered / subject_length) * 100


def calculate_fasta_sum(fasta_file):
    """Total number of sequence bases in a FASTA file (headers excluded)."""
    if not os.path.isfile(fasta_file):
        print(f"  [WARNING] FASTA file not found: {fasta_file}")
        return 0
    try:
        with open(fasta_file) as f:
            seq = ''.join(l.strip() for l in f if not l.startswith('>'))
        return len(seq)
    except Exception as e:
        print(f"  [ERROR] Reading FASTA {fasta_file}: {e}")
        return 0


def get_contig_info(fasta_file):
    """
    Parse a FASTA file and return:
        { contig_id: (length, circular_flag) }
    circular_flag is "true" when the header contains 'circular=true'.
    """
    contig_info = {}
    if not os.path.isfile(fasta_file):
        print(f"  [WARNING] FASTA file not found: {fasta_file}")
        return contig_info
    try:
        current_id = current_header = None
        current_seq = []
        with open(fasta_file) as f:
            for raw in f:
                line = raw.strip()
                if line.startswith('>'):
                    if current_id is not None:
                        seq = ''.join(current_seq)
                        circ = "true" if current_header and "circular=true" in current_header.lower() else "false"
                        contig_info[current_id] = (len(seq), circ)
                    current_header = line
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
        if current_id is not None and current_seq:
            seq = ''.join(current_seq)
            circ = "true" if current_header and "circular=true" in current_header.lower() else "false"
            contig_info[current_id] = (len(seq), circ)
    except Exception as e:
        print(f"  [ERROR] Parsing FASTA {fasta_file}: {e}")
    return contig_info


# Run blastn

def run_blastn(query, subject, blast_out, perc_identity, num_threads):
    print("\n[Step 1] Running blastn …")
    cmd = [
        "blastn",
        "-query", query,
        "-subject", subject,
        "-perc_identity", str(perc_identity),
        "-num_threads", str(num_threads),
        "-outfmt", "6 qseqid sseqid qstart qend sstart send evalue bitscore pident qcovs slen",
        "-out", blast_out,
    ]
    print("  CMD:", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  [ERROR] blastn failed:\n{result.stderr}")
        sys.exit(1)

    cols = "qseqid\tsseqid\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tpident\tqcovs\tslen"
    with open(blast_out) as fh:
        body = fh.read()
    with open(blast_out, 'w') as fh:
        fh.write(cols + "\n" + body)

    lines = body.strip().splitlines()
    print(f"  Done. {len(lines)} alignment(s) written to {blast_out}")
    return blast_out


# Step 2 – Subject coverage 

def add_subject_coverage(input_file, output_file, sample_name, coverage_cutoff=0.0):
    print("\n[Step 2] Calculating subject coverage (new_cov) …")

    with open(input_file) as fh:
        lines = fh.readlines()

    if not lines:
        print("  [ERROR] BLAST output is empty.")
        sys.exit(1)

    header = lines[0].strip().split('\t')
    if len(header) < 11:
        print(f"  [ERROR] Unexpected header: {header}")
        sys.exit(1)

    plasmid_alignments = defaultdict(list)
    plasmid_lengths = {}

    for line in lines[1:]:
        cols = line.strip().split('\t')
        if len(cols) < 11:
            continue
        try:
            sseqid = cols[1]
            sstart, send, slen = int(cols[4]), int(cols[5]), int(cols[10])
            if sstart > send:
                sstart, send = send, sstart
            plasmid_alignments[sseqid].append((sstart, send))
            plasmid_lengths[sseqid] = slen
        except ValueError:
            continue

    plasmid_coverage = {
        sid: calculate_coverage_percentage(merge_intervals(ivs), plasmid_lengths[sid])
        for sid, ivs in plasmid_alignments.items()
    }

    out_rows = []
    for line in lines[1:]:
        cols = line.strip().split('\t')
        if len(cols) < 11:
            continue
        sseqid = cols[1]
        cov_pct = plasmid_coverage.get(sseqid, 0.0)
        try:
            slen = int(cols[10])
        except ValueError:
            slen = 0
        bases_cov = slen * cov_pct / 100
        if cov_pct >= coverage_cutoff:
            out_rows.append((sseqid, line.strip(), cov_pct, bases_cov))

    out_rows.sort(key=lambda x: x[0])

    with open(output_file, 'w') as fh:
        fh.write(f"sample_name\t{lines[0].strip()}\tcoverage_percentage\tbases_covered\n")
        for _, line, cov_pct, bases_cov in out_rows:
            fh.write(f"{sample_name}\t{line}\t{cov_pct:.2f}%\t{bases_cov:.2f}\n")

    print(f"  Done. {len(out_rows)} row(s) written to {output_file}")
    return output_file



#Query-side aligned bases

def add_query_aligned(input_file, output_file):
    print("\n[Step 3] Calculating query-aligned bases (q_prec_cov_alt) …")

    df = pd.read_csv(input_file, sep='\t')
    df["qsseqid_aligned"] = 0

    qseqid_aligned_map = {}
    sum_per_sseqid = []

    for sseqid, group in df.groupby("sseqid"):
        qseqid_aligned = {}
        for qseqid in group["qseqid"].unique():
            rows = group[group["qseqid"] == qseqid]
            max_end = int(rows["qend"].max())
            coverage = np.zeros(max_end + 1, dtype=bool)
            for _, row in rows.iterrows():
                coverage[int(row["qstart"]):int(row["qend"]) + 1] = True
            qseqid_aligned[qseqid] = int(coverage.sum())

        for qseqid, aligned in qseqid_aligned.items():
            qseqid_aligned_map[(sseqid, qseqid)] = aligned

        sum_per_sseqid.append((sseqid, sum(qseqid_aligned.values())))

    df["qsseqid_aligned"] = df.apply(
        lambda r: qseqid_aligned_map.get((r["sseqid"], r["qseqid"]), 0), axis=1
    )

    sum_df = pd.DataFrame(sum_per_sseqid, columns=["sseqid", "sum_qsseqid_aligned"])
    df = df.merge(sum_df, on="sseqid")

    df.to_csv(output_file, sep='\t', index=False)
    print(f"  Done. {len(df)} row(s) written to {output_file}")
    return output_file


#qlen / contig coverage / circularity

def add_qlen_metrics(input_file, output_file, query_fasta):
    print("\n[Step 4] Calculating qlen / contig coverage / circular (multi_qlen_cov) …")

    with open(input_file) as fh:
        lines = fh.readlines()

    header = lines[0].strip().split('\t')
    required = {'sample_name', 'bases_covered', 'sum_qsseqid_aligned', 'qseqid', 'qsseqid_aligned', 'sseqid'}
    if not required.issubset(set(header)):
        missing = required - set(header)
        print(f"  [ERROR] Missing columns: {missing}")
        sys.exit(1)

    sum_qlen     = calculate_fasta_sum(query_fasta)
    contig_info  = get_contig_info(query_fasta)
    print(f"  Query genome total bases : {sum_qlen}")
    print(f"  Query contigs found      : {len(contig_info)}")

    sample_to_sum_qlen    = {}
    sample_to_contig_info = {}


    for line in lines[1:]:
        if line.startswith('#') or not line.strip():
            continue
        cols = line.strip().split('\t')
        if len(cols) < len(header):
            continue
        sname = cols[header.index('sample_name')]
        if sname not in sample_to_sum_qlen:
            sample_to_sum_qlen[sname]    = sum_qlen
            sample_to_contig_info[sname] = contig_info


    qseqid_to_max_contig_cov = defaultdict(float)
    aligned_qseqids = set()

    for line in lines[1:]:
        if line.startswith('#') or not line.strip():
            continue
        cols = line.strip().split('\t')
        if len(cols) < len(header):
            continue
        sname          = cols[header.index('sample_name')]
        qseqid         = cols[header.index('qseqid')]
        qsseqid_aligned = float(cols[header.index('qsseqid_aligned')])
        aligned_qseqids.add(qseqid)

        ci = sample_to_contig_info.get(sname, {})
        qlen, _ = ci.get(qseqid, (0, "false"))
        contig_cov = (qsseqid_aligned / qlen * 100) if qlen > 0 else 0.0
        if contig_cov > qseqid_to_max_contig_cov[qseqid]:
            qseqid_to_max_contig_cov[qseqid] = contig_cov


    updated = [lines[0].strip() + '\tqlen\tsum_qlen\tq_prec_cov\tcontig_cov\tmax_contig_cov\tcircular\n']

    for line in lines[1:]:
        if line.startswith('#') or not line.strip():
            continue
        cols = line.strip().split('\t')
        if len(cols) < len(header):
            continue

        sname              = cols[header.index('sample_name')]
        qseqid             = cols[header.index('qseqid')]
        sum_qsseqid_aligned = float(cols[header.index('sum_qsseqid_aligned')])
        qsseqid_aligned     = float(cols[header.index('qsseqid_aligned')])

        sq   = sample_to_sum_qlen.get(sname, 0)
        ci   = sample_to_contig_info.get(sname, {})
        qlen, circular_flag = ci.get(qseqid, (0, "false"))

        q_prec_cov   = (sum_qsseqid_aligned / sq  * 100) if sq   > 0 else 0.0
        contig_cov   = (qsseqid_aligned      / qlen * 100) if qlen > 0 else 0.0
        max_contig_cov = qseqid_to_max_contig_cov.get(qseqid, 0.0)

        updated.append(
            f"{line.strip()}\t{qlen}\t{sq}\t{q_prec_cov:.4f}\t"
            f"{contig_cov:.4f}\t{max_contig_cov:.4f}\t{circular_flag}\n"
        )


    for sname, ci in sample_to_contig_info.items():
        for qseqid, (qlen, circular_flag) in ci.items():
            if qseqid not in aligned_qseqids:
                empty_vals = []
                for h in header:
                    if h == 'sample_name':
                        empty_vals.append(sname)
                    elif h == 'qseqid':
                        empty_vals.append(qseqid)
                    else:
                        empty_vals.append("0")
                max_contig_cov = qseqid_to_max_contig_cov.get(qseqid, 0.0)
                updated.append(
                    f"{chr(9).join(empty_vals)}\t{qlen}\t0\t0\t0\t{max_contig_cov:.4f}\t{circular_flag}\n"
                )

    with open(output_file, 'w') as fh:
        fh.writelines(updated)

    print(f"  Done. {len(updated)-1} row(s) written to {output_file}")
    return output_file


def main():
    parser = argparse.ArgumentParser(
        description="BLAST pipeline: blastn → subject coverage → query coverage → qlen/contig metrics",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--query",         required=True,  help="Query FASTA file")
    parser.add_argument("--subject",       required=True,  help="Subject FASTA file")
    parser.add_argument("--output_dir",    required=True,  help="Directory for all intermediate and final outputs")
    parser.add_argument("--sample_name",   default=None,   help="Sample name tag added to every row (default: query filename stem)")
    parser.add_argument("--perc_identity", type=float, default=90,  help="blastn -perc_identity")
    parser.add_argument("--num_threads",   type=int,   default=70,  help="blastn -num_threads")
    parser.add_argument("--coverage_cutoff", type=float, default=0.0,
                        help="Minimum subject coverage %% to keep a hit (Step 2)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    sample_name = args.sample_name or os.path.splitext(os.path.basename(args.query))[0]
    base_name   = f"{sample_name}_blast"

    # File paths for each stage
    f_blast   = os.path.join(args.output_dir, f"{base_name}.txt")          # Step 1 raw BLAST
    f_step2   = os.path.join(args.output_dir, f"{base_name}_cov.txt")      # After subject coverage
    f_step3   = os.path.join(args.output_dir, f"{base_name}_qcov.txt")     # After query coverage
    f_final   = os.path.join(args.output_dir, f"{base_name}_final.txt")    # Final output

    run_blastn(args.query, args.subject, f_blast, args.perc_identity, args.num_threads)
    add_subject_coverage(f_blast,  f_step2, sample_name, args.coverage_cutoff)
    add_query_aligned(f_step2, f_step3)
    add_qlen_metrics(f_step3, f_final, args.query)

    print(f"\n{'='*60}")
    print(f"  Pipeline complete!")
    print(f"  Final output : {f_final}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
