#!/usr/bin/env python3
"""
=============================================================================
 Full Validation Pipeline  –  single entry point
=============================================================================
 Steps (run automatically in order):
   1. blastn        – query (assembly) vs subject (complete genome), -subject
   2. new_cov       – subject coverage % and bases covered
   3. q_prec_cov    – query-side aligned bases per contig (qsseqid_aligned)
   4. new_multi_q   – sum_qlen and q_prec_cov
   5. extract       – keep unique sseqid rows, drop positional/score columns
   6. cat           – concatenate all extracted files into one final output

 ── MODE 1: Single sample ────────────────────────────────────────────────────
   python3 validation_pipeline.py \
       --query    assembly.fasta \
       --subject  complete_genome.fasta \
       --output   results/final_output.txt \
       --work_dir results/ \
       [--sample_name my_sample] \
       [--coverage_cutoff 0.0]

 ── MODE 2: Directory of samples ─────────────────────────────────────────────
   python3 validation_pipeline.py \
       --base_dir /path/to/samples/ \
       --output   results/final_output.txt \
       [--coverage_cutoff 0.0]

   Expected folder layout inside base_dir:
     sample_A/
       sample_A_plasmid.fasta    <- subject  (must end with _plasmid.fasta)
       sample_A_assembly.fasta   <- query    (any other *.fasta in folder)
     sample_B/
       ...

 Intermediate directories created inside --work_dir (mode 1) or
 --base_dir (mode 2):
   blast_results/   cov_results/   qprec_results/
   multi_results/   extract_results/
=============================================================================
"""

import os
import glob
import argparse
import subprocess
from collections import defaultdict

import pandas as pd
import numpy as np


# -----------------------------------------------------------------------------
# Shared helpers
# -----------------------------------------------------------------------------

def merge_intervals(intervals):
    sorted_iv = sorted(intervals, key=lambda x: x[0])
    merged = []
    for iv in sorted_iv:
        if not merged or merged[-1][1] < iv[0] - 1:
            merged.append(list(iv))
        else:
            merged[-1][1] = max(merged[-1][1], iv[1])
    return merged


def coverage_pct(merged, subject_length):
    if subject_length == 0:
        return 0.0
    covered = sum(end - start + 1 for start, end in merged)
    return covered / subject_length * 100


def fasta_total_bases(fasta_path):
    """Total sequence bases in a FASTA file (all contigs, headers excluded)."""
    if not os.path.isfile(fasta_path):
        print(f"  [WARNING] FASTA not found: {fasta_path}")
        return 0
    try:
        with open(fasta_path) as fh:
            return sum(len(l.strip()) for l in fh if not l.startswith('>'))
    except Exception as e:
        print(f"  [ERROR] Reading {fasta_path}: {e}")
        return 0


# -----------------------------------------------------------------------------
# Build the list of (sample_name, query_fasta, subject_fasta) jobs
# -----------------------------------------------------------------------------

BLAST_HEADER = (
    "qseqid\tsseqid\tqstart\tqend\tsstart\tsend\t"
    "evalue\tbitscore\tpident\tqcovs\tslen\tsample_name"
)


def collect_jobs_single(query, subject, sample_name):
    """Return one job from explicit --query / --subject arguments."""
    if not os.path.isfile(query):
        raise FileNotFoundError(f"Query file not found: {query}")
    if not os.path.isfile(subject):
        raise FileNotFoundError(f"Subject file not found: {subject}")
    if not sample_name:
        sample_name = os.path.splitext(os.path.basename(query))[0]
    return [(sample_name, query, subject)]


def collect_jobs_directory(base_dir):
    """
    Scan each sub-directory of base_dir for:
      - subject : *_plasmid.fasta
      - queries : every other *.fasta in the same sub-directory
    Returns list of (sample_name, query_fasta, subject_fasta).
    """
    jobs = []
    for sample_dir in sorted(glob.glob(os.path.join(base_dir, "*/"))):
        if not os.path.isdir(sample_dir):
            continue

        sample_name = os.path.basename(sample_dir.rstrip("/"))

        subjects = glob.glob(os.path.join(sample_dir, "*_plasmid.fasta"))
        if not subjects:
            print(f"  [WARNING] No *_plasmid.fasta in {sample_dir} – skipping.")
            continue
        subject_file = subjects[0]

        queries = [
            f for f in glob.glob(os.path.join(sample_dir, "*.fasta"))
            if not f.endswith("_plasmid.fasta")
        ]
        if not queries:
            print(f"  [WARNING] No query FASTA files in {sample_dir} – skipping.")
            continue

        for q in sorted(queries):
            jobs.append((sample_name, q, subject_file))

    return jobs


# -----------------------------------------------------------------------------
# STEP 1 – blastn
# -----------------------------------------------------------------------------

def run_blastn(jobs, blast_dir):
    print("\n" + "="*60)
    print(" STEP 1 – blastn")
    print("="*60)

    for sample_name, query_file, subject_file in jobs:
        plasmid_name = os.path.splitext(os.path.basename(query_file))[0]
        output_file  = os.path.join(blast_dir, f"{plasmid_name}_blast.txt")

        print(f"\n  Sample  : {sample_name}")
        print(f"  Query   : {query_file}")
        print(f"  Subject : {subject_file}")
        print(f"  Output  : {output_file}")
        print(f"  Running blastn ...")

        cmd = [
            "blastn",
            "-query",         query_file,
            "-subject",       subject_file,
            "-perc_identity", "90",
            "-num_threads",   "70",
            "-outfmt",
            "6 qseqid sseqid qstart qend sstart send evalue bitscore pident qcovs slen",
            "-out", output_file,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  [ERROR] blastn failed:\n{result.stderr}")
            continue

        if not os.path.isfile(output_file) or os.path.getsize(output_file) == 0:
            print(f"  [WARNING] No BLAST hits for {plasmid_name} – skipping.")
            if os.path.isfile(output_file):
                os.remove(output_file)
            continue

        hit_count = sum(1 for _ in open(output_file))
        print(f"  Hits    : {hit_count}")

        # Prepend header and append sample_name column
        tmp = output_file + ".tmp"
        with open(tmp, 'w') as out_fh:
            out_fh.write(BLAST_HEADER + "\n")
            with open(output_file) as in_fh:
                for line in in_fh:
                    out_fh.write(line.rstrip("\n") + f"\t{sample_name}\n")
        os.replace(tmp, output_file)
        print(f"  Done    : {output_file}")


# -----------------------------------------------------------------------------
# STEP 2 – new_cov  (subject coverage % and bases covered)
# -----------------------------------------------------------------------------

def step_new_cov(input_dir, output_dir, coverage_cutoff=0.0):
    print("\n" + "="*60)
    print(" STEP 2 – Subject coverage (new_cov)")
    print("="*60)

    processed = 0
    for filename in sorted(os.listdir(input_dir)):
        if not filename.endswith("_blast.txt"):
            continue

        in_path  = os.path.join(input_dir,  filename)
        out_path = os.path.join(output_dir, filename)
        print(f"\n  Processing : {in_path}")

        with open(in_path) as fh:
            lines = fh.readlines()

        if not lines:
            print("  [WARNING] File is empty – skipping.")
            continue

        header = lines[0].strip().split('\t')
        if len(header) < 11:
            print("  [WARNING] Header too short – skipping.")
            continue

        plasmid_alignments = defaultdict(list)
        plasmid_lengths    = {}

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

        plasmid_cov = {
            sid: coverage_pct(merge_intervals(ivs), plasmid_lengths[sid])
            for sid, ivs in plasmid_alignments.items()
        }

        out_rows = []
        for line in lines[1:]:
            cols = line.strip().split('\t')
            if len(cols) < 11:
                continue
            sseqid = cols[1]
            cov    = plasmid_cov.get(sseqid, 0.0)
            try:
                slen = int(cols[10])
            except ValueError:
                slen = 0
            bases = slen * cov / 100
            if cov >= coverage_cutoff:
                out_rows.append((sseqid, line.strip(), cov, bases))

        out_rows.sort(key=lambda x: x[0])

        with open(out_path, 'w') as fh:
            fh.write(lines[0].strip() + "\tcoverage_percentage\tbases_covered\n")
            for _, row, cov, bases in out_rows:
                fh.write(f"{row}\t{cov:.2f}%\t{bases:.2f}\n")

        print(f"  Written    : {out_path}  ({len(out_rows)} rows)")
        processed += 1

    if processed == 0:
        print("  [WARNING] No files processed in step 2.")


# -----------------------------------------------------------------------------
# STEP 3 – q_prec_cov_alt  (query-side aligned bases per contig)
# -----------------------------------------------------------------------------

def step_q_prec_cov(input_dir, output_dir):
    print("\n" + "="*60)
    print(" STEP 3 – Query aligned bases (q_prec_cov_alt)")
    print("="*60)

    for filename in sorted(os.listdir(input_dir)):
        if not filename.endswith("_blast.txt"):
            continue

        in_path  = os.path.join(input_dir,  filename)
        out_path = os.path.join(output_dir, filename)
        print(f"\n  Processing : {in_path}")

        df = pd.read_csv(in_path, sep='\t')
        df["qsseqid_aligned"] = 0

        sum_list = []

        for sseqid, group in df.groupby("sseqid"):
            qseqid_aligned = {}
            for qseqid in group["qseqid"].unique():
                rows    = group[group["qseqid"] == qseqid]
                max_end = int(rows["qend"].max())
                cov_arr = np.zeros(max_end + 1, dtype=bool)
                for _, row in rows.iterrows():
                    cov_arr[int(row["qstart"]):int(row["qend"]) + 1] = True
                qseqid_aligned[qseqid] = int(cov_arr.sum())

            for qseqid, aligned in qseqid_aligned.items():
                df.loc[
                    (df["sseqid"] == sseqid) & (df["qseqid"] == qseqid),
                    "qsseqid_aligned"
                ] = aligned

            sum_list.append((sseqid, sum(qseqid_aligned.values())))

        sum_df = pd.DataFrame(sum_list, columns=["sseqid", "sum_qsseqid_aligned"])
        df = df.merge(sum_df, on="sseqid")
        df.to_csv(out_path, sep='\t', index=False)
        print(f"  Written    : {out_path}  ({len(df)} rows)")


# -----------------------------------------------------------------------------
# STEP 4 – new_multi_q  (sum_qlen and q_prec_cov)
# -----------------------------------------------------------------------------

def step_new_multi_q(input_dir, output_dir, query_fasta_map):
    """
    query_fasta_map : dict  plasmid_stem -> query_fasta_path
                      used to locate the FASTA for sum_qlen calculation.
    """
    print("\n" + "="*60)
    print(" STEP 4 – qlen metrics (new_multi_q)")
    print("="*60)

    processed = 0
    for filename in sorted(os.listdir(input_dir)):
        if not filename.endswith("_blast.txt"):
            continue

        in_path  = os.path.join(input_dir,  filename)
        out_path = os.path.join(output_dir, filename)
        print(f"\n  Processing : {in_path}")

        with open(in_path) as fh:
            lines = fh.readlines()

        header = lines[0].strip().split('\t')
        required = {'sample_name', 'bases_covered', 'sum_qsseqid_aligned'}
        if not required.issubset(set(header)):
            print("  [WARNING] Missing required columns – skipping.")
            continue

        updated            = [lines[0].strip() + "\tsum_qlen\tq_prec_cov\n"]
        sample_to_sum_qlen = {}

        # Derive the query FASTA path from the filename stem
        stem       = filename.replace('_blast.txt', '')
        fasta_path = query_fasta_map.get(stem)

        # Pass 1 – resolve sum_qlen
        for line in lines[1:]:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.strip().split('\t')
            if len(cols) < len(header):
                continue
            sample_name = cols[header.index('sample_name')]
            if sample_name not in sample_to_sum_qlen:
                # Use pre-built map first, fall back to same directory
                resolved = fasta_path or os.path.join(input_dir, f"{stem}.fasta")
                print(f"    Looking for FASTA: {resolved}")
                total = fasta_total_bases(resolved)
                if total == 0:
                    print(f"    [WARNING] Empty or missing FASTA for {sample_name}")
                sample_to_sum_qlen[sample_name] = total

        # Pass 2 – write rows
        for line in lines[1:]:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.strip().split('\t')
            if len(cols) < len(header):
                continue
            sample_name         = cols[header.index('sample_name')]
            sum_qsseqid_aligned = float(cols[header.index('sum_qsseqid_aligned')])
            sum_qlen            = sample_to_sum_qlen.get(sample_name, 0)
            q_prec_cov          = (sum_qsseqid_aligned / sum_qlen * 100) if sum_qlen > 0 else 0.0
            updated.append(f"{line.strip()}\t{sum_qlen}\t{q_prec_cov:.4f}\n")

        with open(out_path, 'w') as fh:
            fh.writelines(updated)

        print(f"  Written    : {out_path}  ({len(updated)-1} rows)")
        processed += 1

    if processed == 0:
        print("  [WARNING] No files processed in step 4.")


# -----------------------------------------------------------------------------
# STEP 5 – extract  (unique sseqid, drop positional/score columns)
# -----------------------------------------------------------------------------

DROP_COLS = {"qstart", "qend", "sstart", "send", "evalue", "bitscore", "pident", "qcovs", "slen"}


def step_extract(input_dir, output_dir):
    print("\n" + "="*60)
    print(" STEP 5 – Extract unique sseqid rows (extract)")
    print("="*60)

    for filename in sorted(os.listdir(input_dir)):
        if not filename.endswith("_blast.txt"):
            continue

        in_path  = os.path.join(input_dir,  filename)
        out_path = os.path.join(output_dir, filename)
        print(f"\n  Processing : {in_path}")

        with open(in_path) as fh:
            lines = fh.readlines()

        if not lines:
            print("  [WARNING] Empty file – skipping.")
            continue

        header     = lines[0].strip().split('\t')
        keep       = [i for i, col in enumerate(header) if col not in DROP_COLS]
        sseqid_col = header.index("sseqid") if "sseqid" in header else 1

        seen_sseqid = set()
        out_lines   = ["\t".join(header[i] for i in keep) + "\n"]

        for line in lines[1:]:
            if line.startswith('#') or not line.strip():
                continue
            cols   = line.strip().split('\t')
            sseqid = cols[sseqid_col] if len(cols) > sseqid_col else ""
            if sseqid not in seen_sseqid:
                seen_sseqid.add(sseqid)
                out_lines.append("\t".join(cols[i] for i in keep) + "\n")

        with open(out_path, 'w') as fh:
            fh.writelines(out_lines)

        print(f"  Written    : {out_path}  ({len(out_lines)-1} unique sseqid rows)")


# -----------------------------------------------------------------------------
# STEP 6 – cat  (concatenate all extracted files into one final output)
# -----------------------------------------------------------------------------

def step_cat(input_dir, output_file):
    print("\n" + "="*60)
    print(" STEP 6 – Concatenate results (cat)")
    print("="*60)

    files = sorted(glob.glob(os.path.join(input_dir, "*.txt")))
    if not files:
        print(f"  [WARNING] No .txt files found in {input_dir}")
        return

    out_dir = os.path.dirname(os.path.abspath(output_file))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(output_file, 'w') as out_fh:
        for idx, filepath in enumerate(files):
            print(f"  Appending : {filepath}")
            with open(filepath) as in_fh:
                for lineno, line in enumerate(in_fh):
                    if idx > 0 and lineno == 0:
                        continue   # skip duplicate headers
                    out_fh.write(line)

    total_rows = sum(1 for _ in open(output_file)) - 1
    print(f"\n  Final output : {output_file}")
    print(f"  Total rows   : {total_rows}  (excluding header)")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Full validation pipeline:\n"
            "  blastn -> new_cov -> q_prec_cov -> new_multi_q -> extract -> cat\n\n"
            "MODE 1 – Single sample:\n"
            "  python3 validation_pipeline.py \\\n"
            "      --query    assembly.fasta \\\n"
            "      --subject  complete_genome.fasta \\\n"
            "      --output   results/final_output.txt \\\n"
            "      --work_dir results/\n\n"
            "MODE 2 – Directory of samples:\n"
            "  python3 validation_pipeline.py \\\n"
            "      --base_dir /data/samples/ \\\n"
            "      --output   results/final_output.txt\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Single-sample mode
    single = parser.add_argument_group("Single-sample mode  (--query + --subject)")
    single.add_argument("--query",
                        help="Query FASTA file (assembly / metagenome)")
    single.add_argument("--subject",
                        help="Subject FASTA file (complete genome / plasmid reference)")
    single.add_argument("--sample_name",
                        help="Label written into every output row (default: query filename stem)")
    single.add_argument("--work_dir",
                        help="Directory for intermediate files (default: parent dir of --output)")

    # Directory mode
    dirmode = parser.add_argument_group("Directory mode  (--base_dir)")
    dirmode.add_argument("--base_dir",
                         help="Root directory containing per-sample sub-directories.\n"
                              "Each sub-dir must contain *_plasmid.fasta (subject) and "
                              "at least one other *.fasta (query). "
                              "The sub-directory name is used as the sample name.")

    # Shared
    parser.add_argument("--output",          required=True,
                        help="Full path for the final concatenated output file.")
    parser.add_argument("--coverage_cutoff", type=float, default=0.0,
                        help="Minimum subject coverage %% to keep a hit (Step 2). Default: 0.0")

    args = parser.parse_args()

    # Detect mode
    single_mode    = bool(args.query or args.subject)
    directory_mode = bool(args.base_dir)

    if single_mode and directory_mode:
        parser.error(
            "Conflicting arguments: use either single-sample mode "
            "(--query / --subject) OR directory mode (--base_dir), not both."
        )
    if not single_mode and not directory_mode:
        parser.error(
            "You must provide either:\n"
            "  Single-sample mode : --query <file> --subject <file>\n"
            "  Directory mode     : --base_dir <dir>"
        )
    if single_mode and not (args.query and args.subject):
        parser.error("Single-sample mode requires both --query and --subject.")

    # Resolve working directory and jobs
    if single_mode:
        work_dir = os.path.realpath(
            args.work_dir or os.path.dirname(os.path.abspath(args.output)) or "."
        )
        jobs = collect_jobs_single(args.query, args.subject, args.sample_name)
    else:
        work_dir = os.path.realpath(args.base_dir)
        jobs     = collect_jobs_directory(work_dir)

    if not jobs:
        print("[ERROR] No valid (query, subject) pairs found. Exiting.")
        return

    print(f"\n  Mode     : {'Single sample' if single_mode else 'Directory'}")
    print(f"  Jobs     : {len(jobs)}")
    print(f"  Work dir : {work_dir}")

    # Intermediate directories
    blast_dir   = os.path.join(work_dir, "blast_results")
    cov_dir     = os.path.join(work_dir, "cov_results")
    qprec_dir   = os.path.join(work_dir, "qprec_results")
    multi_dir   = os.path.join(work_dir, "multi_results")
    extract_dir = os.path.join(work_dir, "extract_results")

    for d in (blast_dir, cov_dir, qprec_dir, multi_dir, extract_dir):
        os.makedirs(d, exist_ok=True)

    # stem -> query_fasta_path  (so step 4 can find the right FASTA)
    query_fasta_map = {
        os.path.splitext(os.path.basename(query))[0]: query
        for _, query, _ in jobs
    }

    # Run pipeline
    run_blastn      (jobs,        blast_dir)
    step_new_cov    (blast_dir,   cov_dir,     args.coverage_cutoff)
    step_q_prec_cov (cov_dir,     qprec_dir)
    step_new_multi_q(qprec_dir,   multi_dir,   query_fasta_map)
    step_extract    (multi_dir,   extract_dir)
    step_cat        (extract_dir, args.output)

    print("\n" + "="*60)
    print(" Pipeline complete!")
    print(f" Final file : {args.output}")
    print("="*60 + "\n")


if __name__ == "__main__":
    main()
