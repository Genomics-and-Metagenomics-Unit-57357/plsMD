# plsMD Step-by-step Tutorial

## Table of Contents

- [Objective](#objective)
- [plsMD Pipeline Overview](#plsmd-pipeline-overview)
- [Setting Up plsMD](#setting-up-plsmd)
- [Step 1: Preprocessing](#step-1-preprocessing)
- [Step 2: Processing](#step-2-processing)
- [Step 3: Annotation](#step-3-annotation)
- [Step 4: Phylogenetic Analysis](#step-4-phylogenetic-analysis)

## Objective

The objective of this tutorial is to illustrate a complete **plsMD pipeline** for short-read whole genome sequencing bacterial samples.

---

## plsMD Pipeline Overview

The plsMD pipeline is divided into the following steps:

- **Preprocessing** – Identifies plasmid replicons to be used as seeds and generates plasmid alignments.
- **Processing** – Uses the seeds and plasmid alignments to create reconstructed plasmids.
- **Annotation** – Reconstructed plasmids are annotated for AMR genes, VF genes, and ISs.
- **Batch Modality** – Reconstructed plasmids with the same replicon type are grouped together for phylogenetic analysis.

---

## Setting Up plsMD

1. **Install plsMD**
Follow the [Installation guide](../README.md#installation) in the main README.md

2. **Download tutorial example files**
```bash
git clone https://github.com/Genomics-and-Metagenomics-Unit-57357/plsMD.git
cd "plsMD/Step-By-Step Tutorial/SRR18543877"
```

---

## Step 1: Preprocessing

Run the preprocessing step to identify plasmid replicons and generate plasmid alignments.
```bash
plsMD --preprocessing \
  --dir <input_directory>/ \
  --output <output_directory> \
  [--db <plsdb_path>] \
  --threads <num_threads>
```

### Preprocessing Output Structure
```
Preprocessing_output/
├── SRR18543877_plasmid.txt
├── SRR18543877_PLSDB.txt
├── SRR18543877_merged.fasta
├── plasmid_list.txt
└── Unicycler_fasta/
    └── SRR18543877.fasta
```

### Output Description

- The FASTA assembly is run on **PlasmidFinder** (`plasmidfinder_results/`) and **rep.mob.typer** (`rep_typer_results/`).
- Both outputs are combined, with redundant genes removed, producing `SRR18543877_plasmid.txt`.
- All contigs are reversed, `_R` is appended to the contig name, and merged with the original contigs into `SRR18543877_merged.fasta`.
- Unique replicon gene names from all samples are concatenated in a file `plasmid_list.txt`

--

## Step 2: Processing

Run the processing step to generate reconstructed plasmids.
```bash
plsMD --processing \
  --dir <input_directory (preprocessing_output)/> \
  --output <output_directory> 
```
### Processing Output Structure
```
Processing_output/
├── SRR18543877_overlap.txt
├── SRR18543877_overlap_filtered.txt
├── extracted_fasta
    └── Col440I_1_pld1/SRR18543877_Col440I_1_pld1.fasta
    └── Col440I_1_pld2/SRR18543877_Col440I_1_pld2.fasta
    └── IncFIB_K__1_Kpn3_pld1/SRR18543877_IncFIB_K__1_Kpn3_pld1.fasta
    └── IncFII_1_pKP91_pld1/SRR18543877_IncFII_1_pKP91_pld1.fasta
    └── IncR_1_pld1/SRR18543877_IncR_1_pld1.fasta
├── gene_directories
    └── Col440I_1_pld1/
        └── SRR18543877_plasmid.txt 
        └── SRR18543877_overlap_filtered.txt
        └── SRR18543877_replicon_contigs.txt
        └── SRR18543877_replicon_filtered.txt
        └── SRR18543877_extracted.fasta   
        └── SRR18543877_NZ_CP143284.1.fasta 
    └── Col440I_1_pld2/
    └── IncFIB_K__1_Kpn3_pld1/
    └── IncFII_1_pKP91_pld1/
    └── IncR_1_pld1/
├── plasmid_files/
    └── SRR18543877_plasmid_contigs.fasta
    └── SRR18543877_plasmid_replicon_filtered.txt
    └── SRR18543877_report.tsv
├── nonplasmid_files/
    └── SRR18543877_nonplasmid_contigs.fasta
```
### Output Description

- The PLSDB alignments are processed as follows:

  1. **Alignments in reverse direction** relative to the reference plasmid are removed.
     Overlapped bases between overlapping alignments are calculated → `SRR18543877_overlap.txt`

  2. **Nested and spurious alignments** are filtered out.
     Aligned bases relative to the reference plasmid are calculated:
     - `bases_covered` — absolute bases covered on the reference plasmid
     - `coverage_percentage` — percentage relative to the length of the reference plasmid
     - Output: `SRR18543877_overlap_filtered.txt`

  3. Directories are created for each gene listed in `plasmid_list.txt` → `gene_directories/`

  4. Within each gene directory, reference plasmid alignments that **do not include the contig
     harboring the respective replicon** are filtered out → `SRR18543877_replicon_contigs.txt`

  5. The **single best reference plasmid** per replicon per sample is selected
     → `SRR18543877_replicon_filtered.txt`

  6. Contig sequences are extracted from `SRR18543877_merged.fasta` and **trimmed for
     overlapped bases** → `SRR18543877_extracted.fasta`

  7. Trimmed sequences are **merged together** into a full reconstructed plasmid
     → `Col440I_1_pld1/SRR18543877_NZ_CP143284.1.fasta`
     The FASTA name (NZ_CP143284.1) is the reference plasmid used as a guide.

  9. All `SRR18543877_replicon_filtered.txt` files across all gene directories concatenated for that sample `plasmid_files/SRR18543877_plasmid_replicon_filtered.txt`. In this step, plasmid names are renamed to include the replicon name instead of the name of the reference plasmid used.
 
  10. FASTA file with all reconstructed plasmid sequences for that sample `plasmid_files/SRR18543877_plasmid_contigs.fasta`.
      In this step, contigs that are tagged as circular in the original unicycler assembly are included regardless of the inclusion of a replicon.
 
  11. Reports for each sample are generated for all reconstructed plasmids and circular non-replicon plasmids. `plasmid_files/SRR18543877_report.tsv`.
 
 
  12. All contigs from the original assembly that were not identified as plasmid-associated are extracted `nonplasmid_files/SRR18543877_nonplasmid_contigs.fasta`.
 
  --

## Step 3: Annotation
Run the annotation step to align the reconstructed and nonplasmid sequences against multiple databases `(amrfinderplus-vfdb-plasmidfinder-ISFinder)`.
```bash
plsMD --annotation \
  --dir <input_directory (processing_output)> \
  --output <output_directory>
  --threads <num_threads>
```
### Annotation Output Structure
```
Annotation_output/
├── plasmid_files/
    └── AMR/SRR18543877_plasmid_contigs_AMR.txt
    └── IS/SRR18543877_plasmid_contigs_IS.txt
    └── PL/SRR18543877_plasmid_contigs_PL.txt
    └── VF/SRR18543877_plasmid_contigs_VF.txt

├── nonplasmid_files/
    └── AMR/SRR18543877_nonplasmid_contigs_AMR.txt
    └── IS/SRR18543877_nonplasmid_contigs_IS.txt
    └── PL/SRR18543877_nonplasmid_contigs_PL.txt
    └── VF/SRR18543877_nonplasmid_contigs_VF.txt
```
--

## Step 4: Phylogenetic analysis  
If multiple samples are processed together, this step allows the phylogenetic analysis of plasmids from different samples harboring the same replicon type.
The user manually needs to copy the desired extracted plasmid directories into a new directory `input_directory`.
```
├── extracted_fasta
    └── Col440I_1_pld1/   →    ├── new_directory/
                                   └── Col440I_1_pld1/
```

```bash
plsMD --phylogenetics \
  --genes_dir <input_directory>/ \
  --output <output_directory> 
  --threads <num_threads>
```
### Phylogenetics Output Structure
```
├── phylogenetics_output/
    └── rotated_sequences/
        └── Col440I_1_pld1/Col440I_1_pld1_concatenated.fasta
        └── Col440I_1_pld1/Col440I_1_pld1_commn_seq.txt
        └── Col440I_1_pld1/Col440I_1_pld1_rotated.fasta
    └── Phylogenetic_analysis/
        └── Col440I_1_pld1/Col440I_1_pld1_aligned.fasta
        └── Col440I_1_pld1/Col440I_1.treefile
 
   ```
### Output Description 

  1. Extracted plasmids from all samples in teh same replicon gene folder are concatenated into a multi-fasta file`Col440I_1_pld1_concatenated.fasta`.

  2. A common sequence is retrieved from all sequences where it exists only once per plasmid sequence to be used as a common stating point for rotation. `Col440I_1_pld1_commn_seq.txt`

  3. The common sequence is used to rotate the sequences and reverse if needed `Col440I_1_pld1_rotated.fasta`.

  4.  Multiple sequence alignment is performed on the rotated sequences `Col440I_1_pld1_aligned.fasta`.
  
  5.  Phylogenetic tree files are generated `Col440I_1.treefile`. 
  




