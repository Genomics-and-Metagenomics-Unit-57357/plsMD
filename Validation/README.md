# plsMD Validation

This folder contains scripts and example results to validate plsMD reconstructions against complete reference genomes.

---

## Table of Contents

- [Input Structure](#input-structure)
- [1. Dataset Validation](#1-dataset-validation)
- [2. Recall and Precision](#2-recall-and-precision)
- [3. Sensitivity and Specificity](#3-sensitivity-and-specificity)
- [4. Synteny Conservation](#4-synteny-conservation)

---

## Input Structure
```
Validation/
├── Dataset-Validation/
│   ├── Scripts/
│   │   └── Dataset_Validation.py
│   ├── failed_sample/
│   │   ├── Input/
│   │   │   ├── SRR1955549_assembly.fasta
│   │   │   └── SRR1955549_complete_genome.fasta
│   │   └── SRR1955549_validation.txt
│   └── passed_sample/
│       ├── Input/
│       │   ├── SRR18543877_assembly.fasta
│       │   └── SRR18543877_complete_genome.fasta
│       └── SRR18543877_validation.txt
├── Recall-Precision/
│   ├── Input/
│   │   └── SRR18543877_recall_precision/
│   └── SRR18543877_recall_precision.txt
├── Sensitivity-Specificity/
│   ├── Input/
│   │   ├── SRR18543877_FN/
│   │   ├── SRR18543877_FP/
│   │   ├── SRR18543877_TN/
│   │   └── SRR18543877_TP/
│   ├── SRR18543877_FN.txt
│   ├── SRR18543877_FP.txt
│   ├── SRR18543877_TN.txt
│   └── SRR18543877_TP.txt
├── Synteny/
│   ├── Reference_prokka/
│   │   └── CP144990.1.fasta/
│   │       └── CP144990.1.fasta.gff
│   └── SRR18543877_prokka/
│       └── SRR18543877_IncFII_1_pKP91/
│           └── SRR18543877_IncFII_1_pKP91.gff
└── Scripts/
    ├── Recall_precision_validation.py
    ├── Sensitivity_specificity_validation.py
    ├── best_matching_plasmids.py
    └── synteny.py
```

---

## 1. Dataset Validation

- **Goal:** Identify assemblies suitable for validation
- **Method:** BLAST original assemblies against complete genomes. Samples with low alignment coverage are excluded
- **Script:** `Dataset-Validation/Scripts/Dataset_Validation.py`
- **Output:** `Dataset-Validation/SRR1955549/SRR1955549_validation.txt`
```bash
python Dataset-Validation/Scripts/Dataset_Validation.py \
  --query Dataset-Validation/failed_sample/Input/SRR1955549_assembly.fasta \
  --subject Dataset-Validation/failed_sample/Input/SRR1955549_complete_genome.fasta \
  --output Dataset-Validation/failed_sample/SRR1955549_validation.txt
```

> The output TSV indicates whether a sample passed or failed. Failed samples are excluded from all downstream steps.

---

## 2. Recall and Precision

- **Goal:** Evaluate how accurately reconstructed plasmids match reference plasmids
- **Method:** BLAST reconstructed plasmids against complete plasmid genomes. The best-matching plasmid per reference is selected
- **Script:** `Scripts/Recall_precision_validation.py`
- **Input:** Passing samples only — per-plasmid FASTA files in `Recall-Precision/Input/SRR18543877_recall_precision/`
- **Output:** `Recall-Precision/SRR18543877_recall_precision.txt`
```bash
python Scripts/Recall_precision_validation.py \
  --query Recall-Precision/Input/SRR18543877_recall_precision/SRR18543877_IncFII_1_pKP91.fasta \
  --subject Recall-Precision/Input/SRR18543877_recall_precision/SRR18543877_plasmid.fasta \
  --output Recall-Precision/SRR18543877_recall_precision.txt
```

---

## 3. Sensitivity and Specificity

- **Goal:** Assess overall plasmid vs. non-plasmid classification accuracy
- **Method:** BLAST reconstructed plasmids and non-plasmid contigs against complete genomes. TP, FP, TN, FN are calculated to derive sensitivity and specificity
- **Script:** `Scripts/Sensitivity_specificity_validation.py`
- **Input:** Passing samples only — classified FASTA files in `Sensitivity-Specificity/Input/`
- **Output:** Per-classification result files in `Sensitivity-Specificity/`

| Input folder | Output file |
|---|---|
| `SRR18543877_TP/` | `SRR18543877_TP.txt` |
| `SRR18543877_FP/` | `SRR18543877_FP.txt` |
| `SRR18543877_TN/` | `SRR18543877_TN.txt` |
| `SRR18543877_FN/` | `SRR18543877_FN.txt` |

```bash
python Scripts/Sensitivity_specificity_validation.py \
  --query Sensitivity-Specificity/Input/SRR18543877_FN/SRR18543877_nonplasmid_contigs.fasta \
  --subject Sensitivity-Specificity/Input/SRR18543877_FN/SRR18543877_plasmid.fasta \
  --output Sensitivity-Specificity/SRR18543877_TN.tx
```

---

## 4. Synteny Conservation

- **Goal:** Compare gene order in reconstructed plasmids vs. complete plasmid genomes
- **Method:** Prokka annotations are generated for both reconstructed and reference plasmids. Gene order is compared for synteny conservation
- **Script:** `Scripts/synteny.py`
- **Input:** Passing samples only — Prokka `.gff` files from `Synteny/SRR18543877_prokka/` and `Synteny/Reference_prokka/`
- **Output:** Per-plasmid synteny comparison
```bash
python Scripts/synteny.py \
  --gff Synteny/SRR18543877_prokka/SRR18543877_IncFII_1_pKP91/SRR18543877_IncFII_1_pKP91.gff \
  --reference Synteny/Reference_prokka/CP144990.1.fasta/CP144990.1.fasta.gff \
  --pairs_file Synteny/sseqid_pairs.txt
  --output Synteny/SRR18543877_synteny.txt
```
