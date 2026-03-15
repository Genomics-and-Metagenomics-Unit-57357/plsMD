# plsMD Validation

This folder contains scripts and example results to validate plsMD reconstructions against complete reference genomes.

---

## Input Structure

Each sample should have **two FASTA files**:

1. `assembly.fasta` – The sample’s assembled genome (from Unicycler or similar).  
2. `genome.fasta` – The corresponding complete reference genome.  

For this tutorial, we include one **failed** and one **passed** sample:
```
Validation/
└── input/
├── failed_sample/
│ ├── assembly.fasta
│ └── genome.fasta
└── passed_sample/
├── assembly.fasta
└── genome.fasta

```

---

## 1. Dataset Validation

- **Goal:** Identify assemblies suitable for validation.  
- **Method:** BLAST original assemblies against complete genomes. Samples with low alignment coverage are excluded.  
- **Script:** `scripts/01_dataset_validation.sh`  
- **Output:** `results/assembly_vs_genomes.tsv`  

```bash
bash scripts/01_dataset_validation.sh --input input/ --output results/assembly_vs_genomes.tsv
```

The TSV file will indicate which sample passed or failed. Failed samples are excluded from downstream steps.


---

## 2. Recall and Precision

- **Goal:** Evaluate how accurately reconstructed plasmids match reference plasmids.  
- **Method:** BLAST reconstructed plasmids against complete plasmid genomes. The best-matching plasmid per reference is selected.  
- **Script:** `scripts/02_recall_precision.sh`
- **Input** Only passing samples.   
- **Output:** `results/recall_precision.tsv`  

```bash
bash scripts/02_recall_precision.sh --input input/passed_sample/ --output results/recall_precision.tsv
```
---

## 3. Sensitivity and Specificity

- **Goal:** Assess overall plasmid vs non-plasmid classification.  
- **Method:** BLAST reconstructed plasmids and non-plasmid contigs against complete genomes. Calculate TP, FP, TN, FN → derive sensitivity & specificity.  
- **Script:** `scripts/03_sensitivity_specificity.sh`
- **Input** Only passing samples.   
- **Output:** `results/sensitivity_specificity.tsv`  

```bash
bash scripts/03_sensitivity_specificity.sh --input input/passed_sample/ --output results/sensitivity_specificity.tsv
```
---

## 4. Synteny Conservation

- **Goal:** Compare gene order in reconstructed plasmids vs complete plasmid genomes.  
- **Method:** Annotate complete genomes using Prokka, extract gene order in reconstructed plasmids, and compare for synteny conservation.  
- **Script:** `scripts/04_synteny_conservation.sh`
- **Input** Only passing samples.   
- **Output:** `results/synteny_comparison.tsv`  

```bash
bash scripts/04_synteny_conservation.sh --input input/passed_sample/ --output results/synteny_comparison.tsv
```
---


### Folder Structure
```
Validation/
├── README.md
├── input/
│   ├── failed_sample/
│   │   ├── assembly.fasta
│   │   └── genome.fasta
│   └── passed_sample/
│       ├── assembly.fasta
│       └── genome.fasta
├── scripts/
│   ├── 01_dataset_validation.sh
│   ├── 02_recall_precision.sh
│   ├── 03_sensitivity_specificity.sh
│   └── 04_synteny_conservation.sh
└── results/
    ├── assembly_vs_genomes.tsv
    ├── recall_precision.tsv
    ├── sensitivity_specificity.tsv
    └── synteny_comparison.tsv
```
