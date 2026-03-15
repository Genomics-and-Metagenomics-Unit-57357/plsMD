# plsMD Validation

This folder contains scripts and example results to validate plsMD reconstructions against complete reference genomes.

---

## 1. Dataset Validation
- **Goal:** Identify assemblies suitable for validation.  
- **Method:** BLAST original assemblies against complete genomes. Samples with low alignment coverage are excluded.  
- **Script:** `scripts/01_dataset_validation.sh`  
- **Outputs:** `results/assembly_vs_genomes.tsv`  

---

## 2. Recall and Precision
- **Goal:** Evaluate how accurately reconstructed plasmids match reference plasmids.  
- **Method:** BLAST reconstructed plasmids against complete plasmid genomes. Select the best-matching plasmid per reference.  
- **Script:** `scripts/02_recall_precision.sh`  
- **Outputs:** `results/recall_precision.tsv`  

---

## 3. Sensitivity and Specificity
- **Goal:** Assess overall plasmid vs non-plasmid classification.  
- **Method:** BLAST reconstructed plasmids and non-plasmid contigs against complete genomes. Calculate TP, FP, TN, FN → derive sensitivity & specificity.  
- **Script:** `scripts/03_sensitivity_specificity.sh`  
- **Outputs:** `results/sensitivity_specificity.tsv`  

---

## 4. Synteny Conservation
- **Goal:** Compare gene order in reconstructed plasmids vs complete plasmid genomes.  
- **Method:** Annotate complete genomes using Prokka, extract gene order in reconstructed plasmids, compare for synteny conservation.  
- **Script:** `scripts/04_synteny_conservation.sh`  
- **Outputs:** `results/synteny_comparison.tsv`  

---

### Folder Structure
Validation/
  ├── README.md
  ├── scripts/
  │ ├── 01_dataset_validation.sh
  │ ├── 02_recall_precision.sh
  │ ├── 03_sensitivity_specificity.sh
  │ └── 04_synteny_conservation.sh
  └── results/
  ├── assembly_vs_genomes.tsv
  ├── recall_precision.tsv
  ├── sensitivity_specificity.tsv
  └── synteny_comparison.tsv
