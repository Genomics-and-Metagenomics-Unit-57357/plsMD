plsMD: A plasmid reconstruction tool for short-read assemblies
-----
Introduction
------------
Whole genome sequencing (WGS) has become a cornerstone of antimicrobial resistance (AMR) surveillance. However, reconstructing plasmid sequences from short-read WGS data is challenging due to repetitive sequences and assembly fragmentation. 
plsMD is a tool designed for comprehensive plasmid reconstruction from short-read assemblies, going beyond simple contig binning. It utilizes Unicycler assemblies with established plasmid databases to reconstruct full plasmid sequences, plsMD offers both single and multi-sample analysis, enabling researchers to perform detailed plasmid characterization and phylogenetic investigations.
## Installation

### Prerequisites

No manual dependency installation is required. All dependencies are pre-installed
inside a dedicated `plsMD` conda environment within the Docker image.
The only requirement on your system is **Docker**.


| Tool | Version | Usage |
|------|---------|-------|
| [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) | 2.13.0 | Sequence alignment against PLSDB plasmid database |
| [Abricate](https://github.com/tseemann/abricate) | latest | Replicon typing and AMR gene detection |
| [AMRFinderPlus](https://github.com/ncbi/amr) | latest | AMR gene annotation |
| [MAFFT](https://mafft.cbrc.jp/alignment/software/) | latest | Multiple sequence alignment for phylogenetics |
| [IQ-TREE](http://www.iqtree.org/) | latest | Phylogenetic tree construction |
| [seqtk](https://github.com/lh3/seqtk) | latest | FASTA/FASTQ sequence processing |
| **rep.mob.typer sequences** | — | Replicon typing sequences, pre-configured in image |
| **IS sequences** | — | Insertion sequence database, pre-configured in image |

#### Python Libraries
- [pandas](https://pandas.pydata.org/)
- [numpy](https://numpy.org/)
- [Biopython](https://biopython.org/)
---

### 1. Clone the Repository
```bash
git clone https://github.com/Genomics-and-Metagenomics-Unit-57357/plsMD.git
cd plsMD
```

---

### 2. Build the Docker Image

Choose **one** of the following options depending on whether you already have the
PLSDB database:

**Option A — You already have the PLSDB database** *(recommended if you have it locally)*
```bash
docker build -t plsmd .
```

> The database will **not** be downloaded during the build. You will reference your
> existing PLSDB database manually when running the pipeline using the `--db` flag:
>
> ```bash
> plsMD --preprocessing --db /path/to/your/plsdb ...
> ```

**Option B — You do not have the PLSDB database** *(downloads automatically during build)*
```bash
docker build --build-arg DOWNLOAD_DB=true -t plsmd .
```

> Requires an active internet connection during the build and may take longer
> depending on your download speed.

---

### 3. Set Up Command Line Access

Make the wrapper script executable and create a system-wide symlink so you can
run `plsMD` from anywhere on your system:
```bash
# Make the wrapper executable
chmod +x plsMD_wrapper.sh

# Create a symlink in your PATH
sudo ln -sf $(pwd)/plsMD_wrapper.sh /usr/local/bin/plsMD
```

> **Important:** Run these commands from inside the `plsMD/` directory,
> otherwise `$(pwd)` will point to the wrong location.

---

### 4. Verify Installation

Confirm that plsMD is correctly installed and accessible:
```bash
plsMD --version
plsMD --help
```

You should see the version number and a list of available options.
If you get a `command not found` error, check that `/usr/local/bin` is in your `$PATH`.

------------
plsMD Pipeline
-----------
--------
Preprocessing 
------
```bash
plsMD --preprocessing --dir <input_directory> --output <preprocessing_output> --threads <num_threads> --db <plsdb_path>
```
| Option | Description | 
|---|---|
| `--dir` | Input directory containing unicylcer assembly FASTAs | 
| `--output` | Output directory containing preprocessing results | 
| `--threads` | Number of threads | 
| `--db` | PLSDB database path | 
------
Processing and reconstruction
----
```bash
plsMD --processing --dir <input_directory> --output <processing_output>
```
| Option | Description | 
|---|---|
| `--dir` | Input directory containing the preprocessed data (output from plsMD --preprocessing) | 
| `--output` | Output directory containing the processsing (reconstructed plasmids) results | 
------
Annotation
----
```bash
plsMD --annotation --dir <input_directory> --output <annotation_output> --threads <num_threads>
```
| Option | Description | 
|---|---|
| `--dir` | Input directory containing the processed data(plasmid_files/nonplasmid_files) (output from plsMD --processing) | 
| `--dir` | Input directory containing the preprocessed data (output from plsMD --preprocessing) | 
| `--threads` | Number of threads | 
-------
plsMD batch modality
------
In addition to the preprocessing and processing steps, run --phylogenetics for batch samples modality.

------
Phylogenetic Analysis
------
```bash
plsMD --phylogenetics --genes_dir <genes_directory> --output <phylogenetics_output> --min_length <min_length> --max_length <max_length> --threads <threads>
```
| Option | Description | 
|---|---|
| `--genes_dir` |Root directory containing extracted fasta folders. Each subdirectory should represent a plasmid replicon | 
| `--output` |Output directory containing phylogenetic results | 
| `--threads` | Number of threads | 
| `--min_length` | Minimum sequence length for common sequence search. |
| `--max_length` | Maximum sequence length for common sequence search. |
---------
Output directories/files
| Directory/File | Description | Generated by |
|---|---|---|
| `<preprocessing_input>/` | **Input directory** — should contain only the Unicycler FASTA assemblies, one per sample | User provided |
| `<preprocessing_output>/<fasta_file>_plasmid.txt` | Copy of merged results for downstream processing | `plsMD --preprocessing` |
| `<preprocessing_output>/<fasta_file>_PLSDB.txt` | BLASTn alignments against the PLSDB database | `plsMD --preprocessing` |
| `<preprocessing_output>/<fasta_file>_merged.fasta` | Original and reversed contigs merged into a single FASTA | `plsMD --preprocessing` |
| `<preprocessing_output>/plasmid_list.txt` | List of identified plasmid replicons | `plsMD --preprocessing` |
| `<preprocessing_output>/unicycler_fasta/` | Copy of input Unicycler FASTA assemblies | `plsMD --preprocessing` |
| `<preprocessing_output>/plasmidfinder_results/<fasta_file>_plasmid.txt` | PlasmidFinder replicon typing results via Abricate | `plsMD --preprocessing` |
| `<preprocessing_output>/rep_typer_results/<fasta_file>_plasmid.txt` | rep.mob.typer replicon typing results via Abricate | `plsMD --preprocessing` |
| `<preprocessing_output>/merged_plasmid_results/<fasta_file>_plasmid.txt` | Merged PlasmidFinder and rep.mob.typer results with redundant entries removed | `plsMD --preprocessing` |
| `<processing_output>/Processed_ColContigs.txt` | Processed Col and Inc plasmid contig assignments | `plsMD --processing` |
| `<processing_output>/gene_directories/` | Per-replicon directories containing alignment and filtering intermediates | `plsMD --processing` |
| `<processing_output>/extracted_fasta/` | Reconstructed FASTA sequences per plasmid replicon | `plsMD --processing` |
| `<processing_output>/plasmid_files/` | Final reconstructed plasmid sequences and reports per sample | `plsMD --processing` |
| `<processing_output>/nonplasmid_files/` | Non-plasmid contigs per sample | `plsMD --processing` |
| `<annotation_output>/` | AMRFinderPlus, Abricate (VFDB, PlasmidFinder) and BLASTn (IS) results for plasmid and non-plasmid sequences | `plsMD --annotation` |
| `<phylogenetics_output>/phylogenetic_tree/` | MSA and phylogenetic analysis output | `plsMD --phylogenetics` |
| `<phylogenetics_output>/rotated_sequences/` | Concatenated and rotated sequences | `plsMD --phylogenetics` |
----------
## Citations

If you use plsMD in your research, please cite the following tools:
- **plsMD** — Lotfi M, Jalal D, et al. bioRxiv 2025.03.17.643493; doi: https://doi.org/10.1101/2025.03.17.643493 **GitHub**. https://github.com/Genomics-and-Metagenomics-Unit-57357/plsMD, 
- **BLAST+** — Camacho et al. (2009) *BMC Bioinformatics* 10:421. https://doi.org/10.1186/1471-2105-10-421
- **Abricate** — Seemann T. GitHub. https://github.com/tseemann/abricate
- **AMRFinderPlus** — Feldgarden et al. (2021) *Sci Rep* 11:12728. https://doi.org/10.1038/s41598-021-91456-0
- **MAFFT** — Katoh & Standley (2013) *Mol Biol Evol* 30:772–780. https://doi.org/10.1093/molbev/mst010
- **IQ-TREE** — Nguyen et al. (2015) *Mol Biol Evol* 32(1):268–274. https://doi.org/10.1093/molbev/msu300
- **ISFinder** — Siguier et al. (2006) *Nucleic Acids Res* 34:D32–D36. https://www.ncbi.nlm.nih.gov/pubmed/16381877 — Database: http://www-is.biotoul.fr
- **mob-suite** — Robertson & Nash (2018) *Microbial Genomics* 4(8):e000206. https://doi.org/10.1099/mgen.0.000206
