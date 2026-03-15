plsMD: A plasmid reconstruction tool for short-read assemblies
-----
Introduction
------------
Whole genome sequencing (WGS) has become a cornerstone of antimicrobial resistance (AMR) surveillance. However, reconstructing plasmid sequences from short-read WGS data is challenging due to repetitive sequences and assembly fragmentation. 
plsMD is a tool designed for comprehensive plasmid reconstruction from short-read assemblies, going beyond simple contig binning. It utilizes Unicycler assemblies with established plasmid databases to reconstruct full plasmid sequences, plsMD offers both single and multi-sample analysis, enabling researchers to perform detailed plasmid characterization and phylogenetic investigations.
## Installation

### Prerequisites

No manual dependency installation is required. All dependencies are bundled and
pre-configured inside the Docker image within a dedicated `plsMD` conda environment.

#### Bioinformatics Tools

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

| Library | Usage |
|---------|-------|
| [pandas](https://pandas.pydata.org/) | Data processing and tabular file handling |
| [numpy](https://numpy.org/) | Numerical operations |
| [Biopython](https://biopython.org/) | FASTA parsing and sequence manipulation |

> All tools above are installed and made compatible within the Docker image.
> The only requirement on your system is **Docker**.

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
plsMD Single modality
-----------
--------
Preprocessing 
------
```bash
plsMD --preprocessing --dir <input_directory> --threads <num_threads> --db <plsdb_path>
```
| Option | Description | 
|---|---|
| `--dir` | Input directory containing unicylcer assembly FASTAs | 
| `--threads` | Number of threads | 
| `--db` | PLSDB database path | 
------
Processing and reconstruction
----
```bash
plsMD --processing --dir <input_directory>
```
| Option | Description | 
|---|---|
| `--dir` | Input directory containing the preprocessed data (output from plsMD --preprocessing) | 
------
Annotation
----
```bash
plsMD --annotation --dir <input_directory> --threads <num_threads> --IS_db <IS_db_path>
```
| Option | Description | 
|---|---|
| `--dir` | Input directory containing the preprocessed data (output from plsMD --preprocessing) | 
| `--threads` | Number of threads | 
| `--IS_db` | Insertion Sequences Database path | 
-------
plsMD batch modality
------
In addition to the preprocessing and processing steps, run --phylogenetics for batch samples modality.

------
Phylogenetic Analysis
------
```bash
plsMD --phylogenetics --genes_dir <genes_directory> --min_length <min_length> --max_length <max_length> --threads <threads>
```
| Option | Description | 
|---|---|
| `--genes_dir` |Root directory containing gene folders. Each subdirectory should represent a plasmid replicon | 
| `--threads` | Number of threads | 
| `--min_length` | Minimum sequence length for common sequence search. |
| `--max_length` | Maximum sequence length for common sequence search. |
---------
Output folders/files
| File/Directory | Description | Generated by |
|---|---|---|
| `<input_directory>/unicycler_fasta/` | Contains the unicycler  FASTA assemblies provided as input. | `plsMD --preprocessing` |
| `<input_directory>/<fasta_file>_plasmid.txt` | ABRICATE output files for each input FASTA file, detailing identified antimicrobial resistance genes. | `plsMD --preprocessing` |
| `<input_directory>/<fasta_file>_PLSDB.txt` | BLASTn output files for each input FASTA file, showing matches against the PLSDB database. | `plsMD --preprocessing` |
| `<input_directory>/plasmid_list.txt` | A list of identified plasmid replicons. | `plsMD --preprocessing` |
| `<genes_directory>/gene_directories/` | Contains the gene directories used in the phylogenetic analysis | `plsMD --process` |
| `<input_directory>/plasmid_files/` | Contains processed plasmid contig files. | `plsMD --processing` |
| `<input_directory>/nonplasmid_files/` | Contains processed non-plasmid contig files. | `plsMD --processing` |
| `<input_directory>/nonreplicon_contigs/` | Contains contigs that does not have a  replicon sequence | `plsMD --processing` |
| `<input_directory>/Processed_ColContigs.txt` | Combined output file for processed Col and Inc plasmids. | `plsMD --processing` |
|  `<input_directory>/plasmid_files/annotation` | Contains Output from AMRfinder, abricate and blastn (IS). | `plsMD --annotation` |
|  `<input_directory>/nonplasmid_files/annotation` | Contains Output from AMRfinder, abricate and blastn (IS). | `plsMD --annotation` |
| `<genes_directory>/phylogenetic_tree/` | Contains the output of the MSA and phylogenetic analysis. | `plsMD --phylogenetics` |
----------
## Citations

If you use plsMD in your research, please cite the following tools:

- **BLAST+** — Camacho et al. (2009) *BMC Bioinformatics* 10:421. https://doi.org/10.1186/1471-2105-10-421
- **Abricate** — Seemann T. GitHub. https://github.com/tseemann/abricate
- **AMRFinderPlus** — Feldgarden et al. (2021) *Sci Rep* 11:12728. https://doi.org/10.1038/s41598-021-91456-0
- **MAFFT** — Katoh & Standley (2013) *Mol Biol Evol* 30:772–780. https://doi.org/10.1093/molbev/mst010
- **IQ-TREE** — Nguyen et al. (2015) *Mol Biol Evol* 32(1):268–274. https://doi.org/10.1093/molbev/msu300
- **ISFinder** — Siguier et al. (2006) *Nucleic Acids Res* 34:D32–D36. https://www.ncbi.nlm.nih.gov/pubmed/16381877 — Database: http://www-is.biotoul.fr
- **mob-suite** — Robertson & Nash (2018) *Microbial Genomics* 4(8):e000206. https://doi.org/10.1099/mgen.0.000206
