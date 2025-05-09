plsMD: A plasmid reconstruction tool for short-read assemblies
-----
Introduction
------------
Whole genome sequencing (WGS) has become a cornerstone of antimicrobial resistance (AMR) surveillance. However, reconstructing plasmid sequences from short-read WGS data is challenging due to repetitive sequences and assembly fragmentation. 
plsMD is a tool designed for comprehensive plasmid reconstruction from short-read assemblies, going beyond simple contig binning. It utilizes Unicycler assemblies with established plasmid databases to reconstruct full plasmid sequences, plsMD offers both single and multi-sample analysis, enabling researchers to perform detailed plasmid characterization and phylogenetic investigations.

-------
Installation
-----------

```bash
git clone https://github.com/Genomics-and-Metagenomics-Unit-57357/plsMD.git
cd plsMD
```
Build the Docker Image (without downloading PLSDB database)
----------
```bash
docker build -t plsmd .
```
Build the Docker Image (with downloading PLSDB database)
--------
```bash
docker build --build-arg DOWNLOAD_DB=true -t plsmd .
```
Set Up Command Line Access
---------
```bash
# Make the wrapper executable
chmod +x plsMD_wrapper.sh

# Create symlink 
sudo ln -sf $(pwd)/plsMD_wrapper.sh /usr/local/bin/plsMD
```
Verify Installation
---------
```bash
plsMD --version
plsMD --help
```
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
| `--min_length` | Maximum sequence length for common sequence search. |
---------
Output folders/files
| File/Directory | Description | Generated by |
|---|---|---|
| `<input_directory>/unicycler_fasta/` | Contains the unicycler  FASTA assemblies provided as input. | `plsMD --preprocessing` |
| `<input_directory>/<fasta_file>_plasmid.tct` | ABRICATE output files for each input FASTA file, detailing identified antimicrobial resistance genes. | `plsMD --preprocessing` |
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
