# A Reproducible Genetic Variant Calling Workflow Across Diverse Species
### Overview
This repository provides a comprehensive workflow, including configuration files and processed outputs, for genetic variant discovery across diverse species. The pipeline implements a standardized, modular, and fully reproducible framework to transform raw Next-Generation Sequencing (NGS) reads into high-quality variant datasets (VCFs). This project aims to advance scientific research by ensuring computational reproducibility, enabling performance benchmarking, and supporting open data sharing.

### Scientific Motivation
This project establishes a transparent pipeline for the genetic variant discovery across diverse species using parameterized workflows and curated datasets. By making all code and results publicly available, we ensure the complete reproducibility of the analyses presented in our associated manuscript.

### Reproducibility
To facilitate ease of use and verification, this repository includes :
-	**Installation Guide**: Comprehensive instructions for setting up the required software environment.
-	**Workflow Scripts**: Fully documented scripts with pre-defined parameters.
-	**Validation Dataset**: Example datasets provided for testing and pipeline verification.

### Contact
For inquiries regarding analytical methods, results, or technical support, please contact:
-	**HyeonJung Lee** (hyeon@kaist.ac.kr)    
*Korea Advanced Institute of Science and Technology (KAIST)*
-	**Dr. Sunhee Kim** (king@kongju.ac.kr)    
*Kongju National University, South Korea*

## Requirements

The [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) package manager is required to manage the pipeline environment. The pipeline and all its dependencies (Python 3, `bwa-mem2`, `samtools`, `picard`, `gatk`) are bundled into a conda environment defined in `pseudoDB_env.yaml'. To find the specific version of Conda compatible with your system, visit the Anaconda Archive at https://repo.anaconda.com/archive/.

[Docker](https://www.docker.com) support is planned as an alternative to conda-based installation. This will allow the pipeline to be run in a fully self-contained container without requiring conda.

## Installation

Run `./install.sh` to prepare the conda environment and install the pseudoDB script. 

The installation directory and conda environment name can be customized with `--install-dir` and `--env-name`.

```bash
# Default installation:
# - User-level installation in $HOME/.local/bin
# - Creates a conda environment named 'gatk3'.
./install.sh

# Customized installation:
./install.sh --install-dir=<INSTALL_DIR> --env-name=<ENV_NAME>
```

After installation, activate the conda environment before running the pipeline:

```bash
conda activate gatk3  # Default conda environment name
conda activate <ENV_NAME>
```

To test the installation, run `./run_examples.sh`.

```bash
./run_examples.sh
```

## Usage

```bash
pseudoDB --species <SPECIES> --fasta <FASTA> --sample-list <SAMPLE_LIST> -output-dir <OUTPUT_DIR> [--database <DATABASE>] [--database-name <DATABASE_NAME>] [--threads <THREADS>] [--softlink]
```

### Required Arguments

- `-sp`, `--species` — Species name used in output file naming
- `-fa`, `--fasta` — Path to reference FASTA file (`.fa`, `.fna`, `.fasta`, optionally `.gz`)
- `-s`, `--sample-list` — Path to file listing paired-end FASTQ input paths (one per line)
- `-o`, `--output-dir` — Path to root output directory

### Optional Arguments

- `-db`, `--database` — Path to database VCF (`.vcf`, `.vcf.gz`). If omitted, the pipeline will run in pseudo-database construction mode
- `-dn`, `--database-name` — Name used in output file naming. If not provided, will be derived by from the database filename
- `-t`, `--threads` — Number of CPU threads passed to `bwa-mem2 mem` (default: `1`)
- `-sl`, `--softlink` — If set, input files are softlinked into the output directory instead of copied

## Examples (taken from `run_examples.sh`)

**Construct a pseudo-database:**
```bash
pseudoDB \
  --species human_chr22 \
  --fasta ./example_data/chr22.fa \
  --sample-list ./example_data/sample_list_chr22.txt \
  --output-dir ./example_out/human_chr22_case1
```

**Variant calling with dbSNP:**
```bash
pseudoDB \
  --species human_chr22 \
  --fasta ./example_data/chr22.fa \
  --sample-list ./example_data/sample_list_chr22.txt \
  --output-dir ./example_out/human_chr22_case2 \
  --database ./example_data/chr22_dbSNP.vcf.gz \
  --database-name dbSNP
```

**Variant calling with a previously generated pseudoDB:**
```bash
pseudoDB \
  --species human_chr22 \
  --fasta ./example_data/chr22.fa \
  --sample-list ./example_data/sample_list_chr22.txt \
  --output-dir ./example_out/human_chr22_case3 \
  --database ./example_data/chr22_pseudoDB.vcf.gz \
  --database-name pseudoDB
```

## Sample List File Format

The file passed to `-s` must contain accessible paths to paired-end FASTQ files, one per line. Files must follow the naming pattern:

```
<SAMPLE_NAME>_1.<EXT>
<SAMPLE_NAME>_2.<EXT>
```

Where `<EXT>` is one of: `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`.

Both R1 and R2 files must be present for every sample. The pipeline will raise an error if any pair is incomplete.

Example:
```
/data/samples/HG00096_1.fastq.gz
/data/samples/HG00096_2.fastq.gz
/data/samples/HG00097_1.fastq.gz
/data/samples/HG00097_2.fastq.gz
```

## Output Directory Structure

<img width="800" height="350" alt="image" src="https://github.com/user-attachments/assets/38e5a96e-2a8c-49be-bc23-f557ea28ded6" />>  <br>
*Fig. 1 : The overall structure of the directories. Here, `--output_dir` is set to \<name of species\>.*

The pipeline creates and populates the following directory tree under `--output-dir`.

User input files, such as reference FASTA, database VCF, and FASTQ files with be copied into the following folders. If the `--softlink` option is used, they will be softlinked instead of copied to save time.

The pipeline output files will be placed in directories under `module/` except for the pseudoDB VCF output, which will be placed in `data/db/`.

```
output_dir/
|-- data/
|   |-- db/        # Database VCF files and pseudoDB output
|   |-- fastq/     # Input FASTQ files (copied or softlinked)
|   `-- ref/       # Reference FASTA and all index files
`-- module/
    |-- align/     # Aligned BAM files (_aligned.bam) and temp files
    |-- error/     # Error rate estimation output (_erate files)
    |-- machine/   # Recalibrated BAM files (_recalibrated.bam)
    |-- model/     # Model-adjusted quality score output (_qs files)
    `-- variants/  # Variant calling VCF output
```

## File Naming Conventions

All output files are named using `<species_name>` (from `-sp`) and `<db_name>` (from `-dn` or derived from `-db`).

- Aligned BAM: `module/align/<SAMPLE_NAME>_aligned.bam`
- Recalibrated BAM: `module/machine/<SAMPLE_NAME>_<DB_NAME>_recalibrated.bam`
- Variant call VCF: `module/variants/<SPECIES>_<DB_NAME>_variant_calling.vcf.gz`
- Error rate file: `module/error/<SAMPLE_NAME>_<DB_NAME>_erate`
- Quality score file: `module/model/<SAMPLE_NAME>_<DB_NAME>_qs`
- PseudoDB VCF: `data/db/<SPECIES>_pseudoDB.vcf.gz`

## Softlinking (`-sl`)

By default, all input files (FASTA, database VCF, FASTQ files) are **copied** into the output directory. Passing `-sl` will create **symbolic links** instead, saving disk space and time.

- If symlink creation fails for any reason, the pipeline automatically falls back to copying the file.
- If source and destination resolve to the same absolute path, the operation is skipped silently.
- Softlinks point to the original source paths, so moving or deleting the source files after running will break them.
- **Please be careful with softlinked files — deleting them from the output directory may affect the original source files on some systems.**

## Resumability

The pipeline is designed to be safely re-run on a partially completed output directory:

- **Reference indexing**: skipped if all index files already exist.
- **Alignment**: per-sample, skipped if `<SAMPLE_NAME>_aligned.bam` already exists in `module/align/`.
- **Base recalibration**: per-sample, skipped if `<SAMPLE_NAME>_<DB_NAME>_recalibrated.bam` already exists in `module/machine/`.
- **Variant calling** and **error rate**: not resumable at the per-sample level — will re-run in full if called.

## Threading

The `-t` / `--threads` argument controls the number of CPU threads passed to `bwa-mem2 mem` only. All other tools are currently not affected by this flag and will use their own defaults.