# A Reproducible Genetic Variant Calling Workflow Across Diverse Species
### Overview
This repository provides a comprehensive workflow, including configuration files and processed outputs, for genetic variant discovery across diverse species. The pipeline implements a standardized, modular, and fully reproducible framework to transform raw Next-Generation Sequencing (NGS) reads into high-quality variant datasets (VCFs). This project aims to advance scientific research by ensuring computational reproducibility, enabling performance benchmarking, and supporting open data sharing.

### Scientific Motivation
This project establishes a transparent pipeline for the genetic variant discovery across diverse species using parameterized workflows and curated datasets. By making all code and results publicly available, we ensure the complete reproducibility of the analyses presented in our associated manuscript.

### Contact
For inquiries regarding analytical methods, results, or technical support, please contact:
-	**HyeonJung Lee** (hyeon@kaist.ac.kr)    
*Korea Advanced Institute of Science and Technology (KAIST)*
-	**Dr. Sunhee Kim** (king@kongju.ac.kr)    
*Kongju National University, South Korea*

### Content List
- [Requirements](#requirements)
  - [Software](#software)
  - [Recommended computing system](#recommended-computing-system)
- [Installation](#installation)
  - [Conda](#conda)
  - [Docker](#docker)
- [Usage](#usage)
  - [Required Arguments](#required-arguments)
  - [Optional Arguments](#optional-arguments)
  - [Sample List File Format](#sample-list-file-format)
- [Examples](#examples-taken-from-run_examplessh)
- [Output Directory Structure](#output-directory-structure)
  - [File Naming Conventions](#file-naming-conventions)
- [Other Details](#other-details)
  - [Threading and Memory](#threading-and-memory)
  - [Using BWA instead of BWA-MEM2](#using-bwa-instead-of-bwa-mem2)
  - [Softlinking (-sl)](#softlinking--sl)
  - [Resumability](#resumability)
- [Links to Example Datasets](#links-to-dataset-examples)

## Requirements

### Software

The [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) package manager is required to manage the pipeline environment. The pipeline and all its dependencies (Python 3, `bwa-mem2`, `samtools`, `picard`, `gatk`) are bundled into a conda environment defined in `pseudoDB_env.yaml`. To find the specific version of Conda compatible with your system, visit the Anaconda Archive at https://repo.anaconda.com/archive/.

Alternatively, the pipeline can also be run via [Docker](https://www.docker.com). Please visit https://docs.docker.com/engine/install for instructions on Docker installation.

### Recommended computing system

The recommended high-performance computing server for the pseudoDB workflow:
- 128 GB RAM
- 1 TB of disk space
- 32 CPU threads @ 2.4 GHz
- A recent stable Unix-based operating system. 

The wall-clock runtime for a typical complete run starting from FASTQ files (e.g., 11.1 GB) is approximately 11.3 hours, including 2.5 hours for sequence alignment, 3.0 hours for pseudoDB construction, 3.1 hours for base quality score recalibration (BQSR), and 2.7 hours for variant calling.

## Installation

### Conda

Clone the repository and `cd` into the main `pseudoDB` directory.

```bash
git clone https://github.com/infoLab204/pseudoDB.git
cd pseudoDB
```

From inside the main `pseudoDB` directory, run `./install.sh` to prepare the conda environment and install the pseudoDB script. The installation directory and conda environment name can be customized with `--install-dir` and `--env-name`.

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

To test the installation, run `./run_examples.sh` from inside the main `pseudoDB` directory.

```bash
./run_examples.sh
```

### Docker

Alternatively, you can also run pseudoDB via Docker. It still uses Conda to install and manage packages, but this approach minimizes dependency conflicts with your host system.

To build the Docker image, run the command below from inside the main `pseudoDB` directory.

```bash
docker build -t pseudodb:latest .
```

To test the installation, run `./run_examples_docker.sh` from inside the main `pseudoDB` directory.

```bash
./run_examples_docker.sh
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
- `-t`, `--threads` — Number of CPU threads passed to compatible tools (default: `8`)
- `-m`, `--memory` - Memory limit (GB) passed to compatible tools (default: `16`)
- `-sl`, `--softlink` — If set, input files are softlinked into the output directory instead of copied.
- `--use-bwa` - Use BWA instead of BWA-MEM2. Takes more time, but consumes less memory. Automatically activated if BWA-MEM2 runs into an error.


### Sample List File Format

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

## Examples (taken from `run_examples.sh`)

This section shows examples with a small human chromosome 22 dataset provided in [example_data/](./example_data).

**Construct a pseudo-database:**
```bash
pseudoDB \
  -sp human_chr22 \
  -fa ./example_data/human_chr22.fa \
  -s ./example_data/human_chr22_sample_list.txt \
  -o ./example_out/human_chr22_case1
```

**Variant calling with dbSNP:**
```bash
pseudoDB \
  -sp human_chr22 \
  -fa ./example_data/human_chr22.fa \
  -s ./example_data/human_chr22_sample_list.txt \
  -o ./example_out/human_chr22_case2 \
  -db ./example_data/human_chr22_dbSNP.vcf.gz \
  -dn dbSNP
```

**Variant calling with a previously generated pseudoDB:**
```bash
pseudoDB \
  -sp human_chr22 \
  -fa ./example_data/human_chr22.fa \
  -s ./example_data/human_chr22_sample_list.txt \
  -o ./example_out/human_chr22_case3 \
  -db ./example_data/human_chr22_pseudoDB.vcf.gz \
  -dn pseudoDB
```

For Docker-based installation, directories where input data files are saved will have to be mounted into the container. All input paths (`-fa`, `-s`, `-db`, etc.) and the output directory (`-o`) must then be under the mounted path so the container can access them. For more information, please visit https://docs.docker.com/engine/storage/bind-mounts.

```bash
docker run \
  -v ./example_data:/data \
  -v ./example_out:/output \
  pseudodb \
  -sp human_chr22 \
  -fa /data/human_chr22.fa \
  -s /data/human_chr22_sample_list_docker.txt \
  -o /output/human_chr22_case1_docker
```

Notes:
- Remember to update the sample list to reflect paths inside the Docker container. For example, instead of `./example_data/data/SAMPLE_1.fastq.gz`, use `/data/SAMPLE_1.fastq.gz` (because `./example_data` is mounted to `/data`).
- Remember to create the output directory (e.g. `./example_out`) before running the Docker container to avoid errors.
- The `--softlink` (`-sl`) option will still work inside the Docker container. The only caveat is that symlinks created inside the container will not resolve correctly from the host after the container exits.
- To use multiple CPU threads, ensure your Docker daemon is configured to allow sufficient CPU resources.

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

### File Naming Conventions

All output files are named using: 
- The FASTQ sample name (`<SN>`)
- The `--species-name` argument (`<SP>`)
- The `--database-name` argument (`<DBN>`)

| Output | Directory | Filename |
|--------|-----------|----------|
| Aligned BAM | `module/align/` | `<SN>_aligned.bam` |
| Recalibrated BAM | `module/machine/` | `<SN>_<DBN>_recalibrated.bam` |
| Variant call VCF | `module/variants/` | `<SP>_<DBN>_variant_calling.vcf.gz` |
| Error rate file | `module/error/` | `<SN>_<DBN>_erate` |
| Quality score file | `module/model/` | `<SN>_<DBN>_qs` |
| PseudoDB VCF | `data/db/` | `<SP>_pseudoDB.vcf.gz` |

## Other Details

### Threading and Memory

The `-t` / `--threads` argument controls the number of CPU threads passed to compatible tools:

- `bwa-mem2 mem` or `bwa mem`
- GATK 3: `UnifiedGenotyper`

The `-m` / `--memory` argument controls the memory limit passed to compatible tools (in GB):

- PICARD: `CreateSequenceDictionary`, `SortSam`, `MarkDuplicates`
- GATK 3: `UnifiedGenotyper`, `BaseRecalibrator`, `PrintReads`

All other tools such are currently not affected by these two arguments and will use their own defaults.

### Using `BWA` instead of `BWA-MEM2`

By default, pseudoDB uses [`BWA-MEM2`](https://github.com/bwa-mem2/bwa-mem2) to index reference FASTA files and align FASTQ files. This is fast, but consumes a lot of memory. 

Using the `--use-bwa` flag allows users to use [`BWA`](https://github.com/lh3/bwa), a slower version of the tool that uses less memory.

Performance on GRCh38_full_analysis_set_plus_decoy_hla.fa indexing and 1 human sample (HG00096; IGSR):

| Step | Tool | Time | RAM | Threads |
|------|------|------|-----|---------|
| FASTA indexing | BWA-MEM2 | ~15 min | ~75 GB | 1 |
| FASTQ alignment | BWA-MEM2 | ~30 min | ~25 GB | 8 |
| FASTA indexing | BWA | ~45 min | ~5 GB | 1 |
| FASTQ alignment | BWA | ~60 min | ~8 GB | 8 |

The default run will also fall back to `BWA` if `BWA-MEM2` fails.

### Softlinking (`-sl`)

By default, all input files (FASTA, database VCF, FASTQ files) are **copied** into the output directory. Passing `-sl` will create **symbolic links** instead, saving disk space and time.

- If symlink creation fails for any reason, the pipeline automatically falls back to copying the file.
- If source and destination resolve to the same absolute path, the operation is skipped silently.
- Softlinks point to the original source paths, so moving or deleting the source files after running will break them.
- **Please be careful with softlinked files — deleting them from the output directory may affect the original source files on some systems.**

### Resumability

The pipeline is designed to be safely re-run on a partially completed output directory:

- **Reference indexing**: Skipped if all index files already exist.
- **Alignment**: Per-sample, skipped if `<SN>_aligned.bam` already exists in `module/align/`.
- **Base recalibration**: Per-sample, skipped if `<SN>_<DBN>_recalibrated.bam` already exists in `module/machine/`.
- **Variant calling** and **error rate**: not resumable at the per-sample level, will re-run in full if called.

## Links to Dataset Examples

This section lists example datasets for human and other species.

1. Human

    - Example FASTQ     : https://www.internationalgenome.org/data-portal/sample (e.g. HG00096) - [Link to download instructions](./legacy/datasets/human.md)
    - Reference FASTA   : http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
    - dbSNP             : https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz
    - pseudoDB          : https://zenodo.org/record/7488070/files/human_pseudoDB.vcf.gz

2. African Oil Palm

    - Example FASTQ     : https://www.ebi.ac.uk/ena/browser/view/PRJEB21246 (e.g. ERR2004436) - [Link to download instructions](./legacy/datasets/african_oil_palm.md)
    - Reference FASTA   : https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/442/705/GCF_000442705.1_EG5/GCF_000442705.1_EG5_genomic.fna.gz
    - dbSNP             : https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_8/by_species/elaeis_guineensis/EG5/51953_GCA_000442705.1_current_ids.vcf.gz
    - pseudoDB          : https://zenodo.org/record/18437285/files/palm_pseudoDB.vcf.gz

3. Brown Bear

    - Example FASTQ     : https://www.ebi.ac.uk/ena/browser/view/PRJNA1139383 (e.g. SRR29938644) - [Link to download instructions](./legacy/datasets/brown_bear.md)
    - Reference FASTA   : https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/065/955/GCF_023065955.2_UrsArc2.0/GCF_023065955.2_UrsArc2.0_genomic.fna.gz
    - dbSNP             : N/A
    - pseudoDB          : https://zenodo.org/record/18437343/files/bear_pseudoDB.vcf.gz

4. Cattle

    - Example FASTQ     : https://www.ebi.ac.uk/ena/browser/view/PRJNA238491 (e.g. SRR1293227) - [Link to download instructions](./legacy/datasets/cattle.md)
    - Reference FASTA   : https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/055/GCF_000003055.6_Bos_taurus_UMD_3.1.1/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fna.gz
    - dbSNP             : https://ftp.ncbi.nih.gov/snp/organisms/archive/cow_9913/VCF/00-All.vcf.gz
    - pseudoDB          : https://zenodo.org/record/18333082/files/cattle_pseudoDB.vcf.gz

5. Chickpea

    - Example FASTQ     : https://db.cngb.org/search/project/CNP0000370/ (e.g. SRR5183095) - [Link to download instructions](./legacy/datasets/chickpea.md)
    - Reference FASTA   : https://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Cicer_arietinum/all_assembly_versions/GCA_000331145.1_ASM33114v1/GCA_000331145.1_ASM33114v1_genomic.fna.gz
    - dbSNP             : https://ftp.ncbi.nih.gov/snp/organisms/archive/chickpea_3827/VCF
    - pseudoDB          : https://zenodo.org/record/7487929/files/chickpea_pseudoDB.vcf.gz

6. Komodo Dragon

    - Example FASTQ     : https://www.ebi.ac.uk/ena/browser/view/PRJNA738464 (e.g. SRR14830854) - [Link to download instructions](./legacy/datasets/komodo_dragon.md)
    - Reference FASTA   : https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/798/865/GCF_004798865.1_ASM479886v1/GCF_004798865.1_ASM479886v1_genomic.fna.gz
    - dbSNP             : N/A
    - pseudoDB          : https://zenodo.org/record/18437361/files/komodo_pseudoDB.vcf.gz

7. Rice

    - Example FASTQ     : https://www.ebi.ac.uk/ena/browser/view/PRJEB6180 (e.g. SAMEA2569416/IRIS_313-10889) - [Link to download instructions](./legacy/datasets/rice.md)
    - Reference FASTA   : https://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Oryza_sativa/all_assembly_versions/GCA_001433935.1_IRGSP-1.0/GCA_001433935.1_IRGSP-1.0_genomic.fna.gz
    - dbSNP             : https://ftp.ncbi.nih.gov/snp/organisms/archive/rice_4530/VCF/00-All.vcf.gz
    - pseudoDB          : https://zenodo.org/record/7488383/files/rice_pseudoDB.vcf.gz

8. Sheep

    - Example FASTQ     : https://www.ebi.ac.uk/ena/browser/view/PRJNA160933 (e.g. SRR501898) - [Link to download instructions](./legacy/datasets/sheep.md)
    - Reference FASTA   : https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Ovis_aries/latest_assembly_versions/GCA_000298735.2_Oar_v4.0/GCA_000298735.2_Oar_v4.0_genomic.fna.gz
    - dbSNP             : https://ftp.ncbi.nih.gov/snp/organisms/archive/sheep_9940/VCF/00-All.vcf.gz
    - pseudoDB          : https://zenodo.org/record/7488425/files/sheep_pseudoDB.vcf.gz

9. Stevia

    - Example FASTQ     : https://www.ebi.ac.uk/ena/browser/view/PRJNA684944 (e.g. SRR13325728) - [Link to download instructions](./legacy/datasets/stevia.md)
    - Reference FASTA   : https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/936/405/GCA_009936405.2_ASM993640v2/GCA_009936405.2_ASM993640v2_genomic.fna.gz
    - dbSNP             : N/A
    - pseudoDB          : https://zenodo.org/record/18437378/files/stevia_pseudoDB.vcf.gz

10. Swan Goose

    - Example FASTQ     : https://www.ebi.ac.uk/ena/browser/view/PRJNA722049 (e.g. SRR14534354) - [Link to download instructions](./legacy/datasets/swan_goose.md)
    - Reference FASTA   : https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/166/845/GCF_002166845.1_GooseV1.0/GCF_002166845.1_GooseV1.0_genomic.fna.gz
    - dbSNP             : https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_8/by_species/anser_cygnoides/GooseV1.0/8845_GCA_002166845.1_current_ids.vcf.gz
    - pseudoDB          : https://zenodo.org/record/18437393/files/goose_pseudoDB.vcf.gz

For more information on how to download these samples, please refer to instructions in [legacy/datasets](./legacy/datasets).