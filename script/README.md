## Script Usage

```bash
pseudoDB.py -sp <species_name> -fa <fasta> -s <sample_list> -o <output_dir> [options]
```

### Required Arguments

| Flag | Long | Description |
|------|------|-------------|
| `-sp` | `--species` | Species name used in output file naming |
| `-fa` | `--fasta` | Path to reference FASTA file (`.fa`, `.fna`, `.fasta`, optionally `.gz`) |
| `-s` | `--sample-list` | Path to file listing FASTQ input paths (one per line) |
| `-o` | `--output-dir` | Path to root output directory (created if it doesn't exist) |

### Optional Arguments

| Flag | Long | Default | Description |
|------|------|---------|-------------|
| `-db` | `--database` | `None` | Path to database VCF (`.vcf`, `.vcf.gz`, `.bcf`, `.bcf.gz`). If omitted, pipeline runs in pseudo-database construction mode |
| `-dn` | `--database-name` | Derived from `-db` filename | Name used in output file naming. If not provided, derived by stripping VCF extension from the database filename |
| `-t` | `--threads` | `1` | Number of CPU threads passed to `bwa-mem2 mem` |
| `-sl` | `--softlink` | `False` | If set, input files are softlinked into the output directory instead of copied |

---

## Examples

**Construct a pseudo-database:**
```bash
pseudoDB.py \
  -sp human \
  -fa /refs/GRCh38.fa \
  -s /data/sample_list.txt \
  -o /results/human_run
```

**Variant calling with dbSNP:**
```bash
pseudoDB.py \
  -sp human \
  -fa /refs/GRCh38.fa \
  -s /data/sample_list.txt \
  -o /results/human_run \
  -db /data/db/dbSNP_v155.vcf.gz \
  -dn dbSNP
```

**Variant calling with a previously generated pseudoDB:**
```bash
pseudoDB.py \
  -sp human \
  -fa /refs/GRCh38.fa \
  -s /data/sample_list.txt \
  -o /results/human_run2 \
  -db /results/human_run/data/db/human_pseudoDB.vcf.gz
  -dn pseudoDB
```

---

## Sample List File Format

The file passed to `-s` must contain accessible paths to FASTQ files, one per line. Files must follow the naming pattern:

```
<sample_name>_1.<ext>
<sample_name>_2.<ext>
```

Where `<ext>` is one of: `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`.

Both R1 and R2 files must be present for every sample. The pipeline will raise an error if any pair is incomplete.

Example:
```
/data/samples/HG00096_1.fastq.gz
/data/samples/HG00096_2.fastq.gz
/data/samples/HG00097_1.fastq.gz
/data/samples/HG00097_2.fastq.gz
```

---

## Output Directory Structure

The pipeline creates and populates the following directory tree under `--output-dir`:

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

---

## File Naming Conventions

All output files are named using `<species_name>` (from `-sp`) and `<db_name>` (from `-dn` or derived from `-db`).

| File | Location | Pattern |
|------|----------|---------|
| Aligned BAM | `module/align/` | `<sample>_aligned.bam` |
| Recalibrated BAM | `module/machine/` | `<sample>_<db_name>_recalibrated.bam` |
| Variant call VCF | `module/variants/` | `<species>_<db_name>_variant_calling.vcf.gz` |
| Error rate file | `module/error/` | `<sample>_<db_name>_erate` |
| Quality score file | `module/model/` | `<sample>_<db_name>_qs` |
| PseudoDB VCF | `data/db/` | `<species>_pseudoDB.vcf.gz` |

---

## Pipeline Modes

### Mode 1: Pseudo-Database Construction (`-db` not provided)

Used to generate a pseudo-database of variants from your own samples.

**Pipeline steps:**
1. Set up output directory and copy/link input files
2. Index reference FASTA
3. Align each sample to reference
4. Call variants across all samples

### Mode 2: Variant Calling with Existing Database (`-db` provided)

Used for variant calling with a known variants database such as dbSNP or a previously generated pseudoDB.

**Pipeline steps:**
1. Set up output directory and copy/link input files
2. Index reference FASTA
3. Align each sample to reference
4. Recalibrate base quality scores
5. Call variants across all samples
6. Estimate sample error rate
7. Estimate model-adjusted base quality scores

---

## Softlinking (`-sl`)

By default, all input files (FASTA, database VCF, FASTQ files) are **copied** into the output directory. Passing `-sl` will create **symbolic links** instead, saving disk space.

- If symlink creation fails for any reason, the pipeline automatically falls back to copying the file.
- If source and destination resolve to the same absolute path, the operation is skipped silently.
- Softlinks point to the original source paths, so moving or deleting the source files after running the pipeline will break them.
- **Please be careful with softlink deletion as the original source file may be affected**

---

## Resumability

The pipeline is designed to be safely re-run on a partially completed output directory:

- **Reference indexing**: skipped if all index files already exist.
- **Alignment**: per-sample, skipped if `<sample>_aligned.bam` already exists in `module/align/`.
- **Base recalibration**: per-sample, skipped if `<sample>_<db_name>_recalibrated.bam` already exists in `module/machine/`.
- **Variant calling** and **error rate**: not resumable at the per-sample level — will re-run in full if called.

---

## Threading

The `-t` / `--threads` argument controls the number of CPU threads passed to `bwa-mem2 mem` only. All other tools (`picard`, `gatk`, `samtools`) are currently not parallelized by this flag and will use their own defaults.