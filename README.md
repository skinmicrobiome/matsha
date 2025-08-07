# MATSHA

**MATSHA** (Microbial Association with Traits from SHotgun Analysis) is a streamlined command-line tool and Python package for detecting genotype-phenotype associations using shotgun metagenomic data.

![Logo](./matsha_logo.png)

## ✨ Features

- 🔍 Supports both single-end and paired-end reads
- 📊 Depth-aware subsampling to reduce confounding
- 🧬 Variant-level and gene-level GWAS
- 🧪 Tested with toy datasets and CI workflows
- ✅ Easy to install and use with a clean CLI interface

## 🚀 Quick Start

### Installation

```bash
git clone https://github.com/skinmicrobiome/matsha.git
cd matsha
conda env create -f environment.yml
conda activate matsha_env
````

Now you can run `matsha` from the command line.

### Example usage

```bash
matsha --input tests/toy_data/input_file.tsv \
       --genome tests/toy_data/reference.fna \
       --genbank tests/toy_data/reference.gbff \
       --output tests/toy_data/output/ \
       --paired \
       --threads 4
```

For full options, run:

```bash
matsha --help
```

## 📁 Input Format

The `input_file.tsv` should be a tab-delimited file with a header row. For **paired-end reads**, the first two columns must be paths to R1 and R2 FASTQ files. For **single-end reads**, only R1 is required. Columns from the third onward can include one or more phenotypes for association testing.

- Both compressed (`.fastq.gz`) and uncompressed (`.fastq`) FASTQ files are supported.
- Phenotypes can be **binary** (e.g., disease status) or **quantitativ** (e.g., BMI). 
- Missing values should be left **blank**. Samples with missing values will be excluded from the analysis for that specific phenotype.

Example:

```tsv
R1                              R2                              disease status   BMI
/path/to/sample1_R1.fastq.gz    /path/to/sample1_R2.fastq.gz    1            27.4
/path/to/sample2_R1.fastq.gz    /path/to/sample2_R2.fastq.gz    0            30.1
/path/to/sample3_R1.fastq.gz    /path/to/sample3_R2.fastq.gz    1            
/path/to/sample4_R1.fastq.gz    /path/to/sample4_R2.fastq.gz    0            22.8
````

## 📁 Outputs
After running `matsha`, the following output structure is generated in the specified output directory:
```text
output_dir/
├── logfile
├── phenotype_subsampled.txt
├── all/
│ ├── combined.g.vcf
│ ├── combined.g.annotated.vcf
│ ├── gwas_output.gene.tsv
│ ├── gwas_output.variant.tsv
│ ├── gwas_output_significant.gene.tsv
│ └── gwas_output_significant.variant.tsv
├── 0/
│ └── (same structure as all/)
├── 1/
└── ...
```
| Path | Description |
|------|-------------|
| `logfile` | A text file recording the runtime log of `matsha`, including command-line arguments, processing steps, and any warnings or errors. Useful for debugging or reproducibility. |
| `phenotype_subsampled.txt` | A tab-delimited file listing the phenotypes for which subsampling was performed. Each listed phenotype corresponds to a numbered folder (e.g., `0/`, `1/`, ...). |
| `all/` | Folder containing results from the full dataset without any subsampling. |
| `0/`, `1/`, ... | Folders for each subsampled phenotype group listed in `phenotype_subsampled.txt`. Each folder contains results specific to that phenotype. |
| `*/combined.g.vcf` | The unannotated merged VCF file containing variant calls across samples. |
| `*/combined.g.annotated.vcf` | The same VCF file, annotated using SnpEff to include gene IDs and predicted impacts. |
| `*/gwas_output.gene.tsv` | All GWAS results for gene-level analysis. |
| `*/gwas_output.variant.tsv` | All GWAS results for variant-level analysis. |
| `*/gwas_output_significant.gene.tsv` | Significantly associated genes based on statistical thresholds. |
| `*/gwas_output_significant.variant.tsv` | Significantly associated variants based on statistical thresholds. |

## 🧪 Testing

To run unit tests using the toy dataset:

```bash
pytest
```

## 🔧 Development

- CLI implemented using [`click`](https://click.palletsprojects.com/)
- Modular structure:
  - `core.py`: main pipeline logic (mapping, coverage, variant calling, etc.)
  - `utils.py`: helper functions
  - `gwas.py`: GWAS-specific functions
- Toy data under `tests/toy_data/` for reproducibility

## 📝 Citation

Coming soon.

## 📜 License

MIT License © 2025 Zeyang Shen
