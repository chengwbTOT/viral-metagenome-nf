# Viral Metagenome Assembly & Annotation Nextflow Pipeline

## Overview
A reproducible end-to-end Nextflow pipeline for viral metagenomic analysis of paired-end short-read sequencing data, adapted from the original bash workflow (version 0.2.4). This pipeline covers raw read quality control, de novo metagenomic assembly, viral sequence identification, taxonomic classification, host prediction, and Pfam functional annotation, optimized for freshwater environmental samples.

## Features
- Raw read trimming and quality control via metaWRAP `read_qc` module
- De novo metagenomic assembly using metaSPAdes
- Viral contig identification, binning and quality assessment via ViWrap
- Viral taxonomic classification and genus-level host prediction
- Automated result curation: filtered high/medium-quality viral genomes, normalized abundance table, and integrated summary report
- Pfam protein domain functional annotation via HMMER hmmscan
- Native support for checkpoint resume, batch parallel sample processing, and customizable resource allocation

## Prerequisites
### Required Conda Environments
- `metawrap-env`: For read quality control and metaSPAdes assembly
- `ViWrap`: For viral identification, classification and host prediction
- `hmmer`: For Pfam functional annotation

### Required Databases
- ViWrap reference database
- Pfam-A v37.1 database

## Quick Start

### Input Preparation
Place all paired-end raw FASTQ files in a `raw_data/` directory, following the naming convention:

```text
raw_data/
├── sample1_R1.fastq
└── sample1_R2.fastq
```

### Basic Run
```bash
nextflow run viral_metagenome_pipeline.nf --threads 40
```

### Resume Interrupted Run
Automatically skip all successfully completed steps (recommended for large datasets):
```bash
nextflow run viral_metagenome_pipeline.nf -resume
```

## Parameters

| Parameter | Default Value | Description |
| :--- | :--- | :--- |
| `--threads` | 40 | Number of CPU threads per process |
| `--assembly_memory` | 500 GB | Memory allocated for metaSPAdes assembly |
| `--min_contig_length` | 5000 | Minimum contig length threshold for viral identification |
| `--pfam_evalue` | 1e-5 | E-value cutoff for Pfam domain annotation |
| `--outdir` | results | Root directory for output files |

---

## Output Structure
Results are organized in a standardized directory structure consistent with the original bash pipeline:

```text
results/
└── <sample_id>/
    ├── 02_READ_QC/
    │   ├── final_pure_reads_1.fastq
    │   └── final_pure_reads_2.fastq
    ├── 03_ASSEMBLY/
    │   └── contigs.fasta
    └── 04_ViWrap/
        └── 08_ViWrap_summary_outdir/
            ├── Virus_summary_simplify.txt   # Final integrated summary: quality, abundance, taxonomy, host prediction
            ├── Virus_genomes_files/         # Renamed viral genomes with sample ID prefix
            └── Annotation/
                └── Pfamv37.1/               # Pfam domain annotation outputs
```

---

## Changelog

* **v0.2.4 (2025-03-05)**:
  - Initial Nextflow implementation ported from the original bash pipeline
  - Added support for batch parallel processing of multiple samples
  - Integrated native checkpoint resume functionality
  - Standardized directory structure and configurable parameters

---

## Citation
If you use this pipeline in your research, please cite the underlying tools:
1. Uritskiy GV, DiRuggiero J, Taylor J. MetaWRAP—a flexible pipeline for genome-resolved metagenomic data analysis. *Microbiome*. 2018;6(1):1-13.
2. Zhou Z, et al. ViWrap: A modular pipeline to identify, bin, classify, and predict viral–host relationships for viruses from metagenomes. *iMeta*. 2023;2(2):e118.

---

## Contact
* **Author**: Wen-Bin Cheng
* **Email**: chengwb@mail.ustc.edu.cn
