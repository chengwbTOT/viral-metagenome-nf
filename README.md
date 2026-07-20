# Viral Metagenome Nextflow Pipeline
A reproducible Nextflow pipeline for viral metagenomic assembly and annotation, ported from the original bash pipeline v0.2.4.

## Overview
This pipeline performs end-to-end viral metagenome analysis for paired-end short-read sequencing data, including quality control, de novo metagenomic assembly, viral sequence identification, taxonomic classification, host prediction and Pfam functional annotation.

## Features
- Read quality control and trimming via metaWRAP `read_qc`
- Metagenomic de novo assembly via metaSPAdes
- Viral sequence identification, binning and host prediction via ViWrap
- Automated result curation: high/medium quality virus filtering, abundance merging, taxonomic and host annotation integration
- Pfam domain functional annotation via HMMER

## Usage
### Basic run
```bash
nextflow run viral_metagenome_pipeline.nf --threads 40
