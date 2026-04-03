<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-funomic_logo_dark.png">
    <img alt="nf-core/funomic" src="docs/images/nf-core-funomic_logo_light.png">
  </picture>
</h1>

## Introduction

A reproducible and scalable Nextflow implementation of the FunOMIC2 pipeline for fungal metagenomic analysis.

This project is currently a work in progress and does not yet meet full nf-core standards. At present, only the taxonomic profiling module has been implemented. Contributions—especially toward completing and improving the functional profiling component—are highly encouraged.

---

## Background

FunOMIC2 is a fungal metagenomics pipeline developed by the Manichanh Lab, designed for taxonomic and functional profiling of fungal communities from metagenomic sequencing data.

Original repository:  
https://github.com/ManichanhLab/FunOMIC2

This project does not reimplement the core algorithms. Instead, it reorganizes and orchestrates the original pipeline within a modern workflow framework using Nextflow.

---

## Limitations

- Modifications to the original script were unavoidable based on methodological considerations described in the referenced study: https://link.springer.com/article/10.1186/s40168-023-01693-w#Sec2  
    As a result, parts of the implementation may diverge from the original design.
    
- During refactoring, several tools and dependencies were upgraded to newer versions, which may lead to differences in behavior or outputs compared to the original pipeline.
    
- Logical errors may have been introduced during adaptation and restructuring. Results should be interpreted with caution.
    
- Functional profiling components were not modified and remain as implemented in the original workflow.
    
- The codebase is not yet fully modularized. Some components remain tightly coupled, which may limit maintainability and extensibility.
    
- Additional validation and benchmarking are required to ensure robustness across diverse datasets and environments.
    

---
## Requirements

- Nextflow (>= 22.x recommended)
- Java (required by Nextflow)
- Docker or Singularity (recommended for reproducibility)

---

## Usage

## 1. Database Setup

## UHGG Database

```sh
#!/usr/bin/env bash
set -euo pipefail

base='https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/species_catalogue'

mkdir -p ./fasta
mkdir -p ../index

curl -fsSL "$base/" \
  | grep -Eo 'MGYG[0-9]+/' \
  | sort -u \
  | xargs -I{} -P 10 sh -c '
      curl -fsSL "$1/$2" | grep -Eo "MGYG[0-9]+/" | sort -u | while read -r species; do
          printf "%s/%s%sgenome/%s.fna\n" "$1" "$2" "$species" "${species%/}"
      done
  ' _ "$base" {} \
  | sort -u > urls.txt

url_count=$(wc -l < urls.txt)

aria2c -i urls.txt --check-certificate=false -x 16 -j 5 -c --dir=./fasta


find ./fasta -name '*.fna' -print0 \
  | sort -z \
  | xargs -0 cat > ./uhgg_catalogue_v2.0.2.fna


bowtie2-build --large-index --threads 16 ./uhgg_catalogue_v2.0.2.fna uhgg_v2.0.2

mv *.bt2l ../index/

```

## FunOMIC-T Database

```
wget https://manichanh.vhir.org/funomic/FunOMIC-T.v2.tar.xz -O FunOMIC-T.v2.tar.xz
```

---

## 2. Configuration

## params.yaml
```yaml             
input  : 'samplesheet.csv'
outdir : 'results'

bacterial_db : '/database/bowtie2/uhgg/index'
taxonomic_db : '/database/FunOMIC/FunOMIC-T.v2'
protein_db : '/database/FunOMIC/FunOMIC.P.v2'

coverage                   : 80.0
min_read_length            : 60
```
---
### samplesheet.csv

```csv
sample,fastq_1,fastq_2
S1,/data/S1_R1.fastq.gz,/data/S1_R2.fastq.gz
S2,/data/S2_R1.fastq.gz,/data/S2_R2.fastq.gz
```

---
## 3. Run Pipeline

```sh
nextflow run "$PIPELINE_DIR/main.nf" \
   -params-file params.yaml \
   -profile docker \
   -resume
```
---
## Credits

This project is based on the original FunOMIC2 pipeline developed by the Manichanh Lab.

---

## Citations

An extensive list of references for the tools used by the pipeline can be found in `CITATIONS.md`.

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

> Xie, Z. and Manichanh, C. (2022). FunOMIC: Pipeline with built-in Fungal Taxonomic and Functional Databases for Human Mycobiome Profiling. [https://doi.org/10.1016/j.csbj.2022.07.010](https://doi.org/10.1016/j.csbj.2022.07.010) 
> Vega-Abellaneda, S., Xie, Z. and Manichanh, C. (to be published). FunOMIC2: Updated pipeline with built-in fungal taxonomic and functional databases for human mycobiome profiling.

