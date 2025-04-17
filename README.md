# ATACseq

This repository provides a step-by-step computational workflow for analyzing ATAC seq data, for use in the Zuniga lab. 

## Contents
- ATACseq_pipeline.md: Shell workflow - from reads to peaks
- DifferentialAnalysis.md: R workflow - peaks to differentially accessible regions

## Requirements: 
**Unix dependencies**
- FastQC (Version 0.12.1)
- Cutadapt (Version 4.5)
- MultiQC
- Bowtie2
- Samtools
- bedtools
- MACS3

**R packages**
- csaw
- GenomicRanges
- ChIPpeakAnni
- edgeR
- ggrepel
- ggrastr
- DiffBind
