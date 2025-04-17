# ATACseq

This repository provides a step-by-step computational workflow for analyzing ATAC seq data, for use in the Zuniga lab. 

## Contents
- ATACseq_pipeline.md: Shell workflow - from reads to peaks
- DifferentialAnalysis.md: R workflow - peaks to differentially accessible regions

## Requirements: 
**Unix dependencies**
- FastQC (Version 0.12.1)
- Cutadapt (Version 4.5)
- MultiQC (Version 1.26)
- Bowtie2 (Version 2.5.2)
- Samtools (Version 1.20)
- bedtools (Version 3.31.1)
- MACS3 (3.0.2)

**R packages**
(Using R 4.1)
- csaw
- GenomicRanges
- ChIPpeakAnni
- edgeR
- ggrepel
- ggrastr
- DiffBind
