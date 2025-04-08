# ATACseq pipeline
**Ulises Gaspar-Eslava**

Experimental design: 3 ATAC seq biological replicates from sorted spleen pDCs for 2 conditions: uninfected and LCMV CL13 chronically infected mice.

Unix dependencies: 
- FastQC
- Cutadapt
- MultiQC
- Bowtie2
- Samtools
- bedtools
- MACS3
- HOMER

## Quality Control

QC raw reads: 
```bash
fastqc CL13-1-ATAC.fastq
fastqc CL13-2-ATAC.fastq
fastqc CL13-3-ATAC.fastq
fastqc NAIVE-1-ATAC.fastq
fastqc NAIVE-2-ATAC.fastq
fastqc NAIVE-3-ATAC.fastq
```
On th
Trim adapter sequences:
```bash
cutadapt -a CTGTCTCTTATACACATCTGACGCTGCCGACGA -o TRIMMED/CL13-1-ATAC_trimmed.fastq RAW_READS/CL13-1-ATAC.fastq --cores 10
cutadapt -a CTGTCTCTTATACACATCTGACGCTGCCGACGA -o TRIMMED/CL13-2-ATAC_trimmed.fastq RAW_READS/CL13-2-ATAC.fastq --cores 10
cutadapt -a CTGTCTCTTATACACATCTGACGCTGCCGACGA -o TRIMMED/CL13-3-ATAC_trimmed.fastq RAW_READS/CL13-3-ATAC.fastq --cores 10
cutadapt -a CTGTCTCTTATACACATCTGACGCTGCCGACGA -o TRIMMED/NAIVE-1-ATAC_trimmed.fastq RAW_READS/NAIVE-1-ATAC.fastq --cores 10
cutadapt -a CTGTCTCTTATACACATCTGACGCTGCCGACGA -o TRIMMED/NAIVE-2-ATAC_trimmed.fastq RAW_READS/NAIVE-2-ATAC.fastq --cores 10
cutadapt -a CTGTCTCTTATACACATCTGACGCTGCCGACGA -o TRIMMED/NAIVE-3-ATAC_trimmed.fastq RAW_READS/NAIVE-3-ATAC.fastq --cores 10
```
Specific adapter sequences can be retreieved from Illumina's [Sequences for Nextera, Illumina Prep, and Illumina PCR Kits](https://support-docs.illumina.com/SHARE/AdapterSequences/Content/SHARE/AdapterSeq/Nextera/SequencesNextera_Illumina.htm). In this case, the reverse complement of the Nextera Transposase Adapters was used. 

QC trimmed reads:
```bash
fastqc CL13-1-ATAC_trimmed.fastq
fastqc CL13-2-ATAC_trimmed.fastq
fastqc CL13-3-ATAC_trimmed.fastq
fastqc NAIVE-1-ATAC_trimmed.fastq
fastqc NAIVE-2-ATAC_trimmed.fastq
fastqc NAIVE-3-ATAC_trimmed.fastq
```
Combine QC data:
```bash
multiqc RAW_READS/
multiqc TRIMMED/
```
This has to be run on both raw and trimmed directories. 

## Alignment
Indexed reference genome can be retrieved from [iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html). For this project, the mouse genome from UCSC build mm10 was used.
```bash
(bowtie2 --no-unal --very-sensitive -p 30 -x REFERENCE/Mus_muscul/UCSC/mm10/Sequence/Bowtie2Index/genome -U TRIMMED/CL13-1-ATAC_trimmed.fastq -S ALIGNMENTS/CL13-1-ATAC.sam) 2>CL13-1-ATAC.log
(bowtie2 --no-unal --very-sensitive -p 30 -x REFERENCE/Mus_muscul/UCSC/mm10/Sequence/Bowtie2Index/genome -U TRIMMED/CL13-2-ATAC_trimmed.fastq -S ALIGNMENTS/CL13-2-ATAC.sam) 2>CL13-2-ATAC.log
(bowtie2 --no-unal --very-sensitive -p 30 -x REFERENCE/Mus_muscul/UCSC/mm10/Sequence/Bowtie2Index/genome -U TRIMMED/CL13-3-ATAC_trimmed.fastq -S ALIGNMENTS/CL13-3-ATAC.sam) 2>CL13-3-ATAC.log
(bowtie2 --no-unal --very-sensitive -p 30 -x REFERENCE/Mus_muscul/UCSC/mm10/Sequence/Bowtie2Index/genome -U TRIMMED/NAIVE-1-ATAC_trimmed.fastq -S ALIGNMENTS/NAIVE-1-ATAC.sam) 2>NAIVE-1-ATAC.log
(bowtie2 --no-unal --very-sensitive -p 30 -x REFERENCE/Mus_muscul/UCSC/mm10/Sequence/Bowtie2Index/genome -U TRIMMED/NAIVE-2-ATAC_trimmed.fastq -S ALIGNMENTS/NAIVE-2-ATAC.sam) 2>NAIVE-2-ATAC.log
(bowtie2 --no-unal --very-sensitive -p 30 -x REFERENCE/Mus_muscul/UCSC/mm10/Sequence/Bowtie2Index/genome -U TRIMMED/NAIVE-3-ATAC_trimmed.fastq -S ALIGNMENTS/NAIVE-3-ATAC.sam) 2>NAIVE-2-ATAC.log
```



