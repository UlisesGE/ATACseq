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

Aditional scripts used:
SamtoolsFiltering (https://github.com/UlisesGE/ATACseq/blob/main/SamtoolsFiltering.sh)

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
This has to be run on the raw file's directories. 

Trim adapter sequences:
```bash
cutadapt -a CTGTCTCTTATACACATCTGACGCTGCCGACGA -o TRIMMED/CL13-1-ATAC_trimmed.fastq RAW_READS/CL13-1-ATAC.fastq --cores 10
cutadapt -a CTGTCTCTTATACACATCTGACGCTGCCGACGA -o TRIMMED/CL13-2-ATAC_trimmed.fastq RAW_READS/CL13-2-ATAC.fastq --cores 10
cutadapt -a CTGTCTCTTATACACATCTGACGCTGCCGACGA -o TRIMMED/CL13-3-ATAC_trimmed.fastq RAW_READS/CL13-3-ATAC.fastq --cores 10
cutadapt -a CTGTCTCTTATACACATCTGACGCTGCCGACGA -o TRIMMED/NAIVE-1-ATAC_trimmed.fastq RAW_READS/NAIVE-1-ATAC.fastq --cores 10
cutadapt -a CTGTCTCTTATACACATCTGACGCTGCCGACGA -o TRIMMED/NAIVE-2-ATAC_trimmed.fastq RAW_READS/NAIVE-2-ATAC.fastq --cores 10
cutadapt -a CTGTCTCTTATACACATCTGACGCTGCCGACGA -o TRIMMED/NAIVE-3-ATAC_trimmed.fastq RAW_READS/NAIVE-3-ATAC.fastq --cores 10
```
Specific adapter sequences can be retreieved from Illumina's [Sequences for Nextera, Illumina Prep, and Illumina PCR Kits](https://support-docs.illumina.com/SHARE/AdapterSequences/Content/SHARE/AdapterSeq/Nextera/SequencesNextera_Illumina.htm). In this case, the **reverse complement** of the Nextera Transposase Adapters was used. 

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
This has to be run on both raw and trimmed directories to combine both reports. 

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
### Alignment adjustments and filtering
It is necessary to edit the paths to alignment files and to output directory for filtered alignments. This script will sort, index, remove misaligned reads, mitochondrial reads, and PCR duplicates. 
```bash
chmod +x SamtoolsFiltering.sh
./SamtoolsFiltering.sh
```
The resultig filtered alignment and index files will be in a new directory, with the suffix *-v2.bam.
## Visualize alignment coverage and pileups
```bash
bamCoverage -b ALIGNMENTS/CL13-1-ATAC-v2.bam -o CL13-1-ATAC-coverage.bigwig -bs 1 -p max/2
bamCoverage -b ALIGNMENTS/CL13-2-ATAC-v2.bam -o CL13-2-ATAC-coverage.bigwig -bs 1 -p max/2
bamCoverage -b ALIGNMENTS/CL13-3-ATAC-v2.bam -o CL13-3-ATAC-coverage.bigwig -bs 1 -p max/2
bamCoverage -b ALIGNMENTS/NAIVE-1-ATAC-v2.bam -o NAIVE-1-ATAC-coverage.bigwig -bs 1 -p max/2
bamCoverage -b ALIGNMENTS/NAIVE-2-ATAC-v2.bam -o NAIVE-2-ATAC-coverage.bigwig -bs 1 -p max/2
bamCoverage -b ALIGNMENTS/NAIVE-3-ATAC-v2.bam -o NAIVE-3-ATAC-coverage.bigwig -bs 1 -p max/2
```
These bigwig files can then be used to visualize signal distribution for the sequencing reads using [IGV](https://igv.org/app/) or [WashU Epigenome Browser](https://epigenomegateway.wustl.edu/legacy/). 
## Peak calling
Call broad signal peaks on each replicate - FDR 0.01 as recommended in the [program's documentation](https://macs3-project.github.io/MACS/docs/callpeak.html). The BAM file must be first converted to a MACS3 compatible BED file.
```bash
macs3 randsample -i ALIGNMENTS/CL13-1-ATAC-v2.bam -f BAM -p 100 -o BED/CL13-1-ATAC.bed
macs3 callpeak -f BED --nomodel --shift -37 --extsize 73 -B --broad -g 2652783500 --keep-dup all -q 0.01 -n CL13-1-ATAC -t CL13-1-ATAC.bed --outdir ./PEAKS 2>> macs3.log
macs3 randsample -i ALIGNMENTS/CL13-2-ATAC-v2.bam -f BAM -p 100 -o BED/CL13-2-ATAC.bed
macs3 callpeak -f BED --nomodel --shift -37 --extsize 73 -B --broad -g 2652783500 --keep-dup all -q 0.01 -n CL13-2-ATAC -t CL13-2-ATAC.bed --outdir ./PEAKS 2>> macs3.log
macs3 randsample -i ALIGNMENTS/CL13-3-ATAC-v2.bam -f BAM -p 100 -o BED/CL13-3-ATAC.bed
macs3 callpeak -f BED --nomodel --shift -37 --extsize 73 -B --broad -g 2652783500 --keep-dup all -q 0.01 -n CL13-3-ATAC -t CL13-3-ATAC.bed --outdir ./PEAKS 2>> macs3.log
macs3 randsample -i ALIGNMENTS/NAIVE-1-ATAC-v2.bam -f BAM -p 100 -o BED/NAIVE-1-ATAC.bed
macs3 callpeak -f BED --nomodel --shift -37 --extsize 73 -B --broad -g 2652783500 --keep-dup all -q 0.01 -n NAIVE-1-ATAC -t NAIVE-1-ATAC.bed --outdir ./PEAKS 2>> macs3.log
macs3 randsample -i ALIGNMENTS/NAIVE-2-ATAC-v2.bam -f BAM -p 100 -o BED/NAIVE-2-ATAC.bed
macs3 callpeak -f BED --nomodel --shift -37 --extsize 73 -B --broad -g 2652783500 --keep-dup all -q 0.01 -n NAIVE-2-ATAC -t NAIVE-2-ATAC.bed --outdir ./PEAKS 2>> macs3.log
macs3 randsample -i ALIGNMENTS/NAIVE-3-ATAC-v2.bam -f BAM -p 100 -o BED/NAIVE-3-ATAC.bed
macs3 callpeak -f BED --nomodel --shift -37 --extsize 73 -B --broad -g 2652783500 --keep-dup all -q 0.01 -n NAIVE-3-ATAC -t NAIVE-3-ATAC.bed --outdir ./PEAKS 2>> macs3.log
```
The resulting *_peaks.broadPeak files can also be visualized on IGV. 
