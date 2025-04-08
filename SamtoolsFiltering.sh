#!/bin/bash

# Ulises Gaspar Eslava
# National Autonomous University of Mexico, 2025
# hectog@lcg.unam.mx

# Performs standard sorting, indexing and filtering steps for ATACseq using bedtools:
# 1. Filter non unique alignments
# 2. Sort and Index alinged reads
# 3. Remove mitochondrial reads
# 4. Remove duplicate reads, sort and index filtered reads

# make sure to first modify path to SAM files and to filtered BAM files
# usage: bash SamtoolsFiltering.sh


#Location of align,ents (sam files)
SAM="/home/hectorg/space3/ATACseq/BOWTIE/CL13-v-Naive"

#Location to save sorted, filtered, indexed files
ALIGNMENTS="/home/hectorg/space3/ATACseq/ALIGNMENTS/CL13-v-Naive"

for file in $SAM/*.sam; do
        #Extract file name (without extension)
        FILE_NAME=$(basename "$file" | cut -d '.' -f1)
        FILE_NAME_TRUE=$(basename "$file")
        echo "Processing file: $FILE_NAME"

        ORG_READS=$(samtools view -c "$SAM/$FILE_NAME_TRUE")
        echo "Original read #: $ORG_READS"

        #Run samtools
        #Filter non-unique alignments, sort and index
        samtools view -b -q 10 -o "$ALIGNMENTS/$FILE_NAME".bam "$SAM/$FILE_NAME_TRUE"
        samtools sort -o "$ALIGNMENTS/$FILE_NAME".bam "$ALIGNMENTS/$FILE_NAME".bam
        samtools index "$ALIGNMENTS/$FILE_NAME".bam
        UNIQUE_ALN=$(samtools view -c "$ALIGNMENTS/$FILE_NAME".bam)
        echo "minus low quality alignments: $UNIQUE_ALN"

        #Remove mitochondrial aligned reads by header ID, save to separate, v2, file
        samtools idxstats "$ALIGNMENTS/$FILE_NAME".bam | cut -f 1 | grep -v MT | xargs samtools view -b "$ALIGNMENTS/$FILE_NAME".bam > "$ALIGNMENTS/$FILE_NAME"-v2.bam
        NO_MT_READS=$(samtools view -c "$ALIGNMENTS/$FILE_NAME"-v2.bam)
        echo "Minus mitochondrial reads: $NO_MT_READS"

        #Sort by name of reads, fill coordinates, sort by coordinates, remove duplicates, save to v2 file
        samtools sort -n "$ALIGNMENTS/$FILE_NAME"-v2.bam -u | \
        samtools fixmate - - -u | \
        samtools sort - -u | \
        samtools markdup -r - \
        "$ALIGNMENTS/$FILE_NAME"-v2.bam
        samtools index "$ALIGNMENTS/$FILE_NAME"-v2.bam
        NO_DUP=$(samtools view -c "$ALIGNMENTS/$FILE_NAME"-v2.bam)
        echo "Minus duplicates: $NO_DUP"

done
