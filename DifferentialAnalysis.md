# Differential accessibility analysis
**Ulises Gaspar Eslava**

[Reske et al. (2020)](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-020-00342-y) recommendations for ATAC seq differential analysis testing:
- Compare different normalization methods
- csaw and edgeR for differential analysis
- DiffBind can also be used

Both csaw and DiffBind workflows are present in this file. 
## csaw workflow
This is an adaptation of Reske et al.'s csaw workflow for 3 biological repeats. 
### Read in data:
```R
library(GenomicRanges)
library(csaw)
library(ChIPpeakAnno)
```
Using the previous pipeline's MACS3 broadPeaks:
```R
#Read in genomic coordinates
Cl13.1.peaks <- read.table("CL13-1-ATAC_peaks.broadPeak", sep="\t")[,1:3]
Cl13.2.peaks <- read.table("CL13-2-ATAC_peaks.broadPeak", sep="\t")[,1:3]
Cl13.3.peaks <- read.table("CL13-3-ATAC_peaks.broadPeak", sep="\t")[,1:3]
Naive.1.peaks <- read.table("NAIVE-1-ATAC_peaks.broadPeak", sep="\t")[,1:3]
Naive.2.peaks <- read.table("NAIVE-2-ATAC_peaks.broadPeak", sep="\t")[,1:3]
Naive.3.peaks <- read.table("NAIVE-3-ATAC_peaks.broadPeak", sep="\t")[,1:3]
#Format appropiately
colnames(Cl13.1.peaks) <- c("chrom", "start", "end")
colnames(Cl13.2.peaks) <- c("chrom", "start", "end")
colnames(Cl13.3.peaks) <- c("chrom", "start", "end")
colnames(Naive.1.peaks) <- c("chrom", "start", "end")
colnames(Naive.2.peaks) <- c("chrom", "start", "end")
colnames(Naive.3.peaks) <- c("chrom", "start", "end")
#Convert to GRanges objects
Cl13.1.peaks <- GRanges(Cl13.1.peaks)
Cl13.2.peaks <- GRanges(Cl13.2.peaks)
Cl13.3.peaks <- GRanges(Cl13.3.peaks)
Naive.1.peaks <- GRanges(Naive.1.peaks)
Naive.2.peaks <- GRanges(Naive.2.peaks)
Naive.3.peaks <- GRanges(Naive.3.peaks)
```
### Define consensus peak set:
Filter genomic regions to test only those with consistently high signal across replicates within each condition.
One of the methods to define this peak set is to intersect between all biological replicates; then use the union between both conditions
```R
Cl13.peaks <- Reduce(intersect, list(Cl13.1.peaks, Cl13.2.peaks, Cl13.3.peaks))
Naive.peaks <- Reduce(intersect, list(Naive.1.peaks, Naive.2.peaks, Naive.3.peaks))
all.peaks <- union(Cl13.peaks, Naive.peaks)
```
### Load reads
```R
#Specify bam files
bams <- c("CL13-1-ATAC-v2.bam", "CL13-2-ATAC-v2.bam", "CL13-3-ATAC-v2.bam", "NAIVE-1-ATAC-v2.bam", "NAIVE-2-ATAC-v2.bam", "NAIVE-3-ATAC-v2.bam")
#Read in blacklist
blacklist <- read.table("ENCFF547MET.bed", sep = "\t")
colnames(blacklist) <- c("chrom", "start", "end")
blacklist <- GRanges(blacklist)
#Parameters to load reads: Only standard chromosmes, discard genomic regions in blacklist
Chromosomes <- paste0("chr", c(1:19, "X", "Y"))
param <- readParam(pe = "none", discard = blacklist, restrict = Chromosomes)
merged.peak.counts <- regionCounts(se.bams, all.peaks, param = param)
```
### Filter peaks by abundace
```R
#This will remove peaks with low read abundance that may mess with downstream analysis
library(edgeR)
peak.abundances <- aveLogCPM(asDGEList(merged.peak.counts))
peak.counts.filt <- merged.peak.counts[peak.abundances > -2.5, ]
min(peak.abundances)  #Usefull for modifying abundance cutoff
```
### csaw normalization
2 normalization methods from the csaw library will be tested: Loess and TMM

TMM generates scaling factors from counts in large genomic bins:
```R
binned <- windowCounts(se.bams, bin = TRUE, width = 10000, param = param)
```
```R
#CSAW normalization
peak.counts.tmm <- peak.counts.filt
peak.counts.tmm <- normFactors(binned, se.out=peak.counts.tmm)

#Loess normalization
peak.counts.loess <- peak.counts.filt
peak.counts.loess <- normOffsets(peak.counts.loess, se.out = TRUE)
```
```R
#Comment out as needed:
#working.windows <- peak.counts.tmm
working.windows <- peak.counts.loess
```
### Differential analysis 
Set up data, design matrix:
```R
#Format data
y <- asDGEList(working.windows)
colnames(y$counts) <- c("Cl13_1", "Cl13_2", "Cl13_3", "Naive_1", "Naive_2", "Naive_3")
rownames(y$samples) <- c("Cl13_1", "Cl13_2", "Cl13_3", "Naive_1", "Naive_2", "Naive_3")
y$samples$group <- c("Cl13", "Cl13", "Cl13", "Naive", "Naive", "Naive")
#set up model matrix
design <- model.matrix(~0+group, data = y$samples)
colnames(design) <- c("Cl13", "Naive")
#design
```
Test for differentially accessible regions:
```R
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust = TRUE)
results <- glmQLFTest(fit, contrast=makeContrasts(Naive-Cl13, levels = design))
#head(results$table)
```
Adjust P-values, combine GRanges with DA statistics:
```R
results$table$FDR <- p.adjust(results$table$PValue, method = "fdr")
rowData(working.windows) <- cbind(rowData(working.windows), results$table)
#working.windows@rowRanges
working.windows.final <- working.windows
working.windows.final@rowRanges$sig <- "n.s."
#Set threshold as desired
working.windows.final@rowRanges$sig[working.windows.final@rowRanges$FDR < 0.05] <- "significant"
```
### Annotation:
```R
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
annotations <- detailRanges(working.windows.final@rowRanges,
                            txdb = TxDb.Mmusculus.UCSC.mm10.knownGene,
                            orgdb = org.Mm.eg.db,
                            promoter = c(1500,500), dist = 10*1e3L)
working.windows.final@rowRanges$overlap <- annotations$overlap
working.windows.final@rowRanges$left <- annotations$left
working.windows.final@rowRanges$right <- annotations$right
#working.windows.final@rowRanges
```
### Visualization
Volcano plot:
```R
library(ggrepel)
library(ggrastr)

ATACdata <- as.data.frame(working.windows.final@rowRanges)
ATACdata$annotated <- ATACdata$overlap != ""
ATACdata$label <- NA
ATACdata$label[ATACdata$annotated == TRUE & ATACdata$FDR < 0.05] <- ATACdata$overlap[ATACdata$annotated == TRUE & ATACdata$FDR < 0.05]


ATACdata$sig <- ifelse(ATACdata$sig=="significant",TRUE,FALSE)
#Set desired threshold
ATACdata$Dir <- "n.s"
ATACdata$Dir[ATACdata$logFC > 0.6 & ATACdata$FDR < 0.05] <- "Naive"
ATACdata$Dir[ATACdata$logFC < -0.6 & ATACdata$FDR < 0.05] <- "CL13"
ATACdata$Dir <- factor(ATACdata$Dir, levels = c("Naive", "CL13", "n.s"))

ggplot(data = ATACdata, aes(x = logFC, y = -log10(FDR), col = Dir)) +
  geom_point() +
  theme_minimal() +
  geom_vline(xintercept = c(-0.6, 0.6), col="gray25", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "gray25", linetype = "dashed") +
  scale_color_manual(values = c("red", "blue", "gray")) +
  labs(title = "DA regions - Loess normalization", color = "Accessibility") +
  geom_text_repel(aes(label = ifelse(!is.na(label), label, "")), label.padding = 0.25, point.padding = 1e-06)
```
![image](https://github.com/user-attachments/assets/9a69504c-2d0c-423b-81c1-114ab41ab8e4)
MA plot (useful to visualize global biases, test normalization methods...):
```R
ggplot(data=ATACdata,aes(x = logCPM, y = logFC, col = factor(sig, levels=c("FALSE", "TRUE")))) + 
  geom_point() + scale_color_manual(values = c("black", "red")) + 
  geom_smooth(inherit.aes=F, aes(x = logCPM, y = logFC), method = "loess") + 
  geom_hline(yintercept = 0, color = "darkgreen", linetype = "dashed") + 
  labs(title = "DA regions - Loess normalization", color = "Significance") +
  geom_text_repel(aes(label = ifelse(!is.na(label), label, "")))
```
![image](https://github.com/user-attachments/assets/d0f62a74-c21c-4814-b4a7-c72bddfbd65d)
### Generate custom BED files
For downstream analysis, where subsets of peaks can be usefull to identify potentially active cis regulatory elements, transcription factor binding motifs, etc.
An example: 
```R
NaiveDApeaks <- DA_regions[DA_regions$logFC > 0,]
df_naive <- data.frame(seqnames = seqnames(NaiveDApeaks),
                 starts=start(NaiveDApeaks)-1,
                 ends=end(NaiveDApeaks),
                 names=elementMetadata(NaiveDApeaks)[,10],
                 scores=elementMetadata(NaiveDApeaks)[,5],
                 strands=strand(NaiveDApeaks))
#head(df_naive)
write.table(df_naive, file="NaivePeaks.bed", quote=F, sep="\t", row.names=F, col.names=F)
```
## DiffBind workflow
```R
library(DiffBind)
```
### Create DBA object 
Using MACS3 excel output.
Define experiment, load corresponding metadata (name, condition, replicate number...):
```R
Cl13_vs_Naive <- dba.peakset(NULL, 
                             peaks = "CL13-1-ATAC_peaks.xls",
                             peak.caller = "macs", sampID = "Cl13-1", condition = "Infected", replicate = 1,
                             bamReads = "CL13-1-ATAC-v2.bam")
Cl13_vs_Naive <- dba.peakset(Cl13_vs_Naive,
                             peaks = "CL13-2-ATAC_peaks.xls",
                             peak.caller = "macs", sampID = "CL13-2", condition = "Infected", replicate = 2,
                             bamReads = "CL13-2-ATAC-v2.bam")
Cl13_vs_Naive <- dba.peakset(Cl13_vs_Naive,
                             peaks = "CL13-3-ATAC_peaks.xls",
                             peak.caller = "macs", sampID = "CL13-3", condition = "Infected", replicate = 3,
                             bamReads = "CL13-3-ATAC-v2.bam")

Cl13_vs_Naive <- dba.peakset(Cl13_vs_Naive,
                             peaks = "NAIVE-1-ATAC_peaks.xls",
                             peak.caller = "macs", sampID = "NAIVE-1", condition = "Uninfected", replicate = 1,
                             bamReads = "NAIVE-1-ATAC-v2.bam")
Cl13_vs_Naive <- dba.peakset(Cl13_vs_Naive,
                             peaks = "NAIVE-2-ATAC_peaks.xls",
                             peak.caller = "macs", sampID = "NAIVE-2", condition = "Uninfected", replicate = 2,
                             bamReads = "NAIVE-2-ATAC-v2.bam")
Cl13_vs_Naive <- dba.peakset(Cl13_vs_Naive,
                             peaks = "NAIVE-3-ATAC_peaks.xls",
                             peak.caller = "macs", sampID = "NAIVE-3", condition = "Uninfected", replicate = 3,
                             bamReads = "NAIVE-3-ATAC-v2.bam")
#Filter out regions in blacklist
dba.blacklist(Cl13_vs_Naive, blacklist = TRUE, greylist = FALSE)
#Overlap within replicates can be visualized with a Venn diagram
#dba.plotVenn(Cl13_vs_Naive, Cl13_vs_Naive$masks$Infected, main = "Peak overlaps in CL13 infected replicates")
#dba.plotVenn(Cl13_vs_Naive, Cl13_vs_Naive$masks$Uninfected, main = "Peak overlaps in naive replicates")
```
Load reads to memory (takes a while): 
```R
Cl13_vs_Naive_counts <- dba.count(Cl13_vs_Naive, bParallel = TRUE, score = DBA_SCORE_READS)
```
A correlation heatmap of all peaks can be ploted. Good for QC and to see if there's correlation within experimental groups:
```R
plot(Cl13_vs_Naive_counts, main = "Correlation heatmap of samples")
```
![image](https://github.com/user-attachments/assets/f671cac7-6976-4e0c-a406-7917e19dba91)
A consensus peak set again has to be defined. Filter to genomic intervals present in all 3 repeats per condition:
```R
Cl13_vs_Naive_occupancy <- dba.peakset(Cl13_vs_Naive, consensus = DBA_CONDITION, minOverlap = 3)
#Cl13_vs_Naive_occupancy
#Overlap of both consensus can be visualized on a Venn map
#dba.plotVenn(Cl13_vs_Naive_occupancy, Cl13_vs_Naive_occupancy$masks$Consensus, main = "Overlap of Cl13 and Naive")
```
### Normalization
There are multiple normalization methods in DiffBind. Reske et al. recommend Reads in peaks for ATAC seq:
```R
Cl13_vs_Naive_counts <- dba.normalize(Cl13_vs_Naive_counts, normalize = DBA_NORM_LIB, library = DBA_LIBSIZE_PEAKREADS)
```
### Differential analysis
```R
#Se up model for differential analysis testing
Cl13_vs_Naive.model <- dba.contrast(Cl13_vs_Naive_counts,
                                    reorderMeta=list(Condition="Infected"))
#Cl13_vs_Naive.model


```



