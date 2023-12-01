
# **Transcriptome and gene regulatory analysis of cotton omics project**

**Author: Li'ang Yu, Andrew Nelson**\
**Date: June 14th, 2022**

<img width="800" alt="image" src="https://user-images.githubusercontent.com/69836931/214942202-6e3ce9af-fea6-43e3-ae6c-6d43a4665c4d.png">


- [**Transcriptome and gene regulatory analysis of cotton omics project**](#transcriptome-and-gene-regulatory-analysis-of-cotton-omics-project)
  - [Quantification of reads counts using featureCounts](#quantification-of-reads-counts-using-featurecounts)
    - [Perform the featureCounts for coding genes at gene level](#perform-the-featurecounts-for-coding-genes-at-gene-level)
    - [check if meta-feature and features counted correctly for transcripts](#check-if-meta-feature-and-features-counted-correctly-for-transcripts)
  - [Normalization of count numbers](#normalization-of-count-numbers)
    - [Normalization using DESEq2](#normalization-using-deseq2)
    - [Import raw count data](#import-raw-count-data)
    - [build matrix for each sample based on tissues, accessions, and replicates](#build-matrix-for-each-sample-based-on-tissues-accessions-and-replicates)
    - [normalization and DEGs with experiemtnal design been specified](#normalization-and-degs-with-experiemtnal-design-been-specified)
    - [Count Normalization](#count-normalization)
  - [Weighted Co-expression analysis](#weighted-co-expression-analysis)
    - [Load packages](#load-packages)
    - [Load expression data and meta-table](#load-expression-data-and-meta-table)
    - [Load phenotype data](#load-phenotype-data)
    - [Calculate the average for each two replicates](#calculate-the-average-for-each-two-replicates)
    - [Load different dataset for analysis](#load-different-dataset-for-analysis)
    - [Calculate the soft thresholding power (beta)](#calculate-the-soft-thresholding-power-beta)
    - [STEP Build co-expression matrix using one-step approach](#step-build-co-expression-matrix-using-one-step-approach)
    - [Visulize the cluster dendrogram](#visulize-the-cluster-dendrogram)
    - [plotDendroAndColors](#plotdendroandcolors)
    - [Plot the eigengenes dendrogram](#plot-the-eigengenes-dendrogram)
    - [Plot Epigengene values from different modules across samples](#plot-epigengene-values-from-different-modules-across-samples)
    - [Module-trait relationship](#module-trait-relationship)
      - [Build correlation between metabolites and expression module](#build-correlation-between-metabolites-and-expression-module)
    - [Plot the heatmap to show metabolites with high-correlation](#plot-the-heatmap-to-show-metabolites-with-high-correlation)
      - [Build correlation between phenotypes and expression module](#build-correlation-between-phenotypes-and-expression-module)
      - [Build correlation between NDVI and expression module](#build-correlation-between-ndvi-and-expression-module)
    - [Correlation of triat-module from certain modules](#correlation-of-triat-module-from-certain-modules)
      - [Calculate the gene-module membership](#calculate-the-gene-module-membership)
    - [GS\_MM based on metabolites](#gs_mm-based-on-metabolites)
    - [GS\_MM based on fiber-yield related traits](#gs_mm-based-on-fiber-yield-related-traits)
    - [Replot the GS-MM for highlight purpose](#replot-the-gs-mm-for-highlight-purpose)
    - [Visulize the gene netowrk](#visulize-the-gene-netowrk)
    - [Export list of genes from each module](#export-list-of-genes-from-each-module)
  - [Plot correlation between gene expression and traits](#plot-correlation-between-gene-expression-and-traits)
    - [Load packages and transformed tables](#load-packages-and-transformed-tables)
    - [Reformat the QTls with TPM data](#reformat-the-qtls-with-tpm-data)
    - [Reformat the QTls with Z-score data](#reformat-the-qtls-with-z-score-data)
    - [Plot QTLs correlation with TPM](#plot-qtls-correlation-with-tpm)
    - [Plot with normalized Z-score for QTLs](#plot-with-normalized-z-score-for-qtls)
    - [Reformat the TFs with TPM data for correlation analysis](#reformat-the-tfs-with-tpm-data-for-correlation-analysis)
    - [Reformat the TFs with Zscore data for correlation analysis](#reformat-the-tfs-with-zscore-data-for-correlation-analysis)
    - [Combine with selected traits data](#combine-with-selected-traits-data)
    - [Plot TFs correlation with TPM](#plot-tfs-correlation-with-tpm)
    - [Plot with normalized Z-score for TFs](#plot-with-normalized-z-score-for-tfs)
  - [Correaltion analysis between Log2FC and phenotype ratio (WL/WW) conditions](#correaltion-analysis-between-log2fc-and-phenotype-ratio-wlww-conditions)
    - [Define correlation function](#define-correlation-function)
    - [Load Log2FC data derived](#load-log2fc-data-derived)
    - [define new dataframe for calculation output (10 traits from df1 to df10)](#define-new-dataframe-for-calculation-output-10-traits-from-df1-to-df10)
    - [Calculate the correlation between gene expression and phenotype (production traits)](#calculate-the-correlation-between-gene-expression-and-phenotype-production-traits)
    - [Calculate the correlation between gene expression and phenotype (vegetative indices)](#calculate-the-correlation-between-gene-expression-and-phenotype-vegetative-indices)
    - [Filter and the correlation](#filter-and-the-correlation)
    - [Statistics of correlation](#statistics-of-correlation)
    - [Plot clustering heatmap of metaboltites assocaited with effects](#plot-clustering-heatmap-of-metaboltites-assocaited-with-effects)
      - [Plot the metaboliets with Genotype variations](#plot-the-metaboliets-with-genotype-variations)
      - [Plot the metaboliets with Contion variations](#plot-the-metaboliets-with-contion-variations)
      - [Plot the metaboliets with Interaction variations](#plot-the-metaboliets-with-interaction-variations)
    - [Pair-wise comparison of metaboltes ratio to phenotype ratio](#pair-wise-comparison-of-metaboltes-ratio-to-phenotype-ratio)
    - [Pair-wise comparison of gene expression to metabolites level](#pair-wise-comparison-of-gene-expression-to-metabolites-level)
    - [Plot effects distribution](#plot-effects-distribution)
    - [Plot effects distribution in venn plot](#plot-effects-distribution-in-venn-plot)
  - [Genomic variants analysis of 22 accessions using WGS data](#genomic-variants-analysis-of-22-accessions-using-wgs-data)
    - [Trim and summarize the reads information](#trim-and-summarize-the-reads-information)
    - [Variants calling of 18 accessions using (DNA-seq data)](#variants-calling-of-18-accessions-using-dna-seq-data)
    - [Reformat the VCF with index information added](#reformat-the-vcf-with-index-information-added)
    - [Filter the variants of multiple samples](#filter-the-variants-of-multiple-samples)
    - [Remove loci based on cutoff derived from plot information](#remove-loci-based-on-cutoff-derived-from-plot-information)
    - [Perform the LD-prune to reduce marker density](#perform-the-ld-prune-to-reduce-marker-density)
    - [Annotate variants using SnpEff](#annotate-variants-using-snpeff)
    - [Annotate variants using VEP (Variant Effect Predictor)](#annotate-variants-using-vep-variant-effect-predictor)
    - [Extract the variants annotation of interested genes](#extract-the-variants-annotation-of-interested-genes)
    - [Build the phylogeny using IQtree](#build-the-phylogeny-using-iqtree)
    - [Extract the SNPs/INDELs data for IGV visulization purpose](#extract-the-snpsindels-data-for-igv-visulization-purpose)
  - [Cistrome Deep learning-based prediction](#cistrome-deep-learning-based-prediction)
    - [Extract the positive and negative sequences for training](#extract-the-positive-and-negative-sequences-for-training)
    - [Extract the promoter sequences (reverse and forward sequences)](#extract-the-promoter-sequences-reverse-and-forward-sequences)
    - [Run DL model tranining and CREs prediction](#run-dl-model-tranining-and-cres-prediction)
  - [DAP-seq data analysis](#dap-seq-data-analysis)
    - [Data preparation](#data-preparation)
    - [Mapping reads to the reference genome](#mapping-reads-to-the-reference-genome)
    - [Calling of DAPseq peaks](#calling-of-dapseq-peaks)
    - [Filter the high-confidence peaks](#filter-the-high-confidence-peaks)
    - [Intersect Variants data with peak profiles](#intersect-variants-data-with-peak-profiles)
    - [DAP-seq data reshape and extraction](#dap-seq-data-reshape-and-extraction)
      - [Install packages](#install-packages)
      - [Load libraries](#load-libraries)
      - [Peak annotation step](#peak-annotation-step)
      - [Check the log2FC of genes been bound by DREB2A and HSFA6B](#check-the-log2fc-of-genes-been-bound-by-dreb2a-and-hsfa6b)
      - [Load the list of genes from two modules](#load-the-list-of-genes-from-two-modules)
  - [Population genomics: SNP analysis of IPS flanking SNPs](#population-genomics-snp-analysis-of-ips-flanking-snps)
    - [BYU panel with deep sequencing data](#byu-panel-with-deep-sequencing-data)
    - [TAMU panel with deep sequencing data](#tamu-panel-with-deep-sequencing-data)
    - [Using the published larger panel for analysis](#using-the-published-larger-panel-for-analysis)
  - [Manuscript preparation plot (Misc information)](#manuscript-preparation-plot-misc-information)
    - [Plot the CV over differnece modules](#plot-the-cv-over-differnece-modules)
    - [Plot the line graph of non-lable probe shift comparison](#plot-the-line-graph-of-non-lable-probe-shift-comparison)


## Quantification of reads counts using featureCounts 
The step is to count the reads per primary transcripts using [**featureCounts**](https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html). 


**Example header of Sorghum.gff_primary_transcripts.gff3**\
exon as "feature"\
Parent as "meta-fature"/"attribute"
```
Chr10	phytozomev12	mRNA	6867058	6894701	.	+	.	ID=Sobic.010G080700.2;Name=Sobic.010G080700.2	pacid=37908434	longest=1	ancestorIdentifier=Sobic.010G080700.2.v2.1
Chr10	phytozomev12	five_prime_UTR.	6867058	6867546	.	+	.	ID=Sobic.010G080700.2.five_prime_UTR.1;Parent=Sobic.010G080700.2	pacid=37908434		
Chr10	phytozomev12	five_prime_UTR.	6867920	6868028	.	+	.	ID=Sobic.010G080700.2.five_prime_UTR.2;Parent=Sobic.010G080700.2	pacid=37908434		
Chr10	phytozomev12	exon	6868029	6868060	.	+	0	ID=Sobic.010G080700.2:exon:1;Parent=Sobic.010G080700.2	pacid=37908434		
Chr10	phytozomev12	exon	6868136	6868238	.	+	1	ID=Sobic.010G080700.2:exon:2;Parent=Sobic.010G080700.2	pacid=37908434		
Chr10	phytozomev12	exon	6868579	6868700	.	+	0	ID=Sobic.010G080700.2:exon:3;Parent=Sobic.010G080700.2	pacid=37908434		
Chr10	phytozomev12	exon	6868787	6868879	.	+	1	ID=Sobic.010G080700.2:exon:4;Parent=Sobic.010G080700.2	pacid=37908434		
Chr10	phytozomev12	exon	6870395	6870461	.	+	1	ID=Sobic.010G080700.2:exon:5;Parent=Sobic.010G080700.2	pacid=37908434		
Chr10	phytozomev12	exon	6870544	6870634	.	+	0	ID=Sobic.010G080700.2:exon:6;Parent=Sobic.010G080700.2	pacid=37908434		
Chr10	phytozomev12	exon	6870733	6870785	.	+	2	ID=Sobic.010G080700.2:exon:7;Parent=Sobic.010G080700.2	pacid=37908434		
```

### Perform the featureCounts for coding genes at gene level
```bash
gff_dir="/mnt/Knives/1_data/9_Cotton_omics/0_data"
BAM_dir="/mnt/Knives/1_data/9_Cotton_omics/0_data/1_Bam" 
output_dir="/mnt/Knives/1_data/9_Cotton_omics/2_featureCounts"

###Perform featureCounts for coding genes
featureCounts -T 64 -s 1 -p \
  -a $gff_dir/Ghirsutum_527_v2.1.gene.gff3 \
  -o $output_dir/Cotton_gene-count \
  -t gene \
  -g ID \
  -O \
  $BAM_dir/*.bam > $output_dir/Cotton_gene-count.log
  
grep "Total reads" $output_dir/Cotton_gene-count.log | cut -d " " -f8 > $output_dir/transcript_gene_Total_reads

###Clean GeneID
grep -v 'featureCounts v1.6.0; Command' $output_dir/Cotton_gene-count | cut -d ":" -f2 > $output_dir/Cotton_gene-count.clean
sed -i 's@/mnt/Knives/1_data/5_MODs_bam/sorghum/@@g' $output_dir/Sorghum_transcript_clean.count
sed -i "s@/mnt/Knives/1_data/5_MODs_bam/sorghum/1_Sample/@@g" Sorghum_transcript.count.summary
```
### check if meta-feature and features counted correctly for transcripts
```
Features : 75376                                                        
Meta-features : 75376                                                   
Chromosomes/contigs : 215  
```
## Normalization of count numbers
### Normalization using DESEq2
Using R to to Normalization based on libraries
Using R to identify DEGs across pair-wise samples
Using vsd-tranformation to prepare the input for WGCNA

```r
library(DESeq2)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(knitr)
library(RColorBrewer)
```

### Import raw count data
```r
#Input list of genes from one dataset been modififed 
getwd()
Cotton_meta <- read.table("C:/Users/Leon/OneDrive - Cornell University/Projects/6_Cotton/0_data/Cotton_FeatureCount_gene_20210318.txt", header = TRUE)

#Place Gene_IDs into appropriate format for DESeq2
names(Cotton_meta) <- sub("X", "", names(Cotton_meta))
row.names(Cotton_meta) <- Cotton_meta[,1]
Cotton <- Cotton_meta[,-1]

#Convert NAs to 0s
Cotton[is.na(Cotton)] = 0 
Cotton <- Cotton[,6:97]

#Check the format
head(Cotton)
``` 

### build matrix for each sample based on tissues, accessions, and replicates
```r
#import Matrix (rep,treatment, accessions)
info.df <- read.csv("C:/Users/Leon/OneDrive - Cornell University/Projects/6_Cotton/0_data/0_sample_info/2019_Cotton_samples_RNA_seq_updated_with_actual_FASTq_IDs.csv", header = T)
info.clean <- info.df[,c(5,6,9,12)]

#Extract rownames for RNA-seq sample
sample.ID <- data.frame(FASTq.ID = colnames(Cotton))

#Join the description of samples to RNAseq samples
samples <- inner_join(info.clean,sample.ID, by = "FASTq.ID")
samples <- samples[c("FASTq.ID", "treatment", "rep","entry")]
samples$rep <- as.character(samples$rep)

write.csv(samples,"C:/Users/Leon/OneDrive - Cornell University/Projects/6_Cotton/4_Cluster/sample_map.csv", row.names = F)
```

### normalization and DEGs with experiemtnal design been specified
```r
#start normalization
dds <- DESeqDataSetFromMatrix(countData = Cotton,
                              colData = samples[,c(2,3,4)],
                              design = ~entry + treatment + rep + entry:treatment)
```

### Count Normalization
```r
#Continuing normalization
dds <- DESeq(dds)

#Check experimental deign sheet
colData(dds)

#Export the count information
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="C:/Users/Leon/OneDrive - Cornell University/Projects/6_Cotton/3_FeatureCounts/Cotton-FC-normalized_counts.txt", sep="\t", quote=F, col.names=NA)
```


## Weighted Co-expression analysis
### Load packages
```r
library(cluster)
library(dplyr)
library(pheatmap)
library(stringr)
library(RColorBrewer)
library(viridis)
library(matrixStats)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(cowplot)
library(knitr)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
```

### Load expression data and meta-table
```r
#Load expression
TPM <- read.csv("/Users/Leon/OneDrive - Cornell University/Projects/6_Cotton/3_FeatureCounts/Gh_WGCNA-Oct10-TPM1.csv", header = T)
#rename the column name
rownames(TPM) <- TPM$X
TPM <- TPM[,-1]

#Reorder samples orders by colnames
TPM <- TPM[ , order(names(TPM))]
head(TPM)
```

### Load phenotype data
```r
#load trait table
datTrait <- read.csv("/Users/Leon/OneDrive - Cornell University/Projects/6_Cotton/10_phenotype/Metabolites_2019_BLUEs.csv", header = T)
rownames(datTrait) <- datTrait$genotype
datTrait

datTrait1 <- read.csv("/Users/Leon/OneDrive - Cornell University/Projects/6_Cotton/10_phenotype/FIBER_yield_and_quality_2019_BLUEs_TIPO_removed.csv", header = T)
rownames(datTrait1) <- datTrait1$genotype
datTrait1 

datTrait2 <- read.csv("/Users/Leon/OneDrive - Cornell University/Projects/6_Cotton/10_phenotype/VIs_2019_BLUEs.csv", header = T)
rownames(datTrait2) <- datTrait2$genotype
datTrait2

#Subtract WW and WL data
TPM_WL <- TPM[,grepl("WL",colnames(TPM))]
TPM_WW <- TPM[,grepl("WW",colnames(TPM))]
```


### Calculate the average for each two replicates

```r
#Build clean data frame for average numbers
#FOR WL
TPM_WL_clean <- data.frame(matrix(ncol = 21, nrow = nrow(TPM_WL)))
vector <- c(datTrait[grepl("WL",rownames(datTrait)),1])
rownames(TPM_WL_clean) <- rownames(TPM_WL)
colnames(TPM_WL_clean) <- vector

#FOR WW
TPM_WW_clean <- data.frame(matrix(ncol = 21, nrow = nrow(TPM_WW)))
vector <- c(datTrait[grepl("WW",rownames(datTrait)),1])
rownames(TPM_WW_clean) <- rownames(TPM_WW)
colnames(TPM_WW_clean) <- vector

#Build WL list to loop 
WL_list <- unique(str_replace(colnames(TPM_WL_clean),"_WL",""))
#Loop calculation for each accesion average
for (i in 1:length(WL_list)){
  TPM_WL_clean[,grepl(WL_list[i],colnames(TPM_WL_clean))] <- rowMeans(TPM_WL[,grepl(WL_list[i],colnames(TPM_WL))])
}

WW_list <- unique(str_replace(colnames(TPM_WW_clean),"_WW",""))
#Loop calculation for each accesion average
for (i in 1:length(WW_list)){
  TPM_WW_clean[,grepl(WW_list[i],colnames(TPM_WW_clean))] <- rowMeans(TPM_WW[,grepl(WW_list[i],colnames(TPM_WW))])
}
```

```r
TPM_WW_clean <- TPM_WW_clean[,1:21]
TPM_WL_clean <- TPM_WL_clean[,1:21]

### Combined data for heatmap analyiss
TPM_bind <- cbind(TPM_WL_clean, TPM_WW_clean)
head(TPM_bind)

write.csv(TPM_bind,"/Users/Leon/OneDrive - Cornell University/Projects/6_Cotton/3_FeatureCounts/Gh_WGCNA-TPM-ave_bind.csv", quote = F, row.names = T)
```

### Load different dataset for analysis
Load top (Genes with high SDs)

```r
# build the dataExpr file (TOP 10% MAD score genes)
datExpr <- as.data.frame(t(TPM_WL_clean[order(apply(TPM_WL_clean,1,mad), decreasing = T)[1:4432],]))
write.csv(data.frame(colnames(datExpr)),"/Users/Leon/OneDrive - Cornell University/Projects/6_Cotton/14_MapMan/MAD_score.list.csv",quote = F )
datExpr[,1:10]

datTrait <- datTrait[datTrait$treatment == "WL",]
datTrait1 <- datTrait1[datTrait1$treatment == "WL",]
datTrait2 <- datTrait2[datTrait2$treatment == "WL",]

sampleNames = rownames(datExpr)
traitRows = match(sampleNames, datTrait$genotype) 
traitRows1 = match(sampleNames, datTrait1$genotype) 
traitRows2 = match(sampleNames, datTrait2$genotype) 

```

### Calculate the soft thresholding power (beta)
```r
allowWGCNAThreads()

powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Pring estimated power
sft$powerEstimate
k <- softConnectivity(datE=datExpr,power=sft$powerEstimate)

sizeGrWindow(10, 5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main="Check Scale free topology\n")
```

```r
# Plot the results:
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
```

### STEP Build co-expression matrix using one-step approach
```r
net = blockwiseModules(
                 datExpr,
                 power = 7,
                 maxBlockSize = 10000,
                 TOMType = "unsigned", minModuleSize = 30,
                 reassignThreshold = 0, mergeCutHeight = 0.2,
                 numericLabels = TRUE, pamRespectsDendro = FALSE,
                 saveTOMs = F, 
                 verbose = 3, deepSplit = 2
 )
table(net$colors) 
```

### Visulize the cluster dendrogram
```r
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
table(mergedColors)

#Count numbers of genes under each of module
Module_count <- as.data.frame(table(mergedColors))
Module_count
#write.csv(Module_count, "/Users/Leon/Desktop/Cotton/1_WGCNA_Module_count.csv", row.names = F)

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
```

### plotDendroAndColors
```r
#check the sample names and gene names
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#cluster tress (check if any outliers existed among samples)
#sample_colors <- numbers2colors(as.numeric(factor(datTrait$Group)), 
 #                               colors = c("blue","red", "white"),signed = FALSE)

par(mar = c(1,4,3,1),cex=0.8)
plotDendroAndColors(datExpr_tree, sample_colors,
                    groupLabels = colnames(sample),
                    cex.dendroLabels = 0.8,
                    marAll = c(1, 4, 3, 1),
                    cex.rowText = 0.01,
                    main = "Sample clustering with phylogenetic group assigned")

```

### Plot the eigengenes dendrogram
```r
# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate weight from the clinical traits
weight = as.data.frame(datTrait1$Lint_yield);
names(weight) = "weight"
# Add the weight to existing module eigengenes
MET = orderMEs(MEs)

# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
plotHeatmaps = FALSE)

```

### Plot Epigengene values from different modules across samples

```r
#extract epigengene values 
datME=moduleEigengenes(datExpr,moduleColors)$eigengenes

#select certain module for analysis 
which.module="turquoise"

sizeGrWindow(8,7);
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))

#Plot heat map of certain genes under module 
plotMat(t(scale(datExpr[,moduleColors==which.module]) ),
nrgcols=30,rlabels=T,rcols=which.module,
main=which.module, cex.main=1)
#Plot barchart of certain genes under module 
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
ylab="eigengene expression",xlab="array sample")
```
```r
#select certain module for analysis 
which.module="salmon"

sizeGrWindow(8,7);
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))

#Plot heat map of certain genes under module 
plotMat(t(scale(datExpr[,moduleColors==which.module]) ),
nrgcols=30,rlabels=T,rcols=which.module,
main=which.module, cex.main=1)
#Plot barchart of certain genes under module 
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
ylab="eigengene expression",xlab="array sample")

```

```r
#select certain module for analysis 
which.module="magenta"

sizeGrWindow(8,7);
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))

#Plot heat map of certain genes under module 
plotMat(t(scale(datExpr[,moduleColors==which.module]) ),
nrgcols=30,rlabels=T,rcols=which.module,
main=which.module, cex.main=1)
#Plot barchart of certain genes under module 
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
ylab="eigengene expression",xlab="array sample")

```

### Module-trait relationship 
#### Build correlation between metabolites and expression module
```r
#module heatmap
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

moduleColors <- labels2colors(net$colors)
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

corType = "pearson"
robustY = ifelse(corType=="pearson",T,F)

##Build correlation between metabolites and expression module
if (corType=="pearson") {
  modTraitCor = cor(MEs, datTrait[,colnames(datTrait)[3:453]], use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, colnames(datTrait[3:453]) , robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}

textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
modTraitCor_table <- as.data.frame(modTraitCor)

#Export correlaion for each phenotype
write.csv(as.data.frame(modTraitCor), "/Users/Leon/OneDrive - Cornell University/Projects/6_Cotton/2_ModTraitCor_metabolites-treeCut0.2-Oct09.csv")
head(modTraitCor_table, 20)  
```
### Plot the heatmap to show metabolites with high-correlation
```r
### chose certain metabolites for analysis 
metabolites <- c(
"X2019_LCMS_154",
"X2019_LCMS_338",
"X2019_LCMS_411",
"X2019_LCMS_391",
"X2019_LCMS_223",
"X2019_LCMS_277",
"X2019_GCMS_32",
"X2019_LCMS_406",
"X2019_GCMS_20",
"X2019_LCMS_453",
"X2019_GCMS_19",
"X2019_GCMS_26",
"X2019_GCMS_32")

##Build correlation between metabolites and expression module
if (corType=="pearson") {
  modTraitCor = cor(MEs,datTrait[,metabolites], use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, metabolites , robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}

textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)


labeledHeatmap(Matrix = modTraitCor, xLabels = metabolites, 
               yLabels = colnames(MEs), 
               cex.lab = 0.5, 
               ySymbols = colnames(MEs), colorLabels = TRUE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = TRUE, 
               cex.text = 0.25, zlim = c(-1,1),
               main = paste("Module-metabolites relationships"))
```


#### Build correlation between phenotypes and expression module
```r
### Under the WL condition 
if (corType=="pearson") {
  modTraitCor1 = cor(MEs, datTrait1[,colnames(datTrait1)[3:12]], use = "p")
  modTraitP1 = corPvalueStudent(modTraitCor1, nSamples)
} else {
  modTraitCorP1 = bicorAndPvalue(MEs_col, datTrait1, robustY=robustY)
  modTraitCor1 = modTraitCorP1$bicor
  modTraitP1   = modTraitCorP1$p
}

textMatrix1 = paste(signif(modTraitCor1, 2), "\n(", signif(modTraitP1, 1), ")", sep = "")
dim(textMatrix1) = dim(modTraitCor1)

####XLabels = colnames(datTrait)[3:6]

labeledHeatmap(Matrix = modTraitCor1, xLabels = colnames(datTrait1)[3:12], 
               yLabels = colnames(MEs), 
               cex.lab = 0.5, 
               ySymbols = colnames(MEs), colorLabels = TRUE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix1, setStdMargins = TRUE, 
               cex.text = 0.25, zlim = c(-1,1),
               main = paste("Module-fiber_traits relationships"))

```
#### Build correlation between NDVI and expression module
```r

if (corType=="pearson") {
  modTraitCor2 = cor(MEs, datTrait2[,colnames(datTrait2)[3:6]], use = "p")
  modTraitP2 = corPvalueStudent(modTraitCor2, nSamples)
} else {
  modTraitCorP2 = bicorAndPvalue(MEs_col, datTrait2 , robustY=robustY)
  modTraitCor2= modTraitCorP2$bicor
  modTraitP2   = modTraitCorP2$p
}

textMatrix2 = paste(signif(modTraitCor2, 2), "\n(", signif(modTraitP2, 1), ")", sep = "")
dim(textMatrix2) = dim(modTraitCor2)

####XLabels = colnames(datTrait)[3:6]

labeledHeatmap(Matrix = modTraitCor2, xLabels = colnames(datTrait2)[3:6], 
               yLabels = colnames(MEs), 
               cex.lab = 0.5, 
               ySymbols = colnames(MEs), colorLabels = TRUE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix2, setStdMargins = TRUE, 
               cex.text = 0.25, zlim = c(-1,1),
               main = paste("Module-NVDI relationships"))
```


### Correlation of triat-module from certain modules
#### Calculate the gene-module membership
```r
# names (colors) of the modules
modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))

MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

head(geneModuleMembership, 20)
#Export the gene-module membership
#write.csv(geneModuleMembership, "/Users/Leon/Desktop/Cotton/6_Gene_module_membership_Oct03-2021.csv")
```

### GS_MM based on metabolites
```r
# Select interested traits from phenotype 
metab <- datTrait[,colnames(datTrait) == "X2019_GCMS_32"]
title_name <- "X2019_GCMS_32"

geneTraitSignificance <- as.data.frame(cor(datExpr, metab, use = "p"))
head(geneTraitSignificance)

GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(metab), sep="");
names(GSPvalue) = paste("p.GS.", names(metab), sep="");

# Select modules 
module <- "salmon"
column <- match(module, modNames)
moduleGenes <- moduleColors==module

#Delimit the size of windows
sizeGrWindow(7, 7)
par(mfrow = c(1,1))

verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = paste0("Gene significance for ",title_name),
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
```

### GS_MM based on fiber-yield related traits
```r
# Select interested traits from phenotype 
phenotype <- datTrait1[,colnames(datTrait1) == "Lint_yield"]
title_name1 <- "Lint_yield"

geneTraitSignificance1 <- as.data.frame(cor(datExpr, phenotype, use = "p"))
head(geneTraitSignificance1)

GSPvalue1 <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance1), nSamples))
names(geneTraitSignificance1) = paste("GS.", names(phenotype), sep="");
names(GSPvalue1) = paste("p.GS.", names(phenotype), sep="");

# Select modules 
module <- "turquoise"
column <- match(module, modNames)
moduleGenes <- moduleColors==module

#Delimit the size of windows
sizeGrWindow(7, 7)
par(mfrow = c(1,1))

verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance1[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = paste0("Gene significance for ",title_name1),
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)



### Export the file for reformat
ExportGS_MM <- data.frame(GeneID = rownames(geneModuleMembership[moduleGenes,]),
                          MM = abs(geneModuleMembership[moduleGenes, column]), 
                          GS = abs(geneTraitSignificance1[moduleGenes, 1]), 
                          Group = "Fill")

write.csv(ExportGS_MM, 
          "/Users/leon/Documents/Project/Cotton/WGCNA/GS-MM_lintModule-replot.csv", quote = F)
```

### Replot the GS-MM for highlight purpose
```r
### Reload the modified table for plot
Replot <- read.csv("/Users/leon/Documents/Project/Cotton/WGCNA/GS-MM_lintModule-PlotNew.csv", header = T)
head(Replot)


### Plot the data
ggplot(Replot, aes(x=MM, y=GS, color=Group, size = Group, alpha = Group, label = GeneID)) +
      geom_point() +
      scale_alpha_manual(values = c(Fill = 0.8, QTLs=1, Traitcorrelated = 0.8, 
                                    TFs= 0.8, HeatResponse= 0.8, StressResponse = 0.8, CoreGene = 0.8)) +
      scale_size_manual(values = c(Fill = 1, QTLs=7, Traitcorrelated = 4, TFs= 4, 
                                   HeatResponse= 4, StressResponse = 4, CoreGene = 3)) +
      scale_colour_manual(values=c(Fill = "grey70", Traitcorrelated = "darkorange", 
                                   TFs = "red", QTLs = "darkgreen", StressResponse = "dodgerblue", 
                                   HeatResponse = "darkorchid1", CoreGene = "bisque")) + 
      geom_text(data= Replot[Replot$GeneID == "IPS" | Replot$GeneID == "ABP" | Replot$GeneID == "L7T" |
                             Replot$GeneID == "API5"| Replot$GeneID == "DBCP"| Replot$GeneID == "DREB2A" | Replot$GeneID == "HSFA6B",], 
                size = 4, color = "black",
                hjust = 1) +
      theme_bw() + 
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+  
      labs(x=paste0("Module Membership"), 
           y=paste0("Gene significance")) +
      theme(axis.title.x =element_text(size=14), 
            axis.title.y=element_text(size=14),
            axis.text = element_text(size = 12),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            plot.title = element_text(hjust = 0.5,size = 16,face = "bold"),
            plot.margin = unit(rep(2,4),'lines')) +
       theme(legend.position = 'none') +
      geom_hline(aes(yintercept=0.7),colour="black",lwd=0.5,linetype=5)+
      geom_vline(aes(xintercept=0.8),colour="black",lwd=0.5,linetype=5)

```

### Visulize the gene netowrk
```r
TOM = TOMsimilarityFromExpr(datExpr, power=sft$powerEstimate, corType=corType, networkType="unsigned")
TOM <- as.matrix(TOM)
dissTOM = 1-TOM
# Transform dissTOM with a power to make moderately strong 
# connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA

TOMPlot_results <- TOMplot(plotTOM, net$dendrograms, moduleColors, 
        main = "Network heatmap plot, all genes")

nSelect = 4000
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")

```

```r
moduleColors_list <- unique(as.data.frame(moduleColors))
moduleColors_list
```

### Export list of genes from each module 
```r
Module_count

for (i in 1:nrow(Module_count)){
  # Select module
  module <- Module_count[i,1]
  # Select module probes
  probes = colnames(datExpr)
  moduleColors = mergedColors
  inModule = (moduleColors==module);
  modProbes = probes[inModule];
  # Select the corresponding Topological Overlap
  modTOM <- TOM[inModule, inModule]
  dimnames(modTOM) <- list(modProbes, modProbes)
  
  cyt <- exportNetworkToCytoscape(
  modTOM,
  edgeFile = 
  paste0("/Users/leon/Documents/edges-MAD10", paste(module), ".csv"),
  nodeFile =
  paste0("/Users/leon/Documents/nodes-MAD", paste(module), ".csv"),
  weighted = TRUE,
  threshold = 0,
  nodeNames = modProbes, 
  nodeAttr = moduleColors[inModule])
}
```

## Plot correlation between gene expression and traits
### Load packages and transformed tables 
```r
library(dplyr)
library(ggplot2)
library(knitr)
library(RColorBrewer)
library(ggpubr)
library(tidyr)
library(cowplot)
library(knitr)
```

**Trait Tables**
```r
Trait <- read.csv("/Users/Leon/OneDrive - Cornell University/Projects/6_Cotton/10_phenotype/FIBER_yield_and_quality_2019_BLUEs_TIPO_removed.csv", header = T)
rownames(Trait) <- Trait$genotype
head(Trait)

Metabolites_BLUE <- read.csv("/Users/Leon/OneDrive - Cornell University/Projects/6_Cotton/10_phenotype/Metabolites_2019_BLUEs.csv", header = T, row.names = 1)
Metabolites_BLUE$genotype <- rownames(Metabolites_BLUE) 
head(Metabolites_BLUE)
```

**Gene expression Table**
```r
TPM_bind <- read.csv("/Users/Leon/OneDrive - Cornell University/Projects/6_Cotton/3_FeatureCounts/Gh_WGCNA-TPM-ave_bind.csv", header = T, row.names = 1)
TPM_corr = as.data.frame(t(TPM_bind))
head(TPM_corr)

Zscore_bind <- read.csv("/Users/Leon/OneDrive - Cornell University/Projects/6_Cotton/14_MapMan/turquoise_Zscore.csv", header = T, row.names = 1)
Zscore_corr = as.data.frame(t(Zscore_bind))
head(Zscore_corr)
```

**Gene information List**

```r
GenePlot <- read.csv("/Users/Leon/OneDrive - Cornell University/Projects/6_Cotton/17_Cistrome/03_DL-prediction/Plot_GeneList.csv", header = T)
head(GenePlot)
```

### Reformat the QTls with TPM data
```r
### import lint related genes
yield_gene <- GenePlot[GenePlot$Class == "QTLs", ]
yield_gene

for (i in 1:nrow(yield_gene)){
  
  yield_df <- TPM_corr[,yield_gene[i,1], drop = F]
  yield_df$Gene <- yield_gene[i,7]
  yield_df$Accessions <- rownames(yield_df)
  rownames(yield_df) <- NULL
  colnames(yield_df) <- c("TPM", "Gene", "genotype")
  
  assign(paste0("yield_QTLs_", yield_gene[i,1]), yield_df)
}

### Merge QTLs into same datafroma
QTLs_combine <- bind_rows(lapply(ls(pattern = "yield_QTLs_"), get))
QTLs_combine$TPM <- as.numeric(QTLs_combine$TPM)
head(QTLs_combine)

#rm(list=ls(pattern="yield_QTLs_"))
```

### Reformat the QTls with Z-score data

```r
### import lint related genes
yield_gene <- GenePlot[GenePlot$Class == "QTLs", ]
yield_gene

for (i in 1:nrow(yield_gene)){
  
  yield_df <- Zscore_corr[,yield_gene[i,1], drop = F]
  yield_df$Gene <- yield_gene[i,7]
  yield_df$Accessions <- rownames(yield_df)
  rownames(yield_df) <- NULL
  colnames(yield_df) <- c("Zscore", "Gene", "genotype")
  
  assign(paste0("Zsore_yieldQTLs-", yield_gene[i,1]), yield_df)
}

### Merge QTLs into same datafroma
QTLs_Zscore <- bind_rows(lapply(ls(pattern = "Zsore_yieldQTLs-"), get))
QTLs_Zscore$Zscore <- as.numeric(QTLs_Zscore$Zscore)
head(QTLs_Zscore)

#rm(list=ls(pattern="yield_QTLs_"))
```

```r
### Extract values of selected traits
Select_traits <- Trait[,c("genotype", "Lint_yield", "treatment")]
rownames(Select_traits) <- NULL

### combine two sets
Select_Comp <- left_join(QTLs_combine, Select_traits, by = "genotype")
head(Select_Comp)

Select_Zscore <- left_join(QTLs_Zscore, Select_traits, by = "genotype")
head(Select_Zscore)
```

### Plot QTLs correlation with TPM 
```r
Plot_WL <- Select_Comp[Select_Comp$treatment == "WL",]

### Plot each QTLs
for (i in 1:nrow(yield_gene)){
  print(ggscatter(Plot_WL[Plot_WL$Gene == yield_gene[i,7],], 
          x = "TPM", y = "Lint_yield", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = paste0("TPM of ", yield_gene[i,7]), ylab = "Lint yield",
          shape = 21, size = 3, 
          add.params = list(color = "blue", fill = "lightgray")))
}
```

### Plot with normalized Z-score for QTLs

```r
Zscore_WL <- Select_Zscore [Select_Zscore$treatment == "WL",]
head(Zscore_WL)

### Plot each QTLs
Zscore_WL %>%
  ggplot(aes(x=Zscore, y=Lint_yield, color=Gene))+
  geom_point(size = 2, alpha = 0.5)+ 
  scale_color_manual(values = c(IPS="#1B9E77", ABP="tomato", L7T="#E7298A", DBC="steelblue", API5="#E6AB02"))+
  theme_bw()+ 
  ggtitle("Correaltion between lint yield and 5 QTLs Zscore") +
  geom_smooth(method = "lm", aes(color= Gene), size = 0.5, level = 0.95) +
  xlab("Zscore of lint yield QTLs") +
  ylab("Lint yieled (Kg)") +
  stat_cor(label.y.npc = c(1,0.5)) +
  theme(legend.position = c(0.8, 0.2))
```

### Reformat the TFs with TPM data for correlation analysis
```r
### import lint related genes
TF_gene <- GenePlot[GenePlot$Class == "TFs", ]
TF_gene 

for (i in 1:nrow(TF_gene)){
  
  TF_df <- TPM_corr[,TF_gene[i,1], drop = F]
  TF_df$Gene <- TF_gene[i,1]
  TF_df$Accessions <- rownames(TF_df)
  rownames(TF_df) <- NULL
  colnames(TF_df) <- c("TPM", "Gene", "genotype")
  
  assign(paste0("yield_TFs_", TF_gene[i,1]), TF_df)
}

### Merge QTLs into same datafroma
TFs_combine <- bind_rows(lapply(ls(pattern = "yield_TFs_"), get))
TFs_combine$TPM <- as.numeric(TFs_combine$TPM)
head(TFs_combine)

```

### Reformat the TFs with Zscore data for correlation analysis
```r
### import lint related genes

for (i in 1:nrow(TF_gene)){
  
  TF_df <- Zscore_corr[,TF_gene[i,1], drop = F]
  TF_df$Gene <- TF_gene[i,1]
  TF_df$Accessions <- rownames(TF_df)
  rownames(TF_df) <- NULL
  colnames(TF_df) <- c("Zscore", "Gene", "genotype")
  
  assign(paste0("Zscore_yield_TFs-", TF_gene[i,1]), TF_df)
}

### Merge QTLs into same datafroma
TFs_Zscore <- bind_rows(lapply(ls(pattern = "Zscore_yield_TFs-"), get))
TFs_Zscore$Zscore <- as.numeric(TFs_Zscore$Zscore)
head(TFs_Zscore)

```

### Combine with selected traits data
```r
### combine two sets
Select_TF <- left_join(TFs_combine, Select_traits, by = "genotype")
head(Select_TF)

Select_Zscore_TF <- left_join(TFs_Zscore, Select_traits, by = "genotype")
head(Select_Zscore_TF)
```

### Plot TFs correlation with TPM 
```r
Plot_WL_TF <- Select_TF[Select_TF$treatment == "WL",] %>% na.omit()
Plot_WL_TF

### Plot each QTLs
for (i in 1:nrow(TF_gene)){
  print(ggscatter(Plot_WL_TF[Plot_WL_TF$Gene == TF_gene[i,1],], 
          x = "TPM", y = "Lint_yield", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = paste0("TPM of ", TF_gene[i,1]), 
          ylab = "Lint yield (Kg)",
          shape = 21, size = 3, 
          add.params = list(color = "blue", fill = "lightgray"),
          title = paste0("Correlation between ", TF_gene[i,7], " and lint yield")))
}
```

### Plot with normalized Z-score for TFs

```r
Zscore_WL_TF <- Select_Zscore_TF[Select_Zscore_TF$treatment == "WL",]
head(Zscore_WL_TF)
```

**MBFC1**
```r
### Plot each QTLs
MBFC_Zscore <- Zscore_WL_TF[Zscore_WL_TF$Gene == "Gohir.A06G061400" | Zscore_WL_TF$Gene == "Gohir.D06G059700", ] 

MBFC_Zscore %>%
  ggplot(aes(x=Zscore, y=Lint_yield, color=Gene))+
  geom_point(size = 2, alpha = 0.7)+ 
  scale_color_manual(values = c(Gohir.A06G061400="steelblue", Gohir.D06G059700="tomato"))+
  theme_bw()+ 
  ggtitle("Correaltion between lint yield and GhMBFC1 Zscore") +
  geom_smooth(method = "lm", aes(color= Gene), size = 0.5, level = 0.95) +
  stat_cor(label.y.npc = "top") +
  xlab("Zscore of TFs") +
  ylab("Lint yieled (Kg)") + 
  theme(legend.position = c(0.8, 0.2))
```

**DREB2C**
```r
### Plot each QTLs
DREB2C_Zscore <- Zscore_WL_TF[Zscore_WL_TF$Gene == "Gohir.A13G021700" | Zscore_WL_TF$Gene == "Gohir.D13G022300", ] 

DREB2C_Zscore %>%
  ggplot(aes(x=Zscore, y=Lint_yield, color=Gene))+
  geom_point(size = 2, alpha = 0.7)+ 
  scale_color_manual(values = c(Gohir.A13G021700="steelblue", Gohir.D13G022300="tomato"))+
  theme_bw()+ 
  ggtitle("Correaltion between lint yield and GhDREB2C Zscore") +
  geom_smooth(method = "lm", aes(color= Gene), size = 0.5, level = 0.95) +
  stat_cor(label.y.npc = "top") +
  xlab("Zscore of TFs") +
  ylab("Lint yieled (Kg)") +
  theme(legend.position = c(0.8, 0.2))
```

**HSFA6B**
```r
### Plot each QTLs
HSFA6B_Zscore <- Zscore_WL_TF[Zscore_WL_TF$Gene == "Gohir.A08G064100" | Zscore_WL_TF$Gene == "Gohir.D08G072600",] 

HSFA6B_Zscore %>%
  ggplot(aes(x=Zscore, y=Lint_yield, color=Gene))+
  geom_point(size = 2, alpha = 0.7)+ 
  scale_color_manual(values = c(Gohir.A08G064100="steelblue", Gohir.D08G072600="tomato"))+
  theme_bw()+ 
  ggtitle("Correaltion between lint yield and GhHSFA6B expression") +
  geom_smooth(method = "lm", aes(color= Gene), size = 0.5, level = 0.95) +
  stat_cor(label.y.npc = "top") +
  xlab("Zscore of GhHSFA6B") +
  ylab("Lint yieled (Kg)") +
  theme(legend.position = c(0.8, 0.2))
```

**HSF7**
```r
### Plot each QTLs
HSF7_Zscore <- Zscore_WL_TF[Zscore_WL_TF$Gene == "Gohir.D03G153000" | Zscore_WL_TF$Gene == "Gohir.D08G198000" |
                              Zscore_WL_TF$Gene == "Gohir.A08G179400",] 

HSF7_Zscore %>%
  ggplot(aes(x=Zscore, y=Lint_yield, color=Gene))+
  geom_point(size = 2, alpha = 0.7)+ 
  scale_color_manual(values = c(Gohir.A08G179400="steelblue", Gohir.D08G198000="tomato", Gohir.D03G153000 = "black"))+
  theme_bw()+ 
  ggtitle("Correaltion between lint yield and GhHSF7 Zscore") +
  geom_smooth(method = "lm", aes(color= Gene), size = 0.5, level = 0.95) +
  stat_cor(label.y.npc = "top") +
  xlab("Zscore of GhHSF7") +
  ylab("Lint yieled (Kg)") +
  theme(legend.position = c(0.8, 0.2))
```

## Correaltion analysis between Log2FC and phenotype ratio (WL/WW) conditions
```r
library(Hmisc)
library(dplyr)
library(reshape2)
library(ggpubr)
library(ggVennDiagram)
```

```r
load("/Users/Leon/OneDrive - Cornell University/Projects/6_Cotton/2_Pipeline/LogFC.RData")
```

### Define correlation function
```r
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}
```

### Load Log2FC data derived
```r
LOG <- read.csv("/Users/Leon/OneDrive - Cornell University/3_Manuscript/Manuscript8/7_LogFC_ratio/Ratio_LogFC.csv", header = T)
row.names(LOG) = LOG[,1]
head(LOG)

Ratio <- read.csv("/Users/Leon/OneDrive - Cornell University/3_Manuscript/Manuscript8/7_LogFC_ratio/LogFC_meta_table.csv", header = T)
row.names(Ratio) = Ratio$Gene
head(Ratio)

```

### define new dataframe for calculation output (10 traits from df1 to df10)
```r
df1 <- data.frame(matrix(ncol = 7 , nrow = 67413))
df2 <- data.frame(matrix(ncol = 7 , nrow = 67413))
df3 <- data.frame(matrix(ncol = 7 , nrow = 67413))
df4 <- data.frame(matrix(ncol = 7 , nrow = 67413))
df5 <- data.frame(matrix(ncol = 7 , nrow = 67413))
df6 <- data.frame(matrix(ncol = 7 , nrow = 67413))
df7 <- data.frame(matrix(ncol = 7 , nrow = 67413))
df8 <- data.frame(matrix(ncol = 7 , nrow = 67413))
df9 <- data.frame(matrix(ncol = 7 , nrow = 67413))
df10 <- data.frame(matrix(ncol = 7 , nrow = 67413))


colnames(df1) <- c("row", "column", "R", "P-value", "Intercept", "Effects", "R.square")
colnames(df2) <- c("row", "column", "R", "P-value", "Intercept", "Effects", "R.square")
colnames(df3) <- c("row", "column", "R", "P-value", "Intercept", "Effects", "R.square")
colnames(df4) <- c("row", "column", "R", "P-value", "Intercept", "Effects", "R.square")
colnames(df5) <- c("row", "column", "R", "P-value", "Intercept", "Effects", "R.square")
colnames(df6) <- c("row", "column", "R", "P-value", "Intercept", "Effects", "R.square")
colnames(df7) <- c("row", "column", "R", "P-value", "Intercept", "Effects", "R.square")
colnames(df8) <- c("row", "column", "R", "P-value", "Intercept", "Effects", "R.square")
colnames(df9) <- c("row", "column", "R", "P-value", "Intercept", "Effects", "R.square")
colnames(df10) <- c("row", "column", "R", "P-value", "Intercept", "Effects", "R.square")
```

### Calculate the correlation between gene expression and phenotype (production traits)
```r
for (i in 1:(nrow(Ratio))){
  res <- rcorr(t(rbind(LOG[1,-1], Ratio[i,-1])))
  df1[i,1:4] <- flattenCorrMatrix(res$r, res$P)
  temp.df = data.frame(t(rbind(LOG[1,-1], Ratio[i,-1])))
  sm = summary(lm(Lint_yield ~ get(Ratio[i,1]), temp.df))
  df1[i,5] = sm$coefficients[1]
  df1[i,6] = sm$coefficients[2]
  df1[i,7] = sm$r.squared
}

for (i in 1:(nrow(Ratio))){
  res <- rcorr(t(rbind(LOG[2,-1], Ratio[i,-1])))
  df2[i,1:4] <- flattenCorrMatrix(res$r, res$P)
  temp.df = data.frame(t(rbind(LOG[2,-1], Ratio[i,-1])))
  sm = summary(lm(Mic ~ get(Ratio[i,1]), temp.df))
  df2[i,5] = sm$coefficients[1]
  df2[i,6] = sm$coefficients[2]
  df2[i,7] = sm$r.squared
}

for (i in 1:(nrow(Ratio))){
  res <- rcorr(t(rbind(LOG[3,-1], Ratio[i,-1])))
  df3[i,1:4] <- flattenCorrMatrix(res$r, res$P)
  temp.df = data.frame(t(rbind(LOG[3,-1], Ratio[i,-1])))
  sm = summary(lm(UHM ~ get(Ratio[i,1]), temp.df))
  df3[i,5] = sm$coefficients[1]
  df3[i,6] = sm$coefficients[2]
  df3[i,7] = sm$r.squared
}

for (i in 1:(nrow(Ratio))){
  res <- rcorr(t(rbind(LOG[4,-1], Ratio[i,-1])))
  df4[i,1:4] <- flattenCorrMatrix(res$r, res$P)
  temp.df = data.frame(t(rbind(LOG[4,-1], Ratio[i,-1])))
  sm = summary(lm(UI ~ get(Ratio[i,1]), temp.df))
  df4[i,5] = sm$coefficients[1]
  df4[i,6] = sm$coefficients[2]
  df4[i,7] = sm$r.squared
}

for (i in 1:(nrow(Ratio))){
  res <- rcorr(t(rbind(LOG[5,-1], Ratio[i,-1])))
  df5[i,1:4] <- flattenCorrMatrix(res$r, res$P)
  temp.df = data.frame(t(rbind(LOG[5,-1], Ratio[i,-1])))
  sm = summary(lm(Str ~ get(Ratio[i,1]), temp.df))
  df5[i,5] = sm$coefficients[1]
  df5[i,6] = sm$coefficients[2]
  df5[i,7] = sm$r.squared
}

for (i in 1:(nrow(Ratio))){
  res <- rcorr(t(rbind(LOG[6,-1], Ratio[i,-1])))
  df6[i,1:4] <- flattenCorrMatrix(res$r, res$P)
  temp.df = data.frame(t(rbind(LOG[6,-1], Ratio[i,-1])))
  sm = summary(lm(Elo ~ get(Ratio[i,1]), temp.df))
  df6[i,5] = sm$coefficients[1]
  df6[i,6] = sm$coefficients[2]
  df6[i,7] = sm$r.squared
}
```

### Calculate the correlation between gene expression and phenotype (vegetative indices)
```r
### for vegetative indices
for (i in 1:(nrow(Ratio))){
  res <- rcorr(t(rbind(LOG[7,-1], Ratio[i,-1])))
  df7[i,1:4] <- flattenCorrMatrix(res$r, res$P)
  temp.df = data.frame(t(rbind(LOG[7,-1], Ratio[i,-1])))
  sm = summary(lm(NDVI ~ get(Ratio[i,1]), temp.df))
  df7[i,5] = sm$coefficients[1]
  df7[i,6] = sm$coefficients[2]
  df7[i,7] = sm$r.squared
}

for (i in 1:(nrow(Ratio))){
  res <- rcorr(t(rbind(LOG[8,-1], Ratio[i,-1])))
  df8[i,1:4] <- flattenCorrMatrix(res$r, res$P)
  temp.df = data.frame(t(rbind(LOG[8,-1], Ratio[i,-1])))
  sm = summary(lm(sPRI ~ get(Ratio[i,1]), temp.df))
  df8[i,5] = sm$coefficients[1]
  df8[i,6] = sm$coefficients[2]
  df8[i,7] = sm$r.squared
}


for (i in 1:(nrow(Ratio))){
  res <- rcorr(t(rbind(LOG[9,-1], Ratio[i,-1])))
  df9[i,1:4] <- flattenCorrMatrix(res$r, res$P)
  temp.df = data.frame(t(rbind(LOG[9,-1], Ratio[i,-1])))
  sm = summary(lm(CRI ~ get(Ratio[i,1]), temp.df))
  df9[i,5] = sm$coefficients[1]
  df9[i,6] = sm$coefficients[2]
  df9[i,7] = sm$r.squared
}

for (i in 1:(nrow(Ratio))){
  res <- rcorr(t(rbind(LOG[10,-1], Ratio[i,-1])))
  df10[i,1:4] <- flattenCorrMatrix(res$r, res$P)
  temp.df = data.frame(t(rbind(LOG[10,-1], Ratio[i,-1])))
  sm = summary(lm(WI.NDVI ~ get(Ratio[i,1]), temp.df))
  df10[i,5] = sm$coefficients[1]
  df10[i,6] = sm$coefficients[2]
  df10[i,7] = sm$r.squared
}
```

### Filter and the correlation
```r
combined <- bind_rows(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10)
head(combine)

combine_filter <- combined[combined$`P-value` < 0.05,]
combine_filter$Type <- "fill"
combine_filter[combine_filter$R < 0,8] = "N correlated"
combine_filter[combine_filter$R > 0,8] = "Pcorrelated"

combine_filter$ABS = abs(combine_filter$Effects)
head(combine_filter)

write.csv(combine_filter,"/Users/Leon/OneDrive - Cornell University/3_Manuscript/Manuscript8/7_LogFC_ratio/Combined_file.csv", row.names = F, quote = F)
write.csv(combined,"/Users/Leon/OneDrive - Cornell University/3_Manuscript/Manuscript8/7_LogFC_ratio/Combined-non-filter_file.csv", row.names = F, quote = F)
```
### Statistics of correlation
```r
combine_stats <- combine_filter %>% group_by(row, Type) %>% mutate(counts = n())
combine_stats <- unique(combine_stats[,c(1,8,10)])
combine_stats <- dcast(combine_stats, row~Type,sum)
combine_stats

write.csv(combine_stats,"/Users/Leon/OneDrive - Cornell University/3_Manuscript/Manuscript8/7_LogFC_ratio/Combined_stats.csv", row.names = F, quote = F)
```

### Plot clustering heatmap of metaboltites assocaited with effects
```r
### Load metabolties Z-score
Zscore <- read.csv("/Users/Leon/OneDrive - Cornell University/3_Manuscript/Manuscript8/11_Metabolties/Zscores_Met_2019_BLUEs_WW_and_WL_plot.csv", header = T)
sample_col <- read.csv("/Users/Leon/OneDrive - Cornell University/3_Manuscript/Manuscript8/11_Metabolties/sample_col.csv", header = T)
sample_row <- read.csv("/Users/Leon/OneDrive - Cornell University/3_Manuscript/Manuscript8/11_Metabolties/sample_row.csv", header = T)

rownames(sample_row) <- sample_row[,1]
rownames(sample_col) <- sample_col[,1]
rownames(Zscore) <- Zscore[,1]

head(Zscore)
head(sample_col)
head(sample_row)
```

#### Plot the metaboliets with Genotype variations
```r
annotation_colors = list(
  treatment = c(WL="red", WW="royalblue"),
  Group = c(Wild="pink", Commercial="slateblue1", Elite="palegreen"))

Metablites_select <- Zscore[Zscore$Class == "Genotype",3:44]
head(Metablites_select)

pheatmap(Metablites_select,
         annotation_col = sample_col[,2:3], 
         annotation_row = sample_row[sample_row$Type == "Genotype",2:3], 
         color = inferno(30),
         annotation_colors = annotation_colors, 
         clustering_distance_r = "manhattan", 
         show_colnames =F,show_rownames = F,
         cluster_rows = T, cluster_cols = T,
         fontsize_row = 4, fontsize_col = 4, 
         angle_col = 90, border_color = NA)
```

#### Plot the metaboliets with Contion variations
```r
annotation_colors = list(
  treatment = c(WL="red", WW="royalblue"),
  Group = c(Wild="pink", Commercial="slateblue1", Elite="palegreen"))

Metablites_select <- Zscore[Zscore$Class == "Env",3:44]
head(Metablites_select)

pheatmap(Metablites_select,
         annotation_col = sample_col[,2:3], 
         annotation_row = sample_row[sample_row$Type == "Env",2:3], 
         color = inferno(30),
         annotation_colors = annotation_colors, 
         clustering_distance_r = "manhattan", 
         show_colnames =F,show_rownames = F,
         cluster_rows = T, cluster_cols = T,
         fontsize_row = 4, fontsize_col = 4, 
         angle_col = 90, border_color = NA)
```
#### Plot the metaboliets with Interaction variations
```r
annotation_colors = list(
  treatment = c(WL="red", WW="royalblue"),
  Group = c(Wild="pink", Commercial="slateblue1", Elite="palegreen"))

Metablites_select <- Zscore[Zscore$Class == "Interaction",3:44]
head(Metablites_select)

pheatmap(Metablites_select,
         annotation_col = sample_col[,2:3], 
         annotation_row = sample_row[sample_row$Type == "Interaction",2:3], 
         color = inferno(30),
         annotation_colors = annotation_colors, 
         clustering_distance_r = "manhattan", 
         show_colnames =F,show_rownames = F,
         cluster_rows = T, cluster_cols = T,
         fontsize_row = 4, fontsize_col = 4, 
         angle_col = 90, border_color = NA)
```


### Pair-wise comparison of metaboltes ratio to phenotype ratio 
```r
Ratio <- read.csv("/Users/Leon/OneDrive - Cornell University/3_Manuscript/Manuscript8/7_LogFC_ratio/Metabolites_WWWL_ratio.csv", header = T)
row.names(Ratio) =Ratio$Metabolites
head(Ratio)

LOG <- read.csv("/Users/Leon/OneDrive - Cornell University/3_Manuscript/Manuscript8/7_LogFC_ratio/Ratio_LogFC.csv", header = T)
row.names(LOG) = LOG[,1]
head(LOG)
```
**PM_combine:** Phenotype-metabolites 
```r
###Define DF first for
df <- data.frame(matrix(ncol = 7 , nrow = 451))
colnames(df) <- c("row", "column", "R", "P-value", "Intercept", "Effects", "R.square")

### Loop analysis nested in metabolites
for (n in 1:nrow(LOG)){
  for (i in 1:(nrow(Ratio))){
    res <- rcorr(t(rbind(LOG[n,-1], Ratio[i,-1])))
    
    df[i,1:4] = flattenCorrMatrix(res$r, res$P)
    temp.df = data.frame(t(rbind(LOG[n,-1], Ratio[i,-1])))
    sm = summary(lm(get(LOG[n,1]) ~ get(Ratio[i,1]), temp.df))
    df[i,5] = sm$coefficients[1]
    df[i,6] = sm$coefficients[2]
    df[i,7] = sm$r.squared
  }
  assign(paste0("Metabolites_df",n), df)
}

merge_list <- lapply(ls(pattern = "Metabolites_df"), get)
PM_combine <- bind_rows(merge_list)
PM_combine_clean <- PM_combine[PM_combine$`P-value` < 0.05, ]
head(PM_combine_clean)

write.csv(PM_combine_clean, "/Users/Leon/OneDrive - Cornell University/3_Manuscript/Manuscript8/7_LogFC_ratio/Phenotype-metabolites.correlation-Feb15-2022.csv", quote = F)
```

### Pair-wise comparison of gene expression to metabolites level 
```r
Ratio <- read.csv("/Users/Leon/OneDrive - Cornell University/3_Manuscript/Manuscript8/7_LogFC_ratio/LogFC_meta_table.csv", header = T)
row.names(Ratio) = Ratio$Gene
head(Ratio)

LOG <- read.csv("/Users/Leon/OneDrive - Cornell University/3_Manuscript/Manuscript8/7_LogFC_ratio/Metabolites_WWWL_ratio.csv", header = T)
row.names(LOG) =LOG$Metabolites
head(LOG)

```
**GM_combine:** Genotype-metabolites 
```r
###Define DF first for
df <- data.frame(matrix(ncol = 7 , nrow = 451))
colnames(df) <- c("row", "column", "R", "P-value", "Intercept", "Effects", "R.square")

### Loop analysis nested in metabolites
for (n in 1:nrow(LOG)){
  for (i in 1:(nrow(Ratio))){
    res <- rcorr(t(rbind(LOG[n,-1], Ratio[i,-1])))
    
    df[i,1:4] = flattenCorrMatrix(res$r, res$P)
    temp.df = data.frame(t(rbind(LOG[n,-1], Ratio[i,-1])))
    sm = summary(lm(get(LOG[n,1]) ~ get(Ratio[i,1]), temp.df))
    df[i,5] = sm$coefficients[1]
    df[i,6] = sm$coefficients[2]
    df[i,7] = sm$r.squared
  }
  assign(paste0("GM_df",n), df)
}

merge_list <- lapply(ls(pattern = "GM_df"), get)
GM_combine <- bind_rows(merge_list)
GM_combine_clean <- GM_combine[GM_combine$`P-value` < 0.05, ]
head(GM_combine_clean)

write.csv(GM_combine_clean, "/Users/Leon/OneDrive - Cornell University/3_Manuscript/Manuscript8/LogFC_ratio/Gene-metabolites.correlation-Feb15-2022.csv", quote = F)
```

### Plot effects distribution 
```r
### PLot the density plot of 10 traits and 451 metabolites
combine_filter = read.csv("/Users/Leon/OneDrive - Cornell University/3_Manuscript/Manuscript8/7_LogFC_ratio/Plot_meta_Feb2022.csv", header = T)
combine_filter$abs.effect = abs(combine_filter$Effects)
head(combine_filter)

### define length function
fun_length <- function(x){
  return(data.frame(y=median(x),label= paste0("N=", length(x))))
}
```

```r
### Compare effects across the three types of phenotype (density plot)
ggplot(data = combine_filter, aes(y=log2(abs.effect), fill = Type))+
  geom_density(alpha = 0.2) +
  facet_wrap(~ Trait, ncol=3) + theme_bw() +
  xlab("Log2 transformed Abosolute value of effects") +
  ylab("Density")
```

```r
### Compare R square across thre types of phenotypes
ggplot(data = combine_filter, aes(x=R.square, color = Type))+
  facet_wrap(~ Trait, ncol=3) +
  geom_density(alpha = 0.5) + theme_bw() +
  xlab("R square of linear regression fit") +
  ylab("Density")
```

```r
### Compare effects across the three types of phenotype (boxplot)
level_order = c("Metabolites", "Production", "Index")
ggplot(data = combine_filter, aes(x = factor(Trait, level = level_order), y = log2(abs.effect)))+
  geom_boxplot(alpha = 0.3, outlier.size = 0.02, fill = "grey") +
  theme_bw() +
  ylab("Log2 transformed value of absolute effect") +
  stat_summary(fun = mean, color = "darkred", position = position_dodge(1),
             geom = "point", shape = 4, size = 6,
             show.legend = F) + 
  stat_summary(fun.data = fun_length, geom = "text", vjust = -18, size = 3, color = "red") +
  stat_compare_means(method = "kruskal") +
  theme(legend.position="none")
```

### Plot effects distribution in venn plot
```r
venn_combine <- list(Meatbolites = combine_filter[combine_filter$Trait == "Metabolites",2] , 
                      Index = combine_filter[combine_filter$Trait == "Index",2],
                     Production = combine_filter[combine_filter$Trait == "Production",2]) 

ggVennDiagram(venn_combine, label_alpha = 2, set_color = "black", color = "red")
```

## Genomic variants analysis of 22 accessions using WGS data
**Using the ncbi sra-tool download**
```bash
#!/bin/bash
for sample in $(cat /mnt/Knives/1_data/9_Cotton_omics/5_WGS/SRR_Acc_List.txt); 
do 
    prefetch ${sample} --max-size 200G
    mv /home/liangyu/ncbi/public/sra/${sample}.sra /mnt/Knives/1_data/9_Cotton_omics/5_WGS/${sample}.sra
    fastq-dump --gzip --split-files /mnt/Knives/1_data/9_Cotton_omics/5_WGS/${sample}.sra
done
```
**Using the ENA**
sample list summary
```
/vol1/fastq/SRR631/001/SRR6311571/SRR6311571_1.fastq.gz	/vol1/fastq/SRR631/001/SRR6311571/SRR6311571_2.fastq.gz
/vol1/fastq/SRR631/007/SRR6311717/SRR6311717_1.fastq.gz	/vol1/fastq/SRR631/007/SRR6311717/SRR6311717_2.fastq.gz
/vol1/fastq/SRR631/002/SRR6311732/SRR6311732_1.fastq.gz	/vol1/fastq/SRR631/002/SRR6311732/SRR6311732_2.fastq.gz
/vol1/fastq/SRR631/000/SRR6311780/SRR6311780_1.fastq.gz	/vol1/fastq/SRR631/000/SRR6311780/SRR6311780_2.fastq.gz
/vol1/fastq/SRR631/000/SRR6311500/SRR6311500_1.fastq.gz	/vol1/fastq/SRR631/000/SRR6311500/SRR6311500_2.fastq.gz
/vol1/fastq/SRR631/007/SRR6311857/SRR6311857_1.fastq.gz	/vol1/fastq/SRR631/007/SRR6311857/SRR6311857_2.fastq.gz
/vol1/fastq/SRR631/006/SRR6311856/SRR6311856_1.fastq.gz	/vol1/fastq/SRR631/006/SRR6311856/SRR6311856_2.fastq.gz
/vol1/fastq/SRR631/000/SRR6311570/SRR6311570_1.fastq.gz	/vol1/fastq/SRR631/000/SRR6311570/SRR6311570_2.fastq.gz
/vol1/fastq/SRR631/004/SRR6311744/SRR6311744_1.fastq.gz	/vol1/fastq/SRR631/004/SRR6311744/SRR6311744_2.fastq.gz
/vol1/fastq/SRR631/007/SRR6311867/SRR6311867_1.fastq.gz	/vol1/fastq/SRR631/007/SRR6311867/SRR6311867_2.fastq.gz
/vol1/fastq/SRR631/003/SRR6311863/SRR6311863_1.fastq.gz	/vol1/fastq/SRR631/003/SRR6311863/SRR6311863_2.fastq.gz
/vol1/fastq/SRR631/004/SRR6311864/SRR6311864_1.fastq.gz	/vol1/fastq/SRR631/004/SRR6311864/SRR6311864_2.fastq.gz
/vol1/fastq/SRR631/000/SRR6311820/SRR6311820_1.fastq.gz	/vol1/fastq/SRR631/000/SRR6311820/SRR6311820_2.fastq.gz
/vol1/fastq/SRR631/003/SRR6311853/SRR6311853_1.fastq.gz	/vol1/fastq/SRR631/003/SRR6311853/SRR6311853_2.fastq.gz
/vol1/fastq/SRR631/002/SRR6311572/SRR6311572_1.fastq.gz /vol1/fastq/SRR631/002/SRR6311572/SRR6311572_2.fastq.gz
/vol1/fastq/SRR631/008/SRR6311538/SRR6311538_1.fastq.gz /vol1/fastq/SRR631/008/SRR6311538/SRR6311538_2.fastq.gz
/vol1/fastq/SRR631/002/SRR6311842/SRR6311842_1.fastq.gz /vol1/fastq/SRR631/002/SRR6311842/SRR6311842_2.fastq.gz
```

```bash
#!/bin/bash
work_dir="/mnt/Knives/1_data/9_Cotton_omics/5_WGS"

IFS=$'\n';
for LINE in $(cat $work_dir/sample.list);
do
	pair1=$(echo ${LINE} | awk '{ print $1}')
	pair2=$(echo ${LINE} | awk '{ print $2}')

    ### forward reads download
	ascp -k 1 -QT -l 500m -P33001 \
    	-i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
    	era-fasp@fasp.sra.ebi.ac.uk:$pair1 \
    	/mnt/Knives/1_data/9_Cotton_omics/5_WGS/

    ### reverse reads download
	ascp -k 1 -QT -l 500m -P33001 \
    	-i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
    	era-fasp@fasp.sra.ebi.ac.uk:$pair2\
    	/mnt/Knives/1_data/9_Cotton_omics/5_WGS/
done
```
**merge reads processed by Duke Diversity panel**
NOTE: 
```bash
seq_dir="/mnt/Knives/1_data/9_Cotton_omics/5_WGS"

cat $seq_dir/21552Dpa_SA-1019_S227_L00?_R1_001.fastq.gz > $seq_dir/SA_1019_R1.fastq.gz
cat $seq_dir/21552Dpa_SA-1019_S227_L00?_R2_001.fastq.gz > $seq_dir/SA_1019_R2.fastq.gz

cat $seq_dir/21552Dpa_SA-0159_S245_L00?_R1_001.fastq.gz > $seq_dir/SA_0159_R1.fastq.gz
cat $seq_dir/21552Dpa_SA-0159_S245_L00?_R2_001.fastq.gz > $seq_dir/SA_0159_R2.fastq.gz

cat 21552Dpa_SA-1212_S222_L00?_R1_001.fastq.gz > SA_1212_R1.fastq.gz
cat 21552Dpa_SA-1212_S222_L00?_R2_001.fastq.gz > SA_1212_R2.fastq.gz
```

### Trim and summarize the reads information
```bash
raw_reads_dir="/mnt/Knives/1_data/9_Cotton_omics/5_WGS"
trim_reads_dir="/mnt/Knives/1_data/9_Cotton_omics/5_WGS/Clean_reads"

for i in SRR6311566 SRR6311571 SRR6311717 SRR6311732 SRR6311780 SRR6311500 SRR6311857 SRR6311856 SRR6311570; 
do
 trimmomatic PE -threads 10 \
    ${raw_reads_dir}/${i}_R1.fastq.gz \
    ${raw_reads_dir}/${i}_R2.fastq.gz \
    ${trim_reads_dir}/${i}_forward_paired.fq.gz ${trim_reads_dir}/${i}_forward_unpaired.fq.gz \
    ${trim_reads_dir}/${i}_reverse_paired.fq.gz ${trim_reads_dir}/${i}_reverse_unpaired.fq.gz \
    ILLUMINACLIP:/home/liangyu/anaconda3/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:2 LEADING:3 TRAILING:3 MINLEN:36

   ### check the numbers of reads within the unpaired reads file
   echo $(zcat ${trim_reads_dir}/${i}_forward_paired.fq.gz | wc -l)/4|bc >> ${trim_reads_dir}/forward_clean.count
   echo $(zcat ${trim_reads_dir}/${i}_reverse_paired.fq.gz | wc -l)/4|bc >> ${trim_reads_dir}/reverse_clean.count

   echo $(zcat ${trim_reads_dir}/${i}_forward_unpaired.fq.gz | wc -l)/4|bc >> ${trim_reads_dir}/forward_unpair.count
   echo $(zcat ${trim_reads_dir}/${i}_reverse_unpaired.fq.gz | wc -l)/4|bc >> ${trim_reads_dir}/reverse_unpair.count
done

### move all untrimmed reads as archives
mv /mnt/Knives/1_data/Sorghum_Pangenome/1_Reads/*.fastq.gz /mnt/Kanta/Sorghum_reads/
```
**Summarize reads info**
```bash
trim_reads_dir="/mnt/Knives/1_data/9_Cotton_omics/5_WGS/Clean_reads"

for i in $(cat ${trim_reads_dir}/clean_sample);
do 
  echo $(zcat ${trim_reads_dir}/${i} | wc -l)/4|bc >> ${trim_reads_dir}/Pannel_reads.count
  seqkit fx2tab -nl ${trim_reads_dir}/${i} | head -1 
done
```

### Variants calling of 18 accessions using (DNA-seq data)
Note from 20220726: The SRR6311717 sample was interrupted and needs to be re-launched
Sample list of all accesions

```
SRR6311566
SRR6311571
SRR6311717
SRR6311732
SRR6311780
SRR6311500
SRR6311857
SRR6311856
SRR6311570
SRR6311744
SRR6311842
SRR6311538
SRR6311572
SRR6311867
SRR6311863
SRR6311864
SRR6311820
SRR6311853
```

```bash
Ref_dir="/mnt/Knives/1_data/9_Cotton_omics/0_data"
workdir="/mnt/Knives/1_data/9_Cotton_omics/6_Varaints" 
BAM_dir="/mnt/Knives/1_data/9_Cotton_omics/5_WGS/Clean_reads" 
GATK_dir="/home/liangyu/anaconda3/pkgs/gatk-3.8-hdfd78af_11/opt/gatk-3.8"

#Build index and dict
samtools faidx $Ref_dir/Ghirsutum_527_v2.0.fa
bwa index $Ref_dir/Ghirsutum_527_v2.0.fa 
picard CreateSequenceDictionary -R $Ref_dir/Ghirsutum_527_v2.0.fa -O $Ref_dir/Ghirsutum_527_v2.0.dict

### set the thread numbers
nt=14

#loop settings
for file in $(cat $workdir/sample.list);
do
  mkdir ${workdir}/${file}

  ### Map reads by bwa piepline
  bwa mem -M -R "@RG\tID:$file\tSM:$file\tPL:ILLUMINA" \
    -t $nt  \
   $Ref_dir/Ghirsutum_527_v2.0.fa \
    ${BAM_dir}/${file}_forward_paired.fq.gz \
    ${BAM_dir}/${file}_reverse_paired.fq.gz | samtools view -@10 -bS - -o ${workdir}/${file}/${file}.reorder.bam 

  samtools view -b -F 4 ${workdir}/${file}/${file}.reorder.bam > ${workdir}/${file}/${file}.cut.bam

  ### Sort bam files
  picard SortSam \
   INPUT=${workdir}/${file}/${file}.cut.bam \
   OUTPUT=${workdir}/${file}/${file}.sort.bam SORT_ORDER=coordinate

  ### Mark and remove duplictaes from Bam
  picard MarkDuplicates -REMOVE_DUPLICATES true \
    -I ${workdir}/${file}/${file}.sort.bam \
    -O ${workdir}/${file}/${file}.sort.dup.bam \
    -M ${workdir}/${file}/${file}.metrics.txt 

  ### add group numbers and infotrmation
  picard AddOrReplaceReadGroups \
   -I ${workdir}/${file}/${file}.sort.dup.bam \
   -O ${workdir}/${file}/${file}.sort.dup.Gr.bam -LB cotton -PL illumina -PU Diversity -SM ${file} &

  ###index bam file
  cd ${workdir}/${file}
  samtools index ${workdir}/${file}/${file}.sort.dup.Gr.bam

  java -Xmx512G -jar $GATK_dir/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -R $Ref_dir/Ghirsutum_527_v2.0.fa \
    -I ${workdir}/${file}/${file}.sort.dup.Gr.bam \
    -nct $nt --emitRefConfidence GVCF \
    -variant_index_type LINEAR \
    -variant_index_parameter 128000 \
    -o ${workdir}/${file}/${file}.gvcf

  mv ${workdir}/${file}/${file}.gvcf ${workdir}/
done

ls ${workdir} | grep ".gvcf" | sort > gvcf.list

### Generate the cohort VCFfiles
java -Xmx512G -jar $GATK_dir/GenomeAnalysisTK.jar \
 -T CombineGVCFs -nct $nt \
 -R $Ref_dir/Ghirsutum_527_v2.0.fa \
 --variant ${workdir}/gvcf.list 
 -o cotton.resequencing.gvcf \

### Transform the VCFs format
java -Xmx512G -jar $GATK_dir/GenomeAnalysisTK.jar \
 -T GenotypeGVCFs -nct $nt \
 -R $Ref_dir/Ghirsutum_527_v2.0.fa \
 -V $workdir/cotton.resequencing.gvcf \
 -o $workdir/cotton_DNA.vcf \

```
### Reformat the VCF with index information added 

```bash
### Extract the header information 
grep "#" > $workdir/header

### reformat the index and combine with the header
grep -v "#" $workdir/cotton_DNA.vcf  | awk '{print $1,$2,$1"_"$2,$0}' | sed 's/ /\t/g' | cut -f1-3,7- | cat header - > cotton_DNA2.vcf
```
Or use the Bcftools to edit:
```bash
bcftools annotate --set-id +'%CHROM\_%POS' cotton_DNA.vcf > cotton_setID.vcf
```
See below as a example before and after:
```
A01     1022    .       C       T       6192.83 .      
A01     1225    .       G       C       3539.09 .     
A01     1331    .       G       A       429.12  .      
A01     1493    .       A       G       3456.76 .      
```
```
A01     1022    A01_1022        C       T       6192.83 .       
A01     1225    A01_1225        G       C       3539.09 .       
A01     1331    A01_1331        G       A       429.12  .      
A01     1493    A01_1493        A       G       3456.76 .      
A01     1535    A01_1535        A       G       4001.20 .       
A01     1704    A01_1704        T       G       840.12  .      
A01     1707    A01_1707        A       C       813.36  . 
```

### Filter the variants of multiple samples
Check link for defination of [**parameter**](https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/QUAL_QD_GQ_Formulation_fDG.htm) 


**DP is determined by average total depth (total depth/population size = individual sample)** \
**minGQ=20**: 1% chance that the call is incorrect** \
**minQUAL=20**: 1 % chance that there is no variant at the site**

```bash
#!/bin/bash
workdir="/data/home/lyu/01_project/02_Cotton/02_Variants/01_Filter"
scripts_dir="/data/home/lyu/02_script/R"

### remove loci with missing data based on cut-off to reduce bias of total mapping depth
#keep loci with up to 10% missing data for plot and minor allel frquency 
vcftools --vcf $workdir/cotton_setID.vcf \
	--max-missing 0.90 \
	--min-alleles 2 \
	--recode \
	--maf 0.05 \
	--recode-INFO-all \
	--out $workdir/cotton_plot.vcf

### PLOT DP, QD, and QUAL for each loci 
###Extrat VCF depth information (This is the combined depth among samples)
egrep -v "^#" $workdir/cotton_plot.vcf.recode.vcf | \
cut -f 8 | \
sed 's/^.*;DP=\([0-9]*\);.*$/\1/' > $workdir/cotton_map_DP.txt

###Extrat VCF QualitybyDepth information 
egrep -v "^#" $workdir/cotton_plot.vcf.recode.vcf | \
cut -f 8 | \
sed 's/^.*;QD=\([0-9]*.[0-9]*\);.*$/\1/' > $workdir/cotton_map_QD.txt

###Extrat VCF Quality information
egrep -v "^#" $workdir/cotton_plot.vcf.recode.vcf | \
cut -f 6 > $workdir/cotton_map_QUAL.txt

### Plot three component
Rscript ${scripts_dir}/Plot_VCF.R \
  --DP $workdir/cotton_map_DP.txt \
  --QUAL $workdir/cotton_map_QUAL.txt \
  --QD $workdir/cotton_map_QD.txt 
```
```
**Results:
SNP number round 1: After filtering, kept 5,849,690 out of a possible 7,245,555 Sites

Density distribution: (percentile summary)
    10% 20% 30% 40% 50% 60% 70% 80% 90%
DP: 243 258 267 275 283 291 299 311 333

        10%      20%      30%      40%      50%      60%      70%      80%
QUAL: 264.300  456.250  628.240  877.276 1303.540 2076.480 3120.650 4335.320

    10%   20%   30%   40%   50%   60%   70%   80%   90%
QD 5.66 18.55 25.57 27.67 29.11 30.32 31.51 32.71 33.96

SNP number round 2: After filtering, kept 5416973 out of a possible 5849690 Sites


```
### Remove loci based on cutoff derived from plot information 

```bash
#!/bin/bash
workdir="/data/home/lyu/01_project/02_Cotton/02_Variants/01_Filter"

### fileter based on GQ score and QUAL
vcftools --vcf cotton_plot.vcf.recode.vcf  \
	--minDP 250 \
	--minGQ 20 \
	--minQ 200 \
	--recode \
	--recode-INFO-all \
	--out cotton_clean.vcf

### filter based qualitybydepth (QD) (QD is the normalized quality based on mapping depth usually set as 2)
java -jar /home/lyu/mambaforge/share/biopet-vcffilter-0.2-1/vcffilter-assembly-0.2.jar \
  -i cotton_clean.vcf.recode.vcf \
  -o output.vcf \
  --filter "QD < 5"

## Line NO. 1079049

### extract genotype information only 
cat $workdir/F1_final_filter.vcf | perl while (<>) { s!:\S+!!g; print; } > cotton_final_filter.clean 
```
### Perform the LD-prune to reduce marker density 
```bash
### build the lib
plink --vcf cotton_clean.vcf.recode.vcf \
  --make-bed \
  --allow-extra-chr \
  --out cotton_bfile

### Use LD and distance to prune
plink --bfile cotton_bfile \
  --indep-pairwise 10 5 0.8 \
  --allow-extra-chr \
  --out Cotton_pruned

### Extract the pruned loci 
bcftools view \
  --include ID==@Cotton_pruned.prune.in \
  cotton_clean.vcf.recode.vcf > cotton.prune.vcf

```

### Annotate variants using SnpEff
**Configure and install Snpeff** \
**NOTE:** modify the configuration file: snpeff.config
```
# SnpEff configuration file
# Databases are stored here
# E.g.: Information for 'hg19' is stored in data.dir/hg19/
# You can use tilde ('~') as first character to refer to your home directory. 
# Also, a non-absolute path will be relative to config's file dir

data.dir = /home/lyu/03_software/snpEff/data
```
**Run Snpeff**
```bash
#!/bin/bash

snp_dir="/home/lyu/03_software/snpEff"
genome_dir="/data/home/lyu/01_project/02_Cotton/01_Genome"
cd $snp_dir

### linke reference genome
mkdir $snp_dir/data/Cotton
ln -s $genome_dir/Ghirsutum_527_v2.1.gene.gff3 $snp_dir/data/Cotton/genes.gff

### build libraries
snpEff -Xmx512g build -gff3 Cotton -c snpEff.config

### perform annotation
snpEff  -Xmx512g \
    Cotton ${vcf_dir}/cotton_clean.vcf.recode.vcf > ${snp_dir}/Cotton.ann.vcf
```

### Annotate variants using VEP (Variant Effect Predictor) 
**NOTE**: (Have to use the gff3 with the exon information)

```bash
genome_dir="/home/liangyu/1_Project/5_Cotton/0_data"

### prepare tehe VCF and ref genome file
mkdir $HOME/vep_data
ln -s $genome_dir/Ghirsutum_527_v2.1.gene.gff3 $HOME/vep_data
scp -P 22 lyu@ghibli.bti.cornell.edu:/home/lyu/01_project/02_Cotton/02_Variants/01_Filter/cotton_clean.vcf.recode.vcf $HOME/vep_data/

### Reformat the gff3
grep -v "#" Ghirsutum_527_v2.1.gene.gff3 | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > Gh.gff.gz
tabix -p gff Gh.gff.gz
```
**Using Singularity**
```bash
### Perform the annotation 
singularity exec $HOME/6_Docker/vep.sif \
    vep --dir $HOME/vep_data \
    -i $HOME/vep_data/cotton_clean.vcf.recode.vcf \
    --gff $HOME/vep_data/Gh.gff.gz \
    --fasta $HOME/vep_data/Ghirsutum_527_v2.0.fa
```
**Using docker**
```bash
docker pull ensemblorg/ensembl-vep:release_106.1

sudo docker run --rm -v $(pwd):/working-dir -w /working-dir ensemblorg/ensembl-vep \
    vep --dir $HOME/vep_data \
    -i $HOME/vep_data/cotton_clean.vcf.recode.vcf \
    --gff $HOME/vep_data/Gh.gff.gz \
    --fasta $HOME/vep_data/Ghirsutum_527_v2.0.fa
```

### Extract the variants annotation of interested genes
```bash
 for i in $(cat /home/liangyu/vep_data/01_Results/GeneList);
 do 
  ### grep genes
  grep $i Cotton_effect_output.txt >> Select_Varaints;

done
```

### Build the phylogeny using IQtree
**NOTE** THE first row of phylip file reprsesents the outgroup of TREE!
```bash
#!/bin/bash

script_dir="/data/home/lyu/03_software/vcf2phylip"
vcf_dir="/data/home/lyu/01_project/02_Cotton/04_Phylogeny"

### Convert vcf to phylip
python $script_dir/vcf2phylip.py -i $vcf_dir/cotton.prune.vcf

### REPLACE "*" into "N" also replace the sequnece ID into certain sample name
IFS=$'\n';
for LINE in $(cat $work_dir/ID_table | grep -v 'ID');do

    #READ DIRECTORIES
    ID=$(echo ${LINE} | awk '{ print $1}')
    Accession=$(echo ${LINE} | awk '{ print $2 }')

    #Repalace ID as accessions names
    sed -i "s/$ID/$Accession/g" $work_dir/cotton.prune.min4.phy.treefile 
done

sed -i 's/*/N/g' cotton.prune.min4.phy 

### IQ-tree construction 
iqtree -s cotton.prune.min4.phy -nt 125 -m MFP
```

### Extract the SNPs/INDELs data for IGV visulization purpose
```bash
### modify the header (sample name)
sed -i 's/SRR6311566/AC_134_CB_4029/g' cotton_setID.vcf
sed -i 's/SRR6311571/CB_4012_KING_KARAJAZSY/g' cotton_setID.vcf
sed -i 's/SRR6311717/FELISTANA_UA_7_18/g' cotton_setID.vcf
sed -i 's/SRR6311732/FUNTUA_FT_5/g' cotton_setID.vcf
sed -i 's/SRR6311780/LOCKETT_BXL/g' cotton_setID.vcf
sed -i 's/SRR6311500/TAM_86III_26/g' cotton_setID.vcf
sed -i 's/SRR6311857/TAM_91C_34/g' cotton_setID.vcf
sed -i 's/SRR6311856/TAM_94WE_37S/g' cotton_setID.vcf
sed -i 's/SRR6311570/COKER_310/g' cotton_setID.vcf
sed -i 's/SRR6311744/PD3/g' cotton_setID.vcf
sed -i 's/SRR6311842/TIPO_CHACO_UA_4_4/g' cotton_setID.vcf
sed -i 's/SRR6311538/WESTERN_STOOMPROOF/g' cotton_setID.vcf
sed -i 's/SRR6311572/MEXICO_910/g' cotton_setID.vcf
sed -i 's/SRR6311867/VIR_7094_COKER_310/g' cotton_setID.vcf
sed -i 's/SRR6311863/VIR_7223/g' cotton_setID.vcf
sed -i 's/SRR6311864/VIR_7153_D_10/g' cotton_setID.vcf
sed -i 's/SRR6311820/PLAINS/g' cotton_setID.vcf
sed -i 's/SRR6311853/VIR_6615_MCU_5/g' cotton_setID.vcf

### compress files
bgzip -c cotton_setID.vcf > cotton_setID.vcf.gz
 
### Index the vcf.gz
bcftools index cotton_setID.vcf.gz

for i in 01 02 03 04 05 06 07 08 09 10 11 12 13;
do
 
  ### Extract the regions of interests
  bcftools view cotton_setID.vcf.gz --regions D$i > D$i.vcf 
  bcftools view cotton_setID.vcf.gz --regions D$i > D$i.vcf 

  ### Index the VCF file
  bcftools index D$i.vcf.gz
  bcftools index A$i.vcf.gz 

done
```


## Cistrome Deep learning-based prediction
### Extract the positive and negative sequences for training 
See R notebook regarding extracting coordinates of each narrowPeaks
Selected TFs used for model traning:

```
HSFA6B_colamp
HSFA6B_col
AT5G60130_col_a
```

```bash
#!/bin/bash
work_dir="/home/lyu/01_project/02_Cotton/05_cistrome/01_PeakIndex"

for FILE in HSF7_colamp_a HSF7_col_a HSF6_col_a bZIP44_col_a HSFA6B_colamp HSFA6B_col;
do
  sed -i 's/chr/Chr/g' ${FILE}_PeakIndex

  IFS=$'\n';
  for LINE in $(cat ${work_dir}/${FILE}_PeakIndex | grep -v 'index');
  do
    seq1=$(echo ${LINE} | awk '{ print $1}')
    seq2=$(echo ${LINE} | awk '{ print $2}')
    seq3=$(echo ${LINE} | awk '{ print $3}')

    ### extract sequences using samtools
    samtools faidx $work_dir/Athaliana_447_TAIR10.fa $seq1 >> ${work_dir}/${FILE}_positive_intact.fasta
    samtools faidx $work_dir/Athaliana_447_TAIR10.fa $seq2 >> ${work_dir}/${FILE}_negative.fasta
    samtools faidx $work_dir/Athaliana_447_TAIR10.fa $seq3 >> ${work_dir}/${FILE}_negative.fasta
  done
   
  ### Generate the reverse complement sequences for peak summit
  revseq -sequence ${work_dir}/${FILE}_positive_intact.fasta -outseq ${work_dir}/${FILE}_positive_RC.fasta
  
  ### Combine sequences
  cat ${work_dir}/${FILE}_positive_intact.fasta ${work_dir}/${FILE}_positive_RC.fasta > ${work_dir}/${FILE}_positive.fasta

  ### reformat sequences
  grep -v ">" ${work_dir}/${FILE}_positive.fasta > ${work_dir}/${FILE}_positive.txt
  grep -v ">" ${work_dir}/${FILE}_negative.fasta > ${work_dir}/${FILE}_negative.txt

  rm ${work_dir}/${FILE}_positive.fasta
  rm ${work_dir}/${FILE}_negative.fasta
done
```
### Extract the promoter sequences (reverse and forward sequences)
**Note: remove the sequences which are less than 2K** \
For reverse sequences:
```bash
### extract the prmoter regions
for i in $(cat /home/liangyu/vep_data/Gh_reverse_promoter);
do
  samtools faidx /home/liangyu/vep_data/Ghirsutum_527_v2.0.fa $i >> Gh_reverse_promoter.fasta
done

### Exclude the short sequences
faSomeRecords Gh_reverse_promoter.fasta -exclude exc_seq Gh_reverse_promoter_exc.fasta
sed -i 's/:/_/g' Gh_reverse_promoter_exc.fasta
sed -i 's/-/_/g' Gh_reverse_promoter_exc.fasta

### reverse complement the sequences for reserve genes 
revseq -sequence Gh_reverse_promoter_exc.fasta -outseq Gh_reverse_promoter_RC.fasta
```

For forward sequences:
```bash
### extract the prmoter regions
for i in $(cat /home/liangyu/vep_data/Gh_forward_promoter);
do
  samtools faidx /home/liangyu/vep_data/Ghirsutum_527_v2.0.fa $i >> Gh_forward_promoter.fasta
done

sed -i 's/:/_/g' Gh_forward_promoter.fasta
sed -i 's/-/_/g' Gh_forward_promoter.fasta

### Combine sequences together
cat Gh_forward_promoter.fasta Gh_reverse_promoter_RC.fasta > Gh_promoter.fasta
```

reformat the pairID before replacement under the Gh_promoterName
```
A01_14115_16115	Gohir.A01G000101
A01_46776_48776	Gohir.A01G000301
A01_50653_52653	Gohir.A01G000401
A01_79810_81810	Gohir.A01G000800
A01_84890_86890	Gohir.A01G000900
A01_87047_89047	Gohir.A01G001000
A01_92702_94702	Gohir.A01G001100
A01_96507_98507	Gohir.A01G001200
A01_103286_105286	Gohir.A01G001400
A01_112107_114107	Gohir.A01G001600
```
```bash
### Using seqkit to replace the ID
seqkit replace -k Gh_promoterName --pattern "^(\w+)" --replacement "{kv}" -w 0 Gh_promoter.fasta > Gh_promoter_clean.fasta

```

### Run DL model tranining and CREs prediction

Note: install the following version of sklearn for compatability issue
```bash
pip install scikit-learn==0.24.2
```
Check the version after installation
```python
import sklearn
print(sklearn.__version__)
```
Reformat the multiple line fasta to one line fasta sequence 
```bash
perl -pe '/^>/ ? print "\n" : chomp' in.fasta > out.fasta
```
See [**tutorial**](https://github.com/Takeshiddd/CisDecoding_cistrome) for details:
```bash
#!/bin/bash
work_dir="/home/lyu/01_project/02_Cotton/05_cistrome/01_PeakIndex"
ROC_dir="/home/lyu/01_project/02_Cotton/05_cistrome/02_Model"
model_dir="/home/lyu/01_project/02_Cotton/05_cistrome/models/home/lyu/01_project/02_Cotton/05_cistrome/02_Model"
output_dir="/home/lyu/01_project/02_Cotton/05_cistrome/03_Results"

for FILE in HSFA6B_colamp HSFA6B_col
do
  ### Training/validation by a fully-connected model
  python FC-cistrome-training.py \
      -p ${work_dir}/${FILE}_positive.txt \
      -n ${work_dir}/${FILE}_negative.txt \
      -o $ROC_dir/${FILE} \
      -e 20 \
      -l 31 

  ### Detection of prediction (of each TF biding) in the target promoter sequences in bin sliding-windows
  python MultiSeq_CREs_prediction_walking.py \
      -f Turquoise_singleLine.fasta \
      -m ${model_dir}/${FILE}.h5 \
      -w 4 \
      -b 31 \
      -o ${output_dir}/${FILE}

  ### Conversion to binary CRE arrays in bins
  python BinIntg2BinaryArray.py \
      -i ${output_dir}/wOTU-${FILE} \
      -b 2 \
      -t 0.8 \
      -o ${output_dir}/${FILE}-bin

done
```

## DAP-seq data analysis


### Data preparation
Create the sample list files for data prep
```
Gohir.A08G064100_A
Gohir.A08G064100_B
Gohir.A11G200600_A
Gohir.A11G200600_B
Gohir.A13G021700_A
Gohir.A13G021700_B
Gohir.D08G072600_A
Gohir.D08G072600_B
Gohir.D11G114200_A
Gohir.D11G114200_B
Gohir.D13G022300_A
Gohir.D13G022300_B
Halo_A
Halo_B
```
```bash
#! /bin/bash
raw_reads_dir="/mnt/Leonof/Xiaodan/DAPseq_cotton/0_rawData/resequencing"
trim_reads_dir="/mnt/Leonof/Xiaodan/DAPseq_cotton/1_trimData"

## when trimming with trim galore, need to specify whether it's paired reads.

for i in $(cat ${raw_reads_dir}/sample.list); 
do
 ~/Software/TrimGalore/trim_galore \
  --phred33 \
  --fastqc \
  --cores 10 \
  --path_to_cutadapt ~/.local/bin/cutadapt \
  --output_dir /mnt/Leonof/Xiaodan/DAPseq_cotton/1_trimData \
  ${raw_reads_dir}/${i}.fastq.gz

done

```

### Mapping reads to the reference genome
```bash
#!/bin/bash
Ref_dir="/mnt/Knives/1_data/9_Cotton_omics/0_data"
trim_reads_dir="/mnt/Leonof/Xiaodan/DAPseq_cotton/1_trimData" 
workdir="/mnt/Leonof/Xiaodan/DAPseq_cotton/2_mapping_bwa" 

#Build index and dict
samtools faidx $Ref_dir/Ghirsutum_527_v2.0.fa
bwa index $Ref_dir/Ghirsutum_527_v2.0.fa 
picard CreateSequenceDictionary -R $Ref_dir/Ghirsutum_527_v2.0.fa -O $Ref_dir/Ghirsutum_527_v2.0.dict

### set the thread numbers
nt=14

#loop settings
for file in $(cat ${trim_reads_dir}/sample.list);
do
  #mkdir ${workdir}/${file}

  ## Map reads by bwa piepline
  bwa mem -M -B 40 -R "@RG\tID:$file\tSM:$file\tPL:ILLUMINA" \
    -t $nt  \
    $Ref_dir/Ghirsutum_527_v2.0.fa \
    ${trim_reads_dir}/${file}_trimmed.fq.gz | samtools view -@10 -bS - -o ${workdir}/${file}.reorder.bam 

  ## uniquely mapped reads 
  sambamba view -t 15 -h -f bam \
    -F "mapping_quality >= 1 and not (unmapped or secondary_alignment) and not ([XA] != null or [SA] != null)" \
    ${workdir}/${file}.reorder.bam \
    -o ${workdir}/${file}.uniqmapped.bam

done
```

### Calling of DAPseq peaks
```bash
#! /bin/bash
mapped_reads_dir="/mnt/Leonof/Xiaodan/DAPseq_cotton/2_mapping_bwa"

for i in $(cat ${mapped_reads_dir}/sample.list); 
do
  
 macs2 callpeak \
 -t ${mapped_reads_dir}/${i}.uniqmapped.bam \
 -c ${mapped_reads_dir}/Halo_A.uniqmapped.bam \
 -n ${i} \
 -f BAM \
 --keep-dup all \
 -g 2.3e9 \

done

```
### Filter the high-confidence peaks
```bash
##
~/Software/idr-2.0.2/bin/idr --samples \
../3_peakcalling/Gohir.D08G072600_A_combined_peaks.narrowPeak \
  ../3_peakcalling/Gohir.D08G072600_B_combined_peaks.narrowPeak \
  --input-file-type narrowPeak --rank q.value --output-file Gohir.D08G072600-idr --plot --log-output-file Gohir.D08G072600.idr.log  

idr --samples \
  ../3_peakcalling/Gohir.A13G021700_A_combined_peaks.narrowPeak \
  ../3_peakcalling/Gohir.A13G021700_B_combined_peaks.narrowPeak \
  --input-file-type narrowPeak --rank q.value --output-file Gohir.A13G021700-idr --plot --log-output-file Gohir.A13G021700.idr.log
```


### Intersect Variants data with peak profiles
1. The peak center where fall into the **Upstream promoter** (-5000bp ~ 0 TSS) and **Downstream promoter** (0 TSS ~ Longest 5‘UTR) regions wers selected as \
potential regulatory elements of two TFs. \
2. Variants calling data (SNPs and INDELs) which associated with selected core genes within modules were slected to intersect with peaks filtered above

```bash
work_dir="/mnt/Knives/1_data/9_Cotton_omics/10_Haplotype"

### use intersection for two bed files
	intersectBed -a ${work_dir}/HSFA6B-Peak.bed \ 
    -b ${work_dir}/CottonVariants.bed -wo  > ${work_dir}/HSFA6B-Variants-Intersect.bed

  intersectBed -a ${work_dir}/HSFA6B-Peak.bed \ 
    -b ${work_dir}/CottonVariants.bed -wo  > ${work_dir}/HSFA6B-Variants-Intersect.bed
```

### DAP-seq data reshape and extraction
#### Install packages
```r
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ChIPseeker")
BiocManager::install("clusterProfiler")
```

#### Load libraries
```r
setwd("/mnt/Leonof/Xiaodan/DAPseq_cotton/5_PeakAnno")

library(ChIPseeker) 
library(clusterProfiler)
library(AnnotationDbi)
require(GenomicFeatures)
```

#### Peak annotation step
```r
TDXB <- makeTxDbFromGFF("/home/liangyu/vep_data/Ghirsutum_527_v2.1.gene_exons.gff3")

Gohir.A13G021700_peaks <- readPeakFile("/mnt/Leonof/Xiaodan/DAPseq_cotton/4_peakIDR/Gohir.A13G021700-idr")
Gohir.A13G021700_R1_peaks <- readPeakFile("/mnt/Leonof/Xiaodan/DAPseq_cotton/3_peakcalling/Gohir.A13G021700_A_peaks.narrowPeak")
Gohir.A13G021700_R2_peaks <- readPeakFile("/mnt/Leonof/Xiaodan/DAPseq_cotton/3_peakcalling/Gohir.A13G021700_B_peaks.narrowPeak")

Gohir.D08G072600_peaks <- readPeakFile("/mnt/Leonof/Xiaodan/DAPseq_cotton/4_peakIDR/Gohir.D08G072600-idr")
Gohir.D08G072600_R1_peaks <- readPeakFile("/mnt/Leonof/Xiaodan/DAPseq_cotton/3_peakcalling/Gohir.D08G072600_A_peaks.narrowPeak")
Gohir.D08G072600_R2_peaks <- readPeakFile("/mnt/Leonof/Xiaodan/DAPseq_cotton/3_peakcalling/Gohir.D08G072600_B_peaks.narrowPeak")

### annotate and set TSS ranges
peakAnno_DREB2C_R1 <- annotatePeak(Gohir.A13G021700_R1_peaks, tssRegion = c(-3000, 3000), TxDb = TDXB)
peakAnno_DREB2C_R2 <- annotatePeak(Gohir.A13G021700_R2_peaks, tssRegion = c(-3000, 3000), TxDb = TDXB)

peakAnno_HSFA6B_R1 <- annotatePeak(Gohir.D08G072600_R1_peaks, tssRegion = c(-3000, 3000), TxDb = TDXB)
peakAnno_HSFA6B_R2 <- annotatePeak(Gohir.D08G072600_R2_peaks, tssRegion = c(-3000, 3000), TxDb = TDXB)

### write the resul
write.table(peakAnno_DREB2C_R1, file = '/home/liangyu/peakAnno_DREB2C_R1_peaks.anno.txt',sep = '\t', quote = FALSE, row.names = FALSE)
write.table(peakAnno_DREB2C_R2, file = '/home/liangyu/peakAnno_DREB2C_R2_peaks.anno.txt',sep = '\t', quote = FALSE, row.names = FALSE)
write.table(peakAnno_HSFA6B_R1, file = '/home/liangyu/peakAnno_HSFA6B_R1_peaks.anno.txt',sep = '\t', quote = FALSE, row.names = FALSE)
write.table(peakAnno_HSFA6B_R2, file = '/home/liangyu/peakAnno_HSFA6B_R2_peaks.anno.txt',sep = '\t', quote = FALSE, row.names = FALSE)

```

#### Check the log2FC of genes been bound by DREB2A and HSFA6B
**Comparison** of two classes of genes will be performed: \
1. Genes bound by 2 TFs from lint yield-related module \
2. Genes bound by 2 TFs from the transcriptome-wide

#### Load the list of genes from two modules
```r
Plot_target <- read.csv("/Users/Leon/Documents/Project/Cotton/TF-targtsGenes.csv", header =T)
LogFC <- read.csv("/Users/Leon/Documents/Project/Cotton/Gh_logFC.csv", header = T)

unique(Plot_target$Class)

logFC_table <- left_join(Plot_target, LogFC, by = "Gene")
logFC_table$log2FC_ave <- rowMeans(logFC_table[,3:23])

### random pick rows from global level 

HSFA6B_Global <- logFC_table[logFC_table$Class == "HSFA6B_Global",]
HSFA6B_select <- HSFA6B_Global[sample(nrow(HSFA6B_Global), 200), ]

DREB2A_Global <- logFC_table[logFC_table$Class == "DREB2A_Global",]
DREB2A_select <- DREB2A_Global[sample(nrow(DREB2A_Global), 100), ]

HSFA6B_local <- logFC_table[logFC_table$Class == "HSFA6B_targets",]
DREB2A_local <- logFC_table[logFC_table$Class == "DREB_Targets",]

Plot_target <- bind_rows(HSFA6B_select, DREB2A_select, DREB2A_local, HSFA6B_local)
```

```r
### Relevel factors
Plot_target $Class <- factor(Plot_target$Class, levels = c("HSFA6B_Global", "HSFA6B_targets", "DREB2A_Global", "DREB_Targets")) 

my_comparisons <- list( c("DREB_Targets", "DREB2A_Global"), c("HSFA6B_targets", "HSFA6B_Global"))

### Plot boxplot over four conditions 
ggplot(Plot_target, 
       (aes(x=Class, y=log2FC_ave, fill = Class))) +
  geom_boxplot(outlier.size = 0.2, outlier.shape = 21) +
  theme_classic() +
  ylab("Log2FC (WL/WW) of TPM")+
  scale_fill_manual(values = c(DREB_Targets="#a38a77", HSFA6B_targets="#a38a77", HSFA6B_Global="#0ea29a", DREB2A_Global="#0ea29a" )) +
  stat_summary(fun = mean, color = "darkred", position = position_dodge(0.75),
               geom = "point", shape = 4, size = 3,
               show.legend = FALSE) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") +
  theme(legend.position = "none")         

```

## Population genomics: SNP analysis of IPS flanking SNPs
### BYU panel with deep sequencing data

```bash
vcf_dir="/data/home/lyu/01_project/03_DIV/07_VCF-BYU"
Ref_dir="/data/home/lyu/01_project/03_DIV/02_Genome"
work_dir="/data/home/lyu/01_project/03_DIV/08_BYU-VCFs-select"

for i in $(cat $work_dir/BYU_list);
do
  gatk SelectVariants \
      -R ${Ref_dir}/GhirsutumCoker_698_v1.0.fa \
      -V ${vcf_dir}/${i}.g.vcf.gz \
      -L A02:98843947-98843832 \
      -O ${work_dir}/${i}.g.select.vcf.gz
done
```

### TAMU panel with deep sequencing data

```bash
vcf_dir="/data/home/lyu/01_project/03_DIV/07_VCF"
Ref_dir="/data/home/lyu/01_project/03_DIV/02_Genome"
work_dir="/data/home/lyu/01_project/03_DIV/09_TAMU-VCFs-select"

for i in $(cat $work_dir/BYU_list);
do
  gatk SelectVariants \
      -R ${Ref_dir}/GhirsutumCoker_698_v1.0.fa \
      -V ${vcf_dir}/${i}.g.vcf.gz \
      -L A02:98843947-98843832 \
      -O ${work_dir}/${i}.g.select.vcf.gz
done
```
Combined gVCFs using the selected gVCFs
```bash
Ref_dir="/data/home/lyu/01_project/03_DIV/02_Genome"
work_dir="/data/home/lyu/01_project/03_DIV/08_BYU-VCFs-select"
sample_dir="/data/home/lyu/01_project/03_DIV/08_BYU-VCFs-select"

gatk CombineGVCFs \
  --java-options  "-Xmx512g -XX:ParallelGCThreads=100" \
  -R ${Ref_dir}/GhirsutumCoker_698_v1.0.fa \
  --variant ${sample_dir}/gVCF.list \
  -O ${work_dir}/Cotton-IPS-BYU-panel.gvcf.gz


gatk GenotypeGVCFs
 -R ${Ref_dir}/GhirsutumCoker_698_v1.0.fa \
 -V ${work_dir}/Cotton-IPS-BYU-panel.gvcf.gz \
 -O ${work_dir}/Cotton-IPS-BYU-panel.vcf
```

### Using the published larger panel for analysis
The [**Advanced Science paper**](https://doi.org/10.1002/advs.202003634) containts the lager population panel for extraction and we downloaded \
the [**genome V1.1**](https://data.jgi.doe.gov/refine-download/phytozome?organism=Ghirsutum&expanded=458&_gl=1*112jyu9*_ga*ODA0ODUzNTk5LjE2ODc4MDQ1MTY.*_ga_YBLMHYR3C2*MTY5NDQ1NDIyMi4xNy4xLjE2OTQ0NTQyNTguMC4wLjA.) for to search corresponded region in the v2.1 genome using blastn


**Extraction of the sequences of IPS=-DAPseq**
```bash
samtools faidx Ghirsutum_527_v2.0.fa A02:98974021-98974133 > IPS-peak-sequences.fasta
```

**Blast these sequences**
```bash
makeblastdb -in /data/home/nelsonlab/cotton_vcf/blastn/Ghirsutum_458_v1.0.fa \
    -dbtype nucl \
    -out /data/home/nelsonlab/cotton_vcf/blastn/Ghirsutum_458 \
    -parse_seqids

### Run blastn with identity of 90 as cut-off

blastn -query /data/home/nelsonlab/cotton_vcf/blastn/IPS-peak-sequences.fasta \
    -out /data/home/nelsonlab/cotton_vcf/blastn/IPS-peak \
    -db  /data/home/nelsonlab/cotton_vcf/blastn/Ghirsutum_458 \
    -word_size 25 -perc_identity 90 -dust no -num_threads 20 -max_target_seqs 10 -evalue 1e-10 \
    -outfmt 6 
```
Check the coordinates from blastn
```
A02:98974021-98974133	A02	100.00	113	0	0	1	113	80977441	80977553	1e-52	209
A02:98974021-98974133	A09	98.23	113	2	0	1	113	54083448	54083560	3e-49	198
A02:98974021-98974133	A09	97.35	113	3	0	1	113	1343207	1343319	1e-47	193
A02:98974021-98974133	A09	97.35	113	3	0	1	113	1351452	1351564	1e-47	193
A02:98974021-98974133	A09	97.35	113	3	0	1	113	4435934	4435822	1e-47	193
A02:98974021-98974133	A09	96.46	113	4	0	1	113	4444193	4444081	7e-46	187
A02:98974021-98974133	A09	96.46	113	4	0	1	113	59657210	59657098	7e-46	187
A02:98974021-98974133	A09	96.46	113	4	0	1	113	59665414	59665302	7e-46	187

```

```bash
picard CreateSequenceDictionary -R Ghirsutum_458_v1.0.fa -O Ghirsutum_458_v1.0.dict
```


## Manuscript preparation plot (Misc information)

```r
library(pheatmap)
library(ggplot2)
library(ggpubr)
library(cowplot)

```

```r
my.data <- read.csv("C:/Users/Leon/OneDrive - Cornell University/3_Manuscript/Manuscript8/1_Figures/3_revision3/Phenotype-Ratio.csv", header = T, row.names = 1)
Phenotype <- my.data[,1:9]
metaFeature <- subset(my.data , select = "Group")

annotation_row = metaFeature 

annotation_colors = list(
    Group = c(Foreign="#D01A5B", Mixed="#2a499d", U.S.="#3ab035"))
   

pheatmap(Phenotype,
         color = inferno(21),
         annotation_colors = annotation_colors, 
         annotation_row = annotation_col,
         clustering_distance_r = "correlation", 
         show_colnames =T,show_rownames = F,cluster_rows = T, cluster_cols = T,
         fontsize_row = 4, fontsize_col = 4, angle_col = 45, border_color = NA)

```

```r
GenotypeTest <- read.csv("C:/Users/Leon/OneDrive - Cornell University/Projects/6_Cotton/10_phenotype/SNPs-LintYield-Test.csv", header = T)

A <- ggplot(GenotypeTest[GenotypeTest$Genotype == "AA" | GenotypeTest$Genotype == "GG",], 
       (aes(x=Genotype, y=Yield.WL, fill = Genotype))) +
  geom_boxplot(outlier.size = 0.2, outlier.shape = 21) +
  theme_bw() +
  ylab("Lint yiled under WL condition")+
  scale_fill_manual(values = c(AA="#a38a77", GG="#0ea29a" )) +
  stat_summary(fun = mean, color = "darkred", position = position_dodge(0.75),
               geom = "point", shape = 4, size = 3,
               show.legend = FALSE) +
  stat_compare_means(aes(group = Genotype), method = "t.test", label = "p.format") +
  geom_jitter(color = "black") + ylim(0,120) +
  theme(legend.position = "none")         



B <- ggplot(GenotypeTest[GenotypeTest$Genotype == "AA" | GenotypeTest$Genotype == "GG",], 
       (aes(x=Genotype, y=Yield.WW, fill = Genotype))) +
  geom_boxplot(outlier.size = 0.2, outlier.shape = 21) +
  theme_bw() +
  ylab("Lint yiled under WW condition")+
  scale_fill_manual(values = c(AA="#a38a77", GG="#0ea29a" )) +
  stat_summary(fun = mean, color = "darkred", position = position_dodge(0.75),
               geom = "point", shape = 4, size = 3,
               show.legend = FALSE) +
  stat_compare_means(aes(group = Genotype), method = "t.test", label = "p.format") +
  geom_jitter(color = "black") + ylim(0,120) +
  theme(legend.position = "none")         

plot_grid(A, B)
```

```r
AME_plot <- read.csv("C:/Users/Leon/OneDrive - Cornell University/3_Manuscript/Manuscript8/1_Figures/3_revision3/AMEcount_plot.csv", header = T)

AME_plot$Module <- factor(AME_plot$Module, levels = c("tan", "red", "purple", "yellow", 
                                    "blue", "lightgreen", "darkgreen", "greenyellow", "cyan", "lightcyan", "green","salmon", "midnightblue", "black", "lightyellow", "royalblue", "turquiose", "brown")) 

ggplot(data=AME_plot, aes(x=Module, y=AME)) +
  geom_col(width = 0.5) + 
  theme_classic2()


   
```

### Plot the CV over differnece modules
```r
### read the CV files
CV_plot <- read.csv("C:/Users/Leon/OneDrive - Cornell University/3_Manuscript/Manuscript8/1_Figures/3_revision3/Gh_WGCNA-TPM-CV-plot.csv", header = T)
head(CV_plot)


### Relevel the order of boxplot


random_plot <- CV_plot[sample(nrow(CV_plot), 193), ]
random_plot$Module <- "random"

Plot_final <- rbind(random_plot, CV_plot)
Plot_final$Module <- factor(Plot_final$Module, levels = c("tan", "red", "purple", "yellow", 
                                    "blue", "lightgreen", "darkgreen", "greenyellow", "cyan", "lightcyan",
                                    "green","salmon", "midnightblue", "black", "lightyellow", "royalblue", 
                                    "turquoise", "brown", "random")) 

### Plot the boxplot
ggplot(Plot_final[Plot_final$Module != "NA",], 
       (aes(x=Module, y=CV))) +
  geom_boxplot(fill = "#0ea29a", alpha = 0.5, outlier.colour = "black", outlier.size = 0.6) +
  theme_classic() +
  ylab("CV across accession over two conditions")+
  stat_summary(fun = mean, color = "darkred", position = position_dodge(0.75),
               geom = "point", shape = 4, size = 3,
               show.legend = FALSE) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45
                                                             , vjust = 0.5, hjust=1)) 

```

### Plot the line graph of non-lable probe shift comparison
```r
EMSA <- read.csv("/Users/Leon/OneDrive - Cornell University/3_Manuscript/Manuscript8/1_Figures/EMSA-probe-complexSignal.csv", header = T)

EMSA$Index <- paste0(EMSA$SNP, "_", EMSA$Concentration)
EMSA
```

```r

EMSA$Concentration <- factor(EMSA$Concentration, levels = c("1X", "2X", "5X", "10X", "20X")) 

EMSA  %>% ggplot(aes(x=Concentration, y = Ratio, color= SNP, group = SNP)) +
  geom_line(alpha = 0.8,size = 0.8)+ 
  geom_point(alpha = 0.8,size = 3)+
  theme_classic()+
  scale_color_manual(values = c("G"="#1B9E77", "C"="#D95F02", "T"="midnightblue")) +
  ylab("(Protein + probe) / probe ratio") +
  xlab("Concentration of non-lable probe")
```


