
# **Transcriptome and gene regulatory analysis of cotton omics project**

**Author: Li'ang Yu, Andrew Nelson**\
**Date: June 14th, 2022**

<img width="800" alt="image" src="https://user-images.githubusercontent.com/69836931/214942202-6e3ce9af-fea6-43e3-ae6c-6d43a4665c4d.png">


- [**Transcriptome and gene regulatory analysis of cotton omics project**](#transcriptome-and-gene-regulatory-analysis-of-cotton-omics-project)
  - [STEP 1 Quantification of reads counts using featureCounts](#step-1-quantification-of-reads-counts-using-featurecounts)
    - [Perform the featureCounts for coding genes at gene level](#perform-the-featurecounts-for-coding-genes-at-gene-level)
    - [check if meta-feature and features counted correctly for transcripts](#check-if-meta-feature-and-features-counted-correctly-for-transcripts)
  - [STEP 2 Normalization of count numbers](#step-2-normalization-of-count-numbers)
    - [Normalization using TPM](#normalization-using-tpm)
    - [Normalization using DESEq2](#normalization-using-deseq2)
  - [output: html\_document](#output-html_document)
    - [Import raw count data](#import-raw-count-data)
    - [build matrix for each sample based on tissues, accessions, and replicates](#build-matrix-for-each-sample-based-on-tissues-accessions-and-replicates)
    - [normalization and DEGs with experiemtnal design been specified](#normalization-and-degs-with-experiemtnal-design-been-specified)
    - [Count Normalization](#count-normalization)
  - [Variants calling of 22 accessions using WGS data](#variants-calling-of-22-accessions-using-wgs-data)
    - [Trim the reads](#trim-the-reads)
    - [Variants calling of 22 accessions using (DNA-seq data)](#variants-calling-of-22-accessions-using-dna-seq-data)
    - [Filter the variants of multiple samples](#filter-the-variants-of-multiple-samples)
  - [Correlation of trait genes TPM and phenotype](#correlation-of-trait-genes-tpm-and-phenotype)
  - [**ATAC-seq chromatin accesibility profiling**](#atac-seq-chromatin-accesibility-profiling)
    - [Trim reads and quality control of reads](#trim-reads-and-quality-control-of-reads)


## STEP 1 Quantification of reads counts using featureCounts 
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
## STEP 2 Normalization of count numbers
### Normalization using TPM


### Normalization using DESEq2
Using R to to Normalization based on libraries
Using R to identify DEGs across pair-wise samples
Using vsd-tranformation to prepare the input for WGCNA

---
title: "Cotton_comp_DESeq2""
author: "A, Nelson; L, Yu"
date: "2021-3-17"
output: html_document
---
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

## Variants calling of 22 accessions using WGS data
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
```bash
seq_dir="/mnt/Knives/1_data/9_Cotton_omics/5_WGS"

cat $seq_dir/21552Dpa_SA-1019_S227_L00?_R1_001.fastq.gz > $seq_dir/SA_1019_R1.fastq.gz
cat $seq_dir/21552Dpa_SA-1019_S227_L00?_R2_001.fastq.gz > $seq_dir/SA_1019_R2.fastq.gz

cat $seq_dir/21552Dpa_SA-0159_S245_L00?_R1_001.fastq.gz > $seq_dir/SA_0159_R1.fastq.gz
cat $seq_dir/21552Dpa_SA-0159_S245_L00?_R2_001.fastq.gz > $seq_dir/SA_0159_R2.fastq.gz

cat 21552Dpa_SA-1212_S222_L00?_R1_001.fastq.gz > SA_1212_R1.fastq.gz
cat 21552Dpa_SA-1212_S222_L00?_R2_001.fastq.gz > SA_1212_R2.fastq.gz
```

### Trim the reads
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
#mv /mnt/Knives/1_data/Sorghum_Pangenome/1_Reads/*.fastq.gz /mnt/Kanta/Sorghum_reads/
```

### Variants calling of 22 accessions using (DNA-seq data)
Note from 20220726: The SRR6311717 sample was interrupted and needs to be re-launched
Sample list of all accesions
```
SA_0159
SA_1019
SA_1212
SRR6311500
SRR6311566
SRR6311570
SRR6311571
SRR6311717
SRR6311732
SRR6311780
SRR6311856
SRR6311857
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

  ## Map reads by bwa piepline
  bwa mem -M -R "@RG\tID:$file\tSM:$file\tPL:ILLUMINA" \
    -t $nt  \
    $Ref_dir/Ghirsutum_527_v2.0.fa \
    ${BAM_dir}/${file}_forward_paired.fq.gz \
    ${BAM_dir}/${file}_reverse_paired.fq.gz | samtools view -@10 -bS - -o ${workdir}/${file}/${file}.reorder.bam 

  samtools view -b -F 4 ${workdir}/${file}/${file}.reorder.bam > ${workdir}/${file}/${file}.cut.bam

  ## Sort bam files
  picard SortSam \
   INPUT=${workdir}/${file}/${file}.cut.bam \
   OUTPUT=${workdir}/${file}/${file}.sort.bam SORT_ORDER=coordinate

    #Mark and remove duplictaes from Bam
   picard MarkDuplicates -REMOVE_DUPLICATES true \
    -I ${workdir}/${file}/${file}.sort.bam \
    -O ${workdir}/${file}/${file}.sort.dup.bam \
    -M ${workdir}/${file}/${file}.metrics.txt 

 #add group numbers and infotrmation
  picard AddOrReplaceReadGroups \
    -I ${workdir}/${file}/${file}.sort.dup.bam \
    -O ${workdir}/${file}/${file}.sort.dup.Gr.bam -LB cotton -PL illumina -PU Diversity -SM ${file} &

   #index bam file
   samtools index ${workdir}/${file}/${file}.sort.dup.Gr.bam

  java -Xmx128G -jar $GATK_dir/GenomeAnalysisTK.jar \
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
java -Xmx128G -jar $GATK_dir/GenomeAnalysisTK.jar \
 -T CombineGVCFs -nct $nt \
 -R $Ref_dir/Ghirsutum_527_v2.0.fa \
 --variant ${workdir}/gvcf.list 
 -o cotton.cohort.gvcf \

### Transform the VCFs format
java -Xmx128G -jar $GATK_dir/GenomeAnalysisTK.jar \
 -T GenotypeGVCFs -nct $nt \
 -R $Ref_dir/Ghirsutum_527_v2.0.fa \
 --variant $workdir/cohort.g.vcf -o cotton_DNA_genotype.vcf
```
### Filter the variants of multiple samples
```bash
#! bin/bash
workdir="/home/liangyu/1_Project/4_spinach/5_VCFs" 

### STEP.1 remove loci with missing data based on cut-off to reduce bias of total mapping depth
#keep loci with up to 10% missing data for plot
vcftools --vcf $workdir/F1_map_variants.vcf \
	--max-missing 0.90 --min-alleles 2 \
	--recode --recode-INFO-all \
	--out $workdir/F1_plot.vcf


### STEP.2 PLOT DP, QD, and QUAL for each loci 
###Extrat VCF depth information (This is the combined depth among samples)
egrep -v "^#" $workdir/F1_plot.vcf.recode.vcf | \
cut -f 8 | \
sed 's/^.*;DP=\([0-9]*\);.*$/\1/' > $workdir/F1_map_DP.txt


###Extrat VCF QualitybyDepth information 
egrep -v "^#" $workdir/F1_plot.vcf.recode.vcf | \
cut -f 8 | \
sed 's/^.*;QD=\([0-9]*.[0-9]*\);.*$/\1/' > $workdir/F1_map_QD.txt

###Extrat VCF Quality information
egrep -v "^#" $workdir/F1_plot.vcf.recode.vcf | \
cut -f 6 > $workdir/F1_map_QUAL.txt

#Plot three component
Rscript 1_Plot_depth_quality.R --DP $workdir/F1_map_DP.txt --QUAL $workdir/F1_map_QUAL.txt --QD $workdir/F1_map_QD.txt


### STEP.3 Remove loci based on cutoff derived from plot information 
## MinDP and MaxDP is determined by average total depth (total depth/population size)
## minGQ=20: 1% chance that the call is incorrect
## minQ=20: 1 % chance that there is no variant at the site

vcftools --vcf $workdir/F1_plot.vcf.recode.vcf  \
	--minDP 5 --maxDP 25 \
	--minGQ 20 --minQ 20 \
	--recode --recode-INFO-all \
	--out $workdir/F1.vcf

#filter based qualitybydepth (QD) (QD is the normalized quality based on mapping depth usually set as 2)
vcffilter -f "QD > 5" $workdir/F1.vcf.recode.vcf > $workdir/F1_final_filter.vcf
## Line NO. 1079049


### STEP 4. extract genotype information only 
cat $workdir/F1_final_filter.vcf | perl while (<>) { s!:\S+!!g; print; } > F1_final_filter.clean 
```








## Correlation of trait genes TPM and phenotype
```r
### Plot WW AND WL together
datTrait1 <- read.csv("C:/Users/Leon/OneDrive - Cornell University/Projects/6_Cotton/10_phenotype/FIBER_yield_and_quality_2019_BLUEs_TIPO_removed.csv", header = T)
rownames(datTrait1) <- datTrait1$genotype

datTrait2 <- read.csv("C:/Users/Leon/OneDrive - Cornell University/Projects/6_Cotton/Results/salmon/Metabolites_2019_BLUEs.csv", header = T, row.names = 1)
datTrait2$genotype <- rownames(datTrait2) 

TPM_bind <- read.csv("C:/Users/Leon/OneDrive - Cornell University/Projects/6_Cotton/3_FeatureCounts/Gh_WGCNA-TPM-ave_bind.csv", header = T, row.names = 1)
TPM_corr = as.data.frame(t(TPM_bind))
```

**select gene of interests**
```r
### import lint related genes
yield_gene <- c("Gohir.D13G022300", 
                "Gohir.A02G132300", 
                "Gohir.A08G064100", 
                "Gohir.D08G072600")

Cor1 <- TPM_corr[,yield_gene]
Cor1$genotype <- row.names(Cor1)
Cor2 = datTrait1[,c("genotype", "Lint_yield", "treatment")]

### combine two sets
Comp = left_join(Cor1, Cor2, by = "genotype")
write.csv(Comp, "C:/Users/Leon/OneDrive - Cornell University/Projects/6_Cotton/14_MapMan/lint_expression.csv")
```

```r
### import sucrose related genes
sucrose_gene <- c("Gohir.D13G022300", 
                  "Gohir.A02G132300", 
                  "Gohir.A08G064100", 
                  "Gohir.D08G072600")

Cor1 <- TPM_corr[,X2019_GCMS_15]
Cor1$genotype <- row.names(Cor1)
Cor2 = datTrait1[,c("genotype", "X2019_GCMS_15", "treatment")]

### combine two sets
Comp = left_join(Cor1, Cor2, by = "genotype")
write.csv(Comp, "C:/Users/Leon/OneDrive - Cornell University/Projects/6_Cotton/14_MapMan/sugar_expression.csv")
```
**make the scatter plot**
```r
Cor1 <- data.frame(cbind(TPM_corr[,"Gohir.D07G180900"],row.names(TPM_corr)))
colnames(Cor1) = c("TPM","genotype")
Cor2 = datTrait2[,c("genotype", "X2019_GCMS_15", "treatment")]

Comp = left_join(Cor1, Cor2, by = "genotype")
Comp$TPM = as.numeric(Comp$TPM)
Comp$X2019_GCMS_15 = as.numeric(Comp$X2019_GCMS_15)

###Plot only the WL condition 
ggscatter(na.omit(Comp[Comp$treatment == "WL",]), x = "TPM", y = "X2019_GCMS_15", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "TPM", ylab = "Glucose content",
          color = "black", shape = 21, size = 3, 
          add.params = list(color = "blue", fill = "lightgray"))

###Plot under the two conditions
ggscatter(na.omit(Comp), x = "TPM", y = "X2019_GCMS_15", 
          add = "reg.line", 
          color = "treatment",
          fullrange = TRUE,
          rug = TRUE,
          cor.method = "spearman",
          shape = "treatment",
          xlab = "TPM", ylab = "X2019_GCMS_15") + 
            stat_cor(aes(color = treatment), label.x = 3)
```

## **ATAC-seq chromatin accessibility profiling**
### Trim reads and quality control of reads
### Reads mapping and alignment

## **DAP-seq data of selected hub transcription factors**
