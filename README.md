
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
    - [Intersect Variants data with peak profiles](#intersect-variants-data-with-peak-profiles)
    - [Population genomics: SNP analysis of IPS flanking SNPs](#population-genomics-snp-analysis-of-ips-flanking-snps)
      - [BYU panel with deep sequencing data](#byu-panel-with-deep-sequencing-data)
      - [TAMU panel with deep sequencing data](#tamu-panel-with-deep-sequencing-data)
      - [Using the published larger panel for ana](#using-the-published-larger-panel-for-analysis)


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


### Population genomics: SNP analysis of IPS flanking SNPs
#### BYU panel with deep sequencing data

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

#### TAMU panel with deep sequencing data

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

#### Using the published larger panel for analysis
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



