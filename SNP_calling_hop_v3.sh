##This document details commands used for the SNP calling for two hop genoes against the Cascade reference genome 

## ----------------------0.Download reference genomes---------------------------------------------------------------------------------------

## Latest 2021 Cascade assembly from Hopbase

##This assembly was generated from PacBio reads using Falcon

mkdir -p projects/hop/raw_data/cascade_2021
wget http://hopbase.org/content/maskedDedupCascade/version5/maskedCascadePrimary.fasta.gz
# unzip data:
gunzip *

 
### Pilgrim and 316/1/10 raw sequences were moved to project directory

## Pilgrim (S1) and 316/1/10 (S3) P150x2 sequences were generaged on illumina novaseq6000 platform 



##------------------------------------------1. Estimate genome coverage of S1 and S3----------------------------------------------------------------------------------

## ~20x coverage of this size genomes is sufficient for SNP calling.

##Change file permissions from archive to executble
chmod +x /projects/hop/raw_data/trim/S*/*/*.fq.gz


for Strain in S1 S3; do
 for DataDir in $(ls -d ../../projects/hop/raw_data/trim/$Strain); do
    F_Read=$(ls $DataDir/F/*fq.gz)
    R_Read=$(ls $DataDir/R/*fq.gz)
    Outdir=/projects/hop/raw_data/trim/$Strain/genome_coverage_2021
    echo $F_Read
    echo $R_Read
    echo $Outdir
    mkdir -p $Outdir
    ProgDir=/home/hajduk/git_repos/scripts/hop_genomics/dna_qc
    sbatch $ProgDir/count_nucl.sh $F_Read $R_Read 3000 $Outdir #Estimated genome size
 done
done



##------------------------------------------2. Carryout alignment of genomes (bowtie2)-----------------------------------------------


## Alignemt of S1 S3 reads vs Cascade genome

##Alignment of reads from a single run using Bowtie:

for Hop_genome in S1 S3; do
  Reference=$(ls ../../projects/hop/cascade_2021/cascadePrimary.fasta)
  for StrainPath in $(ls -d ../../projects/hop/raw_data/trim/"$Hop_genome"); do
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    F_Read=$(ls $StrainPath/F/*_trim.fq.gz)
    R_Read=$(ls $StrainPath/R/*_trim.fq.gz)
    echo "$Organism - $Strain"
    echo $F_Read
    echo $R_Read
    OutDir=$StrainPath/alignment_2021
    mkdir -p $OutDir
    ProgDir=/home/hajduk/git_repos/scripts/hop_genomics/SNP_calling
    sbatch $ProgDir/bowtie.sh $Reference $F_Read $R_Read $OutDir
  done
done


##------------------------------------------3. Rename input mapping files in each folder by prefixing with the strain ID------------------------------------ 
 
for Strain in S1 S3; do
  for filename in /projects/hop/raw_data/trim/$Strain/alignment_2021; do #filename is the pathway specified after "in"
   echo $Strain
     mv "$filename/cascadePrimary.fasta_aligned.sam" "$filename/"$Strain"_unmasked.fa_aligned.sam"
     mv "$filename/cascadePrimary.fasta_aligned.bam" "$filename/"$Strain"_unmasked.fa_aligned.bam"
     mv "$filename/cascadePrimary.fasta_aligned_sorted.bam" "$filename/"$Strain"_unmasked.fa_aligned_sorted.bam"
     mv "$filename/cascadePrimary.fasta.indexed.1.bt2" "$filename/"$Strain"_contigs_unmasked.fa.indexed.1.bt2"
     mv "$filename/cascadePrimary.fasta.indexed.2.bt2" "$filename/"$Strain"_contigs_unmasked.fa.indexed.2.bt2"
     mv "$filename/cascadePrimary.fasta.indexed.3.bt2" "$filename/"$Strain"_contigs_unmasked.fa.indexed.3.bt2"
     mv "$filename/cascadePrimary.fasta.indexed.4.bt2" "$filename/"$Strain"_contigs_unmasked.fa.indexed.4.bt2"
     mv "$filename/cascadePrimary.fasta.indexed.rev.1.bt2" "$filename/"$Strain"_contigs_unmasked.fa.indexed.rev.1.bt2"
     mv "$filename/cascadePrimary.fasta.indexed.rev.2.bt2" "$filename/"$Strain"_contigs_unmasked.fa.indexed.rev.2.bt2"
     mv "$filename/cascadePrimary.fasta_RPK.txt" "$filename/"$Strain"_contigs_unmasked.fa_RPK.txt"
     mv "$filename/cascadePrimary.fasta_aligned_sorted.bam.index" "$filename/"$Strain"_unmasked.fa_aligned_sorted.bam.index"
  done 
done 


##-----------------------------------------4. Remove multimapping reads, discordant reads, optical duplicates (samtools)--------------------------------------------------- 
##and add read group and sample name to each mapped read. 

  for Strain in S1 S3; do 
    for input in ../../projects/hop/raw_data/trim/$Strain/alignment_2021/"$Strain"_unmasked.fa_aligned.sam; do
  #   Jobs=$(squeue -u ${USER} --noheader --array | wc -l)
  #      while [ $Jobs -gt 1 ]; do
  #          sleep 25m
  #          printf "."
  #          Jobs=$(squeue -u ${USER} --noheader --array | wc -l)
  #      done
  #  printf "\n"
    OutDir=../../projects/hop/archive_raw_data/trim/$Strain/alignment_2021/nomulti
  #  mkdir -p $OutDir  
    ProgDir=/home/hajduk/git_repos/scripts/hop_genomics/SNP_calling
    sbatch $ProgDir/pre_SNP_calling.slurm.sh $input $Strain $OutDir # This will add read group and sample name to each mapped read. Preferably, use the shortest ID possible.
    done
  done
 


##this is just copying the reference to all of the query genomes' nomulti folders 
reference=/projects/hop/cascade_2021
for Strain in S1 S3; do
  for filename in /projects/hop/raw_data/trim/$Strain/alignment_2021; do
   echo $Strain
    cp "$reference"/cascadePrimary.fasta "$filename"/nomulti
  done 
done 


##----------------------------------------5. Run RealignerTargetCreator (gatk)---------------------------------------------------------
##to get realigner.intervals used in the next step
## not part of GATK4
Reference=../../projects/hop/cascade_2021/cascadePrimary.fasta  
for Strain in S1 S3; do
  for input in ../../projects/hop/raw_data/trim/$Strain/alignment_2021/nomulti/"$Strain"_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam; do
    echo $Strain
    Outdir=/projects/hop/raw_data/trim/$Strain/alignment/nomulti/pre_indel_realignment
    mkdir -p $Outdir
    ProgDir=/home/hajduk/git_repos/scripts/hop_genomics/SNP_calling
    sbatch $ProgDir/pre_indel_realignment_v2.sh $Reference $input $Strain $Outdir
  done 
done     


##----------------------------------------6. Run IndelRealigner using realigner.intervals file (gatk)----------------------------------------------
##from previous step
##indelrealigner is not part of the GATK 4 pipeline. Realignment based haplotype caller does this step now. 

Reference=../../projects/hop/cascade_2021/cascadePrimary.fasta  
for Strain in S1 S3; do
  for input in ../../projects/hop/raw_data/trim/$Strain/alignment_2021/nomulti/"$Strain"_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam; do
    target_intervals=../../projects/hop/raw_data/trim/$Strain/alignment_2021/nomulti/pre_indel_realignment/realigner.intervals
    echo $Strain
    Outdir=/projects/hop/raw_data/trim/$Strain/alignment_2021/nomulti/corrected_bam
    mkdir -p $Outdir
    ProgDir=/home/hajduk/git_repos/scripts/hop_genomics/SNP_calling
    sbatch $ProgDir/indel_realignment_v2.sh $Reference $Strain $input $target_intervals $Outdir
  done 
done  

## The next step in the GATK SNP calling pileup pipeline is base quality score recalibration (BQSR) which uses known variant sets to mask out positions where variations are commonly found. 
## If no preexisting SNP and INDEL data is available, variant calling may be performed without applying BQSR, by applying hard filtering on the raw call data set 


##----------------------------------------7. Submit bam files from step 6. for SNP calling with haplotype caller in gvcf module (gatk) --------------------------------

Reference=../../projects/hop/cascade_2021/cascadePrimary.fasta  
for Strain in S1 S3; do
  for input in ../../projects/hop/raw_data/trim/$Strain/alignment_2021/nomulti/corrected_bam/*_realigned.bam; do
    echo $Strain
    Outdir=/projects/hop/SNP_calling/haplotype_calls_v2/$Strain
    mkdir -p $Outdir
    ProgDir=/home/hajduk/git_repos/scripts/hop_genomics/SNP_calling
    sbatch $ProgDir/haplotyp_caller_v2.sh $Reference $Strain $input $Outdir
  done 
done  

##copy gvcf from individual haplotype calls to shared folder
##reference genome with .dict and .fa should also be copied over

##-------------------------------------------8. CombineGVCF and GenotypeGVCF (gatk)----------------------------------------------------------------------------------
# submit g.vcf files generated from haplotype caller step
## filter for only biallelic SNPs. 

ProgDir=/home/hajduk/git_repos/scripts/hop_genomics/SNP_calling
sbatch $ProgDir/genotype_gvcf_v2.sh 



#--------------------------------------------9. Hard filter variants (gatk)---------------------------------------------------------------------------------------- 

#with vcftools

#for strain in S1 S3; do 
# for vcf in /projects/hop/SNP_calling/haplotype_calls_v2/$strain/*_genotype.vcf; do 
  # mq=40
  # qual=30
  # dp=10
  # gq=30
  # na=0.95
  # removeindel=Y/N
#  Outdir=/projects/hop/SNP_calling/haplotype_calls_v2/$strain/filtered_vcftools
#  mkdir -p $Outdir #does not create the outdir directories for some reason
#  ProgDir=/home/hajduk/git_repos/scripts/hop_genomics/SNP_calling
#  sbatch $ProgDir/sub_vcf_parser.sh $vcf 40 30 10 30 1 N $Outdir $strain
# done 
#done 


#Hard filtering with GATK tools. Recommended when a variant callset is too small for VQSR or for which truth/training sets are not available.
#https://sites.google.com/a/broadinstitute.org/legacy-gatk-forum-discussions/tutorials/2806-how-to-apply-hard-filters-to-a-call-set

Reference=../../projects/hop/cascade_2021/cascadePrimary.fasta
for strain in cohort_v2; do 
 for vcf in /projects/hop/SNP_calling/haplotype_calls_v2/$strain/*_genotype.vcf; do 
  Outdir=/projects/hop/SNP_calling/haplotype_calls_v2/$strain/filtered_gatk
  mkdir -p $Outdir
  ProgDir=/home/hajduk/git_repos/scripts/hop_genomics/SNP_calling
  sbatch $ProgDir/variant_hard_filter3.sh $Reference $strain $vcf $Outdir
 done 
done 


#-------------------------------------------10. Make table with filtered variants (gatk)---------------------------------------------------------------------------------

##creates tsvs that show filtered/unfiltered variants and reason for filtering. this was used for plotting filtering results.
##creates tsv with snp types and genotype field

Reference=../../projects/hop/cascade_2021/cascadePrimary.fasta
for strain in hop_joint; do 
 for vcf in /projects/hop/SNP_calling/haplotype_calls_v2/$strain/filtered_gatk/*_filtered_*; do 
  Outdir=/projects/hop/SNP_calling/haplotype_calls_v2/$strain/filtered_gatk
  mkdir -p $Outdir
  ProgDir=/home/hajduk/git_repos/scripts/hop_genomics/SNP_calling
  sbatch $ProgDir/variants_to_table.sh $Reference $strain $vcf $Outdir
 done 
done 
 

##-------------------------------------------------11. Plot filtering results (R)------------------------------------------------------------------------

##setwd("C:/Users/hajdk/OneDrive/New folder/Thesis/Chapter4_Hop_genome_analysis")
##library('ggplot2')
##library('tools')
##library('tidyverse')
##library('ggpubr')
##
##
##df1<-read.table("cohort_snps_stats.tsv" , sep = '\t', header = TRUE)
##
##str(df1)
##
##df1$PASS <- df1$FILTER
##df1$PASS[df1$PASS!="PASS"] <- "FAIL"
##df1$PASS <- as.factor(df1$PASS)
#### Histogram with density plot
##geom_histogram(aes(y=..density..), alpha=0.2, bins=50)
##
##p0 <- ggplot(df1, aes(x=QD)) +
##  geom_density(alpha=.2) +
##  theme(legend.position="none")
##print(p0)
##ggsave("../QD_all.pdf", p0)
##
##p1 <- ggplot(df1, aes(x=QD, colour=PASS, fill=PASS)) +
##  geom_density(alpha=.2) +
##  theme(legend.position="right") 
##print(p1)
##
##p2 <- ggplot(df1, aes(x=FS, colour=PASS, fill=PASS)) +
##  geom_density(alpha=.2) +
##  theme(legend.position="") +
##  scale_x_continuous(trans='log10')
##print(p2)
##
##p3 <- ggplot(df1, aes(x=SOR, colour=PASS, fill=PASS)) +
##  geom_density(alpha=.2) +
##  theme(legend.position="none")
##print(p3)
##
##p4 <- ggplot(df1, aes(x=MQ, colour=PASS, fill=PASS)) +
##  geom_density(alpha=.2) +
##  theme(legend.position="none")
##print(p4)
##
##p5 <- ggplot(df1, aes(x=MQRankSum, colour=PASS, fill=PASS)) +
##  geom_density(alpha=.2) +
##  theme(legend.position="none")
##print(p5)
##
##p6 <- ggplot(df1, aes(x=ReadPosRankSum, colour=PASS, fill=PASS)) +
##  geom_density(alpha=.2) +
##  theme(legend.position="none")
##print(p6)
##
##
##df2<-read.table("cohort_indels_stats.tsv" , sep = '\t', header = TRUE)
##
##str(df1)
##
##df2$PASS <- df2$FILTER
##df2$PASS[df2$PASS!="PASS"] <- "FAIL"
##df2$PASS <- as.factor(df1$PASS)
#### Histogram with density plot
##geom_histogram(aes(y=..density..), alpha=0.4, bins=50)
##
##p7 <- ggplot(df2, aes(x=QD, colour=PASS, fill=PASS)) +
## geom_density(alpha=.4) +
##  theme(legend.position="none")
##print(p7)
##
##p8 <- ggplot(df2, aes(x=FS, colour=PASS, fill=PASS)) +
##  geom_density(alpha=.2) +
##  theme(legend.position="none") +
##  scale_x_continuous(trans='log10')
##print(p8)
##
##p9 <- ggplot(df2, aes(x=ReadPosRankSum, colour=PASS, fill=PASS)) +
##  geom_density(alpha=.2) +
##  theme(legend.position="none")
##print(p9)


#alternatively use bash command to execute R script for visualising filtering results (this didn't work for me and just ran the script above directly from R)
#```bash
#srun --partition long --time 0-04:00:00 --mem 10gb --cpus-per-task 10 --pty bash
#SnpTsv=$(ls projects/hop/SNP_calling/haplotype_calls/cohort/filtered_gatk5/*_snps_stats.tsv)
#Outdir=/projects/hop/SNP_calling/haplotype_calls/cohort/filtering_plots_R
#mkdir -p $Outdir
#ProgDir=/home/hajduk/git_repos/scripts/hop_genomics/SNP_calling
#Rscript --vanilla $ProgDir/plot_filter_results_KH.R $SnpTsv
#```
#
#```bash
#srun --partition long --time 0-04:00:00 --mem 10gb --cpus-per-task 10 --pty bash
#IndelTsv=$(ls projects/hop/SNP_calling/haplotype_calls/cohort/filtered_gatk5/*_indels_stats.tsv)
#Outdir=/projects/hop/SNP_calling/haplotype_calls/cohort/filtering_plots_R
#mkdir -p $Outdir
#ProgDir=/home/hajduk/git_repos/scripts/hop_genomics/SNP_calling
#Rscript --vanilla $ProgDir/plot_filter_results_KH.R $IndelTsv
#``` 


##--------------------------------------------12. Annotate filtered VCF file for SNP effects (snpeff)---------------------------------------------------------------------------------- 

## Create custom SnpEff genome database

```bash
SnpEff=/projects/oldhome/sobczm/bin/snpEff
nano $SnpEff/snpEff.config
```

##download and unzip genemodels from hopbase.org
## /projects/hop mkdir gene_predictions
## $gene_predictions wget http://hopbase.org/content/maskedDedupCascade/version5/maskedCascadePrimary.genes.gff.gz
## $gene_predictions gunzip maskedCascadePrimary.genes.gff.gz

## Add the following lines to the section with databases:

```
##---
## EMR Databases
##----
## H. lupulus cascade genome
cascadePrimary.genome : hop_ref2021
```

##Collect input files

```bash
Reference=$(ls /projects/hop/cascade_2021/cascadePrimary.fasta)
Gff=$(ls /projects/hop/gene_prediction/cascade/maskedCascadePrimary.genes.gff)
SnpEff=/projects/oldhome/sobczm/bin/snpEff
mkdir $SnpEff/data/cascadePrimary
cp $Reference $SnpEff/data/cascadePrimary/sequences.fa
cp $Gff $SnpEff/data/cascadePrimary/genes.gff

##Build database using GFF3 annotation
java -jar $SnpEff/snpEff.jar build -gff3 -v cascadePrimary
```

## Annotate VCF files
CurDir=../../projects/hop/SNP_calling
cd $CurDir
for a in $(ls /projects/hop/SNP_calling/haplotype_calls_v2/hop_joint/filtered_gatk/*_filtered_snps.vcf); do  
    echo $a
    filename=$(basename "$a")
    Prefix=${filename%.vcf}
    OutDir=$(ls -d /projects/hop/SNP_calling/haplotype_calls_v2/hop_joint)
    #mkdir -p $Outdir
    SnpEff=/projects/oldhome/sobczm/bin/snpEff
    java -Xmx4g -jar $SnpEff/snpEff.jar -v -ud 0 cascadePrimary $a > $OutDir/"$Prefix"_annotated.vcf
    mv snpEff_genes.txt $OutDir/snpEff_genes_"$Prefix".txt
    mv snpEff_summary.html $OutDir/snpEff_summary_"$Prefix".html
done

  #submitted scripts below from the directory where files are located: /projects/hop/SNP_calling/haplotype_calls/cohort 
  ```bash 
    SnpEff=/projects/oldhome/sobczm/bin/snpEff
    #mv cohort_filtered_snps* $OutDir/.
    #-
    #Create subsamples of SNPs containing those in a given category
    #-
    #genic (includes 5', 3' UTRs)
    java -jar $SnpEff/SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') || (ANN[*].EFFECT has 'nonsense_variant') || (ANN[*].EFFECT has 'synonymous_variant') || (ANN[*].EFFECT has 'intron_variant') || (ANN[*].EFFECT has '5_prime_UTR_variant') || (ANN[*].EFFECT has '3_prime_UTR_variant')" hop_joint_filtered_snps_annotated.vcf > hop_joint_filtered_snps_genic.vcf
    #coding
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'synonymous_variant')" hop_joint_filtered_snps_annotated.vcf > hop_joint_filtered_snps_coding.vcf
    #non-synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant')" hop_joint_filtered_snps_annotated.vcf > hop_joint_filtered_snps_nonsyn.vcf
    #synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'synonymous_variant')" hop_joint_filtered_snps_annotated.vcf > hop_joint_filtered_snps_syn.vcf
    #Four-fold degenrate sites (output file suffix: 4fd)
   ```
  #```bash 
  #  ProgDir=/projects/oldhome/sobczm/bin/popgen/summary_stats
  #  Outdir=/projects/hop/SNP_calling/haplotype_calls/cohort/
  #  python $ProgDir/parse_snpeff_synonymous.py $OutDir/cohort_filtered_snps_syn.vcf
  #  AllSnps=$(cat $Outdir/cohort_filtered_snps_annotated.vcf | grep -v '#' | wc -l)
  #  GeneSnps=$(cat $Outdir/cohort_filtered_snps_gene.vcf | grep -v '#' | wc -l)
  #  CdsSnps=$(cat $Outdir/cohort_filtered_snps_coding.vcf | grep -v '#' | wc -
  #  NonsynSnps=$(cat $Outdir/cohort_filtered_snps_nonsyn.vcf | grep -v '#' | wc -l  )
  #  SynSnps=$(cat $Outdir/cohort_filtered_snps_syn.vcf | grep -v '#' | wc -l)
  #  printf "Comparison\$AllSnps\tGeneSnps\tCdsSnps\tSynSnps\tNonsynSnps\n"
  #  printf "cohort_filtered_snps\t$AllSnps\t$GeneSnps\t$CdsSnps\t$SynSnps\t$NonsynSnps\n"
  #```
grep -v '#' | wc -l hop_joint_filtered* 
  ### Results:
  ### 8850932 cohort_filtered_snps_annotated.vcf
  ### 178284 cohort_filtered_snps_coding.vcf
  ###  585881 cohort_filtered_snps_genic.vcf
  ###  115559 cohort_filtered_snps_nonsyn.vcf
  ###   71507 cohort_filtered_snps_syn.vcf

#manipulate snpeff vcf:
more cohort_annotated_snps_extracted.vcf | grep -v '#' | grep "PASS" | cut -f 1,2,8,9,10 | sed 's/|/\t/g' | cut -f 1,2,4,5,7 | grep -v "intergenic\|stream" > cohort_filtered_snps_annotated_table > 
# grep -v "" takes every line that start with " 
# wc -l counts  
# cut -f will take whole columns from a file , numbers to define which ones
# cut -f 1 | sort | uniq | wc -l  take a column and removes duplicates
# sed 's/|/\t/g' will put things in separate columns that had a pipe (|) between them. after they are separated, information can be extracted
# https://www.youtube.com/watch?v=-rmreyRAbkE very good snpeff interpretation and manipulation video
#want to extract information on where the snp is , impact level, what gene is in 
# venny (bioinfogp.cnb.csic.es) makes vennn diagrams for comparing snps in different genomes

##Pipeline complete 


grep -v '#' | wc -l filtered_cohort_onlyparent_gt.tsv   1877794 


