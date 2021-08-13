#!/usr/bin/env bash
#SBATCH -J variants_to_table
#SBATCH --partition=long
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=30


##########################################################################

#Extract fields from a VCF file to a tab-delimited table. This tool extracts specified fields for each variant in a VCF file to a tab-delimited table, 
#which may be easier to work with than a VCF. By default, the tool only extracts PASS or . (unfiltered) variants in the VCF file. 
#Filtered variants may be included in the output by adding the --show-filtered flag. The tool can extract both INFO (i.e. site-level) 
#fields and FORMAT (i.e. sample-level) fields.

#INPUT:
# 1st argument: Refrence fasta 
# 2nd argument: sample name (prefix) to be used to identify it in the future
# 3nd argument: Input hard-filtered vcf files with SNPs and INDELs 
#OUTPUT:
# tabular format of 


reference=$1
strain=$2
vcf=$3
outdir=$4


WorkDir=/projects/hop/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir

cp $reference $WorkDir
cp $vcf $WorkDir
cd $WorkDir


samtools faidx cascadePrimary.fasta

picard=/home/connellj/miniconda2/share/picard-2.18.29-0/picard.jar
java -jar $picard CreateSequenceDictionary \
     R=cascadePrimary.fasta \
     O=cascadePrimary.dict 


gatk VariantsToTable \
     -R cascadePrimary.fasta \
     -V ${strain}_filtered_snps.vcf \
     -F CHROM -F POS -F TYPE -F FILTER \
     -F QD -F FS -F SOR -F MQ -F MQRankSum -F ReadPosRankSum \
     -GF AD \
     -O ${strain}_snps_stats.tsv \
     -raw  

gatk VariantsToTable \
     -R cascadePrimary.fasta \
     -V ${strain}_filtered_indels.vcf \
     -F CHROM -F POS -F TYPE -F FILTER \
     -F QD -F FS -F ReadPosRankSum \
     -GF AD \
     -O ${strain}_indels_stats.tsv \
     -raw 

# table with snp types and genotype fields
gatk VariantsToTable \
     -R cascadePrimary.fasta \
     -V ${strain}_filtered_snps.vcf \
     -F CHROM -F POS -F REF -F ALT -F HET -F HOM-VAR -F HOM-REF -F NO-CALL -F MULTI-ALLELIC -F FILTER\
     -GF GT\
     -O cohort_filtered_snps_genotypes.tsv \
     -raw  


cp -r $WorkDir/*.tsv $outdir
rm -r $WorkDir
