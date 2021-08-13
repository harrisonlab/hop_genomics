#!/usr/bin/env bash
#SBATCH -J haplotype_caller
#SBATCH --partition=himem 
#SBATCH --mem-per-cpu=3G
#SBATCH --cpus-per-task=6


##########################################################################
#INPUT:
# 1st argument: Refrence fasta 
# 2nd argument: sample name (prefix) to be used to identify it in the future
# 3nd argument: Input BAM 
#OUTPUT:
# VCF out


reference=$1
strain=$2
input_bam=$3
outdir=$4


WorkDir=/projects/hop/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir


cp $reference $WorkDir
cp $input_bam $WorkDir
cd $WorkDir


samtools faidx cascadePrimary.fasta
samtools index "$Strain"_realigned.bam


picard=/home/connellj/miniconda2/share/picard-2.18.29-0/picard.jar
java -jar $picard CreateSequenceDictionary \
	R=cascadePrimary.fasta \
	O=cascadePrimary.dict 


gatk HaplotypeCaller \
     -ploidy 2 \
     -R cascadePrimary.fasta \
     -I "$strain"_realigned.bam \
     -o "$strain"_SNP_calls.g.vcf \
     -ERC GVCF \
     --allow_potentially_misencoded_quality_scores \
     -variant_index_type LINEAR \
     -variant_index_parameter 128000

 

cp $WorkDir/"$strain"_SNP_calls.g.vcf $outdir
rm -r $WorkDir



#the hajdooo guide to snp calling 


#1.) calll snips 



#2.) the end 



#3.) ask john why it dosent work 



#4.) use bowtie2 to make algnments 




#5.) haplotypecaller 



#6.) segg dugo 