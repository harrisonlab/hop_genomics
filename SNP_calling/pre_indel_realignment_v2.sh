#!/usr/bin/env bash
#SBATCH -J pre_indel_realignment
#SBATCH --partition=long
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=30


##########################################################################
#INPUT:
# 1st argument: Refrence fasta 
# 2nd argument:input BAM file from pre_snp_calling file with your mappings with duplicates marked no multimapping sorted 
# 3rd argument: sample name (prefix) to be used to identify bam in the future
#OUTPUT:
# realigner intervals to be used in the next stages of indel realignment 


reference=$1
input_bam=$2
strain=$3
outdir=$4


WorkDir=/projects/hop/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir


cp $reference $WorkDir
cp $input_bam $WorkDir
cd $WorkDir


samtools faidx cascadePrimary.fasta
samtools index *_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam


picard=/home/connellj/miniconda2/share/picard-2.18.29-0/picard.jar
java -jar $picard CreateSequenceDictionary \
	R=cascadePrimary.fasta \
	O=cascadePrimary.dict 


gatk RealignerTargetCreator \
     -R cascadePrimary.fasta \
     -I *_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -o realigner.intervals 


cp -r $WorkDir/realigner.intervals $outdir
rm -r $WorkDir