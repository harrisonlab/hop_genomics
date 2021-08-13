#!/bin/bash
#!/usr/bin/env bash
#SBATCH -J fastq-mcf
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=60G
#SBATCH --cpus-per-task=3

. ~/.profile
# Commands to perform adapter trimming of either RNA of DNA data using trimmomatic
# Please visualise the quality of data before and after this trimming using fastqc.
USAGE="qc_trimmomatic.sh <F_reads.fq/fq.gz> <R_reads.fq/fq.gz> <Illumina_adapters.fa> <Output_directory> <logfile prefix>"
echo "$USAGE"
echo ""
F_IN=$1
R_IN=$2
ADAPTER_FILE=$3
OutDir=$4
Prefix=$5
CUR_PATH=$PWD
WORK_DIR=/projects/hop/${SLURM_JOB_ID}
#WORK_DIR=/mnt/beegfs/scratch/$USER/${SLURM_JOB_ID}
mkdir -p $WORK_DIR
cd $WORK_DIR
LOGFILE=${Prefix}_trim_log.txt
F_NO_ADAPT=${Prefix}_F_trim.fq.gz
R_NO_ADAPT=${Prefix}_R_trim.fq.gz
UNPAIRED_F_NO_ADAPT=${Prefix}_F_trim_unpaired.fq.gz
UNPAIRED_R_NO_ADAPT=${Prefix}_R_trim_unpaired.fq.gz
java -Xmx1g -jar /projects/oldhome/armita/prog/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar \
  PE -phred33 -threads 10 -trimlog $LOGFILE \
  $CUR_PATH/$F_IN $CUR_PATH/$R_IN \
  $F_NO_ADAPT $UNPAIRED_F_NO_ADAPT \
  $R_NO_ADAPT $UNPAIRED_R_NO_ADAPT \
  ILLUMINACLIP:$ADAPTER_FILE:2:30:10 \
  LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36
# Remove Illumina adapters provided in the TruSeq3-PE.fa file (provided). Initially Trimmomatic will look for seed matches (16 bases) allowing maximally 2 mismatches. These seeds will be extended and clipped if in the case of paired end reads a score of 30 is reached (about 50 bases), or in the case of single ended reads a score of 10, (about 17 bases).
# Remove leading low quality or N bases (below quality 5)
# Remove trailing low quality or N bases (below quality 5)
# Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15
# Drop reads which are less than 36 bases long after these steps
OutDirF=$CUR_PATH/$OutDir/F
OutDirR=$CUR_PATH/$OutDir/R
mkdir -p $OutDirF
cp $F_NO_ADAPT $OutDirF/.
cp $UNPAIRED_F_NO_ADAPT $OutDirF/.
mkdir -p $OutDirR
cp $R_NO_ADAPT $OutDirR/.
cp $UNPAIRED_R_NO_ADAPT $OutDirR/.
cp $LOGFILE $CUR_PATH/$OutDir/.
rm -r $WORK_DIR
