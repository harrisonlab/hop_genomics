#!/usr/bin/env bash
#SBATCH -J genome_coverage
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=50G
#SBATCH --cpus-per-task=3
# ---------------
# INPUT:
# 1st argument: Forward read
# 2nd argument: Reverse read
# 3rd argument: Genome size
# OUTPUT:
# x Genome coverage
Read_F=$(basename $1)
Read_R=$(basename $2)
Genome_size=$3
OutDir=$4
#DATA_TYPE=$(echo $Read_F | rev | cut -d "/" -f6 | rev)
#READ_TYPE=$(echo $Read_F | rev | cut -d "/" -f5 | rev)
#ORGANISM=$(echo $Read_F | rev | cut -d "/" -f4 | rev)
#STRAIN=$(echo $Read_F | rev | cut -d "/" -f3 | rev)
# ---------------
# Step 2
# Copy data
# ---------------
WorkDir=/projects/hop/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
#cd $WorkDir

cp $1 $WorkDir
cp $2 $WorkDir
cd $WorkDir
gunzip $Read_F
gunzip $Read_R
# Step 3
# Calculate estimated genome coverage

Sub1=*1_trim.fq
Sub2=*2_trim.fq

/data/scratch/gomeza/prog/count_nucl.pl -i $Sub1 -i $Sub2 -g $3 > estimated_coverage.log
cp $WorkDir/estimated_coverage.log $OutDir
rm -r $WorkDir

