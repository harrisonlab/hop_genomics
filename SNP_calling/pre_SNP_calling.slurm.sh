#!/usr/bin/env bash
#SBATCH -J pre_snp_calling
#SBATCH --partition=long
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=8


##############################################
# Prep mappings from Bowtie2 for SNP calling
### Remove multimapping reads, discordant reads. PCR and optical duplicates, and
### add read group and sample name to each mapped read (preferably, the shortest ID possible)
#INPUT:
# 1st argument: input SAM file with your mappings
# 2nd argument: sample name (prefix) to be used to identify it in the future
#OUTPUT:
# Indexed BAM file with suffix "nodup_rg" to be fed into SNP calling with GATK.
#############################################
input_sam=$1
prefix=$2
OutDir=$3
filename=$(basename "$input_sam")
name="${filename%.*}"


WorkDir=/projects/hop/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
### Prep

cp $input_sam $WorkDir
cd $WorkDir

### Get rid of multimapping reads by filtering out on the XS:i: tag
grep -v "XS:i" $filename >temp && mv temp $filename
samtools view -bS -o $name.bam $filename
samtools sort $name.bam -o $name\_nomulti\_sorted.bam
#For samtools older versions (0.1.18) use: 
#samtools sort $name.bam $name\_nomulti\_sorted
samtools index $name\_nomulti\_sorted.bam

### Keep only reads with "paired reads" and "properly paired reads" flags.
samtools view -b -h -f 3 -o $name\_nomulti\_proper.bam $name\_nomulti\_sorted.bam
### Sort for downstream analyses
samtools sort $name\_nomulti\_proper.bam -o $name\_nomulti\_proper\_sorted.bam
#For samtools older versions (0.1.18) use: 
#samtools sort $name\_proper.bam $name\_proper\_sorted
samtools index $name\_nomulti\_proper\_sorted.bam

### Remove PCR and optical duplicates
picard=/home/connellj/miniconda2/share/picard-2.18.29-0/picard.jar
java -jar $picard MarkDuplicates \
	I=$name\_nomulti\_proper\_sorted.bam \
	O=$name\_nomulti\_proper\_sorted\_nodup.bam \
	METRICS_FILE=$name\_nomulti\_proper\_sorted\_nodup.txt \
	REMOVE_DUPLICATES=TRUE \
	ASSUME_SORTED=TRUE \
	MAX_RECORDS_IN_RAM=500000000 \
	MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
	VALIDATION_STRINGENCY=LENIENT 


### Add group and sample name (prefix)
java -jar $picard AddOrReplaceReadGroups \
	I=$name\_nomulti\_proper\_sorted\_nodup.bam \
	O=$name\_nomulti\_proper\_sorted\_nodup_rg.bam \
	SORT_ORDER=coordinate \
	CREATE_INDEX=true \
	RGID=$prefix \
	RGSM=$prefix \
	RGPL=Illumina \
	RGLB=library \
	RGPU=barcode \
	VALIDATION_STRINGENCY=LENIENT 


samtools index $name\_nomulti\_proper\_sorted\_nodup_rg.bam

### Cleanup
mv $WorkDir/$name\_nomulti\_proper\_sorted\_nodup.txt $OutDir
mv $WorkDir/$name\_nomulti\_proper\_sorted\_nodup_rg.bam $OutDir
mv $WorkDir/$name\_nomulti\_proper\_sorted\_nodup_rg.bam.bai $OutDir
#rm -rf $WorkDir

