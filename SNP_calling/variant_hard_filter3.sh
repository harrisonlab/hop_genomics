#!/usr/bin/env bash
#SBATCH -J filter_variants
#SBATCH --partition=long
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=30


##########################################################################
#steps
# 1. Extract the SNPs from the call set
# 2. Determine parameters for filtering SNPs
# 3. Apply the filter to the SNP call set
# 4. Extract the Indels from the call set
# 5. Determine parameters for filtering indels
# 6. Apply the filter to the Indel call set


#INPUT:
# 1st argument: Refrence fasta 
# 2nd argument: sample name (prefix) to be used to identify it in the future
# 3nd argument: Input VCF with unfiltered SNP and INDEL calls
#OUTPUT:
# filtered VCF


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

gatk SelectVariants \
     -R cascadePrimary.fasta \
     -V ${strain}_SNPs_genotype.vcf \
     -select-type SNP \
     -O ${strain}_raw_snps.vcf

gatk VariantFiltration \
    -V ${strain}_raw_snps.vcf \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    --missing-values-evaluate-as-failing \
    -O ${strain}_filtered_snps.vcf

gatk SelectVariants \
     -R cascadePrimary.fasta \
     -V ${strain}_SNPs_genotype.vcf \
     -select-type INDEL \
     -O ${strain}_raw_indels.vcf

gatk VariantFiltration \
    -V ${strain}_raw_indels.vcf \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    --missing-values-evaluate-as-failing \
    -O ${strain}_filtered_indels.vcf 

gatk CountVariants \
     -V *.vcf

cp -r $WorkDir/*.vcf $outdir
rm -r $WorkDir
