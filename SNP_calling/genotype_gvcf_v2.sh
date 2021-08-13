#!/usr/bin/env bash
#SBATCH -J genotype_gvcf
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=3G
#SBATCH --cpus-per-task=6


##########################################################################
#SNP call files and WT genome should be collected into the same diretory, WT.dict and WT.fai must also be present. 


in_file=/projects/hop/SNP_calling/haplotype_calls_v2
reference=/projects/hop/SNP_calling/haplotype_calls_v2/cascadePrimary.fasta 


filename=$(basename "$reference")
output=/projects/hop/SNP_calling/haplotype_calls_v2/"${filename%.*}_SNPs_genotype.vcf"
output2=/projects/hop/SNP_calling/haplotype_calls_v2/"${filename%.*}_SNPs_genotype_2.vcf"


gatk GenotypeGVCFs \
     -R $reference \
     -V $in_file/S1_SNP_calls.g.vcf \
	-V $in_file/S3_SNP_calls.g.vcf \
 	-o $output	


												#Perform joint genotyping on one or more samples pre-called with HaplotypeCaller diploid organism \
gatk GenotypeGVCFs \
     -T VariantsToAllelicPrimitives \
     -R $reference \
     -V $output \
     -o $output2								#A VCF with alleles broken into primitive type

#Simplify multi-nucleotide variants (MNPs) into more basic/primitive alleles