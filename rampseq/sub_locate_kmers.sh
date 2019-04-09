#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 2
#$ -l virtual_free=1G
#$ -l h=blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace


Usage="sub_locate_kmers <assembly.fasta> <filtered_kmers.tsv> <output_path/prefix>"
echo "$Usage"

Contig=$1
Kmers=$2
Prefix=$3

ProgDir=/home/armita/git_repos/emr_repos/scripts/hop_genomics/rampseq
$ProgDir/locate_kmers.py --assembly $Contig --kmers $Kmers --prefix $Prefix
