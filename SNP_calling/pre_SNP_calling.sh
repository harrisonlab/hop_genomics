


-----------------------------------------------------------------Alignment---------------------------------------------------------------------
# Alignemt of S1 S3 reads vs Cascade genome

 Alignment of reads from a single run using Bowtie:

for Hop_genome in S1 S3; do
  Reference=$(ls ../../projects/hop/raw_data/cascasde/cascadePrimary.fasta)
  for StrainPath in $(ls -d ../../projects/hop/raw_data/trim/"$Hop_genome"); do
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    F_Read=$(ls $StrainPath/F/*_trim.fq.gz)
    R_Read=$(ls $StrainPath/R/*_trim.fq.gz)
    echo "$Organism - $Strain"
    echo $F_Read
    echo $R_Read
    OutDir=$StrainPath/alignment
    mkdir -p $OutDir
    ProgDir=/home/hajduk/gitrepos/scripts
    sbatch $ProgDir/bowtie.sh $Reference $F_Read $R_Read $OutDir
  done
done

#the line below converts bam files to sam files 
for Strain in S1 S3; do 
 for file in /projects/hop/raw_data/trim/alignment/cascadePrimary.fasta_aligned_sorted.bam
  do
     echo $file 
     samtools view -h $file > ${file/.bam/.sam}
 done
done

-------------------------------------------------------------Renaming_files------------------------------------------------------------------
## Rename input mapping files in each folder by prefixing with the strain ID

```bash
for Strain in S3 ; do
  for filename in /projects/hop/raw_data/trim/$Strain/alignment; do
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
```

--------------------------------------------------------------RemoveMultimapReadss----------------------------------------------------------

## Remove multimapping reads, discordant reads. PCR and optical duplicates, and add read group and sample name to each mapped read (preferably, the shortest ID possible)

Convention used:
qsub $ProgDir/sub_pre_snp_calling.sh <INPUT SAM FILE> <SAMPLE_ID>

```bash
  for Strain in S1; do # Replace with the strain name
    for input in ../../projects/hop/raw_data/trim/$Strain/alignment/"$Strain"_unmasked.fa_aligned.sam; do
  #   Jobs=$(squeue -u ${USER} --noheader --array | wc -l)
  #      while [ $Jobs -gt 1 ]; do
  #          sleep 25m
  #          printf "."
  #          Jobs=$(squeue -u ${USER} --noheader --array | wc -l)
  #      done
  #  printf "\n"
    OutDir=../../projects/hop/raw_data/trim/$Strain/alignment/nomulti
  #  mkdir -p $OutDir  
    ProgDir=/home/hajduk/gitrepos/scripts/SNP_calling
    sbatch $ProgDir/pre_SNP_calling.sh $input $Strain $OutDir # This will add read group and sample name to each mapped read. Preferably, use the shortest ID possible.
    done
  done
 ```  

 