

---------------------------------------------------Quality_assessment_pre_trimming------------------------------------------------------------------------------------

Data qc 
programs: fastqc fastq-mcf kmc
Data quality was visualised using fastqc:


#Assess read qulity and genome coverage 

# Run fastqc
for Strain in S1 S3; do
    RawData=$(ls ../../projects/hop/raw_data/$Strain/R/*.fq)
    echo $RawData;
    ProgDir=/home/hajduk/gitrepos/scripts/dna_qc
    sbatch $ProgDir/fastqc.sh $RawData 
done


-------------------------------------------------------------------trimming------------------------------------------------------------------------------------

Trimming was performed on data to trim adapters from sequences and remove poor quality data. This was done with fastq-mcf

  # Run fastq-mcf
for Strain in S3; do
    Read_F=../../projects/hop/raw_data/$Strain/F/*.fq
    Read_R=../../projects/hop/raw_data/$Strain/R/*.fq 
    echo $Read_F;
    echo $Read_R;
    IluminaAdapters=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc/illumina_full_adapters.fa
    ProgDir=/home/hajduk/gitrepos/scripts/dna_qc
    sbatch $ProgDir/fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA 
done

-------------------------------------------------------------------GenomeCoverage------------------------------------------------------------------------------------

# Assess genome sequencing coverage any coverage less then 30x should be excluded from SNP calling

#Change file permissions from archive to executble
chmod +x /projects/hop/raw_data/trim/S3/*/*.fq.gz


for Strain in S3; do
 for DataDir in $(ls -d ../../projects/hop/raw_data/trim/$Strain); do
    F_Read=$(ls $DataDir/F/*fq.gz)
    R_Read=$(ls $DataDir/R/*fq.gz)
    Outdir=../../projects/hop/raw_data/trim/$Strain/genome_coverage
    echo $F_Read
    echo $R_Read
    echo $Outdir
    mkdir -p $Outdir
    ProgDir=/home/hajduk/gitrepos/scripts/dna_qc
    sbatch $ProgDir/count_nucl.sh $F_Read $R_Read 3000 $Outdir #Estimated genome size
 done
done


-------------------------------------------------------------------Visualisation-------------------------------------------------------------------------------

Data qc 
programs: fastqc fastq-mcf kmc
Data quality was visualised using fastqc:


for Strain in S3; do
    RawData=$(ls ../../projects/hop/raw_data/trim/$Strain/R/*_trim.fq.gz)
    echo $RawData;
    ProgDir=/home/hajduk/gitrepos/scripts/dna_qc
    sbatch $ProgDir/fastqc.sh $RawData 
done



