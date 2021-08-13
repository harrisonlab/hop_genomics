for RawData in $(ls /projects/hop/raw_data/S1/*/*.fastq.gz); do
echo $RawData;
ProgDir=/home/hajduk/gitrepos/scripts/dna_qc$ less run_fastqc.sh
sbatch $ProgDir/run_fastqc.sh $RawData; $OutDir
done

