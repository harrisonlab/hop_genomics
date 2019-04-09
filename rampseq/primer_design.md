# Commands used in the design of RAMPseq primers


```bash
  cd /home/groups/harrisonlab/project_files/hop_genomics
```

## Download reference genomes:


### Cascade assembly from Hopbase

This assembly was generated from PacBio reads using Falcon

```bash
CurDir=$PWD
OutDir=assembly/reference/H.lupulus/cascade/hopbase
mkdir -p $OutDir
cd $OutDir
wget http://hopbase.cgrb.oregonstate.edu/content/dedupCascade/cascadePrimary.fasta.gz
wget http://hopbase.cgrb.oregonstate.edu/content/maskedDedupCascade/sorted_modified_combinedMipsREdatRepMask_LTRharvestRepMask_GBrowse.gff.gz
wget http://hopbase.cgrb.oregonstate.edu/content/dedupCascade/cascadeHenningStringtieTranscripts.gff3.gz
wget http://hopbase.cgrb.oregonstate.edu/content/dedupCascade/dedupHopBaseCascade_vs_uniprot.gff3.gz
wget http://hopbase.cgrb.oregonstate.edu/content/cascadeHaplotigs/cascadeHaplotigs.fasta.gz
# unzip data:
gunzip *
```


###




# Kmer filtering:


<!--
Test pipelines using Alternaria (small fungal repeat-sparse genome)

```bash
qlogin -pe smp 8 -l virtual_free=5.9G
cd /home/groups/harrisonlab/project_files/hop_genomics
# Test runs using Alternaria data:

Assembly=$(ls /data/scratch/armita/alternaria/repeat_masked/A.gaisen/650/filtered_contigs/650_contigs_unmasked.fa)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=rampseq/primer_design/kat/$Organism/$Strain
mkdir -p $OutDir

Threads=8
MinCount=1
MaxCount=2000
MinGC=1
MaxGC=100
KmerSz=20

Prefix=Alt_kmers

kat filter kmer --output_prefix  $OutDir/$Prefix --threads $Threads --low_count $MinCount --high_count $MaxCount --low_gc $MinGC --high_gc $MaxGC --non_canonical --mer_len $KmerSz $Assembly

/home/armita/prog/kat/bin/kat_jellyfish info $OutDir/${Prefix}-in.jf20
/home/armita/prog/kat/bin/kat_jellyfish dump $OutDir/${Prefix}-in.jf20 --lower-count $MinCount --upper-count $MaxCount --column > $OutDir/$Prefix.tsv

# K-mers in input   : 33911372 distinct; 34346437 total.
# K-mers to keep    : 33910583 distinct; 34345455 total.

Assembly=$(ls /data/scratch/armita/alternaria/repeat_masked/A.gaisen/650/filtered_contigs/650_contigs_softmasked_repeatmasker_TPSI_appended.fa)

Threads=8
MinCount=1
MaxCount=2000
MinGC=1
MaxGC=100
KmerSz=20

Prefix=Alt_kmers_soft

kat filter kmer --output_prefix  $OutDir/$Prefix --threads $Threads --low_count $MinCount --high_count $MaxCount --low_gc $MinGC --high_gc $MaxGC --non_canonical --mer_len $KmerSz $Assembly

/home/armita/prog/kat/bin/kat_jellyfish info $OutDir/${Prefix}-in.jf20
/home/armita/prog/kat/bin/kat_jellyfish dump $OutDir/${Prefix}-in.jf20 --lower-count $MinCount --upper-count $MaxCount --column > $OutDir/$Prefix.tsv

# K-mers in input   : 33911372 distinct; 34346437 total.
# K-mers to keep    : 33910583 distinct; 34345455 total.

Threads=8
MinCount=100
MaxCount=10000
MinGC=1
MaxGC=100
KmerSz=20

Prefix=Alt_kmers_soft

kat filter kmer --output_prefix  $OutDir/$Prefix --threads $Threads --low_count $MinCount --high_count $MaxCount --low_gc $MinGC --high_gc $MaxGC --non_canonical --mer_len $KmerSz $Assembly

jellyfish info $OutDir/${Prefix}-in.jf20
jellyfish dump $OutDir/${Prefix}-in.jf20 > $OutDir/$Prefix.fa
```
-->

## Identify kmer copy number and perform filtering:

```bash
qlogin -pe smp 8 -l virtual_free=5.9G
cd /home/groups/harrisonlab/project_files/hop_genomics
Assembly=$(ls assembly/reference/H.lupulus/cascade/hopbase/cascadePrimary.fasta)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=/data/scratch/armita/hop_genomics/rampseq/primer_design/kat/$Organism/$Strain
mkdir -p $OutDir

Threads=8
MinCount=1
MaxCount=10000
# MinGC=35
# MaxGC=65
MinGC=1
MaxGC=100
KmerSz=20

Prefix=cascade_kmers


kat filter kmer --output_prefix  $OutDir/$Prefix --threads $Threads --low_count $MinCount --high_count $MaxCount --low_gc $MinGC --high_gc $MaxGC --non_canonical --mer_len $KmerSz $Assembly


MinCount=1000
MaxCount=10000
/home/armita/prog/kat/bin/kat_jellyfish info $OutDir/${Prefix}-in.jf20
/home/armita/prog/kat/bin/kat_jellyfish dump $OutDir/${Prefix}-in.jf20 --lower-count $MinCount --upper-count $MaxCount --column > $OutDir/$Prefix.tsv

ProgDir=/home/armita/git_repos/emr_repos/scripts/hop_genomics/rampseq
# $ProgDir/process_kmers.py --kmer $OutDir/$Prefix.tsv --min_copy 1000 --min_GC 35 --max_GC 65 > $OutDir/${Prefix}_filtered.tsv
$ProgDir/process_kmers.py --kmer $OutDir/$Prefix.tsv --min_copy 1000 --min_GC 35 --max_GC 65 > $OutDir/${Prefix}_filtered.fa
```

## Identify kmer locations

Test scripts:

```bash
qlogin -pe smp 1 -l virtual_free=2G
cd /home/groups/harrisonlab/project_files/hop_genomics
Assembly=$(ls assembly/reference/H.lupulus/cascade/hopbase/cascadePrimary.fasta)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
# OutDir=/data/scratch/armita/hop_genomics/rampseq/primer_design/bwa/$Organism/$Strain
OutDir=/data/scratch/armita/hop_genomics/rampseq/primer_design/kat/$Organism/$Strain
mkdir -p $OutDir

Prefix=${Strain}_kmers
Kmers=$(ls /data/scratch/armita/hop_genomics/rampseq/primer_design/kat/$Organism/$Strain/${Prefix}_filtered.fa)
ProgDir=/home/armita/git_repos/emr_repos/scripts/hop_genomics/rampseq

# $ProgDir/locate_kmers.py --assembly $Assembly --kmers $OutDir/${Prefix}_filtered.fa | less
cat $OutDir/${Prefix}_filtered.fa | head -n2000 > $OutDir/${Prefix}_filtered_test.fa
$ProgDir/locate_kmers.py --assembly $OutDir/test.fa --kmers $OutDir/${Prefix}_filtered_test.fa --prefix $OutDir/${Strain}

cat $OutDir/cascade_000000F_kmer_pairs.tsv | cut -f1 | sort | uniq -c | sort -nr
cat $OutDir/cascade_000000F_kmer_pairs.tsv | grep '159_3798_40.0-909_3912_35.0' | cut -f1,7 | sed "s/^/>/g" | sed "s/\t/\n/g" | less

$ProgDir/locate_kmers.py --assembly $OutDir/000000F.fa --kmers $OutDir/${Prefix}_filtered.fa --prefix $OutDir/${Strain}_2


# cat $OutDir/cascade_000000F_kmer_pairs.tsv | cut -f1 | sort | uniq -c | sort -nr | head
#
# cat $OutDir/cascade_000000F_kmer_pairs.tsv | grep '52241_2216_65.0-69655_2605_65.0' | cut -f1,7 | sed "s/^/>/g" | sed "s/\t/\n/g" | less
#
# cat $OutDir/cascade_000000F_kmer_pairs.tsv | grep '31385_7521_35.0-38330_9906_40.0' | cut -f1,7 | sed "s/^/>/g" | sed "s/\t/\n/g" | less

cat $OutDir/cascade_000000F_2_kmer_pairs.tsv | cut -f1 | sort | uniq -c | sort -nr | head
cat /data/scratch/armita/hop_genomics/rampseq/primer_design/kat/H.lupulus/cascade/cascade_2_000000F_kmer_pairs.tsv | grep '2294_8271_40.0-70797_7373_35.0' | cut -f1,7 | sed "s/^/>/g" | sed "s/\t/\n/g" | less
```

Submit analysis for whole genome:

This work was done in the temporary /data dirve


Make a local copy of the assembly
```bash
cd /data/scratch/armita/hop_genomics
Assembly=$(ls /home/groups/harrisonlab/project_files/hop_genomics/assembly/reference/H.lupulus/cascade/hopbase/cascadePrimary.fasta)
mkdir -p assembly/reference/H.lupulus/cascade/hopbase
cp -s $Assembly assembly/reference/H.lupulus/cascade/hopbase/.
```

Split the assembly into individual contigs:

```bash
qlogin
cd /data/scratch/armita/hop_genomics
Assembly=$(ls assembly/reference/H.lupulus/cascade/hopbase/cascadePrimary.fasta)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
SplitDir=$(dirname $Assembly)/split
mkdir $SplitDir
ProgDir=/home/armita/git_repos/emr_repos/scripts/hop_genomics/rampseq
$ProgDir/split_fasta.py --fasta $Assembly --prefix $SplitDir/${Strain}
```

for each contig identify kmers present:

```bash
  cd /data/scratch/armita/hop_genomics
  for Contig in $(ls assembly/reference/H.lupulus/cascade/hopbase/split/*.fa | tail -n+2); do
    Organism=$(echo $Contig | rev | cut -f5 -d '/' | rev)
    Strain=$(echo $Contig | rev | cut -f4 -d '/' | rev)
    SplitDir=$(dirname $Contig)/split
    OutDir=rampseq/primer_design/locations/$Organism/$Strain/split
    mkdir -p $OutDir
    Prefix=${Strain}_kmers
    Kmers=$(ls rampseq/primer_design/kat/$Organism/$Strain/${Prefix}_filtered.fa)
    Jobs=$(qstat | grep 'sub_locate' | grep 'qw' | wc -l)
    while [ $Jobs -gt 1 ]; do
      sleep 30s
      printf "."
      Jobs=$(qstat | grep 'sub_locate' | grep 'qw' | wc -l)
    done		
    printf "\n"
    echo $(basename $Contig)
    ProgDir=/home/armita/git_repos/emr_repos/scripts/hop_genomics/rampseq
    qsub $ProgDir/sub_locate_kmers.sh $Contig $Kmers $OutDir/$Strain
  done
```


```bash
cd /data/scratch/armita/hop_genomics
ls rampseq/primer_design/locations/H.lupulus/cascade/split/cascade_*_kmer_pairs.tsv | wc -l
for File in $(ls rampseq/primer_design/locations/H.lupulus/cascade/split/cascade_*_kmer_pairs.tsv); do
  basename $File
  OutName=$(echo $File | sed 's/.tsv/_paircounts.txt/g')
  cat $File | cut -f1 | sort | uniq -c | sort -nr > $OutName
done


OutFile=rampseq/primer_design/locations/H.lupulus/cascade/paircounts_min100.txt
ProgDir=/home/armita/git_repos/emr_repos/scripts/hop_genomics/rampseq
$ProgDir/sum_primer_occurence.py --inF $(ls rampseq/primer_design/locations/H.lupulus/cascade/split/cascade_*_kmer_pairs_paircounts.txt) --min 100 > $OutFile

wc -l $OutFile
head -n100 $OutFile
```

```
1184574
```

```
11542_6815_60.0-24312_9039_55.0	10331
24312_9039_55.0-61464_6814_60.0	10324
24312_9039_55.0-53689_6810_60.0	10317
4052_9042_55.0-53689_6810_60.0	10301
18830_7612_35.0-55267_8278_40.0	10296
11084_7253_50.0-53259_7961_55.0	10221
29789_8117_40.0-60968_7670_35.0	10220
54833_9189_45.0-68374_8565_50.0	10209
16293_9961_60.0-49159_6513_40.0	10209
19803_8137_45.0-53259_7961_55.0	10206
14921_8691_55.0-19792_9624_40.0	10198
4052_9042_55.0-49924_6966_55.0	10192
29789_8117_40.0-42424_7681_35.0	10189
11084_7253_50.0-65144_7947_55.0	10185
19803_8137_45.0-65144_7947_55.0	10177
11084_7253_50.0-42508_7943_55.0	10172
19803_8137_45.0-42508_7943_55.0	10162
2294_8271_40.0-60968_7670_35.0	10145
10554_7539_35.0-40267_8424_35.0	10143
38330_9906_40.0-60968_7670_35.0	10143
18830_7612_35.0-29789_8117_40.0	10132
2294_8271_40.0-42424_7681_35.0	10128
49924_6966_55.0-66601_8144_50.0	10125
38330_9906_40.0-42424_7681_35.0	10121
14921_8691_55.0-49159_6513_40.0	10120
68374_8565_50.0-8781_9023_45.0	10096
46424_8016_45.0-60968_7670_35.0	10094
42424_7681_35.0-46424_8016_45.0	10069
55372_9751_40.0-60968_7670_35.0	10066
18830_7612_35.0-2294_8271_40.0	10061
18830_7612_35.0-38330_9906_40.0	10051
```

```bash
Pair="11542_6815_60.0-24312_9039_55.0"
OutDir=rampseq/primer_design/locations/H.lupulus/cascade/$Pair
mkdir $OutDir
for File in $(ls rampseq/primer_design/locations/H.lupulus/cascade/split/cascade_*_kmer_pairs.tsv); do
  echo $(basename $File)
  cat $File | grep -w "$Pair" >> $OutDir/${Pair}_information.tsv
done
cat $OutDir/${Pair}_information.tsv | cut -f7 | sort | uniq -c | sort -n > $OutDir/${Pair}_counts.txt
echo "Number of unique sequences:"
cat $OutDir/${Pair}_counts.txt | grep -w '1' | wc -l
cat $OutDir/${Pair}_information.tsv | cut -f1,7 | sed "s/^/>/g" | sed "s/\t/\n/g" > $OutDir/${Pair}_seqs.fa
```

```
  Number of unique sequences:
  7902
```

If this locus has ~10000 amplicons and a MiSeq can generate up to 25million reads
also a MiSeq can cope with 96 illumina barcodes in a single run.

25000000 / 96 = 260416 reads per sample
260416 / 10000 = 26 reads per amplicon for each sample.

This means we would want to use a minimum 3 MiSeq flow cells for a single population with this locus.
