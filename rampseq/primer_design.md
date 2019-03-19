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
wget http://hopbase.cgrb.oregonstate.edu/content/dedupCascade/cascadeHenningStringtieTranscripts.gff3.gz
wget http://hopbase.cgrb.oregonstate.edu/content/dedupCascade/dedupHopBaseCascade_vs_uniprot.gff3.gz
wget http://hopbase.cgrb.oregonstate.edu/content/cascadeHaplotigs/cascadeHaplotigs.fasta.gz
# unzip data:
gunzip *
```


###




# Kmer filtering:


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

Identify kmer locations:

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

# bwa index $Assembly
#
# bwa mem -M -t 12 $Assembly $Kmers | samtools view -S -b - > $OutDir/"$Prefix".bam
#
# ### Add group and sample name (Prefix)
# bamaddrg -b "$Prefix".bam -s $Prefix -r $Prefix > "$Prefix"_unsorted.bam
# ### Sort the full BAM file.
# samtools sort -@ 12  "$Prefix"_unsorted.bam -o "$Prefix"_sorted.bam
#
# #index
# samtools index "$Prefix"_sorted.bam

```
