#!/usr/bin/env bash
#SBATCH --job-name=bash
#SBATCH --partition=long
#SBATCH --time=0-20:00:00



INFILE=$1
DATA_TYPE=$(echo $INFILE | rev | cut -d "/" -f6 | rev)
READ_TYPE=$(echo $INFILE | rev | cut -d "/" -f5 | rev)
ORGANISM=$(echo $INFILE | rev | cut -d "/" -f4 | rev)
STRAIN=$(echo $INFILE | rev | cut -d "/" -f3 | rev)
DIRECTION=$(echo $INFILE | rev | cut -d "/" -f2 | rev)
READS=$(echo $INFILE | rev | cut -d "/" -f1 | rev)

CUR_PATH=$PWD
WORK_DIR=$TMPDIR/"$STRAIN"_fastqc

mkdir -p $WORK_DIR
cd $WORK_DIR

cp $CUR_PATH/$INFILE .

fastqc --nogroup $READS

mkdir -p $CUR_PATH/$DATA_TYPE/$READ_TYPE/$ORGANISM/$STRAIN/$DIRECTION/

cp -r $WORK_DIR/*/ $CUR_PATH/$DATA_TYPE/$READ_TYPE/$ORGANISM/$STRAIN/$DIRECTION/.

rm -r $TMPDIR

exit
