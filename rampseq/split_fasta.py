#!/usr/bin/python
#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import os
import argparse
from Bio import SeqIO

ap = argparse.ArgumentParser()

ap.add_argument('--fasta', required=True,type=str,help='Input fasta file')
ap.add_argument('--prefix', required=True,type=str,help='filepath and prefix for output files')

conf = ap.parse_args()

seq_records = list(SeqIO.parse(conf.fasta, "fasta"))
out = conf.prefix

#-----------------------------------------------------
# Step 2
# Write sequences
#-----------------------------------------------------

for record in seq_records:
    seq_id = record.id
    o = "_".join([out, seq_id + ".fa"])
    SeqIO.write(record, o, "fasta")
