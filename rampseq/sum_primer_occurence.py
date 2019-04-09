#!/usr/bin/python
#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import os
import argparse
import re
from collections import defaultdict

ap = argparse.ArgumentParser()

ap.add_argument('--inF', required=True,type=str, nargs='+',help='director containing file with _paircounts.txt endings')
ap.add_argument('--min', required=True,type=int,help='minimum occurence of pairs printed to output')

# ap.add_argument('--out', required=True,type=str,help='filepath and filename for the output file')

conf = ap.parse_args()

# out = conf.out

count_dict = defaultdict(int)

inList = conf.inF
minC = conf.min

for file in inList:
    with open(file) as f:
        kmer_lines = f.readlines()
    for line in kmer_lines:
        line = line.rstrip()
        line = re.sub(r"^\s+", "", line)
        line = re.sub(r" +", " ", line)
        # print line
        split_line = line.split(" ")
        count = split_line[0]
        pair = split_line[1]
        count_dict[pair] += int(count)

sorted_keys = sorted(count_dict, reverse=True, key=lambda k: count_dict[k])
for key in sorted_keys:
    count = count_dict[key]
    if count < minC:
        continue
    print "\t".join([key, str(count)])
