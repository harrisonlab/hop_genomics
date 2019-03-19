#!/usr/bin/python
#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import os
import argparse
import re
import math
import numpy as np
from sets import Set
from collections import defaultdict
from collections import Counter
from operator import itemgetter

ap = argparse.ArgumentParser()

ap.add_argument('--kmer', required=True,type=str,help='tsv file output from KAT (kat_jellyfish dump) containing non-canonical kmer sequence and copy number')
ap.add_argument('--min_copy', required=True,type=int,help='Minimum sequence copy number')
# ap.add_argument('--max_copy', required=True,type=int,help='Maximum sequence copy number')
ap.add_argument('--min_GC', required=True,type=int,help='Minimum GC content')
ap.add_argument('--max_GC', required=True,type=int,help='Maximum GC content')
# ap.add_argument('--prefix', required=True,type=str,help='filepath and prefix for output files')

conf = ap.parse_args()

min_copy = conf.min_copy
min_GC = conf.min_GC
max_GC = conf.max_GC
with open(conf.kmer) as f:
    kmer_lines = f.readlines()

#-----------------------------------------------------
# Step 2
# Define classes
#-----------------------------------------------------


def add_info(kmer, copy_num):
    """"""
    GC_count = sum(kmer.count(x) for x in ['G', 'C'])
    # print GC_count
    # print(np.divide(GC_count, float(len(kmer))))
    GC_perc = round(np.divide(GC_count, float(len(kmer))), 4) * 100
    if any(x in kmer[-3:] for x in ['G', 'C']):
        GC_clamp = 'yes'
    else:
        GC_clamp = ''
    return(kmer, str(copy_num), str(GC_perc), GC_clamp)

# def filter_kmer(kmer, copy_num, GC_perc, GC_clamp):
#     """"""
#     GC_count = sum(kmer.count(x) for x in ['G', 'C'])
#     GC_perc = round(GC_count / len(kmer), 4) * 100
#     if any(x in kmer[-3:] for x in ['G', 'C']):
#         GC_clamp = 'yes'
#     else:
#         GC_clamp = ''
#     return(kmer, str(copy_num), str(GC_perc), GC_clamp)


#-----------------------------------------------------
# Step
# Process kmers
#-----------------------------------------------------

i = 0
for line in kmer_lines:
    line = line.rstrip()
    split_line = line.split()
    [kmer, copy_num, GC_perc, GC_clamp] = add_info(split_line[0], split_line[1])
    # out_line = "\t".join([kmer, copy_num, GC_perc, GC_clamp])
    # print out_line
    # if int(copy_num) > 1:
    #     if float(GC_perc) >= min_GC and float(GC_perc) <= max_GC:
    if (int(copy_num) > min_copy and
        float(GC_perc) >= min_GC and
        float(GC_perc) <= max_GC and
        GC_clamp == 'yes'
        ):
        i += 1
        # out_line = "\t".join([kmer, copy_num, GC_perc, GC_clamp])
        # out_list.append(out_line)
        # print out_line
        header = "_".join([">" + str(i), copy_num, GC_perc])
        print "\n".join([header, kmer])
