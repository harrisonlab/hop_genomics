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
from Bio import SeqIO

ap = argparse.ArgumentParser()

ap.add_argument('--assembly', required=True,type=str,help='fasta file of assembly')
ap.add_argument('--kmers', required=True,type=str,help='kmer tsv')
ap.add_argument('--prefix', required=True,type=str,help='path and prefix for output files')
conf = ap.parse_args()


# with open(conf.kmers) as f:
#     kmer_lines = f.readlines()
kmer_records = list(SeqIO.parse(conf.kmers, "fasta"))

seq_records = list(SeqIO.parse(conf.assembly, "fasta"))

out = conf.prefix

#-----------------------------------------------------
# Step 2
# Define classes
#-----------------------------------------------------


class SeqObj(object):
    """
    Object containing location of high abundance kmers and their vicinity
    one another.
    Attributes:
        seq_name: sequence data
        kmer_locations: dictionary of kmer locations of kmer start sites.
        kmer_pairs:
    """
    def __init__(self, id, seq):
        """set the object, defining the dictionary."""
        dict = defaultdict(list)
        self.id = id
        self.seq = seq
        self.kmer_loc_dict = defaultdict(str)
        self.kmer_pair_dict = defaultdict(list)
        self.F_hit_dict = defaultdict(str)
        self.R_hit_dict = defaultdict(str)
    def locate_kmer(self, k_rec):
        """ Create the conversion dictionary from input lines """
        k_id = k_rec.id
        k_seq = str(k_rec.seq)
        # k_rec_rc = k_rec.reverse_complement
        # print str(k_rec_rc)
        k_rec_rc = k_rec.reverse_complement()
        k_seq_rc = str(k_rec_rc.seq)
        # k_seq_rc = str(k_rec_rc.seq)
        # print k_seq_rc
        # Identify forward matches:
        iter = re.finditer(k_seq, self.seq)
        indices = [m.start(0) for m in iter]
        for start_loc in indices:
            self.F_hit_dict[start_loc] = k_id
        # Identify revcomp matches:
        iter = re.finditer(k_seq_rc, self.seq)
        indices = [m.end(0) for m in iter]
        for start_loc in indices:
            self.R_hit_dict[start_loc] = k_id
    def get_locations(self):
        """ Create the conversion dictionary from input lines """
        outlines = []
        outlines.append(self.id + "\t" + ",".join((str(x) for x in self.F_hit_dict.keys())))
        outlines.append(self.id + "_rev)\t" + ",".join((str(x) for x in self.R_hit_dict.keys())))
        return(outlines)
        # return([",".join(self.F_hit_dict.keys()), ",".join(self.R_hit_dict.keys())])
        # return(self.F_hit_dict.keys())
    def identify_pairs_old(self):
        outlines = []
        F_hits = self.F_hit_dict.keys()
        R_hits = self.R_hit_dict.keys()
        for F in sorted(F_hits):
            for R in sorted(R_hits):
                dist = int(R) - int(F)
                if (dist >= 125 and dist <= 200):
                    amplicon = self.seq[F:R]
                    key = "-".join(sorted([self.F_hit_dict[F], self.R_hit_dict[R]]))
                    self.kmer_pair_dict[key].append("\t".join([self.id, str(F), str(R), self.F_hit_dict[F], self.R_hit_dict[R]]))
                    outlines.append("\t".join([key, self.id, str(F), str(R), self.F_hit_dict[F], self.R_hit_dict[R], amplicon]))
        return(outlines)
    def identify_pairs(self):
        outlines = []
        F_hits = self.F_hit_dict.keys()
        # R_hits = self.R_hit_dict.keys()
        for F in sorted(F_hits):
            for i in range(125, 200):
                i += 40
                i += F
                if self.R_hit_dict[i]:
                    R = i
                    amplicon = self.seq[F:R]
                    key = "-".join(sorted([self.F_hit_dict[F], self.R_hit_dict[R]]))
                    self.kmer_pair_dict[key].append("\t".join([self.id, str(F), str(R), self.F_hit_dict[F], self.R_hit_dict[R]]))
                    outlines.append("\t".join([key, self.id, str(F), str(R), self.F_hit_dict[F], self.R_hit_dict[R], amplicon]))
        return(outlines)


        # print(indices)
    #     for line in lines:
    #         line = line.rstrip()
    #         split_line = line.split("\t")
    #         old_gene_id = split_line[0]
    #         new_gene_id = split_line[2]
    #         conv_dict = self.conversion_dict
    #         conv_dict[old_gene_id] = new_gene_id
    #         self.conversion_dict = conv_dict

# def kmer2gff(kmer, seq_records):
#     """"""
#     for record in seq_records:


seq_dict = defaultdict(list)

for record in seq_records:
    record = record.upper()
    seq_id = record.id
    seq = str(record.seq)
    seq_dict[seq_id] = SeqObj(seq_id, seq)
    for k_rec in kmer_records:
        # print k_rec.seq
        # print k_rec.reverse_complement().seq
        seq_dict[seq_id].locate_kmer(k_rec)
    # k_id = '2_2428_60.0'
    # k_seq = 'GTGCTCGCGGCTTTAGTGCA'
    # seq_dict[seq_id].locate_kmer(k_id, k_seq)
    f = open("_".join([out, seq_id, "kmer_locations.tsv"]),"w+")
    f.write("\n".join(seq_dict[seq_id].get_locations()))
    f.close()

    f = open("_".join([out, seq_id, "kmer_pairs.tsv"]),"w+")
    f.write("\n".join(seq_dict[seq_id].identify_pairs()))
    f.close()

    # seq_dict[seq_id].locate_kmer(kmer)
    # if seq_id == "000000F":
    #     print "\n".join([">"+record.id, str(record.seq)])
    # if seq_id == "000001F":
    #     print "\n".join([">"+record.id, str(record.seq)])
