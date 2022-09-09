#!/usr/bin/env python3

import os
import sys
import re
import json
import ast
import argparse
 
parser = argparse.ArgumentParser(description="Get final table for assembly pipeline",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-s", "--stats", nargs="+",  help="List with stats files from different assemblies")
parser.add_argument("-b", "--buscos", nargs="+",  help="List with busco short summaries from different assemblies")
parser.add_argument("-m", "--merqs", nargs="+",  help="List with merqury results from different assemblies")

args = parser.parse_args()
print("assembly\tcN50\tcL50\tsN50\tsL50\ttotal_len\ttotal_seq\tBUSCOv5\tQV\tMerqury_completeness\tFalse_duplications")
buscos_dict = {}
merqs_dict = {}
for file in args.buscos:
  base_elements = os.path.basename(file).split('.')[:-2]
  base = ".".join(base_elements)
  buscos_dict[base] = file
for file in args.merqs:
  base = os.path.basename(os.path.dirname(file))
  merqs_dict[base] = os.path.dirname(file)

assemblies={}
for file in args.stats:
  base_elements = os.path.basename(file).split('.')[:-2]
  base = ".".join(base_elements)
  cn50 = cl50 = sl50 = sn50 = len = seqs = 0
  if base not in assemblies:
    with open (file, 'r') as file:
      for line in file:
        if line.startswith('Contig N50'):
          cn50 = line.split('\t')[-1].rstrip()
        if line.startswith('Contig L50'):
          cl50 = line.split('\t')[-1].rstrip()
        if line.startswith('Scaffold N50'):
          sn50 = line.split('\t')[-1].rstrip()
        if line.startswith('Scaffold L50'):
          sl50 = line.split('\t')[-1].rstrip()
        if line.startswith('Scaffold num_bp'):
          len = line.split('\t')[-1].rstrip()
        if line.startswith('Scaffold num_seq'):
          seqs = line.split('\t')[-1].rstrip()
    with open (buscos_dict[base], 'r') as busco_file:
      for line in busco_file:
        if line.startswith('	C:'):
          busco = line.split('\s+')[-1].rstrip()
    merqdir = merqs_dict[base]
    with open (merqdir + "/" + base + ".qv", 'r') as merq_file:
      for line in merq_file:
        qv = line.split('\t')[3]
    with open (merqdir + "/" + base + ".completeness.stats", 'r') as merq_file:
      for line in merq_file:
        completeness = line.split('\t')[4].rstrip()
    with open (merqdir + "/" + base + ".false_duplications.txt", 'r') as merq_file:
      for line in merq_file:
        dup = line.split('\t')[9].rstrip()
    assemblies[base] = ""
    print(base,cn50,cl50,sn50,sl50,len,seqs,busco, qv, completeness, dup, sep="\t")