#!/usr/bin/env python3

import os
import sys
import re
import json

dir = sys.argv[1]
chunks = sys.argv[2]

for n in range(0, int(chunks)+1):
  seqs = []
  cont = 0
#  print (n)
  with open (dir + "chunk" + str(n) + ".bed", 'r') as bedfile:
    for line in bedfile:
      seq = line.rstrip().split('\t')[0] 
      seqs.append(seq + "_pilon")
  with open (dir + "chunk" + str(n) + ".polished.fasta", 'r') as fastafile:
    for line in fastafile:
      if line.startswith('>'):
        seqid = line.rstrip().replace('>', '')
        if seqid in seqs:
          cont = 1
          print (line.rstrip())
        else: 
          cont = 0
      elif cont == 1:
        print(line.rstrip())
