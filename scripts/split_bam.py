#!/usr/bin/env python3

import os
import sys
import re
import json

lengths = sys.argv[1]
map_all = sys.argv[2]
chunks = sys.argv[3]
threads = sys.argv[4]

n=0

def write_bed(bed):
  bed.write(line.rstrip().split(' ')[1] + "\t0\t" + line.split(' ')[0] + "\n")

def split_bam(bed, bamchunk, bam_all, threads):
   os.system("samtools view -@ " + threads + " " + bam_all + " -L " + bed + " -Sb -o " + bamchunk)
   os.system("samtools index " + bamchunk)

with open (lengths, 'r') as file:
  for line in file:
    if n < int(chunks) -1:
      bed = open("chunk" + str(n) + ".bed", "w")
      write_bed(bed)
      bed.close()
      split_bam("chunk" + str(n) + ".bed","chunk" + str(n) + ".bam", map_all, threads)
      n+=1
    elif n == int(chunks)-1:
      bed = open("chunk" + str(n) + ".bed", "a")
      write_bed(bed)
      bed.close()
  split_bam("chunk" + str(n) + ".bed","chunk" + str(n) + ".bam", map_all, threads)

