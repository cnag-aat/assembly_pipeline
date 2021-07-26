#!/usr/bin/env python3

import os
import sys
import re
import json

dir = sys.argv[1]
base = sys.argv[2]

with open (dir + "/stats/" + base + ".nseries.txt", 'r') as file:
  for line in file:
    if line.startswith('L_ass'):
      length = line.split('\t')[1]
    if line.startswith('num_blocks'):
      blocks = line.split('\t')[1]
    elif line.startswith('50'):
      n50 = line.split('\t')[1]
      l50 = line.split('\t')[2].rstrip()

with open (dir + "/busco/" + base + ".short_summary.txt", 'r') as file:
  for line in file:
    if re.search("C:", line):
      busco = line.rstrip().replace('\t', '') 

if os.path.exists(dir + "/merqury/" + base):
  with open (dir + "/merqury/" + base + "/completeness.stats", 'r') as file:
    for line in file:
      merq = line.split('\t')[4].rstrip()
  with open (dir + "/merqury/" + base + "/" + base + ".qv", 'r') as file:
    for line in file:
      qv = line.split('\t')[3]
else:
  merq = ""
  qv = ""

print (base + "\t" +  n50 + "\t" + l50 + "\t" + length + "\t" + blocks + "\t" + busco + "\t" + qv + "\t" + merq)
