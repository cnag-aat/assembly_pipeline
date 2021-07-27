#!/usr/bin/env python3
import sys
import subprocess

fasta = sys.argv[1]
bam = sys.argv[2]

glen_out = subprocess.check_output("fastalength " + fasta + "|awk \'{len += $1} END{print len}\'", shell=True)
glen=str(glen_out).split("'")[1].replace('\\n','')
cov_out = subprocess.check_output("module purge; module load bedtools; bedtools genomecov -bg -ibam " + bam + " | awk \'{len=$3-$2; cov+=len*$4}END{print cov}\'", shell=True)
bcov = str(cov_out).split("'")[1].replace('\\n','')
gcov = int(bcov)/int(glen);
print (gcov)

