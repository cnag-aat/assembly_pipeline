#!/usr/bin/env python3
import sys
import subprocess
import os

##Script to estimate the illumina coverage before running hypo.
##Author: Jèssica Gómez Garrido
##email: jessica.gomez@cnag.crg.eu


fasta = sys.argv[1]
bam = sys.argv[2]

dir=os.path.dirname(os.path.realpath(__file__))
glen_out = subprocess.check_output(dir + "/fastalength " + fasta + "|awk \'{len += $1} END{print len}\'", shell=True)
glen=str(glen_out).split("'")[1].replace('\\n','')
cov_out = subprocess.check_output("bedtools genomecov -bg -ibam " + bam + " | awk \'{len=$3-$2; cov+=len*$4}END{print cov}\'", shell=True)
bcov = str(cov_out).split("'")[1].replace('\\n','')
gcov = int(bcov)/int(glen);
print (str(gcov).rstrip())

