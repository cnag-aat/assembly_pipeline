from datetime import datetime
import re
import os

date = datetime.now().strftime('%Y%m%d.%H%M%S')

if not os.path.exists("logs"):
  os.makedirs("logs")

rule purge_dups:
  input:
    assembly_in = "assembly.fa",
    reads = os.getcwd() + "/ontreads.fastq.gz",
    mapping = os.getcwd() + "/mappings/ontreads.paf.gz",
  output:
    assembly_out = "assembly.purged.fa"
  params:
    module = "PURGEDUPS/1.2.5",
    base = "assembly",
    dir = os.getcwd()
  threads: 12
  log:
    "logs/" + str(date) + ".purge_dups.out",
    "logs/" + str(date) + ".purge_dups.err",
  shell:
    "module purge; module load {params.module} PIGZ/2.3.3;"
    "cd {params.dir};"
    "split_fa {input.assembly_in} > {params.base}.split;"
    "minimap2 -t {threads} -xasm5 -DP {params.base}.split {params.base}.split | pigz -p {threads} -c > {params.base}.split.self.paf.gz;"
    "pbcstat {input.mapping};"
    "calcuts PB.stat > cutoffs 2>calcults.log;"
    "purge_dups -2 -T cutoffs -c PB.base.cov {params.base}.split.self.paf.gz > dups.bed 2> purge_dups.log;"
    "get_seqs -e dups.bed {input.assembly_in};"
    "ln -s purged.fa {output.assembly_out};"
    "module purge; module load gcc/6.3.0 PYTHON/3.7.1;"
    "python3 /apps/{params.module}/scripts/hist_plot.py -c cutoffs PB.stat PB.cov.png;"
