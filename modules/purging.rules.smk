from datetime import datetime
import re
import os

rule purge_dups:
  input:
    assembly_in = "assembly.fa",
    reads = os.getcwd() + "/ontreads.fastq.gz",
    mapping = os.getcwd() + "/mappings/ontreads.paf.gz",
  output:
    assembly_out = "assembly.purged.fa",
    plot = "PB.cov.png"
  params:
    scripts_dir = "",
    base = "assembly",
    dir = os.getcwd(),
    calcuts_opts = ""
  threads: 12
  conda:
    "../envs/purge_dups1.2.6.yaml"
  shell:
    "cd {params.dir};"
    "if [ -f {params.dir}/{params.base}.split.self.paf.gz ]; then " \
      "echo '{params.dir}/{params.base}.split.self.paf.gz already exists, skipping alignment';" \
    "else " \
      "split_fa {input.assembly_in} > {params.base}.split;" \
      "minimap2 -t {threads} -xasm5 -DP {params.base}.split {params.base}.split | pigz -p {threads} -c > {params.base}.split.self.paf.gz;" \
      "pbcstat {input.mapping}; fi;" \
    "calcuts {params.calcuts_opts} PB.stat > cutoffs 2>calcults.log;"
    "purge_dups -2 -T cutoffs -c PB.base.cov {params.base}.split.self.paf.gz > dups.bed 2> purge_dups.log;"
    "get_seqs -e dups.bed {input.assembly_in};"
    "ln -s purged.fa {output.assembly_out};"
    "python3 {params.scripts_dir}hist_plot.py -c cutoffs PB.stat {output.plot};"
    