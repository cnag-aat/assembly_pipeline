from datetime import datetime
import os
import re

date = datetime.now().strftime('%Y%m%d.%H%M%S')

rule hypo:
    input:
      genome = "assembly.fasta",
      lr_bam = "assembly_minimap2.bam",
      sr_bam = "assembly_bwa.bam",
      reads_file = ["reads.1.fastq.gz", "reads.2.fastq.gz"],
      cov_script = "../scripts/get_cov.py"     
    output:
      polished ="assembly.hypo1.fasta"
    params:
      genome_size = "500000000",
      cov = 0,
      proc = 6,
      opts = " -B "
    wildcard_constraints:
      param="\d+"
    conda:
      "../envs/hypo1.0.3.yaml"
    threads: 6
    shell:
      "cd {wildcards.directory}hypo;"
      "ls -1 {input.reads_file} > {wildcards.directory}hypo/short_reads.list.txt;"
      "if [ {params.cov} == 0 ]; then coverage=$({input.cov_script} {input.genome} {input.sr_bam}); else coverage={params.cov}; fi;"
      "echo 'Illumina coverage is: '$coverage;"
      "echo 'Running: hypo -r @short_reads.list.txt -d {input.genome} -b {input.sr_bam} -c '$coverage' -s {params.genome_size}  -t {threads} -o {output.polished} -p {params.proc} {params.opts} {input.lr_bam}';"
      "hypo -r @short_reads.list.txt -d {input.genome} -b {input.sr_bam} -c $coverage -s {params.genome_size}  -t {threads} -o {output.polished} -p {params.proc} {params.opts} {input.lr_bam};"