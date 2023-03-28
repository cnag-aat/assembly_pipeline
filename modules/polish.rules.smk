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

rule nextpolish_lr:
  input:
    genome = "assembly.fasta",
    bam =  "assembly_minimap2.bam"
  output:
    polished = "assembly.nextpolish_ont1.fasta"
  params:
    lrtype = "ont",
    path = "NEXTPOLISH/v1.4.0/NextPolish"
  envmodules:
    "NextPolish/1.4.1-GCC-11.2.0"
  threads: 12
  shell:
    "cd {wildcards.directory}nextpolish;"
    "echo {input.bam} > lgs.fofn;"
    "python {params.path}/lib/nextpolish2.py -g {input.genome} -p {threads} -l lgs.fofn -r {params.lrtype} > {output.polished};"

rule nextpolish_sr:
  input: 
    genome = "assembly.fasta",
    bam =  "assembly_bwa.bam"
  output:
    polished = "assembly.nextpolish_ont2.nextpolish_ill1.fasta",
  params:
    task = 1,
    path = "NEXTPOLISH/v1.4.0/NextPolish"
  envmodules:
    "NextPolish/1.4.1-GCC-11.2.0"
  threads: 12
  shell:
    "cd {wildcards.directory}nextpolish;"
    "samtools faidx {input.genome};"
    "python {params.path}/lib/nextpolish1.py -g {input.genome}  -p {threads} -s {input.bam} -t {params.task} > {output.polished};"
