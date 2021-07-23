from datetime import datetime
import re
import os

date = datetime.now().strftime('%Y%m%d.%H%M%S')
keepfiles = False

rule trim_galore:
  input:
    read1 = "illumina.1.fastq.gz",
    read2 = "illumina.2.fastq.gz",
  output:
    trim1 = "illumina.trimed.1.fastq.gz",
    trim2 = "illumina.trimed.2.fastq.gz",
  params:
    outdir = "illumina_trim",
    opts = "--gzip -q 20 --paired --retain_unpaired",
  threads: 4
  shell:
    "source ~jgomez/init_shell.sh;"
    "conda activate ~jgomez/conda_environments/preprocess_illumina;"
    "mkdir -p {params.outdir};"
    "cd {params.outdir}; "
    "trim_galore -j {threads} {params.opts} {input.read1} {input.read2} ;"
    "b=`basename {input.read1} .1.fastq.gz`;"
    "ln -s $b.1_val_1.fq.gz {output.trim1};"
    "ln -s $b.1_val_2.fq.gz {output.trim2};"
    "conda deactivate;"

rule concat_reads:
  input:
    fastqs = ["reads1.fastq.gz", "reads2.fastq.gz"]
  output:
    final_fastq = "all_reads.fastq.gz"
  threads: 2
  shell:
    "zcat {input.fastqs} | pigz -p {threads} -c  > {output.final_fastq};"

rule build_meryl_db:
  input:
    fastq = "all_reads.fastq.gz"
  output:
     out_dir = directory("meryl_db.meryl")
  params:
    environment = "~fcruz/.conda/envs/merqury_v1.1/",
    kmer = 21
  threads: 4
  shell:
    "source ~jgomez/init_shell.sh;"
    "conda activate {params.environment};"
    "meryl k={params.kmer} count output {output.out_dir} {input.fastq};"
    "conda deactivate;"

rule concat_meryl:
  input:
    input_run = "meryl_db"
  output:
    meryl_all = directory("meryl_db_all")
  params:
    environment = "~fcruz/.conda/envs/merqury_v1.1/;",
  threads: 4
  shell:
    "source ~jgomez/init_shell.sh;"
    "conda activate {params.environment};"
    "meryl union-sum output {output.meryl_all} {input.input_run};"
    "conda deactivate;"

rule long_ranger:
  input: 
    mkfastq_dir = "input_dir/"
  output:
    fastq_out = "preprocess_10X_linkedreads/10X.barcoded.fastq.gz",
    sum_out = "preprocess_10X_linkedreads/10X.summary.csv"
  params:
    path = "/scratch/project/devel/aateam/src/10X/longranger-2.2.2",
    outdir = "preprocess_10X_linkedreads/",
    sample = "raw_10X"
  threads: 8
  run:
    shell(
      "export PATH={params.path}:$PATH;"
      "mkdir -p {params.outdir};"
      "cd {params.outdir};"
      "longranger basic --id={params.sample} --sample={params.sample} --fastqs={input.mkfastq_dir} --localcores={threads};"
      "cp {params.sample}/outs/barcoded.fastq.gz {output.fastq_out};"
      "cp {params.sample}/outs/summary.csv {output.sum_out};"
    )
    if keepfiles == False:
      shell(
        "echo 'Removing longranger run dir:{params.outdir}{params.sample}';"
        "module load bsc;"
        "lrm -o {params.outdir}{params.sample}.rmlist.sh {params.outdir}{params.sample}; sh {params.outdir}{params.sample}.rmlist.sh;"
      )

rule filtlong:
  input:
    reads = "ont_reads.fastq.gz"
  output:
    outreads = "ont_reads.filtlong.fastq.gz"
  params:
    path = "/scratch/project/devel/aateam/bin/filtlong",
    minlen = 1000,
    min_mean_q = 80,
    opts = ""
  threads: 8
  shell:
    "module purge; module load PIGZ/2.3.3 gcc/4.9.3;"
    "export PATH={params.path}:$PATH;"
    "filtlong --version;"
    "filtlong --min_length {params.minlen} --min_mean_q {params.min_mean_q} {params.opts} {input.reads} | pigz -p {threads} -c > {output.outreads};"
