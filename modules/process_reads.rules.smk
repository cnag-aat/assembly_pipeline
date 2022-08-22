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
    "ln -s $b.2_val_2.fq.gz {output.trim2};"
    "conda deactivate;"

rule concat_reads:
  input:
    fastqs = ["reads1.fastq.gz", "reads2.fastq.gz"]
  output:
    final_fastq = "all_reads.fastq.gz"
  log:
    "{dir}logs/" + str(date) + ".j%j.concat.{ext}.out",
    "{dir}logs/" + str(date) + ".j%j.concat.{ext}.err"
  benchmark:
    "{dir}logs/" + str(date) + ".concat.benchmark.{ext}.txt"
  conda:
    '../envs/ass_base.yaml'
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
    kmer = 21,
    genomescope_path = "/home/devel/jgomez/bin/genomescope2.0/genomescope2.0/",
    ploidy = 2,
  threads: 4
  shell:
    "source ~jgomez/init_shell.sh;"
    "conda activate {params.environment};"
    "meryl union-sum output {output.meryl_all} {input.input_run};"
    "d=`dirname {output.meryl_all}`;"
    "meryl histogram {output.meryl_all} > $d/meryl.hist;"
    "conda deactivate;"
    "module purge; module load gcc/6.3.0 R/3.6.0 mkl/12.1.6 PYTHON/3.6.0;"
    "export PATH={params.genomescope_path}:$PATH;"
    "export R_LIBS_USER=\"/home/devel/jgomez/R_libs\";"
    "genomescope.R -i $d/meryl.hist -o $d/out_k{params.kmer}/ -p {params.ploidy} -k {params.kmer};"

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
#    path = "/scratch/project/devel/aateam/bin",
    minlen = 1000,
    min_mean_q = 80,
    opts = ""
  threads: 8
  log:
    "{dir}logs/" + str(date) + ".j%j.filtlong.out",
    "{dir}logs/" + str(date) + ".j%j.filtlong.err"
  benchmark:
    "{dir}logs/" + str(date) + ".filtlong.benchmark.txt"
  conda:
    '../envs/filtlong0.2.1.yaml'
  shell:
 #   "module purge; module load PIGZ/2.3.3 gcc/4.9.3;"
    "filtlong --version;"
    "filtlong --min_length {params.minlen} --min_mean_q {params.min_mean_q} {params.opts} {input.reads} | pigz -p {threads} -c > {output.outreads};"

rule nanoplot:
  input:
    fastq = "reads.ont.fastq.gz"
  output:
    stats = "nanostats_out/NanoStats.txt",
    yield_len = "nanostats_out/Yield_By_Length.png",
    read_len = "nanostats_out/WeightedHistogramReadlength.png"
  params:
    outdir = "nanostats_out/",
  log:
    "{dir}logs/" + str(date) + ".NanoStats.out",
    "{dir}logs/" + str(date) + ".NanoStats.err"
  benchmark:
    "{dir}logs/" + str(date) + ".NanoStats.benchmark.txt"
  threads: 4
  conda:
    '../envs/nanoplot1.40.0'
  shell:
    "mkdir -p {params.outdir};" 
    "cd {params.outdir}; "
    "NanoPlot -t {threads} --plots dot --fastq {input.fastq} -o .; "

rule Kraken2:
  input:
    read = "reads.fastq.gz",
    database = "/scratch/groups/assembly/shared/databases/kraken2/bacteria_db/",
    kmers = "/scratch/groups/assembly/shared/databases/kraken2/bacteria_db/database100mers.kmer_distrib"
  output:
    report = "kraken2.report",
    readsout = "kraken2.seqs.out",
    abundance =  "bracken_abundance.txt",
  params:
    additional = "",
    prefix = "assembly"
  log:
    "logs/" + str(date) + ".kraken.out",
    "logs/" + str(date) + ".kraken.err",
  benchmark:
    "logs/" + str(date) + ".kraken.benchmark.txt",
  threads: 4
  conda:
    '../envs/kraken2.1.2.yaml'
  shell:
     "kraken2 --threads {threads} --db {input.database}  --use-names --report {output.report} {params.additional} {input.read} > {output.readsout}; "
     "est_abundance.py -i {output.report} -k {input.kmers} -l S -t 10 -o {output.abundance}; "