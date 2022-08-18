from datetime import datetime
import re
import os

date = datetime.now().strftime('%Y%m%d.%H%M%S')

if not os.path.exists("logs"):
  os.makedirs("logs")

rule flye:
  input:
    reads = os.getcwd() + "/ontreads.fastq.gz",
  output:
    assembly = os.getcwd() + "/flye/flye.assembly.fasta",
  params:
    outdir = os.getcwd() + "/flye/",
    readtype = "nano-raw",
    pol_iterations = 2,
    other_flye_opts = ""  #include here genome size in pipeline
  threads: 24
  log:
    "logs/" + str(date) + ".flye.out",
    "logs/" + str(date) + ".flye.err",
  benchmark:
    "logs/" + str(date) + ".flye.benchmark.txt"
  conda:
    '../envs/flye2.9.1.yaml'
  shell:
    "mkdir -p {params.outdir}out;"
    "cd {params.outdir};"
    "echo 'Running command: flye --{params.readtype} {input.reads} -o {params.outdir}out -t {threads} -i {params.pol_iterations} {params.other_flye_opts}';"
    "flye --{params.readtype} {input.reads} -o {params.outdir}out -t {threads} -i {params.pol_iterations} {params.other_flye_opts};"
    "ln -s {params.outdir}out/assembly.fasta {output.assembly};"

rule nextdenovo:
  input:
    reads = os.getcwd() + "/ontreads.fastq.gz"
  output:
    assembly = os.getcwd() + "/nextDenovo/nextdenovo.assembly.fasta"
  params:
    outdir = os.getcwd() + "/nextDenovo/",
    config = "nextdenovo.cfg"
  threads: 24
  log:
    "logs/" + str(date) + ".nextdenovo.out",
    "logs/" + str(date) + ".nextdenovo.err",
  benchmark:
    "logs/" + str(date) + ".nextdenovo.benchmark.txt"
  envmodules:
    "NextDenovo/2.5.0"
  shell:
    "mkdir -p {params.outdir};"
    "ls {input.reads} > {params.outdir}long_reads.fofn;"
    "nextDenovo {params.config};"
    "ln -s {params.outdir}03.ctg_graph/nd.asm.fasta {output.assembly};"


