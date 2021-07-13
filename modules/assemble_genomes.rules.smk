from datetime import datetime
import re
import os

date = datetime.now().strftime('%Y%m%d.%H%M%S')
shell.prefix("echo 'Cluster jobid $SLURM_JOB_ID';")
#keepfiles = False
#scripts_dir = os.path.dirname(sys.argv[0]) + "/../scripts/" 

if not os.path.exists("logs"):
  os.makedirs("logs")

rule flye:
  input:
    reads = os.getcwd() + "/ontreads.fastq.gz",
    env = "/home/devel/talioto/miniconda3/envs/flye-v2.8.3/"
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
  shell:
    "source  ~jgomez/init_shell.sh;"
    "conda activate {input.env};"
    "echo 'Running command: flye --{params.readtype} {input.reads} -o {params.outdir} -t {threads} -i {params.pol_iterations} {params.other_flye_opts}';"
    "flye --{params.readtype} {input.reads} -o {params.outdir} -t {threads} -i {params.pol_iterations} {params.other_flye_opts};"
    "ln -s {params.outdir}assembly.fasta {output.assembly};"
    "conda deactivate;"

rule nextdenovo:
  input:
    reads = os.getcwd() + "/ontreads.fastq.gz",
    config = "nextdenovo.cfg"
  output:
    assembly = os.getcwd() + "/nextDenovo/nextdenovo.assembly.fasta"
  params:
    outdir = os.getcwd() + "/nextDenovo/",
    module = "NEXTDENOVO/2.4.0"
  threads: 24
  log:
    "logs/" + str(date) + ".nextdenovo.out",
    "logs/" + str(date) + ".nextdenovo.err",
  shell:
    "module purge; module load {params.module};"
    "mkdir -p {params.outdir};"
    "ls {input.reads} > {params.outdir}long_reads.fofn;"
    "nextDenovo {input.config};"
    "ln -s {params.outdir}03.ctg_graph/nd.asm.fasta {output.assembly};"


