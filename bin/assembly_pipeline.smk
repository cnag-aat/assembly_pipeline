from datetime import datetime
import os
import re
import subprocess

date = datetime.now().strftime('%Y%m%d.%H%M%S')

module assembly_workflow:
  snakefile: "../modules/assemble_genomes.rules.smk"
module eval_workflow:
  snakefile: "../modules/evaluate_assemblies.rules.smk"
module preprocess_workflow:
  snakefile: "../modules/process_reads.rules.smk"

working_dir = config["Outputs"]["base_dir"]
scripts_dir = config["Inputs"]["scripts_dir"]
shell.prefix("echo 'Cluster jobid $SLURM_JOB_ID'; export PATH=" + scripts_dir + ":$PATH;")
logs_dir = config["Parameters"]["logs_dir"]
if not os.path.exists(logs_dir):
  os.makedirs(logs_dir)

keepfiles = config["Parameters"]["keep_intermediate"]
base = config["Parameters"]["base_name"]

##0. Define path for files and variables

flye_dir = config["Outputs"]["flye_dir"]
nextdenovo_dir = config["Outputs"]["nextdenovo_dir"]
flye_assembly = config["Outputs"]["flye_out"]
nextdenovo_assembly = config["Outputs"]["nextdenovo_out"]
targets = []
if config["Parameters"]["run_flye"] == True:
  targets.append(flye_assembly)
  if not os.path.exists(flye_dir + "logs"):
    os.makedirs(flye_dir + "logs")
if config["Parameters"]["run_nextdenovo"] == True:
  targets.append(nextdenovo_assembly)
  if not os.path.exists(nextdenovo_dir + "logs"):
    os.makedirs(nextdenovo_dir + "logs")

#print(targets)
#1- Define rule all
rule all:
  input:
    targets,
    config["Outputs"]["stats_out"]
  log:
    logs_dir + str(date) + ".rule_all.out",
    logs_dir + str(date) + ".rule_all.err"

##2- Obtain input reads
fastqs = {}
if config["Parameters"]["run_flye"] == True or config["Parameters"]["run_nextdenovo"] == True:
  ONT_filtered = config["Inputs"]["ONT_filtered"]
  if not os.path.exists(ONT_filtered):
    if config["Inputs"]["ONT_reads"] == None:
      ont_dir = config["Inputs"]["ONT_dir"]
      ont_list = config["Wildcards"]["ONT_wildcards"].split(',')
      ont_reads = working_dir + "reads.ont.fastq.gz"
      extensions = ["fastq.gz"]
      for i in extensions:
        fastqs["ont."+i] = []
        for file in ont_list:
          fastqs["ont."+i].append(ont_dir + file + "." + i)
    else:
      ont_reads = config["Inputs"]["ONT_reads"] 

if len(fastqs) > 0:
 use rule concat_reads from preprocess_workflow with:
  input:
    fastqs = lambda wildcards: fastqs[wildcards.ext]
  output:
    final_fastq = working_dir + "reads.{ext}"
  log:
    logs_dir + str(date) + ".concat.{ext}.out",
    logs_dir + str(date) + ".concat.{ext}.err"
  threads: config["Parameters"]["concat_cores"]  

if not os.path.exists(ONT_filtered):
  use rule filtlong from preprocess_workflow with:
    input:
      reads = ont_reads
    output:
      outreads = ONT_filtered
    params:
      path = config["Filtlong"]["Filtlong path"],
      minlen = config["Filtlong"]["Filtlong minlen"],
      min_mean_q = config["Filtlong"]["Filtlong min_mean_q"]
    log:
      logs_dir + str(date) + ".filtlong.out",
      logs_dir + str(date) + ".filtlong.err"
    threads: config["Parameters"]["concat_cores"]  

##Run assemblers
if config["Parameters"]["run_flye"] == True:
  use rule flye from assembly_workflow with:
    input:
      reads = ONT_filtered,
      env = config["Flye"]["Flye environment"]
    output:
      assembly = flye_assembly
    params:
      outdir = flye_dir,
      readtype = config["Parameters"]["lr_type"],
      pol_iterations = config["Flye"]["Flye polishing iterations"],
      other_flye_opts = config["Flye"]["options"]
    log:
      flye_dir + "logs/" + str(date) + ".flye.out",
      flye_dir + "logs/" + str(date) + ".flye.err"
    threads: config["Flye"]["Flye cores"]

if config["Parameters"]["run_nextdenovo"] == True:
 use rule nextdenovo from assembly_workflow with:
  input:
    reads = ONT_filtered,
    config = config["Parameters"]["ndconfFile"]
  output:
    assembly = nextdenovo_assembly
  params:
    outdir = nextdenovo_dir,
    module = config["Nextdenovo"]["Nextdenovo module"],
  log:
    nextdenovo_dir + "logs/" + str(date) + ".nextdenovo.out",
    nextdenovo_dir + "logs/" + str(date) + ".nextdenovo.err"
  threads: config["Nextdenovo"]["Nextdenovo cores"]

##Run polishing pipeline
include: "../modules/polish_assemblies.v03.smk"
