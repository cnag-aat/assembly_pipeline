from datetime import datetime
import os
import re
import subprocess

module polish_workflow:
  snakefile: "../modules/polish.rules.smk"

#0. Main  
nor = config["Parameters"]["nextpolish_ont_rounds"]
nir = config["Parameters"]["nextpolish_ill_rounds"]
hr = config["Parameters"]["hypo_rounds"]
keepfiles = config["Parameters"]["keep_intermediate"]
base = config["Parameters"]["base_name"]

##1. Define path for files and variables
make_logs = []
inputs = config["Inputs"]["Assemblies for polishing"]
polished = []
bwa = {}
minimap2 = {}
hypo_in = {}
hypo_lrmap = {}
hypo_srmap = {}
nor_in = {}
nor_map={}
nir_in = {}
nir_map={}
tasks={}
for assembly_in in inputs:
  basename = os.path.splitext(os.path.basename(assembly_in))[0]
  rundir = inputs[assembly_in]
  mappings_dir = rundir + "mappings/"
  if not os.path.exists(mappings_dir):
    os.makedirs(mappings_dir)
  make_logs.append(rundir)
  if hr > 0:
    hypo_dir = rundir + "hypo/"
    make_logs.append(hypo_dir)
    bwa[basename] = assembly_in
    hypo_in[hypo_dir + basename + ".hypo1.fasta"] = assembly_in
    hypo_srmap[hypo_dir + basename + ".hypo1.fasta"] = mappings_dir + basename + "_bwa.bam"

    if config["Hypo"]["long_reads"]:
      minimap2[basename] = assembly_in
      hypo_lrmap[hypo_dir + basename + ".hypo1.fasta"] = mappings_dir + basename + "_minimap2.bam"
    
    for i in range(1,hr):
      file_name = basename + ".hypo" + str(i) 
      bwa[file_name] = hypo_dir + file_name + ".fasta"    
      n=i+1
      hypo_in[hypo_dir + basename + ".hypo" + str(n) + ".fasta"] = hypo_dir + file_name + ".fasta"     
      hypo_srmap[hypo_dir + basename + ".hypo" + str(n) + ".fasta"] = mappings_dir + file_name + "_bwa.bam"

      if config["Hypo"]["long_reads"]:
        minimap2[file_name] = hypo_dir + file_name + ".fasta"
        hypo_lrmap[hypo_dir + basename + ".hypo" + str(n) + ".fasta"] = mappings_dir + file_name + "_minimap2.bam"
    for i in hypo_in:
      polished.append(i)
  
  if nir > 0 or nor > 0:
    nextpolish_dir = rundir + "nextpolish/"
    make_logs.append(nextpolish_dir)

  if nor > 0:
    minimap2[basename] = assembly_in
    nor_in[nextpolish_dir + basename + ".nextpolish_ont1.fasta"] = assembly_in
    nor_map[nextpolish_dir + basename + ".nextpolish_ont1.fasta"] = mappings_dir + basename + "_minimap2.bam"

    for i in range(1,nor):
      file_name = basename + ".nextpolish_ont" + str(i) 
      minimap2[file_name] = nextpolish_dir + file_name + ".fasta"    
      n=i+1
      nor_in[nextpolish_dir + basename + ".nextpolish_ont" + str(n) + ".fasta"] = nextpolish_dir + file_name + ".fasta"     
      nor_map[nextpolish_dir + basename + ".nextpolish_ont" + str(n) + ".fasta"] = mappings_dir + file_name + "_minimap2.bam"
    for i in nor_in:
      polished.append(i)

  if nir >0:
    nir_base = basename
    if nor >0:
      nir_base += ".nextpolish_ont" + str(nor)
      nir_in[nextpolish_dir + nir_base + ".nextpolish_ill1_tmp.fa"] = nextpolish_dir + nir_base + ".fasta"
    else:
      nir_in[nextpolish_dir + nir_base + ".nextpolish_ill1_tmp.fa"] = assembly_in

    bwa[nir_base] = nextpolish_dir + nir_base + ".fasta"
    bwa[nir_base + ".nextpolish_ill1_tmp"] = nextpolish_dir + nir_base + ".nextpolish_ill1_tmp.fa"
    
    nir_in[nextpolish_dir + nir_base + ".nextpolish_ill1.fasta"] = nextpolish_dir + nir_base + ".nextpolish_ill1_tmp.fa"
    nir_map[nextpolish_dir +  nir_base + ".nextpolish_ill1_tmp.fa"] = mappings_dir + nir_base + "_bwa.bam"
    nir_map[nextpolish_dir + nir_base + ".nextpolish_ill1.fasta"] = mappings_dir +  nir_base + ".nextpolish_ill1_tmp_bwa.bam" 
    tasks[nextpolish_dir + nir_base + ".nextpolish_ill1_tmp.fa"] = 1
    tasks[nextpolish_dir +  nir_base + ".nextpolish_ill1.fasta"] = 2
    polished.append(nextpolish_dir + nir_base + ".nextpolish_ill1.fasta")
    for i in (range(1,nir)):
      pref = nir_base + ".nextpolish_ill"
      n = i+1
      bwa[pref + str(n) + "_tmp"] = nextpolish_dir + pref + str(n) + "_tmp.fa"       
      bwa[pref + str(i)] = nextpolish_dir + pref + str(i) + ".fasta"
      polished.append(nextpolish_dir + pref + str(i) + ".fasta")
      nir_in[nextpolish_dir + pref + str(n) + "_tmp.fa"] = nextpolish_dir + pref + str(i) + ".fasta"
      nir_in[nextpolish_dir + pref + str(n) + ".fasta"] = nextpolish_dir + pref + str(n) + "_tmp.fa"
      nir_map[nextpolish_dir + pref + str(n) + "_tmp.fa"] = mappings_dir + pref + str(i) + "_bwa.bam"
      nir_map[nextpolish_dir + pref + str(n) + ".fasta"] = mappings_dir + pref + str(n) + "_tmp_bwa.bam"
      tasks[nextpolish_dir + pref + str(n) + "_tmp.fa"] = 1
      tasks[nextpolish_dir + pref + str(n) + ".fasta"] = 2

for d in make_logs:
  if not os.path.exists(d + "logs"):
    os.makedirs(d+"logs")

#Polish genomes
if hr > 0:
  additional_hypo_opts = ""
 
  if config["Hypo"]["options"] != None:
    additional_hypo_opts = config["Hypo"]["options"]

  use rule hypo from polish_workflow with:
    input:
      genome = lambda wildcards: hypo_in[wildcards.directory + "hypo/" + wildcards.base + ".hypo" + wildcards.param + ".fasta"],
      lr_bam = lambda wildcards: expand(hypo_lrmap[wildcards.directory + "hypo/" + wildcards.base + ".hypo" + wildcards.param + ".fasta"]) \
               if config["Hypo"]["long_reads"] else [],
      sr_bam = lambda wildcards: hypo_srmap[wildcards.directory + "hypo/" + wildcards.base + ".hypo" + wildcards.param + ".fasta"],
      reads_file = [pe1_reads, pe2_reads] if r10X_reads == None else [r10X_reads],
      cov_script = scripts_dir + "get_cov.py"       
    output:
      polished ="{directory}hypo/{base}.hypo{param}.fasta"
    params:
      genome_size = config["Parameters"]["genome_size"],
      cov = config["Hypo"]["illumina coverage"],
      proc = config["Hypo"]["processes"],
      opts = lambda wildcards: additional_hypo_opts + " -B " \
                               if config["Hypo"]["long_reads"] else additional_hypo_opts 
    wildcard_constraints:
      param="\d+"
    log:
      "{directory}hypo/logs/" + str(date) + ".j%j.{base}.hypo{param}.out",
      "{directory}hypo/logs/" + str(date) + ".j%j.{base}.hypo{param}.err",
    benchmark:
      "{directory}hypo/logs/" + str(date) + ".{base}.hypo{param}.benchmark.txt"
    threads: config["Parameters"]["hypo_cores"]

if nor > 0:
  nextpolish_lrtype= ""
  if re.search("nano",config["Parameters"]["lr_type"],):
    nextpolish_lrtype = "ont"
  elif re.match("pacbio-corr", config["Parameters"]["lr_type"]):
    nextpolish_lrtype = "hifi"
  elif re.match("pacbio-raw", config["Parameters"]["lr_type"], ):
    nextpolish_lrtype = "clr"
  
  use rule nextpolish_lr from polish_workflow with:
    input:
      genome = lambda wildcards: nor_in[wildcards.directory + "nextpolish/" + wildcards.base + ".nextpolish_ont" + wildcards.param + ".fasta"],
      bam =   lambda wildcards: nor_map[wildcards.directory + "nextpolish/" + wildcards.base + ".nextpolish_ont" + wildcards.param + ".fasta"]
    output:
      polished = "{directory}nextpolish/{base}.nextpolish_ont{param}.fasta"
    params:
      lrtype = nextpolish_lrtype,
      path = "/software/assembly/easybuild/software/NextPolish/1.4.1-GCC-11.2.0"
    wildcard_constraints:
      param="\d+"
    log:
      "{directory}nextpolish/logs/" + str(date) + ".j%j.{base}.nextpolish_lr{param}.out",
      "{directory}nextpolish/logs/" + str(date) + ".j%j.{base}.nextpolish_lr{param}.err",
    benchmark:
      "{directory}nextpolish/logs/" + str(date) + ".{base}.nextpolish_lr{param}.benchmark.txt",
    threads:  config["Parameters"]["nextpolish_cores"]

use rule nextpolish_sr from polish_workflow with:
  input: 
    genome = lambda wildcards: nir_in[wildcards.directory + "/nextpolish/" + wildcards.base + ".nextpolish_ill" + wildcards.param ],
    bam =   lambda wildcards: nir_map[wildcards.directory + "/nextpolish/" + wildcards.base + ".nextpolish_ill" + wildcards.param ]
  output:
    polished = "{directory}/nextpolish/{base}.nextpolish_ill{param}"
  params:
    task = lambda wildcards: tasks[wildcards.directory + "/nextpolish/" + wildcards.base + ".nextpolish_ill" + wildcards.param],
    path = "/software/assembly/easybuild/software/NextPolish/1.4.1-GCC-11.2.0"
  log:
    "{directory}/nextpolish/logs/" + str(date) + ".j%j.{base}.nextpolish_sr{param}.out",
    "{directory}/nextpolish/logs/" + str(date) + ".j%j.{base}.nextpolish_sr{param}.err",
  benchmark:
    "{directory}/nextpolish/logs/" + str(date) + ".{base}.nextpolish_sr{param}.benchmark.txt",
  threads:  config["Parameters"]["nextpolish_cores"]