from datetime import datetime
import os
import re
import subprocess

module polish_workflow:
  snakefile: "../modules/polish.rules.smk"

#0. Main  
rr = config["Parameters"]["racon_rounds"]
pr = config["Parameters"]["pilon_rounds"]
mr = config["Parameters"]["medaka_rounds"]
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
      reads_file = [pe1_reads, pe2_reads],
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