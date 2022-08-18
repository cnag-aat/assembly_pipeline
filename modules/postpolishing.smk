from datetime import datetime
import os
import re
import subprocess

date = datetime.now().strftime('%Y%m%d.%H%M%S')

module purging_workflow:
  snakefile: "../modules/purging.rules.smk"
module eval_workflow:
  snakefile: "../modules/evaluate_assemblies.rules.smk"
module preprocess_workflow:
  snakefile: "../modules/process_reads.rules.smk"

working_dir = config["Outputs"]["base_dir"]
eval_dir = config["Outputs"]["eval_dir"]
scripts_dir = config["Inputs"]["scripts_dir"]
logs_dir = working_dir + "logs"
if not os.path.exists(logs_dir):
  os.makedirs(logs_dir)

shell.prefix( "export PATH=" + scripts_dir + ":$PATH;")
keepfiles = config["Parameters"]["keep_intermediate"]
base = config["Parameters"]["base_name"]

##0. Define path for files and variables
targets_purgedups = []
input_assemblies = {}
if config["Parameters"]["run_purgedups"] == True:
  if config["Inputs"]["ONT_reads"] != None:
    ont_reads = config["Inputs"]["ONT_reads"]
  elif config["Outputs"]["filtlong_dir"] != None:
    ont_reads = config["Outputs"]["filtlong_dir"] + "reads.ont.fastq.gz"

  for i in config["Inputs"]["Assemblies for postpolishing"]:
    input_assemblies[working_dir + config["Inputs"]["Assemblies for postpolishing"][i] + "_run_purgedups/"] = i
    base = os.path.splitext(os.path.basename(i))[0]
    targets_purgedups.append( working_dir + config["Inputs"]["Assemblies for postpolishing"][i] + "_run_purgedups/" + base + ".purged.fa")
   # print (targets_purgedups)
    if not os.path.exists(working_dir + config["Inputs"]["Assemblies for postpolishing"][i] + "_run_purgedups/logs"):
      os.makedirs(working_dir + config["Inputs"]["Assemblies for postpolishing"][i] + "_run_purgedups/logs")
    mappings_dir = working_dir + config["Inputs"]["Assemblies for postpolishing"][i] + "_run_purgedups/mappings"
    if not os.path.exists(mappings_dir):
      os.makedirs(mappings_dir)

#1- Define rule all
rule all_postpolishing:
  input:
    targets_purgedups
  log:
    logs_dir + str(date) + ".rule_all.out",
    logs_dir + str(date) + ".rule_all.err"


#1- Run purgedups 
if config["Parameters"]["run_purgedups"] == True:
  calcuts = ""
  if config["Purge_dups"]["calcuts_options"] != None:
    calcuts = config["Purge_dups"]["calcuts_options"]

  use rule purge_dups from purging_workflow with:
    input:
      assembly_in = lambda wildcards: input_assemblies[wildcards.dir + "/"],
      reads = ont_reads,
      mapping = "{dir}/mappings/{base_in}_minimap2.allreads.paf.gz"
    output:
      assembly_out = "{dir}/{base_in}.purged.fa" 
    params:
      module = config["Purge_dups"]["purgedups_module"],
      base = "{base_in}",
      dir = "{dir}",
      calcuts_opts = calcuts
    log:
      "{dir}/logs/" + str(date) + ".j%j.{base_in}.purge_dups.out",
      "{dir}/logs/" + str(date) + ".j%j.{base_in}.purge_dups.err",
    threads: config["Purge_dups"]["purgedups_cores"]

