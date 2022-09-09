from datetime import datetime
import os
import re
import subprocess

module purging_workflow:
  snakefile: "../modules/purging.rules.smk"

##0. Define path for files and variables
input_assemblies = {}
postpolish = []
if config["Parameters"]["run_purgedups"] == True:
  for i in config["Inputs"]["Assemblies for postpolishing"]:
    input_assemblies[working_dir + config["Inputs"]["Assemblies for postpolishing"][i] + "_run_purgedups/"] = i
    base_postpolish = os.path.splitext(os.path.basename(i))[0]
    if not os.path.exists(working_dir + config["Inputs"]["Assemblies for postpolishing"][i] + "_run_purgedups/logs"):
      os.makedirs(working_dir + config["Inputs"]["Assemblies for postpolishing"][i] + "_run_purgedups/logs")
    mappings_dir = working_dir + config["Inputs"]["Assemblies for postpolishing"][i] + "_run_purgedups/mappings"
    if not os.path.exists(mappings_dir):
      os.makedirs(mappings_dir)    
    postpolish.append(i)
    postpolish.append( working_dir + config["Inputs"]["Assemblies for postpolishing"][i] + "_run_purgedups/" + base_postpolish + ".purged.fa")
    minimap2[base_postpolish] = i

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
      assembly_out = "{dir}/{base_in}.purged.fa",
      plot = report ("{dir}/{base_in}.PB.cov.png",
             caption="../report/purgedups.rst",
             category = "Purgedups")
    params:
      scripts_dir = scripts_dir,
      base = "{base_in}",
      dir = "{dir}",
      calcuts_opts = calcuts
    log:
      "{dir}/logs/" + str(date) + ".j%j.{base_in}.purge_dups.out",
      "{dir}/logs/" + str(date) + ".j%j.{base_in}.purge_dups.err"
    benchmark:
      "{dir}/logs/" + str(date) + ".{base_in}.purge_dups.benchmark.txt"
    conda:
      "../envs/purge_dups1.2.6.yaml"
    threads: config["Purge_dups"]["purgedups_cores"]

