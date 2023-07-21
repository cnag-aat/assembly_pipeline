from datetime import datetime
import os
import re
import subprocess

module purging_workflow:
  snakefile: "../modules/purging.rules.smk"
module scaffolding_workflow:
  snakefile: "../modules/scaffolding_rules.smk"

##0. Define path for files and variables
input_assemblies = {}
postpolish = []
scripts_dir = config["Inputs"]["scripts_dir"]
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

tigmint_assemblies = {}
if config["Parameters"]["run_tigmint"] == True:
  for i in config["Inputs"]["Assemblies for postpolishing"]:
    step = config["Inputs"]["Assemblies for postpolishing"][i]
    base = os.path.splitext(os.path.basename(i))[0]
    assembly = i
    if config["Parameters"]["run_purgedups"] == True:
      pstep = step.split('_')[0]
      nstep = pstep.replace('s','')
      cstep = float(nstep) + 1
      step = "s0" + str(cstep) + "_" + step.split('_')[0].replace('s','p')
      assembly =  working_dir + config["Inputs"]["Assemblies for postpolishing"][i] + "_run_purgedups/" + base + ".purged.fa"

    tigmint_assemblies[working_dir + step + "_Scaffolding_tigmint_ARKS_LINKS/"] = assembly
    base = os.path.splitext(os.path.basename(assembly))[0]
    postpolish.append( working_dir + step + "_Scaffolding_tigmint_ARKS_LINKS/" + base + ".10X.scaffolds.fa")
    if not os.path.exists(working_dir +step + "_Scaffolding_tigmint_ARKS_LINKS/logs"):
      os.makedirs(working_dir + step + "_Scaffolding_tigmint_ARKS_LINKS/logs")

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

#2- Run tigmint
if config["Parameters"]["run_tigmint"] == True:
  tigmint_opts = ""
  if config["scaffolding_10X"]["tigmint_options"] != None:
    tigmint_opts = config["scaffolding_10X"]["tigmint_options"]

  if config["Inputs"]["processed_10X"] != None or len(config["Inputs"]["raw_10X"])>0:
    r10X_dir = config["Inputs"]["processed_10X"]
    r10X_reads = r10X_dir + "reads.illumina10X.barcoded.fastq.gz"

  use rule scaffolding_10X from scaffolding_workflow with:
    input:
      assembly = lambda wildcards: tigmint_assemblies[wildcards.dir + "/"],
      reads = r10X_reads
    output:
      scaffolded = "{dir}/{base_in}.10X.scaffolds.fa"
    params:
      base = "{base_in}",
      dir = "{dir}",
      opts = tigmint_opts,
      scripts = scripts_dir,
      rmcmd =  lambda wildcards: "echo 'Removing bam file'; "  \
            # "rm $base_ass.$base_reads.sortbx.bam;" \
             "echo 'Removing tmp dir'; " \
             "rm -r $TMPDIR;"
             if keepfiles == False else "" 
    log:
      "{dir}/logs/" + str(date) + ".j%j.{base_in}.10X_scaffolding.out",
      "{dir}/logs/" + str(date) + ".j%j.{base_in}.10X_scaffolding.err",
    benchmark:
      "{dir}/logs/" + str(date) + ".{base_in}.10X_scaffolding.benchmark.txt",
    threads: config["scaffolding_10X"]["tigmint_cores"]

