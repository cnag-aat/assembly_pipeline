from datetime import datetime
import os
import re
import subprocess

date = datetime.now().strftime('%Y%m%d.%H%M%S')

module eval_workflow:
  snakefile: "../modules/evaluate_assemblies.rules.smk"

working_dir = config["Outputs"]["base_dir"]
scripts_dir = config["Inputs"]["scripts_dir"]
shell.prefix("export PATH=" + scripts_dir + ":$PATH;")

keepfiles = config["Parameters"]["keep_intermediate"]
base = config["Parameters"]["base_name"]
eval_dir = config["Outputs"]["eval_dir"]

##0. Define path for files and variables

inputs = config["Inputs"]["Assemblies for polishing"]

assemblies = []
if config["Finalize"]["final Evaluations"] == True:
  for i in inputs:
    assemblies.append(i)
  if config["Parameters"]["run_flye"] == True:
    flye_assembly = config["Outputs"]["flye_out"]
    assemblies.append(flye_assembly)
  if config["Parameters"]["run_nextdenovo"] == True:
    nextdenovo_assembly = config["Outputs"]["nextdenovo_out"]
    assemblies.append(nextdenovo_assembly)

in_files = {}
evals_dir = {}
BuscoSummaries = []
StatsFiles = []
for file in assemblies:
  ass_base = os.path.splitext(os.path.basename(file))[0]
  basedirname = os.path.basename(os.path.dirname(file))
  #if basedirname == "hypo" or basedirname == "rmp" or basedirname == "nextpolish":
   #   basedirname = os.path.basename(os.path.dirname(os.path.dirname(busco)))
  evalassdir =  eval_dir + basedirname + "/"
  if not os.path.exists(evalassdir + "logs/"):
    os.makedirs(evalassdir + "logs/")
  if not os.path.exists(evalassdir + "stats"):
    os.makedirs(evalassdir + "stats")
  buscodir = evalassdir + "busco/"
  if not os.path.exists(buscodir):
    os.makedirs(buscodir)
  in_files[evalassdir + ass_base] = file
  evals_dir[evalassdir + ass_base] = evalassdir
  BuscoSummaries.append(buscodir + ass_base + ".short_summary.txt")
  StatsFiles.append(evalassdir + "stats/" + ass_base + ".stats.txt")

#1- Run evaluations

use rule get_stats from eval_workflow with:
  input:
    assembly =  lambda wildcards: in_files[eval_dir + wildcards.dir + "/" + wildcards.buscobase],
  output: 
    nseries = report(eval_dir + "{dir}/stats/{buscobase}.nseries.txt",
              caption="../report/stats.rst",
              category = "Evaluate assemblies",
              subcategory = "{dir}"),
    stats = report (eval_dir + "{dir}/stats/{buscobase}.stats.txt",
              caption="../report/stats.rst",
              category = "Evaluate assemblies",
              subcategory = "{dir}"),
    gaps = report(eval_dir + "{dir}/stats/{buscobase}.gaps.txt",
              caption="../report/stats.rst",
              category = "Evaluate assemblies",
              subcategory = "{dir}")
  params:
    outbase = "{buscobase}",
    scripts_dir = scripts_dir
  log:
    eval_dir + "{dir}/logs/" + str(date) + ".j%j.get_stats.{buscobase}.out",
    eval_dir + "{dir}/logs/" + str(date) + ".j%j.get_stats.{buscobase}.err"
  benchmark:  
    eval_dir + "{dir}/logs/" + str(date) + "get_stats.{buscobase}.benchmark.out",
  conda:
    '../envs/ass_base.yaml'
  threads: 1

if config["Finalize"]["BUSCO lineage"] != None:
  use rule run_busco from eval_workflow with:
    input:
      assembly = lambda wildcards: in_files[eval_dir + wildcards.dir + "/" + wildcards.buscobase],
      lineage = config["Finalize"]["BUSCO lineage"]
    output:
      summary = report(eval_dir + "{dir}/busco/{buscobase}.short_summary.txt",
                caption="../report/busco.rst",
                category = "Evaluate assemblies",
                subcategory = "{dir}"),
      full = eval_dir + "{dir}/busco/{buscobase}.full_table.tsv",
    params:
      out_path = lambda wildcards: evals_dir[eval_dir + wildcards.dir +"/" + wildcards.buscobase],
      odb = os.path.basename(config["Finalize"]["BUSCO lineage"]),
      buscobase = lambda wildcards:  wildcards.buscobase,
      rmcmd = "echo 'Removing BUSCO run dir:{params.out_path}{params.buscobase}'; \
            rm -r {params.out_path}{params.buscobase};" if keepfiles == False else "" 
    log:
      eval_dir + "{dir}/logs/" + str(date) + ".j%j.busco.{buscobase}.out",
      eval_dir + "{dir}/logs/" + str(date) + ".j%j.busco.{buscobase}.err",
    benchmark:
      eval_dir + "{dir}/logs/" + str(date) + ".j%j.busco.{buscobase}.benchmark.txt"
    conda:
      '../envs/busco5.4.0.yaml'
    threads: config["Parameters"]["busco_cores"]

use rule finalize from eval_workflow with:
  input:
    #assembly = lambda wildcards: expand(in_files[eval_dir + wildcards.dir + "/" + wildcards.buscobase]),
    assembly = assemblies,
    buscos = BuscoSummaries,
    stats= StatsFiles
    #buscos = rules.run_busco.output.summary,
    #stats = rules.get_stats.output.stats
  output:
    output = config["Outputs"]["stats_out"]
  log:
    logs_dir + str(date) + ".j%j.finalize.out",
    logs_dir + str(date) + ".j%j.finalize.err"
  benchmark:
    logs_dir + str(date) + ".finalize.benchmark.txt"
  conda:
    '../envs/ass_base.yaml'
  threads: 1
 