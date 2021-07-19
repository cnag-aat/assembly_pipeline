from datetime import datetime
import os
import re
import subprocess

date = datetime.now().strftime('%Y%m%d.%H%M%S')

module curate_workflow:
  snakefile: "../modules/curation.rules.smk"
module eval_workflow:
  snakefile: "../modules/evaluate_assemblies.rules.smk"
module preprocess_workflow:
  snakefile: "../modules/process_reads.rules.smk"

working_dir = config["Outputs"]["base_dir"]
eval_dir = config["Outputs"]["eval_dir"]
scripts_dir = config["Inputs"]["scripts_dir"]
logs_dir = config["Parameters"]["logs_dir"]
if not os.path.exists(logs_dir):
  os.makedirs(logs_dir)

shell.prefix("echo 'Cluster jobid $SLURM_JOBID'; export PATH=" + scripts_dir + ":$PATH;")
keepfiles = config["Parameters"]["keep_intermediate"]
base = config["Parameters"]["base_name"]

##0. Define path for files and variables
targets_purgedups = []
input_assemblies = {}
minimap2 = {}
if config["Parameters"]["run_purgedups"] == True:
  for i in config["Inputs"]["Assemblies for curation"]:
    input_assemblies[working_dir + config["Inputs"]["Assemblies for curation"][i] + "_run_purgedups/"] = i
    base = os.path.splitext(os.path.basename(i))[0]
    targets_purgedups.append( working_dir + config["Inputs"]["Assemblies for curation"][i] + "_run_purgedups/" + base + ".purged.fa")
    if not os.path.exists(working_dir + config["Inputs"]["Assemblies for curation"][i] + "_run_purgedups/logs"):
      os.makedirs(working_dir + config["Inputs"]["Assemblies for curation"][i] + "_run_purgedups/logs")
    minimap2[base] = i
    mappings_dir = working_dir + config["Inputs"]["Assemblies for curation"][i] + "_run_purgedups/mappings"
    if not os.path.exists(mappings_dir):
      os.makedirs(mappings_dir)

#1- Define rule all
rule all_curation:
  input:
    targets_purgedups
  log:
    logs_dir + str(date) + ".rule_all.out",
    logs_dir + str(date) + ".rule_all.err"

##2- Obtain input reads
fastqs = {}
if config["Parameters"]["run_purgedups"] == True:
  if config["Inputs"]["ONT_reads"] == None:
    ont_dir = config["Inputs"]["ONT_dir"]
    ont_list = config["Wildcards"]["ONT_wildcards"].split(',')
    ont_reads = config["Outputs"]["filtlong_dir"] + "reads.ont.fastq.gz"
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
    final_fastq = "{dir}reads.{ext}"
  log:
    "{dir}logs/" + str(date) + ".concat.{ext}.out",
    "{dir}logs/" + str(date) + ".concat.{ext}.err"
  threads: config["Parameters"]["concat_cores"] 

  
#3- Perform mappings
if config["Parameters"]["run_purgedups"] == True:
  use rule align_ont from eval_workflow with:
    input:
      genome = lambda wildcards: minimap2[wildcards.name],
      reads = ont_reads,
    output:
      mapping = "{directory}/mappings/{name}_{ext}"
    params:
      env = config["Racon"]["Minimap environment"],
      align_opts = lambda wildcards:"ax map-ont" if wildcards.ext == "minimap2.bam" else "x map-ont",
      tmp = "$TMPDIR/{name}_{ext}.tmp"
    wildcard_constraints:
      ext = "minimap2.(.+)"
    log:
      "{directory}/logs/" + str(date) + ".{ext}.{name}.out",
      "{directory}/logs/" + str(date) + ".{ext}.{name}.err",
    threads: config["Parameters"]["minimap2_cores"]

#4- Run curation rules
if config["Parameters"]["run_purgedups"] == True:
  use rule purge_dups from curate_workflow with:
    input:
      assembly_in = lambda wildcards: input_assemblies[wildcards.dir + "/"],
      reads = ont_reads,
      mapping = "{dir}/mappings/{base_in}_minimap2.paf.gz"
    output:
      assembly_out = "{dir}/{base_in}.purged.fa" 
    params:
      module = config["Parameters"]["purgedups_module"],
      base = "{base_in}",
      dir = "{dir}"
    log:
      "{dir}/logs/" + str(date) + ".{base_in}.purge_dups.out",
      "{dir}/logs/" + str(date) + ".{base_in}.purge_dups.err",
    threads: config["Parameters"]["purgedups_cores"]

