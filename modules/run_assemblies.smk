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
shell.prefix("export PATH=" + scripts_dir + ":$PATH;")

keepfiles = config["Parameters"]["keep_intermediate"]
base = config["Parameters"]["base_name"]

##0. Define path for files and variables

flye_dir = config["Outputs"]["flye_dir"]
nextdenovo_dir = config["Outputs"]["nextdenovo_dir"]
flye_assembly = config["Outputs"]["flye_out"]
nextdenovo_assembly = config["Outputs"]["nextdenovo_out"]

##2- Obtain input reads
fastqs = {}
ont_reads = ""
ONT_filtered = ""
ont_fastqs = {}
kraken_ins = {}
if config["Parameters"]["run_flye"] == True or config["Parameters"]["run_nextdenovo"] == True:
  ONT_filtered = config["Inputs"]["ONT_filtered"]
  ont_fastqs["filtered_ont"] = ONT_filtered
  if config["Parameters"]["run_kraken2"] == True:
    kraken_ins[os.path.dirname(ONT_filtered)] = ONT_filtered
  if not os.path.exists(ONT_filtered):
    if not os.path.exists(config["Outputs"]["filtlong_dir"] + "logs"):
      os.makedirs(config["Outputs"]["filtlong_dir"]  + "logs")
    if config["Inputs"]["ONT_reads"] == None:
      ont_dir = config["Inputs"]["ONT_dir"]
      ont_list = config["Wildcards"]["ONT_wildcards"].split(',')
      ont_reads =config["Outputs"]["filtlong_dir"] + "reads.ont.fastq.gz"
      extensions = ["fastq.gz"]
      for i in extensions:
        fastqs["ont."+i] = []
        for file in ont_list:
          fastqs["ont."+i].append(ont_dir + file + "." + i)
    else:
      ont_reads = config["Inputs"]["ONT_reads"] 
    ont_fastqs["raw_ont"] = ont_reads

if len(fastqs) > 0:
 use rule concat_reads from preprocess_workflow with:
  input:
    fastqs = lambda wildcards: fastqs[wildcards.ext]
  output:
    final_fastq = "{dir}reads.{ext}"
  log:
    "{dir}logs/" + str(date) + ".j%j.concat.{ext}.out",
    "{dir}logs/" + str(date) + ".j%j.concat.{ext}.err"
  benchmark:
    "{dir}logs/" + str(date) + ".concat.benchmark.{ext}.txt"
  conda:
    '../envs/ass_base.yaml'
  threads: config["Parameters"]["concat_cores"]  

if config["Inputs"]["ONT_filtered"] !=None and not os.path.exists(config["Inputs"]["ONT_filtered"]):
  extra_filtlong_opts = config["Filtlong"]["options"]
  if extra_filtlong_opts == None:
    extra_filtlong_opts = ""

  use rule filtlong from preprocess_workflow with:
    input:
      reads = ont_reads
    output:
      outreads = ONT_filtered
    params:
      minlen = config["Filtlong"]["Filtlong minlen"],
      min_mean_q = config["Filtlong"]["Filtlong min_mean_q"],
      opts = extra_filtlong_opts
    log:
      config["Outputs"]["filtlong_dir"] + "logs/" + str(date) + ".j%j.filtlong.out",
      config["Outputs"]["filtlong_dir"] + "logs/" + str(date) + ".j%j.filtlong.err"
    benchmark:
      config["Outputs"]["filtlong_dir"] + "logs/" + str(date) + ".filtlong.benchmark.txt"
    conda:
      '../envs/filtlong0.2.1.yaml'
    threads: config["Parameters"]["concat_cores"]  

#EVALUATE INPUT READS

use rule nanoplot from preprocess_workflow with:
  input:
    fastq = lambda wildcards: ont_fastqs[wildcards.prefix]
  output:
    stats = config["Outputs"]["filtlong_dir"] + "nanostats/{prefix}/NanoStats.txt",
    yield_len = report(config["Outputs"]["filtlong_dir"] + "nanostats/{prefix}/Yield_By_Length.png",
                caption="../report/nanostats.rst",
                category = "Process reads",
                subcategory = "{prefix}"),
    read_len = report(config["Outputs"]["filtlong_dir"] + "nanostats/{prefix}/WeightedHistogramReadlength.png",
               caption="../report/nanostats.rst",
               category = "Process reads",
               subcategory = "{prefix}"),
  params:
    outdir = config["Outputs"]["filtlong_dir"] + "nanostats/{prefix}/"
  log:
    config["Outputs"]["filtlong_dir"] + "logs/" + str(date) + ".j%j.NanoStats.{prefix}.out",
    config["Outputs"]["filtlong_dir"] + "logs/" + str(date) + ".j%j.NanoStats.{prefix}.err"
  benchmark:
    config["Outputs"]["filtlong_dir"] + "logs/" + str(date) + ".NanoStats.{prefix}.benchmark.txt"
  conda:
    '../envs/nanoplot1.40.0.yaml'
  threads: config["Parameters"]["concat_cores"]

use rule Kraken2 from preprocess_workflow with: 
  input:
    read = lambda wildcards: kraken_ins[wildcards.dir],
   # database = config["Kraken2"]["database"]
  output:
    report = "{dir}/kraken2.report",
  params:
    additional = config["Kraken2"]["additional_opts"]
  log:
    "{dir}/logs" + str(date) + ".j%j.kraken.out",
    "{dir}/logs/" + str(date) + ".j%j.kraken.err",
  benchmark:
    "{dir}/logs/" + str(date) + ".j%j.kraken.benchmark.txt",
  conda:
    '../envs/kraken2.1.2.yaml'
  threads: config["Kraken2"]["threads"]


##Run assemblers
if config["Parameters"]["run_flye"] == True:
  use rule flye from assembly_workflow with:
    input:
      reads = ONT_filtered,
    output:
      assembly = flye_assembly
    params:
      outdir = flye_dir,
      readtype = config["Parameters"]["lr_type"],
      pol_iterations = config["Flye"]["Flye polishing iterations"],
      other_flye_opts = config["Flye"]["options"]
    log:
      flye_dir + "logs/" + str(date) + ".j%j.flye.out",
      flye_dir + "logs/" + str(date) + ".j%j.flye.err"
    benchmark:
      flye_dir + "logs/" + str(date) + ".flye.benchmark.txt",
    conda:
      '../envs/flye2.9.1.yaml'
    threads: config["Flye"]["Flye cores"]

if config["Parameters"]["run_nextdenovo"] == True:
 use rule nextdenovo from assembly_workflow with:
  input:
    reads = ONT_filtered
  output:
    assembly = nextdenovo_assembly
  params:
    outdir = nextdenovo_dir,
    config = config["Parameters"]["ndconfFile"]
  log:
    nextdenovo_dir + "logs/" + str(date) + ".j%j.nextdenovo.out",
    nextdenovo_dir + "logs/" + str(date) + ".j%j.nextdenovo.err"
  benchmark:
    nextdenovo_dir + "logs/" + str(date) + ".j%j.nextdenovo.benchmark.txt"
  envmodules:
    "NextDenovo/2.5.0"
  threads: config["Nextdenovo"]["Nextdenovo cores"]