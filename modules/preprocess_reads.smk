from datetime import datetime
import os
import re
import subprocess

date = datetime.now().strftime('%Y%m%d.%H%M%S')

module preprocess_workflow:
  snakefile: "../modules/process_reads.rules.smk"

working_dir = config["Outputs"]["base_dir"]
scripts_dir = config["Inputs"]["scripts_dir"]
#shell.prefix("export PATH=" + scripts_dir + ":$PATH;")

keepfiles = config["Parameters"]["keep_intermediate"]
base = config["Parameters"]["base_name"]

#0. Obtain inputs

fastqs = {}
ont_reads = ""
ont_fastqs = {}
kraken_ins = {}

if config["Inputs"]["ONT_reads"]:
  ont_reads = config["Inputs"]["ONT_reads"]      
ONT_filtered = config["Inputs"]["ONT_filtered"] 
ont_fastqs["filtered_ont"] = ONT_filtered

if config["Inputs"]["ONT_dir"] and config["Wildcards"]["ONT_wildcards"]:
  ont_dir = config["Inputs"]["ONT_dir"]
  ont_list = config["Wildcards"]["ONT_wildcards"].split(',')
  if ont_reads == "":
    ont_reads = config["Outputs"]["filtlong_dir"] + "reads.ont.fastq.gz"  
  extensions = ["fastq.gz"]
  for i in extensions:
    fastqs["ont."+i] = []
    for file in ont_list:
      fastqs["ont."+i].append(ont_dir + file + "." + i)

if ont_reads != "":
  ont_fastqs["raw_ont"] = ont_reads

if config["Parameters"]["run_kraken2"] == True:
  dbname = os.path.basename(config["Kraken2"]["database"])
  kraken_ins[os.path.dirname(ONT_filtered) + "/Kraken/filtered_ont/"+dbname+"/filtlong_"+dbname] = ONT_filtered
  logs = os.path.dirname(ONT_filtered)+"/Kraken/filtered_ont/"+dbname + "/filtlong_"+dbname+"_logs"
  if not os.path.exists(logs):
      os.makedirs(logs)
  if ont_reads != "":
    kraken_ins[os.path.dirname(ONT_filtered) + "/Kraken/raw_ont/"+dbname+"/raw_ont_"+dbname] = ont_reads
    logs = os.path.dirname(ont_reads)+"/Kraken/raw_ont/"+dbname + "/raw_ont_"+dbname+"_logs"
    if not os.path.exists(logs):
      os.makedirs(logs)

extra_filtlong_opts = config["Filtlong"]["options"]
if extra_filtlong_opts == None:
  extra_filtlong_opts = ""

#1. Preprocess reads

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

#2. EVALUATE INPUT READS

use rule nanoplot from preprocess_workflow with:
  input:
    fastq = lambda wildcards: ont_fastqs[wildcards.prefix]
  output:
    stats = report(config["Outputs"]["filtlong_dir"] + "nanostats/{prefix}/NanoStats.txt",
            caption="../report/nanostats.rst",
            category = "Process reads",
            subcategory = "{prefix}"),
    yield_len = report(config["Outputs"]["filtlong_dir"] + "nanostats/{prefix}/Yield_By_Length.png",
                caption="../report/nanostats.rst",
                category = "Process reads",
                subcategory = "{prefix}"),
    read_len = report(config["Outputs"]["filtlong_dir"] + "nanostats/{prefix}/WeightedHistogramReadlength.png",
               caption= "../report/nanostats.rst",
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
    read = lambda wildcards: kraken_ins[wildcards.base],
    database = config["Kraken2"]["database"],
    kmers = config["Kraken2"]["kmer_dist"]
  output:
    report = report("{base}.kraken2.report.txt",
             caption="../report/kraken.rst",
             category = "Process reads",
             subcategory = "Kraken reports"),
    readsout = "{base}.kraken2.seqs.out",
    abundance = "{base}.bracken_abundance.txt"
  params:
    additional = config["Kraken2"]["additional_opts"],
    prefix = lambda wildcards: os.path.basename(wildcards.base)
  log:
    "{base}_logs/" + str(date) + ".j%j.kraken.out",
    "{base}_logs/" + str(date) + ".j%j.kraken.err",
  benchmark:
     "{base}_logs/" + str(date) + ".j%j.kraken.benchmark.txt",
  conda:
    '../envs/kraken2.1.2.yaml'
  threads: config["Kraken2"]["threads"]