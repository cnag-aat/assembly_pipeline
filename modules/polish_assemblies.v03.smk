from datetime import datetime
import os
import re
import subprocess

date = datetime.now().strftime('%Y%m%d.%H%M%S')

module eval_workflow:
  snakefile: "../modules/evaluate_assemblies.rules.smk"
module preprocess_workflow:
  snakefile: "../modules/process_reads.rules.smk"

def get_targets(assembly):
  basename = os.path.splitext(os.path.basename(assembly))[0]
  rundir = inputs[assembly]
 # print (basename)
  mappings_dir = rundir + "mappings/"
  rmp_dir = rundir + "rmp/"
  nextpolish_dir = rundir + "nextpolish/"
  hypo_dir = rundir + "hypo/"
  make_logs = []
  make_logs.append(rundir)
  if rr > 0 or nor > 0 or mr >0 or hr >0:
    minimap2[basename] = assembly  
    if not os.path.exists(mappings_dir):
      os.makedirs(mappings_dir)
  if rr > 0:
    make_logs.append(rmp_dir)
    racon_in[rmp_dir + basename + ".racon1.fasta"] = assembly
    racon_map[rmp_dir + basename + ".racon1.fasta"] = mappings_dir + basename + "_minimap2.paf.gz"
    targets.append(rmp_dir + basename + ".racon" + str(rr) + ".fasta")
    if mr == 0 and pr == 0:
      terminalassemblies.append(rmp_dir + basename + ".racon" + str(rr) + ".fasta")
    else:
      intermediateassemblies.append(rmp_dir + basename + ".racon" + str(rr) + ".fasta")

  if mr > 0:
    make_logs.append(rmp_dir)
    base_medaka = basename
    if rr > 0:
      base_medaka += ".racon" + str(rr)
      medaka_in[rmp_dir + base_medaka + ".medaka1.fasta"] = rmp_dir + base_medaka + ".fasta"
      minimap2[base_medaka] = rmp_dir + base_medaka + ".fasta"
    else:
      medaka_in[rmp_dir + base_medaka + ".medaka1.fasta"] = assembly
    medaka_map[rmp_dir + base_medaka + ".medaka1.fasta"] = mappings_dir + base_medaka + "_minimap2.bam"
    targets.append(rmp_dir + base_medaka + ".medaka" + str(mr) + ".fasta")   
    if pr == 0:
      terminalassemblies.append(rmp_dir + base_medaka + ".medaka" + str(mr) + ".fasta") 
    else:
      intermediateassemblies.append(rmp_dir + base_medaka + ".medaka" + str(mr) + ".fasta")

  if nor > 0:
    make_logs.append(nextpolish_dir)
    nor_in[nextpolish_dir + basename + ".nextpolish_ont1.fasta"] = assembly
    nor_map[nextpolish_dir + basename + ".nextpolish_ont1.fasta"] = mappings_dir + basename + "_minimap2.bam"
    targets.append(nextpolish_dir + basename + ".nextpolish_ont" + str(nor) + ".fasta")
    if nir == 0:
      terminalassemblies.append(nextpolish_dir + basename + ".nextpolish_ont" + str(nor) + ".fasta")
    else:
      intermediateassemblies.append(nextpolish_dir + basename + ".nextpolish_ont" + str(nor) + ".fasta")

  for i in (range(1,rr)):
    pref = basename + ".racon"
    minimap2[pref + str(i)] = rmp_dir + pref + str(i) + ".fasta"
    targets.append(rmp_dir + pref + str(i) + ".fasta")
    n = i+1
    racon_in[rmp_dir + pref + str(n) + ".fasta"]= rmp_dir + pref + str(i) + ".fasta"
    racon_map[rmp_dir + pref + str(n) + ".fasta"] = mappings_dir + pref + str(i) + "_minimap2.paf.gz"
    intermediateassemblies.append(rmp_dir + basename + ".racon" + str(i) + ".fasta")
  for i in (range(1,mr)):
    pref = base_medaka + ".medaka"
    minimap2[pref + str(i)] = rmp_dir + pref + str(i) + ".fasta"
    targets.append(rmp_dir + pref + str(i) + ".fasta")
    intermediateassemblies.append(rmp_dir + base_medaka + ".medaka" + str(i) + ".fasta")
    n=i+1
    medaka_in[rmp_dir + pref + str(n) + ".fasta"] = rmp_dir + pref + str(i) + ".fasta"
    medaka_map[rmp_dir + pref + str(n) + ".fasta"] = mappings_dir + pref + str(i) + "_minimap2.bam"
  for i in (range(1,nor)):
    pref = basename + ".nextpolish_ont"
    minimap2[pref + str(i)] = nextpolish_dir + pref + str(i) + ".fasta"
    targets.append(nextpolish_dir + pref + str(i) + ".fasta")
    intermediateassemblies.append(nextpolish_dir + basename + ".nextpolish_ont" + str(i) + ".fasta")
    n = i+1
    nor_in[nextpolish_dir + pref + str(n) + ".fasta"] = nextpolish_dir + pref + str(i) + ".fasta"
    nor_map[nextpolish_dir + pref + str(n) + ".fasta"] = mappings_dir + pref + str(i) + "_minimap2.bam"

  if (nor == 0 and nir > 0) or (rr == 0 and mr == 0 and pr > 0) or hr > 0:
    bwa[basename] = assembly 
    if not os.path.exists(mappings_dir):
      os.makedirs(mappings_dir) 
  if pr > 0:
    make_logs.append(rmp_dir)
    pilon_base = basename
    if rr > 0:
      pilon_base+=".racon" + str(rr)
    if mr >0:
      pilon_base += ".medaka" + str(mr)
    if rr > 0  or mr >0:
      bwa[pilon_base] = rmp_dir + pilon_base + ".fasta"
      pilon_in[rmp_dir + pilon_base + ".pilon1.fasta"] = rmp_dir + pilon_base + ".fasta"
    else:
      pilon_in[rmp_dir + pilon_base + ".pilon1.fasta"] = assembly 
    pilon_map[rmp_dir + pilon_base + ".pilon1.fasta"] = mappings_dir + pilon_base + "_bwa.bam" 
    targets.append(rmp_dir + pilon_base + ".pilon" + str(pr) + ".fasta")
    terminalassemblies.append(rmp_dir + pilon_base + ".pilon" + str(pr) + ".fasta") 
    for i in (range(1,pr)):    
      pref = pilon_base + ".pilon"
      bwa[pref + str(i)] = rmp_dir +pref+ str(i) + ".fasta"
      targets.append(rmp_dir + pref + str(i) + ".fasta"),
      n=i+1
      pilon_in[rmp_dir + pref + str(n) + ".fasta"] = rmp_dir + pref + str(i) + ".fasta"
      pilon_map[rmp_dir + pref + str(n) + ".fasta"] = mappings_dir + pref + str(i) + "_bwa.bam" 
      intermediateassemblies.append(rmp_dir + pilon_base + ".pilon" + str(i) + ".fasta") 
    
  if nir >0:
    make_logs.append(nextpolish_dir)
    nir_base = basename
    if nor >0:
      nir_base += ".nextpolish_ont" + str(nor)
      bwa[nir_base] = nextpolish_dir + nir_base + ".fasta"
    #  bwa[nir_base + ".nextpolish_ill1_tmp"] = nextpolish_dir + nir_base + ".nextpolish_ill1_tmp.fa"
      nir_in[nextpolish_dir + nir_base + ".nextpolish_ill1_tmp.fa"] = nextpolish_dir + nir_base + ".fasta"
    else:
      nir_in[nextpolish_dir + nir_base + ".nextpolish_ill1_tmp.fa"] = assembly
    bwa[nir_base + ".nextpolish_ill1_tmp"] = nextpolish_dir + nir_base + ".nextpolish_ill1_tmp.fa"
    nir_in[nextpolish_dir + nir_base + ".nextpolish_ill1.fasta"] = nextpolish_dir + nir_base + ".nextpolish_ill1_tmp.fa"
    nir_map[nextpolish_dir +  nir_base + ".nextpolish_ill1_tmp.fa"] = mappings_dir + nir_base + "_bwa.bam"
    nir_map[nextpolish_dir + nir_base + ".nextpolish_ill1.fasta"] = mappings_dir +  nir_base + ".nextpolish_ill1_tmp_bwa.bam" 
    targets.append(nextpolish_dir + nir_base + ".nextpolish_ill" + str(nir) + ".fasta")
    terminalassemblies.append(nextpolish_dir + nir_base + ".nextpolish_ill" + str(nir) + ".fasta")
    tasks[nextpolish_dir + nir_base + ".nextpolish_ill1_tmp.fa"] = 1
    tasks[nextpolish_dir +  nir_base + ".nextpolish_ill1.fasta"] = 2
    for i in (range(1,nir)):
      pref = nir_base + ".nextpolish_ill"
      n = i+1
      bwa[pref + str(n) + "_tmp"] = nextpolish_dir + pref + str(n) + "_tmp.fa"       
      bwa[pref + str(i)] = nextpolish_dir + pref + str(i) + ".fasta"
      targets.append(nextpolish_dir + pref + str(i) + ".fasta")
      intermediateassemblies.append(nextpolish_dir + nir_base + ".nextpolish_ill" + str(i) + ".fasta")
      nir_in[nextpolish_dir + pref + str(n) + "_tmp.fa"] = nextpolish_dir + pref + str(i) + ".fasta"
      nir_in[nextpolish_dir + pref + str(n) + ".fasta"] = nextpolish_dir + pref + str(n) + "_tmp.fa"
      nir_map[nextpolish_dir + pref + str(n) + "_tmp.fa"] = mappings_dir + pref + str(i) + "_bwa.bam"
      nir_map[nextpolish_dir + pref + str(n) + ".fasta"] = mappings_dir + pref + str(n) + "_tmp_bwa.bam"
      tasks[nextpolish_dir + pref + str(n) + "_tmp.fa"] = 1
      tasks[nextpolish_dir + pref + str(n) + ".fasta"] = 2
  if hr >0:
    make_logs.append(hypo_dir)
    hypo_in[hypo_dir + basename + ".hypo1.fasta"] = assembly
    hypo_lrmap[hypo_dir + basename + ".hypo1.fasta"] = mappings_dir + basename + "_minimap2.bam"
    hypo_srmap[hypo_dir + basename + ".hypo1.fasta"] = mappings_dir + basename + "_bwa.bam"
    targets.append(hypo_dir + basename + ".hypo" + str(hr) + ".fasta")
    terminalassemblies.append(hypo_dir + basename + ".hypo" + str(hr) + ".fasta")
    for i in (range(1,hr)):
      bwa[basename + ".hypo" + str(i)] = hypo_dir + basename + ".hypo" + str(i) + ".fasta"
      minimap2[basename + ".hypo" + str(i)] = hypo_dir + basename + ".hypo" + str(i) + ".fasta"
      targets.append(hypo_dir + basename + ".hypo" + str(i) + ".fasta")
      intermediateassemblies.append(hypo_dir + basename + ".hypo" + str(i) + ".fasta")
      n=i+1
      hypo_in[hypo_dir + basename + ".hypo" + str(n) + ".fasta"] = hypo_dir + basename + ".hypo" + str(i) + ".fasta"
      hypo_lrmap[hypo_dir + basename + ".hypo" + str(n) + ".fasta"] = mappings_dir + basename + ".hypo" + str(i) + "_minimap2.bam"
      hypo_srmap[hypo_dir + basename + ".hypo" + str(n) + ".fasta"] = mappings_dir + basename + ".hypo" + str(i) + "_bwa.bam"
    if config["Hypo"]["options"] != None:
      additional_hypo_opts = config["Hypo"]["options"]
  for d in make_logs:
    if not os.path.exists(d + "logs"):
      os.makedirs(d+"logs")
##Main  
working_dir = config["Outputs"]["base_dir"]
eval_dir = config["Outputs"]["eval_dir"]
scripts_dir = config["Inputs"]["scripts_dir"]
logs_dir = config["Parameters"]["logs_dir"]
if not os.path.exists(logs_dir):
  os.makedirs(logs_dir)

rr = config["Parameters"]["racon_rounds"]
pr = config["Parameters"]["pilon_rounds"]
mr = config["Parameters"]["medaka_rounds"]
nor = config["Parameters"]["nextpolish_ont_rounds"]
nir = config["Parameters"]["nextpolish_ill_rounds"]
hr = config["Parameters"]["hypo_rounds"]
keepfiles = config["Parameters"]["keep_intermediate"]
base = config["Parameters"]["base_name"]

##0. Define path for files and variables
targets = []
terminalassemblies = []
intermediateassemblies = []
minimap2={}
bwa = {}
nir_base = "nextpolish_ill"
pilon_base = "pilon"
racon_opts = ""
if rr > 0:
  racon_in = {}
  racon_map = {}
  if config["Racon"]["options"] != None:
    racon_opts = config["Racon"]["options"] 
  racon_shell_command = "module purge; module load GCC/9.3.0 python3; export PATH=" + config["Racon"]["Racon dir"] + "/scripts:$PATH; export PATH=" + config["Racon"]["Racon dir"] + "/build/bin:$PATH;"

medaka_consensus_opts = ""
if mr > 0:
  medaka_in = {}
  medaka_map = {}
  if config["Parameters"]["medaka_cores"] > 8:
    medaka_threads = "8" 
  else:
    medaka_threads = config["Parameters"]["medaka_cores"]

  if config["Medaka"]["consensus options"] != None:
    consensus_opts = config["Medaka"]["consensus options"] 

nextpolish_lrtype= ""
if nor > 0:
  nor_in = {}
  nor_map = {}
  if re.search("nano",config["Parameters"]["lr_type"],):
    nextpolish_lrtype = "ont"
  elif re.match("pacbio-corr", config["Parameters"]["lr_type"]):
    nextpolish_lrtype = "hifi"
  elif re.match("pacbio-raw", config["Parameters"]["lr_type"], ):
    nextpolish_lrtype = "clr"

if pr > 0:
  pilon_in = {}
  pilon_map = {}
  if config["Pilon"]["JAVA options"] == None:
    java_opt = ""
  else:
    java_opt = "-" + config["Pilon"]["JAVA options"] 

if nir >0:
  nir_in = {}
  nir_map = {}
  tasks = {}

additional_hypo_opts = ""
if hr >0:
  hypo_in = {}
  hypo_lrmap = {}
  hypo_srmap = {}

inputs = config["Inputs"]["Assemblies for polishing"]
for assembly_in in inputs:
  get_targets(assembly_in)


assemblies4Busco = []
if config["Finalize"]["intermediate BUSCOs"] == True or config["Finalize"]["final BUSCOs"] == True:
  targets.append(config["Outputs"]["stats_out"])
if config["Finalize"]["intermediate BUSCOs"] == True:
  for i in intermediateassemblies:
    assemblies4Busco.append(i)
if config["Finalize"]["final BUSCOs"] == True:
  for i in terminalassemblies:
    assemblies4Busco.append(i)
  for i in inputs:
    assemblies4Busco.append(i)

busco_in = {}
evals_dir = {}
BuscoSummaries = []
if len(assemblies4Busco) > 0:
  for busco in assemblies4Busco:
    buscobase = os.path.splitext(os.path.basename(busco))[0]
    basedirname = os.path.basename(os.path.dirname(os.path.dirname(busco)))
    evaldir =  eval_dir + basedirname + "/"
    if not os.path.exists(evaldir + "logs/"):
      os.makedirs(evaldir + "logs/")
    buscodir = evaldir + "busco/"
    if not os.path.exists(buscodir):
      os.makedirs(buscodir)
    busco_in[evaldir + buscobase] = busco
    evals_dir[evaldir + buscobase] = evaldir
    BuscoSummaries.append(buscodir + buscobase + ".short_summary.txt")

assemblies4Merqury = []
if config["Finalize"]["Merqury db"] != None:
  if not os.path.exists(config["Finalize"]["Merqury db"]):
    targets.append(config["Finalize"]["Merqury db"])
  if config["Finalize"]["final BUSCOs"] == True:
    for i in terminalassemblies:
      assemblies4Merqury.append(i)
    for i in inputs:
      assemblies4Merqury.append(i)

merq_in = {}
MerqurySummaries = []
MerquryQV = []
if len(assemblies4Merqury) > 0:
  for merq in assemblies4Merqury:
    merqbase = os.path.splitext(os.path.basename(merq))[0]
    basedirname = os.path.basename(os.path.dirname(os.path.dirname(merq)))
    evaldir =  eval_dir + basedirname + "/"
    merqdir = evaldir + "merqury/" + merqbase + "/"
    if not os.path.exists(merqdir):
      os.makedirs(merqdir)
    merq_in[evaldir + merqbase] = merq
    MerqurySummaries.append(merqdir + "completeness.stats")
    MerquryQV.append(merqdir + merqbase + ".qv")

#1- Define rule all
rule all_polishing:
  input:
    targets
  log:
    logs_dir + str(date) + ".rule_all.out",
    logs_dir + str(date) + ".rule_all.err"

##2- Obtain input reads
fastqs = {}
if rr > 0 or mr > 0 or nor or hr >0:
  ONT_filtered = config["Inputs"]["ONT_filtered"]
  if not os.path.exists(ONT_filtered):
    if not os.path.exists(config["Outputs"]["filtlong_dir"] + "logs"):
      os.makedirs(config["Outputs"]["filtlong_dir"]  + "logs")
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

longranger_inputs = {}
pe1_reads = None
pe2_reads = None
r10X_reads = None
if pr > 0 or nir>0 or hr>0:
  if config["Inputs"]["ILLUMINA_pe1"] == None and config["Inputs"]["ILLUMINA_pe2"] == None and config["Inputs"]["ILLUMINA_10X"] == None:
    if config["Inputs"]["illumina_dir"] != None:
      illumina_dir = config["Inputs"]["illumina_dir"]
      illumina_list = config["Wildcards"]["illumina_wildcards"].split(',')
      extensions = ["1.fastq.gz", "2.fastq.gz"]
      pe1_reads = working_dir + "reads.illumina.1.fastq.gz"
      pe2_reads = working_dir + "reads.illumina.2.fastq.gz"
      for i in extensions:
        fastqs["illumina." + i] = []
        for file in illumina_list:
          fastqs["illumina." + i].append(illumina_dir + file + "." + i)
    if config["Inputs"]["processed_10X"] != None or config["Inputs"]["raw_10X"] != None:
      r10X_list = config["Wildcards"]["10X_wildcards"].split(',')
      r10X_dir = config["Inputs"]["processed_10X"]
      if config["Inputs"]["raw_10X"] != None:
        raw_10X_dir = config["Inputs"]["raw_10X"]
        for i in r10X_list:
          longranger_inputs[i] = [i]
      extensions = ["barcoded.fastq.gz"]
      r10X_reads = r10X_dir + "reads.illumina10X.barcoded.fastq.gz"
      if not os.path.exists(r10X_dir + "logs/"):
        os.makedirs(r10X_dir + "logs/")
      for i in extensions:
        fastqs["illumina10X." + i] = []
        for file in r10X_list:
          fastqs["illumina10X." + i].append(r10X_dir + file + ".lr." + i)

  elif config["Inputs"]["ILLUMINA_10X"] != None:
    r10X_reads = config["Inputs"]["ILLUMINA_10X"]
  else:
    pe1_reads = config["Inputs"]["ILLUMINA_pe1"].split(',')[0]
    pe2_reads = config["Inputs"]["ILLUMINA_pe2"].split(',')[0]

if config["Finalize"]["Merqury db"]:
  if not os.path.exists(config["Finalize"]["Merqury db"]):
    meryl_loc = os.path.dirname(config["Finalize"]["Merqury db"]) + "/tmp_meryl/"
    if not os.path.exists(meryl_loc):
      os.makedirs(meryl_loc)
    reads_loc = {}
    meryl_dbs = []
    if config["Inputs"]["processed_10X"]!= None or config["Inputs"]["raw_10X"] != None:
      r10X_list = config["Wildcards"]["10X_wildcards"].split(',')
      r10X_dir = config["Inputs"]["processed_10X"]
      if config["Inputs"]["raw_10X"] != None: 
        raw_10X_dir = config["Inputs"]["raw_10X"]
        for i in r10X_list:
          longranger_inputs[i] = [i]
      if not os.path.exists(r10X_dir + "logs/"):
        os.makedirs(r10X_dir + "logs/")
      extensions = [".barcoded.fastq.gz"]
      for i in r10X_list:
        for e in extensions:
          reads_loc[i + e] = r10X_dir + i + e
          meryl_dbs.append(i + e)
    elif r10X_reads != None:
      reads_loc[os.path.basename(r10X_reads)] = r10X_reads
      meryl_dbs.append(os.path.basename(r10X_reads))
    if pe1_reads != None:
      reads_loc[os.path.basename(pe1_reads)] = pe1_reads
      meryl_dbs.append(os.path.basename(pe1_reads))
    if pe2_reads != None:
      reads_loc[os.path.basename(pe2_reads)] = pe2_reads
      meryl_dbs.append(os.path.basename(pe2_reads))

if len(longranger_inputs) > 0:
  use rule long_ranger from preprocess_workflow with: 
    input: 
      mkfastq_dir = raw_10X_dir
    output:
      fastq_out = r10X_dir + "{bname}.lr.barcoded.fastq.gz",
      sum_out = r10X_dir + "{bname}.lr.barcoded.summary.csv"
    params:
      path = config["Parameters"]["longranger_path"],
      outdir = r10X_dir,
      sample = lambda wildcards: longranger_inputs[wildcards.bname]
    log:
      r10X_dir + "logs/" + str(date) + ".longranger.{bname}.out",
      r10X_dir + "logs/" + str(date) + ".longranger.{bname}.err"
    threads: config["Parameters"]["longranger_cores"] 

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

if not os.path.exists(ONT_filtered):
  extra_filtlong_opts = config["Filtlong"]["options"]
  if extra_filtlong_opts == None:
    extra_filtlong_opts = ""
  use rule filtlong from preprocess_workflow with:
    input:
      reads = ont_reads
    output:
      outreads = ONT_filtered
    params:
      path = config["Filtlong"]["Filtlong path"],
      minlen = config["Filtlong"]["Filtlong minlen"],
      min_mean_q = config["Filtlong"]["Filtlong min_mean_q"],
      opts = extra_filtlong_opts
    log:
      config["Outputs"]["filtlong_dir"] + "logs/" + str(date) + ".filtlong.out",
      config["Outputs"]["filtlong_dir"] + "logs/" + str(date) + ".filtlong.err"
    threads: config["Parameters"]["concat_cores"] 

if config["Finalize"]["Merqury db"]:
  if not os.path.exists(config["Finalize"]["Merqury db"]):
    use rule build_meryl_db from preprocess_workflow with:
      input:
        fastq = lambda wildcards: reads_loc[wildcards.db]
      output:
        out_dir = directory(meryl_loc + "{db}.meryl")  
      params:
        environment = config["Finalize"]["Merqury environment"],
        kmer = config["Finalize"]["Meryl K"],
      log:
        logs_dir + str(date) + ".build_meryl.{db}.out",
        logs_dir + str(date) + ".build_meryl.{db}.err" 

    use rule concat_meryl from preprocess_workflow with:
      input:
        input_run = lambda wildcards: expand(rules.build_meryl_db.output.out_dir, db=meryl_dbs)
      output:
        meryl_all = directory(config["Finalize"]["Merqury db"])
      params:
        environment = config["Finalize"]["Merqury environment"],
      log:
        logs_dir + str(date) + ".concat_meryl.out",
        logs_dir + str(date) + ".concat_meryl.err" 

#3- Perform mappings
if rr > 0 or mr > 0 or nor>0 or hr>0:
  use rule align_ont from eval_workflow with:
    input:
      genome = lambda wildcards: minimap2[wildcards.name],
      reads = ont_reads
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

if pr > 0 or nir>0 or hr>0:
 use rule align_illumina from eval_workflow with:
  input:
    genome = lambda wildcards: bwa[wildcards.name],
    reads = [pe1_reads, pe2_reads] if r10X_reads == None else r10X_reads
  output:
    mapping = "{directory}/mappings/{name}_bwa.bam",
  log:
    "{directory}/logs/" + str(date) + ".bwa.{name}.out",
    "{directory}/logs/" + str(date) + ".bwa.{name}.err",
  threads: config["Parameters"]["BWA_cores"]

#4- Run polishers

if rr > 0:
  rule racon:
    input:
      assembly = lambda wildcards: racon_in[wildcards.directory + "/rmp/" + wildcards.base + ".racon" + wildcards.param + ".fasta"],
      reads = ont_reads, 
      mapping = lambda wildcards: racon_map[wildcards.directory + "/rmp/" + wildcards.base + ".racon" + wildcards.param + ".fasta"], 
    output:
      polished = "{directory}/rmp/{base}.racon{param}.fasta"
    wildcard_constraints:
      param = "\d+"
    params:
      racon_env = config["Racon"]["Racon dir"],
      opts = racon_opts
    log:
      "{directory}/rmp/logs/" + str(date) + "{base}.racon{param}.out",
      "{directory}/rmp/logs/" + str(date) + "{base}.racon{param}.err",
    threads: config["Parameters"]["racon_cores"]
    run:
      shell(
        "{racon_shell_command}"
        "cd {wildcards.directory};"
        "ln -s {input.reads} reads4racon.fastq.gz;"
        "{params.racon_env}/scripts/racon_wrapper.py -u {params.opts} -t {threads} reads4racon.fastq.gz {input.mapping} {input.assembly} > {output.polished};"
        "rm reads4racon.fastq.gz;"
      )
      if os.stat(output.polished).st_size == 0:
        shell (
          "rm {output.polished};"
        )
if mr > 0:
  rule medaka:
    input:
      assembly = lambda wildcards: medaka_in[wildcards.directory + "/rmp/" + wildcards.base + ".medaka" + wildcards.param + ".fasta"],
      reads = ont_reads, 
      mapping = lambda wildcards: medaka_map[wildcards.directory + "/rmp/" + wildcards.base + ".medaka" + wildcards.param + ".fasta"], 
    output:
      polished = "{directory}/rmp/{base}.medaka{param}.fasta"
    wildcard_constraints:
      param = "\d+"
    params:
      medaka_env = config["Medaka"]["environment"],
      medaka_dir = config["Medaka"]["Working_dir"],
      model = config["Medaka"]["Model"],
      consensus_opts = medaka_consensus_opts
    log:
      "{directory}/rmp/logs/" + str(date) + "{base}.medaka{param}.out",
      "{directory}/rmp/logs/" + str(date) + "{base}.medaka{param}.err",
    threads: config["Parameters"]["medaka_cores"]
    run:
      if os.path.exists(wildcards.directory + "/rmp/" + wildcards.base + ".medaka" + wildcards.param + ".hdf"):
        shell(
          "echo '{wildcards.directory}/rmp/{wildcards.base}.medaka{wildcards.param}.hdf already exists, deleting';"
          "rm {wildcards.directory}/rmp/{wildcards.base}.medaka{wildcards.param}.hdf;"
        )        
      shell(
        "module purge;"
        "source ~jgomez/init_shell.sh;"
        "cd {params.medaka_dir};"
        "conda activate {params.medaka_env};"       
        "medaka consensus {input.mapping} {wildcards.directory}/rmp/{wildcards.base}.medaka{wildcards.param}.hdf --threads {medaka_threads} --model {params.model} {params.consensus_opts};"
        "echo 'Consensus has been called, generating polished output now';"
        "medaka stitch --threads {threads} {wildcards.directory}/rmp/{wildcards.base}.medaka{wildcards.param}.hdf {input.assembly} {output.polished};"
        "conda deactivate;"
      )
  

if config["Pilon"]["chunks"] > 1 and pr > 1: 
  rule split_pilon:
    input:
      assembly = lambda wildcards: pilon_in[wildcards.directory + "/rmp/" + wildcards.base + ".pilon" + wildcards.param + ".fasta"],
      mapping = lambda wildcards: pilon_map[wildcards.directory + "/rmp/" + wildcards.base + ".pilon" + wildcards.param + ".fasta"], 
    output:
      bed = expand("{{directory}}/rmp/{{base}}.split{{param}}/chunk{c}.bed", c=range(0,config["Pilon"]["chunks"])),
      bam = expand("{{directory}}/rmp/{{base}}.split{{param}}/chunk{c}.bam", c=range(0,config["Pilon"]["chunks"])),
      bai = expand("{{directory}}/rmp/{{base}}.split{{param}}/chunk{c}.bam.bai", c=range(0,config["Pilon"]["chunks"]))
    wildcard_constraints:
      param = "\d+"
    params:
      chunks = config["Pilon"]["chunks"],
      splitdir = "{directory}/rmp/{base}.split{param}/"
    log:
      "{directory}/rmp/logs/" + str(date) + ".{base}.split_pilon{param}.out",
      "{directory}/rmp/logs/" + str(date) + ".{base}.split_pilon{param}.err",
    threads: config["Parameters"]["BWA_cores"]
    shell:
      "module load SAMTOOLS/1.12 bwa java/1.8.0u31 PILON/1.21;"
      "mkdir -p {params.splitdir};cd {params.splitdir};"
      "{scripts_dir}fastalength {input.assembly} | sort -k1,1nr > assembly.len;"
      "{scripts_dir}split_bam.py assembly.len {input.mapping} {params.chunks} {threads};" 

if pr > 0:
  rule pilon:
    input:
      assembly = lambda wildcards: pilon_in[wildcards.directory + "/rmp/" + wildcards.base + ".pilon" + wildcards.param + ".fasta"], 
      alignment = "{directory}/rmp/{base}.split{param}/chunk{c}.bam" if config["Pilon"]["chunks"] > 1 else lambda wildcards: pilon_map[wildcards.directory + "/rmp/" + wildcards.base + ".pilon" + wildcards.param + ".fasta"], 
    output:
      polished = "{directory}/rmp/{base}.split{param}/chunk{c}.polished.fasta" if config["Pilon"]["chunks"] > 1 else "{directory}/rmp/{base}.pilon{param}.fasta",
      changes =  "{directory}/rmp/{base}.split{param}/chunk{c}.polished.changes" if config["Pilon"]["chunks"] > 1 else "{directory}/rmp/{base}.pilon{param}.changes",
    params:
      path = config["Pilon"]["path"],
      opts = config["Pilon"]["options"],
      java_opts = java_opt,
      dir = "{directory}/rmp/{base}.split{param}/" if config["Pilon"]["chunks"] > 1 else "{directory}/rmp/"
    log:
      "{directory}/rmp/logs/" + str(date) + ".{base}.pilon{param}.chunk{c}.out" if config["Pilon"]["chunks"] > 1 else "{directory}/rmp/logs/" + str(date) + ".{base}.pilon{param}.out",
      "{directory}/rmp/logs/" + str(date) + ".{base}.pilon{param}.chunk{c}.err" if config["Pilon"]["chunks"] > 1 else  "{directory}/rmp/logs/" + str(date) + ".{base}.pilon{param}.err"
    threads: config["Parameters"]["pilon_cores"]
    run:
      basename = os.path.splitext(os.path.basename(output.polished))[0]
      shell(
        "module purge;"
        "module load SAMTOOLS/1.12 bwa java/1.8.0u31 PILON/1.21;"
        "echo 'Running Pilon on {input.assembly} with {input.alignment}';"
        "cd {params.dir};"
        "java {params.java_opts} -jar {params.path} --genome {input.assembly} --frags {input.alignment} {params.opts} --threads {threads} --output {basename};"
      )

if config["Pilon"]["chunks"] > 1 and pr > 1:
  rule join_pilon:
    input:
      chunk_polished =expand("{{directory}}/rmp/{{base}}.split{{param}}/chunk{c}.polished.fasta", c=range(0,config["Pilon"]["chunks"])),
      chunk_changes = expand("{{directory}}/rmp/{{base}}.split{{param}}/chunk{c}.polished.changes", c=range(0,config["Pilon"]["chunks"]))
    output:
      polished = "{directory}/rmp/{base}.pilon{param}.fasta",
      changed = "{directory}/rmp/{base}.pilon{param}.changes",   
    params:
      splitdir ="{directory}/rmp/{base}.split{param}/",
      chunks = config["Pilon"]["chunks"]
    log:
      "{directory}/rmp/logs/" + str(date) + ".{base}.concat_pilon{param}.out",
      "{directory}/rmp/logs/" + str(date) + ".{base}.concat_pilon{param}.err",
    threads: 1
    run:
      shell (
        "{scripts_dir}/concat_pilon.py {params.splitdir} {params.chunks} > {output.polished};"
        "cat {input.chunk_changes} > {output.changed};"
      )
      if keepfiles == False:
        shell(
          "echo 'Deleting {params.splitdir}';"
          "rm -r {params.splitdir};" 
        )   

if nor > 0:
  rule nextpolish_lr:
    input:
      genome = lambda wildcards: nor_in[wildcards.directory + "nextpolish/" + wildcards.base + ".nextpolish_ont" + wildcards.param + ".fasta"],
      bam =   lambda wildcards: nor_map[wildcards.directory + "nextpolish/" + wildcards.base + ".nextpolish_ont" + wildcards.param + ".fasta"]
    output:
      polished = "{directory}nextpolish/{base}.nextpolish_ont{param}.fasta"
    params:
      lrtype = nextpolish_lrtype
    wildcard_constraints:
      param="\d+"
    log:
      "{directory}nextpolish/logs/" + str(date) + ".{base}.nextpolish_lr{param}.out",
      "{directory}nextpolish/logs/" + str(date) + ".{base}.nextpolish_lr{param}.err",
    threads:  config["Parameters"]["nextpolish_cores"]
    shell:
      "module purge; module unload intel; module load NEXTPOLISH/1.3.1;"
      "cd {wildcards.directory}nextpolish;"
      "echo {input.bam} > lgs.fofn;"
      "python /apps/NEXTPOLISH/1.3.1/lib/nextpolish2.py -g {input.genome} -p {threads} -l lgs.fofn -r {params.lrtype} > {output.polished};"

rule nextpolish_sr:
  input: 
    genome = lambda wildcards: nir_in[wildcards.directory + "nextpolish/" + wildcards.base + ".nextpolish_ill" + wildcards.param ],
    bam =   lambda wildcards: nir_map[wildcards.directory + "nextpolish/" + wildcards.base + ".nextpolish_ill" + wildcards.param ]
  output:
    polished = "{directory}nextpolish/{base}.nextpolish_ill{param}"
  params:
    task = lambda wildcards: tasks[wildcards.directory + "nextpolish/" + wildcards.base + ".nextpolish_ill" + wildcards.param]
  log:
    "{directory}nextpolish/logs/" + str(date) + ".{base}.nextpolish_sr{param}.out",
    "{directory}nextpolish/logs/" + str(date) + ".{base}.nextpolish_sr{param}.err",
  threads:  config["Parameters"]["nextpolish_cores"]
  shell:
    "module purge;module unload intel; module load NEXTPOLISH/1.3.1;"
    "cd {wildcards.directory}nextpolish;"
    "python /apps/NEXTPOLISH/1.3.1/lib/nextpolish1.py -g {input.genome}  -p {threads} -s {input.bam} -t {params.task} > {output.polished};"

if hr > 0:
  rule hypo:
    input:
      genome = lambda wildcards: hypo_in[wildcards.directory + "hypo/" + wildcards.base + ".hypo" + wildcards.param + ".fasta"],
      lr_bam = lambda wildcards: hypo_lrmap[wildcards.directory + "hypo/" + wildcards.base + ".hypo" + wildcards.param + ".fasta"],
      sr_bam = lambda wildcards: hypo_srmap[wildcards.directory + "hypo/" + wildcards.base + ".hypo" + wildcards.param + ".fasta"],
      reads_file = [pe1_reads, pe2_reads] if r10X_reads == None else [r10X_reads]
  #  pe1 = pe1_reads,
   # pe2 = pe2_reads
    output:
      polished ="{directory}hypo/{base}.hypo{param}.fasta"
    params:
      hypo_env = config["Hypo"]["environment"],
      genome_size = config["Parameters"]["genome_size"],
      cov = config["Hypo"]["illumina coverage"],
      proc = config["Hypo"]["processes"],
      opts = additional_hypo_opts
    wildcard_constraints:
      param="\d+"
    log:
      "{directory}hypo/logs/" + str(date) + ".{base}.hypo{param}.out",
      "{directory}hypo/logs/" + str(date) + ".{base}.hypo{param}.err",
    threads: config["Parameters"]["hypo_cores"]
    run:
      if (params.cov == 0):
        cov_out = subprocess.check_output(scripts_dir + "get_cov.py " + input.genome + " " + input.sr_bam, shell=True)
        coverage = str(cov_out).split("'")[1].replace('\\n','')
      else:
        coverage = params.cov
      shell(
        "echo 'Illumina coverage is:';"
        "echo {coverage};"
      )
      reads = open(wildcards.directory + "hypo/short_reads.list.txt", 'w')
    #  reads_file = input.reads.split(' ')
      for f in input.reads_file:
        reads.write(f + "\n")
      reads.close()
      shell(
        "cd {wildcards.directory}hypo;"
        "module purge;"
        "source ~jgomez/init_shell.sh;"
        "conda activate {params.hypo_env};"
        "echo 'Running: hypo -r @short_reads.list.txt -d {input.genome} -b {input.sr_bam} -c {coverage} -s {params.genome_size} -B {input.lr_bam} -t {threads} -o {output.polished} -p {params.proc} {params.opts}';"
        "hypo -r @short_reads.list.txt -d {input.genome} -b {input.sr_bam} -c {coverage} -s {params.genome_size} -B {input.lr_bam} -t {threads} -o {output.polished} -p {params.proc} {params.opts};"
        "conda deactivate;"
      )

#4- Run evaluations
if config["Finalize"]["BUSCO lineage"] != None:
  use rule run_busco from eval_workflow with:
    input:
      assembly = lambda wildcards: busco_in[eval_dir + wildcards.dir + "/" + wildcards.buscobase],
      lineage = config["Finalize"]["BUSCO lineage"]
    output:
      summary = eval_dir + "{dir}/busco/{buscobase}.short_summary.txt",
      full = eval_dir + "{dir}/busco/{buscobase}.full_table.tsv",
    params:
      busco_env = config["Finalize"]["BUSCO environment"],
      out_path = lambda wildcards: evals_dir[eval_dir + wildcards.dir +"/" + wildcards.buscobase],
      odb = os.path.basename(config["Finalize"]["BUSCO lineage"]),
      buscobase = lambda wildcards:  wildcards.buscobase
    log:
      eval_dir + "{dir}/logs/" + str(date) + ".busco.{buscobase}.out",
      eval_dir + "{dir}/logs/" + str(date) + ".busco.{buscobase}.err",
    threads: config["Parameters"]["busco_cores"]

if  config["Finalize"]["Merqury db"] != None:
 use rule run_merqury from eval_workflow with:
  input:
    meryl_db = config["Finalize"]["Merqury db"],
    assembly = lambda wildcards: merq_in[eval_dir + wildcards.dir +"/" + wildcards.merqbase],
  output:
    completeness = eval_dir + "{dir}/merqury/{merqbase}/completeness.stats",
    qv = eval_dir + "{dir}/merqury/{merqbase}/{merqbase}.qv",
    hist = eval_dir + "{dir}/merqury/{merqbase}/{merqbase}.{merqbase}.spectra-cn.hist",
    plots = [eval_dir + "{dir}/merqury/{merqbase}/{merqbase}.spectra-cn.ln.png", eval_dir + "{dir}/merqury/{merqbase}/{merqbase}.spectra-cn.fl.png", eval_dir + "{dir}/merqury/{merqbase}/{merqbase}.spectra-cn.st.png"]
  params:
    out_pref = "{merqbase}",
    conda_env = config["Finalize"]["Merqury environment"],
    directory= eval_dir + "{dir}/merqury/{merqbase}",
  log:
    eval_dir + "{dir}/logs/" + str(date) + ".merqury.{merqbase}.out",
    eval_dir + "{dir}/logs/" + str(date) + ".merqury.{merqbase}.err"

rule finalize:
  input:
    assembly = assemblies4Busco,
    buscos = BuscoSummaries,
    merqs= MerqurySummaries,
    qv = MerquryQV 
  output:
    output = config["Outputs"]["stats_out"]
  params:
  log:
    logs_dir + str(date) + ".finalize.out",
    logs_dir + str(date) + ".finalize.err"
  threads: 1
  run:
    shell(
      '''
      echo "assembly\tN50\tL50\ttotal_len\ttotal_seq\tBUSCO\tQV\tMerqury_completeness" > {output.output}
      '''
    )
    for genome in input.assembly:
      outbase = os.path.splitext(os.path.basename(genome))[0]  
      dir = eval_dir + os.path.basename(os.path.dirname(os.path.dirname(genome))) + "/"
      if not os.path.exists(dir + "stats"):
        os.makedirs(dir + "stats")
      shell(
        "{scripts_dir}fastalength {genome} | {scripts_dir}Nseries.pl > {dir}stats/{outbase}.nseries.txt;"
        "{scripts_dir}get_final_tbl.py {dir} {outbase} >> {output.output};"
      )
    if keepfiles == False:
      shell (
        "echo 'Pipeline has been completed succesffully, we are now going to delete temporary files:';"
      )
     # for f in working_dir + "tmp_meryl", hypo_dir + "aux", working_dir + "mappings":
      for i in inputs:
        rundir = inputs[i]
        for f in rundir + "tmp_meryl", rundir + "hypo/aux", rundir + "mappings":
          if (os.path.exists(f)):
            shell(
              "echo 'Deleting {f}';"
              "rm -r {f};"
            )
    shell(
      "echo 'Pipeline completed, bye';"
    )
