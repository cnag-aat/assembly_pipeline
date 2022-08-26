from datetime import datetime
import os
import re
import subprocess
import get_targets

date = datetime.now().strftime('%Y%m%d.%H%M%S')

#0. Main  
working_dir = config["Outputs"]["base_dir"]
eval_dir = config["Outputs"]["eval_dir"]
scripts_dir = config["Inputs"]["scripts_dir"]

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
bwa = {}
minimap2 = {}
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
    minimap2[basename] = assembly_in
    for i in range(1,hr):
      file_name = basename + ".hypo" + str(i) + ".fasta"
      bwa[file_name] = hypo_dir + file_name
      minimap2[file_name] = hypo_dir + file_name

for d in make_logs:
  if not os.path.exists(d + "logs"):
    os.makedirs(d+"logs")

#Perform alignments

if len(bwa) > 0:
  use rule align_illumina from eval_workflow with:
    input:
      genome = lambda wildcards: bwa[wildcards.name],
      reads = [pe1_reads, pe2_reads] 
    output:
      mapping = "{directory}/mappings/{name}_bwa.bam",
      stats = report("{directory}/{name}_bwa.stats.txt",
              caption = "bwa.rst",
              category = "Read Mapping",
              subcategory = "{name}_bwa")
    log:
      "{directory}/logs/" + str(date) + ".j%j.bwa.{name}.out",
      "{directory}/logs/" + str(date) + ".j%j.bwa.{name}.err",
    benchmark:
      "{directory}/logs/" + str(date) + ".bwa.{name}.benchmark.txt",
    conda:
      "../envs/bwa-mem2.2.1.yaml"
    threads: config["Parameters"]["BWA_cores"]

if len(minimap2) > 0:
  use rule align_ont from eval_workflow with:
    input:
      genome = lambda wildcards: minimap2[wildcards.name],
      reads = lambda wildcards: ont_reads if wildcards.ext == "minimap2.allreads.paf.gz" else ONT_filtered,
    output:
      mapping = "{directory}/mappings/{name}_{ext}"
    params:
      align_opts = lambda wildcards:"ax map-ont" if wildcards.ext == "minimap2.bam" else "x map-ont",
      tmp = "{directory}/mappings/{name}_{ext}.tmp",
      compress_cmd = lambda wildcards : "samtools view -Sb " + wildcards.directory + "/mappings/" + wildcards.name + "_" + wildcards.ext + ".tmp | " \
                     "samtools sort -@ " + str(config["Parameters"]["minimap2_cores"]) +" -o " + wildcards.directory + "/mappings/" + wildcards.name + "_" + wildcards.ext +";" +\
                     "samtools index " + wildcards.directory + "/mappings/" + wildcards.name + "_" + wildcards.ext  \
                     if wildcards.ext == "minimap2.bam" else \
                     "gzip -c " + wildcards.directory + "/mappings/" + wildcards.name + "_" + wildcards.ext + ".tmp > " + wildcards.directory + "/mappings/" + wildcards.name + "_" + wildcards.ext         
    wildcard_constraints:
      ext = "minimap2.(.+)"
    log:
      "{directory}/logs/" + str(date) + ".j%j.{ext}.{name}.out",
      "{directory}/logs/" + str(date) + ".j%j.{ext}.{name}.err",
    benchmark:
      "{directory}/logs/" + str(date) + ".{ext}.{name}.benchmark.txt",
    conda:
      "../envs/minimap2.24.yaml"
    threads: config["Parameters"]["minimap2_cores"]