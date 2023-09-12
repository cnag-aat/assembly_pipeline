from datetime import datetime
import os
import re
import subprocess

module eval_workflow:
  snakefile: "../modules/evaluate_assemblies.rules.smk"

##0. Define path for files and variables

inputs = config["Inputs"]["Assemblies for polishing"]

assemblies = []
if config["Finalize"]["final Evaluations"] == True:
  for i in inputs:
    assemblies.append(i)
  for i in polished:
    assemblies.append(i)   
  for i in postpolish:
    assemblies.append(i)
  for i in assembled:
    assemblies.append(i)
  for i in pretext_in:
    assemblies.append(i)

in_files = {}
evals_dir = {}
BuscoSummaries = []
StatsFiles = []
MerqurySummaries = []
MerquryQV = []
MerquryDups = []
lrtype = config["Parameters"]["lr_type"]
stats_loc = {}

for file in assemblies:
  ass_base = os.path.splitext(os.path.basename(file))[0]
  basedirname = os.path.basename(os.path.dirname(file))
  if basedirname == "hypo" or basedirname == "rmp" or basedirname == "nextpolish":
    basedirname = os.path.basename(os.path.dirname(os.path.dirname(file)))
  evalassdir =  eval_dir + basedirname + "/"
  if not os.path.exists(evalassdir + "logs/"):
    os.makedirs(evalassdir + "logs/")
  if not os.path.exists(evalassdir + "stats"):
    os.makedirs(evalassdir + "stats")
  buscodir = evalassdir + "busco/"
  if not os.path.exists(buscodir) and config["Finalize"]["BUSCO lineage"]:
    os.makedirs(buscodir)
  merqdir = evalassdir + "merqury/" + ass_base 
  if not os.path.exists(merqdir) and config["Finalize"]["Merqury db"]:
    os.makedirs(merqdir)
  in_files[evalassdir + ass_base] = file
  evals_dir[evalassdir + ass_base] = evalassdir

  if config["Finalize"]["BUSCO lineage"]:
    BuscoSummaries.append(buscodir + ass_base + ".short_summary.txt")

  if config["Finalize"]["Merqury db"]:
    MerqurySummaries.append(merqdir + "/" + ass_base + ".completeness.stats")
    MerquryQV.append(merqdir + "/" + ass_base + ".qv")
    MerquryDups.append(merqdir + "/" + ass_base + ".false_duplications.txt")

  StatsFiles.append(evalassdir + "stats/" + ass_base + ".stats.txt")
  stats_loc[ass_base] = evalassdir + "stats/" + ass_base + ".gaps.txt"

#1- Perform alignments

if len(bwa) > 0:
  use rule align_illumina from eval_workflow with:
    input:
      genome = lambda wildcards: bwa[wildcards.name],
      reads = [pe1_reads, pe2_reads] if r10X_reads == None else r10X_reads
    output:
      mapping = "{directory}/mappings/{name}_bwa.bam",
      stats = report("{directory}/{name}_bwa.stats.txt",
              caption = "../report/bwa.rst",
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
      align_opts = lambda wildcards:"ax" if wildcards.ext == "minimap2.bam" else "x",
      type = lambda wildcards:"map-ont" if lrtype == "nano-raw" else "map-hifi",
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

if len(hic_assemblies) > 0:

  use rule align_hic from eval_workflow with:
    input:
        ass = lambda wildcards: hic_assemblies[wildcards.name],
        bwt = lambda wildcards: hic_assemblies[wildcards.name] + ".bwt",
        read1 = hic_pe1,
        read2 = hic_pe2
    output:
        mapped = "{directory}/mappings/{name}.full_hic.bam",
        unmapped = "{directory}/mappings/{name}.unmapped_hic.bam"
    params:
        name = '{name}',
        outd = '{directory}/mappings/',
        options = config['HiC']['align_opts']
    log:
        "{directory}/logs/" + str(date) + ".j%j.rule_mapping.{name}.out",
        "{directory}/logs/" + str(date) + ".j%j.rule_mapping.{name}.err"
    benchmark:
        "{directory}/logs/" + str(date) + ".rule_mapping.benchmark.{name}.txt"
    threads: config["Parameters"]["BWA_cores"]

  use rule pairtools_processing from eval_workflow with:
    input: 
      mapped = lambda wildcards: hic_bams[wildcards.name],
      alength = lambda wildcards: asslength[wildcards.name]
    output:
      stats = "{directory}/pairtools_out/stats.mq{mq}.{name}.txt",
      mappedptsort = "{directory}/pairtools_out/mapped.PT.mq{mq}.{name}.name_sorted.bam",
      mappedpt = "{directory}/pairtools_out/mapped.PT.mq{mq}.{name}.bam",
    wildcard_constraints:
      mq="\d+",
    params:
      scripts_dir = scripts_dir,
      mq = '{mq}',
      name = '{name}',
      outd = '{directory}/pairtools_out',
      tmpd = '{directory}/{name}_mq{mq}_tmp',
      rmcmd = "rm -r {directory}/{name}_mq{mq}_tmp;" if keepfiles == False
               else ""
    log:
        "{directory}/logs/" + str(date) + ".j%j.rule_pairtools.mq{mq}.{name}.out",
        "{directory}/logs/" + str(date) + ".j%j.rule_pairtools.mq{mq}.{name}.err"
    benchmark:
        "{directory}/logs/" + str(date) + ".rule_pairtools.benchmark.mq{mq}.{name}.txt"
    threads: config['Parameters']['pairtools_cores']

  use rule qc_statistics from eval_workflow with: 
    input:
      statmq = "{directory}/pairtools_out/stats.mq{mq}.{name}.txt",
      mapbam = "{directory}/pairtools_out/mapped.PT.mq{mq}.{name}.bam"
    output:
      libstats = "{directory}/pairtools_out/HiC_QC_LibraryStats_mq{mq}.{name}.txt",
      pslibstats = "{directory}/pairtools_out/HiC_QC_LibraryStats_extrapolated_mq{mq}.{name}.txt" if not config['HiC']['deepseq'] else []
    wildcard_constraints:
      mq="\d+",
    params:
      scripts_dir = scripts_dir,
      outd = '{directory}/out',
      assemblylength = config['HiC']['qc_assemblylen'],
      deepseq = config['HiC']['deepseq'],
    log:
        "{directory}/logs/" + str(date) + ".j%j.rule_qc_stats.mq{mq}.{name}.out",
        "{directory}/logs/" + str(date) + ".j%j.rule_qc_stats.mq{mq}.{name}.err"
    benchmark:
        "{directory}/logs/" + str(date) + ".rule_qc_stats.benchmark.mq{mq}.{name}.txt"
    threads: 2
  if config['HiC']['get_pretext']:
    use rule get_extension_gaps from eval_workflow with:
      input:
        sla = lambda wildcards: hic_assemblies[wildcards.name],
        gaps_bed = lambda wildcards: stats_loc[wildcards.name]
      output:
        gaps = "{directory}/{name}.gaps.bg"
      params:
        scripts_dir = scripts_dir,
      log:
        "{directory}/logs/" + str(date) + ".j%j.rule_gaps.{name}.out",
        "{directory}/logs/" + str(date) + ".j%j.rule_gaps.{name}.err"
      benchmark:
        "{directory}/logs/" + str(date) + ".rule_gaps.{name}.benchmark.txt"
      threads:  1

    use rule get_extension_telomeres from eval_workflow with:
      input:
        sla = lambda wildcards: hic_assemblies[wildcards.name]
      output:
        tel = "{directory}/{name}.telomeres.bg"
      params:
        outd = "{directory}"
      log:
        "{directory}/logs/" + str(date) + ".j%j.rule_telomeres.{name}.out",
        "{directory}/logs/" + str(date) + ".j%j.rule_telomeres.{name}.err"
      benchmark:
        "{directory}/logs/" + str(date) + ".rule_telomeres.{name}.benchmark.txt"
      threads:  2

    if len(minimap2) > 0:
      use rule get_extension_ont from eval_workflow with:
        input:
          sbam = lambda wildcards: pretext_lrmap[wildcards.name],
        output: 
          ontcov = "{directory}/{name}.ONTcoverage.bg"
        log:
          "{directory}/logs/" + str(date) + ".j%j.rule_ont_bed.{name}.out",
          "{directory}/logs/" + str(date) + ".j%j.rule_ont_bed.{name}.err"
        benchmark:
          "{directory}/logs/" + str(date) + ".rule_ont_bed.{name}.benchmark.txt"
        threads:  1

  if not config['HiC']['deepseq']:  #run for qc as well
    use rule read_screening from eval_workflow with:
      input:
        umapped = "{directory}/mappings/{name}.unmapped_hic.bam",
      output:
        ufasta = "{directory}/blast/unmapped_hic.{name}.{readsblast}.fasta",
        subufasta = "{directory}/blast/unmapped_hic.{name}.{readsblast}_reads.fasta",
        blastoutbscore = "{directory}/blast/unmapped_hic.{name}.{readsblast}_unmapped_reads_vs_nt_25cul1_1e25.megablast.sorted_by_bitscore.out",
        blastoutbscorethits = "{directory}/blast/unmapped_hic.{name}.{readsblast}_unmapped_reads_vs_nt_25cul1_1e25.megablast.sorted_by_bitscore.tophits",
        blastoutorganisms = "{directory}/blast/unmapped_hic.{name}.{readsblast}_unmapped_reads_vs_nt_25cul1_1e25.megablast.organisms.txt" 
      params:
        scripts_dir = scripts_dir,
        outd = "{directory}/blast/",
        readforblast = config['HiC']['reads_for_blast'],
        blastdb = config['HiC']['blastdb']
      log:
        "{directory}/logs/" + str(date) + ".j%j.rule_screening.{name}.{readsblast}_reads.out",
        "{directory}/logs/" + str(date) + ".j%j.rule_screening.{name}.{readsblast}_reads.err"
      benchmark:
        "{directory}/logs/" + str(date) + ".rule_screening.{name}.{readsblast}_reads.benchmark.txt"
      threads:  config['HiC']['blast_cores']

#2- Run evaluations

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
    eval_dir + "{dir}/logs/" + str(date) + ".get_stats.{buscobase}.benchmark.out",
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
      rmcmd =  lambda wildcards: "echo 'Removing BUSCO run dir'; rm -r " + \
             evals_dir[eval_dir + wildcards.dir +"/" + wildcards.buscobase] + wildcards.buscobase + ";" \
             if keepfiles == False else "" 
    log:
      eval_dir + "{dir}/logs/" + str(date) + ".j%j.busco.{buscobase}.out",
      eval_dir + "{dir}/logs/" + str(date) + ".j%j.busco.{buscobase}.err",
    benchmark:
      eval_dir + "{dir}/logs/" + str(date) + ".busco.{buscobase}.benchmark.txt"
    conda:
      '../envs/busco5.4.0.yaml'
    threads: config["Parameters"]["busco_cores"]

if  config["Finalize"]["Merqury db"] != None:
  use rule run_merqury from eval_workflow with:
    input:
      meryl_db = config["Finalize"]["Merqury db"],
      assembly = lambda wildcards: in_files[eval_dir + wildcards.dir +"/" + wildcards.merqbase],
    output:
      completeness = report (eval_dir + "{dir}/merqury/{merqbase}/{merqbase}.completeness.stats",
                     caption="../report/merqury.rst",
                     category = "Evaluate assemblies",
                     subcategory = "{dir}"),
      qv = report(eval_dir + "{dir}/merqury/{merqbase}/{merqbase}.qv",
           caption="../report/merqury.rst",
           category = "Evaluate assemblies",
           subcategory = "{dir}"),
      hist = eval_dir + "{dir}/merqury/{merqbase}/{merqbase}.{merqbase}.spectra-cn.hist",
      false_dups = report(eval_dir + "{dir}/merqury/{merqbase}/{merqbase}.false_duplications.txt",
           caption="../report/merqury.rst",
           category = "Evaluate assemblies",
           subcategory = "{dir}"),
      plots = report ([eval_dir + "{dir}/merqury/{merqbase}/{merqbase}.spectra-cn.ln.png", eval_dir + "{dir}/merqury/{merqbase}/{merqbase}.spectra-cn.fl.png", eval_dir + "{dir}/merqury/{merqbase}/{merqbase}.spectra-cn.st.png"],
              caption="../report/merqury.rst",
              category = "Evaluate assemblies",
              subcategory = "{dir}"),
    params:
      out_pref = "{merqbase}",
      directory= eval_dir + "{dir}/merqury/{merqbase}",
    log:
      eval_dir + "{dir}/logs/" + str(date) + ".j%j.merqury.{merqbase}.out",
      eval_dir + "{dir}/logs/" + str(date) + ".j%j.merqury.{merqbase}.err"
    benchmark:
      eval_dir + "{dir}/logs/" + str(date) + ".merqury.{merqbase}.benchmark.txt"
    conda:
      "../envs/merqury1.3.yaml"
    threads: 
      config["Finalize"]["Meryl threads"]

#Run final job of the pipeline
rmcmd = ""
if keepfiles == False:
  rmcmd = "echo 'Pipeline has been completed succesffully, we are now going to delete temporary files:';"
  if config["Inputs"]["processed_illumina"] != None:
    t = config["Inputs"]["processed_illumina"]
    illumina_list = config["Wildcards"]["illumina_wildcards"].split(',')
    if os.path.exists(config["Inputs"]["processed_illumina"] + illumina_list[0] + ".1_val_1.fq.gz"):
      rmcmd += "echo 'Deleting fastqs in " + t + "'; rm " + t + "*.gz;"
    if config["Finalize"]["Merqury db"] != None:
      t = os.path.dirname(config["Finalize"]["Merqury db"]) + "/tmp_meryl/"
      if (os.path.exists(t)):
        rmcmd += "echo 'Deleting " + t + "'; rm -r " + t + ";"
  for i in inputs:
    rundir = inputs[i]
    rmcmd += "echo 'Deleting mappings in " + rundir + "mappings'; rm -r " + rundir +  "mappings;"
    if os.path.exists(rundir + "hypo/aux"):
      rmcmd += "echo 'Deleting " + rundir + "hypo/aux;'; rm -r " + rundir + "/hypo/aux;"

use rule finalize from eval_workflow with:
  input:
    assembly = assemblies,
    buscos = BuscoSummaries,
    stats= StatsFiles,
    merqs=MerquryQV,
    pretext = pretext_files
  output:
    output = report(config["Outputs"]["stats_out"],
             caption="../report/final.rst",
             category = "Final Table"),
  params:
    scripts_dir = scripts_dir,
    rmcmd = rmcmd,
  log:
    logs_dir + str(date) + ".j%j.finalize.out",
    logs_dir + str(date) + ".j%j.finalize.err"
  benchmark:
    logs_dir + str(date) + ".finalize.benchmark.txt"
  conda:
    '../envs/ass_base.yaml'
  threads: 1
  
# use rule get_report from eval_workflow with:
#   input:
#     stats = config["Outputs"]["stats_out"],
#     config = config["Parameters"]["configFile"]
#   output:
#     report = base + ".report.zip"
#   log:
#     logs_dir + str(date) + ".j%j.get_report.out",
#     logs_dir + str(date) + ".j%j.get_report.err"
#   benchmark:
#     logs_dir + str(date) + ".get_report.benchmark.txt"
#   conda:
#     '../envs/ass_base.yaml'
#   threads: 1
