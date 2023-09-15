from datetime import datetime
import os
import re
import subprocess

module purging_workflow:
  snakefile: "../modules/purging.rules.smk"
module scaffolding_workflow:
  snakefile: "../modules/scaffolding_rules.smk"
module hic_workflow:
    snakefile: "../modules/hic.rules.smk"

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

hic_assemblies = {}
asslength = {}
pretext_files = []
pretext_in = []
tpf_files = []

if config['Inputs']['HiC_dir'] and config["Wildcards"]["HiC_wildcards"]:
  hic_out_dir = {}
  if config['HiC']['deepseq'] == False:
    i = config['HiC']['assembly_qc']
    name = os.path.splitext(os.path.basename(i))[0]
    hic_assemblies[name] = i
    hic_out_dir[name] = config['Outputs']['hic_qc_dir']
  else:
    if config["Parameters"]["run_tigmint"] == False and config["Parameters"]["run_purgedups"] == False:
      for i in config["Inputs"]["Assemblies for postpolishing"]:
        name = os.path.splitext(os.path.basename(i))[0]
        hic_assemblies[name] = i 
    elif config["Parameters"]["run_tigmint"] == True:
      for i in postpolish:
        if "10X" in i:
          name = os.path.splitext(os.path.basename(i))[0]
          hic_assemblies[name] = i
    elif config["Parameters"]["run_purgedups"] == True:
      for i in postpolish:
        if "purgedups" in i:
          name = os.path.splitext(os.path.basename(i))[0]
          hic_assemblies[name] = i 

    for i in config["Inputs"]["Assemblies for postpolishing"]:
      step = config["Inputs"]["Assemblies for postpolishing"][i]
      base = os.path.splitext(os.path.basename(i))[0]
      assembly = i
      pstep = step.split('_')[0]
      nstep = pstep.replace('s','')
      if config["Parameters"]["run_tigmint"] == True and config["Parameters"]["run_purgedups"] == True:
        cstep = float(nstep) + 2
        step = "s0" + str(cstep) + "_" + step.split('_')[0].replace('s','p')
        assembly =  working_dir + config["Inputs"]["Assemblies for postpolishing"][i] + "_run_purgedups/" + base + ".purged.10X.scaffolds.fa"
      elif config["Parameters"]["run_purgedups"] == True:
        cstep = float(nstep) + 1
        step = "s0" + str(cstep) + "_" + step.split('_')[0].replace('s','p')
        assembly =  working_dir + config["Inputs"]["Assemblies for postpolishing"][i] + "_run_purgedups/" + base + ".purged.fa"
      elif config["Parameters"]["run_tigmint"] == True:
        cstep = float(nstep) + 1
        step = "s0" + str(cstep) + "_" + step.split('_')[0].replace('s','p')
        assembly =  working_dir + config["Inputs"]["Assemblies for postpolishing"][i] + "_run_purgedups/" + base + ".10X.scaffolds.fa"
      name = os.path.splitext(os.path.basename(assembly))[0]
      hic_out_dir[name] = working_dir + step + "_HiC_scaffolding/"
      postpolish.append(hic_out_dir[name] + name + ".yahs_scaffolds_final.fa")
      tpf_files.append(hic_out_dir[name] + name + ".yahs_scaffolds_final.fa.tpf")

  hic_bams = {}
  pretext_lrmap = {}
  for i in hic_out_dir:
    name = i
    hic_bams[name] = hic_out_dir[i] + "mappings/" + name + ".full_hic.bam"
    pretext_lrmap[name] = hic_out_dir[i] + "mappings/" + name + "_minimap2.bam"
    if not os.path.exists(hic_out_dir[i] + "/logs"):
      os.makedirs(hic_out_dir[i] + "/logs")
    if config['Parameters']['run_yahs']:
      name = i + ".yahs_scaffolds_final"
      hic_assemblies[name] = hic_out_dir[i] + name + ".fa"
      hic_bams[name] = hic_out_dir[i] + "mappings/" + name + ".full_hic.bam"
      pretext_lrmap[name] = hic_out_dir[i] + "mappings/" + name + "_minimap2.bam"

  for i in hic_assemblies:
    asslength[i] = os.path.dirname(hic_assemblies[i]) + "/"+ i + ".genome"
    pretext_in.append(os.path.dirname(hic_assemblies[i]) + "/" + i + ".fa")
    minimap2[i] = os.path.dirname(hic_assemblies[i]) + "/" + i + ".fa"
    for mq in config['HiC']['MQ']:
      if i in hic_out_dir:
        pretext_files.append(hic_out_dir[i] + "in_pretext/" + i + "_mq" + str(mq) + ".extensions.pretext")
        if not os.path.exists(hic_out_dir[i] + "in_pretext/logs"):
          os.makedirs(hic_out_dir[i] + "in_pretext/logs")
      else:
        pretext_files.append(os.path.dirname(hic_assemblies[i]) + "/out_pretext/" + i + "_mq" + str(mq) + ".extensions.pretext")
        if not os.path.exists(os.path.dirname(hic_assemblies[i]) + "/out_pretext/logs"):
          os.makedirs(os.path.dirname(hic_assemblies[i]) + "/out_pretext/logs")

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

if config['Inputs']['HiC_dir'] and config["Wildcards"]["HiC_wildcards"]:

  use rule assembly_prepare from hic_workflow with:
    input:
      lambda wildcards: hic_assemblies[wildcards.name]
    output:
      glen = "{dir}/{name}.genome",
      bwa = "{dir}/{name}.fa.bwt",
      faidx = "{dir}/{name}.fa.fai"   
    params:
      scripts = scripts_dir,
      workdir = lambda wildcards: os.path.dirname(hic_assemblies[wildcards.name]),
    log:
      "{dir}/logs/" + str(date) + ".j%j.rule_assembly.prepare.{name}.out",
      "{dir}/logs/"  + str(date) + ".j%j.rule_assembly.prepare.{name}.err"  
    benchmark:
      "{dir}/logs/" + str(date) + ".rule_assembly.prepare.benchmark.{name}.txt"
    threads: 2

  if config['Parameters']['run_yahs']:
    use rule run_yahs from hic_workflow with:
      input:
        mappedptsort = "{directory}/in_pretext/pairtools_out/mapped.PT.mq" + str(config['HiC']["yahs_mq"]) + ".{name}.name_sorted.bam",
        sla = lambda wildcards: hic_assemblies[wildcards.name],
        index = lambda wildcards: hic_assemblies[wildcards.name] + ".fai"
      output:
        outyahs = "{directory}/{name}.yahs_scaffolds_final.fa",
        agp = "{directory}/{name}.yahs_scaffolds_final.agp"  
      params:
        name = "{name}",
        yahsdir = "{directory}/run_yahs/",
        mq = config['HiC']["yahs_mq"],
        yahsopts = config["HiC"]['yahs_opts'],
      log:
        "{directory}/logs/" + str(date) + ".j%j.rule_yahs.{name}.out",
        "{directory}/logs/" + str(date) + ".j%j.rule_yahs.{name}.err"
      benchmark:
        "{directory}/logs/" + str(date) + ".rule_yahs.{name}.benchmark.txt"
      threads:  config['HiC']['yahs_cores']

if config['HiC']['get_pretext']:
  use rule generate_pretext from hic_workflow with:
    input:
      mapbam = "{directory}/pairtools_out/mapped.PT.mq{mq}.{name}.bam"
    output:
      pret = "{directory}/{name}_mq{mq}.pretext", 
      ptd = "{directory}/snapshots/PTdone{mq}.{name}.txt"
    wildcard_constraints:
      mq="\d+",
    params:
      scripts_dir = scripts_dir,
      outd = "{directory}",
      name = "{name}",
      mq = "{mq}",
    log:
      "{directory}/logs/" + str(date) + ".j%j.rule_pretext.{mq}.{name}.out",
      "{directory}/logs/" + str(date) + ".j%j.rule_pretext.{mq}.{name}.err"
    benchmark:
      "{directory}/logs/" + str(date) + ".rule_pretext.{mq}.{name}.txt"
    threads: config['Parameters']['pairtools_cores']

  use rule add_extensions_pretext from hic_workflow with:
    input:
      tel = "{directory}/{name}.telomeres.bg",
      gaps = "{directory}/{name}.gaps.bg",
      ontcov = lambda wildcards: "{directory}/{name}.ONTcoverage.bg" if len(minimap2) > 0 else "",
      pret = "{directory}/{name}_mq{mq}.pretext",  
    output:
      pretext = "{directory}/{name}_mq{mq}.extensions.pretext", 
    wildcard_constraints:
      mq="\d+",
    params:
      scripts_dir = scripts_dir,
      outd = "{directory}",
      mq = "{mq}",
    log:
      "{directory}/logs/" + str(date) + ".j%j.rule_add_ext.{mq}.{name}.out",
      "{directory}/logs/" +  str(date) + ".j%j.rule_add_ext.{mq}.{name}.err"
    benchmark:
      "{directory}/logs/" + str(date) + ".rule_add_ext.{mq}.{name}.benchmark.txt"
    threads:  config['Parameters']['pairtools_cores']

  use rule get_tpf from hic_workflow with:
    input:
      fasta = "{directory}/{name}.yahs_scaffolds_final.fa",
    output:
      tpf = "{directory}/{name}.yahs_scaffolds_final.fa.tpf",
    log:
      "{directory}/logs/" + str(date) + ".j%j.get_tpf.{name}.out",
      "{directory}/logs/" + str(date) + ".j%j.get_tpf.{name}.err"
    benchmark:
      "{directory}/logs/" + str(date) + ".get_tpf.{name}.benchmark.txt"
    threads:  2


