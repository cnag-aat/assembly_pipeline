from datetime import datetime
import os
import re
import subprocess

report: "../report/workflow.rst"


date = datetime.now().strftime('%Y%m%d.%H%M%S')

working_dir = config["Outputs"]["base_dir"]
#shell.prefix("TMPDIR=" + working_dir + "tmp/; echo tmpdir is set to $TMPDIR;")
scripts_dir = config["Inputs"]["scripts_dir"]
logs_dir = working_dir + "logs/"
if not os.path.exists(logs_dir):
  os.makedirs(logs_dir)

#if not os.path.exists(working_dir + "tmp/"):
#  os.makedirs(working_dir + "tmp/")

keepfiles = config["Parameters"]["keep_intermediate"]
base = config["Parameters"]["base_name"]
eval_dir = config["Outputs"]["eval_dir"]

##0. Define path for files and variables

flye_dir = config["Outputs"]["flye_dir"]
nextdenovo_dir = config["Outputs"]["nextdenovo_dir"]
flye_assembly = config["Outputs"]["flye_out"]
nextdenovo_assembly = config["Outputs"]["nextdenovo_out"]

targets = []

krakendb = ""
if config["Parameters"]["run_kraken2"] == True:
  krakendb = os.path.basename(config["Kraken2"]["database"])

if config["Parameters"]["run_kraken2"] == True:
  if config["Inputs"]["processed_illumina"]:
    targets.append(os.path.dirname(os.path.dirname(config["Inputs"]["processed_illumina"]))+ "/Kraken/" + krakendb + "/illumina_" + krakendb+".kraken2.report.txt")
  if config["Inputs"]["processed_10X"]:
    targets.append(os.path.dirname(os.path.dirname(config["Inputs"]["processed_10X"]))+ "/Kraken/" + krakendb + "/10X_" + krakendb+".kraken2.report.txt")

if config["Finalize"]["Merqury db"]:
  merqury_db = config["Finalize"]["Merqury db"]
  genomescope_dir = os.path.dirname(merqury_db) + "/genomescope2_k" + str(config["Finalize"]["Meryl K"]),
  targets.append(genomescope_dir)

if config["Parameters"]["run_flye"] == True or config["Parameters"]["run_nextdenovo"] == True:
  ONT_filtered = config["Inputs"]["ONT_filtered"]
  nanostats_dir = os.path.dirname(ONT_filtered)
  targets.append(nanostats_dir + "/nanostats/filtered_ont/NanoStats.txt")
  if config["Parameters"]["run_kraken2"] == True:
    targets.append(nanostats_dir + "/Kraken/filtered_ont/" + krakendb + "/filtlong_"+krakendb+".kraken2.report.txt")
  if config["Inputs"]["ONT_dir"] and config["Wildcards"]["ONT_wildcards"].split(','):
    targets.append(nanostats_dir + "/nanostats/raw_ont/NanoStats.txt")
    if config["Parameters"]["run_kraken2"] == True:
      targets.append(nanostats_dir + "/Kraken/raw_ont/" + krakendb + "/raw_ont_"+krakendb+".kraken2.report.txt")
    
  if config["Parameters"]["run_flye"] == True:
    targets.append(flye_assembly)
    if not os.path.exists(flye_dir + "logs"):
      os.makedirs(flye_dir + "logs")
  if config["Parameters"]["run_nextdenovo"] == True:
    targets.append(nextdenovo_assembly)
    if not os.path.exists(nextdenovo_dir + "logs"):
      os.makedirs(nextdenovo_dir + "logs")

if config["Finalize"]["final Evaluations"] == True:
  targets.append(config["Outputs"]["stats_out"])

 # targets.append(base + ".report.zip")

#1- Define rule all
rule all:
  input:
    targets
  log:
    logs_dir + str(date) + ".j%j.rule_all.out",
    logs_dir + str(date) + ".j%j.rule_all.err"

include: "../modules/preprocess_reads.smk"
include: "../modules/run_assemblies.smk"
include: "../modules/polish_assemblies.v04.smk"
include: "../modules/postpolishing.smk"
include: "../modules/assembly_evaluation.smk"
