from datetime import datetime
import re
import os

date = datetime.now().strftime('%Y%m%d.%H%M%S')
keepfiles = False

scripts_dir = os.path.dirname(sys.argv[0]) + "/../scripts/"

shell.prefix( "export PATH=" + scripts_dir + ":$PATH;")

if not os.path.exists("logs"):
  os.makedirs("logs")

rule scaffolding_10X:
  input:
    assembly = "flye_filtlong.assembly.hypo2.purged.fa",
    reads = "BALDOLAU_08.barcoded.fastq.gz"
  output:
    scaffolded = "flye_filtlong.assembly.hypo2.purged.10X.scaffolds.fa"
  params:
    opts = "",
    dir = os.getcwd(),
    rmcmd = "",
    scripts = "../scripts/"
  threads: 24
  conda:
    "../envs/tigmint1.2.6_arcs1.2.4.yaml"
  shell:
      "base_ass=$(basename {input.assembly} .fa);"
      "base_reads=$(basename {input.reads} .fastq.gz);"
      "mkdir -p {params.dir}/tigmint_with_ARKS/tmp;"
      "cd {params.dir}/tigmint_with_ARKS;"
      "ln -s {input.assembly} $base_ass.fa;"
      "ln -s {input.reads} $base_reads.fq.gz;"
      "TMPDIR={params.dir}/tigmint_with_ARKS/tmp;"
      "echo \"Tmpdir is $TMPDIR\";"
      "tigmint-make arcs draft=$base_ass reads=$base_reads {params.opts} t={threads};"
      "sed \'s/,/ /\' {params.dir}/tigmint_with_ARKS/$base_ass_*.scaffolds.fa | {params.scripts_dir}FastaToTbl | {params.scripts_dir}TblToFasta > {output.scaffolded};"
      "{params.rmcmd}"
