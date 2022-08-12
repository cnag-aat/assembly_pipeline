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
    env = "/scratch/project/devel/aateam/src/TIGMINT/tigmint_conda_env/",
    opts = "",
    dir = os.getcwd()
  log:
    "logs/" + str(date) + ".10X_scaffolding.out",
    "logs/" + str(date) + ".10X_scaffolding.err"
  threads: 24
  run:
    shell(
      "module purge; source ~jgomez/init_shell.sh; conda activate {params.env};"
      "base_ass=$(basename {input.assembly} .fa);"
      "base_reads=$(basename {input.reads} .fastq.gz);"
      "mkdir -p {params.dir}/tigmint_with_ARKS;"
      "cd {params.dir}/tigmint_with_ARKS;"
      "ln -s {input.assembly} $base_ass.fa;"
      "ln -s {input.reads} $base_reads.fq.gz;"
      "arcs-make arks-with-tigmint draft=$base_ass reads=$base_reads {params.opts} t={threads};"
      "sed \'s/,/ /\' {params.dir}/tigmint_with_ARKS/$base_ass_*.scaffolds.fa |FastaToTbl | TblToFasta > {output.scaffolded};"
      "conda deactivate;"
    )
    if keepfiles == False:
      shell(
        "cd {params.dir}/tigmint_with_ARKS;"
        "base_ass=$(basename {input.assembly} .fa);"
        "base_reads=$(basename {input.reads} .fastq.gz);"
        "rm $base_ass.$base_reads.sortbx.bam;"
      )
