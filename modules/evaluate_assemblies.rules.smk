from datetime import datetime
import re
import os

date = datetime.now().strftime('%Y%m%d.%H%M%S')
keepfiles = False
scripts_dir = os.path.dirname(sys.argv[0]) + "/../scripts/"  

shell.prefix("export PATH=" + scripts_dir + ":$PATH;")
logs_dir = "logs/"
if not os.path.exists(logs_dir):
  os.makedirs(logs_dir)

rule finalize:
  input:
    assembly = "assembly.fasta",
    buscos = "busco_short_summary.txt",
    stats = "assembly.stats.txt"
  output:
    output = "stats.txt"
  params:
  log:
    "logs/" + str(date) + ".j%j.finalize.out",
    "logs/" + str(date) + ".j%j.finalize.err"
  benchmark:
    "logs/" + str(date) + ".finalize.benchmark.txt"
  conda:
    '../envs/ass_base.yaml'
  threads: 1
  shell:
    "touch {output.output}"

rule get_stats:
  input:
    assembly = "assembly.fasta"
  output:
    nseries =  "assembly.nseries.txt",
    stats = "assembly.stats.txt",
    gaps = "assembly.gaps.txt"
  params:
    outbase = "assembly",
    scripts_dir = scripts_dir
  log:
    "logs/" + str(date) + ".j%j.assembly_stats.out",
    "logs/" + str(date) + ".j%j.assembly_stats.err"
  benchmark:
    "logs/" + str(date) + ".assembly_stats.benchmark.txt"
  conda:
    "../envs/ass_base.yaml"
  threads: 1
  shell:
    "{params.scripts_dir}fastalength {input.assembly} | {params.scripts_dir}Nseries.pl > {output.nseries};"
    "{params.scripts_dir}fasta-stats.py -f {input.assembly} -s {output.stats} -r {output.gaps};"

rule run_busco:
  input:
    assembly = os.getcwd() + "/assembly.fasta",
    lineage = "scratch/groups/assembly/shared/databases/busco_v5/busco_v5_lineages/lineages/vertebrata_odb10"
  output:
    summary = os.getcwd() + "/busco/short_summary.txt",
    full = os.getcwd() + "/busco/full_table.tsv",
  params:
    out_path = os.getcwd() + "/busco",
    odb = "vertebrata_odb10",
    buscobase = "assembly",
    rmcmd = "echo 'Removing BUSCO run dir:{params.out_path}{params.buscobase}'; \
            rm -r {params.out_path}{params.buscobase};" if keepfiles == True else "" 
  log:
    "logs/" + str(date) + ".busco.out",
    "logs/" + str(date) + ".busco.err",
  benchmark:
    "{dir}/logs/" + str(date) + ".busco.benchmark.txt",
  conda:
    '../envs/busco5.4.0.yaml'
  threads: 8
  shell:
    "cd {params.out_path};" 
    "busco -i {input.assembly} -c {threads} -m genome --offline -o {params.buscobase} --out_path {params.out_path} -l {input.lineage};"
    "mv {params.out_path}{params.buscobase}/run_{params.odb}/short_summary.txt {output.summary};"
    "mv {params.out_path}{params.buscobase}/run_{params.odb}/full_table.tsv {output.full};"
    "{params.rmcmd}"
  
  #if keepfiles == False:
   #   shell(
    #    "echo 'Removing BUSCO run dir:{params.out_path}{params.buscobase}';"
     #   "rm -r {params.out_path}{params.buscobase};" 
     # )

rule run_merqury:
  input:
    meryl_db = os.getcwd() + "/database.meryl",
    assembly = os.getcwd() + "/assembly.fasta"
  output:
    completeness = os.getcwd() + "/merqury_run/completeness.stats",
    hist = os.getcwd() + "/merqury_run/assembly.assembly.spectra-cn.hist",
    plots = [os.getcwd() + "/merqury_run/assembly.spectra-cn.ln.png", os.getcwd() + "/merqury_run/assembly.spectra-cn.fl.png", os.getcwd() + "/merqury_run/assembly.spectra-cn.st.png"]
  params:
    out_pref = "assembly",
    conda_env = "~fcruz/.conda/envs/merqury_v1.1/",
    directory = os.getcwd() + "/merqury_run",
  log:
    "logs/" + date + ".merqury.out",
    "logs/" + date + ".merqury.err",
  threads: 1
  shell:
    "module purge;"
    "source ~jgomez/init_shell.sh;"
    "conda activate {params.conda_env};"
    "mkdir -p {params.directory}; cd {params.directory};"
    "merqury.sh {input.meryl_db} {input.assembly} {params.out_pref};"
    "{params.conda_env}/share/merqury/plot/plot_spectra_cn.R -f {output.hist} -o {params.out_pref}.spectra-cn -z {params.out_pref}.{params.out_pref}.only.hist;"
    "conda deactivate;"

rule align_ont:
  input:
    genome = os.getcwd() + "/assembly.fasta",
    reads = os.getcwd() + "/reads.ont.fastq.gz"
  output:
    mapping = "minimap2.bam"
  params:
    env = "/scratch/project/devel/aateam/src/RACON/v1.4.20_conda_env/",
    align_opts = "ax map-ont",
    tmp = "minimap2.sam"
  log:
    "logs/" + date + ".minimap2.out",
    "logs/" + date + ".minimap2.err",
  threads: 4
  run:
    shell(
      "module purge; source ~jgomez/init_shell.sh; conda activate {params.env};" 
      "minimap2 -{params.align_opts} -t {threads} {input.genome} {input.reads} > {params.tmp};"
      "conda deactivate;"
    )
    if params.align_opts == "ax map-ont":
      shell(
        "module purge;"
        "module load SAMTOOLS/1.12;"
        "samtools view -Sb {params.tmp}| samtools sort -@ {threads} -o {output.mapping} -;"
        "samtools index {output.mapping};"
        "rm {params.tmp};"
      )
    elif params.align_opts == "x map-ont":
      shell (
        "gzip -c {params.tmp} > {output.mapping};"
        "rm {params.tmp};"
      )

rule align_illumina:
  input:
    genome = os.getcwd() + "/assembly.fasta",
    reads = [os.getcwd() + "/reads.illumina.1.fastq.gz", os.getcwd() + "/reads.illumina.2.fastq.gz"]
 #   pe1 = os.getcwd() + "/reads.illumina.1.fastq.gz",
  #  pe2 = os.getcwd() + "/reads.illumina.2.fastq.gz",
  output:
    mapping = "illumina.bam"
  params:
    options = ""
  log:
    "logs/" + date + ".bwa2.out",
    "logs/" + date + ".bwa2.err",
  threads: 4
  shell:
    "module purge;"
    "module load SAMTOOLS/1.12 bwa java/1.8.0u31 PILON/1.21;"
    "bwa index {input.genome};"
    "bwa mem -Y {params.options} -t {threads} {input.genome} {input.reads} | samtools view -Sb - | samtools sort -@ {threads} -o {output.mapping} -;"
    "samtools index {output.mapping};"
