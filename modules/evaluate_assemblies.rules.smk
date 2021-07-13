shell.prefix("source ~jgomez/init_shell.sh;")

from datetime import datetime
import re
import os

date = datetime.now().strftime('%Y%m%d.%H%M%S')
keepfiles = False
scripts_dir = os.path.dirname(sys.argv[0]) + "/../scripts/"  

rule finalize_assembly:
  input: 
    "assembly.fasta"
  output:
    "assembly.scaffolds.fa"
  params:
    name = "assembly",
    s2c = "",
    outdir = os.getcwd(),
    evaldir =os.getcwd() + "evaluations/",
  threads: 2
  shell:
    "mkdir -p {params.outdir};"
    "cd {params.outdir};"
    "{scripts_dir}scaffolds2contigs.pl -i {input} -name {params.name} {params.s2c};"
    "ln -s {params.name}.scaffolds.fa {output};"
    "mkdir -p {params.evaldir};"
    "cd {params.evaldir};"
    "mkdir -p stats; cd stats;"
    "{scripts_dir}fastalength {final_assembly} | {scripts_dir}Nseries.pl > {params.name}.scaffolds.nseries.txt;"
    "{scripts_dir}fastalength {params.outdir}{params.name}.contigs.fa | {scripts_dir}Nseries.pl > {params.name}.contigs.nseries.txt;"

rule run_busco:
  input:
    assembly = os.getcwd() + "/assembly.fasta",
    lineage = "/scratch/project/devel/aateam/bin/busco_envs/lineages/odb10/vertebrata_odb10"
  output:
    summary = os.getcwd() + "/busco/short_summary.txt",
    full = os.getcwd() + "/busco/full_table.tsv",
  params:
    busco_env = "/scratch/project/devel/aateam/bin/busco_envs/busco_v4.0.6/",
    out_path = os.getcwd() + "/busco",
    odb = "vertebrata_odb10",
    buscobase = "assembly"
  threads: 8
  log:
    "logs/" + date + ".busco.out",
    "logs/" + date + ".busco.err",
  run:
    shell(
      "module purge;"
      "source ~jgomez/init_shell.sh;"
      "conda activate {params.busco_env};"    
      "cd {params.out_path};" 
      "busco -i {input.assembly} -c {threads} -m genome --offline -o {params.buscobase} --out_path {params.out_path} -l {input.lineage};"
      "conda deactivate;"
      "rm -r busco_downloads;"
      "mv {params.out_path}{params.buscobase}/run_{params.odb}/short_summary.txt {output.summary};"
      "mv {params.out_path}{params.buscobase}/run_{params.odb}/full_table.tsv {output.full};"
    )
    if keepfiles == False:
      shell(
        "echo 'Removing BUSCO run dir:{params.out_path}{params.buscobase}';"
        "module load bsc;"
        "lrm -o {params.out_path}{params.buscobase}.rmlist.sh {params.out_path}{params.buscobase}; sh {params.out_path}{params.buscobase}.rmlist.sh;"
      )

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
