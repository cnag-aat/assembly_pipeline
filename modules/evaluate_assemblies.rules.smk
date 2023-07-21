from datetime import datetime
import re
import os

date = datetime.now().strftime('%Y%m%d.%H%M%S')
scripts_dir = ""
keepfiles = False

rule finalize:
  input:
    assembly = "assembly.fasta",
    buscos = "busco_short_summary.txt",
    stats = "assembly.stats.txt",
    merqs= "assembly.qv",
  output:
    output = "stats.txt",
  params:
    scripts_dir = "../scripts/",
    rmcmd = "",
  conda:
    '../envs/ass_base.yaml'
  threads: 1
  shell:
    "{params.scripts_dir}get_final_tbl_mult.py -s {input.stats} -b {input.buscos} -m {input.merqs} > {output.output};"
    "{params.rmcmd}"
    "echo 'Pipeline completed';"

# rule get_report:
#   input:
#     stats = "stats,txt",
#     config = "assembly_pipeline.config"
#   output:
#     report = "{base}.report.zip"
#   conda:
#     '../envs/ass_base.yaml'
#   threads: 1
#   shell:
#     "snakemake --snakefile /software/assembly/pipelines/Assembly_pipeline/v2.0/assembly_pipeline/bin/assembly_pipeline.smk --configfile {input.config} --report {output.report};"

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

rule run_merqury:
  input:
    meryl_db = os.getcwd() + "/database.meryl",
    assembly = os.getcwd() + "/assembly.fasta"
  output:
    completeness = os.getcwd() + "/merqury_run/assembly.completeness.stats",
    hist = os.getcwd() + "/merqury_run/assembly.assembly.spectra-cn.hist",
    false_dups = os.getcwd() + "/merqury_run/assembly.false_duplications.txt",
    plots = [os.getcwd() + "/merqury_run/assembly.spectra-cn.ln.png", os.getcwd() + "/merqury_run/assembly.spectra-cn.fl.png", os.getcwd() + "/merqury_run/assembly.spectra-cn.st.png"]
  params:
    out_pref = "assembly",
    directory = os.getcwd() + "/merqury_run",
  log:
    "logs/" + date + ".merqury.out",
    "logs/" + date + ".merqury.err",
  conda:
    "../envs/merqury1.3.yaml"
  threads: 4
  shell:
    "mkdir -p {params.directory}; cd {params.directory};"
    "merqury.sh {input.meryl_db} {input.assembly} {params.out_pref};"
    "$CONDA_PREFIX/share/merqury/plot/plot_spectra_cn.R -f {output.hist} -o {params.out_pref}.spectra-cn -z {params.out_pref}.{params.out_pref}.only.hist;"
    "$CONDA_PREFIX/share/merqury/eval/false_duplications.sh {output.hist} > {output.false_dups};"

rule run_busco:
  input:
    assembly = os.getcwd() + "/assembly.fasta",
    lineage = "/scratch/groups/assembly/shared/databases/busco_v5/busco_v5_lineages/lineages/vertebrata_odb10"
  output:
    summary = os.getcwd() + "/busco/short_summary.txt",
    full = os.getcwd() + "/busco/full_table.tsv",
  params:
    out_path = os.getcwd() + "/busco",
    odb = "vertebrata_odb10",
    buscobase = "assembly",
    rmcmd = "echo 'Removing BUSCO run dir: busco'; \
            rm -r busco;" if keepfiles == True else "" 
  conda:
    '../envs/busco5.4.0.yaml'
  threads: 8
  shell:
    "cd {params.out_path};" 
    "busco -i {input.assembly} -c {threads} -m genome --offline -o {params.buscobase} --out_path {params.out_path} -l {input.lineage};"
    "mv {params.out_path}{params.buscobase}/run_{params.odb}/short_summary.txt {output.summary};"
    "mv {params.out_path}{params.buscobase}/run_{params.odb}/full_table.tsv {output.full};"
    "{params.rmcmd}"
  
    
rule align_ont:
  input:
    genome = os.getcwd() + "/assembly.fasta",
    reads = os.getcwd() + "/reads.ont.fastq.gz"
  output:
    mapping = "minimap2.bam"
  params:
   # align_opts = "ax map-ont",
    align_opts = "ax",
    type = "map-ont",
    tmp = "minimap2.sam",
    compress_cmd = lambda wildcards : "samtools view -Sb " + wildcards.directory + "/mappings/" + wildcards.name + "_" + wildcards.ext + ".tmp | " \
                   "samtools sort -@ " + str(config["Parameters"]["minimap2_cores"]) +" -o " + wildcards.directory + "/mappings/" + wildcards.name + "_" + wildcards.ext +";" +\
                   "samtools index " + wildcards.directory + "/mappings/" + wildcards.name + "_" + wildcards.ext 
  conda:
    "../envs/minimap2.24.yaml"
  threads: 4
  shell:
    "minimap2 -{params.align_opts} {params.type} -t {threads} {input.genome} {input.reads} > {params.tmp};"
    "{params.compress_cmd};"
    "rm {params.tmp};"

rule align_illumina:
  input:
    genome = os.getcwd() + "/assembly.fasta",
    reads = [os.getcwd() + "/reads.illumina.1.fastq.gz", os.getcwd() + "/reads.illumina.2.fastq.gz"]
  output:
    mapping = "illumina.bam",
    stats = "illumina.stats.txt"
  params:
    options = ""
  conda:
    "../envs/bwa-mem2.2.1.yaml"
  threads: 4
  shell:
    "echo 'this should be in the .out file';"
    ">&2 echo 'this should be in the .err file';"
    "env;"
    "bwa-mem2 index {input.genome};"
    "bwa-mem2 mem -Y {params.options} -t {threads} {input.genome} {input.reads} | samtools view -Sb - | samtools sort -@ {threads} -o {output.mapping} -;"
    "samtools index {output.mapping};"
    "samtools stats -@ {threads} {output.mapping} > {output.stats};"