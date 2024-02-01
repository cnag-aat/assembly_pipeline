
from datetime import datetime
import re
import os

rule assembly_prepare:
  input: 
   "assembly_in.fa"
  output:
    glen = "assembly.genome",
    bwa = "assembly.bwt",
    faidx = "assembly.fai"   
  params:
    scripts_dir = "../scripts/",  
    workdir = "assembly/",
  conda:
    "../envs/dovetail_tools.yaml"
  shell:
    """
    export PATH="{params.scripts}:$PATH";
    cd {params.workdir}/; 
    fastalength {input} | gawk '{{print $2"\t"$1}}' > {output.glen};
    bwa index {input}; 
    samtools faidx {input};
    """

rule run_yahs:
  input: 
    mappedptsort = "mapped.PT.mq.name_sorted.bam",
    sla = "assembly.fa",
    index = "assembly.fa.fai"
  output:
    outyahs = "yahs_scaffolds_final.fa",
    agp = "yahs_scaffolds_final.agp"
  params:
    yahsdir = "s06.1_p05.1_HiC_scaffolding/",
    yahsopts = "",
    mq = 40,
    name = "assembly"
  conda:
    "../envs/yahs1.2a.2.yaml"
  shell:
    """
    mkdir -p {params.yahsdir} 
    cd {params.yahsdir}
    yahs {input.sla} {params.yahsopts} {input.mappedptsort} 
    ln -s {params.yahsdir}/yahs.out_scaffolds_final.fa {output.outyahs}
    ln -s {params.yahsdir}/yahs.out_scaffolds_final.agp {output.agp}
    """

rule generate_pretext:
  input:
    mapbam = "mapped.PT.mq40.bam"  
  output:
    pret = "assembly_mq40.pretext", 
    ptd = "snapshots/PTdonemq40.txt"
  params:
    scripts_dir = "../scripts",
    outd = "s06.1_p05.1_HiC_scaffolding",
    name = "assembly",
    mq = 40,
  conda:
    "../envs/pretext-suite0.0.2_plus_samtools1.6.yaml"
  shell:
    """
    export PATH="{params.scripts_dir}:$PATH;"
    mkdir -p {params.outd}/snapshots/three_wave_blue_green_yellow/
    cd {params.outd}/

    samtools view -@ {threads} -h {input.mapbam} | PretextMap -o {output.pret} \
    --sortby length --sortorder descend --mapq {params.mq}

    PretextSnapshot -m {output.pret}  --sequences "=full" \
    -o snapshots/three_wave_blue_green_yellow

    PretextSnapshot -m {output.pret} --sequences "=all" \
    -o snapshots/three_wave_blue_green_yellow

    touch {output.ptd}    
    """

rule add_extensions_pretext:
  input:
    tel = "telomeres.bg",
    gaps = "gaps.bg",
    pret = "assembly_mq.pretext",
    ontcov = "ONTcoverage.bg",
  output:
    pretext = "assembly_mq.extensions.pretext"
  params:
    scripts_dir = "../scripts",
    outd = "s06.1_p05.1_HiC_scaffolding",
    mq = 40,
  conda:
    "../envs/pretext-suite0.0.2_plus_samtools1.6.yaml"
  shell:
    """
    export PATH="{params.scripts_dir}:$PATH;"
    cd {params.outd}
    cp {input.pret} {output.pretext}

    if [[ -s "{input.gaps}" ]]; then
    cat {input.gaps} | PretextGraph -i {output.pretext} -n "GAPs"
    fi

    if [[ -s "{input.tel}" ]]; then 
    cat {input.tel} | PretextGraph -i {output.pretext} -n "TELOMERES"
    fi

    if [[ -s "{input.ontcov}" ]]; then
    cat {input.ontcov} | PretextGraph -i {output.pretext} -n "ONTcov" 
    fi
    """

rule get_tpf:
  input:
    fasta = "yahs.out_scaffolds_final.fa"
  output:
    tpf = "genome.yahs_scaffolded.fa.tpf"
  params:
    scripts_dir = "../scripts/",
    dir = "s06.1_p05.1_HiC_scaffolding"
  conda:
    "../envs/rapidcuration.yaml"
  shell:
    "mkdir -p {params.dir}/run_yahs; cd {params.dir}/run_yahs;"
    "{params.scripts_dir}rapid_split.pl -fa yahs.out_scaffolds_final.fa;"
    "ln -s run_yahs/yahs.out_scaffolds_final.fa.tpf {output.tpf};"



