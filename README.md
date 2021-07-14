# assembly_pipeline
 Snakemake Pipeline used for de novo genome assembly @CNAG. It has been developed using Snakemake v6.0.5, which can be loaded in the cluster by doing:
 ```bash
 conda activate /home/devel/jgomez/conda_environments/snakemake
 ```

Currently it can take ONT reads, illumina paired-end data and illumina 10X data. It does the preprocessing of the reads, assembly, polishing and evaluations. By default it will preprocess the reads, run Flye + Hypo and evaluate the resulting assemblies with BUSCO, MERQURY and Nseries. 
It needs a config file and a spec file (json file with instructions on which resources to use in the CNAG cluster for each of the jobs). Both files are created by the script "create_config_assembly.py" that is located in the bin directory. To check all the options accepted by the script, do:

```
bin/create_config_assembly.py -h
```

# How to provide input data:

There are several ways of providing the reads for assembly and polishing.

1- ONT reads

1.1 Using the option ``--ont-dir {DIR}`` in create_config_assembly.py:

If you do so, it will look for all the files in the directory that end in '.fastq.gz' and will add the basenames to "ONT_wildcards". These wildcards will be processed by the pipeline that will: 

- Concatenate all the files into a single file

- Run filtlong with the specified parameters for minlen and min_mean_q

- Use the resulting file for assembly and/or polishing

You can also specify the basenames of the files that you want to use with the ``--ont-list `` option. In this case, the pipeline will use the wildcards that you're providing instead of merging all the files in the directory.

1.2 Using the option ```--ont-reads {FILE}``` in create_config_assembly.py:

If you do so, it will consider that you already have all the reads in one file and will:  

- Run filtlong with the specified parameters for minlen and min_mean_q

- Use the resulting file for assembly and/or polishing

1.3 Using the option ```--ont-filt {FILE}```. It will use this file as the output from filtlong and will skip the preprocessing steps and directly use it for assembly and/or polishing. 



2-Illumina 10X-linked data

2.1 Using the  ```--raw-10X {DIR}``` and ``` --10X-list``` options:

It will take each basename in the list to get the corresponding fastqs from the directory and run longranger on them. Afterwards, it will build meryldbs for each "barcoded" file. Finally, it will concatenate all the meryldbs and "barcoded" files. Resulting "barcoded" file will be used for polishing. 

2.2 Using the ``--processed-10X {DIR}`` parameter. This directory can already be there or be produced by the pipeline as described in step 2.1. Once all the "barcoded" fastq files are there, meryldbs will be built for each "barcoded" file.  Finally, it will concatenate all the meryldbs and "barcoded" files. Resulting "barcoded" file will be used for polishing. 

2.3 Using the ``--10X`` option. The argument would be the path to the concatenated ".barcoded" file that needs to be used for polishing. If the pre-concatenated files are not given, meryldbs will be directly generated with this file, but it may run out of memory. 

3- Illumina short-read data

3.1 Using the ``--illumina-dir {DIR}`` option, that will look for all the files in the directory that end in '.1.fastq.gz' and will add the basenames to "illumina_wildcards". These wildcards will be processed by the pipeline that will: 

- Concatenate all the *.1.fastq.gz and the *2.fastq.gz in one file per pair. 

- The resulting reads will be used for building meryldbs and polishing. 

3.2 Using the ``--pe1 {FILE} and --pe2 {FILE}`` options. That will consider that these are the paired files containing all the illumina reads to be used and will build meryldbs and polish with them.

4- Input assemblies

If you want to polish an already assembled assembly, you can give it to the pipeline by using the option ``--assembly-in ASSEMBLY_IN [ASSEMBLY_IN ...]
                        Dictionary with assemblies that need to be polished but not assembled and directory where they should
                        be polished. Example: '{"assembly1":"polishing_dir1"}' '{"assembly2"="polishing_dir2"}' ...``
			


# Description of implemented rules

1- Preprocessing:
	
- Read concatenation:

``zcat {input.fastqs} | pigz -p {threads} -c  > {output.final_fastq}``
	
- Longranger for 10X reads: it uses the Longranger version installed in the path specified in the configfile

``longranger basic --id={params.sample} --sample={params.sample} --fastqs={input.mkfastq_dir} --localcores={threads}``

- Filtlong:  it uses the Filtlong version installed in the path specified in the configfile

``filtlong --min_length {params.minlen} --min_mean_q {params.min_mean_q} {params.opts} {input.reads} | pigz -p {threads} -c > {output.outreads}``
	
- Build meryldb (with processed 10X reads or illumina reads): with the merqury conda environment specified in the configfile

``meryl k={params.kmer} count output {output.out_dir} {input.fastq}``
	
- Concat meryldbs: with the merqury conda environment specified in the configfile

``meryl union-sum output {output.meryl_all} {input.input_run}``
	
- Align ONT (Minimap2): it aligns the reads using minimap2 and outputs the alignment either in bam or in paf.gz formats. It uses the minimap2 conda environment specified in the configfile

``minimap2 -{params.align_opts} -t {threads} {input.genome} {input.reads} ``

- Align Illumina (BWA-MEM): it aligns the reads with BWA-mem and outputs a bam file

``bwa mem -Y {params.options} -t {threads} {input.genome} {input.reads} | samtools view -Sb - | samtools sort -@ {threads} -o {output.mapping} -``

2- Assembly

- Flye (default)
	
- Nextdenovo (if turned on)

3- Polishing

- Hypo (default)
	
- Racon (if turned on)
	
- Medaka (if turned on)
	
- Pilon (if turned on)
	
- Nextpolish ont (if turned on)
	
- Nextpolish illumina (if turned on)

4- Evaluations
	
- Merqury
	
- Busco
	
- Nseries

# Description of all options

```
bin/create_config_assembly.py -h
usage: create_configuration_file [-h] [--configFile configFile] [--specFile specFile] [--ndconfFile ndconfFile]
                                 [--logs-dir logs_dir] [--concat-cores concat_cores] [--genome-size genome_size]
                                 [--lr-type lr_type] [--basename base_name] [--keep-intermediate]
                                 [--preprocess-lr-step PREPROCESS_ONT_STEP] [--preprocess-10X-step PREPROCESS_10X_STEP]
                                 [--flye-step FLYE_STEP] [--no-flye] [--nextdenovo-step NEXTDENOVO_STEP] [--run-nextdenovo]
                                 [--racon-cores racon_cores] [--nextpolish-cores nextpolish_cores]
                                 [--minimap2-cores minimap2_cores] [--bwa-cores bwa_cores] [--pilon-cores pilon_cores]
                                 [--medaka-cores medaka_cores] [--hypo-cores hypo_cores] [--busco-cores busco_cores]
                                 [--racon-rounds racon_rounds] [--pilon-rounds pilon_rounds] [--medaka-rounds medaka_rounds]
                                 [--nextpolish-ont-rounds nextpolish_ont_rounds]
                                 [--nextpolish-ill-rounds nextpolish_ill_rounds] [--hypo-rounds hypo_rounds]
                                 [--longranger-cores longranger_cores] [--longranger-path longranger_path]
                                 [--scripts-dir SCRIPTS_DIR] [--ont-reads ONT_READS] [--ont-dir ONT_DIR]
                                 [--ont-filt ONT_FILTERED] [--pe1 PE1] [--pe2 PE2] [--raw-10X RAW_10X]
                                 [--processed-10X PROCESSED_10X] [--10X R10X] [--illumina-dir ILLUMINA_DIR]
                                 [--assembly-in ASSEMBLY_IN [ASSEMBLY_IN ...]] [--pipeline-workdir PIPELINE_WORKDIR]
                                 [--filtlong-dir FILTLONG_DIR] [--flye-dir FLYE_DIR] [--nextdenovo-dir NEXTDENOVO_DIR]
                                 [--flye-polishing-dir POLISH_FLYE_DIR] [--nextdenovo-polishing-dir POLISH_NEXTDENOVO_DIR]
                                 [--eval-dir eval_dir] [--stats-out stats_out] [--filtlong-path filtlong_path]
                                 [--filtlong-minlen filtlong_minlen] [--filtlong-min-mean-q filtlong_min_mean_q]
                                 [--filtlong-opts filtlong_opts] [--flye-env flye_env] [--flye-cores flye_cores]
                                 [--flye-polishing-iterations flye_pol_it] [--other-flye-opts other_flye_opts]
                                 [--nextdenovo-module nextdenovo_module] [--nextdenovo-cores nextdenovo_cores]
                                 [--nextdenovo-task nextdenovo_task] [--nextdenovo-rewrite nextdenovo_rewrite]
                                 [--nextdenovo-parallel_jobs nextdenovo_parallel_jobs]
                                 [--nextdenovo-minreadlen nextdenovo_minreadlen]
                                 [--nextdenovo-seeddepth nextdenovo_seeddepth]
                                 [--nextdenovo-seedcutoff nextdenovo_seedcutoff]
                                 [--nextdenovo-blocksize nextdenovo_blocksize]
                                 [--nextdenovo-pa-correction  nextdenovo_pa_correction]
                                 [--nextdenovo-minimap_raw nextdenovo_minimap_raw]
                                 [--nextdenovo-minimap_cns nextdenovo_minimap_cns]
                                 [--nextdenovo-minimap_map nextdenovo_minimap_map] [--nextdenovo-sort nextdenovo_sort]
                                 [--nextdenovo-correction_opts nextdenovo_correction_opts]
                                 [--nextdenovo-nextgraph_opt nextdenovo_nextgraph_opt] [--hypo-env hypo_env]
                                 [--sr-cov ill_cov] [--hypo-proc hypo_processes] [--hypo-opts hypo_opts]
                                 [--minimap-env minimap_env] [--racon-dir racon_dir] [--racon-opts racon_opts]
                                 [--medaka-env medaka_env] [--medaka-workdir medaka_workdir] [--medaka-model medaka_model]
                                 [--medaka-consensus-opts medaka_consensus_opts] [--pilon-path pilon_path]
                                 [--pilon-opts pilon_opts] [--java-opts java_opts] [--pilon-subs pilon_subsampling]
                                 [--pilon-chunks pilon_chunks] [--intermediate-evals] [--no-final-evals]
                                 [--busco-env busco_env] [--busco-lin busco_lineage] [--merqury-env merqury_env]
                                 [--merqury-db merqury_db] [--meryl-k meryl_k] [--ont-list ONT_wildcards]
                                 [--illumina-list illumina_wildcards] [--r10X-list r10X_wildcards]

Create a configuration json file for the assembly pipeline.

optional arguments:
  -h, --help            show this help message and exit

General Parameters:
  --configFile configFile
                        Configuration JSON to be generated. Default assembly.config
  --specFile specFile   Cluster specifications JSON fileto be generated. Default assembly.spec
  --ndconfFile ndconfFile
                        Name pf the nextdenovo config file. Default nextdenovo.config
  --logs-dir logs_dir   Directory to keep all the log files. Default logs/
  --concat-cores concat_cores
                        Number of threads to concatenate reads and to run filtlong. Default 4
  --genome-size genome_size
                        Approximate genome size. Example: 615m or 2.6g. Default None
  --lr-type lr_type     Type of long reads (options are flye read-type options). Default nano-raw
  --basename base_name  Base name for the project. Default None
  --keep-intermediate   Set this to True if you do not want intermediate files to be removed. Default False
  --preprocess-lr-step PREPROCESS_ONT_STEP
                        Step for preprocessing long-reads. Default 01.1
  --preprocess-10X-step PREPROCESS_10X_STEP
                        Step for preprocessing 10X reads. Default 01.2
  --flye-step FLYE_STEP
                        Step for running flye. Default 02.1
  --no-flye             Give this option if you do not want to run Flye.
  --nextdenovo-step NEXTDENOVO_STEP
                        Step for running nextdenovo. Default 02.2
  --run-nextdenovo      Give this option if you do want to run Nextdenovo.
  --racon-cores racon_cores
                        Number of threads to run the racon step. Default 16
  --nextpolish-cores nextpolish_cores
                        Number of threads to run the nextpolish step. Default 14
  --minimap2-cores minimap2_cores
                        Number of threads to run the alignment with minimap2. Default 16
  --bwa-cores bwa_cores
                        Number of threads to run the alignments with BWA-Mem. Default 16
  --pilon-cores pilon_cores
                        Number of threads to run the pilon step. Default 16
  --medaka-cores medaka_cores
                        Number of threads to run the medaka step. Default 16
  --hypo-cores hypo_cores
                        Number of threads to run the hypo step. Default 24
  --busco-cores busco_cores
                        Number of threads to run BUSCO. Default 16
  --racon-rounds racon_rounds
                        Number of rounds of racon to run. Default 0
  --pilon-rounds pilon_rounds
                        Number of rounds of pilon to run. Default 0
  --medaka-rounds medaka_rounds
                        Number of rounds of medaka to run. Default 0
  --nextpolish-ont-rounds nextpolish_ont_rounds
                        Number of rounds to run the Nextpolish with ONT step. Default 0
  --nextpolish-ill-rounds nextpolish_ill_rounds
                        Number of rounds to run the Nextpolish with illumina step. Default 0
  --hypo-rounds hypo_rounds
                        Number of rounds to run the Hypostep. Default 1
  --longranger-cores longranger_cores
                        Number of threads to run longranger. Default 8
  --longranger-path longranger_path
                        Path to longranger executable. Default /scratch/project/devel/aateam/src/10X/longranger-2.2.2

Inputs:
  --scripts-dir SCRIPTS_DIR
                        Directory with the different scripts for the pipeline. Default
                        /scratch/project/devel/aateam/src/assembly_pipeline/bin/../scripts/
  --ont-reads ONT_READS
                        File with all the ONT reads. Default None
  --ont-dir ONT_DIR     Directory where the ONT fastqs are stored. Default None
  --ont-filt ONT_FILTERED
                        File with the ONT reads after running filtlong on them. Default None
  --pe1 PE1             File with the illumina paired-end fastqs, pair 1.
  --pe2 PE2             File with the illumina paired-end fastqs, pair 2.
  --raw-10X RAW_10X     Directory to mkfastq Raw 10X reads.
  --processed-10X PROCESSED_10X
                        Directory to Processed 10X reads. Already there or to be produced by the pipeline.
  --10X R10X            File with barcoded 10X reads in fastq.gz format, concatenated.
  --illumina-dir ILLUMINA_DIR
                        Directory where the illumina fastqs are stored. Default None
  --assembly-in ASSEMBLY_IN [ASSEMBLY_IN ...]
                        Dictionary with assemblies that need to be polished but not assembled and directory where they should
                        be polished. Example: '{"assembly1":"polishing_dir1"}' '{"assembly2"="polishing_dir2"}' ...

Outputs:
  --pipeline-workdir PIPELINE_WORKDIR
                        Base directory for the pipeline run. Default /scratch/devel/jgomez/test_assembly_pipeline/
  --filtlong-dir FILTLONG_DIR
                        Directory to process the ONT reads with filtlong. Default s01.1_p01.1_Filtlong
  --flye-dir FLYE_DIR   Directory to run flye. Default s02.1_p01.1_flye/
  --nextdenovo-dir NEXTDENOVO_DIR
                        Directory to run nextdenovo. Default s02.2_p01.1_nextdenovo/
  --flye-polishing-dir POLISH_FLYE_DIR
                        Directory to polish the flye assembly. Default s03.1_p02.1_polishing/
  --nextdenovo-polishing-dir POLISH_NEXTDENOVO_DIR
                        Directory to run nextdenovo. Default s03.2_p02.2_polishing/
  --eval-dir eval_dir   Base directory for the evaluations. Default evaluations/
  --stats-out stats_out
                        Path to the file with the final statistics.

Filtlong:
  --filtlong-path filtlong_path
                        Path to the filtlong software. Default /scratch/project/devel/aateam/bin/filtlong
  --filtlong-minlen filtlong_minlen
                        Minimum read length to use with Filtlong. Default 1000
  --filtlong-min-mean-q filtlong_min_mean_q
                        Minimum mean quality to use with Filtlong. Default 80
  --filtlong-opts filtlong_opts
                        Extra options to run Filtlong (eg. -t 4000000000)

Flye:
  --flye-env flye_env   Conda environment to run FLYE. Default /home/devel/talioto/miniconda3/envs/flye-v2.8.3/
  --flye-cores flye_cores
                        Number of threads to run FLYE. Default 24
  --flye-polishing-iterations flye_pol_it
                        Number of polishing iterations to use with FLYE. Default 2
  --other-flye-opts other_flye_opts
                        Additional options to run Flye. Default None

Nextdenovo:
  --nextdenovo-module nextdenovo_module
                        Cluster module to run nextdenovo. Default NEXTDENOVO/2.4.0
  --nextdenovo-cores nextdenovo_cores
                        Number of threads to run nextdenovo. Default 24
  --nextdenovo-task nextdenovo_task
                        Task need to run. Default all
  --nextdenovo-rewrite nextdenovo_rewrite
                        Overwrite existing directory. Default yes
  --nextdenovo-parallel_jobs nextdenovo_parallel_jobs
                        Number of tasks used to run in parallel. Default 4
  --nextdenovo-minreadlen nextdenovo_minreadlen
                        Filter reads with length < minreadlen. Default 1k
  --nextdenovo-seeddepth nextdenovo_seeddepth
                        Expected seed depth, used to calculate seed_cutoff, co-use with genome_size, you can try to set it
                        30-45 to get a better assembly result. Default 45
  --nextdenovo-seedcutoff nextdenovo_seedcutoff
                        Minimum seed length, <=0 means calculate it automatically using bin/seq_stat. Default 0
  --nextdenovo-blocksize nextdenovo_blocksize
                        Block size for parallel running, split non-seed reads into small files, the maximum size of each file
                        is blocksize. Default 1g
  --nextdenovo-pa-correction  nextdenovo_pa_correction
                        number of corrected tasks used to run in parallel, each corrected task requires ~TOTAL_INPUT_BASES/4
                        bytes of memory usage, overwrite parallel_jobs only for this step. Default 4
  --nextdenovo-minimap_raw nextdenovo_minimap_raw
                        minimap2 options, used to find overlaps between raw reads, see minimap2-nd for details. Default -t 6
                        -x ava-ont
  --nextdenovo-minimap_cns nextdenovo_minimap_cns
                        minimap2 options, used to find overlaps between corrected reads. Default -t 6 -x ava-ont -k17 -w17
  --nextdenovo-minimap_map nextdenovo_minimap_map
                        minimap2 options, used to map reads back to the assembly. Default -t 6 -x ava-ont
  --nextdenovo-sort nextdenovo_sort
                        sort options, see ovl_sort for details. Default -m 40g -t 20
  --nextdenovo-correction_opts nextdenovo_correction_opts
                        Correction options. Default -p 6
  --nextdenovo-nextgraph_opt nextdenovo_nextgraph_opt
                        nextgraph options, see nextgraph for details. Default -a 1

Hypo:
  --hypo-env hypo_env   Conda environment to run Hypo. Default /scratch/project/devel/aateam/src/HyPo/HyPov1_conda_env/
  --sr-cov ill_cov      Approximate short read coverage for hypo Default 0
  --hypo-proc hypo_processes
                        Number of contigs to be processed in parallel by HyPo. Default 6
  --hypo-opts hypo_opts
                        Additional options to run Hypo. Default None

Racon:
  --minimap-env minimap_env
                        Conda environment to run Minimap2. Default /scratch/project/devel/aateam/src/RACON/v1.4.20_conda_env
  --racon-dir racon_dir
                        Directory with Racon installation. Default /scratch/project/devel/aateam/src/RACON/v1.4.21_github/
  --racon-opts racon_opts
                        Extra options to run Racon_wrapper(eg. --split 100000000) to split the assembly in chunks of
                        specified size and decrease memory requirements. Do racon_wrapper -h for more info.

Medaka:
  --medaka-env medaka_env
                        Conda environment to run Medaka. Default /scratch/project/devel/aateam/src/MEDAKA/medaka-141
  --medaka-workdir medaka_workdir
                        Directory to run Medaka. Default $TMPDIR/
  --medaka-model medaka_model
                        User needs to specify the model that will be used to run Medaka
  --medaka-consensus-opts medaka_consensus_opts
                        Specify any parameters to change when running medaka consensus

Pilon:
  --pilon-path pilon_path
                        Path to Pilon executable. Default /apps/PILON/1.21/pilon
  --pilon-opts pilon_opts
                        Additional options to run Pilon. Default --fix bases --changes
  --java-opts java_opts
                        Options for the java execution of Pilon. Default None
  --pilon-subs pilon_subsampling
                        Percentage of reads to randomly use when running pilon. Default 0
  --pilon-chunks pilon_chunks
                        Number of chunks to split the pilon polishing jobs. Default 25

Finalize:
  --intermediate-evals  If specified, run evaluations on intermediate assemblies. Default False
  --no-final-evals      If specified, do not run evaluations on final assemblies. Default True
  --busco-env busco_env
                        Conda environment to run BUSCO. Default /scratch/project/devel/aateam/bin/busco_envs/busco_v4.0.6/
  --busco-lin busco_lineage
                        Path to the lineage directory to run Busco with. Default None
  --merqury-env merqury_env
                        Conda environment to run merqury. Default /home/devel/fcruz/.conda/envs/merqury_v1.1/
  --merqury-db merqury_db
                        Meryl database. Default None
  --meryl-k meryl_k     Kmer length to build the meryl database. Default None

Wildcards:
  --ont-list ONT_wildcards
                        List with basename of the ONT fastqs that will be used. Default None
  --illumina-list illumina_wildcards
                        List with basename of the illumina fastqs. Default None
  --r10X-list r10X_wildcards
                        List with basename of the raw 10X fastqs. For raw 10X we need to give this argument, for processed
                        10X reads, the pipeline can obtain it. Default None


```

# Steps that will be included soon
	
Trim_galore (rule done, need to integrate it)


10X scaffolding (pending)
	
Shasta?

