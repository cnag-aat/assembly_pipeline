# assembly_pipeline
 Snakemake Pipeline used for de novo genome assembly @CNAG. It has been developed for Snakemake v6.0.5.

It accepts Oxford Nanopore Technologies (ONT) reads, illumina paired-end data and illumina 10X data. It does the preprocessing of the reads, assembly, polishing, purge_dups, 10X scaffodling and different evaluation steps. By default it will preprocess the reads, run Flye + Hypo + purge_dups and evaluate the resulting assemblies with BUSCO, MERQURY, Nseries and assembly_stats. 
It needs a config file and a spec file (json file with instructions on which resources should slurm use for each of the jobs). Both files are created by the script "create_config_assembly.py" that is located in the bin directory. To check all the options accepted by the script, do:

```
bin/create_config_assembly.py -h
```

Once the 2 config files are produced, the pipeline can be launched using snakemake like this:

``snakemake --notemp -j 999 --snakefile assembly_pipeline.smk --configfile assembly.config --is --cluster-conf assembly.spec  -np ``

If you are using an HPC cluster, please check how should you run snakemake to launch the jobs to the cluster. 

# How to provide input data:

There are several ways of providing the reads.

### 1- ONT reads

1.1 Using the option ``--ont-dir {DIR}`` in create_config_assembly.py.

If you do so, it will look for all the files in the directory that end in '.fastq.gz' and will add the basenames to "ONT_wildcards". These wildcards will be processed by the pipeline that will:

- Concatenate all the files into a single file

- Run filtlong with the default or specified parameters. 

- Use the resulting file for assembly, polishing and/or purging.

You can also specify the basenames of the files that you want to use with the ``--ont-list `` option. In this case, the pipeline will use the wildcards that you're providing instead of merging all the files in the directory.

1.2 Using the option ```--ont-reads {FILE}``` in create_config_assembly.py.

If you do so, it will consider that you already have all the reads in one file and will:  

- Run filtlong with the default or specified parameters.

- Use the resulting file for assembly, polishing and/or purging.

1.3 Using the option ```--ont-filt {FILE}```. It will use this file as the output from filtlong. Hence, it will skip the preprocessing steps and directly use it for assembly, polishing and/or purging. 



### 2-Illumina 10X-linked data

2.1 Using the  ```--raw-10X {DIR:list}``` option. 

Dictionary with 10X raw read directories, it has to be the mkfastq dir. You must specify as well the sampleIDs from this run. Example: '{"mkfastq-                        dir":"sample1,sample2,sample3"}'...

It will take each basename in the list to get the fastqs from the corresponding directory and run longranger on each sample. Afterwards, it will build meryldbs for each "barcoded" file. Finally, it will concatenate all the meryldbs and "barcoded" files. Resulting "barcoded" file will be used for polishing. 

2.2 Using the ``--processed-10X {DIR}`` parameter. 

This directory can already be there or be produced by the pipeline as described in step 2.1. Once all the "barcoded" fastq files are there, meryldbs will be built for each "barcoded" file.  Finally, it will concatenate all the meryldbs and "barcoded" files. Resulting "barcoded" file will be used for polishing. 

2.3 Using the ``--10X`` option. 

The argument to this is the path to the concatenated ".barcoded" file that needs to be used for polishing. If the pre-concatenated files are not given, meryldbs will be directly generated with this file, but it may run out of memory. 

### 3- Illumina short-read data

3.1 Using the ``--illumina-dir {DIR}`` option, that will look for all the files in the directory that end in '.1.fastq.gz' and will add the basenames to "illumina_wildcards". These wildcards will be processed by the pipeline that will: 

- Trim adaptors with Trimgalore

- Concatenate all the trimmed *.1.fastq.gz and the *2.fastq.gz in one file per pair. 

- The resulting reads will be used for building meryldbs and polishing. 

3.2 Using the ``--processed-illumina`` option. If the directory exists and contains files, the pipeline will look for all the files in the directory that end in '.1.fastq.gz' and will add the basenames to "illumina_wildcards". These wildcards will be processed by the pipeline that will:

- Concatenate all the trimmed *.1.fastq.gz and the *2.fastq.gz in one file per pair. 

- The resulting reads will be used for building meryldbs and polishing. 

3.3 Using the ``--pe1 {FILE} and --pe2 {FILE}`` options. That will consider that these are the paired files containing all the illumina reads ready to be used and will build meryldbs and polish with them.

### 4- Input assemblies

If you want to polish an already assembled assembly, you can give it to the pipeline by using the option ``--assembly-in ASSEMBLY_IN [ASSEMBLY_IN ...]
                        Dictionary with assemblies that need to be polished but not assembled and directory where they should
                        be polished. Example: '{"assembly1":"polishing_dir1"}' '{"assembly2"="polishing_dir2"}' ...``
			
If you want to improve an already polished assembly, you can give it to the pipeline by using the option ``--postpolish-assemblies POSTPOLISH_ASSEMBLIES [POSTPOLISH_ASSEMBLIES ...]
                        Dictionary with assemblies for whic postpolishing steps need to be run but that are not assembled and
                        base step for the directory where the first postpolishing step should be run. Example:
                        '{"assembly1":"s04.1_p03.1"}' '{"assembly2"="s04.2_p03.2"}' ...``



```
# Description of implemented rules

1- Preprocessing:
	
- **Read concatenation:**

``zcat {input.fastqs} | pigz -p {threads} -c  > {output.final_fastq}``
	
- **Longranger for 10X reads**: it uses the Longranger version installed in the path specified in the configfile

``longranger basic --id={params.sample} --sample={params.sample} --fastqs={input.mkfastq_dir} --localcores={threads}``

- **Trimgalore:** By default it gives the ``--gzip -q 20 --paired --retain_unpaired`` options, but it can be changed with the ``--trim-galore-opts `` argument. 

``trim_galore -j {threads} {params.opts} {input.read1} {input.read2}``

- **Filtlong:** it uses the Filtlong version installed in the path specified in the configfile. By default it gives the min_length and min_mean_q parameters, but extra parameters can be added with the ``--filtlong-opts`` option.

``filtlong --min_length {params.minlen} --min_mean_q {params.min_mean_q} {params.opts} {input.reads} | pigz -p {threads} -c > {output.outreads}``
	
- **Build meryldb** (with processed 10X reads or illumina reads): it uses the merqury conda environment specified in the configfile. It takes as argument the `--mery-k` value that needs to be estimated first for the genime size. 

``meryl k={params.kmer} count output {output.out_dir} {input.fastq}``
	
- Concat meryldbs: with the merqury conda environment specified in the configfile

``meryl union-sum output {output.meryl_all} {input.input_run}``
	
- **Align ONT (Minimap2):** it aligns the reads using minimap2 and outputs the alignment either in bam or in paf.gz formats. It uses the minimap2 conda environment specified in the configfile

``minimap2 -{params.align_opts} -t {threads} {input.genome} {input.reads} ``

- **Align Illumina (BWA-MEM):** it aligns the reads with BWA-mem and outputs a bam file

``bwa mem -Y {params.options} -t {threads} {input.genome} {input.reads} | samtools view -Sb - | samtools sort -@ {threads} -o {output.mapping} -``

2- Assembly

- **Flye (default)**. It is run by default, if you don't want the pipeline to run it, you can give `--no-flye` option when creating the config. It uses the conda environment specified in the config. By default it is set to 2 polishing iterations and gives the genome-size estimate that has been given when creating the config. Extra options can be provided with the `--flye-opts`.

``flye --{params.readtype} {input.reads} -o {params.outdir}out -t {threads} -i {params.pol_iterations} {params.other_flye_opts} ``
	
- **Nextdenovo (if ``run-nextdenovo``):** It uses the cluster module specified in the config. If nextdenovo option is turned on, the create_config script will also create the nextdenovo config file. Check the create_config help to see which options can be modified on it. 

``nextDenovo {input.config}``

3- Polishing

- **Hypo (default):** It is the polisher that the pipeline uses by default, it can be turned off specifying ``--no-hypo`` when creating the config. If selected, the reads will be aligned in previous rules and then hypo will be run, it requires illumina data. It uses the conda environment specified in the config. 

``hypo -r @short_reads.list.txt -d {input.genome} -b {input.sr_bam} -c {coverage} -s {params.genome_size} -B {input.lr_bam} -t {threads} -o {output.polished} -p {params.proc} {params.opts} ``
	
- **Racon (if turned on):** to run racon, specify ``--racon-rounds `` and the number of rounds of it you want to run. It uses the conda environment specified in the config file. 

``{params.racon_env}/scripts/racon_wrapper.py -u {params.opts} -t {threads} reads4racon.fastq.gz {input.mapping} {input.assembly} > {output.polished} ``
	
- **Medaka (if turned on):** to run medaka, specify ``--medaka-rounds`` and the nummber of rounds of it you want to run. It uses the conda environment specified in the config file. It'll run after racon and before pilon, if they are also selected. 

`` medaka consensus {input.mapping} {wildcards.directory}/rmp/{wildcards.base}.medaka{wildcards.param}.hdf --threads {medaka_threads} --model {params.model} {params.consensus_opts};
medaka stitch --threads {threads} {wildcards.directory}/rmp/{wildcards.base}.medaka{wildcards.param}.hdf {input.assembly} {output.polished}``
	
- **Pilon (if turned on):** to run Pilon, specify ``--pilon-rounds`` and the number of rounds of it you want to run. If it's a big genome, the pipeline will split the consensus step in several jobs, each of them running on certain scaffolds. It uses the version installed in the path specified in the config. 

``{scripts_dir}split_bam.py assembly.len {input.mapping} {params.chunks} {threads};
java {params.java_opts} -jar {params.path} --genome {input.assembly} --frags {input.alignment} {params.opts} --threads {threads} --output {basename}; 
{scripts_dir}/concat_pilon.py {params.splitdir} {params.chunks} > {output.polished}``

	
- **Nextpolish ont (if turned on):** to run nextpolish with ONT reads, specify ``--nextpolish-ont-rounds`` and the number of rounds you want to run of it. 

``"python /apps/NEXTPOLISH/1.3.1/lib/nextpolish2.py -g {input.genome} -p {threads} -l lgs.fofn -r {params.lrtype} > {output.polished}``
	
- **Nextpolish illumina (if turned on):** to run nextpolish with ONT reads, specify ``--nextpolish-ill-rounds`` and the number of rounds you want to run of it. 

``"python /apps/NEXTPOLISH/1.3.1/lib/nextpolish1.py -g {input.genome}  -p {threads} -s {input.bam} -t {params.task} > {output.polished}``

4- Post-assembly

- **Purge_dups (by default):** select ``--no-purgedups`` if you don't want to run it. If no manual cutoffs are given, it'll run purgedups with automatic cutoffs and then will rerun it selecting the mean cutoff as 0.75\*cov. It uses the version installed in the cluster module specified in the config. 

5- Evaluations
	
- **Merqury:** It runs on each 'terminal' assembly. This is, the base assembly and the resulting assembly from each branch of the pipeline. 
	
- **Busco:** It can be run only in the terminal assemblies or on all the assemblies produced by the pipeline. It uses the conda environment specified in the config as well as the parameters specified. 
	
- **Nseries:** This is run during the *finalize* on all the assemblies that are evaluated. After it, that rule combines the statistics produced by all the evaluation rules. 

# Description of all options

```
bin/create_config_assembly.py -h
usage: create_configuration_file [-h] [--configFile configFile] [--specFile specFile] [--ndconfFile ndconfFile] [--concat-cores concat_cores] [--genome-size genome_size] [--lr-type lr_type] [--basename base_name] [--species species] [--keep-intermediate]
                                 [--preprocess-lr-step PREPROCESS_ONT_STEP] [--preprocess-10X-step PREPROCESS_10X_STEP] [--preprocess-illumina-step PREPROCESS_ILLUMINA_STEP] [--flye-step FLYE_STEP] [--no-flye] [--nextdenovo-step NEXTDENOVO_STEP] [--run-nextdenovo]
                                 [--racon-cores racon_cores] [--nextpolish-cores nextpolish_cores] [--minimap2-cores minimap2_cores] [--bwa-cores bwa_cores] [--pilon-cores pilon_cores] [--medaka-cores medaka_cores] [--hypo-cores hypo_cores] [--busco-cores busco_cores]
                                 [--racon-rounds racon_rounds] [--pilon-rounds pilon_rounds] [--medaka-rounds medaka_rounds] [--nextpolish-ont-rounds nextpolish_ont_rounds] [--nextpolish-ill-rounds nextpolish_ill_rounds] [--hypo-rounds hypo_rounds]
                                 [--longranger-cores longranger_cores] [--longranger-path longranger_path] [--genomescope-opts genomescope_additional] [--no-purgedups] [--ploidy ploidy] [--run-tigmint] [--run-kraken2] [--scripts-dir SCRIPTS_DIR] [--ont-reads ONT_READS]
                                 [--ont-dir ONT_DIR] [--ont-filt ONT_FILTERED] [--pe1 PE1] [--pe2 PE2] [--processed-illumina PROCESSED_ILLUMINA] [--raw-10X RAW_10X [RAW_10X ...]] [--processed-10X PROCESSED_10X] [--10X R10X] [--illumina-dir ILLUMINA_DIR]
                                 [--assembly-in ASSEMBLY_IN [ASSEMBLY_IN ...]] [--postpolish-assemblies POSTPOLISH_ASSEMBLIES [POSTPOLISH_ASSEMBLIES ...]] [--pipeline-workdir PIPELINE_WORKDIR] [--filtlong-dir FILTLONG_DIR] [--flye-dir FLYE_DIR]
                                 [--nextdenovo-dir NEXTDENOVO_DIR] [--flye-polishing-dir POLISH_FLYE_DIR] [--nextdenovo-polishing-dir POLISH_NEXTDENOVO_DIR] [--eval-dir eval_dir] [--stats-out stats_out] [--filtlong-minlen filtlong_minlen]
                                 [--filtlong-min-mean-q filtlong_min_mean_q] [--filtlong-opts filtlong_opts] [--kraken2-db kraken2_db] [--kraken2-kmer kraken2_kmers] [--kraken2-opts additional_kraken2_opts] [--kraken2-cores kraken2_threads]
                                 [--trim-galore-opts trim_galore_opts] [--trim-Illumina-cores Trim_Illumina_cores] [--flye-cores flye_cores] [--flye-polishing-iterations flye_pol_it] [--other-flye-opts other_flye_opts] [--nextdenovo-cores nextdenovo_cores]
                                 [--nextdenovo-task nextdenovo_task] [--nextdenovo-rewrite nextdenovo_rewrite] [--nextdenovo-parallel_jobs nextdenovo_parallel_jobs] [--nextdenovo-minreadlen nextdenovo_minreadlen] [--nextdenovo-seeddepth nextdenovo_seeddepth]
                                 [--nextdenovo-seedcutoff nextdenovo_seedcutoff] [--nextdenovo-blocksize nextdenovo_blocksize] [--nextdenovo-pa-correction  nextdenovo_pa_correction] [--nextdenovo-minimap_raw nextdenovo_minimap_raw]
                                 [--nextdenovo-minimap_cns nextdenovo_minimap_cns] [--nextdenovo-minimap_map nextdenovo_minimap_map] [--nextdenovo-sort nextdenovo_sort] [--nextdenovo-correction_opts nextdenovo_correction_opts]
                                 [--nextdenovo-nextgraph_opt nextdenovo_nextgraph_opt] [--sr-cov ill_cov] [--hypo-proc hypo_processes] [--hypo-no-lr] [--hypo-opts hypo_opts] [--racon-dir racon_dir] [--racon-opts racon_opts] [--medaka-env medaka_env]
                                 [--medaka-workdir medaka_workdir] [--medaka-model medaka_model] [--medaka-consensus-opts medaka_consensus_opts] [--pilon-path pilon_path] [--pilon-opts pilon_opts] [--java-opts java_opts] [--pilon-subs pilon_subsampling]
                                 [--pilon-chunks pilon_chunks] [--purgedups-cores purgedups_cores] [--purgedups-calcuts-opts calcuts_opts] [--tigmint-cores tigmint_cores] [--tigmint-opts tigmint_opts] [--no-final-evals] [--busco-lin busco_lineage]
                                 [--merqury-db merqury_db] [--meryl-k meryl_k] [--meryl-threads meryl_threads] [--ont-list ONT_wildcards] [--illumina-list illumina_wildcards] [--r10X-list r10X_wildcards]

Create a configuration json file for the assembly pipeline.

options:
  -h, --help            show this help message and exit

General Parameters:
  --configFile configFile
                        Configuration JSON to be generated. Default assembly.config
  --specFile specFile   Cluster specifications JSON fileto be generated. Default assembly.spec
  --ndconfFile ndconfFile
                        Name pf the nextdenovo config file. Default nextdenovo.config
  --concat-cores concat_cores
                        Number of threads to concatenate reads and to run filtlong. Default 4
  --genome-size genome_size
                        Approximate genome size. Example: 615m or 2.6g. Default None
  --lr-type lr_type     Type of long reads (options are flye read-type options). Default nano-raw
  --basename base_name  Base name for the project. Default None
  --species species     Name of the species to be assembled. Default None
  --keep-intermediate   Set this to True if you do not want intermediate files to be removed. Default False
  --preprocess-lr-step PREPROCESS_ONT_STEP
                        Step for preprocessing long-reads. Default 02.1
  --preprocess-10X-step PREPROCESS_10X_STEP
                        Step for preprocessing 10X reads. Default 02.2
  --preprocess-illumina-step PREPROCESS_ILLUMINA_STEP
                        Step for preprocessing illumina reads. Default 02.2
  --flye-step FLYE_STEP
                        Step for running flye. Default 03.1
  --no-flye             Give this option if you do not want to run Flye.
  --nextdenovo-step NEXTDENOVO_STEP
                        Step for running nextdenovo. Default 03.2
  --run-nextdenovo      Give this option if you do want to run Nextdenovo.
  --racon-cores racon_cores
                        Number of threads to run the racon step. Default 16
  --nextpolish-cores nextpolish_cores
                        Number of threads to run the nextpolish step. Default 14
  --minimap2-cores minimap2_cores
                        Number of threads to run the alignment with minimap2. Default 16
  --bwa-cores bwa_cores
                        Number of threads to run the alignments with BWA-Mem2. Default 16
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
  --genomescope-opts genomescope_additional
                        Additional options to run Genomescope2 with. Default
  --no-purgedups        Give this option if you do not want to run Purgedups.
  --ploidy ploidy       Expected ploidy. Default 2
  --run-tigmint         Give this option if you want to run the scaffolding with 10X reads step.
  --run-kraken2         Give this option if you want to run Kraken2 on the input reads.

Inputs:
  --scripts-dir SCRIPTS_DIR
                        Directory with the different scripts for the pipeline. Default bin/../scripts/
  --ont-reads ONT_READS
                        File with all the ONT reads. Default None
  --ont-dir ONT_DIR     Directory where the ONT fastqs are stored. Default None
  --ont-filt ONT_FILTERED
                        File with the ONT reads after running filtlong on them. Default None
  --pe1 PE1             File with the illumina paired-end fastqs, already trimmed, pair 1.
  --pe2 PE2             File with the illumina paired-end fastqs, already trimmed, pair 2.
  --processed-illumina PROCESSED_ILLUMINA
                        Directory to Processed illumina reads. Already there or to be produced by the pipeline.
  --raw-10X RAW_10X [RAW_10X ...]
                        Dictionary with 10X raw read directories, it has to be the mkfastq dir. You must specify as well the sampleIDs from this run. Example: '{"mkfastq-dir":"sample1,sample2,sample3"}'...
  --processed-10X PROCESSED_10X
                        Directory to Processed 10X reads. Already there or to be produced by the pipeline.
  --10X R10X            File with barcoded 10X reads in fastq.gz format, concatenated.
  --illumina-dir ILLUMINA_DIR
                        Directory where the raw illumina fastqs are stored. Default None
  --assembly-in ASSEMBLY_IN [ASSEMBLY_IN ...]
                        Dictionary with assemblies that need to be polished but not assembled and directory where they should be polished. Example: '{"assembly1":"polishing_dir1"}' '{"assembly2"="polishing_dir2"}' ...
  --postpolish-assemblies POSTPOLISH_ASSEMBLIES [POSTPOLISH_ASSEMBLIES ...]
                        Dictionary with assemblies for whic postpolishing steps need to be run but that are not assembled and base step for the directory where the first postpolishing step should be run. Example: '{"assembly1":"s04.1_p03.1"}'
                        '{"assembly2":"s04.2_p03.2"}' ...

Outputs:
  --pipeline-workdir PIPELINE_WORKDIR
                        Base directory for the pipeline run. Default /software/assembly/pipelines/Assembly_pipeline/v2.0/assembly_pipeline/
  --filtlong-dir FILTLONG_DIR
                        Directory to process the ONT reads with filtlong. Default s02.1_p01.1_Filtlong
  --flye-dir FLYE_DIR   Directory to run flye. Default s03.1_p02.1_flye/
  --nextdenovo-dir NEXTDENOVO_DIR
                        Directory to run nextdenovo. Default s03.2_p02.1_nextdenovo/
  --flye-polishing-dir POLISH_FLYE_DIR
                        Directory to polish the flye assembly. Default s04.1_p03.1_polishing/
  --nextdenovo-polishing-dir POLISH_NEXTDENOVO_DIR
                        Directory to run nextdenovo. Default s04.2_p03.2_polishing/
  --eval-dir eval_dir   Base directory for the evaluations. Default evaluations/
  --stats-out stats_out
                        Path to the file with the final statistics.

Filtlong:
  --filtlong-minlen filtlong_minlen
                        Minimum read length to use with Filtlong. Default 1000
  --filtlong-min-mean-q filtlong_min_mean_q
                        Minimum mean quality to use with Filtlong. Default 80
  --filtlong-opts filtlong_opts
                        Extra options to run Filtlong (eg. -t 4000000000)

Kraken2:
  --kraken2-db kraken2_db
                        Database to be used for running Kraken2. Default None
  --kraken2-kmer kraken2_kmers
                        Database to be used for running Kraken2. Default None
  --kraken2-opts additional_kraken2_opts
                        Optional parameters for the rule Kraken2. Default
  --kraken2-cores kraken2_threads
                        Number of threads to run the Kraken2 step. Default 16

Trim_Galore:
  --trim-galore-opts trim_galore_opts
                        Optional parameters for the rule trim_galore. Default --gzip -q 20 --paired --retain_unpaired
  --trim-Illumina-cores Trim_Illumina_cores
                        Number of threads to run the Illumina trimming step. Default 4

Flye:
  --flye-cores flye_cores
                        Number of threads to run FLYE. Default 72
  --flye-polishing-iterations flye_pol_it
                        Number of polishing iterations to use with FLYE. Default 2
  --other-flye-opts other_flye_opts
                        Additional options to run Flye. Default --scaffold

Nextdenovo:
  --nextdenovo-cores nextdenovo_cores
                        Number of threads to run nextdenovo. Default 72
  --nextdenovo-task nextdenovo_task
                        Task need to run. Default all
  --nextdenovo-rewrite nextdenovo_rewrite
                        Overwrite existing directory. Default yes
  --nextdenovo-parallel_jobs nextdenovo_parallel_jobs
                        Number of tasks used to run in parallel. Default 18
  --nextdenovo-minreadlen nextdenovo_minreadlen
                        Filter reads with length < minreadlen. Default 1k
  --nextdenovo-seeddepth nextdenovo_seeddepth
                        Expected seed depth, used to calculate seed_cutoff, co-use with genome_size, you can try to set it 30-45 to get a better assembly result. Default 45
  --nextdenovo-seedcutoff nextdenovo_seedcutoff
                        Minimum seed length, <=0 means calculate it automatically using bin/seq_stat. Default 0
  --nextdenovo-blocksize nextdenovo_blocksize
                        Block size for parallel running, split non-seed reads into small files, the maximum size of each file is blocksize. Default 1g
  --nextdenovo-pa-correction  nextdenovo_pa_correction
                        number of corrected tasks used to run in parallel, each corrected task requires ~TOTAL_INPUT_BASES/4 bytes of memory usage, overwrite parallel_jobs only for this step. Default 96
  --nextdenovo-minimap_raw nextdenovo_minimap_raw
                        minimap2 options, used to find overlaps between raw reads, see minimap2-nd for details. Default -t 30
  --nextdenovo-minimap_cns nextdenovo_minimap_cns
                        minimap2 options, used to find overlaps between corrected reads. Default -t 30
  --nextdenovo-minimap_map nextdenovo_minimap_map
                        minimap2 options, used to map reads back to the assembly. Default -t 30
  --nextdenovo-sort nextdenovo_sort
                        sort options, see ovl_sort for details. Default -m 40g -t 20
  --nextdenovo-correction_opts nextdenovo_correction_opts
                        Correction options. Default -p 18
  --nextdenovo-nextgraph_opt nextdenovo_nextgraph_opt
                        nextgraph options, see nextgraph for details. Default -a 1

Hypo:
  --sr-cov ill_cov      Approximate short read coverage for hypo Default 0
  --hypo-proc hypo_processes
                        Number of contigs to be processed in parallel by HyPo. Default 6
  --hypo-no-lr          Set this to false if you don¡t want to run hypo with long reads. Default True
  --hypo-opts hypo_opts
                        Additional options to run Hypo. Default None

Racon:
  --racon-dir racon_dir
                        Directory with Racon installation. Default /scratch/project/devel/aateam/src/RACON/v1.4.21_github/
  --racon-opts racon_opts
                        Extra options to run Racon_wrapper(eg. --split 100000000) to split the assembly in chunks of specified size and decrease memory requirements. Do racon_wrapper -h for more info.

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

Purge_dups:
  --purgedups-cores purgedups_cores
                        Number of threads to run purgedups. Default 8
  --purgedups-calcuts-opts calcuts_opts
                        Adjusted values to run calcuts for purgedups. Default None

Scaffold_with_10X:
  --tigmint-cores tigmint_cores
                        Number of threads to run the 10X scaffolding step. Default 12
  --tigmint-opts tigmint_opts
                        Adjusted values to run the scaffolding with 10X reads. Default None

Finalize:
  --no-final-evals      If specified, do not run evaluations on final assemblies. Default True
  --busco-lin busco_lineage
                        Path to the lineage directory to run Busco with. Default None
  --merqury-db merqury_db
                        Meryl database. Default None
  --meryl-k meryl_k     Kmer length to build the meryl database. Default None
  --meryl-threads meryl_threads
                        Number of threads to run meryl and merqury. Default 4

Wildcards:
  --ont-list ONT_wildcards
                        List with basename of the ONT fastqs that will be used. Default None
  --illumina-list illumina_wildcards
                        List with basename of the illumina fastqs. Default None
  --r10X-list r10X_wildcards
                        List with basename of the raw 10X fastqs. Default None
```

