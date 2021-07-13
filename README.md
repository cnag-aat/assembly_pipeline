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
	
- Read concatenation
	
- Longranger for 10X reads

- Filtlong
	
- Build meryldb (with processed 10X reads or illumina reads)
	
- Concat meryldbs
	
- Align ONT (Minimap2)
	
- Align Illumina (BWA-MEM)

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

# Steps that will be included soon
	
Trim_galore (rule done, need to integrate it)


10X scaffolding (pending)
	
Shasta?

