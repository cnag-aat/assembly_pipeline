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

1- ONT reads

There are several ways of providing the ONT reads for assembly and polishing.

1.1 Using the option ``--ont-dir {DIR}`` in create_config_assembly.py:

If you do so, it will look for all the files in the directory that end in '.fastq.gz' and will add the basenames to "ONT_wildcards". These wildcards will be processed by the pipeline that will: 

-- Concatenate all the files into a single file

-- Run filtlong with the specified parameters for minlen and min_mean_q

-- Use the resulting file for assembly and/or polishing

You can also specify the basenames of the files that you want to use with the ``--ont-list `` option. In this case, the pipeline will use the wildcards that you're providing instead of merging all the files in the directory.

1.2 Using the option ```--ont-reads {FILE}``` in create_config_assembly.py:

If you do so, it will consider that you already have all the reads in one file and will:  

-- Run filtlong with the specified parameters for minlen and min_mean_q

-- Use the resulting file for assembly and/or polishing

1.3 Using the option ```--ont-filtered {FILE}```. It will use this file as the output from filtlong and will skip the preprocessing steps and directly use it for assembly and/or polishing. 



2- 10X data


# List of steps currently implemented: 

1- Preprocessing:

	
Read concatenation
	
Longranger for 10X reads

Filtlong
	
Build meryldb (with processed 10X reads or illumina reads)
	
Concat meryldbs
	
Align ONT (Minimap2)
	
Align Illumina (BWA-MEM)

2- Assembly

Flye (default)
	
Nextdenovo (if turned on)

3- Polishing


Hypo (default)
	
Racon (if turned on)
	
Medaka (if turned on)
	
Pilon (if turned on)
	
Nextpolish ont (if turned on)
	
Nextpolish illumina (if turned on)


4- Evaluations
	
Merqury
	
Busco
	
Nseries
# Steps that will be included soon
	
Trim_galore (rule done, need to integrate it)


10X scaffolding (pending)
	
Shasta?

