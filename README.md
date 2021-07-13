# assembly_pipeline
 Snakemake Pipeline used for de novo genome assembly @CNAG. It has been developed using Snakemake v6.0.5, which can be loaded in the cluster by doing:
 ```bash
 conda activate /home/devel/jgomez/conda_environments/snakemake
 ```

Currently it can take ONT reads, illumina paired-end data and illumina 10X data. It does the preprocessing of the reads, assembly, polishing and evaluations. By default it will preprocess the reads, run Flye + Hypo and evaluate the resulting assemblies with BUSCO, MERQURY and Nseries. 
It needs a config file and a spec file (json file with instructions on which resources to use in the CNAG cluster for each of the jobs). Both files are created by the script "create_config_assembly.py" that is located in the bin directory. 

List of steps currently implemented: 

# 1- Preprocessing:

	
Read concatenation
	
Longranger for 10X reads

Filtlong
	
Build meryldb (with processed 10X reads or illumina reads)
	
Concat meryldbs
	
Align ONT (Minimap2)
	
Align Illumina (BWA-MEM)

# 2- Assembly

Flye (default)
	
Nextdenovo (if turned on)

# 3- Polishing


Hypo (default)
	
Racon (if turned on)
	
Medaka (if turned on)
	
Pilon (if turned on)
	
Nextpolish ont (if turned on)
	
Nextpolish illumina (if turned on)


# 4- Evaluations
	
Merqury
	
Busco
	
Nseries

# Steps that will be included soon
	
Trim_galore (rule done, need to integrate it)


10X scaffolding (pending)
	
Shasta?

