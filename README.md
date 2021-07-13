# assembly_pipeline
 Snakemake Pipeline used for de novo genome assembly @CNAG. It has been developed using Snakemake v6.0.5, which can be loaded in the cluster by doing:
 ```bash
 conda activate /home/devel/jgomez/conda_environments/snakemake
 ```

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
	
Filtlong (pending)
	
Trim_galore (rule done, need to integrate it)


10X scaffolding (pending)
	
Shasta?

