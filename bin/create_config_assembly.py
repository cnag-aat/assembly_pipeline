#!/usr/bin/env python3
import os
import json
import argparse
import sys
import re
from configparser import ConfigParser
nextdenovo_config = ConfigParser()

def get_wildcards(dir, wildcards, ext):
  for r, d, f in os.walk(dir):
    for file in f:
      if re.search(ext, file):
        a = file.replace(ext,'')
        if wildcards == None:
          wildcards = a
        else:
          wildcards += "," + a
  return wildcards

class CreateConfigurationFile(object):
    """Class which manages Configuration file Manager"""
      
    def __init__(self):
        """Class constructor"""
        #GENERAL PARAMETERS
        self.configFile = "assembly.config"                                                       #Name of the json configuration file to be created
        self.specFile = "assembly.spec"                                                           #Name of the spec file to be created
        self.ndconfFile = "nextdenovo.config"                                                     #Name of the nextdenovo config file to be created
        self.concat_cores = 4                                                                     #Number of threads to concatenate the reads and to run filtlong
        self.keepintermediate = False                                                             #Set this to True if you do not want intermediate files to be removed
        self.lr_type  = "nano-raw"                                                                #Type of long reads (options are flye read-type options)
        self.base_name = None                                                                     #Base name for the project
        self.species = None                                                                       #Name of the species to be assembled
        self.genome_size = None							                                                      #Estimated genome size
        self.preprocess_ont_step = "02.1"                                                         #Step directory for preprocessing ont
        self.preprocess_10X_step = "02.2"
        self.preprocess_illumina_step = "02.2"
        self.flye_step = "03.1"                                                                   #Step direcotory for running flye
        self.nextdenovo_step = "03.2"                                                             #Step direcotory for running nextdenovo
        self.run_flye = True       
        self.run_nextdenovo = False   
        self.racon_rounds = 0                                                                     #Number of rounds of racon to run
        self.pilon_rounds = 0                                                                     #Number of rounds of pilon to run
        self.medaka_rounds = 0                                                                    #Number of rounds of medaka to run  
        self.nextpolish_ont_rounds = 0                                                            #Number of rounds for running Nexpolish with ONT
        self.nextpolish_ill_rounds = 0                                                            #Number of rounds for running Nexpolish with illumina
        self.hypo_rounds = 1                                                                      #Number of rounds for running hypo
        self.run_purgedups = True 
        self.run_tigmint = False
        self.minimap2_cores = 16                                                                  #Number of threads to run the alignment with minimap2
        self.bwa_cores = 16                                                                       #Number of threads to run the alignment for the pilon step
        self.racon_cores = 16                                                                     #Number of threads to run the racon step
        self.nextpolish_cores = 14                                                                #Number of threads to run the nextpolish step
        self.pilon_cores = 16                                                                     #Number of threads to run the pilon step
        self.medaka_cores = 16                                                                    #Number of threads to tun the medaka step
        self.hypo_cores = 24                                                                      #Number of threads to tun the hypo step
        self.busco_cores = 16                                                                     #Number of threads to tun the BUSCO    
        self.longranger_cores = 8                                                                 #Number of threads to run longranger   
        self.longranger_path = "/scratch/project/devel/aateam/src/10X/longranger-2.2.2" 
        self.genomescope_additional = ""  
        self.ploidy = 2
        self.run_kraken2 = False

        #ALL SPEC PARAMETERS
        self.all_qos = "test"
        self.all_time = "00:05:00"
        self.all_queue = "genD"

        #LONGRANGER SPEC PARAMETERS
        self.longranger_qos = "normal"
        self.longranger_time = "12:00:00"
        self.longranger_queue = "genD"
        self.longranger_mem = "50G"

        #TRIMGALORE PARAMETERS
        self.trim_galore_opts = "--gzip -q 20 --paired --retain_unpaired"
        self.Trim_Illumina_cores = 4                                                              #Number of threads to run the trim Illumina step

        #TRIMGALORE SPEC PARAMETERS
        self.trimgalore_qos = "normal"
        self.trimgalore_time = "3:00:00"
        self.trimgalore_queue = "genD"
        self.trimgalore_mem = "50000"

        #CONCAT READS SPEC PARAMETERS
        self.concat_reads_qos = "normal"
        self.concat_reads_time = "10:00:00"
        self.concat_reads_queue = "genD"
        self.concat_reads_mem = "5000"

        #NANOPLOT SPEC PARAMETERS
        self.nanoplot_qos = "normal"
        self.nanoplot_time = "6:00:00"
        self.nanoplot_queue = "genD"
        self.nanoplot_mem = "10000"

        #KRAKEN PARAMETERS
        self.kraken2_db = None  
        self.kraken2_kmers = None
        self.kraken2_threads = 16
        self.additional_kraken2_opts = ""

        #KRAKEN SPEC PARAMETERS
        self.kraken2_qos = "vlong"
        self.kraken2_time = "48:00:00"
        self.kraken2_queue = "genD"
        self.kraken2_mem = "1000"

        #BUILD MERYL SPEC PARAMETERS
        self.build_meryl_qos = "normal"
        self.build_meryl_time = "6:00:00"
        self.build_meryl_queue = "genD"
        self.build_meryl_mem = "50G"

        #CONCAT MERYL SPEC PARAMETERS
        self.concat_meryl_qos = "normal"
        self.concat_meryl_time = "6:00:00"
        self.concat_meryl_queue = "genD"
        self.concat_meryl_mem = "10000"


        #GENOMESCOPE2 SPEC PARAMETERS
        self.genomescope_qos = "short"
        self.genomescope_time = "1:00:00"
        self.genomescope_queue = "genD"
        self.genomescope_mem = "100G"

        #INPUT PARAMETERS
        self.scripts_dir = os.path.dirname(sys.argv[0]) + "/../scripts/"                          #Directory with the different scripts for the pipeline
        self.ONT_reads = None                                                                     #File with all the ONT reads    
        self.ONT_dir = None                                                                       #Directory with the ONT reads, give this option if you don't have a single fastq file with all the reads
        self.ONT_filtered = None                                                                  #File with the ONT reads after running filtlong         
        self.pe1 = None                                                                           #File with the illumina paired-end fastqs, already trimeed, pair 1
        self.pe2 = None                                                                           #File with the illumina paired-end fastqs, already trimmed, pair 2
        self.r10X = None                                                                          #File with barcoded 10X reads in fastq.gz format, concatenated
        self.raw_10X = {}                                                                         #List with 10X raw read directories, it has to be the mkfastq dir. You must specify as well the sampleIDs from this run. Example: 
        self.processed_10X = None                                                                 #Directory to Processed 10X reads, already there or to be produced by the pipeline
        self.illumina_dir = None                                                                  #Directory with the illumina raw reads, give this option if you don't have a single fastq file with all the reads
        self.processed_illumina = None                                                            #Directory to Processed illumina reads, already there or to be produced by the pipeline
        self.assembly_in = {}                                                                     #List of input assemblies that need to be polished but are not assembled by the pipeline
        self.assemblies = {}
        self.postpolish_assemblies = {}                                                           #List of input assemblies for which postpolishing steps need to be run but are not produced by the pipeline
        self.assemblies_cur = {}
        self.r10X_reads = {}

        #OUTPUT PARAMETERS
        self.pipeline_workdir = os.getcwd() + "/"                                                 #Base directory for the pipeline run
        self.filtlong_dir = "s" + self.preprocess_ont_step + "_p01.1_Filtlong"                    #Directory to process the ONT reads
        self.flye_dir = "s" + self.flye_step + "_p" + self.preprocess_ont_step + "_flye/"         #Directory to run flye 
        self.nextdenovo_dir =  "s" + self.nextdenovo_step + "_p" + self.preprocess_ont_step + "_nextdenovo/"         #Directory to run Nextdenovo 
        self.flye_out = self.flye_dir + "flye.assembly.fasta"
        self.nextdenovo_out = self.nextdenovo_dir + "nextdenovo.asssembly.fasta"
        self.polish_flye_dir = "s04.1_p" + self.flye_step + "_polishing/"                          #Base directory to run polishing pipeline in flye assembly
        self.polish_nextdenovo_dir = "s04.2_p" + self.nextdenovo_step + "_polishing/"              #Base directory to run polishing pipeline in nextdenovo assembly  
        self.eval_dir = "evaluations/"                                                             #Base directory for the evaluations
        self.stats_out = None                                                                      #Path to the file with the final pipeline statistics

        #FILTLONG PARAMETERS
        self.filtlong_minlen = "1000"
        self.filtlong_min_mean_q = "80"
        self.filtlong_opts = None

        #FILTLONG SPEC PARAMETERS
        self.filtlong_qos = "normal"
        self.filtlong_time = "10:00:00"
        self.filtlong_queue = "genD"
        self.filtlong_mem = "10000"

        #FLYE PARAMETERS
        self.flye_cores = 128	                                                                  #Number of threads to run Flye
        self.flye_pol_it = 2				      			                  #Number of polishing iterations to use with FLYE
        self.other_flye_opts = " --scaffold "                                                     #include here genome size in pipeline											

        #FLYE SPEC PARAMETERS
        self.flye_qos = "vlong"
        self.flye_time = "48:00:00"
        self.flye_queue = "genD"
        self.flye_mem = "900G"

        #NEXTDENOVO PARAMETERS
        self.nextdenovo_cores = 128	                                                          #Number of threads to run nextdenovo        
        self.nextdenovo_task = "all"
        self.nextdenovo_rewrite = "yes"
        self.nextdenovo_parallel_jobs = 18
        self.nextdenovo_minreadlen = "1k"
        self.nextdenovo_seeddepth = 45
        self.nextdenovo_seedcutoff = 0
        self.nextdenovo_blocksize = "1g"
        self.nextdenovo_pa_correction = 96
        self.nextdenovo_minimap_raw = "-t 30"
        self.nextdenovo_sort = "-m 40g -t 20"
        self.nextdenovo_correction_opts = "-p 18"                      
        self.nextdenovo_minimap_cns = "-t 30 "
        self.nextdenovo_minimap_map = "-t 30 "              
        self.nextdenovo_nextgraph_opt = "-a 1"

        #NEXTDENOVO SPEC PARAMETERS
        self.nextdenovo_qos = "vlong"
        self.nextdenovo_time = "48:00:00"
        self.nextdenovo_queue = "genD"
        self.nextdenovo_mem = "900G"

        #MINIMAP2 SPEC PARAMETERS
        self.minimap_qos = "normal"
        self.minimap_time = "6:00:00"
        self.minimap_queue = "genD"
        self.minimap_mem = "500G"
        
        #BWA SPEC PARAMETERS
        self.bwa_qos = "normal"
        self.bwa_time = "6:00:00"
        self.bwa_queue = "genD"
        self.bwa_mem = "100G"

        #HYPO PARAMETERS
        self.ill_cov = 0                                                                          #Approximate short read coverage for hypo
        self.hypo_processes = 6                                                                   #Number of contigs to be processed in parallel by hypo
        self.hypo_lr = True                                                                       #Set this to true if you want to run hypo with both long and short reads
        self.hypo_opts = None                                                                     #Extra options to run Hypo 

        #HYPO SPEC PARAMETERS
        self.hypo_qos = "normal"
        self.hypo_time = "6:00:00"
        self.hypo_queue = "genD"
        self.hypo_mem = "250G"

        #RACON PARAMETERS
        self.racon_dir = "/scratch/project/devel/aateam/src/RACON/v1.4.21_github/"
        self.racon_opts = None                                                                    #Extra options to run Racon_wrapper

        #RACON SPEC PARAMETERS
        self.racon_qos = "normal"
        self.racon_time = "12:00:00"
        self.racon_queue = "main"

        #MEDAKA PARAMETERS
        self.medaka_env = "/scratch/project/devel/aateam/src/MEDAKA/medaka-141" 
        self.medaka_workdir = "$TMPDIR/"      
        self.medaka_model = None                                                                  #User needs to specify the model that will be used to run medaka
        self.medaka_consensus_opts = None                                                         #Specify any parameters to change when running medaka consensus

        #MEDAKA SPEC PARAMETERS
        self.medaka_qos = "normal"
        self.medaka_time = "16:00:00"
        self.medaka_queue = "main"

        #PILON PARAMETERS
        self.pilon_path = "/apps/PILON/1.21/pilon"                                                #Path to Pilon executable
        self.pilon_opts = "--fix bases --changes "                                                #Options to run pilon with
        self.java_opts = None                                                                     #Options for the java execution of pilon
        self.pilon_subsampling = 0                                                                #Percentage of reads to randomly use when running pilon      
        self.pilon_chunks = 25                                                                    #Number of chunks to split pilon polishing  

        #SPLIT PILON SPEC PARAMETERS
        self.split_pilon_qos = "normal"
        self.split_pilon_time = "6:00:00"
        self.split_pilon_queue = "genB,main"

        #PILON SPEC PARAMETERS
        self.pilon_qos = "normal"
        self.pilon_time = "12:00:00"
        self.pilon_queue = "main"
        self.pilon_threads = self.pilon_cores

        #JOIN PILON SPEC PARAMETERS
        self.join_pilon_qos = "normal"
        self.join_pilon_time = "6:00:00"
        self.join_pilon_queue = "genB,main"

        #NEXTPOLISH LR SPEC PARAMETERS
        self.nextpolish_lr_qos = "normal"
        self.nextpolish_lr_time = "6:00:00"
        self.nextpolish_lr_queue = "genD"
        self.nextpolish_lr_mem = "500"
        
        #NEXTPOLISH SR SPEC PARAMETERS
        self.nextpolish_sr_qos = "normal"
        self.nextpolish_sr_time = "6:00:00"
        self.nextpolish_sr_queue = "genD"
        self.nextpolish_sr_mem = "2000"
        
        #PURGEDUPS PARAMETERS
        self.purgedups_cores = 8 
        self.calcuts_opts = None                                                                  #Adjusted values to run calcuts for purgedups

        #PURGEDUPS SPEC PARAMETERS
        self.purgedups_qos = "normal"
        self.purgedups_time = "1:00:00"
        self.purgedups_queue = "genD"
        self.purgedups_mem = "200G"

        #10X SCAFFOLDING PARAMETERS
        self.tigmint_opts = None                                                                  #Additional option to give to the 10X scaffolding step
        self.tigmint_cores = 12

        #10X SCAFFOLDING SPEC PARAMETERS
        self.tigmint_qos = "long"
        self.tigmint_time = "24:00:00"
        self.tigmint_queue = "genD"
        self.tigmint_mem = "100G"   

        #FINALIZE PARAMETERS
        self.final_evals = True                                                                  #Set this to true if you want evaluations to be run on each of the final assemblies     
        self.busco_lineage = None                                                                 #Path to the lineage directory to run Busco with
        self.merqury_db = None
        self.meryl_k = None
        self.meryl_threads = 4

        #STATS SPEC PARAMETERS
        self.stats_qos = "test"
        self.stats_time = "0:10:00"
        self.stats_queue = "genD"  
        self.stats_mem = "1000"
        
        #BUSCO SPEC PARAMETERS
        self.busco_qos = "short"
        self.busco_time = "6:00:00"
        self.busco_queue = "genD"  
        self.busco_mem = "50G"

        #MERQURY SPEC PARAMETERS
        self.merq_qos = "normal"
        self.merq_time = "3:00:00"
        self.merq_queue = "genD"
        self.merq_mem = "100G"          

        #FINALIZE SPEC PARAMETERS
        self.fin_qos = "short"
        self.fin_time = "2:00:00"
        self.fin_queue = "genD"
        self.fin_mem = "1000" 

        #WILDCARDS
        self.ONT_wildcards = None
        self.illumina_wildcards = None
        self.r10X_wildcards = None                                                                #For raw 10X we need to give this argument, for processed 10X reads, the pipeline can obtain it

###
        #DICTIONARIES
        self.allParameters = {}
        self.generalParameters = {}
        self.allSpecParameters = {}
        self.inputParameters = {}
        self.outputParameters = {}
        self.longrangerSpecParameters = {}
        self.trimgaloreParameters = {}
        self.trimgaloreSpecParameters = {}
        self.concatreadsSpecParameters = {}
        self.nanoplotSpecParameters = {}
        self.kraken2Parameters = {}
        self.kraken2SpecParameters = {}
        self.buildmerylSpecParameters = {}
        self.concatmerylSpecParameters = {}
        self.genomescopeSpecParameters = {}
        self.filtlongParameters = {}
        self.filtlongSpecParameters = {}
        self.flyeParameters = {}
        self.flyeSpecParameters = {}
        self.nextdenovoParameters = {}
        self.nextdenovoSpecParameters = {}
        self.nextdenovoConfig = {}
        self.minimapSpecParameters = {}
        self.bwaSpecParameters= {}
        self.hypoParameters = {}
        self.hypoSpecParameters = {}
        self.raconParameters = {}
        self.raconSpecParameters = {}
        self.medakaParameters = {}
        self.medakaSpecParameters = {}
        self.pilonParameters = {}
        self.splitpilonSpecParameters = {}
        self.pilonSpecParameters = {}
        self.joinpilonSpecParameters = {}
        self.nextpolishlrSpecParameters = {}
        self.nextpolishsrSpecParameters = {}
        self.purgedupsParameters = {}
        self.purgedupsSpecParameters = {}
        self.scaffold10XParameters = {}
        self.scaffold10XSpecParameters = {}
        self.finalizeParameters = {}
        self.statsSpecParameters = {}
        self.buscoSpecParameters = {}
        self.merqSpecParameters = {}
        self.finSpecParameters = {}
        self.repSpecParameters = {}
        self.wildcardParameters = {}

####

    def register_parameter(self, parser):
        """Register all parameters with the given
        argparse parser"""
        self.register_general(parser)
        self.register_input(parser)
        self.register_output(parser)
        self.register_filtlong(parser)
        self.register_kraken2(parser)
        self.register_trimgalore(parser)
        self.register_flye(parser)
        self.register_nextdenovo(parser)
        self.register_hypo(parser)
        self.register_racon(parser)
        self.register_medaka(parser)
        self.register_pilon(parser)
        self.register_purgedups(parser)
        self.register_scaffold10X(parser)
        self.register_finalize(parser)
        self.register_wildcards(parser)

    def register_general(self, parser):
        """Register all general parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        general_group = parser.add_argument_group('General Parameters')
        general_group.add_argument('--configFile', dest="configFile", metavar="configFile", help='Configuration JSON to be generated. Default %s' % self.configFile)
        general_group.add_argument('--specFile', dest="specFile", metavar="specFile", help='Cluster specifications JSON  fileto be generated. Default %s' % self.specFile)
        general_group.add_argument('--ndconfFile', dest="ndconfFile", metavar="ndconfFile", help='Name pf the nextdenovo config file. Default %s' % self.ndconfFile)
        general_group.add_argument('--concat-cores', type = int, dest="concat_cores", metavar="concat_cores", default=self.concat_cores, help='Number of threads to concatenate reads and to run filtlong. Default %s' % self.concat_cores)
        general_group.add_argument('--genome-size', dest="genome_size", metavar="genome_size", help='Approximate genome size. Example: 615m or 2.6g. Default %s' % self.genome_size)
        general_group.add_argument('--lr-type', dest="lr_type", metavar="lr_type", default=self.lr_type, choices=['pacbio-raw', 'pacbio-corr', 'pacbio-hifi', 'nano-raw', 'nano-corr', 'nano-hq', 'subassemblies'],  help='Type of long reads (options are flye read-type options). Default %s' % self.lr_type)
        general_group.add_argument('--basename', dest="base_name", metavar="base_name", help='Base name for the project. Default %s' % self.base_name)
        general_group.add_argument('--species', dest="species", metavar="species", help='Name of the species to be assembled. Default %s' % self.species)
        general_group.add_argument('--keep-intermediate', dest="keepintermediate", action="store_true", help='Set this to True if you do not want intermediate files to be removed. Default %s' % self.keepintermediate)
        general_group.add_argument('--preprocess-lr-step', dest="preprocess_ont_step", default=self.preprocess_ont_step, help='Step for preprocessing long-reads. Default %s' % self.preprocess_ont_step)
        general_group.add_argument('--preprocess-10X-step', dest="preprocess_10X_step", default=self.preprocess_10X_step, help='Step for preprocessing 10X reads. Default %s' % self.preprocess_10X_step)
        general_group.add_argument('--preprocess-illumina-step', dest="preprocess_illumina_step", default=self.preprocess_illumina_step, help='Step for preprocessing illumina reads. Default %s' % self.preprocess_illumina_step)
        general_group.add_argument('--flye-step', dest="flye_step", default=self.flye_step, help='Step for running flye. Default %s' % self.flye_step)
        general_group.add_argument('--no-flye', dest="run_flye", action="store_false", help='Give this option if you do not want to run Flye.')
        general_group.add_argument('--nextdenovo-step', dest="nextdenovo_step", default=self.nextdenovo_step, help='Step for running nextdenovo. Default %s' % self.nextdenovo_step)
        general_group.add_argument('--run-nextdenovo', dest="run_nextdenovo", action="store_true", help='Give this option if you do want to run Nextdenovo.')
        general_group.add_argument('--racon-cores', type = int, dest="racon_cores", metavar="racon_cores", default=self.racon_cores, help='Number of threads to run the racon step. Default %s' % self.racon_cores)
        general_group.add_argument('--nextpolish-cores', type = int, dest="nextpolish_cores", metavar="nextpolish_cores", default=self.nextpolish_cores, help='Number of threads to run the nextpolish step. Default %s' % self.nextpolish_cores)
        general_group.add_argument('--minimap2-cores', type = int, dest="minimap2_cores", metavar="minimap2_cores", default=self.minimap2_cores, help='Number of threads to run the alignment with minimap2. Default %s' % self.minimap2_cores)
        general_group.add_argument('--bwa-cores', type = int, dest="bwa_cores", metavar="bwa_cores", default=self.bwa_cores, help='Number of threads to run the alignments with BWA-Mem2. Default %s' % self.bwa_cores)
        general_group.add_argument('--pilon-cores', type = int, dest="pilon_cores", metavar="pilon_cores", default=self.pilon_cores, help='Number of threads to run the pilon step. Default %s' % self.pilon_cores)
        general_group.add_argument('--medaka-cores', type = int, dest="medaka_cores", metavar="medaka_cores", default=self.medaka_cores, help='Number of threads to run the medaka step. Default %s' % self.medaka_cores)
        general_group.add_argument('--hypo-cores', type = int, dest="hypo_cores", metavar="hypo_cores", default=self.hypo_cores, help='Number of threads to run the hypo step. Default %s' % self.hypo_cores)
        general_group.add_argument('--busco-cores', type = int, dest="busco_cores", metavar="busco_cores", default=self.busco_cores, help='Number of threads to run BUSCO. Default %s' % self.busco_cores)
        general_group.add_argument('--racon-rounds', type = int, dest="racon_rounds", metavar="racon_rounds", default=self.racon_rounds, help='Number of rounds of racon to run. Default %s' % self.racon_rounds)
        general_group.add_argument('--pilon-rounds', type = int, dest="pilon_rounds", metavar="pilon_rounds", default=self.pilon_rounds, help='Number of rounds of pilon to run. Default %s' % self.pilon_rounds)
        general_group.add_argument('--medaka-rounds', type = int, dest="medaka_rounds", metavar="medaka_rounds", default=self.medaka_rounds, help='Number of rounds of medaka to run. Default %s' % self.medaka_rounds)
        general_group.add_argument('--nextpolish-ont-rounds', type = int, dest="nextpolish_ont_rounds", metavar="nextpolish_ont_rounds", default=self.nextpolish_ont_rounds, help='Number of rounds to run the Nextpolish with ONT step. Default %s' % self.nextpolish_ont_rounds)
        general_group.add_argument('--nextpolish-ill-rounds', type = int, dest="nextpolish_ill_rounds", metavar="nextpolish_ill_rounds", default=self.nextpolish_ill_rounds, help='Number of rounds to run the Nextpolish with illumina step. Default %s' % self.nextpolish_ill_rounds)
        general_group.add_argument('--hypo-rounds', type = int, dest="hypo_rounds", metavar="hypo_rounds", default=self.hypo_rounds, help='Number of rounds to run the Hypostep. Default %s' % self.hypo_rounds)
        general_group.add_argument('--longranger-cores', type = int, dest="longranger_cores", metavar="longranger_cores", default=self.longranger_cores, help='Number of threads to run longranger. Default %s' % self.longranger_cores)
        general_group.add_argument('--longranger-path', dest="longranger_path", metavar="longranger_path", help='Path to longranger executable. Default %s' % self.longranger_path)
        general_group.add_argument('--genomescope-opts', dest="genomescope_additional", metavar="genomescope_additional", help='Additional options to run Genomescope2 with. Default %s' % self.genomescope_additional)
        general_group.add_argument('--no-purgedups', dest="run_purgedups", action="store_false", help='Give this option if you do not want to run Purgedups.')
        general_group.add_argument('--ploidy', type = int, dest="ploidy", metavar="ploidy", default=self.ploidy, help='Expected ploidy. Default %s' % self.ploidy) 
        general_group.add_argument('--run-tigmint', dest="run_tigmint", action="store_true", help='Give this option if you want to run the scaffolding with 10X reads step.')
        general_group.add_argument('--run-kraken2', dest="run_kraken2", action="store_true", help='Give this option if you want to run Kraken2 on the input reads.')

    def register_input(self, parser):
        """Register all input parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        input_group = parser.add_argument_group('Inputs')
        input_group.add_argument('--scripts-dir', dest="scripts_dir", help='Directory with the different scripts for the pipeline. Default %s' % self.scripts_dir)
        input_group.add_argument('--ont-reads', dest="ONT_reads", help='File with all the ONT reads. Default %s' % self.ONT_reads)
        input_group.add_argument('--ont-dir', dest="ONT_dir", help='Directory where the ONT fastqs are stored. Default %s' % self.ONT_dir)
        input_group.add_argument('--ont-filt', dest="ONT_filtered", help='File with the ONT reads after running filtlong on them. Default %s' % self.ONT_filtered)
        input_group.add_argument('--pe1', dest="pe1", help='File with the illumina paired-end fastqs, already trimmed,  pair 1.')
        input_group.add_argument('--pe2', dest="pe2", help='File with the illumina paired-end fastqs, already trimmed, pair 2.')
        input_group.add_argument('--processed-illumina', dest="processed_illumina", help='Directory to Processed illumina reads. Already there or to be produced by the pipeline.')
        input_group.add_argument('--raw-10X', dest="raw_10X",  nargs="+", type=json.loads, default=self.raw_10X, help='Dictionary with 10X raw read directories, it has to be the mkfastq dir. You must specify as well the sampleIDs from this run. Example: \'{\"mkfastq-dir\":\"sample1,sample2,sample3\"}\'...')
        input_group.add_argument('--processed-10X', dest="processed_10X", help='Directory to Processed 10X reads. Already there or to be produced by the pipeline.')
        input_group.add_argument('--10X', dest="r10X", help='File with barcoded 10X reads in fastq.gz format, concatenated.')
        input_group.add_argument('--illumina-dir', dest="illumina_dir", help='Directory where the raw illumina fastqs are stored. Default %s' % self.illumina_dir)
        input_group.add_argument('--assembly-in', dest="assembly_in", nargs="+", type=json.loads, default=self.assembly_in, help='Dictionary with assemblies that need to be polished but not assembled and directory where they should be polished. Example: \'{\"assembly1\":\"polishing_dir1\"}\' \'{\"assembly2\"=\"polishing_dir2\"}\' ...')
        input_group.add_argument('--postpolish-assemblies', dest="postpolish_assemblies", nargs="+", type=json.loads, default=self.postpolish_assemblies, help='Dictionary with assemblies for whic postpolishing steps need to be run but that are not assembled and base step for the directory where the first postpolishing step should be run. Example: \'{\"assembly1\":\"s04.1_p03.1\"}\' \'{\"assembly2\":\"s04.2_p03.2\"}\' ...')

    def register_output(self, parser):
        """Register all output parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        output_group = parser.add_argument_group('Outputs')
        output_group.add_argument('--pipeline-workdir', dest="pipeline_workdir", help='Base directory for the pipeline run. Default %s' % self.pipeline_workdir)
        output_group.add_argument('--filtlong-dir', dest="filtlong_dir",  help='Directory to process the ONT reads with filtlong. Default %s' % self.filtlong_dir)
        output_group.add_argument('--flye-dir', dest="flye_dir",  help='Directory to run flye. Default %s' % self.flye_dir)
        output_group.add_argument('--nextdenovo-dir', dest="nextdenovo_dir",  help='Directory to run nextdenovo. Default %s' % self.nextdenovo_dir)
        output_group.add_argument('--flye-polishing-dir', dest="polish_flye_dir",  help='Directory to polish the flye assembly. Default %s' % self.polish_flye_dir)
        output_group.add_argument('--nextdenovo-polishing-dir', dest="polish_nextdenovo_dir",  help='Directory to run nextdenovo. Default %s' % self.polish_nextdenovo_dir)
        output_group.add_argument('--eval-dir', dest="eval_dir", metavar="eval_dir",  help='Base directory for the evaluations. Default %s' %self.eval_dir)
        output_group.add_argument('--stats-out', dest="stats_out", metavar="stats_out",  help='Path to the file with the final statistics.')

    def register_filtlong(self, parser):
        """Register all filtlong parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        filtlong_group = parser.add_argument_group('Filtlong')
        filtlong_group.add_argument('--filtlong-minlen', dest="filtlong_minlen", metavar="filtlong_minlen", default=self.filtlong_minlen, type = int, help='Minimum read length to use with Filtlong. Default %s' % self.filtlong_minlen)
        filtlong_group.add_argument('--filtlong-min-mean-q', dest="filtlong_min_mean_q", metavar="filtlong_min_mean_q", default=self.filtlong_min_mean_q, type = int, help='Minimum mean quality to use with Filtlong. Default %s' % self.filtlong_min_mean_q)
        filtlong_group.add_argument('--filtlong-opts', dest="filtlong_opts", metavar="filtlong_opts", default=self.filtlong_opts, help='Extra options to run Filtlong (eg. -t 4000000000)')

    def register_trimgalore(self, parser):
        """Register all trimgalore parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        trimgalore_group = parser.add_argument_group('Trim_Galore')
        trimgalore_group.add_argument('--trim-galore-opts', dest="trim_galore_opts", metavar="trim_galore_opts", default=self.trim_galore_opts, help='Optional parameters for the rule trim_galore. Default %s' % self.trim_galore_opts)
        trimgalore_group.add_argument('--trim-Illumina-cores', type = int, dest="Trim_Illumina_cores", metavar="Trim_Illumina_cores", default=self.Trim_Illumina_cores, help='Number of threads to run the Illumina trimming step. Default %s' % self.Trim_Illumina_cores)

    def register_kraken2(self, parser):
        """Register all Kraken2 parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        kraken2_group = parser.add_argument_group('Kraken2')
        kraken2_group.add_argument('--kraken2-db', dest="kraken2_db", metavar="kraken2_db", default=self.kraken2_db, help='Database to be used for running Kraken2. Default %s' % self.kraken2_db)
        kraken2_group.add_argument('--kraken2-kmer', dest="kraken2_kmers", metavar="kraken2_kmers", default=self.kraken2_kmers, help='Database to be used for running Kraken2. Default %s' % self.kraken2_kmers)
        kraken2_group.add_argument('--kraken2-opts', dest="additional_kraken2_opts", metavar="additional_kraken2_opts", default=self.additional_kraken2_opts, help='Optional parameters for the rule Kraken2. Default %s' % self.additional_kraken2_opts)
        kraken2_group.add_argument('--kraken2-cores', type = int, dest="kraken2_threads", metavar="kraken2_threads", default=self.kraken2_threads, help='Number of threads to run the Kraken2 step. Default %s' % self.kraken2_threads)

    def register_flye(self, parser):
        """Register all flye parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        flye_group = parser.add_argument_group('Flye')
        flye_group.add_argument('--flye-cores', dest="flye_cores", metavar="flye_cores", default=self.flye_cores, type = int, help='Number of threads to run FLYE. Default %s' % self.flye_cores)
        flye_group.add_argument('--flye-polishing-iterations', dest="flye_pol_it", metavar="flye_pol_it", default=self.flye_pol_it, type = int, help='Number of polishing iterations to use with FLYE. Default %s' % self.flye_pol_it)
        flye_group.add_argument('--other-flye-opts', dest="other_flye_opts", metavar="other_flye_opts", default=self.other_flye_opts, help='Additional options to run Flye. Default %s' % self.other_flye_opts)

    def register_nextdenovo(self, parser):
        """Register all nextdenovo parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        nextdenovo_group = parser.add_argument_group('Nextdenovo')
        nextdenovo_group.add_argument('--nextdenovo-cores', dest="nextdenovo_cores", metavar="nextdenovo_cores", default=self.nextdenovo_cores, type = int, help='Number of threads to run nextdenovo. Default %s' % self.nextdenovo_cores)
        nextdenovo_group.add_argument('--nextdenovo-task', dest="nextdenovo_task", metavar="nextdenovo_task", choices=['all', 'correct', 'assemble'], default=self.nextdenovo_task, help='Task need to run. Default %s' % self.nextdenovo_task)
        nextdenovo_group.add_argument('--nextdenovo-rewrite', dest="nextdenovo_rewrite", metavar="nextdenovo_rewrite", choices=['yes', 'no'], default=self.nextdenovo_rewrite, help='Overwrite existing directory. Default %s' % self.nextdenovo_rewrite)
        nextdenovo_group.add_argument('--nextdenovo-parallel_jobs', dest="nextdenovo_parallel_jobs", metavar="nextdenovo_parallel_jobs", default=self.nextdenovo_parallel_jobs, type = int, help='Number of tasks used to run in parallel. Default %s' % self.nextdenovo_parallel_jobs)
    #    nextdenovo_group.add_argument('--nextdenovo-tmp', dest="nextdenovo_tmp", metavar="nextdenovo_tmp", default=self.nextdenovo_tmp, help='Temporary directory in compute nodes to avoid high IO wait. Default %s' % self.nextdenovo_tmp) --> it cannot be used with local
        nextdenovo_group.add_argument('--nextdenovo-minreadlen', dest="nextdenovo_minreadlen", metavar="nextdenovo_minreadlen", default=self.nextdenovo_minreadlen, help='Filter reads with length < minreadlen. Default %s' % self.nextdenovo_minreadlen)
        nextdenovo_group.add_argument('--nextdenovo-seeddepth', dest="nextdenovo_seeddepth", metavar="nextdenovo_seeddepth", default=self.nextdenovo_seeddepth, type = int, help='Expected seed depth, used to calculate seed_cutoff, co-use with genome_size, you can try to set it 30-45 to get a better assembly result. Default %s' % self.nextdenovo_seeddepth)
        nextdenovo_group.add_argument('--nextdenovo-seedcutoff', dest="nextdenovo_seedcutoff", metavar="nextdenovo_seedcutoff", default=self.nextdenovo_seedcutoff, type = int, help='Minimum seed length, <=0 means calculate it automatically using bin/seq_stat. Default %s' % self.nextdenovo_seedcutoff)
        nextdenovo_group.add_argument('--nextdenovo-blocksize', dest="nextdenovo_blocksize", metavar="nextdenovo_blocksize", default=self.nextdenovo_blocksize, help='Block size for parallel running, split non-seed reads into small files, the maximum size of each file is blocksize. Default %s' % self.nextdenovo_blocksize)
        nextdenovo_group.add_argument('--nextdenovo-pa-correction ', dest="nextdenovo_pa_correction", metavar="nextdenovo_pa_correction", default=self.nextdenovo_pa_correction, type = int, help='number of corrected tasks used to run in parallel, each corrected task requires ~TOTAL_INPUT_BASES/4 bytes of memory usage, overwrite parallel_jobs only for this step. Default %s' % self.nextdenovo_pa_correction )
        nextdenovo_group.add_argument('--nextdenovo-minimap_raw', dest="nextdenovo_minimap_raw", metavar="nextdenovo_minimap_raw", default=self.nextdenovo_minimap_raw , help='minimap2 options, used to find overlaps between raw reads, see minimap2-nd for details. Default %s' % self.nextdenovo_minimap_raw )
        nextdenovo_group.add_argument('--nextdenovo-minimap_cns', dest="nextdenovo_minimap_cns", metavar="nextdenovo_minimap_cns", default=self.nextdenovo_minimap_cns , help='minimap2 options, used to find overlaps between corrected reads. Default %s' % self.nextdenovo_minimap_cns )
        nextdenovo_group.add_argument('--nextdenovo-minimap_map', dest="nextdenovo_minimap_map", metavar="nextdenovo_minimap_map", default=self.nextdenovo_minimap_map , help='minimap2 options, used to map reads back to the assembly. Default %s' % self.nextdenovo_minimap_map)
        nextdenovo_group.add_argument('--nextdenovo-sort', dest="nextdenovo_sort", metavar="nextdenovo_sort", default=self.nextdenovo_sort, help='sort options, see ovl_sort for details. Default %s' % self.nextdenovo_sort)
        nextdenovo_group.add_argument('--nextdenovo-correction_opts', dest="nextdenovo_correction_opts", metavar="nextdenovo_correction_opts", default=self.nextdenovo_correction_opts, help='Correction options. Default %s' % self.nextdenovo_correction_opts)
        nextdenovo_group.add_argument('--nextdenovo-nextgraph_opt', dest="nextdenovo_nextgraph_opt", metavar="nextdenovo_nextgraph_opt", default=self.nextdenovo_nextgraph_opt, help='nextgraph options, see nextgraph for details. Default %s' % self.nextdenovo_nextgraph_opt)

    def register_hypo(self, parser):
        """Register all hypo parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        hypo_group = parser.add_argument_group('Hypo')
        hypo_group.add_argument('--sr-cov', dest="ill_cov", metavar="ill_cov", default= self.ill_cov, type=int, help='Approximate short read coverage for hypo Default %s' % self.ill_cov)
        hypo_group.add_argument('--hypo-proc', dest="hypo_processes", metavar="hypo_processes", default=self.hypo_processes, type = int, help='Number of contigs to be processed in parallel by HyPo. Default %s' % self.hypo_processes)
        hypo_group.add_argument('--hypo-no-lr', dest="hypo_lr", default=self.hypo_lr, action= "store_false", help='Set this to false if you don¡t want to run hypo with long reads. Default %s' % self.hypo_lr)
        hypo_group.add_argument('--hypo-opts', dest="hypo_opts", metavar="hypo_opts", default=self.hypo_opts, help='Additional options to run Hypo. Default %s' % self.hypo_opts)

    def register_racon(self, parser):
        """Register all racon parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        racon_group = parser.add_argument_group('Racon')
        racon_group.add_argument('--racon-dir', dest="racon_dir", metavar="racon_dir", help='Directory with Racon installation. Default %s' % self.racon_dir)
        racon_group.add_argument('--racon-opts', dest="racon_opts", metavar="racon_opts", default=self.racon_opts, help='Extra options to run Racon_wrapper(eg. --split 100000000) to split the assembly in chunks of specified size and decrease memory requirements. Do racon_wrapper -h for more info.')

    def register_medaka(self, parser):
        """Register all medaka parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        medaka_group = parser.add_argument_group('Medaka')
        medaka_group.add_argument('--medaka-env', dest="medaka_env", metavar="medaka_env", help='Conda environment to run Medaka. Default %s' % self.medaka_env)
        medaka_group.add_argument('--medaka-workdir', dest="medaka_workdir", metavar="medaka_workdir", help='Directory to run Medaka. Default %s' % self.medaka_workdir)
        medaka_group.add_argument('--medaka-model', dest="medaka_model", metavar="medaka_model", help='User needs to specify the model that will be used to run Medaka')
        medaka_group.add_argument('--medaka-consensus-opts', dest="medaka_consensus_opts", metavar="medaka_consensus_opts", help='Specify any parameters to change when running medaka consensus')

    def register_pilon(self, parser):
        """Register all pilon parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        pilon_group = parser.add_argument_group('Pilon')
        pilon_group.add_argument('--pilon-path', dest="pilon_path", metavar="pilon_path", help='Path to Pilon executable. Default %s' % self.pilon_path)
        pilon_group.add_argument('--pilon-opts', dest="pilon_opts", metavar="pilon_opts", default=self.pilon_opts, help='Additional options to run Pilon. Default %s' % self.pilon_opts)
        pilon_group.add_argument('--java-opts', dest="java_opts", metavar="java_opts", default=self.java_opts, help='Options for the java execution of Pilon. Default %s' % self.java_opts)
        pilon_group.add_argument('--pilon-subs', dest="pilon_subsampling", metavar="pilon_subsampling", type = int, default=self.pilon_subsampling, help='Percentage of reads to randomly use when running pilon. Default %s' % self.pilon_subsampling)
        pilon_group.add_argument('--pilon-chunks', dest="pilon_chunks", metavar="pilon_chunks", type = int, default=self.pilon_chunks, help='Number of chunks to split the pilon polishing jobs. Default %s' % self.pilon_chunks)

    def register_purgedups(self, parser):
        """Register all purgedups parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        purgedups_group = parser.add_argument_group('Purge_dups')
        purgedups_group.add_argument('--purgedups-cores', type = int, dest="purgedups_cores", metavar="purgedups_cores", default=self.purgedups_cores, help='Number of threads to run purgedups. Default %s' % self.purgedups_cores)
        purgedups_group.add_argument('--purgedups-calcuts-opts', dest="calcuts_opts", metavar="calcuts_opts", default = self.calcuts_opts, help='Adjusted values to run calcuts for purgedups. Default %s' % self.calcuts_opts)

    def register_scaffold10X(self, parser):
        """Register all 10X scaffolding parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        scaffold10X_group = parser.add_argument_group('Scaffold_with_10X')
        scaffold10X_group.add_argument('--tigmint-cores', type = int, dest="tigmint_cores", metavar="tigmint_cores", default=self.tigmint_cores, help='Number of threads to run the 10X scaffolding step. Default %s' % self.tigmint_cores)
        scaffold10X_group.add_argument('--tigmint-opts', dest="tigmint_opts", metavar="tigmint_opts", default = self.tigmint_opts, help='Adjusted values to run the scaffolding with 10X reads.  Default %s' % self.tigmint_opts)

    def register_finalize(self, parser):
        """Register all finalize parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        finalize_group = parser.add_argument_group('Finalize')
        finalize_group.add_argument('--no-final-evals', dest="final_evals", action="store_false", help='If specified, do not run evaluations on final assemblies. Default %s' % self.final_evals)
        finalize_group.add_argument('--busco-lin', dest="busco_lineage", metavar="busco_lineage", help='Path to the lineage directory to run Busco with. Default %s' % self.busco_lineage)
        finalize_group.add_argument('--merqury-db', dest="merqury_db", metavar="merqury_db", help='Meryl database. Default %s' % self.merqury_db)
        finalize_group.add_argument('--meryl-k', dest="meryl_k", metavar="meryl_k", type = int, help='Kmer length to build the meryl database. Default %s' % self.meryl_k)
        finalize_group.add_argument('--meryl-threads', dest="meryl_threads", metavar="meryl_threads", type = int, default = self.meryl_threads, help='Number of threads to run meryl and merqury. Default %s' % self.meryl_threads)

    def register_wildcards(self, parser):
        """Register all wildcards parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        wildcards_group = parser.add_argument_group('Wildcards')
        wildcards_group.add_argument('--ont-list', dest="ONT_wildcards", metavar="ONT_wildcards", help='List with basename of the ONT fastqs that will be used. Default %s' % self.ONT_wildcards)
        wildcards_group.add_argument('--illumina-list', dest="illumina_wildcards", metavar="illumina_wildcards", help='List with basename of the illumina fastqs. Default %s' % self.illumina_wildcards)
        wildcards_group.add_argument('--r10X-list', dest="r10X_wildcards", metavar="r10X_wildcards", help='List with basename of the raw 10X fastqs. Default %s' % self.r10X_wildcards)
####

    def check_parameters(self,args):
        """Check parameters consistency
            
        args -- set of parsed arguments"""

        if args.pipeline_workdir != None:
          args.pipeline_workdir = os.path.abspath(args.pipeline_workdir) + "/"
        else:
          args.pipeline_workdir = os.getcwd() + "/"

        if args.configFile != None:
          args.configFile = os.path.abspath(args.configFile) 
        else:
          args.configFile = os.path.abspath(self.configFile) 

        if args.specFile != None:
          args.specFile = os.path.abspath(args.specFile) 
        else:
          args.specFile = os.path.abspath(self.specFile) 

        if args.ndconfFile != None:
          args.ndconfFile = os.path.abspath(args.ndconfFile) 
        else:
          args.ndconfFile = os.path.abspath(self.ndconfFile) 

        if args.scripts_dir:
          args.scripts_dir = os.path.abspath(args.scripts_dir) + "/"
        else:
          args.scripts_dir = os.path.abspath(self.scripts_dir) + "/"
        if not os.path.exists(args.scripts_dir):
          print (args.scripts_dir + " not found")

        if args.eval_dir:
          args.eval_dir = os.path.abspath(args.eval_dir) + "/"
        else:
          args.eval_dir = args.pipeline_workdir  + self.eval_dir

        if args.base_name == None:
          parser.print_help()
          print ("You need to provide a base name for the project")
          sys.exit(-1)   

        if args.species == None:
          parser.print_help()
          print ("You need to provide the species name")
          sys.exit(-1)   

        if args.stats_out == None:
          args.stats_out = args.eval_dir + args.base_name + ".stats_summary.txt"
        else:
          args.stats_out = os.path.abspath(args.stats_out)

        if args.longranger_path:
          args.longranger_path = os.path.abspath(args.longranger_path)
        else:
          args.longranger_path =  os.path.abspath(self.longranger_path)
        if not os.path.exists(args.longranger_path):
          print (args.longranger_path + " not found")

        if args.genome_size == None:
          parser.print_help()
          print ("You need to provide a genome size estimate")
          sys.exit(-1)
        else:
          if re.search("m", args.genome_size):
            gsize = float(re.split("m", args.genome_size)[0])
          elif re.search("g", args.genome_size):
            gsize = float(re.split("g", args.genome_size)[0]) * 1000
          else:
            parser.print_help()
            print ("Provided genome size estimate is in unknown format, please check suggestions and give and appropriate value")
            sys.exit(-1)
          print ("Genome size is " + str(gsize) + " megabases")

        if args.racon_opts == None:
          if gsize > 1000:
            args.racon_opts = " --split 100000000 "
        elif not re.search("split", args.racon_opts):
            args.racon_opts += " --split 100000000 "

        args.filtlong_qos =  self.filtlong_qos
        args.filtlong_time = self.filtlong_time 
        args.filtlong_queue = self.filtlong_queue 
        args.filtlong_mem =  self.filtlong_mem        

        args.longranger_qos =  self.longranger_qos
        args.longranger_time = self.longranger_time 
        args.longranger_queue = self.longranger_queue   
        args.longranger_mem = self.longranger_mem       

        args.trimgalore_qos = self.trimgalore_qos
        args.trimgalore_time = self.trimgalore_time
        args.trimgalore_queue = self.trimgalore_queue
        args.trimgalore_mem = self.trimgalore_mem

        args.concat_reads_qos =  self.concat_reads_qos
        args.concat_reads_time = self.concat_reads_time 
        args.concat_reads_queue = self.concat_reads_queue
        args.concat_reads_mem =  self.concat_reads_mem 

        args.nanoplot_qos =  self.nanoplot_qos
        args.nanoplot_time = self.nanoplot_time 
        args.nanoplot_queue = self.nanoplot_queue
        args.nanoplot_mem =  self.nanoplot_mem     

        args.kraken2_qos = self.kraken2_qos
        args.kraken2_time = self.kraken2_time 
        args.kraken2_queue = self.kraken2_queue
        args.kraken2_mem = self.kraken2_mem

        args.build_meryl_qos =  self.build_meryl_qos
        args.build_meryl_time = self.build_meryl_time 
        args.build_meryl_queue = self.build_meryl_queue
        args.build_meryl_mem = self.build_meryl_mem

        args.concat_meryl_qos =  self.concat_meryl_qos
        args.concat_meryl_time = self.concat_meryl_time 
        args.concat_meryl_queue = self.concat_meryl_queue
        args.concat_meryl_mem = self.concat_meryl_mem

        args.genomescope_qos =  self.genomescope_qos
        args.genomescope_time = self.genomescope_time 
        args.genomescope_queue = self.genomescope_queue
        args.genomescope_mem = self.genomescope_mem

        args.flye_qos =  self.flye_qos
        args.flye_time = self.flye_time 
        args.flye_queue = self.flye_queue
        args.flye_mem =  self.flye_mem      

        args.nextdenovo_qos =  self.nextdenovo_qos
        args.nextdenovo_time = self.nextdenovo_time 
        args.nextdenovo_queue = self.nextdenovo_queue
        args.nextdenovo_mem =  self.nextdenovo_mem

        args.minimap_qos =  self.minimap_qos
        args.minimap_time = self.minimap_time 
        args.minimap_queue = self.minimap_queue
        args.minimap_mem = self.minimap_mem

        args.bwa_qos =  self.bwa_qos
        args.bwa_time = self.bwa_time 
        args.bwa_queue = self.bwa_queue
        args.bwa_mem = self.bwa_mem

        args.hypo_qos =  self.hypo_qos
        args.hypo_time = self.hypo_time 
        args.hypo_queue = self.hypo_queue
        args.hypo_mem = self.hypo_mem

        args.racon_qos =  self.racon_qos
        args.racon_time = self.racon_time 
        args.racon_queue = self.racon_queue

        args.medaka_qos =  self.medaka_qos
        args.medaka_time = self.medaka_time 
        args.medaka_queue = self.medaka_queue

        args.split_pilon_qos =  self.split_pilon_qos
        args.split_pilon_time = self.split_pilon_time 
        args.split_pilon_queue = self.split_pilon_queue

        args.pilon_qos =  self.pilon_qos
        args.pilon_time = self.pilon_time 
        args.pilon_queue = self.pilon_queue
        args.pilon_threads = self.pilon_threads

        args.join_pilon_qos =  self.join_pilon_qos
        args.join_pilon_time = self.join_pilon_time 
        args.join_pilon_queue = self.join_pilon_queue

        args.nextpolish_lr_qos =  self.nextpolish_lr_qos
        args.nextpolish_lr_time = self.nextpolish_lr_time 
        args.nextpolish_lr_queue = self.nextpolish_lr_queue
        args.nextpolish_lr_mem = self.nextpolish_lr_mem
        
        args.nextpolish_sr_qos =  self.nextpolish_sr_qos
        args.nextpolish_sr_time = self.nextpolish_sr_time 
        args.nextpolish_sr_queue = self.nextpolish_sr_queue
        args.nextpolish_sr_mem = self.nextpolish_sr_mem

        args.purgedups_qos =  self.purgedups_qos
        args.purgedups_time = self.purgedups_time 
        args.purgedups_queue = self.purgedups_queue
        args.purgedups_mem = self.purgedups_mem

        args.tigmint_qos =  self.tigmint_qos
        args.tigmint_time = self.tigmint_time 
        args.tigmint_queue = self.tigmint_queue
        args.tigmint_mem = self.tigmint_mem


        args.stats_qos =  self.stats_qos
        args.stats_time = self.stats_time 
        args.stats_queue = self.stats_queue
        args.stats_mem = self.stats_mem

        args.busco_qos =  self.busco_qos
        args.busco_time = self.busco_time 
        args.busco_queue = self.busco_queue
        args.busco_mem = self.busco_mem

        args.merq_qos =  self.merq_qos
        args.merq_time = self.merq_time 
        args.merq_queue = self.merq_queue
        args.merq_mem = self.merq_mem

        args.fin_qos =  self.fin_qos
        args.fin_time = self.fin_time 
        args.fin_queue = self.fin_queue
        args.fin_mem = self.fin_mem

        if gsize > 500:
          args.java_opts = "Xmx150g"

        if gsize > 1000:
          args.concat_cores = 16
          args.flye_qos = "marathon"
          args.flye_time = "150:00:00"
          args.flye_cores = 128
          args.nextdenovo_qos = "marathon"
          args.nextdenovo_time = "150:00:00"
          args.bwa_time = "24:00:00"
          args.bwa_qos = "long"
          args.minimap_time = "12:00:00"
          args.hypo_time = "10:00:00"
          args.busco_time = "24:00:00"
          args.busco_qos = "long"
          args.merq_time = "6:00:00"
          args.pilon_cores = 20

        if gsize < 100:
          args.busco_time = "2:00:00"
          args.merq_time = "2:00:00"
          args.fin_time = "0:30:00"
          args.nextdenovo_qos = "long"
          args.nextdenovo_time = "24:00:00"
          args.nextdenovo_mem = "100G"
          args.flye_cores = 20

        if args.run_flye == True or args.run_nextdenovo == True or args.racon_rounds > 0 or args.medaka_rounds > 0 or args.nextpolish_ont_rounds > 0 or args.hypo_rounds > 0 or args.run_purgedups != None:
          if args.ONT_filtered:
            args.ONT_filtered = os.path.abspath(args.ONT_filtered)
            args.filtlong_dir = os.path.dirname(args.ONT_filtered) + "/"
          else:
            if args.filtlong_dir == None:
              args.filtlong_dir = args.pipeline_workdir + "s" + args.preprocess_ont_step + "_p01.1_Filtlong/"
            else:
              args.filtlong_dir = os.path.abspath(args.filtlong_dir) + "/" 
            args.ONT_filtered = args.filtlong_dir + "ont.filtlong.fastq.gz"

          if args.ONT_dir == None and args.ONT_reads == None and not os.path.exists(args.ONT_filtered):
            parser.print_help()
            print ("The long reads are needed")
            sys.exit(-1)
          if args.ONT_dir:
            args.ONT_dir = os.path.abspath(args.ONT_dir) + "/" 
            if not os.path.exists(args.ONT_dir):
              parser.print_help()
              print (args.ONT_dir + " not found")
              sys.exit(-1)
            elif args.ONT_wildcards == None:
              args.ONT_wildcards = get_wildcards(args.ONT_dir, args.ONT_wildcards, '.fastq.gz')

          if args.ONT_reads:
            args.ONT_reads = os.path.abspath(args.ONT_reads) 
            if not os.path.exists(args.ONT_reads):
              if args.ONT_dir == None:
                parser.print_help()
                print (args.ONT_reads + " not found")
                sys.exit(-1)

        if args.merqury_db:
          args.merqury_db = os.path.abspath(args.merqury_db)
        
        args.r10X_reads = {}
        if args.pilon_rounds > 0 or args.nextpolish_ill_rounds > 0 or args.hypo_rounds >0  or args.run_tigmint == True or args.merqury_db:
          if args.illumina_dir == None and args.pe1 == None and args.pe2==None and args.r10X==None and args.processed_illumina == None and len(args.raw_10X) == 0 and args.processed_10X == None:
            parser.print_help()
            print ("The illumina reads are needed")
            sys.exit(-1)
          elif args.pe1 or args.pe2:
            if args.pe1:
              args.pe1 = os.path.abspath(args.pe1) 
            else: 
              parser.print_help()
              print ("Both illumina pairs are needed, missing pe1")
              sys.exit(-1)
            if args.pe2:
              args.pe2 = os.path.abspath(args.pe2) 
            else: 
              parser.print_help()
              print ("Both illumina pairs are needed, missing pe2")
              sys.exit(-1)

            if not os.path.exists(args.pe1) or not os.path.exists(args.pe2):
              if args.illumina_dir == None and args.processed_illumina == None:
                parser.print_help()
                print ("Illumina reads are not found and are needed")
                sys.exit(-1)
          elif args.r10X:
            args.r10X = os.path.abspath(args.r10X) 
            if not os.path.exists(args.r10X):
              if len(raw_10X) == 0 and args.processed_10X == None:
                parser.print_help()
                print ("Illumina reads are not found and are needed")
                sys.exit(-1)
          elif args.illumina_dir != None:
            args.illumina_dir = os.path.abspath(args.illumina_dir) + "/" 
            if args.processed_illumina == None:
              args.processed_illumina =  args.pipeline_workdir + "s" + self.preprocess_illumina_step + "_p01.1_preprocess_illumina_reads/trim/"
            if not os.path.exists(args.illumina_dir):
              parser.print_help()
              print (args.illumina_dir + " not found")
              sys.exit(-1)
            elif args.illumina_wildcards == None:
              args.illumina_wildcards = get_wildcards(args.illumina_dir, args.illumina_wildcards, '.1.fastq.gz')
          elif args.processed_illumina != None:
            args.processed_illumina = os.path.abspath(args.processed_illumina) + "/"
            if not os.path.exists(args.processed_illumina):
              parser.print_help()
              print (args.processed_illumina + " not found and missing raw illumina directory")
              sys.exit(-1) 
            elif args.illumina_wildcards == None:
              args.illumina_wildcards = get_wildcards(args.processed_illumina, args.illumina_wildcards, '.1.fastq.gz')
          else:
            if len(args.raw_10X):
             # print (args.raw_10X)

              args.r10X_wildcards = ""

              for my_dict in args.raw_10X:
                for key in my_dict:
              #    print (key)
                  args.r10X_reads[os.path.abspath(key)] = my_dict[key]
                  if args.r10X_wildcards:
                    args.r10X_wildcards += ","
                  args.r10X_wildcards += my_dict[key]
                  if not os.path.exists(os.path.abspath(key)):
                    parser.print_help()
                    print (os.path.abspath(key) + " not found")
                    sys.exit(-1)
                  elif key == None:
                    parser.print_help()
                    print ('If you want to process the 10X reads with longranger, you need to provide the 10X basenames together with the directory')
                    sys.exit(-1)    
              if args.processed_10X == None:
                args.processed_10X = args.pipeline_workdir + "s" + self.preprocess_10X_step + "_p01.1_preprocess_10X_linkedreads"
            if args.processed_10X != None:
              args.processed_10X = os.path.abspath(args.processed_10X) + "/"
              if args.r10X_wildcards == None:
                if not os.path.exists(args.processed_10X):
                  parser.print_help()
                  print (args.processed_10X + " not found")
                  sys.exit(-1)
                else:
                  args.r10X_wildcards = get_wildcards(args.processed_10X, args.r10X_wildcards, '.barcoded.fastq.gz') 
          if args.merqury_db:
            if not os.path.exists(args.merqury_db):
              if args.meryl_k == None:
                parser.print_help()
                print (args.merqury_db + " not found, the pipeline will create it but you need to provide the kmer value")
                sys.exit(-1)             
              elif args.pe1 != None or args.r10X != None or args.r10X_wildcards != None or args.illumina_wildcards != None:
                print (args.merqury_db + " not found, the pipeline will create it")   
              else:
                print ('We cannot create ' + args.merqury_db + ' without illumina reads')

        if args.run_tigmint == True and args.r10X == None:
         # if args.processed_10X != None:   
         #   args.r10X = args.processed_10X + "reads.illumina10X.barcoded.fastq.gz" 
         # else:   
          if args.processed_10X == None:
            parser.print_help()
            print ("10X reads are required to run the 10X scaffolding step")
            sys.exit(-1) 

        if args.run_kraken2 == True:
          if args.kraken2_db == None or args.kraken2_kmers == None:
            parser.print_help()
            print ("Please, specify a kraken2 database and a kmer distribution file if you want to run Kraken2 and Bracken on the reads. Otherwise, do not use the --run-kraken parameter.")
            sys.exit(-1)
          else:
            args.kraken2_db = os.path.abspath(args.kraken2_db)
            args.kraken2_kmers = os.path.abspath(args.kraken2_kmers)
            if not os.path.exists(args.kraken2_db):
              parser.print_help()
              print ("\n" + args.kraken2_db + " should exist and it does not.")
              sys.exit(-1)
            if not os.path.exists(args.kraken2_kmers):
              parser.print_help()
              print ("\n" + args.kraken2_kmers + " should exist and it does not.")
              sys.exit(-1)

        if args.other_flye_opts == None:
          args.other_flye_opts = self.other_flye_opts + " -g " + args.genome_size + " "
        elif not re.match("-g ",args.other_flye_opts) and not re.match("--genome-size ", args.other_flye_opts):
          args.other_flye_opts += " -g " + args.genome_size + " "

        if args.flye_dir == None:
          args.flye_dir = args.pipeline_workdir + "s" + args.flye_step + "_p" + args.preprocess_ont_step + "_flye/"
        else:
          args.flye_dir = os.path.abspath(args.flye_dir) + "/"  
        args.flye_out = args.flye_dir + "flye.assembly.fasta"
             
        if args.nextdenovo_dir == None:
          args.nextdenovo_dir = args.pipeline_workdir + "s" + args.nextdenovo_step + "_p" + args.preprocess_ont_step + "_nextdenovo/"
        else:
          args.nextdenovo_dir = os.path.abspath(args.nextdenovo_dir) + "/" 
        args.nextdenovo_out = args.nextdenovo_dir + "nextdenovo.assembly.fasta"

        if args.racon_dir:
          args.racon_dir = os.path.abspath(args.racon_dir) + "/"
        else:
          args.racon_dir =  os.path.abspath(self.racon_dir) + "/"
        if not os.path.exists(args.racon_dir):
          print (args.racon_dir + " not found")

        if args.medaka_env:
          args.medaka_env = os.path.abspath(args.medaka_env)
        else:
          args.medaka_env =  os.path.abspath(self.medaka_env)
        if not os.path.exists(args.medaka_env):
          print (args.medaka_env + " not found")
         
        if args.medaka_rounds > 0 and not args.medaka_model:
          parser.print_help()
          print ("You have chosen to run Medaka but no model for it has been specified")
          sys.exit(-1)

        if args.medaka_workdir:
          args.medaka_workdir = os.path.abspath(args.medaka_workdir) + "/"
        else:
          args.medaka_workdir = self.medaka_workdir 

        if args.pilon_path:
          args.pilon_path = os.path.abspath(args.pilon_path)
        else:
          args.pilon_path =  os.path.abspath(self.pilon_path)
        # if not os.path.exists(args.pilon_path):
        #   print (args.pilon_path + " not found")

        if args.busco_lineage:
          args.busco_lineage = os.path.abspath(args.busco_lineage)
          if not os.path.exists(args.busco_lineage):
            print (args.busco_lineage + " not found")
        elif args.final_evals == True:
          print ("busco lineage is needed if you want to run Busco")

        args.assemblies={}
        if len(args.assembly_in):
          for my_dict in args.assembly_in:
            for key in my_dict:
              args.assemblies[os.path.abspath(key)] = os.path.abspath(my_dict[key]) + "/"

        args.assemblies_cur = {}
        if len(args.postpolish_assemblies):
          for my_dict in args.postpolish_assemblies:
            for key in my_dict:
              args.assemblies_cur[os.path.abspath(key)] = my_dict[key]

        if args.pilon_rounds > 0 or args.nextpolish_ill_rounds > 0 or args.hypo_rounds >0 or args.racon_rounds > 0 or args.medaka_rounds > 0 or args.nextpolish_ont_rounds:
          if args.run_flye == True:
            if args.polish_flye_dir != None:
              args.polish_flye_dir = os.path.abspath(args.polish_flye_dir)
            else:
              step = float(args.flye_step) + 1
              args.polish_flye_dir = args.pipeline_workdir + "s0" + str(step) + "_p" + args.flye_step + "_polishing/"
            args.assemblies[args.flye_out] = args.polish_flye_dir
          if args.run_nextdenovo== True:
            if args.polish_nextdenovo_dir != None:
              args.polish_nextdenovo_dir = os.path.abspath(args.polish_nextdenovo_dir)
            else:
              step = float(args.nextdenovo_step) + 1
              args.polish_nextdenovo_dir = args.pipeline_workdir + "s0" + str(step) + "_p" + args.nextdenovo_step + "_polishing/"
            args.assemblies[args.nextdenovo_out] = args.polish_nextdenovo_dir

          if args.run_purgedups == True:          
            if len(args.assemblies) > 0:
              pol_bases = {}
              base_tmp = ""
              if  args.hypo_rounds >0:
                pol_bases["hypo"] = "hypo" + str(args.hypo_rounds)
              if args.racon_rounds > 0:
                base_tmp+= "racon" + str(args.racon_rounds)
              if base_tmp != "":
                base_tmp += "."
              if args.medaka_rounds > 0:
                base_tmp += "medaka" + str(args.medaka_rounds)
              if args.pilon_rounds > 0:
                if base_tmp != "":
                  base_tmp += "."
                base_tmp += "pilon" + str(args.pilon_rounds)
              if base_tmp != "":
                pol_bases["rmp"] = base_tmp
                base_tmp = ""
              if args.nextpolish_ont_rounds > 0:
                base_tmp+= "nextpolish_ont" + str(args.nextpolish_ont_rounds)
              if args.nextpolish_ill_rounds > 0:
                if base_tmp != "":
                  base_tmp += "."
                base_tmp += "nextpolish_ill" + str(args.nextpolish_ill_rounds)
              if base_tmp != "":
                pol_bases["nextpolish"] = base_tmp
              paths = 0
              for m in args.assemblies:
                bpol = os.path.splitext(os.path.basename(m))[0]
                path = args.assemblies[m]
                pstep = path.split('/')[-2].split('_')[0]
                
                nstep = pstep.replace('s','')
                #cstep = float(nstep) + 1 + paths
                for p in pol_bases:
                  cstep = round(float(nstep) + 1 + paths,1)
                  args.assemblies_cur[args.assemblies[m] + p + "/" + bpol + "." +  pol_bases[p] + ".fasta"] = "s0" + str(cstep) + "_p" + nstep
                  paths += 0.1

  
###

    def storeGeneralParameters(self,args):
        """Updates general parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.generalParameters["configFile"] = args.configFile
        self.generalParameters["specFile"] = args.specFile
        self.generalParameters["ndconfFile"] = args.ndconfFile
        self.generalParameters["lr_type"] = args.lr_type
        self.generalParameters["concat_cores"] = args.concat_cores
        self.generalParameters["keep_intermediate"] = args.keepintermediate
        self.generalParameters["base_name"] = args.base_name
        self.generalParameters["species"] = args.species
        self.generalParameters["genome_size"] = args.genome_size
        self.generalParameters["preprocess_ont_step"] = args.preprocess_ont_step
        self.generalParameters["preprocess_illumina_step"] = args.preprocess_illumina_step
        self.generalParameters["preprocess_10X_step"] = args.preprocess_10X_step
        self.generalParameters["flye_step"] = args.flye_step
        self.generalParameters["run_flye"] = args.run_flye
        self.generalParameters["nextdenovo_step"] = args.nextdenovo_step
        self.generalParameters["run_nextdenovo"] = args.run_nextdenovo
        self.generalParameters["racon_rounds"] = args.racon_rounds
        self.generalParameters["pilon_rounds"] = args.pilon_rounds
        self.generalParameters["medaka_rounds"] = args.medaka_rounds
        self.generalParameters["nextpolish_ont_rounds"] = args.nextpolish_ont_rounds
        self.generalParameters["nextpolish_ill_rounds"] = args.nextpolish_ill_rounds
        self.generalParameters["hypo_rounds"] = args.hypo_rounds
        self.generalParameters["racon_cores"] = args.racon_cores
        self.generalParameters["nextpolish_cores"] = args.nextpolish_cores
        self.generalParameters["minimap2_cores"] = args.minimap2_cores
        self.generalParameters["BWA_cores"] = args.bwa_cores
        self.generalParameters["pilon_cores"] = args.pilon_cores
        self.generalParameters["medaka_cores"] = args.medaka_cores
        self.generalParameters["hypo_cores"] = args.hypo_cores
        self.generalParameters["busco_cores"] = args.busco_cores
        self.generalParameters["longranger_cores"] = args.longranger_cores
        self.generalParameters["longranger_path"] = args.longranger_path
        self.generalParameters["genomescope_additional_options"] = args.genomescope_additional
        self.generalParameters["ploidy"] = args.ploidy
        self.generalParameters["run_purgedups"] = args.run_purgedups
        self.generalParameters["run_tigmint"] = args.run_tigmint
        self.generalParameters["run_kraken2"] = args.run_kraken2
        self.allParameters["Parameters"] = self.generalParameters

    def storeallSpecParameters(self,args):
        """Updates rule all cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.allSpecParameters["name"] = "{rule}_{base}_assembly_pipeline"
        self.allSpecParameters["qos"] = self.all_qos
        self.allSpecParameters["time"] = self.all_time
        self.allSpecParameters["queue"] = self.all_queue
        self.allParameters ["all"] = self.allSpecParameters


    def storeInputParameters(self,args):
        """Updates input parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.inputParameters["scripts_dir"] = args.scripts_dir
        self.inputParameters["ONT_reads"] = args.ONT_reads
        self.inputParameters["ONT_dir"] = args.ONT_dir
        self.inputParameters["ONT_filtered"] = args.ONT_filtered
        self.inputParameters["ILLUMINA_pe1"] = args.pe1
        self.inputParameters["ILLUMINA_pe2"] = args.pe2
        self.inputParameters["processed_illumina"] = args.processed_illumina
        self.inputParameters["raw_10X"] = args.r10X_reads
        self.inputParameters["processed_10X"] = args.processed_10X
        self.inputParameters["ILLUMINA_10X"] = args.r10X
        self.inputParameters["illumina_dir"] = args.illumina_dir
        self.inputParameters["Assemblies for polishing"] = args.assemblies
        self.inputParameters["Assemblies for postpolishing"] = args.assemblies_cur
        self.allParameters ["Inputs"] = self.inputParameters

    def storeOutputParameters(self,args):
        """Updates output parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.outputParameters["base_dir"] = args.pipeline_workdir 
        self.outputParameters["filtlong_dir"] = args.filtlong_dir
        self.outputParameters["flye_dir"] = args.flye_dir
        self.outputParameters["nextdenovo_dir"] = args.nextdenovo_dir
        self.outputParameters["flye_out"] = args.flye_out
        self.outputParameters["nextdenovo_out"] = args.nextdenovo_out
        self.outputParameters["eval_dir"] = args.eval_dir
        self.outputParameters["stats_out"] = args.stats_out
        self.allParameters ["Outputs"] = self.outputParameters

    def storelongrangerSpecParameters(self,args):
        """Updates longranger cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.longrangerSpecParameters["name"] = "{rule}_" + args.base_name + "_{base}_{wildcards.bname}"
        self.longrangerSpecParameters["qos"] = args.longranger_qos
        self.longrangerSpecParameters["time"] = args.longranger_time
        self.longrangerSpecParameters["queue"] = args.longranger_queue
        self.longrangerSpecParameters["mem"] = args.longranger_mem
        self.allParameters ["long_ranger"] = self.longrangerSpecParameters

    def storeTrimgaloreParameters(self,args):
        """Updates the Trim_Galore parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.trimgaloreParameters["options"] = args.trim_galore_opts
        self.trimgaloreParameters["Trim_Illumina_cores"] = args.Trim_Illumina_cores
        self.allParameters ["Trim_Galore"] = self.trimgaloreParameters

    def storetrimgaloreSpecParameters(self,args):
        """Updates trimgalore cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.trimgaloreSpecParameters["name"] = "{rule}_{base}_{wildcards.file}"
        self.trimgaloreSpecParameters["qos"] = args.trimgalore_qos
        self.trimgaloreSpecParameters["time"] = args.trimgalore_time
        self.trimgaloreSpecParameters["queue"] = args.trimgalore_queue
        self.trimgaloreSpecParameters["mem"] = args.trimgalore_mem
        self.allParameters ["trim_galore"] = self.trimgaloreSpecParameters

    def storeconcatreadsSpecParameters(self,args):
        """Updates concat reads cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.concatreadsSpecParameters["name"] = "{rule}_{base}_{wildcards.ext}"
        self.concatreadsSpecParameters["qos"] = args.concat_reads_qos
        self.concatreadsSpecParameters["time"] = args.concat_reads_time
        self.concatreadsSpecParameters["queue"] = args.concat_reads_queue
        self.concatreadsSpecParameters["mem"] = args.concat_reads_mem
        self.allParameters ["concat_reads"] = self.concatreadsSpecParameters

    def storenanoplotSpecParameters(self,args):
        """Updates Nanoplot cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.nanoplotSpecParameters["name"] = "{rule}_{base}_{wildcards.prefix}"
        self.nanoplotSpecParameters["qos"] = args.nanoplot_qos
        self.nanoplotSpecParameters["time"] = args.nanoplot_time
        self.nanoplotSpecParameters["queue"] = args.nanoplot_queue
        self.nanoplotSpecParameters["mem"] = args.nanoplot_mem
        self.allParameters ["nanoplot"] = self.nanoplotSpecParameters

    def storekraken2Parameters(self,args):
        """Updates Kraken2 parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.kraken2Parameters["database"] = args.kraken2_db
        self.kraken2Parameters["kmer_dist"] = args.kraken2_kmers
        self.kraken2Parameters["threads"] = args.kraken2_threads
        self.kraken2Parameters["additional_opts"] = args.additional_kraken2_opts
        self.allParameters ["Kraken2"] = self.kraken2Parameters

    def storekraken2SpecParameters(self,args):
        """Updates Kraken2 cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.kraken2SpecParameters["name"] = "{rule}_{base}_{params.prefix}"
        self.kraken2SpecParameters["qos"] = args.kraken2_qos
        self.kraken2SpecParameters["time"] = args.kraken2_time
        self.kraken2SpecParameters["queue"] = args.kraken2_queue
        self.kraken2SpecParameters["mem"] = args.kraken2_mem
        self.allParameters ["Kraken2"] = self.kraken2SpecParameters

    def storebuildmerylSpecParameters(self,args):
        """Updates build meryl cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.buildmerylSpecParameters["name"] = "{rule}_" + args.base_name + "_{base}_{wildcards.db}"
        self.buildmerylSpecParameters["qos"] = args.build_meryl_qos
        self.buildmerylSpecParameters["time"] = args.build_meryl_time
        self.buildmerylSpecParameters["queue"] = args.build_meryl_queue
        self.buildmerylSpecParameters["mem"] = args.build_meryl_mem
        self.allParameters ["build_meryl_db"] = self.buildmerylSpecParameters

    def storeconcatmerylSpecParameters(self,args):
        """Updates concat meryl cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.concatmerylSpecParameters["name"] = "{rule}_" + args.base_name + "_{base}"
        self.concatmerylSpecParameters["qos"] = args.concat_meryl_qos
        self.concatmerylSpecParameters["time"] = args.concat_meryl_time
        self.concatmerylSpecParameters["queue"] = args.concat_meryl_queue
        self.concatmerylSpecParameters["mem"] = args.concat_meryl_mem
        self.allParameters ["concat_meryl"] = self.concatmerylSpecParameters

    def storegenomescopeSpecParameters(self,args):
        """Updates genomescope cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.genomescopeSpecParameters["name"] = "{rule}_" + args.base_name + "_{base}"
        self.genomescopeSpecParameters["qos"] = args.genomescope_qos
        self.genomescopeSpecParameters["time"] = args.genomescope_time
        self.genomescopeSpecParameters["queue"] = args.genomescope_queue
        self.genomescopeSpecParameters["mem"] = args.genomescope_mem
        self.allParameters ["genomescope2"] = self.genomescopeSpecParameters


    def storeFiltlongParameters(self,args):
        """Updates filtlong parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.filtlongParameters["Filtlong minlen"] = args.filtlong_minlen
        self.filtlongParameters["Filtlong min_mean_q"] = args.filtlong_min_mean_q
        self.filtlongParameters["options"] = args.filtlong_opts
        self.allParameters ["Filtlong"] = self.filtlongParameters

    def storefiltlongSpecParameters(self,args):
        """Updates filtlong cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.filtlongSpecParameters["name"] = "{rule}_{base}_s" + args.preprocess_ont_step 
        self.filtlongSpecParameters["qos"] = args.filtlong_qos
        self.filtlongSpecParameters["time"] = args.filtlong_time
        self.filtlongSpecParameters["mem"] = args.filtlong_mem
        self.filtlongSpecParameters["queue"] = args.filtlong_queue
        self.allParameters ["filtlong"] = self.filtlongSpecParameters

    def storeFlyeParameters(self,args):
        """Updates flye parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.flyeParameters["Flye cores"] = args.flye_cores
        self.flyeParameters["Flye polishing iterations"] = args.flye_pol_it
        self.flyeParameters["options"] = args.other_flye_opts
        self.allParameters ["Flye"] = self.flyeParameters

    def storeflyeSpecParameters(self,args):
        """Updates flye cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.flyeSpecParameters["name"] = "{rule}_{base}_s" + args.flye_step 
        self.flyeSpecParameters["qos"] = args.flye_qos
        self.flyeSpecParameters["time"] = args.flye_time
        self.flyeSpecParameters["queue"] = args.flye_queue
        self.flyeSpecParameters["mem"] = args.flye_mem
        self.allParameters ["flye"] = self.flyeSpecParameters

    def storeNextdenovoParameters(self,args):
        """Updates nextdenovo parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.nextdenovoParameters["Nextdenovo cores"] = args.nextdenovo_cores
        self.allParameters ["Nextdenovo"] = self.nextdenovoParameters

    def storenextdenovoSpecParameters(self,args):
        """Updates nextdenovo cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.nextdenovoSpecParameters["name"] = "{rule}_{base}_s" + args.nextdenovo_step 
        self.nextdenovoSpecParameters["qos"] = args.nextdenovo_qos
        self.nextdenovoSpecParameters["time"] = args.nextdenovo_time
        self.nextdenovoSpecParameters["queue"] = args.nextdenovo_queue
        self.nextdenovoSpecParameters["mem"] = args.nextdenovo_mem
        self.allParameters ["nextdenovo"] = self.nextdenovoSpecParameters

    def storeminimapSpecParameters(self,args):
        """Updates minimap2 cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.minimapSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.name}_{wildcards.ext}" 
        self.minimapSpecParameters["qos"] = args.minimap_qos
        self.minimapSpecParameters["time"] = args.minimap_time
        self.minimapSpecParameters["queue"] = args.minimap_queue
        self.minimapSpecParameters["mem"] = args.minimap_mem
        self.allParameters ["align_ont"] = self.minimapSpecParameters

    def storebwaSpecParameters(self,args):
        """Updates BWA cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.bwaSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.name}" 
        self.bwaSpecParameters["qos"] = args.bwa_qos
        self.bwaSpecParameters["time"] = args.bwa_time
        self.bwaSpecParameters["queue"] = args.bwa_queue
        self.bwaSpecParameters["mem"] = args.bwa_mem
        self.allParameters ["align_illumina"] = self.bwaSpecParameters

    def storeHypoParameters(self,args):
        """Updates hypo parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.hypoParameters["illumina coverage"] = args.ill_cov
        self.hypoParameters["processes"] = args.hypo_processes
        self.hypoParameters["long_reads"] = args.hypo_lr
        self.hypoParameters["options"] = args.hypo_opts
        self.allParameters ["Hypo"] = self.hypoParameters

    def storehypoSpecParameters(self,args):
        """Updates hypo cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.hypoSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.base}.{wildcards.param}"
        self.hypoSpecParameters["qos"] = args.hypo_qos
        self.hypoSpecParameters["time"] = args.hypo_time
        self.hypoSpecParameters["queue"] = args.hypo_queue
        self.hypoSpecParameters["mem"] = args.hypo_mem
        self.allParameters ["hypo"] = self.hypoSpecParameters

    def storeRaconParameters(self,args):
        """Updates racon parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.raconParameters["Racon dir"] = args.racon_dir
        self.raconParameters["options"] = args.racon_opts
        self.allParameters ["Racon"] = self.raconParameters

    def storeraconSpecParameters(self,args):
        """Updates racon cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.raconSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.base}.{wildcards.param}"
        self.raconSpecParameters["qos"] = args.racon_qos
        self.raconSpecParameters["time"] = args.racon_time
        self.raconSpecParameters["queue"] = args.racon_queue
        self.allParameters ["racon"] = self.raconSpecParameters

    def storeMedakaParameters(self,args):
        """Updates medaka parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.medakaParameters["environment"] = args.medaka_env
        self.medakaParameters["Model"] = args.medaka_model
        self.medakaParameters["Working_dir"] = args.medaka_workdir
        self.medakaParameters["consensus options"] = args.medaka_consensus_opts
        self.allParameters ["Medaka"] = self.medakaParameters

    def storemedakaSpecParameters(self,args):
        """Updates medaka cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.medakaSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.base}.{wildcards.param}"
        self.medakaSpecParameters["qos"] = args.medaka_qos
        self.medakaSpecParameters["time"] = args.medaka_time
        self.medakaSpecParameters["queue"] = args.medaka_queue
        self.allParameters ["medaka"] = self.medakaSpecParameters

    def storePilonParameters(self,args):
        """Updates pilon parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.pilonParameters["path"] = args.pilon_path
        self.pilonParameters["options"] = args.pilon_opts
        self.pilonParameters["JAVA options"] = args.java_opts
        self.pilonParameters["subsampling"] = args.pilon_subsampling
        self.pilonParameters["chunks"] = args.pilon_chunks
        self.allParameters ["Pilon"] = self.pilonParameters

    def storesplitpilonSpecParameters(self,args):
        """Updates split pilon cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.splitpilonSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.base}.{wildcards.param}"
        self.splitpilonSpecParameters["qos"] = args.split_pilon_qos
        self.splitpilonSpecParameters["time"] = args.split_pilon_time
        self.splitpilonSpecParameters["queue"] = args.split_pilon_queue
        self.allParameters ["split_pilon"] = self.splitpilonSpecParameters

    def storepilonSpecParameters(self,args):
        """Updates pilon cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.pilonSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.base}.{wildcards.param}"
        self.pilonSpecParameters["qos"] = args.pilon_qos
        self.pilonSpecParameters["time"] = args.pilon_time
        self.pilonSpecParameters["queue"] = args.pilon_queue
        self.pilonSpecParameters["threads"] = args.pilon_threads
        self.allParameters ["pilon"] = self.pilonSpecParameters

    def storejoinpilonSpecParameters(self,args):
        """Updates join pilon cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.joinpilonSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.base}.{wildcards.param}"
        self.joinpilonSpecParameters["qos"] = args.join_pilon_qos
        self.joinpilonSpecParameters["time"] = args.join_pilon_time
        self.joinpilonSpecParameters["queue"] = args.join_pilon_queue
        self.allParameters ["join_pilon"] = self.joinpilonSpecParameters

    def storenextpolishlrSpecParameters(self,args):
        """Updates nextpolish lr cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.nextpolishlrSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.base}.{wildcards.param}"
        self.nextpolishlrSpecParameters["qos"] = args.nextpolish_lr_qos
        self.nextpolishlrSpecParameters["time"] = args.nextpolish_lr_time
        self.nextpolishlrSpecParameters["queue"] = args.nextpolish_lr_queue
        self.nextpolishlrSpecParameters["mem"] = args.nextpolish_lr_mem
        self.allParameters ["nextpolish_lr"] = self.nextpolishlrSpecParameters

    def storenextpolishsrSpecParameters(self,args):
        """Updates nextpolish sr cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.nextpolishsrSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.base}.{wildcards.param}"
        self.nextpolishsrSpecParameters["qos"] = args.nextpolish_sr_qos
        self.nextpolishsrSpecParameters["time"] = args.nextpolish_sr_time
        self.nextpolishsrSpecParameters["queue"] = args.nextpolish_sr_queue
        self.nextpolishsrSpecParameters["mem"] = args.nextpolish_sr_mem
        self.allParameters ["nextpolish_sr"] = self.nextpolishsrSpecParameters

    def storePurgedupsParameters(self,args):
        """Updates purge_dups parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.purgedupsParameters["purgedups_cores"] = args.purgedups_cores
        self.purgedupsParameters["calcuts_options"] = args.calcuts_opts
        self.allParameters ["Purge_dups"] = self.purgedupsParameters

    def storepurgedupsSpecParameters(self,args):
        """Updates purgedups cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.purgedupsSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.base_in}"
        self.purgedupsSpecParameters["qos"] = args.purgedups_qos
        self.purgedupsSpecParameters["time"] = args.purgedups_time
        self.purgedupsSpecParameters["queue"] = args.purgedups_queue
        self.purgedupsSpecParameters["mem"] = args.purgedups_mem
        self.allParameters ["purge_dups"] = self.purgedupsSpecParameters


    def storescaffold10XParameters(self,args):
        """Updates 10X scaffolding parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.scaffold10XParameters["tigmint_cores"] = args.tigmint_cores
        self.scaffold10XParameters["tigmint_options"] = args.tigmint_opts
        self.allParameters ["scaffolding_10X"] = self.scaffold10XParameters

    def storescaffold10XSpecParameters(self,args):
        """Updates 10X scaffolding cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.scaffold10XSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.base_in}"
        self.scaffold10XSpecParameters["qos"] = args.tigmint_qos
        self.scaffold10XSpecParameters["time"] = args.tigmint_time
        self.scaffold10XSpecParameters["queue"] = args.tigmint_queue
        self.scaffold10XSpecParameters["mem"] = args.tigmint_mem
        self.allParameters ["scaffolding_10X"] = self.scaffold10XSpecParameters

    def storeFinalizeParameters(self,args):
        """Updates finalize parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.finalizeParameters["final Evaluations"] = args.final_evals
        self.finalizeParameters["BUSCO lineage"] = args.busco_lineage
        self.finalizeParameters["Merqury db"] = args.merqury_db
        self.finalizeParameters["Meryl K"] = args.meryl_k
        self.finalizeParameters["Meryl threads"] = args.meryl_threads
        self.allParameters ["Finalize"] = self.finalizeParameters

    def storestatsSpecParameters(self,args):
        """Updates stats cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.statsSpecParameters["name"] = "{rule}_{base}.{params.outbase}"
        self.statsSpecParameters["qos"] = args.stats_qos
        self.statsSpecParameters["time"] = args.stats_time
        self.statsSpecParameters["queue"] = args.stats_queue
        self.statsSpecParameters["mem"] = args.stats_mem
        self.allParameters ["get_stats"] = self.statsSpecParameters

    def storebuscoSpecParameters(self,args):
        """Updates busco cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.buscoSpecParameters["name"] = "{rule}_{base}.{params.buscobase}"
        self.buscoSpecParameters["qos"] = args.busco_qos
        self.buscoSpecParameters["time"] = args.busco_time
        self.buscoSpecParameters["queue"] = args.busco_queue
        self.buscoSpecParameters["mem"] = args.busco_mem
        self.allParameters ["run_busco"] = self.buscoSpecParameters

    def storemerqurySpecParameters(self,args):
        """Updates merqury cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.merqSpecParameters["name"] = "{rule}_{base}.{wildcards.merqbase}"
        self.merqSpecParameters["qos"] = args.merq_qos
        self.merqSpecParameters["time"] = args.merq_time
        self.merqSpecParameters["queue"] = args.merq_queue
        self.merqSpecParameters["mem"] = args.merq_mem
        self.allParameters ["run_merqury"] = self.merqSpecParameters

    def storefinalizeSpecParameters(self,args):
        """Updates finalize cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.finSpecParameters["name"] = "{rule}_{base}"
        self.finSpecParameters["qos"] = args.fin_qos
        self.finSpecParameters["time"] = args.fin_time
        self.finSpecParameters["queue"] = args.fin_queue
        self.finSpecParameters["mem"] = args.fin_mem
        self.allParameters ["finalize"] = self.finSpecParameters

    def storeget_reportSpecParameters(self,args):
        """Updates finalize cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.repSpecParameters["name"] = "{rule}_{base}"
        self.repSpecParameters["qos"] = args.fin_qos
        self.repSpecParameters["time"] = args.fin_time
        self.repSpecParameters["queue"] = args.fin_queue
        self.allParameters ["get_report"] = self.repSpecParameters
       
    def storeWildcardParameters(self,args):
        """Updates wildcard parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.wildcardParameters["ONT_wildcards"] = args.ONT_wildcards
        self.wildcardParameters["illumina_wildcards"] = args.illumina_wildcards
        self.wildcardParameters["10X_wildcards"] = args.r10X_wildcards
        self.allParameters ["Wildcards"] = self.wildcardParameters

    def create_nextdenovo_config(self, args):
        nextdenovo_config.add_section('General')
        nextdenovo_config.set('General', 'job_type', 'local')
        nextdenovo_config.set('General', 'task', args.nextdenovo_task)
        nextdenovo_config.set('General', 'rewrite', args.nextdenovo_rewrite)
        nextdenovo_config.set('General', 'parallel_jobs', str(args.nextdenovo_parallel_jobs))

        if args.keepintermediate == True:
          nextdenovo_config.set('General', 'deltmp', 'no')
        else:
          nextdenovo_config.set('General', 'deltmp', 'yes')

        if re.search("nano", args.lr_type):
          nextdenovo_config.set('General', 'read_type', 'ont')
          if args.lr_type == "nano-raw":
            nextdenovo_config.set('General', 'input_type', 'raw')
          else:
            nextdenovo_config.set('General', 'input_type', 'corrected')
        elif re.search("hifi", args.lr_type):
          nextdenovo_config.set('General', 'read_type', 'hifi')
        else:
          print ("Need to select proper read type for running Nextdenovo")

        nextdenovo_config.set('General',  'workdir', args.nextdenovo_dir)
        nextdenovo_config.set('General',  'input_fofn', args.nextdenovo_dir + 'long_reads.fofn')

        nextdenovo_config.add_section('correct_option')
        nextdenovo_config.set('correct_option',  'read_cutoff', args.nextdenovo_minreadlen)
        nextdenovo_config.set('correct_option',  'genome_size', args.genome_size)
        nextdenovo_config.set('correct_option',  'seed_depth', str(args.nextdenovo_seeddepth))
        nextdenovo_config.set('correct_option',  'seed_cutoff', str(args.nextdenovo_seedcutoff))
        nextdenovo_config.set('correct_option',  'blocksize', args.nextdenovo_blocksize)
        nextdenovo_config.set('correct_option',  'pa_correction', str(args.nextdenovo_pa_correction))
        nextdenovo_config.set('correct_option',  'minimap2_options_raw', args.nextdenovo_minimap_raw)
        nextdenovo_config.set('correct_option',  'sort_options', args.nextdenovo_sort)
        nextdenovo_config.set('correct_option',  'correction_options', args.nextdenovo_correction_opts)

        nextdenovo_config.add_section('assemble_option')
        nextdenovo_config.set('assemble_option',  'minimap2_options_cns', args.nextdenovo_minimap_cns)
        nextdenovo_config.set('assemble_option',  'minimap2_options_map', args.nextdenovo_minimap_map)
        nextdenovo_config.set('assemble_option',  'nextgraph_options', args.nextdenovo_nextgraph_opt)

#1.Create object class Configuration File
configManager = CreateConfigurationFile()
specManager = CreateConfigurationFile()
NDConfManager = CreateConfigurationFile()

#2.Create object for argument parsinng
parser = argparse.ArgumentParser(prog="create_configuration_file",
                description="Create a configuration json file for the assembly pipeline."
                )     
#2.1 Updates arguments and parsing
configManager.register_parameter(parser)

args = parser.parse_args()

#2.2 Check Parameters
configManager.check_parameters(args)

#3. store arguments to super map structure
configManager.storeGeneralParameters(args)
configManager.storeInputParameters(args)
configManager.storeOutputParameters(args)
configManager.storeTrimgaloreParameters(args)
configManager.storekraken2Parameters(args)
configManager.storeFiltlongParameters(args)
configManager.storeFlyeParameters(args)
configManager.storeNextdenovoParameters(args)
configManager.storeHypoParameters(args)
configManager.storeRaconParameters(args)
configManager.storeMedakaParameters(args)
configManager.storePilonParameters(args)
configManager.storePurgedupsParameters(args)
configManager.storescaffold10XParameters(args)
configManager.storeFinalizeParameters(args)
configManager.storeWildcardParameters(args)

specManager.storeallSpecParameters(args)
if args.r10X_wildcards != None:
  specManager.storelongrangerSpecParameters(args)
if args.illumina_dir != None:
  specManager.storetrimgaloreSpecParameters(args)
if args.illumina_wildcards != None or args.ONT_wildcards != None or args.r10X_wildcards:
  specManager.storeconcatreadsSpecParameters(args)
if args.ONT_wildcards != None or args.ONT_reads != None:
  if not os.path.exists(args.ONT_filtered):
    specManager.storefiltlongSpecParameters(args)
  specManager.storenanoplotSpecParameters(args)
if args.run_kraken2 == True:
  specManager.storekraken2SpecParameters(args)
if args.run_flye == True:
  specManager.storeflyeSpecParameters(args)
if args.run_nextdenovo == True:
  specManager.storenextdenovoSpecParameters(args)
  with open(args.ndconfFile, 'w') as ndconf:
    configManager.create_nextdenovo_config(args)
    nextdenovo_config.write(ndconf)
if args.pilon_rounds > 0 or args.nextpolish_ill_rounds > 0 or args.hypo_rounds >0:
  specManager.storebwaSpecParameters(args)
if args.racon_rounds > 0 or args.medaka_rounds > 0 or args.nextpolish_ont_rounds > 0 or args.hypo_rounds >0 or args.run_purgedups == True:
  specManager.storeminimapSpecParameters(args)
if args.hypo_rounds > 0:
  specManager.storehypoSpecParameters(args)
if args.racon_rounds > 0:
  specManager.storeraconSpecParameters(args)
if args.medaka_rounds > 0:
  specManager.storemedakaSpecParameters(args)
if args.pilon_rounds > 0:
  specManager.storesplitpilonSpecParameters(args)
  specManager.storepilonSpecParameters(args)
  specManager.storejoinpilonSpecParameters(args)
if args.nextpolish_ont_rounds > 0:
  specManager.storenextpolishlrSpecParameters(args)
if args.nextpolish_ill_rounds > 0:
  specManager.storenextpolishsrSpecParameters(args)
if args.run_purgedups == True:
  specManager.storepurgedupsSpecParameters(args)
if args.run_tigmint == True:
  specManager.storescaffold10XSpecParameters(args)
if args.merqury_db:
  if not os.path.exists(args.merqury_db):
    specManager.storebuildmerylSpecParameters(args)
    specManager.storeconcatmerylSpecParameters(args)
    specManager.storegenomescopeSpecParameters(args)
  specManager.storemerqurySpecParameters(args)
if args.final_evals:
  specManager.storestatsSpecParameters(args)
  specManager.storebuscoSpecParameters(args)
specManager.storefinalizeSpecParameters(args)
specManager.storeget_reportSpecParameters(args)
#4. Store JSON file
with open(args.configFile, 'w') as of:
    json.dump(configManager.allParameters, of, indent=2)
with open(args.specFile, 'w') as of:
    json.dump(specManager.allParameters, of, indent=2)
