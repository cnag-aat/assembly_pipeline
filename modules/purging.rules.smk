from datetime import datetime
import re
import os

date = datetime.now().strftime('%Y%m%d.%H%M%S')

if not os.path.exists("logs"):
  os.makedirs("logs")

rule purge_dups:
  input:
    assembly_in = "assembly.fa",
    reads = os.getcwd() + "/ontreads.fastq.gz",
    mapping = os.getcwd() + "/mappings/ontreads.paf.gz",
  output:
    assembly_out = "assembly.purged.fa"
  params:
    module = "PURGEDUPS/1.2.5",
    base = "assembly",
    dir = os.getcwd(),
    calcuts_opts = ""
  threads: 12
  log:
    "logs/" + str(date) + ".purge_dups.out",
    "logs/" + str(date) + ".purge_dups.err",
  run:
    if not os.path.exists(params.dir + "/" +  params.base + ".split.self.paf.gz"):
      shell(
        "module purge; module load {params.module} PIGZ/2.3.3;"
        "cd {params.dir};"
        "split_fa {input.assembly_in} > {params.base}.split;"
        "minimap2 -t {threads} -xasm5 -DP {params.base}.split {params.base}.split | pigz -p {threads} -c > {params.base}.split.self.paf.gz;"
        "pbcstat {input.mapping};"
      )
    else:
      shell(
        "echo '{params.dir}{params.base}.split.self.paf.gz already exists, skipping alignment';"
      )
    shell(
        "module purge; module load {params.module} PIGZ/2.3.3;"
        "cd {params.dir};"
        "calcuts {params.calcuts_opts} PB.stat > cutoffs 2>calcults.log;"
    )
    cov = 0
    for line in open(params.dir + "/calcults.log"):  
      if re.search("M::calcuts] mean: ", line):
        line=line.replace("M::calcuts] mean: ", '')
      #  print(line)
        cov = int(line.split(',')[0].lstrip('['))
    print (cov)

    if params.calcuts_opts != "" or cov == 0:
      shell(
        "module purge; module load {params.module} PIGZ/2.3.3;"
        "cd {params.dir};"
        "purge_dups -2 -T cutoffs -c PB.base.cov {params.base}.split.self.paf.gz > dups.bed 2> purge_dups.log;"
        "get_seqs -e dups.bed {input.assembly_in};"
        "ln -s purged.fa {output.assembly_out};"
        "module purge; module load gcc/6.3.0 PYTHON/3.7.1;"
        "python3 /apps/{params.module}/scripts/hist_plot.py -c cutoffs PB.stat PB.cov.png;"
      )
    else:
      shell(
        "module purge; module load {params.module} PIGZ/2.3.3;"
        "cd {params.dir};"
        "echo 'No manual cutoffs have been given, we will run purgedups first to estimate the coverage and then with the adjusted parameters';"
        #"calcuts PB.stat > cutoffs 2>calcults.log;"
        "purge_dups -2 -T cutoffs -c PB.base.cov {params.base}.split.self.paf.gz > dups.bed 2> purge_dups.log;"
        "cov=$(cat calcults.log | sed -n 's/\[M::calcuts\] mean: //p' | cut -f 1 -d ',');"
        "mean=$(($cov*75/100));"
        "min=$(cut -f 1 cutoffs);"
        "max=$(cut -f 6 cutoffs);"
        "calcuts -l $min -m $mean -u $max PB.stat > cutoffs_adjusted;"
        "purge_dups -2 -T cutoffs_adjusted -c PB.base.cov {params.base}.split.self.paf.gz > adjusted_dups.bed 2> adjusted_purge_dups.log;"
        "get_seqs -e adjusted_dups.bed {input.assembly_in};"
        "ln -s purged.fa {output.assembly_out};"
        "module purge; module load gcc/6.3.0 PYTHON/3.7.1;"
        "python3 /apps/{params.module}/scripts/hist_plot.py -c cutoffs PB.stat PB.cov.png;"
        "python3 /apps/{params.module}/scripts/hist_plot.py -c cutoffs_adjusted PB.stat adjusted_PB.cov.png;"
      )
