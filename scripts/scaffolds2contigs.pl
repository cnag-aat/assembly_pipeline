#!/usr/bin/env perl
use strict;
use lib "/home/devel/talioto/myperlmods";
use lib "/home/devel/talioto/myperlmods/Bio";
use SeqOp;
use Getopt::Long;
use File::Basename qw( fileparse );
my $fname = 0;
my $num_n = 1;
my $species = '';
my $taxid = '';
my $assembly_name = 0;
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my @abbr = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
my $realyear = 1900+$year;
my $date = "$mday-$abbr[$mon]-$realyear";
my $center = "CNAG";
my $description = "";
my $min_scaff_size = 20;
my $min_contig_size = 20;
my $sp = 0;
my $cp = 0;
my $lookuptable = 1;
my $keep = 0;
my $sort = 1;
GetOptions(
	   'f|s|i:s'        => \$fname,
	   'k|keep'         => \$keep,
	   'n:s'            => \$num_n,
	   'name:s'         =>\$assembly_name,
	   'taxid:s'        =>\$taxid,
	   'species:s'      =>\$species,
	   'center:s'       =>\$center,
	   'date:s'         =>\$date,
	   'mins:s'         =>\$min_scaff_size,
	   'minc:s'         =>\$min_contig_size,
	   'sp:s'           => \$sp,
	   'cp:s'           => \$cp,
	   'lookup!'       => \$lookuptable,
	   'sort!'			=> \$sort
	);
my $options = <<HERE;
           -f,-s,-i	name of fasta file, required
           -n		number of Ns on which to break (1)
           -name	assembly name (filename)
           -taxid	taxonomy id from NCBI
           -species	species
           -center	center (CNAG)
           -date	date of assembly (current local time)
           -mins	minimum scaffold size to keep (20)
           -minc	minimum contig size to keep (20)
HERE
die "usage: $0 -f <scaffolds.fa> [options]\n$options\n" if !$fname;
print STDERR "'-n 0' does not make sense. Running with '-n 1'\n" if !$num_n;
my ($base,$path,$ext) = fileparse($fname,qw(\.fa \.fasta));
if ($assembly_name){
  $base = $assembly_name;
}else{
  $assembly_name = $base;
}
my $cout = "$base.contigs.fa";
my $sout = "$base.scaffolds.fa";
my $offsets = "$base.contigs.offsets";
my $soffsets = "$base.scaffolds.offsets";
my $gaps = "$base.scaffolds.gaps.bed";
my $agp = "$base.agp";
my $lup = "$base.lookup.txt";
my @scaffolds;
open IN, "FastaToTbl $fname |";
my $maxc=1;
my $nums=0;
while (<IN>) {
  $nums++;
  my $c=1;
  my @F=split;
  push @scaffolds,{'len'=>length($F[1]),'id'=>$F[0],'seq'=>$F[1]};
    while ($F[1]=~/N+/g) {
      $c++;
    }
  if ($c>$maxc) {
    $maxc=$c;
  }
}
close IN;
if (!$sp) {
  $sp = length($nums);
}
if (!$cp) {
  $cp = length($maxc);
}

#open IN, "FastaToTbl $fname |";
unlink $sout if -e $sout;
unlink $cout if -e $cout;
open (AGP,">$agp");
print AGP '##agp-version 2.0',"\n";
print AGP '# ',"ORGANISM: $species\n";
print AGP '# ',"TAX_ID: $taxid\n";
print AGP '# ',"ASSEMBLY NAME: $assembly_name\n";
print AGP '# ',"ASSEMBLY DATE: $date\n";
print AGP '# ',"GENOME CENTER: $center\n";
print AGP '# ',"DESCRIPTION: $description\n";


open SOUT, "| TblToFasta >> $sout";
open COUT, "| TblToFasta >> $cout";
open OFF, ">$offsets";
open SOFF, ">$soffsets";
open GAP, ">$gaps";
if ($lookuptable){open LU,">$lup";}
my $scount = 1;
print STDERR "Running $0 on $fname\n";
#while (<IN>){
my @sorted_scaffolds=();
if ($sort){@sorted_scaffolds = sort lensort @scaffolds;@scaffolds=@sorted_scaffolds;}
foreach my $scaff (@scaffolds){
  chomp;
  my $sid = $scaff->{id};
  my $seq = $scaff->{seq};
  my $oldid = $sid; 
  my $seq_up = uc($seq);
  $seq = $seq_up;
  #print STDERR "$sid\n";
  my $num = $min_contig_size - 1;
  if ($num>0){
  	$seq=~s/^([^N]{1,$num})N/"N" x (length($1)+1)/ge;  
  	$seq=~s/N([^N]{1,$num})N/"N" x (length($1)+2)/ge; 
  	$seq=~s/N([^N]{1,$num})$/"N" x (length($1)+1)/ge;
  }
  if ($seq=~s/^(N+)//){ ### this creates new  offset
    my $scaffoffset = length($1);
    print SOFF "$oldid\t$scaffoffset\n";
  }else{
    print SOFF "$oldid\t0\n";
  }
  $seq=~s/N+$//;
  next if length($seq)<$min_scaff_size;

  my $outseq = $seq;
  $outseq=~s/[^ACGT]/N/g;
  next if $outseq!~/[ACGTacgt]/;
  if(!$keep){
    $sid = sprintf("%s"."_s%0$sp"."i",$base,$scount++);
  }
  print SOUT "$sid\t$outseq\n";
  if($lookuptable){print LU "$sid\t$oldid\n";}
  my $count = 1;
  # my @contigs = split /N+/,$seq;
  # foreach my $contig (@contigs){
    
  #   my $cid = sprintf("%s.contig%04d",$sid,$count++);
  #   print COUT "$cid\t$contig\n";
  # }
  my $END = 0;
  my $offset = 0;
  my $element = 0;
  my $part_number = 1;
  foreach my $contig_or_gap (split(/(N+)/i,$seq)){# foreach my $contig_or_gap (split(/(N{$num_n,})/i,$seq))
    #my $contig = $1;
    if ($element % 2){
      my $gapsize = length($contig_or_gap);
      print AGP join("\t",($sid,$offset+1,$offset+$gapsize, $part_number++,'N',$gapsize,'scaffold','yes','paired-ends')),"\n";
      $offset+=$gapsize;
    }else{
      $contig_or_gap=~s/[^ACGTacgt]/N/g; #replace IUPAC ambiguity codes;
      #my $offset = $end - length($contig);
      print GAP  "$sid\t$END\t$offset\n" if $count > 1;;
      my $cid = sprintf("%sc%0$cp"."d",$sid,$count++);
      my $contig_len = length($contig_or_gap);
      print COUT "$cid\t$contig_or_gap\n";
      print OFF  "$cid\t$offset\n";
      print AGP join("\t",($sid,$offset+1,$offset + $contig_len, $part_number++,'W',$cid,'1',$contig_len,'+')),"\n";
      $offset+=$contig_len;
      $END = $offset;
    }
    $element++;
  }
}
#close IN;
close COUT;
close SOUT;
close GAP;
close AGP;
close OFF;
if ($lookuptable){close LU;}

sub lensort
  {
    $b->{len} <=> $a->{len}
  }

