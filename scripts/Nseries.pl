#!/usr/bin/env perl
use strict;
use Getopt::Long;
my $genome_size = 0; 
my $ass = '';
my $std = 0;
my $sort = 1;
my $fosmid_sim = 0;
my $aln = 0;
GetOptions(
	   'g|s|l:s'  => \$genome_size,
	   'id|a:s'    => \$ass,
	   'sort!' => \$sort,
	   'fos' => \$fosmid_sim,
	   'aln|alignment' => \$aln
	  );

my $numblocks = 0;
my $sum = 0;
my $first = 1;
my $largest = 0;
my $smallest = 0;
my $assembly_size = 0;
my $cov = 0;
my $ps = 0;
open OUT,"|transpose.pl";
#open OUT,*STDOUT;
if ($fosmid_sim) {
  $ass=~/p(\d+)_(\d+)/;
  $cov = $2;
  $ps = $1;
}
my @nseries = map{5*$_}(0 .. 24);

my %N;
my %NG;
my %nbN;
my %nbNG;


#initialize all values to 0;
foreach my $s (@nseries) {
  $N{$s} = 0;
  $nbN{$s} = 0;
}
foreach my $s (@nseries) {
  $NG{$s} = 0;
  $nbNG{$s} = 0;
}

my @lens;
while (<>) {			#read reverse sorted input
  chomp;
  my ($length,$rest) = split;
  $assembly_size+=$length;
  push @lens, $length;
}
my @sorted;
if ($sort) {
  @sorted = sort numerically @lens;
} else {
  @sorted = @lens;
}
if (!$genome_size) {
  #$genome_size = $assembly_size;
}
#print  "Assembly\t$ass\n" if $ass;
foreach my $size (@sorted) {	#read reverse sorted input
  $numblocks++;
  if ($first) {
    $largest = $size; $first = 0;
    $N{0}=$largest;
    $NG{0}=$largest;
    $nbN{0}=1;
    $nbNG{0}=1;
  }
  $sum+=$size;
  my $L=$size;

  foreach my $s (@nseries) {
    if (!$N{$s} &&  $sum>=$assembly_size*($s/100)) {
      $N{$s} = $L;
      $nbN{$s}=$numblocks;
    }
  }
  if ($genome_size) {
    foreach my $s (@nseries) {
      if (!$NG{$s} &&  $sum>=$genome_size*($s/100)) {
	$NG{$s} = $L;
	$nbNG{$s}=$numblocks;
      }
    }
  }
} 
$smallest = $sorted[-1];
print OUT "metric\tname";
if ($fosmid_sim) {
  print OUT "\tpsize\tcov";
}
print OUT "\tL_gen\tL_ass\tnum_blocks";
foreach my $s (@nseries) {
  print OUT "\t$s";
}
print OUT "\n";
my $end = 0;
my $A = '';
if ($aln) {
  $A = 'A';
}
if ($genome_size) {
  if ($fosmid_sim) {
    print OUT "NG$A\t$ass\t$ps\t$cov\t$genome_size\t$assembly_size\t$numblocks";
  } else {
    print OUT "NG$A\t$ass\t$genome_size\t$assembly_size\t$numblocks";
  }
  
  foreach my $s (@nseries) {
    if (!$NG{$s}) {
      if ($end) {
	print OUT "\t-";
      } else {
	print OUT "\t$smallest";
	$end = 1;
      }
    } else {
      print OUT "\t$NG{$s}";
    }
  }
  print OUT "\n";
  if ($fosmid_sim) {
    print OUT "nbNG$A\t$ass\t$ps\t$cov\t$genome_size\t$assembly_size\t$numblocks";
  } else {
    print OUT "nbNG$A\t$ass\t$genome_size\t$assembly_size\t$numblocks";
  }
  $end = 0;
  foreach my $s (@nseries) {
    if (!$nbNG{$s}) {
      if ($end) {
	print OUT "\t-";
      } else {
	print OUT "\t$numblocks";
	$end = 1;
      }
    } else {
      print OUT "\t$nbNG{$s}";
    }
  }
  print OUT "\n";
}
if ($fosmid_sim) {
  print OUT "N$A\t$ass\t$ps\t$cov\t-\t$assembly_size\t$numblocks";
} else {
  print OUT "N$A\t$ass\t-\t$assembly_size\t$numblocks";
}
foreach my $s (@nseries) {
  if ($s > 100) {
    print OUT "\t-";
  } else {
    print OUT "\t$N{$s}";
  }
}
print OUT "\n";
if ($fosmid_sim) {
  print OUT "nbN$A\t$ass\t$ps\t$cov\t-\t$assembly_size\t$numblocks";
} else {
  print OUT "nbN$A\t$ass\t-\t$assembly_size\t$numblocks";
}
foreach my $s (@nseries) {
  if ($s > 100) {
    print OUT "\t-";
  } else {
    print OUT "\t$nbN{$s}";
  }
}
print OUT "\n";


sub numerically { $b <=> $a }
