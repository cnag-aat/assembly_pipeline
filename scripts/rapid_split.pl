#!/usr/bin/env perl
# kj2 08.10.2020

use strict;
use Getopt::Long;

my $fa;
my $printfa;
my $split;
my $help;

GetOptions (
    "fa:s" => \$fa,
    "printfa" => \$printfa,
    "h" => \$help,
    "help" =>   \$help,
    "split" => \$split,
);

if (($help) || (!$fa)) {
    print "This script takes a fasta file and produces a TPF on contig level, i.e. it splits at all Ns\n";
    print "Usage:\n";
    print "perl split.pl -fa <fasta>\n";
    print "              -printfa \# if you want to print out the split fasta\n";
    print "              -h/help  \# this message\n";
    exit(0);
}

open(OUT,">./${fa}.split") if ($printfa);
open(OUA,">./${fa}.tpf");
open(SPL,"seqtk cutN -n 1 $fa |");
my $lastcoord;
my $last;
while (<SPL>) {
    print OUT if ($printfa);
    /\>((\S+)\:(\d+)-(\d+))/ and do {
        if ($last && ($last eq $2)) {
            my $length = $3-$lastcoord-1;
            print OUA "GAP\tTYPE-2\t$length\n";
        }
        print OUA "?\t$1\t$2\tPLUS\n";
        $lastcoord = $4;
        $last = $2;
    }
}
if ($split) {
    die("you need to also set printfa\n") unless ($printfa);
    system("mkdir temp") unless (-e "temp");
    chdir("temp");
    system("perl ~kj2/bin/mod_fasta.pl -single ../${fa}.split");
}

