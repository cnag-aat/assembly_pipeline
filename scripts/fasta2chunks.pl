#!/usr/bin/env perl
use Getopt::Long;
my $description;
use Bio::SeqIO;
use Bio::Seq;
use File::Basename;
my $format = 'fasta';
#my $chunksize = 2500;
my $numchunks = 100;
my $file = 0;
my $pad = 0;
GetOptions(
	   "n:s"=>\$numchunks,
	   "f:s"=>\$file,
	   'pad'=>\$pad
	  );
my @suffixlist = qw(.gff .gff3 .gff2 .gtf .fa .fsa .fasta .gb);
my ($name,$path,$suffix) = fileparse($file,@suffixlist);
my $ext = 'fa';
die "no file given" if !$file;
if ($suffix eq ".gb") {
    $format = 'Genbank';
}
$in  = Bio::SeqIO->new( '-format' => $format , -file => $file);
#$file =~ s/\.[^.]+$//;

my @out = ();
my $c = 1;
my $n = 0;
my $j = 0;
my $fw = 2;
$fw = (int ((log $numchunks) / (log 10)))+1;

for (my $j=1;$j<=$numchunks;$j++) {
    my $cn = $j;
    if ($pad){$cn = sprintf("%0".$fw."d",$j);}
    push @out, "$name.$cn.$ext";
}

print STDERR "Writing to chunk files...\n";
while ($seqobj = $in->next_seq()) {
    my $seqio = Bio::SeqIO->new('-format' => $format,-file => ">>".$out[$n++%$numchunks]);
    $seqio->write_seq($seqobj);
}
