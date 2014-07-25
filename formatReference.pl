#!/usr/bin/perl

# Read in tab-delited file of ACCEPTOR/DONOR seqIDs and sequences, and output only the # of basepairs specified

### INPUT FILES
my $in = $ARGV[0];    
my $num = $ARGV[1];
my $out = $ARGV[2];
print $num;

open(IN, "< $in") or die "cannot open $in : $!\n";
open(OUT, "> $out") or die "cannot open $out : $!\n";

while($line=<IN>){

  chomp $line;
  next if $line =~ m/seq/i;
  @fields = split("\t",$line);
  
  $mySeq = substr($fields[1],0,$num);
  printf OUT ">%s\n%s\n", $fields[0],$mySeq;

 }
 close(IN);
 close(OUT);
 