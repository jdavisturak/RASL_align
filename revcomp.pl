#!/usr/bin/perl -w

foreach  $arg (@ARGV){
	print( revcomp($arg) . "\n");
}

sub revcomp {
  my $dna = shift;
  my $rev = reverse($dna);
  $rev =~ tr/ACGTacgt/TGCAtgca/;

  return $rev;
}
