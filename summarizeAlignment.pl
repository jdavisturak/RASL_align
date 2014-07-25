#!/usr/bin/perl

### 1) INPUT FILES
my $filebase = $ARGV[0]; 

$singleA = readAndSum($filebase . ".singleAcceptors",1);
$singleD = readAndSum($filebase . ".singleDonors",1);
$mislig = readAndSum($filebase . ".misligations",1);
$matches = readAndSum($filebase . ".matches",3);

printf "Singleton Acceptors:\t%d\nSingletonDonors:\t%d\nMisligations:\t%d\nMatches:\t%d\n",$singleA,$singleD,$mislig,$matches;

sub readAndSum{
  $file = shift;
  $column = shift;
  $counts = 0;
  
  open(IN, "< $file") or die "cannot open $file : $!\n";
  while(defined($line=<IN>)){    
    chomp $line;
    @fields = split("\t",$line);
    $counts += $fields[$column];
#    print $counts;
  }
  return $counts;
}
    
