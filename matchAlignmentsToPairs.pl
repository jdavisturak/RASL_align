#!/usr/bin/perl

## This script takes a sample name as an argument and matches the aligned acceptor and donors, trying to find a pair of legit ligations to match to  <filename>.donors.aligned and <filename>.acceptors.aligned

### Pseudo code:

## 1) Read in data, make pairs:
# read in donor file, put it into a hash

# read in acceptor file:
#    as it is being read in, compare to donor hash and try to find the same reads.  
#      If a pair is found, save that pair to the 'pairs' hash, and delete it from the donor hash.
#      Otherwise, if THIS ACCEPTOR has no pair, add it to the lonely acceptor array...

## 2) Count the singletons
# For donors in the existing hash, create a new hash where the keys are seqIDs.  Increment this when one is found.
# For acceptors in the array, do the same thing.

## 3) Count the pairs.
# Read in the list of legal pairs, for each one, look for its occurence in the 'pairs' hash.  Record it (even if 0), then delete entries in the pair hash.

## OUTPUT:
# Legal pairs hash
# Bad pairs=misligations hash
# Singletons: Acceptors, Donors (1 or two files?)


### 1) INPUT FILES
my $filebase = $ARGV[0]; 

## Read in donor file
my $donor_ref = getDonors($filebase . "_donors.aligned");

## Read in acceptor file, and figure out pairs and singletons
my $acceptor_output_ref = readAndCompareAcceptors($filebase . "_acceptors.aligned",$donor_ref); # returns: \%pairs_hash,\%singleAcceptors]; # 
$pairs_ref = $acceptor_output_ref->[0];
$singleAcceptors_ref = $acceptor_output_ref->[1];

### 2) Count the singletons
# For donors in the existing hash, create a new hash where the keys are seqIDs.  Increment this when one is found.
$singleDonors_ref = tabulateSingleDonors($donor_ref);

### 3) Tabulate the pairs (also includes output, instead of keep it all in memory)
processGoodPairs($ARGV[1],$filebase .".matches", $pairs_ref);

## OUTPUT non-matches:
outputHashRef($filebase .".misligations",$pairs_ref);    # misligations = what's left in pairs_ref
outputHashRef($filebase .".singleAcceptors",$singleAcceptors_ref);
outputHashRef($filebase .".singleDonors",$singleDonors_ref);

  
################################################################################################
################################################################################################
################################     FUNCTIONS     ############################################
################################################################################################
################################################################################################

################################################################################################
# Generic function to print the contents of a hash reference to tab-delimited file
sub outputHashRef{
  $out_file = shift;
  $ref = shift;
  
  open(OUT, "> $out_file") or die "cannot open output file $out_file : $!\n";
  while( ($k,$v) = each(%$ref)){
    print OUT "$k\t$v\n";
  }
  close(OUT);
}
################################################################################################
# This loads a file of legal probe pairs, then looks in 'pairs_ref' to report how often each pair is seen. It also deletes seen pairs
sub processGoodPairs{
  $database_file = shift;
  $out_file = shift;
  $pairs_ref = shift;
 
  open(in_handle, "< $database_file") or die "cannot open database file $database_file : $!\n";
  open(OUT, "> $out_file") or die "cannot open output $out_file  (2): $!\n";

  while(defined($line=<in_handle>)){    
    chomp $line;
    @fields = split("\t",$line);
    next if $fields[0] =~ m/GeneName/;
    
    $pairsID = $fields[1];
    $pairsValue = 0 + $pairs_ref->{$pairsID};

    # Simply output the row here
    print OUT "$line\t$pairsValue\n";

    # Delete this pair from hash:
    delete($pairs_ref->{$pairsID});
  }
  close(OUT);
  close(in_handle);
  
}
################################################################################################

sub tabulateSingleDonors{
  $donor_ref = shift;
  my %singleDonors = ();
  for (keys %$donor_ref){
    $seqID = $donor_ref->{$_};
    $singleDonors{$seqID}++;
  }
  return \%singleDonors;
}
################################################################################################

sub readAndCompareAcceptors{
  $acceptorfile = shift;
  $donor_ref = shift;

  open(acceptor_handle, "< $acceptorfile") or die "cannot open acceptor file $acceptorfile : $!\n";

  #my %donor_hash = %$donor_ref;    # This hash is the readID=>seqID
  my %pairs_hash = ();             # This hash will be pairID=>count

  #print "Reading in acceptors ... ";

  my %singleAcceptors = ();
  while(defined($line=<acceptor_handle>)){
 
    chomp $line;
    @fields = split("\t",$line);

    next if $fields[6] > 0; # This is the number of accepted alignments for this read

    $readID = substr($fields[0],0,length($fields[0])-2); # Name of the read: I appended '_D' to all donor sequence names
    $acceptor_seqID  = $fields[2];	

    ## Retrieve, if possible, the where this donor aligned.
    #$donor_seqID = $donor_hash{$readID};
    $donor_seqID = $donor_ref->{$readID};

    if (defined($donor_seqID)){
      # This means there IS a value in the donor hash for this read: now we have a pair.
      $pairID = $donor_seqID . "_" . $acceptor_seqID;

      $pairs_hash{$pairID}++;    # Increment the count of this pairID
      #delete($donor_hash{$readID}); # get rid of this read in the donors' singleton list.
      delete $donor_ref->{$readID};
    }
    else{
      # This means the acceptor is the singleton
      $singleAcceptors{$acceptor_seqID}++;
    }
  }
  close(acceptor_handle);
  #print "Done.\n";

  ## OUTPUT: references to the hashes - hopefully they won't get destroyed?
  return [\%pairs_hash,\%singleAcceptors]; # this is creating a referenced array; also, it contains only references!  Oy...
}

################################################################################################

sub getDonors{
  $donorfile = shift;
  open(donor_handle, "< $donorfile") or die "cannot open donor file $donorfile : $!\n";

  #print "Reading in donors ... ";

  # read in donor file, put it into a hash
  my %donors = ();
  while(defined($line=<donor_handle>)){
    chomp $line;
    @fields = split("\t",$line);
    
    next if $fields[6] > 0; # This is the number of accepted alignments for this read
    $readID = substr($fields[0],0,length($fields[0])-2); # Name of the read: I appended '_D' to all donor sequence names
    $seqID  = $fields[2];	
    $donors{$readID} = $seqID;
    
  } 
  close(donor_handle);
  #print "Finished.\n";	
  return \%donors; # referenced
}  

################################################################################################

sub revcomp {
  my $dna = shift;
  my $rev = reverse($dna);
  $rev =~ tr/ACGTacgt/TGCAtgca/;

  return $rev;
}

