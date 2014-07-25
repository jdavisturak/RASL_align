#!/usr/bin/perl -w
# This is my master pipeline for aligning RASL data

# Step 1: parse/decode

# Step 2: run bowtie (on each file)

# Step 3: Evaluate matches/mismatches (on each file)

# Step 4: Make summary table 

use strict;
use Getopt::Mixed;
use List::MoreUtils qw/ uniq /;
use Switch;

my $helpString;
$helpString = '
### INPUTS: 
#
## REQUIRED:
# -i/input: - name of the fastq file containing all of the reads
# -t/targets: A tab-delimited file containig these fields : (Well    Sequence        Sample), including header line.
# -A/acceptors:   location of Acceptors library  
# -D/donors:   location of Donors library
# -p: File of legal pairs: a tab-delimited file containig these fields : (GeneName        pairID  pairType), including header line
#
## OPTIONAL:
#
# [[Control of which steps are performed]]
#  -r 1234... runs everything, eg. 1234
#
# -o/output (optional):  string that all output files will be built off of
# -B:   Barcode length - number of characters of the barcode that are being used. (default: the length of Sequence[0]) in targets.
# -b:   offset for finding the barcode of length B in the targets file   (default: 3)
# -a:   length of Acceptor sequence (default: 20)
# -d:   length of Donor sequence (default: 20);
# -h:  help.
# -q/quiet:  Suppress output
';

##############################################################################################################################
## PROCESS INPUTS

my( $INPUT, $TARGETS, $LIB_A, $LIB_D, @STEPS, $fileBase, $barLength, $barOffset,$A_length,$D_length,$VERBOSE, @QUIET,$PAIRS ) =
  ( undef,   undef,   undef,    undef,  undef, undef, undef, 3, 20, 16, 1,undef, undef );
$barOffset = 3;
$A_length = 20;
$D_length = 16;
$VERBOSE = 1;
@QUIET = ();
@STEPS = ();

my $err = undef;


Getopt::Mixed::init( q{i=s input>i fastq>i
                    t=s targets>t
                    A=s acceptors>A
                    D=s donors>D
                    r=i R>r run>r Run>r
                    o=s out>o output>o
                    B=i   b=i  a=i  d=i  h  q
                    p=s P>p pairs>p
                    });

OPTS: {

  while( my( $option, $value, $pretty ) = Getopt::Mixed::nextOption() ) {
        if ($option eq 'h'){
          print STDERR $helpString;
          exit;
          }

        $INPUT = $value if $option eq 'i';
        $TARGETS = $value if $option eq 't';
        $LIB_A = $value if $option eq 'A';
        $LIB_D = $value if $option eq 'D';
        @STEPS = sort {$a <=> $b}  uniq(split("",$value)) if $option eq 'r';
        $fileBase = $value if $option eq 'o';
        $barLength = $value if $option eq 'B';
        $barOffset = $value if $option eq 'b';
        $A_length = $value if $option eq 'a';
        $D_length = $value if $option eq 'd';
        $VERBOSE = 0 if $option eq 'q';
        @QUIET = ("--quiet") if $option eq 'q';
        $PAIRS = $value if $option eq 'p';        
  }

  $err .= "   Must specify a value for targets file (-t)\n" if !defined($TARGETS);               
}

##############################################################################################################################
## PRE-PROCESS

@STEPS = uniq(@STEPS);
$fileBase = $INPUT if !defined($fileBase);

print "Scheduled actions:\n"  if $VERBOSE;
foreach (@STEPS){

  if($VERBOSE){ ### List what we are going to do:
    print "STEP 1: Decode fastq file.\n" if $_ eq 1; 
    print "STEP 2: Run Bowtie.\n" if $_ eq 2; 
    print "STEP 3: Compile Matches.\n" if $_ eq 3;  
    print "STEP 4: Generate summary statistics.\n" if $_ eq 4;  
  }
  
  ### Check that we have the correct info to proceed:
  $err .= "   Must specify a value for fastq file (-i)\n" if !defined($INPUT) and $_ eq 1;  
  $err .= "   Must specify a value for Acceptors index (-A)\n" if !defined($LIB_A) and $_ eq 2;               
  $err .= "   Must specify a value for Donors index (-D)\n" if !defined($LIB_D)  and $_ eq 2;                 
  $err .= "   Must specify a Pairs file (-p)  \n" if !defined($PAIRS) and ($_ eq 3);               
  $err .= "   Must specify an input file (-i) or output string (-o)  \n" if !defined($fileBase) and ($_ eq 2 or $_ eq 3 or $_ eq 4);               

}

die "Error: failed to specify correct command line options:\n$err" if defined($err);

##############################################################################################################################
## PROCESS STEPS 

my ($step, $samples_ref,$samples,$acceptorFiles, $donorFiles,$acceptors_aligned,$donors_aligned, $args, $sampleBase) = 
  (undef, undef,undef, undef, undef, undef, undef, undef, undef);
my $logFile = $fileBase. ".log";

## READ IN THE TARGETS
($samples, $acceptorFiles, $donorFiles, $acceptors_aligned, $donors_aligned) = getSamples($TARGETS,$fileBase); 

## LOOP through all STEPS
foreach $step (@STEPS){
  switch ($step){     # Step 1: parse/decode

    case 1 {
      ## Decode fastq file
      
      ## Set up command:
      $args = "perl /Users/jeremydt/RASL/parseRASL1.pl $INPUT $TARGETS $fileBase $D_length";
      print "Running command: ----> $args \n" if $VERBOSE;
      system($args) == 0  or die "System Command $args failed: $?";
    
    }
    
    case 2 {    # Step 2: run bowtie (on each file)
        
        for(my $i=0;$i < scalar(@$samples);$i ++){          

          # Align donors                 
          $args = sprintf("echo '\n<D> Sample %s donors:' >> %s",$samples->[$i], $logFile);
          system($args) == 0  or die "System Command $args failed: $?";          
          $args = sprintf("bowtie -v 1 -m 1 --best -f -p 5 %s %s %s 2>> %s",  $LIB_D, $donorFiles->[$i] ,  $donors_aligned->[$i],$logFile);          
          print "Running command: ----> $args \n\n" if $VERBOSE;
          system($args) == 0  or die "System Command $args failed: $?";
       
          # Align acceptors  
          $args = sprintf("echo '\n<A> Sample %s acceptors:' >> %s",$samples->[$i], $logFile); ## Write a header in the log files
          system($args) == 0  or die "System Command $args failed: $?";                        
          $args = sprintf("bowtie -v 1 -m 1 --best -f -p 5 %s %s %s 2>> %s",$LIB_A,  $acceptorFiles->[$i] ,  $acceptors_aligned->[$i],$logFile);
          print "Running command: ----> $args \n\n" if $VERBOSE;
          system($args) == 0  or die "System Command $`args failed: $?";
       
        }
    }
    
    case 3 {    # Step 3: Evaluate matches/mismatches (on each file)
           
        for(my $i=0;$i < scalar(@$samples);$i ++){                    
          my $sampleBase =  $fileBase . "_" . $samples->[$i];
          $args = "perl /Users/jeremydt/RASL/matchAlignmentsToPairs.pl $sampleBase $PAIRS";
          print "Running command: ----> $args \n\n" if $VERBOSE;
          system($args) == 0  or die "System Command $args failed: $?" ;
        }
    
    }
    
    case 4{     
        #### Step 4: Make summary table 
        ## Need raw reads, # matches, # mislgations, singleA, singleD;
        # also no match, bad N
         
        my $summaryFile = $fileBase . ".summary.txt";
        my ($Ns,$Wrong,$Good) = (undef,undef,undef); 
        ### ALL OF THIS IS BASED ON THE LOG FILE, WHICH IS BASED ON THE ORDER OF THE HASH KEYS.
        open(LOGFILE, "< $logFile") or die "cannot open $logFile : $!\n";
        open(SUMMARYFILE, "> $summaryFile") or die "cannot open $summaryFile : $!\n";;
        
        my $sampleCount = scalar(@$samples);
        my %samplesHash = ();
        
        while(<LOGFILE>){                   
          chomp;
          if ($sampleCount--){        # for each sample, get the decode counts from the parseRASL.pl output
            my ($sample, $code, $reads) = split("\t",$_);             # Total reads: '$reads';    
            $samplesHash{$sample} = $reads;
          }
          
          else{
            ## Parse the summary output from parseRASL.pl
            my $m = $_;
            $m =~ s/N: ([0-9]*), Wrong code: ([0-9]*) Good Codes: ([0-9]*)/$1 $2 $3/; # This is now Ns, Wrong codes, Good codes
            ($Ns,$Wrong,$Good) = split(" ",$m);
            close (LOGFILE);  last;            # ignore the rest of the lines
          }
        }      
        
       
        #
        
        ## Construct final table
        print SUMMARYFILE "Sample\tReads\tMatches\tPercent Matches\tMisligations\tSingleton Acceptors\tSingleton Donors\n";
        
        my $totalReads = 0;
        my $totalMatches = 0;
        my $samp;        
                
        foreach $samp (@$samples){           
           my $rawReads = $samplesHash{$samp};
           my $matches = readAndSum($fileBase . "_$samp.matches",3);
           my $percentMatch = asPercent($matches, $rawReads );
           my $misLig = readAndSum($fileBase . "_$samp.misLigations",1);
           my $singleA = readAndSum($fileBase . "_$samp.singleDonors",1);
           my $singleD = readAndSum($fileBase . "_$samp.singleAcceptors",1);

           print SUMMARYFILE join("\t",($samp,$rawReads, $matches, $percentMatch, $misLig, $singleA, $singleD)) ."\n";
           
           $totalReads += $rawReads;
           $totalMatches += $matches;
        }
        my $totalRaw = ($totalReads + $Wrong + $Ns);
        
        print SUMMARYFILE "Total decoded\t$totalReads\t$totalMatches\t" . asPercent($totalMatches,$totalReads) . "\n";
        print SUMMARYFILE "Bad decode\t$Wrong\n";
        print SUMMARYFILE "Contain N\t$Ns\n";
        print SUMMARYFILE "Total Raw\t$totalRaw\tDecode Rate\t" . asPercent($totalReads, $totalRaw) ."\tRaw Match Rate\t" . asPercent($totalMatches, $totalRaw) ."\n";                        
    }  
  }  
}

### SUBROUTINES ###

## Read in samples
sub getSamples{
  my $TARGETS = shift;
  my $fileBase = shift;
  
  my (@samples,@acceptorFiles,@donorFiles,@acceptor_aligned,@donor_aligned) = ((),(),(),(),());
  
  open(INPUT_TARGETS, "< $TARGETS") or die "Cannot open $TARGETS \n";
  
  my ($line,$foo,$bar,$SAMPLE);
  while(defined($line= <INPUT_TARGETS>)){
  
  	chomp $line;
  	next if ($line =~ m/Seq/i);
   	($foo,$bar,$SAMPLE) = split("\t",$line);
  	
  	push(@samples,$SAMPLE);
    push(@acceptorFiles, ($fileBase  . "_" . $SAMPLE . "_acceptors.fa") );
    push(@donorFiles, ($fileBase  . "_" . $SAMPLE . "_donors.fa") );      
    push(@acceptor_aligned, ($fileBase  . "_" . $SAMPLE . "_acceptors.aligned") );
    push(@donor_aligned, ($fileBase  . "_" . $SAMPLE . "_donors.aligned") );
  }
  close (INPUT_TARGETS);
  
  return (\@samples,\@acceptorFiles,\@donorFiles,\@acceptor_aligned,\@donor_aligned); # this is creating a referenced array; also, it contains only references!  Oy...
}

## Read in a file and sum one column
sub readAndSum{
  my $file = shift;
  my $column = shift;
  my $counts = 0;
  
  open(IN, "< $file") or die "cannot open $file : $!\n";
  while(<IN>){    
    chomp;
    my @fields = split("\t",$_);
    $counts += $fields[$column];
  }
  return $counts;
  close (IN);
}

## Return a formatted percentage from 2 numbers
sub asPercent{
  my $num = shift;
  my $denom = shift;
  return 0 if $denom eq 0;
  return sprintf("%.2f",$num/$denom*100);
}