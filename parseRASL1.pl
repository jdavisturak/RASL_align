#!/usr/bin/perl -w

#First step:
#keep track of the # of reads:
#	Total
#	Bad:
#		Ns in the way of decoding, acceptor or donor?
#		Otherwise bad decoding
#	Total for each lane...
#	Write to each lane's new files .

### INPUT FILES
my $raslfile = $ARGV[0];          # Data
my $targetsfile = $ARGV[1];       # Samples names/barcodes
my $fileBase = $ARGV[2];          # For naming other files
my $donorLength = $ARGV[3];
$donorLength = 16 if !defined($donorLength);
$acceptorLength = 20;

$fileBase = $raslfile if !defined($fileBase);           

open(INPUT_RASL, "< $raslfile") or die "cannot open $raslfile  : $!\n";
open(INPUT_TARGETS, "< $targetsfile") or die "cannot open $targetsfile : $!\n";;

### OUTPUT FILES
my $outputfile = $fileBase . ".log";    # log, statistics
my $badCodesFile = $fileBase. ".badCodes"; # Reads that are barcoded incorrectly: for future use, if desired.
open(OUTPUT, "> $outputfile");
open(BADCODES, "> $badCodesFile");

### INITIALIZE VARIABLES
my %sampleNames = ();     #Hash of sample names, with key=code
my %targetCounts = ();    #Hash (# reads), with key=sampleName
my %acceptorFiles = ();     #Hash (with key=code) of target files to put the reads, for alignment  
my %donorFiles = ();     #Hash (with key=code) of target files to put the reads, for alignment  

my $totalCount = 0;
my $Ncount = 0;
my $wrongCode = 0;

# Read in targets (Sample), initialize hashes
while(defined($line= <INPUT_TARGETS>)){
	chomp $line;
	next if ($line =~ m/Seq/i);
 	@fields = split("\t",$line);
	
	### Get multiplex code (6 bp of unique code), save sample name
	$code = substr(revcomp($fields[1]),2,6);
	$sampleNames{$code} = $fields[2];
	$targetCounts{$sampleNames{$code}} = 0;
	
  ### Initialize handles to output files	
	$acceptorOutFile =  $fileBase  . "_" . $fields[2] . "_acceptors.fa";
	$donorOutFile    =  $fileBase  . "_" . $fields[2] . "_donors.fa";

  open($acceptorFiles{$code}, "> $acceptorOutFile");
	open($donorFiles{$code}, "> $donorOutFile");
}

#### read in RASL file ...
while(defined($line = <INPUT_RASL>)) 
{
	
	$line = <INPUT_RASL>;  ## This second line is the real sequence
	$garbage = <INPUT_RASL>;
	$garbage = <INPUT_RASL>;

	### retrieve sequence elements:
	$acceptorSeq = substr($line,0,$acceptorLength);	
	$donorSeq = substr($line,$acceptorLength,$donorLength);
	$codeSeq = substr($line,$acceptorLength+$donorLength,6);
  
	### Check for N
	if($acceptorSeq =~m/N/ or $donorSeq =~ m/N/ or $codeSeq =~m/N/){
		$Ncount++;
		next;
	}
	
	### Check code
	if (defined($sampleNames{$codeSeq})){ 
		$goodCodes++;
		### OK, we're in business!
		$targetCounts{$sampleNames{$codeSeq}}++;

		##### Write output to acceptor, donor files 
		## First generate a fastq name for this read: @Sample_num_(A/D)
		$myFastq =  sprintf(">%s_%d_", $sampleNames{$codeSeq},$targetCounts{$sampleNames{$codeSeq}}); 
		
		## Print to acceptor, then donor files
	  print {$acceptorFiles{$codeSeq}} sprintf("%sA\n%s\n" , $myFastq, $acceptorSeq);
    print {$donorFiles{$codeSeq}} sprintf("%sD\n%s\n" , $myFastq, $donorSeq);
	}
	else{
    $wrongCode++;
    print BADCODES $codeSeq . "\n";
    next;
	}
}
	
	
### Print out summary	
# Total number of reads for each sample:
while( my ($code, $sample) = each %sampleNames ) {     
    print OUTPUT "$sample\t$code\t" . $targetCounts{$sample} . "\n"; 
    close($acceptorFiles{$code});
    close($donorFiles{$code});
} 
# Summary stats for the lane:
print OUTPUT "N: " . $Ncount . ", Wrong code: " . $wrongCode . " Good Codes: " , $goodCodes . "\n";  

close(BADCODES);
close(OUTPUT);

################################################################################################

sub revcomp {
  my $dna = shift;
  my $rev = reverse($dna);
  $rev =~ tr/ACGTacgt/TGCAtgca/;

  return $rev;
}

