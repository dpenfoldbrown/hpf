#!/usr/bin/perl -w

# Richard Bonneau , Initl/ auth.
# copywrite 2002 -  Institue for Systems Biology, 
# Seattle, WA 98103


# this combines silent mode files from Rosetta ... very simple!
# these are the result of runnig rosetta -silent
# see the readme that comes with the distribution

$maxFiles = shift @ARGV;
@files = @ARGV;

unless( $maxFiles) {
	die "silentcombine.pl MAX-FILES FILE1 FILE2 ...\n"
}

$first = 1;
$nfiles = 0;

foreach $f (@files) {
	open(INFILE,"$f") || warn "$f not a valid file\n"; 
	$nfiles++;
	if ($nfiles >= $maxFiles) {
		warn "skipping $f ... max-files exceeded...\n";
		next;
	}
	
	if ($first  == 1) {
		$line = <INFILE>; ## seq line
		chomp $line;
		@w = split /\s+/,$line;
		$initialSeq = $w[1]; ## keep initial seq for checking
		print "$line\n";

		$line = <INFILE>; ## score decription
                print $line;
		
		$first = 0;
	} else {
		$line = <INFILE>; # skip and check seq line
		chomp $line;
		@w = split /\s+/,$line;
		$thisSeq = $w[1];
		if ($thisSeq ne $initialSeq) {
			die "SEQUENCE MISMATCH AT $f ... exiting ...\n";
		}
		$line = <INFILE>; # skip
	}
	while ($line = <INFILE>) {
		print $line;
	}

}



