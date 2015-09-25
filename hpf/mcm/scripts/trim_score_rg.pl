#!/usr/bin/perl
use List::Util 'shuffle';
use Getopt::Std;

$useage = "-f hpf2 filename -s number of decoys with top scores, 0 if no filter -r number of decoys with top rg's, 0 if no filter -x number of random decoys to select";

%options=();
getopts("f:s:r:x:",\%options);
# like the shell getopt, "d:" means d takes an argument
if ( defined $options{f} ) { $infile = $options{f} } else { die "$useage" };
if ( defined $options{s} ) { $score_num = $options{s} } else { $score_num = 20000 };
if ( defined $options{r} ) { $rg_num = $options{r} } else { $rg_num = 15000 };
if ( defined $options{x} ) { $rand_num = $options{x} } else { $rand_num = 0 };
print "Unprocessed by Getopt::Std:\n" if $ARGV[0];



$i = 0;
$tmpdir = "./";


my %rg_hash = ();
my %score_hash = ();
my %header_hash = ();
my %decoy_hash = ();

my $sequence = "";
my $currentDesc = "";
#kdrew: create array to store decoys
my $currentDecoy_ref;

print "opening file: $infile \n";
open(INFO, $infile);		# Open the file
#foreach $line (@lines)
while(<INFO>)
{
	my($line) = $_;
	if( $line =~ /^SCORE.+/ )
	{
		#print "score: " . $line;

		@arrayD = @$currentDecoy_ref;
		if($#arrayD >= 0)
		{
			#print "decoy not empty\n";
			#kdrew: store former decoy array
			$decoy_hash{ $currentDesc } = $currentDecoy_ref;
			#print "decoy_hash: @{$decoy_hash{ $currentDesc }} \n";
			my @currentDecoy = ();
			$currentDecoy_ref = \@currentDecoy;
		}


		#kdrew: split line by spaces to extract out score, rg and description
		my @split_line = split(/\s+/, $line);
		my $score = $split_line[1];
		my $rg = $split_line[12];
		
		#patrick hpf1 files have a different header line, but description is always last
		#my $description = $split_line[32];
		$size = @split_line;
    	my $description = $split_line[$size-1];
		
		
		#print "score: " . $score . " rg: " . $rg . " description: " . $description . "\n";
		$currentDesc = $description;

		#kdrew: hash on description
		$header_hash{ $currentDesc } = $line;
		$rg_hash{ $currentDesc } = $rg;
		$score_hash{ $currentDesc } = $score;
		
	}
	elsif ( $line =~ /^SEQUENCE/ )
	{
		#print "sequence: " . $line;
		$sequence = $line;
	}
	elsif ( $line =~ /.+/ )
	{
		#print "decoy: " . $line;
		my @split_line = split(/\s+/, $line);
		my $description = $split_line[13];
		#print "description: " . $description . "\n";
		if($description == $currentDesc)
		{
			#print "descriptions match\n";
			push(@$currentDecoy_ref, $line);
		}
		else
		{
			#print "problem with descriptions\n";
		}
	}
}

@arrayD = @$currentDecoy_ref;
if($#arrayD >= 0)
{
	#print "decoy not empty\n";
	#kdrew: store former decoy array
	$decoy_hash{ $currentDesc } = $currentDecoy_ref;
	#print "decoy_hash: @{$decoy_hash{ $currentDesc }} \n";
	my @currentDecoy = ();
	$currentDecoy_ref = \@currentDecoy;
}

#kdrew: cleanup, remove the header from all the hashes
delete $score_hash{ "description" };
delete $rg_hash{ "description" };
delete $decoy_hash{ "description" };



#kdrew: if the option was passed in as 0 don't filter
if(0 != $rand_num)
{
	#patrick shuffling the keys to randomize them.
	foreach $value (shuffle(keys %decoy_hash))
	{
		#print "rg: $value $rg_hash{$value}\n";
		$rand_num--; 

		if(0 > $rand_num)
		{
			delete $score_hash{ $value };
			delete $rg_hash{ $value };
			delete $decoy_hash{ $value };
		}
	}
}

#kdrew: if the option was passed in as 0 don't filter
if(0 != $score_num)
{
	foreach $value (sort {$score_hash{$b} cmp $score_hash{$a} } keys %score_hash)
	{
		#print "score: $value $score_hash{$value}\n";
		$score_num--;

		if(0 > $score_num)
		{
			delete $score_hash{ $value };
			delete $rg_hash{ $value };
			delete $decoy_hash{ $value };
		}
	}
}

#kdrew: if the option was passed in as 0 don't filter
if(0 != $rg_num)
{
	foreach $value (sort {$rg_hash{$a} cmp $rg_hash{$b} } keys %rg_hash)
	{
		#print "rg: $value $rg_hash{$value}\n";
		$rg_num--; 

		if(0 > $rg_num)
		{
			delete $score_hash{ $value };
			delete $rg_hash{ $value };
			delete $decoy_hash{ $value };
		}
	}
}

open(OUTFILE, ">$infile\.trim");
print OUTFILE $sequence;
print OUTFILE $header_hash{ "description" };
while ( my ($key, $value) = each(%decoy_hash) ) 
{
	#print "$key => $header_hash{ $key }\n";
	#print "$key => @$value\n";
	print OUTFILE $header_hash{ $key };
	print OUTFILE @$value;
}


close(OUTFILE);
close(INFO);			# Close the file
