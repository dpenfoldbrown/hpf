#!/usr/bin/perl 

# Transforms a series of PSI-BLAST checkpoint files from a binary 
# format to an ASCII-readable format, and it corrects errors with
# the way that PSI-BLAST does pseudocounts. To save space, it will
# compress the checkpoint files after they're created.
#
# Copyright 2005, University of Washington
#   This document contains private and confidential information and
#   its disclosure does not constitute publication.  All rights are
#   reserved by University of Washington, except those specifically 
#   granted by license.
#
#
#  Initial Author: James Thompson (tex@u.washington.edu)
#  $Revision: 8054 $
#  $Date: 2006-04-27 22:54:39 -0400 (Thu, 27 Apr 2006) $
#  $Author: tex $
#
###################################################################

use File::Basename;

my $readblastmatrix  = $ENV{HOME} . '/src/msaUtil/readblastmatrix';
my $blosum_matrix    = $ENV{HOME} . '/src/msaUtil/blosum62.qij';
my $check_directory  = $ENV{HOME} . '/vall/corrected/check';
my $output_directory = $ENV{HOME} . '/vall/corrected/profile';

my @check_files      = glob("$check_directory/*.check");
#my @check_files = qw|/users/tex/vall/corrected/check/12ASA.check|;

# Read in BLOSUM62 matrix
if (-f "$blosum_matrix" ) {
    open(B62,"$blosum_matrix") or die "couldnt find blosum table\n";
} else {
    die "Can't find blosum table\n";
}

my @b62;
my @blos_aa = (1,15,12,3,2,14,4,6,7,8,10,9,11,5,13,16,17,19,20,18);
my %aaNum = ( A => 1,
	   C => 2,
	   D => 3,
	   E => 4,
	   F => 5,
	   G => 6,
	   H => 7,
	   I => 8,
	   K => 9,
           L => 10,
           M => 11,
           N => 12,
           P => 13,
           Q => 14,
           R => 15,
           S => 16,
           T => 17,
           V => 18,
           W => 19,
           Y => 20
);

my $trash = <B62>;$trash = <B62>;$trash = <B62>;$trash = <B62>;
for (my $i=0;$i<20;$i++) {
    my $line = <B62>;    
    chomp $line;
    my @words  = split /\s+/,$line;
    for (my $j=0;$j<=$i;$j++) {    
        $b62[$blos_aa[$i]][$blos_aa[$j]] = $words[$j];
        $b62[$blos_aa[$j]][$blos_aa[$i]] = $words[$j];
    }
}

close B62 or die $!;

for (my $i=1;$i<=20;$i++) {
    my $sum = 0.0;
    for (my $j=1;$j<=20;$j++) {
        $sum = $sum + $b62[$i][$j];
    }
    for (my $j=1;$j<=20;$j++) {    ## normalize so each row sums to 1.0
        $b62[$i][$j]    = $b62[$i][$j]    / $sum;
    }
} # End BLOSUM62 matrix processing.

# Process all of the checkpoint files.
foreach my $check_file (@check_files) {
    my $filename = basename $check_file;
    my $output_filename;
    if ( $filename =~ /([\d\w\_]+)\.check/ ) {
        $output_filename = $output_directory . '/' . $1 . '.ASCII';
    }

    # Run readblastmatrix to convert the binary file into an ASCII
    # matrix.
    my $ascii_matrix = `echo $check_file | $readblastmatrix`;

    # Clean up the ASCII checkpoint matrix (using code from cleanCheck.pl)
    open CHECKOUT, ">$output_filename" or die "Error opening $output_filename ($!)";
    foreach my $line ( split /\n/, $ascii_matrix ) {
        chomp($line);

        # Skip blank lines.
        if ( $line =~ /^\s*$/ ) { next };

        my @tokens = split /\s+/,$line;    

        # Skip lines that have no tokens
        if ( scalar(@tokens) == 1 ) { next };
	if ( $tokens[0] eq 'X' ) { next };

        if ($tokens[0] =~ /Enter/) {
            print CHECKOUT "$tokens[3]\n";    
        } elsif ($tokens[1] =~ /\d.\d+/ && $line !~ /0\.0000 0\.0000 0\.0000 0\.0000 0\.0000 0\.0000/) {
            print CHECKOUT "$line p\n";
        } elsif ($line =~ /0\.0000 0\.0000 0\.0000 0\.0000 0\.0000 0\.0000/) {
            ## now we take the lines that psi-blast truncates and write out the blosumn line for that aa    
            printf CHECKOUT "%1s ",$tokens[0];
            for (my $i=1;$i<21;$i++) {
                if ( !defined $b62[$aaNum{$tokens[0]}][$i] ) {
                    print "i = $i\n";
                    print "output_filename = $output_filename\n";
                    print "line = $line\n";
                    print "token = ", $tokens[0], "\n";
                }
                printf CHECKOUT "%6.4f ",$b62[$aaNum{$tokens[0]}][$i];
            }
            print CHECKOUT " t\n";
        }
    }
    close CHECKOUT or die $!;
}
