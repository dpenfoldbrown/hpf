#!/usr/bin/perl -w

use strict;
use File::Basename;

my $dssp_binary  = $ENV{BAKER_HOME} . '/src/shareware/dssp/dssp';
my $dssp_options = ' -na ';
my $pdb_file_dir = './corrected/pdb';
my $data_dir     = './corrected/dssp';


my @pdb_files = glob("$pdb_file_dir/*.pdb");
foreach my $pdb_file (@pdb_files) {
    my $dssp_output = `$dssp_binary $dssp_options $pdb_file`;
    my @output = split /\n/, $dssp_output;

    # Correct DSSP output (shamelessly stolen from dssp_correct.pl)
    my $dssp_correct_buf = '';
    my $count = 0;

    foreach my $line ( @output ) {
        if ($line=~/^\s\s(...)\s\s(...)(\s.\s[^!].*)/ and  $count==0) {
            $dssp_correct_buf .= "  $1  $2$3\n";
        } elsif ($line=~/^\s\s(...)\s\s(...)(\s.\s[^!].*)/ and $count!=0) {
            
            my $formatted = sprintf("%3d",$count);
            $dssp_correct_buf .= "  $formatted  $2$3\n";
            $count++;
        }
        elsif ($line=~/^\s\s(...)\s\s(...)(\s.\s!.*)/ and $count==0) {
            print STDERR "CHAIN BREAK $line";
            $count=scalar($1);
        } elsif ($line=~/^\s\s(...)\s\s(...)\s.\s!.*/ and $count!=0) {
            print STDERR "CHAIN BREAK $line";
            #Do nothing
        } elsif ($line!~/^\s\s(...)\s\s(...)\s.\s[^!].*/) {
            $dssp_correct_buf .= $line;
        }
    }

    my $pdb_id;
    if ( (basename $pdb_file) =~ /^([\d\w]{5})\.pdb/ ) {
        $pdb_id = $1;
    } else {
        print STDERR "Error: don't recognize PDB ID from $pdb_file!\n";
    }

    open  FILE, ">$data_dir/$pdb_id.dssp" or die $!;
    print FILE $dssp_correct_buf;
    close FILE or die $!;
}
