#!/usr/bin/perl -w

# Script for running completePdbCoords.pl and fixMissingDensity.pl on
# a list of PDB files culled from the PISCES server.
# By James Thompson <tex@u.washington.edu> 

use strict;
use File::Basename;
use Data::Dumper::Simple;

$0 = basename $0;

my $usage = <<USAGE;
usage:  $0 pdb_list.txt
USAGE

my $filename = $ARGV[0];
die $usage unless defined $filename && -f $filename;

my $scratchdir = './corrected';
my $ids = get_ids($filename);

# Setup necessary directory structure.
if ( ! -d $scratchdir ) {
    mkdir $scratchdir;
}

foreach my $pdb_id ( @$ids ) {
    if ( ! is_completed( $pdb_id ) ) {
        print "fixing $pdb_id\n";
        fix_pdb_file($pdb_id);
    }
}

# Returns an array reference to a list of PDB ID's that that will be 
# corrected. Takes a single argument that specifies a filename 
# containing PDB IDs, one per line.

sub get_ids {
    my $id_list_file = shift;
    open FILE, "<$id_list_file" or die $!;
    my @lines = <FILE>;
    close FILE or die $!;

    my @id_list;
    # iterate over @lines, grab id's, put them in @id_list
    for (@lines) {
        if ( /^([\d\w_]{5})/ ) {
            push @id_list, $1;
        }
    }

    return \@id_list;
}

# Returns true if there exists a FASTA and PDB file for the given PDB id. Otherwise
# returns false.
sub is_completed {
    my $pdb_id = shift;
    
    if ( -f "$scratchdir/pdb/$pdb_id.pdb" && -f "$scratchdir/fasta/$pdb_id.fasta" ) {
        return 1;
    } else {
        return 0;
    }
}

# Fixes inconsistencies in PDB files by running completePdbCoords.pl and
# removeMissingDensity.pl on a PDB file. Takes a single argument, which is the
# PDB id to be fixed. 

sub fix_pdb_file { 
    my $pdb_chain = shift;

    my $pdb_id = substr($pdb_chain,0,4);
    my $chain  = substr($pdb_chain,4,1);

    my $orig_filename  = &get_pdb_filename($pdb_id);           # original PDB file
    my $int_filename   = "$scratchdir/$pdb_chain.temp";        # intermediate filename
    my $final_filename = "$scratchdir/pdb/$pdb_chain.pdb";     # final, corrected PDB file
    my $fasta_filename = "$scratchdir/fasta/$pdb_chain.fasta"; # fasta sequence file

    if ( ! -d "$scratchdir/pdb" ) {
        mkdir "$scratchdir/pdb";
    }

    if ( ! -d "$scratchdir/fasta" ) {
        mkdir "$scratchdir/fasta";
    }

    # run completePdbCoords.pl
    system("./completePdbCoords.pl -pdbfile $orig_filename -chain $chain -outfile $int_filename -fastaout $fasta_filename ");

    # run removeMissingDensity.pl
    system("./removeMissingDensity.pl -pdbfile $int_filename -outfile $final_filename");

    # get rid of the intermediate filename
    unlink($int_filename);
}

# Takes a PDB id as an argument, returns a path to the PDB file containing 
# that PDB id. 

sub get_pdb_filename {
    my $pdb_id = shift;

    my $rootdir = '/net/pdb';
    my $nextdir = substr($pdb_id,1,2);
    my $pdb_filename = "$rootdir/$nextdir/$pdb_id.pdb";

    return lc $pdb_filename;
}
