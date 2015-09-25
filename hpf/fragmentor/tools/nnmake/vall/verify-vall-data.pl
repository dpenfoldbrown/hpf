#!/usr/bin/perl -w

use strict;
use File::Basename;

$0 = basename $0;
my $usage = <<USAGE;
usage:  $0 pdb_list

USAGE

die $usage unless $ARGV[0] && -f $ARGV[0];

# Configuration Options
my $clean_pdb_dir     = $ENV{HOME} . '/vall/corrected/pdb';
my $idealized_pdb_dir = $ENV{HOME} . '/vall/corrected/ideal';
my $profile_directory = $ENV{HOME} . '/vall/corrected/profile';
my $dssp_dir          = $ENV{HOME} . '/vall/corrected/dssp';
my $fasta_dir         = $ENV{HOME} . '/vall/corrected/fasta';

# Get the list of PDB ids that we'll be working with.
open PDB_LIST, "<$ARGV[0]" or die $!;
my @pdb_list = <PDB_LIST>;
close PDB_LIST or die $!;

@pdb_list = map { chomp $_; $_ } @pdb_list;

# Iterate over the list of PDB files, add appropriate information to 
# the Vall.

my $missed;
foreach my $pdb_id ( @pdb_list ) {
    if ( ! -f "$profile_directory/$pdb_id.ASCII" ) {
        push @{$missed->{profile}}, $pdb_id;
    }
    if ( ! -f "$clean_pdb_dir/$pdb_id.pdb" ) {
        push @{$missed->{clean_pdb}}, $pdb_id;
    }

    if ( ! -f "$dssp_dir/$pdb_id.dssp" ) {
        push @{$missed->{dssp}}, $pdb_id;
    }
    if ( ! -f "$idealized_pdb_dir/$pdb_id"."_0001.pdb" ) {
        push @{$missed->{ideal}}, $pdb_id;
    }
    if ( ! -f "$fasta_dir/$pdb_id.fasta" ) {
        push @{$missed->{fasta}}, $pdb_id;
    }
}

while ( my ($category, $arrayref) = each %$missed ) {
    print "Missing file ", $category, ' (', scalar(@$arrayref), " total)\n";
    print join "\n", @$arrayref;
    print "\n", '-' x 80, "\n";
}
