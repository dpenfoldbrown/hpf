#!/usr/bin/perl -w

use strict;

my $readmatrix = '/users/tex/src/msaUtil/readblastmatrix';
my $scratch    = '/scratch/tex';
my $psiblast   = '/users/tex/blast/blastpgp';
my $dbname     = '/scratch/shared/genomes/nr';
my $hval       = '0.001';
my $eval       = '0.001';
my $num_passes = 3;

my @pdb_ids = qw/1L2P/;
my $blast_arguments = " -j $num_passes -h $hval -e $eval -b 0 -k 0 -d $dbname ";

foreach my $pdb_id (@pdb_ids) {
    # run PSI-BLAST on each query sequence
    my $cmd;
    $cmd .= "$psiblast $blast_arguments -i corrected/fasta/$pdb_id.fasta ";
    $cmd .= "-o $scratch/$pdb_id.blast -C $scratch/$pdb_id.check -Q $scratch/$pdb_id.matrix";
    print $cmd, "\n";
    system($cmd);
}
