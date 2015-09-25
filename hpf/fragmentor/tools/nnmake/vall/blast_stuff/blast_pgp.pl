#!/usr/bin/perl -w

use strict;

my $fasta_file = $ARGV[0];

my $home = "/users/tex/";
my $blast_pgp = './blastpgp';
die "Error: $fasta_file doesn't exist!" unless -f "$home/fasta/$fasta_file";

if ( $fasta_file =~ /([\d\w_]{5})\.fasta/ ) {
    my $id = $1;
    # PSI-BLAST on filtnr to remove low-complexity sequences.
    my $filtr_blast_args = " -i $home/fasta/$fasta_file " .
                           " -o $home/blast/filtnr_hits/$id.filtnr.blast " .
                           #" -C $home/blast/filtnr_check/$id.filtnr.check " .
                           " -C /tmp/$id.filtnr.check " .
                           " -j 3 -h 0.000001 -e 0.000001 -b 0 -k 0 -d $home/filtnr ";

    #print "Running round 1 ... ";
    system("$blast_pgp $filtr_blast_args");
    #print "done.\n";

    my $blast_args = " -i $home/fasta/$fasta_file " .
                     " -C $home/blast/check/$id.check ".
                     #" -R $home/blast/filtnr_check/$id.filtnr.check" .
                     " -R /tmp/$id.filtnr.check" .
                     " -o $home/blast/hits/$id.blast " .
                     #" -j 2 -h 0.01 -e 0.01 -b 0 -k 0 -d $home/nr ";
                     " -t 1 -h 0.01 -e 0.01 -b 0 -k 0 -d $home/nr ";
    #print "Running round 2 ... ";
    system("$blast_pgp $blast_args");
    #print "done.\n";
} else {
    die "Error: don't recognize fasta_file $fasta_file.\n";
}
