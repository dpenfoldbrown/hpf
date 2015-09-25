#!/usr/bin/perl -w

# fragpicker.pl - frontend script to make_fragments.pl with ability to 
# make fragment files of different sizes.
# James Thompson <tex@u.washington.edu>

use lib  '/users/tex/lib/perl5/site_perl';

use strict;
use Getopt::Long;
use File::Basename;
use Parallel::ForkManager;

my $fragdir    = './output';
my $make_frags = './make_fragments.pl';
my $nnmake_dir = './nnmake';

my @fasta_files;
my @sizes;

&GetOptions(
    'fasta=s' => \@fasta_files,
    'size=i'  => \@sizes,
);

$0 = basename $0;
my $usage = <<USAGE
usage: $0 -fasta=<fasta_filename> -size=<integer>
USAGE
;

if ( !@sizes || !@fasta_files ) {
    # put a better usage message in here someday
    warn $usage, "\n";
    exit 1;
}

# handle case when we've been given a glob of fastas
my @entries  = @fasta_files;
@fasta_files = ();
foreach my $entry (@entries) {
    if ( $entry =~ /\*/ ) {
        my @globbed_files = glob($entry); 
        push @fasta_files, @globbed_files;
    } else {
        push @fasta_files, $entry;
    }
}

my $pm = Parallel::ForkManager->new( 5 );
foreach my $fasta (@fasta_files) {
    # Validate the FASTA files given as input by the user
    if ( $fasta !~ /([\d\w_]{4})([\d\w_]{1})\.fasta/ ) {
        die "Error: don't recognize $fasta as a pdb id and chain!\n";
    }
    my $pdb_chain = $1 . $2;
    my $pdbid = $1;
    my $pdb_file = $pdbid . '.pdb';
    
    if ( ! -f $pdb_file ) {
        warn "Can't find $pdb_file, DME values will be meaningless!\n";
    }

    $pm->start and next;
    foreach my $neighbor_count (@sizes) {
        # skip the length if it doesn't look like an integer
        if ( $neighbor_count !~ /^\d+$/ ) {
            print "Error: don't recognize $neighbor_count as an integer!\n";
        }

        print "Making fragments for $fasta with $neighbor_count neighbors.\n";
        # make the appropriate directory structure (if necessary)
        if ( ! -d "$fragdir/$pdb_chain" ) { mkdir "$fragdir/$pdb_chain"; }
        my $working_dir = "$fragdir/$pdb_chain/$neighbor_count/";
        if ( ! -d $working_dir ) { mkdir $working_dir; }
    
        # copy the NNMAKE source code to the proper location
        my $checkout_cmd = "cp $nnmake_dir/* $working_dir/ 2>/dev/null";
        system($checkout_cmd);
    
        # use some trickiness to edit param.h in place and reset max_nn
        my $replace_cmd = "cd $working_dir; perl -p -i -e \'s/max_nn=\\d+/max_nn=$neighbor_count/g\' param.h";
        system($replace_cmd);
    
        # remake NNMAKE binary
        my $make_cmd = "cd $working_dir; make 2>/dev/null;";
        my $output = `$make_cmd`;
    
        # don't forget to copy the fasta file to the correct directory
        my $copy_cmd = "cp $fasta $working_dir; cp -f $pdb_file $working_dir";
        system($copy_cmd);
        
        # make the necessary fragments
        my $make_frags_cmd = "cd $working_dir; $make_frags -id $pdb_chain -nohoms -xx aa $fasta";
	#Took out -nosam!
        print "Executing ", $make_frags_cmd, "\n";
        system($make_frags_cmd);
    }

    $pm->finish;
    print "Finished making fragment files for $fasta with the following sizes:\n";
    print join "\n", @sizes, "\n";
    print '-' x 80, "\n";
}

$pm->wait_all_children;
