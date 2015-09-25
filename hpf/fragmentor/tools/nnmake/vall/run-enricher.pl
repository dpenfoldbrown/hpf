#!/usr/bin/perl -w

use strict;
$|++;

use File::Copy;
use File::Basename;
use Parallel::ForkManager;

my $rosetta     = $ENV{HOME} . '/src/bin/rosetta.gcc';
my $enricher    = $ENV{HOME} . '/src/bin/enrich_fragments.py';
my $frag_folder = '.';

my @pdb_list = map { basename $_ } glob("$frag_folder/*");

$0 = basename $0;
my $pm = Parallel::ForkManager->new( 4 );
foreach my $full_pdb_id ( @pdb_list ) {
    if ( ! -d $full_pdb_id ) { next };
    
    my $pdb_id = substr($full_pdb_id,0,4);
    my $chain  = substr($full_pdb_id,4,1);
    warn "Error: pdb_id not defined for $full_pdb_id\n" unless $pdb_id;
    warn "Error: chain not defined for $full_pdb_id\n"  unless $chain;

    my $working_folder = "$frag_folder/$full_pdb_id/";

    # skip this folder if it's been locked (just to be safe)
    my $lockfile = "$working_folder/.$0.lock";
    if ( -f $lockfile ) { 
        next; 
    } else {
        open FILE, ">>$lockfile" or die "Error creating lockfile $lockfile! ($!)\n";
        close FILE or die "Error creating lockfile $lockfile! ($!)\n";
    }

    $pm->start and next;

    my @fragment_files = map { basename $_ } glob("$full_pdb_id/*v1_3");
    print join "\t", @fragment_files, "\n";

    foreach my $ff (@fragment_files) {
        my $cmd   = "cd $working_folder; james_condense_fragment_file.py $ff -TOP_N 999 -STD_FRAG 150 -DIVERSE_FRAG 50";
        print $cmd, "\n";
        system($cmd);
    }

    # don't forget to remove the lockfile
    unlink $lockfile;

    $pm->finish;
}

$pm->wait_all_children;
