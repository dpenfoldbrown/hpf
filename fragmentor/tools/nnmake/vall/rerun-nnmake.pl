#!/usr/bin/perl -w

use strict;
$|++;

use File::Copy;
use File::Basename;
use Parallel::ForkManager;

my $nnmake      = $ENV{HOME} . '/src/nnmake/pNNMAKE.gnu';
my $make_frags  = $ENV{HOME} . '/src/nnmake/make_fragments.pl';
my $frag_folder = '.';
my $paths_file  = 'paths_defs.txt';

my @pdb_list = map { basename $_ } glob("$frag_folder/*");

$0 = basename $0;
my $pm = Parallel::ForkManager->new( 4 );
foreach my $full_pdb_id ( @pdb_list ) {
    my $pdb_id = substr($full_pdb_id,0,4);
    my $chain  = substr($full_pdb_id,4,1);
    warn "Error: pdb_id not defined for $full_pdb_id\n" unless $pdb_id;
    warn "Error: chain not defined for $full_pdb_id\n"  unless $chain;

    my $working_folder = "$frag_folder/$full_pdb_id/";
    open PATHS, ">$working_folder/path_defs.txt" or die $!;
    print PATHS make_paths_file();
    close PATHS or die $!;

    # skip this folder if it's been locked (just to be safe)
    my $lockfile = "$working_folder/.$0.lock";
    if ( -f $lockfile ) { 
        next; 
    } else {
        open FILE, ">>$lockfile" or die "Error creating lockfile $lockfile! ($!)\n";
        close FILE or die "Error creating lockfile $lockfile! ($!)\n";
    }

    #my $cmd = "cd $working_folder; $nnmake aa $pdb_id $chain";
    $pm->start and next;
    my $fasta = "$full_pdb_id.fasta";
    my $cmd   = "cd $working_folder; $make_frags -nohoms -id $full_pdb_id -xx aa -nosam $fasta";
    print $cmd, "\n";
    system($cmd);

    # don't forget to remove the lockfile
    unlink $lockfile;

    $pm->finish;
}

$pm->wait_all_children;

sub make_paths_file {
    my $output = <<OUTPUT
./
./
./
/net/local/scratch/nnmake_database/
./
./
./
./
./
./
./
./
/scratch/shared/nnmake_database/
./
./
./
./
vall.dat.2005-11-30
OUTPUT
;

    return $output;
}
