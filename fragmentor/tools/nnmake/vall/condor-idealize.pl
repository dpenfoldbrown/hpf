#!/usr/bin/perl -w

###############################################################################
##
##  Copyright 2005, University of Washington, the Baker Lab, and James Thompson.
##   This document contains private and confidential information and its 
##   disclosure does not constitute publication.  All rights are reserved by 
##   University of Washington, the Baker Lab, and James Thompson, except those 
##   specifically granted by license.
##
##  Initial Author: James Thompson <tex@u.washington.edu>
###############################################################################

use strict;
use Getopt::Long;
use File::Basename;

$0 = basename $0; # avoid printing out the full path to this program

# Set up (configurable) options.
my @clusters           = qw/atum/;
my $setup_cluster      = 0;  # don't set up the cluster filesystem by default
my $make_jobs          = 0;  # don't make job files by default
my $submit_jobs        = 0;  # don't submit jobs by default
my $home               = $ENV{HOME};
my $user               = $ENV{USER};
my $pdb_dir            = $home . '/vall/corrected/pdb';
my $rosetta_dir        = $home . '/vall/rosetta_stuff';; #directory with blast binaries and databases
my $max_condor_jobs    = 500;
my $condor_file_prefix = $home . '/vall/condorjob.idealize';
my $sleep_seconds      = $ENV{SLEEP} || 300; # seconds to sleep between checks

my $usage = <<USAGE;
usage:  $0 [options]
    Valid options:
    --setup            setup the cluster filesystem (takes a while)
    --make_jobs        make files containing Condor job listings
    --submit_jobs      submit jobs to Condor clusters that I know about

    --pdb_dir          directory containing .pdb files
    --rosetta_dir      directory containing Rosetta databases and executables

    --max_condor_jobs  maximum number of jobs allowed on a cluster
    --cluster          name of cluster (can specificy more than once)
USAGE

GetOptions( 
    "pdb_dir=s"         => \$pdb_dir,
    "rosetta_dir=s"     => \$rosetta_dir,
    "cluster=s"         => \@clusters,
    "setup"             => \$setup_cluster,
    "make_jobs"         => \$make_jobs,
    "submit_jobs"       => \$submit_jobs,
    "max_condor_jobs=i" => \$max_condor_jobs,
    "help"              => sub { warn "$usage\n"; exit 1; },
);

# Remove duplicates from the list of clusters.
my %hash = map { $_, 1 } @clusters;
@clusters = keys %hash;

# If we haven't been told to $setup_cluster or $make_jobs, then die with a warning.
if ( !$setup_cluster && !$make_jobs && !$submit_jobs ) {
    warn "Error: Must be told to either setup a cluster or make job files!\n";
    warn "$usage\n"; 
    exit 1;
}

# Setup the cluster filesystem if necessary.
if ( $setup_cluster ) {
    foreach my $cluster (@clusters) {
        # Transfer FASTA files and PSI-BLAST files over to the cluster.
        print STDERR "Transmitting data files to cluster $cluster (may take a while) ... ";
        system("rsync --compress --ignore-existing -r -e ssh $pdb_dir/ $cluster:$home/pdb" );
        system("rsync --compress --ignore-existing -r -e ssh $rosetta_dir/ $cluster:$home" );
        print STDERR "done.\n";
        # remove old log files and old idealized files
        system("ssh $cluster 'rm -f *_0001.pdb; rm -rf idealize.error idealize.log idealize.output'");
    }
}

if ( $make_jobs ) {
    # Create some big Condor job submission files. These can only submit a
    # maxmimum of $max_condor_jobs jobs per cluster. We'll make these job
    # files contain $max_condor_jobs/2 jobs per file.

    # round down $max_condor_jobs if necessary
    if ( ($max_condor_jobs % 2) != 0 ) {
        $max_condor_jobs--;
    }
    
    
    my $job_count = 0;
    my $chunk_size = $max_condor_jobs / 2; # number of jobs per job file
    my @pdb_files = map { basename $_ } glob("$pdb_dir/*.pdb");
    
    while ( $job_count < scalar(@pdb_files) - 1 ) {

        my $condor_job = <<TEMPLATE;

universe     = vanilla
Executable   = ./rosetta++/rosetta.gcc
Log          = idealize.log
Output       = idealize.output
Error        = idealize.error
Notification = Error

TEMPLATE

        # Slice the array up, be sure not to run off the edge of the array.
        my $upper_limit = $job_count + $chunk_size - 1;
        if ( $upper_limit > scalar(@pdb_files) - 1 ) {
            $upper_limit = scalar(@pdb_files) - 1;
        }

        my @local_pdb_files = @pdb_files[ $job_count .. $upper_limit ];

        foreach my $pdb_file ( @local_pdb_files ) {
            if ( $pdb_file =~ /([\w\d]{4})\.pdb/ ) {
                my $id = $1; 
        
                $condor_job .= "Arguments = -idealize -s pdb/$pdb_file\n";
                $condor_job .= "Queue\n";
                $job_count++;
            } else {
                print STDERR "Don't recognize $pdb_file!\n";
            }
        }
        
        # put this job file into a new file, named $condor_file_prefix . $job_count
        open FILE, ">$condor_file_prefix.$job_count" or die $!;
        print FILE $condor_job;
        close FILE or die $!;
    }
}

if ( $submit_jobs ) {
    # get a list of the Condor job files we have.
    my @condor_job_files = map { basename $_ } glob("$condor_file_prefix.*");

    # Check on each cluster, if the number of jobs submitted by each user
    # becomes lower than $max_condor_jobs / 2, then submit another file of
    # jobs.

    while ( scalar (@condor_job_files) > 0 ) {
        foreach my $cluster (@clusters) {

            # Check and see how many jobs I have running on this cluster.
            my $output = `ssh $cluster 'condor_q -sub $user | tail -1' 2>&1`;
            my $number_of_jobs;
            if ( $output =~ /^(\d+) jobs\;/ ) {
                $number_of_jobs = $1;
            } elsif ( $output =~ /Error: Collector has no record of schedd\/submitter/ ) {
                $number_of_jobs = 0;
            }

            # If we are under the limit $max_condor_jobs / 2, submit
            # another job file.
            if ( $number_of_jobs < ($max_condor_jobs / 2) ) {
                if ( my $condor_file = shift @condor_job_files ) {
                    print STDERR "Submitting $condor_file to $cluster ... ";
                    system("rsync -e ssh $condor_file $cluster:$home");
                    system("ssh $cluster 'condor_submit $condor_file'");
                    print STDERR "done.\n";
                    print STDERR scalar(@condor_job_files), " job files left.\n";
                }
            }
        }

        print STDERR "Sleeping $sleep_seconds ...\n";
        # Wait $sleep_seconds until checking again.
        sleep $sleep_seconds;
    }

    print STDERR "Finished submitting all jobs!\n";
}
