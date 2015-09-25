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

$|++;
$0 = basename $0; # avoid printing out the full path to this program

# Set up some (configurable) options.
#my @clusters           = qw/dua/;
my @clusters           = ();
my $setup_cluster      = 0;  # don't set up the cluster filesystem by default
my $make_jobs          = 0;  # don't make job files by default
my $submit_jobs        = 0;  # don't submit jobs by default
my $home               = $ENV{HOME};
my $user               = $ENV{USER};
my $fasta_dir          = $home . '/vall/corrected/fasta';
my $blast_dir          = $home . '/vall/blast_stuff'; #directory with blast binaries and databases
my $max_condor_jobs    = 500;
my $condor_file_prefix = $home . '/vall/condorjob.blast';
my $sleep_seconds      = $ENV{SLEEP} || 300; # seconds to sleep between checking clusters

my $nr_blast_arguments = " -j 2 -h 0.001 -e 0.001 -b 0 -k 0 -d nr ";
my $filtrnr_blast_arguments = " -j 3 -h 0.000001 -e 0.000001 -b 0 -k 0 -d nr ";

my $usage = <<USAGE;
usage:  $0 [options]
    Valid options:
    --setup            setup the cluster filesystem (takes a while)
    --make_jobs        make files containing Condor job listings
    --submit_jobs      submit jobs to Condor clusters that I know about

    --fasta_dir        directory containing .fasta files
    --blast_dir        directory containing BLAST databases and executables

    --max_condor_jobs  maximum number of jobs allowed on a cluster
    --cluster          cluster_name
USAGE

GetOptions( 
    "fasta_dir=s"       => \$fasta_dir,
    "blast_dir=s"       => \$blast_dir,
    "cluster=s"         => \@clusters,
    "setup"             => \$setup_cluster,
    "make_jobs"         => \$make_jobs,
    "submit_jobs"       => \$submit_jobs,
    "max_condor_jobs=i" => \$max_condor_jobs,
    "help"              => sub { warn "$usage\n"; exit 1; },
);

# If we haven't been told to $setup_cluster or $make_jobs, then die with a warning.
if ( !$setup_cluster && !$make_jobs && !$submit_jobs ) {
    warn "Error: Must be told to --make_jobs, --submit_jobs or --setup!\n";
    warn "$usage\n"; 
    exit 1;
}

if ( !@clusters ) {
    print "Error: must define a cluster with --cluster!\n";
}

# Remove duplicates from the list of clusters.
my %hash = map { $_, 1 } @clusters;
@clusters = keys %hash;

# Setup the cluster filesystem if necessary.
if ( $setup_cluster ) {
    foreach my $cluster ( @clusters ) {
        # Transfer FASTA files and PSI-BLAST files over to the cluster.
        print STDERR "Transmitting data files to cluster $cluster (may take a while) ... ";
        system("rsync --compress --ignore-existing -r -e ssh $fasta_dir/ $cluster:$home/fasta" );
        system("rsync --compress --ignore-existing -r -e ssh $blast_dir/ $cluster:$home" );
        print STDERR "done.\n";
        # ssh to the given cluster, set up the proper directory structure: 
        # $ENV{HOME}/fasta       - holds FASTA files
        # $ENV{HOME}/blast/check - holds checkpointing files
        # $ENV{HOME}/blast/hits  - holds blast results
        
        print STDERR "Creating directory structure on cluster $cluster ... ";
        system("ssh $cluster 'rm -rf blast; mkdir blast; cd blast; mkdir check hits filtnr_hits filtnr_check' 2>/dev/null");
        print STDERR "done.\n";

        # remove old log files
        system("ssh $cluster 'rm -f blast.error blast.log blast.output'");
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
    my @fasta_files = map { basename $_ } glob("$fasta_dir/*.fasta");
    
    while ( $job_count < scalar(@fasta_files) - 1 ) {

        my $condor_job = <<TEMPLATE;

universe     = vanilla
Executable   = blast_pgp.pl
Log          = blast.log
Output       = blast.output
Error        = blast.error
Notification = Error
Requirements = (Arch == "INTEL") && (OpSys == "LINUX") && (Disk >= 0) && (Memory >= 250)

TEMPLATE

        # Slice the array up, be sure not to run off the edge of the array.
        my $upper_limit = $job_count + $chunk_size - 1;
        if ( $upper_limit > scalar(@fasta_files) - 1 ) {
            $upper_limit = scalar(@fasta_files) - 1;
        }

        my @local_fasta_files = @fasta_files[ $job_count .. $upper_limit ];

        foreach my $fasta_file ( @local_fasta_files ) {
            if ( $fasta_file =~ /([\w\d_]{5})\.fasta/ ) {
                my $id = $1; 
                my $local_blast_args = $nr_blast_arguments .
                                       " -i $home/fasta/$fasta_file " .
                                       " -o $home/blast/hits/$id.blast " .
                                       " -C $home/blast/check/$id.check ";
        
                #$condor_job .= "Arguments = $local_blast_args \n";
                $condor_job .= "Arguments = $fasta_file\n";
                $condor_job .= "Queue\n";
                $job_count++;
            } else {
                print STDERR "Don't recognize $fasta_file!\n";
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
                    #system("ssh $cluster 'rm $condor_file'");
                    print STDERR "done.\n";
                    print STDERR scalar(@condor_job_files), " job files left.\n";
                }
            }
        }

        #print STDERR "Sleeping $sleep_seconds ...\n";
        # Wait $sleep_seconds until checking again.
        sleep $sleep_seconds;
    }


    print STDERR "Finished submitting all jobs!\n";
}
