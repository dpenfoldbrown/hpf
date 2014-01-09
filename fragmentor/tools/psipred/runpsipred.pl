#!/users/robetta/bin/perl

# This is a simple script which will carry out all of the basic steps
# required to make a PSIPRED V2 prediction. Note that it assumes that the
# following programs are in the appropriate directories:
# blastpgp - PSIBLAST executable (from NCBI toolkit)
# makemat - IMPALA utility (from NCBI toolkit)
# psipred - PSIPRED V2 program
# psipass2 - PSIPRED V2 program

###############################################################################
# conf
###############################################################################

# The name of the BLAST data bank
#set dbname = allfilt
$dbname = "/scratch/shared/genomes/filtnr";

# Tmp dir for working
$tmpdir = "/scratch/tmp";

# Where the NCBI programs have been installed
#set ncbidir = /usr/local/bin
$ncbidir = "/home/rbonneau/shareware/blast";

# Where the PSIPRED V2 programs have been installed
#set execdir = ./bin
$home = "/home/rbonneau/shareware/psipred2.01";
$execdir = "$home/bin";

# Where the PSIPRED V2 data files have been installed
#set datadir = ./data
$datadir = "$home/data";

###############################################################################
# init
###############################################################################

if ($#ARGV < 0) {
    print STDERR "usage: $0 <fastafile>\n";
    exit -1;
}
$fastafile = shift @ARGV;

$rootname = $basename = $fastafile;
$basename =~ s!^.*\/!!;
$basename =~ s!\.[^\.]*$!!;
$rootname =~ s!\.[^\.]*$!!;
$outname  = "$basename.chk";

system (qq{cp -f $fastafile psitmp.fasta});

###############################################################################
# run
###############################################################################

print "Running PSI-BLAST with sequence $fastafile...\n";

system (qq{$ncbidir/blastpgp -b 0 -j 3 -h 0.001 -d $dbname -i psitmp.fasta -C psitmp.chk -o $tmpdir/$basename.psipred_blast});

print "Predicting secondary structure...\n";

system (qq{echo psitmp.chk > psitmp.pn});
system (qq{echo psitmp.fasta > psitmp.sn});
system (qq{$ncbidir/makemat -P psitmp});

print "Pass1 ...\n";

system (qq{$execdir/psipred psitmp.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 $datadir/weights.dat4 > $basename.psipred_ss});

print "Pass2 ...\n";

system (qq{$execdir/psipass2 $datadir/weights_p2.dat 1 1.0 1.0 $basename.psipred_ss2 $basename.psipred_ss > $basename.psipred_horiz});

# Remove temporary files

print "Cleaning up ...\n";
system (qq{rm -f psitmp.* error.log $tmpdir/$basename.psipred_blast});

print "Final output files: $basename.psipred_ss2 $basename.psipred_horiz\n";
print "Finished.\n";

exit 0;

###############################################################################
# end
###############################################################################
