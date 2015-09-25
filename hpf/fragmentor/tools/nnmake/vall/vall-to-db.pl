#!/usr/bin/perl -w

use strict;
use File::Basename;

# quick and dirty script for reformatting a Vall .dat file  and putting it in a
# MySQL database. 

$|++;

$0 = basename $0;
my $usage = <<USAGE;
usage: $0 filename tablename
USAGE

# Configuration information
my $database = 'paul2';
my $hostname = 'stool';
my $username = 'murphp';
my $password = 'dhs123';

my $vall_file = $ARGV[0];
my $tablename = $ARGV[1];

die "$usage\n" unless $tablename && -f $vall_file;

print "Creating table $tablename in database ... ";
my $sql = make_create_table_sql( $tablename );
my $create_table_cmd = "mysql -e \'$sql\' -h $hostname --user=$username --password=$password $database";
system($create_table_cmd);
print "done.\n";

print "Reformatting old Vall ... ";
# Reformat the Vall file, put it into an appropriately named temporary file.
open FILE, "<$vall_file" or die "Error opening file $vall_file ($!)";
my @file = <FILE>;
close FILE or die $!;

# Replace spaces with tabs so that mysqlimport is happy.
my @new_file = map { s/\s+/\t/g; $_ } @file;
print "done.\n";

my $outfile = "/tmp/$tablename.txt";
print "Writing Vall to $outfile ... ";
open OUTFILE, ">$outfile" or die $!;
for (@new_file) { print OUTFILE $_, "\n" };
close OUTFILE or die $!;
print "done.\n";

# Run mysqlimport on the newly created file.
my $command = "mysqlimport -h $hostname --local --user=$username --fields-terminated-by='\t' --password=$password $database $outfile";
print "Importing data ...\n";
system($command);
print "done.\n";

# Unlink the old file.
#unlink $outfile or warn $!;

sub make_create_table_sql {
    my $tablename = shift;
    
my $sql = <<SQL;
CREATE TABLE $tablename (
    pdb_chain   CHAR(5),
    ss          CHAR(1),
    begin       INT,
    end         INT,
    CA_x        FLOAT,
    CA_y        FLOAT,
    CA_z        FLOAT,
    phi         FLOAT,
    psi         FLOAT,
    omega       FLOAT,
    chi         FLOAT,
    n_align     INT,
    solv_acc    FLOAT,
    gap_freq    FLOAT,
    profile_A   FLOAT,
    profile_C   FLOAT,
    profile_D   FLOAT,
    profile_E   FLOAT,
    profile_F   FLOAT,
    profile_G   FLOAT,
    profile_H   FLOAT,
    profile_I   FLOAT,
    profile_K   FLOAT,
    profile_L   FLOAT,
    profile_M   FLOAT,
    profile_N   FLOAT,
    profile_P   FLOAT,
    profile_Q   FLOAT,
    profile_R   FLOAT,
    profile_S   FLOAT,
    profile_T   FLOAT,
    profile_V   FLOAT,
    profile_W   FLOAT,
    profile_Y   FLOAT
);
SQL

    return $sql;
}
