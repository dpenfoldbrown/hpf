#!/usr/bin/perl -w


use strict;
use DBI;

#kdrew: library for parsing command line options
use Getopt::Std;

my $useage = "-t tmpdir (local) -e experiment number -u mysql username -p mysql password\n";

my $file_dir; 
my $experiment_id;
my $mysql_user;
my $mysql_pass;

my %options=();
getopts("t:e:u:p:",\%options);
# like the shell getopt, "d:" means d takes an argument
if ( defined $options{t} ) { $file_dir = $options{t} } else  { $file_dir = "/tmp/" };
if ( defined $options{e} ) { $experiment_id = $options{e} } else  { die "$useage" };
if ( defined $options{u} ) { $mysql_user = $options{u} } else  { die "$useage" };
if ( defined $options{p} ) { $mysql_pass = $options{p} } else  { die "$useage" };
print "Unprocessed by Getopt::Std:\n" if $ARGV[0];


print "extract decoys from bddb\n";

my $dbh = DBI->connect('DBI:mysql:kdrew;host=mysql', $mysql_user, $mysql_pass);

		#create table temp_keys select distinct sequence_key from protein where experiment_key = 84;
		my $create_query = "create temporary table temp_keys select distinct sequence_key from bddb.protein where experiment_key = $experiment_id";

		print "$create_query\n";
		my $example = $dbh->prepare($create_query);
		$example->execute;

		#kdrew: change query below to this:
		#select yk.sequence_key, d.domain_sequence_key, md.outfile_key, mdy.decoy_name from domain as d, yeast_keys as yk, mcmData as md, mcmDecoy as mdy where yk.sequence_key = d.parent_sequence_key and md.sequence_key = d.domain_sequence_key and mdy.outfile_key = md.outfile_key limit 10;

		#my $query = "select distinct md.sequence_key, md.mcm_decoy_key, a2s.nr_ac, a2s.ac2 from mcmData as md, mcmDecoy as mdc, ac2sequence as a2s where a2s.sequence_key = md.sequence_key and mdc.id = md.mcm_decoy_key and a2s.db = \"pdb\" ";

		#my $query = "select distinct d.domain_sequence_key, md.structure_key from bddb.domain as d, temp_keys as tk, bddb.mcmData as md, ddbMeta.structure as strt where tk.sequence_key = d.parent_sequence_key and md.sequence_key = d.domain_sequence_key and md.structure_key = strt.id";
		my $query = "select distinct d.domain_sequence_key, md.cluster_center_index from bddb.domain as d, temp_keys as tk, bddb.mcmData as md, ddbMeta.structure as strt where tk.sequence_key = d.parent_sequence_key and md.outfile_key = d.outfile_key and md.structure_key = strt.id";

		print "$query\n";
		$example = $dbh->prepare($query);
		$example->execute;
		while ((my @results)=$example->fetchrow_array)
		{ 
			print("$results[0]\n"); 

			my $filename = $file_dir.$results[0].".".$results[1].".pdb";
			print "filename: $filename\n";
			open(PDBOUT,">$filename");


			my $query2 = "select uncompress(compress_file_content) from ddbMeta.structure where id = $results[1]";
			print "$query2\n";
			my $example2 = $dbh->prepare($query2);
			$example2->execute;
			while ((my @pdb)=$example2->fetchrow_array)
			{ 
				print PDBOUT "@pdb";
			}
		}

		#`cd $file_dir`;
		#`tar czvf experiment_$experiment_id.tgz *.pdb`;
