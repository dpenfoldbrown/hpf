#!/usr/bin/perl

#kdrew: show summary of an organism (by experiment_key) 

use DBI;
use Getopt::Std;

my $useage = "-e experiment -u mysql username -p mysql password \n";

my $mysql_user;
my $mysql_pass;
my $experiment_key;

my %options=();
getopts("e:u:p:",\%options);
# like the shell getopt, "d:" means d takes an argument
if ( defined $options{e} ) { $experiment_key = $options{e} } else  { print "using all experiments\n" };
if ( defined $options{u} ) { $mysql_user = $options{u} } else  { die "$useage" };
if ( defined $options{p} ) { $mysql_pass = $options{p} } else  { die "$useage" };
print "Unprocessed by Getopt::Std:\n" if $ARGV[0];

	$dbh = DBI->connect('dbi:mysql:kdrew;host=mysql',$mysql_user,$mysql_pass);

	$sql = "select name from bddb.experiment as e where e.id = $experiment_key";
	$sth = $dbh->prepare($sql);
	$sth->execute || die "Could not execute SQL statement ... maybe invalid?";
	my ($experiment_name)=$sth->fetchrow_array;

	$sql = "select count(distinct p.sequence_key) from bddb.protein as p where p.experiment_key = $experiment_key";
	$sth = $dbh->prepare($sql);
	$sth->execute || die "Could not execute SQL statement ... maybe invalid?";
	my ($organism_sequence_count)=$sth->fetchrow_array;

	$sql = "select count(distinct p.sequence_key) from bddb.protein as p, ddbCommon.domain as d where p.sequence_key = d.parent_sequence_key and p.experiment_key = $experiment_key";
	$sth = $dbh->prepare($sql);
	$sth->execute || die "Could not execute SQL statement ... maybe invalid?";
	my ($sequence_count_domain)=$sth->fetchrow_array;

	$sql = "select count(distinct d.domain_sequence_key) from bddb.protein as p, ddbCommon.domain as d where p.sequence_key = d.parent_sequence_key and p.experiment_key = $experiment_key";
	$sth = $dbh->prepare($sql);
	$sth->execute || die "Could not execute SQL statement ... maybe invalid?";
	my ($domain_count)=$sth->fetchrow_array;

	#$sql = "select count(distinct d.domain_sequence_key) from bddb.protein as p, ddbCommon.domain as d, bddb.mcmData as md where p.sequence_key = d.parent_sequence_key and d.domain_sequence_key = md.sequence_key and p.experiment_key = $experiment_key";
	$sql = "select count(distinct d.domain_sequence_key) from bddb.protein as p, ddbCommon.domain as d, bddb.mcmData as md where p.sequence_key = d.parent_sequence_key and d.outfile_key = md.outfile_key and p.experiment_key = $experiment_key and d.domain_source = 'ginzu'";
	$sth = $dbh->prepare($sql);
	$sth->execute || die "Could not execute SQL statement ... maybe invalid?";
	my ($mcm_domain_count)=$sth->fetchrow_array;

	#$sql = "select count(distinct d.parent_sequence_key) from bddb.protein as p, bddb.domain as d, bddb.mcmData as md where p.sequence_key = d.parent_sequence_key and d.domain_sequence_key = md.sequence_key and p.experiment_key = $experiment_key";
	$sql = "select count(distinct d.parent_sequence_key) from bddb.protein as p, ddbCommon.domain as d, bddb.mcmData as md where p.sequence_key = d.parent_sequence_key and d.outfile_key = md.outfile_key and p.experiment_key = $experiment_key and d.domain_source = 'ginzu'";
	$sth = $dbh->prepare($sql);
	$sth->execute || die "Could not execute SQL statement ... maybe invalid?";
	my ($mcm_parent_count)=$sth->fetchrow_array;

	print("\n$experiment_key:$experiment_name\n");
	print("total sequence count: $organism_sequence_count\n");
	print("sequence count in domain table: $sequence_count_domain\n");
	print("domain count: $domain_count\n");
	print("domains with mcm scores count: $mcm_domain_count\n");
	print("proteins with mcm scores count: $mcm_parent_count\n");

	$query = "select d.domain_type, count(distinct domain_sequence_key) from bddb.protein as p, ddbCommon.domain as d where p.sequence_key = d.parent_sequence_key and p.experiment_key = $experiment_key group by d.domain_type limit 10";
	$sth = $dbh->prepare($query);
	$sth->execute || die "Could not execute SQL statement ... maybe invalid?";

	print("\nDOMAIN_TYPE:NUM_DOMAINS\n");
	while (my ($domain_type, $num_domain)=$sth->fetchrow_array)
	{ 
		print("$domain_type:$num_domain\n");
	}



