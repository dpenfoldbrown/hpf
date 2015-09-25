#!/users/robetta/bin/perl

my $fasta = shift @ARGV;

my $samt99 = "/users/robetta/src/shareware/sam/bin/target99";

if ($fasta =~ /^\S+\.fasta$/) {

	my $chain = $fasta;


	print "Running samt99\n";
	print `$samt99 -seed $fasta -out $chain`; 


}
