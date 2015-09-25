#!/users/robetta/bin/perl -w
# Condense 6-state sam-ss prediction into 3-state
use strict;
die "USAGE: $0 <file>\n" if $#ARGV != 0; # One argument; rdb6 file
my $param;
$param->{'ext'} = '.rdb'; # Extension of file to output
my $file = shift;
my $rows;
open IN, "<$file" || die "Cannot open file $file: $!\n";
@{$rows} = <IN>;
close (IN);
$rows = &condense_rdb6_rdb($rows);
$file =~ /^(.*\/?[^\/]+?)\.[^\.]+$/;
open OUT, ">".$1.$param->{'ext'} || die "Cannot open file ".$1.$param->{'ext'}.": $!\n";
print OUT @{$rows};
close (OUT);

# *** Subroutines ***
sub condense_rdb6_rdb {
	my($rows)=@_;
	my($outrows,@parts,$line);
	for (@{$rows}) {
		if ( /^#/) {
			push (@{$outrows}, $_);
		} elsif (/^Pos/) {
			push (@{$outrows}, "Pos\tAA\tE\tH\tL\n");
		} elsif (/^10N/) {
			push (@{$outrows}, "10N\t1S\t5N\t5N\t5N\n");
		} elsif (/^[0-9]/) {
			@parts = split;
			die "Wrong number of columns in results row: $_\n" if $#parts != 7;
			$line = sprintf "%8s\t%1s\t%8s\t%8s\t%8s\n",$parts[0]-1,$parts[1],$parts[2]+$parts[3],$parts[4]+$parts[5],$parts[6]+$parts[7];
			push (@{$outrows}, $line);
		} else {
			die "Unrecognized row: $_\n";
		}
	}
	return $outrows;
}
