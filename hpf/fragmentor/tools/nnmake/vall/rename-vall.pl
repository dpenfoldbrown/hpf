#!/usr/bin/perl -w

use strict;

open FILE, "<$ARGV[0]" or die $!;

while (<FILE>) {
    if ( /^([\d\w_]{4})/ ) {
        my $new = lc $1;
        s/$1/$new/g;
        print $_;
    }
}

close FILE or die $!;
