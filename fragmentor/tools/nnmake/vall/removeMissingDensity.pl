#!/usr/bin/perl -w
##
##
## Copyright 2000, University of Washington
##   This document contains private and confidential information and
##   its disclosure does not constitute publication.  All rights are
##   reserved by University of Washington, except those specifically 
##   granted by license.
##
##
##  Initial Author: Dylan Chivian (dylan@lazy8.com)
##  $Revision: 6785 $
##  $Date: 2005-10-07 17:35:49 -0400 (Fri, 07 Oct 2005) $
##  $Author: tex $
##
##
###############################################################################

###############################################################################
package PDButil;
###############################################################################

$| = 1;                                                   # don't buffer stdout

local %opts = &getCommandLineOptions ();
local $pdbfile = $opts{pdbfile};
local $outfile = $opts{outfile};

$pdbID = $pdbfile;
$pdbID =~ s!^.*/?p?d?b?(\w\w\w\w)\.[pe][dn][bt]\.?g?z?Z?$!$1!;
$pdbID = lc $pdbID;

$pdb_basefile = $pdbfile;
$pdb_basefile =~ s!^.*\/!!;

###############################################################################
# main
###############################################################################

# read
#
#@pdbbuf = &fileBufArray ($pdbfile);
#foreach $line (@pdbbuf) {
#    if ($line =~ /^ATOM/) {
#	$occ = substr ($line, 54, 6);
#	if ($occ > 0.0) {
#	    push (@out, $line);
#	}
#    } elsif ($line !~ /^HETATM/) {
#	push (@out, $line);
#    }
#}

@pdbbuf = &fileBufArray ($pdbfile);
$last_chain  = '';
$last_resseq = undef;
$start_discard_resseq = 1;
$stop_discard_resseq  = undef;
@aa_buf = ();
@out = ();
$bb_defined = undef;
for ($i=0; $i <= $#pdbbuf; ++$i) {
    $line = $pdbbuf[$i];
    if ($line =~ /^ATOM/) {

	$occ = substr ($line, 54, 6);
	next if ($occ <= 0.0);
        $resseq = substr ($line, 22, 4);
#        $x      = substr ($line, 30, 8);
        $y      = substr ($line, 38, 8);
        $z      = substr ($line, 46, 8);
#        next if ($x == (-999.000+$resseq) && $y == (-999.000+$resseq) && $z == (-999.000+$resseq));
#        next if ($x == (-999.000+$resseq) && $y == -999.000 && $z == -999.000);
#        next if ($x == -999.000 && $y == -999.000 && $z == -999.000);
#        next if ($x == (-900.000+5*$resseq) && $y == -900.000 && $z == -900.000);
        next if ($y == -999.999 && $z == -999.999);
        next if ($y == -900.000 && $z == -900.000);
	
	$chain = substr ($line, 21, 1);
	if ($chain ne $last_chain) {
	    $last_chain = $chain;
	    $last_resseq = undef;
	    @aa_buf = ();
	    $N_occ      = undef;
	    $CA_occ     = undef;
	    $C_occ      = undef;
	    $O_occ      = undef;
	    $bb_defined = undef;
	}

	if (! defined $last_resseq || $resseq != $last_resseq || $i == $#pdbbuf) {
	    $last_resseq = $resseq  if (! defined $last_resseq);
	    if ($i == $#pdbbuf) {
		push (@aa_buf, $line);
		$atomname = substr ($line, 12, 4);
		if ($atomname eq ' N  ') {
		    $N_occ  = 1;
		} elsif ($atomname eq ' CA ') {
		    $CA_occ = 1;
		} elsif ($atomname eq ' C  ') {
		    $C_occ  = 1;
		} elsif ($atomname eq ' O  ') {
		    $O_occ  = 1;
		}
		$bb_defined = 'true'  if ($N_occ && $CA_occ && $C_occ && $O_occ);
	    }
	    if ($bb_defined) {
		$bb_defined_somewhere = 'true';
		push (@out, @aa_buf);
		$start_discard_resseq = $resseq;        # position lags density
	    } else {
		$stop_discard_resseq = ($i == $#pdbbuf) ? $resseq : $resseq-1; # pos lags density
		if ($resseq != 1) {
		    printf ("DISCARDING from $pdb_basefile chain: '%1s', resseq_range: %4d-%-4d\n", $chain, $start_discard_resseq, $stop_discard_resseq);
		}
	    }
	    $last_resseq = $resseq;
	    @aa_buf = ();
	    $N_occ      = undef;
	    $CA_occ     = undef;
	    $C_occ      = undef;
	    $O_occ      = undef;
	    $bb_defined = undef;
	}
	push (@aa_buf, $line);
	$atomname = substr ($line, 12, 4);
	if ($atomname eq ' N  ') {
	    $N_occ  = 1;
	} elsif ($atomname eq ' CA ') {
	    $CA_occ = 1;
	} elsif ($atomname eq ' C  ') {
	    $C_occ  = 1;
	} elsif ($atomname eq ' O  ') {
	    $O_occ  = 1;
	}
	$bb_defined = 'true'  if ($N_occ && $CA_occ && $C_occ && $O_occ);
    }
    elsif ($line !~ /^HETATM/) {
	if (@aa_buf) {
	    if ($bb_defined) {
		push (@out, @aa_buf);
		$bb_defined_somewhere = 'true';
	    } else {
		$stop_discard_resseq = $resseq; # pos lags density
		printf ("DISCARDING from $pdb_basefile chain: '%1s', resseq_range: %4d-%-4d\n", $chain, $start_discard_resseq, $stop_discard_resseq);
	    }
	    $last_chain  = undef;
	    $last_resseq = undef;
	    @aa_buf = ();
	    $N_occ      = undef;
	    $CA_occ     = undef;
	    $C_occ      = undef;
	    $O_occ      = undef;
	    $bb_defined = undef;
	}
	push (@out, $line);
    }
}
if (! $bb_defined_somewhere) {
    print STDERR ("no backbone in pdb file '$pdbfile'\n");
    exit 0;
}


# output
#
$outbuf = join ("\n", @out)."\n"  if (@out);
if ($outfile) {
#    print "creating $outfile\n";
    open (OUTFILE, '>'.$outfile);
    select (OUTFILE);
}
print $outbuf  if ($outbuf);
if ($outfile) {
    close (OUTFILE);
    select (STDOUT);
}
&abort ("failed to create out file '$outfile'") if ($outfile && ! -s $outfile);


# exit
#
exit 0;

###############################################################################
# subs
###############################################################################

# getCommandLineOptions()
#
#  desc: get the command line options
#
#  args: none
#
#  rets: \%opts  pointer to hash of kv pairs of command line options
#
sub getCommandLineOptions {
    use Getopt::Long;
    local $usage = qq{
usage: $0 
\t-pdbfile          <pdbfile>
\t[-outfile         <outfile>]
};

    # Get args
    #
    local %opts = ();
    &GetOptions (\%opts, "pdbfile=s", "outfile=s");


    # Check for legal invocation
    #
    if (! defined $opts{pdbfile}) {
        print STDERR "$usage\n";
        exit -1;
    }
    &checkExistence ('f', $opts{pdbfile});	

    return %opts;
}

###############################################################################
# util
###############################################################################

sub maxInt {
    local ($v1, $v2) = @_;
    return ($v1 > $v2) ? $v1 : $v2;
}

sub tidyDecimals {
    my ($num, $decimal_places) = @_;
    if ($num !~ /\./) {
	$num .= '.' . '0' x $decimal_places;
	$num =~ s/^0+//;
    }
    else {
	if ($num =~ s/(.*\.\d{$decimal_places})(\d).*$/$1/) {
	    my $nextbit = $2;
	    if ($nextbit >= 5) {
		my $flip = '0.' . '0' x ($decimal_places - 1) . '1'; 
		$num += $flip;
	    }
        }
	$num =~ s/^0//;
	my $extra_places = ($decimal_places + 1) - length $num;
	$num .= '0' x $extra_places  if ($extra_places > 0);
    }

    return $num;
}

sub distsq {
    local @dims = @_;
    local $v = 0;
    foreach $dim (@dims) {
	$v += $dim*$dim;
    }
    return $v;
}

sub logMsg {
    local ($msg, $logfile) = @_;

    if ($logfile) {
        open (LOGFILE, ">".$logfile);
        select (LOGFILE);
    }
    else {
	select (STDERR);
    }
    print $msg, "\n";
    if ($logfile) {
        close (LOGFILE);
    }
    select (STDOUT);

    return 'true';
}

sub checkExistence {
    local ($type, $path) = @_;
    if ($type eq 'd') {
	if (! -d $path) { 
            print STDERR "$0: dirnotfound: $path\n";
            exit -3;
	}
    }
    elsif ($type eq 'f') {
	if (! -f $path) {
            print STDERR "$0: filenotfound: $path\n";
            exit -3;
	}
    }
}

sub abort {
    my $msg = shift;
    my $date = `date +'%Y-%m-%d_%T'`;  chomp $date;
    print STDERR "[$date]:$0:ABORT: $msg\n";
    exit -2;
}

sub writeBufToFile {
    ($file, $bufptr) = @_;
    if (! open (FILE, '>'.$file)) {
	&abort ("$0: unable to open file $file for writing");
    }
    print FILE join ("\n", @{$bufptr}), "\n";
    close (FILE);
    return;
}

sub fileBufString {
    local $file = shift;
    local $oldsep = $/;
    undef $/;
    if ($file =~ /\.gz|\.Z/) {
	if (! open (FILE, "gzip -dc $file |")) {
	    &abort ("$0: unable to open file $file for gzip -dc");
	}
    }
    elsif (! open (FILE, $file)) {
	&abort ("$0: unable to open file $file for reading");
    }
    local $buf = <FILE>;
    close (FILE);
    $/ = $oldsep;
    return $buf;
}

sub fileBufArray {
    local $file = shift;
    local $oldsep = $/;
    undef $/;
    if ($file =~ /\.gz|\.Z/) {
	if (! open (FILE, "gzip -dc $file |")) {
	    &abort ("$0: unable to open file $file for gzip -dc");
	}
    }
    elsif (! open (FILE, $file)) {
	&abort ("$0: unable to open file $file for reading");
    }
    local $buf = <FILE>;
    close (FILE);
    $/ = $oldsep;
    @buf = split (/$oldsep/, $buf);
    pop (@buf)  if ($buf[$#buf] eq '');
    return @buf;
}

###############################################################################
1; # package end
# end
###############################################################################
