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

$debug = 1;

###############################################################################
# init
###############################################################################

$| = 1;                                                   # don't buffer stdout

local %opts = &getCommandLineOptions ();
local $pdbfile          = $opts{pdbfile};
local $chainOI          = $opts{chain};
local $outfile          = $opts{outfile};
local $fastaout         = $opts{fastaout};
local $fastain          = $opts{fastain};
local $resmapout        = $opts{resmapout};

$pdbID = $pdbfile;
$pdbID =~ s!^.*/!!;
$pdbID =~ s!^pdb|\.[pe][dn][bt]\.?g?z?Z?$!!g;
$pdbID =~ s!^(....)(.?)(.?)(.*)!(lc $1).(uc $2).(uc $3).$4!e;

if ($outfile) {
    $outdir = $outfile;
    $outdir =~ s!^(.*)/[^/]*!$1!;
    $outdir = '.' if ($outdir eq $outfile);
    system (qq{mkdir -p $outdir})  if (! -d $outdir);
}
if ($fastaout) {
    $fastaoutdir = $fastaout;
    $fastaoutdir =~ s!^(.*)/[^/]*!$1!;
    $fastaoutdir = '.' if ($fastaoutdir eq $fastaout);
    system (qq{mkdir -p $fastaoutdir})  if (! -d $fastaoutdir);
}
if ($resmapout) {
    $resmapoutdir = $resmapout;
    $resmapoutdir =~ s!^(.*)/[^/]*!$1!;
    $resmapoutdir = '.' if ($resmapoutdir eq $resmapout);
    system (qq{mkdir -p $resmapoutdir})  if (! -d $resmapoutdir);
}

###############################################################################
# main
###############################################################################

# read
#
($pdb, $chainOI) = &readPdbFile ($pdbfile, $chainOI, $fastain);
&abort ("no chain $chainOI for pdb $pdbID") if (!$pdb->{bodybufs}->{"$chainOI"});


# write fasta
#
if ($fastaout) {
    $fasta_buf = ">$pdbID$chainOI  ".($#{$pdb->{seqres}->{"$chainOI"}}+1)." res\n";
    for (my $i=0; $i <= $#{$pdb->{seqres}->{"$chainOI"}}; ++$i) {
	$fasta_buf .= $pdb->{seqres}->{"$chainOI"}->[$i];
	$fasta_buf .= "\n"  if (($i+1) % 50 == 0 && $i != $#{$pdb->{seqres}->{"$chainOI"}});
    }
    $fasta_buf .= "\n";
    open  (FASTA, '>'.$fastaout);
    print  FASTA $fasta_buf;
    close (FASTA);
}


# write resmap
#
if ($resmapout) {
    $resmap_buf = '';
    for (my $i=0; $i <= $#{$pdb->{seqres}->{"$chainOI"}}; ++$i) {
	$coordseq_n = $pdb->{resmap}->{"$chainOI"}->[$i];
	$coordseq_n = -9999  if (! defined $coordseq_n);
	$resmap_buf .= sprintf ("%1s %5d %5s\n", $pdb->{seqres}->{"$chainOI"}->[$i], $i + 1, $coordseq_n);
    }
    open  (RESMAP, '>'.$resmapout);
    print  RESMAP $resmap_buf;
    close (RESMAP);
}


#
# build pdb output buffer
#

# get body
$res_buf = &residueRecs ($pdb, $chainOI);

# HEADER
foreach $line (@{$pdb->{headerbuf}}) {
    $outbuf .= "$line\n";
}

# BODY
my @body = ();
push (@body, split (/\n/, $res_buf));
push (@body, "TER    atmi");

# HETATMs
foreach $line (@{$pdb->{hetatmbuf}}) {
    #push (@body, $line)  if ($line !~ /HOH/);
    push (@body, $line);
}

# finish
for ($i=0; $i <= $#body; ++$i) {
    substr ($body[$i], 6, 5) = sprintf ("%5d", $i+1);
}
push (@body, "END");
$outbuf .= join ("\n", @body)."\n";


# output
#

if ($outfile) {
    if ($outfile =~ s/\.gz$//) {
	$gzip_flag = 'true';
    }
    elsif ($outfile =~ s/\.Z$//) { 
	$compress_flag = 'true';
    }

#    print "creating $outfile\n";
    open (OUTFILE, '>'.$outfile);
    select (OUTFILE);
}
print $outbuf;
if ($outfile) {
    close (OUTFILE);
    select (STDOUT);

    if ($gzip_flag) {
	system (qq{gzip -9f $outfile});
    }
    elsif ($compress_flag) {
	system (qq{compress -f $outfile});
    }
}


# exit
#
exit 0;

###############################################################################
# subs
###############################################################################

# residue rec string
#
sub residueRecs {
    my ($pdb, $chainOI) = @_;
    my $outbuf = '';
    my $res = $pdb->{res};
    my $atom_type = '';
    my $res_i = 0, $line_i = 0;
    my @atom_types = qw (N CA C O CB 
			 CG OG OG1 CG1 CG2 SG
			 CD OD1 ND1 ND2 CD1 CD2 OD2 SD
			 CE NE1 CE1 CE2 CE3 OE1 OE2 NE NE2 
			 CZ CZ2 CZ3 NZ 
			 NH1 NH2 CH2 OH
			);

    for ($res_i=0, $line_i=0; $res_i <= $#{$res->{"$chainOI"}}; ++$res_i) {
	foreach $atom_type (@atom_types) {
	    next if (! defined $res->{"$chainOI"}->[$res_i]->{atoms}->{$atom_type}->{rec});
	    $rec = $res->{"$chainOI"}->[$res_i]->{atoms}->{$atom_type}->{rec};
	    substr ($rec, 6, 5) = sprintf ("%5d", ++$line_i); 
	    $outbuf .= "$rec\n";
	}
    }
    return $outbuf;
}


# readPdbFile()
#
sub readPdbFile {
    ($pdbfile, $chainOI, $fastain) = @_; # note: must be able to update chainOI
    my $pdb = +{};
    my $headerbuf = +[];                                # $headerbuf->[$line_i]
    my $bodybufs = +{};                      # $bodybufs->{"$chain"}->[$line_i]
    my $hetatmbuf = +[];                                # $hetatmbuf->[$line_i]
    my $header = +{};                                       # $header->{$field}
    my $atoms = +{};                   # $atoms->{"$chain"}->[$res_i]->{$field}
    my $seqres = +{};                           # $seqres->{"$chain"}->[$res_i]
    my $coordseq = +{};                       # $coordseq->{"$chain"}->[$res_i]
    my $resmap = +{};                           # $resmap->{"$chain"}->[$res_i]
    my $res = +{};     # $res->{"$chain"}->[$res_i]->{atoms}->{$atom}->{$field}
    my $modres = +{};
    my $modres_reverse = +{};
    my $last_ResSeq_iCode = +{};
    my $last_ResName = +{};
    my $coordseq_i = +{};
    my $altLoc_used = +{};
    my $atomrec = +{};

    my @pdbbuf = &fileBufArray ($pdbfile);


    # split file into header and body
    #
    local $header_done = undef;
    for (my $i=0; $i <= $#pdbbuf && ! $header_done; ++$i) {
	$recType = substr ($pdbbuf[$i], 0, 6);
	if ($recType !~ /^ATOM|^HETATM/) {
	    if ($recType !~ /^HELIX|^SHEET/) {
		push (@$headerbuf, $pdbbuf[$i]);
		$rec_name = substr ($pdbbuf[$i], 0, 6);
		push (@{$header->{$rec_name}}, substr ($pdbbuf[$i], 6));
	    }
	}
	else {
	    $header_done = 'true';
	    $body_start = $i;
	}
    }


    # build modres
    #
    if ($header->{MODRES}) {
	foreach $line (@{$header->{MODRES}}) {
	    my $deriv = substr ($line, 6, 3);
	    my $src = substr ($line, 18, 3);
	    if ($src !~ /^\s*$/ && $deriv !~ /^\s*$/) {
#		if (&mapResCode($deriv) eq 'Z') {
		if (&mapResCode($deriv) eq 'X') {              # haven't seen it
		    $modres->{$deriv} = $src;
		    $modres_reverse->{$src} = $deriv;
		}
	    }
	}
    }


    # build seqres
    #
    if (! $fastain) {
	if (! $header->{SEQRES}) {
	    &abort ("no SEQRES found in header for pdb $pdbID");
	}
	foreach $line (@{$header->{SEQRES}}) {
	    my $chain = substr ($line, 5, 1);
	    $chain = '_'  if ($chain eq ' ');
	    my @seqres3 = split (/\s+/, substr ($line, 13, 51));
	    shift (@seqres3) if (! $seqres3[0]);
	    pop (@seqres3)  if (! $seqres3[$#seqres3]);
	    for (my $i=0; $i <= $#seqres3; ++$i) {
		if(length $seqres3[$i] != 3) {                        # DNA/RNA
		    $seqres3[$i] = (' 'x(3 - length $seqres3[$i])).$seqres3[$i];
		}
		next if ($seqres3[$i] eq 'EXC');
		if ($modres->{$seqres3[$i]}) {
		    push (@{$seqres->{"$chain"}}, &mapResCode($modres->{$seqres3[$i]}));
		} else {
		    push (@{$seqres->{"$chain"}}, &mapResCode($seqres3[$i]));
		}
	    }
	}

	# handle 1xob (SEQRES: chain '_', COORDS: chain 'A')
	if ($pdbID =~ /^1xob/) {
	    $seqres->{'A'} = $seqres->{'_'};
	    $seqres->{'_'} = undef;
	}

	# change $chainOI from '_' to 'A' if that's what there is
	if ($chainOI eq '_' && defined $seqres->{'A'}) {
	    print STDERR "$0: WARNING: changing chain from _ to A for pdbfile $pdbfile\n";
	    $chainOI = 'A';
	}
	# change $chainOI from 'A' to '_' if that's what there is
	if ($chainOI eq 'A' && defined $seqres->{'_'}) {
	    print STDERR "$0: WARNING: changing chain from A to _ for pdbfile $pdbfile\n";
	    $chainOI = '_';
	}


	# handle obnoxious exceptions where seqres differs from coords
	#  trust coords
	#  (BAD xtalographers, BAD!)
	if ($pdbID =~ /^1a6y/) {
	    $seqres->{'A'}->[24] = 'L';
	    $seqres->{'B'}->[24] = 'L';
	}
	elsif ($pdbID =~ /^1e13/) {              
	    $seqres->{'A'}->[0] = 'H'; 
	    $seqres->{'A'}->[1] = 'M'; 
	    $seqres->{'B'}->[0] = 'H'; 
	    $seqres->{'B'}->[1] = 'M'; 
	}
	elsif ($pdbID =~ /^1fna/) {
	    $seqres->{'_'}->[0] = 'G'; 
	    $seqres->{'_'}->[1] = 'S'; 
	}
	elsif ($pdbID =~ /^1hra/) {
	    $seqres->{'_'}->[0] = 'M';
	}
	elsif ($pdbID =~ /^1l3a/) {
	    $seqres->{'A'}->[111] = 'T';
	    $seqres->{'B'}->[111] = 'T';
	    $seqres->{'C'}->[111] = 'T';
	    $seqres->{'D'}->[111] = 'T';
	}
	elsif ($pdbID =~ /^1mcb/ ||
	       $pdbID =~ /^1mcc/ ||
	       $pdbID =~ /^1mcd/ ||
	       $pdbID =~ /^1mce/ ||
	       $pdbID =~ /^1mcf/ ||
	       $pdbID =~ /^1mch/ ||
	       $pdbID =~ /^1mci/ ||
	       $pdbID =~ /^1mcj/ ||
	       $pdbID =~ /^1mck/ ||
	       $pdbID =~ /^1mcl/ ||
	       $pdbID =~ /^1mcn/ ||
	       $pdbID =~ /^1mcq/ ||
	       $pdbID =~ /^1mcr/ ||
	       $pdbID =~ /^1mcs/ ) {                 
	    $seqres->{'A'}->[0] = 'P';
	    $seqres->{'B'}->[0] = 'P';
	}
	elsif ($pdbID =~ /^1mio/) {
	    $seqres->{'A'}->[424] = 'A';
	}
	elsif ($pdbID =~ /^1pda/) {
	    $seqres->{'_'}->[240] = 'G';
	    $seqres->{'_'}->[260] = 'A';
	}
	elsif ($pdbID =~ /^1s01/) {
	    $seqres->{'_'}->[55]  = 'N';
	    $seqres->{'_'}->[56]  = 'P';
	    $seqres->{'_'}->[60]  = 'N';
	    $seqres->{'_'}->[87]  = 'A';
	    $seqres->{'_'}->[88]  = 'S';
	    $seqres->{'_'}->[97]  = 'A';
	    $seqres->{'_'}->[98]  = 'D';
	    $seqres->{'_'}->[157] = 'T';
	    $seqres->{'_'}->[158] = 'S';
	    $seqres->{'_'}->[168] = 'A';
	    $seqres->{'_'}->[181] = 'S';
	    $seqres->{'_'}->[250] = 'E';
	}
	elsif ($pdbID =~ /^1sos/) {
	    $seqres->{'A'}->[0] = 'A';
	    $seqres->{'B'}->[0] = 'A';
	    $seqres->{'C'}->[0] = 'A';
	    $seqres->{'D'}->[0] = 'A';
	    $seqres->{'E'}->[0] = 'A';
	    $seqres->{'F'}->[0] = 'A';
	    $seqres->{'G'}->[0] = 'A';
	    $seqres->{'H'}->[0] = 'A';
	    $seqres->{'I'}->[0] = 'A';
	    $seqres->{'J'}->[0] = 'A';
	}
	elsif ($pdbID =~ /^1tbp/) {
	    $seqres->{'A'}->[0] = 'M';
	    $seqres->{'B'}->[0] = 'M';
	}
	elsif ($pdbID =~ /^1tnf/) {
	    $seqres->{'A'}->[142] = 'L'; 
	    $seqres->{'B'}->[142] = 'L'; 
	    $seqres->{'C'}->[142] = 'L'; 
	}
	elsif ($pdbID =~ /^1tx4/) {
	    $seqres->{'B'}->[48] = 'K'; 
	}
	elsif ($pdbID =~ /^1ysa/) {
	    $seqres->{'C'}->[0] = 'K';
	    $seqres->{'D'}->[0] = 'M';
	    $seqres->{'D'}->[1] = 'K';
	}
	elsif ($pdbID =~ /^2fam/) {
	    $seqres->{'_'}->[146] = 'A';
	}
	elsif ($pdbID =~ /^2bce/) {
	    $seqres->{'_'}->[186] = 'Q';
	    $seqres->{'_'}->[360] = 'Q';
	}
	elsif ($pdbID =~ /^3hvt/) {
	    $seqres->{'B'}->[276] = 'R';
	}
    } 
    else {                                        # just read the fasta in file
	@fasta_buf = &fileBufArray ($fastain);
	foreach $line (@fasta_buf) {
	    next if ($line =~ /^\>/);
	    $line =~ s/\s+//g;
	    push (@{$seqres->{"$chainOI"}}, split (//, $line));
	}
    }
    &abort ("SEQRES header does not have records for pdb '$pdbID' chain '$chainOI'")  if (! $seqres->{"$chainOI"});
    

    # read body
    #
    for (my $i=$body_start; $i <= $#pdbbuf; ++$i) {
	last if ($pdbbuf[$i] =~ /^END/);
	my $recType = substr ($pdbbuf[$i], 0, 6);
	if ($recType =~ /^ATOM|^HETATM/) {
	    if ($recType =~ /^HETATM/) {
		my $resName = substr ($pdbbuf[$i], 17, 3);
#		next if ($resName eq 'HOH');                       # unecessary
#		next if (! $modres->{$resName} && (&mapResCode($resName) eq 'X' || &mapResCode($resName) eq 'Z'));  # wrong
		next if (! $modres->{$resName} && ! &common_modres_map($resName));
	    }

	    # parse atom
	    #
	    $atomrec = &parseAtomRec ($pdbbuf[$i], $modres);

	    next if ($atomrec->{resName} eq 'X' || $atomrec->{resName} eq 'Z');

	    # get chain
	    #
#	    $chain = substr ($pdbbuf[$i], 21, 1);
#	    $chain = '_'  if ($chain eq ' ');
	    $chain = $atomrec->{chainID};

	    # skip chains we're not interested in
	    next if ($chain ne $chainOI);

	    # skip all ATOM records for a chain if TER has been found
	    next  if ($ter_found{"$chain"} && ($recType =~ /^ATOM/ || $recType =~ /^HETATM/));

	    # skip all ATOM and HETATM records for a chain if ENDMDL has been found
	    next  if ($endmdl_found{"$chain"} && ($recType =~ /^ATOM/ || $recType =~ /^HETATM/));


	    # skip problem coordinates (from annoying pdb files)
	    #
	    if ($pdbID =~ /^1ak9/) {
		if ($atomrec->{chainID} eq '_') {
		    if (($atomrec->{resSeq} == 257 && $atomrec->{altLoc} eq 'B') ||
			($atomrec->{resSeq} == 258 && $atomrec->{altLoc} eq 'B')) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1b3o/) {
		if ($atomrec->{chainID} eq 'A' ||
		    $atomrec->{chainID} eq 'B') {
		    if ($atomrec->{resSeq} == 499) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1b89/) {                    
		if ($atomrec->{chainID} eq 'A') {
		    if ($atomrec->{resSeq} >=  1182 && $atomrec->{resSeq} <= 1194) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1cbn/) {
		if ($atomrec->{chainID} eq '_') {
		    if (($atomrec->{resSeq} == 22 && $atomrec->{iCode} eq 'B') ||
			($atomrec->{resSeq} == 25 && $atomrec->{iCode} eq 'B')) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1d2q/) {
		if ($atomrec->{chainID} eq 'A') {
		    if (($atomrec->{resSeq} == 134 && $atomrec->{resName} eq 'N') ||
			($atomrec->{resSeq} == 139 && $atomrec->{resName} eq 'A') ||
			($atomrec->{resSeq} == 140 && $atomrec->{resName} eq 'A') ||
			($atomrec->{resSeq} == 263 && $atomrec->{resName} eq 'A') ||
			($atomrec->{resSeq} == 264 && $atomrec->{resName} eq 'H') ||
			($atomrec->{resSeq} == 265 && $atomrec->{resName} eq 'A')) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1e3h/) {                    
		if ($atomrec->{chainID} eq 'A') {
#		    if ($atomrec->{resSeq} >=  623 && $atomrec->{resSeq} <= 632) {
		    if ($atomrec->{resSeq} >=  605) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1eis/) {                    
		if ($atomrec->{chainID} eq 'A') {
		    if ($atomrec->{resSeq} == 10 && $atomrec->{resName} eq 'G') {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1ejg/) {
		if ($atomrec->{chainID} eq 'A') {
		    if (($atomrec->{resSeq} == 22 && $atomrec->{resName} eq 'S') ||
			($atomrec->{resSeq} == 25 && $atomrec->{resName} eq 'I')) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1en2/) {
		if ($atomrec->{chainID} eq 'A') {
		    if (($atomrec->{resSeq} == 10 && $atomrec->{resName} eq 'S') ||
			($atomrec->{resSeq} == 14 && $atomrec->{resName} eq 'G') ||
			($atomrec->{resSeq} == 16 && $atomrec->{resName} eq 'R') ||
			($atomrec->{resSeq} == 80 && $atomrec->{resName} eq 'G') ||
			($atomrec->{resSeq} == 81 && $atomrec->{resName} eq 'N')) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1enm/) {
		if ($atomrec->{chainID} eq 'A') {
		    if (($atomrec->{resSeq} == 10 && $atomrec->{resName} eq 'S') ||
			($atomrec->{resSeq} == 14 && $atomrec->{resName} eq 'G') ||
			($atomrec->{resSeq} == 16 && $atomrec->{resName} eq 'R') ||
			($atomrec->{resSeq} == 80 && $atomrec->{resName} eq 'G') ||
			($atomrec->{resSeq} == 81 && $atomrec->{resName} eq 'N')) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1eta/) {
		if ($atomrec->{chainID} eq '1' || 
		    $atomrec->{chainID} eq '2') {
		    if ($atomrec->{resSeq} == 30 && $atomrec->{iCode} eq 'B') {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1fcs/) {
		if ($atomrec->{chainID} eq '_') {
		    if ($atomrec->{resSeq} == 0) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1fm5/) {
		if ($atomrec->{chainID} eq 'A') {
		    if ($atomrec->{resSeq} == 78 ||
			$atomrec->{resSeq} == 79 ||
			$atomrec->{resSeq} == 80) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1fh2/) {
		if ($atomrec->{chainID} eq 'A' || $atomrec->{chainID} eq 'B') {
		    if (($atomrec->{resSeq} ==  30 && $atomrec->{resName} eq 'M') ||
			($atomrec->{resSeq} == 119 && $atomrec->{resName} eq 'T')) {
			next;
		    }
		}
	    }
#	    elsif ($pdbID =~ /^1ft8/) {                    
#		if ($atomrec->{chainID} eq 'E') {
###		    &abort ("cannot handle this pdb $pdbID");
##		    if ($atomrec->{resSeq} >= 160) {
##			next;
##		    }
#		}
#	    }
	    elsif ($pdbID =~ /^1gtl/) {                    
		if ($atomrec->{chainID} eq '1' || $atomrec->{chainID} eq '2') {
		    if ($atomrec->{resSeq} >= 353) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1gw8/) {                    
		if ($atomrec->{chainID} eq 'B' || 
		    $atomrec->{chainID} eq 'E' || 
		    $atomrec->{chainID} eq 'H' || 
		    $atomrec->{chainID} eq 'K') {
		    if ($atomrec->{resSeq} <= 2016) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1hb7/) {                    
		if ($atomrec->{chainID} eq 'A' ||

		    $atomrec->{chainID} eq 'D' || 
		    $atomrec->{chainID} eq 'G' || 
		    $atomrec->{chainID} eq 'J') {
		    if (($atomrec->{resSeq} <= 1018) ||
			($atomrec->{resSeq} >= 1353 && $atomrec->{resSeq} <= 1357) ||
			($atomrec->{resSeq} == 1384)) {
			next;
		    }
		}
		elsif ($atomrec->{chainID} eq 'B' ||
		       $atomrec->{chainID} eq 'E' || 
		       $atomrec->{chainID} eq 'H' || 
		       $atomrec->{chainID} eq 'K') {
		    if (($atomrec->{resSeq} <= 2016) ||
			($atomrec->{resSeq} >= 2353 && $atomrec->{resSeq} <= 2357) ||
			($atomrec->{resSeq} == 2384)) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1hhk/) {
		if ($atomrec->{chainID} eq 'A') {
		    if ($atomrec->{resSeq} >= 271) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1hj9/) {                    
		if ($atomrec->{chainID} eq 'A') {
		    if ($atomrec->{resSeq} >= 184 && $atomrec->{resSeq} <= 188) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1i7w/) {                    
		if ($atomrec->{chainID} eq 'D') {
		    if ($atomrec->{resSeq} >= 632 && $atomrec->{resSeq} <= 635) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1jdf/) {                    
		if ($atomrec->{chainID} eq 'D') {
		    if ($atomrec->{resSeq} >= 5 && $atomrec->{resSeq} <= 6) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1jhf/) {
		if ($atomrec->{chainID} eq 'B') {
		    if (($atomrec->{resSeq} == 75 && $atomrec->{resName} eq 'G') ||
			($atomrec->{resSeq} == 76 && $atomrec->{resName} eq 'L') ||
			($atomrec->{resSeq} == 77 && $atomrec->{resName} eq 'P') ||
			($atomrec->{resSeq} == 78 && $atomrec->{resName} eq 'L') ||
			($atomrec->{resSeq} == 79 && $atomrec->{resName} eq 'V')) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1jhr/) {
		if ($atomrec->{chainID} eq 'A') {
		    if ($atomrec->{resSeq} == 346) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1jj2/) {
		if ($atomrec->{chainID} eq 'G') {
		    if ($atomrec->{resSeq} >= 63) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1jni/) {
		if ($atomrec->{chainID} eq 'A') {
		    if ($atomrec->{resSeq} == 37 ||
			$atomrec->{resSeq} == 92) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1jxt/) {
		if ($atomrec->{chainID} eq 'A') {
		    if (($atomrec->{resSeq} == 22 && $atomrec->{resName} eq 'S') ||
			($atomrec->{resSeq} == 25 && $atomrec->{resName} eq 'I')) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1jxu/) {
		if ($atomrec->{chainID} eq 'A') {
		    if (($atomrec->{resSeq} == 22 && $atomrec->{resName} eq 'S') ||
			($atomrec->{resSeq} == 25 && $atomrec->{resName} eq 'I')) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1jxw/) {
		if ($atomrec->{chainID} eq 'A') {
		    if (($atomrec->{resSeq} == 22 && $atomrec->{resName} eq 'S') ||
			($atomrec->{resSeq} == 25 && $atomrec->{resName} eq 'I')) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1jxx/) {
		if ($atomrec->{chainID} eq 'A') {
		    if (($atomrec->{resSeq} == 22 && $atomrec->{resName} eq 'S') ||
			($atomrec->{resSeq} == 25 && $atomrec->{resName} eq 'I')) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1jxy/) {
		if ($atomrec->{chainID} eq 'A') {
		    if (($atomrec->{resSeq} == 22 && $atomrec->{resName} eq 'S') ||
			($atomrec->{resSeq} == 25 && $atomrec->{resName} eq 'I')) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1k8a/) {
		if ($atomrec->{chainID} eq 'I') {
		    if ($atomrec->{resSeq} >= 63) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1k9m/) {
		if ($atomrec->{chainID} eq 'I') {
		    if ($atomrec->{resSeq} >= 63) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1kd1/) {
		if ($atomrec->{chainID} eq 'I') {
		    if ($atomrec->{resSeq} >= 63) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1kqs/) {
		if ($atomrec->{chainID} eq 'G') {
		    if ($atomrec->{resSeq} >= 63) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1lta/) {
		if ($atomrec->{chainID} eq 'C') {
		    if ($atomrec->{resSeq} >= 237) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1ltb/) {
		if ($atomrec->{chainID} eq 'C') {
		    if ($atomrec->{resSeq} >= 233) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1m1k/) {
		if ($atomrec->{chainID} eq 'I') {
		    if ($atomrec->{resSeq} >= 63) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1m90/) {
		if ($atomrec->{chainID} eq 'I') {
		    if ($atomrec->{resSeq} >= 63) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1mdc/) {
		if ($atomrec->{chainID} eq '_') {
		    if ($atomrec->{resSeq} == 132) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1mnd/) {
		if ($atomrec->{chainID} eq '_') {
		    if ($atomrec->{resSeq} >= 685) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^1scf/) {
		if ($atomrec->{chainID} eq 'B') {
		    if ($atomrec->{resSeq} >= 137) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^2hmq/) {
		if ($atomrec->{chainID} eq 'A' || 
		    $atomrec->{chainID} eq 'B' || 
		    $atomrec->{chainID} eq 'C' || 
		    $atomrec->{chainID} eq 'D') {
		    if ($atomrec->{resSeq} == 21 && $atomrec->{iCode} eq 'A') {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^2hmz/) {
		if ($atomrec->{chainID} eq 'A' ||
		    $atomrec->{chainID} eq 'B' ||
		    $atomrec->{chainID} eq 'C' ||
		    $atomrec->{chainID} eq 'D') {
		    if ($atomrec->{resSeq} == 21 && $atomrec->{iCode} eq 'A') {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^2pf2/) {
		if ($atomrec->{chainID} eq '_') {
		    if ($atomrec->{resSeq} == 15 ||
			$atomrec->{resSeq} == 146) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^2plv/) {
		if ($atomrec->{chainID} eq '1') {
		    if (($atomrec->{resSeq} ==  6 && $atomrec->{resName} eq 'G') ||
			($atomrec->{resSeq} ==  7 && $atomrec->{resName} eq 'S') ||
			($atomrec->{resSeq} ==  8 && $atomrec->{resName} eq 'S') ||
			($atomrec->{resSeq} ==  9 && $atomrec->{resName} eq 'S') ||
			($atomrec->{resSeq} == 10 && $atomrec->{resName} eq 'T')) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^3hhr/) {
		if ($atomrec->{chainID} eq 'C') {
		    if ($atomrec->{resSeq} >= 235) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^8gpb/) {
		if ($atomrec->{chainID} eq '_') {
		    if ($atomrec->{resSeq} >= 840) {
			next;
		    }
		}
	    }
	    elsif ($pdbID =~ /^9pti/) {
		if ($atomrec->{chainID} eq '_') {
		    if ($atomrec->{resSeq} == 57 ||
			$atomrec->{resSeq} == 58) {
			next;
		    }
		}
	    }


	    # get coordseq
	    #
	    if (! defined $last_ResSeq_iCode->{"$chain"}) {
		if ($atomrec->{resName} ne 'X' && $atomrec->{resName} ne 'Z') {
		    $coordseq_i->{"$chain"} = 0;              # res_i starts at 0
		    $last_ResSeq_iCode->{"$chain"} = $atomrec->{resSeq}.$atomrec->{iCode};
		    $last_ResName->{"$chain"} = $atomrec->{resName};
		    push (@{$coordseq->{"$chain"}}, $atomrec->{resName});
		}
	    }
	    elsif ($atomrec->{resSeq}.$atomrec->{iCode} ne $last_ResSeq_iCode->{"$chain"} ||
		   $atomrec->{resName} ne $last_ResName->{"$chain"}) {
		if ($atomrec->{resName} ne 'X' && $atomrec->{resName} ne 'Z') {
		    ++$coordseq_i->{"$chain"};
		    $last_ResSeq_iCode->{"$chain"} = $atomrec->{resSeq}.$atomrec->{iCode};
		    $last_ResName->{"$chain"} = $atomrec->{resName};
		    push (@{$coordseq->{"$chain"}}, $atomrec->{resName});
		}
	    }
	    $atomrec->{coordseq_i} = $coordseq_i->{"$chain"};
	    
	    # pick first altLoc seen and stick with it for this residue
	    if ($atomrec->{altLoc} && ! defined $altLoc_used->{"$chain"}->[$coordseq_i->{"$chain"}]) {
		$altLoc_used->{"$chain"}->[$coordseq_i->{"$chain"}] = $atomrec->{altLoc};
	    }

	    # store CA coords for use in alignment
	    if ($atomrec->{name} eq 'CA' 
		&& (! $atomrec->{altLoc} || $atomrec->{altLoc} eq $altLoc_used->{"$chain"}->[$coordseq_i->{"$chain"}])) {
		$CA_coords->{"$chain"}->[$coordseq_i->{"$chain"}] = +[$atomrec->{x},
								      $atomrec->{y},
								      $atomrec->{z}];
	    }

	    # still not totally debugged?
	    if (! $atomrec->{altLoc} || 
		$atomrec->{altLoc} eq $altLoc_used->{"$chain"}->[$coordseq_i->{"$chain"}]) {
		push (@{$bodybufs->{"$chain"}}, $pdbbuf[$i]);
		push (@{$atoms->{"$chain"}}, $atomrec);
	    }
	}
	elsif ($recType =~ /^ANISOU/) {
	    next;
	}
	elsif ($recType =~ /^TER/) {
	    my $atom_n = -1;
	    if (length ($pdbbuf[$i]) >= 11) {
		$atom_n = substr ($pdbbuf[$i], 6, 5);
#		$atom_n =~ s/^\s+|\s+$//g;                       # not necessary
	    }
	    else {
		$atom_n = undef;
	    }

	    # special cases
	    if    ($pdbID =~ /^1acb/ && ($atom_n == 88 || $atom_n == 1069))   { next; }
	    elsif ($pdbID =~ /^1cho/ && ($atom_n == 1069))                    { next; }
	    elsif ($pdbID =~ /^1dv4/)                                         { next; }
	    elsif ($pdbID =~ /^1gha/ && ($atom_n == 70 || $atom_n == 1051))   { next; }
	    elsif ($pdbID =~ /^1ghb/ && ($atom_n == 70 || $atom_n == 1051))   { next; }
	    elsif ($pdbID =~ /^1gmc/ && ($atom_n == 70 || $atom_n == 1051))   { next; }
	    elsif ($pdbID =~ /^1gmd/ && ($atom_n == 70 || $atom_n == 1051))   { next; }
	    elsif ($pdbID =~ /^1gmh/ && ($atom_n == 69 || $atom_n == 1046))   { next; }
#	    elsif ($pdbID =~ /^1h68/ && ($atom_n == 1629))                   { next; }
	    elsif ($pdbID =~ /^1rnu/ && ($atom_n == 122))                     { next; }
	    elsif ($pdbID =~ /^1rnv/ && ($atom_n == 123))                     { next; }
	    elsif ($pdbID =~ /^1qd7/)                                         { next; }
	    elsif ($pdbID =~ /^1tbg/ && ($atom_n == 3125))                    { next; }
	    elsif ($pdbID =~ /^2cha/ && ($atom_n == 1069))                    { next; }
	    elsif ($pdbID =~ /^2csm/ && ($atom_n == 1744))                    { next; }
	    elsif ($pdbID =~ /^2gch/ && ($atom_n == 1069))                    { next; }
	    elsif ($pdbID =~ /^2gct/ && ($atom_n == 1069))                    { next; }
	    elsif ($pdbID =~ /^2rns/ && ($atom_n == 122))                     { next; }
	    elsif ($pdbID =~ /^3gch/ && ($atom_n == 1069))                    { next; }
	    elsif ($pdbID =~ /^3gct/ && ($atom_n == 1069))                    { next; }
	    elsif ($pdbID =~ /^4cha/ && ($atom_n == 1069 || $atom_n == 2841)) { next; }
	    elsif ($pdbID =~ /^4gch/ && ($atom_n == 1069))                    { next; }
	    elsif ($pdbID =~ /^4ovo/ && ($atom_n == 130))                     { next; }
	    elsif ($pdbID =~ /^5cha/ && ($atom_n == 1069 || $atom_n == 2841)) { next; }
	    elsif ($pdbID =~ /^5gch/ && ($atom_n == 1069))                    { next; }
	    elsif ($pdbID =~ /^6cha/ && ($atom_n == 1069 || $atom_n == 2852)) { next; }
	    elsif ($pdbID =~ /^6gch/ && ($atom_n == 1069))                    { next; }
	    elsif ($pdbID =~ /^7gch/ && ($atom_n == 1069))                    { next; }
	    elsif ($pdbID =~ /^8gch/ && ($atom_n == 72 || $atom_n == 1051))   { next; }
	    else {
		$ter_found{"$chain"} = 'true';
	    }
	    next;
	}
	elsif ($recType =~ /^ENDMDL/) {
	    $endmdl_found{"$chain"} = 'true';
	    next;
	}
	else {
	    #last;
	    next;
	}
    }


    # read HETATM records
    #
    $modres_hetatmbuf = +{};
    for (my $i=0; $i <= $#pdbbuf; ++$i) {
	$recType = substr ($pdbbuf[$i], 0, 6);
	if ($recType =~ /^HETATM/) {
	    my $deriv = substr ($pdbbuf[$i], 17, 3);
	    my $chain = substr ($pdbbuf[$i], 21, 1);
	    if ($modres->{$deriv} || &mapResCode($deriv, 'SILENT') ne 'X') {
		push (@{$modres_hetatmbuf->{$chain}}, $pdbbuf[$i]);
	    } else {
		push (@$hetatmbuf, $pdbbuf[$i]);
	    }
	}
    }


    # get register for each chain
    #
#    foreach $chain (keys %$seqres) {      # this is what we ultimately want to use
    foreach $chain ($chainOI) {

#debug
#print "seqres:   ";
#for (my $i=0; $i <= $#{$seqres->{"$chain"}}; ++$i) {
#    print $seqres->{"$chain"}->[$i];
#}
#print "\n";
#print "coordseq: ";
#for (my $i=0; $i <= $#{$coordseq->{"$chain"}}; ++$i) {
#    print $coordseq->{"$chain"}->[$i];
#}
#print "\n";
# end debug

	if (! $seqres->{"$chain"}) {
	    &abort ("no SEQRES found in header for pdb $pdbID, chain $chain");
	}
	if (! $coordseq->{"$chain"}) {
	    &abort ("no COORDSEQ found in body for pdb $pdbID, chain $chain");
	}
	if ($#{$coordseq->{"$chain"}} > $#{$seqres->{"$chain"}}) {
	    &abort ("coordseq longer than SEQRES for id $pdbID, chain $chain (". $#{$coordseq->{"$chain"}} ." > ". $#{$seqres->{"$chain"}}. ")\n"."seqres:   '".join ('', @{$seqres->{"$chain"}})."'\n"."coordseq: '".join('',@{$coordseq->{"$chain"}})."'");
	}

	# the alignment itself
	#
	my $alignment = &alignCoordSeqToSEQRES ('L', $seqres->{"$chain"}, $coordseq->{"$chain"}, $CA_coords->{"$chain"});


	# check the alignment visually
	#
	if ($debug) {
	    $lastAlignment = -1;
	    print STDERR "seqres  : ";
	    for ($i=0; $i <= $#{$seqres->{"$chain"}}; ++$i) {
		if (! defined $alignment->{'1to2'}->[$i]
		    || $alignment->{'1to2'}->[$i] == -1) {
		    print STDERR $seqres->{"$chain"}->[$i];
		} else {
		    print STDERR '-' while ($alignment->{'1to2'}->[$i] > ++$lastAlignment);
		    print STDERR $seqres->{"$chain"}->[$i];
		}
	    }
	    print STDERR "\n";
	    $lastAlignment = -1;
	    print STDERR "coordseq: ";
	    for ($i=0; $i <= $#{$coordseq->{"$chain"}}; ++$i) {
		if (! defined $alignment->{'2to1'}->[$i]
		    || $alignment->{'2to1'}->[$i] == -1) {
		    print STDERR $coordseq->{"$chain"}->[$i];
		} else {
		    print STDERR '-' while ($alignment->{'2to1'}->[$i] > ++$lastAlignment);
		    print STDERR $coordseq->{"$chain"}->[$i];
		}
	    }
	    print STDERR "\n";

#debug
#for ($i=0; $i <= $#{$coordseq->{"$chain"}}; ++$i) {
#    print $coordseq->{"$chain"}->[$i]." $i -> ". $alignment->{'2to1'}->[$i] ."\n";
#}
#print "\n";
	}


	for (my $j=0; $j <= $#{$atoms->{"$chain"}}; ++$j) {
# debug
#print "res: '". $atoms->{"$chain"}->[$j]->{resName}."' coordseq_i: '".$atoms->{"$chain"}->[$j]->{coordseq_i}."' -> map: '". $alignment->{'2to1'}->[$atoms->{"$chain"}->[$j]->{coordseq_i}]."'\n";

	    # fix ATOM identifier (from 'HETATM')
            substr ($bodybufs->{"$chain"}->[$j], 0, 6) = sprintf ("%6s", 'ATOM  ');

	    # fix resName field
	    my $code3 = &mapResCode ($atoms->{"$chain"}->[$j]->{resName});
	    substr ($bodybufs->{"$chain"}->[$j], 17, 3) = sprintf ("%3s", $code3);

	    # get resmap, fix resSeq field, and store resSeqReal
	    $coordseq_n = substr ($bodybufs->{"$chain"}->[$j], 22, 5);  # + insert code
	    $coordseq_n =~ s/\s+//g;
	    $coordseq_i = $atoms->{"$chain"}->[$j]->{coordseq_i};
	    
	    if (! defined $alignment->{'2to1'}->[$coordseq_i]) {
		&abort ("some coords did not align to SEQRES at res $coordseq_n");
	    }
	    $real_res_i = $alignment->{'2to1'}->[$coordseq_i];
	    $resmap->{"$chain"}->[$real_res_i] = $coordseq_n;
	    $atoms->{"$chain"}->[$j]->{resSeqReal} = $real_res_i + 1;
	    substr ($bodybufs->{"$chain"}->[$j], 22, 4) = sprintf ("%4d", $real_res_i + 1);

	    # remove any insertion code
	    substr ($bodybufs->{"$chain"}->[$j], 26, 1) = ' ';

	    # fix MSE/MET 'SE  ' atom type (selenium is really sulfur-delta)
	    if ($atoms->{"$chain"}->[$j]->{resName} eq 'M' && $atoms->{"$chain"}->[$j]->{name} eq 'SE') {
		substr ($bodybufs->{"$chain"}->[$j], 12, 4) = sprintf ("%4s", ' SD ');
		$atoms->{"$chain"}->[$j]->{name} = 'SD';
	    }
	    # fix CSE/CYS 'SE  ' atom type (selenium is really sulfur-gamma)
	    if ($atoms->{"$chain"}->[$j]->{resName} eq 'C' && $atoms->{"$chain"}->[$j]->{name} eq 'SE') {
		substr ($bodybufs->{"$chain"}->[$j], 12, 4) = sprintf ("%4s", ' SG ');
		$atoms->{"$chain"}->[$j]->{name} = 'SG';
	    }
	}
    }


    # create atom records organized by residue and atom type
    #
    #foreach $chain (sort keys %$seqres) {
    foreach $chain ($chainOI) {
	$chain_display = ($chain eq '_') ? ' ' : $chain;

	# zero records
	#
	for ($res_i=0; $res_i <= $#{$seqres->{"$chain"}}; ++$res_i) {
	    $res_type = $seqres->{"$chain"}->[$res_i];
	    @atom_recs = &resAtomRecs ($res_type);
	    if (@atom_recs) {
		foreach $rec (@atom_recs) {
		    substr ($rec, 22, 4) = sprintf ("%4d", $res_i+1);
		    $atom_type = substr ($rec, 12, 4);
		    $atom_type =~ s/^\s+|\s+$//g;
		    $res->{"$chain"}->[$res_i]->{atoms}->{$atom_type}->{rec} = $rec;
		}
	    }
	}

	# get coord data
	#
	for ($i=0; $i <= $#{$atoms->{"$chain"}}; ++$i) {
	    $res_i     = $atoms->{"$chain"}->[$i]->{resSeqReal} - 1;
	    $atom_type = $atoms->{"$chain"}->[$i]->{name};
	    $rec = $bodybufs->{"$chain"}->[$i];
	    $rec_len = length $rec;
	    $rec .= ' 'x(80-$rec_len)  if ($rec_len < 80);
	    $occ = substr ($rec, 54, 6);
	    if ($occ =~ /^\s+$/) {
		substr ($rec, 54, 6) = sprintf ("%6.2f", 1.00);
	    } elsif ($occ == 0) {
		substr ($rec, 54, 6) = sprintf ("%6.2f", 0.10);
	    }
	    $res->{"$chain"}->[$res_i]->{atoms}->{$atom_type}->{rec} = $rec;
	}

	# restore missing modified residues that are available in @$modres_hetatmbuf
	#
	for ($res_i=0; $res_i <= $#{$seqres->{"$chain"}}; ++$res_i) {

	    # skip residues that are already defined
	    $defined_ca = undef;
	    if (defined $res->{"$chain"}->[$res_i]->{atoms}->{'CA'}->{rec}) {
		$res_rec = $res->{"$chain"}->[$res_i]->{atoms}->{'CA'}->{rec};
		$defined_ca = (substr ($res_rec, 54, 6) >= 0);
	    }
	    next if ($defined_ca);

	    # only fill in missing MODRES residues
	    next  if (! defined $res_rec);
	    $res_type = substr ($res_rec, 17, 3);
	    $likely_modres = undef;
	    if (defined $modres_reverse->{$res_type}) {
		$likely_modres = 'true';
	    }
	    else {
		for ($modres_atom_i=0; $modres_atom_i <= $#{$modres_hetatmbuf->{$chain}}; ++$modres_atom_i) {
		    $hetatm_res_type = substr ($modres_hetatmbuf->{$chain}->[$modres_atom_i], 17, 3);
		    if (&mapResCode ($hetatm_res_type) eq &mapResCode ($res_type)) {
			$likely_modres = 'true';
			last;
		    }
		}
	    }
	    if ($likely_modres) {
		for ($modres_atom_i=0; $modres_atom_i <= $#{$modres_hetatmbuf->{$chain}}; ++$modres_atom_i) {
		    $hetatm_atom_type = substr ($modres_hetatmbuf->{$chain}->[$modres_atom_i], 12, 4);
		    $hetatm_res_type  = substr ($modres_hetatmbuf->{$chain}->[$modres_atom_i], 17, 3);
		    $hetatm_res_n     = substr ($modres_hetatmbuf->{$chain}->[$modres_atom_i], 22, 5); # + insertion code

		    # check if right hetatm type
		    if ($hetatm_atom_type eq ' CA ' &&
			((defined $modres_reverse->{$res_type} && $modres_reverse->{$res_type} eq $hetatm_res_type) ||
			(&mapResCode ($hetatm_res_type) eq &mapResCode ($res_type)))) {

			# check to see if the CA is in the right proximity to neighbor residues
			@modres_CA_pos = (substr ($modres_hetatmbuf->{$chain}->[$modres_atom_i], 30, 8),
					  substr ($modres_hetatmbuf->{$chain}->[$modres_atom_i], 38, 8),
					  substr ($modres_hetatmbuf->{$chain}->[$modres_atom_i], 46, 8));
			$defined_Nterm_ca = undef;
			if (defined $res->{"$chain"}->[$res_i-1]->{atoms}->{'CA'}->{rec}) {
			    $Nterm_res_rec = $res->{"$chain"}->[$res_i-1]->{atoms}->{'CA'}->{rec};
			    $defined_Nterm_ca = (substr ($Nterm_res_rec, 54, 6) >= 0);
			    @Nterm_CA_pos = (substr ($Nterm_res_rec, 30, 8),
					     substr ($Nterm_res_rec, 38, 8),
					     substr ($Nterm_res_rec, 46, 8));
			}

			$defined_Cterm_ca = undef;
			if (defined $res->{"$chain"}->[$res_i+1]->{atoms}->{'CA'}->{rec}) {
			    $Cterm_res_rec = $res->{"$chain"}->[$res_i+1]->{atoms}->{'CA'}->{rec};
			    $defined_Cterm_ca = (substr ($Cterm_res_rec, 54, 6) >= 0);
			    @Cterm_CA_pos = (substr ($Cterm_res_rec, 30, 8),
					     substr ($Cterm_res_rec, 38, 8),
					     substr ($Cterm_res_rec, 46, 8));
			}

			$fits_here = 'true';
			if (! $defined_Nterm_ca && ! $defined_Cterm_ca) {
			    $fits_here = undef;
			}
			if ($defined_Nterm_ca && ! &CA_neighbor (\@modres_CA_pos, \@Nterm_CA_pos)) {
			    $fits_here = undef;
			}
			if ($defined_Cterm_ca && ! &CA_neighbor (\@modres_CA_pos, \@Cterm_CA_pos)) {
			    $fits_here = undef;
			}

			# if it fits spatially, then insert into $res structure
			if (! $fits_here) {
			    next;
			}
			else {
                            # DEBUG
			    #print "$res_type $res_i has a HETATM fit\n";			    
			    $found_hetatms = undef;
			    for ($modres_atom_j=0; $modres_atom_j <= $#{$modres_hetatmbuf->{$chain}}; ++$modres_atom_j) {
				$hetatm_res_m = substr ($modres_hetatmbuf->{$chain}->[$modres_atom_j], 22, 5); # + insertion code

				# found our atom records
				if ($hetatm_res_m eq $hetatm_res_n) {

				    $found_hetatms = 'true';

				    # work on $rec instead of $modres_hetatmbuf
				    my $rec = $modres_hetatmbuf->{$chain}->[$modres_atom_j];
				    $rec_len = length $rec;
				    $rec .= ' 'x(80-$rec_len)  if ($rec_len < 80);

				    # we need atom_type now
				    $atom_type = substr ($rec, 12, 4);
				    $atom_type =~ s/^\s+|\s+$//g;

				    # fix MSE/MET 'SE  ' atom type (selenium is really sulfur-delta)
				    if ($res_type eq 'MET' && $atom_type eq 'SE') {
					substr ($rec, 12, 4) = sprintf ("%4s", ' SD ');
					$atom_type = 'SD';
				    }
				    # fix CSE/CYS 'SE  ' atom type (selenium is really sulfur-gamma)
				    if ($res_type eq 'CYS' && $atom_type eq 'SE') {
					substr ($rec, 12, 4) = sprintf ("%4s", ' SG ');
					$atom_type = 'SG';
				    }

				    # only replace regular atoms in protein
				    if (defined $res->{"$chain"}->[$res_i]->{atoms}->{$atom_type}->{rec}) {

					# fix ATOM identifier (from 'HETATM')
					substr ($rec, 0, 6) = sprintf ("%6s", 'ATOM  ');
					
					# fix resName field
					my $code3 = $res_type;
					substr ($rec, 17, 3) = sprintf ("%3s", $code3);
					
					# get resmap, fix resSeq field, and store resSeqReal
					$coordseq_n = $hetatm_res_n;
					$coordseq_n =~ s/\s+//g;
					$real_res_i = $res_i;
					$resmap->{"$chain"}->[$real_res_i] = $coordseq_n;
					substr ($rec, 22, 4) = sprintf ("%4d", $real_res_i + 1);
					
					# remove any insertion code
					substr ($rec, 26, 1) = ' ';
					
					# and finally assign the record to $res structure
					$res->{"$chain"}->[$res_i]->{atoms}->{$atom_type}->{rec} = $rec;
				    }
				}
				elsif ($found_hetatms) {
				    last;
				}
			    }
			    last;
			}
		    }
		}
	    }
	}


	# remove stray solitary residues (must be 3 res long to keep)
	#
	for ($res_i=0; $res_i <= $#{$seqres->{"$chain"}}; ++$res_i) {
	    $defined_ca = undef;
	    if (defined $res->{"$chain"}->[$res_i]->{atoms}->{'CA'}->{rec}) {
		$res_rec = $res->{"$chain"}->[$res_i]->{atoms}->{'CA'}->{rec};
		$defined_ca = (substr ($res_rec, 54, 6) >= 0);
	    }
	    next if (! $defined_ca);
	    
	    $defined_ca_pos1 = undef;
	    $defined_ca_pos2 = undef;
	    $defined_ca_pos3 = undef;
	    
	    # last of three
	    if ($res_i > 1) {
		if (defined $res->{"$chain"}->[$res_i-2]->{atoms}->{'CA'}->{rec}) {
		    $res_rec = $res->{"$chain"}->[$res_i-2]->{atoms}->{'CA'}->{rec};
		    $defined_ca_pos1 = (substr ($res_rec, 54, 6) >= 0);
		}
		if (defined $res->{"$chain"}->[$res_i-1]->{atoms}->{'CA'}->{rec}) {
		    $res_rec = $res->{"$chain"}->[$res_i-1]->{atoms}->{'CA'}->{rec};
		    $defined_ca_pos2 = (substr ($res_rec, 54, 6) >= 0);
		}
		next if ($defined_ca_pos1 && $defined_ca_pos2);
	    }
	    # first of three
	    if ($res_i < $#{$seqres->{"$chain"}} - 1) {
		if (defined $res->{"$chain"}->[$res_i+1]->{atoms}->{'CA'}->{rec}) {
		    $res_rec = $res->{"$chain"}->[$res_i+1]->{atoms}->{'CA'}->{rec};
		    $defined_ca_pos2 = (substr ($res_rec, 54, 6) >= 0);
		}
		if (defined $res->{"$chain"}->[$res_i+2]->{atoms}->{'CA'}->{rec}) {
		    $res_rec = $res->{"$chain"}->[$res_i+2]->{atoms}->{'CA'}->{rec};
		    $defined_ca_pos3 = (substr ($res_rec, 54, 6) >= 0);
		}
		next if ($defined_ca_pos2 && $defined_ca_pos3);
	    } 
	    # middle of three
	    if (defined $res->{"$chain"}->[$res_i-1]->{atoms}->{'CA'}->{rec}) {
		$res_rec = $res->{"$chain"}->[$res_i-1]->{atoms}->{'CA'}->{rec};
		$defined_ca_pos1 = (substr ($res_rec, 54, 6) >= 0);
	    }
	    if (defined $res->{"$chain"}->[$res_i+1]->{atoms}->{'CA'}->{rec}) {
		$res_rec = $res->{"$chain"}->[$res_i+1]->{atoms}->{'CA'}->{rec};
		$defined_ca_pos3 = (substr ($res_rec, 54, 6) >= 0);
	    }
	    next if ($defined_ca_pos1 && $defined_ca_pos3);

	    # not member of three contiguous, zero that bad boy
	    $res_type = $seqres->{"$chain"}->[$res_i];
	    @atom_recs = &resAtomRecs ($res_type);
	    if (@atom_recs) {
		foreach $rec (@atom_recs) {
		    substr ($rec, 22, 4) = sprintf ("%4d", $res_i+1);
		    $atom_type = substr ($rec, 12, 4);
		    $atom_type =~ s/^\s+|\s+$//g;
		    $res->{"$chain"}->[$res_i]->{atoms}->{$atom_type}->{rec} = $rec;
		}
	    }
	}

	# parse important fields into data structure
	#
	for ($res_i=0, $line_i=0; $res_i <= $#{$seqres->{"$chain"}}; ++$res_i) {
	    next if (! defined $res->{"$chain"}->[$res_i]->{atoms}->{'CA'}->{rec});
	    $rec = $res->{"$chain"}->[$res_i]->{atoms}->{'CA'}->{rec};
	    $res->{"$chain"}->[$res_i]->{resName} = &mapResCode (substr ($rec, 17, 3));

	    foreach $atom_type (keys %{$res->{"$chain"}->[$res_i]->{atoms}}) {
		# tidy up $rec
		$rec = $res->{"$chain"}->[$res_i]->{atoms}->{$atom_type}->{rec};
		$rec_len = length $rec;
		$rec .= ' 'x(80-$rec_len)  if ($rec_len < 80);
		substr ($rec, 6, 5)  = sprintf ("%5d", ++$line_i);
		substr ($rec, 21, 1) = sprintf ("%1s", $chain_display);
		$res->{"$chain"}->[$res_i]->{atoms}->{$atom_type}->{rec} = $rec;

		# get field vals
		$x   = substr ($rec, 30, 8);
		$y   = substr ($rec, 38, 8);
		$z   = substr ($rec, 46, 8);
		$occ = substr ($rec, 54, 6);
		$b   = substr ($rec, 60, 6);
		$chg = substr ($rec, 78, 2);
		$x   =~ s/^\s+|\s+$//g;
		$y   =~ s/^\s+|\s+$//g;
		$z   =~ s/^\s+|\s+$//g;
		$occ =~ s/^\s+|\s+$//g;
		$b   =~ s/^\s+|\s+$//g;
		$chg =~ s/^\s+|\s+$//g;
		$occ = "1.00"  if ($occ eq '');
		$b   = "0.00"  if ($occ eq '');
		$chg = "x"     if ($chg eq '');

		$res->{"$chain"}->[$res_i]->{atoms}->{$atom_type}->{x}   = $x;
		$res->{"$chain"}->[$res_i]->{atoms}->{$atom_type}->{y}   = $y;
		$res->{"$chain"}->[$res_i]->{atoms}->{$atom_type}->{z}   = $z;
		$res->{"$chain"}->[$res_i]->{atoms}->{$atom_type}->{occ} = $occ;
		$res->{"$chain"}->[$res_i]->{atoms}->{$atom_type}->{b}   = $b;
		$res->{"$chain"}->[$res_i]->{atoms}->{$atom_type}->{chg} = $chg;
	    }
	}
    }


    # attach data to pdb struct
    #
    $pdb->{headerbuf} = $headerbuf;
    $pdb->{bodybufs}  = $bodybufs;
    $pdb->{hetatmbuf} = $hetatmbuf;
    $pdb->{header}    = $header;
    $pdb->{seqres}    = $seqres;
    $pdb->{resmap}    = $resmap;
    $pdb->{modres}    = $modres;
    $pdb->{atoms}     = $atoms;
    $pdb->{res}       = $res;
    
    return ($pdb, $chainOI);
}


# parseAtomRec()
#
#
# PDB ATOM Record Format 
#
# COLUMNS        DATA TYPE       FIELD         DEFINITION
# -----------------------------------------------------------------------------
#  1 -  6        Record name     "ATOM  "
#  7 - 11        Integer         serial        Atom serial number.
# 13 - 16        Atom            name          Atom name.
# 17             Character       altLoc        Alternate location indicator.
# 18 - 20        Residue name    resName       Residue name.
# 22             Character       chainID       Chain identifier.
# 23 - 26        Integer         resSeq        Residue sequence number.
# 27             AChar           iCode         Code for insertion of residues.
# 31 - 38        Real(8.3)       x             Orthogonal coords for X in Ang
# 39 - 46        Real(8.3)       y             Orthogonal coords for Y in Ang
# 47 - 54        Real(8.3)       z             Orthogonal coords for Z in Ang
# 55 - 60        Real(6.2)       occupancy     Occupancy.
# 61 - 66        Real(6.2)       tempFactor    Temperature factor.
# 73 - 76        LString(4)      segID         Segment id, left-justified.
# 77 - 78        LString(2)      element       Element symbol, right-justified.
# 79 - 80        LString(2)      charge        Charge on the atom.
#
#
sub parseAtomRec {
    local ($line, $modres) = @_;
    my $atomrec = +{};

    my $line_len = length ($line);
    &abort ("short atom line '$line'")  if ($line_len < 54);

    $atomrec->{rec}        = substr ($line, 0, 6);
    $atomrec->{serial}     = substr ($line, 6, 5);
    $atomrec->{name}       = substr ($line, 12, 4);
    $atomrec->{altLoc}     = substr ($line, 16, 1);
    $atomrec->{resName}    = substr ($line, 17, 3);
    $atomrec->{chainID}    = substr ($line, 21, 1);
    $atomrec->{resSeq}     = substr ($line, 22, 4);
    $atomrec->{iCode}      = substr ($line, 26, 1);
    $atomrec->{x}          = substr ($line, 30, 8);
    $atomrec->{y}          = substr ($line, 38, 8);
    $atomrec->{z}          = substr ($line, 46, 8);
    $atomrec->{occupancy}  = ($line_len < 60) ? '1.00' : substr ($line, 54, 6);
    $atomrec->{tempFactor} = ($line_len < 66) ? '0.00' : substr ($line, 60, 6);
    $atomrec->{segID}      = ($line_len < 76) ? ''     : substr ($line, 72, 4);
    $atomrec->{element}    = ($line_len < 78) ? ''     : substr ($line, 76, 2);
    $atomrec->{charge}     = ($line_len < 80) ? ''     : substr ($line, 78, 2);
    
    if (length $atomrec->{resName} == 3) {
	if ($modres->{$atomrec->{resName}}) {
	    $atomrec->{resName} = $modres->{$atomrec->{resName}};
	}
	$atomrec->{resName} = &mapResCode ($atomrec->{resName});
    }
    
    foreach $k (keys %$atomrec) {
	# note: sometimes fields are missing!
	$atomrec->{$k} =~ s/^\s+|\s+$//g  if (defined $atomrec->{$k});
    }

    $atomrec->{chainID} = '_' if ($atomrec->{chainID} =~ /^\s*$/);

    return $atomrec;
}


# alignCoordSeqToSEQRES()
#
sub alignCoordSeqToSEQRES {
    my ($scope, $seq1, $seq2, $CA) = @_;
    my $alignment = +{};

    # straightforward linear ungapped (less expensive)
    #
    $alignment = &noGapLinearAlignment ($seq1, $seq2);

    if (! $alignment) {
	# gapping in seq2 alignment 
	#
	$alignment = &alignSeqsNoGaps1 ($scope, $seq1, $seq2, $CA);

	# alignment must be perfect (that means > 1000*(n-1) residues) 
	#
	if ($alignment->{aligned_residues_cnt} != $#{$seq2}+1) {
#debug
print "lenseq: ".($#{$seq2}+1)."\n";
print "score : ".$alignment->{aligned_residues_cnt}."\n";

#	    &abort ("imperfect alignment between coordseq and SEQRES for pdb $pdbID, chain $chain");
	}
    }

    return $alignment;
}


# noGapLinearAlignment ()
#
sub noGapLinearAlignment {
    my ($seq1, $seq2) = @_;
    my $alignment = undef;

    my $register_found = undef;
    my $register_shift = -1;
    my $seq1_copy = +[];
    for (my $i=0; $i <= $#{$seq1}; ++$i) {
	$seq1_copy->[$i] = $seq1->[$i];
    }

    for (my $base_i=0; $base_i <= $#{$seq1} - $#{$seq2} && ! $register_found; ++$base_i) {

	if (&getIdentity ($seq1_copy, $seq2) == 1) {
	    $register_found = 'true';
	    for (my $i=0; $i <= $#{$seq2}; ++$i) {
		$alignment->{'2to1'}->[$i] = $base_i + $i;
		$alignment->{'1to2'}->[$base_i+$i] = $i;
	    }
	    for (my $i=0; $i < $base_i; ++$i) {
		$alignment->{'1to2'}->[$i] = -1;
	    }
	    for (my $i=$base_i+$#{$seq2}+1; $i <= $#{$seq1}; ++$i) {
		$alignment->{'1to2'}->[$i] = -1;
	    }
	}
	else {
	    shift @{$seq1_copy};
	}
    }
    
    return $alignment;
}


# alignSeqsNoGaps1()
#
sub alignSeqsNoGaps1 {
    my ($scope, $seq1, $seq2, $CA) = @_;
    my $alignment = +{};
    my $V = +[];                       # score for opt Q1...Qi<>T1...Tj
    #NOGAP my $E = +[];                # score for opt Q1...Qi<>T1...Tj, - <>Tj
    my $F = +[];                       # score for opt Q1...Qi<>T1...Tj, Qi<> -
    my $G = +[];                       # score for opt Q1...Qi<>T1...Tj, Qi<>Tj
    my $Vsource = +[];
    my $Fsource = +[];
    #my $Esource = +[];
    my $walkback_state;
# debug
#    my $pair = 1000;
#    my $mispair = -1000;
#    my $gap_init = 11;
#    my $gap_ext = 1;
    my $init_score;
    my $ext_score;
    my $pair = 1000;
    my $mispair = -10000;
    my $base_gap_init = 500;  
#    my $gap_ext = 0;  
    my $gap_ext = 0.1;  # changed from 0 to 0.1 on 2003-03-21
    my $gap_init = $base_gap_init;
    my $INT_MIN = -1000000;
    my ($i, $j);

# debug
#print "seq1: '";
#for ($i=0; $i<=$#{$seq1}; ++$i) {
#    print $seq1->[$i];
#}
#print "'\n";
#print "seq2: '";
#for ($i=0; $i<=$#{$seq2}; ++$i) {
#    print $seq2->[$i];
#}
#print "'\n";
#exit 0;
# end debug

    # basis
    #
    $V->[0]->[0] = 0;
    #$E->[0]->[0] = $INT_MIN;   # never actually accessed
    #$F->[0]->[0] = $INT_MIN;   # never actually accessed
    for ($i=1; $i <= $#{$seq1}+1; ++$i) {
	$V->[$i]->[0] = ($scope eq 'G')  ? -$gap_init - $i*$gap_ext  : 0;
	#NOGAP $E->[$i]->[0] = $INT_MIN;
    }
    for ($j=1; $j <= $#{$seq2}+1; ++$j) {
	$V->[0]->[$j] = ($scope eq 'G')  ? -$gap_init - $j*$gap_ext  : 0;
	$F->[0]->[$j] = $INT_MIN;
    }


    # recurrence
    #
    $maxScore = 0;
    for ($i=1; $i<= $#{$seq1}+1; ++$i) {
	for ($j=1; $j<= $#{$seq2}+1; ++$j) {

	    # note: seq1[i-1]==Qi and seq2[j-1]==Tj (i.e. res 1 stored as 0)

	    # check for closeness of residues to adjust gap_init so that ambiguous
	    #   residue placement on either side of missing density is resolved
	    #
	    $gap_init = $base_gap_init;
	    if ($j >= 1) {
		#if (! &CA_neighbor ($CA->[$j-1], $CA->[$j-2])) {
		if ($j <= $#{$seq2} && 
		    ! &CA_neighbor ($CA->[$j], $CA->[$j-1])) {
#		    $gap_init = -100;
		    $gap_init = 0;     # changed from -100 to 0 on 2003-10-14
		}
	    }


	    # G, F, and E
	    #
	    $G->[$i]->[$j] = $V->[$i-1][$j-1]
		+ &scorePair ($seq1->[$i-1], $seq2->[$j-1], $pair, $mispair);
	    $init_score = $V->[$i-1]->[$j] -$gap_init -$gap_ext;
	    $ext_score =  $F->[$i-1]->[$j] -$gap_ext;
	    if ($init_score > $ext_score) {
		$F->[$i]->[$j] = $init_score;
		$Fsource->[$i]->[$j] = 'V';
	    } else {
		$F->[$i]->[$j] = $ext_score;
		$Fsource->[$i]->[$j] = 'F';
	    }		
	    #NOGAP $E->[$i]->[$j] = &maxInt ($V->[$i][$j-1] -$gap_init -$gap_ext, $E->[$i]->[$j-1] -$gap_ext);


	    # V
	    #
	    # Local scope and null string superior
	    #NOGAP if ($scope eq 'L' && $F->[$i]->[$j] < 0 && $E->[$i]->[$j] < 0 && $G->[$i]->[$j] < 0) {
	    if ($scope eq 'L' && $F->[$i]->[$j] < 0 && $G->[$i]->[$j] < 0) {
		$V->[$i]->[$j] = 0;
		$Vsource->[$i]->[$j] = 'N';
	    }
	    # Global scope or null string inferior
	    else {
		#NOGAP if ($F->[$i]->[$j] >= $G->[$i]->[$j] || $E->[$i]->[$j] >= $G->[$i]->[$j]) {
		if ($F->[$i]->[$j] >= $G->[$i]->[$j]) {
		    #NOGAP if ($F->[$i]->[$j] > $E->[$i]->[$j]) {
			$V->[$i]->[$j] = $F->[$i]->[$j];
			$Vsource->[$i]->[$j] = 'F';
		    #NOGAP } else {
			#NOGAP $V->[$i]->[$j] = $E->[$i]->[$j];
			#NOGAP $Vsource->[$i]->[$j] = 'E';
		    #NOGAP }	
		} else {
		    $V->[$i]->[$j] = $G->[$i]->[$j];
		    $Vsource->[$i]->[$j] = 'G';
		}
	    }
	    
	    # maxScore
	    #
	    if ($V->[$i]->[$j] > $maxScore) {
		$maxScore = $V->[$i]->[$j];
		$maxScore_i = $i;
		$maxScore_j = $j;
	    }
	}
    }
    $alignment->{score} = ($scope eq 'G')
	                      ? $V->[$#{$seq1}+1]->[$#{$seq2}+1]
			      : $maxScore;
    
    # walk back
    #
    $walkback_state = 'V';
    if ($scope eq 'G') {
	$i = $#{$seq1}+1;
	$j = $#{$seq2}+1;
    } else {
	$i = $maxScore_i;
	$j = $maxScore_j;
    }
    $walkback_done = undef;
    while (! $walkback_done) {
	# Global stop condition
	if ($scope eq 'G' && ($i == 0 || $j == 0)) {
	    $walkback_done = 'TRUE';
	    last;
	}
	# Local stop condition
	elsif ($scope eq 'L' &&           # if(Vsource[i][j]=='N'): we hit null
	       ($V->[$i]->[$j] == 0 || $i == 0 || $j == 0)) {
	    $walkback_done = 'TRUE';
	    last;
	}
	
	# walking back in V
	if ($walkback_state eq 'V') {
	    if ($Vsource->[$i]->[$j] eq 'G') {
		#$walkback_state = 'V';                           # unnecessary
		$alignment->{'1to2'}->[$i-1] = $j-1;            # seq1[i-1]==Qi
		$alignment->{'2to1'}->[$j-1] = $i-1;            # seq2[j-1]==Tj
		++$alignment->{aligned_residues_cnt};
		--$i; --$j;
	    }
	    elsif ($Vsource->[$i]->[$j] eq 'F') {
		if ($Fsource->[$i]->[$j] eq 'V') { 
		    #$walkback_state = 'V';                       # unnecessary
		} else {
		    $walkback_state = 'F';
		}
		$alignment->{'1to2'}->[$i-1] = -1;
		--$i;
	    }
	#  NOGAP
	#  NOGAP   else {                          # Vsource->[$i]->[$j] eq 'E'
	#  NOGAP	if ($Esource->[$i]->[$j] eq 'V') {
	#  NOGAP	    $walkback_state = 'V';                # unnecessary
	#  NOGAP	} else {
	#  NOGAP	    $walkback_state = 'E';
	#  NOGAP	}	  
	#  NOGAP	$alignment->{'2to1'}->[$j-1] = -1;
	#  NOGAP	--$j;
	#  NOGAP    }
	}
	
	# walking back in F
	elsif ($walkback_state eq 'F') {
	    if ($Fsource->[$i]->[$j] eq 'V') { 
		$walkback_state = 'V';
	    }
	    #else {
	    #$walkback_state = 'F';                               # unnecessary
	    #}
	    $alignment->{'1to2'}->[$i-1] = -1;
	    --$i;
	}
	# NOGAP
	# NOGAP	# walking back in E
	# NOGAP	else { #if ($walkback_state eq 'E') {
	# NOGAP	    if ($Esource->[$i]->[$j] eq 'V') {
	# NOGAP		$walkback_state = 'V';
	# NOGAP	    }
	# NOGAP	    #else {
	# NOGAP	    #  $walkback_state = 'E';                     # unnecessary
	# NOGAP	    #}	  
	# NOGAP	    $alignment->{'2to1'}->[$j-1] = -1;
	# NOGAP	    --$j;
	# NOGAP	}
    }


    return ($alignment);
}                     


# scorePair()
#
sub scorePair {
    local ($r1, $r2, $pair, $mispair) = @_;
    return $pair  if ($r1 eq $r2);
    return $mispair;
}


# CA_neighbor()
#
sub CA_neighbor {
    my ($CA_res1, $CA_res2) = @_;
    my $i;

    # if not defined at this position, don't penalize
    #   ORIG BUGGY: return undef  if (! defined $CA_res1 || ! defined $CA_res2);
    return 'true'  if (! defined $CA_res1 || ! defined $CA_res2);

    my $x1 = $CA_res1->[0];
    my $x2 = $CA_res2->[0];
    my $y1 = $CA_res1->[1];
    my $y2 = $CA_res2->[1];
    my $z1 = $CA_res1->[2];
    my $z2 = $CA_res2->[2];
    my $dx = $x2 - $x1;
    my $dy = $y2 - $y1;
    my $dz = $z2 - $z1;

    return 'true'  if ($dx*$dx + $dy*$dy + $dz*$dz <= 18.0625);  # 4.25 loose vs 3.8
    return undef;
}


# getIdentity()
#
sub getIdentity {
    my ($a, $b) = @_;
    local $score = 0;
    local $align_len = 0;
    local $a_len = $#{$a} + 1;
    local $b_len = $#{$b} + 1;
    local $len = ($a_len <= $b_len) ? $a_len : $b_len;
    &abort ("attempt to compare a zero length segment!")  if ($len == 0);

    for (my $i=0; $i < $len; ++$i) {
        next      if (! $a->[$i] || ! $b->[$i]);
	next      if ($a->[$i] eq ' ' || $b->[$i] eq ' ');
        next      if ($a->[$i] eq '-' || $b->[$i] eq '-');
        ++$score  if ($a->[$i] eq $b->[$i]);
        ++$align_len;
    }
    &abort ("attempt to compare unaligned segment!")  if ($align_len == 0); 
    return $score / $align_len;
}


sub common_modres_map {
    my ($code3) = shift;

    # all of these are supposed to be handled by MODRES
    #   (but aren't always or aren't done properly)
    my %three_to_one = (
			 '5HP' => 'Q',
			 'ABA' => 'C',
			 'AGM' => 'R',
			 'ALY' => 'K',
			 'ASI' => 'D',
			 'BHD' => 'D',
			 'CEA' => 'C',
			 'CGU' => 'E',
			 'CME' => 'C',
			 'CSB' => 'C',
			 'CSE' => 'C',
			 'CSD' => 'C',
			 'CSO' => 'C',
			 'CSP' => 'C',
			 'CSS' => 'C',
			 'CSW' => 'C',
			 'CSX' => 'C',
			 'CXM' => 'M',
			 'CYM' => 'C',
			 'CYG' => 'C',
			 'DOH' => 'D',
			 'FME' => 'M',
			 'GL3' => 'G',
			 'GLH' => 'E',
			 'HYP' => 'P',
			 'ILG' => 'E',
			 'KCX' => 'K',
			 'LLP' => 'K',
			 'LYZ' => 'K',
			 'MEN' => 'N',
			 'MGN' => 'Q',
			 'MHS' => 'H',
			 'MIS' => 'S',
			 'MLY' => 'K',
			 'MSE' => 'M',
			 'NEP' => 'H',
			 'OCS' => 'C',
			 'PCA' => 'Q',
			 'PTR' => 'Y',
			 'SAC' => 'S',
			 'SEP' => 'S',
			 'SMC' => 'C',
			 'STY' => 'Y',
			 'SVA' => 'S',            
			 'TPO' => 'T',
			 'TPQ' => 'Y',
			 'TRN' => 'W',
			 'TRO' => 'W',
			 'YOF' => 'Y',
			);

    return $three_to_one{$code3};
}


sub mapResCode {
    local ($incode, $silent) = @_;
    $incode = uc $incode;
    my $newcode = undef;

    my %one_to_three = ( 'G' => 'GLY',
			 'A' => 'ALA',
			 'V' => 'VAL',
			 'L' => 'LEU',
			 'I' => 'ILE',
			 'P' => 'PRO',
			 'C' => 'CYS',
			 'M' => 'MET',
			 'H' => 'HIS',
			 'F' => 'PHE',
			 'Y' => 'TYR',
			 'W' => 'TRP',
			 'N' => 'ASN',
			 'Q' => 'GLN',
			 'S' => 'SER',
			 'T' => 'THR',
			 'K' => 'LYS',
			 'R' => 'ARG',
			 'D' => 'ASP',
			 'E' => 'GLU',
			 'X' => 'XXX',
			 '0' => '  A',
			 '1' => '  C',
			 '2' => '  G',
			 '3' => '  T',
			 '4' => '  U'
			);

    my %three_to_one = ( 'GLY' => 'G',
			 'ALA' => 'A',
			 'VAL' => 'V',
			 'LEU' => 'L',
			 'ILE' => 'I',
			 'PRO' => 'P',
			 'CYS' => 'C',
			 'MET' => 'M',
			 'HIS' => 'H',
			 'PHE' => 'F',
			 'TYR' => 'Y',
			 'TRP' => 'W',
			 'ASN' => 'N',
			 'GLN' => 'Q',
			 'SER' => 'S',
			 'THR' => 'T',
			 'LYS' => 'K',
			 'ARG' => 'R',
			 'ASP' => 'D',
			 'GLU' => 'E',
			 
			 '  X' => 'X',
			 '  A' => '0',
			 '  C' => '1',
			 '  G' => '2',
			 '  T' => '3',
			 '  U' => '4',
			 ' +A' => '0',
			 ' +C' => '1',
			 ' +G' => '2',
			 ' +T' => '3',
			 ' +U' => '4',
			 
			 # all of these are supposed to be handled by MODRES
			 #   (but aren't always or aren't done properly)
			 '5HP' => 'Q',
			 'ABA' => 'C',
			 'AGM' => 'R',
			 'ALY' => 'K',
			 'ASI' => 'D',
			 'BHD' => 'D',
			 'CEA' => 'C',
			 'CGU' => 'E',
			 'CME' => 'C',
			 'CSB' => 'C',
			 'CSE' => 'C',
			 'CSD' => 'C',
			 'CSO' => 'C',
			 'CSP' => 'C',
			 'CSS' => 'C',
			 'CSW' => 'C',
			 'CSX' => 'C',
			 'CXM' => 'M',
			 'CYM' => 'C',
			 'CYG' => 'C',
			 'DOH' => 'D',
			 'FME' => 'M',
			 'GL3' => 'G',
			 'GLH' => 'E',
			 'HYP' => 'P',
			 'ILG' => 'E',
			 'KCX' => 'K',
			 'LLP' => 'K',
			 'LYZ' => 'K',
			 'MEN' => 'N',
			 'MGN' => 'Q',
			 'MHS' => 'H',
			 'MIS' => 'S',
			 'MLY' => 'K',
			 'MSE' => 'M',
			 'NEP' => 'H',
			 'OCS' => 'C',
			 'PAS' => 'D',
			 'PCA' => 'Q',
			 'PTR' => 'Y',
			 'SAC' => 'S',
			 'SEP' => 'S',
			 'SMC' => 'C',
			 'STY' => 'Y',
			 'SVA' => 'S',            
			 'TPO' => 'T',
			 'TPQ' => 'Y',
			 'TRN' => 'W',
			 'TRO' => 'W',
			 'YOF' => 'Y',
			 
			 '1MG' => 'X',
			 '2DA' => 'X',
			 '2PP' => 'X',
			 '4SC' => 'X',
			 '4SU' => 'X',
			 '5IU' => 'X',
			 '5MC' => 'X',
			 '5MU' => 'X',
			 'ACB' => 'X',
			 'ACE' => 'X',
			 'ACL' => 'X',
			 'ADD' => 'X',
			 'AHO' => 'X',
			 'AIB' => 'X',
			 'ALS' => 'X',
			 'ARM' => 'X',
			 'ASK' => 'X',
			 'ASX' => 'X',          # NOT B, PREFER TOTAL AMBIGUITY
			 'BAL' => 'X',
			 'BE2' => 'X',
			 'CAB' => 'X',
			 'CBX' => 'X',
			 'CBZ' => 'X',
			 'CCC' => 'X',
			 'CHA' => 'X',
			 'CH2' => 'X',
			 'CH3' => 'X',
			 'CHG' => 'X',
			 'CPN' => 'X',
			 'CRO' => 'X',
			 'DAL' => 'X',
			 'DGL' => 'X',
			 'DOC' => 'X',
			 'DPN' => 'X',
			 'EXC' => 'X',
			 'EYS' => 'X',
			 'FGL' => 'X',
			 'FOR' => 'X',
			 'G7M' => 'X',
			 'GLQ' => 'X',
			 'GLX' => 'X',          # NOT Z, PREFER TOTAL AMBIGUITY
			 'GLZ' => 'X',
			 'GTP' => 'X',
			 'H2U' => 'X',
			 'HAC' => 'X',
			 'HEM' => 'X',
			 'HMF' => 'X',
			 'HPB' => 'X',
			 'IAS' => 'X',
			 'IIL' => 'X',
			 'IPN' => 'X',
			 'LAC' => 'X',
			 'LYT' => 'X',
			 'LYW' => 'X',
			 'MAA' => 'X',
			 'MAI' => 'X',
			 'MHO' => 'X',
			 'MLZ' => 'X',
			 'MYR' => 'X',
			 'NAD' => 'X',
			 'NAL' => 'X',
			 'NH2' => 'X',
			 'NIT' => 'X',
			 'NLE' => 'X',
			 'ODS' => 'X',
			 'OXY' => 'X',
			 'PHD' => 'X',
			 'PHL' => 'X',
			 'PNL' => 'X',
			 'PPH' => 'X',
			 'PPL' => 'X',
			 'PRN' => 'X',
			 'PSS' => 'X',
			 'PSU' => 'X',
			 'PVL' => 'X',
			 'PY2' => 'X',
			 'QND' => 'X',
			 'QUO' => 'X',
			 'SEC' => 'X',
			 'SEG' => 'X',
			 'SEM' => 'X',
			 'SET' => 'X',
			 'SIN' => 'X',
			 'SLE' => 'X',
			 'THC' => 'X',
			 'TPN' => 'X',
			 'TRF' => 'X',
			 'UNK' => 'X',
			 'VAS' => 'X',
			 'YRR' => 'X',
			);

    my %fullname_to_one = ( 'GLYCINE'          => 'G',
			    'ALANINE'          => 'A',
			    'VALINE'           => 'V',
			    'LEUCINE'          => 'L',
			    'ISOLEUCINE'       => 'I',
			    'PROLINE'          => 'P',
			    'CYSTEINE'         => 'C',
			    'METHIONINE'       => 'M',
			    'HISTIDINE'        => 'H',
			    'PHENYLALANINE'    => 'F',
			    'TYROSINE'         => 'Y',
			    'TRYPTOPHAN'       => 'W',
			    'ASPARAGINE'       => 'N',
			    'GLUTAMINE'        => 'Q',
			    'SERINE'           => 'S',
			    'THREONINE'        => 'T',
			    'LYSINE'           => 'K',
			    'ARGININE'         => 'R',
			    'ASPARTATE'        => 'D',
			    'GLUTAMATE'        => 'E',
			    'ASPARTIC ACID'    => 'D',
			    'GLUTAMATIC ACID'  => 'E',
			    'ASPARTIC_ACID'    => 'D',
			    'GLUTAMATIC_ACID'  => 'E',
			    'SELENOMETHIONINE' => 'M',
			    'SELENOCYSTEINE'   => 'M',
			    'ADENINE'          => '0',
			    'CYTOSINE'         => '1',
			    'GUANINE'          => '2',
			    'THYMINE'          => '3',
			    'URACIL'           => '4'
			  );
    
    # map it
    #
    if (length $incode == 1) {
	$newcode = $one_to_three{$incode};
    }
    elsif (length $incode == 3) {
        $newcode = $three_to_one{$incode};
    }
    else {
	$newcode = $fullname_to_one{$incode};
    }


    # check for weirdness
    #
    if (! defined $newcode) {
#	&abort ("unknown residue '$incode'");
#	print STDERR ("unknown residue '$incode' (mapping to 'Z')\n");
#	$newcode = 'Z';
	if (! $silent) {
	    print STDERR ("unknown residue '$incode' (mapping to 'X')\n");
	}
	$newcode = 'X';
    }
    elsif ($newcode eq 'X') {
	if (! $silent) {
	    print STDERR ("strange residue '$incode' (seen code, mapping to 'X')\n");
	}
    }

    return $newcode;
}		      


# resAtomRecs()
#
sub resAtomRecs {
    my $incode = shift;
    my @atom_recs = ();
    my $code1 = (length $incode == 1) ? $incode : &mapResCode ($incode);

    my %atomRecs = ( 
'A' => q{
ATOM   atmi  N   ALA  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  ALA  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   ALA  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   ALA  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  ALA  resi       0.000   0.000   0.000 -1.00  0.00
},
'C' => q{
ATOM   atmi  N   CYS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  CYS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   CYS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   CYS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  CYS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  SG  CYS  resi       0.000   0.000   0.000 -1.00  0.00
},
'D' => q{
ATOM   atmi  N   ASP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  ASP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   ASP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   ASP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  ASP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG  ASP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  OD1 ASP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  OD2 ASP  resi       0.000   0.000   0.000 -1.00  0.00
},
'E' => q{
ATOM   atmi  N   GLU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  GLU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   GLU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   GLU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  GLU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG  GLU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD  GLU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  OE1 GLU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  OE2 GLU  resi       0.000   0.000   0.000 -1.00  0.00
},
'F' => q{
ATOM   atmi  N   PHE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  PHE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   PHE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   PHE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  PHE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG  PHE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD1 PHE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD2 PHE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CE1 PHE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CE2 PHE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CZ  PHE  resi       0.000   0.000   0.000 -1.00  0.00
},
'G' => q{
ATOM   atmi  N   GLY  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  GLY  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   GLY  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   GLY  resi       0.000   0.000   0.000 -1.00  0.00
},
'H' => q{
ATOM   atmi  N   HIS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  HIS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   HIS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   HIS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  HIS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG  HIS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  ND1 HIS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD2 HIS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CE1 HIS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  NE2 HIS  resi       0.000   0.000   0.000 -1.00  0.00
},
'I' => q{
ATOM   atmi  N   ILE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  ILE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   ILE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   ILE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  ILE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG1 ILE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG2 ILE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD1 ILE  resi       0.000   0.000   0.000 -1.00  0.00
},
'K' => q{
ATOM   atmi  N   LYS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  LYS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   LYS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   LYS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  LYS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG  LYS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD  LYS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CE  LYS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  NZ  LYS  resi       0.000   0.000   0.000 -1.00  0.00
},
'L' => q{
ATOM   atmi  N   LEU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  LEU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   LEU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   LEU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  LEU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG  LEU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD1 LEU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD2 LEU  resi       0.000   0.000   0.000 -1.00  0.00
},
'M' => q{
ATOM   atmi  N   MET  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  MET  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   MET  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   MET  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  MET  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG  MET  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  SD  MET  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CE  MET  resi       0.000   0.000   0.000 -1.00  0.00
},
'N' => q{
ATOM   atmi  N   ASN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  ASN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   ASN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   ASN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  ASN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG  ASN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  OD1 ASN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  ND2 ASN  resi       0.000   0.000   0.000 -1.00  0.00
},
'P' => q{
ATOM   atmi  N   PRO  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  PRO  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   PRO  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   PRO  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  PRO  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG  PRO  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD  PRO  resi       0.000   0.000   0.000 -1.00  0.00
},
'Q' => q{
ATOM   atmi  N   GLN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  GLN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   GLN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   GLN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  GLN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG  GLN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD  GLN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  OE1 GLN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  NE2 GLN  resi       0.000   0.000   0.000 -1.00  0.00
},
'R' => q{
ATOM   atmi  N   ARG  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  ARG  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   ARG  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   ARG  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  ARG  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG  ARG  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD  ARG  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  NE  ARG  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CZ  ARG  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  NH1 ARG  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  NH2 ARG  resi       0.000   0.000   0.000 -1.00  0.00
},
'S' => q{
ATOM   atmi  N   SER  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  SER  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   SER  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   SER  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  SER  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  OG  SER  resi       0.000   0.000   0.000 -1.00  0.00
},
'T' => q{
ATOM   atmi  N   THR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  THR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   THR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   THR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  THR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  OG1 THR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG2 THR  resi       0.000   0.000   0.000 -1.00  0.00
},
'V' => q{
ATOM   atmi  N   VAL  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  VAL  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   VAL  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   VAL  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  VAL  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG1 VAL  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG2 VAL  resi       0.000   0.000   0.000 -1.00  0.00
},
'W' => q{
ATOM   atmi  N   TRP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  TRP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   TRP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   TRP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  TRP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG  TRP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD1 TRP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD2 TRP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  NE1 TRP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CE2 TRP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CE3 TRP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CZ2 TRP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CZ3 TRP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CH2 TRP  resi       0.000   0.000   0.000 -1.00  0.00
},
'Y' => q{
ATOM   atmi  N   TYR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  TYR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   TYR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   TYR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  TYR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG  TYR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD1 TYR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD2 TYR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CE1 TYR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CE2 TYR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CZ  TYR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  OH  TYR  resi       0.000   0.000   0.000 -1.00  0.00
},
'X' => q{
ATOM   atmi  N   UNK  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  UNK  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   UNK  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   UNK  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  UNK  resi       0.000   0.000   0.000 -1.00  0.00
}
    );
    
    my $lines = $atomRecs{$code1};
    if (defined $lines) {
	$lines =~ s/^\s+|\s+$//g;
	@atom_recs =  split (/\n/, $lines);
    }
    
    return @atom_recs;
}


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
\t -pdbfile    <pdbfile>
\t[-chain      <chain>]
\t[-outfile    <outfile>]
\t[-fastain    <fastain>]    
\t[-fastaout   <fastaout>]    
\t[-resmapout  <resmapout>]    
};
#\t[-maxpeptidebond  <maxpeptidebond>]  (def: 1.7 Ang)

    # Get args
    #
    local %opts = ();
#    &GetOptions (\%opts, "pdbfile=s", "outfile=s", "chain=s", "maxpeptidebond=f");
    &GetOptions (\%opts, "pdbfile=s", "outfile=s", "chain=s", "fastaout=s", "fastain=s", "resmapout=s");


    # Check for legal invocation
    #
    if (! defined $opts{pdbfile}) {
        print STDERR "$usage\n";
        exit -1;
    }
    &checkExistence ('f', $opts{pdbfile});	


    # Defaults
    #
    #$opts{chain} = '_'  if (! defined $opts{chain} || $opts{chain} eq '0');
    $opts{chain} = '_'  if (! defined $opts{chain});
    $opts{chain} = uc $opts{chain};
#    $opts{maxpeptidebond} = 1.7  if (! defined $opts{maxpeptidebond});

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
