#!/usr/bin/perl -w
#  CVS information:
#  $Revision: 4052 $
#  $Date: 2004-02-26 16:54:18 -0500 (Thu, 26 Feb 2004) $
#  $Author: jeff $
# sorts nnmake_loop library down
# sorts to find best 100 loop closers
# sorts these by meta-score rank and takes top 20
# prints these
# then follows this with up to 40 more ranked soley by meta_score
# duplicates are dropped.
# thus top 20 are best loop closers, next 40 are best meta_score.


if (@ARGV < 1) {
    print STDERR "usage: $0 loop_library\n";
    print STDERR "example: $0 2hntTest.loops_all > 2hntTest.loops\n";
    print STDERR "this script trims a loop library for use in rosetta loop mode\n";
    exit;
}


# nnmake writes out: dme_score,rms_score,rms_score2,seq_score,template_ss_score
#                    ss_score1, ss_score2, ss_score3, ss_type v_point,frag_name
# nnmake ranks by:
#               rank =seq_score+0.5*(ss_score2+ss_score3)+template_ss_score+  0.05*rms_score
#               loop<13
#                   rms<1.4          add 0.05*rms                  
#                   1.4< rms <2.4    add 1.05*rms
#                   2.4<rms          add 2.05*rms
#               and if loop<5        add another 1*rms
# then in rosetta, (CASP) throw out loops with rms >0.5*length for length <=7
# now, for loops < 12discard ones can't close to <0.5rms by small moves
$| = 1;
#$TOPN = 49;
#$TOPX = 249;
$GAP = 5;  #rms
$CLOSE = 6; #rms2
$SEQ = 7;   
$JUFO = 9;
$JONES = 10;
$SAM = 11;
$RANK = $SAM+1;  #line number from nnmake file (ie rank according to nnmake)
$FINAL =$RANK+1; #entire line 
$loop_len = 0;
$last_loop_start = 0;
$last_loop_stop  = 0;

# read lines for this loop region into X
@lines = &fileBufArray ($ARGV[0]);
for ($line_i=0; $line_i <= $#lines; ++$line_i) {
    $line = $lines[$line_i];
    $i |=0 ; push @X,$line ; $i++; 
    
    $loop_len   = substr ($line, 0, 3);
    $loop_start = substr ($line, 3, 5);
    $loop_stop  = substr ($line, 8, 5);
    $rest       = substr ($line, 13);
    if ($rest) {
	$last_loop_start = $loop_start if ($last_loop_start == 0);
	$last_loop_stop  = $loop_stop  if ($last_loop_stop  == 0);
    }

    if (($line_i == $#lines || 
	     $loop_start != $last_loop_start || 
	     $loop_stop != $last_loop_stop) && 
          ($last_loop_start != 0 || $last_loop_stop != 0))  {


	pop @X;
	
	$TOPX = ($i < 250)   ? $i-1       : 249;
	$TOPN = ($i < 250)   ? int ($i/5) : 49;
	$TOPN = ($TOPN < 25) ? $TOPX      : $TOPN;
	warn "$last_loop_start $last_loop_stop I: $i-1 TOPN: $TOPN TOPX: $TOPX\n";
	
        # for each loop line, make @x comprised of first elements of loop line
	# (up to sam score), orig line number (RANK), entire line (FINAL)
        # put all these x arrays into array H
	$i = 0;
	@H = map { @x=split; [ @x[0..$SAM],$i++,$_] } @X;
	
        # Make a list sorted by rms2 score (CLOSE)
        # keep top TOPN of these (but ordered by RANK)
        # of top TOPX best CLOSE, keep top TOPN JUFO, top TOPN SAM and top TOPN JONES (removing redundancies):

        # sort H on the CLOSE score, store in loop_closers
	@loop_closers = sort { $a->[$CLOSE] <=> $b->[$CLOSE] } @H; 
        # take best TOPN CLOSE scores, sort on RANK and put in B 
	@B = sort { $a->[$RANK] <=> $b->[$RANK] } @loop_closers[0..$TOPN]; 
        # print these, and hash them  (ie key is RANK, value is 1 (stored))  
	print  map{$hash{ $_->[$RANK] } = 1  ; $_->[$FINAL]." loop+rank\n" } @B; 

        # sort top TOPX CLOSE scores by JUFO, store in B
	@B = sort { $a->[$JUFO] <=> $b->[$JUFO] } @loop_closers[0..$TOPX]; 
        # of the top TOPN of these, put the non redundant ones in I
	@I = grep { not $hash{ $_->[$RANK]   } } @B[0..$TOPN]; 
        # print them and hash   
	print  map { $hash{ $_->[$RANK] } = 1  ; $_->[$FINAL] ." loop+jufo\n";} @I if @I; 
   
        # repeat for SAM
	@B = sort { $a->[$SAM] <=> $b->[$SAM] } @loop_closers[0..$TOPX]; 
	@I = grep { not  $hash{ $_->[$RANK]  }   } @B[0..$TOPN]; 
	print  map { $hash{ $_->[$RANK] } = 1  ; $_->[$FINAL]." loop+SAM\n"} @I if @I; 
	
        # repeat for JONES
	@B = sort { $a->[$JONES] <=> $b->[$JONES] } @loop_closers[0..$TOPX]; 
	@I = grep {not  $hash{ $_->[$RANK]  }  } @B[0..$TOPN]; 
	print  map { $hash{ $_->[$RANK] } = 1  ; $_->[$FINAL]." loop+JONES\n"; } @I if @I; 
	
	
        # Make a list sorted by JUFO+SEQ score
        # keep top TOPN
        # of to TOPX,  keep top TOPN GAP; 
   
        # B=sorted by JUFO+SEQ; print top TOPN non-redundant and hash
	@B = sort { $a->[$SEQ]+$a->[$JUFO] <=> $b->[$SEQ]+$b->[$JUFO] } @H;
	print  map { $hash{ $_->[$RANK] } = 1  ; $_->[$FINAL] ." jufo\n";} grep { not defined ($hash{ $_->[$RANK] } )  } @B[0..$TOPN]; 
	
        # B=top TOPX SEQ+JUFO sorted by GAP; print top TOPN non-redundant & hash
	@B = sort { $a->[$GAP] <=> $b->[$GAP] } @B[0..$TOPX]; 
	print  map { $hash{ $_->[$RANK] } = 1  ; $_->[$FINAL] ." jufo+loop\n";} grep { not defined ($hash{ $_->[$RANK] } )  } @B[0..$TOPN]; 
	
        # repeat for JONES+SEQ     
	@B = sort { $a->[$SEQ]+$a->[$JONES] <=> $b->[$SEQ]+$b->[$JONES] } @H;
	print  map { $hash{ $_->[$RANK] } = 1  ; $_->[$FINAL] ." jones\n";} grep { not defined ($hash{ $_->[$RANK] } )  } @B[0..$TOPN]; 
	@B = sort { $a->[$GAP] <=> $b->[$GAP] } @loop_closers[0..$TOPX]; 
	print  map { $hash{ $_->[$RANK] } = 1  ; $_->[$FINAL] ." jones+loop\n";} grep { not defined ($hash{ $_->[$RANK] } )  } @B[0..$TOPN]; 
	
        # repeat for SAM+SEQ   
	@B = sort { $a->[$SEQ]+$a->[$SAM] <=> $b->[$SEQ]+$b->[$SAM] } @H;
	print  map { $hash{ $_->[$RANK] } = 1  ; $_->[$FINAL] ." sam\n";} grep { not defined ($hash{ $_->[$RANK] } )  } @B[0..$TOPN]; 
	@B = sort { $a->[$GAP] <=> $b->[$GAP] } @loop_closers[0..$TOPX]; 
	print  map { $hash{ $_->[$RANK] } = 1  ; $_->[$FINAL] ." sam+loop\n";} grep { not defined ($hash{ $_->[$RANK] } )  } @B[0..$TOPN]; 
	
	
        # reinit variables to prepare for next block of 2000; 
	
	$i=1; 
	@X=($line);
	%hash = ();
	$last_loop_start = $loop_start;   
	$last_loop_stop  = $loop_stop ;   
    }
    if (! $rest) {  # no loop data for this one
	warn "skipping loop $loop_len $loop_start $loop_stop\n";
	$i=0; 
	@X=();
	%hash = ();
	print $line."\n";
	$last_loop_start = 0;   
	$last_loop_stop  = 0;   
	next;
    }
}

###############################################################################
# util
###############################################################################

# abort()
#
sub abort {
    my $msg = shift;
    print STDERR "$0: $msg\n";
    exit -2;
}
                 
# fileBufString()
#
sub fileBufString {
    my $file = shift;
    my $oldsep = $/;
    undef $/;
    if ($file =~ /\.gz|\.Z/) {
        if (! open (FILE, "gzip -dc $file |")) {
            &abort ("$0: unable to open file $file for gzip -dc");
        }
    }
    elsif (! open (FILE, $file)) {
        &abort ("$0: unable to open file $file for reading");
    }
    my $buf = <FILE>;
    close (FILE);
    $/ = $oldsep;
    return $buf;
}

# fileBufArray()
#
sub fileBufArray {
    my $file = shift;
    my $oldsep = $/;
    undef $/;
    if ($file =~ /\.gz|\.Z/) {
        if (! open (FILE, "gzip -dc $file |")) {
            &abort ("$0: unable to open file $file for gzip -dc");
        }
    }
    elsif (! open (FILE, $file)) {
        &abort ("$0: unable to open file $file for reading");
    }
    my $buf = <FILE>;
    close (FILE);
    $/ = $oldsep;
    @buf = split (/$oldsep/, $buf);
    pop (@buf)  if ($buf[$#buf] eq '');
    return @buf;
}

# bigFileBufArray()
#
sub bigFileBufArray {
    my $file = shift;
    my $buf = +[];
    if ($file =~ /\.gz|\.Z/) {
        if (! open (FILE, "gzip -dc $file |")) {
            &abort ("$0: unable to open file $file for gzip -dc");
        }
    }
    elsif (! open (FILE, $file)) {
        &abort ("$0: unable to open file $file for reading");
    }
    while (<FILE>) {
        chomp;
        push (@$buf, $_);
    }
    close (FILE);
    return $buf;
}

###############################################################################
# end
1;                                                     # in case it's a package
###############################################################################
