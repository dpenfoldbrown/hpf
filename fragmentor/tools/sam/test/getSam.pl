#!/users/robetta/bin/perl

###############################################################################
# conf
###############################################################################

# paths
$srcdir       = $ENV{'HOME'}."/src";
$samdir       = "$srcdir/shareware/sam-t99";
$samt99dir    = "$srcdir/shareware/sam";
$samdir       =~ s!//!/!g;
$samt99dir    =~ s!//!/!g;
$target99     = "$samt99dir/bin/target99";
$uniqueseq    = "$sam99dir/bin/uniqueseq";
$uniqueseq_percent_id = 0.9;
$predict_2nd  = "$samdir/predict-2nd";
$condense_rdb = "$samdir/condense_rdb.pl";


###############################################################################
# init
###############################################################################

%opts = &getCommandLineOptions ();
$fastafile = $opts{'fastafile'};
$outdir    = $opts{'outdir'};

$outdir =~ s!//!/!g;

$id = $fastafile;
$id =~ s!\.[^\.]*$!!;
$id =~ s!^.*/!!;

###############################################################################
# main
###############################################################################

# sam dies if fasta is missing header line!
$tmp_fastafile = "$outdir/$id-tmp.fasta";
$tmp_id = "$outdir/$id-tmp";
$tmp_uniqueseq_a2m_id = "$id-uniqueseq-tmp";
$tmp_uniqueseq_a2mfile = "$outdir/$id-uniqueseq-tmp.a2m";
$tmp_a2mfile = "$outdir/$id-tmp.a2m";

@tmp_fasta_buf = &fileBufArray ($fastafile);
if ($tmp_fasta_buf[0] !~ /^\s*\>/) {
    unshift (@tmp_fasta_buf, ">$id");
}
open (TMP_FASTA, '>'.$tmp_fastafile);
print TMP_FASTA join ("\n", @tmp_fasta_buf);
close (TMP_FASTA);


# don't want to be in $outdir, rather in $samdir
#chdir $outdir;

# get into $samdir so sam can read file recode3.20comp
chdir $samdir;                                

# create samscript
$sam_6state = "$outdir/$id.sam_6state";
$sam_ebghtl = "$outdir/$id.sam_ebghtl";
$sam_log    = "$outdir/$id.sam_log";
$samscript_txt_buf = 
qq{ReadAlphabet $samdir/std.alphabet
ReadAlphabet $samdir/DSSP.alphabet
ReadNeuralNet $samdir/overrep-3617-IDaa13-7-10-11-10-11-7-7-ebghtl-seeded3-stride-trained.net
ReadA2M $tmp_uniqueseq_a2mfile
PrintRDB $sam_6state
PrintPredictionFasta $sam_ebghtl
};
$samscript_txt = "$outdir/$id.samscript.txt";
open  (SAMSCRIPT, '>'.$samscript_txt);
print  SAMSCRIPT $samscript_txt_buf;
close (SAMSCRIPT);

# run sam

&runCmd (qq{$target99 -seed $tmp_fastafile -out $tmp_id}); 
&runCmd (qq{$uniqueseq $tmp_uniqueseq_a2m_id -percent_id 0.9 -alignfile $tmp_a2mfile});
&runCmd (qq{$predict_2nd -noalph < $samscript_txt >& $sam_log});
&runCmd (qq{$condense_rdb $sam_6state});

# clean up
#unlink $sam_6state;                                  # might want to keep this
#unlink $sam_ebghtl;                                  # might want to keep this
#unlink $tmp_a2mfile;                                  # might want to keep this
unlink $samscript_txt;
unlink $sam_log;
unlink $tmp_fastafile;

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
    local $usage = qq{usage: $0 
\t -fastafile  <fastafile>
\t[-outdir     <outdir>]    (def: .)
};

    # Get args
    #
    local %opts = ();
    &GetOptions (\%opts, "fastafile=s", "outdir=s");

    # Check for legal invocation
    #
    if (!defined $opts{fastafile}) {
        print STDERR "$usage\n";
        exit -1;
    }

    # defaults
    #
    $opts{outdir} = '.'  if (! defined $opts{outdir});


    $opts{outdir}   =~ s!/$!!;
    &checkExist          ('f', $opts{fastafile});
    &checkExistAndCreate ('d', $opts{outdir});


    return %opts;
}

###############################################################################
# util
###############################################################################

# runCmd()
#
sub runCmd {
    local $cmd = shift;
    local $retcode = 0;
#    print "RUN: '$cmd'\n\n";
    if (system ($cmd) != 0) {
#	&abort ("FAIL: $cmd");
	print STDERR ("FAIL: $cmd");
	$retcode = $?>>8;
    }
    return $retcode;
}

# checkExist()
#
sub checkExist {
    local ($type, $path) = @_;
    if ($type eq 'f' && ! -f $path) {
	&abort ("filenotfound $path");
	
    }
    if ($type eq 'd' && ! -d $path) { 
	&abort ("dirnotfound $path");
    }
    return 'true';
}

# checkExistAndCreate()
#
sub checkExistAndCreate {
    local ($type, $path) = @_;
    if ($type eq 'f' && ! -f $path) {
	#print "creating $path...\n";
	open (FILE, '>'.$path);
	close (FILE);
    }
    if ($type eq 'd' && ! -d $path) { 
	#print "creating $path...\n";
	$mode = 0777  if (! $mode);
	mkdir ($path, $mode);
    }
    return 'true';
}

# abort()
#
sub abort {
    local $msg = shift;
    print STDERR "$0: $msg\n";
    exit -2;
}

# fileBufString()
#
sub fileBufString {
    local $file = shift;
    local $oldsep = $/;
    undef $/;
    if (! open (FILE, $file)) {
	&abort ("$0: unable to open file $file for reading");
    }
    local $buf = <FILE>;
    close (FILE);
    $/ = $oldsep;
    return $buf;
}

# fileBufArray()
#
sub fileBufArray {
    local $file = shift;
    local $oldsep = $/;
    undef $/;
    if (! open (FILE, $file)) {
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
# end
###############################################################################
