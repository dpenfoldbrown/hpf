#!/usr/bin/perl
##
## Copyright 2002, University of Washington, the Baker Lab, and Dylan Chivian.
##   This document contains private and confidential information and its 
##   disclosure does not constitute publication.  All rights are reserved by 
##   University of Washington, the Baker Lab, and Dylan Chivian, except those 
##   specifically granted by license.
##
##  Initial Author: Dylan Chivian (dylan@lazy8.com)
##  $Revision: 1.3 $
##  $Date: 2002/05/15 18:36:37 $
##  $Author: rohl $
##
###############################################################################


###############################################################################
# init
###############################################################################

%opts       = &getCommandLineOptions ();
$fastafile  = $opts{fastafile};
$outfile    = $opts{outfile};
$label      = $opts{label};
$method     = $opts{method};
$email      = $opts{email};
$mbox       = $opts{mbox};
$parse_mbox = $opts{parsembox};
$patience   = $opts{patience};

$timestamp = `date +%Y-%m-%d_%T`;  chomp $timestamp;
$host = $ENV{HOST};
$pid = $$;
$runid = $timestamp.'.'.$host.'.'.$pid;
$runid = ($label) ? $label.'.'.$runid : $runid;

###############################################################################
# conf
###############################################################################

#$| = 1;  # disable STDOUT buffering

# urls
#$psipred_server = "http://insulin.brunel.ac.uk/cgi-bin/psipred/psipred.cgi";
#$phd_server     = "http://dodo.cpmc.columbia.edu/cgi/pp/submit";
$psipred_server = "http://bioinf.cs.ucl.ac.uk/cgi-bin/psipred/psipred.cgi";
$phd_server     = "http://cubic.bioc.columbia.edu/cgi/pp/submit";
$sam_server     = "http://www.cse.ucsc.edu/~farmer/cgi-bin/T99-query.pl";

# paths
$tmp   = "./";

###############################################################################
# main
###############################################################################

$fasta = &getFastaFromFile ($fastafile);
$fasta =~ s/\s+//g;

&submitSSpred ($fasta, $method, $runid);

if ($parse_mbox) {
    &parseSSpredFromMbox ($mbox, $method, $runid, $outfile, $patience);
}

exit 0;

###############################################################################
# subs
###############################################################################

# getFastaFromFile ()
#
sub getFastaFromFile {
    my $fastafile = shift;
    my $fasta = '';
    my $line = '';
    foreach $line (&fileBufArray ($fastafile)) {
	next if ($line =~ /^\s*\>/);
	$fasta .= "$line\n";
    }
    $fasta =~ s/^\s+|\s+$//g;
    return $fasta;
}


# submitSSpred ()
#
sub submitSSpred {
    my ($fasta, $method) = @_;
    my $url      = '';
    my $get_args = '';

    if ($method eq 'psipred') {
        $url         = $psipred_server;
        my $Email    = &escapeGetArg ($email);
        my $Sequence = &escapeGetArg ($fasta);
        my $Subject  = &escapeGetArg ($runid);
        my $Program  = 'psipred';
        my $Output   = 'opnone';

        $get_args = qq{
	    Email    = $Email \\&
	    Sequence = $Sequence \\&
	    Subject  = $Subject \\&
	    Program  = $Program \\&
	    Output   = $Output
	};
	$get_args =~ s/\s+//g;
    }

    elsif ($method eq 'phd') {
        $url           = $phd_server;
        #$url           = "http://shampoo.baker:8008/";
        my $usr_email  = &escapeGetArg ($email);
        my $sequence   = &escapeGetArg ($fasta);
        my $seq_name   = &escapeGetArg ($runid);
	my $opt_prd    = &escapeGetArg ('secondary structure only (PHDsec)');
	my $opt_ali    = &escapeGetArg ('no alignment returned');
	my $seq_format = &escapeGetArg ('Default= single sequence - click to select other');
	my $resp_mode  = &escapeGetArg ('batch (results returned by email, if address correct!)');
	my $usr_pwd = '';
	my $opt_ret_concise = '';
	my $sequence_file = '';

        $post_args = qq{
	    usr-email       = $usr_email \&
	    usr-pwd         = $usr_pwd \&
	    opt-prd         = $opt_prd \&
	    opt-ali         = $opt_ali \&
	    seq-name        = $seq_name \&
	    seq-format      = $seq_format \&
	    sequence        = $sequence \&
	    resp-mode       = $resp_mode
        };
	$post_args =~ s/\s+//g;
	$post_args .= '&';
	$post_args .= qq{sequenceFile%22%3b+filename%3d%22=};
    }

    elsif ($method eq 'sam') {
        $url            = $sam_server;
        my $address     = &escapeGetArg ($email);
        my $subjectline = &escapeGetArg ($runid);
        my $sequence    = &escapeGetArg ($fasta);
	my $db          = &escapeGetArg ('PDB and SCOP domains');
	my $RetTargetAlignmentA2M         = 'on';
	my $RetTargetAlignmentPrettyAlign = 'on';
	my $RetTargetAlignmentHtml        = 'on';
	my $SearchDb                      = 'on';
	my $ScoreMethod                   = 'Average';
	my $evalue                        = '1.0';
	my $numalign                      = '40';
	my $PairwisePrettyAlign           = 'on';
	my $PairwiseA2m                   = 'on';
	my $RetSecStrRdb                  = 'on';
	my $RetSecStrFasta                = 'on';

        $get_args = qq{
	    address                       = $address \\&
	    subjectline                   = $subjectline \\&
	    sequence                      = $sequence \\&
	    db                            = $db \\&
	    RetTargetAlignmentA2M         = $RetTargetAlignmentA2M \\&
	    RetTargetAlignmentPrettyAlign = $RetTargetAlignmentPrettyAlign \\&
	    RetTargetAlignmentHtml        = $RetTargetAlignmentHtml \\&
	    SearchDb                      = $SearchDb \\&
	    ScoreMethod                   = $ScoreMethod \\&
	    evalue                        = $evalue \\&
	    numalign                      = $numalign \\&
	    PairwisePrettyAlign           = $PairwisePrettyAlign \\&
	    PairwiseA2m                   = $PairwiseA2m \\&
	    RetSecStrRdb                  = $RetSecStrRdb \\&
	    RetSecStrFasta                = $RetSecStrFasta
	};
	$get_args =~ s/\s+//g;
    }

    else {
        &abort ("unknown ss pred method '$method'");
    }

    # submit it
    #
    my $ret_html = '';
    if ($method eq 'phd') {
	$postargs_file = "$tmp/$runid.post";
	open  (POSTARGS, '>'.$postargs_file);
	print  POSTARGS $post_args;
	close (POSTARGS);
	# multipart broken
	#$ret_html = `$httppost $url $postargs_file multipart`;
	#$ret_html = `$httppost $url $postargs_file flat`;
	$ret_html = &httppost ($url, $postargs_file, "flat");
	unlink ($postargs_file);
	# debug
	#print STDOUT "$RET_HTML: '$ret_html'\n";
    } else {
	#$ret_html = `$httpget $url\\?$get_args`;
	$ret_html = &httpget ("$url\\?$get_args");
    }
    #if ($? != 0) {
    if ($ret_html =~ /server error/i) {
	&abort ("failure getting ss prediction from $method server.\nreturn message: '$ret_html'");
    }

    return;
}


# parseSSpredFromMbox ()
#
sub parseSSpredFromMbox {
    my ($mbox, $method, $runid, $outfile, $patience) = @_;
    my $found       = undef;
    my @mbox_buf    = ();
    my $psipred_msg = '';
    my $interval = 600;   # sec, so 10 min
    my $time_spent  = 0;
    my $i=0, $j=0, $start_i=0;

    #if ($method eq 'psipred' || $method eq 'phd') {
    #if ($method eq 'psipred') {
    #$patience = $short_patience;
    #$interval = 30;     # 30 sec
    #} elsif ($method eq 'sam' || $method eq 'phd') {
    #$patience = $long_patience;
    #$interval = 30;     # 30 sec
    #} else {
    #&abort ("unknown method '$method'");
    #}

    while ($time_spent < $patience) {
	@mbox_buf = &fileBufArray ($mbox);
	for ($i=0; $i <= $#mbox_buf; ++$i) {
	    if ($method eq 'psipred') {
		$found = 'true' if ($mbox_buf[$i] =~ /^Subject\:\s+$runid\s+PSIPRED Results/);
	    } 
	    elsif ($method eq 'phd') {
		$found = 'true' if ($mbox_buf[$i] =~ /description=$runid/);
	    }
	    elsif ($method eq 'sam') {
		$found = 'true' if ($mbox_buf[$i] =~ /^Subject:\s*$runid\s*/);
	    }
	    else {
		&abort ("unknown method '$method'");
	    }
	    next if (! $found);

	    for ($j=$i-1; $j >= 0 && $mbox_buf[$j] !~ /^From\s+/i; --$j) {
		$start_i = $j-1;
	    }
	    $sspred_msg = "$mbox_buf[$start_i]\n";
	    for ($j=$start_i+1; $j <= $#mbox_buf && $mbox_buf[$j] !~ /^From\s+/i; ++$j) {
		$sspred_msg .= "$mbox_buf[$j]\n";
	    }
	    last;
	}
	last if ($found);

	$time_spent += $interval;
	sleep $interval;
    }

    # write it (or if sam, retrieve from server)
    #
    if (! $found) {
	&abort ("patience $patience secs exceeded waiting for response");
    }
    if ($method eq 'psipred' || $method eq 'phd') {
	&writeBuf ($sspred_msg, $outfile);
    }
    elsif ($method eq 'sam') {
	$pickup_url = $sspred_msg;
	$pickup_url =~ s/.*(http\S+).*/$1/s;
	$pickup_url =~ s/index.html$/secstr.rdb.txt/;

	#$rdb_buf = `$httpget $pickup_url`;
	$rdb_buf = &httpget ($pickup_url);
	&writeBuf ($rdb_buf, $outfile);
    }

    return;
}


# httpget()
#
sub httpget {
    my $showhead = 0;
    my $ret_msg = '';
    ($host, $port, $uri) = &getGetParams (@_);
    ($stat, $head, $page) = &SiteSucker::makeHttpRequest ($host, $port, $uri); 
    
    if ($showhead) {
	foreach $k (sort keys %$head) {
	    $ret_msg .= $head->{$k}."\n";
	}
	$ret_msg .= "\n\n";
    }

    $ret_msg .= "$page\n";
    return $ret_msg;
}


# httppost()
#
sub httppost {
    my $showhead = 0;
    my $ret_msg = '';
    ($host, $port, $uri, $postfile, $flat_or_multipart) = &getPostParams (@_);
    $post = &fileBufString ($postfile);
    ($stat, $head, $page) = &SiteSucker::makeHttpRequestPost ($host, $port, $uri, $post, $flat_or_multipart); 
    
    if ($showhead) {
	foreach $k (sort keys %$head) {
	    $ret_msg .= $head->{$k}."\n";
	}
	$ret_msg .= "\n\n";
    }
    $ret_msg .= "$page\n";
    return $ret_msg;
}


# getGetParams()
#
sub getGetParams {
    my @argv = @_;
    my ($host, $port, $uri);

    if ($#argv < 0) {
        print "usage: $0 <url>\n";
        exit 0;
    }
    $url = $argv[0];;

    $url =~ s!^http://!!;
    if ($url =~ m!^([^:/]+):?(\d*)(.*)!) {
        $host = $1;
        $port = $2;
        $uri = $3;
    }
    else {
        print STDERR "$0: malformed url '$url'\n";
        exit -2;
    }

    $port = 80  if (! $port);
    $uri = '/'  if (! $uri);

    #$uri .= '/'  if ($uri !~ /\/$/ && $uri !~ /\./ && $uri !~ /\?/);

    return ($host, $port, $uri);
}


# getPostParams
#
sub getPostParams {
    my @argv = @_;
    my ($host, $port, $uri);

    if ($#argv < 0) {
	print "usage: $0 <url> <postfile> <flat_or_multipart>\n";
	exit 0;
    }
    $url               = $argv[0];
    $postfile          = $argv[1];
    $flat_or_multipart = $argv[2];

    $url =~ s!^http://!!;
    if ($url =~ m!^([^:/]+):?(\d*)(.*)!) {
	$host = $1;
	$port = $2;
	$uri  = $3;
    }
    else {
	print STDERR "$0: malformed url '$url'\n";
	exit -2;
    }

    $port = 80  if (! $port);
    $uri = '/'  if (! $uri);

    #$uri .= '/'  if ($uri !~ /\/$/ && $uri !~ /\./ && $uri !~ /\?/);

    return ($host, $port, $uri, $postfile, $flat_or_multipart);
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
\t -fastafile <fastafile>
\t[-method    <method>]  (psipred, phd, sam) (def: psipred)
\t[-label     <label>]                       (def: date_stamp)
\t[-outfile   <outfile>]                     (def: STDOUT)
\t[-email     <email_addr>                   (def: $ENV{USER}\@mail.bakerlab.org)
\t[-mbox      <mbox_spool>                   (def: $ENV{MAIL})
\t[-parsembox <T/F>]                         (def: F)
\t[-patience  <patience_with_mbox>]  (sec)   (def: psipred: 1800, other: 86400)
};

    # Get args
    #
    local %opts = ();
    &GetOptions (\%opts, 
		 "fastafile=s",
		 "method=s",
		 "label=s",
		 "outfile=s",
		 "email=s",
		 "mbox=s",
		 "parsembox=s",
		 "patience=i");

    # Check for legal invocation
    #
    if (! defined $opts{fastafile}) {
        print STDERR "$usage\n";
        exit -1;
    }

    # defaults
    #
    $opts{method}    = 'psipred'                           if (! defined $opts{method});
    $opts{parsembox} = undef                               if ($opts{parsembox} =~ /^f/i);
    $opts{email}     = $ENV{'USER'}."\@mail.bakerlab.org"  if (! $opts{email});
    $opts{mbox}      = $ENV{'MAIL'}                        if (! $opts{mbox});

    if ($opts{method} eq 'psipred') {
	$opts{patience} = 1800  if (! defined $opts{patience});    # 30 min
    } elsif ($opts{method} eq 'phd' || $opts{method} eq 'sam') {
	#$long_patience  = 36000  if (! defined $opts{patience});  # 10 hours
	$opts{patience} = 86400  if (! defined $opts{patience});   # 24 hours
    } else {
	&abort ("unknown method '".$opts{method}."'");
    }

    # existence checks
    &checkExist ('f', $opts{fastafile});
    &checkExist ('f', $opts{mbox});

    return %opts;
}
# end getCommandLineOptions()

###############################################################################
# util
###############################################################################

# maxInt()
#
sub maxInt {
    local ($v1, $v2) = @_;
    return ($v1 > $v2) ? $v1 : $v2;
}
# end maxInt()


# makeGetArgs (\%a)
#
#
#   desc:  make a string suitable for GET args out of associative array
#
#   args:  \%a    array with key->val pairs
#
#   rets:  $ret   string for GET args
#
#
sub makeGetArgs {
    local ($a) = @_;
    my ($k, $s, $ret);

    foreach $k (keys %$a) {
        $s = $a->{$k};
        $s = &escapeGetArg ($s);
        $ret .= "$k=$s&";
    }
    chop $ret;                                      # Remove trailing ampersand
    return $ret;
}
# end makeGetArgs ()


# escapeGetArg ()
#
sub escapeGetArg {
    local $str = shift;

    $str =~ s/ /+/go;
    $str =~ s/\%/\%25/go;
    $str =~ s/([^0-9a-zA-Z\%+])/"\%".&charToHex($1)/ge;

    return $str;
}
# end escapeGetArg ()


# charToHex ()
#
sub charToHex {
    my $ascii = ord($_[0]);
    my %hexMap = (  0 => '0',
                    1 => '1',
                    2 => '2',
                    3 => '3',
                    4 => '4',
                    5 => '5',
                    6 => '6',
                    7 => '7',
                    8 => '8',
                    9 => '9',
                   10 => 'a',
                   11 => 'b',
                   12 => 'c',
                   13 => 'd',
                   14 => 'e',
                   15 => 'f'
		    );

    return $hexMap{(($ascii & 0xf0) >> 4)} . $hexMap{($ascii & 0x0f)};
}
# end charToHex ()


#  hexToChar ()
#
sub hexToChar {
    my $ascii = hex($_[0]);
    return chr $ascii;
}
# end hexToChar ()


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
        print "creating $path...\n";
        open (FILE, '>'.$path);
        close (FILE);
    }
    if ($type eq 'd' && ! -d $path) {
        print "creating $path...\n";
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
#    print "$0: $msg\n";
    exit -2;
}


# writeBuf ()
#
sub writeBuf {
    my ($buf, $outfile) = @_;
    if ($outfile) {
        #print "creating $outfile\n";
        open (OUTFILE, '>'.$outfile);
        select (OUTFILE);
    }
    print $buf;
    if ($outfile) {
        close (OUTFILE);
        select (STDOUT);
    }
    return;
}


# fileBufString ()
#
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


# fileBufArray ()
#
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
package SiteSucker;
###############################################################################

use Socket;
use FileHandle;


$accept_cookies = 0;
$keep_sid = 0;
$cookie = '';
$debug = 0;

$| = 1  if ($debug);
%uriVisited = ();


sub crawlURL {
    my ($host, $port, $uri, $refuri, $sid) = @_;
    my ($req_uri, $clean_uri, $cwd_uri) = ('','','');
    my %snd_hdrs = ();
    my $head = undef;
    my $page = '';
    my @links = ();
    my ($new_host, $new_port, $link, $path_elt) = ('','','','');
    my @new_uri = ();
    my $nxt_uri = '';
    my $head_end = -1;
    my $set_cookie = '';
    

    if ($uri !~ m!/cgi-bin/!) {
	$uri =~ s!(/[^\./]+)$!$1/!;
    }
    $clean_uri = $uri;
    $clean_uri =~ s!\#[^/]+$!!;
    $clean_uri =~ s!\?[^/]+$!!;
    return  if ($uriVisited{$clean_uri});
    $uriVisited{$clean_uri} = 1;

    $cwd_uri = $uri;
    $cwd_uri =~ s![^/]*$!!;
    $cwd_uri = '/'  if (! $cwd_uri);


    # Read page
    #
    $req_uri = $uri;
    if ($cookie) {
	$snd_hdrs{Cookie} = $cookie;
    }
    else {
	if (! $accept_cookies && $keep_sid && $sid) {
	    $req_uri =  &insertSession ($uri, $sid) ;
	}
    }
    ($stat, $head, $page) = &makeHttpRequest($host,$port,$req_uri,%snd_hdrs);

    if ($stat != 1) {
	print STDERR "$0: failure: unable to makeHttpRequest()\n";
	exit 2;
    }


    # Play with mime and page
    #
    if (defined &main::examineMime) {
	&main::examineMime ($host, $port, $uri, $refuri, $head);
    }
    if (defined &main::examinePage) {    
	&main::examinePage ($host, $port, $uri, $refuri, $page);
    }


    # Crawl links
    #
#    if (($head_end = index ($page, "</head>")) == -1) {
#	$head_end = index ($page, "</HEAD>");
#    }
#    if ($head_end != -1) {
#	$head_end += length ("</head>");
#	$page = substr ($page, $head_end, length ($page) - $head_end);
#    }
    while ($page =~ s/<a[^<>]+href\s*=\s*\"?([^\"\'\s>]+)[^<>]*>//im) {
	push (@links, $1);
    }
    while ($page =~ s/<frame[^<>]+src\s*=\s*\"?([^\"\'\s>]+)[^<>]*>//im) {
	push (@links, $1);
    }
    undef $page;                                                  # save memory
    undef $head;


    foreach $link (@links) {
	next  if ($link =~ /^mailto:/i);
	next  if ($link =~ /^ftp:/i);
	next  if ($link =~ /^javascript:/i);

	next  if ($link =~ /^\#/);
	$link =~ s/\#.*$//;
	if (! $accept_cookies && $keep_sid) {
	    $sid = &getSession ($link);
	}
	$link = &removeSession ($link);

	if ($link =~ s!^http://([^\:/]+)!!i) {
	    $new_host = $1;
	    next  if ($new_host ne $host);
	}

	if ($link =~ s!^:(\d+)!!) {
	    $new_port = $1;
	    next  if ($new_port ne $port);
	}

	if ($link =~ m!^/!) {
	    $nxt_uri = $link;
	}
	else {
	    $cwd_uri =~ s!^/|/$!!g;
	    @new_uri = split (/\/{1,2}/, $cwd_uri);
	    
	    foreach $path_elt (split (/\/{1,2}/, $link)) {
		next if ($path_elt eq '.' || $path_elt eq '');
		if ($path_elt eq '..') {
		    pop (@new_uri);
		}
		else {
		    push (@new_uri, $path_elt);
		}
	    }
	    $nxt_uri = '/' . join ('/', @new_uri);
	    $nxt_uri .= '/'  if ($new_uri[$#new_uri] !~ /\./);
	}
	&crawlURL ($host, $port, $nxt_uri, $uri, $sid);
    }

    return;
}


sub makeHttpRequest {
    my ($host, $port, $uri, %snd_hdrs) = @_;
    my ($nic);
    my ($name, $aliases, $proto, $type, $len, $addr);
    my $head = undef;
    my $data = '';
    my $status = 1;
    my $linesep = $/;

    
    ($name, $aliases, $proto) = getprotobyname('tcp');
    ($name, $aliases, $type, $len, $addr) = gethostbyname($host);
    $nic = pack('S n N x8', AF_INET, $port, unpack("N", $addr));
    
    if (!socket(S, PF_INET, SOCK_STREAM, $proto)) {
	$status = -1;
	goto done;
    }
    if (!connect(S, $nic)) {
	$status = -2;
	goto done;
    }

    S->autoflush(1);
    $request  = "GET $uri HTTP/1.0\r\n";
    foreach $i (keys %snd_hdrs) {
	next if (! $i);
	$request .= "$i: " . $snd_hdrs{$i} . "\r\n";
    }
    $request .= "\r\n";
    print S $request;


    # Get mime header
    #
    $head->{status} = <S>;
    while (<S>) {
	last  if ($_ =~ /^\s*$/);
	$_ =~ /^([^:]+):\s*(.*)/;
	$head->{$1} = $2;
	if ($accept_cookies && $_ =~ /^set-cookie: (.*)/i) {
	    $set_cookie = $1;
	    $set_cookie =~ s/expires=[^\;]+//;
	    $set_cookie =~ s/path=[^\;]+//;
	    $set_cookie =~ s/domain=[^\;]+//;
	    $set_cookie =~ s/secure//;
	    $set_cookie =~ s/\s*\;\s*//g;
	    $cookie .= " $set_cookie\;";
	}
    }
    $cookie =~ s/^\s+|\s*\;$//g;


    # Get data
    #
    undef $/;
    $data = <S>;
    close (S);
    $/ = $linesep;
    
done:
    return ($status, $head, $data);
}


sub makeHttpRequestPost {
    my ($host, $port, $uri, $post, $flat_or_multipart, %snd_hdrs) = @_;
    my ($nic);
    my ($name, $aliases, $proto, $type, $len, $addr);
    my $head = undef;
    my $data = '';
    my $status = 1;
    my $linesep = $/;

    my $boundary;

    if ($flat_or_multipart !~ /^f/i) { 
	($post, $boundary) = &makeMultiPartPost ($post);
    }

    ($name, $aliases, $proto) = getprotobyname('tcp');
    ($name, $aliases, $type, $len, $addr) = gethostbyname($host);
    $nic = pack('S n N x8', AF_INET, $port, unpack("N", $addr));
    
    if (!socket(S, PF_INET, SOCK_STREAM, $proto)) {
	$status = -1;
	goto done;
    }
    if (!connect(S, $nic)) {
	$status = -2;
	goto done;
    }

    S->autoflush(1);
    #$request  = "POST $uri HTTP/1.0\r\n";
    $request  = "POST $uri HTTP/1.1\r\n";
    foreach $i (keys %snd_hdrs) {
	next if (! $i);
	$request .= "$i: " . $snd_hdrs{$i} . "\r\n";
    }
    my $content_length = length ($post);
    $request .= "Connection: " . "close"          . "\r\n";
    #$request .= "Connection: " . "Keep-Alive"     . "\r\n";
    $request .= "User-Agent: " . "HTTPposter/1.0" . "\r\n";
    $request .= "Host: "       . "$host:$port"    . "\r\n";
    $request .= ($flat_or_multipart !~ /^f/i) 
	        ? "Content-Type: " . "multipart/form-data; boundary=$boundary" 
		: "Content-Type: " . "application/x-www-form-urlencoded";
    $request .= "\r\n"; 
    $request .= "Content-Length: " . $content_length . "\r\n";
    $request .= "\r\n";
    $request .= $post . "\r\n";
    print S $request;

    # debug
    #print STDOUT $request;

    # Get mime header
    #
    $head->{status} = <S>;
    while (<S>) {
	last  if ($_ =~ /^\s*$/);
	$_ =~ /^([^:]+):\s*(.*)/;
	$head->{$1} = $2;
	if ($accept_cookies && $_ =~ /^set-cookie: (.*)/i) {
	    $set_cookie = $1;
	    $set_cookie =~ s/expires=[^\;]+//;
	    $set_cookie =~ s/path=[^\;]+//;
	    $set_cookie =~ s/domain=[^\;]+//;
	    $set_cookie =~ s/secure//;
	    $set_cookie =~ s/\s*\;\s*//g;
	    $cookie .= " $set_cookie\;";
	}
    }
    $cookie =~ s/^\s+|\s*\;$//g;


    # Get data
    #
    undef $/;
    $data = <S>;
    close (S);
    $/ = $linesep;
    
done:
    return ($status, $head, $data);
}


sub makeMultiPartPost {
    my $post = shift;
    my ($ret_post, $boundary);
    my ($pair, $k, $v);
    my $rand;

    srand();
    $boundary = '-'x29;
    for (my $i=0; $i < 28; ++$i) {
	$rand = rand 10;
	$rand = int ($rand);
	$boundary .= $rand;
    }

    my @args = split (/\&/, $post);
    foreach $pair (@args) {
	next if ($pair =~ /^\s*$/);
	($k, $v) = split (/=/, $pair);
	$k = &unEscapeGetArg ($k);
	$v = &unEscapeGetArg ($v);
	$ret_post .= qq{$boundary\r\n};
	$ret_post .= qq{Content-Disposition: form-data; name="$k"\r\n\r\n$v\r\n};
    }
    $ret_post .= "$boundary--";

    return ($ret_post, $boundary);
}


sub getSession {
    my $uri = shift;
    $uri =~ /\@SK\@([^\@]+)\@\@/;
    return $1;
}


sub removeSession {
    my $uri = shift;

    $uri =~ s/\@SK\@[^\@]+\@\@//;
    return $uri;
}


sub insertSession {
    my ($uri, $sid) = @_;
    my $new_uri = '';
    my $sid_path = "\@SK\@$sid\@\@";
    my $uri_len = length ($uri);
    my $part1 = '';
    my $part1_len = 0;
    my $marker = 0;
    my ($octothorpe_marker, $questionmark_marker) = (-1,-1);


    if (($octothorpe_marker = index ($uri, "\#")) != -1) {
	$marker = $octothorpe_marker;
    }
    if (($questionmark_marker = index ($uri, "?")) != -1
	&& $questionmark_marker < $marker) {
	$marker = $questionmark_marker;
    }
    if ($marker > 0) {
	--$marker;
    }
    else {
	$marker = $uri_len - 1;
    }


    if (substr ($uri, $marker, 1) eq '/') {
	$new_uri = substr ($uri, 0, $marker + 1);
	$new_uri .= $sid_path;
	++$marker;
	$new_uri .= substr ($uri, $marker)  if ($marker < $uri_len);
    }
    else {
	$part1 = substr ($uri, 0, $marker + 1);
	$part1 =~ s/[^\.\/]*$//;

	if ($part1 =~ /\.$/) {
	    chop $part1;
	}
	# else ($part1 =~ /\/$/)
	
	$new_uri = $part1;
	$new_uri .= $sid_path;
	$part1_len = length ($part1);
	$new_uri .= substr ($uri, $part1_len)  if ($part1_len < $uri_len);
    }

    return $new_uri;
}


# makeGetArgs (\%a)
#
#
#   desc:  make a string suitable for GET args out of associative array
#
#   args:  \%a    array with key->val pairs
#
#   rets:  $ret   string for GET args
#
#
sub makeGetArgs {
    local ($a) = @_;
    my ($k, $s, $ret);
    foreach $k (keys %$a) {
        $s = $a->{$k};
        $s = &escapeGetArg ($s);
        $ret .= "$k=$s&";
    }
    chop $ret;                                      # Remove trailing ampersand
    return $ret;
}
# end makeGetArgs ()


# escapeGetArg ()
#
sub escapeGetArg {
    local $str = shift;
    $str =~ s/ /+/go;
    $str =~ s/\%/\%25/go;
    $str =~ s/([^0-9a-zA-Z\%+])/"\%".&charToHex($1)/ge;
    return $str;
}
# end escapeGetArg ()


# unEscapeGetArg ()
#
sub unEscapeGetArg {
    local $str = shift;
    $str =~ s/\%25/PRESERVE_PERCENT_SYMBOL/go;
    $str =~ s/\%(\d[0-9a-fA-F])/&hexToChar($1)/ge;
    $str =~ s/PRESERVE_PERCENT_SYMBOL/\%/go;
    $str =~ s/\+/ /go;
    return $str;
}
# end escapeGetArg ()


# charToHex ()
#
sub charToHex {
    my $ascii = ord($_[0]);
    my %hexMap = (  0 => '0',
                    1 => '1',
                    2 => '2',
                    3 => '3',
                    4 => '4',
                    5 => '5',
                    6 => '6',
                    7 => '7',
                    8 => '8',
                    9 => '9',
                   10 => 'a',
                   11 => 'b',
                   12 => 'c',
                   13 => 'd',
                   14 => 'e',
                   15 => 'f'
                 );

    return $hexMap{(($ascii & 0xf0) >> 4)} . $hexMap{($ascii & 0x0f)};
}
# end charToHex ()


#  hexToChar ()
#
sub hexToChar {
    my $ascii = hex($_[0]);
    return chr $ascii;
}
# end hexToChar ()


1;

###############################################################################
# END
###############################################################################
