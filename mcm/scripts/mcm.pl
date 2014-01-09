#!/usr/bin/perl -w
#### 
#
# Needs gnu scientific library installed (gsl).
#
####


# libraries - should all in installed on most systems

use strict;
use Carp;
use Getopt::Long;
use File::Basename;

# read command line
my $ar = {};
$ar->{call} = join " ", @ARGV;
my @opt = qw( test=s rundir=s decoydir=s epath=s dpath=s outfile=s ssfile=s debug=i zscore_threshold=s number_of_clustercenters=i mammothDb=s -help -force);
$ar->{debug} = 0 unless defined $ar->{debug};
&GetOptions( $ar, @opt );
&help if $ar->{help} || !-f $ar->{outfile};
sub help {
	print "Takes an rosetta silentmode file (-outfile <file>; required ) and the corresponding secondary structure prediction file (-ssfile <ssfile>), clusters, extracts decoys (number_of_clustercenters; defaults to 25), runs mammoth and estimates the confidence in each mammoth match above a zscore of zscore_threshold (defaults to 4.5);\n";
	exit;
}

# Setup
$ar->{epath} = $ENV{MCM} if $ENV{MCM} && !$ar->{epath};
if (!$ar->{epath} && $0 =~ /^\//) {
	my @parts = split /\//, $0;
	pop @parts;
	$ar->{epath} = join "\/", @parts;
}
confess "Could not detemine where the executables are located. Set either the MCM environment or give absolute path with the -epath option\n" unless $ar->{epath} && -d $ar->{epath};
$ar->{clusterer} = sprintf "%s/robetta_cluster_mod",$ar->{epath}; # clustering software. Needs to output the simulation convergence term
#kdrew: old line with Bill's extractor $ar->{extractor} = sprintf "%s/reconstruct_ROSETTA_pdb_by_index",$ar->{epath}; # Extract decoys from a silentmode file.
$ar->{extractor} = sprintf "%s/rosetta",$ar->{epath}; # Extract decoys from a silentmode file.
$ar->{mammoth} = sprintf "%s/mammoth",$ar->{epath}; # structure-structure alignment program. Needs to calculate the contact order of the matched regions
$ar->{mdbdatafile} = sprintf "%s/data.mammothDb",$ar->{dpath}; # contains data about the domains in the mammoth database
$ar->{mammothDb} = sprintf "%s/mammothDbCEMS",$ar->{dpath} unless $ar->{mammothDb}; # localtion of the mammoth database
$ar->{mammothDbList} = sprintf "%s/list.mammoth",$ar->{dpath}; # mammoth list of domains in the mammoth database

# Defaults
$ar->{zscore_threshold} = 4.5 unless $ar->{zscore_threshold};
$ar->{number_of_clustercenters} = 25 unless $ar->{number_of_clustercenters};

# for testing of individual subroutines
&test() if $ar->{test};

# Main
&init();
&check('start');
&run();
&check('end');
&post();

sub post {
	# output warnings
	unless ($#{ $ar->{warning} } < 0) {
		$ar->{xml} .= "\t<warnings>\n";
		for my $warning (@{ $ar->{warning} }) {
			$ar->{xml} .= sprintf "\t\t<warning>%s</warning>\n", $warning;
		}
		$ar->{xml} .= "\t</warnings>\n";
	}
	# terminate the xml and output to file
	$ar->{xml} .= "</mcm>\n";
	my $logfile;
	my $count = 0;
	while (1==1) {
		$count++;
		$logfile = sprintf "log.%02d.xml", $count;
		last unless -f $logfile;
	}
	confess "File exists\n" if -f $logfile;
	unlink "log.xml" if -l "log.xml";
	open XML, ">$logfile";
	print XML $ar->{xml};
	close XML;
	symlink $logfile, "log.xml";
}
sub run {
	# run the program
	my $ret;
	eval {
		$ret = &cluster(); # cluster the decoys
		print $ret if $ar->{debug} > 0;
	};
	&error( $@ ) if $@;
	eval {
		$ret = &read_convergence(); # read in the convergence
		print $ret if $ar->{debug} > 0;
	};
	&error( $@ ) if $@;
	eval {
		$ret = &extract_rosetta(); # extract cluster centers
		print $ret if $ar->{debug} > 0;
	};
	&error( $@ ) if $@;
	eval {
		$ret = &mammoth(); # mammoth cluster centers 
		print $ret if $ar->{debug} > 0;
	};
	&error( $@ ) if $@;
	eval {
		$ret = &confidence(); # calculate confidence and output to xml file
		print $ret if $ar->{debug} > 0;
	};
	&error( $@ ) if $@;
	eval {
		$ret = &decoys2xml(); # put the decoys in the xml file
		print $ret if $ar->{debug} > 0;
	};
	&error( $@ ) if $@;
}
sub error {
	# terminates the program with error message
	my($error)=@_;
	$ar->{xml} .= sprintf "<error>%s</error>\n", $error;
	warn "$error\n";
	&post();
	exit 1;
}
# testing subroutines
sub test {
	# test indivitual subroutines
	&init();
	&check('start');
	if ($ar->{test} eq 'extract') {
		print &extract();
	} elsif ($ar->{test} eq 'mammoth') {
		print &mammoth();
	} elsif ($ar->{test} eq 'cluster') {
		print &cluster();
	} elsif ($ar->{test} eq 'confidence') {
		print &confidence();
	} elsif ($ar->{test} eq 'read_mdbdata') {
		print &read_mdbdata();
	} else {
		confess "Unknown test: $ar->{test}\n";
	}
	exit 0;
}
sub init {
	# setup and parse essential information. Initialize xml
	&setup('pwd');
	$ar->{outfile} = "$ar->{pwd}/$ar->{outfile}" unless $ar->{outfile} =~ /^\//;
	&setup('mammothDbList');
	&setup('rundir') unless $ar->{rundir};
	&setup('sequence');
	&setup('ss');
	&setup('hostname');
	&check('logfile');
	$ar->{xml} = "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";
	$ar->{xml} .= "<mcm>\n";
	$ar->{xml} .= "<param>\n"; # all parameters needed to reproduce the run
	my $now_time = localtime;
	$ar->{xml} .= sprintf "\t<runtime>%s</runtime>\n",$now_time;
	$ar->{xml} .= sprintf "\t<hostname>%s</hostname>\n",$ar->{hostname};
	$ar->{xml} .= sprintf "\t<number_of_clustercenters>%d</number_of_clustercenters>\n",$ar->{number_of_clustercenters};
	$ar->{xml} .= sprintf "\t<zscore_threshold>%s</zscore_threshold>\n",$ar->{zscore_threshold};
	$ar->{xml} .= sprintf "\t<call>%s</call>\n",join " ", $ar->{call};
	$ar->{xml} .= sprintf "\t<pwd>%s</pwd>\n",join " ", $ar->{pwd};
	$ar->{xml} .= sprintf "\t<rundir>%s</rundir>\n",join " ", $ar->{rundir};
	$ar->{xml} .= sprintf "\t<decoydir>%s</decoydir>\n",join " ", $ar->{decoydir};
	$ar->{xml} .= sprintf "\t<ss>%s</ss>\n",join " ", $ar->{ss};
	$ar->{xml} .= sprintf "\t<sequence>%s</sequence>\n",join " ", $ar->{sequence};
	$ar->{xml} .= sprintf "\t<outfile>%s</outfile>\n",join " ", $ar->{outfile};
	$ar->{xml} .= sprintf "\t<ssfile>%s</ssfile>\n",join " ", $ar->{ssfile};

	$ar->{xml} .= sprintf "\t<epath>%s</epath>\n",join " ", $ar->{epath};
	$ar->{xml} .= sprintf "\t<clusterer>%s</clusterer>\n",join " ", $ar->{clusterer};
	$ar->{xml} .= sprintf "\t<extractor>%s</extractor>\n",join " ", $ar->{extractor};
	$ar->{xml} .= sprintf "\t<mammoth>%s</mammoth>\n",join " ", $ar->{mammoth};
	$ar->{xml} .= sprintf "\t<mdbdatafile>%s</mdbdatafile>\n",join " ", $ar->{mdbdatafile};
	$ar->{xml} .= sprintf "\t<mammothDb>%s</mammothDb>\n",join " ", $ar->{mammothDb};
	$ar->{xml} .= sprintf "\t<mammothDbList>%s</mammothDbList>\n",join " ", $ar->{mammothDbList};
	$ar->{xml} .= "</param>\n";
	print "init ok\n" if $ar->{debug} > 1;
}
sub setup {
	# will setup various aspects, such as rundir and hostname
	my $mode = shift;
	if ($mode eq 'mammothDbList') {
		return '' if -f $ar->{mammothDbList};
		confess "Cannot find mammothDb directory\n" unless -d $ar->{mammothDb};
		chdir $ar->{mammothDb};
		open LIST, ">$ar->{mammothDbList}";
		print LIST "MAMMOTH List\n";
		printf LIST "%s\n", $ar->{mammothDb};
		my @pdbs = glob("*.pdb");
		chomp @pdbs;
		for (@pdbs) {
			printf LIST "%s\n", $_;
		}
		confess "Unsuccessful in producing mammothDbList\n" unless -f $ar->{mammothDbList};
	} elsif ($mode eq 'rundir') {
		# rundir (the work-directory). All outfiles are written out in this directory
		# Get from outfile
		confess "No outfile\n" unless $ar->{outfile};
		confess "Cannot find outfile\n" unless -f $ar->{outfile};
		$ar->{rundir} = $ar->{outfile}.".p";
		mkdir $ar->{rundir} unless -d $ar->{rundir};
		confess "Could not make rundir $ar->{rundir}: $!\n" unless -d $ar->{rundir};
		chdir $ar->{rundir};
		print "setup-rundir ok\n" if $ar->{debug} > 1;
	} elsif ($mode eq 'hostname') {
		# hostname of the machine. Useful when running on clusters
		$ar->{hostname} = $ENV{hostname};
		$ar->{hostname} = `hostname` unless $ar->{hostname};
		chomp $ar->{hostname};
	} elsif ($mode eq 'pwd') {
		# the directory where the silentmode-file and the ss-file are located. the rundir will be put in this directory
		$ar->{pwd} = `pwd`;
		chomp $ar->{pwd};
		confess "Cannot find work-directory\n" unless -d $ar->{pwd};
		print "setup-pwd ok\n" if $ar->{debug} > 1;
	} elsif ($mode eq 'sequence') {
		# get the sequence from the silentmode file
		$ar->{sequence} = `head -1 $ar->{outfile}`;
		chomp $ar->{sequence};
		$ar->{sequence} =~ s/^SEQUENCE: // || confess "Could not remove expected tag from sequence\n";
		confess "Could not parse out sequence\n" unless length($ar->{sequence}) > 0;
		print "setup-sequence ok\n" if $ar->{debug} > 1;
	} elsif ($mode eq 'ss') {
		# parse the ss-prediction
		unless ($ar->{ssfile}) {
			$ar->{ssfile} = $ar->{outfile};
			if ($ar->{outfile} =~ /out$/) {
				$ar->{ssfile} =~ s/out$/ss/ || confess "Could not replace expected extension\n";
			} else {
				$ar->{ssfile} .= '.ss';
			}
		}
		$ar->{ssfile} = "$ar->{pwd}/$ar->{ssfile}" unless $ar->{ssfile} =~ /^\//;
		confess "ss file same as outfile\n" if $ar->{ssfile} eq $ar->{outfile};
		confess "Cannot find ss file ($ar->{ssfile})\n" unless -f $ar->{ssfile};
		&parseSS();
		print "setup-ss ok\n" if $ar->{debug} > 1;
	} else {
		confess "Unknown $mode\n";
	}
}
sub check {
	my $mode = shift;
	# try to detemine if everything needed to run is present
	if ($mode eq 'logfile') {
		my @logfiles = glob("$ar->{rundir}/*.xml");
		if ($#logfiles >= 0 && !$ar->{force}) {
			warn "Has xml-log. Will not run unless forced\n";
			exit;
		}
	} elsif ($mode eq 'start') {
		# setup rundir
		# look for executables, directories, outfiles etc
		confess "Cannot find mammoth executable\n" unless -x $ar->{mammoth};
		confess "Cannot find extractor executable\n" unless -x $ar->{extractor};
		confess "Cannot find clusterer executable\n" unless -x $ar->{clusterer};
		confess "Cannot find mammothDb directory\n" unless -d $ar->{mammothDb};
		confess "Cannot find mammothDb data file ($ar->{mdbdatafile})\n" unless -f $ar->{mdbdatafile};
		confess "Cannot find mammoth list $ar->{mammothDbList} \n" unless -f $ar->{mammothDbList};
		confess "Cannot find file $ar->{outfile} \n" unless -f $ar->{outfile};
		printf "check-init ok\n" if $ar->{debug} > 0;
	} elsif ($mode eq 'end') {
		# make sure all expected files were produced
		confess "Mammoth cmd file missing\n" unless -f "cmd.mammoth";
		confess "Cluster cmd file missing\n" unless -f "cmd.cluster";
		confess "Mammoth data file missing\n" unless -f "data.mammoth";
		confess "cluster data file missing\n" unless -f "data.cluster";
		confess "cluster log file missing\n" unless -f "log.cluster";
		confess "mammoth log file missing\n" unless -f "log.mammoth";
		confess "extract log file missing\n" unless -f "log.extract";
		my @decoys = glob("decoy*.pdb");
		confess "no cluster centers\n" if $#decoys < 0;
		printf "check-end ok\n" if $ar->{debug} > 0;
	} else {
		confess "Unknown $mode\n";
	}
}
sub parseSS {
	# parse the ss prediction
	confess "Needs sequence in parseSS to compare length\n" unless $ar->{sequence};
	open SS, "<$ar->{ssfile}";
	my @content = map{ $_ =~ s/Pred: //; $_ }grep{ /^Pred: / }<SS>;
	close SS;
	chomp @content;
	map{$_ =~ s/Pred: //;}@content ;
	$ar->{ss} = join "", @content;
	my $ss_len = length($ar->{ss});
	my $sequence_len = length($ar->{sequence});

	print "sequence: $ar->{sequence}, $sequence_len\n";
	print "ss: $ar->{ss}, $ss_len\n";
	confess "Lengths differ\n" unless length($ar->{sequence}) == length($ar->{ss});
	my $tmpss = $ar->{ss};
	my $ppa = $tmpss =~ s/H//g;
	my $ppb = $tmpss =~ s/E//g;
	my $ppc = $tmpss =~ s/C//g;
	confess "Still some left $tmpss\n" if $tmpss;
	$ar->{prediction_percent_alpha} = $ppa / length($ar->{sequence});
	$ar->{prediction_percent_beta} = $ppb / length($ar->{sequence});
	$ar->{prediction_percent_coil} = $ppc / length($ar->{sequence});
}
### RUN SUB ###
sub cluster {
	# Generate cmd file
	return "cmd.cluster file exists - assuming clustering has be done. Not overwriting.\n" if -f 'cmd.cluster' && !$ar->{overwrite};
	$ar->{xml} .= sprintf "<clusterer>\n";
	confess "No sequence in cluster\n" unless $ar->{sequence};
	open CLUSTER, ">cmd.cluster";
	printf CLUSTER sprintf "OUTPUT_FILE data.cluster\n";
	printf CLUSTER sprintf "TARGET %s %s\n",$ar->{outfile},$ar->{sequence};
	close CLUSTER;
	my $shell = "$ar->{clusterer} cmd.cluster >& log.cluster";
	$ar->{xml} .= sprintf "<shell>%s</shell>\n", $shell;
	`$shell`;
	$ar->{xml} .= sprintf "</clusterer>\n";
	# Example cmd.cluster files
	# Example 1 - with homologs:
	#OUTPUT_FILE ./aat000_.clusters
	#TARGET aat000.out MATRTPKLVKHTLLTRFKDEITREQIDNYINDYTNLLDLIPSMKSFNWGTDLGMESAELNRGYTHAFESTFESKSGLQEYLDSAALAAFAEGFLPTLSQRLVIDYFLY
	#HOMOLOG aah001.out --------VEHVIFHCWSSDVSD--VAAMQAEGVATLGALPGVRRIAAGVAETP-----NAPYPYLWAMTFATPEALAHFRTHPEHQRFASRFRPFASERITIDFHL-
	#HOMOLOG aah002.out -------MLFHQVFFWLKNPGDKADRDKLIAGLKALKA-IDVIQQLHVGVPAATERDVVDNSYDVSELMVFKSVEDQKRYRDHPLLQKFVADCSHLWSKVVVYD----
	# Example 2 - without homologs:
	#OUTPUT_FILE ./aat000_.clusters_mod
	#TARGET test/aat000_.out MATRTPKLVKHTLLTRFKDEITREQIDNYINDYTNLLDLIPSMKSFNWGTDLGMESAELNRGYTHAFESTFESKSGLQEYLDSAALAAFAEGFLPTLSQRLVIDYFLY
}
sub mammoth {
	my(%param)=@_;
	# Construct prediction list
 	if (-f 'data.mammoth' && !$ar->{overwrite}) {
		my $tail = `tail -1 data.mammoth`; # look at last line
		chomp $tail;
		if ($tail =~ /NORMAL_EXIT/) { # see if mammoth completed as expected
			return "Have data.mammoth file. Assuming mammoth is complete. Not overwriting.\n";
		} else {
			return "Have data.mammoth file, but cannot find completion tag in last line (last line: $tail). Not overwriting.\n";
		}
	}
	$ar->{xml} .= "<mammoth>\n";

	#kdrew: glob seems to get an invalid directory, trying without glob
	#my @files = glob("$ar->{rundir}/decoy_*.pdb");
	#my @files = `ls $ar->{decoydir}/decoy_*.pdb`;
	my @files = `ls $ar->{rundir}/decoy_*.pdb`;
	confess "Could not find any files\n" if $#files < 0;
	$ar->{xml} .= sprintf "<n_decoys>%s</n_decoys>\n",$#files+1;
	open LIST, ">cmd.mammoth";
	print LIST "MAMMOTH List\n";
	#printf LIST "/\n"; # ,$directory;
	printf LIST "\n"; # ,$directory;
	for my $file (@files) {
		# patrick 09-26-2008
		# Printing basename instead of file.  Problems with mammoth paths.  Mammoth will load the decoys relative to the cmd.mammoth file.
		my $basename = basename($file);
		printf LIST "%s\n",$basename;
	}
	close LIST;
	confess "Could not create file...\n" unless -f "cmd.mammoth";
	# Mammoth
	my $shell = "$ar->{mammoth} -p cmd.mammoth -e $ar->{mammothDbList} -o data.mammoth -v 0 -r 0 -c 1 >& log.mammoth";
	$ar->{xml} .= sprintf "<shell>%s</shell>\n",$shell;
	`$shell`;
	confess "No data.mammoth file produced.\n" unless -f 'data.mammoth';
	my $tail = `tail -1 data.mammoth`; # look at last line
	chomp $tail;
	$ar->{xml} .= sprintf "<tail>%s</tail>\n",$tail;
	$ar->{xml} .= "</mammoth>\n";
	confess "Mammoth did not complete as expected. Last line: $tail.\n" unless $tail =~ /NORMAL_EXIT/;
}

#kdrew: new extract subroutine to use rosetta instead of Bill's extractor
sub extract_rosetta {
	my(%param)=@_;
	my @decoys = glob("decoy*.pdb");

	#kdrew: temporary file used as input into rosetta's extraction function
	my $tmp_filename = "tmpList";

	#return "Have decoys. Not overwriting\n" unless $#decoys < 0;
	$ar->{xml} .= "<extract>\n";
	confess "No infile\n" unless -f 'data.cluster';
	open IN, "<data.cluster";
	my $count = 0;
	while (<IN>) {
		last if $_ =~ /CLUSTER MEMBER/;
		#kdrew: changing regex to fit hpf2
		#if ($_ =~ /^\d+\:\s+(\d+),(\w+)\s+(\d+)\s/) {
		#patrick: Allowing both regex's to parse.  Whichever works!
		if ($_ =~ /^\d+\:\s+(\d+),.+(S_.+)\s+(\d+)\s/  || $_ =~ /^\d+\:\s+(\d+),(\w+)\s+(\d+)\s/ ) {
			$ar->{xml} .= sprintf "\t<center%03d>%s</center%03d>\n",$count,$1,$count;
			$ar->{xml} .= sprintf "\t<center%03dsize>%s</center%03dsize>\n",$count,$3,$count;
			$ar->{center}->{$1}->{rank} = $count;
			$ar->{center}->{$1}->{size} = $3;
            
            # From the RE captures (against a cluster center line from a robetta_cluster_mod output file):
            # 1 - index (#) of the silent file record in the silent file given to cluster
            # 2 - the S_1234_1234 decoy ID rosetta gives in the silent file (labelled as "description")
            #     NOTE: in clusterer output, this ID has an underscore appended. Removed below.
            # 3 - the size of the cluster as reported in clustered outfile

			my $tmp1 = $1;
			my $tmp2 = $2;
			my $tmp3 = $3;

			#kdrew: strips ending underscore if present	
			if($tmp2 =~ m/_$/)
			{
				$tmp2 = $`;
			}
			#kdrew: put the identifier of the decoy into a file so rosetta can extract it
				#kdrew: use append and move call to rosetta out of loop for faster processing
			`echo $tmp2 > $tmp_filename`;
			#kdrew: cp paths.txt to current directory for rosetta to work
				#kdrew: bad form, move out of the loop and parameterize location of path.txt
			`cp $ar->{epath}/paths.txt .`;
			my $basename = basename($ar->{outfile});
			`$ar->{extractor} -score -silent_input -s $basename -refold -nstruct 1 -new_reader -l $tmp_filename >> log.extract` if $#decoys < 0;
			#`$ar->{extractor} -score -silent_input -s $ar->{outfile} -refold -nstruct 1 -new_reader -l $tmp_filename >> log.extract` if $#decoys < 0;
			`mv $tmp2*.pdb decoy_$tmp1.pdb`;
			last if ++$count >= $ar->{number_of_clustercenters};
		} else {
			push @{ $ar->{warning} }, "extract: Cannot parse $_";
			next;
		}
	}
	close IN;
	$ar->{xml} .= sprintf "\t<n_decoys>%s</n_decoys>\n",$count;
	$ar->{xml} .= "</extract>\n";
}
#kdrew: old subroutine to use Bill's extractor
sub extract {
	my(%param)=@_;
	my @decoys = glob("decoy*.pdb");
	#return "Have decoys. Not overwriting\n" unless $#decoys < 0;
	$ar->{xml} .= "<extract>\n";
	confess "No infile\n" unless -f 'data.cluster';
	open IN, "<data.cluster";
	my $count = 0;
	while (<IN>) {
		last if $_ =~ /CLUSTER MEMBER/;
		if ($_ =~ /^\d+\:\s+(\d+),\w+\s+(\d+)\s/) {
			$ar->{xml} .= sprintf "\t<center%03d>%s</center%03d>\n",$count,$1,$count;
			$ar->{xml} .= sprintf "\t<center%03dsize>%s</center%03dsize>\n",$count,$2,$count;
			$ar->{center}->{$1}->{rank} = $count;
			$ar->{center}->{$1}->{size} = $2;
			`$ar->{extractor} $ar->{outfile} $1 >> log.extract` if $#decoys < 0;
			last if ++$count >= $ar->{number_of_clustercenters};
		} else {
			push @{ $ar->{warning} }, "extract: Cannot parse $_";
			next;
		}
	}
	close IN;
	$ar->{xml} .= sprintf "\t<n_decoys>%s</n_decoys>\n",$count;
	$ar->{xml} .= "</extract>\n";
}
sub read_mdbdata {
	open MDBDATA, "<$ar->{mdbdatafile}";
	my %data;
	while (<MDBDATA>) {
		chomp;
		my @cols = split /\t/, $_;
		confess sprintf "Wrong number of columns read (%d; %s)\n", $#cols, $_ unless $#cols == 6;
		next if $cols[0] eq 'structure_key';
		$data{$cols[0]} = { sequence_key => $cols[1], length => $cols[2], percent_alpha => $cols[3], percent_beta => $cols[4], sccs => $cols[5],astral_ac => $cols[6] };
	}
	if ($ar->{debug} > 2) {
		for my $structure_key (keys %data) {
			for (keys %{ $data{$structure_key} }) {
				printf "%s => %s\n", $_, $data{$structure_key}->{$_};
			}
		}
	}
	close MDBDATA;
	return \%data;
}
sub read_convergence {
	my $string = `head -2 data.cluster | tail -1`;
	my ($size1,$radius1,$size2,$radius2,$total_n) = $string =~ /standard_thresholds: size1= (\d+) threshold1= ([\.\d]+) size2= (\d+) threshold2= ([\.\d]+) total_decoys= (\d+)/;
	$ar->{xml} .= "<convergence>\n";
	$ar->{xml} .= sprintf "\t<size1>%s</size1>\n", $size1;
	$ar->{xml} .= sprintf "\t<radius1>%s</radius1>\n", $radius1;
	$ar->{xml} .= sprintf "\t<size2>%s</size2>\n", $size2;
	$ar->{xml} .= sprintf "\t<radius2>%s</radius2>\n", $radius2;
	$ar->{xml} .= sprintf "\t<total_n_decoys>%s</total_n_decoys>\n", $total_n;
	confess "Cannot read...\n" unless $size1 && $radius1;
	$ar->{xml} .= "</convergence>\n";
	$ar->{convergence} = $radius1;
	$ar->{n_decoys_in_outfile} = $total_n;
	return '';
}
sub decoys2xml {
	my @decoys = glob("decoy*.pdb");
	undef $/;
	$ar->{xml} .= "<decoys>\n";
	for my $decoy (@decoys) {
		open IN, "<$decoy";
		my $content = <IN>;
		confess "Nothing read form $decoy\n" unless $content;
		$ar->{xml} .= sprintf "<decoy><name>%s</name><atomrecord>\n%s\n</atomrecord></decoy>\n", $decoy,$content;
		close IN;
	}
	$ar->{xml} .= "</decoys>\n";
}
sub confidence {
	confess "Could not find the mammoth data file\n" unless -f "data.mammoth";
	confess "No sequence\n" unless $ar->{sequence};
	confess "No prediction_percent_alpha\n" unless defined($ar->{prediction_percent_alpha});
	confess "No prediction_percent_beta\n" unless defined($ar->{prediction_percent_beta});
	$ar->{xml} .= "<confidence>\n";
	open MDATA, "<data.mammoth";
	my $mdbdata = &read_mdbdata();
	confess "No convergence\n" unless $ar->{convergence};
	confess "No n_decoys_in_outfile\n" unless $ar->{n_decoys_in_outfile};
	#$ar->{xml} .= sprintf "<header>: p e Zscore -ln(E) Evalue score #e #p nsup nss psi1 psi2 CO_p CO_e prediction_file experiment_file convergence prediction_length prediction_percent_alpha prediction_percent_beta experiment_sequence_key experiment_length experiment_percent_alpha experiment_percent_beta experiment_sccs experiment_astral_ac probability</header>\n";
	my $target = (split /\//, $ar->{outfile})[-1];
	$ar->{xml} .= "<data>\n";
	while (<MDATA>) {
		next unless substr($_,0,1) eq ':';
		my %data;
		my $junk;my $rest;
		($junk,$data{prediction_index},$data{experiment_index},$data{zscore},$data{ln_e},$data{evalue},$data{score},$data{prediction_sequence_length},$data{experiment_sequence_length},$data{nsup},$data{nss},$data{psi1},$data{psi2},$data{prediction_contact_order},$data{experiment_contact_order},$data{prediction_file},$data{experiment_file},$rest) = split /\s+/, $_;
		confess "Has rest..\n" if $rest;
		# filter lines with less and threshold
		next if $data{zscore} eq '******';
		next if $data{zscore} < $ar->{zscore_threshold};
		$data{convergence} = $ar->{convergence};
		$data{n_decoys_in_outfile} = $ar->{n_decoys_in_outfile};
		$data{target} = $target;
		$data{prediction_percent_alpha} =   $ar->{prediction_percent_alpha};
		$data{prediction_percent_beta} =   $ar->{prediction_percent_beta};
		my ($center_index)  = $data{prediction_file} =~ /decoy_(\d+).pdb$/;
		confess "No center index read from $data{prediction_file}\n" unless $center_index;
		confess "Not defined.\n" unless defined $ar->{center}->{$center_index}->{rank};
		$data{cluster_center_index} = $center_index;
		$data{cluster_center_rank} = $ar->{center}->{$center_index}->{rank};
		$data{cluster_center_size} = $ar->{center}->{$center_index}->{size};
		my ($structure_key)  = $data{experiment_file} =~ /^(\d+).pdb$/;
		confess "No structure key read from $data{experiment_file}\n" unless $structure_key;
		if ($mdbdata->{$structure_key}) {
			$data{experiment_sequence_key} = $mdbdata->{$structure_key}->{'sequence_key'} || confess "No sequence_key\n";
			$data{experiment_percent_alpha} = $mdbdata->{$structure_key}->{'percent_alpha'};
			$data{experiment_percent_beta} = $mdbdata->{$structure_key}->{'percent_beta'};
			$data{experiment_sccs} = $mdbdata->{$structure_key}->{'sccs'} || confess "No sccs\n";
			$data{experiment_astral_ac} = $mdbdata->{$structure_key}->{'astral_ac'} || confess "No astral_ac\n";
		} else {
			push @{ $ar->{warning} }, "No mdbdata for $structure_key";
			next;
		}
		$data{probability} = &probability(\%data);
		$ar->{xml} .= "\t<entry>\n";
		for my $key (sort keys %data) {
			confess "Something is wrong for $key $data{$key}\n" unless defined( $key ) && defined( $data{$key} );
			$ar->{xml} .= sprintf "\t\t<%s>%s</%s>\n", $key, $data{$key},$key;
		}
		$ar->{xml} .= "\t</entry>\n";
	}
	$ar->{xml} .= "</data>\n";
	close MDATA;
	$ar->{xml} .= "</confidence>\n";
	return '';
}
sub probability {
	my($hash)=@_;
	$hash->{class} = 3; # default
	$hash->{class} = 1 if $hash->{prediction_percent_alpha} >= 0.15 && $hash->{prediction_percent_beta} < 0.15; # alpha proteins
	$hash->{class} = 2 if $hash->{prediction_percent_alpha} < 0.15 && $hash->{prediction_percent_beta} >= 0.15; # beta proteins
	$hash->{ratio} = $hash->{prediction_sequence_length}/$hash->{experiment_sequence_length};
	if ($hash->{prediction_percent_alpha} != 0 && $hash->{experiment_percent_alpha} != 0) {
		$hash->{aratio} = $hash->{prediction_percent_alpha}/$hash->{experiment_percent_alpha};
	} elsif ($hash->{prediction_percent_alpha} == 0 && $hash->{experiment_percent_alpha} == 0) {
		$hash->{aratio} = 1;
	} elsif ($hash->{prediction_percent_alpha} == 0) {
		$hash->{aratio} = 0.0067379;
	} elsif ($hash->{experiment_percent_alpha} == 0) {
		$hash->{aratio} = 148.41;
	} else {
		confess "this cannot happend...\n";
	}
	if ($hash->{prediction_percent_beta} != 0 && $hash->{experiment_percent_beta} != 0) {
		$hash->{bratio} = $hash->{prediction_percent_beta}/$hash->{experiment_percent_beta};
	} elsif ($hash->{prediction_percent_beta} == 0 && $hash->{experiment_percent_beta} == 0) {
		$hash->{bratio} = 1;
	} elsif ($hash->{prediction_percent_beta} == 0) {
		$hash->{bratio} = 0.0067379;
	} elsif ($hash->{experiment_percent_beta} == 0) {
		$hash->{bratio} = 148.41;
	} else {
		confess "this cannot happend...\n";
	}
	my %value;
	# A
	$value{1}->{resp} = -1.833620
		+ $hash->{zscore} * 0.453397
		+ $hash->{prediction_contact_order} * 0.055677
		+ $hash->{convergence} * -0.112188
		+ abs(log($hash->{ratio})) * -2.591992
		+ abs(log($hash->{aratio})) * -1.771192
		+ abs(log($hash->{bratio})) * -0.009326;
	# B
	$value{2}->{resp} = 0.794572
		+ $hash->{zscore} * 0.408039
		+ $hash->{prediction_contact_order} * 0.040630
		+ $hash->{convergence} * -0.303383
		+ abs(log($hash->{ratio})) * -4.636701
		+ abs(log($hash->{aratio})) * -0.033375
		+ abs(log($hash->{bratio})) * -0.910910;
	# CD
	$value{3}->{resp} = -0.078608
		+ $hash->{zscore} * 0.442554
		+ $hash->{prediction_contact_order} * -0.018660
		+ $hash->{convergence} * -0.157854
		+ abs(log($hash->{ratio})) * -4.027780
		+ abs(log($hash->{aratio})) * -0.802676
		+ abs(log($hash->{bratio})) * -0.768912;
	$value{2}->{prob} = 1/(1+1/exp($value{2}->{resp}));
	$value{1}->{prob} = 1/(1+1/exp($value{1}->{resp}));
	$value{3}->{prob} = 1/(1+1/exp($value{3}->{resp}));
	#printf "%s %s %s %s %s %s %s\n", $hash->{id}, $hash->{ratio},$hash->{aratio},$hash->{bratio}, $hash->{prediction_sccs},$hash->{class},$diff;
	#printf "A:  %s %s %s\n",$value{1}->{resp},$value{1}->{prob},$hash->{probability};
	#printf "B:  %s %s %s\n",$value{2}->{resp},$value{2}->{prob},$hash->{probability};
	#printf "CD: %s %s %s\n",$value{3}->{resp},$value{3}->{prob},$hash->{probability};
	return $value{$hash->{class}}->{prob};
}
