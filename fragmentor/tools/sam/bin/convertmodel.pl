#!/usr/local/bin/perl5
#
# convertmodel.pl
# written by Rachel Karchin for UCSC Computational Biology - November 1996
#
# convertmodel.pl converts the SAM model files in directories designated
# on the command line from ascii to binary format or from binary to
# ascii format.  The type of conversion is controlled by selection of
# ascii or binary.  ascii converts a model file from binary
# to ascii and binary converts a model file from ascii to binary.
#
# convertmodel.pl
# 1.   searches for all model files created by SAM 
#      in directories designated by user on the command line.
#      The user has the choice of traversing all subdirectories
#      of the chosen directories by entering the -r option or
#      doing conversions only in the selected directories (omit -r)
# 
# 2. converts each modelfile found from binary to ascii format or
#    vice versa using hmmconvert.
#    The option ascii implements a conversion from binary
#    to ascii, while binary implements a conversion from
#    ascii to binary.
#
# Usage: convertmodel.pl binary -r directory1 directory2
#
#--------------------------------------------------------------------------

if ($#ARGV < 1 ) {&print_usage_exit};

#get requested directories and output format from the command line
#user can enter arguments in any order 
foreach $i ( 0 .. $#ARGV ) {
   $next_arg = @ARGV[$i];
   if ($next_arg eq "binary" ) {
       $output_style = "binary";
   }
   if ($next_arg eq "ascii" ) {
       $output_style = "ascii";
   }
   if ($next_arg eq "-r") { 
       $recurse = "-r";
   }
   if ( -d $next_arg == 1 ) {
       push(@the_dirs, $next_arg); 
   }
   $i++;
}
  if ( ($output_style eq "") || ( @the_dirs[0] eq "" ) ) {
     &print_usage_exit;
  }
     
  $convert = "hmmconvert";
# get the model files
  while (@the_dirs) {
      $working_dir = shift(@the_dirs);
      chdir $working_dir;
      if ( $recurse eq "-r" ) {
	  &recursive_convert;
    } else {
	&non_recursive_convert;
    }
  }
  exit (1);


#----------------------- subroutines ----------------------------

sub print_usage_exit {
print "\n\n";
print "  Usage: convertmodel.pl binary <dir1> [<dir2>...]\n";
print "              Converts models in the directories to binary format\n\n";
print "         convertmodel.pl ascii <dir1> [<dir2>...]\n";
print "              Converts models in the directories to ascii format\n\n";
print "         convertmodel.pl binary -r <dir1> [<dir2>...]\n";
print "              Converts models in directories and recursively in\n";
print "              their subdirectories to binary format\n\n";
print "         convertmodel.pl ascii -r <dir1> [<dir2>...]\n";
print "              Converts models in directories and recursively in\n";
print "              their subdirectories to ascii format\n\n"; 
exit (-1);
}

# file_check checks that the file is not a link or directory, and that
# it has a .mod extension.  It returns -1 if the file is undesirable,
# 1 if it is desirable.
sub file_check {
    if ( -d $next_file ) {
        return (-1);
    }
     if ( -l $next_file ) {
        return (-1);
    }

    if ( !( $next_file =~ /mod$/ )) {
        return (-1);
    }
    return 1;
}

sub recursive_convert {
    open(FIND, "find . -name \\*.mod \\( ! -type l \\)  \\( ! -type d \\) -print |") || die "Couldn't run find: $!\n";

     while ($filename = <FIND> ) {
       chop $filename;
       &do_hmm_convert;
     }
   close(FIND);
}

sub non_recursive_convert {
     open(THEFILES, "ls | ");
     $i = 1;
     while(<THEFILES>) {
	 chop;
         print "Checking :$_:\n";
         $next_file = $_;
         $n = &file_check;
         if ( $n == 1 ) {
	    print "Adding :$_:";
	    @the_files[$i] = $next_file;
	 }
         $i++;
     }
     close(THEFILES);

     while (@the_files){
        $filename = shift(@the_files); 
        if ($filename ne "") {
            &do_hmm_convert;             
	}
    }
 }

sub do_hmm_convert {

    if (  !($filename =~ /hmmer/) && !($filename =~ /genmod/) )
    {
         $l = length($filename);
         $fname = substr($filename, 0, $l-4);
         if ($output_style eq "ascii" ) {
             system "$convert $fname -modelfile $filename -binary_output 0";
	 }
         if ($output_style eq "binary" ) {
             system "$convert $fname -modelfile $filename -binary_output 1";
	 }
     }
}


