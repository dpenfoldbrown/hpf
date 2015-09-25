# This is a PERL5 module for useful commands in SAM-T99 perl scripts.
# It also locates and sources the initialization file.
#
# $Id: SamT99.pm,v 1.5 2001/08/23 20:53:03 rph Exp $
#


package SamT99;
use English;
use Exporter;
use File::Basename;
@ISA = qw(Exporter);
@EXPORT = qw(fatal run_prog absolute absolute_or_in_dir make_directories);


# Load the configuration
if (defined($ENV{SAM_T99_CONF})) {
    $sam_t99_conf = $ENV{SAM_T99_CONF};
} else {
    $sam_t99_conf = dirname($PROGRAM_NAME) . "/sam-t99.conf";
}
print STDERR "Reading $sam_t99_conf";
do $sam_t99_conf || die "Can't read config file: $sam_t99_conf";
$path = $ENV{PATH};
print STDERR "\nPath: $path";



#
# fatal(msg)
#
# Generate a fatal error message.
#
sub fatal($) 
{   my($msg) = @_;
    die("*** Error: $PROGRAM_NAME: $msg\n");
    exit(1);
}

#
# run_prog(cmd, [errmsg])
#
# Run a command using system.  Exit if the command fails.  An optional
# error message for failure maybe specified or one will be generated.
#
sub run_prog($$) 
{   my($cmd,$errmsg) = @_;
    if ($::verbose)
    {   print STDERR "\n\n@@@@ $cmd\n";
    }

    if (system($cmd) != 0) 
    {	if (!defined($errmsg)) 
        {   $errmsg = "command failed: $cmd";
        }
        fatal($errmsg);
    }
}

# prepend relative path names with current working directory
# usage $absolute_file = absolute($file)

sub absolute($)
{    my ($file) = @_;
    return $file if ($file =~ /^\/.*/);
    my $wd = `pwd`;
    chop($wd);
    if ($wd=~/\/auto(.*)/) {return "$1/$file"};
    return "$wd/$file";

}

# prepend relative path names with current working directory
# or specified directory (must already exist in one or the other)
sub absolute_or_in_dir($$)
{
    my ($file, $dir) = @_;
    my($abs_file) = absolute($file);
    return $abs_file if (-e $abs_file);
    my($dir_file) = "$dir/$file";
    return $dir_file if (-e $dir_file);
    fatal("Can't find either $abs_file or $dir_file\n");
}


# make all directories along path to specified directory
# THIS SHOULD PROBABLY BE REPLACED 
#	by standard mkpath() from File::Path,
# 	but we found the group setting to be useful on our system.
sub make_directories ($)
{   my ($dir) = @_;
    my $groups = `groups`;
    my $in_protein_group = ($groups =~ /protein/ );

    @pieces = split /\//, $dir;
    for (my $i = 0; $i <=$#pieces; $i++)
    {	my $path = join ('/',  @pieces[0 .. $i]);
	next if $path eq "";
	if (! -d $path)
	{   mkdir $path, 0775;
	    chmod 0775, $path;
	    if ($in_protein_group )
	    {	run_prog("chgrp protein $path", "chgrp failed");
	    }
	}
    }
}


# packages must end with a true value---a restriction that seems
# to be poorly documented.
1;
