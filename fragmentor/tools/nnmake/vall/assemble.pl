#!/usr/bin/perl

# assemble.pl - originally by Dylan Chivian <dylan@lazy8.com>
# Takes two arguments: a file listing PDB id's, one per line, and
# an output file for making the vall database.

# Maintained by James Thompson <tex@u.washington.edu>

use strict;
use File::Basename;

$0 = basename $0;
my $usage = <<USAGE
$0 input_pdb_list output_filename
    where input_pdb_list is a file containing a list of pdb id's one per line,
    and output_filename specifies the output file for the Vall.
USAGE
;

if ( scalar(@ARGV) != 2 ) {
    print STDERR "Error: incorrect number of arguments!\n";
    print STDERR $usage, "\n";
    exit -1;
}

# Configuration Options
my $clean_pdb_dir     = $ENV{HOME} . '/vall/corrected/pdb';
my $idealized_pdb_dir = $ENV{HOME} . '/vall/corrected/ideal';
my $profile_directory = $ENV{HOME} . '/vall/corrected/profile';
my $dssp_dir          = $ENV{HOME} . '/vall/corrected/dssp';

# Get the list of PDB ids that we'll be working with.
open PDB_LIST, "<$ARGV[0]" or die $!;
my @pdb_list = <PDB_LIST>;
close PDB_LIST or die $!;

#@pdb_list = map { chomp $_; uc $_ } @pdb_list;
@pdb_list = map { chomp $_; $_ } @pdb_list;

# Iterate over the list of PDB files, add appropriate information to 
# the Vall.
my $vall = ''; # $vall is a \n separated representation of the Vall
foreach my $pdb_id ( @pdb_list ) {
    $vall .= make_vall_line( $pdb_id );
}

open  VALL, ">$ARGV[1]" or die $!;
print VALL $vall;
close VALL or die $!;

#### Subs #########

sub make_vall_line {
    my $pdb_id = shift;
    
    my $vall_lines = '';

    my $profile_file  = "$profile_directory/$pdb_id.ASCII";
    my $cleanpdb_file = "$clean_pdb_dir/$pdb_id.pdb";
    my $dssp_file     = "$dssp_dir/$pdb_id.dssp";
    my $ideal_file    = "$idealized_pdb_dir/$pdb_id"."_0001.pdb";

    if ( !-f $profile_file )  { print "Can't find $profile_file\n"; next };
    if ( !-f $cleanpdb_file ) { print "Can't find $cleanpdb_file\n"; next };
    if ( !-f $dssp_file )     { print "Can't find $dssp_file\n";    next };
    if ( !-f $ideal_file )    { print "Can't find $ideal_file\n";   next };

    my @profile_buf  = @{ &slurp_file($profile_file) };
    my @cleanpdb_buf = @{ &slurp_file($cleanpdb_file) };
    my @dssp_buf     = @{ &slurp_file($dssp_file) };
    my @ideal_buf    = @{ &slurp_file($ideal_file) };

    # column order: A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
    my $header_done = undef;
    my $profiles    = []; # Array of arrays containing profiles
    my @res_ids     = (); # Contains a list of the residues such that:
                          # $res_ids[$i] represents the original amino acid
                          # for the profile at $profiles[$i+1]
    my $res_i = 0;                         
    foreach my $line (@profile_buf) {
        my ($res_id, @freqs) = split (/\s+/, $line);

        pop (@freqs); # remove trailing character

        foreach my $freq (@freqs) {
            push (@{$profiles->[$res_i]}, $freq);
        }
        push (@res_ids, $res_id);

        ++$res_i;
    }

    # Get real list of residues from cleaned PDB file.
    my @real_res_i = (); # Contains a list of the integer indexes i,
                         # such that real_res_i[$i] implies that @res_ids[$i]
                         # specifies a residue with molecular weight.

    my $last_res_index;
    foreach my $line (@cleanpdb_buf) {
        # read atom records only
        if ( $line =~ /^ATOM/ )  {     
            my $res_index = substr ($line, 22, 4);
            $res_index =~ s/^\s*|\s*$//g;
            if ( !defined $last_res_index || ($res_index != $last_res_index) ) {
                $last_res_index = $res_index;
                push (@real_res_i, $res_index); # the first res in a chain is '1'
            }
        }
    }

    # Get ideal Ca positions from the idealized PDB files.
    my $ideal_info = [];
    for (my $i=0, $res_i=0; $ideal_buf[$i] =~ /^ATOM/ && $i <= $#ideal_buf; ++$i) {
        my $name = substr($ideal_buf[$i], 12, 4);

        if ($name =~ /CA/) {
            my $chain = substr($ideal_buf[$i],21,1);
            if ( $chain !~ /^\s*$/ ) { 
                $ideal_info->[$res_i]->{chain} = $chain;
            } else {
                $ideal_info->[$res_i]->{chain} = '_';
            }
            $ideal_info->[$res_i]->{CA_x}  =  substr($ideal_buf[$i],30,8);
            $ideal_info->[$res_i]->{CA_y}  =  substr($ideal_buf[$i],38,8);
            $ideal_info->[$res_i]->{CA_z}  =  substr($ideal_buf[$i],46,8);
            $ideal_info->[$res_i]->{CA_x}  =~ s/^\s+|\s+$//g;
            $ideal_info->[$res_i]->{CA_y}  =~ s/^\s+|\s+$//g;
            $ideal_info->[$res_i]->{CA_z}  =~ s/^\s+|\s+$//g;
            ++$res_i;
        }
    }

    # Get phi psi and omega values from the idealized PDB file.
    #
    my $angles_started = undef;
    
    my $begin = 0; #Begin of segment
            
    for (my $i = 0, $res_i = 0; $i <= $#ideal_buf; ++$i) {
        if ($ideal_buf[$i] =~ /^complete/) { #idealized file format was changed
            $angles_started = 'true';
            next;
        }
        next if (! $angles_started);
        next if ($ideal_buf[$i] !~ /^\s*\d+\s.\s*([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\.\-]+)/); #as well here
        my $phi = $1;
        my $psi = $2;
        my $omg = $3;

        $ideal_info->[$res_i]->{phi} = $phi;
        $ideal_info->[$res_i]->{psi} = $psi;
        $ideal_info->[$res_i]->{omg} = $omg;
        $ideal_info->[$res_i]->{chi} = '0.000';
        ++$res_i;

    }

    # get SS information from DSSP output file
    my $i = 1;
    foreach my $line (@dssp_buf) {
        # Skip the DSSP file header
        if ( $line =~ /#  RESIDUE AA STRUCTURE/ ) {
            $header_done = 'true';
            next;
        }
        next if (!$header_done);

        # Extract the 2' structure from this line, convert
        # it to a 3-letter alphabet using &ss3state.
        my $ss_dssp = substr ($line, 16, 1);
        $ideal_info->[$i]->{ss} = &ss3state($ss_dssp);

        ++$i;
    }

    # build vall row    

    for (my $i=0; $i <= $#{$ideal_info}; ++$i) {
        #$res_i = $real_res_i[$i] - 1;       # our lists start from '0', not '1'
        $res_i = $real_res_i[$i];       # our lists start from '0', not '1'

        next if (!$profiles->[$res_i]);         # not all fastas find matches!
        next if ($res_ids[$res_i] eq 'X');
        #next if ($up[$i] eq "discard");     #discard too short and idealized badly fragments.

        my $total_residues = $#{$ideal_info} + 2;

        if ( !$ideal_info->[$i]->{ss} ) {
            $ideal_info->[$i]->{ss} = 'L';
        }

        # Set up some arbitrary numbers for now for now
        my $gap = '0.000';
        my $nalign = '10';
        my $acc = '0.00';
        my $protein_position_begin = $i + 1;
        my $protein_position_end   = $total_residues - ($i + 2);
        #my $protein_position_begin = $i;
        #my $protein_position_end   = $total_residues - $i;
        my $uc_chain = $ideal_info->[$i]->{chain};

        # fix the pdb_id 
        my $tmp = lc(substr($pdb_id,0,4));
        $pdb_id = $tmp . substr($pdb_id,4,1);
        
        # the rest
        $vall_lines .= sprintf ("%5s %1s %1s %5d %4d %4d %8.2f %8.2f %8.2f %8.3f %8.3f %8.3f %8.3f %3d %4.2f %5.3f",
                                    $pdb_id,
                                    $res_ids[$res_i],
                                    $ideal_info->[$i]->{ss},
                                    $real_res_i[$i],
                                    $protein_position_end,
                                    $protein_position_begin,
                                    $ideal_info->[$i]->{CA_x},
                                    $ideal_info->[$i]->{CA_y},
                                    $ideal_info->[$i]->{CA_z},
                                    $ideal_info->[$i]->{phi},
                                    $ideal_info->[$i]->{psi},
                                    $ideal_info->[$i]->{omg},
                                    $ideal_info->[$i]->{chi},
                                    $nalign,
                                    $acc,
                                    $gap
                                   );
        # profiles
        for (my $j=0; $j < 20; ++$j) {
            $vall_lines .= sprintf (" %5.3f", $profiles->[$res_i]->[$j]);
        }
        $vall_lines .= "\n";                                          
    }

    return $vall_lines;
}


# ss3state()
# Changes the six-state DSSP alphabet for 2' structure into a 
# three-state alphabet. Here's how it works:
# Input         Output
#   H             H
#   I             H
#   G             H
#   E             B
#   B             B
#   L             L
#
# In the shortened alphabet, H denotes a residue that is part of
# an alpha helix, B denotes a residue that is part of a beta sheet,
# and L denotes a residue that is part of loop.

sub ss3state {
    my $ss = shift;

    my $letter;
    if ($ss eq 'H' || $ss eq 'I' || $ss eq 'G') {
        $letter = 'H';
    }
    elsif ($ss eq 'E' || $ss eq 'B') {
        $letter = 'E';
    } else {
        $letter = 'L';
    }
        
    if ( !defined $letter ) { 
        abort("Error: can't determine structure for letter ($ss)\n");
    }

    return $letter;
}

sub abort {
    my $msg = shift;
    print STDERR "$0: $msg\n";
    exit -2;
}

sub slurp_file {
    my $file = shift;
    open FILE, "<$file" or die $!;
    my @buffer = <FILE>;
    close FILE or die $!;

    return \@buffer;
}
