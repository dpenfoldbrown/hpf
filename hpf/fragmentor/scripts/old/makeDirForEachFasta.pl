#!/usr/bin/perl -w

# Richard Bonneau , Initl/ auth.
# copywrite 2002 -  Institue for Systems Biology, 
# Seattle, WA 98103

#example use #($opt_hash,$args) = &my_getopts("o:");   ## left with arg, right logical
#example use #my(@file_list) = @$args;
#example use #if ($opt_hash->{opt_o}) {$outfile = $opt_hash->{opt_o} } else {die "-o outfile\n";}

@file = @ARGV;
print "file @file\n";
foreach $f (@file) {
    print "f $f\n";
	$baseName = $f;
	$baseName =~ s/\.fasta//g;
	mkdir "$baseName/";
	`mv $f $baseName/`;	
}








#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
sub my_getopts {
# general purpose command line parser
# operates on global @ARGV
# call with string 'xyz:pdq' plus optional code to designate how to handle errors.
# where x y and z are command line switches that take arguments
# and p,d, &q are ones that do not take args.
# e.g. suppose contents of ARGV resulted from:
#  -x first_arg -pd -y second_arg not_an_arg also_not_arg -q -z third_arg  again_not_arg
# then return value is a hash containing:
# x => first_arg, y=>second_arg, z=>third_arg, p=>1,d=>1,q=>1
# plus an array containing (not_an_arg,also_not_arg,again_not_arg)
# $err flags if error occurs such as missing args
     my($argumentative,$ignore_err) = @_;
    local( $_,$first,$rest);
    my $err = 0;
    my ($with_args,$without_args) = split /:/ , $argumentative;
    my $wal = length $with_args;
    my $i=0;
    my (%arg_hash,@not_args);
     while ( $i<=$#ARGV ) {
      if ( $ARGV[$i] =~ /^-(.)(.*)/ ) # does it start with a dash?
         {
            ($first,$rest) = ($1,$2);
             $pos = index($argumentative,$first);
           if($pos < 0)
             {  # did not match any arg ?
                if ($ignore_err<2)
                 {  #
                        print STDERR "Switch -$first  no such switch \n";
                        $err = -1;  # a soft error
                        last unless ($ignore_err) ;
                 }
                $pos = $wal+1  #  pretend it is a non-argument type flag
             }
           if($pos <= $wal )
             {  # needs an argument
                if (!$rest) #if $rest blank then arg is the next ARGV element
                  {
                        $i++; #will absorb this arg
                        if ($i>$#ARGV)
                           {  # got any more args?
                              print STDERR "Switch -$first requires an argument \n";
                             $err=1;
                             last unless ($ignore_err >1) ;
                             $ARGV[$i] = 1;  # make up the arg for it}
                           }
                          $rest = $ARGV[$i];
                   }
              }
                else # does not need an arg
              {
                   if ($rest) # if more tags then recylce the argument remainder
                     {
                        $ARGV[$i--] = "-$rest";
                     }
                       $rest = 1;
              }
              $arg_hash{"opt_$first"} = $rest;   # remember the arg
         } # did not start with dash
           else # is not an arg, so save it
         {
               push @not_args, $ARGV[$i];
         }
           $i++;  # next arg
       } #next while
     return (\%arg_hash, \@not_args,$err);
} #end sub



