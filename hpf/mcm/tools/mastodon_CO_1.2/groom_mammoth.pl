#!/usr/bin/perl -w
# cems 2001  copyright (C) charlie e.m. strauss,2001
# check files experimental CA traces of sinlge chains for
# consecutive complete unambiguous and properly located calphas
# chains failing this test are not suited to rosetta, nnmake, or mammoth.
# run with single arg on command line: name of pdb file
# errors are marked with the line number and residue number.
my $file;
my (@errs);
use strict;
foreach $file (@ARGV) {
open FH,$file || do { warn "cant open $file .. will skip \n";next;};
print STDERR "Processing $file ";

my @CA = grep { /^ATOM / && substr($_,12,3) eq " CA"  } <FH>;
print STDERR "Number of CA ATOM records: $#CA\n";
close FH;
my @pdb = ();
my $line = shift @CA;
push @pdb,$line;
my %h = ();

my @err=(); 
my $d;

my $alt = substr($line,16,1);
my $chain = substr($line,21,1);
my $seq=substr($line,22,4);
my $acode=substr($line,26,1);
my $x = substr($line,30,8);
my $y = substr($line,38,8);
my $z = substr($line,46,8);
my $too_close;
my $i = 0;
if ($alt ne " ") {chomp($line);
        push @errs,$line,"$i $seq altloc "};

my ($l_x,$l_y,$l_z,$l_chain,$l_seq) = ($x,$y,$z,$chain,$seq);
 $h{$seq} .=$acode;
foreach $line (@CA) {
$i++;
$too_close=0;
$alt = substr($line,16,1);
$chain = substr($line,21,1);
$seq=substr($line,22,4);
$acode=substr($line,26,1);
$x = substr($line,30,8);
$y = substr($line,38,8);
$z = substr($line,46,8);

if ($alt ne " ") {push @err,"$i $seq altloc $alt\n"};
if ($chain ne $l_chain) {

       push @err,"$i $seq chain break :$chain: != :$l_chain: \n";}
else {

 if (  exists  $h{$seq} && $acode =~ $h{$seq} ) {
     
        push @err, "$i $seq repeated residue ";
        }
 else {
   if ($seq != $l_seq+1 && $acode eq " ") {
      push  @err,"$i $seq non-consecutive residue: previous $l_seq\n"; }
   else {
    $d = ($x-$l_x)*($x-$l_x)+
      ($y-$l_y)*($y-$l_y)+ ($z-$l_z)*($z-$l_z);
      $d = sqrt($d) if $d>0;


   if (($d < 1.0) || ($d > 5.5)) {
      push @err, "$i $seq bad chain ca separation $d ".
        ": previous xyz: $l_x $l_y $l_z \n";
        $too_close= $d<1 ? 1:0;};
   }
  }
  }
push  @pdb,$line  unless $too_close ||(exists $h{$seq} && $acode =~ $h{$seq} && $alt ne " ") ;

($l_x,$l_y,$l_z,$l_chain,$l_seq) = ($x,$y,$z,$chain,$seq);
    $h{$seq} .=$acode;
if (@err) { chomp($line);push @errs,"FILE $file line $line \n",@err};
@err=();


}
open FH,">$file" || die "cant open $file";
print FH @pdb;
close FH;
}

if (@errs) {
  print "errors exist  \n";
  print join "\n",@errs,"\n";
} else
{ print "No Errors !!!\n" }

