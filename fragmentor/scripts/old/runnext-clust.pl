#!/usr/bin/perl -w


# needs 2 things:
# the node number or descriptor
# and a list of rosetta format directories
## you should be in the group dir , one up, when
## running this ... the rosetta dir should have a rosetta
## name ... ie. xxxxc/

#open(IFILE,"$list")  || die "can't find dir/file list\n\n";

$orgdir = `pwd`;
chomp $orgdir;

#$up = `uptime`;
#chomp $up;
#@tmp=split /\s+/, $up; 
#$iup=$tmp[$#tmp-2];
#$iup =~ s/,//g;
#$b = sprintf "%.0f", $iup;
#if ($b > 5) { 
#	$checkMSG = "node: $nodeNum  load: $b    ... (not run , load higher than max)\n";
#	`echo $checkMSG >> $check`;
#	die "load to high of this host ( node: $nodeNum  load: $b ) \n";
#}

foreach $dd (@ARGV) {
    $d = substr $dd,0,4;
    $c = substr $dd,4,1;

    `mkdir $d$c `;
    `mv $dd $d$c`;
    chdir "$d$c";    


    `gunzip -f *.out.gz`;
    $dd =~ s/\.gz//g;
        #`/home/rbonneau/bin/silentCombine.pl 15 [0-9]*.out >  ao$d.out ` ;
    system("/bin/nice -19 /home/rbonneau/utils_pred/philClust/cluster_info_silent.out $dd - $d$c 20,75,100 3,7 co,5 0.80 ");
        #system("/bin/nice -19 /home/rbonneau/utils_pred/philClust/cluster_info_silent.out ao$d.out - $d$c 25,75,100 3,7 >& log-noCO-clust");
        #`gzip -f *.subcluster.*.pdb *.cluster*.contacts`;
    `gzip *.pdb`;
    `gzip *.out`;
        #`gzip log-noCO-clust`;
    chdir "$orgdir";

}

