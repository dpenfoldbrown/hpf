#!/bin/tcsh

unalias rm
echo pwd: $PWD
foreach d (?????)

     echo "doing:"
     echo $d
     cd $d
     rm -f ss_blast
     rm -f sstmp.*
     rm -f *.check
     rm -f *.checkpoint
     rm -f *.dat
     rm -f *.jufo_ss
     rm -f names.*
     bzip2 -9 *_v1_3
     cd ../

end


