#!/bin/tcsh

unalias rm
foreach dd (??/)
  cd $dd
  foreach d (?????)

     echo "doing:"
     echo $d
     cd $d
     ~/.local/share/wIBM/tools/make_fragments.local.pl -nosam $d.fasta
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
  cd ../
end


