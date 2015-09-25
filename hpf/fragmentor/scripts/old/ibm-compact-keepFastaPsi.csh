#!/bin/tcsh

unalias rm
foreach dd ( ??    )
  cd $dd
  foreach d (?????/)
     echo "doing:"
     echo $d
     cd $d
     mv *.fasta* ../
     mv *.psipred* ../
     cd ../
     rm -vfr $d
  end
  cd ../
end


