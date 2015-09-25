#!/bin/tcsh

foreach d ( ?? )
  ../tools/move_complete_workUnits.pl -o part_done/ -i $d
end

