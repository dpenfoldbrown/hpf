#!/bin/tcsh

foreach d (  ?? )
  tar cvf $d.000-999.tar $d
  bzip2 -9 $d.000-999.tar
  md5sum $d.000-999.tar.bz2 > $d.000-999.md5
  mv $d.000-999.md5 ~grid/incoming/
  mv $d.000-999.tar ~grid/incoming/
end

