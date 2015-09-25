#!/bin/tcsh

foreach d ( gi gu gw hb hc he )
  tar cvf $d.000-999.tar $d
  bzip2 $d.000-999.tar
  md5sum $d.000-999.tar.bz2 > $d.000-999.tar.md5
  mv $d.000-999.tar.* ~grid/incoming/
end

