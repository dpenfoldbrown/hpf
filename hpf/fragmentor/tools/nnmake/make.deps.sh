#  CVS information:
#  $Revision: 1320 $
#  $Date: 2002-03-26 21:11:52 -0500 (Tue, 26 Mar 2002) $
#  $Author: rohl $

# make.deps.sh: make a list of fortran dependency files for make
# It works by finding all the include statements and picking out the filenames
# usage: make.deps.sh source.f .source.d source.o
# JJG 4/3/1

echo making $2 from $1
echo $2 .'$(COMPILER)'.$3: \\ > $2
grep "include ['\"]" $1 | sed -e 's/^[[:alpha:]].*$//g' -e 's/include //g' -e 's/"//g' -e "s/'//g" -e 's/  */ /g' -e 's/!.*$//g' -e 's/$/ \\/g' >> $2
echo \ $1 >> $2
echo >> $2

#cat $2
