c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 2336 $
c  $Date: 2002-08-09 15:02:53 -0400 (Fri, 09 Aug 2002) $
c  $Author: rohl $

      subroutine file_manager()
      implicit none
      include 'path_defs.h'
c *********
c charlie e m strauss
c
c  reads in file paths and stores them in a common block.
c  the order of pathnames in the file is critical!!!

c  the order is determined by the order the variables in the COMMON blocks in 
c path_defs.h.
c  the variables have plain english names so their interpretatio should be 
c self evident


c  Here we use a clever trick to automate the reading in of the
c  pathnames without having to know the variable names explicitly
c  we use an equivalence statment between an array and the first variable
c in the common block 

c note: it is ASSUMED that the paths do not contain a "space" in them.  
c If you want paths
c that contain spaces you will need to edit the automatice string length
c  calculation below
c note: paths should be less than 132 chars in length and end with a slash
     
         
      integer j
      character*132 names(MAX_PATHS)
      
      equivalence (loop_path,names(1))
      integer name_length(MAX_PATHS)
      equivalence (loop_path_l, name_length(1) )

      integer iunit
      integer vall_l

      iunit=22

      open(unit=iunit, file='path_defs.txt',
     #     status='old',iostat=j)

      if (j.ne.0) then
         write(0,*) "This program needs to be run in the directory "
     #        //"containing path_defs.txt"
         stop
      endif

      do j=1,MAX_PATHS 
         read(iunit,'(a132)',END=999, err=9999) names(j)
         name_length(j) = index(names(j),' ')-1
         if (name_length(j).gt.0) then
            write(0,*) j, names(j)(1:name_length(j)) 
         else 
            write(0,*) j, 'no path '
            name_length(j)=1
         endif
      enddo
      read(iunit,'(a50)',END=998,err=9999)vall_name

      if (index(vall_name,'vall.dat').ne.1) then
         write(0,*)"non-standard vall name: ",vall_name
         write(0,*)"expecting vall.dat.<identifier>"
         stop
      endif

      vall_l=index(vall_name,' ')
      vall_cst_coord_name = 'vall_cst_coord'//vall_name(5:vall_l)

      close(iunit)
      return
 998  continue
      write(0,*)"vall name not included in path_defs.txt"
      stop
 999  continue
      write(0,*)'expecting ', MAX_PATHS, ' paths'
      write(0,*)'only found ',j-1, ' paths in path_defs.txt'
      stop
      
 9999 continue
      write(0,*) "dang, error reading path_defs.txt file"
      stop
      
      end

           


