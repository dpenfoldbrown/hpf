c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 3384 $
c  $Date: 2003-09-10 19:35:07 -0400 (Wed, 10 Sep 2003) $
c  $Author: rohl $

      subroutine make_dat   
      implicit none
      include 'path_defs.h'
c----------------------------------------------------------------------
      include 'param.h'
      include 'structure.h'
      integer i,iunit
c----------------------------------------------------------------------

      iunit=16
      write(0,*) 'make_dat',iunit,make_dat_path(1:make_dat_l)//
     #     lower_name//'.dat'
      open(iunit,file=make_dat_path(1:make_dat_l)//
     #     lower_name//'.dat',status='unknown') 
      do i=1,total_residue
         write(iunit,50)24,total_residue-i,i,i,
     #        0.0,0.0,0.0, 
     #        0.0,0.0,0.0,0.0,  
     #        lower_name,chain_letter,seq(i)
      enddo
      close(iunit)  
      
 50   format(i2,1x,3i4,1x,3f9.3,4f8.3,1x,a4,a1,1x,a1)
      return
      end 


