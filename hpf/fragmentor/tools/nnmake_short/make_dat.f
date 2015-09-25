c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 1.3 $
c  $Date: 2002/06/18 04:41:53 $
c  $Author: rohl $

      subroutine make_dat   
      implicit none
      include 'path_defs.h'
c----------------------------------------------------------------------
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
     #        lower_name,chain_letter,seq(i),nalign
      enddo
      close(iunit)  
      
 50   format(i2,1x,3i4,1x,3f9.3,4f8.3,1x,a4,a1,1x,a1,1x,i3)
      return
      end 


