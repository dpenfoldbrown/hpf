c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 1.5 $
c  $Date: 2002/06/18 04:41:53 $
c  $Author: rohl $

      subroutine all_ss 
      implicit none 
c--------------------------------------------------------------
      include 'structure.h' 
      integer i,j,k  
      character id*6
c--------------------------------------------------------------
      k=total_residue 
      write(60,*) 
      do i=1,k,60 
         id='seq: '
         write(60,1000)id,(seq(j),j=i,min(i+59,k)) 
         write(60,*) 
         id='nativ:'
         write(60,1000)id,(ss_native(j),j=i,min(i+59,k)) 
         id='phd:  '
         write(60,1000)id,(ss_phd(j),j=i,min(i+59,k)) 
         id='jones:'
         write(60,1000)id,(ss_jon(j),j=i,min(i+59,k)) 
         id='rdb:  '
         write(60,1000)id,(ss_rdb(j),j=i,min(i+59,k)) 
         id='jufo: '
         write(60,1000)id,(ss_jufo(j),j=i,min(i+59,k)) 
         id='nn:   '
         write(60,1000)id,(ss_nn(j),j=i,min(i+59,k)) 
         write(60,*) 
      enddo 
  
      return 
 1000 format(a6,1x,60a1)
      end 

