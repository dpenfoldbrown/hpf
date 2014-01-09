c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 3384 $
c  $Date: 2003-09-10 19:35:07 -0400 (Wed, 10 Sep 2003) $
c  $Author: rohl $

      subroutine all_ss 
      implicit none 
c--------------------------------------------------------------
      include 'param.h'
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

