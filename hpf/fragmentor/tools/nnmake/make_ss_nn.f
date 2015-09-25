c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 2338 $
c  $Date: 2002-08-09 18:23:11 -0400 (Fri, 09 Aug 2002) $
c  $Author: rohl $

      subroutine make_ss_nn 

car formerly required 8mers

      implicit none 
c------------------------------------------------------------
      include 'param.h'
      include 'structure.h'
      integer i,j,helix,strand,loop,place,start  
      integer size
c------------------------------------------------------------

      size=sizes(1)    
      do i=1,n_sizes
         if (sizes(i).eq.8) size=8    !8-mers  if available
      enddo
      
      do i=1,total_residue 
         start=i-4 
         if (start.lt.1) start=1 
         if (start.gt.total_residue-size) start=total_residue-size 
         helix=0 
         strand=0 
         loop=0 
         do j=1,min(25,max_nn)    
            place=best_nn(size,start,j)+(i-start)  
            if (ss(place).eq.'H') helix=helix+1 
            if (ss(place).eq.'E') strand=strand+1 
            if (ss(place).eq.'L') loop=loop+1 
         enddo 
         ss_nn(i)='L' 
         if (helix.gt.strand.and.helix.gt.loop) ss_nn(i)='H' 
         if (strand.gt.helix.and.strand.gt.loop) ss_nn(i)='E' 
      enddo 

      return 
      end 

