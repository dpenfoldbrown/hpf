c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 2338 $
c  $Date: 2002-08-09 18:23:11 -0400 (Fri, 09 Aug 2002) $
c  $Author: rohl $

      subroutine convert1_num                
      implicit none 
c-------------------------------------------------------
      include 'param.h'
      include 'structure.h' 
      integer i              
c------------------------------------------------------
      do i=1,total_residue 
         seqnum(i)=0 
      enddo 
      do i=1,total_residue 
         if (seq(i).eq.'A') then 
            seqnum(i)=1 
         endif 
         if (seq(i).eq.'C') then 
            seqnum(i)=2 
         endif 
         if (seq(i).eq.'D'.or.seq(i).eq.'B') then 
            seq(i)='D' 
            seqnum(i)=3 
         endif 
         if (seq(i).eq.'E'.or.seq(i).eq.'Z') then 
            seq(i)='E' 
            seqnum(i)=4 
         endif 
         if (seq(i).eq.'F') then 
            seqnum(i)=5 
         endif 
         if (seq(i).eq.'G') then 
            seqnum(i)=6 
         endif 
         if (seq(i).eq.'H') then 
            seqnum(i)=7 
         endif 
         if (seq(i).eq.'I') then 
            seqnum(i)=8 
         endif 
         if (seq(i).eq.'K') then 
            seqnum(i)=9 
         endif 
         if (seq(i).eq.'L') then 
            seqnum(i)=10
         endif 
         if (seq(i).eq.'M') then 
            seqnum(i)=11
         endif 
         if (seq(i).eq.'N') then 
            seqnum(i)=12
         endif 
         if (seq(i).eq.'P') then 
            seqnum(i)=13
         endif 
         if (seq(i).eq.'Q') then 
            seqnum(i)=14
         endif 
         if (seq(i).eq.'R') then 
            seqnum(i)=15
         endif 
         if (seq(i).eq.'S') then 
            seqnum(i)=16
         endif 
         if (seq(i).eq.'T') then 
            seqnum(i)=17
         endif 
         if (seq(i).eq.'V') then 
            seqnum(i)=18
         endif 
         if (seq(i).eq.'W') then 
            seqnum(i)=19
         endif 
         if (seq(i).eq.'Y') then      
            seqnum(i)=20
         endif 
      enddo 
      return 
      end 
