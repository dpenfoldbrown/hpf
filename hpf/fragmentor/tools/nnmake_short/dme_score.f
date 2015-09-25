c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 1.3 $
c  $Date: 2002/06/18 04:41:53 $
c  $Author: rohl $

      real function dme_score(itarget,ivall,size,firstnum)  
      implicit none 
c-------------------------------------------------------
      include 'structure.h' 
      integer i,j,m,n,place,itarget,ivall,size,tot_pair,firstnum     
      real dist1,dist2(0:max_len,0:max_len),dme
c-------------------------------------------------------
      dme_score=99.0
      if (size.le.1) return 
      tot_pair=0 
      do i=1,size
         tot_pair=tot_pair+i
      enddo
      if (tot_pair.gt.0) then 
         do m=0,size
            do n=m+1,size
               dist2(m,n)=(vall_ca(1,ivall+m)-vall_ca(1,ivall+n))**2  
     #                   +(vall_ca(2,ivall+m)-vall_ca(2,ivall+n))**2  
     #                   +(vall_ca(3,ivall+m)-vall_ca(3,ivall+n))**2  
               dist2(m,n)=sqrt(dist2(m,n)) 
            enddo 
         enddo 
         do j=1,firstnum
            dme=0.0 
            place=best_nn(size,itarget,j)  
            do m=0,size
               do n=m+1,size
                  dist1=(vall_ca(1,place+m)-vall_ca(1,place+n))**2  
     #                 +(vall_ca(2,place+m)-vall_ca(2,place+n))**2  
     #                 +(vall_ca(3,place+m)-vall_ca(3,place+n))**2  
                  dme=dme+(sqrt(dist1)-dist2(m,n))**2   
               enddo 
            enddo 
            dme=sqrt(dme/tot_pair)   
            if (dme.lt.dme_score) dme_score=dme
         enddo 
      endif

      return 
      end 

