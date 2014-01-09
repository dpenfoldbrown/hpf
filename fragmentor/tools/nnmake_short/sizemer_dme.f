c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 1.3 $
c  $Date: 2002/06/18 04:41:53 $
c  $Author: rohl $

      subroutine sizemer_dme(size)  
      implicit none 
c-------------------------------------------------------
      include 'structure.h' 
      integer i,j,m,n,place,good_dme,total_dme,stat_neigh,
     #        coverage,coverage5,new_pos,size,tot_pair    
      real dist1,dist2(0:max_len,0:max_len),dme,cutoff  
      integer dummy
c-------------------------------------------------------
      stat_neigh=25       !compute statistics over 1st stat_neigh in list
      stat_neigh=min(max_nn,stat_neigh) 
      if (size.le.1) return 
      total_dme=0 
      good_dme=0 
      tot_pair=0  
      coverage=0 
      coverage5=0 
      do i=1,size 
         tot_pair=tot_pair+i 
      enddo 
   
      if (size.ge.2)cutoff=0.12
      if (size.ge.4) cutoff=0.5
      if (size.ge.6) cutoff=0.8
      if (size.ge.8) cutoff=1.0
      if (size.ge.12) cutoff=1.3
      if (size.ge.16) cutoff=1.6
      if (size.ge.18) cutoff=1.8
      if (size.ge.20) cutoff=2.1
      write(0,*) " sizemer ",size,cutoff
      do i=1,total_residue-size
         new_pos=0 
         do m=0,size
            do n=m+1,size
               dist2(m,n)=(calpha(1,i+m)-calpha(1,i+n))**2  
     #                   +(calpha(2,i+m)-calpha(2,i+n))**2  
     #                   +(calpha(3,i+m)-calpha(3,i+n))**2  
               dist2(m,n)=sqrt(dist2(m,n)) 
               dist2(n,m)=dist2(m,n) 
            enddo 
         enddo 
         do j=1,stat_neigh  !collect stats over first stat_neigh 
            total_dme=total_dme+1 
            dme=0.0 
            place=best_nn(size,i,j)  
            do m=0,size
               do n=m+1,size
                  dist1=(vall_ca(1,place+m)-vall_ca(1,place+n))**2  
     #                 +(vall_ca(2,place+m)-vall_ca(2,place+n))**2  
     #                 +(vall_ca(3,place+m)-vall_ca(3,place+n))**2  
                  dme=dme+(sqrt(dist1)-dist2(m,n))**2   
               enddo 
            enddo 
            dme=sqrt(dme/tot_pair)   
            if (dme.lt.cutoff) then 
               good_dme=good_dme+1 
               if (new_pos.eq.0) coverage=coverage+1
               new_pos=new_pos+1 
            endif 
         enddo    
         if (new_pos.ge.5) coverage5=coverage5+1
      enddo 
      dummy = size+1
      if (total_dme.gt.0) then 
      write(60,'(I2,A3,a6,2i8,f8.5)')dummy,'mer',' tot: ',good_dme,
     #        total_dme,real(good_dme)/total_dme 
      write(60,'(I2,A3,a6,i8,f8.5)')dummy,'mer',' one: ',coverage,
     #     real(coverage*min(max_nn,stat_neigh))/total_dme
      write(60,'(I2,A3,a7,i8,f8.5)')dummy,'mer',' five: ',coverage5,
     #     real(coverage5*min(max_nn,stat_neigh))/total_dme
      endif 

      return 
      end 

