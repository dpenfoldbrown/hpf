c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 1.3 $
c  $Date: 2002/06/18 04:41:53 $
c  $Author: rohl $

      real  function dme_native(it,ivall,size)  
      implicit none 
c-----------------------------------------------------------------
      include 'structure.h' 
      integer i,m,n,it,ivall,size,tot_pair    
      real dist1,dist2,dme
c-----------------------------------------------------------------
      
      dme_native=99.0
      dme=0.0 
      tot_pair=0 
      do i=1,size
         tot_pair=tot_pair+i
      enddo
      if (tot_pair.gt.0) then 
         do m=0,size
            do n=m+1,size
               dist2=(calpha(1,it+m)-calpha(1,it+n))**2  
     #              +(calpha(2,it+m)-calpha(2,it+n))**2  
     #              +(calpha(3,it+m)-calpha(3,it+n))**2  
               dist1=(vall_ca(1,ivall+m)-vall_ca(1,ivall+n))**2  
     #              +(vall_ca(2,ivall+m)-vall_ca(2,ivall+n))**2  
     #              +(vall_ca(3,ivall+m)-vall_ca(3,ivall+n))**2  
               dme=dme+(sqrt(dist1)-sqrt(dist2))**2   
            enddo 
         enddo 
         dme_native=sqrt(dme/tot_pair)
      endif

      return 
      end 

c-----------------------------------------------------------------

      function dme_full(it,ivall,size)  
      implicit none 
      include 'structure.h' 
      integer i,m,n,it,ivall,count,size,tot_pair    
      real dist1,dist2,dme,dme_full 

      dme_full = 99.0
      dme=0.0 
      count=0
      tot_pair=0 
      do i=1,size
         tot_pair=tot_pair+i
      enddo
      if (tot_pair.gt.0) then 
         do m=0,size
          do n=m+1,size
          count = count+1
          dist2=(calpha(1,it+m)-calpha(1,it+n))**2  
     #         +(calpha(2,it+m)-calpha(2,it+n))**2  
     #         +(calpha(3,it+m)-calpha(3,it+n))**2  
          dist1=(vall_ca(1,ivall+m)-vall_ca(1,ivall+n))**2  
     #         +(vall_ca(2,ivall+m)-vall_ca(2,ivall+n))**2  
     #         +(vall_ca(3,ivall+m)-vall_ca(3,ivall+n))**2  
          dme=dme+(sqrt(dist1)-sqrt(dist2))**2 
c dme is the c_alpha dme,  also include other atoms in dem for dme_full  
           dist2=(native_cnh(1,it+m)-native_cnh(1,it+n))**2  
     #         +(native_cnh(2,it+m)-native_cnh(2,it+n))**2  
     #         +(native_cnh(3,it+m)-native_cnh(3,it+n))**2  
          dist1=(cnh(1,ivall+m)-cnh(1,ivall+n))**2  
     #         +(cnh(2,ivall+m)-cnh(2,ivall+n))**2  
     #         +(cnh(3,ivall+m)-cnh(3,ivall+n))**2 
     
          dme=dme +(sqrt(dist1)-sqrt(dist2))**2
          dist2=(native_cnh(1,it+m)-calpha(1,it+n))**2  
     #         +(native_cnh(2,it+m)-calpha(2,it+n))**2  
     #         +(native_cnh(3,it+m)-calpha(3,it+n))**2  
          dist1=(cnh(1,ivall+m)-vall_ca(1,ivall+n))**2  
     #         +(cnh(2,ivall+m)-vall_ca(2,ivall+n))**2  
     #         +(cnh(3,ivall+m)-vall_ca(3,ivall+n))**2 
     
          dme=dme +(sqrt(dist1)-sqrt(dist2))**2 

          dist2=(native_cnh(1,it+n)-calpha(1,it+m))**2  
     #         +(native_cnh(2,it+n)-calpha(2,it+m))**2  
     #         +(native_cnh(3,it+n)-calpha(3,it+m))**2  
          dist1=(cnh(1,ivall+n)-vall_ca(1,ivall+m))**2  
     #         +(cnh(2,ivall+n)-vall_ca(2,ivall+m))**2  
     #         +(cnh(3,ivall+n)-vall_ca(3,ivall+m))**2 

          dme=dme +(sqrt(dist1)-sqrt(dist2))**2 

            enddo 
         enddo 
           dme_full = sqrt( dme / (4*tot_pair))
c$$$         if (dme_full.gt.30) then
c$$$          do m=0,size
c$$$          do n=m+1,size
c$$$          dist2=(calpha(1,it+m)-calpha(1,it+n))**2  
c$$$     #         +(calpha(2,it+m)-calpha(2,it+n))**2  
c$$$     #         +(calpha(3,it+m)-calpha(3,it+n))**2  
c$$$          dist1=(vall_ca(1,ivall+m)-vall_ca(1,ivall+n))**2  
c$$$     #         +(vall_ca(2,ivall+m)-vall_ca(2,ivall+n))**2  
c$$$     #         +(vall_ca(3,ivall+m)-vall_ca(3,ivall+n))**2  
c$$$          d1 = (sqrt(dist1)-sqrt(dist2))
c$$$c dme is the c_alpha dme,  also include other atoms in dem for dme_full  
c$$$           dist2=(native_cnh(1,it+m)-native_cnh(1,it+n))**2  
c$$$     #         +(native_cnh(2,it+m)-native_cnh(2,it+n))**2  
c$$$     #         +(native_cnh(3,it+m)-native_cnh(3,it+n))**2  
c$$$          dist1=(cnh(1,ivall+m)-cnh(1,ivall+n))**2  
c$$$     #         +(cnh(2,ivall+m)-cnh(2,ivall+n))**2  
c$$$     #         +(cnh(3,ivall+m)-cnh(3,ivall+n))**2 
c$$$     
c$$$          d2= (sqrt(dist1)-sqrt(dist2))
c$$$          dist2=(native_cnh(1,it+m)-calpha(1,it+n))**2  
c$$$     #         +(native_cnh(2,it+m)-calpha(2,it+n))**2  
c$$$     #         +(native_cnh(3,it+m)-calpha(3,it+n))**2  
c$$$          dist1=(cnh(1,ivall+m)-vall_ca(1,ivall+n))**2  
c$$$     #         +(cnh(2,ivall+m)-vall_ca(2,ivall+n))**2  
c$$$     #         +(cnh(3,ivall+m)-vall_ca(3,ivall+n))**2 
c$$$     
c$$$          d3 = (sqrt(dist1)-sqrt(dist2)) 
c$$$ 
c$$$        write(0,*) m,n,it,vall
c$$$        write(0,*) d1,d2,d3
c$$$        write(0,*) dist1,dist2
c$$$     	write(0,*) ( calpha(i,it+m),i=1,3), ( calpha(i,it+n),i=1,3)
c$$$	write(0,*) ( native_cnh(i,it+m),i=1,3), ( native_cnh(i,it+n),i=1,3)
c$$$	write(0,*) ( vall_ca(i,ivall+m),i=1,3), ( vall_ca(i,ivall+n),i=1,3)
c$$$        write(0,*) ( cnh(i,ivall+m),i=1,3), ( cnh(i,ivall+n),i=1,3)
c$$$        enddo
c$$$        enddo
c$$$        write(0,*) '*******************************************************'
c$$$        write(0,*) tot_pair
c$$$        pause
c$$$        endif
      endif
c$$$         if (count.ne.tot_pair) pause
      return 
      end 
    
