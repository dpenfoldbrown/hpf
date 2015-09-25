       
       ! author: copyright 2001 Charlie E. M. Strauss 2001
       
       subroutine rms_setup(n,xe,xp,xm,xre,xrp,xse,xsp)
       !author: charlie E. M. Strauss 2001
       ! precomputes moments and cross moments of two aligned point vectors
       ! this routine is just a middle man calling others to do its work
       !  in particular I note that if you are going to be calling this with 
       ! various pairings of the
       ! same arrays, then its only the xm array that is coing to change, and 
       !  you might want
       ! to rewrite this for better effiency.
       
       
       implicit none
       ! inputs
       integer n
       real*8 xe(3,n),xp(3,n)
       ! outputs
       real*8 xre(0:n),xrp(0:n)
       real*8 xse(3,0:n),xsp(3,0:n)
       real*8 xm(3,3,0:n)

       
       call  prepare_cross_moments(n,xe,xp,xm)
       call  prepare_rsums(n,xe,xse)
       call  prepare_rsums(n,xp,xsp)
       call  prepare_radii(n,xe,xre)
       call  prepare_radii(n,xp,xrp)
       
       
       
       return
       end
       
       
       ! the headache here is dealing with the center of mass in a sliding way.
  
         real*8 function  rms2_fast2(ii,jj,xm,xre,xrp,xse,xsp,det,ev)
c   AUTHOR: charlie strauss 2001
c   computes the rms between two non-weighted point vectors.
c   intended to be applied to a windowed region of a large array of points.
c
c
c   the vectors are not handed in.
c   instead their pre-computed sliding moments are handed in
c   plus the start(ii) and stop(jj) atoms defining the window in the larger 
c   array
c   xm is the sliding cross moments array
c   xre is the sliding radius array for the first set of points
c   xrp is the sliding radius array for the second set of points
c   xse is the sliding sum_coord array for first
c   xsp is the sliding sum_coord array for the second
c
c  assumes caller is not an idiot:  jj>ii>0

c   det is an out value of the determinant of the cross moment matrix
c   returned value is the rms
c
c   most of this is double precision for good reasons.  first there are some 
c   large
c   differences of small numbers.  and second the rsymm_eignen() routine can 
c   internally have numbers
c   larger than the largest real*4 number.  (you could do some fancy foot work 
c   to rescale things if
c   you really had a problem with this.
c
c   
c    
	   implicit none
       integer ii,jj
       real*8 xm(3,3,0:jj)             ! note zero origin
       real*8 xre(0:jj),xrp(0:jj)      ! note zero origin
       real*8 xse(3,0:jj),xsp(3,0:jj)  ! note zero origin

   
       real*8 det

       real*8 ev(3)
       real*8 rms2_fast3
       integer n,kk
       kk = ii-1
       n= jj-kk
       
        rms2_fast2 = rms2_fast3(n,xm(1,1,kk),xre(kk),xrp(kk),xse(1,kk),xsp(1,kk),
     &                   xm(1,1,jj),xre(jj),xrp(jj),xse(1,jj),xsp(1,jj))
       return
       end
       
       real*8 function  rms2_fast3(n,xmi,xrei,xrpi,xsei,xspi,
     &                               xmj,xrej,xrpj,xsej,xspj)
c   AUTHOR: charlie strauss 2001
c   computes the rms between two non-weighted point vectors.
c   intended to be applied to a windowed region of a large array of points.
c
c
c   the vectors are not handed in.
c   instead their pre-computed sliding moments are handed in
c   plus the start(ii) and stop(jj) atoms defining the window in the larger 
c   array
c   xm is the sliding cross moments array
c   xre is the sliding radius array for the first set of points
c   xrp is the sliding radius array for the second set of points
c   xse is the sliding sum_coord array for first
c   xsp is the sliding sum_coord array for the second
c
c  assumes caller is not an idiot:  jj>ii>0

c   det is an out value of the determinant of the cross moment matrix
c   returned value is the rms
c
c   most of this is double precision for good reasons.  first there are some 
c   large
c   differences of small numbers.  and second the rsymm_eignen() routine can 
c   internally have numbers
c   larger than the largest real*4 number.  (you could do some fancy foot work 
c   to rescale things if
c   you really had a problem with this.
c
c   
c    
	   implicit none
       integer n
       real*8 xmi(3,3)             ! note zero origin
       real*8 xrei,xrpi      ! note zero origin
       real*8 xsei(3),xspi(3)  ! note zero origin
       real*8 xmj(3,3)             ! note zero origin
       real*8 xrej,xrpj      ! note zero origin
       real*8 xsej(3),xspj(3)  ! note zero origin



       ! local
       real*8 det
       integer i,j,k
       real*8 temp1,invmass
       real*8 come(3),comp(3)
       real*8 ev(3)
       real*8 m_moment(3,3),rr_moment(3,3)
       real*8 det3
       real*8 rms_ctx,rms2_sum
       real*8 handedness
      
        
c compute center of mass 

        invmass = 1.0d0/(n)    ! 1 over the mass
        come(1) = (xsej(1)-xsei(1))*invmass       !x_com 
        come(2) = (xsej(2)-xsei(2))*invmass       !y_com 
        come(3) = (xsej(3)-xsei(3))*invmass      !z_com 
        
        comp(1) = (xspj(1)-xspi(1))*invmass
        comp(2) = (xspj(2)-xspi(2))*invmass 
        comp(3) = (xspj(3)-xspi(3))*invmass

C	Make cross moments matrix  

        do  k=1,3
        ! do  j = 1,3
             !  Ye Olde moments trick:   sum (x-x0)(y-y0) = sum (xy -x0y0) 
             m_moment(k,1) = (xmj(k,1)-xmi(k,1))*invmass- come(k)*comp(1)  ! flopped  
             m_moment(k,2) = (xmj(k,2)-xmi(k,2))*invmass- come(k)*comp(2) 
             m_moment(k,3) = (xmj(k,3)-xmi(k,3))*invmass- come(k)*comp(3)         
        ! enddo
       enddo

	    det = det3(m_moment) 

            if (Dabs(det) .le. 1e-24) then
               !write(0,*) 'Warning:degenerate cross moments: det=',det
                !returning a zero rms, to avoid any chance of
                ! Floating Point Errors.
            
             rms2_fast3=0.0d0
             return
            
             endif
             ! get handedness  of frame from determinant
         handedness = DSIGN(1.0D0,det) 
             ! weird but documented "feature" of DSIGN(a,b) (but not SIGN) is 
             ! that if fails if a<0
             
	    !  multiply cross moments by itself
         
	    do i = 1,3
	    do j = i,3
	      rr_moment(i,j) = m_moment(1,i)*m_moment(1,j)
     &                   + m_moment(2,i)*m_moment(2,j)
     &                   + m_moment(3,i)*m_moment(3,j)
          rr_moment(j,i) = rr_moment(i,j)
     	  enddo
        enddo
        
        
            !  compute eigen values of cross-cross moments

              call rsym_eigenval(rr_moment,ev)
              
              if (handedness.lt.0) then
              
              ! reorder eigen values  so that ev(3) is the smallest eigenvalue
              
              if (ev(2) .gt. ev(3) ) then
                    if (ev(3).gt.ev(1)) then
                        temp1 = ev(3)
                        ev(3) = ev(1)
                        ev(1) = temp1
                    endif
               else
                    if (ev(2).gt.ev(1)) then
                       temp1 = ev(3)
                       ev(3) = ev(1)
                       ev(1) = temp1
                     else
                       temp1 = ev(3)
                       ev(3) = ev(2)
                       ev(2) = temp1
                 endif
               endif
               endif  
                 ! ev(3) is now the smallest eigen value.  the other two are not
                 !  sorted.
                 
           
            ! now we must catch the special case of the rotation with inversion.
            ! we cannot allow inversion rotations.
            ! fortunatley, and curiously, the optimal non-inverted rotation 
            ! matrix 
            ! will have the similar eigen values.
            ! we just have to make a slight change in how we handle things 
            ! depending on determinant
            
             
            	rms_ctx = sqrt(abs(ev(1)))+sqrt(abs(ev(2)))
     &            + handedness*sqrt(abs(ev(3)))
  
             
                
            ! the abs() are theoretically unneccessary since the eigen values
            ! of a real symmetric
            ! matrix are non-negative.  in practice sometimes small eigen vals 
            ! end up just negative
            temp1 = come(1)**2+come(2)**2+come(3)**2+
     &              comp(1)**2+comp(2)**2+comp(3)**2
            
            rms2_sum = (xrej-xrei + xrpj-xrpi )*invmass -temp1
          
            ! and combine the outer and cross terms into the final calculation.
            !  (the abs() just saves us a headache when the roundoff error 
            ! accidantally makes the sum negative)
            
            rms2_fast3 = sqrt(abs(rms2_sum-2.0D0*rms_ctx))
            
            
          
         return
         end





       
       real*8 function  rms2_fast(ii,jj,xm,xre,xrp,xse,xsp,det,ev)
c   AUTHOR: charlie strauss 2001
c   computes the rms between two non-weighted point vectors.
c   intended to be applied to a windowed region of a large array of points.
c
c
c   the vectors are not handed in.
c   instead their pre-computed sliding moments are handed in
c   plus the start(ii) and stop(jj) atoms defining the window in the larger 
c   array
c   xm is the sliding cross moments array
c   xre is the sliding radius array for the first set of points
c   xrp is the sliding radius array for the second set of points
c   xse is the sliding sum_coord array for first
c   xsp is the sliding sum_coord array for the second
c
c  assumes caller is not an idiot:  jj>ii>0

c   det is an out value of the determinant of the cross moment matrix
c   returned value is the rms
c
c   most of this is double precision for good reasons.  first there are some 
c   large
c   differences of small numbers.  and second the rsymm_eignen() routine can 
c   internally have numbers
c   larger than the largest real*4 number.  (you could do some fancy foot work 
c   to rescale things if
c   you really had a problem with this.
c
c   
c    
	   implicit none
       integer ii,jj
       real*8 xm(3,3,0:jj)             ! note zero origin
       real*8 xre(0:jj),xrp(0:jj)      ! note zero origin
       real*8 xse(3,0:jj),xsp(3,0:jj)  ! note zero origin


       ! local
       real*8 det
       integer i,j,k
       real*8 temp1,invmass
       real*8 come(3),comp(3)
       real*8 ev(3)
       real*8 m_moment(3,3),rr_moment(3,3)
       real*8 det3
       real*8 rms_ctx,rms2_sum
       real*8 handedness
      
        
c compute center of mass 

        invmass = 1.0d0/(jj-ii+1.d0)    ! 1 over the mass
        come(1) = (xse(1,jj)-xse(1,ii-1))*invmass       !x_com 
        come(2) = (xse(2,jj)-xse(2,ii-1))*invmass       !y_com 
        come(3) = (xse(3,jj)-xse(3,ii-1))*invmass      !z_com 
        
        comp(1) = (xsp(1,jj)-xsp(1,ii-1))*invmass
        comp(2) = (xsp(2,jj)-xsp(2,ii-1))*invmass 
        comp(3) = (xsp(3,jj)-xsp(3,ii-1))*invmass

C	Make cross moments matrix  

        do  k=1,3
        ! do  j = 1,3
             !  Ye Olde moments trick:   sum (x-x0)(y-y0) = sum (xy -x0y0) 
             m_moment(k,1) = (xm(k,1,jj)-xm(k,1,ii-1))*invmass- come(k)*comp(1)  ! flopped  
             m_moment(k,2) = (xm(k,2,jj)-xm(k,2,ii-1))*invmass- come(k)*comp(2) 
             m_moment(k,3) = (xm(k,3,jj)-xm(k,3,ii-1))*invmass- come(k)*comp(3)         
        ! enddo
       enddo

	    det = det3(m_moment) 

            if (Dabs(det) .le. 1e-24) then
               !write(0,*) 'Warning:degenerate cross moments: det=',det
                !returning a zero rms, to avoid any chance of
                ! Floating Point Errors.
            
             rms2_fast=0.0d0
             return
            
             endif
             ! get handedness  of frame from determinant
         handedness = DSIGN(1.0D0,det) 
             ! weird but documented "feature" of DSIGN(a,b) (but not SIGN) is 
             ! that if fails if a<0
             
	    !  multiply cross moments by itself
         
	    do i = 1,3
	    do j = i,3
	      rr_moment(i,j) = m_moment(1,i)*m_moment(1,j)
     &                   + m_moment(2,i)*m_moment(2,j)
     &                   + m_moment(3,i)*m_moment(3,j)
          rr_moment(j,i) = rr_moment(i,j)
     	  enddo
        enddo
        
        
            !  compute eigen values of cross-cross moments

              call rsym_eigenval(rr_moment,ev)
              
              if (handedness.lt.0) then
              
              ! reorder eigen values  so that ev(3) is the smallest eigenvalue
              
              if (ev(2) .gt. ev(3) ) then
                    if (ev(3).gt.ev(1)) then
                        temp1 = ev(3)
                        ev(3) = ev(1)
                        ev(1) = temp1
                    endif
               else
                    if (ev(2).gt.ev(1)) then
                       temp1 = ev(3)
                       ev(3) = ev(1)
                       ev(1) = temp1
                     else
                       temp1 = ev(3)
                       ev(3) = ev(2)
                       ev(2) = temp1
                 endif
               endif
               endif  
                 ! ev(3) is now the smallest eigen value.  the other two are not
                 !  sorted.
                 
           
            ! now we must catch the special case of the rotation with inversion.
            ! we cannot allow inversion rotations.
            ! fortunatley, and curiously, the optimal non-inverted rotation 
            ! matrix 
            ! will have the similar eigen values.
            ! we just have to make a slight change in how we handle things 
            ! depending on determinant
            
             
            	rms_ctx = sqrt(abs(ev(1)))+sqrt(abs(ev(2)))
     &            + handedness*sqrt(abs(ev(3)))
  
             
                
            ! the abs() are theoretically unneccessary since the eigen values
            ! of a real symmetric
            ! matrix are non-negative.  in practice sometimes small eigen vals 
            ! end up just negative
            temp1 = come(1)**2+come(2)**2+come(3)**2+
     &              comp(1)**2+comp(2)**2+comp(3)**2
            
            rms2_sum = (xre(jj)-xre(ii-1) + xrp(jj)-xrp(ii-1) )*invmass -temp1
          
            ! and combine the outer and cross terms into the final calculation.
            !  (the abs() just saves us a headache when the roundoff error 
            ! accidantally makes the sum negative)
            
            rms2_fast = sqrt(abs(rms2_sum-2.0D0*rms_ctx))
            
            
          
         return
         end






        subroutine prepare_cross_moments(nover,xe,xp,xm)
        ! author: charlie E. M. Strauss 2001
        ! computes a 3x3 matrix of cross moments between the x,y,z components of
        !  the two input vectors.
        ! the output is the running sum of these matricies
        implicit none
        ! inputs
        integer nover
        real*8 xe(3,nover),xp(3,nover)
        ! outputs
        
        real*8 xm(3,3,0:nover)  ! note zero origin
        
        ! local dummy variables
      
        integer i,j,k
        
        ! init array we will making into a running sum
        do j=1,3
             do k = 1,3
             
                xm(j,k,0) = 0.0d0
             enddo
           enddo
        
        
        
        do i = 1,nover
        
           do j=1,3
             do k = 1,3
             
                xm(k,j,i) = xm(k,j,i-1)+xp(j,i)*xe(k,i)   ! flopped
             enddo
           enddo
        enddo
        
        return
        end
        
        subroutine prepare_radii(n,xu,xr)
        ! this computes the running sum of the x
        ! the output is an array containing the running sum of the radii squared
        implicit none
        integer n
        real*8 xu(3,n),xr(0:n)  ! note zero origin
        
        ! local vars
        integer i

             xr(0) = 0.0d0
      
        do i=1,n
           xr(i)=xr(i-1) +xu(1,i)**2 +xu(2,i)**2 +xu(3,i)**2
        enddo
        
        return
        end
        
        
         subroutine prepare_rsums(n,xu,xs)
        ! this computes the running sum of the x
        ! the output is an array containing the running sum of the radii squared
        implicit none
        integer n
        real*8 xu(3,n),xs(3,0:n)  ! note zero origin
        
        ! local vars
        integer i
        
        
           xs(1,0)=0.0d0
           xs(2,0)=0.0d0
           xs(3,0)=0.0d0
        do i=1,n
         
           xs(1,i)=xs(1,i-1)+xu(1,i)
           xs(2,i)=xs(2,i-1)+xu(2,i)
           xs(3,i)=xs(3,i-1)+xu(3,i)
         enddo
        
        return
        end
        
        
        subroutine clear_rms()
        
        implicit none
            real*8 xre,xrp 
        real*8 xsp(3),xse(3)  ! note zero origin
        integer count
        real*8 xm(3,3)  ! note zero origin
        integer j,k
        common /rms_obj/xm,xre,xrp,xse,xsp,count
         ! init array we will making into a running sum
        do j=1,3
             do k = 1,3
             
                xm(j,k) = 0.0d0
             enddo
           enddo
         xre = 0.0d0
         xrp = 0.0d0
         
           xse(1)=0.0d0
           xse(2)=0.0d0
           xse(3)=0.0d0
           xsp(1)=0.0d0
           xsp(2)=0.0d0
           xsp(3)=0.0d0
           count=0
        return 
        end
        
        subroutine add_rms(i,xp,xe)
        implicit none
        ! computes a 3x3 matrix of cross moments between the x,y,z components of
        !  the two input vectors.
        ! the output is the running sum of these matricies
  
        ! inputs
        integer i
        real*8 xe(3,i),xp(3,i)
        
        ! outputs
        integer count
        real*8 xre,xrp 
        real*8 xsp(3),xse(3)  ! note zero origin
       
        real*8 xm(3,3)  ! note zero origin
        
        ! local dummy variables
      
        integer j,k
        
        common /rms_obj/xm,xre,xrp,xse,xsp,count
        
        count = count+1
           do j=1,3
             do k = 1,3
                xm(k,j) = xm(k,j)+xp(j,i)*xe(k,i)  ! flopped
             enddo
           enddo

           xre=xre +xe(1,i)**2 +xe(2,i)**2 +xe(3,i)**2
           xrp=xrp +xp(1,i)**2 +xp(2,i)**2 +xp(3,i)**2
      

           xse(1)=xse(1)+xe(1,i)
           xse(2)=xse(2)+xe(2,i)
           xse(3)=xse(3)+xe(3,i)
           xsp(1)=xsp(1)+xp(1,i)
           xsp(2)=xsp(2)+xp(2,i)
           xsp(3)=xsp(3)+xp(3,i)
        
        
        return
        end
        
        
       
      subroutine   rmsfitca3(npoints,xx0,xx,yy0,yy,esq)
c   AUTHOR: charlie strauss 2001
c   
c
c   This subroutine gets its alingment info via a common block!!
c   Alingment (rotation matrix) and rms(esq) are computed on the basis
c   of residues previously designated by calls to add_rms().
c   However, the rotation is applied to all Npoints of XX0,yy0 with the 
c   results returned in xx,yy.


c
c   most of this is double precision for good reasons.  
c   first there are some large differences of small numbers. 
c    second the rsymm_eignen() routine can internally have numbers
c   larger than the largest real*4 number.  (you could do some fancy foot work 
c   to rescale m_moment if you really had a problem with this.)
c   
c
c   NOTE: det is a double precision real
c   NOTE: (xx,yy) can be same arrays as (xx_0,yy_0) if desired
c   
c    
	   implicit none
       integer Npoints
       real*8 xx0(3,npoints),xx(3,npoints),yy0(3,npoints),yy(3,npoints)
        real*8 esq
        
        ! outputs
        integer count
        real*8 xre,xrp 
        real*8 xsp(3),xse(3)  ! note zero origin
       
        real*8 xm(3,3)  ! note zero origin

        
        common /rms_obj/xm,xre,xrp,xse,xsp,count
        

       ! local
       real*8 det
       integer i,j,k
       real*8 temp1,mass
       real*8 come(3),comp(3)
       real*8 ev(3)
       real*8 m_moment(3,3),rr_moment(3,3)
       real*8 det3
       real*8 rms_ctx,rms2_sum
       real*8 handedness
      
       real*8 r(3,3) ,t(3)
c compute center of mass 


        mass = count
        come(1) = xse(1)/mass       !x_com 
        come(2) = xse(2)/mass       !y_com 
        come(3) = xse(3)/mass       !z_com 
        
        comp(1) = xsp(1)/mass 
        comp(2) = xsp(2)/mass 
        comp(3) = xsp(3)/mass 

        
C	Make cross moments matrix  

       
        do  k=1,3
         do  j = 1,3
            
             m_moment(k,j) = xm(k,j)/mass- come(k)*comp(j)  ! flopped com
                     
         enddo
       enddo
        

	    det = det3(m_moment) ! get handedness  of frame from determinant
	    

            if (Dabs(det) .le. 1e-24) then
               !write(0,*) 'Warning:degenerate cross moments: det=',det
                ! might think about returning a zero rms, to avoid any chance of
                ! Floating Point Errors?
            
             esq=0.0d0
             return
            
             endif
         handedness = DSIGN(1.0D0,det)
             ! weird but documented "feature" of DSIGN(a,b) (but not SIGN) is 
             ! that if fails if a<0
             
	    !  multiply cross moments by itself
         
	    do i = 1,3
	    do j = i,3
	      rr_moment(i,j) = m_moment(1,i)*m_moment(1,j)
     &                   + m_moment(2,i)*m_moment(2,j)
     &                   + m_moment(3,i)*m_moment(3,j)
          rr_moment(j,i) = rr_moment(i,j)
     	  enddo
        enddo
        
        
            !  compute eigen values of cross-cross moments

              call rsym_eigenval(rr_moment,ev)
              
              ! reorder eigen values  so that ev(3) is the smallest eigenvalue
              
              if (ev(2) .gt. ev(3) ) then
                    if (ev(3).gt.ev(1)) then
                        temp1 = ev(3)
                        ev(3) = ev(1)
                        ev(1) = temp1
                    endif
               else
                    if (ev(2).gt.ev(1)) then
                       temp1 = ev(3)
                       ev(3) = ev(1)
                       ev(1) = temp1
                     else
                       temp1 = ev(3)
                       ev(3) = ev(2)
                       ev(2) = temp1
                 endif
               endif
                 
                 ! ev(3) is now the smallest eigen value.  the other two are not
                 !  sorted.  This is prefered order for computing the rotation matrix
                 
             call rsym_rotation(m_moment,rr_moment,ev,r)   
  
            ! now we rotate and offset all npoints 
                
            
            do i =1,npoints
              do k = 1,3   ! remove center of mass
                 yy(k,i) = yy0(k,i)-come(k)
                 xx(k,i) = xx0(k,i)-comp(k)
              enddo
              do j = 1,3   ! compute rotation
               ! temp1 = 0.0
               ! do k = 1,3 
               !  temp1 = temp1+r(j,k)*yy(k,i)
               ! enddo
               !t(j) = temp1
               t(j)=         r(j,1)*yy(1,i)
     &                      +r(j,2)*yy(2,i)
     &                      +r(j,3)*yy(3,i)

              enddo
              yy(1,i) = t(1) 
              yy(2,i) = t(2)
              yy(3,i) = t(3)
            enddo     
             
            ! now we must catch the special case of the rotation with inversion.
            ! fortunatley, and curiously, the optimal non-inverted rotation 
            ! matrix will have a similar relation between rmsd and the eigen values.
            ! we just have to make a slight change in how we handle things 
            ! depending on determinant
            
             
            rms_ctx = Dsqrt(Dabs(ev(1)))+Dsqrt(Dabs(ev(2))) 
     &                      + handedness*Dsqrt(Dabs(ev(3)))
  
            ! the abs() are theoretically unneccessary since the eigen values
            ! of a real symmetric matrix are non-negative.
            ! in practice sometimes small eigen vals end up as tiny negatives.
            ! 
            temp1 = come(1)**2+come(2)**2+come(3)**2+
     &              comp(1)**2+comp(2)**2+comp(3)**2
            
            rms2_sum = (xre + xrp)/mass -temp1
          
            ! and combine the outer and cross terms into the final calculation.
            !  (the abs() just saves us a headache when the roundoff error 
            ! accidantally makes the sum negative)
            
            esq = sqrt(abs(rms2_sum-2.0D0*rms_ctx))

         return
         end

       
        subroutine  rmsfitca2(npoints,xx,yy,ww,natsel,esq)
c AUTHOR: charlie strauss 2001
c   computes the rms between two weighted point vectors.
c   xx_0,yy_0 are the input vectors of of points and ww is their weights
c   xx,yy are by-product output vectors of the same points offset to remove center of mass
c   det is an out value of the determinant of the cross moment matrix
c   returned value is the rms
c
c   most of this is double precision for good reasons.  first there are some large
c   differences of small numbers.  and second the rsymm_eignen() routine can internally have numbers
c   larger than the largest real*4 number.  (you could do some fancy foot work to rescale things if
c   you really had a problem with this.
c
c   NOTE: det is a double precision real
c   NOTE: (xx,yy) can be same arrays as (xx_0,yy_0) if desired
c   
c    
	   implicit none
       integer Npoints,natsel
       real*8 xx(3,npoints),yy(3,npoints),ww(npoints)
  
       real*8 det
       integer i,j,k
       real*8 temp1,temp3
       real*8 ev(3)
       real*8 m_moment(3,3),rr_moment(3,3)
       real*8 det3
       real*8 rms_ctx
       real*8 rms_sum,esq
       real*8 handedness
       real*8 t(3)
       
       real*8  XPC,YPC,ZPC,XEC,YEC,ZEC,R(3,3)
       COMMON /TRANSFORM/ XPC,YPC,ZPC,XEC,YEC,ZEC,R
        
c align center of mass to origin

      CALL COMAS(xx,Ww,npoints,XPC,YPC,ZPC)
      CALL COMAS(yy,Ww,npoints,XEC,YEC,ZEC)
      
      temp3=0.0d0
      do i = 1,npoints
         temp3 = temp3+ww(i)
      enddo
        
C	Make cross moments matrix   INCLUDE THE WEIGHTS HERE
        do  k=1,3
         do  j = 1,3
            
            TEMP1=0.0
            do i = 1,npoints
               TEMP1 = TEMP1+ ww(i)*YY(k,i)*XX(j,i)
        
            enddo
            m_moment(k,j) =TEMP1    /(temp3)		! rescale by temp3
         enddo
           
       enddo
       

	    det = det3(m_moment) ! will get handedness  of frame from determinant
   
           
            if (Dabs(det) .le. 1e-24) then
                !  write(0,*) 'Warning:degenerate cross moments: det=',det
                ! might think about returning a zero rms, to avoid any chance of Floating Point Errors?
            
               esq=0.0d0
               return
            
             endif
         handedness = DSIGN(1.0D0,det)
             ! weird but documented fortran "feature" of DSIGN(a,b) (but not SIGN) is that if fails if a<0
             
	    !  multiply cross moments by itself
         
	    do i = 1,3
	    do j = i,3
	      rr_moment(i,j) = m_moment(1,i)*m_moment(1,j)
     &                   + m_moment(2,i)*m_moment(2,j)
     &                   + m_moment(3,i)*m_moment(3,j)
          rr_moment(j,i) = rr_moment(i,j)              ! well it is symmetric afterall
     	  enddo
     	  
        enddo
         
            !  compute eigen values of cross-cross moments
          
              call rsym_eigenval(rr_moment,ev)
              
              
              
               ! reorder eigen values  so that ev(3) is the smallest eigenvalue
              
              if (ev(2) .gt. ev(3) ) then
                    if (ev(3).gt.ev(1)) then
                        temp1 = ev(3)
                        ev(3) = ev(1)
                        ev(1) = temp1
                    endif
               else
                    if (ev(2).gt.ev(1)) then
                       temp1 = ev(3)
                       ev(3) = ev(1)
                       ev(1) = temp1
                     else
                       temp1 = ev(3)
                       ev(3) = ev(2)
                       ev(2) = temp1
                 endif
               endif
                              
                 ! ev(3) is now the smallest eigen value.  the other two are not
                 !  sorted.  this is prefered order for rotation matrix
               
              
            call rsym_rotation(m_moment,rr_moment,ev,r)

C$$$             do i =1,npoints
C$$$               do j = 1,3
C$$$                 temp1 = 0.0
C$$$                do k = 1,3 
C$$$                  temp1 = temp1+r(j,k)*yy(k,i)
C$$$                enddo
C$$$                t(j) = temp1
C$$$               enddo
C$$$               yy(1,i) = t(1) 
C$$$               yy(2,i) = t(2)
C$$$               yy(3,i) = t(3)
C$$$             enddo
           
            do i =1,npoints
              do j = 1,3   ! compute rotation
               t(j)=         r(j,1)*yy(1,i)
     &                      +r(j,2)*yy(2,i)
     &                      +r(j,3)*yy(3,i)
              enddo
              yy(1,i) = t(1) 
              yy(2,i) = t(2)
              yy(3,i) = t(3)
            enddo  
            ! now we must catch the special case of the rotation with inversion.
            ! we cannot allow inversion rotations.
            ! fortunatley, and curiously, the optimal non-inverted rotation matrix 
            ! will have the similar eigen values.
            ! we just have to make a slight change in how we handle things depending on determinant
            
             
            	rms_ctx = Dsqrt(Dabs(ev(1)))+Dsqrt(Dabs(ev(2))) + handedness*Dsqrt(Dabs(ev(3)))
  
                rms_ctx = rms_ctx*temp3
                
            ! the abs() are theoretically unneccessary since the eigen values of a real symmetric
            ! matrix are non-negative.  in practice sometimes small eigen vals end up just negative
            rms_sum = 0.0D0
            do  i = 1,npoints
               do  j = 1,3
               rms_sum = rms_sum+ ww(i)*(yy(j,i)**2 + xx(j,i)**2)
              enddo	
            enddo
            rms_sum = rms_sum   !   /temp3   (will use natsel instead)
             
            ! and combine the outer and cross terms into the final calculation.
            !  (the abs() just saves us a headache when the roundoff error accidantally makes the sum negative)
            
            esq = sqrt(dabs(rms_sum-2.0D0*rms_ctx)/float(natsel))
           
         return
         end
            
 
  

        
        
        	
        