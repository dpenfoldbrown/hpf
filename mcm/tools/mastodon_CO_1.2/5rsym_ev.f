      subroutine rsym_rotation(mm,m,ev,rot)
                                ! Author: charlie strauss (cems) 2001
                                ! finds a (proper) rotation matrix that minimizes the rms.
                                ! this computes the rotation matrix based on the eigenvectors of m that gives the
                                ! the mimimum rms.  this is determined using mm the cross moments matrix.
                                !
                                ! for best results the third eigen value should be the 
                                ! smallest eigen value!
      
      implicit none
                                !input
      real*8 m(3,3),mm(3,3)
      real*8 ev(3)
                                !output
      real*8 rot(3,3)
                                ! local 
      integer i,j,k
      
      real*8 temp(3,3),mvec(3,3),norm
      
      
      call rsym_evector(m,ev,mvec)
      
      do i = 1,2                ! dont need no stinkin third component
         norm = 1 /sqrt(abs(ev(i) )) ! abs just fixes boo boos   
                                ! if one was nervous here one could explicitly
                                ! compute the norm of temp.
         do j = 1,3
            temp(j,i) = 0.0d0
            do k = 1,3
               temp(j,i) = temp(j,i)+mvec(k,i)*mm(j,k)
            enddo
            temp(j,i) = temp(j,i)*norm
         enddo
      enddo

      temp(1,3) =  temp(2,1)*temp(3,2) -temp(2,2)*temp(3,1)
      temp(2,3) = -temp(1,1)*temp(3,2) +temp(1,2)*temp(3,1)
      temp(3,3) =  temp(1,1)*temp(2,2) -temp(1,2)*temp(2,1)

      do i = 1,3
         do j = 1,3
            rot(j,i) = 0.0d0
            do k = 1,3
               rot(j,i) = rot(j,i)+ temp(i,k)*mvec(j,k)
            enddo
         enddo
      enddo
      
      return
      end
      
               real*8 function  det3(m)
         implicit none
         real*8 m(3,3)
c AUTHOR: charlie strauss 2001
! determinant of a 3x3 matrix
! cems: cute factoid: det of a 3x3 is the dot product of one row with the cross product of the other two.
! this explains why a right hand coordinate system has a positive determinant. cute huh?
          det3 =    m(1,3)*(  m(2,1)*m(3,2) -m(2,2)*m(3,1)) 
     &               -m(2,3)*( m(1,1)*m(3,2) -m(1,2)*m(3,1))
     &               +m(3,3)*( m(1,1)*m(2,2) -m(1,2)*m(2,1))
           return
           end
           
  

        subroutine  rsym_eigenval(m,ev)
        implicit none
c computes the eigen values of a real symmetric 3x3 matrix
c AUTHOR: charlie strauss 2001
c the method used is a deterministic analytic result that I hand factored.(whew!)
c Amusingly, while I suspect this factorization is not yet optimal in the number of calcs required
c I cannot find a signifcantly better one.  
c  (if it were optimal I suspect I would not have to compute imaginary numbers that
c  I know must eventually cancel out to give a net real result.)
c this method relys on the fact that an analytic factoring of an order 3 polynomial exists.
c m(3,3) is a 3x3 real symmetric matrix: only the upper triangle is actually used
c  ev(3) is a real vector of eigen values, not neccesarily in sorted order.


        real*8 m(3,3)
        real*8 xx,yy,zz,xy,xz,yz
        real*8 a,b,c,s0
        complex*16 f1,f2,f3,f4,f5  ! should be complex*16 but old versions of absoft have a BUG!
        real*8 ev(3)
        real*8 norm
        ! first, for lexical sanity, name some temporary variables
        xx= m(1,1)
        yy= m(2,2)    
        zz=m(3,3)     
        xy=m(1,2)   
        xz=m(1,3)  
        yz=m(2,3)   
        

        ! coefficients of characterisitic polynomial
        a = xx+yy+zz
        b = -xx*zz-xx*yy-yy*zz+xy*xy+xz*xz+yz*yz
        c=xx*yy*zz-xz*xz*yy-xy*xy*zz-yz*yz*xx+2*xy*xz*yz
        
        ! For numerical dynamic range we rescale the variables here
        ! with rare exceptions this is unnessary but it doesn't add much to the calcluation time.
        ! this also allows use of real*4 if desired. cems
        !  for complex*8 the rescaling trigger  should be  1e15  (or less)
        !  for complex*16 the rescaling trigger should be  1e150 (or less)
        ! note we completely ignore the possiblity of needing rescaling to avoid
        ! underflow due to too small numbers. left as an excercise to the reader.

        norm = max(abs(a),max(abs(b),abs(c)))
        if (norm.gt.1.0D50) then   ! rescaling trigger
          a = a/norm
          b = b/(norm**2)
          c = c/(norm**3)
         else
           norm=1.0d0
         endif
        ! we undo the scaling by de-scaling the eigen values at the end
        
       
        ! eigenvals are the roots of the characteristic polymonial  0 = c+b*e +a*e^2  - e^3
        ! solving for the three roots now:
        ! dont try to follow this in detail: its just a tricky
        ! factorization of the formulas for cubic equation roots.

         s0 = -12.0D0*(b**3)-3.0D0*(b**2)*(a**2)+54.0D0*c*b*a+81.0D0*(c**2)+12.0D0*c*(a**3) ! butt ugly term

          f1 = b*a/6.0D0+c/2.0D0+(a**3)/27.0+cdsqrt(s0*(1.0D0,0.0D0))/18.0D0

          f2 =   f1**(1.0/3.0)   ! note f1 is a complex number
                                 ! footnote: Absoft's OSX compiler has a bug on complex*16
                                 ! and computes the incorrect value at this step.  
                                 ! should be a double in exponent?

          f3 = (-b/3.0D0-(a**2)/9.0D0)/f2

         f4 =  f2-f3
        
         f5 = Dsqrt(3.0D0)*(0.0D0,1.0D0)*(f2+f3)   ! just our imaginary friend, mr. i
        
 	    s0 = a/3.0D0
        ev(1) =  f4    ! note implicitly take real part, imag part "should" be zero
        ev(1) = norm*(ev(1)+s0 ) ! do addition after type conversion in previous line.
        ev(2) = -f4+f5  ! note implicitly take real part, imag part is zero
        ev(2) = norm*(ev(2)*0.5D0 +s0 )   
        ev(3) = (-f4-f5) ! note implicitly take real part, imag part is zero
        ev(3) = norm*(ev(3)*0.5D0 +s0)
      
        return
        end



      
      subroutine rsym_evector(m,ev,mvec)
                                ! Author: charlie strauss (cems) 2001
                                ! given original matrix plus its eigen values
                                ! compute the eigen vectors.
                                ! USAGE notice: for computing min rms rotations be sure to call this 
                                ! with the lowest eigen value in Ev(3).
                                !
                                ! The minimal factorization of the eigenvector problem I have derived below has a puzzling
                                ! or, rather, interesting assymetry.  Namely, it doesn't matter what 
                                !  either ZZ=M(3,3) is or what Ev(3) is!
                                !  One reason for this is that of course all we need to  know is contained in the
                                !  M matrix in the firstplace, so the eigen values overdetermine the problem
                                !  We are just exploiting the redundant info in the eigen values to hasten the solution.
                                !  What I dont know is if infact there exists any symmetic form using the eigen values.
                                ! 
                                ! we deliberately introduce another assymetry for numerical stability, and
                                ! to force proper rotations  (since eigen vectors are not unique within a sign change). 
                                !   first, we explicitly numerically norm the vectors to 1 rather than assuming the
                                !   algebra will do it with enough accuracy. ( and its faster to boot!)  
                                !   second, we explicitly compute the third vector as being the cross product of the
                                !   previous two rather than using the M matrix and third eigen value.
                                !   If you arrange the eigen values so that the third eigen value is the lowest
                                !   then this gaurentees a stable result under the case of either very small
                                !   eigen value or degenerate eigen values.
                                !   this norm, and ignoring the third eigen value also gaurentee us the even if the
                                !   eigen vectors are not perfectly accurate that, at least, the matrix
                                !   that results is a pure orthonormal rotation matrix, which for most applications
                                !   is the most important form of accuracy.
                                !
      
      implicit none
                                !input
      real*8 m(3,3)
      real*8 ev(3)
                                !output
      real*8 mvec(3,3)
                                !local
      real*8 xx,yy,xy,zx,yz     !zz
      real*8 e1,e2,e3,znorm
      integer i
      
      
                                ! first, for sanity only, name some temporary variables
      !zz= m(3,3)     
      xy= m(1,2)   
      zx= m(1,3)  
      yz= m(2,3)   

      if (ev(1).ne.ev(2)) then  ! test for degenerate eigen values

         do i = 1,2             ! only computer first two eigen vectors using this method
                                ! note you could compute all three this way if you wanted to,
                                ! but you would run into problems with degenerate eigen values.
            
            xx=m(1,1)-ev(i)
            yy=m(2,2)-ev(i)
                                ! I marvel at how simple this is when you know the eigen values.
            e1 = xy*yz-zx*yy
            e2 = xy*zx-yz*xx
            e3 = xx*yy-xy*xy
            
            znorm= sqrt(e1**2+e2**2+e3**2) 
            
            mvec(1,i) = e1/znorm
            mvec(2,i) = e2/znorm
            mvec(3,i) = e3/znorm
            
         enddo
         
                                ! now compute the third eigenvector
         mvec(1,3) =  mvec(2,1)*mvec(3,2) -mvec(2,2)*mvec(3,1)
         mvec(2,3) = -mvec(1,1)*mvec(3,2) +mvec(1,2)*mvec(3,1)
         mvec(3,3) =  mvec(1,1)*mvec(2,2) -mvec(1,2)*mvec(2,1)
         
                                ! pathologically nervous people would explicitly normalize this vector too.
         
         return
         
      else
         
         if (ev(2).ne.ev(3)) then
            write(0,*) " hey is this the right thing to be doing??? "
            
            do i = 2,3          ! Okay, since 1 and 2 are degenerate we will use 2 and 3 instead.
               
               xx=m(1,1)-ev(i)
               yy=m(2,2)-ev(i)
                                ! I marvel at how simple this is when you know the eigen values.
               e1 = xy*yz-zx*yy
               e2 = xy*zx-yz*xx
               e3 = xx*yy-xy*xy ! yes you sharp eyed person, its not quite symmetric here too. 
                                !                   life is odd.
               
               znorm= sqrt(e1**2+e2**2+e3**2) 
               
               mvec(1,i) = e1/znorm
               mvec(2,i) = e2/znorm
               mvec(3,i) = e3/znorm
               
            enddo
            
                                ! now compute the third eigenvector
            mvec(1,1) =  mvec(2,2)*mvec(3,3) -mvec(2,3)*mvec(3,2)
            mvec(2,1) = -mvec(1,2)*mvec(3,3) +mvec(1,3)*mvec(3,2)
            mvec(3,1) =  mvec(1,2)*mvec(2,3) -mvec(1,3)*mvec(2,2)
            
                                ! pathologically nervous people would explicitly normalize this vector too.
            
            return
            
         else
            
            write(0,*) 'warning: all eigen values are equal'
            
            do i = 1,3
               mvec(1,i) = 0.0d0
               mvec(2,i) = 0.0d0
               mvec(3,i) = 0.0d0
               mvec(i,i) = 1.0d0
            enddo
            return
         endif
      endif
      
      
      end
      
      
      
      
      
