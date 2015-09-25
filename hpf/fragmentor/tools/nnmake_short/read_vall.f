c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 1.12 $
c  $Date: 2002/06/18 04:41:53 $
c  $Author: rohl $

      subroutine read_vall   
      implicit none 
c----------------------------------------------------------------------
      include 'structure.h' 
       include 'path_defs.h'
      integer temps,temph,i,j,total_positions,no_align,res_zero,ss_zero
      integer fail
      integer om0,total_new, total_chains,seq_list   
      integer iunit
      character lastchain*1,lastname*4      
      real chi,junk
      character*500 filename
      character*60 buffer
c----------------------------------------------------------------------
C   modified on 4/19/2000 to input already weighted profiles 
c  (which were created w/ psiblast -C)
C    so all the sequence weighting and psuedo-count adding is out the window!
C     also it accepts the format decided on by DCChivian when the two of us remade the vall 
C	RBonneau , 4/2000
C----------------------------------------------------------------------

      total_positions=0 
      i=0 
      iunit=20

      filename=read_vall_path(1:read_vall_l)//'vall.dat.2001-02-02'
      write(0,*) 'read_vall (unit=',iunit,'): ',
     #     filename(1:index(filename,' '))
      open(iunit,file=filename,status='old',iostat=fail)
      if (fail.ne.0) then
         filename=read_vall2_path(1:read_vall2_l)//'vall.dat.2001-02-02'
         write(0,*) 'read_vall (unit=',iunit,'): ',
     #        filename(1:index(filename,' '))
         open(iunit,file=filename,status='old',iostat=fail)
         if (fail.ne.0) then
            write(0,*)"FAIL: unable to open vall file"
            stop
         endif
      endif

      do i=1,MAX_VALL
         resseq(i) = 'moom'

	 read(iunit,50,end=77)frag_name(i),chain(i),res_list(i),ss(i),
     #        resseq(i),back_seq(i),seq_list,
     #        vall_ca(1,i),vall_ca(2,i),vall_ca(3,i),phi_list(i),
     #        psi_list(i),
     #        omega_list(i),chi,no_align,acc_vall(i),junk,
     #        (vall_pro(j,i),j=0,19)
             
 50   format(a4,a1,1x,a1,1x,a1,1x,a5,2(1x,i4),3(1x,f8.2),4(1x,f8.3),
     #        1x,i3,1x,f4.2,1x,f5.3,20(1x,f5.3))        

      enddo

c if we got here then either the vall array is too small for the file or
c the vall array is exactly the right size.  so we test this by
c trying to read one more line.  if we cant read one more line then
c end test will catch this case and everybody is cool
      read(iunit,"(A60)",end=77) buffer
      write(0,*) "ERROR:Vall file bigger than declared storage, "
      write(0,*)"increace max_vall in structure.h"
      stop
      
 77   continue
      vall_size=i-1
      write(0,*)"vall_size",vall_size

CCC summary part ------------------------------------------------------
CCC not vital - output EVERY TIME stats on the vall!
CCC -------------------------------------------------------------------
CCC in case of mystery bug comment untill end!
c old      total_chains=0 
      j=0 
      om0=0 
      do i=1,vall_size 
         if (psi_list(i).eq.0.0) phi_list(i)=0.0 
         if ((omega_list(i).lt.170.0.or.omega_list(i).gt.190.0)
     #        .and.psi_list(i).ne.0.0) om0=om0+1 
         if (lastname.ne.frag_name(i).or. 
     #       lastchain.ne.chain(i)) then 
c          if (total_new.ge.1.and.i.gt.1) 
c     #        write(41,*)frag_name(i-1),chain(i-1),total_new 
        total_new=0 
            total_chains=total_chains+1 
            lastname=frag_name(i) 
            lastchain=chain(i)
         endif 
c old         if (psi_list(i).eq.0.0.and.i.lt.vall_size) then 
c old           read(resseq(i)(1:4),*)m 
c old           read(resseq(i+1)(1:4),*)n 
c old           if (m+1.ne.n) then 
c old              off_zero=off_zero+1 
c old              total_new=total_new+1 
c old           endif 
c old            phi_list(i)=0.0 
c old         endif 
         if (res_list(i).eq.' ') then 
            phi_list(i)=0.0 
            res_zero=res_zero+1 
         endif 
         if (ss(i).eq.' ') then 
            phi_list(i)=0.0 
            ss_zero=ss_zero+1 
         endif 
         if (ss(i).ne.'H'.and.ss(i).ne.'E'.and.
     #       ss(i).ne.'L') phi_list(i)=0.0 
         if (frag_name(i).eq.'    ') phi_list(i)=0.0 
         if (phi_list(i).ne.0.0) j=j+1 
      enddo 
      write(60,*)"total chains read from vall ",total_chains 
      write(60,*)"total good read from vall ",j," of ",vall_size  
      write(60,*)"omega of not 170-190 ",om0  
      write(60,*)"res zero ",res_zero," ss zero ",ss_zero 

      close(iunit) 

      temph=0    
      temps=0    
      do i=1,vall_size
         if (ss(i).eq.'H') temph=temph+1 
         if (ss(i).eq.'E') temps=temps+1 
      enddo 
      write(60,*)'vall fraction helix: ',real(temph)/vall_size,
     #           ' strand: ',real(temps)/vall_size  

  100 format(20f5.2) 
c----------------------------------
C add blosumn62 psuedocounts !!!
C no longer neccesary ( RB 4/2000) we now use psiblast checkpoint files
C-------------------------------
c      do i=1,vall_size 
c         do j=1,20 
c            temp_pro(j)=0.0 
c         enddo 
c         do j=1,20 
c            do k=1,20 
c               temp_pro(k)=temp_pro(k)
c     #            +vall_pro(j,i)*substitution(j,k)
c            enddo 
c         enddo  
c         do j=1,20 
c            vall_pro(j,i)=temp_pro(j)
c         enddo 
c      enddo 
c---------------------------------------
c  NORMALIZE OUT THE GAP         
c      do i=1,vall_size
c         total=0.0
c         do j=1,20
c            total=total+vall_pro(j,i)
c         enddo
c         if (total.eq.0.0) then
c            do j=1,20
c               vall_pro(j,i)=200.0
c            enddo
c         else
c            do j=1,20
c               vall_pro(j,i)=vall_pro(j,i)/total
c            enddo
c         endif
c         if (total.lt.0.0) then 
c            write(60,*)'total is ',total,' BIG PROBLEM!!' 
c            stop 
c         endif 
c      enddo
c      write(60,*)'profiled vallnn files' 
    
      return 
      end 

c------------------------------------------------------------------------
      subroutine read_vall_constraints_coord()

      include 'structure.h'    !vall_size, HCA,cnh,frag_name,chain,res_list
      include 'path_defs.h'    !vall_path
      include 'dipolar.h'

      real*8 vall_svdcoeff(ORDERSIZE,MAXDIPTYPES,MAX_VALL)  
                       !(data,type,residue)
                       !data: svd coefficients
                       !dipolar types for which data is precomputed:
                       !1:HN-N,2:HA-CA,3:CA-C,4:C(i-1)-N,5:C(i-1)-H   
      common /dipolar_valldata/vall_svdcoeff

      
car internal
      integer i,j,k

      character*1 res
      character*1 ch
      character*4 frag
      character*500 filename
      
      integer iunit,fail

      iunit=20
      
      filename=read_vall_path(1:read_vall_l)//
     #      'vall_cst_coord.dat.2001-02-02'
      write(0,*) 'vall_constraints file: ',filename(1:index(filename,' '))
      write(0,*) 'unit=',iunit
      open(unit=iunit,file=filename,status='old',iostat=fail)
      if (fail.ne.0) then
         filename=read_vall2_path(1:read_vall2_l)//'vallcst_coord.2001-02-02'
         write(0,*) 'read_vall (unit=',iunit,'): ',
     #        filename(1:index(filename,' '))
         open(iunit,file=filename,status='old',iostat=fail)
         if (fail.ne.0) then
            write(0,*)"FAIL: unable to open vall_cst_coord file"
            stop
         endif
      endif
      
      do i=1,vall_size
         read(iunit,1000)frag,ch,res,(cNH(j,i),j=1,3),(HCA(j,i),j=1,3),
     #        ((vall_svdcoeff(j,k,i),j=1,ORDERSIZE),k=1,MAXDIPTYPES),
     #        (xyz_vall(k,1,i),k=1,3),(xyz_vall(k,2,i),k=1,3),    !N,CA,CB,C,O
     #        (xyz_vall(k,3,i),k=1,3),(xyz_vall(k,4,i),k=1,3),
     #        (xyz_vall(k,5,i),k=1,3)

         do k=1,3                    !replace the ca coord from read_vall
            vall_ca(k,i)=xyz_vall(k,2,i)
         enddo
 
         if (frag.ne.frag_name(i) .or. ch.ne.chain(i) .or.
     #        res.ne.res_list(i)) then
            write(0,*)"ERROR:  vall_constraints file doesn't match "
     #           //"vall file at line ",i
            write(0,*)"vall file contains:",frag_name(i)//chain(i)//" "
     #           //res_list(i)
            write(0,*)"vall_constraints file contains:",
     #           frag//ch//" "//res
            stop
         endif
      enddo
      close(iunit)

      return
c 1000 format(a4,a1,1x,a1,6(1x,f8.2),25(1x,f9.6),15(1x,f8.2))
 1000 format(a4,a1,1x,a1,6(1x,f8.2),30(1x,f9.6),15(1x,f8.2))
 2000 format(a4,a1,1x,a1)
      end


