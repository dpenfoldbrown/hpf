c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 3729 $
c  $Date: 2003-11-26 15:06:37 -0500 (Wed, 26 Nov 2003) $
c  $Author: rohl $

      program make_vall_constraints

      implicit none

      include 'param.h'
      include 'structure.h'
      include 'path_defs.h'  
      include 'dipolar.h'

      character*200 filename
      integer i,j,k,l
      integer iunit,fail
      integer iargc

      integer aan(max_res)
      integer nres
      integer first,line
      real*8 angles(3,max_res)
      real*8 xyz_dp(3,5,max_res)
      real xyz(3,5,max_res)
      real HN(3)
      real HA(3)
      real svdcoeff(ORDERSIZE,MAXDIPTYPES)    !(type,coeff) 
           !type:1=HN-N,2=HA-CA,3=CA-C,4=C(i-1)-N,5=c(i-1)-HN
           !coeff:direction cosine coeff for svd

      i=iargc()
      if (i.lt.1) then
         write(0,*)'Usage: pVallCoord <vall.dat.id>'
         stop
      endif
      call getarg(1,vall_name)
      if (index(vall_name,'vall.dat').ne.1) then
         write(0,*)"non-standard vall name: ",vall_name
         write(0,*)"expecting vall.dat.<identifier>"
         stop
      endif

      vall_cst_coord_name = 'vall_cst_coord'//
     #     vall_name(5:index(vall_name,' '))
      read_vall_path='./ '
      read_vall_l =2
      read_vall2_path='./ '
      read_vall2_l=2

      call read_vall()
      iunit=12
      filename=read_vall_path(1:read_vall_l)//vall_cst_coord_name
      open(unit=iunit,file=filename,status='unknown',iostat=fail)
      if (fail.ne.0) then
         write(0,*)"unable to open for output(unit=',iunit,') :",
     #        filename(1:index(filename,' '))
         filename=read_vall2_path(1:read_vall2_l)//vall_cst_coord_name
         open(iunit,file=filename,status='old',iostat=fail)
         if (fail.ne.0) then
            write(0,*)"unable to open for output(unit=',iunit,') :",
     #           filename(1:index(filename,' '))
            write(0,*)"FAIL: unable to open vall_cst_coord file"
            stop
         endif
      endif

      write(0,*)"Output file (unit=',iunit,') :",
     #     filename(1:index(filename,' '))

      if (vall_size+1.gt.max_vall) then
         write(0,*)"increase max_vall to ",vall_size+1
      endif
      phi_list(vall_size+1)=0.0
      nres=1
      do i=1,vall_size

         if (phi_list(i+1).eq.0.0 .and.
     #        psi_list(i).eq.0.0 .and. omega_list(i).eq.0.0) then
            first=i-nres+1
            call phipsi_to_angles(phi_list(first),psi_list(first),
     #           omega_list(first),angles,nres)
            call refold_xyz(angles,xyz_dp,1,nres,nres)
            call xyz_to_position(xyz_dp,xyz,nres)
            do j=1,nres 
               call compute_hxyz(-1,xyz(1,1,j),xyz(1,2,j),
     #              xyz(1,4,j),HA)                 !HA
               call precalculate_dipolar(HA,xyz(1,2,j), !HA-CA
     #              svdcoeff(1,2))  
               call precalculate_dipolar(xyz(1,2,j),xyz(1,4,j), !CA-C
     #              svdcoeff(1,3))

               if (j.gt.1) then
                  call compute_hxyz(0,xyz(1,4,j-1),xyz(1,1,j),
     #                 xyz(1,2,j),HN) !HN
                  call precalculate_dipolar(HN,xyz(1,1,j),      !HN-N
     #                 svdcoeff(1,1))
                  call precalculate_dipolar(xyz(1,1,j),            !C(i-1)-N
     #                 xyz(1,4,j-1),svdcoeff(1,4))  
                  call precalculate_dipolar(HN,xyz(1,4,j-1),    !C(i-1)-HN
     #                 svdcoeff(1,5))
                  call precalculate_dipolar(HN,HA,    !HN(i)-HA(i)
     #                 svdcoeff(1,6))
                  
               else
                  do k=1,ORDERSIZE
                     svdcoeff(k,1)=0.0
                     svdcoeff(k,4)=0.0
                     svdcoeff(k,5)=0.0
                     svdcoeff(k,6)=0.0
                  enddo
                  do k=1,3
                     HN(k)=xyz(k,1,1)    !set HN of residue 1 to N coord
                  enddo
               endif
               line=first+j-1
               write(iunit,200)frag_name(line),chain(line),
     #              res_list(line),(HN(k),k=1,3),(HA(k),k=1,3),
     #              ((svdcoeff(k,l),k=1,ORDERSIZE),l=1,MAXDIPTYPES),
     #              (xyz(k,1,j),k=1,3),(xyz(k,2,j),k=1,3),         !N,CA,CB,C,O
     #              (xyz(k,3,j),k=1,3),(xyz(k,4,j),k=1,3),
     #              (xyz(k,5,j),k=1,3)
              


            enddo               !j
            nres=1              !reset res counter
         else
            nres=nres+1
             if (nres.gt.max_res) then
               write(0,*)"nres exceeds max_res"
               write(0,*)"nres=",nres
               stop
            endif
         endif
      enddo                      !i
      
      close(iunit)
 100  format(a4,a1,1x,a1,6(1x,f8.2),30(1x,f9.6))
 200  format(a4,a1,1x,a1,6(1x,f8.2),30(1x,f9.6),15(1x,f8.2))

      end

 
c----------------------------------------------------------------------------

      subroutine precalculate_dipolar(m,n,svdcoeff)


      implicit none

      real m(3)
      real n(3)
      real svdcoeff(*)  

car internal      
      integer i
      real temp
      real umn(3)              !unit vector between coupled atoms m,n 
      real r                   !distance between atoms m,n

      umn(1)=m(1)-n(1)          !get vector mn
      umn(2)=m(2)-n(2)
      umn(3)=m(3)-n(3)

      temp =umn(1)*umn(1)+umn(2)*umn(2)+umn(3)*umn(3)
      r=sqrt(temp)

      temp=1.0/r       
      umn(1)=umn(1)*temp        !normalize (= cosines) 
      umn(2)=umn(2)*temp
      umn(3)=umn(3)*temp

      temp =umn(1)*umn(1)                   
      svdcoeff(1)= umn(2)*umn(2)-temp  !direction cosine svd coeff
      svdcoeff(2)= umn(3)*umn(3)-temp
      svdcoeff(3)= 2.0*umn(1)*umn(2)
      svdcoeff(4)= 2.0*umn(1)*umn(3)
      svdcoeff(5)= 2.0*umn(2)*umn(3)
c      write(0,*)svdcoeff(1),svdcoeff(2),svdcoeff(3),
c     #     svdcoeff(4),svdcoeff(5)
      return
      end
      
