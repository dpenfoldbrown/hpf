c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 3729 $
c  $Date: 2003-11-26 15:06:37 -0500 (Wed, 26 Nov 2003) $
c  $Author: rohl $

car   5/2/00: made consistent with rosetta conventions:

car     renamed position->xyz to avoid confusion
car       'position'                         'xyz'
*****  Atom numbers                     
**      1  2  3 4  5             **      1  2  3 4  5 
**      NH CA CB C CO .....      **      NH CA C CO CB.....

car  added accessory subroutines to handle format convresions
car  position_to_xyz,xyz_to_position,phipsi_to_angles,angles_to_phipsi 

**   name changes of subroutines:
car    refold->refold_xyz
car    Dcross->Dcros           
car    transform-> transform_xyz  

car  appended functions from vectors.f:
car   Dmxm, Dtranslate,Drotate,DtransposeM,Ddotprod

car   uses functions from bystroff.f:
car     Dunitvec,Dcros,Dsubvec 

car   removed subroutines 
car       transform_sc (transforms with sidechain coordinates)
car       xyz_to_torsion (converts xyz to angles)
c-----------------------------------------------------------------------
      subroutine position_to_xyz(position,xyz,nres)

      implicit none

car   changes ordering of coordinates in xyz array
car    1        1  2  3  4  5             -->      1  2  3  4  5 
car             NH CA CB C  CO            -->      NH CA C  CO CB

      integer nres
      real position(3,5,nres)
      real*8 xyz(3,5,nres)
     
      integer i,j

      do i=1,nres
         do j=1,3
            xyz(j,1,i)=position(j,1,i)
            xyz(j,2,i)=position(j,2,i)
            xyz(j,3,i)=position(j,4,i)
            xyz(j,4,i)=position(j,5,i)
            xyz(j,5,i)=position(j,3,i)
         enddo
      enddo
      
      return
      end
c---------------------------------------------------------------------
      subroutine xyz_to_position(xyz,position,nres)

      implicit none

car   changes ordering of coordinates in xyz array
car    1        1  2  3  4  5             -->      1  2  3  4  5 
car             NH CA C  CO CB            -->      NH CA CB C  CO

      integer nres
      real position(3,5,nres)
      real*8 xyz(3,5,nres)
      

      integer i,j

      do i=1,nres
         do j=1,3
            position(j,1,i)=xyz(j,1,i)
            position(j,2,i)=xyz(j,2,i)
            position(j,3,i)=xyz(j,5,i)
            position(j,4,i)=xyz(j,3,i)
            position(j,5,i)=xyz(j,4,i)
         enddo
      enddo
      
      return
      end
c------------------------------------------------------------------------

      subroutine phipsi_to_angles(phi,psi,omega,angles,nres)

car converts rosetta-style phi,psi,omega arrays in degrees to double
car  precision angles array in radians for refold_xyz

      implicit none

      integer nres
      real phi(nres)
      real psi(nres)
      real omega(nres)
      real*8 angles(3,nres)

      integer i

      real convert
      parameter (convert=3.14159/180.0)
      
      do i=1,nres
         angles(1,i)=phi(i)*convert
         angles(2,i)=psi(i)*convert
         angles(3,i)=omega(i)*convert
      enddo

      return
      end
c--------------------------------------------------------------------------

      subroutine angles_to_phipsi(angles,phi,psi,omega,nres)

car converts rosetta-style phi,psi,omega arrays in degrees to double
car  precision angles array in radians for refold_xyz

      implicit none

      integer nres
      real*8 angles(3,nres)
      real phi(nres)
      real psi(nres)
      real omega(nres)

      integer i

      real convert
      parameter (convert=180./3.14159)
      
      do i=1,nres
         phi(i)=angles(1,i)*convert
         psi(i)=angles(2,i)*convert
         omega(i)=angles(3,i)*convert
      enddo

      return
      end
c--------------------------------------------------------------------------
     
      subroutine refold_xyz(angles, xyz, first, last, total)

car   This routine generates a complete coordinate set for backbone atoms
car   using torsion angles from residues from first to last.
car   Residues outside of this range are assumed to already have coordinates
car   consistent with the backbone torsion angles, and the cordinates
car   of these residues are transformed to be continguous with the 
car   region inside first-last
car   angles are in radians
car   all variables passed in and out via argument lists, no common blocks


      implicit none

car parameter
      integer maxatom
      parameter (maxatom=5)

      real*8 angles(3,*)             !(phi/psi/omega, residue) 
      real*8 xyz(3,maxatom,*)   !(xyz, N/CA/C/O/CB residue)  cartesian coord
      integer first
      integer last
      integer total

car internal
      real*8 axes1(3,3)      !rotation matrix for first atom-
                           !eventually the coord system to place next atom
      real*8 axes2(3,3)      !rotation matrix for last atom
      real*8 N1(3)           !coordinates of N of first
      real*8 N2(3)           !coordinates of N of last+1
      real*8 init_angles(3)  !fake angles for building first residue
      integer start
      integer i

      if (first.eq.1) then  !build first residue
         do i=1,3
            N1(i) = 0.0             !Initial Nitrogen at Origin
            init_angles(i) = 0.0    !omega of "residue" 0,phi,psi of res 1
            axes1(i,1)=0.0
            axes1(i,2)=0.0
            axes1(i,3)=0.0
            axes1(i,i)=1.0          ! Initial axes = Identity matrix
         enddo
         init_angles(3)=angles(2,1)       !psi of res 1
         call torsion_to_xyz(init_angles, xyz(1,1,1), axes1, N1,1) 
         start = 2
      else
         do i=1,3 
            N1(i) = xyz(i,1,first)          !initial N xyz
            axes1(i,1)=xyz(i,1,first)- xyz(i,3,first-1)   !C-N vector
            axes1(i,2)= xyz(i,3,first-1)-xyz(i,2,first-1) !CA-C vector
         enddo
         call norm_axes(axes1)                   !Convert to axes
         start=first
      endif

      if (last.lt.total) then    !rotation matrix before torsion_to_xyz
         do i=1,3 
            N2(i) = xyz(i,1,last+1)          !trailing N xyz
            axes2(i,1)=xyz(i,1,last+1)-xyz(i,3,last) !C-N vector
            axes2(i,2)=xyz(i,3,last)-xyz(i,2,last)   !CA-C vector       
         enddo
         call norm_axes(axes2)            !Convert to axes
      endif      

      call torsion_to_xyz(angles(3,start-1), xyz(1,1,start), axes1, 
     #     N1,last-start+1)    !Refold the residues from "first" to "last"
      
 
car now need to transform coordinates of either 1->last or last->total
car so the two pieces are contiguous
      if (last.eq.total) return             !no transformation required

c      write(0,*)N1(1),N1(2),N1(3)
c      do j=1,3
c         write(0,*)axes1(1,j),axes1(2,j),axes1(3,j)
c      enddo
c      write(0,*)
c      write(0,*)N2(1),N2(2),N2(3)
c      do j=1,3
c         write(0,*)axes2(1,j),axes2(2,j),axes2(3,j)
c      enddo
c      write(0,*)

      if ((total-last).le.last) then        !transform  last->total
         call transform_xyz(xyz(1,1,last+1),maxatom*(total-last),
     #   axes1,axes2,N1,N2)
      else
         call transform_xyz(xyz,maxatom*(last),axes2,axes1,N2,N1)
      endif
      return
      end


c------------------------------------------------------------------------

      subroutine torsion_to_xyz(angles,xyz,axes,Nxyz,nres)

c     Generates coordinates for nres residues from the backbone torsion
c     angles.  
c
c     It assumes the axes have already 
c     been initialized properly, and requires the initial nitrogen position
c     Nxyz, to start the folding process.
c
c     Upon leaving the routine, the axes
c     will be properly initialized for the following residue, and the Nxyz
c     position will be calculated for the following residue.
c     
c     N position depends on psi(i-1)
c     CA determined by omega(i-1)
c     C determined by phi
c
c     axes are unit vectors of coordinate system to place next atom
c     x axis is along previous bond (ie dihedral) 

c     note: angles are shifted relative to normal order in angle array!!
c     this is because angles are used in order omega/phi/psi to build 
c     CA, C and N

      implicit none

car parameters
      integer maxatom
      parameter (maxatom=5)

car input/output      
      real*8 angles(3,*)             !(omega/phi/psi, residue) 
      real*8 xyz(3,maxatom,*)  !(xyz, N/CA/C/O/CB,residue)  cartesian coord
      real*8 axes(3,3)            !rotation matrix for next atom
      real*8 Nxyz(3)            !coord of N of first residue
      integer nres             !num of residues to refold

car internal
      integer i,j
      real*8 CBpos(3)           !xyz, measured from CA, w/C axes
      real*8 Opos(3)            !xyz, measured from C', w/ N+1 axes
      data CBpos /-0.532605943, 0.775186438, 1.19513831 /
      data Opos  /-0.670458753, 1.03241588,   0.00161216482/

c Keith Frost's values
c      data CBpos /-0.532605938, 0.775186444, -1.19513831/
c      data Opos  /-0.670458748, 1.03241588,   0.0/
  
      
      do i=1,3                      
         xyz(i,1,1)=Nxyz(i)      !Place initial N
      enddo

      do i=1,nres
         call put_atom(2, xyz(1,1,i), 
     #        axes,angles(1,i),xyz(1,2,i))         !CA
         call put_atom(3, xyz(1,2,i), 
     #        axes,angles(2,i),xyz(1,3,i))          !C'
         do j=1,3               !CB
            xyz(j,5,i) = xyz(j,2,i) +
     #           CBpos(1) * axes(j,1) + CBpos(2) * axes(j,2) + 
     #           CBpos(3)*axes(j,3)
         enddo
         if (i.eq.nres) then
            call put_atom(1, xyz(1,3,i),axes, angles(3,i),Nxyz)
         else
            call put_atom(1, xyz(1,3,i),axes, !N(i+1)      
     #           angles(3,i),xyz(1,1,i+1))
         endif
         do j=1,3                                              !O (i)
            xyz(j,4,i) = xyz(j,3,i) +   
     #           Opos(1) * axes(j,1) + Opos(2) * axes(j,2)
c     #              +  Opos(3) * axes(i,3) 
         enddo
      enddo
      
 
      return
      end
c----------------------------------------------------------------------------
      subroutine put_atom(atomType, pos1, axes, ang, pos2)

car   given the coord of the previous backbone atom (pos1), the coordinates
car   of the next backbone atom of type (atomType) are returned (pos2)
car   using the rotation matrix (axes) and dihedral angle (angle). The
car   rotation matrix is replaced with the one required for the next atom.

      implicit none
      
      integer atomType   !N=1, CA=2, C'=3
      real*8 pos1(3)       !coord of previous backbone atom
      real*8 axes(3,3)     !rotation matrix for next atom
      real*8 ang          !relevant dihedral angle (omega, phi, or psi)
      real*8 pos2(3)       !coord of next backbone atom
      
car internal 
      integer i
      real*8 bkbZ(3)      ! N, CA, C'      
      real*8 bkbR(3)      ! N, CA, C'     

car axes: 
car x axis is along last bond made
car y axis is such that the 2nd to last bond is in the xy plane
car ie for CA, x axis is from c' to n with previous ca in xy plane

c      data bkbZ /0.5866521260492, 0.76607355607, 0.551474428047/
c      data bkbR /1.1922300926227, 1.24037947222, 1.421801230135/
c Keith's values
      data bkbZ /0.58665212, 0.76607356, 0.55084803/ ! N, CA, C'
      data bkbR /1.19223010, 1.24037948, 1.42017031/ ! N, CA, C'   

car spherical to rectangular coordinates
car x=Rsin(ang1)cos(theta)  where theta is in xy plane measured from x to y  
car y=Rsin(ang1)sin(theta)  ang1 is measured from z toward projection in xy
car z=Rcos(ang1)
car Rcos(ang1) -length of projection on z axis
car Rsin(ang1) -length of projection on xy plane

car here axes are rotated so that z->x, x->y, y->z
car  dihedral is measured from -y to z
car so theta= -dihedral +180  and
car    sin(theta)= -sin(dihedral)
car    cos(theta)= -cos(dihedral)

car bkbZ=Rcos(ang) -length of projection on x axis
car bkbR=Rsin(ang) -length of projection on yz plane
car so x=bkbZ y=bkbR*(-cos(dihedral)) z=bkbR*(-sin(dihedral))
car multiply by axes to transform to proper coordinate system

!      tmp=3.14159*(ang)/180.0 !convert to radians
      do i=1,3
        pos2(i) = pos1(i)+bkbZ(atomType)*axes(i,1)
     #         - dcos(ang) * bkbR(atomType) * axes(i,2)
     #         - dsin(ang) * bkbR(atomType) * axes(i,3)
      enddo
     
car update axes for next atom placement
C      write(0,*)"**",pos1(1),pos2(1)
      do i=1,3
         axes(i,2)=axes(i,1)
         axes(i,1)=pos2(i)-pos1(i)
      enddo
      call norm_axes(axes)
     
      return
      end
  

c-------------------------------------------------------------------------
      subroutine norm_axes(axes)
      
car this subroutine returns a matrix of unit vectors for a coordinate
car system with the x axis along a1 and the xy plane defined by the
car vectors a1 and a2
      implicit none
      real*8 axes(3,3)

car hardcode cross products and unit vectors
      call Dcros(axes(1,1),axes(1,2),axes(1,3))    !fxn in bystroff.f
      call Dcros(axes(1,3),axes(1,1),axes(1,2))
      call Dunitvec(axes(1,1),axes(1,1))           !fxn in bystroff.f
      call Dunitvec(axes(1,3),axes(1,3))
      call Dunitvec(axes(1,2),axes(1,2))

      return
      end
c ---------------------------------------------------------------------------

      subroutine transform_xyz(xyz,natom, axes1,axes2,pt1,pt2) 

car transforms coordinates of natoms in the  xyz array 
car by rotation required to align axes2 on axes1, and the subsequent
car translation required to position pt2 on pt1
      
      implicit none

      real*8 axes1(3,3)  !rotation matrix of final coord system
      real*8 axes2(3,3)  !rotation matix of starting coord system 
      real*8 xyz(3,*) !coordinates to be transformed
      integer natom     !number of atoms to be transformed
      real*8 pt1(3),pt2(3)     !coordinates first atom to be moved to

car internal
      real*8 rot(3,3)    !rotation matrix to rotate axes2 into axes1
      real*8 offset(3)   !offset between atom
      
      integer i

      do i=1,3
         offset(i)=-pt2(i)
      enddo
      call DtransposeM(axes1)
      call Dmxm(axes1, axes2, rot) !Compute rotation matrix from axes
      call Drotate(offset,rot,1)
      call Dtranslate(offset,pt1,1)

      call Drotate(xyz,rot,natom)
      call Dtranslate(xyz,offset,natom)
      
      return
      end

c--------------------------------------------------------------------------
      subroutine DtransposeM(mat)  !F90 has intrinsic function 'transpose'

c      0 1 2        0 3 6
c      3 4 5  ===>  1 4 7 
c      6 7 8        2 5 8

      implicit none
      
      real*8 mat(3,3)
      
      real*8 tmp
      
      tmp=mat(1,2)
      mat(1,2)=mat(2,1)
      mat(2,1)=tmp

      tmp=mat(1,3)
      mat(1,3)=mat(3,1)
      mat(3,1)=tmp

      tmp=mat(2,3)
      mat(2,3)=mat(3,2)
      mat(3,2)=tmp

      return
      end
            
c---------------------------------------------------------------------
      subroutine Dmxm(A, B, AxB)

car multiplies two 3x3 matrices
      
c    0 1 2     a b c     (012).(adg) (012).(beh) (012).(cfi)
c    3 4 5  *  d e f  =  (345).(adg) (345).(beh) (012).(cfi)
c    6 7 8     g h i     (678).(adg) (678).(beh) (678).(cfi)

      implicit none
      
      real*8 A(3,3)
      real*8 B(3,3)
      real*8 AxB(3,3)
      
      integer i
      real*8 col(3)        !temporary holder for column vector

car function calls
      real*8 Ddotprod
      
      do i=1,3
         col(1)=B(i,1)
         col(2)=B(i,2)
         col(3)=B(i,3)
         AxB(i,1)=Ddotprod(A(1,1),col)
         AxB(i,2)=Ddotprod(A(1,2),col)
         AxB(i,3)=Ddotprod(A(1,3),col)
      enddo
      
      return
      end
c------------------------------------------------------------------------
      real*8 function Ddotprod(v1,v2)

      implicit none

      real*8 v1(3),v2(3)

      Ddotprod = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)

      return
      end

c------------------------------------------------------------------------
      subroutine Dtranslate(xyz, vec, num)

car translates num atoms in xyz by vec

      implicit none

      real*8 xyz(3,*)        !coordinates
      real*8 vec(3)
      integer num          !total atoms to translate

      integer i

      
      do i=1,num
         xyz(1,i)=xyz(1,i)+vec(1)
         xyz(2,i)=xyz(2,i)+vec(2)
         xyz(3,i)=xyz(3,i)+vec(3)
      enddo
      
      return
      end

c--------------------------------------------------------------------------
      subroutine Drotate(xyz, mat, num)

car rotates num atoms in xyz by mat

      implicit none

      real*8 xyz(3,*)
      real*8 mat(3,3)
      integer num

      integer i
      real*8 x,y,z

      do i=1,num
         x=mat(1,1)*xyz(1,i)+mat(2,1)*xyz(2,i)+mat(3,1)*xyz(3,i)
         y=mat(1,2)*xyz(1,i)+mat(2,2)*xyz(2,i)+mat(3,2)*xyz(3,i)
         z=mat(1,3)*xyz(1,i)+mat(2,3)*xyz(2,i)+mat(3,3)*xyz(3,i)
        
         xyz(1,i)=x
         xyz(2,i)=y
         xyz(3,i)=z
      enddo

      return
      end


      subroutine Dcros(v1,v2,v3) ! warning there is another sub called cross
** v3 = v1 X v2
      real*8 v1(3), v2(3), v3(3)
      v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
      v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
      v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
      return
      end
      subroutine Dunitvec(v1,v2)
      real*8 v1(3),v2(3),cc
       real*8 invX 
      cc = sqrt(v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3)  )
 
      if (cc.eq.0.) cc=0.0001
! cems changed division to inverse multiply
      invX = 1.0/cc
      v2(1) = v1(1)*invX
      v2(2) = v1(2)*invX
      v2(3) = v1(3)*invX
    
      return
      end


	 
