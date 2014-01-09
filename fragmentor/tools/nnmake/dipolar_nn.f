c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 17461 $
c  $Date: 2007-10-01 16:16:26 -0400 (Mon, 01 Oct 2007) $
c  $Author: dekim $


      logical function get_dipolar_exist()

      implicit none
      include 'dipolar.h'

      get_dipolar_exist=dipolar_exist

      return
      end

c------------------------------------------------------------------------
      subroutine dipolar_set_verbose(flag)

car use to turn out detailed output for next call to eval_dipolar
car reset to false after every call to eval_dipolar

      implicit none

      logical flag
      logical verbose
      common /dipolar_vflag/ verbose
      data verbose /.false./
      verbose=flag
      return
      end
c-------------------------------------------------------------------------
      logical function dipolar_get_verbose()

      implicit none

      logical verbose
      common /dipolar_vflag/ verbose
      dipolar_get_verbose=verbose

      return
      end
c------------------------------------------------------------------------
      subroutine dipolar_set_verbose_unit(iunit)

      implicit none
      include 'dipolar.h'

      integer iunit
      dunit=iunit

      return
      end
c------------------------------------------------------------------------

      subroutine read_dipolar(filename,iunit,exist,dipolar_count)

car reads in dipolar constraint file and stores data in arrays

car note: all coupling constants assumed to be collected under
car the same conditions, such that the same order tensor can be applied
car to all couplings

car file format:
car line1: version     (not yet in use)
car line2: parameters  (not yet in use)
car line3: comment
car line4: # of lines to read   (<=physical number of lines)
car line5 -> beginning in 1st column
car    tag      res1     atom1    res2    atom2      J
car     a1  2x   i4   1x  a4   ix  a4  1x  4x    1x f10.2
car
car 200  format(a1,2x,i4,1x,a4,1x,i4,1x,a4,1x,f10.2)
car
car if the first character of a line is '#', the line is read, but
car    the data is ignored.
car otherwise, the first column of each line contains a tag identifying
car    the set to which the measurement belongs. (ie different experimental
car    conditions with different alignment tensor) tags may be any
car    single character other than '#'
car
car Residues and atoms are stored internally with res1 <= res2,
car   regardless of input order, in addition, the list is sorted in
car   increasing value of res1

car Only five types of dipolar couplings are assumed to be possible
car       1H-1H
car       1H-15N
car       1H-13C
car      13C-15N
car      13C-13C
car If more are added, need to correct function get_invDcnst accordingly
car Also, only bb atoms (including H) are included at present.

      implicit none

      include 'dipolar.h'

car input
      character*(*) filename    !filename with dipolar constraints
      integer iunit      !safe unit number of dipolar constraint file

car output
      logical exist
      integer dipolar_count

*****  Atom numbers start
*****   with the backbone atoms as follows
**      1  2  3 4  5
**      NH CA CB C CO .....

car   0  HN         !all protons should be <= 0
car  -1  HA  or HA1
car  -2  HA2        !glycine

car parameters
      integer maxprint
      parameter (maxprint=50)

car internal
      integer error      !iostatus of constraint file
      character*100 comment
      integer res1,res2
      character*4 atom1,atom2   ! n,ca,cb, c, o,hn,ha
      character*1 flag
      character*1 setmap(MAXDIPSETS)  !maps characters to set numbers

      integer nlines
      integer i,j
      integer type

      real*8 x(ORDERSIZE,MAXDIPSETS)
      real*8 Azz(MAXDIPSETS)
      real*8 eta(MAXDIPSETS)
      real*8  paf_rot(3,3,MAXDIPSETS)
      logical reject(MAXDIPSETS)

      common /align_tensor/ x,Azz,eta,paf_rot,reject


car function calls

      integer atomlookup
      real*8 get_invDcnst
      integer dipolartype

car data: fixed bond lengths in Angstroms
      real fixed_dist(5)   !fixed pair distances for different types
      data fixed_dist /1.01, 1.08,  1.52325877, 1.32874878, 2.032764/
car H-X bond lengths from get_hxyz
car bb bond lengths from refold or template residue
car Clore et al. (1998)JMR 133:216 uses
car  HCA 1.08; HN 1.02;  N-C' 1.341; CA-C' 1.525

      data dunit /0/

      dipolar_count = 0
      open(iunit,file=filename,status='old',iostat=error)
      if (error.ne.0) then
         write(0,*)"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
         write(0,*)"WARNING: DIPOLAR CONSTRAINT FILE NOT FOUND"
         write(0,*)" Searched for: ", filename(1:index(filename,' '))
         write(0,*)" Dipolar constraints will not be used "
         write(0,*)"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
         ndipolar=0
         ndipolar_sets=0
         exist=.false.
         return
      endif

      exist = .true.
      dipolar_exist=.true.
      write(0,*) " dipolar constraints: ",filename(1:maxprint)
      read(iunit,100,err=1000)comment
      write(0,160)comment(1:maxprint)
      read(iunit,100,err=1000)comment   !parameter line, not in use
      write(0,160)comment(1:maxprint)
      read(iunit,100,err=1000)comment   !user comment line
      write(0,160)comment(1:maxprint)
      read(iunit,150,err=2000)nlines

      ndipolar_sets=1
      ndipolar = 1
      do  i=1,nlines
         read(iunit,200,end=500,err=5000)flag,res1,atom1,res2,atom2,Jdipolar(ndipolar)
         if (flag.ne.'#') then
            if (res1.le.res2) then
               pairDipolar(1,ndipolar)=res1
               pairDipolar(3,ndipolar)=res2
               pairDipolar(2,ndipolar)=atomlookup(atom1)
               pairDipolar(4,ndipolar)=atomlookup(atom2)
            else
               pairDipolar(1,ndipolar)=res2
               pairDipolar(3,ndipolar)=res1
               pairDipolar(2,ndipolar)=atomlookup(atom2)
               pairDipolar(4,ndipolar)=atomlookup(atom1)
            endif
            if (ndipolar.eq.1) setmap(1)=flag
            do j=1,nDipolar_sets
               if (flag.eq.setmap(j)) then
                  pairDipolar(5,ndipolar)=j
                  goto 400
               endif
            enddo
            ndipolar_sets=ndipolar_sets+1
            setmap(ndipolar_sets)=flag
            pairDipolar(5,ndipolar)=ndipolar_sets
 400        continue            !jump here if set exists
            ndipolar = ndipolar +1
         endif
c         write(0,*)i,(pairdipolar(j,i),j=1,5),Jdipolar(i)
       enddo
 500  continue
      close(iunit)
      ndipolar = ndipolar -1
      dipolar_count = ndipolar

      if (ndipolar.gt.MAXDIPOLAR)  then
         write(0,*) "ERROR in read_dipolar:Ndipolar exceeds MAXDIPOLAR."
         write(0,*) "Ndipolar = ",ndipolar," MAXDIPOLAR = ", MAXDIPOLAR
         write(0,*) "Increase MAXDIPOLAR in dipolar.h"
         write(0,*)"NNMAKE FAILED: .dpl file has more than ",
     #        MAXDIPOLAR," restraints"
         stop
      endif

      call sortdipolar(ndipolar,pairDipolar,Jdipolar)

car store invDcnst values & distances for bonded atoms
car conventions:  if distance can't be precalculated, Ddist<0
car  and the value of invDmax must be calculated in eval_dipolar

      do i=1,ndipolar
         invDcnst(i)=get_invDcnst(pairdipolar(1,i))
         type=dipolartype(pairdipolar(1,i))
         if (type.eq.0) then
            write(0,*)"WARNING: unrecognized dipolar type in dipolar "
     #           //"restraint file"
            write(0,*)"residue 1: ",pairdipolar(1,i)," atom 1: ",
     #           pairdipolar(2,i)," residue 2: ",pairdipolar(3,i),
     #           " atom 2: ",pairdipolar(4,i)," data set: ",
     #           pairdipolar(4,i)
            write(0,*)"This restraint will not be used"
         elseif (type.le.5) then
            Ddist(i)=fixed_dist(type)
            invDmax(i)=invDcnst(i)*(fixed_dist(type)**3)
            Jdipolar(i)=Jdipolar(i)*invDmax(i)
         else
            Ddist(i)=0.0
            invDmax(i)=invDcnst(i)
            Jdipolar(i)=Jdipolar(i)*invDcnst(i)
          endif
      enddo

c$$$      do i=1,ndipolar
c$$$         write(0,*)i,(pairdipolar(j,i),j=1,5),Jdipolar(i)/invDmax(i)
c$$$      enddo

      write(0,*)"number of dipolar constraints:", ndipolar
      write(0,*)"number of dipolar data sets:", ndipolar_sets

car initialize values of tensor and reject flags:

      do j=1,MAXDIPSETS
         do i=1,ORDERSIZE
            x(i,j)=0.0
         enddo
         reject(j)=.true.
      enddo

      return

 100  format(a)
 150  format(i5)
 160  format(4x,a)
 200  format(a1,2x,i4,1x,a4,1x,i4,1x,a4,1x,f10.2)
 1000 write(0,*)"NNMAKE FAILED: format error in .dpl file lines 1-3"
 2000 write(0,*)"NNMAKE FAILED: format error in .dpl file line 4"
 5000 write(0,*)"NNMAKE FAILED: format error in .dpl file restraint ",i
      stop
      end
c-------------------------------------------------------------------------
      integer function atomlookup(atom)

car differs for rosetta scheme
      implicit none

      character*4 atom

      if (atom(1:1).ne.' ' .and. atom(1:1).ne.'1' .and. atom(1:1).ne.'2'
     #    .and. atom(1:1) .ne. '3' .and. atom(1:1) .ne. '#') then
         atom(2:4)=atom(1:3)      !check that fields are aligned
         atom(1:1)=' '
      endif

      if (atom.eq.' H ') then
         atomlookup = 0
         return
      elseif (atom.eq.' HA ' .or. atom.eq.'1HA ') then
         atomlookup = -1
         return
      elseif (atom.eq.'2HA ') then
         atomlookup = -2
         return
      elseif  (atom.eq.' CA ' .or. atom.eq.'CA  ') then
         atomlookup = 2
         return
      elseif  (atom.eq.' N   ' .or. atom.eq.'N    ') then
         atomlookup = 1
         return
      elseif (atom.eq.' C  '.or.atom.eq.'C ') then !C'=atom 4
         atomlookup = 4
         return
      else
         write(0,*) "Atom type: '",atom,"' unknown in dipolar_nn.f"
         stop
      endif
      return
      end

c--------------------------------------------------------------------------
      character*4 function namelookup(num)

      implicit none

      integer num

      if (num .eq. 0) then
         namelookup = ' H  '
      elseif (num.eq.-1) then
         namelookup=' HA '
      elseif (num.eq.-2) then
         namelookup='2HA '
      elseif (num.eq.2) then
         namelookup=' CA '
      elseif (num.eq.1) then
         namelookup=' N  '
      elseif (num.eq.4) then
         namelookup=' C  '
      else
         namelookup='UNKN'
      endif
      return
      end

c--------------------------------------------------------------------------
      integer function dipolartype(pairdipolar)

car function differs for rosetta numbering scheme

      implicit none

car input
      integer pairdipolar(*)

car internal
      integer low,high

      dipolartype=0       !undefined
      if (pairdipolar(1).eq.pairdipolar(3)) then       !same residue
         low=min(pairdipolar(2),pairdipolar(4))
         high=max(pairdipolar(2),pairdipolar(4))
         if (low.eq.0 .and. high.eq.1) dipolartype=1   !HN-N
         if (low.eq.-1 .and. high.eq.2) dipolartype=2  !HA-CA
         if (low.eq.2 .and. high.eq.4) dipolartype=3    !CA-C'
         if (low.eq.-1 .and. high.eq.0) dipolartype=6  !HN-HA
      elseif (pairdipolar(1).eq.pairdipolar(3)-1) then !adjacent residue
         low=pairdipolar(4)
         high=pairdipolar(2)       !C' must be of i-1
         if (low.eq.1 .and. high.eq.4) dipolartype=4 !C'(i-1)-N(i)
         if (low.eq.0 .and. high.eq.4) dipolartype=5 !C'(i-1)-HN(i)
      endif
      if (dipolartype.gt.0) return

      write(0,*) "undefined dipolar type"
      write(0,*) "residue1 ",pairdipolar(1)," atom1 ",pairdipolar(2),
     #     " residue2 ",pairdipolar(3)," atom2 ",pairdipolar(4),
     #     " data set ",pairdipolar(5)

      return
      end
c---------------------------------------------------------------------------
      real*8 function get_invDcnst(pairdipolar)

car Dmax =-gamma(n)*gamma(m)*h/(2*pi**2*r**3         Losonczi et al.
car      =-mu0*gamma(n)*gamma(m)*h/(16*pi**3*r**3)   Fischer et al.
car where gamma is gyromagnetic ratio of nucleus
car  gamma(1H)     2.6752e8 1/(Ts)   T=Ns/(Cm)
car  gamma(15N)   -2.712e7 1/(Ts)
car  gamma(13C)    6.728e7 1/(Ts)
car  h             6.62608e-34 Js    J=Nm
car  uo            4*pi*1e-7 Ns**2/C**2
car
car from Fischer eq, units come out in Hz (1/s)
car Dmax=-(1e-7)(6.62608e-34)(gamma(m))(gamma(n)/[4(pi**2)(1e-30)(r**3)]
car where 1e-30 converts r3 from cubic meters to cubic angstroms
car invDmax=invDcnst*(r**3) with r in Angstroms
car reduced dipolar coupling = Dobs*invDmax=Jdipolar*r**3*invDcnst
car with r**3 in Angstroms

car dipolar types           Dmax   These are values from ORDERTEN
car         1H-1H        -240200.0        from Prestegard lab
car         1H-15N         24350.0  !note paper has a factor of 2 error
car         1H-13C        -60400.0    !Losonczi et al. JMR 1999 138:334
car         15N-13C         6125.0    !eq 1, 2 in denom should be squared
car         13C-13C       -15200.0

car my calculations
car HH -120118.4     difference 2x
car HN   12177.1
car HC  -30209.2
car NC    3062.48
car CC   -7597.47

car and inverses:
car  HH -0.0000083251
car  HN   0.0000821215
car  HC  -0.0000331025
car  NC    0.000326533
car  CC   -0.000131623

      implicit none
car input
      integer pairDipolar(*)

car internal
      integer low,high      !atom numbers, assumed ordered
                            !note assumptions in possible pairs

car note:  all atoms number >2 are C
car            atom 1 is N
car            atoms <1 are H
car same for rosetta and non-rosetta scheme

      low=min(pairdipolar(2),pairdipolar(4))
      high=max(pairdipolar(2),pairdipolar(4))

      if (high.le.0) then         !H-H
         get_invDcnst= -0.0000083251
         return
      endif
      if (low.le.0) then          !H-X
         if (high.eq.1) then         !HN
            get_invDcnst=   0.0000821215
            return
         else                        !HC  (HX=HC if X.ne.N or H)
            get_invDcnst=  -0.0000331025
            return
         endif
      elseif (low.eq.1) then      !N-C  (NX=NC if X.ne.H)
         get_invDcnst=  0.000326533
         return
      endif
      get_invDcnst=   -0.000131623     !C-C  (X1-X2=CC if X1.ne.H or N)
      return
      end


c------------------------------------------------------------------------
      subroutine calc_orderparam(x,vec,Azz,eta)

      implicit none

      real*8 x(*)
      real*8 vec(3,3)           !eigenvectors
      real*8 Azz
      real*8 eta
car internal
      integer i
      real*8 S(3,3)      !order matrix
      real*8 val(3)      !eigenvalues
      integer nrot
      integer sort(3)   ! sorted index to val: val(sort(1))=largest abs val.
      real*8 temp1,temp2



c assemble ordermatrix
      S(1,1)=-x(1)-x(2)
      S(2,2)=x(1)
      S(3,3)=x(2)
      S(1,2)=x(3)
      S(2,1)=x(3)
      S(1,3)=x(4)
      S(3,1)=x(4)
      S(2,3)=x(5)
      S(3,2)=x(5)


c diagonalize
      call dp_jacobi_xyz(S,val,vec,nrot)

c sort eigenvalues

      sort(1)=1
      sort(2)=2
      sort(3)=3
      if (abs(val(1)).lt.abs(val(2))) then     !largest absolute value
         sort(2)=1
         sort(1)=2
      endif
      if( abs(val(sort(2))).lt.abs(val(3)) ) then
         sort(3)=sort(2)
         sort(2)=3
         if (abs(val(sort(1))).lt.abs(val(3))) then
            sort(2) = sort(1)
            sort(1) = 3
         endif
      endif


      Azz=val(sort(1))
      eta=(2.0/3.0)*abs(val(sort(2))-val(sort(3))/Azz)

c sort eigen values      !largest to smallest : Azz,Ayy,Axx
      temp1 = val(sort(1))
      temp2 = val(sort(2))
      val(3) = val(sort(3))
      val(2) = temp2
      val(1) = temp1

c sort eigen vectors
      do i = 1,3
         temp1 = vec(i,sort(3))
         temp2 = vec(i,sort(2))
         vec(i,3)=vec(i,sort(1))
         vec(i,1) =temp1
         vec(i,2) =temp2
      enddo

      Azz=val(1)
      eta=(2.0/3.0)*(val(3)-val(2))/val(1)

      return
      end

c-------------------------------------------------------------------------
      subroutine calc_dipscore(A,x,b,map,ndata,Azz,score)

      implicit none
      include 'dipolar.h'

      real*8 A(MAXDIPOLAR,ORDERSIZE)
      real*8 x(ORDERSIZE)
      real*8 b(MAXDIPOLAR)
      integer map(MAXDIPOLAR)
      integer ndata
      real*8 Azz
      real score

      logical  verbose

      integer i,j
      real Jcalc

      logical dipolar_get_verbose
      character*4 namelookup

      verbose=dipolar_get_verbose()

      score=0.0
      do i=1,ndata
         Jcalc=0.
         do j=1,ORDERSIZE
            Jcalc=Jcalc+A(i,j)*x(j)
         enddo
         score=score + (b(i)-Jcalc)**2
         if (verbose) write(dunit,1000)pairdipolar(1,map(i)),
     #        b(i)/invDmax(map(i)),Jcalc/invDmax(map(i)),
     #        (b(i)/invDmax(map(i))-Jcalc/invDmax(map(i)))**2,
     #        pairdipolar(1,map(i)),namelookup(pairdipolar(2,map(i))),
     #        pairdipolar(3,map(i)),namelookup(pairdipolar(4,map(i)))
      enddo
      score=score/(ndata*(Azz**2))  !normalize by #data, principal order param.

      return

 1000 format(i3,3(1x,f7.2),1x,i3,1x,a4,1x,i3,1x,a4)

      end


c---------------------------------------------------------------------------
      subroutine get_atomxyz(xyz,res,atom,atomxyz)

car differs for rosetta scheme
car returns coordinates for a specified atom, either copied from the
car position array, or for protons, calculated from the backbone
car coordinates using a template coordinate

      implicit none

      real xyz(3,5,*)     !coordinates
      integer res         !residue (ie third index of xyz)
      integer atom
      real atomxyz(3)

*****  Atom numbers start
*****   with the backbone atoms as follows
**      1  2  3 4  5
**      NH CA CB C CO  .....

car   0  HN     !all protons should be <= 0
car  -1  HA


      if (atom.ge.1) then       !backbone atom
         atomxyz(1)=xyz(1,atom,res)
         atomxyz(2)=xyz(2,atom,res)
         atomxyz(3)=xyz(3,atom,res)
         return
      elseif (atom.eq.0) then
         if (res.eq.1) then
            write(0,*)'WARNING: Error in eval_dipolar!'
            write(0,*)' Attempt to calculate coupling constant for '
     #           //'HN of N-terminal residue.'
            stop
         endif
         call compute_hxyz(atom,xyz(1,4,res-1),xyz(1,1,res),
     #        xyz(1,2,res),atomxyz)     !C' is atom 4 in rosetta scheme
         return
      elseif (atom.eq.-1 .or. atom.eq.-2)then     !HA,HA1,HA2
         call compute_hxyz(atom,xyz(1,1,res),xyz(1,2,res),
     #        xyz(1,4,res),atomxyz)      !C' is atom 4 in rosetta scheme
         return
      endif

      write(0,*)'WARNING: unrecognized atom type in dipolar restraints'
      write(0,*)'atom number', atom, ' is not defined in get_atomxyz'

      return
      end


c-------------------------------------------------------------------------
      subroutine compute_hxyz(type,atom1,atom2,atom3,hxyz)

car differs for rosetta conventions cross->cros

car returns hxyz, the coordinates of a proton of specified type.
car types:
car   0  HN         !all protons should be <= 0
car  -1  HA  or HA1
car  -2  HA2        !glycine
car The coordinates of three atoms are required to place the proton
car atom2=atom to which proton is bound
car atom1=N-terminally adjacent to atom2, (or closer to backbone)
car atom3-C-terminally adjacent to atom2, (or further from backbone)

      implicit none

car input
      integer type    !proton type 0=HN,-1=HA,-2=glyHA2
      real atom1(3)   !coordinates of atom N-term to proton-bound atom
      real atom2(3)   !coordinates of proton-bound atom
      real atom3(3)   !coordinates of atom C-term to proton-bound atom
car output
      real hxyz(3)

car internal
      integer i
      real axes(3,3)

car parameters
      real*8 Hpos(3,-2:0)  !xyz measured from atom2 using coordinate sys where
                           !atom2->atom1 is colinear w/ x axis
                           !atom3 lies in xy plane and
                           !atom3->atom2 has a +y component

      data Hpos /-0.230853993,  0.488168418,  0.935306426, ! -2=HA2
     #           -0.230853993,  0.488168418, -0.935306426, ! -1=HA or HA1
     #           -0.477136907,  0.890191189,  0.0  /       !  0=HN

car  -0.481861031,  0.899004954.  0.0     !HN when bondlength = 1.02 clore
car  -0.477136907,  0.890191189,  0.0     !HN when bondlength = 1.01 cems
car  -0.472412784,  0.881377423,  0.0     !HN when bondlength = 1.00 bk
car  -0.491966549  0.882082136  0.0 !1.01A from N, bisecting the CNCA angle
                                          !w/refold template, 1.01 A

      do i=1,3
         axes(i,1)=atom3(i)-atom2(i)
         axes(i,2)=atom2(i)-atom1(i)
      enddo
      call cros(axes(1,1),axes(1,2),axes(1,3))   !generate rotation matrix
      call cros(axes(1,3),axes(1,1),axes(1,2))  !rosetta function cros
      call unitvec(axes(1,1),axes(1,1))
      call unitvec(axes(1,3),axes(1,3))
      call unitvec(axes(1,2),axes(1,2))

      do i=1,3          !template coord x rotation matrix + base coord
         hxyz(i)=atom2(i)+ Hpos(1,type) * axes(i,1) +
     #        Hpos(2,type) * axes(i,2) +
     #        Hpos(3,type) * axes(i,3)
      enddo

      return
      end

c--------------------------------------------------------------------------
      subroutine calc_ordermatrix(nrow,A,b,x,reject)

c returns the order matrix x and if invalid, sets reject true
c nrow gets increased if the matrix needs to be padded
c no other values are altered

      implicit none
      include 'dipolar.h'

car input/output
      integer nrow
      real*8 A(MAXDIPOLAR,ORDERSIZE)
      real*8 b(MAXDIPOLAR)
      real*8 x(ORDERSIZE)
      logical reject

car parameters
      real factor
      parameter (factor=1e-6)   !cutoff factor for singular values in svd

car internal
      integer i,j
      real*8 U(MAXDIPOLAR,ORDERSIZE)
      real*8 w(ORDERSIZE)              !singular values
      real*8 v(ORDERSIZE,ORDERSIZE)
      real*8 wmin, wmax           !min and max values of w
      integer sing              !number of singular values in w
      real*8 Sxx
      integer ndata

      ndata=nrow
      do i=1,nrow                   !copy A
         U(i,1)= A(i,1)   !copy into U
         U(i,2)= A(i,2)
         U(i,3)= A(i,3)
         U(i,4)= A(i,4)
         U(i,5)= A(i,5)
      enddo

      do while (nrow.lt.ORDERSIZE)   !not enoughs rows, fill to square w/zeros
         nrow=nrow+1
         U(nrow,1)=0.0
         U(nrow,2)=0.0
         U(nrow,3)=0.0
         U(nrow,4)=0.0
         U(nrow,5)=0.0
      enddo
      call svdcmp(U,nrow,ORDERSIZE,MAXDIPOLAR,ORDERSIZE,w,v)

car check w for singular values prior to inversion
      wmax=0.0
      do j=1,ORDERSIZE
         if (w(j).gt.wmax) wmax=w(j)
      enddo
      wmin=wmax*factor
      sing=0
      do j=1,ORDERSIZE
         if (w(j).lt.wmin) then
            w(j)=0.0
            sing=sing+1
         endif
      enddo
      if (sing .gt. abs(ndata-ORDERSIZE))
     #     write(0,*)"SVD yielded a matrix singular above expectation "
     #     //"in get_ordermatrix"


car find solution for exact dipolar values:
      call svbksb(u,w,v,nrow,ORDERSIZE,MAXDIPOLAR,ORDERSIZE,b,x)

car x components: (Syy,Szz,Sxy,Sxz,Syz)
car check for acceptable values

      reject=.false.
      Sxx=-x(1)-x(2)
      if (Sxx.lt.-0.5 .or. Sxx.gt.1.0) reject=.true. !Sxx
      if (x(1).lt.-0.5 .or. x(1).gt.1.0) reject=.true. !Syy
      if (x(2).lt.-0.5 .or. x(2).gt.1.0) reject=.true. !Szz
      if (x(3).lt.-0.75 .or. x(3).gt.0.75) reject=.true. !Sxy
      if (x(4).lt.-0.75 .or. x(4).gt.0.75) reject=.true. !Sxz
      if (x(5).lt.-0.75 .or. x(5).gt.0.75) reject=.true. !Syz

      if (reject) then
         write(0,*)'order matrix not physically meaningful'
         return
                                !try with errors on dipolar values? map error?
                                !score=0.0
                                !return
      endif

      return
      end


c---------------------------------------------------------------------------
c  subroutines from numerical recipes
c---------------------------------------------------------------------------

      SUBROUTINE svdcmp(a,m,n,mp,np,w,v)

car Given a matrix A, with logical dimensions MxN and physical dimensions
car MPxNP, this routine computes its singular value decomposition, A=U.W.Vt
car The matrix U replaces A on output.  The diagonal matrix of singular
car values (W) is output as a vector W, The matrix V (not the transpose Vt)
car is output as V.  M must be greater than or equal to N.  If it is smaller,
car then A should be filled up to square with zero rows.

      INTEGER m,mp,n,np,NMAX
      REAL*8 a(mp,np),v(np,np),w(np)
      PARAMETER (NMAX=500)
CU    USES pythag
      INTEGER i,its,j,jj,k,l,nm
      REAL*8 anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX),pythag
      g=0.0
      scale=0.0
      anorm=0.0
      do 25 i=1,n
        l=i+1
        rv1(i)=scale*g
        g=0.0
        s=0.0
        scale=0.0
        if(i.le.m)then
          do 11 k=i,m
            scale=scale+abs(a(k,i))
11        continue
          if(scale.ne.0.0)then
            do 12 k=i,m
              a(k,i)=a(k,i)/scale
              s=s+a(k,i)*a(k,i)
12          continue
            f=a(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,i)=f-g
            do 15 j=l,n
              s=0.0
              do 13 k=i,m
                s=s+a(k,i)*a(k,j)
13            continue
              f=s/h
              do 14 k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
14            continue
15          continue
            do 16 k=i,m
              a(k,i)=scale*a(k,i)
16          continue
          endif
        endif
        w(i)=scale *g
        g=0.0
        s=0.0
        scale=0.0
        if((i.le.m).and.(i.ne.n))then
          do 17 k=l,n
            scale=scale+abs(a(i,k))
17        continue
          if(scale.ne.0.0)then
            do 18 k=l,n
              a(i,k)=a(i,k)/scale
              s=s+a(i,k)*a(i,k)
18          continue
            f=a(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,l)=f-g
            do 19 k=l,n
              rv1(k)=a(i,k)/h
19          continue
            do 23 j=l,m
              s=0.0
              do 21 k=l,n
                s=s+a(j,k)*a(i,k)
21            continue
              do 22 k=l,n
                a(j,k)=a(j,k)+s*rv1(k)
22            continue
23          continue
            do 24 k=l,n
              a(i,k)=scale*a(i,k)
24          continue
          endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
25    continue
      do 32 i=n,1,-1
        if(i.lt.n)then
          if(g.ne.0.0)then
            do 26 j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
26          continue
            do 29 j=l,n
              s=0.0
              do 27 k=l,n
                s=s+a(i,k)*v(k,j)
27            continue
              do 28 k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
28            continue
29          continue
          endif
          do 31 j=l,n
            v(i,j)=0.0
            v(j,i)=0.0
31        continue
        endif
        v(i,i)=1.0
        g=rv1(i)
        l=i
32    continue
      do 39 i=min(m,n),1,-1
        l=i+1
        g=w(i)
        do 33 j=l,n
          a(i,j)=0.0
33      continue
        if(g.ne.0.0)then
          g=1.0/g
          do 36 j=l,n
            s=0.0
            do 34 k=l,m
              s=s+a(k,i)*a(k,j)
34          continue
            f=(s/a(i,i))*g
            do 35 k=i,m
              a(k,j)=a(k,j)+f*a(k,i)
35          continue
36        continue
          do 37 j=i,m
            a(j,i)=a(j,i)*g
37        continue
        else
          do 38 j= i,m
            a(j,i)=0.0
38        continue
        endif
        a(i,i)=a(i,i)+1.0
39    continue
      do 49 k=n,1,-1
        do 48 its=1,30
          do 41 l=k,1,-1
            nm=l-1
            if((abs(rv1(l))+anorm).eq.anorm)  goto 2
            if((abs(w(nm))+anorm).eq.anorm)  goto 1
41        continue
1         c=0.0
          s=1.0
          do 43 i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            if((abs(f)+anorm).eq.anorm) goto 2
            g=w(i)
            h=pythag(f,g)
            w(i)=h
            h=1.0/h
            c= (g*h)
            s=-(f*h)
            do 42 j=1,m
              y=a(j,nm)
              z=a(j,i)
              a(j,nm)=(y*c)+(z*s)
              a(j,i)=-(y*s)+(z*c)
42          continue
43        continue
2         z=w(k)
          if(l.eq.k)then
            if(z.lt.0.0)then
              w(k)=-z
              do 44 j=1,n
                v(j,k)=-v(j,k)
44            continue
            endif
            goto 3
          endif
          if(its.eq.30) pause 'no convergence in svdcmp'
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
          g=pythag(f,1.0d0)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0
          s=1.0
          do 47 j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=pythag(f,h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do 45 jj=1,n
              x=v(jj,j)
              z=v(jj,i)
              v(jj,j)= (x*c)+(z*s)
              v(jj,i)=-(x*s)+(z*c)
45          continue
            z=pythag(f,h)
            w(j)=z
            if(z.ne.0.0)then
              z=1.0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do 46 jj=1,m
              y=a(jj,j)
              z=a(jj,i)
              a(jj,j)= (y*c)+(z*s)
              a(jj,i)=-(y*s)+(z*c)
46          continue
47        continue
          rv1(l)=0.0
          rv1(k)=f
          w(k)=x
48      continue
3       continue
49    continue
      return
      END
c---------------------------------------------------------------------

      FUNCTION pythag(a,b)
      REAL*8 a,b,pythag
      REAL*8 absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pythag=absa*sqrt(1.+(absb/absa)**2)
      else
        if(absb.eq.0.)then
          pythag=0.
        else
          pythag=absb*sqrt(1.+(absa/absb)**2)
        endif
      endif
      return
      END

c---------------------------------------------------------------------
      SUBROUTINE svbksb(u,w,v,m,n,mp,np,b,x)

car solves A.x=b for a vector X where A is specified by the arrays U,W,V
car as returned by SVDCMP. M and N are the logical dimensions of A. MP
car and NP are the phpyscial dimensions of A.  B is the input right-hand
car side.  X is the output solution vector.  No input quantities are
car destroyed.


      INTEGER m,mp,n,np,NMAX
      REAL*8 b(mp),u(mp,np),v(np,np),w(np),x(np)
      PARAMETER (NMAX=500)
      INTEGER i,j,jj
      REAL*8 s,tmp(NMAX)
      do 12 j=1,n
        s=0.
        if(w(j).ne.0.)then
          do 11 i=1,m
            s=s+u(i,j)*b(i)
11        continue
          s=s/w(j)
        endif
        tmp(j)=s
12    continue
      do 14 j=1,n
        s=0.
        do 13 jj=1,n
          s=s+v(j,jj)*tmp(jj)
13      continue
        x(j)=s
14    continue
      return
      END
c--------------------------------------------------------------------------
      SUBROUTINE dp_jacobi_xyz(a,d,v,nrot)
c------------------------------------------------------------------------
car from numerical recipes, modified by R.Bonneau to handle 3x3 matrices
car specifically
car computes all eigenvalues and eigenvectors of a symmetric matrix a (3x3).
car A is destroyed. D returns the eigenvalues,V contains the normalized
car eigenvectors, by column. NROT returns the number of jacobi rotations
car required.
car single precision

      INTEGER n,np,nrot,NMAX
      REAL*8 a(3,3),d(3),v(3,3)
      PARAMETER (NMAX=500)
      INTEGER i,ip,iq,j
      REAL*8 c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
c-------------------------------------------------------------------------
c solves for eigen-values of real-symetric matrix
c using jacobi method of iterative rotations
c fast for a three by three
c not neccesarily optimal- numerical recipies
c in fortran : Cambridge-University Press
c--------------------------------------------------------------------------
      n=3   !!! hardwired
      np=3  !!! hardwired
      do 12 ip=1,n   !!! coordinates
        do 11 iq=1,n
          v(ip,iq)=0.
11      continue
        v(ip,ip)=1.
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.
13    continue
      nrot=0
      do 24 i=1,50
        sm=0.
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))
14        continue
15      continue
        if(sm.eq.0.)return
        if(i.lt.4)then
          tresh=0.2*sm/n**2
        else
          tresh=0.
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+
     # g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h
              else
                theta=0.5*h/a(ip,iq)
                t=1./(abs(theta)+sqrt(1.+theta**2))
                if(theta.lt.0.)t=-t
              endif
              c=1./sqrt(1+t**2)
              s=t*c
              tau=s/(1.+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.
23      continue
24    continue
      pause 'too many iterations in jacobi'
      return
      END
c------------------------------------------------------------------------
      subroutine sortdipolar(n,ra,rb)

car heapsort algorithm from numerical recipes
car modified here from hpsort to use multidimensional ra and also sort rb
      implicit none
      integer size                !first dim of ra
      parameter (size=5)
      integer key                 !key to sort on
      parameter (key=1)

      INTEGER n
      integer ra(size,n)
      real*8 rb(n)
      INTEGER i,ir,j,l,k
      integer rra(size)
      real*8 rrb

      if (n.lt.2) return
      l=n/2+1
      ir=n
10    continue
      if(l.gt.1)then
         l=l-1
         do k=1,size
            rra(k)=ra(k,l)
         enddo
         rrb=rb(l)
      else
         do k=1,size
            rra(k)=ra(k,ir)
         enddo
         rrb=rb(ir)
         do k=1,size
            ra(k,ir)=ra(k,1)
         enddo
         rb(ir)=rb(1)
         ir=ir-1
         if(ir.eq.1)then
            do k=1,size
               ra(k,1)=rra(k)
            enddo
            rb(1)=rrb
            return
         endif
      endif
      i=l
      j=l+l
 20   if(j.le.ir)then
         if(j.lt.ir)then
      if(ra(key,j).lt.ra(key,j+1)) j=j+1
         endif
         if(rra(key).lt.ra(key,j)) then
            do k=1,size
               ra(k,i)=ra(k,j)
            enddo
            rb(i)=rb(j)
            i=j
            j=j+j
         else
            j=ir+1
         endif
         goto 20
      endif
      do k=1,size
         ra(k,i)=rra(k)
      enddo
      rb(i)=rrb
      goto 10
      END

c---------------------------------------------------------------------------
c---------------------------------------------------------------------------
c  subroutines from bystroff.f required by dipolar.f
c--------------------------------------------------------------------------
c-------------------------------------------------------------------------
      subroutine unitvec(v1,v2)
      real*4 v1(3),v2(3),cc
       real invX
      cc = sqrt(v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3)  )

      if (cc.eq.0.)cc=0.0001
! cems changed division to inverse multiply
      invX = 1.0/cc
      v2(1) = v1(1)*invX
      v2(2) = v1(2)*invX
      v2(3) = v1(3)*invX

      return
      end


      subroutine cros(v1,v2,v3) ! warning there is another sub called cross
** v3 = v1 X v2
      real*4 v1(3), v2(3), v3(3)
      v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
      v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
      v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
      return
      end

cc------------------------------------------------------------------------
cc additional subroutines required for nnmake (not rosetta)
cc------------------------------------------------------------------------

      real function retrieve_frag_dipolar_score(frag_start)

      implicit none
      include 'dipolar.h'

      integer max_frags                  !see also prepar_dipolar_score
      parameter (max_frags=600000)

      real frag_dipolar_score(max_frags)
      common /dipolar_frag_scores/ frag_dipolar_score

      integer frag_start

      if (.not.dipolar_exist) then
	retrieve_frag_dipolar_score=0.0
	return
      endif
      retrieve_frag_dipolar_score= frag_dipolar_score(frag_start)

      return
      end


c------------------------------------------------------------------------
      real function score_frag_dipolar(fragstart)

car calculates the total chi squared between predicted and calculated
car dipolar coupling constants:
car score = sum(Jdipolar_obs**2-Jdipolar_calc**2)/((# constraints)
car normalized by the principal order parameter
car
car Note that dipolar
car couplings involving the HN of residue i or j are not included in the
car score because the orientation of the HN bond vector is not defined
car unless the phi of residue first is defined.
car
car Input dipolar coupling data are obtained via the common block in
car dipolar.h
car
car References:
car  Losonczi et al. (1999) JMR 138:334
car  Fischer et al. (1999) Biochemistry 38:9013
car  Clore et al., (1998) JMR 133:216
car  Bothner-By in Encyclopedia of NMR (Grant & Harris ed) pp2932 Wiley 1995

      implicit none

      include 'dipolar.h'

car input
      integer fragstart         !starting vall residue of fragment

      logical verbose           !print out diagnostics
      common /dipolar_vflag/verbose
car input (via common block)
      integer map(MAXDIPOLAR,MAXDIPSETS)
      integer maplength(MAXDIPSETS)
      integer mapstart
      common /dipolar_map/ map,maplength,mapstart

car internal
      real score              !temporary holder for score
      integer iset
      integer nrow             !rows of matrix A = number of dipolars in frag
      real setscore

      real*8 A(MAXDIPOLAR,ORDERSIZE)  !coeff matrix for SVD; see  Losonzi et al
      real*8 b(MAXDIPOLAR)      !reduced dipolar couplings (Jdipolar/Jmax)

      real*8 x(ORDERSIZE,MAXDIPSETS)
      real*8 Azz(MAXDIPSETS)
      real*8 eta(MAXDIPSETS)
      real*8  paf_rot(3,3,MAXDIPSETS)
      logical reject(MAXDIPSETS)

      common /align_tensor/ x,Azz,eta,paf_rot,reject

      integer i
car ----------------function body ----------------------------------------

car constraints assumed to be stored by read_dipolar such that
car pairDipolar(1,i) .le.  pairDipolar(3,i)    and
car pairDipolar(1,i) .le.  pairDipolar(1,i+1)
      score=0.0
      score_frag_dipolar=0.0
      if (.not.dipolar_exist) return
      if (verbose) write(0,*)"residue, Jobs, Jcalc, residual**2"
      do iset=1,ndipolar_sets
         if (verbose) write(0,*)"dipolar constraint set ",iset
         nrow=maplength(iset)
         if (nrow.lt. 5) goto 300
         call assemble_datamatrix_nn(fragstart,mapstart,map(1,iset),
     #        nrow,A,b)
         call calc_ordermatrix(nrow,A,b,x(1,iset),reject(iset))
         if (reject(iset)) then
            write(0,*)"iset=",iset," fragstart=",fragstart
            write(0,*)'rejected'
            write(0,*)x(1,iset),x(2,iset),x(3,iset),x(4,iset),x(5,iset)
            do i=1,maplength(iset)
               write(0,*)i,A(i,1),A(i,2),A(i,3),A(i,4),A(i,5),
     #              b(i)
            enddo
            goto 300
         endif
         call calc_orderparam(x(1,iset),paf_rot(1,1,iset),Azz(iset),
     #        eta(iset))
         call calc_dipscore(A,x(1,iset),b,map(1,iset),maplength(iset),
     #        Azz(iset),setscore)
         if (verbose) write(0,*)"Dipolar score for constraint set ",
     #        iset," : ",setscore
         score=score+setscore
 300     continue             !next set
      enddo
      if (verbose) write(0,*)"total dipolar score: ",score
      score_frag_dipolar=score
      return
      end
c------------------------------------------------------------------------
      subroutine prepare_dipolar_score(first,size,vall_size)

car selects the appropriate dipolar coupling data to be scored against,
car and generates a map to the data  stored in a common block accessed
car by assemble_datamatrix_nn

      implicit none

      include 'dipolar.h'

      integer max_frags                  !see also retrieve_frag_dipolar_score
      parameter (max_frags=600000)

car input
      integer first      !first protein residue to score, ie frag begin
      integer size
      integer vall_size

car output (to common block)
      integer maplength(MAXDIPSETS)
      integer map(MAXDIPOLAR,MAXDIPSETS)
      integer mapstart     !first query residue  for which map is built
      common /dipolar_map/ map,maplength,mapstart

      real frag_dipolar_score(max_frags)
      common /dipolar_frag_scores/ frag_dipolar_score

car internal

      integer iset,i
      integer last
      integer nfrag

      real mean,sd,weight

car function calls
      integer dipolartype
      real score_frag_dipolar

      if (.not.dipolar_exist) return

car check array size
      if (max_frags.lt.vall_size) then
         write(0,*)"increase max_frags in dipolar_nn.f"
         write(0,*)"max_frags: ",max_frags," vall_size: ",vall_size
         stop
      endif


      last=first+size-1
      mapstart=first        !save this to pass via common block

      do iset=1,Ndipolar_sets
         maplength(iset) =0
      enddo
      do i=1,Ndipolar
         if (pairDipolar(1,i).gt.last) goto 250 !return
         if (pairDipolar(1,i).lt.first) goto 200
         if (pairDipolar(3,i).lt.first.or. pairDipolar(3,i).gt.last)
     #        goto 200
         if ((pairDipolar(1,i).eq.first.and.pairDipolar(2,i).eq.0)
     #        .or.(pairDipolar(3,i).eq.first.and.
     #        pairDipolar(4,i).eq.0))
     #        goto 200          !HN of first not determined
         if (dipolartype(pairDipolar(1,i)).eq.0) goto 200 !don't have
                                !data for undefined types

         iset=pairDipolar(5,i)
         maplength(iset)=maplength(iset)+1 !acceptable constraint
         map(maplength(iset),iset)=i
 200     continue
      enddo

car calculate significance score (zscore)

car score all frags, calculate mean
 250  mean=0.0
      nfrag=0
      do i=1,vall_size-size+1
        frag_dipolar_score(i)=score_frag_dipolar(i)
        if (frag_dipolar_score(i).lt.1.e9 .and.
     #       frag_dipolar_score(i).gt.0.0) then
           nfrag=nfrag+1
           mean=mean+frag_dipolar_score(i)
         endif
      enddo
      if (nfrag.eq.0) return
      mean=mean/nfrag

car calculate std dev (rmsd)

      sd=0.0
      do i=1,vall_size-size+1
         if (frag_dipolar_score(i).lt.1e9. and.
     #       frag_dipolar_score(i).gt.0.0) then
            sd=sd+(mean-frag_dipolar_score(i))**2
         endif
      enddo
      sd=sqrt(sd/nfrag)

car weight so than mean score is 20:   long tail above mean, most
car scores much less than 20
      weight=20.0/mean

      do i=1,vall_size-size+1
         frag_dipolar_score(i)=frag_dipolar_score(i)*weight
         if (frag_dipolar_score(i).eq.0.0)frag_dipolar_score(i)=mean*weight
      enddo


      return
      end
c--------------------------------------------------------------------------
      subroutine assemble_datamatrix_nn(fragstart,mapstart,map,nrow,A,b)

car assembles data matrices A and b for singular value decomposition
car using data in the common blocks as specified by the map generated
car by prepare_dipolar_score
car
car nrow is the number of rows in A and the size of b and length of map
car map is a lookup table to access the data in the global data arrays
car pairDipolar and Jdipolar

      implicit none
      include 'param.h'
      include 'structure.h'    !vall_size
      include 'dipolar.h'

car input
      integer fragstart   !first residue of vall in fragment
      integer mapstart    !first query residue for which map is built
      integer map(MAXDIPOLAR)  !map to data for corrent data set
      integer nrow      !number of data points in map

car common block input
      real*8 vall_svdcoeff(ORDERSIZE,MAXDIPTYPES,MAX_VALL)
                       !(data,type,residue)
                       !data:  1-5:direction cosines products
                       !dipolar types for which data is precomputed:
                       !1:HN-N,2:HA-CA,3:CA-C,4:C(i-1)-N,5:C(i-1)-H
      common /dipolar_valldata/vall_svdcoeff


car output
      real*8 A(MAXDIPOLAR,ORDERSIZE)
      real*8 b(MAXDIPOLAR)

car internal
      integer i
      integer type,res

car function calls
      integer dipolartype

      do i=1,nrow
         type=dipolartype(pairdipolar(1,map(i)))
         b(i) =Jdipolar(map(i)) !already distance normalized
         res=fragstart+(pairDipolar(3,map(i))-mapstart)   !vall residue
         A(i,1)=vall_svdcoeff(1,type,res)
         A(i,2)=vall_svdcoeff(2,type,res)
         A(i,3)=vall_svdcoeff(3,type,res)
         A(i,4)=vall_svdcoeff(4,type,res)
         A(i,5)=vall_svdcoeff(5,type,res)
      enddo
      return
      end

c--------------------------------------------------------------------------
      subroutine delete_dipolar()

cmj function to delete all dipllar couplings
      implicit none
      include 'dipolar.h'

      dipolar_exist = .false.
      Ndipolar_sets = 0
      Ndipolar = 0

      return
      end






