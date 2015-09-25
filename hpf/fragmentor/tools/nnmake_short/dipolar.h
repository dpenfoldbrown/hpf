c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 1.4 $
c  $Date: 2002/06/18 04:41:53 $
c  $Author: rohl $


      integer MAXDIPOLAR
      parameter (MAXDIPOLAR=900)
      integer MAXDIPSETS
      parameter (MAXDIPSETS=2)
      integer ORDERSIZE
      parameter (ORDERSIZE=5)    !# independent elements in ordermatrix
      integer MAXDIPTYPES        !types for which data is precomputed in vall
      parameter (MAXDIPTYPES=6)  !1:HN-N,2:HA-CA,3:CA-C,4:C(i-1)-N,5:C(i-1)-HN,HN(i)-HA(i)

      real dipolar_cutoff
      parameter (dipolar_cutoff=30.0)  !score cutoff for full atom

      integer Ndipolar          !number of dipolar constraints
      integer pairDipolar(5,MAXDIPOLAR)      
                                !(key,index) key:1=res1,2=atom1
                                !                3=res2,4=atom2
                                !                5=set#
      integer  Ndipolar_sets
      real*8 Jdipolar(MAXDIPOLAR) !coupling constant of ith constraint (Hz)
      real*8 invDcnst(MAXDIPOLAR) !1/dipolar interaction constant
      real Ddist(MAXDIPOLAR)    !dist for atom pairs whose distance is fixed
                                !<0 for atoms whose distance is variable
      real*8 invDmax(MAXDIPOLAR)  !dist**3/Dcnst/dist**3 for fixed dist pairs
                                !when Ddist2<0 (ie variable dist) invDmax
                                !is equal to the last calculated value
      integer dunit            !unit number to dump output to

      common /dipolar_private/ Ndipolar,pairDipolar,Ndipolar_sets,Jdipolar,
     #                    invDcnst,Ddist,invDmax,dunit



      logical dipolar_exist 

      common /dipolar_public/ dipolar_exist

