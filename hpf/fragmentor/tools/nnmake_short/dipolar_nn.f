c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 1.7 $
c  $Date: 2002/06/18 04:41:53 $
c  $Author: rohl $

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
      if (.not.dipolar_exist) return
      if (verbose) write(0,*)"residue, Jobs, Jcalc, residual**2"
      do iset=1,ndipolar_sets
         if (verbose) write(0,*)"dipolar constraint set ",iset
         nrow=maplength(iset)
         call assemble_datamatrix_nn(fragstart,mapstart,map(1,iset),
     #        nrow,A,b)
         call calc_ordermatrix(nrow,A,b,x(1,iset),reject(iset))
         if (reject(iset)) then
            write(0,*)"iset=",iset," fragstart=",fragstart
            write(0,*)'rejected' 
            write(0,*)x(1,iset),x(2,iset),x(3,iset),x(4,iset),x(5,iset)
            do i=1,maplength(iset)
               write(0,'(i,6f6.3)')i,A(i,1),A(i,2),A(i,3),A(i,4),A(i,5),
     #              b(i)
            enddo
            score_frag_dipolar=1e9
            return
c            goto 300
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

      real mean,sd,min,weight

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

