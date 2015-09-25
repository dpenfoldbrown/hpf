



c===============================================================================
c
      subroutine maxsub(ires1,ires2,nsec1,nsec2,nvnew,iv1,iv2,
     &                  xca1,xca2,norm,match,lsup,sec1,sec2,zthresh,generate_pdb,
     &                  zscore,evalue,scorel,nssmatch,
     &                  norm_o1,nsup_o1,rms_o1,psi_o1,
     &                  norm_o2,nali_o2,rms_o2,psi_o2,
     &					co_1,co_2,COMPUTE_CO)
    
c-------------------------------------------------------------------------------
c --- Here applies a modification of Dani Fischer's heuristic algorithm for 
c --- finding the largest subset of residues for superimposition within
c --- a threshold. The part that restraints the secondary structure
c --- matching needs to be changed.
c ---
c --- At this point, the algorithm works as follows: first, the residue 
c --- assignment is done on the basis of the global secondary structure
c --- similarity. Then, with this assignment ofr residue pairs, the
c --- heuristic procedure of Fisher is used.
c-------------------------------------------------------------------------------

      implicit REAL*8 (a-h,o-z)
      
      real*8 rms_o1,psi_o1,rms_o2,psi_o2
      integer norm_o1,nsup_o1,norm_o2,nali_o2
	  integer COMPUTE_CO
	  real*8 co_f1,co_f2    ! contact order
        

c --- Vector alignment
      integer maxres,maxlen,maxfrag
      real*8 rmstol,angmax
    
      parameter  (maxres = 3000)
      parameter  (maxlen = 2*maxres)
      parameter  (rmstol = 4.0d0)
      parameter  (maxfrag=  100)
      parameter  (angmax = 60.0d0)
      
      real*8 xca1,xca2,xp,xe,xp0,xe0,w,wmax
      dimension    xca1(maxres,3),xca2(maxres,3)
      dimension    xp(3*maxres),xe(3*maxres),xp0(3*maxres),xe0(3*maxres)
      dimension    w(maxres),wmax(maxres)

      character*1  sec1(maxres), sec2(maxres)
      character*1  match(maxlen),lsup(maxlen)
      character*3  rnam1(maxres),rnam2(maxres),aaname3(0:20)

      logical      assigned1(maxres),assigned2(maxres)
      logical      matched(maxres),logical_w(maxres)

      integer      iv1(maxlen),iv2(maxlen),ial(maxlen),smax
      integer      nsec1(maxlen),nsec2(maxlen),np(maxres)
      integer      ires1(maxres),ires2(maxres)
      integer      ir1(maxres),ir2(maxres)
      integer      block1(maxfrag,2),block2(maxfrag,2)
      integer      pairs(maxres,2),isup1(maxres),isup2(maxres)
      
      integer last
      real*8 min_d2,d2
	  
c -- contact order	  
	  integer N_co1,n_co2  ! count of contacts
	  real*8 co_1,co_2     ! contact order
	  real*8 CO_CUTOFF2
	  parameter (CO_CUTOFF2 = 64.0)
      
	  
c --  output control
      real*8 zthresh
      logical generate_pdb
c -- functions     


      COMMON /TRANSFORM/ XPC,YPC,ZPC,XEC,YEC,ZEC,R(3,3)


      data aaname3 /'ALA','CYS','ASP','GLU','PHE','GLY','HIS',
     &              'ILE','LYS','LEU','MET','ASN','PRO','GLN',
     &              'ARG','SER','THR','VAL','TRP','TYR','XXX'/


c-------------------------------------------------------------------------------
c  First selects the residues allowed to be superimposed
c-------------------------------------------------------------------------------

      nsup = 0
      nssmatch =0
      k    = 0
      do i = 1, nvnew
         if (iv1(i).ne.0 .and. iv2(i).ne.0) then
            nsup = nsup + 1
            if (match(i).eq.'|') then 
               nssmatch = nssmatch+1
               matched(nsup) = .true.
            else 
               matched(nsup) = .false.
            endif
            ial(nsup) = i
            do j = 1, 3
               k=k+1
               xp(k) = xca1(iv1(i),j)
               xe(k) = xca2(iv2(i),j)
               
            enddo
            w(nsup) = 1.0d0
         endif
      enddo
      
      l = 7
      
      if (nsup.gt.2)  call rmsfitca2(nsup,xp,xe,w,nsup,rms) ! allignment of less than 2 residues?
      
      psi = (float(nsup)/float(norm))*100.0d0
      
      psi_o1 = psi
      nsup_o1 = nsup
      norm_o1 = norm
      rms_o1 = rms
      
      
C$$$       write(6,*)
C$$$       write(6,21)' PSI(ini)= ',psi,'NALI= ',nsup,'NORM= ',norm,
C$$$      &     'RMS= ',rms,'NSS= ',nssmatch
     
      if (nsup.le.2) then
         evalue= 1.0
         Zscore= -10.0
         psi_o2 = psi
        nali_o2 = nsup
        norm_o2 = norm
        rms_o2 = rms
         
C$$$          Write(6,*)"E-value=     1.0"    
C$$$          write(6,*) 
C$$$          write(6,*)"Z-score=   -10.0        -ln(E)=      0.000" 
C$$$          Write (6,*)  "Z-score: -10.0  "
       
         write (6,*)  "WARNING: aborting maxsub, too few superimposable residues  ",nsup," out of ",nvnew
         return 
      endif

      do i=1,nsup*3
         xe0(i)=xe(i)
         xp0(i)=xp(i)
      enddo
c-------------------------------------------------------------------------------
c  Now apply Fishers's maxsub algorithm. An heptapeptide is used.
c  The algorithm is modified in that only pairs of residues with
c  similar local secondary structures are allowed to be used as 
c  seed residues to superimpose.
c------------------------------------------------------------------------------
      smax = 0
      do i = 1, nvnew           ! should this not be 1,3*nsup?
                                !w(i)    = 0.0d0
         wmax(i) = 0.0d0
      enddo
      rmsmax=1000
      
      if (nssmatch.lt.2) goto 77 ! basically aborting but outputing some diagnostics too 

      do 5 i = 1, nsup-l+1
         call clear_rms()
         if (matched(i)) then  ! this line is Angel's variation on Danni Fischer's method.
                               ! only seed at points with matched SS.
                               
            ! build up a seed segment of length l                   
            lmax = 0
            do j = 1, nsup 
               if (j.ge.i .and. j.le.i+l-1) then ! could do this without if statement
                  call add_rms(j,xp0,xe0)
                  logical_w(j)=.TRUE. !w(j) = 1.0d0
                  lmax = lmax + 1
               else
                  logical_w(j)=.FALSE. ! w(j) = 0.0d0
               endif
            enddo
            
            
            ! find initial alignment using seed
            ! rmsfitca3 rotates all of the residues but the rotation only alligns the 
            ! residues pushed into add_rms.
            call rmsfitca3(nsup,xp0,xp,xe0,xe,rms)   
            
            ! next we iterate the following algorithm 
            !  1) using current allignment find all atoms within a threshold t of being superimposed
            !  2) add these close ones to set to be aligned
            !  3) aling using the current set, then reorient all atoms.
            !  4) increment t by a small amount.
            !  5) repeat this until theshold= 7 angstroms.                    
            t=0.0d0
            last = lmax
            do while (t.lt.7.0d0)  !do m = 1, 7 
               last = lmax
               t=t+1.0d0        ! increment threshold by one angstrom. (note this is ties to int(min_d))
               t2 = t*t         ! t squared
               min_d2 = 8.0*8.0 
                                !t = float(m)   !*7.0/float(7)  ! huh? must be a relic??? 
               do n = 1, nsup
                  if (.not.logical_w(n)) then ! if (w(n).eq.0.0d0) then
                     
                     k = 3*(n-1)
                     d2 = ((xp(k+1)-xe(k+1))**2+(xp(k+2)-xe(k+2))**2+(xp(k+3)-xe(k+3))**2)
                         ! squared distance
                     if (d2.le.t2) then  ! is this atom within threshold?
                        call add_rms(n,xp0,xe0)  ! if so, add to list
                        logical_w(n)=.TRUE. !w(n) = 1.0d0  ! set membership flag
                        lmax = lmax + 1          ! keep a count of members
                     else 
                        ! if not below threshold, then find the closest atom
                        if (d2.le.min_d2) min_d2=d2 !min_d = min(d,min_d)
                     endif
                  endif
               enddo
                                ! check if we added any residues on this iteration.
                                ! if not then 1) we dont need to refit the calphas 2) we can advance threshold level 
               if (lmax.ne.last) then
                  call rmsfitca3(nsup,xp0,xp,xe0,xe,rms)
               else
                  
                                !write(0,*) i,'skipping',t,int(min_d),lmax,min_d
                  t=int(sqrt(min_d2))  ! advance the threshold
               endif
            enddo
            
                                ! huh? logic here is confusing.
            if ((lmax.gt.smax).and.(rms.le.rmstol)) then
               smax = lmax
               do n = 1, nsup
                  if (logical_w(n) ) then
                     wmax(n) = 1.0d0
                  else 
                     wmax(n) = 0.0d0
                  endif
                                ! wmax(n) = w(n)
               enddo
               rmsmax=rms
            else if ((lmax.eq.smax).and.(rms.lt.rmsmax)) then
               
               smax = lmax
               do n = 1, nsup
                  if (logical_w(n) ) then
                     wmax(n) = 1.0d0
                  else 
                     wmax(n) = 0.0d0
                  endif
                                !wmax(n) = w(n)
               enddo
               rmsmax=rms
            endif
         endif
 5    continue
      
c-------------------------------------------------------------------------------
c     --- Confirm final superimposition.
c     --- first, compile regions without indels. Report rms
c-------------------------------------------------------------------------------
 77   continue
      k  = 0
      nr = 0
      do 10 i = 1, nvnew
         lsup(i)=' '
         if (iv1(i).eq.0.or.iv2(i).eq.0) goto 10
         nr = nr + 1
         rnam1(nr) = aaname3(nsec1(iv1(i)))
         rnam2(nr) = aaname3(nsec2(iv2(i)))
         if (int(wmax(nr)) .eq. 1) lsup(i)='*'
         do j = 1, 3
            k=k+1
            xp(k) = xca1(iv1(i),j)
            xe(k) = xca2(iv2(i),j)
         enddo
		 kk = 0
				 
         ir1(nr) = ires1(iv1(i))
         ir2(nr) = ires2(iv2(i))
 10   continue
      
      if (smax.gt.1) then         
         call rmsfitca2(nr,xp,xe,wmax,smax,rms) ! side effect sets xpc,zpc etc... via common
      else
      ! if smax is less than 2 then basically we failed to find an allignement
      !  to make the best of a bad situation we simply revert to alligning all of the residue
         do i =1,nr
            wmax(i) = 1.0
         enddo
         call rmsfitca2(nr,xp,xe,wmax,nr,rms) ! side effect sets xpc,zpc etc... via common
      endif
c-------------------------------------------------------------------------------
c compute the contact order if requested
c-------------------------------------------------------------------------------



        co_1=0.0
		n_co1=0
		co_2=0.0
		n_co2=0
		
	if (COMPUTE_CO.eq.1) then
	  do 22 i = 1, nvnew

		 if (int(wmax(i)) .ne. 1)  goto 22  ! not included
		

		do 23 j = i+1,nvnew

		! note this is dumb in the sense that we do twice as many distance calc than are required.
		! since each pair gets computed once for first atom and once for the second atom
		! we also loop over more than we need to since it's always over nvnew rather than actual lengths
       
	   ! now we have to do the contact order of the overlap region on both the predicted and experimental structs
		if (int(wmax(j)) .ne. 1)  goto 22  ! not included in overlap
	    if (iv1(i).ne.0.and.iv1(j).ne.0) then  ! when would this be true? does wmax() already contain this check?
		 if (iv1(j)-iv1(i) .ge. 4) then ! only do i to i+4
			d2 = 
     &		        (xca1(iv1(i),1)-xca1(iv1(j),1))**2+
     &				(xca1(iv1(i),2)-xca1(iv1(j),2))**2+
     &				(xca1(iv1(i),3)-xca1(iv1(j),3))**2

			if (d2.le.CO_CUTOFF2) then  ! this is squareof distance cutoff in angstroms
			co_1 = co_1+ abs(iv1(j)-iv1(i))
			n_co1 = n_co1+1

			endif
			endif
		endif
	    if (iv2(i).ne.0.and.iv2(j).ne.0) then  ! when would this be true? does wmax() already contain this check?
          		 if (iv2(j)-iv2(i) .ge. 4) then ! only do residues separated by three residues.
			d2 = 
     &		        (xca2(iv2(i),1)-xca2(iv2(j),1))**2+
     &				(xca2(iv2(i),2)-xca2(iv2(j),2))**2+
     &				(xca2(iv2(i),3)-xca2(iv2(j),3))**2
	 
		   if (d2.le.CO_CUTOFF2) then  ! this is squareof distance cutoff in angstroms
			co_2 = co_2+ abs(iv2(j)-iv2(i))
			n_co2 = n_co2+1

			endif
			endif
		endif
			   
23		continue
22		continue
       if (n_co1 .gt. 0) then
			co_1 = co_1/n_co1
	   endif
	   
		if (n_co2 .gt. 0) then
			co_2 = co_2/n_co2
	   endif
	   
		endif  ! COMPUTE_CO test
		
	  
	  
c     write(6,'(/,a,f7.3,a,i5,/)') 'RMS = ', rms, ' SMAX = ', smax
      psi = (float(smax)/float(norm))*100.0d0

c-------------------------------------------------------------------------------
c --- Transform PSI to Levitt & Gerstein Sstr
c --- NEED TO BE REMOVED
c-------------------------------------------------------------------------------

c compute score without gaps

      score = 0.0d0
      nali  = 0
      ngap  = 0
      do i = 1, nvnew
        if (iv1(i).eq.0.or.iv2(i).eq.0) ngap=ngap+1
        if (wmax(i).eq.1) then
          d = 0.0d0
          k = 3*(i-1)  ! fixed bug.  was 3*(n-1)
          do j = 1, 3
            d = d + (xp(k+j)-xe(k+j))**2
          enddo
          d = dsqrt(d)
          score = score + 1.0d0/(1.0d0+(d/rmstol)**2)
          nali  = nali + 1
        endif
      enddo
c remove N and C terminal gaps

      nend = 0
      do i = 1, nvnew
        if (iv1(i).ne.0.and.iv2(i).ne.0) goto 15
        nend = nend + 1
      enddo
   15 do i = nvnew, 1, -1
        if (iv1(i).ne.0.and.iv2(i).ne.0) goto 20
        nend = nend + 1
      enddo
   20 ngap = ngap - nend

c correct score

      scoren = score + (ngap/2)
      scorel = 20 * scoren

c-------------------------------------------------------------------------------
c --- All done. Report match statistics
c-------------------------------------------------------------------------------
      psi_o2 = psi
      nali_o2 = nali
      norm_o2 = norm
      rms_o2 = rms
C$$$       write(6,21)' PSI(end)= ',psi,'NALI= ',nali,'NORM= ',norm,
C$$$      &           'RMS= ',rms

C$$$ 
C$$$       write(6,21)' Sstr(LG)= ',scorel,'NALI= ',nali,'NORM= ',norm,
C$$$      &           'RMS= ',rms

c     write(6,*) scorel, scoren, ngap, nend, nali, norm
  21  format(a,f7.2,2x,2(a,i3,2x),a,f7.2,2x,a,i3)
c-------------------------------------------------------------------------------
c --- Report probabilities for random matches. 
c --- These are preliminary values for the average and standard deviation
c --- of PSI as a function of NORM. These values are:
c --- m(L) = 759.31 * L**(-0.7545)
c --- s(L) = 393.32 * L**(-0.9009)
c-------------------------------------------------------------------------------

c     am = 759.31 * norm**(-0.7545)
c     as = 393.32 * norm**(-0.9009)
c     am = 695.41 * norm**(-0.7278)
c     as = 340.00 * norm**(-0.9045)
c EV fitting, using N>70
      am = 747.29 * nr**(-0.7971)
      as = 124.99 * nr**(-0.6882)
      zscore = (psi-am)/as
c this is the gaussian approach. Actually, extreme-value is more adequate
c      evalue = 0.5d0 * erfcc(zscore/sqrt(2.0d0))
c here it is the EV approach
      znew   = 0.730*((1.2825755*zscore)+0.5772)
      evalue = 1.0d0-exp(-exp(-znew))
c the following change was made on may 18 2005
C  we used to truncate the e-value at 2.65E-14 because the round-off error was large for better values.
c  now we will use a taylor approximation of the function expanded near zero

c due to numerical errors, the e-value is cutoff at 2.650d-14
c  if (evalue.lt.2.650d-14) evalue=2.650d-14  ! old code
      if (evalue.lt.2.650d-14) then
             evalue = exp(-znew)
             evalue = evalue*(1.0d0 - evalue/2.0d0) 
      endif

      
      
      
 
C$$$     
C$$$       write(6,*)
C$$$       write(6,"(a,1x,g17.8)") 'E-value= ', evalue
C$$$       write(6,*)
C$$$       write(6,"(a,1x,g17.8,1x,a,1x,g17.8)") 'Z-score= ', zscore, '   -ln(E)= ', -log(evalue)

c-------------------------------------------------------------------------------
c --- This now is for the output of the superimposed coordinates
c --- now, apply transformation to indels as well
c --- It is always assumed that xp is the model and xe the experiment,
c --- and that transformation is based upon rotation of xe.
c --- Watch out for last bead!!! (not yet incorporated)
c-------------------------------------------------------------------------------
      if (generate_pdb .and. zscore.ge.zthresh ) then  
     
      nsup1=nsup
      nsup2=nsup
      do i = 1, nvnew

        if (iv2(i).eq.0) then

          l=iv1(i)
          ns3 = 3*nsup1
          xp(ns3+1) = xca1(l,1)-xpc
          xp(ns3+2) = xca1(l,2)-ypc
          xp(ns3+3) = xca1(l,3)-zpc
          nsup1=nsup1+1
          rnam1(nsup1) = aaname3(nsec1(l))
          ir1(nsup1)   = ires1(l)

        else if (iv1(i).eq.0) then

          l=iv2(i)
          ns3 = 3*nsup2
          xt = xca2(l,1)-xec
          yt = xca2(l,2)-yec
          zt = xca2(l,3)-zec
          xe(ns3+1) = r(1,1)*xt+r(1,2)*yt+r(1,3)*zt
          xe(ns3+2) = r(2,1)*xt+r(2,2)*yt+r(2,3)*zt
          xe(ns3+3) = r(3,1)*xt+r(3,2)*yt+r(3,3)*zt
          nsup2=nsup2+1
          rnam2(nsup2) = aaname3(nsec2(l))
          ir2(nsup2) = ires2(l)

        endif

      enddo
 
c-------------------------------------------------------------------------------
c --- finally, write everything in pdb format. Write also a header
c --- with the key numerical results from the calculation. For the pdb,
c --- first sort the residue number to give a better output for graphics
c-------------------------------------------------------------------------------

      open(unit=40,file='maxsub_sup.pdb')

      write(40,'(a)')      'REMARK '
      write(40,'(a,f7.2)') 'REMARK PSI(end)= ',psi
      write(40,'(a,i7)')   'REMARK NALI    = ',nali
      write(40,'(a,i7)')   'REMARK NORM    = ',norm
      write(40,'(a,f7.4)') 'REMARK RMS     = ',rms
      write(40,'(a)')      'REMARK '
      write(40,'(a,f7.2)') 'REMARK Z-score = ',zscore
      write(40,'(a,f7.2)') 'REMARK -ln(E)  = ',-log(evalue)
      write(40,'(a)')      'REMARK '
      write(40,'(a)')      'REMARK Transformation Matrix:'
      write(40,'(a,3(f9.5,1x))') 'REMARK ',r(1,1),r(1,2),r(1,3)
      write(40,'(a,3(f9.5,1x))') 'REMARK ',r(2,1),r(2,2),r(2,3)
      write(40,'(a,3(f9.5,1x))') 'REMARK ',r(3,1),r(3,2),r(3,3)
      write(40,'(a)')      'REMARK '
      write(40,'(a)')      'REMARK Translation vector (Prediction):'
      write(40,'(a,3(f9.5,1x))') 'REMARK ',xpc,ypc,zpc
      write(40,'(a)')      'REMARK '
      write(40,'(a)')      'REMARK Translation vector (Experiment):'
      write(40,'(a,3(f9.5,1x))') 'REMARK ',xec,yec,zec
      write(40,'(a)')      'REMARK '

      call sort(nsup1,ir1,xp,rnam1)
      call putpdb(nsup1,xp,40,rnam1,ir1)

      call sort(nsup2,ir2,xe,rnam2)
      call putpdb(nsup2,xe,40,rnam2,ir2)

c-------------------------------------------------------------------------------
c --- Test for the refinement step
c --- Here also writes the superimposition in a script for rasmol
c --- NEED TO MAKE SUPERIMPOSITION WITH ALL RESIDUES SELECTED
c --- PUT AT THE END OF THE SUBROUTINE
c-------------------------------------------------------------------------------
      
      npairs = 0
      do i = 1, nvnew
        assigned1(i)=.false.
        assigned2(i)=.false.
      enddo
      do i = 1, nvnew
        if (lsup(i).eq.'*') then
          npairs=npairs+1
          pairs(npairs,1)=iv1(i)
          pairs(npairs,2)=iv2(i)
          assigned1(iv1(i))=.true.
          assigned2(iv2(i))=.true.
          np(npairs)=i
        endif
      enddo

c This needs some refinement yet... Not clear that this is the best
c procedure. Perhaps DP using something similar to Levitt&Gerstein...

      npal=npairs
      itimes=0
c      do while (itimes.lt.3)
c        call extendali(nvnew,lsup,nsup1,nsup2,iv1,iv2,
c     &                 assigned1,assigned2,npairs,pairs,xp,xe)
c        call extendcheck(npairs,nvnew,w,xca1,xca2,xp,xe,pairs,iv1,iv2,
c     &                   nsec1,nsec2,ires1,ires2,ir1,ir2,rnam1,rnam2)
c        itimes=itimes+1
c      enddo

c-------------------------------------------------------------------------------
c now we need to update the final structural alignment. We need to
c update iv1 and iv2...
c
c Now executes the new structural alignment...
c-------------------------------------------------------------------------------

      open(30,file='rasmol.tcl')

      psi = (float(npairs)/float(norm))*100.0d0
      write(6,*)
      write(30,'(a,i5)') '# itimes = ', itimes
      write(30,21)'# PSI(end)= ',psi,'NALI= ',npairs,'NORM= ',norm,
     &           'RMS= ',rms

      open(unit=50,file='maxsub_sup2.pdb')

c apply transformation for the first protein.

      nsup1=0
      k=0
      do i = 1, nvnew
        l=iv1(i)
        if (l.ne.0) then
          xp(k+1) = xca1(l,1)-xpc
          xp(k+2) = xca1(l,2)-ypc
          xp(k+3) = xca1(l,3)-zpc
          nsup1=nsup1+1
          rnam1(nsup1) = aaname3(nsec1(l))
          ir1(nsup1) = ires1(l)
          k=k+3
        endif
      enddo
      call sort(nsup1,ir1,xp,rnam1)
      call putpdb(nsup1,xp,50,rnam1,ir1)

c apply transformation for the second protein.

      nsup2=0
      k=0
      do i = 1, nvnew
        l=iv2(i)
        if (l.ne.0) then
          xt = xca2(l,1)-xec
          yt = xca2(l,2)-yec
          zt = xca2(l,3)-zec
          xe(k+1) = r(1,1)*xt+r(1,2)*yt+r(1,3)*zt
          xe(k+2) = r(2,1)*xt+r(2,2)*yt+r(2,3)*zt
          xe(k+3) = r(3,1)*xt+r(3,2)*yt+r(3,3)*zt
          nsup2=nsup2+1
          rnam2(nsup2) = aaname3(nsec2(l))
          ir1(nsup2) = ires2(l)
          k=k+3
        endif
      enddo
      call sort(nsup2,ir2,xe,rnam2)
      call putpdb(nsup2,xe,50,rnam2,ir2)

      close(50)

c-------------------------------------------------------------------------------
c --- calculate the new e-value:
c-------------------------------------------------------------------------------

      psi    = (float(npairs)/float(norm))*100.0d0
      zscore = (psi-am)/as
      znew   = 0.730*((1.2825755*zscore)+0.5772)
      evalue = 1.0d0-exp(-exp(-znew))
      if (evalue.lt.2.650d-14) evalue=2.650d-14
      write(30,'(a,f10.6)') '# -ln(E) = ', -log(evalue)
      write(30,'(a,i4)')    '# npairs = ', npairs
      do i = 1, npairs
        write(30,'(a,2i4)') '# ', pairs(i,1), pairs(i,2)
      enddo

      write(30,'(a)') 'load inline'
      write(30,'(a)') 'select :A'
      write(30,'(a)') 'color [100,149,237]'
      write(30,'(a)') 'backbone 50'
      write(30,'(a)') 'select:B'
      write(30,'(a)') 'colour [255,20,147]'
      write(30,'(a)') 'backbone 50'

      do i = 1, nvnew
        isup1(i)=0
        isup2(i)=0
      enddo
      do i = 1, npairs
        isup1(pairs(i,1)) = 1
        isup2(pairs(i,2)) = 1
      enddo

      nb1=0
      do i = 1, nsup1-1
        if (i.eq.1.and.isup1(i).eq.1) then
          nb1=nb1+1
          block1(nb1,1)=i
        elseif (i.eq.1.and.isup1(i).eq.0.and.i+1.eq.2.and.isup1(i).eq.1)
     &  then
          nb1=nb1+1
          block1(nb1,1)=i+1
        endif
        if (isup1(i).eq.0.and.isup1(i+1).eq.1) then
          nb1=nb1+1 
          block1(nb1,1)=i+1
        endif
        if (isup1(i).eq.1.and.isup1(i+1).eq.0) block1(nb1,2)=i
      enddo
      if (isup1(nsup1).eq.1) block1(nb1,2)=nsup1
      do i = 1, nb1
        write(30,31) 'select ',block1(i,1),'-',block1(i,2),':A'
        write(30,'(a)') 'backbone 150'
      enddo

      nb2=0
      do i = 1, nsup2-1
        if (i.eq.1.and.isup2(1).eq.1) then
          nb2=nb2+1
          block2(nb2,1)=i
        elseif (i.eq.1.and.isup2(i).eq.0.and.i+1.eq.2.and.isup2(i).eq.1) 
     &  then
          nb2=nb2+1
          block2(nb2,1)=i+1
        endif
        if (isup2(i).eq.0.and.isup2(i+1).eq.1) then
          nb2=nb2+1 
          block2(nb2,1)=i+1
        endif
        if (isup2(i).eq.1.and.isup2(i+1).eq.0) block2(nb2,2)=i
      enddo
      if (isup2(nsup2).eq.1) block2(nb2,2)=nsup2
      do i = 1, nb2
        write(30,31) 'select ',block2(i,1),'-',block2(i,2),':B'
        write(30,'(a)') 'backbone 150'
      enddo

      write(30,'(a)') 'select all'
      write(30,'(a)') 'wireframe off'
      write(30,'(a)') 'select none'
      write(30,'(a)') 'exit'
      j=0
      do i = 1, nsup1
c       write(30,500) i,rnam1(i),ir1(i),xp(j+1),xp(j+2),xp(j+3)
        write(30,500) i,rnam1(i),i,xp(j+1),xp(j+2),xp(j+3)
        j=j+3
      enddo
      j=0
      do i = 1, nsup2
c       write(30,501) i,rnam2(i),ir2(i),xe(j+1),xe(j+2),xe(j+3)
        write(30,501) i,rnam2(i),i,xe(j+1),xe(j+2),xe(j+3)
        j=j+3
      enddo
   31 format(a,2(i4,a))
  500 format('ATOM ',i6,1x,' CA ',1x,a3,' A',i4,4x,3f8.3)
  501 format('ATOM ',i6,1x,' CA ',1x,a3,' B',i4,4x,3f8.3)

c-------------------------------------------------------------------------------
c --- End of the program
c-------------------------------------------------------------------------------

      close(30)
      endif
      return
      end
c
c===============================================================================
c
      function angle(x,y)

c-------------------------------------------------------------------------------
c --- calculate the angle between two vectors (given in degrees)
c-------------------------------------------------------------------------------

      implicit REAL*8 (a-h,o-z)

      parameter (radian=57.29577951308232088d0)
      dimension x(3),y(3)

      xmod=0.0d0
      do k = 1, 3
        xmod = xmod + x(k)*x(k)
      enddo
      xmod = dsqrt(xmod)

      ymod=0.0d0
      do k = 1, 3
        ymod = ymod + y(k)*y(k)
      enddo
      ymod = dsqrt(ymod)

      q = xmod*ymod
      if (q.lt.1.0d-7) q=1.0d-7

      xy = 0.0d0
      do k = 1, 3
        xy = xy + x(k)*y(k)
      enddo

      cos_angle = xy/q
      if (abs(cos_angle) .ge. 1.0d0)  cos_angle = sign(1.0d0,cos_angle)

      angle = radian * acos(cos_angle)

      return
      end
c
c===============================================================================
c
      subroutine extendcheck(npairs,nvnew,w,xca1,xca2,xp,xe,pairs,
     &                       iv1,iv2,nsec1,nsec2,ires1,ires2,
     &                       ir1,ir2,rnam1,rnam2)

c-------------------------------------------------------------------------------
c --- Checks the extension of the alignment
c-------------------------------------------------------------------------------

      implicit REAL*8 (a-h,o-z)

      parameter  (maxres = 3000)
      parameter  (maxlen = 2*maxres)
      parameter  (rmstol = 4.0d0)

      dimension    xca1(maxres,3),xca2(maxres,3)
      dimension    xp(3*maxres),xe(3*maxres),w(maxres)

      character*1  lsup(maxlen)
      character*3  rnam1(maxres),rnam2(maxres),aaname3(0:20)

      integer      iv1(maxlen),iv2(maxlen)
      integer      nsec1(maxlen),nsec2(maxlen),np(maxres)
      integer      ires1(maxres),ires2(maxres),nptmp(maxres)
      integer      ir1(maxres),ir2(maxres)
      integer      pairs(maxres,2)

      COMMON /TRANSFORM/ XPC,YPC,ZPC,XEC,YEC,ZEC,R(3,3)

      data aaname3 /'ALA','CYS','ASP','GLU','PHE','GLY','HIS',
     &              'ILE','LYS','LEU','MET','ASN','PRO','GLN',
     &              'ARG','SER','THR','VAL','TRP','TYR','XXX'/

      nptrial=npairs
      do i = 1, nvnew
        w(i) = 0.0d0
      enddo
      k=0
        do i = 1, npairs
          do j = 1, 3
            k=k+1
            xp(k) = xca1(pairs(i,1),j)
            xe(k) = xca2(pairs(i,2),j)
          enddo
          w(i) = 1.0d0
        enddo
       call rmsfitca2(npairs,xp,xe,w,npairs,rms)
 
        nsup1=0
        k=0
        do i = 1, nvnew
          l=iv1(i)
          if (l.ne.0) then
            xp(k+1) = xca1(l,1)-xpc
            xp(k+2) = xca1(l,2)-ypc
            xp(k+3) = xca1(l,3)-zpc
            nsup1=nsup1+1
            rnam1(nsup1) = aaname3(nsec1(l))
            ir1(nsup1) = ires1(l)
            k=k+3
          endif
        enddo
        call sort(nsup1,ir1,xp,rnam1)
        nsup2=0
        k=0
        do i = 1, nvnew
          l=iv2(i)
          if (l.ne.0) then
            xt = xca2(l,1)-xec
            yt = xca2(l,2)-yec
            zt = xca2(l,3)-zec
            xe(k+1) = r(1,1)*xt+r(1,2)*yt+r(1,3)*zt
            xe(k+2) = r(2,1)*xt+r(2,2)*yt+r(2,3)*zt
            xe(k+3) = r(3,1)*xt+r(3,2)*yt+r(3,3)*zt
            nsup2=nsup2+1
            rnam2(nsup2) = aaname3(nsec2(l))
            ir1(nsup2) = ires2(l)
            k=k+3
          endif
        enddo
        call sort(nsup2,ir2,xe,rnam2)
        npafter=0
        do i = 1, npal
          k=(3*pairs(i,1))-1
          xpi=xp(k+1)
          xpj=xp(k+2)
          xpk=xp(k+3)
          l=(3*pairs(i,2))-1
          xei=xe(l+1)
          xej=xe(l+2)
          xek=xe(l+3)
          d=dsqrt((xpi-xei)**2+(xpj-xej)**2+(xpk-xek)**2)
          if (d.lt.rmstol) then
            npafter=npafter+1
            nptmp(npafter)=np(i)
          else
            nptrial=nptrial-1
          endif
        enddo
        if (nptrial.lt.npal) return
        do i = 1, nvnew 
          lsup(i)=' '
        enddo
        do i = 1, npafter
          np(i)=nptmp(i) 
          lsup(np(i))='*'
        enddo
        npal=npafter

        return
        end
c
c===============================================================================
c
      subroutine extendali(nvnew,lsup,nsup1,nsup2,iv1,iv2,
     &                     assigned1,assigned2,npairs,pairs,xp,xe)

c-------------------------------------------------------------------------------
c --- Extends the alignment provided by maxsub.
c --- It has not been checked carefully!!!.
c-------------------------------------------------------------------------------

      implicit REAL*8 (a-h,o-z)

      parameter  (maxres = 3000)
      parameter  (maxlen = 2*maxres)
      parameter  (rmstol = 4.0d0)
      parameter  (angmax = 60.0d0)

      dimension    xp(3*maxres),xe(3*maxres),ai(3),aj(3)
      character*1  lsup(maxlen)
      logical      assigned1(maxres),assigned2(maxres)
      integer      iv1(maxlen),iv2(maxlen),pairs(maxres,2)

      external     angle

c-------------------------------------------------------------------------------
c --- pairs from the initial local backbone alignment
c-------------------------------------------------------------------------------

      npairs = 0
      do i = 1, nvnew
        assigned1(i)=.false.
        assigned2(i)=.false.
      enddo
      do i = 1, nvnew
        if (lsup(i).eq.'*') then
          npairs=npairs+1
          pairs(npairs,1)=iv1(i)
          pairs(npairs,2)=iv2(i)
          assigned1(iv1(i))=.true.
          assigned2(iv2(i))=.true.
        endif
      enddo
      npal=npairs

c-------------------------------------------------------------------------------
c --- extension of the original alignment
c-------------------------------------------------------------------------------

      m21=nsup1-2
      m22=nsup2-2
      ii=0
      do i = 1, nsup1 
        dmin=1000.0
        if (.not.assigned1(i)) then
          xix=xp(ii+1)
          xiy=xp(ii+2)
          xiz=xp(ii+3)
          jj=0
          do j = 1, nsup2  
            if (.not.assigned2(j)) then
              xjx=xe(jj+1)
              xjy=xe(jj+2)
              xjz=xe(jj+3)
              d=dsqrt((xix-xjx)**2+(xiy-xjy)**2+(xiz-xjz)**2)
              if (d.le.dmin) then
                if((i.ge.2.or.i.le.m21).and.(j.ge.2.or.j.le.m22))then
                  ai(1) = xp(ii+4)-xp(ii-2)
                  ai(2) = xp(ii+5)-xp(ii-1)
                  ai(3) = xp(ii+6)-xp(ii  )
                  aj(1) = xe(jj+4)-xe(jj-2)
                  aj(2) = xe(jj+5)-xe(jj-1)
                  aj(3) = xe(jj+6)-xe(jj  )
                  if (angle(ai,aj).lt.angmax) then
                    jas=j
                    dmin=d
                  endif
                endif
              endif
            endif
            jj=jj+3
          enddo
          if (dmin.le.rmstol+0.30) then
            do k = 1, npal-1
              if (pairs(k,1).lt.i.and.pairs(k+1,1).gt.i) then
                if(pairs(k,2).gt.jas+15.or.pairs(k+1,2).lt.jas-15)then
                  goto 115
                endif
              endif
            enddo
            if(i.lt.pairs(1,1).and.jas-15.gt.pairs(1,2))goto 115
            if(i.gt.pairs(npal,1).and.jas.lt.pairs(npal,2)-15)goto 115
            npairs=npairs+1
            assigned1(i)  =.true.
            assigned2(jas)=.true.
            pairs(npairs,1)=i
            pairs(npairs,2)=jas
          endif
        endif
 115    continue
        ii=ii+3
      enddo

c-------------------------------------------------------------------------------
c --- end of the subroutine
c-------------------------------------------------------------------------------

      return
      end

c
c=========================================================================
c
      subroutine putpdb(nat,x,nunit,rnam,ir)

c-------------------------------------------------------------------------
c --- Output a selected fraction of atoms in a pdb file. --
c-------------------------------------------------------------------------

      implicit REAL*8 (a-h,o-z)

      parameter  (maxres = 3000)

      dimension    x(3*maxres)
      integer      ir (maxres)
      character*3  rnam(maxres)

      j=0
      do i = 1, nat
        write(nunit,500) i,rnam(i),ir(i),x(j+1),x(j+2),x(j+3)
        j=j+3
      enddo
      write(nunit,'(a)') 'TER'
  500 format('ATOM ',i6,1x,' CA ',1x,a3,2x,i4,4x,3f8.3)

      return
      end
c
c========================================================================
c
      FUNCTION erfcc(x)

c-------------------------------------------------------------------------
c  (C) Copr. 1986-92 Numerical Recipes Software
c-------------------------------------------------------------------------

      real*8 erfcc,x
      real*8 t,z
      z=abs(x)
      t=1./(1.+0.5*z)
      erfcc=t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+t*
     *(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+t*
     *(1.48851587+t*(-.82215223+t*.17087277)))))))))
      if (x.lt.0.) erfcc=2.-erfcc
      return
      END
c
c========================================================================

         SUBROUTINE COMAS(C,WT,NAT,XC,YC,ZC)
c
c-------------------------------------------------------------------------
c --- Calculate the center of geometry for the selected atoms ---
c-------------------------------------------------------------------------
c
      implicit REAL*8 (a-h,o-z)
c
      DIMENSION C(*),WT(*)
      DATA ZERO/0.0D+00/
c
      SUMX = ZERO
      SUMY = ZERO
      SUMZ = ZERO
      SUM  = ZERO
      I3 = 0
      DO 100 I = 1, NAT
        WTI = WT(I)
        SUMX = SUMX+C(I3+1)*WTI
        SUMY = SUMY+C(I3+2)*WTI
        SUMZ = SUMZ+C(I3+3)*WTI
        SUM = SUM+WTI
        I3 = I3+3
  100 CONTINUE
      SUM = 1.0D0/SUM
      XC = SUMX*SUM
      YC = SUMY*SUM
      ZC = SUMZ*SUM
c
      I3 = 0
      DO 120 I = 1, NAT
        C(I3+1) = C(I3+1)-XC
        C(I3+2) = C(I3+2)-YC
        C(I3+3) = C(I3+3)-ZC
        I3 = I3+3
  120 CONTINUE
c       WRITE(6,9008) XC,YC,ZC
c9008   FORMAT(/5X,'CENTER OF MASS:', 19x, 3F10.4)
      RETURN
      END
c

      
 
