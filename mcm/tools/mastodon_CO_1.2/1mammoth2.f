c --- Version 1.2  may 19 2006
c --- last revision by CEMS
c=========================================================================
c ---
c ---                          ~ TURBO MAMMOTH ~
c ---
c ---         MAtching Molecular Models Obtained from THeory 
c ---               Copyright (C) Angel R. Ortiz, 2001
c ---        Copyright (C) Accelerated Algorithm By Charlie E.M. Strauss,2001
c_______________________________________________________________________________
c
c     This software is copyrighted by Angel R. Ortiz and by Charlie E.M. Strauss.
c     The following terms apply to all files associated with the software 
c     unless explicitly disclaimed in individual files.
c
c     The authors hereby grants permission to use, copy and modify this software
c     and its documentation for any non-commercial purpose, provided that existing 
c     copyright notices are retained in all copies and that this notice is included
c     verbatim in any distributions, and that the source is cited in publications.
c     Redistribution as a part of a commercial distribution is considered commerical use.
c     No written agreement, license, or royalty fee is required for any of the 
c     authorized uses.
c     
c     
c
c     IN NO EVENT SHALL THE AUTHORS OR DISTRIBUTORS BE LIABLE TO ANY PARTY
c     FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
c     ARISING OUT OF THE USE OF THIS SOFTWARE, ITS DOCUMENTATION, OR ANY
c     DERIVATIVES THEREOF, EVEN IF THE AUTHOR HAS BEEN ADVISED OF THE
c     POSSIBILITY OF SUCH DAMAGE.
c
c     THE AUTHORS AND DISTRIBUTORS SPECIFICALLY DISCLAIM ANY WARRANTIES,
c     INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY,
c     FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT.  THIS SOFTWARE
c     IS PROVIDED ON AN "AS IS" BASIS, AND THE AUTHOR AND DISTRIBUTORS HAVE
c     NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
c     MODIFICATIONS.
c
c_______________________________________________________________________________
c 
c ---
c --- Algorithm:
c --- =========
c ---
c --- The method used by this program is as follows:
c ---
c ---  1.- From the Calpha trace, compute the unit-vector U-RMS between
c --- all pairs of heptapeptides of both model and experimental structure.
c --- The U-RMS is described in: Kedem, Chew & Elber (1999) Proteins
c --- 37(4):554-64, and in Chew, Huttenlocher, Kedem & Kleinberg (1999)
c --- J.Comp.Biol. 6, 313-325. This is a measure sensitive to the local 
c --- structure.
c ---
c ---  2.- Use the matrix derived in step 1 to find and alignment of
c --- local structures that maximizes the local similarity of both the
c --- model and the experimental structure. For that, use a global
c --- alignment method with zero end gaps, as described in
c --- Needleman & Wunsch (1970) J.Mol.Biol. 48, 443-453. 
c ---
c ---  3.- Find the maximum subset of similar local structures that have
c --- their corresponding Calphas close in cartesian space. "Close" is
c --- considered here as a distance less or equal than 4.0 A. The method
c --- to find this subset is an accelerated variant of the MaxSub algorithm 
c --- from the Fischer group: Siew, Elofsson, Rychlewski & Fischer (2000)
c --- Bioinformatics, in press. 
c ---
c ---  4.- Obtain the probability of obtaining the given proportion of
c --- aligned residues (with respect to the shortest protein model) by
c --- chance. This metric (E-value) is then used as the final score 
c --- (or the corresponding Z-score, both are equivalent for gaussian
c --- distributions, however the Z-score is a more manegable index).
c --- In order to obtain this value, an approach similar to that of
c --- Levitt & Gerstein (1998) PNAS 95, 5913 is used, as described in
c --- Abagyan & Batalov (1997) J.Mol.Biol. 273, 355-368. The E-value
c --- estimation is based on extreme-value fitting. In a test set
c --- with the SCOP database, it shows rather good performance.
c ---
c ---
c --- Compilation: 
c --- ===========
c --- 
c ---
c --- 1. Under Linux 6.0  with Absoft:
c ---        
c --- f77 -O -lU77 -N109 -W -s *.f
c or  f77 -O -N18 -W -N34 -N109 -N1 -s -lU77 *.f 
c --- 
c --- 2. Under Mac OS X(10.1) with Absoft
c --- 
c --- f77  -O -B18 -W -N34 -N109 -N1 -s -lU77 -lfio -lf90math -lf77math
c --- 
c --- 3. Under  gnu g77
c ---
c --- g77 -O5 -ffixed-line-length-132 *.f
c ---
c --- 4. Under IRIX 6.5 with r10k:
c ---
c --- f77 -IPA -Ofast=ip28 -LNO -n32 -mips4 -r10000 -o *.f
c ---
c --- 5. Under ALPHA with compaq Fort
c ---
c --- fort -O5 -extend_source  *.f
c ---
c --- things consider on other compliers:
c ---     i) extended fortran line lengths
c ---    ii) static variables (most fortran's default except absoft)
c ---   iii) longjumps/datasizes (required on some compilers for 64 bit addresses)
c ---
c --- Authors & Contact information:
c --- =============================
c --- 
c --- Authors:
c ---
c ---   -Angel R. Ortiz coded the initial version of the program.
c ---   -Osvaldo Olmea contributed to several features of the program
c ---    and also corrected some of the bugs...
c ---   -Charlie Strauss accelerated the algorithm, added features and fixed bugs.
c ---
c --- Contact information:
c ---
c ---   Angel Ramirez Ortiz             |  Phone: 
c ---   Assistant Professor             |  Fax:   
c ---   Dep. Physiology & Biophysics    | 
c ---   Mount Sinai School of Medicine  |  e-mail: ortiz@inka.mssm.edu
c ---   One Gustave Levy Pl., Box 1218  |  http://www.mssm.edu
c ---   New York, NY 10029  
c ---
c ---   Charlie E. M. Strauss           | Phone: (505) 665 4838
c ---   Techincal Staff
c ---   Los Alamos National Laboratory
c ---   Los Alamos NM 87545
c ---   e-mail: cems@lanl.gov   or cems at browndogs dot org
c ---   

c ---    CHANGES
c ---    added Contact order calculation CEMS  2006
c ---    may 18 2006 changed cut-off on e-value  CEMS
c ---        previously all e-values smaller than  2.65E-15 were truncated.
c ---        the practical effect was that comparing any structure to itself
c ---        with a length > 258 residues returned the same e-value.
c ---        now use a taylor expansion to compute values below this point. 
c ---   
c ---


c ---
c ---             | 
c ---   
c ---
c ---                          
c ---
c=========================================================================
c
      implicit REAL*8 (a-h,o-z)

      parameter  (maxatom= 3000)   ! need to reconcile maxatom and maxres, a bit confused in use
      parameter  (maxres = 3000)
      parameter  (maxfrag=  100)
      parameter  (maxlen = 2*maxres)
      parameter  (max_pdb_list1=1000)   ! yeah I know its assymetric
      parameter  (max_pdb_list2=10000)
      parameter  (iEOF = -1)

c --- file management

      character*10 pname1,pname2
      character*340 pdb1,pdb2,outfile
      character*80 string,buffer
      logical      pdb1_exists,pdb2_exists,report
      
      character*256 pdb_path1,pdb_path2
      character*80  pdb_list1(max_pdb_list1),pdb_list2(max_pdb_list2)
      integer       pdb_list_len1,pdb_list_len2
 
      
      logical       generate_pdb   

c --- protein descriptions

      integer      iv1(maxlen),iv2(maxlen)
      integer      nsec1(maxlen),nsec2(maxlen)
      integer      ica1(maxres),ica2(maxres)
      integer      ires1(maxres),ires2(maxres)
      integer      isec1(maxres),isec2(maxres)
      integer      ised1(-2:maxres),ised2(-2:maxres)

      character*1  sec1(maxres), sec2(maxres)
      character*1  aseq1(maxlen),aseq2(maxlen)
      character*1  rseq1(maxlen),rseq2(maxlen)
      character*1  match(maxlen),lsup(maxlen)
      character*1  aaname1(0:20)
      character*3  aaname3(0:20),res

      dimension    x1(maxatom,3), x2(maxatom,3)
      dimension    xca1(maxres,3),xca2(maxres,3)
      
            
      real*8 rms_o1,psi_o1,rms_o2,psi_o2
      real*8 evalue,scorel,zscore
      integer norm_o1,nsup_o1,norm_o2,nali_o2
      integer nssmatch_o
      logical terse
      integer      iatom1, resnum1, iatom2, resnum2
      character*4  atnam1, resnam1, atnam2, resnam2

      common /pdbdata/ iatom1(maxatom), resnum1(maxatom),
     &                 atnam1(maxatom), resnam1(maxatom),
     &                 iatom2(maxatom), resnum2(maxatom),
     &                 atnam2(maxatom), resnam2(maxatom)

      common /vectors/ v1(maxres,3),v2(maxres,3),smapori(maxres,maxres)


c -- contact order
	real*8 co_f1,co_f2
	integer COMPUTE_CO

c --- timing stuff
      character*80  part
      parameter        (maxcall=100)
      common /timings/ ftme(maxcall),part(maxcall)


      data aaname3 /'ALA','CYS','ASP','GLU','PHE','GLY','HIS',
     &              'ILE','LYS','LEU','MET','ASN','PRO','GLN',
     &              'ARG','SER','THR','VAL','TRP','TYR','XXX'/

      data aaname1 /'A','C','D','E','F','G','H','I','K','L','M',
     &              'N','P','Q','R','S','T','V','W','Y','X'/

     
     

c-------------------------------------------------------------------------
c --- Read comand line arguments
c --- pdb1 => Predicted conformation
c --- pdb2 => Experimental conformation
c-------------------------------------------------------------------------

      
c ------------------------------------------------------------------------
c     parse args
c ------------------------------------------------------------------------  
      ! defaults   
      generate_pdb = .true.
      z_thresh = -1.0
      terse=.false.
      outfile = 'STDOUT'
	  COMPUTE_CO = 0
      
      ! read args
      narg = iargc()
      if (narg.lt.4) call help
      do i = 1, narg, 2
         call getarg(i,string)
         if (string(1:1) .eq. '-') then
            if (string(2:2).eq.'p') then
               call getarg(i+1,buffer)
               read(buffer,'(a)') pdb1
               
            else if (string(2:2).eq.'e') then
               call getarg(i+1,buffer)
               read(buffer,'(a)') pdb2
               
            else if (string(2:2).eq.'o') then
               call getarg(i+1,buffer)
               read(buffer,'(a)') outfile
            else if (string(2:2).eq.'t') then ! output threshold cutoff
               call getarg(i+1,buffer)
               read(buffer,*) z_thresh
            else if (string(2:2).eq.'v') then ! verbose outupt
               call getarg(i+1,buffer)
               if (buffer(1:1).eq.'0'.or.buffer(1:1).eq.'n') terse=.TRUE.
            else if (string(2:2).eq.'r') then ! dont make pdbs
			    call getarg(i+1,buffer)
               if (buffer(1:1).eq.'0'.or.buffer(1:1).eq.'n') generate_pdb =.false.
			else if (string(2:2).eq.'c') then 
			     call getarg(i+1,buffer)
				 if (buffer(1:1).eq.'1'.or.buffer(1:1).eq.'y') COMPUTE_CO = 1
				 write (0,*) 'COMPUTE_CO FLAG ',buffer(1:1),COMPUTE_CO
			else	 
               call help
            endif
         else
            call help
         endif
      enddo
      write (0,*) 'COMPUTE_CO FLAG ',COMPUTE_CO
      if (terse) zthresh =1000
      if (terse) generate_pdb =.false.
      
c     ------------------------------------------------------------------------          
c     process args
c     ------------------------------------------------------------------------    
      
      
      
      
c     open output file
      if (outfile(1:len_trim(outfile)) .ne. 'STDOUT') then
         open(unit=6,file=outfile)
      endif
c     validate existence of input files, load lists with names and paths
      call getpdb_names(max_pdb_list1,300,pdb1,pdb_path1,pdb_list1,pdb_list_len1)
      write(0,*) " -p MAMMOTH List entries =", pdb_list_len1
      
      call getpdb_names(max_pdb_list2,300,pdb2,pdb_path2,pdb_list2,pdb_list_len2)
      write(0,*) " -e MAMMOTH List entries =", pdb_list_len2
      
      if (pdb_list_len1.gt.1 .or. pdb_list_len2.gt.1) generate_pdb=.false.
      
c     write header
      write(6,*) "Predicted path:   ",pdb_path1(1: lchar(pdb_path1) )
      write(6,*) "Experimental path:",pdb_path2(1: lchar(pdb_path2) )
      
      write(6,'(1x,79a,/)') ('_',i=1,79)
      write(6,'(/,22x,''T U R B O     M  A  M  M  O  T  H'',/)')
      write(6,'(/,15x,''MAtching Molecular Models Obtained '',
     &     ''from THeory'',/)')
      write(6,'(1x,79a,/)') ('_',i=1,79)
      write(6,*) ' '
      
      if (terse) then
         write(6,24) 'p ','e ','Zscore   ','-ln(E)   ','Evalue   ',
     &        'score    ','#e','#p','nsup','nss','psi1    ','psi2    '
	 
	 if (COMPUTE_CO .eq. 1) then
		write(6,"(2(a7,1x),$)")  'CO_p','CO_e'
	 endif
	 write(6,*) ! carriage return
	 
      endif
 24   format(1x,2(a5,1x),1x,4(a10,1x),4(a4,2x),2(a9,1x),$)
      
      
c-----------------------------------------------------------------------
c     Main Loop
c-----------------------------------------------------------------------
      do item2=1,pdb_list_len2
         pdb2 = pdb_path2(1:lchar(pdb_path2))//pdb_list2(item2)
         pdb2_exists = .false.
         inquire(file=pdb2(1:lchar(pdb2)),exist=pdb2_exists)
         
         if (.not.pdb2_exists) then
            write(0,*) "warning: the -e pdb file does not exist '",pdb2(1:lchar(pdb2)),"'"
            goto 723
         endif
         open(unit=30,file=pdb2,status='old')          
         call getpdb(nsel2,x2,30,iatom2,resnum2,resnam2,atnam2,nc2,ica2)
         close(30)
         
         do 15 i = 1, nc2
            ires2(i) = resnum2(ica2(i))
            res      = resnam2(ica2(i))
            do j = 0, 19
               if (res .eq. aaname3(j)) then
                  nsec2(i) = j
                  goto 15
               endif
            enddo
            write(6,*) 'WARNING: Residue not found in -e : ',pdb2(1:lchar(pdb2)), res
 15      continue
                
         do item1=1,pdb_list_len1 
 721        continue
            pdb1 = pdb_path1(1:lchar(pdb_path1))//pdb_list1(item1)
                                ! remove file from list if it does not exist.
            if (item2.eq.1) then ! only need to test if files exist on first pass of outer loop
               pdb1_exists = .false.
               inquire(file=pdb1,exist=pdb1_exists)
               if (.not. pdb1_exists) then
                  write(0,*) "warning: -p pdb file does not exist ",pdb1(1:lchar(pdb1))
                                ! remove slacker from the list
                  pdb_list1(item1) = pdb_list1(pdb_list_len1) ! overwrite slacker with last file in list
                  pdb_list_len1 = pdb_list_len1-1 ! contract list
                  if (pdb_list_len1.lt.1) then ! list now empty?
                     write(0,*) "No valid files found in pdb1 list path", pdb_path1
                     call my_halt()
                  endif
                  goto 721      ! reprocess this item            ! if not retry
               endif
            endif
            open(unit=20,file=pdb1,status='old')
    
c-----------------------------------------------------------------------
c     Start the clock
c-----------------------------------------------------------------------
            
            ntime = 0
            ttime = 0.0d0
            call setime
c-------------------------------------------------------------------------
c --- Read the pdb file to be superimposed
c-------------------------------------------------------------------------
            call getpdb(nsel1,x1,20,iatom1,resnum1,resnam1,atnam1,nc1,ica1)
            close(20)
 
C$$$  write(6,'(/,a,/)')     '==> PREDICTION: '
C$$$  write(6,'(5x,a,a)')    'Filename: ', pdb1(1:lchar(pdb1))
C$$$  write(6,'(5x,a,i5,/)') 'Number of residues: ', nc1
            do 10 i = 1, nc1
               ires1(i) = resnum1(ica1(i))
               res      = resnam1(ica1(i))
               do j = 0, 19
                  if (res .eq. aaname3(j)) then
                     nsec1(i) = j
c     write(6,*) i, nsec1(i), res, ires1(i)
                     goto 10
                  endif
               enddo
               write(6,*) 'WARNING: Residue not found: ', res
 10         continue
            
            
C$$$  
C$$$  
C$$$  write(6,'(/,a,/)')     '==> EXPERIMENT: '
C$$$  write(6,'(5x,a,a)')    'Filename: ', pdb2(1:lchar(pdb2))
C$$$  write(6,'(5x,a,i5,/)') 'Number of residues: ', nc2
C$$$  
            
      
 
C$$$ 
C$$$       write(6,*) ' '
C$$$       write(6,*) '-----------------------------'
C$$$       write(6,*) ' Structural Alignment Scores '
C$$$       write(6,*) '-----------------------------'
C$$$       write(6,*) ' '

            call storetime('Initialization:       ',ttime,ntime)
c-------------------------------------------------------------------------------
c --- Select calfas and assign secondary structure based on r14.
c --- This part needs to be simplified in only one subroutine
c --- secondary structure assignment is not really needed, it's only
c --- used in th eoutput, thus it can be avoided during database
c --- searching
c-------------------------------------------------------------------------------

            call setime

c --- Predicted coordinates. Get Ca coordinates and Ca-Ca vectors.
c --- Assign also the secondary structure

            do i = 1, nc1
               ic1 = ica1(i)
               do k = 1, 3
                  xca1(i,k) = x1(ic1,k)
               enddo
c     write(6,*) i, (xca1(i,k),k=1,3)
            enddo
            
            call assignss(nc1,xca1,isec1,sec1,ised1)
c     write(6,'(50a1)') (sec1(i),i=1,nc1)
            
            nv1 = nc1-1
            do i = 1, nv1
               v1mod = 0.0d0
               do k = 1, 3
                  v1(i,k) = xca1(i+1,k)-xca1(i,k)
                  v1mod = v1mod + v1(i,k)*v1(i,k)
               enddo
               v1mod = sqrt(v1mod)
               do k = 1, 3
                  v1(i,k) = v1(i,k)/v1mod
               enddo
            enddo
            
c --- Experimental coordinates. Get Ca coordinates and Ca-Ca vectors.
c --- Assign also the secondary structure
            
            do i = 1, nc2
               ic2 = ica2(i)
               do k = 1, 3
                  xca2(i,k) = x2(ic2,k)
               enddo
            enddo
            
            call assignss(nc2,xca2,isec2,sec2,ised2)
            
c     write(6,'(50a1)') (sec2(i),i=1,nc2)
            nv2 = nc2-1
            do i = 1, nv2
               v2mod = 0.0d0
               do k = 1, 3
                  v2(i,k) = xca2(i+1,k)-xca2(i,k)
                  v2mod   = v2mod + v2(i,k)*v2(i,k)
               enddo
               v2mod = dsqrt(v2mod)
               do k = 1, 3
                  v2(i,k) = v2(i,k)/v2mod
               enddo
c     write(6,'(i5,5x,3f8.3)') i, (v2(i,k),k=1,3)
            enddo
            
            call storetime('Secondary structure assignment:',ttime,ntime)
            call setime
c-------------------------------------------------------------------------------
c --- Vector alignment using dynamic programming
c-------------------------------------------------------------------------------
       
            
            report = .true.
            
            call alignv(report,nvnew,nv1,nv2,iv1,iv2,sc,pgaps,match)
c     write(6,'(50i2)') (iv1(i),i=1,nvnew)
c     write(6,*)
c     write(6,'(50i2)') (iv2(i),i=1,nvnew)
            
            call storetime('Secondary structure alignment: ',ttime,ntime)
            
c-------------------------------------------------------------------------------
c     --- Now find the maximum subset of residues aligning below 4.0 A (rmstol)
c-------------------------------------------------------------------------------
            
            call setime
            
c     do i = 1, nc1
c     write(6,*) i, aaname3(nsec1(i))
c     enddo
            norm = min(nc1,nc2)
            
            call maxsub(ires1,ires2,nsec1,nsec2,nvnew,iv1,iv2,
     &           xca1,xca2,norm,match,lsup,sec1,sec2,z_thresh,generate_pdb,
     &           zscore,evalue,scorel,nssmatch_o,
     &           norm_o1,nsup_o1,rms_o1,psi_o1,
     &           norm_o2,nali_o2,rms_o2,psi_o2,
     &					co_f1,co_f2,COMPUTE_CO)   
            call storetime('Tertiary structure matching:   ',ttime,ntime)
c-------------------------------------------------------------------------------
c     --- Reports structure-based  alignment
c-------------------------------------------------------------------------------
            call setime
            if (terse) then
               
               write (6,23) item1,item2
     &              ,zscore,-log(evalue),evalue,scorel
     &              ,nc1,nc2,nsup_o1,nssmatch_o
     &              ,psi_o1,psi_o2
	 
 23            format(':',2(i5,1x),1x,4(g10.4,1x),4(i4,2x),2(g9.3,2x),$)	 
	               if (COMPUTE_CO .eq. 1) then
					write(6,"(2(f7.1,x),$)") co_f2,co_f1
				   endif
				   
			write(6,"(2x,2(a,2x))") 
     &              ,pdb_list1(item1)(1:lchar(pdb_list1(item1)))
     &              ,pdb_list2(item2)(1:lchar(pdb_list2(item2)))


             
			  		
            else
      
               write(6,*) '-------------------'
               write(6,*) ' Input information '
               write(6,*) '-------------------'
               write(6,*) ' '
               
               write(6,'(/,a,/)')     '==> PREDICTION: '
               write(6,'(5x,a,a)')    'Filename: ', pdb1(1:lchar(pdb1))
               write(6,'(5x,a,i5,/)') 'Number of residues: ', nc1
               write(6,'(/,a,/)')     '==> EXPERIMENT: '
               write(6,'(5x,a,a)')    'Filename: ', pdb2(1:lchar(pdb2))
               write(6,'(5x,a,i5,/)') 'Number of residues: ', nc2
               
               write(6,*) ' '
               write(6,*) '-----------------------------'
               write(6,*) ' Structural Alignment Scores '
               write(6,*) '-----------------------------'
               write(6,*) ' '
               
               
               write(6,*)
               write(6,21)' PSI(ini)= ',psi_o1,'NALI= ',nsup_o1,'NORM= ',norm_o1,
     &              'RMS= ',rms_o1
	           write(6,"(2x,a,i3)") 'NSS= ',nssmatch_o
	         
	 
               write(6,21)' PSI(end)= ',psi_o2,'NALI= ',nali_o2,'NORM= ',norm_o2,
     &              'RMS= ',rms_o2
               
                  if (COMPUTE_CO .eq. 1) then
				write(6,"(2(2x,a,f6.2),$)") 'CO_p=',co_f2, 'CO_e=',co_f1
			endif
			   write(6,*)
			   
               write(6,21)' Sstr(LG)= ',scorel,'NALI= ',nali_o2,'NORM= ',norm_o2,
     &              'RMS= ',rms_o2
               
               
               write(6,*)
               write(6,"(a,1x,g17.8)") 'E-value= ', evalue
               write(6,*)
               write(6,"(a,1x,g17.8,1x,a,1x,g17.8)") 'Z-score= ', zscore, '   -ln(E)= ', -log(evalue)
               
 21            format(a,f7.2,2x,2(a,i3,2x),a,f7.2,$)   
               if (zscore.ge.z_thresh) then
                  do i = 1, nvnew
                     if (iv1(i).eq.0) then
                        aseq1(i) = '.'
                        rseq1(i) = '.'
                     else
                        aseq1(i) = sec1(iv1(i))
                        rseq1(i) = aaname1(nsec1(iv1(i)))
                     endif
                     if (iv2(i).eq.0) then
                        aseq2(i) = '.'
                        rseq2(i) = '.'
                     else
                        aseq2(i) = sec2(iv2(i))
                        rseq2(i) = aaname1(nsec2(iv2(i)))
                     endif
                  enddo
                  
                  write(6,*) ' '
                  write(6,*) '----------------------------'
                  write(6,*) ' Final Structural Alignment '
                  write(6,*) '----------------------------'
                  write(6,*) ' '
                  ioff = 0
 70               ibeg = ioff + 1
                  if (nvnew .ge. ibeg+49) then
                     iend = ibeg+49
                  else
                     iend = nvnew
                  endif
c     pname1 = pdb1(1:lchar(pdb1))
c     pname2 = pdb2(1:lchar(pdb2))
                  pname1 = 'Prediction'
                  pname2 = 'Experiment'
                  write(6,301)        (lsup(i),i=ibeg,iend)
                  write(6,300) pname1,(rseq1(i),i=ibeg,iend)
                  write(6,300) pname1,(aseq1(i),i=ibeg,iend)
                  write(6,300) ' ',(match(i),i=ibeg,iend)
                  write(6,300) pname2,(aseq2(i),i=ibeg,iend)
                  write(6,300) pname2,(rseq2(i),i=ibeg,iend)
                  write(6,301)        (lsup(i),i=ibeg,iend)
                  write(6,*) ' '
                  ioff = iend
                  if (iend .lt. nvnew) goto 70
c     write(*,500) sc, sc-smax, smax
               endif
               
               call storetime('Text Output',ttime,ntime)
c-------------------------------------------------------------------------------
c     --- End of the program
c-------------------------------------------------------------------------------
               
               call printime(ntime)
        
            if (pdb_list_len1.ne.1 .or. pdb_list_len2.ne.1) then
               write(6,'(a)') "<MAMMOTH> END_OF_RECORD"
            endif
          endif ! end terse mode
         enddo
 723     continue
      enddo
      write(6,'(a)')  "<MAMMOTH> NORMAL_EXIT"      
      close(6)
      
 300  format(a10,1x,5(10a1,1x))
 301  format(11x,5(10a1,1x))
      
      end
c     
c==========================================================================
c
      subroutine assignss(nres,xyz,isec,sec,ised)

      implicit   REAL*8 (a-h,o-z)

      parameter  (maxres = 3000)

      character*1 sec(maxres)
      integer     isec(maxres),ised(-2:maxres)
      dimension   xyz(maxres,3)
c     logical     natlbin(maxres)
c     integer     natlbin(maxres)

      do i = 1, nres
        isec(i) = 5
      enddo
      call getl3g(nres,xyz,isec,ised)


      do i = 1, nres
        if (isec(i).eq.2) then
          sec(i) = 'H'
        else if (isec(i).eq.4) then
          sec(i) = 'S'
        else
          sec(i) = '-'
        endif
      enddo

      return
      end
c
c========================================================================
c
      subroutine getl3g(nres,xyz,natlbin,ised)
C
c-------------------------------------------------------------------------
c
C       secassign. Assigns the dominant secondary secondary elements
c       based on the r14 criterion:
c
C       *********************************************************
C       BASED ON A SMOOTHED DISTRIBTION USES THE NEW 94 database*
C       *               ENER14G.F       3/2/93                  *
C       *       CALCULATES THE R14 DISTRIBUTION FUNCTION IN THE *
C       *       DATABASE                                        *
C       *       FOR THE 311 LATTICE BUT WITH A SMOOTHED         *
C       *       DISTRIBUTION                                    *
C       *                                                       *
C       *********************************************************
c
c-------------------------------------------------------------------------
c
        implicit REAL*8 (a-h,o-z)

        parameter (maxres = 3000)
        PARAMETER(ibins=6)
        DIMENSION xa(3,maxres),pxyz(3,maxres),xyz(maxres,3)
c       logical*1 natlbin(maxres)
        integer   natlbin(maxres),ised(-2:maxres)
        common/is/ibb(-87:91)
c
c --- Patch to keep he data structure of the r14s and getl3g routines,
c --- then transform real space to lattice units (conv. factor = 1/1.22)
c
        do i = 1, nres
          do j = 1, 3
            pxyz(j,i) = xyz(i,j)
            xa(j,i)   = pxyz(j,i) * 0.8197
          enddo
        enddo
c
c --- Generate bin states
c
        do i=-86,91
          if(i .lt. -56)                        ibb(i)=4
          if(i .ge. -56 .and. i .lt. -25)       ibb(i)=3
          if(i .ge. -25 .and. i .lt.   0)       ibb(i)=1
          if(i .ge.   0 .and. i .lt.  26)       ibb(i)=2
          if(i .ge.  26 .and. i .lt.  56)       ibb(i)=3
          if(i .ge.  56)                        ibb(i)=4
        end do
c
c --- Assign r1-4
c
      call r14s(nres,xa,natlbin,ised)
C
      return
      end
C
c=========================================================================
c
        SUBROUTINE R14S(nres,xa,natlbin,ised)
C
c-------------------------------------------------------------------------
C       ALL  DIMENSIONS ARE IN LATTICE UNITS IN THIS PROGRAM
c------------------------------------------------------------------------
c
        implicit REAL*8 (a-h,o-z)

        LOGICAL QSIDE
        parameter(maxres = 3000)
        parameter(maxsec =  200)
        integer sed(-2:maxres),ised(-2:maxres)
        DIMENSION QSIDE(maxres)
c       logical*1 natlbin(maxres)
        integer   natlbin(maxres)
        common/is/ibb(-87:91)
        dimension mbegin(0:maxres)
        dimension mend(0:maxres),lbegin(maxsec),lend(0:maxsec)
        dimension istate(maxsec),lstate(maxsec)
        dimension iconf(0:maxres)
        dimension y(3,-7:maxres)
        dimension rc(3),xa(3,maxres),xb(3,-10:maxres)
        dimension rad(maxres),ased(maxres)
        dimension awt(0:5)
 
c-------------------------------------------------------------------------
c  Define a few parameters and initialize
c-------------------------------------------------------------------------

        awt(5)= 1./243.
        awt(4)= 5./243.
        awt(3)=15./243.
        awt(2)=30./243.
        awt(1)=45./243.
        awt(0)=51./243

        DO I=nres+1,maxres
          ised(i)=0
          sed(i)=0
          qside(i)=.false.
        end do

        DO I=1,NRES
          ised(i) =0
          sed(i)  =0
          qside(i)=.true.
        end do

        CA=3.785/1.22
        CA2=CA*CA
 
c-------------------------------------------------------------------------
c       cos(180-valence_angle) is in the range [-.30 0.901] including.
c       which corresponds to the range of val_ang 72.5 to 154
c-------------------------------------------------------------------------
        DO I=2,NRES-1
          R21=0.
          R22=0.
          DP=0.
          DO J=1,3
             R22=R22+(XA(J,I+1)-XA(J,I))**2
            R21=R21+(XA(J,I)-XA(J,I-1))**2
            DP=DP+(XA(J,I+1)-XA(J,I))*(XA(J,I)-XA(J,I-1))
          END DO
          IF(ABS(R21/CA2)-1. .GT. 0.2)THEN
            QSIDE(I)=.FALSE.
            GO TO 101
          END IF

C       AAA IS (180-VALANCE ANGLE)
          AAA=DP/SQRT(R22*R21)
          IF(AAA.LT.-0.30 .OR. AAA.GT.0.901) THEN
            QSIDE(I)=.FALSE.
          ELSE
            QSIDE(I)=.TRUE.
          END IF
101     CONTINUE
        END DO
 
c-------------------------------------------------------------------------
c       COMPUTE R14 AND CHECK THAT R14 IS NOT UNREASONABLE
c-------------------------------------------------------------------------
        do 140 i=2,nres-2
          ip1=i+1
          if(qside(i).and. qside(ip1))then
            im1=i-1
            ip2=i+2
            r14=(xa(1,im1)-xa(1,ip2))**2 +
     &          (xa(2,im1)-xa(2,ip2))**2 +
     &          (xa(3,im1)-xa(3,ip2))**2
 
            IR14=INT(R14+0.5)
            IF(IR14 .GT. 91)GO TO 12
            WX1=XA(1,I)-XA(1,IM1)
            WY1=XA(2,I)-XA(2,IM1)
            WZ1=XA(3,I)-XA(3,IM1)
            WX2=XA(1,IP1)-XA(1,I)
            WY2=XA(2,IP1)-XA(2,I)
            WZ2=XA(3,IP1)-XA(3,I)
            PX=WY1*WZ2-WY2*WZ1
            PY=WZ1*WX2-WZ2*WX1
            PZ=WX1*WY2-WX2*WY1

            WX3=XA(1,IP2)-XA(1,IP1)
            WY3=XA(2,IP2)-XA(2,IP1)
            WZ3=XA(3,IP2)-XA(3,IP1)
            IHAND=PX*WX3+PY*WY3+PZ*WZ3
            IF(IHAND .LT. 0) IR14=-IR14
            IF(Ir14 .lt. -86)IR14=-86
            I14=ibb(IR14)
            iconf(i)=ir14
            SED(i)=I14
            ised(i)=i14
12          CONTINUE
          end if
140     continue
c
c-------------------------------------------------------------------------
c       set up the location of the secondary elements:
c       renumber if there are less than 3 residues per beta state
c-------------------------------------------------------------------------
c
c       create average tube:
c
        do i=1,nres
          do j=1,3
            xb(j,i)=xa(j,i)
          end do
        end do
        do i=-4,0
          do j=1,3
            xb(j,i)=xb(j,1)
          enddo
        end do
        do i=1,5
          do j=1,3
            xb(j,i+nres)=xb(j,nres)
          end do
        end do
        do i=1,nres
          do j=1,3
            y(j,i)=0.
            do iw=-5,5
              iwa=abs(iw)
              y(j,i)=y(j,i)+awt(iwa)*xb(j,i+iw)
            end do
          end do
        end do
        do i=1,nres
          do j=1,3
            xa(j,i)=y(j,i)
          end do
        end do

        do i=3,nres-3
          rad(i)=0.
          sum=0.
          anorm=0.
          do j=1,3
            anorm=anorm+(y(j,i+1)-y(j,i-1))**2
            rc(j)=y(j,i-2)+y(j,i+2)-2.*y(j,i)
            sum=sum+rc(j)**2/16.
          end do
          anorm=sqrt(anorm)
          sum=sum/(anorm+0.0001)
          rad(i)=sum
          if(sum .lt. .09)then
            sed(i)=0
          else
            sed(i)=1
          end if
        end do
22      format(i4,3f6.4)

        icount=0
C100    continue
        SED(0)=0
        SED(1)=0
        SED(2)=0
        sed(3)=0
        sed(nres-2)=0
        sed(nres-1)=0
        sed(nres)=0
        icount=icount+1
Cif(icount .gt. 3) goto 200
        NLOOP=0
        DO I=3,NRES-3
c       IF(SED(I).EQ.1 .and. SED(I-1) .ne. 1 .and. SED(i-2) .ne.1)THEN
        IF(SED(I).EQ.1 .and. SED(I-1) .eq. 0)then
                nloop=nloop+1
                istate(nloop)=1
                MBEGIN(nloop)=I
                mend(nloop)=I
c       ELSEIF(SED(i).eq.1.and.SED(I+1).ne.4.and.SED(I+2).ne.4)THEN
        ELSEIF(SED(i).eq.1 .and.SED(I+1).eq.0)then
        mend(nloop)=i
        END IF
        END DO
Cwrite(6,*)'no. of loops =', nloop
        do k=1,nloop
Cwrite(6,*)'begin=',mbegin(k),'  end=',mend(k)
        end do


        nelmt=0
        mend(0)=0
        mbegin(nloop+1)=nres
        do k=1,Nloop
        if(Mbegin(k)-mend(K-1).gt. 2)then
                nelmt=nelmt+1
                lbegin(nelmt)=mbegin(K)
                lend(nelmt)=mend(K)
                lstate(nelmt)=istate(k)
                IF(LEND(nelmt) .eq. nres-1)LEND(nelmt)=nres-2
        else
        lend(nelmt)=mend(k)
        end if
        end do
Cwrite(6,*)'the no. of secondary structural elements=',nelmt
        do k=1,nelmt
        ik=lstate(k)
Cwrite(6,*)k,' ',lbegin(k),' ',lend(k)
        end do
        lend(0)=0
        lbegin(nelmt+1)=nres-1
C
        do i = 1,nres
           ased(i) = 3.0
        enddo
C
        do 2 k=1,nelmt+1
        ib=lend(k-1)+1
        ie=lbegin(k)-1
        if (ib.eq.ie) goto 2
        a=0.
        is=0
        do  ii=ib+1,ie
        is=is+1
        do j=1,3
        a=a+(xa(j,ii)-xa(j,ii-1))**2
        end do
        end do
        a=sqrt(a)/is
        if(a .lt. 0.6)then
c         write(6,*)'element=',k,a,ib,ie,' helix'
          lstate(k)=2
        else if(a .lt. 1.1)then
c         write(6,*)'element=',k,a,ib,ie,' beta'
          lstate(k)=4
        else
c         write(6,*)'element=',k,a,ib,ie,'extended  beta/turn'
          lstate(k)=4
        end if
C
        do i = ib,ie
           ased(i) = lstate(k)
        enddo
  2     continue
        !end do
C
        ased(1)=ased(2)
        ased(2)=ased(3)
        ased(nres)=ased(nres-2)
        ased(nres-1)=ased(nres-2)
C
        do i = 1, nres
c          write(6,*) ased(i)
           natlbin(i) = int(ased(i))
        enddo
C
        return
C
        end
c
c===============================================================================
c
      subroutine alignv(report,nvnew,nv1,nv2,iv1,iv2,sc,pgaps,match)
 
c-------------------------------------------------------------------------------
c ---                                                                        ---
c ---              Alignment of two sets of vectors by the                   ---
c ---                     Needleman-Wunsch algorithm                         ---
c ---                                                                        ---
c ---  Needleman, S.B. & Wunsch, C.D. (1970) "A general method applicable    ---
c ---  to the search for similarities in the aminoacid sequence of two       ---
c ---  proteins". J.Mol.Biol. 48, 443-453.                                   ---
c ---                                                                        ---
c ---  This is a global alignment algorithm. No gap ends are used.           ---
c ---                                                                        ---
c ---                               =ARO-00= 
c-------------------------------------------------------------------------------
 
      implicit REAL*8 (a-h,o-z)

      parameter  (maxres = 3000)
      parameter  (maxlen = 2*maxres)
C     
      character*1 match(maxlen)
      dimension   ipath(2,maxres,maxres)
      dimension   smap(maxres,maxres)
      integer     iv1(maxlen),iv2(maxlen)
      logical     report
      
      common /vectors/ v1(maxres,3),v2(maxres,3),smapori(maxres,maxres)
      
c------------------------------------------------------------------------------
c     build the original score map
c------------------------------------------------------------------------------
 
      gapi = 7.00
      gape = 0.45
      
     
      call simat(nv1,nv2,smap,maxl,gapi,gape)
      
c------------------------------------------------------------------------------
c     Build penalty weighted maximum match matrix
c------------------------------------------------------------------------------
 
      do i = nv1-1, 1, -1
         look1= min(nv1,i+2+maxl)
         do j = nv2-1, 1, -1
            ii=i+1
            jj=j+1
            smax = smap(ii,jj)
            imax = ii
            jmax = jj
            do ii = i+2,look1
               value = smap(ii,jj) - gapi - float(ii-i-1)*gape
               if (value .gt. smax) then
                  smax = value
                  imax = ii
                  jmax = jj
               endif
            enddo
            ii=i+1
            look2= min(nv2,j+2+maxl)
            do jj = j+2,look2
               value = smap(ii,jj) - gapi - float(jj-j-1)*gape
               if (value .gt. smax) then
                  smax = value
                  imax = ii
                  jmax = jj
               endif
            enddo
            smap(i,j) = smap(i,j) + smax
            ipath(1,i,j) = imax
            ipath(2,i,j) = jmax
         enddo
      enddo
      
c------------------------------------------------------------------------------
c     Find the path
c------------------------------------------------------------------------------
 
      sc = 0.0
      nvnew = 0
      ngaps = 0

      i_p = 0
      j_p = 0
C
20    imax = i_p + 1
      jmax = j_p + 1
      smax = smap(i_p+1,j_p+1)
C
      do j = j_p+2, nv2
        value = smap(i_p+1,j)
        if (value .gt. smax) then
           smax = value
           imax = i_p+1
           jmax = j
        endif
      enddo
C     
      do i = i_p + 2, nv1
        value = smap(i,j_p+1)
        if (value.gt. smax) then
          smax = value
          imax = i
          jmax = j_p+1
       endif
      enddo
C
      ns1p1 = nv1+1
      ns2p1 = nv2+1
      do while (imax .ne. ns1p1 .or. jmax .ne. ns2p1)

         do i = i_p+1, imax-1
            nvnew = nvnew + 1
            iv1(nvnew) = i
            iv2(nvnew) = 0
            match(nvnew) = ' '
            ngaps = ngaps+1
         enddo
         do j = j_p+1, jmax-1
            nvnew = nvnew + 1
            iv1(nvnew) = 0
            iv2(nvnew) = j
            match(nvnew) = ' '
            ngaps = ngaps+1
         enddo
         
         sc = sc + smapori(imax,jmax)
        nvnew = nvnew + 1
        iv1(nvnew) = imax
        iv2(nvnew) = jmax
        if (smapori(imax,jmax).gt.3.0) then
           match(nvnew) = '|'
        else
           match(nvnew) = ' '
        endif
C     
        if (imax .eq. nv1 .or. jmax .eq. nv2) goto 50
        i_p = imax
        j_p = jmax
        imax = ipath(1,i_p,j_p)
        jmax = ipath(2,i_p,j_p)
        
      enddo
C     
 50   do i = imax+1,nv1
         nvnew = nvnew + 1
         iv1(nvnew) = i
         iv2(nvnew) = 0
         match(nvnew) = ' '
         ngaps = ngaps+1
      enddo
C     
      do j = jmax+1,nv2
         nvnew = nvnew + 1
         iv1(nvnew) = 0
         iv2(nvnew) = j
         match(nvnew) = ' '
         ngaps = ngaps+1
      enddo
      
c     pgaps = (float(ngaps) / float(nseqnew))*100.0
      if (.not.report) return
      
c------------------------------------------------------------------------------
c     Path found. Compute matches and write results
c------------------------------------------------------------------------------
      
c     do i = 1, nseqnew
c     if (aseq1(i) .ne. '.' .and. aseq1(i) .eq. aseq2(i)) then
c         match(i) = '|'
c       else
c         match(i) = ' '
c       endif
c     enddo
c
c     write(*,*)' '
c     write(*,*) 'Original Sequences:'
c     write(*,*) '-------------------'
c     write(*,*)' '

c     numseq = max(nseq1,nseq2)
c     ioff = 0
c 60  ibeg = ioff + 1
c     if (numseq .ge. ibeg+49) then
c       iend = ibeg+49
c     else
c       iend = numseq
c     endif
c     write(*,300) pname1,(cseq1(i),i=ibeg,iend)
c     write(*,300) pname2,(cseq2(i),i=ibeg,iend)
c     write(*,*)' '
c     ioff = iend
c     if (iend .lt. numseq) goto 60
C
c     write(*,*) 'Aligned Sequences:'
c     write(*,*) '-------------------'
c     write(*,*) ' '
c     ioff = 0
c 70  ibeg = ioff + 1
c     if (nseqnew .ge. ibeg+49) then
c       iend = ibeg+49
c     else
c       iend = nseqnew
c     endif
c     write(*,300) pname1,(aseq1(i),i=ibeg,iend)
c     write(*,300) ' ',(match(i),i=ibeg,iend)
c     write(*,300) pname2,(aseq2(i),i=ibeg,iend)
c     write(*,*) ' '
c     ioff = iend
c     if (iend .lt. nseqnew) goto 70
c     write(*,350) nsame
c     write(*,400) float(nsame) / float(min(nseq1,nseq2)) * 100.0
c     write(*,500) sc, sc-smax, smax
C
      return
C
c100  format(1x,100a3)
c120  format(a1,x,100f3.0)
c200  format(2i3)
c300  format(a12,x,5(10a1,x))
c350  format('Aligned Identical Residues: ',i3)
c400  format('Sequence Identity = ',f5.1,'%')
c500  format('Score = ',f6.1,'  Penalty = ',f6.1,'  Net score = ',f7.1)
C
      end
c
c========================================================================
c
      subroutine simat(nv1,nv2,smap,maxl,gapi,gape)
 
c-------------------------------------------------------------------------
c --- Compute similarity matrix between ca-ca vectors ---
c ---                      10x speed enhancement by cems 01            ---
c-------------------------------------------------------------------------

      implicit REAL*8 (a-h,o-z)
      
      parameter  (maxres = 3000)
      parameter  (maxfrag= maxres) ! 100)
      parameter  (maxlen = 2*maxres)
      
       integer npoints
       real*8 xre(0:maxres),xrp(0:maxres)
       real*8 xse(3,0:maxres),xsp(3,0:maxres)
       real*8 xm(3,3,0:maxres)
       real*8 det2,rms2_fast3,ev(3)
      
 
      
C
      dimension  smap(maxres,maxres)
      dimension  xp(3,maxfrag),xe(3,maxfrag)  ! ,w(maxfrag)

      common /vectors/ v1(maxres,3),v2(maxres,3),smapori(maxres,maxres)

    
      
      l = 7
      r = sqrt(float(l))
      r = sqrt(2.0-(2.84/r))

      do i = 1, nv1
         do j = 1, nv2
            smap(i,j) = 0.0d0
            smapori(i,j) = 0.0d0
         enddo
      enddo
      
      smapmin= 1000.0
      smapmax=-1000.0
      
      
                                ! reverse bass-ackward array indexing scheme.
      do n = 1, nv2
         do j = 1, 3
            
            xe(j,n) = v2(n,j)
         enddo
      enddo 
      
      do n = 1, nv1                  
         do j = 1, 3
            
            xp(j,n) = v1(n,j)
         enddo
         
      enddo
      
                                ! index overlap regions
c imagine the two sequences written one above the other so that the last L residues
c overlap.  The tail of 2 overlapping the beginiing of 1.
c if we number these and label the bin at the beginging of 1 as bin zero
c then the begining of 2 is at position minus nv2-L
c if we now slide 2 along, keeping 1 fixed then when we are done so that the
c last L residues of 1 overlap the first L residues of 2 then the tip of
c 2 will lie at bin Nv1-L+1
c  
c let's loop over all those shifts.

c now remember these are the differnce vectors not the residues so nv1 and nv2 are one less.


      do ihh = -(nv2-l),nv1-l+1
         
        iis = max0(1-ihh,1)  ! starting point of overlap for seq 2
        jjs = max0(1,ihh)       ! starting point of overlap for seq 1
        npoints=min0(nv2-iis+1,nv1-jjs+1) ! number of overlapping points
        
        call rms_setup(npoints,xp(1,jjs),xe(1,iis),xm,xre,xrp,xse,xsp)
        
                                ! stride widows through overlap regions 
        do kk=1,npoints-l+1
           ii = iis+kk-1
           jj = jjs+kk-1
           
           kkk = kk-1
           jjj = kk+l-1
              
           urms = rms2_fast3(l,xm(1,1,kkk),xre(kkk),xrp(kkk),xse(1,kkk),xsp(1,kkk),
     &                   xm(1,1,jjj),xre(jjj),xrp(jjj),xse(1,jjj),xsp(1,jjj))
           !urms = rms2_fast2(kk,kk+l-1,xm,xre,xrp,xse,xsp,det2,ev) ! rms of kk->kk+L-1 residues
           
           sim = ((r-urms)/r)*10.0
           if (sim.lt.0.0d0) sim = 0.0d0
           smapmin=min(smapmin,sim)
           smapmax=max(smapmax,sim)
           smap(jj,ii) = sim
           smapori(jj,ii) = sim
           
c           write(0,*) jj,ii,ihh,kk,npoints,iis+kkk-1,jjs+jjj-1
        enddo
        
      enddo
      
      
      maxl = int((smapmax-smapmin-gapi)/gape)+2
      
      
      return
      end
      
c=========================================================================
c
      subroutine getpdb_names(max_files,max_path_length,filename,pdb_path,pdb_list,nfiles)
      implicit none
c-------------------------------------------------------------------------
c --- processes file to see it it is alist file
c------------------------------------------------------------------------- 
      ! input
      integer max_files,max_path_length   
      character*(*)filename
      ! output
      character*(*)pdb_path,pdb_list(max_files)
      integer nfiles
      ! local
      logical name_exists
      integer istatus
      character*(512) line_buffer
      integer i,iEOF
      parameter  (iEOF = -1)
      ! functions
      integer lchar
            
      
            inquire(file=filename,exist=name_exists)
            if (.not.name_exists) then
              write(6,'(/,a,/)') 'ERROR: file ',filename,' does not exist.'
              call my_halt()
            endif
      
c   open file
         
        open(unit=20,file=filename,status='old')
        
        read(20,"(A)") line_buffer
        
         if ( line_buffer(1:12) .eq. "MAMMOTH List" ) then
             !"its a mammoth list file not a pdb file"
             i=1
             read(20,"(A)",iostat=istatus) pdb_path        ! read path header
             if (istatus.eq.iEOF) goto 778
             
             if (lchar(pdb_path).eq.max_path_length) then
                write(0,*) "path is too long",pdb_path
                call my_halt()
             endif
             
              if (lchar(pdb_path).gt.0) pdb_path = pdb_path(1:lchar(pdb_path))//'/'  ! assume unix    
             
  777          continue
             read(20,"(A)",iostat=istatus) pdb_list(i)     ! read list
             
             if (istatus.eq.iEOF) goto 778
             if (lchar(pdb_list(i)).lt.1) goto 777  ! its a blank line in file
             i=i+1
             if (i.gt.max_files) then
               write(0,*) "List file has too many entries max_files=",max_files,i
               call my_halt()
             endif
             goto 777
  778          continue
             i=i-1
             if (i .le. 1 ) then
                 write(0,*) "List file is empty ",filename(1:lchar(filename))
                 call my_halt()
             endif
             nfiles = i
            
         else
            ! its not a mammoth list file
            nfiles=1
            pdb_path = ""
            pdb_list(1) = filename
            
          
         endif
         close(unit=20)
        
         return
         end
      
      
c
c=========================================================================
c
      subroutine getpdb(NR,xyz,nunit,iatom,resnum,resnam,atnam,nc,ica)
c
c-------------------------------------------------------------------------
c --- Read coordinates (all) from the pdb file. ---
c --- Store Calphas as a pointer ---
c-------------------------------------------------------------------------
c
      implicit REAL*8 (a-h,o-z)
c
      parameter  (maxatom= 3000)
      parameter  (maxline=30000)
      parameter  (maxres = 3000)
c
      dimension     xyz(maxatom,3)
      integer       iatom(maxatom), resnum(maxatom)
      integer       ica(maxres)
      character*4   atnam(maxatom), resnam(maxatom)
      character*80  buffer*80
c
      ncalfa = 0
      natom  = 0
      do 5 k = 1, maxline
        read(nunit,'(a80)',end=6) BUFFER
        if (BUFFER(1:6).eq.'ENDMDL') then
          write(6,*) 'WARNING: pdb entry contains multiple models'
          write(6,*) 'WARNING: only the first model will be used '
          goto 6
        endif
        if (BUFFER(1:4).eq.'ATOM') THEN
          read(BUFFER,'(30x,3f8.3)') x,y,z
          if (x.eq.0.0d0.and.y.eq.0.0d0.and.z.eq.0.0d0) goto 5
          natom = natom + 1
          read(BUFFER,'(5x,i6,1x,a4,1x,a3,2x,i4,4x,3f8.3)')
     &    iatom(natom), atnam(natom), resnam(natom),
     &    resnum(natom), (xyz(natom,j), j=1,3)
          if((buffer(14:15).eq.'CA').or.
     &       (buffer(13:14).eq.'CA').or.
     &       (buffer(15:16).eq.'CA'))then
            ncalfa = ncalfa + 1
            ica(ncalfa) = natom
          endif
        else if (BUFFER(1:4).eq.'TER ') then
          goto 6
        endif
    5 continue
c
    6 NR = natom
      nc = ncalfa
      if (nc .eq. 0) then
        write(6,*) 'ERROR: No C-alphas selected'
        call my_halt()
      endif
c
      return
      end
c
c===============================================================================
c
      subroutine setime

c-------------------------------------------------------------------------------
c     ##  subroutine setime  --  initialize elapsed CPU time clock  ##
c         "setime" initializes the elapsed interval CPU timer
c-------------------------------------------------------------------------------

      implicit none
      real*8 cputim
      common /chrono/ cputim
c
c initialize interval at elapsed CPU time for current job
c
      call clock (cputim)

      return
      END
c
c===============================================================================
c
      subroutine getime (elapsed)

c-------------------------------------------------------------------------------
c     ##  subroutine getime  --  get elapsed CPU time in seconds  ##
c        "getime" gets elapsed CPU time in seconds for an interval
c-------------------------------------------------------------------------------

      implicit none
      real*8 cputim
      common /chrono/ cputim
      real*8 elapsed
c
c  elapsed time for the interval is the current total CPU
c  time minus the total time at the start of the interval
c
      call clock (elapsed)
      elapsed = elapsed - cputim

      return
      END
c
c===============================================================================
c
      subroutine clock (seconds)

c-------------------------------------------------------------------------------
c     ##  subroutine clock  --  find elapsed time for current job  ##
c        "clock" determines elapsed CPU time in seconds since the
c         start of the job; only one of the implementations should
c         be activated by removing comment characters from the code
c-------------------------------------------------------------------------------

      implicit none
      real*8 seconds
c
c Unix machines have access to the "etime" intrinsic,
c this code works for Sun, DEC, SGI, Convex and others
c
      real etime,times(2)
      seconds = dble(etime(times))
c
      return
      END
c
c================================================================================
c
      subroutine storetime(subname,ttime,ntime)
c
      implicit     none
      integer      maxcall,ntime
      parameter    (maxcall = 100)
      real*8       time,ttime,ftme
      character*80  part
      character*(*) subname

      common /timings/ ftme(maxcall),part(maxcall)

      call getime(time)

      ntime = ntime + 1
      part(ntime) = subname
      ftme(ntime) = time
      ttime = ttime + time

      return
      END
c
c================================================================================
c
      subroutine printime(ntime)

      implicit      none
      integer       maxcall,ntime,i
      parameter     (maxcall = 100)
      real*8        ftme
      character*80  part

      common /timings/ ftme(maxcall),part(maxcall)

      write(6,*) ' '
      write(6,*) '---------'
      write(6,*) ' Timings '
      write(6,*) '---------'
      write(6,*) ' '

      do i = 1, ntime
        write (6,33)  part(i),ftme(i)
      enddo
      write(6,*) ' '
  33  format ('   < ',a30,f12.3,' sec >')

      return
      END
c
c================================================================================
c
      subroutine help
c
c-------------------------------------------------------------------------------
c --- help messages
c-------------------------------------------------------------------------------
c
      implicit REAL*8 (a-h,o-z)

      write(0,*)
      write(0,*) 'MAMMOTH: The valid command parameters are:'
      write(0,*) 'REQUIRED args'
      write(0,*) '-p <predicted conformation file>'
      write(0,*) '-e <experimental conformation file>'
      write(0,*) '-o <output information>'
      write(0,*) 'OPTIONAL args'
      write(0,*) '-t <threshold z-score> for output of pdbs and SS alignment'
      write(0,*) '-r <1 or 0>   make pdb file  (default = 1,yes  i.e. true)'
      write(0,*) '-v <1 or 0>   verbose output (default = 1,yes) (false also sets -t => inifity)'
	  write(0,*) '-c <1 or 0>   compute contact order (default = 0,no)'

      write(0,*) 'OPTIONAL BULK LIST PROCESSING MODE:'
      write(0,*) 'if the first 12 characters of -e or -p file are "MAMMOTH List"'
      write(0,*) 'then the next line is interpreted as a path to a directory and'
      write(0,*) 'all subsequent lines are intepreted as names of files to process'
      write(0,*)
        call my_halt()

      end
c
c===============================================================================
c
      integer function lchar(string)
c
c-------------------------------------------------------------------------
c --- Returns the position of the last printable character (not a blank)
c --- in the string STRING
c-------------------------------------------------------------------------
c
      implicit REAL*8 (a-h,o-z)

      character*(*)  string
c      integer j
C

c      most compilers have a len_trim intrinsic.  if yours does uncomment these lines

c      lchar = len_trim(string)
c      if (lchar.eq.1) then  ! test stupid  fortran default for empty string
c        if (string(1:1).eq.'') lchar = 0
c      endif
c      return

c      compliers having len_trim(): alpha(compaq), Absoft, G77
c      compilers lacking len_trim(): PGI
      
      ! the following is provided if your compiler does not have len_trim intrinsic
     
      do j=len(string),1,-1
        
        if (string(j:j).ne.' ' ) then
         
          lchar=j
          ! unfortunately fortrans differ on how they treat a blank string
          ! absoft for example will not match the empty string to " " correctly, when
          ! that is the first character of a string. (geeeesh!)  so we trap this here:
          if (j.eq.1.and.string(1:1).eq.'' ) lchar=0
          
          return
         
       end if
       
        
      end do
      lchar=0
      return
      
      
      end
c
c==========================================================================

       subroutine  my_halt()
       ! halts after closing files
       write(6,*) "<MAMMOTH> HALT"
       close(6)
       close(20)
       close(30)
       
       stop
       end
       
