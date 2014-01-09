c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 6809 $
c  $Date: 2005-10-11 10:26:48 -0700 (Tue, 11 Oct 2005) $
c  $Author: dekim $

      subroutine read_barcode(error) 
      implicit none 
c----------------------------------------------------------------------
      include 'param.h'
      include 'structure.h' 
      include 'path_defs.h'


cdk from rosetta++ barcode_classes.cc
cdk torsion and SS csts
c----------------------------------------------------------------------
c   each cst starts with a tag-string which determines the rest of the parsing
c   current tags are:

c   PHI <rsd> <weight> <fvalue>
c   PSI <rsd> <weight> <fvalue>
c   OMEGA <rsd> <weight> <fvalue>
c   BB_BIG_BIN <rsd> <weight> <cvalue>
c   SS <rsd> <weight> <cvalue>
c   BB_SMALL_BIN <RSD> <weight> <fvalue> <fvalue> <fvalue> <fvalue> ***1
c   BB_CLUSTER <RSD> <weight> <fvalue> <fvalue> <fvalue> <fvalue> ***3

c   fvalue means floating point, ivalue means int, cvalue means char

c   DANGER DANGER DANGER DANGER:
c  
c   ***1 the order of the values should be: min_phi,max_phi,min_psi,max_psi
c   ***2 the order of the values should be: min_chi1,max_chi1,min_chi2,max_chi2
c   ***3 the order of the values should be: phi,psi,omega,thresh
c   BUT: if the bin spans the 180 divider, then min_phi should be greater than
c        max_phi,etc; ie the check becomes: (phi >= min_phi || phi <= max_phi)
c        rather than (phi >= min_phi && phi <= max_phi)

c   the format for the constraints file is: (ignore the leading "// ")

c   <feature1-tag> <flavor1-freq> <constraint1> {<constraint2> <constraint3> ... }
c   <feature1-tag> <flavor2-freq> <constraint1> {<constraint2> <constraint3> ... }
c   ...
c   <feature2-tag> <flavor1-freq> <constraint1> {<constraint2> <constraint3> ... }
c   <feature2-tag> <flavor2-freq> <constraint1> {<constraint2> <constraint3> ... }
c   ...

c   here each feature-tag is just a unique, non-white-space string that uniquely
c   identifies a common set of flavors

c   <constraint1> means insert the relevant info for constraint1, eg
c   <constraint1> could be "PHI 12 1.0 46.7" or "SS 24 1.0 E"
c  
c----------------------------------------------------------------------

      integer i
      integer iunit
      logical error
      character*10 tag
      real freq
      character*10 cst
      integer cstpos
      integer cstweight
      character*10 cstvalue
c------------------------------------------------------------------
      iunit=16

      write(0,*) 'read_barcode: ',barcode_name
      open(iunit,file=barcode_name,status='old',iostat=i)
      if (i.ne.0) then
         write(0,*) 'no barcode file'
         error=.true.
         return
      endif

      do
          read (iunit,*,end=101) tag, freq, cst, cstpos, cstweight, cstvalue
          write(0,*) 'barcode feature:',tag,' cst:',cst,' pos:',cstpos,' val:',cstvalue
          if (cst.eq.'BB_BIG_BIN') then 
              force_bb_big_bin(cstpos) = cstvalue(1:1)
          elseif (cst(1:2).eq.'SS') then
              force_ss(cstpos) = cstvalue(1:1)
          endif
      enddo
c 101   do i=1,total_residue
c          write(0,*) 'pos',i,' bin ',force_bb_big_bin(i) 
c          write(0,*) 'pos',i,' ss ',force_ss(i) 
c      enddo
c      return
 
 101   return

      end



