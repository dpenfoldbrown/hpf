c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 1.10 $
c  $Date: 2003/09/10 23:15:00 $
c  $Author: rohl $

      subroutine read_prof_rdb(error)
c     reads the prof prediction for the protein
c     this is read into the "phd" arrays since prof
c     is an updated version of phd
c
c     this reads an "rdb"-style file made by a python conversion script
c     from the original messy prof output
c
c     flags errors if residues dont match fasta.
      implicit none
      include 'path_defs.h'
      include 'param.h'
      include 'structure.h'
      
      character*1 rdb_aa(max_res)
      integer iunit,i,nread
      logical error
      
      iunit=16
      if (chain_letter.eq.' ') chain_letter='_'
      write(0,*) 'reading prof_rdb file',read_rdb_path(1:read_rdb_l) 
     #     //lower_name//chain_letter//'.prof_rdb'
      open(iunit,file=read_rdb_path(1:read_rdb_l) 
     #     //lower_name//chain_letter//'.prof_rdb',status='old',iostat=i) 
      if (i.ne.0) then 
         write(0,*) 'no prof_rdb file'
         error=.true.
         return
      endif 
        
      call read_rdb_single(iunit,rdb_aa,rel_phdH,Rel_phdE,Rel_phdL,
     #     max_res,Nread)
      close(iunit)

 30   continue
      if (nread.ne.total_residue) then
         if (rdb_aa(nread).eq. ' ') then  !maybe have blank lines at end
            nread=nread-1  
            goto 30
         endif
         write(0,*) "STOPPING: rdb file size mismatched to total_residue",nread,total_residue
         stop
      endif
      
      do i=1,total_residue
         if (rdb_aa(i).ne.seq(i)) then
            write(0,*) "STOPPING: residue mismatch in rdb file res "
     #           ,i,rdb_aa(i),' should be ',seq(i)
            stop
         endif
         rel_phdH(i) = rel_phdH(i)
         rel_phdE(i) = rel_phdE(i)
         rel_phdL(i) = rel_phdL(i)
         
         ss_phd(i) = "L"
         if ( rel_phdH(i).gt.rel_phdE(i) .and. 
     #        rel_phdH(i).gt.rel_phdL(i) ) ss_phd(i) = "H"
         if ( rel_phdE(i).gt.rel_phdH(i) .and. 
     #        rel_phdE(i).gt.rel_phdL(i) ) ss_phd(i) = "E"
         if ( rel_phdL(i).gt.rel_phdH(i) .and. 
     #        rel_phdL(i).gt.rel_phdE(i) ) ss_phd(i) = "L"
         
      enddo
      error=.false.
      return
      end
