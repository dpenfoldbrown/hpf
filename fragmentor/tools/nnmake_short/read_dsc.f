c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 1.5 $
c  $Date: 2002/06/18 04:41:53 $
c  $Author: rohl $

      subroutine read_dsc(error)
      implicit none 
      include 'path_defs.h'
c------------------------------------------------------------------
      include 'structure.h' 
      integer i
      character indicator*8,lower2upper*1   
      integer iunit
      logical error
c------------------------------------------------------------------
      iunit=16
      if (chain_letter.eq.' ') chain_letter='_' 
      write(0,*) 'read_rdb',iunit,read_rdb_path(1:read_rdb_l) 
     #    //lower_name//chain_letter//'.rdb'
      open(iunit,file=read_rdb_path(1:read_rdb_l) 
     #    //lower_name//chain_letter//'.rdb',status='old',iostat=i) 
      if (i.ne.0) then 
         write(60,*)'no dsc file'
         error=.true.
         return
      endif 
      indicator='abcdefgh' 
      do while(indicator.ne.'NO.  RES') 
         read(iunit,'(a8)',end=201)indicator
      enddo 
      do i=1,total_residue    
         read(iunit,'(6x,a1,6x,a1,3(5x,f5.3))')seq_rdb(i),
     #      ss_rdb(i),rel_rdbH(i),rel_rdbE(i),rel_rdbL(i)  
      enddo 

 101  continue 
      close(iunit) 

      do i=1,total_residue 
         seq_rdb(i)=lower2upper(seq_rdb(i))  
         rel_rdbH(i)=rel_rdbH(i)
         rel_rdbE(i)=rel_rdbE(i)
         rel_rdbL(i)=rel_rdbL(i)
         if (seq_rdb(i).ne.seq(i)) then  
          write(60,*)'rdb mismatch: ',i,seq(i),' rdb: ',seq_rdb(i) 
          stop 
c        else 
c         write(60,*)i,'seq: ',seq(i,1),' rdb: ',seq_rdb(i) 
         endif 
         if (ss_rdb(i).eq.'C') ss_rdb(i)='L' 
      enddo 

      error=.false.
      return 

 201  write(60,*)'unable to read dsc file' 
      stop 

      end 

