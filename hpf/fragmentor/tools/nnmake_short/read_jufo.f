c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 1.3 $
c  $Date: 2002/06/18 04:41:53 $
c  $Author: rohl $

      subroutine read_jufo(error)
      implicit none 
      include 'path_defs.h'
c------------------------------------------------------------------
      include 'structure.h' 
      integer i
      character indicator*8
      integer iunit
      logical error
c------------------------------------------------------------------
      iunit=16
      if (chain_letter.eq.' ') chain_letter='_' 
      write(0,*) 'read_jufo',iunit,read_jufo_path(1:read_jufo_l) 
     #    //lower_name//chain_letter//'.jufo_ss'
      open(iunit,file=read_jufo_path(1:read_jufo_l) 
     #    //lower_name//chain_letter//'.jufo_ss',status='old',iostat=i) 
      if (i.ne.0) then 
         write(60,*)'no jufo file'
         error=.true.
         return
      endif 
      indicator='abcdefgh' 
      do i=1,total_residue    
         read(iunit,'(5x,a1,1x,a1,1x,3(2x,f5.3))')seq_jufo(i),
     #      ss_jufo(i),rel_jufoL(i),rel_jufoH(i),rel_jufoE(i)  
      enddo 

 101  continue 
      close(iunit) 
      do i=1,total_residue 
         rel_jufoH(i)=rel_jufoH(i)
         rel_jufoE(i)=rel_jufoE(i)
         rel_jufoL(i)=rel_jufoL(i)
         if (seq_jufo(i).ne.seq(i)) then  
          write(60,*)'jufo mismatch: ',i,seq(i),' jufo: ',seq_jufo(i) 
          stop 
c        else 
c         write(60,*)i,'seq: ',seq(i),' jufo: ',seq_jufo(i) 
         endif 
         if (ss_jufo(i).eq.'U') ss_jufo(i)='L'
         if (ss_jufo(i).eq.'S') ss_jufo(i)='E' 
      enddo 

      error=.false.
      return 

 201  write(60,*)'unable to read jufo file' 
      stop 

      end 

