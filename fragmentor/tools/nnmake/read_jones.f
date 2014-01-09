c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 4019 $
c  $Date: 2004-02-17 17:03:17 -0500 (Tue, 17 Feb 2004) $
c  $Author: rohl $

      subroutine read_jones_1state(error)
      implicit none 
      include 'path_defs.h'
      include 'param.h'
      include 'structure.h' 

car output
      logical error
car local
      integer i,j 
      integer iunit
      integer line
      character*7  indicator
      character*200 filename


      iunit=16
c look in two places then die

      if (chain_letter.eq.' ') chain_letter='_' 
      filename=read_jones_path(1:read_jones_l) 
     #     //lower_name//chain_letter//'.jones'
      write(0,*) "read_jones 16 ",filename(1:index(filename,' '))
      open(iunit,file=filename,status='old',iostat=i) 
      if (i.eq.0) goto 10
      
      filename=read_jones_path(1:read_jones_l) 
     #     //lower_name//chain_letter//'.psipred'
      write(0,*) "read_jones 16 ",filename(1:index(filename,' '))
      open(iunit,file=filename,status='old',iostat=i) 
      if (i.eq.0) goto 10

      filename=read_jones_path(1:read_jones_l) 
     #        //lower_name//chain_letter//'.jones2_15_99'
      write(0,*) "read_jones 16 ",filename(1:index(filename,' '))
      open(iunit,file=filename,status='old',iostat=i) 
      if (i.eq.0) goto 10
      
      write(60,*)'no jones file' 
      error=.true.
      return

 10   continue

      line=0
 666  indicator='abcdef' 
      do while(indicator(1:5).ne.'Conf:') 
         line=line+1
         read(iunit,'(a7)',end=201)indicator
      enddo 
      if (llt(indicator(7:7),'0') .or. lgt(indicator(7:7),'9')) goto 666
      rewind(iunit)
      do i=1,line-1
          read(iunit,*)
       enddo      

      do i=1,total_residue,60 
         read(iunit,'(6x,60f1.0)')(rel_jon(j),j=i,min(i+59,total_residue)) 
         read(iunit,'(6x,60a1)')(ss_jon(j),j=i,min(i+59,total_residue)) 
         read(iunit,'(6x,60a1)')(seq_jon(j),j=i,min(i+59,total_residue)) 
         if (min(i+59,total_residue).lt.total_residue) then 
            read(iunit,*) 
            read(iunit,*) 
            read(iunit,*) 
         endif  
      enddo 
 101  continue 
      close(iunit) 
      do i=1,total_residue 
         if (ss_jon(i).eq.'C') ss_jon(i)='L'
         rel_jonH(i)=0.3-0.0333*real(rel_jon(i))
         rel_jonE(i)=0.3-0.0333*real(rel_jon(i))
         rel_jonL(i)=0.3-0.0333*real(rel_jon(i))
         if (ss_jon(i).eq.'H') rel_jonH(i)=0.3+0.0667*real(rel_jon(i)) 
         if (ss_jon(i).eq.'E') rel_jonE(i)=0.3+0.0667*real(rel_jon(i))
         if (ss_jon(i).eq.'L') rel_jonL(i)=0.3+0.0667*real(rel_jon(i))
      enddo 
      do i=1,total_residue 
         if (seq_jon(i).ne.seq(i)) then  
          write(0,*)'jones mismatch: ',i,seq(i),
     #       ' jon: ',seq_jon(i) 
          write(60,*)'jones mismatch: ',i,seq(i),
     #       ' jon: ',seq_jon(i) 
          stop 
         endif 
         if (ss_jon(i).eq.' ') ss_jon(i)='L' 
      enddo 
      error=.false.
      return 

 201  write(60,*)'unable to read jones file' 
      error=.true.
      return

      end 
c------------------------------------------------------------------

      subroutine read_jones(error)
      implicit none 
      include 'path_defs.h'
c------------------------------------------------------------------
      include 'param.h'
      include 'structure.h' 
      integer i
      integer iunit
      logical error
      real jonL,jonH,jonE
      character*80 line
c------------------------------------------------------------------
      iunit=16
      if (chain_letter.eq.' ') chain_letter='_' 
      write(0,*) 'read_jones',iunit,read_jones_path(1:read_jones_l) 
     #    //lower_name//chain_letter//'.psipred_ss2'
      open(iunit,file=read_jones_path(1:read_jones_l) 
     #    //lower_name//chain_letter//'.psipred_ss2',status='old',iostat=i) 
      if (i.ne.0) then 
         write(60,*)'no jones file'
         goto 301
      endif 
      line(1:4)='####'
      do while (line(1:4).ne.'   1')
         read(iunit,'(a80)',end=201)line
      enddo
      read (line,'(5x,a1,1x,a1,1x,3(2x,f5.3))')seq_jon(1),
     #      ss_jon(1),jonL,jonH,jonE
      
      do i=1,total_residue    
         if (i.gt.1) read(iunit,'(5x,a1,1x,a1,1x,3(2x,f5.3))',err=201)
     #      seq_jon(i),ss_jon(i),jonL,jonH,jonE
         rel_jonH(i)=jonH
         rel_jonE(i)=jonE
         rel_jonL(i)=jonL
         if (ss_jon(i).eq.'C') ss_jon(i)='L' 
         if (seq_jon(i).ne.seq(i)) then  
            write(0,*)'jones mismatch: ',i,seq(i),' jones: ',seq_jon(i)
            goto 301
         endif
         
      enddo 

 101  continue 
      close(iunit) 

      error=.false.
      return 

 201  write(0,*)'unable to read jones file' 
 301  continue
      error=.true.
      call read_jones_1state(error)
      return

      end 

