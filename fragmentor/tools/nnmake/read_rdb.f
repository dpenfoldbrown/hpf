c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 3381 $
c  $Date: 2003-09-10 19:15:00 -0400 (Wed, 10 Sep 2003) $
c  $Author: rohl $

      subroutine read_rdb(error)
c reads the rdb file for current protein
c rdb file are relative sequence propensity from the SAM (santa cruz) server
c flags errors if residues dont match fasta.
      implicit none
      include 'path_defs.h'
      include 'param.h'
      include 'structure.h'
      
      character*1 rdb_aa(max_res)
      integer iunit,i,nread
      logical error
      
      iunit=16
      write(0,*) 'reading rdb file',read_rdb_path(1:read_rdb_l) 
     #     //lower_name//chain_letter//'.rdb'
      open(iunit,file=read_rdb_path(1:read_rdb_l) 
     #     //lower_name//chain_letter//'.rdb',status='old',iostat=i) 
      if (i.ne.0) then 
         write(0,*) 'no rdb file'
         error=.true.
         return
      endif 
        
      call read_rdb_single(iunit,rdb_aa,rel_rdbH,Rel_rdbE,Rel_rdbL,
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
         rel_rdbH(i) = rel_rdbH(i)
         rel_rdbE(i) = rel_rdbE(i)
         rel_rdbL(i) = rel_rdbL(i)
         
         ss_rdb(i) = "L"
         if ( rel_rdbH(i).gt.rel_rdbE(i) .and. 
     #        rel_rdbH(i).gt.rel_rdbL(i) ) ss_rdb(i) = "H"
         if ( rel_rdbE(i).gt.rel_rdbH(i) .and. 
     #        rel_rdbE(i).gt.rel_rdbL(i) ) ss_rdb(i) = "E"
         if ( rel_rdbL(i).gt.rel_rdbH(i) .and. 
     #        rel_rdbL(i).gt.rel_rdbE(i) ) ss_rdb(i) = "L"
         
      enddo
      error=.false.
      return
      end
c--------------------------------------------------------------------------
      subroutine read_rdb_single( iunit,aa,relH,RelE,RelL,max_res,i)
      implicit none
c input
      integer iunit
      integer max_res
c output
      character*1 aa(max_res)	
      real relH(max_res),relE(max_res),relL(max_res)
      integer i
	
c local
      character*100 buffer
c	local
      integer ii
      i=0
       buffer='   ' 
      do while(buffer(1:3).ne.'Pos')	
         read(iunit,'(a3)',end = 99,err=100) buffer
         
      enddo
	 
      read(iunit,*,end=99,err=100) ! read one more line of header
	 
      do while(i.lt.max_res)	
         read(iunit,'(A100)',end = 99,err=100) buffer
         i=i+1
c absoft format
         read(buffer,'(i8,8x,a1,f15.3,2f16.3)',err=200)
     #        ii,aa(i),relE(i),relH(i),relL(i)
	goto 300
c pgi format
 200         read(buffer,'(i8,1x,a1,1x,f8.3,1x,f8.3,1x,f8.3)',err=400)
     #        ii,aa(i),relE(i),relH(i),relL(i)
 300	continue
      enddo
        
 99   continue
      return

 100  write(0,*) 'error reading rdb file',buffer
      stop
 400  write(0,*)'error parsing buffer'//buffer
	write(0,*)'*'//buffer
      end

