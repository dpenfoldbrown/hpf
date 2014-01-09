c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 1.5 $
c  $Date: 2002/06/18 04:41:53 $
c  $Author: rohl $

      subroutine read_phd(error) 
      implicit none 
      include 'path_defs.h'
      include 'structure.h' 

car   output
      logical error

car local
      integer i,j 
      integer iunit
      character indicator*12,line_id*3,sec*3  
      integer rel_phd_int( max_res),rel_phdH_int(max_res)
      integer rel_phdE_int(max_res),rel_phdL_int(max_res)


      error=.true.
      iunit=16
      if (chain_letter.eq.' ') chain_letter='_'
       write(0,*) 'read_phd',iunit, read_phd_path(1:read_phd_l)
     #    //lower_name//chain_letter//'.phd'
       open(iunit,file= read_phd_path(1:read_phd_l)
     #    //lower_name//chain_letter//'.phd',status='old',iostat=i) 
      if (i.ne.0) then 
	 write(0,*) 'No file in path'
         write(60,*)'no phd file' 
         return
      endif 
      indicator='    abcdein:' 
      do while(indicator.ne.'protein:    '.and.
     #         indicator.ne.'    protein:') 
         read(iunit,'(1x,a12)',end=201)indicator
      enddo 
      do i=1,total_residue,60 
         line_id='NON' 
         sec='sss' 
         do while (line_id.ne.'AA ') 
            read(iunit,'(9x,a3,a3)',end=101)line_id,sec  
         enddo 
         backspace(iunit) 
         if (sec.eq.'   ') then  
            read(iunit,1000)(seq_phd(j),j=i,min(i+59,total_residue)) 
            read(iunit,1000)(ss_phd(j),j=i,min(i+59,total_residue)) 
            read(iunit,2000)(rel_phd_int(j),j=i,min(i+59,total_residue)) 
            read(iunit,*) 
            read(iunit,2000)(rel_phdH_int(j),j=i,min(i+59,total_residue)) 
            read(iunit,2000)(rel_phdE_int(j),j=i,min(i+59,total_residue)) 
            read(iunit,2000)(rel_phdL_int(j),j=i,min(i+59,total_residue)) 
         else 
            read(iunit,3000)(seq_phd(j),j=i,min(i+59,total_residue)) 
            read(iunit,3000)(ss_phd(j),j=i,min(i+59,total_residue)) 
            read(iunit,4000)(rel_phd_int(j),j=i,min(i+59,total_residue)) 
            read(iunit,*) 
            read(iunit,4000)(rel_phdH_int(j),j=i,min(i+59,total_residue)) 
            read(iunit,4000)(rel_phdE_int(j),j=i,min(i+59,total_residue)) 
            read(iunit,4000)(rel_phdL_int(j),j=i,min(i+59,total_residue)) 
         endif
      enddo 
 101  continue 
      close(iunit) 
      do i=1,total_residue
         rel_phd(i) =0.1*real(rel_phd_int(i))
         rel_phdH(i)=0.1*real(rel_phdH_int(i))
         rel_phdE(i)=0.1*real(rel_phdE_int(i))
         rel_phdL(i)=0.1*real(rel_phdL_int(i))
         if (seq_phd(i).ne.seq(i)) then  
            write(0,*) 'phd mismatch: see status file for details... '
            write(0,*) 'quitting'
            write(60,*)'phd mismatch: ',i,seq(i),' phd: ',seq_phd(i) 
            stop 
         else 
         endif 
         if (ss_phd(i).eq.' ') ss_phd(i)='L' 
      enddo 
      error=.false.
      return 
      
 201  write(60,*)'unable to read phd' 
      stop 
      
 1000 format(18x,60a1)
 2000 format(18x,60i1)
 3000 format(14x,60a1)     
 4000 format(14x,60i1)     

      end 

