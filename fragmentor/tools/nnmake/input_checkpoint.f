c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 2338 $
c  $Date: 2002-08-09 18:23:11 -0400 (Fri, 09 Aug 2002) $
c  $Author: rohl $

      subroutine input_checkpoint()   
      implicit none
c reads blastpgp -C made checkpoint file which contains
c position specific seq-weights and inteligent psuedocounts 
c from blosumn62 (i.e. the longer the alignment less
c the substitution table is used
c-------------------------------------------------------
c   MATCH SEGMENTS TO PROFILE
c-------------------------------------------------------
      include 'param.h'
      include 'structure.h'
      include 'path_defs.h'

      character name*5
      character filename*200
      integer h,i,j,num_res
      integer funit

      funit = 11

c  OPEN AND READ checkpoint FILE
   
      if (chain_letter.eq.' ') chain_letter = "_"
      name=lower_name//chain_letter
      if (chain_letter.eq.'_') chain_letter = " "

      filename=input_chkpt_path(1:input_chkpt_l)//name//'.checkpoint'
      write(0,*) "input checkpoint",funit,filename(1:index(filename,' '))
      open(funit,file= filename,status='old',iostat=h)
      if (h.eq.0) goto 10
 
      filename=input_chkpt_path(1:input_chkpt_l)//name//'.checkPROFILE'
      write(0,*) "input checkpoint",funit,filename(1:index(filename,' '))
      open(funit,file= filename,status='old',iostat=h)
      if (h.eq.0) goto 10
     
      write(0,*)" no checkpoint file for ",lower_name
      stop

 10   continue

C these files are cleaned w/ perl (so e-z!)
      read(funit,*)num_res
      if (num_res.gt.max_res) then
         write(0,*)'increase max_res in structure.h to ', num_res
         stop
      endif
      
      do i=1,num_res
         read(funit,'(a1,20(1x,f6.4))')seq(i),(profile(j,i),j=0,19)
      enddo
      close(funit)
         
c if pdb is not read then we use check point as source of 
c information anbout query sequence
      if (use_pdb.eq.'N') then
         total_residue=num_res
         do i = 1,total_residue
            residue1(i) = seq(i)  
         enddo 
      else
c     otherwise verify that check point and pdb agree
         if (total_residue.ne.num_res) then
            write(0,*) 'check point file length mismatch', total_residue,num_res
            stop
         endif
         do i = 1,total_residue
            if (residue1(i) .ne. seq(i) ) then
               write(0,*) 'check point file residue mismatch',i,residue1(i),seq(i)
            endif
            
         enddo
         
      endif
      if (total_residue.gt.max_res) then
         write(0,*)'increase max_res in structure.h to ', total_residue
         stop
      endif

      call convert1_num
   
      write(0,*)"protein_name, chain",lower_name, chain_letter
      return
      end
