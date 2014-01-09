c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 7305 $
c  $Date: 2006-01-17 16:46:00 -0500 (Tue, 17 Jan 2006) $
c  $Author: pbradley $

      subroutine homologs 
      implicit none 
      include 'path_defs.h'
c---------------------------------------------------------------------
      include 'param.h'
      include 'structure.h' 
      character ch1*1,ch2*1,prot1*4,prot2*4  
      integer i,j,new,damn   
c---------------------------------------------------------------------
c look in two places then die if cant find it
      if (chain_letter.eq.' ') chain_letter='_' 
      num_homologs=1 
      homolog_id(num_homologs)=lower_name 
      homolog_ch(num_homologs)=chain_letter 
      write(0,*) "read homolog 18 ",homolog_path(1:homolog_l)//
     # lower_name//chain_letter//'.'//'homolog_vall'
       open(18,file= homolog_path(1:homolog_l)//
     # lower_name//chain_letter//'.'//'homolog_vall',
     # status='old',iostat=damn)

  
c       if (damn.ne.0) then 
c        write(0,*)  "read homolog 18 ",
c     #  '/usr5/users/ksimons/database/select/cull_7_99.homolog'  
c         write(60,*) 'cant find homologn, will try exclude list.'
c       open(18,file=  
c     #  '/usr5/users/ksimons/database/select/cull_7_99.homolog'
c     #   ,status='old',iostat=damn)
 
 
        if (damn.ne.0) then
 
	write(0,*) " WARNING:: unable to find/open homolog file"
        return
       endif
c      endif

      do while(num_homologs.lt.max_homologs)
         read(18,'(2(a4,a1,1x))',end=51)prot1,ch1,prot2,ch2 
car 0 _ and space  all indicate no chain
car convert all to _ which is used in the vall
         if (ch1.eq.'0') ch1='_'
         if (ch1.eq.' ') ch1='_'
         if (ch2.eq.'0') ch2='_'
         if (ch2.eq.' ') ch2='_'
         
	 !write(0,*)prot1,ch1,prot2,ch2,"-",lower_name,chain_letter
         if ((prot1.eq.lower_name.and.ch1.eq.chain_letter).or.
     1        (prot1.eq.lower_name .and. ch1.eq.'_' 
     2         .and.chain_letter.eq.'0') .or. !since we change it above
     3        (prot2.eq.lower_name .and. ch2.eq.'_' 
     4         .and.chain_letter.eq.'0') .or. !since we change it above
     5        (prot2.eq.lower_name.and.ch2.eq.chain_letter)) then 
            new=1 
            do j=1,num_homologs 
               if (homolog_id(j).eq.prot1.and.
     #             homolog_ch(j).eq.ch1) new=0    
            enddo 
            if (new.eq.1) then 
               num_homologs=num_homologs+1 
               homolog_id(num_homologs)=prot1 
               homolog_ch(num_homologs)=ch1 
            endif 
            new=1 
            do j=1,num_homologs 
               if (homolog_id(j).eq.prot2.and.
     #             homolog_ch(j).eq.ch2) new=0    
            enddo 
            if (new.eq.1) then 
               num_homologs=num_homologs+1 
               homolog_id(num_homologs)=prot2 
               homolog_ch(num_homologs)=ch2 
            endif 
         endif 
      enddo
 51   continue
      close(18)
      if (num_homologs.eq.max_homologs) then 
         write(60,*)'not enough homologs space' 
         stop 
      endif 

      if (num_homologs.eq.0) return 
      write(60,*)'SIMILAR STRUCTURES IN PDB_SELECT' 
      write(60,*)'      removing ',num_homologs,' proteins' 
      do i=1,num_homologs   
         write(60,*)'should remove: ',homolog_id(i),homolog_ch(i)  
c         if (homolog_ch(i).eq.'_') homolog_ch(i)=' ' 
         new=0 
         do j=1,vall_size 
  
            if (frag_name(j).eq.homolog_id(i).and. 
     #          chain(j).eq.homolog_ch(i)) then 
               new=1 
               phi_list(j)=0.0 
               psi_list(j)=0.0 
            endif 
 

         enddo 
         if (new.eq.1) then  
            write(60,*)'does remove: ',homolog_id(i),homolog_ch(i)    
         endif 
      enddo 

 100  format(a8)
 200  format(3x,a4,a1,1x,a4,a1)  
      return 
      end 
