c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 1.3 $
c  $Date: 2002/06/18 04:41:53 $
c  $Author: rohl $

      subroutine write_names()  

      implicit none

      include 'structure.h' 
      include 'path_defs.h'
c-----------------------------------------------------
c------------------------------------------------------ 
      integer i,j,k,size,place,do_it,tot_names,
     #        max_names 
      integer iunit
      parameter (max_names=  75)
      integer tot_list(max_names) 
      character list(max_names)*4 
      character*200 filename
c------------------------------------------------- 

      iunit=10
      tot_names=0
      filename=write_names_path(1:write_names_l)//
     #     'names.'//extension(1:index(extension," ")-1)
     #     //prefix_code//lower_name//chain_letter
      write(0,*)'write_names',iunit,filename(1:index(filename,' '))
      open(iunit,file=filename,status='unknown') 
c     size=8 
      size=sizes(1)    
      do i=1,n_sizes
         if (sizes(i).eq.8) size=8    !8-mers  if available
      enddo
      do i=1,total_residue-size                       
        do j=1,max_nn    
           place=best_nn(size,i,j)      
           do_it=1 
           do k=1,tot_names 
              if (frag_name(place).eq.list(k)) then 
                 do_it=0  
                 tot_list(k)=tot_list(k)+1 
              endif 
           enddo 
           if (do_it.eq.1.and.tot_names.lt.max_names) then 
              tot_names=tot_names+1 
              tot_list(tot_names)=1 
              list(tot_names)=frag_name(place) 
           endif 
         enddo
      enddo 
      do i=1,tot_names 
         write(iunit,'(i6,1x,a4,1x,i6)')i,list(i),tot_list(i)  
      enddo 
      close(iunit) 

      return 
      end 
 

