c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 1.5 $
c  $Date: 2002/06/18 04:41:53 $
c  $Author: rohl $

      subroutine zscore_dme(ss_type)                    
      implicit none
c-------------------------------------------------------
      include 'structure.h' 
      include 'path_defs.h'
c-------------------------------------------------------
      integer ss_type
      integer q,j,size,i
      real dme_native
      external dme_native
c surbroutine prints out stats on to 200 frags for each size and paradigm
c assumes that best_nn contains the fragment in sorted order.  
c if not then use big_best_nn
c print out is arranged so that residue positions are across with lesser 
c ranked alternatives
c underneath.  z_score and dme alternate along row
c arrays for each score are arranged successively.
      integer iunit
      character*1 int_to_str
      character*1 ss_num
      
      real score_fragment
 
      iunit=70

      ss_num=int_to_str(ss_type,1)

      open(iunit,file = setup_path(1:setup_l)//'zscore_'//ss_num//'.'//
     #     extension(1:index(extension," ")-1)//'_'//lower_name//
     #     chain_letter,status='unknown')

      
      do q = 1,n_sizes
         size = sizes(q)
         write(iunit,*) "SIZE:",size
         do j = 1,max_nn
            
c  TEMPORARY CHANGE:  ZSCORE REPLACE BY ACTUAL SCORE.               
c     write(iunit,'(i4,1x,200(f7.3,1x,f7.3,2x))') j+1000*ss_type, 
c     #           (z_score(size,j,i),dme_native(i,best_nn(size,i,j),size)  
c     #           , i=1,total_residue-size)
            write(iunit,'(i4,1x,200(g9.3,1x,f7.3,2x))') 
     #           j+1000*ss_type, (score_fragment(i,best_nn(size,i,j),
     #           size,best_nn_ss_type(size,i,j)),
     #           dme_native(i,best_nn(size,i,j),size),
     #           i=1,total_residue-size)   
         enddo
      enddo
      
      close(iunit)

      return
      end

c--------------------------------------------------

      character*(*) function int_to_str(num,length)

car converts an integer to a string of characters of size "length"
car pads the left side with zeros as required
car truncates from left if num>10**(length-1)

      implicit none

      integer num
      integer length

      integer digit
      integer i
      character*(1) ch

      int_to_str=' '
      do i=1,length
         digit = mod(num/(10**(i-1)),10)
         ch = char(digit+ichar('0'))
         int_to_str= ch//int_to_str(1:index(int_to_str,' ')-1)
      enddo

      return
      end

