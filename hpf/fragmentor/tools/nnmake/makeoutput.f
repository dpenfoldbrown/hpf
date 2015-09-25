c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 11090 $
c  $Date: 2006-11-02 15:17:28 -0500 (Thu, 02 Nov 2006) $
c  $Author: dekim $

      subroutine makeoutput(tag)    
c------------------------------------------------------
c  THIS SUBROUTINE WILL WRiTE  THE *.NN FILE 
c------------------------------------------------------ 
      implicit none 
      include 'path_defs.h'
      include 'param.h'
      include 'structure.h' 

car  tag_dummy should be a reflections of the quota_cap mixture?

      integer tag

      integer i,j,k,q,size,place     
      real dme,dme_native,dme_f,dme_full,score,score_fragment,noe
      real  test_frag_noe
      character*3 tag_str,num2str 
      character tag_string*9 
      integer iunit
      character*2 int_to_str
      character*2 sizenum
      character*100 filename
      real retrieve_frag_dipolar_score,dipolar
c------------------------------------------------- 
      iunit=10
 
      tag_str=num2str(tag)
 
      if (chain_letter.eq.' ') chain_letter='_' 
      
      do q=1,n_sizes
         size=sizes(q) 
         sizenum=int_to_str(size+1,2)
         filename=make_output_path(1:make_output_l)//
     #        prefix_code//lower_name//chain_letter//sizenum
     #         //'_'//tag_str(2:3)//'.'//extension(1:index(extension,' ')-1)
         write(0,*) 'make_output',iunit,filename(1:index(filename,' '))
         
         open (iunit,file=filename,status='unknown',iostat=j) 
         
         if (j.ne.0) then
            write(0,*) "file open in makeoutput failed"
            stop
         endif
         do i=1,total_residue-size 

            write(iunit,2000)i,max_nn
            write(iunit,*) 
            
            call prepare_constraints_score(i,size+1,vall_size)
            
            do j=1,max_nn    
               place=best_nn(size,i,j)  
               
               score = score_fragment(i,place,size,
     #              best_nn_ss_type(size,i,j)) 
               if (score.gt.9999) then
                   score = 999.999
               endif
               dme=dme_native(i,place,size) 
               dme_f= dme_full(i,place,size)
               noe = test_frag_noe(place)
               dipolar= retrieve_frag_dipolar_score(place)
               write(tag_string,'(a1,i3,1x,a1,i3)') 'P',i,'F',j
               do k=0,size 
                  place=best_nn(size,i,j)+k    
                  if (chain(place).eq.' ') chain(place)='_' 
                  write(iunit,1000)
     #                 frag_name(place),chain(place), 
     #                 resseq(place),res_list(place),ss(place),  
     #                 phi_list(place),psi_list(place),
     #                 omega_list(place), score,dme,dme_f,   
     #                 best_nn_ss_type(size,i,j),dipolar+noe,tag_string 
              enddo            ! k frags / residue
               write(iunit,*)  
            enddo
         enddo                  ! i residues
         close(iunit) 
      enddo                     !q sizes

      return 
 1000 format(1x,a4,1x,a1,1x,a5,2(1x,a1),5(f9.3),1x,f9.3,i2,1x,f9.3,1x,
     #     a9)
 2000 format(' position: ',i12,' neighbors:   ',i10)     
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

 

