c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 1.7 $
c  $Date: 2002/06/18 04:41:53 $
c  $Author: rohl $

      subroutine makeoutput(tag)    
c------------------------------------------------------
c  THIS SUBROUTINE WILL WRiTE  THE *.NN FILE 
c------------------------------------------------------ 
      implicit none 
      include 'path_defs.h'
      include 'structure.h' 

car  tag_dummy should be a reflections of the quota_cap mixture?

      integer tag

      integer i,j,k,q,size,place     
      real dme,dme_native,dme_f,dme_full,score,score_fragment,zsm ,noe
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
               zsm = -(score-
     #              z_score_mean(i,q,1+best_nn_ss_type(size,i,j)))/
     #              z_score_deviation(i,q,1+best_nn_ss_type(size,i,j))
               dme=dme_native(i,place,size) 
               dme_f= dme_full(i,place,size)
               noe = test_frag_noe(place)
               dipolar= retrieve_frag_dipolar_score(place)
               write(tag_string,'(a1,i3,1x,a1,i3)') 'P',i,'F',j
               do k=0,size 
                  place=best_nn(size,i,j)+k    
                  if (chain(place).eq.' ') chain(place)='_' 
cRB               write(iunit,1000)
cRB  #                 frag_name(place),chain(place), 
cRB  #                 resseq(place),res_list(place),ss(place),  
cRB  #                 phi_list(place),psi_list(place),
cRB  #                 omega_list(place), score,dme,dme_f,   
cRB  #                 best_nn_ss_type(size,i,j),dipolar+noe,tag_string 
                  write(iunit,1050)
     #                 frag_name(place),chain(place),
     #                 resseq(place),res_list(place),ss(place),
     #                 phi_list(place),psi_list(place),
     #                 omega_list(place)
              enddo            ! k frags / residue
               write(iunit,*)  
            enddo
         enddo                  ! i residues
         close(iunit) 
      enddo                     !q sizes

      return 
 1000 format(1x,a4,1x,a1,1x,a5,2(1x,a1),5(f9.3),1x,f9.3,i2,1x,f9.3,1x,
     #     a9)
 1050 format(1x,a4,1x,a1,1x,a5,2(1x,a1),3(f9.3),1x)
 2000 format(' position: ',i12,' neighbors:   ',i10)     
      end 
 

