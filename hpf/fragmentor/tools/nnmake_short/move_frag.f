c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 1.10 $
c  $Date: 2002/06/18 04:41:53 $
c  $Author: rohl $

      subroutine move_frag_set_quota(quota)
      implicit none
      include 'structure.h'

c *** inputs
      real quota(max_ss_type)
car local   
      integer i
      real sum

car private common 
      integer inserts_by_ss_type(max_ss_type)
      integer pointer_by_ss_type(max_ss_type)
      integer inserts_by_ssH, inserts_by_ssE, inserts_by_ssL
      integer total_inserts
      real ss_type_quota(max_ss_type)

      common/move_frag_static/inserts_by_ss_type, pointer_by_ss_type,
     #     total_inserts, ss_type_quota, inserts_by_ssH, inserts_by_ssE, 
     #     inserts_by_ssL
c     data ss_type_quota/0.10,0.35,0.45,0.25,0.0,0.0/ !cems  
c     data ss_type_quota/0.0,0.3,0.6,0.3,0.0,0.0/ !bonneau, casp4  
c     data ss_type_quota/0.10,0.3,0.45,0.3,0.3,0.0/ !add jufo prediction
      data ss_type_quota/0.00,0.2,0.6,0.2,0.2,0.0,0.0,0.0,0.0/
     
      sum=0.0
      write(0,*)"ss type quotas:"
      do i=1,max_ss_type
         ss_type_quota(i) = quota(i)
         sum=sum+ss_type_quota(i)
         write(0,*) i, ss_type_quota(i)
      enddo
      if (sum.ge.1.0) return
      write(0,*) 'frag picking quotas must sum to 1.0 or preferably larger'
      stop
      end
c---------------------------------------------------------------------
     
      real function move_frag_get_quota(type)
      implicit none
      include 'structure.h'

      integer type

car private common 
      integer inserts_by_ss_type(max_ss_type)
      integer pointer_by_ss_type(max_ss_type)
      integer inserts_by_ssH, inserts_by_ssE, inserts_by_ssL
      integer total_inserts
      real ss_type_quota(max_ss_type)

      common/move_frag_static/inserts_by_ss_type, pointer_by_ss_type,
     #     total_inserts, ss_type_quota, inserts_by_ssH, inserts_by_ssE, 
     #     inserts_by_ssL

      move_frag_get_quota=ss_type_quota(type)
      return
      end
     
c---------------------------------------------------------------------

      subroutine move_frag(bb_nn,b_nn,b_nn_ss_type)                      
c frags are picked in round robin (starting with jones), with lose-your-turn 
c if over your quota cap or if you are redundant to previous pick.
c this strategy prevents an ss-type from contributing (approx) more than 
c its cap. quotas should sum to slightly more than one since they 
c are quota-caps not quotas.

car takes frags from big_best_nn (bb_nn) and assembles best_nn (b_nn)

car:  note cems does 3:1 interleave of nmr_type to all other types
car for final ratio of 70:30 see frag_merg
      implicit none
      include 'structure.h' 

      integer bb_nn(max_ss_type,max_nn,max_res,max_len)
! intended to be an implicit equivalence to big_best_nn
      integer b_nn(max_len,max_res,max_nn) !generally best_nn
      integer b_nn_ss_type(max_len,max_res,max_nn)  !generally best_nn_sstype  
 
      integer size,i,j,q 
      integer ss_type
      logical fragments_exist
      logical reset_done

      integer inserts_by_ss_type(max_ss_type)
      integer pointer_by_ss_type(max_ss_type)
      integer inserts_by_ssH, inserts_by_ssE, inserts_by_ssL
      integer total_inserts
      real ss_type_quota(max_ss_type)

      common/move_frag_static/inserts_by_ss_type, pointer_by_ss_type,
     #     total_inserts, ss_type_quota, inserts_by_ssH, inserts_by_ssE, 
     #     inserts_by_ssL

c frags are picked in round robin (starting with jones), with lose-your-turn 
c if over your quota cap or if you are redundant to previous pick.
c this strategy prevents an ss-type from contributing (approx) more than 
c its cap. quotas should sum to slightly more than one since they 
c are quota-caps not quotas.

c-------------------------------------------------------	 

      do q = 1,n_sizes
         size=sizes(q)
         do i = 1,total_residue-size

cmj   reset parameters for next fragment

            fragments_exist=.false.
            reset_done = .false.
            total_inserts=0
            inserts_by_ssH = 0
            inserts_by_ssE = 0
            inserts_by_ssL = 0
            ss_type_quota(helix_type) = 0.0
            ss_type_quota(sheet_type) = 0.0
            ss_type_quota(loop_type)  = 0.0
            do j=1,max_ss_type
               inserts_by_ss_type(j) = 0
               pointer_by_ss_type(j) = 1
            enddo
            ss_type=jones_type          ! start process with jones frag
          
cmj   insert max_nnfragments

            do while(total_inserts.lt.max_nn)

cmj   check if fragment list is not emty. if so, set fragments_exist to be true 
cmj   if fragments_exist is false at the end of the ss_type loop the quota for the best
cmj   helix, sheet and loop fragments is set to be 1.0 (see before)
cmj   call insertion if list is not emty and quota not equal 0.0

cmj            if(size.eq.8 .and. i.eq.1) write(0,*) ss_type, 
cmj  #              ss_type_quota(ss_type)

               if(ss_type_quota(ss_type).gt.0.0.and.
     #           pointer_by_ss_type(ss_type).lt.max_nn) then
                 fragments_exist = .true.
                 call place_frag(bb_nn,size,ss_type,i,b_nn,b_nn_ss_type)
               endif

cmj   switch to next ss_type

               ss_type = ss_type+1

cmj   if ss_type greater then max_ss_type reset to 1
cmj   if no fragments were found in last iteration, set quotas for best
cmj   helix, sheet and loop fragments to be 1 so that 3x200 new fragments exist
cmj   set fragments_exist to false => it will be set true if one of the lists is not empty

               if (ss_type.gt.max_ss_type) then
                  if(.not.fragments_exist) then
                     ss_type_quota(helix_type) = 1.0
                     ss_type_quota(sheet_type) = 1.0
                     ss_type_quota(loop_type)  = 1.0
                  endif
                  ss_type=1
                  fragments_exist=.false.
               endif

cmj   reset fragment lists for 25..max_nn fragments 

               if (total_inserts.eq.25.and. .not.reset_done) then
                  fragments_exist=.true.
                  reset_done = .true.
                  ss_type_quota(helix_type) = 0.0
                  ss_type_quota(sheet_type) = 0.0
                  ss_type_quota(loop_type)  = 0.0
                  do j=1,max_ss_type 
                     pointer_by_ss_type(j) = 1
                  enddo
               endif
            enddo               ! ss_type/total_inserts next ss_type and next insertion
         enddo                  ! i next residue position
      enddo                     ! q next size
      return
      end


c---------------------------------------------------------------------------
      subroutine place_frag(bb_nn,size,ss_type,i,b_nn,b_nn_ss_type)
      implicit none
      include 'structure.h'
      integer bb_nn(max_ss_type,max_nn,max_res,max_len)
      integer size,ss_type
      integer i

      integer b_nn(max_len,max_res,max_nn)          !generally best_nn
      integer b_nn_ss_type(max_len,max_res,max_nn)  !generally best_nn_sstype  

      integer k,vall_frag
         
      integer inserts_by_ss_type(max_ss_type)
      integer pointer_by_ss_type(max_ss_type)
      integer inserts_by_ssH, inserts_by_ssE, inserts_by_ssL
      integer total_inserts
      real ss_type_quota(max_ss_type)

      common/move_frag_static/inserts_by_ss_type, pointer_by_ss_type,
     #     total_inserts, ss_type_quota, inserts_by_ssH, inserts_by_ssE, 
     #     inserts_by_ssL

      integer quota_number

cmj   this number ensures that the quotas for ss_type and the secondary structure
cmj   are fulfilled within the first 25 and also within all max_nn fragments 

      if(total_inserts.lt.25) then
         quota_number = 25
      else
         quota_number = max_nn
      endif

cmj   check to see if inserting this frag would exceed its overall quota
cmj   if so, set pointer_by_ss_type to max_nn => list is not used
  
      if (real(inserts_by_ss_type(ss_type)).ge.
     #     (ss_type_quota(ss_type)*quota_number)) then
         pointer_by_ss_type(ss_type) = max_nn
         goto 102
      endif

      vall_frag = bb_nn(ss_type,pointer_by_ss_type(ss_type),i,size)

cmj   check to see if the overall ratio of helix, sheet and loop  is
cmj   already fullfilled for the actual secondary type => if so skip it

      if (ss(vall_frag+size/2).eq.'H'.and.real(inserts_by_ssH).gt.
     #     (rel_aveH(i+size/2)*quota_number)) goto 101 
      if (ss(vall_frag+size/2).eq.'E'.and.real(inserts_by_ssE).gt.
     #     (rel_aveE(i+size/2)*quota_number)) goto 101
      if (ss(vall_frag+size/2).eq.'L'.and.real(inserts_by_ssL).gt.
     #     (rel_aveL(i+size/2)*quota_number)) goto 101

cmj   check if the fragment was already inserted => if so skip it

      do k = 1,total_inserts,1                         
         if (b_nn(size,i,k).eq.vall_frag) goto 101
      enddo

cmj   fragment is accepted for insertion into list

      total_inserts = total_inserts+1 
      b_nn(size,i,total_inserts) = vall_frag

cmj   keep track of ss_type and number of insertions from this ss_type

      b_nn_ss_type(size,i,total_inserts) = ss_type
      inserts_by_ss_type(ss_type) = inserts_by_ss_type(ss_type)+1

cmj   keep track of the secondary structure of the middle position

      if (ss(vall_frag+size/2).eq.'H') inserts_by_ssH = inserts_by_ssH+1
      if (ss(vall_frag+size/2).eq.'E') inserts_by_ssE = inserts_by_ssE+1
      if (ss(vall_frag+size/2).eq.'L') inserts_by_ssL = inserts_by_ssL+1

cmj   safety check - this should never happen

      if (inserts_by_ss_type(ss_type).gt.max_nn) then
         write(5,*) "ERROR encountered in move_frag: probably mis-allocated ss_type_quotas in move_frag.f"
         stop
      endif

cmj   go to next fragment in the frgment list

 101  continue
      pointer_by_ss_type(ss_type)=pointer_by_ss_type(ss_type)+1

cmj   end
      
 102  continue
      return
      end



c--------------------------------------------------------------------------

      subroutine merge_frag_lists(list1,sstype1,ratio1,list2,sstype2,
     #     ratio2,nsizes, size_list,max_len,nres,max_res,max_nn,
     #     outlist,outsstype)

      implicit none

      integer max_len,max_res,max_nn
      integer list1(max_len,max_res,max_nn)  
      integer list2(max_len,max_res,max_nn)
      integer sstype1(max_len,max_res,max_nn)
      integer sstype2(max_len,max_res,max_nn)
      integer ratio1,ratio2      !relative amounts from the 2 lists
      integer nsizes
      integer size_list(nsizes)
      integer nres 
      integer outlist(max_len,max_res,max_nn) 
      integer outsstype(max_len,max_res,max_nn) 

car local
      integer ptr1,ptr2,ptr3     !ptrs into list1,list2
      integer next_sstype
      integer next_frag
      integer size
      integer q,position,k
      integer which_list
      integer from_current_list

      do q=1,nsizes
         size=size_list(q)
         do position=1,nres-size
          
            ptr1=1
            ptr2=1
            ptr3=0
            from_current_list=0
            which_list=1
            
            do while (ptr3.lt.max_nn) 
car get next frag from list
               if (which_list.eq.1) then
                  next_frag=list1(size,position,ptr1)
                  next_sstype=sstype1(size,position,ptr1)
                  ptr1=ptr1+1   !don't look at this one again
               else
                  next_frag=list2(size,position,ptr2)
                  next_sstype=sstype2(size,position,ptr2)
                  ptr2=ptr2+1   !don't look at this one again
               endif
               if (next_frag.eq.0) then
                  write(0,*)"WARNING!!! There is a null fragment in "//
     #                 "in the fraglist"
                  write(0,*)'frag0',which_list,size,position,ptr1,ptr2
               endif
car check for redundant frags, but note that frag still counts toward
car from_current_list counter:
               from_current_list=from_current_list+1
               do k = 1,ptr3    !check against current frags
                  if (outlist(size,position,k).eq.next_frag) then
                     goto 101   ! already been inserted, skip value
                  endif
               enddo
               ptr3=ptr3 +1     !add to the list
               outlist(size,position,ptr3)=next_frag
               outsstype(size,position,ptr3)=next_sstype

car check if we need to switch input lists               
 101           continue       !frag rejected
               if (which_list.eq.1) then
                  if (from_current_list.ge.ratio1) then
                     which_list=2
                     from_current_list=0
                  endif
               else
                  if (from_current_list.ge.ratio2) then
                     which_list=1
                     from_current_list=0
                  endif
               endif
car check if there are more frags on the current list:
               if (which_list.eq.1 .and. ptr1.gt.max_nn) then
                  if (ptr2.gt.max_nn) goto 666
                  which_list=2
               elseif (which_list.eq.2 .and. ptr2.gt.max_nn) then
                  if (ptr1.gt.max_nn) goto 666
                  which_list=1
               endif

            enddo   !while less than max_nn
 105        continue
         enddo   !next position
              
      enddo    !next size

      return
 666  write(0,*)'WARNING!!! not enough frags to fill merged list!'
      write(0,*)'ptr1: ', ptr1,' ptr2: ',ptr2,' ptr3: ',ptr3
      write(0,*)'position: ', position, 'size: ',size
      goto 105

      end

c------------------------------------------------------------------------
      subroutine copy_bbnn_to_fraglist(bb_nn,ss_type,max_ss_type,max_nn,
     #     max_res,max_len,list,sstype_list)

      implicit none

car input
      integer max_ss_type,max_nn,max_res,max_len
      integer bb_nn(max_ss_type,max_nn,max_res,max_len)
      integer ss_type

car output
      integer list(max_len,max_res,max_nn)
      integer sstype_list(max_len,max_res,max_nn)

car local
      integer i,j,k

      do i=1,max_len
         do j=1,max_res
            do k = 1,max_nn
               list(i,j,k)= bb_nn(ss_type,k,j,i)
               sstype_list(i,j,k)=ss_type
            enddo
         enddo
      enddo

      return
      end
