c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 1.7 $
c  $Date: 2002/06/18 04:41:53 $
c  $Author: rohl $


      subroutine find_nn(isize)
      implicit none
c-------------------------------------------------------
c     MATCH SEGMENTS TO PROFILE
c-------------------------------------------------------
      include 'structure.h'

c input
      integer isize             ! actual size not the fake size-1
      real proline_phi_score
      logical bad
c local
      integer i,j,k
c heap storage
      integer max_heap
      parameter (max_heap =max_nn)
      integer heap_index(max_heap+2,max_ss_type)
      real heap_value(max_heap+2,max_ss_type)
c ss_score storage
      real ss_scores(max_ss_type)

c some temporary and handshaking variables for hash
      logical h_err
      real h_coval
      integer h_val,h_index
      real base_score,seq_score

      ss_scores(seq_type) = 0.0 ! initialize this to zero

      do i=1,total_residue-isize+1   ! all frag positions

c make a list of noe  data to score, 
c dipolar: score all frags to get chi distribution
         call prepare_constraints_score(i,isize,vall_size)

c initialize heap
         do k = 1,max_ss_type
            call heap_empty(heap_index(1,k),heap_value(1,k),max_heap)
         enddo

         do j=1,vall_size-isize+1 !all vall fragments

            call  all_score(i,j,isize,
     #           seq_score,ss_scores,proline_phi_score,bad)

            if (.not.bad) then  ! flags gaps in vall plus cis omegas
               base_score = seq_score*seq_weight
     #              +proline_phi_score*weight_proline_phi

               
               do k =1,max_ss_type      !store all score in heap
                  call heap_insert(heap_index(1,k),heap_value(1,k),j,
     #                 -(base_score+ss_weight*ss_scores(k)),h_err)
               enddo            ! k (scores)
            endif
         enddo                  ! j (vall)

c extract the contents of all the heaps (in sorted order)
c store in big_best_nn for later merging
cems  I note that there is no reason we have to store all of these frags 
cems  before writing them out; we could just merge and write out results
cems  for each frag position as we generate them

         do k = 1,max_ss_type    ! all the heaps
            call heap_get_size(heap_index(1,k),h_val)
            if (h_val.lt.max_nn) then
               write(0,*) "find_nn: insufficient fragments at position",
     #              i," have ",h_val," but wanted",max_nn
               stop
            endif
            
            do j = 1,h_val
               call heap_extract(heap_index(1,k),heap_value(1,k),
     #              h_index,h_coval,h_err)
car save in reverse order
               big_best_nn(k,h_val-j+1,i,isize-1)=h_index
             enddo
          enddo          
      enddo                     ! i (residue position)
      write(*,*) "find nn done"
      return 
      end 

c----------------------------------------------------------------------------
c----------------------------------------------------------------------------
c     charlie strauss 1999,2000

c     general purpose heap manager.  heaps are a very efficient way of
c     sorting.  more importantly they are the most efficient for the
c     problem of finding the TOP-N of a long list.  they do a just in
c     time sorting scheme that does not actually sort the list until it
c     is actually required to be in sorted order.  and when seeking the
c     top N of something it doe no more sorting than minimially required
c     for finding the top-N.

c     heap and coheap are a pair of arrays; coheap contains the values on
c     which you want things sorted; heap contains keys associated with
c     these value and is sorted identicaally to coheap.  currently heap
c     is heap is an arra of integers with two extra storage elments at
c     the front-- the first is the max dimension of the heap, the second
c     is the current number of entries, the third is the the start of
c     the heap integers (and the key to the smallest value in coheap).

c     variables 
c           integer max_heap = 2000 
c           integer heap(max_heap+2) ! note the +2 VERY VERY IMPORTANT !!!!!  
c           real coheap(max_heap) !note: +2 is not required here 
c     WARNING: heap NEVER checks array ounds even if compiler is set to check.

c     initialize the heap:  call heap_empty(heap,coheap,max_heap) 

c     insert a value, repeat this many times:
c             callheap_insert(heap,coheap,val,coval,err) 
c     this places the val and coval into the heap; when the heap
c     contains max_heap items it is full; at this point new items
c     inserted bump out old items -- the LOWEST value in the heap is
c     discarded to make room for the new item.
c
c     extraction in sorted order from lowest to highest call
c     heap_extract(heap,coheap,val,coval,err) each call extracts one

c     charlie strauss 1999,2000
c----------------------------------------------------------------------        
      subroutine heap_empty(heap,coheap,max_size)
        
c creates (or erases) a heap with no entries
c and sets maximum size of max_size

      implicit none
      integer heap(-2:*),max_size
      real coheap(0:*)
      
      heap(-1) = 0
      heap(-2) = max_size
        
      return
      end

c----------------------------------------------------------------------        

      subroutine heap_extract(heap,coheap,val,coval,err)
      implicit none
c modifes heap and last_val
c return val and err.          
      integer val
      logical err
      integer heap(-2:*)        ! convert to zero offset matrix
      real coval,coheap(0:*)
      
      err = .true.
      if (heap(-1).lt. 1 ) return
      
      val = heap(0)
      coval = coheap(0)
      err = .false.
      heap(-1) = heap(-1)-1
       
      
      if (heap(-1) .eq. 0) return ! special case for singel value in heap
      
      heap(0) = heap(heap(-1))  ! move last value to front
      coheap(0) = coheap(heap(-1))
      
      call heap_down(heap,coheap,1)
      
      return
      end
 
c----------------------------------------------------------------------        
      subroutine heap_insert(heap,coheap,val,coval,err)
      implicit none
      
c modifes heap and last_dummy, inserts val, returns err
c requires heap_max to be previously initialized
          
      integer val
      logical err
      integer heap(-2:*)        ! convert to zero offset matrix
      real coheap(0:*),coval
      
      if (heap(-1).ge. heap(-2)) then ! list is full, use replace instead
         err = .true.
	 
         if (coheap(0).lt.coval) call heap_replace(heap,coheap,val,coval) 
         return
      endif
c     write(0,*) "inserting",heap(-1),val,coval
      err = .false.
      heap(heap(-1)) = val      ! empty spot on end (zero offset)
      coheap(heap(-1)) = coval
      
      heap(-1) = heap(-1)+1
      call heap_up(heap,coheap,heap(-1))
      
      return
      end
c----------------------------------------------------------------------        
      
      subroutine heap_replace(heap,coheap,val,coval)
! overwrite the lowest element 
      implicit none
c modifes heap  
          
      integer val 
      real coval, coheap(0:*) 
      logical err
      integer heap(-2:*)        ! convert to zero offset matrix
      err = .false.
c	write(0,*) 'replacing',val,coval,heap(0),coheap(0) 
      heap(0) = val             ! overwrite the lowest element
      coheap(0) = coval 
      call heap_down(heap,coheap,1)
      return
      end

c----------------------------------------------------------------------        
      subroutine heap_get_size(heap,size)
      implicit none
      integer heap(-2:*)     
      integer size
      size = heap(-1)
      return
      end
c----------------------------------------------------------------------        

      subroutine heap_down( heap,coheap,index_dummy )
      implicit none
      integer index_dummy 
      integer heap(-2:*)        ! convert to zero offset matrix
      real coheap(0:*),coiv,cocv,cocv2
      integer index,child,iv,cv,cv2,last
      index = index_dummy-1     ! convert to zero offset matrix
      last  =  heap(-1)-1       !convert to zero offset matrix
      
      if (last .le. 0 )  return ! empty or single element
      if (index.gt.last) return ! dumbass
      
      iv = heap(index)          ! the inserted value
      coiv = coheap(index)
      
      do while (index .lt. last)
         child = 2*index+1
         
         if (child .gt. last ) goto 20 ! loop escape
         
         cv  = heap(child)
         cocv =  coheap(child)
         
         if (child .lt. last ) then
            cv2 = heap (child+1)
            cocv2 = coheap(child+1)
            
            if ( cocv2 .lt. cocv ) then
               cv = cv2
               cocv = cocv2
               
               child = child+1
            endif
         endif
         
         if (coiv .le. cocv ) goto 20 ! loop escape
         coheap(index) = cocv
         heap(index) = cv
         index = child
      enddo
      
 20   continue                  ! loop escape
      heap(index) = iv
      coheap(index) = coiv
      return
      end
      
c----------------------------------------------------------------------        

      subroutine heap_up(heap,coheap,index_dummy)
      implicit none
      integer index_dummy
      integer heap(-2:*)        ! convert to zero offset matrix
      real coheap(0:*),covalue,copv
      
      integer index,parent,value,pv  
      index = index_dummy-1     ! convert to zero offset matrix
      

      value = heap(index)
      covalue=coheap(index)
      
      do while (index .ne. 0) 
         parent = int ((index-1)/2)
         pv = heap(parent)
         copv = coheap(parent)
         if (copv .lt. covalue)  goto 20 ! loop escape
         coheap(index) = copv
         heap(index) = pv
         index = parent
      enddo
      
 20   continue                  ! loop escape
      coheap(index) = covalue
      heap(index) = value
      return
      end
         
 




