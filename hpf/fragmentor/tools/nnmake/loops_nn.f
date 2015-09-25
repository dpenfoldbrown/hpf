c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 4051 $
c  $Date: 2004-02-26 16:30:08 -0500 (Thu, 26 Feb 2004) $
c  $Author: jeff $

      subroutine setup_template(zone_file,loop_list,nloop)
      
      implicit none
      include 'param.h'
      include 'structure.h'
      include 'path_defs.h'

      character*(*) zone_file

      integer nloop,i
      integer loop_list(5,max_loops)
      
c package globals
      integer max_temp_res
      parameter (max_temp_res=1000)
      real temp_cnitro(3,max_temp_res),temp_calpha(3,max_temp_res)
      real temp_ccarbon(3,max_temp_res),temp_coxygen(3,max_temp_res)
      real temp_cbeta(3,max_temp_res)
      character*1 template_ss(max_temp_res)
      integer template_residue
      common/template_atoms/temp_cnitro,temp_calpha,temp_ccarbon,
     #     temp_coxygen,temp_cbeta,template_ss,template_residue   

car local
      integer iunit
      character*1000 filename

      data iunit /72/

      filename=loop_path(1:loop_path_l)
     #     //zone_file(1:index(zone_file,' ')-1)//'.zones'
      call read_zone_file(filename,iunit,loop_list,nloop,max_loops,
     #     total_residue)

      filename=loop_path(1:loop_path_l)
     #     //zone_file(1:index(zone_file,' ')-1)//'.pdb'
      call input_pdb_template(filename,iunit,temp_cnitro,temp_calpha,
     #     temp_ccarbon,temp_coxygen,temp_cbeta,max_res,template_residue) 
      do i=1,template_residue
         if (seqnum(i).eq.6) call compute_CB_coord(temp_cnitro(1,i),
     #        temp_calpha(1,i),temp_ccarbon(1,i),temp_cbeta(1,i))
      enddo

      filename=loop_path(1:loop_path_l)//
     #     zone_file(1:index(zone_file,' ')-1)//'.ssa'
      call read_ssa(filename,iunit,template_residue,template_ss)
      return
      end

c---------------------------------------------------------------------
      subroutine read_zone_file(filename,iunit,loop_list,nloop,
     #     max_loops,total_residue)

      implicit none
car input
      character*(*) filename
      integer iunit
      integer total_residue
      integer max_loops
c output
      integer loop_list(5,max_loops),nloop

car local
      integer i,j,k
      character*1000 line

car note charlie numbered loops s.t loop_begin is the last aligned residue
car and loop_end is the first aligned residue of the next aligned segment

car zone format
car ZONE a-b:c-d      residues a to b of query align with c to d of parent
      nloop=0
      write(0,*)'.zones file: ',filename(1:index(filename,' '))
      open (unit=iunit,file=filename,status='OLD')

 100  read(iunit,1000,end=999)line
      if (line(1:4).ne.'ZONE' .and. line(1:4).ne.'zone') goto 100
      nloop=nloop+1
      if (nloop.gt.max_loops) then
         write(0,*)'increase max_loops in read_ksync_file'
         write(0,*)'likely this is max_loops in structure.h'
         stop
      endif
      if (nloop.gt.1) then
         read(line,2000)loop_list(3,nloop-1),loop_list(2,nloop),
     #        loop_list(5,nloop-1),loop_list(4,nloop)
      else
         read(line,2000)i,loop_list(2,nloop),j,
     #        loop_list(4,nloop)
         if (i.ne.1) then
            loop_list(2,2)=loop_list(2,1)
            loop_list(4,2)=loop_list(4,1)
            loop_list(2,1)=0
            loop_list(3,1)=i
            loop_list(4,1)=j-loop_list(2,2)+loop_list(3,1)
            loop_list(5,1)=j
            nloop=nloop+1
         endif
      endif
      goto 100
      
 999  continue
      if (loop_list(2,nloop).ne.total_residue) then
         loop_list(3,nloop)=total_residue+1
         loop_list(5,nloop)=loop_list(3,nloop)-
     #     loop_list(2,nloop)+loop_list(4,nloop)
      else
         nloop=nloop-1
      endif

      write(0,*) 'length, query start, query end, parent start, parent end'
      do i=1,nloop
         loop_list(1,i)=loop_list(3,i)-loop_list(2,i)+1
         write(0,*)(loop_list(k,i),k=1,5)
      enddo
      
      return
 1000 format(1000a)
 2000 format(5x,i4,1x,i4,1x,i4,1x,i4)
      end

c-----------------------------------------------------------------------------
    
      subroutine read_ssa(filename,iunit,total_residue,ss)

      implicit none

      character*(*) filename
      integer iunit
      integer total_residue
      character*(1) ss(total_residue)

car local
      integer i
      integer length
      character*1000 linebuffer
      
      write(0,*)'Looking for secondary structure assignment file:',
     #     filename(1:index(filename,' '))
      open(iunit,file=filename,status='old',iostat=i)

      if (i.ne.0) then 
         write(0,*) 'Warning!! ssa file not found'
         return
      endif

      read(iunit,200)linebuffer
      length=index(linebuffer,' ')-1
      if (length.ne.total_residue) then
         write(0,*)'Disagreement in chain length in ssa file'
         write(0,*)'total residue', total_residue,'ssa length', length
         stop
      endif

      close(iunit)
    
      do i=1,total_residue
         ss(i)=linebuffer(i:i)
         if (ss(i).eq.'C') ss(i)= 'L'
         if (ss(i).ne.'L' .and. ss(i).ne.'H' .and. ss(i).ne.'E') 
     #        goto 150
      enddo

      return

 150  write(0,*)'Unrecognized character in ssa file: ',ss(i)  
      stop

 200  format(a1000)
      end


c---------------------------------------------------------------------
      subroutine find_loops(zone_file,loop_list,N_loop_alignments)

      implicit none

c     MATCH SEGMENTS TO PROFILE
car  NOTES::  phd scored but not used;  quotas should be
car  checked to see if ss predictions were read in correctly
car  should add checks for cis omegas and bad proline phi angles

car note no constraint info is ever used

      include 'param.h'
      include 'structure.h'
      include 'path_defs.h'

      integer funit
      parameter (funit=22)

c     input
      integer n_loop_alignments
      integer loop_list(5,N_loop_alignments)
      character*(100) zone_file
      
c     local
      integer size
      integer i,ii,j
  
      real diff
      real rms_score,seq_score
      integer loop_size,loop_begin,loop_end
      real dme,dme_native
      real phd_score,jones_score,rdb_score,jufo_score,template_ss_score
      real splice_rms
      character*(600) filename
      logical terminal,good

c heap storage
      integer heap_index(max_depth+2)
      real heap_value(max_depth+2)
      integer top_loops(max_depth)  !array to reverse order of heap
c some temporary and handshaking variables for heap
      logical heap_err
      integer next    !next frag in heap
      real next_score
      integer heap_depth
     
c-------------------------------------------------------

car open files:
      filename=make_output_path(1:make_output_l)//
     #     zone_file(1:index(zone_file,' ')-1)//'.loops_all'
      write(0,*)'loop file: ',filename(1:index(filename,' '))
      open(unit=funit, file=filename)

      do ii=1, N_loop_alignments !  ii= loop number 
         loop_size = loop_list(1,ii)   
	 loop_begin = loop_list(2,ii) !last aligned
	 loop_end =   loop_list(3,ii) !first aligned
         size = loop_end-loop_begin ! dont forget:  insane definition of size
                                !nnmake numbers frags or loops from 0->'size'
         i = 1                  ! i is query position, fix to 1 here cause we
                                !write out the loop, each loop overwrites last
         terminal=.false.
         if (loop_begin.eq.0 .or. loop_end.eq.total_residue+1) terminal=.true.

         if (size.gt.max_loop_length) then
            call write_empty_loop(funit,loop_begin+1,loop_end-1)
            goto 200
         endif
c display target functions
c$$$         write(0,*) "loop_param ",ii,(loop_list(j,ii),j=1,5)
c$$$         write(0,*)'loop_sequence ',lower_name, loop_size-2,(seq(j),j=loop_begin+1,loop_end-1)
c$$$	write(0,*) 'jones ss query loop ', ii, ':'
c$$$	do j = 0,size
c$$$	    write(0,*) rel_jonH(loop_begin+j), 
c$$$     #          rel_jonE(loop_begin+j) ,
c$$$     #          rel_jonL(loop_begin+j)  
c$$$        enddo

car setup the template for 
         if (.not.terminal) then
            call compute_loop_ss_gap(loop_list(1,ii))
            call set_rms_weight(1,1.0)
         endif

c initialize heap
         call heap_empty(heap_index,heap_value,max_depth)

         do j=1,vall_size-size !all vall fragments

            call score_loops(loop_begin,j,size,rms_score,
     #           seq_score,phd_score,jones_score,rdb_score,jufo_score,
     #           template_ss_score,good)

            if (good) then
CAR Here's the score that's actually used!!
               diff= seq_score+
     #              0.33*(jones_score+rdb_score+jufo_score)+
     #              template_ss_score+0.05*rms_score
               if (loop_size.LE.14) then
                  diff = diff+0.05*rms_score
                  if (rms_score.gt.1.4) then
                     diff = diff+rms_score 
                     if (rms_score.gt.2.4) diff = diff+rms_score 
                  endif
               endif  
               if (loop_size.le.6) diff = diff+rms_score
c  store in heap, change sign of score cause heap keeps biggest values
              call heap_insert(heap_index,heap_value,j,-diff,heap_err)
            endif
         enddo                  ! j (vall)

c check to see that heap is fully populated
         call heap_get_size(heap_index,heap_depth)
         if (heap_depth.lt.max_depth) then
            write(0,*) "find_loops: insufficient loops for loop",
     #           ii," have ",heap_depth," but wanted",max_depth
            stop
         endif

c extract the contents of the heap (in sorted order),
c rescore them and write them out
         do j = 1,heap_depth
            call heap_extract(heap_index,heap_value,next,next_score,heap_err)
            top_loops(heap_depth-j+1) = next !reverse order
         enddo
         do j=1,heap_depth
            next=top_loops(j)
            call score_loops(loop_begin, next,size,rms_score,
     #           seq_score,phd_score,jones_score,rdb_score,jufo_score,
     #           template_ss_score,good)
            if (terminal) then
               splice_rms=0.0
            else
               call calc_splice_rms(xyz_vall(1,1,next),
     #              omega_list(next),phi_list(next+loop_size-1),
     #              loop_size-2,loop_begin+1,loop_end-1,splice_rms)
            endif
            dme=0.0
            if (use_pdb.eq.'Y') dme = dme_native(loop_begin,next,size)
c add 1 to loop position (write out only unaligned positions)
            call write_one(funit,loop_begin+1,loop_end-1,
     #           size-2,next+1,dme,rms_score,splice_rms,
     #           seq_score,template_ss_score,jufo_score,
     #           jones_score,rdb_score)       
         enddo
 200     continue              !escape for longer loops
      enddo                     ! next loop region (i)
      return 
      end 

c----------------------------------------------------------------------
      
      subroutine write_one(funit,loop_begin,loop_end,size,v_point, 
     #     dme_score,rms_score,splice_score,
     #     seq_score,template_ss_score,jufo_score,jones_score,rdb_score)
      
      implicit none
c     wrapper to hide vall common block
      include 'param.h'
      include 'structure.h'
      
c     input      
      integer funit
      integer loop_begin,loop_end,size,v_point
      real dme_score,rms_score,splice_score
      real jufo_score,jones_score,rdb_score,template_ss_score,seq_score
      
c     local
      
      call write_loop(funit,loop_begin,loop_end,size,
     #     phi_list(v_point),psi_list(v_point),omega_list(v_point),
     #     ss(v_point),frag_name(v_point),dme_score,rms_score,
     #     splice_score,template_ss_score,jufo_score,jones_score,
     #     rdb_score,seq_score,v_point)
      return
      end
      
c------------------------------------------------------------------------      
      subroutine write_empty_loop(funit,loop_begin,loop_end)
      
      implicit none

car input
      integer funit,loop_begin,loop_end

car local
      integer size
      
      size=loop_end-loop_begin+1        
      write(funit,90)size,loop_begin,loop_end
      return
 90    format( i3,1x,I4,1x,I4,1x)
      
      end
      
c------------------------------------------------------------------------      
      subroutine write_loop(funit,loop_begin,loop_end,size,
     #     phi_list,psi_list,omega_list,ss,frag_name,dme_score,
     #     rms_score,splice_score,template_ss_score,jufo_score,
     #     jones_score,rdb_score,seq_score,v_point)
      implicit none
c     writes a loop fragment out to disk in standard format
c     inputs
      integer funit
      integer loop_begin,loop_end,size ! frag size and location
      real phi_list(size+1),psi_list(size+1),omega_list(size+1) ! frag angles
      character*1 ss(size+1)    ! ss array of frag residues
      character*(*) frag_name   ! the name of the frag
      real dme_score            ! the tru error in the frag 
      real rms_score,splice_score,template_ss_score
      real jufo_score,jones_score,rdb_score
      real seq_score ! frag scores
      integer v_point   ! source of ss_score, vall location
      
      
c     local
      
      integer i,j,ss_type
      integer isize
      

      character*1000   line

      ss_type=1
      isize = size+1            ! to restore my sanity
      write(line(1:14),90)isize,loop_begin,loop_end
      j=15
! ss
      do i=1,isize
         write(line(j:j),'(A1)')ss(i)
         j=j+1
      enddo
! scores
      write(line(j:j+102),91)dme_score,rms_score,splice_score,seq_score,
     #     template_ss_score,jufo_score,jones_score,rdb_score,
     #     ss_type,v_point,frag_name
      j=j+103
!angles
      do i=1,isize
         write(line(j:j+25),92)phi_list(i),psi_list(i),omega_list(i)
         j=j+26
      enddo
      
      write(funit,'(A)')line(1:j-1)

 90    format( i3,1x,I4,1x,I4,1x)
 91    format(1x,4(f9.4,1x),1x,4(f9.4,1x),1x,i2,1x,i9,3x,a4,1x)
 
 92    format(3(f7.2,1x),2x)
      
      return
      end
c--------------------------------------------------------------------
      subroutine score_loops(i_in,j,size_in,rms_score,
     #     seq_score,phd_score,jones_score,rdb_score,jufo_score,
     #     template_ss_score,good)
      implicit none

      include 'param.h'
      include 'structure.h'
      include 'dipolar.h'

CAR ADD PROLINE SCORES

c inputs
      integer i_in,j,size_in
c outputs
      real rms_score,seq_score,phd_score,jones_score,jufo_score
      real rdb_score,template_ss_score
      logical good
c local
      integer i,size
      integer k,aa
      real ss_score
      logical terminal
c---------------------------------------------------------
CCCC    passed in from find_nn: i-pos in query seq, 
CCCC        j-pos in vall , size= length of frag -1 (insane)
C---------------------------------------------------------

      good=.false.
      seq_score=0.0         
      ss_score=0.0
      phd_score=0.0      
      jones_score=0.0      
      rdb_score=0.0   
      jufo_score=0.0
      rms_score = 0.0
      template_ss_score = 0.0

      terminal=.false.
      if (i_in.eq.0) then
         i=1
         size=size_in-1
         terminal=.true.
      else
         i=i_in
      endif
      if (i_in+size_in.eq.total_residue+1) then
         size=size_in-1
         terminal=.true.
      else
         size=size_in
      endif
c     call score_cis_omega(size+1,omega,seq,cis)        
c     if (cis) return

      do k=0,size   
         if (phi_list(j+k).eq.0.0) return  !chain break
         
c  SEQUENCE MATCHING SCORE 
         do aa=1,20 
            seq_score=seq_score      
     #           +abs(vall_pro(aa,j+k)-profile(aa,i+k)) !city block difference
         enddo 
      enddo 
      
      do k = 0,size
C we sum the ss prediction reliabilities if they 
C for the actual ss at these positions in the vall
        
         if (ss(j+k).eq.'H') then
            phd_score = phd_score-rel_phdH(i+k) 
            jones_score = jones_score-rel_jonH(i+k)
            rdb_score = rdb_score-rel_rdbH(i+k) 
            jufo_score = jufo_score-rel_jufoH(i+k) 
            
         elseif(ss(j+k).eq.'E') then
            phd_score = phd_score-rel_phdE(i+k) 
            jones_score = jones_score-rel_jonE(i+k)
            rdb_score = rdb_score-rel_rdbE(i+k) 
            jufo_score = jufo_score-rel_jufoE(i+k) 
            
         elseif (ss(j+k).eq.'L') then 
            phd_score = phd_score-rel_phdL(i+k) 
            jones_score = jones_score-rel_jonL(i+k)
            rdb_score = rdb_score-rel_rdbL(i+k) 
            jufo_score = jufo_score-rel_jufoL(i+k) 
            
         endif
      enddo 
 

c compute the rms_score for the overlap of fragment with template
c note we will not reach this point if vall contains a bad spot 
c over this region  
      if (.not. terminal) then
         call score_rms(xyz_vall(1,1,j),rms_score)
         call score_template_ss(ss(j),template_ss_score)
      endif

      seq_score = seq_score*seq_weight
      phd_score = phd_score*ss_weight
      jones_score = jones_score*ss_weight
      rdb_score = rdb_score*ss_weight
      jufo_score = jufo_score*ss_weight
      template_ss_score = template_ss_score*ss_weight
      good=.true.

      return 
      end 


c-----------------------------------------------------------------------
      subroutine calc_splice_rms(vall,vall_omega,vall_phi,loop_size,
     #     loop_begin,loop_end,rms)
     

car here loop_begin and loop_end are the first and last UNaligned residues
cr and size is the number of unaligned residues
car vall vall_psi, and vall_omega start at loop_begin-1 (ie last aligned residue)
car vall_phi starts at loop_end+1 (ie first aligned residue)

      implicit none
      integer loop_size,loop_begin,loop_end
      real vall(3,5,*),vall_omega,vall_phi
c ouput
      real rms
       
c local
      integer k
         
      real temp_1(3,5),temp_2(3,5),temp_omega,temp_phi

car functions
      real dihedral

c package globals
      integer max_temp_res
      parameter (max_temp_res=1000)
      real temp_cnitro(3,max_temp_res),temp_calpha(3,max_temp_res)
      real temp_ccarbon(3,max_temp_res),temp_coxygen(3,max_temp_res)
      real temp_cbeta(3,max_temp_res)
      character*1 template_ss(max_temp_res)
      integer template_residue
      common/template_atoms/temp_cnitro,temp_calpha,temp_ccarbon,
     #     temp_coxygen,temp_cbeta,template_ss,template_residue


      temp_omega=  dihedral(              !omega of last aligned res
     #      temp_calpha(1,loop_begin-1),
     #      temp_ccarbon(1,loop_begin-1),
     #	    temp_cnitro(1,loop_begin),
     #      temp_calpha(1,loop_begin))
      temp_phi=  dihedral(                !phi of first aligned res
     #      temp_ccarbon(1,loop_end),
     #	    temp_cnitro(1,loop_end+1),
     #      temp_calpha(1,loop_end+1),
     #      temp_ccarbon(1,loop_end+1))

c copy into temp arrays so atoms are consecutive
      do k= 1,3
car 1/15/01: redefine which atoms rms is calculated over to be consistent
car with rosetta:
car will calc rms between temp_1 (template) and temp_2 (vall) after
car correcting for differences in angles

         temp_1(k,1) = temp_cnitro(k,loop_end+1) !first aligned res
         temp_1(k,2) = temp_calpha(k,loop_end+1)
         temp_1(k,3) = temp_ccarbon(k,loop_end+1)
         temp_1(k,4) = temp_coxygen(k,loop_end+1)
         temp_1(k,5) = temp_cbeta(k,loop_end+1)

         temp_2(k,1) = vall(k,1,loop_size+2)  !N   !first aligned res
         temp_2(k,2) = vall(k,2,loop_size+2)  !CA
         temp_2(k,3) = vall(k,4,loop_size+2)  !C !differ by omega(loop_size+1) 
         temp_2(k,4) = vall(k,5,loop_size+2)  !O
         temp_2(k,5) = vall(k,3,loop_size+2)  !CB                    
      enddo
      

      call get_splice_rms(temp_ccarbon(1,loop_begin-1),
     #     temp_cnitro(1,loop_begin),
     #     temp_calpha(1,loop_begin),temp_omega,temp_phi,
     #     temp_1,vall(1,4,1),vall(1,1,2),vall(1,2,2),vall_omega,
     #     vall_phi,temp_2,rms)
      return
      end
c----------------------------------------------------------------------
      subroutine get_splice_rms(p1,p2,p3,omega_p,phi_p,test_p,
     #     q1,q2,q3,omega_q,phi_q,test_q,rms)

c smoothly grafts q1,q2,q3 onto p1,p2,p3 after applying the difference
c in omega angles for the two residues.

c then computes new position of test_q  (correcting for diffence in phi_p
c and phi_q as well)
c  the position of test_q is then comared with test_p and the rms computed.
c effectively this computes what would happen if fragment were inserted and 
c refolded
c then it compares the downstream positions and compares those to template.

car *_p is from template,  *_q is from the vall
car p1,p2,p3 are C of last aligned residue and 
car              N,CA  of first unaligned residue of loop
car omega_p is the omega of the last aligned residue preceding loop
car phi_p is the phi of the first aligned residue following loop
car test_p is N,CA,C,O,CB of first aligned residue following loop

      implicit none
c inputs
      real p1(3),p2(3),p3(3),omega_p,phi_p
      real test_p(3,5)  
      real q1(3),q2(3),q3(3),omega_q,phi_q
      real test_q(3,5)   
c out
      real rms
	
c local
      integer i,k
      real matb(3,3),vecb(3),mata(3,3),veca(3),matc(3,3),vecc(3)
      real temp(3),q3_prime(3),domega,dphi
      real test_q_prime(3,5)

c  compute amount we need to rotate fragment to make it have an omega angle
c  of omega_p instead of omega_q	
      domega = omega_p-omega_q
	
c  compute rotation and offset that accomplishes this rotation	
      call getrot(q1,q2,domega,mata,veca)

c apply this rotation to q3 to determine where q3 would have been if 
c omega_q = omega_p	
      call rotate(mata,q3,temp)
      call addvec(temp,veca,q3_prime)

c and apply to test_q
      do i = 1,5
         call rotate(mata,test_q(1,i),temp)
         call addvec(temp,veca,test_q_prime(1,i))
      enddo

c now find rotation that lines up q1,q2,q3_prime with p1,p2,p3	
      call refold_align_transform(q1,q2,q3_prime,p1,p2,p3,Matb)

c find offset vector maping q1->p1	
      call rotate(Matb,q1,vecb)
      call subvec(p1,vecb,vecb)

c apply this rotation and offset to each test point in turn	
c to determine where test_q would be if q1,q2,q3' on p1,p2,p3
      do i = 1,5
         call rotate(matb,test_q_prime(1,i),temp)
         call addvec(temp,vecb,test_q_prime(1,i))
      enddo

c compute amout we need to rotate last CO to account for differences in phi
      dphi=phi_p-phi_q
c compute rotation and offset that accomplishes this rotation	
      call getrot(test_q_prime(1,1),test_q_prime(1,2),dphi,matc,vecc)
c apply to C,O   
      do i=3,4       
         call rotate(matc,test_q_prime(1,i),temp)
         call addvec(temp,vecc,test_q_prime(1,i))
      enddo
      rms = 0.0
      do i=1,5
         do k=1,3
            rms = rms+(test_p(k,i)-test_q_prime(k,i))**2
         enddo
      enddo
      rms = sqrt(rms/5.0)
      return
      end
      
c--------------------------------------------------------------------------
      subroutine getrot(a1,a2,chi,mat,vec)
c arguments are in degrees!!!!!!!!
      implicit none
** find the matrix (mat) and vector (vec) that
** rotate chi deg about an axis defined by a1-->a2
** rotation is right-handed. That is, it is a
** clockwise rotation when looking from a1 to a2
** chi in radians
      real a1(3),a2(3),vec(3)
      real mat(3,3)
      real a(3),c(3),da,phi,psi,chi  
      
      integer i,k
**
      if (abs(chi).lt.1.0E-4) goto 100
      phi=0.0
      psi=0.0
** get difference vector: a
      call subvec(a2,a1,a)
** get length of vector A
      da=a(1)*a(1) +a(2)*a(2) +a(3)*a(3) 
      da=sqrt(da)
** get phi, psi
      if (da.ne.0.0) psi=acos(a(3)/da)
      if (a(1).ne.0.0.or.a(2).ne.0.0) phi=atan2(-a(1),a(2))
        
** get matrix
      call getmat(phi,psi,chi,mat)
** get vector, = -Mat*a1 + a1
      call rotate(mat,a1,c)
      call subvec(a1,c,vec)
      return
 100  do I=1,3
         do K=1,3
            mat(K,I)=0.
         enddo
      enddo
      do I=1,3 
         mat(I,I)=1.
         vec(I)=0.
      enddo
      return
      end

c----------------------------------------------------------------------------
      subroutine getmat(phi,psi,kappa,aa)
c arguments are in degrees!!!!!!!!
      implicit none
** get matrix for Tanaka convention polar angles
** CB  14-JAN-1991, subroutine 3-JUN-1991

** >>> phi,psi,kappa in radians <<<
      real phi,psi,kappa,sf,cf,ss,cs,sk,ck,aa(3,3)
        
      integer i,j
****
cb Try to avoid underflows...
      if (abs(phi).lt.1.0E-4) phi = 0.0
      if (abs(psi).lt.1.0E-4) psi = 0.0
      if (abs(kappa).lt.1.0E-4) goto 100
      sf = sin(phi)
      cf = cos(phi)
      ss = sin(psi)
      cs = cos(psi)
      sk = sin(-kappa)
      ck = cos(-kappa)
** now calculate the product matrix of the five rotation matrices
      aa(1,1)=(cf*cf+sf*sf*cs*cs)*ck+sf*sf*ss*ss
      aa(1,2)=ss*ss*sf*cf*(ck-1.0)-cs*sk
      aa(1,3)=sf*cs*ss*(ck-1.0)+cf*sk*ss
      aa(2,1)=ss*ss*cf*sf*(ck-1.0)+cs*sk
      aa(2,2)=(sf*sf+cf*cf*cs*cs)*ck+cf*cf*ss*ss
      aa(2,3)=cf*ss*cs*(1.0-ck)+sf*ss*sk
      aa(3,1)=cs*ss*sf*(ck-1.0)-cf*ss*sk
      aa(3,2)=cf*ss*cs*(1.0-ck)-sf*ss*sk
      aa(3,3)=ss*ss*ck+cs*cs
** done.
      return
 100  do I=1,3
         do J=1,3
            aa(I,J)=0.0
         enddo
      enddo 
      do I=1,3 
         aa(I,I)=1.0
      enddo 
      return
      end
      
c----------------------------------------------------------------------------
c----------------------------------------------------------------------------
       
      subroutine rotate(mat,v1,v2)
      real mat(3,3),v1(3),v2(3)
c input mat,v1
c output: v2 = MAT*V1
c profiler says do loop is expensive! use explicit notation
c$$$        do I=1,3
c$$$        v2(I)=mat(1,I)*v1(1)+mat(2,I)*v1(2)+mat(3,I)*v1(3)
c$$$        enddo
      v2(1)=mat(1,1)*v1(1)+mat(2,1)*v1(2)+mat(3,1)*v1(3)
      v2(2)=mat(1,2)*v1(1)+mat(2,2)*v1(2)+mat(3,2)*v1(3)
      v2(3)=mat(1,3)*v1(1)+mat(2,3)*v1(2)+mat(3,3)*v1(3)
      return
      end
c----------------------------------------------------------------------------
c----------------------------------------------------------------------------

      subroutine refold_align_transform(p1,p2,p3,q1,q2,q3,Mat)
      implicit none

c inputs
c P1, p2 ,p3 are three points in space
c q1,q2,q3 are the same three points in space  but rotated and offset.
c outputs
c Mat is the matrix  needed to transform
c P->Q namely M*P+OFF = Q  
c the offset is an arbitrary translation and is not returned.
c it can be found from q1-Mat*P1 or q2-Mat*p2.
c if for some reason P and Q aren't the same set of points then
c the orientations are defined as follows
c the plane p1,p2,p3 and q1,q2,q3 are aligned
c and p1->p2 direction is aligned with q1->q2
 
      real p1(3),p2(3),p3(3),q1(3),q2(3),q3(3)
      real Mat(3,3) 
      real MatQ(3,3), MatP(3,3) 
      
 
      call refold_coord_sys(p1,p2,p3,MatP)
      call refold_coord_sys(q1,q2,q3,MatQ)
      call mat_multiply_transpose(MatQ,MatP,Mat)
      return
      end
c-----------------------------------------------------------------------
c----------------------------------------------------------------------------



      subroutine refold_coord_sys(p1,p2,p3,Mat)

      implicit none
c creates a unique orientation  matrix given three points
c this matrix uniquely specifies the orientation
c of the 3-d coordinate system given any three reference points.
c (the matrix is of course arbitrary wrt to any fixed rotation)
c however, the relative orientation of two such matrices can be compared
c to determine the rotation needed to transform one 3-d coordinate system into 
c the other.
c infact the transformation matrix taking coordinate system B into A
c  A*transpose(B)
c where A and B are the matrices returned from this routine.
      real p1(3),p2(2),p3(2),temp(3)
      real mat(3,3)
      call subvec(p1,p2,temp)
      call unitvec(temp,mat(1,1))
      call subvec(p3,p2,temp)
      call cros(mat(1,1),temp,mat(1,3))
      call unitvec(mat(1,3),mat(1,3))
      call cros(mat(1,3),mat(1,1),mat(1,2))
      return
      end
c----------------------------------------------------------------------------
c----------------------------------------------------------------------------
      
      subroutine mat_multiply_transpose(a,b,c_out) 
c returns result in c_out cems
      real a(3,3),b(3,3),c_out(3,3)
      integer i,j
      do i = 1,3
         do j = 1,3
            
            c_out(i,j) = b(i,1)*a(j,1)+b(i,2)*a(j,2)+b(i,3)*a(j,3)
         enddo
      enddo
      return
      end
c----------------------------------------------------------------------------
      real*4 function dotprod(v1,v2)
      real*4 v1(3),v2(3)
      dotprod = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)  
      return
      end
  
      subroutine addvec(v1,v2,v3)
      real*4 v1(3),v2(3),v3(3)
**    v3 = v1 - v2
      v3(1) = v1(1) + v2(1)
      v3(2) = v1(2) + v2(2)
      v3(3) = v1(3) + v2(3)
      return
      end

      subroutine subvec(v1,v2,v3)
      real*4 v1(3),v2(3),v3(3)
**    v3 = v1 - v2
      v3(1) = v1(1) - v2(1)
      v3(2) = v1(2) - v2(2)
      v3(3) = v1(3) - v2(3)
      return
      end

c---------------------------------------------------------------
    
      subroutine compute_loop_ss_gap(loop_list)

      implicit none
 
      include 'param.h'
      include 'structure.h'

c loads common blocks with values specified for template used in rms routine.
c in typical use a single template structure is compared against MANY possible
c structures taken from the VALL.  hence we treat this structure slightly
c assymetrically:  setting it here rather than placing it in the call 
c arguments of score_rms()
c allows us to avoid passing the template down through some hierarchy 
c of scoring functions.

car from read_ksync_file, loop_begin is the first UNaligned residue, and loop_end
car is the last UNaligned residue
c loop_begin is the residue number of the template structure
c loop_end is where the loop joins the template structure
c loop_size is NOT loop_end-loop_begin+1  
c loop_size is the length of the loop in the query (not the template) that is replacing
c the loop_begin to loop_end region of the template
c debugging
      real gap 
c input
      integer loop_list(5)


c package globals

car parameters for the current loop
      integer max_template_res,max_template_atoms
      parameter (max_template_res=4,max_template_atoms=3)
      integer template_map(max_template_res)
      real*8 p1(3,MAX_template_res*max_template_atoms),uu(3,3)
      real*8 p2(3,max_template_atoms*max_template_res)
      real*8 ww(max_template_atoms*max_template_res)
      real rel_templateH(max_template_res)
      real rel_templateE(max_template_res)
      real rel_templateL(max_template_res)
      common/rms_template/template_map,p2,ww,p1,uu,
     #     rel_templateH,rel_templateE,rel_templateL
      
car parameters for the entire template      
      integer max_temp_res
      parameter (max_temp_res=1000)
      real temp_cnitro(3,max_temp_res),temp_calpha(3,max_temp_res)
      real temp_ccarbon(3,max_temp_res),temp_coxygen(3,max_temp_res)
      real temp_cbeta(3,max_temp_res)
      character*1 template_ss(max_temp_res)
      integer template_residue
      common/template_atoms/temp_cnitro,temp_calpha,temp_ccarbon,
     #     temp_coxygen,temp_cbeta,template_ss,template_residue
     
      
c local
      integer i,j,k,count
      integer loop_begin,loop_end
      integer loop_size

      real template(3,max_template_res*max_template_atoms)
      character*1 type1,type2
c copy into named variables
      loop_size = loop_list(1)
      loop_begin = loop_list(2)   !last aligned
      loop_end = loop_list(3)     !first aligned in next seg
       
 
car compute gap      
      gap = 0.0
      do k = 1,3
         gap = gap+(temp_cnitro(k,loop_begin+1)
     #        -temp_cnitro(k,loop_end-1))**2
      enddo
      gap = sqrt(gap)

c copy template coordinates at loop ends into template variable so consecutive
c      write(0,*) 'compute_loop_ss_gap query loop',loop_begin,loop_end
      do k = 1,3
         count = 0
         do j = 0,1             
            count = count+1
            template(k,count) = temp_cnitro(k,loop_begin+j)
         enddo
         do j = 1,0,-1             
            count = count+1
            template(k,count) = temp_cnitro(k,loop_end-j)
         enddo
         
         do j = 0,1             
            count = count+1
            template(k,count) = temp_calpha(k,loop_begin+j)
         enddo
         do j = 1,0,-1          
            count = count+1
            template(k,count) = temp_calpha(k,loop_end-j)
         enddo
         
         do j = 0,1             
            count = count+1
            template(k,count) = temp_ccarbon(k,loop_begin+j)
         enddo
         do j = 1,0,-1             
            count = count+1
            template(k,count) = temp_ccarbon(k,loop_end-j)
         enddo
      enddo

car template_map handles the fact that the residues in template are 
car non contiguous. Residue 1 and 2 are the last aligned and first
car unaligned residue at the N-terminal end of the loop
car residues 3 and 4 are the last unaligned and first aligned residue
car at the C-terminal end of the loop
	
      template_map(1) = 1
      template_map(2) = 2
      template_map(3) = loop_size-1
      template_map(4) = loop_size
c      write(0,*) 'template_map',(template_map(k),k=1,4)

c set default weights
	
      call set_rms_weight(1,1.0) 
c copy in to common block and change precision
	
      do i = 1,max_template_atoms*max_template_res 
         do k = 1,3
            p2(k,i) = template(k,i)
         enddo
      enddo
c$$$      write(0,*) ' set weights '
c$$$      write(0,'(4(f9.3,1x))') ((ww(j*4+k),k=1,4),j=0,max_template_atoms-1)
      
c backtrack on template to determine the ss patterns
c idea is as follows, starting from end of template position, we backtrack working back into
c the template from the loop, hunting for the first encountered regular secondary structure.
c then we create a template for this region which is either of type loop, or the observed regular  secondary structure
c found by backtracking.  the concept is that in an insertion type situation either the old secondary structure in the
c template gets extended further, or it is really loop.
c
c for our purposes here we are going to step back one more than we already were stepping back (total of three)

	 
        
      type1 = 'L'
      if (template_ss(loop_begin).ne.'L') 
     #     type1 = template_ss(loop_begin)
      if (template_ss(loop_begin+1).ne.'L') 
     #     type1 = template_ss(loop_begin+1)         
      
      type2='L'
      if (template_ss(loop_end).ne.'L')  
     #     type2 = template_ss(loop_end)
      if (template_ss(loop_end-1).ne.'L') 
     #     type2 = template_ss(loop_end-1)
      
      do i = 1,max_template_res
         rel_templateH(i)=0.0
         rel_templateE(i)=0.0
         rel_templateL(i)=3.0
      enddo
	     
      if (type1.eq.'H') then 
         rel_templateH(1) = 7
         rel_templateH(2) = 7
      elseif (type1.eq.'E') then 
         rel_templateE(1) = 6
         rel_templateE(2) = 6
      endif
      
      if (type2.eq.'H') then
         rel_templateH(3) = 7
         rel_templateH(4) = 7
      elseif (type2.eq.'E') then 
         rel_templateE(3) = 6
         rel_templateE(4) = 6
      endif
      

      do i = 1,max_template_res/2
         type1 = template_ss(loop_begin+i-1)
         if (type1.eq.'H') rel_templateH(i) = rel_templateH(i)+7
         if (type1.eq.'E') rel_templateE(i) = rel_templateE(i)+6.5
         if (type1.eq.'L') rel_templateL(i) = rel_templateL(i)+6
         
         type1 = template_ss(loop_end-2+i)
         if (type1.eq.'H') rel_templateH(2+i) = rel_templateH(i+2)+7
         if (type1.eq.'E') rel_templateE(2+i) = rel_templateE(i+2)+6.5
         if (type1.eq.'L') rel_templateL(2+i) = rel_templateL(i+2)+6
      enddo
      
	
	
c$$$      write(0,*) "template_ss for loop ends  ",
c$$$     #     (template_ss(loop_begin+i),i=0,1),':',
c$$$     #     (template_ss(loop_end+i-2),i=1,2)
c$$$      write(0,*)"rel_template ss for loop ends"
c$$$      do i = 1,max_template_res
c$$$         write(0,*) rel_templateH(i),rel_templateE(i),rel_templateL(i)
c$$$      enddo
      
      return
      end

c-----------------------------------------------------------	
      subroutine score_template_ss(vall_ss,ss_score)
      implicit none
c this is a hardcoded scoring function that compares residues before and 
c after the loop
c with the secondary structure of the template.  a score of some sort is 
c returned
c note this essentially is "peaking" past the end of the loop region.  
c some loops may
c happen to lie on ends of proteins or be adjacent to "bad" spots.  for 
c these we have no 
c basis for comaprison and we just have to fake it.
      
c inputs
      character*1 vall_ss(*)
      real ss_score
        
c package globals
car parameters for the current loop
      integer max_template_res,max_template_atoms
      parameter (max_template_res=4,max_template_atoms=3)
      integer template_map(max_template_res)
      real*8 p1(3,MAX_template_res*max_template_atoms),uu(3,3)
      real*8 p2(3,max_template_atoms*max_template_res)
      real*8 ww(max_template_atoms*max_template_res)
      real rel_templateH(max_template_res)
      real rel_templateE(max_template_res)
      real rel_templateL(max_template_res)
      common/rms_template/template_map,p2,ww,p1,uu,
     #     rel_templateH,rel_templateE,rel_templateL

c locals
      integer k
      character*1 temp
      ss_score = 0.0
      do k = 1,max_template_res

C we sum the ss prediction reliabilities if they 
C for the actual ss at these positions in the vall
         temp = vall_ss(template_map(k))
         if (temp.eq.'H') then
            ss_score = ss_score-rel_templateH(k) 
         elseif(temp.eq.'E') then
            ss_score = ss_score-rel_templateE(k)
         elseif (temp.eq.'L') then 
            ss_score = ss_score-rel_templateL(k) 
         endif

      enddo 	
         

      return
      end
c---------------------------------------------------------------------	 

      subroutine set_rms_weight(scheme,level) 
      implicit none

c loads common blocks with values specified for template used in rms routine.
c in typical use a single template structure is compared against MANY possible
c structures taken from the VALL.  hence we treat this structure slightly
c assymetrically:  setting it here rather than placing it in the call arguments of score_rms()
c allows us to avoid passing the template down through some hierarchy of scoring functions.

      integer scheme
      real level
c package globals
car parameters for the current loop
      integer max_template_res,max_template_atoms
      parameter (max_template_res=4,max_template_atoms=3)
      integer template_map(max_template_res)
      real*8 p1(3,MAX_template_res*max_template_atoms),uu(3,3)
      real*8 p2(3,max_template_atoms*max_template_res)
      real*8 ww(max_template_atoms*max_template_res)
      real rel_templateH(max_template_res)
      real rel_templateE(max_template_res)
      real rel_templateL(max_template_res)
      common/rms_template/template_map,p2,ww,p1,uu,
     #     rel_templateH,rel_templateE,rel_templateL

c local
      integer k
      real*8 norm
      norm = 0.0
      do k = 0,max_template_atoms-1
         ww(1+k*max_template_res) = 0.0 !1e-2*level
         ww(2+k*max_template_res) =  1.0 ! *level*scheme
         ww(3+k*max_template_res) =  1.0
         ww(4+k*max_template_res) = 0.0 !1e-2
      enddo
c normalize to mean weight of 1 per atom
      do k = 1,max_template_res*max_template_atoms
         norm = norm +ww(k)
      enddo
	
      norm = max_template_res*max_template_atoms/norm
      do k = 1,max_template_res*max_template_atoms
         ww(k)=ww(k)*norm
      enddo
      
c now WW is normalized so that average to a mass of 1 per atom
  
      return
      end

c----------------------------------------------------	
	      
      subroutine score_rms(vall_pos,rms_score)

      implicit none
      

c inputs
      real vall_pos(3,5,*)
c output 
      real rms_score

c package globals
car parameters for the current loop
      integer max_template_res,max_template_atoms
      parameter (max_template_res=4,max_template_atoms=3)
      integer template_map(max_template_res)
      real*8 p1(3,MAX_template_res*max_template_atoms),uu(3,3)
      real*8 p2(3,max_template_atoms*max_template_res)
      real*8 ww(max_template_atoms*max_template_res)
      real rel_templateH(max_template_res)
      real rel_templateE(max_template_res)
      real rel_templateL(max_template_res)
      common/rms_template/template_map,p2,ww,p1,uu,
     #     rel_templateH,rel_templateE,rel_templateL
	
c local
      integer npoints
c debug 
      
      real d1,d2

c compares two structures.  this is a wrapper rouitne for calc_rms
c this fills arrays with atoms to be compared then computes rms_overlap
c this value is returned as the score.
c 
      real*8 ctx
      integer k,j,jj 
       
      logical verbose
      data verbose/.false./
  
c filling p1 array
 
      Npoints = max_template_atoms*max_template_res 
      do jj = 1,max_template_res
         j = template_map(jj)
         
	do k = 1,3
           p1(k,jj) = vall_pos(k,1,j) 
           p1(k,jj+max_template_res) = vall_pos(k,2,j)
           p1(k,jj+2*max_template_res) = vall_pos(k,4,j)
 	enddo
      enddo
        
      if (verbose) then 
         do jj=1,npoints
            write(0,*)'points',p1(1,jj),p2(1,jj),ww(jj)
         enddo
         write(0,*) '==========',max_template_res
         
         do jj = 1,npoints,max_template_res
            d1 = 0.0
            d2 = 0.0
            do k = 1,3
               d1 = d1+(p1(k,jj)-p1(k,5))**2
               d2 = d2+(p2(k,jj)-p2(k,5))**2
            enddo
            write(0,*) jj,sqrt(d1),sqrt(d2),sqrt(d1)-sqrt(d2)
         enddo
         do jj = 2,npoints,max_template_res
            d1 = 0.0
            d2 = 0.0
            do k = 1,3
               d1 = d1+(p1(k,jj)-p1(k,5))**2
               d2 = d2+(p2(k,jj)-p2(k,5))**2
            enddo
            
            write(0,*) jj,sqrt(d1),sqrt(d2),sqrt(d1)-sqrt(d2)
         enddo
      endif
c ************************************************
c dont blink!

	call findUU(p1,p2,ww,nPoints,uu,ctx)

c now have 3x3 rotation array UU
c will compute optimal RMS without
c having to actually rotate the arrays (its magic! but it
c does need the value of ctx from findUU )

	call calc_rms(rms_score,p1,p2,ww,npoints,ctx)

c rms_score is now the rms error
c ************************************************

	 
	return
	end
c--------------------------------------------------------------------------
	
c =======================================
c =====   the main event
c========================================		
      subroutine findUU(XX,YY,WW,Npoints,UU,sigma3)
      implicit none
      integer Npoints,sort(3)
      real*8 XX(3,Npoints),YY(3,Npoints),WW(Npoints),UU(3,3)
c purpose: intended to rotate one protein xyz array onto 
c another one such that the point-by-point rms is minimized.
c inputs:  XX,YY  are 2D arrays of x,y,z position of each atom
c 	   WW is a weight matrix for the points 
c	   Npoints is the number of XYZ points (need not be  physical array size)
c output:  UU is 3x3 rotation matrix. 
C          SIGMA3 IS TO BE PASSED TO OPTIONAL FAST_RMS CALC ROUTINE. 
c SIDEEFECTS: the Matrices XX, YY are Modified so that their weighted center of mass 
c				is  moved to (0,0,0).
c CAVEATS: 
c	   
c          1) it is CRITICAL that the first physical dimension of XX and YY is 3
c	   2) an iterative approx algorithm computes the diagonalization of
c              a 3x3 matrix.  if this program needs a speed up this could be
c               made into an analytic but uggggly diagonalization routine.
c	  
c____________________________
c Written by Charlie Strauss 1999
c	
c Mathethematical Basis from paper:
c (Wolfgang Kabsch) acta Cryst (1976) A32 page 922 
C REVISED APRIL 21
C   1) ORIGINAL PAPER HAD ERROR IN HANDEDNESS OF VECTORS, LEADING 
c      TO INVERSION MATRICIES ON OCCASION. OOPS. NOW FIXED.
C       SEE ACTA CRYST(1978) A34 PAGE 827 FOR REVISED MATH
C   2) TRAP DIVIDE BY ZERO ERRORS WHEN NO ROTATIONS REQUIRED.
C   3) ADDED WEIGHTS (WEIGHTS NOW WORK)
C   4) ADDED FAST RMS CALC AUXILIRARY ROUTINE.
c   5) CHANGED TO REAL*8 TO DEAL WITH HIGHLY DISSIMILAR BUT LARGE PROTEINS.
C      
c  Revised april 22
c   switched order of array subscripts so that can use logical array sizes
c_____________________________
c _________
C XX and YY are lists of Npoints  XYZ vectors (3xNpoint matrix) to be co-alligned
c these matrices are returned slightly modified: they are translated so their origins
c are at the center of mass (see Weights below ).
c WW is the weight or importance of each point (weighted RMS) vector of size Npoints
c The center of mass is figured including this variable.  
c UU is a 3x3 symmetric orthonornal rotation matrix that will rotate YY onto XX such that
c the weighted RMS distance is minimized.
c __________
      real*8 sigma3
      real*8 eVec(3,3),bb(3,3),w_w(3)
      real*8 m_moment(3,3),rr_moment(3,3)
      integer i, k,j,NROT
      real*8 temp1
      real*8 temp2
      real*8 temp3
      real*8 Ra(3)
c debugging
      integer err_unit
      parameter (err_unit=384)

c align center of mass to origin
      do 12 k = 1,3
         temp1=0.0
         temp2=0.0
         temp3=0.0
         do 11 j=1,Npoints
            temp1 = temp1+XX(k,j)*WW(j)
            temp2 = temp2+YY(k,j)*WW(j)
            temp3 = temp3+WW(j)
 11      continue
         temp1 = temp1/temp3
         temp2 = temp2/temp3
         do 10 j=1,Npoints
            XX(k,j) = XX(k,j)-temp1
            YY(k,j) = YY(k,j)-temp2
 10      continue
 12   continue

c     Make cross moments matrix   INCLUDE THE WEIGHTS HERE

      do 33 k=1,3
         do 32 j = 1,3
            
            TEMP1=0.0
            do 31 i = 1,npoints
               TEMP1 = TEMP1+ ww(i)*YY(k,i)*XX(j,i)
               
 31         continue
            m_moment(k,j) =TEMP1		
 32      continue
 33   continue
      
c	Multiply CROSS MOMENTS by transpose
      call BlankMatrixMult(m_moment,3,3,1,m_moment,3,0,RR_moment)

c find eigenvalues, eigenvectors of symmetric matrix RR_moment,

      call jacobi_xyz(RR_moment,w_w,eVec,NROT)

c	write(err_unit,*) "avec"
c	call debug(evec,3)
	

c explicitly coded 3 level index sort using eigenvalues
      do i = 1,3
         sort(i) = i
      enddo
      
      if (w_w(1).lt.w_w(2)) then
         sort(2)=1
         sort(1)=2
      endif
      
      if( w_w(sort(2)).lt.w_w(3) ) then
         sort(3)=sort(2)
         sort(2)=3
         
         if (w_w(sort(1)).lt.w_w(3)) then
            sort(2) = sort(1)
            sort(1) = 3
         endif
      endif
      
c sort is now an index to order of eigen values
 
c fix trivial roundoff errors.
      if (w_w(sort(3)).lt.0.0) then !(negative values theoretically  impossible)
 
         if (w_w(sort(3)).lt.-1e-3)  then ! Yikes! non-trivial error--should NEVER happen
            write(0,*) "absurd negative eigen value",w_w(sort(3))
            stop
         endif
         
         w_w(sort(3)) = 0.0     ! fix insignificant roundoff error 
      endif
c fix degenerate rotations
      if (w_w(sort(2)).le.0.0) then ! holy smokes,  two eigen values are zeros or heavens negative
         write(0,*) "**** degenerate rotation you fool!  ****"



c          return identity rotation matrix to moron
         do i=1,3
            do k=1,3
               uu(i,k) = 0.0
               if (i.eq.k) uu(i,i) = 1.0
            enddo
         enddo
         if (w_w(sort(1)).le.0.0) w_w(sort(1)) = 0.0 ! should NEVER happen
         sigma3=sqrt(w_w(sort(1)))
         return                 ! make like a prom dress and slip off
      endif

c sort eigen values
      temp1 = w_w(sort(1))
      temp2 = w_w(sort(2))
      w_w(3) = w_w(sort(3))
      w_w(2) = temp2
      w_w(1) = temp1	

	
c sort first two eigen vectors (dont care about third)
      do i = 1,3
         temp1 = evec(i,sort(1)) 
         temp2 = evec(i,sort(2))
         evec(i,1) =temp1
         evec(i,2) =temp2
      enddo


c april 20: the fix not only fixes bad eigen vectors but solves
c a problem of forcing a right-handed coordinate system
	 
      call fixEigenvector(eVec)  
c at this point we now have three good eigenvectors in a right hand coordinate system.

c make bb basis vectors	= moments*eVec

      call BlankMatrixMult(m_moment,3,3,0,eVec,3,0,BB)
c squirel away a free copy of the third eigenvector before normalization/fix
      do j = 1,3
         Ra(j) = bb(j,3)
      enddo
      
c normalize first two bb-basis vectors
c dont care about third since were going to replace it with b1xb2
c this also avoids problem of possible zero third eigen value
      do  j=1,2
         temp1 = 1/sqrt(w_w(j)) ! zero checked for above
         do k=1,3               !x,y,z
            bb(k,j) = bb(k,j)*temp1
         enddo 
      enddo 

	

c  fix things so that bb eigenvecs are right handed  	
	      
      call fixEigenvector(BB)	! need to fix this one too
c find  product of eVec and BB matrices

      call BlankMatrixMult(eVec,3,3,0,BB,3,1,UU)
C result is returned in UU.
		

c and lastly determine a value used in another routine to compute the 
c rms
      sigma3=0.0
      do j=1,3
         sigma3 = sigma3+BB(j,3)*RA(j)
c     write(0,*) "eigenvalues j=",j,w_w(j)
      enddo
      if (sigma3.lt.0.0) then 
         sigma3 = sqrt(w_w(1))+sqrt(w_w(2))-sqrt(w_w(3)) ! note sign on w_w(3) is negative
      else
         sigma3 = sqrt(w_w(1))+sqrt(w_w(2))+sqrt(w_w(3))
      endif
      end


c *********************************************************************	

      subroutine fixEigenvector(m_v)  
      real*8 m_v(3,3)
      
c 	m_v is a 3x3  matrix of 3 eigen vectors
c       replaces the third  eigenvector by taking cross product of
c 	of the first two eigenvectors

      real*8 norm
      
      m_v(1,3) = m_v(2,1)*m_v(3,2)-m_v(3,1)*m_v(2,2)
      m_v(2,3) = m_v(3,1)*m_v(1,2)-m_v(1,1)*m_v(3,2)
      m_v(3,3) = m_v(1,1)*m_v(2,2)-m_v(2,1)*m_v(1,2)
c     normalize it to 1 (should already be one but lets be safe)
         
      norm = sqrt(1/( m_v(1,3)**2+m_v(2,3)**2+m_v(3,3)**2))
      
      m_v(1,3) = norm*m_v(1,3)
      m_v(2,3) = norm*m_v(2,3)
      m_v(3,3) = norm*m_v(3,3)
      return
      end

c *********************************************************************	
	
      subroutine calc_rms(rms_out,XX,YY,WW,npoints,ctx) 
c companion routine for findUU ( it is optional )
c computes the minimum RMS deviation beteen XX and YY as though
c it rotated the arrays, without actually rotating them.
c NOTE: the XX, YY, WW must be the same as call to findUU
c (remember that findUU offests the XX and YY weighted  COM to the origin!)
c ctx is a magic number computed during find UU that is needed for this calculation
c rms_out is the  REAL*4  output value of the rms deviation.
      implicit none
      integer npoints
      real*4 rms_out
      real*8 xx(3,npoints),yy(3,npoints),ww(npoints)
      real*8 rms,ctx
      real*8 eVec(3,3),w_w(3)
      real*8 m_moment(3,3),rr_moment(3,3)
      common/findUU_public/w_w,evec,m_moment,rr_moment
      
      integer j,k
	 
      rms =  0.0
      
      do k = 1,npoints
         do j = 1,3
            rms = rms + ww(k)*(xx(j,k)**2 +yy(j,k)**2)
            
         enddo
      enddo
      
      rms = rms -2*ctx
c watch roundoff error
      if (rms.lt.0.0) then
         if (rms.lt.-1e-3) then
            write(0,*) "absurd negative rms in ctx calc",rms,ctx
            stop
         endif
c fix roundoff
         rms=0.0
      endif
      
      rms_out = sqrt(rms/real(nPoints)) ! retruurn a real*4
      return
      end

	   
	 
c *********************************************************************	


      subroutine BlankMatrixMult(A,n,np,transposeA,B,m,transposeB,AxB_out)
      integer n,np,transposeA,m,transposeB
      real*8 a(np,n),b(np,m),AxB_out(m,n)
      integer k,j
c     fills output matrix with zeros before calling matrix multiply
      do 12 k=1,m
         do 11 j=1,n
            AxB_out(k,j) = 0.0
 11      continue
 12   continue
      
      
      call matrixMult(A,n,np,transposeA,B,m,transposeB,AxB_out)
      return
      end
c *********************************************************************	

      subroutine MatrixMult(A,n,np,transposeA,B,m,transposeB,AxB_out)
c multiplys matrices A (npXn) and B (npXn). results in AxB.out
c IF THE MATRICES are SQUARE.  you can also multiply the transposes of these matrices
c to do so set the TransposeA or TransposeB flags to 1, otherwise they should be zero.
c NOTE WELL: transpose only works correctly for square matricies!
c ___________
c charlie strauss 1999
c__________________
      integer n,np,transposeA,m,transposeB
      real*8 a(np,n),b(np,m),AxB_out(m,n)
      integer k,j,i
      

      if  (transposeA.eq.0) then
         if (transposeB.eq.0) then
            
            do 23 k=1,m
               do22 j = 1,n
               do 21 i = 1,np
                  AxB_out(k,j)  = AxB_out(k,j)+ A(k,i)*B(i,j)
                  
 21            continue
               
 22         continue
 23      continue
         

         else
            do 13 k=1,m
               do 12 j = 1,n
                  do 11 i = 1,np
                     AxB_out(k,j)  = AxB_out(k,j)+ A(k,i)*B(j,i)
 11               continue
                  
 12            continue
 13         continue
            

         endif
      else   
         if  (TransposeB.eq.0) then
            do 33 k=1,m
               do 32 j = 1,n
                  do 31 i = 1,np
                     AxB_out(k,j)  = AxB_out(k,j)+ A(i,k)*B(i,j)
 31               continue
                  
 32            continue
 33         continue

         else
            do 43 k=1,m
               do 42 j = 1,n
                  do 41 i = 1,np
                     AxB_out(k,j) = AxB_out(k,j)+A(i,k)*B(j,i)
 41               continue
                  
 42            continue
 43         continue
            

         endif
      endif
      return
      
      end


c ***************************************************************
c following routine is from Numerical recipes

      SUBROUTINE jacobi_xyz(a,d,v,nrot)
c-------------------------------------------------------------------------------
      INTEGER n,np,nrot,NMAX
      REAL*8 a(3,3),d(3),v(3,3)
      PARAMETER (NMAX=500)
      INTEGER i,ip,iq,j
      REAL*8 c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
c-------------------------------------------------------------------------------
c solves for eigen-values of REAL*8-symetric matrix
c using jacobi method of iterative rotations
c fast for a three by three
c not neccesarily optimal- numerical recipies
c in fortran : Cambridge-University Press
c-------------------------------------------------------------------------------
      n=3   !!! not noramly , expect inner-product
      np=3  !!! from xyz*xyz of protein hydrophobic
      do 12 ip=1,n   !!! coordinates
        do 11 iq=1,n
          v(ip,iq)=0.
11      continue
        v(ip,ip)=1.
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.
13    continue
      nrot=0
      do 24 i=1,50
        sm=0.
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))
14        continue
15      continue
        if(sm.eq.0.)return
        if(i.lt.4)then
          tresh=0.2*sm/n**2
        else
          tresh=0.
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+
     # g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h
              else
                theta=0.5*h/a(ip,iq)
                t=1./(abs(theta)+sqrt(1.+theta**2))
                if(theta.lt.0.)t=-t
              endif
              c=1./sqrt(1+t**2)
              s=t*c
              tau=s/(1.+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.
23      continue
24    continue
      pause 'too many iterations in jacobi'
      return
      END

