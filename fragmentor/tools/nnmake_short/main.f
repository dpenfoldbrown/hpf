c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 1.17 $
c  $Date: 2002/06/18 04:41:53 $
c  $Author: rohl $

      program main 
      implicit none 
c-----------------------------------------------------------
c  MAKES FRAGMENT FILES
c-----------------------------------------------------------
c  modified from cems nn_make_noe_chsft
c  added R.Bonneau's checkpoint file MSA method, homolog checking
c  added dipolar constraints,updated constraints
c  included frag merging of std 

      include 'structure.h'
      include 'path_defs.h'  
       
      integer i
      character*100 filename

      integer nzones
      character*100 zone_list(max_zones)

      real quota(max_ss_type), sum_quota
      real sum_seq
      real noise
      integer nmr_ratio  

      logical rdb_error,phd_error,jones_error,jufo_error
      logical chsft_exist,dipolar_exist,noe_exist

      data sizes /2, 8/      !must specify n_sizes values
      data nmr_ratio /3/     !nmr_ratio:1=mixing for nmr:ss/seq frags
c------------------------------------------------------------
      call file_manager()
      
      call setup(zone_list,nzones)
      
      call input_pdb_nn()       !read pdb of question structure
      call input_checkpoint()   !overwrites sequence and total_residue
      call read_phd(phd_error)     
      call read_jones(jones_error)
      call read_jufo(jufo_error) 
      call read_rdb(rdb_error)
      if (rdb_error) call read_dsc(rdb_error)

cmj set noise level between 0.00 and 1.00       
      noise = 0.0

car set quotas for different ss_types
      sum_quota = 0.0
      if(phd_error) then
         quota(phd_type) = 0.0
      else
         quota(phd_type) = .2
         sum_quota = sum_quota + quota(phd_type)
         write(0,*)"phd ss prediction in use"
      endif
      if(jones_error) then 
         quota(jones_type) = 0.0
      else
         quota(jones_type) = .6
         sum_quota = sum_quota + quota(jones_type)
         write(0,*)"jones ss prediction in use"
      endif
      if(rdb_error) then
         quota(rdb_type) = 0.0
      else
         quota(rdb_type) = .2
         sum_quota = sum_quota + quota(rdb_type)
         write(0,*)"rdb/dsc ss prediction in use"
      endif
      if(jufo_error) then 
         quota(jufo_type) = 0.0
      else
         quota(jufo_type) = .2
         sum_quota = sum_quota + quota(jufo_type)
         write(0,*)"jufo ss prediction in use"
      endif

      quota(seq_type) = 0.0 !0.105*sum_quota
      sum_quota = sum_quota + quota(seq_type)

      quota(helix_type)=0.0     !these ss_types for maintaining ss balance
      quota(sheet_type)=0.0
      quota(loop_type)=0.0
      
      if (sum_quota.eq.quota(seq_type)) then
         write(0,*)"NNMAKE FAILED: no valid  secondary structure "//
     #        "predictions found"
         stop
      endif


cmj compute average secondary structure predictions

      do i=1,total_residue
         rel_aveH(i) = 0.0
         rel_aveE(i) = 0.0
         rel_aveL(i) = 0.0
         if(quota(phd_type).ne.0.0) then
            rel_aveH(i) = rel_aveH(i) + rel_phdH(i)
            rel_aveE(i) = rel_aveE(i) + rel_phdE(i)
            rel_aveL(i) = rel_aveL(i) + rel_phdL(i)
            write(0,50)"phd: ",i,rel_phdH(i),rel_phdE(i),rel_phdL(i)
         endif
         if(quota(jones_type).ne.0.0) then
            rel_aveH(i) = rel_aveH(i) + rel_jonH(i)
            rel_aveE(i) = rel_aveE(i) + rel_jonE(i)
            rel_aveL(i) = rel_aveL(i) + rel_jonL(i)
            write(0,50)"jon: ",i,rel_jonH(i),rel_jonE(i),rel_jonL(i)
         endif
         if(quota(rdb_type).ne.0.0) then
            rel_aveH(i) = rel_aveH(i) + rel_rdbH(i)
            rel_aveE(i) = rel_aveE(i) + rel_rdbE(i)
            rel_aveL(i) = rel_aveL(i) + rel_rdbL(i)
            write(0,50)"rdb: ",i,rel_rdbH(i),rel_rdbE(i),rel_rdbL(i)
         endif
         if(quota(jufo_type).ne.0.0) then
            rel_aveH(i) = rel_aveH(i) + rel_jufoH(i)
            rel_aveE(i) = rel_aveE(i) + rel_jufoE(i)
            rel_aveL(i) = rel_aveL(i) + rel_jufoL(i)
            write(0,50)"jufo: ",i,rel_jufoH(i),rel_jufoE(i),rel_jufoL(i)
         endif
         write(0,*) !to print a newline
         sum_seq = rel_aveH(i) + rel_aveE(i) + rel_aveL(i)

cmj add 10% bias to the average predicted values 
         rel_aveH(i) = noise + (1. - 3. * noise) * rel_aveH(i)/sum_seq
         rel_aveE(i) = noise + (1. - 3. * noise) * rel_aveE(i)/sum_seq
         rel_aveL(i) = noise + (1. - 3. * noise) * rel_aveL(i)/sum_seq

cmj try to make fragments better by setting only one sstype if rel_ave > 0.8
c        if(rel_aveH(i).gt.0.75) then
c           rel_aveH(i)=1.0
c           rel_aveL(i)=0.0
c           rel_aveE(i)=0.0
c
c        elseif(rel_aveE(i).gt.0.75) then
c           rel_aveH(i)=0.0
c           rel_aveL(i)=0.0
c           rel_aveE(i)=1.0
c
c        elseif(rel_aveL(i).gt.0.75) then
c           rel_aveH(i)=0.0
c           rel_aveL(i)=1.0
c           rel_aveE(i)=0.0
c        endif
c        write(0,*)"ave: ",rel_aveH(i),rel_aveE(i),rel_aveL(i)
      enddo

c read experimental data
      call read_nmr_pred(chsft_exist)     !chsft file
      call read_noe_constraints(noe_exist)
      filename=read_noe_path(1:read_noe_l)//
     #     lower_name//chain_letter//'.dpl'
      call read_dipolar(filename,112,dipolar_exist)
      quota(nmr_type)=0.0       !experimental data- merged in separate routine
      if (.not.noe_exist.and.      !if no data, turn of nmr merging
     #     .not.dipolar_exist.and.
     #     .not.chsft_exist) nmr_ratio=0  

car normalize quotas
      do i=1,max_ss_type
         if(sum_quota.ne.0.0) then
            quota(i) = 1.2 * quota(i) / sum_quota
         endif
      enddo
      call  move_frag_set_quota(quota)
 
car read vall  
      call read_vall()            !read seq/struct files 
      if (nzones.gt.0 .or. noe_exist .or. dipolar_exist) 
     #     call read_vall_constraints_coord()
      
car read homologs
      call homologs()	    	!removes homologues in dunbrack's files 
      call homologs_nr()          !removes homologs -double check

      if (nzones.eq.0) then
         call make_frag_files(nmr_ratio)
car else pick loops
      else

         call make_loop_library(nzones,zone_list)
      endif

c close files opened in setup
      close(60)
      write(0,*)'NNMAKE DONE:: normal exit'
 50   format(2x,a,i3,3f5.2,$)
      end  
c-----------------------------------------------------------------------------

      subroutine make_frag_files(nmr_ratio)

      implicit none
      include 'structure.h'

car input
      integer nmr_ratio
cems local	
	integer q,isize

      integer i
      integer nnmake_tag
     
      integer ss_frags(max_len,max_res,max_nn)
      integer nmr_frags(max_len,max_res,max_nn)
      integer ssfrags_sstype(max_len,max_res,max_nn)
      integer nmrfrags_sstype(max_len,max_res,max_nn)

      if (use_pdb.eq.'N') then
         call make_dat()
         write(0,*)'not reading pdb file for ',lower_name  
      endif
      write(60,*)'total residue: ',total_residue 
      
      
      do q=1,n_sizes
         isize=sizes(q)+1  ! stupid size array is one smaller than actual size
         write(0,*)"size ",isize,"inserts",total_residue-isize
         write(60,*)'doing size ',isize
         call find_nn(isize)
      enddo

      write(0,*)"merging fragments"
car generate a single list of frags picked by secondary structure prediction
car using move_frag which does round robin, lose turn if over cap or frag is
car redundant

car if nmr data exist, the ss frags are put in a temporary list, and
car then merged with nmr-picked frags in a 3:1 nmr:ss ratio. Straight merging
car with redundant frags removed, (but counting for both lists with respect
car to the merge ratio)

car in either case, final frag list ends up in best_nn and best_nn_ss_type

      if (nmr_ratio.gt.0) then
         call move_frag(big_best_nn,ss_frags,ssfrags_sstype)
         call copy_bbnn_to_fraglist(big_best_nn,nmr_type,max_ss_type,max_nn,
     #        max_res,max_len,nmr_frags,nmrfrags_sstype)
         call merge_frag_lists(nmr_frags,nmrfrags_sstype,nmr_ratio,
     #        ss_frags,ssfrags_sstype,1,n_sizes,sizes,max_len,
     #        total_residue,max_res,max_nn,best_nn,best_nn_ss_type)
         nnmake_tag=6
      else
         nnmake_tag=5           !reflects quota mixture in move_frag
         call move_frag(big_best_nn,best_nn,best_nn_ss_type)
      endif

*******************************************************************
c output codes for frag names (nnmake_tag)
c    5 = ss pred only
c    6 = + experimental data
*******************************************************************

c diagnostics to status file:
      call make_ss_nn()           !make ss pred from 50 top best_nn only  
      call all_ss()               !prints all secondary structure info  
         
      do i=1,n_sizes
         call sizemer_dme(sizes(i))  
      enddo
      call write_names()   
      write(0,*)'done writing names...'
      
      call makeoutput(nnmake_tag) !write out final fragment files


      return
      end

c---------------------------------------------------------------------

      subroutine make_loop_library(nzones,zone_list)

      implicit none
      include 'structure.h'

car input
      integer nzones
      character*100 zone_list(nzones)
      
car local
      integer i
      integer N_loop_alignments
      integer loop_list(5,max_loops) 

      do i=1,nzones
         write(0,*)'NEXT ZONE FILE: ', 
     #        zone_list(i)(1:index(zone_list(i),' ')-1)
         call setup_template(zone_list(i),loop_list,n_loop_alignments)
         
         call find_loops(zone_list(i),loop_list,
     #        N_loop_alignments) 
      enddo
      
      return
      end







