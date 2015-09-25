c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 11090 $
c  $Date: 2006-11-02 15:17:28 -0500 (Thu, 02 Nov 2006) $
c  $Author: dekim $

      integer vall_size
      common/vall_read/vall_size         !set in read_vall.f
      integer sizes(n_sizes)             !set in setup.f
      common /frag_sizes/ sizes


car global variables     
      integer nseq,num_homologs   
      integer best_nn(max_len,max_res,max_depth),total_residue,acc(max_res) 
      integer best_nn_ss_type(max_len,max_res,max_nn)
      integer big_best_nn(max_ss_type,max_nn,max_res,max_len)
      integer*1 seqnum(max_res)      
      common/ints1/nseq,num_homologs,seqnum
      common/ints2/best_nn,total_residue,acc,big_best_nn, best_nn_ss_type

      real substitution(0:20,0:20),profile(0:21,max_res),
     #     acc_vall(max_vall),accaa(0:20,0:10)   
      common/everything/substitution,profile,acc_vall,accaa    

      character seq(max_res)*1,lower_name*4, 
     #          extension*10,use_pdb*1,residue3(max_res)*3,
     #          residue1(max_res)*1,ed_file*1,res_list(max_vall)*1,  
     #          chain_letter*1,resseq(max_vall)*5
      character homolog_id(max_homologs)*4,
     #          homolog_ch(max_homologs)*1,prefix_code*2  

      character force_bb_big_bin(max_res)*1
      character force_ss(max_res)*1

      common/chars/seq,lower_name,extension,use_pdb,
     #     residue3,residue1,ed_file,res_list,chain_letter,resseq,
     #     homolog_id,homolog_ch,prefix_code,force_bb_big_bin,force_ss

 
      character*1 ss(max_vall)
      integer back_seq(max_vall)  !only in read_vall
    
      real phi_list(max_vall),psi_list(max_vall), 
     #     omega_list(max_vall), 
     #     vall_pro(0:20,max_vall),avg_pro(0:20),
     #     vall_ca(3,max_vall), calpha(3,max_res),
     #     cnitro(3,max_res),ccarbon(3,max_res),
     #     coxygen(3,max_res),
     #      cnh(3,max_vall),  HCA(3,max_vall),native_cnh(3,max_res), 
     #      native_HCA(3,max_res)   
      character chain(max_vall)*1,frag_name(max_vall)*4 
    
      common/vall_char/chain,frag_name,ss  
      common/vall_total/phi_list,psi_list,
     #     omega_list,vall_pro
      common/others/avg_pro,calpha,vall_ca,cnitro,cnh,HCA,
     #              native_cnh,native_HCA,
     #              ccarbon,coxygen   
      common/vall_ints/ back_seq  
 
      real seq_weight,ss_weight,chsft_weight,weight_proline_cis,
     #     weight_proline_phi
      common/weight/seq_weight,ss_weight,chsft_weight,
     #     weight_proline_cis,weight_proline_phi

      character seq_phd(max_res)*1,ss_phd(max_res)*1 
      real rel_phd(max_res),rel_phdH(max_res),
     #        rel_phdE(max_res),rel_phdL(max_res) 
      common/phd_stuff/seq_phd,ss_phd,rel_phd,rel_phdH,rel_phdE,
     #                 rel_phdL
     
      character seq_jon(max_res)*1,ss_jon(max_res)*1 
      real rel_jon(max_res),rel_jonH(max_res),
     #     rel_jonE(max_res),rel_jonL(max_res) 
      common/jones_stuff/seq_jon,ss_jon,rel_jon,rel_jonH,rel_jonE,
     #                 rel_jonL       
 
      character seq_jufo(max_res)*1,ss_jufo(max_res)*1 
      real rel_jufo(max_res),rel_jufoH(max_res),
     #     rel_jufoE(max_res),rel_jufoL(max_res) 
      common/jufo_stuff/seq_jufo,ss_jufo,rel_jufo,rel_jufoH,
     #                 rel_jufoE,rel_jufoL       

      character seq_rdb(max_res)*1,ss_rdb(max_res)*1,ss_nn(max_res)*1 
      real rel_rdbH(max_res),rel_rdbE(max_res),
     #     rel_rdbL(max_res)
      common/rdb_stuff/seq_rdb,ss_rdb,ss_nn,
     #                 rel_rdbH,rel_rdbE,rel_rdbL
     
     
      real rel_aveH(max_res),rel_aveE(max_res),
     #     rel_aveL(max_res)
      common/ave_stuff/rel_aveH,rel_aveE,rel_aveL

      integer*1 vall_accept(20,max_vall)
      common/intts1/vall_accept

      real xyz_vall(3,5,max_vall),xyz_native(3,5,max_res),
     #       cbeta(3,max_res)

      common /vall_coord/ xyz_vall,xyz_native,cbeta
