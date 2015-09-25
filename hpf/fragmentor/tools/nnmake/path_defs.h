c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 11090 $
c  $Date: 2006-11-02 15:17:28 -0500 (Thu, 02 Nov 2006) $
c  $Author: dekim $

C header file containing path names 
c WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
c WARNING!!!! do no alter the order of the variables!!! the order is crucial:
c see EQUIVELANCE statment in file_manager.f to understand why.  
c NOTE: both common
c blocks have to be in the same order for the same corresponding variables!!!
c WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
c lines of text designating file paths in the file path_defs.txt should 
c correspond to the order
c of the variables in the path_names common block.
      integer MAX_PATHS
      parameter (MAX_PATHS=17)

      character*132 loop_path,homolog_path,input_chkpt_path,
     #     read_vall2_path,
     #     input_pdb_path1,input_pdb_path2,make_dat_path,
     #     make_output_path,
     #     read_jufo_path,read_rdb_path,read_jones_path,read_phd_path,
     #     read_vall_path,setup_path,write_names_path,read_nmr_path,
     #     read_noe_path
      

      common/path_names/ loop_path,homolog_path,input_chkpt_path,
     #     read_vall2_path,
     #     input_pdb_path1,input_pdb_path2,make_dat_path,
     #     make_output_path,
     #     read_jufo_path,read_rdb_path,read_jones_path,read_phd_path,
     #     read_vall_path,setup_path,write_names_path,read_nmr_path,
     #     read_noe_path

c variables to hold the length of the strings.  these need to be in exactly
c the same order as the string variables above.  Do not alter their order in 
c the common block


      integer loop_path_l,homolog_l,input_chkpt_l,read_vall2_l,
     #     input_pdb_l1,input_pdb_l2,make_dat_l,make_output_l,
     #     read_jufo_l,read_rdb_l,read_jones_l,read_phd_l,
     #     read_vall_l,setup_l,write_names_l,read_nmr_l,
     #     read_noe_l
      

      common /path_lengths/loop_path_l,homolog_l,input_chkpt_l,
     #     read_vall2_l,
     #     input_pdb_l1,input_pdb_l2,make_dat_l,make_output_l,
     #     read_jufo_l,read_rdb_l,read_jones_l,read_phd_l,
     #     read_vall_l,setup_l,write_names_l,read_nmr_l,
     #     read_noe_l


      character*50 vall_name,vall_cst_coord_name,barcode_name
      common /file_names/vall_name,vall_cst_coord_name,barcode_name
