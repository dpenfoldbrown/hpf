c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 8261M $
c  $Date: 2008-12-17 13:59:34 -0500 (Wed, 17 Dec 2008) $
c  $Author: (local) $

      character*5 version
c      parameter (version = 'v1_4 ') ! pb 08/02/05 prof+new quotas
      parameter (version = 'v1_3 ') !used for frag filenames5
      integer n_sizes               !the actual number of sizes used, 
      parameter (n_sizes = 2)       !different from array size max_len
      integer max_res,max_vall,max_nn,max_homologs,
     #        max_len ,max_ss_type 
      parameter (max_res=2000,max_nn=75,max_vall=2000000,
     #           max_homologs=20000,max_len=8,max_ss_type=9)  
                                         !max_len is length-1 
					 !where length is what you want 

      integer seq_type
      parameter (seq_type=1)
      integer phd_type
      parameter (phd_type=2)
      integer jones_type
      parameter (jones_type=3)
      integer rdb_type
      parameter (rdb_type=4)
      integer jufo_type
      parameter (jufo_type=5)
      integer helix_type
      parameter (helix_type=6)
      integer sheet_type
      parameter (sheet_type=7)
      integer loop_type
      parameter (loop_type=8)
CAR WARNING WARNING WARNING
car  nmr_type must be the last ss_type
      integer nmr_type
      parameter (nmr_type=max_ss_type)

car parameter for loop picking
      integer max_zones             !number of different zone files
      parameter (max_zones=5000)
      integer max_loops
      parameter (max_loops=25)
      integer max_depth
      parameter (max_depth=2000)
      integer max_loop_length       !dont pick database loops longer than this
      parameter (max_loop_length=17)  
 
