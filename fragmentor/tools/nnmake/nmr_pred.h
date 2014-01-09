c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 2338 $
c  $Date: 2002-08-09 18:23:11 -0400 (Fri, 09 Aug 2002) $
c  $Author: rohl $

        real averphi(max_res),averpsi(max_res)
        real phivarprd(max_res),psivarprd(max_res)
        real log_mean_phi_err(max_res),log_mean_psi_err(max_res)    
       
        common /r_nmrpred/ phivarprd,psivarprd,averphi,averpsi, 
     #       log_mean_phi_err,log_mean_psi_err

	logical chsft_exist
	common  /l_nmrpred/ chsft_exist

