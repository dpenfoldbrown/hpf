c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 1.3 $
c  $Date: 2002/06/18 04:41:53 $
c  $Author: rohl $

	integer max_nmr
	parameter(max_nmr=400)

        real averphi(max_nmr),averpsi(max_nmr)
        real phivarprd(max_nmr),psivarprd(max_nmr)
        real log_mean_phi_err(max_nmr),log_mean_psi_err(max_nmr)    
       
        common /r_nmrpred/ phivarprd,psivarprd,averphi,averpsi, 
     #       log_mean_phi_err,log_mean_psi_err

	logical chsft_exist
	common  /l_nmrpred/ chsft_exist

