c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 4746 $
c  $Date: 2004-06-15 16:04:40 -0400 (Tue, 15 Jun 2004) $
c  $Author: rohl $

      subroutine prepare_constraints_score(i,size,vall_size)

      implicit none
      
      integer i      !starting position in protein sequence
      integer size   !fragment size
      integer vall_size  !frags in vall
     
      call prepare_noe_score(i,size) 
      call prepare_dipolar_score(i,size,vall_size)

      return
      end

c-------------------------------------------------------------------------
      logical function get_chsft_exist()

      implicit none
      include 'param.h'
      include 'nmr_pred.h'
      
      get_chsft_exist=chsft_exist

      return
      end
c-------------------------------------------------------------------------
      logical function get_constraints_exist()

      implicit none
      include 'distConstraint.h'
      
      get_constraints_exist=constraints_exist

      return
      end
c-------------------------------------------------------------------------

      subroutine read_nmr_pred(exist) 

      implicit none           
      include 'param.h'
      include 'structure.h'
      include 'path_defs.h'
      include 'nmr_pred.h'

car output
      logical exist


      real nmr_psi(3,max_res),nmr_phi(3,max_res)
      real actphi(max_res),actpsi(max_res)
      real phisd(max_res),psisd(max_res)
      real nmr_ss(max_res),averss(max_res)
      integer which(max_res),aanum(max_res)
      real psi_nmr_score(0:361),phi_nmr_score(0:361)

      integer jj,kk
      real t
      character*180 file_to_open

      integer iunit

      iunit=98

      file_to_open = read_nmr_path(1:read_nmr_l)
     #     //lower_name//chain_letter//'.chsft'
      write(0,*) "Reading nmr chemical shift file:",
     #     file_to_open(1:index(file_to_open,' '))
      open(unit=iunit,file=file_to_open,status='old',iostat=jj)
      if (jj.ne.0) then
         write(0,*) "WARNING: can't find .chsft file (read_nmr_pred.f)"
c         write(0,*)'making frags using checkpoint and ss prediction only'
         chsft_exist=.false.
         exist=.false.
         do kk=1,total_residue   !set to zero so no attempt to use this data
            phivarprd(kk)=0.0
         enddo
         return
      endif

      write(0,*)"residue,which,averphi,averpsi,phivarprd,psivarprd,"
     #     //"nmr_phi _psi(1),nmr_phi _psi(2),nmr_phi _psi(3),actphi,"
     #     //"actpsi,phisd,psisd,nmr_ss,averss"
      
      do kk=1,total_residue
         
         read(iunit,10,err=200,end=100) aanum(kk),which(kk),averphi(kk),
     #        averpsi(kk),phivarprd(kk),psivarprd(kk),
     #        nmr_phi(1,kk),nmr_psi(1,kk),nmr_phi(2,kk),
     #        nmr_psi(2,kk),nmr_phi(3,kk),nmr_psi(3,kk),      
     #        actphi(kk),actpsi(kk),phisd(kk),psisd(kk),
     #        nmr_ss(kk),averss(kk)                                     
c$$$         write(0,10)   aanum(kk),which(kk),averphi(kk),
c$$$     #        averpsi(kk),phivarprd(kk),psivarprd(kk),
c$$$     #        nmr_phi(1,kk),nmr_psi(1,kk),nmr_phi(2,kk),
c$$$     #        nmr_psi(2,kk),nmr_phi(3,kk),nmr_psi(3,kk),      
c$$$     #        actphi(kk),actpsi(kk),phisd(kk),psisd(kk),
c$$$     #        nmr_ss(kk),averss(kk) 

         if ((averphi(kk).eq.180.00).and.(averpsi(kk).eq.-180.00))
     #        phivarprd(kk)=0.0 ! detect evil  !car should be 0 already
         
c     preprocess data for faster scoring: 
c     compute log of mean delta_phi error
c     compute inverse of mean_delta_phi error
c     if mean_delta_phi error is zero then this means there was no chem 
c     shift prediction
c   for this residue and we want to zero its contribution in scoring function
c note: its important to use natural logarithm       
         if (phivarprd(kk).gt.0.0) then
            log_mean_phi_err(kk)  = log( phivarprd(kk))
            phivarprd(kk)=1/phivarprd(kk)
         else
            log_mean_phi_err(kk)  = log(360.0)
            phivarprd(kk)= 0.0
         endif
         if (psivarprd(kk).gt.0.0) then
            log_mean_psi_err(kk)  = log( psivarprd(kk))
            psivarprd(kk)=1/psivarprd(kk)
         else
            log_mean_psi_err(kk)  = log(360.0)
            psivarprd(kk)= 0.0
         endif
         
      enddo

car -------------      
car as far as I can tell, these parameters are calculated and then
car never used
      t = 1.0
      t = -1.0
      do kk=1,361               ! degree bins
           
         t=kk*1.0
         if (kk.gt.180) t = abs(360.0-t)
         if (abs(t).gt.90.0) then
            psi_nmr_score(kk) = log(2.09)
            phi_nmr_score(kk) = log(2.50)
         else
            psi_nmr_score(kk) = log(2.09+83.46*exp(-(t/26.22 )**2))
            phi_nmr_score(kk) = log(2.50+55.423*exp(-(t/14.46)**2))
         endif
      enddo
car -------------      
      close(iunit)
      chsft_exist=.true.
      exist=.true.
      return

 10   format(i3,1x,i3,1x,16f9.2)

 100  write(0,*)'chain length mismatch in read_nmr_pred'
      write(0,*)'total_residue: ',total_residue
      write(0,*)'.chsft filelength: ', kk-1
      write(0,*)"NNMAKE FAILED: .cshft file does not contain the "
     #     //"correct number of residues; correct the .cshft_in file"
      stop

 200  write(0,*)"NNMAKE FAILED: format error in .cshft file line", kk,
     #     " Check the .cshft_in file format"
      stop
      end
c---------------------------------------------------

      function test_frag_noe(place)
      implicit none                          	
      include 'distConstraint.h'

      integer noe_paired_residue(1000),n_noe,residue_offset
      common/noe_frag/noe_paired_residue,n_noe,residue_offset
      integer i
      integer pair,place
      real dist
      
      real vall_noe_dist             ! function call,res_offset
      real test_frag_noe        ! return value
      test_frag_noe=0.0
      
      if (n_noe.eq.0) return  
      
      do i=1,n_noe
         pair = noe_paired_residue(i)
         dist = vall_noe_dist(pair,place)
         test_frag_noe=  test_frag_noe+ 
     #        ( max(dist- pairRadius(pair),0.0))**2 
      enddo            
      test_frag_noe = test_frag_noe
      return 
      
      end

c-------------------------------------------------------------------
      subroutine prepare_noe_score(start_res,dues)
car assemble a list of noes that can be evaluated given the current
car frag insertion point (start_res) and size (dues)

c both the starting and ending constraints must lie within the 
c bounds of the fragment.
c noe_paired_residue is a list of  the indicies of relevant constraints.
c n_noe is th number of valid constriaints in the list.
c residue_offset is the frag insertion point
      implicit none
      include 'distConstraint.h'
      integer noe_paired_residue(1000),n_noe,residue_offset
      common/noe_frag/noe_paired_residue,n_noe,residue_offset
      integer start_res,dues 
      integer i,lower_constraint,upper_constraint
      n_noe=0
      
      if (.not.constraints_exist) return

      do i =1,allvalid

car only evaluate HA/HN constraints         
c$$$         if (pair_atom_index(1,i).ne.-1 .and. 
c$$$     #       pair_atom_index(1,i).ne.0) goto 100
c$$$         if (pair_atom_index(2,i).ne.-1 .and. 
c$$$     #       pair_atom_index(2,i).ne.0) goto 100
         if (pair_atom_index(1,i).eq.999) goto 100     !check for undefined
         if (pair_atom_index(2,i).eq.999) goto 100 

car both residues must be within fragment
         lower_constraint = min(constraintpair(i,1),constraintpair(i,2))
         upper_constraint = max(constraintpair(i,1),constraintpair(i,2))
         if (lower_constraint.ge.start_res) goto 100
         if (upper_constraint.le.(dues+start_res-1)) goto 100

car add to list:
         n_noe = n_noe+1
         noe_paired_residue(n_noe)=i

 100     continue
      enddo
      residue_offset = start_res
      return
      end 
c------------------------------------------------------------------------      
      real function vall_noe_dist(pair,place)
                               
      implicit none
      include 'distConstraint.h' ! for the pair atom index  

      integer pair,place

      integer noe_paired_residue(1000),n_noe,residue_offset
      common/noe_frag/noe_paired_residue,n_noe,residue_offset
      integer i,j 
      real xyz1(3),xyz2(3)
      
      i = constraintPair(pair,1)-residue_offset+place 
      j = constraintPair(pair,2)-residue_offset+place  
         
      call get_vall_xyz(i,pair_atom_index(1,pair),xyz1)
      call get_vall_xyz(j,pair_atom_index(2,pair),xyz2)

      vall_noe_dist = ( xyz1(1)- xyz2(1) ) **2
     #     + (xyz1(2)-xyz2(2))**2
     #     + (xyz1(3)-xyz2(3))**2
                                
      vall_noe_dist=sqrt(vall_noe_dist)

      return
      end
c---------------------------------------------------------
      subroutine get_vall_xyz(res,type_index,xyz)

      implicit none
      include 'param.h'
      include 'structure.h'
      
      integer res
      integer type_index
      real xyz(3)

      if  (type_index.eq.0) then
         xyz(1)=cnh(1,res)
         xyz(2)=cnh(2,res)
         xyz(3)=cnh(3,res)
      elseif (type_index.eq.-1) then
         xyz(1)=hca(1,res)
         xyz(2)=hca(2,res)
         xyz(3)=hca(3,res)
      elseif (type_index.ge.1 .and. type_index.le.5) then
         xyz(1)=xyz_vall(1,type_index,res)
         xyz(2)=xyz_vall(2,type_index,res)
         xyz(3)=xyz_vall(3,type_index,res)
      else
         write(0,*)'unable to determine coordinates'
         write(0,*)'res ',res, ' type_index ',type_index
         stop
      endif
      return
         
 
      end
c---------------------------------------------------------
c------------------------------------------------------------------------
      subroutine read_noe_constraints(noe_exist,noe_count)

car read NMR_v3.0, calls read_constraints_old to handle older versions.
car mondified version for nnmake based on read_constraints

      implicit none                      
      
      include 'distConstraint.h'
      include 'path_defs.h'
      include 'param.h'
      include 'structure.h'      !for sequence
      

      logical noe_exist
      integer n_pairLines, noe_count

car comment
car comment
car comment
car n_pairLines                 (# of lines of constraints to read)
car tag,res1,atom1,res2,atom2,upperbound,lowerbound,{true distance}
car
car format(a1,2x,i4,1x,a4,1x,i4,1x,a4,1x,f10.2,1x,f10.2)
car tag: '#':ignore,(actually any non-whitespace char)
car true distance is optional and not read
car atom1 and atom2 should follow pdb-style atom names with IUPAC numbering
car also allowed: ' CEN' for centroid constraints (not evaluted in fullatom mode)
car degenerate hydrogen references allowed (no degenerate heavy atoms)

      real degenerate_padding
      parameter (degenerate_padding = 1.3)
      real general_padding
      parameter (general_padding = 1.3)

      character*200 filename
      integer j,k,kk,helpme 
      character*1 tag
      character*4 attype(MAX_CONSTRAINTS,2) 
      character*44 diagnostic_line(MAX_CONSTRAINTS)
      character*64 comment
      
      integer type2index

      integer cnstr_x  

      noe_count = 0
      cnstr_x=112
      filename=read_noe_path(1:read_noe_l)//lower_name//chain_letter
     #     //'.cst'
      open(cnstr_x, file=filename,status='OLD',
     #     iostat=helpme)
      if (helpme.ne.0) then 
         write(0,*)" XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
         write(0,*)"WARNING: CONSTRAINT FILE NOT FOUND"
         write(0,*) "Searched for: ",filename(1:index(filename,' '))
         WRITE(0,*) "  --will use no constraints --"
         write(0,*)" XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
         pc_weight=0.0
         pc_score=0.0
         np_pairs=0
         allValid=0
         n_pairLines=0
         constraints_exist = .false.
         noe_exist=.false.
         return 
      endif
	   
      constraints_exist = .true.
      noe_exist=.true.

         pc_weight=100.0/30.0   !std weight
         read(cnstr_x,'(A64)',err=1000) comment ! comment line
         write(0,'(A64)') comment
         read(cnstr_x,'(A64)',err=1000) comment ! comment line
         write(0,'(A64)') comment
         read(cnstr_x,'(A64)',err=1000) comment ! comment line
         write(0,'(A64)',err=1000) comment 
         read(cnstr_x,*,err=2000) n_pairLines

         k=1                    !current constraints
         j=1                    !current diagnostic constraint
         allValid=0             !total number of constraints
            
         do kk=1,n_pairLines
            read(cnstr_x,'(a1)',end=111,err=3000)tag
            if (tag.ne. ' ') goto 112 !skip line
            allValid=allValid+1
            backspace(cnstr_x)
            if (tag.ne. ' ') then
               read(cnstr_x,'(A44)',err=4000)diagnostic_line(j)
               j=j+1
            else
               read(cnstr_x,36,err=5000)tag,constraintPair(k,1),attype(k,1),
     #              constraintPair(k,2),attype(k,2),pairRadius(k),
     #              pairMinRadius(k)
car check for protons for which heavy atom coordinates will be used
car and pad bounds accordingly
car terminal H constraints,degenerate hydrogens
            if (constraintPair(k,1) .eq.1 .and. attype(k,1)(2:3).eq.'H ') then
               pairRadius(k)= pairRadius(k) + degenerate_padding
            elseif (attype(k,1)(1:1).eq.'#') then
                pairRadius(k)= pairRadius(k) + degenerate_padding
            endif
            if (constraintPair(k,2) .eq.1 .and. attype(k,2)(2:3).eq.'H ') then
              pairRadius(k)= pairRadius(k) + degenerate_padding
            elseif ( attype(k,2)(1:1).eq.'#') then
              pairRadius(k)= pairRadius(k) + degenerate_padding
            endif


car add general padding to all constraints
            pairRadius(k) = pairRadius(k) + general_padding

               k=k+1
            endif
 112        continue
         enddo
         j=j-1

 111     close(cnstr_x)
         np_pairs=k-1
         do kk=1,j
            read(diagnostic_line(kk),36)tag,constraintPair(k,1),
     #           attype(k,1),constraintPair(k,2),attype(k,2),
     #           pairRadius(k),pairMinRadius(k)
            k=k+1
            if (k.gt.MAX_CONSTRAINTS) then
               write(0,*)'WARNING: '//
     #              'increase MAX_CONSTRAINTS in distConstraint.h'
               write(0,*)' MAX_CONSTRAINTS ',
     #              MAX_CONSTRAINTS,' restraints read ',k
               write(0,*)"NNMAKE FAILED: .dpl file has more than ",
     #              MAX_CONSTRAINTS," restraints"
               stop
            endif
         enddo
         k=k-1
         noe_count = k
            
         if (k.ne.allValid .or. j.ne.allValid-np_pairs) then 
            write(0,*)'NNMAKE FAILED: error counting lines in .cst file'
            stop
         endif

         do k=1,allValid
            constraintPair(k,3)=seqnum(constraintPair(k,1))  !for nnmake
            constraintPair(k,4)=seqnum(constraintPair(k,2))
            pair_atom_index(1,k)=type2index(attype(k,1),constraintPair(k,3))
            pair_atom_index(2,k)=type2index(attype(k,2),constraintPair(k,4))
            pairStage(k)=abs(constraintPair(k,1)-constraintPair(k,2))
         enddo

      write(0,*)' stage res1  atom1 res2  atom2     upper     lower'
      do  k=1,allValid
         write(0,77)pairStage(k),constraintPair(k,1),pair_atom_index(1,k),
     #        constraintPair(k,2),pair_atom_index(2,k),pairRadius(k),
     #        pairMinRadius(k)
         if (k.eq.np_pairs) write(0,*)
     #        '------------diagnostic constraints-----------------'
      enddo
   
      write(0,*)'------------------------------------------------------'
      return
 36   format(a1,2x,i4,1x,a4,1x,i4,1x,a4,1x,f10.2,1x,f10.2)
 77   format(5(i5,1x),2(f10.2))
 1000 write(0,*)"NNMAKE FAILED: format error in .cst file comment lines 1-3"
 2000 write(0,*)"NNMAKE FAILED: format error in .cst file line 4"
 3000 write(0,*)"NNMAKE FAILED: format error in .cst file tag field"
 4000 write(0,*)"NNMAKE FAILED: format error in .cst file diagnostic line"
 5000 write(0,*)"NNMAKE FAILED: format error in .cst file data line"
      stop
      end

c----------------------------------------------------------------------
      integer function type2index(atom,restype)

      implicit none

      character*4 atom
      integer restype     !residue type

      logical IUPAC
      parameter (IUPAC = .false.)

      if (atom(1:1).ne.' ' .and. atom(1:1).ne.'1' .and. atom(1:1).ne.'2'
     #    .and. atom(1:1) .ne. '3' .and. atom(1:1) .ne. '#') then
         atom(2:4)=atom(1:3)      !check that fields are aligned
         atom(1:1)=' '
      endif

      type2index=999
car bb protons:
      if     (atom.eq.' H  ') then
         type2index = 0
      elseif (atom.eq.' HA ') then
         type2index = -1
      elseif (atom.eq.'1HA ' ) then
         if (.not. IUPAC) type2index = -1
      elseif (atom.eq.'2HA ' ) then
         type2index = -2
         if (IUPAC)  type2index =-1
      elseif (atom.eq.'3HA ' ) then
         if (IUPAC) type2index = -2


car bb heavy atoms  :
      elseif (atom.eq.' CA ') then
         type2index = 2
      elseif (atom.eq.' N  ') then
         type2index = 1
      elseif (atom.eq.' C  ') then
         type2index = 4
      elseif (atom.eq.' CB ') then
         type2index = 3
      elseif (atom.eq.' O  ') then
         type2index = 5

car differs for nnmake, sc atoms are undefined because we don't have
car these coordinates      

      endif

      return
         
 100  write(0,*)"Atom type: ",atom," unknown to type2index"
      stop
      
      end

c--------------------------------------------------------------------------
      subroutine delete_noe_constraints()

cmj function to delete all dipllar couplings
      implicit none
      include 'distConstraint.h'

      constraints_exist = .false.
      np_pairs = 0
      allValid = 0
      

      return
      end







