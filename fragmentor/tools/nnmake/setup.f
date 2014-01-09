c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 11090 $
c  $Date: 2006-11-02 15:17:28 -0500 (Thu, 02 Nov 2006) $
c  $Author: dekim $

      subroutine setup(zone_list,nzones,ssonly) 
      implicit none 

      include 'path_defs.h'
      include 'param.h' 
      include 'structure.h'

car output
      integer nzones
      character*(*) zone_list(*)
      logical ssonly

car local
      integer iargc,i,iunit,j
      character*3 num2str
      character*30 listname

      data iunit /45/
c----------------------------------------------------------
      i=iargc()
      ssonly =.false.
      if (i.lt.3) then
         write(0,*)'NNMAKE requires three arguments'
         write(0,*)'  prefix code (char*2),'
         write(0,*)'  pdb name in lowercase (char*4)'
         write(0,*)"  chain name (char*1, '_' for none)"
         write(0,*)"e.g. 'pNNMAKE ta 1ubq _ ' "
         write(0,*)
         write(0,*)' 4 arguments required for forcing barcodes:'
         write(0,*)'pNNMAKE code pdbname chain -barcode_file file'
         write(0,*)'   OR tries pdbnamechain.bar'
         write(0,*)
         write(0,*)' 4 arguments required for making loop libraries:'
         write(0,*)'pNNMAKE code pdbname chain zone_file'
         write(0,*)'   OR'
         write(0,*)'pNNMAKE code pdbname chain -zonelist filename'
         write(0,*)
         write(0,*)'NOTE: zone files should be base name for zones, '
         write(0,*)'      template pdb, template ssa, and loop files'
         stop
      else
         call getarg(1,prefix_code)
         call getarg(2,lower_name)
         call getarg(3,chain_letter)
         nzones=0        
         barcode_name = lower_name//chain_letter//'.bar'
         if (i.gt.3) then
            call getarg(4,listname)
            if (listname.eq.'-ssonly') then
                ssonly=.true.
            elseif (listname.eq. '-zonelist') then
               nzones=1
               call getarg(5,listname)
               open(unit=iunit,file=listname,status='old')
               do while (.true.)
                  if (nzones.gt.max_zones) then
                     write(0,*)'NNMAKE FAILED: number of zone files '//
     #                   'exceeds limit'
                     write(0,*)' Increase max_zones in structure.h'
                     stop
                  endif
                  read(iunit,*,end=100)zone_list(nzones)
                  nzones=nzones+1
               enddo
 100           continue
               nzones=nzones-1
            elseif (listname.eq. '-barcode_file') then
               call getarg(5,barcode_name)
            else   !assume fourth arg is a zonefile
               nzones =1
               zone_list(nzones) = listname
            endif
         endif
       endif
       do i=1,nzones
         j = index(zone_list(i),'.zones')
         if (j.gt.0) zone_list(i)=zone_list(i)(1:j-1)
       enddo

c  USER DEFINES LOWER_NAME AND CHAIN_LETTER
      if (chain_letter.eq.'_') chain_letter=' '
      extension=  num2str(max_nn)//"_"//version


      write(0,*) 'setup 60 ',setup_path(1:setup_l)// 'status.'
     #     //extension(1:index(extension," ")-1)//'_'//prefix_code//
     #     lower_name//chain_letter
      open(60,file= setup_path(1:setup_l)//'status.'//
     #     extension(1:index(extension," ")-1)//'_'//prefix_code//
     #     lower_name//chain_letter,status='unknown')

      if (ssonly) then
         seq_weight=0.0
      else
         seq_weight=1.00 
      endif
      ss_weight=0.50
      chsft_weight=0.45
      weight_proline_phi = 1.0/log(10.0)
      write(60,*)'seq_weight: ',seq_weight 
      write(60,*)'ss_weight: ',ss_weight 
      write(60,*)'chsft_weight: ',chsft_weight
      write(0,*)'seq_weight: ',seq_weight 
      write(0,*)'ss_weight: ',ss_weight
      write(0,*)'chsft_weight: ',chsft_weight
      write(0,*) 'weight_proline_phi',weight_proline_phi

      write(60,*)'protein: ',lower_name,' chain ID: ',chain_letter


      return 
      end 
c---------------------------------------------------------------------
      function  num2str(no)
      implicit none
      integer no
      integer no1,no2,no3
      character*3  num2str
      character suf(0:9)*1
       data suf/'0','1','2','3','4','5','6','7','8','9'/

      no1=(mod(no,1000)-mod(no,100))/100      
      no2=(mod(no,100)-mod(no,10))/10      
      no3=mod(no,10) 
      write(0,*) no1,no2,no3
      num2str = suf(no1)//suf(no2)//suf(no3)
      return 
      end
  

