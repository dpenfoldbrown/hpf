c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 3382 $
c  $Date: 2003-09-10 19:18:51 -0400 (Wed, 10 Sep 2003) $
c  $Author: rohl $

      subroutine input_pdb_nn()

car differs from input_pdb, hence nn tag         
c  original author: Kim T Simons 

      implicit none 

      include 'param.h'
      include 'structure.h'
      include 'path_defs.h'
      real x,y,z
      character res_name*3,pdb_card*6,atom_name*4,
     #          chain_id*1,res_pdb*5,res_conf*1,prev_pdb*5    
      integer res_count,i,j 
      integer connect(max_res)
      character save_chain_letter*1
      
      integer iunit
c------------------------------------------------ 
      iunit=11
      save_chain_letter=chain_letter

      total_residue=0 
      res_count=0 
      do i=1,max_res 
         residue3(i)='   ' 
         residue1(i)=' ' 
         connect(i)=0 
         do j=1,3 
            calpha(j,i)=0.0  
            cbeta(j,i)=0.0  
            ccarbon(j,i)=0.0  
            coxygen(j,i)=0.0 
            cnitro(j,i)=0.0  
         enddo 
      enddo 


      use_pdb='Y' 
      write(0,*) 'input_pdb',iunit,input_pdb_path1(1:input_pdb_l1)//
     #     lower_name//'.pdb'
      open(iunit,file=input_pdb_path1(1:input_pdb_l1)//
     #     lower_name//'.pdb',status='old',iostat=i)
      if (i.ne.0) then 
         write(60,*)'using pdb file from database' 
         write(0,*) 'input_pdb',iunit,input_pdb_path2(1:input_pdb_l2)//
     #        lower_name//'.pdb'
         
         open(iunit,file=input_pdb_path2(1:input_pdb_l2)//
     #        lower_name//'.pdb',status='old',iostat=i)
         if (i.ne.0) then 
            write(60,*)'no pdb file for ',lower_name 
            use_pdb='N' 
            return  
         endif 
      endif 

      if (chain_letter.eq.'_') then 
         chain_letter=' ' 
      endif 
 
      atom_name='  '
      pdb_card='      '
      prev_pdb='     '
      chain_id='z' 
      write(0,*) chain_letter
 40   do while (atom_name.ne.' N  '.or.chain_id.ne.chain_letter.or.
     #     (pdb_card.ne.'HETATM'.and.pdb_card.ne.'ATOM  '))
         read(iunit,100,err=40,end=500)pdb_card,atom_name,res_conf,
     #        res_name,chain_id,res_pdb,x,y,z 
         
      enddo 
      do while (chain_id.eq.chain_letter)
         if (res_count.gt.max_res) then
            write(0,*)'need to increase max_res in structure.h'
            write(0,*)'max_res: ',max_res,' res_count: ',res_count
            stop
         endif
         if (pdb_card.eq.'TER   '.or.pdb_card.eq.'ENDMDL') goto 500  
         if (pdb_card.eq.'ATOM  '.or.pdb_card.eq.'HETATM') then   
            if (atom_name.eq.' N  ') then 
               if (res_pdb.ne.prev_pdb) then  
                  res_count=res_count+1
               endif 
               
               prev_pdb=res_pdb    
               cnitro(1,res_count)=x 
               cnitro(2,res_count)=y 
               cnitro(3,res_count)=z 
               residue3(res_count)=res_name 
            elseif (atom_name.eq.' CA ') then 
               calpha(1,res_count)=x 
               calpha(2,res_count)=y 
               calpha(3,res_count)=z  
               
            elseif (atom_name.eq.' C  ') then 
               ccarbon(1,res_count)=x 
               ccarbon(2,res_count)=y 
               ccarbon(3,res_count)=z 
            elseif (atom_name.eq.' CB ') then
               cbeta(1,res_count)=x
               cbeta(2,res_count)=y
               cbeta(3,res_count)=z
            elseif (atom_name.eq.' O  ') then
               coxygen(1,res_count)=x
               coxygen(2,res_count)=y
               coxygen(3,res_count)=z
            endif
         endif 
         read(iunit,100,end=500,err=500)pdb_card,atom_name,
     #        res_conf,res_name,chain_id,res_pdb,x,y,z  
      enddo 
 500  continue 
      total_residue=res_count 
      if (res_count.eq.0) then 
         write(0,*) "XXXXXXXXX***********XXXXXXXXXX********XXXXXXX"
         write(0,*) "WARNING: pdb file exists but cannot be read!!!"
         write(0,*) "XXXXXXXXX***********XXXXXXXXXX********XXXXXXX"
         
         use_pdb='N' 
         return
      endif
      
      write(0,*) " from pdb file, total_residue",total_residue
      if (total_residue.gt.max_res)
     #     write(17,*)'need to increase size of total_cbeta',
     #     ' in parameter file for structure ',lower_name  
      close(iunit) 

      call convert3_num  
      if (residue1(i).eq.'G') call compute_CB_coord(cnitro(1,i),
     #        calpha(1,i),ccarbon(1,i),cbeta(1,i))
	
      chain_letter = save_chain_letter
 
      call build_protons(cnitro,calpha,ccarbon,total_residue,
     #        native_hca,native_cnh)    !note these coord aren't used

      
      return 
 50   format(A4)
 100  FORMAT(A6,6X,A4,A1,A3,1x,A1,A5,3x,3f8.3)
      end 

c----------------------------------------------------------------------------
      subroutine build_protons(N,CA,C,nres,HA,HN)
      implicit none

      integer i,j

      real N(3,*)
      real CA(3,*)
      real C(3,*)
      integer nres
      real HA(3,*)
      real HN(3,*)
      
      do i=1,nres
         call compute_hxyz(-1,N(1,i),CA(1,i),CA(1,i),HA(1,i))
         if (i.gt.1) then
            call compute_hxyz(0,C(1,i-1),N(1,i),CA(1,i),HN(1,i))
         else
            do j=1,3
               HN(j,i)=N(j,i)
            enddo
         endif
      enddo

      return
      end
      

c-----------------------------------------------------------------------
            
      subroutine input_pdb_template(filename,iunit,
     #     cnitro,calpha,ccarbon,coxygen,cbeta,max_res,total_residue)        
 
      implicit none 

      include 'path_defs.h'

car input
      character*(*) filename
      integer iunit
      integer max_res
 
car output
      integer total_residue
      real cnitro(3,max_res),calpha(3,max_res),ccarbon(3,max_res)
      real coxygen(3,max_res),cbeta(3,max_res)
      

car local
      real x,y,z  
      character res_name*3,pdb_card*6,atom_name*4,
     #     chain_letter*1,chain_id*1,res_pdb*5,res_conf*1,prev_pdb*5    
      integer res_count,i,j
	

      total_residue=0 
      res_count=0 
      do i=1,max_res 
         do j=1,3 
            calpha(j,i)=0.0  
            ccarbon(j,i)=0.0  
            coxygen(j,i)=0.0  
            cnitro(j,i)=0.0  
            cbeta(j,i)=0.0
         enddo 
      enddo 

      write(0,*) 'input_pdb_template 11 ',filename(1:index(filename,' '))
      
      open(iunit,file=filename,status='OLD',iostat=i)
      

      if (i.ne.0) then 
         write(0,*)'no pdb file for ',filename(1:index(filename,' '))
         stop
      endif
      
      atom_name='  '
      pdb_card='      '
      prev_pdb='     '
      chain_id='z' 
 40   do while (atom_name.ne.' N  ' .or.
     #     (pdb_card.ne.'HETATM'.and.pdb_card.ne.'ATOM  '))
         read(iunit,100,err=40,end=500)pdb_card,atom_name,res_conf,
     #        res_name,chain_id,res_pdb,x,y,z 
      enddo 
      chain_letter=chain_id
      do while (chain_id.eq.chain_letter.and.res_count.lt.max_res)   
         if (pdb_card.eq.'TER   '.or.pdb_card.eq.'ENDMDL') goto 500  
         if (pdb_card.eq.'ATOM  '.or.pdb_card.eq.'HETATM') then   
            if (atom_name.eq.' N  ') then 
               if (res_pdb.ne.prev_pdb) then  
                  res_count=res_count+1 
               endif 
               prev_pdb=res_pdb  
               cnitro(1,res_count)=x 
               cnitro(2,res_count)=y 
               cnitro(3,res_count)=z 
            elseif (atom_name.eq.' CA ') then 
               calpha(1,res_count)=x 
               calpha(2,res_count)=y 
               calpha(3,res_count)=z  
               
            elseif (atom_name.eq.' C  ') then 
               ccarbon(1,res_count)=x 
               ccarbon(2,res_count)=y 
               ccarbon(3,res_count)=z 
            elseif (atom_name.eq.' CB ') then
               cbeta(1,res_count)=x
               cbeta(2,res_count)=y
               cbeta(3,res_count)=z
            elseif (atom_name.eq.' O  ') then
               coxygen(1,res_count)=x
               coxygen(2,res_count)=y
               coxygen(3,res_count)=z
            endif
         endif 
         read(iunit,100,end=500,err=500)pdb_card,atom_name,
     #        res_conf,res_name,chain_id,res_pdb,x,y,z  
      enddo 
 500  continue 
      total_residue=res_count 
      if (total_residue.gt.max_res) then
         write(0,*)'increase  max_res in structure.h'
         write(0,*)'max_res: ',max_res,' total_residue: ',total_residue
         stop
      endif
      close(iunit) 
       
      write(0,*) ' return from pdb_input with ',total_residue

      return 
 50   format(A4)
 100  FORMAT(A6,6X,A4,A1,A3,1x,A1,A5,3x,3f8.3)
      end 

      
c-----------------------------------------------------------------------------
      subroutine torsion_from_position(nres,xyz,phi,psi,omega)
*** determines phi,psi angles from a protein structure
*** non-existant angles set to 0.0
******************************************************************

      implicit none

      integer nres              !# of residues in the protein, logical
      real phi(nres),psi(nres),omega(nres)
      real xyz(3,5,nres)        !(xyz,atom#,res#) atom order:N CA CB C O 
      
car local variables
      integer i
      real dihedral


      phi(1)=0.0
      do i=1,nres
         if (i.gt.1) phi(i)=dihedral(xyz(1,4,i-1),xyz(1,1,i),xyz(1,2,i),
     +        xyz(1,4,i))
         if (i.lt.nres) then
            psi(i)=dihedral(xyz(1,1,i),xyz(1,2,i),xyz(1,4,i),
     +           xyz(1,1,i+1)) !psi
            omega(i)=dihedral(xyz(1,2,i),xyz(1,4,i),xyz(1,1,i+1),
     +           xyz(1,2,i+1))  !omega
         endif
      enddo
      psi(nres)=0.0
      omega(nres)=0.0

      return
      end

c---------------------------------------------------------------------
      real function dihedral(a1,a2,a3,a4)

car returns dihedral angle from four points.
car  Define a coordinate system in which vector v23 lies on the Z axis
car  vector v21 lies in the xy plane.

car functions from bystroff.f


      implicit none

      real a1(3),a2(3),a3(3),a4(3)
      real v23(3),v21(3),v34(3),vy(3),vx(3),uy(3),ux(3),cx,cy

      real dotprod,cc

      real conversion
      parameter (conversion=180.0/3.14159265359)

      call subvec(a3,a2,v23)
      call subvec(a1,a2,v21)
      call cros(v23,v21,vy)           !vy is vector along y axis
      call unitvec(vy,uy)              !uy is unit vector on y axis
      call cros(uy,v23,vx)            !vx is vector on x axis
      call unitvec(vx,ux)              !ux is unit vector on y axis
      call subvec(a4,a3,v34)
      cx = dotprod(v34,ux)             !cx: component of v34 on x axis
      cy = dotprod(v34,uy)             !cy: component of v34 on y axis
      cc = atan2(cy,cx)
      cc = cc*conversion     !convert to degrees
      dihedral = cc

      return
      end                                                                  

	

c---------------------------------------------------------------------

      subroutine compute_CB_coord(cnitro,calpha,ccarbon,cbeta)

car differs for rosetta conventions cross->cros

car returns hxyz, the coordinates of a proton of specified type.
car types:
car   0  HN         !all protons should be <= 0
car  -1  HA  or HA1
car  -2  HA2        !glycine
car The coordinates of three atoms are required to place the proton
car atom2=atom to which proton is bound
car atom1=N-terminally adjacent to atom2, (or closer to backbone)      
car atom3-C-terminally adjacent to atom2, (or further from backbone)

      implicit none

car input
      real cnitro(3),calpha(3),ccarbon(3)

car output
      real cbeta(3)
car internal
      integer i
      real axes(3,3)

car paramters
      real*8 CBpos(3)  !location of CB atom measured from CA with
                       ! x axis along N-CA bond,  yaxis in plane defined
                       ! by cA->C bond vector
      data CBpos /-0.532605943, 0.775186438, 1.19513831 /


      do i=1,3
         if (cnitro(i).ne.0.0) goto 100
         if (calpha(i).ne.0.0) goto 100
         if (ccarbon(i).ne.0.0) goto 100
      enddo
      do i=1,3
         cbeta(i)=0.0
      enddo
      return

 100  continue

      do i=1,3
         axes(i,1)=ccarbon(i)-calpha(i)
         axes(i,2)=calpha(i)-cnitro(i)
      enddo
      call cros(axes(1,1),axes(1,2),axes(1,3))   !generate rotation matrix
      call cros(axes(1,3),axes(1,1),axes(1,2))  !rosetta function cros
      call unitvec(axes(1,1),axes(1,1))
      call unitvec(axes(1,3),axes(1,3))
      call unitvec(axes(1,2),axes(1,2))

      do i=1,3          !template coord x rotation matrix + base coord
         cbeta(i)=calpha(i) + CBpos(1) * axes(i,1) +
     #        CBpos(2) * axes(i,2) +
     #        CBpos(3) * axes(i,3)
      enddo

      return
      end
             
              
