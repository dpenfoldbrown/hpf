c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 2338 $
c  $Date: 2002-08-09 18:23:11 -0400 (Fri, 09 Aug 2002) $
c  $Author: rohl $

      subroutine convert3_num  
      implicit none 
c-------------------------------------------------------
c  THIS SUBROUTINE CONVERTS THE LETTER CODES INTO
c    1letter CODES 
c  
c   ONE LETTER CODE   THREE LETTER CODE   NUMBER CODE 
c         A                   ALA              1
c         C                   CYS              2
c         D                   ASP              3
c         E                   GLU              4
c         F                   PHE              5
c         G                   GLY              6
c         H                   HIS              7
c         I                   ILE              8
c         K                   LYS              9
c         L                   LEU             10
c         M                   MET             11
c         N                   ASN             12
c         P                   PRO             13
c         Q                   GLN             14
c         R                   ARG             15 
c         S                   SER             16
c         T                   THR             17
c         V                   VAL             18
c         W                   TRP             19
c         Y                   TYR             20 
c------------------------------------------------------
      include 'param.h'
      include 'structure.h' 
      integer i 
c------------------------------------------------------
      do 100 i=1,total_residue  
         if (residue3(i).eq.'ALA') then 
            residue1(i)='A' 
            goto 100 
         endif 
         if (residue3(i).eq.'CYS') then 
            residue1(i)='C' 
            goto 100 
         endif 
         if (residue3(i).eq.'ASP'.or.residue3(i).eq.'ASX') then 
            residue1(i)='D' 
            goto 100 
         endif 
         if (residue3(i).eq.'GLU'.or.residue3(i).eq.'GLX') then 
            residue1(i)='E' 
            goto 100 
         endif 
         if (residue3(i).eq.'PHE') then 
            residue1(i)='F' 
            goto 100 
         endif 
         if (residue3(i).eq.'GLY') then 
            residue1(i)='G' 
            goto 100 
         endif 
         if (residue3(i).eq.'HIS') then 
            residue1(i)='H' 
            goto 100 
         endif 
         if (residue3(i).eq.'ILE') then 
            residue1(i)='I' 
            goto 100 
         endif 
         if (residue3(i).eq.'LYS') then 
            residue1(i)='K' 
            goto 100 
         endif 
         if (residue3(i).eq.'LEU') then 
            residue1(i)='L' 
            goto 100 
         endif 
         if (residue3(i).eq.'MET') then 
            residue1(i)='M' 
            goto 100 
         endif 
         if (residue3(i).eq.'ASN') then 
            residue1(i)='N' 
            goto 100 
         endif 
         if (residue3(i).eq.'PRO') then 
            residue1(i)='P' 
            goto 100 
         endif 
         if (residue3(i).eq.'GLN') then 
            residue1(i)='Q' 
            goto 100 
         endif 
         if (residue3(i).eq.'ARG') then 
            residue1(i)='R' 
            goto 100 
         endif 
         if (residue3(i).eq.'SER') then 
            residue1(i)='S' 
            goto 100 
         endif 
         if (residue3(i).eq.'THR') then 
            residue1(i)='T' 
            goto 100 
         endif 
         if (residue3(i).eq.'VAL') then 
            residue1(i)='V' 
            goto 100 
         endif 
         if (residue3(i).eq.'TRP') then 
            residue1(i)='W' 
            goto 100 
         endif 
         if (residue3(i).eq.'TYR') then      
            residue1(i)='Y' 
            goto 100 
         endif 
         if (residue3(i).ne.'   ') then 
              write(60,*)'cannot handle ',residue3(i)  
         endif 
         residue1(i)=' ' 
 100  continue 
      return 
      end 
