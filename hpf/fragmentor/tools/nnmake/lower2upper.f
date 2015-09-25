c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 2072 $
c  $Date: 2002-06-18 00:41:53 -0400 (Tue, 18 Jun 2002) $
c  $Author: rohl $

      function lower2upper(char1) 
      implicit none 
c-------------------------------------------------------------
      character lower2upper*1,char1*1  
c-------------------------------------------------------------
      lower2upper=' '
      if (char1.eq.'a') lower2upper='A' 
      if (char1.eq.'b') lower2upper='B' 
      if (char1.eq.'c') lower2upper='C' 
      if (char1.eq.'d') lower2upper='D' 
      if (char1.eq.'e') lower2upper='E' 
      if (char1.eq.'f') lower2upper='F' 
      if (char1.eq.'g') lower2upper='G' 
      if (char1.eq.'h') lower2upper='H' 
      if (char1.eq.'i') lower2upper='I' 
      if (char1.eq.'j') lower2upper='J' 
      if (char1.eq.'k') lower2upper='K' 
      if (char1.eq.'l') lower2upper='L' 
      if (char1.eq.'m') lower2upper='M' 
      if (char1.eq.'n') lower2upper='N' 
      if (char1.eq.'o') lower2upper='O' 
      if (char1.eq.'p') lower2upper='P' 
      if (char1.eq.'q') lower2upper='Q' 
      if (char1.eq.'r') lower2upper='R' 
      if (char1.eq.'s') lower2upper='S' 
      if (char1.eq.'t') lower2upper='T' 
      if (char1.eq.'u') lower2upper='U' 
      if (char1.eq.'v') lower2upper='V' 
      if (char1.eq.'w') lower2upper='W' 
      if (char1.eq.'y') lower2upper='Y' 
      if (char1.eq.'z') lower2upper='Z' 

      return 
      end 

