c===============================================================================
c
      SUBROUTINE sort(n,ir,x,nam)

c-------------------------------------------------------------------------------
C  Modified from routine SORT3 from Numerical Recipes Software 
c-------------------------------------------------------------------------------

      implicit none

      integer      maxres

      parameter  (maxres = 3000)

      integer      j,n,k,ktmp,iwksp(maxres)
      integer      ir(maxres),irtmp(maxres)
      character*3  nam(maxres),namtmp(maxres)
      real*8       x(3*maxres),xtmp(3*maxres)
      real         wksp(maxres)
      
      k = 0
      do j = 1, n
        wksp(j)   = float(ir(j))
        irtmp(j)  = ir(j)
        namtmp(j) = nam(j)
        xtmp(k+1) = x(k+1)
        xtmp(k+2) = x(k+2)
        xtmp(k+3) = x(k+3)
        k=k+3
      enddo

      call indexx(n,wksp,iwksp)

      k = 0
      do j = 1, n
        ktmp   = 3*(iwksp(j)-1)
        ir(j)  = irtmp(iwksp(j))
        nam(j) = namtmp(iwksp(j))
        x(k+1) = xtmp(ktmp+1)
        x(k+2) = xtmp(ktmp+2)
        x(k+3) = xtmp(ktmp+3)
        k=k+3
      enddo
        
      return
      end
c
c===============================================================================
c
      SUBROUTINE indexx(n,arr,indx)

c-------------------------------------------------------------------------------
C  (C) Copr. 1986-92 Numerical Recipes Software
c-------------------------------------------------------------------------------

      parameter  (maxres = 3000)
      INTEGER n,indx(maxres),M,NSTACK
      real arr(maxres)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      real a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else

        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
c
c===============================================================================
