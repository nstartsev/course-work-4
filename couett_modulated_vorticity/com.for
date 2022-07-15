
      subroutine com(epsy)
      parameter (Jmax=128)
      parameter (Nmax=32)
      implicit real*8 (a-h,o-z)
        real*8
     >  rwork(10000),dpr(Nmax,Jmax)
      character(len = 50)
     >  tmpchar,bechar, v0char, fname, rechar
      common
     >/dimy/yn(0:1024),yn1(0:1024),ym(0:1025),ym1(0:1025),Jm
     >/pry/apym(256),bpym(256),cpym(256),apyn(256),bpyn(256),cpyn(256)
     >/pyp/apyp(256),bpyp(256),cpyp(256),apzn(256), bpzn(256), cpzn(256)
     >/pw/pwork(10000)
     >/alren/al,be0, Re, Nm, v0
     >/basefl/ub(0:257,0:257), vb(0:257,0:257),wb(0:257,0:257),Mm
         ci=(0.d0,1.d0)
	     alci=al*ci
	     cial=alci
*
*        be = 5.d0/3
        be = be0
        pi = 4*atan(1.)
        hy=2.d0/Jm
        do j=0,Jm
          yn(j)=hy*j-1
          yn1(j)=hy
        end do
        do j=0,Jm+1
          ym(j)=(j-0.5d0)*hy-1
          ym1(j)=hy
        end do


* race coefficients

        do j=1,Jm
         apym(j) = 1.d0/(ym1(j)*yn1(j-1))
         cpym(j) = 1.d0/(ym1(j)*yn1(j))
         bpym(j) = -apym(j)- cpym(j)
        end do

        do j = 1, Jm-1
            apyn(j) = 1.d0/(yn1(j)*ym1(j))
            cpyn(j) = 1.d0/(yn1(j)*ym1(j+1))
            bpyn(j) = -apyn(j) - cpyn(j)
        end do
*

* pressure coefficients
      do j = 1, Jm
        apyp(j) = apym(j)
        bpyp(j) = bpym(j)-al**2
        cpyp(j) = cpym(j)
      end do


        bpyp(1)=   apyp(1)+bpyp(1)
        apyp(1) = 0.
        bpyp(Jm) = bpyp(Jm)+cpyp(Jm)
        cpyp(Jm) = 0.

      do n = 1, Nm
        apzn(n) = 0.d0
        bpzn(n) = -((n-1)*be)**2
        cpzn(n) = 0.d0
      end do
* Base flow


*        do j = 1, Jm
*            ub(0,j) = ym(j)
*        end do
10    format(F10.3)
      write(bechar,10) be0
      bechar = adjustl(bechar)
      write(v0char,10) v0
      v0char = adjustl(v0char)

      write(rechar,10) re
      rechar = adjustl(rechar)

      fname = 'ubase_v0_'//
     > trim(v0char)//'_be_'//trim(bechar)//'_Re_'//trim(rechar)//'.dat'

        open(5, file=fname)
        do m=0, Mm-1
            do j = 0, Jm+1
                read(5,*) ub(m,j)
            end do
        end do
        close(5)
        do m = 0 , 0
            do j = 0 , Jm+1
*                 ub(m,j) = -ub(m, Jm+1-j)
            end do
        end do

        do m = 1 , Mm-1
            do j = 0 , Jm+1
*                 ub(m,j) = ub(m, Jm+1-j)
            end do
        end do
*        do j = 1, Jm
*            ub(0,j) = ym(j)
*        end do
        ub(0,0)    = -ub(0,1) - 2
        ub(0,Jm+1) = -ub(0,Jm) + 2

        do j=0, Jm+1
            wb(0,j) = 0.
        end do

        do  j = 0, Jm
            vb(0,j) = 0.
        end do
*
* fill pwork
        call blktri(0,1,Nm,apzn,bpzn,cpzn,1,Jm,apyp,bpyp,cpyp,Nmax,
     >               dpr,ierror,pwork)

*






* fill pwork
        call blktri(0,1,Nm,apzn,bpzn,cpzn,1,Jm,apyp,bpyp,cpyp,Nmax,
     >               dpr,ierror,pwork)

*

      return
      end
