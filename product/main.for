       implicit real*8 (a-h,o-z)
      parameter (Jmax=64)
      parameter (Kmax=64)

      common
     >/dimy/yn(0:1024),yn1(0:1024),ym(0:1025),ym1(0:1025),Jm
     >/dimz/zn(0:1024),zn1(0:1024),zm(0:1025),zm1(0:1025),Km
     >/alre/al,be, Re,v0
     >/basefl/ub(0:257,0:257),Oxb(0:257,0:257),vb(0:257,0:257)
     >,Oyb(0:257,0:257)
     >,wb(0:257,0:257),Ozb(0:257,0:257)
     >/eigen1/ur(0:257,0:257),vr(0:257,0:257),wr(0:257,0:257),
     > ui(0:257,0:257),vi(0:257,0:257),wi(0:257,0:257)

         ci=(0.d0,1.d0)
	     alci=al*ci
	     cial=alci
       open(5,file='oz.car')
        read(5,*)Jm
        read(5,*)Km
        read(5,*) al
        read(5,*) be
        read(5,*) Re
        read(5,*) eps
        read(5,*) v0

        close(5)

*        call com1(eps)
*        call product()
*        call check1()
*        stop

         open(20, file ='prod(al,be).dat')
         write(20,300) 50,50
        al = 0.0d0
        be = 0.0d0
        do j1 = 1, 50
            al = 0.0d0
            be = be+0.1d0
            do j2 = 1, 50
            al = al+0.1d0

        write(*,*) 'Jm=',Jm,'Km=', Km,'al=', al , 'be=', be ,'v0=',v0
        write(*,*) 'eigen1'
        call Com1(eps)
        call product()
        call check1()
            end do
        end do

        close(20)


*        write(*,*) 'eigen2'
*        call Com2(eps)
*        call product()
*        call check1()
*        write(*,*) 'eigen3'
*        call Com3(eps)
*        call product()


*        write(*,*) 'eigen4'
*        call Com4(eps)
*        call product()
*           write(*,*) 'eigen2'

*           call product()
* проверка векторных полей

300   format('  ZONE  I=',i4,'  J=',i4,'  F=POINT')
10    format(F10.3)
      end
