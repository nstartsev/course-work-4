      subroutine check1()
      parameter (Jmax=128)
      parameter (Kmax=128)
      implicit real*8 (a-h,o-z)
      common
     >/dimy/yn(0:1024),yn1(0:1024),ym(0:1025),ym1(0:1025),Jm
     >/dimz/zn(0:1024),zn1(0:1024),zm(0:1025),zm1(0:1025),Km
     >/alre/al,be, Re, v0
     >/basefl/ub(0:257,0:257),Oxb(0:257,0:257),vb(0:257,0:257)
     >,Oyb(0:257,0:257)
     >,wb(0:257,0:257),Ozb(0:257,0:257)
     >/eigen1/ur(0:257,0:257),vr(0:257,0:257),wr(0:257,0:257),
     > ui(0:257,0:257),vi(0:257,0:257),wi(0:257,0:257)
       character(len = 50)
     > ueig1, veig1, weig1, alchar,bechar, ueig2, veig2, weig2, ceig,
     > urealchar, uimagechar, v0char, fname,
     > vrealchar, vimagechar,
     > wrealchar, wimagechar,
     > uchar, vchar, wchar,
     > prealchar, pimagechar,streamfunchar
         ci=(0.d0,1.d0)
	     alci=al*ci
	     cial=alci



10    format(F10.3)
        write(alchar,10) al
        alchar = adjustl(alchar)
        write(bechar,10) be
        bechar = adjustl(bechar)
        write(v0char,10) v0
        v0char = adjustl(v0char)

      open(4,file='check_ucouett.dat')
      write(4,300) Jm,Km
      do k = 1 , Km
        do j = 1, Jm
            write(4,*) ym(j), zm(k), ym(j)
        end do
      end do
      open(4, File='check_ubase.dat')
      write(4,300)Jm+1,Km+2
        do k = 0, Km+1
            do j = 0, Jm
                write(4,*) yn(j), zm(k), 0.5d0*(ub(j,k)+ub(j+1,k))
            end do
        end do
      close(4)
*

      open(4,File='check_vbase.dat')
      write(4,300) Jm+1, Km
        do k = 1, Km
            do j = 0, Jm
                write(4,*) yn(j), zm(k) , vb(j,k)

            end do
        end do
      close(4)

      open(4,File='check_vwbase.dat')
      write(4,300) Jm/2+1, Km/2+1
        do k = 0, Km,2
            do j = 0, Jm ,2
                tmp1 = 0.5d0*(vb(j,k)+vb(j,k+1))
                tmp2 = 0.5d0*(wb(j,k)+wb(j+1,k))
                write(4,*) yn(j),zn(k), tmp1, tmp2
            end do
        end do
        close(4)
      open(5,File='check_wbase.dat')
      write(5,300) Jm,Km+1
        do k = 0, Km
            do j = 1, Jm

                write(5,*) ym(j), zn(k) , wb(j,k)
            end do
        end do

      close(5)

      open(4, File='check_ureigen.dat')
      write(4,300) Jm+1, Km+2
      do k = 0 , Km+1
        do j = 0 , Jm
            write(4,*) yn(j), zm(k), 0.5d0*(ur(j,k)+ur(j+1,k))
        end do
      end do
      close(4)

      open(4, File='check_uieigen.dat')
      write(4,300) Jm+1, Km+2
      do k = 0 , Km+1
        do j = 0 , Jm
            write(4,*) yn(j), zm(k), 0.5d0*(ui(j,k)+ui(j+1,k))
        end do
      end do
      close(4)


        uchar = 'ueigen11'//'_al_'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'


      open(4, File=uchar)


      write(4,300) Jm+2, Km+2
      do k = 0 , Km+1
        do j = 0 , Jm+1
            tmp = ur(j,k)**2+ui(j,k)**2
            write(4,*) ym(j), zm(k), tmp
        end do
      end do
      close(4)


      open(4, File='check_vreigen.dat')
      write(4,300) Jm+1, Km+2
      do k = 0 , Km+1
        do j = 0 , Jm
            write(4,*) yn(j), zm(k), vr(j,k)
        end do
      end do
      close(4)

      open(4, File='check_vieigen.dat')
      write(4,300) Jm+1, Km+2
      do k = 0 , Km+1
        do j = 0 , Jm
            write(4,*) yn(j), zm(k), vi(j,k)
        end do
      end do
      close(4)


        vchar = 'veigen11'//'_al_'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'


      open(4, File=vchar)
      write(4,300) Jm+1, Km+2
      do k = 0 , Km+1
        do j = 0 , Jm
            tmp = (vi(j,k))**2+(vr(j,k))**2
            write(4,*) yn(j), zm(k), tmp
        end do
      end do
      close(4)


      open(4, File='check_wreigen.dat')
      write(4,300) Jm+2, Km+1
      do k = 0 , Km
        do j = 0 , Jm+1
            write(4,*) ym(j), zn(k), wr(j,k)
        end do
      end do
      close(4)

      open(4, File='check_wieigen.dat')
      write(4,300) Jm+2, Km+1
      do k = 0 , Km
        do j = 0 , Jm+1
            write(4,*) ym(j), zn(k), wi(j,k)
        end do
      end do
      close(4)
        wchar = 'weigen11'//'_al_'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
      open(4, File=wchar)
      write(4,300) Jm+2, Km+1
      do k = 0 , Km
        do j = 0 , Jm+1
            tmp = (wi(j,k))**2+(wr(j,k))**2
            write(4,*) ym(j), zn(k), tmp
        end do
      end do
      close(4)

300   format('  ZONE  I=',i4,'  J=',i4,'  F=POINT')


      return
      end
