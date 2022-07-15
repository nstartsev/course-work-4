*
      subroutine com1(epsy)
      parameter (Jmax=128)
      parameter (Kmax=128)
      implicit real*8 (a-h,o-z)
      common
     >/dimy/yn(0:1024),yn1(0:1024),ym(0:1025),ym1(0:1025),Jm
     >/dimz/zn(0:1024),zn1(0:1024),zm(0:1025),zm1(0:1025),Km
     >/alre/al,be, Re,v0
     >/basefl/ub(0:257,0:257),Oxb(0:257,0:257),vb(0:257,0:257)
     >,Oyb(0:257,0:257)
     >,wb(0:257,0:257),Ozb(0:257,0:257)
     >/eigen1/ur(0:257,0:257),vr(0:257,0:257),wr(0:257,0:257),
     > ui(0:257,0:257),vi(0:257,0:257),wi(0:257,0:257)
      character(len = 50)
     > ueig1, veig1, weig1, alchar,bechar, ueig2, veig2, weig2, ceig,
     > urealchar, uimagechar, v0char, fname, rechar,
     > vrealchar, vimagechar,
     > wrealchar, wimagechar,
     > prealchar, pimagechar,streamfunchar,tmpchar
         ci=(0.d0,1.d0)
	     alci=al*ci
	     cial=alci
*

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

        hz=4*pi/be/Km
        do k=0,Km
          zn(k)=hz*k-2*pi/be
          zn1(k)=hz
        end do
        do k=0,Km+1
          zm(k)=(k-0.5d0)*hz-2*pi/be
          zm1(k)=hz
        end do


* Base flow

10    format(F10.3)
300   format('  ZONE  I=',i4,'  J=',i4,'  F=POINT')
        write(alchar,10) al
        alchar = adjustl(alchar)
        write(bechar,10) be
        bechar = adjustl(bechar)
        write(v0char,10) v0
        v0char = adjustl(v0char)
          write(rechar,10) Re
          rechar = adjustl(rechar)

      fname = 'produbase_v0_'//
     > trim(v0char)//'_be_'//trim(bechar)//'_Re_'//trim(rechar)//'.dat'
      open(4, File = fname)
          do j = 1, Jm
            do k = Km/2+1, Km
                 read(4,*) a,b,temp
                 ub(j,k) = temp
            end do
        end do

        do k = 0, Km/2
            do j = 0, Jm+1
                ub(j,k) = ub(j,Km+1-k)
            end do
        end do
        close(4)

        do k=0,Km+1
            ub(0,k)=-ub(1,k)-2
            ub(Jm+1,k)=-ub(Jm,k)+2
        end do

        do j=0,Jm+1
            ub(j,0)=ub(j,1)
            ub(j,Km+1)=ub(j,Km)
        end do



      fname = 'prodvbase_v0_'//
     > trim(v0char)//'_be_'//trim(bechar)//'_Re_'//trim(rechar)//'.dat'
      open(4, File=fname)
      fname = 'prodwbase_v0_'//
     > trim(v0char)//'_be_'//trim(bechar)//'_Re_'//trim(rechar)//'.dat'
      open(5, File=fname)


      do j = 0, Jm
        do k = Km/2+1, Km
            read(4,*) a,b , vb(j,k)
        end do
      end do


        do k = 0, Km/2
            do j = 0, Jm+1
                vb(j,k) = vb(j,Km+1-k)
            end do
        end do


        do j=0,Jm+1
            vb(j,0)=vb(j,1)
            vb(j,Km+1)=vb(j,Km)
        end do


      do j = 1, Jm
        do k = Km/2,Km
            read(5,*) a,b , wb(j,k)
        end do
      end do

        do k = 0, Km/2
            do j = 0, Jm+1
                wb(j,k) = -wb(j,Km-k)
            end do
        end do

        do k=0,Km+1
            wb(0,k)=   -wb(1,k)
            wb(Jm+1,k)=-wb(Jm,k)
        end do




      close(4)
      close(5)

        do j=0,Jm
            do k=0,Km
                oxb(j,k)=-(vb(j,k+1)-vb(j,k))/zn1(k)
     >          +(wb(j+1,k)-wb(j,k))/yn1(j)
            end do
        end do
        do j=0,Jm+1
            do k=0,Km
                oyb(j,k)=(ub(j,k+1)-ub(j,k))/zn1(k)
     >          -cial*wb(j,k)
            end do
        end do
        do j=0,Jm
            do k=0,Km+1
                ozb(j,k)=-(ub(j+1,k)-ub(j,k))/yn1(j)
     >          +cial*vb(j,k)
            end do
        end do
*
*     writing u, v, w components of first eigenmode

      ueig1 = 'ueigen1_al_'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
      veig1 = 'veigen1_al_'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
      weig1 = 'weigen1_al_'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'


      open(4, File = ueig1)
      open(5, File = veig1)
      open(7, File = weig1)
      open(8, File = 'ch.dat')
      write(8,300)  Km/2, Jm
      do j = 1 , Jm
        do k = Km/2+1, Km/2+Km/2
            read(4,*) a ,b, ur(j,k),ui(j,k)
            write(8,*) zm(k), ym(j), ui(j,k)
        end do
      end do


      do j = 1, Jm-1
        do k = 1+Km/2, Km/2+Km/2
            read(5,*) a,b, vr(j,k),vi(j,k)
        end do
      end do

       do j = 1, Jm
        do k = 0+Km/2 , Km/2+Km/2
            read(7,*)  a,b, wr(j,k),wi(j,k)
        end do
       end do

      close(4)
      close(5)
      close(7)


      do j = 0, Jm+1
        do k = 0, Km/2
            ur(j,k) = -ur(j, Km+1 - k)
            ui(j,k) = -ui(j, Km+1 - k)
        end do
      end do


        do k = 1,Km/2
            do j = 0, Jm
                vr(j,k) = -vr(j,Km+1-k)
                vi(j,k) = -vi(j,Km+1-k)
            end do
        end do


        do k = 0,Km/2
            do j = 0, Jm+1
                wr(j,k) = wr(j,Km-k)
                wi(j,k) = wi(j,Km-k)
            end do
        end do

* boundary conditions
      do k=0,Km
        wr(0,k)   = -wr(1,k)
        wr(Jm+1,k)= -wr(Jm,k)
        wi(0,k)   = -wi(1,k)
        wi(Jm+1,k)= -wi(Jm,k)
      end do

      do j=0,Jm
        vr(j,0)    = -vr(j,1)
        vr(j,Km+1) = -vr(j,Km)
        vi(j,0)    = -vi(j,1)
        vi(j,Km+1) = -vi(j,Km)
      end do

      do k = 0, Km+1
        vr(0,k) = 0
        vi(0,k) = 0
        vr(Jm,k) = 0
        vi(Jm,k) = 0
      end do

      do k=0,Km+1
        ur(0,k)   = -ur(1,k)
        ur(Jm+1,k)= -ur(Jm,k)
        ui(0,k)   = -ui(1,k)
        ui(Jm+1,k)= -ui(Jm,k)
      end do


      do j=0,Jm+1
        ur(j,0)   =  -ur(j,1)
        ur(j,Km+1)=  -ur(j,Km)
        ui(j,0)   =  -ui(j,1)
        ui(j,Km+1)=  -ui(j,Km)
      end do
*


      return
      end
*
      subroutine com2(epsy)
      parameter (Jmax=128)
      parameter (Kmax=128)
      implicit real*8 (a-h,o-z)
      common
     >/dimy/yn(0:1024),yn1(0:1024),ym(0:1025),ym1(0:1025),Jm
     >/dimz/zn(0:1024),zn1(0:1024),zm(0:1025),zm1(0:1025),Km
     >/alre/al,be, Re,v0
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
     > prealchar, pimagechar,streamfunchar
         ci=(0.d0,1.d0)
	     alci=al*ci
	     cial=alci
*

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

        hz=4*pi/be/Km
        do k=0,Km
          zn(k)=hz*k-2*pi/be
          zn1(k)=hz
        end do
        do k=0,Km+1
          zm(k)=(k-0.5d0)*hz-2*pi/be
          zm1(k)=hz
        end do


* Base flow

10    format(F10.3)
        write(alchar,10) al
        alchar = adjustl(alchar)
        write(bechar,10) be
        bechar = adjustl(bechar)
        write(v0char,10) v0
        v0char = adjustl(v0char)

      streamfunchar = 'streamfun'//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
      open(4, File = streamfunchar)
        do k = Km/2+1, Km+1
            do j = 0, Jm+1
                 read(4,*) a,b,temp
                 ub(j,k) = temp
            end do
        end do

        do k = 0, Km/2
            do j = 0, Jm+1
                ub(j,k) = ub(j,Km+1-k)
            end do
        end do
        close(4)

        do k=0,Km+1
            ub(0,k)=-ub(1,k)-2
            ub(Jm+1,k)=-ub(Jm,k)+2
        end do

        do j=0,Jm+1
            ub(j,0)=ub(j,1)
            ub(j,Km+1)=ub(j,Km)
        end do



      fname = 'v_v0_'//trim(v0char)//'_be_'//trim(bechar)//'.dat'
      open(4, File=fname)
      fname = 'w_v0_'//trim(v0char)//'_be_'//trim(bechar)//'.dat'
      open(5, File=fname)

*      write(4,300)Jm+1,Km/2

      do k = 1, Km
        do j = 0, Jm
            read(4,*) a,b , vb(j,k)
        end do
      end do
*      write(5,300)Jm,Km/2+1

      do k = 0,Km
        do j = 1, Jm
            read(5,*) a,b , wb(j,k)
        end do
      end do

      close(4)
      close(5)

        do j=0,Jm
            do k=0,Km
                oxb(j,k)=-(vb(j,k+1)-vb(j,k))/zn1(k)
     >          +(wb(j+1,k)-wb(j,k))/yn1(j)
            end do
        end do
        do j=0,Jm+1
            do k=0,Km
                oyb(j,k)=(ub(j,k+1)-ub(j,k))/zn1(k)
     >          -cial*wb(j,k)
            end do
        end do
        do j=0,Jm
            do k=0,Km+1
                ozb(j,k)=-(ub(j+1,k)-ub(j,k))/yn1(j)
     >          +cial*vb(j,k)
            end do
        end do
*


      urealchar = 'ureal2(y,z)'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
      uimagechar ='uimage2(y,z)'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
      open(4, File = urealchar)
      open(5, File = uimagechar)
      do k = Km/2, Km+1
        do j = 0, Jm+1
            read(4,*) a,b, ur(j,k)
            read(5,*) a,b, ui(j,k)
        end do
      end do

      do j = 0, Jm+1
        do k = 0, Km/2
            ur(j,k) = -ur(j, Km+1 - k)
            ui(j,k) = -ui(j, Km+1 - k)
        end do
      end do

      vrealchar = 'vreal2(y,z)'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
      vimagechar ='vimage2(y,z)'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
      open(4, File = vrealchar)
      open(5, File = vimagechar)

      do k = Km/2+1 , Km
        do j = 0, Jm
            read(4,*) a,b, vr(j,k)
            read(5,*) a,b, vi(j,k)
          end do
      end do
      close(4)
      close(5)

        do k = 1,Km/2
            do j = 0, Jm
                vr(j,k) = -vr(j,Km+1-k)
                vi(j,k) = -vi(j,Km+1-k)
            end do
        end do


      wrealchar = 'wreal2(y,z)'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
      wimagechar ='wimage2(y,z)'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
      open(4, File = wrealchar)
      open(5, File = wimagechar)
        do k = Km/2, Km
            do j = 0, Jm+1
                read(4,*) a,b, wr(j,k)
                read(5,*) a,b, wi(j,k)
            end do
        end do

        do k = 0,Km/2-1
            do j = 0, Jm+1
                wr(j,k) = wr(j,Km-k)
                wi(j,k) = wi(j,Km-k)
            end do
        end do

* boundary conditions
      do k=0,Km
        wr(0,k)   = -wr(1,k)
        wr(Jm+1,k)= -wr(Jm,k)
        wi(0,k)   = -wi(1,k)
        wi(Jm+1,k)= -wi(Jm,k)
      end do

      do j=0,Jm
        vr(j,0)    = -vr(j,1)
        vr(j,Km+1) = -vr(j,Km)
        vi(j,0)    = -vi(j,1)
        vi(j,Km+1) = -vi(j,Km)
      end do

      do k = 0, Km+1
        vr(0,k) = 0
        vi(0,k) = 0
        vr(Jm,k) = 0
        vi(Jm,k) = 0
      end do

      do k=0,Km+1
        ur(0,k)   = -ur(1,k)
        ur(Jm+1,k)= -ur(Jm,k)
        ui(0,k)   = -ui(1,k)
        ui(Jm+1,k)= -ui(Jm,k)
      end do


      do j=0,Jm+1
        ur(j,0)   =  -ur(j,1)
        ur(j,Km+1)=  -ur(j,Km)
        ui(j,0)   =  -ui(j,1)
        ui(j,Km+1)=  -ui(j,Km)
      end do
*


      return
      end
*
      subroutine com3(epsy)
      parameter (Jmax=128)
      parameter (Kmax=128)
      implicit real*8 (a-h,o-z)
      common
     >/dimy/yn(0:1024),yn1(0:1024),ym(0:1025),ym1(0:1025),Jm
     >/dimz/zn(0:1024),zn1(0:1024),zm(0:1025),zm1(0:1025),Km
     >/alre/al,be, Re,v0
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
     > prealchar, pimagechar,streamfunchar
         ci=(0.d0,1.d0)
	     alci=al*ci
	     cial=alci
*

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

        hz=4*pi/be/Km
        do k=0,Km
          zn(k)=hz*k-2*pi/be
          zn1(k)=hz
        end do
        do k=0,Km+1
          zm(k)=(k-0.5d0)*hz-2*pi/be
          zm1(k)=hz
        end do


* Base flow

10    format(F10.3)
        write(alchar,10) al
        alchar = adjustl(alchar)
        write(bechar,10) be
        bechar = adjustl(bechar)
        write(v0char,10) v0
        v0char = adjustl(v0char)

      streamfunchar = 'streamfun'//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
      open(4, File = streamfunchar)
        do k = Km/2+1, Km+1
            do j = 0, Jm+1
                 read(4,*) a,b,temp
                 ub(j,k) = temp
            end do
        end do

        do k = 0, Km/2
            do j = 0, Jm+1
                ub(j,k) = ub(j,Km+1-k)
            end do
        end do
        close(4)

        do k=0,Km+1
            ub(0,k)=-ub(1,k)-2
            ub(Jm+1,k)=-ub(Jm,k)+2
        end do

        do j=0,Jm+1
            ub(j,0)=ub(j,1)
            ub(j,Km+1)=ub(j,Km)
        end do



      fname = 'v_v0_'//trim(v0char)//'_be_'//trim(bechar)//'.dat'
      open(4, File=fname)
      fname = 'w_v0_'//trim(v0char)//'_be_'//trim(bechar)//'.dat'
      open(5, File=fname)

*      write(4,300)Jm+1,Km/2

      do k = 1, Km
        do j = 0, Jm
            read(4,*) a,b , vb(j,k)
        end do
      end do
*      write(5,300)Jm,Km/2+1

      do k = 0,Km
        do j = 1, Jm
            read(5,*) a,b , wb(j,k)
        end do
      end do

      close(4)
      close(5)

        do j=0,Jm
            do k=0,Km
                oxb(j,k)=-(vb(j,k+1)-vb(j,k))/zn1(k)
     >          +(wb(j+1,k)-wb(j,k))/yn1(j)
            end do
        end do
        do j=0,Jm+1
            do k=0,Km
                oyb(j,k)=(ub(j,k+1)-ub(j,k))/zn1(k)
     >          -cial*wb(j,k)
            end do
        end do
        do j=0,Jm
            do k=0,Km+1
                ozb(j,k)=-(ub(j+1,k)-ub(j,k))/yn1(j)
     >          +cial*vb(j,k)
            end do
        end do
*


      urealchar = 'ureal3(y,z)'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
      uimagechar ='uimage3(y,z)'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
      open(4, File = urealchar)
      open(5, File = uimagechar)
      do k = Km/2, Km+1
        do j = 0, Jm+1
            read(4,*) a,b, ur(j,k)
            read(5,*) a,b, ui(j,k)
        end do
      end do

      do j = 0, Jm+1
        do k = 0, Km/2
            ur(j,k) = -ur(j, Km+1 - k)
            ui(j,k) = -ui(j, Km+1 - k)
        end do
      end do

      vrealchar = 'vreal3(y,z)'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
      vimagechar ='vimage3(y,z)'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
      open(4, File = vrealchar)
      open(5, File = vimagechar)

      do k = Km/2+1 , Km
        do j = 0, Jm
            read(4,*) a,b, vr(j,k)
            read(5,*) a,b, vi(j,k)
          end do
      end do
      close(4)
      close(5)

        do k = 1,Km/2
            do j = 0, Jm
                vr(j,k) = -vr(j,Km+1-k)
                vi(j,k) = -vi(j,Km+1-k)
            end do
        end do


      wrealchar = 'wreal3(y,z)'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
      wimagechar ='wimage3(y,z)'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
      open(4, File = wrealchar)
      open(5, File = wimagechar)
        do k = Km/2, Km
            do j = 0, Jm+1
                read(4,*) a,b, wr(j,k)
                read(5,*) a,b, wi(j,k)
            end do
        end do

        do k = 0,Km/2-1
            do j = 0, Jm+1
                wr(j,k) = wr(j,Km-k)
                wi(j,k) = wi(j,Km-k)
            end do
        end do

* boundary conditions
      do k=0,Km
        wr(0,k)   = -wr(1,k)
        wr(Jm+1,k)= -wr(Jm,k)
        wi(0,k)   = -wi(1,k)
        wi(Jm+1,k)= -wi(Jm,k)
      end do

      do j=0,Jm
        vr(j,0)    = -vr(j,1)
        vr(j,Km+1) = -vr(j,Km)
        vi(j,0)    = -vi(j,1)
        vi(j,Km+1) = -vi(j,Km)
      end do

      do k = 0, Km+1
        vr(0,k) = 0
        vi(0,k) = 0
        vr(Jm,k) = 0
        vi(Jm,k) = 0
      end do

      do k=0,Km+1
        ur(0,k)   = -ur(1,k)
        ur(Jm+1,k)= -ur(Jm,k)
        ui(0,k)   = -ui(1,k)
        ui(Jm+1,k)= -ui(Jm,k)
      end do


      do j=0,Jm+1
        ur(j,0)   =  -ur(j,1)
        ur(j,Km+1)=  -ur(j,Km)
        ui(j,0)   =  -ui(j,1)
        ui(j,Km+1)=  -ui(j,Km)
      end do
*


      return
      end
*
      subroutine com4(epsy)
      parameter (Jmax=128)
      parameter (Kmax=128)
      implicit real*8 (a-h,o-z)
      common
     >/dimy/yn(0:1024),yn1(0:1024),ym(0:1025),ym1(0:1025),Jm
     >/dimz/zn(0:1024),zn1(0:1024),zm(0:1025),zm1(0:1025),Km
     >/alre/al,be, Re,v0
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
     > prealchar, pimagechar,streamfunchar
         ci=(0.d0,1.d0)
	     alci=al*ci
	     cial=alci
*

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

        hz=4*pi/be/Km
        do k=0,Km
          zn(k)=hz*k-2*pi/be
          zn1(k)=hz
        end do
        do k=0,Km+1
          zm(k)=(k-0.5d0)*hz-2*pi/be
          zm1(k)=hz
        end do


* Base flow

10    format(F10.3)
        write(alchar,10) al
        alchar = adjustl(alchar)
        write(bechar,10) be
        bechar = adjustl(bechar)
        write(v0char,10) v0
        v0char = adjustl(v0char)

      streamfunchar = 'streamfun'//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
      open(4, File = streamfunchar)
        do k = Km/2+1, Km+1
            do j = 0, Jm+1
                 read(4,*) a,b,temp
                 ub(j,k) = temp
            end do
        end do

        do k = 0, Km/2
            do j = 0, Jm+1
                ub(j,k) = ub(j,Km+1-k)
            end do
        end do
        close(4)

        do k=0,Km+1
            ub(0,k)=-ub(1,k)-2
            ub(Jm+1,k)=-ub(Jm,k)+2
        end do

        do j=0,Jm+1
            ub(j,0)=ub(j,1)
            ub(j,Km+1)=ub(j,Km)
        end do



      fname = 'v_v0_'//trim(v0char)//'_be_'//trim(bechar)//'.dat'
      open(4, File=fname)
      fname = 'w_v0_'//trim(v0char)//'_be_'//trim(bechar)//'.dat'
      open(5, File=fname)

*      write(4,300)Jm+1,Km/2

      do k = 1, Km
        do j = 0, Jm
            read(4,*) a,b , vb(j,k)
        end do
      end do
*      write(5,300)Jm,Km/2+1

      do k = 0,Km
        do j = 1, Jm
            read(5,*) a,b , wb(j,k)
        end do
      end do

      close(4)
      close(5)

        do j=0,Jm
            do k=0,Km
                oxb(j,k)=-(vb(j,k+1)-vb(j,k))/zn1(k)
     >          +(wb(j+1,k)-wb(j,k))/yn1(j)
            end do
        end do
        do j=0,Jm+1
            do k=0,Km
                oyb(j,k)=(ub(j,k+1)-ub(j,k))/zn1(k)
     >          -cial*wb(j,k)
            end do
        end do
        do j=0,Jm
            do k=0,Km+1
                ozb(j,k)=-(ub(j+1,k)-ub(j,k))/yn1(j)
     >          +cial*vb(j,k)
            end do
        end do
*


      urealchar = 'ureal4(y,z)'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
      uimagechar ='uimage4(y,z)'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
      open(4, File = urealchar)
      open(5, File = uimagechar)
      do k = Km/2, Km+1
        do j = 0, Jm+1
            read(4,*) a,b, ur(j,k)
            read(5,*) a,b, ui(j,k)
        end do
      end do

      do j = 0, Jm+1
        do k = 0, Km/2
            ur(j,k) = -ur(j, Km+1 - k)
            ui(j,k) = -ui(j, Km+1 - k)
        end do
      end do

      vrealchar = 'vreal4(y,z)'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
      vimagechar ='vimage4(y,z)'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
      open(4, File = vrealchar)
      open(5, File = vimagechar)

      do k = Km/2+1 , Km
        do j = 0, Jm
            read(4,*) a,b, vr(j,k)
            read(5,*) a,b, vi(j,k)
          end do
      end do
      close(4)
      close(5)

        do k = 1,Km/2
            do j = 0, Jm
                vr(j,k) = -vr(j,Km+1-k)
                vi(j,k) = -vi(j,Km+1-k)
            end do
        end do


      wrealchar = 'wreal4(y,z)'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
      wimagechar ='wimage4(y,z)'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
      open(4, File = wrealchar)
      open(5, File = wimagechar)
        do k = Km/2, Km
            do j = 0, Jm+1
                read(4,*) a,b, wr(j,k)
                read(5,*) a,b, wi(j,k)
            end do
        end do

        do k = 0,Km/2-1
            do j = 0, Jm+1
                wr(j,k) = wr(j,Km-k)
                wi(j,k) = wi(j,Km-k)
            end do
        end do

* boundary conditions
      do k=0,Km
        wr(0,k)   = -wr(1,k)
        wr(Jm+1,k)= -wr(Jm,k)
        wi(0,k)   = -wi(1,k)
        wi(Jm+1,k)= -wi(Jm,k)
      end do

      do j=0,Jm
        vr(j,0)    = -vr(j,1)
        vr(j,Km+1) = -vr(j,Km)
        vi(j,0)    = -vi(j,1)
        vi(j,Km+1) = -vi(j,Km)
      end do

      do k = 0, Km+1
        vr(0,k) = 0
        vi(0,k) = 0
        vr(Jm,k) = 0
        vi(Jm,k) = 0
      end do

      do k=0,Km+1
        ur(0,k)   = -ur(1,k)
        ur(Jm+1,k)= -ur(Jm,k)
        ui(0,k)   = -ui(1,k)
        ui(Jm+1,k)= -ui(Jm,k)
      end do


      do j=0,Jm+1
        ur(j,0)   =  -ur(j,1)
        ur(j,Km+1)=  -ur(j,Km)
        ui(j,0)   =  -ui(j,1)
        ui(j,Km+1)=  -ui(j,Km)
      end do
*


      return
      end
