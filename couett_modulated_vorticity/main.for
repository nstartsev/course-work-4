       implicit real*8 (a-h,o-z)

      parameter (Jmax=128)
      parameter (Nmax=32)

      integer LWMAX
      parameter ( LWMAX = 100000 )
      parameter  (NJmax = 5000)
      real*8
     >rwork(LWMAX), zm(0:100), zn(0:100)
      complex*16
     > u(0:Nmax,0:Jmax),v(0:Nmax,0:Jmax),w(0:Nmax,0:Jmax)
     >,u1(0:Nmax,0:Jmax),v1(0:Nmax,0:Jmax),w1(0:Nmax,0:Jmax)
     >,p(0:Nmax,0:Jmax)
     >,ut(0:Nmax,0:Jmax),vt(0:Nmax,0:Jmax),wt(0:Nmax,0:Jmax)
     >,dp(Nmax,Jmax), tmp1, tmp2 ,tmp3
     >,Amat(Jmax*Nmax,Jmax*Nmax), evecl(Jmax*Nmax,Jmax*Nmax)
     >,evec(Jmax*Nmax,Jmax*Nmax),eval(Jmax*Nmax),work(LWMAX)
     >,Amat1(Jmax*Nmax,Jmax*Nmax)
     >,ci,alci,cial,angle
         real*8
     >  dpr(Nmax,Jmax), dpi(Nmax,Jmax)
      character(len = 50)
     > ueig1, veig1, weig1, alchar,bechar, ueig2, veig2, weig2, ceig,
     > urealchar, uimagechar, v0char, fname, rechar
     > vrealchar, vimagechar,
     > wrealchar, wimagechar,
     > prealchar, pimagechar,streamfunchar,tmpchar
      common
     >/dimy/yn(0:1024),yn1(0:1024),ym(0:1025),ym1(0:1025),Jm
     >/pry/apym(256),bpym(256),cpym(256),apyn(256),bpyn(256),cpyn(256)
     >/pyp/apyp(256),bpyp(256),cpyp(256),apzn(256), bpzn(256), cpzn(256)
     >/pw/pwork(10000)
     >/alren/al,be0, Re, Nm,v0
     >/basefl/ub(0:257,0:257), vb(0:257,0:257),wb(0:257,0:257),Mm
	     integer info, lwork
     	INTRINSIC INT, MIN

       open(5,file='oz.car')
        read(5,*) Jm
        read(5,*) Mm
        read(5,*) Nm
        read(5,*) al
        read(5,*) be0
        read(5,*) Re
        read(5,*) v0
        close(5)

        be0 = 3.8d0

        open(11,file='c(al,be).dat')


        write(11,300) 50, 50



        do kkk = 1, 12
            be0 = be0 + 0.1
            be = be0
*            al = 0.0d0
            al = 0.d0
        do kk = 1, 50
            al = al+0.1d0

       call Com(eps)
       write(*,*) 'Jm=',Jm,'Mm=', Mm, 'Nm', Nm, 'al=', al ,
     >      'be=', be0 ,'Re=', Re,'v0=', v0
         ci=(0.d0,1.d0)
	     alci=al*ci
	     cial=alci

* check ubase
        Km = 32
        pi = 4*atan(1.)
        hz=2*pi/be/Km
        do k=0,Km+1
         zm(k)=(k-0.5d0)*hz
         zn(k) = k*hz

        end do
          do m = 0 , Mm-1, 2
            do j= 0 , Jm/2+1
                ub(m,j) = -ub(m,Jm+1-j)
            end do
          end do
          do m = 1 , Mm-1, 2
            do j= 0 , Jm/2+1
                ub(m,j) = +ub(m,Jm+1-j)
            end do
          end do
        open(5, file='checkubase.dat')
        write(5,300) Km, Jm
          do j = 1, Jm
            do k = 1, Km
                res = 0.d0
                do m = 0, Mm-1
                    res = res + ub(m,j)*cos(m*be*zm(k))
                end do
                write(5,*) zm(k), ym(j), res
            end do
        end do
        close(5)
*        stop

*


*       Main loop
      Nnm=(Jm-1)*(Nm-1)+(Jm)*(Nm)
*      Nnm = (Jm-1)*(Nm-1)
*        do n=1,Nm
        nn = 0
        do n = 1 , Nm-1
            do j = 1, Jm-1
                nn = nn+1
                do j1 = 0 , Jm+1
                    do n1 = 0 , Nm-1
                        v(n1,j1)=0
                        w(n1,j1)=0
                    end do
                end do
                v(n,j) = (1.d0,0.d0)
                call rp(t,u,v,w,u1,v1,w1,p,Jmax, Nmax)

                m=0
                do n1 = 1, Nm-1
                    do j1 = 1, Jm-1
                        m = m+1
                        amat(m,nn) = v1(n1,j1)
                    end do
                end do

                do n1 = 0 , Nm-1
                    do j1 = 1, Jm
                        m=m+1
                        amat(m,nn) = w1(n1,j1)
                    end do
                end do

            end do
        end do
        write(*,*) nn, '!!!!!!!!!!!!!!!!!!!!!!!!!!'

        do n = 0 , Nm-1
            do j = 1, Jm
                nn=nn+1
                do j1 = 0 , Jm+1
                    do n1 = 0 , Nm-1
                        v(n1,j1)=0
                        w(n1,j1)=0
                    end do
                end do
                w(n,j) = (1.d0,0.d0)

                call rp(t,u,v,w,u1,v1,w1,p,Jmax, Nmax)


                m=0
                do n1 = 1, Nm-1
                    do j1 = 1, Jm-1
                        m = m+1
                        amat(m,nn) = v1(n1,j1)
                    end do
                end do

                do n1 = 0 , Nm-1
                    do j1 = 1, Jm
                        m=m+1
                        amat(m,nn) = w1(n1,j1)
                    end do
                end do

            end do
        end do
        write(*,*) nn
        open(5, file = 'amat.dat')
        do k = 1 , Nnm
            do j = 1 , Nnm
                amat1(j,k)= amat(j,k)
                write(5,*) amat1(j,k), j , k
            end do
        end do
        close(5)

* Eigenvectors and eigenvalues

          info = 2
          lwork=-1
          CALL ZGEEV( 'N', 'V', Nnm, amat, Jmax*Nmax, eval
     >, evecl, Jmax*Nmax, evec, Jmax*Nmax, work, lwork, rwork, info )
          lwork = MIN( LWMAX, INT( work( 1 ) ) )
          CALL ZGEEV( 'N', 'V', Nnm, amat, Jmax*Nmax, eval
     >, evecl, Jmax*Nmax, evec, Jmax*Nmax, work, lwork, rwork, info )
*

      do k=1,Nnm-1
         do n=k+1,Nnm
             if(dreal(eval(n)).gt.dreal(eval(k))) then
              eval(Nnm+1)=eval(k)
              eval(k)=eval(n)
              eval(n)=eval(Nnm+1)
              do j=1,Nnm
                evec(j,Nnm+1)=evec(j,k)
                evec(j,k)=evec(j,n)
                evec(j,n)=evec(j,Nnm+1)
              end do
             end if
         end do
      end do

*
      tmp33 = 0
      do l = 1, Nnm
          do j=1, Nnm
              tmp1 = (0.d0, 0.d0)
              tmp22 = 0
              do  k = 1, Nnm
                  tmp1 = tmp1 + amat1(j,k)*evec(k,l)
              end do
              tmp22 = max(abs(eval(l)*evec(j,l) - tmp1), tmp22)
          end do
              tmp33 = max(tmp22, tmp33)
      end do
      write(*,*) 'zgeev' , tmp33

*

        write(11,*) al,be0, dreal(eval(1))

        write(*,*) al,be0, dreal(eval(1))

* m-я  мода

      meval = 1

      do j = 0, Jm+1
        do n = 0, Nm
            v(n,j) = 0
            w(n,j) = 0
            vt(n,j) = 0
            wt(n,j) = 0
        end do
      end do
      nn = 0
      do n = 1, Nm-1
        do j = 1, Jm-1
            nn = nn+1
            v(n,j) = evec(nn,meval)
        end do
      end do
      do n = 0 , Nm-1
        do j = 1, Jm
            nn=nn+1
            w(n,j) = evec(nn,meval)
        end do
       end do


*        angle = (w(0,1)+w(0,Jm))/abs(w(0,1)+w(0,Jm))
*        angle = (1.d0,0)
*        write(*,*) 'angle', angle, abs(angle)
*        do n = 1, Nm-1
*            do j = 1, Jm-1
*                v(n,j) = v(n,j)/angle
*            end do
*        end do

*        do n = 0, Nm-1
*            do j = 1, Jm
*                w(n,j) =w(n,j)/angle
*                write(*,*) (w(n,j)), abs(w(n,j)), n , j
*            end do
*        end do

        do n = 0, Nm - 1
            v(n,0) = (0.d0,0.d0)
            v(n,Jm) = (0.d0,0.d0)
        end do
        do n = 0, Nm - 1
            w(n,0)    = - w(n,1)
            w(n,Jm+1) = - w(n,Jm)
        end do

        do n = 0, Nm - 1
            do j = 1, Jm
                u(n,j) = -((v(n,j)-v(n,j-1))/ym1(j)
     >          - n*be*w(n,j))/cial
            end do
        end do
        do n = 0, Nm -1
            u(n,0)    = - u(n,1)
            u(n,Jm+1) = - u(n,Jm)
        end do



        write(alchar,10) al
        alchar = adjustl(alchar)
        write(bechar,10) be
        bechar = adjustl(bechar)
        write(v0char,10) v0
        v0char = adjustl(v0char)

      ueig1 = 'ueigen1_al_'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
      veig1 = 'veigen1_al_'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
      weig1 = 'weigen1_al_'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'


      open(4, File = ueig1)
      open(5, File = veig1)
      open(7, File = weig1)


          do j = 1, Jm
            do k = 1, Km
                tmp1 = 0.d0
                do n = 0, Nm-1
                    tmp1 = tmp1 + u(n,j)*sin(n*be*zm(k))
                end do
                write(4,*) zm(k), ym(j), dreal(tmp1),dimag(tmp1)
            end do
          end do




          do j = 1, Jm-1
            do k = 1, Km
                tmp1 = 0.d0
                do n = 0, Nm-1
                    tmp1 = tmp1 + v(n,j)*sin(n*be*zm(k))
                end do
                write(5,*) zm(k), yn(j), dreal(tmp1),dimag(tmp1)
            end do
        end do


          do j = 1, Jm
            do k = 0, Km
                tmp1 = 0.d0
                do n = 0, Nm-1
                    tmp1 = tmp1 + w(n,j)*cos(n*be*zn(k))
                end do
                write(7,*) zn(k), ym(j), dreal(tmp1),dimag(tmp1)
            end do
        end do

      close(4)
      close(5)
      close(7)



        end do
        end do

        close(11)

* m-я  мода

      meval = 1

      do j = 0, Jm+1
        do n = 0, Nm
            v(n,j) = 0
            w(n,j) = 0
            vt(n,j) = 0
            wt(n,j) = 0
        end do
      end do
      nn = 0
      do n = 1, Nm-1
        do j = 1, Jm-1
            nn = nn+1
            v(n,j) = evec(nn,meval)
        end do
      end do
      do n = 0 , Nm-1
        do j = 1, Jm
            nn=nn+1
            w(n,j) = evec(nn,meval)
        end do
       end do


*        angle = (w(0,1)+w(0,Jm))/abs(w(0,1)+w(0,Jm))
        angle = (1.d0,0)
        write(*,*) 'angle', angle, abs(angle)
        do n = 1, Nm-1
            do j = 1, Jm-1
                v(n,j) = v(n,j)/angle
            end do
        end do

        do n = 0, Nm-1
            do j = 1, Jm
                w(n,j) =w(n,j)/angle
                write(*,*) (w(n,j)), abs(w(n,j)), n , j
            end do
        end do

        do n = 0, Nm - 1
            v(n,0) = (0.d0,0.d0)
            v(n,Jm) = (0.d0,0.d0)
        end do
        do n = 0, Nm - 1
            w(n,0)    = - w(n,1)
            w(n,Jm+1) = - w(n,Jm)
        end do

        do n = 0, Nm - 1
            do j = 1, Jm
                u(n,j) = -((v(n,j)-v(n,j-1))/ym1(j)
     >          - n*be*w(n,j))/cial
            end do
        end do
        do n = 0, Nm -1
            u(n,0)    = - u(n,1)
            u(n,Jm+1) = - u(n,Jm)
        end do





        do n = 0, Nm - 1
            do j = 1, Jm
                ut(n,j) = -((vt(n,j)-vt(n,j-1))/ym1(j)
     >          - n*be*wt(n,j))/cial
            end do
        end do
        do n = 0, Nm -1
            ut(n,0)    = - ut(n,1)
            ut(n,Jm+1) = - ut(n,Jm)
        end do

        do j1 = 1, Jm-1
            do n1 = 1, Nm-1
*                    write(*,*) vt(n1,j1), n1, j1
            end do
        end do



        call rp(t,u,v,w,u1,v1,w1,p,Jmax,Nmax)




        tmp11 = 0.d0
        do n = 0 , Nm-1
            do j = 1, Jm
                tmp1 = ut(n,j)*cial
     >            +(vt(n,j)-vt(n,j-1))/ym1(j)
     >          - n*be*(wt(n,j))
                tmp11 = max(abs(tmp1), tmp11)
            end do
        end do
        write(*,*) 'cont2', tmp11

        tmp11 = 0.d0
        do n = 0, Nm-1
            do j = 2, Jm-1
                tmp1 = eval(meval)*u(n,j) - u1(n,j)
                tmp11 = max(tmp11, abs(tmp1))

            end do
        end do
        write(*,*) 'ut', tmp11


        tmp11 = 0.d0
        do n = 1, Nm-1
            do j = 1, Jm-1
                tmp1 = eval(meval)*v(n,j) - v1(n,j)
                tmp11 = max(tmp11, abs(tmp1))

            end do
        end do
        write(*,*) 'vt', tmp11

        tmp11 = 0.d0
        do n = 0, Nm-1
            do j = 1, Jm
*                write(*,*) wt(n,j), n, j
                tmp1 = eval(meval)*w(n,j) - w1(n,j)
                tmp11 = max(tmp11, abs(tmp1))
            end do
        end do
        write(*,*) 'wt',tmp11

        tmp11 = 0.d0
        do n = 0 , Nm-1
            do j = 1, Jm
                tmp1 = u(n,j)*cial
     >            +(v(n,j)-v(n,j-1))/ym1(j)
     >          - n*be*(w(n,j))
                tmp11 = max(abs(tmp1), tmp11)
            end do
        end do
        write(*,*) 'cont', tmp11

*******
        open(5, file='checkueigen.dat')
        write(5,300) Km, Jm
          do j = 1, Jm
            do k = 1, Km
                tmp1 = 0.d0
                do n = 0, Nm-1
                    tmp1 = tmp1 + u(n,j)*sin(n*be*zm(k))
                end do
                write(5,*) zm(k), ym(j), dimag(tmp1)
            end do
        end do
        close(5)


      open(4,file='ureal.dat')
      do n = 0, Nm-1
          do j = 1, Jm
                write(4,*) ym(j), dreal(u(n,j))
          end do
      end do
      close(4)
      open(4,file='uimag.dat')
      do n = 0, Nm-1
          do j = 1, Jm
                write(4,*) ym(j), dimag(u(n,j))
          end do
      end do
      close(4)


      open(4,file='vreal.dat')
      do n = 1, 1
          do j = 0, Jm
                write(4,*) yn(j), dreal(v(n,j))
          end do
      end do
      close(4)
      open(4,file='vimag.dat')
      do n = 1, 1
          do j = 0, Jm
                write(4,*) yn(j), dimag(v(n,j))
          end do
      end do
      close(4)
      open(4,file='wreal.dat')
      do n = 0, 0
          do j = 1, Jm
                write(4,*) ym(j), dreal(w(n,j))
          end do
      end do
      close(4)
      open(4,file='wimag.dat')
      do n = 0, 0
          do j = 1, Jm
                write(4,*) ym(j), dimag(w(n,j))
          end do
      end do
      close(4)

      open(4,file='preal.dat')
      do n = 0, Nm-1
          do j = 1, Jm
                write(4,*) n,j, dreal(p(n,j))
          end do
      end do
      close(4)
      open(4,file='pimag.dat')
      do n = 0, Nm-1
          do j = 1, Jm
                write(4,*) n,j, dimag(p(n,j))
          end do
      end do
      close(4)

*


* выписывание массива собственных значений
        open(6, FILE ='c1.dat')

        do n=1, Nnm
             write(6,*) dreal(eval(n)), dimag(eval(n))
        end do
        close(6)
*

300   format('  ZONE  I=',i4,'  J=',i4,'  F=POINT')
10    format(F10.3)


      end
