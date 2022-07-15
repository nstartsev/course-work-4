*
      subroutine rp(t,u,v,w,ut,vt,wt,p,Jmax, Nmax)
      implicit real*8 (a-h,o-z)
	     complex*16
     > u(0:Nmax,0:Jmax),v(0:Nmax,0:Jmax),w(0:Nmax,0:Jmax)
     >,ut(0:Nmax,0:Jmax),vt(0:Nmax,0:Jmax),wt(0:Nmax,0:Jmax)
     >,p(0:Nmax,0:Jmax),dp(0:Nmax,0:Jmax)
     >,ci,alci,cial, tmp , tmp1, tmp2 , tmp3,tmp4,tmp5,tmp6
     >, ox(0:Nmax,0:Jmax), oy(0:Nmax,0:Jmax), oz(0:Nmax,0:Jmax),
     > a(0:Nmax),b(0:Nmax), c1(0:Nmax), c2(0:Nmax), c3(0:Nmax),
     > c4(0:Nmax)
     >, oxb(0:Nmax,0:Jmax), oyb(0:Nmax,0:Jmax), ozb(0:Nmax,0:Jmax)
         real*8
     >  dpr(Nmax,Jmax), dpi(Nmax,Jmax)

      common
     >/dimy/yn(0:1024),yn1(0:1024),ym(0:1025),ym1(0:1025),Jm
     >/pry/apym(256),bpym(256),cpym(256),apyn(256),bpyn(256),cpyn(256)
     >/pyp/apyp(256),bpyp(256),cpyp(256),apzn(256), bpzn(256), cpzn(256)
     >/pw/pwork(10000)
     >/alren/al,be, Re, Nm
     >/basefl/ub(0:257,0:257), vb(0:257,0:257),wb(0:257,0:257),Mm


      ci=(0.d0,1.d0)
      alci=al*ci
      cial=alci
*velocity
        do n = 0, Nm - 1
            v(n,0) = (0.d0,0.d0)
            v(n,Jm) = (0.d0,0.d0)
        end do
        do n = 0, Nm - 1
            w(n,0)    = - w(n,1)
            w(n,Jm+1) = - w(n,Jm)
        end do
        diff = 0.
        do n = 0, Nm - 1
            do j = 1, Jm
                u(n,j) = -((v(n,j)-v(n,j-1))/ym1(j)
     >          - n*be*w(n,j))/cial
                diff1 = abs(u(n,j))
                diff = max(diff1, diff)
            end do
        end do
*        write(*,*) diff
        do n = 0, Nm -1
            u(n,0)    = - u(n,1)
            u(n,Jm+1) = - u(n,Jm)
        end do



*
***

* Vorticities

      do n = 0, Nm-1
            do j = 0, Jm
                oxb(n,j) = 0.
                oyb(n,j) = -(n*be)*ub(n,j)
                ozb(n,j) = -(ub(n,j+1)-ub(n,j))/yn1(j)
            end do
      end do
        diff = 0.
      do n = 0, Nm-1
            do j = 0, Jm
                ox(n,j) = (w(n,j+1)-w(n,j))/yn1(j)
     >                       - (n*be)*v(n,j)
                oy(n,j) =  (n*be)*u(n,j)
     >                    -  cial*w(n,j)
                oz(n,j) = cial*v(n,j)
     >          - (u(n,j+1)-u(n,j))/yn1(j)
*                write(*,*) oz(n,j), n, j
                diff1 = abs(ox(n,j))
                diff = max(diff, diff1)
            end do
      end do
*      write(*,*) diff
****


*   ut(n,j)
        do j = 1, Jm




            do n = 0, Nm - 1

                a(n) = v(n,j)
                b(n) = ozb(n,j)
            end do
            call prodsincos(a,b,c2,Nm,Nm, Nm, Nmax)

            do n = 0, Nm - 1

                a(n) = v(n,j-1)
                b(n) = ozb(n,j-1)
            end do
            call prodsincos(a,b,c3,Nm, Nm, Nm, Nmax)

            do n = 0, Nm - 1

                a(n) = oyb(n,j)
                b(n) = w(n,j)

            end do
            call prodsincos(a,b,c4,Nm, Nm, Nm, Nmax)


            do n = 0 , Nm-1
                tmp1 = 0.
                tmp2 = 0.5d0 *( c2(n)+c3(n))
                tmp3 = - c4(n)
                tmp5 = - (oz(n,j)-oz(n,j-1))/(ym1(j)*Re)
                tmp6 =  - (n*be)*oy(n,j)/Re
                ut(n,j) =  tmp2 + tmp3  + tmp5 + tmp6

            end do
      open(4,file='dp.dat')
      do n = 0, 0
*          do j = 1, Jm
*                write(*,*) n,j, (ut(n,j))
*          end do
      end do
      close(4)

        end do

*        pause
**


*   vt(n,j)
        do j = 0, Jm
            do n = 0, Nm - 1

                a(n) = oz(n,j)
                b(n) = 0.5d0*(ub(n,j+1)+ub(n,j))

            end do
            call prodsincos(a,b,c1,Nm,Nm, Nm, Nmax)


            do n = 0, Nm - 1

                a(n) = 0.5d0*(u(n,j+1)+u(n,j))
                b(n) = ozb(n,j)

            end do
            call prodsincos(a,b,c2,Nm,Nm, Nm, Nmax)



            do n = 0 , Nm-1
            tmp3 = -c1(n)
            tmp4 = -c2(n)
            tmp5 = (n*be)*ox(n,j)/Re
            tmp6 = +cial*oz(n,j)/Re
            vt(n,j) = tmp3 + tmp4 + tmp5 + tmp6
            end do
        end do
****



*   wt(n,j)
        do j = 1, Jm
            do n = 0, Nm-1

                b(n) = ub(n,j)
                a(n) = oy(n,j)
            end do

            call prodcos(a,b,c1,Nm,Nm,Nm, Nmax)

            do n = 0, Nm-1

                a(n) = u(n,j)
                b(n) = oyb(n,j)

            end do

            call prodsin(a,b,c2,Nm,Nm,Nm, Nmax)
            do n = 0, Nm-1
            tmp1 = c1(n)
            tmp2 = c2(n)
            tmp5 = -cial* oy(n,j)/Re
            tmp6 = (ox(n,j)-ox(n,j-1))/(ym1(j)*Re)
            wt(n,j) = tmp1 + tmp2 + tmp5 + tmp6
            end do
        end do
*

*





* pressure
        do n = 0, Nm - 1
            do j = 1, Jm

                tmp1 =   cial*ut(n,j)
                tmp2 = (vt(n,j)-vt(n,j-1))/ym1(j)
                tmp3 = -(n*be)*wt(n,j)
                dp(n,j) = (tmp1+tmp2+tmp3)


            end do
        end do

        do n = 0, Nm-1
            dp(n,1)  = dp(n,1)  + apym(1) *yn1(0)*vt(n,0)
            dp(n,Jm) = dp(n,Jm) - cpym(Jm)*yn1(Jm)*vt(n,Jm)
        end do

        do j = 1 , Jm
            do n = 1 , Nm
                dpr(n,j) = dreal(dp(n-1,j))
                dpi(n,j) = dimag(dp(n-1,j))
            end do
        end do


        call blktri(1,1,Nm,apzn,bpzn,cpzn,1,Jm,apyp,bpyp,cpyp,Nmax,
     >               dpr,ierror,pwork)


        call blktri(1,1,Nm,apzn,bpzn,cpzn,1,Jm,apyp,bpyp,cpyp,Nmax,
     >               dpi,ierror,pwork)


        do j = 1, Jm
            do n = 0, Nm-1
                p(n,j) = dpr(n+1,j) + ci*dpi(n+1,j)
*                write(*,*) n-1, j, p(n,j)
            end do
        end do

*

* check pressure
        do n = 0, Nm-1
            p(n,0) = p(n,1) - yn1(0)*vt(n,0)
            p(n,Jm+1) = p(n,Jm) + yn1(Jm)*vt(n,Jm)
        end do

        diff = 0.d0
        do j = 1, Jm
            do n = 0, Nm-1
                tmp1 = (cial)**2*p(n,j) - (n*be)**2*p(n,j) +
     >          (p(n,j-1)-2*p(n,j)+p(n,j+1))/(ym1(j)**2)
                tmp2 = ut(n,j)*cial +(vt(n,j)-vt(n,j-1))/ym1(j)
     >          - n*be*wt(n,j)
                diff1 = abs(tmp1-tmp2)
                diff = max(diff,diff1)
*                write(*,*) n,j, p(n,j)
            end do
        end do
*        write(*,*) diff
*        stop
        diff = 0.
        do n = 0, Nm - 1
            do j = 1, Jm
                tmp1 = (cial)*(cial*p(n,j)-ut(n,j)) +
     >          n*be*(-(n*be)*p(n,j)+wt(n,j))
                tmp2 = (p(n,j-1)-2*p(n,j)+p(n,j+1))/(ym1(j)**2) -
     >           (vt(n,j)-vt(n,j-1))/ym1(j)
                tmp3 = (p(n,j+1)-p(n,j))/yn1(j) - vt(n,j)
                tmp4 = (p(n,j)-p(n,j-1))/yn1(j-1) - vt(n,j-1)
                tmp5 = (tmp3-tmp4)/ym1(j)

                diff1 = abs(tmp1+tmp5)
                diff = max(diff,diff1)
            end do
        end do
*        write(*,*) diff
*        do j = 1 , Jm
*            p(0,j) = 0
*        end do

        do n = 0, Nm-1
            do j = 1, Jm
                ut(n,j) = ut(n,j) - cial*p(n,j)
            end do
        end do

        do n = 0, Nm -1
            ut(n,0)    = - ut(n,1)
            ut(n,Jm+1) = - ut(n,Jm)
        end do
        diff = 0.
        do n = 0, Nm-1
            do j = 0, Jm
                vt(n,j) = vt(n,j) - (p(n,j+1)-p(n,j))/yn1(j)
                diff1 = abs(vt(n,j))
                diff = max(diff,diff1)
            end do
        end do
*        write(*,*) diff

        do n = 0, Nm-1
            do j = 0, Jm
                wt(n,j) = wt(n,j) - (n*be) * p(n,j)
            end do
        end do
        do n = 0, Nm - 1
            wt(n,0)    = - wt(n,1)
            wt(n,Jm+1) = - wt(n,Jm)
        end do

      return
      end
