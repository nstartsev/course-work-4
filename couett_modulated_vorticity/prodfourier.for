      subroutine prodcos(a,b, c, Nm, Mm, Lm,Nmax)
      implicit real*8 (a-h,o-z)
	     complex*16
     >  a(0:Nmax),b(0:Nmax), c(0:Nmax) ,tmp, tmp1, tmp2, tmp3
       do n = Nm, Nm+Mm
            a(n) = (0.d0,0.d0)
       end do
        do m = Mm, Nm+Mm
            b(m) = (0.d0,0.d0)
        end do

* c(0)
        c(0) = 1.d0/2*a(0)*b(0)
        tmp = (0.d0,0.d0)
        do n = 0, Nm-1
            tmp = tmp + 1.d0/2*a(n)*b(n)
        end do
        c(0) = c(0) + tmp

        do l = 1, Lm-1
            tmp1 = (0.d0,0.d0)
            tmp2 = (0.d0, 0.d0)
            tmp3 = (0.d0, 0.d0)
            do n = 0, l
                tmp1 = 1.d0/2*a(n)*b(l-n) + tmp1
            end do
            do n = 0, Nm-1
                tmp2 = 1.d0/2*a(n)*b(l+n) + tmp2
            end do
            do m = 0, Mm-1
                tmp3 = 1.d0/2*a(l+m)*b(m) + tmp3
            end do
            c(l) = tmp1 + tmp2 + tmp3
        end do

        end

      subroutine prodsin(a,b, c, Nm, Mm, Lm,Nmax)
      implicit real*8 (a-h,o-z)
	     complex*16
     >  a(0:Nmax),b(0:Nmax), c(0:Nmax) , tmp, tmp1, tmp2, tmp3
       do n = Nm, Nm+Mm
            a(n) = (0.d0,0.d0)
       end do
        do m = Mm, Nm+Mm
            b(m) = (0.d0,0.d0)
        end do

* c(0)
        c(0) = 1.d0/2*a(0)*b(0)
        tmp = (0.d0,0.d0)
        do n = 0, Nm-1
            tmp = tmp + 1.d0/2*a(n)*b(n)
        end do
        c(0) = c(0) + tmp

        do l = 1, Lm-1
            tmp1 = (0.d0,0.d0)
            tmp2 = (0.d0, 0.d0)
            tmp3 = (0.d0, 0.d0)
            do n = 0, l
                tmp1 = 1.d0/2*a(n)*b(l-n) + tmp1
            end do
            do n = 0, Nm-1
                tmp2 = 1.d0/2*a(n)*b(l+n) + tmp2
            end do
            do m = 0, Mm-1
                tmp3 = 1.d0/2*a(l+m)*b(m) + tmp3
            end do
            c(l) = -tmp1 + tmp2 + tmp3
        end do

        end
      subroutine prodsincos(a,b, c, Nm, Mm, Lm,Nmax)
      implicit real*8 (a-h,o-z)
	     complex*16
     >  a(0:Nmax),b(0:Nmax), c(0:Nmax) , tmp ,tmp1, tmp2, tmp3
       do n = Nm, Nm+Mm
            a(n) = (0.d0,0.d0)
       end do
        do m = Mm, Nm+Mm
            b(m) = (0.d0,0.d0)
        end do

* c(0)
        c(0) = a(0)*b(0)


        do l = 1, Lm-1
            tmp1 = (0.d0,0.d0)
            tmp2 = (0.d0, 0.d0)
            tmp3 = (0.d0, 0.d0)
            do n = 0, l
                tmp1 = 1.d0/2*a(n)*b(l-n) + tmp1
            end do
            do n = 0, Nm-1
                tmp2 = 1.d0/2*a(n)*b(l+n) + tmp2
            end do
            do m = 0, Mm-1
                tmp3 = 1.d0/2*a(l+m)*b(m) + tmp3
            end do
            c(l) = tmp1 - tmp2 + tmp3
        end do

        end
