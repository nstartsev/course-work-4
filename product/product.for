      subroutine product()
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
      real*8
     >tmp01, tmp02, tmp1 , tmp2 , tmp31, tmp32 , tmp33
     >, tmp34, tmp35, tmp36 , tmp41 , tmp42, tmp43, tmp44
     >, tmp45 , tmp46, tmp11, tmp12, tmp21, tmp22
     >, res, prod
      character(len = 50)
     > ueig1, veig1, weig1, alchar,bechar, ueig2, veig2, weig2, ceig,
     > urealchar, uimagechar, v0char, fname,  prodchar



         ci=(0.d0,1.d0)
	     alci=al*ci
	     cial=alci

300   format('  ZONE  I=',i4,'  J=',i4,'  F=POINT')

* normalize
      Ekin = 0.d0
      do j = 1, Jm
        do k = 1, Km
          tmp11 = 0.5d0*(vr(j,k)*vr(j,k)+vi(j,k)*vi(j,k))
          tmp12 = 0.5d0*(vr(j-1,k)*vr(j-1,k)+vi(j-1,k)*vi(j-1,k))

          tmp1 = 0.5d0*(tmp11+tmp12)

          tmp21 = 0.5d0*(wr(j,k)*wr(j,k)+wi(j,k)*wi(j,k))
          tmp22 = 0.5d0*(wr(j,k-1)**2 + wi(j,k-1)**2)

          tmp2 = 0.5d0*(tmp21 + tmp22)

          tmp3 = 0.5d0*(ur(j,k)*ur(j,k) + ui(j,k)*ui(j,k))
          Ekin = Ekin + (tmp1+tmp2+tmp3)*ym1(j)*zm1(k)
         end do
       end do


       do j = 0, Jm+1
        do k = 0, Km+1
            ur(j,k) = ur(j,k)/(sqrt(Ekin))
            ui(j,k) = ui(j,k)/(sqrt(Ekin))
        end do
       end do

       do j = 0, Jm
          do k = 0, Km+1
              vr(j,k) = vr(j,k)/(sqrt(Ekin))
              vi(j,k) = vi(j,k)/(sqrt(Ekin))
          end do
       end do

       do j = 0, Jm+1
          do k = 0, Km
              wr(j,k) = wr(j,k)/(sqrt(Ekin))
              wi(j,k) = wi(j,k)/(sqrt(Ekin))
          end do
       end do

        write(*,*) 'Ekin=', Ekin
*

*     product

      prod = 0.d0
10    format(F10.3)
        write(alchar,10) al
        alchar = adjustl(alchar)
        write(bechar,10) be
        bechar = adjustl(bechar)
        write(v0char,10) v0
        v0char = adjustl(v0char)

        prodchar = 'prod'//'_al_'//trim(alchar)//'_v0_'//trim(v0char)//
     > '_be_'//trim(bechar)//'.dat'
        open(3,file = prodchar)
      write(3,300) Km, Jm
      do j = 1, Jm
        do k = 1,Km
          tmp01 = 0.5d0*(vb(j,k)+vb(j-1,k))
          tmp02 = 0.5d0*(wb(j,k)+wb(j,k-1))

          tmp11 = 0.5d0*(vr(j,k)*vr(j,k)+vi(j,k)*vi(j,k))
          tmp12 = 0.5d0*(vr(j-1,k)*vr(j-1,k)+vi(j-1,k)*vi(j-1,k))

          tmp1 = (tmp11-tmp12)/ym1(j)

          tmp21 = 0.5d0*(wr(j,k)*wr(j,k)+wi(j,k)*wi(j,k))
          tmp22 = 0.5d0*(wr(j,k-1)**2 + wi(j,k-1)**2)

          tmp2 = (tmp21 - tmp22)/zm1(k)

          tmp31 = 0.25d0*(vr(j,k)+vr(j,k-1)+vr(j-1,k)+vr(j-1,k-1))
          tmp32 = 0.25d0*(vi(j,k)+vi(j,k-1)+vi(j-1,k)+vi(j-1,k-1))

          tmp33 = 0.25d0*(vr(j,k+1)+vr(j,k)+vr(j-1,k+1)+vr(j-1,k))
          tmp34 = 0.25d0*(vi(j,k+1)+vi(j,k)+vi(j-1,k+1)+vi(j-1,k))

          tmp35 = 0.5d0*(tmp31*wr(j,k-1)+tmp32*wi(j,k-1))
          tmp36 = 0.5d0*(tmp33*wr(j,k)  +tmp34*wi(j,k))

          tmp3 = (tmp36 - tmp35)/zm1(k)

          tmp41 = 0.25d0*(wr(j,k-1)+wr(j,k)+wr(j-1,k-1)+wr(j-1,k))
          tmp42 = 0.25d0*(wi(j,k-1)+wi(j,k-1)+wi(j-1,k-1)+wi(j-1,k))

          tmp43 = 0.25d0*(wr(j+1,k-1)+wr(j+1,k)+wr(j,k-1)+wr(j,k))
          tmp44 = 0.25d0*(wi(j+1,k-1)+wi(j+1,k)+wi(j,k-1)+wi(j,k))

          tmp45 = 0.5d0*(tmp41*vr(j-1,k)+tmp42*vi(j-1,k))
          tmp46 = 0.5d0*(tmp43*vr(j,k)+ tmp44*vi(j,k))

          tmp4 = (tmp46 - tmp45)/ym1(j)

          res = (tmp01*(tmp1+tmp3)+tmp02*(tmp2+tmp4))
          prod = prod + res
          write(3,*) ym(j), zm(k) , -res
*          write(*,*) 'res=',res,'prod=', prod
        end do
      end do
      close(3)
      write(*,*) 'product=', -prod
*     end product
        emain = 0.0d0
*     Emain
      do j = 1, Jm
        do k = 1, Km
          tmp11 = (vb(j,k)-vb(j-1,k))/ym1(j)
          tmp1  = tmp11**2

          tmp21 = (wb(j,k)-wb(j,k-1))/zm1(k)
          tmp2  = tmp21**2

          tmp31 = 0.25d0*(vb(j,k-1)+vb(j,k)+vb(j-1,k-1)+vb(j-1,k))
          tmp32 = 0.25d0*(vb(j,k)+vb(j,k+1)+vb(j-1,k)+vb(j-1,k+1))
          tmp33 = (tmp32- tmp31)/zm1(k)
          tmp3 = tmp33**2

          tmp41 = 0.25d0*(wb(j,k-1)+wb(j,k)+wb(j-1,k-1)+wb(j-1,k))
          tmp42 = 0.25d0*(wb(j+1,k-1)+wb(j+1,k)+wb(j,k-1)+wb(j,k))
          tmp43 = (tmp42-tmp41)/ym1(j)
          tmp4 = tmp43**2

          res = (tmp1+tmp2+tmp3+tmp4)
          Emain = Emain+ res
        end do
      end do
      Emain = Emain/Re
      write(*,*) 'Emain=', Emain
      amplitude = -Emain/prod
      write(20,*) al,be , abs(amplitude)**(1./2)*sign(1.d0,amplitude)
      write(*,*) 'amplitude' , amplitude
*
*

      return
      end
