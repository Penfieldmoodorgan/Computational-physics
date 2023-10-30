c subrutina metode implicit
c inputs kappa, vector per a t=0 (tf), numero de punts per x (xpunts), numero de punts per t (npas)
c output matriu psi1(xpunts,npas)

      subroutine implicit(k,tf,xpunts,npas,psi1)
      implicit none
      integer xpunts,npas,i,j
      real*8 k,tf(0:xpunts),psi1(npas,xpunts-1),dt,h
      real*8 alpha,a(xpunts-1),b(xpunts-1),c(xpunts-1),t
      real*8 r(xpunts-1),psi(xpunts-1)
      common/dades/h,dt

c diagonals
      alpha=k*dt/(h**2)      
      a=-alpha
      b=1.d0+2.d0*alpha
      c=-alpha
      
c posem zeros per a igualar dimensions
      a(1)=0.d0
      c(xpunts-1)=0.d0

      do i=1,xpunts-1
        r(i)=tf(i)
      end do

      r(1)=r(1)+alpha*tf(0)
      r(xpunts-1)=r(xpunts-1)+alpha*tf(xpunts)

c càlcul

      do i=1,npas
        t=0.d0+i*dt
        call tridiag(a,b,c,r,psi,xpunts-1)
        do j=1,xpunts-1
          psi1(i,j)=psi(j)
        end do
        r=psi
        r(1)=r(1)+alpha*tf(0)
        r(xpunts-1)=r(xpunts-1)+alpha*tf(xpunts)
      end do

      return
      end subroutine