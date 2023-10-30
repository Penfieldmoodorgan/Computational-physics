c subrutina que implementa el metode de simpson compost
c retorna integral, que es la aproximacio obtinguda per a la integral emprant aquest metode
c la utilitzarem per a calcular norma

      subroutine simpson(a,b,npas,fcn,integral)
      implicit none
      real*8 h,x,a,b,fcn(500),i,integral
      integer j,npas


      i=fcn(1)+fcn(npas)
      h=(b-a)/dble(npas-1)
      do j=1,npas
        x=a+dble(j)*h
        if (mod(j,2).eq.0) i=i+2.d0*fcn(j)
        if (mod(j,2).ne.0) i=i+4.d0*fcn(j)
      end do  
      integral=i*h/3.d0

      return
      end subroutine