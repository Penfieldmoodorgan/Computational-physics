c subrutina que implementa el metode dels trapezis
c nint:nombre intervals, (a,b): interval integracio, fcn:funcio a integrar

      subroutine trapezis(ninter,fcn,a,b,integral)
      implicit none
      real*8 h,x,a,b,fcn,integral,int
      integer ninter,i

      int=fcn(a)+fcn(b)
      h=(b-a)/dble(ninter)
      do i=1,ninter-1
        x=a+dble(i)*h
        int=int+2.d0*fcn(x)
      end do
      integral=int*h/2.d0

      return
      end subroutine

c subrutina que implementa el metode de simpson
C (a,b):interval integracio, fcn: funcio a integrar, m:exponent -> nombre intervals: 2^m
c retorna integral, que es la aproximacio obtinguda per a la integral emprant aquest metode

      subroutine simpson(a,b,m,fcn,integral)
      implicit none
      real*8 h,x,a,b,fcn,i,integral
      integer j,m,n

      i=fcn(a)+fcn(b)
      n=2**m
      h=(b-a)/dble(n)
      do j=1,n-1
        x=a+dble(j)*h
        if (mod(j,2).eq.0) i=i+2.d0*fcn(x)
        if (mod(j,2).ne.0) i=i+4.d0*fcn(x)
      end do  
      integral=i*h/3.d0

      return
      end subroutine

