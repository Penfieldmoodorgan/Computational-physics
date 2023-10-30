c -----------------------------------------------------------------------
C                      PRÀCTICA 3
c-----------------------------------------------------------------------
c  òrbita cometa Kohoutek. Integració numèrica.
C ----------------------------------------------------------------------
      program pract3

      implicit none
      real*8 ykohoutek,xini,xfin,at,as,h,a,pi,exact,errt,errs,b
      real*8 exact2,xi2,xf2,a2t,a2s,h2,h4
      integer k,n
      external ykohoutek

c      write(*,*) ykohoutek(-3.5d0)
      a=508.633
      b=429.074
      xini=-4.d0*a
      xfin=-3.d0*a
      pi=dacos(-1.0d0)
      exact=pi*a*b

c-----------------FORMATS-----------------------------------------------
c ----------------------------------------------------------------------
100   format(9x,a,22x,a,22x,a,22x,a,17x,a,17x,a,17x,a)
200   format(7(e21.14,3x))
300   format(9x,a,22x,a,22x,a,19x,a,17x,a,17x,a,17x,a)
c ----------------------------------------------------------------------

      open(1,file='P3-18P-res1.dat')
      write(1,100) 'h','At','As','Error t','Error s','h^2','h^4'
      do k=2,20
      	n=2**k
        h=(xfin-xini)/dble(n)
        call trapezis(n,ykohoutek,xini,xfin,at)
        call simpson(xini,xfin,k,ykohoutek,as)
        at=4.d0*at
        as=4.d0*as
        errt=dabs(exact-at)
        errs=dabs(exact-as)
        h2= h**2
        h4=h**4
        write(1,200) h,at,as,errt,errs,h2,h4
      end do
      close(1)

c apartat c

      xi2=-4.d0*a
      xf2=-7.d0*a/2.d0
      exact2=a*b*(3.d0*dsqrt(3.d0)+2.d0*pi)/24.d0

      open(2,file='P3-18P-res2.dat')
      write(2,300) 'h','A2t','A2s','Error t','Error s','h^2','h^4'
      do k=2,20
      	n=2**k
        h=(xf2-xi2)/dble(n)
        call trapezis(n,ykohoutek,xi2,xf2,a2t)
        call simpson(xi2,xf2,k,ykohoutek,a2s)
        errt=dabs(exact2-a2t)
        errs=dabs(exact2-a2s)
        h2=h**2
        h4=h**4
        write(2,200) h,a2t,a2s,errt,errs,h2,h4
      end do
      close(2)

      stop
      end program

c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c                         FUNCIONS
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------

c funcio que calcula orbita del cometa
      real*8 function ykohoutek(x)
      implicit none
      real*8 a,b,funcio,x,div

      a=508.633
      b=429.074

      div=((x+4.d0*a)/a)**2
      funcio=b*dsqrt(1.d0-div)
      ykohoutek=funcio

      return
      end function

c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c                         SUBRUTINES
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------

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
c-----------------------------------------------------------------------