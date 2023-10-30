c Resol el problema T*psi=r 
c T = matriu tridiagonal, A = diagonal inferior, afegim cero primera component
c B = diagonal central, C = diagonal superior, afegim cero a la última component
c els zeros afegits a A i C són per a que tinguin la mateixa dimensió que B
c inputs: vectors A,B,C,R , imax = dimensió vectors
c output: vector PSI

      subroutine tridiag(a,b,c,r,psi,imax)
      implicit none
      real*8 bet,gam(4001)
      real*8 a(imax),b(imax),c(imax),r(imax),psi(imax)
      integer imax,j

      psi=0.d0 !inicialitzem a zero

      if(b(1).eq.0.d0) write(*,*)'Atencio! sistema indeterminat'
      bet=b(1)
      psi(1)=r(1)/bet
      do j=2,imax
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.d0) write(*,*) 'Atencio! bet=0'
        psi(j)=(r(j)-a(j)*psi(j-1))/bet
      end do

      do j=imax-1,1,-1
        psi(j)=psi(j)-gam(j+1)*psi(j+1)
      end do

      return
      end subroutine