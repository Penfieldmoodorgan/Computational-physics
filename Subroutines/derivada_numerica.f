c subrutina que crea vectors x i funci
c inputs:ndat: nombre de dades ,fun: funcio ,[a,b] interval on volem avaluar la funcio
c outputs: vectors funci(ndat): valors de la funcio, x(ndat): valors de la variable
      subroutine vectors(ndat,fun,a,b,funci,x)
      implicit none
      real*8 a,b,fun,funci(80),x(80),h,v
      integer ndat,j

      h=(b-a)/dble(ndat-1)
      v=a
      j=1
      do while (v.le.b)
            x(j)=v
            funci(j)=fun(v)
            v=v+h
            j=j+1
      end do

      return
      end subroutine


c subrutina que calcula derivada
c inputs: vector amb valors de la variable equiespaiats (xk+1-xk=h) x(ndat),ndat
c vector amb els valors corresponents a la funcio funci(xk), funci(ndat)
c output: vector amb la derivada calculada numéricament dfunci(ndat)

      subroutine defunT1(ndat,x,funci,dfunci)
      implicit none
      real*8 x(80),funci(80),dfunci(80),h
      integer ndat,k
      
      h=dble(x(2)-x(1))
      dfunci(1)=(funci(2)-funci(1))/h
      do k=2,ndat-1
            dfunci(k)=(funci(k+1)-funci(k-1))/(2.d0*h)
      end do
      dfunci(ndat)=(funci(ndat)-funci(ndat-1))/h

      return
      end subroutine
