C CALCULA UNA INTEGRAL EMPRANT EL METODE DE ROMBERG

      PROGRAM ROMBERG

      IMPLICIT NONE
      REAL*8 f,pi,a,b,exact,Int,factor,precisio,error,x1,x2,e
      REAL*8 T(1:20,1:20)
      INTEGER i,j,k,w

C DEFINEIX ELS LIMITS DE INTEGRACIO

      pi=dacos(-1.0d0)
      a=pi/4.d0
      b=0.75d0*pi

C DEFINEIX EL VALOR EXACTE DE LA INTEGRAL CALCULAT NUMERICAMENT   

      exact=dexp(0.75*pi)/dsqrt(2.0d0)
    
C MENTRE LA DIFERENCIA ENTRE ELS DOS ULTIMS VALORS DE LA DIAGONAL DE LA MATRIU SIGUI MAJOR QUE LA PRECISIO DEFINIDA, APLICA EL ALGORISME I GUARDA ELS VALORS A LA MATRIU
      precisio=1.d-15
      i=1
      e=1

      DO WHILE (e.GT.precisio)
        CALL trapezis(2**(i-1),Int,a,b)
        T(i,1)=Int
        DO w=i+1,20
          T(i,w)=0
        END DO
        DO j=2,i
          factor=1.d0/((4**j)-1.d0)
          T(i,j)=((4**j)*T(i,j-1)-T(i-1,j-1))*factor
          x1=T(i-1,i-1)
          x2=T(i,i)
          e=error(x1,x2)
        END DO
        i=i+1
      END DO
      k=i-1

      OPEN(1,FILE='Romberg.dat')
    
C ESCRIU ELS VALORS OBTINGUTS PER A CADA ITERACIO
C COM FORTRAN ESCRIU LES DADES PER COLUMNES, PER A OBTENIR LA VISUALITZACIÓ DE LES DADES DESITJADA AL FITXER, ESCRIU LA MATRIU TRANSPOSADA. 

      DO i=1,k
        WRITE(1,*) (T(j,i),j=1,k)
      END DO

C ESCRIU AL FITXER LA APROXIMACIO OBTINGUDA, EL VALOR EXACTE DE LA INTEGRAL DEFINIDA I EL ERROR DE LA PRIMERA. 
      WRITE(1,*) ''
      WRITE(1,*) 'Nombre iteracions', k
      WRITE(1,*) 'Valor obtingut:', T(k,k)
      WRITE(1,*) 'Valor exacte:', exact
      WRITE(1,*) 'Error aproximacio:', error(T(k,k),exact)


      CLOSE(1)
      STOP 
      END PROGRAM

C SUBRUTINA QUE IMPLEMENTA EL METODE DELS TRAPEZIS

      SUBROUTINE TRAPEZIS(n,Int,a,b) 
      IMPLICIT NONE
      REAL*8 Int,h,x,a,b,f
      INTEGER n,i 

      x=a
      Int=f(a)+f(b)
      h=(b-a)/dble(n)
      DO i=1,n-1
        x=x+h
        Int=Int+2.d0*f(x)
      END DO
      Int=Int*h/2.d0

      RETURN
      END SUBROUTINE

C FUNCIO PER A LA QUE ES VOL CALCULAR LA INTEGRAL

      REAL*8 FUNCTION f(x)
      IMPLICIT NONE
      REAL*8 x,temp

      temp=dexp(x)*dsin(x)
      f=temp

      RETURN
      END FUNCTION

C FUNCIO QUE CALCULA EL ERROR RELATIU

      REAL*8 FUNCTION error(x1,x2)
      IMPLICIT NONE
      REAL*8 x1,x2,temp1

      temp1=(x2-x1)/x1
      error=abs(temp1)

      RETURN
      END FUNCTION



      
    