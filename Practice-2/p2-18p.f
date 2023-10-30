c -----------------------------------------------------------------------
C                      PRÀCTICA 2
c-----------------------------------------------------------------------
c posició pistons en funció del temps. Interpolació
C ----------------------------------------------------------------------
      PROGRAM PISTO
      IMPLICIT NONE
      REAL*8 radit1,w0,L,phi,t,r,m,n,h,augment,tin,xout,x1
      REAL*8 x(5),posis(0:500),temps(0:500)
      INTEGER j
      COMMON/DADES/temps,posis

c-----------------FORMATS-----------------------------------------------
c ----------------------------------------------------------------------
100   FORMAT(6(E20.14,4X))
200   FORMAT(3(E20.14,4X))
c-----------------------------------------------------------------------

      w0=5.d0
      L=18.5d0
      
C GUARDA EN UN FITXER EL TEMPS, I LA POSICIÓ DELS QUATRE PISTONS EN FUNCIO D'AQUEST

      OPEN(1,file='P2-18P-res1.dat')

      t=0.d0

      DO WHILE (t.LE.5)
        CALL posit1(w0,L,t,x)
        WRITE(1,100) t,x(1:5)
        t=t+0.01d0
      END DO
      CLOSE(1)

  

      OPEN(2,file='P2-18P-res1.dat',status='old')
c      OPEN(4,file='miau.dat')
      DO j=0,500
            READ(2,*)temps(j),r,m,n,posis(j),h
c            WRITE(4,*)temps(j),posis(j)
      END DO 
      CLOSE(2) 
c      CLOSE(4)

C GUARDA EN UN FITXER ELS VALORS DE TEMPS I POSICIO DEL PISTO 4 OBTINGUTS AMB LA INTERPOLACIO

      OPEN(3,file='P2-18P-res2.dat')
      augment=(3.d0/2000.d0)
      tin=0.d0
      DO WHILE(tin.LE.3)
        CALL xinterpol(tin,xout)
        x1=xout
        CALL xinterpol0(tin,xout)
        WRITE(3,200) tin,xout,x1
        tin=tin+augment
      END DO
      CLOSE(3)

      STOP
      END PROGRAM

c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c                         FUNCIONS
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------

C VARIABLES ENTRADA: LONGITUD I NUMERO DEL PISTO
C RETORNA RADI PISTO K EN CM       

      REAL*8 FUNCTION radit1(L,k)
      IMPLICIT NONE
      REAL*8 L,temp
      INTEGER k

      temp=(L/dble(k))-0.5d0
      radit1=temp

      RETURN
      END FUNCTION

C FUNCIO QUE CALCULA LA FASE INICIAL
      REAL*8 FUNCTION phi(k)
      IMPLICIT NONE
      REAL*8 pi,temp1
      INTEGER k

      pi=dacos(-1.0d0)
      temp1=pi*(dble(k)/5.d0)**2
      phi=temp1

      RETURN
      END FUNCTION

c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c                         SUBRUTINES
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------

      SUBROUTINE posit1(w0,l,t,x)
      IMPLICIT NONE
      REAL*8 x(5),phi,L,t,w0,radit1,s1,s2
      INTEGER k

      DO k=1,5
        s1=radit1(L,k)*dcos(w0*t+phi(k))
        s2=dsqrt(L**2-(radit1(L,k)**2)*dsin(w0*t+phi(k))**2)
        x(k)=s1+s2
      END DO

      RETURN
      END SUBROUTINE

C INTERPOLACIO LINEAL
    
      SUBROUTINE xinterpol(tin,xout)
      IMPLICIT NONE
      REAL*8 tin,xout,posis(0:500),temps(0:500),xt,div
      INTEGER a0,a1
      COMMON/DADES/temps,posis

      a0=int(tin/0.01d0)
      a1=a0+1

      div=(posis(a1)-posis(a0))*(tin-temps(a0))/(temps(a1)-temps(a0))
      xt=posis(a0)+div

      xout=xt

      RETURN
      END SUBROUTINE     

C INTERPOLACIO ORDRE ZERO
    
      SUBROUTINE xinterpol0(tin,xout)
      IMPLICIT NONE
      REAL*8 tin,xout,posis(0:500),temp2,temps(0:500),t0,t1
      INTEGER i
      COMMON/DADES/temps,posis

      DO i=0,499
        t0=temps(i)
        t1=temps(i+1)
        IF((t0.LE.tin).AND.(t1.GE.tin)) temp2=posis(i)
      END DO

      xout=temp2

      RETURN
      END SUBROUTINE     

C ----------------------------------------------------------------------