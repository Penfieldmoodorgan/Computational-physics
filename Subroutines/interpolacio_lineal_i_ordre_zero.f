c !!!! temps i posis com a common!

C INTERPOLACIO LINEAL x(t)=posis(temps)
c uneix parelles de punts successius amb una línia recta
    
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

C INTERPOLACIO ORDRE ZERO x(t)=posis(temps)
c dóna a la funcio un valor constant, x(tk), dins de cada subinterval [tk,tk+1]
    
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