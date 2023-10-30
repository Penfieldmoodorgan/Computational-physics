c -----------------------------------------------------------------------
C                      PRÀCTICA 1
c-----------------------------------------------------------------------
c Sumatori de funció Pk
C ----------------------------------------------------------------------

      PROGRAM PRACTICA

      IMPLICIT NONE
      REAL*8 div,j,p,suma,s,asimpt
      INTEGER n
      
C LLEGEIX NUMERO ENTER K I CALCULA PK, MOSTRA VALOR A PANTALLA
C SI EL NOMBRE NO ES TROBA ENTRE 3 I 35 TORNA A PREGUNTAR
C SI EL NOMBRE NO ES UN ENTER TORNA A PREGUNTAR

      j = 0

      DO WHILE ((j.LT.15).OR.(j.GT.221).OR.(INT(j).NE.j))
            WRITE (*,*) 'Introdueixi un nombre enter entre 15 i 221'
            READ  (*,*) j
            IF ((j.LT.15).OR.(j.GT.221)) THEN
                  WRITE(*,*)'Aquest numero no es troba al interval.'
            ELSE IF (INT(j).NE.j) THEN
                  WRITE(*,*) 'Aquest no es un numero enter'
            ENDIF
      END DO
      WRITE(*,*) 'Valor Pk', P(INT(j))

C CALCULA SUMATORI DES DE 28 FINS A 65

      WRITE (*,*) 'Valor del sumatori des de 28 fins a 65', SUMA(28,65)

C CREA FITXER (1) AMB TRES COLUMNES
C COLUMNA 1: VALORS N ENTRE 11 I 311 DE 3 EN 3 
C COLUMNA 2 : VALOR SUMATORI DES DE N FINS A 8
C COLUMNA 3 : COMPORTAMENT ASIMPTOTIC 1/5*N^3

C CREA FITXER (2) AMB DOS COLUMNES: N, DIVISIO ENTRE EL SUMATORI DES DE N FINS A 8 I COMPORTAMENT ASIMPTOTIC
      
      OPEN(1,FILE='P1-18P-res1.dat')
      OPEN(2,FILE='P1-18P-res2.dat')

      DO n = 11,311,3
            s = SUMA(8,n)
            asimpt = (1.0d0/5.0d0)*(n**3)
            div = s/asimpt
            WRITE (1,*) n,s,asimpt 
            WRITE (2,*) n,div     
      END DO
      CLOSE(1)
      CLOSE(2)

      STOP
      END PROGRAM 
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c                         FUNCIONS
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------

C FUNCIO CALCUL Pk  

      REAL*8 FUNCTION p(k)
      IMPLICIT NONE
      INTEGER k
      REAL*8 e

      e = dexp(1.0d0)

      p = (3.0d0/5.0d0)*(k**2)+e+(10.0d0*k)
      RETURN
      END FUNCTION

C FUNCIO CALCUL SUMATORI

      REAL*8 FUNCTION SUMA(n1,n2)
      IMPLICIT NONE
      INTEGER n1,n2,k
      REAL*8 p,temp

      temp = 0.0d0
      DO k = n1,n2
            temp = temp + p(k)
      END DO
      suma=temp
      RETURN
      END FUNCTION      
 
C-----------------------------------------------------------------------
