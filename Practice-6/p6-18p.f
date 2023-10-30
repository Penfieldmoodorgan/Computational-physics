c -----------------------------------------------------------------------
C                      PRÀCTICA 6
c-----------------------------------------------------------------------
c      Monte Carlo
C ----------------------------------------------------------------------

c practica 6

      program pre6
      implicit none
      real*8 fun,f1,f2,pi,f3,f4,f5,g,l,psi
      external fun,f1,f2,f3,f4,f5,g,psi
      COMMON/DADES/pi,l
      pi=dacos(-1.0d0)
      l=16.d0 !longitud en micres
     
      call montecarlop6()
      call mcarloMD()


      stop
      end program

c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c              SUBRUTINES
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------

c subrutina metode del rebuig

      subroutine subaccepta(ndat,xnums,a,b,M,fun)
      implicit none
      integer ndat,iseed,i
      real*8 xnums(1:1000000),a,b,M,fun,x1,x2,x,p,f,mitja,mitj2,des,var

      iseed=16404850 !seed niub
            call srand(iseed)

      do i=1,ndat
            f=-1
            p=0
            do while (f.lt.p)
                  x1=rand()
                  x2=rand()
                  x=(b-a)*x1+a
                  p=M*x2
                  f=fun(x)
            end do
            xnums(i)=x
      end do
 
      mitja=0.d0
      var=0.d0

      do i=1,ndat
        mitja=mitja+xnums(i)
      end do
      mitja=mitja/dble(ndat)
      
      do i=1,ndat
            mitj2=(xnums(i)-mitja)**2
            var=var+mitj2
      end do

      des=dsqrt(var)

      open(4,file='P6-18P-res.dat',status='unknown',position='append')
      write(4,*) 'Mostra 20000 nombres funcio:'
      write(4,*) ' '
      write(4,600) 'Valor mitja','Variancia','Desviacio estandard'
      write(4,500) mitja,var,des
      write(4,*) ' '
      write(4,*) ' '
      close(4)

500   format(2x,e20.12,5x,e20.12,5x,e20.12)
600   format(8x,a,14x,a,12x,a)

      return
      end subroutine

c subrutina montecarlo 1D

      subroutine montecarlop6()
      implicit none
      real*8 fun,xa,xb,pi,f1,f2,i1,sigma1
      real*8 i1exacte,x1,f3,f4,f5
      real*8 xnums(1:1000000),m,a,b
      real*8 s1,s2,div,factor1,l
      integer n,j,iseed,ndat
      external f1,f2,fun
      COMMON/DADES/pi,l

c calcul integrals funcio 1 i 2

      i1exacte=(122.d0/27.d0)-(21.d0*(pi**2)/4.d0)!valor exacte integral funcio 1
c      write(*,*) i1exacte,i2exacte
      xa=0.d0 !limit inferior integral 1
      xb=(3.d0/2.d0)*pi !limit superior integral 1

      iseed=16404850 !seed niub
      call srand(iseed)

100   format(10x,a,16x,a,21x,a)
200   format(2x,i10,2(5x,e20.12))


      open(1,file='P6-18P-res.dat',status='unknown')
      write(1,100) 'N','I1','Sigma1'
      write(1,*) ' '
     
c inicializem a zero
      
      s1=0.d0
      s2=0.d0
      n=1000000
c      write(*,*) (xb-xa),(xd-xc)

      do j=1,n 
            x1=xa+(rand()*(xb-xa))
            s1=s1+(f1(x1)*(xb-xa))
            s2=s2+((f1(x1)*(xb-xa))**2)
            i1=s1/dble(j)
            factor1=dsqrt((s2/dble(j))-((s1/dble(j))**2))
            sigma1=factor1/dsqrt(dble(j))
            if ( mod(j,10000).eq.0 ) then
                  write(1,200) j,i1,sigma1
            end if
      end do
      write(1,*) ' '
      write(1,*) ' '
      close(1)

c genera 1000000 numeros segons p(x)=fun(x)

      ndat=1000000

      a=-l/2.d0
      b=l/2.d0
      M=0.15d0
      call subaccepta(ndat,xnums,a,b,M,fun)

c inicializem a zero

      open(2,file='P6-18P-res.dat',status='unknown',position='append')
      write(2,100) 'N','I2','Sigma2'
      write(2,*) ' '

      s1=0.d0
      s2=0.d0
      n=1000000

      do j=1,n
            s1=s1+(f2(xnums(j)))
            s2=s2+(f2(xnums(j))**2)
            i1=s1/dble(j)
            factor1=dsqrt((s2/dble(j))-((s1/dble(j))**2))
            sigma1=factor1/dsqrt(dble(j))
            if ( mod(j,10000).eq.0 ) then
                 write(2,200) j,i1,sigma1
            end if  
      end do
      write(2,*) ' '
      write(2,*) ' '
      close(2)


      return
      end subroutine

      subroutine mcarloMD()
      implicit none
      real*8 pi,l,psi,s1,s2,x1,x2,x3,xa,xb,fun,sigma1,factor1,i1,div
      real*8 xnums(1000000)
      integer i,j,n
      external psi,fun
      COMMON/DADES/pi,l
      COMMON/RANDOM/xnums

      s1=0.d0
      s2=0.d0
      n=250000

      open(8,file='P6-18P-res.dat',status='unknown',position='append')
      write(8,700) 'N','I6','Sigma'
      write(8,*) ' '

      xa=-l/2.d0
      xb=l/2.d0

      i=1           
      do j=1,n
            x1=xa+xnums(i)*(xb-xa)
            x2=xa+xnums(i+1)*(xb-xa)
            x3=xa+xnums(i+2)*(xb-xa)
            i=i+3
            s1=s1+(psi(x1,x2,x3)*(xb-xa))
            s2=s2+((psi(x1,x2,x3)*(xb-xa))**2)
            div=1.d0/dble(j)
            i1=s1*div
            factor1=dsqrt((s2*div)-((s1*div)**2))
            sigma1=factor1*dsqrt(div)
            if ( mod(j,10000).eq.0 ) then
                  write(8,800) j,i1,sigma1
            end if      
      end do
      write(8,*) ' '
      write(8,*) ' '
      close(8)

700   format(10x,a,16x,a,21x,a)
800   format(2x,i10,2(5x,e20.12))


      return 
      end subroutine


c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c                 FUNCIONS      
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------

c------ densitat probabilitat p(x)

      real*8 function fun(x)
      implicit none
      real*8 x,pi,div,temp,l
      COMMON/DADES/pi,l


      div=pi*(x-(l/2.d0))/l
      temp=(2.d0/l)*(dsin(div)**2)
      fun=temp

      return
      end function

c------- funcio 1

      real*8 function f1(x)
      implicit none
      real*8 x,funcio

      funcio=(x**3)*(dsin(x)**3)
      f1=funcio
  
      return
      end function

c ----funcio 2
 
      real*8 function f2(x)
      implicit none
      real*8 x,funcio,div,pi,l
      COMMON/DADES/pi,l

      div=8.d0*pi*(x-(l/2.d0))/l
      funcio=dsin(div)**2
      f2=funcio
  
      return
      end function

c ---funcio 3

      real*8 function psi(x1,x2,x3)
      implicit none
      real*8 pi,l,x(3),x1,x2,x3,prod1,temp1,temp2,prod2,temp3,temp4
      integer i,j
      COMMON/DADES/pi,l

      x=(/x1,x2,x3/)

      do i=1,3
      	temp1=(pi*(x(i)-(l/2.d0)))/l
      	prod1=prod1*dsin(temp1)
      	do j=1,3
      		if (j.lt.i) then
      			temp2=(pi*(x(j)-(l/2.d0)))/l
      			temp3=(pi*(x(i)-(l/2.d0)))/l
      			prod2=prod2*(dcos(temp2)-dcos(temp3))
      		end if
      	end do
      end do

      temp4=dabs(prod1*prod2)
      psi=temp4**2

      return
      end function

