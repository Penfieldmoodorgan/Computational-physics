c -----------------------------------------------------------------------
C                      PRÀCTICA 5
c-----------------------------------------------------------------------
c      moviment brownià
C ----------------------------------------------------------------------

      program pre5
      implicit none
      real*8 xhisto(1:120),errhisto(1:120),histo(1:120),tf,d(1:250)
      real*8 xgaus(1:120000),mu,sigma,xa,xb,delta,x(1:250),y(1:250),t,z
      real*8 var,inx,iny,suma,suma2,mitj2,variancia(1:240),mitj,fun,mes
      integer ndat,ncaixes,i,j,nm,contador
      external fun

c-----------------FORMATS-----------------------------------------------
c ----------------------------------------------------------------------
100   format(12x,11(a,23x))
200   format(2x,11(e20.12,5x))
300   format(12x,a,20x,a)
400   format(2x,e20.12,5x,e20.12)
c-----------------------------------------------------------------------

c apartat 1

      ndat=120000
      mu=0.d0
      sigma=1.d0
      call subgaussians(ndat,mu,sigma,xgaus)

      ncaixes=120
      xa=-5.d0
      xb=5.d0
      call hist(ndat,xgaus,xa,xb,ncaixes,xhisto,histo,errhisto)

      open(6,file='P5-18P-res1.dat',status='unknown',position='append')
      mes=(xb-xa)/150.d0
      z=0.d0
      do while(z.lt.xb)
      	z=z+mes
      	write(6,400) z, fun(z)
      end do
      close(6)


c apartat 2

      delta=0.02d0 !increment de temps
      nm=250 ! molecules oxigen independents en 2D
      var=delta*2.21d-5 ! variancia distribucio increment x, increment y
      x=0.d0 !x inicial (1,...,nm)
      y=0.d0 !y inicial (1,...,nm)

      open(1,file='P5-18P-res2.dat',status='unknown')
      write(1,100) 't','x1','y1','x2','y2','x3','y3','x4','y4','x5','y5'

      contador = 1
      do i=1,240 !240 passos de temps
      	t=i*delta
        suma=0.d0
        suma2=0.d0
      	do j=1,nm !nm components de x i y
      		inx=dsqrt(var)*xgaus(contador)
      		iny=dsqrt(var)*xgaus(contador+1)
      		x(j)=x(j)+inx
      		y(j)=y(j)+iny
      		suma=suma+y(j)
      		suma2=suma2+y(j)**2
      		contador=contador+2
      	end do
      	write(1,200) t,x(1),y(1),x(2),y(2),x(3),y(3),x(4),y(4),x(5),y(5)
      	mitj=suma/dble(nm)
      	mitj2=suma2/dble(nm)
      	variancia(i)=mitj2-mitj**2
      end do
      write(1,*) ' '
      write(1,*) ' '
      close(1)

c apartat b

      open(2,file='P5-18P-res2.dat',status='unknown',position='append')
      write(2,300) 't','Var(y(t))'
      do i=1,240
      	t=i*delta
      	write(2,400) t,variancia(i)
      end do
      write(2,*) ' '
      write(2,*) ' '
      close(2)

c apartat extra
      
c      open(3,file='P5-18P-extra.dat',status='unknown')

      do j=1,nm
        d(j)=dsqrt((x(j)**2)+(y(j)**2))
c        write(3,*) d(j)
      end do
c      close(3)
      
c FALTA TROS!
      stop
      end program


c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c              SUBRUTINES
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------

c subrutina que crea histograma normalitzat
c inputs: ndat=nombre de dades, xdata=vector amb ndat components dels valor de x,
c [xa,xb]=limits superior i inferior interval, ncaixes=nombre de caixes del histograma
c outputs: xhisto=vector amb ncaixes components que indiquen posicio central caixa
c histo=vector amb ncaixes components que son el nombre de valors de x en cada caixa entre el nombre de valors total
c errhisto=vector amb ncaixes components que son el error de cada caixa

      subroutine hist(ndat,xdata,xa,xb,ncaixes,xhisto,histo,errhisto)
      implicit none
      integer ncaixes,ndat,i,j
      integer inthisto(1:150)
      real*8 xdata(1:120000),xhisto(1:120),errhisto(1:120),histo(1:120)
      real*8 xa,xb,l,valor,linf,lsup,div,sigma,suma
       
      inthisto=0  ! inicialitza totes components histo a zero
      l=dabs(xb-xa)/(ncaixes-1)

      do j=1,ndat
      	valor=xdata(j)
        if ((valor.le.xb).and.(valor.ge.xa)) then !cond valor dins interval
        	do i=1,ncaixes
        		xhisto(i)=xa+(i-1)*l
                linf=xhisto(i)-(l/2.d0)
                lsup=xhisto(i)+(l/2.d0)
                if((valor.gt.linf).and.(valor.le.lsup)) then ! si valor dins de caixa compta +1
                	inthisto(i)=inthisto(i)+1
                end if
             end do
        end if
      end do

c normalitzem nk=nk/wkN 
      
      suma=0.d0
      do i=1,ncaixes
      	div=dble(inthisto(i))/(dble(ndat)*l) ! posar amplada inteval dividint o no?
        histo(i)=div                         ! sense amplada queda suma total =1
        sigma=dsqrt(div*(1.d0-div)/dble(ndat))/l  
c        errhisto(i)=3.d0*sigma ! posar error +-3sigma?o nomes sigma?
        errhisto(i)=sigma 
c        suma=suma+histo(i)
      end do
c      write(*,*)suma

c escriu valors en fitxer

500   format(2x,e20.12,5x,e20.12,5x,e20.12)
600   format(8x,a,16x,a,16x,a)

      open(3,file='P5-18P-res1.dat',status='unknown')
      write(3,*) 'Dades histograma normalitzat:'
c      write(3,*) '(Error de tres desviacions estandard, barres error con
c     +tindran el valor real en un 99.7% dels casos) '
      write(3,*) ' '
      write(3,600) 'Punt mig', 'Numero valors','Error'
c      write(3,*) ' '
      do i=1,ncaixes
      	write(3,500) xhisto(i),histo(i),errhisto(i)
      end do
      write(3,*) ' '
      write(3,*) ' '
      close(3)   

      return
      end subroutine

c subrutina que genera ndat nombres gaussians de valor mitja xmu i variancia igual a xsigma^2

      subroutine subgaussians(ndat,xmu,xsigma,xgaus)
      implicit none
      integer ndat,iseed,j
      real*8 xmu,xsigma,xgaus(1:120000),pi,x1,x2,z1,z2
      
      pi=dacos(-1.0d0)

      iseed=16404850 !seed niub
      call srand(iseed)

      do j=1,ndat-1,2
      	x1=rand()
      	x2=rand()
      	z1=dsqrt(-2.d0*dlog(x1))*dcos(2.d0*pi*x2)
        z2=dsqrt(-2.d0*dlog(x1))*dsin(2.d0*pi*x2)
        xgaus(j)=z1*xsigma+xmu
        xgaus(j+1)=z2*xsigma+xmu
      end do

      return
      end subroutine

c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c                         FUNCIONS
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------

      real*8 function fun(x)
      implicit none
      real*8 x,pi,div,temp

      pi=dacos(-1.0d0)
      div=1/dsqrt(2.d0*pi)
      temp=div*dexp(-(x**2)/2.d0)
      fun=temp

      return
      end function
c-----------------------------------------------------------------------