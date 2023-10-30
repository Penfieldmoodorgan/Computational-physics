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
