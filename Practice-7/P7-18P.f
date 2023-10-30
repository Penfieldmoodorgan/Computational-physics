c -----------------------------------------------------------------------
C                      PRÀCTICA 7
c-----------------------------------------------------------------------
c                  Pèndol simple
C ----------------------------------------------------------------------

      program pre7
      implicit none
      real*8 m,l,g,tn,wn,pi,tmax,y0,f0,tmin,dt,f01,f02
      real*8 fun,epot,ekin,aprox
      real*8 kin,pot,total,kin2,pot2,total2
      real*8 y(2000),f(2000),y2(50001),f2(50001),y22(50001),f22(50001)
      integer j,k,npas,pascon(4)
      external fun,ekin,epot,aprox
      COMMON/DADES/m,l,g

      pi=dacos(-1.0d0) 
      tmin=0.d0
      pascon=(/200,600,4000,50000/)

c-----------------FORMATS-----------------------------------------------
c ----------------------------------------------------------------------
300   format(2x,7(e20.12,5x))
400   format(12x,a,12x,a,6x,a,3x,a,7x,a,5x,a,5x,a)
100   format(2x,5(e20.12,5x))
600   format(a,2x,i10)
700   format(10x,a,2x,e20.12)
200   format(10x,a,19x,a,3(15x,a))
500   format(29x,a,11x,a)
800   format(45x,a,29x,a)
250   format(10x,a,16x,a)
c ---------------------------------------------------------------------

c parametres pendol descrit per eq diferencial l*alfa''=-g*sin alfa
c alfa=angle, alfa''= acceleracio angular

      m=0.95d0 ! massa en kg
      l=1.05d0 !longitud en metres
      g=1.66d0 ! gravetat en metres/segons^2
      wn=dsqrt(g/l) !frequencia angular
      tn=(2.d0*pi)/wn !periode 
      tmax=6.d0*tn 

c apartat a----petites oscilacions

      npas=1500 
      y0=0.15d0 !alfa(0) en rad
      f0=0.d0 ! alfa'(0) en rad/s
      dt=(tmax-tmin)/dble(npas) 

      open(1,file='P7-18P-res.dat',status='unknown')
      write(1,*) ' -------------APARTAT A------------'
      write(1,*) ' '
      write(1,600) ' Numero de passos de temps:',npas 
      write(1,*) ' '
      write(1,500) 'Valors obtinguts amb metode euler cru: ','Valors obt
     +inguts amb metode predictor/corrector:'
      write(1,*) ' '
      write(1,200) 't(s)','alfa(rad)','dalfa(rad/s)','alfa(rad)','dalfa(
     +rad/s)'
      write(1,*) ' '
      
      call euler(npas,fun,tmax,tmin,y0,f0,y,f)
      call corrector(npas,fun,tmax,tmin,y0,f0,y2,f2)

      do j=1,npas+1
      	write(1,100) dble(j-1)*dt,y(j),f(j),y2(j),f2(j)
      end do

      write(1,*) ' '
      write(1,*) ' '

c---------aproximacio petites oscilacions: sin(alfa)->alfa----
      
      write(1,*) 'Aproximació petites oscilacions: sin(alfa) -> alfa'
      write(1,*) ' '
      write(1,500) 'Valors obtinguts amb metode euler cru: ','Valors obt
     +inguts amb metode predictor/corrector:'
      write(1,*) ' '
      write(1,200) 't(s)','alfa(rad)','dalfa(rad/s)','alfa(rad)','dalfa(
     +rad/s)'
      write(1,*) ' '
      
      call euler(npas,aprox,tmax,tmin,y0,f0,y,f)
      call corrector(npas,aprox,tmax,tmin,y0,f0,y2,f2)

      do j=1,npas+1
            write(1,100) dble(j-1)*dt,y(j),f(j),y2(j),f2(j)
      end do

      write(1,*) ' '
      write(1,*) ' '

c apartat b --- oscilacions grans

      y0=pi-0.15d0 !alfa(0) en rad
      f0=0.d0 ! alfa'(0) en rad/s
      npas=1500
      dt=(tmax-tmin)/dble(npas) 

      open(1,file='P7-18P-res.dat',status='unknown')
      write(1,*) ' -------------APARTAT B------------'
      write(1,*) ' '
      write(1,600) ' Numero de passos de temps:',npas 
      write(1,*) ' '
      write(1,500) 'Valors obtinguts amb metode euler cru: ','Valors obt
     +inguts amb metode predictor/corrector:'
      write(1,*) ' '
      write(1,200) 't(s)','alfa(rad)','dalfa(rad/s)','alfa(rad)','dalfa(
     +rad/s)'
      write(1,*) ' '
      
      call euler(npas,fun,tmax,tmin,y0,f0,y,f)
      call corrector(npas,fun,tmax,tmin,y0,f0,y2,f2)

      do j=1,npas+1
      	write(1,100) dble(j-1)*dt,y(j),f(j),y2(j),f2(j)
      end do

      write(1,*) ' '
      write(1,*) ' '

c apartat c----------energia

      y0=pi-0.015d0
      f0=0.1d0
      npas=1500
      dt=(tmax-tmin)/dble(npas) 

      write(1,*) ' -------------APARTAT C------------'
      write(1,*) ' '
      write(1,600) ' Numero de passos de temps:',npas 
      write(1,*) ' '
      

      call euler(npas,fun,tmax,tmin,y0,f0,y,f)
      call corrector(npas,fun,tmax,tmin,y0,f0,y2,f2)

      write(1,700) 'alfa:', y0
      write(1,*) ' '
      write(1,800) 'Valors obtinguts amb metode euler cru: ','Valors obt
     +inguts amb metode predictor/corrector:'
      write(1,*) ' '
      write(1,400)'t(s)','Energia cinetica (J)','Energia potencial (J)',
     +'Energia total (J)','Energia cinetica (J)','Energia potencial (J)'
     +,'Energia total (J)'
      write(1,*)' '
  
      do j=1,npas+1
      	kin=ekin(f(j))
      	kin2=ekin(f2(j))
      	pot=epot(y(j))
      	pot2=epot(y2(j))
      	total=kin+pot
      	total2=kin2+pot2
      	write(1,300) dt*(j-1),kin,pot,total,kin2,pot2,total2
      end do
      write(1,*) ' '
      write(1,*) ' '

c apartat d----transicio

      y0=0.d0
      f01=2.d0*wn+0.1d0
      f02=2.d0*wn-0.1d0
      tmax=12.d0*tn
      npas=5000
      dt=(tmax-tmin)/dble(npas) 

      call corrector(npas,fun,tmax,tmin,y0,f01,y2,f2)
      call corrector(npas,fun,tmax,tmin,y0,f02,y22,f22)

      write(1,*) ' -------------APARTAT D------------'
      write(1,*) ' '
      write(1,600) ' Numero de passos de temps:',npas 
      write(1,*) ' '
      write(1,*) 'Valors obtinguts amb metode predictor/corrector: '
      write(1,*) ' '
      write(1,200) 't(s)','alfa+(rad)','dalfa+(rad/s)','alfa-(rad)','dal
     +fa-(rad/s)'
      write(1,*) ' '

      do j=1,npas+1
      	write(1,100) dt*(j-1),y2(j),f2(j),y22(j),f22(j)
      end do
      
      write(1,*) ' '
      write(1,*) ' '

c apartat e-----------convergencia del metode

      write(1,*) ' -------------APARTAT E------------'
      write(1,*) ' '
      write(1,*) 'Valors obtinguts amb metode predictor/corrector: '
      write(1,*) ' '

      do k=1,4
        y0=3.d0       
        f0=0.d0        
        tmax=10.d0*tn
      	npas=pascon(k)
        dt=(tmax-tmin)/dble(npas)
      	call corrector(npas,fun,tmax,tmin,y0,f0,y2,f2)
      	write(1,600) ' Numero de passos de temps:',npas 
        write(1,*) ' '
        write(1,250) 't(s)','Energia total(J)'
        write(1,*) ' '
      	do j=1,npas+1
      		kin=ekin(f2(j))
      		pot=epot(y2(j))
      		total=kin+pot
      		write(1,100) dt*(j-1),total
      	end do
        write(1,*) ' '
        write(1,*) ' '
      end do

      close(1)

c extra----animacio
      
      open(2,file='anime.dat',status='unknown')

      y0=pi-0.15d0 !alfa(0) en rad
      f0=0.d0 ! alfa'(0) en rad/s
      npas=1500
      dt=(tmax-tmin)/dble(npas) 
     
      call euler(npas,fun,tmax,tmin,y0,f0,y,f)
      call corrector(npas,fun,tmax,tmin,y0,f0,y2,f2)

      do j=1,npas+1
      	write(2,*) dble(j-1)*dt,y(j),f(j),y2(j),f2(j)
      end do

      write(2,*) ' '
      write(2,*) ' '

      close(2)

      return
      end program


c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c              SUBRUTINES
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------

c mètode de euler per a calcular y i derivada primera de y en funcio del temps
c pas de temps=dt, num passos de temps=npas
c yo=fucio en 0, f0= derivada primera en 0

      subroutine euler(npas,fun,tmax,tmin,y0,f0,y,f)
      implicit none
      real*8 y(2000),f(2000),tmax,dt,y0,f0,fun,tmin
      integer npas,i

      y=0.d0  ! inicialitzem a zero
      f=0.d0

      dt=(tmax-tmin)/dble(npas)

      y(1)=y0
      f(1)=f0

      do i=2,npas+1

      	y(i)=y(i-1)+dt*f(i-1)    
      	f(i)=f(i-1)+dt*fun(y(i-1))

      end do
      
      return
      end subroutine

c mètode de euler MILLORAT 

      subroutine corrector(npas,fun,tmax,tmin,y0,f0,y2,f2)
      implicit none
      real*8 y2(50001),f2(50001),tmax,dt,y0,f0,fun,tmin,y1,fa,ya,f1
      integer npas,i

      y2=0.d0 ! inicialitzem a zero
      f2=0.d0

      dt=(tmax-tmin)/dble(npas)

      ya=y0
      fa=f0
c      y1=y0+dt*f0
c      f1=f0+dt*fun(y0)

      do i=1,npas+1

      	y1=ya+dt*fa
        f1=fa+dt*fun(ya)

      	y2(i)=ya+(dt/2.d0)*(fa+f1)
      	f2(i)=fa+(dt/2.d0)*(fun(ya)+fun(y1))

      	ya=y2(i)
      	fa=f2(i)

      end do
      
      return
      end subroutine


c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c                 FUNCIONS      
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------

c ------ funcio pendol
c ------ input x=angle

      real*8 function fun(x)
      implicit none
      real*8 x,g,l,temp,m
      COMMON/DADES/m,l,g

      temp=-(g/l)*dsin(x)
      fun=temp

      return
      end function

c ------ funcio pendol aproximacio petites oscilacions
c ------ input x=angle

      real*8 function aprox(x)
      implicit none
      real*8 x,g,l,temp,m
      COMMON/DADES/m,l,g

      temp=-(g/l)*x
      aprox=temp

      return
      end function

c ------ funcio energia cinetica
c ------ input x=velocitat angular

      real*8 function ekin(x)
      implicit none
      real*8 x,m,temp,l,g
      COMMON/DADES/m,l,g

      temp=0.5d0*m*(l**2)*(x**2)
      ekin=temp

      return
      end function

c ------ funcio energia potencial
c ------ input x=angle

      real*8 function epot(x)
      implicit none
      real*8 x,m,temp,g,l
      COMMON/DADES/m,l,g

      temp=-(m*g*l*dcos(x))
      epot=temp

      return
      end function

c-----------------------------------------------------------------------