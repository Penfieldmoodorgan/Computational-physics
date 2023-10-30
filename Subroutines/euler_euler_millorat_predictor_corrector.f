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

c mètode de euler modificat 

      subroutine emodificat(npas,fun,tmax,tmin,y0,f0,y2,f2)
      implicit none
      real*8 y2(20000),f2(20000),tmax,dt,y0,f0,fun,tmin
      integer npas,i

      y2=0.d0 ! inicialitzem a zero
      f2=0.d0

      dt=(tmax-tmin)/dble(npas-1)

      y2(1)=y0
      f2(1)=f0
      y2(2)=y2(1)+dt*f2(1)
      f2(2)=f2(1)+dt*fun(y2(1)) 

      do i=3,npas

            y2(i)=y2(i-2)+2.d0*dt*f2(i-1)
            f2(i)=f2(i-2)+2.d0*dt*fun(y2(i-1))

      end do
      
      return
      end subroutine


c mètode corrector-predictor (euler millorat)

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
