c subrutina que calcula un pas del metode runge kutta 4 per a un sistema de 
c nequs eq de 1er ordre acoblades

c     ·Runge kutta 4t ordre:
c         k1=h*f(x0+y0)
c         k2=h*f(x0+0.5*h,y0+0.5*k1)
c         k3=h*f(x0+0.5*h,y0+0.5*k2)
c         k4=h*f(x0+h,y0+k3)
c      y(x0+h)=y0+(k1+2*k2+2*k3+k4)

c aqui h=dx,f(x,y)=dy/dx,y=y'',x=x0,yyin=y0 vector

      subroutine myRK4(x,dx,nequs,yyin,yyout)
      implicit none
      real*8 yyin(nequs),yyout(nequs),yin(nequs),x,dx
      real*8 k1(nequs),k2(nequs),k3(nequs),k4(nequs),ksuma(nequs)
      integer nequs,i

c inicialitzem vectors a 0
      yyout=0.d0
      k1=0.d0
      k2=0.d0
      k3=0.d0
      k4=0.d0
      ksuma=0.d0
      yin=0.d0
      
c----obtencio k1--------  
      call derivades (nequs,x,yyin,k1)
c---obtencio k2-----------
      do i=1,nequs
      	yin(i)=yyin(i)+(dx*0.5d0*k1(i))
      end do

      call derivades (nequs,x+0.5d0*dx,yin,k2)
c---obtencio k3----------- 
      do i=1,nequs
      	yin(i)=yyin(i)+(dx*0.5d0*k2(i))
      end do

      call derivades (nequs,x+0.5d0*dx,yin,k3)
c---obtencio k4-----------      
      do i=1,nequs
      	yin(i)=yyin(i)+(dx*k3(i))
      end do

      call derivades (nequs,x+dx,yin,k4)

c-----calcul vector yyout    
      do i=1,nequs
      	ksuma(i)=k1(i) + 2.d0*k2(i) + 2.d0*k3(i)+k4(i)
      	yyout(i)=yyin(i)+(dx/6.d0)*ksuma(i)
      end do

      return
      end subroutine
