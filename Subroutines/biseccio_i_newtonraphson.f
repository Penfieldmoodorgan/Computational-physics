c subrutina que implementa metode biseccio
c inputs: [a,b]:interval on busquem les arrels de la funcio, eps: precissió que volem
c outputs: niter: nombre de iteracions que han calgut per a arribar a la precissió desitjada
c xarrel:arrel (zero de la funció) trobada

      subroutine biT1(a,b,eps,niter,fun,xarrel)
      implicit none
      real*8 a,b,eps,fun,xarrel,x,error
      integer niter,n

      error=1000.d0
      n=0

      if (fun(a)*fun(b).ge.0.d0) then 
            write(*,*) 'La funcio no canvia de signe, no hi ha arrel'
      else
            do while ((error.gt.eps).and.(fun(x).ne.0.d0))
                 x=(a+b)/2.d0
                 n=n+1
                 if (fun(a)*fun(x).gt.0.d0) a=x
                 if (fun(a)*fun(x).lt.0.d0) b=x
                 error=dabs(a-b)
            end do
            niter=n
            xarrel=x
      end if

      return
      end subroutine

c subrutina que implementa metode newton-raphson
c escriu en fitxer x0, arrel, nombre iteració
c inputs: x0:punt del qual partim, eps:precissió,fun:funció,dfun:derivada de la funció
c outputs: niter: nombre de iteracions que han calgut per a arribar a la precissió desitjada
c xarrel:arrel (zero de la funció) trobada

      subroutine nrT1(x0,eps,niter,fun,dfun,xarrel)
      implicit none
      real*8 x0,eps,fun,dfun,xarrel,error,x1,x
      integer niter,n

600   format(e20.12,5x,e20.12,5x,i10)
500   format(9x,a,23x,a,15x,a)

      open(4,file='P4-18P-res.dat',status='unknown',position='append') 
c      write(2,500) 'Convergencia Newton-Raphson x0:',x0
c      write(2,*) ' '
      write(4,500) 'x0','Arrel','Iteracio'
      error=1000.d0
      x1=x0
      n=0
     
      do while(error.gt.eps)
      	x=x1-(fun(x1)/dfun(x1))
      	error=dabs(x-x1)
      	n=n+1
      	write(4,600) x0,x,n
      	x1=x
      end do
      write(4,*) ' '
      write(4,*) ' '
      
      niter=n
      xarrel=x

      close(4)

      return 
      end subroutine