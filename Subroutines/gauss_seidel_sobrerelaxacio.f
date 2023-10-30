c subrutina un pas de gauss-seidel o sobrerelaxacio. Elecció mètode controlada
c per variable icontrol. icontrol=1 -> gauss, icontrol=2 ->sobrerelaxació
c ti és la matriu inicial,tf es la matriu final i sigma és la tolerància
c xpunts és el nombre de punts en que hem dividit x
c ypunts és el nombre de punts en que hem dividit y
c (x0,y0) punt que ens interessa estudiar
c AMB FOGONS!

      subroutine metode(icontrol,ypunts,xpunts,sigma,x0,y0,tini,tf)
      implicit none
      real*8 a(0:xpunts,0:ypunts),b(0:xpunts,0:ypunts) ! a i b son matrius temporals per a operacions intermitjes
      real*8 tini(0:xpunts,0:ypunts),tf(0:xpunts,0:ypunts)
      real*8 ti(0:xpunts,0:ypunts)
      real*8 sigma,epsilon,w,error,h,ro,x,y,x0,y0
      integer icontrol,j,k,xpunts,ypunts,i
      common/dades/h

c     ------------FORMATS-----------------------------------------------
c     ------------------------------------------------------------------
100   format(6x,i10,18x,e20.12)
300   format(12x,a,23x,a)
c     ------------------------------------------------------------------
      
      w=1.35d0

c inicialitzem matrius a zero      
      tf=0.d0
      a=0.d0
      b=0.d0

c passem dades de tini a ti per a treballar amb aquesta última i no sobreescriure la primera
      do i=0,xpunts
      	do j=0,ypunts
      		ti(i,j)=tini(i,j)
      	end do
      end do
      
c escrivim a la matriu les condcions de contorn
      do i=0,xpunts
      	tf(i,0)=ti(i,0)
      	tf(i,ypunts)=ti(i,ypunts)
      	do j=1,ypunts
      		tf(0,j)=ti(0,j)
      		tf(xpunts,j)=ti(xpunts,j)
      	end do
      end do

      epsilon=1.d0 ! asignem valor gran a epsilon per entrar al bucle 
      k=0 !nombre iteració

      open(2,file='P9-18P.dat',status='unknown',position='append')
      if(icontrol.eq.1) write(2,*)'Mètode Gauss-Seidel'
      if(icontrol.eq.2) write(2,*)'Mètode Sobrerelaxació'
      write(2,*) ' '
      write(2,300) 'Iteració','Valor'
      write(2,*) ' '

c calculem punts centrals      
      do while(epsilon.gt.sigma)
        epsilon=0.d0 ! ara que ja estem a dins del bucle inicialitzem epsilon a zero 
        do i=1,ypunts-1 !recorre les y sense passar per valors determinats per cond contorn
          y=0.d0+i*h
          do j=1,xpunts-1 !recorre les x sense passar per valors determinats per cond contorn
            x=0.d0+j*h
            call densitat(x,y,ro) !calculem la densitat de temperatura al punt desitjat
            a(j,i)=(ti(j+1,i)+ti(j,i+1)+(h**2)*ro)
            if(icontrol.eq.1) then
              tf(j,i)=(a(j,i)+ti(j-1,i)+ti(j,i-1))/4.d0
            else if(icontrol.eq.2) then
              b(j,i)=((a(j,i)+tf(j-1,i)+tf(j,i-1))/4.d0)-ti(j,i)    
              tf(j,i)=ti(j,i)+(w*b(j,i))
            end if
            error=dabs(ti(j,i)-tf(j,i))
            if(error.gt.epsilon) epsilon=error
            if((x.eq.x0).and.(y.eq.y0)) write(2,100) k,tf(j,i)
          end do
        end do
        do i=1,ypunts-1
          do j=1,xpunts-1
            ti(j,i)=tf(j,i)
          end do
        end do
        k=k+1
      end do
      write(2,*) ' '
      write(2,*) ' '

      close(2)

      return 
      end subroutine

c sense fogons

      subroutine metode1(icontrol,ypunts,xpunts,sigma,tini,tf)
      implicit none
      real*8 a(0:xpunts,0:ypunts),b(0:xpunts,0:ypunts) ! a i b son matrius temporals per a operacions intermitjes
      real*8 tini(0:xpunts,0:ypunts),tf(0:xpunts,0:ypunts)
      real*8 ti(0:xpunts,0:ypunts)
      real*8 sigma,epsilon,w,error
      integer icontrol,j,k,xpunts,ypunts,i
      
      w=1.35d0

c inicialitzem matrius a zero      
      tf=0.d0
      a=0.d0
      b=0.d0

c passem dades de tini a ti per a treballar amb aquesta última i no sobreescriure la primera
      do i=0,xpunts
        do j=0,ypunts
          ti(i,j)=tini(i,j)
        end do
      end do
      
c escrivim a la matriu les condcions de contorn
      do i=0,xpunts
        tf(i,0)=ti(i,0)
        tf(i,ypunts)=ti(i,ypunts)
        do j=1,ypunts
          tf(0,j)=ti(0,j)
          tf(xpunts,j)=ti(xpunts,j)
        end do
      end do

      epsilon=1.d0 ! asignem valor gran a epsilon per entrar al bucle 
      k=0 !nombre iteració

c calculem punts centrals      
      do while(epsilon.gt.sigma)
        epsilon=0.d0 ! ara que ja estem a dins del bucle inicialitzem epsilon a zero 
        do i=1,ypunts-1 !recorre les y sense passar per valors determinats per cond contorn
          do j=1,xpunts-1 !recorre les x sense passar per valors determinats per cond contorn
            a(j,i)=(ti(j+1,i)+ti(j,i+1))
            if(icontrol.eq.1) then
              tf(j,i)=(a(j,i)+ti(j-1,i)+ti(j,i-1))/4.d0
            else if(icontrol.eq.2) then
              b(j,i)=((a(j,i)+tf(j-1,i)+tf(j,i-1))/4.d0)-ti(j,i)    
              tf(j,i)=ti(j,i)+(w*b(j,i))
            end if
            error=dabs(ti(j,i)-tf(j,i))
            if(error.gt.epsilon) epsilon=error
          end do
        end do
        do i=1,ypunts-1
          do j=1,xpunts-1
            ti(j,i)=tf(j,i)
          end do
        end do
        k=k+1
      end do

      return 
      end subroutine


c     versio 2- gauss,jacobi,sobrerelaxacio

c subrutina que implementa el mètode Gauss-Siedel o sobrerelaxació
c Elecció mètode controlada per variable icontrol. 
c icontrol=1 -> gauss-siedel, icontrol=2 -> jacobi, icontrol=3 -> sobrerelaxació
c ti és la matriu inicial,tf es la matriu final i sigma és la tolerància
c xpunts és el nombre de punts en que hem dividit x

      subroutine Gauss(icontrol,xpunts,ro,sigma,ti,tf)
      implicit none
      real*8 ro,sigma,h,w,epsilon,error,x
      real*8 ti(0:xpunts),tf(0:xpunts),b(0:xpunts)
      integer xpunts,icontrol,j,k
      common/dades/h

      w=1.5d0

c inicialitzem vectors a zero
      tf=0.d0
      b=0.d0

c escrivim al vector les condcions de contorn
      tf(0)=ti(0)
      tf(xpunts)=ti(xpunts)

      epsilon=1.d0 ! asignem valor gran a epsilon per entrar al bucle 
      k=0 !nombre iteració

c calculem punts centrals      
      do while(epsilon.gt.sigma)
        epsilon=0.d0 ! ara que ja estem a dins del bucle inicialitzem epsilon a zero 
        do j=1,xpunts-1 !recorre les x sense passar per valors determinats per cond contorn
          x=0.d0+j*h
          if(icontrol.eq.1) then ! Gauss-Siedel
            tf(j)=(ti(j+1)+ti(j-1)+(h**2)*ro(x))/2.d0
            else if(icontrol.eq.2) then !Jacobi
              tf(j)=(ti(j+1)+tf(j-1)+(h**2)*ro(x))/2.d0
          else if(icontrol.eq.3) then !Sobrerelaxació
            b(j)=((ti(j+1)+tf(j-1)+(h**2)*ro(x))/2.d0)-ti(j)    
            tf(j)=ti(j)+(w*b(j))
          end if                                         
          error=dabs(ti(j)-tf(j))
          if(error.gt.epsilon) epsilon=error
        end do
        do j=1,xpunts-1
          ti(j)=tf(j)
        end do
        k=k+1
      end do

      return
      end subroutine