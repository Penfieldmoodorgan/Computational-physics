      subroutine metropolis(s,l,t,m,time)
      implicit none
      real*8 t,time,valor,prob
      integer l,s(l,l),m,n,dh,canvi,i,j,a,b,c,d

      n=l**2 !nombre de espins

c tirem numeros aleatoris amb distribució U(1,L) per a x i y, com tenim U(0,1) fem canvi y=a+(b-a)*x on x € U(0,1)
c d'aquesta manera tenim els dos índex (x,y) escollits de manera aleatòria a partir de la distribució U(1,L) per triar uniformement a l'atzar un espí 

      i=int(1.d0+(l-1)*rand()) 
      j=int(1.d0+(l-1)*rand())

c fem avançar el temps físic com t->t+(1/N)
      time=time+(1.d0/dble(n))

c porposem canvi espin s(i,j)=-s(i,j)      
      canvi=-s(i,j)

c condicions periòdiques, no tenim malla en 2D plana sino sobre un torus      
      if(j.eq.1) then
            a=s(i,l)
      else
            a=s(i,j-1)
      end if

      if(j.eq.l) then
            b=s(i,1)
      else
            b=s(i,j+1)
      end if

      if(i.eq.l) then
            c=s(1,j)
      else
            c=s(i+1,j)
      end if

      if(i.eq.1) then
            d=s(l,j)
      else
            d=s(i-1,j)
      end if

c variació energia associada a un possible canvi (dh=2*s(i,j)*(s(i,j-1)+s(i,j+1)+s(i+1,j)+s(i-1,j)))
      dh=2*s(i,j)*(a+b+c+d)

c  si canvi energia menor que zero acceptem canvi,si canvi energia major o igual que zero acceptem canvi amb una probabilitat e^-dh/t
c probabilitat acceptar canvi p=e^-dh/t, prob no acceptar 1-p
      if(dh.lt.0) then 
            s(i,j)=canvi
c si l'espí canvia de negatiu a positiu, tindrem un canvi positiu en la magnetització, si l'espí canvia de positiu a negatiu, tindrem un canvi negatiu en la magnetització
c com la magnetització absoluta és igual a n*magnetització de l'espí i n=l^2, si el canvi en la magnetització en l'espí és de m->m +/- (2/l^2), el canvi en la magnetització absoluta serà M->M +/- 2
            if(canvi.gt.0) m=m+2 
            if(canvi.lt.0) m=m-2 
      else              
            prob=dexp(-dh/t)
            valor=rand()
            if(valor.lt.prob) then 
                  s(i,j)=canvi
                  if(canvi.gt.0) m=m+2
                  if(canvi.lt.0) m=m-2
            end if
      end if

      return
      end subroutine