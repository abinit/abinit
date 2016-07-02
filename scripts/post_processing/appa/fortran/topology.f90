!------------------------------------!
!-----PAIR DISTRIBUTION FUNCTION-----!
!------------------------------------!
subroutine pair_distribution(nbtime,position,rprim,f,inc,a,b,c,deltaR,rho,indexAtom1,indexAtom2,r,res1 &
,res2,size1,size2,size3,nba1,nba2)


  !input 
  integer :: size1,size2,size3,nba1,nba2
  integer, intent(in) :: inc,nbtime
  real, intent(in) :: deltaR,rho,f
  !f2py  optional  ,  depend(r)           ::  size1=len(r)
  !f2py  optional  ,  depend(position)    ::  size2=len(position[:,1,1])
  !f2py  optional  ,  depend(position)    ::  size3=len(position[1,:,1])
  !f2py  optional  ,  depend(indexAtom1)  ::  nba1=len(indexAtom1)
  !f2py  optional  ,  depend(indexAtom2)  ::  nba2=len(indexAtom2)
  real, intent(in)  ::  r(size1)
  real, intent(in)  ::  position(size2,size3,3),a(size2),b(size2),c(size2),rprim(3,3)
  integer, intent(in)  ::  indexAtom1(nba1),indexAtom2(nba2)
  !output
  real(8),intent(out)  ::  res1(size1), res2(size1) 
  !local variable
  integer :: bin,j,k,l,nb,xRmax,yRmax,zRmax,xx,yy,zz
  real(8) :: distance,x0,y0,z0,s, Pi,Rmax
  real(8), dimension(3)  ::  p1=(/1,1,1/),p2=(/1,1,1/)
  real(8) :: HIST(size1),G1(size1),G2(size1),G3(size1),Integ(size1) !,aprim(3),bprim(3),cprim(3)
     
  nb = 0

  j  = 1

  s  = min(a(size2),b(size2),c(size2))/2.

  HIST(:) = 0

  Pi = 4*atan(1.) !accurate value of pi

  do while(j<=nbtime) !loop on step

     Rmax = f*min(a(j),b(j),c(j))   !half boxes calculated.
     xRmax = int(2*Rmax/a(j) + 1/2.)  !number of step (box, and not half of boxes) for one side of the supercell outside of the normal cell.
     yRmax = int(2*Rmax/b(j) + 1/2.)
     zRmax = int(2*Rmax/c(j) + 1/2.)
 
!     aprim = rprim(1,:)/a(j)
!     bprim = rprim(2,:)/b(j)
!     cprim = rprim(3,:)/c(j)

     do k=1,size(indexAtom1)
        do l=1,size(indexAtom2)

           if (k /= l) then

              p1 = position(j,indexAtom1(k),:)
              p2 = position(j,indexAtom2(l),:)
              
              x = p1(1)-p2(1)
              y = p1(2)-p2(2)
              z = p1(3)-p2(3)

              x = x - a(j)*int(x/a(j))  !replace atoms in a box [-1:1]

              if (x/a(j) > 1/2.) then
                 x = x - a(j)
              endif

              if (x/a(j) < -1/2.) then
                 x = x + a(j)
              endif

              y = y - b(j)*int(y/b(j))  !replace atoms in a box [-1:1]

              if (y/b(j) > 1/2.) then
                 y = y - b(j)
              endif

              if (y/b(j) < -1/2.) then
                 y = y + b(j)
              endif

              z = z - c(j)*int(z/c(j))  !replace atoms in a box [-1:1]

              if (z/c(j) > 1/2.) then
                 z = z - c(j)
              endif

              if (z/c(j) < -1/2.) then
                 z = z + c(j)
              endif

              !use symetries to have the rdf further the box limit

              x0=x
              y0=y
              z0=z
              
              do xx=-xRmax,xRmax
                 x=x0+xx*a(j)
                 do yy=-yRmax,yRmax
                    y=y0+yy*b(j)
                    do zz=-zRmax,zRmax
                       z=z0+zz*c(j)

                       distance = (x**2+y**2+z**2)**0.5
                       
                       !distance = ((x*aprim(1)+y*aprim(2)+z*aprim(3))**2+ &
                       !     (x*bprim(1)+y*bprim(2)+z*bprim(3))**2+&
                       !     (x*cprim(1)+y*cprim(2)+z*cprim(3))**2)**0.5

                       if (distance < Rmax .and. distance /=0) then
                          
                          bin = int( distance / deltaR ) + 1
                          
                          if(bin <= size1) then                             

                             HIST(bin) = HIST(bin) + 1

                          endif

                       endif

                    enddo
                 enddo
              enddo
              
           endif

        end do
     end do
     j = j + inc
     nb = nb+1
  end do

  G1 =  HIST / (nb*size(indexAtom1)*size(indexAtom2)*rho * (4./3) * Pi  * ( (r+deltaR)**3 - r**3 ) )

  G2 =  HIST / (nb*size(indexAtom1)*size(indexAtom2)*rho * (4./3) * Pi  * ( (r+deltaR)**3 - r**3 ) )

  G3 =  HIST / (nb*size(indexAtom1)*size(indexAtom2)*rho * (4./3) * Pi  * ( (r+deltaR)**3 - r**3 ) )

!  G2 =  HIST / (nb*size(indexAtom1)*size(indexAtom2)*rho * Pi/12. *( (3 - 36*(r+deltaR)**2 + 32*(r+deltaR)**3)&
!       - (3 - 36*r**2 + 32.*r**3) ))
     
!  G3 = HIST / (nb*size(indexAtom1)*size(indexAtom2)*rho * ( (-Pi/4.+3*Pi*(r+deltaR)**2+&
!       (4*(r+deltaR)**2-2)**0.5+(1-12*(r+deltaR)**2)*atan((4*(r+deltaR)**2-2)**0.5)+16/3.*(r+deltaR)**3*&
!       atan(2*(r+deltaR)*(4*(r+deltaR)**2-3)/((4*(r+deltaR)**2-2)**0.5)*(4*(r+deltaR)**2+1))) - &
!       (-Pi/4.+3*Pi*r**2+(4*r**2-2)**0.5+(1-12*r**2)*atan((4*r**2-2)**0.5)+16/3.*r**3*&
!       atan(2*r*(4*r**2-3)/((4*r**2-2)**0.5)*(4*r**2+1))) )) 

  Integ = 3*r**2 / (2*nb*(size(indexAtom1))* ((r + deltaR/10000.)**3 - (r**3) )*10000.)

  do i=2,size1

     distance = (i-1)*deltaR
       
     res2(i) = res2(i-1) + Integ(i)*(HIST(i)+HIST(i-1))*deltaR
 
     if (distance <= s .or. distance > s*3**0.5) then
           
        res1(i) = G1(i)
           
     elseif (distance > s .and. distance <= s*2**0.5) then
        
        res1(i) = G2(i)
           
     elseif (distance > s*2**0.5 .and. distance <= s*3**0.5) then
           
        res1(i) = G3(i)

     endif

  enddo

end subroutine pair_distribution

!---------------------------!
!-----RDF DECONVOLUTION-----!
!-------------------------- !

subroutine rdf_deconvolution(nei,iter1,nbtime,position,inc,a,b,c,deltaR,indexAtom1,indexAtom2,r,res &
,size1,size2,size3,nba1,nba2)

  !nei gives the number of neighbors we consider

  !input 
  integer :: size1,size2,size3,nba1,nba2
  integer, intent(in) :: nbtime,inc,nei
  real, intent(in) :: deltaR
  !f2py  optional  ,  depend(r)           ::  size1=len(r)
  !f2py  optional  ,  depend(position)    ::  size2=len(position[:,1,1])
  !f2py  optional  ,  depend(position)    ::  size3=len(position[1,:,1])
  !f2py  optional  ,  depend(indexAtom1)  ::  nba1=len(indexAtom1)
  !f2py  optional  ,  depend(indexAtom2)  ::  nba2=len(indexAtom2)
  real, intent(in) :: r(size1)
  real, intent(in)  ::  position(size2,size3,3),iter1(size3,3),a(size2),b(size2),c(size2)
  integer, intent(in)  ::  indexAtom1(nba1),indexAtom2(nba2)
  !output
  real(8),intent(out)  ::  res(size1)
  !local variable
  integer :: bin,i,j,k,l,m,nb
  real(8) :: x,y,z,x1,y1,z1,x2,y2,z2,Pi,scal,ang,max,distance
  real(8) :: HIST(size1), Dist(size(indexAtom2)-1), Org(size2,size3-1)
  logical :: test

  HIST(:) = 0

  Pi = 4*atan(1.)

  do j=1,size(indexAtom1)
        
     m = 0
     max = 0
        
     Dist(:) = 0
           
     do k=1,size(indexAtom2)

        if (k /= j) then

           m = m + 1

           x = iter1(indexAtom1(j),1)-position(i,indexAtom2(k),1)
           y = iter1(indexAtom1(j),2)-position(i,indexAtom2(k),2)
           z = iter1(indexAtom1(j),3)-position(i,indexAtom2(k),3)
                 
           x = x - a(i)*int(x/a(i))

           if (x/a(i) > 1/2.) then
              x = x - a(i)
           endif

           if (x/a(i) < -1/2.) then
              x = x + a(i)
           endif

           y = y - b(i)*int(y/b(i))
           
           if (y/b(i) > 1/2.) then
              y = y - b(i)
           endif
              
           if (y/b(i) < -1/2.) then
              y = y + b(i)
           endif

           z = z - c(i)*int(z/c(i))
              
           if (z/c(i) > 1/2.) then
              z = z - c(i)
           endif
              
           if (z/c(i) < -1/2.) then
              z = z + c(i)
           endif

           distance = ( x**2 + y**2 + z**2 )**0.5
              
           !organisation of matrices Dist and Pos to have atoms sort with distances

           if (distance >= max) then
                 
              Dist(m) = distance
                 
              Org(j,m) = m

              max = distance

           else !distance < max

              test = .false. 
              l = m-1
                 
              do while (test .eqv. .false.)
                       
                 Dist(l+1) = Dist(l)
                       
                 Org(j,l+1) = Org(j,l)
                    
                 if (distance >= Dist(l-1)) then
                          
                    test = .true.

                    Dist(l) = distance
                    
                    Org(j,l) = m
                    
                 else

                    l = l-1

                 endif

              enddo

           endif

        endif
          
     enddo

  enddo

  nb = 0
  i  = 1

  do while (i <= nbtime)

     do j=1,size(indexAtom1)

        l = Org(j,nei)    !number of the neith neibor of the atom j (on the first iteration) at the iteration i.

        x = position(i,indexAtom1(j),1)-position(i,l,1)
        y = position(i,indexAtom1(j),2)-position(i,l,2)
        z = position(i,indexAtom1(j),3)-position(i,l,3)

        x = x - a(i)*int(x/a(i))

        if (x/a(i) > 1/2.) then
           x = x - a(i)
        endif

        if (x/a(i) < -1/2.) then
           x = x + a(i)
        endif

        y = y - b(i)*int(y/b(i))
              
        if (y/b(i) > 1/2.) then
           y = y - b(i)
        endif
              
        if (y/b(i) < -1/2.) then
           y = y + b(i)
        endif

        z = z - c(i)*int(z/c(i))
              
        if (z/c(i) > 1/2.) then
           z = z - c(i)
        endif
              
        if (z/c(i) < -1/2.) then
           z = z + c(i)
        endif

        distance = ( x**2 + y**2 + z**2 )**0.5

        bin = int(distance/deltaR) + 1

        if (bin <= size1) then

           HIST(bin) = HIST(bin) + 1
             
        endif

     enddo
     nb = nb + 1
     i = i + inc
  enddo

  res = HIST/nb


end subroutine rdf_deconvolution


!---------------------------------------!
!-----ANGULAR DISTRIBUTION FUNCTION-----!
!---------------------------------------!

subroutine angular_distribution(nei,nbtime,position,rprim,inc,a,b,c,delta,indexAtom1,indexAtom2,theta,res &
,size1,size2,size3,nba1,nba2)

  !nei gives the number of neighbors we consider

  !input 
  integer :: size1,size2,size3,nba1,nba2
  integer, intent(in) :: nbtime,inc
  real, intent(in) :: delta,nei
  !f2py  optional  ,  depend(theta)       ::  size1=len(theta)
  !f2py  optional  ,  depend(position)    ::  size2=len(position[:,1,1])
  !f2py  optional  ,  depend(position)    ::  size3=len(position[1,:,1])
  !f2py  optional  ,  depend(indexAtom1)  ::  nba1=len(indexAtom1)
  !f2py  optional  ,  depend(indexAtom2)  ::  nba2=len(indexAtom2)
  real, intent(in) :: theta(size1)
  real, intent(in)  ::  position(size2,size3,3),a(size2),b(size2),c(size2),rprim(3,3)
  integer, intent(in)  ::  indexAtom1(nba1),indexAtom2(nba2)
  !output
  real(8),intent(out)  ::  res(size1)
  !local variable
  integer :: bin,i,j,k,l,m,nb
  real(8) :: x,y,z,x1,y1,z1,x2,y2,z2,Pi,scal,ang,max,distance,d1,d2
  real(8) :: HIST(size1), Dist(size(indexAtom2)-1), Pos(size(indexAtom2)-1,3),aprim(3),bprim(3),cprim(3)
  logical :: test

  HIST(:) = 0

  nb = 0

  i = 1

  Pi = 4*atan(1.)

  do while(i<=nbtime) !loop on step

     aprim = rprim(1,:)/a(i)
     bprim = rprim(2,:)/b(i)
     cprim = rprim(3,:)/c(i)

     do j=1,size(indexAtom1)
        
        m = 0
        max = 0
        
        Dist(:) = 0
        Pos(:,:) = 0
           
        do k=1,size(indexAtom2)

           if (k /= j) then

              m = m + 1

              x = position(i,indexAtom1(j),1)-position(i,indexAtom2(k),1)
              y = position(i,indexAtom1(j),2)-position(i,indexAtom2(k),2)
              z = position(i,indexAtom1(j),3)-position(i,indexAtom2(k),3)
                 
              x = x - a(i)*int(x/a(i))

              if (x/a(i) > 1/2.) then
                 x = x - a(i)
              endif

              if (x/a(i) < -1/2.) then
                 x = x + a(i)
              endif

              y = y - b(i)*int(y/b(i))
              
              if (y/b(i) > 1/2.) then
                 y = y - b(i)
              endif
              
              if (y/b(i) < -1/2.) then
                 y = y + b(i)
              endif

              z = z - c(i)*int(z/c(i))
              
              if (z/c(i) > 1/2.) then
                 z = z - c(i)
              endif
              
              if (z/c(i) < -1/2.) then
                 z = z + c(i)
              endif

              distance = ((x*aprim(1)+y*aprim(2)+z*aprim(3))**2+ &
                            (x*bprim(1)+y*bprim(2)+z*bprim(3))**2+&
                            (x*cprim(1)+y*cprim(2)+z*cprim(3))**2)**0.5
              
              !organisation of matrices Dist and Pos to have atoms sort with distances

              if (distance >= max) then
                 
                 Dist(m) = distance
                 
                 Pos(m,1) = x
                 Pos(m,2) = y
                 Pos(m,3) = z

                 max = distance

              else !distance < max

                 test = .false. 
                 l = m-1
                 
                 do while (test .eqv. .false.)
                       
                    Dist(l+1) = Dist(l)
                       
                    Pos(l+1,1) = Pos(l,1)
                    Pos(l+1,2) = Pos(l,2)
                    Pos(l+1,3) = Pos(l,3)
                    
                    if (distance >= Dist(l-1)) then
                          
                       test = .true.

                       Dist(l) = distance
                       
                       Pos(l,1) = x
                       Pos(l,2) = y
                       Pos(l,3) = z
                       
                    else

                       l = l-1

                    endif

                 enddo

              endif

           endif
          
        enddo

        k = 0

        d1 = 0
        d2 = 0

        do while (d1 < nei)
                 
           k = k + 1

           x1 = Pos(k,1)
           y1 = Pos(k,2)
           z1 = Pos(k,3)

           d1 = Dist(k)

           l = 0
 
           do while (d2 < nei)

              l = l + 1
 
              x2 = Pos(l,1)
              y2 = Pos(l,2)
              z2 = Pos(l,3)
                    
              d2 = Dist(l)         

              if (k /= l) then
                 
                 scal = x1*x2+y1*y2+z1*z2
                 
                 !obtain angle from scalar product
                 
                 ang = acos(scal/(d1*d2))*180./Pi  
                 
                 bin = int(ang/delta)

                 if (bin >= 0) then

                    bin = bin + 1

                 endif

                 if (bin <= 0) then
              
                    HIST(size1) = HIST(size1) + 1
                    
                 endif

                 if (bin <= size1 .and. bin > 0) then
                    
                    HIST(bin) = HIST(bin) + 1

                 endif
                    
              endif

           enddo

        enddo

     enddo
     nb = nb + 1
     i = i + inc
  enddo

  res = HIST / nb

end subroutine angular_distribution



!----------------------------------------!
!-----NEIGHBOR DISTRIBUTION FUNCTION-----!
!----------------------------------------!

subroutine m_distribution(nei,nbtime,position,rprim,inc,it,a,b,c,indexAtom1,indexAtom2 &
,res_max,res_min,Data,size1,size2,nba1,nba2)

  !nei gives the number of neighbors we consider


  !input 
  integer :: size1,size2,nba1,nba2
  integer, intent(in) :: nbtime,inc,it,nei
  !f2py  optional  ,  depend(position)    ::  size1=len(position[:,1,1])
  !f2py  optional  ,  depend(position)    ::  size2=len(position[1,:,1])
  !f2py  optional  ,  depend(indexAtom1)  ::  nba1=len(indexAtom1)
  !f2py  optional  ,  depend(indexAtom2)  ::  nba2=len(indexAtom2)
  real, intent(in)  ::  position(size1,size2,3),a(size1),b(size1),c(size1)
  integer, intent(in)  ::  indexAtom1(nba1),indexAtom2(nba2),rprim(3,3)
  !output
  real(8),intent(out)  ::  res_max,res_min,Data(it) 
  !local variable
  integer :: i,j,k,l,m,n
  real(8) :: x,y,z,max,distance
  real(8), dimension(3)  ::  p1=(/1,1,1/),p2=(/1,1,1/)
  real(8) :: Dist(nba2-1),aprim(3),bprim(3),cprim(3)
  logical :: test

  res_max = 0
  res_min = 0

  Data(:) = 0

  i = 1

  n = 0

  do while(i<=nbtime) !loop on step

     aprim = rprim(1,:)/a(i)
     bprim = rprim(2,:)/b(i)
     cprim = rprim(3,:)/c(i)

     do j=1,size(indexAtom1)
        
        m = 0
        max = 0
                   
        Dist(:) = 0

        do k=1,size(indexAtom2)

           if (k /= j) then

              m = m + 1

              p1 = position(i,indexAtom1(j),:)
              p2 = position(i,indexAtom2(k),:)

              x = p1(1)-p2(1)
              y = p1(2)-p2(2)
              z = p1(3)-p2(3)
                 
              x = x - a(i)*int(x/a(i))

              if (x/a(i) > 1/2.) then
                 x = x - a(i)
              endif

              if (x/a(i) < -1/2.) then
                 x = x + a(i)
              endif

              y = y - b(i)*int(y/b(i))
              
              if (y/b(i) > 1/2.) then
                 y = y - b(i)
              endif

              if (y/b(i) < -1/2.) then
                 y = y + b(i)
              endif
              
              z = z - c(i)*int(z/c(i))
                 
              if (z/c(i) > 1/2.) then
                 z = z - c(i)
              endif

              if (z/c(i) < -1/2.) then
                 z = z + c(i)
              endif

              distance = ((x*aprim(1)+y*aprim(2)+z*aprim(3))**2+ &
                            (x*bprim(1)+y*bprim(2)+z*bprim(3))**2+&
                            (x*cprim(1)+y*cprim(2)+z*cprim(3))**2)**0.5

              !organisation of matrices Dist and Pos to have atoms sort with distances

              if (distance >= max) then
                    
                 Data((nba2-1)*n+m) = distance
                 Dist(m) = distance
                 max = distance

              else !distance < max

                 test = .false. 
                 l = m-1

                 do while (test .eqv. .false.)

                    Data((nba2-1)*n+l+1) = Data((nba2-1)*n+l)
                    Dist(l+1) = Dist(l)

                    if (distance >= Dist(l-1)) then
                          
                       test = .true.

                       Data((nba2-1)*n+l) = distance
                       Dist(l) = distance

                    else

                       l = l-1

                    endif
                    
                 enddo

              endif

           endif

        enddo

        if (nei /= 0) then

           if (i == 1 .and. j == 1) then

              res_max = Dist(nei)
              res_min = Dist(nei)
           
           else

              if (Dist(nei) > res_max) then

                 res_max = Dist(nei)
                 
              endif
              
              if (Dist(nei) < res_min) then
                 
                 res_min = Dist(nei)

              endif
              
           endif

        endif
        n = n+1
     enddo
     i = i + inc
  enddo

end subroutine m_distribution


subroutine n_distribution(nei,nbtime,data,inc,delta,indexAtom1,indexAtom2,dmin,r,res &
,it,size1,nba1,nba2)

  !nei gives the number of neighbors we consider

  !input 
  integer :: it,size1,nba1,nba2
  integer, intent(in) :: nbtime,inc,nei
  real, intent(in) :: delta,dmin
  !f2py  optional  ,  depend(data)        ::  it=len(data)
  !f2py  optional  ,  depend(r)           ::  size1=len(r)
  !f2py  optional  ,  depend(indexAtom1)  ::  nba1=len(indexAtom1)
  !f2py  optional  ,  depend(indexAtom1)  ::  nba2=len(indexAtom2)
  real, intent(in) :: r(size1),data(it)
  integer, intent(in)  ::  indexAtom1(nba1),indexAtom2(nba2)
  !output
  real(8),intent(out)  ::  res(size1)
  !local variable
  integer :: bin,i,j,k,nb
  real(8) :: HIST(size1)

  HIST(:) = 0

  nb = 0

  i = 1
  k = 0

  do while(i<=nbtime) !loop on step
     do j=1,size(indexAtom1)
        
        bin = int((data((nba2-1)*k+nei)-dmin)/delta) + 1

        if (bin <= size1) then

           HIST(bin) = HIST(bin) + 1
             
        endif

        k= k+1
        
     enddo
     nb = nb + 1
     i = i + inc
  enddo

  res = HIST / nb

end subroutine n_distribution



!---------------------------------------!
!-----NEIGHBOR PROBABILITY FUNCTION-----!
!---------------------------------------!

subroutine probability(nei,nbtime,position,rprim,a,b,c,indexAtom1,indexAtom2,res,size1,size2,nba1,nba2)

  !nei gives the number of neighbors we consider


  !input 
  integer :: size1,size2,nba1,nba2
  integer, intent(in) :: nbtime,nei
  !f2py  optional  ,  depend(position)    ::  size1=len(position[:,1,1])
  !f2py  optional  ,  depend(position)    ::  size2=len(position[1,:,1])
  !f2py  optional  ,  depend(indexAtom1)  ::  nba1=len(indexAtom1)
  !f2py  optional  ,  depend(indexAtom2)  ::  nba2=len(indexAtom2)
  real, intent(in)  ::  position(size1,size2,3),a(size1),b(size1),c(size1),rprim(3,3)
  integer, intent(in)  ::  indexAtom1(nba1),indexAtom2(nba2)
  !output
  real(8),intent(out)  ::  res(nei)
  !local variable
  integer :: i,j,k,l,m,n
  real(8) :: x,y,z,max,distance
  real(8), dimension(3)  ::  p1=(/1,1,1/),p2=(/1,1,1/)
  real(8) :: HIST(nei), Dist(size(indexAtom2)-1), Org(size(indexAtom1),size(indexAtom2)-1,2), & 
Reorg(size(indexAtom1),size(indexAtom2)-1,2),aprim(3),bprim(3),cprim(3)
  logical :: test

  HIST(:) = 0

  do i=1,2

     aprim = rprim(1,:)/a(i)
     bprim = rprim(2,:)/b(i)
     cprim = rprim(3,:)/c(i)

     do j=1,size(indexAtom1)
           
        m = 0
        max = 0
           
        Dist(:) = 0

        do k=1,size(indexAtom2)

           if (k /= j) then

              m = m + 1

              if (i == 1) then

                 p1 = position(1,indexAtom1(j),:)
                 p2 = position(1,indexAtom2(k),:)
              
              else

                 p1 = position(nbtime,indexAtom1(j),:)
                 p2 = position(nbtime,indexAtom2(k),:)

              endif
              
              x = p1(1)-p2(1)
              y = p1(2)-p2(2)
              z = p1(3)-p2(3)
                 
              x = x - a(i)*int(x/a(i))

              if (x/a(i) > 1/2.) then
                 x = x - a(i)
              endif

              if (x/a(i) < -1/2.) then
                 x = x + a(i)
              endif

              y = y - b(i)*int(y/b(i))

              if (y/b(i) > 1/2.) then
                 y = y - b(i)
              endif

              if (y/b(i) < -1/2.) then
                 y = y + b(i)
              endif

              z = z - c(i)*int(z/c(i))
              
              if (z/c(i) > 1/2.) then
                 z = z - c(i)
              endif

              if (z/c(i) < -1/2.) then
                 z = z + c(i)
              endif

              distance = ((x*aprim(1)+y*aprim(2)+z*aprim(3))**2+ &
                            (x*bprim(1)+y*bprim(2)+z*bprim(3))**2+&
                            (x*cprim(1)+y*cprim(2)+z*cprim(3))**2)**0.5

              !organisation of matrices Dist and Pos to have atoms sort with distances

              if (distance >= max) then
                 
                 Dist(m) = distance

                 Org(j,m,i) = m

                 max = distance

              else !distance < max

                 test = .false. 
                 l = m-1
                    
                 do while (test .eqv. .false.)
                       
                    Dist(l+1) = Dist(l)
                       
                    Org(j,l+1,i) = Org(j,l,i)
                       
                    if (distance >= Dist(l-1)) then
                          
                       test = .true.

                       Dist(l) = distance
                          
                       Org(j,l,i) = m

                    else

                       l = l-1

                    endif
                    
                 enddo

              endif

           endif
           
        enddo

        if (i == 2) then

           do k=1,size(indexAtom2)-1

              l = 1
              test = .false.
              do while (test .eqv. .false.)

                 if (Org(j,k,1) == Org(j,l,2)) then

                    Reorg(j,k,1) = k
                    Reorg(j,k,2) = l
                    
                    test = .true.

                 endif
                 l=l+1
              enddo
              
           enddo

        endif

     enddo
  enddo

  
  do i=1,nei
     do j=1,size(indexAtom1)
        n = 0
        do k=1,i
        
           if (Reorg(j,k,1) <= i .and. Reorg(j,k,2) <= i) then

              n = n + 1

           endif

        enddo

        if (n == i) then

           HIST(i) = HIST(i) + 1
           
        endif

     enddo
  enddo

  res = HIST / size(indexAtom1)
     
end subroutine probability
