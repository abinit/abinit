!---------------------------------!
!-----Atomic Position Symetry-----!
!---------------------------------!

subroutine atom_sym(nbtime,position,a,b,c,indexAtom1,Pos1,Pos2,Pos3,size1,size2,nba1)

  !input 
  integer :: size1,size2,nba1
  integer, intent(in) :: nbtime
  !f2py  optional  ,  depend(position)    ::  size1=len(position[:,1,1])
  !f2py  optional  ,  depend(position)    ::  size2=len(position[1,:,1])
  !f2py  optional  ,  depend(indexAtom1)  ::  nba1=len(indexAtom1)
  real, intent(in)  ::  position(size1,size2,3)
  integer, intent(in)  ::  indexAtom1(nba1)
  real, intent(in) :: a(size1),b(size1),c(size1)
  !output
  real(8),intent(out)  ::  Pos1(nbtime*nba1),Pos2(nbtime*nba1),Pos3(nbtime*nba1)
  !local variable
  real(8) :: x,y,z
  integer :: i,j,k

  k = 1

  do  i=1,nbtime
     do j=1,nba1
 
       x = position(i,indexAtom1(j),1)
       y = position(i,indexAtom1(j),2)
       z = position(i,indexAtom1(j),3)

       if (x >= 0) then  !replace atoms in a box [0:1]
         
         Pos1(k) = x - a(i)*int(x/a(i))
         
       else
         
         Pos1(k) = x - a(i)*int(x/a(i)-1)
         
       endif
       
       if (y >= 0) then  !replace atoms in a box [0:1]
         
         Pos2(k) = y - b(i)*int(y/b(i))
         
       else
         
         Pos2(k) = y - b(i)*int(y/b(i)-1)
         
       endif
       
       if (z >= 0) then  !replace atoms in a box [0:1]
         
         Pos3(k) = z - c(i)*int(z/c(i))
         
       else
         
         Pos3(k) = z - c(i)*int(z/c(i)-1)
         
       endif
        
        k = k+1
        
      enddo
    enddo

end subroutine atom_sym


!---------------------------------!
!-----Atomic Position Average-----!
!---------------------------------!

subroutine atom_ave(nbtime,position,indexAtom1,pos1,pos2,pos3,size1,size2,nba1)

  !input 
  integer :: size1,size2,nba1
  integer, intent(in) :: nbtime
  !f2py  optional  ,  depend(position)    ::  size1=len(position[:,1,1])
  !f2py  optional  ,  depend(position)    ::  size2=len(position[1,:,1])
  !f2py  optional  ,  depend(indexAtom1)  ::  nba1=len(indexAtom1)
  real, intent(in)  ::  position(size1,size2,3)
  integer, intent(in)  ::  indexAtom1(nba1)
  !output
  real(8),intent(out)  ::  pos1(nba1),pos2(nba1),pos3(nba1)
  !local variable
  real :: x,y,z,xred,yred,zred
  integer :: i,j

  do  i=1,nbtime
     do j=1,nba1

        x = position(i,indexAtom1(j),1)
        y = position(i,indexAtom1(j),2)
        z = position(i,indexAtom1(j),3)

        pos1(j) = pos1(j) + x/nbtime
        pos2(j) = pos2(j) + y/nbtime
        pos3(j) = pos3(j) + z/nbtime

     enddo
  enddo

end subroutine atom_ave

!---------------------------------!
!-----Atomic Position Standard----!
!---------------------------------!

subroutine atom_std(nbtime,position,indexAtom1,pos1,pos2,pos3,size1,size2,nba1)

  !input 
  integer :: size1,size2,nba1
  integer, intent(in) :: nbtime
  !f2py  optional  ,  depend(position)    ::  size1=len(position[:,1,1])
  !f2py  optional  ,  depend(position)    ::  size2=len(position[1,:,1])
  !f2py  optional  ,  depend(indexAtom1)  ::  nba1=len(indexAtom1)
  real, intent(in)  ::  position(size1,size2,3)
  integer, intent(in)  ::  indexAtom1(nba1)
  !output
  real(8),intent(out)  ::  pos1(nbtime*nba1),pos2(nbtime*nba1),pos3(nbtime*nba1)
  !local variable
  real :: x,y,z,xred,yred,zred
  integer :: i,j,k

  k = 1
  do  i=1,nbtime
     do j=1,nba1

        x = position(i,indexAtom1(j),1)
        y = position(i,indexAtom1(j),2)
        z = position(i,indexAtom1(j),3)

        pos1(k) = x
        pos2(k) = y
        pos3(k) = z

        k=k+1
     enddo
  enddo

end subroutine atom_std

