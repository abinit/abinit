!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_per_cond
!! NAME
!!  m_per_cond
!!
!! FUNCTION
!!  This module contains basic tools for periodic conditions traitement.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2020 ABINIT group (MM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_per_cond

 use defs_basis
 use m_abicore

 implicit none

 private

 public :: per_cond   ! Return an arithmetic progression
 public :: per_dist   ! Return the distance on a periodic grid


 interface per_cond
  module procedure per_cond_re
  module procedure per_cond_int
  module procedure per_cond_int3
 end interface per_cond

 interface per_dist
  module procedure per_dist_int
  module procedure per_dist_int1
  module procedure per_dist_int3
 end interface per_dist

CONTAINS  !===========================================================
!!***

!!****f* m_per_cond/per_cond_re
!! NAME
!!  per_cond_re
!!
!! FUNCTION
!! Given a 2d-array of integer initial(3,nb), it calulates the values
!! of any of the array in the periodic orthogonal discretized
!! grid begining in 0 and of lengths dim_grid(3)
!!
!!
!! INPUTS
!!  initial(1:3,0:nb-1)=initial point
!!  dim_grid(1:3)= box lengths
!!  nb= dimension
!!  metric=is the factor scale of the box axes
!!
!! OUTPUT
!!  per_cond_re= initial in the box
!!
!! SOURCE

function per_cond_re(nb,initial,dim_grid,metric)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nb
 integer,intent(in) :: dim_grid(3)
 integer :: per_cond_re(3,0:nb-1)
 real(dp),intent(in) :: initial(3,0:nb-1)
 real(dp),intent(in) :: metric(1:3)

!Local variables-------------------------------
 integer :: ii,jj

! *********************************************************************
 do jj=0,nb-1
  do ii=1,3
   per_cond_re(ii,jj)=modulo(nint(initial(ii,jj)/metric(ii)),dim_grid(ii))
  end do
 end do


end function per_cond_re
!!***

!!****f* m_per_cond/per_cond_int
!! NAME
!!  per_cond_int
!!
!! FUNCTION
!! Given a 2d-array of integer initial(3,nb), it calulates the values
!! of any of the array in the periodic orthogonal discretized
!! grid begining in 0 and of lengths dim_grid(3)
!!
!!
!! INPUTS
!!  initial(1:3,0:nb-1)=initial point
!!  dim_grid(1:3)= box lengths
!!  nb= dimension
!!  metric=is the factor scale of the box axes
!!
!! OUTPUT
!!  per_cond_int= initial in the box
!!
!! SOURCE

function per_cond_int(nb,initial,dim_grid)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nb
 integer,intent(in) :: dim_grid(3)
 integer :: per_cond_int(3,0:nb-1)
 integer, intent(in) :: initial(3,0:nb-1)

!Local variables-------------------------------
 integer :: ii,jj

! *********************************************************************
 do jj=0,nb-1
  do ii=1,3
   per_cond_int(ii,jj)=modulo(initial(ii,jj),dim_grid(ii))
  end do
 end do


end function per_cond_int
!!***

!!****f* m_per_cond/per_cond_int3
!! NAME
!!  per_cond_int3
!!
!! FUNCTION
!! Given a 2d-array of integer initial(3,nb), it calulates the values
!! of any of the array in the periodic orthogonal discretized
!! grid begining in 0 and of lengths dim_grid(3)
!!
!!
!! INPUTS
!!  initial(1:3,0:nb-1)=initial point
!!  dim_grid(1:3)= box lengths
!!  nb= dimension
!!  metric=is the factor scale of the box axes
!!
!! OUTPUT
!!  per_cond_int3= initial in the box
!!
!! SOURCE

function per_cond_int3(initial,dim_grid)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dim_grid(3)
 integer :: per_cond_int3(3)
 integer,intent(in) :: initial(3)


!Local variables-------------------------------
 integer :: ii

! *********************************************************************

 do ii=1,3
  per_cond_int3(ii) = modulo(initial(ii),dim_grid(ii))
 end do

end function per_cond_int3
!!***

!!****f* m_per_cond/per_dist_int
!! NAME
!!  per_dist_int
!!
!! FUNCTION
!! Given a two 2d-array of integer initA(3,nb),initB(3,nb) in the
!! periodic grid, it calulates the values
!! of any of the distances in the periodic orthogonal discretized
!! grid begining in 0 and of lengths dim_grid(3)
!!
!!
!! INPUTS
!!  initA(1:3,0:nb-1)=initial point
!!  initB(1:3,0:nb-1)=initial point
!!  dim_grid(1:3)= box lengths
!!  nb= dimension
!!
!! OUTPUT
!!  per_dist_int= initial in the box
!!
!! SOURCE

function per_dist_int(nb,initA,initB,dim_grid)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nb
 integer,intent(in) :: dim_grid(3)
 integer :: per_dist_int(3,0:nb-1)
 integer,intent(in) :: initA(3,0:nb-1),initB(3,0:nb-1)


!Local variables-------------------------------
 integer :: ii,jj, bo

! *********************************************************************
 do jj=0,nb-1
  do ii=1,3
   bo = abs(initA(ii,jj)-initB(ii,jj))
   per_dist_int(ii,jj) = minval((/ bo,dim_grid(ii)-bo/))
  end do
 end do


end function per_dist_int
!!***


!!****f* m_per_cond/per_dist_int1
!! NAME
!!  per_dist_int1
!!
!! FUNCTION
!! Given a two scalars of integer initA,initB in the
!! periodic grid, it calulates the values
!! of any of the distances in the periodic orthogonal discretized
!! grid begining in 0 and of lengths dim_grid
!!
!!
!! INPUTS
!!  initA=initial point
!!  initB=initial point
!!  dim_grid= box lengths
!!
!! OUTPUT
!!  per_dist_int1= initial in the box
!!
!! SOURCE

function per_dist_int1(initA,initB,dim_grid)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dim_grid,initA,initB
 integer :: per_dist_int1


!Local variables-------------------------------
 integer :: bo

! *********************************************************************
 bo = abs(initA-initB)
 per_dist_int1 = minval((/ bo,dim_grid-bo/))

end function per_dist_int1
!!***


!!****f* m_per_cond/per_dist_int3
!! NAME
!!  per_dist_int3
!!
!! FUNCTION
!! Given a two 3d-vector of integer initA(3),initB(3) in the
!! periodic grid, it calulates the values
!! of any of the distances in the periodic orthogonal discretized
!! grid begining in 0 and of lengths dim_grid(3)
!!
!!
!! INPUTS
!!  initA(1:3)=initial point
!!  initB(1:3)=initial point
!!  dim_grid(1:3)= box lengths
!!
!! OUTPUT
!!  per_dist_int3= initial in the box
!!
!! SOURCE

function per_dist_int3(initA,initB,dim_grid)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dim_grid(3)
 integer :: per_dist_int3(3)
 integer,intent(in) :: initA(3),initB(3)


!Local variables-------------------------------
 integer :: ii, bo

! *********************************************************************
 do ii=1,3
  bo = abs(initA(ii)-initB(ii))
  per_dist_int3(ii) = minval((/ bo,dim_grid(ii)-bo/))
 end do

end function per_dist_int3
!!***

END MODULE m_per_cond
!!***
