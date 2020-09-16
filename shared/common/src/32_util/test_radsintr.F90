!!****p* ABINIT/check/test_radsintr
!! NAME
!!  test_radsintr
!!
!! FUNCTION
!!  Tests the routine src/32_util/radsintr.F90 that performs
!!  radial Fourier transform.
!!
!! COPYRIGHT
!!  Copyright (C) 2011-2020 ABINIT group (CE)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! NOTES
!!  All required data are generated internally. This program
!!  performs radial FT of a gaussian function and then transforms
!!  it back to real space to check possible differences.
!!
!! PARENTS
!!
!! CHILDREN
!!      radsintr
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program test_radsintr

 use defs_basis
 use m_abicore
 use m_errors

 use m_integrals,     only : radsintr

 implicit none

!Arguments ------------------------------------

!Local variables-------------------------------
!scalars
 integer :: ii,qmesh_size,rmesh_size
 real(dp) :: qstep,rstep,yq1,yqn,yr1,yrn
!character(len=500) :: message ! to be uncommented, if needed
!arrays
 real(dp), allocatable :: funr(:),funr2(:),funq(:),qmesh(:),rmesh(:)

! *************************************************************************

!Parameters for the generation of the grids
 rmesh_size = 150
 rstep = 0.02d0
 qmesh_size = 50
 qstep = 0.04

!Initialization of meshes
 ABI_ALLOCATE(rmesh,(rmesh_size))
 ABI_ALLOCATE(qmesh,(qmesh_size))
 do ii=1,rmesh_size
   rmesh(ii) =rstep*dble(ii-1)
 end do
 do ii=1,qmesh_size
   qmesh(ii) =qstep*dble(ii-1)
 end do

 ABI_ALLOCATE(funr,(rmesh_size))
 ABI_ALLOCATE(funq,(qmesh_size))
 ABI_ALLOCATE(funr2,(rmesh_size))

 do ii=1,rmesh_size
   funr(ii)=exp(-((rmesh(ii)-one)**two)/0.1d0)
 end do

 call radsintr(funr,funq,qmesh_size,rmesh_size,qmesh,rmesh,yq1,yqn)

 write(std_out,*) ch10, 'r  F(r):',ch10
 do ii=1,rmesh_size
   write(std_out,*) rmesh(ii), funr(ii)
 end do

 write(std_out,*) ch10, 'q  F(q):',ch10
 do ii=1,qmesh_size
   write(std_out,*) qmesh(ii), funq(ii)
 end do

 call radsintr(funq,funr2,rmesh_size,qmesh_size,rmesh,qmesh,yr1,yrn)

 write(std_out,*) ch10, 'r  F2(r):',ch10
 do ii=1,rmesh_size
   write(std_out,*) rmesh(ii), funr2(ii)
 end do

 ABI_DEALLOCATE(funq)
 ABI_DEALLOCATE(funr)
 ABI_DEALLOCATE(funr2)
 ABI_DEALLOCATE(rmesh)
 ABI_DEALLOCATE(qmesh)

 end program test_radsintr
!!***
