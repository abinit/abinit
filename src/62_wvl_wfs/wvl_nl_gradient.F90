!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_nl_gradient
!! NAME
!! wvl_nl_gradient
!!
!! FUNCTION
!! Compute the non local part of the wavefunction gradient.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2016 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      afterscfloop,vtorho
!!
!! CHILDREN
!!      nonlocal_forces,wrtout,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_nl_gradient(grnl, mpi_enreg, natom, rprimd, wvl, xcart)

 use defs_basis
 use defs_abitypes
 use defs_wvltypes
 use m_profiling_abi
 use m_errors
 use m_xmpi

#if defined HAVE_BIGDFT
 use BigDFT_API, only: nonlocal_forces
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_nl_gradient'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom
 type(MPI_type),intent(in) :: mpi_enreg
 type(wvl_data),intent(inout) :: wvl
!arrays
 real(dp),intent(in) :: xcart(3,natom),rprimd(3,3)
 real(dp),intent(inout) :: grnl(3,natom)

!Local variables-------------------------------
#if defined HAVE_BIGDFT
!scalars
 integer :: ia,ierr,igeo,me,nproc,spaceComm
 character(len=500) :: message
!arrays
 real(dp),allocatable :: gxyz(:,:)
 real(dp)::strtens(6,4)
#endif

! *************************************************************************

#if defined HAVE_BIGDFT

!Compute forces
 write(message, '(a,a)' ) ' wvl_nl_gradient(): compute non-local part to gradient.'
 call wrtout(std_out,message,'COLL')

!Nullify output arrays.
 grnl(:, :) = zero
 strtens(:,:)=zero
 
 ABI_ALLOCATE(gxyz,(3, natom))
 gxyz(:,:) = zero

!Add the nonlocal part of the forces to grtn (BigDFT routine)
 spaceComm=mpi_enreg%comm_wvl
 me=xmpi_comm_rank(spaceComm)
 nproc=xmpi_comm_size(spaceComm)
 call nonlocal_forces(wvl%descr%Glr, &
& wvl%descr%h(1), wvl%descr%h(2), wvl%descr%h(3), wvl%descr%atoms, &
& xcart, wvl%wfs%ks%orbs, wvl%projectors%nlpsp, wvl%wfs%ks%Lzd%Glr%wfd, &
& wvl%wfs%ks%psi, gxyz, .true.,strtens(1,2), &
& proj_G=wvl%projectors%G,paw=wvl%descr%paw)

 if (nproc > 1) then
   call xmpi_sum(gxyz, spaceComm, ierr)
 end if

!Forces should be in reduced coordinates.
 do ia = 1, natom, 1
   do igeo = 1, 3, 1
     grnl(igeo, ia) = - rprimd(1, igeo) * gxyz(1, ia) - &
&                       rprimd(2, igeo) * gxyz(2, ia) - &
&                       rprimd(3, igeo) * gxyz(3, ia)
   end do
 end do
 ABI_DEALLOCATE(gxyz)

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) natom,mpi_enreg%nproc,wvl%wfs%ks,xcart(1,1),rprimd(1,1),grnl(1,1)
#endif

end subroutine wvl_nl_gradient
!!***
