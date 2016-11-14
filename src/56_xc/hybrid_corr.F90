!{\src2tex{textfont=tt}}
!!****f* ABINIT/hybrid_corr
!! NAME
!! hybrid_corr
!!
!! FUNCTION
!! Compute the correction to the XC potetential and
!! energy due to the hybridation
!!
!! COPYRIGHT
!!  Copyright (C) 2015-2016 ABINIT group (FA)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  ixc = choice of exchange-correlation functional.
!!  mpi_enreg=information about MPI parallelization
!!  nfft = number of fft grid points.
!!  ngfft(1:3) = integer fft box dimensions, see getng for ngfft(4:8).
!!  nspden = number of spin-density components.
!!  rhor(nfft,nspden) = electron density in real space in electrons/bohr**3
!!   (total in first half and spin-up in second half if nspden = 2).
!!  rprimd(3,3) = dimensional primitive translations for real space in Bohr.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  enxc = exchange correlation energy corrected for hybrid functionals
!!  vxc = exchange correlation potential corrected for hybrid functionals
!!
!! WARNINGS
!! Current restrictions are:
!!  a - Spin-polarized case not tested.
!!
!! PARENTS
!!
!! CHILDREN
!!      dtset_copy,dtset_free,libxc_functionals_end,libxc_functionals_init
!!      rhohxc
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine hybrid_corr(dtset,ixc,nkxc,mpi_enreg,nfft,ngfft,nspden,rhor,rprimd,hybrid_mixing,vxc,enxc)

 use defs_basis
 use m_profiling_abi
 use m_errors
 use defs_abitypes, only : MPI_type, dataset_type
 use m_dtset,       only : dtset_copy, dtset_free

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hybrid_corr'
 use interfaces_56_xc, except_this_one => hybrid_corr
!End of the abilint section

 implicit none

!Arguments -------------------------------------------------------------
!scalars
 integer,intent(in) :: ixc,nfft,nspden,nkxc
 real(dp),intent(in) :: hybrid_mixing
 real(dp),intent(inout) :: enxc
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rhor(nfft,nspden),rprimd(3,3)
 real(dp),intent(inout) :: vxc(nfft,nspden)

!Local variables -------------------------------------------------------
!No improved xc quadrature.
!No core correction.
!Dummy here.
!For debugging purposes (see tests below):
!integer :: i1,i2,i3,k1,n1,n2,n3
!real(dp) :: kx,rho,rhomax,ftest
!scalars
 character(len=500) :: message
 type(dataset_type) :: dtLocal
!arrays
 real(dp) :: dum(0)
 real(dp) :: enxc_corr
 real(dp),allocatable :: kxcr(:,:)
 real(dp),allocatable :: xccc3d(:),vxc_corr(:,:)

!***********************************************************************
!For debugging purposes (see tests below):

!ftest(i1,n1,k1) = 0._dp+1._dp*cos(k1*two_pi*float(i1)/float(n1))

!***********************************************************************

!Check input parameters.

 if (nspden > 2) then
   message = ' kxc_alda does not work yet for nspden > 2.'
   MSG_ERROR(message)
 end if

!Copy the input variables from the current dataset to a temporary one
!to tune some parameters
 call dtset_copy(dtLocal, dtset)
 dtLocal%intxc = 0
 dtLocal%ixc   = -101

 ! Reinitialize the libxc module with the overriden values
 if (dtset%ixc<0) then
   call libxc_functionals_end()
 end if
 if (dtLocal%ixc<0) then
   call libxc_functionals_init(dtLocal%ixc,dtLocal%nspden)
 end if

 call rhohxc(dtLocal,enxc_corr,dum,0,kxcr,mpi_enreg,nfft,dum,dum,0,dum,0,nkxc,0,nspden,n3xccc,&
& 0,dum,rhor,rprimd,strsxc,1,dum,vxc_corr,dum,xccc3d)

 vxc(:,:) = vxc(:,:) + hybrid_mixing*vxc_corr(:,:)
 enxc = enxc + hybrid_mixing*enxc_corr

 call rhohxc(dtLocal,enxc_corr,dum,0,kxcr,mpi_enreg,dum,dum,dum,0,dum,0,nkxc,0,nspden,0,&
& 0,dum,rhor,rprimd,strsxc,1,dum,vxc_corr,vxcavg,0)

 vxc(:,:) = vxc(:,:) - hybrid_mixing*vxc_corr(:,:)
 enxc = enxc - hybrid_mixing*enxc_corr


! Revert libxc module to the original settings
 if (dtLocal%ixc<0) then
   call libxc_functionals_end()
 end if
 if (dtset%ixc<0) then
   call libxc_functionals_init(dtset%ixc,dtset%nspden)
 end if

!Free memory.
 call dtset_free(dtLocal)
 ABI_FREE(kxcr)
 ABI_FREE(xccc3d)

end subroutine hybrid_corr
!!***

