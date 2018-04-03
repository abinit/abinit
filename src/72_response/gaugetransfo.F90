!{\src2tex{textfont=tt}}
!!****f* ABINIT/gaugetransfo
!! NAME
!! gaugetransfo
!!
!! FUNCTION
!! This routine allows the passage from the parallel-transport gauge
!! to the diagonal gauge for the first-order wavefunctions
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group (MVeithen)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg_k(2,mpw*nspinor*mband*nsppol)=planewave coefficients of wavefunctions
!!                                   for a particular k point.
!!  cwavef(2,npw1_k*nspinor)=first order wavefunction for a particular k point
!!                           in the parallel gauge
!!  eig_k(mband*nsppol)=GS eigenvalues at k (hartree)
!!  eig1_k(2*nsppol*mband**2)=matrix of first-order eigenvalues (hartree)
!!  iband=band index of the 1WF for which the transformation has to be applied
!!  mband=maximum number of bands
!!  nband_k=number of bands for this k point
!!  npw_k=maximum dimensioned size of npw or wfs at k
!!  npw1_k=number of plane waves at this k+q point
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ_k(nband_k)=occupation number for each band (usually 2) for each k
!!
!! OUTPUT
!!  cwavef_d(2,npw1_k*nspinor)=first order wavefunction for a particular k point
!!                             in the diagonal gauge
!!
!! PARENTS
!!      dfpt_nstwf
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine gaugetransfo(cg_k,cwavef,cwavef_d,eig_k,eig1_k,iband,nband_k, &
&                      mband,npw_k,npw1_k,nspinor,nsppol,occ_k)


 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gaugetransfo'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iband,mband,nband_k,npw1_k,npw_k,nspinor,nsppol
!arrays
 real(dp),intent(in) :: cg_k(2,npw_k*nspinor*nband_k),cwavef(2,npw1_k*nspinor)
 real(dp),intent(in) :: eig1_k(2*nsppol*mband**2),eig_k(mband*nsppol)
 real(dp),intent(in) :: occ_k(nband_k)
 real(dp),intent(out) :: cwavef_d(2,npw1_k*nspinor)

!Local variables-------------------------------
!tolerance for non degenerated levels
!scalars
 integer :: jband
 real(dp),parameter :: etol=1.0d-3
!arrays
 real(dp) :: cwave0(2,npw1_k*nspinor),eig1(2)

! *********************************************************************

!DEBUG
!write(100,*) 'gaugetransfo: ',iband
!ENDDEBUG

 cwavef_d(:,:) = cwavef(:,:)

 do jband = 1,nband_k !loop over bands

   if ((abs(eig_k(iband)-eig_k(jband)) > etol).and.(abs(occ_k(jband)) > tol8 )) then

     cwave0(:,:) = cg_k(:,1+(jband-1)*npw_k*nspinor:jband*npw_k*nspinor)

     eig1(1) = eig1_k(2*jband-1+(iband-1)*2*nband_k)
     eig1(2) = eig1_k(2*jband +(iband-1)*2*nband_k)

     cwavef_d(1,:)=cwavef_d(1,:) &
&     - (eig1(1)*cwave0(1,:)-eig1(2)*cwave0(2,:))/(eig_k(jband)-eig_k(iband))
     cwavef_d(2,:)=cwavef_d(2,:) &
&     - (eig1(1)*cwave0(2,:)+eig1(2)*cwave0(1,:))/(eig_k(jband)-eig_k(iband))

   end if

 end do    !loop over bands

end subroutine gaugetransfo
!!***
