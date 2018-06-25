!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_setngfft
!! NAME
!! wvl_setngfft
!!
!! FUNCTION
!! When wavelets are used, the FFT grid is used to store potentials and
!! density. The size of the grid takes into account the two resolution in wavelet
!! description and also the distribution over processor in the parallel case.
!!
!! The FFT grid is not in strict terms an FFT grid but rather a real space grid.
!! Its dimensions are not directly compatible with FFTs. This is not relevant
!! when using the wavelet part of the code and in the Poisson solver the arrays
!! are extended to match FFT dimensions internally. But for other parts of the
!! code, this must be taken into account.
!!
!! see doc/variables/vargs.html#ngfft for details about ngfft
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SIDE EFFECTS
!!  mpi_enreg=informations about MPI parallelization (description of the
!!            density and potentials scatterring is allocated and updated).
!!  dtset <type(dataset_type)>=the FFT grid is changed.
!!
!! PARENTS
!!      gstate,wvl_wfsinp_reformat
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_setngfft(me_wvl, mgfft, nfft, ngfft, nproc_wvl, n1i, n2i, n3i,n3d)

 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_setngfft'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(out) :: mgfft, nfft
 integer, intent(in)  :: n1i, n2i, n3i,n3d, nproc_wvl, me_wvl
!arrays
 integer, intent(out) :: ngfft(18)

!Local variables-------------------------------
!scalars
#if defined HAVE_BIGDFT
 character(len=500) :: message
#endif

! *************************************************************************

#if defined HAVE_BIGDFT
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_setngfft : Changing the FFT grid definition.'
 call wrtout(std_out,message,'COLL')

!Change nfft and ngfft
!Now ngfft will use the density definition (since the potential size
!is always smaller than the density one). ????
 ngfft(1) = n1i
 ngfft(2) = n2i
 ngfft(3) = n3i

 nfft = n1i*n2i*n3d
!Set up fft array dimensions ngfft(4,5,6) to avoid cache conflicts
!Code paste from getng()
 ngfft(4) = 2 * (ngfft(1) / 2) + 1
 ngfft(5) = 2 * (ngfft(2) / 2) + 1
 ngfft(6) = ngfft(3)
 if (nproc_wvl == 0) then
   ngfft(9)  = 0    ! paral_fft
   ngfft(10) = 1    ! nproc_fft
   ngfft(11) = 0    ! me_fft
   ngfft(12) = 0    ! n2proc
   ngfft(13) = 0    ! n3proc
 else
   ngfft(9)  = 1    ! paral_fft
   ngfft(10) = nproc_wvl
   ngfft(11) = me_wvl
   ngfft(12) = ngfft(2)
   ngfft(13) = n3d
 end if

 write(message, '(a,3I12)' ) &
& '  | ngfft(1:3) is now:    ', ngfft(1:3)
 call wrtout(std_out,message,'COLL')
 write(message, '(a,3I12)' ) &
& '  | ngfft(4:6) is now:    ', ngfft(4:6)
 call wrtout(std_out,message,'COLL')

!Set mgfft
 mgfft= max(ngfft(1), ngfft(2), ngfft(3))

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) mgfft,nfft,n1i,n2i,n3i,n3d,nproc_wvl,me_wvl,ngfft(1)
#endif

end subroutine wvl_setngfft
!!***
