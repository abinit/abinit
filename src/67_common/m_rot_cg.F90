!{\src2tex{textfont=tt}}
!!****f* ABINIT/m_rot_cg
!! NAME
!! m_rot_cg
!!
!! FUNCTION
!!  Rotate the cg coefficient with the rotation matrix obtained from the
!!  diagonalisation of the non-diagonal occupation matrix produced by DMFT.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      mkrho
!!
!! CHILDREN
!!
!! SOURCE
!!
!! TODO /!\ No parallel computing yet !

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_rot_cg

  use defs_basis
  use m_xmpi
  use m_errors
  use m_abicore

  implicit none

  private

  public :: rot_cg

  contains
!!***

!!****f* ABINIT/diag_occ
!! NAME
!! diag_occ
!!
!! FUNCTION
!! Use for DMFT in KGB parallelisation. Diagonalise the occupation matrix
!! and return diagonalised occupations and associated eigen vectors sorted
!! with descending occupation
!!
!! INPUTS
!!   occ_nd_cpx(nband, nband) = matrix of non diagonal occupations for DMFT
!!   nband = number of band to be processed
!!
!!
!! OUTPUT
!!   occ_diag(2, nband) = diagonal occupation in the new band space
!!   cwavef_rot(2, npw, nband) = fourier coefficient of wave functions for all bands rotated in the band space
!!
!! PARENTS
!!      mkrho
!!
!! CHILDREN
!!
!! SOURCE
!!
!! TODO /!\ No parallel computing yet !
!! TODO add the possibility of using ScaLAPACK to do computation in parallel

subroutine diag_occ(occ_nd_cpx, nband, occ_diag)

!Arguments ------------------------------------
!scalars
  integer,intent(in) :: nband
!! type(MPI_type),intent(inout) :: mpi_enreg
!! type(dataset_type),intent(in) :: dtset
!! type(paw_dmft_type), intent(in)  :: paw_dmft
!no_abirules
  complex(kind=dpc), intent(inout) :: occ_nd_cpx(nband, nband)
  real(kind=dp), intent(out) :: occ_diag(nband)

!Local variables-------------------------------

!scalars
  integer :: info, lwork
  character(len=500) :: message

!arrays
  real(kind=dp) :: rwork(3*nband-1)
  complex(kind=dpc), allocatable :: work(:)

! *************************************************************************

  DBG_ENTER("COLL")

!! Initialisation
  rwork = zero
  info = 0

!! use the opposite to have zheev orders the eigenvalues descending (once the
!! opposite have been taken again)
  occ_nd_cpx = -occ_nd_cpx

!! Get diagonal occupations and associeted base

! Compute the optimal working array size
  ABI_ALLOCATE(work,(1))
  work = czero
  call zheev('V', 'U', nband, occ_nd_cpx, nband, occ_diag, work, -1, rwork, info)
  lwork = work(1)
  ABI_DEALLOCATE(work)

! Compute the eigenvalues (occ_diag) and vectors
  ABI_ALLOCATE(work,(lwork))
  work = czero

  call zheev('V', 'U', nband, occ_nd_cpx, nband, occ_diag, work, lwork, rwork, info)

!! obtain the true eigen values of occupation matrix in descending order
  occ_diag = -occ_diag

  ABI_DEALLOCATE(work)

  if (info > 0) then
    message=""
    write(message, "(a,i5)") " something wrong happened with the diagonalisation of &
&the occupation matrix (did't converge), info=",info
    MSG_ERROR(message)
  else if (info < 0) then
    message=""
    write(message, "(a,i5)") " something wrong happened with the diagonalisation of &
&the occupation matrix (bad input argument), info=",info
    MSG_ERROR(message)
  end if

  DBG_EXIT("COLL")

end subroutine diag_occ
!!***

!!****f* ABINIT/rot_cg
!! NAME
!! rot_cg
!!
!! FUNCTION
!! Use for DMFT in KGB parallelisation. Diagonalise the occupation matrix
!! and use the resulting base to represent the wave functions.
!!
!! INPUTS
!!   occ_nd(2, nband, nband) = matrix of non diagonal occupations for DMFT
!!   cwavef(2, npw, nband) = fourier coefficient of wave functions for all bands
!!   npw = number of G vectors computed in this iteration
!!   nband = number of band to be processed
!!   blocksize = size of the block for th LO.. algorithm
!!               still have to be equal to nband
!!   nspinor = number of spinor components
!!   first_bandc = index of the first correlated band
!!   nbandc = number of bands correlated
!!
!!
!! OUTPUT
!!   occ_diag(2, nband) = diagonal occupation in the new band space
!!
!! SIDE EFFECT
!!   cwavef is rotated with the unitary matrix obtained from the diagonalisation
!!   of occupations (occ_nd)
!! PARENTS
!!      mkrho
!!
!! CHILDREN
!!
!! SOURCE
!!
!! TODO /!\ No parallel computing yet !
!! TODO add the possibility of using ScaLAPACK to do computation in parallel
!! TODO Make the computation of the new wf parallel

subroutine rot_cg(occ_nd, cwavef, npw, nband, blocksize, nspinor, first_bandc, nbandc, occ_diag)

!Arguments ------------------------------------
!scalars
  integer,intent(in) :: npw, nband, blocksize, nspinor, first_bandc, nbandc
!! type(MPI_type),intent(inout) :: mpi_enreg
!! type(dataset_type),intent(in) :: dtset
!! type(paw_dmft_type), intent(in)  :: band_in
!no_abirules
  real(kind=dp), intent(in) :: occ_nd(2, blocksize, blocksize)
  real(kind=dp), intent(inout) :: cwavef(2, npw, blocksize, nspinor)
  real(kind=dp), intent(out) :: occ_diag(blocksize)

!Local variables-------------------------------

!scalars
  integer :: n, np, ig
  character(len=500) :: message

!arrays
  real(kind=dp) :: occ_diag_red(nbandc)
  complex(kind=dpc) :: occ_nd_cpx(nbandc, nbandc)
  complex(kind=dpc) :: cwavef_rot_g(nbandc, nspinor)

! *************************************************************************

  DBG_ENTER("COLL")

  if(nband /= blocksize) then
    message = " DMFT in KGB cannot be used with multiple blocks yet. Make sure that bandpp*npband = nband."
    MSG_ERROR(message)
  end if

!! Initialisation

  do n=1,nbandc
    do np=1,nbandc
      occ_nd_cpx(n,np) = cmplx(occ_nd(1,n+first_bandc-1,np+first_bandc-1), occ_nd(2,n+first_bandc-1,np+first_bandc-1), kind=dpc)
    end do
  end do

!! Get diagonal occupations and associeted base

  call diag_occ(occ_nd_cpx, nbandc, occ_diag_red)

  do n=1,nband
    if (n < first_bandc .or. n >= first_bandc+nbandc) then
      occ_diag(n) = occ_nd(1, n, n)
    else
      occ_diag(n) = occ_diag_red(n-first_bandc+1)
    end if
  end do

!! Compute the corresponding wave functions if nothing wrong happened
  ! $c^{rot}_{n,k}(g) =  \sum_{n'} [\bar{f_{n',n}} * c_{n',k}(g)]$
  do ig=1,npw
    cwavef_rot_g(:,:) = czero
    do n=1,nbandc
      do np=1,nbandc
        cwavef_rot_g(n,:) = cwavef_rot_g(n,:) + occ_nd_cpx(np, n) * &
&                           cmplx(cwavef(1,ig,np+first_bandc-1,:), cwavef(2,ig,np+first_bandc-1,:), kind=dpc)
      end do
    end do
    cwavef(1,ig,first_bandc:first_bandc+nbandc-1,:) = dreal(cwavef_rot_g)
    cwavef(2,ig,first_bandc:first_bandc+nbandc-1,:) = dimag(cwavef_rot_g)
  end do

  DBG_EXIT("COLL")

 end subroutine rot_cg

end module m_rot_cg
!!***
