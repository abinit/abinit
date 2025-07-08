!!****f* ABINIT/m_rot_cg
!! NAME
!! m_rot_cg
!!
!! FUNCTION
!!  Rotate the cg coefficient with the rotation matrix obtained from the
!!  diagonalization of the non-diagonal occupation matrix produced by DMFT.
!!
!! INPUTS
!!
!! OUTPUT
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
  use m_abi_linalg, only : abi_xgemm
  use m_abicore
  use m_errors
  use m_xmpi

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
!! and return diagonalised occupations and associated eigenvectors sorted
!! with descending occupation
!!
!! INPUTS
!!   occ_nd_cpx(nband,nband) = matrix of non diagonal occupations for DMFT
!!   nband = number of bands to be processed
!!
!! OUTPUT
!!   occ_diag(2,nband) = diagonal occupations in the new band space
!!
!! SOURCE
!!
!! TODO /!\ No parallel computing yet !
!! TODO add the possibility of using ScaLAPACK to do computation in parallel

subroutine diag_occ(occ_nd_cpx,nband,occ_diag)

!Arguments ------------------------------------
!scalars
  integer, intent(in) :: nband
!! type(MPI_type),intent(inout) :: mpi_enreg
!! type(dataset_type),intent(in) :: dtset
!! type(paw_dmft_type), intent(in)  :: paw_dmft
!no_abirules
  complex(dpc), intent(inout) :: occ_nd_cpx(nband,nband)
  real(dp), intent(inout) :: occ_diag(nband)
!Local variables-------------------------------
  integer :: info,lwork
  character(len=500) :: message
  real(dp) :: rwork(3*nband-1)
  complex(dpc), allocatable :: work(:)
! *************************************************************************

  DBG_ENTER("COLL")

!! Use the opposite to have zheev orders the eigenvalues by descending order.
!! Afterwards, we multiply by -1 once again.
  occ_nd_cpx(:,:) = -occ_nd_cpx(:,:)

!! Get diagonal occupations and associated base

! Compute the optimal working array size
  ABI_MALLOC(work,(1))
  call zheev('V','U',nband,occ_nd_cpx(:,:),nband,occ_diag(:),work(:),-1,rwork(:),info)
  lwork = int(work(1))
  ABI_FREE(work)

! Compute the eigenvalues (occ_diag) and vectors
  ABI_MALLOC(work,(lwork))

  call zheev('V','U',nband,occ_nd_cpx(:,:),nband,occ_diag(:),work(:),lwork,rwork(:),info)

!! Obtain the true eigenvalues of occupation matrix in descending order
  occ_diag(:) = -occ_diag(:)

  ABI_FREE(work)

  if (info > 0) then
    message = ""
    write(message,"(a,i5)") " something wrong happened with the diagonalization of &
       & the occupation matrix (didn't converge), info=",info
    ABI_ERROR(message)
  else if (info < 0) then
    message = ""
    write(message,"(a,i5)") " something wrong happened with the diagonalization of &
       & the occupation matrix (bad input argument), info=",info
    ABI_ERROR(message)
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
!!   occ_nd(2,nband,nband) = matrix of non diagonal occupations for DMFT
!!   cwavef(2,npw,nband) = Fourier coefficients of wave functions for all bands
!!   npw = number of G vectors computed in this iteration
!!   nband = number of bands to be processed
!!   blocksize = size of the block for the LO.. algorithm
!!               still has to be equal to nband
!!   nspinor = number of spinor components
!!   first_bandc = index of the first correlated band
!!   nbandc = number of correlated bands
!!   dmft_test = for backwards compatibility of some tests
!!
!! OUTPUT
!!   occ_diag(nband) = diagonal occupations in the new band space
!!
!! SIDE EFFECT
!!   cwavef is rotated with the unitary matrix obtained from the diagonalization
!!   of occupations (occ_nd)
!! SOURCE
!!
!! TODO /!\ No parallel computing yet !
!! TODO add the possibility of using ScaLAPACK to do computation in parallel
!! TODO Make the computation of the new wf parallel

subroutine rot_cg(occ_nd,cwavef,npw,nband,blocksize,nspinor,first_bandc,nbandc,occ_diag,dmft_test)

!Arguments ------------------------------------
!scalars
  integer, intent(in) :: blocksize,dmft_test,first_bandc,nband,nbandc,npw,nspinor
!! type(MPI_type),intent(inout) :: mpi_enreg
!! type(dataset_type),intent(in) :: dtset
!! type(paw_dmft_type), intent(in)  :: band_in
!no_abirules
  real(dp), intent(in) :: occ_nd(2,blocksize,blocksize)
  real(dp), intent(inout) :: occ_diag(blocksize)
  real(dp), intent(inout) :: cwavef(2,npw,blocksize,nspinor)
!Local variables-------------------------------
!scalars
  integer :: ispinor,n
  character(len=500) :: message
!arrays
  real(dp), allocatable :: occ_diag_red(:)
  complex(dpc), allocatable :: mat_tmp(:,:),mat_tmp2(:,:),occ_nd_cpx(:,:)
  !complex(kind=dpc) :: cwavef_rot_g(nbandc, nspinor)
! *************************************************************************

  DBG_ENTER("COLL")

  if (nband /= blocksize) then
    message = " DMFT in KGB cannot be used with multiple blocks yet. Make sure that bandpp*npband = nband."
    ABI_ERROR(message)
  end if

!! Initialization

  ABI_MALLOC(mat_tmp,(npw,nbandc))
  ABI_MALLOC(mat_tmp2,(npw,nbandc))
  ABI_MALLOC(occ_diag_red,(nbandc))
  ABI_MALLOC(occ_nd_cpx,(nbandc,nbandc))

  occ_nd_cpx(:,:) = cmplx(occ_nd(1,first_bandc:first_bandc+nbandc-1,first_bandc:first_bandc+nbandc-1), &
                        & occ_nd(2,first_bandc:first_bandc+nbandc-1,first_bandc:first_bandc+nbandc-1),kind=dp)

!! Get diagonal occupations and associated base

  call diag_occ(occ_nd_cpx(:,:),nbandc,occ_diag_red(:))

  do n=1,nband
    if (n < first_bandc .or. n >= first_bandc+nbandc) then
      occ_diag(n) = occ_nd(1,n,n)
    else
      occ_diag(n) = occ_diag_red(n-first_bandc+1)
    end if
  end do ! n

!! Compute the corresponding wave functions if nothing wrong happened
  ! $c^{rot}_{n,k}(g) =  \sum_{n'} [\bar{f_{n',n}} * c_{n',k}(g)]$

  if (dmft_test == 1) occ_nd_cpx(:,:) = conjg(occ_nd_cpx(:,:))

  do ispinor=1,nspinor
    mat_tmp(:,:) = cmplx(cwavef(1,1:npw,first_bandc:first_bandc+nbandc-1,ispinor), &
                       & cwavef(2,1:npw,first_bandc:first_bandc+nbandc-1,ispinor),kind=dp)
    call abi_xgemm("n","n",npw,nbandc,nbandc,cone,mat_tmp(:,:),npw, &
                 & occ_nd_cpx(:,:),nbandc,czero,mat_tmp2(:,:),npw)
    cwavef(1,1:npw,first_bandc:first_bandc+nbandc-1,ispinor) = dble(mat_tmp2(:,:))
    cwavef(2,1:npw,first_bandc:first_bandc+nbandc-1,ispinor) = aimag(mat_tmp2(:,:))
  end do ! ispinor

  !do ig=1,npw
  !  cwavef_rot_g(:,:) = czero
  !  do n=1,nbandc
  !    do np=1,nbandc
  !      cwavef_rot_g(n,:) = cwavef_rot_g(n,:) + occ_nd_cpx(np, n) * &
!&                           cmplx(cwavef(1,ig,np+first_bandc-1,:), cwavef(2,ig,np+first_bandc-1,:), kind=dpc)
 !     end do
 !   end do
 !   cwavef(1,ig,first_bandc:first_bandc+nbandc-1,:) = dreal(cwavef_rot_g)
 !   cwavef(2,ig,first_bandc:first_bandc+nbandc-1,:) = dimag(cwavef_rot_g)
 ! end do

  ABI_FREE(mat_tmp)
  ABI_FREE(mat_tmp2)
  ABI_FREE(occ_diag_red)
  ABI_FREE(occ_nd_cpx)

  DBG_EXIT("COLL")

 end subroutine rot_cg
!!***

end module m_rot_cg
!!***
