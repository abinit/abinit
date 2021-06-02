! Experimental tight binding solver
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tight_binding
  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi

  use m_multibinit_dataset, only: multibinit_dtset_type

  implicit none
  private

  type, public :: tight_binding_t
     integer :: nrpt, nkpt
     integer, allocatable :: rpts(:,:)
     integer :: nbasis
     integer :: nspin
     logical :: spinor
     real(dp) :: cell(3,3)
     real(dp), allocatable :: xred(:,:)
     real(dp), allocatable :: T(:,:,:)

   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: set_hopping
     procedure :: set_kmesh
     procedure :: gen_hamk
     procedure :: solve_k

  end type tight_binding_t

contains
  subroutine initialize(self, nbasis, nspin, spinor, cell, xred, nrpt)
    class(tight_binding_t), intent(inout) :: self
    integer, intent(in):: nbasis, nspin, nrpt
    logical, intent(in):: spinor
    real(dp), intent(in) :: cell(3,3)
    real(dp), intent(in) :: xred(3, nbasis)

    self%nbasis=nbasis
    self%nspin=nspin
    self%spinor=spinor

    self%nrpt=nrpt

    ABI_MALLOC(self%xred,(3, nbasis))
    ABI_MALLOC(self%rpts, (3, nrpt))

    self%cell(:,:)=cell(:,:)
    self%xred(:,:)=xred(:,:)

  end subroutine initialize

  subroutine finalize(self)
    class(tight_binding_t), intent(inout) :: self
    if(allocated(self%xred)) then
       ABI_FREE(self%xred)
    end if

    if(allocated(self%rpts)) then
       ABI_FREE(self%rpts)
    end if
  end subroutine finalize

  subroutine set_hopping(self, T)
    class(tight_binding_t), intent(inout) :: self
    real(dp), intent(in) :: T(:,:,:)
    self%T(:,:,:)=T(:,:,:)
  end subroutine set_hopping

  subroutine gen_hamk(self, kpt, ham)
    class(tight_binding_t), intent(inout) :: self
    real(dp), intent(in) :: kpt(3)
    complex(dp), intent(inout) :: ham(self%nbasis,self%nbasis)
    complex(dp) :: phase
    integer :: i

    ham(:,:)=0.0
    do i=1, self%nrpt
       phase= exp(2.0_dp* j_dpc * pi* dot_product(kpt, self%rpts(:,i)))
       ham(:,:)=ham(:,:)+self%T(:,:, i)*phase
    end do
  end subroutine gen_hamk

  subroutine solve_k(self, kpt, evals, evecs)
    class(tight_binding_t), intent(inout) :: self
    real(dp), intent(in) :: kpt(3)
    real(dp), intent(inout) :: evals(self%nbasis)
    complex(dp), intent(inout) :: evecs(self%nbasis,self%nbasis)
    complex(dp)  :: work(3*self%nbasis-2)
    integer ::  info, lwork
    real(dp) ::  Rwork(3*self%nbasis-2)
    external ZHEEV

    call self%gen_hamk(kpt, evecs)

    lwork=3*self%nbasis-2
    call ZHEEV('V', 'U', self%nbasis, evecs, self%nbasis, evals, work, lwork, rwork, info)

  end subroutine solve_k

  subroutine set_kmesh(self, kmesh)
    class(tight_binding_t), intent(inout) :: self
    integer, intent(in) :: kmesh(3)
    ABI_UNUSED_A(self)
    ABI_UNUSED_A(kmesh)
  end subroutine set_kmesh


end module m_tight_binding
