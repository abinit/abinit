!!****m* ABINIT/m_frskerker2
!! NAME
!! m_frskerker2
!!
!! FUNCTION
!! provide the ability to compute the
!! penalty function and its first derivative associated
!! with some residuals and a real space dielectric function
!!
!! COPYRIGHT
!! Copyright (C) 1998-2020 ABINIT group (DCA, XG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! NOTES
!! this is neither a function nor a subroutine. This is a module
!! It is made of two functions and one init subroutine
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

module m_frskerker2

  use defs_basis
  use m_abicore
  use m_dtset

  use defs_abitypes, only : MPI_type
  use m_spacepar, only : laplacian
  use m_numeric_tools, only : dotproduct

  implicit none

  !! common variables copied from input
  integer,save,private                  :: nfft,nspden,ngfft(18)
  real(dp),save,allocatable,private    :: deltaW(:,:),mat(:,:),rdielng(:)
  real(dp),save,private                :: gprimd(3,3)
  type(dataset_type),pointer,save,private  :: dtset_ptr
  type(MPI_type),pointer,save,private  :: mpi_enreg_ptr
  !! common variables computed
  logical,save,private :: ok=.false.
! *************************************************************************

contains
!!***

!!****f* m_frskerker2/frskerker2__init
!! NAME
!! frskerker2__init
!!
!! FUNCTION
!! initialisation subroutine
!! Copy every variables required for the energy calculation
!! Allocate the required memory
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_prcref
!!
!! CHILDREN
!!      laplacian
!!
!! SOURCE

subroutine frskerker2__init(dtset_in,mpi_enreg_in,nfft_in,ngfft_in,nspden_in,rdielng_in,deltaW_in,gprimd_in,mat_in )

!Arguments ------------------------------------
 type(dataset_type),target,intent(in) :: dtset_in
 integer,intent(in)  :: nfft_in,ngfft_in(18),nspden_in
 real(dp),intent(in) :: deltaW_in(nfft_in,nspden_in),mat_in(nfft_in,nspden_in)
 real(dp),intent(in) :: rdielng_in(nfft_in)
 real(dp),intent(in)  :: gprimd_in(3,3)
 type(MPI_type),target,intent(in)  :: mpi_enreg_in

! *************************************************************************
! !allocation and data transfer
! !Thought it would have been more logical to use the privates intrinsic of the module as
! !input variables it seems that it is not possible...
  if(.not.ok) then
   dtset_ptr => dtset_in
   mpi_enreg_ptr => mpi_enreg_in
   nspden=nspden_in
   ngfft=ngfft_in
   nfft=nfft_in
   ABI_ALLOCATE(deltaW,(size(deltaW_in,1),size(deltaW_in,2)))
   ABI_ALLOCATE(mat,(size(mat_in,1),size(mat_in,2)))
   ABI_ALLOCATE(rdielng,(size(rdielng_in)))
   deltaW=deltaW_in
   rdielng=rdielng_in
   mat=mat_in
   gprimd=gprimd_in
   ok = .true.
  end if

 end subroutine frskerker2__init
!!***

!!****f* m_frskerker2/frskerker2__end
!! NAME
!! frskerker2__end
!!
!! FUNCTION
!! ending subroutine
!! deallocate memory areas
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_prcref
!!
!! CHILDREN
!!      laplacian
!!
!! SOURCE

  subroutine frskerker2__end()

! *************************************************************************
  if(ok) then
!  ! set ok to false which prevent using the pf and dpf
   ok = .false.
!  ! free memory
   ABI_DEALLOCATE(deltaW)
   ABI_DEALLOCATE(mat)
   ABI_DEALLOCATE(rdielng)
  end if

 end subroutine frskerker2__end
!!***

!!****f* m_frskerker2/frskerker2__newvres2
!! NAME
!! frskerker2__newvres2
!!
!! FUNCTION
!! affectation subroutine
!! do the required renormalisation when providing a new value for
!! the density after application of the gradient
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      laplacian
!!
!! SOURCE

subroutine frskerker2__newvres2(nv1,nv2,x, grad, vrespc)

!Arguments ------------------------------------
 integer,intent(in)    :: nv1,nv2
 real(dp),intent(in)   :: x
 real(dp),intent(inout)::grad(nv1,nv2)
 real(dp),intent(inout)::vrespc(nv1,nv2)

! *************************************************************************
  grad(:,:)=x*grad(:,:)
  vrespc(:,:)=vrespc(:,:)+grad(:,:)

 end subroutine frskerker2__newvres2
!!***

!!****f* m_frskerker2/frskerker2__pf
!! NAME
!! frskerker2__pf
!!
!! FUNCTION
!! penalty function associated with the preconditionned residuals
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

  function frskerker2__pf(nv1,nv2,vrespc)

!Arguments ------------------------------------
 integer,intent(in)    :: nv1,nv2
 real(dp),intent(in) ::vrespc(nv1,nv2)
 real(dp)            ::frskerker2__pf

!Local variables-------------------------------
 real(dp)            :: buffer1(nv1,nv2),buffer2(nv1,nv2)
 integer             :: ispden
! *************************************************************************

  if(ok) then
   buffer1=vrespc
   call laplacian(gprimd,mpi_enreg_ptr,nfft,nspden,ngfft,rdfuncr=buffer1,laplacerdfuncr=buffer2)
   do ispden=1,nspden
    buffer2(:,ispden)=(vrespc(:,ispden)-((rdielng(:))**2)*buffer2(:,ispden))  &
&    *half  -  deltaW(:,ispden)
   end do
!  pf_rscgres=dotproduct(vrespc,buffer2)*half-dotproduct(vrespc,deltaW)
   frskerker2__pf=dotproduct(nv1,nv2,vrespc,buffer2)
  else
   frskerker2__pf=zero
  end if
 end function frskerker2__pf
!!***

!!****f* m_frskerker2/frskerker2__dpf
!! NAME
!! frskerker2__dpf
!!
!! FUNCTION
!! derivative of the penalty function
!! actually not the derivative but something allowing minimization
!! at constant density
!! formula from the work of rackowski,canning and wang
!! H*phi - int(phi**2H d3r)phi
!! that is the simple projection of the change on a direction
!! normal to the density changes
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

function frskerker2__dpf(nv1,nv2,vrespc)

!Arguments ------------------------------------
 integer,intent(in) :: nv1,nv2
 real(dp),intent(in)::vrespc(nv1,nv2)
 real(dp)           :: frskerker2__dpf(nv1,nv2)

!Local variables-------------------------------
 real(dp):: buffer1(nv1,nv2),buffer2(nv1,nv2)
 integer :: ispden

! *************************************************************************

  if(ok) then
   buffer1=vrespc
   call laplacian(gprimd,mpi_enreg_ptr,nfft,nspden,ngfft,rdfuncr=buffer1,laplacerdfuncr=buffer2)
   do ispden=1,nspden
    frskerker2__dpf(:,ispden)= vrespc(:,ispden)-deltaW(:,ispden)-((rdielng(:))**2)*buffer2(:,ispden)
   end do
  else
   frskerker2__dpf = zero
  end if

end function frskerker2__dpf
!!***

end module m_frskerker2
!!***
