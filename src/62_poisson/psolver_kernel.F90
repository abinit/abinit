!{\src2tex{textfont=tt}}
!!****f* ABINIT/psolver_kernel
!! NAME
!! psolver_kernel
!!
!! FUNCTION
!! Build, get or free the kernel matrix used by the Poisson solver to compute the
!! the convolution between 1/r and rho. The kernel is a saved variable. If
!! this routine is called for building while a kernel already exists, it is not
!! recomputed if all parameters (grid step and data size) are unchanged. Otherwise
!! the kernel is freed and recompute again. The build action has a returned variable
!! which is a pointer on the kernel. The get action also returns the kernel, or
!! NULL if none has been associated.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR, TRangel).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  iaction=0 to free all kernel allocated array,
!!          1 to compute the kernel (parallel case),
!!          2 to get it (parallel case),
!!          3 to compute the kernel (sequential),
!!          4 to get the sequential kernel.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!
!! OUTPUT
!!  kernel= associated kernel on build (iaction = 1) and get action (iaction = 2).
!!
!! PARENTS
!!      gstate,mklocl_realspace,psolver_hartree,psolver_rhohxc
!!      wvl_wfsinp_reformat
!!
!! CHILDREN
!!      deallocate_coulomb_operator,nullify_coulomb_operator,pkernel_set,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psolver_kernel(hgrid, iaction,  icoulomb, &
& iproc, kernel, mpi_comm, ngfft, nproc, nscforder)
  
 use defs_basis
 use m_profiling_abi
 use m_errors

#if defined HAVE_BIGDFT
 use BigDFT_API,     only  : coulomb_operator,nullify_coulomb_operator, &
&                            deallocate_coulomb_operator,mpi_environment
 use poisson_solver, only : pkernel_init,pkernel_set
#else
 use defs_wvltypes,  only : coulomb_operator
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psolver_kernel'
 use interfaces_14_hidewrite
!End of the abilint section

  implicit none

!Arguments ------------------------------------
  !scalars
  integer,intent(in) :: iaction, icoulomb, mpi_comm, nscforder, iproc, nproc
  !arrays
  integer, intent(in) :: ngfft(3)
  type(coulomb_operator),intent(inout)::kernel
  real(dp),intent(in) :: hgrid(3)

!Local variables-------------------------
#if defined HAVE_BIGDFT
  !scalars
  integer,parameter :: igpu=0 !no GPUs
  !arrays
  integer, save :: kernel_scfOrder 
  integer, save :: kernel_icoulomb
  integer, save :: data_size(3) = (/ -2, -2, -2 /) 
  real(dp), save :: kernel_hgrid(3)   ! Grid step used when creating the kernel.
  character(len = 1) :: geocode
  character(len=500) :: message
  integer :: current_size(3)
  type(coulomb_operator),save :: pkernel, kernelseq 
  type(mpi_environment) :: mpi_env
#endif

! *************************************************************************

#if defined HAVE_BIGDFT

 if (icoulomb == 0) then
!  The kernel is built with 'P'eriodic boundary counditions.
   geocode = 'P'
 else if (icoulomb == 1) then
!  The kernel is built with 'F'ree boundary counditions.
   geocode = 'F'
 else if (icoulomb == 2) then
!  The kernel is built with 'S'urface boundary counditions.
   geocode = 'S'
 end if
 current_size(:) = ngfft(1:3)

!Initialise kernel_array pointer.
 if (maxval(data_size) == -2) then
   call nullify_coulomb_operator(pkernel)
   call nullify_coulomb_operator(kernelseq)
 end if

!If iaction == 0, we free the kernel.
 if (iaction == 0) then
   if (associated(pkernel%kernel)) then
     write(message, "(A)") "Psolver_kernel() : deallocating pkernel..."
     call wrtout(std_out, message,'COLL')
     
     call deallocate_coulomb_operator(pkernel)
   end if
   if (associated(kernelseq%kernel)) then
     write(message, "(A)") "Psolver_kernel() : deallocating kernelseq..."
     call wrtout(std_out, message,'COLL')

     call deallocate_coulomb_operator(kernelseq)
   end if
   data_size = (/ -1, -1, -1 /)
   return
 end if


!Action is build or get. We check the sizes before doing anything else.

!!$!Get the size depending on wavelets calculations or not
!!$ if (dtset%usewvl == 0) then
!!$   hgrid(1) = rprimd(1, 1) / ngfft(1)
!!$   hgrid(2) = rprimd(2, 2) / ngfft(2)
!!$   hgrid(3) = rprimd(3, 3) / ngfft(3)
!!$
!!$ else
!!$   hgrid(:) = 0.5d0 * wvl%h(:)
!!$   current_size(1:3) = (/ wvl%Glr%d%n1i, wvl%Glr%d%n2i, wvl%Glr%d%n3i /)
!!$ end if

!Compute a new kernel if grid size has changed or if the kernel
!has never been computed.
 if ((iaction == 1 .and. .not. associated(pkernel%kernel))    .or. &
& (iaction == 3 .and. .not. associated(kernelseq%kernel)) .or. &
& kernel_icoulomb /= icoulomb  .or. &
& data_size(1)    /= current_size(1) .or. &
& data_size(2)    /= current_size(2) .or. &
& data_size(3)    /= current_size(3) .or. &
& kernel_hgrid(1) /= hgrid(1)        .or. &
& kernel_hgrid(2) /= hgrid(2)        .or. &
& kernel_hgrid(3) /= hgrid(3)        .or. &
& kernel_scfOrder /= nscforder) then
   write(message, "(A,A,A,3I6)") "Psolver_kernel() : building kernel...", ch10, &
&   " | data dimensions:", current_size
   call wrtout(std_out, message, 'COLL')

   if (iaction == 1 .or. iaction == 2) then
     if (associated(pkernel%kernel)) then
       call deallocate_coulomb_operator(pkernel)
     end if
     mpi_env%mpi_comm = mpi_comm
     mpi_env%iproc    = iproc
     mpi_env%nproc    = nproc
     mpi_env%igroup   = 0 ! no task groups
     mpi_env%ngroup   = 1 ! no task groups
     pkernel= pkernel_init(.True.,iproc,nproc,igpu,geocode,&
&     current_size,hgrid,nscforder, mpi_env = mpi_env)
     call pkernel_set(pkernel,.True.)
   end if

   if (iaction == 3 .or. iaction == 4) then
     if (associated(kernelseq%kernel)) then
       call deallocate_coulomb_operator(kernelseq)
     end if
     mpi_env%mpi_comm = mpi_comm
     mpi_env%iproc    = 0
     mpi_env%nproc    = 1
     mpi_env%igroup   = 0 ! no task groups
     mpi_env%ngroup   = 1 ! no task groups
     kernelseq= pkernel_init(.True.,iproc,nproc,igpu,geocode,&
&     current_size,hgrid,nscforder, mpi_env = mpi_env)
     call pkernel_set(kernelseq,.True.)
   end if

!  Storing variables which were used to make the kernel
   kernel_icoulomb = icoulomb
   data_size(:)    = current_size(:)
   kernel_hgrid(:) = hgrid(:)
   kernel_scfOrder = nscforder
 end if
 
 ! Shallow copy if kernel has been associated.
 if (iaction == 1 .or. iaction == 2) then
   kernel = pkernel
 end if
 if (iaction == 3 .or. iaction == 4) then
   kernel = kernelseq
 end if

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*)  iaction,icoulomb,mpi_comm,nscforder,iproc,nproc,ngfft(1),kernel%co,hgrid(1)
#endif

end subroutine psolver_kernel
!!***
