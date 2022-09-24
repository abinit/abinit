!!****m* ABINIT/m_gemm_nonlop
!! NAME
!! m_gemm_nonlop
!!
!! FUNCTION
!!  This module provides functions to compute the nonlocal operator by means of the BLAS GEMM
!!  routine. By treating ndat simultaneous wavefunctions, it is able to exploit BLAS3 routines,
!!  which leads to excellent CPU efficiency and OpenMP scalability.
!!
!! COPYRIGHT
!! Copyright (C) 2014-2022 ABINIT group (AL)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

! TODO list :
! Don't allocate the full nkpt structures, only those that are treated by this proc: use same init as in m_bandfft_kpt
! support more options (forces & stresses mostly)
! Support RF/other computations (only GS right now)
! handle the case where nloalg(2) < 0, ie no precomputation of ph3d
! more systematic checking of the workflow (right now, only works if init/make/gemm/destroy, no multiple makes, etc)
! Avoid allocating the complex matrix when istwfk > 1
! Merge with chebfi's invovl


#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_gemm_nonlop

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_abi_linalg  ! copy_on_gpu, copy_from_gpu, alloc_on_gpu, dealloc_on_gpu, gpu_memset, gpu_allocated

 use defs_abitypes, only : MPI_type
 use m_opernlc_ylm, only : opernlc_ylm
 use m_opernlc_ylm_gpu, only : opernlc_ylm_gpu
 use m_opernlc_ylm_allwf_cpu, only : opernlc_ylm_allwf_cpu
 use m_pawcprj, only : pawcprj_type
 use m_geometry, only : strconv
 use m_kg, only : mkkpg

#if defined(HAVE_GPU_CUDA)
  use m_alloc_hamilt_gpu, only : gemm_nonlop_kokkos
#endif

#ifdef HAVE_FC_ISO_C_BINDING
 use, intrinsic :: iso_c_binding, only : c_ptr, c_int32_t, c_int64_t, c_float, c_double, c_size_t, c_loc
#endif


 implicit none

 private

 ! Use these routines in order: first call init, then call make_gemm_nonlop for each k point,
 ! then call gemm_nonlop to do the actual computation, and call destroy when done. See gstate and vtorho.
 public :: init_gemm_nonlop
 public :: destroy_gemm_nonlop
 public :: make_gemm_nonlop
 public :: gemm_nonlop

#if defined HAVE_GPU_CUDA
 public :: init_gemm_nonlop_gpu
 public :: destroy_gemm_nonlop_gpu
 !public :: make_gemm_nonlop_gpu
 public :: gemm_nonlop_gpu
#endif


!!***

!----------------------------------------------------------------------

!!****t* m_gemm_nonlop/gemm_nonlop_type
!! NAME
!! gemm_nonlop_type
!!
!! FUNCTION
!! Contains information needed to apply the nonlocal operator
!!
!! SOURCE
 type,public :: gemm_nonlop_type

   integer :: nprojs
   integer :: ngrads

   real(c_double), allocatable :: projs(:, :, :)
   ! (2, npw, nprojs)

   real(c_double), allocatable :: projs_r(:, :, :)
   ! (1, npw, nprojs)

   real(c_double), allocatable :: projs_i(:, :, :)
   ! (1, npw, nprojs)

   real(dp), allocatable :: dprojs(:, :, :)
   ! (2, npw, nprojs*ngrads)
   real(dp), allocatable :: dprojs_r(:, :, :)
   ! (1, npw, nprojs*ngrads)
   real(dp), allocatable :: dprojs_i(:, :, :)
   ! (1, npw, nprojs*ngrads)

 end type gemm_nonlop_type
!!***

 type(gemm_nonlop_type), save, public, allocatable, target :: gemm_nonlop_kpt(:)
 !(nkpt)

#if defined(HAVE_FC_ISO_C_BINDING) && defined(HAVE_GPU_CUDA)

 type, bind(c), public :: gemm_nonlop_gpu_type

   integer(kind=c_int32_t) :: npw
   integer(kind=c_int32_t) :: nprojs

   ! array of double on GPU, dimensions are (2, npw, nprojs)
   type(c_ptr) :: projs

   ! array of double on GPU, dimensions are (1, npw, nprojs)
   type(c_ptr) :: projs_r

   ! array of double on GPU, dimensions are (1, npw, nprojs)
   type(c_ptr) :: projs_i

 end type gemm_nonlop_gpu_type

 !! array of size nkpt of sobjects of type gemm_nonlop_gpu_type, array size is nkpt
 type(gemm_nonlop_gpu_type), save, public, allocatable :: gemm_nonlop_kpt_gpu(:)
 !(nkpt)

#endif


 integer, save, public :: gemm_nonlop_ikpt_this_proc_being_treated
 !! This is oh so very crude, but I can't find any other way to do it without passing ikpt deep down to nonlop

 logical, save, public :: gemm_nonlop_use_gemm = .false.
 ! Public variable indicating whether we should call gemm_nonlop or fall back to the usual nonlop. Set to false
 ! in order not to interfere with non-GS calls to nonlop.

 logical, save, public :: gemm_nonlop_use_gemm_gpu = .false.
 ! Public variable controlled by input dataset var named use_gemm_nonlop_cuda (0 or 1).
 ! When 0, nonlop is computed by calling the regular nonlop_gpu
 ! When 1, nonlop is computed by calling gemm_nonlop_gpu

 logical, save, public :: gemm_nonlop_use_kokkos = .false.
 ! public variable controlled by dataset variable use_kokkos - probably removed when kokkos version is debuged

 !!
 !! GPU interface
 !!
#if defined(HAVE_FC_ISO_C_BINDING) && defined(HAVE_GPU_CUDA)

 interface

   !> allocate a object of type gemm_nonlop_gpu_type
   subroutine gemm_nonlop_gpu_allocate(gemm_nonlop_gpu_obj, npw, nprojs, istwf_k) bind(c, name='cuda_gemm_nonlop_gpu_allocate')
     use, intrinsic :: iso_c_binding
     import gemm_nonlop_gpu_type
     implicit none
     type(gemm_nonlop_gpu_type),     intent(inout) :: gemm_nonlop_gpu_obj
     integer(kind=c_int32_t), value, intent(in)    :: npw
     integer(kind=c_int32_t), value, intent(in)    :: nprojs
     integer(kind=c_int32_t), value, intent(in)    :: istwf_k
   end subroutine gemm_nonlop_gpu_allocate

   !> deallocate a object of type gemm_nonlop_gpu_type
   subroutine gemm_nonlop_gpu_deallocate(gemm_nonlop_gpu_obj) bind(c, name='cuda_gemm_nonlop_gpu_deallocate')
     use, intrinsic :: iso_c_binding
     import gemm_nonlop_gpu_type
     implicit none
     type(gemm_nonlop_gpu_type),     intent(inout) :: gemm_nonlop_gpu_obj
   end subroutine gemm_nonlop_gpu_deallocate

   !> extract real part of a complex vector
   !> data_in and data_out must be pointers in device memory
   subroutine extract_real_part(data_out, data_in, size) bind(c, name='cuda_extract_real_part')
     use, intrinsic :: iso_c_binding
     implicit none
     type(c_ptr),             value :: data_out
     type(c_ptr),             value :: data_in
     integer(kind=c_int32_t), value :: size
   end subroutine extract_real_part

   !> extract real part of a complex vector
   !> data_in and data_out must be pointers in device memory
   subroutine extract_imag_part(data_out, data_in, size) bind(c, name='cuda_extract_imag_part')
     use, intrinsic :: iso_c_binding
     implicit none
     type(c_ptr),             value :: data_out
     type(c_ptr),             value :: data_in
     integer(kind=c_int32_t), value :: size
   end subroutine extract_imag_part

   !> insert real part of a complex vector
   !> data_in and data_out must be pointers in device memory
   subroutine insert_real_part(data_out, data_in, size) bind(c, name='cuda_insert_real_part')
     use, intrinsic :: iso_c_binding
     implicit none
     type(c_ptr),             value :: data_out
     type(c_ptr),             value :: data_in
     integer(kind=c_int32_t), value :: size
   end subroutine insert_real_part

   !> insert real part of a complex vector
   !> data_in and data_out must be pointers in device memory
   subroutine insert_imag_part(data_out, data_in, size) bind(c, name='cuda_insert_imag_part')
     use, intrinsic :: iso_c_binding
     implicit none
     type(c_ptr),             value :: data_out
     type(c_ptr),             value :: data_in
     integer(kind=c_int32_t), value :: size
   end subroutine insert_imag_part

   !> data_in and data_out must be pointers in device memory
   subroutine fix_realvec(data, npw_in, ndat_nspinor, option) bind(c, name='cuda_fix_realvec')
     use, intrinsic :: iso_c_binding
     implicit none
     type(c_ptr),                    intent(inout) :: data
     integer(kind=c_int32_t), value, intent(in)    :: npw_in
     integer(kind=c_int32_t), value, intent(in)    :: ndat_nspinor
     integer(kind=c_int32_t), value, intent(in)    :: option
   end subroutine fix_realvec

#if defined(HAVE_KOKKOS)

    subroutine opernlc_ylm_allwf_kokkos(cplex, cplex_enl, cplex_fac, &
      &                                 dimenl1, dimenl2, dimekbq, &
      &                                 iatm, itypat, ntypat, nprojs, &
      &                                 natom, nincat, nspinor, &
      &                                 nspinortot, paw_opt, &
      &                                 nlmn, lmnmax, &
      &                                 enl_gpu, &
      &                                 gx_gpu, &
      &                                 gxfac_gpu, &
      &                                 gxfac2_gpu, &
      &                                 gxfac_sij_gpu, &
      &                                 shift_spinor, ndat, &
      &                                 atindx1_gpu, &
      &                                 indlmn_gpu, &
      &                                 lambda_gpu, &
      &                                 sij_typ_gpu, &
      &                                 shift_proj, &
      &                                 nattyp_max) bind(c, name='opernlc_ylm_allwf_kokkos_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int32_t), value, intent(in)    :: cplex, cplex_enl, cplex_fac
      integer(kind=c_int32_t), value, intent(in)    :: dimenl1, dimenl2, dimekbq
      integer(kind=c_int32_t), value, intent(in)    :: iatm, itypat, ntypat, nprojs
      integer(kind=c_int32_t), value, intent(in)    :: natom, nincat, nspinor
      integer(kind=c_int32_t), value, intent(in)    :: nspinortot, paw_opt
      integer(kind=c_int32_t), value, intent(in)    :: nlmn, lmnmax
      type(c_ptr),             value                :: enl_gpu ! (dimenl1, dimenl2, nspinortot**2, dimekbq)
      type(c_ptr),             value                :: gx_gpu ! (cplex,nlmn,nincat,nspinor*ndat)
      type(c_ptr),             value                :: gxfac_gpu ! (cplex_fac,nlmn,nincat,nspinor*ndat)
      type(c_ptr),             value                :: gxfac2_gpu ! (cplex_fac,nlmn,nincat,nspinor*ndat)
      type(c_ptr),             value                :: gxfac_sij_gpu !(cplex,nlmn,nincat,nspinor*ndat*(paw_opt/3))
      integer(kind=c_int32_t), value, intent(in)    :: shift_spinor, ndat
      type(c_ptr),             value                :: atindx1_gpu ! (natom)
      type(c_ptr),             value                :: indlmn_gpu  ! (6,nlmn)
      type(c_ptr),             value                :: lambda_gpu ! (ndat)
      type(c_ptr),             value                :: sij_typ_gpu ! (((paw_opt+1)/3)*nlmn*(nlmn+1)/2)
      integer(kind=c_int32_t), value, intent(in)    :: shift_proj, nattyp_max
    end subroutine opernlc_ylm_allwf_kokkos

#endif

 end interface

#endif

contains

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop/init_gemm_nonlop
!! NAME
!! init_gemm_nonlop
!!
!! FUNCTION
!! Initalization of the gemm_nonlop_kpt array
!!
!! INPUTS
!! nkpt= number of k-points
!!
!! SOURCE
 subroutine init_gemm_nonlop(nkpt)

  integer,intent(in) :: nkpt
  integer :: ikpt

! *************************************************************************

  ! TODO only allocate the number of kpt treated by this proc
  ABI_MALLOC(gemm_nonlop_kpt, (nkpt))
  do ikpt=1,nkpt
    gemm_nonlop_kpt(ikpt)%nprojs = -1
    gemm_nonlop_kpt(ikpt)%ngrads = -1
  end do

  gemm_nonlop_kokkos % allocated = .false.

 end subroutine init_gemm_nonlop
!!***

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop/destroy_gemm_nonlop
!! NAME
!! destroy_gemm_nonlop
!!
!! FUNCTION
!! Destruction of the gemm_nonlop_kpt array
!!
!! INPUTS
!! nkpt= number of k-points
!!
!! SOURCE
 subroutine destroy_gemm_nonlop(nkpt)

  integer,intent(in) :: nkpt
  integer :: ikpt

! *************************************************************************

! TODO add cycling if kpt parallelism
  do ikpt = 1,nkpt
    call free_gemm_nonlop_ikpt(ikpt)
  end do

 ABI_FREE(gemm_nonlop_kpt)

 end subroutine destroy_gemm_nonlop
!!***

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop/free_gemm_nonlop_ikpt
!! NAME
!! free_destroy_gemm_nonlop_ikpt
!!
!! FUNCTION
!! Release memory for one kpt value of the gemm_nonlop_kpt array
!!
!! INPUTS
!! ikpt= index of gemm_nonlop_kptto be released
!!
!! SOURCE
 subroutine free_gemm_nonlop_ikpt(ikpt)

  integer,intent(in) :: ikpt

! *************************************************************************

 if(gemm_nonlop_kpt(ikpt)%nprojs /= -1) then
   if (allocated(gemm_nonlop_kpt(ikpt)%projs)) then
     ABI_FREE(gemm_nonlop_kpt(ikpt)%projs)
   end if
   if (allocated(gemm_nonlop_kpt(ikpt)%projs_r)) then
     ABI_FREE(gemm_nonlop_kpt(ikpt)%projs_r)
   end if
   if (allocated(gemm_nonlop_kpt(ikpt)%projs_i)) then
   ABI_FREE(gemm_nonlop_kpt(ikpt)%projs_i)
   end if
   gemm_nonlop_kpt(ikpt)%nprojs = -1
   if(gemm_nonlop_kpt(ikpt)%ngrads /= -1) then
     if (allocated(gemm_nonlop_kpt(ikpt)%dprojs)) then
       ABI_FREE(gemm_nonlop_kpt(ikpt)%dprojs)
     end if
     if (allocated(gemm_nonlop_kpt(ikpt)%dprojs_r)) then
       ABI_FREE(gemm_nonlop_kpt(ikpt)%dprojs_r)
     end if
     if (allocated(gemm_nonlop_kpt(ikpt)%dprojs_i)) then
       ABI_FREE(gemm_nonlop_kpt(ikpt)%dprojs_i)
     end if
     gemm_nonlop_kpt(ikpt)%ngrads = -1
   end if
 end if

 end subroutine free_gemm_nonlop_ikpt
!!***

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop/make_gemm_nonlop
!! NAME
!! make_gemm_nonlop
!!
!! FUNCTION
!! Build the gemm_nonlop array
!!
!! INPUTS
!!
!! SOURCE
 subroutine make_gemm_nonlop(ikpt,npw,lmnmax,ntypat,indlmn,nattyp,istwf_k,ucvol,ffnl_k, &
&                            ph3d_k,kpt_k,kg_k,kpg_k,use_gemm_nonlop_gpu, &
&                            compute_grad_strain,compute_grad_atom) ! Optional parameters

  integer, intent(in) :: ikpt
  integer, intent(in) :: npw, lmnmax, ntypat
  integer, intent(in) :: indlmn(:,:,:), kg_k(:,:)
  integer, intent(in) :: nattyp(ntypat)
  integer, intent(in) :: istwf_k
  logical, intent(in), optional :: compute_grad_strain,compute_grad_atom
  real(dp), intent(in) :: ucvol
  real(dp), intent(in) :: ffnl_k(:,:,:,:)
  real(dp), intent(in) :: ph3d_k(:,:,:)
  real(dp), intent(in) :: kpt_k(:)
  real(dp), intent(in), target :: kpg_k(:,:)
  integer, intent(in) :: use_gemm_nonlop_gpu

  integer :: nprojs,ndprojs,ngrads

  integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
  integer :: itypat, ilmn, nlmn, ia, iaph3d, igrad, shift, shift_grad
  integer :: il, ipw, idir, idir1, idir2, nkpg_local
  logical :: parity,my_compute_grad_strain,my_compute_grad_atom
  real(dp),allocatable :: atom_projs(:,:,:), atom_dprojs(:,:,:,:), temp(:)
  real(dp),pointer :: kpg(:,:)

! *************************************************************************

  my_compute_grad_strain=.false. ; if (present(compute_grad_strain)) my_compute_grad_strain=compute_grad_strain
  my_compute_grad_atom=.false. ; if (present(compute_grad_atom)) my_compute_grad_atom=compute_grad_atom
  ABI_CHECK(size(ph3d_k)>0,'nloalg(2)<0 not compatible with use_gemm_nonlop=1!')
!  ABI_CHECK((.not.my_compute_grad_strain).or.size(kpg_k)>0,"kpg_k should be allocated to compute gradients!")
!  ABI_CHECK((.not.my_compute_grad_atom).or.size(kpg_k)>0,"kpg_k should be allocated to compute gradients!")

  iaph3d = 1

  ABI_MALLOC(atom_projs, (2, npw, lmnmax))
  if (my_compute_grad_strain) then
    ndprojs = 3
    ABI_MALLOC(atom_dprojs, (2, npw, ndprojs, lmnmax))
  else
    ndprojs = 0
  end if

  ABI_MALLOC(temp, (npw))

  call free_gemm_nonlop_ikpt(ikpt)

  ! build nprojs, ngrads
  nprojs = 0 ; ngrads = 0
  do itypat=1,ntypat
    nprojs = nprojs + count(indlmn(3,:,itypat)>0)*nattyp(itypat)
  end do
  if (my_compute_grad_strain) ngrads=ngrads+6
  if (my_compute_grad_atom) ngrads=ngrads+3
  if (nprojs>0) gemm_nonlop_kpt(ikpt)%nprojs = nprojs
  if (ngrads>0) gemm_nonlop_kpt(ikpt)%ngrads = ngrads

  if(istwf_k <= 1) then
    ABI_MALLOC(gemm_nonlop_kpt(ikpt)%projs, (2, npw, nprojs))
    gemm_nonlop_kpt(ikpt)%projs = zero
    if(ngrads>0) then
      ABI_MALLOC(gemm_nonlop_kpt(ikpt)%dprojs, (2, npw, nprojs*ngrads))
      gemm_nonlop_kpt(ikpt)%dprojs = zero
    end if
  else
    ABI_MALLOC(gemm_nonlop_kpt(ikpt)%projs_r, (1, npw, nprojs))
    ABI_MALLOC(gemm_nonlop_kpt(ikpt)%projs_i, (1, npw, nprojs))
    gemm_nonlop_kpt(ikpt)%projs_r = zero
    gemm_nonlop_kpt(ikpt)%projs_i = zero
    if(ngrads>0) then
      ABI_MALLOC(gemm_nonlop_kpt(ikpt)%dprojs_r, (1, npw, nprojs*ngrads))
      ABI_MALLOC(gemm_nonlop_kpt(ikpt)%dprojs_i, (1, npw, nprojs*ngrads))
      gemm_nonlop_kpt(ikpt)%dprojs_r = zero
      gemm_nonlop_kpt(ikpt)%dprojs_i = zero
    end if
  end if

  ! allocate gemm_nonlop_kpt_gpu
#if defined HAVE_GPU_CUDA
  if (use_gemm_nonlop_gpu == 1) then
    gemm_nonlop_kpt_gpu(ikpt)%npw    = npw
    gemm_nonlop_kpt_gpu(ikpt)%nprojs = nprojs
    !call gemm_nonlop_gpu_allocate(gemm_nonlop_kpt_gpu(ikpt), npw, nprojs, istwf_k)
    call alloc_on_gpu(gemm_nonlop_kpt_gpu(ikpt)%projs, 2*npw*nprojs*dp)
    if (istwf_k > 1) then
      call alloc_on_gpu(gemm_nonlop_kpt_gpu(ikpt)%projs_r, 1*npw*nprojs*dp)
      call alloc_on_gpu(gemm_nonlop_kpt_gpu(ikpt)%projs_i, 1*npw*nprojs*dp)
    end if
  end if
#endif


  ! Compute (k+G) vectors if needed
  nkpg_local=0
  if ((my_compute_grad_strain.or.my_compute_grad_atom).and.size(kpg_k)==0) then
    nkpg_local=3
    ABI_MALLOC(kpg,(npw,nkpg_local))
    call mkkpg(kg_k,kpg,kpt_k,nkpg_local,npw)
  else
    kpg => kpg_k
  end if

  shift = 0 ; shift_grad = 0
  do itypat = 1, ntypat
    nlmn = count(indlmn(3,:,itypat)>0)

    do ia = 1, nattyp(itypat)

      !! build atom_projs, from opernlb
      !! P = 4pi/sqrt(ucvol)* conj(diag(ph3d)) * ffnl * diag(parity), with parity = (-i)^l

      ! start from 4pi/sqrt(ucvol)*ffnl
      ! atom_projs(1, :, 1:nlmn) = four_pi/sqrt(ham%ucvol) * ham%ffnl_k(:, 1, 1:nlmn)
      ! TODO vectorize (DCOPY with stride)
      atom_projs(:,:,:) = zero
      do ipw=1, npw
        atom_projs(1,ipw, 1:nlmn) = four_pi/sqrt(ucvol) * ffnl_k(ipw, 1, 1:nlmn, itypat)
      end do
      if (ndprojs>0) atom_dprojs(:,:,:,:) = zero
      if (my_compute_grad_strain) then
        do ipw=1, npw
          atom_dprojs(1,ipw, 1:3, 1:nlmn) = four_pi/sqrt(ucvol) * ffnl_k(ipw, 2:4, 1:nlmn, itypat)
        end do
      end if

      ! multiply by (-i)^l
      do ilmn=1,nlmn
        il=mod(indlmn(1,ilmn, itypat),4);
        parity=(mod(il,2)==0)
        if (il>1) then
          ! multiply by -1
          atom_projs(:,:,ilmn) = -atom_projs(:,:,ilmn)
          if (ndprojs>0) atom_dprojs(:,:,:,ilmn) = -atom_dprojs(:,:,:,ilmn)
        end if
        if(.not. parity) then
          ! multiply by -i
          temp = atom_projs(2,:,ilmn)
          atom_projs(2,:,ilmn) = -atom_projs(1,:,ilmn)
          atom_projs(1,:,ilmn) =  temp
          if (ndprojs>0) then
            do idir=1,ndprojs
              temp = atom_dprojs(2,:,idir,ilmn)
              atom_dprojs(2,:,idir,ilmn) = -atom_dprojs(1,:,idir,ilmn)
              atom_dprojs(1,:,idir,ilmn) =  temp
            end do
          end if
        end if
      end do

      ! multiply by conj(ph3d)
      do ilmn=1,nlmn
        temp = atom_projs(1, :, ilmn)
        atom_projs(1, :, ilmn) = atom_projs(1, :, ilmn) * ph3d_k(1, :, iaph3d) &
&                              + atom_projs(2, :, ilmn) * ph3d_k(2, :, iaph3d)
        atom_projs(2, :, ilmn) = atom_projs(2, :, ilmn) * ph3d_k(1, :, iaph3d) &
&                              - temp                   * ph3d_k(2, :, iaph3d)
      end do
      if (ndprojs>0) then
        do ilmn=1,nlmn
          do idir=1,ndprojs
            temp = atom_dprojs(1, :, idir,ilmn)
            atom_dprojs(1, :, idir,ilmn) = atom_dprojs(1, :, idir,ilmn) * ph3d_k(1, :, iaph3d) &
&                                        + atom_dprojs(2, :, idir,ilmn) * ph3d_k(2, :, iaph3d)
            atom_dprojs(2, :, idir,ilmn) = atom_dprojs(2, :, idir,ilmn) * ph3d_k(1, :, iaph3d) &
&                                        - temp                         * ph3d_k(2, :, iaph3d)
          end do
        end do
      end if

      !! atom_projs is built, copy to projs / dprojs

      if(istwf_k <= 1) then
        gemm_nonlop_kpt(ikpt)%projs(1:2, :, shift+1:shift+nlmn) = atom_projs(:, :, 1:nlmn)
        if(ngrads>0) then
          igrad=0
          if(my_compute_grad_strain) then
            do idir=1,6
              idir1=alpha(idir);idir2=beta(idir)
              do ilmn=1,nlmn
                do ipw=1,npw
                  gemm_nonlop_kpt(ikpt)%dprojs(1:2, ipw, shift_grad+nlmn*igrad+ilmn) = &
&                  -half*(atom_dprojs(1:2, ipw, idir1, ilmn)*kpg(ipw,idir2) &
&                        +atom_dprojs(1:2, ipw, idir2, ilmn)*kpg(ipw,idir1))
                end do
              end do
              igrad=igrad+1
            end do
          end if
          if(my_compute_grad_atom) then
            do idir=1,3
              do ilmn=1,nlmn
                do ipw=1,npw
                  gemm_nonlop_kpt(ikpt)%dprojs(1, ipw, shift_grad+nlmn*igrad+ilmn) = &
&                  +atom_projs(2, ipw, ilmn)*kpg(ipw,idir)*two_pi
                  gemm_nonlop_kpt(ikpt)%dprojs(2, ipw, shift_grad+nlmn*igrad+ilmn) = &
&                  -atom_projs(1, ipw, ilmn)*kpg(ipw,idir)*two_pi
                end do
              end do
              igrad=igrad+1
            end do
          end if
        end if

      else ! istwf_k>1
        gemm_nonlop_kpt(ikpt)%projs_r(1, :, shift+1:shift+nlmn) = atom_projs(1, :, 1:nlmn)
        gemm_nonlop_kpt(ikpt)%projs_i(1, :, shift+1:shift+nlmn) = atom_projs(2, :, 1:nlmn)
        if(ngrads>0) then
          igrad=0
          if(my_compute_grad_strain) then
            do idir=1,6
              idir1=alpha(idir);idir2=beta(idir)
              do ilmn=1,nlmn
                do ipw=1,npw
                  gemm_nonlop_kpt(ikpt)%dprojs_r(1, ipw, shift_grad+nlmn*igrad+ilmn) = &
&                  -half*(atom_dprojs(1, ipw, idir1, ilmn)*kpg(ipw,idir2) &
&                        +atom_dprojs(1, ipw, idir2, ilmn)*kpg(ipw,idir1))

                  gemm_nonlop_kpt(ikpt)%dprojs_i(1, ipw, shift_grad+nlmn*igrad+ilmn) = &
&                  -half*(atom_dprojs(2, ipw, idir1, ilmn)*kpg(ipw,idir2) &
&                        +atom_dprojs(2, ipw, idir2, ilmn)*kpg(ipw,idir1))
                end do
              end do
              igrad=igrad+1
            end do
          end if
          if(my_compute_grad_atom) then
            do idir=1,3
              do ilmn=1,nlmn
                do ipw=1,npw
                  gemm_nonlop_kpt(ikpt)%dprojs_r(1, ipw, shift_grad+nlmn*igrad+ilmn) = &
&                  +atom_projs(2, ipw, ilmn)*kpg(ipw,idir)*two_pi
                  gemm_nonlop_kpt(ikpt)%dprojs_i(1, ipw, shift_grad+nlmn*igrad+ilmn) = &
&                  -atom_projs(1, ipw, ilmn)*kpg(ipw,idir)*two_pi
                end do
              end do
              igrad=igrad+1
            end do
          end if
        end if
      end if ! istwf_k

      shift = shift + nlmn
      shift_grad = shift_grad + ngrads*nlmn
      iaph3d = iaph3d + 1
    end do
  end do

#if defined HAVE_GPU_CUDA
  if (use_gemm_nonlop_gpu == 1) then
    ! upload data to gpu memory
    call copy_on_gpu(c_loc(gemm_nonlop_kpt(ikpt)%projs(1,1,1)), gemm_nonlop_kpt_gpu(ikpt)%projs, 2*npw*nprojs*dp)
    if(istwf_k > 1) then
      call copy_on_gpu(c_loc(gemm_nonlop_kpt(ikpt)%projs_r(1,1,1)), gemm_nonlop_kpt_gpu(ikpt)%projs_r, 1*npw*nprojs*dp)
      call copy_on_gpu(c_loc(gemm_nonlop_kpt(ikpt)%projs_i(1,1,1)), gemm_nonlop_kpt_gpu(ikpt)%projs_i, 1*npw*nprojs*dp)
    end if
  end if ! use_gemm_nonlop_gpu == 1
#endif

  ABI_FREE(atom_projs)
  ABI_FREE(temp)
  if (allocated(atom_dprojs)) then
    ABI_FREE(atom_dprojs)
  end if
  if (nkpg_local>0) then
    ABI_FREE(kpg)
  end if

 end subroutine make_gemm_nonlop
!!***

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop/gemm_nonlop
!! NAME
!! gemm_nonlop
!!
!! FUNCTION
!! Replacement of nonlop. same prototype as nonlop although not all options are implemented.
!!
!! INPUTS
!!
!! SOURCE

 subroutine gemm_nonlop(atindx1,choice,cpopt,cprjin,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,&
&                 enl,enlout,ffnlin,ffnlout,gmet,gprimd,idir,indlmn,istwf_k,&
&                 kgin,kgout,kpgin,kpgout,kptin,kptout,lambda,lmnmax,matblk,mgfft,&
&                 mpi_enreg,mpsang,mpssoang,natom,nattyp,ndat,ngfft,nkpgin,nkpgout,nloalg,&
&                 nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO,paw_opt,phkxredin,&
&                 phkxredout,ph1d,ph3din,ph3dout,signs,sij,svectout,&
&                 tim_nonlop,ucvol,useylm,vectin,vectout,&
&                 use_gpu_cuda)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: choice,cpopt,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,idir
  integer,intent(in) :: istwf_k,lmnmax,matblk,mgfft,mpsang,mpssoang,natom,ndat,nkpgin
  integer,intent(in) :: nkpgout,nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO
  integer,intent(in) :: paw_opt,signs,tim_nonlop,useylm
  integer,optional,intent(in) :: use_gpu_cuda
  real(dp),intent(in) :: lambda(ndat),ucvol
  type(MPI_type),intent(in) :: mpi_enreg
  !arrays
  integer,intent(in) :: atindx1(natom),indlmn(6,lmnmax,ntypat),kgin(3,npwin)
  integer,intent(in) :: kgout(3,npwout),nattyp(ntypat),ngfft(18),nloalg(3)
  real(dp),intent(in) :: enl(dimenl1,dimenl2,nspinortot**2,dimekbq)
  real(dp),intent(in) :: ffnlin(npwin,dimffnlin,lmnmax,ntypat)
  real(dp),intent(in) :: ffnlout(npwout,dimffnlout,lmnmax,ntypat),gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3),kpgin(npwin,nkpgin*useylm)
  real(dp),intent(in) :: kpgout(npwout,nkpgout*useylm),kptin(3),kptout(3)
  real(dp),intent(in) :: phkxredin(2,natom),phkxredout(2,natom),ph1d(2,*)
  real(dp),intent(in) :: sij(dimenl1,ntypat*((paw_opt+1)/3))
  real(dp),intent(inout) :: ph3din(2,npwin,matblk),ph3dout(2,npwout,matblk)
  real(dp),intent(inout) :: vectin(2,npwin*nspinor*ndat)
  real(dp),intent(inout) :: enlout(nnlout*ndat)
  real(dp),intent(out) :: svectout(2,npwout*nspinor*(paw_opt/3)*ndat)
  real(dp),intent(inout) :: vectout(2,npwout*nspinor*ndat) !vz_i
  type(pawcprj_type),intent(inout) :: cprjin(natom,nspinor*((cpopt+5)/5)*ndat)

  ! locals
  integer :: ii, ia, idat, igrad, nprojs, ngrads, shift, iatom, nlmn, ierr, ibeg, iend
  integer :: cplex, cplex_enl, cplex_fac, proj_shift, grad_shift
  integer :: enlout_shift, force_shift, nnlout_test
  integer :: iatm, ndgxdt, ndgxdtfac, nd2gxdt, nd2gxdtfac, optder, itypat, ilmn
  integer :: cplex_dgxdt(1), cplex_d2gxdt(1)
  real(dp) :: esum
  real(dp) :: work(6)
  real(dp) :: dgxdt_dum_in(1,1,1,1,1), dgxdt_dum_out(1,1,1,1,1),dgxdt_dum_out2(1,1,1,1,1)
  real(dp) :: d2gxdt_dum_in(1,1,1,1,1), d2gxdt_dum_out(1,1,1,1,1),d2gxdt_dum_out2(1,1,1,1,1)
  real(dp), allocatable :: enlk(:),sij_typ(:)
  real(dp), allocatable :: projections(:,:,:), s_projections(:,:,:), vnl_projections(:,:,:)
  real(dp), allocatable :: dprojections(:,:,:), temp_realvec(:)

! *************************************************************************

  ! We keep the same interface as nonlop, but we don't use many of those
  ABI_UNUSED((/ffnlin,ffnlout,gmet,kpgin,kpgout/))
  ABI_UNUSED((/ph1d(1,1),ph3din,ph3dout/))
  ABI_UNUSED((/phkxredin,phkxredout,ucvol/))
  ABI_UNUSED((/mgfft,mpsang,mpssoang/))
  ABI_UNUSED((/kptin,kptout/))
  ABI_UNUSED((/idir,nloalg,ngfft,kgin,kgout,ngfft,only_SO,tim_nonlop,use_gpu_cuda/))

  ! Check supported options
  if (.not.gemm_nonlop_use_gemm) then
    ABI_BUG('computation not prepared for gemm_nonlop use!')
  end if
  if ( (choice>1.and.choice/=7.and.signs==2) .or. &
&      (choice>3.and.choice/=7.and.choice/=23.and.signs==1) .or. &
&      (useylm/=1) ) then
    ABI_BUG('gemm_nonlop option not supported!')
  end if
  if (signs==1) then
    nnlout_test=0
    if (choice==1) nnlout_test=1
    if (choice==2) nnlout_test=3*natom
    if (choice==3) nnlout_test=6
    if (choice==23) nnlout_test=6+3*natom
    if (nnlout<nnlout_test) then
      ABI_BUG('wrong nnlout size!')
    end if
  end if

  cplex=2;if (istwf_k>1) cplex=1
  cplex_enl=1;if (paw_opt>0) cplex_enl=2*dimenl1/(lmnmax*(lmnmax+1)) ! is enl complex?
  cplex_fac=max(cplex,dimekbq)
  if ((nspinortot==2.or.cplex_enl==2).and.paw_opt>0.and.choice/=7) cplex_fac=2 ! is vnl_projections complex?

  nprojs = gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%nprojs
  ngrads = gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%ngrads
  if(choice==1) ngrads=0
  ABI_CHECK(ngrads>=3.or.choice/=2 ,"Bad allocation in gemm_nonlop (2)!")
  ABI_CHECK(ngrads>=6.or.choice/=3 ,"Bad allocation in gemm_nonlop (3)!")
  ABI_CHECK(ngrads>=9.or.choice/=23,"Bad allocation in gemm_nonlop (23)!")

  ! These will store the non-local factors for vectin, svectout and vectout respectively
  ABI_MALLOC(projections,(cplex, nprojs,nspinor*ndat))
  ABI_MALLOC(s_projections,(cplex, nprojs,nspinor*ndat))
  ABI_MALLOC(vnl_projections,(cplex_fac, nprojs,nspinor*ndat))
  projections = zero
  s_projections = zero
  vnl_projections = zero
  if (signs==1.and.ngrads>0) then
    ABI_MALLOC(dprojections,(cplex, ngrads*nprojs,nspinor*ndat))
    dprojections = zero
  end if

  if(nprojs == 0) then
    ! TODO check if this is correct
    if(signs == 1) then
      enlout=zero
      return
    end if
    if(signs == 2) then
      vectout = zero
      if(paw_opt>0) svectout = vectin
      return
    end if
  end if

  ! determine precisely when temp_realvec needs to be allocated
  ! to factorize allocate (resp. deallocate) at the begining (resp. at the end) of subroutine
  ! to avoid multiple allocate/deallocate that can be costly
  if (cplex /= 2) then
    if ( (cpopt < 2) .or. &
      &  (paw_opt == 3 .or. paw_opt == 4) .or. &
      &  (paw_opt == 0 .or. paw_opt == 1 .or. paw_opt == 4)) then
       ABI_MALLOC(temp_realvec,(MAX(npwout,npwin)*nspinor*ndat))
    end if
  end if


  if(cpopt >= 2) then
    ! retrieve from cprjin
    do idat=1, ndat*nspinor
      shift = 0
      do iatom = 1, natom
        nlmn = cprjin(iatom, idat)%nlmn
        projections(1:cplex, shift+1:shift+nlmn, idat) = cprjin(iatom, idat)%cp(1:cplex, 1:nlmn)
        shift = shift + nlmn
      end do
    end do
    if(cpopt==4.and.allocated(dprojections)) then
      ABI_CHECK(cprjin(1,1)%ncpgr>=ngrads,"cprjin%ncpgr not correct! (1)")
      do idat=1, ndat*nspinor
        shift = 0
        do iatom = 1, natom
          nlmn  = cprjin(iatom, idat)%nlmn
          do igrad=1,ngrads
            dprojections(1:cplex, shift+1:shift+nlmn, idat) = &
&                   cprjin(iatom, idat)%dcp(1:cplex,igrad,1:ilmn)
            shift = shift + nlmn
          end do
        end do
      end do
    end if
  else ! cpopt < 2
    ! opernla
    if(cplex == 2) then
      call abi_zgemm_2r('C', 'N', nprojs, ndat*nspinor, npwin, cone, &
&                gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%projs, npwin,&
&                vectin, npwin, czero, projections, nprojs)
      if(signs==1.and.ngrads>0) then
        call abi_zgemm_2r('C', 'N', ngrads*nprojs, ndat*nspinor, npwin, cone, &
                 gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%dprojs, npwin,&
                 vectin, npwin, czero, dprojections, ngrads*nprojs)
      end if
    else

      if (.not. allocated(temp_realvec)) then
        ABI_ERROR("Please provide memory allocation for temp_realvec array")
      end if

      ! only compute real part of projections = P^* psi => projections_r = P_r^T psi_r + P_i^T psi_i
      temp_realvec(1:npwin*nspinor*ndat) = vectin(1,1:npwin*nspinor*ndat)
      if(istwf_k == 2 .and. mpi_enreg%me_g0 == 1) then
        do idat=1, ndat*nspinor
          temp_realvec(1+(idat-1)*npwin) = temp_realvec(1+(idat-1)*npwin)/2
        end do
      end if
      call DGEMM('T', 'N', nprojs, ndat*nspinor, npwin, one, &
&                gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%projs_r, npwin, &
&                temp_realvec, npwin, zero, projections, nprojs)
      if(signs==1.and.ngrads>0) then
        call DGEMM('T', 'N', ngrads*nprojs, ndat*nspinor, npwin, one, &
&                  gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%dprojs_r, npwin, &
&                  temp_realvec, npwin, zero, dprojections, ngrads*nprojs)
      end if
      temp_realvec(1:npwin*nspinor*ndat) = vectin(2,1:npwin*nspinor*ndat)
      if(istwf_k == 2 .and. mpi_enreg%me_g0 == 1) then
        do idat=1, ndat*nspinor
          temp_realvec(1+(idat-1)*npwin) = zero
        end do
      end if
      call DGEMM('T', 'N', nprojs, ndat*nspinor, npwin, one, &
&                gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%projs_i, npwin, &
&                temp_realvec, npwin, one , projections, nprojs)
      projections = projections * 2
      if(signs==1.and.ngrads>0) then
        call DGEMM('T', 'N', ngrads*nprojs, ndat*nspinor, npwin, one, &
&                  gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%dprojs_i, npwin, &
&                  temp_realvec, npwin, one , dprojections, ngrads*nprojs)
        dprojections = dprojections * 2
      end if
    end if ! cplex == 2
    call xmpi_sum(projections,mpi_enreg%comm_fft,ierr)
    if (choice>1) then
      call xmpi_sum(dprojections,mpi_enreg%comm_fft,ierr)
    end if

    if(cpopt >= 0) then
      ! store in cprjin
      do idat=1, ndat*nspinor
        shift = 0
        do iatom = 1, natom
          nlmn = cprjin(iatom, idat)%nlmn
          cprjin(iatom, idat)%cp(1:cplex, 1:nlmn) = projections(1:cplex, shift+1:shift+nlmn, idat)
          shift = shift + nlmn
        end do
      end do
      if(cpopt==3) then
        ABI_CHECK(cprjin(1,1)%ncpgr>=ngrads,"cprjin%ncpgr not correct! (2)")
        shift = 0
        do iatom = 1, natom
          nlmn = cprjin(iatom, idat)%nlmn
          do igrad=1,ngrads
            cprjin(iatom, idat)%dcp(1:cplex,igrad,1:nlmn) = &
&              dprojections(1:cplex, shift+1:shift+nlmn, idat)
            shift = shift + nlmn
          end do
        end do
      end if
    end if ! cpopt >= 0
  end if ! cpopt >= 2

  if(choice > 0) then

    if(choice /= 7) then
      ! opernlc
      iatm = 0
      ndgxdt = 0
      ndgxdtfac = 0
      nd2gxdt = 0
      nd2gxdtfac = 0
      optder = 0

      shift = 0
      ABI_MALLOC(sij_typ,(((paw_opt+1)/3)*lmnmax*(lmnmax+1)/2))
      do itypat=1, ntypat
        nlmn=count(indlmn(3,:,itypat)>0)
        if (paw_opt>=2) then
          if (cplex_enl==1) then
            do ilmn=1,nlmn*(nlmn+1)/2
              sij_typ(ilmn)=sij(ilmn,itypat)
            end do
          else
            do ilmn=1,nlmn*(nlmn+1)/2
              sij_typ(ilmn)=sij(2*ilmn-1,itypat)
            end do
          end if
        end if

        ibeg = shift+1
        iend = shift+nattyp(itypat)*nlmn

        do idat = 1,ndat
          call opernlc_ylm(atindx1,cplex,cplex_dgxdt,cplex_d2gxdt,cplex_enl,cplex_fac,dgxdt_dum_in,dgxdt_dum_out,dgxdt_dum_out2,&
&         d2gxdt_dum_in,d2gxdt_dum_out,d2gxdt_dum_out2,dimenl1,dimenl2,dimekbq,enl,&
&         projections(:, ibeg:iend, 1+nspinor*(idat-1):nspinor*idat),&
&         vnl_projections(:, ibeg:iend,1+nspinor*(idat-1):nspinor*idat),&
&         s_projections(:, ibeg:iend,1+nspinor*(idat-1):nspinor*idat),&
&         iatm,indlmn(:,:,itypat),itypat,lambda(idat),mpi_enreg,natom,ndgxdt,ndgxdtfac,nd2gxdt,nd2gxdtfac,&
&         nattyp(itypat),nlmn,nspinor,nspinortot,optder,paw_opt,sij_typ)
        end do

        shift = shift + nattyp(itypat)*nlmn
        iatm = iatm+nattyp(itypat)
      end do
      ABI_FREE(sij_typ)
    else
      s_projections = projections
    end if ! choice /= 7

    ! opernlb (only choice=1)
    if(signs==2) then
      if(paw_opt == 3 .or. paw_opt == 4) then

        ! Get svectout from s_projections
        if(cplex == 2) then

          call abi_zgemm_2r('N', 'N', npwout, ndat*nspinor, nprojs, cone, &
&                        gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%projs, npwout, &
&                        s_projections, nprojs, czero, svectout, npwout)
        else

          if (.not. allocated(temp_realvec)) then
            ABI_ERROR("Please provide memory allocation for temp_realvec array")
          end if

          call DGEMM('N', 'N', npwout, ndat*nspinor, nprojs, one, &
&                    gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%projs_r, npwout, &
&                    s_projections, nprojs, zero, temp_realvec, npwout)
          svectout(1,1:npwout*nspinor*ndat) = temp_realvec(1:npwout*nspinor*ndat)
          call DGEMM('N', 'N', npwout, ndat*nspinor, nprojs, one, &
&                    gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%projs_i, npwout,&
&                    s_projections, nprojs, zero, temp_realvec, npwout)
          svectout(2,1:npwout*nspinor*ndat) = temp_realvec(1:npwout*nspinor*ndat)

        end if ! cplex = 2
        if(choice /= 7) svectout = svectout + vectin ! TODO understand this

      end if  ! (paw_opt == 3 .or. paw_opt == 4)

      if(paw_opt == 0 .or. paw_opt == 1 .or. paw_opt == 4) then
        ! Get vectout from vnl_projections
        if(cplex_fac == 2) then

          call abi_zgemm_2r('N', 'N', npwout, ndat*nspinor, nprojs, cone, &
&                        gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%projs, npwout, &
&                        vnl_projections, nprojs, czero, vectout, npwout)
        else

          if (.not. allocated(temp_realvec)) then
            ABI_ERROR("Please provide memory allocation for temp_realvec array")
          end if

          call DGEMM('N', 'N', npwout, ndat*nspinor, nprojs, one, &
&                    gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%projs_r, npwout, &
&                    vnl_projections, nprojs, zero, temp_realvec, npwout)
          vectout(1,1:npwout*nspinor*ndat) = temp_realvec(1:npwout*nspinor*ndat)
          call DGEMM('N', 'N', npwout, ndat*nspinor, nprojs, one, &
&                    gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%projs_i, npwout, &
&                    vnl_projections, nprojs, zero, temp_realvec, npwout)
          vectout(2,1:npwout*nspinor*ndat) = temp_realvec(1:npwout*nspinor*ndat)

        end if ! cplex_fac == 2

      end if  ! (paw_opt == 0 .or. paw_opt == 1 .or. paw_opt == 4)
    end if ! opernlb

    ! opernld
    if(signs==1) then
      enlout=zero
      if(choice==1.or.choice==3.or.choice==23) then
        ABI_MALLOC(enlk,(ndat))
        enlk=zero
        do idat=1,ndat*nspinor
          proj_shift=0
          do itypat=1, ntypat
            nlmn=count(indlmn(3,:,itypat)>0)
            do ia=1,nattyp(itypat)
              !Following loops are a [D][Z]DOT
              esum=zero
              do ilmn=1,nlmn
                do ii=1,cplex
                  esum=esum +vnl_projections(ii,proj_shift+ilmn,idat) &
&                           *projections    (ii,proj_shift+ilmn,idat)
                end do
              end do
              proj_shift=proj_shift+nlmn
              enlk(idat) = enlk(idat) + esum
            end do
          end do
        end do
        if (choice==1) enlout(1:ndat)=enlk(1:ndat)
      end if ! choice=1/3/23
      if(choice==2.or.choice==3.or.choice==23) then
        do idat=1,ndat*nspinor
          proj_shift=0 ; grad_shift=0
          enlout_shift=(idat-1)*nnlout
          force_shift=merge(6,0,choice==23)
          do itypat=1, ntypat
            nlmn=count(indlmn(3,:,itypat)>0)
            do ia=1,nattyp(itypat)
              if (choice==3.or.choice==23) then
                do igrad=1,6
                  !Following loops are a [D][Z]DOT
                  esum=zero
                  do ilmn=1,nlmn
                    do ii=1,cplex
                      esum=esum +vnl_projections(ii,proj_shift+ilmn,idat) &
&                               *dprojections   (ii,grad_shift+ilmn,idat)
                    end do
                  end do
                  grad_shift=grad_shift+nlmn
                  enlout(enlout_shift+igrad)=enlout(enlout_shift+igrad) + two*esum
                end do
              end if
              if (choice==2.or.choice==23) then
                do igrad=1,3
                  !Following loops are a [D][Z]DOT
                  esum=zero
                  do ilmn=1,nlmn
                    do ii=1,cplex
                      esum=esum +vnl_projections(ii,proj_shift+ilmn,idat) &
&                               *dprojections   (ii,grad_shift+ilmn,idat)
                    end do
                  end do
                  grad_shift=grad_shift+nlmn
                  enlout(enlout_shift+force_shift+igrad)= &
&                               enlout(enlout_shift+force_shift+igrad) + two*esum
                end do
                force_shift=force_shift+3
              end if
              proj_shift=proj_shift+nlmn
            end do
          end do
        end do
      end if ! choice=2, 3 or 23

    end if !opernld

  end if ! choice>0

! Reduction in case of parallelism
  if (signs==1.and.mpi_enreg%paral_spinor==1) then
    if (size(enlout)>0) then
      call xmpi_sum(enlout,mpi_enreg%comm_spinor,ierr)
    end if
    if (choice==3.or.choice==23) then
      call xmpi_sum(enlk,mpi_enreg%comm_spinor,ierr)
    end if
  end if

! Derivatives wrt strain
!  - Convert from reduced to cartesian coordinates
!  - Substract volume contribution
 if ((choice==3.or.choice==23).and.signs==1.and.paw_opt<=3) then
   do idat=1,ndat
     enlout_shift=(idat-1)*nnlout
     call strconv(enlout(enlout_shift+1:enlout_shift+6),gprimd,work)
     enlout(enlout_shift+1:enlout_shift+3)=(work(1:3)-enlk(idat))
     enlout(enlout_shift+4:enlout_shift+6)= work(4:6)
   end do
 end if

! Release memory
  ABI_FREE(projections)
  ABI_FREE(s_projections)
  ABI_FREE(vnl_projections)
  if (allocated(dprojections)) then
    ABI_FREE(dprojections)
  end if
  if (allocated(enlk)) then
    ABI_FREE(enlk)
  end if
  if (allocated(temp_realvec)) then
    ABI_FREE(temp_realvec)
  end if

 end subroutine gemm_nonlop
!***

!----------------------------------------------------------------------
#if defined HAVE_GPU_CUDA

 !!****f* m_gemm_nonlop/init_gemm_nonlop_gpu
 !! NAME
 !! init_gemm_nonlop_gpu
 !!
 !! FUNCTION
 !! Memory allocation of the gemm_nonlop_kpt_gpu array
 !!
 !! INPUTS
 !! nkpt= number of k-points
 !!
 !! PARENTS
 !!      m_gstate
 !!
 !! CHILDREN
 !!      abi_zgemm_2r,dgemm,opernlc_ylm,xmpi_sum
 !!
 !! SOURCE
 subroutine init_gemm_nonlop_gpu(nkpt)

   integer,intent(in) :: nkpt
   integer :: ikpt

   ! *************************************************************************

   ! TODO only allocate the number of kpt treated by this proc
   ABI_MALLOC(gemm_nonlop_kpt_gpu, (nkpt))
   do ikpt=1,nkpt
     gemm_nonlop_kpt_gpu(ikpt)%npw = -1
     gemm_nonlop_kpt_gpu(ikpt)%nprojs = -1
   end do

 end subroutine init_gemm_nonlop_gpu
 !!***

 !!****f* m_gemm_nonlop/destroy_gemm_nonlop_gpu
 !! NAME
 !! destroy_gemm_nonlop_gpu
 !!
 !! FUNCTION
 !! Initalization of the gemm_nonlop_kpt array
 !!
 !! INPUTS
 !! nkpt= number of k-points
 !!
 !! PARENTS
 !!      m_gstate
 !!
 !! CHILDREN
 !!      abi_zgemm_2r,dgemm,opernlc_ylm,xmpi_sum
 !!
 !! SOURCE
 subroutine destroy_gemm_nonlop_gpu(nkpt)

   integer,intent(in) :: nkpt
   integer :: ikpt

   ! *************************************************************************

   ! TODO add cycling if kpt parallelism

   ! deallocate GPU ressource for each k point
   do ikpt = 1,nkpt
     if(gemm_nonlop_kpt_gpu(ikpt)%nprojs /= -1) then
       ! deallocate arrays projs, projs_r and projs_i
       !call gemm_nonlop_gpu_deallocate(gemm_nonlop_kpt_gpu(ikpt))
       call dealloc_on_gpu(gemm_nonlop_kpt_gpu(ikpt)%projs)
       call dealloc_on_gpu(gemm_nonlop_kpt_gpu(ikpt)%projs_r)
       call dealloc_on_gpu(gemm_nonlop_kpt_gpu(ikpt)%projs_i)
       gemm_nonlop_kpt_gpu(ikpt)%nprojs = -1
     end if
   end do

   ABI_FREE(gemm_nonlop_kpt_gpu)

 end subroutine destroy_gemm_nonlop_gpu
 !!***

 !!****f* m_gemm_nonlop/gemm_nonlop_gpu
!! NAME
!! gemm_nonlop_gpu
!!
!! FUNCTION
!! Replacement of nonlop.
!!
!! INPUTS
!!
!! PARENTS
!!      m_nonlop
!!
!! CHILDREN
!!      abi_zgemm_2r,dgemm,opernlc_ylm,xmpi_sum
!!
!! SOURCE
 subroutine gemm_nonlop_gpu(atindx1,choice,cpopt,cprjin,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,&
   &                        enl,indlmn,istwf_k,&
   &                        lambda,lmnmax,matblk,&
   &                        mpi_enreg,natom,nattyp,ndat,nkpgin,nkpgout,&
   &                        nnlout,npwin,npwout,nspinor,nspinortot,ntypat,paw_opt,&
   &                        sij,svectout,&
   &                        useylm,vectin,vectout,&
   &                        use_gpu_cuda)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: choice,cpopt,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout
  integer,intent(in) :: istwf_k,lmnmax,matblk,natom,ndat,nkpgin
  integer,intent(in) :: nkpgout,nnlout,npwin,npwout,nspinor,nspinortot,ntypat
  integer,intent(in) :: paw_opt,useylm
  integer,optional,intent(in)      :: use_gpu_cuda
  real(dp), target, intent(in)     :: lambda(ndat)
  type(pawcprj_type),intent(inout) :: cprjin(natom,nspinor*((cpopt+5)/5)*ndat)
  type(MPI_type),   intent(in)     :: mpi_enreg

  !arrays
  integer,  target, intent(in)     :: atindx1(natom)
  integer,  target, intent(in)     :: indlmn(6,lmnmax,ntypat)
  integer,          intent(in)     :: nattyp(ntypat)
  real(dp), target, intent(in)     :: enl(dimenl1,dimenl2,nspinortot**2,dimekbq)
  real(dp), target, intent(in)     :: sij(dimenl1,ntypat*((paw_opt+1)/3))
  real(dp), target, intent(inout)  ::  vectin (2,npwin*nspinor*ndat)
  real(dp), target, intent(out)    :: svectout(2,npwout*nspinor*(paw_opt/3)*ndat)
  real(dp), target, intent(inout)  ::  vectout(2,npwout*nspinor*ndat) !vz_i

  ! locals
  integer :: idat, nprojs, shift, iatom, nlmn, ierr, ibeg, iend
  integer :: cplex, cplex_enl, cplex_fac
  integer :: iatm, ndgxdt, ndgxdtfac, nd2gxdt, nd2gxdtfac, optder, itypat, ilmn
  integer :: cplex_dgxdt(1), cplex_d2gxdt(1)
  real(dp) :: dgxdt_dum_in(1,1,1,1,1), dgxdt_dum_out(1,1,1,1,1),dgxdt_dum_out2(1,1,1,1,1)
  real(dp) :: d2gxdt_dum_in(1,1,1,1,1), d2gxdt_dum_out(1,1,1,1,1),d2gxdt_dum_out2(1,1,1,1,1)
  integer :: npw_max
  integer :: nattyp_max

  real(dp), allocatable, target :: sij_typ(:)

  !type(c_ptr)                      :: projections_gpu,        s_projections_gpu,        vnl_projections_gpu
  real(dp),    allocatable, target :: projections_cpu(:,:,:), s_projections_cpu(:,:,:), vnl_projections_cpu(:,:,:)

  ! used inside opernlc_ylm_allwf_kokkos_cpp when iphase > 1
  type(c_ptr)                      :: vnl_projections2_gpu

  type(c_ptr)                      :: temp_realvec_gpu

  ! GPU waveform data are allocated in m_alloc_hamilt_gpu
  !type(c_ptr)                      :: vectin_gpu, vectout_gpu, svectout_gpu

  type(c_ptr)                      :: enl_gpu
  integer                          :: enl_size_bytes

  integer                          :: sizeof_int

  type(c_ptr)                      :: atindx1_gpu
  integer                          :: atindx1_size_bytes

  type(c_ptr)                      :: indlmn_gpu
  integer                          :: indlmn_size_bytes

  type(c_ptr)                      :: lambda_gpu
  integer                          :: lambda_size_bytes

  type(c_ptr)                      :: sij_typ_gpu
  integer                          :: sij_typ_size_bytes

  integer(kind=c_int32_t), parameter :: izero = 0
  integer(kind=c_int32_t), parameter :: fix_realvec_divide_by_2 = 0
  integer(kind=c_int32_t), parameter :: fix_realvec_zero_out    = 1

  integer                          :: shift_spinor

! *************************************************************************

  npw_max = MAX(npwin, npwout)

  cplex = 2; if (istwf_k>1) cplex=1
  cplex_enl = 1; if (paw_opt>0) cplex_enl=2*dimenl1/(lmnmax*(lmnmax+1)) ! is enl complex?
  cplex_fac = max(cplex,dimekbq)
  if ((nspinortot==2.or.cplex_enl==2).and.paw_opt>0.and.choice/=7) cplex_fac=2 ! is vnl_projections complex?

  nprojs = gemm_nonlop_kpt_gpu(gemm_nonlop_ikpt_this_proc_being_treated)%nprojs

  ! These will store the non-local factors for vectin, svectout and vectout respectively
  call gpu_memset(gemm_nonlop_kokkos%    projections_gpu, izero, cplex *     nprojs * nspinor*ndat * dp)
  call gpu_memset(gemm_nonlop_kokkos%  s_projections_gpu, izero, cplex *     nprojs * nspinor*ndat * dp)
  call gpu_memset(gemm_nonlop_kokkos%vnl_projections_gpu, izero, cplex_fac * nprojs * nspinor*ndat * dp)

  if (dimekbq>1) then
    ABI_CHECK(cplex_fac==2,"BUG: invalid cplex_fac==1 when dimekbq=2!")
    ABI_MALLOC_CUDA(vnl_projections2_gpu, (cplex_fac * nprojs * nspinor*ndat * dp))
    call gpu_memset(vnl_projections2_gpu, izero, cplex_fac * nprojs * nspinor*ndat * dp)
  end if

  !ABI_MALLOC_CUDA( vectin_gpu,  2 * npwin *nspinor*ndat * dp)
  !ABI_MALLOC_CUDA( vectout_gpu, 2 * npwout*nspinor*ndat * dp)
  !ABI_MALLOC_CUDA(svectout_gpu, 2 * npwout*nspinor*ndat*(paw_opt/3) * dp)

  call copy_on_gpu(C_LOC(vectin(1,1)), gemm_nonlop_kokkos%vectin_gpu, 2*npwin*nspinor*ndat*dp)

  !! gpu alloc and init : enl_gpu
  enl_size_bytes = dimenl1 * dimenl2 * nspinortot**2 * dimekbq * dp
  ABI_MALLOC_CUDA( enl_gpu, enl_size_bytes )
  call copy_on_gpu( C_LOC(enl(1,1,1,1)) , enl_gpu, enl_size_bytes )

  !! gpu alloc and init atindx1_gpu
  sizeof_int = 4
  atindx1_size_bytes = natom * sizeof_int
  ABI_MALLOC_CUDA( atindx1_gpu,  atindx1_size_bytes )
  call copy_on_gpu( C_LOC(atindx1(1)) , atindx1_gpu, atindx1_size_bytes )

  !! gpu alloc and init indlmn_gpu
  indlmn_size_bytes = 6*lmnmax*ntypat * sizeof_int
  ABI_MALLOC_CUDA( indlmn_gpu,  indlmn_size_bytes )
  call copy_on_gpu( C_LOC(indlmn(1,1,1)) , indlmn_gpu, indlmn_size_bytes )

  !! gpu alloc and init lambda_gpu
  lambda_size_bytes =  ndat * dp
  ABI_MALLOC_CUDA( lambda_gpu, lambda_size_bytes )
  call copy_on_gpu( C_LOC(lambda(1)) , lambda_gpu, lambda_size_bytes )

  if(nprojs == 0) then
    ! TODO check if this is correct
    vectout = zero
    if(paw_opt>0) svectout = vectin
    return
  end if

  ! determine precisely when temp_realvec needs to be allocated
  ! to factorize allocate (resp. deallocate) at the begining (resp. at the end) of subroutine
  ! to avoid multiple allocate/deallocate that can be costly
  if (cplex /= 2) then
    if ( (cpopt < 2) .or. &
      &  (paw_opt == 3 .or. paw_opt == 4) .or. &
      &  (paw_opt == 0 .or. paw_opt == 1 .or. paw_opt == 4)) then
      ABI_MALLOC_CUDA(temp_realvec_gpu, (npw_max * nspinor * ndat * dp))
    end if
  end if

  ABI_MALLOC(    projections_cpu,(cplex,     nprojs,nspinor*ndat))
  ABI_MALLOC(  s_projections_cpu,(cplex,     nprojs,nspinor*ndat)) ! TODO - TO BE REMOVED ONCE CUDA-IZATION IS OK
  ABI_MALLOC(vnl_projections_cpu,(cplex_fac, nprojs,nspinor*ndat)) ! TODO - TO BE REMOVED ONCE CUDA-IZATION IS OK
  projections_cpu = zero
  s_projections_cpu = zero
  vnl_projections_cpu = zero

  if(cpopt >= 2) then

    ! retrieve from cprjin
    do idat=1, ndat*nspinor
      shift = 0
      do iatom = 1, natom
        nlmn = cprjin(iatom, idat)%nlmn
        projections_cpu(1:cplex, shift+1:shift+nlmn, idat) = cprjin(iatom, idat)%cp(1:cplex, 1:nlmn)
        shift = shift + nlmn
      end do
    end do

    ! copy from HOST projections_cpu to GPU projections_gpu
    call copy_on_gpu(C_LOC(projections_cpu(1,1,1)), gemm_nonlop_kokkos%projections_gpu, cplex * nprojs * nspinor*ndat * dp)

  else ! cpopt < 2

     ! opernla
     if(cplex == 2) then

       ! projections_gpu = projs * vectin_gpu
       call gpu_xgemm(cplex, 'C', 'N', &
         &            nprojs, ndat*nspinor, npwin, &                                                ! M,N,K
         &            cone, &                                                                       ! alpha
         &            gemm_nonlop_kpt_gpu(gemm_nonlop_ikpt_this_proc_being_treated)%projs, npwin, & ! A, LDA
         &            gemm_nonlop_kokkos%vectin_gpu, npwin, &                                                          ! B, LDB
         &            czero, &                                                                      ! beta
         &            gemm_nonlop_kokkos%projections_gpu, nprojs)                                   ! C, LDC

     else

       if (.not. gpu_allocated(temp_realvec_gpu)) then
         ABI_ERROR("Please provide memory allocation for temp_realvec_gpu array")
       end if


       ! only compute real part of projections = P^* psi => projections_r = P_r^T psi_r + P_i^T psi_i
       !temp_realvec(1:npwin*nspinor*ndat) = vectin(1,1:npwin*nspinor*ndat)
       call extract_real_part(temp_realvec_gpu, gemm_nonlop_kokkos%vectin_gpu, npwin*nspinor*ndat)

       if(istwf_k == 2 .and. mpi_enreg%me_g0 == 1) then
         ! do idat=1, ndat*nspinor
         !   temp_realvec(1+(idat-1)*npwin) = temp_realvec(1+(idat-1)*npwin)/2
         ! end do
         call fix_realvec(temp_realvec_gpu, npwin, ndat*nspinor, fix_realvec_divide_by_2)
       end if

       call gpu_xgemm(cplex, 'T', 'N', &
         &            nprojs, ndat*nspinor, npwin, &                                                  ! M,N,K
         &            cone, &                                                                         ! alpha
         &            gemm_nonlop_kpt_gpu(gemm_nonlop_ikpt_this_proc_being_treated)%projs_r, npwin, & ! A, LDA
         &            temp_realvec_gpu, npwin, &                                                      ! B, LDB
         &            czero, &                                                                        ! beta
         &            gemm_nonlop_kokkos%projections_gpu, nprojs)                                     ! C, LDC

       !temp_realvec(1:npwin*nspinor*ndat) = vectin(2,1:npwin*nspinor*ndat)
       call extract_imag_part(temp_realvec_gpu, gemm_nonlop_kokkos%vectin_gpu, npwin*nspinor*ndat)

       if(istwf_k == 2 .and. mpi_enreg%me_g0 == 1) then
         ! do idat=1, ndat*nspinor
         !   temp_realvec(1+(idat-1)*npwin) = zero
         ! end do
         call fix_realvec(temp_realvec_gpu, npwin, ndat*nspinor, fix_realvec_zero_out)
       end if
       call gpu_xgemm(cplex, 'T', 'N', &
         &            nprojs, ndat*nspinor, npwin, &
         &            cone, &
         &            gemm_nonlop_kpt_gpu(gemm_nonlop_ikpt_this_proc_being_treated)%projs_i, npwin, &
         &            temp_realvec_gpu, npwin, &
         &            cone, &
         &            gemm_nonlop_kokkos%projections_gpu, nprojs)

       !gemm_nonlop_kokkos%projections_gpu = 2 * gemm_nonlop_kokkos%projections_gpu
       call gpu_xscal(cplex, nprojs*nspinor*ndat, ctwo, gemm_nonlop_kokkos%projections_gpu, 1)

     end if ! cplex == 2

!     call xmpi_sum(projections,mpi_enreg%comm_fft,ierr)

    call copy_from_gpu(C_LOC(  projections_cpu(1,1,1)),   gemm_nonlop_kokkos%projections_gpu, cplex * nprojs * nspinor*ndat * dp)

    if(cpopt >= 0) then
      ! copy from GPU projections_gpu to HOST projections_cpu
      call copy_from_gpu(C_LOC(projections_cpu(1,1,1)), gemm_nonlop_kokkos%projections_gpu, cplex * nprojs * nspinor*ndat * dp)

      ! store in cprjin
      do idat=1, ndat*nspinor
        shift = 0
        do iatom = 1, natom
          nlmn = cprjin(iatom, idat)%nlmn
          cprjin(iatom, idat)%cp(1:cplex, 1:nlmn) = projections_cpu(1:cplex, shift+1:shift+nlmn, idat)
          shift = shift + nlmn
        end do
      end do

    end if ! cpopt >= 0

   end if ! cpopt >= 2

    if(choice > 0) then

      if(choice /= 7) then
        ! opernlc
        iatm = 0
        ndgxdt = 0
        ndgxdtfac = 0
        nd2gxdt = 0
        nd2gxdtfac = 0
        optder = 0

        ! get the maximun of nattyp array
        nattyp_max = maxval(nattyp)

        !! gpu alloc and init sij_typ_size_bytes
        sij_typ_size_bytes = (((paw_opt+1)/3)*lmnmax*(lmnmax+1)/2) * dp
        ABI_MALLOC_CUDA( sij_typ_gpu, sij_typ_size_bytes )

        ABI_MALLOC     ( sij_typ    , (((paw_opt+1)/3)*lmnmax*(lmnmax+1)/2) )

        shift = 0
        do itypat=1, ntypat
          nlmn=count(indlmn(3,:,itypat)>0)
          if (paw_opt>=2) then

            if (cplex_enl==1) then

              do ilmn=1,nlmn*(nlmn+1)/2
                sij_typ(ilmn) = sij(ilmn,itypat)
              end do

            else

              do ilmn=1,nlmn*(nlmn+1)/2
                sij_typ(ilmn) = sij(2*ilmn-1,itypat)
              end do

            end if

            call copy_on_gpu( C_LOC(sij_typ(1)) , sij_typ_gpu, sij_typ_size_bytes )

          end if ! paw_opt>=2

          ibeg = shift+1
          iend = shift+nattyp(itypat)*nlmn

          !   ! TODO - PK - old version - to be removed soon
          ! do idat = 1,ndat
          !   call opernlc_ylm_gpu(atindx1, cplex, cplex_enl, cplex_fac, &
          !     &                  dimenl1, dimenl2, dimekbq, enl, &
          !     &                  projections_cpu(:, ibeg:iend, 1+nspinor*(idat-1):nspinor*idat), &
          !     &                  vnl_projections_cpu(:, ibeg:iend,1+nspinor*(idat-1):nspinor*idat), &
          !     &                  s_projections_cpu(:, ibeg:iend,1+nspinor*(idat-1):nspinor*idat), &
          !     &                  iatm, indlmn(:,:,itypat), itypat, lambda(idat), mpi_enreg, natom, &
          !     &                  nattyp(itypat), nlmn, nspinor, nspinortot, paw_opt, sij_typ)
          ! end do ! idat

#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS)
          ! TODO - PK

          !Parallelization over spinors treatment
          shift_spinor = 0
          if (mpi_enreg%paral_spinor==1) then
            shift_spinor = mpi_enreg%me_spinor
          end if

          if (gemm_nonlop_use_kokkos) then
            !write(*,*) "opernlc_ylm_allwf_kokkos"
            call opernlc_ylm_allwf_kokkos(cplex, cplex_enl, cplex_fac, &
              &                           dimenl1, dimenl2, dimekbq, &
              &                           iatm, itypat, ntypat, nprojs, &
              &                           natom, nattyp(itypat), nspinor, &
              &                           nspinortot, paw_opt, &
              &                           nlmn, lmnmax, &
              &                           enl_gpu, &
              &                           gemm_nonlop_kokkos%projections_gpu, &
              &                           gemm_nonlop_kokkos%vnl_projections_gpu, &
              &                           vnl_projections2_gpu, &
              &                           gemm_nonlop_kokkos%s_projections_gpu, &
              &                           shift_spinor, ndat, &
              &                           atindx1_gpu, &
              &                           indlmn_gpu, &
              &                           lambda_gpu, &
              &                           sij_typ_gpu, &
              &                           shift, nattyp_max)
          else
            ! TODO - PK
            !write(*,*) "opernlc_ylm_allwf_cpu"
            call opernlc_ylm_allwf_cpu(atindx1, cplex, cplex_enl, cplex_fac, &
              &                  dimenl1, dimenl2, dimekbq, enl, &
              &                  projections_cpu(:, ibeg:iend, 1:nspinor*ndat), &
              &                  vnl_projections_cpu(:, ibeg:iend, 1:nspinor*ndat), &
              &                  s_projections_cpu(:, ibeg:iend, 1:nspinor*ndat), &
              &                  iatm, indlmn(:,:,itypat), itypat, ndat, lambda, mpi_enreg, natom, &
              &                  nattyp(itypat), nlmn, nspinor, nspinortot, paw_opt, sij_typ)
          end if

#else

          ! TODO - PK
          call opernlc_ylm_allwf_cpu(atindx1, cplex, cplex_enl, cplex_fac, &
            &                  dimenl1, dimenl2, dimekbq, enl, &
            &                  projections_cpu(:, ibeg:iend, 1:nspinor*ndat), &
            &                  vnl_projections_cpu(:, ibeg:iend, 1:nspinor*ndat), &
            &                  s_projections_cpu(:, ibeg:iend, 1:nspinor*ndat), &
            &                  iatm, indlmn(:,:,itypat), itypat, ndat, lambda, mpi_enreg, natom, &
            &                  nattyp(itypat), nlmn, nspinor, nspinortot, paw_opt, sij_typ)

#endif

          shift = shift + nattyp(itypat)*nlmn
          iatm  = iatm  + nattyp(itypat)
        end do ! itypat
        ABI_FREE(sij_typ)
        ABI_FREE_CUDA( sij_typ_gpu )

        if (gemm_nonlop_use_kokkos) then
          ! nothing to do, s_projections and vnl_projections data are already in GPU memory
        else
          ! TO BE REMOVED LATTER
          ! upload s_projections and vnl_projections to GPU
          call copy_on_gpu(C_LOC(  s_projections_cpu(1,1,1)), gemm_nonlop_kokkos%  s_projections_gpu, cplex     * nprojs * nspinor*ndat * dp)
          call copy_on_gpu(C_LOC(vnl_projections_cpu(1,1,1)), gemm_nonlop_kokkos%vnl_projections_gpu, cplex_fac * nprojs * nspinor*ndat * dp)
        end if


      else ! choice == 7

        ! TO BE REMOVED - DEBUG ONLY
        !s_projections_cpu = projections_cpu
        call copy_gpu_to_gpu(gemm_nonlop_kokkos%s_projections_gpu, &
          &                  gemm_nonlop_kokkos%projections_gpu, cplex * nprojs * nspinor*ndat * dp)

      end if ! choice

      ! opernlb
      if (paw_opt == 3 .or. paw_opt == 4) then

        ! Get svectout from s_projections
        if(cplex == 2) then

          call gpu_xgemm(cplex, 'N', 'N', &
            &            npwout, ndat*nspinor, nprojs, &                                                ! M,N,K
            &            cone, &                                                                        ! alpha
            &            gemm_nonlop_kpt_gpu(gemm_nonlop_ikpt_this_proc_being_treated)%projs, npwout, & ! A, LDA
            &            gemm_nonlop_kokkos%s_projections_gpu, nprojs, &                                ! B, LDB
            &            czero, &                                                                       ! beta
            &            gemm_nonlop_kokkos%svectout_gpu, npwout)                                       ! C, LDC

        else

          if (.not. gpu_allocated(temp_realvec_gpu)) then
            ABI_ERROR("Please provide memory allocation for temp_realvec_gpu array")
          end if

          call gpu_xgemm(cplex, 'N', 'N', &
            &            npwout, ndat*nspinor, nprojs, &                                                ! M,N,K
            &            cone, &                                                                        ! alpha
            &            gemm_nonlop_kpt_gpu(gemm_nonlop_ikpt_this_proc_being_treated)%projs_r, npwout, & ! A, LDA
            &            gemm_nonlop_kokkos%s_projections_gpu, nprojs, &                                ! B, LDB
            &            czero, &                                                                       ! beta
            &            temp_realvec_gpu, npwout)                                                      ! C, LDC
          !svectout(1,1:npwout*nspinor*ndat) = temp_realvec_gpu(1:npwout*nspinor*ndat)
          call insert_real_part(gemm_nonlop_kokkos%svectout_gpu, temp_realvec_gpu, npwout*nspinor*ndat)

          call gpu_xgemm(cplex, 'N', 'N', &
            &            npwout, ndat*nspinor, nprojs, &                                                ! M,N,K
            &            cone, &                                                                        ! alpha
            &            gemm_nonlop_kpt_gpu(gemm_nonlop_ikpt_this_proc_being_treated)%projs_i, npwout,& ! A, LDA
            &            gemm_nonlop_kokkos%s_projections_gpu, nprojs, &                                ! B, LDB
            &            czero, &                                                                       ! beta
            &            temp_realvec_gpu, npwout)                                                      ! C, LDC
          !svectout(2,1:npwout*nspinor*ndat) = temp_realvec(1:npwout*nspinor*ndat)
          call insert_imag_part(gemm_nonlop_kokkos%svectout_gpu, temp_realvec_gpu, npwout*nspinor*ndat)

        end if ! cplex == 2

        if (choice /= 7) then

          ! compute : svectout_gpu = svectout_gpu + vectin_gpu
          ! this a axpy operation with x => vectin_gpu, y => svectout_gpu and alpha=1
          ! please remember that svectout_gpu and vectin_gpu have same size when paw_opt >= 3 and paw_opt<6
          ! this is the case here
          call gpu_xaxpy(1, &                                  ! real
            &            2*npwin*nspinor*ndat, &               ! size
            &            cone, &                               ! alpha
            &            gemm_nonlop_kokkos%vectin_gpu, 1, &   ! X, incrx
            &            gemm_nonlop_kokkos%svectout_gpu, 1)   ! Y, incry

        endif

        ! copy back results on host
        call copy_from_gpu(C_LOC(svectout(1,1)), gemm_nonlop_kokkos%svectout_gpu,&
          & 2*npwout*nspinor*(paw_opt/3)*ndat * dp)

      end if ! (paw_opt == 3 .or. paw_opt == 4)

      if (paw_opt == 0 .or. paw_opt == 1 .or. paw_opt == 4) then

        ! Get vectout from vnl_projections
        if (cplex_fac == 2) then

          call gpu_xgemm(cplex, 'N', 'N', &
            &            npwout, ndat*nspinor, nprojs, &                                                ! M,N,K
            &            cone, &                                                                        ! alpha
            &            gemm_nonlop_kpt_gpu(gemm_nonlop_ikpt_this_proc_being_treated)%projs, npwout, & ! A, LDA
            &            gemm_nonlop_kokkos%vnl_projections_gpu, nprojs, &                              ! B, LDB
            &            czero, &                                                                       ! beta
            &            gemm_nonlop_kokkos%vectout_gpu, npwout)                                        ! C, LDC

        else

          if (.not. gpu_allocated(temp_realvec_gpu)) then
            ABI_ERROR("Please provide memory allocation for temp_realvec_gpu array")
          end if

          call gpu_xgemm(cplex, 'N', 'N', &
            &            npwout, ndat*nspinor, nprojs, &                                               ! M,N,K
            &            cone, &                                                                       ! alpha
            &            gemm_nonlop_kpt_gpu(gemm_nonlop_ikpt_this_proc_being_treated)%projs_r, npwout, &  ! A, LDA
            &            gemm_nonlop_kokkos%vnl_projections_gpu, nprojs, &                             ! B, LDB
            &            czero, &                                                                      ! beta
            &            temp_realvec_gpu, npwout)                                                     ! C, LDC
          ! vectout(1,1:npwout*nspinor*ndat) = temp_realvec(1:npwout*nspinor*ndat)
          call insert_real_part(gemm_nonlop_kokkos%vectout_gpu, temp_realvec_gpu, npwout*nspinor*ndat)

          call gpu_xgemm(cplex, 'N', 'N', &
            &            npwout, ndat*nspinor, nprojs, &                                               ! M,N,K
            &            cone, &                                                                       ! alpha
            &            gemm_nonlop_kpt_gpu(gemm_nonlop_ikpt_this_proc_being_treated)%projs_i, npwout, & ! A, LDA
            &            gemm_nonlop_kokkos%vnl_projections_gpu, nprojs, &                             ! B, LDB
            &            czero, &                                                                      ! beta
            &            temp_realvec_gpu, npwout)                                                     ! C, LDC
          ! vectout(2,1:npwout*nspinor*ndat) = temp_realvec(1:npwout*nspinor*ndat)
          call insert_imag_part(gemm_nonlop_kokkos%vectout_gpu, temp_realvec_gpu, npwout*nspinor*ndat)

        end if  ! cplex_fac == 2

        ! copy back results on host
        call copy_from_gpu(C_LOC(vectout(1,1)), gemm_nonlop_kokkos%vectout_gpu, 2*npwout*nspinor*ndat * dp)

      end if ! (paw_opt == 0 .or. paw_opt == 1 .or. paw_opt == 4)

    end if ! choice > 0

  if (dimekbq>1) then
    ABI_FREE_CUDA(vnl_projections2_gpu)
  end if

  !ABI_FREE_CUDA(      vectin_gpu)
  !ABI_FREE_CUDA(     vectout_gpu)
  !ABI_FREE_CUDA(    svectout_gpu)

  if (gpu_allocated(temp_realvec_gpu)) then
    ABI_FREE_CUDA(temp_realvec_gpu)
  end if

  ABI_FREE_CUDA( enl_gpu )
  ABI_FREE_CUDA( atindx1_gpu )
  ABI_FREE_CUDA( indlmn_gpu )
  ABI_FREE_CUDA( lambda_gpu )

  ! if projections_cpu was allocated, then free it here
  ABI_FREE(    projections_cpu)
  ABI_FREE(  s_projections_cpu) ! TO BE REMOVED
  ABI_FREE(vnl_projections_cpu) ! TO BE REMOVED


 end subroutine gemm_nonlop_gpu
!***

#endif

end module m_gemm_nonlop
!!***
