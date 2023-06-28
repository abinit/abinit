!!****m* ABINIT/m_gemm_nonlop_gpu
!! NAME
!! m_gemm_nonlop_gpu
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


#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_gemm_nonlop_gpu

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_fstrings,    only : itoa, ftoa, sjoin

 use m_abi_linalg  ! copy_on_gpu, copy_from_gpu, alloc_on_gpu, dealloc_on_gpu, gpu_memset, gpu_allocated
 use defs_abitypes, only : MPI_type
 use m_opernlc_ylm_allwf_cpu, only : opernlc_ylm_allwf_cpu
 use m_pawcprj, only : pawcprj_type
 use m_gemm_nonlop, only : gemm_nonlop_type,gemm_nonlop_ikpt_this_proc_being_treated,make_gemm_nonlop,gemm_nonlop_kpt

#if defined(HAVE_GPU_CUDA)
  use m_gpu_toolbox
  use m_alloc_hamilt_gpu, only : gemm_nonlop_kokkos
#endif

#ifdef HAVE_KOKKOS
 use m_manage_kokkos, only : opernlc_ylm_allwf_kokkos
#endif

#ifdef HAVE_FC_ISO_C_BINDING
 use, intrinsic :: iso_c_binding, only : c_ptr, c_int32_t, c_int64_t, c_float, c_double, c_size_t, c_loc
#endif

 implicit none

 private

 public :: init_gemm_nonlop_gpu
 public :: destroy_gemm_nonlop_gpu
 public :: make_gemm_nonlop_gpu
 public :: gemm_nonlop_gpu

 !!***

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

 end interface

#endif

 contains

!----------------------------------------------------------------------

#if defined HAVE_GPU_CUDA

 !!****f* m_gemm_nonlop_gpu/init_gemm_nonlop_gpu
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

#if defined(HAVE_KOKKOS)
   gemm_nonlop_kokkos % allocated = .false.
#endif

 end subroutine init_gemm_nonlop_gpu
 !!***

 !!****f* m_gemm_nonlop_gpu/destroy_gemm_nonlop_gpu
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

   ! This function must be called before destroy_gemm_nonlop so it can
   ! properly figure out which GPU buffer to free.
   if(.not. allocated(gemm_nonlop_kpt)) then
     ABI_BUG("gemm_nonlop is already free, cannot free GPU resources !")
   end if

   ! deallocate GPU ressource for each k point
   do ikpt = 1,nkpt
     if(gemm_nonlop_kpt_gpu(ikpt)%nprojs /= -1) then
       ! deallocate arrays projs, projs_r and projs_i
       if (allocated(gemm_nonlop_kpt(ikpt)%projs)) then
         call dealloc_on_gpu(gemm_nonlop_kpt_gpu(ikpt)%projs)
       end if
       if (allocated(gemm_nonlop_kpt(ikpt)%projs_r)) then
         call dealloc_on_gpu(gemm_nonlop_kpt_gpu(ikpt)%projs_r)
       end if
       if (allocated(gemm_nonlop_kpt(ikpt)%projs_i)) then
         call dealloc_on_gpu(gemm_nonlop_kpt_gpu(ikpt)%projs_i)
       end if
       gemm_nonlop_kpt_gpu(ikpt)%nprojs = -1
     end if
   end do

   ABI_FREE(gemm_nonlop_kpt_gpu)

 end subroutine destroy_gemm_nonlop_gpu
 !!***


!!****f* m_gemm_nonlop_gpu/make_gemm_nonlop_gpu
!! NAME
!! make_gemm_nonlop_gpu
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
 subroutine make_gemm_nonlop_gpu(ikpt,npw,lmnmax,ntypat,indlmn,nattyp,istwf_k,ucvol,ffnl_k, &
&                     ph3d_k,kpt_k,kg_k,kpg_k, &
&                     compute_grad_strain,compute_grad_atom) ! Optional parameters

  !Arguments ------------------------------------
  integer, intent(in) :: ikpt
  integer, intent(in) :: npw, lmnmax,ntypat
  integer, intent(in) :: indlmn(:,:,:), kg_k(:,:)
  integer, intent(in) :: nattyp(ntypat)
  integer, intent(in) :: istwf_k
  logical, intent(in), optional :: compute_grad_strain,compute_grad_atom
  real(dp), intent(in) :: ucvol
  real(dp), intent(in) :: ffnl_k(:,:,:,:)
  real(dp), intent(in) :: ph3d_k(:,:,:)
  real(dp), intent(in) :: kpt_k(:)
  real(dp), intent(in), target :: kpg_k(:,:)

  ! locals
  integer              :: nprojs, itypat

! *************************************************************************

  ABI_UNUSED((/ikpt,npw,lmnmax,ntypat,indlmn,kg_k,nattyp,istwf_k/))
  ABI_UNUSED((/ucvol,ffnl_k,ph3d_k,kpt_k,kpg_k/))
  ABI_UNUSED((/compute_grad_strain,compute_grad_atom/))

  call make_gemm_nonlop(ikpt,npw,lmnmax,ntypat,indlmn,nattyp,istwf_k,ucvol,ffnl_k, &
         ph3d_k,kpt_k,kg_k,kpg_k, &
         compute_grad_strain=compute_grad_strain,compute_grad_atom=compute_grad_atom)

  nprojs = 0
  do itypat=1,ntypat
    nprojs = nprojs + count(indlmn(3,:,itypat)>0)*nattyp(itypat)
  end do
  gemm_nonlop_kpt_gpu(ikpt)%npw    = npw
  gemm_nonlop_kpt_gpu(ikpt)%nprojs = nprojs

#ifdef DEBUG_VERBOSE_GPU
  if(xmpi_comm_rank(xmpi_world) == 0) then
    call check_gpu_mem("make_gemm_nonlop begin")
    call wrtout(std_out,sjoin(" npw    .......", itoa(npw)),    'COLL')
    call wrtout(std_out,sjoin(" nprojs .......", itoa(nprojs)), 'COLL')
  end if
#endif

  if(istwf_k <= 1) then
    call alloc_on_gpu(gemm_nonlop_kpt_gpu(ikpt)%projs, INT(2,c_size_t)*npw*nprojs*dp)
    ! TODO : gradients
  else
    call alloc_on_gpu(gemm_nonlop_kpt_gpu(ikpt)%projs_r, INT(1, c_size_t)*npw*nprojs*dp)
    call alloc_on_gpu(gemm_nonlop_kpt_gpu(ikpt)%projs_i, INT(1, c_size_t)*npw*nprojs*dp)
    ! TODO : gradients
  end if

#ifdef DEBUG_VERBOSE_GPU
  if(xmpi_comm_rank(xmpi_world) == 0) then
    call check_gpu_mem("make_gemm_nonlop end  ")
  end if
#endif

  ! upload data to gpu memory
  if(istwf_k <= 1) then
    call copy_on_gpu(c_loc(gemm_nonlop_kpt(ikpt)%projs(1,1,1)), gemm_nonlop_kpt_gpu(ikpt)%projs, INT(2, c_size_t)*npw*nprojs*dp)
    ! TODO : gradients
  else
    call copy_on_gpu(c_loc(gemm_nonlop_kpt(ikpt)%projs_r(1,1,1)), gemm_nonlop_kpt_gpu(ikpt)%projs_r, &
      &                    INT(1, c_size_t)*npw*nprojs*dp)
    call copy_on_gpu(c_loc(gemm_nonlop_kpt(ikpt)%projs_i(1,1,1)), gemm_nonlop_kpt_gpu(ikpt)%projs_i, &
      &                    INT(1, c_size_t)*npw*nprojs*dp)
  end if

 end subroutine make_gemm_nonlop_gpu
!!***

!!****f* m_gemm_nonlop_gpu/gemm_nonlop_gpu
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
  integer(c_size_t)                :: vectin_size

  type(c_ptr)                      :: enl_gpu
  integer(c_size_t)                :: enl_size_bytes

  integer                          :: sizeof_int

  type(c_ptr)                      :: atindx1_gpu
  integer(c_size_t)                :: atindx1_size_bytes

  type(c_ptr)                      :: indlmn_gpu
  integer(c_size_t)                :: indlmn_size_bytes

  type(c_ptr)                      :: lambda_gpu
  integer(c_size_t)                :: lambda_size_bytes

  type(c_ptr)                      :: sij_typ_gpu
  integer(c_size_t)                :: sij_typ_size_bytes

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
  call gpu_memset(gemm_nonlop_kokkos%    projections_gpu, izero, INT(cplex,     c_size_t) * nprojs * nspinor*ndat * dp)
  call gpu_memset(gemm_nonlop_kokkos%  s_projections_gpu, izero, INT(cplex,     c_size_t) * nprojs * nspinor*ndat * dp)
  call gpu_memset(gemm_nonlop_kokkos%vnl_projections_gpu, izero, INT(cplex_fac, c_size_t) * nprojs * nspinor*ndat * dp)

  if (dimekbq>1) then
    ABI_CHECK(cplex_fac==2,"BUG: invalid cplex_fac==1 when dimekbq=2!")
    ABI_MALLOC_CUDA(vnl_projections2_gpu,        INT(cplex_fac, c_size_t) * nprojs * nspinor*ndat * dp)
    call gpu_memset(vnl_projections2_gpu, izero, INT(cplex_fac, c_size_t) * nprojs * nspinor*ndat * dp)
  end if

  !ABI_MALLOC_CUDA( vectin_gpu,  INT(2, c_size_t) * npwin *nspinor*ndat * dp)
  !ABI_MALLOC_CUDA( vectout_gpu, INT(2, c_size_t) * npwout*nspinor*ndat * dp)
  !ABI_MALLOC_CUDA(svectout_gpu, INT(2, c_size_t) * npwout*nspinor*ndat*(paw_opt/3) * dp)

  !call copy_on_gpu(C_LOC(vectin(1,1)), gemm_nonlop_kokkos%vectin_gpu, INT(2, c_size_t)*npwin*nspinor*ndat*dp)
  vectin_size = 2*npwin*nspinor*ndat*dp
  call gpu_data_prefetch_async(C_LOC(vectin), vectin_size)

  if (choice == 7) then
    call gpu_data_prefetch_async(C_LOC(svectout), vectin_size)
  end if

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
      ABI_MALLOC_CUDA(temp_realvec_gpu, (INT(npw_max, c_size_t) * nspinor * ndat * dp))
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
    call copy_on_gpu(C_LOC(projections_cpu(1,1,1)), gemm_nonlop_kokkos%projections_gpu,&
        INT(cplex, c_size_t) * nprojs * nspinor*ndat * dp)

  else ! cpopt < 2

     ! opernla
     if(cplex == 2) then

       ! projections_gpu = projs * vectin_gpu
       call abi_gpu_xgemm(cplex, 'C', 'N', &
         &            nprojs, ndat*nspinor, npwin, &                                                ! M,N,K
         &            cone, &                                                                       ! alpha
         &            gemm_nonlop_kpt_gpu(gemm_nonlop_ikpt_this_proc_being_treated)%projs, npwin, & ! A, LDA
         &            C_LOC(vectin), npwin, &                                                       ! B, LDB
         &            czero, &                                                                      ! beta
         &            gemm_nonlop_kokkos%projections_gpu, nprojs)                                   ! C, LDC

     else

       if (.not. gpu_allocated(temp_realvec_gpu)) then
         ABI_ERROR("Please provide memory allocation for temp_realvec_gpu array")
       end if


       ! only compute real part of projections = P^* psi => projections_r = P_r^T psi_r + P_i^T psi_i
       !temp_realvec(1:npwin*nspinor*ndat) = vectin(1,1:npwin*nspinor*ndat)
       call extract_real_part(temp_realvec_gpu, C_LOC(vectin), npwin*nspinor*ndat)

       if(istwf_k == 2 .and. mpi_enreg%me_g0 == 1) then
         ! do idat=1, ndat*nspinor
         !   temp_realvec(1+(idat-1)*npwin) = temp_realvec(1+(idat-1)*npwin)/2
         ! end do
         call fix_realvec(temp_realvec_gpu, npwin, ndat*nspinor, fix_realvec_divide_by_2)
       end if

       call abi_gpu_xgemm(cplex, 'T', 'N', &
         &            nprojs, ndat*nspinor, npwin, &                                                  ! M,N,K
         &            cone, &                                                                         ! alpha
         &            gemm_nonlop_kpt_gpu(gemm_nonlop_ikpt_this_proc_being_treated)%projs_r, npwin, & ! A, LDA
         &            temp_realvec_gpu, npwin, &                                                      ! B, LDB
         &            czero, &                                                                        ! beta
         &            gemm_nonlop_kokkos%projections_gpu, nprojs)                                     ! C, LDC

       !temp_realvec(1:npwin*nspinor*ndat) = vectin(2,1:npwin*nspinor*ndat)
       call extract_imag_part(temp_realvec_gpu, C_LOC(vectin), npwin*nspinor*ndat)

       if(istwf_k == 2 .and. mpi_enreg%me_g0 == 1) then
         ! do idat=1, ndat*nspinor
         !   temp_realvec(1+(idat-1)*npwin) = zero
         ! end do
         call fix_realvec(temp_realvec_gpu, npwin, ndat*nspinor, fix_realvec_zero_out)
       end if
       call abi_gpu_xgemm(cplex, 'T', 'N', &
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

    call copy_from_gpu(C_LOC(  projections_cpu(1,1,1)),   gemm_nonlop_kokkos%projections_gpu,&
        INT(cplex, c_size_t) * nprojs * nspinor*ndat * dp)

    if(cpopt >= 0) then
      ! copy from GPU projections_gpu to HOST projections_cpu
      call copy_from_gpu(C_LOC(projections_cpu(1,1,1)), gemm_nonlop_kokkos%projections_gpu,&
          INT(cplex, c_size_t) * nprojs * nspinor*ndat * dp)

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
          call copy_on_gpu(C_LOC(  s_projections_cpu(1,1,1)), gemm_nonlop_kokkos%  s_projections_gpu, &
            &              INT(cplex,     c_size_t) * nprojs * nspinor*ndat * dp)
          call copy_on_gpu(C_LOC(vnl_projections_cpu(1,1,1)), gemm_nonlop_kokkos%vnl_projections_gpu, &
            &              INT(cplex_fac, c_size_t) * nprojs * nspinor*ndat * dp)
        end if


      else ! choice == 7

        ! TO BE REMOVED - DEBUG ONLY
        !s_projections_cpu = projections_cpu
        call copy_gpu_to_gpu(gemm_nonlop_kokkos%s_projections_gpu, &
             &               gemm_nonlop_kokkos%projections_gpu, &
             &               INT(cplex, c_size_t) * nprojs * nspinor * ndat * dp)

      end if ! choice

      ! opernlb
      if (paw_opt == 3 .or. paw_opt == 4) then

        ! Get svectout from s_projections
        if(cplex == 2) then

          call abi_gpu_xgemm(cplex, 'N', 'N', &
            &            npwout, ndat*nspinor, nprojs, &                                                ! M,N,K
            &            cone, &                                                                        ! alpha
            &            gemm_nonlop_kpt_gpu(gemm_nonlop_ikpt_this_proc_being_treated)%projs, npwout, & ! A, LDA
            &            gemm_nonlop_kokkos%s_projections_gpu, nprojs, &                                ! B, LDB
            &            czero, &                                                                       ! beta
            &            C_LOC(svectout(1,1)), npwout)                                                  ! C, LDC

        else

          if (.not. gpu_allocated(temp_realvec_gpu)) then
            ABI_ERROR("Please provide memory allocation for temp_realvec_gpu array")
          end if

          call abi_gpu_xgemm(cplex, 'N', 'N', &
            &            npwout, ndat*nspinor, nprojs, &                                                ! M,N,K
            &            cone, &                                                                        ! alpha
            &            gemm_nonlop_kpt_gpu(gemm_nonlop_ikpt_this_proc_being_treated)%projs_r, npwout, & ! A, LDA
            &            gemm_nonlop_kokkos%s_projections_gpu, nprojs, &                                ! B, LDB
            &            czero, &                                                                       ! beta
            &            temp_realvec_gpu, npwout)                                                      ! C, LDC
          !svectout(1,1:npwout*nspinor*ndat) = temp_realvec_gpu(1:npwout*nspinor*ndat)
          call insert_real_part(C_LOC(svectout(1,1)), temp_realvec_gpu, npwout*nspinor*ndat)

          call abi_gpu_xgemm(cplex, 'N', 'N', &
            &            npwout, ndat*nspinor, nprojs, &                                                ! M,N,K
            &            cone, &                                                                        ! alpha
            &            gemm_nonlop_kpt_gpu(gemm_nonlop_ikpt_this_proc_being_treated)%projs_i, npwout,& ! A, LDA
            &            gemm_nonlop_kokkos%s_projections_gpu, nprojs, &                                ! B, LDB
            &            czero, &                                                                       ! beta
            &            temp_realvec_gpu, npwout)                                                      ! C, LDC
          !svectout(2,1:npwout*nspinor*ndat) = temp_realvec(1:npwout*nspinor*ndat)
          call insert_imag_part(C_LOC(svectout(1,1)), temp_realvec_gpu, npwout*nspinor*ndat)

        end if ! cplex == 2

        if (choice /= 7) then

          ! compute : svectout_gpu = svectout_gpu + vectin_gpu
          ! this a axpy operation with x => vectin_gpu, y => svectout_gpu and alpha=1
          ! please remember that svectout_gpu and vectin_gpu have same size when paw_opt >= 3 and paw_opt<6
          ! this is the case here
          call abi_gpu_xaxpy(1, &                              ! real
            &            2*npwin*nspinor*ndat, &               ! size
            &            cone, &                               ! alpha
            &            C_LOC(vectin), 1, &                   ! X, incrx
            &            C_LOC(svectout), 1)                   ! Y, incry

        endif

        ! copy back results on host
        !call copy_from_gpu(C_LOC(svectout(1,1)), gemm_nonlop_kokkos%svectout_gpu,&
        !  & INT(2, c_size_t)*npwout*nspinor*(paw_opt/3)*ndat * dp)

        ! reminder : a call to cudaDeviceSynchronize is needed so that svectout can be re-used safely on host
        ! don't do it here, only when really needed

      end if ! (paw_opt == 3 .or. paw_opt == 4)

      if (paw_opt == 0 .or. paw_opt == 1 .or. paw_opt == 4) then

        ! Get vectout from vnl_projections
        if (cplex_fac == 2) then

          call abi_gpu_xgemm(cplex, 'N', 'N', &
            &            npwout, ndat*nspinor, nprojs, &                                                ! M,N,K
            &            cone, &                                                                        ! alpha
            &            gemm_nonlop_kpt_gpu(gemm_nonlop_ikpt_this_proc_being_treated)%projs, npwout, & ! A, LDA
            &            gemm_nonlop_kokkos%vnl_projections_gpu, nprojs, &                              ! B, LDB
            &            czero, &                                                                       ! beta
            &            C_LOC(vectout), npwout)                                        ! C, LDC

        else

          if (.not. gpu_allocated(temp_realvec_gpu)) then
            ABI_ERROR("Please provide memory allocation for temp_realvec_gpu array")
          end if

          call abi_gpu_xgemm(cplex, 'N', 'N', &
            &            npwout, ndat*nspinor, nprojs, &                                               ! M,N,K
            &            cone, &                                                                       ! alpha
            &            gemm_nonlop_kpt_gpu(gemm_nonlop_ikpt_this_proc_being_treated)%projs_r, npwout, &  ! A, LDA
            &            gemm_nonlop_kokkos%vnl_projections_gpu, nprojs, &                             ! B, LDB
            &            czero, &                                                                      ! beta
            &            temp_realvec_gpu, npwout)                                                     ! C, LDC
          ! vectout(1,1:npwout*nspinor*ndat) = temp_realvec(1:npwout*nspinor*ndat)
          call insert_real_part(C_LOC(vectout), temp_realvec_gpu, npwout*nspinor*ndat)

          call abi_gpu_xgemm(cplex, 'N', 'N', &
            &            npwout, ndat*nspinor, nprojs, &                                               ! M,N,K
            &            cone, &                                                                       ! alpha
            &            gemm_nonlop_kpt_gpu(gemm_nonlop_ikpt_this_proc_being_treated)%projs_i, npwout, & ! A, LDA
            &            gemm_nonlop_kokkos%vnl_projections_gpu, nprojs, &                             ! B, LDB
            &            czero, &                                                                      ! beta
            &            temp_realvec_gpu, npwout)                                                     ! C, LDC
          ! vectout(2,1:npwout*nspinor*ndat) = temp_realvec(1:npwout*nspinor*ndat)
          call insert_imag_part(C_LOC(vectout), temp_realvec_gpu, npwout*nspinor*ndat)

        end if  ! cplex_fac == 2

        ! copy back results on host
        !call copy_from_gpu(C_LOC(vectout(1,1)), gemm_nonlop_kokkos%vectout_gpu, INT(2, c_size_t)*npwout*nspinor*ndat * dp)

      end if ! (paw_opt == 0 .or. paw_opt == 1 .or. paw_opt == 4)

    end if ! choice > 0

  if (dimekbq>1) then
    ABI_FREE_CUDA(vnl_projections2_gpu)
  end if

  !ABI_FREE_CUDA(      vectin_gpu)
  !ABI_FREE_CUDA(     vectout_gpu)
  !ABI_FREE_CUDA(    svectout_gpu)

  ! device sync before reusing data computed on device
  ! this may not be the best location to have this sync
  call gpu_device_synchronize()

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
 !!***


#else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! stubs for compiling with CUDA disabled.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine init_gemm_nonlop_gpu(nkpt)
  integer,intent(in) :: nkpt

  ABI_UNUSED((/nkpt/))
  ABI_BUG("Unhandled configuration for CUDA immplementation")
 end subroutine init_gemm_nonlop_gpu

 subroutine destroy_gemm_nonlop_gpu(nkpt)
  integer,intent(in) :: nkpt

  ABI_UNUSED((/nkpt/))
  ABI_BUG("Unhandled configuration for CUDA immplementation")
 end subroutine destroy_gemm_nonlop_gpu

 subroutine make_gemm_nonlop_gpu(ikpt,npw,lmnmax,ntypat,indlmn,nattyp,istwf_k,ucvol,ffnl_k, &
&                            ph3d_k,kpt_k,kg_k,kpg_k, &
&                            compute_grad_strain,compute_grad_atom) ! Optional parameters

  integer, intent(in) :: ikpt
  integer, intent(in) :: npw, lmnmax,ntypat
  integer, intent(in) :: indlmn(:,:,:), kg_k(:,:)
  integer, intent(in) :: nattyp(ntypat)
  integer, intent(in) :: istwf_k
  logical, intent(in), optional :: compute_grad_strain,compute_grad_atom
  real(dp), intent(in) :: ucvol
  real(dp), intent(in) :: ffnl_k(:,:,:,:)
  real(dp), intent(in) :: ph3d_k(:,:,:)
  real(dp), intent(in) :: kpt_k(:)
  real(dp), intent(in), target :: kpg_k(:,:)

  ABI_UNUSED((/ikpt,npw,lmnmax,ntypat,indlmn,kg_k,nattyp,istwf_k/))
  ABI_UNUSED((/ucvol,ffnl_k,ph3d_k,kpt_k,kpg_k/))
  ABI_UNUSED((/compute_grad_strain,compute_grad_atom/))
  ABI_BUG("Unhandled configuration for CUDA immplementation")

 end subroutine make_gemm_nonlop_gpu

 subroutine gemm_nonlop_gpu(atindx1,choice,cpopt,cprjin,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,&
&                 enl,enlout,ffnlin,ffnlout,gmet,gprimd,idir,indlmn,istwf_k,&
&                 kgin,kgout,kpgin,kpgout,kptin,kptout,lambda,lmnmax,matblk,mgfft,&
&                 mpi_enreg,mpsang,mpssoang,natom,nattyp,ndat,ngfft,nkpgin,nkpgout,nloalg,&
&                 nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO,paw_opt,phkxredin,&
&                 phkxredout,ph1d,ph3din,ph3dout,signs,sij,svectout,&
&                 tim_nonlop,ucvol,useylm,vectin,vectout,&
&                 vectproj,use_gpu_cuda)

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
  real(dp),intent(inout),target :: vectin(2,npwin*nspinor*ndat)
  real(dp),intent(inout) :: enlout(nnlout*ndat)
  real(dp),intent(out),target :: svectout(2,npwout*nspinor*(paw_opt/3)*ndat)
  real(dp),intent(inout),target :: vectout(2,npwout*nspinor*ndat) !vz_i
  real(dp),intent(inout),optional,target :: vectproj(:,:,:)
  type(pawcprj_type),intent(inout) :: cprjin(natom,nspinor*((cpopt+5)/5)*ndat)

  ABI_UNUSED((/choice,cpopt,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,idir/))
  ABI_UNUSED((/istwf_k,lmnmax,matblk,mgfft,mpsang,mpssoang,natom,ndat,nkpgin/))
  ABI_UNUSED((/nkpgout,nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO/))
  ABI_UNUSED((/paw_opt,signs,tim_nonlop,useylm,use_gpu_cuda/))
  ABI_UNUSED((/atindx1,indlmn,kgin,kgout,nattyp,ngfft,nloalg/))
  ABI_UNUSED((/enl,ffnlin,ffnlout,gmet,gprimd,kpgin,kpgout,kptin,kptout,phkxredin,phkxredout/))
  ABI_UNUSED((/ucvol,lambda,sij,ph1d(1,1),ph3din,ph3dout,vectin,enlout,svectout,vectout,vectproj/))
  ABI_UNUSED_A(cprjin)
  ABI_UNUSED_A(mpi_enreg)
  ABI_BUG("Unhandled configuration for CUDA immplementation")

 end subroutine gemm_nonlop_gpu

#endif

end module m_gemm_nonlop_gpu
!!***

