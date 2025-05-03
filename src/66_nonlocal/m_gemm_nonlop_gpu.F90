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
!! Copyright (C) 2014-2025 ABINIT group (AL,MS)
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
 use m_opernlc_ylm_allwf, only : opernlc_ylm_allwf
 use m_pawcprj, only : pawcprj_type
 use m_gemm_nonlop_projectors
 use m_hamiltonian, only : KPRIME_H_K, K_H_KPRIME, K_H_K, KPRIME_H_KPRIME

#if defined(HAVE_GPU_CUDA)
 use m_gpu_toolbox
 use m_alloc_hamilt_gpu, only : gemm_nonlop_gpu_data
#endif

#ifdef HAVE_KOKKOS
 use m_manage_kokkos, only : opernlc_ylm_allwf_kokkos
#endif

#ifdef HAVE_FC_ISO_C_BINDING
 use, intrinsic :: iso_c_binding, only : c_ptr, c_int32_t, c_int64_t, c_float, c_double, c_size_t, c_loc
#endif

 implicit none

 private

 public :: gemm_nonlop_gpu

 !!***

 !!
 !! GPU interface
 !!
#if defined(HAVE_FC_ISO_C_BINDING) && defined(HAVE_GPU_CUDA)

 interface

   !> extract real part of a complex vector
   !> data_in and data_out must be pointers in device memory
   subroutine extract_real_part(data_out, data_in, size) bind(c, name='cuda_extract_real_part')
     use, intrinsic :: iso_c_binding
     type(c_ptr),             value :: data_out
     type(c_ptr),             value :: data_in
     integer(kind=c_int32_t), value :: size
   end subroutine extract_real_part

   !> extract real part of a complex vector
   !> data_in and data_out must be pointers in device memory
   subroutine extract_imag_part(data_out, data_in, size) bind(c, name='cuda_extract_imag_part')
     use, intrinsic :: iso_c_binding
     type(c_ptr),             value :: data_out
     type(c_ptr),             value :: data_in
     integer(kind=c_int32_t), value :: size
   end subroutine extract_imag_part

   !> insert real part of a complex vector
   !> data_in and data_out must be pointers in device memory
   subroutine insert_real_part(data_out, data_in, size) bind(c, name='cuda_insert_real_part')
     use, intrinsic :: iso_c_binding
     type(c_ptr),             value :: data_out
     type(c_ptr),             value :: data_in
     integer(kind=c_int32_t), value :: size
   end subroutine insert_real_part

   !> insert real part of a complex vector
   !> data_in and data_out must be pointers in device memory
   subroutine insert_imag_part(data_out, data_in, size) bind(c, name='cuda_insert_imag_part')
     use, intrinsic :: iso_c_binding
     type(c_ptr),             value :: data_out
     type(c_ptr),             value :: data_in
     integer(kind=c_int32_t), value :: size
   end subroutine insert_imag_part

   !> data_in and data_out must be pointers in device memory
   subroutine fix_realvec(data, npw_in, ndat_nspinor, option) bind(c, name='cuda_fix_realvec')
     use, intrinsic :: iso_c_binding
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

!!****f* m_gemm_nonlop_gpu/gemm_nonlop_gpu
!! NAME
!! gemm_nonlop_gpu
!!
!! FUNCTION
!! Replacement of nonlop.
!!
!! INPUTS
!! [gpu_option] = GPU implementation to use, i.e. cuda, openMP, ... (0=not using GPU)
!!
!! PARENTS
!!      m_nonlop
!!
!! CHILDREN
!!      abi_zgemm_2r,dgemm,opernlc_ylm,xmpi_sum
!!
!! SOURCE
 subroutine gemm_nonlop_gpu(atindx1,choice,cpopt,cprjin,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,&
&                 enl,ffnlin,ffnlout,indlmn,istwf_k,&
&                 lambda,lmnmax,matblk,&
&                 mpi_enreg,natom,nattyp,ndat,nkpgin,nkpgout,&
&                 nnlout,npwin,npwout,nspinor,nspinortot,ntypat,paw_opt,&
&                 ph3din,ph3dout,sij,svectout,&
&                 ucvol,useylm,vectin,vectout,select_k,&
&                 gpu_option,vectproj)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: choice,cpopt,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout
  integer,intent(in) :: istwf_k,lmnmax,matblk,natom,ndat,nkpgin
  integer,intent(in) :: nkpgout,nnlout,npwin,npwout,nspinor,nspinortot,ntypat
  integer,intent(in) :: paw_opt,useylm,select_k
  integer,intent(in) :: gpu_option
  real(dp),intent(in) :: ucvol
  real(dp),target, intent(in)     :: lambda(ndat)
  type(MPI_type),intent(in) :: mpi_enreg
  !arrays
  integer,intent(in),target :: atindx1(natom),indlmn(6,lmnmax,ntypat)
  integer,intent(in),target :: nattyp(ntypat)
  real(dp),intent(in),target :: enl(dimenl1,dimenl2,nspinortot**2,dimekbq)
  real(dp),intent(in),target :: ffnlin(npwin,dimffnlin,lmnmax,ntypat)
  real(dp),intent(in),target :: ffnlout(npwout,dimffnlout,lmnmax,ntypat)
  real(dp),intent(in),target :: sij(dimenl1,ntypat*((paw_opt+1)/3))
  real(dp),intent(inout),target :: ph3din(2,npwin,matblk),ph3dout(2,npwout,matblk)
  real(dp),intent(inout),target :: vectin(2,npwin*nspinor*ndat)
  real(dp),intent(out),target :: svectout(:,:)
  real(dp),intent(inout),target :: vectout(:,:)
  real(dp),intent(inout),optional, ABI_CONTIGUOUS target :: vectproj(:,:,:)
  type(pawcprj_type),intent(inout) :: cprjin(natom,nspinor*((cpopt+5)/5)*ndat)

  ! locals
  integer :: idat, nprojs, shift, iatom, nlmn, ierr, ibeg, iend, ikin, ikout
  integer :: cplex, cplex_enl, cplex_fac
  integer :: iatm, ndgxdt, ndgxdtfac, nd2gxdt, nd2gxdtfac, optder, itypat, ilmn
  integer :: cplex_dgxdt(9), cplex_d2gxdt(18)
  logical :: local_vectproj
  real(dp) :: dgxdt_dum_in(1,1,1,1,1), dgxdt_dum_out(1,1,1,1,1),dgxdt_dum_out2(1,1,1,1,1)
  real(dp) :: d2gxdt_dum_in(1,1,1,1,1), d2gxdt_dum_out(1,1,1,1,1),d2gxdt_dum_out2(1,1,1,1,1)
  real(dp), allocatable, target :: sij_typ(:)
  integer :: npw_max
  integer :: nattyp_max


  !type(c_ptr)                      :: projections_gpu,        s_projections_gpu,        vnl_projections_gpu
  real(dp), ABI_CONTIGUOUS pointer :: projections_cpu(:,:,:)
  real(dp),    allocatable, target :: s_projections_cpu(:,:,:), vnl_projections_cpu(:,:,:)

  ! used inside opernlc_ylm_allwf_kokkos_cpp when iphase > 1
  type(c_ptr)                      :: vnl_projections2_gpu

  type(c_ptr)                      :: temp_realvec_gpu

  ! GPU waveform data are allocated in m_alloc_hamilt_gpu
  !type(c_ptr)                      :: vectin_gpu, vectout_gpu, svectout_gpu
  ! Pointers to either CUDA allocated or managed memory data
  type(c_ptr)                      :: vectin_ptr, vectout_ptr, svectout_ptr
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

  ! This function should only be called within CUDA legacy or Kokkos code paths
  if(gpu_option/=ABI_GPU_LEGACY .and. gpu_option/=ABI_GPU_KOKKOS) then
    ABI_BUG("Unhandled GPU value !")
  end if

  npw_max = MAX(npwin, npwout)

  ikin=1; ikout=1;
  select case (select_k)
  case (K_H_K)
    ikin=1; ikout=1;
  case (K_H_KPRIME)
    ikin=2; ikout=1;
  case (KPRIME_H_K)
    ikin=1; ikout=2;
  case (KPRIME_H_KPRIME)
    ikin=2; ikout=2;
  end select
  cplex=2;if (istwf_k>1) cplex=1
  cplex_enl=1;if (paw_opt>0) cplex_enl=2*dimenl1/(lmnmax*(lmnmax+1)) ! is enl complex?
  cplex_fac=max(cplex,dimekbq)
  if ((nspinortot==2.or.cplex_enl==2).and.paw_opt>0.and.choice/=7) cplex_fac=2 ! is vnl_projections complex?
  ndgxdt=0; nd2gxdt=0

  ! Compute projectors if need be
  ! FIXME Would be nice to rely on opernl(ab)_gemm instead but some refacto is needed
  nprojs=0
  do itypat=1,ntypat
    nprojs = nprojs + count(indlmn(3,:,itypat)>0)*nattyp(itypat)
  end do

  call refresh_projectors(npwin,istwf_k,nprojs,ndgxdt,nd2gxdt,(ikin==2),&
  &                       gpu_option)
  if(gemm_nonlop_kpt(ikin)%ikpt/=gemm_nonlop_ikpt_this_proc_being_treated) then
    call prep_projectors(npwin,lmnmax,ntypat,indlmn,nattyp,istwf_k,&
    &                    ucvol,ffnlin,ph3din,dimffnlin,matblk,&
    &                    nprojs,(ikin==2),gpu_option,0)
    gemm_nonlop_kpt(ikin)%ikpt=gemm_nonlop_ikpt_this_proc_being_treated
  end if
  if(ikin/=ikout .and. choice>0) then
    call refresh_projectors(npwout,istwf_k,nprojs,ndgxdt,nd2gxdt,(ikout==2),&
    &                       gpu_option)
    if(gemm_nonlop_kpt(ikout)%ikpt/=gemm_nonlop_ikpt_this_proc_being_treated) then
      call prep_projectors(npwout,lmnmax,ntypat,indlmn,nattyp,istwf_k,&
      &                    ucvol,ffnlout,ph3dout,dimffnlout,matblk,&
      &                    nprojs,(ikout==2),gpu_option,0)
      gemm_nonlop_kpt(ikout)%ikpt=gemm_nonlop_ikpt_this_proc_being_treated
    end if
  end if

  ! These will store the non-local factors for vectin, svectout and vectout respectively
  call gpu_memset(gemm_nonlop_gpu_data%    projections_gpu, izero, INT(cplex,     c_size_t) * nprojs * nspinor*ndat * dp)
  call gpu_memset(gemm_nonlop_gpu_data%  s_projections_gpu, izero, INT(cplex,     c_size_t) * nprojs * nspinor*ndat * dp)
  call gpu_memset(gemm_nonlop_gpu_data%vnl_projections_gpu, izero, INT(cplex_fac, c_size_t) * nprojs * nspinor*ndat * dp)

  if (dimekbq>1) then
    ABI_CHECK(cplex_fac==2,"BUG: invalid cplex_fac==1 when dimekbq=2!")
    ABI_MALLOC_CUDA(vnl_projections2_gpu,        INT(cplex_fac, c_size_t) * nprojs * nspinor*ndat * dp)
    call gpu_memset(vnl_projections2_gpu, izero, INT(cplex_fac, c_size_t) * nprojs * nspinor*ndat * dp)
  end if

  vectin_size = 2*npwin*nspinor*ndat*dp

  if(gpu_option==ABI_GPU_LEGACY) then
    vectin_ptr  = gemm_nonlop_gpu_data%vectin_gpu
    vectout_ptr = gemm_nonlop_gpu_data%vectout_gpu
    svectout_ptr= gemm_nonlop_gpu_data%svectout_gpu

    call copy_on_gpu(vectin, gemm_nonlop_gpu_data%vectin_gpu, vectin_size)
    if (choice == 7) then
      call copy_on_gpu(svectout, gemm_nonlop_gpu_data%svectout_gpu, INT(2, c_size_t) * npwout*nspinor*ndat * dp)
    end if

  else if(gpu_option==ABI_GPU_KOKKOS) then
    vectin_ptr  =C_LOC(vectin)
    vectout_ptr =C_LOC(vectout)
    svectout_ptr=C_LOC(svectout)

    call gpu_data_prefetch_async(vectin_ptr, vectin_size)
    if (choice == 7) then
      call gpu_data_prefetch_async(C_LOC(svectout), vectin_size)
    end if

  end if

  !! gpu alloc and init : enl_gpu
  enl_size_bytes = dimenl1 * dimenl2 * nspinortot**2 * dimekbq * dp
  ABI_MALLOC_CUDA( enl_gpu, enl_size_bytes )
  call copy_on_gpu(enl, enl_gpu, enl_size_bytes )

  !! gpu alloc and init atindx1_gpu
  sizeof_int = 4
  atindx1_size_bytes = natom * sizeof_int
  ABI_MALLOC_CUDA( atindx1_gpu,  atindx1_size_bytes )
  call copy_on_gpu(atindx1, atindx1_gpu, atindx1_size_bytes )

  !! gpu alloc and init indlmn_gpu
  indlmn_size_bytes = 6*lmnmax*ntypat * sizeof_int
  ABI_MALLOC_CUDA( indlmn_gpu,  indlmn_size_bytes )
  call copy_on_gpu(indlmn, indlmn_gpu, indlmn_size_bytes )

  !! gpu alloc and init lambda_gpu
  lambda_size_bytes =  ndat * dp
  ABI_MALLOC_CUDA( lambda_gpu, lambda_size_bytes )
  call copy_on_gpu(lambda, lambda_gpu, lambda_size_bytes )

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

  ! If vectproj is provided, use it for further calculations, use allocated array otherwise
  local_vectproj=.false.
  if(PRESENT(vectproj)) then
    if(size(vectproj)>1) local_vectproj=.true.
  end if
  if (local_vectproj) then
    projections_cpu => vectproj
  else
    ABI_MALLOC(    projections_cpu,(cplex,     nprojs,nspinor*ndat))
    projections_cpu = zero
  end if
  ABI_MALLOC(  s_projections_cpu,(cplex,     nprojs,nspinor*ndat)) ! TODO - TO BE REMOVED ONCE CUDA-IZATION IS OK
  ABI_MALLOC(vnl_projections_cpu,(cplex_fac, nprojs,nspinor*ndat)) ! TODO - TO BE REMOVED ONCE CUDA-IZATION IS OK
  s_projections_cpu = zero
  vnl_projections_cpu = zero

  if(cpopt >= 2) then

    if(.not. local_vectproj) then
      !This use-case is extremely painful for GEMM GPU nonlop performance.
      !vectproj paramter should always be provided to avoid it

      ! retrieve from cprjin
      do idat=1, ndat*nspinor
        shift = 0
        do iatom = 1, natom
          nlmn = cprjin(iatom, idat)%nlmn
          projections_cpu(1:cplex, shift+1:shift+nlmn, idat) = cprjin(iatom, idat)%cp(1:cplex, 1:nlmn)
          shift = shift + nlmn
        end do
      end do
    end if

    ! copy from HOST projections_cpu to GPU projections_gpu
    call copy_on_gpu(projections_cpu, gemm_nonlop_gpu_data%projections_gpu,&
        INT(cplex, c_size_t) * nprojs * nspinor*ndat * dp)

  else ! cpopt < 2

    ! opernla
    if(cplex == 2) then

      ! projections_gpu = projs * vectin_gpu
      call abi_gpu_xgemm(cplex, 'C', 'N', &
        &            nprojs, ndat*nspinor, npwin, &                                                ! M,N,K
        &            cone, &                                                                       ! alpha
        &            gemm_nonlop_kpt_gpu(ikin)%projs, npwin, & ! A, LDA
        &            vectin_ptr, npwin, &                                                          ! B, LDB
        &            czero, &                                                                      ! beta
        &            gemm_nonlop_gpu_data%projections_gpu, nprojs)                                 ! C, LDC

    else

      if (.not. gpu_allocated(temp_realvec_gpu)) then
        ABI_ERROR("Please provide memory allocation for temp_realvec_gpu array")
      end if


      ! only compute real part of projections = P^* psi => projections_r = P_r^T psi_r + P_i^T psi_i
      !temp_realvec(1:npwin*nspinor*ndat) = vectin(1,1:npwin*nspinor*ndat)
      call extract_real_part(temp_realvec_gpu, vectin_ptr, npwin*nspinor*ndat)

      if(istwf_k == 2 .and. mpi_enreg%me_g0_fft == 1) then
        ! do idat=1, ndat*nspinor
        !   temp_realvec(1+(idat-1)*npwin) = temp_realvec(1+(idat-1)*npwin)/2
        ! end do
        call fix_realvec(temp_realvec_gpu, npwin, ndat*nspinor, fix_realvec_divide_by_2)
      end if

      call abi_gpu_xgemm(cplex, 'T', 'N', &
        &            nprojs, ndat*nspinor, npwin, &                                                  ! M,N,K
        &            cone, &                                                                         ! alpha
        &            gemm_nonlop_kpt_gpu(ikin)%projs_r, npwin, & ! A, LDA
        &            temp_realvec_gpu, npwin, &                                                      ! B, LDB
        &            czero, &                                                                        ! beta
        &            gemm_nonlop_gpu_data%projections_gpu, nprojs)                                   ! C, LDC

      !temp_realvec(1:npwin*nspinor*ndat) = vectin(2,1:npwin*nspinor*ndat)
      call extract_imag_part(temp_realvec_gpu, vectin_ptr, npwin*nspinor*ndat)

      if(istwf_k == 2 .and. mpi_enreg%me_g0_fft == 1) then
        ! do idat=1, ndat*nspinor
        !   temp_realvec(1+(idat-1)*npwin) = zero
        ! end do
        call fix_realvec(temp_realvec_gpu, npwin, ndat*nspinor, fix_realvec_zero_out)
      end if
      call abi_gpu_xgemm(cplex, 'T', 'N', &
        &            nprojs, ndat*nspinor, npwin, &
        &            cone, &
        &            gemm_nonlop_kpt_gpu(ikin)%projs_i, npwin, &
        &            temp_realvec_gpu, npwin, &
        &            cone, &
        &            gemm_nonlop_gpu_data%projections_gpu, nprojs)

      !gemm_nonlop_gpu_data%projections_gpu = 2 * gemm_nonlop_gpu_data%projections_gpu
      call gpu_xscal(cplex, nprojs*nspinor*ndat, ctwo, gemm_nonlop_gpu_data%projections_gpu, 1)

    end if ! cplex == 2

!    call xmpi_sum(projections,mpi_enreg%comm_fft,ierr)

    if(cpopt >= 0) then
      ! copy from GPU projections_gpu to HOST projections_cpu
      call copy_from_gpu(projections_cpu, gemm_nonlop_gpu_data%projections_gpu,&
          INT(cplex, c_size_t) * nprojs * nspinor*ndat * dp)

      if(.not. local_vectproj) then
        !This use-case is extremely painful for GEMM GPU nonlop performance.
        !vectproj paramter should always be provided to avoid it

        ! store in cprjin
        do idat=1, ndat*nspinor
          shift = 0
          do iatom = 1, natom
            nlmn = cprjin(iatom, idat)%nlmn
            cprjin(iatom, idat)%cp(1:cplex, 1:nlmn) = projections_cpu(1:cplex, shift+1:shift+nlmn, idat)
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

          call copy_on_gpu(sij_typ, sij_typ_gpu, sij_typ_size_bytes )

        end if ! paw_opt>=2

        ibeg = shift+1
        iend = shift+nattyp(itypat)*nlmn

        !Parallelization over spinors treatment
        shift_spinor = 0
        if (mpi_enreg%paral_spinor==1) then
          shift_spinor = mpi_enreg%me_spinor
        end if

        ! Use the Kokkos implementation of opernlc if available
        if (gpu_option == ABI_GPU_KOKKOS) then
#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS)
          call opernlc_ylm_allwf_kokkos(cplex, cplex_enl, cplex_fac, &
            &                           dimenl1, dimenl2, dimekbq, &
            &                           iatm, itypat, ntypat, nprojs, &
            &                           natom, nattyp(itypat), nspinor, &
            &                           nspinortot, paw_opt, &
            &                           nlmn, lmnmax, &
            &                           enl_gpu, &
            &                           gemm_nonlop_gpu_data%projections_gpu, &
            &                           gemm_nonlop_gpu_data%vnl_projections_gpu, &
            &                           vnl_projections2_gpu, &
            &                           gemm_nonlop_gpu_data%s_projections_gpu, &
            &                           shift_spinor, ndat, &
            &                           atindx1_gpu, &
            &                           indlmn_gpu, &
            &                           lambda_gpu, &
            &                           sij_typ_gpu, &
            &                           shift, nattyp_max)
#endif
        ! Fall back on CPU implementation of opernlc
        else

          call copy_from_gpu(  projections_cpu, gemm_nonlop_gpu_data%projections_gpu, &
            &              INT(cplex,     c_size_t) * nprojs * nspinor*ndat * dp)

          call opernlc_ylm_allwf(atindx1, cplex, cplex_dgxdt, cplex_d2gxdt,&
            &                  cplex_enl, cplex_fac, &
            &                  dgxdt_dum_in, dgxdt_dum_out, dgxdt_dum_out2, &
            &                  d2gxdt_dum_in, d2gxdt_dum_out, d2gxdt_dum_out2, &
            &                  dimenl1, dimenl2, dimekbq, enl, &
            &                  projections_cpu, &
            &                  vnl_projections_cpu, &
            &                  s_projections_cpu, &
            &                  iatm, indlmn(:,:,itypat), itypat, lambda, mpi_enreg, natom, &
            &                  ndgxdt,ndgxdtfac,nd2gxdt,nd2gxdtfac,&
            &                  nattyp(itypat), nlmn, nspinor, nspinortot, optder, paw_opt, sij_typ, &
            &                  ndat, ibeg-1, iend, nprojs, 1, ABI_GPU_DISABLED)

          call copy_on_gpu(  s_projections_cpu, gemm_nonlop_gpu_data%  s_projections_gpu, &
            &              INT(cplex,     c_size_t) * nprojs * nspinor*ndat * dp)
          call copy_on_gpu(vnl_projections_cpu, gemm_nonlop_gpu_data%vnl_projections_gpu, &
            &              INT(cplex_fac, c_size_t) * nprojs * nspinor*ndat * dp)
        end if

        shift = shift + nattyp(itypat)*nlmn
        iatm  = iatm  + nattyp(itypat)
      end do ! itypat

      ABI_FREE(sij_typ)
      ABI_FREE_CUDA( sij_typ_gpu )

    else ! choice == 7

      ! TO BE REMOVED - DEBUG ONLY
      call copy_gpu_to_gpu(gemm_nonlop_gpu_data%s_projections_gpu, &
           &               gemm_nonlop_gpu_data%projections_gpu, &
           &               INT(cplex, c_size_t) * nprojs * nspinor * ndat * dp)

    end if ! choice

    ! opernlb
    if (paw_opt == 3 .or. paw_opt == 4) then

      ! Get svectout from s_projections
      if(cplex == 2) then

        call abi_gpu_xgemm(cplex, 'N', 'N', &
          &            npwout, ndat*nspinor, nprojs, &                                                ! M,N,K
          &            cone, &                                                                        ! alpha
          &            gemm_nonlop_kpt_gpu(ikout)%projs, npwout, & ! A, LDA
          &            gemm_nonlop_gpu_data%s_projections_gpu, nprojs, &                              ! B, LDB
          &            czero, &                                                                       ! beta
          &            svectout_ptr, npwout)                                                          ! C, LDC

      else

        if (.not. gpu_allocated(temp_realvec_gpu)) then
          ABI_ERROR("Please provide memory allocation for temp_realvec_gpu array")
        end if

        call abi_gpu_xgemm(cplex, 'N', 'N', &
          &            npwout, ndat*nspinor, nprojs, &                                                ! M,N,K
          &            cone, &                                                                        ! alpha
          &            gemm_nonlop_kpt_gpu(ikout)%projs_r, npwout, & ! A, LDA
          &            gemm_nonlop_gpu_data%s_projections_gpu, nprojs, &                              ! B, LDB
          &            czero, &                                                                       ! beta
          &            temp_realvec_gpu, npwout)                                                      ! C, LDC
        !svectout(1,1:npwout*nspinor*ndat) = temp_realvec_gpu(1:npwout*nspinor*ndat)
        call insert_real_part(svectout_ptr, temp_realvec_gpu, npwout*nspinor*ndat)

        call abi_gpu_xgemm(cplex, 'N', 'N', &
          &            npwout, ndat*nspinor, nprojs, &                                                ! M,N,K
          &            cone, &                                                                        ! alpha
          &            gemm_nonlop_kpt_gpu(ikout)%projs_i, npwout,& ! A, LDA
          &            gemm_nonlop_gpu_data%s_projections_gpu, nprojs, &                              ! B, LDB
          &            czero, &                                                                       ! beta
          &            temp_realvec_gpu, npwout)                                                      ! C, LDC
        !svectout(2,1:npwout*nspinor*ndat) = temp_realvec(1:npwout*nspinor*ndat)
        call insert_imag_part(svectout_ptr, temp_realvec_gpu, npwout*nspinor*ndat)

      end if ! cplex == 2

      if (choice /= 7) then

        ! compute : svectout_gpu = svectout_gpu + vectin_gpu
        ! this a axpy operation with x => vectin_gpu, y => svectout_gpu and alpha=1
        ! please remember that svectout_gpu and vectin_gpu have same size when paw_opt >= 3 and paw_opt<6
        ! this is the case here
        call abi_gpu_xaxpy(1, &                              ! real
          &            2*npwin*nspinor*ndat, &               ! size
          &            cone, &                               ! alpha
          &            vectin_ptr, 1, &                      ! X, incrx
          &            svectout_ptr, 1)                      ! Y, incry

      endif

      if(gpu_option==ABI_GPU_LEGACY) then
        ! copy back results on host
        call copy_from_gpu(svectout, svectout_ptr, INT(2, c_size_t)*npwout*nspinor*(paw_opt/3)*ndat * dp)
      end if

      ! reminder : a call to cudaDeviceSynchronize is needed so that svectout can be re-used safely on host
      ! don't do it here, only when really needed

    end if ! (paw_opt == 3 .or. paw_opt == 4)

    if (paw_opt == 0 .or. paw_opt == 1 .or. paw_opt == 4) then

      ! Get vectout from vnl_projections
      if (cplex_fac == 2) then

        call abi_gpu_xgemm(cplex, 'N', 'N', &
          &            npwout, ndat*nspinor, nprojs, &                                                ! M,N,K
          &            cone, &                                                                        ! alpha
          &            gemm_nonlop_kpt_gpu(ikout)%projs, npwout, & ! A, LDA
          &            gemm_nonlop_gpu_data%vnl_projections_gpu, nprojs, &                            ! B, LDB
          &            czero, &                                                                       ! beta
          &            vectout_ptr, npwout)                                                           ! C, LDC

      else

        if (.not. gpu_allocated(temp_realvec_gpu)) then
          ABI_ERROR("Please provide memory allocation for temp_realvec_gpu array")
        end if

        call abi_gpu_xgemm(cplex, 'N', 'N', &
          &            npwout, ndat*nspinor, nprojs, &                                               ! M,N,K
          &            cone, &                                                                       ! alpha
          &            gemm_nonlop_kpt_gpu(ikout)%projs_r, npwout, &  ! A, LDA
          &            gemm_nonlop_gpu_data%vnl_projections_gpu, nprojs, &                           ! B, LDB
          &            czero, &                                                                      ! beta
          &            temp_realvec_gpu, npwout)                                                     ! C, LDC
        ! vectout(1,1:npwout*nspinor*ndat) = temp_realvec(1:npwout*nspinor*ndat)
        call insert_real_part(vectout_ptr, temp_realvec_gpu, npwout*nspinor*ndat)

        call abi_gpu_xgemm(cplex, 'N', 'N', &
          &            npwout, ndat*nspinor, nprojs, &                                               ! M,N,K
          &            cone, &                                                                       ! alpha
          &            gemm_nonlop_kpt_gpu(ikout)%projs_i, npwout, & ! A, LDA
          &            gemm_nonlop_gpu_data%vnl_projections_gpu, nprojs, &                           ! B, LDB
          &            czero, &                                                                      ! beta
          &            temp_realvec_gpu, npwout)                                                     ! C, LDC
        ! vectout(2,1:npwout*nspinor*ndat) = temp_realvec(1:npwout*nspinor*ndat)
        call insert_imag_part(vectout_ptr, temp_realvec_gpu, npwout*nspinor*ndat)

      end if  ! cplex_fac == 2

      ! copy back results on host
      if(gpu_option==ABI_GPU_LEGACY) then
        call copy_from_gpu(vectout, vectout_ptr, INT(2, c_size_t)*npwout*nspinor*ndat * dp)
      end if

    end if ! (paw_opt == 0 .or. paw_opt == 1 .or. paw_opt == 4)

  end if ! choice > 0

  if (dimekbq>1) then
    ABI_FREE_CUDA(vnl_projections2_gpu)
  end if

  ! device sync before reusing data computed on device
  ! this may not be the best location to have this sync
  call gpu_device_synchronize()

  if (cplex /= 2) then
    if ( (cpopt < 2) .or. &
      &  (paw_opt == 3 .or. paw_opt == 4) .or. &
      &  (paw_opt == 0 .or. paw_opt == 1 .or. paw_opt == 4)) then
      ABI_FREE_CUDA(temp_realvec_gpu)
    end if
  end if

  ABI_FREE_CUDA( enl_gpu )
  ABI_FREE_CUDA( atindx1_gpu )
  ABI_FREE_CUDA( indlmn_gpu )
  ABI_FREE_CUDA( lambda_gpu )

  ! if projections_cpu was allocated, then free it here
  if(.not. local_vectproj) then
    ABI_FREE(projections_cpu)
  end if
  ABI_FREE(  s_projections_cpu) ! TO BE REMOVED
  ABI_FREE(vnl_projections_cpu) ! TO BE REMOVED

 end subroutine gemm_nonlop_gpu
 !!***


#else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! stubs for compiling with CUDA disabled.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine gemm_nonlop_gpu(atindx1,choice,cpopt,cprjin,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,&
&                 enl,enlout,ffnlin,ffnlout,gmet,gprimd,idir,indlmn,istwf_k,&
&                 kgin,kgout,kpgin,kpgout,kptin,kptout,lambda,lmnmax,matblk,mgfft,&
&                 mpi_enreg,mpsang,mpssoang,natom,nattyp,ndat,ngfft,nkpgin,nkpgout,nloalg,&
&                 nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO,paw_opt,phkxredin,&
&                 phkxredout,ph1d,ph3din,ph3dout,signs,sij,svectout,&
&                 tim_nonlop,ucvol,useylm,vectin,vectout,select_k,&
&                 gpu_option,vectproj)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: choice,cpopt,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,idir
  integer,intent(in) :: istwf_k,lmnmax,matblk,mgfft,mpsang,mpssoang,natom,ndat,nkpgin
  integer,intent(in) :: nkpgout,nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO
  integer,intent(in) :: paw_opt,signs,tim_nonlop,useylm,select_k
  integer,intent(in) :: gpu_option
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
  ABI_UNUSED((/istwf_k,lmnmax,matblk,mgfft,mpsang,mpssoang,natom,ndat,nkpgin,select_k/))
  ABI_UNUSED((/nkpgout,nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO/))
  ABI_UNUSED((/paw_opt,signs,tim_nonlop,useylm,gpu_option/))
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

