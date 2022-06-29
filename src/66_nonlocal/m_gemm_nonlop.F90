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
 use m_abi_linalg

 use defs_abitypes, only : MPI_type
 use m_opernlc_ylm, only : opernlc_ylm
 use m_pawcprj, only : pawcprj_type
 use m_geometry, only : strconv
 use m_kg, only : mkkpg

 implicit none

 private

 ! Use these routines in order: first call init, then call make_gemm_nonlop for each k point,
 ! then call gemm_nonlop to do the actual computation, and call destroy when done. See gstate and vtorho.
 public :: init_gemm_nonlop
 public :: destroy_gemm_nonlop
 public :: make_gemm_nonlop
 public :: gemm_nonlop
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

   real(dp), allocatable :: projs(:, :, :)
   ! (2, npw, nprojs)
   real(dp), allocatable :: projs_r(:, :, :)
   ! (1, npw, nprojs)
   real(dp), allocatable :: projs_i(:, :, :)
   ! (1, npw, nprojs)

   real(dp), allocatable :: dprojs(:, :, :)
   ! (2, npw, nprojs*ngrads)
   real(dp), allocatable :: dprojs_r(:, :, :)
   ! (1, npw, nprojs*ngrads)
   real(dp), allocatable :: dprojs_i(:, :, :)
   ! (1, npw, nprojs*ngrads)

 end type gemm_nonlop_type
!!***

 type(gemm_nonlop_type), save, public, allocatable :: gemm_nonlop_kpt(:)
 !(nkpt)

 integer, save, public :: gemm_nonlop_ikpt_this_proc_being_treated
 !! This is oh so very crude, but I can't find any other way to do it without passing ikpt deep down to nonlop

 logical, save, public :: gemm_nonlop_use_gemm = .false.
 ! Public variable indicating whether we should call gemm_nonlop or fall back to the usual nonlop. Set to false
 ! in order not to interfere with non-GS calls to nonlop.

 logical, save, public :: gemm_nonlop_use_gemm_gpu = .false.
 ! Public variable controlled by input dataset var named use_gemm_nonlop_cuda (0 or 1).
 ! When 0, nonlop is computed by calling the regular nonlop_gpu
 ! When 1, nonlop is computed by calling gemm_nonlop_gpu

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
  else
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
      ABI_MALLOC(temp_realvec,(MAX(npwout,npwin)*nspinor*ndat))
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
      ABI_FREE(temp_realvec)
    end if
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
    end if
  end if

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
    end if

    ! opernlb (only choice=1)
    if(signs==2) then
      if(paw_opt == 3 .or. paw_opt == 4) then
        ! Get svectout from s_projections
        if(cplex == 2) then
          call abi_zgemm_2r('N', 'N', npwout, ndat*nspinor, nprojs, cone, &
&                        gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%projs, npwout, &
&                        s_projections, nprojs, czero, svectout, npwout)
        else
          ABI_MALLOC(temp_realvec,(MAX(npwout,npwin)*nspinor*ndat))
          call DGEMM('N', 'N', npwout, ndat*nspinor, nprojs, one, &
&                    gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%projs_r, npwout, &
&                    s_projections, nprojs, zero, temp_realvec, npwout)
          svectout(1,1:npwout*nspinor*ndat) = temp_realvec(1:npwout*nspinor*ndat)
          call DGEMM('N', 'N', npwout, ndat*nspinor, nprojs, one, &
&                    gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%projs_i, npwout,&
&                    s_projections, nprojs, zero, temp_realvec, npwout)
          svectout(2,1:npwout*nspinor*ndat) = temp_realvec(1:npwout*nspinor*ndat)
          ABI_FREE(temp_realvec)
        end if
        if(choice /= 7) svectout = svectout + vectin ! TODO understand this
      end if
      if(paw_opt == 0 .or. paw_opt == 1 .or. paw_opt == 4) then
        ! Get vectout from vnl_projections
        if(cplex_fac == 2) then
          call abi_zgemm_2r('N', 'N', npwout, ndat*nspinor, nprojs, cone, &
&                        gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%projs, npwout, &
&                        vnl_projections, nprojs, czero, vectout, npwout)
        else
          ABI_MALLOC(temp_realvec,(MAX(npwout,npwin)*nspinor*ndat))
          call DGEMM('N', 'N', npwout, ndat*nspinor, nprojs, one, &
&                    gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%projs_r, npwout, &
&                    vnl_projections, nprojs, zero, temp_realvec, npwout)
          vectout(1,1:npwout*nspinor*ndat) = temp_realvec(1:npwout*nspinor*ndat)
          call DGEMM('N', 'N', npwout, ndat*nspinor, nprojs, one, &
&                    gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%projs_i, npwout, &
&                    vnl_projections, nprojs, zero, temp_realvec, npwout)
          vectout(2,1:npwout*nspinor*ndat) = temp_realvec(1:npwout*nspinor*ndat)
          ABI_FREE(temp_realvec)
        end if
      end if
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

 end subroutine gemm_nonlop
!***

!----------------------------------------------------------------------

end module m_gemm_nonlop
!!***
