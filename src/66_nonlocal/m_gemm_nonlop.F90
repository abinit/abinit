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
!! Copyright (C) 2014-2024 ABINIT group (AL)
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
 use m_fstrings,    only : itoa, ftoa, sjoin
 use m_abi_linalg

 use defs_abitypes, only : MPI_type
 use m_opernlc_ylm, only : opernlc_ylm
 use m_opernla_gemm, only : opernla_gemm
 use m_opernlb_gemm, only : opernlb_gemm
 use m_opernld_ylm_allwf_cpu, only : opernld_ylm_allwf_cpu
 use m_pawcprj, only : pawcprj_type
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

   integer :: nprojs_blk
   integer :: nprojs_last_blk

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

 type(gemm_nonlop_type), save, public, allocatable, target :: gemm_nonlop_kpt(:)
 !(nkpt)

 integer, save, public :: gemm_nonlop_ikpt_this_proc_being_treated
 !! This is oh so very crude, but I can't find any other way to do it without passing ikpt deep down to nonlop

 logical, save, public :: gemm_nonlop_use_gemm = .false.
 ! Public variable indicating whether we should call gemm_nonlop or fall back to the usual nonlop. Set to false
 ! in order not to interfere with non-GS calls to nonlop.

 logical, save, public :: gemm_nonlop_is_distributed = .false.
 ! Public variable indicating whether we should gemm_nonlop operated in a distributed manner. Set to false by default
 ! but might be enabled by memory constraints or forced by user through parameters.

 integer, save, public :: gemm_nonlop_nblocks = 1
 ! Public variable indicating in how many blocks of MPI tasks should the projs arrays be ditributed.
 ! Default size 1 indicates no distribution at all.

 integer, save, public :: gemm_nonlop_block_comm = xmpi_comm_null
 ! MPI communicator for MPI tasks processing the same gemm_nonlop block for projs array distribution

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
&                            idir_pert,&
&                            compute_grad_strain,compute_grad_atom) ! Optional parameters

  integer, intent(in) :: ikpt
  integer, intent(in) :: npw, lmnmax, ntypat
  integer, intent(in) :: indlmn(:,:,:), kg_k(:,:)
  integer, intent(in) :: nattyp(ntypat)
  integer, intent(in) :: istwf_k
  logical, intent(in), optional :: compute_grad_strain,compute_grad_atom
  integer, intent(in), optional :: idir_pert
  real(dp), intent(in) :: ucvol
  real(dp), intent(in) :: ffnl_k(:,:,:,:)
  real(dp), intent(in) :: ph3d_k(:,:,:)
  real(dp), intent(in) :: kpt_k(:)
  real(dp), intent(in), target :: kpg_k(:,:)

  integer :: nprojs,ndprojs,ngrads

  integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
  integer :: itypat, ilmn, nlmn, ia, iaph3d, igrad, shift, shift_grad
  integer :: il, ipw, idir, idir1, idir2, nkpg_local, ffnl_dir, dimffnl
  logical :: parity,my_compute_grad_strain,my_compute_grad_atom,compute_pert
  real(dp):: wt
  real(dp),allocatable :: atom_projs(:,:,:), atom_dprojs(:,:,:,:), temp(:)
  real(dp),pointer :: kpg(:,:)
  integer :: rank, nprocs, ierr
  integer :: nprojs_blk, nprojs_last_blk, nprojs_my_blk
  integer :: lmn_beg,ibeg,iend,shift_do,nlmn_o,lmn_grad_beg
  integer :: idir_beg,idir_end,idir_pert_
  logical :: is_last_rank
  real(dp),allocatable :: dprojs_tmp(:,:,:),dprojs_r_tmp(:,:,:),dprojs_i_tmp(:,:,:)

! *************************************************************************

  my_compute_grad_strain=.false. ; if (present(compute_grad_strain)) my_compute_grad_strain=compute_grad_strain
  my_compute_grad_atom=.false. ; if (present(compute_grad_atom)) my_compute_grad_atom=compute_grad_atom
  ABI_CHECK(size(ph3d_k)>0,'nloalg(2)<0 not compatible with use_gemm_nonlop=1!')
!  ABI_CHECK((.not.my_compute_grad_strain).or.size(kpg_k)>0,"kpg_k should be allocated to compute gradients!")
!  ABI_CHECK((.not.my_compute_grad_atom).or.size(kpg_k)>0,"kpg_k should be allocated to compute gradients!")
  idir_pert_=0; if(present(idir_pert)) idir_pert_=idir_pert
  compute_pert=.false.; if(present(idir_pert) .and. (.not. my_compute_grad_atom .and. .not. my_compute_grad_strain)) compute_pert=.true. !FIXME Need refacto

  iaph3d = 1
  wt=four_pi/sqrt(ucvol)
  dimffnl = size(ffnl_k, dim=2)
  ffnl_dir=1; if(dimffnl>2) ffnl_dir=idir_pert_

  ABI_MALLOC(atom_projs, (2, npw, lmnmax))
  if (my_compute_grad_strain .and. .not. present(idir_pert)) then
    ndprojs = 3
    ABI_MALLOC(atom_dprojs, (2, npw, ndprojs, lmnmax))
  else if(compute_pert .or. my_compute_grad_strain) then
    ndprojs = 1
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
  if(present(idir_pert)) then
    ngrads=ngrads+1
  else
    if (my_compute_grad_strain) ngrads=ngrads+6
    if (my_compute_grad_atom) ngrads=ngrads+3
  end if
  if (nprojs>0) gemm_nonlop_kpt(ikpt)%nprojs = nprojs
  if (ngrads>0) gemm_nonlop_kpt(ikpt)%ngrads = ngrads

  if(gemm_nonlop_is_distributed) then
    rank = xmpi_comm_rank(xmpi_world); nprocs = xmpi_comm_size(xmpi_world)
    is_last_rank = (rank==nprocs-1)
    if(gemm_nonlop_block_comm/=xmpi_comm_null) call xmpi_comm_free(gemm_nonlop_block_comm)
    if(gemm_nonlop_nblocks==1) gemm_nonlop_nblocks=nprocs
    write(std_out,*) "Splitting on ", gemm_nonlop_nblocks
    call xmpi_comm_split(xmpi_world, rank/(nprocs/gemm_nonlop_nblocks), rank, gemm_nonlop_block_comm, ierr)
    if(ierr/=0) ABI_BUG("Bug split!")
  end if
  rank = xmpi_comm_rank(gemm_nonlop_block_comm); nprocs = xmpi_comm_size(gemm_nonlop_block_comm)
  is_last_rank = (rank==nprocs-1)

  if(gemm_nonlop_is_distributed) then
    nprojs_blk = nprojs / nprocs
    nprojs_last_blk = nprojs_blk + modulo(nprojs,nprojs_blk)
    gemm_nonlop_kpt(ikpt)%nprojs_blk = nprojs_blk
    gemm_nonlop_kpt(ikpt)%nprojs_last_blk = nprojs_last_blk
    if(is_last_rank) then
      nprojs_my_blk = nprojs_last_blk
    else
      nprojs_my_blk = nprojs_blk
    end if
  else
    nprojs_last_blk = nprojs
    nprojs_my_blk = nprojs
    nprojs_blk = nprojs
  end if

  if(gemm_nonlop_is_distributed .and. ngrads>0) then
    !Temporary buffer to help distribute dprojs correctly (and easily)
    if(istwf_k <= 1) then
      ABI_MALLOC(dprojs_tmp, (2, npw, nprojs_my_blk*ngrads*2))
    else
      ABI_MALLOC(dprojs_r_tmp, (1, npw, nprojs_my_blk*ngrads*2))
      ABI_MALLOC(dprojs_i_tmp, (1, npw, nprojs_my_blk*ngrads*2))
    end if
  end if

  if(istwf_k <= 1) then
    ABI_MALLOC(gemm_nonlop_kpt(ikpt)%projs, (2, npw, nprojs_last_blk))
    gemm_nonlop_kpt(ikpt)%projs = zero
    if(ngrads>0) then
      ABI_MALLOC(gemm_nonlop_kpt(ikpt)%dprojs, (2, npw, nprojs_last_blk*ngrads))
      gemm_nonlop_kpt(ikpt)%dprojs = zero
    end if
  else
    ABI_MALLOC(gemm_nonlop_kpt(ikpt)%projs_r, (1, npw, nprojs_last_blk))
    ABI_MALLOC(gemm_nonlop_kpt(ikpt)%projs_i, (1, npw, nprojs_last_blk))
    gemm_nonlop_kpt(ikpt)%projs_r = zero
    gemm_nonlop_kpt(ikpt)%projs_i = zero
    if(ngrads>0) then
      ABI_MALLOC(gemm_nonlop_kpt(ikpt)%dprojs_r, (1, npw, nprojs_last_blk*ngrads))
      ABI_MALLOC(gemm_nonlop_kpt(ikpt)%dprojs_i, (1, npw, nprojs_last_blk*ngrads))
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

  if(gemm_nonlop_is_distributed) then
    ibeg = rank*nprojs_blk+1
    iend = ibeg+nprojs_my_blk
    shift_do = 0
    lmn_grad_beg = -1
  end if

  shift = 0 ; shift_grad = 0
  lmn_beg = 1

  do itypat = 1, ntypat
    nlmn = count(indlmn(3,:,itypat)>0)
    nlmn_o = nlmn

    do ia = 1, nattyp(itypat)

      ! In distributed mode, loops are skipped until reach the section
      ! of "ilmn" to be stored by local rank
      if(gemm_nonlop_is_distributed) then
        if(shift_do+nlmn < ibeg) then
          shift_do = shift_do + nlmn
          iaph3d = iaph3d + 1
          cycle
        end if

        lmn_beg = max(1,ibeg-shift_do)
        if(lmn_grad_beg==-1) lmn_grad_beg = (lmn_beg-1)*ngrads
        if(shift_do+nlmn > iend) nlmn = iend - shift_do - 1
      end if

      !! build atom_projs, from opernlb
      !! P = 4pi/sqrt(ucvol)* conj(diag(ph3d)) * ffnl * diag(parity), with parity = (-i)^l

      ! start from 4pi/sqrt(ucvol)*ffnl
      ! atom_projs(1, :, lmn_beg:nlmn) = four_pi/sqrt(ham%ucvol) * ham%ffnl_k(:, 1, lmn_beg:nlmn)
      ! TODO vectorize (DCOPY with stride)
      atom_projs(:,:,:) = zero
      do ipw=1, npw
        atom_projs(1,ipw, 1:nlmn_o) = wt * ffnl_k(ipw, 1, 1:nlmn_o, itypat)
      end do
      if (ndprojs>0) atom_dprojs(:,:,:,:) = zero
      if (my_compute_grad_strain .and. .not. present(idir_pert)) then
        do ipw=1, npw
          atom_dprojs(1,ipw, 1:ndprojs, 1:nlmn_o) = wt * ffnl_k(ipw, 2:ndprojs+1, 1:nlmn_o, itypat)
        end do
      end if
      if (compute_pert .or. (my_compute_grad_strain .and. present(idir_pert)) ) then
        do ipw=1, npw
          atom_dprojs(1,ipw, 1, 1:nlmn_o) = wt * ffnl_k(ipw, 1+ffnl_dir, 1:nlmn_o, itypat)
        end do
      end if


      ! multiply by (-i)^l
      do ilmn=1,nlmn_o
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
      do ilmn=1,nlmn_o
        temp = atom_projs(1, :, ilmn)
        atom_projs(1, :, ilmn) = atom_projs(1, :, ilmn) * ph3d_k(1, :, iaph3d) &
&                              + atom_projs(2, :, ilmn) * ph3d_k(2, :, iaph3d)
        atom_projs(2, :, ilmn) = atom_projs(2, :, ilmn) * ph3d_k(1, :, iaph3d) &
&                              - temp                   * ph3d_k(2, :, iaph3d)
      end do
      if (ndprojs>0) then
        do ilmn=1,nlmn_o
          do idir=1,ndprojs
            temp = atom_dprojs(1, :, idir,ilmn)
            atom_dprojs(1, :, idir,ilmn) = atom_dprojs(1, :, idir,ilmn) * ph3d_k(1, :, iaph3d) &
&                                        + atom_dprojs(2, :, idir,ilmn) * ph3d_k(2, :, iaph3d)
            atom_dprojs(2, :, idir,ilmn) = atom_dprojs(2, :, idir,ilmn) * ph3d_k(1, :, iaph3d) &
&                                        - temp                         * ph3d_k(2, :, iaph3d)
          end do
        end do
      end if

      !! atom_projs is built, copy to projs

      if(istwf_k <= 1) then
        gemm_nonlop_kpt(ikpt)%projs(1:2, :, shift+1:shift+(nlmn-(lmn_beg-1))) = atom_projs(:, :, lmn_beg:nlmn)
      else ! istwf_k>1
        gemm_nonlop_kpt(ikpt)%projs_r(1, :, shift+1:shift+(nlmn-(lmn_beg-1))) = atom_projs(1, :, lmn_beg:nlmn)
        gemm_nonlop_kpt(ikpt)%projs_i(1, :, shift+1:shift+(nlmn-(lmn_beg-1))) = atom_projs(2, :, lmn_beg:nlmn)
      end if

      !! Handling dprojs

      ! In distributed case, we compute values for all nlmn*ngrads in dprojs_tmp.
      ! dprojs will be filled by "cutting" relevant values out
      if(ngrads>0) then
        if(gemm_nonlop_is_distributed) then
          if(istwf_k <= 1) then
            igrad=0
            if(my_compute_grad_strain) then
              do idir=1,6
                idir1=alpha(idir);idir2=beta(idir)
                do ilmn=1,nlmn_o
                  do ipw=1,npw
                    dprojs_tmp(1:2, ipw, shift_grad+nlmn_o*igrad+ilmn) = &
  &                  -half*(atom_dprojs(1:2, ipw, idir1, ilmn)*kpg(ipw,idir2) &
  &                        +atom_dprojs(1:2, ipw, idir2, ilmn)*kpg(ipw,idir1))
                  end do
                end do
                igrad=igrad+1
              end do
            end if
            if(my_compute_grad_atom) then
              do idir=1,3
                do ilmn=1,nlmn_o
                  do ipw=1,npw
                    dprojs_tmp(1, ipw, shift_grad+nlmn_o*igrad+ilmn) = &
  &                  +atom_projs(2, ipw, ilmn)*kpg(ipw,idir)*two_pi
                    dprojs_tmp(2, ipw, shift_grad+nlmn_o*igrad+ilmn) = &
  &                  -atom_projs(1, ipw, ilmn)*kpg(ipw,idir)*two_pi
                  end do
                end do
                igrad=igrad+1
              end do
            end if

          else ! istwf_k>1
            igrad=0
            if(my_compute_grad_strain) then
              do idir=1,6
                idir1=alpha(idir);idir2=beta(idir)
                do ilmn=1,nlmn_o
                  do ipw=1,npw
                    dprojs_r_tmp(1, ipw, shift_grad+nlmn_o*igrad+ilmn) = &
  &                  -half*(atom_dprojs(1, ipw, idir1, ilmn)*kpg(ipw,idir2) &
  &                        +atom_dprojs(1, ipw, idir2, ilmn)*kpg(ipw,idir1))

                    dprojs_i_tmp(1, ipw, shift_grad+nlmn_o*igrad+ilmn) = &
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
                    dprojs_r_tmp(1, ipw, shift_grad+nlmn_o*igrad+ilmn) = &
  &                  +atom_projs(2, ipw, ilmn)*kpg(ipw,idir)*two_pi
                    dprojs_i_tmp(1, ipw, shift_grad+nlmn_o*igrad+ilmn) = &
  &                  -atom_projs(1, ipw, ilmn)*kpg(ipw,idir)*two_pi
                  end do
                end do
                igrad=igrad+1
              end do
            end if
          end if ! istwf_k

        else ! non-distributed case

          if(istwf_k <= 1) then
            igrad=0
            if(my_compute_grad_strain .and. .not. present(idir_pert)) then
              do idir=1,6
                idir1=alpha(idir);idir2=beta(idir)
                do ilmn=lmn_beg,nlmn
                  do ipw=1,npw
                    gemm_nonlop_kpt(ikpt)%dprojs(1:2, ipw, shift_grad+nlmn*igrad+ilmn) = &
  &                  -half*(atom_dprojs(1:2, ipw, idir1, ilmn)*kpg(ipw,idir2) &
  &                        +atom_dprojs(1:2, ipw, idir2, ilmn)*kpg(ipw,idir1))
                  end do
                end do
                igrad=igrad+1
              end do
            end if
            if(my_compute_grad_atom .and. .not. present(idir_pert)) then
              do idir=1,3
                do ilmn=lmn_beg,nlmn
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
            if(my_compute_grad_atom .and. present(idir_pert)) then
              !FIXME Ugly workaround so it works on many tasks
              if(rank/=0) write(std_out,*) idir_pert_
              do ilmn=lmn_beg,nlmn
                do ipw=1,npw
                  gemm_nonlop_kpt(ikpt)%dprojs(1, ipw, shift_grad+nlmn*igrad+ilmn) = &
&                  +atom_projs(2, ipw, ilmn)*kpg(ipw,idir_pert_)*two_pi
                  gemm_nonlop_kpt(ikpt)%dprojs(2, ipw, shift_grad+nlmn*igrad+ilmn) = &
&                  -atom_projs(1, ipw, ilmn)*kpg(ipw,idir_pert_)*two_pi
                end do
              end do
              igrad=igrad+1
            end if
            if(my_compute_grad_strain .and. present(idir_pert)) then
              do ilmn=lmn_beg,nlmn
                do ipw=1,npw
                  gemm_nonlop_kpt(ikpt)%dprojs(1, ipw, shift_grad+nlmn*igrad+ilmn) = &
  &                -atom_dprojs(1, ipw, 1, ilmn)
                  gemm_nonlop_kpt(ikpt)%dprojs(2, ipw, shift_grad+nlmn*igrad+ilmn) = &
  &                -atom_dprojs(2, ipw, 1, ilmn)
                end do
              end do
              igrad=igrad+1
            end if
            if(compute_pert) then
              do ilmn=lmn_beg,nlmn
                do ipw=1,npw
                  gemm_nonlop_kpt(ikpt)%dprojs(1, ipw, shift_grad+nlmn*igrad+ilmn) = &
  &                +atom_dprojs(1, ipw, 1, ilmn)
                  gemm_nonlop_kpt(ikpt)%dprojs(2, ipw, shift_grad+nlmn*igrad+ilmn) = &
  &                +atom_dprojs(2, ipw, 1, ilmn)
                end do
              end do
              igrad=igrad+1
            end if

          else ! istwf_k>1
            igrad=0
            if(my_compute_grad_strain) then
              do idir=1,6
                idir1=alpha(idir);idir2=beta(idir)
                do ilmn=lmn_beg,nlmn
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
                do ilmn=lmn_beg,nlmn
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
          end if ! istwf_k
        end if ! gemm_nonlop_is_distributed
      end if ! ngrads

      iaph3d = iaph3d + 1
      shift_grad = shift_grad + ngrads*nlmn_o

      if(gemm_nonlop_is_distributed) then
        shift = shift + nlmn - (lmn_beg-1)
        shift_do = shift_do + nlmn
        if(shift_do >= iend - 1) exit
      else
        shift = shift + nlmn
      end if
    end do
    if(gemm_nonlop_is_distributed .and. shift_do >= iend - 1) exit
  end do

  ! Filling dprojs by extracting values from dprojs_tmp
  if(gemm_nonlop_is_distributed .and. ngrads>0) then
    shift_grad = lmn_grad_beg
    if(istwf_k <= 1) then
      gemm_nonlop_kpt(ikpt)%dprojs(1:2, 1:npw, 1:ngrads*nprojs_my_blk) = &
      &      dprojs_tmp(1:2, 1:npw, shift_grad+1:shift_grad+ngrads*nprojs_my_blk)
    else
      gemm_nonlop_kpt(ikpt)%dprojs_r(1, 1:npw, 1:ngrads*nprojs_my_blk) = &
      &      dprojs_r_tmp(1, 1:npw, shift_grad+1:shift_grad+ngrads*nprojs_my_blk)
      gemm_nonlop_kpt(ikpt)%dprojs_i(1, 1:npw, 1:ngrads*nprojs_my_blk) = &
      &      dprojs_i_tmp(1, 1:npw, shift_grad+1:shift_grad+ngrads*nprojs_my_blk)
    end if
  end if

  if(gemm_nonlop_is_distributed .and. ngrads>0) then
    if(istwf_k <= 1) then
      ABI_FREE(dprojs_tmp)
    else
      ABI_FREE(dprojs_r_tmp)
      ABI_FREE(dprojs_i_tmp)
    end if
  end if
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
!! [gpu_option] = GPU implementation to use, i.e. cuda, openMP, ... (0=not using GPU)
!!
!! SOURCE

 subroutine gemm_nonlop(atindx1,choice,cpopt,cprjin,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,&
&                 enl,enlout,ffnlin,ffnlout,gmet,gprimd,idir,indlmn,istwf_k,&
&                 kgin,kgout,kpgin,kpgout,kptin,kptout,lambda,lmnmax,matblk,mgfft,&
&                 mpi_enreg,mpsang,mpssoang,natom,nattyp,ndat,ngfft,nkpgin,nkpgout,nloalg,&
&                 nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO,paw_opt,phkxredin,&
&                 phkxredout,ph1d,ph3din,ph3dout,signs,sij,svectout,&
&                 tim_nonlop,ucvol,useylm,vectin,vectout,atom_proj_shift,&
&                 vectproj,gpu_option)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: choice,cpopt,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,idir
  integer,intent(in) :: istwf_k,lmnmax,matblk,mgfft,mpsang,mpssoang,natom,ndat,nkpgin
  integer,intent(in) :: nkpgout,nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO
  integer,intent(in) :: paw_opt,signs,tim_nonlop,useylm,atom_proj_shift
  integer,optional,intent(in) :: gpu_option
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
  real(dp),intent(inout),optional, ABI_CONTIGUOUS target :: vectproj(:,:,:)
  type(pawcprj_type),intent(inout) :: cprjin(natom,nspinor*((cpopt+5)/5)*ndat)

  ! locals
  complex(dpc), parameter :: cminusone  = (-1._dp,0._dp)
  integer :: ii, ia, idat, igrad, nprojs, ngrads, shift, iatom, nlmn, ierr, ibeg, iend, ikpt
  integer :: cplex, cplex_enl, cplex_fac, proj_shift, grad_shift
  integer :: projs_beg,projs_end
  integer :: nnlout_test
  integer :: iatm, ndgxdt, ndgxdtfac, nd2gxdt, nd2gxdtfac, optder, itypat, ilmn
  integer :: cplex_dgxdt(1), cplex_d2gxdt(1)
  logical :: local_vectproj
  real(dp) :: dgxdt_dum_in(1,1,1,1,1), dgxdt_dum_out(1,1,1,1,1),dgxdt_dum_out2(1,1,1,1,1)
  real(dp) :: d2gxdt_dum_in(1,1,1,1,1), d2gxdt_dum_out(1,1,1,1,1),d2gxdt_dum_out2(1,1,1,1,1)
  real(dp), allocatable :: sij_typ(:)
  real(dp),  ABI_CONTIGUOUS pointer :: projections(:,:,:)
  real(dp), allocatable :: s_projections(:,:,:), vnl_projections(:,:,:)
  real(dp), allocatable :: dprojections(:,:,:), temp_realvec(:)
  real(dp), allocatable :: s_dprojections(:,:,:), vnl_dprojections(:,:,:)
  integer :: nprojs_my_blk
  integer :: rank, nprocs
  logical :: is_last
  real(dp), allocatable :: projs_local(:,:,:), dprojs_local(:,:,:)
  real(dp), allocatable :: projs_recv(:,:,:), dprojs_recv(:,:,:)
  integer :: ngrads_tmp

! *************************************************************************

  ! We keep the same interface as nonlop, but we don't use many of those
  ABI_UNUSED((/ffnlin,ffnlout,gmet,kpgin,kpgout/))
  ABI_UNUSED((/ph1d(1,1),ph3din,ph3dout/))
  ABI_UNUSED((/phkxredin,phkxredout,ucvol/))
  ABI_UNUSED((/mgfft,mpsang,mpssoang/))
  ABI_UNUSED((/kptin,kptout/))
  ABI_UNUSED((/idir,nloalg,ngfft,kgin,kgout,ngfft,only_SO,tim_nonlop,gpu_option/))

  ! Check supported options
  if (.not.gemm_nonlop_use_gemm) then
    ABI_BUG('computation not prepared for gemm_nonlop use!')
  end if
  if ( (choice>3.and.choice/=7.and.choice/=5.and.choice/=51.and.signs==2) .or. &
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

  ikpt=gemm_nonlop_ikpt_this_proc_being_treated
  cplex=2;if (istwf_k>1) cplex=1
  cplex_enl=1;if (paw_opt>0) cplex_enl=2*dimenl1/(lmnmax*(lmnmax+1)) ! is enl complex?
  cplex_fac=max(cplex,dimekbq)
  if ((nspinortot==2.or.cplex_enl==2).and.paw_opt>0.and.choice/=7) cplex_fac=2 ! is vnl_projections complex?

  ! The number of projectors used for computation may vary among
  ! nonlop calls, from computing on all atoms to a select one for
  ! some perturbations.
  ! In such cases, projs arrays must be recomputed
  nprojs=0
  do itypat=1,ntypat
    nprojs = nprojs + count(indlmn(3,:,itypat)>0)*nattyp(itypat)
  end do
  projs_beg=1; projs_end=nprojs;
  if(choice==2) then
    projs_beg=atom_proj_shift+1
    projs_end=projs_beg+nprojs
  end if

  ngrads = gemm_nonlop_kpt(ikpt)%ngrads
  nprojs_my_blk = gemm_nonlop_kpt(ikpt)%nprojs_blk

  if(gemm_nonlop_is_distributed) then
    rank = xmpi_comm_rank(gemm_nonlop_block_comm); nprocs = xmpi_comm_size(gemm_nonlop_block_comm)
    is_last = (rank==nprocs-1)
    if(is_last) nprojs_my_blk = gemm_nonlop_kpt(ikpt)%nprojs_last_blk
  end if

  if(choice==1) ngrads=0
  if(signs==1) then
    ABI_CHECK(ngrads>=3.or.choice/=2 ,"Bad allocation in gemm_nonlop (2)!")
    ABI_CHECK(ngrads>=6.or.choice/=3 ,"Bad allocation in gemm_nonlop (3)!")
    ABI_CHECK(ngrads>=9.or.choice/=23,"Bad allocation in gemm_nonlop (23)!")
  else if(signs==2) then
    ABI_CHECK(ngrads==1.or.choice==1 ,"Bad allocation in gemm_nonlop !")
  end if

  ! If vectproj is provided, use it for further calculations, use allocated array otherwise
  local_vectproj=.false.
  if(PRESENT(vectproj)) then
    if(size(vectproj)>1) local_vectproj=.true.
  end if
  if (local_vectproj) then
    projections => vectproj
  end if

  ! These will store the non-local factors for vectin, svectout and vectout respectively
  if(.not. local_vectproj) then
    ABI_MALLOC(projections,(cplex, nprojs,nspinor*ndat))
  end if
  ABI_MALLOC(s_projections,(cplex, nprojs,nspinor*ndat))
  ABI_MALLOC(vnl_projections,(cplex_fac, nprojs,nspinor*ndat))

  if(gemm_nonlop_is_distributed) then
    ABI_MALLOC(projs_recv,   (cplex, npwin, gemm_nonlop_kpt(ikpt)%nprojs_last_blk))
    ABI_MALLOC(projs_local,   (cplex, npwin, gemm_nonlop_kpt(ikpt)%nprojs_last_blk))
    if (ngrads>0) then
      ABI_MALLOC(dprojs_recv,   (cplex, npwin, ngrads*gemm_nonlop_kpt(ikpt)%nprojs_last_blk))
      ABI_MALLOC(dprojs_local,   (cplex, npwin, ngrads*gemm_nonlop_kpt(ikpt)%nprojs_last_blk))
    end if
  end if

  if(.not. local_vectproj) projections = zero
  s_projections = zero
  vnl_projections = zero
  !if (ngrads>0) then
  ngrads_tmp=1; if (ngrads>0) ngrads_tmp=ngrads
    ABI_MALLOC(dprojections,(cplex, ngrads_tmp*nprojs,nspinor*ndat))
    ABI_MALLOC(s_dprojections,(cplex, ngrads_tmp*nprojs,nspinor*ndat))
    ABI_MALLOC(vnl_dprojections,(cplex_fac, ngrads_tmp*nprojs,nspinor*ndat))
    dprojections(:,:,:) = zero
    s_dprojections(:,:,:) = zero
    vnl_dprojections(:,:,:) = zero
  !end if

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
    if(.not. local_vectproj) then
      !This use-case is extremely painful for GEMM nonlop performance.
      !vectproj paramter should always be provided to avoid it

      ! retrieve from cprjin
      do idat=1, ndat*nspinor
        shift = 0
        do iatom = 1, natom
          nlmn = cprjin(iatom, idat)%nlmn
          projections(1:cplex, shift+1:shift+nlmn, idat) = cprjin(iatom, idat)%cp(1:cplex, 1:nlmn)
          shift = shift + nlmn
        end do
      end do
    end if
    if(cpopt==4.and.allocated(dprojections)) then
      ABI_CHECK(cprjin(1,1)%ncpgr>=ngrads,"cprjin%ncpgr not correct! (1)")
      do idat=1, ndat*nspinor
        shift = 0
        do iatom = 1, natom
          nlmn  = cprjin(iatom, idat)%nlmn
          do igrad=1,ngrads
            dprojections(1:cplex, shift+1:shift+nlmn, idat) = &
&                   cprjin(iatom, idat)%dcp(1:cplex,igrad,1:nlmn)
            shift = shift + nlmn
          end do
        end do
      end do
    end if
  end if ! cpopt

  if(cpopt<=1.or.(cpopt<=3.and.(choice==2 .or. choice==3 .or. choice==5.or.choice==51))) then

    call opernla_gemm(choice,cplex,cplex_dgxdt,cplex_d2gxdt,d2gxdt_dum_in,dprojections,projections,&
    &       idir,istwf_k,mpi_enreg,nd2gxdt,ngrads,&
    &       npwin,nspinor,signs,ndat,rank,&
    &       cpopt,nprocs,&
    &       nprojs,gemm_nonlop_kpt(ikpt)%nprojs_blk,nprojs_my_blk,gemm_nonlop_kpt(ikpt)%nprojs_last_blk,&
    &       vectin,gemm_nonlop_kpt(ikpt)%projs,gemm_nonlop_kpt(ikpt)%dprojs,&
    &       gemm_nonlop_kpt(ikpt)%projs_r,gemm_nonlop_kpt(ikpt)%projs_i,&
    &       gemm_nonlop_kpt(ikpt)%dprojs_r,gemm_nonlop_kpt(ikpt)%dprojs_i,&
    &       temp_realvec,&
    &       projs_local,projs_recv,dprojs_local,dprojs_recv,&
    &       gemm_nonlop_is_distributed)


    if(cpopt >= 0) then
      if(.not. local_vectproj) then
        !This use-case is extremely painful for GEMM nonlop performance.
        !vectproj paramter should always be provided to avoid it

        ! store in cprjin
        do idat=1, ndat*nspinor
          shift = 0
          do iatom = 1, natom
            nlmn = cprjin(iatom, idat)%nlmn
            cprjin(iatom, idat)%cp(1:cplex, 1:nlmn) = projections(1:cplex, shift+1:shift+nlmn, idat)
            shift = shift + nlmn
          end do
        end do
      end if
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
      ndgxdt = ngrads
      ndgxdtfac = ngrads
      nd2gxdt = 0
      nd2gxdtfac = 0
      optder = 0;if (ndgxdtfac>0 .and. signs == 2) optder=1
      cplex_dgxdt(:) = 1;
      if(ndgxdt > 0) then
       if (choice==5.or.choice==51) cplex_dgxdt(:) = 2
      end if

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
          call opernlc_ylm(atindx1,cplex,cplex_dgxdt,cplex_d2gxdt,cplex_enl,cplex_fac,&
&         dprojections(:, ibeg:iend, 1+nspinor*(idat-1):nspinor*idat),&
&         vnl_dprojections(:, ibeg:iend, 1+nspinor*(idat-1):nspinor*idat),&
&         s_dprojections(:, ibeg:iend, 1+nspinor*(idat-1):nspinor*idat),&
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

    ! opernlb
    if(signs==2) then
      call opernlb_gemm(choice,cplex,cplex_dgxdt,cplex_d2gxdt,cplex_fac,&
      &       d2gxdt_dum_out,d2gxdt_dum_out,&
      &       vnl_dprojections,s_dprojections,&
      &       vnl_projections,s_projections,&
      &       idir,istwf_k,mpi_enreg,nd2gxdt,ngrads,&
      &       npwin,nspinor,signs,ndat,rank,&
      &       cpopt,nprocs,paw_opt,&
      &       nprojs,gemm_nonlop_kpt(ikpt)%nprojs_blk,nprojs_my_blk,gemm_nonlop_kpt(ikpt)%nprojs_last_blk,&
      &       vectin,vectout,svectout,&
      &       gemm_nonlop_kpt(ikpt)%projs(:,:,projs_beg:projs_end),&
      &       gemm_nonlop_kpt(ikpt)%dprojs(:,:,projs_beg:projs_end),&
      &       gemm_nonlop_kpt(ikpt)%projs_r,gemm_nonlop_kpt(ikpt)%projs_i,&
      &       gemm_nonlop_kpt(ikpt)%dprojs_r,gemm_nonlop_kpt(ikpt)%dprojs_i,&
      &       temp_realvec,&
      &       projs_local,projs_recv,dprojs_local,dprojs_recv,&
      &       gemm_nonlop_is_distributed)
    end if

    ! opernld
    if(signs==1) then
      call opernld_ylm_allwf_cpu(choice,cplex,cplex_fac,&
      &       dprojections,vnl_dprojections,s_dprojections,d2gxdt_dum_in,&
      &       enlout,projections,vnl_projections,s_projections,&
      &       ndat,nd2gxdt,ndgxdt,&
      &       ndgxdtfac,indlmn,ntypat,lmnmax,nprojs,nnlout,nspinor,paw_opt,&
      &       gprimd,nattyp,mpi_enreg)
    end if !opernld

  end if ! choice>0

! Release memory
  if(.not. local_vectproj) then
    ABI_FREE(projections)
  end if
  ABI_FREE(s_projections)
  ABI_FREE(vnl_projections)
  if(gemm_nonlop_is_distributed) then
    ABI_FREE(projs_local)
    ABI_FREE(projs_recv)
    if (signs==1.and.ngrads>0) then
      ABI_FREE(dprojs_local)
      ABI_FREE(dprojs_recv)
    end if
  end if
  if (allocated(dprojections)) then
    ABI_FREE(dprojections)
  end if
  if (allocated(s_dprojections)) then
    ABI_FREE(s_dprojections)
  end if
  if (allocated(vnl_dprojections)) then
    ABI_FREE(vnl_dprojections)
  end if
  if (allocated(temp_realvec)) then
    ABI_FREE(temp_realvec)
  end if

 end subroutine gemm_nonlop
!***

end module m_gemm_nonlop
!!***
