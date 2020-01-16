!!****m* ABINIT/m_invovl
!! NAME
!!  m_invovl
!!
!! FUNCTION
!!  Provides functions to invert the overlap matrix S. Used by Chebyshev in PAW
!!  See paper by A. Levitt and M. Torrent for details
!!  S = 1 + projs * s_projs * projs'
!!  S^-1 = 1 + projs * inv_s_projs * projs', with
!!  inv_s_projs = - (s_projs^-1 + projs'*projs)^-1
!!
!! COPYRIGHT
!! Copyright (C) 2013-2019 ABINIT group (AL)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_invovl

 use defs_basis
 use m_errors
 use m_xmpi
 use m_abicore

 use defs_abitypes, only : mpi_type
 use m_time,        only : timab
 use m_hamiltonian, only : gs_hamiltonian_type
 use m_bandfft_kpt, only : bandfft_kpt_get_ikpt
 use m_pawcprj,     only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_axpby
 use m_nonlop,      only : nonlop
 use m_prep_kgb,    only : prep_nonlop

 implicit none

 private

!public procedures.
 public :: init_invovl
 public :: make_invovl
 public :: apply_invovl
 public :: destroy_invovl

!!***

!!****t* m_invovl/invovl_kpt_type
!! NAME
!! invovl_kpt_type
!!
!! FUNCTION
!! Contains information needed to invert the overlap matrix S
!!
!! SOURCE

 type, public :: invovl_kpt_type

   integer :: nprojs
   !  total number of projectors
   !    nlmn for a specific atom = count(indlmn(3,:,itypat)>0)
   !  A value of -1 means that the following arrays are not allocated

   real(dp), allocatable :: gram_projs(:,:,:)
   ! gram_projs(cplx, nprojs, nprojs)
   ! projs' * projs

   real(dp), allocatable :: inv_sij(:, :,:,:)
   ! inv_sij(cplx, lmnmax, lmnmax, ntypat)
   ! inverse of ham%sij

   real(dp), allocatable :: inv_s_approx(:,:,:,:)
   ! inv_s_approx(cplx, lmnmax, lmnmax, ntypat)
   ! preconditionner

end type invovl_kpt_type
!!***

 type(invovl_kpt_type), public,save,allocatable, target :: invovl_kpt(:)

CONTAINS

!!****f* m_invovl/init_invovl
!! NAME
!! init_invovl
!!
!! FUNCTION
!! Initalization of the invovl_kpt array
!!
!! INPUTS
!! nkpt= number of k-points
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      dsymm,zhemm
!!
!! SOURCE

 subroutine init_invovl(nkpt)

  integer, intent(in) :: nkpt
  integer :: ikpt

! *************************************************************************

  ABI_DATATYPE_ALLOCATE(invovl_kpt, (nkpt))
  ! TODO add cycling if kpt parallelism
  do ikpt=1,nkpt
    invovl_kpt(ikpt)%nprojs = -1
  end do

 end subroutine init_invovl
!!***

!!****f* m_invovl/destroy_invovl
!! NAME
!! destroy_invovl
!!
!! FUNCTION
!! Destruction of the invovl_kpt array
!!
!! INPUTS
!! nkpt= number of k-points
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      dsymm,zhemm
!!
!! SOURCE
 subroutine destroy_invovl(nkpt)

  integer, intent(in) :: nkpt
  integer :: ikpt

! *************************************************************************

  ! TODO add cycling if kpt parallelism
  do ikpt=1,nkpt
    if(invovl_kpt(ikpt)%nprojs == -1) then
      ! write(0, *) 'ERROR invovl_kpt is unallocated'
      cycle
    end if

    ABI_DEALLOCATE(invovl_kpt(ikpt)%gram_projs)
    ABI_DEALLOCATE(invovl_kpt(ikpt)%inv_sij)
    ABI_DEALLOCATE(invovl_kpt(ikpt)%inv_s_approx)
    invovl_kpt(ikpt)%nprojs = -1
  end do

  ABI_DATATYPE_DEALLOCATE(invovl_kpt)

 end subroutine destroy_invovl
!!***

!!****f* m_invovl/make_invovl
!! NAME
!! make_invovl
!!
!! FUNCTION
!! Builds of the invovl structure
!!
!! INPUTS
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!      dsymm,zhemm
!!
!! SOURCE

subroutine make_invovl(ham, dimffnl, ffnl, ph3d, mpi_enreg)

 use m_abi_linalg
 implicit none

 type(gs_hamiltonian_type),intent(in) :: ham
 integer, intent(in) :: dimffnl
 real(dp),intent(in) :: ffnl(ham%npw_k,dimffnl,ham%lmnmax,ham%ntypat)
 real(dp),intent(in) :: ph3d(2,ham%npw_k,ham%matblk)
 type(mpi_type) :: mpi_enreg

 real(dp) :: atom_projs(2, ham%npw_k, ham%lmnmax)
 real(dp) :: temp(ham%npw_k)
 complex(dpc), allocatable :: work(:)
 real(dp), allocatable :: projs(:,:,:)
 real(dp), allocatable :: gram_proj(:,:,:)
 integer, allocatable :: ipiv(:)

 integer :: itypat, ilmn, nlmn, jlmn, ia, iaph3d, shift
 integer :: il, ilm, jlm, ipw, info, ierr, cplx
 integer :: ikpt_this_proc
 logical :: parity
 real(dp) :: tsec(2)
 character(len=500) :: message
 character :: blas_transpose

 type(invovl_kpt_type), pointer :: invovl
 integer :: array_nprojs_pp(mpi_enreg%nproc_fft)
 integer :: iproc, slice_size
 real(dp), allocatable :: gramwork(:,:,:)

 integer, parameter :: timer_mkinvovl = 1620, timer_mkinvovl_build_d = 1621, timer_mkinvovl_build_ptp = 1622

! *************************************************************************

 !! S = 1 + PDP', so
 !! S^-1 = 1 + P inv_s_projs P', with
 !! inv_s_projs = - (D^-1 + P'*P)^-1

 if(ham%istwf_k == 1) then
   cplx = 2
   blas_transpose = 'c'
 else
   cplx = 1
   blas_transpose = 't'
 end if

 ikpt_this_proc=bandfft_kpt_get_ikpt()
 invovl => invovl_kpt(ikpt_this_proc)

 if(invovl%nprojs /= -1) then
   ! We have been here before, cleanup before remaking
   ABI_DEALLOCATE(invovl%gram_projs)
   ABI_DEALLOCATE(invovl%inv_sij)
   ABI_DEALLOCATE(invovl%inv_s_approx)
   invovl%nprojs = -1
 end if

 iaph3d = 1

 call timab(timer_mkinvovl,1,tsec)
 call timab(timer_mkinvovl_build_d,1,tsec)

 ! build nprojs
 invovl%nprojs = 0
 do itypat=1,ham%ntypat
   invovl%nprojs = invovl%nprojs + count(ham%indlmn(3,:,itypat)>0)*ham%nattyp(itypat)
 end do

 ABI_ALLOCATE(projs, (2, ham%npw_k, invovl%nprojs))
 ABI_ALLOCATE(invovl%inv_sij, (cplx, ham%lmnmax, ham%lmnmax, ham%ntypat))
 ABI_ALLOCATE(invovl%inv_s_approx, (cplx, ham%lmnmax, ham%lmnmax, ham%ntypat))
 ! workspace for inversion
 ABI_ALLOCATE(ipiv, (ham%lmnmax))
 ABI_ALLOCATE(work, (64*ham%lmnmax))

 projs = zero
 invovl%inv_sij = zero
 invovl%inv_s_approx = zero

 shift = 0
 do itypat = 1, ham%ntypat
   nlmn = count(ham%indlmn(3,:,itypat)>0)
   !! unpack ham%sij into inv_sij
   do jlmn = 1, nlmn
     invovl%inv_sij(1, jlmn, jlmn, itypat) = ham%sij(jlmn*(jlmn-1)/2 + jlmn, itypat)
     jlm=ham%indlmn(4,jlmn, itypat)
     do ilmn = 1, jlmn-1
       ilm=ham%indlmn(4,ilmn, itypat)
       if (ilm == jlm) then ! sparsity check
         invovl%inv_sij(1, ilmn, jlmn, itypat) = ham%sij(jlmn*(jlmn-1)/2 + ilmn, itypat)
         invovl%inv_sij(1, jlmn, ilmn, itypat) = ham%sij(jlmn*(jlmn-1)/2 + ilmn, itypat)
       end if
     end do
   end do

   ! Invert sij
   if(cplx == 2) then
     call ZHETRF('U', nlmn, invovl%inv_sij(:,:,:,itypat), ham%lmnmax, ipiv, work, (64*nlmn), info)
     call ZHETRI('U', nlmn, invovl%inv_sij(:,:,:,itypat), ham%lmnmax, ipiv, work, info)
   else
     call DSYTRF('U', nlmn, invovl%inv_sij(:,:,:,itypat), ham%lmnmax, ipiv, work, (64*nlmn), info)
     call DSYTRI('U', nlmn, invovl%inv_sij(:,:,:,itypat), ham%lmnmax, ipiv, work, info)
   end if
   ! complete the matrix
   do ilm=1, nlmn
     do jlm=1, ilm-1
       invovl%inv_sij(1,ilm,jlm,itypat) =  invovl%inv_sij(1,jlm,ilm,itypat)
     end do
   end do

   !! loop on atoms to build atom_projs and fill projs, s_projs
   do ia = 1, ham%nattyp(itypat)

     !! build atom_projs, from opernlb
     !! P = 4pi/sqrt(ucvol)* conj(diag(ph3d)) * ffnl * diag(parity), with parity = (-i)^l
     atom_projs(:,:,:) = 0

     ! start from 4pi/sqrt(ucvol)*ffnl
     ! atom_projs(1, :, 1:nlmn) = four_pi/sqrt(ham%ucvol) * ffnl(:, 1, 1:nlmn)
     ! TODO vectorize (DCOPY with stride)
     do ipw=1, ham%npw_k
       atom_projs(1,ipw, 1:nlmn) = four_pi/sqrt(ham%ucvol) * ffnl(ipw, 1, 1:nlmn, itypat)
     end do

     ! multiply by (-i)^l
     do ilmn=1,nlmn
       il=mod(ham%indlmn(1,ilmn, itypat),4);
       parity=(mod(il,2)==0)
       if (il>1) then
         ! multiply by -1
         atom_projs(:,:,ilmn) = -atom_projs(:,:,ilmn)
       end if
       if(.not. parity) then
         ! multiply by -i
         temp = atom_projs(2,:,ilmn)
         atom_projs(2,:,ilmn) = -atom_projs(1,:,ilmn)
         atom_projs(1,:,ilmn) =  temp
       end if
     end do

     ! multiply by conj(ph3d)
     do ilmn=1,nlmn
       temp = atom_projs(1, :, ilmn)
       atom_projs(1, :, ilmn) = atom_projs(1, :, ilmn) * ph3d(1, :, iaph3d) + atom_projs(2, :, ilmn) * ph3d(2, :, iaph3d)
       atom_projs(2, :, ilmn) = atom_projs(2, :, ilmn) * ph3d(1, :, iaph3d) - temp                   * ph3d(2, :, iaph3d)
     end do

     ! me_g0 trick
     if(ham%istwf_k == 2 .and. mpi_enreg%me_g0 == 1) then
       atom_projs(1,1,:) = atom_projs(1,1,:) / sqrt2
       atom_projs(2,1,:) = zero
     end if
     if(ham%istwf_k == 2) then
       atom_projs(:,:,:) = atom_projs(:,:,:)*sqrt2
     end if


     !! atom_projs and typat_s_projs are built, copy them to projs and inv_s_projs
     projs(:, :, shift+1:shift+nlmn) = atom_projs(:, :, 1:nlmn)
     shift = shift + nlmn

     ! build inv_s_approx = (D^-1+PtP)^-1 restricted to a single atom block
     ! can be optimized (real, build directly from ffnl)
     if(ia == 1) then
       ! D^-1
       invovl%inv_s_approx(1, :, :, itypat) = invovl%inv_sij(1, :, :, itypat)
       ! + PtP
       ABI_ALLOCATE(gram_proj, (cplx, nlmn, nlmn))
       call abi_xgemm(blas_transpose,'N', nlmn, nlmn, (3-cplx)*ham%npw_k, cone, atom_projs(:,:,1), (3-cplx)*ham%npw_k, &
&                     atom_projs(:,:,1), (3-cplx)*ham%npw_k, czero, gram_proj(:,:,1), nlmn,x_cplx=cplx)
       call xmpi_sum(gram_proj,mpi_enreg%comm_bandspinorfft,ierr)
       invovl%inv_s_approx(:,1:nlmn,1:nlmn,itypat) = invovl%inv_s_approx(:,1:nlmn,1:nlmn,itypat) + gram_proj(:,:,:)
       ABI_DEALLOCATE(gram_proj)
       ! ^-1
       if(cplx == 2) then
         call ZHETRF('U', nlmn, invovl%inv_s_approx(:,:,:,itypat), ham%lmnmax, ipiv, work, (64*nlmn), info)
         call ZHETRI('U', nlmn, invovl%inv_s_approx(:,:,:,itypat), ham%lmnmax, ipiv, work, info)
       else
         call DSYTRF('U', nlmn, invovl%inv_s_approx(:,:,:,itypat), ham%lmnmax, ipiv, work, (64*nlmn), info)
         call DSYTRI('U', nlmn, invovl%inv_s_approx(:,:,:,itypat), ham%lmnmax, ipiv, work, info)
       end if
       ! complete lower triangle of matrix
       do ilm=1, nlmn
         do jlm=1, ilm-1
           invovl%inv_s_approx(1, ilm, jlm, itypat) =  invovl%inv_s_approx(1, jlm, ilm, itypat)
           if(cplx == 2) then
             invovl%inv_s_approx(2, ilm, jlm, itypat) =  -invovl%inv_s_approx(2, jlm, ilm, itypat)
           end if
         end do
       end do
     end if

     iaph3d = iaph3d + 1
   end do
 end do
 ABI_DEALLOCATE(ipiv)
 ABI_DEALLOCATE(work)

 call timab(timer_mkinvovl_build_d, 2, tsec)
 call timab(timer_mkinvovl_build_ptp, 1, tsec)

 ! Compute P'P one column slice at a time (might be too big to fit in one proc)
 if(mpi_enreg%paral_kgb == 1) then
   ! Split the work evenly the fft processors
   array_nprojs_pp(:) = invovl%nprojs / mpi_enreg%nproc_fft
   ! not enough work, there's MOD(nprojs,mpi_enreg%nproc_fft) tasks left
   ! assign them to the first ones
   array_nprojs_pp(1:MOD(invovl%nprojs,mpi_enreg%nproc_fft)) = array_nprojs_pp(1:MOD(invovl%nprojs,mpi_enreg%nproc_fft)) + 1
 else
   array_nprojs_pp = invovl%nprojs
 end if
 ABI_ALLOCATE(invovl%gram_projs, (cplx,invovl%nprojs,array_nprojs_pp(mpi_enreg%me_fft+1)))
 shift = 0
 do iproc = 1, mpi_enreg%nproc_fft
   ! compute local contribution to slice iproc of gram_projs
   slice_size = array_nprojs_pp(iproc)
   ABI_ALLOCATE(gramwork, (cplx,invovl%nprojs,slice_size))
   call abi_xgemm(blas_transpose,'N', invovl%nprojs, slice_size, (3-cplx)*ham%npw_k, cone, projs(:,:,1), (3-cplx)*ham%npw_k, &
&                      projs(:, :, shift+1), (3-cplx)*ham%npw_k, czero, gramwork(:,:,1), invovl%nprojs,x_cplx=cplx)
   shift = shift + slice_size
   ! reduce on proc i
   call xmpi_sum_master(gramwork, iproc-1, mpi_enreg%comm_fft, ierr)
   if(iproc == mpi_enreg%me_fft+1) then
     invovl%gram_projs = gramwork
   end if
   ABI_DEALLOCATE(gramwork)
 end do
 call xmpi_sum(invovl%gram_projs,mpi_enreg%comm_band,ierr)

 call timab(timer_mkinvovl_build_ptp, 2, tsec)
 call timab(timer_mkinvovl,2,tsec)

 ABI_DEALLOCATE(projs)

 write(message,*) 'Invovl built'
 call wrtout(std_out,message,'COLL')

end subroutine make_invovl
!!***

!!****f* m_invovl/apply_invovl
!! NAME
!! apply_invovl
!!
!! FUNCTION
!! Applies the inverse of the overlap matrix to cwavef
!!
!! INPUTS
!!
!! PARENTS
!!      chebfi
!!
!! CHILDREN
!!      dsymm,zhemm
!!
!! SOURCE

 subroutine apply_invovl(ham, cwavef, sm1cwavef, cwaveprj, npw, ndat, mpi_enreg, nspinor)

 implicit none

 ! args
 type(gs_hamiltonian_type),intent(in) :: ham
 integer,intent(in) :: npw, ndat
 integer,intent(in) :: nspinor
 real(dp), intent(inout) :: cwavef(2, npw*nspinor*ndat) ! TODO should be in, fix nonlop
 type(mpi_type) :: mpi_enreg
 real(dp), intent(inout) :: sm1cwavef(2, npw*nspinor*ndat)
 type(pawcprj_type), intent(inout) :: cwaveprj(ham%natom,nspinor*ndat)

 real(dp),allocatable :: proj(:,:,:),sm1proj(:,:,:),PtPsm1proj(:,:,:)
 integer :: idat, iatom, nlmn, shift
 real(dp) :: tsec(2)

 integer :: choice, cpopt, paw_opt , cplx
 character :: blas_transpose
 type(pawcprj_type) :: cwaveprj_in(ham%natom,nspinor*ndat)

 integer :: ikpt_this_proc
 integer, parameter :: tim_nonlop = 13
 ! dummies
 real(dp) :: enlout(ndat), lambda_block(1), gvnlxc(1,1)
 integer, parameter :: nnlout = 0, idir = 0, signs = 2

 type(invovl_kpt_type), pointer :: invovl
 integer, parameter :: timer_apply_inv_ovl_opernla = 1630, timer_apply_inv_ovl_opernlb = 1631, &
&                      timer_apply_inv_ovl_inv_s = 1632

! *************************************************************************

 ikpt_this_proc=bandfft_kpt_get_ikpt()
 invovl => invovl_kpt(ikpt_this_proc)

 if(ham%istwf_k == 1) then
   cplx = 2
   blas_transpose = 'c'
 else
   cplx = 1
   blas_transpose = 't'
 end if

 ABI_ALLOCATE(proj, (cplx,invovl%nprojs,nspinor*ndat))
 ABI_ALLOCATE(sm1proj, (cplx,invovl%nprojs,nspinor*ndat))
 ABI_ALLOCATE(PtPsm1proj, (cplx,invovl%nprojs,nspinor*ndat))
 proj = zero
 sm1proj = zero
 PtPsm1proj = zero

 call timab(timer_apply_inv_ovl_opernla, 1, tsec)

 call pawcprj_alloc(cwaveprj_in,0,ham%dimcprj)

 ! get the cprj
 choice = 0 ! only compute cprj, nothing else
 cpopt = 0 ! compute and save cprj
 paw_opt = 3 ! S nonlocal operator
 if (mpi_enreg%paral_kgb==1) then
   call prep_nonlop(choice,cpopt,cwaveprj_in,enlout,ham,idir,lambda_block,ndat,mpi_enreg,&
&                   nnlout,paw_opt,signs,sm1cwavef,tim_nonlop,cwavef,gvnlxc,already_transposed=.true.)
 else
   call nonlop(choice,cpopt,cwaveprj_in,enlout,ham,idir,lambda_block,mpi_enreg,ndat,nnlout,&
&              paw_opt,signs,sm1cwavef,tim_nonlop,cwavef,gvnlxc)
 end if

 call timab(timer_apply_inv_ovl_opernla, 2, tsec)
 call timab(timer_apply_inv_ovl_inv_s, 1, tsec)

 ! copy cwaveprj_in to proj(:,:)
 do idat=1, ndat*nspinor
   shift = 0
   do iatom = 1, ham%natom
     nlmn = cwaveprj_in(iatom, idat)%nlmn
     proj(1:cplx, shift+1:shift+nlmn, idat) = cwaveprj_in(iatom, idat)%cp(1:cplx, 1:nlmn)
     shift = shift + nlmn
   end do
 end do

 !multiply by S^1
 call solve_inner(invovl, ham, cplx, mpi_enreg, proj, ndat*nspinor, sm1proj, PtPsm1proj)
 sm1proj = - sm1proj
 PtPsm1proj = - PtPsm1proj

 ! copy sm1proj to cwaveprj(:,:)
 do idat=1, ndat*nspinor
   shift = 0
   do iatom = 1, ham%natom
     nlmn = cwaveprj(iatom, idat)%nlmn
     cwaveprj(iatom, idat)%cp(1:cplx, 1:nlmn) = sm1proj(1:cplx, shift+1:shift+nlmn, idat)
     shift = shift + nlmn
   end do
 end do
 call timab(timer_apply_inv_ovl_inv_s, 2, tsec)
 call timab(timer_apply_inv_ovl_opernlb, 1, tsec)

 ! get the corresponding wf
 cpopt = 2 ! reuse cprj
 choice = 7 ! get wf from cprj, without the application of S
 paw_opt = 3
 if (mpi_enreg%paral_kgb==1) then
   call prep_nonlop(choice,cpopt,cwaveprj,enlout,ham,idir,lambda_block,ndat,mpi_enreg,nnlout,&
&                   paw_opt,signs,sm1cwavef,tim_nonlop,cwavef,gvnlxc,already_transposed=.true.)
 else
   call nonlop(choice,cpopt,cwaveprj,enlout,ham,idir,lambda_block,mpi_enreg,ndat,nnlout,paw_opt,&
&              signs,sm1cwavef,tim_nonlop,cwavef,gvnlxc)
 end if

 call timab(timer_apply_inv_ovl_opernlb, 2, tsec)

 ! copy PtPSm1proj to cwaveprj(:,:)
 do idat=1, ndat*nspinor
   shift = 0
   do iatom = 1, ham%natom
     nlmn = cwaveprj(iatom, idat)%nlmn
     cwaveprj(iatom, idat)%cp(1:cplx, 1:nlmn) = PtPsm1proj(1:cplx, shift+1:shift+nlmn, idat)
     shift = shift + nlmn
   end do
 end do

 sm1cwavef = cwavef + sm1cwavef
 call pawcprj_axpby(one, one, cwaveprj_in, cwaveprj)
 call pawcprj_free(cwaveprj_in)

 ABI_DEALLOCATE(proj)
 ABI_DEALLOCATE(sm1proj)
 ABI_DEALLOCATE(PtPsm1proj)

end subroutine apply_invovl
!!***

!!****f* m_invovl/solve_inner
!! NAME
!! solve_inner
!!
!! FUNCTION
!! Helper function: iteratively solves the inner system
!!
!! INPUTS
!!
!! PARENTS
!!      m_invovl
!!
!! CHILDREN
!!      dsymm,zhemm
!!
!! SOURCE
subroutine solve_inner(invovl, ham, cplx, mpi_enreg, proj, ndat, sm1proj, PtPsm1proj)

 use m_abi_linalg
 implicit none

 integer,intent(in) :: ndat,cplx
 type(invovl_kpt_type), intent(in) :: invovl
 real(dp), intent(inout) :: proj(cplx, invovl%nprojs,ndat), sm1proj(cplx, invovl%nprojs,ndat), PtPsm1proj(cplx, invovl%nprojs,ndat)
 real(dp), allocatable :: temp_proj(:,:,:)
 type(mpi_type), intent(in) :: mpi_enreg
 type(gs_hamiltonian_type),intent(in) :: ham

 integer :: array_nlmntot_pp(mpi_enreg%nproc_fft)
 integer :: nlmntot_this_proc, ibeg, iend, ierr, i, nprojs
 real(dp) :: resid(cplx, invovl%nprojs,ndat), precondresid(cplx, invovl%nprojs,ndat)
 real(dp) :: normprojs(ndat), errs(ndat), maxerr, previous_maxerr
 character(len=500) :: message

 real(dp), parameter :: precision = 1e-16 ! maximum relative error. TODO: use tolwfr ?
 real(dp) :: convergence_rate
 integer :: additional_steps_to_take

! *************************************************************************

 nprojs = invovl%nprojs
 normprojs = SUM(SUM(proj**2, 1),1)
 call xmpi_sum(normprojs, mpi_enreg%comm_fft, ierr)

 ! Compute work distribution : split nprojs evenly between the fft processors
 if(mpi_enreg%paral_kgb == 1) then
   array_nlmntot_pp(:) = nprojs / mpi_enreg%nproc_fft
   ! not enough work, there's MOD(nprojs,mpi_enreg%nproc_fft) tasks left
   ! assign them to the first ones
   array_nlmntot_pp(1:MOD(nprojs,mpi_enreg%nproc_fft)) = array_nlmntot_pp(1:MOD(nprojs,mpi_enreg%nproc_fft)) + 1
   ibeg = SUM(array_nlmntot_pp(1:mpi_enreg%me_fft)) + 1
   iend = ibeg + array_nlmntot_pp(1+mpi_enreg%me_fft) - 1
   nlmntot_this_proc = iend - ibeg + 1
 else
   ibeg = 1
   iend = nprojs
   nlmntot_this_proc = nprojs
 end if

 ABI_ALLOCATE(temp_proj, (cplx, nlmntot_this_proc,ndat))

 ! first guess for sm1proj
 call apply_block(ham, cplx, invovl%inv_s_approx, nprojs, ndat, proj, sm1proj)

 ! Iterative refinement
 ! TODO use a more efficient iterative algorithm than iterative refinement, use locking
 additional_steps_to_take = -1
 do i=1, 30
  ! compute resid = proj - (D^-1 + PtP)sm1proj
   call apply_block(ham, cplx, invovl%inv_sij, nprojs, ndat, sm1proj, resid)
   temp_proj = sm1proj(:,ibeg:iend,:)
   call abi_xgemm('N', 'N', nprojs, ndat, nlmntot_this_proc, cone, invovl%gram_projs(:,:,1), nprojs, &
&            temp_proj(:,:,1), nlmntot_this_proc, czero, PtPsm1proj(:,:,1), nprojs, x_cplx=cplx)
   call xmpi_sum(PtPsm1proj, mpi_enreg%comm_fft, ierr)
   resid = proj - resid - Ptpsm1proj
   ! exit check
   errs = SUM(SUM(resid**2, 1),1)
   call xmpi_sum(errs, mpi_enreg%comm_fft, ierr)

   maxerr = sqrt(MAXVAL(errs/normprojs))
   if(maxerr < precision .or. additional_steps_to_take == 1) then
     exit
     ! We might stall and never get to the specified precision because of machine errors.
     ! If we got to 1e-10, extrapolate convergence rate and determine the number of additional
     ! steps to take to reach precision
   else if(maxerr < 1e-10 .and. additional_steps_to_take == -1) then
     convergence_rate = -LOG(1e-10) / i
     additional_steps_to_take = CEILING(-LOG(precision/1e-10)/convergence_rate) + 1
   else if(additional_steps_to_take > 0) then
     if(previous_maxerr<maxerr)exit
     additional_steps_to_take = additional_steps_to_take - 1
   end if
   previous_maxerr=maxerr

   ! add preconditionned residual
   call apply_block(ham, cplx, invovl%inv_s_approx, nprojs, ndat, resid, precondresid)
   sm1proj = sm1proj + precondresid
 end do

 if(maxerr >= precision .and. maxerr >= 1e-10) then
   write(message, *) 'In invovl, max error was', maxerr, ' after 30 iterations'
   MSG_WARNING(message)
 else
   ! write(message,'(a,i2,a,es13.5)') 'Iterative solver in invovl finished in ', i, ' iterations, error', maxerr
   ! call wrtout(std_out,message,'COLL')
 end if

 ABI_DEALLOCATE(temp_proj)

end subroutine solve_inner
!!***

!!****f* m_invovl/apply_block
!! NAME
!! apply_block
!!
!! FUNCTION
!! Helper function: applies a block-diagonal matrix mat(lmnmax, lmnmax, ntypat)
!!
!! INPUTS
!!
!! PARENTS
!!      m_invovl
!!
!! CHILDREN
!!      dsymm,zhemm
!!
!! SOURCE
subroutine apply_block(ham, cplx, mat, nprojs, ndat, x, y)

  use m_abi_linalg
  implicit none

  integer,intent(in) :: ndat, nprojs, cplx
  real(dp), intent(inout) :: x(cplx, nprojs, ndat), y(cplx, nprojs, ndat)
  type(gs_hamiltonian_type),intent(in) :: ham
  real(dp), intent(in) :: mat(cplx, ham%lmnmax, ham%lmnmax, ham%ntypat)

  integer :: nlmn, shift, itypat, idat

! *************************************************************************

  do idat = 1, ndat
    shift = 1
    do itypat=1, ham%ntypat
      nlmn = count(ham%indlmn(3,:,itypat)>0)
      !! apply mat to all atoms at once
      ! perform natom multiplications of size nlmn
      if(cplx == 2) then
        call ZHEMM('L','U', nlmn, ham%nattyp(itypat), cone, mat(:, :, :, itypat), ham%lmnmax, &
&                x(:, shift:shift+nlmn*ham%nattyp(itypat)-1, idat), nlmn, czero, &
&                y(:,shift:shift+nlmn*ham%nattyp(itypat)-1,idat), nlmn)
      else
        call DSYMM('L','U', nlmn, ham%nattyp(itypat), one, mat(:, :, :, itypat), ham%lmnmax, &
&                x(:, shift:shift+nlmn*ham%nattyp(itypat)-1, idat), nlmn, zero, &
&                y(:,shift:shift+nlmn*ham%nattyp(itypat)-1,idat), nlmn)
      end if
      shift = shift + ham%nattyp(itypat)*nlmn
    end do
  end do

end subroutine apply_block
!!***

end MODULE m_invovl
!!***
