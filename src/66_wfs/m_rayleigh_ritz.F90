!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_rayleigh_ritz
!! NAME
!!  m_rayleigh_ritz
!!
!! FUNCTION
!! This file contains routines that perform the Rayleigh-Ritz,
!! either by forming the full matrix (_subdiago) or by forming the
!! distributed matrix in block-cyclic form (_distributed)
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2019 ABINIT group (AL)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_rayleigh_ritz

 use defs_basis
 use m_errors
 use m_cgtools
 use m_xmpi
 use m_abicore
 use m_abi_linalg
 use m_slk

 use defs_abitypes,   only : mpi_type
 use m_time,          only : timab
 use m_numeric_tools, only : pack_matrix

 implicit none

 private
!!***

 public :: rayleigh_ritz_subdiago
#if defined HAVE_LINALG_SCALAPACK
 public :: rayleigh_ritz_distributed
#endif
!!***

contains
!!***

!!****f* ABINIT/rayleigh_ritz_subdiago
!! NAME
!! rayleigh_ritz_subdiago
!!
!! FUNCTION
!! Performs a rayleigh-ritz procedure (subspace rotation), building the
!! hamiltonian/overlap matrices in full and calling the subdiago method
!!
!! INPUTS
!!  mpi_enreg=informations about MPI parallelization
!!  nband=number of bands at this k point for that spin polarization
!!  npw=number of plane waves at this k point
!!  nspinor=number of plane waves at this k point
!!  usepaw=if 1 we use the PAW method
!!
!! OUTPUT
!!  eig(nband)=array for holding eigenvalues (hartree)
!!
!! SIDE EFFECTS
!!  cg(2,*)=updated wavefunctions
!!  ghc(2,*)=updated ghc
!!  gsc(2,*)=updated gsc
!!  gvnlxc(2,*)=updated gvnlxc
!!
!! PARENTS
!!      chebfi
!!
!! CHILDREN
!!
!! NOTES
!!  TODO choose generalized eigenproblem or ortho + diago (see #if 1)
!!
!! SOURCE

subroutine rayleigh_ritz_subdiago(cg,ghc,gsc,gvnlxc,eig,has_fock,istwf_k,mpi_enreg,nband,npw,nspinor,usepaw)

 ! Arguments
 type(mpi_type),intent(inout) :: mpi_enreg
 integer,intent(in) :: nband,npw,nspinor,usepaw,istwf_k
 real(dp),intent(inout) :: cg(2,npw*nspinor*nband),gsc(2,npw*nspinor*nband),ghc(2,npw*nspinor*nband),gvnlxc(2,npw*nspinor*nband)
 real(dp),intent(out) :: eig(nband)
 logical :: has_fock

 ! Locals
 real(dp), allocatable :: subham(:), totham(:,:)
 real(dp), allocatable :: subovl(:), totovl(:,:)
 real(dp), allocatable :: evec(:,:), edummy(:,:)
 integer :: ierr, cplx, vectsize
 real(dp) :: gtempc(2,npw*nspinor*nband), tsec(2)
 character :: blas_transpose

 integer, parameter :: timer_chebfi = 1600, timer_alltoall = 1601, timer_apply_inv_ovl = 1602, timer_rotation = 1603
 integer, parameter :: timer_subdiago = 1604, timer_subham = 1605, timer_ortho = 1606, timer_getghc = 1607
 integer, parameter :: timer_residuals = 1608, timer_update_eigen = 1609, timer_sync = 1610

 ! *************************************************************************

 if(istwf_k == 1) then
   cplx = 2
   vectsize = npw*nspinor
   blas_transpose = 'c'
 else
   cplx = 1
   vectsize = 2*npw*nspinor
   blas_transpose = 't'
 end if

#if 1
 call timab(timer_subham, 1, tsec)

 ! Transform cg, ghc and maybe gsc, according to istwf_k
 if(istwf_k == 2) then
   cg = cg * sqrt2
   if(mpi_enreg%me_g0 == 1) cg(:, 1:npw*nspinor*nband:npw) = cg(:, 1:npw*nspinor*nband:npw) / sqrt2
   ghc = ghc * sqrt2
   if(mpi_enreg%me_g0 == 1) ghc(:, 1:npw*nspinor*nband:npw) = ghc(:, 1:npw*nspinor*nband:npw) / sqrt2
   if(usepaw == 1) then
     gsc = gsc * sqrt2
     if(mpi_enreg%me_g0 == 1) gsc(:, 1:npw*nspinor*nband:npw) = gsc(:, 1:npw*nspinor*nband:npw) / sqrt2
   end if
 end if

 ! Build, pack and sum suham
 ABI_ALLOCATE(subham, (cplx*nband*(nband+1)/2))
 ABI_ALLOCATE(totham, (cplx, nband*nband))
 call abi_xgemm(blas_transpose,'n',nband,nband,vectsize,cone,ghc,vectsize,&
& cg,vectsize,czero,totham,nband, x_cplx=cplx)
 call pack_matrix(totham, subham, nband, cplx)
 ABI_DEALLOCATE(totham)
 call xmpi_sum(subham,mpi_enreg%comm_bandspinorfft,ierr)


 ! Same for subovl
 ABI_ALLOCATE(subovl, (cplx*nband*(nband+1)/2))
 ABI_ALLOCATE(totovl, (cplx, nband*nband))
 if(usepaw == 1) then
   call abi_xgemm(blas_transpose,'n',nband,nband,vectsize,cone,gsc,vectsize,&
&   cg,vectsize,czero,totovl,nband, x_cplx=cplx)
 else
   call abi_xgemm(blas_transpose,'n',nband,nband,vectsize,cone,cg,vectsize,&
&   cg,vectsize,czero,totovl,nband, x_cplx=cplx)
 end if
 call pack_matrix(totovl, subovl, nband, cplx)
 ABI_DEALLOCATE(totovl)
 call xmpi_sum(subovl,mpi_enreg%comm_bandspinorfft,ierr)


 ! Transform back
 if(istwf_k == 2) then
   cg = cg / sqrt2
   if(mpi_enreg%me_g0 == 1) cg(:, 1:npw*nspinor*nband:npw) = cg(:, 1:npw*nspinor*nband:npw) * sqrt2
   ghc = ghc / sqrt2
   if(mpi_enreg%me_g0 == 1) ghc(:, 1:npw*nspinor*nband:npw) = ghc(:, 1:npw*nspinor*nband:npw) * sqrt2
   if(usepaw == 1) then
     gsc = gsc / sqrt2
     if(mpi_enreg%me_g0 == 1) gsc(:, 1:npw*nspinor*nband:npw) = gsc(:, 1:npw*nspinor*nband:npw) * sqrt2
   end if
 end if
 call timab(timer_subham, 2, tsec)


 call timab(timer_subdiago, 1, tsec)
 ABI_ALLOCATE(evec, (cplx*nband, nband))

 call abi_xhpgv(1,'V','U',nband,subham,subovl,eig,evec,nband,istwf_k=istwf_k,use_slk=mpi_enreg%paral_kgb)

 ABI_DEALLOCATE(subham)
 ABI_DEALLOCATE(subovl)

! Fix the phase (this is because of the simultaneous diagonalisation of this
! matrix by different processors, allowing to get different unitary transforms, thus breaking the
! coherency of parts of cg stored on different processors).
! call cg_normev(evec,nband,nband)  ! Unfortunately, for cg_normev to work, one needs the vectors to be normalized, so uses fxphas_seq
 ABI_ALLOCATE(edummy, (cplx*nband, nband))
 call fxphas_seq(evec,edummy,0,0,1,nband*nband,nband*nband,nband,nband,0)
 ABI_DEALLOCATE(edummy)



 ! Rotate
 call abi_xgemm('n','n',vectsize,nband, nband,cone,cg , vectsize, evec, nband, czero, gtempc, vectsize, x_cplx=cplx)
 cg = gtempc
 call abi_xgemm('n','n',vectsize,nband, nband,cone,ghc, vectsize, evec, nband, czero, gtempc, vectsize, x_cplx=cplx)
 ghc = gtempc
 if(usepaw == 1) then
   call abi_xgemm('n','n',vectsize,nband, nband,cone,gsc, vectsize, evec, nband, czero, gtempc, vectsize, x_cplx=cplx)
   gsc = gtempc
 endif
 if(usepaw==0 .or. has_fock)then
   call abi_xgemm('n','n',vectsize,nband, nband,cone,gvnlxc, vectsize, evec, nband, czero, gtempc, vectsize, x_cplx=cplx)
   gvnlxc = gtempc
 end if
 ABI_DEALLOCATE(evec)
 call timab(timer_subdiago, 2, tsec)

#else
 !! TODO non-functional, should be rewritten. Possibly faster (tests needed)
 write(message, *) 'Transposed, orthogonalizing'
 call wrtout(std_out,message,'COLL')

 ! orthonormalization
 call timab(timer_ortho, 1, tsec)
 if (usepaw==1) then
   call abi_xorthonormalize(cg, gsc,nband, mpi_enreg%comm_bandspinorfft, sqgram, npw*nspinor, 2)
 else
   call abi_xorthonormalize(cg, cg, nband, mpi_enreg%comm_bandspinorfft, sqgram, npw*nspinor, 2)
 end if
 call timab(timer_ortho, 2, tsec)

 ! rotate ghc, gsc and gvnlxc
 call timab(timer_rotation, 1, tsec)
 call abi_xtrsm('r','u','n','n',npw*nspinor,nband,cone,sqgram,nband, ghc,npw*nspinor,x_cplx=2)
 if(usepaw==1) then
   call abi_xtrsm('r','u','n','n',npw*nspinor,nband,cone,sqgram,nband, gsc,npw*nspinor,x_cplx=2)
 endif
 if(usepaw==0 .or has_fock)then
   call abi_xtrsm('r','u','n','n',npw*nspinor,nband,cone,sqgram,nband, gvnlxc,npw*nspinor,x_cplx=2)
 end if
 call timab(timer_rotation, 2, tsec)

 write(message, *) 'Orthogonalized, building subham'
 call wrtout(std_out,message,'COLL')

 ! build hamiltonian  in subspace
 call timab(timer_subham, 1, tsec)
 call abi_xgemm(blas_transpose,'n',nband,nband,npw*nspinor,cone,ghc,npw*nspinor,&
& cg,npw*nspinor,czero,totham,nband, x_cplx=2)
 ! pack for compatibility with subdiago
 call pack_matrix(totham, subham, nband)
 call xmpi_sum(subham,mpi_enreg%comm_bandspinorfft,ierr)
 call timab(timer_subham, 2, tsec)

 write(message, *) 'Subham built, diagonalizing'
 call wrtout(std_out,message,'COLL')

 ! Rayleigh-Ritz
 call timab(timer_subdiago,1,tsec)
 call subdiago(cg,eig,evec,gsc,0,0,gs_hamk%istwf_k,&
& mcg,mcg,nband,npw,nspinor,dtset%paral_kgb,&
& subham,dummy,0,gs_hamk%usepaw,mpi_enreg%me_g0)
 call timab(timer_subdiago,2,tsec)

 write(message, *) 'Diagonalization done'
 call wrtout(std_out,message,'COLL')

 ! Rotate ghc and gvnlxc according to evecs
 call timab(timer_rotation, 1, tsec)
 call abi_xgemm('n','n',npw*nspinor,nband, nband,cone,ghc, npw*nspinor, evec, nband, czero, gtempc, npw*nspinor, x_cplx=2)
 ghc = gtempc
 if(usepaw==0 .or has_fock)then
   call abi_xgemm('n','n',npw*nspinor,nband, nband,cone,gvnlxc, npw*nspinor, evec, nband, czero, gtempc, npw*nspinor, x_cplx=2)
   gvnlxc = gtempc
 end if
 call timab(timer_rotation, 2, tsec)

#endif

end subroutine rayleigh_ritz_subdiago
!!***

#if defined HAVE_LINALG_SCALAPACK
!!****f* ABINIT/rayleigh_ritz_distributed
!! NAME
!! rayleigh_ritz_distributed
!!
!! FUNCTION
!! Performs a rayleigh-ritz procedure (subspace rotation), building the distributed
!! hamiltonian/overlap matrices directly, and calling the ScaLapack routines
!!
!! INPUTS
!!  mpi_enreg=informations about MPI parallelization
!!  nband=number of bands at this k point for that spin polarization
!!  npw=number of plane waves at this k point
!!  nspinor=number of plane waves at this k point
!!  usepaw=do we use the PAW method
!!
!! OUTPUT
!!  eig(nband)=array for holding eigenvalues (hartree)
!!
!! SIDE EFFECTS
!!  cg(2,*)=updated wavefunctions
!!  ghc(2,*)=updated ghc
!!  gsc(2,*)=updated gsc
!!  gvnlxc(2,*)=updated gvnlxc
!!
!! PARENTS
!!      chebfi
!!
!! CHILDREN
!!
!! NOTES
!!  Should profile for large test cases and see where the bottleneck is.
!!  Is it the copies? Should we do partial GEMMs?
!!  Is it the latency? Should we buffer more?
!!  Should we overlap computations and communications? (easy in theory, tedious in practice)
!!
!! SOURCE

subroutine rayleigh_ritz_distributed(cg,ghc,gsc,gvnlxc,eig,has_fock,istwf_k,mpi_enreg,nband,npw,nspinor,usepaw)

 integer,external :: NUMROC

 ! Arguments
 type(mpi_type),intent(inout) :: mpi_enreg
 integer,intent(in) :: nband,npw,nspinor,usepaw,istwf_k
 real(dp),intent(inout) :: cg(2,npw*nspinor*nband),gsc(2,npw*nspinor*nband),ghc(2,npw*nspinor*nband),gvnlxc(2,npw*nspinor*nband)
 real(dp),intent(out) :: eig(nband)
 logical :: has_fock

 ! Locals
 integer :: blocksize,nbproc,iproc,ierr,cplx,vectsize
 integer :: buffsize_iproc(2), coords_iproc(2), grid_dims(2)
 real(dp) :: cg_new(2,npw*nspinor*nband),gsc_or_vnlxc_new(2,npw*nspinor*nband),ghc_new(2,npw*nspinor*nband)
 real(dp), allocatable :: ham_iproc(:,:), ovl_iproc(:,:), evec_iproc(:,:), left_temp(:,:), right_temp(:,:)
 real(dp) :: tsec(2)
 type(matrix_scalapack) :: sca_ham, sca_ovl, sca_evec
 character(len=500) :: message
 character :: blas_transpose

 integer, parameter :: timer_chebfi = 1600, timer_alltoall = 1601, timer_apply_inv_ovl = 1602, timer_rotation = 1603
 integer, parameter :: timer_subdiago = 1604, timer_subham = 1605, timer_ortho = 1606, timer_getghc = 1607
 integer, parameter :: timer_residuals = 1608, timer_update_eigen = 1609, timer_sync = 1610

 ! *************************************************************************

 if(istwf_k == 1) then
   cplx = 2
   vectsize = npw*nspinor
   blas_transpose = 'c'
 else
   cplx = 1
   vectsize = 2*npw*nspinor
   blas_transpose = 't'
 end if

 !write(message, *) 'RR: init'
 !call wrtout(std_out,message,'COLL')
 !======================================================================================================
 ! Init Scalapack matrices
 !======================================================================================================
 call init_matrix_scalapack(sca_ham ,nband,nband,slk_processor,istwf_k,10)
 call init_matrix_scalapack(sca_ovl ,nband,nband,slk_processor,istwf_k,10)
 call init_matrix_scalapack(sca_evec,nband,nband,slk_processor,istwf_k,10)

 ! Get info
 blocksize = sca_ham%sizeb_blocs(1) ! Assume square blocs
 nbproc = slk_processor%grid%nbprocs
 grid_dims = slk_processor%grid%dims

 !======================================================================================================
 ! Build hamiltonian and overlap matrices
 !======================================================================================================
 ! TODO maybe we should avoid copies at the price of less BLAS efficiency (when blocksize is small)? must profile.

 call timab(timer_subham, 1, tsec)

 ! Transform cg, ghc and maybe gsc, according to istwf_k
 if(istwf_k == 2) then
   cg = cg * sqrt2
   if(mpi_enreg%me_g0 == 1) cg(:, 1:npw*nspinor*nband:npw) = cg(:, 1:npw*nspinor*nband:npw) / sqrt2
   ghc = ghc * sqrt2
   if(mpi_enreg%me_g0 == 1) ghc(:, 1:npw*nspinor*nband:npw) = ghc(:, 1:npw*nspinor*nband:npw) / sqrt2
   if(usepaw == 1) then
     gsc = gsc * sqrt2
     if(mpi_enreg%me_g0 == 1) gsc(:, 1:npw*nspinor*nband:npw) = gsc(:, 1:npw*nspinor*nband:npw) / sqrt2
   end if
 end if

 do iproc=0,nbproc-1
   ! Build the local matrix belonging to processor iproc
   !write(message, *) 'RR: build', iproc
   !call wrtout(std_out,message,'COLL')

   ! Get coordinates of iproc
   coords_iproc(1) = INT(iproc / grid_dims(2))
   coords_iproc(2) = MOD(iproc,  grid_dims(2))

   ! Get buffersize of iproc
   buffsize_iproc(1) = NUMROC(nband,blocksize,coords_iproc(1),0,slk_processor%grid%dims(1))
   buffsize_iproc(2) = NUMROC(nband,blocksize,coords_iproc(2),0,slk_processor%grid%dims(2))

   ! Allocate matrices_iproc, that will gather the contribution of this proc to the block owned by iproc
   ABI_ALLOCATE(ham_iproc, (cplx*buffsize_iproc(1), buffsize_iproc(2)))
   ABI_ALLOCATE(ovl_iproc, (cplx*buffsize_iproc(1), buffsize_iproc(2)))

   ! Build them
   ABI_ALLOCATE(left_temp,  (2, npw*nspinor*buffsize_iproc(1)))
   ABI_ALLOCATE(right_temp, (2, npw*nspinor*buffsize_iproc(2)))

   ! ovl
   call from_mat_to_block_cyclic(cg, npw*nspinor, nband, left_temp, &
&   buffsize_iproc(1), blocksize, coords_iproc(1), grid_dims(1))
   if(usepaw == 1) then
     call from_mat_to_block_cyclic(gsc, npw*nspinor, nband, right_temp, &
&     buffsize_iproc(2), blocksize, coords_iproc(2), grid_dims(2))
   else
     call from_mat_to_block_cyclic(cg, npw*nspinor, nband, right_temp, &
&     buffsize_iproc(2), blocksize, coords_iproc(2), grid_dims(2))
   end if
   call abi_xgemm(blas_transpose,'n',buffsize_iproc(1),buffsize_iproc(2),vectsize,cone,left_temp,vectsize,&
&   right_temp,vectsize,czero,ovl_iproc,buffsize_iproc(1), x_cplx=cplx)

   ! ham
   call from_mat_to_block_cyclic(ghc, npw*nspinor, nband, right_temp, &
&   buffsize_iproc(2), blocksize, coords_iproc(2), grid_dims(2))
   call abi_xgemm(blas_transpose,'n',buffsize_iproc(1),buffsize_iproc(2),vectsize,cone,left_temp,vectsize,&
&   right_temp,vectsize,czero,ham_iproc,buffsize_iproc(1), x_cplx=cplx)

   ! Sum to iproc, and fill sca_ matrices
   call xmpi_sum_master(ham_iproc, iproc, slk_communicator, ierr)
   call xmpi_sum_master(ovl_iproc, iproc, slk_communicator, ierr)
   if(iproc == slk_processor%myproc) then
     ! DCOPY to bypass the real/complex issue
     if(cplx == 2) then
       call DCOPY(cplx*buffsize_iproc(1)*buffsize_iproc(2), ham_iproc, 1, sca_ham%buffer_cplx, 1)
       call DCOPY(cplx*buffsize_iproc(1)*buffsize_iproc(2), ovl_iproc, 1, sca_ovl%buffer_cplx, 1)
     else
       call DCOPY(cplx*buffsize_iproc(1)*buffsize_iproc(2), ham_iproc, 1, sca_ham%buffer_real, 1)
       call DCOPY(cplx*buffsize_iproc(1)*buffsize_iproc(2), ovl_iproc, 1, sca_ovl%buffer_real, 1)
     end if
   end if

   ABI_DEALLOCATE(ham_iproc)
   ABI_DEALLOCATE(ovl_iproc)
   ABI_DEALLOCATE(left_temp)
   ABI_DEALLOCATE(right_temp)
 end do

 ! Final sum
 if(cplx == 2) then
   call xmpi_sum(sca_ham%buffer_cplx, slk_complement_communicator, ierr)
   call xmpi_sum(sca_ovl%buffer_cplx, slk_complement_communicator, ierr)
 else
   call xmpi_sum(sca_ham%buffer_real, slk_complement_communicator, ierr)
   call xmpi_sum(sca_ovl%buffer_real, slk_complement_communicator, ierr)
 end if

 ! Transform back
 if(istwf_k == 2) then
   cg = cg / sqrt2
   if(mpi_enreg%me_g0 == 1) cg(:, 1:npw*nspinor*nband:npw) = cg(:, 1:npw*nspinor*nband:npw) * sqrt2
   ghc = ghc / sqrt2
   if(mpi_enreg%me_g0 == 1) ghc(:, 1:npw*nspinor*nband:npw) = ghc(:, 1:npw*nspinor*nband:npw) * sqrt2
   if(usepaw == 1) then
     gsc = gsc / sqrt2
     if(mpi_enreg%me_g0 == 1) gsc(:, 1:npw*nspinor*nband:npw) = gsc(:, 1:npw*nspinor*nband:npw) * sqrt2
   end if
 end if
 call timab(timer_subham, 2, tsec)

 !======================================================================================================
 ! Do the diagonalization
 !======================================================================================================
 !write(message, *) 'RR: diag'
 !call wrtout(std_out,message,'COLL')
 call timab(timer_subdiago, 1, tsec)
 call compute_generalized_eigen_problem(slk_processor,sca_ham,sca_ovl,&
& sca_evec,eig,slk_communicator,istwf_k)
 call timab(timer_subdiago, 2, tsec)

 !======================================================================================================
 ! Perform rotation
 !======================================================================================================
 call timab(timer_rotation, 1, tsec)
 cg_new = zero
 gsc_or_vnlxc_new = zero
 ghc_new = zero
 do iproc=0,nbproc-1
   ! Compute the contribution to the rotated matrices from this block
   !write(message, *) 'RR: rot', iproc
   !call wrtout(std_out,message,'COLL')

   ! Get coordinates of iproc
   coords_iproc(1) = INT(iproc / grid_dims(2))
   coords_iproc(2) = MOD(iproc,  grid_dims(2))

   ! Get buffersize of iproc
   buffsize_iproc(1) = NUMROC(nband,blocksize,coords_iproc(1),0,slk_processor%grid%dims(1))
   buffsize_iproc(2) = NUMROC(nband,blocksize,coords_iproc(2),0,slk_processor%grid%dims(2))

   ! Get data from iproc
   ABI_ALLOCATE(evec_iproc, (cplx*buffsize_iproc(1), buffsize_iproc(2)))
   if(iproc == slk_processor%myproc) then
     if(cplx == 2) then
       call DCOPY(cplx*buffsize_iproc(1)*buffsize_iproc(2), sca_evec%buffer_cplx, 1, evec_iproc, 1)
     else
       call DCOPY(cplx*buffsize_iproc(1)*buffsize_iproc(2), sca_evec%buffer_real, 1, evec_iproc, 1)
     end if
   end if
   call xmpi_bcast(evec_iproc,iproc,slk_communicator,ierr)

   ! Compute contribution to the rotated matrices from iproc
   ABI_ALLOCATE(left_temp,  (2,npw*nspinor*buffsize_iproc(1)))
   ABI_ALLOCATE(right_temp,  (2,npw*nspinor*buffsize_iproc(2)))

   ! cg
   call from_mat_to_block_cyclic(cg, npw*nspinor, nband, left_temp, &
&   buffsize_iproc(1), blocksize, coords_iproc(1), grid_dims(1))
   call abi_xgemm('n','n',vectsize,buffsize_iproc(2),buffsize_iproc(1),cone,left_temp,vectsize,&
&   evec_iproc, buffsize_iproc(1), czero, right_temp, vectsize, x_cplx=cplx)
   call from_block_cyclic_to_mat(cg_new, npw*nspinor, nband, right_temp, &
&   buffsize_iproc(2), blocksize, coords_iproc(2), grid_dims(2))

   ! ghc
   call from_mat_to_block_cyclic(ghc, npw*nspinor, nband, left_temp, &
&   buffsize_iproc(1), blocksize, coords_iproc(1), grid_dims(1))
   call abi_xgemm('n','n',vectsize,buffsize_iproc(2),buffsize_iproc(1),cone,left_temp,vectsize,&
&   evec_iproc, buffsize_iproc(1), czero, right_temp, vectsize,x_cplx=cplx)
   call from_block_cyclic_to_mat(ghc_new, npw*nspinor, nband, right_temp,&
&   buffsize_iproc(2), blocksize, coords_iproc(2), grid_dims(2))

   ! gsc or vnlc
   if(usepaw == 1) then
     call from_mat_to_block_cyclic(gsc, npw*nspinor, nband, left_temp, &
&     buffsize_iproc(1), blocksize, coords_iproc(1), grid_dims(1))
   endif
   if(usepaw==0 .or. has_fock)then
     call from_mat_to_block_cyclic(gvnlxc, npw*nspinor, nband, left_temp, &
&     buffsize_iproc(1), blocksize, coords_iproc(1), grid_dims(1))
   end if
   call abi_xgemm('n','n',vectsize,buffsize_iproc(2),buffsize_iproc(1),cone,left_temp,vectsize,&
&   evec_iproc, buffsize_iproc(1), czero, right_temp, vectsize, x_cplx=cplx)
   call from_block_cyclic_to_mat(gsc_or_vnlxc_new, npw*nspinor, nband, right_temp, &
&   buffsize_iproc(2), blocksize, coords_iproc(2), grid_dims(2))

   ABI_DEALLOCATE(evec_iproc)
   ABI_DEALLOCATE(left_temp)
   ABI_DEALLOCATE(right_temp)
 end do

 ! Overwrite with _new
 cg = cg_new
 ghc = ghc_new
 if(usepaw == 1) then
   gsc = gsc_or_vnlxc_new
 endif
 if(usepaw==0 .or. has_fock)then
   gvnlxc = gsc_or_vnlxc_new
 end if
 call timab(timer_rotation, 2, tsec)

end subroutine rayleigh_ritz_distributed
!!***

!!****f* ABINIT/from_mat_to_block_cyclic
!! NAME
!! from_mat_to_block_cyclic
!!
!! FUNCTION
!! Fills block_cyclic_mat with the columns of full_mat owned by iproc, using a 1D block-cyclic distribution
!!
!! INPUTS
!! full_mat(2,vectsize*nband)=full input matrix
!! vectsize=number of rows
!! nband=number of columns of the full matrix
!! buffsize=size of the local part of the matrix
!! blocksize=size of block for the block-cyclic distribution
!! iproc=local processor number
!! nprocs=total number of processors
!!
!! OUTPUT
!! block_cyclic_mat(2,vectsize*buffsize)=local part of full_mat
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      rayleigh_ritz
!!
!! CHILDREN
!!
!! SOURCE
subroutine from_mat_to_block_cyclic(full_mat, vectsize, nband, block_cyclic_mat, buffsize, blocksize, iproc, nprocs)

 integer, intent(in) :: vectsize, nband, buffsize, blocksize, iproc, nprocs
 real(dp), intent(in) :: full_mat(2, vectsize*nband)
 real(dp), intent(inout) :: block_cyclic_mat(2, vectsize*buffsize)

 integer :: shift_fullmat, shift_block_cyclic_mat
 integer :: cur_bsize

 ! *************************************************************************

 shift_fullmat = iproc*blocksize
 shift_block_cyclic_mat = 0

 do while(.true.)
   cur_bsize = MIN(blocksize, nband-shift_fullmat)
   if(cur_bsize <= 0) exit
   block_cyclic_mat(:, shift_block_cyclic_mat*vectsize+1 : shift_block_cyclic_mat*vectsize + cur_bsize*vectsize) = &
&   full_mat(:, shift_fullmat*vectsize+1 : shift_fullmat*vectsize + cur_bsize*vectsize)

   shift_block_cyclic_mat = shift_block_cyclic_mat + blocksize
   shift_fullmat = shift_fullmat + blocksize*nprocs
 end do

end subroutine from_mat_to_block_cyclic
!!***

!!****f* ABINIT/from_block_cyclic_to_mat
!! NAME
!! from_block_cyclic_to_mat
!!
!! FUNCTION
!! Fills the columns of full_mat owned by iproc with block_cyclic_mat, using a 1D block-cyclic distribution
!!
!! INPUTS
!! block_cyclic_mat(2,vectsize*buffsize)=local part of full_mat
!! vectsize=number of rows
!! nband=number of columns of the full matrix
!! buffsize=size of the local part of the matrix
!! blocksize=size of block for the block-cyclic distribution
!! iproc=local processor number
!! nprocs=total number of processors
!!
!! OUTPUT
!! full_mat(2,vectsize*nband)=full matrix
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      rayleigh_ritz
!!
!! CHILDREN
!!
!! SOURCE

subroutine from_block_cyclic_to_mat(full_mat, vectsize, nband, block_cyclic_mat, buffsize, blocksize, iproc, nprocs)

 integer, intent(in) :: vectsize, nband, buffsize, blocksize, iproc, nprocs
 real(dp), intent(inout) :: full_mat(2, vectsize*nband)
 real(dp), intent(in) :: block_cyclic_mat(2, vectsize*buffsize)

 integer :: shift_fullmat, shift_block_cyclic_mat
 integer :: cur_bsize

 ! *************************************************************************

 shift_fullmat = iproc*blocksize
 shift_block_cyclic_mat = 0

 do while(.true.)
   cur_bsize = MIN(blocksize, nband-shift_fullmat)
   if(cur_bsize <= 0) exit
   full_mat(:, shift_fullmat*vectsize+1 : shift_fullmat*vectsize + cur_bsize*vectsize) =&
&   full_mat(:, shift_fullmat*vectsize+1 : shift_fullmat*vectsize + cur_bsize*vectsize) +&
&   block_cyclic_mat(:, shift_block_cyclic_mat*vectsize+1 : shift_block_cyclic_mat*vectsize + cur_bsize*vectsize)

   shift_block_cyclic_mat = shift_block_cyclic_mat + blocksize
   shift_fullmat = shift_fullmat + blocksize*nprocs
 end do

end subroutine from_block_cyclic_to_mat
!!***

#endif
!ScaLapack

end module m_rayleigh_ritz
!!***
