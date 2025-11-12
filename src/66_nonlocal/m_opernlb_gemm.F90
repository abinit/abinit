!!****m* ABINIT/m_opernlb_gemm
!! NAME
!!  m_opernlb_gemm
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2022 ABINIT group (MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_opernlb_gemm

 use defs_basis
 use m_abicore
 use m_errors
 USE_MPI
 use m_xmpi
 use m_abi_linalg
 use m_gemm_nonlop_projectors

 use defs_abitypes, only : MPI_type
 use m_time,        only : timab

#ifdef HAVE_FC_ISO_C_BINDING
 use, intrinsic :: iso_c_binding, only : c_ptr,c_loc,c_size_t
#endif

 implicit none

 private
!!***

 public :: opernlb_gemm

contains
!!***


!----------------------------------------------------------------------

!!****f* m_opernlb_gemm/opernlb_gemm_distributed
!! NAME
!! opernlb_gemm_distributed
!!
!! FUNCTION
!! Distributed version of "opernlb" GEMM called in gemm_nonlop.
!!
!! INPUTS
!!
!! SOURCE
subroutine opernlb_gemm_distributed(rank,nprocs,npw,ndat,&
&                                   nprojs,nprojs_blk,nprojs_last_blk,cplex,&
&                                   projs_local,projections,vectout,gpu_option)
 integer,  intent(in)     :: rank,nprocs,npw,ndat,gpu_option
 integer,  intent(in)     :: nprojs,nprojs_blk,nprojs_last_blk,cplex
 real(dp), intent(in),  target    :: projs_local(cplex,npw,nprojs_last_blk)
 real(dp), intent(in),  target    :: projections(cplex,nprojs,ndat)
 real(dp), intent(out), target    :: vectout(cplex,npw*ndat)

 !Local variables
 integer :: iblock,ibeg,iend,req(2),ierr,nprojs_cur_blk,rank_prev,rank_next
 complex(dp) :: beta
 real(dp), ABI_CONTIGUOUS pointer :: recv_buf(:,:,:), work_buf(:,:,:)
 real(dp), allocatable, target  :: projs_recv(:,:,:)
! *************************************************************************

 ABI_MALLOC(projs_recv, (cplex, npw, nprojs_last_blk))
#ifdef HAVE_OPENMP_OFFLOAD
 !$OMP TARGET ENTER DATA MAP(alloc:projs_recv) IF(gpu_option==ABI_GPU_OPENMP)
#endif

 rank_next=modulo(rank + 1,nprocs)
 rank_prev=rank - 1
 if(rank_prev == -1) rank_prev = nprocs - 1

 beta = czero

 do iblock=1,nprocs

   if(rank+iblock == nprocs) then
     nprojs_cur_blk = nprojs_last_blk
   else
     nprojs_cur_blk = nprojs_blk
   end if

   if(modulo(iblock,2)==1) then
!    XG20241028 : This coding confused the gnu_8.5 compiler of buda2_gnu_8.5_cuda, wrt the contiguous character of the target.
!    It declared an error. Make it simple !
!    work_buf => projs_local(1:cplex,1:npw,1:nprojs_last_blk)
!    recv_buf => projs_recv(1:cplex,1:npw,1:nprojs_last_blk)
     work_buf => projs_local
     recv_buf => projs_recv
   else
!    XG20241028 : Same as above
!    work_buf => projs_recv(1:cplex,1:npw,1:nprojs_last_blk)
!    recv_buf => projs_local(1:cplex,1:npw,1:nprojs_last_blk)
     work_buf => projs_recv
     recv_buf => projs_local
   end if

   if(gpu_option == ABI_GPU_DISABLED) then
     call xmpi_isend(work_buf,rank_prev,iblock,gemm_nonlop_block_comm,req(1),ierr)
     call xmpi_irecv(recv_buf,rank_next,iblock,gemm_nonlop_block_comm,req(2),ierr)
   else if(gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
#ifndef HAVE_GPU_MPI

     ! GPU-aware MPI not available : perform MPI comms on CPU
     !$OMP TARGET UPDATE FROM(work_buf) if(iblock==1)
     call xmpi_isend(work_buf,rank_prev,iblock,gemm_nonlop_block_comm,req(1),ierr)
     call xmpi_irecv(recv_buf,rank_next,iblock,gemm_nonlop_block_comm,req(2),ierr)

#else

     ! GPU-aware MPI available : pass GPU buffers to MPI
     call xmpi_isend(work_buf,rank_prev,iblock,gemm_nonlop_block_comm,req(1),ierr,use_omp_map=.true.)
     call xmpi_irecv(recv_buf,rank_next,iblock,gemm_nonlop_block_comm,req(2),ierr,use_omp_map=.true.)

#endif
#endif
   end if

   ibeg = 1 + modulo(rank+iblock-1,nprocs)*nprojs_blk
   iend = ibeg+nprojs_cur_blk-1

   if(gpu_option == ABI_GPU_DISABLED) then
     if(cplex==2) then
       call abi_zgemm_2r('N', 'N', npw, ndat, nprojs_cur_blk, cone, &
       &                 work_buf, npw, &
       &                 projections(:,ibeg:iend,:), nprojs_cur_blk,&
       &                 beta, &
       &                 vectout, npw)
     else
       call DGEMM('N', 'N', npw, ndat, nprojs_cur_blk, one, &
       &          work_buf, npw, &
       &          projections(:,ibeg:iend,:), nprojs_cur_blk, &
       &          real(beta), &
       &          vectout, npw)
     end if
   else if(gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET DATA USE_DEVICE_ADDR(work_buf,vectout,projections)
     call abi_gpu_xgemm(cplex, 'N','N', &
     &                  npw, ndat, nprojs_cur_blk, cone, &
     &                  c_loc(work_buf), npw, &
     &                  c_loc(projections(1,ibeg,1)), nprojs, &
     &                  beta, &
     &                  c_loc(vectout), npw)
     !$OMP END TARGET DATA
#endif
   end if

   beta = cone

   call xmpi_wait(req(1),ierr)
   call xmpi_wait(req(2),ierr)
   !call xmpi_waitall(req,ierr)

#ifdef HAVE_OPENMP_OFFLOAD
#ifndef HAVE_GPU_MPI
   ! If MPI is not GPU-aware, push received data to GPU
   !$OMP TARGET UPDATE TO(recv_buf) IF(gpu_option==ABI_GPU_OPENMP)
#endif
#endif

 end do

 if(modulo(iblock,2)==1) then
   if(gpu_option == ABI_GPU_DISABLED) then
     call DCOPY(cplex*npw*nprojs_cur_blk, recv_buf, 1, work_buf, 1)
   else if(gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET DATA USE_DEVICE_ADDR(work_buf,recv_buf)
     call copy_gpu_to_gpu(c_loc(work_buf), c_loc(recv_buf), INT(cplex, c_size_t)*npw*nprojs_last_blk*dp)
     !$OMP END TARGET DATA
#endif
   end if
 end if

#ifdef HAVE_OPENMP_OFFLOAD
 !$OMP TARGET EXIT DATA MAP(delete:projs_recv) IF(gpu_option==ABI_GPU_OPENMP)
#endif
 ABI_FREE(projs_recv)

end subroutine opernlb_gemm_distributed
!!***

!----------------------------------------------------------------------

subroutine opernlb_xgemm(cplex,transa,transb,npw,ndat,nprojs,alpha,a,lda,b,ldb,beta,c,ldc,&
&                       rank, nprocs,&
&                       nprojs_blk, nprojs_last_blk,&
&                       iblock,&
&                       gpu_option,use_distrib,use_sliced_gemms)

!Arguments ------------------------------------
 integer,intent(in) :: cplex,lda,ldb,ldc,npw,ndat,nprojs,gpu_option
 integer,intent(in) :: rank,nprocs,nprojs_blk,nprojs_last_blk
 integer,intent(in) :: iblock
 logical,intent(in) :: use_distrib,use_sliced_gemms
 complex(dp),intent(in) :: alpha,beta
 character(len=1),intent(in) :: transa,transb
 real(dp),target,intent(in) :: a(cplex,lda,nprojs), b(cplex,ldb,ndat)
 real(dp),target,intent(inout) :: c(cplex,ldc,ndat)

 integer :: ibeg
! *********************************************************************
 ibeg = 1
 if(use_sliced_gemms) ibeg = 1 + (iblock-1)*nprojs_blk

 if(use_distrib) then
   call opernlb_gemm_distributed(rank,nprocs,npw,ndat,&
   &                             nprojs,&
   &                             nprojs_blk,&
   &                             nprojs_last_blk,&
   &                             cplex,&
   &                             a,&
   &                             b,c,gpu_option)
 else
   if (gpu_option == ABI_GPU_DISABLED) then
     if(cplex==2) then
       call abi_zgemm_2r(transa,transb,npw,ndat,nprojs,alpha,&
       &    a,lda,b,ldb,beta,c,ldc)
     else ! cplex==1
       call DGEMM(transa,transb,npw,ndat,nprojs,real(alpha),&
       &    a,lda,b,ldb,real(beta),c,ldc)
     end if
   else if (gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET DATA USE_DEVICE_ADDR(a,b,c)
     call abi_gpu_xgemm(cplex,transa,transb,npw,ndat,nprojs,alpha,&
     &    c_loc(a),lda,c_loc(b(1,ibeg,1)),ldb,beta,c_loc(c),ldc)
     !$OMP END TARGET DATA
#endif
   end if
 end if

 end subroutine opernlb_xgemm

!----------------------------------------------------------------------

!!****f* m_opernlb_gemm/opernlb_gemm
!! NAME
!! opernlb_gemm
!!
!! FUNCTION
!! For a given wave-function |c>, get all projected scalars
!! <p_lmn|c> where |p_lmn> are non-local projectors
!!   With:
!!   <p_lmn|c>=4pi/sqrt(vol) (i)^l Sum_g[c(g).f_nl(g).Y_lm(g).exp(2pi.i.g.R)]
!!
!! INPUTS
!!  choice=chooses possible output:
!!         if choice>=0: compute projected scalars
!!         if choice<0: same as choice>0 but use already computed projected scalars
!!         if ABS(choice)>1, then compute additional quantities:
!!           2: compute projected scalars and derivatives wrt atm pos.
!!           3: compute projected scalars and derivatives wrt strains
!!           22: compute projected scalars and 2nd derivatives wrt atm pos. and q-vector.
!!           23: compute projected scalars, derivatives wrt atm pos. and derivatives wrt strains
!!           25: compute projected scalars and 3rd derivatives wrt atm pos. and two q-vectors.
!!           4, 24: compute projected scalars, derivatives wrt atm pos.
!!                  and 2nd derivatives wrt atm pos.
!!           33: compute projected scalars and 2nd derivatives wrt strain and q-vector.
!!           5,51,52: compute projected scalars and derivatives wrt wave vector k
!!           53: compute projected scalars and derivatives wrt wave vector k in direction idir+1 and idir+2 mod 3
!!           54: compute projected scalars, deriv. wrt atm pos., deriv. wrt wave vector k
!!               and 2nd derivatives wrt right wave vector k and atm pos.
!!           55: compute projected scalars, deriv. strains, deriv. wrt wave vector k
!!               and 2nd derivatives wrt right wave vector k and strain
!!           6: compute projected scalars, derivatives wrt atm pos., derivatives wrt strains,
!!              2nd derivatives wrt 2 strains and derivatives wrt strain and atm pos.
!!           7: not available
!!           8: compute projected scalars, derivatives wrt wave vector k
!!              and 2nd derivatives wrt 2 wave vectors k
!!  cplex=1 if <p_lmn|c> scalars are real or pure imaginary (equivalent to istwfk>1)
!!        2 if <p_lmn|c> scalars are complex
!!  dimffnl=second dimension of ffnl
!!  ffnl(npw,dimffnl,nlmn)= nonlocal quantities containing nonlocal form factors
!!  ia3=gives the number of the first atom in the subset presently treated
!!  idir=direction of the - atom to be moved in the case (choice=2,signs=2) or (choice=22,signs=2)
!!                        - k point direction in the case (choice=5,signs=2)
!!                        - strain component (1:6) in the case (choice=3,signs=2) or (choice=6,signs=1)
!!                        - strain component (1:9) in the case (choice=33,signs=2)
!!                        - (1:9) components to specify the atom to be moved and the second q-gradient
!!                          direction in the case (choice=25,signs=2)
!!  indlmn(6,nlmn)= array giving l,m,n,lm,ln,s for i=lmn
!!  istwf_k=option parameter that describes the storage of wfs
!!  kpg(npw,nkpg)=(k+G) components          for ikpg=1...3   (if nkpg=3 or 9)
!!       [(k+G)_a].[(k+G)_b] quantities for ikpg=4...9   (if nkpg=9)
!!       (k+G) Cartesian components for choice==33
!!  matblk=dimension of the array ph3d
!!  mpi_enreg=information about MPI parallelization
!!  ndgxdt=second dimension of dgxdt
!!  nd2gxdt=second dimension of d2gxdt
!!  nincat=number of atoms in the subset here treated
!!  nkpg=second dimension of array kpg (0, 3 or 9)
!!  nlmn=number of (l,m,n) numbers for current type of atom
!!  nloalg(3)=governs the choice of the algorithm for non-local operator.
!!  npw=number of plane waves in reciprocal space
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  ph3d(2,npw,matblk)=three-dimensional phase factors
!!  [qdir]= optional, direction of the q-gradient (only for choice=22 choice=25 and choice=33)
!!  signs=chooses possible output:
!!   signs=1: compute derivatives in all directions
!!   signs=2: compute derivative in direction IDIR only
!!            compatible only with 1st-order derivatives and "single" derivatives
!!  ucvol=unit cell volume (bohr^3)
!!  vect(2,npw*my_nspinor)=starting vector in reciprocal space
!!
!! OUTPUT
!!  if (choice>1) dgxdt(cplex,ndgxdt,nlmn,nincat,nspinor)=
!!     gradients of projected scalars wrt coords  (choice=2, 23, 4, 54, 6)
!!                                    wrt strains (choice=3, 23, 55)
!!                                    wrt k wave vect. (choice=5, 51, 52, 53, 54, 55, 8)
!!                                    wrt coords and q vect (choice=22)
!!                                    wrt coords and two q vects (choice=25)
!!                                    wrt strains and q vect (choice=33)
!!  if (choice=4, 24, 33, 54, 55, 6, 8) d2gxdt(cplex,nd2gxdt,nlmn,nincat,nspinor)=
!!     2nd grads of projected scalars wrt 2 coords (choice=4 or 24)
!!                                    wrt coords & k wave vect. (choice=54)
!!                                    wrt strains & k wave vect. (choice=55)
!!                                    wrt coords & strains (choice=6)
!!                                    wrt 2 strains (choice=6)
!!                                    wrt 2 k wave vect. (choice=8)
!!                                    wrt strains and q vect (choice=33)
!!     only compatible with signs=1
!!  cplex_dgxdt(ndgxdt) = used only when cplex = 1
!!             cplex_dgxdt(i) = 1 if dgxdt(1,i,:,:)   is real, 2 if it is pure imaginary
!!  cplex_d2gxdt(nd2gxdt) = used only when cplex = 1
!!             cplex_d2gxdt(i) = 1 if d2gxdt(1,i,:,:) is real, 2 if it is pure imaginary
!!
!! SIDE EFFECTS
!!  gx(cplex,nlmn,nincat,nspinor)= projected scalars - input if choice<0, output if choice>=0
!!
!! NOTES
!! 1-The openMP version is different from the standard version:
!!   the standard version is more effifient on one CPU core.
!! 2-Operate for one type of atom, and within this given type of atom,
!!   for a subset of at most nincat atoms.
!!
!! SOURCE
subroutine opernlb_gemm(choice,cplex,cplex_dgxdt,cplex_d2gxdt,cplex_fac,&
&       d2gxdtfac,d2gxdtfac_sij,dgxdtfac,dgxdtfac_sij,&
&       dimffnl,ffnl,gxfac,gxfac_sij,&
&       idir,indlmn,kpg,matblk,istwf_k,&
&       nd2gxdt,nd2gxdtfac,ndgxdt,ndgxdtfac,&
&       nkpg,npw,nspinor,signs,ucvol,ndat,ntypat,lmnmax,nattyp,&
&       is_kprime,iatom_only,atom_proj_shift,&
&       paw_opt,ph3d,&
&       nprojs,&
&       vectin,vectout,svectout,&
&       temp_realvec_r,temp_realvec_i,&
&       gpu_option,use_distrib)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,cplex,cplex_fac,idir,istwf_k,nd2gxdt,nd2gxdtfac
 integer,intent(in) :: ndgxdt,dimffnl,nkpg,lmnmax,ntypat,ndgxdtfac,matblk,npw,nspinor
 integer,intent(in) :: paw_opt,signs,ndat,iatom_only,atom_proj_shift
 integer,intent(in) :: nprojs
 real(dp),intent(in) :: ucvol
 integer,intent(in) :: gpu_option
 logical,intent(in) :: use_distrib,is_kprime
!arrays
 integer,intent(in)  :: indlmn(6,lmnmax,ntypat),nattyp(ntypat)
 integer,intent(in)  :: cplex_dgxdt(ndgxdt),cplex_d2gxdt(nd2gxdt)
 real(dp),intent(in) :: ffnl(npw,dimffnl,lmnmax,ntypat),kpg(npw,nkpg)
 real(dp),intent(in) :: ph3d(2,npw,matblk)
 real(dp),target,intent(in)  :: vectin(:,:)
 real(dp),target,intent(inout) :: vectout(:,:),svectout(:,:)
 real(dp),target,intent(in) :: d2gxdtfac(cplex_fac,nd2gxdtfac,nprojs,ndat*nspinor)
 real(dp),target,intent(in) :: d2gxdtfac_sij(cplex,nd2gxdtfac,nprojs,ndat*nspinor)
 real(dp),target,intent(inout) :: dgxdtfac(cplex_fac,ndgxdtfac*nprojs,ndat*nspinor)
 real(dp),target,intent(inout) :: dgxdtfac_sij(cplex,ndgxdtfac*nprojs,ndat*nspinor)
 real(dp),target,intent(in) :: gxfac(cplex_fac,nprojs,ndat*nspinor)
 real(dp),target,intent(in) :: gxfac_sij(cplex,nprojs,ndat*nspinor)
 real(dp),target,intent(out) :: temp_realvec_r(:),temp_realvec_i(:)

!Local variables-------------------------------
 integer :: idat,i,ik,nprojs_all,iproj,iplex
 integer :: projs_beg,projs_end,dprojs_beg,dprojs_end
 integer :: nprojs_blk,nprojs_last_blk,nprojs_cur_blk,rank,nprocs,iblock,nblocks
 logical :: use_sliced_gemms
 complex(dp) :: beta
 real(dp), ABI_CONTIGUOUS pointer :: projs(:,:,:),projs_r(:,:,:),projs_i(:,:,:)
 real(dp), ABI_CONTIGUOUS pointer :: dprojs(:,:,:),dprojs_r(:,:,:),dprojs_i(:,:,:)

 ik=1; if(is_kprime) ik=2
#ifndef HAVE_OPENMP_OFFLOAD
 ABI_UNUSED((/iproj,idat,iplex/))
#endif
 ABI_UNUSED(cplex_dgxdt)
 ABI_UNUSED(cplex_d2gxdt)
 ABI_UNUSED(d2gxdtfac)
 ABI_UNUSED(d2gxdtfac_sij)

 nprojs_all=nprojs
 if(iatom_only>0) then
   nprojs_all=0
   do i=1,ntypat
     nprojs_all = nprojs_all + count(indlmn(3,:,i)>0)*nattyp(i)
   end do
 end if
 nprojs_last_blk=nprojs_all
 iblock=1; nblocks=1
 use_sliced_gemms=.false.
 if(gemm_nonlop_block_size>1 .and. .not. use_distrib) then
   nblocks=gemm_nonlop_block_size
   use_sliced_gemms=.true.
 end if

 call refresh_projectors(npw,istwf_k,nprojs_all,ndgxdt,nd2gxdt,&
 &                       is_kprime,gpu_option)
 if(nprojs_all/=gemm_nonlop_kpt(ik)%nprojs) ABI_BUG("Problem")
 nprojs_blk = nprojs
 nprojs_last_blk = nprojs
 if(use_distrib) then
   rank=xmpi_comm_rank(gemm_nonlop_block_comm); nprocs=xmpi_comm_size(gemm_nonlop_block_comm)
   nprojs_blk      = gemm_nonlop_kpt(ik)%nprojs_blk
   nprojs_last_blk = gemm_nonlop_kpt(ik)%nprojs_last_blk
   iblock=rank+1
 else if(gemm_nonlop_block_size>1) then
    nprojs_blk = nprojs / gemm_nonlop_block_size
    nprojs_last_blk = nprojs_blk + modulo(nprojs,nprojs_blk)
 end if

 projs_beg=1; projs_end=nprojs;
 dprojs_beg=1; dprojs_end=max(1,nprojs*ndgxdt)
 if((choice==2 .and. signs==2)) then
   projs_beg=atom_proj_shift+1
   projs_end=projs_beg+nprojs-1
   dprojs_beg=atom_proj_shift*ndgxdt+1
   dprojs_end=dprojs_beg+nprojs*ndgxdt-1
 end if

 if(gemm_nonlop_block_size>1) then
   projs_beg=1; projs_end=nprojs_last_blk;
   dprojs_beg=1; dprojs_end=max(1,nprojs_last_blk*ndgxdt)
 end if

 if(istwf_k == 1) then
   projs => gemm_nonlop_kpt(ik)%projs(:,:,projs_beg:projs_end)
   if(ndgxdt>0)  dprojs => gemm_nonlop_kpt(ik)%dprojs(:,:,dprojs_beg:dprojs_end)
 else
   projs_r  => gemm_nonlop_kpt(ik)%projs_r(:,:,projs_beg:projs_end)
   projs_i  => gemm_nonlop_kpt(ik)%projs_i(:,:,projs_beg:projs_end)
   if(ndgxdt>0)  dprojs_r => gemm_nonlop_kpt(ik)%dprojs_r(:,:,dprojs_beg:dprojs_end)
   if(ndgxdt>0)  dprojs_i => gemm_nonlop_kpt(ik)%dprojs_i(:,:,dprojs_beg:dprojs_end)
 end if

 if(gemm_nonlop_kpt(ik)%ikpt/=gemm_nonlop_ikpt_this_proc_being_treated .or. use_sliced_gemms) then
   call prep_projectors(npw,lmnmax,ntypat,indlmn,nattyp,istwf_k,&
   &                    ucvol,ffnl,ph3d,dimffnl,matblk,&
   &                    nprojs_last_blk,is_kprime,gpu_option,iblock)
   gemm_nonlop_kpt(ik)%ikpt=gemm_nonlop_ikpt_this_proc_being_treated
 end if
 if(choice>1 .and. ndgxdt>0) then
   if(     nd2gxdt/=gemm_nonlop_kpt(ik)%ngrads2 &
   &   .or. ndgxdt/=gemm_nonlop_kpt(ik)%ngrads &
   &   .or. choice/=gemm_nonlop_kpt(ik)%choice &
   &   .or.   idir/=gemm_nonlop_kpt(ik)%idir &
   &   .or. use_sliced_gemms) then
     call prep_dprojectors(npw,lmnmax,ntypat,indlmn,nattyp,istwf_k,&
     &                    ucvol,ffnl,ph3d,kpg,nkpg,dimffnl,matblk,&
     &                    nprojs_last_blk,ndgxdt,nd2gxdt,choice,signs,idir,&
     &                    is_kprime,gpu_option,iblock)
     gemm_nonlop_kpt(ik)%choice = choice
     gemm_nonlop_kpt(ik)%idir = idir
   end if
 end if

 if(paw_opt == 3 .or. paw_opt == 4) then

   ! Get svectout from gxfac_sij
   if(cplex == 2) then

     ! With many blocks, GEMM results will be summed using beta=cone
     beta = czero

     do i=1,nblocks
       if(use_sliced_gemms .and. i>1) then
         call prep_projectors(npw,lmnmax,ntypat,indlmn,nattyp,istwf_k,&
         &                    ucvol,ffnl,ph3d,dimffnl,matblk,&
         &                    nprojs_last_blk,is_kprime,gpu_option,i)
         gemm_nonlop_kpt(ik)%ikpt=gemm_nonlop_ikpt_this_proc_being_treated
         if(choice>1 .and. ndgxdt>0) then
           call prep_dprojectors(npw,lmnmax,ntypat,indlmn,nattyp,istwf_k,&
           &                    ucvol,ffnl,ph3d,kpg,nkpg,dimffnl,matblk,&
           &                    nprojs_last_blk,ndgxdt,nd2gxdt,choice,signs,idir,&
           &                    is_kprime,gpu_option,i)
           gemm_nonlop_kpt(ik)%choice = choice
           gemm_nonlop_kpt(ik)%idir = idir
         end if
       end if

       nprojs_cur_blk=nprojs
       if(use_sliced_gemms) then
         if(i<nblocks) then
           nprojs_cur_blk=nprojs_blk
         else
           nprojs_cur_blk=nprojs_last_blk
         end if
       end if

       if(choice==1 .or. choice==7) then
         call opernlb_xgemm(cplex, 'N', 'N', npw, ndat*nspinor, nprojs_cur_blk, cone, &
         &    projs, npw,&
         &    gxfac_sij, nprojs, beta, svectout, npw,&
         &    rank, nprocs,&
         &    nprojs_blk, nprojs_last_blk, i,&
         &    gpu_option, use_distrib, use_sliced_gemms)
       else if(choice==2) then
         call opernlb_xgemm(cplex, 'N', 'N', npw, ndat*nspinor, nprojs_cur_blk, cone, &
         &    dprojs, npw, &
         &    gxfac_sij, nprojs, beta, svectout, npw,&
         &    rank, nprocs,&
         &    nprojs_blk, nprojs_last_blk, i,&
         &    gpu_option, use_distrib, use_sliced_gemms)
         call opernlb_xgemm(cplex, 'N', 'N', npw, ndat*nspinor, nprojs_cur_blk, cone, &
         &    projs, npw, &
         &    dgxdtfac_sij, nprojs, cone, svectout, npw,&
         &    rank, nprocs,&
         &    nprojs_blk, nprojs_last_blk, i,&
         &    gpu_option, use_distrib, use_sliced_gemms)
       else if(choice==3) then
         if(idir<=3) then
           if(gpu_option == ABI_GPU_DISABLED) then
             dgxdtfac_sij(:,:,:) = dgxdtfac_sij(:,:,:) - gxfac_sij(:,:,:)
           else if(gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
             !$OMP TARGET TEAMS DISTRIBUTE &
             !$OMP& MAP(to:gxfac_sij,dgxdtfac_sij) PRIVATE(idat)
             do idat=1,ndat
               !$OMP PARALLEL DO PRIVATE(iproj,iplex) COLLAPSE(2)
               do iproj=1,nprojs
                 do iplex=1,cplex
                   dgxdtfac_sij(iplex,iproj,idat) = dgxdtfac_sij(iplex,iproj,idat) - gxfac_sij(iplex,iproj,idat)
                 end do
               end do
             end do
#endif
           end if
         end if
         call opernlb_xgemm(cplex, 'N', 'N', npw, ndat*nspinor, nprojs_cur_blk, cone, &
         &    projs, npw, &
         &    dgxdtfac_sij, nprojs, beta, svectout, npw,&
         &    rank, nprocs,&
         &    nprojs_blk, nprojs_last_blk, i,&
         &    gpu_option, use_distrib, use_sliced_gemms)
         call opernlb_xgemm(cplex, 'N', 'N', npw, ndat*nspinor, nprojs_cur_blk, cone, &
         &    dprojs, npw, &
         &    gxfac_sij, nprojs, cone, svectout, npw,&
         &    rank, nprocs,&
         &    nprojs_blk, nprojs_last_blk, i,&
         &    gpu_option, use_distrib, use_sliced_gemms)
       else if(choice==5) then
         call opernlb_xgemm(cplex, 'N', 'N', npw, ndat*nspinor, nprojs_cur_blk, cone, &
         &    projs, npw, &
         &    dgxdtfac_sij, nprojs, beta, svectout, npw,&
         &    rank, nprocs,&
         &    nprojs_blk, nprojs_last_blk, i,&
         &    gpu_option, use_distrib, use_sliced_gemms)
         call opernlb_xgemm(cplex, 'N', 'N', npw, ndat*nspinor, nprojs_cur_blk, cone, &
         &    dprojs, npw, &
         &    gxfac_sij, nprojs, cone, svectout, npw,&
         &    rank, nprocs,&
         &    nprojs_blk, nprojs_last_blk, i,&
         &    gpu_option, use_distrib, use_sliced_gemms)
       else if(choice==51) then
         call opernlb_xgemm(cplex, 'N', 'N', npw, ndat*nspinor, nprojs_cur_blk, cone, &
         &    projs, npw, &
         &    dgxdtfac_sij, nprojs, beta, svectout, npw,&
         &    rank, nprocs,&
         &    nprojs_blk, nprojs_last_blk, i,&
         &    gpu_option, use_distrib, use_sliced_gemms)
       end if

       beta = cone
     end do
   else

     ! With many blocks, GEMM results will be summed using beta=cone
     beta = czero

     do i=1,nblocks
       if(use_sliced_gemms .and. i>1) then
         call prep_projectors(npw,lmnmax,ntypat,indlmn,nattyp,istwf_k,&
         &                    ucvol,ffnl,ph3d,dimffnl,matblk,&
         &                    nprojs_last_blk,is_kprime,gpu_option,i)
         gemm_nonlop_kpt(ik)%ikpt=gemm_nonlop_ikpt_this_proc_being_treated
         if(choice>1 .and. ndgxdt>0) then
           call prep_dprojectors(npw,lmnmax,ntypat,indlmn,nattyp,istwf_k,&
           &                    ucvol,ffnl,ph3d,kpg,nkpg,dimffnl,matblk,&
           &                    nprojs_last_blk,ndgxdt,nd2gxdt,choice,signs,idir,&
           &                    is_kprime,gpu_option,i)
           gemm_nonlop_kpt(ik)%choice = choice
           gemm_nonlop_kpt(ik)%idir = idir
         end if
       end if

       nprojs_cur_blk=nprojs
       if(use_sliced_gemms) then
         if(i<nblocks) then
           nprojs_cur_blk=nprojs_blk
         else
           nprojs_cur_blk=nprojs_last_blk
         end if
       end if

       call opernlb_xgemm(cplex, 'N', 'N', npw, ndat*nspinor, nprojs_cur_blk, cone, &
       &    projs_r, npw, &
       &    gxfac_sij, nprojs, beta, temp_realvec_r, npw,&
       &    rank, nprocs,&
       &    nprojs_blk, nprojs_last_blk, i,&
       &    gpu_option, use_distrib, use_sliced_gemms)
       call opernlb_xgemm(cplex, 'N', 'N', npw, ndat*nspinor, nprojs_cur_blk, cone, &
       &    projs_i, npw,&
       &    gxfac_sij, nprojs, beta, temp_realvec_i, npw,&
       &    rank, nprocs,&
       &    nprojs_blk, nprojs_last_blk, i,&
       &    gpu_option, use_distrib, use_sliced_gemms)

       beta=cone

     end do

     if(gpu_option == ABI_GPU_DISABLED) then
       svectout(1,1:npw*nspinor*ndat) = temp_realvec_r(1:npw*nspinor*ndat)
       svectout(2,1:npw*nspinor*ndat) = temp_realvec_i(1:npw*nspinor*ndat)
     else if(gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
       !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
       !$OMP& MAP(to:temp_realvec_r,temp_realvec_i,svectout) PRIVATE(i)
       do i=1, npw*nspinor*ndat
         svectout(1,i) = temp_realvec_r(i)
         svectout(2,i) = temp_realvec_i(i)
       end do
#endif
     end if

   end if ! cplex = 2

   if(choice /= 7 .and. choice /= 5 .and. choice/=51 .and. choice/=2 .and. choice/=3) then
     if(gpu_option == ABI_GPU_DISABLED) then
       svectout = svectout + vectin ! TODO understand this
     else if(gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
       !$OMP TARGET DATA USE_DEVICE_ADDR(vectin,svectout)
       call abi_gpu_xaxpy(1, 2*npw*nspinor*ndat, cone, &
       &    c_loc(vectin), 1, c_loc(svectout), 1)
       !$OMP END TARGET DATA
#endif
     end if
   end if
 end if  ! (paw_opt == 3 .or. paw_opt == 4)

 if(paw_opt == 0 .or. paw_opt == 1 .or. paw_opt == 4) then
   ! Get vectout from gxfac
   if(cplex_fac == 2) then

     ! With many blocks, GEMM results will be summed using beta=cone
     beta = czero

     do i=1,nblocks
       if(use_sliced_gemms .and. (i>1 .or. paw_opt == 3 .or. paw_opt == 4)) then
         call prep_projectors(npw,lmnmax,ntypat,indlmn,nattyp,istwf_k,&
         &                    ucvol,ffnl,ph3d,dimffnl,matblk,&
         &                    nprojs_last_blk,is_kprime,gpu_option,i)
         gemm_nonlop_kpt(ik)%ikpt=gemm_nonlop_ikpt_this_proc_being_treated
         if(choice>1 .and. ndgxdt>0) then
           call prep_dprojectors(npw,lmnmax,ntypat,indlmn,nattyp,istwf_k,&
           &                    ucvol,ffnl,ph3d,kpg,nkpg,dimffnl,matblk,&
           &                    nprojs_last_blk,ndgxdt,nd2gxdt,choice,signs,idir,&
           &                    is_kprime,gpu_option,i)
           gemm_nonlop_kpt(ik)%choice = choice
           gemm_nonlop_kpt(ik)%idir = idir
         end if
       end if

       nprojs_cur_blk=nprojs
       if(use_sliced_gemms) then
         if(i<nblocks) then
           nprojs_cur_blk=nprojs_blk
         else
           nprojs_cur_blk=nprojs_last_blk
         end if
       end if

       if(choice==1 .or. choice==7) then
         call opernlb_xgemm(cplex_fac, 'N', 'N', npw, ndat*nspinor, nprojs_cur_blk, cone, &
         &    projs, npw, &
         &    gxfac, nprojs, beta, vectout, npw,&
         &    rank, nprocs,&
         &    nprojs_blk, nprojs_last_blk, i,&
         &    gpu_option, use_distrib, use_sliced_gemms)
       else if(choice==2) then
         call opernlb_xgemm(cplex_fac, 'N', 'N', npw, ndat*nspinor, nprojs_cur_blk, cone, &
         &    dprojs, npw, &
         &    gxfac, nprojs, beta, vectout, npw,&
         &    rank, nprocs,&
         &    nprojs_blk, nprojs_last_blk, i,&
         &    gpu_option, use_distrib, use_sliced_gemms)
         call opernlb_xgemm(cplex_fac, 'N', 'N', npw, ndat*nspinor, nprojs_cur_blk, cone, &
         &    projs, npw, &
         &    dgxdtfac, nprojs, cone, vectout, npw,&
         &    rank, nprocs,&
         &    nprojs_blk, nprojs_last_blk, i,&
         &    gpu_option, use_distrib, use_sliced_gemms)
       else if(choice==3) then
         if(idir<=3) then
           if(gpu_option == ABI_GPU_DISABLED) then
             dgxdtfac(:,:,:) = dgxdtfac(:,:,:) - gxfac(:,:,:)
           else if(gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
             !$OMP TARGET TEAMS DISTRIBUTE &
             !$OMP& MAP(to:gxfac,dgxdtfac) PRIVATE(idat)
             do idat=1,ndat
               !$OMP PARALLEL DO PRIVATE(iproj,iplex) COLLAPSE(2)
               do iproj=1,nprojs
                 do iplex=1,cplex_fac
                   dgxdtfac(iplex,iproj,idat) = dgxdtfac(iplex,iproj,idat) - gxfac(iplex,iproj,idat)
                 end do
               end do
             end do
#endif
           end if
         end if
         call opernlb_xgemm(cplex_fac, 'N', 'N', npw, ndat*nspinor, nprojs_cur_blk, cone, &
         &    projs, npw, &
         &    dgxdtfac, nprojs, beta, vectout, npw,&
         &    rank, nprocs,&
         &    nprojs_blk, nprojs_last_blk, i,&
         &    gpu_option, use_distrib, use_sliced_gemms)
         call opernlb_xgemm(cplex_fac, 'N', 'N', npw, ndat*nspinor, nprojs_cur_blk, cone, &
         &    dprojs, npw, &
         &    gxfac, nprojs, cone, vectout, npw,&
         &    rank, nprocs,&
         &    nprojs_blk, nprojs_last_blk, i,&
         &    gpu_option, use_distrib, use_sliced_gemms)
       else if(choice==5) then
         call opernlb_xgemm(cplex_fac, 'N', 'N', npw, ndat*nspinor, nprojs_cur_blk, cone, &
         &    projs, npw, &
         &    dgxdtfac, nprojs, beta, vectout, npw,&
         &    rank, nprocs,&
         &    nprojs_blk, nprojs_last_blk, i,&
         &    gpu_option, use_distrib, use_sliced_gemms)
         call opernlb_xgemm(cplex_fac, 'N', 'N', npw, ndat*nspinor, nprojs_cur_blk, cone, &
         &    dprojs, npw, &
         &    gxfac, nprojs, cone, vectout, npw,&
         &    rank, nprocs,&
         &    nprojs_blk, nprojs_last_blk, i,&
         &    gpu_option, use_distrib, use_sliced_gemms)
       else if(choice==51) then
         call opernlb_xgemm(cplex_fac, 'N', 'N', npw, ndat*nspinor, nprojs_cur_blk, cone, &
         &    projs, npw, &
         &    dgxdtfac, nprojs, beta, vectout, npw,&
         &    rank, nprocs,&
         &    nprojs_blk, nprojs_last_blk, i,&
         &    gpu_option, use_distrib, use_sliced_gemms)
       end if
       beta = cone
     end do
   else

     ! With many blocks, GEMM results will be summed using beta=cone
     beta = czero

     do i=1,nblocks
       if(use_sliced_gemms .and. (i>1 .or. paw_opt == 3 .or. paw_opt == 4)) then
         call prep_projectors(npw,lmnmax,ntypat,indlmn,nattyp,istwf_k,&
         &                    ucvol,ffnl,ph3d,dimffnl,matblk,&
         &                    nprojs_last_blk,is_kprime,gpu_option,i)
         gemm_nonlop_kpt(ik)%ikpt=gemm_nonlop_ikpt_this_proc_being_treated
         if(choice>1 .and. ndgxdt>0) then
           call prep_dprojectors(npw,lmnmax,ntypat,indlmn,nattyp,istwf_k,&
           &                    ucvol,ffnl,ph3d,kpg,nkpg,dimffnl,matblk,&
           &                    nprojs_last_blk,ndgxdt,nd2gxdt,choice,signs,idir,&
           &                    is_kprime,gpu_option,i)
           gemm_nonlop_kpt(ik)%choice = choice
           gemm_nonlop_kpt(ik)%idir = idir
         end if
       end if

       nprojs_cur_blk=nprojs
       if(use_sliced_gemms) then
         if(i<nblocks) then
           nprojs_cur_blk=nprojs_blk
         else
           nprojs_cur_blk=nprojs_last_blk
         end if
       end if

       call opernlb_xgemm(cplex_fac, 'N', 'N', npw, ndat*nspinor, nprojs_cur_blk, cone, &
       &        projs_r, npw, &
       &        gxfac, nprojs, beta, temp_realvec_r, npw,&
       &        rank, nprocs,&
       &        nprojs_blk, nprojs_last_blk, i,&
       &        gpu_option, use_distrib, use_sliced_gemms)
       call opernlb_xgemm(cplex_fac, 'N', 'N', npw, ndat*nspinor, nprojs_cur_blk, cone, &
       &        projs_i, npw, &
       &        gxfac, nprojs, beta, temp_realvec_i, npw,&
       &        rank, nprocs,&
       &        nprojs_blk, nprojs_last_blk, i,&
       &        gpu_option, use_distrib, use_sliced_gemms)

       beta=cone

     end do

     if(gpu_option == ABI_GPU_DISABLED) then
       vectout(1,1:npw*nspinor*ndat) = temp_realvec_r(1:npw*nspinor*ndat)
       vectout(2,1:npw*nspinor*ndat) = temp_realvec_i(1:npw*nspinor*ndat)
     else if(gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
       !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
       !$OMP& MAP(to:temp_realvec_r,temp_realvec_i,vectout) PRIVATE(i)
       do i=1, npw*nspinor*ndat
         vectout(1,i) = temp_realvec_r(i)
         vectout(2,i) = temp_realvec_i(i)
       end do
#endif
     end if
   end if ! cplex_fac == 2
 end if  ! (paw_opt == 0 .or. paw_opt == 1 .or. paw_opt == 4)

end subroutine opernlb_gemm
!!***

end module m_opernlb_gemm
!!***
