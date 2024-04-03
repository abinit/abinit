!!****m* ABINIT/m_opernla_gemm
!! NAME
!!  m_opernla_gemm
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

module m_opernla_gemm

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_abi_linalg

 use defs_abitypes, only : MPI_type
 use m_time,        only : timab

#ifdef HAVE_FC_ISO_C_BINDING
 use, intrinsic :: iso_c_binding, only : c_ptr,c_loc
#endif

 implicit none

 private
!!***

 public :: opernla_gemm

contains
!!***


!!****f* m_gemm_nonlop/gemm_nonlop_distributed_gemm_opernla
!! NAME
!! gemm_nonlop_distributed_gemm_opernla
!!
!! FUNCTION
!! Distributed version of "opernla" GEMM called in gemm_nonlop.
!!
!! INPUTS
!!
!! SOURCE
subroutine gemm_nonlop_distributed_gemm_opernla(rank,nprocs,npwin,ndat,nspinor,&
&                                               nprojs_blk,nprojs_last_blk,nprojs_my_blk,cplex,beta,&
&                                               projs_local,projs_recv,projs,vectin,projections)
  integer,  intent(in)     :: rank,nprocs,npwin,ndat,nspinor
  integer,  intent(in)     :: nprojs_blk,nprojs_last_blk,nprojs_my_blk,cplex
  real(dp), intent(in)     :: projs(:,:,:),vectin(*)
  complex(dpc), intent(in) :: beta
  real(dp), intent(inout)  :: projs_local(:,:,:),projs_recv(:,:,:)
  real(dp), intent(out)    :: projections(:,:,:)

  !Local variables
  integer :: iblock,ibeg,iend,req(2),ierr,nprojs_cur_blk,rank_prev,rank_next
  integer :: gemm_nonlop_block_comm

! *************************************************************************

  rank_next=modulo(rank + 1,nprocs)
  rank_prev=rank - 1
  if(rank_prev == -1) rank_prev = nprocs - 1
  do iblock=1,nprocs

    if(rank+iblock == nprocs) then
      nprojs_cur_blk = nprojs_last_blk
    else
      nprojs_cur_blk = nprojs_blk
    end if

    if(iblock == 1) then
      call DCOPY(cplex*npwin*nprojs_my_blk , projs,      1, projs_local, 1)
    else
      call DCOPY(cplex*npwin*nprojs_cur_blk, projs_recv, 1, projs_local, 1)
    end if

    if(iblock < nprocs) then
      ABI_BUG("FIX gemm_nonlop_block_comm damn it !")
      call xmpi_isend(projs_local,rank_prev,iblock,gemm_nonlop_block_comm,req(1),ierr)
      call xmpi_irecv(projs_recv,rank_next,iblock,gemm_nonlop_block_comm,req(2),ierr)
    end if

    ibeg = 1 + modulo(rank+iblock-1,nprocs)*nprojs_blk
    iend = ibeg+nprojs_cur_blk-1
    if(cplex==2) then
      call abi_zgemm_2r('C', 'N', nprojs_cur_blk, ndat*nspinor, npwin, cone, &
      &                projs_local, npwin,&
      &                vectin, npwin, beta, projections(:,ibeg:iend,:), nprojs_cur_blk)
    else
      call DGEMM('T', 'N', nprojs_cur_blk, ndat*nspinor, npwin, one, &
      &          projs_local, npwin, &
      &          vectin, npwin, real(beta), projections(:,ibeg:iend,:), nprojs_cur_blk)
    end if

    if(iblock < nprocs) then
      call xmpi_waitall(req,ierr)
    end if

  end do
end subroutine gemm_nonlop_distributed_gemm_opernla

!----------------------------------------------------------------------

#ifdef HAVE_OPENMP_OFFLOAD
!!****f* m_gemm_nonlop_ompgpu/gemm_nonlop_ompgpu_distributed_gemm_opernla
!! NAME
!! gemm_nonlop_ompgpu_distributed_gemm_opernla
!!
!! FUNCTION
!! Distributed version of "opernla" GEMM called in gemm_nonlop.
!!
!! INPUTS
!!
!! SOURCE
 subroutine gemm_nonlop_ompgpu_distributed_gemm_opernla(rank,nprocs,npwin,ndat,nspinor,&
 &                                               nprojs,nprojs_blk,nprojs_last_blk,nprojs_my_blk,cplex,beta,&
 &                                               projs_local,projs_recv,vectin,cprojs)
   integer,  intent(in)     :: rank,nprocs,npwin,ndat,nspinor
   integer,  intent(in)     :: nprojs,nprojs_blk,nprojs_last_blk,nprojs_my_blk,cplex
   real(dp), intent(in), target     :: vectin(*)
   complex(dpc), intent(in) :: beta
   real(dp), intent(inout), target  :: projs_local(cplex,npwin,nprojs_last_blk),projs_recv(cplex,npwin,nprojs_last_blk)
   real(dp), intent(out), target  :: cprojs(:,:,:)

   !Local variables
   integer :: iblock,ibeg,iend,req(2),ierr,nprojs_cur_blk,rank_prev,rank_next,idat
   real(dp), ABI_CONTIGUOUS pointer :: recv_buf(:,:,:), work_buf(:,:,:)
   real(dp), ABI_CONTIGUOUS pointer :: recv_buf_f(:), work_buf_f(:)


   integer :: gemm_nonlop_block_comm
! *************************************************************************

   ABI_BUG("broken")
   rank_next=modulo(rank + 1,nprocs)
   rank_prev=rank - 1
   if(rank_prev == -1) rank_prev = nprocs - 1
   do iblock=1,nprocs

     call xmpi_barrier(gemm_nonlop_block_comm)
     if(rank+iblock == nprocs) then
       nprojs_cur_blk = nprojs_last_blk
     else
       nprojs_cur_blk = nprojs_blk
     end if

     if(modulo(iblock,2)==1) then
       work_buf => projs_local(1:cplex,1:npwin,1:nprojs_last_blk)
       recv_buf => projs_recv(1:cplex,1:npwin,1:nprojs_last_blk)
     else
       work_buf => projs_recv(1:cplex,1:npwin,1:nprojs_last_blk)
       recv_buf => projs_local(1:cplex,1:npwin,1:nprojs_last_blk)
     end if

#ifndef HAVE_GPU_MPI

       ! GPU-aware MPI not available : perform MPI comms on CPU
       call xmpi_isend(work_buf,rank_prev,iblock,gemm_nonlop_block_comm,req(1),ierr)
       call xmpi_irecv(recv_buf,rank_next,iblock,gemm_nonlop_block_comm,req(2),ierr)
       !$OMP TARGET UPDATE TO(work_buf)

#else

       ! GPU-aware MPI available : pass GPU buffers to MPI
       !$OMP TARGET DATA USE_DEVICE_PTR(work_buf,recv_buf)
       call c_f_pointer(c_loc(work_buf), work_buf_f, [cplex*npwin*nprojs_last_blk])
       call c_f_pointer(c_loc(recv_buf), recv_buf_f, [cplex*npwin*nprojs_last_blk])
       call MPI_ISEND(work_buf_f,cplex*npwin*nprojs_cur_blk,MPI_DOUBLE_PRECISION,&
       &    rank_prev,iblock,gemm_nonlop_block_comm,req(1),ierr)
       call MPI_IRECV(recv_buf_f,cplex*npwin*nprojs_cur_blk,MPI_DOUBLE_PRECISION,&
       &    rank_next,iblock,gemm_nonlop_block_comm,req(2),ierr)
       !$OMP END TARGET DATA

#endif


     ibeg = 1 + modulo(rank+iblock-1,nprocs)*nprojs_blk
     iend = ibeg+nprojs_cur_blk-1
     if(cplex == 2) then
       !$OMP TARGET DATA USE_DEVICE_PTR(work_buf,vectin,cprojs)
       call abi_gpu_xgemm(cplex, 'C','N', &
               nprojs_cur_blk, ndat*nspinor, npwin, cone, &
               c_loc(work_buf), npwin, &
               c_loc(vectin), npwin, &
               beta, &
               c_loc(cprojs(1,ibeg,1)), nprojs)
       !$OMP END TARGET DATA
     else
       call DGEMM('T', 'N', nprojs_cur_blk, ndat*nspinor, npwin, one, &
       &          work_buf, npwin, &
       &          vectin, npwin, real(beta), cprojs(:,ibeg:iend,:), nprojs_cur_blk)
     end if

     call xmpi_waitall(req,ierr)

     call xmpi_barrier(gemm_nonlop_block_comm)
   end do

   if(modulo(iblock,2)==1) then
#ifndef HAVE_GPU_MPI
       call DCOPY(cplex*npwin*nprojs_cur_blk, recv_buf, 1, work_buf, 1)
#else
       !$OMP TARGET DATA USE_DEVICE_PTR(work_buf,recv_buf)
       call copy_gpu_to_gpu(c_loc(work_buf), c_loc(recv_buf), INT(cplex, c_size_t)*npwin*nprojs_last_blk*dp)
       !$OMP END TARGET DATA
#endif
   end if

 end subroutine gemm_nonlop_ompgpu_distributed_gemm_opernla
!!***
#endif

!----------------------------------------------------------------------

 subroutine opernl_xgemm(cplx,transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc,&
 &                      gpu_option,use_distrib)

!Arguments ------------------------------------
 integer,intent(in) :: cplx,lda,ldb,ldc,m,n,k,gpu_option
 logical,intent(in) :: use_distrib
 complex(dpc),intent(in) :: alpha,beta
 character(len=1),intent(in) :: transa,transb
 real(dp),target,intent(in) :: a(*),b(*)
 real(dp),target,intent(inout) :: c(*)

! *********************************************************************

 if (gpu_option == ABI_GPU_DISABLED) then
   if(cplx==2) then
     call abi_zgemm_2r(transa,transb,m,n,k,alpha,&
     &    a,lda,b,ldb,beta,c,ldc)
   else ! cplx==1
     call DGEMM(transa,transb,m,n,k,real(alpha),&
     &    a,lda,b,ldb,real(beta),c,ldc)
   end if
 else if (gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET DATA USE_DEVICE_PTR(a,b,c)
   call abi_gpu_xgemm(cplx,transa,transb,m,n,k,alpha,&
   &    c_loc(a),lda,c_loc(b),ldb,beta,c_loc(c),ldc)
   !$OMP END TARGET DATA
#endif
 end if

 end subroutine opernl_xgemm

!----------------------------------------------------------------------

!!***
!!****f* ABINIT/opernla_gemm
!! NAME
!! opernla_gemm
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

subroutine opernla_gemm(choice,cplex,cplex_dgxdt,cplex_d2gxdt,d2gxdt,dgxdt,gx,&
&       idir,istwf_k,mpi_enreg,nd2gxdt,ndgxdt,&
&       npw,nspinor,signs,ndat,rank,&
&       cpopt,nprocs,&
&       nprojs,nprojs_blk,nprojs_my_blk,nprojs_last_blk,&
&       vectin,projs,dprojs,d2projs,&
&       projs_r,projs_i,dprojs_r,dprojs_i,d2projs_r,d2projs_i,temp_realvec,&
&       projs_local,projs_recv,&
&       gpu_option,use_distrib)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,cplex,idir,istwf_k,nd2gxdt
 integer,intent(in) :: ndgxdt,npw,nspinor,signs,ndat,rank
 integer,intent(in) :: cpopt,nprocs
 integer,intent(in) :: nprojs,nprojs_blk,nprojs_my_blk,nprojs_last_blk
 type(MPI_type),intent(in) :: mpi_enreg
 integer,intent(in) :: gpu_option
 logical,intent(in) :: use_distrib
!arrays
 integer,intent(out) :: cplex_dgxdt(ndgxdt),cplex_d2gxdt(nd2gxdt)
 real(dp),target,intent(in) :: vectin(:,:)
 real(dp),target,intent(inout) :: d2gxdt(cplex,nd2gxdt,nprojs,ndat*nspinor)
 real(dp),target,intent(inout) :: dgxdt(cplex,ndgxdt*nprojs,ndat*nspinor)
 real(dp),target,intent(inout) :: gx(cplex,nprojs,ndat*nspinor)
 real(dp),target,intent(inout) :: projs(:,:,:),projs_r(:,:,:),projs_i(:,:,:)
 real(dp),target,intent(inout) :: dprojs(:,:,:),dprojs_r(:,:,:),dprojs_i(:,:,:)
 real(dp),target,intent(inout) :: d2projs(:,:,:),d2projs_r(:,:,:),d2projs_i(:,:,:)
 real(dp),target,intent(inout) :: projs_local(:,:,:),projs_recv(:,:,:),temp_realvec(:)

!Local variables-------------------------------
 integer :: idat,ierr,i,i1,i2,i3

 ABI_UNUSED(d2gxdt)
 ABI_UNUSED(cplex_d2gxdt)
 ABI_UNUSED(idir)

 if(cplex == 2) then
   if(.not. use_distrib) then
     if(cpopt<=1) then
       call opernl_xgemm(cplex, 'C', 'N', nprojs, ndat*nspinor, npw, cone, &
       &    projs, npw,&
       &    vectin, npw, czero, gx, nprojs,&
       &    gpu_option, use_distrib)
     end if

     if(ndgxdt>0 .and. cpopt<=3) then
       call opernl_xgemm(cplex, 'C', 'N', ndgxdt*nprojs, ndat*nspinor, npw, cone, &
       &    dprojs, npw,&
       &    vectin, npw, czero, dgxdt, ndgxdt*nprojs,&
       &    gpu_option, use_distrib)
     end if

     if(nd2gxdt>0) then
       call opernl_xgemm(cplex, 'C', 'N', nd2gxdt*nprojs, ndat*nspinor, npw, cone, &
       &    d2projs, npw,&
       &    vectin, npw, czero, d2gxdt, nd2gxdt*nprojs,&
       &    gpu_option, use_distrib)
     end if
   else
     if(gpu_option == ABI_GPU_DISABLED) then
       if(cpopt<=1) then
         call gemm_nonlop_distributed_gemm_opernla(rank,nprocs,npw,ndat,nspinor,&
         &                                         nprojs_blk,&
         &                                         nprojs_last_blk,&
         &                                         nprojs_my_blk,cplex,czero,&
         &                                         projs_local,projs_recv,&
         &                                         projs,&
         &                                         vectin,gx)
       end if
       if(ndgxdt>0 .and. cpopt<=3) then
         call gemm_nonlop_distributed_gemm_opernla(rank,nprocs,npw,ndat,nspinor,&
         &                                         ndgxdt*nprojs_blk,&
         &                                         ndgxdt*nprojs_last_blk,&
         &                                         ndgxdt*nprojs_my_blk,cplex,czero,&
         &                                         projs_local,projs_recv,&
         &                                         dprojs,&
         &                                         vectin,dgxdt)
       end if
     else if(gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
       if(cpopt<=1) then
         call gemm_nonlop_ompgpu_distributed_gemm_opernla(rank,nprocs,npw,ndat,nspinor,&
         &                                         nprojs,&
         &                                         nprojs_blk,&
         &                                         nprojs_last_blk,&
         &                                         nprojs_my_blk,cplex,czero,&
         &                                         projs,projs_recv,&
         &                                         vectin,gx)
       end if
       if(ndgxdt>0 .and. cpopt<=3) then
         call gemm_nonlop_ompgpu_distributed_gemm_opernla(rank,nprocs,npw,ndat,nspinor,&
         &                                         ndgxdt*nprojs,&
         &                                         ndgxdt*nprojs_blk,&
         &                                         ndgxdt*nprojs_last_blk,&
         &                                         ndgxdt*nprojs_my_blk,cplex,czero,&
         &                                         dprojs,projs_recv,&
         &                                         vectin,dgxdt)
       end if
#endif
     end if
   end if

 else ! cplex==1

   ! only compute real part of gx = P^* psi => gx_r = P_r^T psi_r + P_i^T psi_i
   if(gpu_option == ABI_GPU_DISABLED) then
     temp_realvec(1:npw*nspinor*ndat) = vectin(1,1:npw*nspinor*ndat)
   else if(gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
     !$OMP& MAP(to:temp_realvec,vectin) PRIVATE(i)
     do i=1, npw*nspinor*ndat
       temp_realvec(i) = vectin(1,i)
     end do
#endif
   end if

   if(istwf_k == 2 .and. mpi_enreg%me_g0_fft == 1) then
     if(gpu_option == ABI_GPU_DISABLED) then
       do idat=1, ndat*nspinor
         temp_realvec(1+(idat-1)*npw) = temp_realvec(1+(idat-1)*npw)/2
       end do
     else if(gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
       !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
       !$OMP& MAP(to:temp_realvec) PRIVATE(idat)
       do idat=1, ndat*nspinor
         temp_realvec(1+(idat-1)*npw) = temp_realvec(1+(idat-1)*npw)/2
       end do
#endif
     end if
   end if
   if(.not. use_distrib) then
     if(cpopt<=1) then
       call opernl_xgemm(cplex, 'T', 'N', nprojs, ndat*nspinor, npw, cone, &
       &    projs_r, npw, &
       &    temp_realvec, npw, czero, gx, nprojs,&
       &    gpu_option, use_distrib)
     end if
     if(ndgxdt>0 .and. cpopt<=3) then
       call opernl_xgemm(cplex, 'T', 'N', ndgxdt*nprojs, ndat*nspinor, npw, cone, &
       &    dprojs_r, npw, &
       &    temp_realvec, npw, czero, dgxdt, ndgxdt*nprojs,&
       &    gpu_option, use_distrib)
     end if
   else
     if(gpu_option == ABI_GPU_DISABLED) then
       if(cpopt<=1) then
         call gemm_nonlop_distributed_gemm_opernla(rank,nprocs,npw,ndat,nspinor,&
         &                                         nprojs_blk,&
         &                                         nprojs_last_blk,&
         &                                         nprojs_my_blk,cplex,czero,&
         &                                         projs_local,projs_recv,&
         &                                         projs_r,&
         &                                         temp_realvec,gx)
       end if
       if(ndgxdt>0 .and. cpopt<=3) then
         call gemm_nonlop_distributed_gemm_opernla(rank,nprocs,npw,ndat,nspinor,&
         &                                         ndgxdt*nprojs_blk,&
         &                                         ndgxdt*nprojs_last_blk,&
         &                                         ndgxdt*nprojs_my_blk,cplex,czero,&
         &                                         projs_local,projs_recv,&
         &                                         dprojs_r,&
         &                                         temp_realvec,dgxdt)
       end if
     else if(gpu_option == ABI_GPU_OPENMP) then
       ABI_BUG("Case not handled wit OpenMP GPU")
     end if
   end if

   if(gpu_option == ABI_GPU_DISABLED) then
     temp_realvec(1:npw*nspinor*ndat) = vectin(2,1:npw*nspinor*ndat)
   else if(gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
     !$OMP& MAP(to:temp_realvec,vectin) PRIVATE(i)
     do i=1, npw*nspinor*ndat
       temp_realvec(i) = vectin(2,i)
     end do
#endif
   end if

   if(istwf_k == 2 .and. mpi_enreg%me_g0_fft == 1) then
     if(gpu_option == ABI_GPU_DISABLED) then
       do idat=1, ndat*nspinor
         temp_realvec(1+(idat-1)*npw) = zero
       end do
     else if(gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
       !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
       !$OMP& MAP(to:temp_realvec) PRIVATE(idat)
       do idat=1, ndat*nspinor
         temp_realvec(1+(idat-1)*npw) = zero
       end do
#endif
     end if
   end if
   if(.not. use_distrib) then
     if(cpopt<=1) then
       call opernl_xgemm(cplex, 'T', 'N', nprojs, ndat*nspinor, npw, cone, &
       &    projs_i, npw, &
       &    temp_realvec, npw, cone , gx, nprojs,&
       &    gpu_option, use_distrib)

       if(gpu_option == ABI_GPU_DISABLED) then
         gx = gx * 2
       else if(gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
         !$OMP TARGET DATA USE_DEVICE_PTR(gx)
         call abi_gpu_xscal(cplex, nprojs*nspinor*ndat, ctwo, c_loc(gx), 1)
         !$OMP END TARGET DATA
#endif
       end if

     end if
     if(ndgxdt>0) then
       call opernl_xgemm(cplex, 'T', 'N', ndgxdt*nprojs, ndat*nspinor, npw, cone, &
       &    dprojs_i, npw, &
       &    temp_realvec, npw, cone , dgxdt, ndgxdt*nprojs,&
       &    gpu_option, use_distrib)

       if(gpu_option == ABI_GPU_DISABLED) then
         dgxdt = dgxdt * 2
       else if(gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
         !$OMP TARGET DATA USE_DEVICE_PTR(dgxdt)
         call abi_gpu_xscal(cplex, ndgxdt*nprojs*nspinor*ndat, ctwo, c_loc(dgxdt), 1)
         !$OMP END TARGET DATA
#endif
       end if
     end if
   else ! use_distrib
     if(gpu_option == ABI_GPU_DISABLED) then
       if(cpopt<=1) then
         call gemm_nonlop_distributed_gemm_opernla(rank,nprocs,npw,ndat,nspinor,&
         &                                         nprojs_blk,&
         &                                         nprojs_last_blk,&
         &                                         nprojs_my_blk,cplex,cone,&
         &                                         projs_local,projs_recv,&
         &                                         projs_i,&
         &                                         temp_realvec,gx)
         gx = gx * 2
       end if
       if(ndgxdt>0 .and. cpopt<=3) then
         call gemm_nonlop_distributed_gemm_opernla(rank,nprocs,npw,ndat,nspinor,&
         &                                         nprojs_blk,&
         &                                         nprojs_last_blk,&
         &                                         ndgxdt*nprojs_my_blk,cplex,cone,&
         &                                         projs_local,projs_recv,&
         &                                         dprojs_i,&
         &                                         temp_realvec,dgxdt)
         dgxdt = dgxdt * 2
       end if
     end if
   end if
 end if ! cplex == 2
 if(gpu_option == ABI_GPU_DISABLED) then
   call xmpi_sum(gx,mpi_enreg%comm_fft,ierr)
   if (choice>1) then
     call xmpi_sum(dgxdt,mpi_enreg%comm_fft,ierr)
   end if
 end if

end subroutine opernla_gemm
end module m_opernla_gemm
