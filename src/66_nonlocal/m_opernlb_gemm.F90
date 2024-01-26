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
 use m_xmpi
 use m_abi_linalg

 use defs_abitypes, only : MPI_type
 use m_time,        only : timab

 implicit none

 private
!!***

 public :: opernlb_gemm

contains
!!***


!----------------------------------------------------------------------

!!****f* m_gemm_nonlop/gemm_nonlop_distributed_gemm_opernlb
!! NAME
!! gemm_nonlop_distributed_gemm_opernlb
!!
!! FUNCTION
!! Distributed version of "opernlb" GEMM called in gemm_nonlop.
!!
!! INPUTS
!!
!! SOURCE
 subroutine gemm_nonlop_distributed_gemm_opernlb(rank,nprocs,npwout,ndat,nspinor,&
 &                                               nprojs_blk,nprojs_last_blk,nprojs_my_blk,cplex,&
 &                                               projs_local,projs_recv,projs,projections,vectout)
   integer,  intent(in)     :: rank,nprocs,npwout,ndat,nspinor
   integer,  intent(in)     :: nprojs_blk,nprojs_last_blk,nprojs_my_blk,cplex
   real(dp), intent(in)     :: projs(:,:,:),projections(:,:,:)
   real(dp), intent(inout)  :: projs_local(:,:,:),projs_recv(:,:,:)
   real(dp), intent(out)    :: vectout(*)

   !Local variables
   integer :: iblock,ibeg,iend,req(2),ierr,nprojs_cur_blk,rank_prev,rank_next
   complex(dpc) :: beta
   integer :: gemm_nonlop_block_comm

! *************************************************************************

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

     if(iblock == 1) then
       call DCOPY(cplex*npwout*nprojs_my_blk, projs, 1, projs_local, 1)
     else
       call DCOPY(cplex*npwout*nprojs_cur_blk, projs_recv, 1, projs_local, 1)
     end if

     if(iblock < nprocs) then
       ABI_BUG("FIX gemm_nonlop_block_comm damn it !")
       call xmpi_isend(projs_local,rank_prev,iblock,gemm_nonlop_block_comm,req(1),ierr)
       call xmpi_irecv(projs_recv,rank_next,iblock,gemm_nonlop_block_comm,req(2),ierr)
     end if

     ibeg = 1 + modulo(rank+iblock-1,nprocs)*nprojs_blk
     iend = ibeg+nprojs_cur_blk-1

     if(cplex==2) then
       call abi_zgemm_2r('N', 'N', npwout, ndat*nspinor, nprojs_cur_blk, cone, &
       &                 projs_local, npwout, &
       &                 projections(:,ibeg:iend,:), nprojs_cur_blk, beta, vectout, npwout)
     else
       call DGEMM('N', 'N', npwout, ndat*nspinor, nprojs_cur_blk, one, &
       &          projs_local, npwout, &
       &          projections(:,ibeg:iend,:), nprojs_cur_blk, real(beta), vectout, npwout)
     end if

     beta = cone

     if(iblock < nprocs) then
       call xmpi_waitall(req,ierr)
     end if

   end do
 end subroutine gemm_nonlop_distributed_gemm_opernlb
!!***

!----------------------------------------------------------------------

!!***
!!****f* ABINIT/opernlb_gemm
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
&       d2gxdtfac,d2gxdtfac_sij,dgxdtfac,dgxdtfac_sij,gxfac,gxfac_sij,&
&       idir,istwf_k,mpi_enreg,nd2gxdt,ndgxdt,&
&       npw,nspinor,signs,ndat,rank,&
&       cpopt,nprocs,paw_opt,&
&       nprojs,nprojs_blk,nprojs_my_blk,nprojs_last_blk,&
&       vectin,vectout,svectout,projs,dprojs,&
&       projs_r,projs_i,dprojs_r,dprojs_i,temp_realvec,&
&       projs_local,projs_recv,dprojs_local,dprojs_recv,&
&       use_distrib)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,cplex,cplex_fac,idir,istwf_k,nd2gxdt
 integer,intent(in) :: ndgxdt,npw,nspinor,signs,ndat,rank
 integer,intent(in) :: cpopt,nprocs,paw_opt
 integer,intent(in) :: nprojs,nprojs_blk,nprojs_my_blk,nprojs_last_blk
 type(MPI_type),intent(in) :: mpi_enreg
 logical :: use_distrib
!arrays
 integer,intent(in)  :: cplex_dgxdt(ndgxdt),cplex_d2gxdt(nd2gxdt)
 real(dp),intent(in)  :: vectin(:,:)
 real(dp),intent(out) :: vectout(:,:),svectout(:,:)
 real(dp),intent(in) :: d2gxdtfac(cplex,nd2gxdt,nprojs,ndat*nspinor)
 real(dp),intent(in) :: d2gxdtfac_sij(cplex_fac,nd2gxdt,nprojs,ndat*nspinor)
 real(dp),intent(inout) :: dgxdtfac(cplex,ndgxdt*nprojs,ndat*nspinor)
 real(dp),intent(inout) :: dgxdtfac_sij(cplex_fac,ndgxdt*nprojs,ndat*nspinor)
 real(dp),intent(in) :: gxfac(cplex,nprojs,ndat*nspinor)
 real(dp),intent(in) :: gxfac_sij(cplex_fac,nprojs,ndat*nspinor)
 real(dp),intent(in) :: projs(:,:,:),projs_r(:,:,:),projs_i(:,:,:)
 real(dp),intent(in) :: dprojs(:,:,:),dprojs_r(:,:,:),dprojs_i(:,:,:)
 real(dp),intent(out) :: projs_local(:,:,:),projs_recv(:,:,:),dprojs_local(:,:,:),dprojs_recv(:,:,:),temp_realvec(:)

!Local variables-------------------------------
 complex(dpc), parameter :: cminusone  = (-1._dp,0._dp)
 integer :: idat,ierr


 if(paw_opt == 3 .or. paw_opt == 4) then

   ! Get svectout from gxfac_sij
   if(cplex == 2) then

     if(.not. use_distrib) then
       if(choice==1) then
         call abi_zgemm_2r('N', 'N', npw, ndat*nspinor, nprojs, cone, &
         &    projs, npw, &
         &    gxfac_sij, nprojs, czero, svectout, npw)
       else if(choice==2) then
         call abi_zgemm_2r('N', 'N', npw, ndat*nspinor, nprojs, cone, &
         &    dprojs, npw, &
         &    gxfac_sij, nprojs, czero, svectout, npw)
         call abi_zgemm_2r('N', 'N', npw, ndat*nspinor, nprojs, cone, &
         &    projs, npw, &
         &    dgxdtfac_sij, nprojs, cone, svectout, npw)
       else if(choice==3) then
         if(idir<=3) then
           dgxdtfac_sij(:,:,:) = dgxdtfac_sij(:,:,:) - gxfac_sij(:,:,:)
         end if
         call abi_zgemm_2r('N', 'N', npw, ndat*nspinor, nprojs, cone, &
         &    projs, npw, &
         &    dgxdtfac_sij, nprojs, czero, svectout, npw)
         call abi_zgemm_2r('N', 'N', npw, ndat*nspinor, nprojs, cminusone, &
         &    dprojs, npw, &
         &    gxfac_sij, nprojs, cone, svectout, npw)
       else if(choice==5) then
         call abi_zgemm_2r('N', 'N', npw, ndat*nspinor, nprojs, cone, &
         &    projs, npw, &
         &    dgxdtfac_sij, nprojs, czero, svectout, npw)
         call abi_zgemm_2r('N', 'N', npw, ndat*nspinor, nprojs, cone, &
         &    dprojs, npw, &
         &    gxfac_sij, nprojs, cone, svectout, npw)
       else if(choice==51) then
         call abi_zgemm_2r('N', 'N', npw, ndat*nspinor, nprojs, cone, &
         &    projs, npw, &
         &    dgxdtfac_sij, nprojs, czero, svectout, npw)
       end if
     else
       call gemm_nonlop_distributed_gemm_opernlb(rank,nprocs,npw,ndat,nspinor,&
       &                                         nprojs_blk,&
       &                                         nprojs_last_blk,&
       &                                         nprojs_my_blk,cplex,&
       &                                         projs_local,projs_recv,&
       &                                         projs,&
       &                                         gxfac_sij,svectout)
     end if
   else

     if(.not. use_distrib) then
       call DGEMM('N', 'N', npw, ndat*nspinor, nprojs, one, &
       &          projs_r, npw, &
       &          gxfac_sij, nprojs, zero, temp_realvec, npw)
       svectout(1,1:npw*nspinor*ndat) = temp_realvec(1:npw*nspinor*ndat)
       call DGEMM('N', 'N', npw, ndat*nspinor, nprojs, one, &
       &          projs_i, npw,&
       &          gxfac_sij, nprojs, zero, temp_realvec, npw)
       svectout(2,1:npw*nspinor*ndat) = temp_realvec(1:npw*nspinor*ndat)
     else
       call gemm_nonlop_distributed_gemm_opernlb(rank,nprocs,npw,ndat,nspinor,&
       &                                         nprojs_blk,&
       &                                         nprojs_last_blk,&
       &                                         nprojs_my_blk,cplex,&
       &                                         projs_local,projs_recv,&
       &                                         projs_r,&
       &                                         gxfac_sij,temp_realvec)
       svectout(1,1:npw*nspinor*ndat) = temp_realvec(1:npw*nspinor*ndat)
       call gemm_nonlop_distributed_gemm_opernlb(rank,nprocs,npw,ndat,nspinor,&
       &                                         nprojs_blk,&
       &                                         nprojs_last_blk,&
       &                                         nprojs_my_blk,cplex,&
       &                                         projs_local,projs_recv,&
       &                                         projs_i,&
       &                                         gxfac_sij,temp_realvec)
       svectout(2,1:npw*nspinor*ndat) = temp_realvec(1:npw*nspinor*ndat)
     end if

   end if ! cplex = 2
   if(choice /= 7 .and. choice /= 5 .and. choice/=51 .and. choice/=2 .and. choice/=3) svectout = svectout + vectin ! TODO understand this

 end if  ! (paw_opt == 3 .or. paw_opt == 4)

 if(paw_opt == 0 .or. paw_opt == 1 .or. paw_opt == 4) then
   ! Get vectout from vnl_projections
   if(cplex_fac == 2) then

     if(.not. use_distrib) then
       if(choice==1) then
         call abi_zgemm_2r('N', 'N', npw, ndat*nspinor, nprojs, cone, &
         &    projs, npw, &
         &    gxfac, nprojs, czero, vectout, npw)
       else if(choice==2) then
         call abi_zgemm_2r('N', 'N', npw, ndat*nspinor, nprojs, cone, &
         &    dprojs, npw, &
         &    gxfac, nprojs, czero, vectout, npw)
         call abi_zgemm_2r('N', 'N', npw, ndat*nspinor, nprojs, cone, &
         &    projs, npw, &
         &    dgxdtfac, nprojs, cone, vectout, npw)
       else if(choice==3) then
         if(idir<=3) then
           dgxdtfac(:,:,:) = dgxdtfac(:,:,:) - gxfac(:,:,:)
         end if
         call abi_zgemm_2r('N', 'N', npw, ndat*nspinor, nprojs, cone, &
         &    projs, npw, &
         &    dgxdtfac, nprojs, czero, vectout, npw)
         call abi_zgemm_2r('N', 'N', npw, ndat*nspinor, nprojs, cminusone, &
         &    dprojs, npw, &
         &    gxfac, nprojs, cone, vectout, npw)
       else if(choice==5) then
         call abi_zgemm_2r('N', 'N', npw, ndat*nspinor, nprojs, cone, &
         &    projs, npw, &
         &    dgxdtfac, nprojs, czero, vectout, npw)
         call abi_zgemm_2r('N', 'N', npw, ndat*nspinor, nprojs, cone, &
         &    dprojs, npw, &
         &    gxfac, nprojs, cone, vectout, npw)
       else if(choice==51) then
         call abi_zgemm_2r('N', 'N', npw, ndat*nspinor, nprojs, cone, &
         &    projs, npw, &
         &    dgxdtfac, nprojs, czero, vectout, npw)
       end if
     else
       call gemm_nonlop_distributed_gemm_opernlb(rank,nprocs,npw,ndat,nspinor,&
       &                                         nprojs_blk,&
       &                                         nprojs_last_blk,&
       &                                         nprojs_my_blk,cplex,&
       &                                         projs_local,projs_recv,&
       &                                         projs,&
       &                                         gxfac,vectout)
     end if
   else

     if(.not. use_distrib) then
       call DGEMM('N', 'N', npw, ndat*nspinor, nprojs, one, &
       &          projs_r, npw, &
       &          gxfac, nprojs, zero, temp_realvec, npw)
       vectout(1,1:npw*nspinor*ndat) = temp_realvec(1:npw*nspinor*ndat)
       call DGEMM('N', 'N', npw, ndat*nspinor, nprojs, one, &
       &          projs_i, npw, &
       &          gxfac, nprojs, zero, temp_realvec, npw)
       vectout(2,1:npw*nspinor*ndat) = temp_realvec(1:npw*nspinor*ndat)
     else
       call gemm_nonlop_distributed_gemm_opernlb(rank,nprocs,npw,ndat,nspinor,&
       &                                         nprojs_blk,&
       &                                         nprojs_last_blk,&
       &                                         nprojs_my_blk,cplex,&
       &                                         projs_local,projs_recv,&
       &                                         projs_r,&
       &                                         gxfac,temp_realvec)
       vectout(1,1:npw*nspinor*ndat) = temp_realvec(1:npw*nspinor*ndat)
       call gemm_nonlop_distributed_gemm_opernlb(rank,nprocs,npw,ndat,nspinor,&
       &                                         nprojs_blk,&
       &                                         nprojs_last_blk,&
       &                                         nprojs_my_blk,cplex,&
       &                                         projs_local,projs_recv,&
       &                                         projs_i,&
       &                                         gxfac,temp_realvec)
       vectout(2,1:npw*nspinor*ndat) = temp_realvec(1:npw*nspinor*ndat)
     end if

   end if ! cplex_fac == 2

 end if  ! (paw_opt == 0 .or. paw_opt == 1 .or. paw_opt == 4)

end subroutine opernlb_gemm

end module m_opernlb_gemm

