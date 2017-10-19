!{\src2tex{textfont=tt}}
!!****f* ABINIT/fock_ACE_getghc
!! NAME
!!  fock_ACE_getghc
!!
!! FUNCTION
!!  Compute the matrix elements <G|Vx|psi> of the Fock operator in the ACE context.
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2017 ABINIT group (CMartins,FJ,MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cwavef(2,npw*nspinor*ndat)= planewave coefficients of wavefunctions on which Fock operator is applied.
!!  gs_ham <type(gs_hamiltonian_type)>=all data for the Hamiltonian to be applied
!!  mpi_enreg= information about MPI parallelization
!!
!! SIDE EFFECTS
!!  ghc(2,npw*ndat)= matrix elements <G|H|C> or <G|H-lambda.S|C> (if sij_opt>=0 or =-1 in getghc)
!!                   contains the fock exchange term for cwavef at the end.
!!
!! NOTES
!!  The current version assumes that :
!!   * nspinor = 1
!!   * no "my_nspinor"
!!   * no restriction to the value of istwfk_bz (but must be tested in all case)
!!   * all the data for the occupied states (cgocc_bz) are the same as those for the current states (cg)
!!
!! PARENTS
!!      getghc
!!
!! CHILDREN
!!      dotprod_g
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine fock_ACE_getghc(cwavef,ghc,gs_ham,mpi_enreg)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_xmpi
 use m_cgtools, only :dotprod_g
 use m_fock
 use m_hamiltonian, only : gs_hamiltonian_type
 use defs_datatypes, only: pseudopotential_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fock_ACE_getghc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
! Scalars
 type(MPI_type),intent(in) :: mpi_enreg
 type(gs_hamiltonian_type),target,intent(inout) :: gs_ham
! Arrays
 real(dp),intent(inout) :: cwavef(:,:)!,ghc(2,gs_ham%npw_k)
 real(dp),intent(inout) :: ghc(:,:)

!Local variables-------------------------------
! Scalars
 integer :: iband, ier,ikpt,ipw,my_nspinor,nband_k,npw
 real(dp) :: doti,dotr,eigen
 type(fock_common_type),pointer :: fockcommon
! Arrays
 real(dp), allocatable :: ghc1(:,:),xi(:,:)


! *************************************************************************

 ABI_CHECK(associated(gs_ham%fockcommon),"fock must be associated!")
 fockcommon => gs_ham%fockcommon

 ABI_CHECK(gs_ham%nspinor==1,"only allowed for nspinor=1!")
 ABI_CHECK(gs_ham%npw_k==gs_ham%npw_kp,"only allowed for npw_k=npw_kp (ground state)!")

 ikpt=fockcommon%ikpt
 npw=gs_ham%npw_k
 nband_k=fockcommon%nband(ikpt)
 my_nspinor=max(1,gs_ham%nspinor/mpi_enreg%nproc_spinor)
!*Initialization of the array ghc1
!*ghc1 will contain the exact exchange contribution to the Hamiltonian
 ABI_ALLOCATE(ghc1,(2,npw*my_nspinor))
 ghc1=zero
 ABI_ALLOCATE(xi,(2,npw*my_nspinor))

 do iband=1, nband_k
   xi(1,:)=gs_ham%fockACE_k%xi(1,:,iband)
   xi(2,:)=gs_ham%fockACE_k%xi(2,:,iband)

   call dotprod_g(dotr,doti,gs_ham%istwf_k,npw*my_nspinor,2,xi,cwavef,mpi_enreg%me_g0,mpi_enreg%comm_fft)

   ghc1(1,:)=ghc1(1,:)-(dotr*gs_ham%fockACE_k%xi(1,:,iband)-doti*gs_ham%fockACE_k%xi(2,:,iband))
   ghc1(2,:)=ghc1(2,:)-(dotr*gs_ham%fockACE_k%xi(2,:,iband)+doti*gs_ham%fockACE_k%xi(1,:,iband))
 end do
 ABI_DEALLOCATE(xi)

!* If the calculation is parallelized, perform an MPI_allreduce to sum all the contributions in the array ghc
! ghc(:,:)=ghc(:,:)/mpi_enreg%nproc_kpt + ghc1(:,:)
 ghc(:,:)=ghc(:,:) + ghc1(:,:)

! call xmpi_sum(ghc,mpi_enreg%comm_kpt,ier)

! ============================================
! === Calculate the contribution to energy ===
! ============================================
!* Only the contribution when cwavef=cgocc_bz are calculated, in order to cancel exactly the self-interaction
!* at each convergence step. (consistent definition with the defintion of hartree energy)
 if (fockcommon%ieigen/=0) then
   eigen=zero
!* Dot product of cwavef and ghc
!* inspired from the routine 53_spacepar/meanvalue_g but without the reference to parallelism and filtering
   if(gs_ham%istwf_k==2) then
     eigen=half*cwavef(1,1)*ghc1(1,1)
   else
     eigen=cwavef(1,1)*ghc1(1,1)+cwavef(2,1)*ghc1(2,1)
   end if
   do ipw=2,npw
     eigen=eigen+cwavef(1,ipw)*ghc1(1,ipw)+cwavef(2,ipw)*ghc1(2,ipw)
   end do
   if(gs_ham%istwf_k>=2) eigen=two*eigen
!   call xmpi_sum(eigen,mpi_enreg%comm_kpt,ier)
   fockcommon%eigen_ikpt(fockcommon%ieigen)= eigen
   fockcommon%ieigen = 0
 end if

! ===============================
! === Deallocate local arrays ===
! ===============================

 ABI_DEALLOCATE(ghc1)

end subroutine fock_ACE_getghc
!!***
