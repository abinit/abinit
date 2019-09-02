!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ephtk
!! NAME
!!  m_ephtk
!!
!! FUNCTION
!!  Helper functions common to e-ph calculations.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2018 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_ephtk

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_dtset
 use m_crystal
 use m_krank
 use m_xmpi

 use m_fstrings,     only : itoa, sjoin, ltoa, ktoa
 use m_bz_mesh,      only : isamek

 implicit none

 private

 public :: ephtk_set_phmodes_ship     ! Setup a mask to skip accumulating the contribution of certain phonon modes.
 public :: ephtk_set_pertables        ! Set tables for parallelism over perturbations from my_npert and comm
 public :: ephtk_mkqtabs              ! Build tables with correspondence between q-points as needed by complete_gamma.
 public :: ephtk_gam_atm2qnu          ! Compute phonon linewidths from gamma matrix in reduced coordinates.
!!***

contains  !=====================================================
!!***

!!****f* m_ephtk/ephtk_set_phmodes_ship
!! NAME
!!  ephtk_set_phmodes_ship
!!
!! FUNCTION
 ! Setup a mask to skip accumulating the contribution of certain phonon modes.
!!
!! INPUT
!!  dtset<dataset_type>=All input variables for this dataset.
!!
!! OUTPUT
!!   phmodes_skip(natom3) For each mode: 1 to skip this contribution else 0
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ephtk_set_phmodes_ship(dtset, phmodes_skip)

!Arguments ------------------------------------
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,allocatable,intent(out) :: phmodes_skip(:)

!Local variables ------------------------------
!scalars
 integer :: natom3

! *************************************************************************

 ! Setup a mask to skip accumulating the contribution of certain phonon modes.
 ! By default do not skip, if set skip all but specified
 natom3 = dtset%natom * 3
 ABI_MALLOC(phmodes_skip, (natom3))
 phmodes_skip = 0

 if (all(dtset%eph_phrange /= 0)) then
   if (minval(dtset%eph_phrange) < 1 .or. &
       maxval(dtset%eph_phrange) > natom3 .or. &
       dtset%eph_phrange(2) < dtset%eph_phrange(1)) then
     MSG_ERROR('Invalid range for eph_phrange. Should be between [1, 3*natom] and eph_modes(2) > eph_modes(1)')
   end if
   call wrtout(std_out, sjoin(" Including phonon modes between [", &
               itoa(dtset%eph_phrange(1)), ',', itoa(dtset%eph_phrange(2)), "]"))
   phmodes_skip = 0
   phmodes_skip(dtset%eph_phrange(1):dtset%eph_phrange(2)) = 1
 end if

end subroutine ephtk_set_phmodes_ship
!!***

!!****f* m_ephtk/ephtk_set_pertables
!! NAME
!!  ephtk_set_pertables
!!
!! FUNCTION
!!  Build tables for parallelism over perturbations from my_npert and comm
!!
!! INPUT
!!  natom: Number of atoms
!!  my_npert: Number of atomic perturbations or phonon modes treated by this MPI rank.
!!  comm: MPI communicator for parallelism over atomic perturbations.
!!
!! OUTPUT
!!  integer,allocatable :: my_pinfo(:,:)
!!     my_pinfo(3, my_npert)
!!     my_pinfo(1, ip) gives the `idir` index of the ip-th perturbation.
!!     my_pinfo(2, ip) gives the `ipert` index of the ip-th perturbation.
!!     my_pinfo(3, ip) gives `pertcase`=idir + (ipert-1)*3
!!  integer,allocatable :: pert_table(:,:)
!!     pert_table(2, natom3)
!!     pert_table(1, npert): rank of the processor treating this atomic perturbation.
!!     pert_table(2, npert): imyp index in my_pinfo table, -1 if this rank is not treating ipert.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ephtk_set_pertables(natom, my_npert, pert_table, my_pinfo, comm)

!Arguments ------------------------------------
 integer,intent(in) :: natom, my_npert, comm
!arrays
 integer,allocatable :: pert_table(:,:)
 integer,allocatable :: my_pinfo(:,:)

!Local variables ------------------------------
!scalars
 integer :: iatom, idir, pertcase, bstart, bstop, ii, ip, natom3, my_rank, nproc
!arrays
 integer :: all_pinfo(3, natom*3)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

 ! Build table with list of perturbations treated by this CPU.
 natom3 = natom * 3
 ABI_MALLOC(my_pinfo, (3, my_npert))
 ABI_MALLOC(pert_table, (2, natom3))

 do iatom=1,natom
   do idir=1,3
     pertcase = idir + (iatom-1) * 3
     all_pinfo(:, pertcase) = [idir, iatom, pertcase]
     pert_table(1, pertcase) = (pertcase - 1) / (natom3 / nproc)
   end do
 end do
 bstart = (natom3 / nproc) * my_rank + 1
 bstop = bstart + my_npert - 1
 my_pinfo = all_pinfo(:, bstart:bstop)

 pert_table(2, :) = -1
 do ii=1,my_npert
   ip = my_pinfo(3, ii)
   pert_table(2, ip) = ii
 end do
 !write(std_out,*)"my_npert", my_npert, "nproc", nproc; write(std_out,*)"my_pinfo", my_pinfo

end subroutine ephtk_set_pertables
!!***

!!****f* m_ephtk/ephtk_mkqtabs
!! NAME
!!  ephtk_mkqtabs
!!
!! FUNCTION
!!  Build tables with correspondence between q-points in the IBZ/BZ as needed by complete_gamma.
!!
!! INPUT
!!  cryst<crystal_t>=Crystal structure.
!!  nqibz, qibz = Points in the IBZ
!!  nqbz, qbz = Points in the BZ
!!
!! OUTPUT
!! qirredtofull(nqibz) = mapping irred to full qpoints
!! qpttoqpt(2, cryst%nsym, nqbz)) = qpoint index mapping under symops.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ephtk_mkqtabs(cryst, nqibz, qibz, nqbz, qbz, qirredtofull, qpttoqpt)

!Arguments ------------------------------------
 type(crystal_t),intent(in) :: cryst
 integer,intent(in) :: nqibz, nqbz
!arrays
 real(dp),intent(in) :: qibz(3, nqibz), qbz(3, nqbz)
 integer,allocatable :: qirredtofull(:),qpttoqpt(:,:,:)

!Local variables ------------------------------
!scalars
 integer :: iq_bz, iq_ibz, isq_bz, isym
 type(krank_t) :: qrank
!arrays
 integer :: g0(3)
 real(dp) :: qirr(3), tmp_qpt(3)

! *************************************************************************

 qrank = krank_new(nqbz, qbz)

 ! Compute index of IBZ q-point in the BZ array
 ABI_CALLOC(qirredtofull, (nqibz))

 do iq_ibz=1,nqibz
   qirr = qibz(:,iq_ibz)
   iq_bz = qrank%get_index(qirr)
   if (iq_bz /= -1) then
     ABI_CHECK(isamek(qirr, qbz(:,iq_bz), g0), "isamek")
     qirredtofull(iq_ibz) = iq_bz
   else
     MSG_ERROR(sjoin("Full BZ does not contain IBZ q-point:", ktoa(qirr)))
   end if
 end do

 ! Build qpttoqpt table. See also mkqptequiv
 ABI_MALLOC(qpttoqpt, (2, cryst%nsym, nqbz))
 qpttoqpt = -1
 do iq_bz=1,nqbz
   do isym=1,cryst%nsym
     tmp_qpt = matmul(cryst%symrec(:,:,isym), qbz(:,iq_bz))

     isq_bz = qrank%get_index(tmp_qpt)
     if (isq_bz == -1) then
       MSG_ERROR("Looks like no kpoint equiv to q by symmetry without time reversal!")
     end if
     qpttoqpt(1,isym,isq_bz) = iq_bz

     ! q --> -q
     tmp_qpt = -tmp_qpt
     isq_bz = qrank%get_index(tmp_qpt)
     if (isq_bz == -1) then
       MSG_ERROR("Looks like no kpoint equiv to q by symmetry with time reversal!")
     end if
     qpttoqpt(2,isym,isq_bz) = iq_bz
   end do
 end do

 call qrank%free()

end subroutine ephtk_mkqtabs
!!***

!----------------------------------------------------------------------

!!****f* m_ephtk/ephtk_gam_atm2qnu
!! NAME
!! ephtk_gam_atm2qnu
!!
!! FUNCTION
!! This routine takes the gamma matrices in the atom representation and
!! multiplies them by the displ_red matrices. Based on gam_mult_displ
!!
!! INPUTS
!!   natom3 = number of phonon branches (3*natom)
!!   displ_red = phonon mode displacement vectors in reduced coordinates.
!!   gam_bare = bare gamma matrices before multiplication
!!
!! OUTPUT
!!   gam_now = output gamma matrices multiplied by displacement matrices
!!
!! PARENTS
!!
!! CHILDREN
!!      zgemm
!!
!! SOURCE

subroutine ephtk_gam_atm2qnu(natom3, displ_red, gam_atm, gam_qnu)

!Arguments -------------------------------
 integer, intent(in)  :: natom3
 real(dp), intent(in)  :: displ_red(2,natom3,natom3)
 real(dp), intent(in)  :: gam_atm(2,natom3,natom3)
 real(dp), intent(out) :: gam_qnu(natom3)

!Local variables -------------------------
 integer,save :: enough = 0
 integer :: nu
 character(len=500) :: msg
 real(dp) :: zgemm_tmp_mat(2,natom3,natom3), gam_now(2,natom3,natom3)

! *********************************************************************

 call zgemm('c','n',natom3, natom3, natom3, cone, displ_red, natom3, gam_atm, natom3, czero, zgemm_tmp_mat, natom3)

 gam_now = zero
 call zgemm('n','n',natom3,natom3,natom3,cone,zgemm_tmp_mat,natom3,displ_red,natom3,czero,gam_now,natom3)

 ! Extract gamma(q,nu)
 do nu=1,natom3
   gam_qnu(nu) = gam_now(1, nu, nu)
   if (abs(gam_now(2, nu, nu)) > tol8) then
     enough = enough + 1
     if (enough <= 30) then
       write (msg,'(a,i0,a,es16.8)')' non-zero imaginary part for branch: ',nu,', img: ',gam_now(2, nu, nu)
       MSG_WARNING(msg)
     end if
   end if
 end do

end subroutine ephtk_gam_atm2qnu
!!***

end module m_ephtk
!!***
