!!****m* ABINIT/m_ephtk
!! NAME
!!  m_ephtk
!!
!! FUNCTION
!!  Helper functions common to e-ph calculations.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2025 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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
 use m_ebands
 use m_crystal
 use m_krank
 use m_xmpi

 use m_fstrings,     only : itoa, sjoin, ltoa, ftoa, ktoa
 use m_bz_mesh,      only : isamek
 use m_fftcore,      only : get_kg

 implicit none

 private

 public :: ephtk_set_phmodes_skip     ! Setup a mask to skip accumulating the contribution of certain phonon modes.
 public :: ephtk_set_pertables        ! Set tables for parallelism over perturbations from my_npert and comm
 public :: ephtk_mkqtabs              ! Build tables with correspondence between q-points as needed by complete_gamma.
 public :: ephtk_gam_atm2qnu          ! Compute phonon linewidths from gamma matrix in reduced coordinates.
 public :: ephtk_gkknu_from_atm       ! Transform the gkk matrix elements from (atom, red_direction) basis to phonon-mode basis.
 public :: ephtk_update_ebands        ! Update ebands according to dtset%occopt, tsmear, mbpt_sciss, eph_fermie, eph_extrael
 public :: ephtk_get_mpw_gmax
 public :: ephtk_v1atm_to_vqnu        !  Receive potentials in atomic representation and return potential in phonon representation
!!***

 real(dp),public,parameter :: EPHTK_WTOL = tol6
  ! Tolerance for phonon frequencies to be ignored.
  ! Lambda coefficients are set to zero when abs(w) < EPHTK_WTOL
  ! This tolerance is also used in the integrals of a2F(w).

contains  !=====================================================
!!***

!!****f* m_ephtk/ephtk_set_phmodes_skip
!! NAME
!!  ephtk_set_phmodes_skip
!!
!! FUNCTION
!! Setup a mask to skip accumulating the contribution of certain phonon modes.
!!
!! INPUT
!!  eph_phrange=Abinit input variable.
!!
!! OUTPUT
!!   phmodes_skip(natom3) For each mode: 1 to skip the contribution given by this phonon branch else 0
!!
!! SOURCE

subroutine ephtk_set_phmodes_skip(natom, eph_phrange, phmodes_skip)

!Arguments ------------------------------------
 integer,intent(in) :: natom
!arrays
 integer,intent(in) :: eph_phrange(2)
 integer,allocatable,intent(out) :: phmodes_skip(:)

!Local variables ------------------------------
 integer :: natom3
! *************************************************************************

 ! Setup a mask to skip accumulating the contribution of certain phonon modes.
 ! By default do not skip, if set skip all but specified
 natom3 = natom * 3
 ABI_MALLOC(phmodes_skip, (natom3))
 phmodes_skip = 0

 if (all(eph_phrange /= 0)) then
   if (minval(abs(eph_phrange)) < 1 .or. &
       maxval(abs(eph_phrange)) > natom3 .or. &
       abs(eph_phrange(2)) < abs(eph_phrange(1))) then
     ABI_ERROR('Invalid range for eph_phrange. Should be between [1, 3*natom] and eph_modes(2) > eph_modes(1)')
   end if
   if (all(eph_phrange > 0)) then
      call wrtout(std_out, sjoin(" Including phonon modes between [", itoa(eph_phrange(1)), ',', itoa(eph_phrange(2)), "]"))
      phmodes_skip = 1
      phmodes_skip(eph_phrange(1):eph_phrange(2)) = 0
   else if (all(eph_phrange < 0)) then
      call wrtout(std_out, sjoin(" Excluding phonon modes between [", &
                   itoa(abs(eph_phrange(1))), ',', itoa(abs(eph_phrange(2))), "]"))
      phmodes_skip = 0
      phmodes_skip(abs(eph_phrange(1)):abs(eph_phrange(2))) = 1
   else
      ABI_ERROR(sjoin("Invalid eph_phrange: ", itoa(eph_phrange(1)), ',', itoa(eph_phrange(2))))
   end if
 end if

end subroutine ephtk_set_phmodes_skip
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
!! SOURCE

subroutine ephtk_set_pertables(natom, my_npert, pert_table, my_pinfo, comm)

!Arguments ------------------------------------
 integer,intent(in) :: natom, my_npert, comm
!arrays
 integer,allocatable,intent(out) :: pert_table(:,:), my_pinfo(:,:)

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
     ABI_ERROR(sjoin("Full BZ does not contain IBZ q-point:", ktoa(qirr)))
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
       ABI_ERROR("Looks like no kpoint equiv to q by symmetry without time reversal!")
     end if
     qpttoqpt(1,isym,isq_bz) = iq_bz

     ! q --> -q
     tmp_qpt = -tmp_qpt
     isq_bz = qrank%get_index(tmp_qpt)
     if (isq_bz == -1) then
       ABI_ERROR("Looks like no kpoint equiv to q by symmetry with time reversal!")
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
!! This routine takes the gamma matrices in the atomic representation and
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
       ABI_WARNING(msg)
     end if
   end if
 end do

end subroutine ephtk_gam_atm2qnu
!!***

!----------------------------------------------------------------------

!!****f* m_ephtk/ephtk_gkknu_from_atm
!! NAME
!!  ephtk_gkknu_from_atm
!!
!! FUNCTION
!!  Transform the gkk matrix elements from (atom, red_direction) basis to phonon-mode basis.
!!
!! INPUTS
!!  nb1,nb2=Number of bands in gkq_atm matrix.
!!  nk=Number of k-points (usually 1)
!!  natom=Number of atoms.
!!  gkq_atm(2,nb1,nb2,3*natom)=EPH matrix elements in the atomic basis.
!!  phfrq(3*natom)=Phonon frequencies in Ha
!!  displ_red(2,3*natom,3*natom)=Phonon displacement in reduced coordinates.
!!
!! OUTPUT
!!  gkq_nu(2,nb1,nb2,3*natom)=EPH matrix elements in the phonon-mode basis.
!!
!! SOURCE

subroutine ephtk_gkknu_from_atm(nb1, nb2, nk, natom, gkq_atm, phfrq, displ_red, gkq_nu)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nb1, nb2, nk, natom
!arrays
 real(dp),intent(in) :: phfrq(3*natom),displ_red(2,3*natom,3*natom)
 real(dp),intent(in) :: gkq_atm(2,nb1,nb2,nk,3*natom)
 real(dp),intent(out) :: gkq_nu(2,nb1,nb2,nk,3*natom)

!Local variables-------------------------
 integer :: nu,ipc
! *************************************************************************

 gkq_nu = zero

 ! Loop over phonon branches.
 do nu=1,3*natom
   ! Ignore negative or too small frequencies
   if (phfrq(nu) < EPHTK_WTOL) cycle

   ! Transform the gkk from (atom, reduced direction) basis to phonon mode representation
   do ipc=1,3*natom
     gkq_nu(1,:,:,:,nu) = gkq_nu(1,:,:,:,nu) &
       + gkq_atm(1,:,:,:,ipc) * displ_red(1,ipc,nu) &
       - gkq_atm(2,:,:,:,ipc) * displ_red(2,ipc,nu)
     gkq_nu(2,:,:,:,nu) = gkq_nu(2,:,:,:,nu) &
       + gkq_atm(1,:,:,:,ipc) * displ_red(2,ipc,nu) &
       + gkq_atm(2,:,:,:,ipc) * displ_red(1,ipc,nu)
   end do

   gkq_nu(:,:,:,:,nu) = gkq_nu(:,:,:,:,nu) / sqrt(two * phfrq(nu))

   ! Perform the transformation using array operations
   !gkq_nu(1,:,:,:,nu) = sum(gkq_atm(1,:,:,:,:) * displ_red(1,:,nu) - gkq_atm(2,:,:,:,:) * displ_red(2,:,nu), dim=5)
   !gkq_nu(2,:,:,:,nu) = sum(gkq_atm(1,:,:,:,:) * displ_red(2,:,nu) + gkq_atm(2,:,:,:,:) * displ_red(1,:,nu), dim=5)
   !! Apply the normalization factor
   !factor = one / sqrt(two * phfrq(nu))
   !gkq_nu(:,:,:,:,nu) = gkq_nu(:,:,:,:,nu) * factor
 end do

end subroutine ephtk_gkknu_from_atm
!!***

!----------------------------------------------------------------------

!!****f* m_ephtk/ephtk_update_ebands
!! NAME
!!  ephtk_update_ebands
!!
!! FUNCTION
!!  Update ebands according to dtset%occopt, tsmear, mbpt_sciss, eph_fermie, eph_extrael
!!
!! INPUTS
!!  dtset<dataset_type>=All input variables for this dataset.
!!
!! SOURCE

subroutine ephtk_update_ebands(dtset, ebands, header)

!Arguments ------------------------------------
!scalars
 type(dataset_type),intent(in) :: dtset
 type(ebands_t),intent(inout) :: ebands
 character(len=*),intent(in) :: header

!Local variables-------------------------
!scalars
 real(dp),parameter :: nholes = zero
 character(len=500) :: msg
 integer :: unts(2)
! *************************************************************************

 unts = [std_out, ab_out]

 if (dtset%occopt /= ebands%occopt .or. abs(dtset%tsmear - ebands%tsmear) > tol12) then
 !if (.True.) then
   write(msg,"(2a,2(a,i0,a,f14.6,a))")&
   " Changing occupation scheme as input occopt and tsmear differ from those read from WFK file.",ch10,&
   "   From WFK file: occopt = ",ebands%occopt,", tsmear = ",ebands%tsmear,ch10,&
   "   From input:    occopt = ",dtset%occopt,", tsmear = ",dtset%tsmear,ch10
   call wrtout(unts, msg)
   call ebands%set_scheme(dtset%occopt, dtset%tsmear, dtset%spinmagntarget, dtset%prtvol)

   if (abs(dtset%mbpt_sciss) > tol6) then
     ! Apply the scissor operator
     call wrtout(unts, sjoin(" Applying scissors operator to the conduction states with value: ", &
                 ftoa(dtset%mbpt_sciss * Ha_eV, fmt="(f6.2)"), " (eV)"))
     call ebands%apply_scissors(dtset%mbpt_sciss)
   end if
 end if

 ! Default value of eph_fermie is zero hence no tolerance is used!
 if (dtset%eph_fermie /= zero) then
   ABI_CHECK(dtset%eph_extrael == zero, "eph_fermie and eph_extrael are mutually exclusive")
   call wrtout(unts, sjoin(" Fermi level set by the user at:", ftoa(dtset%eph_fermie * Ha_eV, fmt="(f6.2)"), " (eV)"))
   call ebands%set_fermie(dtset%eph_fermie, msg)
   call wrtout(unts, msg)

 else if (abs(dtset%eph_extrael) > zero) then
   call wrtout(unts, sjoin(" Adding eph_extrael:", ftoa(dtset%eph_extrael), "to input nelect:", ftoa(ebands%nelect)))
   call ebands%set_scheme(dtset%occopt, dtset%tsmear, dtset%spinmagntarget, dtset%prtvol, update_occ=.False.)
   call ebands%set_extrael(dtset%eph_extrael, nholes, dtset%spinmagntarget, msg)
   call wrtout(unts, msg)
 end if

 ! Recompute occupations. This is needed if WFK files have been produced in a NSCF run
 ! since occ are set to zero, and fermie is taken from the previous density.
 if (dtset%kptopt > 0) then
   call ebands%update_occ(dtset%spinmagntarget, prtvol=dtset%prtvol)
   call ebands%print([std_out], header=header, prtvol=dtset%prtvol)
 end if

end subroutine ephtk_update_ebands
!!***

!----------------------------------------------------------------------

!!****f* m_ephtk/ephtk_get_mpw_gmax
!! NAME
!!  ephtk_get_mpw_gmax
!!
!! FUNCTION
!! mpw is the maximum number of plane-waves over k and k+q where k and k+q are in the BZ.
!! we also need the max components of the G-spheres (k, k+q) in order to allocate the workspace array work
!! used to symmetrize the wavefunctions in G-space.
!! Note that we loop over the full BZ instead of the IBZ(k)
!! This part is slow for very dense meshes, should try to use a geometrical approach...
!!
!! INPUTS
!!
!! SOURCE

subroutine ephtk_get_mpw_gmax(nkpt, kpts, ecut, gmet, mpw, gmax, comm, init_with_zero)

!Arguments ------------------------------------
 integer,intent(in) :: nkpt
 integer,intent(out) :: mpw, gmax(3)
 real(dp),intent(in) :: ecut, kpts(3,nkpt), gmet(3,3)
 integer,intent(in) :: comm
 logical,optional,intent(in) :: init_with_zero

!Local variables ------------------------------
 integer,parameter :: istwfk1 = 1
 integer :: ik,i1,i2,i3,cnt,ipw,ii,onpw,my_mpw,my_gmax(3),ierr, my_rank, nprocs
 real(dp) :: kk(3), kq(3)
 integer,allocatable :: gtmp(:,:)
 logical :: init_with_zero__
! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 init_with_zero__ = .True.; if (present(init_with_zero)) init_with_zero__ = init_with_zero

 if (init_with_zero__) then
   mpw = 0; gmax = 0
 end if

 cnt = 0
 do ik=1,nkpt
   kk = kpts(:, ik)
   do i3=-1,1
     do i2=-1,1
       do i1=-1,1
         cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle ! MPI parallelism inside comm
         kq = kk + half * [i1, i2, i3]
         ! TODO: g0 umklapp here can enter into play gmax may not be large enough!
         call get_kg(kq, istwfk1, 1.1_dp * ecut, gmet, onpw, gtmp)
         mpw = max(mpw, onpw)
         do ipw=1,onpw
           do ii=1,3
             gmax(ii) = max(gmax(ii), abs(gtmp(ii, ipw)))
           end do
         end do
         ABI_FREE(gtmp)
       end do
     end do
   end do
 end do

 my_mpw = mpw; call xmpi_max(my_mpw, mpw, comm, ierr)
 my_gmax = gmax; call xmpi_max(my_gmax, gmax, comm, ierr)

end subroutine ephtk_get_mpw_gmax
!!***

!!****f* m_epthk/ephtk_v1atm_to_vqnu
!! NAME
!!  ephtk_v1atm_to_vqnu
!!
!! FUNCTION
!!  Receive potentials in atomic representation and return potential in phonon representation
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

pure subroutine ephtk_v1atm_to_vqnu(cplex, nfft, nspden, natom3, v1_atm, displ_red, v1_qnu)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex, nfft, nspden, natom3
!arrays
 real(dp),intent(in) :: v1_atm(cplex, nfft, nspden, natom3)
 real(dp),intent(out) :: v1_qnu(2, nfft, nspden, natom3)
 real(dp),intent(in) :: displ_red(2, natom3, natom3)

!Local variables-------------------------------
!scalars
 integer :: nu, ip, ispden
!************************************************************************

 do nu=1,natom3
   ! v1_qnu = \sum_{ka} phdispl{ka}(q,nu) D_{ka,q} V_scf(r)
   ! NOTE: prefactor 1/sqrt(2 w(q,nu)) is not included in the potentials.
   ! v1_qnu(2, nfft, nspden, natom3), v1_atm(cplex, nfft, nspden, natom3)
   v1_qnu(:, :, :, nu) = zero
   do ip=1,natom3
     do ispden=1,nspden
       if (cplex == 2) then
         v1_qnu(1, :, ispden, nu) = v1_qnu(1, :, ispden, nu) + &
           displ_red(1,ip,nu) * v1_atm(1,:,ispden,ip) - displ_red(2,ip,nu) * v1_atm(2,:,ispden,ip)
         v1_qnu(2, :, ispden, nu) = v1_qnu(2, :, ispden, nu) + &
           displ_red(2,ip,nu) * v1_atm(1,:,ispden,ip) + displ_red(1,ip,nu) * v1_atm(2,:,ispden,ip)
       else
         ! Gamma point. d(q) = d(-q)* --> d is real.
         v1_qnu(1, :, ispden, nu) = v1_qnu(1, :, ispden, nu) + displ_red(1,ip,nu) * v1_atm(1,:,ispden,ip)
       end if
     end do
   end do
 end do

end subroutine ephtk_v1atm_to_vqnu
!!***

end module m_ephtk
!!***
