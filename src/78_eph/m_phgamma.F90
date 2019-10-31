!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_phgamma
!! NAME
!!
!! FUNCTION
!!  Computation of phonon linewidths, isotropic superconducting properties and transport in metals with the LOVA.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group (MG)
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

module m_phgamma

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_krank
 use m_htetra
 use m_tetrahedron
 use libtetrabz
 use m_ifc
 use m_ebands
 use m_fstab
 use iso_c_binding
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_wfk
 use m_ddb
 use m_ddk
 use m_dvdb
 use m_crystal
 use m_fft
 use m_hamiltonian
 use m_pawcprj
 use m_dtset
 use m_dtfil
 use m_wfd
 use m_ephtk

 use defs_abitypes,    only : mpi_type
 use m_time,           only : cwtime, cwtime_report
 use m_fstrings,       only : toupper, itoa, sjoin, ktoa, ltoa, strcat
 use m_numeric_tools,  only : arth, wrap2_pmhalf, simpson_int, simpson, mkherm, get_diag, isdiagmat
 use m_io_tools,       only : open_file, iomode_from_fname
 use m_symtk,          only : littlegroup_q
 use m_geometry,       only : normv
 use m_special_funcs,  only : gaussian
 use m_fftcore,        only : ngfft_seq, get_kg
 use m_cgtools,        only : cg_zdotc
 use m_kg,             only : getph
 use m_dynmat,         only : symdyma, ftgam_init, ftgam, asrif9
 use defs_datatypes,   only : ebands_t, pseudopotential_type
 use m_bz_mesh,        only : kpath_t, kpath_new
 use m_special_funcs,  only : fermi_dirac
 use m_kpts,           only : kpts_ibz_from_kptrlatt, tetra_from_kptrlatt, listkk
 use defs_elphon,      only : complete_gamma !, complete_gamma_tr
 use m_getgh1c,        only : getgh1c, rf_transgrid_and_pack, getgh1c_setup
 use m_pawang,         only : pawang_type
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type
 use m_pawfgr,         only : pawfgr_type

 implicit none

 private
!!***

 public :: eph_phgamma

!----------------------------------------------------------------------

!!****t* m_phgamma/phgamma_t
!! NAME
!! phgamma_t
!!
!! FUNCTION
!! Provides methods for computing phonon linewidths, interpolating the results
!! in q-space and evaluate superconducting properties.
!!
!! SOURCE

 type,public :: phgamma_t

  integer :: natom
  ! Number of atoms per unit cell.

  integer :: natom3
  ! Number of phonon branches i.e. 3*natom.

  integer :: nsppol
  ! Number of independent spin polarizations.

  integer :: nspinor
  ! Number of spinorial components.

  integer :: nqibz
  ! Number of q-points in the IBZ.

  integer :: my_nqibz
  ! Number of q-points from the IBZ for current processor

  integer :: nqbz
  ! Number of q-points in the BZ.

  integer :: eph_scalprod = 0
  ! This to call anaddb routines. Note that eph_scalprod 1 is not supported in eph.

  integer :: bcorr = 0
   ! 1 to include Blochl correction in the tetrahedron method else 0.

  integer :: prteliash = 0
  ! This flaf activates the computation of the Eliashberg function.

  integer :: nrpt
  ! Number of points in the real space representation of the gamma matrices.

  integer :: symgamma
  ! 1 if gamma matrices should be symmetrized by symdyma when using Fourier interpolation

  integer :: asr
  ! If "Acoustic rule at gamma should be enforced.

  integer :: ndir_transp
  ! 0 if no transport, otherwise 3

  integer :: ngqpt(3)
  ! Number of divisions in the Q mesh.

  integer :: nene
  ! Number of chemical potential values used for inelastic integration

  !integer :: my_nqpt
  !integer,allocatable :: my_iqpt(:)

  integer :: my_nfsk_q
  ! Number of k-points in the FS treated by this processor for a given q
  ! Computed in phgamma_setup_qpoint

  integer,allocatable :: my_ifsk_q(:)
  ! Index of the FS k-points treated by this processor for a given q
  ! Computed in phgamma_setup_qpoint

  integer :: my_nspins
   ! Number of spins treated by this MPI rank

  integer,allocatable :: my_spins(:)
   ! my_spins(my_nspins)
   ! Indirect table giving the spin indices treated by this rank.
   ! Used only the collinear case with nspinor == 1

  !integer :: my_npert
   ! Number of atomic perturbations or phonon modes treated by this MPI rank.

  !integer(i1b),allocatable :: itreat_qibz(:)
   ! itreat_qibz(nqibz)
   ! Table used to distribute potentials over q-points in the IBZ.
   ! The loop over qpts in the IBZ(k) is MPI distributed inside qpt_comm accordinging to this table.
   ! 0 if this IBZ point is not treated by this proc.
   ! 1 if this IBZ is treated.

  !integer,allocatable :: my_pinfo(:,:)
   ! my_pinfo(3, my_npert)
   ! my_pinfo(1, ip) gives the `idir` index of the ip-th perturbation.
   ! my_pinfo(2, ip) gives the `ipert` index of the ip-th perturbation.
   ! my_pinfo(3, ip) gives `pertcase`=idir + (ipert-1)*3

  !integer,allocatable :: pert_table(:,:)
   ! pert_table(2, natom3)
   ! pert_table(1, npert): rank of the processor treating this atomic perturbation.
   ! pert_table(2, npert): imyp index in my_pinfo table, -1 if this rank is not treating ipert.

  integer, allocatable :: my_iqibz(:)
  ! indices of ibz iq in local array. -1 if iq does not belong to current proc

  integer, allocatable :: my_iqbz(:)
  ! List of iq_ibz indices treated by current proc

  real(dp) :: enemin
  ! Minimal chemical potential value used for inelastic integration Copied from fstab

  real(dp) :: deltaene
  ! Chemical potential increment for inelastic integration Copied from fstab
  !   for simplicity could be made equal to phonon frequency step

  real(dp) :: gprim(3,3)
  ! Needed for Fourier interpolation.
  ! NOTE: gprim (not gprimd) is used for all FT interpolations,
  ! to be consistent with the dimensions of the rpt, which come from anaddb.

  real(dp),allocatable :: n0(:)
   ! n0(gams%nsppol)=Density of states at the Fermi level.

  real(dp),allocatable :: qibz(:,:)
  ! qibz(3,nqibz)
  ! Reduced coordinates of the q-points in the IBZ.

  real(dp),allocatable :: wtq(:)
  ! wtq(nqibz)
  ! Weights of the q-points in the IBZ (normalized to one)

  real(dp),allocatable :: qbz(:,:)
  ! qbz(3,nqbz)
  ! Reduced coordinates of the q-points in the BZ.

  real(dp),allocatable :: rpt(:,:)
  ! rpt(3,nrpt)
  !  Reduced coordinates ***in terms of rprim*** of the lattice points used
  !  for the Fourier transform of the phonon linewidths.

  real(dp),allocatable :: wghatm(:,:,:)
  ! wghatm(natom,natom,nrpt)
  ! Weights used in the FT of the phonon linewidths.

  real(dp),allocatable :: vals_qibz(:,:,:,:,:)
  ! vals_qibz(2,natom3,natom3,nqibz,nsppol)) in reduced coordinates for each q-point in the IBZ.
  ! matii_qibz {\tau'\alpha',\tau\alpha} = sum over k, ib1, ib2
  !   <psi_{k+q,ib1} | H(1)_{\tau'\alpha'} | psi_{k,ib2}>*  \cdot
  !   <psi_{k+q,ib1} | H(1)_{\tau \alpha } | psi_{k,ib2}>

  !NOTE: choice to put nsppol before or after nqbz is a bit arbitrary
  !   abinit uses nband,nkpt,nsppol, but here for convenience nkpt_phon,nsppol,nqbz interpolation is on qpt
  !MG: I think that nsppol should be the last dimension set call to ftgam.

  real(dp),allocatable :: vals_rpt(:,:,:,:)
  ! vals_rpt(2,natom3**2,nrpt,nsppol)
  ! tgamma matrices in real space in reduced coordinates. Used for the Fourier interpolation.

  ! transport stuff with velocity factors
  real(dp),allocatable :: vals_in_qibz(:,:,:,:,:,:)
  real(dp),allocatable :: vals_out_qibz(:,:,:,:,:,:)
  ! vals_XX_qibz(2,9,natom3,natom3,nqibz,nsppol)) in reduced coordinates for each q-point in the IBZ.

  real(dp),allocatable :: vals_in_rpt(:,:,:,:,:)
  real(dp),allocatable :: vals_out_rpt(:,:,:,:,:)
  ! vals_XX_rpt(2,9,natom3**2,nrpt,nsppol)
  ! tgamma matrices in real space in reduced coordinates.

  ! gamma matrices keeping full electron energy dependency
  real(dp),allocatable :: vals_ee(:,:,:,:,:,:,:)
  ! vals_eew(2, nene, nene, natom3, natom3, nqibz, nsppol)

  contains

   procedure :: free => phgamma_free
    ! Free memory.

   procedure :: interp => phgamma_interp
    ! Interpolates the phonon linewidths.

   procedure :: eval_qibz => phgamma_eval_qibz
    ! Evaluate phonon linewidths without Fourier interpolation.

   procedure :: interp_setup => phgamma_interp_setup
    ! Compute the tables used for the interpolation in q-space.

   procedure :: linwid => phgamma_linwid
    ! Interpolate linewidths along an arbitrary q-path.

 end type phgamma_t

 public :: phgamma_init          ! Creation method.
!!***

!----------------------------------------------------------------------

!!****t* m_phgamma/a2fw_t
!! NAME
!! a2fw_t
!!
!! FUNCTION
!! Store the Eliashberg function a2F(w).
!!
!! SOURCE

 type,public :: a2fw_t

  integer :: nomega
  ! Number of frequency points in a2f(w).

  integer :: nsppol
  ! Number of independent spin polarizations.

  integer :: natom3
  ! Number of phonon modes.

  integer :: nene
  ! Number of chemical potential values used for inelastic integration

  real(dp) :: enemin
  ! Minimal chemical potential value used for inelastic integration Copied from fstab

  real(dp) :: deltaene
  ! Chemical potential increment for inelastic integration Copied from fstab
  !   for simplicity could be made equal to phonon frequency step

  real(dp) :: omega_min,omega_max
  ! min and Max frequency (Ha) in the linear mesh.

  real(dp) :: wstep
  ! Step of the linear mesh

  real(dp) :: smear
  ! Gaussian broadening used to approximated the Dirac distribution.

  integer :: nqshift
  ! Number of shifts in the q-mesh

  integer :: ngqpt(3)
  ! The q-mesh used for calculating vals(w).

  real(dp),allocatable :: qshift(:,:)
  ! qshift(3,nqshift)
  ! The shifts used to generate the q-mesh.

  real(dp),allocatable :: n0(:)
  ! n0(nsppol)
  ! Electronic DOS at the Fermi level.

  real(dp),allocatable :: omega(:)
  ! omega(nomega)
  ! Frequency mesh in Hartree (linear).

  real(dp),allocatable :: vals(:,:,:)
  ! vals(nomega,0:natom3,nsppol)
  ! Eliashberg function
  !   vals(w,1:natom3,1:nsppol): a2f(w) decomposed per phonon branch and spin
  !   vals(w,0,1:nsppol): a2f(w) summed over phonons modes, decomposed in spin

  real(dp),allocatable :: vals_ee(:,:,:,:)
  ! vals_ee(nene,nene,nomega,nsppol)
  ! Eliashberg function
  !   vals(e,e',w,0,1:nsppol): a2f(e,e',w) summed over phonons modes, decomposed in spin

  real(dp),allocatable :: lambdaw(:,:,:)
  ! lambda(nomega,0:natom3,nsppol)

  contains

   procedure :: free => a2fw_free
    ! Free the memory allocated in the structure.

   procedure :: write => a2fw_write
    ! Write alpha^2F(w) to an external file in text/netcdf format

   procedure :: get_moment => a2fw_get_moment
    ! Compute moments of alpha^2F(w)/w .

 end type a2fw_t

 public :: a2fw_init            ! Calculates the FS averaged alpha^2F(w) function.
 !public :: a2fw_solve_gap      ! DEPRECATED
!!***

!!****t* m_phgamma/a2fw_tr_t
!! NAME
!! a2fw_tr_t
!!
!! FUNCTION
!! Store the Eliashberg transport spectral functions:
!!    a2F_trin(w, x, x')
!!    a2F_trout(w, x, x')
!!    a2F_tr(w, x, x') = in - out
!!    a2F_tr_gen(e, e', w, x, x')
!!
!! SOURCE

 type,public :: a2fw_tr_t

  integer :: nomega
  ! Number of frequency points in a2f_tr(w).

  integer :: nene
  ! Number of chemical potential values used for inelastic integration
  ! Number of electron points in a2f_tr_gen(e,e',w).

  integer :: nsppol
  ! Number of independent spin polarizations.

  integer :: natom3
  ! Number of phonon modes.

  real(dp) :: enemin
  ! Minimal chemical potential value used for inelastic integration

  real(dp) :: deltaene
  ! Chemical potential increment for inelastic integration

  real(dp) :: omega_min,omega_max
  ! min and Max frequency (Ha) in the linear mesh.

  real(dp) :: wstep
  ! Step of the linear mesh

  real(dp) :: smear
  ! Gaussian broadening used to approximated the Dirac distribution.

  integer :: nqshift
  ! Number of shifts in the q-mesh

  integer :: ngqpt(3)
  ! The q-mesh used for calculating vals(w).

  real(dp),allocatable :: qshift(:,:)
  ! qshift(3,nqshift)
  ! The shifts used to generate the q-mesh.

  real(dp),allocatable :: n0(:)
  ! n0(nsppol)
  ! Electronic DOS at the Fermi level.

  real(dp),allocatable :: omega(:)
  ! omega(nomega)
  ! Frequency mesh in Hartree (linear).

  real(dp),allocatable :: vals_in(:,:,:,:,:)
  real(dp),allocatable :: vals_out(:,:,:,:,:)
  ! vals_in(nomega,3,3,0:natom3,nsppol)
  ! Eliashberg transport functions for in and out scattering
  !   vals_in(w,3,3,1:natom3,1:nsppol): a2f_tr(w) decomposed per phonon branch and spin
  !   vals_in(w,3,3,0,1:nsppol): a2f_tr(w) summed over phonons modes, decomposed in spin

  real(dp),allocatable :: vals_tr(:,:,:,:,:)
  ! vals_tr(nomega,3,3,0:natom3,nsppol)
  ! transport spectral function = in-out

  real(dp),allocatable :: vals_tr_gen(:,:,:,:,:,:,:)
  ! vals(nene,nene,nomega,3,3,nsppol)
  ! generalized transport spectral function from PB Allen Phys. Rev. Lett. 59, 1460 (1987) [[cite:Allen1987]]

  real(dp),allocatable :: lambdaw_tr(:,:,:,:,:)
  ! lambda(nomega,3,3,0:natom3,nsppol)

  contains

  procedure :: free => a2fw_tr_free
   ! Free the memory allocated in the structure.

  procedure :: write => a2fw_tr_write
   ! Write alpha^2F(w) to an external file in text/netcdf format

 end type a2fw_tr_t

 public :: a2fw_tr_init            ! Calculates the FS averaged alpha^2F_tr,in,out(w, x, x') functions.
!!***

 real(dp),private,parameter :: EPH_Q0TOL = 0.01_dp

contains  !=====================================================
!!***

!----------------------------------------------------------------------

!!****f* m_phgamma/phgamma_free
!! NAME
!! phgamma_free
!!
!! FUNCTION
!!  Free the dynamic memory in a <phgamma_t> datatype
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine phgamma_free(gams)

!Arguments ------------------------------------
 class(phgamma_t),intent(inout) :: gams

! *************************************************************************

 !real
 ABI_SFREE(gams%n0)
 ABI_SFREE(gams%qibz)
 ABI_SFREE(gams%wtq)
 ABI_SFREE(gams%qbz)
 ABI_SFREE(gams%rpt)
 ABI_SFREE(gams%wghatm)
 ABI_SFREE(gams%vals_qibz)
 ABI_SFREE(gams%vals_rpt)
 ABI_SFREE(gams%vals_in_qibz)
 ABI_SFREE(gams%vals_in_rpt)
 ABI_SFREE(gams%vals_out_qibz)
 ABI_SFREE(gams%vals_out_rpt)
 ABI_SFREE(gams%vals_ee)
 ABI_SFREE(gams%my_iqibz)
 ABI_SFREE(gams%my_ifsk_q)
 ABI_SFREE(gams%my_spins)

end subroutine phgamma_free
!!***

!----------------------------------------------------------------------

!!****f* m_phgamma/phgamma_init
!! NAME
!! phgamma_init
!!
!! FUNCTION
!!  Creation method for the phgamma_t datatype.
!!
!! INPUTS
!! cryst<crystal_t>
!! ifc<ifc_type>=Interatomic force constants.
!! symdynmat=1 to activa symmetrization of gamma matrices.
!! ngqpt(3)=Q-mesh divisions
!! n0(dtset%nsppol)=Density of states at the Fermi level.
!!
!! OUTPUT
!! gams<phgamma_t>
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine phgamma_init(gams, cryst, ifc, fstab, dtset, eph_scalprod, ngqpt, n0, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: eph_scalprod
 integer,intent(in) :: comm
 type(crystal_t),intent(in) :: cryst
 type(ifc_type),intent(in) :: ifc
 type(phgamma_t),intent(out) :: gams
 type(fstab_t), intent(in) :: fstab
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: ngqpt(3)
 real(dp),intent(in) :: n0(dtset%nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: qptopt1 = 1
 integer :: my_rank, nproc, ierr, nsppol, natom3
!arrays
 integer :: qptrlatt(3,3)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

 ! Set basic dimensions.
 nsppol = dtset%nsppol; gams%nsppol = nsppol; gams%nspinor = dtset%nspinor

 gams%natom = cryst%natom; gams%natom3 = 3*cryst%natom; natom3 = gams%natom3

 gams%symgamma = dtset%symdynmat; gams%eph_scalprod = eph_scalprod
 !gams%asr = ifc%asr
 gams%asr = 0
 gams%prteliash = dtset%prteliash

 gams%ndir_transp = 0; if (dtset%eph_transport > 0) gams%ndir_transp = 3

 ABI_MALLOC(gams%n0, (nsppol))
 gams%n0 = n0

 gams%nene = fstab%nene
 gams%enemin = fstab%enemin
 gams%deltaene = fstab%deltaene

 ! Setup IBZ, weights and BZ. Always use q --> -q symmetry for phonons even in systems wo inversion
 gams%ngqpt = ngqpt
 qptrlatt = 0; qptrlatt(1,1) = ngqpt(1); qptrlatt(2,2) = ngqpt(2); qptrlatt(3,3) = ngqpt(3)
 call kpts_ibz_from_kptrlatt(cryst, qptrlatt, qptopt1, 1, [zero, zero, zero], &
   gams%nqibz, gams%qibz, gams%wtq, gams%nqbz, gams%qbz)

 ! Allocate matrices in the IBZ.
 ABI_MALLOC_OR_DIE(gams%vals_qibz, (2, natom3, natom3, gams%nqibz, nsppol), ierr)
 gams%vals_qibz = zero

 if (dtset%eph_transport > 0) then
   ABI_MALLOC_OR_DIE(gams%vals_in_qibz, (2, 9, natom3, natom3, gams%nqibz, nsppol), ierr)
   gams%vals_in_qibz = zero
   ABI_MALLOC_OR_DIE(gams%vals_out_qibz, (2, 9, natom3, natom3, gams%nqibz, nsppol), ierr)
   gams%vals_out_qibz = zero
 end if

 if (gams%prteliash == 3) then
   ABI_MALLOC_OR_DIE(gams%vals_ee, (2, gams%nene, gams%nene, natom3, natom3, gams%nqibz, nsppol), ierr)
   gams%vals_ee = zero
 end if

 ! Prepare Fourier interpolation.
 gams%gprim = ifc%gprim
 gams%nrpt  = ifc%nrpt
 ABI_MALLOC(gams%rpt, (3, gams%nrpt))
 gams%rpt = ifc%rpt
 ABI_MALLOC(gams%wghatm, (gams%natom, gams%natom, gams%nrpt))
 gams%wghatm = ifc%wghatm

end subroutine phgamma_init
!!***

!----------------------------------------------------------------------

!!****f* m_phgamma/phgamma_ncwrite
!! NAME
!! phgamma_ncwrite
!!
!! FUNCTION
!!  Write the results stored in the phgamma_t datatype to netcdf file.
!!
!! INPUTS
!! cryst<crystal_t>=Crystalline structure.
!! ifc<ifc_type>=Interatomic force constants.
!! ncid=Netcdf file handler (already open in the caller).
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine phgamma_ncwrite(gams, cryst, ifc, ncid)

!Arguments ------------------------------------
!scalars
 type(phgamma_t),intent(inout) :: gams
 type(crystal_t),intent(in) :: cryst
 type(ifc_type),intent(in) :: ifc
 integer,intent(in) :: ncid

!Local variables-------------------------------
!scalars
 integer :: iq_ibz,spin,mu
 real(dp) :: lambda_tot
 character(len=500) :: msg
!arrays
 real(dp) :: phfrq(3*cryst%natom),gamma_ph(3*cryst%natom),lambda_ph(3*cryst%natom)
 real(dp) :: displ_cart(2,3*cryst%natom,3*cryst%natom)

! *************************************************************************

 ! ==========================================================
 ! write data to files for each q point
 ! ==========================================================
 ! Compute total lambda
 lambda_tot = zero
 do spin=1,gams%nsppol
   do iq_ibz=1,gams%nqibz

     ! Get phonon frequencies and eigenvectors.
     call gams%eval_qibz(cryst, ifc, iq_ibz, spin, phfrq, gamma_ph, lambda_ph, displ_cart)

     do mu=1,gams%natom3
       lambda_tot = lambda_tot + lambda_ph(mu) * gams%wtq(iq_ibz)
     end do

#ifdef HAVE_NETCDF
   ! Write data to netcdf file
   if (ncid /= nctk_noid) then
     if (spin == 1) then
       NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "phfreq_qibz"), phfrq, start=[1, iq_ibz]))
       NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'phdispl_cart_qibz'), displ_cart, start=[1, 1, 1, iq_ibz]))
     end if
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'phgamma_qibz'), gamma_ph, start=[1, iq_ibz, spin]))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'phlambda_qibz'), lambda_ph, start=[1, iq_ibz, spin]))
   end if
#endif

     ! Output to the main output file
     if (gams%nsppol == 2) then
       write(msg,'(2a,3es16.6,a,i1,a,a)')ch10,&
         ' q-point =',gams%qibz(:, iq_ibz),'   spin = ',spin,ch10,&
         ' Mode number    Frequency (Ha)  Linewidth (Ha)  Lambda(q,n)'
     else
       write(msg,'(2a,3es16.6,a,a)')ch10,&
         ' q-point =',gams%qibz(:, iq_ibz),ch10,&
         ' Mode number    Frequency (Ha)  Linewidth (Ha)  Lambda(q,n)'
     end if
     call wrtout([std_out, ab_out], msg)

     do mu=1,gams%natom3
       write(msg,'(i5,es20.6,2es16.6)')mu, phfrq(mu), gamma_ph(mu), lambda_ph(mu)
       call wrtout([std_out, ab_out], msg)
     end do

   end do
   ! Add blank lines to output files between spins
   call wrtout([std_out, ab_out], "")
 end do

 write(ab_out,"(a,f8.4)")" lambda= ",lambda_tot

end subroutine phgamma_ncwrite
!!***

!----------------------------------------------------------------------

!!****f* m_phgamma/tgamma_symm
!! NAME
!! tgamma_symm
!!
!! FUNCTION
!!  Symmetrize the tgamma matrix
!!
!! INPUTS
!! qpt(3)=phonon wavevector in reduced coordinates.
!! cryst<crystal_t>=Crystalline structure.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine tgamma_symm(cryst, qpt, tgamma)

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: cryst
!arrays
 real(dp),intent(in) :: qpt(3)
 real(dp),intent(inout) :: tgamma(2,3*cryst%natom,3*cryst%natom)

!Local variables-------------------------------
!scalars
 integer :: ii,natom3,k
!arrays
 real(dp) :: tgcart(2,3*cryst%natom,3*cryst%natom)
 real(dp) :: umat(2,3*cryst%natom,3*cryst%natom),tmp_mat(2,3*cryst%natom,3*cryst%natom)

! *********************************************************************

 ! Build U matrix.
 umat = zero; k = 1
 do ii=1,cryst%natom
   umat(1,k:k+2, k:k+2) = cryst%gprimd
   k = k + 3
 end do

 natom3 = 3 * cryst%natom

 ! Reduced --> Cartesian
 call zgemm('N','N', natom3, natom3, natom3, cone, tgamma, natom3, umat, natom3, czero, tmp_mat, natom3)
 call zgemm('T','N', natom3, natom3, natom3, cone, umat, natom3, tmp_mat, natom3, czero, tgcart, natom3)

 ! Make the matrix hermitian
 call mkherm(tgcart, 3*cryst%natom)

 ! Symmetrize tgamma matrix.
 call symdyma(tgcart, cryst%indsym, cryst%natom, cryst%nsym, qpt, cryst%rprimd, cryst%symrel, cryst%symafm)

 umat = zero; k = 1
 do ii=0,cryst%natom-1
   umat(1,k:k+2, k:k+2) = cryst%rprimd
   k = k + 3
 end do

 ! Cartesian --> Reduced
 call zgemm('N','N', natom3, natom3, natom3, cone, tgcart, natom3, umat, natom3, czero, tmp_mat, natom3)
 call zgemm('T','N', natom3, natom3, natom3, cone, umat, natom3, tmp_mat, natom3, czero, tgamma, natom3)

end subroutine tgamma_symm
!!***

!----------------------------------------------------------------------

!!****f* m_phgamma/phgamma_eval_qibz
!! NAME
!! phgamma_eval_qibz
!!
!! FUNCTION
!! Compute the phonon linewidths for q-points in the IBZ without performing interpolation.
!!
!! INPUTS
!!  gams<phgamma_t>
!!  cryst<crystal_t>=Crystal structure.
!!  ifc<ifc_type>=Interatomic force constants.
!!  iq_ibz=Index of the q-point in the IBZ array.
!!  spin=Spin index
!!
!! OUTPUT
!!  phfrq(gams%natom3)=Phonon frequencies
!!  gamma_ph(gams%natom3)=Phonon linewidths.
!!  lambda_ph(gams%natom3)=Phonon linewidths.
!!  displ_cart(2,3,cry%natom,3*cryst%natom)=Phonon displacement in cartesian coordinates.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine phgamma_eval_qibz(gams, cryst, ifc, iq_ibz, spin, phfrq, gamma_ph, lambda_ph, displ_cart, gamma_ph_ee)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq_ibz,spin
 class(phgamma_t),intent(inout) :: gams
 type(crystal_t),intent(in) :: cryst
 type(ifc_type),intent(in) :: ifc
!arrays
 real(dp),intent(out) :: phfrq(gams%natom3),gamma_ph(gams%natom3),lambda_ph(gams%natom3)
 real(dp),intent(out),optional :: gamma_ph_ee(gams%nene,gams%nene,gams%natom3)
 real(dp),intent(out) :: displ_cart(2,3,cryst%natom,3*cryst%natom)

!Local variables-------------------------------
!scalars
 integer,parameter :: qtor0 = 0
 integer :: natom3, nu1, iene, jene
 real(dp) :: spinfact
 !character(len=500) :: msg
 !arrays
 real(dp) :: displ_red(2,gams%natom3,gams%natom3)
 real(dp) :: work_qnu(gams%natom3)
 real(dp) :: gam_atm(2,gams%natom3,gams%natom3)

! *************************************************************************

 !@phgamma_t
 natom3 = gams%natom3

 ! Get phonon frequencies and eigenvectors.
 call ifc%fourq(cryst, gams%qibz(:,iq_ibz), phfrq, displ_cart, out_displ_red=displ_red)

 ! Scalar product with the displ_red vectors.
 gam_atm = reshape(gams%vals_qibz(:,:,:,iq_ibz,spin), [2, natom3, natom3])
 call ephtk_gam_atm2qnu(natom3, displ_red, gam_atm, gamma_ph)

 if (present (gamma_ph_ee)) then
   do iene = 1, gams%nene
     do jene = 1, gams%nene
       gam_atm = reshape(gams%vals_ee(:,jene,iene,:,:,iq_ibz,spin), [2, natom3, natom3])
       call ephtk_gam_atm2qnu(natom3, displ_red, gam_atm, work_qnu)
       gamma_ph_ee(jene, iene, :) = work_qnu
     end do
   end do
 end if

 ! Compute lambda
 ! TODO: check this - looks like a factor of 2 wrt the inline documentation!
 ! NB: one factor of 2 comes from the phonon propagator and BE factor,
 ! then you have to be careful with the convention for the Fermi level DOS
 !
 ! spinfact should be 1 for a normal non sppol calculation without spinorbit
 ! for spinors it should also be 1 as bands are twice as numerous but n0 has been divided by 2
 ! for nsppol 2 it should be 0.5 as we have 2 spin channels to sum
 spinfact = two / (gams%nsppol * gams%nspinor)

 do nu1=1,gams%natom3
   gamma_ph(nu1) =  gamma_ph(nu1) * pi * spinfact
   lambda_ph(nu1) = zero
   if (abs(phfrq(nu1)) > EPHTK_WTOL) lambda_ph(nu1) = gamma_ph(nu1) / (two * pi * gams%n0(spin) * phfrq(nu1)**2)
   if (present(gamma_ph_ee)) gamma_ph_ee(:,:,nu1) =  gamma_ph_ee(:,:,nu1) * pi * spinfact
 end do

 ! This to avoid spurious results for the acoustic modes.
 ! In principle, gamma(q) --> 0 and lambda(q) --> 0 for q --> 0
 ! but the Fourier interpolated gammas do not fulfill this property so we set everything
 ! to zero when we are inside a sphere or radius
 if (normv(gams%qibz(:,iq_ibz), cryst%gmet, "G") < EPH_Q0TOL) then
   gamma_ph(1:3) = zero
   lambda_ph(1:3) = zero
 end if

end subroutine phgamma_eval_qibz
!!***

!----------------------------------------------------------------------

!!****f* m_phgamma/phgamma_interp
!! NAME
!! phgamma_interp
!!
!! FUNCTION
!!  Interpolate the phonon linewidths at a given q-point.
!!
!! INPUTS
!!  gams<phgamma_t>
!!  cryst<crystal_t>=crystalline structure.
!!  ifc<ifc_type>=Interatomic force constants.
!!  spin=Spin index
!!  qpt(3)=q-point in reduced coordinates
!!  gamma_ph(3*natom)=Phonon linewidths
!!
!! OUTPUT
!!  gamma_ph(gams%natom3)=Interpolated Phonon linewidths.
!!  lamda_ph(3*natom)=Lambda coefficients for the different phonon modes.
!!  phfrq(3*natom)=phonon frequencies at current q
!!  displ_cart(2,3,natom,3*natom) = Phonon displacement in Cartesian coordinates
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine phgamma_interp(gams, cryst, ifc, spin, qpt, phfrq, gamma_ph, lambda_ph, displ_cart, gamma_ph_ee)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spin
 class(phgamma_t),intent(inout) :: gams
 type(crystal_t),intent(in) :: cryst
 type(ifc_type),intent(in) :: ifc
!arrays
 real(dp),intent(in) :: qpt(3)
 real(dp),intent(out) :: phfrq(gams%natom3),gamma_ph(gams%natom3),lambda_ph(gams%natom3)
 real(dp),intent(out),optional :: gamma_ph_ee(gams%nene,gams%nene,gams%natom3)
 real(dp),intent(out) :: displ_cart(2,3,cryst%natom,3*cryst%natom)

!Local variables-------------------------------
!scalars
 integer,parameter :: qtor0 = 0
 integer, save :: icall=0
 integer :: natom3,nu1
 real(dp) :: spinfact
 character(len=500) :: msg
 !arrays
 real(dp) :: displ_red(2,gams%natom3,gams%natom3)
 real(dp) :: gam_now(2,gams%natom3**2) !,work_qnu(gams%natom3)
 real(dp) :: gam_atm(2,gams%natom3,gams%natom3)
 real(dp),allocatable :: coskr(:,:),sinkr(:,:)

! *************************************************************************

 ! Compute internal tables used for Fourier interpolation.
 if (.not. allocated(gams%vals_rpt)) call gams%interp_setup(cryst)

 if (present(gamma_ph_ee) .and. icall == 0) then
   gamma_ph_ee = zero
   write (msg,'(3a)')&
     " For the moment gams_ee matrix elements are not FT interpolated wrt q,",ch10,&
     " only evaluated on the electron k grid. The resulting a2feew will be 0"
   MSG_WARNING(msg)
   icall = 1
 end if

 natom3 = gams%natom3

 ! Taken from mkph_linwid
 ! This reduced version of ftgkk supposes the kpoints have been integrated
 ! in integrate_gamma. Do FT from real-space gamma grid to 1 qpt.
 ABI_MALLOC(coskr, (1, gams%nrpt))
 ABI_MALLOC(sinkr, (1, gams%nrpt))
 ! TODO: This is not optimal
 call ftgam_init(gams%gprim, 1, gams%nrpt, qpt, gams%rpt, coskr, sinkr)

 call ftgam(gams%wghatm, gam_now, gams%vals_rpt(:,:,:,spin), gams%natom, 1, gams%nrpt, qtor0, coskr, sinkr)

 ! This call is not executed in elphon!
 if (gams%symgamma == 1) call tgamma_symm(cryst, qpt, gam_now)

 ABI_FREE(coskr)
 ABI_FREE(sinkr)

 ! Get phonon frequencies and eigenvectors.
 call ifc%fourq(cryst, qpt, phfrq, displ_cart, out_displ_red=displ_red)

 ! Scalar product with the displ_red
 gam_atm = reshape(gam_now, [2, natom3, natom3])
 call ephtk_gam_atm2qnu(natom3, displ_red, gam_atm, gamma_ph)

 ! Compute lambda
 ! spinfact should be 1 for a normal non sppol calculation without spinorbit
 ! for spinors it should also be 1 as bands are twice as numerous but n0 has been divided by 2
 ! for nsppol 2 it should be 0.5 as we have 2 spin channels to sum
 spinfact = two / (gams%nsppol * gams%nspinor)

 ! Compute lambda
 do nu1=1,gams%natom3
   gamma_ph(nu1) = gamma_ph(nu1) * pi * spinfact
   lambda_ph(nu1) = zero
   if (abs(phfrq(nu1)) > EPHTK_WTOL) lambda_ph(nu1) = gamma_ph(nu1) / (two * pi * gams%n0(spin) * phfrq(nu1)**2)
 end do

 ! This to avoid spurious results for the acoustic modes.
 ! In principle, gamma(q) --> 0 and lambda(q) --> 0 for q --> 0
 ! but the Fourier interpolated gammas do not fulfill this property so we set everything
 ! to zero when we are inside a sphere or radius
 if (normv(qpt, cryst%gmet, "G") < EPH_Q0TOL) then
   write(std_out,*)"Setting values to zero."
   gamma_ph(1:3) = zero
   lambda_ph(1:3) = zero
 end if

end subroutine phgamma_interp
!!***

!----------------------------------------------------------------------

!!****f* m_phgamma/phgamma_interp_setup
!! NAME
!! phgamma_interp_setup
!!
!! FUNCTION
!!  This routines prepares the internal tables used to interpolate the linewidths in q-space
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine phgamma_interp_setup(gams, cryst)

!Arguments ------------------------------------
!scalars
 class(phgamma_t),intent(inout) :: gams
 type(crystal_t),intent(in) :: cryst

!Local variables-------------------------------
!scalars
 integer,parameter :: qtor1 = 1
 integer :: iq_bz,iq_ibz,spin,ierr,ii
 !character(len=500) :: msg
!arrays
 integer,allocatable :: qirredtofull(:),qpttoqpt(:,:,:)
 real(dp),allocatable :: coskr(:,:),sinkr(:,:)
 real(dp),allocatable :: gamma_qpt(:,:,:,:),atmfrc(:,:)
 real(dp),allocatable :: vals_bz(:,:,:,:)

! *************************************************************************

 ABI_MALLOC_OR_DIE(vals_bz, (2, gams%natom3**2, gams%nqbz, gams%nsppol), ierr)
 vals_bz = zero

 ! Build tables needed by complete_gamma.
 call ephtk_mkqtabs(cryst, gams%nqibz, gams%qibz, gams%nqbz, gams%qbz, qirredtofull, qpttoqpt)

 ! Fill BZ array with IBZ data.
 do spin=1,gams%nsppol
   do iq_ibz=1,gams%nqibz
     iq_bz = qirredtofull(iq_ibz)
     vals_bz(:,:,iq_bz,spin) = reshape(gams%vals_qibz(:,:,:,iq_ibz,spin), [2, gams%natom3**2])
   end do
 end do

 ! Complete vals_bz in the full BZ.
 ! FIXME: Change complete_gamma API to pass (..., nsppol)
 ABI_MALLOC(gamma_qpt, (2, gams%natom3**2, gams%nsppol, gams%nqbz))
 do spin=1,gams%nsppol
   gamma_qpt(:, :, spin, :) = vals_bz(:, :, :, spin)
 end do

 call complete_gamma(cryst, gams%natom3, gams%nsppol, gams%nqibz, gams%nqbz, &
   gams%eph_scalprod, qirredtofull, qpttoqpt, gamma_qpt)

 do spin=1,gams%nsppol
   vals_bz(:, :, :, spin) = gamma_qpt(:, :, spin, :)
 end do

 ABI_FREE(gamma_qpt)
 ABI_FREE(qirredtofull)
 ABI_FREE(qpttoqpt)

 ! This call is not executed in elphon!
 if (gams%symgamma == 1) then
   do spin=1,gams%nsppol
     do iq_bz=1,gams%nqbz
       call tgamma_symm(cryst, gams%qbz(:,iq_bz), vals_bz(:,:,iq_bz,spin))
     end do
   end do
 end if

 ! Now FT to real space too
 ! NOTE: gprim (not gprimd) is used for all FT interpolations,
 ! to be consistent with the dimensions of the rpt, which come from anaddb.
 ABI_MALLOC_OR_DIE(gams%vals_rpt, (2, gams%natom3**2, gams%nrpt, gams%nsppol), ierr)
 gams%vals_rpt = zero

 ! q --> r
 ABI_MALLOC(coskr, (gams%nqbz, gams%nrpt))
 ABI_MALLOC(sinkr, (gams%nqbz, gams%nrpt))
 call ftgam_init(gams%gprim, gams%nqbz, gams%nrpt, gams%qbz, gams%rpt, coskr, sinkr)

 do spin=1,gams%nsppol
   call ftgam(gams%wghatm, vals_bz(:,:,:,spin), gams%vals_rpt(:,:,:,spin), gams%natom, gams%nqbz,&
      gams%nrpt, qtor1, coskr, sinkr)

   ! Enforce "acoustic" rule on vals_rpt
   ! This call is not executed in elphon!
   if (gams%asr /= 0) then
     ABI_MALLOC(atmfrc, (3*gams%natom*3*gams%natom, gams%nrpt))
     do ii=1,2
       atmfrc = gams%vals_rpt(ii,:,:,spin)
       !gals%vals_rpt(2,:,,spin) = zero
       call asrif9(gams%asr, atmfrc, gams%natom, gams%nrpt, gams%rpt, gams%wghatm)
       gams%vals_rpt(ii,:,:,spin) = atmfrc
     end do
     ABI_FREE(atmfrc)
   end if
 end do

 ABI_FREE(vals_bz)
 ABI_FREE(coskr)
 ABI_FREE(sinkr)

end subroutine phgamma_interp_setup
!!***

!----------------------------------------------------------------------

!!****f* m_phgamma/phgamma_vv_eval_qibz
!! NAME
!! phgamma_vv_eval_qibz
!!
!! FUNCTION
!! Compute the phonon linewidths times velocity squared, for q-points in the IBZ without performing interpolation.
!!
!! INPUTS
!!  gams<phgamma_t>
!!  cryst<crystal_t>=Crystal structure.
!!  ifc<ifc_type>=Interatomic force constants.
!!  iq_ibz=Index of the q-point in the IBZ array.
!!  spin=Spin index
!!
!! OUTPUT
!!  phfrq(gams%natom3)=Phonon frequencies
!!  gamma_in_ph(gams%natom3)=Phonon linewidths.
!!  gamma_out_ph(gams%natom3)=Phonon linewidths.
!!  lambda_in_ph(gams%natom3)=Phonon linewidths.
!!  lambda_out_ph(gams%natom3)=Phonon linewidths.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine phgamma_vv_eval_qibz(gams, cryst, ifc, iq_ibz, spin, phfrq, gamma_in_ph, gamma_out_ph, lambda_in_ph, lambda_out_ph)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq_ibz,spin
 type(phgamma_t),intent(inout) :: gams
 type(crystal_t),intent(in) :: cryst
 type(ifc_type),intent(in) :: ifc
!arrays
 real(dp),intent(out) :: phfrq(gams%natom3)
 real(dp),intent(out) :: gamma_in_ph(3,3,gams%natom3),lambda_in_ph(3,3,gams%natom3)
 real(dp),intent(out) :: gamma_out_ph(3,3,gams%natom3),lambda_out_ph(3,3,gams%natom3)

!Local variables-------------------------------
!scalars
 integer,parameter :: qtor0 = 0
 integer :: natom3, nu1, idir, jdir, ii
 real(dp) :: spinfact
 !character(len=500) :: msg
!arrays
 real(dp) :: displ_cart(2,3,cryst%natom,3*cryst%natom), displ_red(2,gams%natom3,gams%natom3)
 real(dp) :: work_qnu(gams%natom3), gam_atm(2,gams%natom3,gams%natom3)

! *************************************************************************

 !@phgamma_t
 natom3 = gams%natom3

 ! Get phonon frequencies and eigenvectors.
 call ifc%fourq(cryst, gams%qibz(:,iq_ibz), phfrq, displ_cart, out_displ_red=displ_red)

 do jdir=1,3
   do idir=1,3
     ii = idir + gams%ndir_transp * (jdir - 1)
     ! Scalar product with the displ_red vectors.
     gam_atm = reshape(gams%vals_in_qibz(:,ii,:,:,iq_ibz,spin), [2, natom3, natom3])
     call ephtk_gam_atm2qnu(natom3, displ_red, gam_atm, work_qnu)
     gamma_in_ph(idir,jdir,:) = work_qnu

     gam_atm = reshape(gams%vals_out_qibz(:,ii,:,:,iq_ibz,spin), [2, natom3, natom3])
     call ephtk_gam_atm2qnu(natom3, displ_red, gam_atm, work_qnu)
     gamma_out_ph(idir,jdir,:) = work_qnu
   end do ! idir
 end do ! jdir

 ! Compute lambda
 do jdir =1,3
   do idir =1,3
     ! TODO : check this - looks like a factor of 2 wrt the inline documentation!
     !spinfact should be 1 for a normal non sppol calculation without spinorbit
     !for spinors it should also be 1 as bands are twice as numerous but n0 has been divided by 2
     !for sppol 2 it should be 0.5 as we have 2 spin channels to sum
     spinfact = two / (gams%nsppol*gams%nspinor)

     do nu1=1,gams%natom3
       gamma_in_ph(idir,jdir,nu1)  =  gamma_in_ph(idir,jdir,nu1)  * pi * spinfact
       gamma_out_ph(idir,jdir,nu1) =  gamma_out_ph(idir,jdir,nu1) * pi * spinfact
       lambda_in_ph(idir,jdir,nu1)  = zero
       lambda_out_ph(idir,jdir,nu1) = zero
       if (abs(phfrq(nu1)) > EPHTK_WTOL) then
         lambda_in_ph(idir,jdir,nu1)  = gamma_in_ph(idir,jdir,nu1)  / (two * pi * gams%n0(spin) * phfrq(nu1)**2)
         lambda_out_ph(idir,jdir,nu1) = gamma_out_ph(idir,jdir,nu1) / (two * pi * gams%n0(spin) * phfrq(nu1)**2)
       end if
     end do
   end do ! idir
 end do ! jdir

end subroutine phgamma_vv_eval_qibz
!!***

!----------------------------------------------------------------------

!!****f* m_phgamma/phgamma_vv_interp
!! NAME
!! phgamma_interp
!!
!! FUNCTION
!!  Interpolate the linewidths at a given q-point.
!!
!! INPUTS
!!  gams<phgamma_t>
!!  cryst<crystal_t>=crystalline structure.
!!  ifc<ifc_type>=Interatomic force constants.
!!  spin=Spin index
!!  qpt(3)=q-point in reduced coordinates
!!  gamma_ph(3*natom)=Phonon linewidths
!!  displ_cart(2,3,cryst%natom,3*cryst%natom)=Phonon displacement in cartesian coordinates.
!!
!! OUTPUT
!!  gamma_ph(gams%natom3)=Interpolated Phonon linewidths.
!!  lamda_ph(3*natom)=Lambda coefficients for the different phonon modes.
!!  phfrq(3*natom)=phonon frequencies at current q
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine phgamma_vv_interp(gams, cryst, ifc, spin, qpt, phfrq, gamma_in_ph, gamma_out_ph, lambda_in_ph, lambda_out_ph)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spin
 type(phgamma_t),intent(inout) :: gams
 type(crystal_t),intent(in) :: cryst
 type(ifc_type),intent(in) :: ifc
!arrays
 real(dp),intent(in) :: qpt(3)
 real(dp),intent(out) :: phfrq(gams%natom3)
 real(dp),intent(out) :: gamma_in_ph(3,3,gams%natom3),lambda_in_ph(3,3,gams%natom3)
 real(dp),intent(out) :: gamma_out_ph(3,3,gams%natom3),lambda_out_ph(3,3,gams%natom3)

!Local variables-------------------------------
!scalars
 integer,parameter :: qtor0 = 0
 integer :: natom3, nu1, idir,jdir,ii
 real(dp) :: spinfact
 !character(len=500) :: msg
 !arrays
 real(dp) :: displ_cart(2,3,cryst%natom,3*cryst%natom)
 real(dp) :: displ_red(2,gams%natom3,gams%natom3),work_qnu(gams%natom3)
 real(dp) :: gam_in_now(2,3,3,gams%natom3**2)
 real(dp) :: gam_out_now(2,3,3,gams%natom3**2)
 real(dp) :: gam_atm(2,gams%natom3,gams%natom3)
 real(dp),allocatable :: coskr(:,:),sinkr(:,:)

! *************************************************************************

 ! Compute internal tables used for Fourier interpolation.
 if (.not.allocated(gams%vals_in_rpt)) call phgamma_vv_interp_setup(gams, cryst)

 natom3 = gams%natom3

 ! Get phonon frequencies and eigenvectors.
 call ifc%fourq(cryst, qpt, phfrq, displ_cart, out_displ_red=displ_red)

 ! Taken from mkph_linwid
 ! This reduced version of ftgkk supposes the kpoints have been integrated
 ! in integrate_gamma. Do FT from real-space gamma grid to 1 qpt.
 ABI_MALLOC(coskr, (1,gams%nrpt))
 ABI_MALLOC(sinkr, (1,gams%nrpt))
 ! TODO: This is not optimal
 call ftgam_init(gams%gprim, 1, gams%nrpt, qpt, gams%rpt, coskr, sinkr)

 do idir=1,3
   do jdir=1,3
     ii = idir+gams%ndir_transp*(jdir-1)

     call ftgam(gams%wghatm, gam_in_now, gams%vals_in_rpt(:,ii,:,:,spin), gams%natom, 1, gams%nrpt, qtor0, coskr, sinkr)
     call ftgam(gams%wghatm, gam_out_now, gams%vals_out_rpt(:,ii,:,:,spin), gams%natom, 1, gams%nrpt, qtor0, coskr, sinkr)

     ! This call is not executed in elphon!
     ! TODO: needs to take into account the matrix nature of _in_ and _out_ gammas
     if (gams%symgamma == 1) then
       call tgamma_symm(cryst, qpt, gam_in_now)
       call tgamma_symm(cryst, qpt, gam_out_now)
     end if

     ! Scalar product with the displ_red
     gam_atm = reshape(gam_in_now, [2, natom3, natom3])
     call ephtk_gam_atm2qnu(natom3, displ_red, gam_atm, work_qnu)
     gamma_in_ph(idir,jdir,:) = work_qnu

     gam_atm = reshape(gam_out_now, [2, natom3, natom3])
     call ephtk_gam_atm2qnu(natom3, displ_red, gam_atm, work_qnu)
     gamma_out_ph(idir,jdir,:) = work_qnu
   end do
 end do

 ABI_FREE(coskr)
 ABI_FREE(sinkr)

 ! Compute lambda
 ! spinfact should be 1 for a normal non sppol calculation without spinorbit
 ! for spinors it should also be 1 as bands are twice as numerous but n0 has been divided by 2
 ! for sppol 2 it should be 0.5 as we have 2 spin channels to sum
 spinfact = two / (gams%nsppol * gams%nspinor)

 ! Compute lambda
 do nu1=1,gams%natom3
   do idir=1,3
     do jdir=1,3
       gamma_in_ph(idir,jdir,nu1) = gamma_in_ph(idir,jdir,nu1) * pi * spinfact
       gamma_out_ph(idir,jdir,nu1) = gamma_out_ph(idir,jdir,nu1) * pi * spinfact
       lambda_in_ph(idir,jdir,nu1) = zero
       lambda_out_ph(idir,jdir,nu1) = zero
       if (abs(phfrq(nu1)) > EPHTK_WTOL) then
         lambda_in_ph(idir,jdir,nu1)  = gamma_in_ph(idir,jdir,nu1) / (two * pi * gams%n0(spin) * phfrq(nu1)**2)
         lambda_out_ph(idir,jdir,nu1) = gamma_out_ph(idir,jdir,nu1) / (two * pi * gams%n0(spin) * phfrq(nu1)**2)
       end if
     end do
   end do
 end do

end subroutine phgamma_vv_interp
!!***

!----------------------------------------------------------------------

!!****f* m_phgamma/phgamma_vv_interp_setup
!! NAME
!! phgamma_vv_interp_setup
!!
!! FUNCTION
!!  This routines preparse the internal tables used to interpolate the vv_linewidths in q-space
!!
!! INPUTS
!!  action =
!!    "INIT" to allocate and compute the internal tables (default)
!!    "FREE" to deallocate the internal tables.
!!
!! SIDE EFFECTS
!!  gams<phgamma_t>= gams%vals_in_rpt, etc... depending on action.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine phgamma_vv_interp_setup(gams, cryst)

!Arguments ------------------------------------
!scalars
 type(phgamma_t),intent(inout) :: gams
 type(crystal_t),intent(in) :: cryst

!Local variables-------------------------------
!scalars
 integer,parameter :: qtor1 = 1
 integer :: iq_bz,iq_ibz,spin,ierr, ii, idir, jdir
 !character(len=500) :: msg
!arrays
 integer,allocatable :: qirredtofull(:),qpttoqpt(:,:,:)
 real(dp),allocatable :: coskr(:,:),sinkr(:,:)
 real(dp),allocatable :: vals_in_bz(:,:,:,:,:)
 real(dp),allocatable :: vals_out_bz(:,:,:,:,:)

! *************************************************************************

 ABI_MALLOC_OR_DIE(vals_in_bz,(2, 9, gams%natom3**2, gams%nqbz, gams%nsppol), ierr)
 ABI_MALLOC_OR_DIE(vals_out_bz,(2, 9, gams%natom3**2, gams%nqbz, gams%nsppol), ierr)
 vals_in_bz = zero
 vals_out_bz = zero

 ! Build tables needed by complete_gamma.
 call ephtk_mkqtabs(cryst, gams%nqibz, gams%qibz, gams%nqbz, gams%qbz, qirredtofull, qpttoqpt)

 ! Fill BZ array with IBZ data.
 do spin=1,gams%nsppol
   do iq_ibz=1,gams%nqibz
     iq_bz = qirredtofull(iq_ibz)
     vals_in_bz(:,:,:,iq_bz,spin) = reshape(gams%vals_in_qibz(:,:,:,:,iq_ibz,spin), [2, 9, gams%natom3**2])
     vals_out_bz(:,:,:,iq_bz,spin) = reshape(gams%vals_out_qibz(:,:,:,:,iq_ibz,spin), [2, 9, gams%natom3**2])
   end do
 end do

 ! Complete vals_bz in the full BZ.
 !TODO!!! rotate the vv in and out matrices, according to the symmetry operation, instead of just copying them
 !call complete_gamma_vv(cryst, gams%natom3, gams%nsppol, gams%nqibz, gams%nqbz,&
 !  gams%eph_scalprod, qirredtofull, qpttoqpt, vals_in_bz(:,:,:,:,spin))
 !call complete_gamma_vv(cryst, gams%natom3, gams%nsppol, gams%nqibz, gams%nqbz,&
 !  gams%eph_scalprod, qirredtofull, qpttoqpt, vals_out_bz(:,:,:,:,spin))

 ! TODO: replace the above with these calls from anaddb
 !call complete_gamma_tr(cryst,elph_ds%ep_scalprod,elph_ds%nbranch,elph_ds%nqptirred,&
 !&   elph_ds%nqpt_full,elph_ds%nsppol,elph_tr_ds%gamma_qpt_trout,elph_ds%qirredtofull,qpttoqpt)

 ! TODO: idem for vv_vals 3x3 matrices

 ABI_FREE(qirredtofull)
 ABI_FREE(qpttoqpt)

 !! This call is not executed in elphon!
 !if (gams%symgamma == 1) then
 !  do spin=1,gams%nsppol
 !    do iq_bz=1,gams%nqbz
 !      call tgamma_symm_vv(cryst, gams%qbz(:,iq_bz), vals_in_bz(:,:,:,iq_bz,spin))
 !      call tgamma_symm_vv(cryst, gams%qbz(:,iq_bz), vals_out_bz(:,:,:,iq_bz,spin))
 !    end do
 !  end do
 !end if

 ! Now FT to real space
 ! NOTE: gprim (not gprimd) is used for all FT interpolations,
 ! to be consistent with the dimensions of the rpt, which come from anaddb.
 ! TODO: this is needed only if FT is used, not when the linear interpolation is employed.
 !if (.not. allocated(gams%vals_in_rpt)) then
 ABI_MALLOC_OR_DIE(gams%vals_in_rpt, (2, 9, gams%natom3**2, gams%nrpt, gams%nsppol), ierr)
 ABI_MALLOC_OR_DIE(gams%vals_out_rpt, (2, 9, gams%natom3**2, gams%nrpt, gams%nsppol), ierr)
 gams%vals_in_rpt = zero
 gams%vals_out_rpt = zero

 ! q --> r
 ABI_MALLOC(coskr, (gams%nqbz,gams%nrpt))
 ABI_MALLOC(sinkr, (gams%nqbz,gams%nrpt))
 call ftgam_init(gams%gprim, gams%nqbz, gams%nrpt, gams%qbz, gams%rpt, coskr, sinkr)

 do spin=1,gams%nsppol
   do idir=1,3
     do jdir=1,3
       ii = idir+gams%ndir_transp*(jdir-1)
       ! TODO: this is no contiguous in memory and will be slow. Make adapted ftgam?
       call ftgam(gams%wghatm, vals_in_bz(:,ii,:,:,spin), gams%vals_in_rpt(:,ii,:,:,spin), gams%natom, gams%nqbz,&
         gams%nrpt, qtor1, coskr, sinkr)
       call ftgam(gams%wghatm, vals_out_bz(:,ii,:,:,spin), gams%vals_out_rpt(:,ii,:,:,spin), gams%natom, gams%nqbz,&
         gams%nrpt, qtor1, coskr, sinkr)
     end do
   end do
 end do

 ABI_FREE(vals_in_bz)
 ABI_FREE(vals_out_bz)
 ABI_FREE(coskr)
 ABI_FREE(sinkr)

end subroutine phgamma_vv_interp_setup
!!***
!----------------------------------------------------------------------

!!****f* m_phgamma/phgamma_linwid
!! NAME
!! phgamma_linwid
!!
!! FUNCTION
!!  Interpolate the phonon linewidths along an arbitrary q-path (use Fourier interpolation).
!!
!! INPUTS
!!  cryst<crystal_t>=Info on the unit cell and symmetries.
!!  ifc<ifc_type>=Interatomic force constants.
!!  ndivsm=Number of points used to sample the smallest segment
!!  nvert = Number of extrema in qverts
!!  qverts(3,nvert) = vertices of reciprocal space trajectory
!!  basename=name used to create the different output files (text format).
!!  ncid=Netcdf file handler (already open in the caller).
!!  comm=MPI communicator
!!
!! OUTPUT
!!  wminmax=Minimum and max phonon frequency obtained on the path (Hartree units)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine phgamma_linwid(gams, cryst, ifc, ndivsm, nvert, qverts, basename, ncid, wminmax, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nvert,ndivsm,comm,ncid
 type(crystal_t),intent(in) :: cryst
 type(ifc_type),intent(in) :: ifc
 class(phgamma_t),intent(inout) :: gams
 character(len=*),intent(in) :: basename
!arrays
 real(dp),intent(in) :: qverts(3,nvert)
 real(dp),intent(out) :: wminmax(2)

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: natom,ii,mu,iqpt,natom3,nsppol,ierr
 integer :: spin,unt,nqpt,nrpt,cnt,nproc,my_rank
#ifdef HAVE_NETCDF
 integer :: ncerr
#endif
 real(dp) :: omega_min,omega_max,wtmp,omega
 character(len=500) :: msg
 type(kpath_t) :: qpath
!arrays
 real(dp) :: gamma_spin(gams%nsppol),lambda_spin(gams%nsppol)
 real(dp) :: displ_cart(2,3*cryst%natom,3*cryst%natom)
 real(dp) :: phfrq(3*cryst%natom),gamma_ph(3*cryst%natom),lambda_ph(3*cryst%natom)
 real(dp) :: qpt(3),shift(3)
 real(dp),allocatable :: all_phfreq(:,:),all_gammaq(:,:,:),all_lambdaq(:,:,:),all_displ_cart(:,:,:,:)

! *********************************************************************

 DBG_ENTER("COLL")

 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 natom = cryst%natom; natom3 = gams%natom3; nsppol = gams%nsppol; nrpt = gams%nrpt

 ! Define the q-path along which phonon linwid will be interpolated.
 qpath = kpath_new(qverts, cryst%gprimd, ndivsm)
 nqpt = qpath%npts

 ! Allocate workspace arrays for MPI.
 ABI_CALLOC(all_phfreq, (natom3, nqpt))
 ABI_CALLOC(all_gammaq, (natom3, nqpt, nsppol))
 ABI_CALLOC(all_lambdaq, (natom3, nqpt, nsppol))
 ABI_CALLOC(all_displ_cart, (2, natom3, natom3, nqpt))

 ! initialize the minimum and maximum phonon frequency
 omega_min = huge(one); omega_max = -huge(one)

 ! Interpolation along specified path in q space (keep spin dep.)
 cnt = 0
 do spin=1,nsppol
   do iqpt=1,nqpt
     cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle
     call wrap2_pmhalf(qpath%points(:,iqpt), qpt, shift)

     ! Get phgamma
     call gams%interp(cryst, ifc, spin, qpt, phfrq, gamma_ph, lambda_ph, displ_cart)
     all_gammaq(:, iqpt, spin) = gamma_ph
     all_lambdaq(:, iqpt, spin) = lambda_ph
     if (spin == 1) then
       all_phfreq(:, iqpt) = phfrq
       all_displ_cart(:, :, :, iqpt) = displ_cart
     end if

     ! Find max/min phonon frequency along path chosen
     ! presumed to be representative of full BZ to within 10 percent
     omega_min = min(omega_min, phfrq(1))
     omega_max = max(omega_max, phfrq(natom3))
   end do ! end iqpt do
 end do ! spin

 ! Collect results
 if (omega_min > tol12) omega_min = zero
 wtmp = omega_min; call xmpi_min(wtmp, omega_min, comm, ierr)
 wtmp = omega_max; call xmpi_max(wtmp, omega_max, comm, ierr)
 wminmax = [omega_min, omega_max]

 call xmpi_sum_master(all_gammaq, master, comm, ierr)
 call xmpi_sum_master(all_lambdaq, master, comm, ierr)
 call xmpi_sum_master(all_phfreq, master, comm, ierr)
 call xmpi_sum_master(all_displ_cart, master, comm, ierr)

 ! Master writes text file with final results.
 ! 3 * natom blocks, one block for each phonon mode.
 ! Each block contains:
 !   iqpt  omega_(q) gamma(q) lambda(q) nesting(q) ....

 if (xmpi_comm_rank(comm) == master) then
   if (open_file(strcat(basename, '_PHGAMMA'), msg, newunit=unt, action="write", form="formatted", status="unknown") /= 0) then
     MSG_ERROR(msg)
   end if

   write(unt,'(a)')     '#'
   write(unt,'(a)')     '# ABINIT package: E-PH band structure file. Hartree units'
   write(unt,'(a)')     '#'
   write(unt,'(a,i0,a)')'# Phonon frequencies, ph linewidths and lambda calculated on ',nqpt,' q-points'
   call qpath%print(header="Description of the q-path:", unit=unt, pre="#")
   do ii=1,2; write(unt,'(a)')     "# "; end do

   write(unt,'(a,e16.6)')"# Total DOS at Fermi level ",sum(gams%n0)
   do spin=1,nsppol
     write(unt,"(a,i0,a,e16.6)")"# The DOS at Fermi level for spin ",spin," is ",gams%n0(spin)
   end do

   do mu=1,natom3
     write(unt,'(a)')"#"
     if (nsppol == 1) write(unt,'(a,i0,a)')"# phonon mode [",mu,"] q-index omega gamma lambda"
     if (nsppol == 2) write(unt,'(a,i0,a)')&
       "# phonon mode [",mu,"] q-index omega gamma_tot lambda_tot gamma[spin=1] lambda[spin=1] gamma[2] lambda[2]"
     write(unt,'(a)')"#"
     do iqpt=1,nqpt
       omega = all_phfreq(mu, iqpt)
       gamma_spin = all_gammaq(mu, iqpt, :)
       lambda_spin = all_lambdaq(mu, iqpt, :)
       if (nsppol == 1) then
         write(unt,'(i8,3es16.6)' )iqpt,omega,gamma_spin(1),lambda_spin(1)
       else
         write(unt,'(i8,es20.6,6es16.6)' )iqpt,omega,&
            sum(gamma_spin),sum(lambda_spin),&
            gamma_spin(1),lambda_spin(1),&
            gamma_spin(2),lambda_spin(2)
       end if
     end do
   end do

   close(unt)

   ! Write data to netcdf file
   if (ncid /= nctk_noid) then
#ifdef HAVE_NETCDF
     ncerr = nctk_def_dims(ncid, [&
       nctkdim_t("natom3", 3*natom), nctkdim_t("nqpath", nqpt), nctkdim_t("number_of_spins", nsppol) &
     ], defmode=.True.)
     NCF_CHECK(ncerr)

     ncerr = nctk_def_arrays(ncid, [&
       nctkarr_t('qpath', "dp", "number_of_reduced_dimensions, nqpath"), &
       nctkarr_t('phfreq_qpath', "dp", "natom3, nqpath"), &
       nctkarr_t('phdispl_cart_qpath', "dp", "two, natom3, natom3, nqpath"), &
       nctkarr_t('phgamma_qpath', "dp", "natom3, nqpath, number_of_spins"), &
       nctkarr_t('phlambda_qpath', "dp", "natom3, nqpath, number_of_spins")])
     NCF_CHECK(ncerr)

     NCF_CHECK(nctk_set_datamode(ncid))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "qpath"), qpath%points))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "phfreq_qpath"), all_phfreq))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "phdispl_cart_qpath"), all_displ_cart))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "phgamma_qpath"), all_gammaq))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "phlambda_qpath"), all_lambdaq))
#endif
   end if
 end if ! master

 ABI_FREE(all_phfreq)
 ABI_FREE(all_gammaq)
 ABI_FREE(all_lambdaq)
 ABI_FREE(all_displ_cart)

 call qpath%free()

 DBG_EXIT("COLL")

end subroutine phgamma_linwid
!!***

!----------------------------------------------------------------------

!!****f* m_phgamma/a2fw_free
!! NAME
!! a2fw_free
!!
!! FUNCTION
!!  Free the memory allocated in a2f
!!
!! SIDE EFFECTS
!!  a2f<a2fw_t>=Structure storing the Eliashberg function a2F.
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine a2fw_free(a2f)

!Arguments ------------------------------------
 class(a2fw_t),intent(inout) :: a2f

! *********************************************************************

 ! integer
 ABI_SFREE(a2f%qshift)

 ! real
 ABI_SFREE(a2f%n0)
 ABI_SFREE(a2f%omega)
 ABI_SFREE(a2f%vals)
 ABI_SFREE(a2f%vals_ee)
 ABI_SFREE(a2f%lambdaw)

end subroutine a2fw_free
!!***

!----------------------------------------------------------------------

!!****f* m_phgamma/a2fw_init
!! NAME
!! a2fw_init
!!
!! FUNCTION
!!  Calculates the FS averaged alpha^2F(w) function
!!
!! INPUTS
!!  cryst<crystal_t>=Info on the unit cell.
!!  ifc<ifc_type>=Interatomic force constants.
!!  gams<phgamma_t>=Structure storing the phonon linewidths.
!!  wstep=Step for linear frequency mesh in Ha.
!!  wminmax(2)=Minimum and maximum phonon frequency. Used to construct the linear mesh for A2F(w).
!!  intmeth=Integration method: 1 for gaussian, 2 for tetrahedra
!!  smear=Gaussian broadening used to approximate the Dirac delta.
!!  ngqpt(3)=Divisions of the Q-mesh used for interpolating the phonon linewidths (see also nqshift and qshift).
!!  nqshift=Number of shifts used to generated the Q-mesh.
!!  qshift(3,nqshift)=The shifts.
!!  comm=MPI communicator
!!  [qintp]=If set to False, ngqgpt, nqshift, qshift and qptop  are ignored and
!!     A2F(w) is computed from the IBZ values stored in gams. Default: True i.e use Fourier interpolation.
!!  [qptopt]=Controls the generation of the q-points. If not specified, the routine takes fully into account
!!    the symmetries of the system to generate the q points in the IBZone i.e. qptopt=1
!!    Other values of qptopt can be used for debugging purpose.
!!
!! OUTPUT
!!  a2f<a2fw_t>=Structure storing the Eliashberg function a2F(w).
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine a2fw_init(a2f, gams, cryst, ifc, intmeth, wstep, wminmax, smear, ngqpt, nqshift, qshift, comm, &
  qintp, qptopt) ! optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: intmeth,nqshift,comm
 integer,intent(in),optional :: qptopt
 real(dp),intent(in) :: wstep,smear
 logical,optional,intent(in) :: qintp
 type(phgamma_t),intent(inout) :: gams
 type(ifc_type),intent(in) :: ifc
 type(crystal_t),intent(in) :: cryst
 type(a2fw_t),target,intent(out) :: a2f
!arrays
 integer,intent(in) :: ngqpt(3)
 real(dp),intent(in) :: wminmax(2),qshift(3,nqshift)

!Local variables -------------------------
!scalars
 integer,parameter :: master = 0
 integer :: my_qptopt,iq_ibz,nqibz,ount,my_rank,nproc,cnt
 integer :: mu,iw,natom3,nsppol,spin,ierr,nomega,nqbz
 integer :: iene, jene, itemp, ntemp, jene_jump
 real(dp) :: cpu,wall,gflops
 real(dp) :: lambda_iso,omega,omega_log,xx,omega_min,omega_max,ww,mustar,tc_macmill
 real(dp) :: temp_el, min_temp, delta_temp, chempot, ene1, ene2, G0
 logical :: do_qintp
 character(len=500) :: msg
 type(htetra_t) :: qtetra
!arrays
 integer :: qptrlatt(3,3),new_qptrlatt(3,3)
 real(dp) :: displ_cart(2,3,cryst%natom,3*cryst%natom)
 real(dp) :: phfrq(gams%natom3),gamma_ph(gams%natom3),lambda_ph(gams%natom3)
 real(dp) :: invphfrq(gams%natom3)
 real(dp),allocatable :: my_qshift(:,:), gamma_ph_ee(:,:,:,:), tmp_a2f(:)
 real(dp), ABI_CONTIGUOUS pointer :: a2f_1d(:)
 real(dp),allocatable :: qibz(:,:),wtq(:),qbz(:,:)
 real(dp),allocatable :: a2f_1mom(:),a2flogmom(:),a2flogmom_int(:),wdt(:,:)
 real(dp),allocatable :: lambda_tetra(:,:,:),phfreq_tetra(:,:), tmp_gaussian(:,:)
 real(dp), allocatable :: a2feew_partial(:), a2feew_partial_int(:), a2feew_w(:), a2feew_w_int(:)

! *********************************************************************

 call cwtime(cpu, wall, gflops, "start")
 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

 my_qptopt = 1; if (present(qptopt)) my_qptopt = qptopt
 do_qintp = .True.; if (present(qintp)) do_qintp = qintp
 nsppol = gams%nsppol; natom3 = gams%natom3

 if (gams%prteliash == 3) then
   ABI_MALLOC(gamma_ph_ee, (gams%nene, gams%nene, gams%natom3, gams%nsppol))
 end if

 if (do_qintp) then
   ! Generate fine q-mesh, findthe IBZ and the corresponding weights.
   qptrlatt = 0; qptrlatt(1,1) = ngqpt(1); qptrlatt(2,2) = ngqpt(2); qptrlatt(3,3) = ngqpt(3)

   call kpts_ibz_from_kptrlatt(cryst, qptrlatt, my_qptopt, nqshift, qshift, nqibz, qibz, wtq, nqbz, qbz, &
      new_kptrlatt=new_qptrlatt, new_shiftk=my_qshift)
   ABI_FREE(qbz)

   ! Store quantities that cannot be easily (and safely) calculated if we only know the IBZ.
   a2f%ngqpt = ngqpt; a2f%nqshift = size(my_qshift, dim=2)
   ABI_MALLOC(a2f%qshift, (3, a2f%nqshift))
   a2f%qshift = my_qshift
   ABI_FREE(my_qshift)

 else
   ! No interpolation. Use q-mesh parameters from gams%
   a2f%ngqpt = gams%ngqpt; a2f%nqshift = 1
   nqibz = gams%nqibz
   ABI_MALLOC(qibz, (3, nqibz))
   ABI_MALLOC(wtq, (nqibz))
   ABI_MALLOC(a2f%qshift, (3, a2f%nqshift))
   qibz = gams%qibz; wtq = gams%wtq
   a2f%qshift = zero ! Note: assuming q-mesh centered on Gamma.
 end if

 call cwtime_report(" a2fw_init, q-setup:", cpu, wall, gflops)

 ! Define Min and max frequency for the mesh (enlarge it a bit)
 omega_min = wminmax(1); omega_max = wminmax(2)
 omega_min = omega_min - 0.1*abs(omega_min)
 if (omega_min >= zero) omega_min = one/Ha_meV
 omega_max = omega_max + 0.1*abs(omega_max)

 a2f%nsppol = nsppol; a2f%natom3 = gams%natom3; a2f%smear = smear
 a2f%omega_min = omega_min; a2f%omega_max = omega_max
 nomega = int((omega_max - omega_min) / wstep); a2f%nomega = nomega; a2f%wstep = wstep
 a2f%nene = gams%nene
 a2f%enemin = gams%enemin
 a2f%deltaene = gams%deltaene

 ABI_MALLOC(a2f%n0, (nsppol))
 a2f%n0 = gams%n0
 ! Build linear mesh.
 ABI_MALLOC(a2f%omega, (nomega))
 a2f%omega = arth(omega_min, wstep, nomega)
 ABI_CALLOC(a2f%vals, (nomega,0:natom3, nsppol))
 ABI_CALLOC(a2f%lambdaw, (nomega,0:natom3, nsppol))

 if (gams%prteliash == 3) then
   ABI_CALLOC(a2f%vals_ee, (gams%nene, gams%nene, nomega, nsppol))
 end if

 ABI_MALLOC(tmp_a2f, (nomega))

 if (intmeth == 2) then
   call cwtime(cpu, wall, gflops, "start")

   ! Prepare tetrahedron integration.
   qptrlatt = 0; qptrlatt(1, 1) = a2f%ngqpt(1); qptrlatt(2, 2) = a2f%ngqpt(2); qptrlatt(3, 3) = a2f%ngqpt(3)
   qtetra = tetra_from_kptrlatt(cryst, my_qptopt, qptrlatt, a2f%nqshift, a2f%qshift, nqibz, qibz, comm, msg, ierr)
   if (ierr /= 0) MSG_ERROR(msg)

   ABI_CALLOC(lambda_tetra, (nqibz, natom3, nsppol))
   ABI_CALLOC(phfreq_tetra, (nqibz, natom3))
   cnt = 0
   do iq_ibz = 1, nqibz
     cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle
     call ifc%fourq(cryst, qibz(:,iq_ibz), phfrq, displ_cart)
     ! save for tetrahedron interpolation
     phfreq_tetra(iq_ibz,:) = phfrq(:)
   end do
   call xmpi_sum(phfreq_tetra, comm, ierr)

   call cwtime_report(" a2fw_init%tetra", cpu, wall, gflops)
 end if

 call cwtime(cpu, wall, gflops, "start")

 ! DEV_MJV
 !open (unit=900, file="a2fvals_ee.dat")
 !write (900,*) '# do_qintp ', do_qintp

 ! Loop over spins and qpoints in the IBZ. For the moment parallelize over iq_ibz
 do spin=1,nsppol
   cnt = 0
   do iq_ibz=1,nqibz
     ! TODO: for the moment the memory is not distributed, only the calculation exception is vals_ee
     cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle

     ! Interpolate or evaluate gamma directly.
     if (do_qintp) then
       if (gams%prteliash == 3) then
         call gams%interp(cryst, ifc, spin, qibz(:,iq_ibz), phfrq, gamma_ph, lambda_ph, displ_cart, gamma_ph_ee=gamma_ph_ee)
       else
         call gams%interp(cryst, ifc, spin, qibz(:,iq_ibz), phfrq, gamma_ph, lambda_ph, displ_cart)
       end if
     else
       if (gams%prteliash == 3) then
         call gams%eval_qibz(cryst, ifc, iq_ibz, spin, phfrq, gamma_ph, lambda_ph, displ_cart, gamma_ph_ee=gamma_ph_ee)
       else
         call gams%eval_qibz(cryst, ifc, iq_ibz, spin, phfrq, gamma_ph, lambda_ph, displ_cart)
       end if
     end if

     select case (intmeth)
     case (1)
       ! Gaussian: Add all contributions from the phonon modes at this qpoint to a2f
       ! (note that unstable modes are included).
       ABI_MALLOC(tmp_gaussian, (nomega, natom3))
       do mu=1,natom3
         tmp_a2f = zero
         do iw=1,nomega
           xx = a2f%omega(iw) - phfrq(mu)
           tmp_gaussian(iw, mu) = gaussian(xx, smear)
           tmp_a2f(iw) = tmp_a2f(iw) + tmp_gaussian(iw,mu) * lambda_ph(mu) * abs(phfrq(mu))
         end do
         a2f%vals(:,mu,spin) = a2f%vals(:,mu,spin) + tmp_a2f * wtq(iq_ibz)
       end do

       if (gams%prteliash == 3) then
         ! reset phfrq for low freq modes
         invphfrq = zero
         do mu=1,natom3
           if (abs(phfrq(mu)) > EPHTK_WTOL) invphfrq(mu) = one / abs(phfrq(mu))
         end do

         do iene= 1, gams%nene
           do jene= 1, gams%nene
             tmp_a2f = zero
             ! TODO: following block is just a GEMM
             do mu=1,natom3
               do iw=1,nomega
                 tmp_a2f(iw) = tmp_a2f(iw) + tmp_gaussian(iw,mu) * gamma_ph_ee(jene,iene,mu,spin) * invphfrq(mu)
               end do
             end do

             a2f%vals_ee(jene, iene, :, spin) = a2f%vals_ee(jene, iene, :, spin) + tmp_a2f(:) * wtq(iq_ibz)
             !if (iene == gams%nene/2 .and. jene == gams%nene/2) then
             !  write (900, '(a,E20.10,2x,2x,I6,3E20.10)') '#', wtq(iq_ibz), iq_ibz, invphfrq(1:3)
             !  do iw=1,nomega
             !    write (900, '(i6,2x,E20.10,2x,3E20.10,2x,3E20.10)') iw, a2f%vals_ee(jene, iene, iw, spin), &
             !        gamma_ph_ee(jene,iene,:,spin), tmp_gaussian(iw,1:3)
             !  end do
             !  write (900,*)
             !end if
           end do
         end do
       end if

       ABI_FREE(tmp_gaussian)

     case (2)
       ! Tetra: store data.
       do mu=1,natom3
         lambda_tetra(iq_ibz, mu, spin) = lambda_ph(mu) * abs(phfrq(mu))
       end do

       if (gams%prteliash == 3) then
         MSG_ERROR("Eliashberg function with tetra not coded")
       end if

     case default
       MSG_ERROR(sjoin("Wrong intmeth:", itoa(intmeth)))
     end select

   end do ! iq_ibz
 end do ! spin

 !DEV_MJV
 !close(900)

 if (intmeth == 2) then
   ! workspace for tetra.
   ABI_MALLOC(wdt, (nomega, 2))

   ! For each mode get its contribution
   do spin=1,nsppol
     do mu=1,natom3
       cnt = 0
       do iq_ibz=1,nqibz
         ! NB: if we are interpolating the gamma, nqibz > gams%nqibz
         cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle

         call qtetra%get_onewk(iq_ibz, gams%bcorr, nomega, nqibz, phfreq_tetra(:,mu), omega_min, omega_max, one, wdt)
         wdt = wdt * wtq(iq_ibz)

         ! Accumulate (Integral of a2F is computed afterwards)
         a2f%vals(:,mu,spin) = a2f%vals(:,mu,spin) + wdt(:,1) * lambda_tetra(iq_ibz,mu,spin)
         !a2f%lambdaw(:,mu,spin) = a2f%lambdaw(:,mu,spin) + wdt(:,2) * lambda_tetra(iq_ibz, mu, spin)
      end do
     end do
   end do

   ! Free memory allocated for tetra.
   ABI_FREE(wdt)
   ABI_FREE(lambda_tetra)
   ABI_FREE(phfreq_tetra)
   call qtetra%free()
 end if

 ! Collect final results on each node
 call xmpi_sum(a2f%vals, comm, ierr)
 do spin=1,nsppol
   a2f%vals(:,0,spin) = sum(a2f%vals(:,1:natom3,spin), dim=2)
   ! previously would divide by g(eF, spin)
   !a2f%vals(:,:,spin) = a2f%vals(:,:,spin) / (two_pi*a2f%n0(spin))
 end do

 if (gams%prteliash == 3) then
   ! For the moment vals_ee only works with gaussians
   call xmpi_sum(a2f%vals_ee, comm, ierr)
   do spin=1,nsppol
     a2f%vals_ee(:,:,:,spin) = a2f%vals_ee(:,:,:,spin) / (two * pi * gams%n0(spin))
   end do
 end if

 call cwtime_report(" a2fw_init, a2f_eval:", cpu, wall, gflops)

 ABI_MALLOC(a2f_1mom, (nomega))
 ABI_MALLOC(a2flogmom, (nomega))
 ABI_MALLOC(a2flogmom_int, (nomega))

 !call a2fw_print_info()

 ! Compute lambda(w) (resolved in spin and phonon mode).
 do spin=1,nsppol
   do mu=0,natom3

     do iw=1,nomega
       ww = a2f%omega(iw)
       if (abs(ww) > EPHTK_WTOL) then
          a2f_1mom(iw) = a2f%vals(iw,mu,spin) / abs(ww)
       else
          a2f_1mom(iw) = zero
       end if
       ! TODO: Strange that the first value of int(f) is not set to zero!
       ! FIXME: The frequency integration overestimates lambda if there are negative phonon frequencies!
     end do

     call simpson_int(nomega, wstep, a2f_1mom, a2f%lambdaw(:,mu,spin))
   end do
 end do

 ! print log moments of the alpha 2 F functions
 do spin=1,nsppol
   a2f_1d => a2f%vals(:,0,spin)

   lambda_iso = a2f%get_moment(0, spin)

   ! Get log moment of alpha^2F.
   a2flogmom = zero
   do iw=1,nomega
     omega = a2f%omega(iw)
     if (abs(omega) > EPHTK_WTOL) then
       a2flogmom(iw) = (two / lambda_iso) * a2f_1d(iw) * log(abs(omega)) / abs(omega)
       ! I think this is the correct expression.
       !a2flogmom(iw) = (one / lambda_iso) * a2f_1d(iw) * log(abs(omega)) / abs(omega)
     end if
   end do
   call simpson_int(nomega, wstep, a2flogmom, a2flogmom_int)
   omega_log = exp(a2flogmom_int(nomega))

   mustar = 0.12
   tc_macmill = omega_log/1.2_dp * exp((-1.04_dp*(one+lambda_iso)) / (lambda_iso-mustar*(one+0.62_dp*lambda_iso)))

   if (my_rank == master) then
     ount = std_out
     if (nsppol > 1) then
       write(msg,'(3a)') ch10,'Warning: some of the following quantities should be integrated over spin', ch10
       call wrtout(ount, msg)
     end if

     write(ount,'(a)')' Superconductivity: isotropic evaluation of parameters from electron-phonon coupling.'
     write(ount,'(a,es16.6)')' isotropic lambda = ',lambda_iso
     write(ount,'(a,es16.6,a,es16.6,a)' )' omegalog  = ',omega_log,' (Ha) ', omega_log * Ha_K, ' (Kelvin) '
     write(ount,'(a,es16.6,a,es16.6,a)')' MacMillan Tc = ',tc_macmill,' (Ha) ', tc_macmill * Ha_K, ' (Kelvin) '
     write(ount,"(a)")'    positive moments of alpha2F:'
     write(ount,'(a,es16.6)' )' lambda <omega^2> = ',a2f%get_moment(2, spin)
     write(ount,'(a,es16.6)' )' lambda <omega^3> = ',a2f%get_moment(3, spin)
     write(ount,'(a,es16.6)' )' lambda <omega^4> = ',a2f%get_moment(4, spin)
     write(ount,'(a,es16.6)' )' lambda <omega^5> = ',a2f%get_moment(5, spin)
   end if
 end do

 ! Calculate the temperature dependence of the a2f(e,e',w) integrals (G_0(T_e) in PRL 110 016405 (2013) [[cite:Arnaud2013]])
 if (gams%prteliash == 3 .and. my_rank == master) then
   ntemp = 100
   min_temp = zero
   delta_temp = 40._dp ! Kelvin
   ABI_MALLOC(a2feew_partial, (a2f%nene))
   ABI_MALLOC(a2feew_partial_int, (a2f%nene))
   ABI_MALLOC(a2feew_w, (nomega))
   ABI_MALLOC(a2feew_w_int, (nomega))

   if (open_file("EPC_strength_aafo_T.dat", msg, newunit=ount, form="formatted", action="write", status="unknown") /= 0) then
     MSG_ERROR(msg)
   end if

   write(ount, "(a)")"# temp_el, G_0(T_e) in W/m^3/K, spin"
   do spin=1,nsppol
     do itemp = 1, ntemp
       temp_el = min_temp + (itemp-1)*delta_temp

       ! TODO: need to evolve the chemical potential with T, but I do not have access to the full DOS here!!
       ! possible fix using local information on DOS variation near E_F...
       chempot = a2f%enemin + half * a2f%nene * a2f%deltaene

       do iw=1,nomega
         omega = a2f%omega(iw)
         jene_jump = nint(omega/a2f%deltaene)
         a2feew_partial = zero
         do iene=1,a2f%nene
           ene1 = a2f%enemin + (iene-1)*a2f%deltaene
           ene2 = ene1 + omega
           a2feew_partial(iene) = a2f%vals_ee(min(a2f%nene,iene+jene_jump), iene, iw, spin) * &
              (fermi_dirac(ene1, chempot, temp_el/Ha_K) - fermi_dirac(ene2, chempot, temp_el/Ha_K))
         end do
         call simpson_int(a2f%nene, a2f%deltaene, a2feew_partial, a2feew_partial_int)
         a2feew_w(iw) = a2feew_partial_int(a2f%nene)
       end do
       call simpson_int(nomega, wstep, a2feew_w, a2feew_w_int)
       G0 = a2feew_w_int(nomega) * two_pi * a2f%n0(spin) / cryst%ucvol
       ! conversion factor for G0 to SI units =  Ha_J / Time_Sec / (Bohr_meter)**3 ~ 1.2163049915755545e+30
       write(ount, "(2(e20.10,2x))") temp_el, G0  * kb_HaK / Time_Sec / (Bohr_meter)**3, spin !* Ha_J???
     end do
   end do
   close(ount)

   ABI_FREE(a2feew_partial)
   ABI_FREE(a2feew_partial_int)
   ABI_FREE(a2feew_w)
   ABI_FREE(a2feew_w_int)
   ABI_FREE(gamma_ph_ee)
 end if

 ABI_FREE(tmp_a2f)
 ABI_FREE(a2f_1mom)
 ABI_FREE(a2flogmom)
 ABI_FREE(a2flogmom_int)
 ABI_FREE(qibz)
 ABI_FREE(wtq)

end subroutine a2fw_init
!!***

!----------------------------------------------------------------------

!!****f* m_phgamma/a2fw_get_moment
!! NAME
!! a2fw_get_moment
!!
!! FUNCTION
!!  Compute \int dw [a2F(w)/w] w^n
!!  From Allen PRL 59 1460 [[cite:Allen1987]] (See also [[cite:Grimvall1981]], Eq 6.72 page 175)
!!
!! INPUTS
!!  a2f<a2fw_t>=Structure storing the Eliashberg function.
!!  nn=Value of n
!!  spin=The spin component
!!
!! OUTPUT
!!  a2fw_get_moment = \int dw [a2F(w)/w] w^n
!!  [out_int(x)] = \int^{x} dw [a2F(w)/w] w^n
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

real(dp) function a2fw_get_moment(a2f, nn, spin, out_int)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spin,nn
 class(a2fw_t),intent(in) :: a2f
!arrays
 real(dp),intent(out),optional :: out_int(a2f%nomega)

!Local variables -------------------------
!scalars
 integer :: iw
 real(dp) :: omg, omg_nm1
!arrays
 real(dp) :: ff(a2f%nomega),int_ff(a2f%nomega)

! *********************************************************************

 ! Construct the integrand function. [a2F(w)/w] w^n
 ff = zero; int_ff = zero

 if (nn - 1 >= 0) then
   do iw=1,a2f%nomega
     omg = a2f%omega(iw)
     omg_nm1 = omg ** (nn-1)
     ff(iw) = a2f%vals(iw,0,spin) * omg_nm1
   end do
 else
   do iw=1,a2f%nomega
     omg = a2f%omega(iw)
     omg_nm1 = zero; if (abs(omg) > EPHTK_WTOL) omg_nm1 = omg**(nn-1)
     ff(iw) = a2f%vals(iw,0,spin) * omg_nm1
   end do
 end if

 ! Integration with simpson rule on a linear mesh.
 call simpson_int(a2f%nomega, a2f%wstep, ff, int_ff)

 a2fw_get_moment = int_ff(a2f%nomega)
 if (present(out_int)) out_int = int_ff

end function a2fw_get_moment
!!***

!----------------------------------------------------------------------

!!****f* m_phgamma/a2fw_tr_moment
!! NAME
!! a2fw_tr_moment
!!
!! FUNCTION
!!  Compute \int dw [a2F_tr(w)/w] w^n
!!  From Allen PRL 59 1460 [[cite:Allen1987]] and later PRB papers (See also [[cite:Grimvall1981]] book)
!!
!! INPUTS
!!  a2f_tr<a2fw_tr_t>=Structure storing the Eliashberg function.
!!  nn=Value of n
!!  spin=The spin component
!!
!! OUTPUT
!!  a2fw_tr_moment = \int dw [a2F_tr(w)/w] w^n
!!  [out_int(x)] = \int^{x} dw [a2F_tr(w)/w] w^n
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function a2fw_tr_moment(a2f_tr, nn, spin, out_int)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spin,nn
 type(a2fw_tr_t),intent(in) :: a2f_tr
!arrays
 real(dp),intent(out),optional :: out_int(a2f_tr%nomega,3,3)
 real(dp), dimension(3,3) :: a2fw_tr_moment

!Local variables -------------------------
!scalars
 integer :: iw, idir, jdir
 real(dp) :: omg,omg_nm1
!arrays
 real(dp) :: ff(a2f_tr%nomega),int_ff(a2f_tr%nomega)

! *********************************************************************

 do jdir = 1, 3
   do idir = 1, 3
     ! Construct the integrand function. [a2F_tr(w)/w] w^n
     ff = zero; int_ff = zero

     if (nn-1 >= 0) then
       do iw=1,a2f_tr%nomega
         omg = a2f_tr%omega(iw)
         omg_nm1 = omg ** (nn-1)
         ff(iw) = a2f_tr%vals_tr(iw,idir,jdir,0,spin) * omg_nm1
       end do
     else
       do iw=1,a2f_tr%nomega
         omg = a2f_tr%omega(iw)
         omg_nm1 = zero; if (abs(omg) > EPHTK_WTOL) omg_nm1 = omg**(nn-1)
         ff(iw) = a2f_tr%vals_tr(iw,idir,jdir,0,spin) * omg_nm1
       end do
     end if

     ! Integration with simpson rule on a linear mesh.
     call simpson_int(a2f_tr%nomega, a2f_tr%wstep, ff, int_ff)

     a2fw_tr_moment(idir, jdir) = int_ff(a2f_tr%nomega)
     if (present(out_int)) out_int(:,idir,jdir) = int_ff(:)
   end do
 end do

end function a2fw_tr_moment
!!***

!----------------------------------------------------------------------

!!****f* m_phgamma/a2fw_logmoment
!! NAME
!! a2fw_logmoment
!!
!! FUNCTION
!!
!! INPUTS
!!  a2f<a2fw_t>=Structure storing the Eliashberg function.
!!  nn
!!  spin=Spin index
!!
!! OUTPUT
!!  a2fw_logmoment =
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

real(dp) function a2fw_logmoment(a2f,spin) result(res)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spin
 type(a2fw_t),intent(in) :: a2f

!Local variables -------------------------
!scalars
 integer :: iw
 real(dp) :: omg,lambda_iso
!arrays
 real(dp) :: a2flogmom(a2f%nomega) !,a2flogmom_int(a2f%nomega)

! *********************************************************************

 ! Get log moment of alpha^2F.
 a2flogmom = zero
 lambda_iso = a2f%get_moment(0, spin)

 do iw=1,a2f%nomega
   omg = a2f%omega(iw)
   if (abs(omg) > EPHTK_WTOL) then
     !a2flogmom(iw) = (two/lambda_iso) * a2f%vals(iw,spin) * log(abs(omg))/abs(omg)
   end if
 end do

 !call simpson_int(nomega, a2f%wstep, a2flogmom, a2flogmom_int)
 !res = exp(a2flogmom_int(nomega))
 res = zero

end function a2fw_logmoment
!!***

!----------------------------------------------------------------------

!!****f* m_phgamma/a2fw_lambda_wij
!! NAME
!! a2fw_lambda_wij
!!
!! FUNCTION
!!  Compute \int dw [w x a2F(w)] (w_i-w_j)**2 + w**2)
!!
!! INPUTS
!!  a2f<a2fw_t>=Structure storing the Eliashberg function.
!!  w1,wj=Matsubara frequencies (real) in Hartree
!!  spin=Spin index
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

real(dp) function a2fw_lambda_wij(a2f,wi,wj,spin) result(res)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spin
 real(dp),intent(in) :: wi,wj
 type(a2fw_t),intent(in) :: a2f

!Local variables -------------------------
!scalars
 integer :: iw
 real(dp) :: wd2,vv,omg
 logical :: iszero
!arrays
 real(dp) :: values(a2f%nomega)

! *********************************************************************

 wd2 = (wi - wj)** 2
 iszero = (abs(wi - wj) < EPHTK_WTOL)

 do iw=1,a2f%nomega
   omg = a2f%omega(iw)
   vv = a2f%vals(iw,0,spin)
   if (abs(a2f%omega(iw)) > EPHTK_WTOL) then
     values(iw) = vv * omg / (wd2 + omg**2)
   else
     if (iszero) then ! TODO
       values(iw) = zero
     else
       values(iw) = vv * omg / (wd2 + omg**2)
     end if
   end if
 end do

 res = simpson(a2f%wstep, values)
 if (res < zero) res = zero

end function a2fw_lambda_wij
!!***

!----------------------------------------------------------------------

!!****f* m_phgamma/a2fw_write
!! NAME
!! a2fw_write
!!
!! FUNCTION
!!  Write alpha^2F(w) to an external file in text form
!!
!! INPUTS
!!  a2f<a2fw_t>=Container storing the Eliashberg functions.
!!  basename=Filename for output.
!!  post=String appended to netcdf variables e.g. _qcoarse, _qintp
!!  ncid=Netcdf file handler. Set it to nctk_noid to disable output.
!!
!! OUTPUT
!!  Output is written to file. This routine should be called by one MPI proc.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine a2fw_write(a2f, basename, post, ncid)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 character(len=*),intent(in) :: basename, post
 class(a2fw_t),intent(in) :: a2f

!Local variables -------------------------
!scalars
 integer :: iw,spin,unt,ii,mu
#ifdef HAVE_NETCDF
 integer :: ncerr
 character(len=500) :: dim1_name
#endif
 character(len=500) :: msg
 character(len=fnlen) :: path

! *********************************************************************

 ! Write spin-resolved a2F(w)
 path = strcat(basename, "_A2FW")
 if (open_file(path, msg, newunit=unt, form="formatted", action="write", status="unknown") /= 0) then
   MSG_ERROR(msg)
 end if

 call write_a2fw_header()

 if (a2f%nsppol == 1) then
   write(unt,'(a)')"# Frequency, a2F(w), lambda(w)"
   do iw=1,a2f%nomega
     write(unt,*) a2f%omega(iw), a2f%vals(iw,0,1), a2f%lambdaw(iw,0,1)
   end do

 else
   write(unt,'(a)')"# Frequency, a2F_tot(w), lambda_tot(w)dw, a2F_spin1(w), lambda_spin1(w) ..."
   do iw=1,a2f%nomega
     write(unt,*) a2f%omega(iw), &
                  sum(a2f%vals(iw,0,:)) , sum(a2f%lambdaw(iw,0,:)), &  ! TOT
                  a2f%vals(iw,0,1)      , a2f%lambdaw(iw,0,1),      &  ! UP
                  a2f%vals(iw,0,2)      , a2f%lambdaw(iw,0,2)          ! DOWN
   end do
 end if

 close(unt)

 ! Write phonon contributions to a2F(w)
 path = strcat(basename, "_PH_A2FW")
 if (open_file(path, msg, newunit=unt, form="formatted", action="write", status="unknown") /= 0) then
   MSG_ERROR(msg)
 end if
 call write_a2fw_header()

 if (a2f%nsppol == 1) then
   do mu=0,a2f%natom3
     write(unt,'(a,i0)')"# Phonon mode ",mu
     write(unt,'(a)')"# Frequency, a2F(w), lambda(w)"
     do iw=1,a2f%nomega
       write(unt,*) a2f%omega(iw), a2f%vals(iw,mu,1), a2f%lambdaw(iw,mu,1)
     end do
     do ii=1,2; write(unt,'(a)')""; end do
   end do

 else
   do mu=0,a2f%natom3
     write(unt,'(a,i0)')"# Phonon mode ",mu
     write(unt,'(a)')"# Frequency, a2F_tot(w), lambda_tot(w)dw, a2F_spin1(w), lambda_spin1(w) ..."
     do iw=1,a2f%nomega
       write(unt,*) a2f%omega(iw), &
                    sum(a2f%vals(iw,mu,:)), sum(a2f%lambdaw(iw,mu,:)), &  ! TOT
                    a2f%vals(iw,mu,1)     , a2f%lambdaw(iw,mu,1),      &  ! UP
                    a2f%vals(iw,mu,2)     , a2f%lambdaw(iw,mu,2)          ! DOWN
     end do
   end do
 end if

 close(unt)

 if (ncid /= nctk_noid) then
#ifdef HAVE_NETCDF
   ! Define dimensions.
   dim1_name = strcat("a2f_nomega", post)
   ncerr = nctk_def_dims(ncid, [ &
     nctkdim_t(dim1_name, a2f%nomega), &
     nctkdim_t("natom3p1", a2f%natom3 + 1), nctkdim_t("number_of_spins", a2f%nsppol) &
     ], defmode=.True.)
   NCF_CHECK(ncerr)

   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t(strcat('a2f_mesh', post), "dp", dim1_name), &
     nctkarr_t(strcat('a2f_values', post), "dp", strcat(dim1_name, ", natom3p1, number_of_spins")), &
     nctkarr_t(strcat('a2f_lambdaw', post), "dp", strcat(dim1_name, ", natom3p1, number_of_spins")) &
     ])
   NCF_CHECK(ncerr)

   NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, strcat("a2f_mesh", post)), a2f%omega))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, strcat("a2f_values", post)), a2f%vals))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, strcat("a2f_lambdaw", post)), a2f%lambdaw))
#endif
 end if

contains

subroutine write_a2fw_header()

 ! Output the header.
 write(unt,'(a)')              '#'
 write(unt,'(a)')              '# ABINIT package: a2F(w) file'
 write(unt,'(a)')              '#'
 write(unt,'(a)')              '# a2F(w) function integrated over the FS. omega in a.u.'
 write(unt,'(2a)')             '# ngqpt: ',trim(ltoa(a2f%ngqpt))
 write(unt,'(a,i0,3(a,e16.6))')'# number of frequencies: ',a2f%nomega," between omega_min: ",a2f%omega_min,&
                               ' Ha and omega_max: ',a2f%omega_max,' Ha with step:',a2f%wstep
 write(unt,'(a,e16.6)')         '#  the smearing width for gaussians is ',a2f%smear
 write(unt,'(a,e16.6)')"# Total DOS at Fermi level ",sum(a2f%n0)
 do spin=1,a2f%nsppol
   write(unt,"(a,i0,a,e16.6)")"# The DOS at Fermi level for spin ",spin," is ",a2f%n0(spin)
 end do
 do ii=1,2; write(unt,'(a)')     "# "; end do

end subroutine write_a2fw_header

end subroutine a2fw_write
!!***

!----------------------------------------------------------------------
!!****f* m_phgamma/a2fw_ee_write
!! NAME
!! a2fw_ee_write
!!
!! FUNCTION
!!  Write alpha^2F(e,e',w) to an external file in text form
!!
!! INPUTS
!!  a2f<a2fw_t>=Container storing the Eliashberg functions.
!!  basename=Filename for output.
!!
!! OUTPUT
!!  Output is written to file. This routine should be called by one MPI proc.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine a2fw_ee_write(a2f, basename)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: basename
 type(a2fw_t),intent(in) :: a2f

!Local variables -------------------------
!scalars
 integer :: iw, spin, unt, ii, iene, jene
 real(dp) :: ene1, ene2
 character(len=500) :: msg
 character(len=fnlen) :: path

! *********************************************************************

 ! Write spin-resolved a2F(e,e',w)
 path = strcat(basename, "_A2FEEW")
 if (open_file(path, msg, newunit=unt, form="formatted", action="write", status="unknown") /= 0) then
   MSG_ERROR(msg)
 end if

 call write_a2fw_header()

 write(unt,'(a)')"# en2, en1, Frequency, a2F_tot(w) (presently summed over spin)"
 do iw=1,a2f%nomega
   do iene=1,a2f%nene
     ene1 = a2f%enemin + (iene-1)*a2f%deltaene
     do jene=1,a2f%nene
       ene2 = a2f%enemin + (jene-1)*a2f%deltaene
       write(unt,'(3E20.10,2x,E20.10)') ene2, ene1, a2f%omega(iw), sum(a2f%vals_ee(jene,iene,iw,:))   ! TOT
     end do
   end do
 end do

 close(unt)

 ! Write spin-resolved a2F(ef, ef, w)
 path = strcat(basename, "_A2FW_reference")
 if (open_file(path, msg, newunit=unt, form="formatted", status="unknown") /= 0) then
   MSG_ERROR(msg)
 end if

 call write_a2fw_header()

 write(unt,'(a)')"# Frequency, a2F_tot(ef,ef,w) for comparison with normal a2F(w) (presently summed over spin)"
 iene = int(a2f%nene/2)
 jene = int(a2f%nene/2)
 do iw=1,a2f%nomega
   write(unt,'(E20.10,2x,E20.10)') a2f%omega(iw), sum(a2f%vals_ee(jene,iene,iw,:))   ! TOT
 end do

 close(unt)

contains

subroutine write_a2fw_header()

 ! Output the header.
 write(unt,'(a)')              '#'
 write(unt,'(a)')              '# ABINIT package: a2F(e,eprime,w) file'
 write(unt,'(a)')              '#'
 write(unt,'(a)')              '# a2F(e,eprime,w) function integrated over the FS. omega in a.u.'
 write(unt,'(2a)')             '# ngqpt: ',trim(ltoa(a2f%ngqpt))
 write(unt,'(a,i0,2(a,e16.6))')'# number of energies: ',a2f%nene," from enemin: ",a2f%enemin,' Ha with step ', a2f%deltaene
 write(unt,'(a,i0,3(a,e16.6))')'# number of frequencies: ',a2f%nomega," between omega_min: ",a2f%omega_min, &
                               ' Ha and omega_max: ',a2f%omega_max,' Ha with step:',a2f%wstep
 write(unt,'(a,e16.6)')         '#  the smearing width for gaussians is ',a2f%smear
 write(unt,'(a,e16.6)')"# Total DOS at Fermi level ",sum(a2f%n0)
 do spin=1,a2f%nsppol
   write(unt,"(a,i0,a,e16.6)")"# The DOS at Fermi level for spin ",spin," is ",a2f%n0(spin)
 end do
 do ii=1,2
   write(unt,'(a)')     "# "
 end do

end subroutine write_a2fw_header

end subroutine a2fw_ee_write
!!***

!----------------------------------------------------------------------

!!****f* m_phgamma/a2fw_solve_gap
!! NAME
!! a2fw_solve_gap
!!
!! FUNCTION
!!
!! INPUTS
!!  a2f<a2fw_t>=Container storing the Eliashberg functions.
!!  ntemp=Number of temperatures
!!  temp_range = min and max temperatures
!!  wcut = frequency cutoff for Matsubara sums
!!  mustar= mustar parameter
!!  nstep=Max number of SCF steps
!!  reltol = relative tolerance accepted for exit of main Eliashberg loop
!!  comm=MPI communicator
!!
!! TODO
!!  Implement python version in AbiPy
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine a2fw_solve_gap(a2f,cryst,ntemp,temp_range,wcut,mustar,nstep,reltol,prefix,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ntemp,nstep,comm
 real(dp),intent(in) :: wcut,mustar,reltol
 character(len=*),intent(in) :: prefix
 type(a2fw_t),intent(inout) :: a2f
 type(crystal_t),intent(in) :: cryst
!arrays
 real(dp),intent(in) :: temp_range(2)

!Local variables -------------------------
!scalars
 integer,parameter :: master = 0
 integer :: istep,ii,jj,it,spin,iwn,nwm,conv !iw,
 integer :: my_rank,nproc,ncid,ncerr
 real(dp) :: summ,kT,tstep,abs_delta,rel_delta,alpha,gap
 !character(len=500) :: msg
!arrays
 real(dp),allocatable :: wmts(:),lambda_ij(:,:),tmesh(:)
 real(dp),allocatable :: din(:),dout(:),zin(:),zout(:)

! *********************************************************************

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

 where (a2f%vals < zero) a2f%vals = zero

 ! Build linear mesh of temperatures.
 ABI_CHECK(ntemp > 1, "ntemp cannot be 1")
 tstep = (temp_range(2) - temp_range(1)) / (ntemp - 1)
 ABI_MALLOC(tmesh, (ntemp))
 tmesh = arth(temp_range(1), tstep, ntemp)

 ! Matsubara frequencies: i w_n = i (2n+1) pi T
 kT = kb_HaK * tmesh(1)
 nwm = 0
 do
   if ((2*nwm + 1) * pi * kT > wcut) exit
   nwm = nwm + 1
 end do
 ABI_CHECK(nwm /= 0, "Empy list of Matsubara frequencies, increase wcut")

#ifdef HAVE_NETCDF
 ! Open the netcdf file used to store the results of the calculation.
 if (my_rank == master) then
   NCF_CHECK(nctk_open_create(ncid, strcat(prefix, "_ELIASHBERG.nc"), xmpi_comm_self))
   NCF_CHECK(cryst%ncwrite(ncid))

   ! Define dimensions.
   ncerr = nctk_def_dims(ncid, [&
     nctkdim_t("maxnum_matsubara_frequencies", nwm), nctkdim_t("num_temperatures", ntemp) &
     ], defmode=.True.)
   NCF_CHECK(ncerr)

   ncerr = nctk_def_arrays(ncid, [&
     nctkarr_t('temperatures', "dp", "num_temperatures"),&
     nctkarr_t('delta_imag_axis', "dp", "maxnum_matsubara_frequencies, num_temperatures"),&
     nctkarr_t('zeta_imag_axis', "dp", "maxnum_matsubara_frequencies, num_temperatures"), &
     nctkarr_t('delta_real_axis', "dp", "maxnum_matsubara_frequencies, num_temperatures"),&
     nctkarr_t('zeta_real_axis', "dp", "maxnum_matsubara_frequencies, num_temperatures")  &
     ])
   NCF_CHECK(ncerr)

   NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "temperatures"), tmesh))
 end if
#endif

 do it=1,ntemp
   kT = kb_HaK * tmesh(it)

   ! Matsubara frequencies: i w_n = i (2n+1) pi T
   nwm = 0
   do
     if ((2*nwm + 1) * pi * kT > wcut) exit
     nwm = nwm + 1
   end do
   ABI_CHECK(nwm /= 0, "Empy list of Matsubara frequencies, increase wcut")
   write(std_out,*)"Number of matsubara frequencies",nwm

   ABI_MALLOC(wmts, (nwm))
   do ii=1,nwm
     wmts(ii) = (2*(ii-1) + 1) * pi * kT
   end do

   ! Compute lambda(n-n') kernel (symmetric)
   ABI_MALLOC(lambda_ij, (nwm, nwm))
   do spin=1,a2f%nsppol
     do jj=1,nwm
       do ii=1,jj
          lambda_ij(ii, jj) = a2fw_lambda_wij(a2f,  wmts(ii), wmts(jj), spin)
          if (ii /= jj) lambda_ij(jj, ii) = lambda_ij(ii, jj)
       end do
     end do
   end do
   !lambda_ij = lambda_ij * half

   ABI_MALLOC(din, (nwm))
   ABI_MALLOC(dout, (nwm))
   ABI_MALLOC(zin, (nwm))
   ABI_MALLOC(zout, (nwm))

   ! Initalize din
   if (it == 1) then
     din = 3.4*10-4 * eV_Ha
   else
     din = gap
   end if
   !dout = din

   conv = 0
   do istep=1,nstep
     !where (din < zero) din = zero
     do iwn=1,nwm
       summ = zero
       do jj=1,nwm
         summ = summ + wmts(jj) * lambda_ij(jj, iwn) / sqrt(wmts(jj)**2 + din(jj)**2)
       end do
       zout(iwn) = one + (pi * kT / wmts(iwn)) * summ
     end do

     do iwn=1,nwm
       summ = zero
       do jj=1,nwm
         summ = summ + din(jj) * (lambda_ij(jj, iwn) - mustar) / sqrt(wmts(jj)**2 + din(jj)**2)
       end do
       dout(iwn) = (pi * kT / zout(iwn)) * summ
     end do

     ! Test for convergence
     abs_delta = sum(abs(din))
     rel_delta = sum(abs(dout - din))
     !write(std_out,*)"rel_delta / abs_delta", (rel_delta / abs_delta)

     if ((rel_delta / abs_delta) < reltol) then
       conv = conv + 1
     else
       conv = 0
     end if
     if (conv == 2) exit

     ! TODO: Broyden mixing
     alpha = 0.2
     !alpha = 0.4
     !alpha = one
     din = alpha * dout + (one-alpha) * din
     zin = zout
   end do

   gap = dout(1)
   if (conv == 2) then
     write(std_out,*)"Converged at iteration: ",istep
   else
     write(std_out,*)"Not converged",rel_delta / abs_delta
   end if
   write(std_out,*)"T=",tmesh(it)," [K], gap ",gap*Ha_eV," [eV]"

   ! Write data to netcd file.
#ifdef HAVE_NETCDF
   if (my_rank == master) then
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "zeta_imag_axis"), zin, start=[1,it]))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "delta_imag_axis"), din, start=[1,it]))
   end if
#endif

   !if (it == 1) then
   !  do iwn=1,nwm
   !    write(999,*)wmts(iwn),dout(iwn),zout(iwn)
   !  end do
   !end if

   ABI_FREE(din)
   ABI_FREE(dout)
   ABI_FREE(zin)
   ABI_FREE(zout)

   ABI_FREE(wmts)
   ABI_FREE(lambda_ij)
 end do ! it

 ABI_FREE(tmesh)

#ifdef HAVE_NETCDF
 if (my_rank == master) then
   NCF_CHECK(nf90_close(ncid))
 end if
#endif

end subroutine a2fw_solve_gap
!!***

!----------------------------------------------------------------------

!!****f* m_phgamma/a2fw_tr_free
!! NAME
!! a2fw_tr_free
!!
!! FUNCTION
!!  Free the memory allocated in a2f
!!
!! SIDE EFFECTS
!!  a2f<a2fw_tr_t>=Structure storing the Eliashberg function a2F.
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine a2fw_tr_free(a2f_tr)

!Arguments ------------------------------------
 class(a2fw_tr_t),intent(inout) :: a2f_tr

! *********************************************************************

 ! @a2fw_tr_t

 ! integer
 ABI_SFREE(a2f_tr%qshift)

 ! real
 ABI_SFREE(a2f_tr%n0)
 ABI_SFREE(a2f_tr%omega)
 ABI_SFREE(a2f_tr%vals_in)
 ABI_SFREE(a2f_tr%vals_out)
 ABI_SFREE(a2f_tr%vals_tr)
 ABI_SFREE(a2f_tr%vals_tr_gen)
 ABI_SFREE(a2f_tr%lambdaw_tr)

end subroutine a2fw_tr_free
!!***

!----------------------------------------------------------------------

!!****f* m_phgamma/a2fw_tr_init
!! NAME
!! a2fw_tr_init
!!
!! FUNCTION
!!  Calculates the FS averaged alpha^2F_tr,in,out(w) functions
!!
!! INPUTS
!!  cryst<crystal_t>=Info on the unit cell.
!!  ifc<ifc_type>=Interatomic force constants.
!!  gams<phgamma_t>=Structure storing the phonon linewidths.
!!  wstep=Step for linear frequency mesh in Ha.
!!  wminmax(2)=Minimum and maximum phonon frequency. Used to construct the linear mesh for A2F(w).
!!  intmeth=Integration method: 1 for gaussian, 2 for tetrahedra
!!  smear=Gaussian broadening used to approximate the Dirac delta.
!!  ngqpt(3)=Divisions of the Q-mesh used for interpolating the phonon linewidths (see also nqshift and qshift).
!!  nqshift=Number of shifts used to generated the Q-mesh.
!!  qshift(3,nqshift)=The shifts.
!!  comm=MPI communicator
!!  [qintp]=If set to False, ngqgpt, nqshift, qshift and qptop  are ignored and
!!     A2F(w) is computed from the IBZ values stored in gams. Default: True i.e use Fourier interpolation.
!!  [qptopt]=Controls the generation of the q-points. If not specified, the routine takes fully into account
!!    the symmetries of the system to generate the q points in the IBZone i.e. qptopt=1
!!    Other values of qptopt can be used for debugging purpose.
!!
!! OUTPUT
!!  a2f_tr<a2fw_tr_t>=Structure storing the Eliashberg transport function a2F_tr(w).
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine a2fw_tr_init(a2f_tr, gams, cryst, ifc, intmeth, wstep, wminmax, smear, ngqpt, nqshift, qshift, comm,&
  qintp, qptopt) ! optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: intmeth,nqshift,comm
 integer,intent(in),optional :: qptopt
 real(dp),intent(in) :: wstep,smear
 logical,optional,intent(in) :: qintp
 type(phgamma_t),intent(inout) :: gams
 type(ifc_type),intent(in) :: ifc
 type(crystal_t),intent(in) :: cryst
 type(a2fw_tr_t),target,intent(out) :: a2f_tr
!arrays
 integer,intent(in) :: ngqpt(3)
 real(dp),intent(in) :: wminmax(2),qshift(3,nqshift)

!Local variables -------------------------
!scalars
 integer,parameter :: master = 0
 integer :: my_qptopt,iq_ibz,nqibz,ount,my_rank,nproc,cnt
 integer :: mu,iw,natom3,nsppol,spin,ierr,nomega,nqbz, idir, jdir
 real(dp) :: cpu, wall, gflops
 real(dp) :: omega,xx,omega_min,omega_max,ww
 logical :: do_qintp
 character(len=500) :: msg
 type(htetra_t) :: qtetra
!arrays
 integer :: qptrlatt(3,3),new_qptrlatt(3,3)
 real(dp),allocatable :: my_qshift(:,:)
 real(dp) :: lambda_iso(3,3), omega_log(3,3)
 real(dp) :: phfrq(gams%natom3)
 real(dp) :: gamma_in_ph(3,3,gams%natom3), gamma_out_ph(3,3,gams%natom3)
 real(dp) :: lambda_in_ph(3,3,gams%natom3), lambda_out_ph(3,3,gams%natom3)
 real(dp),allocatable :: tmp_a2f_in(:,:,:), tmp_a2f_out(:,:,:)
 real(dp), ABI_CONTIGUOUS pointer :: a2f_tr_1d(:)
 real(dp),allocatable :: qibz(:,:),wtq(:),qbz(:,:)
 real(dp),allocatable :: a2f_tr_1mom(:),a2f_tr_logmom(:),a2f_tr_logmom_int(:),wdt(:,:)
 real(dp),allocatable :: lambda_in_tetra(:,:,:,:,:),phfreq_tetra(:,:,:)
 real(dp),allocatable :: lambda_out_tetra(:,:,:,:,:)

! *********************************************************************

 !@a2fw_tr_t
 my_qptopt = 1; if (present(qptopt)) my_qptopt = qptopt
 do_qintp = .True.; if (present(qintp)) do_qintp = qintp

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)
 nsppol = gams%nsppol; natom3 = gams%natom3

 call cwtime(cpu, wall, gflops, "start")

 if (do_qintp) then
   ! Generate the q-mesh by finding the IBZ and the corresponding weights.
   qptrlatt = 0; qptrlatt(1,1) = ngqpt(1); qptrlatt(2,2) = ngqpt(2); qptrlatt(3,3) = ngqpt(3)

   call kpts_ibz_from_kptrlatt(cryst, qptrlatt, my_qptopt, nqshift, qshift, nqibz, qibz, wtq, nqbz, qbz, &
      new_kptrlatt=new_qptrlatt, new_shiftk=my_qshift)
   ABI_FREE(qbz)

   ! Store quantities that cannot be easily (and safely) calculated if we only know the IBZ.
   a2f_tr%ngqpt = ngqpt; a2f_tr%nqshift = size(my_qshift, dim=2)
   ABI_MALLOC(a2f_tr%qshift, (3, a2f_tr%nqshift))
   a2f_tr%qshift = my_qshift
   ABI_FREE(my_qshift)

 else
   ! No interpolation. Use q-mesh parameters from gams%
   a2f_tr%ngqpt = gams%ngqpt; a2f_tr%nqshift = 1
   nqibz = gams%nqibz
   ABI_MALLOC(qibz, (3,nqibz))
   ABI_MALLOC(wtq, (nqibz))
   ABI_MALLOC(a2f_tr%qshift, (3, a2f_tr%nqshift))
   qibz = gams%qibz; wtq = gams%wtq
   a2f_tr%qshift = zero ! Note: assuming q-mesh centered on Gamma.
 end if

 call cwtime_report(" a2fw_tr_init, q-setup", cpu, wall, gflops)

 ! Define Min and max frequency for the mesh (enlarge it a bit)
 omega_min = wminmax(1); omega_max = wminmax(2)
 omega_min = omega_min - 0.1*abs(omega_min)
 if (omega_min >= zero) omega_min = one/Ha_meV
 omega_max = omega_max + 0.1*abs(omega_max)

 a2f_tr%nsppol = nsppol; a2f_tr%natom3 = gams%natom3; a2f_tr%smear = smear
 a2f_tr%omega_min = omega_min; a2f_tr%omega_max = omega_max
 nomega = int((omega_max - omega_min) / wstep); a2f_tr%nomega = nomega; a2f_tr%wstep = wstep

 ABI_MALLOC(a2f_tr%n0, (nsppol))
 a2f_tr%n0 = gams%n0
 ! Build linear mesh.
 ABI_MALLOC(a2f_tr%omega, (nomega))
 a2f_tr%omega = arth(omega_min, wstep, nomega)
 ABI_CALLOC(a2f_tr%vals_in, (nomega, 3, 3, 0:natom3, nsppol))
 ABI_CALLOC(a2f_tr%vals_out, (nomega,3, 3, 0:natom3, nsppol))
 ABI_CALLOC(a2f_tr%vals_tr, (nomega,3, 3, 0:natom3, nsppol))
 ABI_CALLOC(a2f_tr%lambdaw_tr, (nomega, 3, 3, 0:natom3, nsppol))
 ABI_MALLOC(tmp_a2f_in, (nomega, 3, 3))
 ABI_MALLOC(tmp_a2f_out, (nomega, 3, 3))

 if (intmeth == 2) then
   call cwtime(cpu, wall, gflops, "start")

   ! Prepare tetrahedron integration.
   qptrlatt = 0; qptrlatt(1, 1) = a2f_tr%ngqpt(1); qptrlatt(2, 2) = a2f_tr%ngqpt(2); qptrlatt(3, 3) = a2f_tr%ngqpt(3)
   qtetra = tetra_from_kptrlatt(cryst, my_qptopt, qptrlatt, a2f_tr%nqshift, a2f_tr%qshift, nqibz, qibz, comm, msg, ierr)
   if (ierr/=0) MSG_ERROR(msg)

   ABI_MALLOC_OR_DIE(lambda_in_tetra, (nqibz, 3, 3, natom3, nsppol), ierr)
   ABI_MALLOC_OR_DIE(lambda_out_tetra, (nqibz, 3, 3, natom3, nsppol), ierr)
   lambda_in_tetra = zero
   lambda_out_tetra = zero

   ABI_MALLOC_OR_DIE(phfreq_tetra, (nqibz, natom3, nsppol), ierr)
   phfreq_tetra = zero

   call cwtime_report(" a2fw_tr_init%tetra,", cpu, wall, gflops)
 end if

 call cwtime(cpu, wall, gflops, "start")

 ! Loop over spins and qpoints in the IBZ
 cnt = 0
 do spin=1,nsppol
   do iq_ibz=1,nqibz
     cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle

     ! Interpolate or evaluate gamma directly.
     if (do_qintp) then
       call phgamma_vv_interp(gams,cryst,ifc,spin,qibz(:,iq_ibz),phfrq,gamma_in_ph,gamma_out_ph,lambda_in_ph,lambda_out_ph)
     else
       call phgamma_vv_eval_qibz(gams,cryst,ifc,iq_ibz,spin,phfrq,gamma_in_ph,gamma_out_ph,lambda_in_ph,lambda_out_ph)
     end if

     select case (intmeth)
     case (1)
       ! Gaussian: Add all contributions from the phonon modes at this qpoint to a2f_tr
       ! (note that unstable modes are included).
       do mu=1,natom3
         tmp_a2f_in = zero
         tmp_a2f_out = zero
         do iw=1,nomega
           xx = a2f_tr%omega(iw) - phfrq(mu)
           tmp_a2f_in(iw,:,:)  = tmp_a2f_in(iw,:,:)  + gaussian(xx, smear) * lambda_in_ph(:,:,mu)  * abs(phfrq(mu))
           tmp_a2f_out(iw,:,:) = tmp_a2f_out(iw,:,:) + gaussian(xx, smear) * lambda_out_ph(:,:,mu) * abs(phfrq(mu))
         end do
         a2f_tr%vals_in(:,:,:,mu,spin)  = a2f_tr%vals_in(:,:,:,mu,spin)  + tmp_a2f_in(:,:,:) * wtq(iq_ibz)
         a2f_tr%vals_out(:,:,:,mu,spin) = a2f_tr%vals_out(:,:,:,mu,spin) + tmp_a2f_out(:,:,:) * wtq(iq_ibz)
       end do

     case (2)
       ! Tetra: store data.
       do mu=1,natom3
         lambda_in_tetra(iq_ibz, :,:, mu, spin) = lambda_in_ph(:,:,mu) * abs(phfrq(mu))
         lambda_out_tetra(iq_ibz, :,:, mu, spin) = lambda_out_ph(:,:,mu) * abs(phfrq(mu))
         phfreq_tetra(iq_ibz, mu, spin) = phfrq(mu)
       end do
     end select

   end do ! iq_ibz
 end do ! spin

 if (intmeth == 2) then
   ! Collect results on each node.
   call xmpi_sum(lambda_in_tetra, comm, ierr)
   call xmpi_sum(lambda_out_tetra, comm, ierr)
   call xmpi_sum(phfreq_tetra, comm, ierr)

   ! workspace for tetra.
   ABI_MALLOC(wdt, (nomega, 2))

   ! TODO: with the tetra_get_onewk call we can integrate this above
   ! and avoid allocating all of lambda_in_tetra and phfreq_tetra!!!
   ! For each mode get its contribution
   cnt = 0
   do spin=1,nsppol
     do mu=1,natom3
       do iq_ibz=1,nqibz
         cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle ! mpi-parallelism

         call qtetra%get_onewk(iq_ibz, gams%bcorr, nomega, nqibz, phfreq_tetra(:,mu,spin), omega_min, omega_max, one, wdt)
         wdt = wdt * wtq(iq_ibz)

         ! Accumulate (Integral of a2F_tr is computed afterwards)
         do idir=1,3
           do jdir=1,3
             a2f_tr%vals_in(:,idir,jdir,mu,spin)  = a2f_tr%vals_in(:,idir,jdir,mu,spin)  &
                + wdt(:,1) * lambda_in_tetra(iq_ibz,idir,jdir,mu,spin)
             a2f_tr%vals_out(:,idir,jdir,mu,spin) = a2f_tr%vals_out(:,idir,jdir,mu,spin) &
                + wdt(:,1) * lambda_out_tetra(iq_ibz,idir,jdir,mu,spin)
             !a2f_tr%lambdaw(:,mu,spin) = a2f_tr%lambdaw(:,mu,spin) + wdt(:,2) * lambda_XX_tetra(iq_ibz, mu, spin)
           end do
         end do
       end do
     end do
   end do

   ! Free memory allocated for qtetra.
   ABI_FREE(wdt)
   ABI_FREE(lambda_in_tetra)
   ABI_FREE(lambda_out_tetra)
   ABI_FREE(phfreq_tetra)
   call qtetra%free()
 end if

 ! Collect final results on each node and divide by g(eF, spin)
 ! TODO: check whether this has to be reinstated, in particular for a2f_tr
 call xmpi_sum(a2f_tr%vals_in, comm, ierr)
 call xmpi_sum(a2f_tr%vals_out, comm, ierr)

 ! TODO: check normalization with N(0) and <v^2> from Savrasov eq 20
 do spin=1,nsppol
   a2f_tr%vals_in(:,:,:,0,spin)  = sum(a2f_tr%vals_in(:,:,:,1:natom3,spin), dim=4)
   a2f_tr%vals_out(:,:,:,0,spin) = sum(a2f_tr%vals_out(:,:,:,1:natom3,spin), dim=4)
   !a2f_tr%vals(:,:,spin) = a2f_tr%vals(:,:,spin) / (two_pi*a2f_tr%n0(spin))
 end do

 a2f_tr%vals_tr = a2f_tr%vals_out - a2f_tr%vals_in

 call cwtime_report(" a2fw_tr_init, a2f_eval", cpu, wall, gflops)

 ! NOTE: for the moment calculate the log moms etc on just the trace of the a2f functions
 ABI_MALLOC(a2f_tr_1mom, (nomega))
 ABI_MALLOC(a2f_tr_logmom, (nomega))
 ABI_MALLOC(a2f_tr_logmom_int, (nomega))

 !call a2fw_tr_print_info()

 ! Compute lambda(w) (resolved in spin mode and directions)
 do spin=1,nsppol
   do mu=0,natom3
     do idir=1,3
       do jdir=1,3

         do iw=1,nomega
           ww = a2f_tr%omega(iw)
           if (abs(ww) > EPHTK_WTOL) then
              a2f_tr_1mom(iw)  = a2f_tr%vals_tr(iw,idir,jdir,mu,spin) / abs(ww)
           else
              a2f_tr_1mom(iw) = zero
           end if
           ! TODO: Strange that the first value of int(f) is not set to zero!
           ! FIXME: The frequency integration overestimates lambda if there are negative phonon frequencies!
         end do

         call simpson_int(nomega, wstep, a2f_tr_1mom, a2f_tr%lambdaw_tr(:,idir,jdir,mu,spin))
       end do ! jdir
     end do ! idir
   end do ! mu
 end do ! spin

 if (my_rank == master) then
   ount = std_out
   if (nsppol > 1) then
     write(msg,'(3a)') ch10,'Comment: some of the following quantities should be integrated over spin', ch10
     call wrtout(ount, msg)
   end if
 end if

 do spin=1,nsppol
   lambda_iso = a2fw_tr_moment(a2f_tr, 0, spin)

   do idir=1,3
     do jdir=1,3
       a2f_tr_1d => a2f_tr%vals_tr(:,idir,jdir,0,spin)

       ! Get log moment of alpha^2F.
       a2f_tr_logmom = zero
       do iw=1,nomega
         omega = a2f_tr%omega(iw)
         if (abs(omega) > EPHTK_WTOL .and. abs(lambda_iso(idir,jdir)) > EPHTK_WTOL) then
           a2f_tr_logmom(iw) = (two/lambda_iso(idir,jdir)) * a2f_tr_1d(iw)*log(abs(omega))/abs(omega)
         end if
       end do
       call simpson_int(nomega, wstep, a2f_tr_logmom, a2f_tr_logmom_int)
       !write(std_out,*)' iw,nomega,greatest_real,a2f_tr_logmom_int(nomega)=',& iw,nomega,greatest_real,a2f_tr_logmom_int(nomega)
       if(abs(a2f_tr_logmom_int(nomega)) < log(greatest_real*tol6)) then
         omega_log(idir,jdir) = exp(a2f_tr_logmom_int(nomega))
       else
         omega_log(idir,jdir)=greatest_real*tol6
       endif
     end do
   end do

   ! TODO: make output only for irred values xx yy zz and top half of matrix
   if (my_rank == master) then
     write(ount,'(a)')' Evaluation of parameters analogous to electron-phonon coupling for 3x3 directions '
     write(ount,'(a,3(3es10.3,2x))') ' lambda = ',lambda_iso
     write(ount,'(a,3(3es10.3,2x),a)' )' omegalog  = ',omega_log,' (Ha) '
     write(ount,'(a,3(3es10.3,2x),a)' )'             ',omega_log*Ha_K, ' (Kelvin) '
     write(ount,"(a)")'    positive moments of alpha2Ftr:'
     write(ount,'(a,9(es10.3,2x))' )' lambda <omega^2> = ',a2fw_tr_moment(a2f_tr, 2, spin)
     write(ount,'(a,9(es10.3,2x))' )' lambda <omega^3> = ',a2fw_tr_moment(a2f_tr, 3, spin)
     write(ount,'(a,9(es10.3,2x))' )' lambda <omega^4> = ',a2fw_tr_moment(a2f_tr, 4, spin)
     write(ount,'(a,9(es10.3,2x))' )' lambda <omega^5> = ',a2fw_tr_moment(a2f_tr, 5, spin)
   end if
 end do

 ABI_FREE(tmp_a2f_in)
 ABI_FREE(tmp_a2f_out)
 ABI_FREE(a2f_tr_1mom)
 ABI_FREE(a2f_tr_logmom)
 ABI_FREE(a2f_tr_logmom_int)
 ABI_FREE(qibz)
 ABI_FREE(wtq)

end subroutine a2fw_tr_init
!!***

!----------------------------------------------------------------------

!!****f* m_phgamma/a2fw_tr_write
!! NAME
!! a2fw_tr_write
!!
!! FUNCTION
!!  Write alpha^2F_tr(w) to an external file in text form
!!
!! INPUTS
!!  a2f_tr<a2fw_tr_t>=Container storing the Eliashberg transport functions.
!!  basename=Filename for output.
!!  post=String appended to netcdf variables e.g. _qcoarse, _qintp
!!  ncid=Netcdf file handler. Set it to nctk_noid to disable output.
!!
!! OUTPUT
!!  Output is written to file. This routine should be called by one MPI proc.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine a2fw_tr_write(a2f_tr, basename, post, ncid)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 character(len=*),intent(in) :: basename, post
 class(a2fw_tr_t),intent(in) :: a2f_tr

!Local variables -------------------------
!scalars
 integer :: iw,spin,unt,ii,mu,idir, jdir
#ifdef HAVE_NETCDF
 integer :: ncerr
 character(len=500) :: dim1_name
#endif
 character(len=500) :: msg
 character(len=fnlen) :: path

! *********************************************************************

 ! Write spin-resolved a2F_tr(w)
 path = strcat(basename, "_A2FW_tr")
 if (open_file(path, msg, newunit=unt, form="formatted", action="write", status="unknown") /= 0) then
   MSG_ERROR(msg)
 end if

 call write_a2fw_tr_header()

 if (a2f_tr%nsppol == 1) then
   write(unt,'(a)')"# Frequency, a2F_tr(w,i,j)"
   do iw=1,a2f_tr%nomega
     write(unt,'(e16.6, 2x, 3(3e16.6,1x))') a2f_tr%omega(iw), a2f_tr%vals_tr(iw,:,:,0,1)
   end do
   write(unt,'(3a)')ch10, ch10, "# Frequency, lambda(w,i,j)"
   do iw=1,a2f_tr%nomega
     write(unt,'(e16.6, 2x, 3(3e16.6,1x))') a2f_tr%omega(iw), a2f_tr%lambdaw_tr(iw,:,:,0,1)
   end do

 else
   write(unt,'(a)')"# Frequency, a2F_tr_tot(w,i,j), a2F_tr_spin1(w,i,j) ..."
   do iw=1,a2f_tr%nomega
     write(unt,'(e16.6, 2x, 3(3e16.6,1x), 2x,  3(3e16.6,1x), 2x,  3(3e16.6,1x))') a2f_tr%omega(iw), &
                  ((sum(a2f_tr%vals_tr(iw,idir,jdir,0,:)), idir=1,3), jdir=1,3) , &
                  a2f_tr%vals_tr(iw,:,:,0,1), &  ! UP
                  a2f_tr%vals_tr(iw,:,:,0,2)     ! DOWN
   end do
   write(unt,'(3a)') ch10, ch10, "# Frequency, lambda_tr_tot(w,i,j)dw, lambda_tr_spin1(w,i,j) ..."
   do iw=1,a2f_tr%nomega
     write(unt,'(e16.6, 2x, 3(3e16.6,1x), 2x, 3(3e16.6,1x), 2x,      3(3e16.6,1x), 2x,  3(3e16.6,1x))') a2f_tr%omega(iw), &
                  ((sum(a2f_tr%lambdaw_tr(iw,idir,jdir,0,:)), idir=1,3), jdir=1,3) , &  ! TOT
                  a2f_tr%lambdaw_tr(iw,:,:,0,1), &  ! UP
                  a2f_tr%lambdaw_tr(iw,:,:,0,2)     ! DOWN
   end do
 end if

 close(unt)

 ! Write phonon mode contributions to a2F_tr(w,i,j)
 path = strcat(basename, "_PH_A2FW_tr")
 if (open_file(path, msg, newunit=unt, form="formatted", action="write", status="unknown") /= 0) then
   MSG_ERROR(msg)
 end if
 call write_a2fw_tr_header()

 if (a2f_tr%nsppol == 1) then
   do mu=0,a2f_tr%natom3
     write(unt,'(a,i0)')"# Phonon mode ",mu
     write(unt,'(a)')"# Frequency, a2F_tr(w)"
     do iw=1,a2f_tr%nomega
       write(unt,'(e16.6,2x,3(3e16.6,1x))') a2f_tr%omega(iw), a2f_tr%vals_tr(iw,:,:,mu,1)
     end do
     do ii=1,2; write(unt,'(a)')""; end do
     write(unt,'(a)')"# Frequency, lambda_tr(w)"
     do iw=1,a2f_tr%nomega
       write(unt,'(e16.6,2x,3(3e16.6,1x))') a2f_tr%omega(iw), a2f_tr%lambdaw_tr(iw,:,:,mu,1)
     end do
     do ii=1,2; write(unt,'(a)')""; end do
   end do

 else
   do mu=0,a2f_tr%natom3
     write(unt,'(a,i0)')"# Phonon mode ",mu
     write(unt,'(a)')"# Frequency, a2F_tot(w), a2F_spin1(w) ..."
     do iw=1,a2f_tr%nomega
       write(unt,*) a2f_tr%omega(iw), &
                    sum(a2f_tr%vals_tr(iw,:,:,mu,:)), &  ! TOT
                    a2f_tr%vals_tr(iw,:,:,mu,1)     , &  ! UP
                    a2f_tr%vals_tr(iw,:,:,mu,2)          ! DOWN
     end do
     write(unt,'(3a)')ch10,ch10,"# Frequency, lambda_tot(w)dw, lambda_spin1(w) ..."
     do iw=1,a2f_tr%nomega
       write(unt,*) a2f_tr%omega(iw), &
                     sum(a2f_tr%lambdaw_tr(iw,:,:,mu,:)), &  ! TOT
                     a2f_tr%lambdaw_tr(iw,:,:,mu,1),      &  ! UP
                     a2f_tr%lambdaw_tr(iw,:,:,mu,2)          ! DOWN
     end do
   end do
 end if

 close(unt)

 if (ncid /= nctk_noid) then
#ifdef HAVE_NETCDF
   ! Define dimensions.
   dim1_name = strcat("a2ftr_nomega", post)
   ncerr = nctk_def_dims(ncid, [ &
     nctkdim_t(dim1_name, a2f_tr%nomega), &
     nctkdim_t("natom3p1", a2f_tr%natom3 + 1), nctkdim_t("number_of_spins", a2f_tr%nsppol) &
   ], defmode=.True.)
   NCF_CHECK(ncerr)

   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t(strcat('a2ftr_mesh', post), "dp", dim1_name), &
     nctkarr_t(strcat('a2ftr_values', post), "dp", strcat(dim1_name, ", three, three, natom3p1, number_of_spins")), &
     nctkarr_t(strcat('a2ftr_lambdaw', post), "dp", strcat(dim1_name, ", three, three, natom3p1, number_of_spins")) &
   ])
   NCF_CHECK(ncerr)

   NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, strcat("a2ftr_mesh", post)), a2f_tr%omega))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, strcat("a2ftr_values", post)), a2f_tr%vals_tr))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, strcat("a2ftr_lambdaw", post)), a2f_tr%lambdaw_tr))
#endif
 end if

contains

subroutine write_a2fw_tr_header()

 ! Output the header.
 write(unt,'(a)')              '#'
 write(unt,'(a)')              '# ABINIT package: a2F_tr(w) file'
 write(unt,'(a)')              '#'
 write(unt,'(a)')              '# a2F_tr(w) function integrated over the FS. omega in a.u.'
 write(unt,'(2a)')             '# ngqpt: ',trim(ltoa(a2f_tr%ngqpt))
 write(unt,'(a,i0,3(a,e16.6))')'# number of frequencies: ',a2f_tr%nomega," between omega_min: ",a2f_tr%omega_min,&
                               ' Ha and omega_max: ',a2f_tr%omega_max,' Ha with step:',a2f_tr%wstep
 write(unt,'(a,e16.6)')         '#  the smearing width for gaussians is ',a2f_tr%smear
 write(unt,'(a,e16.6)')"# Total DOS at Fermi level ",sum(a2f_tr%n0)
 do spin=1,a2f_tr%nsppol
   write(unt,"(a,i0,a,e16.6)")"# The DOS at Fermi level for spin ",spin," is ",a2f_tr%n0(spin)
 end do
 do ii=1,2; write(unt,'(a)')     "# "; end do

end subroutine write_a2fw_tr_header

end subroutine a2fw_tr_write
!!***

!----------------------------------------------------------------------

!!****f* m_phgamma/eph_phgamma
!! NAME
!!  eph_phgamma
!!
!! FUNCTION
!!  Compute phonon linewidths in metals.
!!
!! INPUTS
!! wk0_path=String with the path to the GS unperturbed WFK file.
!! ngfft(18),ngfftf(18)=Coarse and Fine FFT meshes.
!! dtset<dataset_type>=All input variables for this dataset.
!! ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!! dvdb<dbdb_type>=Database with the DFPT SCF potentials.
!! ifc<ifc_type>=interatomic force constants and corresponding real space grid info.
!! pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!! pawang<pawang_type)>=PAW angular mesh and related data.
!! pawrad(ntypat*usepaw)<pawrad_type>=Paw radial mesh and related data.
!! pawtab(ntypat*usepaw)<pawtab_type>=Paw tabulated starting data.
!! psps<pseudopotential_type>=Variables related to pseudopotentials.
!! comm=MPI communicator.
!!
!! OUTPUT
!!
!! PARENTS
!!      eph
!!
!! CHILDREN
!!
!! SOURCE

subroutine eph_phgamma(wfk0_path, dtfil, ngfft, ngfftf, dtset, cryst, ebands, dvdb, ifc, &
                       pawfgr, pawang, pawrad, pawtab, psps, mpi_enreg, comm)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: wfk0_path
 integer,intent(in) :: comm
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 type(dvdb_t),intent(inout) :: dvdb
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
 type(pawfgr_type),intent(in) :: pawfgr
 type(ifc_type),intent(in) :: ifc
 type(mpi_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18),ngfftf(18)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: tim_getgh1c = 1, berryopt0 = 0, timrev1 = 1
 integer,parameter :: useylmgr = 0, useylmgr1 = 0, master = 0, ndat1 = 1, eph_scalprod0 = 0
 integer :: my_rank,nproc,mband,nsppol,nkibz,idir,ipert,iq_ibz
 integer :: cplex,db_iqpt,natom,natom3,ipc,ipc1,ipc2,nspinor,onpw
 integer :: bstart_k,bstart_kq,nband_k,nband_kq,ib1,ib2, band_k, band_kq
 integer :: ik_ibz,ik_bz,ikq_bz,ikq_ibz,isym_k,isym_kq,trev_k,trev_kq,timerev_q
 integer :: ik_fs, myik, mys !, myiq
 integer :: spin,istwf_k,istwf_kq,npw_k,npw_kq
 integer :: ii,ipw,mpw,my_mpw,mnb,ierr,cnt,ncid
 integer :: n1,n2,n3,n4,n5,n6,nspden,do_ftv1q, ltetra
 integer :: sij_opt,usecprj,usevnl,optlocal,optnl,opt_gvnlx1
 integer :: nfft,nfftf,mgfft,mgfftf,kqcount,nkpg,nkpg1,edos_intmeth
 integer :: jene, iene, comm_rpt, method, nesting
 integer :: my_npert, imyp, imyq
#ifdef HAVE_NETCDF
 integer :: ncerr
#endif
 real(dp) :: cpu, wall, gflops, cpu_q, wall_q, gflops_q, cpu_k, wall_k, gflops_k, cpu_all, wall_all, gflops_all
 real(dp) :: edos_step, edos_broad, sigma
 real(dp) :: ecut,eshift,eig0nk,dksqmax
 logical :: gen_eigenpb, need_velocities
 type(wfd_t) :: wfd
 type(fstab_t),pointer :: fs
 type(gs_hamiltonian_type) :: gs_hamkq
 type(rf_hamiltonian_type) :: rf_hamkq
 type(edos_t) :: edos
 type(phgamma_t) :: gams
 type(a2fw_t) :: a2fw
 type(a2fw_tr_t) :: a2fw_tr
 type(ddkop_t) :: ddkop
 type(xcomm_t) :: pert_comm, qs_comm, qpt_comm, bsum_comm, kpt_comm, spin_comm, pkb_comm !, ncwrite_comm
 character(len=500) :: msg
 character(len=fnlen) :: path
!arrays
 integer :: g0_k(3),g0bz_kq(3),g0_kq(3),symq(4,2,cryst%nsym)
 integer :: work_ngfft(18),gmax(3),my_gmax(3),gamma_ngqpt(3) !g0ibz_kq(3),
 integer :: indkk_kq(1,6)
 integer,allocatable :: kg_k(:,:),kg_kq(:,:),gtmp(:,:),nband(:,:),wfd_istwfk(:)
 integer,allocatable :: my_pinfo(:,:), pert_table(:,:)
 real(dp) :: kk(3),kq(3),kk_ibz(3),kq_ibz(3),qpt(3), lf(2),rg(2),res(2)
 real(dp) :: wminmax(2), n0(ebands%nsppol)
 real(dp) :: resvv_in(2,9), resvv_out(2,9), phfrq(3*cryst%natom)
 real(dp),allocatable :: displ_cart(:,:,:,:),displ_red(:,:,:,:)
 real(dp),allocatable :: grad_berry(:,:),kinpw1(:),kpg1_k(:,:),kpg_k(:,:),dkinpw(:)
 real(dp),allocatable :: ffnlk(:,:,:,:),ffnl1(:,:,:,:),ph3d(:,:,:),ph3d1(:,:,:)
 real(dp),allocatable :: v1scf(:,:,:,:),tgam(:,:,:),gkk_atm(:,:,:,:)
 real(dp),allocatable :: bras_kq(:,:,:),kets_k(:,:,:),h1kets_kq(:,:,:)
 real(dp),allocatable :: ph1d(:,:),vlocal(:,:,:,:),vlocal1(:,:,:,:,:)
 real(dp),allocatable :: ylm_kq(:,:),ylm_k(:,:),ylmgr_kq(:,:,:)
 real(dp),allocatable :: dummy_vtrial(:,:),gvnlx1(:,:),work(:,:,:,:)
 real(dp),allocatable :: gs1c(:,:), v1_work(:,:,:,:)
 real(dp),allocatable :: wt_ek(:,:),wt_ekq(:,:), dbldelta_wts(:,:)
 real(dp), allocatable :: tgamvv_in(:,:,:,:),  vv_kk(:,:,:)
 real(dp), allocatable :: tgamvv_out(:,:,:,:), vv_kkq(:,:,:)
 real(dp), allocatable :: tmp_vals_ee(:,:,:,:,:), emesh(:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:)
 type(fstab_t),target,allocatable :: fstab(:)
 type(pawcprj_type),allocatable  :: cwaveprj0(:,:)
 ! for the Eliashberg solver
 !integer :: ntemp
 !real(dp) :: reltol,wcut
 !real(dp) :: temp_range(2)
#ifdef HAVE_MPI
 integer :: ndims, comm_cart, me_cart
 logical :: reorder
 integer :: coords(5)
 integer,allocatable :: dims(:)
 logical,allocatable :: periods(:), keepdim(:)
#endif

!************************************************************************

 if (psps%usepaw == 1) then
   MSG_ERROR("PAW not implemented")
   ABI_UNUSED((/pawang%nsym, pawrad(1)%mesh_size/))
 end if

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)
 call cwtime(cpu_all, wall_all, gflops_all, "start")

 ! Copy important dimensions
 natom = cryst%natom; natom3 = 3 * natom; nsppol = ebands%nsppol; nspinor = ebands%nspinor; nspden = dtset%nspden
 nkibz = ebands%nkpt; mband = ebands%mband

 ! FFT meshes
 nfftf = product(ngfftf(1:3)); mgfftf = maxval(ngfftf(1:3))
 nfft = product(ngfft(1:3)) ; mgfft = maxval(ngfft(1:3))
 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3)
 n4 = ngfft(4); n5 = ngfft(5); n6 = ngfft(6)

 ! Compute electron DOS.
 edos_intmeth = 2; if (dtset%prtdos == 1) edos_intmeth = 1
 edos_step = dtset%dosdeltae; edos_broad = dtset%tsmear
 edos_step = 0.01 * eV_Ha; edos_broad = 0.3 * eV_Ha
 edos = ebands_get_edos(ebands, cryst, edos_intmeth, edos_step, edos_broad, comm)

 ! Store DOS per spin channels
 n0(:) = edos%gef(1:edos%nsppol)
 if (my_rank == master) then
   call edos%print(unit=ab_out)
   path = strcat(dtfil%filnam_ds(4), "_EDOS")
   call wrtout(ab_out, sjoin("- Writing electron DOS to file:", path, ch10))
   call edos%write(path)
 end if

 ! Find Fermi surface k points
 ! TODO: support kptopt, change setup of k-points if tetra: fist tetra weights then k-points on the Fermi surface!
 ABI_MALLOC(fstab, (nsppol))
 call fstab_init(fstab, ebands, cryst, dtset, comm)
 if (my_rank == master) then
   call fstab_print(fstab, unit=std_out)
   call fstab_print(fstab, unit=ab_out)
 end if

 ! Define q-mesh. eph_ngqpt_fine activates interpolation of potentials.
 gamma_ngqpt = ifc%ngqpt; if (all(dtset%eph_ngqpt_fine /= 0)) gamma_ngqpt = dtset%eph_ngqpt_fine

 call phgamma_init(gams, cryst, ifc, fstab(1), dtset, eph_scalprod0, gamma_ngqpt, n0, comm)
 call wrtout(std_out, sjoin("Will compute", itoa(gams%nqibz), "q-points in the IBZ"))

 ltetra = 2 !; ltetra = dtset%useria + 1
 ! Compute optimal energy window for tetra.
 !call find_ewin(gams%nqibz, gams%qibz, cryst, ebands, ltetra, sigma, comm)

 ! ========================
 ! === MPI DISTRIBUTION ===
 ! ========================
 ! Init for sequential execution.
 my_npert = natom3

 if (any(dtset%eph_np_pqbks /= 0)) then
   ! Use parameters from input file.
   pert_comm%nproc = dtset%eph_np_pqbks(1)
   qpt_comm%nproc  = dtset%eph_np_pqbks(2)
   bsum_comm%nproc = dtset%eph_np_pqbks(3)
   kpt_comm%nproc = dtset%eph_np_pqbks(4)
   spin_comm%nproc = dtset%eph_np_pqbks(5)
   my_npert = natom3 / pert_comm%nproc
   ABI_CHECK(my_npert > 0, "pert_comm_nproc cannot be greater than 3 * natom.")
   ABI_CHECK(mod(natom3, pert_comm%nproc) == 0, "pert_comm_nproc must divide 3 * natom.")
 else
   ! Automatic grid generation over k-points and spins.
   if (nsppol == 2 .and. mod(nproc, 2) == 0) then
     spin_comm%nproc = 2
     kpt_comm%nproc = nproc / 2
   else
     kpt_comm%nproc = nproc
   end if
 end if

 ! Consistency check.
 if (pert_comm%nproc * qpt_comm%nproc * bsum_comm%nproc * kpt_comm%nproc * spin_comm%nproc /= nproc) then
   write(msg, "(a,i0,3a, 6(a,1x,i0))") &
     "Cannot create 5d Cartesian grid with total nproc: ", nproc, ch10, &
     "Idle processes are not supported. The product of the `nproc_*` vars should be equal to nproc.", ch10, &
     "pert_nproc (", pert_comm%nproc, ") x qpt_nproc (", qpt_comm%nproc, ") x bsum_nproc (", bsum_comm%nproc, &
     ") x kcalc_nproc (", kpt_comm%nproc, ") x spin_nproc (", spin_comm%nproc, ") != ", &
     pert_comm%nproc * qpt_comm%nproc * bsum_comm%nproc * kpt_comm%nproc * spin_comm%nproc
   MSG_ERROR(msg)
 end if

 ABI_CHECK(bsum_comm%nproc == 1, "Band parallelism not implemented in m_phgamma")

#ifdef HAVE_MPI
 ! Create 5d cartesian communicator: 3*natom perturbations, q-points in IBZ, bands in FS, kpoints in FS, spins
 ndims = 5
 ABI_MALLOC(dims, (ndims))
 ABI_MALLOC(periods, (ndims))
 ABI_MALLOC(keepdim, (ndims))
 periods(:) = .False.; reorder = .False.
 dims = [pert_comm%nproc, qpt_comm%nproc, bsum_comm%nproc, kpt_comm%nproc, spin_comm%nproc]

 call MPI_CART_CREATE(comm, ndims, dims, periods, reorder, comm_cart, ierr)
 ! Find the index and coordinates of the current processor
 call MPI_COMM_RANK(comm_cart, me_cart, ierr)
 call MPI_CART_COORDS(comm_cart, me_cart, ndims, coords, ierr)

 ! Create communicator to distribute natom3 perturbations.
 keepdim = .False.; keepdim(1) = .True.
 call MPI_CART_SUB(comm_cart, keepdim, pert_comm%value, ierr); pert_comm%me = xmpi_comm_rank(pert_comm%value)

 ! Create communicator for qpoints in self-energy integration.
 keepdim = .False.; keepdim(2) = .True.
 call MPI_CART_SUB(comm_cart, keepdim, qpt_comm%value, ierr); qpt_comm%me = xmpi_comm_rank(qpt_comm%value)

 ! Create communicator for bands for self-energy summation
 keepdim = .False.; keepdim(3) = .True.
 call MPI_CART_SUB(comm_cart, keepdim, bsum_comm%value, ierr); bsum_comm%me = xmpi_comm_rank(bsum_comm%value)

 ! Create communicator for kpoints.
 keepdim = .False.; keepdim(4) = .True.
 call MPI_CART_SUB(comm_cart, keepdim, kpt_comm%value, ierr); kpt_comm%me = xmpi_comm_rank(kpt_comm%value)

 ! Create communicator for spins.
 keepdim = .False.; keepdim(5) = .True.
 call MPI_CART_SUB(comm_cart, keepdim, spin_comm%value, ierr); spin_comm%me = xmpi_comm_rank(spin_comm%value)

 ! Create communicator for the (qpoint, spin) loops
 keepdim = .False.; keepdim(2) = .True.; keepdim(5) = .True.
 call MPI_CART_SUB(comm_cart, keepdim, qs_comm%value, ierr); qs_comm%me = xmpi_comm_rank(qs_comm%value)

 ! Create communicator for the (perturbation, k-point, band_sum)
 keepdim = .False.; keepdim(1) = .True.; keepdim(3:4) = .True.
 call MPI_CART_SUB(comm_cart, keepdim, pkb_comm%value, ierr); pkb_comm%me = xmpi_comm_rank(pkb_comm%value)

 ABI_FREE(dims)
 ABI_FREE(periods)
 ABI_FREE(keepdim)
 call xmpi_comm_free(comm_cart)
#endif

 ! Distribute k-points and create mapping to ikcalc index.
 call xmpi_split_cyclic(nsppol, spin_comm%value, gams%my_nspins, gams%my_spins)
 ABI_CHECK(gams%my_nspins > 0, sjoin("nsppol (", itoa(nsppol), ") < spin_comm_nproc (", itoa(spin_comm%nproc), ")"))

 call xmpi_split_cyclic(gams%nqibz, qpt_comm%value, gams%my_nqibz, gams%my_iqibz)
 ABI_CHECK(gams%my_nqibz > 0, sjoin("nqibz (", itoa(gams%nqibz), ") < qpt_comm_nproc (", itoa(qpt_comm%nproc), ")"))

 path = strcat(dtfil%filnam_ds(4), "_A2F.nc")
 ncid = nctk_noid
#ifdef HAVE_NETCDF
 if (my_rank == master) then
   ! Master creates the netcdf file used to store the results of the calculation.
   NCF_CHECK(nctk_open_create(ncid, path, xmpi_comm_self))
   NCF_CHECK(cryst%ncwrite(ncid))
   NCF_CHECK(ebands_ncwrite(ebands, ncid))
   NCF_CHECK(edos%ncwrite(ncid))

   ! Add dimensions.
   ncerr = nctk_def_dims(ncid, [nctkdim_t("nqibz", gams%nqibz), nctkdim_t("natom3", natom3)], defmode=.True.)
   NCF_CHECK(ncerr)
   ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "eph_intmeth", "eph_transport", "symdynmat"])
   NCF_CHECK(ncerr)
   ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "ph_intmeth"])
   NCF_CHECK(ncerr)
   ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "eph_fsewin", "eph_fsmear", "eph_extrael", "eph_fermie"])
   NCF_CHECK(ncerr)
   ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "ph_wstep", "ph_smear"])
   NCF_CHECK(ncerr)

   ! Define arrays for results.
   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("ngqpt", "int", "three"), &
     nctkarr_t("eph_ngqpt_fine", "int", "three"), &
     nctkarr_t("ddb_ngqpt", "int", "three"), &
     nctkarr_t("ph_ngqpt", "int", "three"), &
     nctkarr_t('qibz', "dp", "number_of_reduced_dimensions, nqibz"), &
     nctkarr_t('wtq', "dp", "nqibz"), &
     nctkarr_t('phfreq_qibz', "dp", "natom3, nqibz"), &
     nctkarr_t('phdispl_cart_qibz', "dp", "two, natom3, natom3, nqibz"), &
     nctkarr_t('phgamma_qibz', "dp", "natom3, nqibz, number_of_spins"), &
     nctkarr_t('phlambda_qibz', "dp", "natom3, nqibz, number_of_spins") &
   ])
   NCF_CHECK(ncerr)

   ! ======================================================
   ! Write data that do not depend on the (kpt, spin) loop.
   ! ======================================================
   NCF_CHECK(nctk_set_datamode(ncid))

   ncerr = nctk_write_iscalars(ncid, &
       [character(len=nctk_slen) :: "eph_intmeth", "eph_transport", "symdynmat", "ph_intmeth"], &
       [dtset%eph_intmeth, dtset%eph_transport, dtset%symdynmat, dtset%ph_intmeth])
   NCF_CHECK(ncerr)
   ncerr = nctk_write_dpscalars(ncid, &
     [character(len=nctk_slen) :: "eph_fsewin", "eph_fsmear", "eph_extrael", "eph_fermie", "ph_wstep", "ph_smear"], &
     [dtset%eph_fsewin, dtset%eph_fsmear, dtset%eph_extrael, dtset%eph_fermie, dtset%ph_wstep, dtset%ph_smear])
   NCF_CHECK(ncerr)

   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'qibz'), gams%qibz))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'wtq'), gams%wtq))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ngqpt"), gamma_ngqpt))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "eph_ngqpt_fine"), dtset%eph_ngqpt_fine))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ph_ngqpt"), dtset%ph_ngqpt))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ddb_ngqpt"), dtset%ddb_ngqpt))
   !NCF_CHECK(nf90_close(ncid))
 end if

 !call xmpi_barrier(comm)

 ! Now reopen the file inside ncwrite_comm to perform pararallel-IO (required for q-point parallelism).
 !if (self%ncwrite_comm%value /= xmpi_comm_null) then
 !  NCF_CHECK(nctk_open_modify(self%ncid, path, self%ncwrite_comm%value))
 !  NCF_CHECK(nctk_set_datamode(self%ncid))
 !end if
#endif
 call edos%free()

 ! Open the DVDB file
 call dvdb%open_read(ngfftf, xmpi_comm_self)

 if (pert_comm%nproc > 1) then
   ! Build table with list of perturbations treated by this CPU inside pert_comm
   call ephtk_set_pertables(cryst%natom, my_npert, pert_table, my_pinfo, pert_comm%value)
   ! Activate parallelism over perturbations
   call dvdb%set_pert_distrib(my_npert, natom3, my_pinfo, pert_table, pert_comm%value)
   ABI_FREE(my_pinfo)
   ABI_FREE(pert_table)
 end if

 ! Check whether all q-points are available in the DVDB file
 do_ftv1q = 0
 do iq_ibz=1,gams%nqibz
   qpt = gams%qibz(:,iq_ibz)
   if (dvdb%findq(qpt) == -1) do_ftv1q = do_ftv1q + 1
   do spin=1,nsppol
     fs => fstab(spin)
     kqcount = 0
     do ik_bz=1,fs%nkfs
       kk = fs%kpts(:, ik_bz); kq = kk + qpt
       if (fs%findkg0(kq, g0bz_kq) == -1) cycle
       kqcount = kqcount + 1
     end do
     write(std_out,"((a,i0,2a,a,i0))")" For spin: ",spin,", qpt: ",trim(ktoa(qpt)),", number of (k,q) pairs: ",kqcount
   end do
 end do
 call wrtout(std_out, " ", do_flush=.True.)

 if (do_ftv1q /= 0) then
   call wrtout([std_out, ab_out], " Cannot find eph_ngqpt_fine q-points in DVDB --> Activating Fourier interpolation.")
   ! Prepare Fourier interpolation of DFPT potentials.
   comm_rpt = xmpi_comm_self
   !comm_rpt = bqs_comm%value
   method = dtset%userid
   call dvdb%ftinterp_setup(dtset%dvdb_ngqpt, dtset%ddb_qrefine, 1, dtset%ddb_shiftq, nfftf, ngfftf, method, comm_rpt)
 end if

 ! Initialize the wave function descriptor.
 ABI_MALLOC(nband, (nkibz, nsppol))
 ABI_MALLOC(bks_mask, (mband, nkibz, nsppol))
 ABI_MALLOC(keep_ur, (mband, nkibz ,nsppol))
 nband = mband; bks_mask = .False.; keep_ur = .False.

 ! Only wavefunctions on the FS are stored in wfd.
 ! Need all k-points on the FS because of k+q, spin is not distributed for the time being.
 ! One could reduce the memory allocated per MPI-rank via MPI-FFT or OpenMP.
 do mys=1,gams%my_nspins
   spin = gams%my_spins(mys)
   fs => fstab(spin)
   do ik_bz=1,fs%nkfs
     ik_ibz = fs%indkk_fs(1, ik_bz)
     bstart_k = fs%bstart_cnt_ibz(1, ik_ibz); nband_k = fs%bstart_cnt_ibz(2, ik_ibz)
     bks_mask(bstart_k:bstart_k+nband_k-1, ik_ibz, spin) = .True.
   end do
 end do
 !bks_mask(1:mband,:,:) = .True. ! no memory distribution, each node has the full set of states.

 ! Impose istwfk = 1 for all k points. This is also done in respfn (see inkpts)
 ! wfd_read_wfk will handle a possible conversion if WFK contains istwfk /= 1.
 ABI_MALLOC(wfd_istwfk, (nkibz))
 wfd_istwfk = 1

 ecut = dtset%ecut ! dtset%dilatmx
 call wfd_init(wfd, cryst, pawtab, psps, keep_ur, mband, nband, nkibz, nsppol, bks_mask,&
   nspden, nspinor, ecut, dtset%ecutsm, dtset%dilatmx, wfd_istwfk, ebands%kptns, ngfft,&
   dtset%nloalg, dtset%prtvol, dtset%pawprtvol, comm)

 call wfd%print(header="Wavefunctions on the Fermi surface")

 ABI_FREE(nband)
 ABI_FREE(bks_mask)
 ABI_FREE(keep_ur)
 ABI_FREE(wfd_istwfk)

 call wfd%read_wfk(wfk0_path, iomode_from_fname(wfk0_path))

 ! one-dimensional structure factor information on the coarse grid.
 ABI_MALLOC(ph1d, (2, 3*(2*mgfft+1)*natom))
 call getph(cryst%atindx, natom, n1, n2, n3, ph1d, cryst%xred)

 ! mpw is the maximum number of plane-waves over k and k+q where k and k+q are in the BZ.
 ! we also need the max components of the G-spheres (k, k+q) in order to allocate the workspace array work
 ! that will be used to symmetrize the wavefunctions in G-space.
 call cwtime(cpu, wall, gflops, "start")
 mpw = 0; gmax = 0; cnt = 0
 do iq_ibz=1,gams%nqibz
   qpt = gams%qibz(:,iq_ibz)
   do spin=1,nsppol
     fs => fstab(spin)
     do ik_bz=1,fs%nkfs
       cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle ! MPI parallelism inside comm
       kk = fs%kpts(:, ik_bz); kq = kk + qpt
       ! Compute G sphere, returning npw. Note istwfk == 1.
       call get_kg(kk, 1, ecut, cryst%gmet, onpw, gtmp)
       mpw = max(mpw, onpw)
       do ipw=1,onpw
         do ii=1,3
          gmax(ii) = max(gmax(ii), abs(gtmp(ii,ipw)))
         end do
       end do
       ABI_FREE(gtmp)

       ! TODO: g0 umklapp here can enter into play!
       ! fstab should contains the max of the umlapp G-vectors.
       ! gmax could not be large enough!
       call get_kg(kq, 1, ecut, cryst%gmet, onpw, gtmp)
       mpw = max(mpw, onpw)
       do ipw=1,onpw
         do ii=1,3
          gmax(ii) = max(gmax(ii), abs(gtmp(ii,ipw)))
         end do
       end do
       ABI_FREE(gtmp)
     end do
   end do
 end do
 call cwtime_report(" gmax and mpw", cpu, wall, gflops)

 my_mpw = mpw; call xmpi_max(my_mpw, mpw, comm, ierr)
 my_gmax = gmax; call xmpi_max(my_gmax, gmax, comm, ierr)
 call wrtout(std_out, sjoin(' Optimal value of mpw: ', itoa(mpw)))

 ! Init work_ngfft
 gmax = gmax + 4 ! FIXME: this is to account for umklapp
 gmax = 2*gmax + 1
 call ngfft_seq(work_ngfft, gmax)
 !write(std_out,*)"work_ngfft(1:3): ",work_ngfft(1:3)
 ABI_MALLOC(work, (2, work_ngfft(4), work_ngfft(5), work_ngfft(6)))

 ! Allow PW-arrays dimensioned with mpw
 ABI_MALLOC(kg_k, (3, mpw))
 ABI_MALLOC(kg_kq, (3, mpw))

 ! Spherical Harmonics for useylm == 1.
 ABI_MALLOC(ylm_k, (mpw, psps%mpsang*psps%mpsang*psps%useylm))
 ABI_MALLOC(ylm_kq, (mpw, psps%mpsang*psps%mpsang*psps%useylm))
 ABI_MALLOC(ylmgr_kq, (mpw, 3, psps%mpsang*psps%mpsang*psps%useylm*useylmgr1))

 ! TODO FOR PAW
 usecprj = 0
 ABI_MALLOC(cwaveprj0, (natom, nspinor*usecprj))

 ! Prepare call to getgh1c
 usevnl = 0
 optlocal = 1    ! local part of H^(1) is computed in gh1c=<G|H^(1)|C>
 optnl = 2       ! non-local part of H^(1) is totally computed in gh1c=<G|H^(1)|C>
 opt_gvnlx1 = 0  ! gvnlx1 is output
 ABI_MALLOC(gvnlx1, (2, usevnl))
 ABI_MALLOC(grad_berry, (2, nspinor*(berryopt0/4)))

 ! This part is taken from dfpt_vtorho
 !==== Initialize most of the Hamiltonian (and derivative) ====
 !1) Allocate all arrays and initialize quantities that do not depend on k and spin.
 !2) Perform the setup needed for the non-local factors:
 ! Norm-conserving: Constant kleimann-Bylander energies are copied from psps to gs_hamk.
 ! PAW: Initialize the overlap coefficients and allocate the Dij coefficients.

 call init_hamiltonian(gs_hamkq, psps, pawtab, nspinor, nsppol, nspden, natom, &
   dtset%typat, cryst%xred, nfft, mgfft, ngfft, cryst%rprimd, dtset%nloalg, &
   comm_atom=mpi_enreg%comm_atom, mpi_atmtab=mpi_enreg%my_atmtab, mpi_spintab=mpi_enreg%my_isppoltab, &
   usecprj=usecprj, ph1d=ph1d, nucdipmom=dtset%nucdipmom, use_gpu_cuda=dtset%use_gpu_cuda)

 ! Allocate vlocal. Note nvloc
 ! I set vlocal to huge to trigger possible bugs (DFPT routines should not access the data)
 ABI_MALLOC(vlocal, (n4, n5, n6, gs_hamkq%nvloc))
 vlocal = huge(one)

 ! Allocate work space arrays.
 ABI_MALLOC(tgam, (2, natom3, natom3))
 ABI_MALLOC(displ_cart, (2, 3, cryst%natom, natom3))
 ABI_MALLOC(displ_red, (2, 3, cryst%natom, natom3))

 ABI_CALLOC(dummy_vtrial, (nfftf, nspden))
 ! TODO: Save data to netcdf file for each q in IBZ and then read data to build a2Fw once all big
 ! datastructures (wfd, dvdb) have been deallocated.
 ! As a side effect, one can also implement restart over q-points

 ddkop = ddkop_new(dtset, cryst, pawtab, psps, wfd%mpi_enreg, mpw, wfd%ngfft)

 if (dtset%eph_transport > 0) then
   ABI_MALLOC(tgamvv_in, (2, 9, natom3, natom3))
   ABI_MALLOC(tgamvv_out, (2, 9, natom3, natom3))
 end if

 ! Loop over my q-point in the IBZ.
 do imyq=1,gams%my_nqibz
   iq_ibz = gams%my_iqibz(imyq)
   call cwtime(cpu_q, wall_q, gflops_q, "start")

   qpt = gams%qibz(:, iq_ibz)
   msg = sjoin("[", itoa(iq_ibz), "/", itoa(gams%nqibz), "]")
   call wrtout(std_out, sjoin(" Computing phonon linewidths for IBZ q-point:", ktoa(qpt), msg))

   tgam = zero
   if (dtset%eph_transport > 0) then
     tgamvv_in = zero; tgamvv_out = zero
   end if

   if (do_ftv1q == 0) then
     ! No interpolation --> Find the index of the q-point in the DVDB.
     db_iqpt = dvdb%findq(qpt)

     if (db_iqpt /= -1) then
       ! Read and reconstruct the dvscf potentials for all 3*natom perturbations.
       ! This call allocates v1scf(cplex, nfftf, nspden, 3*natom))
       call dvdb%readsym_allv1(db_iqpt, cplex, nfftf, ngfftf, v1scf, pkb_comm%value)

       if (dvdb%my_npert /= natom3) then
         ! Extract my npert from v1scf
         ABI_MALLOC(v1_work, (cplex, nfftf, nspden, dvdb%my_npert))
         do imyp=1,dvdb%my_npert
           v1_work(:,:,:,imyp) = v1scf(:,:,:,dvdb%my_pinfo(3, imyp))
         end do
         ABI_MOVE_ALLOC(v1_work, v1scf)
       end if

     else
       MSG_ERROR(sjoin("Cannot find q-point:", ktoa(qpt), "in the DVDB file."))
     end if

   else
     ! Use Fourier interpolation of DFPT potentials to get my_npert potentials.
     cplex = 2
     ABI_MALLOC(v1scf, (cplex, nfft, nspden, dvdb%my_npert))
     call dvdb%ftinterp_qpt(qpt, nfftf, ngfftf, v1scf, dvdb%comm_rpt)
   end if

   ! Examine the symmetries of the q wavevector.
   call littlegroup_q(cryst%nsym, qpt, symq, cryst%symrec, cryst%symafm, timerev_q, prtvol=dtset%prtvol)

   ! Get phonon frequencies and eigenvectors.
   call ifc%fourq(cryst, qpt, phfrq, displ_cart, out_displ_red=displ_red)

   ! Allocate vlocal1 with correct cplex. Note nvloc and my_npert.
   ABI_MALLOC(vlocal1, (cplex*n4, n5, n6, gs_hamkq%nvloc, my_npert))

   ! Loop over my spins.
   do mys=1,gams%my_nspins
     spin = gams%my_spins(mys)
     fs => fstab(spin)

     if (dtset%prteliash == 3) then
       ABI_MALLOC_OR_DIE(tmp_vals_ee, (2, gams%nene, gams%nene, gams%natom3, gams%natom3), ierr)
       tmp_vals_ee = zero
       ! Energy mesh for electrons in a2F(e,e',w) NB: It depends on the spin.
       ABI_MALLOC(emesh, (gams%nene))
       emesh = arth(fs%enemin, fs%deltaene, gams%nene)
     end if

     ! Set up local potential vlocal1 with proper dimensioning from vtrial1 taking into account the spin.
     do imyp=1,my_npert
       call rf_transgrid_and_pack(spin, nspden, psps%usepaw, cplex, nfftf, nfft, ngfft, gs_hamkq%nvloc,&
                 pawfgr, mpi_enreg, dummy_vtrial, v1scf(:,:,:,imyp), vlocal, vlocal1(:,:,:,:,imyp))
     end do

     ! Continue to initialize the GS Hamiltonian
     call gs_hamkq%load_spin(spin, vlocal=vlocal, with_nonlocal=.true.)

     ! Allocate workspace for wavefunctions. Make npw larger than expected.
     ! maxnb is the maximum number of bands crossing the FS, used to dimension arrays.
     mnb = fs%maxnb
     ABI_MALLOC(bras_kq, (2, mpw*nspinor, mnb))
     ABI_MALLOC(kets_k, (2, mpw*nspinor, mnb))
     ABI_MALLOC(h1kets_kq, (2, mpw*nspinor, mnb))
     ABI_MALLOC(gkk_atm, (2, mnb, mnb, natom3))

     ! The weights for FS integration.
     ABI_MALLOC(dbldelta_wts, (mnb, mnb))

     if (dtset%eph_transport > 0) then
       ABI_MALLOC(vv_kk, (9, mnb, mnb))
       ABI_MALLOC(vv_kkq, (9, mnb, mnb))
     end if

     if (dtset%prteliash == 3) then
       ABI_MALLOC(wt_ek, (gams%nene, mnb))
       ABI_MALLOC(wt_ekq, (gams%nene, mnb))
     end if

     ! =====================================
     ! Integration over the FS for this spin
     ! =====================================
     ! Compute integration weights and distribute k-points (my_nfsk_q)
     call phgamma_setup_qpoint(gams, fs, cryst, ebands, spin, ltetra, qpt, nesting, kpt_comm%value)

     do myik=1,gams%my_nfsk_q
       call cwtime(cpu_k, wall_k, gflops_k, "start")

       ! The k-point and the symmetries relating the BZ k-point to the IBZ.
       ik_fs = gams%my_ifsk_q(myik)
       kk = fs%kpts(:, ik_fs)
       ik_ibz = fs%indkk_fs(1, ik_fs); isym_k = fs%indkk_fs(2, ik_fs)
       trev_k = fs%indkk_fs(6, ik_fs); g0_k = fs%indkk_fs(3:5,ik_fs)
       kk_ibz = ebands%kptns(:,ik_ibz)

       ! Number of bands crossing the Fermi level at k
       bstart_k = fs%bstart_cnt_ibz(1, ik_ibz); nband_k = fs%bstart_cnt_ibz(2, ik_ibz)

       ! Find k+q in the extended zone and extract symmetry info. cycle if k+q not in FS.
       ! Be careful here because there are two umklapp vectors to be considered:
       !
       !   k + q = k_bz + g0_bz = IS(k_ibz) + g0_ibz + g0_bz
       !
       kq = kk + qpt; ikq_bz = fs%findkg0(kq, g0bz_kq)
       ! Skip this point if kq does not belong to the FS window.
       if (ikq_bz == -1) cycle

       call listkk(dksqmax, cryst%gmet, indkk_kq, ebands%kptns, kq, ebands%nkpt, 1, cryst%nsym,&
          1, cryst%symafm, cryst%symrel, timrev1, xmpi_comm_self, use_symrec=.False.)

       if (dksqmax > tol12) then
         write(msg, '(3a,es16.6,7a)' ) &
          "The WFK file cannot be used to compute phonon linewidths.",ch10, &
          "At least one of the k-points on the FS could not be generated from a symmetrical one. dksqmax: ", dksqmax,ch10, &
          "Q-mesh: ", ltoa(gamma_ngqpt), ", K-mesh (from kptrlatt) ", ltoa(get_diag(dtset%kptrlatt)), &
          'Action: check your WFK file and the (k, q) point input variables.'
          MSG_ERROR(msg)
       end if

       ikq_ibz = indkk_kq(1, 1); isym_kq = indkk_kq(1, 2)
       trev_kq = indkk_kq(1, 6); g0_kq = indkk_kq(1, 3:5)
       kq_ibz = ebands%kptns(:, ikq_ibz)

       ! Number of bands crossing the Fermi level at k+q
       bstart_kq = fs%bstart_cnt_ibz(1, ikq_ibz); nband_kq = fs%bstart_cnt_ibz(2, ikq_ibz)
       ABI_CHECK(nband_k <= mnb .and. nband_kq <= mnb, "wrong nband")

       ! Get npw_k, kg_k and symmetrize wavefunctions from IBZ (if needed).
       call wfd%sym_ug_kg(ecut, kk, kk_ibz, bstart_k, nband_k, spin, mpw, fs%indkk_fs(:,ik_fs), cryst, &
                          work_ngfft, work, istwf_k, npw_k, kg_k, kets_k)

       ! Get npw_kq, kg_kq and symmetrize wavefunctions from IBZ (if needed).
       call wfd%sym_ug_kg(ecut, kq, kq_ibz, bstart_kq, nband_kq, spin, mpw, indkk_kq(1,:), cryst, &
                          work_ngfft, work, istwf_kq, npw_kq, kg_kq, bras_kq)

       ! if PAW, one has to solve a generalized eigenproblem
       ! Be careful here because I will need sij_opt==-1
       gen_eigenpb = psps%usepaw == 1; sij_opt = 0; if (gen_eigenpb) sij_opt = 1
       ABI_MALLOC(gs1c, (2, npw_kq*nspinor*((sij_opt+1)/2)))

       ! Set up the spherical harmonics (Ylm) at k and k+q. See also dfpt_looppert
       !if (psps%useylm == 1) then
       !   optder = 0; if (useylmgr == 1) optder = 1
       !   call initylmg(cryst%gprimd, kg_k, kk, mkmem1, mpi_enreg, psps%mpsang, mpw, nband, mkmem1,&
       !     [npw_k], dtset%nsppol, optder, cryst%rprimd, ylm_k, ylmgr)
       !   call initylmg(cryst%gprimd, kg_kq, kq, mkmem1, mpi_enreg, psps%mpsang, mpw, nband, mkmem1,&
       !     [npw_kq], dtset%nsppol, optder, cryst%rprimd, ylm_kq, ylmgr_kq)
       !end if

       ! Loop over all my perturbations.
       gkk_atm = zero
       do imyp=1,my_npert
         idir = dvdb%my_pinfo(1, imyp); ipert = dvdb%my_pinfo(2, imyp); ipc = dvdb%my_pinfo(3, imyp)

         ! Prepare application of the NL part.
         call init_rf_hamiltonian(cplex, gs_hamkq, ipert, rf_hamkq, has_e1kbsc=.true.)

         call rf_hamkq%load_spin(spin, vlocal1=vlocal1(:,:,:,:,imyp), with_nonlocal=.true.)

         ! This call is not optimal because there are quantities in out that do not depend on idir,ipert
         call getgh1c_setup(gs_hamkq, rf_hamkq, dtset, psps, kk, kq, idir, ipert, &  ! In
           cryst%natom, cryst%rmet, cryst%gprimd, cryst%gmet, istwf_k, &             ! In
           npw_k, npw_kq, useylmgr1, kg_k, ylm_k, kg_kq, ylm_kq, ylmgr_kq, &         ! In
           dkinpw, nkpg, nkpg1, kpg_k, kpg1_k, kinpw1, ffnlk, ffnl1, ph3d, ph3d1)    ! Out

         ! Calculate dvscf * psi_k, results stored in h1kets_kq on the k+q sphere.
         ! Compute H(1) applied to GS wavefunction Psi(0)
         do ib2=1,nband_k
           band_k = ib2 + bstart_k - 1
           eig0nk = ebands%eig(band_k, ik_ibz, spin)
           ! Use scissor shift on 0-order eigenvalue
           eshift = eig0nk - dtset%dfpt_sciss

           call getgh1c(berryopt0, kets_k(:,:,ib2), cwaveprj0, h1kets_kq(:,:,ib2), &
                        grad_berry, gs1c, gs_hamkq, gvnlx1, idir, ipert, eshift, mpi_enreg, optlocal, &
                        optnl, opt_gvnlx1, rf_hamkq, sij_opt, tim_getgh1c, usevnl)
         end do

         call rf_hamkq%free()

         ABI_FREE(kinpw1)
         ABI_FREE(kpg1_k)
         ABI_FREE(kpg_k)
         ABI_FREE(dkinpw)
         ABI_FREE(ffnlk)
         ABI_FREE(ffnl1)
         ABI_FREE(ph3d)
         ABI_SFREE(ph3d1)

         ! Calculate elphmat(j,i) = <psi_{k+q,j}|dvscf_q*psi_{k,i}> for this perturbation.
         ! No need to handle istwf_kq because it's always 1.
         ! The array eig1_k contains:
         !
         ! <u_(band,k+q)^(0)|H_(k+q,k)^(1)|u_(band,k)^(0)>                           (NC psps)
         ! <u_(band,k+q)^(0)|H_(k+q,k)^(1)-(eig0_k+eig0_k+q)/2.S^(1)|u_(band,k)^(0)> (PAW)
         do ib2=1,nband_k
           do ib1=1,nband_kq
             gkk_atm(:, ib1, ib2, ipc) = cg_zdotc(npw_kq*nspinor, bras_kq(1,1,ib1), h1kets_kq(1,1,ib2))
           end do
         end do

       end do ! imyp (loop over my_npert atomic perturbations)

       ABI_FREE(gs1c)

       ! Collect gkk_atm inside pert_comm so that all procs can operate on data.
       if (pert_comm%nproc > 1) call xmpi_sum(gkk_atm, pert_comm%value, ierr)

       !call ephtk_gkknu_from_atm(1, nbcalc_ks, 1, natom, gkq_atm, phfrq, displ_red, gkq_nu)

       need_velocities = dtset%eph_transport > 0 .or. fs%eph_fsmear < zero .or. nesting /= 0

       if (need_velocities) then
         ! Compute diagonal matrix elements of velocity operator with DFPT routines
         ! Velocities are in Cartesian coordinates.
         call ddkop%setup_spin_kpoint(dtset, cryst, psps, spin, kk, istwf_k, npw_k, kg_k)
         do ib2=1,nband_k
           band_k = ib2 + bstart_k - 1
           fs%vk(:,ib2) = ddkop%get_vdiag(ebands%eig(band_k, ik_ibz, spin), &
                                          istwf_k, npw_k, wfd%nspinor, kets_k(:,:,ib2), cwaveprj0)
         end do

         call ddkop%setup_spin_kpoint(dtset, cryst, psps, spin, kq, istwf_kq, npw_kq, kg_kq)
         do ib1=1,nband_kq
           band_kq = ib1 + bstart_kq - 1
           fs%vkq(:,ib1) = ddkop%get_vdiag(ebands%eig(band_kq, ikq_ibz, spin), &
                                           istwf_kq, npw_kq, wfd%nspinor, bras_kq(:,:,ib1), cwaveprj0)
         end do
       end if

       ! Compute weights for double delta integration at the Fermi level.
       call fs%get_dbldelta_weights(ebands, ik_fs, ik_ibz, ikq_ibz, spin, nesting, dbldelta_wts)

       ! Accumulate results in tgam (sum over FS and bands for this spin).
       do ipc2=1,natom3
         do ipc1=1,natom3
           do ib2=1,nband_k
             do ib1=1,nband_kq
               lf = gkk_atm(:, ib1, ib2, ipc1)
               rg = gkk_atm(:, ib1, ib2, ipc2)
               res(1) = lf(1) * rg(1) + lf(2) * rg(2)
               res(2) = lf(1) * rg(2) - lf(2) * rg(1)
               tgam(:,ipc1,ipc2) = tgam(:,ipc1,ipc2) + res(:) * dbldelta_wts(ib1, ib2)
             end do
           end do
         end do
       end do

       if (dtset%eph_transport > 0) then
         ! TODO: could almost make this a BLAS call plus a reshape...
         do ib2 = 1,nband_k
           do ipc2 = 1,3
             do ib1 = 1,nband_kq
               do ipc1 = 1,3
                 ! vk x vk
                 vv_kk(ipc1 + (ipc2-1)*3, ib1,ib2) = fs%vk(ipc1, ib1) * fs%vk(ipc2, ib2)
                    !ddk%velocity(ipc1,ib1,ik_fs,spin) * ddk%velocity(ipc2,ib2,ik_fs,spin)
                 ! vk x vk+q
                 vv_kkq(ipc1 + (ipc2-1)*3, ib1,ib2) = fs%vkq(ipc1, ib1) * fs%vk(ipc2, ib2)
                    !ddk%velocity(ipc1,ib1,ikq_bz,spin) * ddk%velocity(ipc2,ib2,ik_fs,spin)
               end do
             end do
           end do
         end do

         ! Accumulate results (sum over FS and bands).
         do ipc2=1,natom3
           do ipc1=1,natom3
              do ib2=1,nband_k
                do ib1=1,nband_kq
                  lf = gkk_atm(:, ib1, ib2, ipc1)
                  rg = gkk_atm(:, ib1, ib2, ipc2)
                  ! res was missing in MJV version!
                  res(1) = lf(1) * rg(1) + lf(2) * rg(2)
                  res(2) = lf(1) * rg(2) - lf(2) * rg(1)
                  resvv_in(1,:) = res(1) * vv_kkq(:,ib1,ib2)
                  resvv_in(2,:) = res(2) * vv_kkq(:,ib1,ib2)
                  resvv_out(1,:) = res(1) * vv_kk(:,ib1,ib2)
                  resvv_out(2,:) = res(2) * vv_kk(:,ib1,ib2)
                  tgamvv_in(:,:,ipc1,ipc2)  = tgamvv_in(:,:,ipc1,ipc2) + resvv_in * dbldelta_wts(ib1, ib2)
                  tgamvv_out(:,:,ipc1,ipc2) = tgamvv_out(:,:,ipc1,ipc2) + resvv_out * dbldelta_wts(ib1, ib2)
                end do
              end do
           end do
         end do
       end if ! add transport things

       if (dtset%prteliash == 3) then

         ! Precompute deltas with gaussian (tetra is not supported here
         ! also because one should allocate (nene, nene, mband, mband)
         !call get_dbldelta_weights_emesh()
         do ib2=1,nband_k
           band_k = ib2 + bstart_k - 1
           sigma = fs%eph_fsmear
           if (fs%eph_fsmear < zero) then
             sigma = max(maxval([(abs(dot_product(fs%vk(:, ib2), fs%kmesh_cartvec(:,ii))), ii=1,3)]), fs%min_smear)
           end if
           wt_ek(:, ib2) = gaussian(emesh - ebands%eig(band_k, ik_ibz, spin), sigma) / sqrt(one * fs%nktot)
         end do

         do ib1=1,nband_kq
           band_kq = ib1 + bstart_kq - 1
           sigma = fs%eph_fsmear
           if (fs%eph_fsmear < zero) then
             sigma = max(maxval([(abs(dot_product(fs%vkq(:, ib1), fs%kmesh_cartvec(:,ii))), ii=1,3)]), fs%min_smear)
           end if
           wt_ekq(:, ib1) = gaussian(emesh - ebands%eig(band_kq, ikq_ibz, spin), sigma) / sqrt(one * fs%nktot)
         end do

         do ib2=1,nband_k
           do ib1=1,nband_kq
             do ipc2=1,natom3
               do ipc1=1,natom3
                 lf = gkk_atm(:, ib1, ib2, ipc1)
                 rg = gkk_atm(:, ib1, ib2, ipc2)
                 res(1) = lf(1) * rg(1) + lf(2) * rg(2)
                 res(2) = lf(1) * rg(2) - lf(2) * rg(1)
                 do jene = 1, gams%nene
                   do iene = 1, gams%nene
                     ! TODO: distribute this in procs over q. Make a temp array here for 1 q
                     ! then mpi sync it and save it only on 1 processor below after mpisum over k
                     tmp_vals_ee(:,iene,jene,ipc1,ipc2) = tmp_vals_ee(:,iene,jene,ipc1,ipc2) + &
                       res(:) * wt_ekq(iene, ib1) * wt_ek(jene, ib2)
                   end do
                 end do
               end do
             end do
           end do
         end do
       end if

       if (myik < 20 .or. (fs%nkfs > 100 .and. mod(myik, 200) == 0)) then
         write(msg,'(4(a,i0),a,f8.2)')" q-point [", iq_ibz, "/", gams%nqibz, "] k-point [", myik, "/", gams%my_nfsk_q, "]"
         call cwtime_report(msg, cpu_k, wall_k, gflops_k)
       end if

     end do ! ik_fs: sum over k-points on the BZ FS for this spin.

     ABI_FREE(dbldelta_wts)
     ABI_FREE(bras_kq)
     ABI_FREE(kets_k)
     ABI_FREE(h1kets_kq)
     ABI_FREE(gkk_atm)
     ABI_SFREE(wt_ek)
     ABI_SFREE(wt_ekq)

     ! Collect tgam values inside kpt comm and divide by the total number of k-points in the full mesh.
     call xmpi_sum(tgam, kpt_comm%value, ierr)
     ! Symmetrize gamma matrices. Note that this call is not executed in elphon!
     if (gams%symgamma == 1) call tgamma_symm(cryst, gams%qibz(:,iq_ibz), tgam)
     gams%vals_qibz(:,:,:,iq_ibz, spin) = tgam

     if (dtset%eph_transport > 0) then
       ABI_FREE(vv_kk)
       ABI_FREE(vv_kkq)
       call xmpi_sum(tgamvv_in, kpt_comm%value, ierr)
       call xmpi_sum(tgamvv_out, kpt_comm%value, ierr)
       gams%vals_in_qibz(:,:,:,:,iq_ibz,spin)  = tgamvv_in
       gams%vals_out_qibz(:,:,:,:,iq_ibz,spin) = tgamvv_out
     end if

     if (dtset%prteliash == 3) then
       call xmpi_sum(tmp_vals_ee, kpt_comm%value, ierr)
       !TODO: write arrays to ncfile.
       gams%vals_ee(:,:,:,:,:, iq_ibz, spin) = tmp_vals_ee
       ABI_FREE(tmp_vals_ee)
       ABI_FREE(emesh)
     end if

   end do ! spin

   ABI_FREE(v1scf)
   ABI_FREE(vlocal1)

   write(msg,'(2(a,i0),a)')" Computation of q-point [", iq_ibz, "/", gams%nqibz, "]"
   call cwtime_report(msg, cpu_q, wall_q, gflops_q, pre_str=ch10, end_str=ch10)
 end do ! iq_ibz

 call cwtime_report(" Full computation of phonon linewidths", cpu_all, wall_all, gflops_all, pre_str=ch10, end_str=ch10)

 ! Free memory
 ABI_FREE(gvnlx1)
 ABI_FREE(grad_berry)
 ABI_FREE(dummy_vtrial)
 ABI_FREE(work)
 ABI_FREE(ph1d)
 ABI_FREE(vlocal)
 ABI_FREE(kg_k)
 ABI_FREE(kg_kq)
 ABI_FREE(ylm_k)
 ABI_FREE(ylm_kq)
 ABI_FREE(ylmgr_kq)
 ABI_FREE(tgam)
 ABI_FREE(displ_cart)
 ABI_FREE(displ_red)

 if (dtset%eph_transport > 0) then
   ABI_FREE(tgamvv_in)
   ABI_FREE(tgamvv_out)
 end if

 call pawcprj_free(cwaveprj0)
 ABI_FREE(cwaveprj0)
 call ddkop%free()
 call gs_hamkq%free()
 call wfd%free()
 do spin=1,ebands%nsppol
   call fstab(spin)%free()
 end do
 ABI_FREE(fstab)

 ! Collect results on each node
 call xmpi_sum(gams%vals_qibz, qs_comm%value, ierr)
 if (dtset%eph_transport > 0) then
   call xmpi_sum(gams%vals_in_qibz, qs_comm%value, ierr)
   call xmpi_sum(gams%vals_out_qibz, qs_comm%value, ierr)
 end if
 if (dtset%prteliash == 3) then
   call xmpi_sum(gams%vals_ee, qs_comm%value, ierr)
 end if

#ifdef HAVE_NETCDF
 ! Close the netcdf file then master reopens it
 !if (ncwrite_comm%value /= xmpi_comm_null) then
 !  NCF_CHECK(nf90_close(ncid))
 !end if
 !call xmpi_barrier(comm)

 !ncid = nctk_noid
 !if (my_rank == master) then
 !  NCF_CHECK(nctk_open_modify(ncid, path, xmpi_comm_self))
 !  NCF_CHECK(nctk_set_datamode(ncid))
 !end if
#endif

 ! Deallocate MPI communicators.
 call pert_comm%free()
 call qpt_comm%free()
 call bsum_comm%free()
 call qs_comm%free()
 call kpt_comm%free()
 call spin_comm%free()
 call pkb_comm%free()
 !call ncwrite_comm%free()

 ! Print gamma(IBZ) to ab_out and ncid
 if (my_rank == master) call phgamma_ncwrite(gams, cryst, ifc, ncid)

 ! Interpolate linewidths along the q-path.
 if (dtset%ph_nqpath <= 0) then
   write(msg, '(7a,es16.6,4a)' ) &
    'You have not specified a path for the linewidth calculation - no interpolation or output will be done ',ch10,&
    'Action: check your input variables ph_nqpath and ph_qpath'
   MSG_WARNING(msg)
 else
   call gams%linwid(cryst, ifc, dtset%ph_ndivsm, dtset%ph_nqpath, dtset%ph_qpath, dtfil%filnam_ds(4), ncid, wminmax, comm)
 end if

 ! Compute a2Fw using ab-initio q-points (no interpolation here)
 call a2fw_init(a2fw, gams, cryst, ifc, dtset%ph_intmeth, dtset%ph_wstep, wminmax, dtset%ph_smear, &
   dtset%ph_ngqpt, dtset%ph_nqshift, dtset%ph_qshift, comm, qintp=.False., qptopt=1)

 if (my_rank == master) then
   call a2fw%write(strcat(dtfil%filnam_ds(4), "_NOINTP"), "_qcoarse", ncid)
   if (dtset%prteliash == 3) call a2fw_ee_write(a2fw, strcat(dtfil%filnam_ds(4), "_NOINTP"))
 end if
 call a2fw%free()

 ! Compute a2Fw using Fourier interpolation (R -> q) and ph_ngqpt grid.
 call a2fw_init(a2fw, gams, cryst, ifc, dtset%ph_intmeth, dtset%ph_wstep, wminmax, dtset%ph_smear, &
   dtset%ph_ngqpt, dtset%ph_nqshift, dtset%ph_qshift, comm, qptopt=1)

 if (my_rank == master) then
   call a2fw%write(dtfil%filnam_ds(4), "_qintp", ncid)
   if (dtset%prteliash == 3) call a2fw_ee_write(a2fw, dtfil%filnam_ds(4))
 end if

 ! TODO: Use KT mesh instead of T but read T from input.
 !ntemp = 6
 !temp_range = [0.6_dp, 1.2_dp]
 !wcut = 10 * wminmax(2); reltol = 0.001
 !call a2fw_solve_gap(a2fw,cryst,dtset%tmesh,wcut,dtset%eph_mustar,dtset%nstep,reltol,dtfil%filnam_ds(4),comm)

 call a2fw%free()

 ! Compute A2fw using Fourier interpolation and full BZ for debugging purposes.
 !call a2fw_init(a2fw, gams, cryst, ifc, dtset%ph_intmeth, dtset%ph_wstep, wminmax, dtset%ph_smear,&
 !  dtset%ph_ngqpt, dtset%ph_nqshift, dtset%ph_qshift, comm, qptopt=3)
 !if (my_rank == master) call a2fw%write(strcat(dtfil%filnam_ds(4), "_A2FW_QPTOPT3"), "fake", nctk_noid)
 !call a2fw%free()

 if (dtset%eph_transport == 1) then
   ! Calculate and output transport quantities

   ! Compute a2Fw_tr using ab-initio q-points (no interpolation)
   call a2fw_tr_init(a2fw_tr, gams, cryst, ifc, dtset%ph_intmeth, dtset%ph_wstep, wminmax, dtset%ph_smear,&
     dtset%ph_ngqpt, dtset%ph_nqshift, dtset%ph_qshift, comm, qintp=.False., qptopt=1)

   if (my_rank == master) call a2fw_tr%write(strcat(dtfil%filnam_ds(4), "_NOINTP"), "_qcoarse", ncid)
   call a2fw_tr%free()

   ! Compute a2Fw_tr using Fourier interpolation (R --> q) and ph_ngqpt grid
   call a2fw_tr_init(a2fw_tr, gams, cryst, ifc, dtset%ph_intmeth, dtset%ph_wstep, wminmax, dtset%ph_smear,&
     dtset%ph_ngqpt, dtset%ph_nqshift, dtset%ph_qshift, comm, qptopt=1)

   if (my_rank == master) call a2fw_tr%write(dtfil%filnam_ds(4), "_qintp", ncid)
   call a2fw_tr%free()

   ! Compute A2fw_tr using Fourier interpolation and full BZ for debugging purposes.
   !call a2fw_tr_init(a2fw_tr, gams, cryst, ifc, dtset%ph_intmeth, dtset%ph_wstep, wminmax, dtset%ph_smear,&
   !  dtset%ph_ngqpt, dtset%ph_nqshift, dtset%ph_qshift, comm, qptopt=3)
   !if (my_rank == master) call a2fw_tr%write(strcat(dtfil%filnam_ds(4), "_A2FWTR_QPTOPT3"), "fake", nctk_noid)
   !call a2fw_tr%free()
 end if

#ifdef HAVE_NETCDF
 if (my_rank == master) then
   NCF_CHECK(nf90_close(ncid))
 end if
#endif

 call gams%free()

end subroutine eph_phgamma
!!***

!----------------------------------------------------------------------

!!****f* m_phgamma/phgamma_setup_qpoint
!! NAME
!! phgamma_setup_qpoint
!!
!! FUNCTION
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine phgamma_setup_qpoint(gams, fs, cryst, ebands, spin, ltetra, qpt, nesting, comm)

!Arguments ------------------------------------
 type(phgamma_t),intent(inout) :: gams
 type(fstab_t),intent(inout) :: fs
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 integer,intent(in) :: spin, ltetra, comm
 integer,intent(out) :: nesting
!arrays
 real(dp),intent(in) :: qpt(3)

!Local variables-------------------------------
!scalars
 integer,parameter :: enough = 50
 integer :: nkbz, ierr, nb, ik_bz, ik_ibz, ikq_ibz, ikq_fs, ik_fs, i1, i2, i3, nkfs_q, nene
 integer :: ib1, ib2 ! band_k, band_kq
 real(dp),parameter :: max_occ1 = one
 real(dp) :: cpu, wall, gflops, enemin, enemax
 !character(len=500) :: msg
 type(krank_t) :: ibz_krank
 type(t_tetrahedron) :: tetra
 character(len=80) :: errorstring
!arrays
 integer :: nge(3), ngw(3), g0bz_kq(3)
 integer,allocatable :: select_ikfs(:), indkpt(:), kbz2fs(:) !, symrecfm(:,:,:)
 real(dp) :: kk(3), kq(3)
 real(dp),allocatable :: eig_k(:,:), eig_kq(:,:), wght_bz(:,:,:), kbz(:,:)
 real(dp),allocatable :: work_k(:), work_kq(:), dtweightde(:,:,:), tweight(:,:,:)

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 ! The double delta is ill-defined for q == 0. Set nesting to 1 and return
 nesting = 0
 if (fs%eph_intmeth == 2 .and. all(abs(qpt) < tol12)) then
   MSG_COMMENT("Tetrahedron for double grid with q = 0 is ill-defined. Returning nesting 1 to activate adaptive gaussian.")
   nesting = 1
 end if

 if (fs%eph_intmeth == 1 .or. nesting == 1) then
   ! Distribute k-points within the FS window inside comm.
   ! 1) Select k-points such that k+q is stil inside the FS window
   ! 2) Distribute effective k-points assuming all procs in comm have all FS k-points (no filtering)
   ABI_MALLOC(select_ikfs, (fs%nkfs))
   nkfs_q = 0
   do ik_fs=1,fs%nkfs
     kq = fs%kpts(:, ik_fs) + qpt
     ikq_fs = fs%findkg0(kq, g0bz_kq); if (ikq_fs == -1) cycle
     nkfs_q = nkfs_q + 1
     select_ikfs(nkfs_q) = ik_fs
   end do
   ABI_SFREE(gams%my_ifsk_q)
   call xmpi_split_list(nkfs_q, select_ikfs, comm, gams%my_nfsk_q, gams%my_ifsk_q)
   ABI_FREE(select_ikfs)
   call cwtime_report(" phgamma_setup_qpoint", cpu, wall, gflops)
   return
 end if

 ! Tetrahedron method:
 ! 1) Compute weights for double delta integration.
 ! 2) Filter k-points according to the weights (k must be in FS window with non-zero weight)
 ! 2) Distribute effective k-points assuming all procs have all FS k-points.

 fs%dbldelta_tetra_weights_kfs = zero

 ABI_CHECK(isdiagmat(ebands%kptrlatt), "kptrlatt must be diagonal when tetra is used.")
 ABI_CHECK(ebands%nshiftk == 1, "nshiftk must be 1 when tetra is used")
 nge = get_diag(ebands%kptrlatt); ngw = nge
 nkbz = product(nge(1:3))

 ! TODO: Handle symmetries in a cleaner way. Change API of krank_new to pass symafm and kptopt
 !timrev = 0; if (use_tr) timrev=1
 ibz_krank = krank_new(ebands%nkpt, ebands%kptns, nsym=cryst%nsym, symrec=cryst%symrec, time_reversal=.True.)

 nb = fs%maxnb
 ABI_MALLOC(kbz, (3, nkbz))
 ABI_MALLOC(indkpt, (nkbz))
 ABI_MALLOC(eig_k, (nb, nkbz))
 ABI_MALLOC(eig_kq, (nb, nkbz))
 ABI_MALLOC(kbz2fs, (nkbz))

 kbz2fs = -1; ierr = 0; ik_bz = 0
 do i3=0,nge(3) - 1
   do i2=0,nge(2) - 1
     do i1=0,nge(1) - 1
       ik_bz = ik_bz + 1
       indkpt(ik_bz) = ik_bz

       ! Find correspondence between the grid and the IBZ
       kk = ([i1, i2, i3] + ebands%shiftk(:, 1)) / nge(:)
       kbz(:, ik_bz) = kk

       ik_fs = fs%krank%get_index(kk)
       if (ik_fs /= -1) kbz2fs(ik_bz) = ik_fs

       ik_ibz = ibz_krank%get_index(kk)
       if (ik_ibz < 1) then
         if (ierr <= enough) then
           MSG_WARNING(sjoin('kpt:', trim(ktoa(kk)), 'has no symmetric among the k-points!'))
         end if
         ierr = ierr + 1; cycle
       end if

       eig_k(:, ik_bz) = ebands%eig(fs%bmin:fs%bmax, ik_ibz, spin)

       ! Find correspondence between the BZ grid and the IBZ.
       kq = kk + qpt
       ikq_ibz = ibz_krank%get_index(kq)

       if (ikq_ibz < 1) then
         if (ierr <= enough) then
           MSG_WARNING(sjoin('kpt + qpt:', trim(ktoa(kq)), 'has no symmetric among the k-points!'))
         end if
         ierr = ierr + 1; cycle
       end if

       eig_kq(:, ik_bz) = ebands%eig(fs%bmin:fs%bmax, ikq_ibz, spin)
     end do
   end do
 end do

 ABI_CHECK(ierr == 0, "See above warnings")
 call ibz_krank%free()

 if (any(ltetra == [1, 2])) then
   ! Compute weights for double delta integration. Note that libtetra assumes Ef set to zero.
   ! TODO: Average weights over degenerate states?
   write(std_out,"(a,i0,2a)")" Calling libtetrabz_dbldelta with ltetra: ", ltetra, " for q-point:", trim(ktoa(qpt))
   eig_k = eig_k - ebands%fermie; eig_kq = eig_kq - ebands%fermie
   ABI_MALLOC(wght_bz, (nb, nb, nkbz))
   call libtetrabz_dbldelta(ltetra, cryst%gprimd, nb, nge, eig_k, eig_kq, ngw, wght_bz) !, comm=comm)
   eig_k = eig_k + ebands%fermie; eig_kq = eig_kq + ebands%fermie

   ! Reindex from full BZ to fs% kpoints.
   do ik_bz=1,nkbz
     ik_fs = kbz2fs(ik_bz)
     if (ik_fs /= -1) then
       fs%dbldelta_tetra_weights_kfs(:,:,ik_fs) = wght_bz(:,:,ik_bz)
     else
       !write(std_out,*)"should be zero :", wght_bz(:,:,ik_bz)
     end if
   end do
   ABI_FREE(wght_bz)

 else if (ltetra == 3) then
   ! Tetrahedron method with Allen's approach for double delta.
   write(std_out,"(2a)")" Calling Allen's version for q-point: ", trim(ktoa(qpt))
   nene = 3
   enemin = ebands%fermie - tol6
   enemax = ebands%fermie + tol6

   ABI_MALLOC(dtweightde, (nkbz, nene, nene))
   ABI_MALLOC(tweight, (nkbz, nene, nene))
   ABI_MALLOC(work_k, (nkbz))
   ABI_MALLOC(work_kq, (nkbz))

   call init_tetra(indkpt, cryst%gprimd, fs%klatt, kbz, nkbz, tetra, ierr, errorstring, comm)
   ABI_CHECK(ierr == 0, "init_tetra returned ierr /= 0")

   fs%dbldelta_tetra_weights_kfs = zero
   do ib2=1,nb
     work_k = eig_k(ib2, :)
     do ib1=1,nb
      work_kq = eig_kq(ib1, :)
      ! TODO: Average weights over degenerate states?
      call get_dbl_tetra_weight(work_k, work_kq, enemin, enemax, enemin, enemax, &
                                max_occ1, nene, nene, nkbz, tetra, tweight, dtweightde, ierr)
      ABI_CHECK(ierr == 0, "get_dbldelta_weights returned ierr /= 0")

      ! Reindex from full BZ to fs% kpoints.
      do ik_bz=1,nkbz
        ik_fs = kbz2fs(ik_bz)
        if (ik_fs /= -1) then
          fs%dbldelta_tetra_weights_kfs(ib1,ib2,ik_fs) = dtweightde(ik_bz, 2, 2) ! / nkbz
        else
          !write(std_out,*)"should be zero :", wght_bz(:,:,ik_bz)
        end if
      end do
     end do
   end do

   call destroy_tetra(tetra)
   ABI_FREE(work_k)
   ABI_FREE(work_kq)
   ABI_FREE(dtweightde)
   ABI_FREE(tweight)
 else
   MSG_ERROR(sjoin("Invalid value of ltetra:", itoa(ltetra)))
 end if

 ! Now we can filter the k-points according to the tetra weights and distribute inside comm.
 ! Assuming all procs in comm have all k-points in the IBZ.
 ABI_MALLOC(select_ikfs, (fs%nkfs))
 nkfs_q = 0
 do ik_fs=1,fs%nkfs
   if (any(abs(fs%dbldelta_tetra_weights_kfs(:,:,ik_fs)) > zero)) then
     nkfs_q = nkfs_q + 1
     select_ikfs(nkfs_q) = ik_fs
   end if
 end do

 ABI_SFREE(gams%my_ifsk_q)
 call xmpi_split_list(nkfs_q, select_ikfs, comm, gams%my_nfsk_q, gams%my_ifsk_q)
 ABI_FREE(select_ikfs)

 write(std_out,"(2(a,i0),/)") &
   " Treating ", gams%my_nfsk_q, " k-points in the FS window over total nkfs: ", fs%nkfs

 ABI_FREE(kbz2fs)
 ABI_FREE(indkpt)
 ABI_FREE(kbz)
 ABI_FREE(eig_k)
 ABI_FREE(eig_kq)

 call cwtime_report(" phgamma_setup_qpoint", cpu, wall, gflops)

end subroutine phgamma_setup_qpoint
!!***

!----------------------------------------------------------------------

!!****f* m_phgamma/find_ewin
!! NAME
!! find_ewing
!!
!! FUNCTION
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine find_ewin(nqibz, qibz, cryst, ebands, ltetra, fs_ewin, comm)

!Arguments ------------------------------------
 integer,intent(in) :: nqibz
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 integer,intent(in) :: ltetra, comm
 real(dp),intent(out) :: fs_ewin
!arrays
 real(dp),intent(in) :: qibz(3, nqibz)

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: my_rank, iq_ibz, spin, ii
 real(dp) :: cpu, wall, gflops
 integer :: bstarts(3), bstops(3)
 real(dp) :: elows(3), ehighs(3), ewins(3), qsums(3)
 real(dp), allocatable :: wtqs(:,:,:)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm)

 call cwtime(cpu, wall, gflops, "start")

 ABI_MALLOC(wtqs, (nqibz, ebands%nsppol, 3))

 ! 1 is low, 2 is high, 3 is for the workspace
 ewins(1) = half * eV_Ha; elows(1) = ebands%fermie - ewins(1); ehighs(1) = ebands%fermie + ewins(1)
 call get_bands_from_erange(ebands, elows(1), ehighs(1), bstarts(1), bstops(1))
 call calc_dbldelta(cryst, ebands, ltetra, bstarts(1), bstops(1), nqibz, qibz, wtqs(:,:,1), comm)

 ewins(2) = five * eV_Ha; elows(2) = ebands%fermie - ewins(2); ehighs(2) = ebands%fermie + ewins(2)
 call get_bands_from_erange(ebands, elows(2), ehighs(2), bstarts(2), bstops(2))
 call calc_dbldelta(cryst, ebands, ltetra, bstarts(2), bstops(2), nqibz, qibz, wtqs(:,:,2), comm)

 if (abs(sum(abs(wtqs(:,:,1)) - sum(abs(wtqs(:,:,2))))) / sum(abs(wtqs(:,:,2))) < tol2) then
   fs_ewin = ewins(1)
   call print_weights(1)
   goto 100
 end if

 do
   ewins(3) = (ewins(1) + ewins(2)) / two; elows(3) = ebands%fermie - ewins(3); ehighs(3) = ebands%fermie + ewins(3)
   call get_bands_from_erange(ebands, elows(3), ehighs(3), bstarts(3), bstops(3))
   call calc_dbldelta(cryst, ebands, ltetra, bstarts(3), bstops(3), nqibz, qibz, wtqs(:,:,3), comm)

   do ii=1,3
     qsums(ii) = sum(abs(wtqs(:,:,ii))) / nqibz
   end do

   write(std_out,*)qsums(1), ", for ewins(1)_eV", ewins(1) * Ha_eV, " bstart(1)", bstarts(1), " bstopt(1)", bstops(1)
   write(std_out,*)qsums(2), ", for ewins(2)_eV", ewins(2) * Ha_eV, " bstart(1)", bstarts(2), " bstopt(1)", bstops(2)
   write(std_out,*)qsums(3), ", for ewins(3)_eV", ewins(3) * Ha_eV, " bstart(1)", bstarts(3), " bstopt(1)", bstops(3)

   !if (abs(sum(abs(wtqs(:,:,3)) - sum(abs(wtqs(:,:,2))))) / sum(abs(wtqs(:,:,2))) < tol3) then
   !  !  If 3 is close to 2, try to reduce the window by moving towards 1
   !  ewins(2) = ewins(3); elows(2) = elows(3); ehighs(2) = ehighs(3)

   if (abs(qsums(3) - qsums(2)) / qsums(2) < tol2) then
     fs_ewin = ewins(3)
     call print_weights(3)
     exit

   else
    ewins(1) = ewins(3); elows(1) = elows(3); ehighs(1) = ehighs(3)
   end if

 end do

100 continue
 ABI_FREE(wtqs)

 call cwtime_report(" find_ewin", cpu, wall, gflops)

contains

subroutine print_weights(ind)

 integer,intent(in) :: ind

 if (my_rank == master) then
   write(std_out, "(a, f8.3, a)")" Optimal FS energy window: ", fs_ewin * Ha_eV, " (eV)"
   do iq_ibz=1,nqibz
     write(std_out, "(2a)")" For q-point:", trim(ktoa(qibz(:, iq_ibz)))
     do spin=1,ebands%nsppol
       write(std_out, *) wtqs(iq_ibz, spin, ind)
     end do
   end do
 end if

end subroutine print_weights

end subroutine find_ewin
!!***

!----------------------------------------------------------------------

!!****f* m_phgamma/calc_dbldelta
!! NAME
!! calc_dbldelta
!!
!! FUNCTION
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine calc_dbldelta(cryst, ebands, ltetra, bstart, bstop, nqibz, qibz, wtqs, comm)

!Arguments ------------------------------------
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 integer,intent(in) :: ltetra, bstart, bstop, nqibz, comm
!arrays
 real(dp),intent(in) :: qibz(3,nqibz)
 real(dp),intent(out) :: wtqs(nqibz, ebands%nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: enough = 50
 integer :: nkbz, ierr, nb, ik_bz, ik_ibz, ikq_ibz, i1, i2, i3, nene, spin, iq_ibz
 integer :: ib1, ib2, my_rank, nproc, cnt
 real(dp),parameter :: max_occ1 = one
 real(dp) :: enemin, enemax !cpu, wall, gflops,
 !character(len=500) :: msg
 type(krank_t) :: ibz_krank
 type(t_tetrahedron) :: tetra
 !character(len=80) :: errorstring
!arrays
 integer :: nge(3), ngw(3)
 integer,allocatable :: indkpt(:) !, symrecfm(:,:,:)
 real(dp) :: qpt(3), kk(3), kq(3)
 real(dp),allocatable :: eig_k(:,:), eig_kq(:,:), wght_bz(:,:,:), kbz(:,:)
 real(dp),allocatable :: work_k(:), work_kq(:), dtweightde(:,:,:), tweight(:,:,:)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)
 !call cwtime(cpu, wall, gflops, "start")

 ABI_CHECK(nqibz > 1, "Need more that 1 q-point")
 ABI_CHECK(isdiagmat(ebands%kptrlatt), "kptrlatt must be diagonal when tetra is used.")
 ABI_CHECK(ebands%nshiftk == 1, "nshiftk must be 1 when tetra is used")
 nge = get_diag(ebands%kptrlatt); ngw = nge
 nkbz = product(nge(1:3))

 ! TODO: Handle symmetries in a cleaner way. Change API of krank_new to pass symafm and kptopt
 !timrev = 0; if (use_tr) timrev=1
 ibz_krank = krank_new(ebands%nkpt, ebands%kptns, nsym=cryst%nsym, symrec=cryst%symrec, time_reversal=.True.)

 nb = bstop - bstart + 1
 ABI_MALLOC(kbz, (3, nkbz))
 ABI_MALLOC(indkpt, (nkbz))
 ABI_MALLOC(eig_k, (nb, nkbz))
 ABI_MALLOC(eig_kq, (nb, nkbz))

 wtqs = zero

 cnt = 0
 do spin=1,ebands%nsppol
   do iq_ibz=1,nqibz
     cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle ! MPI parallelism.
     qpt = qibz(:, iq_ibz)

     ! The double delta is ill-defined for q == 0
     if (all(abs(qpt) < tol12)) cycle

     ierr = 0; ik_bz = 0
     do i3=0,nge(3) - 1
       do i2=0,nge(2) - 1
         do i1=0,nge(1) - 1
           ik_bz = ik_bz + 1
           indkpt(ik_bz) = ik_bz

           ! Find correspondence between the grid and the IBZ
           kk = ([i1, i2, i3] + ebands%shiftk(:, 1)) / nge(:)
           kbz(:, ik_bz) = kk

           ik_ibz = ibz_krank%get_index(kk)
           if (ik_ibz < 1) then
             if (ierr <= enough) then
               MSG_WARNING(sjoin('kpt:', trim(ktoa(kk)), 'has no symmetric among the k-points!'))
             end if
             ierr = ierr + 1; cycle
           end if

           eig_k(:, ik_bz) = ebands%eig(bstart:bstop, ik_ibz, spin)

           ! Find correspondence between the BZ grid and the IBZ.
           kq = kk + qpt
           ikq_ibz = ibz_krank%get_index(kq)

           if (ikq_ibz < 1) then
             if (ierr <= enough) then
               MSG_WARNING(sjoin('kpt + qpt:', trim(ktoa(kq)), 'has no symmetric among the k-points!'))
             end if
             ierr = ierr + 1; cycle
           end if

           eig_kq(:, ik_bz) = ebands%eig(bstart:bstop, ikq_ibz, spin)
         end do
       end do
     end do

     ABI_CHECK(ierr == 0, "See above warnings")

     if (any(ltetra == [1, 2])) then
       ! Compute weights for double delta integration. Note that libtetra assumes Ef set to zero.
       ! TODO: Average weights over degenerate states?
       !write(std_out,"(a,i0,2a)")" Calling libtetrabz_dbldelta with ltetra: ", ltetra, " for q-point:", trim(ktoa(qpt))
       eig_k = eig_k - ebands%fermie; eig_kq = eig_kq - ebands%fermie
       ABI_MALLOC(wght_bz, (nb, nb, nkbz))
       call libtetrabz_dbldelta(ltetra, cryst%gprimd, nb, nge, eig_k, eig_kq, ngw, wght_bz) !, comm=comm)
       eig_k = eig_k + ebands%fermie; eig_kq = eig_kq + ebands%fermie

       wtqs(iq_ibz, spin) = sum(wght_bz)
       ABI_FREE(wght_bz)

     else if (ltetra == 3) then
       ! Tetrahedron method with Allen's approach for double delta.
       !write(std_out,"(2a)")" Calling Allen's version for q-point: ", trim(ktoa(qpt))
       nene = 3
       enemin = ebands%fermie - tol6
       enemax = ebands%fermie + tol6

       ABI_MALLOC(dtweightde, (nkbz, nene, nene))
       ABI_MALLOC(tweight, (nkbz, nene, nene))
       ABI_MALLOC(work_k, (nkbz))
       ABI_MALLOC(work_kq, (nkbz))

       ! TODO
       MSG_ERROR("Not Implemented")
       !call init_tetra(indkpt, cryst%gprimd, fs%klatt, kbz, nkbz, tetra, ierr, errorstring, comm)
       ABI_CHECK(ierr == 0, "init_tetra returned ierr /= 0")

       !fs%dbldelta_tetra_weights_kfs = zero
       do ib2=1,nb
         work_k = eig_k(ib2, :)
         do ib1=1,nb
          work_kq = eig_kq(ib1, :)
          ! TODO: Average weights over degenerate states?
          call get_dbl_tetra_weight(work_k, work_kq, enemin, enemax, enemin, enemax, &
                                    max_occ1, nene, nene, nkbz, tetra, tweight, dtweightde, ierr)
          ABI_CHECK(ierr == 0, "get_dbldelta_weights returned ierr /= 0")

          wtqs(iq_ibz, spin) = wtqs(iq_ibz, spin) + sum(dtweightde(:, 2, 2))
         end do
       end do

       call destroy_tetra(tetra)
       ABI_FREE(work_k)
       ABI_FREE(work_kq)
       ABI_FREE(dtweightde)
       ABI_FREE(tweight)
     else
       MSG_ERROR(sjoin("Invalid value of ltetra:", itoa(ltetra)))
     end if
   end do ! iq_ibz
 end do ! spin

 ABI_FREE(indkpt)
 ABI_FREE(kbz)
 ABI_FREE(eig_k)
 ABI_FREE(eig_kq)

 call ibz_krank%free()

 call xmpi_sum(wtqs, comm, ierr)

 !call cwtime_report(" calc_dbldelta", cpu, wall, gflops)

end subroutine calc_dbldelta
!!***

end module m_phgamma
!!***
