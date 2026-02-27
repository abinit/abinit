!!****m* ABINIT/m_sigma
!! NAME
!!  m_sigma
!!
!! FUNCTION
!!  This module provides the definition of the sigma_t data type
!!  used to store results of the GW calculation.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2026 ABINIT group (MG, FB, GMR, VO, LR, RWG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_sigma

 use, intrinsic :: iso_c_binding
 use defs_basis
 use m_gwdefs
 use m_xmpi
 use m_abicore
 use m_errors
 use m_nctk
 use m_yaml
 use m_melemts
 use netcdf
 use m_wfd

 use defs_abitypes,    only : MPI_type
 !use m_gwdefs,         only : unt_gw, unt_sig, unt_sgr, unt_sgm, unt_gwdiag, sigparams_t, unt_sigc
 use m_fstrings,       only : itoa, sjoin
 use m_numeric_tools,  only : c2r
 use m_crystal,        only : crystal_t
 use m_ebands,         only : ebands_t
 use m_bz_mesh,        only : kmesh_t, littlegroup_t, findqg0
 use m_screening,      only : epsm1_t

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_sigma/sigma_t
!! NAME
!! sigma_t
!!
!! FUNCTION
!! The sigma_t structured datatype gathers the results of a GW calculation.
!!
!! TODO
!!   ragged arrays (nk,nsppol) --> values ?
!!
!! SOURCE

 type,public :: sigma_t

  integer :: b1gw, b2gw     ! min and Max gw band indices over spin and k-points (used to dimension arrays)
  integer :: gwcalctyp      ! Flag defining the calculation type.
  integer :: nkptgw         ! No. of points calculated
  integer :: nkibz          ! No. of irreducible k-points.
  integer :: nbnds          ! Total number of bands
  integer :: nomega_r       ! No. of real frequencies for the spectral function.
  integer :: nomega_i       ! No. of frequencies along the imaginary axis.
  integer :: nomega4sd      ! No. of real frequencies to evaluate the derivative of $\Sigma(E)$.
  integer :: nsig_ab        ! 1 if nspinor=1,4 for noncollinear case.
  integer :: nsppol         ! No. of spin polarizations.
  integer :: usepawu        ! 1 if we are using DFT+U as starting point (only for PAW)

  real(dp) :: deltae       ! Frequency step for the calculation of d\Sigma/dE
  real(dp) :: maxomega4sd  ! Max frequency around E_ks for d\Sigma/dE.
  real(dp) :: maxomega_r   ! Max frequency for spectral function.
  real(dp) :: scissor_ene  ! Scissor energy value. zero for None.

  integer,allocatable :: maxbnd(:,:)
  ! (nkptgw, nsppol)
  ! Max band index considered in GW for this k-point.

  integer,allocatable :: minbnd(:,:)
  ! (nkptgw, nsppol)
  ! Min band index considered in GW for this k-point.

  real(dp),allocatable :: degwgap(:,:)
  ! (nkibz, nsppol)
  ! Difference btw the QP and the KS direct gap.

  real(dp),allocatable :: egwgap(:,:)
  ! (nkibz, nsppol))
  ! QP direct gap at each k-point and spin.

  real(dp),allocatable :: en_qp_diago(:,:,:)
  ! (nbnds, nkibz, nsppol))
  ! QP energies obtained from the diagonalization of the Hermitian approximation to Sigma (QPSCGW)

  real(dp),allocatable :: e0(:,:,:)
  ! (nbnds, nkibz, nsppol)
  ! KS eigenvalues for each band, k-point and spin. In case of self-consistent?

  real(dp),allocatable :: e0gap(:,:)
  ! (nkibz, nsppol),
  ! KS gap at each k-point, for each spin.

  real(dp),allocatable :: omega_r(:)
  ! (nomega_r)
  ! real frequencies used for the self energy.

  real(dp),allocatable :: kptgw(:,:)
  ! (3, nkptgw)
  ! ! TODO there is a similar array in sigparams_t
  ! List of calculated k-points.

  real(dp),allocatable :: sigxme(:,:,:)
  ! (b1gw:b2gw, nkibz, nsppol*nsig_ab))
  ! Diagonal matrix elements $\<nks|\Sigma_x|nks\>$

  real(dp),allocatable :: sigxcnofme(:,:,:)
  ! (b1gw:b2gw, nkibz, nsppol*nsig_ab))
  ! Diagonal matrix elements $\<nks|\Sigma_xc|nks\>$ taking sqrt(occs) in \Sigma_x, occs in [0,1]

  complex(dp),allocatable :: x_mat(:,:,:,:)
  ! (b1gw:b2gw, b1gw:b2gw, nkibz, nsppol*nsig_ab)
  ! Matrix elements of $\<nks|\Sigma_x|mks\>$

  real(dp),allocatable :: vxcme(:,:,:)
  ! (b1gw:b2gw, nkibz, nsppol*nsig_ab))
  ! $\<nks|v_{xc}[n_val]|nks\>$ matrix elements of vxc
  ! NB: valence-only contribution i.e. computed without model core charge

  real(dp),allocatable :: vUme(:,:,:)
  ! (b1gw:b2gw, nkibz, nsppol*nsig_ab))
  ! $\<nks|v_{U}|nks\>$ for DFT+U.

  complex(dp),allocatable :: degw(:,:,:)
  ! (b1gw:b2gw, nkibz, nsppol))
  ! Difference between the QP and the KS energies.

  complex(dp),allocatable :: dsigmee0(:,:,:)
  ! (b1gw:b2gw, nkibz, nsppol*nsig_ab))
  ! Derivative of $\Sigma_c(E)$ calculated at the KS eigenvalue.

  complex(dp),allocatable :: egw(:,:,:)
  ! (nbnds, nkibz, nsppol))
  ! QP energies, $\epsilon_{nks}^{QP}$.

  logical :: needs_eigvec_qp = .True.

 ! FIXME: These arrays are huge and should be allocated only if self-consistent
  complex(dp),allocatable :: eigvec_qp(:,:,:,:)
  ! (nbnds, nbnds, nkibz, nsppol))
  ! Expansion of the QP amplitudes in the QP basis set of the previous iteration.

  complex(dp),allocatable :: m_ks_to_qp(:,:,:,:)
  ! (nbnds, nbnds, nkibz, nsppol))
  ! m_ks_to_qp(ib,jb,k,s) := <\psi_{ib,k,s}^{KS}|\psi_{jb,k,s}^{QP}>

  complex(dp),allocatable :: hhartree(:,:,:,:)
  ! (b1gw:b2gw, b1gw:b2gw, nkibz, nsppol*nsig_ab)
  ! $\<nks|T+v_H+v_{loc}+v_{nl}|mks\>$
  ! Note that v_{loc} does not include the contribution to vxc(r) given by the model core charge.

  complex(dp),allocatable :: sigcme(:,:,:,:)
  ! (b1gw:b2gw, nkibz, nomega_r, nsppol*nsig_ab))
  ! $\<nks|\Sigma_{c}(E)|nks\>$ at each nomega_r frequency

  complex(dp),allocatable :: sigmee(:,:,:)
  ! (b1gw:b2gw, nkibz, nsppol*nsig_ab))
  ! $\Sigma_{xc}E_{KS} + (E_{QP}- E_{KS})*dSigma/dE_KS

  complex(dp),allocatable :: sigcmee0(:,:,:)
  ! (b1gw:b2gw, nkibz, nsppol*nsig_ab))
  ! Diagonal matrix elements of $\Sigma_c(E)$ calculated at the KS energy $E_{KS}$

  complex(dp),allocatable :: sigcmesi(:,:,:,:)
  ! (b1gw:b2gw, nkibz, nomega_i, nsppol*nsig_ab))
  ! Matrix elements of $\Sigma_c$ along the imaginary axis.
  ! Only used in case of analytical continuation.

  complex(dp),allocatable :: sigcme4sd(:,:,:,:)
  ! (b1gw:b2gw, nkibz, nomega4sd, nsppol*nsig_ab))
  ! Diagonal matrix elements of \Sigma_c around the zeroth order eigenvalue (usually KS).

  complex(dp),allocatable :: sigxcme(:,:,:,:)
  ! (b1gw:b2gw, nkibz, nomega_r, nsppol*nsig_ab))
  ! $\<nks|\Sigma_{xc}(E)|nks\>$ at each real frequency frequency.

  complex(dp),allocatable :: sigxcmesi(:,:,:,:)
  ! (b1gw:b2gw, nkibz, nomega_i, nsppol*nsig_ab))
  ! Matrix elements of $\Sigma_{xc}$ along the imaginary axis.
  ! Only used in case of analytical continuation.

  complex(dp),allocatable :: sigxcme4sd(:,:,:,:)
  ! (b1gw:b2gw, nkibz, nomega4sd, nsppol*nsig_ab))
  ! Diagonal matrix elements of \Sigma_xc for frequencies around the zeroth order eigenvalues.

  complex(dp),allocatable :: ze0(:,:,:)
  ! (b1gw:b2gw, nkibz, nsppol))
  ! renormalization factor. $(1-\dfrac{\partial\Sigma_c} {\partial E_{KS}})^{-1}$

  complex(dp),allocatable :: omega_i(:)
  ! (nomega_i)
  ! Frequencies along the imaginary axis used for the analytical continuation.

  complex(dp),allocatable :: omega4sd(:,:,:,:)
  ! (b1gw:b2gw, nkibz, nomega4sd, nsppol).
  ! Frequencies used to evaluate the Derivative of Sigma.

 contains
   procedure :: init => sigma_init
    ! Initialize the object.

   procedure :: free => sigma_free
    ! Deallocate memory.

   procedure :: get_exene => sigma_get_exene
    ! Compute exchange energy.

    procedure :: get_excene => sigma_get_excene
    ! Compute exchange-correlation MBB (Nat. Orb. Funct. Approx.) energy.

    procedure :: get_haene => sigma_get_haene
     ! Compute hartree energy.

    procedure :: get_kiene => sigma_get_kiene
     ! Compute kinetic energy.

    procedure :: ncwrite => sigma_ncwrite
     ! Write data in netcdf format.

    procedure :: write_results => sigma_write_results
    procedure :: print_perturbative => sigma_print_pertubative
    procedure :: print_qpsc => sigma_print_qpsc
 end type sigma_t

 public  :: sigma_distribute_bks
 public ::  write_sigma_header
!!***

contains  !========================================================================================
!!***

!!****f* m_sigma/write_sigma_header
!! NAME
!! write_sigma_header
!!
!! FUNCTION
!!  Write basic info and dimensions used during the calculation
!!  of the QP correctoions (optdriver==4).
!!
!! INPUTS
!!  Sigp=sigparams_t
!!  Cryst<crystal_t>= Info on the Crystal structure
!!  Kmesh<kmesh_t>= Description of the BZ sampling.
!!
!! OUTPUT
!!  (for writing routines, no output) otherwise, should be described
!!
!! NOTES
!!
!! SOURCE

subroutine write_sigma_header(Sigp, epsm1, Cryst, Kmesh, Qmesh)

!Arguments ------------------------------------
!scalars
 class(sigparams_t),intent(in) :: Sigp
 type(kmesh_t),intent(in) :: Kmesh,Qmesh
 type(crystal_t),intent(in) :: Cryst
 type(epsm1_t),intent(in) :: epsm1

!Local variables-------------------------------
!scalars
 integer :: gwcalctyp, mod10
 character(len=500) :: msg
 integer :: units(2)
! *************************************************************************

 units = [std_out, ab_out]
 call wrtout(units, ' SIGMA fundamental parameters:')

 gwcalctyp = Sigp%gwcalctyp
 mod10 = MOD(Sigp%gwcalctyp,10)

 SELECT CASE (mod10)
 CASE (SIG_GW_PPM)
   write(msg,'(a,i2)')' PLASMON POLE MODEL ',Sigp%ppmodel
 CASE (SIG_GW_AC)
   write(msg,'(a)')' ANALYTIC CONTINUATION'
 CASE (SIG_GW_CD)
   write(msg,'(a)')' CONTOUR DEFORMATION'
 CASE (SIG_HF)
   write(msg,'(a)')' Hartree-Fock'
 CASE (SIG_SEX)
   write(msg,'(a)')' Screened Exchange'
 CASE (SIG_COHSEX)
   write(msg,'(a)')' COHSEX'
 CASE (SIG_QPGW_PPM)
   write(msg,'(a,i2)')' MODEL GW with PLASMON POLE MODEL ',Sigp%ppmodel
 CASE (SIG_QPGW_CD)
   write(msg,'(a)')' MODEL GW without PLASMON POLE MODEL'
 CASE DEFAULT
   write(msg,'(a,i0)')' Wrong value for Sigp%gwcalctyp = ',Sigp%gwcalctyp
   ABI_BUG(msg)
 END SELECT
 call wrtout(units, msg)

 write(msg,'(a,i12)')' number of plane-waves for SigmaX         ',Sigp%npwx
 call wrtout(units, msg)
 write(msg,'(a,i12)')' number of plane-waves for SigmaC and W   ',Sigp%npwc
 call wrtout(units, msg)
 write(msg,'(a,i12)')' number of plane-waves for wavefunctions  ',Sigp%npwwfn
 call wrtout(units, msg)
 write(msg,'(a,i12)')' number of bands                          ',Sigp%nbnds
 call wrtout(units, msg)
 write(msg,'(a,i12)')' number of independent spin polarizations ',Sigp%nsppol
 call wrtout(units, msg)
 write(msg,'(a,i12)')' number of spinorial components           ',Sigp%nspinor
 call wrtout(units, msg)
 write(msg,'(a,i12)')' number of k-points in IBZ                ',Kmesh%nibz
 call wrtout(units, msg)
 write(msg,'(a,i12)')' number of q-points in IBZ                ',Qmesh%nibz
 call wrtout(units, msg)
 write(msg,'(a,i12)')' number of symmetry operations            ',Cryst%nsym
 call wrtout(units, msg)
 write(msg,'(a,i12)')' number of k-points in BZ                 ',Kmesh%nbz
 call wrtout(units, msg)
 write(msg,'(a,i12)')' number of q-points in BZ                 ',Qmesh%nbz
 call wrtout(units, msg)
 write(msg,'(a,i12)')' number of frequencies for dSigma/dE      ',Sigp%nomegasrd
 call wrtout(units, msg)
 write(msg,'(a,f12.2)')' frequency step for dSigma/dE [eV]        ',Sigp%deltae*Ha_eV
 call wrtout(units, msg)
 write(msg,'(a,i12)')' number of omega for Sigma on real axis   ',Sigp%nomegasr
 call wrtout(units, msg)
 write(msg,'(a,f12.2)')' max omega for Sigma on real axis  [eV]   ',Sigp%maxomega_r*Ha_eV
 call wrtout(units, msg)
 write(msg,'(a,f12.2)')' zcut for avoiding poles [eV]             ',Sigp%zcut*Ha_eV
 call wrtout(units, msg)

 if (Sigp%mbpt_sciss>0.1d-4) then
   write(msg,'(a,f12.2)')' scissor energy [eV]                      ',Sigp%mbpt_sciss*Ha_eV
   call wrtout(units, msg)
 end if

 if (mod10 == SIG_GW_AC) then
   write(msg,'(a,i12)')' number of imaginary frequencies for Sigma',Sigp%nomegasi
   call wrtout(units, msg)
   ! MRM not needed for GW 1RDM
   if (gwcalctyp/=21) then
     write(msg,'(a,f12.2)')' max omega for Sigma on imag axis  [eV]   ',Sigp%omegasimax*Ha_eV
     call wrtout(units, msg)
   endif
 end if

 if (Sigp%needs_w()) then
   write(msg,'(2a)')ch10,' EPSILON^-1 parameters (SCR file):'
   call wrtout(units, msg)
   write(msg,'(a,i12)')' dimension of the eps^-1 matrix on file   ',epsm1%Hscr%npwe
   call wrtout(units, msg)
   write(msg,'(a,i12)')' dimension of the eps^-1 matrix used      ',epsm1%npwe
   call wrtout(units, msg)
   write(msg,'(a,i12)')' number of plane-waves for wavefunctions  ',epsm1%Hscr%npwwfn_used
   call wrtout(units, msg)
   write(msg,'(a,i12)')' number of bands                          ',epsm1%Hscr%nbnds_used
   call wrtout(units, msg)
   write(msg,'(a,i12)')' number of q-points in IBZ                ',Qmesh%nibz
   call wrtout(units, msg)
   write(msg,'(a,i12)')' number of frequencies                    ',epsm1%nomega
   call wrtout(units, msg)
   write(msg,'(a,i12)')' number of real frequencies               ',epsm1%nomega_r
   call wrtout(units, msg)
   write(msg,'(a,i12)')' number of imag frequencies               ',epsm1%nomega_i
   call wrtout(units, msg)
 end if

  ! MRM not needed for GW 1RDM
  if (gwcalctyp /= 21) then
    write(msg,'(3a)')ch10,' matrix elements of self-energy operator (all in [eV])',ch10
    call wrtout(units, msg)
    !call wrtout(units, "(a)")" Notations:"
    !call wrtout(units, "(a)")"E0: KS eigenvalue.")
    !call wrtout(units, "(a)")"VxcDFT: KS exchange-correlation potential expectation value.")
    !call wrtout(units, "(a)")"SigX: exchange part of the self-energy.")
    !call wrtout(units, "(a)")"SigC(E0) correlation part of the self-energy, evaluated at the KS eigenenergy.")
    !call wrtout(units, "(a)")"Z: renormalization factor.")
    !call wrtout(units, "(a)")"dSigC/dE: energy derivative of SigC with respect to the energy.")
    !call wrtout(units, "(a)")"SigC(E): correlation part of the self-energy, evaluated at the QP energy.")
    !call wrtout(units, "(a)")"E-E0: difference between QP energy and KS eigenenergy.")
    !call wrtout(units, "(a)")"E: quasiparticle energy.")
    !if (mod10 == SIG_GW_AC) then
    !  call wrtout(units, "For AC calculations, the KS Fermi level has been set to zero.")
    !  call wrtout(units, "KS and QP energies are shifted accordingly.")
    !  call wrtout(units, "IMPORTANT: In AC calculations, the QP energies are obtained by solving the non-linear QP equation along the real-axis")
    !end if
  end if

 if (gwcalctyp < 10) then
   write(msg,'(a)')' Perturbative Calculation'
 else if (gwcalctyp < 20) then
   write(msg,'(a)')' Self-Consistent on Energies only'
 else
   write(msg,'(a)')' Self-Consistent on Energies and Wavefunctions'
 end if
 call wrtout(units, msg)

end subroutine write_sigma_header
!!***

!----------------------------------------------------------------------

!!****f* m_sigma/sigma_write_results
!! NAME
!! sigma_write_results
!!
!! FUNCTION
!!  Write the final results of the GW calculation.
!!
!! INPUTS
!!  Sigp=sigparams_t datatype
!!  ik_ibz= index of the k-point in the array kibz, where GW corrections are calculated
!!  ikcalc= index of the k-point in the array Sigp%kptgw2bz
!!  ks_ebands<ebands_t>=Info on the KS band structure energies.
!!
!! SOURCE

subroutine sigma_write_results(sigma, ikcalc, ik_ibz, Sigp, ks_ebands)

!Arguments ------------------------------------
!scalars
 class(sigma_t),intent(in) :: sigma
 integer,intent(in) :: ikcalc,ik_ibz
 type(ebands_t),intent(in) :: ks_ebands
 type(sigparams_t),intent(in) :: Sigp

!Local variables-------------------------------
!scalars
 integer :: ib,io,is,gwcalctyp,mod10
 character(len=500) :: msg
 type(yamldoc_t) :: ydoc
!arrays
 character(len=12) :: tag_spin(2)
! *************************************************************************

 gwcalctyp = Sigp%gwcalctyp
 mod10 = MOD(Sigp%gwcalctyp,10)

 ! unt_gw:  File with GW corrections.
 ! unt_sig: Self-energy as a function of frequency.
 ! unt_sgr: Derivative wrt omega of the Self-energy.
 ! unt_sigc: Sigma_c(eik) MRM
 ! unt_sgm: Sigma on the Matsubara axis (imag axis)

 tag_spin=(/'            ','            '/); if (sigma%nsppol==2) tag_spin=(/',  SPIN UP  ',',  SPIN DOWN'/)

 do is=1,sigma%nsppol
   write(msg,'(2a,3f8.3,a)')ch10,' k = ',Sigp%kptgw(:,ikcalc),tag_spin(is)
   call wrtout(std_out,msg)
   !call wrtout(ab_out,msg)

   msg = ' Band     E0 <VxcDFT>   SigX SigC(E0)      Z dSigC/dE  Sig(E)    E-E0       E'
   if (sigma%usepawu/=0) then
     msg = ' Band     E0 <VxcDFT>   <H_U>  SigX SigC(E0)      Z dSigC/dE  Sig(E)    E-E0       E'
   end if

   if (gwcalctyp>=10) then
     write(msg,'(2a)')&
     ' Band     E_DFT   <VxcDFT>   E(N-1)  <Hhartree>   SigX  SigC[E(N-1)]',&
     '    Z     dSigC/dE  Sig[E(N)]  DeltaE  E(N)_pert E(N)_diago'
   end if
   call wrtout(std_out,msg)

   ydoc = yamldoc_open('SelfEnergy_ee', width=11, real_fmt='(3f8.3)')
   call ydoc%add_real1d('kpoint', Sigp%kptgw(:,ikcalc))
   call ydoc%add_int('spin', is, int_fmt="(i1)")
   call ydoc%add_real('KS_gap', sigma%e0gap(ik_ibz,is)*Ha_eV)
   call ydoc%add_real('QP_gap', sigma%egwgap(ik_ibz,is)*Ha_eV)
   call ydoc%add_real('Delta_QP_KS', sigma%degwgap(ik_ibz,is)*Ha_eV)
   call ydoc%open_tabular('data', tag='SigmaeeData')
   call ydoc%add_tabular_line(msg)

   write(unt_gw,'(3f10.6)')Sigp%kptgw(:,ikcalc)
   write(unt_gw,'(i4)')Sigp%maxbnd(ikcalc,is)-Sigp%minbnd(ikcalc,is)+1

   write(unt_gwdiag,'(3f10.6)')Sigp%kptgw(:,ikcalc)
   write(unt_gwdiag,'(i4)')Sigp%maxbnd(ikcalc,is)-Sigp%minbnd(ikcalc,is)+1

   write(unt_sig,'("# k = ",3f10.6)')Sigp%kptgw(:,ikcalc)
   write(unt_sig,'("# b = ",2i10)')Sigp%minbnd(ikcalc,is),Sigp%maxbnd(ikcalc,is)

   write(unt_sgr,'("# k = ",3f10.6)')Sigp%kptgw(:,ikcalc)
   write(unt_sgr,'("# b = ",2i10)')Sigp%minbnd(ikcalc,is),Sigp%maxbnd(ikcalc,is)

   write(unt_sigc,'("# k = ",3f10.6)')Sigp%kptgw(:,ikcalc)
   write(unt_sigc,'("# b = ",2i10)')Sigp%minbnd(ikcalc,is),Sigp%maxbnd(ikcalc,is)

   do ib=Sigp%minbnd(ikcalc,is),Sigp%maxbnd(ikcalc,is)
     if (gwcalctyp >= 10) then
       call sigma%print_QPSC(ik_ibz, ib, is, ks_ebands, units=[dev_null], ydoc=ydoc)
       call sigma%print_QPSC(ik_ibz, ib, is, ks_ebands, units=[std_out], prtvol=1)

       write(unt_gwdiag,'(i6,3f9.4)')                                 &
        ib,                                                           &
        sigma%en_qp_diago(ib,ik_ibz,is)*Ha_eV,                            &
        (sigma%en_qp_diago(ib,ik_ibz,is) - ks_ebands%eig(ib,ik_ibz,is))*Ha_eV,&
        zero

     else
       ! If not ppmodel, write out also the imaginary part in ab_out
       select case(mod10)
       case (SIG_GW_AC, SIG_GW_CD)
         call sigma%print_perturbative(ik_ibz, ib, is, units=[dev_null], ydoc=ydoc, prtvol=1)
       case default
         call sigma%print_perturbative(ik_ibz, ib, is, units=[dev_null], ydoc=ydoc)
       end select
       call sigma%print_perturbative(ik_ibz, ib, is, units=[std_out], prtvol=1)
     end if

     write(unt_gw,'(i6,3f9.4)')         &
      ib,                               &
      REAL (sigma%egw (ib,ik_ibz,is))*Ha_eV,&
      REAL (sigma%degw(ib,ik_ibz,is))*Ha_eV,&
      AIMAG(sigma%egw (ib,ik_ibz,is))*Ha_eV
   end do !ib

   if (sigma%e0gap(ik_ibz,is)**2+sigma%egwgap(ik_ibz,is)**2+sigma%degwgap(ik_ibz,is)**2 > tol10) then
     ! Output the direct gap for each spin
     ! If all the gaps are zero, this means that they could not be computed in the calling routine
     write(msg,'(2a,f8.3)')ch10,' E^0_gap       ',sigma%e0gap(ik_ibz,is)*Ha_eV
     call wrtout(std_out,msg)
     write(msg,'(a,f8.3)')      ' E^GW_gap      ',sigma%egwgap(ik_ibz,is)*Ha_eV
     call wrtout(std_out,msg)
     write(msg,'(a,f8.3,a)')    ' DeltaE^GW_gap ',sigma%degwgap(ik_ibz,is)*Ha_eV,ch10
     call wrtout(std_out,msg)
   end if

   call ydoc%write_and_free(ab_out)

   ! Output of the spectral function.
   do io=1,sigma%nomega_r
     write(unt_sig,'(100(e12.5,2x))')&
      REAL(sigma%omega_r(io))*Ha_eV,&
      (REAL(sigma%sigxcme(ib,ik_ibz,io,is))*Ha_eV,&
      AIMAG(sigma%sigxcme(ib,ik_ibz,io,is))*Ha_eV,&
      gw_spectral_function(sigma,io,ib,ik_ibz,is),&
      ib=Sigp%minbnd(ikcalc,is),Sigp%maxbnd(ikcalc,is))
   end do

   do ib=Sigp%minbnd(ikcalc,is),Sigp%maxbnd(ikcalc,is)
     write(unt_sgr,'("# ik, ib",2i5)')ik_ibz,ib
     do io=1,sigma%nomega4sd
       write(unt_sgr,'(100(e12.5,2x))')              &
         REAL (sigma%omega4sd  (ib,ik_ibz,io,is)) * Ha_eV,&
         REAL (sigma%sigxcme4sd(ib,ik_ibz,io,is)) * Ha_eV,&
         AIMAG(sigma%sigxcme4sd(ib,ik_ibz,io,is)) * Ha_eV
     end do
   end do

   !MRM
   do ib=Sigp%minbnd(ikcalc,is),Sigp%maxbnd(ikcalc,is)
     write(unt_sigc,'("# ik, ib",2i5)')ik_ibz,ib
     do io=1,sigma%nomega4sd
       write(unt_sigc,'(100(e12.5,2x))')              &
         REAL (sigma%omega4sd  (ib,ik_ibz,io,is)) * Ha_eV,&
         REAL (sigma%sigcme4sd(ib,ik_ibz,io,is))  * Ha_eV,&
         AIMAG(sigma%sigcme4sd(ib,ik_ibz,io,is))  * Ha_eV
     end do
   end do

   if (mod10 == SIG_GW_AC) then
     ! For AC, write sigma matrix elements along the imaginary axis
     do ib=Sigp%minbnd(ikcalc,is),Sigp%maxbnd(ikcalc,is)
       write(unt_sgm,'("# ik, ib",2i5)')ik_ibz,ib
       do io=1,sigma%nomega_i
         write(unt_sgm,'(3(e12.5,2x))')             &
          AIMAG(sigma%omega_i(io))                * Ha_eV,&
          REAL (sigma%sigxcmesi(ib,ik_ibz,io,is)) * Ha_eV,&
          AIMAG(sigma%sigxcmesi(ib,ik_ibz,io,is)) * Ha_eV
       end do
     end do
   end if

 end do !is

end subroutine sigma_write_results
!!***

!----------------------------------------------------------------------

!!****f* m_sigma/gw_spectral_function
!! NAME
!! gw_spectral_function
!!
!! FUNCTION
!!  Compute the spectral function
!!
!! INPUTS
!!  io,ib,ik_ibz,spin=Frequency, band, k-point, spin index
!!
!! SOURCE

real(dp) pure function gw_spectral_function(sigma, io, ib, ik_ibz, spin) result(aw)

!Arguments ------------------------------------
 class(sigma_t),intent(in) :: sigma
 integer,intent(in) :: io,ib,ik_ibz,spin
! *********************************************************************

 aw = one / pi * abs(aimag(sigma%sigcme(ib,ik_ibz,io,spin))) &
   /( (real(sigma%omega_r(io) - sigma%hhartree(ib,ib,ik_ibz,spin) - sigma%sigxcme(ib,ik_ibz,io,spin)))**2 &
     +(aimag(sigma%sigcme(ib,ik_ibz,io,spin))) ** 2) / Ha_eV

end function gw_spectral_function
!!***

!----------------------------------------------------------------------

!!****f* m_sigma/sigma_print_pertubative
!! NAME
!! sigma_print_pertubative
!!
!! FUNCTION
!!  Write the results of the GW calculation done with the perturbative approach
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine sigma_print_pertubative(sigma, ik_ibz, band, spin, units, &
                                   prtvol, with_header, ydoc) ! Optional

!Arguments ------------------------------------
!scalars
 class(sigma_t),intent(in) :: sigma
 integer,intent(in) :: band,ik_ibz,spin
 integer,optional,intent(in) :: prtvol,units(:)
 logical,optional,intent(in) :: with_header
 type(yamldoc_t),intent(inout),optional :: ydoc

!Local variables-------------------------------
!scalars
 integer :: verbose
 character(len=500) :: msg
! *********************************************************************

 verbose=0      ; if (PRESENT(prtvol)) verbose=prtvol

 if (present(with_header)) then
   if (with_header) then
     call wrtout(units,' Band     E0 <VxcDFT>   SigX SigC(E0)      Z dSigC/dE  Sig(E)    E-E0       E ')
   end if
 end if

 if (sigma%usepawu == 0) then

   if (sigma%nsig_ab /= 1) then
     write(msg,'(i5,9f8.3)')                       &
           band,                                  &
           sigma%e0          (band,ik_ibz,1)*Ha_eV,  &
           SUM(sigma%vxcme   (band,ik_ibz,:))*Ha_eV, &
           SUM(sigma%sigxme  (band,ik_ibz,:))*Ha_eV, &
      REAL(SUM(sigma%sigcmee0(band,ik_ibz,:)))*Ha_eV,&
      REAL(sigma%ze0         (band,ik_ibz,1)),       &
      REAL(SUM(sigma%dsigmee0(band,ik_ibz,:))),      &
      REAL(SUM(sigma%sigmee  (band,ik_ibz,:)))*Ha_eV,&
      REAL(sigma%degw        (band,ik_ibz,1))*Ha_eV, &
      REAL(sigma%egw         (band,ik_ibz,1))*Ha_eV
     call wrtout(units, msg)
     if (present(ydoc)) call ydoc%add_tabular_line(msg)
     if (verbose /= 0) then
       write(msg,'(i5,9f8.3)')                        &
              band,                                   &
              zero,                                   &
              zero,                                   &
              zero,                                   &
        AIMAG(SUM(sigma%sigcmee0(band,ik_ibz,:)))*Ha_eV,&
        AIMAG(sigma%ze0         (band,ik_ibz,1)),       &
        AIMAG(SUM(sigma%dsigmee0(band,ik_ibz,:))),      &
        AIMAG(SUM(sigma%sigmee  (band,ik_ibz,:)))*Ha_eV,&
        AIMAG(sigma%degw        (band,ik_ibz,1))*Ha_eV, &
        AIMAG(sigma%egw         (band,ik_ibz,1))*Ha_eV
       call wrtout(units, msg)
       if(present(ydoc)) call ydoc%add_tabular_line(msg)
     end if
  else
    write(msg,'(i5,9f8.3)')                    &
          band,                                &
          sigma%e0      (band,ik_ibz,spin)*Ha_eV, &
          sigma%vxcme   (band,ik_ibz,spin)*Ha_eV, &
          sigma%sigxme  (band,ik_ibz,spin)*Ha_eV, &
     REAL(sigma%sigcmee0(band,ik_ibz,spin))*Ha_eV,&
     REAL(sigma%ze0     (band,ik_ibz,spin)),      &
     REAL(sigma%dsigmee0(band,ik_ibz,spin)),      &
     REAL(sigma%sigmee  (band,ik_ibz,spin))*Ha_eV,&
     REAL(sigma%degw    (band,ik_ibz,spin))*Ha_eV,&
     REAL(sigma%egw     (band,ik_ibz,spin))*Ha_eV
    call wrtout(units, msg)
    if (present(ydoc)) call ydoc%add_tabular_line(msg)

    if (verbose /= 0) then
      write(msg,'(i5,9f8.3)')                      &
              band,                                &
              zero,                                &
              zero,                                &
              zero,                                &
        AIMAG(sigma%sigcmee0(band,ik_ibz,spin))*Ha_eV,&
        AIMAG(sigma%ze0     (band,ik_ibz,spin)),      &
        AIMAG(sigma%dsigmee0(band,ik_ibz,spin)),      &
        AIMAG(sigma%sigmee  (band,ik_ibz,spin))*Ha_eV,&
        AIMAG(sigma%degw    (band,ik_ibz,spin))*Ha_eV,&
        AIMAG(sigma%egw     (band,ik_ibz,spin))*Ha_eV
       call wrtout(units,msg)
       if (present(ydoc)) call ydoc%add_tabular_line(msg)
    end if
  end if

 else
   ! PAW+U+GW calculation.
   ABI_CHECK(sigma%nsig_ab==1, 'DFT+U with spinor not implemented')
   write(msg,'(i5,10f8.3)')                   &
         band,                                &
         sigma%e0      (band,ik_ibz,spin)*Ha_eV, &
         sigma%vxcme   (band,ik_ibz,spin)*Ha_eV, &
         sigma%vUme    (band,ik_ibz,spin)*Ha_eV, &
         sigma%sigxme  (band,ik_ibz,spin)*Ha_eV, &
    REAL(sigma%sigcmee0(band,ik_ibz,spin))*Ha_eV,&
    REAL(sigma%ze0     (band,ik_ibz,spin)),      &
    REAL(sigma%dsigmee0(band,ik_ibz,spin)),      &
    REAL(sigma%sigmee  (band,ik_ibz,spin))*Ha_eV,&
    REAL(sigma%degw    (band,ik_ibz,spin))*Ha_eV,&
    REAL(sigma%egw     (band,ik_ibz,spin))*Ha_eV
   call wrtout(units,msg)
   if(present(ydoc)) call ydoc%add_tabular_line(msg)

   if (verbose/=0) then
     write(msg,'(i5,10f8.3)')                   &
           band,                               &
           zero,                                &
           zero,                                &
           zero,                                &
           zero,                                &
     AIMAG(sigma%sigcmee0(band,ik_ibz,spin))*Ha_eV,&
     AIMAG(sigma%ze0     (band,ik_ibz,spin)),      &
     AIMAG(sigma%dsigmee0(band,ik_ibz,spin)),      &
     AIMAG(sigma%sigmee  (band,ik_ibz,spin))*Ha_eV,&
     AIMAG(sigma%degw    (band,ik_ibz,spin))*Ha_eV,&
     AIMAG(sigma%egw     (band,ik_ibz,spin))*Ha_eV
     call wrtout(units, msg)
     if(present(ydoc)) call ydoc%add_tabular_line(msg)
   end if
 end if

end subroutine sigma_print_pertubative
!!***

!----------------------------------------------------------------------

!!****f* m_sigma/sigma_print_qpsc
!! NAME
!!  sigma_print_qpsc
!!
!! FUNCTION
!!  Write the results of the GW calculation in case of self-consistency
!!
!! SOURCE

subroutine sigma_print_qpsc(sigma, ik_ibz, band, spin, ks_ebands, units, &
                            prtvol, ydoc) ! Optional

!Arguments ------------------------------------
!scalars
 class(sigma_t),intent(in) :: sigma
 integer,intent(in) :: band,ik_ibz,spin
 integer,intent(in),optional :: prtvol, units(:)
 type(ebands_t),intent(in) :: ks_ebands
 type(yamldoc_t),intent(inout),optional :: ydoc

!Local variables-------------------------------
!scalars
 integer :: verbose
 character(len=500) :: msg
! *********************************************************************

 verbose=0      ; if (PRESENT(prtvol    )) verbose=prtvol

! write(msg,'(a)')&
!&   ' Band     E_DFT   <VxcDFT>   E(N-1)  <Hhartree>   SigX  SigC[E(N-1)]',&
!&   '    Z     dSigC/dE  Sig[E(N)]  DeltaE  E(N)_pert E(N)_diago'

 if (sigma%usepawu==0 .or. .TRUE.) then
   if (sigma%nsig_ab/=1) then
     write(msg,'(i5,12(2x,f8.3))')                       &
           band,                                         &
           ks_ebands%eig     (band,ik_ibz,1)*Ha_eV,        &
           SUM(sigma%vxcme   (band,ik_ibz,:))*Ha_eV,       &
           sigma%e0          (band,ik_ibz,1)*Ha_eV,        &
      REAL(SUM(sigma%hhartree(band,band,ik_ibz,:)))*Ha_eV,&
           SUM(sigma%sigxme  (band,ik_ibz,:))*Ha_eV,       &
      REAL(SUM(sigma%sigcmee0(band,ik_ibz,:)))*Ha_eV,      &
      REAL(sigma%ze0         (band,ik_ibz,1)),             &
      REAL(SUM(sigma%dsigmee0(band,ik_ibz,:))),            &
      REAL(SUM(sigma%sigmee  (band,ik_ibz,:)))*Ha_eV,      &
      REAL(sigma%degw        (band,ik_ibz,1))*Ha_eV,       &
      REAL(sigma%egw         (band,ik_ibz,1))*Ha_eV,       &
           sigma%en_qp_diago (band,ik_ibz,1)*Ha_eV
     call wrtout(units, msg)
     if (present(ydoc)) call ydoc%add_tabular_line(msg)

     write(msg,'(i5,12(2x,f8.3))')                        &
            band,                                         &
            zero,                                         &
            zero,                                         &
            zero,                                         &
      AIMAG(SUM(sigma%hhartree(band,band,ik_ibz,:)))*Ha_eV,&
            zero,                                           &
      AIMAG(SUM(sigma%sigcmee0(band,ik_ibz,:)))*Ha_eV,      &
      AIMAG(sigma%ze0         (band,ik_ibz,1)),             &
      AIMAG(SUM(sigma%dsigmee0(band,ik_ibz,:))),            &
      AIMAG(SUM(sigma%sigmee  (band,ik_ibz,:)))*Ha_eV,      &
      AIMAG(sigma%degw        (band,ik_ibz,1))*Ha_eV,       &
      AIMAG(sigma%egw         (band,ik_ibz,1))*Ha_eV,       &
            zero
     if (verbose/=0) then
       call wrtout(units, msg)
       if (present(ydoc)) call ydoc%add_tabular_line(msg)
     end if
   else
     write(msg,'(i5,12(2x,f8.3))')                          &
           band,                                            &
           ks_ebands%eig    (band,ik_ibz,spin)*Ha_eV,       &
           sigma%vxcme      (band,ik_ibz,spin)*Ha_eV,       &
           sigma%e0         (band,ik_ibz,spin)*Ha_eV,       &
      REAL(sigma%hhartree   (band,band,ik_ibz,spin))*Ha_eV, &
           sigma%sigxme     (band,ik_ibz,spin)*Ha_eV,       &
      REAL(sigma%sigcmee0   (band,ik_ibz,spin))*Ha_eV,      &
      REAL(sigma%ze0        (band,ik_ibz,spin)),            &
      REAL(sigma%dsigmee0   (band,ik_ibz,spin)),            &
      REAL(sigma%sigmee     (band,ik_ibz,spin))*Ha_eV,      &
      REAL(sigma%degw       (band,ik_ibz,spin))*Ha_eV,      &
      REAL(sigma%egw        (band,ik_ibz,spin))*Ha_eV,      &
           sigma%en_qp_diago(band,ik_ibz,spin)*Ha_eV
     call wrtout(units, msg)
     if (present(ydoc)) call ydoc%add_tabular_line(msg)

     write(msg,'(i5,12(2x,f8.3))')                       &
            band,                                        &
            zero,                                        &
            zero,                                        &
            zero,                                        &
      AIMAG(sigma%hhartree  (band,band,ik_ibz,spin))*Ha_eV,&
            zero,                                           &
      AIMAG(sigma%sigcmee0   (band,ik_ibz,spin))*Ha_eV,     &
      AIMAG(sigma%ze0        (band,ik_ibz,spin)),           &
      AIMAG(sigma%dsigmee0   (band,ik_ibz,spin)),           &
      AIMAG(sigma%sigmee     (band,ik_ibz,spin))*Ha_eV,     &
      AIMAG(sigma%degw       (band,ik_ibz,spin))*Ha_eV,     &
      AIMAG(sigma%egw        (band,ik_ibz,spin))*Ha_eV,     &
            zero
     if (verbose/=0) then
       call wrtout(units, msg)
       if (present(ydoc)) call ydoc%add_tabular_line(msg)
     end if
   end if

 else
   ! PAW+U+GW calculation.
   ABI_ERROR("PAW+U+GW not yet implemented")
 end if

end subroutine sigma_print_qpsc
!!***

!----------------------------------------------------------------------

!!****f* m_sigma/sigma_init
!! NAME
!! sigma_init
!!
!! FUNCTION
!! Main creation method for the sigma_t data type.
!!
!! INPUTS
!! usepawu= /=0 if we used DFT+U as starting point (only for PAW)
!!
!! SOURCE

subroutine sigma_init(sigma, Sigp, nkibz, usepawu)

!Arguments ------------------------------------
 class(sigma_t),intent(inout) :: sigma
 type(sigparams_t),intent(in) :: Sigp
 integer,intent(in) :: nkibz, usepawu

!Local variables-------------------------------
 integer :: b1gw,b2gw,mod10
! *************************************************************************

 mod10 = MOD(Sigp%gwcalctyp,10)

 ! Copy important dimensions
 sigma%nkptgw     =Sigp%nkptgw
 sigma%gwcalctyp  =Sigp%gwcalctyp
 sigma%deltae     =Sigp%deltae
 sigma%maxomega4sd=Sigp%maxomega4sd
 sigma%maxomega_r =Sigp%maxomega_r
 sigma%scissor_ene=Sigp%mbpt_sciss

 !FIXME this should be done in sigma_allocate
 ABI_MALLOC(sigma%minbnd, (sigma%nkptgw,Sigp%nsppol))
 ABI_MALLOC(sigma%maxbnd, (sigma%nkptgw,Sigp%nsppol))
 sigma%minbnd=Sigp%minbnd; sigma%maxbnd=Sigp%maxbnd
 ABI_MALLOC(sigma%kptgw, (3,sigma%nkptgw))
 sigma%kptgw=Sigp%kptgw

 sigma%b1gw     =Sigp%minbdgw ! min and Max GW band index over k and spin.
 sigma%b2gw     =Sigp%maxbdgw ! Used to dimension arrays.
 sigma%nbnds    =Sigp%nbnds
 sigma%nkibz    =nkibz
 sigma%nsppol   =Sigp%nsppol
 sigma%nsig_ab  =Sigp%nsig_ab
 sigma%nomega_r =Sigp%nomegasr  !FIXME change name
 sigma%nomega_i =Sigp%nomegasi
 sigma%nomega4sd=Sigp%nomegasrd
 sigma%usepawu  =usepawu

 !================================================
 ! === Allocate arrays in the sigma_t datatype ===
 !================================================
 b1gw=sigma%b1gw
 b2gw=sigma%b2gw

 ! hhartree(b1,b2,k,s)= <b1,k,s|T+v_{loc}+v_{nl}+v_{H}|b2,k,s>
 ABI_CALLOC(sigma%hhartree, (b1gw:b2gw,b1gw:b2gw,sigma%nkibz,sigma%nsppol*sigma%nsig_ab))

 ! QP amplitudes and energies.
 ABI_CALLOC(sigma%en_qp_diago, (sigma%nbnds,sigma%nkibz,sigma%nsppol))

 sigma%needs_eigvec_qp = sigp%gwcalctyp >= 10
 if (sigma%needs_eigvec_qp) then
   ABI_CALLOC(sigma%eigvec_qp, (sigma%nbnds,sigma%nbnds,sigma%nkibz,sigma%nsppol))
 end if

 ABI_CALLOC(sigma%vxcme, (b1gw:b2gw, sigma%nkibz, sigma%nsppol*sigma%nsig_ab))
 ABI_CALLOC(sigma%vUme, (b1gw:b2gw, sigma%nkibz, sigma%nsppol*sigma%nsig_ab))
 ABI_CALLOC(sigma%sigxme, (b1gw:b2gw, sigma%nkibz, sigma%nsppol*sigma%nsig_ab))
 ABI_CALLOC(sigma%sigxcnofme, (b1gw:b2gw, sigma%nkibz, sigma%nsppol*sigma%nsig_ab))
 ABI_CALLOC(sigma%x_mat, (b1gw:b2gw, b1gw:b2gw, sigma%nkibz, sigma%nsppol*sigma%nsig_ab))
 ABI_CALLOC(sigma%sigcme, (b1gw:b2gw, sigma%nkibz, sigma%nomega_r, sigma%nsppol*sigma%nsig_ab))
 ABI_CALLOC(sigma%sigxcme, (b1gw:b2gw, sigma%nkibz, sigma%nomega_r, sigma%nsppol*sigma%nsig_ab))
 ABI_CALLOC(sigma%sigcmee0, (b1gw:b2gw, sigma%nkibz, sigma%nsppol*sigma%nsig_ab))
 ABI_CALLOC(sigma%ze0, (b1gw:b2gw, sigma%nkibz, sigma%nsppol))
 ABI_CALLOC(sigma%dsigmee0, (b1gw:b2gw, sigma%nkibz, sigma%nsppol*sigma%nsig_ab))
 ABI_CALLOC(sigma%sigmee, (b1gw:b2gw, sigma%nkibz, sigma%nsppol*sigma%nsig_ab))
 ABI_CALLOC(sigma%degw, (b1gw:b2gw, sigma%nkibz, sigma%nsppol))
 ABI_CALLOC(sigma%e0, (sigma%nbnds, sigma%nkibz, sigma%nsppol))
 ABI_CALLOC(sigma%egw, (sigma%nbnds, sigma%nkibz, sigma%nsppol))
 ABI_CALLOC(sigma%e0gap, (sigma%nkibz, sigma%nsppol))
 ABI_CALLOC(sigma%degwgap, (sigma%nkibz, sigma%nsppol))
 ABI_CALLOC(sigma%egwgap, (sigma%nkibz, sigma%nsppol))

 ! These quantities are used to evaluate $\Sigma(E)$ around the KS\QP eigenvalue
 ABI_CALLOC(sigma%omega4sd, (b1gw:b2gw, sigma%nkibz, sigma%nomega4sd, sigma%nsppol))
 ABI_CALLOC(sigma%sigcme4sd, (b1gw:b2gw, sigma%nkibz, sigma%nomega4sd, sigma%nsppol*sigma%nsig_ab))
 ABI_CALLOC(sigma%sigxcme4sd, (b1gw:b2gw, sigma%nkibz, sigma%nomega4sd, sigma%nsppol*sigma%nsig_ab))

 ! Mesh along the real axis.
 !TODO Find  better treatment
 if (sigma%nomega_r > 0) then
   ABI_MALLOC(sigma%omega_r, (sigma%nomega_r))
   sigma%omega_r(:)=Sigp%omega_r(:)
 end if

 ! Analytic Continuation
 ! FIXME omegasi should not be in Sigp% here we should construct the mesh
 if (mod10 == SIG_GW_AC) then
   ABI_MALLOC(sigma%omega_i, (sigma%nomega_i))
   sigma%omega_i = Sigp%omegasi
   ABI_CALLOC(sigma%sigcmesi,  (b1gw:b2gw, sigma%nkibz, sigma%nomega_i, sigma%nsppol*sigma%nsig_ab))
   ABI_CALLOC(sigma%sigxcmesi, (b1gw:b2gw, sigma%nkibz, sigma%nomega_i, sigma%nsppol*sigma%nsig_ab))
 end if

end subroutine sigma_init
!!***

!----------------------------------------------------------------------

!!****f* m_sigma/sigma_free
!! NAME
!! sigma_free
!!
!! FUNCTION
!!  Deallocate all associated pointers defined in the sigma_t data type.
!!
!! SOURCE

subroutine sigma_free(sigma)

!Arguments ------------------------------------
 class(sigma_t),intent(inout) :: sigma
! *************************************************************************

 ! integer
 ABI_SFREE(sigma%maxbnd)
 ABI_SFREE(sigma%minbnd)

 ! real
 ABI_SFREE(sigma%degwgap)
 ABI_SFREE(sigma%egwgap)
 ABI_SFREE(sigma%en_qp_diago)
 ABI_SFREE(sigma%e0)
 ABI_SFREE(sigma%e0gap)
 ABI_SFREE(sigma%omega_r)
 ABI_SFREE(sigma%kptgw)
 ABI_SFREE(sigma%sigxme)
 ABI_SFREE(sigma%sigxcnofme)
 ABI_SFREE(sigma%x_mat)
 ABI_SFREE(sigma%vxcme)
 ABI_SFREE(sigma%vUme)

 ! complex
 ABI_SFREE(sigma%degw)
 ABI_SFREE(sigma%dsigmee0)
 ABI_SFREE(sigma%egw)
 ABI_SFREE(sigma%eigvec_qp)
 ABI_SFREE(sigma%m_ks_to_qp)
 ABI_SFREE(sigma%hhartree)
 ABI_SFREE(sigma%sigcme)
 ABI_SFREE(sigma%sigmee)
 ABI_SFREE(sigma%sigcmee0)
 ABI_SFREE(sigma%sigcmesi)
 ABI_SFREE(sigma%sigcme4sd)
 ABI_SFREE(sigma%sigxcme)
 ABI_SFREE(sigma%sigxcmesi)
 ABI_SFREE(sigma%sigxcme4sd)
 ABI_SFREE(sigma%ze0)
 ABI_SFREE(sigma%omega_i)
 ABI_SFREE(sigma%omega4sd)

end subroutine sigma_free
!!***

!----------------------------------------------------------------------

!!****f* m_sigma/sigma_get_exene
!! NAME
!!  sigma_get_exene
!!
!! FUNCTION
!!  Compute exchange energy.
!!
!! INPUTS
!!  sigma<sigma_t>=Sigma results
!!  kmesh<kmesh_t>=BZ sampling.
!!  bands<band_t>=Bands with occupation factors
!!
!! SOURCE

real(dp) pure function sigma_get_exene(sigma, kmesh, bands) result(ex_energy)

!Arguments ------------------------------------
 class(sigma_t),intent(in) :: sigma
 type(kmesh_t),intent(in) :: kmesh
 type(ebands_t),intent(in) :: bands

!Local variables-------------------------------
!scalars
 integer :: ik,ib,spin
 real(dp) :: wtk,occ_bks
! *************************************************************************

 ex_energy = zero

 do spin=1,sigma%nsppol
   do ik=1,sigma%nkibz
     wtk = kmesh%wt(ik)
     do ib=sigma%b1gw,sigma%b2gw
       occ_bks = bands%occ(ib,ik,spin)
       if (sigma%nsig_ab == 1) then
         ex_energy = ex_energy + half * occ_bks * wtk * sigma%sigxme(ib,ik,spin)
       else
         ex_energy = ex_energy + half * occ_bks * wtk * SUM(sigma%sigxme(ib,ik,:))
       end if
     end do
   end do
 end do

end function sigma_get_exene
!!***

!----------------------------------------------------------------------

!!****f* m_sigma/sigma_get_excene
!! NAME
!!  sigma_get_excene
!!
!! FUNCTION
!!  Compute exchange correlation energy using MBB (nat. orb. functional approx.).
!!
!! INPUTS
!!  sigma<sigma_t>=Sigma results
!!  kmesh<kmesh_t>=BZ sampling.
!!  bands<band_t>=Bands with occupation factors
!!
!! SOURCE

real(dp) pure function sigma_get_excene(sigma, kmesh, bands) result(exc_energy)

!Arguments ------------------------------------
 class(sigma_t),intent(in) :: sigma
 type(kmesh_t),intent(in) :: kmesh
 type(ebands_t),intent(in) :: bands

!Local variables-------------------------------
!scalars
 integer :: ik, ib, spin
 real(dp) :: wtk, occ_bks
! *************************************************************************

 exc_energy = zero

 do spin=1,sigma%nsppol
   do ik=1,sigma%nkibz
     wtk = kmesh%wt(ik)
     do ib=sigma%b1gw,sigma%b2gw
       occ_bks = bands%occ(ib,ik,spin)
       if (sigma%nsig_ab==1) then
         if (sigma%nsppol==1) then
           exc_energy = exc_energy + sqrt( abs( half * occ_bks ) ) * wtk * sigma%sigxcnofme(ib,ik,spin)   ! 2*sqrt(occ_i), occ in [0,2] -> [0,1].
         else
           exc_energy = exc_energy + half * sqrt( abs( occ_bks ) ) * wtk * sigma%sigxcnofme(ib,ik,spin)   ! 2*sqrt(occ_i), occ in [0,1] -> [0,1].
         end if
       else
         exc_energy = exc_energy + half * sqrt( abs( occ_bks ) ) * wtk * SUM(sigma%sigxcnofme(ib,ik,:)) ! 2*sqrt(occ_i), occ in [0,1].
       end if
     end do
   end do
 end do

end function sigma_get_excene
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/sigma_get_haene
!! NAME
!! sigma_get_haene
!!
!! FUNCTION
!! Compute the Hartree energy
!!
!! INPUTS
!! Kmesh <kmesh_t>=Structure describing the k-point sampling.
!! bands=<ebands_t>=Datatype gathering info on the QP energies (KS if one shot)
!!  eig(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=KS or QP energies for k-points, bands and spin
!!  occ(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=occupation numbers, for each k point in IBZ, each band and spin
!! Mels
!!  %vhartr=matrix elements of $v_H$.
!!
!! OUTPUT
!! Compute the Hartree energy on eh_energy
!!
!! SOURCE

real(dp) pure function sigma_get_haene(sigma, Mels, kmesh, bands) result(eh_energy)

!Arguments ------------------------------------
!scalars
 class(sigma_t),intent(in) :: sigma
 type(kmesh_t),intent(in) :: kmesh
 type(ebands_t),intent(in) :: bands
 type(melements_t),intent(in) :: Mels

!Local variables-------------------------------
 integer :: ik,ib,spin
 real(dp) :: wtk,occ_bks
! *************************************************************************

 eh_energy = zero

 do spin=1,sigma%nsppol
   do ik=1,sigma%nkibz
     wtk = kmesh%wt(ik)
     do ib=sigma%b1gw,sigma%b2gw
       occ_bks = bands%occ(ib,ik,spin)
       if (sigma%nsig_ab == 1) then ! Only closed-shell restricted is programed
         eh_energy = eh_energy + occ_bks * wtk * Mels%vhartree(ib,ib,ik,spin)
       end if
     end do
   end do
 end do

 eh_energy = half * eh_energy

end function sigma_get_haene
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/sigma_get_kiene
!! NAME
!! sigma_get_kiene
!!
!! FUNCTION
!! Compute the kinetic energy
!!
!! INPUTS
!! Kmesh <kmesh_t>=Structure describing the k-point sampling.
!! bands=<ebands_t>=Datatype gathering info on the QP energies (KS if one shot)
!!  eig(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=KS or QP energies for k-points, bands and spin
!!  occ(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=occupation numbers, for each k point in IBZ, each band and spin
!! Mels
!!  %kinetic=matrix elements of $T$.
!!
!! OUTPUT
!! Compute the kinetic energy on ek_energy
!!
!! SOURCE

real(dp) pure function sigma_get_kiene(sigma, Mels, kmesh, bands) result(ek_energy)

!Arguments ------------------------------------
 class(sigma_t),intent(in) :: sigma
 type(kmesh_t),intent(in) :: kmesh
 type(ebands_t),intent(in) :: bands
 type(melements_t),intent(in) :: Mels

!Local variables-------------------------------
!scalars
 integer :: ik, ib, spin
 real(dp) :: wtk, occ_bks
! *************************************************************************

 ek_energy = zero

 do spin=1,sigma%nsppol
   do ik=1,sigma%nkibz
     wtk = kmesh%wt(ik)
     do ib=sigma%b1gw,sigma%b2gw
       occ_bks = bands%occ(ib,ik,spin)
       if (sigma%nsig_ab==1) then ! Only closed-shell restricted is programed
         ek_energy = ek_energy + occ_bks * wtk * Mels%kinetic(ib,ib,ik,spin)
       end if
     end do
   end do
 end do

end function sigma_get_kiene
!!***

!----------------------------------------------------------------------

!!****f* m_sigma/find_wpoles_for_cd
!! NAME
!!  find_wpoles_for_cd
!!
!! FUNCTION
!!  Find the max frequency needed to account for all the poles of the GW self-energy
!!  in the case of contour deformation technique.
!!
!! INPUTS
!!  Sigp=sigparams_t
!!
!! OUTPUT
!!  omega_max
!!
!! SOURCE

subroutine find_wpoles_for_cd(Sigp, sigma, Kmesh, ebands, omega_max)

!Arguments ------------------------------------
!scalars
 class(sigparams_t),intent(in) :: Sigp
 type(sigma_t),intent(in) :: sigma
 type(ebands_t),intent(in) :: ebands
 type(kmesh_t),intent(in) :: Kmesh
 real(dp),intent(out) :: omega_max

!Local variables-------------------------------
!scalars
 integer :: spin,ik_ibz,band_gr,bgw_start,bgw_stop,io,ioe0j
 integer :: ikgw,ikgw_ibz,ikgw_bz,band_gw,nomega_tot
 real(dp) :: e_green,e_screen,theta_mu_minus_e0i,e_qp,fact_sp
 !character(len=500) :: msg
!arrays
 real(dp),allocatable :: omegame0i(:)
! *************************************************************************

 omega_max = smallest_real
 !
 ! === Normalization of theta_mu_minus_e0i ===
 ! * If nsppol==2, qp_occ $\in [0,1]$
 fact_sp=one
 if (ebands%nsppol==1) then
   fact_sp=half; if (ebands%nspinor==2) fact_sp=one
 end if
 !
 ! Total number of frequencies for sigma (Spectral function + mesh for the derivative).
 nomega_tot=sigma%nomega_r+sigma%nomega4sd
 ABI_MALLOC(omegame0i,(nomega_tot))

 ioe0j=sigma%nomega4sd/2+1
 !
 ! Loop over bands used to construct the Green function.
 do spin=1,ebands%nsppol
   do ik_ibz=1,ebands%nkpt
     do band_gr=1,ebands%nband(ik_ibz+(spin-1)*ebands%nkpt)
       e_green           = ebands%eig(band_gr,ik_ibz,spin)
       theta_mu_minus_e0i= ebands%occ(band_gr,ik_ibz,spin)*fact_sp
       !
       ! Loop over GW states.
       do ikgw=1,Sigp%nkptgw
         bgw_start=Sigp%minbnd(ikgw,spin)
         bgw_stop =Sigp%minbnd(ikgw,spin)
         ikgw_bz  =Sigp%kptgw2bz(ikgw_bz)
         ikgw_ibz =Kmesh%tab(ikgw_bz)

         do band_gw=bgw_start,bgw_stop
           e_qp      = ebands%eig(band_gw,ikgw_ibz,spin)
           !
           ! Get frequencies $\omega$-\epsilon_in$ to evaluate  $d\Sigma/dE$, note the spin
           ! subtract e_KS since we have stored e_KS+ Delta \omega in sigma%omega4sd, not required for AC
           if (sigma%nomega_r>0) omegame0i(1:sigma%nomega_r)=DBLE(Sigp%omega_r(1:sigma%nomega_r))-e_green
           do io=sigma%nomega_r+1,nomega_tot
             !omegame0i(io)=DBLE(sigma%omega4sd(band_gw,ikgw_ibz,io-sigma%nomega_r,spin)) - e_green
             !sigma%omega4sd(jb,ik_ibz,io,spin)=sigma%egw(jb,ik_ibz,spin)+Sigp%deltae*(io-ioe0j)
             omegame0i(io) = e_qp + Sigp%deltae*(io-ioe0j) - e_green
           end do

           do io=1,nomega_tot
             e_screen =  ABS(omegame0i(io))
             if (omegame0i(io)>tol12) then
               !ket(spadc+ig,ios)=ket(spadc+ig,ios)+ct*(one-theta_mu_minus_e0i)
               if ( (one-theta_mu_minus_e0i) > tol12 ) omega_max = MAX(omega_max, e_screen)
             end if
             if (omegame0i(io)<-tol12) then
               !ket(spadc+ig,ios)=ket(spadc+ig,ios)-ct*theta_mu_minus_e0i
               if ( theta_mu_minus_e0i > tol12) omega_max = MAX(omega_max, e_screen)
             end if
           end do

         end do
       end do
       !
     end do
   end do
 end do

 ABI_FREE(omegame0i)

end subroutine find_wpoles_for_cd
!!***

!----------------------------------------------------------------------

!!****f* m_sigma/sigma_ncwrite
!! NAME
!! sigma_ncwrite
!!
!! FUNCTION
!!  Save the data stored in the sigma_t data type on a NETCDF file.
!!
!! INPUTS
!!  filename
!!
!! OUTPUT
!!
!! SOURCE

integer function sigma_ncwrite(sigma, Sigp, epsm1, ncid) result (ncerr)

!Arguments ------------------------------------
!scalars
 class(sigma_t),target,intent(in) :: sigma
 class(sigparams_t),target,intent(in) :: Sigp
 integer,intent(in) :: ncid
 type(epsm1_t),target,intent(in) :: epsm1

!Local variables ---------------------------------------
!scalars
 integer :: nbgw,ndim_sig,b1gw,b2gw,cplex
 !character(len=500) :: msg
!arrays
 real(dp),allocatable :: rdata2(:,:),rdata4(:,:,:,:),rdata5(:,:,:,:,:)
! *************************************************************************

 !@sigma_t
 cplex=2; b1gw=sigma%b1gw; b2gw=sigma%b2gw; nbgw=b2gw-b1gw+1
 ndim_sig=sigma%nsppol*sigma%nsig_ab

 ncerr = nctk_def_dims(ncid, [&
   nctkdim_t("cplex", cplex), nctkdim_t("b1gw", sigma%b1gw), nctkdim_t("b2gw", sigma%b2gw),&
   nctkdim_t("nbgw", nbgw), nctkdim_t("nkptgw", sigma%nkptgw), nctkdim_t("ndim_sig", ndim_sig), &
   nctkdim_t("nomega4sd", sigma%nomega4sd), nctkdim_t("nsig_ab", sigma%nsig_ab), &
   nctkdim_t("scr_nomega", epsm1%nomega) &
 ], defmode=.True.)
 NCF_CHECK(ncerr)

 ! No. of real frequencies, might be zero.
 if (sigma%nomega_r > 0) then
   NCF_CHECK(nctk_def_dims(ncid, nctkdim_t("nomega_r", sigma%nomega_r)))
 end if

 ! No. of imaginary frequencies, might be zero.
 if (sigma%nomega_i > 0) then
   NCF_CHECK(nctk_def_dims(ncid, nctkdim_t("nomega_i", sigma%nomega_i)))
 end if

 ! =======================
 ! == Define variables ===
 ! =======================
 ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: &
   'sigma_nband', 'scr_nband', 'gwcalctyp', 'usepawu', "nfreqre", "nfreqim", "nfreqim_conv"])
 NCF_CHECK(ncerr)

 ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: &
   'ecutwfn', 'ecuteps', 'ecutsigx', 'omegasrdmax', 'deltae', 'omegasrmax', 'scissor_ene'])
 NCF_CHECK(ncerr)

 ! TODO: Decrease size of file: Remove arrays whose size scale as mband ** 2
 ! especially those that are not commonly used e.g. hhartree.
 ncerr = nctk_def_arrays(ncid, [ &
   nctkarr_t("kptgw", "dp", "number_of_reduced_dimensions, nkptgw"),&
   nctkarr_t("minbnd", "i", "nkptgw, number_of_spins"),&
   nctkarr_t("maxbnd", "i", "nkptgw, number_of_spins"), &
   nctkarr_t('degwgap', "dp", 'number_of_kpoints, number_of_spins'),&
   nctkarr_t('egwgap', "dp", 'number_of_kpoints, number_of_spins'),&
   nctkarr_t('en_qp_diago', "dp",'max_number_of_states, number_of_kpoints, number_of_spins'),&
   nctkarr_t('e0', "dp", 'max_number_of_states, number_of_kpoints, number_of_spins'),&
   nctkarr_t('e0gap', "dp", 'number_of_kpoints, number_of_spins'),&
   nctkarr_t('sigxme', "dp", 'nbgw, number_of_kpoints, ndim_sig'),&
   nctkarr_t('vxcme', "dp", 'nbgw, number_of_kpoints, ndim_sig'),&
   nctkarr_t('degw', "dp", 'cplex, nbgw, number_of_kpoints, number_of_spins'),&
   nctkarr_t('dsigmee0', "dp", 'cplex, nbgw, number_of_kpoints, ndim_sig'),&
   nctkarr_t('egw', "dp",'cplex, max_number_of_states, number_of_kpoints, number_of_spins'),&
   nctkarr_t('hhartree', "dp",'cplex, nbgw, nbgw, number_of_kpoints, ndim_sig'),&
   nctkarr_t('sigmee', "dp", 'cplex, nbgw, number_of_kpoints, ndim_sig'),&
   nctkarr_t('sigcmee0', "dp",'cplex, nbgw, number_of_kpoints, ndim_sig'),&
   nctkarr_t('sigcme4sd', "dp",'cplex, nbgw, number_of_kpoints, nomega4sd, ndim_sig'),&
   nctkarr_t('sigxcme4sd', "dp", 'cplex, nbgw, number_of_kpoints, nomega4sd, ndim_sig'),&
   nctkarr_t('ze0',"dp", 'cplex, nbgw, number_of_kpoints, number_of_spins'),&
   nctkarr_t('omega4sd', "dp", 'cplex, nbgw, number_of_kpoints, nomega4sd, number_of_spins') &
 ])
 NCF_CHECK(ncerr)

 if (sigma%needs_eigvec_qp) then
   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t('eigvec_qp', "dp",'cplex, max_number_of_states, max_number_of_states, number_of_kpoints, number_of_spins')])
   NCF_CHECK(ncerr)
 end if

 if (sigma%usepawu == 0) then
   ncerr = nctk_def_arrays(ncid, nctkarr_t("vUme", "dp", 'nbgw, number_of_kpoints, ndim_sig'))
   NCF_CHECK(ncerr)
 end if

 if (epsm1%nomega > 0) then
   ncerr = nctk_def_arrays(ncid, nctkarr_t('scr_omega', "dp", 'cplex, scr_nomega'))
   NCF_CHECK(ncerr)
   ABI_MALLOC(rdata2, (2, epsm1%nomega))
   rdata2 = c2r(epsm1%omega)
   NCF_CHECK(nf90_put_var(ncid, vid('scr_omega'), rdata2 * Ha_eV))
   ABI_FREE(rdata2)
  end if

 if (sigma%nomega_r > 0) then
   ncerr = nctk_def_arrays(ncid, [&
     nctkarr_t('omega_r', "dp", "nomega_r"),&
     nctkarr_t('sigcme', "dp", 'cplex, nbgw, number_of_kpoints, nomega_r, ndim_sig'),&
     nctkarr_t('sigxcme', "dp", 'cplex, nbgw, number_of_kpoints, nomega_r, ndim_sig')])
   NCF_CHECK(ncerr)
 end if

 if (sigma%nomega_i > 0) then
   ncerr = nctk_def_arrays(ncid, [&
     nctkarr_t('sigxcmesi', "dp", 'cplex, nbgw, number_of_kpoints, nomega_i, ndim_sig'),&
     nctkarr_t('sigcmesi', "dp",'cplex, nbgw, number_of_kpoints, nomega_i, ndim_sig'),&
     nctkarr_t('omega_i', "dp", 'cplex, nomega_i')])
   NCF_CHECK(ncerr)
 end if

 if (allocated(sigma%m_ks_to_qp)) then
   ncerr = nctk_def_arrays(ncid, [nctkarr_t('m_ks_to_qp', "dp", &
       "cplex, max_number_of_states, max_number_of_states, number_of_kpoints, number_of_spins")])
   NCF_CHECK(ncerr)
 end if

 ! =====================
 ! === Start writing ===
 ! =====================
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, vid('ecutwfn'), Sigp%ecutwfn))
 NCF_CHECK(nf90_put_var(ncid, vid('ecuteps'), Sigp%ecuteps))
 NCF_CHECK(nf90_put_var(ncid, vid('ecutsigx'), Sigp%ecutsigx))
 NCF_CHECK(nf90_put_var(ncid, vid('sigma_nband'), Sigp%nbnds))
 NCF_CHECK(nf90_put_var(ncid, vid('scr_nband'), epsm1%Hscr%nbnds_used))
 NCF_CHECK(nf90_put_var(ncid, vid('gwcalctyp'), sigma%gwcalctyp))
 NCF_CHECK(nf90_put_var(ncid, vid('usepawu'), sigma%usepawu))
 NCF_CHECK(nf90_put_var(ncid, vid('nfreqre'), epsm1%nomega_r))
 NCF_CHECK(nf90_put_var(ncid, vid('nfreqim'), epsm1%nomega_i))
 NCF_CHECK(nf90_put_var(ncid, vid('nfreqim_conv'), epsm1%nomega_i_conv))
 NCF_CHECK(nf90_put_var(ncid, vid('kptgw'), sigma%kptgw))
 NCF_CHECK(nf90_put_var(ncid, vid('minbnd'), sigma%minbnd))
 NCF_CHECK(nf90_put_var(ncid, vid('maxbnd'),sigma%maxbnd))
 NCF_CHECK(nf90_put_var(ncid, vid('omegasrdmax'), sigma%maxomega4sd * Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('deltae'), sigma%deltae * Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('omegasrmax'), sigma%maxomega_r * Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('scissor_ene'), sigma%scissor_ene * Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('degwgap'), sigma%degwgap * Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('egwgap'), sigma%egwgap * Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('en_qp_diago'), sigma%en_qp_diago * Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('e0'), sigma%e0 * Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('e0gap'), sigma%e0gap * Ha_eV))

 if (sigma%nomega_r > 0) then
   NCF_CHECK(nf90_put_var(ncid, vid('omega_r'), sigma%omega_r * Ha_eV))
 end if

 NCF_CHECK(nf90_put_var(ncid, vid('sigxme'), sigma%sigxme * Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('vxcme'), sigma%vxcme * Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('vUme'), sigma%vUme * Ha_eV))

 ! Have to transfer complex arrays
 ABI_MALLOC(rdata4,(cplex, b1gw:b2gw, sigma%nkibz, sigma%nsppol))
 rdata4=c2r(sigma%degw)
 NCF_CHECK(nf90_put_var(ncid, vid('degw'), rdata4 * Ha_eV))
 ABI_FREE(rdata4)

 ABI_MALLOC(rdata4,(cplex, b1gw:b2gw, sigma%nkibz, sigma%nsppol*sigma%nsig_ab))
 rdata4 = c2r(sigma%dsigmee0)
 NCF_CHECK(nf90_put_var(ncid, vid('dsigmee0'), rdata4))
 ABI_FREE(rdata4)

 ABI_MALLOC(rdata4, (cplex, sigma%nbnds, sigma%nkibz, sigma%nsppol))
 rdata4 = c2r(sigma%egw)
 NCF_CHECK(nf90_put_var(ncid, vid('egw'), rdata4 *Ha_eV))
 ABI_FREE(rdata4)

 if (sigma%needs_eigvec_qp) then
   ABI_MALLOC(rdata5, (cplex, sigma%nbnds, sigma%nbnds, sigma%nkibz, sigma%nsppol))
   rdata5 = c2r(sigma%eigvec_qp)
   NCF_CHECK(nf90_put_var(ncid, vid('eigvec_qp'), rdata5))
   ABI_FREE(rdata5)
 end if

 ABI_MALLOC(rdata5,(cplex, nbgw, nbgw, sigma%nkibz, sigma%nsppol * sigma%nsig_ab))
 rdata5 = c2r(sigma%hhartree)
 NCF_CHECK(nf90_put_var(ncid, vid('hhartree'), rdata5 * Ha_eV))
 ABI_FREE(rdata5)

 if (sigma%nomega_r > 0) then
   ABI_MALLOC(rdata5,(cplex, nbgw, sigma%nkibz, sigma%nomega_r, sigma%nsppol*sigma%nsig_ab))
   rdata5 = c2r(sigma%sigcme)
   NCF_CHECK(nf90_put_var(ncid, vid('sigcme'), rdata5 * Ha_eV))
   ABI_FREE(rdata5)
 end if

 ABI_MALLOC(rdata4, (cplex, nbgw, sigma%nkibz, sigma%nsppol*sigma%nsig_ab))
 rdata4 = c2r(sigma%sigmee)
 NCF_CHECK(nf90_put_var(ncid, vid('sigmee'), rdata4 * Ha_eV))
 ABI_FREE(rdata4)

 ABI_MALLOC(rdata4, (cplex, nbgw, sigma%nkibz, sigma%nsppol*sigma%nsig_ab))
 rdata4 = c2r(sigma%sigcmee0)
 NCF_CHECK(nf90_put_var(ncid, vid('sigcmee0'), rdata4 * Ha_eV))
 ABI_FREE(rdata4)

 if (sigma%nomega_i > 0) then
  ABI_MALLOC(rdata5, (cplex, nbgw, sigma%nkibz, sigma%nomega_i, sigma%nsppol*sigma%nsig_ab))
  rdata5 = c2r(sigma%sigcmesi)
  NCF_CHECK(nf90_put_var(ncid, vid('sigcmesi'), rdata5*Ha_eV))
  ABI_FREE(rdata5)
 end if

 ABI_MALLOC(rdata5, (cplex, nbgw, sigma%nkibz, sigma%nomega4sd, sigma%nsppol*sigma%nsig_ab))
 rdata5 = c2r(sigma%sigcme4sd)
 NCF_CHECK(nf90_put_var(ncid, vid('sigcme4sd'), rdata5 * Ha_eV))
 ABI_FREE(rdata5)

 if (sigma%nomega_r > 0) then
   ABI_MALLOC(rdata5,(cplex, nbgw, sigma%nkibz, sigma%nomega_r, sigma%nsppol*sigma%nsig_ab))
   rdata5 = c2r(sigma%sigxcme)
   NCF_CHECK(nf90_put_var(ncid, vid('sigxcme'), rdata5 * Ha_eV))
   ABI_FREE(rdata5)
 end if

 if (sigma%nomega_i > 0) then
   ABI_MALLOC(rdata5,(cplex, nbgw, sigma%nkibz, sigma%nomega_i, sigma%nsppol*sigma%nsig_ab))
   rdata5 = c2r(sigma%sigxcmesi)
   NCF_CHECK(nf90_put_var(ncid, vid('sigxcmesi'), rdata5 * Ha_eV))
   ABI_FREE(rdata5)
 end if

 if (allocated(sigma%m_ks_to_qp)) then
   ABI_MALLOC(rdata5,(cplex, sigma%nbnds, sigma%nbnds, sigma%nkibz, sigma%nsppol))
   rdata5 = c2r(sigma%m_ks_to_qp)
   NCF_CHECK(nf90_put_var(ncid, vid('m_ks_to_qp'), rdata5))
   ABI_FREE(rdata5)
 end if

 ABI_MALLOC(rdata5, (cplex, nbgw, sigma%nkibz, sigma%nomega4sd, sigma%nsppol*sigma%nsig_ab))
 rdata5 = c2r(sigma%sigxcme4sd)
 NCF_CHECK(nf90_put_var(ncid, vid('sigxcme4sd'), rdata5 * Ha_eV))
 ABI_FREE(rdata5)

 ABI_MALLOC(rdata4, (cplex, nbgw, sigma%nkibz, sigma%nsppol))
 rdata4 = c2r(sigma%ze0)
 NCF_CHECK(nf90_put_var(ncid, vid('ze0'), rdata4))
 ABI_FREE(rdata4)

 if (sigma%nomega_i > 0) then
   ABI_MALLOC(rdata2, (cplex, sigma%nomega_i))
   rdata2 = c2r(sigma%omega_i)
   NCF_CHECK(nf90_put_var(ncid, vid('omega_i'), rdata2 * Ha_eV))
   ABI_FREE(rdata2)
 end if

 ABI_MALLOC(rdata5, (cplex, nbgw, sigma%nkibz, sigma%nomega4sd, sigma%nsppol))
 rdata5 = c2r(sigma%omega4sd)
 NCF_CHECK(nf90_put_var(ncid, vid('omega4sd'), rdata5 * Ha_eV))
 ABI_FREE(rdata5)


contains
 integer function vid(vname)
   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
 end function vid

end function sigma_ncwrite
!!***

!----------------------------------------------------------------------

!!****f* m_sigma/sigma_distribute_bks
!! NAME
!!  sigma_distribute_bks
!!
!! FUNCTION
!!  Distribute the loop over (b,k,s) used to calculate the self-energy matrix elements
!!  taking into account the MPI distribution of the wavefunctions and the use of
!!  symmetries to reduce the BZ sum to an appropriate irreducible wedge.
!!
!! INPUTS
!! nsppol
!! can_symmetrize(nsppol)=.TRUE if symmetries can be used to reduce the number of k-points to be summed.
!! Kmesh<kmesh_t>
!! Qmesh<kmesh_t>
!! Ltg_kgw<littlegroup_t>
!! Wfd(wfdgw_t)
!! mg0(3)
!! kptgw(3)
!! [bks_mask(Wfd%mband,Kmesh%nbz,nsppol)]
!! [got(Wfd%nproc)]=The number of tasks already assigned to the nodes.
!! [global]=If true, an MPI global communication is performed such that each node will have the same table. Useful
!!   if for implementing algorithms in which each node needs to know the global distribution of the tasks, not only
!!   the task it has to complete. Defaults to .FALSE.
!!
!! OUTPUT
!!  my_nbks
!!  proc_distrb(Wfd%mband,Kmesh%nbz,nsppol)
!!
!! SIDE EFFECTS
!!  Wfd%bks_tab
!!
!! SOURCE

subroutine sigma_distribute_bks(Wfd,Kmesh,Ltg_kgw,Qmesh,nsppol,can_symmetrize,kptgw,mg0,my_nbks,proc_distrb,got,bks_mask,global)

!Arguments ------------------------------------
!scalars
 class(wfdgw_t),intent(inout) :: Wfd
 integer,intent(in) :: nsppol
 integer,intent(out) :: my_nbks
 logical,optional,intent(in) :: global
 type(kmesh_t),intent(in) :: Kmesh,Qmesh
 type(littlegroup_t),intent(in) :: Ltg_kgw
!arrays
 integer,intent(in) :: mg0(3)
 integer,optional,intent(inout) :: got(Wfd%nproc)
 integer,intent(out) :: proc_distrb(Wfd%mband,Kmesh%nbz,nsppol)
 real(dp),intent(in) :: kptgw(3)
 logical,intent(in) :: can_symmetrize(Wfd%nsppol)
 logical,optional,intent(in) :: bks_mask(Wfd%mband,Kmesh%nbz,nsppol)

!Local variables-------------------------------
!scalars
 integer :: ierr,ik_bz,ik_ibz,spin,iq_bz,my_nband
 !character(len=500) :: msg
!arrays
 integer :: g0(3)
 real(dp) :: kgwmk(3)
 integer :: get_more(Wfd%nproc),my_band_list(Wfd%mband)
 logical :: bmask(Wfd%mband)
!************************************************************************

 call wfd%update_bkstab()

 get_more=0; if (PRESENT(got)) get_more=got

 ! Different distribution of tasks depending whether symmetries can be used or not.
 proc_distrb= xmpi_undefined_rank

 do spin=1,Wfd%nsppol

   if (can_symmetrize(spin)) then
     do ik_bz=1,Kmesh%nbz
       ik_ibz = Kmesh%tab(ik_bz)
       kgwmk = kptgw-Kmesh%bz(:,ik_bz) ! kptgw must be inside the BZ
       call findqg0(iq_bz,g0,kgwmk,Qmesh%nbz,Qmesh%bz,mG0) ! <- (mg0=mG0) Identify q_bz and G0 where q_bz+G0=k_gw-k_bz
       if (Ltg_kgw%ibzq(iq_bz)==1) then
         bmask=.FALSE.; bmask(1:Wfd%nband(ik_ibz,spin))=.TRUE.
         if (PRESENT(bks_mask)) bmask = bks_mask(:,ik_bz,spin)
         call wfd%distribute_bands(ik_ibz,spin,my_nband,my_band_list,got=get_more,bmask=bmask)
         if (my_nband>0) proc_distrb(my_band_list(1:my_nband),ik_bz,spin)=Wfd%my_rank
       end if
     end do

   else
     ! No symmetries for this spin. Divide the full BZ among procs.
     do ik_bz=1,Kmesh%nbz
       ik_ibz = Kmesh%tab(ik_bz)
       bmask=.FALSE.; bmask(1:Wfd%nband(ik_ibz,spin))=.TRUE.
       if (PRESENT(bks_mask)) bmask = bks_mask(:,ik_bz,spin)
       call wfd%distribute_bands(ik_ibz,spin,my_nband,my_band_list,got=get_more,bmask=bmask)
       if (my_nband>0) proc_distrb(my_band_list(1:my_nband),ik_bz,spin)=Wfd%my_rank
     end do
   end if
 end do ! spin

 if (PRESENT(global)) then
   if (global) then ! Each node will have the same table so that it will know how the tasks are distributed.
     proc_distrb = proc_distrb + 1
     where (proc_distrb == xmpi_undefined_rank + 1)
       proc_distrb = 0
     end where
     call xmpi_sum(proc_distrb,Wfd%comm,ierr)
     where (proc_distrb == 0)
       proc_distrb = xmpi_undefined_rank
     elsewhere
       proc_distrb = proc_distrb - 1
     end where
     !where (proc_distrb /= xmpi_undefined_rank)
     !  ltest = (ANY(proc_distrb == (/(ii,ii=0,Wfd%nproc-1)/)))
     !end where
     !if (.not.ltest) then
     !  write(std_out,*)proc_distrb
     !  ABI_BUG("Bug in the generation of proc_distrb table")
     !end if
   end if
 end if

 my_nbks = COUNT(proc_distrb==Wfd%my_rank)
 if (PRESENT(got)) got=get_more

end subroutine sigma_distribute_bks
!!***

!----------------------------------------------------------------------

end module m_sigma
!!***
