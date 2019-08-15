!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_sigma
!! NAME
!!  m_sigma
!!
!! FUNCTION
!!  This module provides the definition of the sigma_t data type
!!  used to store results of the GW calculation as well as as
!!  methods bound to the object.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2019 ABINIT group (MG, FB, GMR, VO, LR, RWG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_sigma

 use defs_basis
 use defs_datatypes
 use m_xmpi
 use m_abicore
 use m_errors
 use iso_c_binding
 use m_nctk
 use m_yaml
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_wfd

 use defs_abitypes,    only : MPI_type
 use m_numeric_tools,  only : c2r
 use m_gwdefs,         only : unt_gw, unt_sig, unt_sgr, unt_sgm, unt_gwdiag, sigparams_t, sigma_needs_w
 use m_crystal,        only : crystal_t
 use m_bz_mesh,        only : kmesh_t, littlegroup_t, findqg0
 use m_screening,      only : epsilonm1_results
 use m_stream_string,  only : stream_string

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_sigma/sigma_t
!! NAME
!! sigma_t
!!
!! FUNCTION
!! For the GW part of ABINIT, the sigma_t structured datatype
!! gather the results of a sigma calculation.
!!
!! TODO
!!   ragged arrays (nk,nsppol) --> values ?
!!
!! SOURCE

 type,public :: sigma_t

  integer :: b1gw,b2gw      ! min and Max gw band indeces over spin and k-points (used to dimension arrays)
  integer :: gwcalctyp      ! Flag defining the calculation type.
  integer :: nkptgw         ! No. of points calculated
  integer :: nkibz          ! No. of irreducible k-points.
  integer :: nbnds          ! Total number of bands
  integer :: nomega_r       ! No. of real frequencies for the spectral function.
  integer :: nomega_i       ! No. of frequencies along the imaginary axis.
  integer :: nomega4sd      ! No. of real frequencies to evaluate the derivative of $\Sigma(E)$.
  integer :: nsig_ab        ! 1 if nspinor=1,4 for noncollinear case.
  integer :: nsppol         ! No. of spin polarizations.
  integer :: usepawu        ! 1 if we are using LDA+U as starting point (only for PAW)

  real(dp) :: deltae       ! Frequency step for the calculation of d\Sigma/dE
  real(dp) :: maxomega4sd  ! Max frequency around E_ks for d\Sigma/dE.
  real(dp) :: maxomega_r   ! Max frequency for spectral function.
  real(dp) :: scissor_ene  ! Scissor energy value. zero for None.

  integer,allocatable :: maxbnd(:,:)
  ! maxbnd(nkptgw,nsppol)
  ! Max band index considered in GW for this k-point.

  integer,allocatable :: minbnd(:,:)
  ! minbnd(nkptgw,nsppol)
  ! Min band index considered in GW for this k-point.

  !real(dp),allocatable :: ame(:,:,:)
  ! ame(nbnds,nkibz,nomega))
  ! Diagonal matrix elements of the spectral function.
  ! Commented out, it can be calculated from the other quantities

  real(dp),allocatable :: degwgap(:,:)
  ! degwgap(nkibz,nsppol)
  ! Difference btw the QP and the KS optical gap.

  real(dp),allocatable :: egwgap(:,:)
  ! egwgap(nkibz,nsppol))
  ! QP optical gap at each k-point and spin.

  real(dp),allocatable :: en_qp_diago(:,:,:)
  ! en_qp_diago(nbnds,nkibz,nsppol))
  ! QP energies obtained from the diagonalization of the Hermitian approximation to Sigma (QPSCGW)

  real(dp),allocatable :: e0(:,:,:)
  ! e0(nbnds,nkibz,nsppol)
  ! KS eigenvalues for each band, k-point and spin. In case of self-consistent?

  real(dp),allocatable :: e0gap(:,:)
  ! e0gap(nkibz,nsppol),
  ! KS gap at each k-point, for each spin.

  real(dp),allocatable :: omega_r(:)
  ! omega_r(nomega_r)
  ! real frequencies used for the self energy.

  real(dp),allocatable :: kptgw(:,:)
  ! kptgw(3,nkptgw)
  ! ! TODO there is a similar array in sigparams_t
  ! List of calculated k-points.

  real(dp),allocatable :: sigxme(:,:,:)
  ! sigxme(b1gw:b2gw,nkibz,nsppol*nsig_ab))
  ! Diagonal matrix elements $\<nks|\Sigma_x|nks\>$

  complex(dp),allocatable :: x_mat(:,:,:,:)
  ! x_mat(b1gw:b2gw,b1gw:b2gw,nkibz,nsppol*nsig_ab))
  ! Matrix elements of $\<nks|\Sigma_x|nk's\>$

  real(dp),allocatable :: vxcme(:,:,:)
  ! vxcme(b1gw:b2gw,nkibz,nsppol*nsig_ab))
  ! $\<nks|v_{xc}[n_val]|nks\>$ matrix elements of vxc (valence-only contribution).

  real(dp),allocatable :: vUme(:,:,:)
  ! vUme(b1gw:b2gw,nkibz,nsppol*nsig_ab))
  ! $\<nks|v_{U}|nks\>$ for LDA+U.

  complex(dpc),allocatable :: degw(:,:,:)
  ! degw(b1gw:b2gw,nkibz,nsppol))
  ! Difference between the QP and the KS energies.

  complex(dpc),allocatable :: dsigmee0(:,:,:)
  ! dsigmee0(b1gw:b2gw,nkibz,nsppol*nsig_ab))
  ! Derivative of $\Sigma_c(E)$ calculated at the KS eigenvalue.

  complex(dpc),allocatable :: egw(:,:,:)
  ! degw(nbnds,nkibz,nsppol))
  ! QP energies, $\epsilon_{nks}^{QP}$.

  complex(dpc),allocatable :: eigvec_qp(:,:,:,:)
  ! eigvec_qp(nbnds,nbnds,nkibz,nsppol))
  ! Expansion of the QP amplitudes in the QP basis set of the previous iteration.

  complex(dpc),allocatable :: m_lda_to_qp(:,:,:,:)
  ! (%nbnds,%nbnds,%nibz,%nsppol))
  !  m_lda_to_qp(ib,jb,k,s) := <\psi_{ib,k,s}^{KS}|\psi_{jb,k,s}^{QP}>

  complex(dpc),allocatable :: hhartree(:,:,:,:)
  ! hhartree(b1gw:b2gw,b1gw:b2gw,nkibz,nsppol*nsig_ab)
  ! $\<nks|T+v_H+v_{loc}+v_{nl}|mks\>$

  complex(dpc),allocatable :: sigcme(:,:,:,:)
  ! sigcme(b1gw:b2gw,nkibz,nomega_r,nsppol*nsig_ab))
  ! $\<nks|\Sigma_{c}(E)|nks\>$ at each nomega_r frequency

  complex(dpc),allocatable :: sigmee(:,:,:)
  ! sigmee(b1gw:b2gw,nkibz,nsppol*nsig_ab))
  ! $\Sigma_{xc}E_{KS} + (E_{QP}- E_{KS})*dSigma/dE_KS

  complex(dpc),allocatable :: sigcmee0(:,:,:)
  ! sigcmee0(b1gw:b2gw,nkibz,nsppol*nsig_ab))
  ! Diagonal mat. elements of $\Sigma_c(E)$ calculated at the KS energy $E_{KS}$

  complex(dpc),allocatable :: sigcmesi(:,:,:,:)
  ! sigcmesi(b1gw:b2gw,nkibz,nomega_i,nsppol*nsig_ab))
  ! Matrix elements of $\Sigma_c$ along the imaginary axis.
  ! Only used in case of analytical continuation.

  complex(dpc),allocatable :: sigcme4sd(:,:,:,:)
  ! sigcme4sd(b1gw:b2gw,nkibz,nomega4sd,nsppol*nsig_ab))
  ! Diagonal matrix elements of \Sigma_c around the zeroth order eigenvalue (usually KS).

  complex(dpc),allocatable :: sigxcme(:,:,:,:)
  ! sigxcme(b1gw:b2gw,nkibz,nomega_r,nsppol*nsig_ab))
  ! $\<nks|\Sigma_{xc}(E)|nks\>$ at each real frequency frequency.

  complex(dpc),allocatable :: sigxcmesi(:,:,:,:)
  ! sigxcmesi(b1gw:b2gw,nkibz,nomega_i,nsppol*nsig_ab))
  ! Matrix elements of $\Sigma_{xc}$ along the imaginary axis.
  ! Only used in case of analytical continuation.

  complex(dpc),allocatable :: sigxcme4sd(:,:,:,:)
  ! sigxcme4sd(b1gw:b2gw,nkibz,nomega4sd,nsppol*nsig_ab))
  ! Diagonal matrix elements of \Sigma_xc for frequencies around the zeroth order eigenvalues.

  complex(dpc),allocatable :: ze0(:,:,:)
  ! ze0(b1gw:b2gw,nkibz,nsppol))
  ! renormalization factor. $(1-\dfrac{\partial\Sigma_c} {\partial E_{KS}})^{-1}$

  complex(dpc),allocatable :: omega_i(:)
  ! omegasi(nomega_i)
  ! Frequencies along the imaginary axis used for the analytical continuation.

  complex(dpc),allocatable :: omega4sd(:,:,:,:)
  ! omega4sd(b1gw:b2gw,nkibz,nomega4sd,nsppol).
  ! Frequencies used to evaluate the Derivative of Sigma.

 end type sigma_t

 public  :: sigma_init                  ! Initialize the object
 public  :: sigma_free                  ! Deallocate memory
 public  :: sigma_get_exene             ! Compute exchange energy.
 public  :: sigma_ncwrite               ! Write data in netcdf format.
 public  :: write_sigma_header
 public  :: write_sigma_results
 public  :: print_Sigma_perturbative
 public  :: print_Sigma_QPSC
 public  :: sigma_distribute_bks
!!***


CONTAINS  !========================================================================================
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
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      findqg0,wfd_distribute_bands,wfd_update_bkstab,xmpi_sum
!!
!! SOURCE

subroutine write_sigma_header(Sigp,Er,Cryst,Kmesh,Qmesh)

!Arguments ------------------------------------
!scalars
 type(kmesh_t),intent(in) :: Kmesh,Qmesh
 type(crystal_t),intent(in) :: Cryst
 type(Epsilonm1_results),intent(in) :: Er
 type(sigparams_t),intent(in) :: Sigp

!Local variables-------------------------------
!scalars
 integer :: gwcalctyp,mod10
 character(len=500) :: msg

! *************************************************************************

 call wrtout([std_out, ab_out], ' SIGMA fundamental parameters:')

 gwcalctyp=Sigp%gwcalctyp
 mod10=MOD(Sigp%gwcalctyp,10)

 SELECT CASE (mod10)
 CASE (0)
   write(msg,'(a,i2)')' PLASMON POLE MODEL ',Sigp%ppmodel
 CASE (1)
   write(msg,'(a)')' ANALYTIC CONTINUATION'
 CASE (2)
   write(msg,'(a)')' CONTOUR DEFORMATION'
 CASE (5)
   write(msg,'(a)')' Hartree-Fock'
 CASE (6)
   write(msg,'(a)')' Screened Exchange'
 CASE (7)
   write(msg,'(a)')' COHSEX'
 CASE (8)
   write(msg,'(a,i2)')' MODEL GW with PLASMON POLE MODEL ',Sigp%ppmodel
 CASE (9)
   write(msg,'(a)')' MODEL GW without PLASMON POLE MODEL'
 CASE DEFAULT
   write(msg,'(a,i3)')' Wrong value for Sigp%gwcalctyp = ',Sigp%gwcalctyp
   MSG_BUG(msg)
 END SELECT
 call wrtout([std_out, ab_out], msg)

 write(msg,'(a,i12)')' number of plane-waves for SigmaX         ',Sigp%npwx
 call wrtout([std_out, ab_out], msg)
 write(msg,'(a,i12)')' number of plane-waves for SigmaC and W   ',Sigp%npwc
 call wrtout([std_out, ab_out], msg)
 write(msg,'(a,i12)')' number of plane-waves for wavefunctions  ',Sigp%npwwfn
 call wrtout([std_out, ab_out], msg)
 write(msg,'(a,i12)')' number of bands                          ',Sigp%nbnds
 call wrtout([std_out, ab_out], msg)
 write(msg,'(a,i12)')' number of independent spin polarizations ',Sigp%nsppol
 call wrtout([std_out, ab_out], msg)
 write(msg,'(a,i12)')' number of spinorial components           ',Sigp%nspinor
 call wrtout([std_out, ab_out], msg)
 write(msg,'(a,i12)')' number of k-points in IBZ                ',Kmesh%nibz
 call wrtout([std_out, ab_out], msg)
 write(msg,'(a,i12)')' number of q-points in IBZ                ',Qmesh%nibz
 call wrtout([std_out, ab_out], msg)
 write(msg,'(a,i12)')' number of symmetry operations            ',Cryst%nsym
 call wrtout([std_out, ab_out], msg)
 write(msg,'(a,i12)')' number of k-points in BZ                 ',Kmesh%nbz
 call wrtout([std_out, ab_out], msg)
 write(msg,'(a,i12)')' number of q-points in BZ                 ',Qmesh%nbz
 call wrtout([std_out, ab_out], msg)
 write(msg,'(a,i12)')' number of frequencies for dSigma/dE      ',Sigp%nomegasrd
 call wrtout([std_out, ab_out], msg)
 write(msg,'(a,f12.2)')' frequency step for dSigma/dE [eV]        ',Sigp%deltae*Ha_eV
 call wrtout([std_out, ab_out], msg)
 write(msg,'(a,i12)')' number of omega for Sigma on real axis   ',Sigp%nomegasr
 call wrtout([std_out, ab_out], msg)
 write(msg,'(a,f12.2)')' max omega for Sigma on real axis  [eV]   ',Sigp%maxomega_r*Ha_eV
 call wrtout([std_out, ab_out], msg)
 write(msg,'(a,f12.2)')' zcut for avoiding poles [eV]             ',Sigp%zcut*Ha_eV
 call wrtout([std_out, ab_out], msg)

 if (Sigp%mbpt_sciss>0.1d-4) then
   write(msg,'(a,f12.2)')' scissor energy [eV]                      ',Sigp%mbpt_sciss*Ha_eV
   call wrtout([std_out, ab_out], msg)
 end if

 if (mod10==1) then
   write(msg,'(a,i12)')' number of imaginary frequencies for Sigma',Sigp%nomegasi
   call wrtout([std_out, ab_out], msg)
   write(msg,'(a,f12.2)')' max omega for Sigma on imag axis  [eV]   ',Sigp%omegasimax*Ha_eV
   call wrtout([std_out, ab_out], msg)
 end if

 if (sigma_needs_w(Sigp)) then
   write(msg,'(2a)')ch10,' EPSILON^-1 parameters (SCR file):'
   call wrtout([std_out, ab_out], msg)
   write(msg,'(a,i12)')' dimension of the eps^-1 matrix on file   ',Er%Hscr%npwe
   call wrtout([std_out, ab_out], msg)
   write(msg,'(a,i12)')' dimension of the eps^-1 matrix used      ',Er%npwe
   call wrtout([std_out, ab_out], msg)
   write(msg,'(a,i12)')' number of plane-waves for wavefunctions  ',Er%Hscr%npwwfn_used
   call wrtout([std_out, ab_out], msg)
   write(msg,'(a,i12)')' number of bands                          ',Er%Hscr%nbnds_used
   call wrtout([std_out, ab_out], msg)
   write(msg,'(a,i12)')' number of q-points in IBZ                ',Qmesh%nibz
   call wrtout([std_out, ab_out], msg)
   write(msg,'(a,i12)')' number of frequencies                    ',Er%nomega
   call wrtout([std_out, ab_out], msg)
   write(msg,'(a,i12)')' number of real frequencies               ',Er%nomega_r
   call wrtout([std_out, ab_out], msg)
   write(msg,'(a,i12)')' number of imag frequencies               ',Er%nomega_i
   call wrtout([std_out, ab_out], msg)
 end if

 write(msg,'(3a)')ch10,' matrix elements of self-energy operator (all in [eV])',ch10
 call wrtout([std_out, ab_out], msg)

 if (gwcalctyp<10) then
   write(msg,'(a)')' Perturbative Calculation'
 else if (gwcalctyp<20) then
   write(msg,'(a)')' Self-Consistent on Energies only'
 else
   write(msg,'(a)')' Self-Consistent on Energies and Wavefunctions'
 end if
 call wrtout([std_out, ab_out], msg)

end subroutine write_sigma_header
!!***

!----------------------------------------------------------------------

!!****f* m_sigma/write_sigma_results
!! NAME
!! write_sigma_results
!!
!! FUNCTION
!!  Write the final results of the GW calculation.
!!
!! INPUTS
!!  KS_BSt<ebands_t>=Info on the KS band structure energies.
!!     %eig(mband,nkibz,nsppol)= KS energies
!!  ikibz= index of the k-point in the array kibz, where GW corrections are calculated
!!  ikcalc= index of the k-point in the array Sigp%kptgw2bz
!!  Sigp=sigparams_t datatype
!!  sr=sigma results datatype
!!
!! OUTPUT
!!  (for writing routines, no output) otherwise, should be described
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      findqg0,wfd_distribute_bands,wfd_update_bkstab,xmpi_sum
!!
!! SOURCE
!!

subroutine write_sigma_results(ikcalc,ikibz,Sigp,Sr,KS_BSt,use_yaml)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikcalc,ikibz, use_yaml
 type(ebands_t),intent(in) :: KS_BSt
 type(sigparams_t),intent(in) :: Sigp
 type(sigma_t),intent(in) :: Sr

!Local variables-------------------------------
!scalars
 integer :: ib,io,is
 integer :: gwcalctyp,mod10
 character(len=500) :: msg
 type(yamldoc_t) :: ydoc
!arrays
 character(len=12) :: tag_spin(2)

! *************************************************************************

 gwcalctyp=Sigp%gwcalctyp
 mod10=MOD(Sigp%gwcalctyp,10)

 !unt_gw  File with GW corrections.
 !unt_sig Self-energy as a function of frequency.
 !unt_sgr Derivative wrt omega of the Self-energy.
 !unt_sgm Sigma on the Matsubara axis.

 tag_spin=(/'            ','            '/); if (Sr%nsppol==2) tag_spin=(/',  SPIN UP  ',',  SPIN DOWN'/)

 do is=1,Sr%nsppol
   write(msg,'(2a,3f8.3,a)')ch10,' k = ',Sigp%kptgw(:,ikcalc),tag_spin(is)
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')

   msg = ' Band     E0 <VxcLDA>   SigX SigC(E0)      Z dSigC/dE  Sig(E)    E-E0       E'
   if (Sr%usepawu/=0) then
     msg = ' Band     E0 <VxcLDA>   <H_U>  SigX SigC(E0)      Z dSigC/dE  Sig(E)    E-E0       E'
   end if

   if (gwcalctyp>=10) then
     write(msg,'(2a)')&
     ' Band     E_lda   <Vxclda>   E(N-1)  <Hhartree>   SigX  SigC[E(N-1)]',&
     '    Z     dSigC/dE  Sig[E(N)]  DeltaE  E(N)_pert E(N)_diago'
   end if
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')

   ydoc = yamldoc_open('SelfEnergy_ee', "", width=11, real_fmt='(3f8.3)')
   ydoc%use_yaml = use_yaml
   call ydoc%add_real1d('kpoint', Sigp%kptgw(:,ikcalc))
   call ydoc%add_int('spin', is, int_fmt="(i1)")
   call ydoc%add_real('KS_gap', Sr%e0gap(ikibz,is)*Ha_eV)
   call ydoc%add_real('QP_gap', Sr%egwgap(ikibz,is)*Ha_eV)
   call ydoc%add_real('Delta_QP_KS', Sr%degwgap(ikibz,is)*Ha_eV)
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

   do ib=Sigp%minbnd(ikcalc,is),Sigp%maxbnd(ikcalc,is)
     if (gwcalctyp>=10) then
       call print_Sigma_QPSC(Sr,ikibz,ib,is,KS_BSt,unit=ab_out, ydoc=ydoc)
       call print_Sigma_QPSC(Sr,ikibz,ib,is,KS_BSt,unit=std_out,prtvol=1)

       write(unt_gwdiag,'(i6,3f9.4)')                                 &
        ib,                                                           &
        Sr%en_qp_diago(ib,ikibz,is)*Ha_eV,                            &
        (Sr%en_qp_diago(ib,ikibz,is) - KS_BSt%eig(ib,ikibz,is))*Ha_eV,&
        zero

     else
       ! If not ppmodel, write out also the imaginary part in ab_out
       SELECT CASE(mod10)
       CASE(1,2)
         call print_Sigma_perturbative(Sr,ikibz,ib,is,unit=ab_out,ydoc=ydoc,prtvol=1)
       CASE DEFAULT
         call print_Sigma_perturbative(Sr,ikibz,ib,is,unit=ab_out,ydoc=ydoc)
       END SELECT
       call print_Sigma_perturbative(Sr,ikibz,ib,is,unit=std_out,prtvol=1)
     end if

     write(unt_gw,'(i6,3f9.4)')          &
      ib,                               &
      REAL (Sr%egw (ib,ikibz,is))*Ha_eV,&
      REAL (Sr%degw(ib,ikibz,is))*Ha_eV,&
      AIMAG(Sr%egw (ib,ikibz,is))*Ha_eV
   end do !ib

   if (Sr%e0gap(ikibz,is)**2+Sr%egwgap(ikibz,is)**2+Sr%degwgap(ikibz,is)**2 > tol10) then
     ! Output the direct gap for each spin
     ! If all the gaps are zero, this means that it could not be computed in the calling routine
     write(msg,'(2a,f8.3)')ch10,' E^0_gap       ',Sr%e0gap(ikibz,is)*Ha_eV
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
     write(msg,'(a,f8.3)')      ' E^GW_gap      ',Sr%egwgap(ikibz,is)*Ha_eV
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
     write(msg,'(a,f8.3,a)')    ' DeltaE^GW_gap ',Sr%degwgap(ikibz,is)*Ha_eV,ch10
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
   end if

   call ydoc%write_and_free(ab_out)

   ! Output of the spectral function
   do io=1,Sr%nomega_r
     write(unt_sig,'(100(e12.5,2x))')&
      REAL(Sr%omega_r(io))*Ha_eV,&
      (REAL(Sr%sigxcme(ib,ikibz,io,is))*Ha_eV,&
      AIMAG(Sr%sigxcme(ib,ikibz,io,is))*Ha_eV,&
      gw_spectral_function(Sr,io,ib,ikibz,is),&
      ib=Sigp%minbnd(ikcalc,is),Sigp%maxbnd(ikcalc,is))
   end do
   !
   do ib=Sigp%minbnd(ikcalc,is),Sigp%maxbnd(ikcalc,is)
     write(unt_sgr,'("# ik, ib",2i5)')ikibz,ib
     do io=1,Sr%nomega4sd
       write(unt_sgr,'(100(e12.5,2x))')              &
         REAL (Sr%omega4sd  (ib,ikibz,io,is)) *Ha_eV,&
         REAL (Sr%sigxcme4sd(ib,ikibz,io,is)) *Ha_eV,&
         AIMAG(Sr%sigxcme4sd(ib,ikibz,io,is)) *Ha_eV
     end do
   end do
   !
   if (mod10==1) then ! For AC, write sigma matrix elements along the imaginary axis
     do ib=Sigp%minbnd(ikcalc,is),Sigp%maxbnd(ikcalc,is)
       write(unt_sgm,'("# ik, ib",2i5)')ikibz,ib
       do io=1,Sr%nomega_i
         write(unt_sgm,'(3(e12.5,2x))')             &
          AIMAG(Sr%omega_i(io))              *Ha_eV,&
          REAL (Sr%sigxcmesi(ib,ikibz,io,is))*Ha_eV,&
          AIMAG(Sr%sigxcmesi(ib,ikibz,io,is))*Ha_eV
       end do
     end do
   end if

 end do !is

end subroutine write_sigma_results
!!***

!----------------------------------------------------------------------

!!****f* m_sigma/gw_spectral_function
!! NAME
!! gw_spectral_function
!!
!! FUNCTION
!!  Compute the spectral function in the GW approximation
!!
!! INPUTS
!!  io,ib,ikibz,is=Frequency, band, k-point, spin index
!!  Sr=sigma results datatype
!!
!! OUTPUT
!!
!! PARENTS
!!
!! SOURCE

function gw_spectral_function(Sr,io,ib,ikibz,is) result(aw)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: io,ib,ikibz,is
 real(dp) :: aw
 type(sigma_t),intent(in) :: Sr

! *********************************************************************

 aw = one/pi*ABS(AIMAG(Sr%sigcme(ib,ikibz,io,is)))&
   /( (REAL(Sr%omega_r(io)-Sr%hhartree(ib,ib,ikibz,is)-Sr%sigxcme(ib,ikibz,io,is)))**2&
     +(AIMAG(Sr%sigcme(ib,ikibz,io,is)))**2) /Ha_eV

end function gw_spectral_function
!!***

!----------------------------------------------------------------------

!!****f* m_sigma/print_Sigma_perturbative
!! NAME
!! print_Sigma_perturbative
!!
!! FUNCTION
!!  Write the results of the GW calculation done with the perturbative approach
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_sigma
!!
!! CHILDREN
!!      findqg0,wfd_distribute_bands,wfd_update_bkstab,xmpi_sum
!!
!! SOURCE

subroutine print_Sigma_perturbative(Sr,ik_ibz,iband,isp,unit,prtvol,mode_paral,witheader,ydoc)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iband,ik_ibz,isp
 integer,optional,intent(in) :: prtvol,unit
 character(len=*),optional,intent(in) :: mode_paral
 logical,optional,intent(in) :: witheader
 type(sigma_t),intent(in) :: Sr
 type(yamldoc_t),intent(inout),optional :: ydoc

!Local variables-------------------------------
!scalars
 integer :: my_unt,verbose
 character(len=4) :: my_mode
 character(len=500) :: msg

! *********************************************************************

 my_unt =std_out; if (PRESENT(unit      )) my_unt =unit
 verbose=0      ; if (PRESENT(prtvol    )) verbose=prtvol
 my_mode='COLL' ; if (PRESENT(mode_paral)) my_mode=mode_paral

 if (PRESENT(witheader)) then
   if (witheader) then
     call wrtout(my_unt,' Band     E0 <VxcLDA>   SigX SigC(E0)      Z dSigC/dE  Sig(E)    E-E0       E ',my_mode)
   end if
 end if

 if (Sr%usepawu==0) then

   if (Sr%nsig_ab/=1) then
     write(msg,'(i5,9f8.3)')                       &
           iband,                                  &
           Sr%e0          (iband,ik_ibz,1)*Ha_eV,  &
           SUM(Sr%vxcme   (iband,ik_ibz,:))*Ha_eV, &
           SUM(Sr%sigxme  (iband,ik_ibz,:))*Ha_eV, &
      REAL(SUM(Sr%sigcmee0(iband,ik_ibz,:)))*Ha_eV,&
      REAL(Sr%ze0         (iband,ik_ibz,1)),       &
      REAL(SUM(Sr%dsigmee0(iband,ik_ibz,:))),      &
      REAL(SUM(Sr%sigmee  (iband,ik_ibz,:)))*Ha_eV,&
      REAL(Sr%degw        (iband,ik_ibz,1))*Ha_eV, &
      REAL(Sr%egw         (iband,ik_ibz,1))*Ha_eV
       call wrtout(my_unt,msg,my_mode)
     if (present(ydoc)) call ydoc%add_tabular_line(msg)
     if (verbose/=0) then
       write(msg,'(i5,9f8.3)')                        &
              iband,                                  &
              zero,                                   &
              zero,                                   &
              zero,                                   &
        AIMAG(SUM(Sr%sigcmee0(iband,ik_ibz,:)))*Ha_eV,&
        AIMAG(Sr%ze0         (iband,ik_ibz,1)),       &
        AIMAG(SUM(Sr%dsigmee0(iband,ik_ibz,:))),      &
        AIMAG(SUM(Sr%sigmee  (iband,ik_ibz,:)))*Ha_eV,&
        AIMAG(Sr%degw        (iband,ik_ibz,1))*Ha_eV, &
        AIMAG(Sr%egw         (iband,ik_ibz,1))*Ha_eV
       call wrtout(my_unt,msg,my_mode)
       if(present(ydoc)) call ydoc%add_tabular_line(msg)
     end if
  else
    write(msg,'(i5,9f8.3)')                    &
          iband,                               &
          Sr%e0      (iband,ik_ibz,isp)*Ha_eV, &
          Sr%vxcme   (iband,ik_ibz,isp)*Ha_eV, &
          Sr%sigxme  (iband,ik_ibz,isp)*Ha_eV, &
     REAL(Sr%sigcmee0(iband,ik_ibz,isp))*Ha_eV,&
     REAL(Sr%ze0     (iband,ik_ibz,isp)),      &
     REAL(Sr%dsigmee0(iband,ik_ibz,isp)),      &
     REAL(Sr%sigmee  (iband,ik_ibz,isp))*Ha_eV,&
     REAL(Sr%degw    (iband,ik_ibz,isp))*Ha_eV,&
     REAL(Sr%egw     (iband,ik_ibz,isp))*Ha_eV
    call wrtout(my_unt,msg,my_mode)
    if (present(ydoc)) call ydoc%add_tabular_line(msg)

    if (verbose/=0) then
      write(msg,'(i5,9f8.3)')                      &
              iband,                               &
              zero,                                &
              zero,                                &
              zero,                                &
        AIMAG(Sr%sigcmee0(iband,ik_ibz,isp))*Ha_eV,&
        AIMAG(Sr%ze0     (iband,ik_ibz,isp)),      &
        AIMAG(Sr%dsigmee0(iband,ik_ibz,isp)),      &
        AIMAG(Sr%sigmee  (iband,ik_ibz,isp))*Ha_eV,&
        AIMAG(Sr%degw    (iband,ik_ibz,isp))*Ha_eV,&
        AIMAG(Sr%egw     (iband,ik_ibz,isp))*Ha_eV
       call wrtout(my_unt,msg,my_mode)
       if (present(ydoc)) call ydoc%add_tabular_line(msg)
    end if
  end if

 else  ! PAW+U+GW calculation.
   ABI_CHECK(Sr%nsig_ab==1,'LDA+U with spinor not implemented')
   write(msg,'(i5,10f8.3)')                   &
         iband,                               &
         Sr%e0      (iband,ik_ibz,isp)*Ha_eV, &
         Sr%vxcme   (iband,ik_ibz,isp)*Ha_eV, &
         Sr%vUme    (iband,ik_ibz,isp)*Ha_eV, &
         Sr%sigxme  (iband,ik_ibz,isp)*Ha_eV, &
    REAL(Sr%sigcmee0(iband,ik_ibz,isp))*Ha_eV,&
    REAL(Sr%ze0     (iband,ik_ibz,isp)),      &
    REAL(Sr%dsigmee0(iband,ik_ibz,isp)),      &
    REAL(Sr%sigmee  (iband,ik_ibz,isp))*Ha_eV,&
    REAL(Sr%degw    (iband,ik_ibz,isp))*Ha_eV,&
    REAL(Sr%egw     (iband,ik_ibz,isp))*Ha_eV
   call wrtout(my_unt,msg,my_mode)
   if(present(ydoc)) call ydoc%add_tabular_line(msg)

   if (verbose/=0) then
     write(msg,'(i5,10f8.3)')                   &
           iband,                               &
           zero,                                &
           zero,                                &
           zero,                                &
           zero,                                &
     AIMAG(Sr%sigcmee0(iband,ik_ibz,isp))*Ha_eV,&
     AIMAG(Sr%ze0     (iband,ik_ibz,isp)),      &
     AIMAG(Sr%dsigmee0(iband,ik_ibz,isp)),      &
     AIMAG(Sr%sigmee  (iband,ik_ibz,isp))*Ha_eV,&
     AIMAG(Sr%degw    (iband,ik_ibz,isp))*Ha_eV,&
     AIMAG(Sr%egw     (iband,ik_ibz,isp))*Ha_eV
      call wrtout(my_unt,msg,my_mode)
     if(present(ydoc)) call ydoc%add_tabular_line(msg)
   end if
 end if

end subroutine print_Sigma_perturbative
!!***

!----------------------------------------------------------------------

!!****f* m_sigma/print_Sigma_QPSC
!! NAME
!!  print_Sigma_QPSC
!!
!! FUNCTION
!!  Write the results of the GW calculation in case of self-consistency
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_sigma
!!
!! CHILDREN
!!      findqg0,wfd_distribute_bands,wfd_update_bkstab,xmpi_sum
!!
!! SOURCE

subroutine print_Sigma_QPSC(Sr,ik_ibz,iband,isp,KS_BSt,unit,prtvol,mode_paral,ydoc)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iband,ik_ibz,isp
 integer,intent(in),optional :: prtvol,unit
 character(len=*),intent(in),optional :: mode_paral
 type(sigma_t),intent(in) :: Sr
 type(ebands_t),intent(in) :: KS_BSt
 type(yamldoc_t),intent(inout),optional :: ydoc

!Local variables-------------------------------
!scalars
 integer :: my_unt,verbose
 character(len=4) :: my_mode
 character(len=500) :: msg

! *********************************************************************

 my_unt =std_out; if (PRESENT(unit      )) my_unt =unit
 verbose=0      ; if (PRESENT(prtvol    )) verbose=prtvol
 my_mode='COLL' ; if (PRESENT(mode_paral)) my_mode=mode_paral

! write(msg,'(a)')&
!&   ' Band     E_lda   <Vxclda>   E(N-1)  <Hhartree>   SigX  SigC[E(N-1)]',&
!&   '    Z     dSigC/dE  Sig[E(N)]  DeltaE  E(N)_pert E(N)_diago'

 if (Sr%usepawu==0 .or. .TRUE.) then
   if (Sr%nsig_ab/=1) then
     write(msg,'(i5,12(2x,f8.3))')                       &
           iband,                                        &
           KS_BSt%eig     (iband,ik_ibz,1)*Ha_eV,        &
           SUM(Sr%vxcme   (iband,ik_ibz,:))*Ha_eV,       &
           Sr%e0          (iband,ik_ibz,1)*Ha_eV,        &
      REAL(SUM(Sr%hhartree(iband,iband,ik_ibz,:)))*Ha_eV,&
           SUM(Sr%sigxme  (iband,ik_ibz,:))*Ha_eV,       &
      REAL(SUM(Sr%sigcmee0(iband,ik_ibz,:)))*Ha_eV,      &
      REAL(Sr%ze0         (iband,ik_ibz,1)),             &
      REAL(SUM(Sr%dsigmee0(iband,ik_ibz,:))),            &
      REAL(SUM(Sr%sigmee  (iband,ik_ibz,:)))*Ha_eV,      &
      REAL(Sr%degw        (iband,ik_ibz,1))*Ha_eV,       &
      REAL(Sr%egw         (iband,ik_ibz,1))*Ha_eV,       &
           Sr%en_qp_diago (iband,ik_ibz,1)*Ha_eV
     call wrtout(my_unt,msg,my_mode)
     if (present(ydoc)) call ydoc%add_tabular_line(msg)

     write(msg,'(i5,12(2x,f8.3))')                        &
            iband,                                        &
            zero,                                         &
            zero,                                         &
            zero,                                         &
      AIMAG(SUM(Sr%hhartree(iband,iband,ik_ibz,:)))*Ha_eV,&
            zero,                                         &
      AIMAG(SUM(Sr%sigcmee0(iband,ik_ibz,:)))*Ha_eV,      &
      AIMAG(Sr%ze0         (iband,ik_ibz,1)),             &
      AIMAG(SUM(Sr%dsigmee0(iband,ik_ibz,:))),            &
      AIMAG(SUM(Sr%sigmee  (iband,ik_ibz,:)))*Ha_eV,      &
      AIMAG(Sr%degw        (iband,ik_ibz,1))*Ha_eV,       &
      AIMAG(Sr%egw         (iband,ik_ibz,1))*Ha_eV,       &
            zero
     if (verbose/=0) then
       call wrtout(my_unt,msg,my_mode)
       if( present(ydoc)) call ydoc%add_tabular_line(msg)
     end if
   else
     write(msg,'(i5,12(2x,f8.3))')                       &
           iband,                                        &
           KS_BSt%eig    (iband,ik_ibz,isp)*Ha_eV,       &
           Sr%vxcme      (iband,ik_ibz,isp)*Ha_eV,       &
           Sr%e0         (iband,ik_ibz,isp)*Ha_eV,       &
      REAL(Sr%hhartree   (iband,iband,ik_ibz,isp))*Ha_eV,&
           Sr%sigxme     (iband,ik_ibz,isp)*Ha_eV,       &
      REAL(Sr%sigcmee0   (iband,ik_ibz,isp))*Ha_eV,      &
      REAL(Sr%ze0        (iband,ik_ibz,isp)),            &
      REAL(Sr%dsigmee0   (iband,ik_ibz,isp)),            &
      REAL(Sr%sigmee     (iband,ik_ibz,isp))*Ha_eV,      &
      REAL(Sr%degw       (iband,ik_ibz,isp))*Ha_eV,      &
      REAL(Sr%egw        (iband,ik_ibz,isp))*Ha_eV,      &
           Sr%en_qp_diago(iband,ik_ibz,isp)*Ha_eV
     call wrtout(my_unt,msg,my_mode)
     if (present(ydoc)) call ydoc%add_tabular_line(msg)

     write(msg,'(i5,12(2x,f8.3))')                       &
            iband,                                       &
            zero,                                        &
            zero,                                        &
            zero,                                        &
      AIMAG(Sr%hhartree  (iband,iband,ik_ibz,isp))*Ha_eV,&
            zero,                                        &
      AIMAG(Sr%sigcmee0   (iband,ik_ibz,isp))*Ha_eV,     &
      AIMAG(Sr%ze0        (iband,ik_ibz,isp)),           &
      AIMAG(Sr%dsigmee0   (iband,ik_ibz,isp)),           &
      AIMAG(Sr%sigmee     (iband,ik_ibz,isp))*Ha_eV,     &
      AIMAG(Sr%degw       (iband,ik_ibz,isp))*Ha_eV,     &
      AIMAG(Sr%egw        (iband,ik_ibz,isp))*Ha_eV,     &
            zero
     if (verbose/=0) then
       call wrtout(my_unt,msg,my_mode)
       if (present(ydoc)) call ydoc%add_tabular_line(msg)
     end if
   end if

 else ! PAW+U+GW calculation.
   MSG_ERROR("PAW+U+GW not yet implemented")
 end if

end subroutine print_Sigma_QPSC
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
!! usepawu= /=0 if we used LDA+U as starting point (only for PAW)
!!
!! OUTPUT
!!
!! TODO
!!  Write documentation.
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      findqg0,wfd_distribute_bands,wfd_update_bkstab,xmpi_sum
!!
!! SOURCE

subroutine sigma_init(Sigp,nkibz,usepawu,Sr)

!Arguments ------------------------------------
 integer,intent(in) :: nkibz,usepawu
!scalars
 type(sigparams_t),intent(in) :: Sigp
 type(sigma_t),intent(inout) :: Sr

!Local variables-------------------------------
!scalars
 integer :: b1gw,b2gw,mod10

! *************************************************************************

 !@sigma_t
 mod10=MOD(Sigp%gwcalctyp,10)

 ! Copy important dimensions
 Sr%nkptgw     =Sigp%nkptgw
 Sr%gwcalctyp  =Sigp%gwcalctyp
 Sr%deltae     =Sigp%deltae
 Sr%maxomega4sd=Sigp%maxomega4sd
 Sr%maxomega_r =Sigp%maxomega_r
 Sr%scissor_ene=Sigp%mbpt_sciss

 !FIXME this should be done in sigma_allocate
 ABI_MALLOC(Sr%minbnd,(Sr%nkptgw,Sigp%nsppol))
 ABI_MALLOC(Sr%maxbnd,(Sr%nkptgw,Sigp%nsppol))
 Sr%minbnd=Sigp%minbnd; Sr%maxbnd=Sigp%maxbnd
 ABI_MALLOC(Sr%kptgw,(3,Sr%nkptgw))
 Sr%kptgw=Sigp%kptgw

 Sr%b1gw     =Sigp%minbdgw ! * min and Max GW band index over k and spin.
 Sr%b2gw     =Sigp%maxbdgw !   Used to dimension arrays.
 Sr%nbnds    =Sigp%nbnds
 Sr%nkibz    =nkibz
 Sr%nsppol   =Sigp%nsppol
 Sr%nsig_ab  =Sigp%nsig_ab
 Sr%nomega_r =Sigp%nomegasr  !FIXME change name
 Sr%nomega_i =Sigp%nomegasi
 Sr%nomega4sd=Sigp%nomegasrd
 Sr%usepawu  =usepawu

 !======================================================
 ! === Allocate arrays in the sigma_t datatype ===
 !======================================================
 b1gw=Sr%b1gw
 b2gw=Sr%b2gw

 !TODO write routine to allocate all this stuff

 ! hhartree(b1,b2,k,s)= <b1,k,s|T+v_{loc}+v_{nl}+v_{H}|b2,k,s>
 ABI_CALLOC(Sr%hhartree,(b1gw:b2gw,b1gw:b2gw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))

 ! QP amplitudes and energies
 ABI_CALLOC(Sr%en_qp_diago,(Sr%nbnds,Sr%nkibz,Sr%nsppol))
 ABI_CALLOC(Sr%eigvec_qp,(Sr%nbnds,Sr%nbnds,Sr%nkibz,Sr%nsppol))

 ! Dont know if it is better to do this here or in the sigma
 ! * Initialize with KS wavefunctions and energies
 !do ib=1,Sr%nbnds
 ! Sr%en_qp_diago(ib,:,:)=en(:,ib,:)
 ! Sr%eigvec_qp(ib,ib,:,:)=cone
 !end do

 ABI_CALLOC(Sr%vxcme, (b1gw:b2gw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))
 ABI_CALLOC(Sr%vUme, (b1gw:b2gw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))
 ABI_CALLOC(Sr%sigxme, (b1gw:b2gw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))
 ABI_CALLOC(Sr%x_mat, (b1gw:b2gw,b1gw:b2gw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))
 ABI_CALLOC(Sr%sigcme, (b1gw:b2gw,Sr%nkibz,Sr%nomega_r,Sr%nsppol*Sr%nsig_ab))
 ABI_CALLOC(Sr%sigxcme, (b1gw:b2gw,Sr%nkibz,Sr%nomega_r,Sr%nsppol*Sr%nsig_ab))
 ABI_CALLOC(Sr%sigcmee0, (b1gw:b2gw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))
 ABI_CALLOC(Sr%ze0, (b1gw:b2gw,Sr%nkibz,Sr%nsppol))
 ABI_CALLOC(Sr%dsigmee0, (b1gw:b2gw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))
 ABI_CALLOC(Sr%sigmee, (b1gw:b2gw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))
 ABI_CALLOC(Sr%degw, (b1gw:b2gw,Sr%nkibz,Sr%nsppol))
 ABI_CALLOC(Sr%e0, (Sr%nbnds,Sr%nkibz,Sr%nsppol))
 ABI_CALLOC(Sr%egw, (Sr%nbnds,Sr%nkibz,Sr%nsppol))
 ABI_CALLOC(Sr%e0gap, (Sr%nkibz,Sr%nsppol))
 ABI_CALLOC(Sr%degwgap, (Sr%nkibz,Sr%nsppol))
 ABI_CALLOC(Sr%egwgap, (Sr%nkibz,Sr%nsppol))

 ! These quantities are used to evaluate $\Sigma(E)$ around the KS\QP eigenvalue
 ABI_CALLOC(Sr%omega4sd,(b1gw:b2gw,Sr%nkibz,Sr%nomega4sd,Sr%nsppol))
 ABI_CALLOC(Sr%sigcme4sd,(b1gw:b2gw,Sr%nkibz,Sr%nomega4sd,Sr%nsppol*Sr%nsig_ab))
 ABI_CALLOC(Sr%sigxcme4sd,(b1gw:b2gw,Sr%nkibz,Sr%nomega4sd,Sr%nsppol*Sr%nsig_ab))

 !TODO Find  better treatment
 ! Mesh along the real axis.
 if (Sr%nomega_r>0) then
   ABI_MALLOC(Sr%omega_r,(Sr%nomega_r))
   Sr%omega_r(:)=Sigp%omega_r(:)
 end if

 ! Analytical Continuation
 ! FIXME omegasi should not be in Sigp% here we should construct the mesh
 if (mod10==1) then
   ABI_MALLOC(Sr%omega_i,(Sr%nomega_i))
   Sr%omega_i=Sigp%omegasi
   ABI_MALLOC(Sr%sigcmesi ,(b1gw:b2gw,Sr%nkibz,Sr%nomega_i,Sr%nsppol*Sr%nsig_ab))
   ABI_MALLOC(Sr%sigxcmesi,(b1gw:b2gw,Sr%nkibz,Sr%nomega_i,Sr%nsppol*Sr%nsig_ab))
   Sr%sigcmesi=czero; Sr%sigxcmesi=czero
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
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      findqg0,wfd_distribute_bands,wfd_update_bkstab,xmpi_sum
!!
!! SOURCE

subroutine sigma_free(Sr)

!Arguments ------------------------------------
!scalars
 type(sigma_t),intent(inout) :: Sr

! *************************************************************************

 !@sigma_t
!integer
 ABI_SFREE(Sr%maxbnd)
 ABI_SFREE(Sr%minbnd)

!real
 ABI_SFREE(Sr%degwgap)
 ABI_SFREE(Sr%egwgap)
 ABI_SFREE(Sr%en_qp_diago)
 ABI_SFREE(Sr%e0)
 ABI_SFREE(Sr%e0gap)
 ABI_SFREE(Sr%omega_r)
 ABI_SFREE(Sr%kptgw)
 ABI_SFREE(Sr%sigxme)
 ABI_SFREE(Sr%x_mat)
 ABI_SFREE(Sr%vxcme)
 ABI_SFREE(Sr%vUme)

!complex
 ABI_SFREE(Sr%degw)
 ABI_SFREE(Sr%dsigmee0)
 ABI_SFREE(Sr%egw)
 ABI_SFREE(Sr%eigvec_qp)
 ABI_SFREE(Sr%m_lda_to_qp)
 ABI_SFREE(Sr%hhartree)
 ABI_SFREE(Sr%sigcme)
 ABI_SFREE(Sr%sigmee)
 ABI_SFREE(Sr%sigcmee0)
 ABI_SFREE(Sr%sigcmesi)
 ABI_SFREE(Sr%sigcme4sd)
 ABI_SFREE(Sr%sigxcme)
 ABI_SFREE(Sr%sigxcmesi)
 ABI_SFREE(Sr%sigxcme4sd)
 ABI_SFREE(Sr%ze0)
 ABI_SFREE(Sr%omega_i)
 ABI_SFREE(Sr%omega4sd)

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
!! PARENTS
!!
!! SOURCE

pure function sigma_get_exene(sigma,kmesh,bands) result(ex_energy)

!Arguments ------------------------------------
!scalars
 real(dp) :: ex_energy
 type(sigma_t),intent(in) :: sigma
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
       if (sigma%nsig_ab==1) then
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

!!****f* m_sigma/find_wpoles_for_cd
!! NAME
!!  find_wpoles_for_cd
!!
!! FUNCTION
!!  Find the max frequency needed to account for all the poles
!!  of GW used in the contour deformation technique.
!!
!! INPUTS
!!  Sigp=sigparams_t
!!
!! OUTPUT
!!  omega_max
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!      findqg0,wfd_distribute_bands,wfd_update_bkstab,xmpi_sum
!!
!! SOURCE

subroutine find_wpoles_for_cd(Sigp,Sr,Kmesh,BSt,omega_max)

!Arguments ------------------------------------
!scalars
 type(sigparams_t),intent(in) :: Sigp
 type(sigma_t),intent(in) :: Sr
 type(ebands_t),intent(in) :: Bst
 type(kmesh_t),intent(in) :: Kmesh
 real(dp),intent(out) :: omega_max

!Local variables-------------------------------
!scalars
 integer :: spin,ik_ibz,band_gr,bgw_start,bgw_stop,io,ioe0j
 integer :: ikgw,ikgw_ibz,ikgw_bz,band_gw,nomega_tot
 real(dp) :: e_green,e_screen,theta_mu_minus_e0i,e_qp
 real(dp) :: fact_sp
 !character(len=500) :: msg
!arrays
 real(dp),allocatable :: omegame0i(:)

! *************************************************************************

 omega_max = smallest_real
 !
 ! === Normalization of theta_mu_minus_e0i ===
 ! * If nsppol==2, qp_occ $\in [0,1]$
 fact_sp=one
 if (Bst%nsppol==1) then
   fact_sp=half; if (Bst%nspinor==2) fact_sp=one
 end if
 !
 ! Total number of frequencies for sigma (Spectral function + mesh for the derivative).
 nomega_tot=Sr%nomega_r+Sr%nomega4sd
 ABI_MALLOC(omegame0i,(nomega_tot))

 ioe0j=Sr%nomega4sd/2+1
 !
 ! Loop over bands used to construct the Green function.
 do spin=1,Bst%nsppol
   do ik_ibz=1,Bst%nkpt
     do band_gr=1,Bst%nband(ik_ibz+(spin-1)*Bst%nkpt)
       e_green           = Bst%eig(band_gr,ik_ibz,spin)
       theta_mu_minus_e0i= Bst%occ(band_gr,ik_ibz,spin)*fact_sp
       !
       ! Loop over GW states.
       do ikgw=1,Sigp%nkptgw
         bgw_start=Sigp%minbnd(ikgw,spin)
         bgw_stop =Sigp%minbnd(ikgw,spin)
         ikgw_bz  =Sigp%kptgw2bz(ikgw_bz)
         ikgw_ibz =Kmesh%tab(ikgw_bz)

         do band_gw=bgw_start,bgw_stop
           e_qp      = Bst%eig(band_gw,ikgw_ibz,spin)
           !
           ! Get frequencies $\omega$-\epsilon_in$ to evaluate  $d\Sigma/dE$, note the spin
           ! subtract e_KS since we have stored e_KS+ Delta \omega in Sr%omega4sd, not required for AC
           if (Sr%nomega_r>0) omegame0i(1:Sr%nomega_r)=DBLE(Sigp%omega_r(1:Sr%nomega_r))-e_green
           do io=Sr%nomega_r+1,nomega_tot
             !omegame0i(io)=DBLE(Sr%omega4sd(band_gw,ikgw_ibz,io-Sr%nomega_r,spin)) - e_green
             !Sr%omega4sd(jb,ik_ibz,io,spin)=Sr%egw(jb,ik_ibz,spin)+Sigp%deltae*(io-ioe0j)
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
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!
!! SOURCE

integer function sigma_ncwrite(Sigp,Er,Sr,ncid) result (ncerr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 type(sigparams_t),target,intent(in) :: Sigp
 type(Epsilonm1_results),target,intent(in) :: Er
 type(sigma_t),target,intent(in) :: Sr

!Local variables ---------------------------------------
#ifdef HAVE_NETCDF
!scalars
 integer :: nbgw,ndim_sig,b1gw,b2gw,cplex
 !character(len=500) :: msg
!arrays
 real(dp),allocatable :: rdata2(:,:),rdata4(:,:,:,:),rdata5(:,:,:,:,:)

! *************************************************************************

 !@sigma_t
 cplex=2; b1gw=Sr%b1gw; b2gw=Sr%b2gw; nbgw=b2gw-b1gw+1
 ndim_sig=Sr%nsppol*Sr%nsig_ab

 ncerr = nctk_def_dims(ncid, [nctkdim_t("cplex", cplex), nctkdim_t("b1gw", sr%b1gw), nctkdim_t("b2gw", sr%b2gw),&
   nctkdim_t("nbgw", nbgw), nctkdim_t("nkptgw", sr%nkptgw), nctkdim_t("ndim_sig", ndim_sig), &
   nctkdim_t("nomega4sd", sr%nomega4sd), nctkdim_t("nsig_ab", sr%nsig_ab)], defmode=.True.)
 NCF_CHECK(ncerr)

 ! No. of real frequencies, might be zero.
 if (Sr%nomega_r>0) then
   NCF_CHECK(nctk_def_dims(ncid, nctkdim_t("nomega_r", Sr%nomega_r)))
 end if

 ! No. of imaginary frequencies, might be zero.
 if (Sr%nomega_i>0) then
   NCF_CHECK(nctk_def_dims(ncid, nctkdim_t("nomega_i", Sr%nomega_i)))
 end if

 ! =======================
 ! == Define variables ===
 ! =======================
 ! parameters of the calculation.

 ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: 'sigma_nband', 'scr_nband', 'gwcalctyp', 'usepawu'])
 NCF_CHECK(ncerr)

 ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: &
&  'ecutwfn', 'ecuteps', 'ecutsigx', 'omegasrdmax', 'deltae', 'omegasrmax', 'scissor_ene'])
 NCF_CHECK(ncerr)

 ! TODO: Decrease size of file: Remove arrays whose size scale as mband ** 2
 ! especially those that are not commonly used e.g. hhartree.
 ncerr = nctk_def_arrays(ncid, [&
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
   nctkarr_t('eigvec_qp', "dp",'cplex, max_number_of_states, max_number_of_states, number_of_kpoints, number_of_spins'),&
   nctkarr_t('hhartree', "dp",'cplex, nbgw, nbgw, number_of_kpoints, ndim_sig'),&
   nctkarr_t('sigmee', "dp", 'cplex, nbgw, number_of_kpoints, ndim_sig'),&
   nctkarr_t('sigcmee0', "dp",'cplex, nbgw, number_of_kpoints, ndim_sig'),&
   nctkarr_t('sigcmesi', "dp",'cplex, nbgw, number_of_kpoints, ndim_sig'),&
   nctkarr_t('sigcme4sd', "dp",'cplex, nbgw, number_of_kpoints, nomega4sd, ndim_sig'),&
   nctkarr_t('sigxcme4sd', "dp", 'cplex, nbgw, number_of_kpoints, nomega4sd, ndim_sig'),&
   nctkarr_t('ze0',"dp", 'cplex, nbgw, number_of_kpoints, number_of_spins'),&
   nctkarr_t('omega4sd', "dp", 'cplex, nbgw, number_of_kpoints, nomega4sd, number_of_spins')])
 NCF_CHECK(ncerr)

 if (Sr%usepawu==0) then
   ncerr = nctk_def_arrays(ncid, nctkarr_t("vUme", "dp", 'nbgw, number_of_kpoints, ndim_sig'))
   NCF_CHECK(ncerr)
 end if

 if (Sr%nomega_r>0) then
   ncerr = nctk_def_arrays(ncid, [&
     nctkarr_t('omega_r', "dp", "nomega_r"),&
     nctkarr_t('sigcme', "dp", 'cplex, nbgw, number_of_kpoints, nomega_r, ndim_sig'),&
     nctkarr_t('sigxcme', "dp", 'cplex, nbgw, number_of_kpoints, nomega_r, ndim_sig')])
   NCF_CHECK(ncerr)
 end if

 if (Sr%nomega_i>0) then
   ncerr = nctk_def_arrays(ncid, [&
     nctkarr_t('sigxcmesi', "dp", 'cplex, nbgw, number_of_kpoints, nomega_i, ndim_sig'),&
     nctkarr_t('omega_i', "dp", 'cplex, nomega_i')])
   NCF_CHECK(ncerr)
 end if

 if (allocated(sr%m_lda_to_qp)) then
   ncerr = nctk_def_arrays(ncid, [nctkarr_t('m_lda_to_qp', "dp", &
&       "cplex, max_number_of_states, max_number_of_states, number_of_kpoints, number_of_spins")])
   NCF_CHECK(ncerr)
 end if

 ! =====================
 ! === Start writing ===
 ! =====================
 ! vid is a small helper function to get the varid from ncid
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, vid('ecutwfn'), Sigp%ecutwfn))
 NCF_CHECK(nf90_put_var(ncid, vid('ecuteps'), Sigp%ecuteps))
 NCF_CHECK(nf90_put_var(ncid, vid('ecutsigx'), Sigp%ecutsigx))
 NCF_CHECK(nf90_put_var(ncid, vid('sigma_nband'), Sigp%nbnds))
 NCF_CHECK(nf90_put_var(ncid, vid('scr_nband'), Er%Hscr%nbnds_used))
 NCF_CHECK(nf90_put_var(ncid, vid('gwcalctyp'), Sr%gwcalctyp))
 NCF_CHECK(nf90_put_var(ncid, vid('usepawu'), Sr%usepawu))
 NCF_CHECK(nf90_put_var(ncid, vid('kptgw'), Sr%kptgw))
 NCF_CHECK(nf90_put_var(ncid, vid('minbnd'), Sr%minbnd))
 NCF_CHECK(nf90_put_var(ncid, vid('maxbnd'),Sr%maxbnd))

 NCF_CHECK(nf90_put_var(ncid, vid('omegasrdmax'), Sr%maxomega4sd*Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('deltae'), Sr%deltae*Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('omegasrmax'), Sr%maxomega_r*Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('scissor_ene'), Sr%scissor_ene*Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('degwgap'), Sr%degwgap*Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('egwgap'), Sr%egwgap*Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('en_qp_diago'), Sr%en_qp_diago*Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('e0'), Sr%e0*Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('e0gap'), Sr%e0gap*Ha_eV))

 if (Sr%nomega_r>0) then
   NCF_CHECK(nf90_put_var(ncid, vid('omega_r'), Sr%omega_r*Ha_eV))
 end if

 NCF_CHECK(nf90_put_var(ncid, vid('sigxme'), Sr%sigxme*Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('vxcme'), Sr%vxcme*Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('vUme'), Sr%vUme*Ha_eV))

 ! * Have to transfer complex arrays
 ABI_MALLOC(rdata4,(cplex,b1gw:b2gw,Sr%nkibz,Sr%nsppol))
 rdata4=c2r(Sr%degw)
 NCF_CHECK(nf90_put_var(ncid, vid('degw'), rdata4*Ha_eV))
 ABI_FREE(rdata4)

 ABI_MALLOC(rdata4,(cplex,b1gw:b2gw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))
 rdata4=c2r(Sr%dsigmee0)
 NCF_CHECK(nf90_put_var(ncid, vid('dsigmee0'), rdata4))
 ABI_FREE(rdata4)

 ABI_MALLOC(rdata4,(cplex,Sr%nbnds,Sr%nkibz,Sr%nsppol))
 rdata4=c2r(Sr%egw)
 NCF_CHECK(nf90_put_var(ncid, vid('egw'), rdata4*Ha_eV))
 ABI_FREE(rdata4)

 ABI_MALLOC(rdata5,(cplex,Sr%nbnds,Sr%nbnds,Sr%nkibz,Sr%nsppol))
 rdata5=c2r(Sr%eigvec_qp)
 NCF_CHECK(nf90_put_var(ncid, vid('eigvec_qp'), rdata5))
 ABI_FREE(rdata5)

 ABI_MALLOC(rdata5,(cplex,nbgw,nbgw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))
 rdata5=c2r(Sr%hhartree)
 NCF_CHECK(nf90_put_var(ncid, vid('hhartree'), rdata5*Ha_eV))
 ABI_FREE(rdata5)

 if (Sr%nomega_r>0) then
   ABI_MALLOC(rdata5,(cplex,nbgw,Sr%nkibz,Sr%nomega_r,Sr%nsppol*Sr%nsig_ab))
   rdata5=c2r(Sr%sigcme)
   NCF_CHECK(nf90_put_var(ncid, vid('sigcme'), rdata5*Ha_eV))
   ABI_FREE(rdata5)
 end if

 ABI_MALLOC(rdata4,(cplex,nbgw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))
 rdata4=c2r(Sr%sigmee)
 NCF_CHECK(nf90_put_var(ncid, vid('sigmee'), rdata4*Ha_eV))
 ABI_FREE(rdata4)

 ABI_MALLOC(rdata4,(cplex,nbgw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))
 rdata4=c2r(Sr%sigcmee0)
 NCF_CHECK(nf90_put_var(ncid, vid('sigcmee0'), rdata4*Ha_eV))
 ABI_FREE(rdata4)

 if (Sr%nomega_i>0) then
  ABI_MALLOC(rdata5,(cplex,nbgw,Sr%nkibz,Sr%nomega_i,Sr%nsppol*Sr%nsig_ab))
  rdata5=c2r(Sr%sigcmesi)
  NCF_CHECK(nf90_put_var(ncid, vid('sigcmesi'), rdata5*Ha_eV))
  ABI_FREE(rdata5)
 end if

 ABI_MALLOC(rdata5,(cplex,nbgw,Sr%nkibz,Sr%nomega4sd,Sr%nsppol*Sr%nsig_ab))
 rdata5=c2r(Sr%sigcme4sd)
 NCF_CHECK(nf90_put_var(ncid, vid('sigcme4sd'), rdata5*Ha_eV))
 ABI_FREE(rdata5)

 if (Sr%nomega_r>0) then
   ABI_MALLOC(rdata5,(cplex,nbgw,Sr%nkibz,Sr%nomega_r,Sr%nsppol*Sr%nsig_ab))
   rdata5=c2r(Sr%sigxcme)
   NCF_CHECK(nf90_put_var(ncid, vid('sigxcme'), rdata5*Ha_eV))
   ABI_FREE(rdata5)
 end if

 if (Sr%nomega_i>0) then
   ABI_MALLOC(rdata5,(cplex,nbgw,Sr%nkibz,Sr%nomega_i,Sr%nsppol*Sr%nsig_ab))
   rdata5=c2r(Sr%sigxcmesi)
   NCF_CHECK(nf90_put_var(ncid, vid('sigxcmesi'), rdata5*Ha_eV))
   ABI_FREE(rdata5)
 end if

 if (allocated(sr%m_lda_to_qp)) then
   ABI_MALLOC(rdata5,(cplex,Sr%nbnds,Sr%nbnds,Sr%nkibz,Sr%nsppol))
   rdata5=c2r(Sr%m_lda_to_qp)
   NCF_CHECK(nf90_put_var(ncid, vid('m_lda_to_qp'), rdata5))
   ABI_FREE(rdata5)
 end if

 ABI_MALLOC(rdata5,(cplex,nbgw,Sr%nkibz,Sr%nomega4sd,Sr%nsppol*Sr%nsig_ab))
 rdata5=c2r(Sr%sigxcme4sd)
 NCF_CHECK(nf90_put_var(ncid, vid('sigxcme4sd'), rdata5*Ha_eV))
 ABI_FREE(rdata5)

 ABI_MALLOC(rdata4,(cplex,nbgw,Sr%nkibz,Sr%nsppol))
 rdata4=c2r(Sr%ze0)
 NCF_CHECK(nf90_put_var(ncid, vid('ze0'), rdata4))
 ABI_FREE(rdata4)

 if (Sr%nomega_i>0) then
   ABI_MALLOC(rdata2,(cplex,Sr%nomega_i))
   rdata2=c2r(Sr%omega_i)
   NCF_CHECK(nf90_put_var(ncid, vid('omega_i'), rdata2*Ha_eV))
   ABI_FREE(rdata2)
 end if

 ABI_MALLOC(rdata5,(cplex,nbgw,Sr%nkibz,Sr%nomega4sd,Sr%nsppol))
 rdata5=c2r(Sr%omega4sd)
 NCF_CHECK(nf90_put_var(ncid, vid('omega4sd'), rdata5*Ha_eV))
 ABI_FREE(rdata5)

#else
  MSG_ERROR('netcdf support is not activated.')
#endif

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
!! Wfd(wfd_t),intent(inout) ::
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
!! PARENTS
!!      calc_sigc_me,calc_sigx_me,cohsex_me
!!
!! CHILDREN
!!      findqg0,wfd_distribute_bands,wfd_update_bkstab,xmpi_sum
!!
!! SOURCE

subroutine sigma_distribute_bks(Wfd,Kmesh,Ltg_kgw,Qmesh,nsppol,can_symmetrize,kptgw,mg0,my_nbks,proc_distrb,got,bks_mask,global)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsppol
 integer,intent(out) :: my_nbks
 logical,optional,intent(in) :: global
 type(kmesh_t),intent(in) :: Kmesh,Qmesh
 type(littlegroup_t),intent(in) :: Ltg_kgw
 type(wfd_t),intent(inout) :: Wfd
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
 !integer :: test(Wfd%mband,Kmesh%nbz,nsppol)
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
       kgwmk= kptgw-Kmesh%bz(:,ik_bz) ! kptgw must be inside the BZ
       call findqg0(iq_bz,g0,kgwmk,Qmesh%nbz,Qmesh%bz,mG0) ! Identify q_bz and G0 where q_bz+G0=k_gw-k_bz
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
     where (proc_distrb == xmpi_undefined_rank+1)
       proc_distrb = 0
     end where
     call xmpi_sum(proc_distrb,Wfd%comm,ierr)
     where (proc_distrb == 0)
       proc_distrb = xmpi_undefined_rank
     elsewhere
       proc_distrb = proc_distrb -1
     end where
     !where (proc_distrb /= xmpi_undefined_rank)
     !  ltest = (ANY(proc_distrb == (/(ii,ii=0,Wfd%nproc-1)/)))
     !end where
     !if (.not.ltest) then
     !  write(std_out,*)proc_distrb
     !  MSG_BUG("Bug in the generation of proc_distrb table")
     !end if
   end if
 end if

 my_nbks = COUNT(proc_distrb==Wfd%my_rank)
 if (PRESENT(got)) got=get_more

end subroutine sigma_distribute_bks
!!***

!----------------------------------------------------------------------

END MODULE m_sigma
!!***
