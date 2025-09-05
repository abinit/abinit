!!****m* ABINIT/m_ppmodel
!! NAME
!! m_ppmodel
!!
!! FUNCTION
!!  Module containing the definition of the ppmodel_t used to deal with the plasmonpole technique.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2025 ABINIT group (MG, GMR, VO, LR, RWG, RS)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_ppmodel

 use defs_basis
 use m_errors
 use m_abicore
 use m_array
 use m_linalg_interfaces
 use m_distribfft
 use m_yaml

 use defs_abitypes,    only : MPI_type
 use m_fstrings,       only : sjoin, itoa, ktoa
 use m_hide_lapack,    only : xhegv
 use m_gwdefs,         only : GW_Q0_DEFAULT, czero_gw
 use m_crystal,        only : crystal_t
 use m_bz_mesh,        only : kmesh_t
 use m_gsphere,        only : gsphere_t
 use m_vcoul,          only : vcoul_t
 use m_qplusg,         only : cmod_qpg
 use m_fft_mesh,       only : g2ifft
 use m_fft,            only : fourdp
 use m_mpinfo,         only : destroy_mpi_enreg, initmpi_seq
 use m_pstat,          only : pstat_proc

 implicit none

 private
!!***

 integer,public,parameter :: PPM_NONE            = 0
 integer,public,parameter :: PPM_GODBY_NEEDS     = 1
 integer,public,parameter :: PPM_HYBERTSEN_LOUIE = 2
 integer,public,parameter :: PPM_LINDEN_HORSH    = 3
 integer,public,parameter :: PPM_ENGEL_FARID     = 4

 ! Flags giving the status of the plasmon-pole tables
 integer,public,parameter :: PPM_NOTAB         = 0
 integer,public,parameter :: PPM_TAB_ALLOCATED = 1
 integer,public,parameter :: PPM_TAB_STORED    = 2

!----------------------------------------------------------------------

!!****t* m_ppmodel/ppmodel_t
!! NAME
!! ppmodel_t
!!
!! FUNCTION
!!  This datatype gathers all the information on the Plasmonpole technique used in the calculations
!!
!! SOURCE

 type,public :: ppmodel_t

   integer :: dm2_botsq
   ! =npwc if ppmodel=1,2
   ! =1    if ppmodel=3,4

   integer :: dm_eig
   ! =0    if ppmodel=1,2,4
   ! =npwc if ppmodel=3

   integer :: dm2_otq
   ! =npwc if ppmodel=1,2
   ! =1    if ppmodel=3,4

   integer :: invalid_freq
   ! what to do when PPM frequencies are invalid.

   integer :: model
   ! The type of Plasmonpole model.

   integer :: mqmem
   ! =nqibz if in-core solution.
   ! =0 for out-of-core for which the last dimension in the ppm arrays has size 1.

   integer :: nqibz
   ! Number of q-points in the IBZ

   integer :: npwc
   ! Number of G vectors in $\tilde \epsilon $

   integer :: userho
   ! 1 if the ppmodel requires rho(G).

   integer :: iq_bz = 0
   ! The index of the q-point in the BZ that is referenced by the internal pointer.
   ! when we perform a symmetrization from the IBZ to the BZ.

   real(dp) :: drude_plsmf
   ! Drude plasma frequency

   real(dp) :: force_plsmf = zero
   ! Force plasma frequency to that set by ppmfreq (not used if set zero).

   ! arrays
   logical,allocatable :: keep_qibz(:)
   ! (nqibz)
   ! .TRUE. if the ppmodel tables for this q in the IBZ are kept in memory.

   integer,allocatable :: has_qibz(:)
   ! Flag defining the status of the tables for the different q. See the PPM_TAB flags.

   complex(gwpc),allocatable :: bigomegatwsq_qbz_vals(:,:)
   ! (Points|Stores) the symmetrized plasmon pole parameters $\tilde\Omega^2_{G Gp}(q_bz)$.

   complex(gwpc),allocatable :: omegatw_qbz_vals(:,:)
   ! (Points|Stores) the symmetrized plasmon pole parameters $\tilde\omega_{G Gp}(q_bz)$.

   complex(gwpc),allocatable :: eigpot_qbz_vals(:,:)
   ! (Points|Stores) the eigvectors of the symmetrized inverse dielectric matrix.

   type(array2_gwpc_t),allocatable :: bigomegatwsq(:)
   ! (nqibz)%value(npwc,dm2_botsq)
   ! Plasmon pole parameters $\tilde\Omega^2_{G Gp}(q)$.

   type(array2_gwpc_t),allocatable :: omegatw(:)
   ! (nqibz)%value(npwc,dm2_otq)
   ! Plasmon pole parameters $\tilde\omega_{G Gp}(q)$.

   type(array2_gwpc_t),allocatable :: eigpot(:)
   ! (nqibz)%value(dm_eig,dm_eig)
   ! Eigvectors of the symmetrized inverse dielectric matrix.

contains

   procedure :: get_qbz => ppm_get_qbz
     ! Symmetrize the ppm parameters in the BZ.

   procedure :: init => ppm_init
     ! Initialize dimensions and pointers

   procedure :: free => ppm_free
     ! Free dynamic memory.

   procedure :: setup => ppm_setup
     ! Main Driver

   procedure :: new_setup => ppm_new_setup
     ! New Main Driver

   procedure :: print => ppm_print
     ! Print info on object

   procedure :: calc_sigc => ppm_calc_sigc
     ! Matrix elements of the correlated self-energy with ppmodel.

   procedure :: rotate_iqbz => ppm_rotate_iqbz

   procedure :: malloc_iqibz => ppm_malloc_iqibz

   procedure :: table_free_iqibz => ppm_table_free_iqibz

   procedure :: get_eigenvalues => ppm_get_eigenvalues

   procedure :: getem1 => ppm_getem1
     ! Reconstruct e^{-1}(w) from ppm.

   procedure :: getem1_one_ggp => ppm_getem1_one_ggp
     ! Reconstruct e^{-1}(w) from ppm for one (G,G') pair

 end type ppmodel_t

 public :: cqratio
!!***

contains
!!***

!!****f* m_ppmodel/ppm_get_qbz
!! NAME
!!  ppm_get_qbz
!!
!! FUNCTION
!!  Compute plasmonpole matrix elements for q in the BZ from the symmetrical image in the IBZ
!!
!! INPUTS
!!  Gsph<gsphere_t>=data related to the G-sphere
!!  Qmesh<kmesh_t>=Info on the q-mesh
!!  iq_bz=Index of the q-point in the BZ where ppmodel parameters have to be symmetrized
!!
!! OUTPUT
!!  botsq
!!  otq
!!  eig (only if ppm%ppmodel==3)
!!
!! NOTES
!!  In the present implementation we are not considering a possible umklapp vector G0.
!!  In this case,indeed, the equation is different since we have to consider G-G0.
!!  There is however a check in sigma
!!
!!  * Remember the symmetry properties of \tilde\espilon^{-1}
!!
!!    If q_bz=Sq_ibz+G0:
!!
!!      $\epsilon^{-1}_{SG1-G0,SG2-G0}(q_bz) = e^{+iS(G2-G1).\tau}\epsilon^{-1}_{G1,G2)}(q)
!!
!!    If time-reversal symmetry can be used then:
!!
!!       $\epsilon^{-1}_{G1,G2}(-q_bz) = e^{+i(G1-G2).\tau}\epsilon^{-1}_{-S^{-1}(G1+Go),-S^{-1}(G2+G0)}^*(q)
!!
!! * Note that eig is used only if ppm%model==3
!!
!! SOURCE

subroutine ppm_get_qbz(ppm, Gsph, Qmesh, iq_bz, botsq, otq, eig)

!Arguments ------------------------------------
!scalars
 class(ppmodel_t),target,intent(inout) :: ppm
 integer,intent(in) :: iq_bz
 type(gsphere_t),target,intent(in) :: Gsph
 type(kmesh_t),intent(in) :: Qmesh
 complex(gwpc),allocatable,intent(out) :: botsq(:,:),otq(:,:),eig(:,:)

!Local variables-------------------------------
!scalars
 integer :: ii,jj,iq_ibz,itim_q,isym_q,iq_curr,isg1,isg2
!arrays
 integer, ABI_CONTIGUOUS pointer :: grottb(:)
 complex(gwpc),ABI_CONTIGUOUS pointer :: phsgt(:),bigomegatwsq(:,:),omegatw(:,:)
! *********************************************************************

 ! Save the index of the q-point for checking purpose.
 ppm%iq_bz = iq_bz

 ABI_MALLOC(botsq, (ppm%npwc, ppm%dm2_botsq))
 ABI_MALLOC(otq, (ppm%npwc, ppm%dm2_otq))
 ABI_MALLOC(eig, (ppm%dm_eig, ppm%dm_eig))

 ! Here there is a problem with the small q, still cannot use BZ methods
 iq_ibz = Qmesh%tab(iq_bz); isym_q = Qmesh%tabo(iq_bz); itim_q = (3-Qmesh%tabi(iq_bz))/2

 !call Qmesh%get_bz_item(iq_bz,qbz,iq_ibz,isym_q,itim_q,isirred=q_isirred)
 iq_curr = iq_ibz; if (ppm%mqmem == 0) iq_curr = 1

 grottb => Gsph%rottb (1:ppm%npwc, itim_q, isym_q)
 phsgt  => Gsph%phmSGt(1:ppm%npwc, isym_q)
 bigomegatwsq => ppm%bigomegatwsq(iq_curr)%vals
 omegatw      => ppm%omegatw(iq_curr)%vals

 ! Symmetrize the PPM parameters
 select case (ppm%model)
 case (PPM_GODBY_NEEDS, PPM_HYBERTSEN_LOUIE)
   ! Plasmon pole frequencies otq are invariant under symmetry
!$omp parallel do private(isg1, isg2)
   do jj=1,ppm%npwc
     isg2 = grottb(jj)
     do ii=1,ppm%npwc
       isg1 = grottb(ii)
       botsq(isg1,isg2) = bigomegatwsq(ii,jj)*phsgt(ii)*CONJG(phsgt(jj))
       otq  (isg1,isg2) = omegatw(ii,jj)
     end do
   end do

 case (PPM_LINDEN_HORSH)
   ! For notations see page 22 of Quasiparticle Calculations in solid (Aulbur et al)
   !  If q_bz = Sq_ibz + G0 then:
   !
   ! $\omega^2_{ii}(q_bz) = \omega^2_{ii}(q)$        (otq array)
   ! $\alpha_{ii}(q_bz)   = \alpha_{ii}(q)$          (botq array
   ! $\Phi_{SG-G0}(q_bz)  = \Phi_{G}(q) e^{-iSG.t}$  (eigenvectors of e^{-1}, eig array)
   !
   do ii=1,ppm%npwc ! DM bands index
     botsq(ii,1) = bigomegatwsq(ii,1)
     otq  (ii,1) = omegatw     (ii,1)
     do jj=1,ppm%npwc
       eig(grottb(jj),ii) = ppm%eigpot(iq_curr)%vals(jj,ii) * phsgt(jj)
     end do
   end do
   if (itim_q==2) eig=CONJG(eig) ! Time-reversal

 case (PPM_ENGEL_FARID)
   ! For notations see page 23 of Quasiparticle Calculations in solid (Aulbur et al.)
   ! If q_bz = Sq_ibz + G0 then:
   !
   ! $\omega^2_{ii}(q_bz) = \omega^2_{ii}(q)$        (otq array)
   ! $y_{SG-G0}(q_bz)     = y_{G}(q) e^{-iSG.t}$     (y=Lx)
   !
   do ii=1,ppm%npwc ! DM bands index
     otq(ii,1) = omegatw(ii,1)
     do jj=1,ppm%npwc
       botsq(grottb(jj),ii) = bigomegatwsq(jj,ii)*phsgt(jj)
     end do
   end do

 case default
   ABI_BUG(sjoin('Wrong ppm%model:',itoa(ppm%model)))
 end select

 ! Take into account time-reversal symmetry.
 if (itim_q == 2) then
!$omp parallel workshare
   botsq=CONJG(botsq)
!$omp end parallel workshare
 end if

end subroutine ppm_get_qbz
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/ppm_free
!! NAME
!!  ppm_free
!!
!! FUNCTION
!!  Deallocate all associated pointers defined in a variable of type ppmodel_t.
!!
!! SOURCE

subroutine ppm_free(ppm)

!Arguments ------------------------------------
 class(ppmodel_t),intent(inout) :: ppm

!Local variables-------------------------------
 integer :: dim_q,iq_ibz
! *********************************************************************

 ABI_SFREE(ppm%bigomegatwsq_qbz_vals)
 ABI_SFREE(ppm%omegatw_qbz_vals)
 ABI_SFREE(ppm%eigpot_qbz_vals)

 dim_q = ppm%nqibz; if (ppm%mqmem==0) dim_q=1

 if (allocated(ppm%bigomegatwsq)) then
   do iq_ibz=1,dim_q
     call ppm%bigomegatwsq(iq_ibz)%free()
   end do
   ABI_FREE(ppm%bigomegatwsq)
 end if
 if (allocated(ppm%omegatw)) then
   do iq_ibz=1,dim_q
     call ppm%omegatw(iq_ibz)%free()
   end do
   ABI_FREE(ppm%omegatw)
 end if
 if (allocated(ppm%eigpot)) then
   do iq_ibz=1,dim_q
     call ppm%eigpot(iq_ibz)%free()
   end do
   ABI_FREE(ppm%eigpot)
 end if

 ! logical flags must be deallocated here.
 ABI_SFREE(ppm%keep_qibz)
 ABI_SFREE(ppm%has_qibz)

end subroutine ppm_free
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/ppm_malloc_iqibz
!! NAME
!!  ppm_malloc_iqibz
!!
!! FUNCTION
!!  Allocate the ppmodel tables for the selected q-point in the IBZ.
!!
!! INPUT
!!  iq_ibz=Index of the q-point in the IBZ.
!!
!! SOURCE

subroutine ppm_malloc_iqibz(ppm, iq_ibz)

!Arguments ------------------------------------
 class(ppmodel_t),intent(inout) :: ppm
 integer,intent(in) :: iq_ibz

!Local variables-------------------------------
 integer :: ierr
! *********************************************************************

 ABI_CHECK(allocated(ppm%bigomegatwsq), "bigomegatwsq is not allocated")
 ABI_CHECK(allocated(ppm%omegatw), "omegatwsq is not allocated")
 ABI_CHECK(allocated(ppm%eigpot), "eigpot is not allocated")

 ABI_CHECK_IGEQ(size(ppm%bigomegatwsq), iq_ibz, "bigomegatwsq too small")
 ABI_CHECK_IGEQ(size(ppm%omegatw), iq_ibz, "omegatwsq too small")
 ABI_CHECK_IGEQ(size(ppm%eigpot), iq_ibz, "eigpot too small")

 ABI_MALLOC_OR_DIE(ppm%bigomegatwsq(iq_ibz)%vals, (ppm%npwc, ppm%dm2_botsq), ierr)
 ABI_MALLOC_OR_DIE(ppm%omegatw(iq_ibz)%vals, (ppm%npwc, ppm%dm2_otq), ierr)
 ABI_MALLOC_OR_DIE(ppm%eigpot(iq_ibz)%vals, (ppm%dm_eig, ppm%dm_eig), ierr)

 ppm%has_qibz(iq_ibz) = PPM_TAB_ALLOCATED

end subroutine ppm_malloc_iqibz
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/ppm_table_free_iqibz
!! NAME
!!  ppm_table_free_iqibz
!!
!! FUNCTION
!!  Free the ppmodel tables for the selected q-point in the IBZ.
!!
!! INPUT
!!  iq_ibz = Index of the q-point in the IBZ.
!!
!! SOURCE

subroutine ppm_table_free_iqibz(ppm, iq_ibz)

!Arguments ------------------------------------
 class(ppmodel_t),intent(inout) :: ppm
 integer,intent(in) :: iq_ibz
! *********************************************************************

 if (allocated(ppm%bigomegatwsq)) call ppm%bigomegatwsq(iq_ibz)%free()
 if (allocated(ppm%omegatw)) call ppm%omegatw(iq_ibz)%free()
 if (allocated(ppm%eigpot)) call ppm%eigpot(iq_ibz)%free()

 ppm%has_qibz(iq_ibz) = PPM_NOTAB

end subroutine ppm_table_free_iqibz
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/ppm_init
!! NAME
!!  ppm_init
!!
!! FUNCTION
!!  Initialize dimensions and other important variables related to the ppmodel
!!
!! INPUTS
!! ppmodel=
!! drude_plsmf=
!!
!! SOURCE

subroutine ppm_init(ppm, mqmem, nqibz, npwe, ppmodel, drude_plsmf, invalid_freq)

!Arguments ------------------------------------
 class(ppmodel_t),intent(out) :: ppm
 integer,intent(in) :: mqmem, nqibz, npwe, ppmodel, invalid_freq
 real(dp),intent(in) :: drude_plsmf

!Local variables-------------------------------
!scalars
 integer :: dim_q,iq_ibz
 logical :: ltest
 !character(len=500) :: msg
! *********************************************************************

 ppm%nqibz = nqibz; ppm%mqmem = mqmem; ppm%invalid_freq = invalid_freq
 !call wrtout(std_out, sjoin(' ppm%mqmem:', itoa(ppm%mqmem), 'ppm%nqibz:', itoa(ppm%nqibz)))
 ltest = (ppm%mqmem == 0 .or. ppm%mqmem == ppm%nqibz)
 ! ABI_CHECK(ltest,'Wrong value for mqmem')

 ppm%npwc        = npwe
 ppm%model       = ppmodel
 ppm%drude_plsmf = drude_plsmf
 ppm%userho      = 0
 if (any(ppmodel == [PPM_HYBERTSEN_LOUIE, PPM_LINDEN_HORSH, PPM_ENGEL_FARID])) ppm%userho = 1

 ABI_MALLOC(ppm%keep_qibz, (nqibz))
 ppm%keep_qibz = .FALSE.; if (ppm%mqmem > 0) ppm%keep_qibz = .TRUE.

 ABI_MALLOC(ppm%has_qibz, (nqibz))
 ppm%has_qibz = PPM_NOTAB

 ! Full q-mesh is stored or out-of-memory solution.
 dim_q = ppm%nqibz; if (ppm%mqmem == 0) dim_q=1

 ABI_MALLOC(ppm%bigomegatwsq, (dim_q))
 ABI_MALLOC(ppm%omegatw, (dim_q))
 ABI_MALLOC(ppm%eigpot, (dim_q))

 select case (ppm%model)
 case (PPM_NONE)
   ABI_WARNING("Called with ppmodel == 0")
   ppm%dm2_botsq = 0
   ppm%dm2_otq   = 0
   ppm%dm_eig    = 0
   RETURN

 case (PPM_GODBY_NEEDS, PPM_HYBERTSEN_LOUIE)
   ppm%dm2_botsq = ppm%npwc
   ppm%dm2_otq   = ppm%npwc
   ppm%dm_eig    = 1 ! Should be set to 0, but g95 does not like zero-sized arrays

 case (PPM_LINDEN_HORSH)
   ppm%dm2_botsq = 1
   ppm%dm2_otq   = 1
   ppm%dm_eig    = ppm%npwc

 case (PPM_ENGEL_FARID)
   ppm%dm2_botsq = ppm%npwc
   ppm%dm2_otq   = 1
   ppm%dm_eig    = 1 ! Should be set to 0, but g95 does not like zero-sized arrays

 case default
   ABI_BUG(sjoin('Wrong ppm%model:', itoa(ppm%model)))
 end select

 ! Allocate tables depending on the value of keep_qibz.
 do iq_ibz=1,dim_q
   call ppm%malloc_iqibz(iq_ibz)
 end do

 call pstat_proc%print(_PSTAT_ARGS_)

end subroutine ppm_init
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/ppm_setup
!! NAME
!! ppm_setup
!!
!! FUNCTION
!!  Initialize some values of several arrays of the ppm datastructure
!!  that are used in case of plasmonpole calculations
!!  This is a wrapper around different plasmonpole routines.
!!
!! INPUTS
!!  Cryst<crystal_t>=Info on the unit cell and crystal symmetries.
!!  Qmesh<kmesh_t>=the q-mesh used for the inverse dielectric matrix
!!    %nibz=number of irreducible q-points
!!    %ibz(3,%nibz)=the irred q-point
!!  npwe=number of G vectors for the correlation part
!!  nomega=number of frequencies in $\epsilon^{-1}$
!!  omega=frequencies in epsm1
!!  epsm1=the inverse dielctric matrix
!!  ngfftf(18)=contain all needed information about the 3D fine FFT mesh, see ~abinit/doc/variables/vargs.htm#ngfft
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!  nfftf=the number of points in the FFT mesh (for this processor)
!!  rhor_tot(nfftf)=the total charge in real space
!!
!! SIDE EFFECTS
!!  == if ppmodel 1 or 2 ==
!!   %omegatw and %bigomegatwsq
!!  == if ppmodel 3 ==
!!   %omegatw, %bigomegatwsq and %eigpot
!!  == if ppmodel 4 ==
!!   %omegatw and %bigomegatwsq
!!
!! NOTES
!! * FFT parallelism not implemented.
!! * TODO: rhor_tot should be replaced by rhog_tot
!!
!! SOURCE

subroutine ppm_setup(ppm, Cryst, Qmesh, npwe, nomega, omega, epsm1, nfftf, gvec, ngfftf, rhor_tot, &
                     iqiA) ! Optional

!Arguments ------------------------------------
!scalars
 class(ppmodel_t),intent(inout) :: ppm
 integer,intent(in) :: nfftf,npwe,nomega
 integer,intent(in),optional :: iqiA
 type(kmesh_t),intent(in) :: Qmesh
 type(crystal_t),intent(in) :: Cryst
!arrays
 integer,intent(in) :: gvec(3,npwe),ngfftf(18)
 real(dp),intent(in) :: rhor_tot(nfftf)
 complex(dpc),intent(in) :: omega(nomega)
 complex(gwpc),intent(in) :: epsm1(:,:,:,:)

!Local variables-------------------------------
!scalars
 integer :: nqiA,iq_ibz
 real(dp) :: n_at_G_zero
 logical :: single_q
 character(len=500) :: msg
!scalars
 real(dp) :: qpt(3)
! *************************************************************************

 !@ppmodel_t
 !
 ! === if iqiA is present, then consider only one qpoint to save memory ===
 ! * This means the object has been already initialized
 nqiA = Qmesh%nibz; single_q = .FALSE.
 if (PRESENT(iqiA)) then
   nqiA = 1; single_q = .TRUE.
 end if

 ! Allocate plasmonpole parameters
 ! TODO ppmodel==1 by default, should be set to 0 if AC and CD
 select case (ppm%model)

 case (PPM_NONE)
   ABI_COMMENT(' Skipping Plasmompole model calculation')

 case (PPM_GODBY_NEEDS)
   ! Note: the q-dependency enters only through epsilon^-1.
   do iq_ibz=1,nqiA
     call cppm1par(npwe,nomega,omega,ppm%drude_plsmf,&
                   epsm1(:,:,:,iq_ibz),ppm%omegatw(iq_ibz)%vals,ppm%bigomegatwsq(iq_ibz)%vals)
   end do

 case (PPM_HYBERTSEN_LOUIE)
   do iq_ibz=1,nqiA
     qpt = Qmesh%ibz(:,iq_ibz); if (single_q) qpt=Qmesh%ibz(:,iqiA)

     call cppm2par(qpt,npwe,epsm1(:,:,1,iq_ibz),ngfftf,gvec,Cryst%gprimd,rhor_tot,nfftf,Cryst%gmet,&
                   ppm%bigomegatwsq(iq_ibz)%vals,ppm%omegatw(iq_ibz)%vals,ppm%invalid_freq)
   end do

   ! Quick-and-dirty change of the plasma frequency. Never executed in standard runs.
   if (ppm%force_plsmf>tol6) then ! Integrate the real-space density
      n_at_G_zero = SUM(rhor_tot(:))/nfftf
      ! Change the prefactor
      write(msg,'(2(a,es16.8))') 'Forced ppmfreq:',ppm%force_plsmf*Ha_eV,' nelec/ucvol:',n_at_G_zero
      ABI_WARNING(msg)
      ppm%force_plsmf = (ppm%force_plsmf**2)/(four_pi*n_at_G_zero)
      do iq_ibz=1,ppm%nqibz
        ppm%bigomegatwsq(iq_ibz)%vals = ppm%force_plsmf * ppm%bigomegatwsq(iq_ibz)%vals
        ppm%omegatw(iq_ibz)%vals      = ppm%force_plsmf * ppm%omegatw(iq_ibz)%vals
      end do
      write(msg,'(a,es16.8)') 'Plasma frequency forced in HL ppmodel, new prefactor is:',ppm%force_plsmf
      ABI_WARNING(msg)
   end if

 case (PPM_LINDEN_HORSH) ! TODO Check better double precision, this routine is in a messy state
   do iq_ibz=1,nqiA
     qpt = Qmesh%ibz(:,iq_ibz); if (single_q) qpt=Qmesh%ibz(:,iqiA)
     call cppm3par(qpt,npwe,epsm1(:,:,1,iq_ibz),ngfftf,gvec,Cryst%gprimd,rhor_tot,nfftf,&
                   ppm%bigomegatwsq(iq_ibz)%vals,ppm%omegatw(iq_ibz)%vals(:,1),ppm%eigpot(iq_ibz)%vals)
   end do

 case (PPM_ENGEL_FARID)  ! TODO Check better double precision, this routine is in a messy state
   do iq_ibz=1,nqiA
     qpt = Qmesh%ibz(:,iq_ibz); if (single_q) qpt=Qmesh%ibz(:,iqiA)
     if ((ALL(ABS(qpt)<1.0e-3))) qpt = GW_Q0_DEFAULT ! FIXME
     call cppm4par(qpt,npwe,epsm1(:,:,1,iq_ibz),ngfftf,gvec,Cryst%gprimd,rhor_tot,nfftf,&
                   ppm%bigomegatwsq(iq_ibz)%vals,ppm%omegatw(iq_ibz)%vals(:,1))
   end do

 case default
   ABI_BUG(sjoin('Wrong ppm%model:',itoa(ppm%model)))
 end select

end subroutine ppm_setup
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/ppm_getem1
!! NAME
!!  ppm_getem1
!!
!! FUNCTION
!!  Calculate the symmetrized inverse dielectric matrix from the parameters of the plasmon-pole model.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine ppm_getem1(ppm, mpwc, iqibz, zcut, nomega, omega, Vcp, em1q, &
                      only_ig1, only_ig2) ! Optional

!Arguments ------------------------------------
!scalars
 class(ppmodel_t),intent(in) :: ppm
 integer,intent(in) :: mpwc,iqibz,nomega
 type(vcoul_t),intent(in) :: Vcp
 real(dp),intent(in) :: zcut
 integer,optional,intent(in) :: only_ig1,only_ig2
!arrays
 complex(dpc),intent(in) :: omega(nomega)
 complex(dpc),intent(out) :: em1q(mpwc,mpwc,nomega)

!Local variables-------------------------------
!scalars
 integer :: ig1,ig2,io,idm,ig1_min,ig2_min,ig1_max,ig2_max
 real(dp) :: den
 complex(dpc) :: qpg1,qpg2,ug1,ug2
 complex(dpc) :: delta,em1ggp,otw,zzpq,yg1,yg2,bot1,bot2,chig1g2
 !character(len=500) :: msg
! *************************************************************************

 ABI_CHECK(ppm%mqmem/=0,'mqmem==0 not implemented')

 !TODO zcut should be an entry in ppm
 delta=CMPLX(zero,zcut)

 ! To save memory, a particular combination of
 ! ig1 and ig2 can be selected
 ig1_min = 1
 ig2_min = 1
 ig1_max = ppm%npwc
 ig2_max = ppm%npwc
 if (present(only_ig1)) then
   ig1_min = only_ig1
   ig1_max = only_ig1
 end if
 if (present(only_ig2)) then
   ig2_min = only_ig2
   ig2_max = only_ig2
 end if

 select case (ppm%model)
 case (PPM_GODBY_NEEDS, PPM_HYBERTSEN_LOUIE)
   do io=1,nomega
     do ig2=ig2_min,ig2_max
       do ig1=ig1_min,ig1_max
        !den = omega(io)**2-REAL(ppm%omegatw(iqibz)%vals(ig1,ig2)**2)
        !if (den**2<zcut**2) den = omega(io)**2-REAL( (ppm%omegatw(iqibz)%vals(ig1,ig2)-delta)**2 )
        den = omega(io)**2 - REAL( (ppm%omegatw(iqibz)%vals(ig1,ig2)-delta)**2 )
        em1ggp = ppm%bigomegatwsq(iqibz)%vals(ig1,ig2)/den
        if (ig1==ig2) em1ggp=em1ggp+one
        em1q(ig1,ig2,io)=em1ggp
        !em1q(ig1,ig2,io)=em1ggp*Vcp%vc_sqrt(ig1,iqibz)*Vcp%vc_sqrt(ig2,iqibz)
       end do
     end do
     !
   end do !io

 case (PPM_LINDEN_HORSH)
   !TODO Check coefficients
   do io=1,nomega
     do ig2=ig2_min,ig2_max
       do ig1=ig1_min,ig1_max
         !
         em1ggp=czero
         do idm=1,ppm%npwc
           !den=omega(io)**2-(ppm%omegatw(iqibz)%vals(ig1,ig2)-delta)**2
           !em1w(io)=em1w(io)+eigvec(ig1,idm,iqibz)*conjg(eigvec(ig2,idm,iqibz))*bigomegatwsq(ig1,ig2,iqibz)/den
           ug1 = ppm%eigpot(iqibz)%vals(ig1,idm)
           ug2 = ppm%eigpot(iqibz)%vals(ig2,idm)
           otw = ppm%bigomegatwsq(iqibz)%vals(idm,1)*ppm%omegatw(iqibz)%vals(idm,1)
           zzpq=ppm%bigomegatwsq(iqibz)%vals(idm,1)
           den=half*REAL(zzpq*otw*( one/(omega(io)-otw+delta) - one/(omega(io)+otw-delta) ))
           em1ggp=em1ggp+ug1*CONJG(ug2)*den
           !eigenvalues(idm,io)=one + half*REAL(zzpq*otw*( one/(omega(io)-otw+delta) - one/(omega(io)+otw-delta) ))
         end do
         if (ig2==ig1) em1ggp=em1ggp+one
         em1q(ig1,ig2,io)=em1ggp
       end do !ig1
     end do !ig2
     !
   end do !iomega

 case (PPM_ENGEL_FARID)
   ! Make e^-1
   do io=1,nomega
     do ig2=ig2_min,ig2_max
       qpg2=one/Vcp%vc_sqrt(ig2,iqibz)
       do ig1=ig1_min,ig1_max
         qpg1=one/Vcp%vc_sqrt(ig1,iqibz)

         chig1g2=czero
         do idm=1,ppm%npwc
           otw =ppm%omegatw(iqibz)%vals(idm,1)
           bot1=ppm%bigomegatwsq(iqibz)%vals(ig1,idm)
           bot2=ppm%bigomegatwsq(iqibz)%vals(ig2,idm)
           yg1=SQRT(otw/four_pi)*qpg1*bot1
           yg2=SQRT(otw/four_pi)*qpg2*bot2
           chig1g2=chig1g2 + yg1*CONJG(yg2)/(omega(io)**2-(otw-delta)**2)
         end do

         em1ggp=four_pi*chig1g2/(qpg1*qpg2)
         if (ig1==ig2) em1ggp=em1ggp+one
         em1q(ig1,ig2,io)=em1ggp !*Vcp%vc_sqrt(ig1,iqibz)*Vcp%vc_sqrt(ig2,iqibz)
       end do !ig1
     end do !ig2
   end do !iomega

 case default
   ABI_BUG(sjoin('Wrong ppm%model:',itoa(ppm%model)))
 end select

end subroutine ppm_getem1
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/ppm_getem1_one_ggp
!! NAME
!!  ppm_getem1_one_ggp
!!
!! FUNCTION
!!  Same as ppm_getem1, but does it for a single set of G,G' vectors
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine ppm_getem1_one_ggp(ppm, iqibz, zcut, nomega, omega, Vcp, em1q, ig1, ig2)

!Arguments ------------------------------------
!scalars
 class(ppmodel_t),intent(in) :: ppm
 integer,intent(in) :: iqibz,nomega
 type(vcoul_t),intent(in) :: Vcp
 real(dp),intent(in) :: zcut
 integer, intent(in) :: ig1,ig2
!arrays
 complex(dpc),intent(in) :: omega(nomega)
 complex(dpc),intent(out) :: em1q(nomega)

!Local variables-------------------------------
!scalars
 integer :: io,idm !,ig1_min,ig2_min,ig2_max
 real(dp) :: den
 complex(dpc) :: qpg1,qpg2,ug1,ug2
 complex(dpc) :: delta,em1ggp,otw,zzpq,yg1,yg2,bot1,bot2,chig1g2
 !character(len=500) :: msg
! *************************************************************************

 ABI_CHECK(ppm%mqmem /= 0, 'mqmem==0 not implemented')

 !TODO zcut should be an entry in ppm
 delta=CMPLX(zero,zcut)

 select case (ppm%model)

 case (PPM_GODBY_NEEDS, PPM_HYBERTSEN_LOUIE)
   do io=1,nomega
     !den = omega(io)**2-REAL(ppm%omegatw(iqibz)%vals(ig1,ig2)**2)
     !if (den**2<zcut**2) den = omega(io)**2-REAL( (ppm%omegatw(iqibz)%value(ig1,ig2)-delta)**2 )
     den = omega(io)**2-REAL( (ppm%omegatw(iqibz)%vals(ig1,ig2)-delta)**2 )
     em1ggp = ppm%bigomegatwsq(iqibz)%vals(ig1,ig2)/den
     if (ig1==ig2) em1ggp=em1ggp+one
     em1q(io)=em1ggp
     !em1q(io)=em1ggp*Vcp%vc_sqrt(ig1,iqibz)*Vcp%vc_sqrt(ig2,iqibz)
   end do !io

 case (PPM_LINDEN_HORSH)
   !TODO Check coefficients
   do io=1,nomega
     em1ggp=czero
     do idm=1,ppm%npwc
       !den=omega(io)**2-(ppm%omegatw(iqibz)%vals(ig1,ig2)-delta)**2
       !em1w(io)=em1w(io)+eigvec(ig1,idm,iqibz)*conjg(eigvec(ig2,idm,iqibz))*bigomegatwsq(ig1,ig2,iqibz)/den
       ug1 =ppm%eigpot(iqibz)%vals(ig1,idm)
       ug2 =ppm%eigpot(iqibz)%vals(ig2,idm)
       otw =ppm%bigomegatwsq(iqibz)%vals(idm,1)*ppm%omegatw(iqibz)%vals(idm,1)
       zzpq=ppm%bigomegatwsq(iqibz)%vals(idm,1)
       den=half*REAL(zzpq*otw*( one/(omega(io)-otw+delta) - one/(omega(io)+otw-delta) ))
       em1ggp=em1ggp+ug1*CONJG(ug2)*den
       !eigenvalues(idm,io)=one + half*REAL(zzpq*otw*( one/(omega(io)-otw+delta) - one/(omega(io)+otw-delta) ))
     end do

     if (ig2==ig1) em1ggp=em1ggp+one
     em1q(io)=em1ggp
   end do !iomega

 case (PPM_ENGEL_FARID)
   ! Make e^-1
   do io=1,nomega
     qpg2=one/Vcp%vc_sqrt(ig2,iqibz)
     qpg1=one/Vcp%vc_sqrt(ig1,iqibz)

     chig1g2=czero
     do idm=1,ppm%npwc
       otw =ppm%omegatw(iqibz)%vals(idm,1)
       bot1=ppm%bigomegatwsq(iqibz)%vals(ig1,idm)
       bot2=ppm%bigomegatwsq(iqibz)%vals(ig2,idm)
       yg1=SQRT(otw/four_pi)*qpg1*bot1
       yg2=SQRT(otw/four_pi)*qpg2*bot2
       chig1g2=chig1g2 + yg1*CONJG(yg2)/(omega(io)**2-(otw-delta)**2)
     end do

     em1ggp=four_pi*chig1g2/(qpg1*qpg2)
     if (ig1==ig2) em1ggp=em1ggp+one
     em1q(io)=em1ggp !*Vcp%vc_sqrt(ig1,iqibz)*Vcp%vc_sqrt(ig2,iqibz)
   end do ! io

 case default
   ABI_BUG(sjoin('Wrong ppm%model:',itoa(ppm%model)))
 end select

end subroutine ppm_getem1_one_ggp
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/ppm_get_eigenvalues
!! NAME
!!  ppm_get_eigenvalues
!!
!! FUNCTION
!!  Constructs the inverse dielectri matrixc starting from the plasmon-pole
!!  parameters and calculates the frequency-dependent eigenvalues for each
!!  of the nomega frequencies specified in the array omega.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine ppm_get_eigenvalues(ppm, iqibz, zcut, nomega, omega, Vcp, eigenvalues)

!Arguments ------------------------------------
!scalars
 class(ppmodel_t),intent(in) :: ppm
 integer,intent(in) :: iqibz,nomega
 type(vcoul_t),intent(in) :: Vcp
 real(dp),intent(in) :: zcut
!arrays
 complex(dpc),intent(in) :: omega(nomega)
 complex(dpc),intent(out) :: eigenvalues(ppm%npwc,nomega)

!Local variables-------------------------------
!scalars
 integer :: info,lwork,negw,ig1,ig2,idx,sdim,iomega,ierr
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: ww(:),rwork(:)
 complex(dpc),allocatable :: work(:),Adpp(:),eigvec(:,:),wwc(:),vs(:,:),Afull(:,:)
 complex(dpc),allocatable :: em1q(:,:,:)
 logical,allocatable :: bwork(:)
 logical :: sortcplx !BUG in abilint
! *************************************************************************

 ABI_CHECK(ppm%mqmem/=0,'mqmem==0 not implemented')

 ABI_MALLOC(em1q, (ppm%npwc,ppm%npwc,nomega))

 call ppm%getem1(ppm%npwc,iqibz,zcut,nomega,omega,Vcp,em1q)

 do iomega=1,nomega
   if (ABS(REAL(omega(iomega)))>0.00001) then
     ! Eigenvalues for a generic complex matrix.

     lwork=4*2*ppm%npwc
     ABI_MALLOC(wwc,(ppm%npwc))
     ABI_MALLOC(work,(lwork))
     ABI_MALLOC(rwork,(ppm%npwc))
     ABI_MALLOC(bwork,(ppm%npwc))
     ABI_MALLOC(vs,(ppm%npwc,ppm%npwc))
     ABI_MALLOC(Afull,(ppm%npwc,ppm%npwc))
     Afull=em1q(:,:,iomega)

     !for the time being, no sorting. Maybe here I should sort using the real part?
     call ZGEES('V','N',sortcplx,ppm%npwc,Afull,ppm%npwc,sdim,wwc,vs,ppm%npwc,work,lwork,rwork,bwork,info)
     if (info/=0) then
      write(msg,'(2a,i10)')' ppm_get_eigenvalues: Error in ZGEES, diagonalizing complex matrix, info = ',info
      call wrtout(std_out,msg)
     end if

     eigenvalues(:,iomega)=wwc(:)

     ABI_FREE(wwc)
     ABI_FREE(work)
     ABI_FREE(rwork)
     ABI_FREE(bwork)
     ABI_FREE(vs)
     ABI_FREE(Afull)

   else
     ! === Hermitian Case ===
     lwork=2*ppm%npwc-1
     ABI_MALLOC(ww,(ppm%npwc))
     ABI_MALLOC(work,(lwork))
     ABI_MALLOC(rwork,(3*ppm%npwc-2))
     ABI_MALLOC(eigvec,(ppm%npwc,ppm%npwc))

     ABI_MALLOC_OR_DIE(Adpp,(ppm%npwc*(ppm%npwc+1)/2), ierr)
     !write(std_out,*) 'in hermitian'

     idx=0
     do ig2=1,ppm%npwc
       do ig1=1,ig2
         idx=idx+1
         Adpp(idx)=em1q(ig1,ig2,iomega)
       end do
     end do

     ! Require eigenvectors as well
     call ZHPEV('V','U',ppm%npwc,Adpp,ww,eigvec,ppm%npwc,work,rwork,info)

     ABI_CHECK(info == 0, sjoin('Error diagonalizing matrix, info: ', itoa(info)))

     negw = (COUNT((REAL(ww)<tol6)))
     if (negw /= 0) then
       write(msg,'(a,i0,a,i0,a,f8.4)')'Found negative eigenvalues. No. ',negw,' at iqibz= ',iqibz,' minval= ',MINVAL(REAL(ww))
        ABI_WARNING(msg)
     end if

     eigenvalues(:,iomega)=ww(:)

     ABI_FREE(ww)
     ABI_FREE(work)
     ABI_FREE(rwork)
     ABI_FREE(eigvec)
     ABI_FREE(Adpp)
   end if
 end do !iomega

 ABI_FREE(em1q)

end subroutine ppm_get_eigenvalues
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/cppm1par
!! NAME
!! cppm1par
!!
!! FUNCTION
!! Calculate the plasmon-pole parameters big-omega-twiddle-squared and omega-twiddle from
!! epsilon-twiddle^-1 calculated for nomega (usually 2) frequencies omega=0 and omega=iE0.
!!
!! INPUTS
!!  epsm1(npwc,npwc,nomega)=dielectric matrix at nomega frequencies.
!!  npwc=number of plane waves
!!  nomega=number of frequencies (usually 2)
!!  omega(nomega)=frequencies
!!  omegaplasma=input variable or Drude plasma frequency
!!
!! OUTPUT
!!  bigomegatwsq(npwc,npwc)=parameter of the plasmon-pole model (see gwa.pdf file)
!!  omegatw(npwc,npwc)=parameter of the plasmon-pole model (see gwa.pdf file)
!!
!! TODO
!!  Calculation can be done in place.
!!
!! SOURCE

subroutine cppm1par(npwc, nomega, omega, omegaplasma, epsm1, omegatw, bigomegatwsq)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nomega,npwc
 real(dp),intent(in) :: omegaplasma
!arrays
 complex(dpc),intent(in) :: omega(nomega)
 complex(gwpc),intent(in) :: epsm1(npwc,npwc,nomega)
 complex(gwpc),intent(out) :: omegatw(npwc,npwc), bigomegatwsq(npwc,npwc)

!Local variables-------------------------------
!scalars
 integer :: ig,igp,io,io0,ioe0
 real(dp) :: e0,minomega
 character(len=500) :: msg
 complex(gwpc) :: AA,omegatwsq,diff,ratio,epsm1_io0,epsm1_ioe0
! *************************************************************************

 ! Find omega=0 and omega=imag (closest to omegaplasma) to fit the ppm parameters
 minomega=1.0d-3; io0=0
 do io=1,nomega
   if (ABS(omega(io))<minomega) then
     io0=io; minomega=ABS(omega(io))
   end if
 end do
 ABI_CHECK(io0 /= 0, "omega=0 not found")

 minomega=1.0d-3; e0=200.0; ioe0=0
 do io=1,nomega
   if (REAL(omega(io))<minomega.and.AIMAG(omega(io))>minomega) then
     if (ABS(AIMAG(omega(io))-omegaplasma)<ABS(e0-omegaplasma)) then
       ioe0=io; e0=AIMAG(omega(io))
     end if
   end if
 end do

 write(msg,'(a,f9.4,a)')' Imaginary frequency for fit located at: ',e0*Ha_eV,' [eV] '
 call wrtout(std_out, msg)
 ABI_CHECK(ioe0 /= 0,"Imaginary omega not found")

 ! ================================================================
 ! === Calculate plasmon-pole A parameter A=epsilon^-1(0)-delta ===
 ! ================================================================
 do ig=1,npwc
   do igp=1,npwc
     epsm1_io0  = epsm1(ig,igp,io0)
     epsm1_ioe0 = epsm1(ig,igp,ioe0)

     AA=epsm1_io0
     if (ig==igp) AA=AA-one

     ! === Calculate plasmon-pole omega-twiddle-square parameter ===
     ! XG201009 Strangely, the next formula does not work with gcc43-debug
     ! omegatwsq=(AA/(epsm1_io0-epsm1_ioe0)-one)*e0**2
     ! This seems to be due to precision issue at the level of division by a complex whose norm squared
     ! is below the smallest representable number.
     ! After many trials, I have decided to shift the difference by a small number ... well, not so small ...
     ! for numerical issues
     diff=epsm1_io0-epsm1_ioe0
     diff=diff+cmplx(tol10,tol10)
     ratio=AA/diff
     omegatwsq=(ratio-cone)*e0**2
     !
     ! If omega-twiddle-squared is negative,set omega-twiddle-squared to 1.0 (a reasonable way of treating
     ! such terms, in which epsilon**-1 was originally increasing along this part of the imaginary axis)
     ! (note: originally these terms were ignored in Sigma; this was changed on 6 March 1990.)

     if (REAL(omegatwsq)<=zero) omegatwsq=one
     !
     ! Get omega-twiddle. Neglect the imag part (if any) in omega-twiddle-squared
     omegatw(ig,igp)=SQRT(REAL(omegatwsq))

     ! Get big-omega-twiddle-squared=-omega-twiddle-squared AA
     bigomegatwsq(ig,igp)=-AA*omegatw(ig,igp)**2

   end do !igp
 end do !ig

 write(msg,'(2a,f15.12,2a,2i5,a)')ch10,&
   ' cppm1par : omega twiddle minval [eV]  = ',MINVAL(ABS(omegatw))*Ha_eV,ch10,&
   '            omega twiddle min location = ',MINLOC(ABS(omegatw)),ch10
 call wrtout(std_out,msg)

end subroutine cppm1par
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/cppm2par
!! NAME
!! cppm2par
!!
!! FUNCTION
!!  Calculate plasmon-pole parameters of the Hybertsen and Louie model (PRB 34, 5390 (1986) [[cite:Hybertsen1986]])
!!
!! INPUTS
!!  qpt(3)=The coordinates of the q-point in the IBZ.
!!  epsm1(npwc,npwc)=symmetrized inverse dielectric (static limit is used)
!!  gmet(3,3)=metric in reciprocal space
!!  ngfftf(18)=contain all needed information about the 3D fine FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  npwc=number of plane waves in epsm1
!!  rhor(nfftf)=charge density on the real space FFT grid
!!  nfftf= total number of points in the fine FFT mesh  (for this processor)
!!  invalid_freq: what to do when PPM omega is found negative or imaginary
!!     0) drop it (default as specified in Hybersen-Louie original GW paper)
!!     1) set to 1 hartree
!!     2) set to infinity
!!
!! OUTPUT
!!  bigomegatwsq(npwc,npwc)= squared bare plasma frequencies
!!   \Omega^2_{G1 G2}(q) = 4\pi \frac {(q+G1).(q+G2)}/{|q+G1|^2} n(G1-G2)
!!  omegatw(npwc,npwc)= plasmon frequencies \tilde\omega_{G1 G2}(q) where:
!!  \tilde\omega^2_{G1 G2}(q) =
!!    \frac {\Omega^2_{G1 G2}(q)} {\delta_{G1 G2}-\tilde\epsilon^{-1}_{G1 G2} (q, \omega=0)}
!!
!! SOURCE

subroutine cppm2par(qpt, npwc, epsm1, ngfftf, gvec, gprimd, rhor, nfftf, gmet, bigomegatwsq, omegatw, invalid_freq)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwc,nfftf,invalid_freq
!arrays
 integer,intent(in) :: gvec(3,npwc), ngfftf(18)
 real(dp),intent(in) :: qpt(3),gmet(3,3),gprimd(3,3), rhor(nfftf)
 complex(gwpc),intent(in) :: epsm1(npwc,npwc)
 complex(gwpc),intent(out) :: bigomegatwsq(npwc,npwc), omegatw(npwc,npwc)

!Local variables-------------------------------
!scalars
 integer :: ig,igp,nimwp,ngfft1,ngfft2,ngfft3,gmgp_idx,ierr
 real(dp) :: lambda,phi,AA
 logical,parameter :: use_symmetrized=.TRUE., check_imppf=.FALSE.
 character(len=500) :: msg
 type(MPI_type) :: MPI_enreg_seq
!arrays
 real(dp) :: qlist(3,1)
 real(dp),allocatable :: tmp_rhor(:),qratio(:,:),qplusg(:),rhog_dp(:,:)
 complex(gwpc),allocatable :: omegatwsq(:,:)
 complex(gwpc),allocatable :: rhog(:),rhogg(:,:),temp(:,:)  !MG these should be double precision TODO
!*************************************************************************

 call initmpi_seq(MPI_enreg_seq)
 call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfftf(2),ngfftf(3),'all')

 ! Calculate qratio(npwec,npvec) = (q+G).(q+Gp)/|q+G|^2 ===
 ABI_MALLOC_OR_DIE(qratio,(npwc,npwc), ierr)

 call cqratio(npwc,gvec,qpt,gmet,gprimd,qratio)
 !
 ! Compute the density in G space rhor(R)--> rhog(G)
 ABI_MALLOC(rhog_dp,(2,nfftf))
 ABI_MALLOC(rhog,(nfftf))
 ngfft1=ngfftf(1); ngfft2=ngfftf(2); ngfft3=ngfftf(3)

 ABI_MALLOC(tmp_rhor,(nfftf))
 tmp_rhor = rhor ! To avoid having to use intent(inout).
 call fourdp(1,rhog_dp,tmp_rhor,-1,MPI_enreg_seq,nfftf,1,ngfftf,0)
 ABI_FREE(tmp_rhor)

 rhog(1:nfftf)=CMPLX(rhog_dp(1,1:nfftf),rhog_dp(2,1:nfftf))

 ! Calculate the FFT index of each (G-Gp) vector and assign
 ! the value of the correspondent density simultaneously
 ABI_MALLOC_OR_DIE(rhogg,(npwc, npwc), ierr)

 ierr=0
 do ig=1,npwc
   do igp=1,npwc
     gmgp_idx = g2ifft(gvec(:,ig)-gvec(:,igp),ngfftf)
     if (gmgp_idx/=0) then
       rhogg(ig,igp)=rhog(gmgp_idx)
     else
       ierr=ierr+1
       rhogg(ig,igp)=czero
     end if
   end do
 end do

 if (ierr /= 0) then
   write(msg,'(a,i0,1x,3a)')&
    'Found ',ierr,' G1-G2 vectors falling outside the FFT box. ',ch10,&
    'Enlarge the FFT mesh to get rid of this problem. '
   ABI_WARNING(msg)
 end if

 rhogg=four_pi*rhogg
 ABI_FREE(rhog_dp)
 ABI_FREE(rhog)

 ! Calculate GPP parameters
 ! unsymmetrized epsm1 -> epsm1=|q+Gp|/|q+G|*epsm1
 ABI_MALLOC(qplusg,(npwc))
 ABI_MALLOC(temp,(npwc,npwc))
 ABI_MALLOC_OR_DIE(omegatwsq,(npwc,npwc), ierr)

 temp = -epsm1(:,:)
 !
 ! RS still not obvious for me whether one shall use the symmetrized inverse DM or the unsymmetrized one
 ! the default here is to use the symmetrized one, I must discuss this with XG
 !
 ! MG it turns out that using the symmetrized inverse DM in the plasmon-pole
 ! equations give the same results for the squared plasmon frequencies omegatwsq while the
 ! squared bare plasma frequencies bigomegatwsq related to the symmetrized dielectric matrix
 ! are obtained multiplying by |q+G1|/|q+G2|.
 !
 if (.not.use_symmetrized) then
   qlist(:,1) = qpt
   call cmod_qpg(1,1,qlist,npwc,gvec,gprimd,qplusg) !MG TODO here take care of small q
   do ig=1,npwc
     do igp=1,npwc
       temp(ig,igp)=qplusg(igp)/qplusg(ig)*temp(ig,igp)
     end do
   end do
 end if

 nimwp=0
 do ig=1,npwc
   temp(ig,ig)=temp(ig,ig)+one
   do igp=1,npwc
     bigomegatwsq(ig,igp) = rhogg(ig,igp)*qratio(ig,igp)
     omegatwsq(ig,igp)=bigomegatwsq(ig,igp)/temp(ig,igp)
     !
     ! Set to an arbitrary value the omegawsq which become negative or imaginary
     ! in principle these correspond to cases where the imaginary part of epsm1 does not have
     ! a well defined peak. The imaginary part of epsm1 in these cases oscillates  with a small amplitude
     ! since the amplitude A_GGpr=-pi/2*bigomegatwsq/omegatw,
     ! it follows that bigomegatwsq shall be set to zero for these cases
     if ( REAL(omegatwsq(ig,igp))<= tol12 .or. AIMAG(omegatwsq(ig,igp))**2*tol12> REAL(omegatwsq(ig,igp))**2) then
       nimwp=nimwp+1

       if ( invalid_freq == 1 ) then
        ! set omegatwsq to 1 hartree
         omegatwsq(ig,igp)=cone
         AA = epsm1(ig,igp)
         if ( ig == igp ) AA = AA - one
         omegatw(ig,igp)=SQRT(REAL(omegatwsq(ig,igp)))
         bigomegatwsq(ig,igp)=-AA*omegatw(ig,igp)**2
       elseif ( invalid_freq == 2 ) then
         ! set omegatwsq to infinity
         omegatwsq(ig,igp)=cone/tol6
         AA = epsm1(ig,igp)
         if ( ig == igp ) AA = AA - one
         omegatw(ig,igp)=SQRT(REAL(omegatwsq(ig,igp)))
         bigomegatwsq(ig,igp)=-AA*omegatw(ig,igp)**2
       else
         ! simply ignore all cases of omegatw with imaginary values
         bigomegatwsq(ig,igp)=(0.,0.)
         omegatw(ig,igp)=(ten,0.)
       end if
       if (check_imppf) then
         write(msg,'(a,2(i0,1x))')' Imaginary plasmon frequency at : ',ig,igp
         call wrtout(std_out,msg)
       end if
     else
       ! this part has been added to deal with systems without inversion symmetry
       ! this new implementation gives the same results as the previous one if
       ! omegatwsq is a pure real number and has the advantage of being an improved
       ! approach for systems without an inversion center.
       lambda=ABS(omegatwsq(ig,igp))
       phi=ATAN(AIMAG(omegatwsq(ig,igp))/REAL(omegatwsq(ig,igp)))
       omegatw(ig,igp)=SQRT(lambda/COS(phi))
       bigomegatwsq(ig,igp)=bigomegatwsq(ig,igp)*(1.-(0.,1.)*TAN(phi))
       ! Uncomment the following line and comment the previous to restore the old version.
       !omegatw(ig,igp)=sqrt(real(omegatwsq(ig,igp)))
     end if
   end do
 end do

 write(msg,'(3a,i0,a,i0)')' At q-point : ',trim(ktoa(qpt)), ' # imaginary plasmonpole frequencies: ',nimwp,' / ',npwc**2
 call wrtout(std_out, msg)
 write(msg,'(a,f12.8,a,3(i0,1x))') &
  " omega twiddle minval: ", MINVAL(ABS(omegatw))*Ha_eV, "[eV], min location: ",MINLOC(ABS(omegatw))
 call wrtout(std_out, msg)

 call destroy_mpi_enreg(MPI_enreg_seq)

 ABI_FREE(omegatwsq)
 ABI_FREE(rhogg)
 ABI_FREE(temp)
 ABI_FREE(qplusg)
 ABI_FREE(qratio)

end subroutine cppm2par
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/cppm3par
!! NAME
!! cppm3par
!!
!! FUNCTION
!! Calculate the plasmon-pole parameters using the von Linden-Horsh model (PRB 37, 8351, 1988) [[cite:vonderLinden1988]]
!! (see also Pag 22 of Quasiparticle Calculations in Solids [[cite:Aulbur2001]].
!!
!! INPUTS
!! epsm1(npwc,npwc))= symmetrized inverse dielectric
!! ngfftf(18)=contain all needed information about 3D fine FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! npwc=number of plane waves in epsm1
!! qratio=(q+G1).(q+G2)/(|q+G1|.|q+G2|)
!! rhor(nfftf)=charge density on the real space FFT grid
!! nfftf=number of points in the FFT grid (for this processor)
!! gvec(3,npwc)= G vectors in reduced coordinates
!!
!! OUTPUT
!!  omegatw(npwc,npwc)= plasmon pole positions
!!  bigomegatwsq(npwc,npwc)=(E_{q,ii}^{-1}-1)*omegatw
!!   where E^{-1} is the eigenvalue of the inverse dielectric matrix
!!  eigtot(npwc,npwc)=the eigvectors of the symmetrized inverse dielectric matrix
!!   (first index for G, second index for bands)
!!
!! SOURCE

subroutine cppm3par(qpt,npwc,epsm1,ngfftf,gvec,gprimd,rhor,nfftf,bigomegatwsq,omegatw,eigtot)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftf,npwc
!arrays
 integer,intent(in) :: gvec(3,npwc),ngfftf(18)
 real(dp),intent(in) :: qpt(3),gprimd(3,3),rhor(nfftf)
 complex(gwpc),intent(in) :: epsm1(npwc,npwc)
 complex(gwpc),intent(out) :: bigomegatwsq(npwc,1),omegatw(npwc) ,eigtot(npwc,npwc)

!Local variables-------------------------------
!TODO these should be dp
!scalars
 integer :: idx,ierr,ig,igp,ii,jj,ngfft1,ngfft2,ngfft3,gmgp_idx
 real(dp) :: num,qpg_dot_qpgp
 complex(dpc) :: conjg_eig
 logical :: qiszero
 character(len=500) :: msg
 type(MPI_type) :: MPI_enreg_seq
!arrays
 real(dp) :: b1(3),b2(3),b3(3),gppq(3),gpq(3),qlist(3,1)
 real(dp),allocatable :: eigval(:),qplusg(:),rhog_dp(:,:),zhpev2(:),tmp_rhor(:)
 complex(dpc),allocatable :: eigvec(:,:),matr(:),mm(:,:),rhog(:),rhogg(:,:), zhpev1(:),zz(:)
!*************************************************************************

 ! Fake MPI_type for the sequential part.
 call initmpi_seq(MPI_enreg_seq)
 call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfftf(2),ngfftf(3),'all')

 qiszero = (ALL(ABS(qpt)<1.0e-3))

 b1 = two_pi*gprimd(:,1); b2 = two_pi*gprimd(:,2); b3 = two_pi*gprimd(:,3)

 ngfft1=ngfftf(1); ngfft2=ngfftf(2); ngfft3=ngfftf(3)

 ABI_MALLOC(rhog_dp,(2,nfftf))
 ABI_MALLOC(rhog,(nfftf))
 ABI_MALLOC_OR_DIE(rhogg,(npwc, npwc), ierr)
 !
 ! === Compute the density in G space rhog(r)--> rho(G) ===
 ! FIXME this has to be fixed, rho(G) should be passed instead of doing FFT for each q

 ABI_MALLOC(tmp_rhor,(nfftf))
 tmp_rhor=rhor ! To avoid having to use intent(inout).
 call fourdp(1,rhog_dp,tmp_rhor,-1,MPI_enreg_seq,nfftf,1,ngfftf,0)
 ABI_FREE(tmp_rhor)

 rhog(1:nfftf)=CMPLX(rhog_dp(1,1:nfftf),rhog_dp(2,1:nfftf))
 !
 ! Calculate the FFT index of each (G-Gp) vector and assign the value
 ! of the correspondent density simultaneously
 ierr=0
 do ig=1,npwc
   do igp=1,npwc
     gmgp_idx = g2ifft(gvec(:,ig)-gvec(:,igp),ngfftf)
     if (gmgp_idx/=0) then
       rhogg(ig,igp)=rhog(gmgp_idx)
     else
       ierr=ierr+1
       rhogg(ig,igp)=czero
     end if
   end do
 end do

 if (ierr /= 0) then
   write(msg,'(a,i0,3a)')&
   'Found ',ierr,' G1-G2 vectors falling outside the FFT box. ',ch10,&
   'Enlarge the FFT mesh to get rid of this problem. '
   ABI_WARNING(msg)
 end if

 ! mm(G,Gp) = (q+G) \cdot (q+Gp) n(G-Gp)
 ABI_MALLOC_OR_DIE(mm, (npwc,npwc), ierr)

 do ig=1,npwc
   if (qiszero) then
     ! To be discussed with Riad, here we should use the small q
     ! to be consistent and consider the limit q-->0
     gpq(:)=gvec(:,ig)
   else
     gpq(:)=gvec(:,ig)+qpt
   end if
   do igp=1,npwc
     if (qiszero) then
       gppq(:)=gvec(:,igp)
     else
       gppq(:)=gvec(:,igp)+qpt
     end if
     qpg_dot_qpgp=zero
     do ii=1,3
       qpg_dot_qpgp=qpg_dot_qpgp+&
        ( gpq(1)*b1(ii) +gpq(2)*b2(ii) +gpq(3)*b3(ii))*&
        (gppq(1)*b1(ii)+gppq(2)*b2(ii)+gppq(3)*b3(ii))
     end do
     mm(ig,igp)=rhogg(ig,igp)*qpg_dot_qpgp
   end do !igp
 end do !ig

 ABI_FREE(rhog_dp)
 ABI_FREE(rhog)
 ! === Now we have rhogg,rho0 ===
 !
 ! Calculate the dielectric matrix eigenvalues and vectors
 ! Use only the static epsm1 i.e., only the w=0 part (eps(:,:,1,:))
 ABI_MALLOC(eigval,(npwc))
 ABI_MALLOC_OR_DIE(eigvec, (npwc, npwc), ierr)

 ABI_MALLOC(zz,(npwc))
 zz=czero

 ABI_MALLOC(qplusg,(npwc))

 ! Store the susceptibility matrix in upper mode before calling zhpev.
 ABI_MALLOC_OR_DIE(matr,(npwc*(npwc+1)/2), ierr)

 idx=1
 do ii=1,npwc
   do jj=1,ii
     matr(idx)=epsm1(jj,ii); idx=idx+1
   end do
 end do

 ABI_MALLOC(zhpev2,(3*npwc-2))
 ABI_MALLOC(zhpev1,(2*npwc-1))

 call ZHPEV('V','U',npwc,matr,eigval,eigvec,npwc,zhpev1,zhpev2,ierr)
 ABI_FREE(matr)
 ABI_FREE(zhpev2)
 ABI_FREE(zhpev1)

 if (ierr < 0) then
   write (msg,'(2a,i0,a)')&
    ' Failed to calculate the eigenvalues and eigenvectors of the dielectric matrix ',ch10,&
    ierr*(-1),'-th argument in the matrix has an illegal value. '
   ABI_ERROR(msg)
 end if

 if (ierr > 0) then
   write(msg,'(3a,i0,2a)')&
    ' Failed to calculate the eigenvalues and eigenvectors of the dielectric matrix ',ch10,&
    ' the algorithm failed to converge; ierr = ', ierr,ch10,&
    ' off-diagonal elements of an intermediate tridiagonal form did not converge to zero. '
   ABI_ERROR(msg)
 end if

 ! Calculate the PPM parameters and the eigenpotentials needed for
 ! the calculation of the generalized overlap matrix
 ! Note: the eigenpotentials has to be calculated on the FFT (G-Gp) index
 !
 ! Save eigenvectors of \tilde\epsilon^{-1}
 ! MG well it is better to save \Theta otherwise
 ! we have to calculare \Theta for each band, spin, k-point but oh well
 eigtot=eigvec

 qlist(:,1) = qpt
 call cmod_qpg(1,1,qlist,npwc,gvec,gprimd,qplusg) !MG TODO here take care of small q
 !
 ! Basic Equation:
 !
 ! \Theta_{q,ii}(G)=\Psi_{q,ii}(G)/|q+G|
 ! where \Psi_{q,ii}(G) is the eigenvector of \tilde\epsilon^{-1}

 ! \tilde\omega_{ii,q}^2= 4\pi (1-eigenval(ii,q)))
 ! \sum_{G,Gp} \Theta^*_{q,ii}(G) (q+G)\cdot(q+Gp) n(G-Gp) \Theta_{q,ii}(Gp)

 do ii=1,npwc !DM band
   ! Calculate \Theta_{q,ii}(G)
   ! why the first element is not modified? if the problem is the small value of qplusg(1)
   ! we could multiply by sqrt(mod((q+G)(q+G'))) and then add the sing at the end
   if (qiszero)then
     eigvec(2:,ii)=eigvec(2:,ii)/qplusg(2:)
   else
     eigvec(:,ii)=eigvec(:,ii)/qplusg(:)
   end if
   do ig=1,npwc
     conjg_eig=CONJG(eigvec(ig,ii))
     do igp=1,npwc
       if(qiszero .and. ig==1 .and. igp==1)then
         zz(ii)=zz(ii)+conjg_eig*rhogg(ig,igp)*eigvec(igp,ii)
       else
         zz(ii)=zz(ii)+conjg_eig*mm(ig,igp)*eigvec(igp,ii)
       end if
     end do
   end do

   num=one-eigval(ii)
   if (num<=zero) then
     ! here I think we should set bigomegatwsq=0 and omegatw to an arbitrary value
     ! maybe we can output a warning TO BE discussed with Riad
     if (ABS(num)<1.0d-4) then
       num=1.0d-5
     else
       ABI_ERROR("One or more imaginary plasmon pole energies")
     end if
   end if

   omegatw(ii)=SQRT(4*pi*REAL(zz(ii))/num)
   ! this should be \alpha = 2\pi omegatw * (1-eigenval)
   ! MG check this, in the review I found a factor 2\pi, maybe it is reintroduced later
   bigomegatwsq(ii,1)=num*omegatw(ii)
 end do

 ABI_FREE(rhogg)
 ABI_FREE(mm)
 ABI_FREE(eigval)
 ABI_FREE(zz)
 ABI_FREE(eigvec)
 ABI_FREE(qplusg)

 call destroy_mpi_enreg(MPI_enreg_seq)

 write(msg,'(2a,f12.8,2a,3i5)')ch10,&
  ' cppm3par : omega twiddle minval [eV]  = ',MINVAL(ABS(omegatw))*Ha_eV,ch10,&
  '            omega twiddle min location = ',MINLOC(ABS(omegatw))
 call wrtout(std_out,msg)

end subroutine cppm3par
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/cppm4par
!! NAME
!! cppm4par
!!
!! FUNCTION
!! Calculate the plasmon-pole parameters using Engel-Farid model (PRB47,15931,1993) [[cite:Engel1993]].
!! See also Quasiparticle Calculations in Solids [[cite:Aulbur2001]] page. 23.
!!
!! INPUTS
!!  qpt(3)=Reduced coordinates of the q-point.
!!  npwc=number of plane waves in epsm1
!!  epsm1(npwc,npwc)=symmetrized inverse dielectric matrix.
!!  ngfftf(18)=contain all needed information about 3D fine FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  gvec(3,npwc)=G vectors in reduced coordinated
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!  rhor(nfftf)=charge density on the real space FFT grid
!!  nfftf=Number of FFT points.
!!
!! OUTPUT
!!  bigomegatwsq(npwc,npwc)=plasmon-pole strength.
!!  omegatw(npwc)=plasmon-pole frequencies.
!!
!! SOURCE

subroutine cppm4par(qpt, npwc, epsm1, ngfftf, gvec, gprimd, rhor, nfftf, bigomegatwsq, omegatw)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftf,npwc
!arrays
 integer,intent(in) :: gvec(3,npwc),ngfftf(18)
 real(dp),intent(in) :: gprimd(3,3),qpt(3), rhor(nfftf)
 complex(gwpc),intent(in) :: epsm1(npwc,npwc)
 complex(gwpc),intent(out) :: bigomegatwsq(npwc,npwc),omegatw(npwc)

!Local variables-------------------------------
!scalars
 integer :: ierr,ig,igp,ii,ngfft1,ngfft2,ngfft3,gmgp_idx
 real(dp) :: qpg_dot_qpgp
 character(len=500) :: msg
 character(len=80) :: bar
 type(MPI_type) :: MPI_enreg_seq
!arrays
 real(dp) :: b1(3),b2(3),b3(3),gppq(3),gpq(3),qlist(3,1)
 real(dp),allocatable :: eigval(:),qplusg(:),rhog_dp(:,:),tmp_rhor(:)
 complex(dpc),allocatable :: chi(:,:),chitmps(:,:), mm(:,:),mtemp(:,:),rhog(:), tmp1(:),zz2(:,:)
!*************************************************************************

 ! Calculate density in G space rhog(G)
 ! FIXME this has to be fixed, rho(G) should be passed instead of doing FFT for each q
 call initmpi_seq(MPI_enreg_seq)
 call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfftf(2),ngfftf(3),'all')

 ABI_MALLOC(rhog_dp, (2,nfftf))

 ! Conduct FFT tho(r)-->rhog(G)
 ABI_MALLOC(tmp_rhor,(nfftf))
 tmp_rhor = rhor ! To avoid having to use intent(inout).
 call fourdp(1,rhog_dp,tmp_rhor,-1,MPI_enreg_seq,nfftf,1,ngfftf,0)

 ABI_FREE(tmp_rhor)
 call destroy_mpi_enreg(MPI_enreg_seq)

 ABI_MALLOC(rhog, (nfftf))
 rhog(1:nfftf)=CMPLX(rhog_dp(1,1:nfftf),rhog_dp(2,1:nfftf))
 ABI_FREE(rhog_dp)

 ! Calculate the FFT index of each (G-Gp) vector and assign the value
 ! of the correspondent density simultaneously
 ngfft1=ngfftf(1)
 ngfft2=ngfftf(2)
 ngfft3=ngfftf(3)

 ABI_MALLOC_OR_DIE(mm, (npwc,npwc), ierr)

 ierr = 0
 do ig=1,npwc
   do igp=1,npwc
     gmgp_idx = g2ifft(gvec(:,ig)-gvec(:,igp),ngfftf)
     if (gmgp_idx /= 0) then
       mm(ig,igp) = rhog(gmgp_idx)
     else
       ierr = ierr + 1
       mm(ig,igp) = czero
     end if
   end do
 end do

 if (ierr /= 0) then
   write(msg,'(a,i0,3a)')&
    'Found ',ierr,' G1-G2 vectors falling outside the FFT box. ',ch10,&
    'Enlarge the FFT mesh to get rid of this problem. '
   ABI_WARNING(msg)
 end if

 ABI_FREE(rhog)

 ! Now we have rhogg, calculate the M matrix (q+G1).(q+G2) n(G1-G2)
 b1=two_pi*gprimd(:,1); b2=two_pi*gprimd(:,2); b3=two_pi*gprimd(:,3)
 do ig=1,npwc
   gpq(:)=gvec(:,ig)+qpt
   do igp=1,npwc
     gppq(:)=gvec(:,igp)+qpt
     qpg_dot_qpgp=zero
     do ii=1,3
       qpg_dot_qpgp = qpg_dot_qpgp + &
         ( gpq(1)*b1(ii) +gpq(2)*b2(ii) +gpq(3)*b3(ii))*&
         (gppq(1)*b1(ii)+gppq(2)*b2(ii)+gppq(3)*b3(ii))
     end do
     mm(ig,igp) = mm(ig,igp)*qpg_dot_qpgp
   end do ! igp
 end do ! ig

 ! Extract the reducible polarizability chi: e^{-1} = 1 + v chi
 ! \tilde\epsilon^{-1}_{G1 G2} = \delta_{G1 G2} + 4\pi \frac{\chi_{G1 G2}}{|q+G1| |q+G2|}
 !MG TODO too much memory in chi, we can do all this stuff inside a loop
 ABI_MALLOC_OR_DIE(chi, (npwc,npwc), ierr)
 ABI_MALLOC(qplusg, (npwc))

 chi(:,:)=epsm1(:,:)
 qlist(:,1) = qpt
 call cmod_qpg(1,1,qlist,npwc,gvec,gprimd,qplusg) !MG TODO here take care of small q

 do ig=1,npwc
   chi(ig,ig)=chi(ig,ig) - one
 end do

 do ig=1,npwc
   do igp=1,npwc
     chi(ig,igp) = chi(ig,igp) * qplusg(ig) * qplusg(igp) / four_pi
   end do
 end do

 ! Solve chi(w=)*X = Lambda M*X where Lambda=-1/em(q)**2
 ABI_MALLOC(eigval, (npwc))
 ABI_MALLOC_OR_DIE(mtemp, (npwc,npwc), ierr)

 ! Copy mm into working array as xhegv changes input matrices
 mtemp(:,:) = mm(:,:)

 call xhegv(1,"Vectors","Upper",npwc,chi,mtemp,eigval)
 ABI_FREE(mtemp)

 ! Now chi contains the eigenvectors.
 ! Eigenvectors are normalized as: X_i^* M X_j = \delta_{ij}

 ! Calculate the plasmon pole parameters
 ! good check: the lowest plasmon energy on gamma should be
 ! close to experimental plasma energy within an error of 10%
 ! this error can be reduced further if one includes the non local
 ! commutators in the calculation of the polarizability at q==0

 ABI_MALLOC(tmp1,(npwc))
 ABI_MALLOC_OR_DIE(zz2, (npwc, npwc), ierr)
 zz2(:,:)= zero

 ! Caller is responsible for handing small q case
 qlist(:,1) = qpt
 call cmod_qpg(1,1,qlist,npwc,gvec,gprimd,qplusg)

 do ii=1,npwc
   ! keeping in mind that the above matrix is negative definite
   ! we might have a small problem with the eigvals corresponding to large G vectors
   ! i.e. DM band index, where the eigevalues become very small with
   ! possibility of being small positive numbers (due to numerical problems)
   ! thus as a caution one can use the following condition
   ! this will not affect the result since such a huge plasmon energy give almost zero
   ! contribution to the self-energy correlation energy.

   if (eigval(ii)>=zero) then
     !write(msg,'(a,i0,a,es16.6)')' Imaginary plasmon pole eigenenergy, eigenvector number ',ii,' with eigval',eigval(ii),ch10
     !ABI_ERROR(msg)
     eigval(ii) = -1.0d-4
   end if

   ! Save plasmon energies omega_p(q)
   omegatw(ii) = SQRT(-one/eigval(ii))

   ! Calculate and save scaled plasmon-pole eigenvectors
   ! defined as \sqrt{4\pi} \frac{Mx}{\sqrt{\tilde\omega} |q+G|}
   tmp1(:)=chi(:,ii)

   do ig=1,npwc
     do igp=1,npwc
       zz2(ig,ii)=zz2(ig,ii)+mm(ig,igp)*tmp1(igp) ! z --> y
     end do
     bigomegatwsq(ig,ii)= SQRT(four_pi) * zz2(ig,ii) / SQRT(omegatw(ii)) / qplusg(ig)
   end do

 end do ! ii

 ABI_FREE(tmp1)
 ABI_FREE(eigval)
 ABI_FREE(zz2)
 ABI_FREE(qplusg)
 ABI_FREE(chi)
 ABI_FREE(mm)

 bar = repeat('-', 80)
 write(msg,'(3a)')bar,ch10,' plasmon energies in eV vs q vector shown for the lowest 10 bands'
 call wrtout(std_out,msg)
 write(msg,'(2x,5x,10f7.3)')(REAL(omegatw(ig))*Ha_eV, ig=1,min(10, npwc))
 call wrtout(std_out,msg)
 write(msg,'(a)')bar
 call wrtout(std_out,msg)

 write(msg,'(2a,f12.8,2a,3i5)')ch10,&
  ' cppm4par: omega twiddle minval [eV]  = ',MINVAL(ABS(omegatw))*Ha_eV,ch10,&
  '           omega twiddle min location = ',MINLOC(ABS(omegatw))
 call wrtout(std_out,msg)

end subroutine cppm4par
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/cqratio
!! NAME
!! cqratio
!!
!! FUNCTION
!!  Calculate qratio(G,Gp,q)= (q+G)\cdot(q+Gp) / |q+G|^2 needed for Hybertsen-Louie and Plasmonpole model
!!
!! INPUTS
!!  npwc=number of planewaves considered (used for the correlation part)
!!  gvec(3,npwc)=reduced coordinates of the plane waves
!!  q(3)=coordinates of q points
!!  gmet(3,3)=metric in reciprocal space
!!  gprimd(3,3)=reciprocal lattice vectors
!!
!! OUTPUT
!!  qratio(npwc,npwc)=(q+G).(q+Gp)
!!
!! SOURCE

subroutine cqratio(npwc, gvec, q, gmet, gprimd, qratio)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwc
!arrays
 integer,intent(in) :: gvec(3,npwc)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),q(3)
 real(dp),intent(out) :: qratio(npwc,npwc)

!Local variables ------------------------------
!scalars
 integer :: ig,igp,ii
 real(dp),parameter :: tol = 0.001_dp
 real(dp) :: qpg_dot_qpgp
!arrays
 real(dp) :: b1(3),b2(3),b3(3),gppq(3),gpq(3),norm(npwc)
!************************************************************************

 b1=two_pi*gprimd(:,1); b2=two_pi*gprimd(:,2); b3=two_pi*gprimd(:,3)

 norm(:)=zero; qratio=zero

 !FIXME this loops have to be rewritten!!!!
 do ig=1,npwc
   gpq(:)=gvec(:,ig)+q
   norm(ig)=two_pi*SQRT(DOT_PRODUCT(gpq,MATMUL(gmet,gpq)))
   !norm(ig)=normv(gpq,gmet,'g')
 end do

 do ig=1,npwc
   gpq(:)=gvec(:,ig)+q
   do igp=1,npwc
     gppq(:)=gvec(:,igp)+q
     qpg_dot_qpgp=zero
     !qpg_dot_qpgp=vdotw(gpq,gppq,gmet,'g')
     do ii=1,3
       qpg_dot_qpgp=qpg_dot_qpgp+&
        ( gpq(1)*b1(ii) +  gpq(2)*b2(ii) + gpq(3)*b3(ii))*&
        (gppq(1)*b1(ii) + gppq(2)*b2(ii) +gppq(3)*b3(ii))
     end do

     ! Now calculate qratio = (q+G).(q+Gp)/|q+G|^2
     ! when |q+G|^2 and (q+G).(q+Gp) are both zero set (q+G).(q+Gp)/|q+G|^2 = 1
     ! when |q+G|^2 is zero and |q+Gp| is not zero set (q+G).(q+Gp)/|q+G|^2 = 0
     if (norm(ig) < tol) then
       if (norm(igp) < tol) then     ! Case q=0 and G=Gp=0
         qratio(ig,igp) = one
       else                          ! Case q=0 and G=0 and Gp !=0
         qratio(ig,igp) = zero
       end if
     else if (norm(igp) < tol) then  ! Case q=0 and G= !0 and Gp=0
       qratio(ig,igp)=zero
     else
       qratio(ig,igp)=qpg_dot_qpgp / norm(ig)**2
     end if

   end do
 end do

end subroutine cqratio
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/ppm_calc_sigc
!!
!! NAME
!! ppm_calc_sigc
!!
!! FUNCTION
!!  Calculate the contribution to self-energy operator for a single band s in the band sum
!!  using a plasmon-pole model.
!!
!! INPUTS
!!  nspinor=Number of spinorial components.
!!  npwc=Number of G vectors in the plasmon pole.
!!  nomega=Number of frequencies.
!!  rhotwgp(npwx)=oscillator matrix elements divided by |q+G| i.e. $\frac{\langle b1 k-q s | e^{-i(q+G)r | b2 k s \rangle}{|q+G|}$.
!!  botsq(npwc,dm2_botsq)=Plasmon pole parameters for this q-point.
!!  otq(npwc,dm2_otq)=Plasmon pole parameters for this q-point.
!!  omegame0i(nomega)=Frequencies used to evaluate \Sigma_c ($\omega$ - $\epsilon_i)$
!!  zcut=Small imaginary part to avoid the divergence. (see related input variable)
!!  theta_mu_minus_e0i= $\theta(\mu-\epsilon_{k-q,b1,s}), defines if the state is occupied or not.
!!  eig(dm_eig,dm_eig)=The eigvectors of the symmetrized inverse dielectric matrix for this q point
!!    (first index for G, second index for bands).
!!  npwx=number of G vectors in rhotwgp
!!
!! OUTPUT
!!  ket(npwc,nomega):
!!
!!  i/two_pi * convolution between G and W ...
!!
!!  === model==1,2 ====
!!
!!    ket(G,omega) += Sum_G2                 Omega(G,G2) * rhotw(G2)
!!                            ---------------------------------------------------
!!                             2 omegatw(G,G2) (omega-E_i + omegatw(G,G2)(2f-1))
!!
!!  sigcme(nomega) (to be described), only relevant if ppm3 or ppm4
!!
!! TODO:
!!  Use BLAS for better efficiency
!!
!! SOURCE

subroutine ppm_calc_sigc(ppm, nspinor, npwc, nomega, rhotwgp, botsq, otq, &
                         omegame0i, zcut, theta_mu_minus_e0i, eig, npwx, ket, sigcme)

!Arguments ------------------------------------
!scalars
 class(ppmodel_t),intent(in) :: ppm
 integer,intent(in) :: nomega, npwc, npwx, nspinor
 real(dp),intent(in) :: theta_mu_minus_e0i, zcut
!arrays
 real(dp),intent(in) :: omegame0i(nomega)
 complex(gwpc),intent(in) :: botsq(npwc,ppm%dm2_botsq), eig(ppm%dm_eig,ppm%dm_eig), otq(npwc,ppm%dm2_otq)
 complex(gwpc),intent(in) :: rhotwgp(npwx*nspinor)
 complex(gwpc),intent(inout) :: ket(npwc*nspinor, nomega)
 complex(gwpc),intent(out) :: sigcme(nomega)

!Local variables-------------------------------
!scalars
 integer :: ig,igp,ii,ios,ispinor,spadc,spadx
 real(dp) :: den,ff,inv_den,omegame0i_io,otw,twofm1,twofm1_zcut
 complex(gwpc) :: ct, num, numf, rhotwgdp_igp
 logical :: fully_occupied,totally_empty
 !character(len=500) :: msg
!arrays
 complex(gwpc),allocatable :: rhotwgdpcc(:)
!*************************************************************************

 select case (ppm%model)

 case (PPM_GODBY_NEEDS, PPM_HYBERTSEN_LOUIE)
   fully_occupied = (abs(theta_mu_minus_e0i-one) < 0.001)
   totally_empty  = (abs(theta_mu_minus_e0i    ) < 0.001)

   do ispinor=1,nspinor
     spadx = (ispinor-1)*npwx; spadc = (ispinor-1)*npwc

     if (.not. totally_empty) then
       ! \Bomega^2_{G1G2}/\omegat_{G1G2} M_{G1,G2}. \theta(\mu-e_s) / (\omega+\omegat_{G1G2}-e_s-i\delta)
       twofm1_zcut = zcut
!$omp parallel do private(omegame0i_io, rhotwgdp_igp, otw, num, den)
       do ios=1,nomega
         omegame0i_io = omegame0i(ios)
         do igp=1,npwc
           rhotwgdp_igp = rhotwgp(spadx+igp)
           do ig=1,npwc
             otw = DBLE(otq(ig,igp)) !in principle otw -> otw - ieta
             num = botsq(ig,igp) * rhotwgdp_igp
             den = omegame0i_io + otw
             if (den**2 > zcut**2) then
               ket(spadc+ig,ios) = ket(spadc+ig,ios) + num/(den*otw) * theta_mu_minus_e0i
             else
               ket(spadc+ig,ios) = ket(spadc+ig,ios) + &
                 num * CMPLX(den,twofm1_zcut) / ((den**2+twofm1_zcut**2)*otw) * theta_mu_minus_e0i
             end if
           end do ! ig
         end do ! igp
       end do ! ios
     end if ! not totally empty

     if (.not. fully_occupied) then
       ! \Bomega^2_{G1G2}/\omegat_{G1G2} M_{G1,G2}. \theta(e_s-\mu) / (\omega-\omegat_{G1G2}-e_s+i\delta)
       twofm1_zcut = -zcut
!$omp parallel do private(omegame0i_io, rhotwgdp_igp, otw, num, den)
       do ios=1,nomega
         omegame0i_io = omegame0i(ios)
         do igp=1,npwc
           rhotwgdp_igp = rhotwgp(spadx+igp)
           do ig=1,npwc
             otw = DBLE(otq(ig,igp)) !in principle otw -> otw - ieta
             num = botsq(ig,igp) * rhotwgdp_igp
             den = omegame0i_io-otw
             if (den**2 > zcut**2) then
               ket(spadc+ig,ios) = ket(spadc+ig,ios) + num / (den*otw)*(one-theta_mu_minus_e0i)
             else
               ket(spadc+ig,ios) = ket(spadc+ig,ios) + &
                 num * CMPLX(den,twofm1_zcut) / ((den**2+twofm1_zcut**2)*otw)*(one-theta_mu_minus_e0i)
             end if
           end do ! ig
         end do ! igp
       end do ! ios
     end if ! not fully occupied

   end do ! ispinor

   ket=ket*half

 case (PPM_LINDEN_HORSH, PPM_ENGEL_FARID)
   ABI_CHECK(nspinor == 1, "nspinor/=1 not allowed")

   ! rho-twiddle(G) is formed, introduce rhotwgdpcc, for speed reason
   ABI_MALLOC(rhotwgdpcc, (npwx))

   ff = theta_mu_minus_e0i      ! occupation number f (include poles if ...)
   twofm1 = two*ff-one          ! 2f-1
   twofm1_zcut = twofm1*zcut
   rhotwgdpcc(:) = CONJG(rhotwgp(:))

   do ios=1,nomega
     omegame0i_io = omegame0i(ios)
     ct = czero_gw
     do ii=1,npwc ! Loop over the DM bands
       num = czero_gw

       select case (ppm%model)
       case (PPM_LINDEN_HORSH)
         ! Calculate \beta (eq. 106 pag 47)
         do ig=1,npwc
           num = num + rhotwgdpcc(ig)*eig(ig,ii)
         end do
         numf=num*CONJG(num) !MG this means that we cannot do SCGW
         numf=numf*botsq(ii,1)

       case (PPM_ENGEL_FARID)
         do ig=1,npwc
           num = num + rhotwgdpcc(ig)*botsq(ig,ii)
         end do
         numf = num*CONJG(num) !MG this means that we cannot do SCGW

       case default
         ABI_ERROR("Wrong ppm%model")
       end select

       otw=DBLE(otq(ii,1)) ! in principle otw -> otw - ieta
       den=omegame0i_io+otw*twofm1

       if (den**2 > zcut**2) then
         inv_den=one/den
         ct=ct+numf*inv_den
       else
         inv_den = one/((den**2+twofm1_zcut**2))
         ct = ct + numf*CMPLX(den,twofm1_zcut)*inv_den
       end if

     end do ! ii DM bands
     sigcme(ios) = ct*half

     !if (ppm%model == PPM_ENGEL_FARID) then
     !ct = dot_product(ket(:, ios), ket(:, ios))
     !if (abs(sigcme(ios) - ct) > tol12) then
     !  ABI_ERROR("foo bar")
     !end if
     !end if
   end do ! ios

   ABI_FREE(rhotwgdpcc)

 case default
   ABI_BUG(sjoin('Wrong ppm%model:',itoa(ppm%model)))
 end select

end subroutine ppm_calc_sigc
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/ppm_rotate_iqbz
!! NAME
!!  ppm_rotate_iqbz
!!
!! FUNCTION
!!  Symmetrize the plasmonpole parameters in the full BZ.
!!
!! INPUTS
!!  iq_bz=Index of the q-point in the BZ where the ppmodel parameters are wanted.
!!  Gsph<gsphere_t>=data related to the G-sphere.
!!  Cryst<crystal_t>=Info on the unit cell and crystal symmetries.
!!  Qmesh<kmesh_t>=the q-mesh used for the inverse dielectric matrix
!!  iq_ibz=Index of the q-point in the BZ.
!!  npwe=number of G vectors for the correlation part
!!  nomega=number of frequencies in $\epsilon^{-1}$
!!  omega=frequencies in epsm1_ggw
!!  epsm1_ggw(npwe,npwe,nomega)=the inverse dielctric matrix
!!  ngfftf(18)=contain all needed information about the 3D fine FFT mesh, see ~abinit/doc/variables/vargs.htm#ngfft
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!  nfftf=the number of points in the FFT mesh (for this processor)
!!  rhor_tot(nfftf)=the total charge in real space
!!
!! SIDE EFFECTS
!!  ppm<ppmodel_t>=data type containing information on the plasmonpole technique.
!!  Internal tables are modified so that they (point|store) the plasmon-pole parameters
!!  for the specified q-point in the BZ.
!!
!! SOURCE

subroutine ppm_rotate_iqbz(ppm, iq_bz, Cryst, Qmesh, Gsph, npwe, nomega, omega, epsm1_ggw, &
                           nfftf, ngfftf, rhor_tot)

!Arguments ------------------------------------
!scalars
 class(ppmodel_t),target,intent(inout) :: ppm
 integer,intent(in) :: nfftf,npwe,nomega,iq_bz
 type(crystal_t),intent(in) :: Cryst
 type(gsphere_t),intent(in) :: Gsph
 type(kmesh_t),intent(in) :: Qmesh
!arrays
 integer,intent(in) :: ngfftf(18)
 real(dp),intent(in) :: rhor_tot(nfftf)
 complex(dpc),intent(in) :: omega(nomega)
 complex(gwpc),intent(in) :: epsm1_ggw(npwe,npwe,nomega)

!Local variables-------------------------------
!scalars
 integer :: iq_ibz,itim_q,isym_q,iq_curr
 logical :: q_isirred
 !character(len=500) :: msg
!arrays
 real(dp) :: qbz(3)
! *********************************************************************

 ! Save the index of the q-point in the BZ for checking purpose.
 ppm%iq_bz = iq_bz

 call qmesh%get_bz_item(iq_bz, qbz, iq_ibz, isym_q, itim_q, isirred=q_isirred)
 iq_curr = iq_ibz; if (ppm%mqmem == 0) iq_curr = 1

 ! =======================================================
 ! ==== Branching for in-core or out-of-core solution ====
 ! =======================================================

 ! Allocate the tables for this q_ibz
 !print *, "ppm%has_qibz(iq_ibz)", ppm%has_qibz(iq_ibz), "q_isirred:", q_isirred
 if (ppm%has_qibz(iq_ibz) == PPM_NOTAB) call ppm%malloc_iqibz(iq_ibz)

 if (ppm%has_qibz(iq_ibz) == PPM_TAB_ALLOCATED) then
   ! Calculate the ppmodel tables for this q_ibz
   call ppm%new_setup(iq_ibz, Cryst, Qmesh, npwe, nomega, omega, epsm1_ggw, nfftf, Gsph%gvec, ngfftf, rhor_tot)
 end if

  ! Allocate memory if not done yet.
#ifdef FC_LLVM
  !FIXME I don't understand why LLVM fails here...
  !I put preproc so others know extra spaces are on purpose
  ABI_REMALLOC(ppm%bigomegatwsq_qbz_vals, (ppm%npwc, ppm%dm2_botsq) )
  ABI_REMALLOC(ppm%omegatw_qbz_vals, (ppm%npwc, ppm%dm2_otq) )
  ABI_REMALLOC(ppm%eigpot_qbz_vals, (ppm%dm_eig, ppm%dm_eig) )
#else
  ABI_REMALLOC(ppm%bigomegatwsq_qbz_vals, (ppm%npwc, ppm%dm2_botsq))
  ABI_REMALLOC(ppm%omegatw_qbz_vals, (ppm%npwc, ppm%dm2_otq))
  ABI_REMALLOC(ppm%eigpot_qbz_vals, (ppm%dm_eig, ppm%dm_eig))
#endif

 if (q_isirred) then
   ! Symmetrization is not needed. Copy the data in memory and change the status.
   ppm%bigomegatwsq_qbz_vals = ppm%bigomegatwsq(iq_ibz)%vals
   ppm%omegatw_qbz_vals = ppm%omegatw(iq_ibz)%vals
   ppm%eigpot_qbz_vals = ppm%eigpot(iq_ibz)%vals

 else
   ! q-point in the BZ. Calculate new table for this q-point in the BZ. Beware: Dimensions should not change.
   call ppm%get_qbz(Gsph, Qmesh, iq_bz, ppm%bigomegatwsq_qbz_vals, ppm%omegatw_qbz_vals, ppm%eigpot_qbz_vals)

   ! Release the table in the IBZ if required.
   if (.not. ppm%keep_qibz(iq_ibz)) call ppm%table_free_iqibz(iq_ibz)
 end if

end subroutine ppm_rotate_iqbz
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/ppm_new_setup
!! NAME
!! ppm_new_setup
!!
!! FUNCTION
!!  Initialize some values of several arrays of the ppm datastructure
!!  that are used in case of plasmonpole calculations
!!  Just a wrapper around different plasmonpole routines.
!!
!! INPUTS
!!  iq_ibz=Index of the q-point in the BZ.
!!  Cryst<crystal_t>=Info on the unit cell and crystal symmetries.
!!  Qmesh<kmesh_t>=the q-mesh used for the inverse dielectric matrix
!!  npwe=number of G vectors for the correlation part
!!  nomega=number of frequencies in $\epsilon^{-1}$
!!  omega=frequencies in epsm1_ggw
!!  epsm1_ggw(npwe,npwe,nomega)=the inverse dielctric matrix
!!  nfftf=the number of points in the FFT mesh (for this processor)
!!  ngfftf(18)=contain all needed information about the 3D fine FFT mesh, see ~abinit/doc/variables/vargs.htm#ngfft
!!  rhor(nfftf)=the total charge in real space.
!!
!! SIDE EFFECTS
!!  == if ppmodel 1 or 2 ==
!!   %omegatw and %bigomegatwsq
!!  == if ppmodel 3 ==
!!   %omegatw, %bigomegatwsq and %eigpot
!!  == if ppmodel 4 ==
!!   %omegatw and %bigomegatwsq
!!
!! NOTES
!! * FFT parallelism not implemented.
!! * TODO: rhor_tot should be replaced by rhog_tot to avoid nq_ibz FFTs.
!!
!! SOURCE

subroutine ppm_new_setup(ppm, iq_ibz, Cryst, Qmesh, npwe, nomega, omega, epsm1_ggw, nfftf, gvec, ngfftf, rhor_tot)

!Arguments ------------------------------------
!scalars
 class(ppmodel_t),intent(inout) :: ppm
 integer,intent(in) :: nfftf,npwe,nomega,iq_ibz
 type(kmesh_t),intent(in) :: Qmesh
 type(crystal_t),intent(in) :: Cryst
!arrays
 integer,intent(in) :: gvec(3,npwe),ngfftf(18)
 real(dp),intent(in) :: rhor_tot(nfftf)
 complex(dpc),intent(in) :: omega(nomega)
 complex(gwpc),intent(in) :: epsm1_ggw(npwe,npwe,nomega)

!Local variables-------------------------------
!scalars
 real(dp) :: n_at_G_zero
 character(len=500) :: msg
!scalars
 real(dp) :: qpt(3)
! *************************************************************************

 if (ppm%has_qibz(iq_ibz) /= PPM_TAB_ALLOCATED) then
   ABI_ERROR(sjoin("ppmodel tables for iq_ibz:", itoa(iq_ibz), "are not allocated! has_qibz=", itoa(ppm%has_qibz(iq_ibz))))
 end if

 qpt = Qmesh%ibz(:,iq_ibz)
 ppm%has_qibz(iq_ibz) = PPM_TAB_STORED

 ! Calculate plasmonpole parameters
 select case (ppm%model)

 case (PPM_NONE)
   ABI_COMMENT('Skipping plasmonpole model calculation')

 case (PPM_GODBY_NEEDS)
   ! Note: the q-dependence enters only through epsilon^-1.
   call cppm1par(npwe, nomega, omega, ppm%drude_plsmf, epsm1_ggw, ppm%omegatw(iq_ibz)%vals, ppm%bigomegatwsq(iq_ibz)%vals)

 case (PPM_HYBERTSEN_LOUIE)
   call cppm2par(qpt, npwe, epsm1_ggw(:,:,1), ngfftf, gvec, Cryst%gprimd, rhor_tot, nfftf, Cryst%gmet, &
                 ppm%bigomegatwsq(iq_ibz)%vals, ppm%omegatw(iq_ibz)%vals, ppm%invalid_freq)

   ! Quick-and-dirty change of the plasmon frequency. Never executed in standard runs.
   if (ppm%force_plsmf > tol6) then
      ! Integrate the real-space density
      n_at_G_zero = SUM(rhor_tot(:))/nfftf
      ! Change the prefactor
      write(msg,'(2(a,es16.8))') 'Forced ppmfreq: ',ppm%force_plsmf*Ha_eV,' nelect/ucvol: ',n_at_G_zero
      ABI_WARNING(msg)

      ppm%force_plsmf = (ppm%force_plsmf**2)/(four_pi*n_at_G_zero)
      ppm%bigomegatwsq(iq_ibz)%vals = ppm%force_plsmf * ppm%bigomegatwsq(iq_ibz)%vals
      ppm%omegatw(iq_ibz)%vals      = ppm%force_plsmf * ppm%omegatw(iq_ibz)%vals
      write(msg,'(a,es16.8)') 'Plasma frequency forced in HL ppmodel, new prefactor is: ',ppm%force_plsmf
      ABI_WARNING(msg)
   end if

 case (PPM_LINDEN_HORSH)
   call cppm3par(qpt, npwe,epsm1_ggw(:,:,1), ngfftf,gvec, Cryst%gprimd, rhor_tot, nfftf, &
                 ppm%bigomegatwsq(iq_ibz)%vals, ppm%omegatw(iq_ibz)%vals(:,1), ppm%eigpot(iq_ibz)%vals)

 case (PPM_ENGEL_FARID)
   if ((ALL(ABS(qpt)<1.0e-3))) qpt = GW_Q0_DEFAULT ! FIXME

   call cppm4par(qpt, npwe,epsm1_ggw(:,:,1), ngfftf, gvec, Cryst%gprimd, rhor_tot, nfftf, &
                 ppm%bigomegatwsq(iq_ibz)%vals, ppm%omegatw(iq_ibz)%vals(:,1))

 case default
   ABI_BUG(sjoin('Wrong ppm%model:', itoa(ppm%model)))
 end select

end subroutine ppm_new_setup
!!***

!!****f* m_ppmodel/ppm_print
!! NAME
!! ppm_print
!!
!! FUNCTION
!!  Print info on object
!!
!! SOURCE

subroutine ppm_print(ppm, units, header)

!Arguments ------------------------------------
 class(ppmodel_t),intent(in) :: ppm
 integer,intent(in) :: units(:)
 character(len=*),optional,intent(in) :: header

!Local variables-------------------------------
 character(len=500) :: msg
 type(yamldoc_t) :: ydoc
!*************************************************************************

 msg = ' ==== Info on the ppm_t object ==== '; if (present(header)) msg=' ==== '//trim(adjustl(header))//' ==== '
 call wrtout(units, msg)

 ydoc = yamldoc_open('Plasmonpole_params') !, width=11, real_fmt='(3f8.3)')
 !call ydoc%add_string("gwr_task", )
 call ydoc%add_int("dm2_botsq", ppm%dm2_botsq)
 call ydoc%add_int("dm_eig", ppm%dm_eig)
 call ydoc%add_int("dm2_otq", ppm%dm2_otq)
 call ydoc%add_int("invalid_freq", ppm%invalid_freq)
 call ydoc%add_int("model", ppm%model)
 call ydoc%add_int("mqmem", ppm%mqmem)
 call ydoc%add_int("nqibz", ppm%nqibz)
 call ydoc%add_int("npwc", ppm%npwc)
 call ydoc%add_int("userho", ppm%userho)
 call ydoc%add_int("iq_bz", ppm%iq_bz)
 call ydoc%add_real("drude_plsmf", ppm%drude_plsmf)
 call ydoc%add_real("force_plsmf", ppm%force_plsmf)
 !call ydoc%add_int1d("keep_qibz", ppm%keep_qibz)
 !call ydoc%add_int1d("has_qibz", ppm%has_qibz)

 call ydoc%write_units_and_free(units)

end subroutine ppm_print
!!***

!----------------------------------------------------------------------

end module m_ppmodel
!!***
