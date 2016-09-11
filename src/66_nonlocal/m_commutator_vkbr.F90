!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_commutator_vkbr
!! NAME
!!  m_commutator_vkbr
!!
!! FUNCTION
!!  This module provides objects and methods used to calculate the matrix elements 
!!  of the commutator [H,r] needed for the correct treatment of the optical limit q-->0
!!  in the matrix elements <k-q,b1|e^{-iqr}|k,b2> when non-local pseudopotentials are used.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2016 ABINIT group (MG, FB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_commutator_vkbr

 use defs_basis
 use defs_datatypes
 use m_blas
 use m_errors
 use m_profiling_abi

 use m_gwdefs,        only : czero_gw
 use m_fstrings,      only : sjoin, itoa
 use m_paw_sphharm,   only : ylmc, ylmcd
 use m_geometry,      only : normv
 use m_crystal,       only : crystal_t

 implicit none

 private
!!***

!!****t* m_commutator_vkbr/kb_form_factors
!! NAME
!! kb_form_factors
!!
!! FUNCTION
!!  The kb_form_factors data type stores basic dimensions and quantities 
!!  used in the GW part for the treatment of the non-analytic behavior of the 
!!  heads and wings of the irreducible polarizability in the long wave-length limit (i.e. q-->0).
!!
!! SOURCE

 type,public :: kb_form_factors

  integer :: npw       
  ! Number of plane waves for this  k-point.

  integer :: lnmax     
  ! Max. number of (l,n) components over all type of pseudos.

  integer :: ntypat    
  ! Number of type of atoms.

  integer,allocatable :: sign_dyad(:,:)
  ! sign_dyad(lnmax,ntypat). 
  ! sign of the KB dyadic product.

  real(dp),allocatable :: ff(:,:,:) 
  ! ff(npw,lnmax,ntypat) 
  ! KB form factor.

  real(dp),allocatable :: ffd(:,:,:)
  ! ffd(npw,lnmax,ntypat) 
  ! Derivative of ff wrt k+G.

 end type kb_form_factors

 !public :: correct_init_kb_ffactors
 !public :: ks_ffactors_free
!!***

!----------------------------------------------------------------------

!!****t* m_commutator_vkbr/kb_potential
!! NAME
!!
!! FUNCTION
!!  Matrix elements in |k+G> space needed for the 
!!  evaluation of the matrix elements of the commutator [Vnl,r] for the 
!!  optical limit in <kb1|e^{-iqr}|kb2>.
!!
!! SOURCE                                                                                   

 type,public :: kb_potential

  integer :: istwfk
  ! Storage mode of the G vectors for this k-point.

  integer :: ntypat
  ! Number of type of atoms

  integer :: natom
  ! Number of atoms

  integer :: mpsang
  ! Max l+1 over atoms

  integer :: npw
  ! Number of G-vectors.

  integer :: inclvkb
  ! Option for calculating the matrix elements of [Vnl,r].

  real(dp) :: kpoint(3)
  ! The k-point in reduced coordinates.

  complex(gwpc),allocatable :: fnl(:,:,:)
  ! fnl(npw,mpsang*mpsang,natom)

  complex(gwpc),allocatable :: fnld(:,:,:,:) 
  ! fnld(3,npw,mpsang*mpsang,natom)

 end type kb_potential

 public :: kb_potential_init
 public :: kb_potential_free
!!***

 interface kb_potential_free
   module procedure kp_potential_free_0D
   module procedure kp_potential_free_1D
 end interface kb_potential_free

 public ::  nc_ihr_comm

CONTAINS  !========================================================================================

!----------------------------------------------------------------------

!!****f* m_commutator_vkbr/correct_init_kb_ffactors
!! NAME
!!  correct_init_kb_ffactors
!!
!! FUNCTION
!!  Calculate KB form factors and derivatives required to evalute
!!  the matrix elements of the commutator [Vnl,r]-. 
!!  This term enters the expression for the oscillator strengths in 
!!  the optical limit q-->0. Pseudopotentials with more than one
!!  projector per angular channel are supported.
!!
!! INPUTS
!!  npw_k=Number of planewaves for this k-point.
!!  kpoint(3)=The kpoint in reduced coordinates.
!!  psps<pseudopotential_type>Structure gathering info on the pseudopotentials.
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  kg_k(3,npw_k)=Set of G-vectors for this k-point.
!!
!! OUTPUT
!!  KBff_k<KB_form_factor>=Structure storing the Kleynmann-Bylander form factor and derivatives for a single k-point.
!!
!! TODO 
!!  Replace old implementation with this new routine. Matrix elements
!!  of the commutator should be calculated on-the-fly in screening only
!!  if really needed. This is the first step toward the elimination
!!  of the KSS file. Modifications in cchi0q0 are needed.
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine correct_init_kb_ffactors(KBff_k,cryst,psps,kpoint,npw_k,kg_k)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'correct_init_kb_ffactors'
 use interfaces_41_geometry
 use interfaces_66_nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_k
 type(kb_form_factors),intent(inout) :: KBff_k
 type(crystal_t),intent(in) :: Cryst
 type(Pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: kg_k(3,npw_k)
 real(dp),intent(in) :: kpoint(3)

!Local variables ------------------------------
!scalars
 integer :: dimffnl,ider,idir,itypat,nkpg
 integer :: il,ilmn,iln,iln0,nlmn,ig
 real(dp) :: fact
 !real(dp) :: ,effmass,ecutsm
 !character(len=500) :: msg
!arrays
 real(dp) :: kpg(3) 
 real(dp),allocatable :: ffnl(:,:,:,:),kpg_dum(:,:),modkplusg(:)
 real(dp),allocatable :: ylm(:,:),ylm_gr(:,:,:),ylm_k(:,:)

! *************************************************************************

 DBG_ENTER("COLL")
 ABI_CHECK(psps%usepaw==0,"You should not be here!")

 ! === Save KB dyadic sign (integer-valued) ===
 ! * Notice how the ordering is chosen correctly unlike in outkss.
 ! * More than one projector per angular channel is allowed but changes in cchi0q0 are needed
 !  allocate(KBff_ksign(psps%mpsang,ntypat))  THIS THE OLD IMPLEMENTATION
 ABI_MALLOC(KBff_k%sign_dyad, (psps%lnmax,Psps%ntypat))
 KBff_k%sign_dyad(:,:) = 0

 do itypat=1,psps%ntypat
   nlmn = count(psps%indlmn(3,:,itypat)>0)
   iln0=0 
   do ilmn=1,nlmn
     iln=psps%indlmn(5,ilmn,itypat)
     if (iln>iln0) then
       iln0=iln
       KBff_k%sign_dyad(iln,itypat)=NINT(DSIGN(one,psps%ekb(ilmn,itypat)))
     end if
   end do
 end do

 KBff_k%npw    = npw_k
 KBff_k%lnmax  = psps%lnmax
 KBff_k%ntypat = psps%ntypat

 ! === Allocate KB form factor and derivative wrt k+G ===
 ! * Also here we use correct ordering for dimensions
 ABI_MALLOC(KBff_k%ff ,(npw_k,psps%lnmax,Psps%ntypat))
 ABI_MALLOC(KBff_k%ffd,(npw_k,psps%lnmax,Psps%ntypat))
 KBff_k%ff(:,:,:)=zero ; KBff_k%ffd(:,:,:)=zero
 
 ider=1 ; dimffnl=2 ! To retrieve the first derivative.
 idir=0 ; nkpg=0
 !
 ! Quantities used only if useylm==1
 ABI_MALLOC(ylm,(npw_k,psps%mpsang**2*Psps%useylm))
 ABI_MALLOC(ylm_gr,(npw_k,3+6*(ider/2),psps%mpsang**2*Psps%useylm))
 ABI_MALLOC(ylm_k,(npw_k,psps%mpsang**2*Psps%useylm))
 ABI_MALLOC(kpg_dum,(npw_k,nkpg))

 ABI_MALLOC(ffnl,(npw_k,dimffnl,psps%lmnmax,Psps%ntypat))

 call mkffnl(psps%dimekb,dimffnl,Psps%ekb,ffnl,Psps%ffspl,cryst%gmet,cryst%gprimd,ider,idir,Psps%indlmn,&
   kg_k,kpg_dum,kpoint,psps%lmnmax,Psps%lnmax,Psps%mpsang,Psps%mqgrid_ff,nkpg,npw_k,& 
   psps%ntypat,Psps%pspso,Psps%qgrid_ff,cryst%rmet,Psps%usepaw,Psps%useylm,ylm_k,ylm_gr)

 ABI_FREE(kpg_dum)
 ABI_FREE(ylm)
 ABI_FREE(ylm_gr)
 ABI_FREE(ylm_k)

 ABI_MALLOC(modkplusg,(npw_k))

 !effmass=one; ecutsm=zero
 !call mkkin(ecut,ecutsm,effmass,cryst%gmet,kg_k,modkplusg,kpoint,npw_k, 0, 0)
 !modkplusg(:)=SQRT(half/pi**2*modkplusg(:))
 !modkplusg(:)=MAX(modkplusg(:),tol10)

 do ig=1,npw_k
   kpg(:)= kpoint(:)+kg_k(:,ig)
   modkplusg(ig) = normv(kpg,cryst%gmet,"G")
 end do

 do itypat=1,psps%ntypat
   nlmn=COUNT(psps%indlmn(3,:,itypat)>0)
   iln0=0 
   do ilmn=1,nlmn
     il= psps%indlmn(1,ilmn,itypat)+1
     iln=psps%indlmn(5,ilmn,itypat)

     if (iln>iln0) then
       iln0=iln

       if (ABS(psps%ekb(ilmn,itypat))>1.0d-10) then
         SELECT CASE (il)
         CASE (1)
           KBff_k%ff (1:npw_k,iln,itypat) = ffnl(:,1,ilmn,itypat)
           KBff_k%ffd(1:npw_k,iln,itypat) = ffnl(:,2,ilmn,itypat)*modkplusg(:)/two_pi

         CASE (2)
           KBff_k%ff (1:npw_k,iln,itypat) =   ffnl(:,1,ilmn,itypat)*modkplusg(:)
           KBff_k%ffd(1:npw_k,iln,itypat) = ((ffnl(:,2,ilmn,itypat)*modkplusg(:)**2)+&
                                              ffnl(:,1,ilmn,itypat))/two_pi
         CASE (3)
           KBff_k%ff (1:npw_k,iln,itypat) =  ffnl(:,1,ilmn,itypat)*modkplusg(:)**2
           KBff_k%ffd(1:npw_k,iln,itypat) = (ffnl(:,2,ilmn,itypat)*modkplusg(:)**3+&
                                           2*ffnl(:,1,ilmn,itypat)*modkplusg(:))/two_pi
         CASE (4)
           KBff_k%ff (1:npw_k,iln,itypat) =  ffnl(:,1,ilmn,itypat)*modkplusg(:)**3
           KBff_k%ffd(1:npw_k,iln,itypat) = (ffnl(:,2,ilmn,itypat)*modkplusg(:)**4+&
                                           3*ffnl(:,1,ilmn,itypat)*modkplusg(:)**2)/two_pi
         CASE DEFAULT
           MSG_ERROR('l greater than g not implemented.')
         END SELECT

         fact = SQRT(four_pi/cryst%ucvol*(2*il-1)*ABS(psps%ekb(ilmn,itypat)))
         KBff_k%ff (:,iln,itypat) = fact * KBff_k%ff (:,iln,itypat)
         KBff_k%ffd(:,iln,itypat) = fact * KBff_k%ffd(:,iln,itypat)

       else ! ekb==0
         KBff_k%ff (:,iln,itypat)=zero
         KBff_k%ffd(:,iln,itypat)=zero
       end if

     end if
   end do
 end do

 ABI_FREE(ffnl)
 ABI_FREE(modkplusg)

 DBG_EXIT("COLL")

end subroutine correct_init_kb_ffactors
!!***

!----------------------------------------------------------------------

!!****f* m_commutator_vkbr/ks_ffactors_free
!! NAME
!!  ks_ffactors_free
!!
!! FUNCTION
!!  Free the memory allocated in a structure of type kb_form_factors
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine ks_ffactors_free(KBff_k)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ks_ffactors_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(kb_form_factors),intent(inout) :: KBff_k

!************************************************************************

 !@kb_form_factors

!integer 
 if (allocated(KBff_k%sign_dyad)) then
   ABI_FREE(KBff_k%sign_dyad)
 end if

!real
 if (allocated(KBff_k%ff )) then
   ABI_FREE(KBff_k%ff)
 end if
 if (allocated(KBff_k%ffd)) then
   ABI_FREE(KBff_k%ffd)
 end if

end subroutine ks_ffactors_free
!!***

!----------------------------------------------------------------------

!!****f* m_commutator_vkbr/kb_potential_init
!! NAME
!!  kb_potential_init
!!
!! FUNCTION
!!  Creation method the the kb_potential structures datatype.
!!
!! INPUTS
!!  cryst<crystal_t>=Datatype gathering info on the crystal structure.
!!  psps<pseudopotential_type>Structure gathering info on the pseudopotentials.
!!  inclvkb=Option defining the algorithm used for the application of [Vnl,r].
!!    2 for Spherical harmonics
!!  istwfk=Storage mode for the wavefunctions at this k-point.
!!  npw=Number of planewaves in <k+G1|[Vnl,r]|k+G2>
!!  kpoint(3)=K-point of interest in reduced coordinates.
!!  gvec(3,npw)=Reduced coordinates of the G-vectors.
!!
!! OUTPUT
!!  KBgrad_k<kb_potential>=Structure containing arrays needed for calculating <\psi_1|[Vnl,r]\psi_2>.
!!    Completely initialized in output.
!!
!! PARENTS
!!      calc_optical_mels,cchi0q0,cchi0q0_intraband
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine kb_potential_init(KBgrad_k,cryst,psps,inclvkb,istwfk,npw,kpoint,gvec)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'kb_potential_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw,inclvkb,istwfk
 type(crystal_t),intent(in) :: cryst
 type(kb_potential),intent(inout) :: KBgrad_k
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: gvec(3,npw)
 real(dp),intent(in) :: kpoint(3)

!Local variables-------------------------------
!scalars
 integer :: ierr
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: vkb(:,:,:),vkbd(:,:,:),vkbsign(:,:)

!************************************************************************

 !@kb_potential
 KBgrad_k%istwfk  = istwfk
 KBgrad_k%ntypat  = cryst%ntypat
 KBgrad_k%natom   = cryst%natom
 KBgrad_k%mpsang  = psps%mpsang
 KBgrad_k%npw  = npw
 KBgrad_k%inclvkb = inclvkb
 KBgrad_k%kpoint  = kpoint
 !
 ! Calculate KB form factors and derivatives.
 !
 ! The following arrays should be allocated using lnmax instead of mpsang
 ! in order to support pseudos with more than projector.
 ! Moreover they should be calculated on-the-fly using calc_vkb
 ! For the moment, we opt for a quick an dirty implementation.
 ABI_MALLOC(vkbsign,(psps%mpsang, cryst%ntypat))
 ABI_MALLOC(vkb ,(npw, psps%mpsang, cryst%ntypat))
 ABI_MALLOC(vkbd,(npw, psps%mpsang, cryst%ntypat))

 call calc_vkb(cryst,psps,kpoint,npw,gvec,vkbsign,vkb,vkbd)
 
 SELECT CASE (inclvkb)
 CASE (2) ! * Complex spherical harmonics (CPU and mem \propto npw).

   write(msg,'(a,f12.1)')'out-of-memory in fnl; Mb= ',one*npw*psps%mpsang**2*cryst%natom*2*gwpc*b2Mb
   ABI_STAT_MALLOC(KBgrad_k%fnl,(npw,psps%mpsang**2,cryst%natom), ierr)
   ABI_CHECK(ierr==0, msg)

   write(msg,'(a,f12.1)')'out-of-memory in fnld; Mb= ',three*npw*psps%mpsang**2*cryst%natom*2*gwpc*b2Mb
   ABI_STAT_MALLOC(KBgrad_k%fnld,(3,npw,psps%mpsang**2,cryst%natom), ierr)
   ABI_CHECK(ierr==0, msg)

   call ccgradvnl_ylm(cryst,psps,npw,gvec,kpoint,vkbsign,vkb,vkbd,KBgrad_k%fnl,KBgrad_k%fnld)

 CASE DEFAULT
   MSG_ERROR(sjoin("Wrong inclvkb= ",itoa(inclvkb)))
 END SELECT 

 ABI_FREE(vkbsign)
 ABI_FREE(vkb)
 ABI_FREE(vkbd)

end subroutine kb_potential_init 
!!***

!----------------------------------------------------------------------

!!****f* m_commutator_vkbr/kp_potential_free_0D
!! NAME
!!  kp_potential_free_0D
!!
!! FUNCTION
!!  Free all memory allocated in a structure of type kb_potential
!!
!! PARENTS
!!      m_commutator_vkbr
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine kp_potential_free_0D(KBgrad_k)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'kp_potential_free_0D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(kb_potential),intent(inout) :: KBgrad_k

!************************************************************************

 !@kb_potential

!complex
 if (allocated(KBgrad_k%fnl)) then
   ABI_FREE(KBgrad_k%fnl)
 end if
 if (allocated(KBgrad_k%fnld)) then
   ABI_FREE(KBgrad_k%fnld)
 end if

end subroutine kp_potential_free_0D
!!***

!----------------------------------------------------------------------

!!****f* m_commutator_vkbr/kp_potential_free_1D
!! NAME
!!  kp_potential_free_1D
!!
!! FUNCTION
!!  Free all memory allocated in a structure of type kb_potential
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine kp_potential_free_1D(KBgrad_k)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'kp_potential_free_1D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(kb_potential),intent(inout) :: KBgrad_k(:)

!Local variables ------------------------------
!scalars
 integer :: ii

!************************************************************************

 do ii=1,SIZE(KBgrad_k)
   call kp_potential_free_0D(KBgrad_k(ii))
 end do

end subroutine kp_potential_free_1D
!!***

!----------------------------------------------------------------------

!!****f* m_commutator_vkbr/add_vnlr_commutator
!! NAME
!!  add_vnlr_commutator
!!
!! FUNCTION
!!  Calculate the matrix elements of the dipole operator <phi1|r|phi2>.
!!  For norm conserving potentials the commutator [Vnl,r] is included according to inclvkb. 
!!
!! INPUTS
!!  KBgrad_k<kb_potential>
!!  cryst<crystal_t>=Datatype gathering info on the crystal structure.
!!  npw=Number of G for wavefunctions.
!!  nspinor=Number of spinorial components.
!!  ug1(npw*nspinor)=Left wavefunction.
!!  ug2(npw*nspinor)=Right wavefunction
!!
!! SIDE EFFECTS
!!  rhotwx(3,nspinor**2)= Updated. Matrix elements in reduced coordinates, see NOTES below.
!!
!! NOTES
!!   1) <k b1|e^{-iq.r}|k b2> = \delta_{b1 b2} -iq <k b1|r|k b2> =  \delta_{b1 b2} -iq ( <k b1| [H,r] |k b2> / (e1-e2) ).
!!
!!      This routine calculates the matrix elements of ir*(e1-e2) 
!!      Remember that [H,r] = -\nabla + [V_nl,r]
!!
!!  2) The Fourier transform of a two-point real function f(r1,r2) satisfies:
!!      a) f_{\Gamma}(G1,G2) = f_{\Gamma}(-G1,-G2)^*
!!      b) f_{G0/2}  (G1,G2) = f_{G0/2}(-G1-G0,-G2-G0)^*
!!
!! TODO
!!  *) Spinorial case is not implemented.
!!
!! PARENTS
!!      m_commutator_vkbr
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine add_vnlr_commutator(KBgrad_k,cryst,npw,nspinor,ug1,ug2,rhotwx)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'add_vnlr_commutator'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw,nspinor
 type(kb_potential),intent(in) :: KBgrad_k
 type(crystal_t),intent(in) :: cryst
!arrays
 complex(gwpc),target,intent(in) :: ug1(npw*nspinor),ug2(npw*nspinor)
 complex(gwpc),intent(inout) :: rhotwx(3,nspinor**2)

!Local variables ------------------------------
!scalars
 integer :: ig1,ig2,iat,ig,ilm,itypat
 complex(gwpc) :: cta1,cta4,ct
!arrays
 complex(gwpc) :: dum(3),cta2(3),cta3(3),gamma_term(3)

!************************************************************************

 ! === Adding term i <c,k|[Vnl,r]|v,k> === 
 ! * Two different algorithms are coded, the second one is much faster

 SELECT CASE (KBgrad_k%inclvkb)
 CASE (2)  ! Complex spherical harmonics (much faster!).

  dum=czero_gw; gamma_term=czero
!TODO this section causes segmentation faults
  !!!!!!!!$OMP PARALLEL DO PRIVATE(cta1,cta2,cta3,cta4) COLLAPSE(2) REDUCTION(+:dum,gamma_term)
  do iat=1,KBgrad_k%natom
    itypat = cryst%typat(iat)
    !psps%indlmn(6, lmnmax, itypat)
    ! array giving l,m,n,lm,ln,spin for i=ln  (if useylm=0)
    do ilm=1,KBgrad_k%mpsang**2
      cta1 = czero_gw; cta2(:) = czero_gw
      cta4 = czero_gw; cta3(:) = czero_gw
      do ig=1,npw 
        ! Here we take advantage of the property Y_(l-m)= (-i)^m Y_lm^*.
        cta1   = cta1    + ug1(ig) * KBgrad_k%fnl (ig,ilm,iat)
        cta2(:)= cta2(:) + ug2(ig) * KBgrad_k%fnld(:,ig,ilm,iat)
        cta3(:)= cta3(:) + ug1(ig) * KBgrad_k%fnld(:,ig,ilm,iat)
        cta4   = cta4    + ug2(ig) * KBgrad_k%fnl (ig,ilm,iat)
        if (ig==1) gamma_term = gamma_term + CONJG(cta1)*cta2(:) +CONJG(cta3(:))*cta4
      end do
      dum(:)= dum(:) +CONJG(cta1)*cta2(:) +CONJG(cta3(:))*cta4
    end do
  end do

  if (KBgrad_k%istwfk>1) then
    dum = two * j_dpc * AIMAG(dum); if (KBgrad_k%istwfk==2) dum = dum - j_dpc * AIMAG(gamma_term)
  end if

  rhotwx(:,1) = rhotwx(:,1) + dum(:)

 CASE DEFAULT
   MSG_ERROR(sjoin("Wrong inclvkb:", itoa(KBgrad_k%inclvkb)))
 END SELECT

end subroutine add_vnlr_commutator
!!***

!----------------------------------------------------------------------

!!****f* m_commutator_vkbr/calc_vkb
!! NAME
!!  calc_vkb
!!
!! FUNCTION
!!  This routine calculates the Kleynman-Bylander form factors and its derivatives 
!!  needed for the evaluation of the matrix elements of the dipole operator <phi1|r|phi2>.
!!
!! INPUTS
!!  cryst<crystal_t>=Crystalline structure
!!  psps<pseudopotential_type>=Structured datatype gathering information on the pseudopotentials.
!!  kpoint(3)=The k-point in reduced coordinates.
!!  npw_k=Number of plane waves for this k-point.
!!  kg_k(3,npw_k)=Reduced coordinates of the G-vectors.
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!
!! OUTPUT
!!  vkb (npw_k,Psps%mpsang,psps%ntypat)=KB form factors.
!!  vkbd(npw_k,Psps%mpsang,psps%ntypat)=KB form factor derivatives.
!!  vkbsign(psps%mpsang,Psps%ntypat)   =KS dyadic sign.
!!
!! NOTES
!!  This piece of code has been extracted from outkss.F90. The implementation is consistent
!!  with the KSS file formata (Fortran version) but it presents two design flaws.
!!
!!   1) Pseudo with more that one projector per l-channel are not supported.
!!
!! TODO
!!  *) Spinorial case is not implemented.
!!  *) Fix the above mentioned programming sins (KSS FORTRAN fileformat has to be modified though)
!!
!! PARENTS
!!      m_commutator_vkbr,m_io_kss
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine calc_vkb(cryst,psps,kpoint,npw_k,kg_k,vkbsign,vkb,vkbd)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_vkb'
 use interfaces_41_geometry
 use interfaces_56_recipspace
 use interfaces_66_nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_k
 type(crystal_t),intent(in) :: cryst
 type(Pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: kg_k(3,npw_k)
 real(dp),intent(in) :: kpoint(3) 
 real(dp),intent(out) :: vkb (npw_k,Psps%mpsang,psps%ntypat)
 real(dp),intent(out) :: vkbd(npw_k,Psps%mpsang,psps%ntypat)
 real(dp),intent(out) :: vkbsign(psps%mpsang,Psps%ntypat)

!Local variables ------------------------------
!scalars
 integer :: dimffnl,ider,idir,itypat,nkpg,il0,in,il,ilmn,ig
 real(dp) :: effmass,ecutsm,ecut
!arrays
 real(dp),allocatable :: ffnl(:,:,:,:),kpg_dum(:,:),modkplusg(:)
 real(dp),allocatable :: ylm_gr(:,:,:),ylm_k(:,:)

! *************************************************************************

 DBG_ENTER("COLL")
 ABI_CHECK(psps%usepaw==0,"You should not be here!")
 ABI_CHECK(psps%useylm==0,"useylm/=0 not considered!")
 !
 ! === Save KB dyadic sign (integer-valued) ===
 vkbsign=zero
 do itypat=1,psps%ntypat
   il0=0 
   do ilmn=1,psps%lmnmax
     il=1+psps%indlmn(1,ilmn,itypat)
     in=psps%indlmn(3,ilmn,itypat)
     if (il/=il0 .and. in==1) then
       il0=il
       vkbsign(il,itypat)=DSIGN(one,psps%ekb(ilmn,itypat))
     end if
   end do
 end do

 ! === Allocate KB form factor and derivative wrt k+G ===
 ! * Here we do not use correct ordering for dimensions
 
 ider=1; dimffnl=2 ! To retrieve the first derivative.
 idir=0; nkpg=0
 !
 ! Quantities used only if useylm==1
 ABI_MALLOC(ylm_k, (npw_k, psps%mpsang**2*Psps%useylm))
 ABI_MALLOC(ylm_gr,(npw_k, 3+6*(ider/2),psps%mpsang**2*Psps%useylm))
 ABI_MALLOC(kpg_dum, (npw_k, nkpg))
 ABI_MALLOC(ffnl, (npw_k,dimffnl, psps%lmnmax, psps%ntypat))

 call mkffnl(psps%dimekb,dimffnl,Psps%ekb,ffnl,Psps%ffspl,cryst%gmet,cryst%gprimd,ider,idir,Psps%indlmn,&
   kg_k,kpg_dum,kpoint,psps%lmnmax,Psps%lnmax,Psps%mpsang,Psps%mqgrid_ff,nkpg,npw_k,& 
   psps%ntypat,Psps%pspso,Psps%qgrid_ff,cryst%rmet,Psps%usepaw,Psps%useylm,ylm_k,ylm_gr)

 ABI_FREE(ylm_k)
 ABI_FREE(ylm_gr)
 ABI_FREE(kpg_dum)

 ABI_MALLOC(modkplusg,(npw_k))
 effmass=one; ecutsm=zero; ecut=HUGE(one)
 call mkkin(ecut,ecutsm,effmass,cryst%gmet,kg_k,modkplusg,kpoint,npw_k,0,0)
 modkplusg(:) = SQRT(half/pi**2*modkplusg(:))
 modkplusg(:) = MAX(modkplusg(:),tol10)

 !do ig=1,npw_k
 ! kpg(:)= kpoint(:)+kg_k(:,ig)
 ! modkplusg(ig) = normv(kpg,cryst%gmet,"G")
 !end do

 ! Calculate matrix elements.
 vkb=zero; vkbd=zero

 do itypat=1,psps%ntypat
   il0=0
   do ilmn=1,psps%lmnmax
     il = 1+psps%indlmn(1,ilmn,itypat)
     in = psps%indlmn(3,ilmn,itypat)
     if ((il/=il0).and.(in==1)) then
       il0=il
       if (ABS(psps%ekb(ilmn,itypat))>1.0d-10) then
         if (il==1) then
           vkb (1:npw_k,il,itypat) = ffnl(:,1,ilmn,itypat)
           vkbd(1:npw_k,il,itypat) = ffnl(:,2,ilmn,itypat)*modkplusg(:)/two_pi
         else if (il==2) then
           vkb(1:npw_k,il,itypat)  = ffnl(:,1,ilmn,itypat)*modkplusg(:)
           do ig=1,npw_k
             vkbd(ig,il,itypat) = ((ffnl(ig,2,ilmn,itypat)*modkplusg(ig)*modkplusg(ig))+&
              ffnl(ig,1,ilmn,itypat) )/two_pi
           end do
         else if (il==3) then
           vkb (1:npw_k,il,itypat) =  ffnl(:,1,ilmn,itypat)*modkplusg(:)**2
           vkbd(1:npw_k,il,itypat) = (ffnl(:,2,ilmn,itypat)*modkplusg(:)**3+&
            2*ffnl(:,1,ilmn,itypat)*modkplusg(:) )/two_pi
         else if (il==4) then
           vkb (1:npw_k,il,itypat) =  ffnl(:,1,ilmn,itypat)*modkplusg(:)**3
           vkbd(1:npw_k,il,itypat) = (ffnl(:,2,ilmn,itypat)*modkplusg(:)**4+&
            3*ffnl(:,1,ilmn,itypat)*modkplusg(:)**2 )/two_pi
         end if
         vkb (:,il,itypat) = SQRT(4*pi/cryst%ucvol*(2*il-1)*ABS(psps%ekb(ilmn,itypat)))*vkb (:,il,itypat)
         vkbd(:,il,itypat) = SQRT(4*pi/cryst%ucvol*(2*il-1)*ABS(psps%ekb(ilmn,itypat)))*vkbd(:,il,itypat)
       else
         vkb (:,il,itypat)=zero
         vkbd(:,il,itypat)=zero
       end if
     end if
   end do
 end do

 ABI_FREE(ffnl)
 ABI_FREE(modkplusg)

 DBG_EXIT("COLL")

end subroutine calc_vkb
!!***

!----------------------------------------------------------------------

!!****f* m_commutator_vkbr/nc_ihr_comm
!! NAME
!!  nc_pwihr_comm
!!
!! FUNCTION
!!  Calculate the matrix elements of the commutator i[H,r]  
!!  For norm conserving potentials the commutator i[Vnl,r] is included depending on inclvkb. 
!!
!! INPUTS
!!  KBgrad_k<kb_potential>
!!  cryst<crystal_t>=Unit cell and symmetries
!!  nspinor=Number of spinorial components.
!!  npw=Number of G for wavefunctions.
!!  istwfk=Storage mode for wavefunctions.
!!  inclvkb=Option defining whether [Vnl,r] is added or not.
!!  kpoint(3)=k-point in reduced coordinates.
!!  ug1(npw*nspinor)=Left wavefunction.
!!  ug2(npw*nspinor)=Right wavefunction
!!  gvec(3,npw)=Planes waves for wavefunctions.
!!
!! OUTPUT
!!  ihr_comm(3,nspinor**2)= Matrix elements of the commutator i[H,r] between the input states.
!!   Result is in reduced coordinates. ug1 and ug2 are supposed to be orthogonal. 
!!
!! NOTES
!!  <k b1|e^{-iq.r}|k b2> = \delta_{b1 b2} -iq <k b1|r|k b2> =  \delta_{b1 b2} -iq ( <k b1| [H,r] |k b2> / (e1-e2) ).
!!  Remember that [H,r] = -\nabla + [V_nl,r]
!!
!! TODO
!!  *) Spinorial case is not implemented.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function nc_ihr_comm(Kbgrad_k,cryst,npw,nspinor,istwfk,inclvkb,kpoint,ug1,ug2,gvec) result(ihr_comm)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nc_ihr_comm'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw,nspinor,inclvkb,istwfk
 type(kb_potential),intent(in) :: KBgrad_k
 type(crystal_t),intent(in) :: Cryst
!arrays
 integer,intent(in) :: gvec(3,npw)
 real(dp),intent(in) :: kpoint(3)
 complex(gwpc),intent(in) :: ug1(npw*nspinor),ug2(npw*nspinor)
 complex(gwpc) :: ihr_comm(3,nspinor**2)

!Local variables ------------------------------
!scalars
 integer :: ig,iab,spad1,spad2
 complex(dpc) :: c_tmp
!arrays
 integer :: spinorwf_pad(2,4)

!************************************************************************

 ! [H, r] = -\nabla + [V_{nl}, r]
 ! V_nl is present only in the case of NC pseudos but 
 ! not in PAW unless even the AE Hamiltonian in non-local e.g. LDA+U or LEXX.

 ! -i <c,k|\nabla_r|v,k> in reduced coordinates is always included.
 ! -i <c,k|\nabla_r|v,k> = \sum_G u_{ck}^*(G) [k+G] u_{vk}(G)
 ! Note that here we assume c/=v, moreover the ug are supposed to be orthonormal and 
 ! hence k+G can be replaced by G.

 spinorwf_pad=RESHAPE((/0,0,npw,npw,0,npw,npw,0/),(/2,4/))

 ihr_comm = czero

 ! -i <c,k|\nabla_r|v,k> in reduced coordinates.
 ! FIXME: This term is  spin diagonal!
 if (istwfk==1) then
   do iab=1,nspinor**2
     spad1 = spinorwf_pad(1,iab)
     spad2 = spinorwf_pad(2,iab)
     do ig=1,npw 
       c_tmp = CONJG(ug1(ig+spad1)) * ug2(ig+spad2)
       ihr_comm(:,iab) = ihr_comm(:,iab) + c_tmp*gvec(:,ig)
     end do
   end do
 else
   ! Symmetrized expression: \sum_G  (k+G) 2i Im [ u_a^*(G) u_b(G) ]. (k0,G0) term is null.
   do ig=1,npw 
     c_tmp = CONJG(ug1(ig)) * ug2(ig)
     ihr_comm(:,1) = ihr_comm(:,1) + two*j_dpc * AIMAG(c_tmp) * (kpoint + gvec(:,ig))
   end do  
 end if

 ! Add second term $i <c,k|[Vnl,r]|v,k> in$ reduced cordinates.
 if (inclvkb/=0) then 
   ABI_CHECK(nspinor == 1, "nspinor/=1 not coded")
   ABI_CHECK(istwfk == KBgrad_k%istwfk, "input istwfk /= KBgrad_k%istwfk")
   !ABI_CHECK(istwfk == 1, "istwfk /=1 not coded")
   call add_vnlr_commutator(KBgrad_k,cryst,npw,nspinor,ug1,ug2,ihr_comm)
 end if

end function nc_ihr_comm
!!***

!----------------------------------------------------------------------

!!****f* m_commutator_vkbr/ccgradvnl_ylm
!! NAME
!! ccgradvnl_ylm
!!
!! FUNCTION
!!  Compute Vnl(K) and grad_K Vnl(K) three reciprocal lattice units components
!!  using spherical harmonics instead of Legendre polynomials
!!  Needed for chi0(q=0)
!!
!! INPUTS
!!  cryst<crystal_t>=Unit cell and symmetries
!!  psps<pseudopotential_type>Structure gathering info on the pseudopotentials.
!!  npw=number of planewaves for wavefunctions
!!  gvec(3,npw)=integer coordinates of each plane wave in reciprocal space
!!  kpoint(3)=K-point in reduced coordinates.
!!  vkbsign(mpsang,ntypat)=sign of each KB dyadic product
!!  vkb(npw,mpsang,ntypat)=KB projector function
!!  vkbd(npw,mpsang,ntypat)=derivative of the KB projector function in reciprocal space
!!
!! OUTPUT
!!  l_fnl(npw,mpsang*2,natom),
!!  l_fnld(3,npw,mpsang*2,natom)
!!
!! NOTES
!!  Subroutine taken from the EXC code  
!!  All the calculations are done in double precision, but the output arrays l_fnl and l_fnld 
!!  are in single precision, should use double precision after modification of the other subroutines 
!!
!!  TODO
!!  the subroutine does not work wity pseudo with more that one projector per angular state 
!!
!! PARENTS
!!      m_commutator_vkbr
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine ccgradvnl_ylm(cryst,psps,npw,gvec,kpoint,vkbsign,vkb,vkbd,l_fnl,l_fnld)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ccgradvnl_ylm'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw
 type(crystal_t),intent(in) :: cryst
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: gvec(3,npw)
 real(dp),intent(in) :: kpoint(3)
 real(dp),intent(in) :: vkb(npw,psps%mpsang,cryst%ntypat)
 real(dp),intent(in) :: vkbd(npw,psps%mpsang,cryst%ntypat) 
 real(dp),intent(in) :: vkbsign(psps%mpsang,Cryst%ntypat)
 complex(gwpc),intent(out) :: l_fnl(npw,psps%mpsang**2,cryst%natom)
 complex(gwpc),intent(out) :: l_fnld(3,npw,psps%mpsang**2,cryst%natom)

!Local variables-------------------------------
!scalars
 integer,parameter :: nlx=4
 integer :: i,iat,ig,il,im,iml,ityp,lmax
 real(dp),parameter :: ppad=tol8
 real(dp) :: cosphi,costh,factor,mkg,mkg2,sinphi,sinth,sq,xdotg
 complex(dpc) :: dphi,dth,sfac
 character(len=500) :: msg
!arrays
 real(dp) :: gcart(3),kcart(3),kg(3)
 real(dp) :: b1(3),b2(3),b3(3),a1(3),a2(3),a3(3)
 complex(dpc) :: dylmcart(3),dylmcrys(3),gradphi(3),gradth(3)
!************************************************************************

 DBG_ENTER("COLL")

 lmax = psps%mpsang
 if (psps%mpsang > nlx) then
   write(msg,'(3a)')&
    'Number of angular momentum components bigger than programmed.',ch10,&
    'Taking into account only s p d f ' 
   MSG_WARNING(msg)
   lmax = nlx
 end if

 a1=cryst%rprimd(:,1); b1=two_pi*Cryst%gprimd(:,1)
 a2=cryst%rprimd(:,2); b2=two_pi*Cryst%gprimd(:,2)
 a3=cryst%rprimd(:,3); b3=two_pi*Cryst%gprimd(:,3)

 ! Calculate Kleiman-Bylander factor and first derivative.
 l_fnl=czero_gw; l_fnld=czero_gw

 do ig=1,npw
   ! Get kcart = k+G in Cartesian coordinates.
   kg(:)= kpoint(:) + REAL(gvec(:,ig))
   kcart(:) = kg(1)*b1(:) + kg(2)*b2(:) + kg(3)*b3(:)
   ! Solve the problem with sinth=0. or sinphi=0
   if (ABS(kcart(2))<ppad) kcart(2) = kcart(2) + ppad

   mkg2 = kcart(1)**2+kcart(2)**2+kcart(3)**2
   mkg = SQRT(mkg2)
   ! The next to solve the problem with k=Gamma.
   !if (mkg < 0.0001) cycle

   sq=SQRT(kcart(1)**2+kcart(2)**2)

   gcart(:)=  REAL(gvec(1,ig))*b1(:)&
&            +REAL(gvec(2,ig))*b2(:)&
&            +REAL(gvec(3,ig))*b3(:)

   ! === Calculate spherical coordinates (th,phi) ===
   costh = kcart(3)/mkg
   sinth = sq/mkg
   cosphi= kcart(1)/sq
   sinphi= kcart(2)/sq
   
   gradth(1)  = kcart(1)*kcart(3)/mkg**3/sinth
   gradth(2)  = kcart(2)*kcart(3)/mkg**3/sinth
   gradth(3)  = -(one/mkg-kcart(3)**2/mkg**3)/sinth
   gradphi(1) = -(one/sq - kcart(1)**2/sq**3)/sinphi
   gradphi(2) = kcart(2)*kcart(1)/sq**3/sinphi
   gradphi(3) = czero
   
   do iat=1,cryst%natom
     ityp = cryst%typat(iat)
     xdotg = gcart(1)*cryst%xcart(1,iat)+gcart(2)*Cryst%xcart(2,iat)+gcart(3)*Cryst%xcart(3,iat)
     ! Remember that in the GW code the reciprocal vectors 
     ! are defined such as a_i*b_j = 2pi delta_ij, no need to introduce 2pi
     sfac=CMPLX(COS(xdotg),SIN(xdotg)) 

     do il=1,lmax
       factor = SQRT(four_pi/REAL(2*(il-1)+1))
       do im=1,2*(il-1)+1
         ! Index of im and il
         iml=im+(il-1)*(il-1)

         ! Calculate the first KB factor, note that l_fnl is simple precision complex
         l_fnl(ig,iml,iat) = factor*sfac*ylmc(il-1,im-il,kcart) * vkb(ig,il,ityp) * vkbsign(il,ityp)

         ! Calculate the second KB factor (involving first derivatives)
         ! dYlm/dK = dYlm/dth * grad_K th + dYlm/dphi + grad_K phi
         call ylmcd(il-1,im-il,kcart,dth,dphi)
         dylmcart(:) = dth*gradth(:) + dphi*gradphi(:)

         ! Cartesian to crystallographic axis
         ! Notice: a bug was discovered by Marco Cazzaniga, december 2009
         ! the transformation matrix A=(a1,a2,a3) must act *on its left* on the
         ! covariant vector dylmcart (a *row* vector). The previous implementation assumed A 
         ! acting on its right on a column vector, yielding wrong results for the (small)
         ! non local contributions to the spectra, such as a spurious anisotropy in isotropic systems.
         ! This is the correct version:
         dylmcrys(1) = (a1(1)*dylmcart(1)+a1(2)*dylmcart(2)+a1(3)*dylmcart(3))/(two_pi)
         dylmcrys(2) = (a2(1)*dylmcart(1)+a2(2)*dylmcart(2)+a2(3)*dylmcart(3))/(two_pi)
         dylmcrys(3) = (a3(1)*dylmcart(1)+a3(2)*dylmcart(2)+a3(3)*dylmcart(3))/(two_pi)

         ! Note that l_fnld is simple precision complex, it could be possible to use double precision
         do i=1,3
           l_fnld(i,ig,iml,iat) = factor*sfac* &
            ( kg(i)/mkg*ylmc(il-1,im-il,kcart)*vkbd(ig,il,ityp) + dylmcrys(i)*vkb(ig,il,ityp) )
         end do 

       end do !im
     end do !il
   end do !iat
 end do !ig

 DBG_EXIT("COLL")

end subroutine ccgradvnl_ylm
!!***

!----------------------------------------------------------------------

END MODULE m_commutator_vkbr
!!***

