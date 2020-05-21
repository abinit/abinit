!!****m* ABINIT/m_vkbr
!! NAME
!!  m_vkbr
!!
!! FUNCTION
!!  This module provides objects and methods used to calculate the matrix elements
!!  of the commutator [H,r] needed for the correct treatment of the optical limit q-->0
!!  in the matrix elements <k-q,b1|e^{-iqr}|k,b2> when non-local pseudopotentials are used.
!!
!! NOTES
!!  This module is deprecated. Use ddkop_t in m_ddk.F90
!!
!! COPYRIGHT
!! Copyright (C) 2008-2020 ABINIT group (MG, FB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_vkbr

 use defs_basis
 use m_hide_blas
 use m_errors
 use m_abicore

 use defs_datatypes,  only : pseudopotential_type
 use m_gwdefs,        only : czero_gw
 use m_fstrings,      only : sjoin, itoa
 use m_paw_sphharm,   only : ylmc, ylmcd
 use m_geometry,      only : normv
 use m_crystal,       only : crystal_t
 use m_kg,            only : mkkin
 use m_mkffnl,        only : mkffnl

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_vkbr/vkbr_t
!! NAME
!!
!! FUNCTION
!!  Matrix elements in |k+G> space needed for the
!!  evaluation of the matrix elements of the commutator [Vnl,r] for the
!!  optical limit in <kb1|e^{-iqr}|kb2>.
!!
!! SOURCE

 type,public :: vkbr_t

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
  ! 0 to exclude commutator, 2 to include it

  real(dp) :: kpoint(3)
  ! The k-point in reduced coordinates.

  complex(gwpc),allocatable :: fnl(:,:,:,:)
  ! fnl(npw,mpsang**2,mproj,natom)

  complex(gwpc),allocatable :: fnld(:,:,:,:,:)
  ! fnld(3,npw,mpsang**2,mproj,natom)

 end type vkbr_t

 public :: vkbr_init       ! vkbr_t Constructor
 public :: vkbr_free       ! Free memory
 public :: nc_ihr_comm     ! Compute matrix elements of the commutator i[H,r] for NC pseudos
 public :: calc_vkb        ! Kleynman-Bylander form factors and derivatives.
!!***

 interface vkbr_free
   module procedure vkbr_free_0D
   module procedure vkbr_free_1D
 end interface vkbr_free

CONTAINS  !========================================================================================

!----------------------------------------------------------------------

!!****f* m_vkbr/vkbr_init
!! NAME
!!  vkbr_init
!!
!! FUNCTION
!!  Creation method the the vkbr_t structures datatype.
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
!!  vkbr<vkbr_t>=Structure containing arrays needed for calculating <\psi_1|[Vnl,r]\psi_2>.
!!    Completely initialized in output.
!!
!! PARENTS
!!      calc_optical_mels,cchi0q0,cchi0q0_intraband
!!
!! CHILDREN
!!      ylmcd
!!
!! SOURCE

subroutine vkbr_init(vkbr,cryst,psps,inclvkb,istwfk,npw,kpoint,gvec)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw,inclvkb,istwfk
 type(crystal_t),intent(in) :: cryst
 type(vkbr_t),intent(inout) :: vkbr
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

 !@vkbr_t
 vkbr%istwfk = istwfk
 vkbr%ntypat = cryst%ntypat
 vkbr%natom = cryst%natom
 vkbr%mpsang  = psps%mpsang
 vkbr%npw = npw
 vkbr%inclvkb = inclvkb
 vkbr%kpoint = kpoint

 ! Calculate KB form factors and derivatives.
 ! The arrays are allocated with lnmax to support pseudos with more than projector.
 ! Note that lnmax takes into account lloc hence arrays are in packed form and one should be
 ! accessed with the indices provided by indlmn.
 ! TODO: they should be calculated on-the-fly using calc_vkb
 !       For the moment, we opt for a quick an dirty implementation.

 ABI_MALLOC(vkbsign, (psps%lnmax, cryst%ntypat))
 ABI_MALLOC(vkb, (npw, psps%lnmax, cryst%ntypat))
 ABI_MALLOC(vkbd, (npw, psps%lnmax, cryst%ntypat))
 call calc_vkb(cryst,psps,kpoint,npw,npw,gvec,vkbsign,vkb,vkbd)

 select case (inclvkb)
 case (2)
   ! Complex spherical harmonics (CPU and mem \propto npw).
   write(msg,'(a,f12.1)')'out-of-memory in fnl; Mb= ',one*npw*psps%mpsang**2*psps%mproj*cryst%natom*2*gwpc*b2Mb
   ABI_STAT_MALLOC(vkbr%fnl,(npw,psps%mpsang**2,psps%mproj,cryst%natom), ierr)
   ABI_CHECK(ierr==0, msg)

   write(msg,'(a,f12.1)')'out-of-memory in fnld; Mb= ',three*npw*psps%mpsang**2*psps%mproj*cryst%natom*2*gwpc*b2Mb
   ABI_STAT_MALLOC(vkbr%fnld,(3,npw,psps%mpsang**2,psps%mproj,cryst%natom), ierr)
   ABI_CHECK(ierr==0, msg)

   call ccgradvnl_ylm(cryst,psps,npw,gvec,kpoint,vkbsign,vkb,vkbd,vkbr%fnl,vkbr%fnld)

 case default
   MSG_ERROR(sjoin("Wrong inclvkb= ",itoa(inclvkb)))
 end select

 ABI_FREE(vkbsign)
 ABI_FREE(vkb)
 ABI_FREE(vkbd)

end subroutine vkbr_init
!!***

!----------------------------------------------------------------------

!!****f* m_vkbr/vkbr_free_0D
!! NAME
!!  vkbr_free_0D
!!
!! FUNCTION
!!  Free all memory allocated in a structure of type vkbr_t
!!
!! PARENTS
!!      m_vkbr
!!
!! CHILDREN
!!      ylmcd
!!
!! SOURCE

subroutine vkbr_free_0D(vkbr)

!Arguments ------------------------------------
!scalars
 type(vkbr_t),intent(inout) :: vkbr

!************************************************************************

!complex
 ABI_SFREE(vkbr%fnl)
 ABI_SFREE(vkbr%fnld)

end subroutine vkbr_free_0D
!!***

!----------------------------------------------------------------------

!!****f* m_vkbr/vkbr_free_1D
!! NAME
!!  vkbr_free_1D
!!
!! FUNCTION
!!  Free all memory allocated in a structure of type vkbr_t
!!
!! PARENTS
!!
!! CHILDREN
!!      ylmcd
!!
!! SOURCE

subroutine vkbr_free_1D(vkbr)

!Arguments ------------------------------------
!arrays
 type(vkbr_t),intent(inout) :: vkbr(:)

!Local variables ------------------------------
!scalars
 integer :: ii

!************************************************************************

 do ii=1,SIZE(vkbr)
   call vkbr_free_0D(vkbr(ii))
 end do

end subroutine vkbr_free_1D
!!***

!----------------------------------------------------------------------

!!****f* m_vkbr/add_vnlr_commutator
!! NAME
!!  add_vnlr_commutator
!!
!! FUNCTION
!!  Calculate the matrix elements of the dipole operator <phi1|r|phi2>.
!!  For norm conserving potentials the commutator [Vnl,r] is included according to inclvkb.
!!
!! INPUTS
!!  vkbr<vkbr_t>
!!  cryst<crystal_t>=Datatype gathering info on the crystal structure.
!!  psps<pseudopotential_type>Structure gathering info on the pseudopotentials.
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
!!      m_vkbr
!!
!! CHILDREN
!!      ylmcd
!!
!! SOURCE

subroutine add_vnlr_commutator(vkbr,cryst,psps,npw,nspinor,ug1,ug2,rhotwx)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw,nspinor
 type(vkbr_t),intent(in) :: vkbr
 type(crystal_t),intent(in) :: cryst
 type(pseudopotential_type),intent(in) :: psps
!arrays
 complex(gwpc),target,intent(in) :: ug1(npw*nspinor),ug2(npw*nspinor)
 complex(gwpc),intent(inout) :: rhotwx(3,nspinor**2)

!Local variables ------------------------------
!scalars
 integer :: iat,ig,ilm,itypat,nlmn,ilmn,iln0,iln,il,in,im
 complex(gwpc) :: cta1,cta4
!arrays
 complex(gwpc) :: dum(3),cta2(3),cta3(3),gamma_term(3)

!************************************************************************

 ABI_CHECK(nspinor == 1, "inclvkb > 0 with nspinor == 2 is not coded")

 ! Adding term i <c,k|[Vnl,r]|v,k> ===
 select case (vkbr%inclvkb)
 case (2)
  ! Complex spherical harmonics (much faster!).
  dum=czero_gw; gamma_term=czero

  do iat=1,vkbr%natom
    itypat = cryst%typat(iat)
    nlmn = count(psps%indlmn(3,:,itypat) > 0)
    iln0 = 0
    do ilmn=1,nlmn
      il = 1 + psps%indlmn(1,ilmn,itypat)
      in = psps%indlmn(3,ilmn,itypat)
      iln = psps%indlmn(5,ilmn,itypat)
      if (iln <= iln0) cycle
      iln0 = iln
      !if (indlmn(6,ilmn,itypat) /= 1 .or. vkbsign(iln,itypat) == zero) cycle
      !in = 1
      do im=1,2*(il-1)+1
        ! Index of im and il
        ilm = im + (il-1)*(il-1)
        cta1 = czero_gw; cta2(:) = czero_gw
        cta4 = czero_gw; cta3(:) = czero_gw
        do ig=1,npw
          ! Here we take advantage of the property Y_(l-m)= (-i)^m Y_lm^*.
          cta1   = cta1    + ug1(ig) * vkbr%fnl (ig,ilm,in,iat)
          cta2(:)= cta2(:) + ug2(ig) * vkbr%fnld(:,ig,ilm,in,iat)
          cta3(:)= cta3(:) + ug1(ig) * vkbr%fnld(:,ig,ilm,in,iat)
          cta4   = cta4    + ug2(ig) * vkbr%fnl (ig,ilm,in,iat)
          if (ig==1) gamma_term = gamma_term + CONJG(cta1)*cta2(:) +CONJG(cta3(:))*cta4
        end do
        dum(:)= dum(:) + CONJG(cta1)*cta2(:) + CONJG(cta3(:))*cta4
      end do

    end do
  end do

  if (vkbr%istwfk>1) then
    dum = two * j_dpc * AIMAG(dum); if (vkbr%istwfk==2) dum = dum - j_dpc * AIMAG(gamma_term)
  end if
  rhotwx(:,1) = rhotwx(:,1) + dum(:)

 case default
   MSG_ERROR(sjoin("Wrong inclvkb:", itoa(vkbr%inclvkb)))
 end select

end subroutine add_vnlr_commutator
!!***

!----------------------------------------------------------------------

!!****f* m_vkbr/calc_vkb
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
!!  vkb (npw_k, %lnmax, %ntypat)=KB form factors.
!!  vkbd(npw_k, %lnmax, %ntypat)=KB form factor derivatives.
!!  vkbsign(%lnmax, %ntypat)   =KS dyadic sign.
!!
!! TODO
!!  SOC not implemented.
!!
!! PARENTS
!!      m_vkbr
!!
!! CHILDREN
!!      ylmcd
!!
!! SOURCE

subroutine calc_vkb(cryst,psps,kpoint,npw_k,mpw,kg_k,vkbsign,vkb,vkbd)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_k, mpw
 type(crystal_t),intent(in) :: cryst
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: kg_k(3,npw_k)
 real(dp),intent(in) :: kpoint(3)
 real(dp),intent(out) :: vkb (mpw,psps%lnmax,psps%ntypat)
 real(dp),intent(out) :: vkbd(mpw,psps%lnmax,psps%ntypat)
 real(dp),intent(out) :: vkbsign(psps%lnmax,psps%ntypat)

!Local variables ------------------------------
!scalars
 integer :: dimffnl,ider,idir,itypat,nkpg,in,il,ilmn,ig,iln,iln0,nlmn
 real(dp) :: effmass_free,ecutsm,ecut
!arrays
 real(dp),allocatable :: ffnl(:,:,:,:),kpg_dum(:,:),modkplusg(:),ylm_gr(:,:,:),ylm_k(:,:)

! *************************************************************************

 DBG_ENTER("COLL")
 ABI_CHECK(psps%usepaw==0, "You should not be here!")
 ABI_CHECK(psps%useylm==0, "useylm/=0 not considered!")

 ! Compute KB dyadic sign.
 vkbsign=zero
 do itypat=1,psps%ntypat
   iln0 = 0
   nlmn = count(psps%indlmn(3,:,itypat) > 0)
   do ilmn=1,nlmn
     iln = psps%indlmn(5,ilmn,itypat)
     if (iln <= iln0) cycle
     iln0 = iln
     if (abs(psps%ekb(iln,itypat)) > 1.0d-10) vkbsign(iln,itypat) = dsign(one, psps%ekb(iln,itypat))
   end do
 end do

 ! Allocate KB form factor and derivative wrt k+G
 ! Here we do not use correct ordering for dimensions
 idir=0; nkpg=0; ider=1; dimffnl=2 ! To retrieve the first derivative.

 ! Quantities used only if useylm==1
 ABI_MALLOC(ylm_k, (npw_k, psps%mpsang**2*Psps%useylm))
 ABI_MALLOC(ylm_gr, (npw_k, 3+6*(ider/2),psps%mpsang**2*Psps%useylm))
 ABI_MALLOC(kpg_dum, (npw_k, nkpg))
 ABI_MALLOC(ffnl, (npw_k,dimffnl, psps%lmnmax, psps%ntypat))

 call mkffnl(psps%dimekb,dimffnl,Psps%ekb,ffnl,Psps%ffspl,cryst%gmet,cryst%gprimd,ider,idir,Psps%indlmn,&
   kg_k,kpg_dum,kpoint,psps%lmnmax,Psps%lnmax,Psps%mpsang,Psps%mqgrid_ff,nkpg,npw_k,&
   psps%ntypat,Psps%pspso,Psps%qgrid_ff,cryst%rmet,Psps%usepaw,Psps%useylm,ylm_k,ylm_gr)

 ABI_FREE(ylm_k)
 ABI_FREE(ylm_gr)
 ABI_FREE(kpg_dum)

 ABI_MALLOC(modkplusg, (npw_k))
 effmass_free = one; ecutsm = zero; ecut = huge(one)
 call mkkin(ecut,ecutsm,effmass_free,cryst%gmet,kg_k,modkplusg,kpoint,npw_k,0,0)
 modkplusg(:) = SQRT(half/pi**2*modkplusg(:))
 modkplusg(:) = MAX(modkplusg(:),tol10)

 ! Calculate matrix elements.
 vkb=zero; vkbd=zero

 do itypat=1,psps%ntypat
   iln0 = 0
   nlmn = count(psps%indlmn(3,:,itypat) > 0)
   do ilmn=1,nlmn
     il = 1 + psps%indlmn(1,ilmn,itypat)
     in = psps%indlmn(3,ilmn,itypat)
     iln = psps%indlmn(5,ilmn,itypat)
     !write(*,*)ilmn, iln, il, in
     if (iln <= iln0) cycle
     iln0 = iln
     !if (vkbsign(iln,itypat) == zero) cycle
     if (ABS(psps%ekb(iln,itypat)) > 1.0d-10) then
       ABI_CHECK(iln == ilmn, "iln != ilmn")
       !ABI_CHECK(il == iln, "il != iln")
       if (il==1) then
         vkb (1:npw_k,iln,itypat) = ffnl(:,1,iln,itypat)
         vkbd(1:npw_k,iln,itypat) = ffnl(:,2,iln,itypat)*modkplusg(:)/two_pi
       else if (il==2) then
         vkb(1:npw_k,iln,itypat)  = ffnl(:,1,iln,itypat)*modkplusg(:)
         do ig=1,npw_k
           vkbd(ig,iln,itypat) = ((ffnl(ig,2,iln,itypat)*modkplusg(ig)*modkplusg(ig))+&
            ffnl(ig,1,iln,itypat) )/two_pi
         end do
       else if (il==3) then
         vkb (1:npw_k,iln,itypat) =  ffnl(:,1,iln,itypat)*modkplusg(:)**2
         vkbd(1:npw_k,iln,itypat) = (ffnl(:,2,iln,itypat)*modkplusg(:)**3+&
          2*ffnl(:,1,iln,itypat)*modkplusg(:) )/two_pi
       else if (il==4) then
         vkb (1:npw_k,iln,itypat) =  ffnl(:,1,iln,itypat)*modkplusg(:)**3
         vkbd(1:npw_k,iln,itypat) = (ffnl(:,2,iln,itypat)*modkplusg(:)**4+&
          3*ffnl(:,1,iln,itypat)*modkplusg(:)**2 )/two_pi
       end if
       vkb (:,iln,itypat) = SQRT(4*pi/cryst%ucvol*(2*il-1)*ABS(psps%ekb(iln,itypat)))*vkb (:,iln,itypat)
       vkbd(:,iln,itypat) = SQRT(4*pi/cryst%ucvol*(2*il-1)*ABS(psps%ekb(iln,itypat)))*vkbd(:,iln,itypat)
     end if
   end do
 end do

 ABI_FREE(ffnl)
 ABI_FREE(modkplusg)

 DBG_EXIT("COLL")

end subroutine calc_vkb
!!***

!----------------------------------------------------------------------

!!****f* m_vkbr/nc_ihr_comm
!! NAME
!!  nc_pwihr_comm
!!
!! FUNCTION
!!  Calculate the matrix elements of the commutator i[H,r]
!!  For norm conserving potentials the commutator i[Vnl,r] is included depending on inclvkb.
!!
!! INPUTS
!!  vkbr<vkbr_t>
!!  cryst<crystal_t>=Unit cell and symmetries
!!  psps<pseudopotential_type>Structure gathering info on the pseudopotentials.
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

function nc_ihr_comm(vkbr,cryst,psps,npw,nspinor,istwfk,inclvkb,kpoint,ug1,ug2,gvec) result(ihr_comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw,nspinor,inclvkb,istwfk
 type(vkbr_t),intent(in) :: vkbr
 type(crystal_t),intent(in) :: cryst
 type(Pseudopotential_type),intent(in) :: psps
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
 ! not in PAW unless even the AE Hamiltonian in non-local e.g. DFT+U or LEXX.

 ! -i <c,k|\nabla_r|v,k> in reduced coordinates is always included.
 ! -i <c,k|\nabla_r|v,k> = \sum_G u_{ck}^*(G) [k+G] u_{vk}(G)
 ! Note that here we assume c/=v, moreover the ug are supposed to be orthonormal and
 ! hence k+G can be replaced by G.
 ! HM 03/08/2018: we need band velocities so we don't assume c/=v anymore and we use k+G.

 spinorwf_pad = RESHAPE([0, 0, npw, npw, 0, npw, npw, 0], [2, 4])
 ihr_comm = czero

 ! -i <c,k|\nabla_r|v,k> in reduced coordinates.
 ! This term is spin diagonal if nspinor == 2
 if (istwfk == 1) then
   do iab=1,nspinor
     spad1 = spinorwf_pad(1,iab); spad2 = spinorwf_pad(2,iab)
     do ig=1,npw
       c_tmp = GWPC_CONJG(ug1(ig+spad1)) * ug2(ig+spad2)
       ihr_comm(:,iab) = ihr_comm(:,iab) + c_tmp * (kpoint + gvec(:,ig))
     end do
   end do
 else
   ! Symmetrized expression: \sum_G  (k+G) 2i Im [ u_a^*(G) u_b(G) ]. (k0,G0) term is null.
   ABI_CHECK(nspinor == 1, "nspinor != 1")
   do ig=1,npw
     c_tmp = GWPC_CONJG(ug1(ig)) * ug2(ig)
     ihr_comm(:,1) = ihr_comm(:,1) + two*j_dpc * AIMAG(c_tmp) * (kpoint + gvec(:,ig))
   end do
 end if

 ! Add second term $i <c,k|[Vnl,r]|v,k> $ in reduced cordinates.
 if (inclvkb /= 0) then
   ABI_CHECK(istwfk == vkbr%istwfk, "input istwfk /= vkbr%istwfk")
   call add_vnlr_commutator(vkbr,cryst,psps,npw,nspinor,ug1,ug2,ihr_comm)
 end if

end function nc_ihr_comm
!!***

!----------------------------------------------------------------------

!!****f* m_vkbr/ccgradvnl_ylm
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
!!  vkbsign(lnmax,ntypat)=sign of each KB dyadic product
!!  vkb(npw,lnmax,ntypat)=KB projector function
!!  vkbd(npw,lnmax,ntypat)=derivative of the KB projector function in reciprocal space
!!
!! OUTPUT
!!  fnl(npw,mpsang*2,natom),
!!  fnld(3,npw,mpsang*2,natom)
!!
!! NOTES
!!  Subroutine taken from the EXC code
!!  All the calculations are done in double precision, but the output arrays fnl and fnld
!!  are in single precision, should use double precision after modification of the other subroutines
!!
!! PARENTS
!!      m_vkbr
!!
!! CHILDREN
!!      ylmcd
!!
!! SOURCE

subroutine ccgradvnl_ylm(cryst,psps,npw,gvec,kpoint,vkbsign,vkb,vkbd,fnl,fnld)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw
 type(crystal_t),intent(in) :: cryst
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: gvec(3,npw)
 real(dp),intent(in) :: kpoint(3)
 real(dp),intent(in) :: vkb(npw,psps%lnmax,cryst%ntypat)
 real(dp),intent(in) :: vkbd(npw,psps%lnmax,cryst%ntypat)
 real(dp),intent(in) :: vkbsign(psps%lnmax,cryst%ntypat)
 complex(gwpc),intent(out) :: fnl(npw,psps%mpsang**2,psps%mproj,cryst%natom)
 complex(gwpc),intent(out) :: fnld(3,npw,psps%mpsang**2,psps%mproj,cryst%natom)

!Local variables-------------------------------
!scalars
 integer :: ii,iat,ig,il,im,ilm,itypat,nlmn,iln0,iln,ilmn,in
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

 if (psps%mpsang > 4) then
   write(msg,'(3a)')&
    'Number of angular momentum components bigger than programmed.',ch10,&
    'Taking into account only s p d f '
   MSG_ERROR(msg)
 end if

 a1=cryst%rprimd(:,1); b1=two_pi*Cryst%gprimd(:,1)
 a2=cryst%rprimd(:,2); b2=two_pi*Cryst%gprimd(:,2)
 a3=cryst%rprimd(:,3); b3=two_pi*Cryst%gprimd(:,3)

 ! Calculate Kleiman-Bylander factor and first derivative.
 fnl=czero_gw; fnld=czero_gw

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

   ! Calculate spherical coordinates (th, phi).
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
     itypat = cryst%typat(iat)
     xdotg = gcart(1)*cryst%xcart(1,iat)+gcart(2)*Cryst%xcart(2,iat)+gcart(3)*Cryst%xcart(3,iat)
     ! Remember that in the GW code the reciprocal vectors
     ! are defined such as a_i*b_j = 2pi delta_ij, no need to introduce 2pi
     sfac=CMPLX(COS(xdotg), SIN(xdotg), kind=dpc)

     iln0 = 0
     nlmn = count(psps%indlmn(3,:,itypat) > 0)
     do ilmn=1,nlmn
       il = 1 + psps%indlmn(1,ilmn,itypat)
       in = psps%indlmn(3,ilmn,itypat)
       iln = psps%indlmn(5,ilmn,itypat)
       ! spin = 1 if scalar term (spin diagonal), 2 if SOC term.
       !spin = psps%indlmn(6, ilmn, itypat)
       if (iln <= iln0) cycle
       iln0 = iln
       if (vkbsign(iln,itypat) == zero) cycle
       !if (spin /= 1 .or. vkbsign(iln,itypat) == zero) cycle
       factor = SQRT(four_pi/REAL(2*(il-1)+1))
       do im=1,2*(il-1)+1
         ! Index of im and il
         ilm = im + (il-1)*(il-1)

         ! Calculate the first KB factor, note that fnl is simple precision complex
         fnl(ig,ilm,in,iat) = factor*sfac*ylmc(il-1,im-il,kcart) * vkb(ig,iln,itypat) * vkbsign(iln,itypat)

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

         ! Note that fnld is simple precision complex, it could be possible to use double precision
         do ii=1,3
           fnld(ii,ig,ilm,in,iat) = factor*sfac* &
            ( kg(ii)/mkg*ylmc(il-1,im-il,kcart)*vkbd(ig,iln,itypat) + dylmcrys(ii)*vkb(ig,iln,itypat) )
         end do

       end do !im
     end do !il
   end do !iat
 end do !ig

 DBG_EXIT("COLL")

end subroutine ccgradvnl_ylm
!!***

!----------------------------------------------------------------------

END MODULE m_vkbr
!!***
