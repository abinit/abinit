!!****m* ABINIT/m_paw_hr
!! NAME
!!  m_paw_hr
!!
!! FUNCTION
!!  This module provides objects and methods to calculate the matrix elements 
!!  of the commutator PAW [H,r] needed for the correct treatment of the optical limit q-->0
!!  in the matrix elements <k-q,b1|e^{-iqr}|k,b2>. As PAW is a full potential method 
!!  the commutator reduces to the contribution given by the velocity operator. 
!!  However, when the all-electron Hamiltonian is non-local (e.g. LDA+U or 
!!  LEXX) additional on-site terms have to be considered in the calculation of the
!!  matrix elements of [H.r].
!!
!! COPYRIGHT
!! Copyright (C) 2008-2020 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_paw_hr

 use defs_basis
 use m_abicore
 use m_errors

 use m_crystal,        only : crystal_t
 use m_pawang,         only : pawang_type
 use m_pawrad,         only : pawrad_type, simp_gen
 use m_pawtab,         only : pawtab_type
 use m_paw_ij,         only : paw_ij_type
 use m_pawfgrtab,      only : pawfgrtab_type
 use m_pawcprj,        only : pawcprj_type
 use m_pawdij,         only : pawpupot
 use m_paw_pwaves_lmn, only : paw_pwaves_lmn_t
 
 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_paw_hr/pawhur_t
!! NAME
!!  pawhur_t
!!
!! FUNCTION
!!  The pawhur_t data type stores basic dimensions and quantities 
!!  used in the GW part for the treatment of the non-analytic behavior of the 
!!  heads and wings of the irreducible polarizability in the long wave-length limit (i.e. q-->0).
!!  Note that, within the PAW formalism, a standard KS Hamiltonian has a semi-local contribution 
!!  arising from the kinetic operator (if we work in the AE representation). 
!!  When LDA+U is used, a fully non-local term is added to the Hamiltonian whose commutator with the position operator 
!!  has to be considered during the calculation of the heads and wings of the polarizability in the optical limit
!!
!! SOURCE

 type,public :: pawhur_t

  integer :: lmn_size
  integer :: lmn2_size
  integer :: nsppol
  !integer :: nsel

  integer,allocatable :: ij_select(:,:,:)
  ! ijselect(lmn_size,lmn_size,nsppol) 
  ! Selection rules of ij matrix elements
  ! Do not take into account selection on x-y-x for the time being.

  real(dp),allocatable :: commutator(:,:,:)
  ! commutator(3,nsel,nsppol)
 end type pawhur_t

 public ::  pawhur_init          ! Init object
 public ::  pawhur_free          ! Deallocate memory
 public ::  paw_ihr              
 public ::  paw_cross_ihr_comm

!!***

CONTAINS  !========================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_paw_hr/pawhur_free
!! NAME
!! pawhur_free
!!
!! FUNCTION
!!  Deallocate memory 
!!
!! PARENTS
!!      bethe_salpeter,cchi0q0,cchi0q0_intraband
!!
!! CHILDREN
!!      simp_gen
!!
!! SOURCE

subroutine pawhur_free(Hur)

 implicit none

!Arguments ------------------------------------
 type(pawhur_t),intent(inout) :: Hur(:)

!Local variables-------------------------------
 integer :: iat
! *************************************************************************

 do iat=1,SIZE(Hur)
   if (allocated(Hur(iat)%ij_select)) then
     ABI_FREE(Hur(iat)%ij_select)
   end if
   if (allocated(Hur(iat)%commutator)) then
     ABI_FREE(Hur(iat)%commutator)
   end if
 end do

end subroutine pawhur_free
!!***

!----------------------------------------------------------------------

!!****f* m_paw_hr/paw_ihr
!! NAME
!! paw_ihr
!!
!! FUNCTION
!!  Calculate the PAW onsite contribution to the matrix elements of the i\nabla operator.
!!  in cartesian coordinates. Take also into account the contribution arising from the U 
!!  part of the Hamiltonian (if any)
!!
!! INPUTS
!!  isppol=Spin index.
!!  nspinor=Number of spinori components.
!!  npw=Number of planewaves for this k-point.
!!  istwfk=Storage mode for the wavefunctions.
!!  kpoint(3)=k-point in reduced coordinates.
!!  Cryst<crystal_t>=Info on the crystal structure.
!!    %natom=Number of atoms in unit cell
!!    %typat(natom)
!!  Pawtab(ntypat)=Only for PAW, TABulated data initialized at start
!!    %lmn_size Number of (l,m,n) elements for the paw basis
!!    %nabla_ij(3,lmn_size,lmn_size)) Onsite contribution
!!      <phi_i|nabla|phi_j>-<tphi_i|nabla|tphi_j> for each type
!!  ug1(nspinor*npwwfn)=Left wavefunction.
!!  ug2(nspinor*npwwfn)=Right wavefunction
!!  HUr(natom)=Commutator of the LDA+U part of the Hamiltonian with the position operator.
!!  Cprj_kb1(natom,nspinor),Cprj_kb2(natom,nspinor) <type(pawcprj_type)>=
!!   projected input wave functions <Proj_i|Cnk> with all NL projectors corresponding to 
!!   wavefunctions (k,b1,s) and (k,b2,s), respectively.
!!
!! OUTPUT
!!  onsite(2,3)=Onsite contribution to  $i<ug1|\nabla|ug2>$
!!
!! PARENTS
!!      cchi0q0,debug_tools,spectra
!!
!! SOURCE

function paw_ihr(isppol,nspinor,npw,istwfk,kpoint,Cryst,Pawtab,ug1,ug2,gvec,Cprj_kb1,Cprj_kb2,HUr) result(ihr_comm)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: isppol,nspinor,npw,istwfk
 complex(gwpc) :: ihr_comm(3,nspinor**2)
 type(crystal_t),intent(in) :: Cryst
!arrays
 integer,intent(in) :: gvec(3,npw)
 real(dp),intent(in) :: kpoint(3)
 complex(gwpc),intent(in) :: ug1(nspinor*npw),ug2(nspinor*npw)
 type(Pawtab_type),target,intent(in) :: Pawtab(Cryst%ntypat)
 type(pawcprj_type),intent(in) :: Cprj_kb1(Cryst%natom,nspinor),Cprj_kb2(Cryst%natom,nspinor)
 type(pawhur_t),intent(in) :: Hur(Cryst%natom)

!Local variables-------------------------------
 integer :: iatom,itypat,lmn_size,ilmn,jlmn,isel
 integer :: ig,iab,spad1,spad2
 real(dp) :: re_p,im_p
 complex(dpc) :: ctemp
!arrays
 integer :: spinorwf_pad(2,4)
 real(dp) :: hurc_ij(3),ons_cart(2,3) !,ons_comm_red(2,3)
 real(dp) :: gspace_cart2red(3,3) !rs_cart2red(3,3),
 real(dp), ABI_CONTIGUOUS pointer :: nabla_ij(:,:,:)
 complex(gwpc) :: ihr_comm_cart(3,nspinor**2)

! *************************************************************************

 ! [H, r] = -\nabla + [V_{nl}, r] 
 ! Note that V_nl is present only if the AE-Hamiltonian is non-local e.g. LDA+U or LEXX.
 spinorwf_pad=RESHAPE((/0,0,npw,npw,0,npw,npw,0/),(/2,4/))
 ihr_comm=zero

 ! -i <c,k|\nabla_r|v,k> = \sum_G u_{ck}^*(G) [k+G] u_{vk}(G) in reduced coordinates.
 if (istwfk==1) then
   do iab=1,nspinor**2
     spad1 = spinorwf_pad(1,iab)
     spad2 = spinorwf_pad(2,iab)
     do ig=1,npw 
       ctemp = CONJG(ug1(ig+spad1)) * ug2(ig+spad2)
       ihr_comm(:,iab) = ihr_comm(:,iab) + ctemp* ( kpoint + gvec(:,ig))
     end do
   end do
 else 
   ! Symmetrized expression: \sum_G  (k+G) 2i Im [ u_a^*(G) u_b(G) ]. (k0,G0) term is null.
   do ig=1,npw 
     ctemp = CONJG(ug1(ig)) * ug2(ig)
     ihr_comm(:,1) = ihr_comm(:,1) + two*j_dpc * AIMAG(ctemp) * (kpoint + gvec(:,ig))
   end do
 end if
 !
 ! Add on-site terms.
 ons_cart=zero
 ABI_CHECK(nspinor==1,"nspinor/=1 not coded")

 do iatom=1,Cryst%natom 
   itypat=Cryst%typat(iatom)
   lmn_size=Pawtab(itypat)%lmn_size
   nabla_ij => Pawtab(itypat)%nabla_ij(:,:,:) 
   !
   !=== Unpacked loop over lmn channels ====
   do jlmn=1,lmn_size
     do ilmn=1,lmn_size
       re_p =  Cprj_kb1(iatom,1)%cp(1,ilmn)*Cprj_kb2(iatom,1)%cp(1,jlmn) &
&             +Cprj_kb1(iatom,1)%cp(2,ilmn)*Cprj_kb2(iatom,1)%cp(2,jlmn) 
 
       im_p =  Cprj_kb1(iatom,1)%cp(1,ilmn)*Cprj_kb2(iatom,1)%cp(2,jlmn) &
&             -Cprj_kb1(iatom,1)%cp(2,ilmn)*Cprj_kb2(iatom,1)%cp(1,jlmn)

       ! Onsite contribution given by -i\nabla.
       ons_cart(1,1)=ons_cart(1,1) + im_p*nabla_ij(1,ilmn,jlmn)
       ons_cart(1,2)=ons_cart(1,2) + im_p*nabla_ij(2,ilmn,jlmn)
       ons_cart(1,3)=ons_cart(1,3) + im_p*nabla_ij(3,ilmn,jlmn)

       ons_cart(2,1)=ons_cart(2,1) - re_p*nabla_ij(1,ilmn,jlmn)
       ons_cart(2,2)=ons_cart(2,2) - re_p*nabla_ij(2,ilmn,jlmn)
       ons_cart(2,3)=ons_cart(2,3) - re_p*nabla_ij(3,ilmn,jlmn)
       !
       if (Pawtab(itypat)%usepawu/=0) then ! Add i[V_u, r] 
         isel=Hur(iatom)%ij_select(ilmn,jlmn,isppol)
         if (isel>0) then
           hurc_ij(:)=Hur(iatom)%commutator(:,isel,isppol)

           ons_cart(1,1)=ons_cart(1,1) - im_p*hurc_ij(1)
           ons_cart(1,2)=ons_cart(1,2) - im_p*hurc_ij(2)
           ons_cart(1,3)=ons_cart(1,3) - im_p*hurc_ij(3)

           ons_cart(2,1)=ons_cart(2,1) + re_p*hurc_ij(1)
           ons_cart(2,2)=ons_cart(2,2) + re_p*hurc_ij(2)
           ons_cart(2,3)=ons_cart(2,3) + re_p*hurc_ij(3)
         end if
       end if

     end do !ilmn
   end do !jlmn
 end do !iatom

 ! ons_cart is in Cartesian coordinates in real space 
 ! while ihr_comm is in reduced coordinates in reciprocal space in terms of gprimd.
 !rs_cart2red = TRANSPOSE(Cryst%gprimd) ! if <r> is in terms of real space vectors
 gspace_cart2red = TRANSPOSE(Cryst%rprimd)

 !ons_comm_red(1,:)=MATMUL(rs_cart2red,ons_comm(1,:))
 !ons_comm_red(2,:)=MATMUL(rs_cart2red,ons_comm(2,:))
 !ihr_comm(:,1) = ihr_comm(:,1) + CMPLX(ons_comm_red(1,:),ons_comm_red(2,:),kind=gwpc)

 ihr_comm_cart(:,1) = two_pi*MATMUL(Cryst%gprimd,ihr_comm(:,1))
 ihr_comm_cart(:,1) = ihr_comm_cart(:,1) + CMPLX(ons_cart(1,:),ons_cart(2,:),kind=gwpc)

 ! Final result is in reduced coordinates, in terms of gprimd.
 ihr_comm(:,1) = MATMUL(gspace_cart2red, ihr_comm_cart(:,1))/two_pi

end function paw_ihr
!!***

!----------------------------------------------------------------------

!!****f* m_paw_hr/paw_cross_ihr_comm
!! NAME
!! paw_cross_ihr_comm
!!
!! FUNCTION
!!  Adds the PAW cross term contribution to the matrix elements of the  i\nabla operator.
!!  in cartesian coordinates. Should take also into account the contribution arising from the U 
!!  part of the Hamiltonian (if any)
!!
!! INPUTS
!!  ihr_comm = the commutator [H,r] evaluated between states i and j, with only the plane-wave and 
!!              the onsite parts included
!!  isppol=Spin index.
!!  nspinor=Number of spinori components.
!!  nr=Number real-space points on the fine fft grid for the ae wavefunctions
!!  kpoint(3)=k-point in reduced coordinates.
!!  Cryst<crystal_t>=Info on the crystal structure.
!!    %natom=Number of atoms in unit cell
!!    %typat(natom)
!!  Pawfgrtab(ntypat)= PAW tabulated data on the fine grid
!!    %lmn_size Number of (l,m,n) elements for the paw basis
!!    %nfgr Number of points on the fine grid
!!    %ifftsph Indexes of the fine-grid points on the fft mesh
!!  Paw_onsite(ntypat)= PAW tabulated data on the fine grid points inside the sphere
!!    %phi_gr(3,nfgr,lmn_size) gradient of phi in cartesian coordinates
!!    %tphi_gr(3,nfgr,lmn_size) gradient of tphi in cartesian coordinates
!!  ur_ae1(nr),ur_ae2(nr)=Left and right AE wavefunction.
!!  ur_ae_onsite1(nr),ur_ae_onsite2(nr)=Left and right AE onsite wavefunction.
!!  Cprj_kb1(natom,nspinor),Cprj_kb2(natom,nspinor) <type(pawcprj_type)>=
!!   projected input wave functions <Proj_i|Cnk> with all NL projectors corresponding to 
!!   wavefunctions (k,b1,s) and (k,b2,s), respectively.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  The cross-term contribution is added to the commutator 
!!
!! PARENTS
!!      cchi0q0
!!
!! CHILDREN
!!      simp_gen
!!
!! SOURCE

subroutine paw_cross_ihr_comm(ihr_comm,nspinor,nr,Cryst,Pawfgrtab,Paw_onsite,&
& ur_ae1,ur_ae2,ur_ae_onsite1,ur_ae_onsite2,Cprj_kb1,Cprj_kb2)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nspinor,nr
 type(crystal_t),intent(in) :: Cryst
!arrays
 complex(gwpc),intent(inout) :: ihr_comm(3,nspinor**2)
 complex(gwpc),intent(in) :: ur_ae1(nr),ur_ae2(nr)
 complex(gwpc),intent(in) :: ur_ae_onsite1(nr),ur_ae_onsite2(nr)
 type(pawfgrtab_type),intent(in) :: Pawfgrtab(Cryst%natom)
 type(paw_pwaves_lmn_t),intent(in) :: Paw_onsite(Cryst%natom)
 type(pawcprj_type),intent(in) :: Cprj_kb1(Cryst%natom,nspinor),Cprj_kb2(Cryst%natom,nspinor)

!Local variables-------------------------------
 integer :: iatom,lmn_size,ilmn,ifgd,ifftsph,nfgd
 complex(dpc) :: cp1, cp2
 complex(dpc) :: cross1,cross2
!arrays
 real(dp) :: gspace_cart2red(3,3)
 complex(gwpc) :: ihr_comm_cart(3,nspinor**2)
 complex(dpc) :: dphigr(3), dphigr1(3),dphigr2(3)

! *************************************************************************

 ABI_CHECK(nspinor==1,"nspinor + pawcross not implemented")

 ! [H, r] = -\nabla + [V_{nl}, r]
 ! The V_nl part, present in case of LDA+U, is omitted for the cross terms contribution
 ! Recall that delta_rho_tw_ij = (psi_i - phi_i)* (phi_j - tphi_j) + (phi_i - tphi_i)* (psi_j - phi_j)
 ihr_comm_cart(:,1) = czero

 do iatom=1,Cryst%natom
   lmn_size = Paw_onsite(iatom)%lmn_size
   nfgd = Pawfgrtab(iatom)%nfgd

   do ifgd=1,nfgd

     ifftsph = Pawfgrtab(iatom)%ifftsph(ifgd)

     cross1 = ur_ae1(ifftsph) - ur_ae_onsite1(ifftsph)
     cross2 = ur_ae2(ifftsph) - ur_ae_onsite2(ifftsph)

     do ilmn=1,lmn_size

       dphigr(1:3) = Paw_onsite(iatom)%phi_gr(1:3,ifgd,ilmn) - Paw_onsite(iatom)%tphi_gr(1:3,ifgd,ilmn)

       cp1 = CMPLX(Cprj_kb1(iatom,1)%cp(1,ilmn),Cprj_kb1(iatom,1)%cp(2,ilmn)) * sqrt(Cryst%ucvol) ! that damn magic factor
       cp2 = CMPLX(Cprj_kb2(iatom,1)%cp(1,ilmn),Cprj_kb2(iatom,1)%cp(2,ilmn)) * sqrt(Cryst%ucvol)

       dphigr1(1:3) = cp1 * dphigr(1:3)
       dphigr2(1:3) = cp2 * dphigr(1:3)

       ihr_comm_cart(1,1) = ihr_comm_cart(1,1) - j_dpc * (CONJG(cross1) * dphigr2(1) - CONJG(dphigr1(1)) * cross2) / nr
       ihr_comm_cart(2,1) = ihr_comm_cart(2,1) - j_dpc * (CONJG(cross1) * dphigr2(2) - CONJG(dphigr1(2)) * cross2) / nr
       ihr_comm_cart(3,1) = ihr_comm_cart(3,1) - j_dpc * (CONJG(cross1) * dphigr2(3) - CONJG(dphigr1(3)) * cross2) / nr

     end do
   end do
 end do

 ! Go to reduced coordinate
 gspace_cart2red = TRANSPOSE(Cryst%rprimd)
 ihr_comm(:,1) = ihr_comm(:,1) +  MATMUL(gspace_cart2red, ihr_comm_cart(:,1)) / two_pi

end subroutine paw_cross_ihr_comm
!!***

!----------------------------------------------------------------------

!!****f* m_paw_hr/pawhur_init
!! NAME
!! pawhur_init
!!
!! FUNCTION
!!  Creation method for the pawhur_t data type.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      bethe_salpeter,cchi0q0,cchi0q0_intraband
!!
!! CHILDREN
!!      simp_gen
!!
!! SOURCE

subroutine pawhur_init(hur,nsppol,pawprtvol,Cryst,Pawtab,Pawang,Pawrad,Paw_ij)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsppol,pawprtvol
 type(crystal_t),intent(in) :: Cryst
 type(Pawang_type),intent(in) :: Pawang
!arrays
 type(Pawtab_type),target,intent(in) :: Pawtab(Cryst%ntypat)
 type(Pawrad_type),intent(in) :: Pawrad(Cryst%ntypat)
 type(Paw_ij_type),intent(in) :: Paw_ij(Cryst%natom)
 type(pawhur_t),intent(inout) :: Hur(Cryst%natom)

!Local variables-------------------------------
!scalars 
 integer :: iatom,ij_idx,isel,itypat,isppol,lmn2_size_max,lmn2_size,lmn_size,lpawu
 integer :: jlmn,jl,jm,jlm,jln,k0lmn,k0lm,k0ln,ilmn,il,im,ilm,iln
 integer :: m2,m1,left_lmn,right_lmn,tot_lmn,nmax
!arrays
 integer :: nsel(3,nsppol)
 integer, ABI_CONTIGUOUS pointer :: indlmn(:,:)
 real(dp) :: sumr_ij(3)
 real(dp),allocatable :: rcart_onsite(:,:,:)
 real(dp),allocatable :: rij_tmp(:,:,:),vpawu(:,:,:,:)

! *************************************************************************

 ! Get onsite matrix elements of the position operator.
 lmn2_size_max=MAXVAL(Pawtab(:)%lmn2_size) 
 ABI_MALLOC(rcart_onsite,(3,lmn2_size_max,Cryst%natom))

 call pawr(Pawtab,Pawrad,Pawang,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xcart,lmn2_size_max,rcart_onsite)

 do iatom=1,Cryst%natom
   itypat=Cryst%typat(iatom)
   if (Pawtab(itypat)%usepawu==0) CYCLE
   lmn2_size=Pawtab(itypat)%lmn2_size
   lmn_size =Pawtab(itypat)%lmn_size
   lpawu=Pawtab(itypat)%lpawu
   Hur(iatom)%lmn2_size=lmn2_size
   Hur(iatom)%lmn_size =lmn_size
   Hur(iatom)%nsppol   =nsppol
   indlmn => Pawtab(itypat)%indlmn

   ABI_MALLOC(rij_tmp,(3,lmn_size**2,nsppol))
   rij_tmp=zero

   ! Get Vpawu^{\sigma}_{m1,m2}
   ABI_MALLOC(vpawu,(Paw_ij(iatom)%cplex_dij,2*lpawu+1,2*lpawu+1,Paw_ij(iatom)%ndij))
   call pawpupot(Paw_ij(iatom)%cplex_dij,Paw_ij(iatom)%ndij,&
&                Paw_ij(iatom)%noccmmp,Paw_ij(iatom)%nocctot,&
&                pawprtvol,Pawtab(itypat),vpawu)

   do isppol=1,nsppol ! spinor not implemented

     ! === Loop on (jl,jm,jn) channels ===
     ij_idx=0
     do jlmn=1,lmn_size
       jl =indlmn(1,jlmn)
       jm =indlmn(2,jlmn)
       jlm=indlmn(4,jlmn)
       jln=indlmn(5,jlmn)

       k0lmn=jlmn*(jlmn-1)/2 
       k0lm =jlm *(jlm -1)/2
       k0ln =jln *(jln -1)/2
       !
       ! === Loop on (il,im,in) channels === 
       ! * Looping over all ij components. Elements are not symmetric.
       do ilmn=1,lmn_size
         il =indlmn(1,ilmn)
         im =indlmn(2,ilmn)
         ilm=indlmn(4,ilmn)
         iln=indlmn(5,ilmn)

         ij_idx=ij_idx+1

         ! === Selection rules ===
         if (il/=lpawu.and.jl/=lpawu) CYCLE 

         sumr_ij(:)=zero 
         do m2=1,2*lpawu+1
           do m1=1,2*lpawu+1
             if (m1==(im-lpawu-1).and.il==lpawu) then 
               left_lmn =ilmn-(il+im+1)+m2
               right_lmn=jlmn
               if (right_lmn>=left_lmn) then 
                 tot_lmn=right_lmn*(right_lmn-1)/2 + left_lmn
               else 
                 tot_lmn=left_lmn*(left_lmn-1)/2 + right_lmn
               end if
               sumr_ij=sumr_ij+vpawu(1,m1,m2,isppol)*rcart_onsite(:,tot_lmn,iatom)
             end if

             if (m2==(jm-lpawu-1).and.jl==lpawu) then 
               left_lmn =ilmn
               right_lmn=jlmn-(jl+jm+1)+m1
               if (right_lmn>=left_lmn) then 
                 tot_lmn=right_lmn*(right_lmn-1)/2 + left_lmn
               else 
                 tot_lmn=left_lmn*(left_lmn-1)/2 + right_lmn
               end if
               sumr_ij=sumr_ij+vpawu(1,m1,m2,isppol)*rcart_onsite(:,tot_lmn,iatom)
             end if
           end do !m1
         end do !m2

         rij_tmp(:,ij_idx,isppol)=sumr_ij(:)

       end do !ilmn
     end do !jlmn
   end do !isppol

   ABI_FREE(vpawu)

   ! === Save values in packed form ===
   ABI_MALLOC(Hur(iatom)%ij_select,(lmn_size,lmn_size,nsppol))
   Hur(iatom)%ij_select=0
   nsel(:,:)=COUNT(ABS(rij_tmp)>tol6,DIM=2)
   nmax=MAXVAL(nsel)
   ABI_MALLOC(Hur(iatom)%commutator,(3,nmax,nsppol))
   do isppol=1,nsppol
     ij_idx=0
     isel  =0
     do jlmn=1,lmn_size
       do ilmn=1,lmn_size
         ij_idx=ij_idx+1
         if (ANY (ABS(rij_tmp(:,ij_idx,isppol))>tol6) ) then
           isel=isel+1
           Hur(iatom)%ij_select(ilmn,jlmn,isppol)=isel
           Hur(iatom)%commutator(:,isel,isppol)=rij_tmp(:,ij_idx,isppol)
         end if
       end do
     end do
   end do

   ABI_FREE(rij_tmp)
 end do !iatom

 ABI_FREE(rcart_onsite)

end subroutine pawhur_init
!!***

!----------------------------------------------------------------------

!!****f* m_paw_hr/pawr
!! NAME
!! pawr
!!
!! FUNCTION
!! Evaluate matrix elements of the position operator between PAW AE partial waves.
!!
!! INPUTS
!!  Pawtab(ntypat) <type(pawtab_type)>=paw tabulated data read at start:
!!     %lmn_size
!!     %lmn2_size
!!     %indklmn
!!     %phiphj
!!  Pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data:
!!     %mesh_size=Dimension of radial mesh
!!     %rad(mesh_size)=The coordinates of all the points of the radial mesh
!!  Pawang <type(pawang_type)>=paw angular mesh and related data
!!     %lmax=Maximum value of angular momentum l+1
!!     %gntselect((2*l_max-1)**2,l_max**2,l_max**2)= selection rules for Gaunt coefficients
!!     %realgnt
!!  natom=number of atoms in unit cell
!!  ntypat=number of types of atom
!!  typat(natom)=type of each atom
!!  xcart(3,natom)=cartesian coordinates
!!
!! OUTPUT
!!  rcart_onsite(3,lmn2_size_max,natom)
!!
!! PARENTS
!!      m_paw_hr
!!
!! CHILDREN
!!      simp_gen
!!
!! SOURCE

subroutine pawr(Pawtab,Pawrad,Pawang,natom,ntypat,typat,xcart,lmn2_size_max,rcart_onsite)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lmn2_size_max,natom,ntypat
 type(Pawang_type),intent(in) :: Pawang

!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: xcart(3,natom)
 real(dp),intent(inout) :: rcart_onsite(3,lmn2_size_max,natom)
 type(Pawrad_type),intent(in) :: Pawrad(ntypat)
 type(Pawtab_type),target,intent(in) :: Pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer,parameter :: ll1=1
 integer :: iatom,idir,ignt,il,ilm,ilm_G,ilmn,iln,im,itypat,jl,jlm,jlmn,jln,jm,k0lm
 integer :: k0lmn,k0ln,klm,klmn,kln,lmn_size,mesh_size,mm_G,lmn2_size
 real(dp) :: fact,intff,rgnt
!arrays
 integer,ABI_CONTIGUOUS pointer :: indlmn(:,:)
 real(dp),allocatable :: ff(:),rad(:),rc_tmp(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 fact=two*SQRT(pi/three)
 rcart_onsite(:,:,:)=zero

 do itypat=1,ntypat
   lmn_size  =Pawtab(itypat)%lmn_size
   lmn2_size =Pawtab(itypat)%lmn2_size
   mesh_size =Pawtab(itypat)%mesh_size
   indlmn =>  Pawtab(itypat)%indlmn

   ABI_MALLOC(ff,(mesh_size))
   ABI_MALLOC(rad,(mesh_size))
   rad(1:mesh_size)=Pawrad(itypat)%rad(1:mesh_size)

   ABI_MALLOC(rc_tmp,(3,lmn2_size))
   rc_tmp=zero
   !
   ! === Loop on (jl,jm,jn) channels
   do jlmn=1,lmn_size
     jl =indlmn(1,jlmn)
     jm =indlmn(2,jlmn)
     jlm=indlmn(4,jlmn)
     jln=indlmn(5,jlmn)

     k0lmn=jlmn*(jlmn-1)/2
     k0lm =jlm *(jlm -1)/2
     k0ln =jln *(jln -1)/2
     !
     ! === Loop on (il,im,in) channels; klmn is the index for packed form ===
     do ilmn=1,jlmn
       il =indlmn(1,ilmn)
       im =indlmn(2,ilmn)
       ilm=indlmn(4,ilmn)
       iln=indlmn(5,ilmn)

       klmn=k0lmn+ilmn
       klm =k0lm +ilm
       kln =k0ln +iln
       !
       ! === For each cartesian direction, use expansion in terms of RSH ===
       ! TODO Add a check if l=1 is in the set
       do idir=1,3
         mm_G=0
         if (idir==1) mm_G= 1
         if (idir==2) mm_G=-1
         if (idir==3) mm_G= 0
         ilm_G=1+ll1**2+ll1+mm_G
         ignt=Pawang%gntselect(ilm_G,klm)
         if (ignt/=0) then
           rgnt=Pawang%realgnt(ignt)
           ff(1)=zero
           !ff(2:mesh_size)=(Pawtab(itypat)%phiphj(2:mesh_size,kln)-Pawtab(itypat)%tphitphj(2:mesh_size,kln))*rad(2:mesh_size)
           ff(2:mesh_size)=Pawtab(itypat)%phiphj(2:mesh_size,kln)*rad(2:mesh_size)
           call simp_gen(intff,ff,Pawrad(itypat))
           rc_tmp(idir,klmn)=fact*intff*rgnt
         end if
       end do !idir

     end do !ilmn
   end do !jllmn

   ! === Make matrix elements for each atom of this type ===
   do jlmn=1,lmn_size
     jl =indlmn(1,jlmn)
     jm =indlmn(2,jlmn)
     jln=indlmn(5,jlmn)

     k0lmn=jlmn*(jlmn-1)/2
     k0ln =jln *(jln -1)/2
     do ilmn=1,jlmn
       il =indlmn(1,ilmn)
       im =indlmn(2,ilmn)
       iln=indlmn(5,ilmn)

       klmn=k0lmn+ilmn
       kln =k0ln +iln

       intff=zero
       if (il==jl.and.jm==im) then
         ff(1:mesh_size)=Pawtab(itypat)%phiphj(1:mesh_size,kln)
         call simp_gen(intff,ff,Pawrad(itypat))
       end if
       do iatom=1,natom
         if (typat(iatom)/=itypat) CYCLE
         rcart_onsite(:,klmn,iatom)=rc_tmp(:,klmn) + xcart(:,iatom)*intff
       end do

     end do ! ilmn
   end do !jlmn

   ABI_FREE(ff)
   ABI_FREE(rad)
   ABI_FREE(rc_tmp)
 end do !itypat

 DBG_EXIT("COLL")

end subroutine pawr
!!***

!----------------------------------------------------------------------

END MODULE m_paw_hr
!!***
