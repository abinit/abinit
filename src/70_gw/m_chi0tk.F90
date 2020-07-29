!!****m* ABINIT/m_chi0tk
!! NAME
!!  m_chi0tk
!!
!! FUNCTION
!!  This module provides tools for the computation of the irreducible polarizability.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2020 ABINIT group (MG, FB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_chi0tk

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_xomp
 use m_sort
 use m_wfd

 use defs_datatypes, only : ebands_t
 use m_gwdefs,   only : GW_TOL_DOCC, czero_gw, cone_gw, one_gw, em1params_t, j_gw
 use m_fstrings, only : sjoin, itoa
 use m_hide_blas,only : xgerc, xgemm, xherk, xher
 use m_crystal,  only : crystal_t
 use m_gsphere,  only : gsphere_t, gsph_gmg_idx, gsph_gmg_fftidx
 use m_bz_mesh,  only : littlegroup_t, kmesh_t, has_BZ_item

 implicit none

 private

 public :: assemblychi0_sym
 public :: symmetrize_afm_chi0
 public :: accumulate_chi0_q0
 public :: accumulate_sfchi0_q0
 public :: assemblychi0sf
 public :: approxdelta
 public :: setup_spectral
 public :: hilbert_transform
 public :: hilbert_transform_headwings
 public :: completechi0_deltapart
 public :: output_chi0sumrule
 public :: accumulate_chi0sumrule
 public :: make_transitions
 public :: chi0_bbp_mask
!!***

CONTAINS  !=======================================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_chi0tk/assemblychi0_sym
!! NAME
!! assemblychi0_sym
!!
!! FUNCTION
!! Update the independent particle susceptibility for the contribution
!! of one pair of occupied-unoccupied band, for each frequency.
!! If symchi=1 the expression is symmetrized taking into account the symmetries
!! of the little group associated to the external q-point.
!! Compute chi0(G1,G2,io)=chi0(G1,G2,io)+\sum_S \hat S (rhotwg(G1)*rhotwg*(G2))*green_w(io)
!! where S are the symmetries of the little group associated to the external q-point.
!!
!! INPUTS
!!  nspinor=Number of spinorial components.
!!  ik_bz=Index of the k-point in the BZ array whose contribution has to be symmetrized and added to cchi0
!!  npwepG0=Maximum number of G vectors taking into account possible umklapp G0, ie enlarged sphere G-G0
!!  rhotwg(npwe)=Oscillator matrix elements for this k-point and the transition that has to be summed
!!  green_w(nomega)=frequency dependent part coming from the green function
!!  Gsph_epsG0<gsphere_t> Information on the "enlarged" G-sphere used for chi0, it contains umklapp G0 vectors
!!    %ng=number of G vectors in the enlarged sphere, actually MUST be equal to the size of rhotwg
!!    %rottbm1(ng,2,nsym)=index of (IR)^{-1} G where I is the identity or the inversion
!!    %phmGt(ng,nsym)=phase factors associated to non-simmorphic operations
!!  Ltg_q<littlegroup_t_type>=Info on the little group associated to the external q-point.
!!    %timrev=2 it time-reversal is used, 1 otherwise
!!    %nsym_sg=Number of space group symmetries
!!    %wtksym(2,nsym,nkbz)=1 if the symmetry (with or without time-reversal) must be considered for this k-point
!!    %flag_umklp(timrev,nsym)= flag for umklapp processes
!!      if 1 that the particular operation (IS) requires a G_o to preserve Q, 0 otherwise
!!    %igmG0(npwepG0,timrev,nsym) index of G-G0 in the array gvec
!!  Ep<em1params_t>=Parameters related to the calculation of chi0/epsilon^-1
!!    %symchi
!!    %nomega=number of frequencies
!!    %npwe=number of plane waves for epsilon (input variable)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  chi0(npwe,npwe,nomega)=independent-particle susceptibility matrix in reciprocal space
!!
!! PARENTS
!!      cchi0,cchi0q0_intraband
!!
!! CHILDREN
!!
!! SOURCE


subroutine assemblychi0_sym(ik_bz,nspinor,Ep,Ltg_q,green_w,npwepG0,rhotwg,Gsph_epsG0,chi0)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz,npwepG0,nspinor
 type(gsphere_t),intent(in) :: Gsph_epsG0
 type(littlegroup_t),intent(in) :: Ltg_q
 type(em1params_t),intent(in) :: Ep
!arrays
 complex(gwpc),intent(in) :: rhotwg(npwepG0)
 complex(dpc),intent(in) :: green_w(Ep%nomega)
 complex(gwpc),intent(inout) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)

!Local variables-------------------------------
!scalars
 integer :: itim,io,isym,ig1,ig2,nthreads
 integer :: isymop,nsymop
 real(gwp) :: dr
 complex(gwpc) :: dd
 !character(len=500) :: msg
!arrays
 integer :: Sm1_gmG0(Ep%npwe)
 complex(gwpc),allocatable :: rhotwg_sym(:,:)

! *************************************************************************

 ABI_UNUSED(nspinor)

 nthreads = xomp_get_max_threads()

 SELECT CASE (Ep%symchi)
 CASE (0)
   ! Do not use symmetries

   ! note that single precision is faster (sometimes factor ~2).
   ! Rely on MKL threads for OPENMP parallelization

   do io=1,Ep%nomega
     ! Check if green_w(io) is real (=> pure imaginary omega)
     ! if yes, the corresponding chi0(io) is hermitian
     if( ABS(AIMAG(green_w(io))) < 1.0e-6_dp ) then
       dr=green_w(io)
       call xher('U',Ep%npwe,dr,rhotwg,1,chi0(:,:,io),Ep%npwe)
     else
       dd=green_w(io)
       call xgerc(Ep%npwe,Ep%npwe,dd,rhotwg,1,rhotwg,1,chi0(:,:,io),Ep%npwe)
     endif
   end do

 CASE (1)
   ! Use symmetries to reconstruct the integrand in the BZ.
   !
   ! Notes on the symmetrization of the oscillator matrix elements
   !  If  Sq = q then  M_G^( Sk,q)= e^{-i(q+G).t} M_{ S^-1G}  (k,q)
   !  If -Sq = q then  M_G^(-Sk,q)= e^{-i(q+G).t} M_{-S^-1G}^*(k,q)
   !
   ! In case of an umklapp process
   !  If  Sq = q+G0 then  M_G( Sk,q)= e^{-i(q+G).t} M_{ S^-1(G-G0}   (k,q)
   !  If -Sq = q+G0 then  M_G(-Sk,q)= e^{-i(q+G).t} M_{-S^-1(G-G0)}^*(k,q)
   !
   ! Ltg_q%igmG0(ig,itim,isym) contains the index of G-G0 where ISq=q+G0
   ! Note that there is no need to take into account the phases due to q,
   ! They cancel in the scalar product ==> phmGt(G,isym)=e^{-iG\cdot t}
   !
   ! Mind the slicing of %rottbm1(npwepG0,timrev,nsym) and %phmGt(npwepG0,nsym) as
   ! these arrays, usually, do not conform to rho_twg_sym(npw) !
   !
   ! Loop over symmetries of the space group and time-reversal.
   nsymop = count(Ltg_q%wtksym(:,:,ik_bz)==1)
   ABI_MALLOC(rhotwg_sym,(Ep%npwe,nsymop))
   isymop = 0

   ! Prepare all the rhotwg at once to use BLAS level 3 routines
   do isym=1,Ltg_q%nsym_sg
     do itim=1,Ltg_q%timrev
       if (Ltg_q%wtksym(itim,isym,ik_bz)==1) then
         ! This operation belongs to the little group and has to be used to reconstruct the BZ ===
         ! * In the following 3 lines mind the slicing (1:npwe)
         ! TODO this is a hot-spot, should add a test on the umklapp
         !
         !gmG0 => Ltg_q%igmG0(1:Ep%npwe,itim,isym)
         Sm1_gmG0(1:Ep%npwe) = Gsph_epsG0%rottbm1( Ltg_q%igmG0(1:Ep%npwe,itim,isym), itim,isym)

         isymop = isymop + 1
         SELECT CASE (itim)
         CASE (1)
           rhotwg_sym(1:Ep%npwe,isymop) = rhotwg(Sm1_gmG0) * Gsph_epsG0%phmGt(1:Ep%npwe,isym)
         CASE (2)
           rhotwg_sym(1:Ep%npwe,isymop) = GWPC_CONJG(rhotwg(Sm1_gmG0))*Gsph_epsG0%phmGt(1:Ep%npwe,isym)
         CASE DEFAULT
           MSG_BUG(sjoin('Wrong itim:', itoa(itim)))
         END SELECT
       end if
     end do
   end do

   ! Multiply rhotwg_sym by green_w(io) and accumulate in chi0(G,Gp,io)
   ! note that single precision is faster (sometimes factor ~2).
   ! Rely on MKL threads for OPENMP parallelization
   do io=1,Ep%nomega
     ! Check if green_w(io) is real (=> pure imaginary omega)
     ! if yes, the corresponding chi0(io) is hermitian
     if( ABS(AIMAG(green_w(io))) < 1.0e-6_dp ) then
       dr=green_w(io)
       call xherk('U','N',Ep%npwe,nsymop,dr,rhotwg_sym,Ep%npwe,one_gw,chi0(:,:,io),Ep%npwe)
     else
       dd=green_w(io)
       call xgemm('N','C',Ep%npwe,Ep%npwe,nsymop,dd,rhotwg_sym,Ep%npwe,rhotwg_sym,Ep%npwe,cone_gw,chi0(:,:,io),Ep%npwe)
     endif
   end do

   ABI_FREE(rhotwg_sym)

 CASE DEFAULT
   MSG_BUG(sjoin('Wrong symchi:', itoa(Ep%symchi)))
 END SELECT

end subroutine assemblychi0_sym
!!***

!----------------------------------------------------------------------

!!****f* m_chi0tk/mkrhotwg_sigma
!! NAME
!! mkrhotwg_sigma
!!
!! FUNCTION
!!  Helper function used to calculate selected linear combination
!!  of the oscillator matrix elements in the case of non-collinear magnetism.
!!
!! INPUTS
!!  ii=Index selecting the particolar combination of spin components.
!!  npw=Number of plane-waves in the oscillators.
!!  nspinor=Number of spinorial components.
!!  rhotwg(npw*nspinor**2)=OScillator matrix elements.
!!
!! OUTPUT
!!  rhotwg_I(npw)=Required linear combination of the oscillator matrix elements.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkrhotwg_sigma(ii,nspinor,npw,rhotwg,rhotwg_I)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ii,npw,nspinor
!arrays
 complex(gwpc),intent(in) :: rhotwg(npw*nspinor**2)
 complex(gwpc),intent(out) :: rhotwg_I(npw)

! *************************************************************************

 SELECT CASE (ii)
 CASE (1)
   ! $ M_0 = M_{\up,\up} + M_{\down,\down} $
   rhotwg_I(:) = rhotwg(1:npw) + rhotwg(npw+1:2*npw)
 CASE (2)
   ! $ M_z = M_{\up,\up} - M_{\down,\down} $
   rhotwg_I(:) = rhotwg(1:npw) - rhotwg(npw+1:2*npw)
 CASE (3)
   ! $ M_x = M_{\up,\down} + M_{\down,\up} $
   rhotwg_I(:) = ( rhotwg(2*npw+1:3*npw) + rhotwg(3*npw+1:4*npw) )
 CASE (4)
   ! $ M_y = i * (M_{\up,\down} -M_{\down,\up}) $
   rhotwg_I(:) = (rhotwg(2*npw+1:3*npw) - rhotwg(3*npw+1:4*npw) )*j_gw
 CASE DEFAULT
   MSG_BUG(sjoin('Wrong ii value:', itoa(ii)))
 END SELECT

end subroutine mkrhotwg_sigma
!!***

!----------------------------------------------------------------------

!!****f* m_chi0tk/symmetrize_afm_chi0tk
!! NAME
!! symmetrize_afm_chi0
!!
!! FUNCTION
!!  Reconstruct the (down, down) component of the irreducible polarizability
!!  starting from the (up,up) element in case of systems with AFM symmetries
!!  (i.e nspden==2 and nsppol=1). Return the trace (up,up)+(down,down) of the
!!  matrix as required by GW calculations.
!!
!! INPUTS
!!  Cryst<crystal_t>= Information on symmetries and unit cell.
!!  Gsph<gsphere_t>= The G-sphere used to descrive chi0.
!!  npwe=Number of G-vectors in chi0.
!!  nomega=number of frequencies.
!!  Ltg_q<littlegroup_t>=Structure with useful table describing the little group of the q-point.
!!
!! SIDE EFFECTS
!! chi0(npwe,npwe,nomega)= In input the up-up component, in output the trace of chi0.
!!   The value of matrix elements that should be zero due to AFM symmetry properties are
!!   forced to be zero (see NOTES below).
!! [chi0_lwing(npwe,nomega,3)] = Lower wings, symmetrized in output.
!! [chi0_uwing(npwe,nomega,3)] = Upper wings, symmetrized in output.
!! [chi0_head(3,3,nomega)    ] = Head of chi0, symmetrized  in output.
!!
!! NOTES
!! In the case of magnetic group Shubnikov type III:
!!   For each set of paired FM-AFM symmetries, the down-down component of
!!   a generic response function in reciprocal space can be obtained according to:
!!
!!    chi^{down,down}_{G1,G2}(q) = chi^{up up}_{G1,G2}(q) e^{iS(G1-G2).(tnonsFM - tnonsAFM)}
!!
!!   where S is the rotational part common to the FM-AFM pair, tnonsFM and tnonsAFM
!!   are the fractional translations associated to the ferromagnetic and antiferromagnetic symmetry, respectively.
!!   Note that, if for a given G1-G2 pair, the phase e^{iS(G1-G2).(tnonsFM - tnonsAFM) depends
!!   on the FM-AFM symmetry pair, then the corresponding matrix element of chi0 must be zero.
!!   Actually this is manually enforced in the code because this property might not be
!!   perfectly satisfied due to round-off errors.
!!
!! In the case of magnetic group Shubnikov type III:
!!   Only the AFM symmetries that preserve the external q-point (with or without time-reversal)
!!   are used to get the (down, down) component using the fact that:
!!
!!    chi^{down,down}_{G1,G2}(Sq) = chi^{up up}_{S^{-1}G1,S^{-1}G2}(q) e^{i(G2-G1).tnons_S }
!!
!!   Actually we perform an average over subset of the little group of q with AFM character in
!!   order to reduce as much as possible errors due to round off errors. In brief we evaluate:
!!
!!    1/N_{Ltq} \sum_{S\in Ltg AFM} chi^{up up}_{S^{-1}G1,S^{-1}G2}(q) e^{i(G2-G1).tnons_S }
!!
!!   where N_{Ltg} is the number of AFM operation in the little group (time reversal included)
!!
!! TODO
!!  It is possible to symmetrize chi0 without any the extra allocation for afm_mat.
!!  More CPU demanding but safer in case of a large chi0 matrix. One might loop over G1 and G2 shells ...
!!
!! PARENTS
!!      cchi0,cchi0q0,cchi0q0_intraband
!!
!! CHILDREN
!!
!! SOURCE

subroutine symmetrize_afm_chi0(Cryst,Gsph,Ltg_q,npwe,nomega,chi0,chi0_head,chi0_lwing,chi0_uwing)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwe,nomega
 type(gsphere_t),intent(in) :: Gsph
 type(crystal_t),intent(in) :: Cryst
 type(littlegroup_t),intent(in) :: Ltg_q
!arrays
 complex(gwpc),intent(inout) :: chi0(npwe,npwe,nomega)
 complex(dpc),optional,intent(inout) :: chi0_lwing(npwe,nomega,3)
 complex(dpc),optional,intent(inout) :: chi0_uwing(npwe,nomega,3)
 complex(dpc),optional,intent(inout) :: chi0_head(3,3,nomega)

!Local variables ------------------------------
!scalars
 integer :: io,ig1,ig2,isymf,isyma,isym,ipair,k0g,kg,npairs,nonzero
 integer :: iSmg1,iSmg2,itim,shubnikov,ntest
 complex(gwpc) :: phase,phase_old,sumchi,ctmp
 logical :: found
 !character(len=500) :: msg
!arrays
 integer :: rotfm(3,3),rotafm(3,3),pairs2sym(2,Cryst%nsym/2)
 real(dp) :: tfm(3),tafm(3)
 complex(gwpc),allocatable :: afm_mat(:),chi0_afm(:,:)

!************************************************************************

 ABI_CHECK(ANY(Cryst%symafm==-1),'Not magnetic space group')
 !
 ! ==== Find shubnikov type ====
 ! TODO This info should be stored in Cryst%
 shubnikov=4; npairs=0

 do isymf=1,Cryst%nsym
   if (Cryst%symafm(isymf)==-1) CYCLE
   rotfm = Cryst%symrec(:,:,isymf)
   tfm = Cryst%tnons(:,isymf)
   found = .FALSE.

   do isyma=1,Cryst%nsym
     if (Cryst%symafm(isyma)==1) CYCLE
     rotafm = Cryst%symrec(:,:,isyma)

     if (ALL(rotfm==rotafm)) then
       found=.TRUE.
       tafm = Cryst%tnons(:,isyma)
       npairs=npairs+1
       ABI_CHECK(npairs<=Cryst%nsym/2,'Wrong AFM group')
       pairs2sym(1,npairs)=isymf
       pairs2sym(2,npairs)=isyma
     end if
   end do !isyma

   if (.not.found) then
    shubnikov=3; EXIT !isymf
   end if
 end do !isymf

 select case (shubnikov)

 case (4)
   call wrtout(std_out,' Found Magnetic group Shubnikov type IV','COLL')
   ABI_CHECK(npairs==Cryst%nsym/2,'Wrong AFM space group')

   ABI_MALLOC(afm_mat,(npwe*(npwe+1)/2))

! jmb
   phase_old=zero

   do ig2=1,npwe
     k0g=ig2*(ig2-1)/2
     do ig1=1,ig2
       kg=k0g+ig1
       nonzero=1

       do ipair=1,Cryst%nsym/2
         isymf = pairs2sym(1,ipair)
         isyma = pairs2sym(2,ipair)
         phase = ( Gsph%phmSGt(ig1,isymf)*CONJG(Gsph%phmSGt(ig1,isyma)) ) * &
&                ( Gsph%phmSGt(ig2,isymf)*CONJG(Gsph%phmSGt(ig2,isyma)) )
         if (ipair>1 .and. (ABS(phase_old-phase) > tol6)) then
           nonzero=0; EXIT
         end if
         phase_old=phase
       end do !ipair

       afm_mat(kg)=nonzero*(cone_gw + phase)
     end do !ig1
   end do !ig2
   !
   ! =======================================================================
   ! ==== Symmetrize chi0 constructing chi0^{\up\up} + chi0{\down\down} ====
   ! =======================================================================
   !
   !  head^{\down\down} = head^{\up\up}
   if (PRESENT(chi0_head)) chi0_head = two * chi0_head
   !
   ! w^{\down\down}_{0 G'} =  e^{-iSG'.(tFM-tAFM)} w^{\up\up}_{0 G'}.
   ! w^{\down\down}_{G 0 } =  e^{+iSG .(tFM-tAFM)} w^{\up\up}_{G 0 }.
   if (PRESENT(chi0_uwing)) then
     do io=1,nomega
       do ig2=1,npwe
         k0g=ig2*(ig2-1)/2
         kg=k0g+1
         chi0_uwing(ig2,io,:)=afm_mat(kg)*chi0_uwing(ig2,io,:)
       end do
     end do
   end if

   if (PRESENT(chi0_lwing)) then
     do io=1,nomega
       do ig1=1,npwe
         k0g=ig1*(ig1-1)/2
         kg=k0g+1
         chi0_lwing(ig1,io,:)=CONJG(afm_mat(kg))*chi0_lwing(ig1,io,:)
       end do
     end do
   end if

   do io=1,nomega
     ! Take care of diagonal.
     do ig1=1,npwe
       chi0(ig1,ig1,io)=two*chi0(ig1,ig1,io)
     end do

     ! Upper and lower triangle are treated differently:
     ! We took advantage of the fact the afm_mat is hermitian to reduce memory.
     do ig2=2,npwe
       k0g=ig2*(ig2-1)/2
       do ig1=1,ig2-1
         kg=k0g+ig1
         chi0(ig1,ig2,io)=afm_mat(kg)*chi0(ig1,ig2,io)
       end do
     end do

     do ig1=2,npwe
       k0g=ig1*(ig1-1)/2
       do ig2=1,ig1-1
         kg=k0g+ig2
         chi0(ig1,ig2,io)=CONJG(afm_mat(kg))*chi0(ig1,ig2,io)
       end do
     end do

   end do !io

   ABI_FREE(afm_mat)

 case (3)
   call wrtout(std_out,' Found Magnetic group Shubnikov type III',"COLL")
   MSG_ERROR('Shubnikov type III not implemented')

   ntest=0
   do itim=1,ltg_q%timrev
     do isym=1,ltg_q%nsym_sg
       ! use only afm sym preserving q with and without time-reversal
       if ( cryst%symafm(isym)==-1 .and. ltg_q%preserve(itim,isym)==1 ) ntest=ntest+1
     end do
   end do

   if (ntest==0) then
       MSG_WARNING("no symmetry can be used!")
   end if
   !RETURN
   ABI_MALLOC(chi0_afm,(npwe,npwe))

   do io=1,nomega

     do ig2=1,npwe
       do ig1=1,npwe
         sumchi=czero_gw

         do itim=1,ltg_q%timrev
           do isym=1,ltg_q%nsym_sg
             ! use only afm sym preserving q with and without time-reversal
             if ( cryst%symafm(isym)==-1 .and. ltg_q%preserve(itim,isym)==1 ) then
               phase =  Gsph%phmGt(ig1,isym)*CONJG(Gsph%phmGt(ig2,isym))
               iSmg1=Gsph%rottbm1(ig1,itim,isym)
               iSmg2=Gsph%rottbm1(ig2,itim,isym)
               ctmp=chi0(iSmg1,iSmg2,io)*phase !; if (itim==2) ctmp=CONJG(ctmp) !check this
               sumchi=sumchi+ctmp !chi0(iSmg1,iSmg2,io)*phase
             end if
           end do ! isym
         end do !itim

         chi0_afm(ig1,ig2)=sumchi/Ltg_q%nsym_ltg  !has to be changed in case of time-reversal
       end do !ig1
     end do !ig2

     ! We want chi_{up,up} +chi_{dwn,dwn}.
     chi0(:,:,io)=chi0(:,:,io)+chi0_afm(:,:)
   end do !iomega

   ABI_FREE(chi0_afm)

 case default
   MSG_BUG(sjoin('Wrong value for shubnikov= ', itoa(shubnikov)))
 end select

end subroutine symmetrize_afm_chi0
!!***

!----------------------------------------------------------------------

!!****f* m_chi0tk/accumulate_chi0_q0
!! NAME
!! accumulate_chi0_q0
!!
!! FUNCTION
!! Update the independent particle susceptibility at q==0 for the contribution
!! of one pair of occupied-unoccupied band, for each frequency.
!! This routine takes advantage of the symmetries of the little group of the external q-point
!! to symmetrize the contribution arising from the input k-point located in the IBZ_q.
!! It computes:
!!
!!   $ \chi_0(G1,G2,io) = \chi_0(G1,G2,io)+\sum_S (rhotwg(G1)*rhotwg^\dagger(G2))*green_w(io) $
!!
!! where S is a symmetry in reciprocal space.
!! The matrix elements of the gradient operator and [V_{nl},r] are symmetrized as well.
!!
!! INPUTS
!!  ik_bz=Index of the k-point whose contribution has to be added to chi0.
!!  isym_kbz=Index of the symmetry such that k_bz = IS k_ibz
!!  itim_kbz=2 if time-reversal has to be used to obtain k_bz, 1 otherwise.
!!  npwepG0=Maximum number of G vectors
!!  rhotwg(npwepG0)=Oscillator matrix elements corresponding to an occupied-unoccupied pair of states.
!!  rhotwx(3,nspinor**2)=Matrix element of the operator $-i[H,r]/(e1-e2) = -i r$ in reciprocal lattice units.
!!  green_w(nomega)=Frequency dependent part of the Green function.
!!  Ltg_q<littlegroup_t_type>=Info on the little group associated to the external q-point.
!!    %timrev=2 it time-reversal is used, 1 otherwise.
!!    %nsym_sg=Number of space group symmetries.
!!    %wtksym(2,nsym,nkbz)=1 if the symmetry (with or without time-reversal) must be considered for this k-point.
!!  Gsph_epsG0<gsphere_t> Information on the "enlarged" G-sphere used for chi0, it contains umklapp G0 vectors
!!    %ng=number of G vectors in the enlarged sphere, actually MUST be equal to the size of rhotwg.
!!    %rottbm1(ng,2,nsym)=index of (IR)^{-1} G where I is the identity or the inversion.
!!    %phmGt(ng,nsym)=phase factors associated to non-simmorphic operations.
!!  Cryst<crystal_t>=Structure defining the unit cell and its symmetries
!!    %nsym=Number of symmetries.
!!    %symrec(3,3,nsym)=Symmetry operation in reciprocal space (reduced coordinates)
!!  Ep<em1params_t>=Parameters of the chi0 calculation.
!!     %npwe=number of plane waves in chi0.
!!     %symchi=1 if symmetrization has to be performed.
!!     %nomega=number of frequencies in chi0.
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  chi0(npwe,npwe,nomega)= Updated independent-particle susceptibility matrix in reciprocal space at q==0.
!!  chi0_head(3,3,Ep%nomega)=Head.
!!  chi0_lwing(Ep%npwe*Ep%nI,Ep%nomega,3)=Lower wing.
!!  chi0_uwing(Ep%npwe*Ep%nJ,Ep%nomega,3)=Upper wing.
!!
!! NOTES
!!
!!  1) Symmetrization of the oscilator matrix elements.
!!    If  Sq = q then  M_G( Sk,q)= e^{-i(q+G).\tau} M_{ S^-1G}  (k,q)
!!    If -Sq = q then  M_G(-Sk,q)= e^{-i(q+G).\tau} M_{-S^-1G}^*(k,q)
!!
!!    In the case of umklapps:
!!    If  Sq = q+G0 then  M_G( Sk,q)= e^{-i(q+G).\tau} M_{ S^-1(G-G0}   (k,q)
!!    If -Sq = q+G0 then  M_G(-Sk,q)= e^{-i(q+G).\tau} M_{-S^-1(G-G0)}^*(k,q)
!!
!!  In the equation below there is no need to take into account the phases due to q.t
!!  as they cancel each other in the scalar product ==> only phmGt(G,isym)=e^{-iG.\tau} is needed.
!!
!!  2) Symmetrization of the matrix elements of the position operator.
!!
!!    <Sk,b|\vec r| Sk,b'> = <k b| R\vec r + \tau|k b'>
!!
!!     where S is one of the symrec operation, R and \tau is the corresponding
!!     operation in real space. The term involving the fractional translation is zero provided that b /= b'.
!!
!! PARENTS
!!      cchi0q0
!!
!! CHILDREN
!!
!! SOURCE

subroutine accumulate_chi0_q0(ik_bz,isym_kbz,itim_kbz,gwcomp,nspinor,npwepG0,Ep,Cryst,Ltg_q,Gsph_epsG0,&
& chi0,rhotwx,rhotwg,green_w,green_enhigh_w,deltaf_b1b2,chi0_head,chi0_lwing,chi0_uwing)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz,isym_kbz,itim_kbz,npwepG0,nspinor,gwcomp
 real(dp),intent(in) :: deltaf_b1b2
 type(littlegroup_t),intent(in) :: Ltg_q
 type(gsphere_t),target,intent(in) :: Gsph_epsG0
 type(crystal_t),intent(in) :: Cryst
 type(em1params_t),intent(in) :: Ep
!arrays
 complex(gwpc),intent(in) :: rhotwg(npwepG0)
 complex(gwpc),intent(in) :: rhotwx(3,nspinor**2)
 complex(gwpc),intent(inout) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)
 complex(dpc),intent(in) :: green_w(Ep%nomega),green_enhigh_w(Ep%nomega)
 complex(dpc),intent(inout) :: chi0_head(3,3,Ep%nomega)
 complex(dpc),intent(inout) :: chi0_lwing(Ep%npwe*Ep%nI,Ep%nomega,3)
 complex(dpc),intent(inout) :: chi0_uwing(Ep%npwe*Ep%nJ,Ep%nomega,3)

!Local variables-------------------------------
!scalars
 integer :: itim,io,isym,idir,jdir
 complex(gwpc) :: dd
 !character(len=500) :: msg
!arrays
 integer,ABI_CONTIGUOUS pointer :: Sm1G(:)
 complex(dpc) :: mir_kbz(3)
 complex(gwpc),allocatable :: rhotwg_sym(:)
 complex(gwpc), ABI_CONTIGUOUS pointer :: phmGt(:)

!************************************************************************

 ABI_UNUSED(deltaf_b1b2)

 SELECT CASE (Ep%symchi)
 CASE (0)
   ! Do not use symmetries.
   ! Symmetrize rhotwg in the full BZ and accumulate over the full BZ i.e.
   !   chi0(G1,G2,io) = chi0(G1,G2,io) + (rhotwg(G1)*CONJG(rhotwg(G2)))*green_w(io)
   !
   ! The non-analytic term is symmetrized for this k-point in the BZ according to:
   !    rhotwg(1)= S^-1q * rhotwx_ibz
   !    rhotwg(1)=-S^-1q * CONJG(rhotwx_ibz) if time-reversal is used.

   ! Multiply elements G1,G2 of rhotwg_sym by green_w(io) and accumulate in chi0(G1,G2,io)
!$omp parallel do private(dd)
   do io=1,Ep%nomega
     dd=green_w(io)
     call XGERC(Ep%npwe,Ep%npwe,dd,rhotwg,1,rhotwg,1,chi0(:,:,io),Ep%npwe)
   end do

   ! === Accumulate heads and wings for each small q ===
   ! FIXME extrapolar method should be checked!!
   ! Symmetrize <r> in full BZ: <Sk b|r|Sk b'> = R <k b|r|k b'> + \tau \delta_{bb'}
   if (nspinor == 1) then
     mir_kbz = (3-2*itim_kbz) * MATMUL(Cryst%symrec(:,:,isym_kbz),rhotwx(:,1))
   else
     mir_kbz = (3-2*itim_kbz) * MATMUL(Cryst%symrec(:,:,isym_kbz), sum(rhotwx(:,1:2), dim=2))
   end if
   if (itim_kbz==2) mir_kbz=CONJG(mir_kbz)

   ! here we might take advantage of Hermiticity along Im axis in RPA (see mkG0w)
   do idir=1,3
     do io=1,Ep%nomega
       chi0_uwing(:,io,idir) = chi0_uwing(:,io,idir) + green_w(io) * mir_kbz(idir)*CONJG(rhotwg)
       chi0_lwing(:,io,idir) = chi0_lwing(:,io,idir) + green_w(io) * rhotwg*CONJG(mir_kbz(idir))
       ! Add contribution due to extrapolar technique.
       !if (gwcomp==1.and.ABS(deltaf_b1b2) >= GW_TOL_DOCC) then
       if (gwcomp==1) then
         chi0_uwing(:,io,idir) = chi0_uwing(:,io,idir) + green_enhigh_w(io) * mir_kbz(idir)*CONJG(rhotwg)
         chi0_lwing(:,io,idir) = chi0_lwing(:,io,idir) + green_enhigh_w(io) * rhotwg*CONJG(mir_kbz(idir))
       end if
     end do
   end do

   ! Accumulate the head.
   do io=1,Ep%nomega
     do jdir=1,3
       do idir=1,3
         chi0_head(idir,jdir,io) = chi0_head(idir,jdir,io) + green_w(io) * mir_kbz(idir)*CONJG(mir_kbz(jdir))
         ! Add contribution due to extrapolar technique.
         !if (gwcomp==1.and.ABS(deltaf_b1b2) >= GW_TOL_DOCC) then
         if (gwcomp==1) then
           chi0_head(idir,jdir,io) = chi0_head(idir,jdir,io) + green_enhigh_w(io) * mir_kbz(idir)*CONJG(mir_kbz(jdir))
         end if
       end do
     end do
   end do

 CASE (1)
   ! Use symmetries to reconstruct the integrand.
   ABI_MALLOC(rhotwg_sym, (Ep%npwe))

   ! Loop over symmetries of the space group and time-reversal.
   do isym=1,Ltg_q%nsym_sg
     do itim=1,Ltg_q%timrev

       if (Ltg_q%wtksym(itim,isym,ik_bz)==1) then
         ! This operation belongs to the little group and has to be considered to reconstruct the BZ.
         phmGt => Gsph_epsG0%phmGt  (1:Ep%npwe,isym) ! In the 2 lines below note the slicing (1:npwe)
         Sm1G  => Gsph_epsG0%rottbm1(1:Ep%npwe,itim,isym)

         SELECT CASE (itim)
         CASE (1)
           rhotwg_sym(1:Ep%npwe)=rhotwg(Sm1G(1:Ep%npwe))*phmGt(1:Ep%npwe)
         CASE (2)
           rhotwg_sym(1:Ep%npwe)=CONJG(rhotwg(Sm1G(1:Ep%npwe)))*phmGt(1:Ep%npwe)
         CASE DEFAULT
           MSG_BUG(sjoin('Wrong value of itim:', itoa(itim)))
         END SELECT

         ! Multiply elements G1,G2 of rhotwg_sym by green_w(io) and accumulate in chi0(G,Gp,io)
!$omp parallel do private(dd)
         do io=1,Ep%nomega
           dd=green_w(io)
           call XGERC(Ep%npwe,Ep%npwe,dd,rhotwg_sym,1,rhotwg_sym,1,chi0(:,:,io),Ep%npwe)
         end do

         ! === Accumulate heads and wings for each small q ===
         ! FIXME extrapolar method should be checked!!

         ! Symmetrize <r> in full BZ: <Sk b|r|Sk b'> = R <k b|r|k b'> + \tau \delta_{bb'}
         if (nspinor == 1) then
           mir_kbz =(3-2*itim) * MATMUL(Cryst%symrec(:,:,isym),rhotwx(:,1))
         else
           mir_kbz = (3-2*itim) * MATMUL(Cryst%symrec(:,:,isym), sum(rhotwx(:,1:2), dim=2))
         end if

         if (itim==2) mir_kbz=CONJG(mir_kbz)

         ! here we might take advantage of Hermiticity along Im axis in RPA (see mkG0w)
         do idir=1,3
           do io=1,Ep%nomega
             chi0_uwing(:,io,idir) = chi0_uwing(:,io,idir) + green_w(io) * mir_kbz(idir)*CONJG(rhotwg_sym)
             chi0_lwing(:,io,idir) = chi0_lwing(:,io,idir) + green_w(io) * rhotwg_sym*CONJG(mir_kbz(idir))
             ! Add contribution due to extrapolar technique.
             !if (gwcomp==1.and.ABS(deltaf_b1b2)>=GW_TOL_DOCC) then
             if (gwcomp==1) then
               chi0_uwing(:,io,idir) = chi0_uwing(:,io,idir) + green_enhigh_w(io) * mir_kbz(idir)*CONJG(rhotwg_sym)
               chi0_lwing(:,io,idir) = chi0_lwing(:,io,idir) + green_enhigh_w(io) * rhotwg_sym*CONJG(mir_kbz(idir))
             end if
           end do
         end do

         ! Accumulate the head.
         do io=1,Ep%nomega
           do jdir=1,3
             do idir=1,3
                chi0_head(idir,jdir,io) = chi0_head(idir,jdir,io) +  green_w(io) * mir_kbz(idir)*CONJG(mir_kbz(jdir))
                ! Add contribution due to extrapolar technique.
                !if (gwcomp==1.and.ABS(deltaf_b1b2) >= GW_TOL_DOCC) then
                if (gwcomp==1) then
                  chi0_head(idir,jdir,io) = chi0_head(idir,jdir,io) + green_enhigh_w(io)*mir_kbz(idir)*CONJG(mir_kbz(jdir))
                end if
             end do
           end do
         end do

       end if !wtksym
     end do !itim
   end do !isym

   ABI_FREE(rhotwg_sym)

 CASE DEFAULT
   MSG_BUG(sjoin('Wrong value of symchi ',itoa(Ep%symchi)))
 END SELECT

end subroutine accumulate_chi0_q0
!!***

!----------------------------------------------------------------------

!!****f* m_chi0tk/accumulate_sfchi0_q0
!! NAME
!! accumulate_sfchi0_q0
!!
!! FUNCTION
!! Update the spectral function of the independent particle susceptibility at q==0 for the contribution
!! of one pair of occupied-unoccupied band, for each frequency.
!! If symchi==1, the symmetries belonging to the little group of the external point q are used
!! to reconstruct the contributions in the full Brillouin zone. In this case, the equation implented is:
!!
!!  $ chi0(G1,G2,io)=chi0(G1,G2,io)+\sum_S (rhotwg(G1)*rhotwg^\dagger(G2))* \delta(\omega -trans) $
!!
!! where S is a symmetry belonging to the little group of q.
!! The subroutine also performs the symmetrization of the matrix elements of the
!! gradient operator and of the commutator [V_{nl},r] with the position operator.
!!
!! INPUTS
!!  ikbz=Index in the BZ of the k-point whose contribution to chi0 has to be added,
!!   if we use symmetries, the contribution to chi0 by this k-point has to be symmetrized.
!!  isym_kbz=Index of the symmetry such as k_bz = IS k_ibz
!!  itim_kbz=2 if time-reversal has to be used to obtain k_bz, 1 otherwise.
!!  my_wl,my_wr=min and Max frequency index treated by this processor.
!!  npwe=Number of plane waves used to describe chi0.
!!  npwepG0=Maximum number of G vectors to account for umklapps.
!!  nomega=Number of frequencies in the imaginary part.
!!  rhotwg(npwepG0)=Oscillator matrix elements corresponding to an occupied-unoccupied pair of states.
!!  rhotwx(3,nspinor**2)=Matrix elements of the gradient and of the commutator of the non-local operator with
!!    the position operator. The second term is present only if inclvkb=1,2.
!!  Gsph_epsG0<gsphere_t> Information on the "enlarged" G-sphere used for chi0, it contains umklapp G0 vectors
!!    %ng=number of G vectors in the enlarged sphere. It MUST be equal to the size of rhotwg
!!    %rottbm1(ng,2,nsym)=index of (IR)^{-1} G where I is the identity or the inversion
!!    %phmGt(ng,nsym)=phase factors associated to non-symmorphic operations
!!  Ltg_q<littlegroup_t_type>=Info on the little group associated to the external q-point.
!!    %timrev=2 it time-reversal is used, 1 otherwise
!!    %nsym_sg=Number of space group symmetries
!!    %wtksym(2,nsym,nkbz)=1 if the symmetry (with or without time-reversal) must be considered for this k-point
!!    %flag_umklp(timrev,nsym)= flag for umklapp processes
!!     if 1 that the particular operation (IS) requires a G_o to preserve Q, 0 otherwise
!! Cryst<crystal_t>=Info on unit cell and it symmetries
!!    %nsym=Number of symmetry operations.
!!    %symrec(3,3,nsym)=Symmetry operations in reciprocal space (reduced coordinates).
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  sf_chi0(npwe,npwe,my_wl:my_wr)=Updated spectral function at q==0.
!!  sf_lwing(npwe,my_wl:my_wr,3)=Updated lower wing of the spectral function.
!!  sf_uwing(npwe,mw_wl:my_wr,3)=Updated upper wing of the spectral function.
!!  sf_head(3,3,my_wl:my_wr)=Updated head of the spectral function.
!!
!! PARENTS
!!      cchi0q0
!!
!! CHILDREN
!!
!! SOURCE

subroutine accumulate_sfchi0_q0(ikbz,isym_kbz,itim_kbz,nspinor,symchi,npwepG0,npwe,Cryst,Ltg_q,Gsph_epsG0,&
& factocc,my_wl,iomegal,wl,my_wr,iomegar,wr,rhotwx,rhotwg,nomegasf,sf_chi0,sf_head,sf_lwing,sf_uwing)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikbz,my_wl,my_wr,nomegasf,npwe,npwepG0,nspinor
 integer,intent(in) :: isym_kbz,itim_kbz,symchi,iomegal,iomegar
 real(dp),intent(in) :: factocc,wl,wr
 type(littlegroup_t),intent(in) :: Ltg_q
 type(gsphere_t),target,intent(in) :: Gsph_epsG0
 type(crystal_t),intent(in) :: Cryst
!arrays
 complex(gwpc),intent(in) :: rhotwg(npwepG0)
 complex(gwpc),intent(in) :: rhotwx(3,nspinor**2)
 complex(gwpc),intent(inout) :: sf_chi0(npwe,npwe,my_wl:my_wr)
 complex(dpc),intent(inout) :: sf_head(3,3,my_wl:my_wr)
 complex(dpc),intent(inout) :: sf_lwing(npwe,my_wl:my_wr,3)
 complex(dpc),intent(inout) :: sf_uwing(npwe,my_wl:my_wr,3)

!Local variables-------------------------------
!scalars
 integer :: itim,isym,idir,jdir
 complex(gwpc) :: num
 character(len=500) :: msg
!arrays
 integer, ABI_CONTIGUOUS pointer :: Sm1G(:)
 complex(dpc) :: mir_kbz(3)
 complex(gwpc), ABI_CONTIGUOUS pointer :: phmGt(:)
 complex(gwpc),allocatable :: rhotwg_sym(:)

!************************************************************************

 if (iomegal<my_wl .or. iomegar>my_wr) then
   write(msg,'(3a,2(a,i0,a,i0,a))')ch10,&
&    'Indices out of boundary ',ch10,&
&    '  my_wl = ',my_wl,' iomegal = ',iomegal,ch10,&
&    '  my_wr = ',my_wr,' iomegar = ',iomegar,ch10
   MSG_BUG(msg)
 end if

 SELECT CASE (symchi)
 CASE (0)
   !
   ! Calculation without symmetries
   ! rhotwg(1)= R^-1q*rhotwx_ibz
   ! rhotwg(1)=-R^-1q*conjg(rhotwx_ibz) for inversion
   if (wl<huge(0.0_dp)*1.d-11) then
     !this is awful but it is still a first coding
     ! Num is single precision needed for cgerc check factocc
     num=-wl*factocc
     call XGERC(npwe,npwe,num,rhotwg,1,rhotwg,1,sf_chi0(:,:,iomegal),npwe)
   end if
   ! Last point, must accumulate left point but not the right one
   if (iomegar/=nomegasf+1 .and. wr<huge(0.0_dp)*1.d-11) then
     num=-wr*factocc
     call XGERC(npwe,npwe,num,rhotwg,1,rhotwg,1,sf_chi0(:,:,iomegar),npwe)
   end if

   ! ================================
   ! ==== Update heads and wings ====
   ! ================================

   ! Symmetrize <r> in full BZ: <Sk b|r|Sk b'> = R <k b|r|k b'> + \tau \delta_{bb'}
   if (nspinor == 1) then
      mir_kbz =(3-2*itim_kbz) * MATMUL(Cryst%symrec(:,:,isym_kbz),rhotwx(:,1))
    else
      mir_kbz = (3-2*itim_kbz) * MATMUL(Cryst%symrec(:,:,isym_kbz), sum(rhotwx(:,1:2), dim=2))
    end if
   if (itim_kbz==2) mir_kbz=CONJG(mir_kbz)

   do jdir=1,3
     if (wl<huge(0.0_dp)*1.d-11) then
       ! this is awful but it is still a first coding
       ! Num is single precision needed for cgerc check factocc
       num=-wl*factocc
       sf_uwing(:,iomegal,jdir) = sf_uwing(:,iomegal,jdir) + num * mir_kbz(jdir) * CONJG(rhotwg(1:npwepG0))
       sf_lwing(:,iomegal,jdir) = sf_lwing(:,iomegal,jdir) + num * rhotwg(1:npwepG0) * CONJG(mir_kbz(jdir))
       do idir=1,3
         sf_head(idir,jdir,iomegal) = sf_head(idir,jdir,iomegal) + num * mir_kbz(idir) * CONJG(mir_kbz(jdir))
       end do
     end if

     ! Last point, must accumulate left point but not the right one
     if (iomegar/=nomegasf+1 .and. wr<huge(0.0_dp)*1.d-11) then
       num=-wr*factocc
       sf_uwing(:,iomegar,jdir) = sf_uwing(:,iomegar,jdir) + num * mir_kbz(jdir) * CONJG(rhotwg(1:npwepG0))
       sf_lwing(:,iomegar,jdir) = sf_lwing(:,iomegar,jdir) + num * rhotwg(1:npwepG0) * CONJG(mir_kbz(jdir))
       do idir=1,3
         sf_head(idir,jdir,iomegar) = sf_head(idir,jdir,iomegar) + num * mir_kbz(idir) * CONJG(mir_kbz(jdir))
       end do
     end if
   end do ! jdir

 CASE (1)
   ! === Notes on the symmetrization of oscillator matrix elements ===
   ! If  Sq = q then  M_G( Sk,q)= e^{-i(q+G)\cdot t} M_{ S^-1G}  (k,q)
   ! If -Sq = q then  M_G(-Sk,q)= e^{-i(q+G)\cdot t} M_{-S^-1G}^*(k,q)
   !
   ! In case of an umklapp process
   ! If  Sq = q+G_o then  M_G( Sk,q)= e^{-i(q+G)\cdot t} M_{ S^-1(G-G_o}   (k,q)
   ! If -Sq = q+G_o then  M_G(-Sk,q)= e^{-i(q+G)\cdot t} M_{-S^-1(G-G-o)}^*(k,q)
   !
   ! rhotwg(1)= R^-1q*rhotwx_ibz
   ! rhotwg(1)=-R^-1q*conjg(rhotwx_ibz) for inversion

   ABI_MALLOC(rhotwg_sym,(npwe))

   ! Loop over symmetries of the space group and time-reversal
   do isym=1,Ltg_q%nsym_sg
     do itim=1,Ltg_q%timrev

       if (Ltg_q%wtksym(itim,isym,ikbz)==1) then
         ! This operation belongs to the little group and has to be considered to reconstruct the BZ
         phmGt => Gsph_epsG0%phmGt(1:npwe,isym) ! In these 2 lines mind the slicing (1:npwe)
         Sm1G  => Gsph_epsG0%rottbm1(1:npwe,itim,isym)

         SELECT CASE (itim)
         CASE (1)
           rhotwg_sym(1:npwe)=rhotwg(Sm1G(1:npwe))*phmGt(1:npwe)
         CASE (2)
           rhotwg_sym(1:npwe)=CONJG(rhotwg(Sm1G(1:npwe)))*phmGt(1:npwe)
         CASE DEFAULT
           MSG_BUG(sjoin('Wrong value of itim:', itoa(itim)))
         END SELECT

         ! Multiply elements G,Gp of rhotwg_sym*num and accumulate in sf_chi0(G,Gp,io)
         if (wl<huge(0.0_dp)*1.d-11) then
           num=-wl*factocc
           call XGERC(npwe,npwe,num,rhotwg_sym,1,rhotwg_sym,1,sf_chi0(:,:,iomegal),npwe)
         end if

         ! Last point, must accumulate left point but not the right one
         if (iomegar/=nomegasf+1 .and. wr<huge(0.0_dp)*1.d-11) then
           num=-wr*factocc
           call XGERC(npwe,npwe,num,rhotwg_sym,1,rhotwg_sym,1,sf_chi0(:,:,iomegar),npwe)
         end if

         ! Accumulate heads and wings.
         ! Symmetrize <r> in full BZ: <Sk b|r|Sk b'> = R <k b|r|k b'> + \tau \delta_{bb'}
         if (nspinor == 1) then
           mir_kbz =(3-2*itim) * MATMUL(Cryst%symrec(:,:,isym),rhotwx(:,1))
         else
           mir_kbz = (3-2*itim) * MATMUL(Cryst%symrec(:,:,isym), sum(rhotwx(:,1:2), dim=2))
         end if
         if (itim==2) mir_kbz=CONJG(mir_kbz)

         do jdir=1,3
           if (wl<huge(0.0_dp)*1.d-11) then
             ! this is awful but it is still a first coding
             ! Num is single precision needed for cgerc check factocc
             num=-wl*factocc
             sf_uwing(:,iomegal,jdir) = sf_uwing(:,iomegal,jdir) + num * mir_kbz(jdir) * CONJG(rhotwg_sym(1:npwe))
             sf_lwing(:,iomegal,jdir) = sf_lwing(:,iomegal,jdir) + num * rhotwg_sym(1:npwe) * CONJG(mir_kbz(jdir))
             do idir=1,3
               sf_head(idir,jdir,iomegal) = sf_head(idir,jdir,iomegal) + num * mir_kbz(idir) * CONJG(mir_kbz(jdir))
             end do
           end if

           ! Last point, must accumulate left point but not the right one
           if (iomegar/=nomegasf+1 .and. wr<huge(0.0_dp)*1.d-11) then
             num=-wr*factocc
             sf_uwing(:,iomegar,jdir) = sf_uwing(:,iomegar,jdir) + num * mir_kbz(jdir) * CONJG(rhotwg_sym(1:npwe))
             sf_lwing(:,iomegar,jdir) = sf_lwing(:,iomegar,jdir) + num * rhotwg_sym(1:npwe) * CONJG(mir_kbz(jdir))
             do idir=1,3
               sf_head(idir,jdir,iomegar) = sf_head(idir,jdir,iomegar) + num * mir_kbz(idir) * CONJG(mir_kbz(jdir))
             end do
           end if
         end do ! jdir

       end if !wtksym
     end do !inv
   end do !isym
   ABI_FREE(rhotwg_sym)

 CASE DEFAULT
   MSG_BUG(sjoin('Wrong value of symchi:', itoa(symchi)))
 END SELECT

end subroutine accumulate_sfchi0_q0
!!***

!----------------------------------------------------------------------

!!****f* m_chi0tk/assemblychi0sf
!! NAME
!! assemblychi0sf
!!
!! FUNCTION
!! Update the spectral function of the irreducible polarizability for the contribution
!! of one pair of occupied-unoccupied states, for each frequency.
!! If symchi==1, the symmetries of the little group of the external q-point are used
!! to symmetrize the contribution in the full Brillouin zone. In this case, the routine computes:
!!
!!   $ chi0(G1,G2,io)=chi0(G1,G2,io)+\sum_S (rhotwg(G1)*rhotwg^\dagger(G2))*\delta(w - trans) $
!!
!! where S are the symmetries of the little group of the external q-point.
!!
!! INPUTS
!!  ik_bz=Index of the k-point in the BZ whose contribution has to be added to the spectral function of chi0
!!    If symchi=1, the contribution is symmetrized.
!!  my_wl,my_wr=min and Max frequency index treated by this processor.
!!  npwe=Number of plane waves used to describe chi0.
!!  npwepG0=Maximum number of G vectors taking into account umklapp vectors.
!!  nomegasf=Number of frequencies for the spectral function.
!!  nspinor=Number of spinorial components.
!!  symchi=1 if symmetries are used, 0 otherwise
!!  rhotwg(npwepG0)=Oscillator matrix elements corresponding to an occupied-unoccupied pair of states.
!!  timrev=if 2, time reversal has to be used to obtain k_bz; 1 otherwise.
!!  Gsph_epsG0<gsphere_t> Information on the "enlarged" G-sphere used for chi0, it contains umklapp G0 vectors
!!    %ng=number of G vectors in the enlarged sphere, actually MUST be equal to the size of rhotwg
!!    %rottbm1(ng,2,nsym)=index of (IR)^{-1} G where I is the identity or the inversion
!!    %phmGt(ng,nsym)=phase factors associated to non-simmorphic operations
!!  Ltg_q<littlegroup_t_type>=Info on the little group associated to the external q-point.
!!    %timrev=2 it time-reversal is used, 1 otherwise
!!    %nsym_sg=Number of space group symmetries
!!    %wtksym(2,nsym,nkbz)=1 if the symmetry (with or without time-reversal) must be considered for this k-point
!!    %flag_umklp(timrev,nsym)= flag for umklapp processes
!!      if 1 that the particular operation (IS) requires a G_o to preserve Q, 0 otherwise
!!    %igmG0(npwepG0,timrev,nsym) index of G-G0 in the array gvec
!!  factocc=occupation factor=f_occ*(ockp-occk) (see cchi0.F90)
!!  wl,wr=Weights used to approximate the delta function.
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  chi0sf(npwe,npwe,my_wl:my_wr)= updated spectral function.
!!
!! NOTES
!!  Umklapp processes are not yet implemented
!!
!! PARENTS
!!      cchi0
!!
!! CHILDREN
!!
!! SOURCE

subroutine assemblychi0sf(ik_bz,symchi,Ltg_q,npwepG0,npwe,rhotwg,Gsph_epsG0,&
& factocc,my_wl,iomegal,wl,my_wr,iomegar,wr,nomegasf,chi0sf)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz,iomegal,iomegar,my_wl,my_wr,nomegasf,npwe,npwepG0
 integer,intent(in) :: symchi
 real(dp),intent(in) :: factocc,wl,wr
 type(gsphere_t),intent(in) :: Gsph_epsG0
 type(littlegroup_t),intent(in) :: Ltg_q
!arrays
 complex(gwpc),intent(in) :: rhotwg(npwepG0)
 complex(gwpc),intent(inout) :: chi0sf(npwe,npwe,my_wl:my_wr)

!Local variables-------------------------------
!scalars
 integer :: isym,itim,ig1,ig2
 complex(gwpc) :: num
 character(len=500) :: msg
!arrays
 integer :: Sm1_gmG0(npwe)
 complex(gwpc) :: rhotwg_sym(npwe)


! *************************************************************************

 if (iomegal < my_wl .or. iomegar > my_wr) then
   write(msg,'(3a,2(a,i0,a,i0,a))')ch10,&
&    ' Indices out of boundary ',ch10,&
&    '  my_wl = ',my_wl,' iomegal = ',iomegal,ch10,&
&    '  my_wr = ',my_wr,' iomegar = ',iomegar,ch10
   MSG_BUG(msg)
 end if

 SELECT CASE (symchi)
 CASE (0)
    ! Do not use symmetries.

! MG: This is the best I can do for this part.
!$omp PARALLEL private(num)
!$omp SECTIONS
!$omp SECTION
   if (wl<huge(0.0_dp)*1.d-11) then !FIXME this is awful
     num=-wl*factocc
     call XGERC(npwe,npwe,num,rhotwg,1,rhotwg,1,chi0sf(:,:,iomegal),npwe)
   end if

   ! Last point, must accumulate left point but not the right one
!$omp SECTION
   if (iomegar/=nomegasf+1 .and. wr<huge(0.0_dp)*1.d-11) then
     num=-wr*factocc
     call XGERC(npwe,npwe,num,rhotwg,1,rhotwg,1,chi0sf(:,:,iomegar),npwe)
   end if
!$omp end SECTIONS
!$omp end PARALLEL

 CASE (1)
   ! Use symmetries to reconstruct oscillator matrix elements
   ! Notes on the symmetrization of the oscillator maxtri elements:
   !
   ! If  Sq=q then  M_G^( Sk,q)= e^{-i(q+G)\cdot t} M_{ S^-1G}  (k,q)
   ! If -Sq=q then  M_G^(-Sk,q)= e^{-i(q+G)\cdot t} M_{-S^-1G}^*(k,q)
   !
   ! In case of an umklapp process
   ! If  Sq=q+G_o then  M_G( Sk,q)= e^{-i(q+G)\cdot t} M_{ S^-1(G-G_o}   (k,q)
   ! If -Sq=q+G_o then  M_G(-Sk,q)= e^{-i(q+G)\cdot t} M_{-S^-1(G-G_o)}^*(k,q)
   !
   ! Ltg_q%igmG0(ig,itim,isym) contains the index of G-G0 where ISq=q+G0
   ! Note that there is no need to take into account the phases due to q,
   ! They cancel in the scalar product ==> phmGt(G,isym)=e^{-iG\cdot t}
   !
   ! Mind the slicing of %rottbm1(npwepG0,timrev,nsym) and %phgt(npwepG0,nsym) as
   ! these arrays, usually, do not conform to rho_twg_sym(npw) !
   !
   !ABI_MALLOC(rhotwg_sym,(npwe))
   !
   ! Loop over symmetries of the space group and time-reversal
   do isym=1,Ltg_q%nsym_sg
     do itim=1,Ltg_q%timrev

       if (Ltg_q%wtksym(itim,isym,ik_bz)==1) then
         ! This operation belongs to the little group and has to be used to reconstruct BZ.
         ! TODO this is a hot-spot, should add a test on the umklapp
         !
         ! In these 3 lines mind the slicing (1:npwe)
         Sm1_gmG0(1:npwe)=Gsph_epsG0%rottbm1( Ltg_q%igmG0(1:npwe,itim,isym), itim,isym)

         SELECT CASE (itim)
         CASE (1)
           rhotwg_sym(1:npwe)=rhotwg(Sm1_gmG0(1:npwe)) * Gsph_epsG0%phmGt(1:npwe,isym)
         CASE (2)
           rhotwg_sym(1:npwe)=CONJG(rhotwg(Sm1_gmG0(1:npwe))) * Gsph_epsG0%phmGt(1:npwe,isym)
         CASE DEFAULT
           MSG_BUG(sjoin('Wrong value for itim:', itoa(itim)))
         END SELECT

#if 0
!! MG: This is the best I can do, at present.
!$omp PARALLEL private(num)
!$omp SECTIONS

!$omp SECTION
         if (wl<huge(0.0_dp)*1.d-11) then
           num=-wl*factocc
           call XGERC(npwe,npwe,num,rhotwg_sym,1,rhotwg_sym,1,chi0sf(:,:,iomegal),npwe)
         end if
!$omp SECTION
         !
         ! Last point, must accumulate left point but not the right one
         if (iomegar/=nomegasf+1 .and. wr<huge(0.0_dp)*1.d-11) then
           num=-wr*factocc
           call XGERC(npwe,npwe,num,rhotwg_sym,1,rhotwg_sym,1,chi0sf(:,:,iomegar),npwe)
         end if
!$omp end SECTIONS
!$omp end PARALLEL
#else

         if (wl<huge(0.0_dp)*1.d-11) then
           !call XGERC(npwe,npwe,num,rhotwg_sym,1,rhotwg_sym,1,chi0sf(:,:,iomegal),npwe)
           num=-wl*factocc
!$omp parallel do
           do ig2=1,npwe
             do ig1=1,npwe
               chi0sf(ig1,ig2,iomegal) = chi0sf(ig1,ig2,iomegal) + num * rhotwg_sym(ig1) * CONJG(rhotwg_sym(ig2))
             end do
           end do
         end if

         ! Last point, must accumulate left point but not the right one
         if (iomegar/=nomegasf+1 .and. wr<huge(0.0_dp)*1.d-11) then
           !call XGERC(npwe,npwe,num,rhotwg_sym,1,rhotwg_sym,1,chi0sf(:,:,iomegar),npwe)
           num=-wr*factocc
!$omp parallel do
           do ig2=1,npwe
             do ig1=1,npwe
               !call XGERC(npwe,npwe,num,rhotwg_sym,1,rhotwg_sym,1,chi0sf(:,:,iomegal),npwe)
               chi0sf(ig1,ig2,iomegar) = chi0sf(ig1,ig2,iomegar) + num * rhotwg_sym(ig1) * CONJG(rhotwg_sym(ig2))
             end do
           end do
         end if
#endif
       end if !wtksym

     end do !inv
   end do !isym
   !ABI_FREE(rhotwg_sym)

 CASE DEFAULT
   MSG_BUG(sjoin('Wrong value for symchi:', itoa(symchi)))
 END SELECT

end subroutine assemblychi0sf
!!***

!----------------------------------------------------------------------

!!****f* m_chi0tk/approxdelta
!! NAME
!!  approxdelta
!!
!! FUNCTION
!!  Approximate the Dirac function using two methods:
!!  method 1) a triangular funtion centered at the value egwdiff_re, Eq 17 of PRB 74, 035101 (2006) [[cite:Shishkin2006]]
!!  method 2) a gaussian of witdth ep%spsmear expandended in Taylor series
!!  (at the moment only the 0-th moments)
!!
!!  Subroutine needed to implement the calculation
!!  of the polarizability using the spectral representation as proposed in:
!!  PRB 74, 035101 (2006) [[cite:Shishkin2006]]
!!  and PRB 61, 7172 (2000) [[cite:Miyake2000]]
!!
!! INPUTS
!!  nomegasf=number of frequencies in the grid for Im \chi_0
!!  omegasf(0:nomega+1)= frequencies (real)
!!  egwdiff_re = transition energy where the delta function is centered
!!
!!  method= 1: a triangular shaped function used to approximated the delta
!!          2: gaussian approximation with standard deviation (smear)
!! smear= used only in case of method==2, defines the width of the gaussian
!!
!! OUTPUT
!!  wl = weight associated to omegal (last omega wich is smaller than egwdiff_re
!!  wr = weight associate to omegar  (first omega larger than egwdff_re
!!  iomegal= index in the array omegasf of the last frequency < egwdiff
!!  iomegar= index in the array omegasf of the first frequency > egwdiff
!!
!! PARENTS
!!      cchi0,cchi0q0
!!
!! CHILDREN
!!
!! SOURCE

subroutine approxdelta(nomegasf,omegasf,egwdiff_re,smear,iomegal,iomegar,wl,wr,spmeth)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nomegasf,spmeth
 integer,intent(out) :: iomegal,iomegar
 real(dp),intent(in) :: egwdiff_re,smear
 real(dp),intent(out) :: wl,wr
!arrays
 real(dp),intent(in) :: omegasf(nomegasf)

!Local variables-------------------------------
 integer :: io,iomega
 real(dp) :: omegal,omegar,deltal,deltar
 !character(len=500) :: msg
! *************************************************************************

 iomega=-999
 do io=nomegasf,1,-1
   if (omegasf(io)<egwdiff_re) then
    iomega=io; EXIT
   end if
 end do

 iomegal=iomega   ; omegal=omegasf(iomegal)
 iomegar=iomegal+1; omegar=omegasf(iomegar)

 SELECT CASE (spmeth)
 CASE (1)
   ! Weights for triangular shaped function
   wr=  (egwdiff_re-omegal)/(omegar-omegal)
   wl= -(egwdiff_re-omegar)/(omegar-omegal)

 CASE (2)
   ! Weights for gaussian method (0-th moment)
   deltal=(egwdiff_re-omegal)/smear
   deltar=(omegar-egwdiff_re)/smear
   if (deltar>=deltal) then
     wl=EXP(-deltal*deltal)
     ! this value is used to avoid double counting and speed-up
     wr=huge(one)*1.d-10
   else
     wl=huge(one)*1.d-10
     wr=exp(-deltal*deltal)
   end if

 CASE DEFAULT
   MSG_BUG(sjoin('Wrong value for spmeth:', itoa(spmeth)))
 END SELECT

end subroutine approxdelta
!!***

!----------------------------------------------------------------------

!!****f* m_chi0tk/calc_kkweight
!! NAME
!!  calc_kkweight
!!
!! FUNCTION
!!  Calculate frequency dependent weights needed to perform the Hilbert transform
!!
!!  Subroutine needed to implement the calculation
!!  of the polarizability using the spectral representation as proposed in:
!!  PRB 74, 035101 (2006) [[cite:Shishkin2006]]
!!  and PRB 61, 7172 (2000) [[cite:Miyake2000]]
!!
!! INPUTS
!! nsp=number of frequencies where the imaginary part of the polarizability is evaluated
!! ne=number of frequencies for the polarizability (same as in epsilon^-1)
!! omegasp(nsp)=real frequencies for the imaginary part of the polarizability
!! omegae(ne)= imaginary frequencies for the polarizability
!! delta=small imaginary part used to avoid poles, input variables
!!
!! OUTPUT
!! kkweight(nsp,ne)=frequency dependent weights Eq A1 PRB 74, 035101 (2006) [[cite:Shishkin2006]]
!!
!! PARENTS
!!      m_chi0
!!
!! CHILDREN
!!
!! SOURCE
!!

subroutine calc_kkweight(ne,omegae,nsp,omegasp,delta,omegamax,kkw)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ne,nsp
 real(dp),intent(in) :: delta,omegamax
!arrays
 real(dp),intent(in) :: omegasp(nsp)
 complex(dpc),intent(in) :: omegae(ne)
 complex(dpc),intent(out) :: kkw(nsp,ne)

!Local variables-------------------------------
!scalars
 integer :: isp,je
 real(dp) :: eta,xx1,xx2,den1,den2
 complex(dpc) :: c1,c2,wt
!************************************************************************

 DBG_ENTER("COLL")

 kkw(:,:)=czero

 do je=1,ne
   eta=delta
   wt=omegae(je)
   ! Not include shift at omega==0, what about metallic systems?
   if (abs(real(omegae(je)))<tol6 .and. abs(aimag(wt))<tol6) eta=tol12
   !  Not include shift along the imaginary axis
   if (abs(aimag(wt))>tol6) eta=zero
   do isp=1,nsp
     if (isp==1) then
       ! Skip negative point, should check that this would not lead to spurious effects
       c1=czero
       den1=one
     else
       xx1=omegasp(isp-1)
       xx2=omegasp(isp)
       den1= xx2-xx1
       c1= -(wt-xx1+j_dpc*eta)*log( (wt-xx2+j_dpc*eta)/(wt-xx1+j_dpc*eta) )&
&          +(wt+xx1-j_dpc*eta)*log( (wt+xx2-j_dpc*eta)/(wt+xx1-j_dpc*eta) )
       c1= c1/den1
     end if
     xx1=omegasp(isp)
     if (isp==nsp) then
       ! Skip last point should check that this would not lead to spurious effects
       xx2=omegamax
     else
       xx2=omegasp(isp+1)
     end if
     den2=xx2-xx1
     c2=  (wt-xx2+j_dpc*eta)*log( (wt-xx2+j_dpc*eta)/(wt-xx1+j_dpc*eta) )&
&        -(wt+xx2-j_dpc*eta)*log( (wt+xx2-j_dpc*eta)/(wt+xx1-j_dpc*eta) )
     c2= c2/den2
     kkw(isp,je)=  c1/den1 + c2/den2
   end do
 end do

 DBG_EXIT("COLL")

end subroutine calc_kkweight
!!***

!----------------------------------------------------------------------

!!****f* m_chi0tk/setup_spectral
!! NAME
!!  setup_spectral
!!
!! FUNCTION
!! Calculation of \chi_o based on the spectral method as proposed in PRB 74, 035101 (2006) [[cite:Shishkin2006]]
!! and PRB 61, 7172 (2000) [[cite:Miyake2000]].
!! Setup of the real frequency mesh for $\Im\chi_o$ and of the frequency-dependent weights for
!! Hilbert transform. Note that CPU time does not depend dramatically on nomegasf unlike memory.
!! spmeth defines the approximant for the delta function:
!!  ==1 : use Triangular approximant (Kresse method)
!!  ==2 : use Gaussian method, requiring smearing (Miyake method)
!!
!! INPUTS
!! nomegasf=number of points for the imaginary part of $\chi0(q,\omega)$
!! nomega=number of frequencies in $\chi0(q,\omega)$.
!! max_rest,min_res=max and min resonant transition energy (for this q-point)
!! my_max_rest,my_min_rest=max and min resonant transition energy treated by this processor
!! method=integer flag defining the type of frequency mesh used for $\Im chi0$
!!  | 0 for a linear mesh
!!  | 1 for a mesh densified around omegaplasma
!! omegaplasma=frequency around which the mesh is densifies (usually Drude plasma frequency)
!!  used only in case of method==1
!! zcut=small imaginary shift to avoid pole in chi0
!!
!! OUTPUT
!!  kkweight(nomegasf,nomega)=Frequency dependent weight for Hilber transform.
!!  omegasf(nomegasf+1)=frequencies for imaginary part.
!!
!! PARENTS
!!      cchi0,cchi0q0
!!
!! CHILDREN
!!
!! SOURCE

subroutine setup_spectral(nomega,omega,nomegasf,omegasf,max_rest,min_rest,my_max_rest,my_min_rest,&
&  method,zcut,omegaplasma,my_wl,my_wr,kkweight)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: method,nomega,nomegasf
 integer,intent(out) :: my_wl,my_wr
 real(dp),intent(in) :: max_rest,min_rest,omegaplasma,zcut
 real(dp),intent(in) :: my_max_rest,my_min_rest
!arrays
 real(dp),intent(out) :: omegasf(nomegasf)
 complex(dpc),intent(in) :: omega(nomega)
 complex(dpc),intent(out) :: kkweight(nomegasf,nomega)

!Local variables-------------------------------
!scalars
 integer :: io,ii
 real(dp) :: nu_min,nu_max,nu1,nu2,dd,domegasf,wp,deltat
 character(len=500) :: msg
!arrays
 integer,allocatable :: insort(:)
!************************************************************************

 ! The mesh must enclose the entire range of transitions.
 dd=(max_rest-min_rest)/(nomegasf-1)
 domegasf=(max_rest-min_rest+2*dd)/(nomegasf-1)

 write(msg,'(4a,f8.3,3a,i5,2a,f8.5,a)')ch10,&
&  ' === Info on the real frequency mesh for spectral method === ',ch10,&
&  '  maximum frequency = ',max_rest*Ha_eV,' [eV]',ch10,&
&  '  nomegasf = ',nomegasf,ch10,&
&  '  domegasf = ',domegasf*Ha_eV,' [eV]'
 call wrtout(std_out,msg,'COLL')

 if (min_rest<tol6) then
   MSG_WARNING("System seems to be metallic")
 end if

 ! ======================================================
 ! === Setup of the w-mesh for the spectral function ====
 ! ======================================================
 SELECT CASE (method)
 CASE (0)
   ! Linear mesh.
   call wrtout(std_out,' Using linear mesh for Im chi0','COLL')
   do io=1,nomegasf
     omegasf(io)=(io-1)*domegasf+min_rest-dd
   end do

 CASE (1)
   ! Non-homogeneous mesh densified around omega_plasma, do not improve results ===
   ! WARNING_ this part has to be checked since I modified omegasf
   write(msg,'(a,f7.4,a)')' Using mesh densified around ',omegaplasma*Ha_eV,' [eV] '
   call wrtout(std_out,msg,'COLL')
   wp=omegaplasma ; deltat=max_rest-min_rest
   nu_min=zero
   if (deltat<wp ) then
     nu_max = wp/sqrt2 *   ATAN(sqrt2*deltat*wp/(-deltat**2+wp**2))
   else
     nu_max = wp/sqrt2 * ( ATAN(sqrt2*deltat*wp/(-deltat**2+wp**2)) + pi)
   end if
   domegasf=(nu_max-nu_min)/(nomegasf+1)
   !write(std_out,*)  -(wp/sqrt2) * atan(sqrt2*deltat*wp/(deltat**2-wp**2))
   omegasf(1)=zero ; omegasf(nomegasf+1)=deltat
   ii=0
   do io=2,nomegasf
     nu1=domegasf*(io-1) ; nu2=TAN(-sqrt2*nu1/wp)
     if (nu2<0) then
       omegasf(io) = wp * (one - SQRT(1+2*nu2**2))/(sqrt2*nu2)
     else
       omegasf(io) = wp * (one + SQRT(1+2*nu2**2))/(sqrt2*nu2)
     end if
     if (omegasf(io)> deltat ) then
       omegasf(io)= deltat-0.1*ii
       ii=ii+1
     end if
     ! write(102,'(i4,2x,3(f9.4,2x))')io,nu1,nu2,ep%omegasf(io)*Ha_eV
   end do

   ! Reorder frequencies in ascending order
   ABI_MALLOC(insort,(nomegasf+1))
   insort(:)=(/ (io,io=1,nomegasf+1) /)
   call sort_dp(nomegasf+1,omegasf,insort,tol14)
   ABI_FREE(insort)

 CASE DEFAULT
   MSG_BUG(sjoin('Wrong value for method:', itoa(method)))
 END SELECT
 !write(std_out,*)omegasf(1)*Ha_eV,omegasf(nomegasf)*Ha_eV

 ! Find min and max index in omegasf treated by this processor.
 my_wr=-999
 do io=1,nomegasf
   if (omegasf(io)>my_max_rest) then
     my_wr=io; EXIT
   end if
 end do
 if (my_wr==nomegasf+2) my_wr=nomegasf+1
 my_wl=-999
 do io=nomegasf,1,-1
   if (omegasf(io)< my_min_rest) then ! Check metals
     my_wl=io; EXIT
   end if
 end do

 write(msg,'(a,2(1x,i0))')' my_wl and my_wr:',my_wl,my_wr
 call wrtout(std_out,msg,'PERS')

 if (my_wl==-999 .or. my_wr==-999) then
   write(msg,'(a,2i6)')' wrong value in my_wl and/or my_wr ',my_wl,my_wr
   MSG_ERROR(msg)
 end if

 ! Calculate weights for Hilbert transform.
 call calc_kkweight(nomega,omega,nomegasf,omegasf,zcut,max_rest,kkweight)

end subroutine setup_spectral
!!***

!----------------------------------------------------------------------

!!****f* m_chi0tk/hilbert_transform
!! NAME
!!  hilbert_transform
!!
!! FUNCTION
!!  Compute the hilbert transform.
!!
!! INPUTS
!! nomegasf=number of points for the imaginary part of $\chi0(q,\omega)$
!! nomega=number of frequencies in $\chi0(q,\omega)$.
!! max_rest,min_res=max and min resonant transition energy (for this q-point)
!! my_max_rest,my_min_rest=max and min resonant transition energy treated by this processor
!!
!! OUTPUT
!!
!! PARENTS
!!      cchi0,cchi0q0
!!
!! CHILDREN
!!
!! SOURCE

subroutine hilbert_transform(npwe,nomega,nomegasf,my_wl,my_wr,kkweight,sf_chi0,chi0,spmeth)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spmeth,nomega,nomegasf,my_wl,my_wr,npwe
!arrays
 complex(dpc),intent(in) :: kkweight(nomegasf,nomega)
 complex(gwpc),intent(inout) :: sf_chi0(npwe,npwe,my_wl:my_wr)
 complex(gwpc),intent(inout) :: chi0(npwe,npwe,nomega)

!Local variables-------------------------------
!scalars
 integer :: ig2,my_nwp
 character(len=500) :: msg
!arrays
 complex(gwpc),allocatable :: A_g1wp(:,:),H_int(:,:),my_kkweight(:,:)

!************************************************************************

#ifdef HAVE_OPENMP
 write(msg,'(2a,i3,a)')ch10,' Performing Hilbert transform (with OpenMP) using method ',spmeth,' It might take some time...'
#else
 write(msg,'(2a,i3,a)')ch10,' Performing Hilbert transform using method ',spmeth,' It might take some time...'
#endif
 call wrtout(std_out, msg, do_flush=.True.)

 my_nwp = my_wr - my_wl +1

!$omp parallel private(my_kkweight, A_g1wp, H_int, ig2)
 ABI_MALLOC(my_kkweight, (my_wl:my_wr,nomega))
 my_kkweight = kkweight(my_wl:my_wr,:)

 ABI_MALLOC(A_g1wp, (npwe, my_nwp))
 ABI_MALLOC(H_int, (npwe, nomega))

!$omp do
 do ig2=1,npwe
   A_g1wp = sf_chi0(:,ig2,:)

   ! Compute H_int = MATMUL(A_g1wp,my_kkweight)
   call XGEMM('N','N',npwe,nomega,my_nwp,cone_gw,A_g1wp,npwe,my_kkweight,my_nwp,czero_gw,H_int,npwe)
   chi0(:,ig2,:) = H_int
 end do

 ABI_FREE(my_kkweight)
 ABI_FREE(A_g1wp)
 ABI_FREE(H_int)
!$omp end parallel

end subroutine hilbert_transform
!!***

!----------------------------------------------------------------------

!!****f* m_chi0tk/hilbert_transform_headwings
!! NAME
!!  hilbert_transform_headwings
!!
!! FUNCTION
!!  Compute the hilbert transform the heads and wings of the polarizability.
!!
!! INPUTS
!! nomegasf=number of points for the imaginary part of $\chi0(q,\omega)$
!! nomega=number of frequencies in $\chi0(q,\omega)$.
!! max_rest,min_res=max and min resonant transition energy (for this q-point)
!! my_max_rest,my_min_rest=max and min resonant transition energy treated by this processor
!!
!! OUTPUT
!!
!! PARENTS
!!      cchi0q0
!!
!! CHILDREN
!!
!! SOURCE

subroutine hilbert_transform_headwings(npwe,nomega,nomegasf,my_wl,my_wr,kkweight, &
& sf_lwing,sf_uwing,sf_head,chi0_lwing,chi0_uwing,chi0_head,spmeth)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spmeth,nomega,nomegasf,my_wl,my_wr,npwe
!arrays
 complex(dpc),intent(in) :: kkweight(nomegasf,nomega)
 complex(dpc),intent(inout) :: sf_lwing(npwe,my_wl:my_wr,3)
 complex(dpc),intent(inout) :: sf_uwing(npwe,my_wl:my_wr,3)
 complex(dpc),intent(inout) :: sf_head(3,3,my_wl:my_wr)
 complex(dpc),intent(inout) :: chi0_lwing(npwe,nomega,3)
 complex(dpc),intent(inout) :: chi0_uwing(npwe,nomega,3)
 complex(dpc),intent(inout) :: chi0_head(3,3,nomega)

!Local variables-------------------------------
!scalars
 integer :: ig1,idir,io,iw
 complex(dpc) :: kkw
 character(len=500) :: msg
!************************************************************************

#ifdef HAVE_OPENMP
 write(msg,'(2a,i3,a)')ch10,' Performing Hilbert transform (with OpenMP) using method ',spmeth,' It might take some time...'
#else
 write(msg,'(2a,i3,a)')ch10,' Performing Hilbert transform using method ',spmeth,' It might take some time...'
#endif
 call wrtout(std_out,msg,'COLL',do_flush=.True.)

 ! Hilbert transform of the head.
 do io=1,nomega
   chi0_head(1,1,io) = SUM(kkweight(my_wl:my_wr,io)*sf_head(1,1,my_wl:my_wr))
   chi0_head(2,1,io) = SUM(kkweight(my_wl:my_wr,io)*sf_head(2,1,my_wl:my_wr))
   chi0_head(3,1,io) = SUM(kkweight(my_wl:my_wr,io)*sf_head(3,1,my_wl:my_wr))
   chi0_head(1,2,io) = SUM(kkweight(my_wl:my_wr,io)*sf_head(1,2,my_wl:my_wr))
   chi0_head(2,2,io) = SUM(kkweight(my_wl:my_wr,io)*sf_head(2,2,my_wl:my_wr))
   chi0_head(3,2,io) = SUM(kkweight(my_wl:my_wr,io)*sf_head(3,2,my_wl:my_wr))
   chi0_head(1,3,io) = SUM(kkweight(my_wl:my_wr,io)*sf_head(1,3,my_wl:my_wr))
   chi0_head(2,3,io) = SUM(kkweight(my_wl:my_wr,io)*sf_head(2,3,my_wl:my_wr))
   chi0_head(3,3,io) = SUM(kkweight(my_wl:my_wr,io)*sf_head(3,3,my_wl:my_wr))
 end do

 ! Hilbert transform for wings.
 ! Partial contributions to chi0 will be summed afterwards.
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(kkw)
 do idir=1,3
   do io=1,nomega
     do iw=my_wl,my_wr
       kkw = kkweight(iw,io)
       do ig1=1,npwe
         chi0_lwing(ig1,io,idir) = chi0_lwing(ig1,io,idir) + kkw*sf_lwing(ig1,iw,idir)
         chi0_uwing(ig1,io,idir) = chi0_uwing(ig1,io,idir) + kkw*sf_uwing(ig1,iw,idir)
       end do
     end do
   end do
 end do  !idir

end subroutine hilbert_transform_headwings
!!***

!----------------------------------------------------------------------

!!****f* m_chi0tk/completechi0_deltapart
!! NAME
!! completechi0_deltapart
!!
!! FUNCTION
!!  Apply the delta part of the completeness correction to chi0
!!
!! INPUTS
!!  ik_bz=Index of the k-point in the full BZ whose contribution has to be added and symmetrized.
!!  qzero=.TRUE. is long wave-length limit.
!!  symchi=1 if we are summing over IBZ_q and symmetrization has to be performed.
!!  npwe=Number of G vectors in chi0.
!!  npwvec=MAX number of G.
!!  nomega=Number of frequencies.
!!  nspinor=Number of spinorial components.
!!  nfftot=Total Number of points in the FFT
!!  ngfft(18)=Info on the FFT.
!!  igfft0(npwvec)=Index of each G in the FFT array.
!!  Gsph_FFT=<gsphere_t>=Info on the largest G-sphere contained in the FFT box used for wavefunctions.
!!  Ltg_q=<littlegroup_t>= Structure gathering information on the little group of the external q.
!!  green_enhigh_w=Approximated frequency dependent part of the Green function entering equation (TODO put reference)
!!  wfwfg=Fourier components of u_{kb1}.u_{kb2}
!!
!! OUTPUT
!!  See SIDES EFFECTS
!!
!! SIDES EFFECTS
!!  chi0(npwe,npwe,nomega)= In input chi0 calculated so far,
!!  In output the "delta part" of the completeness correction is added.
!!
!! PARENTS
!!      cchi0,cchi0q0
!!
!! CHILDREN
!!
!! SOURCE

subroutine completechi0_deltapart(ik_bz,qzero,symchi,npwe,npwvec,nomega,nspinor,&
& nfftot,ngfft,igfft0,Gsph_FFT,Ltg_q,green_enhigh_w,wfwfg,chi0)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz,nfftot,nomega,npwe,npwvec,nspinor,symchi
 logical,intent(in) :: qzero
 type(gsphere_t),intent(in) :: Gsph_FFT
 type(littlegroup_t),intent(in) :: Ltg_q
!arrays
 integer,intent(in) :: igfft0(npwvec),ngfft(18)
 complex(dpc),intent(in) :: green_enhigh_w(nomega)
 complex(gwpc),intent(in) :: wfwfg(nfftot*nspinor**2)
 complex(gwpc),intent(inout) :: chi0(npwe,npwe,nomega)

!Local variables ------------------------------
!scalars
 integer,save :: enough=0
 integer :: iSm1_g1mg2,iSm1_g1mg2_fft,ig,gmg_sph,gmg_fft
 integer :: igp,igstart,isym,itim,outofbox_wfn
 complex(gwpc) :: phmGt
 !character(len=500) :: msg

!************************************************************************

 igstart=1; if (qzero) igstart=2
 outofbox_wfn=0

 SELECT CASE (symchi)

 CASE (0) ! Do not use symmetries.
   ! MG: One has to make sure G1-G2 is still in the FFT mesh for each G1 and G2 in chi0 (not always true)
   ! MODULO wraps G1-G2 in the FFT box but the Fourier components are not periodic!
   do igp=igstart,npwe
     do ig=igstart,npwe
       gmg_fft = gsph_gmg_fftidx(Gsph_FFT,ig,igp,ngfft)
       if (gmg_fft==0) then
         outofbox_wfn=outofbox_wfn+1; CYCLE
       end if
       chi0(ig,igp,:) = chi0(ig,igp,:) + wfwfg(gmg_fft)*green_enhigh_w(:)
     end do
   end do

 CASE (1)
   ! Symmetrize the integrand in the full BZ.
   ! * <Sk b|e^{-i(G1-G2}.r}|b Sk> = e^{-i(G1-G2).\tau} <k b|e^{-i(S^{-1}(G1-G2).r)|b k>
   ! * green_enhigh_w in invariant under symmetry
   ! * We symmetrize using the operations of the little group of q since this routine
   !   is called inside a sum over IBZ_q, it would be possible to symmetrize
   !   this term by just summing over the IBZ and rotating the matrix elements.
   ! * Time-reversal does not lead to a complex conjugated since bra and ket are the same.
   !
   do igp=igstart,npwe
     do ig=igstart,npwe

      ! Get the index of G1-G2.
      gmg_sph = gsph_gmg_idx(Gsph_FFT,ig,igp)
      if (gmg_sph==0) then
        outofbox_wfn=outofbox_wfn+1; CYCLE
      end if

      do itim=1,Ltg_q%timrev
        do isym=1,Ltg_q%nsym_sg
          if (Ltg_q%wtksym(itim,isym,ik_bz)==1) then
            ! * This operation belongs to the little group and has to be used to reconstruct the BZ.
            ! * Time-reversal in not used to rotate (G1-G2) see comment above.
            phmGt          = Gsph_FFT%phmGt  (gmg_sph,isym)
            iSm1_g1mg2     = Gsph_FFT%rottbm1(gmg_sph,1,isym)
            iSm1_g1mg2_fft = igfft0(iSm1_g1mg2)

            chi0(ig,igp,:) = chi0(ig,igp,:) + phmGt*wfwfg(iSm1_g1mg2_fft)*green_enhigh_w(:)
          end if
        end do !isym
      end do !itim

     end do !igp
   end do !ig

 CASE DEFAULT
   MSG_BUG("Wrong value of symchi")
 END SELECT

 if (outofbox_wfn/=0) then
   enough=enough+1
   if (enough<=50) then
     MSG_WARNING(sjoin(' Number of G1-G2 pairs outside the G-sphere for Wfns: ', itoa(outofbox_wfn)))
     if (enough==50) then
       call wrtout(std_out,' ========== Stop writing Warnings ==========','COLL')
     end if
   end if
 end if

end subroutine completechi0_deltapart
!!***

!----------------------------------------------------------------------

!!****f* m_chi0tk/output_chi0sumrule
!! NAME
!! output_chi0sumrule
!!
!! FUNCTION
!!  Calculate and output the value of the sum rule for
!!  the non-interacting polarizability chi0
!!
!! INPUTS
!!
!! OUTPUT
!!  (for writing routines, no output)
!!  otherwise, should be described
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine output_chi0sumrule(qeq0,iq,npwe,omegaplasma,chi0sumrule,epsm1_w0,vc_sqrt)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq,npwe
 real(dp),intent(in) :: omegaplasma
 logical,intent(in) :: qeq0
!arrays
 real(dp),intent(inout) :: chi0sumrule(npwe)
 complex(gwpc),intent(in) :: epsm1_w0(npwe,npwe),vc_sqrt(npwe)

!Local variables ------------------------------
!scalars
 integer :: ig,igstart
 real(dp) :: average,norm
 character(len=500) :: msg

!************************************************************************

 igstart=1; if (qeq0) igstart=2
 !
 ! The sumrule reads:
 ! $ \int d\omega \omega v * Im[ \chi_0(\omega) ] = \pi/2 * w_p^2 $.
 chi0sumrule(igstart:npwe) = chi0sumrule(igstart:npwe) * vc_sqrt(igstart:npwe)**2
 !
 ! Calculate a weighted average of the fulfilment of the sumrule on epsilon
 ! The weight is given according to the significance of each q+G in the
 ! subsequent GW calculation: It is proportional to v * (epsm1 -1 )
 average = zero; norm = zero
 do ig=igstart,npwe
   average = average + chi0sumrule(ig) * real( vc_sqrt(ig)**2 * (epsm1_w0(ig,ig) - 1.0_dp ) )
   norm    = norm    +                   real( vc_sqrt(ig)**2 * (epsm1_w0(ig,ig) - 1.0_dp ) )
   !average = average + chi0sumrule(ig) * real(  (epsm1_w0(ig,ig) - 1.0_dp ) )
   !norm    = norm    +                   real(  (epsm1_w0(ig,ig) - 1.0_dp ) )
   !write(203,'(i4,8(2x,e12.6))') ig,1.0_dp/vc_sqrt(ig),chi0sumrule(ig)/ (0.5d0*omegaplasma**2*pi)
 end do

 if (abs(norm)>tol8) then
   write(msg,'(1x,a,i4,a,f10.2,2x,a)')&
    ' Average fulfillment of the sum rule on Im[epsilon] for q-point ',&
    iq,' :',average/norm/(0.5_dp*omegaplasma**2*pi)*100.0_dp,'[%]'
   call wrtout(std_out,msg,'COLL'); call wrtout(ab_out, msg,'COLL')
 end if

end subroutine output_chi0sumrule
!!***

!----------------------------------------------------------------------

!!****f* m_chi0tk/accumulate_chi0sumrule
!! NAME
!! accumulate_chi0sumrule
!!
!! FUNCTION
!!  Accumulate the contribution to the sum rule for Im chi0
!!  arising from a single transition. Eventually symmetrize
!!  it using the symmetry operations of the little group of q.
!!
!! INPUTS
!!  ik_bz=Index of the k-point in the full BZ whose contribution has to be added and symmetrized.
!!  symchi=1 if we are summing over IBZ_q and symmetrization has to be performed.
!!  npwe=Number of G vectors in chi0.
!!  npwepG0=Number of G vectors in the "enlarged" sphere to treat umklapp.
!!  factor=factor entering the expression.
!!  delta_ene=Transition energy.
!!  Ltg_q=<littlegroup_t>= Structure gathering information on the little group of the external q.
!!  Gsph_epsG0=<gsphere_t>=Info on the G-sphere for chi0.
!!  rhotwg(npwepG0)=Fouriet transform of u_{b1 k-q} u_{b2 k} in the "enlarged" sphere.
!!
!! OUTPUT
!!  See SIDES EFFECTS
!!
!! SIDES EFFECTS
!!  chi0sumrule(npwe)= In input the sum rule calculated so far,
!!  In output the contribution of this transition is accounted for, and, eventually, symmetrized.
!!  using the symmetry operations of the little group of the external q.
!!
!! PARENTS
!!      cchi0,cchi0q0
!!
!! CHILDREN
!!
!! SOURCE

subroutine accumulate_chi0sumrule(ik_bz,symchi,npwe,factor,delta_ene,&
& Ltg_q,Gsph_epsG0,npwepG0,rhotwg,chi0sumrule)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz,npwe,npwepG0,symchi
 real(dp),intent(in) :: delta_ene,factor
 type(gsphere_t),intent(in) :: Gsph_epsG0
 type(littlegroup_t),target,intent(in) :: Ltg_q
!arrays
 real(dp),intent(inout) :: chi0sumrule(npwe)
 complex(gwpc),intent(in) :: rhotwg(npwepG0)

!Local variables-------------------------------
!scalars
 integer :: isym,itim
 !character(len=500) :: msg
!arrays
 integer,allocatable :: Sm1_gmG0(:)
 integer, ABI_CONTIGUOUS pointer :: gmG0(:)
 complex(gwpc),allocatable :: rhotwg_sym(:)

!************************************************************************

 ! Accumulating the sum rule on chi0.
 ! Eq.(5.284) in G. D. Mahan Many-Particle Physics 3rd edition [[cite:Mahan2000]]

 SELECT CASE (symchi)
 CASE (0)
   ! Do not use symmetries, sum is performed in the full BZ.
   chi0sumrule(:)=chi0sumrule(:) + factor*delta_ene*ABS(rhotwg(1:npwe))**2

 CASE (1)
   ! Symmetrize the contribution in the full BZ.
   ABI_ALLOCATE(rhotwg_sym,(npwe))
   ABI_ALLOCATE(Sm1_gmG0,(npwe))

   do itim=1,Ltg_q%timrev
     do isym=1,Ltg_q%nsym_sg
       if (Ltg_q%wtksym(itim,isym,ik_bz)==1) then
        ! This operation belongs to the little group and has to be used to reconstruct the BZ ===
        ! In the following 2 lines mind the slicing (1:npwe)
        gmG0  => Ltg_q%igmG0(1:npwe,itim,isym)
        Sm1_gmG0(1:npwe)=Gsph_epsG0%rottbm1(gmG0(1:npwe),itim,isym)
        rhotwg_sym(1:npwe)=rhotwg(Sm1_gmG0)

        chi0sumrule(:)=chi0sumrule(:) + factor*delta_ene*ABS(rhotwg_sym(1:npwe))**2
       end if
     end do !isym
   end do !itim

   ABI_DEALLOCATE(rhotwg_sym)
   ABI_DEALLOCATE(Sm1_gmG0)

 CASE DEFAULT
   MSG_BUG(sjoin('Wrong value for symchi:', itoa(symchi)))
 END SELECT

end subroutine accumulate_chi0sumrule
!!***

!!****f* m_chi0tk/make_transitions
!! NAME
!! make_transitions
!!
!! FUNCTION
!!  Calculate transition energies entering the espression for the irreducible polarizability.
!!
!! INPUTS
!!  nsspol=1 for spin unpolarized, 2 for spin polarized calculations
!!  nbnds=total number of bands
!!  kmesh<kmesh_t>=datatype gathering info on the k-mesh:
!!   | %nbz=number of k-points in the full BZ
!!   | %nibz=number of k-points in the IBZ
!!   | %tab(nkbz)=table giving for each k-point in the BZ, the corresponding irreducible point in the IBZ array
!!   | %bz(3,nkbz)=reduced coordinated of k-points
!!  TOL_DELTA_OCC=tolerance on the difference of the occupation numbers
!!  gw_energy(nbnds,kmesh%nkibz,nsppol)=quasi-particle energies energies
!!  occ(nbnds,kmesh%nkibz,nsppol)=occupation numbers
!!  chi0alg=integer defining the method used to calculate chi0
!!   0 ==> calculate chi0 using the Adler-Wiser expression
!!   1 ==> use spectral method
!!  timrev=if 2, time-reversal symmetry is considered; 1 otherwise
!!
!! OUTPUT
!! my_max_rest,my_min_rest=Maximum and minimum resonant (posite) transition energy.
!! max_rest,min_rest=Maximun and minimum resonant (posite) transition energy treated by this node.
!!
!! PARENTS
!!      cchi0,cchi0q0
!!
!! CHILDREN
!!
!! SOURCE

subroutine make_transitions(Wfd,chi0alg,nbnds,nbvw,nsppol,symchi,timrev,TOL_DELTA_OCC,&
& max_rest,min_rest,my_max_rest,my_min_rest,Kmesh,Ltg_q,gw_energy,occ,qpoint,bbp_ks_distrb)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: chi0alg,nbnds,nbvw,nsppol,symchi,timrev
 real(dp),intent(in) :: TOL_DELTA_OCC
 real(dp),intent(out) :: max_rest,min_rest
 real(dp),intent(out) :: my_max_rest,my_min_rest
 type(kmesh_t),intent(in) :: Kmesh
 type(littlegroup_t),intent(in) :: Ltg_q
 type(wfd_t),intent(in) :: Wfd
!arrays
 real(dp),intent(in) :: gw_energy(nbnds,Kmesh%nibz,nsppol)
 real(dp),intent(in) :: occ(nbnds,Kmesh%nibz,nsppol),qpoint(3)
 integer,intent(in) :: bbp_ks_distrb(Wfd%mband,Wfd%mband,Kmesh%nbz,Wfd%nsppol)

!Local variables-------------------------------
!scalars
 integer :: ib1,ib2,ii,ik_bz,ik_ibz,ikmq_bz,ikmq_ibz,is,nt,ntrans,my_ntrans,iloop
 real(dp) :: delta_ene,delta_occ,spin_fact
 character(len=500) :: msg
!arrays
 integer :: G0(3)
 real(dp) :: kmq(3)

!************************************************************************

 DBG_ENTER("COLL")

 if (chi0alg<0 .or. chi0alg>=2) then
   MSG_BUG(sjoin('chi0alg:', itoa(chi0alg),' not allowed'))
 end if
 if (timrev/=1 .and. timrev/=2) then
   MSG_BUG(sjoin('timrev:', itoa(timrev),' not allowed'))
 end if

 ABI_UNUSED(nbvw)
 !
 ! In the first loop calculate total number of transitions for this q-point
 ! as well min and max transition without taking into account distribution of bands.
 ! In the second iteration calculate min and Max transition for this processor.
 !
 spin_fact=half; if (nsppol==2) spin_fact=one
 my_max_rest=smallest_real; my_min_rest=greatest_real
    max_rest=smallest_real;    min_rest=greatest_real

 do iloop=1,2
   nt=0
   do ik_bz=1,Kmesh%nbz
     ik_ibz=Kmesh%tab(ik_bz)
     kmq(:)=Kmesh%bz(:,ik_bz)-qpoint(:)

     if (symchi==1) then
       if (Ltg_q%ibzq(ik_bz)/=1) cycle ! This point does not belong to the IBZ defined by the little group
     end if

     ! Find kp=k-q-G0 and also G0 where kp is in the first BZ
     if (.not.has_BZ_item(Kmesh,kmq,ikmq_bz,g0)) then ! Stop as the weight 1.0/nkbz is wrong.
       write(msg,'(4a,2(2a,3f12.6),2a)')ch10,&
&        ' make_transitions : ERROR - ',ch10,&
&        ' kp  = k-q-G0 not found in the BZ mesh',ch10,&
&        ' k   = ',(Kmesh%bz(ii,ik_bz),ii=1,3),ch10,&
&        ' k-q = ',(kmq(ii),ii=1,3),ch10,&
&        ' weight in cchi0/cchi0q is wrong '
       MSG_ERROR(msg)
     end if

     ikmq_ibz=Kmesh%tab(ikmq_bz)
     do is=1,nsppol
       do ib1=1,nbnds
         do ib2=1,nbnds

           if (iloop==2) then
             if (bbp_ks_distrb(ib1,ib2,ik_bz,is)/=Wfd%my_rank) cycle
           end if

           if (timrev==2 .and. ib1<ib2) cycle ! Thanks to time-reversal we gain a factor ~2.

           delta_occ=spin_fact*(occ(ib1,ikmq_ibz,is)-occ(ib2,ik_ibz,is))
           delta_ene=gw_energy(ib1,ikmq_ibz,is)-gw_energy(ib2,ik_ibz,is)

           if (chi0alg==0)  then
             ! Adler-Wiser expression. Skip only if factor due to occupation number is smaller than TOL_DELTA_OCC
             if (abs(delta_occ) < abs(TOL_DELTA_OCC)) cycle
           else if (chi0alg==1) then
             ! Spectral method with time-reversal, only resonant transitions
             ! This has to changed to include spectral method without time-reversal
             if (delta_ene < -abs(TOL_DELTA_OCC) .or. abs(delta_occ) < abs(TOL_DELTA_OCC)) cycle
           end if

           ! We have a new transition
           nt=nt+1

           if (iloop==1) then
             max_rest=MAX(max_rest,zero,delta_ene)
             if (delta_ene>=-tol6) min_rest=MIN(min_rest,delta_ene)
           end if
           if (iloop==2) then
             my_max_rest=MAX(my_max_rest,zero,delta_ene)
             if (delta_ene>=-tol6) my_min_rest=MIN(my_min_rest,delta_ene)
           end if

         end do
       end do
     end do
   end do
   if (iloop==1) ntrans=nt
   if (iloop==2) my_ntrans=nt
 end do !iloop

 write(msg,'(2a,i9,2a,f8.3,3a,f8.3,a)')ch10,&
&  ' Total number of transitions = ',ntrans,ch10,&
&  ' min resonant     = ',min_rest*Ha_eV,' [eV] ',ch10,&
&  ' Max resonant     = ',max_rest*Ha_eV,' [eV] '
 call wrtout(std_out,msg,'COLL')

 if (Wfd%nproc/=1) then
   write(msg,'(2a,i9,2a,f8.3,3a,f8.3,a)')ch10,&
&    ' Total number of transitions for this processor= ',my_ntrans,ch10,&
&    ' min resonant     = ',my_min_rest*Ha_eV,' [eV] ',ch10,&
&    ' Max resonant     = ',my_max_rest*Ha_eV,' [eV] '
   call wrtout(std_out,msg,'PERS')
 end if

 DBG_EXIT("COLL")

end subroutine make_transitions
!!***

!----------------------------------------------------------------------

!!****f* m_chi0tk/chi0_bbp_mask
!! NAME
!!  chi0_bbp_mask
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      cchi0,cchi0q0
!!
!! CHILDREN
!!
!! SOURCE

subroutine chi0_bbp_mask(Ep,use_tr,QP_BSt,mband,ikmq_ibz,ik_ibz,spin,spin_fact,bbp_mask)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spin,ik_ibz,ikmq_ibz,mband
 real(dp),intent(in) :: spin_fact
 logical,intent(in) :: use_tr
 type(em1params_t),intent(in) :: Ep
 type(ebands_t),target,intent(in) :: QP_BSt
!arrays
 logical,intent(out) :: bbp_mask(mband,mband)

!Local variables-------------------------------
!scalars
 integer :: ib1,ib2
 real(dp) :: deltaeGW_b1kmq_b2k,deltaf_b1kmq_b2k,e_b1_kmq,f_b1_kmq
!arrays
 real(dp), ABI_CONTIGUOUS pointer :: qp_energy(:,:,:),qp_occ(:,:,:)

!************************************************************************

 qp_energy => QP_BSt%eig; qp_occ => QP_BSt%occ
 bbp_mask=.FALSE.

 !use_tr = (Ep%awtr==1)
 SELECT CASE (Ep%gwcomp)
 CASE (0)
   do ib1=1,Ep%nbnds
     ! Loop over "conduction" states.
     e_b1_kmq = qp_energy(ib1,ikmq_ibz,spin)
     f_b1_kmq =    qp_occ(ib1,ikmq_ibz,spin)

     do ib2=1,Ep%nbnds ! Loop over "valence" states.
       deltaf_b1kmq_b2k   = spin_fact*(f_b1_kmq-qp_occ(ib2,ik_ibz,spin))
       deltaeGW_b1kmq_b2k = e_b1_kmq-qp_energy(ib2,ik_ibz,spin)

       SELECT CASE (Ep%spmeth)
       CASE (0)
         ! Standard Adler-Wiser expression.
         if (ABS(deltaf_b1kmq_b2k) >= GW_TOL_DOCC) then
           bbp_mask(ib1,ib2)=.TRUE.
           if (use_tr .and. ib1<ib2) bbp_mask(ib1,ib2)=.FALSE. ! GAIN a factor ~2 thanks to time-reversal.
         end if

       CASE (1,2)
         ! Spectral method, WARNING time-reversal here is always assumed!
         if (ABS(deltaf_b1kmq_b2k) >= GW_TOL_DOCC) then
           bbp_mask(ib1,ib2)=.TRUE.
           if (deltaeGW_b1kmq_b2k<zero) bbp_mask(ib1,ib2)=.FALSE. ! Only positive frequencies are needed for the Hilbert transform.
           !$if (use_tr .and. ib1<ib2) bbp_mask(ib1,ib2)=.FALSE. ! GAIN a factor ~2 thanks to time-reversal.
         end if

       CASE DEFAULT
         MSG_ERROR(sjoin(" Wrong value for spmeth:", itoa(Ep%spmeth)))
       END SELECT
       !write(std_out,*) "bbp_mask(ib1,ib2)",bbp_mask(ib1,ib2)
     end do !ib2
   end do !ib1

 CASE (1)
   ! Extrapolar technique
   ABI_CHECK(Ep%spmeth==0,"Hilbert transform and extrapolar method are not compatible")

   ! Loop over "conduction" states.
   do ib1=1,Ep%nbnds
     e_b1_kmq=qp_energy(ib1,ikmq_ibz,spin)
     f_b1_kmq=   qp_occ(ib1,ikmq_ibz,spin)

     ! Loop over "valence" states.
     do ib2=1,Ep%nbnds
       deltaf_b1kmq_b2k  =spin_fact*(f_b1_kmq-qp_occ(ib2,ik_ibz,spin))
       deltaeGW_b1kmq_b2k=e_b1_kmq-qp_energy(ib2,ik_ibz,spin)

       ! When the completeness correction is used,
       ! we need to also consider transitions with vanishing deltaf
       ! Rangel: This is to compute chi in metals correctly with the extrapolar method.
       bbp_mask(ib1,ib2)=.TRUE.
       !if (qp_occ(ib2,ik_ibz,is) < GW_TOL_DOCC) CYCLE
       if (qp_occ(ib2,ik_ibz,spin) < GW_TOL_DOCC .and. (ABS(deltaf_b1kmq_b2k) < GW_TOL_DOCC .or. ib1<ib2)) then
         bbp_mask(ib1,ib2)=.FALSE.
       end if

     end do
   end do

  CASE DEFAULT
    MSG_ERROR(sjoin("Wrong value of gwcomp:", itoa(Ep%gwcomp)))
  END SELECT

end subroutine chi0_bbp_mask
!!***

END MODULE m_chi0tk
!!***
