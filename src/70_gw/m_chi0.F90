!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_chi0
!! NAME
!!  m_chi0
!!
!! FUNCTION
!!  This module provides tools for the computation of the irreducible polarizability.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (MG)
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

MODULE m_chi0

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xomp
 use m_sort

 use m_gwdefs,   only : GW_TOL_DOCC, czero_gw, cone_gw, em1params_t, j_gw
 use m_fstrings, only : sjoin, itoa
 use m_blas,     only : xgerc, xgemm
 use m_crystal,  only : crystal_t
 use m_gsphere,  only : gsphere_t
 use m_bz_mesh,  only : littlegroup_t

 implicit none

 private

 public :: assemblychi0q0_sym
 public :: assemblychi0_sym
 !public :: mkrhotwg_sigma
 public :: assemblychi0sfq0
 public :: symmetrize_afm_chi0
 public :: accumulate_chi0_q0
 public :: accumulate_sfchi0_q0
 public :: assemblychi0sf
 public :: approxdelta
 public :: setup_spectral
 public :: hilbert_transform
 public :: hilbert_transform_headwings
!!***

CONTAINS  !=======================================================================================================
!!***

!!****f* m_chi0/assemblychi0q0_sym
!! NAME
!! assemblychi0q0_sym
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
!!  nqlwl=Number of small q-points used to deal with the optical limit.
!!  qlwl(3,nqlwl)=reciprocal space coordinates of the q-point for long-wavelength limit treatment.
!!  rhotwg(npwepG0*nspinor**2)=Oscillator matrix elements corresponding to an occupied-unoccupied pair of states.
!!  rhotwx(3,nspinor**2)=Matrix element of the operator -i[H,r]/(e1-e2) in reciprocal lattice units.
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
!!  lwing(Ep%npwe*Ep%nI,Ep%nomega,3)=Lower wing (calculated only if nqlwl > 1 )
!!  uwing(Ep%npwe*Ep%nJ,Ep%nomega,3)=Upper wing (calculated only if nqlwl > 1 )
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
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine assemblychi0q0_sym(nqlwl,qlwl,ik_bz,isym_kbz,itim_kbz,gwcomp,nspinor,npwepG0,Ep,Cryst,Ltg_q,Gsph_epsG0,&
& chi0,rhotwx,rhotwg,green_w,green_enhigh_w,deltaf_b1b2,lwing,uwing)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'assemblychi0q0_sym'
 use interfaces_32_util
 use interfaces_70_gw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz,isym_kbz,itim_kbz
 integer,intent(in) :: npwepG0,nqlwl,nspinor,gwcomp
 real(dp),intent(in) :: deltaf_b1b2
 type(littlegroup_t),intent(in) :: Ltg_q
 type(gsphere_t),target,intent(in) :: Gsph_epsG0
 type(crystal_t),intent(in) :: Cryst
 type(em1params_t),intent(in) :: Ep
!arrays
 real(dp),intent(in) :: qlwl(3,nqlwl)
 complex(gwpc),intent(inout) :: rhotwg(npwepG0*nspinor**2)
 complex(gwpc),intent(in) :: rhotwx(3,nspinor**2)
 complex(gwpc),intent(inout) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)
 complex(dpc),intent(in) :: green_w(Ep%nomega),green_enhigh_w(Ep%nomega)
 complex(dpc),intent(inout) :: lwing(Ep%npwe*Ep%nI,Ep%nomega,3)
 complex(dpc),intent(inout) :: uwing(Ep%npwe*Ep%nJ,Ep%nomega,3)

!Local variables-------------------------------
!scalars
 integer :: itim,io,isym,igp,ig,iqlwl,nthreads
 integer :: jj,ii,s_jj,pad_jj,pad_ii,ig1,ig2
 complex(gwpc) :: dd,mqg0,mqg0_sym,rhotwg0_bkp
 !character(len=500) :: msg
!arrays
 integer,ABI_CONTIGUOUS pointer :: Sm1G(:)
 real(dp) :: opinv(3,3),qrot(3)
 real(dp) :: b1(3),b2(3),b3(3)
 complex(gwpc),allocatable :: rhotwg_sym(:),rhotwg_I(:),rhotwg_J(:)
 complex(gwpc),allocatable :: rhotwg_sym_star(:),rhotwg_star(:)
 complex(gwpc), ABI_CONTIGUOUS pointer :: phmGt(:)

!************************************************************************

 b1=two_pi*Gsph_epsG0%gprimd(:,1)
 b2=two_pi*Gsph_epsG0%gprimd(:,2)
 b3=two_pi*Gsph_epsG0%gprimd(:,3)

 nthreads = xomp_get_max_threads()

 SELECT CASE (Ep%symchi)

 CASE (0) ! Do not use symmetries.

   if (nspinor==1) then
    ! * Accumulate over the full BZ i.e
    !    chi0(G1,G2,io) = chi0(G1,G2,io) + (rhotwg(G1)*CONJG(rhotwg(G2)))*green_w(io)
    ! * The non-analytic term is symmetrized for this k-point in the BZ according to:
    !    rhotwg(1)= S^-1q * rhotwx_ibz
    !    rhotwg(1)=-S^-1q * CONJG(rhotwx_ibz) if time-reversal is used.
    opinv(:,:)=REAL(Cryst%symrec(:,:,isym_kbz),dp)
    call matrginv(opinv,3,3)
    qrot = (3-2*itim_kbz) * MATMUL(opinv,qlwl(:,1))

    rhotwg(1)=dotproductqrc(qrot,rhotwx(:,1),b1,b2,b3) !TODO get rid of this
    if (itim_kbz==2) rhotwg(1) = GWPC_CONJG(rhotwg(1))

    if (gwcomp==1) then ! Leave the head and wings uncorrected (does not matter much)
      if (ABS(deltaf_b1b2) < GW_TOL_DOCC) rhotwg(1)=czero_gw
      do igp=1,Ep%npwe
        chi0(1,igp,:) = chi0(1,igp,:) + rhotwg(1) * GWPC_CONJG(rhotwg(igp))*green_enhigh_w(:)
      end do
      do ig=2,Ep%npwe
        chi0(ig,1,:)  = chi0(ig,1,:)  + rhotwg(ig) * GWPC_CONJG(rhotwg(1))  *green_enhigh_w(:)
      end do
    end if

    ! Multiply elements G1,G2 of rhotwg_sym by green_w(io) and accumulate in chi0(G1,G2,io)
     if (nthreads==1 .or. Ep%nomega>=nthreads) then
!$omp parallel do private(dd)
       do io=1,Ep%nomega
         dd=green_w(io)
         call XGERC(Ep%npwe,Ep%npwe,dd,rhotwg,1,rhotwg,1,chi0(:,:,io),Ep%npwe)
       end do
    else
!$omp parallel private(dd)
       do io=1,Ep%nomega
         dd=green_w(io)
!$omp do
         do ig2=1,Ep%npwe
           do ig1=1,Ep%npwe
             chi0(ig1,ig2,io) = chi0(ig1,ig2,io) + dd * rhotwg(ig1) * GWPC_CONJG(rhotwg(ig2))
           end do
         end do
!$omp end do NOWAIT
       end do
!$omp end parallel
    end if

    ! === Accumulate heads and wings for each small q ===
    ! * For better performance, this part is not done if nqlwl==1
    !   lwing and uwing will be filled in cchi0q0 after the MPI collective sum
    ! FIXME extrapolar method should be checked!!
    !do io=1,Ep%nomega
    ! lwing(:,io,1) =  chi0(:,1,io)
    ! uwing(:,io,1) =  chi0(1,:,io)
    !end do

    if (nqlwl>1.and..FALSE.) then
      rhotwg0_bkp = rhotwg(1) ! Save G=0 value of the first q
      ABI_MALLOC(rhotwg_star,(Ep%npwe))
      rhotwg_star = GWPC_CONJG(rhotwg(1:Ep%npwe))

      do iqlwl=2,nqlwl
        qrot =  (3-2*itim_kbz) * MATMUL(opinv,qlwl(:,iqlwl))
        mqg0 = dotproductqrc(qrot,rhotwx(:,1),b1,b2,b3) !TODO get rid of this
        if (itim_kbz==2) mqg0 = GWPC_CONJG(mqg0)
        rhotwg     (1) = mqg0
        rhotwg_star(1) = GWPC_CONJG(mqg0)
        !
        ! here we might take advantage of Hermiticity along Im axis in RPA (see mkG0w)
        do io=1,Ep%nomega
          lwing(:,io,iqlwl) = lwing(:,io,iqlwl) + rhotwg     (1:Ep%npwe) * GWPC_CONJG(mqg0) * green_w(io)
          uwing(:,io,iqlwl) = uwing(:,io,iqlwl) + rhotwg_star(1:Ep%npwe) *            mqg0  * green_w(io)
        end do
      end do ! iqlwl

      ABI_FREE(rhotwg_star)
      rhotwg(1) = rhotwg0_bkp ! Reinstate previous value of rhotwg(1).
    end if ! nqlwl > 1

  else ! spinorial case
    ABI_MALLOC(rhotwg_I,(Ep%npwe))
    ABI_MALLOC(rhotwg_J,(Ep%npwe))

    ABI_CHECK(nqlwl==1,"nqlwl/=1 Not implemented")

    ! I can use symmetries to loop over the upper triangle but
    ! this makes using BLAS more difficult
    ! Important NOTE: treatment of q-->0 limit is correct only
    ! for i=j=0. Other components require additional terms.

    do jj=1,Ep%nJ
      s_jj=1 ; if (jj==4) s_jj=-1
      pad_jj=(jj-1)*Ep%npwe
      call mkrhotwg_sigma(jj,nspinor,Ep%npwe,rhotwg,rhotwg_J)

      rhotwg_J(1) = q0limit(jj,qlwl(:,1),nspinor,rhotwx,b1,b2,b3)
      !TODO RECHECK this
      if (itim_kbz==2) rhotwg_J(1)=-GWPC_CONJG(rhotwg_J(1))

      do ii=1,Ep%nI
        pad_ii=(ii-1)*Ep%npwe

        if (ii/=jj) then
          call mkrhotwg_sigma(ii,nspinor,Ep%npwe,rhotwg,rhotwg_I)
          rhotwg_I(1) = q0limit(ii,qlwl(:,1),nspinor,rhotwx,b1,b2,b3)
          if (itim_kbz==2) rhotwg_I(1)=-GWPC_CONJG(rhotwg_I(1))
        else
          rhotwg_I(:)=rhotwg_J(:)
        end if

        do io=1,Ep%nomega
          dd = s_jj*green_w(io)
          call XGERC(Ep%npwe,Ep%npwe,dd,rhotwg_I,1,rhotwg_J,1,&
&           chi0(pad_ii+1:pad_ii+Ep%npwe,pad_jj+1:pad_jj+Ep%npwe,io),Ep%npwe)
        end do

      end do !ii
    end do !jj

    ABI_FREE(rhotwg_I)
    ABI_FREE(rhotwg_J)
  end if

 CASE (1) ! Use symmetries to reconstruct the integrand.
   if (nspinor==1) then
     ABI_MALLOC(rhotwg_sym,(Ep%npwe))

     ! === Loop over symmetries of the space group and time-reversal ===
     do isym=1,Ltg_q%nsym_sg
       do itim=1,Ltg_q%timrev

         if (Ltg_q%wtksym(itim,isym,ik_bz)==1) then
           ! === This operation belongs to the little group and has to be considered to reconstruct the BZ ===
           ! TODO this is a hot-spot, should add a test on the umklapp
           phmGt => Gsph_epsG0%phmGt  (1:Ep%npwe,isym) ! In the 2 lines below note the slicing (1:npwe)
           Sm1G  => Gsph_epsG0%rottbm1(1:Ep%npwe,itim,isym)

           opinv(:,:)=REAL(Cryst%symrec(:,:,isym),dp)
           call matrginv(opinv,3,3)
           qrot = (3-2*itim) * MATMUL(opinv,qlwl(:,1))

           SELECT CASE (itim)

           CASE (1)
             rhotwg_sym(1:Ep%npwe)=rhotwg(Sm1G(1:Ep%npwe))*phmGt(1:Ep%npwe)
             rhotwg_sym(1)=dotproductqrc(qrot,rhotwx(:,1),b1,b2,b3)

           CASE (2)
             rhotwg_sym(1:Ep%npwe)=GWPC_CONJG(rhotwg(Sm1G(1:Ep%npwe)))*phmGt(1:Ep%npwe)
             rhotwg_sym(1)=GWPC_CONJG(dotproductqrc(qrot,rhotwx(:,1),b1,b2,b3))

           CASE DEFAULT
             MSG_BUG(sjoin('Wrong value of itim:', itoa(itim)))
           END SELECT

           if (gwcomp==1) then ! Leave the head and wings uncorrected (does not matter much)
             if (ABS(deltaf_b1b2) < GW_TOL_DOCC) rhotwg_sym(1)=czero_gw
             do igp=1,Ep%npwe
               chi0(1,igp,:) = chi0(1,igp,:) + rhotwg_sym(1) *GWPC_CONJG(rhotwg_sym(igp))*green_enhigh_w(:)
             end do
             do ig=2,Ep%npwe
               chi0(ig,1,:)  = chi0(ig,1,:)  + rhotwg_sym(ig)*GWPC_CONJG(rhotwg_sym(1))  *green_enhigh_w(:)
             end do
           end if

           ! Multiply elements G1,G2 of rhotwg_sym by green_w(io) and accumulate in chi0(G,Gp,io)
           if (nthreads==1 .or. Ep%nomega>=nthreads) then
! MKL_10.3 xgerc is not threaded.
!$omp parallel do private(dd)
             do io=1,Ep%nomega
               dd=green_w(io)
               call XGERC(Ep%npwe,Ep%npwe,dd,rhotwg_sym,1,rhotwg_sym,1,chi0(:,:,io),Ep%npwe)
             end do
           else
!$omp parallel private(dd)
             do io=1,Ep%nomega
               dd=green_w(io)
!$omp do
               do ig2=1,Ep%npwe
                 do ig1=1,Ep%npwe
                   chi0(ig1,ig2,io) = chi0(ig1,ig2,io) + dd * rhotwg_sym(ig1) * GWPC_CONJG(rhotwg_sym(ig2))
                 end do
               end do
!$omp end do NOWAIT
             end do
!$omp end parallel
           end if
           !
           ! === Accumulate heads and wings for each small q ===
           ! * For better performance, this part is not done if nqlwl==1
           !   lwing and uwing will be filled in cchi0q0 after the MPI collective sum
           ! FIXME extrapolar method should be checked!!
           if (nqlwl>1.and..FALSE.) then
             ABI_MALLOC(rhotwg_sym_star,(Ep%npwe))
             rhotwg_sym_star = GWPC_CONJG(rhotwg_sym)

             do iqlwl=2,nqlwl
               qrot = (3-2*itim) * MATMUL(opinv,qlwl(:,iqlwl))
               mqg0_sym = dotproductqrc(qrot,rhotwx(:,1),b1,b2,b3)
               if (itim==2) mqg0_sym = GWPC_CONJG(mqg0_sym)

               rhotwg_sym     (1) =       mqg0_sym
               rhotwg_sym_star(1) = GWPC_CONJG(mqg0_sym)

               ! here we might take advantage of Hermiticity along Im axis in RPA (see mkG0w)
               do io=1,Ep%nomega
                 lwing(:,io,iqlwl) = lwing(:,io,iqlwl) + rhotwg_sym     (1:Ep%npwe) * GWPC_CONJG(mqg0_sym) * green_w(io)
                 uwing(:,io,iqlwl) = uwing(:,io,iqlwl) + rhotwg_sym_star(1:Ep%npwe) *       mqg0_sym  * green_w(io)
               end do
             end do !iqlwl

             ABI_FREE(rhotwg_sym_star)
           end if !nqlwl>1

         end if !wtksym
       end do !itim
     end do !isym

     ABI_FREE(rhotwg_sym)

   else  !spinorial case
     MSG_BUG('symchi=1 with spinor not implemented ')
     ABI_CHECK(nqlwl==1,"nqlwl/=1 Not implemented")
   end if

 CASE DEFAULT
   MSG_BUG(sjoin('Wrong value of symchi:', itoa(Ep%symchi)))
 END SELECT

end subroutine assemblychi0q0_sym
!!***

!----------------------------------------------------------------------

!!****f* m_chi0/assemblychi0_sym
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
!!  rhotwg(npwe*nspinor**2)=Oscillator matrix elements for this k-point and the transition that has to be summed
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
!!      wrtout
!!
!! SOURCE


subroutine assemblychi0_sym(ik_bz,nspinor,Ep,Ltg_q,green_w,npwepG0,rhotwg,Gsph_epsG0,chi0)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'assemblychi0_sym'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz,npwepG0,nspinor
 type(gsphere_t),intent(in) :: Gsph_epsG0
 type(littlegroup_t),intent(in) :: Ltg_q
 type(em1params_t),intent(in) :: Ep
!arrays
 complex(gwpc),intent(in) :: rhotwg(npwepG0*nspinor**2)
 complex(dpc),intent(in) :: green_w(Ep%nomega)
 complex(gwpc),intent(inout) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)

!Local variables-------------------------------
!scalars
 integer :: itim,io,isym,ig1,ig2,nthreads
 integer :: jj,ii,s_jj,pad_jj,pad_ii
 complex(gwpc) :: dd
 !character(len=500) :: msg
!arrays
 integer :: Sm1_gmG0(Ep%npwe)
 complex(gwpc) :: rhotwg_sym(Ep%npwe)
 complex(gwpc),allocatable :: rhotwg_I(:),rhotwg_J(:)

! *************************************************************************

 nthreads = xomp_get_max_threads()

 SELECT CASE (Ep%symchi)

 CASE (0) ! Do not use symmetries

   if (nspinor==1) then
     if (nthreads==1 .or. Ep%nomega>=nthreads) then
! MKL_10.3 xgerc is not threaded. BTW: single precision is faster (sometimes factor ~2).
!$omp parallel do private(dd)
       do io=1,Ep%nomega
         dd=green_w(io)
         call XGERC(Ep%npwe,Ep%npwe,dd,rhotwg,1,rhotwg,1,chi0(:,:,io),Ep%npwe)
       end do
     else
!$omp parallel private(dd)
       do io=1,Ep%nomega
         dd=green_w(io)
!$omp do
         do ig2=1,Ep%npwe
           do ig1=1,Ep%npwe
             chi0(ig1,ig2,io) = chi0(ig1,ig2,io) + dd * rhotwg(ig1) * GWPC_CONJG(rhotwg(ig2))
           end do
         end do
!$omp end do NOWAIT
       end do
!$omp end parallel
    end if

   else ! spinorial case
     ABI_MALLOC(rhotwg_I,(Ep%npwe))
     ABI_MALLOC(rhotwg_J,(Ep%npwe))

     ! I can use symmetries to loop over the upper triangle but
     ! this makes using BLAS more difficult

     do jj=1,Ep%nJ
       s_jj=1 ; if (jj==4) s_jj=-1
       pad_jj=(jj-1)*Ep%npwe
       call mkrhotwg_sigma(jj,nspinor,Ep%npwe,rhotwg,rhotwg_J)
       !
       do ii=1,Ep%nI
         pad_ii=(ii-1)*Ep%npwe

         if (ii/=jj) then
          call mkrhotwg_sigma(ii,nspinor,Ep%npwe,rhotwg,rhotwg_I)
         else
          rhotwg_I(:)=rhotwg_J(:)
         end if

         do io=1,Ep%nomega
          dd = s_jj*green_w(io)
          call XGERC(Ep%npwe,Ep%npwe,dd,rhotwg_I,1,rhotwg_J,1,chi0(pad_ii+1:pad_ii+Ep%npwe,pad_jj+1:pad_jj+Ep%npwe,io),Ep%npwe)
         end do
         !
       end do !ii
     end do !jj

     ABI_FREE(rhotwg_I)
     ABI_FREE(rhotwg_J)
   end if

 CASE (1) ! Use symmetries to reconstruct the integrand in the BZ.
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
   !
   ! === Loop over symmetries of the space group and time-reversal ===
   do isym=1,Ltg_q%nsym_sg
     do itim=1,Ltg_q%timrev
       !
       if (Ltg_q%wtksym(itim,isym,ik_bz)==1) then
         ! === This operation belongs to the little group and has to be used to reconstruct the BZ ===
         ! * In the following 3 lines mind the slicing (1:npwe)
         ! TODO this is a hot-spot, should add a test on the umklapp
         !
         !gmG0 => Ltg_q%igmG0(1:Ep%npwe,itim,isym)
         Sm1_gmG0(1:Ep%npwe) = Gsph_epsG0%rottbm1( Ltg_q%igmG0(1:Ep%npwe,itim,isym), itim,isym)

         SELECT CASE (itim)
         CASE (1)
           rhotwg_sym(1:Ep%npwe) = rhotwg(Sm1_gmG0) * Gsph_epsG0%phmGt(1:Ep%npwe,isym)
         CASE (2)
           rhotwg_sym(1:Ep%npwe) = GWPC_CONJG(rhotwg(Sm1_gmG0))*Gsph_epsG0%phmGt(1:Ep%npwe,isym)
         CASE DEFAULT
           MSG_BUG(sjoin('Wrong itim:', itoa(itim)))
         END SELECT
         !
         ! Multiply rhotwg_sym by green_w(io) and accumulate in chi0(G,Gp,io)


if (nthreads==1 .or. Ep%nomega>=nthreads) then
! MKL_10.3 xgerc is not threaded. BTW: single precision is faster (sometimes factor ~2).
!$omp parallel do private(dd)
         do io=1,Ep%nomega
           dd=green_w(io)
           call XGERC(Ep%npwe,Ep%npwe,dd,rhotwg_sym,1,rhotwg_sym,1,chi0(:,:,io),Ep%npwe)
         end do

else
!write(std_out,*)"in new with nthreads = ",nthreads
!$omp parallel private(dd)
         do io=1,Ep%nomega
           dd=green_w(io)
!$omp do
           do ig2=1,Ep%npwe
             do ig1=1,Ep%npwe
               chi0(ig1,ig2,io) = chi0(ig1,ig2,io) + dd * rhotwg_sym(ig1) * GWPC_CONJG(rhotwg_sym(ig2))
             end do
           end do
!$omp end do NOWAIT
         end do
!$omp end parallel
         !
end if
       end if
       !
     end do
   end do

 CASE DEFAULT
   MSG_BUG(sjoin('Wrong symchi:', itoa(Ep%symchi)))
 END SELECT

end subroutine assemblychi0_sym
!!***

!----------------------------------------------------------------------

!!****f* m_chi0/mkrhotwg_sigma
!! NAME
!! mkrhotwg_sigma
!!
!! FUNCTION
!!  Helper function used to calculate selected linear combination
!!  of the oscillator matrix elements in the case of noncollinear magnetism.
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
!!      m_chi0
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine mkrhotwg_sigma(ii,nspinor,npw,rhotwg,rhotwg_I)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkrhotwg_sigma'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ii,npw,nspinor
!arrays
 complex(gwpc),intent(in) :: rhotwg(npw*nspinor**2)
 complex(gwpc),intent(out) :: rhotwg_I(npw)

! *************************************************************************

 SELECT CASE (ii)
 CASE (1) ! $ M_0 = M_{\up,\up} + M_{\down,\down} $
   rhotwg_I(:) = rhotwg(1:npw) + rhotwg(npw+1:2*npw)
 CASE (2) ! $ M_z = M_{\up,\up} - M_{\down,\down} $
   rhotwg_I(:) = rhotwg(1:npw) - rhotwg(npw+1:2*npw)
 CASE (3) ! $ M_x = M_{\up,\down} + M_{\down,\up} $
   rhotwg_I(:) = ( rhotwg(2*npw+1:3*npw) + rhotwg(3*npw+1:4*npw) )
 CASE (4) ! $ M_y = i * (M_{\up,\down} -M_{\down,\up}) $
   rhotwg_I(:) = (rhotwg(2*npw+1:3*npw) - rhotwg(3*npw+1:4*npw) )*j_gw
 CASE DEFAULT
   MSG_BUG(sjoin('Wrong ii value:', itoa(ii)))
 END SELECT

end subroutine mkrhotwg_sigma
!!***

!----------------------------------------------------------------------

!!****f* m_chi0/assemblychi0sfq0
!! NAME
!! assemblychi0sfq0
!!
!! FUNCTION
!! Update the spectral function of the independent particle susceptibility at q==0 for the contribution
!! of one pair of occupied-unoccupied band, for each frequency.
!! If symchi==1, the symmetries belonging to the little group of the external point q are used
!! to reconstrunct the contributions in the full Brillouin zone. In this case, the equation implented is:
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
!!  nqlwl=Number of q-points used for the optical limit.
!!  qlwl(3,nqlwl)=Reciprocal space coordinates of the q-points for the long-wavelength limit treatment.
!!  rhotwg(npwepG0*nspinor**2)=Oscillator matrix elements corresponding to an occupied-unoccupied pair of states.
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
!!  chi0sf(npwe,npwe,my_wl:my_wr)=Updated spectral function at q==0.
!!  lwing_sf(npwe,nomega,nqlwl)=Updated lower wing of the spectral function.
!!  uwing_sf(npwe,nomega,nqlwl)=Updated Upper wing of the spectral function.
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine assemblychi0sfq0(nqlwl,qlwl,ikbz,isym_kbz,itim_kbz,nspinor,symchi,npwepG0,npwe,Cryst,Ltg_q,Gsph_epsG0,&
& factocc,my_wl,iomegal,wl,my_wr,iomegar,wr,rhotwx,rhotwg,nomegasf,chi0sf,lwing_sf,uwing_sf)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'assemblychi0sfq0'
 use interfaces_32_util
 use interfaces_70_gw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikbz,my_wl,my_wr,nomegasf,npwe,npwepG0,nqlwl,nspinor
 integer,intent(in) :: isym_kbz,itim_kbz,symchi,iomegal,iomegar
 real(dp),intent(in) :: factocc,wl,wr
 type(littlegroup_t),intent(in) :: Ltg_q
 type(gsphere_t),target,intent(in) :: Gsph_epsG0
 type(crystal_t),intent(in) :: Cryst
!arrays
 real(dp),intent(in) :: qlwl(3,nqlwl)
 complex(gwpc),intent(inout) :: rhotwg(npwepG0*nspinor**2)
 complex(gwpc),intent(in) :: rhotwx(3)
 complex(gwpc),intent(inout) :: chi0sf(npwe,npwe,my_wl:my_wr)
 complex(dpc),intent(inout) :: lwing_sf(npwe,my_wl:my_wr,3)
 complex(dpc),intent(inout) :: uwing_sf(npwe,my_wl:my_wr,3)

!Local variables-------------------------------
!scalars
 integer :: itim,isym,iqlwl
 complex(gwpc) :: num
 complex(gwpc) :: mqg0,mqg0_sym,rhotwg0_bkp
 character(len=500) :: msg
!arrays
 integer, ABI_CONTIGUOUS pointer :: Sm1G(:)
 real(dp) :: opinv(3,3),qrot(3),b1(3),b2(3),b3(3)
 complex(gwpc) :: rhotwg_sym(npwe)
 complex(gwpc),allocatable :: rhotwg_sym_star(:),rhotwg_star(:)
 complex(gwpc),ABI_CONTIGUOUS pointer :: phmGt(:)
!************************************************************************

 if (iomegal<my_wl .or. iomegar>my_wr) then
   write(msg,'(3a,2(a,i0,a,i0))')ch10,&
&    'Indices out of boundary ',ch10,&
&    '  my_wl = ',my_wl,' iomegal = ',iomegal,ch10,&
&    '  my_wr = ',my_wr,' iomegar = ',iomegar,ch10
   MSG_BUG(msg)
 end if

 b1(:)=two_pi*Gsph_epsG0%gprimd(:,1)
 b2(:)=two_pi*Gsph_epsG0%gprimd(:,2)
 b3(:)=two_pi*Gsph_epsG0%gprimd(:,3)

 SELECT CASE (symchi)

 CASE (0)
   !
   ! === Calculation without symmetries ===
   ! * rhotwg(1)= R^-1q*rhotwx_ibz
   ! * rhotwg(1)=-R^-1q*conjg(rhotwx_ibz) for inversion
   ! FIXME My equation reads  -iq* <cSk|\nabla|vSk> = -i \transpose S <ck_i|\nabla\|vk_i>
   if (nspinor==1) then
     opinv(:,:)=REAL(Cryst%symrec(:,:,isym_kbz),dp)
     call matrginv(opinv,3,3)
     qrot =  (3-2*itim_kbz) * MATMUL(opinv,qlwl(:,1))
     rhotwg(1)=dotproductqrc(qrot,rhotwx,b1,b2,b3)
     if (itim_kbz==2) rhotwg(1) = GWPC_CONJG(rhotwg(1))

     if (wl<huge(0.0_dp)*1.d-11) then !this is awful but it is still a first coding
       num=-wl*factocc ! Num is single precision needed for cgerc check factocc
       call XGERC(npwe,npwe,num,rhotwg,1,rhotwg,1,chi0sf(:,:,iomegal),npwe)
     end if
     ! Last point, must accumulate left point but not the right one
     if (iomegar/=nomegasf+1 .and. wr<huge(0.0_dp)*1.d-11) then
       num=-wr*factocc
       call XGERC(npwe,npwe,num,rhotwg,1,rhotwg,1,chi0sf(:,:,iomegar),npwe)
     end if
     !
     ! === Accumulate heads and wings for each small q ===
     ! * For better performance, this part is not done if nqlwl==1
     !   lwing and uwing will be filled in cchi0q0 after the MPI collective sum
     !
     if (nqlwl>1.and..FALSE.) then
       rhotwg0_bkp = rhotwg(1) ! Save G=0 value of the first q
       ABI_MALLOC(rhotwg_star,(npwepG0))
       rhotwg_star = GWPC_CONJG(rhotwg(1:npwepG0))

       do iqlwl=2,nqlwl
         qrot =  (3-2*itim_kbz) * MATMUL(opinv,qlwl(:,iqlwl))
         mqg0 = dotproductqrc(qrot,rhotwx,b1,b2,b3) !TODO get rid of this
         if (itim_kbz==2) mqg0 = GWPC_CONJG(mqg0)
         rhotwg     (1) =mqg0
         rhotwg_star(1) =GWPC_CONJG(mqg0)
         !
         if (wl<huge(0.0_dp)*1.d-11) then !this is awful but it is still a first coding
           num=-wl*factocc ! Num is single precision needed for cgerc check factocc
           lwing_sf(:,iomegal,iqlwl) = lwing_sf(:,iomegal,iqlwl) + rhotwg     (1:npwepG0) * GWPC_CONJG(mqg0) * num
           uwing_sf(:,iomegal,iqlwl) = uwing_sf(:,iomegal,iqlwl) + rhotwg_star(1:npwepG0) *            mqg0  * num
         end if
         !
         ! Last point, must accumulate left point but not the right one
         if (iomegar/=nomegasf+1 .and. wr<huge(0.0_dp)*1.d-11) then
           num=-wr*factocc
           lwing_sf(:,iomegar,iqlwl) = lwing_sf(:,iomegar,iqlwl) + rhotwg     (1:npwepG0) * GWPC_CONJG(mqg0) * num
           uwing_sf(:,iomegar,iqlwl) = uwing_sf(:,iomegar,iqlwl) + rhotwg_star(1:npwepG0) *            mqg0  * num
         end if
       end do ! iqlwl

       ABI_FREE(rhotwg_star)
       rhotwg(1) = rhotwg0_bkp ! Reinstate previous value of rhotwg(1).
     end if !nqlwl

   else ! spinorial case
     MSG_BUG("Spectral method + nspinor==2 not implemented")
   end if


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
   !
   if (nspinor==1) then
     !ABI_MALLOC(rhotwg_sym,(npwe))
     !
     ! === Loop over symmetries of the space group and time-reversal ===
     do isym=1,Ltg_q%nsym_sg
       do itim=1,Ltg_q%timrev

         if (Ltg_q%wtksym(itim,isym,ikbz)==1) then
           ! === This operation belongs to the little group and has to be considered to reconstruct the BZ ===
           ! TODO this is a hot-spot, should add a test on the umklapp
           !
           phmGt => Gsph_epsG0%phmGt(1:npwe,isym) ! In these 2 lines mind the slicing (1:npwe)
           Sm1G  => Gsph_epsG0%rottbm1(1:npwe,itim,isym)

           opinv(:,:)=REAL(Cryst%symrec(:,:,isym),dp)
           call matrginv(opinv,3,3)
           qrot = (3-2*itim) * MATMUL(opinv,qlwl(:,1))

           SELECT CASE (itim)

           CASE (1)
             rhotwg_sym(1:npwe)=rhotwg(Sm1G(1:npwe))*phmGt(1:npwe)
             rhotwg_sym(1)=dotproductqrc(qrot,rhotwx,b1,b2,b3)

           CASE (2)
             rhotwg_sym(1:npwe) = GWPC_CONJG(rhotwg(Sm1G(1:npwe)))*phmGt(1:npwe)
             rhotwg_sym(1) = GWPC_CONJG(dotproductqrc(qrot,rhotwx,b1,b2,b3))

           CASE DEFAULT
             MSG_BUG(sjoin('Wrong value of itim= ', itoa(itim)))
           END SELECT
           !
           ! === Multiply elements G,Gp of rhotwg_sym*num and accumulate in chi0sf(G,Gp,io) ===
           if (wl<huge(0.0_dp)*1.d-11) then
             num=-wl*factocc
             call XGERC(npwe,npwe,num,rhotwg_sym,1,rhotwg_sym,1,chi0sf(:,:,iomegal),npwe)
           end if
           !
           ! Last point, must accumulate left point but not the right one
           if (iomegar/=nomegasf+1 .and. wr<huge(0.0_dp)*1.d-11) then
             num=-wr*factocc
             call XGERC(npwe,npwe,num,rhotwg_sym,1,rhotwg_sym,1,chi0sf(:,:,iomegar),npwe)
           end if

           ! === Accumulate heads and wings for each small q ===
           ! * For better performance, this part is not done if nqlwl==1
           !   lwing and uwing will be filled in cchi0q0 after the MPI collective sum
           if (nqlwl>1.and..FALSE.) then
             ABI_MALLOC(rhotwg_sym_star,(npwe))
             rhotwg_sym_star = GWPC_CONJG(rhotwg_sym(1:npwe))

             do iqlwl=2,nqlwl
               qrot =  (3-2*itim_kbz) * MATMUL(opinv,qlwl(:,iqlwl))
               mqg0_sym = dotproductqrc(qrot,rhotwx,b1,b2,b3) !TODO get rid of this
               if (itim_kbz==2) mqg0_sym=GWPC_CONJG(mqg0_sym)
               rhotwg_sym     (1) =mqg0_sym
               rhotwg_sym_star(1) =GWPC_CONJG(mqg0_sym)
               !
               if (wl<huge(0.0_dp)*1.d-11) then !this is awful but it is still a first coding
                 num=-wl*factocc ! Num is single precision needed for cgerc check factocc
                 lwing_sf(:,iomegal,iqlwl) = lwing_sf(:,iomegal,iqlwl) + rhotwg_sym_star(1:npwe) * GWPC_CONJG(mqg0_sym) * num
                 uwing_sf(:,iomegal,iqlwl) = uwing_sf(:,iomegal,iqlwl) + rhotwg_sym_star(1:npwe) *       mqg0_sym  * num
               end if
               ! Last point, must accumulate left point but not the right one
               if (iomegar/=nomegasf+1 .and. wr<huge(0.0_dp)*1.d-11) then
                 num=-wr*factocc
                 lwing_sf(:,iomegar,iqlwl) = lwing_sf(:,iomegar,iqlwl) + rhotwg_sym_star(1:npwe) * GWPC_CONJG(mqg0_sym) * num
                 uwing_sf(:,iomegar,iqlwl) = uwing_sf(:,iomegar,iqlwl) + rhotwg_sym_star(1:npwe) *            mqg0_sym  * num
               end if
             end do ! iqlwl

             ABI_FREE(rhotwg_sym_star)
           end if !nqlwl

         end if !wtksym
       end do !inv
     end do !isym

     !ABI_FREE(rhotwg_sym)

   else ! spinorial case
     MSG_BUG("Spectral method + nspinor==2 not implemented")
   end if

 CASE DEFAULT
   MSG_BUG(sjoin('Wrong value of symchi= ',itoa(symchi)))
 END SELECT

end subroutine assemblychi0sfq0
!!***

!----------------------------------------------------------------------

!!****f* m_chi0/symmetrize_afm_chi0
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
!!      wrtout
!!
!! SOURCE

subroutine symmetrize_afm_chi0(Cryst,Gsph,Ltg_q,npwe,nomega,chi0,chi0_head,chi0_lwing,chi0_uwing)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symmetrize_afm_chi0'
 use interfaces_14_hidewrite
!End of the abilint section

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
   !
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
     !
     do ig1=1,npwe ! Take care of diagonal.
       chi0(ig1,ig1,io)=two*chi0(ig1,ig1,io)
     end do
     !
     ! Upper and lower triangle are treated differently:
     ! We took advantage of the fact the afm_mat is hermitian to reduce memory.
     do ig2=2,npwe
       k0g=ig2*(ig2-1)/2
       do ig1=1,ig2-1
         kg=k0g+ig1
         chi0(ig1,ig2,io)=afm_mat(kg)*chi0(ig1,ig2,io)
       end do
     end do
     !
     do ig1=2,npwe
       k0g=ig1*(ig1-1)/2
       do ig2=1,ig1-1
         kg=k0g+ig2
         chi0(ig1,ig2,io)=CONJG(afm_mat(kg))*chi0(ig1,ig2,io)
       end do
     end do
     !
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

!!****f* m_chi0/accumulate_chi0_q0
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
!!  rhotwg(npwepG0*nspinor**2)=Oscillator matrix elements corresponding to an occupied-unoccupied pair of states.
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
!!      wrtout
!!
!! SOURCE

subroutine accumulate_chi0_q0(ik_bz,isym_kbz,itim_kbz,gwcomp,nspinor,npwepG0,Ep,Cryst,Ltg_q,Gsph_epsG0,&
& chi0,rhotwx,rhotwg,green_w,green_enhigh_w,deltaf_b1b2,chi0_head,chi0_lwing,chi0_uwing)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'accumulate_chi0_q0'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz,isym_kbz,itim_kbz,npwepG0,nspinor,gwcomp
 real(dp),intent(in) :: deltaf_b1b2
 type(littlegroup_t),intent(in) :: Ltg_q
 type(gsphere_t),target,intent(in) :: Gsph_epsG0
 type(crystal_t),intent(in) :: Cryst
 type(em1params_t),intent(in) :: Ep
!arrays
 complex(gwpc),intent(in) :: rhotwg(npwepG0*nspinor**2)
 complex(gwpc),intent(in) :: rhotwx(3,nspinor**2)
 complex(gwpc),intent(inout) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)
 complex(dpc),intent(in) :: green_w(Ep%nomega),green_enhigh_w(Ep%nomega)
 complex(dpc),intent(inout) :: chi0_head(3,3,Ep%nomega)
 complex(dpc),intent(inout) :: chi0_lwing(Ep%npwe*Ep%nI,Ep%nomega,3)
 complex(dpc),intent(inout) :: chi0_uwing(Ep%npwe*Ep%nJ,Ep%nomega,3)

!Local variables-------------------------------
!scalars
 integer :: itim,io,isym,idir,jdir,jj,ii,s_jj,pad_jj,pad_ii
 complex(gwpc) :: dd
 !character(len=500) :: msg
!arrays
 integer,ABI_CONTIGUOUS pointer :: Sm1G(:)
 complex(dpc) :: mir_kbz(3)
 complex(gwpc),allocatable :: rhotwg_sym(:),rhotwg_I(:),rhotwg_J(:)
 complex(gwpc), ABI_CONTIGUOUS pointer :: phmGt(:)

!************************************************************************

 ABI_UNUSED(deltaf_b1b2)

 SELECT CASE (Ep%symchi)

 CASE (0) ! Do not use symmetries.

   if (nspinor==1) then
    ! * Symmetrize rhotwg in the full BZ and accumulate over the full BZ i.e.
    !     chi0(G1,G2,io) = chi0(G1,G2,io) + (rhotwg(G1)*CONJG(rhotwg(G2)))*green_w(io)
    ! * The non-analytic term is symmetrized for this k-point in the BZ according to:
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
    mir_kbz =(3-2*itim_kbz) * MATMUL(Cryst%symrec(:,:,isym_kbz),rhotwx(:,1))
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
    !
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

  else ! spinorial case
    ABI_MALLOC(rhotwg_I,(Ep%npwe))
    ABI_MALLOC(rhotwg_J,(Ep%npwe))

    MSG_WARNING("Check tensor")

    ! I can use symmetries to loop over the upper triangle but
    ! this makes using BLAS more difficult
    ! Important NOTE: treatment of q-->0 limit is correct only
    ! for i=j=0. Other components require additional terms.

    do jj=1,Ep%nJ
      s_jj=1 ; if (jj==4) s_jj=-1
      pad_jj=(jj-1)*Ep%npwe
      call mkrhotwg_sigma(jj,nspinor,Ep%npwe,rhotwg,rhotwg_J)

      !rhotwg_J(1) = q0limit(jj,qlwl(:,1),nspinor,rhotwx,b1,b2,b3)
      !TODO RECHECK this
      if (itim_kbz==2) rhotwg_J(1)=-CONJG(rhotwg_J(1))

      do ii=1,Ep%nI
        pad_ii=(ii-1)*Ep%npwe

        if (ii/=jj) then
          call mkrhotwg_sigma(ii,nspinor,Ep%npwe,rhotwg,rhotwg_I)
          !rhotwg_I(1) = q0limit(ii,qlwl(:,1),nspinor,rhotwx,b1,b2,b3)
          if (itim_kbz==2) rhotwg_I(1)=-CONJG(rhotwg_I(1))
        else
          rhotwg_I(:)=rhotwg_J(:)
        end if

        do io=1,Ep%nomega
          dd = s_jj*green_w(io)
          call XGERC(Ep%npwe,Ep%npwe,dd,rhotwg_I,1,rhotwg_J,1,&
&           chi0(pad_ii+1:pad_ii+Ep%npwe,pad_jj+1:pad_jj+Ep%npwe,io),Ep%npwe)
        end do

      end do !ii
    end do !jj

    ABI_FREE(rhotwg_I)
    ABI_FREE(rhotwg_J)
  end if

 CASE (1) ! Use symmetries to reconstruct the integrand.
   if (nspinor==1) then
     ABI_MALLOC(rhotwg_sym,(Ep%npwe))

     ! === Loop over symmetries of the space group and time-reversal ===
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
           !
           ! Multiply elements G1,G2 of rhotwg_sym by green_w(io) and accumulate in chi0(G,Gp,io)
!$omp parallel do private(dd)
           do io=1,Ep%nomega
             dd=green_w(io)
             call XGERC(Ep%npwe,Ep%npwe,dd,rhotwg_sym,1,rhotwg_sym,1,chi0(:,:,io),Ep%npwe)
           end do

           ! === Accumulate heads and wings for each small q ===
           ! FIXME extrapolar method should be checked!!

           ! Symmetrize <r> in full BZ: <Sk b|r|Sk b'> = R <k b|r|k b'> + \tau \delta_{bb'}
           mir_kbz =(3-2*itim) * MATMUL(Cryst%symrec(:,:,isym),rhotwx(:,1))
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
           !
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

   else  !spinorial case
     MSG_ERROR('symchi=1 with spinor not implemented ')
     MSG_ERROR("Check tensor")
   end if

 CASE DEFAULT
   MSG_BUG(sjoin('Wrong value of symchi ',itoa(Ep%symchi)))
 END SELECT

end subroutine accumulate_chi0_q0
!!***

!----------------------------------------------------------------------

!!****f* assemblychi0q0_sym/q0limit
!! NAME
!!  q0limit
!!
!! FUNCTION
!!  Return an appropriate linear combination of the spin-dependent oscilator matrix elements at G=0, taking
!!  into account the small q used for the treatment of the optical limit. Used when nspinor=2.
!!
!! INPUTS
!! ii=Index defining the required linear combination of the oscillator strengths.
!! nspinor=Number of spinorial components.
!! qlwl(3)=Reduced components of the small q around Gamma for the treatment of the long wave-length limit.
!! b1(3),b2(3),b3(3)=Lattice vectore of the reciprocal lattice.
!! rhotwx(3,nspinor**2)=Oscillator matrix elements at G=0, for each possible combination of the left and right spin component.
!!
!! OUTPUT
!!  q0limit=Linear combination, see doc below.
!!
!! TODO
!!  This function should be "contained" to facilitate inlining but abilint crashes, dont know why!
!!
!! PARENTS
!!
!! SOURCE

function q0limit(ii,qlwl,nspinor,rhotwx,b1,b2,b3)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'q0limit'
 use interfaces_70_gw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ii,nspinor
 complex(gwpc) :: q0limit
!arrays
 real(dp),intent(in) :: qlwl(3)
 real(dp),intent(in) :: b1(3),b2(3),b3(3)
 complex(gwpc),intent(in) :: rhotwx(3,nspinor**2)

! *********************************************************************

 SELECT CASE (ii)
 CASE (1) ! M_0(q-->0) = Lim M(up,up)+M(dwn,dwn). Exact, neglecting Vnl
   q0limit =  dotproductqrc(qlwl,rhotwx(:,1),b1,b2,b3) &
&            +dotproductqrc(qlwl,rhotwx(:,2),b1,b2,b3)

 CASE (2) ! M_z(q-->0) = Lim M(up,up)-M(dwn,dwn).
   ! WARNING off-diagonal elements of rV12 and rV12 are neglected
   q0limit =  dotproductqrc(qlwl,rhotwx(:,1),b1,b2,b3) &
&            -dotproductqrc(qlwl,rhotwx(:,2),b1,b2,b3)

 CASE (3) ! M_x(q-->0) = M(up,dwn)+M(dwn,up).
   ! Both diagonal elements of the form v12r-rv21 and similiar terms in 12 and 21 are neglected
   q0limit =  dotproductqrc(qlwl,rhotwx(:,3),b1,b2,b3) &
&            +dotproductqrc(qlwl,rhotwx(:,4),b1,b2,b3)

 CASE (4)
   ! Both diagonal elements of the form v12r-rv21 and similiar terms in 12 and 21 are neglected
   q0limit =( dotproductqrc(qlwl,rhotwx(:,3),b1,b2,b3) &
&            -dotproductqrc(qlwl,rhotwx(:,4),b1,b2,b3) )*j_gw

 CASE DEFAULT
   MSG_BUG(sjoin('Wrong value for ii= ',itoa(ii)))
 END SELECT

end function q0limit
!!***

!----------------------------------------------------------------------

!!****f* m_chi0/accumulate_sfchi0_q0
!! NAME
!! accumulate_sfchi0_q0
!!
!! FUNCTION
!! Update the spectral function of the independent particle susceptibility at q==0 for the contribution
!! of one pair of occupied-unoccupied band, for each frequency.
!! If symchi==1, the symmetries belonging to the little group of the external point q are used
!! to reconstrunct the contributions in the full Brillouin zone. In this case, the equation implented is:
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
!!  rhotwg(npwepG0*nspinor**2)=Oscillator matrix elements corresponding to an occupied-unoccupied pair of states.
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
!!      wrtout
!!
!! SOURCE

subroutine accumulate_sfchi0_q0(ikbz,isym_kbz,itim_kbz,nspinor,symchi,npwepG0,npwe,Cryst,Ltg_q,Gsph_epsG0,&
& factocc,my_wl,iomegal,wl,my_wr,iomegar,wr,rhotwx,rhotwg,nomegasf,sf_chi0,sf_head,sf_lwing,sf_uwing)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'accumulate_sfchi0_q0'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikbz,my_wl,my_wr,nomegasf,npwe,npwepG0,nspinor
 integer,intent(in) :: isym_kbz,itim_kbz,symchi,iomegal,iomegar
 real(dp),intent(in) :: factocc,wl,wr
 type(littlegroup_t),intent(in) :: Ltg_q
 type(gsphere_t),target,intent(in) :: Gsph_epsG0
 type(crystal_t),intent(in) :: Cryst
!arrays
 complex(gwpc),intent(in) :: rhotwg(npwepG0*nspinor**2)
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
   ! === Calculation without symmetries ===
   ! * rhotwg(1)= R^-1q*rhotwx_ibz
   ! * rhotwg(1)=-R^-1q*conjg(rhotwx_ibz) for inversion
   if (nspinor==1) then

     if (wl<huge(0.0_dp)*1.d-11) then !this is awful but it is still a first coding
       num=-wl*factocc ! Num is single precision needed for cgerc check factocc
       call XGERC(npwe,npwe,num,rhotwg,1,rhotwg,1,sf_chi0(:,:,iomegal),npwe)
     end if
     ! Last point, must accumulate left point but not the right one
     if (iomegar/=nomegasf+1 .and. wr<huge(0.0_dp)*1.d-11) then
       num=-wr*factocc
       call XGERC(npwe,npwe,num,rhotwg,1,rhotwg,1,sf_chi0(:,:,iomegar),npwe)
     end if
     !
     ! === Accumulate heads and wings for each small q ===
     ! Symmetrize <r> in full BZ: <Sk b|r|Sk b'> = R <k b|r|k b'> + \tau \delta_{bb'}
     mir_kbz =(3-2*itim_kbz) * MATMUL(Cryst%symrec(:,:,isym_kbz),rhotwx(:,1))
     if (itim_kbz==2) mir_kbz=CONJG(mir_kbz)
     !
     ! ================================
     ! ==== Update heads and wings ====
     ! ================================
     do jdir=1,3
       !
       if (wl<huge(0.0_dp)*1.d-11) then !this is awful but it is still a first coding
         num=-wl*factocc ! Num is single precision needed for cgerc check factocc
         sf_uwing(:,iomegal,jdir) = sf_uwing(:,iomegal,jdir) + num * mir_kbz(jdir) * CONJG(rhotwg(1:npwepG0))
         sf_lwing(:,iomegal,jdir) = sf_lwing(:,iomegal,jdir) + num * rhotwg(1:npwepG0) * CONJG(mir_kbz(jdir))
         do idir=1,3
           sf_head(idir,jdir,iomegal) = sf_head(idir,jdir,iomegal) + num * mir_kbz(idir) * CONJG(mir_kbz(jdir))
         end do
       end if
       !
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

   else ! spinorial case
     MSG_BUG("Spectral method + nspinor==2 not implemented")
   end if


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
   !
   if (nspinor==1) then
     ABI_MALLOC(rhotwg_sym,(npwe))
     !
     ! === Loop over symmetries of the space group and time-reversal ===
     do isym=1,Ltg_q%nsym_sg
       do itim=1,Ltg_q%timrev

         if (Ltg_q%wtksym(itim,isym,ikbz)==1) then
           ! === This operation belongs to the little group and has to be considered to reconstruct the BZ ===
           !
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
           !
           ! === Multiply elements G,Gp of rhotwg_sym*num and accumulate in sf_chi0(G,Gp,io) ===
           if (wl<huge(0.0_dp)*1.d-11) then
             num=-wl*factocc
             call XGERC(npwe,npwe,num,rhotwg_sym,1,rhotwg_sym,1,sf_chi0(:,:,iomegal),npwe)
           end if
           !
           ! Last point, must accumulate left point but not the right one
           if (iomegar/=nomegasf+1 .and. wr<huge(0.0_dp)*1.d-11) then
             num=-wr*factocc
             call XGERC(npwe,npwe,num,rhotwg_sym,1,rhotwg_sym,1,sf_chi0(:,:,iomegar),npwe)
           end if

           ! === Accumulate heads and wings for each small q ===
           ! Symmetrize <r> in full BZ: <Sk b|r|Sk b'> = R <k b|r|k b'> + \tau \delta_{bb'}
           mir_kbz =(3-2*itim) * MATMUL(Cryst%symrec(:,:,isym),rhotwx(:,1))
           if (itim==2) mir_kbz=CONJG(mir_kbz)

           do jdir=1,3
             !
             if (wl<huge(0.0_dp)*1.d-11) then !this is awful but it is still a first coding
               num=-wl*factocc ! Num is single precision needed for cgerc check factocc
               sf_uwing(:,iomegal,jdir) = sf_uwing(:,iomegal,jdir) + num * mir_kbz(jdir) * CONJG(rhotwg_sym(1:npwe))
               sf_lwing(:,iomegal,jdir) = sf_lwing(:,iomegal,jdir) + num * rhotwg_sym(1:npwe) * CONJG(mir_kbz(jdir))
               do idir=1,3
                 sf_head(idir,jdir,iomegal) = sf_head(idir,jdir,iomegal) + num * mir_kbz(idir) * CONJG(mir_kbz(jdir))
               end do
             end if
             !
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

   else ! spinorial case
     MSG_BUG("Spectral method + nspinor==2 not implemented")
   end if

 CASE DEFAULT
   MSG_BUG(sjoin('Wrong value of symchi:', itoa(symchi)))
 END SELECT

end subroutine accumulate_sfchi0_q0
!!***

!----------------------------------------------------------------------

!!****f* m_chi0/assemblychi0sf
!! NAME
!! assemblychi0sf
!!
!! FUNCTION
!! Update the spectral function of the irreducible polarizability for the contribution
!! of one pair of occupied-unoccupied states, for each frequenciy.
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
!!  rhotwg(npwepG0*nspinor**2)=Oscillator matrix elements corresponding to an occupied-unoccupied pair of states.
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
!!      wrtout
!!
!! SOURCE

subroutine assemblychi0sf(ik_bz,nspinor,symchi,Ltg_q,npwepG0,npwe,rhotwg,Gsph_epsG0,&
& factocc,my_wl,iomegal,wl,my_wr,iomegar,wr,nomegasf,chi0sf)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'assemblychi0sf'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz,iomegal,iomegar,my_wl,my_wr,nomegasf,npwe,npwepG0
 integer,intent(in) :: nspinor,symchi
 real(dp),intent(in) :: factocc,wl,wr
 type(gsphere_t),intent(in) :: Gsph_epsG0
 type(littlegroup_t),intent(in) :: Ltg_q
!arrays
 complex(gwpc),intent(in) :: rhotwg(npwepG0*nspinor**2)
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

 CASE (0)  ! Do not use symmetries.

! MG: This is the best I can do for this part.
!$omp PARALLEL private(num)
!$omp SECTIONS

!$omp SECTION
   if (wl<huge(0.0_dp)*1.d-11) then !FIXME this is awful
     num=-wl*factocc
     call XGERC(npwe,npwe,num,rhotwg,1,rhotwg,1,chi0sf(:,:,iomegal),npwe)
   end if
   !
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
   ! === Loop over symmetries of the space group and time-reversal ===
   do isym=1,Ltg_q%nsym_sg
     do itim=1,Ltg_q%timrev

       if (Ltg_q%wtksym(itim,isym,ik_bz)==1) then
         ! === This operation belongs to the little group and has to be used to reconstruct BZ ===
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
         !
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

!!****f* m_chi0/approxdelta
!! NAME
!!  approxdelta
!!
!! FUNCTION
!!  Approximate the Dirac function using two methods:
!!  method 1) a triangular funtion centered at the value egwdiff_re (Eq 17 of PRB 74, 035101 (2006)
!!  method 2) a gaussian of witdth ep%spsmear expandended in Taylor series
!!  (at the moment only the 0-th moments)
!!
!!  Subroutine needed to implement the calculation
!!  of the polarizability using the spectral representation as proposed in :
!!  PRB 74, 035101 (2006) and PRB 61, 7172 (1999)
!!
!! INPUTS
!!  nomegasf=number of frequencies in the grid for Im \chi_0
!!  omegasf(0:nomega+1)= frequencies (real)
!!  egwdiff_re = transition energy where the delta function is centered
!!
!!  method= 1 : a triangular shaped function used to approximated the delta
!!          2 : gaussian approximation with standard deviation (smear)
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
!!      wrtout
!!
!! SOURCE

subroutine approxdelta(nomegasf,omegasf,egwdiff_re,smear,iomegal,iomegar,wl,wr,spmeth)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'approxdelta'
!End of the abilint section

 implicit none

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

 iomegal=iomega    ; omegal=omegasf(iomegal)
 iomegar=iomegal+1 ; omegar=omegasf(iomegar)

 SELECT CASE (spmeth)

 CASE (1) ! Weights for triangular shaped function
   wr=  (egwdiff_re-omegal)/(omegar-omegal)
   wl= -(egwdiff_re-omegar)/(omegar-omegal)

 CASE (2) ! Weights for gaussian method (0-th moment)
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

!!****f* m_chi0/calc_kkweight
!! NAME
!!  calc_kkweight
!!
!! FUNCTION
!!  Calculate frequency dependent weights needed to perform the Hilbert transform
!!
!!  Subroutine needed to implement the calculation
!!  of the polarizability using the spectral representation as proposed in :
!!  PRB 74, 035101 (2006) and PRB 61, 7172 (1999)
!!
!! INPUTS
!! nsp=number of frequencies where the imaginary part of the polarizability is evaluated
!! ne=number of frequencies for the polarizability (same as in epsilon^-1)
!! omegasp(nsp)=real frequencies for the imaginary part of the polarizability
!! omegae(ne)= imaginary frequencies for the polarizability
!! delta=small imaginary part used to avoid poles, input variables
!!
!! OUTPUT
!! kkweight(nsp,ne)=frequency dependent weights (Eq A1 PRB 74, 035101 (2006)
!!
!! PARENTS
!!      m_chi0
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE
!!

subroutine calc_kkweight(ne,omegae,nsp,omegasp,delta,omegamax,kkw)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_kkweight'
!End of the abilint section

 implicit none

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
 real(dp) :: eta,xx1,xx2
 real(dp) :: den1,den2
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
     if (isp==nsp) then ! Skip last point should check that this would not lead to spurious effects
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

!!****f* m_chi0/setup_spectral
!! NAME
!!  setup_spectral
!!
!! FUNCTION
!! Calculation of \chi_o based on the spectral method as proposed in PRB 74, 035101 (2006) and PRB 61, 7172 (1999).
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
!!      wrtout
!!
!! SOURCE

subroutine setup_spectral(nomega,omega,nomegasf,omegasf,max_rest,min_rest,my_max_rest,my_min_rest,&
&  method,zcut,omegaplasma,my_wl,my_wr,kkweight)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'setup_spectral'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

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

 ! === The mesh has to enclose the entire range of transitions ===
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
 !
 ! ======================================================
 ! === Setup of the w-mesh for the spectral function ====
 ! ======================================================
 SELECT CASE (method)

 CASE (0) ! Linear mesh.
   call wrtout(std_out,' Using linear mesh for Im chi0','COLL')
   do io=1,nomegasf
     omegasf(io)=(io-1)*domegasf+min_rest-dd
   end do

 CASE (1) ! Non-homogeneous mesh densified around omega_plasma, do not improve results ===
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
 !
 ! === Find min and max index in omegasf treated by this processor ===
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
 !
 ! Calculate weights for Hilbert transform.
 call calc_kkweight(nomega,omega,nomegasf,omegasf,zcut,max_rest,kkweight)

end subroutine setup_spectral
!!***

!----------------------------------------------------------------------

!!****f* m_chi0/hilbert_transform
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
!!      wrtout
!!
!! SOURCE

subroutine hilbert_transform(npwe,nomega,nomegasf,my_wl,my_wr,kkweight,sf_chi0,chi0,spmeth)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hilbert_transform'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spmeth,nomega,nomegasf
 integer,intent(in) :: my_wl,my_wr,npwe
!arrays
 complex(dpc),intent(in) :: kkweight(nomegasf,nomega)
 complex(gwpc), intent(inout) :: sf_chi0(npwe,npwe,my_wl:my_wr)
 complex(gwpc), intent(inout) :: chi0(npwe,npwe,nomega)

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
 call wrtout(std_out,msg,'COLL',do_flush=.True.)

 my_nwp = my_wr-my_wl +1

!$omp parallel private(my_kkweight, A_g1wp, H_int, ig2)
 ABI_MALLOC(my_kkweight,(my_wl:my_wr,nomega))
 my_kkweight = kkweight(my_wl:my_wr,:)

 ABI_MALLOC(A_g1wp,(npwe,my_nwp))
 ABI_MALLOC(H_int,(npwe,nomega))

!$omp do
 do ig2=1,npwe
   A_g1wp = sf_chi0(:,ig2,:)
   !
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

!!****f* m_chi0/hilbert_transform_headwings
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
!!      wrtout
!!
!! SOURCE

subroutine hilbert_transform_headwings(npwe,nomega,nomegasf,my_wl,my_wr,kkweight, &
& sf_lwing,sf_uwing,sf_head,chi0_lwing,chi0_uwing,chi0_head,spmeth)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hilbert_transform_headwings'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spmeth,nomega,nomegasf
 integer,intent(in) :: my_wl,my_wr,npwe
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
 !
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

 ! * Hilbert transform for wings.
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

END MODULE m_chi0
!!***
