!{\src2tex{textfont=tt}}
!!****f* ABINIT/completechi0_deltapart
!! NAME
!! completechi0_deltapart
!!
!! FUNCTION
!!  Apply the delta part of the completeness correction to chi0
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (FB, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see
!! ~abinit/doc/developers/contributors.txt .
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine completechi0_deltapart(ik_bz,qzero,symchi,npwe,npwvec,nomega,nspinor,&
& nfftot,ngfft,igfft0,Gsph_FFT,Ltg_q,green_enhigh_w,wfwfg,chi0)

 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_gsphere,        only : gsphere_t, gsph_gmg_idx, gsph_gmg_fftidx
 use m_bz_mesh,        only : littlegroup_t

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'completechi0_deltapart'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

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
 character(len=500) :: msg

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
     write(msg,'(a,i0)')' Number of G1-G2 pairs outside the G-sphere for Wfns = ',outofbox_wfn
     MSG_WARNING(msg)
     if (enough==50) then
       call wrtout(std_out,' ========== Stop writing Warnings ==========','COLL') 
     end if
   end if
 end if

end subroutine completechi0_deltapart
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/output_chi0sumrule
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

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'output_chi0sumrule'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

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

!!****f* ABINIT/accumulate_chi0sumrule
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

 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_gsphere,  only : gsphere_t
 use m_bz_mesh,  only : littlegroup_t

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'accumulate_chi0sumrule'
!End of the abilint section

 implicit none

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
 character(len=500) :: msg
!arrays
 integer,allocatable :: Sm1_gmG0(:)
 integer, ABI_CONTIGUOUS pointer :: gmG0(:)
 complex(gwpc),allocatable :: rhotwg_sym(:)

!************************************************************************

 ! Accumulating the sum rule on chi0.
 ! Eq.(5.284) in G. D. Mahan Many-Particle Physics 3rd edition

 SELECT CASE (symchi)

 CASE (0) ! Do not use symmetries, sum is performed in the full BZ.
   chi0sumrule(:)=chi0sumrule(:) + factor*delta_ene*ABS(rhotwg(1:npwe))**2

 CASE (1) ! Symmetrize the contribution in the full BZ.
   ABI_ALLOCATE(rhotwg_sym,(npwe))
   ABI_ALLOCATE(Sm1_gmG0,(npwe))

   do itim=1,Ltg_q%timrev
     do isym=1,Ltg_q%nsym_sg
       if (Ltg_q%wtksym(itim,isym,ik_bz)==1) then 
        ! === This operation belongs to the little group and has to be used to reconstruct the BZ ===
        ! * In the following 2 lines mind the slicing (1:npwe)
        !
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
   write(msg,'(a,i3)')'Wrong value for symchi= ',symchi
   MSG_BUG(msg)
 END SELECT

end subroutine accumulate_chi0sumrule
!!***


!subroutine check_completeness(ik_bz,jb,iq,ig,igp,npwepG0)
!!scalars
! integer,intent(in) :: ik_bz,jb,iq,ig,igp,npwepG0
!!arrays
! complex(gwpc),intent(in) :: rhotwg(npwepG0)
!
!! Delta part
!
!call rho_tw_g(nspinor,npwvec,nr,ngfft,ndat1,map2sphere,use_padfft,igfftg0,gbound,&
!& wfn1,i1,ktabr1,ktabp1,spinrot1,&
!& wfn2,i2,ktabr2,ktabp2,spinrot2,&
!& dim_rtwg,rhotwg)
!
!call paw_rho_tw_g(npw,dim_rtwg,nspinor,natom,ntypat,typat,xred,gvec,Cprj_kmqb1,Cprj_kb2,Pwij,rhotwg)
!
!call delta_paw_rho_tw_g()  ! This doesn't exist...yet
!
!! Copied from cchi0
!! maybe the whole thing should be in cchi0...
!
!         ! ==== Form rho-twiddle(r)=u^*_{b1,kmq_bz}(r) u_{b2,kbz}(r) and its FFT transform ====
!         call rho_tw_g(nspinor,Ep%npwepG0,nfftot_gw,ndat1,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound,&
!&          wfr1,itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,&
!&          wfr2,itim_k  ,tabr_k  ,ph_mkt  ,spinrot_k,&
!&          dim_rtwg,rhotwg)
!
!         if (Psps%usepaw==1) then! Add PAW on-site contribution, projectors are already in the BZ.
!           call paw_rho_tw_g(Ep%npwepG0,dim_rtwg,nspinor,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xred,Gsph_epsG0%gvec,&
!&           Cprj1_kmq,Cprj2_k,Pwij,rhotwg)
!
!          if (deltapaw==1) then  ! Add the term arising from the incompleteness of the basis
!            call delta_paw_rho_tw_g()  ! This doesn't exist...yet
!          end if
!         end if
!
!end subroutine check_completeness
