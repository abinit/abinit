!{\src2tex{textfont=tt}}
!!****f* ABINIT/wfd_mkrho
!! NAME
!! wfd_mkrho
!!
!! FUNCTION
!! Calculate the charge density on the fine FFT grid in real space.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (MG,GMR, VO, LR, RWG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ngfftf(18)=array containing all the information for the "fine" FFT.
!!  Cryst<crystal_t> Info on the crystalline structure
!!     %nsym=number of symmetry operations.
!!     %ucvol=unit cell volume.
!!  optcalc=option for calculation. If =0 (default value) perform calculation
!!    of electronic density. If =1, perform calculation of kinetic energy density.
!!    In both cases, the result is returned in rhor.
!!  Psps<type(pseudopotential_type)>=variables related to pseudopotentials
!!  nfftf=Total number of points on the fine FFT grid (for this processor)
!!  Kmesh<kmesh_t>= Info on the k-sampling:
!!     %nibz=number of irreducible k-points.
!!     %nbz=number of k-points in the full Brillouin zone.
!!     %wt(nibz)=irreducible k-points weights.
!!     %timrev=2 if time-reversal symmetry can be used, 1 otherwise.
!!  Wfd<wfd_t)=datatype gathering info on the wavefunctions.
!!    %npwwfn=Number of plane waves used to describe the wave functions.
!!    %nspinor=number of spinorial components.
!!    %nsppol=1 for unpolarized, 2 for spin-polarized calculations.
!!    %nspden=number of spin-density components.
!! [optcalc]=Optional option used to calculate the kinetic energy density. Defaults to 0.
!!
!! OUTPUT
!!  rhor(nfftf,nspden)=The density in the real space on the fine FFT grid.
!!   If nsppol==2, total charge in first half, spin-up component in second half.
!!   (for non-collinear magnetism, first element: total density, 3 next ones: mx,my,mz in units of hbar/2)
!!   If optcalc==1 (optional argument, default value is 0), then rhor will actually
!!   contains kinetic energy density (taur) instead of electronic density.
!!
!! NOTES
!! In the case of PAW calculations:
!!    All computations are done on the fine FFT grid.
!!    All variables (nfftf,ngfftf,mgfftf) refer to this fine FFT grid.
!!    All arrays (densities/potentials...) are computed on this fine FFT grid.
!!    Developers have to be careful when introducing others arrays:
!!      they have to be stored on the fine FFT grid.
!! In the case of norm-conserving calculations:
!!    The mesh is the usual augmented FFT grid to treat correctly the convolution.
!!
!! PARENTS
!!      bethe_salpeter,screening,sigma
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wfd_mkrho(Wfd,Cryst,Psps,Kmesh,Bands,ngfftf,nfftf,rhor,&
&                    optcalc) ! optional arguments


 use defs_basis
 use defs_datatypes
 use m_profiling_abi
 use m_xmpi
 use m_errors

 use m_iterators, only : iter2_t, iter_yield, iter_len, iter_free
 use m_crystal,   only : crystal_t
 use m_bz_mesh,   only : kmesh_t
 use m_fft,       only : fft_ug
 use m_wfd,       only : wfd_t, wfd_get_ur, wfd_update_bkstab, wfd_iterator_bks, wfd_change_ngfft

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_mkrho'
 use interfaces_14_hidewrite
 use interfaces_56_recipspace
 use interfaces_67_common
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftf
 integer,intent(in),optional :: optcalc
 type(ebands_t),intent(in) :: Bands
 type(kmesh_t),intent(in) :: Kmesh
 type(crystal_t),intent(in) :: Cryst
 type(Pseudopotential_type),intent(in) :: Psps
 type(wfd_t),intent(inout) :: Wfd
!arrays
 integer,intent(in) :: ngfftf(18)
 real(dp),intent(out) :: rhor(nfftf,Wfd%nspden)

!Local variables ------------------------------
!scalars
 integer,parameter :: ndat1=1
 integer :: cplex,ib,ib_iter,ierr,ik,ir,is,n1,n2,n3,nfftotf
 integer :: alpha,nalpha,ipw!,ipwsp
 integer :: myoptcalc
 real(dp) :: kpt_cart,kg_k_cart,gp2pi1,gp2pi2,gp2pi3,cwftmp,bks_weight
 character(len=500) :: msg
!arrays
 integer,allocatable :: irrzon(:,:,:)
 real(dp),allocatable :: phnons(:,:,:),rhog(:,:),rhor_down(:),rhor_mx(:),rhor_my(:)
 real(dp),allocatable :: cwavef(:,:)
 complex(dpc),allocatable :: wfr_x(:),wfr_y(:)
 complex(gwpc),allocatable :: gradug(:),work(:)
 complex(gwpc),allocatable,target :: wfr(:)
 complex(gwpc), ABI_CONTIGUOUS pointer :: cwavef1(:),cwavef2(:)
 type(iter2_t) :: Iter_bks

!*************************************************************************

 DBG_ENTER("COLL")

 ! Consistency check.
 if (Wfd%nspden==1.and.Wfd%nspinor==2) then
   MSG_ERROR('nspden==1 and nspinor==2 not implemented')
 end if

 ABI_CHECK(Wfd%nsppol==Bands%nsppol,"mismatch in nspppol")

 if ( ANY(ngfftf(1:3) /= Wfd%ngfft(1:3)) ) then
   call wfd_change_ngfft(Wfd,Cryst,Psps,ngfftf)
 end if
 !
 ! === Calculate IBZ contribution to the charge density ===
 ABI_MALLOC(wfr,(nfftf*Wfd%nspinor))

 if (Wfd%nspinor==2) then
   ABI_MALLOC(wfr_x,(nfftf))
   ABI_MALLOC(wfr_y,(nfftf))
   if (Wfd%nspden==4) then
     ABI_MALLOC(rhor_down,(nfftf))
     ABI_MALLOC(rhor_mx,(nfftf))
     ABI_MALLOC(rhor_my,(nfftf))
     rhor_down=zero
     rhor_mx  =zero
     rhor_my  =zero
   else !TODO
     MSG_ERROR('nspden/=4 and nspinor=2 not implemented')
   end if
 end if

 ! Update the (b,k,s) distribution table.
 call wfd_update_bkstab(Wfd)

 ! Calculate the unsymmetrized density.
 rhor=zero
 myoptcalc=0;  if (present(optcalc)) myoptcalc=optcalc
 nalpha=1;  if (myoptcalc==1) nalpha=3
 if (myoptcalc == 1 .and. wfd%nspinor == 2) then
   MSG_ERROR("kinetic energy density with nspinor == 2 not implemented")
 end if

 ! Build the iterator that will distribute the states in an automated way.
 Iter_bks = wfd_iterator_bks(Wfd,bks_mask=ABS(Bands%occ)>=tol8)

 do alpha=1,nalpha
   do is=1,Wfd%nsppol
     do ik=1,Wfd%nkibz
       !
       do ib_iter=1,iter_len(Iter_bks,ik,is)
         ib = iter_yield(Iter_bks,ib_iter,ik,is)
         bks_weight = Bands%occ(ib,ik,is) * Kmesh%wt(ik) / Cryst%ucvol

         call wfd_get_ur(Wfd,ib,ik,is,wfr)

         cwavef1 => wfr(1:nfftf)
         if(myoptcalc==1) then
           ABI_MALLOC(gradug,(Wfd%Kdata(ik)%npw))
           ABI_MALLOC(cwavef,(2,Wfd%Kdata(ik)%npw))
           ABI_MALLOC(work,(nfftf))
           cwavef(1,:)= REAL(Wfd%Wave(ib,ik,is)%ug(:))
           cwavef(2,:)=AIMAG(Wfd%Wave(ib,ik,is)%ug(:))
!          Multiplication by 2pi i (k+G)_alpha
           gp2pi1=Cryst%gprimd(alpha,1)*two_pi
           gp2pi2=Cryst%gprimd(alpha,2)*two_pi
           gp2pi3=Cryst%gprimd(alpha,3)*two_pi
           kpt_cart=gp2pi1*Wfd%kibz(1,ik)+gp2pi2*Wfd%kibz(2,ik)+gp2pi3*Wfd%kibz(3,ik)
           do ipw=1,Wfd%Kdata(ik)%npw
             kg_k_cart= gp2pi1*Wfd%Kdata(ik)%kg_k(1,ipw) + &
&                       gp2pi2*Wfd%Kdata(ik)%kg_k(2,ipw) + &
&                       gp2pi3*Wfd%Kdata(ik)%kg_k(3,ipw)+kpt_cart
!             ipwsp=ipw!+(ispinor-1)*Wfd%Kdata(ik)%npw
             cwftmp=-cwavef(2,ipw)*kg_k_cart
             cwavef(2,ipw)=cwavef(1,ipw)*kg_k_cart
             cwavef(1,ipw)=cwftmp
           end do
           gradug(:)=CMPLX(cwavef(1,:),cwavef(2,:),gwpc)
           call fft_ug(Wfd%npwarr(ik),nfftf,Wfd%nspinor,ndat1,Wfd%mgfft,Wfd%ngfft,&
&            Wfd%istwfk(ik),Wfd%Kdata(ik)%kg_k,Wfd%Kdata(ik)%gbound,gradug,work)
           cwavef1(:)=work(:)
           ABI_FREE(work)
           ABI_FREE(cwavef)
           ABI_FREE(gradug)
         end if

!$OMP PARALLEL DO
         do ir=1,nfftf
           rhor(ir,is)=rhor(ir,is) + CONJG(cwavef1(ir))*cwavef1(ir)*bks_weight
         end do
         !call cplx_addtorho(n1,n2,n3,n4,n5,n6,ndat,weight_r,ur,rho)

         if (Wfd%nspinor==2) then
           cwavef2 => wfr(1+nfftf:2*nfftf)
           wfr_x(:)=cwavef1(:)+cwavef2(:)       ! $(\Psi^{1}+\Psi^{2})$
           wfr_y(:)=cwavef1(:)-j_dpc*cwavef2(:) ! $(\Psi^{1}-i\Psi^{2})$
!$OMP PARALLEL DO
           do ir=1,nfftf
             rhor_down(ir)=rhor_down(ir)+CONJG(cwavef2(ir))*cwavef2(ir)*bks_weight
             rhor_mx  (ir)=rhor_mx  (ir)+CONJG(wfr_x  (ir))*wfr_x  (ir)*bks_weight
             rhor_my  (ir)=rhor_my  (ir)+CONJG(wfr_y  (ir))*wfr_y  (ir)*bks_weight
           end do
         end if

       end do
     end do
   end do

 end do ! enddo alpha

 if (myoptcalc==1) then ! convention for taur = 1/2 Sum_i |grad phi_i|^2
   rhor(:,:)=half*rhor(:,:)
 end if

 call iter_free(Iter_bks)

 if (wfd%nspinor == 2) then
   rhor(:, 2) = rhor_mx
   rhor(:, 3) = rhor_my
   rhor(:, 4) = rhor_down
 end if

 call xmpi_sum(rhor,Wfd%comm,ierr)

 ! === Symmetrization in G-space implementing also the AFM case ===
 n1=ngfftf(1)
 n2=ngfftf(2)
 n3=ngfftf(3)
 nfftotf=n1*n2*n3

 ABI_MALLOC(irrzon,(nfftotf**(1-1/Cryst%nsym),2,(Wfd%nspden/Wfd%nsppol)-3*(Wfd%nspden/4)))
 ABI_MALLOC(phnons,(2,nfftotf,(Wfd%nspden/Wfd%nsppol)-3*(Wfd%nspden/4)))

 if (Cryst%nsym/=1) then
   call irrzg(irrzon,Wfd%nspden,Wfd%nsppol,Cryst%nsym,n1,n2,n3,phnons,Cryst%symafm,Cryst%symrel,Cryst%tnons)
 end if

 cplex=1
 ABI_MALLOC(rhog,(2,cplex*nfftf))

 call symrhg(cplex,Cryst%gprimd,irrzon,Wfd%MPI_enreg,nfftf,nfftotf,ngfftf,Wfd%nspden,Wfd%nsppol,&
&  Cryst%nsym,Wfd%paral_kgb,phnons,rhog,rhor,Cryst%rprimd,Cryst%symafm,Cryst%symrel)

 ABI_FREE(rhog)
 ABI_FREE(phnons)
 ABI_FREE(irrzon)

 write(msg,'(a,f9.4)')' planewave contribution to nelect: ',SUM(rhor(:,1))*Cryst%ucvol/nfftf
 call wrtout(std_out,msg,'COLL')

 if (Wfd%nspden==4) then
   write(msg,'(a,3f9.4)')&
&     ' mx, my, mz: ',SUM(rhor(:,2))*Cryst%ucvol/nfftf,SUM(rhor(:,3))*Cryst%ucvol/nfftf,SUM(rhor(:,4))*Cryst%ucvol/nfftf
   call wrtout(std_out,msg,'COLL')
 end if

 ABI_FREE(wfr)

 if (Wfd%nspinor==2) then
   ABI_FREE(wfr_x)
   ABI_FREE(wfr_y)
   if (Wfd%nspden==4)  then
     ABI_FREE(rhor_down)
     ABI_FREE(rhor_mx)
     ABI_FREE(rhor_my)
   end if
 end if

 DBG_EXIT("COLL")

end subroutine wfd_mkrho
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/test_charge
!! NAME
!! test_charge
!!
!! FUNCTION
!!  Reports info on the electronic charge as well as Drude plasma frequency.
!!  Mainly used in the GW part.
!!
!! INPUTS
!!  nelectron_exp=Expected total number of electrons (used to normalize the charge)
!!
!! OUTPUT
!!
!! PARENTS
!!      bethe_salpeter,mrgscr,screening,sigma
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine test_charge(nfftf,nelectron_exp,nspden,rhor,ucvol,&
& usepaw,usexcnhat,usefinegrid,compch_sph,compch_fft,omegaplasma)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'test_charge'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftf,nspden,usefinegrid,usepaw,usexcnhat
 real(dp),intent(in) :: compch_fft,compch_sph,ucvol,nelectron_exp
 real(dp),intent(out) :: omegaplasma
!arrays
 real(dp),intent(inout) :: rhor(nfftf,nspden)

!Local variables ------------------------------
!scalars
 real(dp) :: nelectron_tot,nelectron_fft
 real(dp) :: nelectron_pw,nelectron_sph,rhoav,rs,nratio
 character(len=500) :: msg

!*************************************************************************

! ABI_UNUSED(usexcnhat)
if (usexcnhat==0)then
end if

 ! === For PAW output of compensation charges ===
 if (usepaw==1) then
!if (usepaw==1.and.usexcnhat>0) then ! TODO I still dont understand this if!
   write(msg,'(4a)')ch10,' PAW TEST:',ch10,' ==== Compensation charge inside spheres ============'
   if (compch_sph<greatest_real.and.compch_fft<greatest_real) &
&    write(msg,'(3a)')TRIM(msg),ch10,' The following values must be close...'
   if (compch_sph<greatest_real) &
&    write(msg,'(3a,f22.15)')TRIM(msg),ch10,' Compensation charge over spherical meshes = ',compch_sph
   if (compch_fft<greatest_real) then
     if (usefinegrid==1) then
       write(msg,'(3a,f22.15)')TRIM(msg),ch10,' Compensation charge over fine fft grid    = ',compch_fft
     else
       write(msg,'(3a,f22.15)')TRIM(msg),ch10,' Compensation charge over fft grid         = ',compch_fft
     end if
   end if
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,msg,'COLL')
   write(msg,'(a)')ch10
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,msg,'COLL')
 end if !PAW

 nelectron_pw =SUM(rhor(:,1))*ucvol/nfftf
 nelectron_tot=nelectron_pw
 nratio       =nelectron_exp/nelectron_tot

 if (usepaw==1) then
   nelectron_sph=nelectron_pw+compch_sph
   nelectron_fft=nelectron_pw+compch_fft
   nelectron_tot=nelectron_sph
   nratio=(nelectron_exp-nelectron_sph)/nelectron_pw
 end if

 rhoav=nelectron_tot/ucvol ; rs=(three/(four_pi*rhoav))**third
 if (usepaw==0) then
  write(msg,'(2(a,f9.4))')&
&   ' Number of electrons calculated from density = ',nelectron_tot,'; Expected = ',nelectron_exp
 else
   write(msg,'(2(a,f9.4),a)')&
&   ' Total number of electrons per unit cell = ',nelectron_sph,' (Spherical mesh), ',nelectron_fft,' (FFT mesh)'
 end if
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

!$write(msg,'(a,f9.4)')' Renormalizing smooth charge density using nratio = ',nratio
!! rhor(:,:)=nratio*rhor(:,:)

 write(msg,'(a,f9.6)')' average of density, n = ',rhoav
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,f9.4)')' r_s = ',rs
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')
 omegaplasma=SQRT(four_pi*rhoav)
 write(msg,'(a,f9.4,2a)')' omega_plasma = ',omegaplasma*Ha_eV,' [eV]',ch10
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

end subroutine test_charge
!!***
