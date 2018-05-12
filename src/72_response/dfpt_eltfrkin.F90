!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_eltfrkin
!! NAME
!! dfpt_eltfrkin
!!
!! FUNCTION
!! Compute the frozen-wavefunction kinetic enegy contribution to the
!! elastic tensor
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DRH, DCA, XG, GM, AR, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=<G|Cnk>=Fourier coefficients of wavefunction
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  ecutsm=smearing energy for plane wave kinetic energy (Ha) (NOT NEEDED !)
!!  effmass_free=effective mass for electrons (1. in common case)
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=work array for coordinates of G vectors in basis
!!  kptns(3,nkpt)=coordinates of k points in terms of reciprocal space
!!   primitive translations
!!  mband=maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mkmem=number of k points treated by this node.
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum dimension for number of planewaves
!!  nband(nkpt*nsppol)=number of bands being considered per k point
!!  nkpt=number of k points
!!  ngfft(18)=contain all needed information about 3D FFT, i
!!    see ~abinit/doc/variables/vargs.htm#ngfft
!!  npwarr(nkpt)=number of planewaves at each k point, and boundary
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for polarized
!!  occ(mband*nkpt*nsppol)=occupation numbers of bands (usually 2)
!!    at each k point
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  wtk(nkpt)=k point weights
!!
!! OUTPUT
!!  eltfrkin(6,6)=non-symmetrized kinetic energy contribution to the
!!                    elastic tensor
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      d2kindstr2,metric,sphereboundary,timab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfpt_eltfrkin(cg,eltfrkin,ecut,ecutsm,effmass_free,&
&  istwfk,kg,kptns,mband,mgfft,mkmem,mpi_enreg,&
&  mpw,nband,nkpt,ngfft,npwarr,nspinor,nsppol,occ,rprimd,wtk)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_xmpi
 use m_profiling_abi

 use m_time,         only : timab
 use m_geometry,     only : metric
 use m_fftcore,      only : sphereboundary

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_eltfrkin'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mgfft,mkmem,mpw,nkpt,nspinor,nsppol
 real(dp),intent(in) :: ecut,ecutsm,effmass_free
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: istwfk(nkpt),kg(3,mpw*mkmem),nband(nkpt*nsppol)
 integer,intent(in) :: ngfft(18),npwarr(nkpt)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol),kptns(3,nkpt)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),rprimd(3,3),wtk(nkpt)
 real(dp),intent(out) :: eltfrkin(6,6)

!Local variables-------------------------------
!scalars
 integer :: bdtot_index,iband,icg,ierr,ii,ikg
 integer :: ikpt,index,ipw,isppol,istwf_k,jj,master,me,n1,n2
 integer :: n3,nband_k,nkinout,npw_k,spaceComm
 real(dp) :: ucvol
!arrays
 integer,allocatable :: gbound(:,:),kg_k(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),kpoint(3),rmet(3,3),tsec(2)
 real(dp),allocatable :: cwavef(:,:),ekinout(:)
 real(dp),allocatable :: eltfrkink(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

!Default for sequential use
 master=0
!Init mpi_comm
 spaceComm=mpi_enreg%comm_cell
 me=mpi_enreg%me_kpt

!Compute gmet, gprimd and ucvol from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 eltfrkin(:,:)=0.0_dp
 bdtot_index=0
 icg=0

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 ABI_ALLOCATE(kg_k,(3,mpw))
 ABI_ALLOCATE(cwavef,(2,mpw*nspinor))
 ABI_ALLOCATE(eltfrkink,(6,6))

!Define k-points distribution

!LOOP OVER SPINS
 do isppol=1,nsppol
   ikg=0

!  Loop over k points
   do ikpt=1,nkpt

     nband_k=nband(ikpt+(isppol-1)*nkpt)
     istwf_k=istwfk(ikpt)
     npw_k=npwarr(ikpt)

!    Skip this k-point if not the proper processor
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) then
       bdtot_index=bdtot_index+nband_k
       cycle
     end if

     ABI_ALLOCATE(gbound,(2*mgfft+8,2))
     kpoint(:)=kptns(:,ikpt)

     kg_k(:,:) = 0

!$OMP PARALLEL DO PRIVATE(ipw) SHARED(ikg,kg,kg_k,npw_k)
     do ipw=1,npw_k
       kg_k(1,ipw)=kg(1,ipw+ikg)
       kg_k(2,ipw)=kg(2,ipw+ikg)
       kg_k(3,ipw)=kg(3,ipw+ikg)
     end do

     call sphereboundary(gbound,istwf_k,kg_k,mgfft,npw_k)

     index=1+icg

     eltfrkink(:,:)=0.0_dp

     nkinout=6*6
     ABI_ALLOCATE(ekinout,(nkinout))
     ekinout(:)=zero

     do iband=1,nband_k

       if(mpi_enreg%proc_distrb(ikpt,iband,isppol) /= me) cycle

       cwavef(:,1:npw_k*nspinor)=cg(:,1+(iband-1)*npw_k*nspinor+icg:iband*npw_k*nspinor+icg)

       call d2kindstr2(cwavef,ecut,ecutsm,effmass_free,ekinout,gmet,gprimd,&
&       istwf_k,kg_k,kpoint,npw_k,nspinor)

       eltfrkink(:,:)=eltfrkink(:,:)+ occ(iband+bdtot_index)* reshape(ekinout(:), (/6,6/) )

     end do !iband

     ABI_DEALLOCATE(ekinout)

     eltfrkin(:,:)=eltfrkin(:,:)+wtk(ikpt)*eltfrkink(:,:)

     ABI_DEALLOCATE(gbound)

     bdtot_index=bdtot_index+nband_k

     if (mkmem/=0) then
!      Handle case in which kg, cg, are kept in core
       icg=icg+npw_k*nspinor*nband_k
       ikg=ikg+npw_k
     end if

   end do
 end do  ! End loops on isppol and ikpt

!Fill in lower triangle
 do jj=2,6
   do ii=1,jj-1
     eltfrkin(jj,ii)=eltfrkin(ii,jj)
   end do
 end do

!Accumulate eltfrkin on all proc.
 call timab(48,1,tsec)
 call xmpi_sum(eltfrkin,spaceComm,ierr)
 call timab(48,2,tsec)

 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(eltfrkink)
 ABI_DEALLOCATE(kg_k)

 DBG_EXIT("COLL")

contains
!!***

!!****f* ABINIT/d2kindstr2
!! NAME
!! d2kindstr2
!!
!! FUNCTION
!! compute expectation value of the second derivatives of the kinetic energy
!! wrt strain for one band and kpoint
!!
!! INPUTS
!!  cwavef(2,npw*nspinor)=wavefunction for current band
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  ecutsm=smearing energy for plane wave kinetic energy (Ha)
!!  effmass_free=effective mass for electrons (1. in common case)
!!  gmet(3,3)=reciprocal lattice metric tensor ($\textrm{Bohr}^{-2}$)
!!  gprimd(3,3)=primitive vectors in reciprocal space
!!  istwfk=information about wavefunction storage
!!  kg_k(3,npw)=integer coordinates of planewaves in basis sphere.
!!  kpt(3)=reduced coordinates of k point
!!  npw=number of plane waves at kpt.
!!  nspinor=number of spinorial components of the wavefunction
!!
!! OUTPUT
!!  ekinout(36)=expectation values of the second strain derivatives
!!   of the (modified) kinetic energy
!!
!! NOTES
!! Usually, the kinetic energy expression is $(1/2) (2 \pi)^2 (k+G)^2 $
!! However, the present implementation allows for a modification
!! of this kinetic energy, in order to obtain smooth total energy
!! curves with respect to the cut-off energy or the cell size and shape.
!! Thus the usual expression is kept if it is lower then ecut-ecutsm,
!! zero is returned beyond ecut, and in between, the kinetic
!! energy is DIVIDED by a smearing factor (to make it infinite at the
!! cut-off energy). The smearing factor is $x^2 (3-2x)$, where
!! x = (ecut- unmodified energy)/ecutsm.
!!
!! PARENTS
!!      dfpt_eltfrkin
!!
!! CHILDREN
!!
!! SOURCE

subroutine d2kindstr2(cwavef,ecut,ecutsm,effmass_free,ekinout,gmet,gprimd,&
&            istwfk,kg_k,kpt,npw,nspinor)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'd2kindstr2'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwfk,npw,nspinor
 real(dp),intent(in) :: ecut,ecutsm,effmass_free
!arrays
 integer,intent(in) :: kg_k(3,npw)
 real(dp),intent(in) :: cwavef(2,npw*nspinor),gmet(3,3),gprimd(3,3),kpt(3)
 real(dp),intent(inout) :: ekinout(36) !vz_i

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: ig,igs,ii,ispinor,istr1,istr2,ka,kb,kd,kg
 real(dp) :: d2fkin,d2fsm,d2kinacc,d2kpg2,dfkin,dfsm,dkpg21,dkpg22,ecutsm_inv
 real(dp) :: fsm,gpk1,gpk2,gpk3,htpisq,kpg2,term,xx
!arrays
 integer,save :: idx(12)=(/1,1,2,2,3,3,3,2,3,1,2,1/)
 real(dp) :: d2gm(3,3),dgm01(3,3),dgm10(3,3)

! *************************************************************************
!
!htpisq is (1/2) (2 Pi) **2:
 htpisq=0.5_dp*(two_pi)**2

 ecutsm_inv=0.0_dp
 if(ecutsm>1.0d-20)ecutsm_inv=1/ecutsm

!Loop over 2nd strain index
 do istr2=1,6
!  Loop over 1st strain index, upper triangle only
   do istr1=1,istr2

     ka=idx(2*istr1-1);kb=idx(2*istr1);kg=idx(2*istr2-1);kd=idx(2*istr2)

     do ii = 1,3
       dgm01(:,ii)=-(gprimd(ka,:)*gprimd(kb,ii)+gprimd(kb,:)*gprimd(ka,ii))
       dgm10(:,ii)=-(gprimd(kg,:)*gprimd(kd,ii)+gprimd(kd,:)*gprimd(kg,ii))
     end do

     d2gm(:,:)=0._dp
     do ii = 1,3
       if(ka==kg) d2gm(:,ii)=d2gm(:,ii)&
&       +gprimd(kb,:)*gprimd(kd,ii)+gprimd(kd,:)*gprimd(kb,ii)
       if(ka==kd) d2gm(:,ii)=d2gm(:,ii)&
&       +gprimd(kb,:)*gprimd(kg,ii)+gprimd(kg,:)*gprimd(kb,ii)
       if(kb==kg) d2gm(:,ii)=d2gm(:,ii)&
&       +gprimd(ka,:)*gprimd(kd,ii)+gprimd(kd,:)*gprimd(ka,ii)
       if(kb==kd) d2gm(:,ii)=d2gm(:,ii)&
&       +gprimd(ka,:)*gprimd(kg,ii)+gprimd(kg,:)*gprimd(ka,ii)
     end do
     d2gm(:,:)=0.5_dp*d2gm(:,:)

     d2kinacc=0._dp

!    loop on spinor index
     do ispinor=1,nspinor
       igs=(ispinor-1)*npw
!      loop on plane waves
       do ig=1,npw
         gpk1=dble(kg_k(1,ig))+kpt(1)
         gpk2=dble(kg_k(2,ig))+kpt(2)
         gpk3=dble(kg_k(3,ig))+kpt(3)
         kpg2=htpisq*&
&         ( gmet(1,1)*gpk1**2+         &
&         gmet(2,2)*gpk2**2+         &
&         gmet(3,3)*gpk3**2          &
&         +2.0_dp*(gpk1*gmet(1,2)*gpk2+  &
&         gpk1*gmet(1,3)*gpk3+  &
&         gpk2*gmet(2,3)*gpk3 )  )
         dkpg21=htpisq*&
&         ( dgm01(1,1)*gpk1**2+         &
&         dgm01(2,2)*gpk2**2+         &
&         dgm01(3,3)*gpk3**2          &
&         +2.0_dp*(gpk1*dgm01(1,2)*gpk2+  &
&         gpk1*dgm01(1,3)*gpk3+  &
&         gpk2*dgm01(2,3)*gpk3 )  )
         dkpg22=htpisq*&
&         ( dgm10(1,1)*gpk1**2+         &
&         dgm10(2,2)*gpk2**2+         &
&         dgm10(3,3)*gpk3**2          &
&         +2.0_dp*(gpk1*dgm10(1,2)*gpk2+  &
&         gpk1*dgm10(1,3)*gpk3+  &
&         gpk2*dgm10(2,3)*gpk3 )  )
         d2kpg2=htpisq*&
&         ( d2gm(1,1)*gpk1**2+         &
&         d2gm(2,2)*gpk2**2+         &
&         d2gm(3,3)*gpk3**2          &
&         +2.0_dp*(gpk1*d2gm(1,2)*gpk2+  &
&         gpk1*d2gm(1,3)*gpk3+  &
&         gpk2*d2gm(2,3)*gpk3 )  )

         if(kpg2>ecut-tol12)then
           dfkin=0._dp
           d2fkin=0._dp
         elseif(kpg2>ecut-ecutsm)then
!          This kinetic cutoff smoothing function and its xx derivatives
!          were produced with Mathematica and the fortran code has been
!          numerically checked against Mathematica.
           xx=(ecut-kpg2)*ecutsm_inv
           fsm=1.0_dp/(xx**2*(3+xx*(1+xx*(-6+3*xx))))
           dfsm=-3.0_dp*(-1+xx)**2*xx*(2+5*xx)*fsm**2
           d2fsm=6.0_dp*xx**2*(9+xx*(8+xx*(-52+xx*(-3+xx*(137+xx*&
&           (-144+45*xx))))))*fsm**3
           dfkin=fsm-ecutsm_inv*kpg2*dfsm
           d2fkin=ecutsm_inv*(-2.0_dp*dfsm+ecutsm_inv*kpg2*d2fsm)
         else
           dfkin=1._dp
           d2fkin=0._dp
         end if

!        accumulate kinetic energy 2nd derivative with wavefunction components
         term=d2fkin*dkpg21*dkpg22 + dfkin*d2kpg2
         if(istwfk==2 .and. ig/=1)term=2.0_dp*term
         if(istwfk>2)term=2.0_dp*term
         d2kinacc=d2kinacc + term*(cwavef(re,ig+igs)**2 + cwavef(im,ig+igs)**2)

       end do  !ig
     end do !ispinor

     ekinout(istr1+6*(istr2-1))=d2kinacc/effmass_free

   end do !istr1
 end do !istr2

end subroutine d2kindstr2
!!***

end subroutine dfpt_eltfrkin
!!***
