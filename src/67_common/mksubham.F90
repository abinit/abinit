!{\src2tex{textfont=tt}}
!!****f* ABINIT/mksubham
!! NAME
!! mksubham
!!
!! FUNCTION
!! Build the Hamiltonian matrix in the eigenfunctions subspace,
!! for one given band (or for one given block of bands)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cg(2,mcg)=wavefunctions
!!  gsc(2,mgsc)=<g|S|c> matrix elements (S=overlap)
!!  iblock=index of block of bands
!!  icg=shift to be applied on the location of data in the array cg
!!  igsc=shift to be applied on the location of data in the array cg
!!  istwf_k=input parameter that describes the storage of wfs
!!  mcg=second dimension of the cg array
!!  mgsc=second dimension of the gsc array
!!  nband_k=number of bands at this k point for that spin polarization
!!  nbdblock=number of bands in a block
!!  npw_k=number of plane waves at this k point
!!  nspinor=number of spinorial components of the wavefunctions
!!  use_subovl=1 if the overlap matrix is not identity in WFs subspace
!!  use_vnl= 1 if <C band,k|H|C band_prime,k> has to be computed
!!  me_g0=1 if this processors has G=0, 0 otherwise
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  ghc(2,npw_k*nspinor)=<G|H|C band,k> for the current state
!!                       This is an input in non-blocked algorithm
!!                               an output in blocked algorithm
!!  gvnlc(2,npw_k*nspinor)=<G|Vnl|C band,k> for the current state
!!                       This is an input in non-blocked algorithm
!!                               an output in blocked algorithm
!!  isubh=index of current state in array subham
!!  isubo=index of current state in array subovl
!!  subham(nband_k*(nband_k+1))=Hamiltonian expressed in the WFs subspace
!!  subovl(nband_k*(nband_k+1)*use_subovl)=overlap matrix expressed in the WFs subspace
!!  subvnl(nband_k*(nband_k+1)*use_vnl)=non-local Hamiltonian expressed in the WFs subspace
!!
!! PARENTS
!!      cgwf
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mksubham(cg,ghc,gsc,gvnlc,iblock,icg,igsc,istwf_k,&
&                    isubh,isubo,mcg,mgsc,nband_k,nbdblock,npw_k,&
&                    nspinor,subham,subovl,subvnl,use_subovl,use_vnl,me_g0)

 use defs_basis
 use m_profiling_abi
 use m_cgtools

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mksubham'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iblock,icg,igsc,istwf_k,mcg,mgsc,nband_k
 integer,intent(in) :: nbdblock,npw_k,nspinor,use_subovl,use_vnl,me_g0
 integer,intent(inout) :: isubh,isubo
!arrays
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(in) :: gsc(2,mgsc)
 real(dp),intent(inout) :: ghc(2,npw_k*nspinor),gvnlc(2,npw_k*nspinor)
 real(dp),intent(inout) :: subham(nband_k*(nband_k+1))
 real(dp),intent(inout) :: subovl(nband_k*(nband_k+1)*use_subovl)
 real(dp),intent(inout) :: subvnl(nband_k*(nband_k+1)*use_vnl)

!Local variables-------------------------------
!scalars
 integer :: iband,ibdblock,ii,ipw,ipw1,isp,iwavef,jwavef
 real(dp) :: cgimipw,cgreipw,chcim,chcre,cscim,cscre,cvcim,cvcre
!real(dp) :: chc(2),cvc(2),csc(2)

! *********************************************************************

!Loop over bands in a block This loop can be parallelized
 do iband=1+(iblock-1)*nbdblock,min(iblock*nbdblock,nband_k)
   ibdblock=iband-(iblock-1)*nbdblock

!  Compute elements of subspace Hamiltonian <C(i)|H|C(n)> and <C(i)|Vnl|C(n)>
   if(istwf_k==1)then

     do ii=1,iband
       iwavef=(ii-1)*npw_k*nspinor+icg
       chcre=zero ; chcim=zero
       if (use_vnl==0) then
         do ipw=1,npw_k*nspinor
           cgreipw=cg(1,ipw+iwavef)
           cgimipw=cg(2,ipw+iwavef)
           chcre=chcre+cgreipw*ghc(1,ipw)+cgimipw*ghc(2,ipw)
           chcim=chcim+cgreipw*ghc(2,ipw)-cgimipw*ghc(1,ipw)
         end do
!        chc = cg_zdotc(npw_k*nspinor,cg(1,1+iwavef),ghc)
       else
#if 1
         do ipw=1,npw_k*nspinor
           cgreipw=cg(1,ipw+iwavef)
           cgimipw=cg(2,ipw+iwavef)
           chcre=chcre+cgreipw*ghc(1,ipw)+cgimipw*ghc(2,ipw)
           chcim=chcim+cgreipw*ghc(2,ipw)-cgimipw*ghc(1,ipw)
         end do
         cvcre=zero ; cvcim=zero
         do ipw=1,npw_k*nspinor
           cgreipw=cg(1,ipw+iwavef)
           cgimipw=cg(2,ipw+iwavef)
           cvcre=cvcre+cgreipw*gvnlc(1,ipw)+cgimipw*gvnlc(2,ipw)
           cvcim=cvcim+cgreipw*gvnlc(2,ipw)-cgimipw*gvnlc(1,ipw)
         end do
         subvnl(isubh  )=cvcre
         subvnl(isubh+1)=cvcim
#else
!        New version with BLAS1, will require some update of the refs.
         cvc = cg_zdotc(npw_k*nspinor,cg(1,1+iwavef),gvnlc)
         subvnl(isubh  )=cvc(1)
         subvnl(isubh+1)=cvc(2)
         chc = cg_zdotc(npw_k*nspinor,cg(1,1+iwavef),ghc)
         chcre = chc(1)
         chcim = chc(2)
#endif
!        Store real and imag parts in Hermitian storage mode:
       end if
       subham(isubh  )=chcre
       subham(isubh+1)=chcim
!      subham(isubh  )=chc(1)
!      subham(isubh+1)=chc(2)
       isubh=isubh+2
     end do

   else if(istwf_k>=2)then
     do ii=1,iband
       iwavef=(ii-1)*npw_k+icg
!      Use the time-reversal symmetry, but should not double-count G=0
       if(istwf_k==2 .and. me_g0==1) then
         chcre = half*cg(1,1+iwavef)*ghc(1,1)
         if (use_vnl==1) cvcre=half*cg(1,1+iwavef)*gvnlc(1,1)
         ipw1=2
       else
         chcre=zero; ipw1=1
         if (use_vnl==1) cvcre=zero
       end if
       if (use_vnl==0) then
         do isp=1,nspinor
           do ipw=ipw1+(isp-1)*npw_k,npw_k*isp
             cgreipw=cg(1,ipw+iwavef)
             cgimipw=cg(2,ipw+iwavef)
             chcre=chcre+cgreipw*ghc(1,ipw)+cgimipw*ghc(2,ipw)
           end do
         end do
         chcre=two*chcre
       else
         do isp=1,nspinor
           do ipw=ipw1+(isp-1)*npw_k,npw_k*isp
             cgreipw=cg(1,ipw+iwavef)
             cgimipw=cg(2,ipw+iwavef)
             chcre=chcre+cgreipw*ghc(1,ipw)+cgimipw*ghc(2,ipw)
             cvcre=cvcre+cgreipw*gvnlc(1,ipw)+cgimipw*gvnlc(2,ipw)
           end do
         end do
         chcre=two*chcre
         cvcre=two*cvcre
!        Store real and imag parts in Hermitian storage mode:
         subvnl(isubh  )=cvcre
         subvnl(isubh+1)=zero
       end if
       subham(isubh  )=chcre
       subham(isubh+1)=zero
       isubh=isubh+2
     end do
   end if

!  Compute elements of subspace <C(i)|S|C(n)> (S=overlap matrix)
!  <C(i)|S|C(n)> should be closed to Identity.
   if (use_subovl==1) then
     jwavef=(iband-1)*npw_k*nspinor+igsc
     if(istwf_k==1)then
       do ii=1,iband
         iwavef=(ii-1)*npw_k*nspinor+icg
         cscre=zero ; cscim=zero
         do ipw=1,npw_k*nspinor
           cgreipw=cg(1,ipw+iwavef)
           cgimipw=cg(2,ipw+iwavef)
           cscre=cscre+cgreipw*gsc(1,ipw+jwavef)+cgimipw*gsc(2,ipw+jwavef)
           cscim=cscim+cgreipw*gsc(2,ipw+jwavef)-cgimipw*gsc(1,ipw+jwavef)
         end do
!        csc = cg_zdotc(npw_k*nspinor,cg(1,1+iwavef),gsc)
!        subovl(isubo  )=csc(1)
!        subovl(isubo+1)=csc(2)
!        Store real and imag parts in Hermitian storage mode:
         subovl(isubo  )=cscre
         subovl(isubo+1)=cscim
         isubo=isubo+2
       end do
     else if(istwf_k>=2)then
       do ii=1,iband
         iwavef=(ii-1)*npw_k*nspinor+icg
         if(istwf_k==2 .and. me_g0==1)then
           cscre=half*cg(1,1+iwavef)*gsc(1,1+jwavef)
           ipw1=2
         else
           cscre=zero; ipw1=1
         end if
         do isp=1,nspinor
           do ipw=ipw1+(isp-1)*npw_k,npw_k*isp
             cgreipw=cg(1,ipw+iwavef)
             cgimipw=cg(2,ipw+iwavef)
             cscre=cscre+cg(1,ipw+iwavef)*gsc(1,ipw+jwavef)+cg(2,ipw+iwavef)*gsc(2,ipw+jwavef)
           end do
         end do
         cscre=two*cscre
!        Store real and imag parts in Hermitian storage mode:
         subovl(isubo  )=cscre
         subovl(isubo+1)=zero
         isubo=isubo+2
       end do
     end if
   end if

 end do ! iband in a block

end subroutine mksubham
!!***
