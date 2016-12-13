!{\src2tex{textfont=tt}}
!!****f* ABINIT/fxphas_seq
!!
!! NAME
!! fxphas_seq
!!
!! FUNCTION
!! Fix phase of all bands. Keep normalization but maximize real part
!! (minimize imag part). Also fix the sign of real part
!! by setting the first non-zero element to be positive.
!!
!! This version has been stripped of all the mpi_enreg junk by MJV
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR, MT, MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mcg)= contains the wavefunction |c> coefficients.
!!  gsc(2,mgsc)= if useoverlap==1, contains the S|c> coefficients,
!!               where S is an overlap matrix.
!!  icg=shift to be applied on the location of data in the array cg
!!  igsc=shift to be applied on the location of data in the array gsc
!!  istwfk=input option parameter that describes the storage of wfs
!!    (set to 1 if usual complex vectors)
!!  mcg=size of second dimension of cg
!!  mgsc=size of second dimension of gsc
!!  nband_k=number of bands
!!  npw_k=number of planewaves
!!  useoverlap=describe the overlap of wavefunctions:
!!               0: no overlap (S=Identi0,ty_matrix)
!!               1: wavefunctions are overlapping
!!
!! OUTPUT
!!  cg(2,mcg)=same array with altered phase.
!!  gsc(2,mgsc)= same array with altered phase.
!!
!! NOTES
!! When the sign of the real part was fixed (modif v3.1.3g.6), the
!! test Tv3#5 , dataset 5, behaved differently than previously.
!! This should be cleared up.
!!
!! PARENTS
!!      dfpt_phfrq,rayleigh_ritz
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine fxphas_seq(cg,gsc,icg,igsc,istwfk,mcg,mgsc,nband_k,npw_k,useoverlap)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fxphas_seq'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icg,igsc,istwfk,mcg,mgsc,nband_k,npw_k,useoverlap
!arrays
 real(dp),intent(inout) :: cg(2,mcg),gsc(2,mgsc*useoverlap)

!Local variables-------------------------------
!scalars
 integer :: iband,ii,indx
 real(dp) :: cim,cre,gscim,gscre,quotient,root1,root2,saa,sab,sbb,theta
 real(dp) :: thppi,xx,yy
 character(len=500) :: message
!arrays
 real(dp),allocatable :: cimb(:),creb(:),saab(:),sabb(:),sbbb(:) !,sarr(:,:)

! *************************************************************************

!The general case, where a complex phase indeterminacy is present
 if(istwfk==1)then

   ABI_ALLOCATE(cimb,(nband_k))
   ABI_ALLOCATE(creb,(nband_k))
   ABI_ALLOCATE(saab,(nband_k))
   ABI_ALLOCATE(sabb,(nband_k))
   ABI_ALLOCATE(sbbb,(nband_k))
   cimb(:)=zero ; creb(:)=zero

!  Loop over bands
!  TODO: MG store saa arrays in sarr(3,nband_k) to reduce false sharing.
   do iband=1,nband_k
     indx=icg+(iband-1)*npw_k

!    Compute several sums over Re, Im parts of c
     saa=0.0_dp ; sbb=0.0_dp ; sab=0.0_dp
     do ii=1+indx,npw_k+indx
       saa=saa+cg(1,ii)*cg(1,ii)
       sbb=sbb+cg(2,ii)*cg(2,ii)
       sab=sab+cg(1,ii)*cg(2,ii)
     end do
     saab(iband)=saa
     sbbb(iband)=sbb
     sabb(iband)=sab
   end do ! iband


   do iband=1,nband_k

     indx=icg+(iband-1)*npw_k

     saa=saab(iband)
     sbb=sbbb(iband)
     sab=sabb(iband)

!    Get phase angle theta
     if (sbb+saa>tol8)then
       if(abs(sbb-saa)>tol8*(sbb+saa) .or. 2*abs(sab)>tol8*(sbb+saa))then
         if (abs(sbb-saa)>tol8*abs(sab)) then
           quotient=sab/(sbb-saa)
           theta=0.5_dp*atan(2.0_dp*quotient)
         else
!          Taylor expansion of the atan in terms of inverse of its argument. Correct up to 1/x2, included.
           theta=0.25_dp*(pi-(sbb-saa)/sab)
         end if
!        Check roots to get theta for max Re part
         root1=cos(theta)**2*saa+sin(theta)**2*sbb-2.0_dp*cos(theta)*sin(theta)*sab
         thppi=theta+0.5_dp*pi
         root2=cos(thppi)**2*saa+sin(thppi)**2*sbb-2.0_dp*cos(thppi)*sin(thppi)*sab
         if (root2>root1) theta=thppi
       else
!        The real part vector and the imaginary part vector are orthogonal, and of same norm. Strong indeterminacy.
!        Will determine the first non-zero coefficient, and fix its phase 
         do ii=1+indx,npw_k+indx
           cre=cg(1,ii)
           cim=cg(2,ii)
           if(cre**2+cim**2>tol8**2*(saa+sbb))then
             if(cre**2>tol8**2**cim**2)then
               theta=atan(cim/cre)
             else
!              Taylor expansion of the atan in terms of inverse of its argument. Correct up to 1/x2, included.
               theta=pi/2-cre/cim
             end if
             exit
           end if
         end do
       end if
     else
       write(message,'(a,i0,5a)')&
&       'The eigenvector with band ',iband,' has zero norm.',ch10,&
&       'This usually happens when the number of bands (nband) is comparable to the number of planewaves (mpw)',ch10,&
&       'Action: Check the parameters of the calculation. If nband ~ mpw, then decrease nband or, alternatively, increase ecut'
       MSG_ERROR(message)
     end if

     xx=cos(theta)
     yy=sin(theta)

!    Here, set the first non-zero element to be positive
     do ii=1+indx,npw_k+indx
       cre=cg(1,ii)
       cim=cg(2,ii)
       cre=xx*cre-yy*cim
       if(abs(cre)>tol8)exit
     end do
     if(cre<zero)then
       xx=-xx ; yy=-yy
     end if

     creb(iband)=xx
     cimb(iband)=yy

   end do

   do iband=1,nband_k

     indx=icg+(iband-1)*npw_k

     xx=creb(iband)
     yy=cimb(iband)
     do ii=1+indx,npw_k+indx
       cre=cg(1,ii)
       cim=cg(2,ii)
       cg(1,ii)=xx*cre-yy*cim
       cg(2,ii)=xx*cim+yy*cre
     end do

!    Alter phase of array S|cg>
     if (useoverlap==1) then
       indx=igsc+(iband-1)*npw_k
       do ii=1+indx,npw_k+indx
         gscre=gsc(1,ii)
         gscim=gsc(2,ii)
         gsc(1,ii)=xx*gscre-yy*gscim
         gsc(2,ii)=xx*gscim+yy*gscre
       end do
     end if

   end do ! iband

   ABI_DEALLOCATE(cimb)
   ABI_DEALLOCATE(creb)
   ABI_DEALLOCATE(saab)
   ABI_DEALLOCATE(sabb)
   ABI_DEALLOCATE(sbbb)

!  ====================================================================

!  Storages that take into account the time-reversal symmetry : the freedom is only a sign freedom
 else  ! if istwfk/=1

   ABI_ALLOCATE(creb,(nband_k))
   creb(:)=zero
!  Loop over bands
   do iband=1,nband_k

     indx=icg+(iband-1)*npw_k

!    Here, set the first non-zero real element to be positive
     do ii=1+indx,npw_k+indx
       cre=cg(1,ii)
       if(abs(cre)>tol8)exit
     end do
     creb(iband)=cre

   end do ! iband

   do iband=1,nband_k

     cre=creb(iband)
     if(cre<zero)then
       indx=icg+(iband-1)*npw_k
       do ii=1+indx,npw_k+indx
         cg(1,ii)=-cg(1,ii)
         cg(2,ii)=-cg(2,ii)
       end do
       if(useoverlap==1)then
         indx=igsc+(iband-1)*npw_k
         do ii=1+indx,npw_k+indx
           gsc(1,ii)=-gsc(1,ii)
           gsc(2,ii)=-gsc(2,ii)
         end do
       end if
     end if

   end do ! iband

   ABI_DEALLOCATE(creb)

 end if ! istwfk

end subroutine fxphas_seq
!!***
