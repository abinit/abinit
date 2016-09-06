!{\src2tex{textfont=tt}}
!!****f* ABINIT/smatrix
!! NAME
!! smatrix
!!
!! FUNCTION
!! Compute the overlap matrix between the k-points k and k + dk.
!! Depending on the value of job and ddkflag, compute also its determinant,
!! its inverse and the product of its inverse with the wavefunctions at k.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2016 ABINIT  group (MVeithen)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cg(2,mcg_k) = planewave coefficients of wavefunctions at k
!! cgq(2,mcg_q) = planewave coefficients of wavefunctions at q = k + dk
!! ddkflag = 1 : compute product of the inverse overlap matrix
!!               with the wavefunction at k (job = 1 or 11)
!!           0 : do not compute the product of the inverse overlap matrix
!!               with the wavefunction at k
!! icg = shift applied to the wavefunctions in the array cg
!! icg1 = shift applied to the wavefunctions in the array cgq
!! itrs = variable that governs the use of time-reversal symmetry
!!        when converting the wavefunctions from the iBZ to the fBZ
!! job = type of calculation
!!        0 : update overlap matrix only
!!        1 : like 0 but also compute inverse of the overlap matrix
!!       10 : like 0 but also compute determinant of the overlap matrix
!!       11 : like 0 but also compute determinant and inverse of the overlap matrix
!!       20 : like 0 but also transfer cgq to cg1_k without multiplication by S^{-1}
!!       21 : like 1 but also transfer cgq to cg1_k without multiplication by S^{-1}
!! maxbd = used in case ddkflag = 1, defines the highest band for
!!         which the ddk will be computed
!! mcg_k = second dimension of cg
!! mcg_q = second dimension of cg_q
!! mcg1_k = second dimension of cg1_k, should be equal to
!!          mpw*nsppol*nspinor*(maxbd - minbd + 1)
!! minbd = used in case ddkflag = 1, defines the lowest band for
!!         which the ddk will be computed
!! mpw = maximum dimensioned size of npw
!! mband_occ = max number of occupied valence bands for both spins
!! nband_occ = number of (occupied) valence bands
!! npw_k1 = number of plane waves at k
!! npw_k2 = number of plane waves at k + dk
!! nspinor = number of spinorial components of the wavefunctions
!! nsppol = 1 for unpolarized, 2 for spin-polarized
!! pwind_k = array used to compute the overlap matrix
!! pwnsfac = phase factors for non-symmorphic translations
!! shifrbd = shift applied to the location of the WF in the cg-array
!!           after each loop over bands
!!           0 : apply no shift, this is allowed in case cg
!!               contains only the wf of one single band
!!           1 : apply a shift of npw_k1*nspinor, this is the usual option
!!               when cg contains the wf for all occupied bands
!! smat_k_paw : overlap matrix due to on-site PAW terms between different bands
!!              at k and k+b. Only relevant when usepaw = 1, and is to be computed
!!              previously by smatrix_k_paw.F90
!! usepaw = flag governing use of PAW: 0 means no PAW, 1 means PAW is used
!!
!! OUTPUT
!! cg1_k(2,mcg1_k) = product of the inverse overlap matrix with the
!!                   wavefunctions at k; computed in case job = 1 or 11;
!!                   or just cgq in case of job = 20 or 21
!! dtm_k(2) = determinant of the overlap matrix between k and k + dk;
!!            computed in case job = 10 or 11
!! smat_inv = inverse of the overlap matrix
!!
!! SIDE EFFECTS
!! Input/Output
!! sflag_k(iband) = 1 if the elements smat_k(:,iband,:) are up to date
!!                    -> they will not be recomputed
!!                  0 the elements smat_k(:,iband,:) will be recomputed
!!      at the end of the routine, sflag_k(1:mband_occ) = 1
!!      (the whole overlap matrix is up to date)
!! smat_k = overlap matrix between k, k + dk
!!          only the lines for which sflag_k = 0 are computed
!!          smat_k(:,n,m) = < u_{n,k} | u_{m,k+dk} >
!!
!! NOTES
!! This routine is quite flexible in the way it deals with the wavefunctions:
!!  - cg (WF at k) can contain either the whole WF array (all k-points
!!    and bands), in which case the location of the WF at k is specified
!!    by icg and shiftbd = 1, or the WF of a single k-point/band, in which case
!!    shiftbd = 0 and icg = 0.
!!  - cgq (WF at k + dk) can contain either the whole WF array (all k-points
!!    and bands), in which case the location of the WF at k is specified
!!    by icg1, or the WF of a single k-point, in which case
!!    icg1 = 0. cgq must contain the WF of ALL occupied bands.
!!  - cg1_k can either be computed for all valence bands or
!!    for a group of valence bands defined by minbd and maxbd.
!!
!! PARENTS
!!      berryphase_new,cgwf,getcgqphase,make_grad_berry
!!
!! CHILDREN
!!      dzgedi,dzgefa,overlap_g
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine smatrix(cg,cgq,cg1_k,ddkflag,dtm_k,icg,icg1,itrs,job,maxbd,&
&  mcg_k,mcg_q,mcg1_k,minbd,mpw,mband_occ,nband_occ,npw_k1,npw_k2,nspinor,&
&  pwind_k,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_k,smat_k_paw,usepaw)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'smatrix'
 use interfaces_28_numeric_noabirule
 use interfaces_32_util, except_this_one => smatrix
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ddkflag,icg,icg1,itrs,job,maxbd,mcg1_k,mcg_k,mcg_q
 integer,intent(in) :: minbd,mpw,mband_occ,npw_k1,npw_k2,nspinor,shiftbd
 integer,intent(in) :: nband_occ
 integer,intent(in) :: usepaw
!arrays
 integer,intent(in) :: pwind_k(mpw)
 integer,intent(inout) :: sflag_k(mband_occ)
 real(dp),intent(in) :: cg(2,mcg_k),cgq(2,mcg_q),pwnsfac_k(4,mpw)
 real(dp),intent(in) :: smat_k_paw(2,usepaw*mband_occ,usepaw*mband_occ)
 real(dp),intent(inout) :: smat_k(2,mband_occ,mband_occ)
 real(dp),intent(out) :: cg1_k(2,mcg1_k),dtm_k(2)
 real(dp),intent(out) :: smat_inv(2,mband_occ,mband_occ)

!Local variables -------------------------
!scalars
 integer :: count,dzgedi_job,iband,info,ipw,ispinor,jband,jband1,jpw,pwmax,pwmin
 integer :: spnshft_k1, spnshft_k2
 real(dp) :: doti,dotr,fac,wfi,wfr
 character(len=500) :: message
!arrays
 integer,allocatable :: ipvt(:)
! integer,allocatable :: my_pwind_k(:) ! used in debugging below
 real(dp) :: det(2,2)
 real(dp),allocatable :: vect1(:,:),vect2(:,:),zgwork(:,:)

! ***********************************************************************

!DEBUG
!write(std_out,*)'smatrix : enter'
!write(std_out,*)'sflag_k = ',sflag_k
!write(std_out,'(a,4i4)' )'job, ddkflag, shiftbd, itrs = ',job,ddkflag,shiftbd,itrs
!write(std_out,'(a,2i6)')' JWZ smatrix.F90 debug : npw_k1, npw_k2 ',npw_k1,npw_k2
!stop
!ENDDEBUG

 ABI_ALLOCATE(ipvt,(nband_occ))
 ABI_ALLOCATE(zgwork,(2,nband_occ))
 ABI_ALLOCATE(vect1,(2,0:mpw*nspinor))
 ABI_ALLOCATE(vect2,(2,0:mpw*nspinor))
 vect1(:,0) = zero ; vect2(:,0) = zero

!Check if the values of ddkflag and job are compatible

 if ((job /= 0).and.(job /= 1).and.(job /= 10).and.(job /= 11).and.(job/=20).and.(job/=21)) then
   write(message,'(a,i3,a,a)')&
&   ' job is equal to ',job,ch10,&
&   ' while only the values job = 0, 1, 10, 11, 20, or 21 are allowed.'
   MSG_ERROR(message)
 end if

 if (ddkflag == 1) then
   if ((job/=1).and.(job/=11)) then
     write(message,'(a,i0,a,a)')&
&     ' job is equal to ',job,ch10,&
&     ' while ddkflag = 1. This is not allowed.'
     MSG_ERROR(message)
   end if
 end if

!Check the values of sflag_k
 do iband=1,nband_occ
   if (sflag_k(iband)/=0 .and. sflag_k(iband)/=1)then
     write(message,'(3a,i4,a,i4)')&
&     '  The content of sflag_k must be 0 or 1.',ch10,&
&     '  However, for iband=',iband,', sflag_k(iband)=',sflag_k(iband)
     MSG_ERROR(message)
   end if
 end do

!Check if shiftbd is consistent with sflag_k
 if (shiftbd == 0) then
   count = 0
   do iband = 1, nband_occ
     if (sflag_k(iband) == 0) count = count + 1
   end do
   if (count > 1) then
     message = 'in case shiftbd = 0, only 1 element of sflag can be 0'
     MSG_ERROR(message)
   end if
 end if

!Update the lines of the overlap matrix for which sflag = 0
!MVeithen: because of sflag, it is more efficient to perform
!the loop over jband inside the loop over iband

!DEBUG
!write(std_out,*)' smatrix : smat_k(1,1,1)=',smat_k(1,1,1)
!write(std_out,*)' smatrix : sflag_k=',sflag_k
!ENDDEBUG

!!
!! debugging based on norm of k1 vector
!
!ABI_ALLOCATE(my_pwind_k,(mpw))
!!
!do iband = 1, nband_occ
!!
!pwmin = (iband-1)*npw_k1*nspinor*shiftbd
!pwmax = pwmin + npw_k1*nspinor
!!
!!!    Multiply the bra wave function by the phase factor
!if (itrs==1.or.itrs==11) then  ! take complex conjugate of bra
!do ispinor = 1, nspinor
!spnshft_k1=(ispinor-1)*npw_k1
!do ipw = 1,npw_k1
!vect1(1,spnshft_k1+ipw) = cg(1,icg+spnshft_k1+ipw+pwmin)*pwnsfac_k(1,ipw) &
!&           +cg(2,icg+spnshft_k1+ipw+pwmin)*pwnsfac_k(2,ipw)
!vect1(2,spnshft_k1+ipw) = cg(1,icg+spnshft_k1+ipw+pwmin)*pwnsfac_k(2,ipw) &
!&           -cg(2,icg+spnshft_k1+ipw+pwmin)*pwnsfac_k(1,ipw)
!end do
!end do
!else
!do ispinor=1, nspinor
!spnshft_k1=(ispinor-1)*npw_k1
!do ipw = 1,npw_k1
!vect1(1,spnshft_k1+ipw) = cg(1,icg+spnshft_k1+ipw+pwmin)*pwnsfac_k(1,ipw) &
!&           -cg(2,icg+spnshft_k1+ipw+pwmin)*pwnsfac_k(2,ipw)
!vect1(2,spnshft_k1+ipw) = cg(1,icg+spnshft_k1+ipw+pwmin)*pwnsfac_k(2,ipw) &
!&           +cg(2,icg+spnshft_k1+ipw+pwmin)*pwnsfac_k(1,ipw)
!end do
!end do
!end if
!
!!    
!if (npw_k1*nspinor < mpw*nspinor) vect1(:,npw_k1*nspinor+1:mpw*nspinor) = zero
!
!do jband = 1, nband_occ
!
!pwmin = (jband-1)*npw_k1*nspinor
!pwmax = pwmin + npw_k1*nspinor
!
!if (itrs==10.or.itrs==11) then ! take complex conjugate of ket
!do ispinor=1, nspinor
!spnshft_k2=(ispinor-1)*npw_k1
!do ipw = 1, npw_k1
!my_pwind_k(ipw) = ipw
!vect2(1,spnshft_k2+ipw) = cg(1,icg+spnshft_k2+ipw+pwmin)*pwnsfac_k(1,ipw) &
!&             +cg(2,icg+spnshft_k2+ipw+pwmin)*pwnsfac_k(2,ipw)
!vect2(2,spnshft_k2+ipw) = cg(1,icg+spnshft_k2+ipw+pwmin)*pwnsfac_k(2,ipw) &
!&             -cg(2,icg+spnshft_k2+ipw+pwmin)*pwnsfac_k(1,ipw)
!end do
!end do
!else
!do ispinor=1, nspinor
!spnshft_k2=(ispinor-1)*npw_k1
!do ipw = 1, npw_k1
!my_pwind_k(ipw) = ipw
!vect2(1,spnshft_k2+ipw) = cg(1,icg+spnshft_k2+ipw+pwmin)*pwnsfac_k(1,ipw) &
!&             -cg(2,icg+spnshft_k2+ipw+pwmin)*pwnsfac_k(2,ipw)
!vect2(2,spnshft_k2+ipw) = cg(1,icg+spnshft_k2+ipw+pwmin)*pwnsfac_k(2,ipw) &
!&             +cg(2,icg+spnshft_k2+ipw+pwmin)*pwnsfac_k(1,ipw)
!end do
!end do
!end if
!!
!if (npw_k1*nspinor < mpw*nspinor) vect2(:,npw_k1*nspinor+1:mpw*nspinor) = zero
!!
!call overlap_g(doti,dotr,mpw,npw_k1,npw_k1,nspinor,my_pwind_k,vect1,vect2)
!!
!smat_k(1,iband,jband) = dotr
!smat_k(2,iband,jband) = doti
!!
!if (usepaw == 1) then
!smat_k(1,iband,jband) = smat_k(1,iband,jband)+smat_k_paw(1,iband,jband)
!smat_k(2,iband,jband) = smat_k(2,iband,jband)+smat_k_paw(2,iband,jband)
!end if
!!
!if( (smat_k(1,iband,jband)**2+smat_k(2,iband,jband)**2)>tol8) then
!write(std_out,'(a,i2,a,i2,a,es16.8,a,es16.8)')' JWZ Debug: <',iband,'|',jband,'> = ',&
!&                     smat_k(1,iband,jband),' + i ',smat_k(2,iband,jband)
!end if
!!
!end do   ! jband
!!
!end do   ! iband
!!
!ABI_DEALLOCATE(my_pwind_k)
!!
 do iband = 1, nband_occ

   if (sflag_k(iband) == 0) then

     pwmin = (iband-1)*npw_k1*nspinor*shiftbd
     pwmax = pwmin + npw_k1*nspinor
!    
!    old version  (*** multiply by nspinor missing??? ***)
!    vect1(:,1:npw_k1) = cg(:,icg + 1 + pwmin:icg + pwmax)
!    

!    Multiply the bra wave function by the phase factor
     if (itrs==1.or.itrs==11) then  ! take complex conjugate of bra
       do ispinor = 1, nspinor
         spnshft_k1=(ispinor-1)*npw_k1
         do ipw = 1,npw_k1
           vect1(1,spnshft_k1+ipw) = cg(1,icg+spnshft_k1+ipw+pwmin)*pwnsfac_k(1,ipw) &
&           +cg(2,icg+spnshft_k1+ipw+pwmin)*pwnsfac_k(2,ipw)
           vect1(2,spnshft_k1+ipw) = cg(1,icg+spnshft_k1+ipw+pwmin)*pwnsfac_k(2,ipw) &
&           -cg(2,icg+spnshft_k1+ipw+pwmin)*pwnsfac_k(1,ipw)
         end do
       end do
     else
       do ispinor=1, nspinor
         spnshft_k1=(ispinor-1)*npw_k1
         do ipw = 1,npw_k1
           vect1(1,spnshft_k1+ipw) = cg(1,icg+spnshft_k1+ipw+pwmin)*pwnsfac_k(1,ipw) &
&           -cg(2,icg+spnshft_k1+ipw+pwmin)*pwnsfac_k(2,ipw)
           vect1(2,spnshft_k1+ipw) = cg(1,icg+spnshft_k1+ipw+pwmin)*pwnsfac_k(2,ipw) &
&           +cg(2,icg+spnshft_k1+ipw+pwmin)*pwnsfac_k(1,ipw)
         end do
       end do
     end if

!    
     if (npw_k1*nspinor < mpw*nspinor) vect1(:,npw_k1*nspinor+1:mpw*nspinor) = zero

     do jband = 1, nband_occ

       pwmin = (jband-1)*npw_k2*nspinor
       pwmax = pwmin + npw_k2*nspinor

       if (itrs==10.or.itrs==11) then ! take complex conjugate of ket
         do ispinor=1, nspinor
           spnshft_k2=(ispinor-1)*npw_k2
           do ipw = 1, npw_k2
             vect2(1,spnshft_k2+ipw) = cgq(1,icg1+spnshft_k2+ipw+pwmin)*pwnsfac_k(3,ipw) &
&             +cgq(2,icg1+spnshft_k2+ipw+pwmin)*pwnsfac_k(4,ipw)
             vect2(2,spnshft_k2+ipw) = cgq(1,icg1+spnshft_k2+ipw+pwmin)*pwnsfac_k(4,ipw) &
&             -cgq(2,icg1+spnshft_k2+ipw+pwmin)*pwnsfac_k(3,ipw)
           end do
         end do
       else
         do ispinor=1, nspinor
           spnshft_k2=(ispinor-1)*npw_k2
           do ipw = 1, npw_k2
             vect2(1,spnshft_k2+ipw) = cgq(1,icg1+spnshft_k2+ipw+pwmin)*pwnsfac_k(3,ipw) &
&             -cgq(2,icg1+spnshft_k2+ipw+pwmin)*pwnsfac_k(4,ipw)
             vect2(2,spnshft_k2+ipw) = cgq(1,icg1+spnshft_k2+ipw+pwmin)*pwnsfac_k(4,ipw) &
&             +cgq(2,icg1+spnshft_k2+ipw+pwmin)*pwnsfac_k(3,ipw)
           end do
         end do
       end if

       if (npw_k2*nspinor < mpw*nspinor) vect2(:,npw_k2*nspinor+1:mpw*nspinor) = zero

!      DEBUG
!      if(iband==1 .and. jband==1 .and. itrs==0 .and. npw_k1==68 .and. npw_k2==74)then
!      write(std_out,'(a)' )' smatrix : ii,vect1,cg,pwnsfac='
!      do ii=1,npw_k1
!      write(std_out,'(i4,6es16.6)' )ii,vect1(:,ii),cg(:,icg+ii+pwmin),pwnsfac_k(1:2,ii)
!      end do
!      do ii=1,npw_k2
!      write(std_out,'(i4,6es16.6)' )ii,vect2(:,ii),cg(:,icg1+ii+pwmin),pwnsfac_k(3:4,ii)
!      end do
!      end if
!      ENDDEBUG

       call overlap_g(doti,dotr,mpw,npw_k1,npw_k2,nspinor,pwind_k,vect1,vect2)

       smat_k(1,iband,jband) = dotr
       smat_k(2,iband,jband) = doti

       if (usepaw == 1) then
         smat_k(1,iband,jband) = smat_k(1,iband,jband)+smat_k_paw(1,iband,jband)
         smat_k(2,iband,jband) = smat_k(2,iband,jband)+smat_k_paw(2,iband,jband)
       end if

!      DEBUG
!      if(iband==1 .and. jband==1)then
!      write(std_out,'(a,2es16.6,3i4)' )' smatrix : dotr,smat_k(1,iband,jband),mpw,npw_k1,npw_k2',dotr,smat_k(1,iband,jband),mpw,npw_k1,npw_k2
!      end if
!      ENDDEBUG

     end do   ! jband

   end if    ! sflag_k(iband) == 0

 end do   ! iband

!DEBUG
!do iband=1,nband_occ
!do jband=1,nband_occ
!write(std_out,'(a,2i4,2e20.10)') 'smat',iband,jband,smat_k(1,iband,jband),smat_k(2,iband,jband)
!end do
!end do
!write(std_out,*)' smatrix : smat_k(1,1,1)=',smat_k(1,1,1)
!ENDDEBUG

!Update sflag_k
 sflag_k(:) = 1

!Depending on the value of job, compute the determinant of the
!overlap matrix, its inverse or the product of the inverse
!overlap matrix with the WF at k.

 if ((job==1).or.(job==10).or.(job==11).or.(job==21)) then

   smat_inv(:,:,:) = smat_k(:,:,:)

!  DEBUG
!  write(std_out,*)' smatrix : smat_inv=',smat_inv
!  ENDDEBUG

   dzgedi_job=job; if(job==21) dzgedi_job=1
! TODO: should this be over nband_occ(isppol)?
   call dzgefa(smat_inv,mband_occ,nband_occ,ipvt,info)
   call dzgedi(smat_inv,mband_occ,nband_occ,ipvt,det,zgwork,dzgedi_job)

!  DEBUG
!  write(std_out,*)' smatrix : det=',det
!  ENDDEBUG

!  Compute the determinant of the overlap matrix
   dtm_k(:) = zero
   if (job==10 .or. job==11) then
     fac = exp(log(10._dp)*det(1,2))
     dtm_k(1) = fac*(det(1,1)*cos(log(10._dp)*det(2,2)) - &
&     det(2,1)*sin(log(10._dp)*det(2,2)))
     dtm_k(2) = fac*(det(1,1)*sin(log(10._dp)*det(2,2)) + &
&     det(2,1)*cos(log(10._dp)*det(2,2)))
   end if

!  Compute the product of the inverse overlap matrix with the WF

   if (ddkflag == 1) then

     cg1_k(:,:) = zero
     jband1 = 0

     if (itrs == 10 .or. itrs == 11) then

       do jband = minbd, maxbd
         jband1 = jband1 + 1
         do iband = 1, nband_occ

           do ispinor = 1, nspinor
             spnshft_k1 = (ispinor-1)*npw_k1
             spnshft_k2 = (ispinor-1)*npw_k2
             do ipw = 1, npw_k1

               jpw = pwind_k(ipw)
               
               if (jpw > 0) then
                 
                 wfr = cgq(1,icg1+(iband-1)*npw_k2*nspinor+spnshft_k2+jpw)*pwnsfac_k(3,jpw)&
&                 -cgq(2,icg1+(iband-1)*npw_k2*nspinor+spnshft_k2+jpw)*pwnsfac_k(4,jpw)
                 wfi = cgq(1,icg1+(iband-1)*npw_k2*nspinor+spnshft_k2+jpw)*pwnsfac_k(4,jpw)&
&                 +cgq(2,icg1+(iband-1)*npw_k2*nspinor+spnshft_k2+jpw)*pwnsfac_k(3,jpw)
                 
                 cg1_k(1,(jband1-1)*npw_k1*nspinor + spnshft_k1 + ipw) = &
&                 cg1_k(1,(jband1-1)*npw_k1*nspinor + spnshft_k1 + ipw) + &
&                 smat_inv(1,iband,jband)*wfr + smat_inv(2,iband,jband)*wfi

                 cg1_k(2,(jband1-1)*npw_k1*nspinor + spnshft_k1 + ipw) = &
&                 cg1_k(2,(jband1-1)*npw_k1*nspinor + spnshft_k1 + ipw) - &
&                 smat_inv(1,iband,jband)*wfi + smat_inv(2,iband,jband)*wfr

               end if

             end do ! end loop over npw_k1
           end do ! end loop over nspinor

         end do
       end do

     else

       do jband = minbd, maxbd
         jband1 = jband1 + 1
         do iband = 1, nband_occ

           do ispinor = 1, nspinor
             spnshft_k1 = (ispinor-1)*npw_k1
             spnshft_k2 = (ispinor-1)*npw_k2
             do ipw = 1, npw_k1

               jpw = pwind_k(ipw)

               if (jpw > 0) then

                 wfr = cgq(1,icg1+(iband-1)*npw_k2*nspinor+spnshft_k2+jpw)*pwnsfac_k(3,jpw)&
&                 -cgq(2,icg1+(iband-1)*npw_k2*nspinor+spnshft_k2+jpw)*pwnsfac_k(4,jpw)
                 wfi = cgq(1,icg1+(iband-1)*npw_k2*nspinor+spnshft_k2+jpw)*pwnsfac_k(4,jpw)&
&                 +cgq(2,icg1+(iband-1)*npw_k2*nspinor+spnshft_k2+jpw)*pwnsfac_k(3,jpw)

                 cg1_k(1,(jband1-1)*npw_k1*nspinor + spnshft_k1 + ipw) = &
&                 cg1_k(1,(jband1-1)*npw_k1*nspinor + spnshft_k1 + ipw) + &
&                 smat_inv(1,iband,jband)*wfr - smat_inv(2,iband,jband)*wfi

                 cg1_k(2,(jband1-1)*npw_k1*nspinor + spnshft_k1 + ipw) = &
&                 cg1_k(2,(jband1-1)*npw_k1*nspinor + spnshft_k1 + ipw) + &
&                 smat_inv(1,iband,jband)*wfi + smat_inv(2,iband,jband)*wfr

               end if

             end do ! end loop over npw_k1
           end do ! end loop over nspinor

         end do
       end do

     end if     ! itrs

   end if

 end if         ! 

 if(job == 20 .or. job == 21) then ! special case transfering cgq to cg1_k without use of S^{-1}, used in 
!  magnetic field case

   cg1_k(:,:) = zero
   jband1 = 0

   if (itrs == 10 .or. itrs == 11) then

     do jband = minbd, maxbd
       jband1 = jband1 + 1
       do ispinor = 1, nspinor
         spnshft_k1 = (ispinor-1)*npw_k1
         spnshft_k2 = (ispinor-1)*npw_k2
         do ipw = 1, npw_k1
           jpw = pwind_k(ipw)
           
           if (jpw > 0) then
             wfr = cgq(1,icg1+(jband1-1)*npw_k2*nspinor+spnshft_k2+jpw)*pwnsfac_k(3,jpw)&
&             -cgq(2,icg1+(jband1-1)*npw_k2*nspinor+spnshft_k2+jpw)*pwnsfac_k(4,jpw)
             wfi = cgq(1,icg1+(jband1-1)*npw_k2*nspinor+spnshft_k2+jpw)*pwnsfac_k(4,jpw)&
&             +cgq(2,icg1+(jband1-1)*npw_k2*nspinor+spnshft_k2+jpw)*pwnsfac_k(3,jpw)

             cg1_k(1,(jband1-1)*npw_k1*nspinor + spnshft_k1 + ipw) = wfr
             cg1_k(2,(jband1-1)*npw_k1*nspinor + spnshft_k1 + ipw) = wfi

           end if

         end do ! end loop over npw_k1
       end do ! end loop over nspinor

     end do

   else

     do jband = minbd, maxbd
       jband1 = jband1 + 1
       do ispinor = 1, nspinor
         spnshft_k1 = (ispinor-1)*npw_k1
         spnshft_k2 = (ispinor-1)*npw_k2
         do ipw = 1, npw_k1
           jpw = pwind_k(ipw)

           if (jpw > 0) then

             wfr = cgq(1,icg1+(jband1-1)*npw_k2*nspinor+spnshft_k2+jpw)*pwnsfac_k(3,jpw)&
&             -cgq(2,icg1+(jband1-1)*npw_k2*nspinor+spnshft_k2+jpw)*pwnsfac_k(4,jpw)
             wfi = cgq(1,icg1+(jband1-1)*npw_k2*nspinor+spnshft_k2+jpw)*pwnsfac_k(4,jpw)&
             +cgq(2,icg1+(jband1-1)*npw_k2*nspinor+spnshft_k2+jpw)*pwnsfac_k(3,jpw)

             cg1_k(1,(jband1-1)*npw_k1*nspinor + spnshft_k1 + ipw) = wfr
             cg1_k(2,(jband1-1)*npw_k1*nspinor + spnshft_k1 + ipw) = wfi

           end if

         end do ! end loop over npw_k1
       end do ! end loop over nspinor

     end do

   end if     ! itrs

 end if ! end job == 20 .or. job == 21 case

 ABI_DEALLOCATE(ipvt)
 ABI_DEALLOCATE(zgwork)
 ABI_DEALLOCATE(vect1)
 ABI_DEALLOCATE(vect2)

!DEBUG
!write(std_out,*)' dtm_k=',dtm_k(:)
!write(std_out,*)' smatrix : exit '
!ENDDEBUG

end subroutine smatrix
!!***
