!!****m* ABINIT/m_cgtk
!! NAME
!!  m_cgtk
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_cgtk

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_time

 use defs_abitypes, only : MPI_type
 use m_symtk,     only : mati3inv
 use m_geometry,  only : getspinrot
 use m_crystal,   only : crystal_t
 use m_fftcore,   only : sphere
 use m_kg,        only : ph1d3d, getph

 implicit none

 private
!!***

 public :: cgtk_rotate
 public :: cgtk_change_gsphere
 public :: cgtk_fixphase
!!***

contains
!!***

!!****f* ABINIT/cgtk_rotate
!! NAME
!!  cgtk_rotate
!!
!! FUNCTION
!!
!! INPUTS
!!  isym
!!  itimrev
!!  nspinor
!!  ndat
!!  npw1,npw2
!!  istwf1,istwf2
!!  cryst<crystal_t>
!!  shiftg(3)
!!  kg1(3,npw1), kg2(3,npw2)
!!  work_ngfft(18)
!!  kpt1(3)
!!  cg1(2,npw1,nspinor,ndat)
!!
!! OUTPUT
!!  cg2(2,npw2,nspinor,ndat)
!!  work(2,work_ngfft(4),work_ngfft(5),work_ngfft(6)) !*ndat),
!!
!! PARENTS
!!      m_sigmaph,m_wfd,m_wfk
!!
!! NOTES
!!  Inspired to wfconv.
!!
!! CHILDREN
!!      sphere
!!
!! SOURCE

subroutine cgtk_rotate(cryst, kpt1, isym, itimrev, shiftg, nspinor, ndat, &
  npw1, kg1, npw2, kg2, istwf1, istwf2, cg1, cg2, work_ngfft, work)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: isym,itimrev,nspinor,ndat,npw1,npw2,istwf1,istwf2
 type(crystal_t),intent(in) :: cryst
!arrays
 integer,intent(in) :: shiftg(3),kg1(3,npw1),kg2(3,npw2)
 integer,intent(in) :: work_ngfft(18)
 real(dp),intent(in) :: kpt1(3)
 real(dp),intent(in) :: cg1(2,npw1,nspinor,ndat)
 real(dp),intent(out) :: cg2(2,npw2,nspinor,ndat)
 real(dp),intent(out) :: work(2,work_ngfft(4),work_ngfft(5),work_ngfft(6)) !*ndat) for threads?

!Local variables ------------------------------
!scalars
 integer,parameter :: tobox=1,tosph=-1,me_g0=1
 integer :: n1,n2,n3,n4,n5,n6,ipw,idat,isp
 real(dp) :: arg,ar,ai,bi,br,spinrots,spinrotx,spinroty,spinrotz
 logical :: have_phase
!arrays
 integer,parameter :: no_shift(3)=0,atindx(1)=1
 integer :: symm(3,3),symrel_conv(3,3)
 real(dp) :: phktnons(2,1),tnons_conv(3),spinrot(4),tsec(2)
 real(dp),allocatable :: phase1d(:,:),phase3d(:,:),wavef1(:,:)

!************************************************************************

 ! Keep track of total time spent.
 call timab(1780, 1, tsec)

 n1=work_ngfft(1); n2=work_ngfft(2); n3=work_ngfft(3)
 n4=work_ngfft(4); n5=work_ngfft(5); n6=work_ngfft(6)

 symrel_conv = cryst%symrel(:,:,isym)
 call mati3inv(symrel_conv, symm)
 tnons_conv = cryst%tnons(:,isym) !
 have_phase = sum(tnons_conv**2) > tol8

 ! Compute rotation in spinor space
 if (nspinor == 2) call getspinrot(cryst%rprimd,spinrot,symrel_conv)
 if (itimrev == 1) symm = -symm

 ! Need to compute phase factors associated with nonsymmorphic translations?
 if (have_phase) then
   ABI_MALLOC(phase3d,(2,npw1))
   ABI_MALLOC(phase1d,(2,(2*n1+1)+(2*n2+1)+(2*n3+1)))
   ! Although the routine getph is originally written for atomic phase factors, it does precisely what we want
   call getph(atindx,1,n1,n2,n3,phase1d,tnons_conv)

   arg = two_pi*(kpt1(1)*tnons_conv(1)+ kpt1(2)*tnons_conv(2)+ kpt1(3)*tnons_conv(3))
   phktnons(1,1)=cos(arg)
   phktnons(2,1)=sin(arg)
   ! Convert 1D phase factors to 3D phase factors exp(i 2 pi (k+G).tnons )
   call ph1d3d(1,1,kg1,1,1,npw1,n1,n2,n3,phktnons,phase1d,phase3d)
 end if

 ABI_MALLOC(wavef1, (2,npw1))

 do idat=1,ndat
   do isp=1,nspinor
     wavef1 = cg1(:,:,isp,idat)

     if (have_phase) then
       ! Multiply by phase factors due to nonsymmorphic translations.
       do ipw=1,npw1
         ar = phase3d(1,ipw)*wavef1(1,ipw)-phase3d(2,ipw)*wavef1(2,ipw)
         ai = phase3d(2,ipw)*wavef1(1,ipw)+phase3d(1,ipw)*wavef1(2,ipw)
         wavef1(1,ipw)=ar
         wavef1(2,ipw)=ai
       end do
     end if

     ! Take into account time-reversal symmetry, if needed, in the scalar case
     if (itimrev == 1 .and. nspinor == 1) wavef1(2,:npw1) = -wavef1(2,:npw1)

     ! Insert wavef1 in work array.
     call sphere(wavef1,1,npw1,work,n1,n2,n3,n4,n5,n6,kg1,istwf1,tobox,me_g0,no_shift,identity_3d,one)

     ! Apply rotation and extract data.
     call sphere(cg2(:,:,isp,idat),1,npw2,work,n1,n2,n3,n4,n5,n6,kg2,istwf2,tosph,me_g0,shiftg,symm,one)
   end do ! isp

   if (nspinor == 2) then
     if (itimrev == 1) then
       ! Take care of time-reversal symmetry, if needed
       ! Exchange spin-up and spin-down. Make complex conjugate of one component,
       ! and change sign of other component
       do ipw=1,npw2
         ! Here, change sign of real part
         ar = -cg2(1,ipw,1,idat)
         ai =  cg2(2,ipw,1,idat)
         ! Here, change sign of imaginary part
         cg2(1,ipw,1,idat) =  cg2(1,ipw,2,idat)
         cg2(2,ipw,1,idat) = -cg2(2,ipw,2,idat)
         cg2(1,ipw,2,idat) = ar
         cg2(2,ipw,2,idat) = ai
       end do
     end if ! itimrev==1

     ! Rotation in spinor space (see also wfconv)
     spinrots = spinrot(1)
     spinrotx = spinrot(2)
     spinroty = spinrot(3)
     spinrotz = spinrot(4)
     do ipw=1,npw2
       ar = cg2(1,ipw,1,idat)
       ai = cg2(2,ipw,1,idat)
       br = cg2(1,ipw,2,idat)
       bi = cg2(2,ipw,2,idat)
       cg2(1,ipw,1,idat) =  spinrots*ar-spinrotz*ai +spinroty*br-spinrotx*bi
       cg2(2,ipw,1,idat) =  spinrots*ai+spinrotz*ar +spinroty*bi+spinrotx*br
       cg2(1,ipw,2,idat) = -spinroty*ar-spinrotx*ai +spinrots*br+spinrotz*bi
       cg2(2,ipw,2,idat) = -spinroty*ai+spinrotx*ar +spinrots*bi-spinrotz*br
     end do
   end if

 end do ! idat

 ABI_FREE(wavef1)

 if (have_phase) then
   ABI_FREE(phase3d)
   ABI_FREE(phase1d)
 end if

 call timab(1780, 2, tsec)

end subroutine cgtk_rotate
!!***

!!****f* ABINIT/cgtk_change_gsphere
!! NAME
!!  cgtk_change_gsphere
!!
!! FUNCTION
!!  Transfer the G components of the wavefunctions from one sphere to another one.
!!  Can be used to change the value of istwfk e.g. 2 --> 1
!!
!! INPUTS
!!  ndat = Number of wavefunctions to transform.
!!  npw1, npw2 = Number of plane-waves in (input, output) G-sphere
!!  istwf1, istwf2 = Storage mode of (input, output) wavefunctions.
!!  kg1(3,npw1), kg2(3,npw2) = Input/Output G-sphere
!!  cg1(2,npw1,ndat) = Input wavefunctions on kg1 sphere with istwf1 mode.
!!  work_ngfft(18)=Specify work dimensions. Must be large enough to accomodate kg1 and kg2
!!
!! OUTPUT
!!  cg2(2,npw2,ndat) = Output wavefunctions on kg2 sphere with istwf2 mode.
!!  work(2,work_ngfft(4),work_ngfft(5),work_ngfft(6)) = Workspace array
!!
!! PARENTS
!!      m_wfd
!!
!! CHILDREN
!!      sphere
!!
!! SOURCE

subroutine cgtk_change_gsphere(ndat, npw1, istwf1, kg1, cg1, npw2, istwf2, kg2, cg2, work_ngfft, work)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndat,npw1,npw2,istwf1,istwf2
!arrays
 integer,intent(in) :: kg1(3,npw1),kg2(3,npw2)
 integer,intent(in) :: work_ngfft(18)
 real(dp),intent(inout) :: cg1(2,npw1,ndat)  ! TODO: Should be intent(in) but need to change sphere
 real(dp),intent(out) :: cg2(2,npw2,ndat)
 real(dp),intent(out) :: work(2,work_ngfft(4),work_ngfft(5),work_ngfft(6))

!Local variables ------------------------------
!scalars
 integer,parameter :: tobox=1,tosph=-1,me_g0=1
 integer :: n1,n2,n3,n4,n5,n6,idat
!arrays
 integer,parameter :: no_shift(3)=0

!************************************************************************

 n1 = work_ngfft(1); n2 = work_ngfft(2); n3 = work_ngfft(3)
 n4 = work_ngfft(4); n5 = work_ngfft(5); n6 = work_ngfft(6)

 do idat=1,ndat
   ! Insert cg1 in work array taking into account istwf1 (intent in)
   call sphere(cg1(:,:,idat),1,npw1,work,n1,n2,n3,n4,n5,n6,kg1,istwf1,tobox,me_g0,no_shift,identity_3d,one)

   ! Extract cg2 from work array taking into account istwf2
   call sphere(cg2(:,:,idat),1,npw2,work,n1,n2,n3,n4,n5,n6,kg2,istwf2,tosph,me_g0,no_shift,identity_3d,one)
 end do

end subroutine cgtk_change_gsphere
!!***

!!****f* ABINIT/cgtk_fixphase
!! NAME
!! cgtk_fixphase
!!
!! FUNCTION
!! Fix phase of all bands. Keep normalization but maximize real part
!! (minimize imag part). Also fix the sign of real part
!! by setting the first non-zero element to be positive.
!! See also fxphas_seq in m_cgtools
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
!!  mpi_enreg=information about MPI parallelization
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
!! PARENTS
!!      m_vtowfk
!!
!! CHILDREN
!!      timab,xmpi_sum
!!
!! SOURCE

subroutine cgtk_fixphase(cg, gsc, icg, igsc, istwfk, mcg, mgsc, mpi_enreg, nband_k, npw_k, useoverlap)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icg,igsc,istwfk,mcg,mgsc,nband_k,npw_k,useoverlap
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 real(dp),intent(inout) :: cg(2,mcg),gsc(2,mgsc*useoverlap)

!Local variables-------------------------------
!scalars
 integer :: iband,ierr,ii,indx
 real(dp) :: cim,cre,gscim,gscre,quotient,root1,root2,saa,sab,sbb,theta
 real(dp) :: thppi,xx,yy
 character(len=500) :: msg
!arrays
 real(dp) :: buffer2(nband_k,2),buffer3(nband_k,3),tsec(2)
 real(dp),allocatable :: cimb(:),creb(:),saab(:),sabb(:),sbbb(:) !,sarr(:,:)

! *************************************************************************

!The general case, where a complex phase indeterminacy is present
 if(istwfk==1)then

   ABI_MALLOC(cimb,(nband_k))
   ABI_MALLOC(creb,(nband_k))
   ABI_MALLOC(saab,(nband_k))
   ABI_MALLOC(sabb,(nband_k))
   ABI_MALLOC(sbbb,(nband_k))
   cimb(:)=zero ; creb(:)=zero

!  Loop over bands
!  TODO: MG store saa arrays in sarr(3,nband_k) to reduce false sharing.
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nband_k,icg,npw_k,cg,saab,sbbb,sabb)
   do iband=1,nband_k
     indx=icg+(iband-1)*npw_k

!    Compute several sums over Re, Im parts of c
     saa=zero; sbb=zero; sab=zero
     do ii=1+indx,npw_k+indx
       saa=saa+cg(1,ii)*cg(1,ii)
       sbb=sbb+cg(2,ii)*cg(2,ii)
       sab=sab+cg(1,ii)*cg(2,ii)
     end do
     saab(iband)=saa
     sbbb(iband)=sbb
     sabb(iband)=sab
   end do

!  XG030513 : MPIWF : should transmit saab,sbbb,sabb from the procs
!  of the WF group to the master processor of the WF group
   if (mpi_enreg%paral_kgb == 1) then
     buffer3(:,1)=saab(:)
     buffer3(:,2)=sbbb(:)
     buffer3(:,3)=sabb(:)
     call timab(48,1,tsec)
     call xmpi_sum(buffer3,mpi_enreg%comm_fft,ierr)
     if (mpi_enreg%paral_spinor==1) then
       call xmpi_sum(buffer3,mpi_enreg%comm_spinor,ierr)
     end if
     call timab(48,2,tsec)
     saab(:)=buffer3(:,1)
     sbbb(:)=buffer3(:,2)
     sabb(:)=buffer3(:,3)
   end if

!  XG030513 : MPIWF this loop should only be executed by the master of the WF group

   if (mpi_enreg%paral_kgb==0.or.mpi_enreg%me_fft==0) then
     do iband=1,nband_k
       indx=icg+(iband-1)*npw_k

       saa=saab(iband)
       sbb=sbbb(iband)
       sab=sabb(iband)

!      Get phase angle theta
       if (sbb+saa>tol8)then
         if(abs(sbb-saa)>tol8*(sbb+saa) .or. 2*abs(sab)>tol8*(sbb+saa))then
           if (abs(sbb-saa)>tol8*abs(sab)) then
             quotient=sab/(sbb-saa)
             theta=0.5_dp*atan(2.0_dp*quotient)
           else
!            Taylor expansion of the atan in terms of inverse of its argument. Correct up to 1/x2, included.
             theta=0.25_dp*(pi-(sbb-saa)/sab)
           end if
!          Check roots to get theta for max Re part
           root1=cos(theta)**2*saa+sin(theta)**2*sbb-2.0_dp*cos(theta)*sin(theta)*sab
           thppi=theta+0.5_dp*pi
           root2=cos(thppi)**2*saa+sin(thppi)**2*sbb-2.0_dp*cos(thppi)*sin(thppi)*sab
           if (root2>root1) theta=thppi
         else
!          The real part vector and the imaginary part vector are orthogonal, and of same norm. Strong indeterminacy.
!          Will determine the first non-zero coefficient, and fix its phase
!          Hypothesis : there is at least one non-zero element on the master node ...
           do ii=1+indx,npw_k+indx
             cre=cg(1,ii)
             cim=cg(2,ii)
             if(cre**2+cim**2>tol8**2*(saa+sbb))then
               if(cre**2>tol8**2**cim**2)then
                 theta=atan(cim/cre)
               else
!                Taylor expansion of the atan in terms of inverse of its argument. Correct up to 1/x2, included.
                 theta=pi/2-cre/cim
               end if
               exit
             end if
           end do
         end if
       else
         write(msg,'(a,i0,5a)')&
&         'The eigenvector with band ',iband,' has zero norm.',ch10,&
&         'This usually happens when the number of bands (nband) is comparable to the number of planewaves (mpw)',ch10,&
&         'Action: Check the parameters of the calculation. If nband ~ mpw, then decrease nband or, alternatively, increase ecut'
         ABI_ERROR(msg)
       end if

       xx=cos(theta)
       yy=sin(theta)

!      Here, set the first non-zero element to be positive
!      Comment the next nine lines to recover the behaviour of pre v3.1.3g
!      Hypothesis : there is at least one non-zero element on the master node ...
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
   end if

!  XG030513 : MPIWF : should transmit creb(:),cimb(:) of the master
!  processor of the WF group to the others procs of the WF group
   if (mpi_enreg%paral_kgb == 1) then
     call timab(48,1,tsec)
     buffer2(:,1)=creb(:)
     buffer2(:,2)=cimb(:)
     call xmpi_sum(buffer2,mpi_enreg%comm_fft,ierr)
     if (mpi_enreg%paral_spinor==1) then
       call xmpi_sum(buffer2,mpi_enreg%comm_spinor,ierr)
     end if
     call timab(48,2,tsec)
     creb(:)=buffer2(:,1)
     cimb(:)=buffer2(:,2)
   end if

!  MG TODO: Scaling can be done with zscal
!$OMP PARALLEL DO PRIVATE(indx,xx,yy,cre,cim,gscre,gscim)
   do iband=1,nband_k
     indx=icg+(iband-1)*npw_k

     xx=creb(iband)
     yy=cimb(iband)
!    Alter phase of array |cg>
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

   ABI_FREE(cimb)
   ABI_FREE(creb)
   ABI_FREE(saab)
   ABI_FREE(sabb)
   ABI_FREE(sbbb)

 else  ! if istwfk/=1.  Storages that take into account the time-reversal symmetry : the freedom is only a sign freedom

   ABI_MALLOC(creb,(nband_k))
   creb(:)=zero
!  XG030513 : MPIWF : this loop should be done only by the master processor of the WF group

   if (mpi_enreg%paral_kgb==0.or.mpi_enreg%me_fft==0) then

!    Loop over bands
     do iband=1,nband_k

       indx=icg+(iband-1)*npw_k

!      Here, set the first non-zero real element to be positive
       do ii=1+indx,npw_k+indx
         cre=cg(1,ii)
         if(abs(cre)>tol8)exit
       end do
       creb(iband)=cre

     end do ! iband

   end if
!  XG030513 : MPIWF : should transmit cre(:) of the master processor of the WF group to the others
   if (mpi_enreg%paral_kgb == 1) then
     call timab(48,1,tsec)
     call xmpi_sum(creb,mpi_enreg%comm_fft,ierr)
     if (mpi_enreg%paral_spinor==1) then
       call xmpi_sum(creb,mpi_enreg%comm_spinor,ierr)
     end if
     call timab(48,2,tsec)
   end if

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

   ABI_FREE(creb)
 end if ! istwfk

end subroutine cgtk_fixphase
!!***

end module m_cgtk
!!***
