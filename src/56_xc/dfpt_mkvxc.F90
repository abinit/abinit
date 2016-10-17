!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_mkvxc
!! NAME
!! dfpt_mkvxc
!!
!! FUNCTION
!! Compute the first-order change of exchange-correlation potential
!! due to atomic displacement : assemble the first-order
!! density change with the frozen-core density change, then use
!! the exchange-correlation kernel.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2016 ABINIT group (XG, DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cplex= if 1, real space 1-order functions on FFT grid are REAL,
!!         if 2, COMPLEX
!!  ixc= choice of exchange-correlation scheme
!!  kxc(nfft,nkxc)=exchange and correlation kernel (see rhohxc.f)
!!  mpi_enreg=information about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!     see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nhat1(cplex*nfft,2nspden*nhat1dim)= -PAW only- 1st-order compensation density
!!  nhat1dim= -PAW only- 1 if nhat1 array is used ; 0 otherwise
!!  nhat1gr(cplex*nfft,nspden,3*nhat1grdim)= -PAW only- gradients of 1st-order compensation density
!!  nhat1grdim= -PAW only- 1 if nhat1gr array is used ; 0 otherwise
!!  nkxc=second dimension of the kxc array
!!  nspden=number of spin-density components
!!  n3xccc=dimension of xccc3d1 ; 0 if no XC core correction is used, otherwise, nfft
!!  option=if 0, work only with the XC core-correction,
!!         if 1, treat both density change and XC core correction
!!         if 2, treat only density change
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  rhor1(cplex*nfft,nspden)=array for electron density in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  usexcnhat= -PAW only- 1 if nhat density has to be taken into account in Vxc
!!  xccc3d1(cplex*n3xccc)=3D change in core charge density, see n3xccc
!!
!! OUTPUT
!!  vxc1(cplex*nfft,nspden)=change in exchange-correlation potential (including
!!   core-correction, if applicable)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      dfpt_dyxc1,dfpt_mkvxc_noncoll,dfpt_nstdy,dfpt_nstpaw,dfpt_rhotov
!!      dfptnl_loop,m_kxc,nres2vres
!!
!! CHILDREN
!!      dfpt_mkvxcgga,matr3inv,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine dfpt_mkvxc(cplex,ixc,kxc,mpi_enreg,nfft,ngfft,nhat1,nhat1dim,nhat1gr,nhat1grdim,&
&          nkxc,nspden,n3xccc,option,paral_kgb,qphon,rhor1,rprimd,usexcnhat,vxc1,xccc3d1)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_mkvxc'
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_56_xc, except_this_one => dfpt_mkvxc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,ixc,n3xccc,nfft,nhat1dim,nhat1grdim
 integer,intent(in) :: nkxc,nspden,option,paral_kgb,usexcnhat
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: nhat1(cplex*nfft,nspden*nhat1dim)
 real(dp),intent(in) :: nhat1gr(cplex*nfft,nspden,3*nhat1grdim)
 real(dp),intent(in) :: kxc(nfft,nkxc),qphon(3)
 real(dp),intent(in),target :: rhor1(cplex*nfft,nspden)
 real(dp),intent(in) :: rprimd(3,3),xccc3d1(cplex*n3xccc)
 real(dp),intent(out) :: vxc1(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: ii,ir,nhat1dim_tmp,nhat1grdim_tmp
 real(dp) :: rho1_dn,rho1_up,rho1im_dn,rho1im_up,rho1re_dn,rho1re_up
 !character(len=500) :: message
!arrays
 real(dp) :: gprimd(3,3),tsec(2)
 real(dp),allocatable :: nhat1tmp(:,:),nhat1grtmp(:,:,:),rhor1tmp(:,:)
 real(dp),pointer :: rhor1_(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(181,1,tsec)

 if(nspden/=1 .and. nspden/=2) then
   MSG_BUG('For nspden==4 please use dfpt_mkvxc_noncoll!')
 end if

!Special case: no XC applied
 if (ixc==0.or.nkxc==0) then
   MSG_WARNING('Note that no xc is applied (ixc=0)')
   vxc1=zero
   return
 end if

!Treat first LDA
 if(nkxc/=23)then

!  PAW: eventually substract compensation density
   if (option/=0) then
     if (usexcnhat==0.and.nhat1dim==1) then
       ABI_ALLOCATE(rhor1_,(cplex*nfft,nspden))
       rhor1_(:,:)=rhor1(:,:)-nhat1(:,:)
     else
       rhor1_ => rhor1
     end if
   end if

!  Case without non-linear core correction
   if(n3xccc==0 .or. option==2)then

     if(option==0)then  ! No straight XC to compute

       vxc1(:,:)=zero

     else               ! XC, without non-linear XC correction

!      Non-spin-polarized
       if(nspden==1)then
         if(cplex==1)then
           do ir=1,nfft
             vxc1(ir,1)=kxc(ir,1)*rhor1_(ir,1)
           end do
         else
           do ir=1,nfft
             vxc1(2*ir-1,1)=kxc(ir,1)*rhor1_(2*ir-1,1)
             vxc1(2*ir  ,1)=kxc(ir,1)*rhor1_(2*ir  ,1)
           end do
         end if ! cplex==1

!        Spin-polarized
       else
         if(cplex==1)then
           do ir=1,nfft
             rho1_dn=rhor1_(ir,1)-rhor1_(ir,2)
             vxc1(ir,1)=kxc(ir,1)*rhor1_(ir,2)+kxc(ir,2)*rho1_dn
             vxc1(ir,2)=kxc(ir,2)*rhor1_(ir,2)+kxc(ir,3)*rho1_dn
           end do
         else
           do ir=1,nfft
             rho1re_dn=rhor1_(2*ir-1,1)-rhor1_(2*ir-1,2)
             rho1im_dn=rhor1_(2*ir  ,1)-rhor1_(2*ir  ,2)
             vxc1(2*ir-1,1)=kxc(ir,1)*rhor1_(2*ir-1,2)+kxc(ir,2)*rho1re_dn
             vxc1(2*ir  ,1)=kxc(ir,1)*rhor1_(2*ir  ,2)+kxc(ir,2)*rho1im_dn
             vxc1(2*ir-1,2)=kxc(ir,2)*rhor1_(2*ir-1,2)+kxc(ir,3)*rho1re_dn
             vxc1(2*ir  ,2)=kxc(ir,2)*rhor1_(2*ir  ,2)+kxc(ir,3)*rho1im_dn
           end do
         end if ! cplex==1
       end if ! nspden==1

     end if ! option==0

!    Treat case with non-linear core correction
   else

     if(option==0)then

       if(nspden==1)then
         if(cplex==1)then
           do ir=1,nfft
             vxc1(ir,1)=kxc(ir,1)*xccc3d1(ir)
           end do
         else
           do ir=1,nfft
             vxc1(2*ir-1,1)=kxc(ir,1)*xccc3d1(2*ir-1)
             vxc1(2*ir  ,1)=kxc(ir,1)*xccc3d1(2*ir  )
           end do
         end if ! cplex==1
       else
         if(cplex==1)then
           do ir=1,nfft
             vxc1(ir,1)=(kxc(ir,1)+kxc(ir,2))*xccc3d1(ir)*half
             vxc1(ir,2)=(kxc(ir,2)+kxc(ir,3))*xccc3d1(ir)*half
           end do
         else
           do ir=1,nfft
             vxc1(2*ir-1,1)=(kxc(ir,1)+kxc(ir,2))*xccc3d1(2*ir-1)*half
             vxc1(2*ir  ,1)=(kxc(ir,1)+kxc(ir,2))*xccc3d1(2*ir  )*half
             vxc1(2*ir-1,2)=(kxc(ir,2)+kxc(ir,3))*xccc3d1(2*ir-1)*half
             vxc1(2*ir  ,2)=(kxc(ir,2)+kxc(ir,3))*xccc3d1(2*ir  )*half
           end do
         end if ! cplex==1
       end if ! nspden==1

     else ! option/=0

       if(nspden==1)then
         if(cplex==1)then
           do ir=1,nfft
             vxc1(ir,1)=kxc(ir,1)*(rhor1_(ir,1)+xccc3d1(ir))
           end do
         else
           do ir=1,nfft
             vxc1(2*ir-1,1)=kxc(ir,1)*(rhor1_(2*ir-1,1)+xccc3d1(2*ir-1))
             vxc1(2*ir  ,1)=kxc(ir,1)*(rhor1_(2*ir  ,1)+xccc3d1(2*ir  ))
           end do
         end if ! cplex==1
       else
         if(cplex==1)then
           do ir=1,nfft
             rho1_dn=rhor1_(ir,1)-rhor1_(ir,2) + xccc3d1(ir)*half
             rho1_up=rhor1_(ir,2)             + xccc3d1(ir)*half
             vxc1(ir,1)=kxc(ir,1)*rho1_up+kxc(ir,2)*rho1_dn
             vxc1(ir,2)=kxc(ir,2)*rho1_up+kxc(ir,3)*rho1_dn
           end do
         else
           do ir=1,nfft
             rho1re_dn=rhor1_(2*ir-1,1)-rhor1_(2*ir-1,2) + xccc3d1(2*ir-1)*half
             rho1im_dn=rhor1_(2*ir  ,1)-rhor1_(2*ir  ,2) + xccc3d1(2*ir  )*half
             rho1re_up=rhor1_(2*ir-1,2)                 + xccc3d1(2*ir-1)*half
             rho1im_up=rhor1_(2*ir  ,2)                 + xccc3d1(2*ir  )*half
             vxc1(2*ir-1,1)=kxc(ir,1)*rho1re_up+kxc(ir,2)*rho1re_dn
             vxc1(2*ir  ,1)=kxc(ir,1)*rho1im_up+kxc(ir,2)*rho1im_dn
             vxc1(2*ir-1,2)=kxc(ir,2)*rho1re_up+kxc(ir,3)*rho1re_dn
             vxc1(2*ir  ,2)=kxc(ir,2)*rho1im_up+kxc(ir,3)*rho1im_dn
           end do
         end if ! cplex==1
       end if ! nspden==1

     end if ! option==0

   end if ! n3xccc==0

   if (option/=0.and.usexcnhat==0.and.nhat1dim==1) then
     ABI_DEALLOCATE(rhor1_)
   end if

!  Treat GGA
 else

   ABI_ALLOCATE(rhor1tmp,(cplex*nfft,2))

!  Generates gprimd
   call matr3inv(rprimd,gprimd)

!  First transfer the data to spin-polarized storage
   if(option==1 .or. option==2)then   ! Treat the density change
     if(nspden==1)then
       do ir=1,cplex*nfft
         rhor1tmp(ir,1)=rhor1(ir,1)*half
         rhor1tmp(ir,2)=rhor1(ir,1)*half
       end do
     else
       do ir=1,cplex*nfft
         rho1_dn=rhor1(ir,1)-rhor1(ir,2)
         rhor1tmp(ir,1)=rhor1(ir,2)
         rhor1tmp(ir,2)=rho1_dn
       end do
     end if ! nspden==1
   else
     do ir=1,cplex*nfft
       rhor1tmp(ir,1)=zero
       rhor1tmp(ir,2)=zero
     end do
   end if
   if( (option==0 .or. option==1) .and. n3xccc/=0)then
     do ir=1,cplex*nfft
       rhor1tmp(ir,1)=rhor1tmp(ir,1)+xccc3d1(ir)*half
       rhor1tmp(ir,2)=rhor1tmp(ir,2)+xccc3d1(ir)*half
     end do
   end if

!  PAW: treat also compensation density (and gradients)
   if (option/=0.and.nhat1dim==1) then
     ABI_ALLOCATE(nhat1tmp,(cplex*nfft,2))
     nhat1dim_tmp=1
     if (nspden==1)then
       nhat1tmp(:,1)=nhat1(:,1)*half
       nhat1tmp(:,2)=nhat1(:,1)*half
     else
       do ir=1,cplex*nfft
         rho1_dn=nhat1(ir,1)-nhat1(ir,2)
         nhat1tmp(ir,1)=nhat1(ir,2)
         nhat1tmp(ir,2)=rho1_dn
       end do
     end if
   else
     ABI_ALLOCATE(nhat1tmp,(0,0))
     nhat1dim_tmp=0
   end if
   if (option/=0.and.nhat1grdim==1) then
     ABI_ALLOCATE(nhat1grtmp,(cplex*nfft,2,3))
     nhat1grdim_tmp=1
     if (nspden==1)then
       do ii=1,3
         nhat1grtmp(:,1,ii)=nhat1gr(:,1,ii)*half
         nhat1grtmp(:,2,ii)=nhat1gr(:,1,ii)*half
       end do
     else
       do ii=1,3
         do ir=1,cplex*nfft
           rho1_dn=nhat1gr(ir,1,ii)-nhat1gr(ir,2,ii)
           nhat1grtmp(ir,1,ii)=nhat1gr(ir,2,ii)
           nhat1grtmp(ir,2,ii)=rho1_dn
         end do
       end do
     end if
   else
     ABI_ALLOCATE(nhat1grtmp,(0,0,0))
     nhat1grdim_tmp=0
   end if

   call dfpt_mkvxcgga(cplex,gprimd,kxc,mpi_enreg,nfft,ngfft,nhat1tmp,nhat1dim_tmp,&
&   nhat1grtmp,nhat1grdim_tmp,nkxc,nspden,paral_kgb,qphon,rhor1tmp,usexcnhat,vxc1)

   ABI_DEALLOCATE(rhor1tmp)
   ABI_DEALLOCATE(nhat1tmp)
   ABI_DEALLOCATE(nhat1grtmp)

 end if

 call timab(181,2,tsec)

 DBG_EXIT("COLL")

end subroutine dfpt_mkvxc
!!***
