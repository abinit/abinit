!{\src2tex{textfont=tt}}
!!****f* ABINIT/testcode_ctqmc
!! NAME
!! testcode_ctqmc
!!
!! FUNCTION
!! Setup ultra simple hybridization to test CTQMC in simple situations.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! gtmp_nd
!! gw_tmp_nd
!! temp = temperature
!! dmftqmc_l = number of times slices
!! nflavor = number of flavor
!! testrot = 0/1 if rotation of hybridization is tested or not
!! testcode = 1 if tests are activated.
!! opt = 1/2 if pre or postprocessing of CTQMC data.
!!
!! OUTPUT
!! fw1_nd = non diagonal hybridization
!! fw1 = hybridization
!! umod = value of U 
!!  
!!
!! SIDE EFFECTS
!!  gtmp_nd  
!!  gw_tmp_nd
!!
!! NOTES
!!
!! PARENTS
!!      qmc_prep_ctqmc
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine testcode_ctqmc(dmftqmc_l,fw1_nd,fw1,gtmp_nd,gw_tmp_nd,levels_ctqmc,hybri_limit,&
&   nflavor,opt,temp,testrot,testcode,umod)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_ctqmc
 use m_CtqmcInterface
 use m_greenhyb
 !use m_self, only : self_type,initialize_self,destroy_self,print_self,rw_self
 use m_io_tools, only : flush_unit

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'testcode_ctqmc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: dmftqmc_l,nflavor,testrot,testcode,opt
 real(dp), intent(in) :: temp
 real(dp), intent(out) :: umod(2,2)
 complex(dpc), intent(inout) :: gw_tmp_nd(:,:,:)
 real(dp),  intent(inout) :: gtmp_nd(:,:,:)
 complex(dpc), intent(out) :: fw1(:,:)
 complex(dpc), intent(out) :: fw1_nd(:,:,:)
 real(dp),  intent(inout) :: levels_ctqmc(:)
 complex(dpc),  intent(inout) :: hybri_limit(:,:)

!Local variables ------------------------------
 character(len=500) :: message
 integer :: ifreq, itau,realrot,simplehyb
 real(dp) :: omega
 real(dp) :: tbi1,tbi2,e2,tbi3,tbi4,e3,e4,tbi21,tbi12,e3b,e4b,tbi21b,tbi12b
 complex(dpc) :: e1
! arrays
 complex(dpc) :: RR(2,2)
 complex(dpc) :: RR1(2,2)
 complex(dpc) :: RRi(2,2)
 complex(dpc) :: RRt(2,2)
! ************************************************************************
 if (testcode==0) return
 if (nflavor/=2) then
   write(message,'(2a)') ch10,' testcode nflavor.ne.2'
   MSG_ERROR(message)
 end if

 simplehyb=2
 simplehyb=1
 simplehyb=3
 !=========================
 ! Built rotation matrix
 !=========================
 realrot=0
 realrot=2
 if (realrot==1) then
   ! Real rotation
   !=========================
   RR(1,1)  =  SQRT(3.d0)/2.d0
   RR(1,2)  = -1.d0/2.d0
   RR(2,1)  =  1.d0/2.d0
   RR(2,2)  =  SQRT(3.d0)/2.d0
 else if (realrot==2) then
   ! Real rotation
   !=========================
   RR(1,1)  =  SQRT(1.d0/2.d0)
   RR(1,2)  = -SQRT(1.d0/2.d0)
   RR(2,1)  =  SQRT(1.d0/2.d0)
   RR(2,2)  =  SQRT(1.d0/2.d0)
 else
   ! Complex rotation
   !=========================
   RR(1,1)  =  CMPLX(one,two)
   RR(1,2)  =  CMPLX(one,one)
   RR(2,1)  =  CMPLX(one,-one)
   RR(2,2)  =  CMPLX(-one,two)
   RR=RR/sqrt(seven)
 end if
 ! Check rotation is unitary
 !==========================
 RRi(1,1) =  conjg(RR(1,1))
 RRi(1,2) =  conjg(RR(2,1))
 RRi(2,1) =  conjg(RR(1,2))
 RRi(2,2) =  conjg(RR(2,2))
 RR1(:,:)  = MATMUL ( RR(:,:) , RRi(:,:)          )
 !write(6,*) "RR1",RR1
 if(abs(RR1(1,1)-one).gt.tol7.or.abs(RR1(1,2)).gt.tol7.or.abs(RR1(2,2)-one).gt.tol7.or.abs(RR1(2,1)).gt.tol7) then
   write(message,'(2a)') ch10,' testcode error in rotation matrix'
   MSG_ERROR(message)
 end if


 !=================================
 ! Built hybridization  for CTQMC
 !=================================
 if (opt==1) then

 !  Parameters: tight-binding + U
 !  firt test of the code try umod=0, and (tbi1,tbi2,e1,e2)=(2,1,0.5,0.0) testrot=1
 !  second test of the code try umod=four, and (tbi1,tbi2,e1,e2)=(2,1,0.0,0.0) testrot=1
 !=======================================================================================
   fw1_nd(:,:,:)= czero
   tbi1=2.0_dp
   tbi2=1.0_dp
   tbi3=1.0_dp
   tbi4=1.0_dp
   tbi12=2.5_dp
   tbi12b=2.5_dp
   tbi21=2.5_dp
   tbi21b=2.5_dp
   e1=cmplx(0.0,0.0,8)
   e2=zero
   e3=0.2
   e4=0.3
   e3b=0.3
   e4b=-0.2
   umod(:,:)=0.d0

   if(testrot==1.and.(abs(tbi1-tbi2)<tol6)) then
     write(message,'(3a)') ch10,' testrot=1 with tbi1=tbi2 is equivalent' &
     ,'to testrot=0: change testrot'
     MSG_WARNING(message)
   end if
   ! Built fw1_nd
   !==============
   do ifreq=1,dmftqmc_l

     omega=pi*temp*(two*float(ifreq)-1)

     if(simplehyb==1) then
       fw1_nd(ifreq,1,1) =  -umod(1,1)/two+tbi1**2/(dcmplx(0.d0,omega)-e1)
       fw1_nd(ifreq,2,2) =  -umod(1,1)/two+tbi2**2/(dcmplx(0.d0,omega)-e2)
       fw1(ifreq,1)      =  -umod(1,1)/two+tbi1**2/(dcmplx(0.d0,omega)-e1)
       fw1(ifreq,2)      =  -umod(1,1)/two+tbi2**2/(dcmplx(0.d0,omega)-e2)
       hybri_limit(1,1)=tbi1**2
       hybri_limit(2,2)=tbi2**2
       hybri_limit(1,2)=0.d0
       hybri_limit(2,1)=0.d0
     else if(simplehyb==2) then
       fw1_nd(ifreq,1,1) =  -umod(1,1)/two+tbi1**2/(dcmplx(0.d0,omega)-e1)+tbi3**2/(dcmplx(0.d0,omega)-e3)
       fw1_nd(ifreq,2,2) =  -umod(1,1)/two+tbi2**2/(dcmplx(0.d0,omega)-e2)+tbi4**2/(dcmplx(0.d0,omega)-e4)
       fw1(ifreq,1)      =  -umod(1,1)/two+tbi1**2/(dcmplx(0.d0,omega)-e1)
       fw1(ifreq,2)      =  -umod(1,1)/two+tbi2**2/(dcmplx(0.d0,omega)-e2)
     else if(simplehyb==3) then
       fw1_nd(ifreq,1,1) =  -umod(1,1)/two+tbi1**2/(dcmplx(0.d0,omega)-e1)
       fw1_nd(ifreq,2,2) =  -umod(1,1)/two+tbi2**2/(dcmplx(0.d0,omega)-e2)
       fw1_nd(ifreq,1,2) =  tbi12**2/(dcmplx(0.d0,omega)-e3)+tbi12b**2/(dcmplx(0.d0,omega)-e3b)
       fw1_nd(ifreq,2,1) =  tbi21**2/(dcmplx(0.d0,omega)-e4)+tbi21b**2/(dcmplx(0.d0,omega)-e4b)
       fw1(ifreq,1)      =  -umod(1,1)/two+tbi1**2/(dcmplx(0.d0,omega)-e1)
       fw1(ifreq,2)      =  -umod(1,1)/two+tbi2**2/(dcmplx(0.d0,omega)-e2)
       hybri_limit(1,1)=tbi1**2
       hybri_limit(2,2)=tbi2**2
       hybri_limit(1,2)=tbi12**2+tbi12b**2
       hybri_limit(2,1)=tbi21**2+tbi21b**2
     end if
     write(132,*) omega,real(fw1_nd(ifreq,1,1)),aimag(fw1_nd(ifreq,1,1))
     write(133,*) omega,real(fw1_nd(ifreq,1,2)),aimag(fw1_nd(ifreq,1,2))
     write(134,*) omega,real(fw1_nd(ifreq,2,1)),aimag(fw1_nd(ifreq,2,1))
     write(135,*) omega,real(fw1_nd(ifreq,2,2)),aimag(fw1_nd(ifreq,2,2))
     write(1234,*) omega, real(fw1(ifreq,1)),aimag(fw1(ifreq,1))
   end do
   ! Built level and limit of hybridization
   !=======================================
   levels_ctqmc(1:nflavor)=-umod(1,1)/two

   write(std_out,*) "fw1_nd"
   write(std_out,*) fw1_nd(1,1,1), fw1_nd(1,1,2)
   write(std_out,*) fw1_nd(1,2,1), fw1_nd(1,2,2)
   write(std_out,*) "fw1"
   write(std_out,*) fw1(1,1), fw1(1,2)
   write(std_out,*) fw1(2,1), fw1(2,2)

 ! Rotate hybridization if testrot=1
 !==================================
   if(testrot==1) then

     do ifreq=1,dmftqmc_l
       RRt(:,:)  = MATMUL ( RR(:,:)  , fw1_nd(ifreq,:,:) )
   !write(6,*) "RRt"
   !write(6,*) RRt(1,1), RRt(1,2)
   !write(6,*) RRt(2,1), RRt(2,2)
       RR1(:,:)  = MATMUL ( RRt(:,:) , RRi(:,:)          )
   !write(6,*) "RR1"
   !write(6,*) RR1(1,1), RR1(1,2)
   !write(6,*) RR1(2,1), RR1(2,2)
       fw1_nd(ifreq,:,:)=RR1(:,:)
       omega=pi*temp*(two*float(ifreq)+1)
       write(3322,*) omega,real(fw1_nd(ifreq,1,1)),aimag(fw1_nd(ifreq,1,1))
       write(232,*) omega,real(fw1_nd(ifreq,1,1)),aimag(fw1_nd(ifreq,1,1))
       write(233,*) omega,real(fw1_nd(ifreq,1,2)),aimag(fw1_nd(ifreq,1,2))
       write(234,*) omega,real(fw1_nd(ifreq,2,1)),aimag(fw1_nd(ifreq,2,1))
       write(235,*) omega,real(fw1_nd(ifreq,2,2)),aimag(fw1_nd(ifreq,2,2))
     end do

     ! Rotate limit of hybridization
     !=======================================
     RRt(:,:)  = MATMUL ( RR(:,:)  , hybri_limit(:,:)  )
     RR1(:,:)  = MATMUL ( RRt(:,:) , RRi(:,:)          )
     hybri_limit(:,:)=RR1(:,:)

   end if
   ! rajouter test real(fw1_nd(1,:,:)) doit etre diagonale

 !======================================
 ! Rotate Green's function from CTQMC
 !======================================
 else if(opt==2) then

   write(std_out,*) "gw_tmp_nd"
   write(std_out,*) gw_tmp_nd(1,1,1), gw_tmp_nd(1,1,2)
   write(std_out,*) gw_tmp_nd(1,2,1), gw_tmp_nd(1,2,2)
   ! Rotate Green's function back
   !==============================
   if(testrot==1) then
     do ifreq=1,dmftqmc_l
       RRt(1:nflavor,1:nflavor) = MATMUL ( RRi(1:nflavor,1:nflavor),gw_tmp_nd(ifreq,1:nflavor,1:nflavor) )
       RR1(1:nflavor,1:nflavor) = MATMUL ( RRt(1:nflavor,1:nflavor),RR(1:nflavor,1:nflavor) )
       gw_tmp_nd(ifreq,1:nflavor,1:nflavor)=RR1(1:nflavor,1:nflavor)
     end do

     write(std_out,*) "gw_tmp_nd after rotation"
     write(std_out,*) gw_tmp_nd(1,1,1), gw_tmp_nd(1,1,2)
     write(std_out,*) gw_tmp_nd(1,2,1), gw_tmp_nd(1,2,2)

     do itau=1,dmftqmc_l
       RRt(1:nflavor,1:nflavor) = MATMUL ( RRi(1:nflavor,1:nflavor),gtmp_nd(itau,1:nflavor,1:nflavor) )
       RR1(1:nflavor,1:nflavor)  = MATMUL ( RRt(1:nflavor,1:nflavor),RR(1:nflavor,1:nflavor) )
       gtmp_nd(itau,1:nflavor,1:nflavor)=real(RR1(1:nflavor,1:nflavor))
     end do

   ! Rotate Green's function for comparison with testrot=1
   !======================================================
   else if (testrot==0) then ! produce rotated green's function to compare to testrot=1 case

     do itau=1,dmftqmc_l
       RRt(1:nflavor,1:nflavor) = MATMUL ( RR(1:nflavor,1:nflavor),gtmp_nd(itau,1:nflavor,1:nflavor) )
       RR1(1:nflavor,1:nflavor)  = MATMUL ( RRt(1:nflavor,1:nflavor),RRi(1:nflavor,1:nflavor) )
       write(444,*) real(itau-1)/(temp*real(dmftqmc_l)),real(RR1(1,1)),real(RR1(2,2)),real(RR1(1,2)),real(RR1(2,1))
     end do

   end if 

   ! Print out rotated Green's function
   !=====================================
   do itau=1,dmftqmc_l
     write(555,'(e14.5,4(2e14.5,3x))') real(itau-1)/(temp*real(dmftqmc_l)),gtmp_nd(itau,1,1),&
&     gtmp_nd(itau,2,2),gtmp_nd(itau,1,2),gtmp_nd(itau,2,1)
   end do
   
   write(message,'(2a)') ch10,' testcode end of test calculation'
   MSG_ERROR(message)

 end if
 close(444)
 close(555)

end subroutine testcode_ctqmc
!!***
