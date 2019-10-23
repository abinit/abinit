!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ddb_flexo
!! NAME
!!  m_ddb_flexo
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2019 ABINIT group (MR,MS)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
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

module m_ddb_flexo
    
 use defs_basis
 use m_abicore
 use m_profiling_abi
 use m_errors

 use m_fstrings,       only : sjoin
 use m_ddb
 use m_crystal,        only : crystal_t

 implicit none

 private

 public :: ddb_flexo
! *************************************************************************

contains 
!!***

!!****f* ABINIT/m_ddb_flexo/ddb_flexo
!! NAME
!!  ddb_flexo  
!!
!! FUNCTION
!! Get all the contributions to the flexoelectric tensor
!!
!! INPUTS
!!  ddb<type(ddb_type)>=2nd order derivative database.
!!  ddb<type(ddb_type)>=Long wave 3rd order derivative database.
!!  Crystal<type(crystal_t)>=Crystal structure parameters
!!  filnamddb = name of the ddb file
!!  flexoflg=  1 -> Computes all contributions to FxE
!!             2 -> Computes electronic (clamped ion) contribution to FxE
!!             3 -> Computes mixed contribution to FxE
!!             4 -> Computes lattice contribution to FxE
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  anaddb
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddb_flexo(ddb,ddb_lw,crystal,filnamddb,flexoflg)
    
 implicit none

!Arguments ------------------------------------
 integer , intent(in)  :: flexoflg
 class(ddb_type),intent(in) :: ddb,ddb_lw
 type(crystal_t),intent(in) :: crystal
 character(len=fnlen) :: filnamddb

!Local variables-------------------------------
 integer :: iblok,lwsym
! real(dp) :: 
 character(len=500) :: msg

!arrays
 integer :: rfelfd(4),rfphon(4),rfstrs(4)
 integer :: rfqvec(4)
 real(dp) :: qphnrm(3),qphon(3,3)
 real(dp) :: ciflexo(3,3,3,3)
 real(dp) :: mixflexo(3,3,3,3)
 real(dp) :: pol1(3,3,3,ddb%natom)
 
! *************************************************************************

 DBG_ENTER("COLL")
 
! First get the clamped-ion flexoelectric tensor
 if (flexoflg==1.or.flexoflg==2) then
  
   ! Look for the Gamma Block in the DDB
   qphon(:,:)=zero
   qphnrm(:)=one
   rfphon(3)=0
   rfelfd(3)=1
   rfstrs(3)=3
   rfqvec(3)=1
  
   call ddb_lw%get_block(iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,33,rfqvec=rfqvec)
  
   if (iblok == 0) then
     call wrtout(std_out, "--- !WARNING")
     call wrtout(std_out, sjoin("- Cannot find clamped ion FxE tensor in DDB file:", filnamddb))
     call wrtout(std_out, "flexoflag=1 or 2 requires the DDB file to include the corresponding long wave 3rd derivatives")
   end if

   call dtciflexo(ddb_lw%val(:,:,iblok),ddb%mpert,ddb%natom,ciflexo,crystal%ucvol)

 end if

! Then get the mixed contribution to the flexoelectric tensor
 if (flexoflg==1.or.flexoflg==3) then

   ! Extract the P^(1) tensor from the DDB
   lwsym=0
   iblok = ddb_lw%get_quadrupoles(crystal,lwsym,33,pol1)

   ! Look for the Gamma Block of the Phi^(1) tensor in the DDB
   qphon(:,:)=zero
   qphnrm(:)=one
   rfphon(3)=1
   rfelfd(3)=0
   rfstrs(3)=0
   rfqvec(3)=1
  
   call ddb_lw%get_block(iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,33,rfqvec=rfqvec)

   if (iblok == 0) then
     call wrtout(std_out, "--- !WARNING")
     call wrtout(std_out, sjoin("- Cannot find P^(1) tensor in DDB file:", filnamddb))
     call wrtout(std_out, "flexoflag=1 or 3 requires the DDB file to include the corresponding long wave 3rd derivatives")
   end if

   call dtmixflexo(ddb_lw%val(:,:,iblok),mixflexo,ddb%mpert,ddb%natom,pol1,crystal%ucvol)

 end if

 DBG_EXIT("COLL")

end subroutine ddb_flexo
!!***

!!****f* m_ddb/dtciflexo
!! NAME
!! dtciflexo
!!
!! FUNCTION
!! Reads the Clamped Ion Flexoelectric Tensor
!! in the Gamma Block coming from the Derivative Data Base
!! (long wave third-order derivatives).
!!
!! INPUTS
!! blkval(2,3*mpert*3*mpert*3*mpert)= matrix of third-order energies
!! mpert =maximum number of ipert
!! natom= number of atoms in unit cell
!! ucvol= unit cell volume
!!
!! OUTPUT
!! ciflexo(3,3,3,3) = type-II Clamped Ion Flexoelectric Tensor
!!
!! PARENTS
!!      m_ddb
!!
!! CHILDREN
!!
!! SOURCE

subroutine dtciflexo(blkval,mpert,natom,ciflexo,ucvol)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,natom
 real(dp),intent(in) :: ucvol
!arrays
 real(dp),intent(in) :: blkval(2,3*mpert*3*mpert*3*mpert)
 real(dp),intent(out) :: ciflexo(3,3,3,3)

!Local variables -------------------------
!scalars
 integer :: elfd,istrs,ivarA,strsd,strsd1,strsd2,strst,qvecd
 logical :: iwrite
 real(dp) :: fac
 character(len=500) :: msg
!arrays
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 real(dp) :: d3cart(2,3,mpert,3,mpert,3,mpert)

! *********************************************************************

 d3cart(1,:,:,:,:,:,:) = reshape(blkval(1,:),shape = (/3,mpert,3,mpert,3,mpert/))
 d3cart(2,:,:,:,:,:,:) = reshape(blkval(2,:),shape = (/3,mpert,3,mpert,3,mpert/))

!Extraction of the macroscopic polarization induced by an atomic displacement (P^(1) tensor) 
 fac=two/ucvol
 do qvecd=1,3
   do istrs=1,6
     strsd1=alpha(istrs)
     strsd2=beta(istrs)
     strst=natom+3; if (istrs>3) strst=natom+4
     strsd=istrs; if (istrs>3) strsd=istrs-3
     do elfd=1,3
       ciflexo(elfd,qvecd,strsd1,strsd2)=-fac*d3cart(2,elfd,natom+2,strsd,strst,qvecd,natom+8)
       if (istrs>3) ciflexo(elfd,qvecd,strsd2,strsd1)=ciflexo(elfd,qvecd,strsd1,strsd2)
     end do
   end do
 end do

!Print results
 iwrite = ab_out > 0
 if (iwrite) then
   write(msg,'(3a)')ch10,' Type-II electronic flexoelectric tensor (clamped ion) ',ch10
   call wrtout(ab_out,msg,'COLL')
   write(ab_out,*)'      xx          yy          zz          yz          xz          xy          zy          zx          yx'
   do ivarA=1,6
     elfd=alpha(ivarA)
     qvecd=beta(ivarA)
     write(ab_out,'(9f12.6)') ciflexo(elfd,qvecd,1,1),ciflexo(elfd,qvecd,2,2),ciflexo(elfd,qvecd,3,3),&
                              ciflexo(elfd,qvecd,2,3),ciflexo(elfd,qvecd,1,3),ciflexo(elfd,qvecd,1,2),&
                              ciflexo(elfd,qvecd,3,2),ciflexo(elfd,qvecd,3,1),ciflexo(elfd,qvecd,2,1)
   end do
   do ivarA=4,6
     elfd=beta(ivarA)
     qvecd=alpha(ivarA)
     write(ab_out,'(9f12.6)') ciflexo(elfd,qvecd,1,1),ciflexo(elfd,qvecd,2,2),ciflexo(elfd,qvecd,3,3),&
                              ciflexo(elfd,qvecd,2,3),ciflexo(elfd,qvecd,1,3),ciflexo(elfd,qvecd,1,2),&
                              ciflexo(elfd,qvecd,3,2),ciflexo(elfd,qvecd,3,1),ciflexo(elfd,qvecd,2,1)
   end do
 end if

 end subroutine dtciflexo
!!***

!!****f* m_ddb/dtmixflexo
!! NAME
!! dtmixflexo
!!
!! FUNCTION
!! Reads the P^(1) and Phi^(1) tensors 
!! in the Gamma Block coming from the Derivative Data Base
!! (long wave third-order derivatives). And computes the mixed
!! contribution to the flexoelectric tensor.
!!
!! INPUTS
!! blkval(2,3*mpert*3*mpert*3*mpert)= matrix of third-order energies for Phi^(1) tensor
!! mpert =maximum number of ipert
!! natom= number of atoms in unit cell
!! pol1(3,3,3,natom)= tensor with the polarization induced by an atomic displacement (P^(1))
!! ucvol= unit cell volume
!!
!! OUTPUT
!! mixflexo(3,3,3,3) = type-II mixed contribution to the Flexoelectric Tensor
!!
!! PARENTS
!!      m_ddb
!!
!! CHILDREN
!!
!! SOURCE

subroutine dtmixflexo(blkval,mixflexo,mpert,natom,pol1,ucvol)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,natom
 real(dp),intent(in) :: ucvol
!arrays
 real(dp),intent(in) :: blkval(2,3*mpert*3*mpert*3*mpert)
 real(dp),intent(inout) :: pol1(3,3,3,natom)
 real(dp),intent(out) :: mixflexo(3,3,3,3)

!Local variables -------------------------
!scalars
 integer :: elfd,iat,iatd,jat,jatd,ivarA,qvecd
 logical :: iwrite
 real(dp) :: fac
 character(len=500) :: msg
!arrays
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 real(dp) :: d3cart(2,3,mpert,3,mpert,3,mpert)
 real(dp) :: phi1(3,natom,3,natom,3)

! *********************************************************************

 d3cart(1,:,:,:,:,:,:) = reshape(blkval(1,:),shape = (/3,mpert,3,mpert,3,mpert/))
 d3cart(2,:,:,:,:,:,:) = reshape(blkval(2,:),shape = (/3,mpert,3,mpert,3,mpert/))

!P^(1) lacks the 1/ucvol factor
 pol1=pol1/ucvol

!Extraction of Phi^(1) tensor
 do qvecd=1,3
   do jat=1,natom
     do jatd=1,3
       do iat=1,natom
         do iatd=1,3
           phi1(iatd,iat,jatd,jat,qvecd)=-two*d3cart(2,iatd,iat,jatd,jat,qvecd,natom+8)
         end do
       end do
     end do
   end do
 end do

!Print results
! iwrite = ab_out > 0
! if (iwrite) then
!   write(msg,'(3a)')ch10,' Type-II electronic flexoelectric tensor (clamped ion) ',ch10
!   call wrtout(ab_out,msg,'COLL')
!   write(ab_out,*)'      xx          yy          zz          yz          xz          xy          zy          zx          yx'
!   do ivarA=1,6
!     elfd=alpha(ivarA)
!     qvecd=beta(ivarA)
!     write(ab_out,'(9f12.6)') ciflexo(elfd,qvecd,1,1),ciflexo(elfd,qvecd,2,2),ciflexo(elfd,qvecd,3,3),&
!                              ciflexo(elfd,qvecd,2,3),ciflexo(elfd,qvecd,1,3),ciflexo(elfd,qvecd,1,2),&
!                              ciflexo(elfd,qvecd,3,2),ciflexo(elfd,qvecd,3,1),ciflexo(elfd,qvecd,2,1)
!   end do
!   do ivarA=4,6
!     elfd=beta(ivarA)
!     qvecd=alpha(ivarA)
!     write(ab_out,'(9f12.6)') ciflexo(elfd,qvecd,1,1),ciflexo(elfd,qvecd,2,2),ciflexo(elfd,qvecd,3,3),&
!                              ciflexo(elfd,qvecd,2,3),ciflexo(elfd,qvecd,1,3),ciflexo(elfd,qvecd,1,2),&
!                              ciflexo(elfd,qvecd,3,2),ciflexo(elfd,qvecd,3,1),ciflexo(elfd,qvecd,2,1)
!   end do
! end if

 end subroutine dtmixflexo
!!***

end module m_ddb_flexo
!!***
