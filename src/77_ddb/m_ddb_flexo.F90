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

 use m_fstrings,       only : itoa,sjoin
 use m_ddb
 use m_crystal,        only : crystal_t
 use m_dynmat,       only : asria_corr

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
!!  asr= if /=0 acustic sume rule is imposed on the dynamical matrix
!!  d2asr(2,3,natom,3,natom)=ASR-correction
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

subroutine ddb_flexo(asr,d2asr,ddb,ddb_lw,crystal,filnamddb,flexoflg)
    
 implicit none

!Arguments ------------------------------------
!scalars
 integer , intent(in)  :: asr,flexoflg
 class(ddb_type),intent(in) :: ddb,ddb_lw
 type(crystal_t),intent(in) :: crystal
 character(len=fnlen) :: filnamddb
!arrays
 real(dp),intent(in) :: d2asr(2,3,ddb%natom,3,ddb%natom)

!Local variables-------------------------------
 integer :: fblok,iblok,jblok,lwsym
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
     call wrtout(std_out, "   flexoflag=1 or 3 requires the DDB file to include the corresponding long wave 3rd derivatives")
   end if

   ! Look for th block that contains the forces 
   qphon(:,:)=zero
   qphnrm(:)=one
   rfphon(4)=1
   rfelfd(:)=0
   rfstrs(:)=0
   call ddb%get_block(fblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,4)

   if (fblok == 0) then
     call wrtout(std_out, "--- !WARNING")
     call wrtout(std_out, sjoin("- Cannot find forces in DDB file:", filnamddb))
     call wrtout(std_out, "   If there are nonzero residual atomic forces on the structure") 
     call wrtout(std_out, "   flexoflag=1 or 3 will produce an improper piezoelectric force response tensor")
     call wrtout(std_out, "   and a wrong value for the  mixed contribution to the flexoelectric tensor")
   end if

   ! Look for the Gamma Block of the dynamical matrix in the DDB
   qphon(:,:)=zero
   qphnrm(:)=one
   rfphon(:)=0;rfphon(1:2)=1
   rfelfd(:)=0
   rfstrs(:)=0
  
   call ddb%get_block(jblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,1)

   if (jblok == 0) then
     call wrtout(std_out, "--- !WARNING")
     call wrtout(std_out, sjoin("- Cannot find Gamma point Dynamical Matrix in DDB file:", filnamddb))
     call wrtout(std_out, "   flexoflag=1 or 3 requires the DDB file to include the corresponding 2nd derivatives")
   end if

   call dtmixflexo(asr,d2asr,ddb%val(:,:,fblok),ddb%val(:,:,jblok),ddb_lw%val(:,:,iblok),mixflexo,ddb%mpert,ddb%natom,pol1,crystal%ucvol)

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
!! asr= if /=0 acustic sume rule is imposed on the dynamical matrix
!! d2asr(2,3,natom,3,natom)=ASR-correction
!! blkval1d(2,3,mpert,3,mpert)= 1st derivative wrt atom displacements (at least)
!! blkval2d(2,3,mpert,3,mpert)= 2nd derivatives wrt two atom displacements (at least)
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

subroutine dtmixflexo(asr,d2asr,blkval1d,blkval2d,blkval,mixflexo,mpert,natom,pol1,ucvol)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: asr,mpert,natom
 real(dp),intent(in) :: ucvol
!arrays
 real(dp),intent(in) :: d2asr(2,3,mpert,3,mpert)
 real(dp),intent(in) :: blkval1d(1,3,mpert,3,mpert)
 real(dp),intent(in) :: blkval2d(2,3,mpert,3,mpert)
 real(dp),intent(in) :: blkval(2,3*mpert*3*mpert*3*mpert)
 real(dp),intent(inout) :: pol1(3,3,3,natom)
 real(dp),intent(out) :: mixflexo(3,3,3,3)

!Local variables -------------------------
!scalars
 integer :: elfd,iat,iatd,ivar,jat,jatd,jvar,kat,katd,qvecd,qvecd2
 logical :: iwrite
 character(len=500) :: msg
!arrays
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 real(dp) :: d3cart(2,3,mpert,3,mpert,3,mpert)
 real(dp) :: forces(3,natom)
 real(dp) :: intstrn(3,3,3,natom)
 real(dp) :: phi1(3,natom,3,natom,3)
 real(dp) :: piezofr(3,natom,3,3)
 real(dp) :: psinvdm(3*natom,3*natom)

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

!Extraction of the forces to acount for the improper contribution
 do iatd=1,3
   do iat=1,natom
     forces(iatd,iat)=-blkval1d(1,iatd,iat,1,1)
   end do
 end do

!Calculate the piezoelectric force-response tensor including the improper contribution
 piezofr(:,:,:,:)=zero
 do qvecd=1,3
   do iat=1,natom
     do iatd=1,3
       do jatd=1,3
         do jat=1,natom
           piezofr(iatd,iat,jatd,qvecd)=piezofr(iatd,iat,jatd,qvecd) + phi1(iatd,iat,jatd,jat,qvecd)
         end do
         if (iatd==qvecd) piezofr(iatd,iat,jatd,qvecd)=piezofr(iatd,iat,jatd,qvecd) + forces(jatd,iat)
       end do
     end do
   end do
 end do

!Calculate the ion-relaxed internal strain tensor
 !First we need to obtain the pseudo-inverse of the dynamical matrix 
 call dm_psinv(asr,blkval2d,d2asr,ab_out,psinvdm,mpert,natom)

 !Perfom the product with the piezo force-response
 intstrn(:,:,:,:)=zero
 do qvecd=1,3
   do katd=1,3
     do iatd=1,3
       do iat=1,natom
         ivar=(iat-1)*3+iatd
         do jatd=1,3
           do jat=1,natom
             jvar=(jat-1)*3+jatd

             intstrn(qvecd,katd,iatd,iat)= intstrn(qvecd,katd,iatd,iat) + &
           psinvdm(ivar,jvar)*piezofr(jatd,jat,katd,qvecd)

           end do
         end do
       end do
     end do
   end do
 end do

 !Finally calculate the mixed contribution to the FxE tensor
  mixflexo(:,:,:,:)=zero
  do elfd=1,3
    do qvecd=1,3
      do katd=1,3
        do qvecd2=1,3
          do iatd=1,3
            do iat=1,natom

              mixflexo(elfd,qvecd,katd,qvecd2)=mixflexo(elfd,qvecd,katd,qvecd2) - &
            pol1(elfd,qvecd,iatd,iat)*intstrn(qvecd2,katd,iatd,iat)

           end do
         end do
       end do
     end do
   end do
 end do

!Print results
 iwrite = ab_out > 0
 if (iwrite) then
   write(msg,'(3a)')ch10,' Displacement-response internal strain tensor from long-wave magnitudes ',ch10
   call wrtout(ab_out,msg,'COLL')
   write(ab_out,*)' atom   dir        xx          yy          zz          yz          xz          xy'
   do iat=1,natom
     write(ab_out,'(2x,i3,3x,a3,2x,6f12.6)') iat, 'x', intstrn(1,1,1,iat),intstrn(2,2,1,iat),intstrn(3,3,1,iat),&
                                                         intstrn(2,3,1,iat),intstrn(1,3,1,iat),intstrn(1,2,1,iat)
     write(ab_out,'(2x,i3,3x,a3,2x,6f12.6)') iat, 'y', intstrn(1,1,2,iat),intstrn(2,2,2,iat),intstrn(3,3,2,iat),&
                                                         intstrn(2,3,2,iat),intstrn(1,3,2,iat),intstrn(1,2,2,iat)
     write(ab_out,'(2x,i3,3x,a3,2x,6f12.6)') iat, 'z', intstrn(1,1,3,iat),intstrn(2,2,3,iat),intstrn(3,3,3,iat),&
                                                         intstrn(2,3,3,iat),intstrn(1,3,3,iat),intstrn(1,2,3,iat)
   end do 

   write(msg,'(3a)')ch10,' Type-II mixed contribution to flexoelectric tensor ',ch10
   call wrtout(ab_out,msg,'COLL')
   write(ab_out,*)'      xx          yy          zz          yz          xz          xy          zy          zx          yx'
   do ivar=1,6
     elfd=alpha(ivar)
     qvecd=beta(ivar)
     write(ab_out,'(9f12.6)') mixflexo(elfd,qvecd,1,1),mixflexo(elfd,qvecd,2,2),mixflexo(elfd,qvecd,3,3),&
                              mixflexo(elfd,qvecd,2,3),mixflexo(elfd,qvecd,1,3),mixflexo(elfd,qvecd,1,2),&
                              mixflexo(elfd,qvecd,3,2),mixflexo(elfd,qvecd,3,1),mixflexo(elfd,qvecd,2,1)
   end do
   do ivar=4,6
     elfd=beta(ivar)
     qvecd=alpha(ivar)
     write(ab_out,'(9f12.6)') mixflexo(elfd,qvecd,1,1),mixflexo(elfd,qvecd,2,2),mixflexo(elfd,qvecd,3,3),&
                              mixflexo(elfd,qvecd,2,3),mixflexo(elfd,qvecd,1,3),mixflexo(elfd,qvecd,1,2),&
                              mixflexo(elfd,qvecd,3,2),mixflexo(elfd,qvecd,3,1),mixflexo(elfd,qvecd,2,1)
   end do
 end if

 end subroutine dtmixflexo
!!***

!!****f* m_ddb/dm_psinv
!! NAME
!! dm_psinv
!!
!! FUNCTION
!! Computes the pseudoinverse of the dynamical matrix (PRB 72, 035105 (2005)).
!! This piece of code is copied from m_ddb_internalstr.F90.
!!
!! INPUTS
!! asr= if /=0 acustic sume rule is imposed on the dynamical matrix
!! blkval(2,3,mpert,3,mpert)= 2nd derivatives wrt two atom displacements (at least)
!! d2asr(2,3,natom,3,natom)=ASR-correction
!! iout=out file number
!! mpert=maximum number of ipert
!! natom=number of atoms in unit cell
!! 
!! OUTPUT
!! kmatrix(3*natom,3*natom) = array with the pseudo-inverse of dynamical matrix
!!
!! PARENTS
!!      m_ddb_flexo
!!
!! CHILDREN
!!      asria_corr,wrtout,zhpev
!!
!! SOURCE

 subroutine dm_psinv(asr,blkval,d2asr,iout,kmatrix,mpert,natom)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: asr,iout,mpert,natom
!arrays
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert)
 real(dp),intent(in) :: d2asr(2,3,natom,3,natom)
 real(dp),intent(out) :: kmatrix(3*natom,3*natom)

!Local variables -------------------------
!scalars
 integer :: ii1,ipertA,ivarA
 integer :: ii2,ipertB,ivarB
 integer :: ier
 character(len=500) :: message
!arrays
 real(dp) :: Amatr(3*natom-3,3*natom-3),Apmatr(3*natom,3*natom)
 real(dp) :: Bmatr(2,((3*natom-3)*(3*natom-2))/2)
 real(dp) :: Bpmatr(2,(3*natom*(3*natom+1))/2),Cmatr(3*natom-3,3*natom-3)
 real(dp) :: Cpmatr(3*natom,3*natom),Nmatr(3*natom,3*natom)
 real(dp) :: eigval(3*natom-3),eigvalp(3*natom),eigvec(2,3*natom-3,3*natom-3)
 real(dp) :: eigvecp(2,3*natom,3*natom)
 real(dp) :: zhpev1(2,2*3*natom-4)
 real(dp) :: zhpev1p(2,2*3*natom-1),zhpev2(3*3*natom-5),zhpev2p(3*3*natom-2)
 real(dp) :: d2cart(2,3*natom,3*natom)

! *********************************************************************

 d2cart = zero
 do ipertA=1,natom
   do ii1=1,3
     ivarA=ii1+3*(ipertA-1)
     do ipertB=1,natom
       do ii2=1,3
         ivarB=ii2+3*(ipertB-1)
         d2cart(1,ivarA,ivarB)=blkval(1,ii1,ipertA,ii2,ipertB)
       end do
     end do
   end do
 end do

!Eventually impose the acoustic sum rule
!FIXME: this might depend on ifcflag: impose that it is 0 or generalize
 call asria_corr(asr,d2asr,d2cart,natom,natom)
 !call asrq0_apply(asrq0, natom, mpert, msize, crystal%xcart, d2cart)
 kmatrix = d2cart(1,:,:)
 Apmatr(:,:)=kmatrix(:,:)

!DEBUG
!write(std_out,'(/,a,/)')'the force constant matrix'
!do ivarA=1,3*natom
!write(std_out,'(/)')
!do ivarB=1,3*natom
!write(std_out,'(es16.6)')kmatrix(ivarB,ivarA)
!end do
!end do
!ENDDEBUG

 Nmatr(:,:)=0.0_dp
 do ivarA=1,3*natom
   do ivarB=1,3*natom
     if (mod(ivarA,3)==0 .and. mod(ivarB,3)==0)then
       Nmatr(ivarA,ivarB)=one
     end if
     if (mod(ivarA,3)==1 .and. mod(ivarB,3)==1)then
       Nmatr(ivarA,ivarB)=one
     end if
     if (mod(ivarA,3)==2 .and. mod(ivarB,3)==2)then
       Nmatr(ivarA,ivarB)=one
     end if
   end do
 end do

!DEBUG
!do ivarA=1,3*natom
!write(std_out,'(/)')
!do ivarB=1,3*natom
!write(std_out,'(es16.6)')Nmatr(ivarB,ivarA)
!end do
!end do
!ENDDEBUG

!starting the pseudoinverting processes
!then get the eigenvectors of the big matrix,give values to matrixBp
 Bpmatr=0.0_dp
 ii1=1
 do ivarA=1,3*natom
   do ivarB=1,ivarA
     Bpmatr(1,ii1)=Nmatr(ivarB,ivarA)
     ii1=ii1+1
   end do
 end do

!Bpmatr(2,:) is the imaginary part of the force matrix
!then call the subroutines CHPEV and ZHPEV to get the eigenvectors
 call ZHPEV ('V','U',3*natom,Bpmatr,eigvalp,eigvecp,3*natom,zhpev1p,zhpev2p,ier)
 ABI_CHECK(ier == 0, sjoin("ZHPEV returned:", itoa(ier)))

!DEBUG
!the eigenval and eigenvec
!write(std_out,'(/,a,/)')'the eigenvalues and eigenvectors'
!do ivarA=1,3*natom
!write(std_out,'(/)')
!write(std_out,'(es16.6)')eigvalp(ivarA)
!end do
!do ivarA=1,3*natom
!write(std_out,'(/)')
!do ivarB=1,3*natom
!write(std_out,'(es16.6)')eigvecp(1,ivarB,ivarA)
!end do
!end do
!ENDDEBUG

!Then do the multiplication to get the reduced matrix,in two steps
!After this the force constant matrix is decouple in two bloks,
!acoustic and optical ones
 Cpmatr(:,:)=0.0_dp
 do ivarA=1,3*natom
   do ivarB=1,3*natom
     do ii1=1,3*natom
       Cpmatr(ivarA,ivarB)=Cpmatr(ivarA,ivarB)+eigvecp(1,ii1,ivarA)*Apmatr(ii1,ivarB)
     end do
   end do
 end do

 Apmatr(:,:)=0.0_dp
 do ivarA=1,3*natom
   do ivarB=1,3*natom
     do ii1=1,3*natom
       Apmatr(ivarA,ivarB)=Apmatr(ivarA,ivarB)+Cpmatr(ivarA,ii1)*eigvecp(1,ii1,ivarB)
     end do
   end do
 end do

!DEBUG
!the blok diago
!write(std_out,'(/,a,/)')'matrixAp'
!do ivarA=1,3*natom
!write(std_out,'(/)')
!do ivarB=1,3*natom
!write(std_out,'(es16.6)')Apmatr(ivarA,ivarB)
!end do
!end do
!ENDDEBUG

!Check the last three eigenvalues whether too large or not
 ivarB=0
 do ivarA=3*natom-2,3*natom
   if (ABS(Apmatr(ivarA,ivarA))>tol6)then
     ivarB=1
   end if
 end do

 if(ivarB==1)then
   write(message,'(a,a,a,a,a,a,a,a,3es16.6)')ch10,&
&   '  Acoustic sum rule violation met : the eigenvalues of accoustic mode',ch10,&
&   '  are too large at Gamma point.',ch10,&
&   '  Increase cutoff energy or k-points sampling.',ch10,&
&   '  The three eigenvalues are:',Apmatr(3*natom-2,3*natom-2),Apmatr(3*natom-1,natom-1),Apmatr(3*natom,3*natom)
   MSG_WARNING(message)
   call wrtout(iout,message,'COLL')
 end if

!Give the value of reduced matrix form Apmatr to Amatr
 do ivarA=1,3*natom-3
   do ivarB=1,3*natom-3
     Amatr(ivarA,ivarB)=Apmatr(ivarA,ivarB)
   end do
 end do

!Now the reduced matrix is in the matrixA, the convert it
!first give the give the value of matixB from matrixA
 ii1=1
 do ivarA=1,3*natom-3
   do ivarB=1,ivarA
     Bmatr(1,ii1)=Amatr(ivarB,ivarA)
     ii1=ii1+1
   end do
 end do
 Bmatr(2,:)=0.0_dp

!Call the subroutines CHPEV and ZHPEV to get the eigenvectors and the eigenvalues
 call ZHPEV ('V','U',3*natom-3,Bmatr,eigval,eigvec,3*natom-3,zhpev1,zhpev2,ier)
 ABI_CHECK(ier == 0, sjoin("ZHPEV returned:", itoa(ier)))

!Check the unstable phonon modes, if the first is negative then print
!warning message
 if(eigval(1)<-1.0*tol8)then
   write(message,'(9a)') ch10,&
&   ' --- !WARNING',ch10,&
&   '     Unstable eigenvalue detected in force constant matrix at Gamma point',ch10,&
&   '     The system under calculation is physically unstable.',ch10,&
&   ' ---',ch10
   call wrtout(std_out,message,'COLL')
 end if

!Do the matrix mutiplication to get pseudoinverse inverse matrix
 Cmatr(:,:)=0.0_dp
 Amatr(:,:)=0.0_dp
 do ivarA=1,3*natom-3
   Cmatr(ivarA,ivarA)=1.0_dp/eigval(ivarA)
 end do

 do ivarA=1,3*natom-3
   do ivarB=1,3*natom-3
     do ii1=1,3*natom-3
       Amatr(ivarA,ivarB)=Amatr(ivarA,ivarB)+eigvec(1,ivarA,ii1)*Cmatr(ii1,ivarB)
     end do
   end do
 end do


!The second multiplication
 Cmatr(:,:)=0.0_dp
 do ivarA=1,3*natom-3
   do ivarB=1,3*natom-3
     do ii1=1,3*natom-3
       Cmatr(ivarA,ivarB)=Cmatr(ivarA,ivarB)+ Amatr(ivarA,ii1)*eigvec(1,ivarB,ii1)
     end do
   end do
 end do

!DEBUG
!write(std_out,'(/,a,/)')'the pseudo inverse of the force matrix'
!do ivarA=1,3*natom
!write(std_out,'(/)')
!do ivarB=1,3*natom
!write(std_out,'(es16.6)')Cmatr(ivarA,ivarB)
!end do
!end do
!ENDDEBUG

!So now the inverse of the reduced matrix is in the matrixC
!now do another mutilplication to get the pseudoinverse of the original
 Cpmatr(:,:)=0.0_dp
 Apmatr(:,:)=0.0_dp
 do ivarA=1,3*natom-3
   do ivarB=1,3*natom-3
     Cpmatr(ivarA,ivarB)=Cmatr(ivarA,ivarB)
   end do
 end do

!Now times the eigvecp
 do ivarA=1,3*natom
   do ivarB=1,3*natom
     do ii1=1,3*natom
       Apmatr(ivarA,ivarB)=Apmatr(ivarA,ivarB)+eigvecp(1,ivarA,ii1)*&
&       Cpmatr(ii1,ivarB)
     end do
   end do
 end do
 Cpmatr(:,:)=0.0_dp
 do ivarA=1,3*natom
   do ivarB=1,3*natom
     do ii1=1,3*natom
       Cpmatr(ivarA,ivarB)=Cpmatr(ivarA,ivarB)+ Apmatr(ivarA,ii1)*eigvecp(1,ivarB,ii1)
     end do
   end do
 end do

!Now the inverse is in Cpmatr
 kmatrix(:,:)=Cpmatr(:,:)
!transfer the inverse of k-matrix back to the k matrix
!so now the inverse of k matrix is in the kmatrix
!ending the part for pseudoinversing the K matrix

 end subroutine dm_psinv
!!***
end module m_ddb_flexo
!!***
