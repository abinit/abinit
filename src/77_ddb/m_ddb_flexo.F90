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
 use m_dynmat,       only : asria_corr,cart39

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
!!  zeff(3,3,natom)= Born Effective charges
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

subroutine ddb_flexo(asr,d2asr,ddb,ddb_lw,crystal,filnamddb,flexoflg,zeff)
    
 implicit none

!Arguments ------------------------------------
!scalars
 integer , intent(in)  :: asr,flexoflg
 class(ddb_type),intent(in) :: ddb,ddb_lw
 type(crystal_t),intent(in) :: crystal
 character(len=fnlen) :: filnamddb
!arrays
 real(dp),intent(in) :: d2asr(2,3,ddb%natom,3,ddb%natom)
 real(dp),intent(in) :: zeff(3,3,ddb%natom)

!Local variables-------------------------------
 integer :: elfd,iblok,ivar,jblok,kblok,lblok,lwsym,qvecd
 logical :: intstrn_only,iwrite
 character(len=500) :: msg

!arrays
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 integer :: rfelfd(4),rfphon(4),rfstrs(4)
 integer :: rfqvec(4)
 real(dp) :: qphnrm(3),qphon(3,3)
 real(dp) :: ciflexo(3,3,3,3)
 real(dp) :: intstrn(3,3,3,ddb%natom)
 real(dp) :: lattflexo(3,3,3,3)
 real(dp) :: mixflexo(3,3,3,3)
 real(dp) :: pol1(3,3,3,ddb%natom)
 real(dp) :: psinvdm(3*ddb%natom,3*ddb%natom)
 real(dp) :: totflexo(3,3,3,3)
 character(len=2) :: voigt(9)=(/'xx','yy','zz','yz','xz','xy','zy','zx','yx'/)
 
! *************************************************************************

 DBG_ENTER("COLL")
 
! First get the clamped-ion flexoelectric tensor
 ciflexo(:,:,:,:)=zero
 if (flexoflg==1.or.flexoflg==2) then
  
   rfphon(:)=0
   rfelfd(:)=0
   rfstrs(:)=0
   rfqvec(:)=0

   ! Look for the Gamma Block in the DDB
   qphon(:,:)=zero
   qphnrm(:)=one
   rfphon(3)=0
   rfelfd(3)=1
   rfstrs(3)=3
   rfqvec(3)=1
  
   write(msg, '(2a)' ) ch10," Extract the electronic flexoelectric coeficients from 3DTE"
   call wrtout(std_out,msg,'COLL') 
   call ddb_lw%get_block(iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,33,rfqvec=rfqvec)
  
   if (iblok == 0) then
     call wrtout(std_out, "  ")
     call wrtout(std_out, "--- !WARNING")
     call wrtout(std_out, sjoin("- Cannot find clamped ion FxE tensor in DDB file:", filnamddb))
     call wrtout(std_out, "  flexoflag=1 or 2 requires the DDB file to include the corresponding long wave 3rd derivatives")
   end if

   call dtciflexo(ddb_lw%val(:,:,iblok),ddb%mpert,ddb%natom,ciflexo,crystal%ucvol)

 end if

! Then get the mixed contribution to the flexoelectric tensor
 !Activate the calculation of internal strain necessary for lattice mediated contribution.
 intstrn_only=.false.;if (flexoflg==4) intstrn_only=.true.
 mixflexo(:,:,:,:)=zero
 if (flexoflg==1.or.flexoflg==3.or.intstrn_only) then

   ! Extract the P^(1) tensor from the DDB
   if (.not.intstrn_only) then
     lwsym=0
     iblok = ddb_lw%get_quadrupoles(lwsym,33,pol1)
   end if

   rfphon(:)=0
   rfelfd(:)=0
   rfstrs(:)=0
   rfqvec(:)=0

   ! Look for the Gamma Block of the Phi^(1) tensor in the DDB
   qphon(:,:)=zero
   qphnrm(:)=one
   rfphon(3)=1
   rfelfd(3)=0
   rfstrs(3)=0
   rfqvec(3)=1
  
   write(msg, '(2a)' ) ch10," Extract the Phi^(1) coeficients from 3DTE"
   call wrtout(std_out,msg,'COLL') 
   call ddb_lw%get_block(iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,33,rfqvec=rfqvec)

   if (iblok == 0) then
     call wrtout(std_out, "  ")
     call wrtout(std_out, "--- !WARNING")
     call wrtout(std_out, sjoin("- Cannot find Phi^(1) tensor in DDB file:", filnamddb))
     call wrtout(std_out, "  flexoflag=1 or 3 requires the DDB file to include the corresponding long wave 3rd derivatives")
   end if

   ! Look for th block that contains the forces 
   qphon(:,:)=zero
   qphnrm(:)=one
   rfphon(4)=1
   rfelfd(4)=0
   rfstrs(4)=0

   write(msg, '(2a)' ) ch10," Extract the forces from 1DTE"
   call wrtout(std_out,msg,'COLL') 
   call ddb%get_block(kblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,4)

   if (kblok == 0) then
     call wrtout(std_out, "  ")
     call wrtout(std_out, "--- !WARNING")
     call wrtout(std_out, sjoin("- Cannot find forces in DDB file:", filnamddb))
     call wrtout(std_out, "  If there are nonzero residual atomic forces on the structure") 
     call wrtout(std_out, "  flexoflag=1 or 3 will produce an improper piezoelectric force response tensor")
     call wrtout(std_out, "  and a wrong value for the mixed contribution to the flexoelectric tensor")
   end if

   ! Look for the Gamma Block of the dynamical matrix in the DDB
   qphon(:,:)=zero
   qphnrm(:)=one
   rfphon(1:2)=1
   rfelfd(1:2)=0
   rfstrs(1:2)=0
  
   write(msg, '(2a)' ) ch10," Extract the Dynamical Matrix from 2DTE"
   call wrtout(std_out,msg,'COLL') 
   call ddb%get_block(jblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,1)

   if (jblok == 0) then
     call wrtout(std_out, "  ")
     call wrtout(std_out, "--- !WARNING")
     call wrtout(std_out, sjoin("- Cannot find Gamma point Dynamical Matrix in DDB file:", filnamddb))
     call wrtout(std_out, "  flexoflag=1 or 3 requires the DDB file to include the corresponding 2nd derivatives")
   end if

   call dtmixflexo(asr,d2asr,ddb%val(:,:,kblok),ddb%val(:,:,jblok),ddb_lw%val(:,:,iblok),crystal%gprimd,&
  & intstrn,intstrn_only,mixflexo,ddb%mpert,ddb%natom,pol1,psinvdm,crystal%rprimd,crystal%ucvol)

 end if

! Finally get the lattice mediated contribution to the flexoelectric tensor
 lattflexo(:,:,:,:)=zero
 if (flexoflg==1.or.flexoflg==4) then

   rfphon(:)=0
   rfelfd(:)=0
   rfstrs(:)=0
   rfqvec(:)=0

   ! Look for the Gamma Block of the Phi^(1) tensor in the DDB
   qphon(:,:)=zero
   qphnrm(:)=one
   rfphon(3)=1
   rfelfd(3)=0
   rfstrs(3)=0
   rfqvec(3)=1

   write(msg, '(2a)' ) ch10," Extract the Phi^(1) coeficients from 3DTE"
   call wrtout(std_out,msg,'COLL')
   call ddb_lw%get_block(iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,33,rfqvec=rfqvec)

   if (iblok == 0) then
     call wrtout(std_out, "  ")
     call wrtout(std_out, "--- !WARNING")
     call wrtout(std_out, sjoin("- Cannot find Phi^(1) tensor in DDB file:", filnamddb))
     call wrtout(std_out, "  flexoflag=1 or 4 requires the DDB file to include the corresponding long wave 3rd derivatives")
   end if

   ! Look for the Gamma Block of the flexoelectric force response tensor in the DDB
   qphon(:,:)=zero
   qphnrm(:)=one
   rfphon(3)=1
   rfelfd(3)=0
   rfstrs(3)=3
   rfqvec(3)=1

   write(msg, '(2a)' ) ch10," Extract the FxE force response coeficients from 3DTE"
   call wrtout(std_out,msg,'COLL')
   call ddb_lw%get_block(jblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,33,rfqvec=rfqvec)

   if (jblok == 0) then
     call wrtout(std_out, "  ")
     call wrtout(std_out, "--- !WARNING")
     call wrtout(std_out, sjoin("- Cannot find FxE force response tensor in DDB file:", filnamddb))
     call wrtout(std_out, "  flexoflag=1 or 4 requires the DDB file to include the corresponding long wave 3rd derivatives")
   end if

   ! Look for the stress tensor in the DDB
   qphon(:,:)=zero
   qphnrm(:)=one
   rfphon(4)=0
   rfelfd(4)=0
   rfstrs(4)=3

   write(msg, '(2a)' ) ch10," Extract the stress tensor from 1DTE"
   call wrtout(std_out,msg,'COLL')
   call ddb%get_block(lblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,4)

   if (lblok == 0) then
     call wrtout(std_out, "  ")
     call wrtout(std_out, "--- !WARNING")
     call wrtout(std_out, sjoin("- Cannot find stress tensor in DDB file:", filnamddb))
     call wrtout(std_out, "  flexoflag=1 or 4 requires the DDB file to include the corresponding 2nd derivatives")
     call wrtout(std_out, "  to compute the Lagrange Elastic Tensor ")
   end if

   call dtlattflexo(ddb%amu,ddb%val(:,:,lblok),ddb_lw%val(:,:,jblok),ddb_lw%val(:,:,iblok),&
  & intstrn,lattflexo,ddb%mpert,ddb%natom,crystal%ntypat,psinvdm,crystal%typat,crystal%ucvol,zeff)

 end if

!Merge the three contributions and print the total FxE tensor
 totflexo(:,:,:,:)=ciflexo(:,:,:,:)+mixflexo(:,:,:,:)+lattflexo(:,:,:,:)

 iwrite = ab_out > 0
 if (iwrite) then
   write(msg,'(3a)')ch10,' TOTAL flexoelectric tensor (units= nC/m) ',ch10
   call wrtout([ab_out,std_out],msg,'COLL')
   write(msg,*)'           xx          yy          zz          yz          xz          xy'
   call wrtout([ab_out,std_out],msg,'COLL')
   do ivar=1,6
     elfd=alpha(ivar)
     qvecd=beta(ivar)
     write(msg,'(3x,a2,6f12.6)') voigt(ivar),totflexo(elfd,qvecd,1,1),totflexo(elfd,qvecd,2,2),totflexo(elfd,qvecd,3,3),&
                           totflexo(elfd,qvecd,2,3),totflexo(elfd,qvecd,1,3),totflexo(elfd,qvecd,1,2)
     call wrtout([ab_out,std_out],msg,'COLL')
   end do
   do ivar=4,6
     elfd=beta(ivar)
     qvecd=alpha(ivar)
     write(msg,'(3x,a2,6f12.6)') voigt(ivar+3),totflexo(elfd,qvecd,1,1),totflexo(elfd,qvecd,2,2),totflexo(elfd,qvecd,3,3),&
                           totflexo(elfd,qvecd,2,3),totflexo(elfd,qvecd,1,3),totflexo(elfd,qvecd,1,2)
     call wrtout([ab_out,std_out],msg,'COLL')
   end do
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
 real(dp),parameter :: confac=e_Cb/Bohr_meter*1.d9
!arrays
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 real(dp) :: d3cart(2,3,mpert,3,mpert,3,mpert)
 character(len=2) :: voigt(9)=(/'xx','yy','zz','yz','xz','xy','zy','zx','yx'/)

! *********************************************************************

 DBG_ENTER("COLL")

 d3cart(1,:,:,:,:,:,:) = reshape(blkval(1,:),shape = (/3,mpert,3,mpert,3,mpert/))
 d3cart(2,:,:,:,:,:,:) = reshape(blkval(2,:),shape = (/3,mpert,3,mpert,3,mpert/))

!Extraction of the clamped-ion flexoelectric coeficients 
 fac=two/ucvol
 do qvecd=1,3
   do istrs=1,6
     strsd1=alpha(istrs)
     strsd2=beta(istrs)
     strst=natom+3; if (istrs>3) strst=natom+4
     strsd=istrs; if (istrs>3) strsd=istrs-3
     do elfd=1,3
       ciflexo(elfd,qvecd,strsd1,strsd2)=-fac*d3cart(2,elfd,natom+2,strsd,strst,qvecd,natom+8)*confac
       if (istrs>3) ciflexo(elfd,qvecd,strsd2,strsd1)=ciflexo(elfd,qvecd,strsd1,strsd2)
     end do
   end do
 end do

!Print results
 iwrite = ab_out > 0
 if (iwrite) then
   write(msg,'(3a)')ch10,' Type-II electronic (clamped ion) flexoelectric tensor (units= nC/m) ',ch10
   call wrtout([ab_out,std_out],msg,'COLL')
   write(msg,*)'           xx          yy          zz          yz          xz          xy'
   call wrtout([ab_out,std_out],msg,'COLL')
   do ivarA=1,6
     elfd=alpha(ivarA)
     qvecd=beta(ivarA)
     write(msg,'(3x,a2,6f12.6)') voigt(ivarA), ciflexo(elfd,qvecd,1,1),ciflexo(elfd,qvecd,2,2), &
                                 ciflexo(elfd,qvecd,3,3), ciflexo(elfd,qvecd,2,3), &
                                 ciflexo(elfd,qvecd,1,3),ciflexo(elfd,qvecd,1,2)
     call wrtout([ab_out,std_out],msg,'COLL')
   end do
   do ivarA=4,6
     elfd=beta(ivarA)
     qvecd=alpha(ivarA)
     write(msg,'(3x,a2,6f12.6)') voigt(ivarA+3), ciflexo(elfd,qvecd,1,1),ciflexo(elfd,qvecd,2,2), &
                                 ciflexo(elfd,qvecd,3,3), ciflexo(elfd,qvecd,2,3), &
                                 ciflexo(elfd,qvecd,1,3),ciflexo(elfd,qvecd,1,2)
     call wrtout([ab_out,std_out],msg,'COLL')
   end do
 end if

 DBG_EXIT("COLL")

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
!! gprimd(3,3)= basis vectors in the reciprocal space
!! intstrn_only= activates only the calculation of the internal strain tensor
!! mpert =maximum number of ipert
!! natom= number of atoms in unit cell
!! pol1(3,3,3,natom)= tensor with the polarization induced by an atomic displacement (P^(1))
!! rprimd(3,3)= basis vectors in the real space
!! ucvol= unit cell volume
!!
!! OUTPUT
!! mixflexo(3,3,3,3) = type-II mixed contribution to the Flexoelectric Tensor
!! intstrn(3,3,3,natom) = relaxed-ion internal strain tensor
!! psinvdm(3*natom,3*natom) = pseudo inverse of dynamical matrix
!!
!! PARENTS
!!      m_ddb
!!
!! CHILDREN
!!
!! SOURCE

subroutine dtmixflexo(asr,d2asr,blkval1d,blkval2d,blkval,gprimd,intstrn,intstrn_only,mixflexo,mpert,natom,pol1,psinvdm,rprimd,ucvol)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: asr,mpert,natom
 real(dp),intent(in) :: ucvol
 logical,intent(in) :: intstrn_only
!arrays
 real(dp),intent(in) :: d2asr(2,3,mpert,3,mpert)
 real(dp),intent(in) :: blkval1d(2,3,mpert,3,mpert)
 real(dp),intent(in) :: blkval2d(2,3,mpert,3,mpert)
 real(dp),intent(in) :: blkval(2,3*mpert*3*mpert*3*mpert)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(out) :: intstrn(3,3,3,natom)
 real(dp),intent(inout) :: pol1(3,3,3,natom)
 real(dp),intent(out) :: psinvdm(3*natom,3*natom)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(out) :: mixflexo(3,3,3,3)

!Local variables -------------------------
!scalars
 integer :: elfd,iat,iatd,ivar,jat,jatd,jvar,katd,qvecd,qvecd2
 logical :: iwrite
 real(dp),parameter :: confac=e_Cb/Bohr_meter*1.d9
 character(len=500) :: msg
!arrays
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 real(dp) :: d3cart(2,3,mpert,3,mpert,3,mpert)
 real(dp) :: redforces(3,natom),forces(3,natom)
 real(dp) :: phi1(3,natom,3,natom,3)
 real(dp) :: piezofr(3,natom,3,3)
 integer :: flg1(3),flg2(3)
 real(dp) :: vec1(3),vec2(3)
 character(len=2) :: voigt(9)=(/'xx','yy','zz','yz','xz','xy','zy','zx','yx'/)

! *********************************************************************
 
 DBG_ENTER("COLL")

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

!Extraction of the forces and conversion to cartesian coordinates 
!to acount for the improper contribution
 do iat=1,natom
   do iatd=1,3
     redforces(iatd,iat)=-blkval1d(1,iatd,iat,1,1)
   end do
 end do
 forces(:,:)=redforces(:,:)
 flg1(:)=1
 do iat=1,natom
   vec1(:)=forces(:,iat)
   call cart39(flg1,flg2,gprimd,iat,natom,rprimd,vec1,vec2)
   forces(:,iat)=vec2(:)
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

 if (.not.intstrn_only) then
 !Finally calculate the mixed contribution to the FxE tensor
   mixflexo(:,:,:,:)=zero
   do elfd=1,3
     do qvecd=1,3
       do katd=1,3
         do qvecd2=1,3
           do iatd=1,3
             do iat=1,natom
 
               mixflexo(elfd,qvecd,katd,qvecd2)=mixflexo(elfd,qvecd,katd,qvecd2) - &
             pol1(elfd,qvecd,iatd,iat)*intstrn(qvecd2,katd,iatd,iat)*confac
 
            end do
          end do
        end do
      end do
    end do
  end do
 end if

!Print results
 iwrite = ab_out > 0
 if (iwrite) then
   write(msg,'(3a)')ch10,' Force-response internal strain tensor from long-wave magnitudes (units: Hartree/Bohr)',ch10
   call wrtout([ab_out,std_out],msg,'COLL')
   write(msg,*)' atom   dir        xx          yy          zz          yz          xz          xy'
   call wrtout([ab_out,std_out],msg,'COLL')
   do iat=1,natom
     write(msg,'(2x,i3,3x,a3,2x,6f12.6)') iat, 'x', piezofr(1,iat,1,1),piezofr(1,iat,2,2),piezofr(1,iat,3,3),&
                                                    piezofr(1,iat,2,3),piezofr(1,iat,1,3),piezofr(1,iat,1,2)
     call wrtout([ab_out,std_out],msg,'COLL')
     write(msg,'(2x,i3,3x,a3,2x,6f12.6)') iat, 'y', piezofr(2,iat,1,1),piezofr(2,iat,2,2),piezofr(2,iat,3,3),&
                                                    piezofr(2,iat,2,3),piezofr(2,iat,1,3),piezofr(2,iat,1,2)
     call wrtout([ab_out,std_out],msg,'COLL')
     write(msg,'(2x,i3,3x,a3,2x,6f12.6)') iat, 'z', piezofr(3,iat,1,1),piezofr(3,iat,2,2),piezofr(3,iat,3,3),&
                                                    piezofr(3,iat,2,3),piezofr(3,iat,1,3),piezofr(3,iat,1,2)
     call wrtout([ab_out,std_out],msg,'COLL')
   end do 

   write(msg,'(3a)')ch10,' Displacement-response internal strain tensor from long-wave magnitudes (units: Bohr)',ch10
   call wrtout([ab_out,std_out],msg,'COLL')
   write(msg,*)' atom   dir        xx          yy          zz          yz          xz          xy'
   call wrtout([ab_out,std_out],msg,'COLL')
   do iat=1,natom
     write(msg,'(2x,i3,3x,a3,2x,6f12.6)') iat, 'x', intstrn(1,1,1,iat),intstrn(2,2,1,iat),intstrn(3,3,1,iat),&
                                                         intstrn(2,3,1,iat),intstrn(1,3,1,iat),intstrn(1,2,1,iat)
     call wrtout([ab_out,std_out],msg,'COLL')
     write(msg,'(2x,i3,3x,a3,2x,6f12.6)') iat, 'y', intstrn(1,1,2,iat),intstrn(2,2,2,iat),intstrn(3,3,2,iat),&
                                                         intstrn(2,3,2,iat),intstrn(1,3,2,iat),intstrn(1,2,2,iat)
     call wrtout([ab_out,std_out],msg,'COLL')
     write(msg,'(2x,i3,3x,a3,2x,6f12.6)') iat, 'z', intstrn(1,1,3,iat),intstrn(2,2,3,iat),intstrn(3,3,3,iat),&
                                                         intstrn(2,3,3,iat),intstrn(1,3,3,iat),intstrn(1,2,3,iat)
     call wrtout([ab_out,std_out],msg,'COLL')
   end do 

   if (.not.intstrn_only) then
     write(msg,'(3a)')ch10,' Type-II mixed contribution to flexoelectric tensor (units: nC/m)',ch10
     call wrtout([ab_out,std_out],msg,'COLL')
     write(msg,*)'           xx          yy          zz          yz          xz          xy'
     call wrtout([ab_out,std_out],msg,'COLL')
     do ivar=1,6
       elfd=alpha(ivar)
       qvecd=beta(ivar)
       write(msg,'(3x,a2,6f12.6)') voigt(ivar),mixflexo(elfd,qvecd,1,1),mixflexo(elfd,qvecd,2,2),mixflexo(elfd,qvecd,3,3),&
                             mixflexo(elfd,qvecd,2,3),mixflexo(elfd,qvecd,1,3),mixflexo(elfd,qvecd,1,2)
       call wrtout([ab_out,std_out],msg,'COLL')
     end do
     do ivar=4,6
       elfd=beta(ivar)
       qvecd=alpha(ivar)
       write(msg,'(3x,a2,6f12.6)') voigt(ivar+3),mixflexo(elfd,qvecd,1,1),mixflexo(elfd,qvecd,2,2),mixflexo(elfd,qvecd,3,3),&
                             mixflexo(elfd,qvecd,2,3),mixflexo(elfd,qvecd,1,3),mixflexo(elfd,qvecd,1,2)
       call wrtout([ab_out,std_out],msg,'COLL')
     end do
   end if
 end if

 DBG_EXIT("COLL")

 end subroutine dtmixflexo
!!***

!!****f* m_ddb/dtlattflexo
!! NAME
!! dtlattflexo
!!
!! FUNCTION
!! Reads Phi^(1), flexoelectric force response and internal strain tensors 
!! in the Gamma Block coming from the Derivative Data Base
!! (long wave third-order derivatives). And computes the lattice mediated 
!! contribution to the flexoelectric tensor.
!! It also computes and writes the Lagrangian Elastic tensor.
!!
!! INPUTS
!! amu(ntypat)=mass each atom type in the unit cell 
!! blkval1d(2,3,mpert,3,mpert)= 1st derivative wrt stress (at least)
!! blkval2d(2,3,mpert,3,mpert)= 2nd derivatives wrt atom displacements and electric field (at least)
!! blkvalA(2,3*mpert*3*mpert*3*mpert)= matrix of third-order energies for FxE force response tensor
!! blkvalB(2,3*mpert*3*mpert*3*mpert)= matrix of third-order energies for Phi^(1) tensor
!! intstrn(3,3,3,natom) = relaxed-ion internal strain tensor
!! mpert =maximum number of ipert
!! natom= number of atoms in unit cell
!! psinvdm(3*natom,3*natom) = pseudo inverse of dynamical matrix 
!! typat(natom)= Type of each atom in the unit cell
!! ucvol= unit cell volume
!!
!! OUTPUT
!! lattflexo(3,3,3,3) = type-II lattice contribution to the Flexoelectric Tensor
!!
!! PARENTS
!!      m_ddb
!!
!! CHILDREN
!!
!! SOURCE

subroutine dtlattflexo(amu,blkval1d,blkvalA,blkvalB,intstrn,lattflexo,mpert,natom,&
                     & ntypat,psinvdm,typat,ucvol,zeff)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,natom,ntypat
 real(dp),intent(in) :: ucvol

!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: amu(ntypat)
 real(dp),intent(in) :: blkval1d(2,3,mpert,3,mpert)
 real(dp),intent(in) :: blkvalA(2,3*mpert*3*mpert*3*mpert)
 real(dp),intent(in) :: blkvalB(2,3*mpert*3*mpert*3*mpert)
 real(dp),intent(in) :: intstrn(3,3,3,natom)
 real(dp),intent(in) :: psinvdm(3*natom,3*natom)
 real(dp),intent(in) :: zeff(3,3,natom)
 real(dp),intent(out) :: lattflexo(3,3,3,3)

!Local variables -------------------------
!scalars
 integer :: elfd,iat,iatd,istrs,ivar,jat,jatd,jvar,kat,katd,strsd
 integer :: strsd1,strsd2,strst,qvecd,qvecd2
 logical :: iwrite
 real(dp) :: fac,mtot
 real(dp),parameter :: confac=e_Cb/Bohr_meter*1.d9
 character(len=500) :: msg
!arrays
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 real(dp) :: frcelast_t2(3,3,3,3)
 real(dp) :: Csupkap(3,natom,3,3,3)
 real(dp) :: d3cart(2,3,mpert,3,mpert,3,mpert)
 real(dp) :: flexois(3,natom,3,3,3)
 real(dp) :: flexofr(3,natom,3,3,3)
 real(dp) :: hatCsupkap(3,natom,3,3,3)
 real(dp) :: phi1(3,natom,3,natom,3)
 real(dp) :: ricelast_t2(3,3,3,3)
 real(dp) :: roundbkt(3,3,3,3),roundbkt_k(3,3,3,3,natom)
 real(dp) :: sqrbkt_t1(3,3,3,3)
 real(dp) :: stress(3,3)
 character(len=2) :: voigt(9)=(/'xx','yy','zz','yz','xz','xy','zy','zx','yx'/)

! MR: Kept for testing 
! integer :: i,j,k,l
! real(dp) :: delik,deljk,delil,deljl

! *********************************************************************

 DBG_ENTER("COLL")

!Extraction of the stress tensor
 stress(1,1)=blkval1d(1,1,natom+3,1,1)
 stress(2,2)=blkval1d(1,2,natom+3,1,1)
 stress(3,3)=blkval1d(1,3,natom+3,1,1)
 stress(2,3)=blkval1d(1,1,natom+4,1,1);stress(3,2)=stress(2,3)
 stress(1,3)=blkval1d(1,2,natom+4,1,1);stress(3,1)=stress(1,3)
 stress(1,2)=blkval1d(1,3,natom+4,1,1);stress(2,1)=stress(1,2)
 
!Calculate the sublattice-dependent round bracket tensor of PRB 88,174106 (2013)
!First we need to extract the Phi^(1) tensor
 d3cart(1,:,:,:,:,:,:) = reshape(blkvalB(1,:),shape = (/3,mpert,3,mpert,3,mpert/))
 d3cart(2,:,:,:,:,:,:) = reshape(blkvalB(2,:),shape = (/3,mpert,3,mpert,3,mpert/))

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

!Now perform the multiplication with the internal strain
 roundbkt_k(:,:,:,:,:)=zero
 do iat=1,natom
   do iatd=1,3
     do qvecd2=1,3
       do jatd=1,3
         do qvecd=1,3
           do katd=1,3
             do kat=1,natom
               roundbkt_k(iatd,qvecd,jatd,qvecd2,iat)=roundbkt_k(iatd,qvecd,jatd,qvecd2,iat) &
             + phi1(iatd,iat,katd,kat,qvecd)*intstrn(qvecd2,jatd,katd,kat)
             end do
           end do
!         write(100,'(5i3,1x,f12.6)') iat, iatd,qvecd,jatd,qvecd2,roundbkt_k(iatd,qvecd,jatd,qvecd2,iat)
         end do
       end do
     end do
   end do
 end do

!Calculate now the Lagrange elastic tensors 
!First we need to extract the tensor with the q-gradient of the internal strain
!(a.k.a flexoelectric force response tensor) 
 d3cart(1,:,:,:,:,:,:) = reshape(blkvalA(1,:),shape = (/3,mpert,3,mpert,3,mpert/))
 d3cart(2,:,:,:,:,:,:) = reshape(blkvalA(2,:),shape = (/3,mpert,3,mpert,3,mpert/))

 do istrs=1,6
   strsd1=alpha(istrs)
   strsd2=beta(istrs)
   strst=natom+3; if (istrs>3) strst=natom+4
   strsd=istrs; if (istrs>3) strsd=istrs-3
   do qvecd=1,3
     do iat=1,natom
       do iatd=1,3
         flexofr(iatd,iat,qvecd,strsd1,strsd2)=-two*d3cart(1,iatd,iat,strsd,strst,qvecd,natom+8)
         if (istrs>3) flexofr(iatd,iat,qvecd,strsd2,strsd1)=flexofr(iatd,iat,qvecd,strsd1,strsd2)
       end do
     end do
   end do
 end do

!Now compute the type-II frozen-ion elastic tensor (without stress corrected)
 frcelast_t2(:,:,:,:)=zero
 fac=one/ucvol
 do strsd2=1,3
   do strsd1=1,3
     do qvecd=1,3
       do iatd=1,3
         do iat=1,natom
           frcelast_t2(iatd,qvecd,strsd1,strsd2)=frcelast_t2(iatd,qvecd,strsd1,strsd2) + &
        &  flexofr(iatd,iat,qvecd,strsd1,strsd2)*fac
         end do
!         write(100,'(4i3,1x,f12.6)') iatd,qvecd,strsd1,strsd2, frcelast_t2(iatd,qvecd,strsd1,strsd2)
       end do
     end do
   end do
 end do

!Now convert to type-I to obtain the square bracketed tensor of Born and Huang
 do qvecd=1,3
   do strsd2=1,3
     do strsd1=1,3
       do iatd=1,3
         sqrbkt_t1(iatd,strsd1,strsd2,qvecd)=half*(frcelast_t2(iatd,qvecd,strsd1,strsd2) + &
       & frcelast_t2(iatd,strsd2,strsd1,qvecd))
       end do
     end do
   end do
 end do

!Now correct the stress in the square bracketed tesnor tensor
 do qvecd=1,3
   do strsd2=1,3
     do strsd1=1,3
       do iatd=1,3
         if (iatd==strsd1) then
           sqrbkt_t1(iatd,strsd1,strsd2,qvecd)=sqrbkt_t1(iatd,strsd1,strsd2,qvecd) - stress(strsd2,qvecd)
         endif
       end do
     end do
   end do
 end do

!Now convert back to type-II in order to obtain the frozen ion Lagrange elastic tensor.
 frcelast_t2(:,:,:,:)=zero
 do iatd=1,3
   do qvecd=1,3
     do strsd1=1,3
       do strsd2=1,3
         frcelast_t2(iatd,qvecd,strsd1,strsd2)=sqrbkt_t1(iatd,strsd1,qvecd,strsd2) + &
       & sqrbkt_t1(iatd,strsd2,strsd1,qvecd)-sqrbkt_t1(iatd,qvecd,strsd2,strsd1)
       end do
     end do
   end do
 end do

!MR: kept for testing only. If uncommented the resulting elastic tensors must agree with HWRV's ones
! do i=1,3
!   do j=1,3
!     do k=1,3
!       delik=0.d0; if(i==k) delik=1.d0
!       deljk=0.d0; if(j==k) deljk=1.d0
!       do l=1,3
!         delil=0.d0; if(i==l) delil=1.d0
!         deljl=0.d0; if(j==l) deljl=1.d0

!         frcelast_t2(i,j,k,l)=frcelast_t2(i,j,k,l) + &
!       & 0.5d0*( deljk*stress(i,l) + delik*stress(j,l) + delil*stress(j,k) &
!       & + deljl*stress(i,k) )

!       end do
!     end do
!   end do
! end do

!Now compute the rount bracketed tensor of Born and Huang 
!and sum with the clamped ion elastic tensor to obtain the relaxed ion one
 roundbkt(:,:,:,:)=zero
 do strsd2=1,3
   do strsd1=1,3
     do qvecd=1,3
       do iatd=1,3
         do iat=1,natom
           roundbkt(iatd,qvecd,strsd1,strsd2)=roundbkt(iatd,qvecd,strsd1,strsd2) + &
         & roundbkt_k(iatd,qvecd,strsd1,strsd2,iat)*fac
         end do
         ricelast_t2(iatd,qvecd,strsd1,strsd2)=frcelast_t2(iatd,qvecd,strsd1,strsd2) + &
       & roundbkt(iatd,qvecd,strsd1,strsd2)
       end do
     end do
   end do
 end do

!In last place compute the lattice contribution to the FxE tensor
!First obtain the C^{\kappa} tensor of Eq. 59 of PRB 88,174106 (2013) 
 do strsd2=1,3
   do strsd1=1,3
     do qvecd=1,3
       do iatd=1,3
         do iat=1,natom
           Csupkap(iatd,iat,qvecd,strsd1,strsd2)=flexofr(iatd,iat,qvecd,strsd1,strsd2) + &
         & roundbkt_k(iatd,qvecd,strsd1,strsd2,iat) 
         end do
       end do
     end do
   end do
 end do

!Then separate the mass-dependent part using the ion relaxed Lagrange elastic tensor
 mtot=zero
 do iat=1,natom
   mtot=mtot + amu(typat(iat))
 end do

 do strsd2=1,3
   do strsd1=1,3
     do qvecd=1,3
       do iatd=1,3
         do iat=1,natom
           hatCsupkap(iatd,iat,qvecd,strsd1,strsd2)=Csupkap(iatd,iat,qvecd,strsd1,strsd2) - &
         & amu(typat(iat))/mtot*ucvol*ricelast_t2(iatd,qvecd,strsd1,strsd2)
         end do
       end do
     end do
   end do
 end do

!Now compute the type-II flexoelectric internal strain tensor
 flexois(:,:,:,:,:)=zero
 do strsd2=1,3
   do strsd1=1,3
     do qvecd=1,3
       do iat=1,natom
         do iatd=1,3
           ivar=(iat-1)*3+iatd
           do jat=1,natom
             do jatd=1,3
               jvar=(jat-1)*3+jatd
               flexois(iatd,iat,qvecd,strsd1,strsd2)=flexois(iatd,iat,qvecd,strsd1,strsd2) + &
             & psinvdm(ivar,jvar)*hatCsupkap(jatd,jat,qvecd,strsd1,strsd2)
             end do
           end do
         end do
       end do
     end do
   end do
 end do
             
!Finally multiply by the effective charges to obtain the FxE tensor
 lattflexo(:,:,:,:)=zero
 do strsd2=1,3
   do strsd1=1,3
     do qvecd=1,3
       do elfd=1,3
         do iat=1,natom
           do iatd=1,3
             lattflexo(elfd,qvecd,strsd1,strsd2)=lattflexo(elfd,qvecd,strsd1,strsd2) + &
           & zeff(elfd,iatd,iat)*flexois(iatd,iat,qvecd,strsd1,strsd2)/ucvol*confac
           end do
         end do
       end do
     end do
   end do
 end do


 iwrite = ab_out > 0
 if (iwrite) then
   write(msg,'(3a)')ch10,' Lagrange elastic tensor from long wave magnitudes (clamped ion) (units= 10^2 GPa) ',ch10
   call wrtout([ab_out,std_out],msg,'COLL')
   write(msg,*)'      xx          yy          zz          yz          xz          xy'
   call wrtout([ab_out,std_out],msg,'COLL')
   frcelast_t2(:,:,:,:)=frcelast_t2(:,:,:,:)*HaBohr3_GPa/100.00_dp
   do ivar=1,6
     iatd=alpha(ivar)
     qvecd=beta(ivar)
     write(msg,'(9f12.6)') frcelast_t2(iatd,qvecd,1,1), frcelast_t2(iatd,qvecd,2,2),frcelast_t2(iatd,qvecd,3,3), &
                           frcelast_t2(iatd,qvecd,2,3), frcelast_t2(iatd,qvecd,1,3),frcelast_t2(iatd,qvecd,1,2)

     call wrtout([ab_out,std_out],msg,'COLL')
   end do

   write(msg,'(3a)')ch10,' Lagrange elastic tensor from long wave magnitudes (relaxed ion) (units= 10^2 GPa) ',ch10
   call wrtout([ab_out,std_out],msg,'COLL')
   write(msg,*)'      xx          yy          zz          yz          xz          xy'
   call wrtout([ab_out,std_out],msg,'COLL')
   ricelast_t2(:,:,:,:)=ricelast_t2(:,:,:,:)*HaBohr3_GPa/100.00_dp
   do ivar=1,6
     iatd=alpha(ivar)
     qvecd=beta(ivar)
     write(msg,'(9f12.6)') ricelast_t2(iatd,qvecd,1,1), ricelast_t2(iatd,qvecd,2,2),ricelast_t2(iatd,qvecd,3,3), &
                           ricelast_t2(iatd,qvecd,2,3), ricelast_t2(iatd,qvecd,1,3),ricelast_t2(iatd,qvecd,1,2)

     call wrtout([ab_out,std_out],msg,'COLL')
   end do

   write(msg,'(3a)')ch10,' Displacement-response flexoelectric force response tensor (units: eV)',ch10
   call wrtout([ab_out,std_out],msg,'COLL')
   write(msg,*)' atom   dir        xx           yy           zz           yz           xz           xy'
   call wrtout([ab_out,std_out],msg,'COLL')
   Csupkap(:,:,:,:,:)=Csupkap(:,:,:,:,:)*Ha_eV
   do iat=1,natom
     write(msg,'(2x,i3,3x,a3,2x,6(f12.6,1x))') iat, 'xx', Csupkap(1,iat,1,1,1),Csupkap(1,iat,1,2,2),Csupkap(1,iat,1,3,3),&
                                                   & Csupkap(1,iat,1,2,3),Csupkap(1,iat,1,1,3),Csupkap(1,iat,1,1,2)  
     call wrtout([ab_out,std_out],msg,'COLL')
     write(msg,'(2x,i3,3x,a3,2x,6(f12.6,1x))') iat, 'yy', Csupkap(2,iat,2,1,1),Csupkap(2,iat,2,2,2),Csupkap(2,iat,2,3,3),&
                                                   & Csupkap(2,iat,2,2,3),Csupkap(2,iat,2,1,3),Csupkap(2,iat,2,1,2)  
     call wrtout([ab_out,std_out],msg,'COLL')
     write(msg,'(2x,i3,3x,a3,2x,6(f12.6,1x))') iat, 'zz', Csupkap(3,iat,3,1,1),Csupkap(3,iat,3,2,2),Csupkap(3,iat,3,3,3),&
                                                   & Csupkap(3,iat,3,2,3),Csupkap(3,iat,3,1,3),Csupkap(3,iat,3,1,2)  
     call wrtout([ab_out,std_out],msg,'COLL')
     write(msg,'(2x,i3,3x,a3,2x,6(f12.6,1x))') iat, 'yz', Csupkap(2,iat,3,1,1),Csupkap(2,iat,3,2,2),Csupkap(2,iat,3,3,3),&
                                                   & Csupkap(2,iat,3,2,3),Csupkap(2,iat,3,1,3),Csupkap(2,iat,3,1,2)  
     call wrtout([ab_out,std_out],msg,'COLL')
     write(msg,'(2x,i3,3x,a3,2x,6(f12.6,1x))') iat, 'xz', Csupkap(1,iat,3,1,1),Csupkap(1,iat,3,2,2),Csupkap(1,iat,3,3,3),&
                                                   & Csupkap(1,iat,3,2,3),Csupkap(1,iat,3,1,3),Csupkap(1,iat,3,1,2)  
     call wrtout([ab_out,std_out],msg,'COLL')
     write(msg,'(2x,i3,3x,a3,2x,6(f12.6,1x))') iat, 'xy', Csupkap(1,iat,2,1,1),Csupkap(1,iat,2,2,2),Csupkap(1,iat,2,3,3),&
                                                   & Csupkap(1,iat,2,2,3),Csupkap(1,iat,2,1,3),Csupkap(1,iat,2,1,2)  
     call wrtout([ab_out,std_out],msg,'COLL')
   end do


   write(msg,'(3a)')ch10,' Displacement-response flexoelectric internal strain tensor (units: Bohr^2)',ch10
   call wrtout([ab_out,std_out],msg,'COLL')
   write(msg,*)' atom   dir        xx           yy           zz           yz           xz           xy'
   call wrtout([ab_out,std_out],msg,'COLL')
   do iat=1,natom
     write(msg,'(2x,i3,3x,a3,2x,6(f12.6,1x))') iat, 'xx', flexois(1,iat,1,1,1),flexois(1,iat,1,2,2),flexois(1,iat,1,3,3),&
                                                   & flexois(1,iat,1,2,3),flexois(1,iat,1,1,3),flexois(1,iat,1,1,2)  
     call wrtout([ab_out,std_out],msg,'COLL')
     write(msg,'(2x,i3,3x,a3,2x,6(f12.6,1x))') iat, 'yy', flexois(2,iat,2,1,1),flexois(2,iat,2,2,2),flexois(2,iat,2,3,3),&
                                                   & flexois(2,iat,2,2,3),flexois(2,iat,2,1,3),flexois(2,iat,2,1,2)  
     call wrtout([ab_out,std_out],msg,'COLL')
     write(msg,'(2x,i3,3x,a3,2x,6(f12.6,1x))') iat, 'zz', flexois(3,iat,3,1,1),flexois(3,iat,3,2,2),flexois(3,iat,3,3,3),&
                                                   & flexois(3,iat,3,2,3),flexois(3,iat,3,1,3),flexois(3,iat,3,1,2)  
     call wrtout([ab_out,std_out],msg,'COLL')
     write(msg,'(2x,i3,3x,a3,2x,6(f12.6,1x))') iat, 'yz', flexois(2,iat,3,1,1),flexois(2,iat,3,2,2),flexois(2,iat,3,3,3),&
                                                   & flexois(2,iat,3,2,3),flexois(2,iat,3,1,3),flexois(2,iat,3,1,2)  
     call wrtout([ab_out,std_out],msg,'COLL')
     write(msg,'(2x,i3,3x,a3,2x,6(f12.6,1x))') iat, 'xz', flexois(1,iat,3,1,1),flexois(1,iat,3,2,2),flexois(1,iat,3,3,3),&
                                                   & flexois(1,iat,3,2,3),flexois(1,iat,3,1,3),flexois(1,iat,3,1,2)  
     call wrtout([ab_out,std_out],msg,'COLL')
     write(msg,'(2x,i3,3x,a3,2x,6(f12.6,1x))') iat, 'xy', flexois(1,iat,2,1,1),flexois(1,iat,2,2,2),flexois(1,iat,2,3,3),&
                                                   & flexois(1,iat,2,2,3),flexois(1,iat,2,1,3),flexois(1,iat,2,1,2)  
     call wrtout([ab_out,std_out],msg,'COLL')
   end do

   write(msg,'(3a)')ch10,' Type-II lattice contribution to flexoelectric tensor (units= nC/m) ',ch10
   call wrtout([ab_out,std_out],msg,'COLL')
   write(msg,*)'           xx          yy          zz          yz          xz          xy'
   call wrtout([ab_out,std_out],msg,'COLL')
   do ivar=1,6
     elfd=alpha(ivar)
     qvecd=beta(ivar)
     write(msg,'(3x,a2,6f12.6)') voigt(ivar),lattflexo(elfd,qvecd,1,1),lattflexo(elfd,qvecd,2,2),lattflexo(elfd,qvecd,3,3),&
                           lattflexo(elfd,qvecd,2,3),lattflexo(elfd,qvecd,1,3),lattflexo(elfd,qvecd,1,2)
     call wrtout([ab_out,std_out],msg,'COLL')
   end do
   do ivar=4,6
     elfd=beta(ivar)
     qvecd=alpha(ivar)
     write(msg,'(3x,a2,6f12.6)') voigt(ivar+3),lattflexo(elfd,qvecd,1,1),lattflexo(elfd,qvecd,2,2),lattflexo(elfd,qvecd,3,3),&
                           lattflexo(elfd,qvecd,2,3),lattflexo(elfd,qvecd,1,3),lattflexo(elfd,qvecd,1,2)
     call wrtout([ab_out,std_out],msg,'COLL')
   end do

 end if
 DBG_EXIT("COLL")

 end subroutine dtlattflexo
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

 DBG_ENTER("COLL")

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

 DBG_EXIT("COLL")

 end subroutine dm_psinv
!!***
end module m_ddb_flexo
!!***
