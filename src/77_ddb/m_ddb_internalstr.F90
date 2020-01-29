!!****m* ABINIT/m_ddb_internalstr
!! NAME
!!  m_ddb_internalstr
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1999-2020 ABINIT group (XW)
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

module m_ddb_internalstr

 use defs_basis
 use m_abicore
 use m_errors
 use m_crystal
 use m_ddb

 use m_fstrings,     only : itoa, sjoin
 use m_dynmat,       only : asria_corr

 implicit none

 private
!!***

 public :: ddb_internalstr
!!***

contains
!!***

!!****f* ABINIT/ddb_internalstr
!!
!! NAME
!! ddb_internalstr
!!
!! FUNCTION
!! Get the insternal strain tensors,both force response and displacement response ones.
!!
!! INPUTS
!! asrq0<asrq0_t>=Object for the treatment of the ASR based on the q=0 block found in the DDB file.
!! blkval(2,3,mpert,3,mpert,nblok)=
!!   second derivatives of total energy with respect to electric fields
!!   atom displacements,strain,...... all in cartesian coordinates
!! crystal<crystal_t>=Crystalline structure info.
!! iblok= blok number in DDB file
!! iout=out file number
!! mpert=maximum number of ipert
!! msize=Maximum size of dynamical matrices and other perturbations (ddk, dde...)
!! natom=number of atoms in unit cell
!! nblok=number of total bloks in DDB file
!! prt_internalstr=if 2 or higher, print force and displacement internal strain,
!!                 if 1, print only force internal strain,
!!                 if 0, do not print internal strain.
!!
!! OUTPUT
!! instrain=force response internal strain tensor
!!
!! NOTES
!! In output of internal strain tensor,column runs from strain1 to
!! strain6(in Voigt notation),row runs from atom1x,atom1y,atom1z,atom2x,.......
!! sum rule is applied on the internal strain tensor
!!
!! PARENTS
!!      anaddb,m_effective_potential_file
!!
!! CHILDREN
!!      asria_corr,wrtout,zhpev
!!
!! SOURCE

subroutine ddb_internalstr(asr,&
!&crystal,&
& blkval,&
!&asrq0,&
& d2asr,iblok,instrain,iout,mpert,&
!&msize,&
natom,nblok,prt_internalstr)

!Arguments----------------------------------------------
!scalars
 integer,intent(in) :: asr,iblok,iout,mpert,natom,nblok,prt_internalstr
!integer,intent(in) :: msize
!type(crystal_t),intent(in) :: crystal
!type(asrq0_t),intent(inout) :: asrq0
!arrays
 real(dp),intent(in) :: d2asr(2,3,natom,3,natom)
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
 real(dp),intent(out) :: instrain(3*natom,6)

!Local variables------------------------------------
!scalars
 integer :: idirA,idirB,ier,ii1,ii2,ipertA,ipertB,ivarA,ivarB
 character(len=500) :: direction,message
!arrays
 real(dp) :: Amatr(3*natom-3,3*natom-3),Apmatr(3*natom,3*natom)
 real(dp) :: Bmatr(2,((3*natom-3)*(3*natom-2))/2)
 real(dp) :: Bpmatr(2,(3*natom*(3*natom+1))/2),Cmatr(3*natom-3,3*natom-3)
 real(dp) :: Cpmatr(3*natom,3*natom),Nmatr(3*natom,3*natom),deviation(3,6)
 real(dp) :: eigval(3*natom-3),eigvalp(3*natom),eigvec(2,3*natom-3,3*natom-3)
 real(dp) :: eigvecp(2,3*natom,3*natom),instrain_dis(6,3*natom)
 real(dp) :: kmatrix(3*natom,3*natom),zhpev1(2,2*3*natom-4)
 real(dp) :: zhpev1p(2,2*3*natom-1),zhpev2(3*3*natom-5),zhpev2p(3*3*natom-2)
 real(dp) :: d2cart(2,3*natom,3*natom)

!***************************************************************

!extract internal strain from DDB matrix data

 do ipertA=1,natom
   do idirA=1,3
     ivarA=idirA+3*(ipertA-1)
     do ivarB=1,6
       if(ivarB<=3) then
         idirB=ivarB
         ipertB=natom+3
!        for the diagonal modulus
       else if(ivarB>3) then
         idirB=ivarB-3
         ipertB=natom+4
!        for the shear modulus
       end if
       instrain(ivarA,ivarB)=(-1.0_dp)*blkval(1,idirA,ipertA,idirB,ipertB,iblok)
!      write(std_out,'(es16.6)')blkval(1,idirA,ipertA,idirB,ipertB,iblok)
     end do
   end do
 end do
!according to the definition there is a minus sign before the second derivative

!apply sum rule to the internal strain tensor
 deviation(:,:)=zero
 do ivarB=1,6
   do ivarA=1,3*natom
     if(mod(ivarA,3)==0)then
       deviation(1,ivarB)=deviation(1,ivarB)+instrain(ivarA,ivarB)
     end if
     if(mod(ivarA,3)==1)then
       deviation(2,ivarB)=deviation(2,ivarB)+instrain(ivarA,ivarB)
     end if
     if(mod(ivarA,3)==2)then
       deviation(3,ivarB)=deviation(3,ivarB)+instrain(ivarA,ivarB)
     end if
   end do
 end do

 do ivarB=1,6
   do ivarA=1,3*natom
     if(mod(ivarA,3)==0)then
       instrain(ivarA,ivarB)=instrain(ivarA,ivarB)-deviation(1,ivarB)/natom
     end if
     if(mod(ivarA,3)==1)then
       instrain(ivarA,ivarB)=instrain(ivarA,ivarB)-deviation(2,ivarB)/natom
     end if
     if(mod(ivarA,3)==2)then
       instrain(ivarA,ivarB)=instrain(ivarA,ivarB)-deviation(3,ivarB)/natom
     end if
   end do
 end do
!ending the sum rule

!print the force response internal strain constants into the output file
 if(prt_internalstr>0)then
   write(message,'(a,a,a,a)')ch10,&
&   ' Force-response internal strain tensor','(Unit:Hartree/bohr)',ch10
   call wrtout(std_out,message,'COLL')
   call wrtout(iout,message,'COLL')

   write(message,'(a5,a4,a11,a12,a12,a12,a12,a12)')' Atom',' dir','strainxx',&
&   'strainyy','strainzz','strainyz','strainxz','strainxy'
   call wrtout(std_out,message,'COLL')
   do ii1=1,3*natom
     if(mod(ii1,3)==1)then
       direction='x'
     elseif(mod(ii1,3)==2)then
       direction='y'
     elseif(mod(ii1,3)==0)then
       direction='z'
     end if
     write(message,'(a1,i2,a2,a3,6f12.7)')' ',int((ii1-1)/3)+1,'  ',direction,&
&     instrain(ii1,1),instrain(ii1,2),instrain(ii1,3),&
&     instrain(ii1,4),instrain(ii1,5),instrain(ii1,6)
     call wrtout(std_out,message,'COLL')
   end do
 endif

!now write into the ddb output file
 write(message,'(a5,a4,a11,a12,a12,a12,a12,a12)')' Atom',' dir','strainxx',&
& 'strainyy','strainzz','strainyz','strainxz','strainxy'
 call wrtout(iout,message,'COLL')
 do ii1=1,3*natom
   if(mod(ii1,3)==1)then
     direction='x'
   elseif(mod(ii1,3)==2)then
     direction='y'
   elseif(mod(ii1,3)==0)then
     direction='z'
   end if
   write(message,'(a1,i2,a2,a3,6f12.7)')' ',int((ii1-1)/3)+1,'  ',direction,&
&   instrain(ii1,1),&
&   instrain(ii1,2),instrain(ii1,3),&
&   instrain(ii1,4),instrain(ii1,5),instrain(ii1,6)
   call wrtout(iout,message,'COLL')
 end do

! ----------------------------------------------------------------------------------------

!Try to get the displacement response internal strain tensor
!first need the inverse of force constant matrix
 d2cart = zero
 do ipertA=1,natom
   do ii1=1,3
     ivarA=ii1+3*(ipertA-1)
     do ipertB=1,natom
       do ii2=1,3
         ivarB=ii2+3*(ipertB-1)
         d2cart(1,ivarA,ivarB)=blkval(1,ii1,ipertA,ii2,ipertB,iblok)
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

!Now do simple mulplication to obtain the displacement response
!internal strain tensor
 instrain_dis(:,:)=0.0_dp
 do ivarA=1,6
   do ivarB=1,3*natom
     do ii1=1,3*natom
       instrain_dis(ivarA,ivarB)=instrain_dis(ivarA,ivarB)+&
&       instrain(ii1,ivarA)*kmatrix(ii1,ivarB)
     end do
   end do
 end do

!Print out the results
 if(prt_internalstr>1)then
   write(message,'(a,a,a,a)')ch10,&
&   ' Displacement-response internal strain ', 'tensor (Unit:Bohr)',ch10
   call wrtout(std_out,message,'COLL')
   call wrtout(iout,message,'COLL')
   write(message,'(a5,a4,a11,a12,a12,a12,a12,a12)')' Atom',' dir','strainxx',&
&   'strainyy','strainzz','strainyz','strainxz','strainxy'
   call wrtout(std_out,message,'COLL')
   call wrtout(iout,message,'COLL')
   do ivarA=1,3*natom
     if(mod(ivarA,3)==1)then
       direction='x'
     elseif(mod(ivarA,3)==2)then
       direction='y'
     elseif(mod(ivarA,3)==0)then
       direction='z'
     end if
     write(message,'(a1,i2,a2,a3,6f12.7)')' ',int((ivarA-1)/3)+1,'  ',direction,&
&     instrain_dis(1,ivarA),instrain_dis(2,ivarA),&
&     instrain_dis(3,ivarA),instrain_dis(4,ivarA),instrain_dis(5,ivarA),&
&     instrain_dis(6,ivarA)
     call wrtout(std_out,message,'COLL')
     call wrtout(iout,message,'COLL')
   end do
 endif

end subroutine ddb_internalstr
!!***

end module m_ddb_internalstr
!!***
