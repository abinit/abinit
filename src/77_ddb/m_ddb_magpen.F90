!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ddb_magpen
!! NAME
!!  m_ddb_magpen
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2023 ABINIT group (FIXME: add author)
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

module m_ddb_magpen
    
 use defs_basis
 use m_abicore
 use m_profiling_abi
 use m_errors
 use m_ddb
 use m_fstrings,        only : itoa, sjoin

 implicit none

 public :: ddb_magpen       ! Convert the derivatives calculated with the magnetic penalty into the physically relevant ones.

 private

! *************************************************************************

contains 
!!***

!!****f* m_ddb_magpen/ddb_magpen
!! NAME
!! ddb_magpen
!!
!! FUNCTION
!! Convert the second-order derivatives calculated with the magnetic penalty
!! into physically relevant quantities.
!!
!! INPUTS
!! ddb (INOUT) = ddb block datastructure
!! magpen = amplitude (in Ha) of the applied magnetic penalty 
!! mpatpol(2) = Atoms on which the magnetic penalty has been applied
!! mpdir(3) = Directions along which the spin-degrees of freedom have been stiffened 
!! mpert = maximum number of ipert
!! natom= number of atoms in unit cell
!! ntypat= number of atom types
!!
!! OUTPUT
!! ddb_lw= ddb block datastructure
!!
!! SOURCE

 subroutine ddb_magpen(ddb,magpen,mpatpol,mpdir,mpert,natom,ntypat,rftyp)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,natom,ntypat,rftyp
 real(dp),intent(in) :: magpen
!arrays
 type(ddb_type),intent(inout) :: ddb
 integer,intent(in) :: mpatpol(2),mpdir(3)

!Local variables -------------------------
!scalars
 integer :: iblok,ii,jblok,nblok,nsize,cnt,nmat,nmdir
 character(len=500) :: msg
!arrays
 integer :: rfelfd(4),rfphon(4),rfstrs(4),rfmagn(4)
 real(dp) :: qphnrm(3),qphon(3,3)
 real(dp), allocatable :: magsus(:,:,:,:,:)

! *********************************************************************
 if (magpen<zero) then
   nmat= 1
 else if (magpen>zero) then
   nmat= mpatpol(2) - mpatpol(1) + 1
 end if
 nmdir=sum(mpdir(:))
 ABI_MALLOC(magsus,(2,nmdir,nmat,nmdir,nmat))

 do iblok=1,ddb%nblok

   ! Look for the spin-susceptibility block in the DDB
   qphon=zero
   qphon(:,1)=ddb%qpt(:,iblok)
   qphnrm(1)=ddb%nrm(1,iblok)
   rfphon(1:2)=0
   rfelfd(1:2)=0
   rfstrs(1:2)=0
   if (magpen<zero) then
     rfmagn(1:2)= 1
   else if (magpen>zero) then
     rfmagn(1:2)= 2
   end if

   call ddb%get_block(jblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, &
  & mpatpol=mpatpol,mpdir=mpdir,rfmagn=rfmagn)

   ! Calculate and write the spin-susceptibility matrices
   if (jblok /= 0) then
     if (magpen<zero) then
       write(msg, '(2a,(80a),4a)' ) ch10,('=',ii=1,80),ch10,ch10,&
       ' Spin susceptibility (Uniform Zeeman) ',ch10
       call wrtout([std_out, ab_out], msg)
     else if (magpen>zero) then
       write(msg, '(2a,(80a),4a)' ) ch10,('=',ii=1,80),ch10,ch10,&
       ' Spin susceptibility (Local Zeeman) ',ch10
       call wrtout([std_out, ab_out], msg)
     end if

     call spinsus(ddb%val,iblok,magpen,magsus,mpatpol,mpdir,mpert,natom,nblok,nmat,nmdir,unit=dev_null)

   end if

 end do

!Deallocations
 ABI_FREE(magsus)

 end subroutine ddb_magpen
!!***

!!****f* m_ddb_magpen/spinsus
!! NAME
!! spinsus
!!
!! FUNCTION
!! Calculate the spin-susceptibility matrix and its inverse
!!
!! INPUTS
!! blkval(2,3*mpert*3*mpert,nblok)=  Second-order derivative matrices
!!  In our case, the nblok is restricted to iblok
!! iblok= index of the current block
!! magpen = amplitude (in Ha) of the applied magnetic penalty 
!! mpatpol(2) = Atoms on which the magnetic penalty has been applied
!! mpdir(3) = Directions along which the spin-degrees of freedom have been stiffened 
!! mpert =maximum number of ipert
!! natom= number of atoms in unit cell
!! nblok= number of blocks in the DDB
!! nmat= number of atoms on which the magnetic penalty was applied
!!  (=1 if uniform Zeeman case)
!! nmdir= number of directions along which the magnetic penalty was applied
!! [unit]=Output unit number
!!
!! OUTPUT
!! magsus(2,nmdir,nmat,nmdir,nmat)= Spin-sussceptibility tensor
!!
!! SOURCE

 subroutine spinsus(blkval,iblok,magpen,magsus,mpatpol,mpdir,mpert,natom,nblok,nmat,nmdir,unit)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,mpert,natom,nblok,nmat,nmdir
 real(dp),intent(in) :: magpen
 integer,intent(in),optional :: unit
!arrays
 integer,intent(in) :: mpatpol(2),mpdir(3)
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
 real(dp),intent(out) :: magsus(2,nmdir,nmat,nmdir,nmat)

!Local variables -------------------------
!scalars
 integer :: iat1,iat2,icol,idir1,idir2,info,ipert1,ipert2,irow,ndim,unt
 integer :: ipert1_red,ipert2_red,idir1_red,idir2_red
!arrays
 real(dp) :: barmagsus(2,nmdir,nmat,nmdir,nmat)
 real(dp) :: idty(2,nmdir,nmat,nmdir,nmat)
 complex(dp),allocatable :: work1(:,:),work2(:,:)

! *********************************************************************

 unt = std_out; if (present(unit)) unt = unit

 !Extract the penalized susceptibility
 idty=zero
 ipert2_red= 0
 do iat2= mpatpol(1), mpatpol(2)
   ipert2= natom + 11 + iat2
   ipert2_red= ipert2_red + 1
   idir2_red= 0
   do idir2= 1, 3
     if (mpdir(idir2)==0) cycle
     idir2_red= idir2_red + 1
     idty(1,idir2_red,ipert2_red,idir2_red,ipert2_red)=one
     ipert1_red=0
     do iat1= mpatpol(1), mpatpol(2)
       ipert1= natom + 11 + iat1
       ipert1_red= ipert1_red + 1
       idir1_red= 0
       do idir1= 1, 3
         if (mpdir(idir1)==0) cycle
         idir1_red=idir1_red+1

         !TODO: the two factor needs to be applied in ABINIT when
         !passing the magnetic moments to d2etot
         barmagsus(:,idir1_red,ipert1_red,idir2_red,ipert2_red)= &
       & two*blkval(:,idir1,ipert1,idir2,ipert2,iblok)

       end do
     end do
   end do
 end do

 !Use magsus to store the intermediate array
 magsus=idty-magpen*barmagsus
     
 !Invert the arrays
 write(*,*) nmat, nmdir
 ndim=nmat*nmdir
 ABI_MALLOC(work1,(ndim,ndim))
 ABI_MALLOC(work2,(ndim,ndim))
 work1= zero; work2= zero
 do iat1=1,nmat
   do idir1=1,nmdir
     irow=idir1+(iat1-1)*nmdir
     iat2=iat1
     do idir2=idir1,nmdir
       icol=idir2+(iat2-1)*nmdir
       work1(irow,icol)=cmplx(barmagsus(1,idir1,iat1,idir2,iat2),barmagsus(2,idir1,iat1,idir2,iat2),16)
       work2(irow,icol)=cmplx(magsus(1,idir1,iat1,idir2,iat2),magsus(2,idir1,iat1,idir2,iat2),16)
     end do
     do iat2=iat1+1,nmat
       do idir2=1,nmdir
         icol=idir2+(iat2-1)*nmdir
         work1(irow,icol)=cmplx(barmagsus(1,idir1,iat1,idir2,iat2),barmagsus(2,idir1,iat1,idir2,iat2),16)
         work2(irow,icol)=cmplx(magsus(1,idir1,iat1,idir2,iat2),magsus(2,idir1,iat1,idir2,iat2),16)
       end do
     end do
   end do
 end do

 call zpotrf( 'U', ndim, work1, ndim, info )
 ABI_CHECK(info == 0, sjoin('zpotrf returned:', itoa(info)))

 call zpotri( 'U', ndim, work1, ndim, info )
 ABI_CHECK(info == 0, sjoin('zpotri returned:', itoa(info)))

 call zpotrf( 'U', ndim, work2, ndim, info )
 ABI_CHECK(info == 0, sjoin('zpotrf returned:', itoa(info)))

 call zpotri( 'U', ndim, work2, ndim, info )
 ABI_CHECK(info == 0, sjoin('zpotri returned:', itoa(info)))

 ABI_FREE(work1)
 ABI_FREE(work2)

 end subroutine spinsus
!!***

end module m_ddb_magpen
!!***
