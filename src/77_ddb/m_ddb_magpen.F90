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
!! ddb= ddb%val updated with the corrected second-order derivatives
!!
!! SOURCE

 subroutine ddb_magpen(ddb,magpen,mpatpol,mpdir,mpert,natom,ntypat,prtvol,rftyp)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,natom,ntypat,prtvol,rftyp
 real(dp),intent(in) :: magpen
!arrays
 type(ddb_type),intent(inout) :: ddb
 integer,intent(in) :: mpatpol(2),mpdir(3)

!Local variables -------------------------
!scalars
 integer :: iblok,ii,jblok,kblok,nblok,nsize,ndim 
 integer :: nmat,nmdir
 character(len=500) :: msg
!arrays
 integer :: rfelfd(4),rfphon(4),rfstrs(4),rfmagn(4)
 real(dp) :: qphnrm(3),qphon(3,3)
 complex(dp), allocatable :: invbarmagsus(:,:), magsus(:,:)
 complex(dp), allocatable :: mmom(:,:), zfield(:,:)

! *********************************************************************
 if (magpen<zero) then
   nmat= 1
 else if (magpen>zero) then
   nmat= mpatpol(2) - mpatpol(1) + 1
 end if
 nmdir=sum(mpdir(:))
 ndim=nmat*nmdir
 ABI_MALLOC(magsus,(ndim,ndim))
 ABI_MALLOC(invbarmagsus,(ndim,ndim))
 ABI_MALLOC(mmom,(ndim,(natom+2)*3))
 ABI_MALLOC(zfield,(ndim,(natom+2)*3))

 do kblok=1,ddb%nblok

   ! Look for the spin-susceptibility block in the DDB
   qphon=zero
   qphon(:,1)=ddb%qpt(:,kblok)
   qphnrm(1)=ddb%nrm(1,kblok)
   rfphon(1:2)=0
   rfelfd(1:2)=0
   rfstrs(1:2)=0
   if (magpen<zero) then
     rfmagn(1:2)= 1
   else if (magpen>zero) then
     rfmagn(1:2)= 2
   end if

   call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, &
  & mpatpol=mpatpol,mpdir=mpdir,rfmagn=rfmagn)

   ! Calculate and write the spin-susceptibility matrices
   if (iblok /= 0) then
     if (magpen<zero) then
       write(msg, '(2a,(80a),4a)' ) ch10,('=',ii=1,80),ch10,ch10,&
       ' Spin susceptibility (Uniform Zeeman) ',ch10
       call wrtout([std_out, ab_out], msg)
     else if (magpen>zero) then
       write(msg, '(2a,(80a),4a)' ) ch10,('=',ii=1,80),ch10,ch10,&
       ' Spin susceptibility (Local Zeeman) ',ch10
       call wrtout([std_out, ab_out], msg)
     end if

     call spinsus(ddb%val,iblok,invbarmagsus,magpen,magsus,&
   & mpatpol,mpdir,mpert,natom,nblok,ndim,nmat,nmdir,prtvol)

   end if


   ! Calculate and write the induced magnetic moments 
   write(msg, '(2a,(80a),4a)' ) ch10,('=',ii=1,80),ch10,ch10,&
   ' First-order magnetic moments ',ch10
   call wrtout([std_out, ab_out], msg)

   ! First atomic-displacement
   ! Look for the induced magnetic moments block in the DDB
   ! TODO: here we are reading the transpose conjugate, 
   ! this has been de
   rfphon(2)=1
   rfelfd(1:2)=0
   rfstrs(1:2)=0
   rfmagn(:)=0
   if (magpen<zero) then
     rfmagn(1)= 1
   else if (magpen>zero) then
     rfmagn(1)= 2
   end if

   call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, &
  & mpatpol=mpatpol,mpdir=mpdir,rfmagn=rfmagn)

   ! Then electric field
   ! Look for the induced magnetic moments block in the DDB
   rfphon(:)=0
   rfelfd(2)=2
   rfstrs(1:2)=0
   rfmagn(:)=0
   if (magpen<zero) then
     rfmagn(1)= 1
   else if (magpen>zero) then
     rfmagn(1)= 2
   end if

   call ddb%get_block(jblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, &
  & mpatpol=mpatpol,mpdir=mpdir,rfmagn=rfmagn)

   if (iblok /= 0 .or. jblok /=0) then
     call magmom(ddb%val,invbarmagsus,iblok,jblok,magpen,magsus,mmom,&
   & mpatpol,mpdir,mpert,natom,nblok,ndim,nmat,nmdir,prtvol,zfield)
   end if

 end do

!Deallocations
 ABI_FREE(magsus)
 ABI_FREE(invbarmagsus)
 ABI_FREE(mmom)
 ABI_FREE(zfield)

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
!! ndim= dimension of the square susceptibilities 
!! prtvol= control the volume of information written on output
!!
!! OUTPUT
!! invbarmagsus(ndim,ndim)= Inverse of the penalized spin-sussceptibility tensor
!! magsus(ndim,ndim)= Spin-sussceptibility tensor
!!
!! SOURCE

 subroutine spinsus(blkval,iblok,invbarmagsus,magpen,magsus,&
& mpatpol,mpdir,mpert,natom,nblok,ndim,nmat,nmdir,prtvol)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,mpert,natom,nblok,ndim,nmat,nmdir,prtvol
 real(dp),intent(in) :: magpen
!arrays
 integer,intent(in) :: mpatpol(2),mpdir(3)
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
 complex(dp),intent(out) :: invbarmagsus(ndim,ndim)
 complex(dp),intent(out) :: magsus(ndim,ndim)

!Local variables -------------------------
!scalars
 integer :: iat1,iat2,icol,idir1,idir2,info,ipert1,ipert2,irow,lwork
 integer :: ipert1_red,ipert2_red,idir1_red,idir2_red
 character(len=1000) :: msg
!arrays
 complex(dp) :: barmagsus(ndim,ndim)
 complex(dp) :: idty(ndim,ndim)
 integer(dp) :: indexat(ndim),indexdir(ndim)
 integer, allocatable :: ipiv(:)
 complex(dp),allocatable :: work(:),work1(:,:),work2(:,:)
 character(len=1) :: cart(3)=(/'x','y','z'/)

! *********************************************************************

 !Extract the penalized susceptibility
 idty=(zero,zero)
 ipert2_red= 0
 do iat2= mpatpol(1), mpatpol(2)
   ipert2= natom + 11 + iat2
   ipert2_red= ipert2_red + 1
   idir2_red= 0
   do idir2= 1, 3
     if (mpdir(idir2)==0) cycle
     idir2_red= idir2_red + 1
     icol=idir2_red+(ipert2_red-1)*nmdir
     indexat(icol)=iat2
     indexdir(icol)=idir2
     idty(icol,icol)=(one,zero)
     ipert1_red=0
     do iat1= mpatpol(1), mpatpol(2)
       ipert1= natom + 11 + iat1
       ipert1_red= ipert1_red + 1
       idir1_red= 0
       do idir1= 1, 3
         if (mpdir(idir1)==0) cycle
         idir1_red=idir1_red+1
         irow=idir1_red+(ipert1_red-1)*nmdir

         !TODO: the four factor needs to be adapted in ABINIT: when
         !passing the magnetic moments to d2etot a 0.5 factor has to 
         !be applied to be consistent with the Bohr magneton in a.u. 
         !present in the Zeeman field perturbation
         barmagsus(irow,icol)= &
       & four*cmplx(blkval(1,idir1,ipert1,idir2,ipert2,iblok), &
       & blkval(2,idir1,ipert1,idir2,ipert2,iblok),16)

       end do
     end do
   end do
 end do

!Use magsus to store the intermediate array
 magsus=idty-magpen*barmagsus
     
!Invert the arrays
 ABI_MALLOC(work1,(ndim,ndim))
 ABI_MALLOC(work2,(ndim,ndim))
 work1=barmagsus
 work2=magsus

 ABI_MALLOC(ipiv,(ndim))
 call zgetrf( ndim, ndim, work1, ndim, ipiv, info )
 ABI_CHECK(info == 0, sjoin('zgetrf returned:', itoa(info)))

 ABI_MALLOC(work,(2))
 call zgetri( ndim, work1, ndim, ipiv, work, -1, info )
 ABI_CHECK(info == 0, sjoin('zgetri returned:', itoa(info)))
 lwork=int(work(1))

 ABI_REMALLOC(work,(lwork))
 call zgetri( ndim, work1, ndim, ipiv, work, lwork, info )
 ABI_CHECK(info == 0, sjoin('zgetri returned:', itoa(info)))


 call zgetrf( ndim, ndim, work2, ndim, ipiv, info )
 ABI_CHECK(info == 0, sjoin('zgetrf returned:', itoa(info)))

 call zgetri( ndim, work2, ndim, ipiv, work, -1, info )
 ABI_CHECK(info == 0, sjoin('zgetri returned:', itoa(info)))
 lwork=int(work(1))

 ABI_REMALLOC(work,(lwork))
 call zgetri( ndim, work2, ndim, ipiv, work, lwork, info )
 ABI_CHECK(info == 0, sjoin('zgetri returned:', itoa(info)))
 ABI_FREE(work)

 !Write the results in meaningfull arrays
 !use idty here to store an intermediate array
 invbarmagsus=work1
 
 !At last, calculate the susceptibility
 magsus=matmul(work2,barmagsus)

 ABI_FREE(ipiv)
 ABI_FREE(work1)
 ABI_FREE(work2)

 !Write results in output
 if (magpen > zero) then
   call wrtout([ab_out,std_out], ' Local spin susceptibility ')
   call wrtout([ab_out,std_out], '  atom1  dir  atom2  dir        Real              Imag')
   do irow=1, ndim
     do icol=1, ndim
       write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)' ) &
     & indexat(irow), cart(indexdir(irow)), indexat(icol), cart(indexdir(icol)), &
     & real(magsus(irow,icol)), aimag(magsus(irow,icol))
       call wrtout([ab_out,std_out], msg)
     end do
   end do
   call wrtout([ab_out,std_out], '   ')

   if (prtvol > 1) then
     call wrtout([ab_out,std_out], ' Penalized local spin susceptibility ')
     call wrtout([ab_out,std_out], '  atom1  dir  atom2  dir        Real              Imag')
     do irow=1, ndim
       do icol=1, ndim
         write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)' ) &
       & indexat(irow), cart(indexdir(irow)), indexat(icol), cart(indexdir(icol)), &
       & real(barmagsus(irow,icol)), aimag(barmagsus(irow,icol))
         call wrtout([ab_out,std_out], msg)
       end do
     end do
     call wrtout([ab_out,std_out], '   ')
  
     call wrtout([ab_out,std_out], ' Inverse of penalized local spin susceptibility ')
     call wrtout([ab_out,std_out], '  atom1  dir  atom2  dir        Real              Imag')
     do irow=1, ndim
       do icol=1, ndim
         write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)' ) &
       & indexat(irow), cart(indexdir(irow)), indexat(icol), cart(indexdir(icol)), &
       & real(invbarmagsus(irow,icol)), aimag(invbarmagsus(irow,icol))
         call wrtout([ab_out,std_out], msg)
       end do
     end do
     call wrtout([ab_out,std_out], '   ')
   end if

 end if
 
 end subroutine spinsus
!!***

!!****f* m_ddb_magpen/magmom
!! NAME
!! magmom
!!
!! FUNCTION
!! Calculate the first-order magnetic moments and the constrained
!! Zeeman fields
!!
!! INPUTS
!! blkval(2,3*mpert*3*mpert,nblok)=  Second-order derivative matrices
!!  In our case, the nblok is restricted to iblok
!! invbarmagsus(ndim,ndim)= Inverse of the penalized spin-sussceptibility tensor
!! iblok= index of the atomic displacement block
!! jblok= index of the electric field block
!! ip_from= lowest ipert to treat
!! ip_to= highest ipert to treat
!! magsus(ndim,ndim)= Spin-sussceptibility tensor
!! magpen = amplitude (in Ha) of the applied magnetic penalty 
!! mpatpol(2) = Atoms on which the magnetic penalty has been applied
!! mpdir(3) = Directions along which the spin-degrees of freedom have been stiffened 
!! mpert =maximum number of ipert
!! natom= number of atoms in unit cell
!! nblok= number of blocks in the DDB
!! nmat= number of atoms on which the magnetic penalty was applied
!!  (=1 if uniform Zeeman case)
!! nmdir= number of directions along which the magnetic penalty was applied
!! ndim= dimension of the square susceptibilities 
!! prtvol= control the volume of information written on output
!!
!! OUTPUT
!! mmom(ndim,(natom+2)*3)= first order magnetic moments on the atoms and 
!!  directions of the penalty induced by atomic displacements and/or electric fields.
!! zfield(ndim,(natom+2)*3)= Zeeman fields at constrained magnetic moments.
!!
!! SOURCE

 subroutine magmom(blkval,invbarmagsus,iblok,jblok,magpen,magsus,mmom,&
& mpatpol,mpdir,mpert,natom,nblok,ndim,nmat,nmdir,prtvol,zfield)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,jblok,mpert,natom,nblok,ndim,nmat,nmdir,prtvol
 real(dp),intent(in) :: magpen
!arrays
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
 integer,intent(in) :: mpatpol(2),mpdir(3)
 complex(dp),intent(in) :: invbarmagsus(ndim,ndim)
 complex(dp),intent(in) :: magsus(ndim,ndim)
 complex(dp),intent(inout) :: mmom(ndim,(natom+2)*3)
 complex(dp),intent(inout) :: zfield(ndim,(natom+2)*3)
!Local variables -------------------------
!scalars
 integer :: iat1,iat2,icol,idir1,idir2,ipert1,ipert2,irow
 integer :: ipert1_red,ipert2_red,idir1_red,idir2_red
 character(len=1000) :: msg
!arrays
 integer(dp) :: indexat1(ndim),indexdir1(ndim)
 integer(dp) :: indexat2((natom+2)*3),indexdir2((natom+2)*3)
 complex(dp) :: barmmom(ndim,(natom+2)*3)
 character(len=1) :: cart(3)=(/'x','y','z'/)

! *********************************************************************

!Extract the penalized moments
 do ipert2=1,natom+2
   do idir2=1,3
     icol=idir2+(ipert2-1)*3
     indexat2(icol)=ipert2
     indexdir2(icol)=idir2

     ipert1_red= 0
     do iat1= mpatpol(1), mpatpol(2)
       ipert1= natom + 11 + iat1
       ipert1_red= ipert1_red + 1
       idir1_red= 0
       do idir1= 1, 3
         if (mpdir(idir1)==0) cycle
         idir1_red= idir1_red + 1
         irow=idir1_red+(ipert1_red-1)*nmdir
         indexat1(irow)=iat1
         indexdir1(irow)=idir1
         
         if (iblok /=0 .and. ipert2 <= natom) then
           barmmom(irow,icol)=two* &
         & cmplx(blkval(1,idir1,ipert1,idir2,ipert2,iblok), &
         & blkval(2,idir1,ipert1,idir2,ipert2,iblok),16)
         else if (jblok /=0 .and. ipert2 == natom+2) then
           barmmom(irow,icol)=two* &
         & cmplx(blkval(1,idir1,ipert1,idir2,ipert2,jblok), &
         & blkval(2,idir1,ipert1,idir2,ipert2,jblok),16)
         end if

       end do
     end do
   end do
 end do

!Compute the Zeeman fields 
 zfield=matmul(invbarmagsus,barmmom)

!Compute the moments
 mmom=matmul(magsus,zfield)

!Write the results
 if (magpen > zero) then
   if (iblok /= 0) then
     call wrtout([ab_out,std_out], ' Local magnetic moments induced by atomic displacements ')
     call wrtout([ab_out,std_out], '  atom1  dir  atom2  dir        Real              Imag')
     do irow=1, ndim
       do icol=1, natom*3
         write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)' ) &
       & indexat1(irow), cart(indexdir1(irow)), indexat2(icol), cart(indexdir2(icol)), &
       & real(mmom(irow,icol)), aimag(mmom(irow,icol))
         call wrtout([ab_out,std_out], msg)
       end do
     end do
     call wrtout([ab_out,std_out], '   ')
   end if
   if (jblok /= 0) then
     call wrtout([ab_out,std_out], ' Local magnetic moments induced by electric field ')
     call wrtout([ab_out,std_out], '  atom1  dir  efld. dir        Real              Imag')
     do irow=1, ndim
       do icol=(natom+1)*3+1, (natom+2)*3
         write(msg,'(i4,4x,a2,4x,a2,2x,2es18.9)' ) &
       & indexat1(irow), cart(indexdir1(irow)), cart(indexdir2(icol)), &
       & real(mmom(irow,icol)), aimag(mmom(irow,icol))
         call wrtout([ab_out,std_out], msg)
       end do
     end do
     call wrtout([ab_out,std_out], '   ')
   end if
 end if
  
 end subroutine magmom

!!***
end module m_ddb_magpen
!!***
