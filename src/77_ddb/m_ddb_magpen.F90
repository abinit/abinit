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
!! Convert the second- and possibly third-order derivatives calculated with the magnetic penalty
!! into physically relevant quantities.
!!
!! INPUTS
!! ddb (INOUT) = ddb block datastructure
!! ddb_lw (INOUT) = ddb_lw block datastructure
!! magpen = amplitude (in Ha) of the applied magnetic penalty 
!! mpatpol(2) = Atoms on which the magnetic penalty has been applied
!! mpdir(3) = Directions along which the spin-degrees of freedom have been stiffened 
!! mpert = maximum number of ipert
!! mpopt = 1 calculate the frozen-magnetic second-order quantities
!!         2 calculate the spin-relaxed second-order quantities 
!! natom= number of atoms in unit cell
!! ntypat= number of atom types
!! timdisp= 1 calculate the third-order Berry curvatures
!! ucvol= unit cell volume
!!
!! OUTPUT
!! ddb= ddb%val updated with the corrected second-order derivatives
!!
!! SOURCE

 subroutine ddb_magpen(ddb,ddb_lw,magpen,mpatpol,mpdir,mpert,mpopt,natom, &
& ntypat,prtvol,rftyp,ucvol,timdisp,xred)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,mpopt,natom,ntypat,prtvol,rftyp,timdisp
 real(dp),intent(in) :: magpen,ucvol
!arrays
 type(ddb_type),intent(inout) :: ddb,ddb_lw
 integer,intent(in) :: mpatpol(2),mpdir(3)
 real(dp),intent(in) :: xred(3,natom)

!Local variables -------------------------
!scalars
 integer :: iblok,ii,ipert1,ipert2,jblok,kblok,lblok,nblok,nsize,ndim 
 integer :: nmat,nmdir
 character(len=500) :: msg
!arrays
 integer :: rfelfd(4),rfphon(4),rfstrs(4),rfmagn(4),rffreq(4)
 real(dp) :: qphnrm(3),qphon(3,3)
 complex(dpc), allocatable :: barmagsus(:,:),invbarmagsus(:,:)
 complex(dpc), allocatable :: invmagsus(:,:), magsus(:,:), invhmat(:,:)
 complex(dpc), allocatable :: barmmom(:,:),mmom(:,:), zfield(:,:)
 complex(dpc), allocatable :: bc_barmagsus(:,:),bc_ss(:,:),bc_sp(:,:)

! *********************************************************************
 write(msg, '(2a,(80a),4a)' ) ch10,('=',ii=1,80),ch10,ch10,&
 ' Linear-response magnetic penalty section ',ch10
 call wrtout([std_out, ab_out], msg)

 if (magpen<zero) then
   nmat= 1
 else if (magpen>zero) then
   nmat= mpatpol(2) - mpatpol(1) + 1
 end if
 nmdir=sum(mpdir(:))
 ndim=nmat*nmdir
 ABI_MALLOC(barmagsus,(ndim,ndim))
 ABI_MALLOC(magsus,(ndim,ndim))
 ABI_MALLOC(invbarmagsus,(ndim,ndim))
 ABI_MALLOC(invmagsus,(ndim,ndim))
 ABI_MALLOC(invhmat,(ndim,ndim))
 ABI_MALLOC(barmmom,(ndim,(natom+2)*3))
 ABI_MALLOC(mmom,(ndim,(natom+2)*3))
 ABI_MALLOC(zfield,(ndim,(natom+2)*3))

 nblok=ddb%nblok
 do kblok=1,nblok

   ! Look for the spin-susceptibility block in the DDB
   qphon=zero
   qphon(:,1)=ddb%qpt(1:3,kblok)
   qphnrm(:)=ddb%nrm(1,kblok)
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
       write(msg, '(2a,(80a),4a)' ) ch10,('-',ii=1,80),ch10,ch10,&
       ' Spin susceptibility (Uniform Zeeman) ',ch10
       call wrtout([std_out, ab_out], msg)
     else if (magpen>zero) then
       write(msg, '(2a,(80a),4a)' ) ch10,('-',ii=1,80),ch10,ch10,&
       ' Spin susceptibility (Local Zeeman) ',ch10
       call wrtout([std_out, ab_out], msg)
     end if

     call spinsus(barmagsus,ddb%val,iblok,invbarmagsus,invmagsus,invhmat,magpen,magsus,&
   & mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,prtvol)

   end if

   ! Calculate and write the induced magnetic moments 
   write(msg, '(2a,(80a),4a)' ) ch10,('-',ii=1,80),ch10,ch10,&
   ' First-order magnetic moments ',ch10
   call wrtout([std_out, ab_out], msg)

   ! First atomic-displacement
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
   jblok=0
   if (qphnrm(1)<tol8) then
     rfphon(2)=0
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
   end if

   if (iblok /= 0 .or. jblok /=0) then
     call magmom(barmmom,ddb%val,invbarmagsus,invhmat,iblok,jblok,magpen,magsus,mmom,&
   & mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,prtvol,qphon,xred,zfield)
   end if

   !Now calculate the non-magnetic second-order quantities
   write(msg, '(2a,(80a),4a)' ) ch10,('-',ii=1,80),ch10,ch10,&
   ' Second-order linear-response tensors ',ch10
   call wrtout([std_out, ab_out], msg)

   rfmagn(:)=0
   rfelfd(:)=0
   rfphon(:)=0
   
   !IFCs block
   rfphon(1:2)=1
   call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp)
   if (iblok /= 0 ) then
     call mp_ifc(barmagsus,ddb%val,iblok,magsus,mpert,mpopt,&
   & natom,nblok,ndim,prtvol,qphon,xred,zfield)
   end if

   !Born effective charges block
   jblok=0
   if (qphnrm(1)<tol8) then
     rfphon(1:2)=1
     rfelfd(1:2)=2
     call ddb%get_block(jblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp)
   end if
   if (jblok /= 0 ) then
     call mp_zeff(barmagsus,ddb%val,jblok,magsus,mpert,mpopt,&
   & natom,nblok,ndim,prtvol,ucvol,zfield)
   end if

   !Dielectric susceptibility block
   lblok=0
   if (qphnrm(1)<tol8) then
     rfphon(:)=0
     rfelfd(1:2)=2
     call ddb%get_block(lblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp)
   end if
   if (lblok /= 0 ) then
     call mp_diel(barmagsus,ddb%val,lblok,magsus,mpert,mpopt,&
   & natom,nblok,ndim,prtvol,ucvol,zfield)
   end if

 end do

 if (timdisp==1) then

   ABI_MALLOC(bc_barmagsus,(ndim,ndim))
   ABI_MALLOC(bc_ss,(ndim,ndim))
   ABI_MALLOC(bc_sp,(ndim,(natom+2)*3))

   write(msg, '(2a,(80a),4a)' ) ch10,('=',ii=1,80),ch10,ch10,&
   ' Frequency-derivatives magnetic penalty section ',ch10
   call wrtout([std_out, ab_out], msg)

   rffreq(:)=0
   do kblok=1,nblok

     !Berry curvature of the penalized spin-susceptibility
     qphon=zero
     qphon(:,1)=ddb_lw%qpt(1:3,kblok)
     qphnrm(:)=ddb_lw%nrm(1,kblok)
     rfphon(1:3)=0
     rfelfd(1:3)=0
     rfstrs(1:3)=0
     rffreq(3)=1
     iblok=0
     if (magpen<zero) then
       rfmagn(1:2)= 1
     else if (magpen>zero) then
       rfmagn(1:2)= 2
     end if
     call ddb_lw%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, 33, &
   & mpatpol=mpatpol,mpdir=mpdir,rfmagn=rfmagn,rffreq=rffreq)

     if (iblok /= 0) then
       call berrycurv_ss(bc_barmagsus,bc_ss,ddb_lw%val,iblok,invbarmagsus,mpatpol,mpdir,mpert,&
     & natom,nblok,ndim,nmdir,prtvol)
     end if

     !Berry curvature of the induced Zeeman fields

     !First atomic-displacement
     rfphon(2)=1
     rfmagn(:)=0
     if (magpen<zero) then
       rfmagn(1)= 1
     else if (magpen>zero) then
       rfmagn(1)= 2
     end if

     call ddb_lw%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, 33, &
   & mpatpol=mpatpol,mpdir=mpdir,rfmagn=rfmagn,rffreq=rffreq)

     ! Then electric field
     jblok=0
     if (qphnrm(1)<tol8) then
       rfphon(2)=0
       rfelfd(2)=2

       call ddb_lw%get_block(jblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, 33, &
     & mpatpol=mpatpol,mpdir=mpdir,rfmagn=rfmagn,rffreq=rffreq)
     end if

     if (iblok /= 0 .or. jblok /=0) then
       call berrycurv_sp(barmmom,bc_sp,bc_ss,ddb_lw%val,iblok,invbarmagsus,jblok, &
     & mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,prtvol,qphon,xred)
     end if

     !Berry curvature of other second-order quantites

     !IFCs
     rfphon(1:2)=1
     rfelfd(:)=0
     rfmagn(:)=0
     call ddb_lw%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, 33, &
   & rffreq=rffreq)
     if (iblok /= 0 ) then
       call berrycurv_pp(barmagsus,bc_barmagsus,bc_sp,ddb_lw%val,iblok, &
     & mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,prtvol,qphon,xred,zfield)
     end if
   end do 

   ABI_FREE(bc_ss)
 end if

!Deallocations
 ABI_FREE(barmagsus)
 ABI_FREE(magsus)
 ABI_FREE(invbarmagsus)
 ABI_FREE(invmagsus)
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
!! nmdir= number of directions along which the magnetic penalty was applied
!! ndim= dimension of the square susceptibilities 
!! prtvol= control the volume of information written on output
!!
!! OUTPUT
!! barmagsus(ndim,ndim)= Penalized spin-sussceptibility tensor
!! invbarmagsus(ndim,ndim)= Inverse of the penalized spin-sussceptibility tensor
!! magsus(ndim,ndim)= Spin-sussceptibility tensor
!! invmagsus(ndim,ndim)= Inverse of spin-sussceptibility tensor
!!
!! SOURCE

 subroutine spinsus(barmagsus,blkval,iblok,invbarmagsus,invmagsus,invhmat,magpen,magsus,&
& mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,prtvol)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,mpert,natom,nblok,ndim,nmdir,prtvol
 real(dp),intent(in) :: magpen
!arrays
 integer,intent(in) :: mpatpol(2),mpdir(3)
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
 complex(dpc),intent(out) :: barmagsus(ndim,ndim)
 complex(dpc),intent(out) :: invbarmagsus(ndim,ndim)
 complex(dpc),intent(out) :: invhmat(ndim,ndim)
 complex(dpc),intent(out) :: magsus(ndim,ndim)
 complex(dpc),intent(out) :: invmagsus(ndim,ndim)

!Local variables -------------------------
!scalars
 integer :: iat1,iat2,icol,idir1,idir2,info,ipert1,ipert2,irow,lwork
 integer :: ipert1_red,ipert2_red,idir1_red,idir2_red
 real(dp) :: fac
 character(len=1000) :: msg
!arrays
 complex(dpc) :: idty(ndim,ndim)
 integer(dp) :: indexat(ndim),indexdir(ndim)
 integer, allocatable :: ipiv(:)
 complex(dpc),allocatable :: work(:),work1(:,:),work2(:,:)
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

         barmagsus(irow,icol)= &
       & cmplx(blkval(1,idir1,ipert1,idir2,ipert2,iblok), &
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
 
 !At last, calculate the susceptibility and its inverse
 invhmat=work2
 magsus=matmul(work2,barmagsus)
 invmagsus=invbarmagsus-magpen*idty

 ABI_FREE(ipiv)
 ABI_FREE(work1)
 ABI_FREE(work2)

 fac=2.714943600699**2*27.2114/four
 open(10,file='k_ss.txt')
   do irow=1, ndim
     write(10,*) invmagsus(irow,1:ndim)*fac
   end do 
 close(10)

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


   call wrtout([ab_out,std_out], ' Inverse of local spin susceptibility ')
   call wrtout([ab_out,std_out], '  atom1  dir  atom2  dir        Real              Imag')
   do irow=1, ndim
     do icol=1, ndim
       write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)' ) &
     & indexat(irow), cart(indexdir(irow)), indexat(icol), cart(indexdir(icol)), &
     & real(invmagsus(irow,icol)), aimag(invmagsus(irow,icol))
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

     call wrtout([ab_out,std_out], ' Inverse of H matrix (I-\alpha \barchi)^{-1} ')
     call wrtout([ab_out,std_out], '  atom1  dir  atom2  dir        Real              Imag')
     do irow=1, ndim
       do icol=1, ndim
         write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)' ) &
       & indexat(irow), cart(indexdir(irow)), indexat(icol), cart(indexdir(icol)), &
       & real(invhmat(irow,icol)), aimag(invhmat(irow,icol))
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
!! magsus(ndim,ndim)= Spin-sussceptibility tensor
!! magpen = amplitude (in Ha) of the applied magnetic penalty 
!! mpatpol(2) = Atoms on which the magnetic penalty has been applied
!! mpdir(3) = Directions along which the spin-degrees of freedom have been stiffened 
!! mpert =maximum number of ipert
!! natom= number of atoms in unit cell
!! nblok= number of blocks in the DDB
!! nmdir= number of directions along which the magnetic penalty was applied
!! ndim= dimension of the square susceptibilities 
!! prtvol= control the volume of information written on output
!!
!! OUTPUT
!! barmmom(ndim,(natom+2)*3)= penalized first order magnetic moments on the atoms and 
!!  directions of the penalty induced by atomic displacements and/or electric fields.
!! mmom(ndim,(natom+2)*3)= first order magnetic moments on the atoms and 
!!  directions of the penalty induced by atomic displacements and/or electric fields.
!! zfield(ndim,(natom+2)*3)= Zeeman fields at constrained magnetic moments.
!!
!! SOURCE

 subroutine magmom(barmmom,blkval,invbarmagsus,invhmat,iblok,jblok,magpen,magsus,mmom,&
& mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,prtvol,qphon,xred,zfield)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,jblok,mpert,natom,nblok,ndim,nmdir,prtvol
 real(dp),intent(in) :: magpen
!arrays
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
 real(dp),intent(in) :: qphon(3),xred(3,natom)
 integer,intent(in) :: mpatpol(2),mpdir(3)
 complex(dpc),intent(out) :: barmmom(ndim,(natom+2)*3)
 complex(dpc),intent(in) :: invbarmagsus(ndim,ndim)
 complex(dpc),intent(in) :: invhmat(ndim,ndim)
 complex(dpc),intent(in) :: magsus(ndim,ndim)
 complex(dpc),intent(out) :: mmom(ndim,(natom+2)*3)
 complex(dpc),intent(out) :: zfield(ndim,(natom+2)*3)
!Local variables -------------------------
!scalars
 integer :: iat1,iat2,icol,idir1,idir2,ipert1,ipert2,irow
 integer :: ipert1_red,ipert2_red,idir1_red,idir2_red
 real(dp) :: fac
 character(len=1000) :: msg
!arrays
 integer(dp) :: indexat1(ndim),indexdir1(ndim)
 integer(dp) :: indexat2((natom+2)*3),indexdir2((natom+2)*3)
 complex(dpc) :: mmom_alt(ndim,(natom+2)*3)
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
           barmmom(irow,icol)= &
         & cmplx(blkval(1,idir1,ipert1,idir2,ipert2,iblok), &
         & blkval(2,idir1,ipert1,idir2,ipert2,iblok),16)
         else if (jblok /=0 .and. ipert2 == natom+2) then
           barmmom(irow,icol)= &
         & cmplx(blkval(1,idir1,ipert1,idir2,ipert2,jblok), &
         & blkval(2,idir1,ipert1,idir2,ipert2,jblok),16)
         end if

       end do
     end do
   end do
 end do

!Compute the Zeeman fields 
 zfield=-matmul(invbarmagsus,barmmom)

!Compute the moments
 mmom=-matmul(magsus,zfield)
 mmom_alt=matmul(invhmat,barmmom)

 fac=2.714943600699/two*27.2114/0.529177

 open(10,file='k_ps.txt')
 do iat1= 1, natom
   do idir1= 1, 3
     icol= (iat1-1)*3 + idir1
     write(10,*) conjg(zfield(1:ndim,icol)*fac* &
   & exp(two_pi*(0.d0,1.d0)* dot_product(qphon,xred(:,iat1))))
   end do
 end do 
 close(10)

!Write the results
 if (magpen > zero) then
   if (iblok /= 0) then
     call wrtout([ab_out,std_out], ' Local Zeeman fields induced by atomic displacements (at constrained magnetic moments)')
     call wrtout([ab_out,std_out], '  atom1  dir  atom2  dir        Real              Imag')
     do irow=1, ndim
       do icol=1, natom*3
         write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)' ) &
       & indexat1(irow), cart(indexdir1(irow)), indexat2(icol), cart(indexdir2(icol)), &
       & real(zfield(irow,icol)), aimag(zfield(irow,icol))
         call wrtout([ab_out,std_out], msg)
       end do
     end do
     call wrtout([ab_out,std_out], '   ')

     call wrtout([ab_out,std_out], ' Local magnetic moments induced by atomic displacements (from induced Zeeman fields)')
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

     call wrtout([ab_out,std_out], ' Local magnetic moments induced by atomic displacements (from induced penalized moments)')
     call wrtout([ab_out,std_out], '  atom1  dir  atom2  dir        Real              Imag')
     do irow=1, ndim
       do icol=1, natom*3
         write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)' ) &
       & indexat1(irow), cart(indexdir1(irow)), indexat2(icol), cart(indexdir2(icol)), &
       & real(mmom_alt(irow,icol)), aimag(mmom_alt(irow,icol))
         call wrtout([ab_out,std_out], msg)
       end do
     end do
     call wrtout([ab_out,std_out], '   ')
   end if
   if (jblok /= 0) then
     call wrtout([ab_out,std_out], ' Local Zeeman fields induced by electric field (at constrained magnetic moments)')
     call wrtout([ab_out,std_out], '  atom1  dir  efld.dir         Real              Imag')
     do irow=1, ndim
       do icol=(natom+1)*3+1, (natom+2)*3
         write(msg,'(i4,4x,a2,4x,a2,6x,2es18.9)' ) &
       & indexat1(irow), cart(indexdir1(irow)), cart(indexdir2(icol)), &
       & real(zfield(irow,icol)), aimag(zfield(irow,icol))
         call wrtout([ab_out,std_out], msg)
       end do
     end do
     call wrtout([ab_out,std_out], '   ')
     call wrtout([ab_out,std_out], ' Local magnetic moments induced by electric field (from induced Zeeman fields)')
     call wrtout([ab_out,std_out], '  atom1  dir  efld.dir         Real              Imag')
     do irow=1, ndim
       do icol=(natom+1)*3+1, (natom+2)*3
         write(msg,'(i4,4x,a2,4x,a2,6x,2es18.9)' ) &
       & indexat1(irow), cart(indexdir1(irow)), cart(indexdir2(icol)), &
       & real(mmom(irow,icol)), aimag(mmom(irow,icol))
         call wrtout([ab_out,std_out], msg)
       end do
     end do
     call wrtout([ab_out,std_out], '   ')
     call wrtout([ab_out,std_out], ' Local magnetic moments induced by electric field (from induced penalized moments)')
     call wrtout([ab_out,std_out], '  atom1  dir  efld.dir         Real              Imag')
     do irow=1, ndim
       do icol=(natom+1)*3+1, (natom+2)*3
         write(msg,'(i4,4x,a2,4x,a2,6x,2es18.9)' ) &
       & indexat1(irow), cart(indexdir1(irow)), cart(indexdir2(icol)), &
       & real(mmom_alt(irow,icol)), aimag(mmom_alt(irow,icol))
         call wrtout([ab_out,std_out], msg)
       end do
     end do
     call wrtout([ab_out,std_out], '   ')
   end if
   if (prtvol > 1) then
     if (iblok /= 0) then
       call wrtout([ab_out,std_out], ' Penalized local magnetic moments induced by atomic displacements ')
       call wrtout([ab_out,std_out], '  atom1  dir  atom2  dir        Real              Imag')
       do irow=1, ndim
         do icol=1, natom*3
           write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)' ) &
         & indexat1(irow), cart(indexdir1(irow)), indexat2(icol), cart(indexdir2(icol)), &
         & real(barmmom(irow,icol)), aimag(barmmom(irow,icol))
           call wrtout([ab_out,std_out], msg)
         end do
       end do
       call wrtout([ab_out,std_out], '   ')
     end if
     if (jblok /= 0) then
       call wrtout([ab_out,std_out], ' Penalized local magnetic moments induced by electric field ')
       call wrtout([ab_out,std_out], '  atom1  dir  efld.dir         Real              Imag')
       do irow=1, ndim
         do icol=(natom+1)*3+1, (natom+2)*3
           write(msg,'(i4,4x,a2,4x,a2,6x,2es18.9)' ) &
         & indexat1(irow), cart(indexdir1(irow)), cart(indexdir2(icol)), &
         & real(barmmom(irow,icol)), aimag(barmmom(irow,icol))
           call wrtout([ab_out,std_out], msg)
         end do
       end do
       call wrtout([ab_out,std_out], '   ')
     end if
   end if
 end if
  
 end subroutine magmom
!!***

!!****f* m_ddb_magpen/mp_ifc
!! NAME
!! mp_ifc
!!
!! FUNCTION
!! Calculate the different flavor (see mpopt below) of the second-
!! order linear-response quantities
!!
!! INPUTS
!! barmagsus(ndim,ndim)= Penalized spin-sussceptibility tensor
!! blkval(2,3*mpert*3*mpert,nblok)=  Second-order derivative matrices
!!  In our case, the nblok is restricted to iblok
!! iblok= index of the IFCs block
!! magsus(ndim,ndim)= Spin-sussceptibility tensor
!! mpert =maximum number of ipert
!! mpopt = 1 calculate the frozen-magnetic second-order quantities
!!         2 calculate the spin-relaxed second-order quantities 
!! natom= number of atoms in unit cell
!! nblok= number of blocks in the DDB
!! ndim= dimension of the square susceptibilities 
!! prtvol= control the volume of information written on output
!!
!! OUTPUT
!!  directions of the penalty induced by atomic displacements and/or electric fields.
!! zfield(ndim,(natom+2)*3)= Zeeman fields at constrained magnetic moments.
!!
!! SOURCE

 subroutine mp_ifc(barmagsus,blkval,iblok,magsus,mpert,mpopt,&
& natom,nblok,ndim,prtvol,qphon,xred,zfield)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,mpert,mpopt,natom,nblok,ndim,prtvol
!arrays
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
 real(dp),intent(in) :: qphon(3),xred(3,natom)
 complex(dpc),intent(in) :: barmagsus(ndim,ndim)
 complex(dpc),intent(in) :: magsus(ndim,ndim)
 complex(dpc),intent(in) :: zfield(ndim,(natom+2)*3)
!Local variables -------------------------
!scalars
 integer :: icol,idir1,idir2,ipert1,ipert2,irow
 real(dp) :: fac
 character(len=1000) :: msg
!arrays
 character(len=1) :: cart(3)=(/'x','y','z'/)
 complex(dpc) :: barifc(natom*3,natom*3)
 complex(dpc) :: fmifc(natom*3,natom*3)
 complex(dpc) :: fmifc_sf(natom*3,natom*3)
 complex(dpc) :: srifc(natom*3,natom*3)
 complex(dpc) :: ifc_zfield(ndim,natom*3)

! *********************************************************************

!Extract the penalized IFCs 
 do ipert2= 1, natom
   do idir2= 1, 3
     icol=( ipert2-1)*3 + idir2
     do ipert1= 1, natom
       do idir1= 1, 3
         irow=( ipert1-1)*3 + idir1
         barifc(irow,icol)= &
       & cmplx(blkval(1,idir1,ipert1,idir2,ipert2,iblok), &
       &       blkval(2,idir1,ipert1,idir2,ipert2,iblok),16)
       end do
     end do
   end do
 end do 

 !Calculate the frozen-magnetic flavor
 ifc_zfield(:,:)=zfield(:,1:natom*3)
 fmifc= matmul(transpose(conjg(ifc_zfield)),matmul(barmagsus,ifc_zfield))
 fmifc= barifc + fmifc

 !Calculate the spin-relaxed flavor
 if (mpopt==2) then
   srifc= matmul(transpose(conjg(ifc_zfield)),matmul(magsus,ifc_zfield))
   srifc= fmifc - srifc
 end if

 do ipert2= 1, natom
   do idir2= 1, 3
     icol=( ipert2-1)*3 + idir2
     do ipert1= 1, natom
       do idir1= 1, 3
         irow=( ipert1-1)*3 + idir1
         fmifc_sf(irow,icol)=fmifc(irow,icol)* &
       & exp(two_pi*(0.d0,1.d0)* dot_product(qphon,xred(:,ipert2)-xred(:,ipert1)))
       end do
     end do
   end do
 end do 
 
 fac=27.2114/(0.52917)**2
 open(10,file='k_pp.txt')
 do irow=1,natom*3
    write(10,*) fmifc_sf(irow,1:natom*3)*fac
 end do 
 close(10)
 
 !Write the results
 call wrtout([ab_out,std_out], ' Frozen-magnetic interatomic force constants')
 call wrtout([ab_out,std_out], '  atom1  dir  atom2  dir        Real              Imag')
 do ipert1= 1, natom
   do idir1= 1, 3
     irow=( ipert1-1)*3 + idir1
     do ipert2= 1, natom
       do idir2= 1, 3
         icol=( ipert2-1)*3 + idir2
         write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)') &
       & ipert1, cart(idir1), ipert2, cart(idir2), &
       & real(fmifc(irow,icol)), aimag(fmifc(irow,icol))
         call wrtout([ab_out,std_out], msg)
       end do
     end do
     call wrtout([ab_out,std_out], '   ')
   end do
 end do 

 if (mpopt==2) then
   call wrtout([ab_out,std_out], ' Spin-relaxed interatomic force constants')
   call wrtout([ab_out,std_out], '  atom1  dir  atom2  dir        Real              Imag')
   do ipert1= 1, natom
     do idir1= 1, 3
       irow=( ipert1-1)*3 + idir1
       do ipert2= 1, natom
         do idir2= 1, 3
           icol=( ipert2-1)*3 + idir2
           write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)') &
         & ipert1, cart(idir1), ipert2, cart(idir2), &
         & real(srifc(irow,icol)), aimag(srifc(irow,icol))
           call wrtout([ab_out,std_out], msg)
         end do
       end do
       call wrtout([ab_out,std_out], '   ')
     end do
   end do 
 end if

 end subroutine mp_ifc
!!***

!!****f* m_ddb_magpen/mp_diel
!! NAME
!! mp_diel
!!
!! FUNCTION
!! Calculate the different flavor (see mpopt below) of the second-
!! order linear-response quantities
!!
!! INPUTS
!! barmagsus(ndim,ndim)= Penalized spin-sussceptibility tensor
!! blkval(2,3*mpert*3*mpert,nblok)=  Second-order derivative matrices
!!  In our case, the nblok is restricted to iblok
!! iblok= index of the IFCs block
!! magsus(ndim,ndim)= Spin-sussceptibility tensor
!! mpert =maximum number of ipert
!! mpopt = 1 calculate the frozen-magnetic second-order quantities
!!         2 calculate the spin-relaxed second-order quantities 
!! natom= number of atoms in unit cell
!! nblok= number of blocks in the DDB
!! ndim= dimension of the square susceptibilities 
!! prtvol= control the volume of information written on output
!! ucvol= unit cell volume
!!
!! OUTPUT
!!  directions of the penalty induced by atomic displacements and/or electric fields.
!! zfield(ndim,(natom+2)*3)= Zeeman fields at constrained magnetic moments.
!!
!! SOURCE

 subroutine mp_diel(barmagsus,blkval,iblok,magsus,mpert,mpopt,&
& natom,nblok,ndim,prtvol,ucvol,zfield)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,mpert,mpopt,natom,nblok,ndim,prtvol
 real(dp), intent(in) :: ucvol
!arrays
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
 complex(dpc),intent(in) :: barmagsus(ndim,ndim)
 complex(dpc),intent(in) :: magsus(ndim,ndim)
 complex(dpc),intent(in) :: zfield(ndim,(natom+2)*3)

!Local variables -------------------------
!scalars
 integer :: idir1,idir2,ipert1,ipert2
 character(len=1000) :: msg
!arrays
 complex(dpc) :: barepsilon(3,3), bardielsus(3,3)
 complex(dpc) :: fmepsilon(3,3), fmdielsus(3,3)
 complex(dpc) :: srepsilon(3,3), srdielsus(3,3)
 complex(dpc) :: diel_zfield(ndim,3)
 character(len=1) :: cart(3)=(/'x','y','z'/)

! *********************************************************************

!Extract the penalized dielectric tensor
 ipert2= natom + 2
 ipert1= natom + 2
 do idir2= 1, 3
   do idir1= 1, 3
     barepsilon(idir1,idir2)= &
   & cmplx(blkval(1,idir1,ipert1,idir2,ipert2,iblok), &
   &       blkval(2,idir1,ipert1,idir2,ipert2,iblok),16)
     if (idir1==idir2) then
       bardielsus(idir1,idir2)= one - barepsilon(idir1,idir2)
     else
       bardielsus(idir1,idir2)= -barepsilon(idir1,idir2)
     end if
   end do
 end do
 bardielsus= ucvol/four_pi * bardielsus

!Calculate the frozen-magnetic flavor
 diel_zfield(:,:)=zfield(:,(natom+2)*3-2:(natom+2)*3)
 fmdielsus= matmul(transpose(conjg(diel_zfield)),matmul(barmagsus,diel_zfield))
 fmdielsus= bardielsus + fmdielsus
 fmepsilon= -four_pi/ucvol*fmdielsus
 fmepsilon(1,1)= one + fmepsilon(1,1)
 fmepsilon(2,2)= one + fmepsilon(2,2)
 fmepsilon(3,3)= one + fmepsilon(3,3)

!Calculate the spin-relaxed flavor
 if (mpopt==2) then
   srdielsus= matmul(transpose(conjg(diel_zfield)),matmul(magsus,diel_zfield))
   srdielsus= fmdielsus - srdielsus
   srepsilon= -four_pi/ucvol*srdielsus
   srepsilon(1,1)= one + srepsilon(1,1)
   srepsilon(2,2)= one + srepsilon(2,2)
   srepsilon(3,3)= one + srepsilon(3,3)
 end if

 !Write the results
 call wrtout([ab_out,std_out], ' Frozen-magnetic dielectric tensor (clamped ion)')
 call wrtout([ab_out,std_out], '  dir  dir        Real              Imag')
 do idir1= 1, 3
   do idir2= 1, 3
     write(msg,'(2x,a2,3x,a2,2x,2es18.9)') cart(idir1), cart(idir2), &
   & real(fmepsilon(idir1,idir2)), aimag(fmepsilon(idir1,idir2))
     call wrtout([ab_out,std_out], msg)
   end do
 end do
 call wrtout([ab_out,std_out], '   ')

 if (mpopt==2) then
   call wrtout([ab_out,std_out], ' Spin-relaxed dielectric tensor (clamped ion)')
   call wrtout([ab_out,std_out], '  dir  dir        Real              Imag')
   do idir1= 1, 3
     do idir2= 1, 3
       write(msg,'(2x,a2,3x,a2,2x,2es18.9)' ) cart(idir1), cart(idir2), &
     & real(srepsilon(idir1,idir2)), aimag(srepsilon(idir1,idir2))
       call wrtout([ab_out,std_out], msg)
     end do
   end do
   call wrtout([ab_out,std_out], '   ')
 end if
 
 if (prtvol > 1) then
   call wrtout([ab_out,std_out], ' Penalized dielectric tensor (clamped ion)')
   call wrtout([ab_out,std_out], '  dir  dir        Real              Imag')
   do idir1= 1, 3
     do idir2= 1, 3
       write(msg,'(2x,a2,3x,a2,2x,2es18.9)') cart(idir1), cart(idir2), &
     & real(barepsilon(idir1,idir2)), aimag(barepsilon(idir1,idir2))
       call wrtout([ab_out,std_out], msg)
     end do
   end do
   call wrtout([ab_out,std_out], '   ')
 end if

 end subroutine mp_diel
!!***

!!****f* m_ddb_magpen/mp_zeff
!! NAME
!! mp_zeff
!!
!! FUNCTION
!! Calculate the different flavor (see mpopt below) of the second-
!! order linear-response quantities
!!
!! INPUTS
!! barmagsus(ndim,ndim)= Penalized spin-sussceptibility tensor
!! blkval(2,3*mpert*3*mpert,nblok)=  Second-order derivative matrices
!!  In our case, the nblok is restricted to iblok
!! iblok= index of the IFCs block
!! magsus(ndim,ndim)= Spin-sussceptibility tensor
!! mpert =maximum number of ipert
!! mpopt = 1 calculate the frozen-magnetic second-order quantities
!!         2 calculate the spin-relaxed second-order quantities 
!! natom= number of atoms in unit cell
!! nblok= number of blocks in the DDB
!! ndim= dimension of the square susceptibilities 
!! prtvol= control the volume of information written on output
!! ucvol= unit cell volume
!!
!! OUTPUT
!!  directions of the penalty induced by atomic displacements and/or electric fields.
!! zfield(ndim,(natom+2)*3)= Zeeman fields at constrained magnetic moments.
!!
!! SOURCE

 subroutine mp_zeff(barmagsus,blkval,iblok,magsus,mpert,mpopt,&
& natom,nblok,ndim,prtvol,ucvol,zfield)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,mpert,mpopt,natom,nblok,ndim,prtvol
 real(dp), intent(in) :: ucvol
!arrays
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
 complex(dpc),intent(in) :: barmagsus(ndim,ndim)
 complex(dpc),intent(in) :: magsus(ndim,ndim)
 complex(dpc),intent(in) :: zfield(ndim,(natom+2)*3)

!Local variables -------------------------
!scalars
 integer :: icol,idir1,idir2,ipert1,ipert2
 character(len=1000) :: msg
!arrays
 complex(dpc) :: barzeff(3,natom*3)
 complex(dpc) :: fmzeff(3,natom*3)
 complex(dpc) :: srzeff(3,natom*3)
 complex(dpc) :: diel_zfield(ndim,3)
 complex(dpc) :: ifc_zfield(ndim,natom*3)
 character(len=1) :: cart(3)=(/'x','y','z'/)

! *********************************************************************

!Extract the penalized Born charges
 ipert1= natom + 2
 do ipert2= 1, natom
   do idir2= 1, 3
     icol=( ipert2-1)*3 + idir2
     do idir1= 1, 3
       barzeff(idir1,icol)=half * &
     & (cmplx(blkval(1,idir1,ipert1,idir2,ipert2,iblok), &
     &        blkval(2,idir1,ipert1,idir2,ipert2,iblok),16) + &
     &  conjg(cmplx(blkval(1,idir2,ipert2,idir1,ipert1,iblok), &
     &        blkval(2,idir2,ipert2,idir1,ipert1,iblok),16)))
     end do
   end do
 end do

 !Calculate the frozen-magnetic flavor
 ifc_zfield(:,:)=zfield(:,1:natom*3)
 diel_zfield(:,:)=zfield(:,(natom+2)*3-2:(natom+2)*3)
 fmzeff= matmul(transpose(conjg(diel_zfield)),matmul(barmagsus,ifc_zfield))
 fmzeff= barzeff + fmzeff

 !Calculate the spin-relaxed flavor
 if (mpopt==2) then
   srzeff= matmul(transpose(conjg(diel_zfield)),matmul(magsus,ifc_zfield))
   srzeff= fmzeff - srzeff
 end if

 !Write the results
 call wrtout([ab_out,std_out], ' Frozen-magnetic Born effective charges')
 call wrtout([ab_out,std_out], ' efld.dir   atom   dir        Real              Imag')
 do idir1= 1, 3
   do ipert2= 1, natom
     do idir2= 1, 3
       icol=( ipert2-1)*3 + idir2
       write(msg,'(3x,a2,7x,i3,4x,a2,2x,2es18.9)') &
     & cart(idir1), ipert2, cart(idir2), &
     & real(fmzeff(idir1,icol)), aimag(fmzeff(idir1,icol))
       call wrtout([ab_out,std_out], msg)
     end do
   end do
   call wrtout([ab_out,std_out], '   ')
 end do 
 if (mpopt==2) then
   call wrtout([ab_out,std_out], ' Spin-relaxed Born effective charges')
   call wrtout([ab_out,std_out], ' efld.dir   atom   dir        Real              Imag')
   do idir1= 1, 3
     do ipert2= 1, natom
       do idir2= 1, 3
         icol=( ipert2-1)*3 + idir2
         write(msg,'(3x,a2,7x,i3,4x,a2,2x,2es18.9)') &
       & cart(idir1), ipert2, cart(idir2), &
       & real(srzeff(idir1,icol)), aimag(srzeff(idir1,icol))
         call wrtout([ab_out,std_out], msg)
       end do
     end do
     call wrtout([ab_out,std_out], '   ')
   end do 
 end if

 if (prtvol > 1) then
   call wrtout([ab_out,std_out], ' Penalized Born effective charges')
   call wrtout([ab_out,std_out], ' efld.dir   atom   dir        Real              Imag')
   do idir1= 1, 3
     do ipert2= 1, natom
       do idir2= 1, 3
         icol=( ipert2-1)*3 + idir2
         write(msg,'(3x,a2,7x,i3,4x,a2,2x,2es18.9)') &
       & cart(idir1), ipert2, cart(idir2), &
       & real(barzeff(idir1,icol)), aimag(barzeff(idir1,icol))
         call wrtout([ab_out,std_out], msg)
       end do
     end do
     call wrtout([ab_out,std_out], '   ')
   end do 
 end if

 end subroutine mp_zeff
!!***

!!****f* m_ddb_magpen/berrycurv_ss
!! NAME
!! berrycurv_ss
!!
!! FUNCTION
!! Calculate the Berry curvature of the inverse magnetic susceptibility
!! (this is equivalent to the G^(ss) matrix of S.Ren et al.)
!!
!! INPUTS
!! blkval(2,3*mpert*3*mpert*3*mpert,nblok)=  Third-order derivative matrices
!!  In our case, the nblok is restricted to iblok
!! iblok= index of the current block
!! invbarmagsus(ndim,ndim)= Inverse of the penalized spin-sussceptibility tensor
!! mpatpol(2) = Atoms on which the magnetic penalty has been applied
!! mpdir(3) = Directions along which the spin-degrees of freedom have been stiffened 
!! mpert =maximum number of ipert
!! natom= number of atoms in unit cell
!! nblok= number of blocks in the DDB
!! nmdir= number of directions along which the magnetic penalty was applied
!! ndim= dimension of the square susceptibilities 
!! prtvol= control the volume of information written on output
!!
!! OUTPUT
!! bc_ss(ndim,ndim)= Berry-curvature of the inverse local-spin susceptibility
!!
!! SOURCE

 subroutine berrycurv_ss(bc_barmagsus,bc_ss,blkval,iblok,invbarmagsus,&
& mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,prtvol)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,mpert,natom,nblok,ndim,nmdir,prtvol
!arrays
 integer,intent(in) :: mpatpol(2),mpdir(3)
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,3,mpert,nblok)
 complex(dpc),intent(in) :: invbarmagsus(ndim,ndim)
 complex(dpc),intent(out) :: bc_ss(ndim,ndim)
 complex(dpc),intent(out) :: bc_barmagsus(ndim,ndim)
!Local variables -------------------------
!scalars
 integer :: iat1,iat2,icol,idir1,idir2,idir3,info,ipert1,ipert2,ipert3,irow,lwork
 integer :: ipert1_red,ipert2_red,idir1_red,idir2_red
 real(dp) :: fac
 complex(dpc), parameter :: ione=(0.d0,1.d0)
 character(len=1000) :: msg
!arrays
 complex(dpc) :: idty(ndim,ndim)
 integer(dp) :: indexat(ndim),indexdir(ndim)
 integer, allocatable :: ipiv(:)
 complex(dpc),allocatable :: work(:),work1(:,:),work2(:,:)
 character(len=1) :: cart(3)=(/'x','y','z'/)

! *********************************************************************

 !Extract the Berry-curvature of the penalized susceptibility
 idty=(zero,zero)
 ipert2_red= 0
 ipert3= natom + 9
 idir3= 1
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

         bc_barmagsus(irow,icol)= -one* &
       & cmplx(blkval(1,idir1,ipert1,idir2,ipert2,idir3,ipert3,iblok), &
       & blkval(2,idir1,ipert1,idir2,ipert2,idir3,ipert3,iblok),16)

       end do
     end do
   end do
 end do

!Calculate the Berry-curvature of the inverse magnetic susceptibility
 bc_ss=-matmul(invbarmagsus,matmul(bc_barmagsus,invbarmagsus)) 

 fac=2.714943600699**2/four !TMP
 open(10,file='g_ss.txt')
   do irow=1, ndim
     write(10,*) -ione*bc_ss(irow,1:ndim)*fac
   end do 
 close(10)


 call wrtout([ab_out,std_out], ' Berry curvature of the inverse spin susceptibility ')
 call wrtout([ab_out,std_out], '  atom1  dir  atom2  dir        Real              Imag')
 do irow=1, ndim
   do icol=1, ndim
     write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)' ) &
   & indexat(irow), cart(indexdir(irow)), indexat(icol), cart(indexdir(icol)), &
   & real(bc_ss(irow,icol)), aimag(bc_ss(irow,icol))
     call wrtout([ab_out,std_out], msg)
   end do
 end do
 call wrtout([ab_out,std_out], '   ')


 end subroutine berrycurv_ss
!!***

!!****f* m_ddb_magpen/berrycurv_sp
!! NAME
!! berrycurv_sp
!!
!! FUNCTION
!! Calculate the Berry curvature of the spin-phonon Hessian
!! (equivalent to the Berry curvature of the induced Zeeman fields at constrained
!! magnetic moments)
!!
!! INPUTS
!! blkval(2,3*mpert*3*mpert*3*mpert,nblok)=  Third-order derivative matrices
!!  In our case, the nblok is restricted to iblok
!! iblok= index of the current block
!! invbarmagsus(ndim,ndim)= Inverse of the penalized spin-sussceptibility tensor
!! mpatpol(2) = Atoms on which the magnetic penalty has been applied
!! mpdir(3) = Directions along which the spin-degrees of freedom have been stiffened 
!! mpert =maximum number of ipert
!! natom= number of atoms in unit cell
!! nblok= number of blocks in the DDB
!! nmdir= number of directions along which the magnetic penalty was applied
!! ndim= dimension of the square susceptibilities 
!! prtvol= control the volume of information written on output
!!
!! OUTPUT
!!
!! SOURCE

 subroutine berrycurv_sp(barmmom,bc_sp,bc_ss,blkval,iblok,invbarmagsus,&
& jblok,mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,prtvol,qphon,xred)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,jblok,mpert,natom,nblok,ndim,nmdir,prtvol
!arrays
 integer,intent(in) :: mpatpol(2),mpdir(3)
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,3,mpert,nblok)
 real(dp),intent(in) :: qphon(3),xred(3,natom)
 complex(dpc),intent(in) :: barmmom(ndim,(natom+2)*3)
 complex(dpc),intent(in) :: bc_ss(ndim,ndim)
 complex(dpc),intent(in) :: invbarmagsus(ndim,ndim)
 complex(dpc),intent(out) :: bc_sp(ndim,(natom+2)*3)
!Local variables -------------------------
!scalars
 integer :: iat1,iat2,icol,idir1,idir2,idir3,ipert1,ipert2,ipert3,irow
 integer :: ipert1_red,ipert2_red,idir1_red,idir2_red
 real(dp) :: fac,re,im
 complex(dpc), parameter :: ione=(0.d0,1.d0)
 character(len=1000) :: msg
!arrays
 integer(dp) :: indexat1(ndim),indexdir1(ndim)
 integer(dp) :: indexat2((natom+2)*3),indexdir2((natom+2)*3)
 complex(dpc) :: bc_barsp(ndim,(natom+2)*3)
 complex(dpc) :: bc_ps((natom+2)*3,ndim)
 character(len=1) :: cart(3)=(/'x','y','z'/)

! *********************************************************************

!Extract the berry curvature of the penalized moments
 bc_barsp=(zero,zero)
 ipert3= natom + 9
 idir3= 1
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
           bc_barsp(irow,icol)= -one* &
         & cmplx(blkval(1,idir1,ipert1,idir2,ipert2,idir3,ipert3,iblok), &
         & blkval(2,idir1,ipert1,idir2,ipert2,idir3,ipert3,iblok),16)
         else if (jblok /=0 .and. ipert2 == natom+2) then
           bc_barsp(irow,icol)= -one* &
         & cmplx(blkval(1,idir1,ipert1,idir2,ipert2,idir3,ipert3,jblok), &
         & blkval(2,idir1,ipert1,idir2,ipert2,idir3,ipert3,jblok),16)
         end if

       end do
     end do
   end do
 end do

 !Calculate the Berry curvature of the induced Zeeman fields
 bc_sp= -matmul(bc_ss,barmmom) - matmul(invbarmagsus,bc_barsp)

 do irow=1,ndim
   do iat1= 1, natom
     do idir1= 1, 3
       icol= (iat1-1)*3 + idir1
       bc_ps(icol,irow)=conjg(bc_sp(irow,icol)*exp(two_pi*(0.d0,1.d0)* dot_product(qphon,xred(:,iat1))))
     end do
   end do
 end do 

 fac=2.714943600699/two/0.52917 !TMP
 open(10,file='g_ps.txt')
 do irow=1,natom*3
   write(10,*) (0.d0,-1.d0)*bc_ps(irow,1:ndim)*fac
 end do 
 close(10)

 if (iblok /= 0) then
   call wrtout([ab_out,std_out], ' Berry curvature of the Zeeman fields induced by atomic displacements (at constrained magnetic moments)')
   call wrtout([ab_out,std_out], '  atom1  dir  efld.dir         Real              Imag')
   do irow=1, ndim
     do icol=1, natom*3
       write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)' ) &
     & indexat1(irow), cart(indexdir1(irow)), indexat2(icol), cart(indexdir2(icol)), &
     & real(bc_sp(irow,icol)), aimag(bc_sp(irow,icol))
       call wrtout([ab_out,std_out], msg)
     end do
   end do
   call wrtout([ab_out,std_out], '   ')
 end if
 if (jblok /= 0) then
   call wrtout([ab_out,std_out], ' Berry curvature of the Zeeman fields induced by electric field (at constrained magnetic moments)')
   call wrtout([ab_out,std_out], '  atom1  dir  efld.dir         Real              Imag')
   do irow=1, ndim
     do icol=(natom+1)*3+1, (natom+2)*3
       write(msg,'(i4,4x,a2,4x,a2,6x,2es18.9)' ) &
     & indexat1(irow), cart(indexdir1(irow)), cart(indexdir2(icol)), &
     & real(bc_sp(irow,icol)), aimag(bc_sp(irow,icol))
       call wrtout([ab_out,std_out], msg)
     end do
   end do
   call wrtout([ab_out,std_out], '   ')
 end if

 end subroutine berrycurv_sp
!!***

!!****f* m_ddb_magpen/berrycurv_pp
!! NAME
!! berrycurv_pp
!!
!! FUNCTION
!! Calculate the Berry curvature of the phonon-phonon Hessian
!! (at constrained magnetic moments)
!!
!! INPUTS
!! blkval(2,3*mpert*3*mpert*3*mpert,nblok)=  Third-order derivative matrices
!!  In our case, the nblok is restricted to iblok
!! iblok= index of the current block
!! invbarmagsus(ndim,ndim)= Inverse of the penalized spin-sussceptibility tensor
!! mpatpol(2) = Atoms on which the magnetic penalty has been applied
!! mpdir(3) = Directions along which the spin-degrees of freedom have been stiffened 
!! mpert =maximum number of ipert
!! natom= number of atoms in unit cell
!! nblok= number of blocks in the DDB
!! nmdir= number of directions along which the magnetic penalty was applied
!! ndim= dimension of the square susceptibilities 
!! prtvol= control the volume of information written on output
!!
!! OUTPUT
!!
!! SOURCE

 subroutine berrycurv_pp(barmagsus,bc_barmagsus,bc_sp,blkval,iblok,&
& mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,prtvol,qphon,xred,zfield)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,mpert,natom,nblok,ndim,nmdir,prtvol
!arrays
 integer,intent(in) :: mpatpol(2),mpdir(3)
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,3,mpert,nblok)
 real(dp),intent(in) :: qphon(3),xred(3,natom)
 complex(dpc),intent(in) :: barmagsus(ndim,ndim)
 complex(dpc),intent(in) :: bc_barmagsus(ndim,ndim)
 complex(dpc),intent(in) :: bc_sp(ndim,(natom+2)*3)
 complex(dpc),intent(in) :: zfield(ndim,(natom+2)*3)
!Local variables -------------------------
!scalars
 integer :: iat1,iat2,icol,idir1,idir2,idir3,ipert1,ipert2,ipert3,irow
 integer :: ipert1_red,ipert2_red,idir1_red,idir2_red
 real(dp) :: fac
 complex(dpc), parameter :: ione=(0.d0,1.d0)
 character(len=1000) :: msg
!arrays
 integer(dp) :: indexat1(ndim),indexdir1(ndim)
 integer(dp) :: indexat2((natom+2)*3),indexdir2((natom+2)*3)
 complex(dpc) :: bc_barpp(natom*3,natom*3), bc_pp(natom*3,natom*3) 
 complex(dpc) :: bc_pp_sf(natom*3,natom*3), term(natom*3,natom*3,3)
 complex(dpc) :: ifc_bc_sp(ndim,natom*3), ifc_zfield(ndim,natom*3)
 character(len=1) :: cart(3)=(/'x','y','z'/)

! *********************************************************************

!Extract the frequency derivative of the penalized IFCs 
 ipert3= natom + 9
 idir3= 1
 do ipert2= 1, natom
   do idir2= 1, 3
     icol=( ipert2-1)*3 + idir2
     do ipert1= 1, natom
       do idir1= 1, 3
         irow=( ipert1-1)*3 + idir1
         bc_barpp(irow,icol)= &
       & cmplx(blkval(1,idir1,ipert1,idir2,ipert2,idir3,ipert3,iblok), &
       &       blkval(2,idir1,ipert1,idir2,ipert2,idir3,ipert3,iblok),16)
       end do
     end do
   end do
 end do 

!Calculate the different terms entering the Berry curvature
 ifc_zfield(:,:)=zfield(:,1:natom*3)
 ifc_bc_sp(:,:)=bc_sp(:,1:natom*3)
 term(:,:,1)=matmul(transpose(conjg(ifc_bc_sp)),matmul(barmagsus,ifc_zfield))
 term(:,:,2)=matmul(transpose(conjg(ifc_zfield)),matmul(bc_barmagsus,ifc_zfield))
 term(:,:,3)=matmul(transpose(conjg(ifc_zfield)),matmul(barmagsus,bc_sp))

 bc_pp(:,:)= bc_barpp(:,:) + term(:,:,1) + term(:,:,2) + term(:,:,3)

 do ipert2= 1, natom
   do idir2= 1, 3
     icol=( ipert2-1)*3 + idir2
     do ipert1= 1, natom
       do idir1= 1, 3
         irow=( ipert1-1)*3 + idir1
         bc_pp_sf(irow,icol)=bc_pp(irow,icol)* &
       & exp(two_pi*(0.d0,1.d0)* dot_product(qphon,xred(:,ipert2)-xred(:,ipert1)))
       end do
     end do
   end do
 end do 

 fac=one/(0.52917)**2
 open(10,file='g_pp.txt')
 do irow=1,natom*3
   write(10,*) -ione*bc_pp_sf(irow,:)*fac
 end do 
 close(10)

 !Write the results
 call wrtout([ab_out,std_out], ' Berry curvature of interatomic force constants (at constrained magnetic moments)')
 call wrtout([ab_out,std_out], '  atom1  dir  atom2  dir        Real              Imag')
 do ipert1= 1, natom
   do idir1= 1, 3
     irow=( ipert1-1)*3 + idir1
     do ipert2= 1, natom
       do idir2= 1, 3
         icol=( ipert2-1)*3 + idir2
         write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)') &
       & ipert1, cart(idir1), ipert2, cart(idir2), &
       & real(bc_pp(irow,icol)), aimag(bc_pp(irow,icol))
         call wrtout([ab_out,std_out], msg)
       end do
     end do
     call wrtout([ab_out,std_out], '   ')
   end do
 end do 

 end subroutine berrycurv_pp

!!****f* m_ddb_magpen/berrycurv_tt
!! NAME
!! berrycurv_tt
!!
!! FUNCTION
!! Calculate the Berry curvature of the spin-spin Hessian
!!
!! INPUTS
!! blkval(2,3*mpert*3*mpert*3*mpert,nblok)=  Third-order derivative matrices
!!  In our case, the nblok is restricted to iblok
!! iblok= index of the current block
!! invbarmagsus(ndim,ndim)= Inverse of the penalized spin-sussceptibility tensor
!! mpatpol(2) = Atoms on which the magnetic penalty has been applied
!! mpdir(3) = Directions along which the spin-degrees of freedom have been stiffened 
!! mpert =maximum number of ipert
!! natom= number of atoms in unit cell
!! nblok= number of blocks in the DDB
!! nmdir= number of directions along which the magnetic penalty was applied
!! ndim= dimension of the square susceptibilities 
!! prtvol= control the volume of information written on output
!!
!! OUTPUT
!!
!! SOURCE

 subroutine berrycurv_tt(blkval,iblok,invbarmagsus,&
& mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,prtvol)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,mpert,natom,nblok,ndim,nmdir,prtvol
!arrays
 integer,intent(in) :: mpatpol(2),mpdir(3)
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,3,mpert,nblok)
 complex(dpc),intent(in) :: invbarmagsus(ndim,ndim)

 end subroutine berrycurv_tt
!!***

end module m_ddb_magpen
!!***
