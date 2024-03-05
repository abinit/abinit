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
 use m_macroave,        only : POLINT
 use m_io_tools,        only : open_file

 implicit none

 public :: ddb_magpen       ! Convert the derivatives calculated with the magnetic penalty into the physically relevant ones.
 public :: ddb_omega_interpol ! Perform an interpolation of the he derivatives calculated with the magnetic
                              ! penalty and later on convert them into the physically relevant ones at each value of interpolated omega.


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
 integer :: iblok,ii,ipert1,ipert2,jblok,kblok,lblok,nblok,ndim 
 integer :: nmat,nmdir,prtopt
 character(len=500) :: msg
 logical :: qeq0
!arrays
 integer :: rfelfd(4),rfphon(4),rfstrs(4),rfmagn(4),rffreq(4)
 real(dp) :: omega(3),qphnrm(3),qphon(3,3)
 complex(dpc) :: epsilon(3,3)
 complex(dpc), allocatable :: barmagsus(:,:),invbarmagsus(:,:)
 complex(dpc), allocatable :: invmagsus(:,:), magsus(:,:), invhmat(:,:)
 complex(dpc), allocatable :: barmmom(:,:),mmom(:,:), zfield(:,:)
 complex(dpc), allocatable :: bc_barmagsus(:,:),bc_ss(:,:),bc_sp(:,:)
 complex(dpc), allocatable :: ifcmat(:,:),zeff(:,:)

! *********************************************************************
 write(msg, '(2a,(80a),4a)' ) ch10,('=',ii=1,80),ch10,ch10,&
 ' Linear-response magnetic penalty section ',ch10
 call wrtout([std_out, ab_out], msg)

 prtopt=1
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
 ABI_MALLOC(ifcmat,(3*natom,3*natom))
 ABI_MALLOC(zeff,(3,3*natom))

 nblok=ddb%nblok
 do kblok=1,nblok

   ! Look for the spin-susceptibility block in the DDB
   omega=zero
   qphon=zero
   qphon(:,1)=ddb%qpt(1:3,kblok)
   qeq0=(sqrt(sum(qphon(:,1)**2))<tol8)
   qphnrm(:)=ddb%nrm(1,kblok)
   omega(1)=ddb%omega(1,kblok)
   rfphon(1:2)=0
   rfelfd(1:2)=0
   rfstrs(1:2)=0
   if (magpen<zero) then
     rfmagn(1:2)= 1
   else if (magpen>zero) then
     rfmagn(1:2)= 2
   end if

   write(msg, '(1a,(80a),2a,3f16.8,2a,f16.8,a)' ) ch10,('-',ii=1,80),ch10, &
   ' q point  ', qphon(:,1),ch10,&
   ' frequency', ddb%omega(1,kblok), ch10
   call wrtout([std_out, ab_out], msg)

   call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, &
  & mpatpol=mpatpol,mpdir=mpdir,omega=omega,rfmagn=rfmagn)

   ! Calculate and write the spin-susceptibility matrices
   if (iblok /= 0) then
     if (magpen<zero) then
       write(msg, '(2a)' ) ' Spin susceptibility (Uniform Zeeman) ',ch10
       call wrtout([std_out, ab_out], msg)
     else if (magpen>zero) then
       write(msg, '(2a)' ) ' Spin susceptibility (Local Zeeman) ',ch10
       call wrtout([std_out, ab_out], msg)
     end if

     call spinsus(barmagsus,ddb%val,iblok,invbarmagsus,invmagsus,invhmat,magpen,magsus,&
   & mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,prtopt,prtvol)

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
  & mpatpol=mpatpol,mpdir=mpdir,omega=omega,rfmagn=rfmagn)

   ! Then electric field
   ! Look for the induced magnetic moments block in the DDB
   jblok=0
   if (qeq0) then
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
   & mpatpol=mpatpol,mpdir=mpdir,omega=omega,rfmagn=rfmagn)
   end if

   if (iblok /= 0 .or. jblok /=0) then
     call magmom(barmmom,ddb%val,invbarmagsus,invhmat,iblok,jblok,magpen,magsus,mmom,&
   & mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,prtopt,prtvol,qphon,xred,zfield)
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
   call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, omega=omega)
   if (iblok /= 0 ) then
     call mp_ifc(barmagsus,ddb%val,iblok,ifcmat,magsus,mpert,mpopt,&
   & natom,nblok,ndim,prtopt,prtvol,qphon,xred,zfield)
   end if

   !Born effective charges block
   jblok=0
   if (qeq0) then
     rfphon(1:2)=1
     rfelfd(1:2)=2
     call ddb%get_block(jblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, omega=omega)
   end if
   if (jblok /= 0 ) then
     call mp_zeff(barmagsus,ddb%val,jblok,magsus,mpert,mpopt,&
   & natom,nblok,ndim,prtopt,prtvol,ucvol,zeff,zfield)
   end if

   !Dielectric susceptibility block
   lblok=0
   if (qeq0) then
     rfphon(:)=0
     rfelfd(1:2)=2
     call ddb%get_block(lblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, omega=omega)
   end if
   if (lblok /= 0 ) then
     call mp_diel(barmagsus,ddb%val,epsilon,lblok,magsus,mpert,mpopt,&
   & natom,nblok,ndim,prtopt,prtvol,ucvol,zfield)
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
     omega(:)=ddb_lw%omega(:,kblok)
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

     write(msg, '(a,(80a),a,3(a,3f16.8,a),3(a,f16.8,a))' ) ch10,('-',ii=1,80),ch10, &
     ' q point 1  ', qphon(:,1),ch10,&
     ' q point 2  ', qphon(:,2),ch10,&  
     ' q point 3  ', qphon(:,3),ch10,&  
     ' frequency 1', ddb_lw%omega(1,kblok), ch10,&
     ' frequency 2', ddb_lw%omega(2,kblok), ch10,&
     ' frequency 3', ddb_lw%omega(3,kblok), ch10
     call wrtout([std_out, ab_out], msg)

     call ddb_lw%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, 33, &
   & mpatpol=mpatpol,mpdir=mpdir,omega=omega,rfmagn=rfmagn,rffreq=rffreq)

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
   & mpatpol=mpatpol,mpdir=mpdir,omega=omega,rfmagn=rfmagn,rffreq=rffreq)

     ! Then electric field
     jblok=0
     if (qeq0) then
       rfphon(2)=0
       rfelfd(2)=2

       call ddb_lw%get_block(jblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, 33, &
     & mpatpol=mpatpol,mpdir=mpdir,omega=omega,rfmagn=rfmagn,rffreq=rffreq)
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
   & omega=omega,rffreq=rffreq)
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
 ABI_FREE(ifcmat)
 ABI_FREE(zeff)

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
& mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,prtopt,prtvol)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,mpert,natom,nblok,ndim,nmdir,prtopt,prtvol
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

 if (prtopt==1) then
  ! fac=2.714943600699**2*27.2114/four
   fac=27.2114/four
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
  
   end if !magpen>zero

 end if !prtopt
 
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
& mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,prtopt,prtvol,qphon,xred,zfield)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,jblok,mpert,natom,nblok,ndim,nmdir,prtopt,prtvol
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

 if (prtopt==1) then
  ! fac=2.714943600699/two*27.2114/0.529177
   fac=27.2114/0.529177/two
  
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
!! ifcmat(3*natom,3*natom)= IFC matrix calculated at the level of mpopt.
!!
!! SOURCE

 subroutine mp_ifc(barmagsus,blkval,iblok,ifcmat,magsus,mpert,mpopt,&
& natom,nblok,ndim,prtopt,prtvol,qphon,xred,zfield)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,mpert,mpopt,natom,nblok,ndim,prtopt,prtvol
!arrays
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
 real(dp),intent(in) :: qphon(3),xred(3,natom)
 complex(dpc),intent(in) :: barmagsus(ndim,ndim)
 complex(dpc),intent(in) :: magsus(ndim,ndim)
 complex(dpc),intent(in) :: zfield(ndim,(natom+2)*3)
 complex(dpc),intent(out) :: ifcmat(3*natom,3*natom)
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

 ifcmat= fmifc

 !Calculate the spin-relaxed flavor
 if (mpopt==2) then
   srifc= matmul(transpose(conjg(ifc_zfield)),matmul(magsus,ifc_zfield))
   srifc= fmifc - srifc

   ifcmat= srifc
 end if

 !Adopt the same phase convention as for the local Zeeman perturbation
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
 
 if (prtopt==1) then
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
!!  
!! epsilon(3,3)= Complex dielectric tensor calculated at the level of mpopt.
!!
!! SOURCE

 subroutine mp_diel(barmagsus,blkval,epsilon,iblok,magsus,mpert,mpopt,&
& natom,nblok,ndim,prtopt,prtvol,ucvol,zfield)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,mpert,mpopt,natom,nblok,ndim,prtopt,prtvol
 real(dp), intent(in) :: ucvol
!arrays
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
 complex(dpc),intent(in) :: barmagsus(ndim,ndim)
 complex(dpc),intent(out) :: epsilon(3,3)
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

 epsilon=fmepsilon

!Calculate the spin-relaxed flavor
 if (mpopt==2) then
   srdielsus= matmul(transpose(conjg(diel_zfield)),matmul(magsus,diel_zfield))
   srdielsus= fmdielsus - srdielsus
   srepsilon= -four_pi/ucvol*srdielsus
   srepsilon(1,1)= one + srepsilon(1,1)
   srepsilon(2,2)= one + srepsilon(2,2)
   srepsilon(3,3)= one + srepsilon(3,3)

   epsilon=srepsilon
 end if

 if (prtopt==1) then
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
!! zeff(3,natom*3)= Born effective charges calculated at the level of mpopt
!!
!! SOURCE

 subroutine mp_zeff(barmagsus,blkval,iblok,magsus,mpert,mpopt,&
& natom,nblok,ndim,prtopt,prtvol,ucvol,zeff,zfield)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,mpert,mpopt,natom,nblok,ndim,prtopt,prtvol
 real(dp), intent(in) :: ucvol
!arrays
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
 complex(dpc),intent(in) :: barmagsus(ndim,ndim)
 complex(dpc),intent(in) :: magsus(ndim,ndim)
 complex(dpc),intent(out) :: zeff(3,natom*3)
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

 zeff= fmzeff

 !Calculate the spin-relaxed flavor
 if (mpopt==2) then
   srzeff= matmul(transpose(conjg(diel_zfield)),matmul(magsus,ifc_zfield))
   srzeff= fmzeff - srzeff
 
   zeff= srzeff
 end if

 if (prtopt==1) then
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

! fac=2.714943600699**2/four !TMP
 fac=one/four !TMP
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

! fac=2.714943600699/two/0.52917 !TMP
 fac=one/two/0.52917 !TMP
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

!!****f* m_ddb_omega_interpol/ddb_omega_interpol
!! NAME
!! ddb_omega_interpol
!!
!! FUNCTION
!! Interpolate over frequency the secon-order derivatives calculated with 
!! the magnetic penalty and latter on convert them into physically relevant 
!! quantities.
!!
!! INPUTS
!! ddb (INOUT) = ddb block datastructure
!! magpen = amplitude (in Ha) of the applied magnetic penalty 
!! mpatpol(2) = Atoms on which the magnetic penalty has been applied
!! mpdir(3) = Directions along which the spin-degrees of freedom have been stiffened 
!! mpert = maximum number of ipert
!! mpopt = 1 calculate the frozen-magnetic second-order quantities
!!         2 calculate the spin-relaxed second-order quantities 
!! natom= number of atoms in unit cell
!! ntypat= number of atom types
!! ucvol= unit cell volume
!!
!! OUTPUT
!! ddb= ddb%val updated with the corrected second-order derivatives
!!
!! SOURCE

 subroutine ddb_omega_interpol(amu,ddb,outfilename_radix,magpen,mpatpol,mpdir,mpert,mpopt,natom, &
& nomega,ntypat,omegamax,omegamin,prtvol,rftyp,typat,ucvol,xred)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,mpopt,natom,nomega,ntypat,prtvol,rftyp
 real(dp),intent(in) :: magpen,omegamax,omegamin,ucvol
 character(len=*),intent(in) :: outfilename_radix
!arrays
 type(ddb_type),intent(inout) :: ddb
 integer,intent(in) :: mpatpol(2),mpdir(3),typat(natom)
 real(dp),intent(in) :: amu(ntypat)
 real(dp),intent(in) :: xred(3,natom)

!Local variables -------------------------
!scalars
 integer :: diel_unit,i,iblok,ii,imode,ipert1,ipert2,iw,j,jblok,kblok,lblok,mmom_unit,nblok,ndim 
 integer :: nmat,nmdir,nwcalc,phon_unit,prtopt,spin_unit,zeff_unit
 real(dp) :: omegastp
 character(len=5000) :: msg,pfmt
 character(len=fnlen) :: diel_filename,spin_filename,mmom_filename
 character(len=fnlen) :: phon_filename,zeff_filename
!arrays
 real(dp) :: qphnrm(3),qphon(3,3)
 real(dp), allocatable :: dint_barddb(:,:),int_barddb(:,:,:),omega(:),omegacalc(:)
 real(dp), allocatable :: eigvec(:,:,:,:,:),phfrq(:,:),phonspec(:)
 complex(dpc), allocatable :: barmagsus(:,:,:),invbarmagsus(:,:,:)
 complex(dpc), allocatable :: invmagsus(:,:,:), magsus(:,:,:), invhmat(:,:,:)
 complex(dpc), allocatable :: barmmom(:,:,:),mmom(:,:,:), zfield(:,:,:)
 complex(dpc), allocatable :: bc_barmagsus(:,:),bc_ss(:,:),bc_sp(:,:)
 complex(dpc), allocatable :: epsilon(:,:,:),ifcmat(:,:,:)
 complex(dpc), allocatable :: modemm(:,:,:),zeff(:,:,:),modezeff(:,:,:)

! *********************************************************************

 write(msg, '(2a,(80a),4a)' ) ch10,('=',ii=1,80),ch10,ch10,&
 ' Omega interpolation of magnetic penalty quantities section ',ch10
 call wrtout([std_out, ab_out], msg)

!Identify the calculated omegas
 nwcalc=ddb%nblok
 ABI_MALLOC(omegacalc,(nwcalc))
 omegacalc(:)=ddb%omega(1,:)

!Define the omega discretization
 omegastp=(omegamax-omegamin)/(nomega-1) 

 prtopt=0
 if (magpen<zero) then
   nmat= 1
 else if (magpen>zero) then
   nmat= mpatpol(2) - mpatpol(1) + 1
 end if
 nmdir=sum(mpdir(:))
 ndim=nmat*nmdir
 ABI_MALLOC(omega,(nomega))
 ABI_MALLOC(phfrq,(3*natom,nomega))
 ABI_MALLOC(phonspec,(nomega))
 ABI_MALLOC(eigvec,(2,3,natom,3,natom))
 ABI_MALLOC(modemm,(ndim,3*natom,nomega))
 ABI_MALLOC(zeff,(3,3*natom,nomega))
 ABI_MALLOC(modezeff,(3,3*natom,nomega))
 ABI_MALLOC(barmagsus,(ndim,ndim,nomega))
 ABI_MALLOC(magsus,(ndim,ndim,nomega))
 ABI_MALLOC(invbarmagsus,(ndim,ndim,nomega))
 ABI_MALLOC(invmagsus,(ndim,ndim,nomega))
 ABI_MALLOC(invhmat,(ndim,ndim,nomega))
 ABI_MALLOC(barmmom,(ndim,(natom+2)*3,nomega))
 ABI_MALLOC(mmom,(ndim,(natom+2)*3,nomega))
 ABI_MALLOC(zfield,(ndim,(natom+2)*3,nomega))
 ABI_MALLOC(epsilon,(3,3,nomega))
 ABI_MALLOC(ifcmat,(3*natom,3*natom,nomega))
 ABI_MALLOC(dint_barddb,(2,ddb%msize))
 ABI_MALLOC(int_barddb,(2,ddb%msize,1))

!Loop over the frequency
 do iw=1,nomega
   omega(iw)=omegamin+omegastp*(iw-1)

   do ii=1,ddb%msize
     if (all(ddb%flg(ii,:)==1)) then
       call POLINT(omegacalc,ddb%val(1,ii,:),nwcalc,omega(iw),int_barddb(1,ii,1),dint_barddb(1,ii)) 
       call POLINT(omegacalc,ddb%val(2,ii,:),nwcalc,omega(iw),int_barddb(2,ii,1),dint_barddb(2,ii)) 
     else if (count(ddb%flg(ii,:)==0)/=nwcalc) then
       write(msg,'(a,a,a)')&
       'ddb_omega_interpol detects differences between the DDB bloks for each frequency.',ch10,&
     & ' The interpolation has been stopped.' 
       ABI_ERROR(msg)
     end if
   end do 

   !Calculate the spin susceptibilities
   call spinsus(barmagsus(:,:,iw),int_barddb,1,invbarmagsus(:,:,iw),invmagsus(:,:,iw),&
 & invhmat(:,:,iw),magpen,magsus(:,:,iw),mpatpol,mpdir,mpert,natom,1,ndim,nmdir,prtopt,prtvol)

   !Calculate the magnetic moments
   call magmom(barmmom(:,:,iw),int_barddb,invbarmagsus(:,:,iw),invhmat(:,:,iw),1,1,magpen,magsus(:,:,iw),mmom(:,:,iw),&
 & mpatpol,mpdir,mpert,natom,1,ndim,nmdir,prtopt,prtvol,qphon,xred,zfield(:,:,iw))

   !Calculate the dielectric susceptibility
   call mp_diel(barmagsus(:,:,iw),int_barddb,epsilon(:,:,iw),1,magsus(:,:,iw),mpert,mpopt,&
 & natom,1,ndim,prtopt,prtvol,ucvol,zfield(:,:,iw))

   !Calculate the interatomic force constants
   call mp_ifc(barmagsus(:,:,iw),int_barddb,1,ifcmat(:,:,iw),magsus(:,:,iw),mpert,mpopt,&
 & natom,1,ndim,prtopt,prtvol,qphon,xred,zfield(:,:,iw))

   !Calculate the phonon Green's function and spectral function
   call phonon_green(amu,eigvec,ifcmat(:,:,iw),natom,ntypat,omega(iw),phfrq(:,iw),phonspec(iw),typat)

   !Calculate the mode-resolved magnetic moments
   call mode_mmom(amu,eigvec,mmom(:,:,iw),modemm(:,:,iw),natom,ndim,ntypat,typat)

   !Calculate the Born effective charges
   call mp_zeff(barmagsus(:,:,iw),int_barddb,1,magsus(:,:,iw),mpert,mpopt,&
 & natom,1,ndim,prtopt,prtvol,ucvol,zeff(:,:,iw),zfield(:,:,iw))

   !Calculate the mode-resolved Born effective charges
   call mode_zeff(amu,eigvec,natom,ntypat,typat,zeff(:,:,iw),modezeff(:,:,iw))

 end do

!!!  Print results of interpolation
!Spin susceptibilities
 spin_filename=trim(outfilename_radix)//"_SPINSUS"
 if (open_file(spin_filename, msg, newunit=spin_unit) /= 0) then
   ABI_ERROR(msg)
 end if

 write(spin_unit,*) '#'
 write(spin_unit,*) '#  Spin susceptibilities calculated and interpolated by ANADDB'
 write(spin_unit,*) '#'
 write(pfmt, '( "(es15.7, ", I4, "(es15.7))" )' )  ndim**2

 write(spin_unit,*) '#  Real part of local spin susceptibility tensor (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
 call wrtout(spin_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(magsus(i,j,iw)),j=1,ndim),i=1,ndim)
    call wrtout(spin_unit,msg,'COLL')
 end do

 write(spin_unit,*) ' '
 write(spin_unit,*) '#  Imaginary part of local spin susceptibility tensor (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
 call wrtout(spin_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(magsus(i,j,iw)),j=1,ndim),i=1,ndim)
    call wrtout(spin_unit,msg,'COLL')
 end do

 write(spin_unit,*) ' '
 write(spin_unit,*) '#  Real part of the inverse of the penalized local spin susceptibility tensor (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X^{-1}_11     X^{-1}_12     ...     X^{-1}_21     X^{-1}_22     ...'
 call wrtout(spin_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(invbarmagsus(i,j,iw)),j=1,ndim),i=1,ndim)
    call wrtout(spin_unit,msg,'COLL')
 end do

 write(spin_unit,*) ' '
 write(spin_unit,*) '#  Imaginary part of the inverse of the penalized local spin susceptibility tensor (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X^{-1}_11     X^{-1}_12     ...     X^{-1}_21     X^{-1}_22     ...'
 call wrtout(spin_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(invbarmagsus(i,j,iw)),j=1,ndim),i=1,ndim)
    call wrtout(spin_unit,msg,'COLL')
 end do

 close (spin_unit)

!Magnetic moments
 mmom_filename=trim(outfilename_radix)//"_MAGMOM"
 if (open_file(mmom_filename, msg, newunit=mmom_unit) /= 0) then
   ABI_ERROR(msg)
 end if

 write(mmom_unit,*) '#'
 write(mmom_unit,*) '#  Magnetic moments calculated and interpolated by ANADDB'
 write(mmom_unit,*) '#'

 write(pfmt, '( "(es15.7, ", I2, "(es17.7))" )' )  ndim
 do imode= 1, 3*natom
   write(mmom_unit,*) ' '
   write(mmom_unit,'(a,i3)') '#  Real part of magnetic moments (at. units) induced by phonon mode:', imode
   write(msg,'(a,a,a)') ch10,&
 &           ' # At  hw     m_{mat_1,1}     m_{mat_1,2}     ...     m_{mat_2,1}     m_{mat_2,2}'
   call wrtout(mmom_unit,msg,'COLL')
   do iw=1,nomega
     write(msg,pfmt) &
   & omega(iw), (real(modemm(i,imode,iw)),i=1,ndim)
     call wrtout(mmom_unit,msg,'COLL')
   end do
   write(mmom_unit,*) ' '
   write(mmom_unit,'(a,i3)') '#  Imaginary part of magnetic moments (at. units) induced by phonon mode:', imode
   write(msg,'(a,a,a)') ch10,&
 &           ' # At  hw     m_{mat_1,1}     m_{mat_1,2}     ...     m_{mat_2,1}     m_{mat_2,2}'
   call wrtout(mmom_unit,msg,'COLL')
   do iw=1,nomega
     write(msg,pfmt) &
   & omega(iw), (aimag(modemm(i,imode,iw)),i=1,ndim)
     call wrtout(mmom_unit,msg,'COLL')
   end do
 end do

 write(pfmt, '( "(es15.7, ", I2, "(es17.7))" )' )  ndim*3
 write(mmom_unit,*) ' '
 write(mmom_unit,*) '#  Real part of magnetic moments induced by electric field (at. units)'
 write(msg,'(a,a,a)') ch10,&
&           ' # At  hw     m_{mat_1,1}^{Ex}     m_{mat_1,1}^{Ey}',&
&           '      ...     m_{mat_1,2}^{Ex}     ...     m_{mat_2,1}^{Ex}     ...'
 call wrtout(mmom_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(mmom(i,j,iw)),j=3*(natom+1)+1,3*(natom+2)),i=1,ndim)
    call wrtout(mmom_unit,msg,'COLL')
 end do

 write(mmom_unit,*) ' '
 write(mmom_unit,*) '#  Imaginary part of magnetic moments induced by electric field (at. units)'
 write(msg,'(a,a,a)') ch10,&
&           ' # At  hw     m_{mat_1,1}^{Ex}     m_{mat_1,1}^{Ey}',&
&           '      ...     m_{mat_1,2}^{Ex}     ...     m_{mat_2,1}^{Ex}     ...'
 call wrtout(mmom_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(mmom(i,j,iw)),j=3*(natom+1)+1,3*(natom+2)),i=1,ndim)
    call wrtout(mmom_unit,msg,'COLL')
 end do

 close(mmom_unit)

!Dielectric susceptibility
 diel_filename=trim(outfilename_radix)//"_DIELSUS"
 if (open_file(diel_filename, msg, newunit=diel_unit) /= 0) then
   ABI_ERROR(msg)
 end if

 write(diel_unit,*) '#'
 if (mpopt==1) then
   write(diel_unit,*) '#  Frozen-magnetic dielectric tensor calculated and interpolated by ANADDB'
 else if (mpopt==2) then
   write(diel_unit,*) '#  Spin-relaxed dielectric tensor calculated and interpolated by ANADDB'
 else
   write(msg,'(a)') 'ddb_omega_interpol: variable mpopt just can be 1 or 2'
   ABI_ERROR(msg)
 end if
 
 write(diel_unit,*) '#'
 write(pfmt, '( "(es15.7, ", I4, "(es15.7))" )' ) 9 

 write(diel_unit,*) '#  Real part of dielectric tensor'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     eps_11     eps_12     ...     eps_21     eps_22     ...'
 call wrtout(diel_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(epsilon(i,j,iw)),j=1,3),i=1,3)
    call wrtout(diel_unit,msg,'COLL')
 end do

 write(diel_unit,*) ' '
 write(diel_unit,*) '#  Imaginary part of dielectric tensor'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     eps_11     eps_12     ...     eps_21     eps_22     ...'
 call wrtout(diel_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(epsilon(i,j,iw)),j=1,3),i=1,3)
    call wrtout(diel_unit,msg,'COLL')
 end do

 close(diel_unit)

!Phonon spectral function
 phon_filename=trim(outfilename_radix)//"_SPECTRAL_PHONON"
 if (open_file(phon_filename, msg, newunit=phon_unit) /= 0) then
   ABI_ERROR(msg)
 end if

 write(phon_unit,*) '#'
 if (mpopt==1) then
   write(phon_unit,*) '#  Frozen-magnetic phonon spectral function calculated and interpolated by ANADDB'
 else if (mpopt==2) then
   write(phon_unit,*) '#  Spin-relaxed phonon spectral function calculated and interpolated by ANADDB'
 else
   write(msg,'(a)') 'ddb_omega_interpol: variable mpopt just can be 1 or 2'
   ABI_ERROR(msg)
 end if
 
 write(msg,'(a,a)') ch10, ' # At  hw'
 call wrtout(phon_unit,msg,'COLL')
 do iw=1,nomega
   write(msg,*) omega(iw), phonspec(iw)
   call wrtout(phon_unit,msg,'COLL')
 end do

 close(phon_unit)

!Phonon frequencies
 phon_filename=trim(outfilename_radix)//"_PHFRQ"
 if (open_file(phon_filename, msg, newunit=phon_unit) /= 0) then
   ABI_ERROR(msg)
 end if

 write(phon_unit,*) '#'
 if (mpopt==1) then
   write(phon_unit,*) '#  Frozen-magnetic phonon frequencies calculated and interpolated by ANADDB'
 else if (mpopt==2) then
   write(phon_unit,*) '#  Spin-relaxed phonon frequencies calculated and interpolated by ANADDB'
 else
   write(msg,'(a)') 'ddb_omega_interpol: variable mpopt just can be 1 or 2'
   ABI_ERROR(msg)
 end if

 write(pfmt, '( "(es15.7, ", I4, "(es15.7))" )' ) natom*3
 write(msg,'(a,a)') ch10,&
&           ' # At  hw    eval(1)     eval(2) ...'
 call wrtout(phon_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) omega(iw), phfrq(:,iw)
    call wrtout(phon_unit,msg,'COLL')
 end do
 
 close(phon_unit)

!Born effective charges
 zeff_filename=trim(outfilename_radix)//"_ZEFF"
 if (open_file(zeff_filename, msg, newunit=zeff_unit) /= 0) then
   ABI_ERROR(msg)
 end if

 write(zeff_unit,*) '#'
 write(zeff_unit,*) '#  Born effective charges calculated and interpolated by ANADDB'
 write(zeff_unit,*) '#'

 write(pfmt, '( "(es15.7, ", I2, "(es17.7))" )' ) 3 
 do imode= 1, 3*natom
   write(zeff_unit,*) ' '
   write(zeff_unit,'(a,i3)') '#  Real part of Born charge (at. units) induced by phonon mode:', imode
   write(msg,'(a,a)') ch10,&
 &           ' # At  hw     Z^x_{n}     Z^y_{n}     Z^z_{n}'
   call wrtout(zeff_unit,msg,'COLL')
   do iw=1,nomega
     write(msg,pfmt) &
   & omega(iw), (real(modezeff(i,imode,iw)),i=1,3)
     call wrtout(zeff_unit,msg,'COLL')
   end do
   write(zeff_unit,*) ' '
   write(zeff_unit,'(a,i3)') '#  Imaginary part of Born charge (at. units) induced by phonon mode:', imode
   write(msg,'(a,a)') ch10,&
 &           ' # At  hw     Z^x_{n}     Z^y_{n}     Z^z_{n}'
   call wrtout(zeff_unit,msg,'COLL')
   do iw=1,nomega
     write(msg,pfmt) &
   & omega(iw), (aimag(modezeff(i,imode,iw)),i=1,3)
     call wrtout(zeff_unit,msg,'COLL')
   end do
 end do

 close(zeff_unit)

 ABI_FREE(dint_barddb)
 ABI_FREE(int_barddb)
 ABI_FREE(barmagsus)
 ABI_FREE(magsus)
 ABI_FREE(invbarmagsus)
 ABI_FREE(invmagsus)
 ABI_FREE(invhmat)
 ABI_FREE(barmmom)
 ABI_FREE(mmom)
 ABI_FREE(zfield)
 ABI_FREE(epsilon)
 ABI_FREE(ifcmat)
 ABI_FREE(omega)
 ABI_FREE(phfrq)
 ABI_FREE(phonspec)
 ABI_FREE(eigvec)
 ABI_FREE(modemm)
 ABI_FREE(zeff)
 ABI_FREE(modezeff)

 end subroutine ddb_omega_interpol
!!***

!!****f* ABINIT/phonon_green
!! NAME
!!  phonon_green
!!
!! FUNCTION
!!  Computes the phonon Green's function and spectral function at 
!!  a given value of frequency and imaginary damping
!!
!! COPYRIGHT
!!  Copyright (C) 2024 ABINIT group (FIXME: add author)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  amu(ntypat)= atomic masses
!!  ifc(3*natom,3*natom)= Interatomic-force constants calculated at a given omega
!!  natom= number of atoms in the cell
!!  ntypat= number of atom types in the cell
!!  omega= frequency at which the IFCs have been calculated
!!  typat(natom)= array with the type of atoms in the cell
!!
!! OUTPUT
!!  phonspec= phonon spectral function at the input omega
!!  phfrq= phonon frequencies calculated with the IFCs at the input omega
!!  phfrq= phonon eigenvectors calculated with the IFCs at the input omega
!!  
!!
!! SIDE EFFECTS
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


subroutine phonon_green(amu,eigvec,ifc,natom,ntypat,omega,phfrq,phonspec,typat)

 use defs_basis
 use m_errors
 use m_profiling_abi

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in)  :: natom,ntypat 
 real(dp), intent(in) :: omega
 real(dp), intent(out) :: phonspec
!arrays
 integer, intent(in) :: typat(natom)
 real(dp), intent(in) :: amu(ntypat)
 real(dp),intent(out) :: eigvec(2*3*natom*3*natom)
 real(dp),intent(out) :: phfrq(3*natom)
 complex(dpc), intent(in) :: ifc(3*natom,3*natom)

!Local variables-------------------------------
!scalars
 integer :: iat1,iat2,idir1,idir2,icol,ier,imode,info,irow,lwork,ndim
 complex(dpc) :: eta
!arrays
 integer, allocatable :: ipiv(:)
 real(dp), allocatable :: massfac(:,:)
 real(dp), allocatable :: matrx(:,:),zhpev1(:,:),zhpev2(:)
 real(dp), allocatable :: eigval(:)
 complex(dpc), allocatable :: dynmat(:,:),w2dynmat(:,:)
 complex(dpc),allocatable :: work(:),work1(:,:)
!character(len=500) :: msg                   

! *************************************************************************

 DBG_ENTER("COLL")

!Build an array with the mass factors
 ABI_MALLOC(massfac,(natom,natom))
 do iat2= 1, natom
   do iat1= 1, natom
     massfac(iat1,iat2)=one/sqrt(amu(typat(iat1))*amu(typat(iat2)))
   end do
 end do

!Build the ((w+eta)**2 - D(w)) matrix 
 ndim=3*natom
 ABI_MALLOC(dynmat,(ndim,ndim))
 ABI_MALLOC(w2dynmat,(ndim,ndim))
 eta=(0.0_dp,0.00001_dp)
 do iat2= 1, natom
   do idir2= 1, 3
     icol= (iat2-1)*3 + idir2
     do iat1= 1, natom
       do idir1= 1, 3
         irow= (iat1-1)*3 + idir1
         dynmat(irow,icol)= massfac(iat1,iat2)*ifc(irow,icol)
         w2dynmat(irow,icol)= -one*dynmat(irow,icol)
         if (irow==icol) then
           w2dynmat(irow,icol)= (omega+eta)**2 + w2dynmat(irow,icol)
         end if
       end do
     end do
   end do
 end do

!Invert to obtain the phonon Green's function
 ABI_MALLOC(work1,(ndim,ndim))
 work1=w2dynmat

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


!Finally extract the spectral function from the trace
 phonspec= zero
 do irow= 1, ndim
   phonspec= phonspec + aimag(work1(irow,irow))
 end do
 phonspec= -two*omega/pi * phonspec

!Diagonalize the Dynamical matrix
 ABI_MALLOC(matrx,(2,(3*natom*(3*natom+1))/2))
 ABI_MALLOC(eigval,(ndim))
 do icol= 1, ndim
   do irow= 1, icol
     matrx(1,irow + (icol-1)*icol/2)=real(dynmat(irow,icol))
     matrx(2,irow + (icol-1)*icol/2)=aimag(dynmat(irow,icol))
   end do
 end do 

 ABI_MALLOC(zhpev1,(2,2*3*natom-1))
 ABI_MALLOC(zhpev2,(3*3*natom-2))

 call ZHPEV ('V','U',3*natom,matrx,eigval,eigvec,3*natom,zhpev1,zhpev2,ier)
 ABI_CHECK(ier == 0, sjoin('zhpev returned:', itoa(ier)))

 ABI_FREE(matrx)
 ABI_FREE(zhpev1)
 ABI_FREE(zhpev2)

 ! Get the phonon frequencies (negative by convention, if the eigenvalue of the dynamical matrix is negative)
 do imode=1,3*natom
   if(eigval(imode)>=1.0d-16)then
     phfrq(imode)=sqrt(eigval(imode))
   else if(eigval(imode)>=-1.0d-16)then
     phfrq(imode)=zero
   else
     phfrq(imode)=-sqrt(-eigval(imode))
   end if
 end do


 ABI_FREE(ipiv)
 ABI_FREE(work1)
 ABI_FREE(massfac)
 ABI_FREE(dynmat)

 DBG_EXIT("COLL")

end subroutine phonon_green
!!***

!!****f* ABINIT/mode_mmom
!! NAME
!!  mode_mmom
!!
!! FUNCTION
!!  Projects the magnetic moments on the eigenmodes of the dynamical 
!!  matrix calculated at each value of omega
!!
!! COPYRIGHT
!!  Copyright (C) 2024 ABINIT group (FIXME: add author)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  amu(ntypat)= atomic masses
!!  eigvec(2,3,natom,3,natom)= dynamical matrix eigenvectors
!!  mmom(ndim,(natom+2)*3)= first-order magnetic moments
!!  natom= number of atoms in the cell
!!  ndim= dimension of the penalized degrees of freedom
!!  ntypat= number of atom types in the cell
!!  typat(natom)= array with the type of atoms in the cell
!!
!! OUTPUT
!!  modemm(ndim,3*natom)= mode-resolved magnetic moments
!!
!! SIDE EFFECTS
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


subroutine mode_mmom(amu,eigvec,mmom,modemm,natom,ndim,ntypat,typat)

 use defs_basis
 use m_errors
 use m_profiling_abi

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in)  :: natom,ndim,ntypat 
!arrays
 integer, intent(in) :: typat(natom)
 real(dp), intent(in) :: amu(ntypat)
 real(dp), intent(in) :: eigvec(2,3,natom,3,natom)
 complex(dpc), intent(in) :: mmom(ndim,(natom+2)*3)
 complex(dpc), intent(out) :: modemm(ndim,3*natom)

!Local variables-------------------------------
!scalars
 integer :: iat1,iat2,idir1,idir2,icol,im,imode,irow
 real(dp) :: mcell
!arrays
 real(dp), allocatable :: mass(:)
!character(len=500) :: msg                   

! *************************************************************************

 DBG_ENTER("COLL")

!Define the mass factors
 ABI_MALLOC(mass,(natom))
 mcell=zero
 do iat1= 1, natom
   mass(iat1)= amu(typat(iat1))
   mcell= mcell + amu(typat(iat1))
 end do
 mass(:)=sqrt(mcell/mass(:))

!Compute the mode-resolved moments
 modemm(:,:)=(zero,zero)
 do im= 1, ndim
   do iat2= 1, natom
     do idir2= 1, 3
       imode= (iat2-1)*3 + idir2
       do iat1= 1, natom
         do idir1= 1, 3
           irow= (iat1-1)*3 + idir1
           modemm(im,imode)= modemm(im,imode) +  mass(iat1)*mmom(im,irow)* &
         & cmplx(eigvec(1,idir1,iat1,idir2,iat2),eigvec(2,idir1,iat1,idir2,iat2),16)
         end do
       end do
     end do
   end do
 end do

 ABI_FREE(mass)

 DBG_EXIT("COLL")

end subroutine mode_mmom
!!***

!!****f* ABINIT/mode_zeff
!! NAME
!!  mode_zeff
!!
!! FUNCTION
!!  Projects the Born effective charges on the eigenmodes of the dynamical 
!!  matrix calculated at each value of omega
!!
!! COPYRIGHT
!!  Copyright (C) 2024 ABINIT group (FIXME: add author)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  amu(ntypat)= atomic masses
!!  eigvec(2,3,natom,3,natom)= dynamical matrix eigenvectors
!!  zeff(ndim,(natom+2)*3)= first-order magnetic moments
!!  natom= number of atoms in the cell
!!  ntypat= number of atom types in the cell
!!  typat(natom)= array with the type of atoms in the cell
!!  zeff(3,3*natom)= atomic Born effective charges at a given omega
!!
!! OUTPUT
!!  modezeff(3,3*natom)= mode-resolved magnetic moments
!!
!! SIDE EFFECTS
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


subroutine mode_zeff(amu,eigvec,natom,ntypat,typat,zeff,modezeff)

 use defs_basis
 use m_errors
 use m_profiling_abi

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in)  :: natom,ntypat 
!arrays
 integer, intent(in) :: typat(natom)
 real(dp), intent(in) :: amu(ntypat)
 real(dp), intent(in) :: eigvec(2,3,natom,3,natom)
 complex(dpc), intent(in) :: zeff(3,natom*3)
 complex(dpc), intent(out) :: modezeff(3,3*natom)

!Local variables-------------------------------
!scalars
 integer :: iat1,iat2,idir1,idir2,icol,im,imode,irow
 real(dp) :: mcell
!arrays
 real(dp), allocatable :: mass(:)
!character(len=500) :: msg                   

! *************************************************************************

 DBG_ENTER("COLL")

!Define the mass factors
 ABI_MALLOC(mass,(natom))
 mcell=zero
 do iat1= 1, natom
   mass(iat1)= amu(typat(iat1))
   mcell= mcell + amu(typat(iat1))
 end do
 mass(:)=sqrt(mcell/mass(:))

!Compute the mode-resolved Born charges
 modezeff(:,:)=(zero,zero)
 do im= 1, 3
   do iat2= 1, natom
     do idir2= 1, 3
       imode= (iat2-1)*3 + idir2
       do iat1= 1, natom
         do idir1= 1, 3
           irow= (iat1-1)*3 + idir1
           modezeff(im,imode)= modezeff(im,imode) +  mass(iat1)*zeff(im,irow)* &
         & cmplx(eigvec(1,idir1,iat1,idir2,iat2),eigvec(2,idir1,iat1,idir2,iat2),16)
         end do
       end do
     end do
   end do
 end do

 ABI_FREE(mass)

 DBG_EXIT("COLL")

end subroutine mode_zeff
!!***

end module m_ddb_magpen
!!***
