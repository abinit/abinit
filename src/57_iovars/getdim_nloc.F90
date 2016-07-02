!{\src2tex{textfont=tt}}
!!****f* ABINIT/getdim_nloc
!! NAME
!! getdim_nloc
!!
!! FUNCTION
!! Determine the dimensions of arrays that contain
!! the definition of non-local projectors : ekb, ffspl, indlmn
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  mixalch(npspalch,ntypalch,nimage)=alchemical mixing coefficients
!!  nimage=number of images
!!  npsp=number of pseudopotentials
!!  npspalch=number of pseudopotentials for alchemical purposes
!!  ntypat=number of types of pseudo atoms
!!  ntypalch=number of types of alchemical pseudo atoms
!!  pspheads(npsp)=<type pspheader_type>all the important information from the
!!   pseudopotential file headers, as well as the psp file names
!!
!! OUTPUT
!!  lmnmax=maximum number of l,m,n projectors, not taking into account the spin-orbit
!!  lmnmaxso=maximum number of l,m,n projectors, taking into account the spin-orbit
!!  lnmax=maximum number of l,n projectors, not taking into account the spin-orbit
!!  lnmaxso=maximum number of l,n projectors, taking into account the spin-orbit
!!
!! PARENTS
!!      m_psps,memory_eval
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine getdim_nloc(lmnmax,lmnmaxso,lnmax,lnmaxso,mixalch,nimage,npsp,npspalch,&
& ntypat,ntypalch,pspheads)

 use defs_basis
 use defs_datatypes
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getdim_nloc'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nimage,npsp,npspalch,ntypalch,ntypat
 integer,intent(out) :: lmnmax,lmnmaxso,lnmax,lnmaxso
!arrays
 real(dp),intent(in) :: mixalch(npspalch,ntypalch,nimage)
 type(pspheader_type),intent(in) :: pspheads(npsp)

!Local variables-------------------------------
!scalars
 integer :: ilang,ipsp,ipspalch,itypalch,itypat,ntyppure
!integer :: llmax
 character(len=500) :: message
!arrays
 integer,allocatable :: lmnproj_typat(:),lmnprojso_typat(:),lnproj_typat(:)
 integer,allocatable :: lnprojso_typat(:),nproj_typat(:,:),nprojso_typat(:,:)

! *************************************************************************

!write(std_out,*)' getdim_nloc: 'pspheads(1)%nproj(0:3)=',pspheads(1)%nproj(0:3)

 ABI_ALLOCATE(lmnproj_typat,(ntypat))
 ABI_ALLOCATE(lmnprojso_typat,(ntypat))
 ABI_ALLOCATE(lnproj_typat,(ntypat))
 ABI_ALLOCATE(lnprojso_typat,(ntypat))
 ABI_ALLOCATE(nproj_typat,(0:3,ntypat))
 ABI_ALLOCATE(nprojso_typat,(3,ntypat))
 lmnproj_typat(:)=0 ; lmnprojso_typat(:)=0
 lnproj_typat(:)=0 ; lnprojso_typat(:)=0
 nproj_typat(:,:)=0 ; nprojso_typat(:,:)=0

 ntyppure=ntypat-ntypalch

!For each type of pseudo atom, compute the number of projectors
!First, pure pseudo atoms
 if(ntyppure>0)then
   do itypat=1,ntyppure
     nproj_typat(0:3,itypat)=pspheads(itypat)%nproj(0:3)
     nprojso_typat(:,itypat)=pspheads(itypat)%nprojso(:)
   end do
 end if

!Then, alchemical pseudo atoms
 if(ntypalch>0)then
   do itypat=ntyppure+1,ntypat
     itypalch=itypat-ntyppure
     do ipsp=ntyppure+1,npsp
       ipspalch=ipsp-ntyppure
!      If there is some mixing, must accumulate the projectors
       if(sum(abs(mixalch(ipspalch,itypalch,:)))>tol10)then
         nproj_typat(0:3,itypat)=nproj_typat(0:3,itypat)+pspheads(ipsp)%nproj(0:3)
         nprojso_typat(:,itypat)=nprojso_typat(:,itypat)+pspheads(ipsp)%nprojso(:)
       end if
     end do
   end do
 end if

!Now that the number of projectors is known, accumulate the dimensions
 do itypat=1,ntypat
   do ilang=0,3
     lnproj_typat(itypat)=lnproj_typat(itypat)+nproj_typat(ilang,itypat)
     lmnproj_typat(itypat)=lmnproj_typat(itypat)+nproj_typat(ilang,itypat)*(2*ilang+1)
   end do
   lnprojso_typat(itypat)=lnproj_typat(itypat)
   lmnprojso_typat(itypat)=lmnproj_typat(itypat)
   do ilang=1,3
     lnprojso_typat(itypat)=lnprojso_typat(itypat)+nprojso_typat(ilang,itypat)
     lmnprojso_typat(itypat)=lmnprojso_typat(itypat)+nprojso_typat(ilang,itypat)*(2*ilang+1)
   end do
 end do

!Compute the maximal bounds, at least equal to 1, even for local psps
 lmnmax=1;lmnmaxso=1;lnmax=1;lnmaxso=1
 do itypat=1,ntypat
   lmnmax  =max(lmnmax  ,lmnproj_typat  (itypat))
   lmnmaxso=max(lmnmaxso,lmnprojso_typat(itypat))
   lnmax   =max(lnmax   ,lnproj_typat   (itypat))
   lnmaxso =max(lnmaxso ,lnprojso_typat (itypat))
 end do
!The initial coding (below) was not totally portable (MT 110215)
!lmnmax=max(maxval(lmnproj_typat(1:ntypat)),1)
!lmnmaxso=max(maxval(lmnprojso_typat(1:ntypat)),1)
!lnmax=max(maxval(lnproj_typat(1:ntypat)),1)
!lnmaxso=max(maxval(lnprojso_typat(1:ntypat)),1)

 if(maxval(lmnproj_typat(1:ntypat))==0)then
   write(message, '(3a)' )&
&   'Despite there is only a local part to pseudopotential(s),',ch10,&
&   'lmnmax and lnmax are set to 1.'
   MSG_COMMENT(message)
 end if

!XG040806 : These lines make modifications of lnmax and lmnmax
!that are unjustified in many cases, according to the many tests cases
!where they produce a changes, while the test case was working properly.
!One should understand better the needs, and code more appropriate changes ...
!lnmax/lmnmax has to be bigger than 1+lmax (for compatibility reasons)
!llmax=maxval(pspheads(1:ntypat)%lmax)+1 ! And this line might have trouble with HP compiler
!if (lnmax   <llmax) lnmax=llmax
!if (lnmaxso <llmax) lnmaxso=llmax
!if (lmnmax  <llmax) lmnmax=llmax
!if (lmnmaxso<llmax) lmnmaxso=llmax

 write(message, '(a,a,i4,a,i4,3a,i4,a,i4,a)' ) ch10,&
& ' getdim_nloc : deduce lmnmax  =',lmnmax,', lnmax  =',lnmax,',',ch10,&
& '                      lmnmaxso=',lmnmaxso,', lnmaxso=',lnmaxso,'.'
 call wrtout(std_out,message,'COLL')

 ABI_DEALLOCATE(lmnproj_typat)
 ABI_DEALLOCATE(lmnprojso_typat)
 ABI_DEALLOCATE(lnproj_typat)
 ABI_DEALLOCATE(lnprojso_typat)
 ABI_DEALLOCATE(nproj_typat)
 ABI_DEALLOCATE(nprojso_typat)

end subroutine getdim_nloc
!!***
