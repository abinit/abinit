!{\src2tex{textfont=tt}}
!!****f* ABINIT/paw2wvl
!! NAME
!!  paw2wvl
!!
!! FUNCTION
!!  Points WVL objects to PAW objects
!!
!! COPYRIGHT
!!  Copyright (C) 2011-2016 ABINIT group (T. Rangel)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine paw2wvl(pawtab,proj,wvl)

 use m_profiling_abi
 use m_errors
 use defs_basis
 use m_pawtab, only : pawtab_type
 use defs_wvltypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'paw2wvl'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(pawtab_type),intent(in)::pawtab(:)
 type(wvl_internal_type), intent(inout) :: wvl
 type(wvl_projectors_type),intent(inout)::proj

!Local variables-------------------------------
#if defined HAVE_DFT_BIGDFT
 integer :: ib,ig
 integer :: itypat,jb,ll,lmnmax,lmnsz,nn,ntypat,ll_,nn_
 integer ::maxmsz,msz1,ptotgau,max_lmn2_size
 logical :: test_wvl
 real(dp) :: a1
 integer,allocatable :: msz(:)
#endif

!extra variables, use to debug
!  integer::nr,unitp,i_shell,ii,ir,ng
!  real(dp)::step,rmax
!  real(dp),allocatable::r(:), y(:)
!  complex::fac,arg
!  complex(dp),allocatable::f(:),g(:,:)

! *********************************************************************

!DEBUG
!write (std_out,*) ' paw2wvl : enter'
!ENDDEBUG

#if defined HAVE_DFT_BIGDFT

 ntypat=size(pawtab)

!wvl object must be allocated in pawtab
 if (ntypat>0) then
   test_wvl=.true.
   do itypat=1,ntypat
     if (pawtab(itypat)%has_wvl==0.or.(.not.associated(pawtab(itypat)%wvl))) test_wvl=.false.
   end do
   if (.not.test_wvl) then
     MSG_BUG('pawtab%wvl must be allocated!')
   end if
 end if

!Find max mesh size
 ABI_ALLOCATE(msz,(ntypat))
 do itypat=1,ntypat
   msz(itypat)=pawtab(itypat)%wvl%rholoc%msz
 end do
 maxmsz=maxval(msz(1:ntypat))

 ABI_ALLOCATE(wvl%rholoc%msz,(ntypat))
 ABI_ALLOCATE(wvl%rholoc%d,(maxmsz,4,ntypat))
 ABI_ALLOCATE(wvl%rholoc%rad,(maxmsz,ntypat))
 ABI_ALLOCATE(wvl%rholoc%radius,(ntypat))

 do itypat=1,ntypat
   msz1=pawtab(itypat)%wvl%rholoc%msz
   msz(itypat)=msz1
   wvl%rholoc%msz(itypat)=msz1
   wvl%rholoc%d(1:msz1,:,itypat)=pawtab(itypat)%wvl%rholoc%d(1:msz1,:)
   wvl%rholoc%rad(1:msz1,itypat)=pawtab(itypat)%wvl%rholoc%rad(1:msz1)
   wvl%rholoc%radius(itypat)=pawtab(itypat)%rpaw
 end do
 ABI_DEALLOCATE(msz)
!
!Now fill projectors type:
!
 ABI_DATATYPE_ALLOCATE(proj%G,(ntypat))
!
!nullify all for security:
!do itypat=1,ntypat
!call nullify_gaussian_basis(proj%G(itypat))
!end do
 proj%G(:)%nat=1 !not used
 proj%G(:)%ncplx=2 !Complex gaussians
!
!Obtain dimensions:
!jb=0
!ptotgau=0
 do itypat=1,ntypat
!  do ib=1,pawtab(itypat)%basis_size
!  jb=jb+1
!  end do
!  ptotgau=ptotgau+pawtab(itypat)%wvl%ptotgau
   ptotgau=pawtab(itypat)%wvl%ptotgau
   proj%G(itypat)%nexpo=ptotgau
   proj%G(itypat)%nshltot=pawtab(itypat)%basis_size
 end do
!proj%G%nexpo=ptotgau
!proj%G%nshltot=jb
!
!Allocations
 do itypat=1,ntypat
   ABI_ALLOCATE(proj%G(itypat)%ndoc ,(proj%G(itypat)%nshltot))
   ABI_ALLOCATE(proj%G(itypat)%nam  ,(proj%G(itypat)%nshltot))
   ABI_ALLOCATE(proj%G(itypat)%xp   ,(proj%G(itypat)%ncplx,proj%G(itypat)%nexpo))
   ABI_ALLOCATE(proj%G(itypat)%psiat,(proj%G(itypat)%ncplx,proj%G(itypat)%nexpo))
 end do

!jb=0
 do itypat=1,ntypat
   proj%G(itypat)%ndoc(:)=pawtab(itypat)%wvl%pngau(:)
 end do
!
 do itypat=1,ntypat
   proj%G(itypat)%xp(:,:)=pawtab(itypat)%wvl%parg(:,:)
   proj%G(itypat)%psiat(:,:)=pawtab(itypat)%wvl%pfac(:,:)
 end do

!Change the real part of psiat to the form adopted in BigDFT:
!Here we use exp{(a+ib)x^2}, where a is a negative number.
!In BigDFT: exp{-0.5(x/c)^2}exp{i(bx)^2}, and c is positive
!Hence c=sqrt(0.5/abs(a))
 do itypat=1,ntypat
   do ig=1,proj%G(itypat)%nexpo
     a1=proj%G(itypat)%xp(1,ig)
     a1=(sqrt(0.5/abs(a1)))
     proj%G(itypat)%xp(1,ig)=a1
   end do
 end do

!debug
!write(*,*)'paw2wvl 178: erase me set gaussians real equal to hgh (for Li)'
!ABI_DEALLOCATE(proj%G(1)%xp)
!ABI_DEALLOCATE(proj%G(1)%psiat)
!ABI_ALLOCATE(proj%G(1)%psiat,(2,2))
!ABI_ALLOCATE(proj%G(1)%xp,(2,2))
!proj%G(1)%ndoc(:)=1 !1 gaussian per shell
!proj%G(1)%nexpo=2 ! two gaussians in total
!proj%G(1)%xp=zero
!proj%G(1)%xp(1,1)=0.666375d0
!proj%G(1)%xp(1,2)=1.079306d0
!proj%G(1)%psiat(1,:)=1.d0
!proj%G(1)%psiat(2,:)=zero

!
!begin debug
!
!write(*,*)'paw2wvl, comment me'
!and comment out variables 
!rmax=2.d0
!nr=rmax/(0.0001d0)
!ABI_ALLOCATE(r,(nr))
!ABI_ALLOCATE(f,(nr))
!ABI_ALLOCATE(y,(nr))
!step=rmax/real(nr-1,dp)
!do ir=1,nr
!r(ir)=real(ir-1,dp)*step
!end do
!!
!unitp=400
!do itypat=1,ntypat
!ig=0
!do i_shell=1,proj%G(itypat)%nshltot
!unitp=unitp+1
!f(:)=czero
!ng=proj%G(itypat)%ndoc(i_shell)
!ABI_ALLOCATE(g,(nr,ng))
!do ii=1,ng
!ig=ig+1
!fac=cmplx(proj%G(itypat)%psiat(1,ig),proj%G(itypat)%psiat(2,ig))
!a1=-0.5d0/(proj%G(itypat)%xp(1,ig)**2)
!arg=cmplx(a1,proj%G(itypat)%xp(2,ig))
!g(:,ii)=fac*exp(arg*r(:)**2)
!f(:)=f(:)+g(:,ii)
!end do
!do ir=1,nr
!write(unitp,'(9999f16.7)')r(ir),real(f(ir)),(real(g(ir,ii)),&
!&    ii=1,ng)
!end do
!ABI_DEALLOCATE(g)
!end do
!end do
!ABI_DEALLOCATE(r)
!ABI_DEALLOCATE(f)
!ABI_DEALLOCATE(y)
!stop 
!end debug


!jb=0
 do itypat=1,ntypat
   jb=0
   ll_ = 0
   nn_ = 0
   do ib=1,pawtab(itypat)%lmn_size 
     ll=pawtab(itypat)%indlmn(1,ib)
     nn=pawtab(itypat)%indlmn(3,ib)
!    write(*,*)ll,pawtab(itypat)%indlmn(2,ib),nn
     if(ib>1 .and. ll == ll_ .and. nn==nn_) cycle
     jb=jb+1
!    proj%G%nam(jb)= pawtab(itypat)%indlmn(1,ib)+1!l quantum number
!    write(*,*)jb,proj%G%nshltot,proj%G%nam(jb)
     proj%G(itypat)%nam(jb)= pawtab(itypat)%indlmn(1,ib)+1 !l quantum number
!    1 is added due to BigDFT convention
!    write(*,*)jb,pawtab(itypat)%indlmn(1:3,ib)
     ll_=ll
     nn_=nn
   end do
 end do
!Nullify remaining objects
 do itypat=1,ntypat
   nullify(proj%G(itypat)%rxyz)
   nullify(proj%G(itypat)%nshell)
 end do
!nullify(proj%G%rxyz)

!now the index l,m,n objects:
 lmnmax=maxval(pawtab(:)%lmn_size)
 ABI_ALLOCATE(wvl%paw%indlmn,(6,lmnmax,ntypat))
 wvl%paw%lmnmax=lmnmax
 wvl%paw%ntypes=ntypat
 wvl%paw%indlmn=0
 do itypat=1,ntypat
   lmnsz=pawtab(itypat)%lmn_size
   wvl%paw%indlmn(1:6,1:lmnsz,itypat)=pawtab(itypat)%indlmn(1:6,1:lmnsz)
 end do
 
!allocate and copy sij
!max_lmn2_size=max(pawtab(:)%lmn2_size,1)
 max_lmn2_size=lmnmax*(lmnmax+1)/2
 ABI_ALLOCATE(wvl%paw%sij,(max_lmn2_size,ntypat))
!sij is not yet calculated here.
!We copy this in gstate after call to pawinit
!do itypat=1,ntypat
!wvl%paw%sij(1:pawtab(itypat)%lmn2_size,itypat)=pawtab(itypat)%sij(:)
!end do

 ABI_ALLOCATE(wvl%paw%rpaw,(ntypat))
 do itypat=1,ntypat
   wvl%paw%rpaw(itypat)= pawtab(itypat)%rpaw
 end do

 if(allocated(wvl%npspcode_paw_init_guess)) then
   ABI_DEALLOCATE(wvl%npspcode_paw_init_guess)
 end if
 ABI_ALLOCATE(wvl%npspcode_paw_init_guess,(ntypat))
 do itypat=1,ntypat
   wvl%npspcode_paw_init_guess(itypat)=pawtab(itypat)%wvl%npspcode_init_guess
 end do

#else
 if (.false.) write(std_out) pawtab(1)%mesh_size,wvl%h(1),proj%nlpsp
#endif

!DEBUG
!write (std_out,*) ' paw2wvl : exit'
!stop
!ENDDEBUG

end subroutine paw2wvl
!!***
