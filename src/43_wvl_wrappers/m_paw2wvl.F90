!!****m* ABINIT/m_paw2wvl
!! NAME
!!  paw2wvl
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2011-2020 ABINIT group (T. Rangel, MT)
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

module m_paw2wvl

 use defs_basis
 use defs_wvltypes
 use m_abicore
 use m_errors

 use m_pawtab, only : pawtab_type
 use m_paw_ij, only : paw_ij_type

 implicit none

 private
!!***

 public :: paw2wvl
 public :: paw2wvl_ij
 public :: wvl_paw_free
 public :: wvl_cprjreorder
!!***

contains
!!***

!!****f* ABINIT/paw2wvl
!! NAME
!!  paw2wvl
!!
!! FUNCTION
!!  Points WVL objects to PAW objects
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

subroutine paw2wvl(pawtab,proj,wvl)

 implicit none

!Arguments ------------------------------------
 type(pawtab_type),intent(in)::pawtab(:)
 type(wvl_internal_type), intent(inout) :: wvl
 type(wvl_projectors_type),intent(inout)::proj

!Local variables-------------------------------
#if defined HAVE_BIGDFT
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

#if defined HAVE_BIGDFT

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

!!****f* ABINIT/paw2wvl_ij
!! NAME
!!  paw2wvl_ij
!!
!! FUNCTION
!!  FIXME: add description.
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
!!      scfcv
!!
!! CHILDREN
!!      nullify_paw_ij_objects
!!
!! SOURCE

subroutine paw2wvl_ij(option,paw_ij,wvl)

#if defined HAVE_BIGDFT
 use BigDFT_API, only : nullify_paw_ij_objects
#endif
 implicit none

!Arguments ------------------------------------
 integer,intent(in)::option
 type(wvl_internal_type), intent(inout)::wvl
 type(paw_ij_type),intent(in) :: paw_ij(:)
!Local variables-------------------------------
#if defined HAVE_BIGDFT
 integer :: iatom,iaux,my_natom
 character(len=500) :: message
#endif

! *************************************************************************

 DBG_ENTER("COLL")

#if defined HAVE_BIGDFT
 my_natom=size(paw_ij)

!Option==1: allocate and copy
 if(option==1) then
   ABI_DATATYPE_ALLOCATE(wvl%paw%paw_ij,(my_natom))
   do iatom=1,my_natom
     call nullify_paw_ij_objects(wvl%paw%paw_ij(iatom))
     wvl%paw%paw_ij(iatom)%cplex          =paw_ij(iatom)%qphase
     wvl%paw%paw_ij(iatom)%cplex_dij      =paw_ij(iatom)%cplex_dij
     wvl%paw%paw_ij(iatom)%has_dij        =paw_ij(iatom)%has_dij
     wvl%paw%paw_ij(iatom)%has_dijfr      =0
     wvl%paw%paw_ij(iatom)%has_dijhartree =0
     wvl%paw%paw_ij(iatom)%has_dijhat     =0
     wvl%paw%paw_ij(iatom)%has_dijso      =0
     wvl%paw%paw_ij(iatom)%has_dijU       =0
     wvl%paw%paw_ij(iatom)%has_dijxc      =0
     wvl%paw%paw_ij(iatom)%has_dijxc_val  =0
     wvl%paw%paw_ij(iatom)%has_exexch_pot =0
     wvl%paw%paw_ij(iatom)%has_pawu_occ   =0
     wvl%paw%paw_ij(iatom)%lmn_size       =paw_ij(iatom)%lmn_size
     wvl%paw%paw_ij(iatom)%lmn2_size      =paw_ij(iatom)%lmn2_size
     wvl%paw%paw_ij(iatom)%ndij           =paw_ij(iatom)%ndij
     wvl%paw%paw_ij(iatom)%nspden         =paw_ij(iatom)%nspden
     wvl%paw%paw_ij(iatom)%nsppol         =paw_ij(iatom)%nsppol
     if (paw_ij(iatom)%has_dij/=0) then
       iaux=paw_ij(iatom)%cplex_dij*paw_ij(iatom)%lmn2_size
       ABI_ALLOCATE(wvl%paw%paw_ij(iatom)%dij,(iaux,paw_ij(iatom)%ndij))
       wvl%paw%paw_ij(iatom)%dij(:,:)=paw_ij(iatom)%dij(:,:)
     end if
   end do

!  Option==2: deallocate
 elseif(option==2) then
   do iatom=1,my_natom
     wvl%paw%paw_ij(iatom)%has_dij=0
     if (associated(wvl%paw%paw_ij(iatom)%dij)) then
       ABI_DEALLOCATE(wvl%paw%paw_ij(iatom)%dij)
     end if
   end do
   ABI_DATATYPE_DEALLOCATE(wvl%paw%paw_ij)

!  Option==3: only copy
 elseif(option==3) then
   do iatom=1,my_natom
     wvl%paw%paw_ij(iatom)%cplex     =paw_ij(iatom)%qphase
     wvl%paw%paw_ij(iatom)%cplex_dij =paw_ij(iatom)%cplex_dij
     wvl%paw%paw_ij(iatom)%lmn_size  =paw_ij(iatom)%lmn_size
     wvl%paw%paw_ij(iatom)%lmn2_size =paw_ij(iatom)%lmn2_size
     wvl%paw%paw_ij(iatom)%ndij      =paw_ij(iatom)%ndij
     wvl%paw%paw_ij(iatom)%nspden    =paw_ij(iatom)%nspden
     wvl%paw%paw_ij(iatom)%nsppol    =paw_ij(iatom)%nsppol
     wvl%paw%paw_ij(iatom)%dij(:,:)  =paw_ij(iatom)%dij(:,:)
   end do

 else
   message = 'paw2wvl_ij: option should be equal to 1, 2 or 3'
   MSG_ERROR(message)
 end if

#else
 if (.false.) write(std_out,*) option,wvl%h(1),paw_ij(1)%ndij
#endif

 DBG_EXIT("COLL")

end subroutine paw2wvl_ij
!!***

!!****f* ABINIT/wvl_paw_free
!! NAME
!!  wvl_paw_free
!!
!! FUNCTION
!!  Frees memory for WVL+PAW implementation
!!
!! INPUTS
!!  ntypat = number of atom types
!!  wvl= wvl type
!!  wvl_proj= wvl projector type
!!
!! OUTPUT
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


subroutine wvl_paw_free(wvl)

 implicit none

!Arguments ------------------------------------
 type(wvl_internal_type),intent(inout) :: wvl

!Local variables-------------------------------

! *************************************************************************

#if defined HAVE_BIGDFT

!PAW objects
 if( associated(wvl%paw%spsi)) then
   ABI_DEALLOCATE(wvl%paw%spsi)
 end if
 if( associated(wvl%paw%indlmn)) then
   ABI_DEALLOCATE(wvl%paw%indlmn)
 end if
 if( associated(wvl%paw%sij)) then
   ABI_DEALLOCATE(wvl%paw%sij)
 end if
 if( associated(wvl%paw%rpaw)) then
   ABI_DEALLOCATE(wvl%paw%rpaw)
 end if

!rholoc
 if( associated(wvl%rholoc%msz )) then
   ABI_DEALLOCATE(wvl%rholoc%msz)
 end if
 if( associated(wvl%rholoc%d )) then
   ABI_DEALLOCATE(wvl%rholoc%d)
 end if
 if( associated(wvl%rholoc%rad)) then
   ABI_DEALLOCATE(wvl%rholoc%rad)
 end if
 if( associated(wvl%rholoc%radius)) then
   ABI_DEALLOCATE(wvl%rholoc%radius)
 end if

#else
 if (.false.) write(std_out,*) wvl%h(1)
#endif

!paw%paw_ij and paw%cprj are allocated and deallocated inside vtorho

end subroutine wvl_paw_free
!!***


!!****f* ABINIT/wvl_cprjreorder
!! NAME
!!  wvl_cprjreorder
!!
!! FUNCTION
!! Change the order of a wvl-cprj datastructure
!!   From unsorted cprj to atom-sorted cprj (atm_indx=atindx)
!!   From atom-sorted cprj to unsorted cprj (atm_indx=atindx1)
!!
!! INPUTS
!!  atm_indx(natom)=index table for atoms
!!   From unsorted wvl%paw%cprj to atom-sorted wvl%paw%cprj (atm_indx=atindx)
!!   From atom-sorted wvl%paw%cprj to unsorted wvl%paw%cprj (atm_indx=atindx1)
!!
!! OUTPUT
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      cprj_clean,cprj_paw_alloc
!!
!! SOURCE

subroutine wvl_cprjreorder(wvl,atm_indx)

#if defined HAVE_BIGDFT
 use BigDFT_API,only : cprj_objects,cprj_paw_alloc,cprj_clean
 use dynamic_memory
#endif
 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 integer,intent(in) :: atm_indx(:)
 type(wvl_internal_type),intent(inout),target :: wvl

!Local variables-------------------------------
#if defined HAVE_BIGDFT
!scalars
 integer :: iexit,ii,jj,kk,n1atindx,n1cprj,n2cprj,ncpgr
 character(len=100) :: msg
!arrays
 integer,allocatable :: nlmn(:)
 type(cprj_objects),pointer :: cprj(:,:)
 type(cprj_objects),allocatable :: cprj_tmp(:,:)
#endif

! *************************************************************************

 DBG_ENTER("COLL")

#if defined HAVE_BIGDFT
 cprj => wvl%paw%cprj

 n1cprj=size(cprj,dim=1);n2cprj=size(cprj,dim=2)
 n1atindx=size(atm_indx,dim=1)
 if (n1cprj==0.or.n2cprj==0.or.n1atindx<=1) return
 if (n1cprj/=n1atindx) then
   msg='wrong sizes!'
   MSG_BUG(msg)
 end if

!Nothing to do when the atoms are already sorted
 iexit=1;ii=0
 do while (iexit==1.and.ii<n1atindx)
   ii=ii+1
   if (atm_indx(ii)/=ii) iexit=0
 end do
 if (iexit==1) return

 ABI_ALLOCATE(nlmn,(n1cprj))
 do ii=1,n1cprj
   nlmn(ii)=cprj(ii,1)%nlmn
 end do
 ncpgr=cprj(1,1)%ncpgr

 ABI_DATATYPE_ALLOCATE(cprj_tmp,(n1cprj,n2cprj))
 call cprj_paw_alloc(cprj_tmp,ncpgr,nlmn)
 do jj=1,n2cprj
   do ii=1,n1cprj
     cprj_tmp(ii,jj)%nlmn=nlmn(ii)
     cprj_tmp(ii,jj)%ncpgr=ncpgr
     cprj_tmp(ii,jj)%cp(:,:)=cprj(ii,jj)%cp(:,:)
     if (ncpgr>0) cprj_tmp(ii,jj)%dcp(:,:,:)=cprj(ii,jj)%dcp(:,:,:)
   end do
 end do

 call cprj_clean(cprj)

 do jj=1,n2cprj
   do ii=1,n1cprj
     kk=atm_indx(ii)
     cprj(kk,jj)%nlmn=nlmn(ii)
     cprj(kk,jj)%ncpgr=ncpgr
     cprj(kk,jj)%cp=f_malloc_ptr((/2,nlmn(ii)/),id='cprj%cp')
     cprj(kk,jj)%cp(:,:)=cprj_tmp(ii,jj)%cp(:,:)
     if (ncpgr>0) then
       cprj(kk,jj)%dcp=f_malloc_ptr((/2,ncpgr,nlmn(ii)/),id='cprj%dcp')
       cprj(kk,jj)%dcp(:,:,:)=cprj_tmp(kk,jj)%dcp(:,:,:)
     end if
   end do
 end do

 call cprj_clean(cprj_tmp)
 ABI_DATATYPE_DEALLOCATE(cprj_tmp)
 ABI_DEALLOCATE(nlmn)

#else
 if (.false.) write(std_out,*) atm_indx(1),wvl%h(1)
#endif

 DBG_EXIT("COLL")

end subroutine wvl_cprjreorder
!!***

end module m_paw2wvl
!!***
