!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawgylmg
!! NAME
!! pawgylmg
!!
!! FUNCTION
!! PAW: Compute Fourier transform of each g_l(r).Y_lm(r) function
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  kg(3,npw)=integer coordinates of planewaves in basis sphere for this k point.
!!  kpg(npw,nkpg)= (k+G) components (only if useylm=1)
!!  kpt(3)=reduced coordinates of k point
!!  lmax=1+max. value of l angular momentum
!!  nkpg=second dimension of kpg_k (0 if useylm=0)
!!  npw=number of planewaves in basis sphere
!!  ntypat=number of types of atoms
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ylm(npw,lmax**2)=real spherical harmonics for each G and k point
!!
!! OUTPUT
!!  gylmg(npw,lmax**2,ntypat)=Fourier transform of each g_l(r).Y_lm(r) function
!!
!! PARENTS
!!      suscep_stat
!!
!! CHILDREN
!!      splfit
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawgylmg(gprimd,gylmg,kg,kpg,kpt,lmax,nkpg,npw,ntypat,pawtab,ylm)

 use m_profiling_abi

 use defs_basis
 use m_errors
 use m_splines
 use m_pawtab, only : pawtab_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawgylmg'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lmax,nkpg,npw,ntypat
!arrays
 integer,intent(in) :: kg(3,npw)
 real(dp),intent(in) :: gprimd(3,3),kpg(npw,nkpg),kpt(3)
 real(dp),intent(in) :: ylm(npw,lmax**2)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

 real(dp),intent(out) :: gylmg(npw,lmax**2,ntypat)

!Local variables-------------------------------
!scalars
 integer :: ig,ilm,itypat,ll,l0,mm,mqgrid
 real(dp) :: kpg1,kpg2,kpg3,kpgc1,kpgc2,kpgc3
!arrays
 real(dp),allocatable :: glg(:),qgrid(:),kpgnorm(:),shpf(:,:),work(:)

! *************************************************************************

 DBG_ENTER("COLL")

!Get |k+G|:
 ABI_ALLOCATE(kpgnorm,(npw))
 if (nkpg<3) then
   do ig=1,npw
     kpg1=kpt(1)+dble(kg(1,ig));kpg2=kpt(2)+dble(kg(2,ig));kpg3=kpt(3)+dble(kg(3,ig))
     kpgc1=kpg1*gprimd(1,1)+kpg2*gprimd(1,2)+kpg3*gprimd(1,3)
     kpgc2=kpg1*gprimd(2,1)+kpg2*gprimd(2,2)+kpg3*gprimd(2,3)
     kpgc3=kpg1*gprimd(3,1)+kpg2*gprimd(3,2)+kpg3*gprimd(3,3)
     kpgnorm(ig)=sqrt(kpgc1*kpgc1+kpgc2*kpgc2+kpgc3*kpgc3)
   end do
 else
   do ig=1,npw
     kpgc1=kpg(ig,1)*gprimd(1,1)+kpg(ig,2)*gprimd(1,2)+kpg(ig,3)*gprimd(1,3)
     kpgc2=kpg(ig,1)*gprimd(2,1)+kpg(ig,2)*gprimd(2,2)+kpg(ig,3)*gprimd(2,3)
     kpgc3=kpg(ig,1)*gprimd(3,1)+kpg(ig,2)*gprimd(3,2)+kpg(ig,3)*gprimd(3,3)
     kpgnorm(ig)=sqrt(kpgc1*kpgc1+kpgc2*kpgc2+kpgc3*kpgc3)
   end do
 end if

 ABI_ALLOCATE(glg,(npw))
 ABI_ALLOCATE(work,(npw))

 write(std_out,*) ' lmax, pawtab(:)%lcut_size ', lmax, pawtab(:)%lcut_size

!Loop over types of atoms
 do itypat=1,ntypat

   mqgrid=pawtab(itypat)%mqgrid_shp
   ABI_ALLOCATE(qgrid,(mqgrid))
   ABI_ALLOCATE(shpf,(mqgrid,2))
   qgrid(1:mqgrid)=pawtab(itypat)%qgrid_shp(1:mqgrid)

!  Loops over (l,m) values
   do ll=0,pawtab(itypat)%lcut_size-1
     l0=ll**2+ll+1

     shpf(1:mqgrid,1:2)=pawtab(itypat)%shapefncg(1:mqgrid,1:2,1+ll)
     call splfit(qgrid,work,shpf,0,kpgnorm,glg,mqgrid,npw)

     do mm=-ll,ll
       ilm=l0+mm

       gylmg(1:npw,ilm,itypat)=ylm(1:npw,ilm)*glg(1:npw)

!      End loops over (l,m) values
     end do
   end do

!  End loop over atom types
   ABI_DEALLOCATE(qgrid)
   ABI_DEALLOCATE(shpf)
 end do

 ABI_DEALLOCATE(kpgnorm)
 ABI_DEALLOCATE(glg)
 ABI_DEALLOCATE(work)

 DBG_EXIT("COLL")

end subroutine pawgylmg
!!***
