!{\src2tex{textfont=tt}}
!!****f* ABINIT/mksupercell 
!! NAME
!!  mksupercell 
!!
!! FUNCTION 
!!  computes atomic positons, magnetic ordering of supercell
!!  
!! INPUTS
!!  magv_org (optional) magnetic ordering of atoms in primitive cell, 
!!   ordering of atoms given als 1 and -1, if not given fm is assumed     
!!  xred_org relative position of atoms in primitive cell
!!  rprimd_org unit cell dimensions of primitive cell
!!  natom=number of atoms in unit cell
!!  option= 1 output ion-ion distances / 2 output ordering of ion-ion distances / 3 output variables in varlist 
!!           according to ion-ion distances * magnetic ordering
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2017 ABINIT group (DJA)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! OUTPUT
!!  magv_sc magnetic ordering of atoms in supercell
!!  xred_sc relative position of atoms in supercell 
!!  rprimd_sc unit cell dimensions of supercell
!!  
!! SIDE EFFECTS
!!
!! PARENTS
!!      pawuj_det
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mksupercell(xred_org,magv_org,rprimd_org,nat_org,nat_sc,xred_sc,magv_sc,rprimd_sc,ext,prtvol) 

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mksupercell'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)              :: nat_org,nat_sc
 integer,intent(in),optional     :: prtvol
!arrays
 real(dp),intent(in)             :: rprimd_org(3,3)
 integer,intent(in)              :: ext(3)
 real(dp),intent(in)             :: xred_org(3,nat_org)
 real(dp),intent(out)            :: xred_sc(3,nat_sc) 
 real(dp),intent(out)            :: magv_sc(nat_sc)
 real(dp),intent(out)            :: rprimd_sc(3,3)
 integer,intent(in),optional     :: magv_org(nat_org)


!Local variables-------------------------------
!scalars
 integer                      :: prtvoll,ix,iy,iz,nprcl,iprcl,jdim,iatom
!arrays
 real(dp)                     :: magvv_org(nat_org)
 real(dp),allocatable         :: transv(:,:,:)

! *************************************************************************

!DEBUG
!write(std_out,*)'mksupercell: enter'
!END DEBUG

 if (present(magv_org)) then
   magvv_org=magv_org
 else
   magvv_org=(/ (1, iatom=1,nat_org)  /) 
 end if

 if (present(prtvol)) then
   prtvoll=prtvol
 else
   prtvoll=1
 end if

 rprimd_sc=reshape((/ (rprimd_org(ix,:)*ext(ix) ,ix=1,3) /),(/3,3 /))
 nprcl=product(ext)
 ABI_ALLOCATE(transv,(3,nat_org,nprcl))

 transv=reshape((/ (((((/ ix,iy,iz /),iatom=1,nat_org),ix=0,ext(1)-1),iy=0,ext(2)-1),iz=0,ext(3)-1) /), (/ 3, nat_org,nprcl/) )

!DEBUG
!write(std_out,*)'mksupercell: xred_org ' ,xred_org
!END DEBUG

 do iprcl=1,nprcl
   xred_sc(:,1+(iprcl-1)*nat_org:iprcl*nat_org)=xred_org+transv(:,:,iprcl)
   magv_sc(1+(iprcl-1)*nat_org:iprcl*nat_org)=magv_org
 end do


 do jdim=1,3
   xred_sc(jdim,:)=xred_sc(jdim,:)/ext(jdim)
 end do

!DEBUG
!write(std_out,*)'mksupercell: xred_sc ', xred_sc
!write(std_out,*)'mksupercell: magv_sc ', magv_sc
!END DEBUG

 ABI_DEALLOCATE(transv)
 
!DEBUG
!write(std_out,*)'mksupercell: leave'
!END DEBUG
end subroutine mksupercell 
!!***
