!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_xcdata
!! NAME
!!  m_xcdata
!!
!! FUNCTION
!!  This module provides the definition of
!!  the xcdata_type used to drive the computation of the XC energy, potential, kernel, etc.
!!
!! COPYRIGHT
!!  Copyright (C) 2017 ABINIT group (XG)
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

module m_xcdata

 use defs_basis

 implicit none

 private
!!***

!!****t* m_xcdata/xcdata_type
!! NAME
!!  xcdata_type
!!
!! FUNCTION
!!   This object stores the input variables (and derived parameters) needed to compute the exchange-correlation functional, 
!!   not simply to define it. 
!!
!! NOTES
!!
!! SOURCE

 type, public :: xcdata_type

! Integer scalars

  integer :: intxc
    ! 1 if the XC functional has to be interpolated on a more refined mesh than the FFT one.
    ! 0 stick to the original FFT mesh

  integer :: ixc
    ! Choice of exchange-correlation functional. See input variable documentation

  integer :: usekden
    ! 1 if the XC functional depends on the kinetic energy density

  integer :: vdw_xc
    ! Choice of van-der-Waals density functional. See input variable documentation

  integer :: xclevel
    ! Determined from ixc
    ! 0 if no XC functional
    ! 1 if LDA-type XC functional
    ! 2 if GGA-type XC functional
    ! 3 if for TDDFT kernel

  real(dp) :: tphysel
    ! Physical temperature

  real(dp) :: xc_tb09_c
    ! Parameter for Tran-Blaha functional

 end type xcdata_type
!!***

!!****f* m_xc/xcdata_init
!! NAME
!!  xcdata_init
!!
!! FUNCTION
!!  Init all scalars, arrays and pointers in the structure.
!!
!! INPUTS
!!  intxc = 1 if the XC functional has to be interpolated on a more refined mesh than the FFT one
!!  ixc= index of exchange-correlation functional
!!  usekden = 1 if the XC functional depends on the kinetic energy density
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  xcdata <type(xcdata_type)>= the data to calculate exchange-correlation are initialized
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xcdata_init(intxc,ixc,usekden,tphysel,vdw_xc,xc_tb09_c,xcdata)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xcdata_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: intxc,ixc,usekden,vdw_xc
 real(dp),intent(in) :: tphysel,xc_tb09_c
!Local variables-------------------------------
 integer :: xclevel

! *************************************************************************

 DBG_ENTER("COLL")
 xcdata%intxc=intxc
 xcdata%ixc=ixc
 xcdata%usekden=usekden
 xcdata%vdwxc=vdw_xc

 xcdata%tphysel=tphysel
 xcdata%xc_tb09_c=xc_tb09_c

!Compute xclevel (XG 20170925 : what happens for libxc functionals ?)
 xclevel=0
 if ( ( 1<=ixc .and. ixc<=10).or.(30<=ixc .and. ixc<=39) ) xclevel=1 ! LDA
 if ( (11<=ixc .and. ixc<=19).or.(23<=ixc .and. ixc<=29) ) xclevel=2 ! GGA
 if ( 20<=ixc .and. ixc<=22 ) xclevel=3 ! ixc for TDDFT kernel tests
 xcdata%xclevel=xclevel

 DBG_EXIT("COLL")

end subroutine xcdata_init
!!***
