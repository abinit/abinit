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
 use m_errors
 use libxc_functionals

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

  real(dp) :: nelect
    ! Number of electrons in the cell (for Fermi-Amaldi only)

  real(dp) :: tphysel
    ! Physical temperature (for temperature-dependent functional)

  real(dp) :: xc_denpos
    ! density positivity value

  real(dp) :: xc_tb09_c
    ! Parameter for Tran-Blaha functional

 end type xcdata_type
!----------------------------------------------------------------------

 public :: xcdata_init                ! Initialize the object.
 public :: get_xclevel                ! Get the xclevel from ixc (as well as usefock)

contains
!!***

!!****f* m_xcdata/xcdata_init
!! NAME
!!  xcdata_init
!!
!! FUNCTION
!!  Init the structure. Mostly copy input variables, except compute xclevel.
!!
!! INPUTS
!!  intxc = 1 if the XC functional has to be interpolated on a more refined mesh than the FFT one
!!  ixc= index of exchange-correlation functional
!!  nelect = Number of electrons in the cell (for Fermi-Amaldi only)
!!  tphysel = Physical temperature (for temperature-dependent functional)
!!  usekden = 1 if the XC functional depends on the kinetic energy density
!!  vdw_xc = Choice of van-der-Waals density functional
!!  xc_tb09_c = Parameter for Tran-Blaha functional
!!
!! OUTPUT
!!  xcdata <type(xcdata_type)>= the data to calculate exchange-correlation are initialized
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xcdata_init(intxc,ixc,nelect,tphysel,usekden,vdw_xc,xc_tb09_c,xc_denpos,xcdata)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xcdata_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: intxc,ixc,usekden,vdw_xc
 real(dp),intent(in) :: nelect,tphysel,xc_denpos,xc_tb09_c
 type(xcdata_type), intent(out) :: xcdata
!Local variables-------------------------------
 integer :: xclevel

! *************************************************************************

 xcdata%intxc=intxc
 xcdata%ixc=ixc
 xcdata%usekden=usekden
 xcdata%vdw_xc=vdw_xc

 xcdata%nelect=nelect
 xcdata%tphysel=tphysel
 xcdata%xc_denpos=xc_denpos
 xcdata%xc_tb09_c=xc_tb09_c

!Compute xclevel 
 call get_xclevel(ixc,xclevel)
 xcdata%xclevel=xclevel

end subroutine xcdata_init
!!***


!!****f* m_xcdata/get_xclevel
!! NAME
!!  get_xclevel
!!
!! FUNCTION
!!  Compute xclevel.
!!
!! INPUTS
!!  ixc= index of exchange-correlation functional
!!
!! OUTPUT
!!  [usefock = 1 if the XC functional needs the Fock operator]
!!  xclevel= 0 if no XC functional except possibly Fock; 1 if LDA; 2 if GGA ; 3 for TDDFT kernel tests
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_xclevel(ixc,xclevel,usefock)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_xclevel'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ixc
 integer, intent(out) :: xclevel
 integer, intent(out), optional :: usefock

!Local variables-------------------------------
 integer :: ii,isiz,jj
 character(len=500) :: message

! *************************************************************************

 xclevel=0
 if( ( 1<=ixc .and. ixc<=10).or.(30<=ixc .and. ixc<=39).or.(ixc==50) )xclevel=1 ! LDA
 if( (11<=ixc .and. ixc<=19).or.(23<=ixc .and. ixc<=29) )xclevel=2 ! GGA
 if( 20<=ixc .and. ixc<=22 )xclevel=3 ! ixc for TDDFT kernel tests
 if( ixc>=40 .and. ixc<=42 )usefock=1 ! Hartree-Fock or internal hybrid functionals
 if( ixc>=41 .and. ixc<=42)xclevel=2 ! ixc for internal hybrids using GGA
 if (ixc<0) then                                  ! libXC: metaGGA and hybrid functionals
   xclevel=1
   do isiz=1,2
     if (isiz==1) ii=-ixc/1000
     if (isiz==2) ii=-ixc-ii*1000
     jj=libxc_functionals_family_from_id(ii)
     if (jj==XC_FAMILY_GGA    .or.jj==XC_FAMILY_MGGA) xclevel=2
     if (jj==XC_FAMILY_HYB_GGA.or.jj==XC_FAMILY_HYB_MGGA) then
       xclevel=2 ; usefock=1
       if (.not.libxc_functionals_gga_from_hybrid(hybrid_id=ii)) then
         write(message, '(a,i8,3a)' )&
&         'ixc=',ixc,' (libXC hybrid functional) is presently not allowed.',ch10,&
&         'Action: try another hybrid functional or use PAW.'
         MSG_ERROR(message)
       end if
     end if
   end do
 end if

end subroutine get_xclevel

end module m_xcdata
!!***
