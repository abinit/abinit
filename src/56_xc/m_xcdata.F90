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
!!  Copyright (C) 2017-2019 ABINIT group (XG)
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
 use m_dtset, only : dataset_type
 use m_drivexc, only : size_dvxc

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

  integer :: auxc_ixc
    ! Choice of auxiliary exchange-correlation functional. See input variable documentation
    ! If 0, there is no auxiliary xc functional, the one corresponding to ixc has to be used.

  integer :: intxc
    ! 1 if the XC functional has to be interpolated on a more refined mesh than the FFT one.
    ! 0 stick to the original FFT mesh

  integer :: ixc
    ! Choice of exchange-correlation functional. See input variable documentation

  integer :: nspden
    ! Number of spin components of the density

  integer :: usefock
    ! 1 if the XC functional includes a (possibly screened) Fock contribution

  integer :: usegradient
    ! 1 if the XC functional depends on the density gradient

  integer :: uselaplacian
    ! 1 if the XC functional depends on the density laplacian

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

! Real scalars

  real(dp) :: hyb_mixing
    ! Parameter for mixing Fock exchange in native PBEx functionals

  real(dp) :: nelect
    ! Number of electrons in the cell (for Fermi-Amaldi only)

  real(dp) :: tphysel
    ! Physical temperature (for temperature-dependent functional)

  real(dp) :: xc_denpos
    ! density positivity value

 end type xcdata_type

!----------------------------------------------------------------------

 public :: xcdata_init                ! Initialize the object.
 public :: get_xclevel                ! Get the xclevel from ixc (as well as usefock)
 public :: get_auxc_ixc               ! Get the auxiliary xc functional (if it exists)

contains
!!***

!!****f* m_xcdata/xcdata_init
!! NAME
!!  xcdata_init
!!
!! FUNCTION
!!  Init the structure. Mostly copy input variables, except compute and usefock and xclevel.
!!
!! INPUTS
!!  [dtset = the dataset from which the other input variables are taken, if they are not present]
!!  [auxc_ixc = possibly the index of the auxiliary xc functional, otherwise 0.]
!!  [hyb_mixing = parameter for mixing Fock exchange in native PBEx functionals]
!!  [intxc = 1 if the XC functional has to be interpolated on a more refined mesh than the FFT one]
!!  [ixc= index of exchange-correlation functional]
!!  [nelect = Number of electrons in the cell (for Fermi-Amaldi only)]
!!  [tphysel = Physical temperature (for temperature-dependent functional)]
!!  [vdw_xc = Choice of van-der-Waals density functional]
!!
!! OUTPUT
!!  xcdata <type(xcdata_type)>= the data to calculate exchange-correlation are initialized
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      calc_vhxc_me,energy,m_kxc,nonlinear,nres2vres,odamix,prcref,prcref_PMA
!!      respfn,rhotov,scfcv,setvtr,xchybrid_ncpp_cc
!!
!! CHILDREN
!!      get_xclevel
!!
!! SOURCE

subroutine xcdata_init(xcdata,auxc_ixc,dtset,hyb_mixing,intxc,ixc,nelect,nspden,tphysel,&
&                      vdw_xc,xc_denpos)

!Arguments ------------------------------------
!scalars
 integer, intent(in),optional :: auxc_ixc,intxc,ixc,nspden,vdw_xc
 real(dp),intent(in),optional :: hyb_mixing,nelect,tphysel,xc_denpos
 type(dataset_type), intent(in),optional :: dtset
 type(xcdata_type), intent(out) :: xcdata
!Local variables-------------------------------
 integer :: nspden_updn
 character(len=500) :: msg

! *************************************************************************

 if(present(dtset))then
   xcdata%auxc_ixc=dtset%auxc_ixc
   xcdata%intxc=dtset%intxc
   xcdata%ixc=dtset%ixc
   xcdata%nspden=dtset%nspden
   xcdata%vdw_xc=dtset%vdw_xc

   xcdata%hyb_mixing=abs(dtset%hyb_mixing) ! Warning : the absolute value is needed, because of the singular way
                                           ! to define the default values for this input variable.
   xcdata%nelect=dtset%nelect
   xcdata%tphysel=dtset%tphysel
   xcdata%xc_denpos=dtset%xc_denpos

 else
   if(.not.(present(auxc_ixc).and.present(intxc).and.present(ixc).and.&
&           present(vdw_xc).and.present(hyb_mixing).and.&
&           present(nelect).and.present(nspden).and.&
&           present(tphysel).and.present(xc_denpos)))then
     msg='If dtset is not provided, all the other optional arguments must be provided, which is not the case!'
     MSG_BUG(msg)
   endif
 endif

 if(present(auxc_ixc))  xcdata%auxc_ixc=auxc_ixc
 if(present(intxc))     xcdata%intxc=intxc
 if(present(ixc))       xcdata%ixc=ixc
 if(present(nspden))    xcdata%nspden=nspden
 if(present(vdw_xc))    xcdata%vdw_xc=vdw_xc

 if(present(hyb_mixing))xcdata%hyb_mixing=hyb_mixing
 if(present(nelect))    xcdata%nelect=nelect
 if(present(tphysel))   xcdata%tphysel=tphysel
 if(present(xc_denpos)) xcdata%xc_denpos=xc_denpos

!Compute xclevel
 call get_xclevel(xcdata%ixc,xcdata%xclevel,usefock=xcdata%usefock)

!Compute usegradient,uselaplacian,usekden
 nspden_updn=min(xcdata%nspden,2)
 call size_dvxc(xcdata%ixc,1,nspden_updn,usegradient=xcdata%usegradient,&
&               uselaplacian=xcdata%uselaplacian,usekden=xcdata%usekden)

end subroutine xcdata_init
!!***

!----------------------------------------------------------------------

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
!!      invars2,m_xcdata,setup_sigma
!!
!! CHILDREN
!!      get_xclevel
!!
!! SOURCE

subroutine get_xclevel(ixc,xclevel,usefock)

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
 if( (11<=ixc .and. ixc<=19).or.(23<=ixc .and. ixc<=29).or. ixc==1402000)xclevel=2 ! GGA
 if( 20<=ixc .and. ixc<=22 )xclevel=3 ! ixc for TDDFT kernel tests
 if(present(usefock))then
   usefock=0
   if( ixc>=40 .and. ixc<=42 )usefock=1 ! Hartree-Fock or internal hybrid functionals
 endif
 if( ixc>=41 .and. ixc<=42)xclevel=2 ! ixc for internal hybrids using GGA
 if (ixc<0) then                                  ! libXC: metaGGA and hybrid functionals
   xclevel=1
   do isiz=1,2
!    ixc has ABINIT sign convention
!    ii has Libxc sign convention
     if (isiz==1) ii=-ixc/1000
     if (isiz==2) ii=-ixc-ii*1000
     jj=libxc_functionals_family_from_id(ii)
     if (jj==XC_FAMILY_GGA    .or.jj==XC_FAMILY_MGGA) xclevel=2
     if (jj==XC_FAMILY_HYB_GGA.or.jj==XC_FAMILY_HYB_MGGA) then
       xclevel=2
       if(present(usefock))then
         usefock=1
       endif
       if (.not.libxc_functionals_gga_from_hybrid(hybrid_id=ii)) then
         write(message, '(a,i8,3a,i8,2a,2i8,2a)' )&
&         'ixc=',ixc,' (libXC hybrid functional) is presently not allowed.',ch10,&
&         'XC_FAMILY_HYB_GGA=',XC_FAMILY_HYB_GGA,ch10,&
&         'ii,jj=',ii,jj,ch10,&
&         'Action: try another hybrid functional.'
         MSG_ERROR(message)
       end if
     end if
   end do
 end if

end subroutine get_xclevel
!!***

!----------------------------------------------------------------------

!!****f* m_xcdata/get_auxc_ixc
!! NAME
!!  get_auxc_ixc
!!
!! FUNCTION
!!  Returns the ixc of an auxiliary XC functional to be used instead of the input ixc
!!  For most of the functionals, there is no need of an auxiliary functional, in which case auxc_ixc=0
!!  For hybrid functionals, on the contrary, some speedup can be achieved by using such an auxiliary functional
!!  Note that this XC functional intend to replace the whole ixc functional. Generally speakin, it should be
!!  mistaken for the GGA part of the hybrid functional (that is for exchange only, actually).
!!
!!  At present, always return ixc=1, but this might change in the future ...
!!
!! INPUTS
!!  ixc= index of exchange-correlation functional
!!
!! OUTPUT
!!  auxc_ixc= 0 if no need of an auxiliary functional, otherwise, returns the ixc of an auxiliary functional.
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      calc_vhxc_me,invars2
!!
!! CHILDREN
!!      get_xclevel
!!
!! SOURCE

subroutine get_auxc_ixc(auxc_ixc,ixc)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ixc
 integer, intent(out) :: auxc_ixc

!Local variables-------------------------------
 integer :: usefock,xclevel
!integer :: gga_id(2)

! *************************************************************************

 auxc_ixc=11
!Native hybrid functionals from ABINIT
 if (ixc==40.or.ixc==41.or.ixc==42) then
   auxc_ixc = 11
!Hybrid functionals from libxc
 else if (ixc<0) then
   call get_xclevel(ixc,xclevel,usefock)
   if(usefock==1)then
     auxc_ixc=11
!    if (libxc_functionals_gga_from_hybrid(hybrid_id=ixc,gga_id=gga_id)) then
!      auxc_ixc=-gga_id(1)*1000-gga_id(2)
!    endif
   end if
 end if

end subroutine get_auxc_ixc

end module m_xcdata
!!***
