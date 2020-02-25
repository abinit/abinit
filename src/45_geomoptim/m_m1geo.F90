!!****m* ABINIT/m_m1geo
!! NAME
!! m_m1geo
!!
!! FUNCTION
!! This module contains definition the m1geo type (="move unique geometry")
!! and its related init and destroy routines
!!
!! COPYRIGHT
!! Copyright (C) 2018-2020 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_m1geo

 use defs_basis
 use m_profiling_abi
 use m_abimover
 use m_abihist
 use m_dtset
 use m_dtfil

 implicit none

 private

 public :: m1geo_init
 public :: m1geo_destroy

!!***

!----------------------------------------------------------------------

!!****t* m_m1geo/m1geo
!! NAME
!! m1geo type
!!
!! FUNCTION
!! This datatype has the purpose to store all the data taken
!! usually from dtset (but not only) needed to move a single geometry
!! opbtained by the weighting of different images,
!! using the algorithms defined by ionmov.
!!
!! SOURCE

type, public :: m1geo_type
! Scalars
  integer  :: dt_chkdilatmx  ! chkdilatmx from dtset
  integer  :: hmctt          ! hmctt from dtset
  integer  :: icycle         ! cycle index
  integer  :: iexit          ! if nonzero, will have to exit
  integer  :: itime          ! time index for the move1geo algorithms
  integer  :: nctime         !
  integer  :: ncycle         ! foreseen number of cycles
  integer  :: nerr_dilatmx   ! number of dilatmx trespass errors, if chkdilatmx
  integer  :: nimage         ! nimage from dtset
  integer  :: npsp           ! npsp from dtset
  integer  :: ntime          ! number of time steps for the move1geo algorithms
  integer  :: usewvl         ! usewvl from dtset
  real(dp) :: dilatmx        ! dilatmx from dtset
  logical  :: skipcycle      ! .TRUE. when the remaining of the cycle has to be skipped. .FALSE. otherwise.
! Arrays
  real(dp) :: rprimd_orig(3,3)   ! rprimd from dtset
  real(dp),allocatable :: mixesimgf(:)   ! mixesimgf from dtset
! Character string
  character(len=fnlen) :: filnam_ds4
! Datatypes
  type(abimover) :: ab_mover
  type(abimover_specs) :: specs
  type(ab_xfh_type) :: ab_xfh_1geo
  type(delocint) :: deloc
  type(abihist) :: hist_1geo
  type(mttk_type) :: mttk_vars

 end type m1geo_type

contains  !=============================================================
!!***

!!****f* m_gstateimg/m1geo_init
!! NAME
!! m1geo_init
!!
!! FUNCTION
!! Initializes the m1geo structure
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

 subroutine m1geo_init(dtfil,dtset,m1geo_param)

!Arguments ------------------------------------
!scalars
 type(datafiles_type),target,intent(in) :: dtfil
 type(dataset_type),target,intent(in) :: dtset
 type(m1geo_type),intent(inout) :: m1geo_param

!Local variables
!integer
 integer :: ionmov,natom,ncycle,nhisttot,nimage,ntimimage,ntypat
 real(dp), allocatable :: amu_curr(:)

!************************************************************************


 ionmov=dtset%ionmov
 natom=dtset%natom
 nimage=dtset%nimage
 ntimimage=dtset%ntimimage
 ntypat=dtset%ntypat

!Straight initialization from dtset
 m1geo_param%dt_chkdilatmx  =dtset%chkdilatmx
 m1geo_param%dilatmx        =dtset%dilatmx
 m1geo_param%hmctt          =dtset%hmctt
 m1geo_param%nctime         =dtset%nctime
 m1geo_param%nimage         =dtset%nimage
 m1geo_param%npsp           =dtset%npsp
 m1geo_param%usewvl         =dtset%usewvl

 ABI_MALLOC(m1geo_param%mixesimgf,(nimage))
 m1geo_param%mixesimgf(:)=dtset%mixesimgf(1:nimage)


!From dtset, with change or adjustment of meaning
 m1geo_param%ntime          =dtset%ntimimage  ! Beware ntime vs ntimimage
 m1geo_param%rprimd_orig    =dtset%rprimd_orig(:,:,1)

 ABI_MALLOC(amu_curr,(ntypat))
 amu_curr(:)=dtset%amu_orig(:,1)

!From dtfil
 m1geo_param%filnam_ds4     =dtfil%filnam_ds(4)

!###########################################################
!Init sub-datastructures

!1) ab_mover
 call abimover_ini(m1geo_param%ab_mover,amu_curr,dtfil,dtset,m1geo_param%specs)
 ABI_DEALLOCATE(amu_curr)

!2) deloc
 if(ionmov==10 .or. ionmov==11)then
   call delocint_ini(m1geo_param%deloc)
 end if

!3) hist
!In a simple approach, hist_1geo is initialized here.
!Likely one should also make an interface with the hist(nimage) from the calling routine, and average it over the images...

 ncycle=m1geo_param%specs%ncycle
 if(ionmov==25.and.dtset%hmctt>=0) then
   ncycle=dtset%hmctt
   if(dtset%hmcsst>0.and.dtset%optcell/=0) then
     ncycle=ncycle+dtset%hmcsst
   endif
 endif
 nhisttot=ncycle*ntimimage  ! WARNING, here ntime is used instead of ntimimage
 if (dtset%nctime>0) nhisttot=nhisttot+1

!We just store the needed history step not all of them.
 if(m1geo_param%specs%nhist/=-1) nhisttot = m1geo_param%specs%nhist ! We don't need to store all the history

 call abihist_init(m1geo_param%hist_1geo,natom,ntimimage,m1geo_param%specs%isVused,m1geo_param%specs%isARused)

!4) ab_xfh
 m1geo_param%ab_xfh_1geo%mxfh=ntimimage
 m1geo_param%ab_xfh_1geo%nxfh=0
 m1geo_param%ab_xfh_1geo%nxfhr=0
 ABI_MALLOC(m1geo_param%ab_xfh_1geo%xfhist,(3,natom+4,2,ntimimage))

!5) mttk
 if (ionmov==13)then
   call mttk_ini(m1geo_param%mttk_vars,m1geo_param%ab_mover%nnos)
 end if

!###########################################################
!Init other values

 m1geo_param%ncycle=m1geo_param%specs%ncycle
 m1geo_param%icycle=0
 m1geo_param%iexit=0
 m1geo_param%itime=0
 m1geo_param%nerr_dilatmx=0
 m1geo_param%skipcycle=.FALSE.

 end subroutine m1geo_init
!!***

!----------------------------------------------------------------------

!!****f* m_gstateimg/m1geo_destroy
!! NAME
!! m1geo_destroy
!!
!! FUNCTION
!! Destroy the m1geo structure
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

 subroutine m1geo_destroy(m1geo_param)

!Arguments ------------------------------------
!scalars
 type(m1geo_type),intent(inout) :: m1geo_param

!************************************************************************

 ABI_SFREE(m1geo_param%mixesimgf)
 ABI_SFREE(m1geo_param%ab_xfh_1geo%xfhist)

 call abimover_destroy(m1geo_param%ab_mover)
 call delocint_fin(m1geo_param%deloc)
 call abihist_free(m1geo_param%hist_1geo)
 call mttk_fin(m1geo_param%mttk_vars)

 end subroutine m1geo_destroy

!----------------------------------------------------------------------

end module m_m1geo
!!***
