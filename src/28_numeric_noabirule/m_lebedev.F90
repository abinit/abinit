!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_lebedev
!! NAME
!!  m_lebedev
!!
!! FUNCTION
!!  This module contains routines for performing integrations
!!  on the sphere using lebedev-laikov angular grids.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2016 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_lebedev

 use defs_basis
 use m_profiling_abi
 use m_errors

 implicit none

 private 
!!***

!----------------------------------------------------------------------

!!****t* m_lebedev/lebedev_grid_t
!! NAME
!! lebedev_grid_t
!! 
!! FUNCTION
!! Structure storing the knots and the weights of the lebedev-laikov grid.
!! 
!! SOURCE

 type,public :: lebedev_grid_t
   integer :: npts=0
   real(dp),allocatable :: versor(:,:) 
   real(dp),allocatable :: weight(:)  
 end type lebedev_grid_t

!----------------------------------------------------------------------

 integer,private,parameter :: lebedev_ngrids=32
 ! The number of grids available.

 logical,save,private :: gridset_is_init=.FALSE.
 ! Flag defining whether the set of angular grids have been already computed and stored.

 ! The set of grids stored here so that we do not have to recompute them for each integration.
 type(lebedev_grid_t),save,public :: Lgridset(lebedev_ngrids)  

 ! The number of points in each grid.
 integer,public,parameter :: lebedev_npts(lebedev_ngrids)=(/ &
0006,&  
0014,&
0026,&
0038,&
0050,&
0074,&
0086,&
0110,&
0146,&
0170,&
0194,&
0230,&
0266,&
0302,&
0350,&
0434,&
0590,&
0770,&
0974,&
1202,&
1454,&
1730,&
2030,&
2354,&
2702,&
3074,&
3470,&
3890,&
4334,&
4802,&
5294,&
5810/)

 public :: init_lebedev_gridset       ! Calculate and save the 32 angular grids.
 public :: destroy_lebedev_gridset    ! Free the grids stored in this module.
 !public :: lebedev_quadrature         ! Integrate a given function on the sphere.
 public :: m_lebedev_is_init          ! Returns true if lebedev-laikov grids are already stored and computed.

 ! commented because it causes problems to the new version of abilint
 !interface lebedev_quadrature
 !  module procedure lebedev_quadrature_cplx
 !end interface lebedev_quadrature

CONTAINS  !===========================================================
!!***

!!****f* m_lebedev/init_lebedev_grid
!! NAME
!!  init_lebedev_grid
!!
!! FUNCTION
!!  Initialize a lebedev grid. 
!!
!! INPUTS
!!  seq_idx=Sequential index comprised between 1 and 32 defining the order of the mesh.
!!
!! OUTPUT
!!  Lgrid<lebedev_grid_t>=The grid fully initialized.
!!
!! PARENTS
!!      m_lebedev
!!
!! CHILDREN
!!      ld0006,ld0014,ld0026,ld0038,ld0050,ld0074,ld0086,ld0110,ld0146,ld0170
!!      ld0194,ld0230,ld0266,ld0302,ld0350,ld0434,ld0590,ld0770,ld0974,ld1202
!!      ld1454,ld1730,ld2030,ld2354,ld2702,ld3074,ld3470,ld3890,ld4334,ld4802
!!      ld5294,ld5810
!!
!! SOURCE

subroutine init_lebedev_grid(Lgrid,seq_idx)

 use defs_basis

!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_lebedev_grid'
!End of the abilint section

 integer,intent(in) :: seq_idx
 type(lebedev_grid_t),intent(out) :: Lgrid

!Local variables-------------------------------
 integer :: npts,ii
 real(dp),allocatable :: xx(:),yy(:),zz(:)
! *********************************************************************

 if (seq_idx<1.or.seq_idx>lebedev_ngrids) then
   MSG_ERROR("seq_idx out of range")
 end if

 npts = lebedev_npts(seq_idx)
 ABI_MALLOC(xx,(npts))
 ABI_MALLOC(yy,(npts))
 ABI_MALLOC(zz,(npts))
 ABI_MALLOC(Lgrid%weight,(npts))

 call build_lebedev_grid(seq_idx,xx,yy,zz,Lgrid%weight)

 Lgrid%npts = lebedev_npts(seq_idx)
 ABI_MALLOC(Lgrid%versor,(3,Lgrid%npts))

 do ii=1,Lgrid%npts
   Lgrid%versor(:,ii) = (/xx(ii),yy(ii),zz(ii)/)
 end do

 ABI_FREE(xx)
 ABI_FREE(yy)
 ABI_FREE(zz)

end subroutine init_lebedev_grid
!!***

!----------------------------------------------------------------------

!!****f* m_lebedev/destroy_lebedev_grid
!! NAME
!!  destroy_lebedev_grid
!!
!! FUNCTION
!!  Free an instance of lebedev_grid_t
!!
!! PARENTS
!!      m_lebedev
!!
!! CHILDREN
!!      ld0006,ld0014,ld0026,ld0038,ld0050,ld0074,ld0086,ld0110,ld0146,ld0170
!!      ld0194,ld0230,ld0266,ld0302,ld0350,ld0434,ld0590,ld0770,ld0974,ld1202
!!      ld1454,ld1730,ld2030,ld2354,ld2702,ld3074,ld3470,ld3890,ld4334,ld4802
!!      ld5294,ld5810
!!
!! SOURCE

subroutine destroy_lebedev_grid(Lgrid)

 use defs_basis

!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_lebedev_grid'
!End of the abilint section

 type(lebedev_grid_t),intent(inout) :: Lgrid

! *********************************************************************

 Lgrid%npts=0
 if (allocated(Lgrid%versor)) then
   ABI_FREE(Lgrid%versor)
 end if
 if (allocated(Lgrid%weight)) then
   ABI_FREE(Lgrid%weight)
 end if

end subroutine destroy_lebedev_grid
!!***

!----------------------------------------------------------------------

!!****f* m_lebedev/init_lebedev_gridset
!! NAME
!!  init_lebedev_gridset
!!
!! FUNCTION
!!  Initialize the 32 lebedev-laikov angular grids. 
!! 
!! SIDE EFFECTS 
!!  Lgridset(1:32) are initialized and saved in this module.
!!
!! PARENTS
!!      m_screening
!!
!! CHILDREN
!!      ld0006,ld0014,ld0026,ld0038,ld0050,ld0074,ld0086,ld0110,ld0146,ld0170
!!      ld0194,ld0230,ld0266,ld0302,ld0350,ld0434,ld0590,ld0770,ld0974,ld1202
!!      ld1454,ld1730,ld2030,ld2354,ld2702,ld3074,ld3470,ld3890,ld4334,ld4802
!!      ld5294,ld5810
!!
!! SOURCE

subroutine init_lebedev_gridset()

 use defs_basis

!Local variables-------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_lebedev_gridset'
!End of the abilint section

 integer :: igr
! *********************************************************************

 if (.not.gridset_is_init) then
   do igr=1,lebedev_ngrids
     call init_lebedev_grid(Lgridset(igr),igr)
   end do
   gridset_is_init=.TRUE.
 end if

end subroutine init_lebedev_gridset
!!***

!----------------------------------------------------------------------

!!****f* m_lebedev/destroy_lebedev_gridset
!! NAME
!!  destroy_lebedev_gridset
!!
!! FUNCTION
!!  Free the set of grids stored in this module. 
!!
!! PARENTS
!!      m_screening
!!
!! CHILDREN
!!      ld0006,ld0014,ld0026,ld0038,ld0050,ld0074,ld0086,ld0110,ld0146,ld0170
!!      ld0194,ld0230,ld0266,ld0302,ld0350,ld0434,ld0590,ld0770,ld0974,ld1202
!!      ld1454,ld1730,ld2030,ld2354,ld2702,ld3074,ld3470,ld3890,ld4334,ld4802
!!      ld5294,ld5810
!!
!! SOURCE

subroutine destroy_lebedev_gridset()

 use defs_basis

!Local variables-------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_lebedev_gridset'
!End of the abilint section

 integer :: igr
! *********************************************************************

 if (gridset_is_init) then
   do igr=1,lebedev_ngrids
     call destroy_lebedev_grid(Lgridset(igr))
   end do
   gridset_is_init=.FALSE.
 end if

end subroutine destroy_lebedev_gridset
!!***

!----------------------------------------------------------------------

!!****f* m_lebedev/m_lebedev_is_init
!! NAME
!!  m_lebedev_is_init
!!
!! FUNCTION
!!  Returns true if lebedev-laikov grids are already stored and computed.
!! 
!! SOURCE

function m_lebedev_is_init() result(ans)

 use defs_basis

!Local variables-------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'm_lebedev_is_init'
!End of the abilint section

 logical :: ans

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section
! *********************************************************************

 ans = gridset_is_init

end function m_lebedev_is_init
!!***

!----------------------------------------------------------------------

!!****f* m_lebedev/build_lebedev_grid
!! NAME
!!  build_lebedev_grid
!!
!! FUNCTION
!!  Helper subroutine returning the knots and weights of the angular grid
!!  from its sequential index
!!
!! INPUTS
!!  seq_idx=sequential index (must be in [1,32]
!!
!! OUTPUT
!!  xx(:),yy(:),zz(:)=The Cartesian coordinates of the knots.
!!  ww(:)=The weights.
!!
!! PARENTS
!!      m_lebedev
!!
!! CHILDREN
!!      ld0006,ld0014,ld0026,ld0038,ld0050,ld0074,ld0086,ld0110,ld0146,ld0170
!!      ld0194,ld0230,ld0266,ld0302,ld0350,ld0434,ld0590,ld0770,ld0974,ld1202
!!      ld1454,ld1730,ld2030,ld2354,ld2702,ld3074,ld3470,ld3890,ld4334,ld4802
!!      ld5294,ld5810
!!
!! SOURCE

subroutine build_lebedev_grid(seq_idx,xx,yy,zz,ww)

 use defs_basis

!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'build_lebedev_grid'
 use interfaces_28_numeric_noabirule
!End of the abilint section

 integer,intent(in) :: seq_idx
 real(dp) :: xx(:),yy(:),zz(:),ww(:)

!Local variables-------------------------------
 integer :: npts,ii
 character(len=500) :: msg
! *********************************************************************

 ii = assert_eq(SIZE(xx),SIZE(yy),SIZE(zz),SIZE(ww),"Wrong size")
 if (ii /= lebedev_npts(seq_idx)) then
   MSG_ERROR("wrong size in xx,yy,zz,ww")
 end if

 SELECT CASE (seq_idx)

 CASE (1)
   call LD0006(xx,yy,zz,ww,npts) 
 CASE (2)
   call LD0014(xx,yy,zz,ww,npts)
 CASE (3)
   call LD0026(xx,yy,zz,ww,npts)
 CASE (4)
   call LD0038(xx,yy,zz,ww,npts)
 CASE (5)
   call LD0050(xx,yy,zz,ww,npts)
 CASE (6)
   call LD0074(xx,yy,zz,ww,npts)
 CASE (7)
   call LD0086(xx,yy,zz,ww,npts)
 CASE (8)
   call LD0110(xx,yy,zz,ww,npts)
 CASE (9)
   call LD0146(xx,yy,zz,ww,npts)
 CASE (10)
   call LD0170(xx,yy,zz,ww,npts)
 CASE (11)
   call LD0194(xx,yy,zz,ww,npts)
 CASE (12)
   call LD0230(xx,yy,zz,ww,npts)
 CASE (13)
   call LD0266(xx,yy,zz,ww,npts)
 CASE (14)
   call LD0302(xx,yy,zz,ww,npts)
 CASE (15)
   call LD0350(xx,yy,zz,ww,npts)
 CASE (16)
   call LD0434(xx,yy,zz,ww,npts)
 CASE (17)
   call LD0590(xx,yy,zz,ww,npts)
 CASE (18)
   call LD0770(xx,yy,zz,ww,npts)
 CASE (19)
   call LD0974(xx,yy,zz,ww,npts)
 CASE (20)
   call LD1202(xx,yy,zz,ww,npts)
 CASE (21)
   call LD1454(xx,yy,zz,ww,npts)
 CASE (22)
   call LD1730(xx,yy,zz,ww,npts)
 CASE (23)
   call LD2030(xx,yy,zz,ww,npts)
 CASE (24)
   call LD2354(xx,yy,zz,ww,npts)
 CASE (25)
   call LD2702(xx,yy,zz,ww,npts)
 CASE (26)
   call LD3074(xx,yy,zz,ww,npts)
 CASE (27)
   call LD3470(xx,yy,zz,ww,npts)
 CASE (28)
   call LD3890(xx,yy,zz,ww,npts)
 CASE (29)
   call LD4334(xx,yy,zz,ww,npts)
 CASE (30)
   call LD4802(xx,yy,zz,ww,npts)
 CASE (31)
   call LD5294(xx,yy,zz,ww,npts)
 CASE (32)
   call LD5810(xx,yy,zz,ww,npts)

 CASE DEFAULT
   write(msg,'(a,i0)')"Wrong value for seq_idx: ",seq_idx
   MSG_ERROR(msg)
 END SELECT

end subroutine build_lebedev_grid
!!***

!----------------------------------------------------------------------

!!****f* m_lebedev/lebedev_quadrature 
!! NAME
!!  lebedev_quadrature
!!
!! FUNCTION
!!  Perform the integration of the complex function cplx_func on the sphere
!!  using lebedev-laikov angular grid in Cartesian coordinates.
!!  The routine improves the resolution of the grid until the required accuracy is reached
!!
!! INPUTS
!!  cplx_func(external)=Tthe function to be integrated
!!  int_pars(:)=integer parameter passed to cplx_func.
!!  real_pars(:)=real parameter passed to cplx_func.
!!  cplx_pars(:)=complex parameter passed to cplx_func.
!!  [accuracy]=fractional accuracy required. tol12 if not specified.
!!
!! OUTPUT
!!  quad=The integral 1/(4\pi) \int_\Omega func(\Omega) d\Omega
!!  ierr=0 if the quadrature converged. 
!!
!! NOTES
!!   commented because it causes problems to the new version of abilint
!!   Should work 
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

!!!!    subroutine lebedev_quadrature_cplx(cplx_func,int_pars,real_pars,cplx_pars,quad,ierr,accuracy)
!!!!    
!!!!     use defs_basis
!!!!    
!!!!    !Arguments ------------------------------------
!!!!    !scalars
!!!!    
!!!!    !This section has been created automatically by the script Abilint (TD).
!!!!    !Do not modify the following lines by hand.
!!!!    #undef ABI_FUNC
!!!!    #define ABI_FUNC 'lebedev_quadrature_cplx'
!!!!    !End of the abilint section
!!!!    
!!!!     integer,intent(out) :: ierr
!!!!     real(dp),optional,intent(in) :: accuracy
!!!!     complex(dpc),intent(out) :: quad
!!!!    !arrays
!!!!     integer,intent(in) :: int_pars(:)
!!!!     real(dp),intent(in) :: real_pars(:)
!!!!     complex(dpc),intent(in) :: cplx_pars(:)
!!!!     !complex(dpc) :: cplx_func
!!!!    
!!!!     interface
!!!!       function cplx_func(vers,int_pars,real_pars,cplx_pars)
!!!!         use defs_basis
!!!!         real(dp),intent(in) :: vers(3)
!!!!         integer,intent(in) :: int_pars(:)
!!!!         real(dp),intent(in) :: real_pars(:)
!!!!         complex(dpc),intent(in) :: cplx_pars(:)
!!!!         complex(dpc) ::  cplx_func ! abilint removes this line thus causing a compilation error
!!!!       end function cplx_func       ! for the time being, double complex is used to have the correct interface.
!!!!     end interface
!!!!    
!!!!    !Local variables-------------------------------
!!!!    !scalars
!!!!     integer :: igr,ipt
!!!!     real(dp) :: EPS,TOL
!!!!     complex(dpc) :: old_quad,abs_err
!!!!     character(len=500) :: msg      
!!!!    !arrays
!!!!     real(dp) :: vers(3)
!!!!    ! *************************************************************************
!!!!    
!!!!     ierr=0; TOL=tol12; EPS=tol6; if (PRESENT(accuracy)) EPS=accuracy
!!!!    
!!!!     do igr=1,lebedev_ngrids
!!!!       quad = czero
!!!!       do ipt=1,Lgridset(igr)%npts
!!!!         vers = Lgridset(igr)%versor(:,ipt)
!!!!         quad = quad + cplx_func(vers,int_pars,real_pars,cplx_pars) * Lgridset(igr)%weight(ipt)
!!!!       end do
!!!!       !write(std_out,'(a,i2,a,2es14.6)')" Lebedeb-Laikov grid # ",igr," quad= ",quad
!!!!       if (igr>1) then 
!!!!          if (ABS(quad-old_quad)<EPS*ABS(old_quad).or.(ABS(quad)<TOL.and.ABS(old_quad)<TOL)) RETURN
!!!!       end if
!!!!       abs_err = quad-old_quad
!!!!       old_quad = quad
!!!!     end do
!!!!    
!!!!     if (ABS(abs_err)<EPS*ABS(old_quad).or.(ABS(quad)<TOL.and.ABS(old_quad)<TOL)) then
!!!!       write(msg,'(2(a,es14.6),a,2(a,2es14.6))')&
!!!!    &    " Results are not converged within the given accuracy. EPS= ",EPS,"; TOL= ",TOL,ch10,&
!!!!    &    " Estimated absolute error= ",abs_err,"; relative error= ",abs_err/(ABS(quad)+tol16)
!!!!       MSG_WARNING(msg)
!!!!       ierr = -1
!!!!     end if
!!!!    
!!!!    end subroutine lebedev_quadrature_cplx
!!!!    !!***
!!!!    !----------------------------------------------------------------------

END MODULE m_lebedev
!!***
