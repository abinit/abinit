!{\src2tex{textfont=tt}}
!!****f* ABINIT/xfh_recover_deloc
!! NAME
!! xfh_recover_deloc
!!
!! FUNCTION
!! Update the contents of the history xfhist taking values
!! from xred, acell, rprim, fred_corrected and strten
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR, JCC, SE)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,
!! see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!!
!! SIDE EFFECTS
!!
!!
!! PARENTS
!!      pred_delocint
!!
!! CHILDREN
!!      fred2fdeloc,hessupdt,xcart2deloc,xfpack_f2vout,xfpack_x2vin,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine xfh_recover_deloc(ab_xfh,ab_mover,acell,acell0,cycl_main,&
& fred,hessin,ndim,rprim,rprimd0,strten,ucvol,ucvol0,vin,vin_prev,&
& vout,vout_prev,xred,deloc,deloc_int,deloc_force,btinv,gprimd,prim_int,&
& u_matrix)

 use m_profiling_abi
 use defs_basis
 use m_abimover

 use m_results_gs , only : results_gs_type
 use m_bfgs,        only : hessupdt

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xfh_recover_deloc'
 use interfaces_41_geometry
 use interfaces_45_geomoptim, except_this_one => xfh_recover_deloc
!End of the abilint section

implicit none

!Arguments ------------------------------------
!scalars

integer,intent(in) :: ndim
integer,intent(out) :: cycl_main
real(dp),intent(inout) :: ucvol,ucvol0
type(ab_xfh_type),intent(inout) :: ab_xfh
type(abimover),intent(in) :: ab_mover
! DELOCINT specials
type(delocint),intent(inout) :: deloc

!arrays
real(dp),intent(inout) :: acell(3)
real(dp),intent(in) :: acell0(3)
real(dp),intent(inout) :: hessin(:,:)
real(dp),intent(inout) :: xred(3,ab_mover%natom)
real(dp),intent(inout) :: rprim(3,3)
real(dp),intent(inout) :: rprimd0(3,3)
real(dp),intent(inout) :: fred(3,ab_mover%natom)
real(dp),intent(inout) :: strten(6)
real(dp),intent(inout) :: vin(:)
real(dp),intent(inout) :: vin_prev(:)
real(dp),intent(inout) :: vout(:)
real(dp),intent(inout) :: vout_prev(:)
! DELOCINT specials
real(dp),intent(inout) :: deloc_force(3*(ab_mover%natom-1))
real(dp),intent(inout) :: deloc_int(3*(ab_mover%natom-1))
real(dp),intent(inout) :: btinv(3*(ab_mover%natom-1),3*ab_mover%natom)
real(dp),intent(inout) :: prim_int(:),u_matrix(:,:),gprimd(3,3)

!Local variables-------------------------------
!scalars
integer :: ixfh
real(dp) :: xcart(3,ab_mover%natom)

!*********************************************************************

 if(ab_xfh%nxfh/=0)then
!  Loop over previous time steps
   do ixfh=1,ab_xfh%nxfh

!    For that time step, get new (x,f) from xfhist
     xred(:,:)     =ab_xfh%xfhist(:,1:ab_mover%natom        ,1,ixfh)
     rprim(1:3,1:3)=ab_xfh%xfhist(:,ab_mover%natom+2:ab_mover%natom+4,1,ixfh)
     acell(:)      =ab_xfh%xfhist(:,ab_mover%natom+1,1,ixfh)
     fred(:,:)     =ab_xfh%xfhist(:,1:ab_mover%natom,2,ixfh)
!    This use of results_gs is unusual
     strten(1:3)   =ab_xfh%xfhist(:,ab_mover%natom+2,2,ixfh)
     strten(4:6)   =ab_xfh%xfhist(:,ab_mover%natom+3,2,ixfh)

!    !DEBUG
!    write (ab_out,*) '---READ FROM XFHIST---'

!    write (ab_out,*) 'XRED'
!    do kk=1,ab_mover%natom
!    write (ab_out,*) xred(:,kk)
!    end do
!    write (ab_out,*) 'FRED'
!    do kk=1,ab_mover%natom
!    write (ab_out,*) fred(:,kk)
!    end do
!    write(ab_out,*) 'RPRIM'
!    do kk=1,3
!    write(ab_out,*) rprim(:,kk)
!    end do
!    write(ab_out,*) 'ACELL'
!    write(ab_out,*) acell(:)
!    !DEBUG

!    Convert input xred (reduced coordinates) to xcart (cartesian)
     call xred2xcart(ab_mover%natom,rprimd0,xcart,xred)
!    Convert input coordinates in Delocalized internals
     call xcart2deloc(deloc,ab_mover%natom,rprimd0,xcart,&
&     btinv,u_matrix,deloc_int,prim_int)
!    Convert forces to delocalized coordinates for next step
     call fred2fdeloc(btinv,deloc_force,fred,ab_mover%natom,gprimd)

!    Transfer it in vin, vout
     call xfpack_x2vin(acell,acell0,ab_mover%natom-1,&
&     ndim,ab_mover%nsym,ab_mover%optcell,rprim,rprimd0,&
&     ab_mover%symrel,ucvol,ucvol0,vin,deloc_int)
     call xfpack_f2vout(deloc_force,ab_mover%natom-1,&
&     ndim,ab_mover%optcell,ab_mover%strtarget,strten,&
&     ucvol,vout)
!    Get old time step, if any, and update inverse hessian
     if(ixfh/=1)then
       xred(:,:)     =ab_xfh%xfhist(:,1:ab_mover%natom,1,ixfh-1)
       rprim(1:3,1:3)=&
&       ab_xfh%xfhist(:,ab_mover%natom+2:ab_mover%natom+4,1,ixfh-1)
       acell(:)=ab_xfh%xfhist(:,ab_mover%natom+1,1,ixfh-1)
       fred(:,:)=ab_xfh%xfhist(:,1:ab_mover%natom,2,ixfh-1)
!      This use of results_gs is unusual
       strten(1:3)=ab_xfh%xfhist(:,ab_mover%natom+2,2,ixfh-1)
       strten(4:6)=ab_xfh%xfhist(:,ab_mover%natom+3,2,ixfh-1)

!      Convert input xred (reduced coordinates) to xcart (cartesian)
       call xred2xcart(ab_mover%natom,rprimd0,xcart,xred)
!      Convert input coordinates in Delocalized internals
       call xcart2deloc(deloc,ab_mover%natom,rprimd0,xcart,&
&       btinv,u_matrix,deloc_int,prim_int)
!      Convert forces to delocalized coordinates for next step
       call fred2fdeloc(btinv,deloc_force,fred,ab_mover%natom,gprimd)

!      Tranfer it in vin_prev, vout_prev
       call xfpack_x2vin(acell,acell0,ab_mover%natom-1,&
&       ndim,ab_mover%nsym,ab_mover%optcell,rprim,rprimd0,&
&       ab_mover%symrel,ucvol,ucvol0,vin_prev,deloc_int)
       call xfpack_f2vout(deloc_force,ab_mover%natom-1,&
&       ndim,ab_mover%optcell,ab_mover%strtarget,strten,&
&       ucvol,vout_prev)

!      write(ab_out,*) 'Hessian matrix before update',ndim,'x',ndim
!      write(ab_out,*) 'ixfh=',ixfh
!      do kk=1,ndim
!      do jj=1,ndim,3
!      if (jj+2<=ndim)then
!      write(ab_out,*) jj,hessin(jj:jj+2,kk)
!      else
!      write(ab_out,*) jj,hessin(jj:ndim,kk)
!      end if
!      end do
!      end do

       call hessupdt(hessin,ab_mover%iatfix,ab_mover%natom-1,ndim,&
&       vin,vin_prev,vout,vout_prev)

!      !DEBUG
!      write(ab_out,*) 'Hessian matrix after update',ndim,'x',ndim
!      do kk=1,ndim
!      do jj=1,ndim,3
!      if (jj+2<=ndim)then
!      write(ab_out,*) jj,hessin(jj:jj+2,kk)
!      else
!      write(ab_out,*) jj,hessin(jj:ndim,kk)
!      end if
!      end do
!      end do
!      !DEBUG

     end if !if(ab_xfh%nxfh/=0)

!    End loop over previous time steps
   end do

!  The hessian has been generated,
!  as well as the latest vin and vout
!  so will cycle the main loop
   cycl_main=1

 end if

end subroutine xfh_recover_deloc
!!***
