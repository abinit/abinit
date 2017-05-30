!{\src2tex{textfont=tt}}
!!****f* ABINIT/spline_paw_fncs
!! NAME
!! spline_paw_fncs
!!
!! FUNCTION
!! Compute radial PAW functions and their derivatives on a set of points in the PAW sphere.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2017 ABINIT group (JJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! integer :: nnl : number of nl PAW basis functions in set
!! integer :: npts : number of points to perform fits on
!! real(dp) :: points(npts) : SORTED vector of points to perform fits on
!! type(pawrad_type) :: pawrad : paw radial mesh data
!! type(pawtab_type) :: pawtab : paw wavefunctions around each type of atom
!!
!! OUTPUT
!! real(dp) :: phi(npts,nnl), dphi(npts,nnl), tphi(npts,nnl), dtphi(npts,nnl) : PAW functions
!!             phi, tphi and their radial derivatives evaluated at the input points by spline fits
!!
!! NOTES
!! The PAW basis functions are defined by $<r|\phi_i>=(u_i(r)/r)S_{lm}(\hat{r})$, evaluated on a radial
!! grid. This subroutine computes $u(r)$ and $d u(r)/dr$ on the set of points provided on input. Typically
!! the input points will be the fine grid points in the PAW sphere. They are presumed to be sorted
!! already on input.
!!
!! PARENTS
!!
!! CHILDREN
!!      nderiv_gen,spline,splint
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine spline_paw_fncs(dphi,dtphi,nnl,npts,pawrad,pawtab,points,phi,tphi)

 use m_profiling_abi

 use defs_basis
 use m_errors
 use m_splines
 use m_pawrad, only : pawrad_type, nderiv_gen
 use m_pawtab, only : pawtab_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spline_paw_fncs'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nnl,npts
!arrays
 real(dp),intent(in) :: points(npts)
 real(dp),intent(out) :: phi(npts,nnl),dphi(npts,nnl),tphi(npts,nnl),dtphi(npts,nnl)
 type(pawrad_type),intent(in) :: pawrad
 type(pawtab_type),intent(in) :: pawtab

!Local variables-------------------------------
!scalars
 integer :: inl
 real(dp) :: ybcbeg, ybcend
!arrays
 real(dp),allocatable :: der(:),diag(:),ypp(:)

! ************************************************************************

 DBG_ENTER("COLL")

 ABI_ALLOCATE(der,(pawtab%mesh_size))
 ABI_ALLOCATE(ypp,(pawtab%mesh_size))
 ABI_ALLOCATE(diag,(pawtab%mesh_size))
 do inl = 1, nnl

!  spline phi onto points
   ypp(:) = zero; diag(:) = zero; ybcbeg = zero; ybcend = zero
   call spline(pawrad%rad,pawtab%phi(:,inl),pawtab%mesh_size,ybcbeg,ybcend,ypp)
   call splint(pawtab%mesh_size,pawrad%rad,pawtab%phi(:,inl),ypp,npts,points,phi(:,inl))

!  next spline d phi/dr onto points
!  need derivative of phi with respect to radius
   der(:) = zero
   call nderiv_gen(der,pawtab%phi(:,inl),pawrad)
   ypp(:) = zero; diag(:) = zero; ybcbeg = zero; ybcend = zero
   call spline(pawrad%rad,der,pawtab%mesh_size,ybcbeg,ybcend,ypp)
   call splint(pawtab%mesh_size,pawrad%rad,der,ypp,npts,points,dphi(:,inl))

!  next splint tphi onto points
   ypp(:) = zero; diag(:) = zero;
   call spline(pawrad%rad,pawtab%tphi(:,inl),pawtab%mesh_size,ybcbeg,ybcend,ypp)
   call splint(pawtab%mesh_size,pawrad%rad,pawtab%tphi(:,inl),ypp,npts,points,tphi(:,inl))

!  finally spline d tphi/dr onto points
!  need derivative of tphi with respect to radius
   der(:) = zero
   call nderiv_gen(der,pawtab%tphi(:,inl),pawrad)
   ypp(:) = zero; diag(:) = zero; ybcbeg = zero; ybcend = zero
   call spline(pawrad%rad,der,pawtab%mesh_size,ybcbeg,ybcend,ypp)
   call splint(pawtab%mesh_size,pawrad%rad,der,ypp,npts,points,dtphi(:,inl))

 end do ! end loop over nnl basis functions

 ABI_DEALLOCATE(der)
 ABI_DEALLOCATE(ypp)
 ABI_DEALLOCATE(diag)

 DBG_EXIT("COLL")

 end subroutine spline_paw_fncs
!!***
