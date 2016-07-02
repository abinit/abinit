!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp11nl
!! NAME
!! psp11nl
!!
!! FUNCTION
!! Fourier transform the real space UPF projector functions to reciprocal space
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  lmax=maximum ang momentum for which nonlocal form factor is desired.
!!   Usually lmax=1, sometimes = 0 (e.g. for oxygen); lmax <= 2 allowed.
!!  mmax=number of radial grid points for atomic grid
!!  lnmax= maximum index for all l channel projectors, dimension of ffspl
!!  lmnmax= maximum index for all projectors, dimension of indlmn
!!  mqgrid=number of grid points for q grid
!!  n_proj = total number of NL projectors read in
!!  proj = projector data times r, on a real space grid
!!  proj_l = ang mom channel for each projector
!!  proj_np = max number of points used for each projector
!!  qgrid(mqgrid)=values at which form factors are returned
!!  r(mmax)=radial grid values
!!  drdi=derivative of grid point wrt index
!!  useylm = input to use m dependency of NL part, or only Legendre polynomials
!!
!! OUTPUT
!!  ffspl(mqgrid,2,mpsang)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum
!!  indlmn = indexing of each projector, for n, l, m, s, ln, lmn (see pspatm.F90)
!!
!! NOTES
!!
!! PARENTS
!!      upf2abinit
!!
!! CHILDREN
!!      ctrap,jbessel,spline
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psp11nl(ffspl,indlmn,mmax,lnmax,lmnmax,mqgrid,n_proj,&
&                  proj, proj_l, proj_np, qgrid, r, drdi, useylm)

 use defs_basis
 use m_splines
 use m_profiling_abi
 use m_errors

 use m_paw_numeric, only: jbessel=>paw_jbessel

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp11nl'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mmax, lnmax, lmnmax, mqgrid, useylm, n_proj
!arrays
 integer, intent(in) :: proj_l(n_proj)
 integer, intent(in) :: proj_np(n_proj)
 integer, intent(out) :: indlmn(6,lmnmax)
 real(dp),intent(in) :: r(mmax)
 real(dp),intent(in) :: drdi(mmax)
 real(dp),intent(in) :: proj(mmax,n_proj)
 real(dp),intent(in) :: qgrid(mqgrid)
 real(dp),intent(inout) :: ffspl(mqgrid,2,lnmax) !vz_i

!Local variables-------------------------------
!scalars
 integer :: iproj, np, ll, llold, ipsang, i_indlmn
 integer :: iproj_1l, ir, iq, mm
 integer :: bessorder
 real(dp) :: res, arg, besfact, dummy, dummy2
 real(dp), allocatable :: work(:)
 character(len=500) :: message

!*************************************************************************
 bessorder = 0 ! never calculate derivatives of bessel functions

 ffspl = zero
 indlmn = 0
 i_indlmn = 0
 llold = -1
 iproj_1l = 1
!big loop over all projectors
 do iproj = 1, n_proj
   
   if (iproj > lmnmax) then
     write(message,'(a,2i0)') ' Too many projectors found. n_proj, lmnmax =  ',n_proj, lmnmax 
     MSG_ERROR(message)
   end if

   np = proj_np(iproj)
   ABI_ALLOCATE(work,(np))
   ll = proj_l(iproj)
   if (ll < llold) then
     message = 'psp11nl : Error: UPF projectors are not in order of increasing ll'
     MSG_ERROR(message)
   else if (ll == llold) then
     iproj_1l = iproj_1l + 1
   else
     iproj_1l = 1
     llold = ll
   end if
!  determine indlmn for this projector (keep in UPF order and enforce that they are in
!  increasing ll)
   do mm = 1, 2*ll*useylm+1
     i_indlmn = i_indlmn + 1
     indlmn(1,i_indlmn) = ll
     indlmn(2,i_indlmn) = mm-ll*useylm-1
     indlmn(3,i_indlmn) = iproj_1l
     indlmn(4,i_indlmn) = ll*ll+(1-useylm)*ll+mm
     indlmn(5,i_indlmn) = iproj
     indlmn(6,i_indlmn) = 1 !spin? FIXME: to get j for relativistic cases
   end do

!  FT projectors to reciprocal space q
   do iq = 1, mqgrid
     arg = two_pi*qgrid(iq)

!    FIXME: add semianalytic form for integral from 0 to first point
     do ir = 1, np
       call jbessel(besfact, dummy, dummy2, ll, bessorder, arg*r(ir))
!      besfact = sin(arg*r(ir))
       work(ir) = drdi(ir) * besfact * proj(ir, iproj) * r(ir) !* r(ir)
     end do
     call ctrap (np, work, one, res)

     ffspl(iq, 1, iproj) = res
   end do
   ABI_DEALLOCATE(work)
 end do  ! iproj

!add derivative of ffspl(:,1,:) for spline interpolation later
 ABI_ALLOCATE(work,(mqgrid))
 do ipsang = 1, lnmax
   call spline(qgrid,ffspl(:,1,ipsang),mqgrid,zero,zero,ffspl(:,2,ipsang))
 end do
 ABI_DEALLOCATE(work)

end subroutine psp11nl
!!***
