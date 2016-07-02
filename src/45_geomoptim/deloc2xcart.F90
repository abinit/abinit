!{\src2tex{textfont=tt}}
!!****f* ABINIT/deloc2xcart
!! NAME
!! deloc2xcart
!!
!! FUNCTION
!!  Determine the cartesian coordinates which correspond to the
!!  given values of the delocalized coordinates. The relationship
!!  is non-linear, so use an iterative scheme, as in Baker
!!  JCP .105. 192 (1996).
!!  Older reference: Pulay and co. JACS 101 2550 (1979)
!!
!! COPYRIGHT
!! Copyright (C) 2003-2016 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!   deloc <type(delocint)>=Important variables for
!!   |                           pred_delocint
!!   |
!!   | nang     = Number of angles
!!   | nbond    = Number of bonds
!!   | ncart    = Number of cartesian directions
!!   |             (used for constraints)
!!   | ndihed   = Number of dihedrals
!!   | nrshift  = Dimension of rshift
!!   | ninternal= Number of internal coordinates
!!   |            ninternal=nbond+nang+ndihed+ncart
!!   |
!!   | angs(2,3,nang)  = Indexes to characterize angles
!!   | bonds(2,2,nbond)= For a bond between iatom and jatom
!!   |                   bonds(1,1,nbond) = iatom
!!   |                   bonds(2,1,nbond) = icenter
!!   |                   bonds(1,2,nbond) = jatom
!!   |                   bonds(2,2,nbond) = irshift
!!   | carts(2,ncart)  = Index of total primitive internal,
!!   |                   and atom (carts(2,:))
!!   | dihedrals(2,4,ndihed)= Indexes to characterize dihedrals
!!   |
!!   | rshift(3,nrshift)= Shift in xred that must be done to find
!!   |                    all neighbors of a given atom within a
!!   |                    given number of neighboring shells
!! natom = Number of atoms (dtset%natom)
!! rprimd(3,3)=dimensional real space primitive translations (bohr)
!!
!! OUTPUT
!! bt_inv_matrix(3*(natom-1),3*natom)=inverse of transpose of B matrix
!!
!! SIDE EFFECTS
!! u_matrix(ninternal,3*(natom-1))=eigenvectors of G = BB^T matrix
!! xcart(3,natom)=cartesian coordinates of atoms (bohr)
!!
!! NOTES
!!
!! PARENTS
!!      pred_delocint
!!
!! CHILDREN
!!      dgemv,wrtout,xcart2deloc
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine deloc2xcart(deloc,natom,rprimd,xcart,deloc_int,btinv,u_matrix)

 use defs_basis
 use m_abimover
 use m_errors
 use m_profiling_abi
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deloc2xcart'
 use interfaces_14_hidewrite
 use interfaces_45_geomoptim, except_this_one => deloc2xcart
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom
 type(delocint),intent(inout) :: deloc
!arrays
 real(dp),intent(in) :: deloc_int(3*(natom-1)),rprimd(3,3)
 real(dp),intent(inout) :: u_matrix(deloc%ninternal,3*(natom-1))
 real(dp),intent(inout) :: xcart(3,natom)
 real(dp),intent(out) :: btinv(3*(natom-1),3*natom)

!Local variables-------------------------------
!scalars
 integer :: iiter,iprim,niter
 integer :: ii
 real(dp) :: minmix, maxmix
 real(dp) :: mix,tot_diff, toldeloc
 real(dp) :: lntoldeloc
 logical  :: DEBUG=.FALSE.
!arrays
 real(dp) :: btinv_tmp(3*(natom-1),3*natom)
 real(dp) :: cgrad(3*natom),cgrad_old(3*natom)
 real(dp) :: deloc_int_now(3*(natom-1)),prim_int(deloc%ninternal)
 real(dp) :: tmpxcart(3*natom)
 real(dp) :: xdeloc_diff(3*(natom-1))

 character(len=500) :: message

! ******************************************************************

 if (DEBUG) then
   write(ab_out,*) 'ENTERING DELOC2XCART'

   write (message,*) 'BONDS=',deloc%nbond
   call wrtout(ab_out,message,'COLL')
   do ii = 1, deloc%nbond
     write (message,*) ii, deloc%bonds(:,:,ii)
     call wrtout(ab_out,message,'COLL')
   end do

   write (message,*) 'ANGS=',deloc%nang
   call wrtout(ab_out,message,'COLL')
   do ii = 1, deloc%nang
     write (message,*) ii, deloc%angs(:,:,ii)
     call wrtout(ab_out,message,'COLL')
   end do

   write (message,*) 'DIHEDRALS=',deloc%ndihed
   call wrtout(ab_out,message,'COLL')
   do ii = 1, deloc%ndihed
     write (message,*) ii, deloc%dihedrals(:,:,ii)
     call wrtout(ab_out,message,'COLL')
   end do

   write (message,*) 'CARTS=',deloc%ncart
   call wrtout(ab_out,message,'COLL')
   do ii = 1, deloc%ncart
     write (message,*) ii, deloc%carts(:,ii)
     call wrtout(ab_out,message,'COLL')
   end do

   write (ab_out,*) 'xcart (input)'
   do ii=1,natom
     write (ab_out,*) xcart(:,ii)
   end do

 end if

 niter = 200
 tmpxcart = reshape(xcart,(/3*natom/))

 cgrad_old(:) = zero
 cgrad(:) = zero
 maxmix = 0.9_dp
 minmix = 0.2_dp
 toldeloc = tol10
 lntoldeloc = log(toldeloc)

 do iiter=1,niter
   if (iiter==1) then
     mix= minmix
   else
     mix = minmix + (maxmix-minmix)*(log(tot_diff)-lntoldeloc) / lntoldeloc
   end if
   if (mix < minmix) mix = minmix
   if (mix > maxmix) mix = maxmix

   tmpxcart(:) = tmpxcart(:) + mix*cgrad(:)
   xcart = reshape(tmpxcart,(/3,natom/))
   call xcart2deloc(deloc,natom,rprimd,xcart,&
&   btinv_tmp,u_matrix,deloc_int_now,prim_int)
!  update the BT^{-1} matrix?
   btinv(:,:) = btinv_tmp(:,:)

   xdeloc_diff(:) = deloc_int(:) - deloc_int_now(:)

   tot_diff = sum(abs(xdeloc_diff))
   if (tot_diff < toldeloc) exit

   cgrad_old(:) = cgrad(:)

!  gradient vector = btinv^{T} * xdeloc_diff
   call dgemv('T',3*(natom-1),3*natom,one,&
&   btinv,3*(natom-1),xdeloc_diff,1,zero,cgrad,1)
 end do
!end iiter do

 call xcart2deloc(deloc,natom,rprimd,xcart,&
& btinv,u_matrix,deloc_int_now,prim_int)
 write (message,'(3a)') 'delocalized internals, after convergence of xcart = ', ch10
 call wrtout(std_out,message,'COLL')
 do ii = 1, 3*(natom-1)
   write (message,'(I6,E20.10,2x)') ii, deloc_int_now(ii)
   call wrtout(std_out,message,'COLL')
 end do

 xdeloc_diff(:) = deloc_int(:) - deloc_int_now(:)

 write (message,'(a)') 'Primitive internal coordinate values:'
 call wrtout(std_out,message,'COLL')
 do iprim = 1, deloc%nbond
   write (message,'(i6,E20.10)') iprim, prim_int(iprim)
   call wrtout(std_out,message,'COLL')
 end do
 do iprim = deloc%nbond+1, deloc%nbond+deloc%nang+deloc%ndihed
   write (message,'(i6,2E20.10)') iprim, prim_int(iprim), prim_int(iprim)/pi*180.0_dp
   call wrtout(std_out,message,'COLL')
 end do
 do iprim = deloc%nbond+deloc%nang+deloc%ndihed+1, deloc%ninternal
   write (message,'(i6,E20.10)') iprim, prim_int(iprim)
   call wrtout(std_out,message,'COLL')
 end do

 if (iiter == niter+1) then
   write (message,'(a,i6,a,E20.10)') 'deloc2xcart : Error, xcart not converged in ', niter, 'iterations ', tot_diff
   MSG_ERROR(message)
 end if

 if(DEBUG)then
   write (ab_out,*) 'xcart (output)'
   do ii=1,natom
     write (ab_out,*) xcart(:,ii)
   end do
   write(ab_out,*) 'EXITING DELOC2XCART'
 end if

end subroutine deloc2xcart
!!***
