!{\src2tex{textfont=tt}}
!!****f* ABINIT/pred_verlet
!! NAME
!! pred_verlet
!!
!! FUNCTION
!! Ionmov predictors (15) FIRE algorithm
!!
!! IONMOV 15:
!! Given a starting point xred that is a vector of length 3*natom
!! (reduced nuclei coordinates), a velocity vector (in cartesian
!! coordinates), and unit cell parameters (acell and rprimd -
!! without velocities in the present implementation), the Verlet
!! dynamics is performed, using the gradient of the energy
!! (atomic forces and stresses) as calculated by the routine scfcv.
!!
!! At each step, the dot product of velocity and force is calculated.
!! If the v.dot.f is positive for min_downhill consecutive steps, the
!! time step dtion is increased by dtinc until it reaches dtmax. If
!! v.dot.f is negative, dtion is mulitiplied by dtdec.
!!
!! Some atoms can be kept fixed, while the propagation of unit cell
!! parameters is only performed if optcell/=0.
!! No more than ab_mover%ntime steps are performed.
!! The time step is governed by dtion, contained in ab_mover
!! (coming from dtset).
!! Returned quantities are xred, and eventually acell and rprimd
!! (new ones!).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (XHe)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,
!! see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! ab_mover <type(abimover)> : Datatype with all the information
!!                                needed by the preditor
!! itime  : Index of the present iteration
!! ntime  : Maximal number of iterations
!! ionmov : (6 or 7) Specific kind of VERLET
!! zDEBUG : if true print some debugging information
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! hist <type(abihist)> : History of positions,forces
!!                               acell, rprimd, stresses
!!
!! NOTES
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      fcart2fred,hist2var,metric,mkrdim,var2hist,wrtout,xcart2xred
!!      xfpack_f2vout,xfpack_vin2x,xfpack_x2vin,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"
subroutine pred_fire(ab_mover, ab_xfh,preconforstr, hist, itime, ntime, DEBUG, iexit)

 use defs_basis
 use m_profiling_abi
 use m_abimover
 use m_abihist

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pred_fire'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
 use interfaces_45_geomoptim, except_this_one => pred_fire
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(abimover),intent(in) :: ab_mover
 type(ab_xfh_type),intent(inout)    :: ab_xfh
 type(abiforstr),intent(in) :: forstr
 type(abihist),intent(inout) :: hist
 integer,intent(in) :: itime
 integer,intent(in) :: ntime
 integer,intent(in) :: ionmov
 integer,intent(in) :: iexit
 logical,intent(in) :: zDEBUG


!Local variables-------------------------------
!scalars
integer  :: ihist_prev,ndim,
integer, parameter :: npul=0
integer  :: ierr,ii,jj,kk
real(dp),save :: ucvol0
real(dp) :: ucvol,det
real(dp) :: etotal,etotal_prev
real(dp) :: favg
real(dp),save :: hh
real (dp) :: dtinc, dtdec, dtmax, alphadec
integer, save :: ndownhill
! time step is dtion
real(dp), save :: hh

!arrays
real(dp) :: acell(3),strten(6)
real(dp) :: rprim(3,3),rprimd(3,3)
real(dp) :: xred(3,ab_mover%natom),xcart(3,ab_mover%natom)
real(dp), save :: vel(3,ab_mover%natom)
real(dp) :: residual(3,ab_mover%natom),residual_corrected(3,ab_mover%natom)
real(dp),allocatable,save :: vin(:),vout(:)

! PaRaMeters
dtinc=1.1_dp
dtdec=0.5_dp
dtmax=1.0_dp
alphadec=0.99_dp



!***************************************************************************
!Beginning of executable session
!***************************************************************************

 if(iexit/=0)then
   if (allocated(vin))           then
     ABI_DEALLOCATE(vin)
   end if
   if (allocated(vout))          then
     ABI_DEALLOCATE(vout)
   end if
   return
 end if

!##########################################################
!### 01. Compute the dimension of vectors (ndim)

 ndim=3*ab_mover%natom
 if(ab_mover%optcell==1 .or.&
& ab_mover%optcell==4 .or.&
& ab_mover%optcell==5 .or.&
& ab_mover%optcell==6) ndim=ndim+1
 if(ab_mover%optcell==2 .or.&
& ab_mover%optcell==3) ndim=ndim+6
 if(ab_mover%optcell==7 .or.&
& ab_mover%optcell==8 .or.&
& ab_mover%optcell==9) ndim=ndim+3

!##########################################################
!### 02. Allocate the vectors vin

!Notice that vin, vout, etc could be allocated
!From a previous dataset with a different ndim
 if(itime==1)then
   if (allocated(vin))           then
     ABI_DEALLOCATE(vin)
   end if
   if (allocated(vout))          then
     ABI_DEALLOCATE(vout)
   end if
   ABI_ALLOCATE(vin,(ndim))
   ABI_ALLOCATE(vout,(ndim))

! initialize values
   ndownhill=0
   if (ab_mover%dtion>0)then
     hh = ab_mover%dtion
   end if
 end if

!##########################################################
!### 03. Obtain the present values from the history

 call hist2var(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 do ii=1,3
   rprim(ii,1:3)=rprimd(ii,1:3)/acell(1:3)
 end do

 call xred2xcart(ab_mover%natom,rprimd,xcart,xred)
 fcart(:,:)  =hist%fcart(:,:,hist%ihist)
 strten(:)  =hist%strten(:,hist%ihist)
 vel(:,:)   =hist%vel(:,:,hist%ihist)
 etotal     =hist%etot(hist%ihist)


!Fill the residual with forces (No preconditioning)
!Or the preconditioned forces
 if (ab_mover%goprecon==0)then
   call fcart2fred(hist%fcart(:,:,hist%ihist),residual,rprimd,ab_mover%natom)
 else
   residual(:,:)=forstr%fred(:,:)
 end if

!Save initial values
 if (itime==1)then
   acell0(:)=acell(:)
   rprimd0(:,:)=rprimd(:,:)
   ucvol0=ucvol
 end if


 if(zDEBUG)then
   write (std_out,*) 'fcart:'
   do kk=1,ab_mover%natom
     write (std_out,*) fcart(:,kk)
   end do
   write (std_out,*) 'vel:'
   do kk=1,ab_mover%natom
     write (std_out,*) vel(:,kk)
   end do
   write (std_out,*) 'strten:'
   write (std_out,*) strten(1:3),ch10,strten(4:6)
   write (std_out,*) 'etotal:'
   write (std_out,*) etotal
 end if

!Get rid of mean force on whole unit cell, but only if no
!generalized constraints are in effect

 residual_corrected(:,:)=residual(:,:)
 if(ab_mover%nconeq==0)then
   amass_tot=sum(ab_mover%amass(:))
   do kk=1,3
     if (kk/=3.or.ab_mover%jellslab==0) then
       favg=sum(fred_corrected(kk,:))/dble(ab_mover%natom)
       fred_corrected(kk,:)=fred_corrected(kk,:)-favg*ab_mover%amass(:)/amass_tot
     end if
   end do
 end if

!##########################################################
!### 04. Fill the vectors vin and vout

!Initialize input vectors : first vin, then vout
 call xfpack_x2vin(acell, acell0, ab_mover%natom, ndim,&
& ab_mover%nsym, ab_mover%optcell, rprim, rprimd,&
& ab_mover%symrel, ucvol, ucvol0, vin, xred)
 call xfpack_f2vout(fred_corrected, ab_mover%natom, ndim,&
& ab_mover%optcell, ab_mover%strtarget, strten, ucvol,&
& vout)

end subroutine pred_fire
!!***
