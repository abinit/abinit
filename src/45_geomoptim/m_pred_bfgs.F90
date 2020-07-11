!!****m* ABINIT/m_pred_bfgs
!! NAME
!!  m_pred_bfgs
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, JCC, SE, FB)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_pred_bfgs

 use defs_basis
 use m_abicore
 use m_abimover
 use m_abihist
 use m_xfpack
 use m_lbfgs
 use m_errors

 use m_geometry,    only : mkrdim, fcart2fred, metric
 use m_bfgs,        only : hessinit, hessupdt, brdene

 implicit none

 private
!!***

 public :: pred_bfgs
 public :: pred_lbfgs
!!***

contains
!!***

!!****f* ABINIT/pred_bfgs
!! NAME
!! pred_bfgs
!!
!! FUNCTION
!! Ionmov predictors (2 & 3) Broyden-Fletcher-Goldfarb-Shanno
!!
!! IONMOV 2:
!! Given a starting point xred that is a vector of length 3*natom
!! (reduced nuclei coordinates), and unit cell parameters
!! (acell and rprimd) the Broyden-Fletcher-Goldfarb-Shanno
!! minimization is performed on the total energy function, using
!! its gradient (atomic forces and stresse) as calculated
!! by the routine scfcv. Some atoms can be kept fixed,
!! while the optimization of unit cell parameters is
!! only performed if optcell/=0. The convergence requirement on
!! the atomic forces, dtset%tolmxf,  allows an early exit.
!! Otherwise no more than dtset%ntime steps are performed.
!! Returned quantities are xred, and eventually acell and rprim (new ones!).
!! Could see Numerical Recipes (Fortran), 1986, page 307.
!!
!! IONMOV 3:
!! Conduct structural optimization using the Broyden-Fletcher-
!! Goldfarb-Shanno minimization (BFGS), modified to take into
!! account the total energy as well as the gradients (as in usual
!! BFGS). See the paper by Schlegel, J. Comp. Chem. 3, 214 (1982) [[cite:Schlegel1982]].
!! Might be better than ionmov=2 for few degrees of freedom (less than 3 or 4)
!!
!! INPUTS
!! ab_mover <type(abimover)> : Datatype with all the information
!!                                needed by the preditor
!! itime  : Index of the present iteration
!! ntime  : Maximal number of iterations
!! ionmov : (2 or 3) Specific kind of BFGS
!! zDEBUG : if true print some debugging information
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! hist <type(abihist)> : History of positions,forces acell, rprimd, stresses
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      brdene,dgetrf,dgetri,fcart2fred,hessinit,hessupdt,hist2var,metric
!!      mkrdim,var2hist,xfh_recover_new,xfpack_f2vout,xfpack_vin2x,xfpack_x2vin
!!
!! SOURCE

subroutine pred_bfgs(ab_mover,ab_xfh,forstr,hist,ionmov,itime,zDEBUG,iexit)

implicit none

!Arguments ------------------------------------
!scalars
type(abimover),intent(in)       :: ab_mover
type(ab_xfh_type),intent(inout)    :: ab_xfh
type(abihist),intent(inout) :: hist
type(abiforstr),intent(in) :: forstr
integer,intent(in) :: itime
integer,intent(in) :: ionmov
integer,intent(in) :: iexit
logical,intent(in) :: zDEBUG

!Local variables-------------------------------
!scalars
integer  :: ihist_prev,ndim,cycl_main
integer, parameter :: npul=0
integer  :: ierr,ii,jj,kk,nitpul
real(dp),save :: ucvol0
real(dp) :: ucvol,det
real(dp) :: etotal,etotal_prev
real(dp) :: favg,alpha0

!arrays
integer,allocatable :: ipiv(:)
real(dp),allocatable,save :: hessin(:,:),vin(:),vin_prev(:)
real(dp),allocatable,save :: vout(:),vout_prev(:)
real(dp),allocatable,save ::vinres(:,:),vin1(:,:)
real(dp),allocatable ::  amat(:,:),amatinv(:,:),alpha(:,:)
real(dp),allocatable :: rwork(:)
real(dp),save :: acell0(3) ! Initial acell
real(dp),save :: rprimd0(3,3) ! Initial rprimd
real(dp) :: acell(3)
real(dp) :: rprimd(3,3),rprim(3,3)
real(dp) :: gprimd(3,3)
real(dp) :: gmet(3,3)
real(dp) :: rmet(3,3)
real(dp) :: residual(3,ab_mover%natom),residual_corrected(3,ab_mover%natom)
real(dp) :: xred(3,ab_mover%natom),strten(6)

!***************************************************************************
!Beginning of executable session
!***************************************************************************

 if(iexit/=0)then
   if (allocated(vin))        then
     ABI_DEALLOCATE(vin)
   end if
   if (allocated(vout))       then
     ABI_DEALLOCATE(vout)
   end if
   if (allocated(vin_prev))   then
     ABI_DEALLOCATE(vin_prev)
   end if
   if (allocated(vout_prev))  then
     ABI_DEALLOCATE(vout_prev)
   end if
   if (allocated(vinres))  then
     ABI_DEALLOCATE(vinres)
   end if
   if (allocated(vin1))  then
     ABI_DEALLOCATE(vin1)
   end if
   if (allocated(hessin))     then
     ABI_DEALLOCATE(hessin)
   end if
   return
 end if

!write(std_out,*) 'bfgs 01'
!##########################################################
!### 01. Debugging and Verbose

 if(zDEBUG)then
   write(std_out,'(a,3a,35a,42a)') ch10,('-',kk=1,3),&
&   'Debugging and Verbose for pred_bfgs',('-',kk=1,42)
   write(std_out,*) 'ionmov: ',ionmov
   write(std_out,*) 'itime:  ',itime
 end if

!write(std_out,*) 'bfgs 02'
!##########################################################
!### 02. Compute the dimension of vectors (ndim)

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

 if(zDEBUG) write(std_out,*) 'Dimension of vin, vout and hessian (ndim): ',ndim

!write(std_out,*) 'bfgs 03'
!##########################################################
!### 03. Allocate the vectors vin, vout and hessian matrix

!Notice that vin, vout, etc could be allocated
!From a previous dataset with a different ndim
 if(itime==1)then
   if (allocated(vin))        then
     ABI_DEALLOCATE(vin)
   end if
   if (allocated(vout))       then
     ABI_DEALLOCATE(vout)
   end if
   if (allocated(vin_prev))   then
     ABI_DEALLOCATE(vin_prev)
   end if
   if (allocated(vout_prev))  then
     ABI_DEALLOCATE(vout_prev)
   end if
   if (allocated(vinres))  then
     ABI_DEALLOCATE(vinres)
   end if
   if (allocated(vin1))  then
     ABI_DEALLOCATE(vin1)
   end if
   if (allocated(hessin))     then
     ABI_DEALLOCATE(hessin)
   end if
   if(npul>1) then
     ABI_ALLOCATE(vinres,(npul+1,ndim))
     ABI_ALLOCATE(vin1,(npul+1,ndim))
   end if
   ABI_ALLOCATE(vin,(ndim))
   ABI_ALLOCATE(vout,(ndim))
   ABI_ALLOCATE(vin_prev,(ndim))
   ABI_ALLOCATE(vout_prev,(ndim))
   ABI_ALLOCATE(hessin,(ndim,ndim))
 end if

!write(std_out,*) 'bfgs 04'
!##########################################################
!### 04. Obtain the present values from the history

 call hist2var(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)
 do ii=1,3
   rprim(ii,1:3)=rprimd(ii,1:3)/acell(1:3)
 end do

 strten(:)=hist%strten(:,hist%ihist)
 etotal   =hist%etot(hist%ihist)

!Fill the residual with forces (No preconditioning)
!Or the preconditioned forces
 if (ab_mover%goprecon==0)then
   call fcart2fred(hist%fcart(:,:,hist%ihist),residual,rprimd,ab_mover%natom)
 else
   residual(:,:)=forstr%fred(:,:)
 end if

 if(zDEBUG)then
   write (std_out,*) 'residual:'
   do kk=1,ab_mover%natom
     write (std_out,*) residual(:,kk)
   end do
   write (std_out,*) 'strten:'
   write (std_out,*) strten(1:3),ch10,strten(4:6)
   write (std_out,*) 'etotal:'
   write (std_out,*) etotal
 end if

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Save initial values
 if (itime==1)then
   acell0(:)=acell(:)
   rprimd0(:,:)=rprimd(:,:)
   ucvol0=ucvol
 end if

!zDEBUG (UCVOL)
 if(zDEBUG)then
   write(std_out,*) 'Volume of cell (ucvol):',ucvol
 end if

!Get rid of mean force on whole unit cell, but only if no
!generalized constraints are in effect
 residual_corrected(:,:)=residual(:,:)
 if(ab_mover%nconeq==0)then
   do ii=1,3
     if (ii/=3.or.ab_mover%jellslab==0) then
       favg=sum(residual_corrected(ii,:))/dble(ab_mover%natom)
       residual_corrected(ii,:)=residual_corrected(ii,:)-favg
     end if
   end do
 end if

!write(std_out,*) 'bfgs 05'
!##########################################################
!### 05. Fill the vectors vin and vout

!Initialize input vectors : first vin, then vout
!The values of vin from the previous iteration
!should be the same
!if (itime==1)then
 call xfpack_x2vin(acell, acell0, ab_mover%natom, ndim,&
& ab_mover%nsym, ab_mover%optcell, rprim, rprimd0,&
& ab_mover%symrel, ucvol, ucvol0, vin, xred)
!end if

 call xfpack_f2vout(residual_corrected, ab_mover%natom, ndim,&
& ab_mover%optcell, ab_mover%strtarget, strten, ucvol,&
& vout)

!write(std_out,*) 'bfgs 06'
!##########################################################
!### 06. Initialize or update the hessian matrix

!Initialise the Hessian matrix using gmet
 if (itime==1)then

   call hessinit(ab_mover, hessin, gmet, ndim, ucvol)

!  ! Initialize inverse hessian with identity matrix
!  ! in cartesian coordinates, which makes use of metric tensor gmet
!  ! in reduced coordinates.
!  hessin(:,:)=zero
!  do ii=1,ab_mover%natom
!  do kk=1,3
!  do jj=1,3
!  ! Warning : implemented in reduced coordinates
!  if (ab_mover%iatfix(kk,ii)==0 .and.&
!  & ab_mover%iatfix(jj,ii)==0 )then
!  hessin(kk+3*(ii-1),jj+3*(ii-1))=gmet(kk,jj)
!  end if
!  end do
!  end do
!  end do
!  if(ab_mover%optcell/=0)then
!  ! These values might lead to too large changes in some cases
!  diag=ab_mover%strprecon*30.0_dp/ucvol
!  if(ab_mover%optcell==1) diag=diag/three
!  do ii=3*ab_mover%natom+1,ndim
!  hessin(ii,ii)=diag
!  end do
!  end if

   if (ab_mover%restartxf/=0) then

     call xfh_recover_new(ab_xfh,ab_mover,acell,acell0,cycl_main,residual,&
&     hessin,ndim,rprim,rprimd0,strten,ucvol,ucvol0,vin,&
&     vin_prev,vout,vout_prev,xred)

   end if

 end if

 if(itime>1)then
!  Update the hessian matrix, by taking into account the
!  current pair (x,f) and the previous one.
   call hessupdt(hessin,ab_mover%iatfix,ab_mover%natom,ndim,vin,&
&   vin_prev,vout,vout_prev)

 end if

!zDEBUG (vin,vout and hessin before prediction)
 if(zDEBUG)then
   write(std_out,*) 'Vectors vin and vout and inverse of Hessian (hessin) [before]'
   write(std_out,*) 'vin:'
   do ii=1,ndim,3
     if (ii+2<=ndim)then
       write(std_out,*) ii,vin(ii:ii+2)
     else
       write(std_out,*) ii,vin(ii:ndim)
     end if
   end do
   write(std_out,*) 'vout:'
   do ii=1,ndim,3
     if (ii+2<=ndim)then
       write(std_out,*) ii,vout(ii:ii+2)
     else
       write(std_out,*) ii,vout(ii:ndim)
     end if
   end do
   write(std_out,*) 'Inverse Hessian (hessin): ',ndim,'x',ndim
   do kk=1,ndim
     do jj=1,ndim,3
       if (jj+2<=ndim)then
         write(std_out,*) jj,hessin(jj:jj+2,kk)
       else
         write(std_out,*) jj,hessin(jj:ndim,kk)
       end if
     end do
   end do
 end if

!write(std_out,*) 'bfgs 07'
!##########################################################
!### 07. Compute the next values

 if(ionmov==2 .or. itime==1)then

!  Previous cartesian coordinates
   vin_prev(:)=vin(:)

!  New atomic cartesian coordinates are obtained from vin, hessin
!  and vout

   do ii=1,ndim
     vin(:)=vin(:)-hessin(:,ii)*vout(ii)
   end do

!Pulay mixing for vin
   nitpul=0
   if (npul>1) then
     alpha0=1.0_dp
     nitpul=min(itime, npul)
     if (itime>npul) then
       do jj=1,npul-1
         vinres(jj,:)=vinres(jj+1,:)
         vin1(jj,:)=vin1(jj+1,:)
       end do
     end if
     vinres(nitpul,:)=vin(:)-vin_prev(:)
     vin1(nitpul,:)=vin_prev(:)
   end if

   if (nitpul>1) then
     ABI_ALLOCATE(alpha,(nitpul,ndim))
     alpha=zero
     do kk=1,ndim
       ABI_ALLOCATE(amat,(nitpul,nitpul))
       ABI_ALLOCATE(amatinv,(nitpul,nitpul))
       amat=zero;amatinv=zero
       do ii=1,nitpul
         do jj=ii,nitpul
           amat(ii,jj)=vinres(jj,kk)*vinres(ii,kk)
           amat(jj,ii)=amat(ii,jj)
         end do
       end do
       amatinv=amat
       if (abs(vin(kk)-vin_prev(kk))<tol10) then
         alpha(:,kk)=zero
       else
         ABI_ALLOCATE(ipiv,(nitpul))
         ABI_ALLOCATE(rwork,(nitpul))
!          amatinv=1.d5*amatinv
         call dgetrf(nitpul,nitpul,amatinv,nitpul,ipiv,ierr)
         call dgetri(nitpul,amatinv,nitpul,ipiv,rwork,nitpul,ierr)
!          amatinv=1.d5*amatinv
         ABI_DEALLOCATE(ipiv)
         ABI_DEALLOCATE(rwork)
         det=zero
         do ii=1,nitpul
           do jj=1,nitpul
             alpha(ii,kk)=alpha(ii,kk)+amatinv(jj,ii)
             det=det+amatinv(jj,ii)
           end do
         end do
         alpha(:,kk)=alpha(:,kk)/det
       end if
     end do
     ABI_DEALLOCATE(amat)
     ABI_DEALLOCATE(amatinv)
     vin(:)=vin1(nitpul,:)+alpha0*(vin1(nitpul+1,:)-vin1(nitpul,:))
     vin=zero
     do ii=1,nitpul
       vin(:)=vin(:)+ alpha(ii,:)*(vin1(ii,:))
     end do
     ABI_DEALLOCATE(alpha)
   end if


!  Previous atomic forces
   vout_prev(:)=vout(:)

 else
   if(ionmov==3)then
     ihist_prev = abihist_findIndex(hist,-1)
     etotal_prev=hist%etot(ihist_prev)
!    Here the BFGS algorithm, modified to take into account the energy
     call brdene(etotal,etotal_prev,hessin,ndim,vin,vin_prev,vout,vout_prev)
   end if

!  zDEBUG (vin,vout and hessin after prediction)
   if(zDEBUG)then
     write(std_out,*) 'Vectors vin and vout [after prediction]'
     write(std_out,*) 'vin:'
     do ii=1,ndim,3
       if (ii+2<=ndim)then
         write(std_out,*) ii,vin(ii:ii+2)
       else
         write(std_out,*) ii,vin(ii:ndim)
       end if
     end do
     write(std_out,*) 'vout:'
     do ii=1,ndim,3
       if (ii+2<=ndim)then
         write(std_out,*) ii,vout(ii:ii+2)
       else
         write(std_out,*) ii,vout(ii:ndim)
       end if
     end do
   end if


!  FIXME: this should be outside the if clause on ionmov!
!  Implement fixing of atoms : put back old values for fixed
!  components
   do kk=1,ab_mover%natom
     do jj=1,3
!      Warning : implemented in reduced coordinates
       if ( ab_mover%iatfix(jj,kk)==1) then
         vin(jj+(kk-1)*3)=vin_prev(jj+(kk-1)*3)
       end if
     end do
   end do
 end if

!write(std_out,*) 'bfgs 08'
!##########################################################
!### 08. Update the history with the prediction

!Increase indexes
 hist%ihist = abihist_findIndex(hist,+1)

!Transfer vin  to xred, acell and rprim
 call xfpack_vin2x(acell, acell0, ab_mover%natom, ndim,&
& ab_mover%nsym, ab_mover%optcell, rprim, rprimd0,&
& ab_mover%symrel, ucvol, ucvol0,&
& vin, xred)

 if(ab_mover%optcell/=0)then
   call mkrdim(acell,rprim,rprimd)
   call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 end if

!Fill the history with the variables
!xred, acell, rprimd, vel
 call var2hist(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)
 ihist_prev = abihist_findIndex(hist,-1)
 hist%vel(:,:,hist%ihist)=hist%vel(:,:,ihist_prev)

 if(zDEBUG)then
   write (std_out,*) 'residual:'
   do kk=1,ab_mover%natom
     write (std_out,*) residual(:,kk)
   end do
   write (std_out,*) 'strten:'
   write (std_out,*) strten(1:3),ch10,strten(4:6)
   write (std_out,*) 'etotal:'
   write (std_out,*) etotal
 end if

end subroutine pred_bfgs
!!***

!!****f* ABINIT/pred_lbfgs
!! NAME
!! pred_lbfgs
!!
!! FUNCTION
!! Ionmov predictors (22) Limited-memory Broyden-Fletcher-Goldfarb-Shanno
!!
!! IONMOV 22:
!! Given a starting point xred that is a vector of length 3*natom
!! (reduced nuclei coordinates), and unit cell parameters
!! (acell and rprim) the L-Broyden-Fletcher-Goldfarb-Shanno
!! minimization is performed on the total energy function, using
!! its gradient (atomic forces and stress : fred or fcart and
!! stress) as calculated by the routine scfcv. Some atoms can be
!! kept fixed, while the optimization of unit cell parameters is
!! only performed if optcell/=0. The convergence requirement on
!! the atomic forces, dtset%tolmxf,  allows an early exit.
!! Otherwise no more than dtset%ntime steps are performed.
!! Returned quantities are xred, and eventually acell and rprim (new ones!).
!! Could see MinPack on netlib.org
!!
!! INPUTS
!! ab_mover <type(abimover)> : Datatype with all the information needed by the preditor
!! itime  : Index of the present iteration
!! ntime  : Maximal number of iterations
!! ionmov : (22) Specific kind of BFGS
!! zDEBUG : if true print some debugging information
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! hist <type(abihist)> : History of positions,forces acell, rprimd, stresses
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      fcart2fred,hist2var,lbfgs_destroy,lbfgs_init,metric,mkrdim,var2hist
!!      xfh_recover_new,xfpack_f2vout,xfpack_vin2x,xfpack_x2vin
!!
!! SOURCE

subroutine pred_lbfgs(ab_mover,ab_xfh,forstr,hist,ionmov,itime,zDEBUG,iexit)

implicit none

!Arguments ------------------------------------
!scalars
type(abimover),intent(in)       :: ab_mover
type(ab_xfh_type),intent(inout)    :: ab_xfh
type(abihist),intent(inout) :: hist
type(abiforstr),intent(in) :: forstr
integer,intent(in) :: itime
integer,intent(in) :: ionmov
integer,intent(in) :: iexit
logical,intent(in) :: zDEBUG
character(len=500) :: ionmov22_errmsg

!Local variables-------------------------------
!scalars
integer :: info,ihist_prev
integer  :: ndim,cycl_main
integer, parameter :: npul=0
integer  :: ii,jj,kk
real(dp),save :: ucvol0
real(dp) :: ucvol
real(dp) :: etotal
real(dp) :: favg

!arrays
real(dp),allocatable :: diag(:)
real(dp),allocatable,save :: hessin(:,:),vin(:),vin_prev(:)
real(dp),allocatable,save :: vout(:),vout_prev(:)
real(dp),allocatable,save ::vinres(:,:),vin1(:,:)
real(dp),save :: acell0(3) ! Initial acell
real(dp),save :: rprimd0(3,3) ! Initial rprimd
real(dp) :: acell(3)
real(dp) :: rprimd(3,3),rprim(3,3)
real(dp) :: gprimd(3,3)
real(dp) :: gmet(3,3)
real(dp) :: rmet(3,3)
real(dp) :: residual(3,ab_mover%natom),residual_corrected(3,ab_mover%natom)
real(dp) :: xred(3,ab_mover%natom)
real(dp) :: strten(6)

!***************************************************************************
!Beginning of executable session
!***************************************************************************

 if(iexit/=0)then
   call lbfgs_destroy()

   if (allocated(vin))        then
     ABI_DEALLOCATE(vin)
   end if
   if (allocated(vout))       then
     ABI_DEALLOCATE(vout)
   end if
   if (allocated(vin_prev))   then
     ABI_DEALLOCATE(vin_prev)
   end if
   if (allocated(vout_prev))  then
     ABI_DEALLOCATE(vout_prev)
   end if
   if (allocated(vinres))  then
     ABI_DEALLOCATE(vinres)
   end if
   if (allocated(vin1))  then
     ABI_DEALLOCATE(vin1)
   end if
   if (allocated(hessin))     then
     ABI_DEALLOCATE(hessin)
   end if
   return
 end if

!write(std_out,*) 'bfgs 01'
!##########################################################
!### 01. Debugging and Verbose

 if(zDEBUG)then
   write(std_out,'(a,3a,35a,42a)') ch10,('-',kk=1,3),&
&   'Debugging and Verbose for pred_bfgs',('-',kk=1,42)
   write(std_out,*) 'ionmov: ',ionmov
   write(std_out,*) 'itime:  ',itime
 end if

!write(std_out,*) 'bfgs 02'
!##########################################################
!### 02. Compute the dimension of vectors (ndim)

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

 if(zDEBUG) write(std_out,*) 'Dimension of vin, vout and hessian (ndim): ',ndim

!write(std_out,*) 'bfgs 03'
!##########################################################
!### 03. Allocate the vectors vin, vout and hessian matrix

!Notice that vin, vout, etc could be allocated
!From a previous dataset with a different ndim
 if(itime==1)then
   if (allocated(vin))        then
     ABI_DEALLOCATE(vin)
   end if
   if (allocated(vout))       then
     ABI_DEALLOCATE(vout)
   end if
   if (allocated(vin_prev))   then
     ABI_DEALLOCATE(vin_prev)
   end if
   if (allocated(vout_prev))  then
     ABI_DEALLOCATE(vout_prev)
   end if
   if (allocated(vinres))  then
     ABI_DEALLOCATE(vinres)
   end if
   if (allocated(vin1))  then
     ABI_DEALLOCATE(vin1)
   end if
   if (allocated(hessin))     then
     ABI_DEALLOCATE(hessin)
   end if
   if(npul>1) then
     ABI_ALLOCATE(vinres,(npul+1,ndim))
     ABI_ALLOCATE(vin1,(npul+1,ndim))
   end if
   ABI_ALLOCATE(vin,(ndim))
   ABI_ALLOCATE(vout,(ndim))
   ABI_ALLOCATE(vin_prev,(ndim))
   ABI_ALLOCATE(vout_prev,(ndim))
   ABI_ALLOCATE(hessin,(ndim,ndim))
 end if

!write(std_out,*) 'bfgs 04'
!##########################################################
!### 04. Obtain the present values from the history

 call hist2var(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)
 do ii=1,3
   rprim(ii,1:3)=rprimd(ii,1:3)/acell(1:3)
 end do

 strten(:)=hist%strten(:,hist%ihist)
 etotal   =hist%etot(hist%ihist)

!Fill the residual with forces (No preconditioning)
!Or the preconditioned forces
 if (ab_mover%goprecon==0)then
   call fcart2fred(hist%fcart(:,:,hist%ihist),residual,rprimd,ab_mover%natom)
 else
   residual(:,:)= forstr%fred(:,:)
 end if

 if(zDEBUG)then
   write (std_out,*) 'residual:'
   do kk=1,ab_mover%natom
     write (std_out,*) residual(:,kk)
   end do
   write (std_out,*) 'strten:'
   write (std_out,*) strten(1:3),ch10,strten(4:6)
   write (std_out,*) 'etotal:'
   write (std_out,*) etotal
 end if

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Save initial values
 if (itime==1)then
   acell0(:)=acell(:)
   rprimd0(:,:)=rprimd(:,:)
   ucvol0=ucvol
 end if

!zDEBUG (UCVOL)
 if(zDEBUG)then
   write(std_out,*) 'Volume of cell (ucvol):',ucvol
 end if

!Get rid of mean force on whole unit cell, but only if no
!generalized constraints are in effect
 residual_corrected(:,:)=residual(:,:)
 if(ab_mover%nconeq==0)then
   do ii=1,3
     if (ii/=3.or.ab_mover%jellslab==0) then
       favg=sum(residual_corrected(ii,:))/dble(ab_mover%natom)
       residual_corrected(ii,:)=residual_corrected(ii,:)-favg
     end if
   end do
 end if

!write(std_out,*) 'bfgs 05'
!##########################################################
!### 05. Fill the vectors vin and vout

!Initialize input vectors : first vin, then vout
!The values of vin from the previous iteration
!should be the same
!if (itime==1)then
 call xfpack_x2vin(acell, acell0, ab_mover%natom, ndim,&
& ab_mover%nsym, ab_mover%optcell, rprim, rprimd0,&
& ab_mover%symrel, ucvol, ucvol0, vin, xred)
!end if

 call xfpack_f2vout(residual_corrected, ab_mover%natom, ndim,&
& ab_mover%optcell, ab_mover%strtarget, strten, ucvol,&
& vout)

!write(std_out,*) 'bfgs 06'
!##########################################################
!### 06. Initialize or update the hessian matrix

!Initialise the Hessian matrix using gmet
 if (itime==1)then

   ABI_ALLOCATE(diag,(ndim))
   do ii=1,3*ab_mover%natom
!      diag(ii) = 1.00_dp / rprimd(MODULO(ii-1,3)+1,MODULO(ii-1,3)+1)**2
     diag(ii) = gmet(MODULO(ii-1,3)+1,MODULO(ii-1,3)+1)
   end do
   if(ab_mover%optcell/=0)then
!     These values might lead to too large changes in some cases ...
     do ii=3*ab_mover%natom+1,ndim
       diag(ii) = ab_mover%strprecon*30.0_dp/ucvol
       if(ab_mover%optcell==1) diag(ii) = diag(ii) / three
     end do
   end if

   call lbfgs_init(ndim,5,diag)
   ABI_DEALLOCATE(diag)

   if (ab_mover%restartxf/=0) then

     call xfh_recover_new(ab_xfh,ab_mover,acell,acell0,cycl_main,residual,&
&     hessin,ndim,rprim,rprimd0,strten,ucvol,ucvol0,vin,&
&     vin_prev,vout,vout_prev,xred)

   end if

 end if

!zDEBUG (vin,vout and hessin before prediction)
 if(zDEBUG)then
   write(std_out,*) 'Vectors vin and vout and inverse of Hessian (hessin) [before]'
   write(std_out,*) 'vin:'
   do ii=1,ndim,3
     if (ii+2<=ndim)then
       write(std_out,*) ii,vin(ii:ii+2)
     else
       write(std_out,*) ii,vin(ii:ndim)
     end if
   end do
   write(std_out,*) 'vout:'
   do ii=1,ndim,3
     if (ii+2<=ndim)then
       write(std_out,*) ii,vout(ii:ii+2)
     else
       write(std_out,*) ii,vout(ii:ndim)
     end if
   end do
   write(std_out,*) 'Inverse Hessian (hessin): ',ndim,'x',ndim
   do kk=1,ndim
     do jj=1,ndim,3
       if (jj+2<=ndim)then
         write(std_out,*) jj,hessin(jj:jj+2,kk)
       else
         write(std_out,*) jj,hessin(jj:ndim,kk)
       end if
     end do
   end do
 end if

!write(std_out,*) 'bfgs 07'
!##########################################################
!### 07. Compute the next values

 vin_prev(:) = vin
 vout_prev(:) = vout
 info = lbfgs_execute(vin,etotal,vout)

 if (info /= -1) then
   write (ionmov22_errmsg, '(a,i0,3a)') &
    'Lbfgs routine failed. Returned value: ', info,ch10, &
    'Restart your calculation from last step or try a different ionmov'
   MSG_ERROR_CLASS(ionmov22_errmsg, "Ionmov22Error")
 end if

!zDEBUG (vin,vout after prediction)
 if(zDEBUG)then
   write(std_out,*) 'Vectors vin and vout [after prediction]'
   write(std_out,*) 'vin_prev:'
   do ii=1,ndim,3
     if (ii+2<=ndim)then
       write(std_out,*) ii,vin_prev(ii:ii+2)
     else
       write(std_out,*) ii,vin_prev(ii:ndim)
     end if
   end do
   write(std_out,*) 'vin:'
   do ii=1,ndim,3
     if (ii+2<=ndim)then
       write(std_out,*) ii,vin(ii:ii+2)
     else
       write(std_out,*) ii,vin(ii:ndim)
     end if
   end do
   write(std_out,*) 'vout:'
   do ii=1,ndim,3
     if (ii+2<=ndim)then
       write(std_out,*) ii,vout(ii:ii+2)
     else
       write(std_out,*) ii,vout(ii:ndim)
     end if
   end do
 end if


!Implement fixing of atoms : put back old values for fixed
!components
 do kk=1,ab_mover%natom
   do jj=1,3
!    Warning : implemented in reduced coordinates
     if ( ab_mover%iatfix(jj,kk)==1) then
       vin(jj+(kk-1)*3)=vin_prev(jj+(kk-1)*3)
     end if
   end do
 end do


!write(std_out,*) 'bfgs 08'
!##########################################################
!### 08. Update the history with the prediction

!Increase indexes
 hist%ihist = abihist_findIndex(hist,+1)

!Transfer vin  to xred, acell and rprim
 call xfpack_vin2x(acell, acell0, ab_mover%natom, ndim,&
& ab_mover%nsym, ab_mover%optcell, rprim, rprimd0,&
& ab_mover%symrel, ucvol, ucvol0,&
& vin, xred)

 if(ab_mover%optcell/=0)then
   call mkrdim(acell,rprim,rprimd)
   call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 end if

!Fill the history with the variables
!xcart, xred, acell, rprimd
 call var2hist(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)
 ihist_prev = abihist_findIndex(hist,-1)
 hist%vel(:,:,hist%ihist)=hist%vel(:,:,ihist_prev)

 if(zDEBUG)then
   write (std_out,*) 'residual:'
   do kk=1,ab_mover%natom
     write (std_out,*) residual(:,kk)
   end do
   write (std_out,*) 'strten:'
   write (std_out,*) strten(1:3),ch10,strten(4:6)
   write (std_out,*) 'etotal:'
   write (std_out,*) etotal
 end if

end subroutine pred_lbfgs
!!***

end module m_pred_bfgs
!!***
