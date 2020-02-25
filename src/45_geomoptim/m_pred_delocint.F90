!!****m* ABINIT/m_pred_delocint
!! NAME
!! m_pred_delocint
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (MVer, DCA, XG, GMR, JCC, SE)
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

module m_pred_delocint

 use defs_basis
 use m_errors
 use m_abicore
 use m_abimover
 use m_abihist
 use m_xfpack
 use m_linalg_interfaces

 use m_geometry,   only : fcart2fred, xcart2xred, xred2xcart, metric, acrossb
 use m_bfgs,       only : hessinit, hessupdt, brdene
 use m_results_gs, only : results_gs_type

 implicit none

 private
!!***

 public :: pred_delocint

contains
!!***

!!****f* ABINIT/pred_delocint
!! NAME
!! pred_delocint
!!
!! FUNCTION
!! Ionmov predictors (10) BFGS with delocalized internal coordinates
!!
!! IONMOV 10:
!! Given a starting point xred that is a vector of length 3*(natom-1)
!! (reduced nuclei coordinates),
!! and unit cell parameters (acell and rprimd) the
!! Broyden-Fletcher-Goldfarb-Shanno minimization is performed on the
!! total energy function, using its gradient (atomic forces and stresses)
!  as calculated by the routine scfcv. Some atoms can be kept fixed,
!! while the optimization of unit cell
!! parameters is only performed if optcell/=0.
!! The convergence requirement on
!! the atomic forces, 'tolmxf',  allows an early exit.
!! Otherwise no more than 'ntime' steps are performed.
!! Returned quantities are xred, and eventually acell and rprimd (new ones!).
!! Could see Numerical Recipes (Fortran), 1986, page 307.
!!
!!  Implements the delocalized internal coordinate scheme
!!  of Andzelm et al. in CPL .335. 321 (2001) \
!!  and Baker et al. JCP .105. 192 (1996)
!!
!!    B matrix is derivative of delocalized internals wrt cartesian coordinates
!!    U matrix is eigenvectors of G = B*B^{T}
!!    S matrix is eigenvectors of F = B^{T}B
!!
!! INPUTS
!! ab_mover <type(abimover)> : Datatype with all the information needed by the preditor
!! itime  : Index of the present iteration
!! ntime  : Maximal number of iterations
!! ionmov : (10 or 11) Specific kind of BFGS
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
!!      brdene,calc_prim_int,deloc2xcart,fcart2fred,fred2fdeloc,hessupdt
!!      hist2var,make_prim_internals,metric,var2hist,wrtout,xcart2deloc
!!      xcart2xred,xfh_recover_deloc,xfpack_f2vout,xfpack_vin2x,xfpack_x2vin
!!      xred2xcart
!!
!! SOURCE

subroutine pred_delocint(ab_mover,ab_xfh,deloc,forstr,hist,ionmov,itime,zDEBUG,iexit)

 implicit none

!Arguments ------------------------------------
!scalars
 type(abimover),intent(in)       :: ab_mover
 type(ab_xfh_type),intent(inout)    :: ab_xfh
 type(abihist),intent(inout) :: hist
 type(abiforstr),intent(in) :: forstr
 type(delocint),intent(inout) :: deloc
 integer,intent(in) :: itime
 integer,intent(in) :: ionmov
 integer,intent(in) :: iexit
 logical,intent(in) :: zDEBUG

!Local variables-------------------------------
!scalars
 integer  :: ndim,cycl_main
 integer  :: ihist_prev,ii,jj,kk
 real(dp),save :: ucvol0
 real(dp) :: ucvol
 real(dp) :: etotal,etotal_prev
 logical  :: DEBUG=.TRUE.
 !integer,save :: icenter,irshift ! DELOCINT indexes
 integer,save :: ndeloc ! DELOCINT number of
 character(len=500) :: message

!arrays
 real(dp),allocatable,save :: hessin(:,:),vin(:),vin_prev(:)
 real(dp),allocatable,save :: vout(:),vout_prev(:)
 real(dp),save :: acell0(3) ! Initial acell
 real(dp),save :: rprimd0(3,3) ! Initial rprimd
 real(dp),allocatable :: prim_int(:)
 real(dp),allocatable,save :: u_matrix(:,:) ! DELOCINT this may need to be added to type inside ab_mover
 real(dp) :: acell(3)
 real(dp) :: rprimd(3,3)
 real(dp) :: gprimd(3,3)
 real(dp) :: gmet(3,3)
 real(dp) :: rmet(3,3)
 real(dp) :: residual(3,ab_mover%natom)
!real(dp) :: residual_corrected(3,ab_mover%natom)
 real(dp) :: xred(3,ab_mover%natom),xcart(3,ab_mover%natom)
 real(dp) :: strten(6)
 real(dp) :: deloc_force(3*(ab_mover%natom-1))
 real(dp) :: deloc_int(3*(ab_mover%natom-1))
 real(dp) :: bt_inv_matrix(3*(ab_mover%natom-1),3*ab_mover%natom)

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
   if (allocated(hessin))     then
     ABI_DEALLOCATE(hessin)
   end if
   if (allocated(u_matrix))     then
     ABI_DEALLOCATE(u_matrix)
   end if
   return
 end if

!write(std_out,*) 'delocint 01'
!##########################################################
!### 01. Debugging and Verbose

 if(DEBUG)then
   write(std_out,'(a,3a,38a,39a)') ch10,('-',kk=1,3),&
&   'Debugging and Verbose for pred_deloint',('-',kk=1,39)
   write(std_out,*) 'ionmov: ',ionmov
   write(std_out,*) 'itime:  ',itime
 end if

!write(std_out,*) 'delocint 02'
!##########################################################
!### 02. Compute the dimension of vectors (ndim)

!With internal we have 1 coordinate less
 ndeloc = 3*(ab_mover%natom-1)
 ndim=ndeloc
 deloc_int(:)=zero
 deloc_force(:)=zero
 if(ab_mover%optcell==1 .or.&
& ab_mover%optcell==4 .or.&
& ab_mover%optcell==5 .or.&
& ab_mover%optcell==6) ndim=ndim+1
 if(ab_mover%optcell==2 .or.&
& ab_mover%optcell==3) ndim=ndim+6
 if(ab_mover%optcell==7 .or.&
& ab_mover%optcell==8 .or.&
& ab_mover%optcell==9) ndim=ndim+3

 if(DEBUG) write(std_out,*) 'Dimension of vin, vout and hessian (ndim): ',ndim

!write(std_out,*) 'delocint 03'
!##########################################################
!### 03. Allocate the vectors vin, vout and hessian matrix

!Notice thqt vin, vout, etc could be allocated
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
   if (allocated(hessin))     then
     ABI_DEALLOCATE(hessin)
   end if
   ABI_ALLOCATE(vin,(ndim))
   ABI_ALLOCATE(vout,(ndim))
   ABI_ALLOCATE(vin_prev,(ndim))
   ABI_ALLOCATE(vout_prev,(ndim))
   ABI_ALLOCATE(hessin,(ndim,ndim))

 end if


!write(std_out,*) 'delocint 04'
!##########################################################
!### 04. Obtain the present values from the history

 call hist2var(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)
 call xred2xcart(ab_mover%natom,rprimd,xcart,xred)

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

!DEBUG (UCVOL)
 if(DEBUG)then
   write(std_out,*) 'Volume of cell (ucvol):',ucvol
 end if

!Get rid of mean force on whole unit cell, but only if no
!generalized constraints are in effect
!  residual_corrected(:,:)=residual(:,:)
!  if(ab_mover%nconeq==0)then
!    do ii=1,3
!      if (ii/=3.or.ab_mover%jellslab==0) then
!        favg=sum(residual_corrected(ii,:))/dble(ab_mover%natom)
!        residual_corrected(ii,:)=residual_corrected(ii,:)-favg
!      end if
!    end do
!  end if

!write(std_out,*) 'delocint 05'
!##########################################################
!### 05. Compute internals for first time

 if (itime==1)then
   call make_prim_internals(deloc,ab_mover%natom,&
&   ab_mover%ntypat,rprimd,ab_mover%typat,xcart,ab_mover%znucl)

   ABI_ALLOCATE(prim_int,(deloc%ninternal))

   if(DEBUG)then
     write (message,'(a,i6)') 'Number of primitive internal coordinates (ninternal): ',deloc%ninternal
     call wrtout(std_out,  message,'COLL')
   end if

   if (allocated(u_matrix))  then
     ABI_DEALLOCATE(u_matrix)
   end if
   ABI_ALLOCATE(u_matrix,(deloc%ninternal,ndeloc))

   call calc_prim_int(deloc,ab_mover%natom,rprimd,xcart,prim_int)

   if(DEBUG)then
     write (message,'(a)') 'Primitive internal coordinate values:'
     call wrtout(std_out,  message,'COLL')
     write (message,'(a)') ' Bonds:'
     call wrtout(std_out,  message,'COLL')
     do ii = 1, deloc%nbond
       write (message,'(i6,E20.10)') ii, prim_int(ii)
       call wrtout(std_out,  message,'COLL')
     end do

     write (message,'(a)') ' Angles:'
     call wrtout(std_out,  message,'COLL')
     do ii = deloc%nbond+1, deloc%nbond+deloc%nang
       write (message,'(i6,2(E20.10,2x))') ii, prim_int(ii), prim_int(ii)/pi*180.0_dp
       call wrtout(std_out,  message,'COLL')
     end do

     write (message,'(a)') ' Dihedrals:'
     call wrtout(std_out,  message,'COLL')
     do ii = deloc%nbond+deloc%nang+1, deloc%nbond+deloc%nang+deloc%ndihed
       write (message,'(i6,2(E20.10,2x))') ii, prim_int(ii), prim_int(ii)/pi*180.0_dp
       call wrtout(std_out,  message,'COLL')
     end do

     write (message,'(a)') ' Cartesian auxiliary coordinates for constraints:'
     call wrtout(std_out,  message,'COLL')
     do ii = deloc%nbond+deloc%nang+deloc%ndihed+1, deloc%ninternal
       write (message,'(i6,E20.10)') ii, prim_int(ii)
       call wrtout(std_out,  message,'COLL')
     end do
   end if

   ABI_DEALLOCATE(prim_int)

!  equal weight on all internal coordinates as a starting point.
   u_matrix(:,:) = one / dble (ndeloc)

!  Zero the arrays before first use
   deloc_force(:) = zero

 end if

 ABI_ALLOCATE(prim_int,(deloc%ninternal))

!write(std_out,*) 'delocint 06'
!##########################################################
!### 06. Compute delocalized coordinates and forces

!xcart ---> deloc_int

!Convert positions to delocalized coordinates for next step
 call xcart2deloc(deloc,ab_mover%natom,rprimd,xcart,&
& bt_inv_matrix,u_matrix,deloc_int,prim_int)

!fred ---> deloc_force

!Convert forces to delocalized coordinates for next step
 call fred2fdeloc(bt_inv_matrix,deloc_force,residual,ab_mover%natom,gprimd)

!write(std_out,*) 'delocint 07'
!##########################################################
!### 07. Fill the vectors vin and vout

!DEBUG deloc_int and deloc_force before pack
 if(DEBUG)then
   write (std_out,*) 'Delocalized internals and forces (ndeloc):',ndeloc
   write(std_out,*) 'deloc_int'
   do ii=1,ndeloc,3
     if (ii+2<=ndeloc)then
       write(std_out,*) ii,deloc_int(ii:ii+2)
     else
       write(std_out,*) ii,deloc_int(ii:ndeloc)
     end if
   end do
   write(std_out,*) 'deloc_force'
   do ii=1,ndeloc,3
     if (ii+2<=ndeloc)then
       write(std_out,*) ii,deloc_force(ii:ii+2)
     else
       write(std_out,*) ii,deloc_force(ii:ndeloc)
     end if
   end do
 end if

!DELOCINT
!Instead of fred_corrected we use deloc_force
!Instead of xred e use deloc_int
!
!Initialize input vectors : first vin, then vout
!The values of vin from the previous iteration
!should be the same
 call xfpack_x2vin(acell, acell0, ab_mover%natom-1, ndim,&
& ab_mover%nsym, ab_mover%optcell, rprimd, rprimd0,&
& ab_mover%symrel, ucvol, ucvol0, vin, deloc_int)
!end if

 call xfpack_f2vout(deloc_force, ab_mover%natom-1, ndim,&
& ab_mover%optcell, ab_mover%strtarget, strten, ucvol, &
& vout)

!write(std_out,*) 'delocint 08'
!##########################################################
!### 08. Initialize or update the hessian matrix

!Initialise the Hessian matrix using gmet
 if (itime==1)then

!  Initialise the Hessian matrix with ab_mover%userrc.
!  this has become unusable because it imposes ndim >= 3 natom
!  ident = 3x3 identity matrix
!  call hessinit(ab_mover, hessin, gmet, ndim, ucvol)
   hessin = zero
   do ii=1, ndim
     hessin (ii,ii) = one
   end do

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

     call xfh_recover_deloc(ab_xfh,ab_mover,acell,acell0,cycl_main,&
&     residual,hessin,ndim,rprimd,rprimd0,strten,ucvol,ucvol0,vin,vin_prev,&
&     vout,vout_prev,xred,deloc,deloc_int,deloc_force,bt_inv_matrix,gprimd,prim_int,&
&     u_matrix)

   end if

 end if

 ABI_DEALLOCATE(prim_int)

 if(itime>1)then
!  Update the hessian matrix, by taking into account the
!  current pair (x,f) and the previous one.
   call hessupdt(hessin,ab_mover%iatfix,ab_mover%natom,ndim,vin,&
&   vin_prev,vout,vout_prev)

 end if

!DEBUG (vin,vout and hessin before prediction)
 if(DEBUG)then
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

!write(std_out,*) 'delocint 09'
!##########################################################
!### 09. Compute the next values

 if(ionmov==10 .or. itime==1)then

!  Previous cartesian coordinates
   vin_prev(:)=vin(:)

!  New atomic cartesian coordinates are obtained from vin, hessin
!  and vout
   do ii=1,ndim
     vin(:)=vin(:)-hessin(:,ii)*vout(ii)
   end do
!  Previous atomic forces
   vout_prev(:)=vout(:)

 else
   if(ionmov==11)then
     ihist_prev = abihist_findIndex(hist,-1)
     etotal_prev=hist%etot(ihist_prev)
!    Here the BFGS algorithm, modified to take into account the
!    energy
     call brdene(etotal,etotal_prev,hessin,&
&     ndim,vin,vin_prev,vout,vout_prev)

   end if

!  DEBUG (vin,vout and hessin after prediction)
   if(DEBUG)then
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

!write(std_out,*) 'delocint 10'
!##########################################################
!### 10. Convert from delocalized to xcart and xred

!Transfer vin  to deloc_int, acell and rprimd
 call xfpack_vin2x(acell, acell0, ab_mover%natom-1, ndim,&
& ab_mover%nsym, ab_mover%optcell, rprimd, rprimd0,&
& ab_mover%symrel, ucvol, ucvol0,&
& vin, deloc_int)

 if(DEBUG)then
   write (std_out,*) 'Delocalized internals (deloc_int) [after prediction]:'
   write(std_out,*) 'deloc_int:'
   do ii=1,ndeloc,3
     if (ii+2<=ndeloc)then
       write(std_out,*) ii,deloc_int(ii:ii+2)
     else
       write(std_out,*) ii,deloc_int(ii:ndeloc)
     end if
   end do
   write(std_out,*) 'BT Inverse Matrix:'
   do ii=1,3*ab_mover%natom
     write(std_out,*) bt_inv_matrix(:,ii)
   end do
   write (std_out,*) 'xcart (before deloc2xcart):'
   do ii=1,ab_mover%natom
     write (std_out,*) xcart(:,ii)
   end do
 end if

!this routine contains an iterative scheme to find xcart
!from the non-linear relations between deloc and xcart
!SIGNIFICANTLY DIFFERENT FROM xcart2deloc
 call deloc2xcart(deloc,ab_mover%natom,rprimd,xcart,&
& deloc_int,bt_inv_matrix,u_matrix)

!Convert new xcart (cartesian) to xred (reduced coordinates)
 call xcart2xred(ab_mover%natom,rprimd,xcart,xred)


!write(std_out,*) 'delocint 11'
!##########################################################
!### 11. Update the history with the prediction

!Increase indexes
 hist%ihist = abihist_findIndex(hist,+1)

 if(ab_mover%optcell/=0)then
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

end subroutine pred_delocint
!!***

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

subroutine deloc2xcart(deloc,natom,rprimd,xcart,deloc_int,btinv,u_matrix)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom
 type(delocint),intent(in) :: deloc
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

!!****f* ABINIT/fred2fdeloc
!! NAME
!! fred2fdeloc
!!
!! FUNCTION
!!  calculate delocalized forces from reduced coordinate ones
!!
!! INPUTS
!! btinv(3*(natom-1),3*natom)= inverse transpose of B matrix (see delocint)
!! natom = number of atoms
!! gprimd(3,3)=dimensional translations in reciprocal space (bohr-1)
!!
!! OUTPUT
!! deloc_force(3*(natom-1))=delocalized forces from reduced coordinate ones
!! fred(3,natom)=delocalized forces in reduced coordinates
!!
!! PARENTS
!!      pred_delocint,xfh_recover_deloc
!!
!! CHILDREN
!!      dgemm,dgemv,wrtout
!!
!! SOURCE

subroutine fred2fdeloc(btinv,deloc_force,fred,natom,gprimd)

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom
!arrays
 real(dp),intent(in) :: btinv(3*(natom-1),3*natom),gprimd(3,3),fred(3,natom)
 real(dp),intent(out) :: deloc_force(3*(natom-1))

!Local variables-------------------------------
 integer :: ii
!arrays
 real(dp) :: fcart(3,natom)
 character(len=500) :: message

! ******************************************************************

!make cartesian forces

 call dgemm('N','N',3,natom,3,one,&
& gprimd,3,fred,3,zero,fcart,3)

!turn cartesian to delocalized forces
 call dgemv('N',3*(natom-1),3*natom,one,&
& btinv,3*(natom-1),fcart,1,zero,deloc_force,1)

 write (message,'(a)') 'fred2fdeloc : deloc_force = '
 call wrtout(std_out,message,'COLL')

 do ii = 1, 3*(natom-1)
   write (message,'(I6,E16.6)') ii, deloc_force(ii)
   call wrtout(std_out,message,'COLL')
 end do

end subroutine fred2fdeloc
!!***

!!****f* ABINIT/calc_b_matrix
!! NAME
!! calc_b_matrix
!!
!! FUNCTION
!!  calculate values of derivatives of internal coordinates as a function of
!!  cartesian ones =  B matrix
!!
!! INPUTS
!! angs= number of angles
!! bonds(2,2,nbond)=for a bond between iatom and jatom
!!              bonds(1,1,nbond) = iatom
!!              bonds(2,1,nbond) = icenter
!!              bonds(1,2,nbond) = jatom
!!              bonds(2,2,nbond) = irshift
!! carts(2,ncart)= index of total primitive internal, and atom (carts(2,:))
!! dihedrals(2,4,ndihed)=indexes to characterize dihedrals
!! nang(2,3,nang)=indexes to characterize angles
!! nbond=number of bonds
!! ncart=number of auxiliary cartesian atom coordinates (used for constraints)
!! ndihed= number of dihedrals
!! ninternal=nbond+nang+ndihed+ncart: number of internal coordinates
!! nrshift= dimension of rshift
!! rprimd(3,3)=dimensional real space primitive translations (bohr)
!! rshift(3,nrshift)=shift in xred that must be done to find all neighbors of
!!                   a given atom within a given number of neighboring shells
!! xcart(3,natom)=cartesian coordinates of atoms (bohr)
!!
!! OUTPUT
!! b_matrix(ninternal,3*natom)=matrix of derivatives of internal coordinates
!!   wrt cartesians
!!
!! PARENTS
!!      xcart2deloc
!!
!! CHILDREN
!!      acrossb
!!
!! SOURCE

subroutine calc_b_matrix(deloc,natom,rprimd,xcart,b_matrix)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom
 type(delocint),intent(in) :: deloc

!arrays
 real(dp),intent(in) :: rprimd(3,3),xcart(3,natom)
 real(dp),intent(out) :: b_matrix(deloc%ninternal,3*natom)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,i4,iang,ibond,icart,idihed,iprim,s1,s2,s3,s4
!arrays
 real(dp) :: bb(3),r1(3),r2(3),r3(3),r4(3)

! *************************************************************************

 iprim=0
 b_matrix(:,:) = zero

 do ibond=1,deloc%nbond
   i1 = deloc%bonds(1,1,ibond)
   s1 = deloc%bonds(2,1,ibond)
   r1(:) = xcart(:,i1)+deloc%rshift(1,s1)*rprimd(:,1)&
&   +deloc%rshift(2,s1)*rprimd(:,2)&
&   +deloc%rshift(3,s1)*rprimd(:,3)
   i2 = deloc%bonds(1,2,ibond)
   s2 = deloc%bonds(2,2,ibond)
   r2(:) = xcart(:,i2)+deloc%rshift(1,s2)*rprimd(:,1)&
&   +deloc%rshift(2,s2)*rprimd(:,2)&
&   +deloc%rshift(3,s2)*rprimd(:,3)
   iprim=iprim+1
   call dbond_length_d1(r1,r2,bb)
   b_matrix(iprim,3*(i1-1)+1:3*i1) = b_matrix(iprim,3*(i1-1)+1:3*i1) + bb(:)
   call dbond_length_d1(r2,r1,bb)
   b_matrix(iprim,3*(i2-1)+1:3*i2) = b_matrix(iprim,3*(i2-1)+1:3*i2) + bb(:)
 end do

!second: angle values (ang)
 do iang=1,deloc%nang
   i1 = deloc%angs(1,1,iang)
   s1 = deloc%angs(2,1,iang)
   r1(:) = xcart(:,i1)+deloc%rshift(1,s1)*rprimd(:,1)&
&   +deloc%rshift(2,s1)*rprimd(:,2)&
&   +deloc%rshift(3,s1)*rprimd(:,3)
   i2 = deloc%angs(1,2,iang)
   s2 = deloc%angs(2,2,iang)
   r2(:) = xcart(:,i2)+deloc%rshift(1,s2)*rprimd(:,1)&
&   +deloc%rshift(2,s2)*rprimd(:,2)&
&   +deloc%rshift(3,s2)*rprimd(:,3)
   i3 = deloc%angs(1,3,iang)
   s3 = deloc%angs(2,3,iang)
   r3(:) = xcart(:,i3)+deloc%rshift(1,s3)*rprimd(:,1)&
&   +deloc%rshift(2,s3)*rprimd(:,2)&
&   +deloc%rshift(3,s3)*rprimd(:,3)
   iprim=iprim+1
   call dang_d1(r1,r2,r3,bb)
   b_matrix(iprim,3*(i1-1)+1:3*i1) = b_matrix(iprim,3*(i1-1)+1:3*i1) + bb(:)
   call dang_d2(r1,r2,r3,bb)
   b_matrix(iprim,3*(i2-1)+1:3*i2) = b_matrix(iprim,3*(i2-1)+1:3*i2) + bb(:)
   call dang_d1(r3,r2,r1,bb)
   b_matrix(iprim,3*(i3-1)+1:3*i3) = b_matrix(iprim,3*(i3-1)+1:3*i3) + bb(:)
 end do

!third: dihedral values
 do idihed=1,deloc%ndihed
   i1 = deloc%dihedrals(1,1,idihed)
   s1 = deloc%dihedrals(2,1,idihed)
   r1(:) = xcart(:,i1)+deloc%rshift(1,s1)*rprimd(:,1)&
&   +deloc%rshift(2,s1)*rprimd(:,2)&
&   +deloc%rshift(3,s1)*rprimd(:,3)
   i2 = deloc%dihedrals(1,2,idihed)
   s2 = deloc%dihedrals(2,2,idihed)
   r2(:) = xcart(:,i2)+deloc%rshift(1,s2)*rprimd(:,1)&
&   +deloc%rshift(2,s2)*rprimd(:,2)&
&   +deloc%rshift(3,s2)*rprimd(:,3)
   i3 = deloc%dihedrals(1,3,idihed)
   s3 = deloc%dihedrals(2,3,idihed)
   r3(:) = xcart(:,i3)+deloc%rshift(1,s3)*rprimd(:,1)&
&   +deloc%rshift(2,s3)*rprimd(:,2)&
&   +deloc%rshift(3,s3)*rprimd(:,3)
   i4 = deloc%dihedrals(1,4,idihed)
   s4 = deloc%dihedrals(2,4,idihed)
   r4(:) = xcart(:,i4)+deloc%rshift(1,s4)*rprimd(:,1)&
&   +deloc%rshift(2,s4)*rprimd(:,2)&
&   +deloc%rshift(3,s4)*rprimd(:,3)
!  write(std_out,*) 'dihed ',idihed
!  write(std_out,*) r1
!  write(std_out,*) r2
!  write(std_out,*) r3
!  write(std_out,*) r4

   iprim=iprim+1
   call ddihedral_d1(r1,r2,r3,r4,bb)
   b_matrix(iprim,3*(i1-1)+1:3*i1) = b_matrix(iprim,3*(i1-1)+1:3*i1) + bb(:)
   call ddihedral_d2(r1,r2,r3,r4,bb)
   b_matrix(iprim,3*(i2-1)+1:3*i2) = b_matrix(iprim,3*(i2-1)+1:3*i2) + bb(:)
   call ddihedral_d2(r4,r3,r2,r1,bb)
   b_matrix(iprim,3*(i3-1)+1:3*i3) = b_matrix(iprim,3*(i3-1)+1:3*i3) + bb(:)
   call ddihedral_d1(r4,r3,r2,r1,bb)
   b_matrix(iprim,3*(i4-1)+1:3*i4) = b_matrix(iprim,3*(i4-1)+1:3*i4) + bb(:)
 end do

 do icart=1,deloc%ncart
   iprim=iprim+1
   b_matrix(iprim,3*(deloc%carts(2,icart)-1)+deloc%carts(1,icart)) = &
&   b_matrix(iprim,3*(deloc%carts(2,icart)-1)+deloc%carts(1,icart)) + one
 end do

!DEBUG
! write (200,*) 'calc_b_matrix : b_matrix = '
! do iprim=1,deloc%ninternal
!   do i1=1, 3*natom
!     write (200,'(E16.6,2x)',ADVANCE='NO') b_matrix(iprim,i1)
!   end do
!   write (200,*)
! end do
!ENDDEBUG

end subroutine calc_b_matrix
!!***

!!****f* ABINIT/dbond_length_d1
!! NAME
!! dbond_length_d1
!!
!! FUNCTION
!!
!! PARENTS
!!      calc_b_matrix
!!
!! CHILDREN
!!      acrossb
!!
!! SOURCE
!!

subroutine dbond_length_d1(r1,r2,bb)

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: r1(3),r2(3)
 real(dp),intent(out) :: bb(3)

!Local variables ------------------------------
!arrays
 real(dp) :: rpt(3)

!************************************************************************
 rpt(:) = r1(:)-r2(:)
 bb(:) = rpt(:)/bond_length(r1,r2)

end subroutine dbond_length_d1
!!***


!!****f* ABINIT/dang_d1
!! NAME
!! dang_d1
!!
!! FUNCTION
!!
!! PARENTS
!!      calc_b_matrix
!!
!! CHILDREN
!!      acrossb
!!
!! SOURCE
!!

subroutine dang_d1(r1,r2,r3,bb)

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: r1(3),r2(3),r3(3)
 real(dp),intent(out) :: bb(3)

!Local variables ------------------------------
!scalars
 real(dp) :: cos_ang,n1,n1232,n2,tmp
!arrays
 real(dp) :: cp1232(3),rpt(3),rpt12(3),rpt32(3)

!************************************************************************
 n1=bond_length(r1,r2)
 n2=bond_length(r3,r2)

 rpt12(:) = r1(:)-r2(:)
 rpt32(:) = r3(:)-r2(:)

 cos_ang = (rpt12(1)*rpt32(1)+rpt12(2)*rpt32(2)+rpt12(3)*rpt32(3))/n1/n2
 if (cos_ang > one - epsilon(one)*two) then
   cos_ang = one
 else if(cos_ang < -one + epsilon(one)*two) then
   cos_ang = -one
 end if

 rpt(:) = rpt32(:)/n1/n2 - rpt12(:)*cos_ang/n1/n1

 tmp = sqrt(one-cos_ang**2)
 bb(:) = zero
 if (tmp > epsilon(one)) then
   bb(:) = rpt(:) * (-one)/tmp
 end if

!TEST: version from MOLECULAR VIBRATIONS EB Wilson
 call acrossb(rpt12,rpt32,cp1232)
 n1232 = sqrt(cp1232(1)**2+cp1232(2)**2+cp1232(3)**2)
 rpt(:) = (cos_ang*rpt12(:)*n2/n1 - rpt32(:))/n1232
 if (abs(bb(1)-rpt(1))+abs(bb(2)-rpt(2))+abs(bb(3)-rpt(3)) > tol10) then
   write(std_out,*) 'Compare bb ang 1 : '
   write(std_out,*) bb(:), rpt(:), bb(:)-rpt(:)
 end if
 bb(:) = rpt(:)

end subroutine dang_d1
!!***


!!****f* ABINIT/dang_d2
!! NAME
!! dang_d2
!!
!! FUNCTION
!!
!! PARENTS
!!      calc_b_matrix
!!
!! CHILDREN
!!      acrossb
!!
!! SOURCE
!!

subroutine dang_d2(r1,r2,r3,bb)

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: r1(3),r2(3),r3(3)
 real(dp),intent(out) :: bb(3)

!Local variables ------------------------------
!scalars
 real(dp) :: cos_ang,n1,n1232,n2,tmp
!arrays
 real(dp) :: cp1232(3),rpt(3),rpt12(3),rpt32(3)

!************************************************************************
 n1=bond_length(r1,r2)
 n2=bond_length(r3,r2)

 rpt12(:) = r1(:)-r2(:)
 rpt32(:) = r3(:)-r2(:)

 cos_ang = (rpt12(1)*rpt32(1)+rpt12(2)*rpt32(2)+rpt12(3)*rpt32(3))/n1/n2
 if (cos_ang > one - epsilon(one)*two) then
   cos_ang = one
 else if(cos_ang < -one + epsilon(one)*two) then
   cos_ang = -one
 end if

 rpt(:) = -rpt32(:)/n1/n2 - rpt12(:)/n1/n2 &
& + rpt12(:)*cos_ang/n1/n1 + rpt32(:)*cos_ang/n2/n2

 tmp = sqrt(one-cos_ang**2)
 bb(:) = zero
 if (tmp > tol12) then
   bb(:) = rpt(:) * (-one)/tmp
 end if

!TEST: version from MOLECULAR VIBRATIONS EB Wilson
 call acrossb(rpt12,rpt32,cp1232)
 n1232 = sqrt(cp1232(1)**2+cp1232(2)**2+cp1232(3)**2)
 rpt(:) = ((n1-n2*cos_ang)*rpt12(:)/n1 + (n2-n1*cos_ang)*rpt32(:)/n2) / n1232
 if (abs(bb(1)-rpt(1))+abs(bb(2)-rpt(2))+abs(bb(3)-rpt(3))  > tol10) then
   write(std_out,*) 'Compare bb ang 2 : '
   write(std_out,*) bb(:), rpt(:), bb(:)-rpt(:)
 end if
 bb(:) = rpt(:)

end subroutine dang_d2
!!***

!!****f* ABINIT/ddihedral_d1
!! NAME
!! ddihedral_d1
!!
!! FUNCTION
!!
!! PARENTS
!!      calc_b_matrix
!!
!! CHILDREN
!!      acrossb
!!
!! SOURCE
!!

subroutine ddihedral_d1(r1,r2,r3,r4,bb)

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: r1(3),r2(3),r3(3),r4(3)
 real(dp),intent(out) :: bb(3)

!Local variables ------------------------------------
!scalars
 real(dp) :: cos_dihedral,dih_sign,n1,n2,n23,sin_dihedral,tmp
!arrays
 real(dp) :: cp1232(3),cp32_1232(3),cp32_3432(3),cp3432(3),cpcp(3),rpt(3)
 real(dp) :: rpt12(3),rpt32(3),rpt34(3)

!******************************************************************
 rpt12(:) = r1(:)-r2(:)
 rpt32(:) = r3(:)-r2(:)
 rpt34(:) = r3(:)-r4(:)

 call acrossb(rpt12,rpt32,cp1232)
 call acrossb(rpt34,rpt32,cp3432)

!DEBUG
!write(std_out,*) ' cos_dihedral : cp1232 = ', cp1232
!write(std_out,*) ' cos_dihedral : cp3432 = ', cp3432
!ENDDEBUG

 n1 = sqrt(cp1232(1)**2+cp1232(2)**2+cp1232(3)**2)
 n2 = sqrt(cp3432(1)**2+cp3432(2)**2+cp3432(3)**2)

 cos_dihedral = (cp1232(1)*cp3432(1)+cp1232(2)*cp3432(2)+cp1232(3)*cp3432(3))/n1/n2
 if (cos_dihedral > one - epsilon(one)*two) then
   cos_dihedral = one
 else if(cos_dihedral < -one + epsilon(one)*two) then
   cos_dihedral = -one
 end if
!we use complementary of standard angle, so
!cos_dihedral = -cos_dihedral

 call acrossb(cp1232,cp3432,cpcp)
 cpcp(:) = cpcp(:)/n1/n2
!we use complementary of standard angle, but sin is invariant
 sin_dihedral = -(cpcp(1)*rpt32(1)+cpcp(2)*rpt32(2)+cpcp(3)*rpt32(3))&
& /sqrt(rpt32(1)**2+rpt32(2)**2+rpt32(3)**2)
 dih_sign = one
 if (sin_dihedral < -epsilon(one)) then
   dih_sign = -one
 end if

!DEBUG
!write(std_out,'(a,3E16.6)') 'ddihedral_d1 : cos abs(sin) dih_sign= ',&
!&    cos_dihedral,sin_dihedral,dih_sign
!ENDDEBUG

!ddihedral_d1 = dih_sign* acos(cos_dihedral)
 call acrossb(rpt32,cp1232,cp32_1232)
 call acrossb(rpt32,cp3432,cp32_3432)

 rpt(:) = cp32_3432(:)/n1/n2 - cp32_1232(:)/n1/n1 * cos_dihedral
 bb(:) = zero

!DEBUG
!write(std_out,*) 'ddihedral_d1 cp1232 cp3432 = ',cp1232,cp3432,rpt32
!write(std_out,*) 'ddihedral_d1 cp32_1232 cp32_3432 = ',cp32_1232,cp32_3432,cos_dihedral,n1,n2
!write(std_out,*) 'ddihedral_d1 rpt = ',rpt
!ENDDEBUG

 tmp = sqrt(one-cos_dihedral**2)
 if (tmp > tol12) then
!  we use complementary of standard angle, so cosine in acos has - sign,
!  and it appears for the derivative
   bb(:) = -dih_sign * rpt(:) * (-one) / tmp
 else
   bb(:) = dih_sign * cp32_3432(:) / n1 / n2 / &
&   sqrt(cp32_3432(1)**2+cp32_3432(2)**2+cp32_3432(3)**2)
 end if

!TEST: version from MOLECULAR VIBRATIONS EB Wilson

 n23 = sqrt(rpt32(1)*rpt32(1)+rpt32(2)*rpt32(2)+rpt32(3)*rpt32(3))
 rpt(:) = cp1232(:)*n23/n1/n1
!if (abs(bb(1)-rpt(1))+abs(bb(2)-rpt(2))+abs(bb(3)-rpt(3))  > tol10) then
!write(std_out,*) 'Compare bb1 : '
!write(std_out,*) bb(:), rpt(:), bb(:)-rpt(:)
!end if
 bb(:) = rpt(:)

end subroutine ddihedral_d1
!!***

!!****f* ABINIT/ddihedral_d2
!! NAME
!! ddihedral_d2
!!
!! FUNCTION
!!
!! PARENTS
!!      calc_b_matrix
!!
!! CHILDREN
!!      acrossb
!!
!! SOURCE
!!

subroutine ddihedral_d2(r1,r2,r3,r4,bb)

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: r1(3),r2(3),r3(3),r4(3)
 real(dp),intent(out) :: bb(3)

!Local variables
!scalars
 real(dp) :: cos_dihedral,dih_sign,n1,n2,n23,sin_dihedral,sp1232,sp3432,tmp
!arrays
 real(dp) :: cp1232(3),cp1232_12(3),cp1232_34(3),cp32_1232(3),cp32_3432(3)
 real(dp) :: cp3432(3),cp3432_12(3),cp3432_34(3),cpcp(3),rpt(3),rpt12(3)
 real(dp) :: rpt32(3),rpt34(3)

! *************************************************************************
 rpt12(:) = r1(:)-r2(:)
 rpt32(:) = r3(:)-r2(:)
 rpt34(:) = r3(:)-r4(:)

 call acrossb(rpt12,rpt32,cp1232)
 call acrossb(rpt34,rpt32,cp3432)

!DEBUG
!write(std_out,*) ' cos_dihedral : cp1232 = ', cp1232
!write(std_out,*) ' cos_dihedral : cp3432 = ', cp3432
!ENDDEBUG

 n1 = sqrt(cp1232(1)**2+cp1232(2)**2+cp1232(3)**2)
 n2 = sqrt(cp3432(1)**2+cp3432(2)**2+cp3432(3)**2)

 cos_dihedral = (cp1232(1)*cp3432(1)+cp1232(2)*cp3432(2)+cp1232(3)*cp3432(3))/n1/n2
 if (cos_dihedral > one - epsilon(one)*two) then
   cos_dihedral = one
 else if(cos_dihedral < -one + epsilon(one)*two) then
   cos_dihedral = -one
 end if
!we use complementary of standard angle, so
!cos_dihedral = -cos_dihedral

 call acrossb(cp1232,cp3432,cpcp)
 cpcp(:) = cpcp(:)/n1/n2
!we use complementary of standard angle, but sin is invariant
 sin_dihedral = -(cpcp(1)*rpt32(1)+cpcp(2)*rpt32(2)+cpcp(3)*rpt32(3))&
& /sqrt(rpt32(1)**2+rpt32(2)**2+rpt32(3)**2)
 dih_sign = one
 if (sin_dihedral <  -tol12) then
   dih_sign = -one
 end if

!DEBUG
!write(std_out,'(a,3E16.6)') 'ddihedral_d2 : cos abs(sin) dih_sign= ',&
!&    cos_dihedral,sin_dihedral,dih_sign
!ENDDEBUG

!ddihedral_d2 = dih_sign* acos(cos_dihedral)
 call acrossb(rpt32,cp3432,cp32_3432)
 call acrossb(cp3432,rpt12,cp3432_12)
 call acrossb(cp1232,rpt34,cp1232_34)

 call acrossb(rpt32,cp1232,cp32_1232)
 call acrossb(cp1232,rpt12,cp1232_12)
 call acrossb(cp3432,rpt34,cp3432_34)

 rpt(:) = -(cp32_3432(:) + cp3432_12(:) + cp1232_34(:))/n1/n2 &
& +cos_dihedral*(cp32_1232(:)/n1/n1 + cp1232_12(:)/n1/n1 + cp3432_34(:)/n2/n2)
 bb(:) = zero
 tmp = sqrt(one-cos_dihedral**2)
 if (tmp > tol12) then
!  we use complementary of standard angle, so cosine in acos has - sign,
!  and it appears for derivative
   bb(:) = -dih_sign * rpt(:) * (-one) / tmp
 else
   bb(:) = dih_sign * cos_dihedral * &
&   ( cp32_1232(:)/n1/n1/sqrt(cp32_1232(1)**2+cp32_1232(2)**2+cp32_1232(3)**2) &
&   +cp1232_12(:)/n1/n1/sqrt(cp1232_12(1)**2+cp1232_12(2)**2+cp1232_12(3)**2) &
&   +cp3432_34(:)/n2/n2/sqrt(cp3432_34(1)**2+cp3432_34(2)**2+cp3432_34(3)**2) )
 end if

!TEST: version from MOLECULAR VIBRATIONS EB Wilson p. 61
 n23 = sqrt(rpt32(1)*rpt32(1)+rpt32(2)*rpt32(2)+rpt32(3)*rpt32(3))
 sp1232 = rpt12(1)*rpt32(1)+rpt12(2)*rpt32(2)+rpt12(3)*rpt32(3)
 sp3432 = rpt34(1)*rpt32(1)+rpt34(2)*rpt32(2)+rpt34(3)*rpt32(3)

 rpt(:) = -cp1232(:)*(n23-sp1232/n23)/n1/n1 - cp3432(:)*sp3432/n23/n2/n2
!if (abs(bb(1)-rpt(1))+abs(bb(2)-rpt(2))+abs(bb(3)-rpt(3))  > tol10) then
!write(std_out,*) 'Compare bb2 : '
!write(std_out,*) bb(:), rpt(:), bb(:)-rpt(:)
!write(std_out,*) -cp1232(:)*(n23-sp1232/n23)/n1/n1, -cp3432(:)*sp3432/n23/n2/n2
!end if
 bb(:) = rpt(:)

end subroutine ddihedral_d2
!!***

!!****f* ABINIT/xcart2deloc
!! NAME
!! xcart2deloc
!!
!! FUNCTION
!!  Calculate values of delocalized coordinates as a function of
!!  cartesian ones. First primitive internals, then B matrix,
!!  then F, then U then delocalized internals.
!!
!! INPUTS
!! deloc <type(delocint)>=Important variables for pred_delocint
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
!! natom = Number of atoms
!! rprimd(3,3) = Dimensional real space primitive translations
!!               (bohr)
!! xcart(3,natom) = Cartesian coordinates of atoms (bohr)
!!
!! OUTPUT
!! bt_inv_matrix(3*(natom-1),3*natom) = Inverse of B^{T} matrix
!! deloc_int(3*(natom-1)) = Delocalized internal coordinates
!! prim_int(ninternal) = Primitive internal coordinates
!!
!! SIDE EFFECTS
!! u_matrix(ninternal,3*(natom-1)) = Eigenvectors of BB^T matrix
!!
!! NOTES
!!
!! PARENTS
!!      deloc2xcart,pred_delocint,xfh_recover_deloc
!!
!! CHILDREN
!!
!! SOURCE

subroutine xcart2deloc(deloc,natom,rprimd,xcart,bt_inv_matrix,u_matrix,deloc_int,prim_int)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom
 type(delocint),intent(in) :: deloc
!arrays
 real(dp),intent(in) :: rprimd(3,3),xcart(3,natom)
 real(dp),intent(inout) :: u_matrix(deloc%ninternal,3*(natom-1))
 real(dp),intent(out) :: bt_inv_matrix(3*(natom-1),3*natom)
 real(dp),intent(out) :: deloc_int(3*(natom-1))
 real(dp),intent(out) :: prim_int(deloc%ninternal)

!Local variables-------------------------------
!scalars
integer :: ii
logical :: DEBUG=.FALSE.
!arrays
 real(dp) :: b_matrix(deloc%ninternal,3*natom)

! ******************************************************************

 call calc_prim_int(deloc,natom,rprimd,xcart,prim_int)
 if (DEBUG)then
   write(std_out,*) 'Primitive Internals'
   do ii=1,deloc%ninternal
     write(std_out,*) prim_int(ii)
   end do
 end if

 call calc_b_matrix(deloc,natom,rprimd,xcart,b_matrix)
 if (DEBUG)then
   write(std_out,*) 'B Matrix'
   do ii=1,deloc%ninternal
     write(std_out,*) b_matrix(:,ii)
   end do
 end if

 call calc_btinv_matrix(b_matrix,natom,deloc%ninternal,&
& bt_inv_matrix,u_matrix)
 if (DEBUG)then
   write(std_out,*) 'BT Inverse Matrix'
   do ii=1,3*natom
     write(std_out,*) bt_inv_matrix(:,ii)
   end do
 end if

!calculate value of delocalized internals

 call dgemv('T',deloc%ninternal,3*(natom-1),one,&
& u_matrix,deloc%ninternal,prim_int,1,zero,deloc_int,1)

end subroutine xcart2deloc
!!***


!!****f* ABINIT/calc_btinv_matrix
!! NAME
!! calc_btinv_matrix
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      xcart2deloc
!!
!! NOTES
!!   bt_inv_matrix is inverse transpose of the delocalized
!!    coordinate B matrix. b_matrix is the primitive internal B matrix
!!
!! CHILDREN
!!
!! SOURCE

 subroutine calc_btinv_matrix(b_matrix,natom,ninternal,bt_inv_matrix,u_matrix)

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: ninternal,natom
 real(dp),intent(in) :: b_matrix(ninternal,3*natom)
 real(dp),intent(out) :: bt_inv_matrix(3*(natom-1),3*natom)
 real(dp),intent(inout) :: u_matrix(ninternal,3*(natom-1))

!Local variables ------------------------------------
!scalars
 integer :: ii,info,lwork
!arrays
 real(dp) :: f_eigs(3*natom),f_matrix(3*natom,3*natom)
 real(dp) :: s_matrix(3*natom,3*natom)
 real(dp) :: s_red(3*natom,3*(natom-1))
 real(dp) :: u_matrix_old(ninternal,3*(natom-1))
 real(dp),allocatable :: work(:)

!******************************************************************

!f matrix = B^{T} B
 call dgemm('T','N',3*natom,3*natom,ninternal,one,&
& b_matrix,ninternal,b_matrix,ninternal,zero,f_matrix,3*natom)

 lwork = max(1,3*3*natom-1)
 ABI_ALLOCATE(work,(lwork))
 s_matrix(:,:) = f_matrix(:,:)

 call dsyev('V','L',3*natom,s_matrix,3*natom,f_eigs,work,lwork,info)

 ABI_DEALLOCATE(work)

 if (abs(f_eigs(1)) + abs(f_eigs(2)) + abs(f_eigs(3)) > tol10 ) then
   write(std_out,*) 'Error: 3 lowest eigenvalues are not zero'
   write(std_out,*) '  internal coordinates do NOT span the full degrees of freedom !'
   write(std_out,'(6E16.6)') f_eigs
   MSG_ERROR("Aborting now")
 end if
 if ( abs(f_eigs(4)) < tol10 ) then
   write(std_out,*) 'Error: fourth eigenvalue is zero'
   write(std_out,*) '  internal coordinates do NOT span the full degrees of freedom !'
   write(std_out,'(6E16.6)') f_eigs
   MSG_ERROR("Aborting now")
 end if

!calculate U matrix from U = B * S_red * lambda^{-1/2}
 do ii=1,3*(natom-1)
   s_red(:,ii) = s_matrix(:,ii+3)/sqrt(f_eigs(ii+3))
 end do

 u_matrix_old(:,:) = u_matrix(:,:)

 call dgemm('N','N',ninternal,3*(natom-1),3*natom,one,&
& b_matrix,ninternal,s_red,3*natom,zero,u_matrix,ninternal)


!align eigenvectors, to preserve a form of continuity in convergences
!!!! eigenvalues are no longer in increasing order!!! but only s_red is reordered
!so that btinv is correct.
 call align_u_matrices(natom,ninternal,u_matrix,u_matrix_old,s_matrix,f_eigs)

!calculate B_deloc^{-1} matrix for transformation of forces to deloc coord.
!(B^{T}_deloc)^{-1} = (B_deloc B^{T}_deloc)^{-1} B_deloc = lambda^{-3/2} S^{T} F
!= ( S lambda^{3/2} )^{T} F

!! DEFINITION
!! real(dp),intent(out) :: bt_inv_matrix(3*(natom-1),3*natom)

!even better: B_deloc^{-1} = lambda^{-1/2} S^{T}
 do ii=1,3*(natom-1)
!  s_red(:,ii) = s_matrix(:,ii+3)*sqrt(f_eigs(ii+3))
   bt_inv_matrix(ii,:) = s_matrix(:,ii+3)/sqrt(f_eigs(ii+3))
 end do

end subroutine calc_btinv_matrix
!!***

!!****f* ABINIT/align_u_matrices
!! NAME
!! align_u_matrices
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      xcart2deloc
!!
!! CHILDREN
!!
!! SOURCE

 subroutine align_u_matrices(natom,ninternal,u_matrix,u_matrix_old,s_matrix,f_eigs)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ninternal,natom
!arrays
 real(dp),intent(in) :: u_matrix_old(ninternal,3*(natom-1))
 real(dp),intent(inout) :: f_eigs(3*natom)
 real(dp),intent(inout) :: s_matrix(3*natom,3*natom)
 real(dp),intent(inout) :: u_matrix(ninternal,3*(natom-1))

!Local variables ------------------------------
!scalars
 integer :: ii,iint1,imax
 real(dp) :: ss
!arrays
 integer :: eigv_flag(3*(natom-1)),eigv_ind(3*(natom-1))
 real(dp) :: tmps(3*natom,3*natom)
 real(dp) :: tmpu(ninternal,3*(natom-1))
 real(dp) :: tmpf(3*natom)

!******************************************************************

 eigv_flag(:) = 0
 eigv_ind(:) = 0

!just permit a change in sign
 do iint1=1,3*(natom-1)
   ss = zero
   do ii=1,ninternal
     ss = ss + u_matrix_old(ii,iint1)*u_matrix(ii,iint1)
   end do
   if (ss < -tol12) then
     imax = -iint1
   else
     imax = iint1
   end if
   eigv_ind(iint1) = imax
   eigv_flag(abs(imax)) = 1
 end do

 tmpu(:,:) = u_matrix
 tmps(:,:) = s_matrix
 tmpf(:) = f_eigs
!exchange eigenvectors...
 do iint1=1,3*(natom-1)
   ss = one
   if (eigv_ind(iint1) < 0) ss = -one

   imax = abs(eigv_ind(iint1))

   tmpu(:,imax) = ss*u_matrix(:,iint1)

   tmps(:,imax+3) = ss*s_matrix(:,iint1+3)

   tmpf(imax+3) = f_eigs(iint1+3)
 end do

 u_matrix(:,:) = tmpu(:,:)
 s_matrix(:,:) = tmps(:,:)
 f_eigs(:) = tmpf(:)

end subroutine align_u_matrices
!!***

!!****f* ABINIT/xfh_recover_deloc
!! NAME
!! xfh_recover_deloc
!!
!! FUNCTION
!! Update the contents of the history xfhist taking values
!! from xred, acell, rprim, fred_corrected and strten
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      pred_delocint
!!
!! CHILDREN
!!      fred2fdeloc,hessupdt,xcart2deloc,xfpack_f2vout,xfpack_x2vin,xred2xcart
!!
!! SOURCE

subroutine xfh_recover_deloc(ab_xfh,ab_mover,acell,acell0,cycl_main,&
& fred,hessin,ndim,rprim,rprimd0,strten,ucvol,ucvol0,vin,vin_prev,&
& vout,vout_prev,xred,deloc,deloc_int,deloc_force,btinv,gprimd,prim_int,&
& u_matrix)

implicit none

!Arguments ------------------------------------
!scalars

integer,intent(in) :: ndim
integer,intent(out) :: cycl_main
real(dp),intent(inout) :: ucvol,ucvol0
type(ab_xfh_type),intent(inout) :: ab_xfh
type(abimover),intent(in) :: ab_mover
! DELOCINT specials
type(delocint),intent(in) :: deloc

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

end module m_pred_delocint
!!***
