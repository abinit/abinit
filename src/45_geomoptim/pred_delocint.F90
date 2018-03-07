!{\src2tex{textfont=tt}}
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
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR, JCC, SE)
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
!! ionmov : (10 or 11) Specific kind of BFGS
!! zDEBUG : if true print some debugging information
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! hist <type(abihist)> : History of positions,forces
!!                               acell, rprimd, stresses
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pred_delocint(ab_mover,ab_xfh,forstr,hist,ionmov,itime,zDEBUG,iexit)

 use defs_basis
 use m_profiling_abi
 use m_abimover
 use m_abihist

 use m_bfgs, only : hessinit, hessupdt, brdene

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pred_delocint'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
 use interfaces_45_geomoptim, except_this_one => pred_delocint
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(abimover),intent(inout)       :: ab_mover
 type(ab_xfh_type),intent(inout)    :: ab_xfh
 type(abihist),intent(inout) :: hist
 type(abiforstr),intent(in) :: forstr
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
 integer,save :: icenter,irshift ! DELOCINT indexes
 integer,save :: nshell,ndeloc ! DELOCINT number of
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
!  DELOCINT
!  Allocate all the variables of deloc,
!  Needed for compilers such as Pathscale that fails if
!  unallocated variables are passed as arguments
!  TODO:
!  all of the following should go into an init routine in m_abimover

   nshell=3
   ab_mover%deloc%nrshift=(2*nshell+1)**3
   icenter = nshell*(2*nshell+1)**2 + nshell*(2*nshell+1) + nshell + 1

   ABI_ALLOCATE(ab_mover%deloc%rshift,(3,ab_mover%deloc%nrshift))
   irshift=0
   do ii=-nshell,nshell
     do jj=-nshell,nshell
       do kk=-nshell,nshell
         irshift=irshift+1
         ab_mover%deloc%rshift(:,irshift) = (/dble(ii),dble(jj),dble(kk)/)
       end do
     end do
   end do

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
   call make_prim_internals(ab_mover%deloc,icenter,ab_mover%natom,&
&   ab_mover%ntypat,rprimd,ab_mover%typat,xcart,ab_mover%znucl)

   ABI_ALLOCATE(prim_int,(ab_mover%deloc%ninternal))

   if(DEBUG)then
     write (message,'(a,i6)') 'Number of primitive internal coordinates (ninternal): ',ab_mover%deloc%ninternal
     call wrtout(std_out,  message,'COLL')
   end if

   if (allocated(u_matrix))  then
     ABI_DEALLOCATE(u_matrix)
   end if
   ABI_ALLOCATE(u_matrix,(ab_mover%deloc%ninternal,ndeloc))

   call calc_prim_int(ab_mover%deloc,ab_mover%natom,rprimd,xcart,prim_int)

   if(DEBUG)then
     write (message,'(a)') 'Primitive internal coordinate values:'
     call wrtout(std_out,  message,'COLL')
     write (message,'(a)') ' Bonds:'
     call wrtout(std_out,  message,'COLL')
     do ii = 1, ab_mover%deloc%nbond
       write (message,'(i6,E20.10)') ii, prim_int(ii)
       call wrtout(std_out,  message,'COLL')
     end do

     write (message,'(a)') ' Angles:'
     call wrtout(std_out,  message,'COLL')
     do ii = ab_mover%deloc%nbond+1, ab_mover%deloc%nbond+ab_mover%deloc%nang
       write (message,'(i6,2(E20.10,2x))') ii, prim_int(ii), prim_int(ii)/pi*180.0_dp
       call wrtout(std_out,  message,'COLL')
     end do

     write (message,'(a)') ' Dihedrals:'
     call wrtout(std_out,  message,'COLL')
     do ii = ab_mover%deloc%nbond+ab_mover%deloc%nang+1, ab_mover%deloc%nbond+ab_mover%deloc%nang+ab_mover%deloc%ndihed
       write (message,'(i6,2(E20.10,2x))') ii, prim_int(ii), prim_int(ii)/pi*180.0_dp
       call wrtout(std_out,  message,'COLL')
     end do

     write (message,'(a)') ' Cartesian auxiliary coordinates for constraints:'
     call wrtout(std_out,  message,'COLL')
     do ii = ab_mover%deloc%nbond+ab_mover%deloc%nang+ab_mover%deloc%ndihed+1, ab_mover%deloc%ninternal
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

 ABI_ALLOCATE(prim_int,(ab_mover%deloc%ninternal))

!write(std_out,*) 'delocint 06'
!##########################################################
!### 06. Compute delocalized coordinates and forces

!xcart ---> deloc_int

!Convert positions to delocalized coordinates for next step
 call xcart2deloc(ab_mover%deloc,ab_mover%natom,rprimd,xcart,&
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
&     vout,vout_prev,xred,ab_mover%deloc,deloc_int,deloc_force,bt_inv_matrix,gprimd,prim_int,&
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
 call deloc2xcart(ab_mover%deloc,ab_mover%natom,rprimd,xcart,&
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
