!!****m* ABINIT/m_bfgs
!! NAME
!!  m_bfgs
!!
!! FUNCTION
!!  This module provides several routines for the application of a
!!  Broyden-Fletcher-Goldfarb-Shanno (BFGS) minimization algorithm.
!!
!! COPYRIGHT
!! Copyright (C) 2012-2020 ABINIT group (XG,JCC)
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

MODULE m_bfgs

 use defs_basis
 use m_abicore
 use m_errors
 use m_abimover

 use m_io_tools,       only : open_file
 use m_numeric_tools,  only : findmin

 implicit none

 private

!public procedures
 public :: hessinit   ! Initialize Hessian matrix
 public :: hessupdt   ! Update the hessian matrix
 public :: brdene
!!***

CONTAINS  !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_bfgs/hessinit
!! NAME
!! hessinit
!!
!! FUNCTION
!! Initiliase an Hessian matrix, either from disk or using init_matrix.
!! The size ndim must be greater or equal than 3 * ab_mover%natom.
!!
!! INPUTS
!!  fnameabi_hes=filename for Hessian matrix
!!  ab_mover = the input variables relevant for moving ions
!!  init_matrix(3,3)=matrix used for each atom (if iatfix = 0) for initialisation.
!!  ndim=size of the hessian and vectors
!!  ucvol=volume of the box (used when ab_mover%optcell is not null).
!!
!! OUTPUT
!!  hessin(ndim,ndim)=hessian matrix, initialised at output.
!!
!! PARENTS
!!      m_pred_bfgs,m_pred_diisrelax
!!
!! CHILDREN
!!      findmin
!!
!! SOURCE

subroutine hessinit(ab_mover, hessin, init_matrix, ndim, ucvol)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndim
 real(dp),intent(in) :: ucvol
 type(abimover),intent(in) :: ab_mover
!arrays
 real(dp),intent(in) :: init_matrix(3,3)
 real(dp),intent(out) :: hessin(ndim,ndim)

!Local variables-------------------------------
!scalars
 integer :: hess_ok,iatom,idim,idir1,idir2,ii,ios,jj,ndim0,temp_unit
 real(dp) :: diag
 logical :: ex
 character(len=500) :: message

! *********************************************************************

!Initialization of the inverse Hessian to the unit matrix
!Much better choices are possible--this simply corresponds to
!taking the first minimization step as the negative of the
!gradient, with the full length of the gradient vector as
!the step size.  Any spring type model would probably be a better starting guess.

 if (ndim < 3 * ab_mover%natom) then
   write(message, '(a,a,a)' )&
&   'the size of the given hessian matrix is too small.', ch10, &
&   'This is an internal error, contact ABINIT developers.'
   ABI_ERROR(message)
 end if

!Special arrangement: if input hessian file exists, read data from there
 inquire (file=ab_mover%fnameabi_hes,iostat=ios,exist=ex)
 hess_ok=0

 if (ex) then
   ! Read inverse hessian data from file; format is
   if (open_file(ab_mover%fnameabi_hes,message,newunit=temp_unit,form='formatted',status='old') /= 0) then
     ABI_ERROR(message)
   end if
   read (temp_unit,*)
   read (temp_unit,*) ndim0
   if (ndim0/=ndim) then
!    Cannot read data because data file natom disagrees with current job
     write(message,'(5a,i10,a,i10,2a)')&
&     'Tried to read inverse hessian from file',trim(ab_mover%fnameabi_hes),' but',ch10,&
&     'ndim of that file =',ndim0,' , is not equal to input ndim =',ndim,ch10,&
&     ' => initialize inverse hessian with identity matrix.'
     ABI_WARNING(message)
     close(unit=temp_unit)
   else
     ! Read inverse hessian
     do jj=1,ndim
       read (temp_unit,*)
       read (temp_unit,*) (hessin(ii,jj),ii=1,ndim)
     end do
     close (unit=temp_unit)
     write(message,*)' Inverse hessian has been input from input hessian file',trim(ab_mover%fnameabi_hes)
     call wrtout(std_out,message,'COLL')
     hess_ok=1
   end if
 end if

!If hessin was not read, initialize inverse hessian with identity matrix
!in cartesian coordinates, which makes use of metric tensor gmet in reduced coordinates.
 if(hess_ok==0)then
   hessin(:,:)=zero
   do iatom=1,ab_mover%natom
     do idir1=1,3
       do idir2=1,3
!        Warning : implemented in reduced coordinates
         if ( ab_mover%iatfix(idir1,iatom) ==0 .and. ab_mover%iatfix(idir2,iatom) ==0 )then
           hessin(idir1+3*(iatom-1),idir2+3*(iatom-1))=init_matrix(idir1,idir2)
         end if
       end do
     end do
   end do
   if(ab_mover%optcell/=0)then
!    These values might lead to too large changes in some cases ...
     diag=ab_mover%strprecon*30.0_dp/ucvol
     if(ab_mover%optcell==1)diag=diag/three
     do idim=3*ab_mover%natom+1,ndim
       hessin(idim,idim)=diag
     end do
   end if
   call wrtout(std_out,'Inverse hessian has been initialized.','COLL')
 end if

end subroutine hessinit
!!***

!----------------------------------------------------------------------

!!****f* m_bfgs/hessupdt
!! NAME
!! hessupdt
!!
!! FUNCTION
!! Update of the hessian matrix according to the Broyden formula.
!! Could see Numerical Recipes (Fortran), 1986, page 307.
!!
!! INPUTS
!!  iatfix(3,natom)=1 for each atom fixed along specified direction, else 0
!!  natom=number of atoms in unit cell
!!  ndim=size of the hessian and vectors
!!  nimage= -- optional, default=1 --
!!         Number of images of the system described in
!!         vin, vin_prev, vout, vout_prev
!!  vin(ndim)=new input vector
!!  vin_prev(ndim)=previous input vector
!!  vout(ndim)=new output vector
!!  vout_prev(ndim)=previous output vector
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  hessin(ndim,ndim)=hessian matrix, updated at output.
!!
!! PARENTS
!!      m_mep,m_pred_bfgs,m_pred_delocint,m_pred_diisrelax,m_xfpack
!!
!! CHILDREN
!!      findmin
!!
!! SOURCE

subroutine hessupdt(hessin,iatfix,natom,ndim,vin,vin_prev,vout,vout_prev, &
&                   nimage) ! optional argument

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ndim
 integer,intent(in),optional :: nimage
!arrays
 integer,intent(in) :: iatfix(3,natom)
 real(dp),intent(in) :: vin(ndim),vin_prev(ndim),vout(ndim),vout_prev(ndim)
 real(dp),intent(inout) :: hessin(ndim,ndim)

!Local variables-------------------------------
!scalars
 integer :: iatom,idir,ii,jj,nimage_
 real(dp) :: den1,den2,den3
 !character(len=500) :: msg
!arrays
 real(dp) :: bfgs(ndim),din(ndim),dout(ndim),hdelta(ndim)

!***************************************************************************

 nimage_=1;if (present(nimage)) nimage_=nimage

!write(ab_out,*) 'VECTOR INPUT (vin)'
!do ii=1,ndim,3
!if (ii+2<=ndim)then
!write(ab_out,*) ii,vin(ii:ii+2)
!else
!write(ab_out,*) ii,vin(ii:ndim)
!end if
!end do
!write(ab_out,*) 'VECTOR OUTPUT (vout)'
!do ii=1,ndim,3
!if (ii+2<=ndim)then
!write(ab_out,*) ii,vout(ii:ii+2)
!else
!write(ab_out,*) ii,vout(ii:ndim)
!end if
!end do

!write(ab_out,*) 'VECTOR INPUT (vin_prev)'
!do ii=1,ndim,3
!if (ii+2<=ndim)then
!write(ab_out,*) ii,vin(ii:ii+2)
!else
!write(ab_out,*) ii,vin(ii:ndim)
!end if
!end do
!write(ab_out,*) 'VECTOR OUTPUT (vout_prev)'
!do ii=1,ndim,3
!if (ii+2<=ndim)then
!write(ab_out,*) ii,vout(ii:ii+2)
!else
!write(ab_out,*) ii,vout(ii:ndim)
!end if
!end do

 if (mod(ndim,nimage_)/=0) then
   ABI_BUG('nimage must be a dividor of ndim !')
 end if

!Difference between new and previous vectors
 din(:) =vin(:) -vin_prev(:)
 dout(:)=vout(:)-vout_prev(:)

!Implement fixing of atoms; must discard the change of forces on fixed atoms
 do ii=1,nimage_
   jj=3*natom*(ii-1)
   do iatom=1,natom
     do idir=1,3
       if (iatfix(idir,iatom)==1) dout(idir+jj)=zero
     end do
     jj=jj+3
   end do
 end do

!Compute approximate inverse Hessian times delta fcart
!hdelta=hessin*deltaf
 hdelta(:)=zero
 do ii=1,ndim
   hdelta(:)=hdelta(:)+hessin(:,ii)*dout(ii)
 end do

!Calculation of dot products for the denominators
 den1=zero ; den2=zero
 do ii=1,ndim
   den1=den1+dout(ii)*din(ii)
   den2=den2+dout(ii)*hdelta(ii)
 end do

!DEBUG
!write(std_out,*)' hessupdt : den1,den2',den1,den2
!write(std_out,*)' din ',din
!write(std_out,*)' dout ',dout
!write(std_out,*)' hdelta ',hdelta
!ENDDEBUG

!Denominators are multiplicative
 den1=one/den1
 den3=one/den2

!Vectors which make a difference between the BROYDEN and
!the DAVIDON scheme.
 bfgs(:)=den1*din(:)-den3*hdelta(:)

!B.F.G.S. updating formula
 do ii=1,ndim
   do jj=1,ndim
     hessin(ii,jj)=hessin(ii,jj) +den1*din(ii)*din(jj) &
&     -den3*hdelta(ii)*hdelta(jj) +den2*bfgs(ii)*bfgs(jj)
   end do
 end do

end subroutine hessupdt
!!***

!----------------------------------------------------------------------

!!****f* m_bfgs/brdene
!! NAME
!! brdene
!!
!! FUNCTION
!! Update vin according to the Broyden formula, combined
!! with a line minimisation that take into account the total energies.
!! Also transfer vin to vin_prev, vout to vout_prev, and etotal to etotal_prev
!! Could see Numerical Recipes (Fortran), 1986, page 307,
!! as well as Schlegel, J. Comp. Chem. 3, 214 (1982) [[cite:Schlegel1982]].
!!
!! INPUTS
!!  etotal=new total energy (no meaning at output)
!!  hessin(ndim,ndim)=hessian matrix
!!  ndim=size of the hessian and vectors
!!  vout(ndim)=new output vector (no meaning at output)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  etotal_prev=previous total energy; contains input etotal at output
!!  vin(ndim)=new input vector; updated at output
!!  vin_prev(ndim)=previous input vector; contains input vin at output
!!  vout_prev(ndim)=previous output vector; contains input vout at output
!!
!! PARENTS
!!      m_pred_bfgs,m_pred_delocint
!!
!! CHILDREN
!!      findmin
!!
!! SOURCE

subroutine brdene(etotal,etotal_prev,hessin,ndim,vin,vin_prev,vout,vout_prev)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndim
 real(dp),intent(in) :: etotal
 real(dp),intent(inout) :: etotal_prev
!arrays+
 real(dp),intent(in) :: hessin(ndim,ndim),vout(ndim)
 real(dp),intent(inout) :: vin(ndim),vin_prev(ndim),vout_prev(ndim)

!Local variables-------------------------------
!scalars
 integer :: idim,brd_status
 real(dp) :: d2edv2_1,d2edv2_2,d2edv2_predict,dedv_1,dedv_2,dedv_min
 real(dp) :: dedv_predict,etotal_1,etotal_2,etotal_predict,lambda_1,lambda_2
 real(dp) :: lambda_predict
!arrays
 real(dp),allocatable :: dvin(:),vin_min(:),vout_min(:)

!***************************************************************************

 ABI_MALLOC(dvin,(ndim))
 ABI_MALLOC(vin_min,(ndim))
 ABI_MALLOC(vout_min,(ndim))

 lambda_1=1.0_dp       ; lambda_2=0.0_dp
 etotal_1=etotal      ; etotal_2=etotal_prev
 dvin(:)=vin(:)-vin_prev(:)
 dedv_1=dot_product(vout,dvin)
 dedv_2=dot_product(vout_prev,dvin)
 call findmin(dedv_1,dedv_2,dedv_predict,&
& d2edv2_1,d2edv2_2,d2edv2_predict,&
& etotal_1,etotal_2,etotal_predict,&
& lambda_1,lambda_2,lambda_predict,brd_status)

!DEBUG : comes back to usual BFGS !
!lambda_predict=1.0_dp
!dedv_predict=dedv_1
!ENDDEBUG

!Generates vin at the minimum, and an interpolated vout, modified
!to have the right value of dedv_predict, from findmin.
 vin_min(:)=vin_prev(:)+lambda_predict*dvin(:)
 vout_min(:)=vout_prev(:)+lambda_predict*(vout(:)-vout_prev(:))
 dedv_min=dedv_2+lambda_predict*(dedv_1-dedv_2)
!Modify vout_min in order to impose dedv_predict
 vout_min(:)=vout_min(:)+dvin(:)*(dedv_predict-dedv_min)/dot_product(dvin,dvin)

!Previous cartesian coordinates
 etotal_prev=etotal
 vin_prev(:)=vin(:)

!New atomic cartesian coordinates are obtained from vin, hessin and vout
 vin(:)=vin_min(:)
 do idim=1,ndim
   vin(:)=vin(:)-hessin(:,idim)*vout_min(idim)
 end do

!Previous atomic forces
 vout_prev(:)=vout(:)

 ABI_FREE(dvin)
 ABI_FREE(vin_min)
 ABI_FREE(vout_min)

end subroutine brdene
!!***

END MODULE m_bfgs
!!***
