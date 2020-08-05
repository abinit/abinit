!!****m* ABINIT/m_dynmat
!! NAME
!!  m_dynmat
!!
!! FUNCTION
!!  This module provides low-level tools to operate on the dynamical matrix
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2020 ABINIT group (XG, JCC, MJV, NH, RC, MVeithen, MM, MG, MT, DCA)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!!  Use more explicative names for the procedures!
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

module m_dynmat

 use defs_basis
 use m_abicore
 use m_errors
 use m_linalg_interfaces
 use m_xmpi

 use m_fstrings,        only : itoa, sjoin
 use m_numeric_tools,   only : wrap2_pmhalf, mkherm
 use m_symtk,           only : mati3inv, matr3inv, littlegroup_q
 use m_cgtools,         only : fxphas_seq
 use m_ewald,           only : ewald9
 use m_time,            only : timab

 implicit none

 private

 public :: asria_calc           ! Calculate the correction for the Acoustic sum rule on
                                !   the InterAtomic Forces or on the dynamical matrix directly
 public :: asria_corr           ! Imposition of the Acoustic sum rule on the InterAtomic Forces
                                !   or on the dynamical matrix directly from the previously calculated d2asr
 public :: asrprs               ! Imposition of the Acoustic sum rule on the InterAtomic Forces Plus Rotational Symmetry
 public :: cart29               ! Transform a second-derivative matrix from reduced  coordinates to cartesian coordinates, and also
                                !   1) add the ionic part of the effective charges,
                                !   2) normalize the electronic dielectric tensor, and add the vacuum polarisation
 public :: cart39               ! Transform a vector from reduced coordinates to cartesian coordinates,
                                !   taking into account the perturbation from which it was derived,
                                !   and also check the existence of the new values.
 public :: d2cart_to_red        ! Transform a second-derivative matrix
                                ! from cartesian to reduced coordinate.
 public :: chkph3               ! Check the completeness of the dynamical matrix
 public :: chneu9               ! Imposition of the Acoustic sum rule on the Effective charges
 public :: d2sym3               ! Build (nearly) all the other matrix elements that can be build using symmetries.
 public :: q0dy3_apply          ! Takes care of the inclusion of the ewald q=0 term in the dynamical matrix
 public :: q0dy3_calc           ! Calculate the q=0 correction term to the dynamical matrix
 ! TODO: 3 routines to symmetrize. Clarify different use cases
 public :: symdyma              ! Symmetrize the dynamical matrices
 public :: dfpt_sygra           ! Symmetrize derivatives of energy with respect to coordinates,
 public :: dfpt_sydy            ! Symmetrize dynamical matrix (eventually diagonal wrt to the atoms)
 public :: wings3               ! Suppress the wings of the cartesian 2DTE for which the diagonal element is not known
 public :: asrif9               ! Imposes the Acoustic Sum Rule to Interatomic Forces
 public :: get_bigbox_and_weights ! Compute
 public :: make_bigbox          ! Generates a Big Box of R points for the Fourier Transforms the dynamical matrix
 public :: bigbx9               ! Helper functions that faciliates the generation  of a Big Box containing
 public :: canat9               ! From reduced to canonical coordinates
 public :: canct9               ! Convert from canonical coordinates to cartesian coordinates
 public :: chkrp9               ! Check if the rprim used for the definition of the unit cell (in the
                                ! inputs) are consistent with the rprim used in the routine generating  the Big Box
 public :: dist9                ! Compute the distance between atoms in the big box
 public :: ftifc_q2r            ! Fourier transform of the dynamical matrices to obtain interatomic forces (real space).
 private :: ftifc_r2q            ! Fourier transform of the interatomic forces to obtain dynamical matrices (reciprocal space).
 public :: dynmat_dq            ! Compute the derivative D(q)/dq via Fourier transform of the interatomic forces
 public :: ifclo9               ! Convert from cartesian coordinates to local coordinates
 public :: wght9                ! Generates a weight to each R points of the Big Box and for each pair of atoms
 public :: d3sym                ! Given a set of calculated elements of the 3DTE matrix,
                                ! build (nearly) all the other matrix elements that can be build using symmetries.
 public :: sytens               ! Determines the set of irreductible elements of the non-linear optical susceptibility
                                ! and Raman tensors
 public :: axial9               ! Generates the local coordinates system from the  knowledge of the first vector (longitudinal) and
                                !   the ifc matrix in cartesian coordinates
 public :: dymfz9               ! Multiply the dynamical matrix by a phase shift to account for normalized canonical coordinates.
 public :: nanal9               ! Subtract/Add the non-analytical part from one dynamical matrix with number iqpt.
 public :: gtdyn9               ! Generates a dynamical matrix from interatomic force constants and
                                ! long-range electrostatic interactions.
 public :: dfpt_phfrq           ! Diagonalize IFC(q), return phonon frequencies and eigenvectors.
                                ! If q is Gamma, the non-analytical behaviour can be included.
 public :: pheigvec_normalize   ! Normalize input eigenvectors in cartesian coordinates.
 public :: phdispl_from_eigvec  ! Phonon displacements from eigenvectors
 public :: dfpt_prtph           ! Print phonon frequencies
 public :: massmult_and_breaksym  ! Multiply IFC(q) by atomic masses.

 ! TODO: Change name,
 public :: ftgam
 public :: ftgam_init


! *************************************************************************

contains
!!***

!!****f* m_dynmat/asria_calc
!! NAME
!! asria_calc
!!
!! FUNCTION
!! Calculate the correction for the Acoustic sum rule on the InterAtomic Forces
!! or on the dynamical matrix directly
!!
!! INPUTS
!! asr=(0 => no ASR, 1 or 2=> the diagonal element is modified to give the ASR,
!!      5 => impose hermitian solution using lapack call)
!! d2cart=matrix of second derivatives of total energy, in cartesian coordinates
!! mpert =maximum number of ipert
!! natom=number of atom
!!
!! OUTPUT
!! d2asr=matrix used to store the correction needed to fulfill
!! the acoustic sum rule.
!!
!! PARENTS
!!      dfpt_gatherdy,m_ddb,m_effective_potential_file
!!
!! CHILDREN
!!
!! SOURCE

subroutine asria_calc(asr,d2asr,d2cart,mpert,natom)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: asr,mpert,natom
!arrays
 real(dp),intent(in) :: d2cart(2,3,mpert,3,mpert)
 real(dp),intent(out) :: d2asr(2,3,natom,3,natom)

!Local variables-------------------------------
!scalars
 integer :: idir1,idir2,ii,ipert1,ipert2
 integer :: constrank, imatelem, iconst, nconst, nd2_packed, info
 !character(len=500) :: msg
!arrays
 integer, allocatable :: packingindex(:,:,:,:)
 real(dp), allocatable :: constraints(:,:,:)
 real(dp), allocatable :: d2cart_packed(:,:)
 real(dp), allocatable :: singvals(:)
 real(dp), allocatable :: constr_rhs(:,:)
 real(dp), allocatable :: work(:,:),rwork(:)

! *********************************************************************

 d2asr = zero

 if (asr==0) return

 !call wrtout(std_out,' asria_calc: calculation of the correction to the ASR for the interatomic forces.')
 do ipert1=1,natom
   do idir1=1,3
     do idir2=1,3

!      Compute d2asr
       do ipert2=1,natom
         d2asr(:,idir1,ipert1,idir2,ipert1)=&
&         d2asr(:,idir1,ipert1,idir2,ipert1)+&
&         d2cart(:,idir1,ipert1,idir2,ipert2)
       end do
     end do
   end do
 end do

!holistic method: overwrite d2asr with hermitian solution
 if (asr == 5) then
   nconst = 9*natom
   nd2_packed = 3*natom*(3*natom+1)/2
   ABI_MALLOC(constraints,(2,nconst, nd2_packed))
   ABI_MALLOC(d2cart_packed,(2,nd2_packed))
   ABI_MALLOC(constr_rhs,(2,nd2_packed))
   ABI_MALLOC(singvals,(nconst))
   ABI_MALLOC(work,(2,3*nd2_packed))
   ABI_MALLOC(rwork,(5*nd2_packed))
   ABI_MALLOC(packingindex,(3,natom,3,natom))
   ii=1
   packingindex=-1
   do ipert2=1,natom
     do idir2=1,3
       do ipert1=1,ipert2-1
         do idir1=1,3
           packingindex(idir1,ipert1,idir2,ipert2) = ii
           ii = ii+1
         end do
       end do
       do idir1=1,idir2
         packingindex(idir1,ipert2,idir2,ipert2) = ii
         ii = ii+1
       end do
     end do
   end do
!  setup constraint matrix
   constraints = zero
   do ipert1=1,natom
     do idir1=1,3
       do idir2=1,3
         iconst = idir2+3*(idir1-1 + 3*(ipert1-1))
!        set all atom forces, this component
         do ipert2=1,natom
           imatelem = packingindex(idir1,ipert1,idir2,ipert2)
           if (imatelem == -1) then
             imatelem = packingindex(idir2,ipert2,idir1,ipert1)
           end if
           constraints(1,iconst,imatelem) = one
         end do
       end do
     end do
   end do

   d2cart_packed = -999.0d0
   do ipert2=1,natom
     do idir2=1,3
       do ipert1=1,natom
         do idir1=1,3
           imatelem = packingindex(idir1,ipert1,idir2,ipert2)
           if (imatelem == -1) cycle
           d2cart_packed(:,imatelem) = d2cart(:,idir1,ipert1,idir2,ipert2)
         end do
       end do
     end do
   end do
   constr_rhs = zero
   constr_rhs(1,1:nconst) = matmul(constraints(1,:,:),d2cart_packed(1,:))
   constr_rhs(2,1:nconst) = matmul(constraints(1,:,:),d2cart_packed(2,:))

!  lwork = 3*nd2_packed
   call zgelss (nconst,nd2_packed,1,constraints,nconst,constr_rhs,nd2_packed,&
&   singvals,-one,constrank,work,3*nd2_packed,rwork,info)
   ABI_CHECK(info == 0, sjoin('zgelss returned:', itoa(info)))

!  unpack
   do ipert2=1,natom
     do idir2=1,3
       do ipert1=1,natom
         do idir1=1,3
           imatelem = packingindex(idir1,ipert1,idir2,ipert2)
           if (imatelem == -1) then
             imatelem = packingindex(idir2,ipert2,idir1,ipert1)
!            NOTE: should complex conjugate the correction below.
           end if
           d2asr(:,idir1,ipert1,idir2,ipert2) = constr_rhs(:,imatelem)
         end do
       end do
     end do
   end do

   ABI_FREE(constraints)
   ABI_FREE(d2cart_packed)
   ABI_FREE(singvals)
   ABI_FREE(constr_rhs)
   ABI_FREE(work)
   ABI_FREE(rwork)
   ABI_FREE(packingindex)
 end if

end subroutine asria_calc
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/asria_corr
!! NAME
!! asria_corr
!!
!! FUNCTION
!! Imposition of the Acoustic sum rule on the InterAtomic Forces
!! or on the dynamical matrix directly from the previously calculated d2asr
!!
!! INPUTS
!! asr=(0 => no ASR, 1 or 2=> the diagonal element is modified to give the ASR,
!!      5 => impose hermitian solution using lapack call)
!! d2asr=matrix used to store the correction needed to fulfill
!! the acoustic sum rule.
!! mpert =maximum number of ipert
!! natom=number of atom
!!
!! OUTPUT
!! Input/Output:
!! d2cart=matrix of second derivatives of total energy, in cartesian coordinates
!!
!! PARENTS
!!      ddb_elast,ddb_internalstr,dfpt_gatherdy,m_ddb,thmeig
!!
!! CHILDREN
!!
!! SOURCE

subroutine asria_corr(asr,d2asr,d2cart,mpert,natom)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: asr,mpert,natom
!arrays
 real(dp),intent(in) :: d2asr(2,3,natom,3,natom)
 real(dp),intent(inout) :: d2cart(2,3,mpert,3,mpert)

!Local variables-------------------------------
!scalars
 integer :: idir1,idir2,ipert1,ipert2

! *********************************************************************

 if (asr==0) return
 !call wrtout(std_out,' asria_corr: imposition of the ASR for the interatomic forces.')

 ! Remove d2asr
 do ipert2=1,natom
   do idir2=1,3
     do ipert1=1,natom
       do idir1=1,3
         d2cart(:,idir1,ipert1,idir2,ipert2)= d2cart(:,idir1,ipert1,idir2,ipert2) - d2asr(:,idir1,ipert1,idir2,ipert2)
       end do
     end do
   end do
 end do

end subroutine asria_corr
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/asrprs
!! NAME
!! asrprs
!!
!! FUNCTION
!! Imposition of the Acoustic sum rule on the InterAtomic Forces Plus Rotational Symmetry
!!
!! INPUTS
!!  asr=(3 => 1D systems, all elements are modified to give ASR and
!!            rotational symmetry)
!!      (4 => 0D systems, all elements are modified to give ASR and
!!            rotational symmetry)
!!  asrflg=(1 => the correction to enforce asr is computed from
!!           d2cart, but NOT applied;
!!          2 => one uses the previously determined correction)
!!  minvers=previously calculated inverted coefficient matrix
!!  mpert =maximum number of ipert
!!  natom=number of atom
!!  rotinv=(1,2,3 => for linear systems along x,y,z
!!          4 => non-linear molecule
!!  xcart=cartesian coordinates of the ions
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output:
!! d2cart=matrix of second derivatives of total energy, in cartesian coordinates
!! minvers=inverse of the supermatrix for future application of the corrections
!!
!! PARENTS
!!      m_ddb
!!
!! CHILDREN
!!
!! SOURCE

subroutine asrprs(asr,asrflag,rotinv,uinvers,vtinvers,singular,d2cart,mpert,natom,xcart)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: asr,asrflag,mpert,natom,rotinv
!arrays
 real(dp),intent(in) :: xcart(3,natom)
 real(dp),intent(inout) :: d2cart(2,3,mpert,3,mpert)
 real(dp),intent(inout) :: singular(1:3*natom*(3*natom-1)/2)
 real(dp),intent(inout) :: uinvers(1:3*natom*(3*natom-1)/2,1:3*natom*(3*natom-1)/2)
 real(dp),intent(inout) :: vtinvers(1:3*natom*(3*natom-1)/2,1:3*natom*(3*natom-1)/2)

!Local variables-------------------------------
!scalars
 integer :: column,idir1,idir2,ii,info,ipert1,ipert2,jj,n3,row,superdim
 real(dp) :: rcond,test
! real(dp) :: tau ! tau is present but commented out in this routine
 character(len=500) :: msg
!arrays
 integer :: check(3,natom,3)
 real(dp) :: tmp(natom,3,3),weightf(1:natom,1:natom)
 real(dp),allocatable :: d2cartold(:,:,:,:,:),d2vecc(:),d2veccnew(:),d2vecr(:)
 real(dp),allocatable :: d2vecrnew(:),superm(:,:),umatrix(:,:),vtmatrix(:)
 real(dp),allocatable :: work(:)

! *********************************************************************

 if(asr/=3 .and. asr/=4)then
   write(msg,'(3a,i0)')&
   'The argument asr should be 3 or 4,',ch10, 'however, asr = ',asr
   MSG_BUG(msg)
 end if

 if (asr==3.or.asr==4)then
   write(msg, '(a,a)' ) ch10, &
   'asrprs: imposition of the ASR for the interatomic forces and rotational invariance'
   call wrtout(std_out,msg)
 end if

 write(msg,'(a,i0)')' asrflag is ', asrflag
 call wrtout([std_out, ab_out], msg)

!variables for the dimensions of the matrices

!n1=3*natom*(3*natom-1)/2
!n2=9*natom
 n3=3*natom

 superdim=9*natom*(natom-1)/2+n3

 ABI_MALLOC(d2vecr,(1:superdim))
 ABI_MALLOC(d2vecc,(1:superdim))
 d2vecr=0d0
 d2vecc=0d0

!should be changed set to delta function for debugging
 weightf=1d0
!tau=1d-10
 do ii=1, natom
!  do jj=1, ii-1
!  weightf(ii,jj)= &
!  &     ((xcart(1,ii)-xcart(1,jj))**2+(xcart(2,ii)-xcart(2,jj))**2+(xcart(3,ii)-xcart(3,jj))**2)**tau
!  enddo
   weightf(ii,ii)=0d0
 end do

 ABI_MALLOC(d2cartold,(2,3,mpert,3,mpert))

 d2cartold=d2cart

!setup vector with uncorrected derivatives

 do ipert1=1, natom
   do ipert2=1, ipert1-1
     do idir1=1,3
       do idir2=1,3
         row=n3+9*(ipert1-1)*(ipert1-2)/2+9*(ipert2-1)+3*(idir1-1)+idir2
         if(abs(d2cart(1,idir1,ipert1,idir2,ipert2))<1d-6)then
           d2cart(1,idir1,ipert1,idir2,ipert2)=0d0
         else
           d2vecr(row)=4*weightf(ipert1,ipert2)*d2cart(1,idir1,ipert1,idir2,ipert2)
         end if
         if(abs(d2cart(2,idir1,ipert1,idir2,ipert2))<1d-6) then
           d2cart(2,idir1,ipert1,idir2,ipert2)=0d0
         else
           d2vecc(row)=4*weightf(ipert1,ipert2)*d2cart(2,idir1,ipert1,idir2,ipert2)
         end if
       end do
     end do
   end do
 end do

 if(asrflag==1) then !calculate the pseudo-inverse of the supermatrix
   ABI_MALLOC(superm,(1:superdim,1:superdim))

   superm=0d0

!  Setting up the supermatrix containing G, A, D

   do ipert1=1, natom
     do idir1=1, 3
!      Setting up G
       idir2=mod(idir1,3)+1
       row=3*(ipert1-1)+idir1
       do ipert2=1, ipert1-1
         column=9*(ipert1-1)*(ipert1-2)/2+9*(ipert2-1)+3*(rotinv-1)+idir1
         superm(column,row)=xcart(idir2,ipert2)-xcart(idir2,ipert1)
         column=9*(ipert1-1)*(ipert1-2)/2+9*(ipert2-1)+3*(rotinv-1)+idir2
         superm(column,row)=xcart(idir1,ipert1)-xcart(idir1,ipert2)
       end do
       do ipert2=ipert1+1, natom
         column=9*(ipert2-1)*(ipert2-2)/2+9*(ipert1-1)+3*(idir1-1)+rotinv
         superm(column,row)=xcart(idir2,ipert2)-xcart(idir2,ipert1)
         column=9*(ipert2-1)*(ipert2-2)/2+9*(ipert1-1)+3*(idir2-1)+rotinv
         superm(column,row)=xcart(idir1,ipert1)-xcart(idir1,ipert2)
       end do
     end do
     do idir1=1, 3
!      Setting up D
       idir2=mod(idir1,3)+1
       ii=mod(idir1+1,3)+1
       do ipert2=1, ipert1-1
         row=n3+9*(ipert1-1)*(ipert1-2)/2+9*(ipert2-1)+3*(rotinv-1)+idir1
         column=9*natom*(natom-1)/2+3*(ipert1-1)+idir1
         superm(column,row)=superm(column,row)+xcart(idir2,ipert2)-xcart(idir2,ipert1)
         column=9*natom*(natom-1)/2+3*(ipert1-1)+ii
         superm(column,row)=superm(column,row)+xcart(ii,ipert1)-xcart(ii,ipert2)
         row=n3+9*(ipert1-1)*(ipert1-2)/2+9*(ipert2-1)+3*(idir1-1)+rotinv
         column=9*natom*(natom-1)/2+3*(ipert2-1)+idir1
         superm(column,row)=superm(column,row)+xcart(idir2,ipert1)-xcart(idir2,ipert2)
         column=9*natom*(natom-1)/2+3*(ipert2-1)+ii
         superm(column,row)=superm(column,row)+xcart(ii,ipert2)-xcart(ii,ipert1)
       end do
!      Setting up A
       do idir2=1, 3
         do ipert2=1, ipert1-1
           column=9*(ipert1-1)*(ipert1-2)/2+9*(ipert2-1)+3*(idir1-1)+idir2
           row=n3+column
           superm(column,row)=4*weightf(ipert1,ipert2)
         end do
       end do
     end do
   end do

!  calculate the pseudo-inverse of the supermatrix

   ABI_MALLOC(work,(1:6*superdim))
   ABI_MALLOC(vtmatrix,(1:superdim))
   ABI_MALLOC(umatrix,(1:superdim,1:superdim))

!  singular value decomposition of superm

   call dgesvd('A','O',superdim,superdim,superm,superdim,singular,umatrix,superdim, &
               vtmatrix, 1, work,6*superdim,info)
   ABI_CHECK(info == 0, sjoin('dgesvd returned:', itoa(info)))

   ABI_FREE(vtmatrix)
   ABI_FREE(work)

   write(msg, '(a,es16.8,es16.8)' )' Largest and smallest values from svd', singular(1), singular(superdim)
   call wrtout([std_out, ab_out], msg)

!  Invert U and V**T, orthogonal matrices

   do ii=1, superdim
     do jj=1, superdim
       uinvers(ii,jj)=umatrix(jj,ii)
       vtinvers(ii,jj)=superm(jj,ii)
     end do
   end do

   ABI_FREE(umatrix)
   ABI_FREE(superm)

   write(msg,'(a,a)')' asrprs: done with asrflag 1', ch10
   call wrtout(std_out,msg)

 end if !asrflag=1

 if(asrflag==2) then

   ABI_MALLOC(d2vecrnew,(1:superdim))
   ABI_MALLOC(d2veccnew,(1:superdim))

!  Calculate V**T**-1 Sigma**-1 U**-1 *rhs

   do ii=1, superdim
     d2vecrnew(ii)=0d0
     d2veccnew(ii)=0d0
     do jj=1, superdim
       d2vecrnew(ii)=d2vecrnew(ii)+uinvers(ii,jj)*d2vecr(jj)
       d2veccnew(ii)=d2veccnew(ii)+uinvers(ii,jj)*d2vecc(jj)
     end do
   end do

   rcond=1d-10*singular(1)
   do ii=1, superdim
     if(singular(ii)>rcond) then
       d2vecrnew(ii)=d2vecrnew(ii)/singular(ii)
       d2veccnew(ii)=d2veccnew(ii)/singular(ii)
     else
       d2vecrnew(ii)=0d0
       d2veccnew(ii)=0d0
     end if
   end do

   do ii=1, superdim
     d2vecr(ii)=0d0
     d2vecc(ii)=0d0
     do jj=1, superdim
       d2vecr(ii)=d2vecr(ii)+vtinvers(ii,jj)*d2vecrnew(jj)
       d2vecc(ii)=d2vecc(ii)+vtinvers(ii,jj)*d2veccnew(jj)
     end do
   end do

!  Store vector back into the matrix of 2nd order derivates

   do ipert1=1, natom
     do ipert2=1, ipert1-1
       do idir1=1,3
         do idir2=1,3
           row=9*(ipert1-1)*(ipert1-2)/2+9*(ipert2-1)+3*(idir1-1)+idir2
           d2cart(1,idir1,ipert1,idir2,ipert2)=d2vecr(row)
           d2cart(2,idir1,ipert1,idir2,ipert2)=d2vecc(row)
           d2cart(1,idir2,ipert2,idir1,ipert1)=d2vecr(row)
           d2cart(2,idir2,ipert2,idir1,ipert1)=d2vecc(row)
         end do
       end do
     end do
   end do

   ABI_FREE(d2vecrnew)
   ABI_FREE(d2veccnew)

   check=0

   do ipert1=1, natom
     do idir1=1, 3
       do idir2=1, 3
         d2cart(1,idir1,ipert1,idir2,ipert1)=0d0
         d2cart(2,idir1,ipert1,idir2,ipert1)=0d0
         tmp(ipert1,idir1,idir2)=0d0
         do ipert2=1, natom
           if(ipert2/=ipert1) then
             tmp(ipert1,idir1,idir2)=tmp(ipert1,idir1,idir2) &
&             -d2cart(1,idir1,ipert1,idir2,ipert2) &
&             -d2cart(1,idir2,ipert2,idir1,ipert1)
           end if
         end do
       end do
     end do
   end do

   do ipert1=1, natom
     do idir1=1, 3
       do idir2=1, 3
         d2cart(1,idir1,ipert1,idir2,ipert1)=tmp(ipert1,idir1,idir2)/2
         d2cart(1,idir2,ipert1,idir1,ipert1)=d2cart(1,idir1,ipert1,idir2,ipert1)
       end do
     end do
   end do

   write(std_out,*) 'this should all be zero'

   do ipert1=1, natom
     do idir1=1, 3
       do idir2=1, 3
         test=0d0
         do ipert2=1, natom
           test=test+d2cart(1,idir1,ipert1,idir2,ipert2)+d2cart(1,idir2,ipert2,idir1,ipert1)
         end do
         write(std_out,'(i3,i3,i3,es11.3)') idir1,ipert1,idir2,test

         write(msg, '(i3,i3,i3,es11.3)' ) idir1,ipert1,idir2,test
         call wrtout(ab_out,msg)
       end do
     end do
   end do

   write(std_out,*) 'these as well'
   do ipert2=1, natom
     do idir1=1, 3
       do idir2=1, 3
         test=0d0
         do ipert1=1, natom
           test=test+d2cart(1,idir1,ipert1,idir2,ipert2)
         end do
         write(std_out,'(i3,i3,i3,i3,es11.3)') idir1,ipert1,idir2,ipert2,test
       end do
     end do
   end do

   write(msg,'(a,a)')' asrprs: done with asrflag 2', ch10
   call wrtout(std_out,msg)

 end if !ends asrflag=2

 ABI_FREE(d2vecr)
 ABI_FREE(d2vecc)
 ABI_FREE(d2cartold)

end subroutine asrprs
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/cart29
!! NAME
!! cart29
!!
!!
!! FUNCTION
!! Transform a second-derivative matrix from reduced
!! coordinates to cartesian coordinates, and also
!! 1) add the ionic part of the effective charges,
!! 2) normalize the electronic dielectric tensor, and
!!    add the vacuum polarisation
!!
!! INPUTS
!!  blkflg(3,mpert,3,mpert,nblok)=
!!   ( 1 if the element of the dynamical matrix has been calculated ;
!!     0 otherwise )
!!  blkval(2,3,mpert,3,mpert,nblok)=DDB values
!!  gprimd(3,3)=basis vector in the reciprocal space
!!  iblok=number of the blok that will be transformed
!!  mpert =maximum number of ipert
!!  natom=number of atom
!!  nblok=number of blocks (dimension of blkflg and blkval)
!!  ntypat=number of atom types
!!  rprimd(3,3)=basis vector in the real space
!!  typat(natom)=integer label of each type of atom (1,2,...)
!!  ucvol=unit cell volume
!!  zion(ntypat)=charge corresponding to the atom type
!!
!! OUTPUT
!!  carflg(3,mpert,3,mpert)= ( 1 if the element of the cartesian
!!  2DTE matrix has been calculated correctly ; 0 otherwise )
!!  d2cart(2,3,mpert,3,mpert)=
!!    dynamical matrix, effective charges, dielectric tensor,....
!!    all in cartesian coordinates
!!
!! PARENTS
!!      dfpt_gatherdy,m_ddb
!!
!! CHILDREN
!!
!! SOURCE

subroutine cart29(blkflg,blkval,carflg,d2cart,&
& gprimd,iblok,mpert,natom,nblok,ntypat,rprimd,typat,ucvol,zion)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,mpert,natom,nblok,ntypat
 real(dp),intent(in) :: ucvol
!arrays
 integer,intent(in) :: blkflg(3,mpert,3,mpert,nblok),typat(natom)
 integer,intent(out) :: carflg(3,mpert,3,mpert)
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok),gprimd(3,3),rprimd(3,3)
 real(dp),intent(in) :: zion(ntypat)
 real(dp),intent(out) :: d2cart(2,3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer :: idir1,idir2,ii,ipert1,ipert2
!arrays
 integer :: flg1(3),flg2(3)
 real(dp) :: vec1(3),vec2(3)

! *********************************************************************

!First, copy the data blok in place.
 d2cart(:,:,:,:,:)=blkval(:,:,:,:,:,iblok)

!Cartesian coordinates transformation (in two steps)
!First step
 do ipert1=1,mpert
   do ipert2=1,mpert
     do ii=1,2
       do idir1=1,3
         do idir2=1,3
           vec1(idir2)=d2cart(ii,idir1,ipert1,idir2,ipert2)
!          Note here blkflg
           flg1(idir2)=blkflg(idir1,ipert1,idir2,ipert2,iblok)
         end do
         call cart39(flg1,flg2,gprimd,ipert2,natom,rprimd,vec1,vec2)
         do idir2=1,3
           d2cart(ii,idir1,ipert1,idir2,ipert2)=vec2(idir2)
!          And here carflg
           carflg(idir1,ipert1,idir2,ipert2)=flg2(idir2)
         end do
       end do
     end do
   end do
 end do

!Second step
 do ipert1=1,mpert
   do ipert2=1,mpert
     do ii=1,2
       do idir2=1,3
         do idir1=1,3
           vec1(idir1)=d2cart(ii,idir1,ipert1,idir2,ipert2)
!          Note here carflg
           flg1(idir1)=carflg(idir1,ipert1,idir2,ipert2)
         end do
         call cart39(flg1,flg2,gprimd,ipert1,natom,rprimd,vec1,vec2)
         do idir1=1,3
           d2cart(ii,idir1,ipert1,idir2,ipert2)=vec2(idir1)
!          And here carflg again
           carflg(idir1,ipert1,idir2,ipert2)=flg2(idir1)
         end do
       end do
     end do
   end do
 end do

!For the dielectric tensor, takes into account the volume
!of the unit cell, and add the unit matrix (polarization of the vacuum)
 do idir1=1,3
   do idir2=1,3
     do ii=1,2
       d2cart(ii,idir1,natom+2,idir2,natom+2)=&
&       -four_pi/ucvol*d2cart(ii,idir1,natom+2,idir2,natom+2)
     end do
   end do
 end do

 do idir1=1,3
   d2cart(1,idir1,natom+2,idir1,natom+2)=&
&   1.0_dp+d2cart(1,idir1,natom+2,idir1,natom+2)
 end do

!Add the ionic charges to delta z to get the effective charges
 do ipert1=1,natom
   do idir1=1,3
     d2cart(1,idir1,ipert1,idir1,natom+2)=&
&     zion(typat(ipert1))+d2cart(1,idir1,ipert1,idir1,natom+2)
   end do
 end do
 do ipert2=1,natom
   do idir2=1,3
     d2cart(1,idir2,natom+2,idir2,ipert2)=&
&     zion(typat(ipert2))+d2cart(1,idir2,natom+2,idir2,ipert2)
   end do
 end do

!For the piezoelectric tensor, takes into account the volume of the unit cell
 do ipert2=natom+3,natom+4
   do idir1=1,3
     do idir2=1,3
       do ii=1,2
         d2cart(ii,idir1,natom+2,idir2,ipert2)=&
&         (1.0_dp/ucvol)*d2cart(ii,idir1,natom+2,idir2,ipert2)
         d2cart(ii,idir2,ipert2,idir1,natom+2)=&
&         (1.0_dp/ucvol)*d2cart(ii,idir2,ipert2,idir1,natom+2)
       end do
     end do
   end do
 end do

end subroutine cart29
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/cart39
!! NAME
!! cart39
!!
!!
!! FUNCTION
!! Transform a vector from reduced coordinates to cartesian coordinates,
!! taking into account the perturbation from which it was derived,
!! and also check the existence of the new values.
!!
!! INPUTS
!!  flg1(3)=tell if information of each component of vec1 is valid
!!  gprimd(3,3)=basis vector in the reciprocal space
!!  ipert=number of the perturbation
!!  natom=number of atom
!!  rprimd(3,3)=basis vector in the real space
!!  vec1(3)=input vector, in reduced coordinates
!!
!! OUTPUT
!!  flg2(3)=tell if information of each component of vec2 is valid
!!  vec2(3)=output vector, in cartesian coordinates
!!
!! PARENTS
!!      dfpt_gatherdy,m_ddb,m_dynmat
!!
!! CHILDREN
!!
!! SOURCE

subroutine cart39(flg1,flg2,gprimd,ipert,natom,rprimd,vec1,vec2)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: ipert,natom
!arrays
 integer,intent(in) :: flg1(3)
 integer,intent(out) :: flg2(3)
 real(dp),intent(in) :: gprimd(3,3),rprimd(3,3),vec1(3)
 real(dp),intent(out) :: vec2(3)

!Local variables -------------------------
!scalars
 integer :: idir,ii

! *********************************************************************

!Treat phonon-type perturbation
 if(ipert>=1.and.ipert<=natom)then

   do idir=1,3
     vec2(idir)=zero
     flg2(idir)=1
     do ii=1,3
       if(abs(gprimd(idir,ii))>1.0d-10)then
         if(flg1(ii)==1)then
           vec2(idir)=vec2(idir)+gprimd(idir,ii)*vec1(ii)
         else
           flg2(idir)=0
         end if
       end if
     end do
     if(flg2(idir)==0)vec2(idir)=zero
   end do

!  Treat electric field and qvec perturbations
 else if(ipert==natom+2.or.ipert==natom+8) then
!  OCL SCALAR
   do idir=1,3
     vec2(idir)=zero
     flg2(idir)=1
!    OCL SCALAR
     do ii=1,3
       if(abs(rprimd(idir,ii))>1.0d-10)then
         if(flg1(ii)==1)then
           vec2(idir)=vec2(idir)+rprimd(idir,ii)*vec1(ii)/two_pi
         else
           flg2(idir)=0
         end if
       end if
     end do
     if(flg2(idir)==0)vec2(idir)=zero
   end do

!  Treat other perturbations
 else
   do idir=1,3
     vec2(idir)=vec1(idir)
     flg2(idir)=flg1(idir)
   end do
 end if

end subroutine cart39
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/d2cart_to_red
!! NAME
!! d2cart_to_red
!!
!!
!! FUNCTION
!! Transform a second-derivative matrix from cartesian
!! coordinates to reduced coordinate. Also,
!! 1) remove the ionic part of the effective charges,
!! 2) remove the vacuum polarisation from the dielectric tensor
!!    and scale it with the unit cell volume
!! In short, does the inverse operation of cart29.
!!
!! INPUTS
!!  d2cart(2,3,mpert,3,mpert)=
!!    second-derivative matrix in cartesian coordinates
!!  gprimd(3,3)=basis vector in the reciprocal space
!!  rprimd(3,3)=basis vector in the real space
!!  mpert =maximum number of ipert
!!  natom=number of atom
!!
!! OUTPUT
!!  d2red(2,3,mpert,3,mpert)=
!!    second-derivative matrix in reduced coordinates
!!
!! PARENTS
!!      ddb_interpolate
!!
!! CHILDREN
!!
!! SOURCE

subroutine d2cart_to_red(d2cart, d2red, gprimd, rprimd, mpert, natom, &
&                        ntypat,typat,ucvol,zion)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,natom,ntypat
 real(dp),intent(in) :: ucvol
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: d2cart(2,3,mpert,3,mpert)
 real(dp),intent(out) :: d2red(2,3,mpert,3,mpert)
 real(dp),intent(in) :: gprimd(3,3),rprimd(3,3)
 real(dp),intent(in) :: zion(ntypat)

!Local variables -------------------------
!scalars
 integer :: idir1,idir2,ii,ipert1,ipert2
 real(dp) :: fac
!arrays
 integer :: flg1(3),flg2(3)
 real(dp) :: vec1(3),vec2(3)

! *********************************************************************

 flg1 = one
 flg2 = one

 d2red = d2cart

!Remove the ionic charges to z to get the change in effective charges
 do ipert1=1,natom
   do idir1=1,3
     d2red(1,idir1,ipert1,idir1,natom+2)=&
&     d2red(1,idir1,ipert1,idir1,natom+2) - zion(typat(ipert1))
   end do
 end do
 do ipert2=1,natom
   do idir2=1,3
     d2red(1,idir2,natom+2,idir2,ipert2)=&
&     d2red(1,idir2,natom+2,idir2,ipert2) - zion(typat(ipert2))
   end do
 end do

 ! Remove the vacuum polarizability from the dielectric tensor
 do idir1=1,3
   d2red(1,idir1,natom+2,idir1,natom+2)=&
&   d2red(1,idir1,natom+2,idir1,natom+2) - 1.0_dp
 end do

! Scale the dielectric tensor with the volue of the unit cell
 do idir1=1,3
   do idir2=1,3
     do ii=1,2
       d2red(ii,idir1,natom+2,idir2,natom+2)=&
&       - (ucvol / four_pi) * d2red(ii,idir1,natom+2,idir2,natom+2)
     end do
   end do
 end do

!For the piezoelectric tensor, takes into account the volume of the unit cell
 do ipert2=natom+3,natom+4
   do idir1=1,3
     do idir2=1,3
       do ii=1,2
         d2red(ii,idir1,natom+2,idir2,ipert2)=&
&         (ucvol)*d2red(ii,idir1,natom+2,idir2,ipert2)
         d2red(ii,idir2,ipert2,idir1,natom+2)=&
&         (ucvol)*d2red(ii,idir2,ipert2,idir1,natom+2)
       end do
     end do
   end do
 end do

! Reduced coordinates transformation (in two steps)
! Note that rprimd and gprimd are swapped, compared to what cart39 expects
! A factor of (2pi) ** 2 is added to transform the electric field perturbations

!First step
 do ipert1=1,mpert
   fac = one; if (ipert1==natom+2) fac = two_pi ** 2

   do ipert2=1,mpert
     do ii=1,2
       do idir1=1,3
         do idir2=1,3
           vec1(idir2)=d2red(ii,idir1,ipert1,idir2,ipert2)
         end do
         ! Transform vector from cartesian to reduced coordinates
         call cart39(flg1,flg2,transpose(rprimd),ipert1,natom,transpose(gprimd),vec1,vec2)
         do idir2=1,3
           d2red(ii,idir1,ipert1,idir2,ipert2)=vec2(idir2) * fac
         end do
       end do
     end do
   end do
 end do

!Second step
 do ipert1=1,mpert
   do ipert2=1,mpert
     fac = one; if (ipert2==natom+2) fac = two_pi ** 2

     do ii=1,2
       do idir2=1,3
         do idir1=1,3
           vec1(idir1)=d2red(ii,idir1,ipert1,idir2,ipert2)
         end do
         ! Transform vector from cartesian to reduced coordinates
         call cart39(flg1,flg2,transpose(rprimd),ipert2,natom,transpose(gprimd),vec1,vec2)
         do idir1=1,3
           d2red(ii,idir1,ipert1,idir2,ipert2)=vec2(idir1) * fac
         end do
       end do
     end do
   end do
 end do


end subroutine d2cart_to_red
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/chkph3
!! NAME
!! chkph3
!!
!! FUNCTION
!! Check the completeness of the dynamical matrix and eventually send a warning
!!
!! INPUTS
!!  carflg(3,mpert,3,mpert)= ( 1 if the element of the cartesian
!!  2DTE matrix has been calculated correctly ; 0 otherwise )
!!  idir = direction of the eventual electric field
!!  mpert =maximum number of ipert
!!  natom=number of atoms in unit cell
!!
!! OUTPUT
!!  eventually send a warning message
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!
!! SOURCE

subroutine chkph3(carflg,idir,mpert,natom)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: idir,mpert,natom
!arrays
 integer,intent(in) :: carflg(3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer :: idir1,idir2,ipert1,ipert2,send
 character(len=500) :: msg

! *********************************************************************

 send=0

!Check the elements of the analytical part of the dynamical matrix
 do ipert2=1,natom
   do idir2=1,3
     do ipert1=1,natom
       do idir1=1,3
         if(carflg(idir1,ipert1,idir2,ipert2)==0)then
           send=1
         end if
       end do
     end do
   end do
 end do

!If some electric field is present
 if(idir/=0)then

!  Check the dielectric constant
   if(carflg(idir,natom+2,idir,natom+2)==0)then
     send=1
   end if

!  Check the effective charges
   do ipert1=1,natom
     do idir1=1,3
       if(carflg(idir1,ipert1,idir,natom+2)==0)then
         send=1
       end if
     end do
   end do

 end if

 ! If needed, send the message
 if(send==1)then
   write(msg, '(a,a,a,a)' )&
   ' chkph3 : WARNING -',ch10,&
   '  The dynamical matrix was incomplete :',&
   ' phonon frequencies may be wrong ...'
   call wrtout([std_out, ab_out],msg)
 end if

end subroutine chkph3
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/chneu9
!! NAME
!! chneu9
!!
!! FUNCTION
!! Imposition of the Acoustic sum rule on the Effective charges
!! and suppress the imaginary part of the dynamical matrix
!!
!! INPUTS
!!  chneut=(0 => no ASR, 1 => equal repartition,2 => weighted repartition )
!!  mpert =maximum number of ipert
!!  natom=number of atom
!!  ntypat=number of types of atoms in unit cell
!!  selectz=selection of some parts of the effective charge tensor attached to one atom.
!!    (0=> no selection, 1=> trace only, 2=> symmetric part only)
!!  typat(natom)=type of the atom
!!  zion(ntypat)=atomic charge for every type of atom
!!
!! SIDE EFFECTS
!!  Input/Output
!!  d2cart=matrix of second derivatives of total energy, in cartesian
!!       coordinates
!!
!! PARENTS
!!      dfpt_gatherdy,m_ddb
!!
!! CHILDREN
!!
!! SOURCE

subroutine chneu9(chneut,d2cart,mpert,natom,ntypat,selectz,typat,zion)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: chneut,mpert,natom,ntypat,selectz
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: zion(ntypat)
 real(dp),intent(inout) :: d2cart(2,3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer :: idir1,idir2,ii,ipert1,ipert2
 character(len=500) :: msg
!arrays
 real(dp) :: sumwght(2)
 real(dp),allocatable :: wghtat(:)

! *********************************************************************

 ABI_MALLOC(wghtat,(natom))

!In case of acoustic sum rule imposition, compute the weights on each atom.
 if (chneut==1)then

!  The weight is the same for all atom
   do ipert1=1,natom
     wghtat(ipert1)=1./natom
   end do

 else if (chneut==2) then

!  The weight is proportional to the diagonal electronic screening charge of the atom
   sumwght(1)=zero
   do ipert1=1,natom
     wghtat(ipert1)=zero
     do idir1=1,3
       wghtat(ipert1)=wghtat(ipert1)+&
&       d2cart(1,idir1,ipert1,idir1,natom+2)+&
&       d2cart(1,idir1,natom+2,idir1,ipert1)-2*zion(typat(ipert1))
     end do
     sumwght(1)=sumwght(1)+wghtat(ipert1)
   end do

!  Normalize the weights to unity
   do ipert1=1,natom
     wghtat(ipert1)=wghtat(ipert1)/sumwght(1)
   end do
 end if

!Calculation of the violation of the charge neutrality
!and imposition of the charge neutrality condition
 if (chneut/=0)then
   write(msg, '(a,a,a,a,a,a,a)' )&
    ' The violation of the charge neutrality conditions',ch10,&
    ' by the effective charges is as follows :',ch10,&
    '    atom        electric field',ch10,&
    ' displacement     direction   '
   call wrtout(ab_out,msg)
   do idir1=1,3
     do idir2=1,3
       do ii=1,2
         sumwght(ii)=zero
         do ipert1=1,natom
           sumwght(ii)=sumwght(ii)+d2cart(ii,idir1,ipert1,idir2,natom+2)
         end do
         do ipert1=1,natom
           d2cart(ii,idir1,ipert1,idir2,natom+2)=&
           d2cart(ii,idir1,ipert1,idir2,natom+2)-sumwght(ii)*wghtat(ipert1)
         end do
       end do
       write(msg, '(i8,i16,2f16.6)' ) idir1,idir2,sumwght(1),sumwght(2)
       call wrtout(ab_out,msg)
     end do
   end do
   write(msg, '(a)' )' '
   call wrtout(ab_out,msg)

!  The same for the symmetrical part
   do idir1=1,3
     do idir2=1,3
       do ii=1,2
         sumwght(ii)=zero
         do ipert2=1,natom
           sumwght(ii)=sumwght(ii)+d2cart(ii,idir1,natom+2,idir2,ipert2)
         end do
         do ipert2=1,natom
           d2cart(ii,idir1,natom+2,idir2,ipert2)=&
           d2cart(ii,idir1,natom+2,idir2,ipert2)-sumwght(ii)*wghtat(ipert2)
         end do
       end do
     end do
   end do
 end if

!Selection of the trace of the effective charge tensor attached to each atom
 if(selectz==1)then
   do ipert1=1,natom
     do ii=1,2
       sumwght(ii)=zero
       do idir1=1,3
         sumwght(ii)=sumwght(ii)+d2cart(ii,idir1,ipert1,idir1,natom+2)
       end do
       do idir1=1,3
         do idir2=1,3
           d2cart(ii,idir1,ipert1,idir2,natom+2)=zero
         end do
       end do
       do idir1=1,3
         d2cart(ii,idir1,ipert1,idir1,natom+2)=sumwght(ii)/3.0_dp
       end do
     end do
   end do
!  Do the same for the symmetrical part of d2cart
   do ipert2=1,natom
     do ii=1,2
       sumwght(ii)=zero
       do idir1=1,3
         sumwght(ii)=sumwght(ii)+d2cart(ii,idir1,natom+2,idir1,ipert2)
       end do
       do idir1=1,3
         do idir2=1,3
           d2cart(ii,idir1,natom+2,idir2,ipert2)=zero
         end do
       end do
       do idir1=1,3
         d2cart(ii,idir1,natom+2,idir1,ipert2)=sumwght(ii)/3.0_dp
       end do
     end do
   end do
 end if

!Selection of the symmetric part of the effective charge tensor attached to each atom
 if(selectz==2)then
   do ipert1=1,natom
     do ii=1,2
       do idir1=1,3
         do idir2=1,3
           sumwght(ii)=(d2cart(ii,idir1,ipert1,idir2,natom+2)&
&           +d2cart(ii,idir2,ipert1,idir1,natom+2))/2.0_dp
           d2cart(ii,idir1,ipert1,idir2,natom+2)=sumwght(ii)
           d2cart(ii,idir2,ipert1,idir1,natom+2)=sumwght(ii)
         end do
       end do
     end do
   end do
!  Do the same for the symmetrical part of d2cart
   do ipert1=1,natom
     do ii=1,2
       do idir1=1,3
         do idir2=1,3
           sumwght(ii)=(d2cart(ii,idir1,ipert1,idir2,natom+2)&
&           +d2cart(ii,idir2,ipert1,idir1,natom+2))/2.0_dp
           d2cart(ii,idir1,ipert1,idir2,natom+2)=sumwght(ii)
           d2cart(ii,idir2,ipert1,idir1,natom+2)=sumwght(ii)
         end do
       end do
     end do
   end do
 end if

!Write the effective charge tensor
 write(msg, '(a,a,a,a,a,a,a)' )&
   ' Effective charge tensors after ',ch10,&
   ' imposition of the charge neutrality,',ch10,&
   ' and eventual restriction to some part :',ch10,&
  '   atom    displacement  '
 call wrtout(ab_out,msg)

 do ipert1=1,natom
   do idir1=1,3
     write(msg, '(2i10,3es16.6)' )ipert1,idir1,(d2cart(1,idir1,ipert1,idir2,natom+2),idir2=1,3)
     call wrtout(ab_out,msg)
   end do
 end do

!Zero the imaginary part of the dynamical matrix
 write(msg, '(a)' )' Now, the imaginary part of the dynamical matrix is zeroed '
 call wrtout(ab_out,msg)
 call wrtout(std_out,msg)

 do ipert1=1,natom
   do ipert2=1,natom
     do idir1=1,3
       do idir2=1,3
         d2cart(2,idir1,ipert1,idir2,ipert2)=zero
       end do
     end do
   end do
 end do

 ABI_FREE(wghtat)

end subroutine chneu9
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/d2sym3
!! NAME
!! d2sym3
!!
!! FUNCTION
!! Given a set of calculated elements of the 2DTE matrix d2,
!! build (nearly) all the other matrix elements that can be build using symmetries.
!!
!! 1. Perform first some completion by symmetrisation (exchange)
!!    over the two defining perturbations
!! 2. For each element, uses every symmetry, and build the element, in case
!!    EITHER all the needed elements are available,
!!    OR the only missing is itself
!!    OR the perturbation is the electric field, in a diamond
!!    symmetry (the last case was coded rather dirty)
!!
!! INPUTS
!!  indsym(4,nsym,natom)=indirect indexing array : for each
!!   isym,iatom, fourth element is label of atom into which iatom is sent by
!!   INVERSE of symmetry operation isym; first three elements are the primitive
!!   translations which must be subtracted after the transformation to get back
!!   to the original unit cell.
!!  mpert =maximum number of ipert
!!  natom= number of atoms
!!  nsym=number of space group symmetries
!!  qpt(3)=wavevector of the perturbation
!!  symq(4,2,nsym)= (integer) three first numbers define the G vector ;
!!   fourth number is zero if the q-vector is not preserved,
!!              is 1 otherwise
!!   second index is one without time-reversal symmetry,
!!                two with time-reversal symmetry
!!  symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!!  symrel(3,3,nsym)=3x3 matrices of the group symmetries (real space)
!!  timrev=1 if the time-reversal symmetry preserves the wavevector,
!!     modulo a reciprocal lattice vector, timrev=0 otherwise
!!  zero_by_symm= if 1, set blkflg to 1 for the elements that must be zero by symmetry, and zero them.
!!    This has the indirect effect of being able to resymmetrize the whole matrix, thus
!!    enforcing better the symmetry for the 2DTE.
!!
!! SIDE EFFECTS
!!  Input/Output
!!  d2(2,3,mpert,3,mpert)= matrix of the 2DTE
!!  blkflg(3,mpert,3,mpert)= ( 1 if the element of the dynamical
!!     matrix has been calculated ; 0 otherwise)
!!
!! NOTES
!!   The complete search would be to have the possibility
!!   of a set of missing elements. See notes of July 2, 1994,
!!   in the blue notebook 'computer codes'
!!   The partial solution adopted here takes into
!!   account some mirror symmetries
!!   as well as the tetrahedral symmetry of the diamond lattice
!!   On 010331, replaced the loops up to mpert by loops up to
!!   natom+2, because of a crash bug under Windows. However,
!!   the problem lies likely in the use of the indsym array.
!!
!! PARENTS
!!      completeperts,m_ddb,respfn
!!
!! CHILDREN
!!
!! SOURCE

subroutine d2sym3(blkflg,d2,indsym,mpert,natom,nsym,qpt,symq,symrec,symrel,timrev,zero_by_symm)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,natom,nsym,timrev,zero_by_symm
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symq(4,2,nsym)
 integer,intent(in),target :: symrec(3,3,nsym),symrel(3,3,nsym)
 integer,intent(inout) :: blkflg(3,mpert,3,mpert)
 real(dp),intent(in) :: qpt(3)
 real(dp),intent(inout) :: d2(2,3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 logical, parameter :: do_final_sym=.true.
 logical :: qzero
 integer :: exch12,found,idir1,idir2,idisy1,idisy2,ipert1,ipert2
 integer :: ipesy1,ipesy2,isgn,isym,ithree,itirev,nblkflg_is_one,noccur,nsym_used,quit,quit1
 real(dp) :: arg1,arg2,im,norm,re,sumi,sumr,xi,xr
!arrays
 integer,pointer :: sym1_(:,:,:),sym2_(:,:,:)
 real(dp),allocatable :: d2tmp1(:,:,:),d2tmp2(:,:,:),d2work(:,:,:,:,:)

! *********************************************************************

 qzero=(qpt(1)**2+qpt(2)**2+qpt(3)**2<tol16)

!Here look after exchange of 1 and 2 axis,
!for electric field in diamond symmetry
 exch12=0
 if (qzero) then
   do isym=1,nsym
     exch12=1
     if(symrel(1,1,isym)/=0)exch12=0
     if(symrel(1,2,isym)/=1)exch12=0
     if(symrel(1,3,isym)/=0)exch12=0
     if(symrel(2,1,isym)/=1)exch12=0
     if(symrel(2,2,isym)/=0)exch12=0
     if(symrel(2,3,isym)/=0)exch12=0
     if(symrel(3,1,isym)/=0)exch12=0
     if(symrel(3,2,isym)/=0)exch12=0
     if(symrel(3,3,isym)/=1)exch12=0
!    if(exch12==1) write(std_out,*)' d2sym3 : found exchange 1 2 =',isym
     if(exch12==1)exit
   end do
 end if

!Exchange of perturbations

!Consider two cases : either time-reversal symmetry
!conserves the wavevector, or not
 if(timrev==0)then

!  do ipert1=1,mpert  See notes
   do ipert1=1,min(natom+2,mpert)
     do idir1=1,3

!      Since the matrix is hermitian, the diagonal elements are real
       d2(2,idir1,ipert1,idir1,ipert1)=zero

!      do ipert2=1,mpert See notes
       do ipert2=1,min(natom+2,mpert)
         do idir2=1,3

!          If an element exists
           if(blkflg(idir1,ipert1,idir2,ipert2)==1)then

!            Either complete the symmetric missing element
             if(blkflg(idir2,ipert2,idir1,ipert1)==0)then

               d2(1,idir2,ipert2,idir1,ipert1)= d2(1,idir1,ipert1,idir2,ipert2)
               d2(2,idir2,ipert2,idir1,ipert1)=-d2(2,idir1,ipert1,idir2,ipert2)

               blkflg(idir2,ipert2,idir1,ipert1)=1

!              Or symmetrize (the matrix is hermitian) in case both exists
!              (Note : this opportunity has been disabled for more
!              obvious search for bugs in the code )
!              else
!              sumr=d2(1,idir2,ipert2,idir1,ipert1)+d2(1,idir1,ipert1,idir2,ipert2)
!              sumi=d2(1,idir2,ipert2,idir1,ipert1)-d2(1,idir1,ipert1,idir2,ipert2)
!              d2(1,idir2,ipert2,idir1,ipert1)=half*sumr
!              d2(1,idir1,ipert1,idir2,ipert2)=half*sumr
!              d2(2,idir2,ipert2,idir1,ipert1)=half*sumi
!              d2(2,idir1,ipert1,idir2,ipert2)=-half*sumi
             end if
           end if

         end do
       end do

     end do
   end do

!  Here, case with time-reversal symmetry
 else

!  do ipert1=1,mpert See notes
   do ipert1=1,min(natom+2,mpert)
     do idir1=1,3
!      do ipert2=1,mpert See notes
       do ipert2=1,min(natom+2,mpert)
         do idir2=1,3
           d2(2,idir1,ipert1,idir2,ipert2)=zero

!          If an element exists
           if(blkflg(idir1,ipert1,idir2,ipert2)==1)then

!            Either complete the symmetric missing element
             if(blkflg(idir2,ipert2,idir1,ipert1)==0)then

               d2(1,idir2,ipert2,idir1,ipert1)=d2(1,idir1,ipert1,idir2,ipert2)
               blkflg(idir2,ipert2,idir1,ipert1)=1

!              Or symmetrize (the matrix is hermitian) in case both exists
!              (Note : this opportunity has been disabled for more
!              obvious search for bugs in the code )
!              else
!              sumr=d2(1,idir2,ipert2,idir1,ipert1)+d2(1,idir1,ipert1,idir2,ipert2)
!              d2(1,idir2,ipert2,idir1,ipert1)=half*sumr
!              d2(1,idir1,ipert1,idir2,ipert2)=half*sumr
             end if

           end if
         end do
       end do
     end do
   end do
 end if

!Big Big Loop : symmetrize three times, because
!of some cases in which one element is not yet available
!at the first pass, and even at the second one !
 do ithree=1,3

!  Big loop on all elements
!  do ipert1=1,mpert See notes
   do ipert1=1,min(natom+2,mpert)

!    Select the symmetries according to pertubation 1
     if (ipert1<=natom)then
       sym1_ => symrec
     else
       sym1_ => symrel
     end if

     do idir1=1,3
!      do ipert2=1,mpert See notes
       do ipert2=1,min(natom+2,mpert)

    !    Select the symmetries according to pertubation 2
         if (ipert2<=natom)then
           sym2_ => symrec
         else
           sym2_ => symrel
         end if

         do idir2=1,3

!          Will get element (idir1,ipert1,idir2,ipert2)
!          so this element should not yet be present ...
           if(blkflg(idir1,ipert1,idir2,ipert2)/=1)then

             d2(1,idir1,ipert1,idir2,ipert2)=zero
             d2(2,idir1,ipert1,idir2,ipert2)=zero

!            Loop on all symmetries, including time-reversal
             quit1=0
             do isym=1,nsym
               do itirev=1,2
                 isgn=3-2*itirev

                 if(symq(4,itirev,isym)/=0)then
                   found=1

!                  Here select the symmetric of ipert1
                   if(ipert1<=natom)then
                     ipesy1=indsym(4,isym,ipert1)
                   else if(ipert1==(natom+2).and.qzero)then
                     ipesy1=ipert1
                   else
                     found=0
                   end if

!                  Here select the symmetric of ipert2
                   if(ipert2<=natom)then
                     ipesy2=indsym(4,isym,ipert2)
                   else if(ipert2==(natom+2).and.qzero)then
                     ipesy2=ipert2
                   else
                     found=0
                   end if

!                  Now that a symmetric perturbation has been obtained,
!                  including the expression of the symmetry matrix, see
!                  if the symmetric values are available
                   if( found==1 ) then

                     sumr=zero
                     sumi=zero
                     noccur=0
                     nblkflg_is_one=0
                     quit=0
                     do idisy1=1,3
                       do idisy2=1,3
                         if(sym1_(idir1,idisy1,isym)/=0 .and. sym2_(idir2,idisy2,isym)/=0 )then
                           if(blkflg(idisy1,ipesy1,idisy2,ipesy2)==1)then
                             nblkflg_is_one=nblkflg_is_one+1
                             sumr=sumr+sym1_(idir1,idisy1,isym)*sym2_(idir2,idisy2,isym)*&
&                             d2(1,idisy1,ipesy1,idisy2,ipesy2)
                             sumi=sumi+sym1_(idir1,idisy1,isym)*sym2_(idir2,idisy2,isym)*&
&                             d2(2,idisy1,ipesy1,idisy2,ipesy2)

!                            Here, in case the symmetric of the element
!                            is the element, or the symmetric with
!                            respect to permutation of perturbations
!                            (some more conditions on the time-reversal
!                            symmetry must be fulfilled although)
                           else if(  idisy1==idir1 .and. ipesy1==ipert1&
&                             .and. idisy2==idir2 .and. ipesy2==ipert2&
&                             .and.(isgn==1 .or. timrev==1 &
&                             .or. (idir1==idir2 .and. ipert1==ipert2)))&
&                             then
                             noccur=noccur+sym1_(idir1,idisy1,isym)*sym2_(idir2,idisy2,isym)
                           else if(  idisy1==idir2 .and. ipesy1==ipert2&
&                             .and. idisy2==idir1 .and. ipesy2==ipert1&
&                             .and.(isgn==-1 .or. timrev==1&
&                             .or. (idir1==idir2 .and. ipert1==ipert2)))&
&                             then
                             noccur=noccur+sym1_(idir1,idisy1,isym)*sym2_(idir2,idisy2,isym)

!                            Here, electric field case
                           else if( exch12==1 .and. &
&                             ipert1==natom+2 .and. ipert2==natom+2&
&                             .and.(( idisy1+idir1 ==3 &
&                             .and. idisy2==3 .and. idir2==3)&
&                             .or. ( idisy1+idir2 ==3&
&                             .and. idisy2==3 .and. idir1==3)&
&                             .or. ( idisy2+idir2 ==3&
&                             .and. idisy1==3 .and. idir1==3)&
&                             .or. ( idisy2+idir1 ==3&
&                             .and. idisy1==3 .and. idir2==3)))&
&                             then
                             noccur=noccur+sym1_(idir1,idisy1,isym)*sym2_(idir2,idisy2,isym)

                           else
!                            Not found
                             found=0
                             quit=1
                             exit
                           end if

                         end if
                       end do
                       if(quit==1)exit
                     end do
                   end if

!                  In case zero_by_symm==0, the computed materix element must be associated to at least one really computed matrix element
                   if(zero_by_symm==0 .and. nblkflg_is_one==0)then
                     found=0
                   endif

!                  Now, if still found and associated to at least one really computed matrix element, put the correct value into array d2
                   if(found==1)then

!                    In case of phonons, need to take into account the
!                    time-reversal symmetry, and the shift back to the unit cell
!
!                    XG990712 : I am not sure this must be kept for electric field ...
!                    1) Consider time-reversal symmetry
                     sumi=isgn*sumi

                     if(ipert1<=natom .and. ipert2<=natom)then
!                      2) Shift the atoms back to the unit cell.
                       arg1=two_pi*( qpt(1)*indsym(1,isym,ipert1)&
&                       +qpt(2)*indsym(2,isym,ipert1)&
&                       +qpt(3)*indsym(3,isym,ipert1) )
                       arg2=two_pi*( qpt(1)*indsym(1,isym,ipert2)&
&                       +qpt(2)*indsym(2,isym,ipert2)&
&                       +qpt(3)*indsym(3,isym,ipert2) )
                       re=cos(arg1)*cos(arg2)+sin(arg1)*sin(arg2)
!                      XG010117 Must use isgn
                       im=isgn*(cos(arg2)*sin(arg1)-cos(arg1)*sin(arg2))
                     else
                       re=one
                       im=zero
                     end if

!                    Final check, could still fail if the
!                    element was its own symmetric
                     if (abs(1.0_dp-re*noccur)< 1.0d-6.and.abs(im*noccur) <1.0d-6) then
                       found=0
                     end if

                   end if

                   if(found==1)then

                     if(noccur==0)then
                       d2(1,idir1,ipert1,idir2,ipert2)=re*sumr-im*sumi
                       d2(2,idir1,ipert1,idir2,ipert2)=re*sumi+im*sumr
                     else
!                      See page July 2, 1994 in computer codes notebook
                       xr=re*sumr-im*sumi
                       xi=re*sumi+im*sumr
                       norm=one+noccur**2-two*re*noccur
                       xr=xr/norm
                       xi=xi/norm
                       d2(1,idir1,ipert1,idir2,ipert2)=&
&                       (one-re*noccur)*xr-im*noccur*xi
                       d2(2,idir1,ipert1,idir2,ipert2)=&
&                       (one-re*noccur)*xi+im*noccur*xr
                     end if

!                    The element has been constructed !
                     blkflg(idir1,ipert1,idir2,ipert2)=1

                     quit1=1
                     exit ! Exit loop on symmetry operations
                   end if

!                  End loop on all symmetries + time-reversal
                 end if
               end do
               if(quit1==1)exit
             end do

           end if
         end do ! End big loop on all elements
       end do
     end do
   end do

 end do !  End Big Big Loop

!MT oct. 20, 2014:
!Once the matrix has been built, it does not necessarily fulfill the correct symmetries.
!It has just been filled up from rows or columns that only fulfill symmetries preserving
!one particular perturbation.
!An additional symmetrization might solve this (do not consider TR-symmetry)
 if (do_final_sym) then
   ABI_MALLOC(d2tmp1,(2,3,3))
   ABI_MALLOC(d2tmp2,(2,3,3))
   ABI_MALLOC(d2work,(2,3,mpert,3,mpert))
   d2work(:,:,:,:,:)=d2(:,:,:,:,:)
   do ipert1=1,min(natom+2,mpert)
     if ((ipert1==natom+1.or.ipert1==natom+10.or.ipert1==natom+11).or.(ipert1==natom+2.and.(.not.qzero))) cycle
     if (ipert1<=natom)then
       sym1_ => symrec
     else
       sym1_ => symrel
     end if
     do ipert2=1,min(natom+2,mpert)
!      if (any(blkflg(:,ipert1,:,ipert2)==0)) cycle
       if ((ipert2==natom+1.or.ipert2==natom+10.or.ipert2==natom+11).or.(ipert2==natom+2.and.(.not.qzero))) cycle
       if (ipert2<=natom)then
         sym2_ => symrec
       else
         sym2_ => symrel
       end if
       nsym_used=0
       d2tmp2(:,:,:)=zero
       do isym=1,nsym
         if (symq(4,1,isym)==1) then
           ipesy1=ipert1;if (ipert1<=natom) ipesy1=indsym(4,isym,ipert1)
           ipesy2=ipert2;if (ipert2<=natom) ipesy2=indsym(4,isym,ipert2)
!          The condition on next line is too severe, since some elements of sym1_ or sym2_ might be zero,
!          which means not all blkflg(:,ipesy1,:,ipesy2) would need to be 1 to symmetrize the matrix.
!          However, coding something more refined is really more difficult.
!          This condition then has the side effect that more symmetries can be applied when zero_by_symm==1,
!          since blkflg can be set to 1 when the symmetries guarantee the matrix element to be zero.
           if (all(blkflg(:,ipesy1,:,ipesy2)==1)) then
             nsym_used=nsym_used+1
             re=one;im=zero
             if (ipert1<=natom.and.ipert2<=natom.and.(.not.qzero)) then
               arg1=two_pi*(qpt(1)*(indsym(1,isym,ipert1)-indsym(1,isym,ipert2)) &
&                          +qpt(2)*(indsym(2,isym,ipert1)-indsym(2,isym,ipert2)) &
&                          +qpt(3)*(indsym(3,isym,ipert1)-indsym(3,isym,ipert2)))
               re=cos(arg1);im=sin(arg1)
             end if
             d2tmp1(:,:,:)=zero
             do idir2=1,3 !kappa
               do idir1=1,3 !mu
                 do idisy1=1,3 !nu
                   d2tmp1(:,idir1,idir2)=d2tmp1(:,idir1,idir2) &
&                     +sym1_(idir1,idisy1,isym)*d2(:,idisy1,ipesy1,idir2,ipesy2)
                 end do
               end do
             end do
             do idir2=1,3 !mu
               do idir1=1,3 !kappa
                 do idisy2=1,3 !nu
                   d2tmp2(1,idir1,idir2)=d2tmp2(1,idir1,idir2) &
&                  +sym2_(idir2,idisy2,isym)*(d2tmp1(1,idir1,idisy2)*re-d2tmp1(2,idir1,idisy2)*im)
                   d2tmp2(2,idir1,idir2)=d2tmp2(2,idir1,idir2) &
&                  +sym2_(idir2,idisy2,isym)*(d2tmp1(1,idir1,idisy2)*im+d2tmp1(2,idir1,idisy2)*re)
                 end do
               end do
             end do
           end if
         end if
       end do ! isym
       if (nsym_used>0) d2work(:,1:3,ipert1,1:3,ipert2)=d2tmp2(:,1:3,1:3)/dble(nsym_used)
     end do !ipert2
   end do !ipert1
   if (mpert>=natom)   d2(:,1:3,1:natom,1:3,1:natom)=d2work(:,1:3,1:natom,1:3,1:natom)
   if (mpert>=natom+2) then
     d2(:,1:3,natom+2,1:3,1:natom)=d2work(:,1:3,natom+2,1:3,1:natom)
     d2(:,1:3,1:natom,1:3,natom+2)=d2work(:,1:3,1:natom,1:3,natom+2)
     d2(:,1:3,natom+2,1:3,natom+2)=d2work(:,1:3,natom+2,1:3,natom+2)
   end if
   ABI_FREE(d2tmp1)
   ABI_FREE(d2tmp2)
   ABI_FREE(d2work)
 end if

end subroutine d2sym3
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/q0dy3_apply
!! NAME
!! q0dy3_apply
!!
!! FUNCTION
!! Takes care of the inclusion of the ewald q=0 term in the dynamical
!! matrix - corrects the dyew matrix provided as input
!!
!! INPUTS
!!  dyewq0(3,3,natom) = part needed to correct
!!    the dynamical matrix for atom self-interaction.
!!  natom= number of atom in the unit cell
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  dyew(2,3,natom,3,natom)= dynamical matrix corrected on output
!!
!! NOTES
!! Should be used just after each call to dfpt_ewald, for both
!! q==0 and the real wavelength.
!!
!! The q0dy3_apply should be used in conjunction with the subroutine dfpt_ewald (or ewald9):
!! First, the call of dfpt_ewald with q==0 should be done,
!!   then the call to q0dy3_calc will produce
!!   the dyewq0 matrix from the (q=0) dyew matrix
!! Second, the call of dfpt_ewald with the real q (either =0 or diff 0)
!!   should be done, then the call to q0dy3_apply
!!   will produce the correct dynamical matrix dyew starting from
!!   the previously calculated dyewq0 and the bare(non-corrected)
!!   dyew matrix
!!
!! PARENTS
!!      m_dynmat,m_ifc,respfn
!!
!! CHILDREN
!!
!! SOURCE

subroutine q0dy3_apply(natom,dyewq0,dyew)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom
!arrays
 real(dp),intent(in) :: dyewq0(3,3,natom)
 real(dp),intent(inout) :: dyew(2,3,natom,3,natom)

!Local variables -------------------------
!scalars
 integer :: ia,mu,nu

! *********************************************************************

 do mu=1,3
   do nu=1,3
     do ia=1,natom
       dyew(1,mu,ia,nu,ia)=dyew(1,mu,ia,nu,ia)-dyewq0(mu,nu,ia)
     end do
   end do
 end do

end subroutine q0dy3_apply
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/q0dy3_calc
!! NAME
!! q0dy3_calc
!!
!! FUNCTION
!! Calculate the q=0 correction term to the dynamical matrix
!!
!! INPUTS
!!  dyew(2,3,natom,3,natom)= dynamical matrix
!!    input, non-corrected, for q=0 if option=1 or 2
!!  natom= number of atom in the unit cell
!!  option= either 1 or 2:
!!     1: use dyew to calculate dyewq0 symmetrical form
!!     2: use dyew to calculate dyewq0 symmetrical form
!!
!! OUTPUT
!!  dyewq0(3,3,natom) = part needed to correct
!!    the dynamical matrix for atom self-interaction.
!!
!! NOTES
!! Should be used just after each call to dfpt_ewald, for both
!! q==0 and the real wavelength.
!!
!! If option=1 or 2, q0dy3_calc uses an Ewald dynamical matrix at q=0,
!! called dyew, to produce a contracted form called dyewq0 :
!! either:
!!   in an unsymmetrical form (if option=1), or
!!   in a symmetrical form (if option=2).
!!
!! The q0dy3_calc should be used in conjunction with the subroutine dfpt_ewald (or ewald9).
!! First, the call of dfpt_ewald with q==0 should be done ,
!!   then the call to q0dy3_calc will produce
!!   the dyewq0 matrix from the (q=0) dyew matrix
!! Second, the call of dfpt_ewald with the real q (either =0 or diff 0)
!!   should be done, then the call to q0dy3_apply
!!   will produce the correct dynamical matrix dyew starting from
!!   the previously calculated dyewq0 and the bare(non-corrected)
!!   dyew matrix
!!
!! PARENTS
!!      ddb_hybrid,m_ifc,respfn
!!
!! CHILDREN
!!
!! SOURCE

subroutine q0dy3_calc(natom,dyewq0,dyew,option)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom,option
!arrays
 real(dp),intent(in) :: dyew(2,3,natom,3,natom)
 real(dp),intent(out) :: dyewq0(3,3,natom)

!Local variables -------------------------
!scalars
 integer :: ia,ib,mu,nu
 character(len=500) :: msg

! *********************************************************************

 if(option==1.or.option==2)then
   do mu=1,3
     do nu=1,3
       do ia=1,natom
         dyewq0(mu,nu,ia)=zero
         do ib=1,natom
           dyewq0(mu,nu,ia)=dyewq0(mu,nu,ia)+dyew(1,mu,ia,nu,ib)
         end do
       end do
     end do
   end do
 else
   write (msg, '(3a)')&
&   'option should be 1 or 2.',ch10,&
&   'action: correct calling routine'
   MSG_BUG(msg)
 end if

 if(option==2)then
   do ia=1,natom
     do mu=1,3
       do nu=mu,3
         dyewq0(mu,nu,ia)=(dyewq0(mu,nu,ia)+dyewq0(nu,mu,ia))/2
         dyewq0(nu,mu,ia)=dyewq0(mu,nu,ia)
       end do
     end do
   end do
 end if

end subroutine q0dy3_calc
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/symdyma
!!
!! NAME
!! symdyma
!!
!! FUNCTION
!! Symmetrize the dynamical matrices
!!
!! INPUTS
!! indsym(4,nsym*natom)=indirect indexing array : for each
!!   isym,iatom, fourth element is label of atom into which iatom is sent by
!!   INVERSE of symmetry operation isym; first three elements are the primitive
!!   translations which must be subtracted after the transformation to get back
!!   to the original unit cell.
!! natom=number of atoms in unit cell
!! nsym=number of space group symmetries
!! qptn(3)=normalized phonon wavevector
!! rprimd(3,3)=dimensional primitive translations (bohr)
!! symrel(3,3,nsym)=3x3 matrices of the group symmetries (real space)
!!
!! SIDE EFFECTS
!! Input/Output
!! dmati(2*3*natom*3*natom)=dynamical matrices in cartesian coordinates relative to the q
!!  points of the B.Z. sampling
!!
!! NOTES
!! the procedure of the symmetrization of the dynamical matrix follows the
!! equations in: Hendrikse et al., Computer Phys. Comm. 86, 297 (1995) [[cite:Hendrikse1995]]
!!
!! TODO
!! A full description of the equations should be included
!!
!! PARENTS
!!      m_dynmat,m_phgamma,relaxpol
!!
!! CHILDREN
!!
!! SOURCE

subroutine symdyma(dmati,indsym,natom,nsym,qptn,rprimd,symrel,symafm)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom,nsym
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symrel(3,3,nsym)
 integer,intent(in) :: symafm(nsym)
 real(dp),intent(in) :: qptn(3),rprimd(3,3)
 real(dp),intent(inout) :: dmati(2*3*natom*3*natom)

!Local variables -------------------------
!scalars
 integer :: i1,i2,iat,idir,ii,index,isgn,isym,itirev,jat,jdir,jj,kk,ll
 integer :: niat,njat,timrev
 real(dp) :: arg1,arg2,dmint,im,re,sumi,sumr
!arrays
 integer :: indij(natom,natom),symq(4,2,nsym),symrec(3,3,nsym)
 real(dp) :: TqR(3,3),TqS_(3,3),dynmat(2,3,natom,3,natom)
 real(dp) :: dynmatint(2*nsym,2,3,natom,3,natom),gprimd(3,3)
 real(dp) :: symcart(3,3,nsym)

! *********************************************************************

 ! 0) initializations
 call matr3inv(rprimd,gprimd)
 do isym=1,nsym
   call mati3inv(symrel(:,:,isym),symrec(:,:,isym))
 end do

 TqR=zero
 TqS_=zero
 dynmat=zero

!Note: dynmat is used as work space here
 i1=0
 do iat=1,natom
   do idir=1,3
     i1=i1+1
     i2=0
     do jat=1,natom
       do jdir=1,3
         i2=i2+1
         index=i1+3*natom*(i2-1)
         dynmat(1,idir,iat,jdir,jat)=dmati(2*index-1)
         dynmat(2,idir,iat,jdir,jat)=dmati(2*index  )
       end do
     end do
   end do
 end do

!Transform symrel to cartesian coordinates (RC coding)
!do isym=1,nsym
!symcart(:,:,isym)=matmul(rprimd,matmul(dble(symrel(:,:,isym)),gprimd))
!end do

!Coding from symdm9
 do isym=1,nsym
   do jj=1,3
     symcart(:,jj,isym)=zero
     do kk=1,3
       do ll=1,3
         symcart(:,jj,isym)=symcart(:,jj,isym)+rprimd(:,kk)*gprimd(jj,ll)*symrel(kk,ll,isym)
       end do
     end do
   end do
 end do

 ! Get the symq of the CURRENT Q POINT
 ! mjv: set prtvol=0 for production runs.
 call littlegroup_q(nsym,qptn,symq,symrec,symafm,timrev,prtvol=0)

 indij(:,:)=0
 dynmatint=zero

 do isym=1,nsym  ! loop over all the symmetries
   ! write(std_out,*) 'current symmetry',isym
   do itirev=1,2  ! loop over the time-reversal symmetry
     isgn=3-2*itirev
     ! write(std_out,*) 'timereversal',isgn

     if (symq(4,itirev,isym)==1) then ! isym belongs to the wave vector point group
       ! write(std_out,*) 'isym belongs to the wave vector point group'
       do iat=1,natom
         do jat=1,natom
           niat=indsym(4,isym,iat)  ! niat={R|t}iat
           njat=indsym(4,isym,jat)  ! njat={R|t}jat
           indij(niat,njat)=indij(niat,njat)+1
           ! write(std_out,'(a,5i5)') 'current status:',iat,jat,niat,njat,indij(niat,njat)
           ! phase calculation, arg1 and arg2 because of two-atom derivative
           arg1=two_pi*( qptn(1)*indsym(1,isym,iat)+&
             qptn(2)*indsym(2,isym,iat)+&
             qptn(3)*indsym(3,isym,iat) )
           arg2=two_pi*( qptn(1)*indsym(1,isym,jat)+&
             qptn(2)*indsym(2,isym,jat)+&
             qptn(3)*indsym(3,isym,jat) )

           re=cos(arg1)*cos(arg2)+sin(arg1)*sin(arg2)
           im=isgn*(cos(arg2)*sin(arg1)-cos(arg1)*sin(arg2))

           do idir=1,3     ! loop over displacements
             do jdir=1,3   ! loop over displacements
               ! we pick the (iat,jat) (3x3) block of the dyn.mat.
               sumr=zero
               sumi=zero
               do ii=1,3
                 do jj=1,3
                   sumr=sumr+symcart(idir,ii,isym)*dynmat(1,ii,niat,jj,njat)*symcart(jdir,jj,isym)
                   sumi=sumi+symcart(idir,ii,isym)*dynmat(2,ii,niat,jj,njat)*symcart(jdir,jj,isym)
                 end do
               end do
               sumi=isgn*sumi

               dynmatint(nsym*(itirev-1)+isym,1,idir,iat,jdir,jat)=re*sumr-im*sumi
               dynmatint(nsym*(itirev-1)+isym,2,idir,iat,jdir,jat)=re*sumi+im*sumr
             end do
           end do
         end do
       end do ! end treatment of the (iat,jat) (3x3) block of dynmat
     end if ! symmetry check
   end do ! time-reversal
 end do ! symmetries

 !4) make the average, get the final symmetric dynamical matrix
 do iat=1,natom
   do jat=1,natom
     do idir=1,3
       do jdir=1,3
         dmint=zero
         do isym=1,2*nsym
           dmint=dmint+dynmatint(isym,1,idir,iat,jdir,jat)
         end do
         dynmat(1,idir,iat,jdir,jat)=dmint/dble(indij(iat,jat))
         dmint=zero
         do isym=1,2*nsym
           dmint=dmint+dynmatint(isym,2,idir,iat,jdir,jat)
         end do
         dynmat(2,idir,iat,jdir,jat)=dmint/dble(indij(iat,jat))
       end do
     end do
   end do
 end do

 i1=0
 do iat=1,natom
   do idir=1,3
     i1=i1+1
     i2=0
     do jat=1,natom
       do jdir=1,3
         i2=i2+1
         index=i1+3*natom*(i2-1)
         dmati(2*index-1)=dynmat(1,idir,iat,jdir,jat)
         dmati(2*index  )=dynmat(2,idir,iat,jdir,jat)
       end do
     end do
   end do
 end do

end subroutine symdyma
!!***

!!****f* m_dynmat/dfpt_sygra
!!
!! NAME
!! dfpt_sygra
!!
!! FUNCTION
!! Symmetrize derivatives of energy with respect to coordinates,
!! as appearing in phonon calculations.
!! Unsymmetrized gradients are input as deunsy; symmetrized grads are then placed in desym.
!! If nsym=1 simply copy deunsy into desym (only symmetry is identity).
!! The index of the initial perturbation is needed, in case there is a change
!! of atom position (moved in another cell) due to the symmetry operation.
!!
!! INPUTS
!!  natom=number of atoms in cell
!!  deunsy(2,3,natom)=unsymmetrized gradients wrt dimensionless tn (hartree)
!!  note: there is a real and a imaginary part ...
!!  indsym(4,nsym,natom)=label given by subroutine symatm, indicating atom
!!   label which gets rotated into given atom by given symmetry
!!   (first three elements are related primitive translation--
!!   see symatm where this is computed)
!!  nsym=number of symmetry operators in group
!!  ipert=index of the initial perturbation
!!  qpt(3)= wavevector of the phonon, in reduced coordinates
!!  symrec(3,3,nsym)=symmetries of group in terms of operations on
!!    reciprocal space primitive translations--see comments below
!!
!! OUTPUT
!! desym(2,3,natom)=symmetrized gradients wrt dimensionless tn (hartree)
!!
!! NOTES
!! Written by X. Gonze starting from sygrad, written by D.C. Allan:
!!    introduction of the q vector for phonon symmetrization
!! This routine should once be merged with sygrad...
!!
!! PARENTS
!!      dfpt_nstdy,dfpt_nstpaw
!!
!! CHILDREN
!!
!! SOURCE

subroutine dfpt_sygra(natom,desym,deunsy,indsym,ipert,nsym,qpt,symrec)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: ipert,natom,nsym
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symrec(3,3,nsym)
 real(dp),intent(in) :: deunsy(2,3,natom),qpt(3)
 real(dp),intent(out) :: desym(2,3,natom)

!Local variables -------------------------
!scalars
 integer :: ia,ind,isym,mu
 real(dp) :: arg,im,re,sumi,sumr

! *********************************************************************

!DEBUG
!write(std_out,*)' dfpt_sygra : enter '
!write(std_out,*)' dfpt_sygra : qpt(:)',qpt(:)
!do ia=1,natom
!do mu=1,3
!write(std_out,*)' dfpt_sygra : deunsy(:2,mu,ia)',deunsy(:2,mu,ia)
!enddo
!enddo
!ENDDEBUG

 if (nsym==1) then

!  Only symmetry is identity so simply copy
   desym(:,:,:)=deunsy(:,:,:)

 else

!  Actually conduct symmetrization
!  write(std_out,*)' dfpt_sygra : desym(:2,:3,:natom),qpt(:)',desym(:2,:3,:natom),qpt(:)
   do ia=1,natom
!    write(std_out,*)' dfpt_sygra : ia=',ia
     do mu=1,3
       sumr=zero
       sumi=zero
!      write(std_out,*)' dfpt_sygra : mu=',mu
       do isym=1,nsym
         ind=indsym(4,isym,ia)
!        Must shift the atoms back to the unit cell.
!        arg=two_pi*( qpt(1)*indsym(1,isym,ia)&
!        &         +qpt(2)*indsym(2,isym,ia)&
!        &         +qpt(3)*indsym(3,isym,ia) )
!        Selection of non-zero q point, to avoid ipert being outside the 1 ... natom range
         if(qpt(1)**2+qpt(2)**2+qpt(3)**2 > tol16)then
           arg=two_pi*( qpt(1)*(indsym(1,isym,ia)-indsym(1,isym,ipert))&
&           +qpt(2)* (indsym(2,isym,ia)-indsym(2,isym,ipert))&
&           +qpt(3)* (indsym(3,isym,ia)-indsym(3,isym,ipert)))
         else
           arg=zero
         end if

         re=dble(symrec(mu,1,isym))*deunsy(1,1,ind)+&
&         dble(symrec(mu,2,isym))*deunsy(1,2,ind)+&
&         dble(symrec(mu,3,isym))*deunsy(1,3,ind)
         im=dble(symrec(mu,1,isym))*deunsy(2,1,ind)+&
&         dble(symrec(mu,2,isym))*deunsy(2,2,ind)+&
&         dble(symrec(mu,3,isym))*deunsy(2,3,ind)
         sumr=sumr+re*cos(arg)-im*sin(arg)
         sumi=sumi+re*sin(arg)+im*cos(arg)
!        sumr=sumr+re
!        sumi=sumi+im
!        write(std_out,*)' dfpt_sygra : isym,indsym(4,isym,ia),arg,re,im,sumr,sumi',&
!        &      isym,indsym(4,isym,ia),arg,re,im,sumr,sumi
       end do
       desym(1,mu,ia)=sumr/dble(nsym)
       desym(2,mu,ia)=sumi/dble(nsym)
!      write(std_out,*)' dfpt_sygra : desym(:,mu,ia)',desym(:,mu,ia)
     end do
   end do
 end if

end subroutine dfpt_sygra
!!***

!!****f* m_dynmat/dfpt_sydy
!! NAME
!! dfpt_sydy
!!
!! FUNCTION
!! Symmetrize dynamical matrix (eventually diagonal wrt to the atoms)
!! Unsymmetrized dynamical matrix is input as dyfrow;
!! symmetrized dynamical matrix is then placed in sdyfro.
!! If nsym=1 simply copy dyfrow into sdyfro.
!!
!! INPUTS
!!  cplex=1 if dynamical matrix is real, 2 if it is complex
!!  dyfrow(3,3,natom,1+(natom-1)*nondiag)=unsymmetrized dynamical matrix
!!  indsym(4,msym*natom)=indirect indexing array: for each
!!   isym,iatom, fourth element is label of atom into which iatom is sent by
!!   INVERSE of symmetry operation isym; first three elements are the primitive
!!   translations which must be subtracted after the transformation to get back
!!   to the original unit cell.
!!  natom=number of atoms in cell.
!!  nondiag=0 if dynamical matrix is     diagonal with respect to atoms
!           1 if dynamical matrix is non diagonal with respect to atoms
!!  nsym=number of symmetry operators in group.
!!  qphon(3)=wavevector of the phonon
!!  symq(4,2,nsym)=1 if symmetry preserves present qpoint. From littlegroup_q
!!  symrec(3,3,nsym)=symmetries of group in terms of operations on real
!!    space primitive translations (integers).
!!
!! OUTPUT
!!  sdyfro(3,3,natom,1+(natom-1)*nondiag)=symmetrized dynamical matrix
!!
!! NOTES
!! Symmetrization of gradients with respect to reduced
!! coordinates tn is conducted according to the expression
!! $[d(e)/d(t(n,a))]_{symmetrized} = (1/Nsym)*Sum(S)*symrec(n,m,S)*
!!              [d(e)/d(t(m,b))]_{unsymmetrized}$
!! where $t(m,b)= (symrel^{-1})(m,n)*(t(n,a)-tnons(n))$ and tnons
!! is a possible nonsymmorphic translation.  The label "b" here
!! refers to the atom which gets rotated into "a" under symmetry "S".
!! symrel is the symmetry matrix in real space, which is the inverse
!! transpose of symrec.  symrec is the symmetry matrix in reciprocal
!! space.  $sym_{cartesian} = R * symrel * R^{-1} = G * symrec * G^{-1}$
!! where the columns of R and G are the dimensional primitive translations
!! in real and reciprocal space respectively.
!! Note the use of "symrec" in the symmetrization expression above.
!!
!! PARENTS
!!      dfpt_dyfro
!!
!! CHILDREN
!!
!! SOURCE

subroutine dfpt_sydy(cplex,dyfrow,indsym,natom,nondiag,nsym,qphon,sdyfro,symq,symrec)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: cplex,natom,nondiag,nsym
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symq(4,2,nsym),symrec(3,3,nsym)
 real(dp),intent(in) :: dyfrow(cplex,3,3,natom,1+(natom-1)*nondiag),qphon(3)
 real(dp),intent(out) :: sdyfro(cplex,3,3,natom,1+(natom-1)*nondiag)

!Local variables -------------------------
!scalars
 integer :: ia,indi,indj,isym,ja,kappa,mu,natom_nondiag,nsym_used,nu
 logical :: qeq0
 real(dp) :: arg,div,phasei,phaser
!arrays
 real(dp) :: work(cplex,3,3)

! *********************************************************************

 if (nsym==1) then

!  Only symmetry is identity so simply copy
   sdyfro(:,:,:,:,:)=dyfrow(:,:,:,:,:)

 else

!  Actually carry out symmetrization
   sdyfro(:,:,:,:,:)=zero
   qeq0=(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-14)
!  === Diagonal dyn. matrix OR q=0
   if (nondiag==0.or.qeq0) then
     natom_nondiag=1;if (nondiag==1) natom_nondiag=natom
     do ja=1,natom_nondiag
       do ia=1,natom
         do isym=1,nsym
           indi=indsym(4,isym,ia)
           indj=1;if (nondiag==1) indj=indsym(4,isym,ja)
           work(:,:,:)=zero
           do mu=1,3
             do nu=1,3
               do kappa=1,3
                 work(:,mu,kappa)=work(:,mu,kappa)+symrec(mu,nu,isym)*dyfrow(:,nu,kappa,indi,indj)
               end do
             end do
           end do
           do mu=1,3
             do nu=1,3
               do kappa=1,3
                 sdyfro(:,kappa,mu,ia,ja)=sdyfro(:,kappa,mu,ia,ja)+symrec(mu,nu,isym)*work(:,kappa,nu)
               end do
             end do
           end do
         end do
       end do
     end do
     div=one/dble(nsym)
     sdyfro(:,:,:,:,:)=div*sdyfro(:,:,:,:,:)
!    === Non diagonal dyn. matrix AND q<>0
   else
     do ja=1,natom
       do ia=1,natom
         nsym_used=0
         do isym=1,nsym
           if (symq(4,1,isym)==1) then
             arg=two_pi*(qphon(1)*(indsym(1,isym,ia)-indsym(1,isym,ja)) &
&             +qphon(2)*(indsym(2,isym,ia)-indsym(2,isym,ja)) &
&             +qphon(3)*(indsym(3,isym,ia)-indsym(3,isym,ja)))
             phaser=cos(arg);phasei=sin(arg)
             nsym_used=nsym_used+1
             indi=indsym(4,isym,ia)
             indj=indsym(4,isym,ja)
             work(:,:,:)=zero
             do mu=1,3
               do nu=1,3
                 do kappa=1,3
                   work(:,mu,kappa)=work(:,mu,kappa)+symrec(mu,nu,isym)*dyfrow(:,nu,kappa,indi,indj)
                 end do
               end do
             end do
             do mu=1,3
               do nu=1,3
                 do kappa=1,3
                   sdyfro(1,kappa,mu,ia,ja)=sdyfro(1,kappa,mu,ia,ja) &
&                   +symrec(mu,nu,isym)*(work(1,kappa,nu)*phaser-work(2,kappa,nu)*phasei)
                 end do
               end do
             end do
             if (cplex==2) then
               do mu=1,3
                 do nu=1,3
                   do kappa=1,3
                     sdyfro(2,kappa,mu,ia,ja)=sdyfro(2,kappa,mu,ia,ja) &
&                     +symrec(mu,nu,isym)*(work(1,kappa,nu)*phasei+work(2,kappa,nu)*phaser)
                   end do
                 end do
               end do
             end if
           end if
         end do
         div=one/dble(nsym_used)
         sdyfro(:,:,:,ia,ja)=div*sdyfro(:,:,:,ia,ja)
       end do
     end do
   end if

 end if

end subroutine dfpt_sydy
!!***

!    CODE TO BE EVENTUALLY REUSED
!    Sym preserves direction and atom
!    if (symq(1,1,isym)==0.and.symq(2,1,isym)==0.and.symq(3,1,isym)==0.and.symq(4,1,isym)==1)then
!      if (ipert==indsym(4,isym,ipert)) then
!        tok=1
!        do idir1=1,3
!          if ((idir1==idir.and.symrec(idir,idir1,isym)/=1).or.&
! &            (idir1/=idir.and.symrec(idir,idir1,isym)/=0)) tok=0
!        end do
!      end if
!    end if
!    div=one/dble(count(symq(4,1,:)==1))

!----------------------------------------------------------------------

!!****f* m_dynmat/wings3
!! NAME
!! wings3
!!
!! FUNCTION
!!  Suppress the wings of the cartesian 2DTE for which
!!  the diagonal element is not known
!!
!! INPUTS
!!  carflg(3,mpert,3,mpert)= ( 1 if the element of the cartesian
!!  2DTE matrix has been calculated correctly ; 0 otherwise )
!!  d2cart(2,3,mpert,3,mpert)=
!!   dynamical matrix, effective charges, dielectric tensor,....
!!   all in cartesian coordinates
!!  mpert =maximum number of ipert
!!
!! OUTPUT
!!  d2cart(2,3,mpert,3,mpert) without the wings
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!
!! SOURCE

subroutine wings3(carflg,d2cart,mpert)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert
!arrays
 integer,intent(inout) :: carflg(3,mpert,3,mpert)
 real(dp),intent(inout) :: d2cart(2,3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer :: idir,idir1,ipert,ipert1

! *********************************************************************

 do ipert=1,mpert
   do idir=1,3
     if(carflg(idir,ipert,idir,ipert)==0)then
       do ipert1=1,mpert
         do idir1=1,3
           carflg(idir,ipert,idir1,ipert1)=0
           carflg(idir1,ipert1,idir,ipert)=0
           d2cart(1,idir,ipert,idir1,ipert1)=zero
           d2cart(2,idir,ipert,idir1,ipert1)=zero
           d2cart(1,idir1,ipert1,idir,ipert)=zero
           d2cart(2,idir1,ipert1,idir,ipert)=zero
         end do
       end do
     end if
   end do
 end do

end subroutine wings3
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/asrif9
!!
!! NAME
!! asrif9
!!
!! FUNCTION
!! Imposes the Acoustic Sum Rule to Interatomic Forces
!!
!! INPUTS
!! asr= Option for the imposition of the ASR
!!  0 => no ASR,
!!  1 => modify "asymmetrically" the diagonal element
!!  2 => modify "symmetrically" the diagonal element
!! natom= Number of atoms in the unit cell
!! nrpt= Number of R points in the Big Box
!! rpt(3,nprt)= Canonical coordinates of the R points in the unit cell
!!  These coordinates are normalized (=> * acell(3)!!)
!! wghatm(natom,natom,nrpt)= Weight associated to the couple of atoms and the R vector
!! atmfrc(3,natom,3,natom,nrpt)= Interatomic Forces
!!
!! OUTPUT
!! atmfrc(3,natom,3,natom,nrpt)= ASR-imposed Interatomic Forces
!!
!! TODO
!! List of ouput should be included.
!!
!! PARENTS
!!      ddb_hybrid,m_ifc,m_tdep_abitypes
!!
!! CHILDREN
!!
!! SOURCE

subroutine asrif9(asr,atmfrc,natom,nrpt,rpt,wghatm)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: asr,natom,nrpt
!arrays
 real(dp),intent(in) :: rpt(3,nrpt),wghatm(natom,natom,nrpt)
 real(dp),intent(inout) :: atmfrc(3,natom,3,natom,nrpt)

!Local variables -------------------------
!scalars
 integer :: found,ia,ib,irpt,izero,mu,nu
 real(dp) :: sumifc

! *********************************************************************

 if(asr==1.or.asr==2)then
   found=0
   ! Search for the R vector which is equal to ( 0 , 0 , 0 )
   ! This vector leaves the atom a on itself !
   do irpt=1,nrpt
     if (all(abs(rpt(:,irpt))<=1.0d-10)) then
       found=1
       izero=irpt
     end if
     if (found==1) exit
   end do

   if(found==0)then
     MSG_BUG('Not able to find the vector R=(0,0,0).')
   end if

   do mu=1,3
     do nu=1,3
       do ia=1,natom
         sumifc=zero
         do ib=1,natom

           ! Get the sumifc of interatomic forces acting on the atom ia,
           ! either in a symmetrical manner, or an unsymmetrical one.
           if(asr==1)then
             do irpt=1,nrpt
               sumifc=sumifc+wghatm(ia,ib,irpt)*atmfrc(mu,ia,nu,ib,irpt)
             end do
           else if(asr==2)then
             do irpt=1,nrpt
               sumifc=sumifc+&
                (wghatm(ia,ib,irpt)*atmfrc(mu,ia,nu,ib,irpt)+&
                 wghatm(ia,ib,irpt)*atmfrc(nu,ia,mu,ib,irpt))/2
             end do
           end if
         end do

         ! Correct the self-interaction in order to fulfill the ASR
         atmfrc(mu,ia,nu,ia,izero)=atmfrc(mu,ia,nu,ia,izero)-sumifc
         if (asr==2) atmfrc(nu,ia,mu,ia,izero)=atmfrc(mu,ia,nu,ia,izero)
       end do
     end do
   end do
 end if

end subroutine asrif9
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/get_bigbox_and_weights
!! NAME
!! get_bigbox_and_weights
!!
!! FUNCTION
!! Compute the Big Box containing the R points in the cartesian real space needed to Fourier Transform
!! the dynamical matrix into its corresponding interatomic force.
!!
!! INPUTS
!! brav= Bravais Lattice (1 or -1=S.C.;2=F.C.C.;3=BCC;4=Hex.)
!! natom= Number of atoms
!! nqbz= Number of q-points in BZ.
!! ngqpt(3)= Numbers used to generate the q points to sample the Brillouin zone using an homogeneous grid
!! nqshift= number of shifts in q-mesh
!! qshift(3, nqshift) = Q-mesh shifts
!! rprim(3,3)= Normalized coordinates in real space.
!! rprimd, gprimd
!! rcan(3,natom)  = Atomic position in canonical coordinates
!! cutmode=Define the cutoff used to filter the output R-points according to their weights.
!!   0 --> No cutoff (mainly for debugging)
!!   1 --> Include only those R-points for which sum(abs(wg(:,:,irpt)) < tol20
!!         This is the approach used for the dynamical matrix.
!!   2 --> Include only those R-points for which the trace over iatom of abs(wg(iat,iat,irpt)) < tol20
!!         This option is used for objects that depend on a single atomic index.
!! comm= MPI communicator
!!
!! OUTPUT
!! nrpt= Total Number of R points in the Big Box
!! cell(3,nrpt) Give the index of the the cell and irpt
!! rpt(3,nrpt)= Canonical coordinates of the R points in the unit cell. These coordinates are normalized (=> * acell(3)!!)
!! r_inscribed_sphere
!! wghatm(natom,natom,nrpt)= Weights associated to a pair of atoms and to a R vector
!!
!! PARENTS
!!      m_ifc
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_bigbox_and_weights(brav, natom, nqbz, ngqpt, nqshift, qshift, rprim, rprimd, gprim, rcan, &
                                  cutmode, nrpt, rpt, cell, wghatm, r_inscribed_sphere, comm)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: brav, natom, nqbz, nqshift, cutmode, comm
 integer,intent(out) :: nrpt
 real(dp),intent(out) :: r_inscribed_sphere
!arrays
 integer,intent(in) :: ngqpt(3)
 real(dp),intent(in) :: gprim(3,3),rprim(3,3),rprimd(3,3), rcan(3, natom)
 real(dp),intent(in) :: qshift(3, nqshift)
 integer,allocatable,intent(out) :: cell(:,:)
 real(dp),allocatable,intent(out) :: rpt(:,:), wghatm(:,:,:)

!Local variables -------------------------
!scalars
 integer :: my_ierr, ierr, ii, irpt, all_nrpt
 real(dp) :: toldist
 integer :: ngqpt9(9)
 character(len=500*4) :: msg
!arrays
 integer,allocatable :: all_cell(:,:)
 real(dp),allocatable :: all_rpt(:,:), all_wghatm(:,:,:)

! *********************************************************************

 ABI_CHECK(any(cutmode == [0, 1, 2]), "cutmode should be in [0, 1, 2]")

 ! Create the Big Box of R vectors in real space and compute the number of points (cells) in real space
 call make_bigbox(brav, all_cell, ngqpt, nqshift, rprim, all_nrpt, all_rpt)

 ! Weights associated to these R points and to atomic pairs
 ABI_MALLOC(all_wghatm, (natom, natom, all_nrpt))

 ! HM: this tolerance is highly dependent on the compilation/architecture
 !     numeric errors in the DDB text file. Try a few tolerances and check whether all the weights are found.
 ngqpt9 = 0; ngqpt9(1:3) = ngqpt(1:3)
 toldist = tol8
 do while (toldist <= tol6)
   ! Note ngqpt(9) with intent(inout)!
   call wght9(brav, gprim, natom, ngqpt9, nqbz, nqshift, all_nrpt, qshift, rcan, &
              all_rpt, rprimd, toldist, r_inscribed_sphere, all_wghatm, my_ierr)
   call xmpi_max(my_ierr, ierr, comm, ii)
   if (ierr > 0) toldist = toldist * 10
   if (ierr == 0) exit
 end do

 if (ierr > 0) then
   write(msg, '(3a,es14.4,2a,i0, 14a)' ) &
    'The sum of the weight is not equal to nqpt.',ch10,&
    'The sum of the weights is: ',sum(all_wghatm),ch10,&
    'The number of q points is: ',nqbz, ch10, &
    'This might have several sources.',ch10,&
    'If toldist is larger than 1.0e-8, the atom positions might be loose.',ch10,&
    'and the q point weights not computed properly.',ch10,&
    'Action: make input atomic positions more symmetric.',ch10,&
    'Otherwise, you might increase "buffer" in m_dynmat.F90 see bigbx9 subroutine and recompile.',ch10,&
    'Actually, this can also happen when ngqpt is 0 0 0,',ch10,&
    'if abs(brav) /= 1, in this case you should change brav to 1. If brav is already set to 1 (default) try -1.'
   MSG_ERROR(msg)
 end if

 ! Only conserve the necessary points in rpt.
 nrpt = 0
 do irpt=1,all_nrpt
   if (filterw(all_wghatm(:,:,irpt))) cycle
   nrpt = nrpt + 1
 end do

 ! Allocate output arrays and transfer data.
 ABI_MALLOC(rpt, (3, nrpt))
 ABI_MALLOC(cell, (3, nrpt))
 ABI_MALLOC(wghatm, (natom, natom, nrpt))

 ii = 0
 do irpt=1,all_nrpt
   if (filterw(all_wghatm(:,:,irpt))) cycle
   ii = ii + 1
   rpt(:, ii) = all_rpt(:,irpt)
   wghatm(:,:,ii) = all_wghatm(:,:,irpt)
   cell(:,ii) = all_cell(:,irpt)
 end do

 ABI_FREE(all_rpt)
 ABI_FREE(all_wghatm)
 ABI_FREE(all_cell)

contains

logical pure function filterw(wg)

 real(dp),intent(in) :: wg(natom,natom)
 integer :: iat
 real(dp) :: trace

 select case (cutmode)
 case (1)
   filterw = sum(abs(wg)) < tol20
 case (2)
   trace = zero
   do iat=1,natom
     trace = trace + abs(wg(iat,iat))
   end do
   filterw = trace < tol20
 case default
   filterw = .False.
 end select

end function filterw

end subroutine get_bigbox_and_weights
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/make_bigbox
!! NAME
!! make_bigbox
!!
!! FUNCTION
!! Helper functions to faciliate the generation of a Big Box containing
!! all the R points in the cartesian real space needed to Fourier Transform
!! the dynamical matrix into its corresponding interatomic force.
!! See bigbx9 for the algorithm.
!!
!! INPUTS
!! brav= Bravais Lattice (1 or -1=S.C.;2=F.C.C.;3=BCC;4=Hex.)
!! ngqpt(3)= Numbers used to generate the q points to sample the
!!   Brillouin zone using an homogeneous grid
!! nqshft= number of q-points in the repeated cell for the Brillouin zone sampling
!!  When nqshft is not 1, but 2 or 4 (only other allowed values),
!!  the limits for the big box have to be extended by a factor of 2.
!! rprim(3,3)= Normalized coordinates in real space  !!! IS THIS CORRECT?
!!
!! OUTPUT
!! cell= (3,nrpt) Give the index of the the cell and irpt
!! nprt= Number of R points in the Big Box
!! rpt(3,nrpt)= Canonical coordinates of the R points in the unit cell
!!  These coordinates are normalized (=> * acell(3)!!)
!!  The array is allocated here with the proper dimension. Client code is responsible
!!  for the deallocation.
!!
!! PARENTS
!!      m_ifc
!!
!! CHILDREN
!!
!! SOURCE

subroutine make_bigbox(brav, cell, ngqpt, nqshft, rprim, nrpt, rpt)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: brav,nqshft
 integer,intent(out) :: nrpt
!arrays
 integer,intent(in) :: ngqpt(3)
 real(dp),intent(in) :: rprim(3,3)
 real(dp),allocatable,intent(out) :: rpt(:,:)
 integer,allocatable,intent(out) :: cell(:,:)

!Local variables -------------------------
!scalars
 integer :: choice,mrpt
!arrays
 real(dp) :: dummy_rpt(3,1)
 integer:: dummy_cell(1,3)

! *********************************************************************

 ! Compute the number of points (cells) in real space
 choice=0
 call bigbx9(brav,dummy_cell,choice,1,ngqpt,nqshft,mrpt,rprim,dummy_rpt)

 ! Now we can allocate and calculate the points and the weights.
 nrpt = mrpt
 ABI_MALLOC(rpt,(3,nrpt))
 ABI_MALLOC(cell,(3,nrpt))

 choice=1
 call bigbx9(brav,cell,choice,mrpt,ngqpt,nqshft,nrpt,rprim,rpt)

end subroutine make_bigbox
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/bigbx9
!! NAME
!! bigbx9
!!
!! FUNCTION
!! Generation of a Big Box containing all the R points in the
!! cartesian real space needed to Fourier Transforms the dynamical
!! matrix into its corresponding interatomic force.
!!
!! INPUTS
!! brav= Bravais Lattice (1 or -1=S.C.;2=F.C.C.;3=BCC;4=Hex.)
!! choice= if 0, simply count nrpt ; if 1, checks that the input mrpt
!!   is the same as nrpt, and generate rpt(3,mrpt)
!! mrpt=dimension of rpt
!! ngqpt(3)= Numbers used to generate the q points to sample the
!!  Brillouin zone using an homogeneous grid
!! nqshft= number of q-points in the repeated cell for the Brillouin zone sampling
!!  When nqshft is not 1, but 2 or 4 (only other allowed values),
!!  the limits for the big box have to be extended by a factor of 2.
!! rprim(3,3)= Normalized coordinates in real space  !!! IS THIS CORRECT?
!!
!! OUTPUT
!! cell= (3,nrpt) Give the index of the the cell and irpt
!! nprt= Number of R points in the Big Box
!! rpt(3,mrpt)= Canonical coordinates of the R points in the unit cell
!!  These coordinates are normalized (=> * acell(3)!!)
!!  (output only if choice=1)
!!
!! PARENTS
!!      m_dynmat
!!
!! CHILDREN
!!
!! SOURCE

subroutine bigbx9(brav,cell,choice,mrpt,ngqpt,nqshft,nrpt,rprim,rpt)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: brav,choice,mrpt,nqshft
 integer,intent(out) :: nrpt
!arrays
 integer,intent(in) :: ngqpt(3)
 real(dp),intent(in) :: rprim(3,3)
 real(dp),intent(out) :: rpt(3,mrpt)
 integer,intent(out) :: cell(3,mrpt)

!Local variables -------------------------
!In some cases, the atoms coordinates are not packed in the
! [0,1]^3 cube. Then, the parameter "buffer" might be increased,
!to search relevant pairs of atoms in bigger boxes than usual.
!scalars
 integer,parameter :: buffer=1
 integer :: irpt,lim1,lim2,lim3,lqshft,r1,r2,r3
 character(len=500) :: msg

! *********************************************************************

 lqshft=1
 if(nqshft/=1)lqshft=2


!Simple Cubic Lattice
 if (abs(brav)==1) then
   lim1=((ngqpt(1))+1)*lqshft+buffer
   lim2=((ngqpt(2))+1)*lqshft+buffer
   lim3=((ngqpt(3))+1)*lqshft+buffer
   nrpt=(2*lim1+1)*(2*lim2+1)*(2*lim3+1)
   if(choice/=0)then
     if (nrpt/=mrpt) then
       write(msg,'(2(a,i0))')' nrpt=',nrpt,' is not equal to mrpt= ',mrpt
       MSG_BUG(msg)
     end if
     irpt=0
     do r1=-lim1,lim1
       do r2=-lim2,lim2
         do r3=-lim3,lim3
           irpt=irpt+1
           rpt(1,irpt)=r1*rprim(1,1)+r2*rprim(1,2)+r3*rprim(1,3)
           rpt(2,irpt)=r1*rprim(2,1)+r2*rprim(2,2)+r3*rprim(2,3)
           rpt(3,irpt)=r1*rprim(3,1)+r2*rprim(3,2)+r3*rprim(3,3)
           cell(1,irpt)=r1;cell(2,irpt)=r2;cell(3,irpt)=r3
         end do
       end do
     end do
   end if

!  Face Centered Cubic Lattice
 else if (brav==2) then
   lim1=((ngqpt(1)+3)/4)*lqshft+buffer
   lim2=((ngqpt(2)+3)/4)*lqshft+buffer
   lim3=((ngqpt(3)+3)/4)*lqshft+buffer
   nrpt=(2*lim1+1)*(2*lim2+1)*(2*lim3+1)*4
   if(choice/=0)then
     if (nrpt/=mrpt) then
       write(msg,'(2(a,i0))')' nrpt=',nrpt,' is not equal to mrpt= ',mrpt
       MSG_BUG(msg)
     end if
     irpt=0
     do r1=-lim1,lim1
       do r2=-lim2,lim2
         do r3=-lim3,lim3
           irpt=irpt+4
           rpt(1,irpt-3)=r1
           rpt(2,irpt-3)=r2
           rpt(3,irpt-3)=r3
           rpt(1,irpt-2)=r1
           rpt(2,irpt-2)=r2+0.5
           rpt(3,irpt-2)=r3+0.5
           rpt(1,irpt-1)=r1+0.5
           rpt(2,irpt-1)=r2
           rpt(3,irpt-1)=r3+0.5
           rpt(1,irpt)=r1+0.5
           rpt(2,irpt)=r2+0.5
           rpt(3,irpt)=r3
!TEST_AM
!           cell(irpt-3,1)=r1;cell(irpt-3,2)=r2;cell(irpt-3,3)=r3
           cell(1,irpt)=r1;cell(2,irpt)=r2;cell(3,irpt)=r3
         end do
       end do
     end do
   end if

!  Body Centered Cubic Lattice
 else if (brav==3) then
   lim1=((ngqpt(1)+3)/4)*lqshft+buffer
   lim2=((ngqpt(2)+3)/4)*lqshft+buffer
   lim3=((ngqpt(3)+3)/4)*lqshft+buffer
   nrpt=(2*lim1+1)*(2*lim2+1)*(2*lim3+1)*2
   if(choice/=0)then
     if(nrpt/=mrpt) then
       write(msg,'(2(a,i0))')' nrpt= ',nrpt,' is not equal to mrpt= ',mrpt
       MSG_BUG(msg)
     end if
     irpt=0
     do r1=-lim1,lim1
       do r2=-lim2,lim2
         do r3=-lim3,lim3
           irpt=irpt+2
           rpt(1,irpt-1)=r1
           rpt(2,irpt-1)=r2
           rpt(3,irpt-1)=r3
           rpt(1,irpt)=r1+0.5
           rpt(2,irpt)=r2+0.5
           rpt(3,irpt)=r3+0.5
!TEST_AM
!           cell(irpt-1,1)=r1;cell(irpt-1,2)=r2;cell(irpt-1,3)=r3
           cell(1,irpt)=r1;cell(2,irpt)=r2;cell(3,irpt)=r3
         end do
       end do
     end do
   end if

!  Hexagonal Lattice
 else if (brav==4) then
   lim1=(ngqpt(1)+1)*lqshft+buffer
   lim2=(ngqpt(2)+1)*lqshft+buffer
   lim3=((ngqpt(3)/2)+1)*lqshft+buffer
   nrpt=(2*lim1+1)*(2*lim2+1)*(2*lim3+1)
   if(choice/=0)then
     if(nrpt/=mrpt)then
       write(msg,'(2(a,i0))')' nrpt=',nrpt,' is not equal to mrpt=',mrpt
       MSG_BUG(msg)
     end if
     irpt=0
     do r1=-lim1,lim1
       do r2=-lim2,lim2
         do r3=-lim3,lim3
           irpt=irpt+1
           rpt(1,irpt)=r1*rprim(1,1)+r2*rprim(1,2)+r3*rprim(1,3)
           rpt(2,irpt)=r1*rprim(2,1)+r2*rprim(2,2)+r3*rprim(2,3)
           rpt(3,irpt)=r1*rprim(3,1)+r2*rprim(3,2)+r3*rprim(3,3)
           cell(1,irpt)=r1;cell(2,irpt)=r2;cell(3,irpt)=r3
         end do
       end do
     end do
   end if

 else
   write(msg,'(a,i0,a)')' The value of brav= ',brav,' is not allowed (should be -1, 1, 2 or 4).'
   MSG_BUG(msg)
 end if

end subroutine bigbx9
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/canat9
!! NAME
!! canat9
!!
!! FUNCTION
!! Transforms an atom whose coordinates (xred*rprim) would not be
!! in the chosen unit cell used to generate the interatomic forces
!! to its correspondent (rcan) in canonical coordinates.
!!
!! INPUTS
!! brav= Bravais Lattice (1 or -1=S.C.;2=F.C.C.;3=BCC;4=Hex.)
!! natom= Number of atoms in the unit cell
!! rprim(3,3)= Normalized coordinates  of primitive vectors
!!
!! OUTPUT
!! rcan(3,natom)  = Atomic position in canonical coordinates
!! trans(3,natom) = Atomic translations : xred = rcan + trans
!!
!! PARENTS
!!      m_ifc
!!
!! CHILDREN
!!
!! SOURCE


subroutine canat9(brav,natom,rcan,rprim,trans,xred)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: brav,natom
!arrays
 real(dp),intent(in) :: rprim(3,3),xred(3,natom)
 real(dp),intent(out) :: rcan(3,natom),trans(3,natom)

!Local variables -------------------------
!scalars
 integer :: found,iatom,ii
 character(len=500) :: msg
!arrays
 real(dp) :: dontno(3,4),rec(3),rok(3),shift(3),tt(3)

! *********************************************************************

!Normalization of the cartesian atomic coordinates
!If not normalized : rcan(i) <- rcan(i) * acell(i)
 do iatom=1,natom
   rcan(:,iatom)=xred(1,iatom)*rprim(:,1)+xred(2,iatom)*rprim(:,2)+xred(3,iatom)*rprim(:,3)
 end do

 !Study of the different cases for the Bravais lattice:
 if (abs(brav)==1) then
   !Simple Cubic Lattice

   do iatom=1,natom
     ! Canon will produces these coordinate transformations
     ! (Note: here we still use reduced coordinates )
     call wrap2_pmhalf(xred(1,iatom),rok(1),shift(1))
     call wrap2_pmhalf(xred(2,iatom),rok(2),shift(2))
     call wrap2_pmhalf(xred(3,iatom),rok(3),shift(3))

!    New coordinates : rcan
     rcan(:,iatom)=rok(1)*rprim(:,1)+rok(2)*rprim(:,2)+rok(3)*rprim(:,3)
!    Translations between New and Old coordinates
     tt(:)=xred(1,iatom)*rprim(:,1)+xred(2,iatom)*rprim(:,2)+xred(3,iatom)*rprim(:,3)
     trans(:,iatom)=tt(:)-rcan(:,iatom)
   end do

 else if (brav==2) then
   ! Face Centered Lattice
   ! Special possible translations in the F.C.C. case
   dontno(:,:)=zero
   dontno(2,2)=0.5_dp
   dontno(3,2)=0.5_dp
   dontno(1,3)=0.5_dp
   dontno(3,3)=0.5_dp
   dontno(1,4)=0.5_dp
   dontno(2,4)=0.5_dp
   do iatom=1,natom
     found=0
     do ii=1,4
       if (found==1) exit
       ! Canon will produces these coordinate transformations
       call wrap2_pmhalf(rcan(1,iatom)+dontno(1,ii),rok(1),shift(1))
       call wrap2_pmhalf(rcan(2,iatom)+dontno(2,ii),rok(2),shift(2))
       call wrap2_pmhalf(rcan(3,iatom)+dontno(3,ii),rok(3),shift(3))
       ! In the F.C.C., ABS[ Ri ] + ABS[ Rj ] < or = 1/2
       ! The equal signs has been treated using a tolerance parameter
       ! not to have twice the same point in the unit cell !
       rok(1)=rok(1)-1.0d-10
       rok(2)=rok(2)-2.0d-10
       rok(3)=rok(3)-5.0d-10
       if (abs(rok(1))+abs(rok(2))<=0.5_dp) then
         if (abs(rok(1))+abs(rok(3))<=0.5_dp) then
           if (abs(rok(2))+abs(rok(3))<=0.5_dp) then
             tt(:)=rcan(:,iatom)
             ! New coordinates : rcan
             rcan(1,iatom)=rok(1)+1.0d-10
             rcan(2,iatom)=rok(2)+2.0d-10
             rcan(3,iatom)=rok(3)+5.0d-10
             ! Translations between New and Old coordinates
             trans(:,iatom)=tt(:)-rcan(:,iatom)
             found=1
           end if
         end if
       end if
     end do
   end do

 else if (brav==3) then
   ! Body Centered Cubic Lattice
   ! Special possible translations in the B.C.C. case
   dontno(:,1)=zero
   dontno(:,2)=0.5_dp
   do iatom=1,natom
     found=0
     do ii=1,2
       if (found==1) exit
       ! Canon will produce these coordinate transformations
       call wrap2_pmhalf(rcan(1,iatom)+dontno(1,ii),rok(1),shift(1))
       call wrap2_pmhalf(rcan(2,iatom)+dontno(2,ii),rok(2),shift(2))
       call wrap2_pmhalf(rcan(3,iatom)+dontno(3,ii),rok(3),shift(3))
       ! In the F.C.C., ABS[ Ri ] < or = 1/2
       ! and    ABS[ R1 ] + ABS[ R2 ] + ABS[ R3 ] < or = 3/4
       ! The equal signs has been treated using a tolerance parameter
       ! not to have twice the same point in the unit cell !
       rok(1)=rok(1)-1.0d-10
       rok(2)=rok(2)-2.0d-10
       rok(3)=rok(3)-5.0d-10
       if(abs(rok(1))+abs(rok(2))+abs(rok(3))<=0.75_dp)then
         if ( abs(rok(1))<=0.5_dp .and. abs(rok(2))<=0.5_dp .and. abs(rok(3))<=0.5_dp) then
           tt(:)=rcan(:,iatom)
           ! New coordinates : rcan
           rcan(1,iatom)=rok(1)+1.0d-10
           rcan(2,iatom)=rok(2)+2.0d-10
           rcan(3,iatom)=rok(3)+5.0d-10
           ! Translations between New and Old coordinates
           trans(:,iatom)=tt(:)-rcan(:,iatom)
           found=1
         end if
       end if
     end do
   end do

 else if (brav==4) then
   ! Hexagonal Lattice
   ! In this case, it is easier first to work in reduced coordinates space !
   do iatom=1,natom
     ! Passage from the reduced space to the "lozenge" cell
     rec(1)=xred(1,iatom)-0.5_dp
     rec(2)=xred(2,iatom)-0.5_dp
     rec(3)=xred(3,iatom)
     ! Canon will produces these coordinate transformations
     call wrap2_pmhalf(rec(1),rok(1),shift(1))
     call wrap2_pmhalf(rec(2),rok(2),shift(2))
     call wrap2_pmhalf(rec(3),rok(3),shift(3))
     rec(1)=rok(1)+0.5_dp
     rec(2)=rok(2)+0.5_dp
     rec(3)=rok(3)
     ! Passage in Cartesian Normalized Coordinates
     rcan(:,iatom)=rec(1)*rprim(:,1)+rec(2)*rprim(:,2)+rec(3)*rprim(:,3)
     ! Use of a tolerance parameter not to have twice the same pointin the unit cell !
     rcan(1,iatom)=rcan(1,iatom)-1.0d-10
     rcan(2,iatom)=rcan(2,iatom)-2.0d-10
     ! Passage to the honey-com hexagonal unit cell !
     if (rcan(1,iatom)>0.5_dp) then
       rcan(1,iatom)=rcan(1,iatom)-1.0_dp
     end if
     if (rcan(1,iatom)>zero.and.rcan(1,iatom)+sqrt(3.0_dp)*rcan(2,iatom)>1.0_dp) then
       rcan(1,iatom)=rcan(1,iatom)-0.5_dp
       rcan(2,iatom)=rcan(2,iatom)-sqrt(3.0_dp)*0.5_dp
     end if
     if (rcan(1,iatom)<=zero.and.sqrt(3.0_dp)*rcan(2,iatom)-rcan(1,iatom)>1.0_dp) then
       rcan(1,iatom)=rcan(1,iatom)+0.5_dp
       rcan(2,iatom)=rcan(2,iatom)-sqrt(3.0_dp)*0.5_dp
     end if
     ! Translations between New and Old coordinates
     tt(:)=xred(1,iatom)*rprim(:,1)+xred(2,iatom)*rprim(:,2)+xred(3,iatom)*rprim(:,3)
     trans(:,iatom)=tt(:)-rcan(:,iatom)
   end do

   ! End of the possible cases for brav : -1, 1, 2, 4.
 else
   write(msg, '(a,i0,a,a,a)' )&
   'The required value of brav=',brav,' is not available.',ch10,&
   'It should be -1, 1,2 or 4 .'
   MSG_BUG(msg)
 end if

 call wrtout(std_out,' Canonical Atomic Coordinates ')
 do iatom=1,natom
   write(msg, '(a,i5,3es18.8)' )' atom',iatom,rcan(1,iatom),rcan(2,iatom),rcan(3,iatom)
   call wrtout(std_out,msg)
 end do

end subroutine canat9
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/canct9
!!
!! NAME
!! canct9
!!
!! FUNCTION
!! Convert from canonical coordinates to cartesian coordinates
!! a vector defined by its index=ib+natom*(irpt-1)
!!
!! INPUTS
!! acell(3)=length scales by which rprim is to be multiplied
!! gprim(3,3)=dimensionless primitive translations in reciprocal space
!! index= index of the atom
!! natom=number of atoms in unit cell
!! nrpt= Number of R points in the Big Box
!! rcan(3,natom)=canonical coordinates of atoms
!! rprim(3,3)=dimensionless primitive translations in real space
!! rpt(3,nrpt)=canonical coordinates of the points in the BigBox.
!!
!! OUTPUT
!! ib=number of the atom in the unit cell
!! irpt= number of the unit cell to which belong the atom
!! rcart(3)=cartesian coordinate of the atom indexed by index.
!!
!! PARENTS
!!      ddb_hybrid,m_ifc
!!
!! CHILDREN
!!
!! SOURCE

subroutine canct9(acell,gprim,ib,index,irpt,natom,nrpt,rcan,rcart,rprim,rpt)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: index,natom,nrpt
 integer,intent(out) :: ib,irpt
!arrays
 real(dp),intent(in) :: acell(3),gprim(3,3),rcan(3,natom),rprim(3,3)
 real(dp),intent(in) :: rpt(3,nrpt)
 real(dp),intent(out) :: rcart(3)

!Local variables -------------------------
!scalars
 integer :: jj
!arrays
 real(dp) :: xred(3)

! *********************************************************************

 irpt=(index-1)/natom+1
 ib=index-natom*(irpt-1)

!Transform the canonical coordinates to reduced coord.
 do jj=1,3
   xred(jj)=gprim(1,jj)*(rpt(1,irpt)+rcan(1,ib))&
&   +gprim(2,jj)*(rpt(2,irpt)+rcan(2,ib))&
&   +gprim(3,jj)*(rpt(3,irpt)+rcan(3,ib))
 end do

!Then to cartesian coordinates (here the position of the atom b)
 do jj=1,3
   rcart(jj)=xred(1)*acell(1)*rprim(jj,1)+&
&   xred(2)*acell(2)*rprim(jj,2)+&
&   xred(3)*acell(3)*rprim(jj,3)
 end do

end subroutine canct9
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/chkrp9
!! NAME
!! chkrp9
!!
!! FUNCTION
!! Check if the rprim used for the definition of the unit cell (in the
!! inputs) are consistent with the rprim used in the routine generating
!! the Big Box needed to generate the interatomic forces.
!!
!! INPUTS
!! brav=bravais lattice (1 or -1=simple lattice,2=face centered lattice,
!!  3=centered lattice,4=hexagonal lattice)
!! rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!
!! OUTPUT
!!  (only checking)
!!
!! PARENTS
!!      m_ifc
!!
!! CHILDREN
!!
!! SOURCE

subroutine chkrp9(brav,rprim)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: brav
!arrays
 real(dp),intent(in) :: rprim(3,3)

!Local variables -------------------------
!scalars
 integer :: ii,jj
 character(len=500) :: msg

! *********************************************************************

 if (abs(brav)==1) then
!  Simple Cubic Lattice No condition in this case !
   continue

 else if (brav==2) then
!  Face Centered Lattice
   do ii=1,3
     do jj=1,3
       if (  ( ii==jj .and. abs(rprim(ii,jj))>tol10) .or. (ii/=jj .and. abs(rprim(ii,jj)-.5_dp)>tol10) ) then
         write(msg, '(a,a,a,a,a,a,a,a,a,a,a)' )&
         'The input variable rprim does not correspond to the',ch10,&
         'fixed rprim to be used with brav=2 and ifcflag=1 :',ch10,&
         '   0  1/2  1/2',ch10,&
         '  1/2  0   1/2',ch10,&
         '  1/2 1/2   0 ',ch10,&
         'Action: rebuild your DDB by using the latter rprim.'
         MSG_ERROR(msg)
       end if
     end do
   end do

 else if (brav==3) then
!  Body Centered Cubic Lattice
   do ii=1,3
     do jj=1,3
       if (  ( ii==jj .and. abs(rprim(ii,jj)+.5_dp)>tol10) .or. (ii/=jj .and. abs(rprim(ii,jj)-.5_dp)>tol10) ) then
         write(msg, '(a,a,a,a,a,a,a,a,a,a,a)' )&
         'The input variable rprim does not correspond to the',ch10,&
         'fixed rprim to be used with brav=3 and ifcflag=1 :',ch10,&
         '  -1/2  1/2  1/2',ch10,&
         '   1/2 -1/2  1/2',ch10,&
         '   1/2  1/2 -1/2',ch10,&
         'Action: rebuild your DDB by using the latter rprim.'
         MSG_ERROR(msg)
       end if
     end do
   end do

 else if (brav==4) then
!  Hexagonal Lattice
   if (abs(rprim(1,1)-1.0_dp)>tol10 .or. &
       abs(rprim(3,3)-1.0_dp)>tol10 .or. &
       abs(rprim(2,1)      )>tol10 .or. &
       abs(rprim(3,1)      )>tol10 .or. &
       abs(rprim(1,3)      )>tol10 .or. &
       abs(rprim(2,3)      )>tol10 .or. &
       abs(rprim(3,2)      )>tol10 .or. &
       abs(rprim(1,2)+0.5_dp)>tol10 .or. &
       abs(rprim(2,2)-0.5_dp*sqrt(3.0_dp))>tol10 ) then
     write(msg, '(a,a,a,a,a,a,a,a,a,a,a)' )&
      'The input variable rprim does not correspond to the',ch10,&
      'fixed rprim to be used with brav=4 and ifcflag=1 :',ch10,&
      '   1      0      0',ch10,&
      '  -1/2 sqrt[3]/2 0',ch10,&
      '   0      0      1',ch10,&
      'Action: rebuild your DDB by using the latter rprim.'
     MSG_ERROR(msg)
   end if

 else
   write(msg, '(a,i4,a,a,a,a,a)' )&
   'The value of brav=',brav,' is not allowed.',ch10,&
   'Only  -1, 1,2,3 or 4 are allowed.',ch10,&
   'Action: change the value of brav in your input file.'
   MSG_ERROR(msg)
 end if

end subroutine chkrp9
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/dist9
!! NAME
!! dist9
!!
!! FUNCTION
!! Compute the distance between atoms
!!
!! INPUTS
!! acell(3)=length scales by which rprim is to be multiplied
!! dist(natom,natom,nrpt)=distances between atoms
!! gprim(3,3)=dimensionless primitive translations in reciprocal space
!! natom=number of atoms in unit cell
!! nrpt= Number of R points in the Big Box
!! rcan(3,natom)=canonical coordinates of atoms
!! rprim(3,3)=dimensionless primitive translations in real space
!! rpt(3,nrpt)=cartesian coordinates of the points in the BigBox.
!!
!! OUTPUT
!! dist(natom,natom,nrpt)=distances between atoms
!!
!! PARENTS
!!      m_ifc
!!
!! CHILDREN
!!
!! SOURCE

subroutine dist9(acell,dist,gprim,natom,nrpt,rcan,rprim,rpt)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom,nrpt
!arrays
 real(dp),intent(in) :: acell(3),gprim(3,3),rcan(3,natom),rprim(3,3)
 real(dp),intent(in) :: rpt(3,nrpt)
 real(dp),intent(out) :: dist(natom,natom,nrpt)

!Local variables -------------------------
!scalars
 integer :: ia,ib,ii,irpt
!arrays
 real(dp) :: ra(3),rb(3),rdiff(3),red(3),rptcar(3),xred(3)

! *********************************************************************

!BIG loop on all generic atoms
 do ia=1,natom
   ! First transform canonical coordinates to reduced coordinates
   do ii=1,3
     xred(ii)=gprim(1,ii)*rcan(1,ia)+gprim(2,ii)*rcan(2,ia)+gprim(3,ii)*rcan(3,ia)
   end do
   ! Then to cartesian coordinates
   ra(:)=xred(1)*acell(1)*rprim(:,1)+xred(2)*acell(2)*rprim(:,2)+xred(3)*acell(3)*rprim(:,3)
   do ib=1,natom
     do ii=1,3
       xred(ii)=gprim(1,ii)*rcan(1,ib)+gprim(2,ii)*rcan(2,ib)+gprim(3,ii)*rcan(3,ib)
     end do
     do ii=1,3
       rb(ii)=xred(1)*acell(1)*rprim(ii,1)+xred(2)*acell(2)*rprim(ii,2)+xred(3)*acell(3)*rprim(ii,3)
     end do
     do irpt=1,nrpt
       ! First transform it to reduced coordinates
       do ii=1,3
         red(ii)=gprim(1,ii)*rpt(1,irpt)+gprim(2,ii)*rpt(2,irpt)+gprim(3,ii)*rpt(3,irpt)
       end do
       ! Then to cartesian coordinates
       do ii=1,3
         rptcar(ii)=red(1)*acell(1)*rprim(ii,1)+red(2)*acell(2)*rprim(ii,2)+red(3)*acell(3)*rprim(ii,3)
       end do
       do ii=1,3
         rdiff(ii)=-rptcar(ii)+ra(ii)-rb(ii)
       end do
       dist(ia,ib,irpt)=(rdiff(1)**2+rdiff(2)**2+rdiff(3)**2)**0.5
     end do
   end do
 end do

end subroutine dist9
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/ftifc_q2r
!!
!! NAME
!! ftifc_q2r
!!
!! FUNCTION
!!  Generates the Fourier transform of the dynamical matrices
!!  to obtain interatomic forces (real space).
!!
!! INPUTS
!! dynmat(2,3,natom,3,natom,nqpt)= Dynamical matrices coming from the Derivative Data Base
!! gprim(3,3)= Normalized coordinates in reciprocal space
!! natom= Number of atoms in the unit cell
!! nqpt= Number of q points in the Brillouin zone
!! nrpt= Number of R points in the Big Box
!! rpt(3,nprt)= Canonical coordinates of the R points in the unit cell
!!           These coordinates are normalized (=> * acell(3)!!)
!! spqpt(3,nqpt)= Reduced coordinates of the q vectors in reciprocal space
!! comm=MPI communicator.
!!
!! OUTPUT
!! atmfrc(3,natom,3,natom,nrpt)= Interatomic Forces in real space.
!!  We used the imaginary part just for debugging!
!!
!! PARENTS
!!      m_ifc
!!
!! CHILDREN
!!
!! SOURCE

subroutine ftifc_q2r(atmfrc,dynmat,gprim,natom,nqpt,nrpt,rpt,spqpt,comm)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom,nqpt,nrpt,comm
!arrays
 real(dp),intent(in) :: gprim(3,3),rpt(3,nrpt),spqpt(3,nqpt)
 real(dp),intent(out) :: atmfrc(3,natom,3,natom,nrpt)
 real(dp),intent(in) :: dynmat(2,3,natom,3,natom,nqpt)

!Local variables -------------------------
!scalars
 integer :: ia,ib,iqpt,irpt,mu,nu,nprocs,my_rank,ierr
 real(dp) :: im,kr,re
!arrays
 real(dp) :: kk(3)

! *********************************************************************

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 ! Interatomic Forces from Dynamical Matrices
 atmfrc = zero
 do irpt=1,nrpt
   if (mod(irpt, nprocs) /= my_rank) cycle ! mpi-parallelism
   do iqpt=1,nqpt

     ! Calculation of the k coordinates in Normalized Reciprocal coordinates
     kk(1)=spqpt(1,iqpt)*gprim(1,1)+spqpt(2,iqpt)*gprim(1,2)+spqpt(3,iqpt)*gprim(1,3)
     kk(2)=spqpt(1,iqpt)*gprim(2,1)+spqpt(2,iqpt)*gprim(2,2)+spqpt(3,iqpt)*gprim(2,3)
     kk(3)=spqpt(1,iqpt)*gprim(3,1)+spqpt(2,iqpt)*gprim(3,2)+spqpt(3,iqpt)*gprim(3,3)

     ! Product of k and r
     kr=kk(1)*rpt(1,irpt)+kk(2)*rpt(2,irpt)+kk(3)*rpt(3,irpt)

     ! Get the phase factor
     re=cos(two_pi*kr)
     im=sin(two_pi*kr)

     ! Now, big inner loops on atoms and directions
     ! The indices are ordered to give better speed
     do ib=1,natom
       do nu=1,3
         do ia=1,natom
           do mu=1,3
             ! Real and imaginary part of the interatomic forces
             atmfrc(mu,ia,nu,ib,irpt)=atmfrc(mu,ia,nu,ib,irpt) &
              +re*dynmat(1,mu,ia,nu,ib,iqpt)&
              +im*dynmat(2,mu,ia,nu,ib,iqpt)
             !The imaginary part should be equal to zero !!!!!!
             !atmfrc(2,mu,ia,nu,ib,irpt)=atmfrc(2,mu,ia,nu,ib,irpt) &
             !          +re*dynmat(2,mu,ia,nu,ib,iqpt) &
             !          -im*dynmat(1,mu,ia,nu,ib,iqpt)
           end do
         end do
       end do
     end do

   end do
 end do

 call xmpi_sum(atmfrc, comm, ierr)
 !The sumifc has to be weighted by a normalization factor of 1/nqpt
 atmfrc = atmfrc/nqpt

end subroutine ftifc_q2r
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/ftifc_r2q
!!
!! NAME
!! ftifc_r2q
!!
!! FUNCTION
!! Generates the Fourier transform of the interatomic forces
!! to obtain dynamical matrices in reciprocal space: R --> q.
!!
!! INPUTS
!! atmfrc(3,natom,3,natom,nrpt)= Interatomic Forces in real space
!! gprim(3,3)= Normalized coordinates in reciprocal space
!! natom= Number of atoms in the unit cell
!! nqpt= Number of q points in the Brillouin zone
!! nrpt= Number of R points in the Big Box
!! rpt(3,nprt)= Canonical coordinates of the R points in the unit cell
!!   These coordinates are normalized (=> * acell(3)!!)
!! spqpt(3,nqpt)= Reduced coordinates of the q vectors in reciprocal space
!! wghatm(natom,natom,nrpt)= Weights associated to a pair of atoms and to a R vector
!! comm: MPI communicator
!!
!! OUTPUT
!! dynmat(2,3,natom,3,natom,nqpt)= Dynamical matrices coming from the Derivative Data Base
!!
!! PARENTS
!!      m_dynmat
!!
!! CHILDREN
!!
!! SOURCE

subroutine ftifc_r2q(atmfrc, dynmat, gprim, natom, nqpt, nrpt, rpt, spqpt, wghatm, comm)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom,nqpt,nrpt,comm
!arrays
 real(dp),intent(in) :: gprim(3,3),rpt(3,nrpt),spqpt(3,nqpt)
 real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 real(dp),intent(in) :: atmfrc(3,natom,3,natom,nrpt)
 real(dp),intent(out) :: dynmat(2,3,natom,3,natom,nqpt)

!Local variables -------------------------
!scalars
 integer :: ia,ib,iqpt,irpt,mu,nu,cnt,my_rank,nprocs, ierr
 real(dp) :: facti,factr,im,kr,re
!arrays
 real(dp) :: kk(3)

! *********************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 dynmat = zero; cnt = 0

 ! MG: This is an hotspot. I don'tknow whether one should rewrite with BLAS1 dot or not.
 ! Note, however, that simply removing the check on the weights inside the loop over atoms.
 ! leads to a non-negligible speedup with intel (~30% if dipdip -1 is used)
 do iqpt=1,nqpt

   ! Calculation of the k coordinates in Normalized Reciprocal coordinates
   kk(1)=spqpt(1,iqpt)*gprim(1,1)+spqpt(2,iqpt)* gprim(1,2)+spqpt(3,iqpt)*gprim(1,3)
   kk(2)=spqpt(1,iqpt)*gprim(2,1)+spqpt(2,iqpt)* gprim(2,2)+spqpt(3,iqpt)*gprim(2,3)
   kk(3)=spqpt(1,iqpt)*gprim(3,1)+spqpt(2,iqpt)* gprim(3,2)+spqpt(3,iqpt)*gprim(3,3)

   do irpt=1,nrpt
     cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle ! MPI parallelism.

     ! k.R
     kr = kk(1)*rpt(1,irpt)+kk(2)*rpt(2,irpt)+kk(3)*rpt(3,irpt)
     ! Get phase factor
     re = cos(two_pi*kr); im = sin(two_pi*kr)

     ! Inner loop on atoms and directions
     do ib=1,natom
       do ia=1,natom
         !if (abs(wghatm(ia,ib,irpt)) > tol10) then  ! Commented by MG
           factr = re * wghatm(ia,ib,irpt)
           facti = im * wghatm(ia,ib,irpt)
           do nu=1,3
             do mu=1,3
               ! Real and imaginary part of the dynamical matrices
               ! Atmfrc should be real
               dynmat(1,mu,ia,nu,ib,iqpt) = dynmat(1,mu,ia,nu,ib,iqpt) + factr * atmfrc(mu,ia,nu,ib,irpt)
               dynmat(2,mu,ia,nu,ib,iqpt) = dynmat(2,mu,ia,nu,ib,iqpt) + facti * atmfrc(mu,ia,nu,ib,irpt)
             end do
           end do
         !end if
       end do
     end do
   end do
 end do

 if (nprocs > 1) call xmpi_sum(dynmat, comm, ierr)

end subroutine ftifc_r2q
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/dynmat_dq
!!
!! NAME
!! dynmat_dq
!!
!! FUNCTION
!!  Compute the derivative D(q)/dq of the dynamical matrix via Fourier transform
!!  of the interatomic forces
!!
!! INPUTS
!! qpt(3)= Reduced coordinates of the q vector in reciprocal space
!! natom= Number of atoms in the unit cell
!! gprim(3,3)= Normalized coordinates in reciprocal space
!! nrpt= Number of R points in the Big Box
!! rpt(3,nprt)= Canonical coordinates of the R points in the unit cell
!!   These coordinates are normalized (=> * acell(3)!!)
!! atmfrc(3,natom,3,natom,nrpt)= Interatomic Forces in real space
!! wghatm(natom,natom,nrpt)= Weights associated to a pair of atoms and to a R vector
!!
!! OUTPUT
!! dddq(2,3,natom,3,natom,3)= Derivate of the dynamical matrix in cartesian coordinates.
!!  The tree directions are stored in the last dimension.
!!  These coordinates are normalized (=> * acell(3)!!)
!!
!! PARENTS
!!      m_ifc
!!
!! CHILDREN
!!
!! SOURCE

subroutine dynmat_dq(qpt,natom,gprim,nrpt,rpt,atmfrc,wghatm,dddq)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom,nrpt
!arrays
 real(dp),intent(in) :: gprim(3,3),rpt(3,nrpt),qpt(3)
 real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 real(dp),intent(in) :: atmfrc(3,natom,3,natom,nrpt)
 real(dp),intent(out) :: dddq(2,3,natom,3,natom,3)

!Local variables -------------------------
!scalars
 integer :: ia,ib,irpt,mu,nu,ii
 real(dp) :: im,kr,re
!arrays
 real(dp) :: kk(3),fact(2,3)

! *********************************************************************

 dddq = zero

 do irpt=1,nrpt
   ! Calculation of the k coordinates in Normalized Reciprocal coordinates
   kk(1) = qpt(1)*gprim(1,1)+ qpt(2)*gprim(1,2) + qpt(3)*gprim(1,3)
   kk(2) = qpt(1)*gprim(2,1)+ qpt(2)*gprim(2,2) + qpt(3)*gprim(2,3)
   kk(3) = qpt(1)*gprim(3,1)+ qpt(2)*gprim(3,2) + qpt(3)*gprim(3,3)

   ! Product of k and r
   kr=kk(1)*rpt(1,irpt)+kk(2)*rpt(2,irpt)+kk(3)*rpt(3,irpt)

   ! Get phase factor
   re=cos(two_pi*kr); im=sin(two_pi*kr)

   ! Inner loop on atoms and directions
   do ib=1,natom
     do ia=1,natom
       if (abs(wghatm(ia,ib,irpt))>1.0d-10) then
         ! take into account rotation due to i.
         fact(1,:) = -im * wghatm(ia,ib,irpt) * rpt(:,irpt)
         fact(2,:) =  re * wghatm(ia,ib,irpt) * rpt(:,irpt)
         do nu=1,3
           do mu=1,3
             ! Real and imaginary part of the dynamical matrices
             ! Atmfrc should be real
             do ii=1,3
               dddq(1,mu,ia,nu,ib,ii) = dddq(1,mu,ia,nu,ib,ii) + fact(1,ii) * atmfrc(mu,ia,nu,ib,irpt)
               dddq(2,mu,ia,nu,ib,ii) = dddq(2,mu,ia,nu,ib,ii) + fact(2,ii) * atmfrc(mu,ia,nu,ib,irpt)
             end do
           end do
         end do
       end if
     end do
   end do
 end do

end subroutine dynmat_dq
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/ifclo9
!! NAME
!! ifclo9
!!
!! FUNCTION
!! Convert from cartesian coordinates to local coordinates
!! the 3*3 interatomic force constant matrix
!!
!! INPUTS
!! ifccar(3,3)= matrix of interatomic force constants in cartesian
!!  coordinates
!! vect1(3)= cartesian coordinates of the first local vector
!! vect2(3)= cartesian coordinates of the second local vector
!! vect3(3)= cartesian coordinates of the third local vector
!!
!! OUTPUT
!! ifcloc(3,3)= matrix of interatomic force constants in local coordinates
!!
!! PARENTS
!!      m_ifc
!!
!! CHILDREN
!!
!! SOURCE

subroutine ifclo9(ifccar,ifcloc,vect1,vect2,vect3)

!Arguments -------------------------------
!arrays
 real(dp),intent(in) :: ifccar(3,3),vect1(3),vect2(3),vect3(3)
 real(dp),intent(out) :: ifcloc(3,3)

!Local variables -------------------------
!scalars
 integer :: ii,jj
!arrays
 real(dp) :: work(3,3)

! *********************************************************************

 do jj=1,3
   do ii=1,3
     work(jj,ii)=zero
   end do
   do ii=1,3
     work(jj,1)=work(jj,1)+ifccar(jj,ii)*vect1(ii)
     work(jj,2)=work(jj,2)+ifccar(jj,ii)*vect2(ii)
     work(jj,3)=work(jj,3)+ifccar(jj,ii)*vect3(ii)
   end do
 end do

 do jj=1,3
   do ii=1,3
     ifcloc(ii,jj)=zero
   end do
   do ii=1,3
     ifcloc(1,jj)=ifcloc(1,jj)+vect1(ii)*work(ii,jj)
     ifcloc(2,jj)=ifcloc(2,jj)+vect2(ii)*work(ii,jj)
     ifcloc(3,jj)=ifcloc(3,jj)+vect3(ii)*work(ii,jj)
   end do
 end do

end subroutine ifclo9
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/wght9
!! NAME
!! wght9
!!
!! FUNCTION
!! Generates a weight to each R points of the Big Box and for each pair of atoms
!! For each R points included in the space generates by moving
!! the unit cell around each atom; the weight will be one.
!! Border conditions are provided.
!! The R points outside the chosen space will have a 0 weight.
!!
!! INPUTS
!! brav = Bravais lattice (1 or -1=S.C.;2=F.C.C.;4=Hex. -1 is for old algo to find weights, =1 is for Wigner-Seitz algo)
!! gprim(3,3)= Normalized coordinates in reciprocal space
!! natom= Number of atoms in the unit cell
!! ngqpt(6)= Numbers used to sample the Brillouin zone
!! nqpt= Number of q points used in the homogeneous grid
!!  sampling the Brillouin zone
!! nqshft=number of shift vectors in the repeated cell
!! nrpt=Number of R points in the Big Box
!! qshft(3,nqshft)=vectors that will be used to determine
!!  the shifts from (0. 0. 0.)
!! rcan(3,natom)=Atomic position in canonical coordinates
!! rpt(3,nprt)=Canonical coordinates of the R points in the unit cell
!!  These coordinates are normalized (=> * acell(3))
!! rprimd(3,3)=dimensional primitive translations for real space (bohr)
!! toldist= Tolerance on the distance between two R points.
!!
!! OUTPUT
!! wghatm(natom,natom,nrpt)= Weight associated to the couple of atoms and the R vector
!!  The vector r(atom2)-r(atom1)+rpt should be inside the moving box
!! ngqpt(6)= can be modified
!!
!! PARENTS
!!      m_ifc
!!
!! CHILDREN
!!
!! SOURCE

subroutine wght9(brav,gprim,natom,ngqpt,nqpt,nqshft,nrpt,qshft,rcan,rpt,rprimd,toldist,r_inscribed_sphere,wghatm,ierr)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: brav,natom,nqpt,nqshft,nrpt
 integer,intent(out) :: ierr
 real(dp),intent(out) :: r_inscribed_sphere
 real(dp),intent(in) :: toldist
!arrays
 integer,intent(inout) :: ngqpt(9)
 real(dp),intent(in) :: gprim(3,3),qshft(3,4),rcan(3,natom),rpt(3,nrpt),rprimd(3,3)
 real(dp),intent(out) :: wghatm(natom,natom,nrpt)

!Local variables -------------------------
!scalars
 integer :: ia,ib,ii,jj,kk,iqshft,irpt,jqshft,nbordh,tok,nptws,nreq,idir
 real(dp) :: factor,sumwght,normsq,proj
 character(len=500) :: msg
!arrays
 integer :: nbord(9)
 real(dp) :: rdiff(9),red(3,3),ptws(4, 729),pp(3),rdiff_tmp(3)

! *********************************************************************

 ierr = 0

 ! First analyze the vectors qshft
 if (nqshft /= 1) then

   if (brav == 4) then
     write(msg,'(3a,i0,3a)' )&
     'For the time being, only nqshft=1',ch10,&
     'is allowed with brav=4, while it is nqshft=',nqshft,'.',ch10,&
     'Action: in the input file, correct either brav or nqshft.'
     MSG_ERROR(msg)
   end if

   if (nqshft == 2) then
     ! Make sure that the q vectors form a BCC lattice
     do ii=1,3
       if(abs(abs(qshft(ii,1)-qshft(ii,2))-.5_dp)>1.d-10)then
         write(msg, '(a,a,a,a,a,a,a)' )&
         'The test of the q1shft vectors shows that they',ch10,&
         'do not generate a body-centered lattice, which',ch10,&
         'is mandatory for nqshft=2.',ch10,&
         'Action: change the q1shft vectors in your input file.'
         MSG_ERROR(msg)
       end if
     end do
   else if (nqshft == 4) then
     ! Make sure that the q vectors form a FCC lattice
     do iqshft=1,3
       do jqshft=iqshft+1,4
         tok=0
         do ii=1,3
           ! Test on the presence of a +-0.5 difference
           if(abs(abs(qshft(ii,iqshft)-qshft(ii,jqshft))-.5_dp) <1.d-10) tok=tok+1
           ! Test on the presence of a 0 or +-1.0 difference
           if(abs(abs(qshft(ii,iqshft)-qshft(ii,jqshft))-1._dp) <1.d-10  .or.&
              abs(qshft(ii,iqshft)-qshft(ii,jqshft)) < 1.d-10) tok=tok+4
         end do
         ! Test 1 should be satisfied twice, and test 2 once
         if(tok/=6)then
           write(msg, '(7a)' )&
           'The test of the q1shft vectors shows that they',ch10,&
           'do not generate a face-centered lattice, which',ch10,&
           'is mandatory for nqshft=4.',ch10,&
           'Action: change the q1shft vectors in your input file.'
           MSG_ERROR(msg)
         end if
       end do
     end do
   else
     write(msg, '(a,i4,3a)' )&
     'nqshft must be 1, 2 or 4. It is nqshft=',nqshft,'.',ch10,&
     'Action: change nqshft in your input file.'
     MSG_ERROR(msg)
   end if
 end if

 factor=0.5_dp
 if(brav==2 .or. brav==3) factor=0.25_dp
 if(nqshft/=1)factor=factor*2

 if (brav==1) then
   ! Does not support multiple shifts
   if (nqshft/=1) then
     MSG_ERROR('This version of the weights does not support nqshft/=1.')
   end if

   ! Find the points of the lattice given by ngqpt*acell. These are used to define
   ! a Wigner-Seitz cell around the origin. The origin is excluded from the list.
   ! TODO : in principle this should be only -1 to +1 for ii jj kk!
   nptws=0
   do ii=-2,2
     do jj=-2,2
       do kk=-2,2
         do idir=1,3
           pp(idir)=ii*ngqpt(1)*rprimd(idir,1)+ jj*ngqpt(2)*rprimd(idir,2)+ kk*ngqpt(3)*rprimd(idir,3)
         end do
         normsq = pp(1)*pp(1)+pp(2)*pp(2)+pp(3)*pp(3)
         if (normsq > tol6) then
           nptws = nptws + 1
           ptws(:3,nptws) = pp(:)
           ptws(4,nptws) = half*normsq
         end if
       end do
     end do
   end do
 end if ! end new_wght
 !write(std_out,*)'factor,ngqpt',factor,ngqpt(1:3)

 r_inscribed_sphere = sum((matmul(rprimd(:,:),ngqpt(1:3)))**2)
 do ii=-1,1
   do jj=-1,1
     do kk=-1,1
       if (ii==0 .and. jj==0 .and. kk==0) cycle
       do idir=1,3
         pp(idir)=ii*ngqpt(1)*rprimd(idir,1)+ jj*ngqpt(2)*rprimd(idir,2)+ kk*ngqpt(3)*rprimd(idir,3)
       end do
       normsq = pp(1)*pp(1)+pp(2)*pp(2)+pp(3)*pp(3)
       r_inscribed_sphere = min(r_inscribed_sphere, normsq)
     end do
   end do
 end do
 r_inscribed_sphere = sqrt(r_inscribed_sphere)


!Begin the big loop on ia and ib
 do ia=1,natom
   do ib=1,natom

     ! Simple Lattice
     if (abs(brav)==1) then
       ! In this case, it is better to work in reduced coordinates
       ! As rcan is in canonical coordinates, => multiplication by gprim
       do ii=1,3
         red(1,ii)=  rcan(1,ia)*gprim(1,ii) +rcan(2,ia)*gprim(2,ii) +rcan(3,ia)*gprim(3,ii)
         red(2,ii)=  rcan(1,ib)*gprim(1,ii) +rcan(2,ib)*gprim(2,ii) +rcan(3,ib)*gprim(3,ii)
       end do
     end if

     do irpt=1,nrpt

       ! Initialization of the weights to 1.0
       wghatm(ia,ib,irpt)=1.0_dp

       ! Compute the difference vector

       ! Simple Cubic Lattice
       if (abs(brav)==1) then
         ! Change of rpt to reduced coordinates
         do ii=1,3
           red(3,ii)=  rpt(1,irpt)*gprim(1,ii) +rpt(2,irpt)*gprim(2,ii) +rpt(3,irpt)*gprim(3,ii)
           rdiff(ii)=red(2,ii)-red(1,ii)+red(3,ii)
         end do
         if (brav==1) then
           ! rdiff in cartesian coordinates
           do ii=1,3
             rdiff_tmp(ii)=rdiff(1)*rprimd(ii,1)+rdiff(2)*rprimd(ii,2)+rdiff(3)*rprimd(ii,3)
           end do
           rdiff(1:3)=rdiff_tmp(1:3)
         end if

       else
         ! Other lattices
         do ii=1,3
           rdiff(ii)=rcan(ii,ib)-rcan(ii,ia)+rpt(ii,irpt)
         end do
       end if

       ! Assignement of weights

       if(nqshft==1 .and. brav/=4)then

         if (brav/=1) then
           do ii=1,3
             ! If the rpt vector is greater than the allowed space => weight = 0.0
             if (abs(rdiff(ii))-tol10>factor*ngqpt(ii)) then
               wghatm(ia,ib,irpt)=zero
             else if (abs(abs(rdiff(ii))-factor*ngqpt(ii)) <=1.0d-10) then
               ! If the point is in a boundary position => weight/2
               wghatm(ia,ib,irpt)=wghatm(ia,ib,irpt)/2
             end if
           end do
         else
           ! new weights
           wghatm(ia,ib,irpt)=zero
           nreq = 1
           do ii=1,nptws
             proj = rdiff(1)*ptws(1,ii)+rdiff(2)*ptws(2,ii)+rdiff(3)*ptws(3,ii)
             ! if rdiff closer to ptws than the origin the weight is zero
             ! if rdiff close to the origin with respect to all the other ptws the weight is 1
             ! if rdiff is equidistant from the origin and N other ptws the weight is 1/(N+1)
             if (proj - ptws(4,ii) > toldist) then
               nreq = 0
               EXIT
             else if(abs(proj-ptws(4,ii)) <= toldist) then
               nreq=nreq+1
             end if
           end do
           if (nreq>0) then
             wghatm(ia,ib,irpt)=one/DBLE(nreq)
           end if
         end if

       else if(brav==4)then
         ! Hexagonal
         ! Examination of the X and Y boundaries in order to form an hexagon
         ! First generate the relevant boundaries
         rdiff(4)=0.5_dp*( rdiff(1)+sqrt(3.0_dp)*rdiff(2) )
         ngqpt(4)=ngqpt(1)
         rdiff(5)=0.5_dp*( rdiff(1)-sqrt(3.0_dp)*rdiff(2) )
         ngqpt(5)=ngqpt(1)

         ! Test the four inequalities
         do ii=1,5
           if(ii/=2)then

             nbord(ii)=0
             ! If the rpt vector is greater than the allowed space => weight = 0.0
             if (abs(rdiff(ii))-1.0d-10>factor*ngqpt(ii)) then
               wghatm(ia,ib,irpt)=zero
             else if (abs(abs(rdiff(ii))-factor*ngqpt(ii)) <=1.0d-10) then
               ! If the point is in a boundary position increment nbord(ii)
               nbord(ii)=1
             end if

           end if
         end do

         ! Computation of weights
         nbordh=nbord(1)+nbord(4)+nbord(5)
         if (nbordh==1) then
           wghatm(ia,ib,irpt)=wghatm(ia,ib,irpt)/2
         else if (nbordh==2) then
           wghatm(ia,ib,irpt)=wghatm(ia,ib,irpt)/3
         else if (nbordh/=0) then
           MSG_BUG('There is a problem of borders and weights (hex).')
         end if
         if (nbord(3)==1)then
           wghatm(ia,ib,irpt)=wghatm(ia,ib,irpt)/2
         end if

       else if(nqshft==2 .and. brav/=4)then

         ! BCC packing of k-points
         ! First, generate the relevant boundaries
         rdiff(4)= rdiff(1)+rdiff(2)
         rdiff(5)= rdiff(1)-rdiff(2)
         rdiff(6)= rdiff(1)+rdiff(3)
         rdiff(7)= rdiff(1)-rdiff(3)
         rdiff(8)= rdiff(3)+rdiff(2)
         rdiff(9)= rdiff(3)-rdiff(2)
         if(ngqpt(2)/=ngqpt(1) .or. ngqpt(3)/=ngqpt(1))then
           write(msg, '(a,a,a,3i6,a,a,a,a)' )&
           'In the BCC case, the three ngqpt numbers ',ch10,&
           '    ',ngqpt(1),ngqpt(2),ngqpt(3),ch10,&
           'should be equal.',ch10,&
           'Action: use identical ngqpt(1:3) in your input file.'
           MSG_ERROR(msg)
         end if
         do ii=4,9
           ngqpt(ii)=ngqpt(1)
         end do

         ! Test the relevant inequalities
         nbord(1)=0
         do ii=4,9
           ! If the rpt vector is greater than the allowed space => weight = 0.0
           if (abs(rdiff(ii))-1.0d-10>factor*ngqpt(ii)) then
             wghatm(ia,ib,irpt)=zero
           else if (abs(abs(rdiff(ii))-factor*ngqpt(ii)) <=1.0d-10) then
             ! If the point is in a boundary position increment nbord(1)
             nbord(1)=nbord(1)+1
           end if
         end do

         ! Computation of weights
         if (nbord(1)==1) then
           wghatm(ia,ib,irpt)=wghatm(ia,ib,irpt)/2
         else if (nbord(1)==2) then
           wghatm(ia,ib,irpt)=wghatm(ia,ib,irpt)/3
         else if (nbord(1)==3) then
           wghatm(ia,ib,irpt)=wghatm(ia,ib,irpt)/4
         else if (nbord(1)==4) then
           wghatm(ia,ib,irpt)=wghatm(ia,ib,irpt)/6
         else if (nbord(1)/=0) then
           MSG_ERROR(' There is a problem of borders and weights (BCC).')
         end if

       else if(nqshft==4 .and. brav/=4)then

         ! FCC packing of k-points
         ! First, generate the relevant boundaries
         rdiff(4)= (rdiff(1)+rdiff(2)+rdiff(3))*2._dp/3._dp
         rdiff(5)= (rdiff(1)-rdiff(2)+rdiff(3))*2._dp/3._dp
         rdiff(6)= (rdiff(1)+rdiff(2)-rdiff(3))*2._dp/3._dp
         rdiff(7)= (rdiff(1)-rdiff(2)-rdiff(3))*2._dp/3._dp
         if(ngqpt(2)/=ngqpt(1) .or. ngqpt(3)/=ngqpt(1))then
           write(msg, '(a,a,a,3i6,a,a,a,a)' )&
           'In the FCC case, the three ngqpt numbers ',ch10,&
           '    ',ngqpt(1),ngqpt(2),ngqpt(3),ch10,&
           'should be equal.',ch10,&
           'Action: use identical ngqpt(1:3) in your input file.'
           MSG_ERROR(msg)
         end if
         do ii=4,7
           ngqpt(ii)=ngqpt(1)
         end do

         ! Test the relevant inequalities
         nbord(1)=0
         do ii=1,7
           ! If the rpt vector is greater than the allowed space => weight = 0.0
           if (abs(rdiff(ii))-1.0d-10>factor*ngqpt(ii)) then
             wghatm(ia,ib,irpt)=zero
             ! If the point is in a boundary position increment nbord(1)
           else if (abs(abs(rdiff(ii))-factor*ngqpt(ii)) <=1.0d-10) then
             nbord(1)=nbord(1)+1
           end if
         end do

         ! Computation of weights
         if (nbord(1)==1) then
           wghatm(ia,ib,irpt)=wghatm(ia,ib,irpt)/2
         else if (nbord(1)==2) then
           wghatm(ia,ib,irpt)=wghatm(ia,ib,irpt)/3
         else if (nbord(1)==3) then
           wghatm(ia,ib,irpt)=wghatm(ia,ib,irpt)/4
         else if (nbord(1)/=0 .and. wghatm(ia,ib,irpt)>1.d-10) then
           ! Interestingly nbord(1)==4 happens for some points outside of the volume
           MSG_BUG(' There is a problem of borders and weights (FCC).')
         end if

       else
         write(msg, '(3a,i0,a)' )&
         'One should not arrive here ... ',ch10,&
         'The value nqshft ',nqshft,' is not available'
         MSG_BUG(msg)
       end if
     end do ! Assignement of weights is done
   end do ! End of the double loop on ia and ib
 end do

 ! Check the results
 do ia=1,natom
   do ib=1,natom
     sumwght=zero
     do irpt=1,nrpt
       ! Check if the sum of the weights is equal to the number of q points
       sumwght=sumwght+wghatm(ia,ib,irpt)
       !write(std_out,'(a,3(i0,1x))' )' atom1, atom2, irpt ; rpt ; wghatm ',ia,ib,irpt
       !write(std_out,'(3es16.6,es18.6)' )rpt(1,irpt),rpt(2,irpt),rpt(3,irpt),wghatm(ia,ib,irpt)
     end do
     if (abs(sumwght-nqpt)>tol10) ierr = 1
   end do
 end do

end subroutine wght9
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/d3sym
!! NAME
!! d3sym
!!
!! FUNCTION
!! Given a set of calculated elements of the 3DTE matrix,
!! build (nearly) all the other matrix elements that can be build using symmetries.
!!
!! INPUTS
!!  indsym(4,nsym,natom)=indirect indexing array : for each
!!   isym,iatom, fourth element is label of atom into which iatom is sent by
!!   INVERSE of symmetry operation isym; first three elements are the primitive
!!   translations which must be subtracted after the transformation to get back
!!   to the original unit cell.
!!  mpert =maximum number of ipert
!!  natom= number of atoms
!!  nsym=number of space group symmetries
!!  symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!!  symrel(3,3,nsym)=3x3 matrices of the group symmetries (real space)
!!
!! SIDE EFFECTS
!!  Input/Output
!!  blkflg(3,mpert,3,mpert,3,mpert)= matrix that indicates if an
!!   element of d3 is available (1 if available, 0 otherwise)
!!  d3(2,3,mpert,3,mpert,3,mpert)= matrix of the 3DTE
!!
!! PARENTS
!!      m_ddb,nonlinear
!!
!! CHILDREN
!!
!! SOURCE

subroutine d3sym(blkflg,d3,indsym,mpert,natom,nsym,symrec,symrel)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,natom,nsym
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symrec(3,3,nsym),symrel(3,3,nsym)
 integer,intent(inout) :: blkflg(3,mpert,3,mpert,3,mpert)
 real(dp),intent(inout) :: d3(2,3,mpert,3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer :: found,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert,idisy1,idisy2,idisy3
 integer :: ipesy1,ipesy2,ipesy3,isym,ithree
 real(dp) :: sumi,sumr
!arrays
 integer :: sym1(3,3),sym2(3,3),sym3(3,3)

! *********************************************************************

!DEBUG
!write(std_out,*)'d3sym : enter'
!do i1dir = 1, 3
!do i2dir = 1, 3
!do i3dir = 1, 3
!write(std_out,*)i1dir,i2dir,i3dir,blkflg(i1dir,natom+2,i2dir,natom+2,i3dir,natom+2)
!end do
!end do
!end do
!stop
!ENDDEBUG

!First, take into account the permutations symmetry of
!(i1pert,i1dir) and (i3pert,i3dir)
 do i1pert = 1, mpert
   do i2pert = 1, mpert
     do i3pert = 1, mpert

       do i1dir = 1, 3
         do i2dir = 1, 3
           do i3dir = 1, 3

             if ((blkflg(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==1).and. &
&             (blkflg(i3dir,i3pert,i2dir,i2pert,i1dir,i1pert)/=1)) then

               d3(:,i3dir,i3pert,i2dir,i2pert,i1dir,i1pert) = &
&              d3(:,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)

               blkflg(i3dir,i3pert,i2dir,i2pert,i1dir,i1pert) = 1

             end if

           end do
         end do
       end do

     end do
   end do
 end do

!Big Big Loop : symmetrize three times, because
!of some cases in which one element is not yet available
!at the first pass, and even at the second one !

 do ithree=1,3

!  Loop over perturbations
   do i1pert = 1, mpert
     do i2pert = 1, mpert
       do i3pert = 1, mpert

         do i1dir = 1, 3
           do i2dir = 1, 3
             do i3dir = 1, 3

!              Will get element (idir1,ipert1,idir2,ipert2)
!              so this element should not yet be present ...
               if(blkflg(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)/=1)then

                 d3(:,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = 0_dp

                 do isym = 1, nsym

                   found = 1

                   if (i1pert <= natom) then
                     ipesy1 = indsym(4,isym,i1pert)
                     sym1(:,:) = symrec(:,:,isym)
                   else if (i1pert == natom + 2) then
                     ipesy1 = i1pert
                     sym1(:,:) = symrel(:,:,isym)
                   else
                     found = 0
                   end if

                   if (i2pert <= natom) then
                     ipesy2 = indsym(4,isym,i2pert)
                     sym2(:,:) = symrec(:,:,isym)
                   else if (i2pert == natom + 2) then
                     ipesy2 = i2pert
                     sym2(:,:) = symrel(:,:,isym)
                   else
                     found = 0
                   end if

                   if (i3pert <= natom) then
                     ipesy3 = indsym(4,isym,i3pert)
                     sym3(:,:) = symrec(:,:,isym)
                   else if (i3pert == natom + 2) then
                     ipesy3 = i3pert
                     sym3(:,:) = symrel(:,:,isym)
                   else
                     found = 0
                   end if

                   sumr = 0_dp ; sumi = 0_dp;
                   do idisy1 = 1, 3
                     do idisy2 = 1, 3
                       do idisy3 = 1, 3

                         if ((sym1(i1dir,idisy1) /=0).and.(sym2(i2dir,idisy2) /=0) &
&                         .and.(sym3(i3dir,idisy3) /=0)) then

                           if (blkflg(idisy1,ipesy1,idisy2,ipesy2,idisy3,ipesy3) == 1) then

                             sumr = sumr + sym1(i1dir,idisy1)*sym2(i2dir,idisy2)*&
&                             sym3(i3dir,idisy3)*d3(1,idisy1,ipesy1,idisy2,ipesy2,idisy3,ipesy3)
                             sumi = sumi + sym1(i1dir,idisy1)*sym2(i2dir,idisy2)*&
&                             sym3(i3dir,idisy3)*d3(2,idisy1,ipesy1,idisy2,ipesy2,idisy3,ipesy3)

                           else

                             found = 0

                           end if

                         end if

                       end do
                     end do
                   end do

                   if (found == 1) then
                     d3(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = sumr
                     d3(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = sumi
                     blkflg(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = 1
                   end if

                 end do  ! isym

               end if  ! blkflg

!              Close loop over perturbations
             end do
           end do
         end do
       end do
     end do
   end do

 end do  ! close loop over ithree

end subroutine d3sym
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/sytens
!!
!! NAME
!! sytens
!!
!! FUNCTION
!! Determines the set of irreductible elements of the non-linear
!! optical susceptibility and Raman tensors
!!
!! INPUTS
!!  indsym(4,nsym,natom)=indirect indexing array described above: for each
!!   isym,iatom, fourth element is label of atom into which iatom is sent by
!!   INVERSE of symmetry operation isym; first three elements are the primitive
!!   translations which must be subtracted after the transformation to get back
!!   to the original unit cell.
!!  mpert =maximum number of ipert
!!  natom= number of atoms
!!  nsym=number of space group symmetries
!!  symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!!  symrel(3,3,nsym)=3x3 matrices of the group symmetries (real space)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  rfpert(3,mpert,3,mpert,3,mpert) = array defining the type of perturbations
!!       that have to be computed
!!    At the input :
!!       1   ->   element has to be computed explicitely
!!    At the output :
!!       1   ->   element has to be computed explicitely
!!      -1   ->   use symmetry operations to obtain the corresponding element
!!      -2   ->   element is zero by symmetry
!!
!! PARENTS
!!      m_ddb,nonlinear,respfn
!!
!! CHILDREN
!!
!! SOURCE

subroutine sytens(indsym,mpert,natom,nsym,rfpert,symrec,symrel)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,natom,nsym
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symrec(3,3,nsym),symrel(3,3,nsym)
 integer,intent(inout) :: rfpert(3,mpert,3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer :: flag,found,i1dir,i1dir_,i1pert,i1pert_,i2dir,i2dir_,i2pert,i2pert_
 integer :: i3dir,i3dir_,i3pert,i3pert_,idisy1,idisy2,idisy3,ipesy1,ipesy2
 integer :: ipesy3,isym
!arrays
 integer :: sym1(3,3),sym2(3,3),sym3(3,3)
 integer,allocatable :: pertsy(:,:,:,:,:,:)

!***********************************************************************

 ABI_MALLOC(pertsy,(3,mpert,3,mpert,3,mpert))
 pertsy(:,:,:,:,:,:) = 0

!Loop over perturbations

 do i1pert_ = 1, mpert
   do i2pert_ = 1, mpert
     do i3pert_ = 1, mpert

       do i1dir_ = 1, 3
         do i2dir_ = 1, 3
           do i3dir_ = 1, 3

             i1pert = (mpert - i1pert_ + 1)
             if (i1pert <= natom) i1pert = natom + 1 - i1pert
             i2pert = (mpert - i2pert_ + 1)
             if (i2pert <= natom) i2pert = natom + 1 - i2pert
             i3pert = (mpert - i3pert_ + 1)
             if (i3pert <= natom) i3pert = natom + 1 - i3pert

             if (i1pert <= natom) then
               i1dir = i1dir_ ; i2dir = i2dir_ ; i3dir = i3dir_
             else if (i2pert <= natom) then
               i1dir = i2dir_ ; i2dir = i1dir_ ; i3dir = i3dir_
             else if (i3pert <= natom) then
               i1dir = i3dir_ ; i2dir = i2dir_ ; i3dir = i1dir_
             else
               i1dir = i1dir_ ; i2dir = i2dir_ ; i3dir = i3dir_
             end if

             if (rfpert(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) /= 0) then

!              Loop over all symmetries

               flag = 0
               do isym = 1, nsym

                 found = 1

!                Select the symmetric element of i1pert,i2pert,i3pert

                 if (i1pert <= natom) then
                   ipesy1 = indsym(4,isym,i1pert)
                   sym1(:,:) = symrec(:,:,isym)
                 else if (i1pert == natom + 2) then
                   ipesy1 = i1pert
                   sym1(:,:) = symrel(:,:,isym)
                 else
                   found = 0
                 end if

                 if (i2pert <= natom) then
                   ipesy2 = indsym(4,isym,i2pert)
                   sym2(:,:) = symrec(:,:,isym)
                 else if (i2pert == natom + 2) then
                   ipesy2 = i2pert
                   sym2(:,:) = symrel(:,:,isym)
                 else
                   found = 0
                 end if

                 if (i3pert <= natom) then
                   ipesy3 = indsym(4,isym,i3pert)
                   sym3(:,:) = symrec(:,:,isym)
                 else if (i3pert == natom + 2) then
                   ipesy3 = i3pert
                   sym3(:,:) = symrel(:,:,isym)
                 else
                   found = 0
                 end if

!                See if the symmetric element is available and check if some
!                of the elements may be zeor. In the latter case, they do not need
!                to be computed.


                 if ((flag /= -1).and.&
&                 (ipesy1==i1pert).and.(ipesy2==i2pert).and.(ipesy3==i3pert)) then
                   flag = sym1(i1dir,i1dir)*sym2(i2dir,i2dir)*sym3(i3dir,i3dir)
                 end if


                 do idisy1 = 1, 3
                   do idisy2 = 1, 3
                     do idisy3 = 1, 3

                       if ((sym1(i1dir,idisy1) /= 0).and.(sym2(i2dir,idisy2) /= 0).and.&
&                       (sym3(i3dir,idisy3) /= 0)) then
                         if (pertsy(idisy1,ipesy1,idisy2,ipesy2,idisy3,ipesy3) == 0) then
                           found = 0
!                          exit      ! exit loop over symmetries
                         end if
                       end if


                       if ((flag == -1).and.&
&                       ((idisy1/=i1dir).or.(idisy2/=i2dir).or.(idisy3/=i3dir))) then
                         if ((sym1(i1dir,idisy1)/=0).and.(sym2(i2dir,idisy2)/=0).and.&
&                         (sym3(i3dir,idisy3)/=0)) then
                           flag = 0
                         end if
                       end if

                     end do
                   end do
                 end do

                 if (found == 1) then
                   pertsy(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = -1
                 end if

!                In case a symmetry operation only changes the sign of an
!                element, this element has to be equal to zero

                 if (flag == -1) then
                   pertsy(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = -2
                   exit
                 end if

               end do    ! close loop on symmetries

!              If the elemetn i1pert,i2pert,i3pert is not symmetric
!              to a basis element, it is a basis element

               if (pertsy(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) > -1) then
                 pertsy(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = 1
               end if

             end if ! rfpert /= 0

           end do        ! close loop over perturbations
         end do
       end do
     end do
   end do
 end do

!Now, take into account the permutation of (i1pert,i1dir)
!and (i3pert,i3dir)

 do i1pert = 1, mpert
   do i2pert = 1, mpert
     do i3pert = 1, mpert

       do i1dir = 1, 3
         do i2dir = 1, 3
           do i3dir = 1, 3

             if ((i1pert /= i3pert).or.(i1dir /= i3dir)) then

               if ((pertsy(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) == 1).and.&
&               (pertsy(i3dir,i3pert,i2dir,i2pert,i1dir,i1pert) == 1)) then
                 pertsy(i3dir,i3pert,i2dir,i2pert,i1dir,i1pert) = -1
               end if

             end if

           end do
         end do
       end do

     end do
   end do
 end do

 rfpert(:,:,:,:,:,:) = pertsy(:,:,:,:,:,:)

 ABI_FREE(pertsy)

end subroutine sytens
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/axial9
!!
!! NAME
!! axial9
!!
!! FUNCTION
!! Generates the local coordinates system from the
!! knowledge of the first vector (longitudinal) and
!! the ifc matrix in cartesian coordinates
!!
!! INPUTS
!! ifccar(3,3)= matrix of interatomic force constants in cartesian coordinates
!! vect1(3)= cartesian coordinates of the first local vector
!!
!! OUTPUT
!! vect2(3)= cartesian coordinates of the second local vector
!! vect3(3)= cartesian coordinates of the third local vector
!!
!! PARENTS
!!      m_ifc
!!
!! CHILDREN
!!
!! SOURCE

subroutine axial9(ifccar,vect1,vect2,vect3)

!Arguments -------------------------------
!arrays
 real(dp),intent(in) :: ifccar(3,3),vect1(3)
 real(dp),intent(out) :: vect2(3),vect3(3)

!Local variables -------------------------
!scalars
 integer :: flag,ii,itrial,jj
 real(dp) :: innorm,scprod
!arrays
 real(dp) :: work(3)

! *********************************************************************

 do jj=1,3
   work(jj)=zero
   do ii=1,3
     work(jj)=work(jj)+ifccar(jj,ii)*vect1(ii)
   end do
 end do

 flag=0
 do itrial=1,4
   scprod=zero
   do ii=1,3
     scprod=scprod+work(ii)*vect1(ii)
   end do

   do ii=1,3
     work(ii)=work(ii)-vect1(ii)*scprod
   end do

   scprod=zero
   do ii=1,3
     scprod=scprod+work(ii)**2
   end do

   if(scprod<1.0d-10)then
     work(1:3)=zero
     if(itrial>1)work(itrial-1)=1.0_dp
   else
     flag=1
   end if

   if(flag==1)exit
 end do

 innorm=scprod**(-0.5_dp)
 do ii=1,3
   vect2(ii)=work(ii)*innorm
 end do

 vect3(1)=vect1(2)*vect2(3)-vect1(3)*vect2(2)
 vect3(2)=vect1(3)*vect2(1)-vect1(1)*vect2(3)
 vect3(3)=vect1(1)*vect2(2)-vect1(2)*vect2(1)

end subroutine axial9
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/dymfz9
!!
!! NAME
!! dymfz9
!!
!! FUNCTION
!! As the subroutine canatm has transformed the coordinates of the
!! atoms in normalized canonical coordinates, the corresponding
!! dynamical matrix should be multiplied by a phase shift corresponding
!! to the translation between New and Old coordinates of its two
!! corresponding atoms.
!!
!! INPUTS
!! dynmat = non-phase shifted dynamical matrices
!! natom = number of atoms
!! nqpt = number of qpoints
!! gprim = reciprocal lattice vectors (cartesian but dimensionless)
!! option=1 : the matrices are transformed from the old (tn)
!!  coordinate system to the new (normalized canonical)
!!        2 : the matrices are restored from the normalized
!!  canonical coordinate system to the usual (tn) one...
!! rcan = canonical coordinates of atoms
!! spqpt = qpoint coordinates (reduced reciprocal)
!! trans = Atomic translations : xred = rcan + trans
!!
!! OUTPUT
!! dynmat = phase shifted dynamical matrices
!!
!! PARENTS
!!      m_dynmat,m_ifc
!!
!! CHILDREN
!!
!! SOURCE

subroutine dymfz9(dynmat,natom,nqpt,gprim,option,spqpt,trans)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom,nqpt,option
!arrays
 real(dp),intent(in) :: gprim(3,3),spqpt(3,nqpt),trans(3,natom)
 real(dp),intent(inout) :: dynmat(2,3,natom,3,natom,nqpt)

!Local variables -------------------------
!scalars
 integer :: ia,ib,iqpt,mu,nu
 real(dp) :: im,ktrans,re
!arrays
 real(dp) :: kk(3)

! *********************************************************************

 do iqpt=1,nqpt
   ! Definition of q in normalized reciprocal space
   kk(1)=spqpt(1,iqpt)*gprim(1,1)+spqpt(2,iqpt)*gprim(1,2)+spqpt(3,iqpt)*gprim(1,3)
   kk(2)=spqpt(1,iqpt)*gprim(2,1)+spqpt(2,iqpt)*gprim(2,2)+spqpt(3,iqpt)*gprim(2,3)
   kk(3)=spqpt(1,iqpt)*gprim(3,1)+spqpt(2,iqpt)*gprim(3,2)+spqpt(3,iqpt)*gprim(3,3)

   if(option==1)then
     kk(1)=-kk(1)
     kk(2)=-kk(2)
     kk(3)=-kk(3)
   end if

   do ia=1,natom
     do ib=1,natom
       ! Product of q with the differences between the two atomic translations
       ktrans=kk(1)*(trans(1,ia)-trans(1,ib))+kk(2)*(trans(2,ia)-trans(2,ib))+kk(3)*(trans(3,ia)-trans(3,ib))
       do mu=1,3
         do nu=1,3
           re=dynmat(1,mu,ia,nu,ib,iqpt)
           im=dynmat(2,mu,ia,nu,ib,iqpt)
           ! Transformation of the Old dynamical matrices by New ones by multiplication by a phase shift
           dynmat(1,mu,ia,nu,ib,iqpt)=re*cos(two_pi*ktrans)-im*sin(two_pi*ktrans)
           dynmat(2,mu,ia,nu,ib,iqpt)=re*sin(two_pi*ktrans)+im*cos(two_pi*ktrans)
         end do
       end do
     end do
   end do
 end do

end subroutine dymfz9
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/nanal9
!!
!! NAME
!! nanal9
!!
!! FUNCTION
!! If plus=0 then substracts the non-analytical part from one dynamical
!!           matrices, with number iqpt.
!! If plus=1 then adds the non-analytical part to the dynamical
!!           matrices, with number iqpt.
!!
!! INPUTS
!! dyew(2,3,natom,3,natom)= Non-analytical part
!! natom= Number of atoms in the unit cell
!! iqpt= Referenced q point for the dynamical matrix
!! nqpt= Number of q points
!! plus= (see above)
!!
!! OUTPUT
!! dynmat(2,3,natom,3,natom,nqpt)= Dynamical matrices coming from the Derivative Data Base
!!
!! PARENTS
!!      m_dynmat,m_ifc
!!
!! CHILDREN
!!
!! SOURCE

subroutine nanal9(dyew,dynmat,iqpt,natom,nqpt,plus)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iqpt,natom,nqpt,plus
!arrays
 real(dp),intent(in) :: dyew(2,3,natom,3,natom)
 real(dp),intent(inout) :: dynmat(2,3,natom,3,natom,nqpt)

!Local variables -------------------------
!scalars
 integer :: ia,ib,mu,nu
 character(len=500) :: msg

! *********************************************************************

 if (plus==0) then

   do ia=1,natom
     do ib=1,natom
       do mu=1,3
         do nu=1,3
           ! The following four lines are OK
           dynmat(1,mu,ia,nu,ib,iqpt)=dynmat(1,mu,ia,nu,ib,iqpt) - dyew(1,mu,ia,nu,ib)
           dynmat(2,mu,ia,nu,ib,iqpt)=dynmat(2,mu,ia,nu,ib,iqpt) - dyew(2,mu,ia,nu,ib)
         end do
       end do
     end do
   end do

 else if (plus==1) then
   do ia=1,natom
     do ib=1,natom
       do mu=1,3
         do nu=1,3
           dynmat(1,mu,ia,nu,ib,iqpt)=dynmat(1,mu,ia,nu,ib,iqpt) + dyew(1,mu,ia,nu,ib)
           dynmat(2,mu,ia,nu,ib,iqpt)=dynmat(2,mu,ia,nu,ib,iqpt) + dyew(2,mu,ia,nu,ib)
         end do
       end do
     end do
   end do

 else
   write(msg,'(3a,i0,a)' )&
    'The argument "plus" must be equal to 0 or 1.',ch10,&
    'The value ',plus,' is not available.'
   MSG_BUG(msg)
 end if

end subroutine nanal9
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/gtdyn9
!!
!! NAME
!! gtdyn9
!!
!! FUNCTION
!! Generates a dynamical matrix from interatomic force
!! constants and long-range electrostatic interactions.
!!
!! INPUTS
!! acell(3)=length scales by which rprim is to be multiplied
!! atmfrc(3,natom,3,natom,nrpt) = Interatomic Forces in real space
!!  (imaginary part only for debugging)
!! dielt(3,3) = dielectric tensor
!! dipdip= if 0, no dipole-dipole interaction was subtracted in atmfrc
!!  if 1, atmfrc has been build without dipole-dipole part
!! dyewq0(3,3,natom)= Ewald part of the dynamical matrix, at q=0.
!! gmet(3,3)= metric tensor in reciprocal space.
!! gprim(3,3)= Normalized coordinates in reciprocal space
!! mpert =maximum number of ipert
!! natom= Number of atoms in the unit cell
!! nrpt= Number of R points in the Big Box
!! qphnrm= Normalisation coefficient for qpt
!! qpt(3)= Reduced coordinates of the q vectors in reciprocal space
!! rmet(3,3)= Metric tensor in real space.
!! rprim(3,3)= dimensionless primitive translations in real space
!! rpt(3,nprt)= Canonical coordinates of the R points in the unit cell
!!  These coordinates are normalized (=> * acell(3)!!)
!! trans(3,natom)= Atomic translations : xred = rcan + trans
!! ucvol= unit cell volume
!! wghatm(natom,natom,nrpt)= Weights associated to a pair of atoms and to a R vector
!! xred(3,natom)= relative coords of atoms in unit cell (dimensionless)
!! zeff(3,3,natom)=effective charge on each atom, versus electric field and atomic displacement
!! comm=MPI communicator.
!! [dipquad] = if 1, atmfrc has been build without dipole-quadrupole part
!! [quadquad] = if 1, atmfrc has been build without quadrupole-quadrupole part
!!
!! OUTPUT
!! d2cart(2,3,mpert,3,mpert)=dynamical matrix obtained for the wavevector qpt (normalized using qphnrm)
!!
!! PARENTS
!!      anaddb,ddb_interpolate,m_effective_potential_file,m_ifc,m_phonons
!!
!! CHILDREN
!!
!! SOURCE

subroutine gtdyn9(acell,atmfrc,dielt,dipdip,dyewq0,d2cart,gmet,gprim,mpert,natom,&
                  nrpt,qphnrm,qpt,rmet,rprim,rpt,trans,ucvol,wghatm,xred,zeff,qdrp_cart,ewald_option,comm,&
                  dipquad,quadquad)  ! optional

!Arguments -------------------------------
!scalars
 integer,intent(in) :: dipdip,mpert,natom,nrpt,ewald_option,comm
 real(dp),intent(in) :: qphnrm,ucvol
 integer,optional,intent(in) :: dipquad, quadquad
!arrays
 real(dp),intent(in) :: acell(3),dielt(3,3),gmet(3,3),gprim(3,3),qpt(3)
 real(dp),intent(in) :: rmet(3,3),rprim(3,3),rpt(3,nrpt)
 real(dp),intent(in) :: trans(3,natom),wghatm(natom,natom,nrpt),xred(3,natom)
 real(dp),intent(in) :: zeff(3,3,natom)
 real(dp),intent(in) :: qdrp_cart(3,3,3,natom)
 real(dp),intent(in) :: atmfrc(3,natom,3,natom,nrpt)
 real(dp),intent(in) :: dyewq0(3,3,natom)
 real(dp),intent(out) :: d2cart(2,3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer,parameter :: nqpt1 = 1, option2 = 2, sumg0 = 0, plus1 = 1, iqpt1 = 1
 integer :: i1, i2, ib, nsize, dipquad_, quadquad_
!arrays
 real(dp) :: qphon(3) !, tsec(2)
 real(dp),allocatable :: dq(:,:,:,:,:),dyew(:,:,:,:,:)

! *********************************************************************

 ! Keep track of time spent in gtdyn9
 !call timab(1750, 1, tsec)

 ABI_MALLOC(dq,(2,3,natom,3,natom))

 ! Define quadrupolar options
 dipquad_=0; if(present(dipquad)) dipquad_=dipquad
 quadquad_=0; if(present(quadquad)) quadquad_=quadquad

 ! Get the normalized wavevector
 if(abs(qphnrm)<1.0d-7)then
   qphon(1:3)=zero
 else
   qphon(1:3)=qpt(1:3)/qphnrm
 end if

 ! Generate the analytical part from the interatomic forces
 call ftifc_r2q(atmfrc, dq, gprim, natom, nqpt1, nrpt, rpt, qphon, wghatm, comm)

 ! The analytical dynamical matrix dq has been generated
 ! in the normalized canonical coordinate system.
 ! Now, the phase is modified, in order to recover the usual (xred) coordinate of atoms.
 call dymfz9(dq,natom,nqpt1,gprim,option2,qphon,trans)

 if (dipdip==1.or.dipquad_==1.or.quadquad_==1) then
   ! Add the non-analytical part
   ! Compute dyew(2,3,natom,3,natom)= Ewald part of the dynamical matrix,
   ! second energy derivative wrt xred(3,natom) in Hartrees (Denoted A-bar in the notes)
   ABI_MALLOC(dyew,(2,3,natom,3,natom))

   call ewald9(acell,dielt,dyew,gmet,gprim,natom,qphon,rmet,rprim,sumg0,ucvol,xred,zeff,&
      qdrp_cart,option=ewald_option,dipquad=dipquad_,quadquad=quadquad_)

   call q0dy3_apply(natom,dyewq0,dyew)
   call nanal9(dyew,dq,iqpt1,natom,nqpt1,plus1)

   ABI_FREE(dyew)
 end if

 ! Copy the dynamical matrix in the proper location
 ! First zero all the elements
 nsize=2*(3*mpert)**2
 d2cart = zero

 ! Copy the elements from dq to d2cart
 d2cart(:,:,1:natom,:,1:natom)=dq(:,:,1:natom,:,1:natom)

 ! In case we have the gamma point,
 if(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-14)then
   ! Copy the effective charge and dielectric constant in the final array
   do i1=1,3
     do i2=1,3
       d2cart(1,i1,natom+2,i2,natom+2)=dielt(i1,i2)
       do ib=1,natom
         d2cart(1,i1,natom+2,i2,ib)=zeff(i1,i2,ib)
         d2cart(1,i2,ib,i1,natom+2)=zeff(i1,i2,ib)
       end do
     end do
   end do
 end if

 ABI_FREE(dq)

 !call timab(1750, 2, tsec)

end subroutine gtdyn9
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/dfpt_phfrq
!!
!! NAME
!! dfpt_phfrq
!!
!! FUNCTION
!! Get the phonon frequencies and eigenvectors (as well as the corresponding displacements)
!! If q is at Gamma, the non-analytical behaviour can be included.
!! Then, the effective dielectric tensor, the effective charges
!! and oscillator strengths for the limiting direction are also returned
!!
!! INPUTS
!!  amu(ntypat)=mass of the atoms (atomic mass unit) matrix (diagonal in the atoms)
!!  d2cart(2,3,mpert,3,mpert)=dynamical matrix, effective charges, dielectric tensor,.... all in cartesian coordinates
!!  indsym(4,msym*natom)=indirect indexing array : for each
!!   isym,iatom, fourth element is label of atom into which iatom is sent by
!!   INVERSE of symmetry operation isym; first three elements are the primitive
!!   translations which must be subtracted after the transformation to get back to the original unit cell.
!!  mpert =maximum number of ipert
!!  msym=maximum number of symmetries
!!  natom=number of atoms in unit cell
!!  nsym=number of space group symmetries
!!  ntypat=number of atom types
!!  qphnrm=(described above)
!!  qphon(3)= to be divided by qphnrm, give the phonon wavevector;
!!     if qphnrm==0.0_dp, then the wavevector is zero (Gamma point)
!!     and qphon gives the direction of the induced electric field in **CARTESIAN** coordinates.
!!     in the latter case, if qphon is zero, no non-analytical contribution is included.
!!  rprimd(3,3)=dimensional primitive translations (bohr)
!!  symdynmat=if 1, (re)symmetrize the dynamical matrix, except if Gamma wavevector with electric field added.
!!  symrel(3,3,nsym)=matrices of the group symmetries (real space)
!!  typat(natom)=integer label of each type of atom (1,2,...)
!!  ucvol=unit cell volume
!!
!! OUTPUT
!!  displ(2*3*natom*3*natom)= at the end, contains the displacements of atoms in cartesian coordinates.
!!    The first index means either the real or the imaginary part,
!!    The second index runs on the direction and the atoms displaced
!!    The third index runs on the modes.
!!  eigval(3*natom)=contains the eigenvalues of the dynamical matrix
!!  eigvec(2*3*natom*3*natom)= at the end, contains the eigenvectors of the dynamical matrix in cartesian coordinates.
!!  phfrq(3*natom)=phonon frequencies (square root of the dynamical matrix eigenvalues,
!!    except if these are negative, and in this case, give minus the square root of the absolute value
!!    of the matrix eigenvalues). Hartree units.
!!
!! NOTES
!!   1) One makes the dynamical matrix hermitian...
!!   2) In case of q=Gamma, only the real part is used.
!!
!! PARENTS
!!      anaddb,m_ddb,m_effective_potential_file,m_ifc,m_phonons,respfn,thmeig
!!
!! CHILDREN
!!
!! SOURCE

subroutine dfpt_phfrq(amu,displ,d2cart,eigval,eigvec,indsym,&
& mpert,msym,natom,nsym,ntypat,phfrq,qphnrm,qphon,rprimd,&
& symdynmat,symrel,symafm,typat,ucvol)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,msym,natom,nsym,ntypat,symdynmat
 real(dp),intent(in) :: qphnrm,ucvol
!arrays
 integer,intent(in) :: indsym(4,msym*natom),symrel(3,3,nsym),typat(natom)
 integer,intent(in) :: symafm(nsym)
 real(dp),intent(in) :: amu(ntypat),d2cart(2,3,mpert,3,mpert),rprimd(3,3)
 real(dp),intent(inout) :: qphon(3)
 real(dp),intent(out) :: displ(2*3*natom*3*natom),eigval(3*natom)
 real(dp),intent(out) :: eigvec(2*3*natom*3*natom),phfrq(3*natom)

!Local variables -------------------------
!scalars
 integer :: analyt,i1,i2,idir1,idir2,ier,ii,imode,ipert1,ipert2
 integer :: jmode,indexi,indexj,index
 real(dp) :: epsq,qphon2
 logical,parameter :: debug = .False.
 real(dp) :: sc_prod
!arrays
 real(dp) :: qptn(3),dum(2,0) !, tsec(2)
 real(dp),allocatable :: matrx(:,:),zeff(:,:),zhpev1(:,:),zhpev2(:)

! *********************************************************************

 ! Keep track of time spent in dfpt_phfrq
 !call timab(1751, 1, tsec)

 ! Prepare the diagonalisation: analytical part.
 ! Note: displ is used as work space here
 i1=0
 do ipert1=1,natom
   do idir1=1,3
     i1=i1+1
     i2=0
     do ipert2=1,natom
       do idir2=1,3
         i2=i2+1
         index=i1+3*natom*(i2-1)
         displ(2*index-1)=d2cart(1,idir1,ipert1,idir2,ipert2)
         displ(2*index  )=d2cart(2,idir1,ipert1,idir2,ipert2)
       end do
     end do
   end do
 end do

 ! Determine the analyticity of the matrix.
 analyt=1; if(abs(qphnrm)<tol8) analyt=0
 if(abs(qphon(1))<tol8.and.abs(qphon(2))<tol8.and.abs(qphon(3))<tol8) analyt=2

 ! In case of q=Gamma, only the real part is used
 if(analyt==0 .or. analyt==2)then
   do i1=1,3*natom
     do i2=1,3*natom
       index=i1+3*natom*(i2-1)
       displ(2*index)=zero
     end do
   end do
 end if

 ! In the case the non-analyticity is required:
 ! the tensor is in cartesian coordinates and this means that qphon must be in given in Cartesian coordinates.
 if(analyt==0)then

   ! Normalize the limiting direction
   qphon2=qphon(1)**2+qphon(2)**2+qphon(3)**2
   qphon(:)=qphon(:)/sqrt(qphon2)

   ! Get the dielectric constant for the limiting direction
   epsq=zero
   do idir1=1,3
     do idir2=1,3
       epsq= epsq + qphon(idir1)*qphon(idir2) * d2cart(1,idir1,natom+2,idir2,natom+2)
     end do
   end do

   ABI_MALLOC(zeff,(3,natom))

   ! Get the effective charges for the limiting direction
   do idir1=1,3
     do ipert1=1,natom
       zeff(idir1,ipert1)=zero
       do idir2=1,3
         zeff(idir1,ipert1) = zeff(idir1,ipert1) + qphon(idir2)* d2cart(1,idir1,ipert1,idir2,natom+2)
       end do
     end do
   end do

   ! Get the non-analytical part of the dynamical matrix, and suppress its imaginary part.
   i1=0
   do ipert1=1,natom
     do idir1=1,3
       i1=i1+1
       i2=0
       do ipert2=1,natom
         do idir2=1,3
           i2=i2+1
           index=i1+3*natom*(i2-1)
           displ(2*index-1)=displ(2*index-1)+four_pi/ucvol*zeff(idir1,ipert1)*zeff(idir2,ipert2)/epsq
           displ(2*index  )=zero
         end do
       end do
     end do
   end do

   ABI_FREE(zeff)
 end if !  End of the non-analyticity treatment

 ! Multiply IFC(q) by masses
 call massmult_and_breaksym(natom, ntypat, typat, amu, displ)

 ! ***********************************************************************
 ! Diagonalize the dynamical matrix

 !Symmetrize the dynamical matrix
 !FIXME: swap the next 2 lines and update test files to include symmetrization
 !       for Gamma point too (except in non-analytic case)
 !if (symdynmat==1 .and. analyt > 0) then
 if (symdynmat==1 .and. analyt == 1) then
   qptn(:)=qphon(:)
   if (analyt==1) qptn(:)=qphon(:)/qphnrm
   call symdyma(displ,indsym,natom,nsym,qptn,rprimd,symrel,symafm)
 end if

 ii=1
 ABI_MALLOC(matrx,(2,(3*natom*(3*natom+1))/2))
 do i2=1,3*natom
   do i1=1,i2
     matrx(1,ii)=displ(1+2*(i1-1)+2*(i2-1)*3*natom)
     matrx(2,ii)=displ(2+2*(i1-1)+2*(i2-1)*3*natom)
     ii=ii+1
   end do
 end do

 ABI_MALLOC(zhpev1,(2,2*3*natom-1))
 ABI_MALLOC(zhpev2,(3*3*natom-2))

 call ZHPEV ('V','U',3*natom,matrx,eigval,eigvec,3*natom,zhpev1,zhpev2,ier)
 ABI_CHECK(ier == 0, sjoin('zhpev returned:', itoa(ier)))

 ABI_FREE(matrx)
 ABI_FREE(zhpev1)
 ABI_FREE(zhpev2)

 if (debug) then
   ! Check the orthonormality of the eigenvectors
   do imode=1,3*natom
     do jmode=imode,3*natom
       indexi=2*3*natom*(imode-1)
       indexj=2*3*natom*(jmode-1)
       sc_prod=sum(eigvec(indexi+1:indexi+6*natom)*eigvec(indexj+1:indexj+6*natom))
       write(std_out,'(a,2i4,a,es16.6)')' imode,jmode=',imode,jmode,' real scalar product =',sc_prod
     end do
   end do
 end if

 !***********************************************************************

 ! Get the phonon frequencies (negative by convention, if the eigenvalue of the dynamical matrix is negative)
 do imode=1,3*natom
   if(eigval(imode)>=1.0d-16)then
     phfrq(imode)=sqrt(eigval(imode))
   else if(eigval(imode)>=-1.0d-16)then
     phfrq(imode)=zero
   else
     phfrq(imode)=-sqrt(-eigval(imode))
   end if
 end do

 ! Fix the phase of the eigenvectors
 call fxphas_seq(eigvec,dum, 0, 0, 1, 3*natom*3*natom, 0, 3*natom, 3*natom, 0)

 ! Normalise the eigenvectors
 call pheigvec_normalize(natom, eigvec)

 ! Get the phonon displacements
 call phdispl_from_eigvec(natom, ntypat, typat, amu, eigvec, displ)

 if (debug) then
   write(std_out,'(a)')' Phonon eigenvectors and displacements '
   do imode=1,3*natom
     indexi=2*3*natom*(imode-1)
     write(std_out,'(a,i4,a,12es16.6)')' imode=',imode,' eigvec(1:6*natom)=',eigvec(indexi+1:indexi+6*natom)
     write(std_out,'(a,i4,a,12es16.6)')' imode=',imode,' displ(1:6*natom)=',displ(indexi+1:indexi+6*natom)
   end do

   ! Check the orthonormality of the eigenvectors
   do imode=1,3*natom
     do jmode=imode,3*natom
       indexi=2*3*natom*(imode-1)
       indexj=2*3*natom*(jmode-1)
       sc_prod=sum(eigvec(indexi+1:indexi+6*natom)*eigvec(indexj+1:indexj+6*natom))
       write(std_out,'(a,2i4,a,es16.6)')' imode,jmode=',imode,jmode,' real scalar product =',sc_prod
     end do
   end do
 end if

 !call timab(1751, 2, tsec)

end subroutine dfpt_phfrq
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/pheigvec_normalize
!!
!! NAME
!! pheigvec_normalize
!!
!! FUNCTION
!!  Normalize input eigenvectors in cartesian coordinates
!!
!! INPUTS
!!  natom: number of atoms in unit cell
!!
!! SIDE EFFECTS
!!  eigvec(2*3*natom*3*natom)=in output the normalized eigenvectors in cartesian coordinates.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure subroutine pheigvec_normalize(natom, eigvec)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom
!arrays
 real(dp),intent(inout) :: eigvec(2*3*natom*3*natom)

!Local variables -------------------------
!scalars
 integer :: i1,idir1,imode,ipert1,index
 real(dp) :: norm

! *********************************************************************

 do imode=1,3*natom

   norm=zero
   do idir1=1,3
     do ipert1=1,natom
       i1=idir1+(ipert1-1)*3
       index=i1+3*natom*(imode-1)
       norm=norm+eigvec(2*index-1)**2+eigvec(2*index)**2
     end do
   end do
   norm=sqrt(norm)

   do idir1=1,3
     do ipert1=1,natom
       i1=idir1+(ipert1-1)*3
       index=i1+3*natom*(imode-1)
       eigvec(2*index-1)=eigvec(2*index-1)/norm
       eigvec(2*index)=eigvec(2*index)/norm
     end do
   end do

 end do

end subroutine pheigvec_normalize
!!***

!----------------------------------------------------------------------

!!****f* m_dynmat/phdispl_from_eigvec
!!
!! NAME
!! phdispl_from_eigvec
!!
!! FUNCTION
!!  Phonon displacements from eigenvectors
!!
!! INPUTS
!!  natom: number of atoms in unit cell
!!  ntypat=number of atom types
!!  typat(natom)=integer label of each type of atom (1,2,...)
!!  amu(ntypat)=mass of the atoms (atomic mass unit) matrix (diagonal in the atoms)
!!  eigvec(2*3*natom*3*natom)= eigenvectors of the dynamical matrix in cartesian coordinates.
!!
!! OUTPUT
!!  displ(2*3*natom*3*natom)=displacements of atoms in cartesian coordinates.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure subroutine phdispl_from_eigvec(natom, ntypat, typat, amu, eigvec, displ)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom, ntypat
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: amu(ntypat)
 real(dp),intent(in) :: eigvec(2*3*natom*3*natom)
 real(dp),intent(out) :: displ(2*3*natom*3*natom)

!Local variables -------------------------
!scalars
 integer :: i1,idir1,imode,ipert1, index

! *********************************************************************

 do imode=1,3*natom

   do idir1=1,3
     do ipert1=1,natom
       i1=idir1+(ipert1-1)*3
       index=i1+3*natom*(imode-1)
       displ(2*index-1)=eigvec(2*index-1) / sqrt(amu(typat(ipert1))*amu_emass)
       displ(2*index  )=eigvec(2*index  ) / sqrt(amu(typat(ipert1))*amu_emass)
     end do
   end do

 end do

end subroutine phdispl_from_eigvec
!!!***

!!****f* m_dynmat/dfpt_prtph
!! NAME
!! dfpt_prtph
!!
!! FUNCTION
!! Print the phonon frequencies, on unit 6 as well as the printing
!! unit (except if the associated number -iout- is negative),
!! and for the latter, in Hartree, meV, Thz, Kelvin or cm-1.
!! If eivec==1,2, also print the eigenmodes : displacements in cartesian coordinates.
!! If eivec==4, generate output files for band2eps (drawing tool for the phonon band structure
!!
!! INPUTS
!!  displ(2,3*natom,3*natom)= contains the displacements of atoms in cartesian coordinates.
!!  The first index means either the real or the imaginary part,
!!  The second index runs on the direction and the atoms displaced
!!  The third index runs on the modes.
!!  eivec=(if eivec==0, the eigendisplacements are not printed,
!!    if eivec==1,2, the eigendisplacements are printed,
!!    if eivec==4, files for band2eps
!!  enunit=units for output of the phonon frequencies :
!!    0=> Hartree and cm-1, 1=> eV and Thz, other=> Ha,Thz,eV,cm-1 and K
!!  iout= unit for long print (if negative, the routine only print on unit 6, and in Hartree only).
!!  natom= number of atom
!!  phfreq(3*natom)= phonon frequencies in Hartree
!!  qphnrm=phonon wavevector normalisation factor
!!  qphon(3)=phonon wavevector
!!
!! OUTPUT
!!  Only printing
!!
!! NOTES
!! called by one processor only
!!
!! PARENTS
!!      anaddb,m_effective_potential_file,m_ifc,m_phonons,respfn
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine dfpt_prtph(displ,eivec,enunit,iout,natom,phfrq,qphnrm,qphon)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: eivec,enunit,iout,natom
 real(dp),intent(in) :: qphnrm
!arrays
 real(dp),intent(in) :: displ(2,3*natom,3*natom),phfrq(3*natom),qphon(3)

!Local variables -------------------------
!scalars
 integer :: i,idir,ii,imode,jj
 real(dp) :: tolerance
 logical :: t_degenerate
 character(len=500) :: msg
!arrays
 real(dp) :: vecti(3),vectr(3)
 character(len=1) :: metacharacter(3*natom)

! *********************************************************************

!Check the value of eivec
 if (all(eivec /= [0,1,2,4])) then
   write(msg, '(a,i0,a,a)' )&
   'In the calling subroutine, eivec is',eivec,ch10,&
   'but allowed values are between 0 and 4.'
   MSG_BUG(msg)
 end if

!write the phonon frequencies on unit std_out
 write(msg,'(4a)' )' ',ch10,' phonon wavelength (reduced coordinates) , ','norm, and energies in hartree'
 call wrtout(std_out,msg)

!The next format should be rewritten
 write(msg,'(a,4f5.2)' )' ',(qphon(i),i=1,3),qphnrm
 call wrtout(std_out,msg)
 do jj=1,3*natom,5
   if (3*natom-jj<5) then
     write(msg,'(5es17.9)') (phfrq(ii),ii=jj,3*natom)
   else
     write(msg,'(5es17.9)') (phfrq(ii),ii=jj,jj+4)
   end if
   call wrtout(std_out,msg)
 end do
 write(msg,'(a,a,es17.9)') ch10,' Zero Point Motion energy (sum of freqs/2)=',sum(phfrq(1:3*natom))/2
 call wrtout(std_out,msg)

!Put the wavevector in nice format
 if(iout>=0)then
   call wrtout(iout,' ')
   if(qphnrm/=0.0_dp)then
     write(msg, '(a,3f9.5)' )&
     '  Phonon wavevector (reduced coordinates) :',(qphon(i)/qphnrm+tol10,i=1,3)
   else
     write(msg, '(3a,3f9.5)' )&
     '  Phonon at Gamma, with non-analyticity in the',ch10,&
     '  direction (cartesian coordinates)',qphon(1:3)+tol10
   end if
   call wrtout(iout,msg)

!  Write it, in different units.
   if(enunit/=1)then
     write(iout, '(a)' )' Phonon energies in Hartree :'
     do jj=1,3*natom,5
       if (3*natom-jj<5) then
         write(msg, '(1x,5es14.6)') (phfrq(ii),ii=jj,3*natom)
       else
         write(msg, '(1x,5es14.6)') (phfrq(ii),ii=jj,jj+4)
       end if
       call wrtout(iout,msg)
     end do
   end if
   if(enunit/=0)then
     write(iout, '(a)' )' Phonon energies in meV     :'
     do jj=1,3*natom,5
       if (3*natom-jj<5) then
         write(msg, '("-",5es14.6)') (phfrq(ii)*Ha_eV*1.0d3,ii=jj,3*natom)
       else
         write(msg, '("-",5es14.6)') (phfrq(ii)*Ha_eV*1.0d3,ii=jj,jj+4)
       end if
       call wrtout(iout,msg)
     end do
   end if
   if(enunit/=1)then
     write(iout, '(a)' )' Phonon frequencies in cm-1    :'
     do jj=1,3*natom,5
       if (3*natom-jj<5) then
         write(msg, '("-",5es14.6)') (phfrq(ii)*Ha_cmm1,ii=jj,3*natom)
       else
         write(msg, '("-",5es14.6)') (phfrq(ii)*Ha_cmm1,ii=jj,jj+4)
       end if
       call wrtout(iout,msg)
     end do
   end if
   if(enunit/=0)then
     write(iout, '(a)' )' Phonon frequencies in Thz     :'
     do jj=1,3*natom,5
       if (3*natom-jj<5) then
         write(msg, '("-",5es14.6)') (phfrq(ii)*Ha_THz,ii=jj,3*natom)
       else
         write(msg, '("-",5es14.6)') (phfrq(ii)*Ha_THz,ii=jj,jj+4)
       end if
       call wrtout(iout,msg)
     end do
   end if
   if(enunit/=0.and.enunit/=1)then
     write(iout, '(a)' )' Phonon energies in Kelvin  :'
     do jj=1,3*natom,5
       if (3*natom-jj<5) then
         write(msg, '("-",5es14.6)') (phfrq(ii)/kb_HaK,ii=jj,3*natom)
       else
         write(msg, '("-",5es14.6)') (phfrq(ii)/kb_HaK,ii=jj,jj+4)
       end if
       call wrtout(iout,msg)
     end do
   end if
 end if

!Take care of the eigendisplacements
 if(eivec==1 .or. eivec==2)then
   write(msg, '(a,a,a,a,a,a,a,a)' ) ch10,&
   ' Eigendisplacements ',ch10,&
   ' (will be given, for each mode : in cartesian coordinates',ch10,&
   '   for each atom the real part of the displacement vector,',ch10,&
   '   then the imaginary part of the displacement vector - absolute values smaller than 1.0d-7 are set to zero)'
   call wrtout(std_out,msg)
   if(iout>=0) then
     call wrtout(iout,msg)
   end if

!  Examine the degeneracy of each mode. The portability of the echo of the eigendisplacements
!  is very hard to obtain, and has not been attempted.
   do imode=1,3*natom
!    The degenerate modes are not portable
     t_degenerate=.false.
     if(imode>1)then
       if(phfrq(imode)-phfrq(imode-1)<tol6)t_degenerate=.true.
     end if
     if(imode<3*natom)then
       if(phfrq(imode+1)-phfrq(imode)<tol6)t_degenerate=.true.
     end if
     metacharacter(imode)=';'; if(t_degenerate)metacharacter(imode)='-'
   end do

   do imode=1,3*natom
     write(msg,'(a,i4,a,es16.6)' )'  Mode number ',imode,'   Energy',phfrq(imode)
     call wrtout(std_out,msg)
     if(iout>=0)then
       write(msg, '(a,i4,a,es16.6)' )'  Mode number ',imode,'   Energy',phfrq(imode)
       call wrtout(iout,msg)
     end if
     tolerance=1.0d-7
     if(abs(phfrq(imode))<1.0d-5)tolerance=2.0d-7
     if(phfrq(imode)<1.0d-5)then
       write(msg,'(3a)' )' Attention : low frequency mode.',ch10,&
       '   (Could be unstable or acoustic mode)'
       call wrtout(std_out,msg)
       if(iout>=0)then
         write(iout, '(3a)' )' Attention : low frequency mode.',ch10,&
         '   (Could be unstable or acoustic mode)'
       end if
     end if
     do ii=1,natom
       do idir=1,3
         vectr(idir)=displ(1,idir+(ii-1)*3,imode)
         if(abs(vectr(idir))<tolerance)vectr(idir)=0.0_dp
         vecti(idir)=displ(2,idir+(ii-1)*3,imode)
         if(abs(vecti(idir))<tolerance)vecti(idir)=0.0_dp
       end do
       write(msg,'(i4,3es16.8,a,4x,3es16.8)' ) ii,vectr(:),ch10,vecti(:)
       call wrtout(std_out,msg)
       if(iout>=0)then
         write(msg,'(a,i3,3es16.8,2a,3x,3es16.8)') metacharacter(imode),ii,vectr(:),ch10,&
           metacharacter(imode), vecti(:)
         call wrtout(iout,msg)
       end if
     end do
   end do
 end if

end subroutine dfpt_prtph
!!***

!!****f* m_dynmat/massmult_and_breaksym
!!
!! NAME
!!  mult_masses_and_break_symms
!!
!! FUNCTION
!!  Multiply the IFC(q) by the atomic masses, slightly break symmetry to make tests more
!!  portable and make the matrix hermitian before returning.
!!
!! INPUTS
!!  amu(ntypat)=mass of the atoms (atomic mass unit) matrix (diagonal in the atoms)
!!  natom=number of atoms in unit cell
!!  ntypat=number of atom types
!!  typat(natom)=integer label of each type of atom (1,2,...)
!!
!! SIDE EFFECTS
!!  mat(2*3*natom*3*natom)=Multiplies by atomic masses in output.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine massmult_and_breaksym(natom, ntypat, typat, amu, mat)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom,ntypat
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: amu(ntypat)
 real(dp),intent(inout) :: mat(2*3*natom*3*natom)

!Local variables -------------------------
!scalars
 integer :: i1,i2,idir1,idir2,index,ipert1,ipert2
 real(dp),parameter :: break_symm=1.0d-12
 !real(dp),parameter :: break_symm=zero
 real(dp) :: fac
!arrays
 real(dp) :: nearidentity(3,3)

! *********************************************************************

 ! This slight breaking of the symmetry allows the results to be more portable between machines
 nearidentity(:,:)=one
 nearidentity(1,1)=one+break_symm
 nearidentity(3,3)=one-break_symm

 ! Include the masses in the dynamical matrix
 do ipert1=1,natom
   do ipert2=1,natom
     fac=1.0_dp/sqrt(amu(typat(ipert1))*amu(typat(ipert2)))/amu_emass
     do idir1=1,3
       do idir2=1,3
         i1=idir1+(ipert1-1)*3
         i2=idir2+(ipert2-1)*3
         index=i1+3*natom*(i2-1)
         mat(2*index-1)=mat(2*index-1)*fac*nearidentity(idir1,idir2)
         mat(2*index  )=mat(2*index  )*fac*nearidentity(idir1,idir2)
         ! This is to break slightly the translation invariance, and make the automatic tests more portable
         if(ipert1==ipert2 .and. idir1==idir2)then
           mat(2*index-1)=mat(2*index-1)+break_symm*natom/amu_emass/idir1*0.01_dp
         end if
       end do
     end do
   end do
 end do

 ! Make the dynamical matrix hermitian
 call mkherm(mat,3*natom)

end subroutine massmult_and_breaksym
!!***

!!****f* m_dynmat/ftgam
!!
!! NAME
!! ftgam
!!
!! FUNCTION
!! If qtor=1 (q->r):
!!  Generates the Fourier transform of the recip space gkk matrices
!!  to obtain the real space ones.
!! If qtor=0 (r->q):
!!  Generates the Fourier transform of the real space gkk matrices
!!  to obtain the reciprocal space ones.
!!
!! INPUTS
!! natom= Number of atoms in the unit cell
!! nqpt= Number of q points in the Brillouin zone
!!           if qtor=0 this number is read in the input file
!! nrpt= Number of R points in the Big Box
!! qtor= ( q to r : see above )
!! rpt(3,nprt)= Canonical coordinates of the R points in the unit cell
!!           These coordinates are normalized (=> * acell(3)!!)
!! qpt_full(3,nqpt)= Reduced coordinates of the q vectors in reciprocal space
!!           if qtor=0 these vectors are read in the input file
!! wghatm(natom,natom,nrpt)= Weights associated to a pair of atoms and to a R vector
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/output
!! gam_qpt(2,3*natom*3*natom,nqpt)
!!  = gamma matrices in recip space coming from the Derivative Data Base
!! gam_rpt(2,3*natom*3*natom,nrpt)
!!  = gamma matrices in real space stored in file unit_gkk_rpt
!!
!! PARENTS
!!      elphon,get_tau_k,integrate_gamma_alt,m_phgamma,mka2f,mka2f_tr
!!      mka2f_tr_lova,mkph_linwid
!!
!! CHILDREN
!!
!! NOTES
!!   copied from ftiaf9.f
!!   recip to real space: real space is forced to disk file unit_gkk_rpt
!!                        recip space depends on gkqwrite and unitgkq3
!!   real to recip space: real space is forced to disk file unit_gkk_rpt
!!                        recip space is necessarily in memory in gkk_qpt
!!
!!    real space elements are complex, but could be reduced, as (-r) = (+r)*
!!
!! SOURCE

subroutine ftgam (wghatm,gam_qpt,gam_rpt,natom,nqpt,nrpt,qtor,coskr, sinkr)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom,nqpt,nrpt,qtor
!arrays
 real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 real(dp),intent(inout) :: gam_qpt(2,3*natom*3*natom,nqpt)
 real(dp),intent(inout) :: gam_rpt(2,3*natom*3*natom,nrpt)
 real(dp),intent(in) :: coskr(nqpt,nrpt)
 real(dp),intent(in) :: sinkr(nqpt,nrpt)

!Local variables -------------------------
!scalars
 integer :: iatom,idir,ip,iqpt,irpt,jatom,jdir
 real(dp) :: im,re
 character(len=500) :: msg

! *********************************************************************

 select case (qtor)
 case (1)
   ! Recip to real space
   gam_rpt(:,:,:) = zero
   do irpt=1,nrpt
     do iqpt=1,nqpt
       ! Get the phase factor with normalization!
       re=coskr(iqpt,irpt)
       im=sinkr(iqpt,irpt)
       do ip=1,3*natom*3*natom
         ! Real and imaginary part of the real-space gam matrices
         gam_rpt(1,ip,irpt) = gam_rpt(1,ip,irpt) + re*gam_qpt(1,ip,iqpt) + im*gam_qpt(2,ip,iqpt)
         gam_rpt(2,ip,irpt) = gam_rpt(2,ip,irpt) + re*gam_qpt(2,ip,iqpt) - im*gam_qpt(1,ip,iqpt)
       end do
     end do
   end do
   gam_rpt = gam_rpt/nqpt

 case (0)
   ! Recip space from real space
   gam_qpt(:,:,:)=zero

   do irpt=1,nrpt
     do iqpt=1,nqpt

       do iatom=1,natom
         do jatom=1,natom
           re = coskr(iqpt,irpt)*wghatm(iatom,jatom,irpt)
           im = sinkr(iqpt,irpt)*wghatm(iatom,jatom,irpt)

           do idir=1,3
             do jdir=1,3
               ! Get phase factor

               ip= jdir + (jatom-1)*3 + (idir-1)*3*natom + (iatom-1)*9*natom
               ! Real and imaginary part of the interatomic forces
               gam_qpt(1,ip,iqpt) = gam_qpt(1,ip,iqpt) + re*gam_rpt(1,ip,irpt) - im*gam_rpt(2,ip,irpt)
               !DEBUG
               gam_qpt(2,ip,iqpt) = gam_qpt(2,ip,iqpt) + im*gam_rpt(1,ip,irpt) + re*gam_rpt(2,ip,irpt)
               !ENDDEBUG
             end do ! end jdir
           end do ! end idir
         end do
       end do ! end iatom

     end do ! end iqpt
   end do ! end irpt

 case default
   write(msg,'(a,i0,a)' )'The only allowed values for qtor are 0 or 1, while qtor= ',qtor,' has been required.'
   MSG_BUG(msg)
 end select

end subroutine ftgam
!!***

!!****f* m_dynmat/ftgam_init
!!
!! NAME
!! ftgam_init
!!
!! FUNCTION
!!  Generates the sin and cos phases for the Fourier transform of the gkk matrices
!!
!! INPUTS
!! gprim = reciprocal space vectors to get cartesian coord for qpt
!! nqpt= Number of q points in the Brillouin zone
!! nrpt= Number of R points in the Big Box
!! rpt(3,nprt)= Canonical coordinates of the R points in the unit cell
!!           These coordinates are normalized (=> * acell(3)!!)
!! qpt_full(3,nqpt)= Reduced coordinates of the q vectors in reciprocal space
!!           if qtor=0 these vectors are read in the input file
!!
!! OUTPUT
!! coskr, sinkr = cosine and sine of phase factors for given r and q points
!!
!! PARENTS
!!      elphon,get_tau_k,integrate_gamma_alt,m_phgamma,mka2f,mka2f_tr
!!      mka2f_tr_lova,mkph_linwid
!!
!! CHILDREN
!!
!! SOURCE

subroutine ftgam_init (gprim,nqpt,nrpt,qpt_full,rpt,coskr, sinkr)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nqpt,nrpt
!arrays
 real(dp),intent(in) :: gprim(3,3),rpt(3,nrpt),qpt_full(3,nqpt)
 real(dp),intent(out) :: coskr(nqpt,nrpt)
 real(dp),intent(out) :: sinkr(nqpt,nrpt)

!Local variables -------------------------
!scalars
 integer :: iqpt,irpt
 real(dp) :: kr
!arrays
 real(dp) :: kk(3)

! *********************************************************************

! Prepare the phase factors
 do iqpt=1,nqpt
   ! Calculation of the k coordinates in Normalized Reciprocal coordinates
   kk(1) = qpt_full(1,iqpt)*gprim(1,1) + qpt_full(2,iqpt)*gprim(1,2) + qpt_full(3,iqpt)*gprim(1,3)
   kk(2) = qpt_full(1,iqpt)*gprim(2,1) + qpt_full(2,iqpt)*gprim(2,2) + qpt_full(3,iqpt)*gprim(2,3)
   kk(3) = qpt_full(1,iqpt)*gprim(3,1) + qpt_full(2,iqpt)*gprim(3,2) + qpt_full(3,iqpt)*gprim(3,3)
   do irpt=1,nrpt
     ! Product of k and r
     kr = kk(1)*rpt(1,irpt)+ kk(2)*rpt(2,irpt)+ kk(3)*rpt(3,irpt)
     coskr(iqpt,irpt)=cos(two_pi*kr)
     sinkr(iqpt,irpt)=sin(two_pi*kr)
   end do
 end do

end subroutine ftgam_init
!!***

end module m_dynmat
