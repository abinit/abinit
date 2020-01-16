!!****m* ABINIT/m_nesting
!! NAME
!!  m_nesting
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2008-2019 ABINIT group (MG, MJV)
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

MODULE m_nesting

 use defs_basis
 use m_errors
 use m_abicore
 use m_krank
 use m_sort

 use m_numeric_tools,  only : wrap2_zero_one, interpol3d
 use m_io_tools,       only : open_file
 use m_bz_mesh,        only : make_path
 use m_pptools,        only : printxsf

 implicit none

 private

 public :: bfactor
 public :: mknesting
 public :: outnesting
!!***

!----------------------------------------------------------------------

CONTAINS  !=============================================================================
!!***

!!****f* m_nesting/bfactor
!! NAME
!! bfactor
!!
!! FUNCTION
!! Calculate the nesting factor
!!
!! INPUTS
!!  nkptfull = number of k-points in full grid
!!  kptfull(3,nkptfull) = k-point grid
!!  nqpt = number of qpoints
!!  qpt(3,nqpt) = q-point grid (must be a subgrid of the k grid),
!!                the nesting factor will be calculated for each q point in this array
!!  nkpt = eventually reduced number of k-points
!!  weight(nband,nkpt) =  integration weights for each k-point and band (NOT NORMALISED!!!)
!!  nband = number of bands
!!
!! OUTPUT
!!  nestfactor(nqpt) = array containing the nesting factor values
!!
!! NOTES
!! Inspired to nmsq_gam_sumfs and mkqptequiv
!!  TODO : better use of symmetries to reduce the computational effort
!! Must be called with kpt = full grid! Reduction by symmetry is not possible for q-dependent quantities (or not easy :)
!!
!! PARENTS
!!      m_nesting,outelph
!!
!! CHILDREN
!!      make_path,printxsf,wrap2_zero_one
!!
!! SOURCE

subroutine bfactor(nkptfull,kptfull,nqpt,qpt,krank,nkpt,weight,nband,nestfactor)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nband,nkptfull,nqpt,nkpt
!arrays
 real(dp),intent(in) :: kptfull(3,nkptfull),qpt(3,nqpt),weight(nband,nkpt)
 real(dp),intent(out) :: nestfactor(nqpt)
 type(krank_t), intent(in) :: krank

!Local variables-------------------------------
!scalars
 integer :: ib1,ib2,ikplusq_irr,ikpt
 integer :: irank_kpt,ikpt_irr,iqpt,symrank_kpt
 real(dp) :: w1,w2
 !character(len=500) :: message
!arrays
 real(dp) :: kptpq(3)

! *************************************************************************

 nestfactor(:)=zero

 do iqpt=1,nqpt
   do ikpt=1,nkptfull
     irank_kpt = krank%get_rank(kptfull(:,ikpt))
     ikpt_irr = krank%invrank(irank_kpt)

     kptpq(:) = kptfull(:,ikpt) + qpt(:,iqpt)
     symrank_kpt = krank%get_rank(kptpq)

     ikplusq_irr = krank%invrank(symrank_kpt)
     if (ikplusq_irr == -1) then
       MSG_ERROR('It looks like no kpoint equiv to k+q!')
     end if

     do ib1=1,nband
       w1 = weight(ib1,ikpt_irr) !weight for distance from the Fermi surface
       if (w1 < tol6 ) cycle
       do ib2=1,nband
         w2 = weight(ib2,ikplusq_irr) !weight for distance from the Fermi surface
         if (w1 < tol6 ) cycle
         nestfactor(iqpt) = nestfactor(iqpt) + w1*w2
       end do !ib2
     end do !ib1

   end do !ikpt
 end do !iqpt

!need prefactor of (1/nkptfull) for normalisation of integration
 nestfactor(:) = (one/nkptfull) * nestfactor(:)

end subroutine bfactor
!!***

!----------------------------------------------------------------------

!!****f* m_nesting/mknesting
!! NAME
!! mknesting
!!
!! FUNCTION
!!  Calculate the nesting factor over the dense k-grid, interpolate the values along a given q path
!!  and write the data on file in the X-Y format or in the XCrysden format (XSF)
!!
!! INPUTS
!!  nkpt = number of k points
!!  kpt(3,nkpt) = k points
!!  nkx, nky, nkz = number of k-point along each direction
!!  nband = number of bands to be considered in the calculation
!!  weight(nband,nkpt) =  integration weights for each k-point and band
!!  nqpath = number of points requested along the trajectory
!!  qpath_vertices = vertices of the reciprocal space trajectory
!!  base_name = prefix of the output file
!!  gprimd(3,3) dimensional reciprocal lattice vectors
!!  gmet = metric in reciprocal space
!!  prtnest = flags governing the format of the output file
!!
!! OUTPUT
!!   Write data to file.
!!
!! PARENTS
!!      elphon,m_ebands
!!
!! CHILDREN
!!      make_path,printxsf,wrap2_zero_one
!!
!! SOURCE

subroutine mknesting(nkpt,kpt,kptrlatt,nband,weight,nqpath,&
& qpath_vertices,nqptfull,qptfull,base_name,gprimd,gmet,prtnest,qptrlatt,&
& nsym,symrec) ! optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nband,nkpt,nqpath,prtnest
 integer, intent(in) :: nqptfull
 integer, intent(in), optional :: nsym
 character(len=*),intent(in) :: base_name
!arrays
 integer,intent(in) :: kptrlatt(3,3)
 integer,intent(in),optional :: symrec(3,3,*)
 real(dp),intent(in) :: gprimd(3,3),kpt(3,nkpt)
 real(dp),intent(in) :: qptfull(3,nqptfull)
 real(dp),intent(in) :: gmet(3,3)
 real(dp),intent(in) :: qpath_vertices(3,nqpath)
 real(dp),intent(in) :: weight(nband,nkpt)
 integer,intent(in)  :: qptrlatt(3,3)

!Local variables-------------------------------
!scalars
 integer :: ikpt,jkpt
 integer :: ik1, ik2, ik3, nkptfull
 character(len=500) :: message
 type(krank_t) :: krank
!arrays
 integer,allocatable :: tmprank(:),ktable(:)
 character(len=fnlen) :: tmpname
 real(dp),allocatable :: nestfactor(:),nestordered(:)
 real(dp), allocatable :: kptfull(:,:)

! *************************************************************************

 if (kptrlatt(1,2) /= 0 .or. kptrlatt(1,3) /= 0 .or. kptrlatt(2,1) /= 0 .or. &
    kptrlatt(2,3) /= 0 .or. kptrlatt(3,1) /= 0 .or. kptrlatt(3,2) /= 0 ) then
   write (message,'(4a)')&
    'kptrlatt should be diagonal in order to calculate the nesting factor,',ch10,&
    'skipping the nesting factor calculation ',ch10
   MSG_WARNING(message)
   return
 end if

 if (prtnest /= 1 .and. prtnest /= 2) then
   MSG_BUG('prtnest should be 1 or 2')
 end if

 !write(message,'(a,9(i0,1x))')' mknesting : kptrlatt = ',kptrlatt
 !call wrtout(std_out,message,'COLL')

 nkptfull = kptrlatt(1,1)*kptrlatt(2,2)*kptrlatt(3,3)
 ABI_MALLOC(nestordered,(nkptfull))
 nestordered(:)=zero
 ABI_MALLOC(kptfull,(3,nkptfull))

 ikpt = 0
 do ik3 = 0, kptrlatt(3,3)-1
   do ik2 = 0, kptrlatt(2,2)-1
     do ik1 = 0, kptrlatt(1,1)-1
       ikpt = ikpt+1
       kptfull(:,ikpt) = (/dble(ik1)/dble(kptrlatt(1,1)), dble(ik2)/dble(kptrlatt(2,2)),dble(ik3)/dble(kptrlatt(3,3))/)
     end do
   end do
 end do

!NOTE: input weights are not normalised, the normalisation factor in introduced in bfactor
!new version now puts kptfull in correct order before bfactor, so no need to re-order...
 if (present(symrec)) then
   ABI_CHECK(present(nsym), "error - provide nsym and symrec arguments together")
   krank = krank_new(nkpt, kpt, nsym=nsym, symrec=symrec)
 else
   krank = krank_new(nkpt, kpt)
 end if

 call bfactor(nkptfull,kptfull,nkptfull,kptfull,krank,nkpt,weight,nband,nestordered)

!================================================================================================
!use linear interpolation to plot the bfactor along the given q-path
!1) order the kpoints of the grid putting them in increasing x, then y, then z (FORTRAN convention)
!2) make table from input kpts to ordered kpts
!3) perform interpolation
!================================================================================================

 call outnesting(base_name,gmet,gprimd,kptrlatt,nestordered,nkptfull,nqpath,prtnest,qpath_vertices)
 ABI_FREE(nestordered)
!
!now do the same, but for the nesting factor over the phonon qpoints only
!
 ABI_MALLOC(nestfactor,(nqptfull))
 call bfactor(nkptfull,kptfull,nqptfull,qptfull,krank,nkpt,weight,nband,nestfactor)

 call krank%free()
 ABI_FREE(kptfull)

 krank = krank_new(nqptfull, qptfull)

 ABI_MALLOC(ktable,(nqptfull))
 do ikpt=1,nqptfull
   ktable(ikpt) = ikpt
 end do

 ABI_MALLOC(tmprank, (nqptfull))
 do ikpt=1,nqptfull
   tmprank(ikpt) = krank%get_rank(qptfull(:,ikpt))
 end do
 call sort_int(nqptfull, tmprank, ktable)
 ABI_FREE(tmprank)
 call krank%free()

!fill the datagrid for the nesting factor using the Fortran convention and the conventional unit cell
!NOTE: the Fortran convention is a must if we want to plot the data
!in the BXSF format, useful for the linear interpolation since we use interpol3d.F90
 ABI_MALLOC(nestordered,(nqptfull))
 nestordered(:)=zero
 do jkpt=1,nqptfull
   ikpt = ktable(jkpt)
   nestordered(ikpt)=nestfactor(jkpt)
 end do
 ABI_FREE(nestfactor)
 ABI_FREE(ktable)

 tmpname = trim(base_name)//"kplusq"
 call outnesting(tmpname,gmet,gprimd,qptrlatt,nestordered,nqptfull,nqpath,prtnest,qpath_vertices)

 ABI_FREE(nestordered)

end subroutine mknesting
!!***

!----------------------------------------------------------------------

!!****f* m_nesting/outnesting
!! NAME
!! outnesting
!!
!! FUNCTION
!!  Write ou the nesting factors calculated in mknesting
!!  Data on file in the X-Y format (prtnest 1) or
!!  in the XCrysden format (XSF)   (prtnest 2)
!!
!! INPUTS
!!  base_name = prefix of the output file
!!  gmet = metric in reciprocal space
!!  gprimd(3,3) dimensional reciprocal lattice vectors
!!  kptrlatt(3,3) basis vectors for k-grid
!!  nestordered = nesting function on full grid, points ordered in x, then y, then z
!!  nkpt = number of k points
!!  nqpath = number of points requested along the trajectory
!!  prtnest = flags governing the format of the output file
!!  qpath_vertices = vertices of the reciprocal space trajectory
!!
!! OUTPUT
!!  only write to file
!!
!! PARENTS
!!      m_nesting
!!
!! CHILDREN
!!      make_path,printxsf,wrap2_zero_one
!!
!! SOURCE

subroutine outnesting(base_name,gmet,gprimd,kptrlatt,nestordered,nkpt,nqpath,prtnest,qpath_vertices)

!Arguments ------------------------------------
 integer,intent(in) :: nqpath,prtnest,nkpt
 character(len=*),intent(in) :: base_name
!arrays
 integer,intent(in) :: kptrlatt(3,3)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(in) :: gmet(3,3)
 real(dp),intent(in) :: qpath_vertices(3,nqpath)
 real(dp),intent(in) :: nestordered(nkpt)

!Local variables-------------------------------
!scalars
 integer :: unit_nest,nkx,nky,nkz
 integer :: indx,ii,ipoint,npt_tot,realrecip
 character(len=fnlen) :: fname
 character(len=500) :: message
 real(dp) :: res, kval
!arrays
 integer :: ndiv(nqpath-1)
 real(dp),allocatable :: finepath(:,:)
 real(dp) :: tmpkpt(3)
 real(dp) :: origin(3),qpt(3)
! dummy variables for call to printxsf
 integer :: natom, ntypat, typat(1)
 real(dp) :: xcart (3,1), znucl(1)

! *************************************************************************

!===================================================================
!Definition of the q path along which ph linwid will be interpolated
!===================================================================
 call make_path(nqpath,qpath_vertices,gmet,'G',20,ndiv,npt_tot,finepath)

 nkx=kptrlatt(1,1)
 nky=kptrlatt(2,2)
 nkz=kptrlatt(3,3)

 if (nkpt /= nkx*nky*nkz) then
   write(message,'(a,9(i0,1x),2x,i0)')' Wrong input value for kptrlatt  ',kptrlatt, nkpt
   MSG_BUG(message)
 end if

 ! open output file and write header
 if (open_file(base_name,message,newunit=unit_nest,status="unknown",form="formatted",action="write") /= 0) then
    MSG_ERROR(message)
 end if

 write (unit_nest,'(a)')'#'
 write (unit_nest,'(a)')'# ABINIT package : Nesting factor file'
 write (unit_nest,'(a)')'#'
 write (unit_nest,'(a,i10,a)')'# Nesting factor calculated on ',npt_tot,' Q-points'
 write (unit_nest,'(a)')'# Description of the Q-path :'
 write (unit_nest,'(a,i10)')'# Number of line segments = ',nqpath-1
 write (unit_nest,'(a)')'# Vertices of the Q-path and corresponding index = '
 indx=1
 do ii=1,nqpath
   write (unit_nest,'(a,3(E16.6,1x),i8)')'#  ',qpath_vertices(:,ii),indx
   if(ii<nqpath) indx=indx+ndiv(ii)
 end do
 write (unit_nest,'(a)')'#'

!Get qpoint along the q-path from finepath and interpolate the nesting factor
 indx=1

 do ipoint=1, npt_tot
   qpt(:) = finepath(:,ipoint)
   call wrap2_zero_one(qpt(1),tmpkpt(1),res)
   call wrap2_zero_one(qpt(2),tmpkpt(2),res)
   call wrap2_zero_one(qpt(3),tmpkpt(3),res)

   kval = interpol3d(tmpkpt,nkx,nky,nkz,nestordered)

   write(unit_nest,'(i5,18e16.5)')indx,kval
   indx = indx+1
 end do

 close (unit_nest)
 ABI_FREE(finepath)

 if (prtnest==2) then !write also the nest factor in the XSF format
   fname=trim(base_name) // '_NEST_XSF'

   if (open_file(fname,message,newunit=unit_nest,status="unknown",form="formatted",action="write") /= 0) then
      MSG_ERROR(message)
   end if

   origin(:)=zero
   realrecip=1 !reciprocal space
   natom = 1
   ntypat = 1
   typat = (/1/)
   xcart = reshape ((/zero, zero, zero/), (/3,1/))
   znucl = (/one/)
   call printxsf(nkx,nky,nkz,nestordered,gprimd,origin,natom, ntypat, typat, xcart, znucl, unit_nest,realrecip)

   close (unit_nest)
 end if

end subroutine outnesting
!!***

END MODULE m_nesting
!!***
