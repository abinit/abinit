!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_crystal_io
!! NAME
!! m_crystal_io
!!
!! FUNCTION
!! Module containing the methods used to do IO on crystal objects.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2018 ABINIT group (MG, YP, DC)
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

MODULE m_crystal_io

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_crystal
 use m_atomdata
 use m_nctk
 use iso_c_binding
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use defs_abitypes,    only : hdr_type

 implicit none

 private
!!***

 public :: crystal_from_hdr      ! Initialize the object from the abinit header.
 !public :: crystal_from_dtset   ! Initialize the object from the abinit dataset.
 !public :: crystal_ncwrite       ! Dump the object in a netcdf file associated to a ncid.
 !public :: crystal_ncwrite_path  ! Dump the object to file.


CONTAINS

!!****f* m_crystal_io/crystal_from_hdr
!! NAME
!!  crystal_from_hdr
!!
!! FUNCTION
!!  Initializes a crystal_t data type starting from the abinit header.
!!
!! INPUTS
!!  hdr<hdr_type>=the abinit header
!!  timrev ==2 => take advantage of time-reversal symmetry
!!         ==1 ==> do not use time-reversal symmetry
!!  remove_inv [optional]= if .TRUE. the inversion symmetry is removed from the set of operations
!!  even if it is present in the header
!!
!! OUTPUT
!!  cryst<crystal_t>= the data type filled with data reported in the abinit header
!!
!! TODO
!!  Add information on the use of time-reversal in the Abinit header.
!!
!! PARENTS
!!      cut3d,eph,fold2Bloch,gstate,m_ddk,m_dvdb,m_ioarr,m_iowf,m_wfd,m_wfk
!!      mlwfovlp_qp,mrgscr,setup_bse,setup_screening,setup_sigma,wfk_analyze
!!
!! CHILDREN
!!      atomdata_from_znucl,crystal_init
!!
!! SOURCE

subroutine crystal_from_hdr(cryst,hdr,timrev,remove_inv)

 implicit none

!Arguments ------------------------------------
 type(hdr_type),intent(in) :: hdr
 type(crystal_t),intent(out) :: cryst
 integer,intent(in) :: timrev
 logical,optional,intent(in) :: remove_inv

!Local variables-------------------------------
 integer :: space_group
 logical :: rinv,use_antiferro
! *********************************************************************

 rinv=.FALSE.; if (PRESENT(remove_inv)) rinv=remove_inv
 use_antiferro = (hdr%nspden==2.and.hdr%nsppol==1)

 ! Consistency check
 ABI_CHECK(any(timrev == [1, 2]),"timrev should be in (1|2)")
 if (use_antiferro) then
   ABI_CHECK(ANY(hdr%symafm==-1),"Wrong nspden, nsppol, symafm.")
 end if

 space_group=0 !FIXME not known

 call crystal_init(hdr%amu,cryst,space_group,hdr%natom,hdr%npsp,hdr%ntypat,hdr%nsym,hdr%rprimd,hdr%typat,hdr%xred,&
& hdr%zionpsp,hdr%znuclpsp,timrev,use_antiferro,rinv,hdr%title,&
& symrel=hdr%symrel,tnons=hdr%tnons,symafm=hdr%symafm) ! Optional

end subroutine crystal_from_hdr
!!***

end module m_crystal_io
!!***
