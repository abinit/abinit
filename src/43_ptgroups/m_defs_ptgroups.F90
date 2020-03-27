!!****m* ABINIT/m_defs_ptgroups
!! NAME
!! m_defs_ptgroups
!!
!! FUNCTION
!!  This module contains the definition of the point_group_t datatype used to
!!  represent one of the 32 different point groups. It also provides the definition
!!  of the irrept_t structure that is used to store one of the irreducible 
!!  representations of the group
!!
!! NOTES
!!  Methods "bound" to the point_group_t datatype are defined in the
!!  separate module m_ptgroups in order to avoid cyclic dependencies. 
!!  The automatically generated ptg_* routines indeed use irrept_t and 
!!  cannot be #included here.  The size of the final module 
!!  is indeed huge (>4000 lines) and several compilers are not able to
!!  compile the module in reasonable time when -O2 is used.
!!
!! COPYRIGHT
!! Copyright (C) 2010-2020 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_defs_ptgroups

 use defs_basis
 use m_abicore

 implicit none

 private

 character(len=5),public :: ptgroup_names(32) =  (/ &
&  "1    ",&   
&  "-1   ",&
&  "2    ",&
&  "m    ",&
&  "2/m  ",&
&  "222  ",&
&  "mm2  ",&
&  "mmm  ",&
&  "4    ",&
&  "-4   ",&
&  "4/m  ",&
&  "422  ",&
&  "4mm  ",&
&  "-42m ",&
&  "4/mmm",&
&  "3    ",&
&  "-3   ",&
&  "32   ",&
&  "3m   ",&
&  "-3m  ",&
&  "6    ",&
&  "-6   ",&
&  "6/m  ",& 
&  "622  ",&
&  "6mm  ",&
&  "-62m ",&
&  "6/mmm",&
&  "23   ",&
&  "m-3  ",&
&  "432  ",&
&  "-43m ",&
&  "m-3m " /)

 integer,parameter,public :: IRREPNAME_LEN=10

!Structures
!!***

!!****t* m_defs_ptgroups/irrep_t
!! NAME
!! irrep_t
!! 
!! FUNCTION
!!  Datatype representing an (irreducible) representation of the group.
!! 
!! SOURCE

 type,public :: irrep_t

   integer :: dim
   ! The dimension of the irreducible representation.

   integer :: nsym
   ! The number of symmetries in the group. 

   character(len=IRREPNAME_LEN) :: name="???"
   ! The name of the irreducible representation.

   complex(dpc),allocatable :: mat(:,:,:)
   ! mat(dim,dim,nsym)
   ! The irreducible representations of the group.

   complex(dpc),allocatable :: trace(:) 
   ! trace(nsym) 
   ! The trace of each matrix.

 end type irrep_t
!!***

!----------------------------------------------------------------------

!!****t* m_defs_ptgroups/point_group_t
!! NAME
!! point_group_t
!! 
!! FUNCTION
!!  Datatype used to collect data concerning one of the 32 different point groups.
!! 
!! SOURCE

 type,public :: point_group_t

   !integer :: numspg 
   ! Number of the space group.

   integer :: nsym
   ! Number of symmetries in the point group.

   integer :: nclass
   ! Number of classes (equals the number of irreducible representations).

   character(len=5) :: gname
   ! Point group name.

   !character(len=10) :: setting
   ! Information about the space group setting (Table, Standard)

   integer,allocatable :: class_ids(:,:) 
   ! class_ids(2,nclass)
   ! (1,icl) = index of the first symmetry of class icl
   ! (2,icl) = index of the last symmetry of class icl
   ! Note that symmetries in the sym array are packed in classes.

   integer,allocatable :: sym(:,:,:) 
   ! The symmetry operations packed in classes.
   ! NB: operations are referred to the standard coordinate system.
   ! Page 815-816 of International Tables for crystallography Vol.A.

   !$integer,allocatable :: symafm(:)
   ! symafm(nsym)
   ! AFM part of the symmetry operation

   !$real(dp),allocatable :: tnons(:,:) 
   ! tnons(3,nsym)
   ! fractional translations.

   character(len=5),allocatable :: class_names(:) 
   ! class_names(nclass) 
   ! The name of each class.

   type(irrep_t),allocatable :: Irreps(:) 
   ! Irreps(nclass)
   ! Array storing the irreducible representations of the point group.
   ! Initialized from the tables downloaded from the Bilbao server.

 end type point_group_t
!!***

!----------------------------------------------------------------------

!!****t* m_defs_ptgroups/group_k_t
!! NAME
!! group_k_t
!! 
!! FUNCTION
!!  Datatype used to collect data on the little group.
!! 
!! SOURCE

 type,public :: group_k_t

   integer :: spgroup
   ! ITA space group number.

   integer :: nsym
   ! Number of symmetries in the little group.

   integer :: nclass
   ! Number of classes (equals the number of irreducible representations).

   integer,allocatable :: class_ids(:,:) 
   ! class_ids(2,nclass)
   ! (1,icl) = index of the first symmetry of class icl
   ! (2,icl) = index of the last symmetry of class icl
   ! Note that symmetries in sym are packed in classes.

   integer,allocatable :: sym(:,:,:) 
   ! sym(3,3,nsym)
   ! The symmetry operations of the little group packed in classes.
   ! NB: operations are referred to the standard coordinate system.
   ! Page 815-816 of Internationat Tables for crystallography Vol.A.

   real(dp) :: point(3)
   ! The point referred to the standard coordinate system.

   real(dp),allocatable :: tnons(:,:)
   ! tnons(3,nsym)
   ! Fractional translations of the little group.
   ! NB: operations are referred to the standard coordinate system.
   ! Page 815-816 of Internationat Tables for crystallography Vol.A.

   character(len=5),allocatable :: class_names(:) 
   ! class_names(nclass) 
   ! The name of each class.

   type(irrep_t),allocatable :: Irreps(:)
   ! Irreps(nclass)
   ! The set of irreducible representations of the point group.

 end type group_k_t

!CONTAINS  !===========================================================

END MODULE m_defs_ptgroups
!!***
