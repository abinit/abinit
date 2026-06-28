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
!! Copyright (C) 2010-2026 ABINIT group (MG)
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
 use m_errors

 use m_numeric_tools,  only : get_trace

 implicit none

 private

 character(len=5),public :: ptgroup_names(32) =  (/ &
  "1    ",  &
  "-1   ",  &
  "2    ",  &
  "m    ",  &
  "2/m  ",  &
  "222  ",  &
  "mm2  ",  &
  "mmm  ",  &
  "4    ",  &
  "-4   ",  &
  "4/m  ",  &
  "422  ",  &
  "4mm  ",  &
  "-42m ",  &
  "4/mmm",  &
  "3    ",  &
  "-3   ",  &
  "32   ",  &
  "3m   ",  &
  "-3m  ",  &
  "6    ",  &
  "-6   ",  &
  "6/m  ",  &
  "622  ",  &
  "6mm  ",  &
  "-62m ",  &
  "6/mmm",  &
  "23   ",  &
  "m-3  ",  &
  "432  ",  &
  "-43m ",  &
  "m-3m " /)

 integer,parameter,public :: IRREPNAME_LEN=10
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

   complex(dp),allocatable :: mat(:,:,:)
   ! mat(dim,dim,nsym)
   ! The irreducible representations of the group.

   complex(dp),allocatable :: trace(:)
   ! trace(nsym)
   ! The trace of each matrix.

   contains
   procedure :: init => init_irrep
   procedure :: free => irrep_free_0d

 end type irrep_t
!!***

public :: copy_irrep

interface irrep_free
  module procedure irrep_free_0d
  module procedure irrep_free_1d
end interface irrep_free
public :: irrep_free

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

   contains

   procedure :: free => point_group_free
   procedure :: print => point_group_print
   procedure :: locate_sym => locate_sym

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

   contains
   procedure :: free => groupk_free
 end type group_k_t

contains

!----------------------------------------------------------------------

!!****f* m_ptgroups/point_group_free
!! NAME
!! point_group_free
!!
!! FUNCTION
!!  Deallocate all memory allocated in the point_group_t datatype.
!!
!! SOURCE

subroutine point_group_free(Ptg)

!Arguments ------------------------------------
 class(point_group_t),intent(inout) :: Ptg
! *********************************************************************

 ABI_SFREE(Ptg%class_ids)
 ABI_SFREE(Ptg%sym)
 ABI_SFREE(Ptg%class_names)

 if (allocated(Ptg%Irreps)) then
   call irrep_free(Ptg%Irreps)
   ABI_FREE(Ptg%Irreps)
 end if

end subroutine point_group_free
!!***

!----------------------------------------------------------------------

!!****f* m_ptgroups/point_group_print
!! NAME
!! point_group_print
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine point_group_print(Ptg, header, unit, mode_paral, prtvol)

!Arguments ------------------------------------
!scalars
 class(point_group_t),target,intent(in) :: Ptg
 integer,optional,intent(in) :: unit,prtvol
 character(len=4),optional,intent(in) :: mode_paral
 character(len=*),optional,intent(in) :: header

!Local variables-------------------------------
 integer :: my_unt,my_prtvol,irp,icls,sidx
 complex(dp) :: trace
 character(len=4) :: my_mode
 character(len=500) :: msg
 type(irrep_t),pointer :: Row
! *********************************************************************

 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 msg=' ==== Point Group Table ==== '
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 write(std_out,*)REPEAT("=",80)
 write(std_out,*)" Point group : ",TRIM(Ptg%gname)," Number of symmetries ",Ptg%nsym," Number of classes    ",Ptg%nclass

 write(std_out,"(a6)",advance="no")"Class "
 do icls=1,Ptg%nclass
   write(std_out,"('|',a10)",advance="no")Ptg%class_names(icls)
 end do
 write(std_out,"('|')",advance="no")
 write(std_out,*)" "

 write(std_out,"(a6)",advance="no")"Mult  "
 do icls=1,Ptg%nclass
   write(std_out,"('|',i10)",advance="no")Ptg%class_ids(2,icls)-Ptg%class_ids(1,icls) + 1
 end do
 write(std_out,"('|')",advance="no")
 write(std_out,*)" "

 do irp=1,SIZE(Ptg%Irreps)
   Row =>  Ptg%Irreps(irp)
   write(std_out,'(a6)',advance="no")TRIM(Row%name)

   do icls=1,Ptg%nclass
     sidx = Ptg%class_ids(1,icls)
     trace = Row%trace(sidx)
     if (ABS(AIMAG(trace)) > tol6) then
        write(std_out,"('|',(2f5.2))",advance="no")trace
      else
        write(std_out,"('|',(f10.2))",advance="no")REAL(trace)
      end if
   end do

   write(std_out,"('|')",advance="no")
   write(std_out,*)" "
 end do

 write(std_out,*)REPEAT("=",80)

end subroutine point_group_print
!!***

!----------------------------------------------------------------------

!!****f* m_defs_ptgroups/locate_sym
!! NAME
!!  locate_sym
!!
!! FUNCTION
!!  Given a symmetry operation asym, this routine returns its index in the Ptg%sym
!!  array and the index of the class it belongs to.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine locate_sym(Ptg, asym, sym_idx, cls_idx, ierr)

!Arguments ------------------------------------
!scalars
 class(point_group_t),intent(in) :: Ptg
 integer,intent(out) :: sym_idx,cls_idx
 integer,optional,intent(out) :: ierr
!arrays
 integer,intent(in) :: asym(3,3)

!Local variables-------------------------------
 integer :: isym,icls
 character(len=500) :: msg
! *********************************************************************

 sym_idx = 0
 do isym=1,Ptg%nsym
   if (ALL(asym == Ptg%sym(:,:,isym) )) then
     sym_idx = isym
     EXIT
   end if
 end do

 cls_idx = 0
 do icls=1,Ptg%nclass
   if (sym_idx >= Ptg%class_ids(1,icls) .and. sym_idx <= Ptg%class_ids(2,icls) ) then
     cls_idx = icls
     EXIT
   end if
 end do

 if (PRESENT(ierr)) ierr=0
 if (sym_idx==0 .or. cls_idx==0) then
   write(msg,'(a,9(i0,1x),3a,i1,a,i1)')&
     " Symmetry: ",asym," not found in point group table ",ch10,&
     " sym_idx= ",sym_idx, " and cls_idx= ",cls_idx
   if (PRESENT(ierr)) then
     ierr=1
     ABI_WARNING(msg)
   else
     ABI_ERROR(msg)
   end if
 end if

end subroutine locate_sym
!!***

!----------------------------------------------------------------------

!!****f* m_defs_ptgroups/irrep_free_0d
!! NAME
!! irrep_free_0d
!!
!! FUNCTION
!!  Deallocate all memory allocated in the irrep_t datatype.
!!
!! SOURCE

subroutine irrep_free_0d(Irrep)

!Arguments ------------------------------------
 class(irrep_t),intent(inout) :: Irrep
! *********************************************************************

 ABI_SFREE(Irrep%trace)
 ABI_SFREE(Irrep%mat)

end subroutine irrep_free_0d
!!***

!----------------------------------------------------------------------

!!****f* m_defs_ptgroups/irrep_free_1d
!! NAME
!! irrep_free_1d
!!
!! FUNCTION
!!  Deallocate all memory allocated in the irrep_t datatype.
!!
!! SOURCE

subroutine irrep_free_1d(Irrep)

!Arguments ------------------------------------
 class(irrep_t),intent(inout) :: Irrep(:)

!Local variables-------------------------------
 integer :: irp
! *********************************************************************

 do irp=1,SIZE(Irrep)
   call irrep_free_0d(Irrep(irp))
 end do

end subroutine irrep_free_1d
!!***

!----------------------------------------------------------------------

!!****f* m_defs_ptgroups/copy_irrep
!! NAME
!!  copy_irrep
!!
!! FUNCTION
!!  Perform a copy of a set of irrep_t datatypes. Optionally one can multiply
!!  by a phase factor.
!!
!! SOURCE

subroutine copy_irrep(In_irreps, Out_irreps, phase_fact)

!Arguments ------------------------------------
 class(irrep_t),intent(in) :: In_irreps(:)
 class(irrep_t),intent(inout) :: Out_irreps(:)
 complex(dp),optional,intent(in) :: phase_fact(:)

!Local variables-------------------------------
!scalars
 integer :: irp,dim1,dim2,in_nsym,in_dim,isym
!arrays
 complex(dp) :: my_phase_fact(In_irreps(1)%nsym)
! *********************************************************************

 !@irrep_t
 dim1 = SIZE( In_irreps)
 dim2 = SIZE(Out_irreps)
 if (dim1 /= dim2) then
   ABI_ERROR("irreps to be copied have different dimension")
 end if

 my_phase_fact=cone
 if (PRESENT(phase_fact)) then
   my_phase_fact=phase_fact
   if (SIZE(phase_fact) /= In_irreps(1)%nsym) then
     ABI_ERROR("irreps to be copied have different dimension")
   end if
 end if

 do irp=1,dim1
   in_dim  = In_irreps(irp)%dim
   in_nsym = In_irreps(irp)%nsym
   call init_irrep(Out_irreps(irp),in_nsym,in_dim)
   Out_irreps(irp)%name = In_irreps(irp)%name
   do isym=1,in_nsym
     Out_irreps(irp)%mat(:,:,isym) = In_irreps(irp)%mat(:,:,isym) * my_phase_fact(isym)
     Out_irreps(irp)%trace(isym) = get_trace(Out_irreps(irp)%mat(:,:,isym))
   end do
 end do

end subroutine copy_irrep
!!***

!----------------------------------------------------------------------

!!****f* m_defs_ptgroups/init_irrep
!! NAME
!!  alloc_irrep
!!
!! FUNCTION
!!  Initialize an instance of the irrep_t datatype.
!!
!! INPUTS
!!  nsym=The number of symmetries.
!!  irr_dim=The dimension of the irrep.
!!  [irr_name]=The name of theirrep. "???" is used if not given
!!
!! OUTPUT
!!  Irrep<irrep_t>=
!!
!! SOURCE

subroutine init_irrep(Irrep, nsym, irr_dim, irr_name)

!Arguments ------------------------------------
!scalars
 class(irrep_t),intent(inout) :: Irrep
 integer,intent(in) :: nsym
!arrays
 integer,intent(in) :: irr_dim
 character(len=*),optional,intent(in) :: irr_name

!Local variables-------------------------------
 !character(len=500) :: msg
! *********************************************************************

 !@irrep_t
 Irrep%dim  = irr_dim
 Irrep%nsym = nsym
 Irrep%name = "???"
 if (present(irr_name)) Irrep%name = irr_name

 ABI_CALLOC(Irrep%mat,(irr_dim,irr_dim,nsym))
 ABI_CALLOC(Irrep%trace,(nsym))

end subroutine init_irrep
!!***

!----------------------------------------------------------------------

!!****f* m_defs_ptgroups/groupk_free
!! NAME
!!  groupk_free
!!
!! FUNCTION
!!  Deallocate dynamic memory.
!!
!! SOURCE

subroutine groupk_free(Gk)

!Arguments ------------------------------------
 class(group_k_t),intent(inout) :: Gk
! *************************************************************************

! integer
 ABI_SFREE(Gk%class_ids)
 ABI_SFREE(Gk%sym)

!real
 ABI_SFREE(Gk%tnons)

!character
 ABI_SFREE(Gk%class_names)

!type
 if (allocated(Gk%Irreps)) then
   call irrep_free(Gk%Irreps)
   ABI_FREE(Gk%Irreps)
 end if

end subroutine groupk_free
!!***

end module m_defs_ptgroups
!!***
