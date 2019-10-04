!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_paw_an
!! NAME
!!  m_paw_an
!!
!! FUNCTION
!!  This module contains the definition of the paw_an_type structured datatype,
!!  as well as related functions and methods.
!!  paw_an_type variables contain various arrays given on ANgular mesh or ANgular moments
!!  for a given atom.
!!
!! COPYRIGHT
!! Copyright (C) 2013-2019 ABINIT group (MT, FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  FOR DEVELOPPERS: in order to preserve the portability of libPAW library,
!!  please consult ~abinit/src/??_libpaw/libpaw-coding-rules.txt
!!
!! SOURCE

#include "libpaw.h"

MODULE m_paw_an

 USE_DEFS
 USE_MSG_HANDLING
 USE_MPI_WRAPPERS
 USE_MEMORY_PROFILING

 use m_paral_atom, only : get_my_atmtab, free_my_atmtab, get_my_natom
 use m_pawang,     only : pawang_type
 use m_pawtab,     only : pawtab_type

 implicit none

 private

!public procedures.
 public :: paw_an_init
 public :: paw_an_free
 public :: paw_an_nullify
 public :: paw_an_copy
 public :: paw_an_print
 public :: paw_an_gather
 public :: paw_an_redistribute
 public :: paw_an_reset_flags

!private procedures.
 private :: paw_an_isendreceive_getbuffer
 private :: paw_an_isendreceive_fillbuffer
!!***

!----------------------------------------------------------------------

!!****t* m_paw_an/paw_an_type
!! NAME
!! paw_an_type
!!
!! FUNCTION
!! For PAW, various arrays given on ANgular mesh or ANgular moments
!!
!! SOURCE

 type,public :: paw_an_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

!Integer scalars

  integer :: angl_size
   ! Dimension of paw angular mesh (angl_size=ntheta*nphi)

  integer :: cplex
   ! cplex=1 if potentials/densities are real, 2 if they are complex

  integer :: has_kxc
   ! set to 1 if xc kernels kxc1 and kxct1 are allocated and used
   !        2 if they are already computed

  integer :: has_k3xc
   ! set to 1 if xc kernel derivatives k3xc1 and k3xct1 are allocated and used
   !        2 if it is already computed

  integer :: has_vhartree
   ! set to 1 if vh1 and vht1 are allocated and used
   !        2 if they are already computed

  integer :: has_vxc
   ! set to 1 if vxc1 and vxct1 are allocated and used
   !        2 if they are already computed

  integer :: has_vxctau
   ! set to 1 if vxctau1 and vxcttau1 are allocated and used
   !        2 if they are already computed

  integer :: has_vxcval
   ! set to 1 if vxc1_val and vxct1_val are allocated and used
   !        2 if they are already computed

  integer :: has_vxc_ex
   ! set to 1 if vxc_ex and is allocated and used
   !        2 if it is already computed

  integer :: itypat
   ! itypat=type of the atom

  integer :: lm_size
   ! lm_size=(l_size)**2
   ! l is Maximum value of l+1 leading to non zero Gaunt coeffs (l_size=2*l_max+1)

  integer :: mesh_size
   ! Dimension of radial mesh for arrays contained in this paw_an datastructure
   ! May be different from pawrad%mesh_size

  integer :: nkxc1
   ! number of independent components of Kxc1 and Kxct1
   ! (usually 3 for LDA, 23 for GGA)

  integer :: nk3xc1
   ! number of independent components of K3xc1 and K3xct1
   ! (usually 4 for LDA, not available for GGA)

  integer :: nspden
   ! Number of spin-density components

!Logical arrays

  logical, allocatable :: lmselect(:)
   ! lmselect(lm_size)
   ! lmselect(ilm)=select the non-zero LM-moments of "one-center" densities/potentials

!Real (real(dp)) arrays

  real(dp), allocatable :: kxc1 (:,:,:)
   ! kxc1(cplex*mesh_size,lm_size or angl_size,nkxc1)
   ! Gives xc kernel inside the sphere
   !   (theta,phi) values of kernel if pawxcdev=0
   !   LM-moments of kernel if pawxcdev/=0

  real(dp), allocatable :: kxct1 (:,:,:)
   ! kxct1(cplex*mesh_size,lm_size or angl_size,nkxc1)
   ! Gives xc pseudo kernel inside the sphere
   !   (theta,phi) values of kernel if pawxcdev=0
   !   LM-moments of kernel if pawxcdev/=0

  real(dp), allocatable :: k3xc1 (:,:,:)
   ! k3xc1(cplex*mesh_size,lm_size or angl_size,nk3xc1)
   ! Gives xc kernel derivative inside the sphere
   !   (theta,phi) values of kernel derivative if pawxcdev=0
   !   LM-moments of kernel derivative if pawxcdev/=0 => NOT AVAILABLE YET

  real(dp), allocatable :: k3xct1 (:,:,:)
   ! k3xct1(cplex*mesh_size,lm_size or angl_size,nk3xc1)
   ! Gives xc pseudo kernel derivative inside the sphere
   !   (theta,phi) values of kernel derivative if pawxcdev=0
   !   LM-moments of kernel derivative if pawxcdev/=0 => NOT AVAILABLE YET

  real(dp), allocatable :: vh1 (:,:,:)
   ! vh1(cplex*mesh_size,lm_size,nspden)
   ! Gives Hartree potential LM-moments inside the sphere

  real(dp), allocatable :: vht1 (:,:,:)
   ! vht1(cplex*mesh_size,lm_size,nspden)
   ! Gives Hartree tilde potential LM-moments inside the sphere

  real(dp), allocatable :: vxc1 (:,:,:)
   ! vxc1(cplex*mesh_size,lm_size or angl_size,nspden)
   ! Gives xc potential inside the sphere
   !   (theta,phi) values of potential if pawxcdev=0
   !   LM-moments of potential if pawxcdev/=0

  real(dp), allocatable :: vxctau1 (:,:,:)
   ! vxctau1(cplex*mesh_size,lm_size or angl_size,nspden)
   ! Gives xc potential inside the sphere
   !   (theta,phi) values of potential if pawxcdev=0
   !   LM-moments of potential if pawxcdev/=0

  real(dp), allocatable :: vxc1_val (:,:,:)
   ! vxc1_val(cplex*mesh_size,lm_size or angl_size,nspden) (Usually real, Mainly used for GW)
   ! Gives xc potential inside the sphere arising from valence only electrons
   !   (theta,phi) values of potential if pawxcdev=0
   !   LM-moments of potential if pawxcdev/=0

  real(dp), allocatable :: vxct1 (:,:,:)
   ! vxct1(cplex*mesh_size,angl_size,nspden)
   ! Gives xc pseudo potential inside the sphere
   !   (theta,phi) values of potential if pawxcdev=0
   !   LM-moments of potential if pawxcdev/=0

 real(dp), allocatable :: vxcttau1 (:,:,:)
   ! vxcttau1(cplex*mesh_size,angl_size,nspden)
   ! Gives xc pseudo potential inside the sphere
   !   (theta,phi) values of potential if pawxcdev=0
   !   LM-moments of potential if pawxcdev/=0

  real(dp), allocatable :: vxct1_val (:,:,:)
   ! vxct1_val(cplex*mesh_size,angl_size,nspden) (Usually real, Mainly used for GW)
   ! Gives xc pseudo potential inside the sphere
   !   (theta,phi) values of potential if pawxcdev=0
   !   LM-moments of potential if pawxcdev/=0

  real(dp), allocatable :: vxc_ex (:,:,:)
   ! vxc_ex(cplex*mesh_size,angl_size,nspden)
   ! Gives xc  potential for local exact exchange inside the sphere
   !   (theta,phi) values of potential if pawxcdev=0
   !   LM-moments of potential if pawxcdev/=0

 end type paw_an_type
!!***

CONTAINS

!===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_paw_an/paw_an_init
!! NAME
!!  paw_an_init
!!
!! FUNCTION
!!  Initialize a paw_an data type.
!!
!! INPUTS
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!
!! SIDE EFFECTS
!!  Paw_an(:)<type(paw_an_type)>=PAW arrays given on ANgular mesh or ANgular moments.
!!                               Initialized in output
!!
!! PARENTS
!!      bethe_salpeter,dfpt_scfcv,paw_qpscgw,respfn,scfcv,screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine paw_an_init(Paw_an,natom,ntypat,nkxc1,nk3xc1,nspden,cplex,pawxcdev,typat,Pawang,Pawtab,&
&          has_vhartree,has_vxc,has_vxctau,has_vxcval,has_kxc,has_k3xc,has_vxc_ex, & ! optional arguments
&          mpi_atmtab,comm_atom) ! optional arguments (parallelism)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nkxc1,nk3xc1,ntypat,cplex,nspden,pawxcdev
 integer,optional,intent(in) :: has_vhartree,has_vxc,has_vxctau,has_vxcval,has_kxc,has_k3xc,has_vxc_ex
 integer,optional,intent(in) :: comm_atom
!arrays
 integer,intent(in) :: typat(natom)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 type(Pawang_type),intent(in) :: Pawang
 type(Pawtab_type),intent(in) :: Pawtab(ntypat)
 type(Paw_an_type),intent(inout) :: Paw_an(:)

!Local variables-------------------------------
!scalars
 integer :: iat,iat1,itypat,lm_size,my_comm_atom,my_natom,v_size
 logical :: my_atmtab_allocated,paral_atom
!arrays
 integer,pointer :: my_atmtab(:)

! *************************************************************************

!@Paw_an_type

!Set up parallelism over atoms
 my_natom=size(Paw_an);if (my_natom==0) return
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

 do iat=1,my_natom
  iat1=iat;if (paral_atom) iat1=my_atmtab(iat)
  itypat=typat(iat1)

  lm_size                =Pawtab(itypat)%lcut_size**2
  Paw_an(iat)%angl_size  =Pawang%angl_size
  Paw_an(iat)%cplex      =cplex
  Paw_an(iat)%itypat     =itypat
  Paw_an(iat)%lm_size    =lm_size
  Paw_an(iat)%mesh_size  =Pawtab(itypat)%mesh_size
  Paw_an(iat)%nkxc1      =nkxc1
  Paw_an(iat)%nk3xc1     =nk3xc1
  Paw_an(iat)%nspden     =nspden

  ! === Non-zero LM-moments of "one-center" densities/potentials ===
  ! * Filled in pawdenpot.
  LIBPAW_ALLOCATE(Paw_an(iat)%lmselect,(lm_size))

  v_size=Paw_an(iat)%lm_size ; if (pawxcdev==0) v_size=Paw_an(iat)%angl_size

 ! === XC potential inside the sphere ===
 ! * LM-moments of potential if pawxcdev/=0
 ! * (theta,phi) values of potential if pawxcdev=0
  Paw_an(iat)%has_vxc=0
  if (PRESENT(has_vxc)) then
   if (has_vxc>0) then
    Paw_an(iat)%has_vxc=1
    LIBPAW_ALLOCATE(Paw_an(iat)%vxc1 ,(cplex*Paw_an(iat)%mesh_size,v_size,nspden))
    LIBPAW_ALLOCATE(Paw_an(iat)%vxct1,(cplex*Paw_an(iat)%mesh_size,v_size,nspden))
    Paw_an(iat)%vxc1=zero;Paw_an(iat)%vxct1=zero
   end if
  end if

 Paw_an(iat)%has_vxctau=0
 if (PRESENT(has_vxctau)) then
   if (has_vxctau>0) then
    Paw_an(iat)%has_vxctau=1
    LIBPAW_ALLOCATE(Paw_an(iat)%vxctau1 ,(cplex*Paw_an(iat)%mesh_size,v_size,nspden))
    LIBPAW_ALLOCATE(Paw_an(iat)%vxcttau1,(cplex*Paw_an(iat)%mesh_size,v_size,nspden))
    Paw_an(iat)%vxctau1=zero;Paw_an(iat)%vxcttau1=zero
   end if
  end if

  ! ==========================
  ! === Optional arguments ===
  ! ==========================

  ! * XC potential inside PAW spheres generated by valence electrons
  Paw_an(iat)%has_vxcval=0
  if (PRESENT(has_vxcval)) then
   if (has_vxcval>0) then
    Paw_an(iat)%has_vxcval=1
    LIBPAW_ALLOCATE(Paw_an(iat)%vxc1_val ,(cplex*Paw_an(iat)%mesh_size,v_size,nspden))
    LIBPAW_ALLOCATE(Paw_an(iat)%vxct1_val,(cplex*Paw_an(iat)%mesh_size,v_size,nspden))
    Paw_an(iat)%vxc1_val=zero;Paw_an(iat)%vxct1_val=zero
   end if
  end if

  ! * Hartree potential LM-moments inside the sphere
  Paw_an(iat)%has_vhartree=0
  if (PRESENT(has_vhartree)) then
   if (has_vhartree>0) then
    Paw_an(iat)%has_vhartree=1
    LIBPAW_ALLOCATE(Paw_an(iat)%vh1,(cplex*Paw_an(iat)%mesh_size,v_size,nspden))
    LIBPAW_ALLOCATE(Paw_an(iat)%vht1,(cplex*Paw_an(iat)%mesh_size,v_size,nspden))
    Paw_an(iat)%vh1=zero;Paw_an(iat)%vht1=zero
   end if
  end if

  ! xc kernels inside the sphere
  Paw_an(iat)%has_kxc=0
  if (PRESENT(has_kxc)) then
   if (has_kxc>0) then
    Paw_an(iat)%has_kxc=1
    LIBPAW_ALLOCATE(Paw_an(iat)%kxc1 ,(cplex*Paw_an(iat)%mesh_size,v_size,nkxc1))
    LIBPAW_ALLOCATE(Paw_an(iat)%kxct1,(cplex*Paw_an(iat)%mesh_size,v_size,nkxc1))
    if (nkxc1>0) then
      Paw_an(iat)%kxc1=zero;Paw_an(iat)%kxct1=zero
    end if
   end if
  end if

  ! xc kernel derivatives inside the sphere
  Paw_an(iat)%has_k3xc=0
  if (PRESENT(has_k3xc)) then
   if (has_k3xc>0) then
    Paw_an(iat)%has_k3xc=1
    LIBPAW_ALLOCATE(Paw_an(iat)%k3xc1 ,(cplex*Paw_an(iat)%mesh_size,v_size,nk3xc1))
    LIBPAW_ALLOCATE(Paw_an(iat)%k3xct1,(cplex*Paw_an(iat)%mesh_size,v_size,nk3xc1))
    if (nk3xc1>0) then
      Paw_an(iat)%k3xc1=zero;Paw_an(iat)%k3xct1=zero
    end if
   end if
  end if

  ! local exact-exchange potential inside the sphere
  Paw_an(iat)%has_vxc_ex=0
  if (PRESENT(has_vxc_ex)) then
   if (has_vxc_ex>0.and.Pawtab(itypat)%useexexch/=0) then
    Paw_an(iat)%has_vxc_ex=1
    LIBPAW_ALLOCATE(Paw_an(iat)%vxc_ex,(cplex*Paw_an(iat)%mesh_size,v_size,nspden))
    Paw_an(iat)%vxc_ex=zero
   end if
  end if

 end do !iat

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

end subroutine paw_an_init
!!***

!----------------------------------------------------------------------

!!****f* m_paw_an/paw_an_free
!! NAME
!! paw_an_free
!!
!! FUNCTION
!!  Deallocate pointers and nullify flags in a paw_an structure
!!
!! SIDE EFFECTS
!!  Paw_an(:)<type(Paw_an_type)>=various arrays given on ANgular mesh or ANgular moments
!!
!!  All associated pointers in Paw_an(:) are deallocated
!!
!! PARENTS
!!      bethe_salpeter,dfpt_scfcv,m_paral_pert,m_paw_an,respfn,scfcv,screening
!!      sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine paw_an_free(Paw_an)

!Arguments ------------------------------------
!arrays
 type(Paw_an_type),intent(inout) :: Paw_an(:)

!Local variables-------------------------------
 integer :: iat,natom

! *************************************************************************

 !@Paw_an_type

 natom=SIZE(Paw_an);if (natom==0) return

 do iat=1,natom
  if (allocated(Paw_an(iat)%lmselect ))  then
    LIBPAW_DEALLOCATE(Paw_an(iat)%lmselect)
  end if
  if (allocated(Paw_an(iat)%vh1      ))  then
    LIBPAW_DEALLOCATE(Paw_an(iat)%vh1)
  end if
  if (allocated(Paw_an(iat)%vht1     ))  then
    LIBPAW_DEALLOCATE(Paw_an(iat)%vht1)
  end if
  if (allocated(Paw_an(iat)%vxc1     ))  then
    LIBPAW_DEALLOCATE(Paw_an(iat)%vxc1)
  end if
  if (allocated(Paw_an(iat)%vxctau1     ))  then
    LIBPAW_DEALLOCATE(Paw_an(iat)%vxctau1)
  end if
  if (allocated(Paw_an(iat)%vxc1_val ))  then
    LIBPAW_DEALLOCATE(Paw_an(iat)%vxc1_val)
  end if
  if (allocated(Paw_an(iat)%vxct1    ))  then
    LIBPAW_DEALLOCATE(Paw_an(iat)%vxct1)
  end if
  if (allocated(Paw_an(iat)%vxcttau1    ))  then
    LIBPAW_DEALLOCATE(Paw_an(iat)%vxcttau1)
  end if
  if (allocated(Paw_an(iat)%vxct1_val))  then
    LIBPAW_DEALLOCATE(Paw_an(iat)%vxct1_val)
  end if
  if (allocated(Paw_an(iat)%kxc1     ))  then
    LIBPAW_DEALLOCATE(Paw_an(iat)%kxc1)
  end if
  if (allocated(Paw_an(iat)%kxct1    ))  then
    LIBPAW_DEALLOCATE(Paw_an(iat)%kxct1)
  end if
  if (allocated(Paw_an(iat)%k3xc1     ))  then
    LIBPAW_DEALLOCATE(Paw_an(iat)%k3xc1)
  end if
  if (allocated(Paw_an(iat)%k3xct1    ))  then
    LIBPAW_DEALLOCATE(Paw_an(iat)%k3xct1)
  end if
  if (allocated(Paw_an(iat)%vxc_ex   ))  then
    LIBPAW_DEALLOCATE(Paw_an(iat)%vxc_ex)
  end if

  ! === Reset all has_* flags ===
  Paw_an(iat)%has_kxc     =0
  Paw_an(iat)%has_k3xc    =0
  Paw_an(iat)%has_vhartree=0
  Paw_an(iat)%has_vxc     =0
  Paw_an(iat)%has_vxctau  =0
  Paw_an(iat)%has_vxcval  =0
  Paw_an(iat)%has_vxc_ex  =0
 end do !iat

end subroutine paw_an_free
!!***

!----------------------------------------------------------------------

!!****f* m_paw_an/paw_an_nullify
!! NAME
!!  paw_an_nullify
!!
!! FUNCTION
!!  Nullify pointers and flags in a paw_an structure
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!  Paw_an(:)<type(paw_an_type)>=PAW arrays given on ANgular mesh or ANgular moments.
!!                               Nullified in output
!!
!! PARENTS
!!      bethe_salpeter,dfpt_scfcv,m_paw_an,paw_qpscgw,respfn,scfcv,screening
!!      sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine paw_an_nullify(Paw_an)

!Arguments ------------------------------------
 type(Paw_an_type),intent(inout) :: Paw_an(:)

!Local variables-------------------------------
 integer :: iat,natom

! *************************************************************************

 !@Paw_an_type
 ! MGPAW: This one could be removed/renamed,
 ! variables can be initialized in the datatype declaration
 ! Do we need to expose this in the public API?

 natom=SIZE(Paw_an(:));if (natom==0) return

 do iat=1,natom
  ! Set all has_* flags to zero.
  Paw_an(iat)%has_kxc      =0
  Paw_an(iat)%has_k3xc     =0
  Paw_an(iat)%has_vhartree =0
  Paw_an(iat)%has_vxc      =0
  Paw_an(iat)%has_vxctau   =0
  Paw_an(iat)%has_vxcval   =0
  Paw_an(iat)%has_vxc_ex   =0
 end do

end subroutine paw_an_nullify
!!***

!----------------------------------------------------------------------

!!****f* m_paw_an/paw_an_copy
!! NAME
!!  paw_an_copy
!!
!! FUNCTION
!!  Copy one paw_an datastructure into another
!!  Can take into accound changes of dimensions
!!  Can copy a shared paw_an into distributed ones (when parallelism is activated)
!!
!! INPUTS
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  paw_an_in(:)<type(paw_an_type)>= input paw_an datastructure
!!
!! SIDE EFFECTS
!!  paw_an_cpy(:)<type(paw_an_type)>= output paw_an datastructure
!!
!! NOTES
!!  paw_an_cpy must have been allocated in the calling function.
!!
!! PARENTS
!!      m_paw_an
!!
!! CHILDREN
!!
!! SOURCE

subroutine paw_an_copy(paw_an_in,paw_an_cpy,&
&                      mpi_atmtab,comm_atom)  ! optional arguments (parallelism)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: comm_atom
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 type(Paw_an_type),intent(in),target :: paw_an_in(:)
 type(Paw_an_type),intent(inout),target :: paw_an_cpy(:)

!Local variables-------------------------------
!scalars
 integer :: cplx_mesh_size,ij,ij1,lm_size,my_comm_atom,my_natom,nkxc1,nk3xc1,npaw_an_in
 integer :: npaw_an_max,npaw_an_out,nspden,paral_case,v_size,sz1
 logical :: my_atmtab_allocated,paral_atom
 character(len=500) :: msg
 type(Paw_an_type),pointer :: paw_an_in1, paw_an_out1
!arrays
 integer,pointer :: my_atmtab(:)
 type(Paw_an_type), pointer :: paw_an_out(:)

! *************************************************************************

!@Paw_an_type

!Retrieve sizes
 npaw_an_in=size(paw_an_in);npaw_an_out=size(paw_an_cpy)

!Set up parallelism over atoms
 paral_atom=(present(comm_atom));if (paral_atom) paral_atom=(xmpi_comm_size(comm_atom)>1)
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 my_atmtab_allocated=.false.

!Determine in which case we are (parallelism, ...)
!No parallelism: a single copy operation
 paral_case=0;npaw_an_max=npaw_an_in
 paw_an_out => paw_an_cpy
 if (paral_atom) then
   if (npaw_an_out<npaw_an_in) then ! Parallelism: the copy operation is a scatter
     call get_my_natom(my_comm_atom,my_natom,npaw_an_in)
     if (my_natom==npaw_an_out) then
       call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,npaw_an_in)
       paral_case=1;npaw_an_max=npaw_an_out
       paw_an_out => paw_an_cpy
     else
       msg=' npaw_an_out should be equal to my_natom !'
       MSG_BUG(msg)
     end if
   else                            ! Parallelism: the copy operation is a gather
     call get_my_natom(my_comm_atom,my_natom,npaw_an_out)
     if (my_natom==npaw_an_in) then
       paral_case=2;npaw_an_max=npaw_an_in
     else
       msg=' npaw_ij_in should be equal to my_natom !'
       MSG_BUG(msg)
     end if
   end if
 end if

!First case: a simple copy or a scatter
 if (npaw_an_max>0.and.((paral_case==0).or.(paral_case==1))) then
   call paw_an_free(paw_an_cpy)
   call paw_an_nullify(paw_an_cpy)

!  Loop on paw_ij components
   do ij1=1,npaw_an_max
     ij=ij1; if (paral_case==1) ij=my_atmtab(ij1)

     paw_an_in1=>paw_an_in(ij)
     paw_an_out1=>paw_an_out(ij1)
     paw_an_out1%angl_size =paw_an_in1%angl_size
     paw_an_out1%cplex =paw_an_in1%cplex
     paw_an_out1%has_kxc =paw_an_in1%has_kxc
     paw_an_out1%has_k3xc =paw_an_in1%has_k3xc
     paw_an_out1%has_vhartree =paw_an_in1%has_vhartree
     paw_an_out1%has_vxc =paw_an_in1%has_vxc
     paw_an_out1%has_vxctau =paw_an_in1%has_vxctau
     paw_an_out1%has_vxcval =paw_an_in1%has_vxcval
     paw_an_out1%has_vxc_ex =paw_an_in1%has_vxc_ex
     paw_an_out1%itypat =paw_an_in1%itypat
     paw_an_out1%lm_size =paw_an_in1%lm_size
     paw_an_out1%mesh_size =paw_an_in1%mesh_size
     paw_an_out1%nkxc1 =paw_an_in1%nkxc1
     paw_an_out1%nk3xc1 =paw_an_in1%nk3xc1
     paw_an_out1%nspden =paw_an_in1%nspden
     if (allocated(paw_an_in1%lmselect)) then
       sz1=size(paw_an_in1%lmselect)
       LIBPAW_ALLOCATE(paw_an_out1%lmselect,(sz1))
       paw_an_out1%lmselect(:)=paw_an_in1%lmselect(:)
     end if
     v_size=0
     if (paw_an_in1%has_vxc>0) then
       v_size=size(paw_an_in1%vxc1,2)
     else if (paw_an_in1%has_vxctau>0) then
       v_size=size(paw_an_in1%vxctau1,2)
     else if (paw_an_in1%has_kxc>0) then
       v_size=size(paw_an_in1%kxc1,2)
     else if (paw_an_in1%has_k3xc>0) then
       v_size=size(paw_an_in1%k3xc1,2)
     else if (paw_an_in1%has_vxcval>0) then
       v_size=size(paw_an_in1%vxc1_val,2)
     else if (paw_an_in1%has_vxc_ex>0) then
       v_size=size(paw_an_in1%vxc_ex,2)
     else if (paw_an_in1%has_vhartree>0) then
       v_size=size(paw_an_in1%vh1,2)
     end if
     nspden=paw_an_in1%nspden
     lm_size=paw_an_in1%lm_size
     cplx_mesh_size=paw_an_in1%cplex*paw_an_in1%mesh_size
     nkxc1=paw_an_in1%nkxc1
     if (paw_an_in1%has_kxc>0) then
       LIBPAW_ALLOCATE(paw_an_out1%kxc1,(cplx_mesh_size,v_size,nkxc1))
       LIBPAW_ALLOCATE(paw_an_out1%kxct1,(cplx_mesh_size,v_size,nkxc1))
       if (paw_an_in1%has_kxc==2.and.nkxc1>0) then
         paw_an_out1%kxc1(:,:,:)=paw_an_in1%kxc1(:,:,:)
         paw_an_out1%kxct1(:,:,:)=paw_an_in1%kxct1(:,:,:)
       end if
     end if
     nk3xc1=paw_an_in1%nk3xc1
     if (paw_an_in1%has_k3xc>0) then
       LIBPAW_ALLOCATE(paw_an_out1%k3xc1,(cplx_mesh_size,v_size,nk3xc1))
       LIBPAW_ALLOCATE(paw_an_out1%k3xct1,(cplx_mesh_size,v_size,nk3xc1))
       if (paw_an_in1%has_k3xc==2.and.nk3xc1>0) then
         paw_an_out1%k3xc1(:,:,:)=paw_an_in1%k3xc1(:,:,:)
         paw_an_out1%k3xct1(:,:,:)=paw_an_in1%k3xct1(:,:,:)
       end if
     end if
     if (paw_an_in1%has_vhartree>0) then
       LIBPAW_ALLOCATE(paw_an_out1%vh1,(cplx_mesh_size,lm_size,nspden))
       LIBPAW_ALLOCATE(paw_an_out1%vht1,(cplx_mesh_size,lm_size,nspden))
       if (paw_an_in1%has_vhartree==2) then
         paw_an_out1%vh1(:,:,:)=paw_an_in1%vh1(:,:,:)
         paw_an_out1%vht1(:,:,:)=paw_an_in1%vht1(:,:,:)
       end if
     end if
     if (paw_an_in1%has_vxc>0) then
       LIBPAW_ALLOCATE(paw_an_out1%vxc1,(cplx_mesh_size,v_size,nspden))
       LIBPAW_ALLOCATE(paw_an_out1%vxct1,(cplx_mesh_size,v_size,nspden))
       if (paw_an_in1%has_vxc==2) then
         paw_an_out1%vxc1(:,:,:)=paw_an_in1%vxc1(:,:,:)
         paw_an_out1%vxct1(:,:,:)=paw_an_in1%vxct1(:,:,:)
       end if
     end if
     if (paw_an_in1%has_vxctau>0) then
       LIBPAW_ALLOCATE(paw_an_out1%vxctau1,(cplx_mesh_size,v_size,nspden))
       LIBPAW_ALLOCATE(paw_an_out1%vxcttau1,(cplx_mesh_size,v_size,nspden))
       if (paw_an_in1%has_vxc==2) then
         paw_an_out1%vxctau1(:,:,:)=paw_an_in1%vxctau1(:,:,:)
         paw_an_out1%vxcttau1(:,:,:)=paw_an_in1%vxcttau1(:,:,:)
       end if
     end if
     if (paw_an_in1%has_vxcval>0) then
       LIBPAW_ALLOCATE(paw_an_out1%vxc1_val,(cplx_mesh_size,v_size,nspden))
       LIBPAW_ALLOCATE(paw_an_out1%vxct1_val,(cplx_mesh_size,v_size,nspden))
       if (paw_an_in1%has_vxcval==2) then
         paw_an_out1%vxc1_val(:,:,:)=paw_an_in1%vxc1_val(:,:,:)
         paw_an_out1%vxct1_val(:,:,:)=paw_an_in1%vxct1_val(:,:,:)
       end if
     end if
     if (paw_an_in1%has_vxc_ex>0) then
       LIBPAW_ALLOCATE(paw_an_out1%vxc_ex,(cplx_mesh_size,v_size,nspden))
       if (paw_an_in1%has_vxc_ex==2) then
         paw_an_out1%vxc_ex(:,:,:)=paw_an_in1%vxc_ex(:,:,:)
       end if
     end if
   end do
 end if

!Second case: a gather
 if (paral_case==2) then
   call paw_an_gather(paw_an_in,paw_an_cpy,-1,my_comm_atom,my_atmtab)
 end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 end subroutine paw_an_copy
!!***

!----------------------------------------------------------------------

!!****f* m_paw_an/paw_an_print
!! NAME
!!  paw_an_print
!!
!! FUNCTION
!!  Reports basic info on a paw_an datastructure
!!
!! INPUTS
!!  [unit]=the unit number for output
!!  [mode_paral]=either "COLL" or "PERS"
!!  [mpi_atmtab(:)]=indexes of the atoms treated by current proc (can be computed here)
!!  [comm_atom]=MPI communicator over atoms (needed if parallelism over atoms is activated)
!!  [natom]=total number of atom (needed if parallelism over atoms is activated)
!!          if Paw_an is distributed, natom is different from size(Paw_an).
!!
!! OUTPUT
!! (only writing)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine paw_an_print(Paw_an,unit,mode_paral, &
&                       mpi_atmtab,comm_atom,natom)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: comm_atom,natom,unit
 character(len=4),optional,intent(in) :: mode_paral
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 type(Paw_an_type),intent(in) :: Paw_an(:)

!Local variables-------------------------------
!scalars
 integer :: iatom,iatom_tot,my_comm_atom,my_natom,my_unt,size_paw_an
 logical :: my_atmtab_allocated,paral_atom
 character(len=4) :: my_mode
 character(len=500) :: msg
!arrays
 integer,pointer :: my_atmtab(:)

! *************************************************************************

!@Paw_an_type

 size_paw_an=SIZE(Paw_an)
 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_mode  ='PERS' ; if (PRESENT(mode_paral)) my_mode  =mode_paral
 my_natom=size_paw_an; if (PRESENT(natom))      my_natom=natom

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.my_natom/=size_paw_an)
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,my_natom,my_natom_ref=size_paw_an)

 write(msg,'(3a)')ch10,' === Content of the pawfgrtab datatype === ',ch10
 call wrtout(my_unt,msg,my_mode)

 do iatom=1,my_natom
   iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)
   write(msg,'(a)')'                                 '
   call wrtout(my_unt,msg,my_mode)
   write(msg,'(a,i4)')'  ****************************** iatom= ' , iatom_tot
   call wrtout(my_unt,msg,my_mode)
   write(msg,'(a,i4)')'  Dimension of paw angular mesh= ',paw_an(iatom)%angl_size
   call wrtout(my_unt,msg,my_mode)
   write(msg,'(a,i4)')'  cplex (1 if potentials/densities are real, 2 if they are complex)= ',&
&        paw_an(iatom)%cplex
   call wrtout(my_unt,msg,my_mode)
   write(msg,'(a,i4)')'  has_kxc     = ',paw_an(iatom)%has_kxc
   call wrtout(my_unt,msg,my_mode)
   write(msg,'(a,i4)')'  has_k3xc    = ',paw_an(iatom)%has_k3xc
   call wrtout(my_unt,msg,my_mode)
   write(msg,'(a,i4)')'  has_vhartree= ',paw_an(iatom)%has_vhartree
   call wrtout(my_unt,msg,my_mode)
   write(msg,'(a,i4)')'  has_vxc     = ',paw_an(iatom)%has_vxc
   call wrtout(my_unt,msg,my_mode)
   write(msg,'(a,i4)')'  has_vxctau  = ',paw_an(iatom)%has_vxctau
   call wrtout(my_unt,msg,my_mode)
   write(msg,'(a,i4)')'  has_vxcval  = ',paw_an(iatom)%has_vxcval
   call wrtout(my_unt,msg,my_mode)
   write(msg,'(a,i4)')'  has_vxc_ex  = ',paw_an(iatom)%has_vxc_ex
   call wrtout(my_unt,msg,my_mode)
   write(msg,'(a,i4)')'  Atome type  = ',paw_an(iatom)%itypat
   call wrtout(my_unt,msg,my_mode)
   write(msg,'(a,i4)')'  lm_size     = ',paw_an(iatom)%lm_size
   call wrtout(my_unt,msg,my_mode)
   write(msg,'(a,i4)')'  mesh_size   = ',paw_an(iatom)%mesh_size
   call wrtout(my_unt,msg,my_mode)
   write(msg,'(a,i4)')'  nkxc1       = ',paw_an(iatom)%nkxc1
   call wrtout(my_unt,msg,my_mode)
   write(msg,'(a,i4)')'  nk3xc1      = ',paw_an(iatom)%nk3xc1
   call wrtout(my_unt,msg,my_mode)
   write(msg,'(a,i4)')'  nspden      = ',paw_an(iatom)%nspden
   call wrtout(my_unt,msg,my_mode)
 end do

end subroutine paw_an_print
!!***

!----------------------------------------------------------------------

!!****f* m_paw_an/paw_an_gather
!! NAME
!!  paw_an_gather
!!
!! FUNCTION
!!  (All)Gather paw_an datastructures
!!
!! INPUTS
!!  master=master communicator receiving data ; if -1 do a ALLGATHER
!!  comm_atom= communicator over atom
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current calling proc
!!  paw_an_in(:)<type(paw_an_type)>= input paw_an datastructures on every process
!!
!! OUTPUT
!!  paw_an_gathered(:)<type(paw_an_type)>= output paw_an datastructure
!!
!! PARENTS
!!      m_paw_an
!!
!! CHILDREN
!!
!! SOURCE

subroutine paw_an_gather(Paw_an_in,paw_an_gathered,master,comm_atom,mpi_atmtab)

!Arguments ------------------------------------
 integer,intent(in) :: master,comm_atom
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 type(Paw_an_type),target,intent(in) :: Paw_an_in(:)
 type(Paw_an_type),target,intent(inout) :: Paw_an_gathered(:)

!Local variables-------------------------------
!scalars
 integer :: buf_dp_size,buf_dp_size_all,buf_int_size,buf_int_size_all,cplx_mesh_size
 integer :: iat,iatot,ierr,has_lm_select,i1,i2,ij,indx_int,indx_dp
 integer :: lm_size,me_atom
 integer :: my_natom,natom,nkxc1,nk3xc1,npaw_an_in_sum,nproc_atom,nspden,v_size,sz1,sz2,sz3
 logical :: my_atmtab_allocated,paral_atom
 character(len=500) :: msg
 type(Paw_an_type),pointer :: paw_an_in1,paw_an_gathered1
!arrays
 integer :: bufsz(2)
 integer,allocatable :: buf_int(:),buf_int_all(:)
 integer,allocatable :: count_dp(:),count_int(:),count_tot(:),displ_dp(:),displ_int(:)
 integer,pointer :: my_atmtab(:)
 real(dp),allocatable :: buf_dp(:),buf_dp_all(:)

! *************************************************************************

!@Paw_an_type

 if (master/=-1) then
   msg='simple gather (master/=-1) not yet implemented !'
   MSG_BUG(msg)
 end if

 my_natom=size(paw_an_in);natom=size(paw_an_gathered)

!Set up parallelism over atoms
 paral_atom=(my_natom/=natom)
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 call get_my_atmtab(comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)
 nproc_atom=xmpi_comm_size(comm_atom)
 me_atom=xmpi_comm_rank(comm_atom)

!Special case: one process (simple copy)
 if (nproc_atom==1) then
   if (master==-1.or.me_atom==master) then
     call paw_an_free(paw_an_gathered)
     call paw_an_nullify(paw_an_gathered)
     do iat=1,my_natom
       paw_an_in1=>paw_an_in(iat)
       paw_an_gathered1%itypat =paw_an_in1%itypat
       paw_an_gathered1%nspden =paw_an_in1%nspden
       paw_an_gathered1%cplex =paw_an_in1%cplex
       paw_an_gathered1%mesh_size =paw_an_in1%mesh_size
       paw_an_gathered1%angl_size =paw_an_in1%angl_size
       paw_an_gathered1%lm_size =paw_an_in1%lm_size
       paw_an_gathered1%nkxc1 =paw_an_in1%nkxc1
       paw_an_gathered1%nk3xc1 =paw_an_in1%nk3xc1
       paw_an_gathered1%has_vxc =paw_an_in1%has_vxc
       paw_an_gathered1%has_vxctau =paw_an_in1%has_vxctau
       paw_an_gathered1%has_kxc =paw_an_in1%has_kxc
       paw_an_gathered1%has_k3xc =paw_an_in1%has_k3xc
       paw_an_gathered1%has_vxcval =paw_an_in1%has_vxcval
       paw_an_gathered1%has_vxc_ex =paw_an_in1%has_vxc_ex
       paw_an_gathered1%has_vhartree =paw_an_in1%has_vhartree
       if (allocated(paw_an_in1%lmselect)) then
         sz1=size(paw_an_in1%lmselect)
         LIBPAW_ALLOCATE(paw_an_gathered1%lmselect,(sz1))
         paw_an_gathered1%lmselect(:)=paw_an_in1%lmselect(:)
       end if
       if (allocated(paw_an_in1%vxc1)) then
         sz1=size(paw_an_in1%vxc1,1);sz2=size(paw_an_in1%vxc1,2)
         sz3=size(paw_an_in1%vxc1,3)
         LIBPAW_ALLOCATE(paw_an_gathered1%vxc1,(sz1,sz2,sz3))
         paw_an_gathered1%vxc1(:,:,:)=paw_an_in1%vxc1(:,:,:)
       end if
       if (allocated(paw_an_in1%vxctau1)) then
         sz1=size(paw_an_in1%vxctau1,1);sz2=size(paw_an_in1%vxctau1,2)
         sz3=size(paw_an_in1%vxctau1,3)
         LIBPAW_ALLOCATE(paw_an_gathered1%vxctau1,(sz1,sz2,sz3))
         paw_an_gathered1%vxctau1(:,:,:)=paw_an_in1%vxctau1(:,:,:)
       end if
       if (allocated(paw_an_in1%vxcttau1)) then
         sz1=size(paw_an_in1%vxcttau1,1);sz2=size(paw_an_in1%vxcttau1,2)
         sz3=size(paw_an_in1%vxcttau1,3)
         LIBPAW_ALLOCATE(paw_an_gathered1%vxcttau1,(sz1,sz2,sz3))
         paw_an_gathered1%vxcttau1(:,:,:)=paw_an_in1%vxcttau1(:,:,:)
       end if
       if (allocated(paw_an_in1%kxc1)) then
         sz1=size(paw_an_in1%kxc1,1);sz2=size(paw_an_in1%kxc1,2)
         sz3=size(paw_an_in1%kxc1,3)
         LIBPAW_ALLOCATE(paw_an_gathered1%kxc1,(sz1,sz2,sz3))
         if (sz3>0) paw_an_gathered1%kxc1(:,:,:)=paw_an_in1%kxc1(:,:,:)
       end if
       if (allocated(paw_an_in1%k3xc1)) then
         sz1=size(paw_an_in1%k3xc1,1);sz2=size(paw_an_in1%k3xc1,2)
         sz3=size(paw_an_in1%k3xc1,3)
         LIBPAW_ALLOCATE(paw_an_gathered1%k3xc1,(sz1,sz2,sz3))
         if (sz3>0) paw_an_gathered1%k3xc1(:,:,:)=paw_an_in1%k3xc1(:,:,:)
       end if
       if (allocated(paw_an_in1%kxct1)) then
         sz1=size(paw_an_in1%kxct1,1);sz2=size(paw_an_in1%kxct1,2)
         sz3=size(paw_an_in1%kxct1,3)
         LIBPAW_ALLOCATE(paw_an_gathered1%kxct1,(sz1,sz2,sz3))
         if (sz3>0) paw_an_gathered1%kxct1(:,:,:)=paw_an_in1%kxct1(:,:,:)
       end if
       if (allocated(paw_an_in1%k3xct1)) then
         sz1=size(paw_an_in1%k3xct1,1);sz2=size(paw_an_in1%k3xct1,2)
         sz3=size(paw_an_in1%k3xct1,3)
         LIBPAW_ALLOCATE(paw_an_gathered1%k3xct1,(sz1,sz2,sz3))
         if (sz3>0) paw_an_gathered1%k3xct1(:,:,:)=paw_an_in1%k3xct1(:,:,:)
       end if
       if (allocated(paw_an_in1%vxc1_val)) then
         sz1=size(paw_an_in1%vxc1_val,1);sz2=size(paw_an_in1%vxc1_val,2)
         sz3=size(paw_an_in1%vxc1_val,3)
         LIBPAW_ALLOCATE(paw_an_gathered1%vxc1_val,(sz1,sz2,sz3))
         paw_an_gathered1%vxc1_val(:,:,:)=paw_an_in1%vxc1_val(:,:,:)
       end if
       if (allocated(paw_an_in1%vxct1_val)) then
         sz1=size(paw_an_in1%vxct1_val,1);sz2=size(paw_an_in1%vxct1_val,2)
         sz3=size(paw_an_in1%vxct1_val,3)
         LIBPAW_ALLOCATE(paw_an_gathered1%vxct1_val,(sz1,sz2,sz3))
         paw_an_gathered1%vxct1_val(:,:,:)=paw_an_in1%vxct1_val(:,:,:)
       end if
       if (allocated(paw_an_in1%vxc_ex)) then
         sz1=size(paw_an_in1%vxc_ex,1);sz2=size(paw_an_in1%vxc_ex,2)
         sz3=size(paw_an_in1%vxc_ex,3)
         LIBPAW_ALLOCATE(paw_an_gathered1%vxc_ex,(sz1,sz2,sz3))
         paw_an_gathered1%vxc_ex(:,:,:)=paw_an_in1%vxc_ex(:,:,:)
       end if
       if (allocated(paw_an_in1%vh1)) then
         sz1=size(paw_an_in1%vh1,1);sz2=size(paw_an_in1%vh1,2)
         sz3=size(paw_an_in1%vh1,3)
         LIBPAW_ALLOCATE(paw_an_gathered1%vh1,(sz1,sz2,sz3))
         paw_an_gathered1%vh1(:,:,:)=paw_an_in1%vh1(:,:,:)
       end if
       if (allocated(paw_an_in1%vht1)) then
         sz1=size(paw_an_in1%vht1,1);sz2=size(paw_an_in1%vht1,2)
         sz3=size(paw_an_in1%vht1,3)
         LIBPAW_ALLOCATE(paw_an_gathered1%vht1,(sz1,sz2,sz3))
         paw_an_gathered1%vht1(:,:,:)=paw_an_in1%vht1(:,:,:)
       end if
     end do
   end if
  return
 end if

!Test on sizes
 npaw_an_in_sum=my_natom
 call xmpi_sum(npaw_an_in_sum,comm_atom,ierr)
 if (master==-1) then
   if (natom/=npaw_an_in_sum) then
     msg='Wrong sizes sum[npaw_an_in]/=natom !'
     MSG_BUG(msg)
   end if
 else
   if (me_atom==master.and.natom/=npaw_an_in_sum) then
     msg='(2) paw_an_gathered wrongly allocated !'
     MSG_BUG(msg)
   end if
 end if

!Compute sizes of buffers
 buf_int_size=0;buf_dp_size=0
 do ij=1,my_natom
   buf_int_size=buf_int_size+17+size(paw_an_in(ij)%lmselect)
 end do
 do ij=1,my_natom
   paw_an_in1=>paw_an_in(ij)
   if (paw_an_in1%has_vxc==2) then
     buf_dp_size=buf_dp_size+size(paw_an_in1%vxc1)
     buf_dp_size=buf_dp_size+size(paw_an_in1%vxct1)
   end if
   if (paw_an_in1%has_vxctau==2) then
     buf_dp_size=buf_dp_size+size(paw_an_in1%vxctau1)
     buf_dp_size=buf_dp_size+size(paw_an_in1%vxcttau1)
   end if
   if (paw_an_in1%has_kxc==2) then
     buf_dp_size=buf_dp_size+size(paw_an_in1%kxc1)
     buf_dp_size=buf_dp_size+size(paw_an_in1%kxct1)
   end if
   if (paw_an_in1%has_k3xc==2) then
     buf_dp_size=buf_dp_size+size(paw_an_in1%k3xc1)
     buf_dp_size=buf_dp_size+size(paw_an_in1%k3xct1)
   end if
   if (paw_an_in1%has_vxcval==2) then
     buf_dp_size=buf_dp_size+size(paw_an_in1%vxc1_val)
     buf_dp_size=buf_dp_size+size(paw_an_in1%vxct1_val)
   end if
   if (paw_an_in1%has_vxc_ex==2) then
     buf_dp_size=buf_dp_size+size(paw_an_in1%vxc_ex)
   end if
   if (paw_an_in1%has_vhartree==2) then
     buf_dp_size=buf_dp_size+size(paw_an_in1%vh1)
     buf_dp_size=buf_dp_size+size(paw_an_in1%vht1)
   end if
 end do

!Fill in input buffers
 LIBPAW_ALLOCATE(buf_int,(buf_int_size))
 LIBPAW_ALLOCATE(buf_dp ,(buf_dp_size))
 indx_int=1;indx_dp=1
 do ij=1, my_natom
   paw_an_in1=>paw_an_in(ij)
   buf_int(indx_int)=my_atmtab(ij); indx_int=indx_int+1
   buf_int(indx_int)=paw_an_in1%itypat; indx_int=indx_int+1
   buf_int(indx_int)=paw_an_in1%nspden; indx_int=indx_int+1
   buf_int(indx_int)=paw_an_in1%cplex; indx_int=indx_int+1
   buf_int(indx_int)=paw_an_in1%mesh_size; indx_int=indx_int+1
   buf_int(indx_int)=paw_an_in1%angl_size; indx_int=indx_int+1
   buf_int(indx_int)=paw_an_in1%lm_size; indx_int=indx_int+1
   buf_int(indx_int)=paw_an_in1%nkxc1; indx_int=indx_int+1
   buf_int(indx_int)=paw_an_in1%nk3xc1; indx_int=indx_int+1
   buf_int(indx_int)=paw_an_in1%has_vxc; indx_int=indx_int+1
   buf_int(indx_int)=paw_an_in1%has_vxctau; indx_int=indx_int+1
   buf_int(indx_int)=paw_an_in1%has_kxc; indx_int=indx_int+1
   buf_int(indx_int)=paw_an_in1%has_k3xc; indx_int=indx_int+1
   buf_int(indx_int)=paw_an_in1%has_vxcval; indx_int=indx_int+1
   buf_int(indx_int)=paw_an_in1%has_vxc_ex; indx_int=indx_int+1
   buf_int(indx_int)=paw_an_in1%has_vhartree; indx_int=indx_int+1
   v_size=0
   if (paw_an_in1%has_vxc>0) then
     v_size=size(paw_an_in1%vxc1,2)
   else if (paw_an_in1%has_vxctau>0) then
     v_size=size(paw_an_in1%vxctau1,2)
   else if (paw_an_in1%has_kxc>0) then
     v_size=size(paw_an_in1%kxc1,2)
   else if (paw_an_in1%has_k3xc>0) then
     v_size=size(paw_an_in1%k3xc1,2)
   else if (paw_an_in1%has_vxcval>0) then
     v_size=size(paw_an_in1%vxc1_val,2)
   else if (paw_an_in1%has_vxc_ex>0) then
     v_size=size(paw_an_in1%vxc_ex,2)
   else if (paw_an_in1%has_vhartree>0) then
     v_size=size(paw_an_in1%vh1,2)
   end if
   buf_int(indx_int)=v_size;indx_int=indx_int+1
   if (allocated(paw_an_in1%lmselect)) then
     buf_int(indx_int)=1;indx_int=indx_int+1
   else
     buf_int(indx_int)=0;indx_int=indx_int+1
   end if
   nspden=paw_an_in1%nspden
   lm_size=paw_an_in1%lm_size
   cplx_mesh_size=paw_an_in1%cplex*paw_an_in1%mesh_size
   if (lm_size>0) then
     if (allocated(paw_an_in1%lmselect)) then
       do i1=1,lm_size
         if (paw_an_in1%lmselect(i1)) then
           buf_int(indx_int)=1
         else
           buf_int(indx_int)=0
         end if
         indx_int=indx_int+1
       end do
     end if
   end if
   if (paw_an_in1%has_vxc==2) then
     do i1=1,nspden
       do i2=1,v_size
         buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)=paw_an_in1%vxc1(:,i2,i1)
         indx_dp=indx_dp+cplx_mesh_size
       end do
     end do
     do i1=1,nspden
       do i2=1,v_size
         buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)=paw_an_in1%vxct1(:,i2,i1)
         indx_dp=indx_dp+cplx_mesh_size
       end do
     end do
   end if
   if (paw_an_in1%has_vxctau==2) then
     do i1=1,nspden
       do i2=1,v_size
         buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)=paw_an_in1%vxctau1(:,i2,i1)
         indx_dp=indx_dp+cplx_mesh_size
       end do
     end do
     do i1=1,nspden
       do i2=1,v_size
         buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)=paw_an_in1%vxcttau1(:,i2,i1)
         indx_dp=indx_dp+cplx_mesh_size
       end do
     end do
   end if
   if (paw_an_in1%has_kxc==2.and.paw_an_in1%nkxc1>0) then
     do i1=1,paw_an_in1%nkxc1
       do i2=1,v_size
         buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)=paw_an_in1%kxc1(:,i2,i1)
         indx_dp=indx_dp+cplx_mesh_size
       end do
     end do
     do i1=1,paw_an_in1%nkxc1
       do i2=1,v_size
         buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)=paw_an_in1%kxct1(:,i2,i1)
         indx_dp=indx_dp+cplx_mesh_size
       end do
     end do
   end if
   if (paw_an_in1%has_k3xc==2.and.paw_an_in1%nk3xc1>0) then
     do i1=1,paw_an_in1%nk3xc1
       do i2=1,v_size
         buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)=paw_an_in1%k3xc1(:,i2,i1)
         indx_dp=indx_dp+cplx_mesh_size
       end do
     end do
     do i1=1,paw_an_in1%nk3xc1
       do i2=1,v_size
         buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)=paw_an_in1%k3xct1(:,i2,i1)
         indx_dp=indx_dp+cplx_mesh_size
       end do
     end do
   end if
   if (paw_an_in1%has_vxcval==2) then
     do i1=1,nspden
       do i2=1,v_size
         buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)=paw_an_in1%vxc1_val(:,i2,i1)
         indx_dp=indx_dp+cplx_mesh_size
       end do
     end do
     do i1=1,nspden
       do i2=1,v_size
         buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)=paw_an_in1%vxct1_val(:,i2,i1)
         indx_dp=indx_dp+cplx_mesh_size
       end do
     end do
   end if
   if (paw_an_in1%has_vxc_ex==2) then
     do i1=1,nspden
       do i2=1,v_size
         buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)=paw_an_in1%vxc_ex(:,i2,i1)
         indx_dp=indx_dp+cplx_mesh_size
       end do
     end do
   end if
   if (paw_an_in1%has_vhartree==2) then
     do i1=1,nspden
       do i2=1,lm_size
         buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)=paw_an_in1%vh1(:,i2,i1)
         indx_dp=indx_dp+cplx_mesh_size
       end do
     end do
     do i1=1,nspden
       do i2=1,lm_size
         buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)=paw_an_in1%vht1(:,i2,i1)
         indx_dp=indx_dp+cplx_mesh_size
       end do
     end do
   end if
 end do
 if (indx_int/=1+buf_int_size) then
   msg='Error (1) in paw_an_gather: wrong buffer sizes !'
   MSG_BUG(msg)
 end if
 if (indx_dp/=1+buf_dp_size) then
   msg='Error (2) in paw_an_gather: wrong buffer sizes !'
   MSG_BUG(msg)
 end if

!Communicate (1 gather for integers, 1 gather for reals)
 LIBPAW_ALLOCATE(count_int,(nproc_atom))
 LIBPAW_ALLOCATE(displ_int,(nproc_atom))
 LIBPAW_ALLOCATE(count_dp ,(nproc_atom))
 LIBPAW_ALLOCATE(displ_dp ,(nproc_atom))
 LIBPAW_ALLOCATE(count_tot,(2*nproc_atom))
 bufsz(1)=buf_int_size; bufsz(2)=buf_dp_size
 call xmpi_allgather(bufsz,2,count_tot,comm_atom,ierr)
 do ij=1,nproc_atom
   count_int(ij)=count_tot(2*ij-1)
   count_dp (ij)=count_tot(2*ij)
 end do
 displ_int(1)=0;displ_dp(1)=0
 do ij=2,nproc_atom
   displ_int(ij)=displ_int(ij-1)+count_int(ij-1)
   displ_dp (ij)=displ_dp (ij-1)+count_dp (ij-1)
 end do
 buf_int_size_all=sum(count_int)
 buf_dp_size_all =sum(count_dp)
 LIBPAW_DEALLOCATE(count_tot)
 LIBPAW_ALLOCATE(buf_int_all,(buf_int_size_all))
 LIBPAW_ALLOCATE(buf_dp_all ,(buf_dp_size_all))
 call xmpi_allgatherv(buf_int,buf_int_size,buf_int_all,count_int,displ_int,comm_atom,ierr)
 call xmpi_allgatherv(buf_dp ,buf_dp_size ,buf_dp_all ,count_dp ,displ_dp ,comm_atom,ierr)
 LIBPAW_DEALLOCATE(count_int)
 LIBPAW_DEALLOCATE(displ_int)
 LIBPAW_DEALLOCATE(count_dp)
 LIBPAW_DEALLOCATE(displ_dp)

!Fill in output datastructure
 indx_int=1; indx_dp=1
 call paw_an_free(paw_an_gathered)
 call paw_an_nullify(paw_an_gathered)
 do iat=1,natom
   iatot=buf_int_all(indx_int); indx_int=indx_int+1
   paw_an_gathered1=>paw_an_gathered(iatot)
   paw_an_gathered1%itypat=buf_int_all(indx_int); indx_int=indx_int+1
   paw_an_gathered1%nspden=buf_int_all(indx_int); indx_int=indx_int+1
   paw_an_gathered1%cplex=buf_int_all(indx_int); indx_int=indx_int+1
   paw_an_gathered1%mesh_size=buf_int_all(indx_int); indx_int=indx_int+1
   paw_an_gathered1%angl_size=buf_int_all(indx_int); indx_int=indx_int+1
   paw_an_gathered1%lm_size=buf_int_all(indx_int); indx_int=indx_int+1
   paw_an_gathered1%nkxc1=buf_int_all(indx_int); indx_int=indx_int+1
   paw_an_gathered1%nk3xc1=buf_int_all(indx_int); indx_int=indx_int+1
   paw_an_gathered1%has_vxc=buf_int_all(indx_int); indx_int=indx_int+1
   paw_an_gathered1%has_vxctau=buf_int_all(indx_int); indx_int=indx_int+1
   paw_an_gathered1%has_kxc=buf_int_all(indx_int); indx_int=indx_int+1
   paw_an_gathered1%has_k3xc=buf_int_all(indx_int); indx_int=indx_int+1
   paw_an_gathered1%has_vxcval=buf_int_all(indx_int); indx_int=indx_int+1
   paw_an_gathered1%has_vxc_ex=buf_int_all(indx_int); indx_int=indx_int+1
   paw_an_gathered1%has_vhartree=buf_int_all(indx_int); indx_int=indx_int+1
   v_size=buf_int_all(indx_int); indx_int=indx_int+1
   has_lm_select=buf_int_all(indx_int); indx_int=indx_int+1
   nspden=paw_an_gathered1%nspden
   lm_size=paw_an_gathered1%lm_size
   nkxc1=paw_an_gathered1%nkxc1
   nk3xc1=paw_an_gathered1%nk3xc1
   cplx_mesh_size=paw_an_gathered1%cplex*paw_an_gathered1%mesh_size
   if (has_lm_select==1) then
     LIBPAW_ALLOCATE(paw_an_gathered1%lmselect,(lm_size))
     if (lm_size>0) then
       do i1=1,lm_size
         if (buf_int_all(indx_int)==1) then
           paw_an_gathered1%lmselect(i1)=.TRUE.;indx_int=indx_int+1
         else
           paw_an_gathered1%lmselect(i1)=.FALSE.;indx_int=indx_int+1
         end if
       end do
     end if
   end if
   if (paw_an_gathered1%has_vxc>0) then
     LIBPAW_ALLOCATE(paw_an_gathered1%vxc1,(cplx_mesh_size,v_size,nspden))
     LIBPAW_ALLOCATE(paw_an_gathered1%vxct1,(cplx_mesh_size,v_size,nspden))
     if (paw_an_gathered1%has_vxc==2) then
       do i1=1,nspden
         do i2=1,v_size
           paw_an_gathered1%vxc1(:,i2,i1)=buf_dp_all(indx_dp:indx_dp+cplx_mesh_size-1)
           indx_dp=indx_dp+cplx_mesh_size
         end do
       end do
       do i1=1,nspden
         do i2=1,v_size
           paw_an_gathered1%vxct1(:,i2,i1)=buf_dp_all(indx_dp:indx_dp+cplx_mesh_size-1)
           indx_dp=indx_dp+cplx_mesh_size
         end do
       end do
     end if
   end if
   if (paw_an_gathered1%has_vxctau>0) then
     LIBPAW_ALLOCATE(paw_an_gathered1%vxctau1,(cplx_mesh_size,v_size,nspden))
     LIBPAW_ALLOCATE(paw_an_gathered1%vxcttau1,(cplx_mesh_size,v_size,nspden))
     if (paw_an_gathered1%has_vxctau==2) then
       do i1=1,nspden
         do i2=1,v_size
           paw_an_gathered1%vxctau1(:,i2,i1)=buf_dp_all(indx_dp:indx_dp+cplx_mesh_size-1)
           indx_dp=indx_dp+cplx_mesh_size
         end do
       end do
       do i1=1,nspden
         do i2=1,v_size
           paw_an_gathered1%vxcttau1(:,i2,i1)=buf_dp_all(indx_dp:indx_dp+cplx_mesh_size-1)
           indx_dp=indx_dp+cplx_mesh_size
         end do
       end do
     end if
   end if
   if (paw_an_gathered1%has_kxc>0) then
     LIBPAW_ALLOCATE(paw_an_gathered1%kxc1,(cplx_mesh_size,v_size,nkxc1))
     LIBPAW_ALLOCATE(paw_an_gathered1%kxct1,(cplx_mesh_size,v_size,nkxc1))
     if (paw_an_gathered1%has_kxc==2.and.nkxc1>0) then
       do i1=1,nkxc1
         do i2=1,v_size
           paw_an_gathered1%kxc1(:,i2,i1)=buf_dp_all(indx_dp:indx_dp+cplx_mesh_size-1)
           indx_dp=indx_dp+cplx_mesh_size
         end do
       end do
       do i1=1,nkxc1
         do i2=1,v_size
           paw_an_gathered1%kxct1(:,i2,i1)=buf_dp_all(indx_dp:indx_dp+cplx_mesh_size-1)
           indx_dp=indx_dp+cplx_mesh_size
         end do
       end do
     end if
    end if
   if (paw_an_gathered1%has_k3xc>0) then
     LIBPAW_ALLOCATE(paw_an_gathered1%k3xc1,(cplx_mesh_size,v_size,nk3xc1))
     LIBPAW_ALLOCATE(paw_an_gathered1%k3xct1,(cplx_mesh_size,v_size,nk3xc1))
     if (paw_an_gathered1%has_k3xc==2.and.nk3xc1>0) then
       do i1=1,nk3xc1
         do i2=1,v_size
           paw_an_gathered1%k3xc1(:,i2,i1)=buf_dp_all(indx_dp:indx_dp+cplx_mesh_size-1)
           indx_dp=indx_dp+cplx_mesh_size
         end do
       end do
       do i1=1,nk3xc1
         do i2=1,v_size
           paw_an_gathered1%k3xct1(:,i2,i1)=buf_dp_all(indx_dp:indx_dp+cplx_mesh_size-1)
           indx_dp=indx_dp+cplx_mesh_size
         end do
       end do
     end if
    end if
   if (paw_an_gathered1%has_vxcval>0) then
     LIBPAW_ALLOCATE(paw_an_gathered1%vxc1_val,(cplx_mesh_size,v_size,nspden))
     LIBPAW_ALLOCATE(paw_an_gathered1%vxct1_val,(cplx_mesh_size,v_size,nspden))
     if (paw_an_gathered1%has_vxcval==2) then
       do i1=1,nspden
         do i2=1,v_size
           paw_an_gathered1%vxc1_val(:,i2,i1)=buf_dp_all(indx_dp:indx_dp+cplx_mesh_size-1)
           indx_dp=indx_dp+cplx_mesh_size
         end do
       end do
       do i1=1,nspden
         do i2=1,v_size
           paw_an_gathered1%vxct1_val(:,i2,i1)=buf_dp_all(indx_dp:indx_dp+cplx_mesh_size-1)
           indx_dp=indx_dp+cplx_mesh_size
         end do
       end do
     end if
   end if
   if (paw_an_gathered1%has_vxc_ex>0) then
     LIBPAW_ALLOCATE(paw_an_gathered1%vxc_ex,(cplx_mesh_size,v_size,nspden))
     if (paw_an_gathered1%has_vxc_ex==2) then
       do i1=1,nspden
         do i2=1,v_size
           paw_an_gathered1%vxc_ex(:,i2,i1)=buf_dp_all(indx_dp:indx_dp+cplx_mesh_size-1)
           indx_dp=indx_dp+cplx_mesh_size
         end do
       end do
     end if
   end if
   if (paw_an_gathered1%has_vhartree>0) then
     LIBPAW_ALLOCATE(paw_an_gathered1%vh1,(cplx_mesh_size,lm_size,nspden))
     LIBPAW_ALLOCATE(paw_an_gathered1%vht1,(cplx_mesh_size,lm_size,nspden))
     if (paw_an_gathered1%has_vhartree==2) then
       do i1=1,nspden
         do i2=1,lm_size
           paw_an_gathered1%vh1(:,i2,i1)=buf_dp_all(indx_dp:indx_dp+cplx_mesh_size-1)
           indx_dp=indx_dp+cplx_mesh_size
         end do
       end do
       do i1=1,nspden
         do i2=1,lm_size
           paw_an_gathered1%vht1(:,i2,i1)=buf_dp_all(indx_dp:indx_dp+cplx_mesh_size-1)
           indx_dp=indx_dp+cplx_mesh_size
         end do
       end do
     end if
   end if
 end do ! iat

!Free buffers
 LIBPAW_DEALLOCATE(buf_int)
 LIBPAW_DEALLOCATE(buf_int_all)
 LIBPAW_DEALLOCATE(buf_dp)
 LIBPAW_DEALLOCATE(buf_dp_all)

!Destroy atom table
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

end subroutine paw_an_gather
!!***

!----------------------------------------------------------------------

!!****f* m_paw_an/paw_an_redistribute
!! NAME
!! paw_an_redistribute
!!
!! FUNCTION
!!   Redistribute an array of paw_an datastructures
!!   Input paw_an is given on a MPI communicator
!!   Output paw_an is redistributed on another MPI communicator
!!
!! INPUTS
!!  mpi_comm_in= input MPI (atom) communicator
!!  mpi_comm_out= output MPI (atom) communicator
!!  mpi_atmtab_in= --optional-- indexes of the input paw_an treated by current proc
!!                 if not present, will be calculated in the present routine
!!  mpi_atmtab_out= --optional-- indexes of the output paw_an treated by current proc
!!                  if not present, will be calculated in the present routine
!!  natom= --optional-- total number of atoms
!!  ----- Optional arguments used only for asynchronous communications -----
!!    RecvAtomProc(:)= rank of processor from which I expect atom (in mpi_comm_in)
!!    RecvAtomList(:)= indexes of atoms to be received by me
!!      RecvAtomList(irecv) are the atoms I expect from RecvAtomProc(irecv)
!!    SendAtomProc(:)= ranks of process destination of atom (in mpi_comm_in)
!!    SendAtomList(:)= indexes of atoms to be sent by me
!!      SendAtomList(isend) are the atoms sent to SendAtomProc(isend)
!!
!! OUTPUT
!!  [paw_an_out(:)]<type(paw_an_type)>= --optional--
!!                    if present, the redistributed datastructure does not replace
!!                    the input one but is delivered in paw_an_out
!!                    if not present, input and output datastructure are the same.
!!
!! SIDE EFFECTS
!!  paw_an(:)<type(paw_an_type)>= input (and eventually output) paw_an datastructures
!!
!! PARENTS
!!      m_paral_pert
!!
!! CHILDREN
!!
!! SOURCE

subroutine paw_an_redistribute(paw_an,mpi_comm_in,mpi_comm_out,&
&                 natom,mpi_atmtab_in,mpi_atmtab_out,paw_an_out,&
&                 SendAtomProc,SendAtomList,RecvAtomProc,RecvAtomList)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpi_comm_in,mpi_comm_out
 integer,optional,intent(in) :: natom
!arrays
 integer,intent(in),optional,target :: mpi_atmtab_in(:),mpi_atmtab_out(:)
 type(paw_an_type),allocatable,intent(inout) :: paw_an(:)
 type(paw_an_type),pointer,optional :: paw_an_out(:)   !vz_i
 integer,intent(in),optional :: SendAtomProc(:),SendAtomList(:),RecvAtomProc(:),RecvAtomList(:)

!Local variables-------------------------------
!scalars

 integer :: algo_option,i1,iat_in,iat_out,iatom,ierr,iircv,iisend,imsg,imsg_current,imsg1
 integer :: iproc_rcv,iproc_send,ireq,me_exch,mpi_comm_exch,my_natom_in,my_natom_out,my_tag,natom_tot,nb_msg
 integer :: nb_dp,nb_int,nbmsg_incoming,nbrecvmsg,nbsend,nbsendreq,nbsent,nbrecv,next,npaw_an_sent
 integer :: nproc_in,nproc_out
 logical :: flag,in_place,message_yet_prepared,my_atmtab_in_allocated,my_atmtab_out_allocated,paral_atom
!arrays
 integer :: buf_size(3),request1(3)
 integer,pointer :: my_atmtab_in(:),my_atmtab_out(:)
 integer,allocatable :: atmtab_send(:),atm_indx_in(:),atm_indx_out(:),buf_int1(:),From(:),request(:)
 integer,allocatable,target:: buf_int(:)
 integer,pointer :: buf_ints(:)
 logical, allocatable :: msg_pick(:)
 real(dp),allocatable :: buf_dp1(:)
 real(dp),allocatable,target :: buf_dp(:)
 real(dp),pointer :: buf_dps(:)
 type(coeffi1_type),target,allocatable :: tab_buf_int(:),tab_buf_atom(:)
 type(coeff1_type),target,allocatable :: tab_buf_dp(:)
 type(paw_an_type),allocatable :: paw_an_all(:)
 type(paw_an_type),pointer :: paw_an_out1(:)

! *************************************************************************

!@paw_an_type

 in_place=(.not.present(paw_an_out))
 my_natom_in=size(paw_an)

!If not "in_place", destroy the output datastructure
 if (.not.in_place) then
   if (associated(paw_an_out)) then
     call paw_an_free(paw_an_out)
     LIBPAW_DATATYPE_DEALLOCATE(paw_an_out)
   end if
 end if

!Special sequential case
 if (mpi_comm_in==xmpi_comm_self.and.mpi_comm_out==xmpi_comm_self) then
   if ((.not.in_place).and.(my_natom_in>0)) then
     LIBPAW_DATATYPE_ALLOCATE(paw_an_out,(my_natom_in))
     call paw_an_nullify(paw_an_out)
     call paw_an_copy(paw_an,paw_an_out)
   end if
   return
 end if

!Get total natom
 if (present(natom)) then
   natom_tot=natom
 else
   natom_tot=my_natom_in
   call xmpi_sum(natom_tot,mpi_comm_in,ierr)
 end if

!Select input distribution
 if (present(mpi_atmtab_in)) then
   my_atmtab_in => mpi_atmtab_in
   my_atmtab_in_allocated=.false.
 else
   call get_my_atmtab(mpi_comm_in,my_atmtab_in,my_atmtab_in_allocated,&
&                     paral_atom,natom_tot,my_natom_in)
 end if

!Select output distribution
 if (present(mpi_atmtab_out)) then
   my_natom_out=size(mpi_atmtab_out)
   my_atmtab_out => mpi_atmtab_out
   my_atmtab_out_allocated=.false.
 else
   call get_my_natom(mpi_comm_out,my_natom_out,natom_tot)
   call get_my_atmtab(mpi_comm_out,my_atmtab_out,my_atmtab_out_allocated,&
&                     paral_atom,natom_tot)
 end if

!Select algo according to optional input arguments
 algo_option=1
 if (present(SendAtomProc).and.present(SendAtomList).and.&
&    present(RecvAtomProc).and.present(RecvAtomList)) algo_option=2


!Brute force algorithm (allgather + scatter)
!---------------------------------------------------------
 if (algo_option==1) then

   LIBPAW_DATATYPE_ALLOCATE(paw_an_all,(natom_tot))
   call paw_an_nullify(paw_an_all)
   call paw_an_copy(paw_an,paw_an_all,comm_atom=mpi_comm_in,mpi_atmtab=my_atmtab_in)
   if (in_place) then
     call paw_an_free(paw_an)
     LIBPAW_DATATYPE_DEALLOCATE(paw_an)
     LIBPAW_DATATYPE_ALLOCATE(paw_an,(my_natom_out))
     call paw_an_nullify(paw_an)
     call paw_an_copy(paw_an_all,paw_an,comm_atom=mpi_comm_out,mpi_atmtab=my_atmtab_out)
   else
     LIBPAW_DATATYPE_ALLOCATE(paw_an_out,(my_natom_out))
     call paw_an_nullify(paw_an_out)
     call paw_an_copy(paw_an_all,paw_an_out,comm_atom=mpi_comm_out,mpi_atmtab=my_atmtab_out)
   end if
   call paw_an_free(paw_an_all)
   LIBPAW_DATATYPE_DEALLOCATE(paw_an_all)


!Asynchronous algorithm (asynchronous communications)
!---------------------------------------------------------
 else if (algo_option==2) then

   nbsend=size(SendAtomProc) ; nbrecv=size(RecvAtomProc)

   if (in_place) then
     if (my_natom_out > 0) then
       LIBPAW_DATATYPE_ALLOCATE(paw_an_out1,(my_natom_out))
       call paw_an_nullify(paw_an_out1)
     else
       LIBPAW_DATATYPE_ALLOCATE(paw_an_out1,(0))
     end if
   else
     LIBPAW_DATATYPE_ALLOCATE(paw_an_out,(my_natom_out))
     call paw_an_nullify(paw_an_out)
     paw_an_out1=>paw_an_out
   end if

   nproc_in=xmpi_comm_size(mpi_comm_in)
   nproc_out=xmpi_comm_size(mpi_comm_out)
   if (nproc_in<=nproc_out) mpi_comm_exch=mpi_comm_out
   if (nproc_in>nproc_out) mpi_comm_exch=mpi_comm_in
   me_exch=xmpi_comm_rank(mpi_comm_exch)

!  Dimension put to the maximum to send
   LIBPAW_ALLOCATE(atmtab_send,(nbsend))
   LIBPAW_ALLOCATE(atm_indx_in,(natom_tot))
   atm_indx_in=-1
   do iatom=1,my_natom_in
     atm_indx_in(my_atmtab_in(iatom))=iatom
   end do
   LIBPAW_ALLOCATE(atm_indx_out,(natom_tot))
   atm_indx_out=-1
   do iatom=1,my_natom_out
     atm_indx_out(my_atmtab_out(iatom))=iatom
   end do

   LIBPAW_DATATYPE_ALLOCATE(tab_buf_int,(nbsend))
   LIBPAW_DATATYPE_ALLOCATE(tab_buf_dp,(nbsend))
   LIBPAW_DATATYPE_ALLOCATE(tab_buf_atom,(nbsend))
   LIBPAW_ALLOCATE(request,(3*nbsend))

!  A send buffer in an asynchrone communication couldn't be deallocate before it has been receive
   nbsent=0 ; ireq=0 ; iisend=0 ; nbsendreq=0 ; nb_msg=0
   do iisend=1,nbsend
     iproc_rcv=SendAtomProc(iisend)
     next=-1
     if (iisend < nbsend) next=SendAtomProc(iisend+1)
     if (iproc_rcv /= me_exch) then
       nbsent=nbsent+1
       atmtab_send(nbsent)=SendAtomList(iisend) ! we groups the atoms sends to the same process
       if (iproc_rcv /= next) then
         if (nbsent > 0) then
!          Check if message has been yet prepared
           message_yet_prepared=.false.
           do imsg=1,nb_msg
             if (size(tab_buf_atom(imsg)%value) /= nbsent) then
               cycle
             else
               do imsg1=1,nbsent
                 if (tab_buf_atom(imsg)%value(imsg1)/=atmtab_send(imsg1)) exit
                 message_yet_prepared=.true.
                 imsg_current=imsg
               end do
             end if
           end do
!          Create the message
           if (.not.message_yet_prepared) then
             nb_msg=nb_msg+1
             call paw_an_isendreceive_fillbuffer( &
&                   paw_an,atmtab_send,atm_indx_in,nbsent,buf_int,nb_int,buf_dp,nb_dp)
             LIBPAW_ALLOCATE(tab_buf_int(nb_msg)%value,(nb_int))
             LIBPAW_ALLOCATE(tab_buf_dp(nb_msg)%value,(nb_dp))
             tab_buf_int(nb_msg)%value(1:nb_int)=buf_int(1:nb_int)
             tab_buf_dp(nb_msg)%value(1:nb_dp)=buf_dp(1:nb_dp)
             LIBPAW_DEALLOCATE(buf_int)
             LIBPAW_DEALLOCATE(buf_dp)
             LIBPAW_ALLOCATE(tab_buf_atom(nb_msg)%value, (nbsent))
             tab_buf_atom(nb_msg)%value(1:nbsent)=atmtab_send(1:nbsent)
             imsg_current=nb_msg
           end if
!          Communicate
           buf_size(1)=size(tab_buf_int(imsg_current)%value)
           buf_size(2)=size(tab_buf_dp(imsg_current)%value)
           buf_size(3)=nbsent
           buf_ints=>tab_buf_int(imsg_current)%value
           buf_dps=>tab_buf_dp(imsg_current)%value
           my_tag=300
           ireq=ireq+1
           call xmpi_isend(buf_size,iproc_rcv,my_tag,mpi_comm_exch,request(ireq),ierr)
           my_tag=301
           ireq=ireq+1
           call xmpi_isend(buf_ints,iproc_rcv,my_tag,mpi_comm_exch,request(ireq),ierr)
           my_tag=302
           ireq=ireq+1
           call xmpi_isend(buf_dps,iproc_rcv,my_tag,mpi_comm_exch,request(ireq),ierr)
           nbsendreq=ireq
           nbsent=0
         end if
       end if
     else ! Just a renumbering, not a sending
       iat_in=atm_indx_in(SendAtomList(iisend))
       iat_out=atm_indx_out(my_atmtab_in(iat_in))
       call paw_an_copy(paw_an(iat_in:iat_in),paw_an_out1(iat_out:iat_out))
       nbsent=0
     end if
   end do

   LIBPAW_ALLOCATE(From,(nbrecv))
   From(:)=-1 ; nbrecvmsg=0
   do iircv=1,nbrecv
     iproc_send=RecvAtomProc(iircv) !receive from  ( RcvAtomProc is sorted by growing process )
     next=-1
     if (iircv < nbrecv) next=RecvAtomProc(iircv+1)
     if (iproc_send /= me_exch .and. iproc_send/=next) then
       nbrecvmsg=nbrecvmsg+1
       From(nbrecvmsg)=iproc_send
     end if
   end do

   LIBPAW_ALLOCATE(msg_pick,(nbrecvmsg))
   msg_pick=.false.
   nbmsg_incoming=nbrecvmsg
   do while (nbmsg_incoming > 0)
     do i1=1,nbrecvmsg
       if (.not.msg_pick(i1)) then
         iproc_send=From(i1)
         flag=.false.
         my_tag=300
         call xmpi_iprobe(iproc_send,my_tag,mpi_comm_exch,flag,ierr)
         if (flag) then
           msg_pick(i1)=.true.
           call xmpi_irecv(buf_size,iproc_send,my_tag,mpi_comm_exch,request1(1),ierr)
           call xmpi_wait(request1(1),ierr)
           nb_int=buf_size(1)
           nb_dp=buf_size(2)
           npaw_an_sent=buf_size(3)
           LIBPAW_ALLOCATE(buf_int1,(nb_int))
           LIBPAW_ALLOCATE(buf_dp1 ,(nb_dp))
           my_tag=301
           call xmpi_irecv(buf_int1,iproc_send,my_tag,mpi_comm_exch,request1(2),ierr)
           my_tag=302
           call xmpi_irecv(buf_dp1,iproc_send,my_tag,mpi_comm_exch,request1(3),ierr)
           call xmpi_waitall(request1(2:3),ierr)
           call paw_an_isendreceive_getbuffer(paw_an_out1,npaw_an_sent,atm_indx_out,buf_int1,buf_dp1)
           nbmsg_incoming=nbmsg_incoming-1
           LIBPAW_DEALLOCATE(buf_int1)
           LIBPAW_DEALLOCATE(buf_dp1)
         end if
       end if
     end do
   end do
   LIBPAW_DEALLOCATE(msg_pick)

   if (in_place) then
     call paw_an_free(paw_an)
     LIBPAW_DATATYPE_DEALLOCATE(paw_an)
     LIBPAW_DATATYPE_ALLOCATE(paw_an,(my_natom_out))
     call paw_an_nullify(paw_an)
     call paw_an_copy(paw_an_out1,paw_an)
     call paw_an_free(paw_an_out1)
    LIBPAW_DATATYPE_DEALLOCATE(paw_an_out1)
   end if

!  Wait for deallocating arrays that all sending operations has been realized
   if (nbsendreq > 0) then
     call xmpi_waitall(request(1:nbsendreq),ierr)
   end if

!  Deallocate buffers
   do i1=1,nb_msg
     LIBPAW_DEALLOCATE(tab_buf_int(i1)%value)
     LIBPAW_DEALLOCATE(tab_buf_dp(i1)%value)
     LIBPAW_DEALLOCATE(tab_buf_atom(i1)%value)
   end do
   LIBPAW_DATATYPE_DEALLOCATE(tab_buf_int)
   LIBPAW_DATATYPE_DEALLOCATE(tab_buf_dp)
   LIBPAW_DATATYPE_DEALLOCATE(tab_buf_atom)
   LIBPAW_DEALLOCATE(From)
   LIBPAW_DEALLOCATE(request)
   LIBPAW_DEALLOCATE(atmtab_send)
   LIBPAW_DEALLOCATE(atm_indx_in)
   LIBPAW_DEALLOCATE(atm_indx_out)

 end if !algo_option

!Eventually release temporary pointers
 call free_my_atmtab(my_atmtab_in,my_atmtab_in_allocated)
 call free_my_atmtab(my_atmtab_out,my_atmtab_out_allocated)

end subroutine paw_an_redistribute
!!***

!----------------------------------------------------------------------

!!****f* m_paw_an/paw_an_reset_flags
!! NAME
!! paw_an_reset_flags
!!
!! FUNCTION
!!  Set all paw_an flags to 1 (force the recomputation of all arrays)
!!
!! SIDE EFFECTS
!!  Paw_an<type(Paw_an_type)>=paw_an datastructure
!!
!! PARENTS
!!      dfpt_nstpaw,dfpt_scfcv,scfcv
!!
!! CHILDREN
!!
!! SOURCE

subroutine paw_an_reset_flags(Paw_an)

!Arguments ------------------------------------
!arrays
 type(Paw_an_type),intent(inout) :: Paw_an(:)

!Local variables-------------------------------
 integer :: iat,natom

! *************************************************************************

!@Paw_an_type

 natom=SIZE(Paw_an);if (natom==0) return
 do iat=1,natom
   if (Paw_an(iat)%has_kxc     >0) Paw_an(iat)%has_kxc     =1
   if (Paw_an(iat)%has_k3xc    >0) Paw_an(iat)%has_k3xc    =1
   if (Paw_an(iat)%has_vhartree>0) Paw_an(iat)%has_vhartree=1
   if (Paw_an(iat)%has_vxc     >0) Paw_an(iat)%has_vxc     =1
   if (Paw_an(iat)%has_vxctau  >0) Paw_an(iat)%has_vxctau  =1
   if (Paw_an(iat)%has_vxcval  >0) Paw_an(iat)%has_vxcval  =1
   if (Paw_an(iat)%has_vxc_ex  >0) Paw_an(iat)%has_vxc_ex  =1
 end do

end subroutine paw_an_reset_flags
!!***

!----------------------------------------------------------------------

!!****f* m_paw_an/paw_an_isendreceive_getbuffer
!! NAME
!!  paw_an_isendreceive_getbuffer
!!
!! FUNCTION
!!  Fill a paw_an structure with the buffers received in a receive operation
!!  This buffer should have been first extracted by a call to paw_an_isendreceive_fillbuffer
!!
!! INPUTS
!!  atm_indx_recv(1:total number of atoms)= array for receive operation
!!                 Given an index of atom in global numbering, give its index
!!                 in the table of atoms treated by current processor
!!                 or -1 if the atoms is not treated by current processor
!!  buf_int= buffer of receive integers
!!  buf_dp= buffer of receive double precision numbers
!!  npaw_an_send= number of sent atoms
!!
!! OUTPUT
!!  paw_an= output datastructure filled with buffers receive in a receive operation
!!
!! PARENTS
!!      m_paw_an
!!
!! CHILDREN
!!
!! SOURCE

subroutine paw_an_isendreceive_getbuffer(paw_an,npaw_an_send,atm_indx_recv,buf_int,buf_dp)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npaw_an_send
!arrays
 integer,intent(in) ::atm_indx_recv(:),buf_int(:)
 real(dp),intent(in):: buf_dp(:)
 type(paw_an_type),target,intent(inout) :: paw_an(:)

!Local variables-------------------------------
!scalars
 integer :: buf_int_size,buf_dp_size,cplx_mesh_size,has_lm_select,i1,i2
 integer :: iat,iatot,ij,indx_int,indx_dp,lm_size,nkxc1,nk3xc1,nspden,v_size
 character(len=500) :: msg
 type(paw_an_type),pointer :: paw_an1
!arrays

! *********************************************************************

 buf_int_size=size(buf_int)
 buf_dp_size=size(buf_dp)
 indx_int=1; indx_dp=1

 do ij=1,npaw_an_send
   iatot=buf_int(indx_int); indx_int=indx_int+1
   iat= atm_indx_recv(iatot)
   paw_an1=>paw_an(iat)
   paw_an1%itypat=buf_int(indx_int); indx_int=indx_int+1
   paw_an1%nspden=buf_int(indx_int); indx_int=indx_int+1
   paw_an1%cplex=buf_int(indx_int); indx_int=indx_int+1
   paw_an1%mesh_size=buf_int(indx_int); indx_int=indx_int+1
   paw_an1%angl_size=buf_int(indx_int); indx_int=indx_int+1
   paw_an1%lm_size=buf_int(indx_int); indx_int=indx_int+1
   paw_an1%nkxc1=buf_int(indx_int); indx_int=indx_int+1
   paw_an1%nk3xc1=buf_int(indx_int); indx_int=indx_int+1
   paw_an1%has_vxc=buf_int(indx_int); indx_int=indx_int+1
   paw_an1%has_vxctau=buf_int(indx_int); indx_int=indx_int+1
   paw_an1%has_kxc=buf_int(indx_int); indx_int=indx_int+1
   paw_an1%has_k3xc=buf_int(indx_int); indx_int=indx_int+1
   paw_an1%has_vxcval=buf_int(indx_int); indx_int=indx_int+1
   paw_an1%has_vxc_ex=buf_int(indx_int); indx_int=indx_int+1
   paw_an1%has_vhartree=buf_int(indx_int); indx_int=indx_int+1
   v_size=buf_int(indx_int); indx_int=indx_int+1
   has_lm_select=buf_int(indx_int); indx_int=indx_int+1
   nspden=paw_an1%nspden
   lm_size=paw_an1%lm_size
   nkxc1=paw_an1%nkxc1
   nk3xc1=paw_an1%nk3xc1
   cplx_mesh_size=paw_an1%cplex*paw_an1%mesh_size
   if (has_lm_select==1) then
     LIBPAW_ALLOCATE(paw_an1%lmselect,(lm_size))
     if (lm_size>0) then
       do i1=1,lm_size
         if (buf_int(indx_int)==1) then
           paw_an1%lmselect(i1)=.TRUE.;indx_int=indx_int+1
         else
           paw_an1%lmselect(i1)=.FALSE.;indx_int=indx_int+1
         end if
       end do
     end if
   end if
   if (paw_an1%has_vxc>0) then
     LIBPAW_ALLOCATE(paw_an1%vxc1,(cplx_mesh_size,v_size,nspden))
     LIBPAW_ALLOCATE(paw_an1%vxct1,(cplx_mesh_size,v_size,nspden))
     if (paw_an1%has_vxc==2) then
       do i1=1,nspden
         do i2=1,v_size
           paw_an1%vxc1(:,i2,i1)=buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)
           indx_dp=indx_dp+cplx_mesh_size
         end do
       end do
       do i1=1,nspden
         do i2=1,v_size
           paw_an1%vxct1(:,i2,i1)=buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)
           indx_dp=indx_dp+cplx_mesh_size
         end do
       end do
     end if
   end if
   if (paw_an1%has_vxctau>0) then
     LIBPAW_ALLOCATE(paw_an1%vxctau1,(cplx_mesh_size,v_size,nspden))
     LIBPAW_ALLOCATE(paw_an1%vxcttau1,(cplx_mesh_size,v_size,nspden))
     if (paw_an1%has_vxctau==2) then
       do i1=1,nspden
         do i2=1,v_size
           paw_an1%vxctau1(:,i2,i1)=buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)
           indx_dp=indx_dp+cplx_mesh_size
         end do
       end do
       do i1=1,nspden
         do i2=1,v_size
           paw_an1%vxcttau1(:,i2,i1)=buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)
           indx_dp=indx_dp+cplx_mesh_size
         end do
       end do
     end if
   end if
   if (paw_an1%has_kxc>0) then
     LIBPAW_ALLOCATE(paw_an1%kxc1,(cplx_mesh_size,v_size,nkxc1))
     LIBPAW_ALLOCATE(paw_an1%kxct1,(cplx_mesh_size,v_size,nkxc1))
     if (paw_an1%has_kxc==2.and.nkxc1>0) then
       do i1=1,nkxc1
         do i2=1,v_size
           paw_an1%kxc1(:,i2,i1)=buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)
           indx_dp=indx_dp+cplx_mesh_size
         end do
       end do
       do i1=1,nkxc1
         do i2=1,v_size
           paw_an1%kxct1(:,i2,i1)=buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)
           indx_dp=indx_dp+cplx_mesh_size
         end do
       end do
     end if
   end if
   if (paw_an1%has_k3xc>0) then
     LIBPAW_ALLOCATE(paw_an1%k3xc1,(cplx_mesh_size,v_size,nk3xc1))
     LIBPAW_ALLOCATE(paw_an1%k3xct1,(cplx_mesh_size,v_size,nk3xc1))
     if (paw_an1%has_k3xc==2.and.nk3xc1>0) then
       do i1=1,nk3xc1
         do i2=1,v_size
           paw_an1%k3xc1(:,i2,i1)=buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)
           indx_dp=indx_dp+cplx_mesh_size
         end do
       end do
       do i1=1,nk3xc1
         do i2=1,v_size
           paw_an1%k3xct1(:,i2,i1)=buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)
           indx_dp=indx_dp+cplx_mesh_size
         end do
       end do
     end if
   end if
   if (paw_an1%has_vxcval>0) then
     LIBPAW_ALLOCATE(paw_an1%vxc1_val,(cplx_mesh_size,v_size,nspden))
     LIBPAW_ALLOCATE(paw_an1%vxct1_val,(cplx_mesh_size,v_size,nspden))
     if (paw_an1%has_vxcval==2) then
       do i1=1,nspden
         do i2=1,v_size
           paw_an1%vxc1_val(:,i2,i1)=buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)
           indx_dp=indx_dp+cplx_mesh_size
         end do
       end do
       do i1=1,nspden
         do i2=1,v_size
           paw_an1%vxct1_val(:,i2,i1)=buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)
           indx_dp=indx_dp+cplx_mesh_size
         end do
       end do
     end if
   end if
   if (paw_an1%has_vxc_ex>0) then
     LIBPAW_ALLOCATE(paw_an1%vxc_ex,(cplx_mesh_size,v_size,nspden))
     if (paw_an1%has_vxc_ex==2) then
       do i1=1,nspden
         do i2=1,v_size
           paw_an1%vxc_ex(:,i2,i1)=buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)
           indx_dp=indx_dp+cplx_mesh_size
         end do
       end do
     end if
   end if
   if (paw_an1%has_vhartree>0) then
      LIBPAW_ALLOCATE(paw_an1%vh1,(cplx_mesh_size,lm_size,nspden))
      LIBPAW_ALLOCATE(paw_an1%vht1,(cplx_mesh_size,lm_size,nspden))
      if (paw_an1%has_vhartree==2) then
        do i1=1,nspden
          do i2=1,lm_size
            paw_an1%vh1(:,i2,i1)=buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)
            indx_dp=indx_dp+cplx_mesh_size
          end do
        end do
        do i1=1,nspden
          do i2=1,lm_size
            paw_an1%vht1(:,i2,i1)=buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)
            indx_dp=indx_dp+cplx_mesh_size
          end do
        end do
     end if
   end if
 end do ! iat
 if ((indx_int/=1+buf_int_size).or.(indx_dp/=1+buf_dp_size)) then
   write(msg,'(a,i10,a,i10)') 'Wrong buffer sizes: buf_int_size=',buf_int_size,' buf_dp_size=',buf_dp_size
   MSG_BUG(msg)
 end if

end subroutine paw_an_isendreceive_getbuffer
!!***

!----------------------------------------------------------------------


!!****f* m_paw_an/paw_an_isendreceive_fillbuffer
!! NAME
!!  paw_an_isendreceive_fillbuffer
!!
!! FUNCTION
!!  Extract from paw_an and from the global index of atoms
!!  the buffers to send in a sending operation
!!  This function has to be coupled with a call to paw_ij_isendreceive_getbuffer
!! INPUTS
!!  atm_indx_send(1:total number of atoms)= array for send operation,
!!                 Given an index of atom in global numbering, give its index
!!                 in the table of atoms treated by current processor
!!                 or -1 if the atoms is not treated by current processor
!!  npaw_an_send= number of sent atoms
!!  paw_an= data structure from which are extract buffer int and buffer dp
!!
!! OUTPUT
!! buf_int : buffer of integers to be send in a send operation
!! buf_int_size : size of buffer of integers to be send in a send operation
!! buf_dp : buffer of double precision numbers to be send in a send operation
!! buf_dp_size : size of buffer of double precision numbers to be send in a send operation
!!
!! PARENTS
!!      m_paw_an
!!
!! CHILDREN
!!
!! SOURCE

subroutine paw_an_isendreceive_fillbuffer(paw_an, atmtab_send,atm_indx_send,npaw_an_send,&
&                                         buf_int,buf_int_size,buf_dp,buf_dp_size)

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: buf_int_size,buf_dp_size
 integer,intent(in) :: npaw_an_send
!arrays
 integer,intent(in) :: atmtab_send(:),atm_indx_send(:)
 integer,allocatable,intent(out) :: buf_int(:)
 real(dp),allocatable,intent(out):: buf_dp(:)
 type(paw_an_type),target,intent(in) :: paw_an(:)

!Local variables-------------------------------
!scalars
 integer :: cplx_mesh_size,i1,i2,iatom_tot,ij,indx_int,indx_dp
 integer :: ipaw_an_send,lm_size,nspden,v_size
 character(len=500) :: msg
 type(Paw_an_type),pointer :: paw_an1
!arrays

! *********************************************************************

!Compute sizes of buffers
 buf_int_size=0 ; buf_dp_size=0
 do ipaw_an_send=1,npaw_an_send
   iatom_tot=atmtab_send(ipaw_an_send)
   ij = atm_indx_send(iatom_tot)
   paw_an1=>paw_an(ij)
   buf_int_size=buf_int_size+17+size(paw_an1%lmselect)
   if (paw_an1%has_vxc==2) then
     buf_dp_size=buf_dp_size+size(paw_an1%vxc1)
     buf_dp_size=buf_dp_size+size(paw_an1%vxct1)
   end if
   if (paw_an1%has_vxctau==2) then
     buf_dp_size=buf_dp_size+size(paw_an1%vxctau1)
     buf_dp_size=buf_dp_size+size(paw_an1%vxcttau1)
   end if
   if (paw_an1%has_kxc==2) then
     buf_dp_size=buf_dp_size+size(paw_an1%kxc1)
     buf_dp_size=buf_dp_size+size(paw_an1%kxct1)
   end if
   if (paw_an1%has_k3xc==2) then
     buf_dp_size=buf_dp_size+size(paw_an1%k3xc1)
     buf_dp_size=buf_dp_size+size(paw_an1%k3xct1)
   end if
   if (paw_an1%has_vxcval==2) then
     buf_dp_size=buf_dp_size+size(paw_an1%vxc1_val)
     buf_dp_size=buf_dp_size+size(paw_an1%vxct1_val)
   end if
   if (paw_an1%has_vxc_ex==2) then
     buf_dp_size=buf_dp_size+size(paw_an1%vxc_ex)
   end if
   if (paw_an1%has_vhartree==2) then
     buf_dp_size=buf_dp_size+size(paw_an1%vh1)
     buf_dp_size=buf_dp_size+size(paw_an1%vht1)
   end if
 end do

!Fill in input buffers
 LIBPAW_ALLOCATE(buf_int,(buf_int_size))
 LIBPAW_ALLOCATE(buf_dp ,(buf_dp_size))
 indx_int=1;indx_dp=1
 do ipaw_an_send=1,npaw_an_send
   iatom_tot=atmtab_send(ipaw_an_send)
   ij = atm_indx_send(iatom_tot)
   paw_an1=>paw_an(ij)
   buf_int(indx_int)=iatom_tot; indx_int=indx_int+1
   buf_int(indx_int)=paw_an1%itypat; indx_int=indx_int+1
   buf_int(indx_int)=paw_an1%nspden; indx_int=indx_int+1
   buf_int(indx_int)=paw_an1%cplex; indx_int=indx_int+1
   buf_int(indx_int)=paw_an1%mesh_size; indx_int=indx_int+1
   buf_int(indx_int)=paw_an1%angl_size; indx_int=indx_int+1
   buf_int(indx_int)=paw_an1%lm_size; indx_int=indx_int+1
   buf_int(indx_int)=paw_an1%nkxc1; indx_int=indx_int+1
   buf_int(indx_int)=paw_an1%nk3xc1; indx_int=indx_int+1
   buf_int(indx_int)=paw_an1%has_vxc; indx_int=indx_int+1
   buf_int(indx_int)=paw_an1%has_vxctau; indx_int=indx_int+1
   buf_int(indx_int)=paw_an1%has_kxc; indx_int=indx_int+1
   buf_int(indx_int)=paw_an1%has_k3xc; indx_int=indx_int+1
   buf_int(indx_int)=paw_an1%has_vxcval; indx_int=indx_int+1
   buf_int(indx_int)=paw_an1%has_vxc_ex; indx_int=indx_int+1
   buf_int(indx_int)=paw_an1%has_vhartree; indx_int=indx_int+1
   v_size=0
   if (paw_an1%has_vxc>0) then
     v_size=size(paw_an1%vxc1,2)
   else if (paw_an1%has_vxctau>0) then
     v_size=size(paw_an1%vxctau1,2)
   else if (paw_an1%has_kxc>0) then
     v_size=size(paw_an1%kxc1,2)
   else if (paw_an1%has_k3xc>0) then
     v_size=size(paw_an1%k3xc1,2)
   else if (paw_an1%has_vxcval>0) then
     v_size=size(paw_an1%vxc1_val,2)
   else if (paw_an1%has_vxc_ex>0) then
     v_size=size(paw_an1%vxc_ex,2)
   else if (paw_an1%has_vhartree>0) then
     v_size=size(paw_an1%vh1,2)
   end if
   buf_int(indx_int)=v_size;indx_int=indx_int+1
   if (allocated(paw_an1%lmselect)) then
     buf_int(indx_int)=1;indx_int=indx_int+1
   else
     buf_int(indx_int)=0;indx_int=indx_int+1
   end if
   nspden=paw_an1%nspden
   lm_size=paw_an1%lm_size
   cplx_mesh_size=paw_an1%cplex*paw_an1%mesh_size
   if (lm_size>0) then
     if (allocated(paw_an1%lmselect)) then
       do i1=1,lm_size
         if (paw_an1%lmselect(i1)) then
           buf_int(indx_int)=1
         else
           buf_int(indx_int)=0
         end if
         indx_int=indx_int+1
       end do
     end if
   end if
   if (paw_an1%has_vxc==2) then
     do i1=1,nspden
       do i2=1,v_size
         buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)=paw_an1%vxc1(:,i2,i1)
         indx_dp=indx_dp+cplx_mesh_size
       end do
     end do
     do i1=1,nspden
       do i2=1,v_size
         buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)=paw_an1%vxct1(:,i2,i1)
         indx_dp=indx_dp+cplx_mesh_size
       end do
     end do
   end if
   if (paw_an1%has_vxctau==2) then
     do i1=1,nspden
       do i2=1,v_size
         buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)=paw_an1%vxctau1(:,i2,i1)
         indx_dp=indx_dp+cplx_mesh_size
       end do
     end do
     do i1=1,nspden
       do i2=1,v_size
         buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)=paw_an1%vxcttau1(:,i2,i1)
         indx_dp=indx_dp+cplx_mesh_size
       end do
     end do
   end if
   if (paw_an1%has_kxc==2.and.paw_an1%nkxc1>0) then
     do i1=1,paw_an1%nkxc1
       do i2=1,v_size
         buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)=paw_an1%kxc1(:,i2,i1)
         indx_dp=indx_dp+cplx_mesh_size
       end do
     end do
     do i1=1,paw_an1%nkxc1
       do i2=1,v_size
         buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)=paw_an1%kxct1(:,i2,i1)
         indx_dp=indx_dp+cplx_mesh_size
       end do
     end do
   end if
   if (paw_an1%has_k3xc==2.and.paw_an1%nk3xc1>0) then
     do i1=1,paw_an1%nk3xc1
       do i2=1,v_size
         buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)=paw_an1%k3xc1(:,i2,i1)
         indx_dp=indx_dp+cplx_mesh_size
       end do
     end do
     do i1=1,paw_an1%nk3xc1
       do i2=1,v_size
         buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)=paw_an1%k3xct1(:,i2,i1)
         indx_dp=indx_dp+cplx_mesh_size
       end do
     end do
   end if

   if (paw_an1%has_vxcval==2) then
     do i1=1,nspden
       do i2=1,v_size
         buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)=paw_an1%vxc1_val(:,i2,i1)
         indx_dp=indx_dp+cplx_mesh_size
       end do
     end do
     do i1=1,nspden
       do i2=1,v_size
         buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)=paw_an1%vxct1_val(:,i2,i1)
         indx_dp=indx_dp+cplx_mesh_size
       end do
     end do
   end if
   if (paw_an1%has_vxc_ex==2) then
     do i1=1,nspden
       do i2=1,v_size
         buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)=paw_an1%vxc_ex(:,i2,i1)
         indx_dp=indx_dp+cplx_mesh_size
       end do
     end do
   end if
   if (paw_an1%has_vhartree==2) then
     do i1=1,nspden
       do i2=1,lm_size
         buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)=paw_an1%vh1(:,i2,i1)
         indx_dp=indx_dp+cplx_mesh_size
       end do
     end do
     do i1=1,nspden
       do i2=1,lm_size
         buf_dp(indx_dp:indx_dp+cplx_mesh_size-1)=paw_an1%vht1(:,i2,i1)
         indx_dp=indx_dp+cplx_mesh_size
       end do
     end do
   end if
 end do
 if ((indx_int-1/=buf_int_size).or.(indx_dp-1/=buf_dp_size)) then
   write(msg,'(4(a,i10))') 'Wrong buffer sizes: buf_int =',buf_int_size,'/',indx_int-1,&
   &                                           ' buf_dp =',buf_dp_size ,'/',indx_dp-1
   MSG_BUG(msg)
 end if

end subroutine paw_an_isendreceive_fillbuffer
!!***

!----------------------------------------------------------------------

END MODULE m_paw_an
!!***
