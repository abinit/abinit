!!****m* ABINIT/m_results_out
!! NAME
!!  m_results_out
!!
!! FUNCTION
!!  This module provides the definition of the results_out_type used
!!  to store results from GS calculations.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2020 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! TODO
!! One should replace the 'pointer' by 'allocatable'. This was tried, in October 2014,
!! but Petrus_nag complained (test v67mbpt t31...t34), and also max2 (paral#08 np=10).
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_results_out

 use defs_basis
 use m_dtset
 use m_errors
 use m_abicore
 use m_xmpi

 use defs_abitypes, only : MPI_type

 implicit none

 private

! public procedures.
 public :: init_results_out
 public :: destroy_results_out
 public :: copy_results_out
 public :: gather_results_out
!!***

!!****t* m_results_out/results_out_type
!! NAME
!! results_out_type
!!
!! FUNCTION
!! This structured datatype contains a subset of the results of a GS
!! calculation, needed to perform the so-called "internal tests", and
!! to perform the timing analysis
!!
!! SOURCE

 type, public :: results_out_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! Integer scalar

  integer :: natom
   ! The number of atoms for this dataset
  integer :: nimage
   ! The number of images of the cell for this dataset (treated by current proc)
  integer :: nkpt
   ! The number of k-pints for this dataset
  integer :: nocc
   ! The number of occupations for this dataset
  integer :: npsp
   ! The number of pseudopotentials
  integer :: ntypat
   ! The number of types of atoms

! Integer arrays

  integer, pointer :: npwtot(:,:)
   ! npw(mxnkpt,nimage) Full number of plane waves for each
   ! k point, computed with the "true" rprimd
   ! Not taking into account the decrease due to istwfk
   ! Not taking into account the spread of pws on different procs

! Real (real(dp)) arrays

  real(dp), pointer :: acell(:,:)
   ! acell(3,nimage)
   ! Length of primitive vectors

  real(dp), pointer :: amu(:,:)
   ! amu(ntypat,nimage)
   ! Mass of the atomic type

  real(dp), pointer :: etotal(:)
   ! etotal(nimage)
   ! Total energy (Hartree)

  real(dp), pointer :: fcart(:,:,:)
   ! fcart(3,natom,nimage) Cartesian forces (Hartree/Bohr)
   ! Forces in cartesian coordinates (Hartree)

  real(dp), pointer :: fred(:,:,:)
   ! fred(3,natom,nimage)
   ! Forces in reduced coordinates (Hartree)
   ! Actually, gradient of the total energy with respect
   ! to change of reduced coordinates

  real(dp), pointer :: intgres(:,:,:)
   ! intgres(4,natom,nimage)   ! 4 is for nspden
   ! Gradient of the total energy wrt constraints (Hartree)

  real(dp), pointer :: mixalch(:,:,:)
   ! mixalch(npsp,ntypat,nimage)   [note that in psps datastructure, the dimensioning is npspalch,ntypalch]
   ! Mixing coefficients going from the input pseudopotentials (those for alchemical mixing) to the alchemical atoms

  real(dp), pointer :: occ(:,:)
   ! occ(mxmband_upper*mxnkpt*mxnsppol,nimage)
   ! Electronic occupations

  real(dp), pointer :: rprim(:,:,:)
   ! rprim(3,3,nimage)
   ! Dimensionless real space primitive translations

  real(dp), pointer :: strten(:,:)
   ! strten(6,nimage)
   ! Stress tensor

  real(dp), pointer :: vel(:,:,:)
   ! vel(3,natom,nimage)
   ! Atomic velocities

  real(dp), pointer :: vel_cell(:,:,:)
   ! vel_cell(3,3,nimage)
   ! Cell velocities
   ! Time derivatives of dimensional primitive translations

  real(dp), pointer :: xred(:,:,:)
   ! xred(3,natom,nimage)
   ! Atomic positions in reduced coordinates

 end type results_out_type
!!***

CONTAINS

!===========================================================
!!***

!!****f* m_results_out/init_results_out
!! NAME
!!  init_results_out
!!
!! FUNCTION
!!  Init all scalars and pointers in an array of results_out datastructures
!!
!! INPUTS
!!  dtsets(:)= <type datafiles_type> contains all input variables,
!!  option_alloc=0: only allocate datastructure
!!               1: allocate and initialize the whole datastructure
!!               2: allocate datastructure and initialize only first member
!!  option_size=0: allocate results_out with a global number images
!!                  (use mxnimage=max(dtset%nimage))
!!              1: allocate results_out with a number of images per processor
!!                  (use mxnimage=max(mpi_enreg%my_nimage))
!!  mpi_enregs=information about MPI parallelization
!!  mxnimage=-optional- maximal value of nimage over datasets
!!            if this argument is present, it is used for allocations
!!            if it is not present, allocations are automatic
!!  natom= number of atoms
!!  nband= number of bands
!!  nkpt= number of k-points
!!  nsppol= number of independant spin components
!!
!! SIDE EFFECTS
!!  results_out(:)=<type(results_out_type)>=results_out datastructure array
!!
!! PARENTS
!!      abinit,m_results_out
!!
!! CHILDREN
!!      copy_results_out,init_results_out,xmpi_allgather,xmpi_allgatherv
!!      xmpi_gatherv
!!
!! SOURCE

subroutine init_results_out(dtsets,option_alloc,option_size,mpi_enregs,&
&          mxnatom,mxnband,mxnkpt,mxnpsp,mxnsppol,mxntypat,results_out)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: option_alloc,option_size
 integer,intent(in) :: mxnatom,mxnband,mxnkpt,mxnpsp,mxnsppol,mxntypat
!arrays
 type(dataset_type),intent(in) :: dtsets(:)
 type(results_out_type),intent(inout) :: results_out(:)
 type(MPI_type), intent(in) :: mpi_enregs(:)
!Local variables-------------------------------
!scalars
 integer :: dtsets_size,idt1,idt2,idt3,ii,jj,kk
 integer :: mpi_enregs_size,mxnimage_,natom_,nkpt_,nocc_
 integer :: results_out_size
! type(MPI_type) :: mpi_img
!arrays
 integer,allocatable :: img(:,:),nimage(:)
 real(dp),allocatable :: tmp(:,:)

!************************************************************************

 !@results_out_type

 dtsets_size=size(dtsets)
 results_out_size=size(results_out)
 mpi_enregs_size=size(mpi_enregs)
 if (dtsets_size/=mpi_enregs_size .or. dtsets_size/=results_out_size) then
   MSG_ERROR("init_results_out: wrong sizes (2)!")
 endif

 if (results_out_size>0) then

   idt1=lbound(results_out,1);idt2=ubound(results_out,1)
   idt3=idt2;if (option_alloc==2) idt3=idt1
   ABI_ALLOCATE(nimage,(idt1:idt2))
   nimage=0
   mxnimage_=1
   if (option_size==0) then
     do ii=idt1,idt2
       nimage(ii)=dtsets(ii)%nimage
       if (nimage(ii)>mxnimage_) mxnimage_=nimage(ii)
     end do
     if (option_alloc>0) then
       ABI_ALLOCATE(img,(mxnimage_,idt1:idt3))
       img=0
       do ii=idt1,idt3
         do jj=1,nimage(ii)
           img(jj,ii)=jj
         end do
       end do
     end if
   else
     do ii=idt1,idt2
       nimage(ii)=mpi_enregs(ii)%my_nimage
       if (nimage(ii)>mxnimage_) mxnimage_=nimage(ii)
     end do
     if (option_alloc>0) then
       ABI_ALLOCATE(img,(mxnimage_,idt1:idt3))
       img=0
       do ii=idt1,idt3
         do jj=1,nimage(ii)
           img(jj,ii)=mpi_enregs(ii)%my_imgtab(jj)
         end do
       end do
     end if
   end if

   do ii=idt1,idt2

     ABI_ALLOCATE(results_out(ii)%acell,(3,mxnimage_))
     ABI_ALLOCATE(results_out(ii)%amu,(mxntypat,mxnimage_))
     ABI_ALLOCATE(results_out(ii)%etotal,(mxnimage_))
     ABI_ALLOCATE(results_out(ii)%fcart,(3,mxnatom,mxnimage_))
     ABI_ALLOCATE(results_out(ii)%fred,(3,mxnatom,mxnimage_))
     ABI_ALLOCATE(results_out(ii)%intgres,(4,mxnatom,mxnimage_))
     ABI_ALLOCATE(results_out(ii)%mixalch,(mxnpsp,mxntypat,mxnimage_))
     ABI_ALLOCATE(results_out(ii)%npwtot,(mxnkpt,mxnimage_))
     ABI_ALLOCATE(results_out(ii)%occ,(mxnband*mxnkpt*mxnsppol,mxnimage_))
     ABI_ALLOCATE(results_out(ii)%rprim,(3,3,mxnimage_))
     ABI_ALLOCATE(results_out(ii)%strten,(6,mxnimage_))
     ABI_ALLOCATE(results_out(ii)%vel,(3,mxnatom,mxnimage_))
     ABI_ALLOCATE(results_out(ii)%vel_cell,(3,3,mxnimage_))
     ABI_ALLOCATE(results_out(ii)%xred,(3,mxnatom,mxnimage_))

     if ((option_alloc==1).or.(option_alloc==2.and.ii==idt3)) then
       results_out(ii)%nimage=nimage(ii)
       results_out(ii)%natom =mxnatom
       results_out(ii)%nkpt  =mxnkpt
       results_out(ii)%npsp  =mxnpsp
       results_out(ii)%ntypat =mxntypat
       results_out(ii)%nocc  =mxnband*mxnkpt*mxnsppol
       natom_=dtsets(ii)%natom
       nkpt_=dtsets(ii)%nkpt;if(ii==0) nkpt_=mxnkpt
       nocc_=mxnband*dtsets(ii)%nkpt*dtsets(ii)%nsppol
       results_out(ii)%nimage=nimage(ii)
       results_out(ii)%natom=natom_
       results_out(ii)%nkpt=nkpt_
       results_out(ii)%nocc=nocc_
       results_out(ii)%acell=zero
       results_out(ii)%amu=zero
       results_out(ii)%etotal(:)=zero
       results_out(ii)%fcart(:,:,:)=zero
       results_out(ii)%fred(:,:,:)=zero
       results_out(ii)%intgres(:,:,:)=zero
       results_out(ii)%mixalch(:,:,:)=zero
       results_out(ii)%occ=zero
       results_out(ii)%rprim=zero
       results_out(ii)%strten(:,:)=zero
       results_out(ii)%vel=zero
       results_out(ii)%vel_cell=zero
       results_out(ii)%xred=zero
       results_out(ii)%npwtot(:,:)=0
       if (nimage(ii)>0) then
         do jj=1,nimage(ii)
           kk=img(jj,ii)
           results_out(ii)%acell(:,jj)     =dtsets(ii)%acell_orig(:,kk)
           results_out(ii)%amu(:,jj)       =dtsets(ii)%amu_orig(:,kk)
           results_out(ii)%rprim(:,:,jj)   =dtsets(ii)%rprim_orig(:,:,kk)
           results_out(ii)%vel_cell(:,:,jj)=dtsets(ii)%vel_cell_orig(:,:,kk)
           results_out(ii)%mixalch(:,:,jj) =dtsets(ii)%mixalch_orig(:,:,kk)
           if (natom_>0) then
             ABI_ALLOCATE(tmp,(3,natom_))
             tmp(1:3,1:natom_)=dtsets(ii)%vel_orig(1:3,1:natom_,kk)
             results_out(ii)%vel(1:3,1:natom_,jj)=tmp(1:3,1:natom_)
             tmp(1:3,1:natom_)=dtsets(ii)%xred_orig(1:3,1:natom_,kk)
             results_out(ii)%xred(1:3,1:natom_,jj)=tmp(1:3,1:natom_)
             ABI_DEALLOCATE(tmp)
           end if
           if (nocc_>0) then
             results_out(ii)%occ(1:nocc_,jj)=dtsets(ii)%occ_orig(1:nocc_,kk)
           end if
         end do
       end if
     end if

   end do
   ABI_DEALLOCATE(nimage)
   !if (option_size/=0.and.option_alloc==1)  then
   if (allocated(img))  then
     ABI_DEALLOCATE(img)
   end if
 end if

end subroutine init_results_out
!!***

!----------------------------------------------------------------------

!!****f* m_results_out/destroy_results_out
!! NAME
!!  destroy_results_out
!!
!! FUNCTION
!!  Clean and destroy an array of results_out datastructures
!!
!! SIDE EFFECTS
!!  results_out(:)=<type(results_out_type)>=results_out datastructure array
!!
!! PARENTS
!!      abinit,driver
!!
!! CHILDREN
!!      copy_results_out,init_results_out,xmpi_allgather,xmpi_allgatherv
!!      xmpi_gatherv
!!
!! SOURCE

subroutine destroy_results_out(results_out)

!Arguments ------------------------------------
!arrays
 type(results_out_type),intent(inout) :: results_out(:)
!Local variables-------------------------------
!scalars
 integer :: idt1,idt2,ii,results_out_size

!************************************************************************

 !@results_out_type

 results_out_size=size(results_out)
 if (results_out_size>0) then

   idt1=lbound(results_out,1);idt2=ubound(results_out,1)
   do ii=idt1,idt2
     results_out(ii)%nimage=0
     results_out(ii)%natom=0
     results_out(ii)%nkpt=0
     results_out(ii)%nocc=0
     if (associated(results_out(ii)%acell))   then
       ABI_DEALLOCATE(results_out(ii)%acell)
     end if
     if (associated(results_out(ii)%amu))   then
       ABI_DEALLOCATE(results_out(ii)%amu)
     end if
     if (associated(results_out(ii)%etotal))  then
       ABI_DEALLOCATE(results_out(ii)%etotal)
     end if
     if (associated(results_out(ii)%fcart))   then
       ABI_DEALLOCATE(results_out(ii)%fcart)
     end if
     if (associated(results_out(ii)%fred))    then
       ABI_DEALLOCATE(results_out(ii)%fred)
     end if
     if (associated(results_out(ii)%intgres))    then
       ABI_DEALLOCATE(results_out(ii)%intgres)
     end if
     if (associated(results_out(ii)%mixalch))  then
       ABI_DEALLOCATE(results_out(ii)%mixalch)
     end if
     if (associated(results_out(ii)%npwtot))  then
       ABI_DEALLOCATE(results_out(ii)%npwtot)
     end if
     if (associated(results_out(ii)%occ))     then
       ABI_DEALLOCATE(results_out(ii)%occ)
     end if
     if (associated(results_out(ii)%rprim))   then
       ABI_DEALLOCATE(results_out(ii)%rprim)
     end if
     if (associated(results_out(ii)%strten))  then
       ABI_DEALLOCATE(results_out(ii)%strten)
     end if
     if (associated(results_out(ii)%vel))     then
       ABI_DEALLOCATE(results_out(ii)%vel)
     end if
     if (associated(results_out(ii)%vel_cell))  then
       ABI_DEALLOCATE(results_out(ii)%vel_cell)
     end if
     if (associated(results_out(ii)%xred))    then
       ABI_DEALLOCATE(results_out(ii)%xred)
     end if
   end do

 end if

end subroutine destroy_results_out
!!***

!----------------------------------------------------------------------

!!****f* m_results_out/copy_results_out
!! NAME
!!  copy_results_out
!!
!! FUNCTION
!!  Copy a results_out datastructure into another
!!
!! INPUTS
!!  results_out_in=<type(results_out_type)>=input results_out datastructure
!!
!! OUTPUT
!!  results_out_out=<type(results_out_type)>=output results_out datastructure
!!
!! PARENTS
!!      m_results_out
!!
!! CHILDREN
!!      copy_results_out,init_results_out,xmpi_allgather,xmpi_allgatherv
!!      xmpi_gatherv
!!
!! SOURCE

subroutine copy_results_out(results_out_in,results_out_out)

!Arguments ------------------------------------
!arrays
 type(results_out_type),intent(in) :: results_out_in
 type(results_out_type),intent(out) :: results_out_out
!Local variables-------------------------------
!scalars
 integer :: natom_,natom_out,nimage_,nimage_out,nkpt_,nkpt_out,npsp_,npsp_out,nocc_,nocc_out,ntypat_,ntypat_out

!************************************************************************

 !@results_out_type

 nimage_=size(results_out_in%etotal)
 natom_ =size(results_out_in%fcart,2)
 nkpt_  =size(results_out_in%npwtot,1)
 nocc_  =size(results_out_in%occ,1)
 npsp_  =size(results_out_in%mixalch,1)
 ntypat_=size(results_out_in%mixalch,2)
 nimage_out=0;if (associated(results_out_out%etotal))nimage_out=size(results_out_out%etotal)
 natom_out =0;if (associated(results_out_out%fcart)) natom_out =size(results_out_out%fcart,2)
 nkpt_out  =0;if (associated(results_out_out%npwtot))nkpt_out  =size(results_out_out%npwtot,1)
 nocc_out  =0;if (associated(results_out_out%occ))   nocc_out  =size(results_out_out%occ,1)
 npsp_out  =0;if (associated(results_out_out%mixalch))npsp_out  =size(results_out_out%mixalch,1)
 ntypat_out=0;if (associated(results_out_out%mixalch))ntypat_out=size(results_out_out%mixalch,2)

 if (nimage_>nimage_out) then
   if (associated(results_out_out%acell))   then
     ABI_DEALLOCATE(results_out_out%acell)
   end if
   if (associated(results_out_out%etotal))  then
     ABI_DEALLOCATE(results_out_out%etotal)
   end if
   if (associated(results_out_out%rprim))   then
     ABI_DEALLOCATE(results_out_out%rprim)
   end if
   if (associated(results_out_out%strten))  then
     ABI_DEALLOCATE(results_out_out%strten)
   end if
   if (associated(results_out_out%vel_cell))  then
     ABI_DEALLOCATE(results_out_out%vel_cell)
   end if
   ABI_ALLOCATE(results_out_out%acell,(3,nimage_))
   ABI_ALLOCATE(results_out_out%etotal,(nimage_))
   ABI_ALLOCATE(results_out_out%rprim,(3,3,nimage_))
   ABI_ALLOCATE(results_out_out%strten,(6,nimage_))
   ABI_ALLOCATE(results_out_out%vel_cell,(3,3,nimage_))
 end if
 if (nimage_>nimage_out.or.natom_>natom_out) then
   if (associated(results_out_out%fcart))   then
     ABI_DEALLOCATE(results_out_out%fcart)
   end if
   if (associated(results_out_out%fred))    then
     ABI_DEALLOCATE(results_out_out%fred)
   end if
   if (associated(results_out_out%intgres))    then
     ABI_DEALLOCATE(results_out_out%intgres)
   end if
   if (associated(results_out_out%vel))     then
     ABI_DEALLOCATE(results_out_out%vel)
   end if
   if (associated(results_out_out%xred))    then
     ABI_DEALLOCATE(results_out_out%xred)
   end if
   ABI_ALLOCATE(results_out_out%fcart,(3,natom_,nimage_))
   ABI_ALLOCATE(results_out_out%fred,(3,natom_,nimage_))
   ABI_ALLOCATE(results_out_out%intgres,(4,natom_,nimage_))
   ABI_ALLOCATE(results_out_out%vel,(3,natom_,nimage_))
   ABI_ALLOCATE(results_out_out%xred,(3,natom_,nimage_))
 end if
 if (nimage_>nimage_out.or.nkpt_>nkpt_out) then
   if (associated(results_out_out%npwtot))  then
     ABI_DEALLOCATE(results_out_out%npwtot)
   end if
   ABI_ALLOCATE(results_out_out%npwtot,(nkpt_,nimage_))
 end if
 if (nimage_>nimage_out.or.nocc_>nocc_out) then
   if (associated(results_out_out%occ))     then
     ABI_DEALLOCATE(results_out_out%occ)
   end if
   ABI_ALLOCATE(results_out_out%occ,(nocc_,nimage_))
 end if
 if (ntypat_>ntypat_out) then
   if (associated(results_out_out%amu))     then
     ABI_DEALLOCATE(results_out_out%amu)
   end if
   ABI_ALLOCATE(results_out_out%amu,(ntypat_,nimage_))
 end if

 if (npsp_>npsp_out.or.ntypat_>ntypat_out) then
   if (associated(results_out_out%mixalch))     then
     ABI_DEALLOCATE(results_out_out%mixalch)
   end if
   ABI_ALLOCATE(results_out_out%mixalch,(npsp_,ntypat_,nimage_))
 end if

 results_out_out%nimage=results_out_in%nimage
 results_out_out%natom =results_out_in%natom
 results_out_out%nkpt  =results_out_in%nkpt
 results_out_out%nocc  =results_out_in%nocc
 results_out_out%acell(1:3,1:nimage_)         =results_out_in%acell(1:3,1:nimage_)
 results_out_out%amu(1:ntypat_,1:nimage_)      =results_out_in%amu(1:ntypat_,1:nimage_)
 results_out_out%etotal(1:nimage_)            =results_out_in%etotal(1:nimage_)
 results_out_out%fcart(1:3,1:natom_,1:nimage_)=results_out_in%fcart(1:3,1:natom_,1:nimage_)
 results_out_out%fred(1:3,1:natom_,1:nimage_) =results_out_in%fred(1:3,1:natom_,1:nimage_)
 results_out_out%intgres(1:4,1:natom_,1:nimage_) =results_out_in%intgres(1:4,1:natom_,1:nimage_)
 results_out_out%mixalch(1:npsp_,1:ntypat_,1:nimage_)=results_out_in%mixalch(1:npsp_,1:ntypat_,1:nimage_)
 results_out_out%npwtot(1:nkpt_,1:nimage_)    =results_out_in%npwtot(1:nkpt_,1:nimage_)
 results_out_out%occ(1:nocc_,1:nimage_)       =results_out_in%occ(1:nocc_,1:nimage_)
 results_out_out%rprim(1:3,1:3,1:nimage_)     =results_out_in%rprim(1:3,1:3,1:nimage_)
 results_out_out%strten(1:6,1:nimage_)        =results_out_in%strten(1:6,1:nimage_)
 results_out_out%xred(1:3,1:natom_,1:nimage_) =results_out_in%xred(1:3,1:natom_,1:nimage_)
 results_out_out%vel(1:3,1:natom_,1:nimage_)  =results_out_in%vel(1:3,1:natom_,1:nimage_)
 results_out_out%vel_cell(1:3,1:3,1:nimage_)  =results_out_in%vel_cell(1:3,1:3,1:nimage_)

end subroutine copy_results_out
!!***

!----------------------------------------------------------------------

!!****f* m_results_out/gather_results_out
!! NAME
!!  gather_results_out
!!
!! FUNCTION
!!  Gather results_out datastructure array using communicator over images (replicas) of the cell.
!!  Each contribution of single processor is gathered into a big array on master processor
!!
!! INPUTS
!!  allgather= --optional, default=false--  if TRUE do ALL_GATHER instead of GATHER
!!  dtsets(:)= <type datafiles_type> contains all input variables,
!!  master= --optional, default=0-- index of master proc receiving gathered data (if allgather=false)
!!  mpi_enregs=information about MPI parallelization
!!  only_one_per_img= --optional, default=true--  if TRUE, the gather operation
!!                    is only done by one proc per image (master of the comm_cell)
!!  results_out(:)=<type(results_out_type)>=results_out datastructure array on each proc
!!  use_results_all=true if results_out_all datastructure is allocated for current proc
!!
!! SIDE EFFECTS
!!  === f use_results_all=true ===
!!  results_out_all(:)=<type(results_out_type)>=global (gathered) results_out datastructure array
!!
!! PARENTS
!!      abinit,driver
!!
!! CHILDREN
!!      copy_results_out,init_results_out,xmpi_allgather,xmpi_allgatherv
!!      xmpi_gatherv
!!
!! SOURCE

subroutine gather_results_out(dtsets,mpi_enregs,results_out,results_out_all,use_results_all,&
&                             master,allgather,only_one_per_img) ! optional arguments

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: master
 logical,optional,intent(in) :: allgather,only_one_per_img
 logical,intent(in) :: use_results_all
!arrays
 type(dataset_type),intent(in) :: dtsets(:)
 type(results_out_type),intent(in) :: results_out(:)
 type(results_out_type),intent(inout) :: results_out_all(:)
 type(MPI_type), intent(inout) :: mpi_enregs(:)
!Local variables-------------------------------
!scalars
 integer :: dtsets_size
 integer :: ibufi,ibufr
 integer :: idt1,idt2,ierr,ii,iproc,jj
 integer :: isize,isize_img
 integer :: master_all,master_img,master_one_img
 integer :: mpi_enregs_size,mxnatom,mxnband,mxnkpt,mxnpsp,mxnsppol,mxntypat
 integer :: natom_,nkpt_,nocc_,npsp_,ntypat_,nimage,nimagetot
 integer :: results_out_size,results_out_all_size
 integer :: rsize,rsize_img
 logical :: do_allgather,one_per_img
 character(len=500) :: msg
! type(MPI_type):: mpi_img
!arrays
 integer,allocatable :: ibuffer(:),ibuffer_all(:),ibufshft(:)
 integer,allocatable :: iimg(:),isize_img_all(:),nimage_all(:)
 integer,allocatable :: rbufshft(:),rsize_img_all(:)
 real(dp),allocatable :: rbuffer(:),rbuffer_all(:)

!************************************************************************

 !@results_out_type

 one_per_img=.true.;if (present(only_one_per_img)) one_per_img=only_one_per_img
 do_allgather=.false.;if (present(allgather)) do_allgather=allgather
 master_all=0;if (present(master)) master_all=master

! call init_mpi_enreg(mpi_img,init_mpi=.false.)
 master_img=0;master_one_img=0
! i_am_master=(mpi_img%me==master_all)
! use_results_all= &
!&  (((     do_allgather).and.(     one_per_img).and.(mpi_img%me_cell==master_one_img)) .or. &
!&   ((     do_allgather).and.(.not.one_per_img))                                          .or. &
!&   ((.not.do_allgather).and.(     one_per_img).and.(mpi_img%me==master_all))             .or. &
!&   ((.not.do_allgather).and.(.not.one_per_img).and.(mpi_img%me_img==master_img)))

 dtsets_size=size(dtsets);results_out_size=size(results_out)
 mpi_enregs_size=size(mpi_enregs)
 if (dtsets_size/=results_out_size) then
   msg='  Wrong sizes for dtsets and results_out datastructures !'
   MSG_BUG(msg)
 end if
 if (mpi_enregs_size/=results_out_size) then
   msg='  Wrong sizes for dtsets and results_out datastructures !'
   MSG_BUG(msg)
 end if

 if (use_results_all) then
   results_out_all_size=size(results_out_all)
   if (results_out_size/=results_out_all_size) then
     msg='  Wrong size for results_out_all datastructure !'
     MSG_BUG(msg)
   end if
 end if

 if (results_out_size>0) then

   idt1=lbound(results_out,1);idt2=ubound(results_out,1)

!  Create global results_out_all datastructure
   if (use_results_all) then
     mxnatom=1;mxnband=1;mxnkpt=1;mxnpsp=1;mxntypat=1
     do ii=idt1,idt2
       isize=size(results_out(ii)%fcart,2) ;if (isize>mxnatom) mxnatom=isize
       isize=size(results_out(ii)%occ,1)   ;if (isize>mxnband) mxnband=isize
       isize=size(results_out(ii)%mixalch,1);if(isize>mxnpsp) mxnpsp=isize
       isize=size(results_out(ii)%npwtot,1);if (isize>mxnkpt)  mxnkpt=isize
       isize=size(results_out(ii)%mixalch,2);if(isize>mxntypat)  mxntypat=isize
     end do
     mxnband=mxnband/mxnkpt;mxnsppol=1
     call init_results_out(dtsets,2,0,mpi_enregs,mxnatom,mxnband,mxnkpt,mxnpsp,mxnsppol,mxntypat,results_out_all)
   end if

!  Loop over results_out components (datasets)
   do ii=idt1,idt2

!    Simple copy in case of 1 image
     if (dtsets(ii)%npimage<=1) then
       if (use_results_all) then
         call copy_results_out(results_out(ii),results_out_all(ii))
       end if
     else

!      Retrieve MPI information for this dataset

       if ((.not.one_per_img).or.(mpi_enregs(ii)%me_cell==master_one_img)) then

!        Gather number of images treated by each proc
         ABI_ALLOCATE(nimage_all,(mpi_enregs(ii)%nproc_img))
         nimage_all=0
         nimage=results_out(ii)%nimage
         call xmpi_allgather(nimage,nimage_all,mpi_enregs(ii)%comm_img,ierr)
         nimagetot=sum(nimage_all)

!        Copy scalars from distributed results_out to gathered one
         if (use_results_all) then
           results_out_all(ii)%nimage=nimagetot
           results_out_all(ii)%natom =results_out(ii)%natom
           results_out_all(ii)%nkpt  =results_out(ii)%nkpt
           results_out_all(ii)%nocc  =results_out(ii)%nocc
           results_out_all(ii)%npsp  =results_out(ii)%npsp
           results_out_all(ii)%ntypat=results_out(ii)%ntypat
         end if

!        Compute number of integers/reals needed by current
!        results_out structure for current proc
         isize=results_out(ii)%nkpt
         rsize=28+16*results_out(ii)%natom+results_out(ii)%nocc+results_out(ii)%npsp*results_out(ii)%ntypat+results_out(ii)%ntypat
         isize_img=results_out(ii)%nimage*isize
         rsize_img=results_out(ii)%nimage*rsize
         ABI_ALLOCATE(isize_img_all,(mpi_enregs(ii)%nproc_img))
         ABI_ALLOCATE(rsize_img_all,(mpi_enregs(ii)%nproc_img))
         isize_img_all(:)=isize*nimage_all(:)
         rsize_img_all(:)=rsize*nimage_all(:)
         ABI_DEALLOCATE(nimage_all)

!        Compute shifts in buffer arrays for each proc
         ABI_ALLOCATE(ibufshft,(mpi_enregs(ii)%nproc_img))
         ibufshft(1)=0
         ABI_ALLOCATE(rbufshft,(mpi_enregs(ii)%nproc_img))
         rbufshft(1)=0
         do jj=2,mpi_enregs(ii)%nproc_img
           ibufshft(jj)=ibufshft(jj-1)+isize_img_all(jj-1)
           rbufshft(jj)=rbufshft(jj-1)+rsize_img_all(jj-1)
         end do

!        Load buffers
         ABI_ALLOCATE(ibuffer,(isize_img))
         ABI_ALLOCATE(rbuffer,(rsize_img))
         ibufi=0;ibufr=0
         natom_=results_out(ii)%natom
         nkpt_ =results_out(ii)%nkpt
         nocc_ =results_out(ii)%nocc
         npsp_ =results_out(ii)%npsp
         ntypat_ =results_out(ii)%ntypat
         do jj=1,results_out(ii)%nimage
           ibuffer(ibufi+1:ibufi+nkpt_)=results_out(ii)%npwtot(1:nkpt_,jj)
           ibufi=ibufi+nkpt_
           rbuffer(ibufr+1:ibufr+3)=results_out(ii)%acell(1:3,jj)
           ibufr=ibufr+3
           rbuffer(ibufr+1:ibufr+ntypat_)=results_out(ii)%amu(1:ntypat_,jj)
           ibufr=ibufr+ntypat_
           rbuffer(ibufr+1)=results_out(ii)%etotal(jj)
           ibufr=ibufr+1
           rbuffer(ibufr+1:ibufr+3*natom_)=reshape(results_out(ii)%fcart(1:3,1:natom_,jj),(/3*natom_/))
           ibufr=ibufr+3*natom_
           rbuffer(ibufr+1:ibufr+3*natom_)=reshape(results_out(ii)%fred(1:3,1:natom_,jj),(/3*natom_/))
           ibufr=ibufr+3*natom_
           rbuffer(ibufr+1:ibufr+4*natom_)=reshape(results_out(ii)%intgres(1:4,1:natom_,jj),(/4*natom_/))
           ibufr=ibufr+4*natom_
           rbuffer(ibufr+1:ibufr+npsp_*ntypat_)=&
&               reshape(results_out(ii)%mixalch(1:npsp_,1:ntypat_,jj),(/npsp_*ntypat_/) )
           ibufr=ibufr+npsp_*ntypat_
           rbuffer(ibufr+1:ibufr+nocc_)=results_out(ii)%occ(1:nocc_,jj)
           ibufr=ibufr+nocc_
           rbuffer(ibufr+1:ibufr+9)=reshape(results_out(ii)%rprim(1:3,1:3,jj),(/9/))
           ibufr=ibufr+9
           rbuffer(ibufr+1:ibufr+9)=reshape(results_out(ii)%vel_cell(1:3,1:3,jj),(/9/))
           ibufr=ibufr+9
           rbuffer(ibufr+1:ibufr+6)=results_out(ii)%strten(1:6,jj)
           ibufr=ibufr+6
           rbuffer(ibufr+1:ibufr+3*natom_)=reshape(results_out(ii)%vel(1:3,1:natom_,jj),(/3*natom_/))
           ibufr=ibufr+3*natom_
           rbuffer(ibufr+1:ibufr+3*natom_)=reshape(results_out(ii)%xred(1:3,1:natom_,jj),(/3*natom_/))
           ibufr=ibufr+3*natom_
         end do
         if (ibufi/=isize_img.or.ibufr/=rsize_img) then
           msg='  wrong buffer sizes !'
           MSG_BUG(msg)
         end if

!        Gather all data
         if (use_results_all)  then
           ABI_ALLOCATE(ibuffer_all,(isize*nimagetot))
           ABI_ALLOCATE(rbuffer_all,(rsize*nimagetot))
         end if
         if (.not.use_results_all)  then
           ABI_ALLOCATE(ibuffer_all,(0))
           ABI_ALLOCATE(rbuffer_all,(0))
         end if
         if (do_allgather) then
           call xmpi_allgatherv(ibuffer,isize_img,ibuffer_all,isize_img_all,ibufshft,&
&                               mpi_enregs(ii)%comm_img,ierr)
           call xmpi_allgatherv(rbuffer,rsize_img,rbuffer_all,rsize_img_all,rbufshft,&
&                               mpi_enregs(ii)%comm_img,ierr)
         else
           call xmpi_gatherv(ibuffer,isize_img,ibuffer_all,isize_img_all,ibufshft,&
&                            master_img,mpi_enregs(ii)%comm_img,ierr)
           call xmpi_gatherv(rbuffer,rsize_img,rbuffer_all,rsize_img_all,rbufshft,&
&                            master_img,mpi_enregs(ii)%comm_img,ierr)
         end if
         ABI_DEALLOCATE(isize_img_all)
         ABI_DEALLOCATE(rsize_img_all)
         ABI_DEALLOCATE(ibuffer)
         ABI_DEALLOCATE(rbuffer)

!        Transfer buffers into gathered results_out_all (master proc only)
         if (use_results_all) then
           ABI_ALLOCATE(iimg,(mpi_enregs(ii)%nproc_img))
           iimg=0
           natom_=results_out_all(ii)%natom
           nkpt_=results_out_all(ii)%nkpt
           nocc_=results_out_all(ii)%nocc
           npsp_ =results_out_all(ii)%npsp
           ntypat_ =results_out_all(ii)%ntypat
           do jj=1,nimagetot
!            The following line supposes that images are sorted by increasing index
             iproc=mpi_enregs(ii)%distrb_img(jj)+1;iimg(iproc)=iimg(iproc)+1
             ibufi=ibufshft(iproc)+(iimg(iproc)-1)*isize
             ibufr=rbufshft(iproc)+(iimg(iproc)-1)*rsize
             results_out_all(ii)%npwtot(1:nkpt_,jj)=ibuffer_all(ibufi+1:ibufi+nkpt_)
             ibufi=ibufi+nkpt_
             results_out_all(ii)%acell(1:3,jj)=rbuffer_all(ibufr+1:ibufr+3)
             ibufr=ibufr+3
             results_out_all(ii)%amu(1:ntypat_,jj)=rbuffer_all(ibufr+1:ibufr+ntypat_)
             ibufr=ibufr+ntypat_
             results_out_all(ii)%etotal(jj)=rbuffer_all(ibufr+1)
             ibufr=ibufr+1
             results_out_all(ii)%fcart(1:3,1:natom_,jj)= &
&                   reshape(rbuffer_all(ibufr+1:ibufr+3*natom_),(/3,natom_/))
             ibufr=ibufr+3*natom_
             results_out_all(ii)%fred(1:3,1:natom_,jj)= &
&                   reshape(rbuffer_all(ibufr+1:ibufr+3*natom_),(/3,natom_/))
             ibufr=ibufr+3*natom_
             results_out_all(ii)%intgres(1:4,1:natom_,jj)= &
&                   reshape(rbuffer_all(ibufr+1:ibufr+4*natom_),(/4,natom_/))
             ibufr=ibufr+4*natom_
             results_out_all(ii)%mixalch(1:npsp_,1:ntypat_,jj)= &
&                   reshape(rbuffer_all(ibufr+1:ibufr+npsp_*ntypat_),(/npsp_,ntypat_/))
             ibufr=ibufr+npsp_*ntypat_
             results_out_all(ii)%occ(1:nocc_,jj)=rbuffer_all(ibufr+1:ibufr+nocc_)
             ibufr=ibufr+nocc_
             results_out_all(ii)%rprim(1:3,1:3,jj)=reshape(rbuffer_all(ibufr+1:ibufr+9),(/3,3/))
             ibufr=ibufr+9
             results_out_all(ii)%vel_cell(1:3,1:3,jj)=reshape(rbuffer_all(ibufr+1:ibufr+9),(/3,3/))
             ibufr=ibufr+9
             results_out_all(ii)%strten(1:6,jj)=rbuffer_all(ibufr+1:ibufr+6)
             ibufr=ibufr+6
             results_out_all(ii)%vel(1:3,1:natom_,jj)= &
&                   reshape(rbuffer_all(ibufr+1:ibufr+3*natom_),(/3,natom_/))
             ibufr=ibufr+3*natom_
             results_out_all(ii)%xred(1:3,1:natom_,jj)= &
&                   reshape(rbuffer_all(ibufr+1:ibufr+3*natom_),(/3,natom_/))
             ibufr=ibufr+3*natom_
           end do
           ABI_DEALLOCATE(iimg)
         end if

!        Free memory
         ABI_DEALLOCATE(ibufshft)
         ABI_DEALLOCATE(rbufshft)
         ABI_DEALLOCATE(ibuffer_all)
         ABI_DEALLOCATE(rbuffer_all)

       end if
     end if
   end do
 end if

end subroutine gather_results_out
!!***

END MODULE m_results_out
!!***
