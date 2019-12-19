!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_results_img
!! NAME
!!  m_results_img
!!
!! FUNCTION
!!  This module provides the definition of the results_img_type used
!!  to store results from GS calculations for a given image of the cell.
!!
!! COPYRIGHT
!! Copyright (C) 2011-2019 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
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

MODULE m_results_img

 use defs_basis
 use m_energies
 use m_abicore
 use m_results_gs
 use m_errors
 use m_xmpi

 use defs_abitypes, only : mpi_type
 use m_geometry,   only : mkrdim, xred2xcart

 implicit none

 private

! public procedures.
 public :: init_results_img
 public :: destroy_results_img
 public :: nullify_results_img
 public :: copy_results_img
 public :: gather_results_img
 public :: gather_array_img
 public :: scatter_array_img
 public :: get_geometry_img

 interface gather_array_img
   module procedure gather_array_img_1D
   module procedure gather_array_img_2D
 end interface gather_array_img
!!***

!!****t* m_results_img/results_img_type
!! NAME
!! results_img_type
!!
!! FUNCTION
!! This structured datatype contains the results of a GS calculation
!! for a given image of the cell:
!!   energy, forces, stresses,positions, velocities, cell parameter, alchemical mixing parameters

!!
!! SOURCE

 type, public :: results_img_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! Integer scalar

  integer :: natom
   ! The number of atoms for this image

  integer :: npspalch
   ! The number of pseudopotentials for alchemical mixing  for this image

  integer :: nspden
   ! The number of spin-density channels for this image

  integer :: nsppol
   ! The number of spin channels for this image

  integer :: ntypat
   ! The number of types of atoms for this image

  integer :: ntypalch
   ! The number of alchemical pseudoatoms  for this image

! Real (real(dp)) arrays

  real(dp), allocatable :: acell(:)
   ! acell(3)
   ! Dimensions of the cell

  real(dp), allocatable :: amu(:)
   ! amu(ntypat)
   ! Atomic masses of the different types of atoms

  real(dp), allocatable :: mixalch(:,:)
   ! mixalch(npspalch,ntypalch)
   ! Alchemical mixing factors

  real(dp), allocatable :: rprim(:,:)
   ! rprim(3,3)
   ! Primitive translations of the cell

  real(dp), allocatable :: vel(:,:)
   ! vel(3,natom)
   ! Velocities of the atoms

  real(dp), allocatable :: vel_cell(:,:)
   ! vel_cell(3,3)
   ! Time derivative of the cell parameters

  real(dp), allocatable :: xred(:,:)
   ! xred(3,natom)
   ! Reduced coordinates of the atoms

! Other types of data
  type(results_gs_type), pointer :: results_gs => null()
   ! Energies, forces and stresses from the GS calculation

 end type results_img_type
!!***

CONTAINS

!===========================================================
!!***

!!****f* m_results_img/init_results_img
!! NAME
!!  init_results_img
!!
!! FUNCTION
!!  Init all scalars and allocatables in an array of results_img datastructures
!!
!! INPUTS
!!  natom=number of atoms in cell
!!  npspalch = number of pseudopotentials for alchemical mixing  for this image
!!  nspden   = number of spin-density channels
!!  nsppol   = number of spin channels
!!  ntypalch = The number of alchemical pseudoatoms  for this image
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  results_img(:)=<type(results_img_type)>=results_img datastructure array
!!
!! PARENTS
!!      gstateimg,m_results_img
!!
!! CHILDREN
!!      mkrdim,xred2xcart
!!
!! SOURCE

subroutine init_results_img(natom,npspalch,nspden,nsppol,ntypalch,ntypat,results_img)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,npspalch,nspden,nsppol,ntypalch,ntypat
!arrays
 type(results_img_type),intent(inout) :: results_img(:)
!Local variables-------------------------------
!scalars
 integer :: ii,results_img_size
!arrays

!************************************************************************

 !@results_img_type

 results_img_size=size(results_img)

 if (results_img_size>0) then

   do ii=1,results_img_size

     results_img(ii)%natom  =natom
     results_img(ii)%npspalch  =npspalch
     results_img(ii)%nspden    =nspden
     results_img(ii)%nsppol    =nsppol
     results_img(ii)%ntypalch  =ntypalch
     results_img(ii)%ntypat    =ntypat

     ABI_DATATYPE_ALLOCATE(results_img(ii)%results_gs,)
     call init_results_gs(natom,nspden,nsppol,results_img(ii)%results_gs)

     ABI_ALLOCATE(results_img(ii)%acell,(3))
     results_img(ii)%acell=zero
     ABI_ALLOCATE(results_img(ii)%amu,(ntypat))
     results_img(ii)%amu=zero
     ABI_ALLOCATE(results_img(ii)%mixalch,(npspalch,ntypalch))
     results_img(ii)%mixalch=zero
     ABI_ALLOCATE(results_img(ii)%rprim,(3,3))
     results_img(ii)%rprim=zero
     ABI_ALLOCATE(results_img(ii)%xred,(3,natom))
     results_img(ii)%xred =zero
     ABI_ALLOCATE(results_img(ii)%vel,(3,natom))
     results_img(ii)%vel  =zero
     ABI_ALLOCATE(results_img(ii)%vel_cell,(3,3))
     results_img(ii)%vel_cell=zero

   end do
 end if

end subroutine init_results_img
!!***

!----------------------------------------------------------------------

!!****f* m_results_img/destroy_results_img
!! NAME
!!  destroy_results_img
!!
!! FUNCTION
!!  Clean and destroy an array of results_img datastructures
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  results_img(:)=<type(results_img_type)>=results_img datastructure array
!!
!! PARENTS
!!      gstateimg,prtimg
!!
!! CHILDREN
!!      mkrdim,xred2xcart
!!
!! SOURCE

subroutine destroy_results_img(results_img)

!Arguments ------------------------------------
!arrays
 type(results_img_type),intent(inout) :: results_img(:)
!Local variables-------------------------------
!scalars
 integer :: ii,results_img_size

!************************************************************************

 !@results_img_type

 results_img_size=size(results_img)
 if (results_img_size>0) then

   do ii=1,results_img_size

     results_img(ii)%natom=0
     results_img(ii)%npspalch=0
     results_img(ii)%nspden  =0
     results_img(ii)%nsppol  =0
     results_img(ii)%ntypalch=0
     results_img(ii)%ntypat=0

     if (associated(results_img(ii)%results_gs)) then
       call destroy_results_gs(results_img(ii)%results_gs)
       ABI_DATATYPE_DEALLOCATE(results_img(ii)%results_gs)
     end if

     if (allocated(results_img(ii)%acell))  then
       ABI_DEALLOCATE(results_img(ii)%acell)
     end if
     if (allocated(results_img(ii)%amu))  then
       ABI_DEALLOCATE(results_img(ii)%amu)
     end if
     if (allocated(results_img(ii)%mixalch))  then
       ABI_DEALLOCATE(results_img(ii)%mixalch)
     end if
     if (allocated(results_img(ii)%rprim))  then
       ABI_DEALLOCATE(results_img(ii)%rprim)
     end if
     if (allocated(results_img(ii)%xred))   then
       ABI_DEALLOCATE(results_img(ii)%xred)
     end if
     if (allocated(results_img(ii)%vel))    then
       ABI_DEALLOCATE(results_img(ii)%vel)
     end if
     if (allocated(results_img(ii)%vel_cell))    then
       ABI_DEALLOCATE(results_img(ii)%vel_cell)
     end if
   end do

 end if

end subroutine destroy_results_img
!!***

!----------------------------------------------------------------------

!!****f* m_results_img/nullify_results_img
!! NAME
!!  nullify_results_img
!!
!! FUNCTION
!!  Nullify an array of results_img datastructures
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  results_img(:)=<type(results_img_type)>=results_img datastructure array
!!
!! PARENTS
!!
!! CHILDREN
!!      mkrdim,xred2xcart
!!
!! SOURCE

subroutine nullify_results_img(results_img)

!Arguments ------------------------------------
!arrays
 type(results_img_type),intent(inout) :: results_img(:)
!Local variables-------------------------------
!scalars
 integer :: ii,results_img_size

!************************************************************************

 !@results_img_type

 results_img_size=size(results_img)
 if (results_img_size>0) then

   do ii=1,results_img_size
     results_img(ii)%natom=0
     results_img(ii)%npspalch=0
     results_img(ii)%nspden=0
     results_img(ii)%nsppol=0
     results_img(ii)%ntypalch=0
     results_img(ii)%ntypat=0
     nullify(results_img(ii)%results_gs)
   end do

 end if

end subroutine nullify_results_img
!!***

!----------------------------------------------------------------------

!!****f* m_results_img/copy_results_img
!! NAME
!!  copy_results_img
!!
!! FUNCTION
!!  Copy a results_img datastructure into another
!!
!! INPUTS
!!  results_img_in=<type(results_img_type)>=input results_img datastructure
!!
!! OUTPUT
!!  results_img_out=<type(results_img_type)>=output results_img datastructure
!!
!! PARENTS
!!      gstateimg,m_results_img
!!
!! CHILDREN
!!      mkrdim,xred2xcart
!!
!! SOURCE

subroutine copy_results_img(results_img_in,results_img_out)

!Arguments ------------------------------------
!arrays
 type(results_img_type),intent(in) :: results_img_in
 type(results_img_type),intent(inout) :: results_img_out !vz_i
!Local variables-------------------------------
!scalars
 integer :: natom_in,natom_out,npspalch_in,npspalch_out,ntypalch_in,ntypalch_out
 integer :: nspden_in, nspden_out, nsppol_in, nsppol_out,ntypat_in,ntypat_out

!************************************************************************

 !@results_img_type

 natom_in =results_img_in%natom
 natom_out=results_img_out%natom
 npspalch_in =results_img_in%npspalch
 npspalch_out=results_img_out%npspalch
 nspden_in =results_img_in%nspden
 nspden_out=results_img_out%nspden
 nsppol_in =results_img_in%nsppol
 nsppol_out=results_img_out%nsppol
 ntypalch_in =results_img_in%ntypalch
 ntypalch_out=results_img_out%ntypalch
 ntypat_in =results_img_in%ntypat
 ntypat_out=results_img_out%ntypat

 if (natom_in>natom_out) then
   if (allocated(results_img_out%xred))   then
     ABI_DEALLOCATE(results_img_out%xred)
   end if
   if (allocated(results_img_out%vel))    then
     ABI_DEALLOCATE(results_img_out%vel)
   end if
   if (allocated(results_img_out%vel_cell))    then
     ABI_DEALLOCATE(results_img_out%vel_cell)
   end if

   if (allocated(results_img_in%xred))   then
     ABI_ALLOCATE(results_img_out%xred,(3,natom_in))
   end if
   if (allocated(results_img_in%vel))    then
     ABI_ALLOCATE(results_img_out%vel,(3,natom_in))
   end if
   if (allocated(results_img_in%vel_cell))    then
     ABI_ALLOCATE(results_img_out%vel_cell,(3,3))
   end if
 end if

 if (npspalch_in>npspalch_out .or. ntypalch_in>ntypalch_out) then
   if (allocated(results_img_out%mixalch))   then
     ABI_DEALLOCATE(results_img_out%mixalch)
   end if
 endif

 if (ntypat_in>ntypat_out) then
   if (allocated(results_img_in%amu))    then
     ABI_ALLOCATE(results_img_out%amu,(ntypat_in))
   endif
 endif

 results_img_out%natom  =results_img_in%natom
 results_img_out%npspalch  =results_img_in%npspalch
 results_img_out%nspden =results_img_in%nspden
 results_img_out%nsppol =results_img_in%nsppol
 results_img_out%ntypalch  =results_img_in%ntypalch
 results_img_out%ntypat  =results_img_in%ntypat

 call copy_results_gs(results_img_in%results_gs,results_img_out%results_gs)

 if (allocated(results_img_in%acell)) results_img_out%acell(:)=results_img_in%acell(:)
 if (allocated(results_img_in%amu))   results_img_out%amu(1:ntypat_in)=results_img_in%amu(1:ntypat_in)
 if (allocated(results_img_in%mixalch)) &
&   results_img_out%mixalch(1:npspalch_in,1:ntypalch_in)=results_img_in%mixalch(1:npspalch_in,1:ntypalch_in)
 if (allocated(results_img_in%rprim))   results_img_out%rprim(:,:)=results_img_in%rprim(:,:)
 if (allocated(results_img_in%xred))    results_img_out%xred(:,1:natom_in)=results_img_in%xred(:,1:natom_in)
 if (allocated(results_img_in%vel))     results_img_out%vel(:,1:natom_in)=results_img_in%vel(:,1:natom_in)
 if (allocated(results_img_in%vel_cell))results_img_out%vel_cell(:,:)=results_img_in%vel_cell(:,:)

end subroutine copy_results_img
!!***

!----------------------------------------------------------------------

!!****f* m_results_img/gather_results_img
!! NAME
!!  gather_results_img
!!
!! FUNCTION
!!  Gather results_img datastructures using communicator over images (replicas) of the cell.
!!  Each contribution of single processor is gathered into a big array on master processor
!!
!! INPUTS
!!  allgather= --optional, default=false--  if TRUE do ALL_GATHER instead of GATHER
!!  master= --optional, default=0-- index of master proc receiving gathered data (if allgather=false)
!!  mpi_enreg=information about MPI parallelization
!!  only_one_per_img= --optional, default=true--  if TRUE, the gather operation
!!                    is only done by one proc per image (master of the comm_cell)
!!  results_img(:)=<type(results_img_type)>=results_img datastructure array on each proc
!!
!! SIDE EFFECTS
!!  results_img_all(:)=<type(results_img_type)>=global (gathered) results_img datastructure array
!!
!! PARENTS
!!      prtimg
!!
!! CHILDREN
!!      mkrdim,xred2xcart
!!
!! SOURCE

subroutine gather_results_img(mpi_enreg,results_img,results_img_all,&
&                 master,allgather,only_one_per_img) ! optional arguments

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: master
 logical,optional,intent(in) :: allgather,only_one_per_img
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 type(results_img_type),intent(inout) :: results_img(:)
 type(results_img_type),intent(inout) :: results_img_all(:)

!Local variables-------------------------------
!scalars
 integer :: ibufr,ierr,iproc,jj
 integer :: master_all,master_img,master_one_img
 integer :: natom,ngrvdw,nimage,nimagetot,npspalch,nspden,nsppol,ntypalch,ntypat
 integer :: rsize,rsize_img
 logical :: do_allgather,i_am_master,one_per_img,use_results_all
 !character(len=500) :: msg
!arrays
 integer,allocatable :: iimg(:),nimage_all(:),rbufshft(:),rsize_img_all(:)
 real(dp),allocatable :: rbuffer(:),rbuffer_all(:)

!************************************************************************

 !@results_img_type

 one_per_img=.true.;if (present(only_one_per_img)) one_per_img=only_one_per_img
 do_allgather=.false.;if (present(allgather)) do_allgather=allgather
 master_all=0;if (present(master)) master_all=master
 i_am_master=(mpi_enreg%me==master_all)

 master_img=0;master_one_img=0
 use_results_all= &
&  (((     do_allgather).and.(     one_per_img).and.(mpi_enreg%me_cell==master_one_img)) .or. &
&   ((     do_allgather).and.(.not.one_per_img))                                            .or. &
&   ((.not.do_allgather).and.(     one_per_img).and.(mpi_enreg%me==master_all))             .or. &
&   ((.not.do_allgather).and.(.not.one_per_img).and.(mpi_enreg%me_img==master_img)))

!Create global results_img_all datastructure
 if (use_results_all) then
   call init_results_img(results_img(1)%natom,results_img(1)%npspalch,results_img(1)%nspden,results_img(1)%nsppol,&
&    results_img(1)%ntypalch,results_img(1)%ntypat,results_img_all)
 end if

 if ((.not.one_per_img).or.(mpi_enreg%me_cell==master_one_img)) then

!  Simple copy in case of 1 image
   if (use_results_all) then
     if (size(results_img_all,1)<=1) then
       call copy_results_img(results_img(1),results_img_all(1))
       return
     end if
   endif

!  Gather number of images treated by each proc
   nimage=size(results_img,1)
   ABI_ALLOCATE(nimage_all,(mpi_enreg%nproc_img))
   call xmpi_allgather(nimage,nimage_all,mpi_enreg%comm_img,ierr)
   nimagetot=sum(nimage_all)
   if (use_results_all) then
     if (size(results_img_all,1)/=nimagetot) then
       MSG_BUG('Wrong results_img_all size !')
     end if
   end if

!  Copy natom from distributed results_img to gathered one
   natom =results_img(1)%natom
   npspalch =results_img(1)%npspalch
   nspden =results_img(1)%nspden
   nsppol =results_img(1)%nsppol
   ntypalch =results_img(1)%ntypalch
   ntypat   =results_img(1)%ntypat
   ngrvdw=results_img(1)%results_gs%ngrvdw
   if (use_results_all) then
     do jj=1,nimagetot
       results_img_all(jj)%natom =natom
       results_img_all(jj)%npspalch =npspalch
       results_img_all(jj)%nspden =nspden
       results_img_all(jj)%nsppol =nsppol
       results_img_all(jj)%ntypalch =ntypalch
       results_img_all(jj)%ntypat   =ntypat
       results_img_all(jj)%results_gs%ngrvdw=ngrvdw
     enddo
   end if

!  Compute number of data
   rsize=38+n_energies+(30+nspden)*natom+3*nsppol+ntypat+npspalch*ntypalch+3*ngrvdw
   rsize_img=nimage*rsize
   ABI_ALLOCATE(rsize_img_all,(mpi_enreg%nproc_img))
   rsize_img_all(:)=rsize*nimage_all(:)
   ABI_DEALLOCATE(nimage_all)

!  Compute shifts in buffer arrays for each proc
   ABI_ALLOCATE(rbufshft,(mpi_enreg%nproc_img))
   rbufshft(1)=0
   do jj=2,mpi_enreg%nproc_img
     rbufshft(jj)=rbufshft(jj-1)+rsize_img_all(jj-1)
   end do

!  Load buffers
   ABI_ALLOCATE(rbuffer,(rsize_img))
   ibufr=0
   do jj=1,nimage
     rbuffer(ibufr+1)  =results_img(jj)%results_gs%deltae
     rbuffer(ibufr+2)  =results_img(jj)%results_gs%diffor
     rbuffer(ibufr+3)  =results_img(jj)%results_gs%entropy
     rbuffer(ibufr+4)  =results_img(jj)%results_gs%etotal
     rbuffer(ibufr+5)  =results_img(jj)%results_gs%fermie
     rbuffer(ibufr+6)  =results_img(jj)%results_gs%residm
     rbuffer(ibufr+7)  =results_img(jj)%results_gs%res2
     rbuffer(ibufr+8)  =results_img(jj)%results_gs%vxcavg
     rbuffer(ibufr+9:ibufr+11) =results_img(jj)%results_gs%pel(1:3)
     rbuffer(ibufr+12:ibufr+17)=results_img(jj)%results_gs%strten(1:6)
     rbuffer(ibufr+18:ibufr+20)=results_img(jj)%acell(1:3)
     rbuffer(ibufr+21:ibufr+23)=results_img(jj)%rprim(1:3,1)
     rbuffer(ibufr+24:ibufr+26)=results_img(jj)%rprim(1:3,2)
     rbuffer(ibufr+27:ibufr+29)=results_img(jj)%rprim(1:3,3)
     rbuffer(ibufr+30:ibufr+32)=results_img(jj)%vel_cell(1:3,1)
     rbuffer(ibufr+33:ibufr+35)=results_img(jj)%vel_cell(1:3,2)
     rbuffer(ibufr+36:ibufr+38)=results_img(jj)%vel_cell(1:3,3)
     ibufr=ibufr+38
     call energies_to_array(results_img(jj)%results_gs%energies,&
&                           rbuffer(ibufr+1:ibufr+n_energies),1)
     ibufr=ibufr+n_energies
     rbuffer(ibufr+1:ibufr+3*natom)=reshape(results_img(jj)%results_gs%fcart(1:3,1:natom),(/3*natom/))
     ibufr=ibufr+3*natom
     rbuffer(ibufr+1:ibufr+3*natom)=reshape(results_img(jj)%results_gs%fred(1:3,1:natom),(/3*natom/))
     ibufr=ibufr+3*natom
     rbuffer(ibufr+1:ibufr+3*nsppol)=reshape(results_img(jj)%results_gs%gaps(1:3,1:nsppol),(/3*nsppol/))
     ibufr=ibufr+3*nsppol
     rbuffer(ibufr+1:ibufr+3*natom)=reshape(results_img(jj)%results_gs%grchempottn(1:3,1:natom),(/3*natom/))
     ibufr=ibufr+3*natom
     rbuffer(ibufr+1:ibufr+3*natom)=reshape(results_img(jj)%results_gs%grcondft(1:3,1:natom),(/3*natom/))
     ibufr=ibufr+3*natom
     rbuffer(ibufr+1:ibufr+3*natom)=reshape(results_img(jj)%results_gs%gresid(1:3,1:natom),(/3*natom/))
     ibufr=ibufr+3*natom
     rbuffer(ibufr+1:ibufr+3*natom)=reshape(results_img(jj)%results_gs%grewtn(1:3,1:natom),(/3*natom/))
     ibufr=ibufr+3*natom
     rbuffer(ibufr+1:ibufr+3*natom)=reshape(results_img(jj)%results_gs%grxc(1:3,1:natom),(/3*natom/))
     ibufr=ibufr+3*natom
     rbuffer(ibufr+1:ibufr+nspden*natom)=reshape(results_img(jj)%results_gs%intgres(1:nspden,1:natom),(/nspden*natom/))
     ibufr=ibufr+nspden*natom
     rbuffer(ibufr+1:ibufr+3*natom)=reshape(results_img(jj)%results_gs%synlgr(1:3,1:natom),(/3*natom/))
     ibufr=ibufr+3*natom
     rbuffer(ibufr+1:ibufr+ntypat)=results_img(jj)%amu(1:ntypat)
     ibufr=ibufr+ntypat
     rbuffer(ibufr+1:ibufr+npspalch*ntypalch)=reshape(results_img(jj)%mixalch(1:npspalch,1:ntypalch),(/npspalch*ntypalch/))
     ibufr=ibufr+npspalch*ntypalch
     rbuffer(ibufr+1:ibufr+3*natom)=reshape(results_img(jj)%xred(1:3,1:natom),(/3*natom/))
     ibufr=ibufr+3*natom
     rbuffer(ibufr+1:ibufr+3*natom)=reshape(results_img(jj)%vel(1:3,1:natom),(/3*natom/))
     ibufr=ibufr+3*natom
     if (ngrvdw>0) then
       rbuffer(ibufr+1:ibufr+3*ngrvdw)= &
&        reshape(results_img(jj)%results_gs%grvdw(1:3,1:ngrvdw),(/3*ngrvdw/))
       ibufr=ibufr+3*ngrvdw
     end if
   end do
   if (ibufr/=rsize_img) then
     MSG_BUG('wrong buffer size !')
   end if

!  Gather all data
   if (use_results_all)  then
     ABI_ALLOCATE(rbuffer_all,(rsize*nimagetot))
   end if
   if (.not.use_results_all)  then
     ABI_ALLOCATE(rbuffer_all,(0))
   end if
   if (do_allgather) then
     call xmpi_allgatherv(rbuffer,rsize_img,rbuffer_all,rsize_img_all,rbufshft,&
&                         mpi_enreg%comm_img,ierr)
   else
     call xmpi_gatherv(rbuffer,rsize_img,rbuffer_all,rsize_img_all,rbufshft,&
&                          master_img,mpi_enreg%comm_img,ierr)
   end if
   ABI_DEALLOCATE(rbuffer)
   ABI_DEALLOCATE(rsize_img_all)

!  Transfer buffers into gathered results_img_all (master proc only)
   if (use_results_all) then
     ABI_ALLOCATE(iimg,(mpi_enreg%nproc_img))
     iimg=0
     do jj=1,nimagetot
!      The following line supposes that images are sorted by increasing index
       iproc=mpi_enreg%distrb_img(jj)+1;iimg(iproc)=iimg(iproc)+1
       ibufr=rbufshft(iproc)+(iimg(iproc)-1)*rsize
       results_img_all(jj)%results_gs%deltae     =rbuffer_all(ibufr+1)
       results_img_all(jj)%results_gs%diffor     =rbuffer_all(ibufr+2)
       results_img_all(jj)%results_gs%entropy    =rbuffer_all(ibufr+3)
       results_img_all(jj)%results_gs%etotal     =rbuffer_all(ibufr+4)
       results_img_all(jj)%results_gs%fermie     =rbuffer_all(ibufr+5)
       results_img_all(jj)%results_gs%residm     =rbuffer_all(ibufr+6)
       results_img_all(jj)%results_gs%res2       =rbuffer_all(ibufr+7)
       results_img_all(jj)%results_gs%vxcavg     =rbuffer_all(ibufr+8)
       results_img_all(jj)%results_gs%pel(1:3)   =rbuffer_all(ibufr+9:ibufr+11)
       results_img_all(jj)%results_gs%strten(1:6)=rbuffer_all(ibufr+12:ibufr+17)
       results_img_all(jj)%acell(1:3)   =rbuffer_all(ibufr+18:ibufr+20)
       results_img_all(jj)%rprim(1:3,1)=rbuffer_all(ibufr+21:ibufr+23)
       results_img_all(jj)%rprim(1:3,2)=rbuffer_all(ibufr+24:ibufr+26)
       results_img_all(jj)%rprim(1:3,3)=rbuffer_all(ibufr+27:ibufr+29)
       results_img_all(jj)%vel_cell(1:3,1)=rbuffer_all(ibufr+30:ibufr+32)
       results_img_all(jj)%vel_cell(1:3,2)=rbuffer_all(ibufr+33:ibufr+35)
       results_img_all(jj)%vel_cell(1:3,3)=rbuffer_all(ibufr+36:ibufr+38)
       ibufr=ibufr+38
       call energies_to_array(results_img_all(jj)%results_gs%energies,&
&                             rbuffer_all(ibufr+1:ibufr+n_energies),-1)
       ibufr=ibufr+n_energies
       results_img_all(jj)%amu(1:ntypat)=rbuffer_all(ibufr+1:ibufr+ntypat)
       results_img_all(jj)%results_gs%fcart(1:3,1:natom)= &
&             reshape(rbuffer_all(ibufr+1:ibufr+3*natom),(/3,natom/))
       ibufr=ibufr+3*natom
       results_img_all(jj)%results_gs%fred(1:3,1:natom)= &
&             reshape(rbuffer_all(ibufr+1:ibufr+3*natom),(/3,natom/))
       ibufr=ibufr+3*natom
       results_img_all(jj)%results_gs%gaps(1:3,1:nsppol)= &
&             reshape(rbuffer_all(ibufr+1:ibufr+3*nsppol),(/3,nsppol/))
       ibufr=ibufr+3*nsppol
       results_img_all(jj)%results_gs%grchempottn(1:3,1:natom)= &
  &            reshape(rbuffer_all(ibufr+1:ibufr+3*natom),(/3,natom/))
       ibufr=ibufr+3*natom
       results_img_all(jj)%results_gs%grcondft(1:3,1:natom)= &
  &            reshape(rbuffer_all(ibufr+1:ibufr+3*natom),(/3,natom/))
       ibufr=ibufr+3*natom
       results_img_all(jj)%results_gs%gresid(1:3,1:natom)= &
  &            reshape(rbuffer_all(ibufr+1:ibufr+3*natom),(/3,natom/))
       ibufr=ibufr+3*natom
       results_img_all(jj)%results_gs%grewtn(1:3,1:natom)= &
  &            reshape(rbuffer_all(ibufr+1:ibufr+3*natom),(/3,natom/))
       ibufr=ibufr+3*natom
       results_img_all(jj)%results_gs%grxc(1:3,1:natom)= &
  &            reshape(rbuffer_all(ibufr+1:ibufr+3*natom),(/3,natom/))
       ibufr=ibufr+3*natom
       results_img_all(jj)%results_gs%intgres(1:nspden,1:natom)= &
  &            reshape(rbuffer_all(ibufr+1:ibufr+nspden*natom),(/nspden,natom/))
       ibufr=ibufr+nspden*natom
       results_img_all(jj)%results_gs%synlgr(1:3,1:natom)= &
  &            reshape(rbuffer_all(ibufr+1:ibufr+3*natom),(/3,natom/))
       ibufr=ibufr+3*natom
       results_img_all(jj)%amu(1:ntypat)=rbuffer_all(ibufr+1:ibufr+ntypat)
       ibufr=ibufr+ntypat
       results_img_all(jj)%mixalch(1:npspalch,1:ntypalch)= &
  &            reshape(rbuffer_all(ibufr+1:ibufr+npspalch*ntypalch),(/npspalch,ntypalch/))
       ibufr=ibufr+npspalch*ntypalch
       results_img_all(jj)%xred(1:3,1:natom)= &
  &            reshape(rbuffer_all(ibufr+1:ibufr+3*natom),(/3,natom/))
       ibufr=ibufr+3*natom
       results_img_all(jj)%vel(1:3,1:natom)= &
  &            reshape(rbuffer_all(ibufr+1:ibufr+3*natom),(/3,natom/))
       ibufr=ibufr+3*natom
       if (ngrvdw>0) then
         results_img_all(jj)%results_gs%grvdw(1:3,1:ngrvdw)= &
  &              reshape(rbuffer_all(ibufr+1:ibufr+3*ngrvdw),(/3,ngrvdw/))
         ibufr=ibufr+3*ngrvdw
       end if
     end do
     ABI_DEALLOCATE(iimg)
   end if

!  Free memory
   ABI_DEALLOCATE(rbufshft)
   ABI_DEALLOCATE(rbuffer_all)

 end if

end subroutine gather_results_img
!!***

!----------------------------------------------------------------------

!!****f* m_results_img/gather_array_img_1D
!! NAME
!!  gather_array_img_1D
!!
!! FUNCTION
!!  Gather an real 1D-array (part of a results_img datastructure) using communicator
!!  over images (replicas) of the cell.
!!  Each contribution of single processor is gathered into a big array on master processor
!!
!! INPUTS
!!  allgather= --optional, default=false--  if TRUE do ALL_GATHER instead of GATHER
!!  master= --optional, default=0-- index of master proc receiving gathered data (if allgather=false)
!!  mpi_enreg=information about MPI parallelization
!!  only_one_per_img= --optional, default=true--  if TRUE, the gather operation
!!                    is only done by one proc per image (master of the comm_cell)
!!  array_img(:,:)= (real) 1D-array distributed (has 2 dimensions; the 2nd one is nimage)
!!
!! SIDE EFFECTS
!!  array_img_all(:,:)= (real) global (gathered) 1D-array
!!                      (has 2 dimensions; the 2nd one is nimagetot)
!!
!! PARENTS
!!
!! CHILDREN
!!      mkrdim,xred2xcart
!!
!! SOURCE

subroutine gather_array_img_1D(array_img,array_img_all,mpi_enreg,&
&                              master,allgather,only_one_per_img) ! optional arguments

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: master
 logical,optional,intent(in) :: allgather,only_one_per_img
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 real(dp),intent(in) :: array_img(:,:)
 real(dp),intent(inout) :: array_img_all(:,:)

!Local variables-------------------------------
!scalars
 integer :: ibufr,ierr,iproc,jj
 integer :: master_all,master_img,master_one_img,nimage,nimagetot
 integer :: rsize,rsize_img,size1
 logical :: do_allgather,i_am_master,one_per_img,use_array_all
 !character(len=500) :: msg
!arrays
 integer,allocatable :: iimg(:),nimage_all(:),rbufshft(:),rsize_img_all(:)
 real(dp),allocatable :: rbuffer(:),rbuffer_all(:)

!************************************************************************

 !@results_img_type

 one_per_img=.true.;if (present(only_one_per_img)) one_per_img=only_one_per_img
 do_allgather=.false.;if (present(allgather)) do_allgather=allgather
 master_all=0;if (present(master)) master_all=master
 i_am_master=(mpi_enreg%me==master_all)

 master_img=0;master_one_img=0
 use_array_all= &
&  (((     do_allgather).and.(     one_per_img).and.(mpi_enreg%me_cell==master_one_img)) .or. &
&   ((     do_allgather).and.(.not.one_per_img))                                            .or. &
&   ((.not.do_allgather).and.(     one_per_img).and.(mpi_enreg%me==master_all))             .or. &
&   ((.not.do_allgather).and.(.not.one_per_img).and.(mpi_enreg%me_img==master_img)))

 size1=size(array_img,1)
 if (use_array_all) then
   if (size(array_img_all,1)/=size1) then
     MSG_BUG('Wrong array_img_all size (1)')
   end if
 end if

 if ((.not.one_per_img).or.(mpi_enreg%me_cell==master_one_img)) then

!  Simple copy in case of 1 image
   if (use_array_all) then
     if (size(array_img_all,2)<=1) then
       array_img_all(:,1)=array_img(:,1)
       return
     end if
   endif

!  Gather number of images treated by each proc
   nimage=size(array_img,2)
   ABI_ALLOCATE(nimage_all,(mpi_enreg%nproc_img))
   call xmpi_allgather(nimage,nimage_all,mpi_enreg%comm_img,ierr)
   nimagetot=sum(nimage_all)
   if (use_array_all) then
     if (size(array_img_all,2)/=nimagetot) then
       MSG_BUG('Wrong array_img_all size (2)!')
     endif
   end if

!  Compute number of data
   rsize=size1;rsize_img=nimage*rsize
   ABI_ALLOCATE(rsize_img_all,(mpi_enreg%nproc_img))
   rsize_img_all(:)=rsize*nimage_all(:)
   ABI_DEALLOCATE(nimage_all)

!  Compute shifts in buffer arrays for each proc
   ABI_ALLOCATE(rbufshft,(mpi_enreg%nproc_img))
   rbufshft(1)=0
   do jj=2,mpi_enreg%nproc_img
     rbufshft(jj)=rbufshft(jj-1)+rsize_img_all(jj-1)
   end do

!  Load buffers
   ABI_ALLOCATE(rbuffer,(rsize_img))
   ibufr=0
   do jj=1,nimage
     rbuffer(ibufr+1:ibufr+rsize)=reshape(array_img(1:size1,jj),(/rsize/))
     ibufr=ibufr+rsize
   end do

!  Gather all data
   if (use_array_all)  then
     ABI_ALLOCATE(rbuffer_all,(rsize*nimagetot))
   end if
   if (.not.use_array_all)  then
     ABI_ALLOCATE(rbuffer_all,(0))
   end if
   if (do_allgather) then
     call xmpi_allgatherv(rbuffer,rsize_img,rbuffer_all,rsize_img_all,rbufshft,&
&                         mpi_enreg%comm_img,ierr)
   else
     call xmpi_gatherv(rbuffer,rsize_img,rbuffer_all,rsize_img_all,rbufshft,&
&                      master_img,mpi_enreg%comm_img,ierr)
   end if
   ABI_DEALLOCATE(rbuffer)
   ABI_DEALLOCATE(rsize_img_all)

!  Transfer buffers into gathered array_img_all (master proc only)
   if (use_array_all) then
     ABI_ALLOCATE(iimg,(mpi_enreg%nproc_img))
     iimg=0
     do jj=1,nimagetot
!      The following line supposes that images are sorted by increasing index
       iproc=mpi_enreg%distrb_img(jj)+1;iimg(iproc)=iimg(iproc)+1
       ibufr=rbufshft(iproc)+(iimg(iproc)-1)*rsize
       array_img_all(1:size1,jj)=reshape(rbuffer_all(ibufr+1:ibufr+rsize),(/size1/))
     end do
     ABI_DEALLOCATE(iimg)
   end if

!  Free memory
   ABI_DEALLOCATE(rbufshft)
   ABI_DEALLOCATE(rbuffer_all)

 end if

end subroutine gather_array_img_1D
!!***

!----------------------------------------------------------------------

!!****f* m_results_img/gather_array_img_2D
!! NAME
!!  gather_array_img_2D
!!
!! FUNCTION
!!  Gather an real 2D-array (part of a results_img datastructure) using communicator
!!  over images (replicas) of the cell.
!!  Each contribution of single processor is gathered into a big array on master processor
!!
!! INPUTS
!!  allgather= --optional, default=false--  if TRUE do ALL_GATHER instead of GATHER
!!  master= --optional, default=0-- index of master proc receiving gathered data (if allgather=false)
!!  mpi_enreg=information about MPI parallelization
!!  only_one_per_img= --optional, default=true--  if TRUE, the gather operation
!!                    is only done by one proc per image (master of the comm_cell)
!!  array_img(:,:,:)= (real) 2D-array distributed (has 3 dimensions; the 3rd one is nimage)
!!
!! SIDE EFFECTS
!!  array_img_all(:,:,:)= (real) global (gathered) 2D-array
!!                        (has 3 dimensions; the 3rd one is nimagetot)
!!
!! PARENTS
!!
!! CHILDREN
!!      mkrdim,xred2xcart
!!
!! SOURCE

subroutine gather_array_img_2D(array_img,array_img_all,mpi_enreg,&
&                              master,allgather,only_one_per_img) ! optional arguments

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: master
 logical,optional,intent(in) :: allgather,only_one_per_img
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 real(dp),intent(in) :: array_img(:,:,:)
 real(dp),intent(inout) :: array_img_all(:,:,:)

!Local variables-------------------------------
!scalars
 integer :: ibufr,ierr,iproc,jj
 integer :: master_all,master_img,master_one_img,nimage,nimagetot
 integer :: rsize,rsize_img,size1,size2
 logical :: do_allgather,i_am_master,one_per_img,use_array_all
 !character(len=500) :: msg
!arrays
 integer,allocatable :: iimg(:),nimage_all(:),rbufshft(:),rsize_img_all(:)
 real(dp),allocatable :: rbuffer(:),rbuffer_all(:)

!************************************************************************

 !@results_img_type

 one_per_img=.true.;if (present(only_one_per_img)) one_per_img=only_one_per_img
 do_allgather=.false.;if (present(allgather)) do_allgather=allgather
 master_all=0;if (present(master)) master_all=master
 i_am_master=(mpi_enreg%me==master_all)

 master_img=0;master_one_img=0
 use_array_all= &
&  (((     do_allgather).and.(     one_per_img).and.(mpi_enreg%me_cell==master_one_img)) .or. &
&   ((     do_allgather).and.(.not.one_per_img))                                            .or. &
&   ((.not.do_allgather).and.(     one_per_img).and.(mpi_enreg%me==master_all))             .or. &
&   ((.not.do_allgather).and.(.not.one_per_img).and.(mpi_enreg%me_img==master_img)))

 size1=size(array_img,1);size2=size(array_img,2)
 if (use_array_all) then
   if (size(array_img_all,1)/=size1.or.size(array_img_all,2)/=size2) then
     MSG_BUG('Wrong array_img_all size (1)!')
   end if
 end if

 if ((.not.one_per_img).or.(mpi_enreg%me_cell==master_one_img)) then

!  Simple copy in case of 1 image
   if (use_array_all) then
     if (size(array_img_all,3)<=1) then
       array_img_all(:,:,1)=array_img(:,:,1)
       return
     end if
   endif

!  Gather number of images treated by each proc
   nimage=size(array_img,3)
   ABI_ALLOCATE(nimage_all,(mpi_enreg%nproc_img))
   call xmpi_allgather(nimage,nimage_all,mpi_enreg%comm_img,ierr)
   nimagetot=sum(nimage_all)
   if (use_array_all) then
     if (size(array_img_all,3)/=nimagetot) then
       MSG_BUG('Wrong array_img_all size (2)!')
     endif
   end if

!  Compute number of data
   rsize=size1*size2;rsize_img=nimage*rsize
   ABI_ALLOCATE(rsize_img_all,(mpi_enreg%nproc_img))
   rsize_img_all(:)=rsize*nimage_all(:)
   ABI_DEALLOCATE(nimage_all)

!  Compute shifts in buffer arrays for each proc
   ABI_ALLOCATE(rbufshft,(mpi_enreg%nproc_img))
   rbufshft(1)=0
   do jj=2,mpi_enreg%nproc_img
     rbufshft(jj)=rbufshft(jj-1)+rsize_img_all(jj-1)
   end do

!  Load buffers
   ABI_ALLOCATE(rbuffer,(rsize_img))
   ibufr=0
   do jj=1,nimage
     rbuffer(ibufr+1:ibufr+rsize)=reshape(array_img(1:size1,1:size2,jj),(/rsize/))
     ibufr=ibufr+rsize
   end do

!  Gather all data
   if (use_array_all)  then
     ABI_ALLOCATE(rbuffer_all,(rsize*nimagetot))
   end if
   if (.not.use_array_all)  then
     ABI_ALLOCATE(rbuffer_all,(0))
   end if
   if (do_allgather) then
     call xmpi_allgatherv(rbuffer,rsize_img,rbuffer_all,rsize_img_all,rbufshft,&
&                         mpi_enreg%comm_img,ierr)
   else
     call xmpi_gatherv(rbuffer,rsize_img,rbuffer_all,rsize_img_all,rbufshft,&
&                      master_img,mpi_enreg%comm_img,ierr)
   end if
   ABI_DEALLOCATE(rbuffer)
   ABI_DEALLOCATE(rsize_img_all)

!  Transfer buffers into gathered array_img_all (master proc only)
   if (use_array_all) then
     ABI_ALLOCATE(iimg,(mpi_enreg%nproc_img))
     iimg=0
     do jj=1,nimagetot
!      The following line supposes that images are sorted by increasing index
       iproc=mpi_enreg%distrb_img(jj)+1;iimg(iproc)=iimg(iproc)+1
       ibufr=rbufshft(iproc)+(iimg(iproc)-1)*rsize
       array_img_all(1:size1,1:size2,jj)=reshape(rbuffer_all(ibufr+1:ibufr+rsize),(/size1,size2/))
     end do
     ABI_DEALLOCATE(iimg)
   end if

!  Free memory
   ABI_DEALLOCATE(rbufshft)
   ABI_DEALLOCATE(rbuffer_all)

 end if

end subroutine gather_array_img_2D
!!***

!----------------------------------------------------------------------

!!****f* m_results_img/scatter_array_img
!! NAME
!!  scatter_array_img
!!
!! FUNCTION
!!  Scatter an real 2D-array (part of a results_img datastructure) using communicator
!!  over images (replicas) of the cell.
!!  A big array on master processor is scattered into each contribution
!!  of single processor is gathered
!!
!! INPUTS
!!  master= --optional, default=0-- index of master proc sending data
!!  mpi_enreg=information about MPI parallelization
!!  only_one_per_img= --optional, default=true--  if TRUE, the scatter operation
!!                    is only done by one proc per image (master of the comm_cell)
!!  array_img_all(:,:,:)= (real) global 2D-array (has 3 dimensions; the 3rd one is nimagetot)
!!
!! SIDE EFFECTS
!!  array_img(:,:,:)= (real) distributed 2D-array
!!                    (has 3 dimensions; the 3rd one is nimage)
!!
!! PARENTS
!!      predict_pimd
!!
!! CHILDREN
!!      mkrdim,xred2xcart
!!
!! SOURCE

subroutine scatter_array_img(array_img,array_img_all,mpi_enreg,&
&                            master,only_one_per_img) ! optional arguments

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: master
 logical,optional,intent(in) :: only_one_per_img
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 real(dp),intent(inout) :: array_img(:,:,:)
 real(dp),intent(in) :: array_img_all(:,:,:)

!Local variables-------------------------------
!scalars
 integer :: ibufr,ierr,iproc,jj
 integer :: master_all,master_img,master_one_img,nimage,nimagetot
 integer :: rsize,rsize_img,size1,size2
 logical :: i_am_master,one_per_img,use_array_all
 !character(len=500) :: msg
!arrays
 integer,allocatable :: iimg(:),nimage_all(:),rbufshft(:),rsize_img_all(:)
 real(dp),allocatable :: rbuffer(:),rbuffer_all(:)

!************************************************************************

 !@results_img_type

 one_per_img=.true.;if (present(only_one_per_img)) one_per_img=only_one_per_img
 master_all=0;if (present(master)) master_all=master
 i_am_master=(mpi_enreg%me==master_all)

 use_array_all=i_am_master
 master_img=0;master_one_img=0

 size1=size(array_img,1);size2=size(array_img,2)
 if (use_array_all) then
   if (size(array_img_all,1)/=size1.or.size(array_img_all,2)/=size2) then
     MSG_BUG('Wrong array_img_all size (1)!')
   end if
 end if

 if ((.not.one_per_img).or.(mpi_enreg%me_cell==master_one_img)) then

!  Compute (by gather operation) total number of images
   nimage=size(array_img,3)
   ABI_ALLOCATE(nimage_all,(mpi_enreg%nproc_img))
   call xmpi_allgather(nimage,nimage_all,mpi_enreg%comm_img,ierr)
   nimagetot=sum(nimage_all)
   if (use_array_all) then
     if (size(array_img_all,3)/=nimagetot) then
       MSG_BUG('Wrong array_img_all size (2)!')
     endif
   end if

!  Simple copy in case of 1 image
   if (nimagetot<=1) then
     if (use_array_all) array_img(:,:,1)=array_img_all(:,:,1)

   else

!    Compute number of data
     rsize=size1*size2;rsize_img=nimage*rsize
     ABI_ALLOCATE(rsize_img_all,(mpi_enreg%nproc_img))
     rsize_img_all(:)=rsize*nimage_all(:)

!    Compute shifts in buffer arrays for each proc
     ABI_ALLOCATE(rbufshft,(mpi_enreg%nproc_img))
     rbufshft(1)=0
     do jj=2,mpi_enreg%nproc_img
       rbufshft(jj)=rbufshft(jj-1)+rsize_img_all(jj-1)
     end do

!    Load buffer
     if (use_array_all)  then
       ABI_ALLOCATE(rbuffer_all,(rsize*nimagetot))
     end if
     if (.not.use_array_all)  then
       ABI_ALLOCATE(rbuffer_all,(0))
     end if
     if (use_array_all) then
       ABI_ALLOCATE(iimg,(mpi_enreg%nproc_img))
       iimg=0
       do jj=1,nimagetot
!        The following line supposes that images are sorted by increasing index
         iproc=mpi_enreg%distrb_img(jj)+1;iimg(iproc)=iimg(iproc)+1
         ibufr=rbufshft(iproc)+(iimg(iproc)-1)*rsize
         rbuffer_all(ibufr+1:ibufr+rsize)=reshape(array_img_all(1:size1,1:size2,jj),(/rsize/))
       end do
       ABI_DEALLOCATE(iimg)
      end if

!    Scatter all data
     ABI_ALLOCATE(rbuffer,(rsize_img))
     call xmpi_scatterv(rbuffer_all,rsize_img_all,rbufshft,rbuffer,rsize_img,&
&                       master_img,mpi_enreg%comm_img,ierr)
     ABI_DEALLOCATE(rbuffer_all)
     ABI_DEALLOCATE(rbufshft)
     ABI_DEALLOCATE(rsize_img_all)

!    Transfered distributed buffers into array_img (master proc only)
     ibufr=0
     do jj=1,nimage
       array_img(1:size1,1:size2,jj)=reshape(rbuffer(ibufr+1:ibufr+rsize),(/size1,size2/))
       ibufr=ibufr+rsize
     end do
     ABI_DEALLOCATE(rbuffer)

   end if ! nimagetot<=1
   ABI_DEALLOCATE(nimage_all)

!  Now, is requested dispatch data inside each image
   if (.not.one_per_img) then
     call xmpi_bcast(array_img,0,mpi_enreg%comm_cell,ierr)
   end if

 end if

end subroutine scatter_array_img
!!***

!----------------------------------------------------------------------

!!****f* m_results_img/get_geometry_img
!! NAME
!!  get_geometry_img
!!
!! FUNCTION
!!  From a results_img datastructure, get geometry related data
!!
!! INPUTS
!!  natom=number of atoms
!!  nimage=number of images (including static ones)
!!  results_img(nimage)=datastructure that hold data for each image
!!                      (positions, forces, energy, ...)
!!
!! OUTPUT
!!  etotal(3,natom,nimage)=total energy of each image
!!  fcart(3,natom,nimage)=cartesian forces in each image
!!  rprimd(3,3,nimage)   =dimensional primitive translations for each image
!!  xcart(3,natom,nimage)=cartesian coordinates of atoms in each image
!!  xred(3,natom,nimage) =reduced coordinates of atoms in each image
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      predict_neb,predict_steepest,predict_string
!!
!! CHILDREN
!!      mkrdim,xred2xcart
!!
!! SOURCE

subroutine get_geometry_img(etotal,natom,nimage,results_img,fcart,rprimd,xcart,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nimage
!arrays
 real(dp),intent(out) :: etotal(nimage),fcart(3,natom,nimage),rprimd(3,3,nimage)
 real(dp),intent(out) :: xcart(3,natom,nimage),xred(3,natom,nimage)
 type(results_img_type),intent(in) :: results_img(nimage)
!Local variables-------------------------------
!scalars
 integer :: iimage
!arrays
 real(dp) :: acell(3),rprim(3,3)

!************************************************************************

 do iimage=1,nimage
   acell(:)  =results_img(iimage)%acell(:)
   rprim(:,:)=results_img(iimage)%rprim(:,:)
   xred (:,:,iimage)=results_img(iimage)%xred(:,:)
   fcart(:,:,iimage)=results_img(iimage)%results_gs%fcart(:,:)
   call mkrdim(acell,rprim,rprimd(:,:,iimage))
   call xred2xcart(natom,rprimd(:,:,iimage),xcart(:,:,iimage),xred(:,:,iimage))
   etotal(iimage)=results_img(iimage)%results_gs%etotal
 end do

end subroutine get_geometry_img
!!***

END MODULE m_results_img
!!***
