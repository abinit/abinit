!!****m* ABINIT/m_tetrahedron
!! NAME
!! m_tetrahedron
!!
!! FUNCTION
!!  module for tetrahedron interpolation of DOS and similar quantities
!!  depends on sort_tetra and on m_kpt_rank
!!
!! COPYRIGHT
!!  Copyright (C) 2010-2020 ABINIT group (MJV)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!!  1) Test carefully the case of degenerate tethraedron
!!  2) Change API so that we can pass the energy mesh instead of omega_min and omega_max
!!  3) Add table ik_ibz --> tetra_list to avoid cycling inside big loop over ntetra
!!  4) Add options to get only delta and/or theta ?
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "libtetra.h"

module m_tetrahedron

 USE_MEMORY_PROFILING
 USE_MSG_HANDLING
 use m_krank
#ifdef HAVE_MPI2
 use mpi
#endif
#ifdef HAVE_LIBTETRA_ABINIT
 use m_io_tools, only : open_file
 use m_xmpi
#endif

implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

private
!!***

integer, parameter :: dp = kind(1.0d0)

real(dp),parameter  :: tol6 = 1.d-14, tol14 = 1.d-14, zero = 0.d0, one = 1.d0

real(dp), parameter :: sqrtpi = 1.7724538509055159d0


!!****t* m_tetrahedron/t_tetrahedron
!! NAME
!! t_tetrahedron
!!
!! FUNCTION
!! tetrahedron geometry object
!!
!! SOURCE

type, public :: t_tetrahedron

  integer :: ntetra = 0
  ! Number of tetrahedra

  real(dp)  :: vv
  ! volume of the tetrahedra

  real(dp) :: klatt(3, 3)
  ! reciprocal of lattice vectors for full kpoint grid

  integer,allocatable :: tetra_full(:,:,:)
  !(4,2,ntetra)
  ! For each tetra
  !   (:,1,itetra) indices of the vertex in IBZ (symmetrical image)
  !   (:,2,itetra) indices of the vertexes in the BZ

  integer,allocatable :: tetra_mult(:)
  !(ntetra)
  ! multiplicity of each irred tetrahedron

  integer,allocatable :: tetra_wrap(:,:,:)
  !(3,4,ntetra)
  ! flag to wrap tetrahedron summit into IBZ

  integer,allocatable :: ibz_tetra_count(:)
  ! ibz_tetra_mapping(nkpt_ibz)
  ! Number of tetrahedra associated to a point in the IBZ.

  integer,allocatable :: ibz_tetra_mapping(:,:)
  ! ibz_tetra_mapping(nkpt_ibz, maxval(tetra%ibz_tetra_count)))
  ! map ikbz to tetra index.

end type t_tetrahedron

public :: init_tetra               ! Initialize the object
                                   ! See also the high-level interface tetra_from_kptrlatt provided by m_kpts.
public :: get_tetra_weight         ! Calculate integration weights and their derivatives. shape (nkpt, nene).
public :: tetra_blochl_weights     ! Same as in get_tetra_weight but weights have shape (nene, nkpt).
public :: get_dbl_tetra_weight     ! Calculate integration weights for double tetrahedron integration of delta functions.
                                   ! (NB these correspond to the derivative terms in normal tetrahedron).
public :: destroy_tetra            ! Free memory.
public :: tetra_write              ! Write text file (XML format) with tetra info.
public :: tetralib_has_mpi         ! Return True if the library has been compiled with MPI support.
public :: tetra_get_onewk          ! Calculate integration weights and their derivatives for a single k-point in the IBZ.
public :: tetra_get_onewk_wvals    ! Similar to tetra_get_onewk_wvalsa but reveives arbitrary list of frequency points.
public :: tetra_get_onetetra_wvals ! Get weights for one tetrahedra with arbitrary list of frequency points
!!***

contains
!!***

!----------------------------------------------------------------------

!!****f* m_tetrahedron/destroy_tetra
!! NAME
!! destroy_tetra
!!
!! FUNCTION
!! deallocate tetrahedra pointers if needed
!!
!! PARENTS
!!      m_epweights,m_phgamma,m_thmeig,m_unittests
!!
!! CHILDREN
!!      get_onetetra_,sort_tetra
!!
!! SOURCE

subroutine destroy_tetra (tetra)

 type(t_tetrahedron), intent(inout) :: tetra

 if (allocated(tetra%tetra_full)) then
   TETRA_DEALLOCATE(tetra%tetra_full)
 end if
 if (allocated(tetra%tetra_mult)) then
   TETRA_DEALLOCATE(tetra%tetra_mult)
 end if
 if (allocated(tetra%tetra_wrap)) then
   TETRA_DEALLOCATE(tetra%tetra_wrap)
 end if
 if (allocated(tetra%ibz_tetra_count)) then
   TETRA_DEALLOCATE(tetra%ibz_tetra_count)
 end if
 if (allocated(tetra%ibz_tetra_mapping)) then
   TETRA_DEALLOCATE(tetra%ibz_tetra_mapping)
 end if

end subroutine destroy_tetra
!!***

!----------------------------------------------------------------------

!!****f* m_tetrahedron/init_tetra
!! NAME
!! init_tetra
!!
!! FUNCTION
!! get tetrahedra characterized by apexes
!!
!! INPUTS
!!  indkpt(nkpt_fullbz)=indexes of irred kpoints equivalent to kpt_fullbz
!!  gprimd(3,3) = reciprocal space vectors
!!  klatt(3,3)=reciprocal of lattice vectors for full kpoint grid
!!  kpt_fullbz(3,nkpt_fullbz)=kpoints in full brillouin zone
!!  nkpt_fullbz=number of kpoints in full brillouin zone
!!  comm= MPI communicator
!!
!! OUTPUT
!!  tetra%tetra_full(4,2,ntetra)=for each tetrahedron,
!!     the different instances of the tetrahedron (fullbz kpoints)
!!  tetra%tetra_mult(ntetra) = store multiplicity of each irred tetrahedron
!!  tetra%tetra_wrap(3,4,ntetra) = store flag to wrap tetrahedron summit into IBZ
!!  tetra%ntetra = final number of irred tetra (dimensions of tetra_* remain larger)
!!  tetra%vv = tetrahedron volume divided by full BZ volume
!!
!! PARENTS
!!      m_epweights,m_phgamma,m_thmeig,m_unittests
!!
!! CHILDREN
!!      get_onetetra_,sort_tetra
!!
!! SOURCE

subroutine init_tetra(indkpt, gprimd, klatt, kpt_fullbz, nkpt_fullbz, tetra, ierr, errorstring, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt_fullbz, comm
 integer, intent(out) :: ierr
 character(len=80), intent(out) :: errorstring
 type(t_tetrahedron),intent(out) :: tetra
!arrays
 integer,intent(in) :: indkpt(nkpt_fullbz)
 real(dp) ,intent(in) :: gprimd(3,3),klatt(3,3),kpt_fullbz(3,nkpt_fullbz)

!Local variables-------------------------------
!scalars
 integer :: ialltetra,ikpt2,ikpt_full,isummit,itetra,jalltetra,jsummit
 integer :: ii,jj,ikibz,nkpt_ibz, my_rank, nprocs
 integer :: symrankkpt,mtetra,itmp,ntetra_irred
 real(dp) :: shift1,shift2,shift3, rcvol,hashfactor
 !real :: cpu_start, cpu_stop
 type(krank_t) :: krank
!arrays
 integer :: ind_ibz(4), tetra_shifts(3,4,6)  ! 3 dimensions, 4 summits, and 6 tetrahedra / kpoint box
 real(dp)  :: k1(3),k2(3),k3(3)
 integer,allocatable :: tetra_full_(:,:,:)
 integer,allocatable :: tetra_mult_(:)
 integer,allocatable :: tetra_wrap_(:,:,:)
 integer, allocatable :: reforder(:)
 integer, allocatable :: irred_itetra(:)
 real(dp), allocatable :: tetra_hash(:)

! *********************************************************************

 !call cpu_time(cpu_start)

 my_rank = 0; nprocs = 1
#ifdef HAVE_LIBTETRA_ABINIT
 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
#endif

 ierr = 0
 errorstring = ""
!jmb
 shift1 = zero
 shift2 = zero
 shift3 = zero

 tetra%klatt = klatt

 mtetra = 6 * nkpt_fullbz
 TETRA_ALLOCATE(tetra_full_, (4,2,mtetra))
 TETRA_ALLOCATE(tetra_mult_, (mtetra))
 TETRA_ALLOCATE(tetra_wrap_, (3,4,mtetra))

 tetra_mult_ = 1
 tetra_full_ = 0
 tetra_wrap_ = 0

! tetra_shifts(:,1,1) = (/0,0,0/)
! tetra_shifts(:,2,1) = (/0,1,0/)
! tetra_shifts(:,3,1) = (/0,1,1/)
! tetra_shifts(:,4,1) = (/1,1,0/)
! tetra_shifts(:,1,2) = (/0,0,0/)
! tetra_shifts(:,2,2) = (/0,1,1/)
! tetra_shifts(:,3,2) = (/1,1,0/)
! tetra_shifts(:,4,2) = (/1,1,1/)
! tetra_shifts(:,1,3) = (/0,0,0/)
! tetra_shifts(:,2,3) = (/1,0,0/)
! tetra_shifts(:,3,3) = (/1,1,0/)
! tetra_shifts(:,4,3) = (/1,1,1/)
! tetra_shifts(:,1,4) = (/0,0,0/)
! tetra_shifts(:,2,4) = (/0,0,1/)
! tetra_shifts(:,3,4) = (/1,0,0/)
! tetra_shifts(:,4,4) = (/1,1,1/)
! tetra_shifts(:,1,5) = (/0,0,1/)
! tetra_shifts(:,2,5) = (/1,0,0/)
! tetra_shifts(:,3,5) = (/1,0,1/)
! tetra_shifts(:,4,5) = (/1,1,1/)
! tetra_shifts(:,1,6) = (/0,0,0/)
! tetra_shifts(:,2,6) = (/0,0,1/)
! tetra_shifts(:,3,6) = (/0,1,1/)
! tetra_shifts(:,4,6) = (/1,1,1/)

 ! bxu, the following division scheme is according to Bloechl's paper
 tetra_shifts(:,1,1) = (/0,0,0/)
 tetra_shifts(:,2,1) = (/1,0,0/)
 tetra_shifts(:,3,1) = (/0,1,0/)
 tetra_shifts(:,4,1) = (/1,0,1/)
 tetra_shifts(:,1,2) = (/1,0,0/)
 tetra_shifts(:,2,2) = (/1,1,0/)
 tetra_shifts(:,3,2) = (/0,1,0/)
 tetra_shifts(:,4,2) = (/1,0,1/)
 tetra_shifts(:,1,3) = (/0,1,0/)
 tetra_shifts(:,2,3) = (/1,1,0/)
 tetra_shifts(:,3,3) = (/1,0,1/)
 tetra_shifts(:,4,3) = (/1,1,1/)
 tetra_shifts(:,1,4) = (/0,0,0/)
 tetra_shifts(:,2,4) = (/0,1,0/)
 tetra_shifts(:,3,4) = (/0,0,1/)
 tetra_shifts(:,4,4) = (/1,0,1/)
 tetra_shifts(:,1,5) = (/0,0,1/)
 tetra_shifts(:,2,5) = (/1,0,1/)
 tetra_shifts(:,3,5) = (/0,1,0/)
 tetra_shifts(:,4,5) = (/0,1,1/)
 tetra_shifts(:,1,6) = (/0,1,0/)
 tetra_shifts(:,2,6) = (/1,0,1/)
 tetra_shifts(:,3,6) = (/0,1,1/)
 tetra_shifts(:,4,6) = (/1,1,1/)

 ! Make full k-point rank arrays
 ! TODO: Lot of memory allocated here if dense mesh e.g ~ 300 ** 3
 krank = krank_new(nkpt_fullbz, kpt_fullbz)

 ialltetra = 1
 do ikpt_full=1,nkpt_fullbz
   do itetra=1,6
     !ialltetra = itetra + (ikpt_full -1) * 6
     !if (mod(ialltetra, nprocs) /= my_rank) cycle ! MPI parallelism.
     do isummit=1,4
       k1(:) = kpt_fullbz(:,ikpt_full) &
        + tetra_shifts(1,isummit,itetra)*klatt(:,1) &
        + tetra_shifts(2,isummit,itetra)*klatt(:,2) &
        + tetra_shifts(3,isummit,itetra)*klatt(:,3)

       ! Find full kpoint which is summit isummit of tetrahedron itetra around full kpt ikpt_full !
       symrankkpt =  krank%get_rank(k1)
       ikpt2 = krank%invrank(symrankkpt)
       if (ikpt2 < 1) then
         errorstring='Error in ranking k-points - exiting with un-initialized tetrahedra.'
         ierr = 2
         call krank%free()
         TETRA_ALLOCATE(tetra%tetra_full, (4,2,1))
         TETRA_ALLOCATE(tetra%tetra_mult, (1))
         TETRA_ALLOCATE(tetra%tetra_wrap, (3,4,1))
         TETRA_DEALLOCATE(tetra_full_)
         TETRA_DEALLOCATE(tetra_mult_)
         TETRA_DEALLOCATE(tetra_wrap_)
         return
       end if

       ! Store irreducible kpoint equivalent to kpt_fullbz(:,ikpt2)
       tetra_full_(isummit,1,ialltetra) = indkpt(ikpt2)
       tetra_full_(isummit,2,ialltetra) = ikpt2
       shift1 = k1(1)-kpt_fullbz(1,ikpt2)
       shift2 = k1(2)-kpt_fullbz(2,ikpt2)
       shift3 = k1(3)-kpt_fullbz(3,ikpt2)
       if (shift1>0.5d0) then
         tetra_wrap_(1,isummit,ialltetra) = 1
       else if (shift1<-0.5d0) then
         tetra_wrap_(1,isummit,ialltetra) = -1
       end if
       if (shift2>0.5d0) then
         tetra_wrap_(2,isummit,ialltetra) = 1
       else if (shift2<-0.5d0) then
         tetra_wrap_(2,isummit,ialltetra) = -1
       end if
       if (shift3>0.5d0) then
         tetra_wrap_(3,isummit,ialltetra) = 1
       else if (shift3<-0.5d0) then
         tetra_wrap_(3,isummit,ialltetra) = -1
       end if

       ! sort itetra summits
       ! TODO: replace with sort_int
       do jsummit=isummit,2,-1
         if ( tetra_full_(jsummit,1,ialltetra)  <  tetra_full_(jsummit-1,1,ialltetra) ) then
           itmp = tetra_full_(jsummit,1,ialltetra)
           tetra_full_(jsummit,1,ialltetra) = tetra_full_(jsummit-1,1,ialltetra)
           tetra_full_(jsummit-1,1,ialltetra) = itmp
           itmp = tetra_full_(jsummit,2,ialltetra)
           tetra_full_(jsummit,2,ialltetra) = tetra_full_(jsummit-1,2,ialltetra)
           tetra_full_(jsummit-1,2,ialltetra) = itmp
           ! keep fullbz_kpt tetrahedra points in same order
           itmp = tetra_wrap_(1,jsummit,ialltetra)
           tetra_wrap_(1,jsummit,ialltetra) = tetra_wrap_(1,jsummit-1,ialltetra)
           tetra_wrap_(1,jsummit-1,ialltetra) = itmp
           itmp = tetra_wrap_(2,jsummit,ialltetra)
           tetra_wrap_(2,jsummit,ialltetra) = tetra_wrap_(2,jsummit-1,ialltetra)
           tetra_wrap_(2,jsummit-1,ialltetra) = itmp
           itmp = tetra_wrap_(1,jsummit,ialltetra)
           tetra_wrap_(3,jsummit,ialltetra) = tetra_wrap_(3,jsummit-1,ialltetra)
           tetra_wrap_(3,jsummit-1,ialltetra) = itmp
         end if
       end do ! jsummit

     end do ! isummit

     if (ialltetra > mtetra) then
       write (errorstring, '(3a,i0,a,i0)' ) &
        'init_tetra: BUG - ',&
        ' ialltetra > mtetra ',&
        ' ialltetra=  ',ialltetra,', mtetra= ',mtetra
       ierr = 1
       return
     end if
     ialltetra = ialltetra+1
   end do ! itetra
 end do ! ikpt_full

 !call cpu_time(cpu_stop)
 !write(*,*)"tetra_init ikpt_loop:", cpu_stop - cpu_start
 !cpu_start = cpu_stop

 call krank%free()

 rcvol = abs (gprimd(1,1)*(gprimd(2,2)*gprimd(3,3)-gprimd(3,2)*gprimd(2,3)) &
& -gprimd(2,1)*(gprimd(1,2)*gprimd(3,3)-gprimd(3,2)*gprimd(1,3)) &
& +gprimd(3,1)*(gprimd(1,2)*gprimd(2,3)-gprimd(2,2)*gprimd(1,3)))

 ! Volume of all tetrahedra should be the same as that of tetra 1
 ! this is the volume of 1 tetrahedron, should be coherent with notation in Lehmann & Taut
 k1(:) = gprimd(:,1)*klatt(1,1) +  gprimd(:,2)*klatt(2,1) +  gprimd(:,3)*klatt(3,1)
 k2(:) = gprimd(:,1)*klatt(1,2) +  gprimd(:,2)*klatt(2,2) +  gprimd(:,3)*klatt(3,2)
 k3(:) = gprimd(:,1)*klatt(1,3) +  gprimd(:,2)*klatt(2,3) +  gprimd(:,3)*klatt(3,3)
 tetra%vv  = abs (k1(1)*(k2(2)*k3(3)-k2(3)*k3(2)) &
& -k1(2)*(k2(1)*k3(3)-k2(3)*k3(1)) &
& +k1(3)*(k2(1)*k3(2)-k2(2)*k3(1))) / 6.d0 / rcvol

 ! eliminate equivalent tetrahedra by symmetry and account for them in multiplicity tetra_mult
 tetra%ntetra = mtetra

 ! FIXME: could we replace this with a ranking algorithm to avoid the O(tetra%ntetra^2) step? For example:
 ! get tetrahedron rank - problem too many combinations in principle = nkpt_irred^4 - only a few used in practice
 ! sort ranks and keep indices

 ! make hash table = tetra_full_(1)*nkptirred**3+tetra_full_(2)*nkptirred**2+tetra_full_(3)*nkptirred**1+tetra_full_(4)

 hashfactor = 100.d0 ! *acos(-1.d0) ! 100 pi should be far from an integer...
 TETRA_ALLOCATE(tetra_hash, (tetra%ntetra))
 TETRA_ALLOCATE(reforder, (tetra%ntetra))

 !MG: In principle the order of the indices should not matter.
 do ialltetra=1, tetra%ntetra
   tetra_hash(ialltetra) = tetra_full_(1,1,ialltetra)*hashfactor**3+&
&      tetra_full_(2,1,ialltetra)*hashfactor**2+&
&      tetra_full_(3,1,ialltetra)*hashfactor**1+&
&      tetra_full_(4,1,ialltetra)
   reforder(ialltetra) = ialltetra
 end do

 call sort_tetra(tetra%ntetra, tetra_hash, reforder, tol6)
 ! Most of the wall-time is spent in the  preamble of this routine (up to this point).
 ! sort_tetra is not easy to parallelize...

 ! determine number of tetra after reduction
 TETRA_ALLOCATE(irred_itetra, (tetra%ntetra))
 jalltetra = 1
 irred_itetra(1) = 1
 do ialltetra=2, tetra%ntetra
   if (abs(tetra_hash(ialltetra)-tetra_hash(ialltetra-1)) > tol6) then
     ! found a new series
     jalltetra = jalltetra + 1
   end if
   irred_itetra(ialltetra) = jalltetra
 end do

 ! reset number of tetra
 ntetra_irred = jalltetra

 ! allocate definitive tetra arrays and transfer to new arrays
 TETRA_ALLOCATE(tetra%tetra_full, (4,2,ntetra_irred))
 TETRA_ALLOCATE(tetra%tetra_mult, (ntetra_irred))
 TETRA_ALLOCATE(tetra%tetra_wrap, (3,4,ntetra_irred))

 ! eliminate equal rank tetrahedra and accumulate multiplicity into first one
 tetra%tetra_full = 0
 tetra%tetra_mult = 0
 tetra%tetra_wrap = 0
 jalltetra = 1
 tetra%tetra_full(:,:,1) = tetra_full_(:,:,reforder(1))
 tetra%tetra_mult(1) = 1
 tetra%tetra_wrap(:,:,1) = tetra_wrap_(:,:,reforder(1))
 do ialltetra=2, tetra%ntetra
   ! TODO: check if tolerance is adapted
   if (abs(tetra_hash(ialltetra)-tetra_hash(ialltetra-1)) > tol6) then
     ! found a new series
     jalltetra = jalltetra + 1
     tetra%tetra_full(:,:,jalltetra) = tetra_full_(:,:,reforder(ialltetra))
     tetra%tetra_wrap(:,:,jalltetra) = tetra_wrap_(:,:,reforder(ialltetra))
     tetra%tetra_mult(jalltetra) = 1
   else
     ! TODO: add real check that the tetra are equivalent...
     ! otherwise increment jalltetra here as well, generate new series?
     tetra%tetra_mult(jalltetra) = tetra%tetra_mult(jalltetra) + tetra_mult_(reforder(ialltetra))
     !tetra_mult_(reforder(ialltetra)) = 0
   end if
 end do

 ! reset of ntetra for final version after checks and debu
 tetra%ntetra = ntetra_irred

 TETRA_DEALLOCATE(tetra_hash)
 TETRA_DEALLOCATE(reforder)
 TETRA_DEALLOCATE(irred_itetra)
 TETRA_DEALLOCATE(tetra_full_)
 TETRA_DEALLOCATE(tetra_mult_)
 TETRA_DEALLOCATE(tetra_wrap_)

 ! Create mapping between the irreducible k-points
 ! and all the tetrahedron contributing with some weight
 nkpt_ibz = maxval(indkpt)

 ! 1. First we count what is the maximum number of distinct tetrahedra that each k-point contains
 TETRA_ALLOCATE(tetra%ibz_tetra_count,(nkpt_ibz))
 tetra%ibz_tetra_count(:) = 0

 ! Count max tetra contributing
 do ii=1,tetra%ntetra
   ! Here we need the original ordering to reference the correct irred kpoints
   ind_ibz(:) = tetra%tetra_full(:,1,ii)
   ! count max tetra contributing
   do jj=1,4
     ikibz = ind_ibz(jj)
     if (ikibz > nkpt_ibz) cycle
     tetra%ibz_tetra_count(ikibz) = tetra%ibz_tetra_count(ikibz) + 1
   end do
 end do

 ! 2. Then we build mapping of ikbz to tetra
 TETRA_ALLOCATE(tetra%ibz_tetra_mapping,(nkpt_ibz,maxval(tetra%ibz_tetra_count)))
 tetra%ibz_tetra_count(:) = 0
 do ii=1,tetra%ntetra
   ! Here we need the original ordering to reference the correct irred kpoints
   ind_ibz(:) = tetra%tetra_full(:,1,ii)
   ! Use the counter to move pointer and then fill index
   do jj=1,4
     ikibz = ind_ibz(jj)
     if (ikibz > nkpt_ibz) cycle
     ! avoid putting the same index twice
     if (tetra%ibz_tetra_count(ikibz) > 0) then
       if (tetra%ibz_tetra_mapping(ikibz,tetra%ibz_tetra_count(ikibz)) == ii) cycle
     end if
     tetra%ibz_tetra_count(ikibz) = tetra%ibz_tetra_count(ikibz) + 1
     tetra%ibz_tetra_mapping(ikibz,tetra%ibz_tetra_count(ikibz)) = ii
   end do
 end do

 !call cpu_time(cpu_stop)
 !write(*,*)"tetra_init 2nd part:", cpu_stop - cpu_start
 !cpu_start = cpu_stop

end subroutine init_tetra
!!***

!----------------------------------------------------------------------

!!****f* m_tetrahedron/tetra_write
!! NAME
!! tetra_write
!!
!! FUNCTION
!!  Write text file with tetra info.
!!
!! INPUTS
!!  tetra<t_tetrahedron>=tetrahedron geometry object
!!  nkibz=Number of k-points in the IBZ used to generate tetra
!!  kibz(3,nkibz)=Reduced coordinates of the IBZ
!!  path=Name of output file
!!
!! OUTPUT
!!  Output is written to file.
!!
!! PARENTS
!!
!! CHILDREN
!!      get_onetetra_,sort_tetra
!!
!! SOURCE

subroutine tetra_write(tetra, nkibz, kibz, path)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkibz
 character(len=*),intent(in) :: path
 type(t_tetrahedron),intent(in) :: tetra
!arrays
 real(dp),intent(in) :: kibz(3,nkibz)

!Local variables-------------------------------
 integer,parameter :: version=1
 integer :: ik,it,unt
#ifdef HAVE_LIBTETRA_ABINIT
 character(len=500) :: msg
#endif

! *********************************************************************

#ifdef HAVE_LIBTETRA_ABINIT
 if (open_file(file=trim(path),iomsg=msg,newunit=unt,form="formatted",status="unknown",action="write")/=0) then
   TETRA_ERROR(msg)
 end if
#else
 open(file=trim(path),newunit=unt,form="formatted",status="unknown",action="write")
#endif

 write(unt,*)version, " # version number"

 ! Write IBZ
 write(unt,*)nkibz, " # number of k-points in the IBZ"
 write(unt,"(a)")"<irreducible_zone>"
 do ik=1,nkibz
   write(unt,"(3es22.12)") kibz(:,ik)
 end do
 write(unt,"(a)")"</irreducible_zone>"

 ! Write tetra info
 write(unt,"(i0,a)")tetra%ntetra, " # number of tetrahedra"
 write(unt,"(es22.12,a)")tetra%vv, " # tetrahedron volume"

 write(unt,"(a)")"<tetra_full>"
 do it=1,tetra%ntetra
   write(unt,"(8(i0,1x))")tetra%tetra_full(:,:,it)
 end do
 write(unt,"(a)")"</tetra_full>"

 write(unt,"(a)")"<tetra_mult>"
 do it=1,tetra%ntetra
   write(unt,"(i0)")tetra%tetra_mult(it)
 end do
 write(unt,"(a)")"</tetra_mult>"

 write(unt,"(a)")"<tetra_wrap>"
 do it=1,tetra%ntetra
   write(unt,"(12(i0,1x))")tetra%tetra_wrap(:,:,it)
 end do
 write(unt,"(a)")"</tetra_wrap>"

 close(unt)

end subroutine tetra_write
!!***

!----------------------------------------------------------------------

!!****f* m_tetrahedron/get_tetra_weight
!! NAME
!! get_tetra_weight
!!
!! FUNCTION
!! calculate integration weights and their derivatives from Blochl et al PRB 49 16223 [[cite:Bloechl1994a]]
!!
!! INPUTS
!! eigen_in(nkpt)=eigenenergies for each k point
!! enemin=minimal energy for DOS
!! enemax=maximal energy for DOS
!! max_occ=maximal occupation number (2 for nsppol=1, 1 for nsppol=2)
!! nene=number of energies for DOS
!! nkpt=number of irreducible kpoints
!! tetra<t_tetrahedron>
!!   %ntetra=number of tetrahedra
!!   %tetra_full(4,2,ntetra)=for each irred tetrahedron, the list of k point vertices
!!     1 -> irred kpoint   2 -> fullkpt
!!   %tetra_mult(ntetra)=for each irred tetrahedron, its multiplicity
!!   %vv = ratio of volume of one tetrahedron in reciprocal space to full BZ volume
!! bcorr=1 to include Blochl correction else 0.
!! comm=MPI communicator
!!
!! OUTPUT
!!  tweight(nkpt,nene) = integration weights for each irred kpoint from all adjacent tetrahedra
!!  dtweightde(nkpt,nene) = derivative of tweight wrt energy
!!
!! PARENTS
!!      m_epweights,m_thmeig
!!
!! CHILDREN
!!      get_onetetra_,sort_tetra
!!
!! SOURCE

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! THIS FUNCTION IS DEPRECATED, USE tetra_blochl_weights
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine get_tetra_weight(eigen_in,enemin,enemax,max_occ,nene,nkpt,tetra,&
  bcorr,tweight,dtweightde,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nene,nkpt,bcorr,comm
 type(t_tetrahedron), intent(in) :: tetra
 real(dp) ,intent(in) :: enemax,enemin,max_occ
!arrays
 real(dp) ,intent(in) :: eigen_in(nkpt)
 real(dp) ,intent(out) :: dtweightde(nkpt,nene),tweight(nkpt,nene)

!Local variables-------------------------------
 real(dp), allocatable :: dtweightde_ek(:, :), tweight_ek(:, :)

! *********************************************************************

 TETRA_ALLOCATE(dtweightde_ek, (nene, nkpt))
 TETRA_ALLOCATE(tweight_ek, (nene, nkpt))

 call tetra_blochl_weights(tetra,eigen_in,enemin,enemax,max_occ,nene,nkpt,bcorr,tweight_ek,dtweightde_ek,comm)

 ! transpose: otherwise the data access is crap and the code slows by an order of magnitude
 tweight    = transpose(tweight_ek)
 dtweightde = transpose(dtweightde_ek)

 TETRA_DEALLOCATE(dtweightde_ek)
 TETRA_DEALLOCATE(tweight_ek)

end subroutine get_tetra_weight
!!***

!----------------------------------------------------------------------

!!****f* m_tetrahedron/tetra_blochl_weights
!! NAME
!! tetra_blochl_weights
!!
!! FUNCTION
!! calculate integration weights and their derivatives from Blochl et al PRB 49 16223 [[cite:Bloechl1994a]]
!! Same API as get_tetra_weight but output weights here have shape (nene, nkpt)
!!
!! PARENTS
!!      m_tetrahedron,m_unittests
!!
!! CHILDREN
!!      get_onetetra_,sort_tetra
!!
!! SOURCE

subroutine tetra_blochl_weights(tetra,eigen_in,enemin,enemax,max_occ,nene,nkpt,&
  bcorr,tweight_t,dtweightde_t,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nene,nkpt,bcorr,comm
 type(t_tetrahedron), intent(in) :: tetra
 real(dp) ,intent(in) :: enemax,enemin,max_occ
!arrays
 real(dp) ,intent(in) :: eigen_in(nkpt)
 real(dp) ,intent(out) :: dtweightde_t(nene,nkpt),tweight_t(nene,nkpt)

!Local variables-------------------------------
!scalars
 integer :: itetra,nprocs,my_start,my_stop,ierr,ii
!arrays
 integer :: ind_ibz(4)
 real(dp) :: eigen_1tetra(4)
 real(dp), allocatable :: tweight_tmp(:,:),dtweightde_tmp(:,:),buffer(:,:)

! *********************************************************************

 TETRA_ALLOCATE(tweight_tmp, (nene, 4))
 TETRA_ALLOCATE(dtweightde_tmp, (nene, 4))
 tweight_t = zero; dtweightde_t = zero

 call split_work(tetra%ntetra, comm, nprocs, my_start, my_stop, ierr)
 if (ierr /= 0) TETRA_ERROR("Error in MPI layer")

 ! for each tetrahedron
 do itetra=my_start,my_stop
   tweight_tmp = zero
   dtweightde_tmp = zero

   ! Here we need the original ordering to reference the correct irred kpoints
   ind_ibz(:) = tetra%tetra_full(:,1,itetra)

   ! Sort energies before calling get_onetetra_
   eigen_1tetra(:) = eigen_in(ind_ibz(:))
   call sort_tetra(4, eigen_1tetra, ind_ibz, tol14)

   call get_onetetra_(tetra,itetra,eigen_1tetra,enemin,enemax,max_occ,nene,bcorr,tweight_tmp,dtweightde_tmp)

   ! NOTE: the following blas calls are not working systematically, or do not give speed ups, strange...
   !if (nene > 100) then
   !  do ii=1,4
   !    call daxpy (nene, 1.d0, tweight_tmp(:,ii), 1, tweight_t(:,ind_ibz(ii)), 1)
   !  end do
   !  do ii=1,4
   !    call daxpy (nene, 1.d0, dtweightde_tmp(:,ii), 1, dtweightde_t(:,ind_ibz(ii)), 1)
   !  end do
   !else
   do ii=1,4
     tweight_t(:,ind_ibz(ii)) = tweight_t(:,ind_ibz(ii)) + tweight_tmp(:,ii)
   end do
   do ii=1,4
     dtweightde_t(:,ind_ibz(ii)) = dtweightde_t(:,ind_ibz(ii)) + dtweightde_tmp(:,ii)
   end do
   !end if
 end do ! itetra

 TETRA_DEALLOCATE(tweight_tmp)
 TETRA_DEALLOCATE(dtweightde_tmp)

 if (nprocs > 1) then
#ifdef HAVE_MPI
   TETRA_ALLOCATE(buffer, (nene, nkpt))
   call MPI_ALLREDUCE(tweight_t,buffer,nene*nkpt,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
   tweight_t = buffer

   call MPI_ALLREDUCE(dtweightde_t,buffer,nene*nkpt,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
   dtweightde_t = buffer
   TETRA_DEALLOCATE(buffer)
#endif
 end if

end subroutine tetra_blochl_weights
!!***

!----------------------------------------------------------------------

!!****f* m_tetrahedron/get_dbl_tetra_weight
!! NAME
!! get_dbl_tetra_weight
!!
!! FUNCTION
!! calculate integration weights and their derivatives
!! for double tetrahedron method from Allen Phys Stat Sol B 120 529 (1983) [[cite:Allen1983b]]
!! the k-points and tetrahedra must be the same for both grids, of course,
!! but the range of energies is arbitrary
!!
!! Omega is called eigen1 here
!! E is called eigen2 here
!! indexing goes from 1 to 4 for the tetrahedron corners, in order of increasing eigen1
!!  in Allen, from 0 to 3...
!!
!! INPUTS
!! eigen1_in(nkpt)=eigenenergies for each k point
!! eigen2_in(nkpt)=eigenenergies for each k point
!! enemin1=minimal energy for DOS in energy 1
!! enemax1=maximal energy for DOS
!! enemin2=minimal energy for DOS in energy 2
!! enemax2=maximal energy for DOS
!! max_occ=maximal occupation number (2 for nsppol=1, 1 for nsppol=2)
!! nene1=number of energies for DOS in energy 1
!! nene2=number of energies for DOS in energy 2
!! nkpt=number of irreducible kpoints
!! tetra%ntetra=number of tetra
!! tetra%tetra_full(4,2,ntetra)=for each irred tetrahedron, the list of k point vertices
!!   1 -> irred kpoint   2 -> fullkpt
!! tetra%tetra_mult(ntetra)=for each irred tetrahedron, its multiplicity
!! tetra%vv = ratio of volume of one tetrahedron in reciprocal space to full BZ volume
!! ierr = error code on exit
!!
!! OUTPUT
!!  tweight(nkpt,nene1,nene2) = integration weights for each irred kpoint from all adjacent tetrahedra
!!  dtweightde(nkpt,nene1,nene2) = derivative of tweight wrt energy
!!
!! PARENTS
!!      m_phgamma
!!
!! CHILDREN
!!      get_onetetra_,sort_tetra
!!
!! SOURCE

subroutine get_dbl_tetra_weight(eigen1_in,eigen2_in,enemin1,enemax1,enemin2,enemax2,&
&    max_occ,nene1,nene2,nkpt,tetra,tweight,dtweightde, ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nene1,nene2,nkpt
 integer,intent(out) :: ierr
 type(t_tetrahedron), intent(in) :: tetra
 real(dp),intent(in) :: enemax1,enemin1
 real(dp),intent(in) :: enemax2,enemin2
 real(dp),intent(in) :: max_occ
!arrays
 real(dp),intent(in) :: eigen1_in(nkpt)
 real(dp),intent(in) :: eigen2_in(nkpt)
 real(dp),intent(out) :: dtweightde(nkpt,nene1,nene2),tweight(nkpt,nene1,nene2)

!Local variables-------------------------------
!  needed for gaussian replacement of Dirac functions
!  the three coefficients of the DOS as quadratic form,
!    in the interval [eig(ikpt-1), eig(ikpt)]
!    for ikpt = 1 we add a point below eigen(1) which doesnt
!    contribute to the DOS in any tetrahedron
!scalars
 integer :: ieps1,ieps2,itetra
 integer :: nn1_1,nn1_2,nn1_3,nn1_4
 integer :: nn2_1,nn2_2,nn2_3
 integer :: ind_a(3), ind_b(3), ind_c(3)
 real(dp)  :: deltaene1,eps1
 real(dp)  :: deltaene2,eps2
! real(dp)  :: gau_prefactor,gau_width,gau_width2
 real(dp)  :: epsilon1(4,4)
 real(dp)  :: epsilon2(4,4)
 real(dp)  :: inv_epsilon1(4,4)
 real(dp)  :: aa(3),bb(3),cc(3)
 real(dp)  :: delaa(3),delbb(3),delcc(3)
 real(dp)  :: delaa0,delbb0,delcc0
 real(dp)  :: inv_delaa(3),inv_delbb(3),inv_delcc(3)
 real(dp)  :: deleps1, deleps2
 real(dp)  :: inv_deleps1
 real(dp)  :: dccde1, dccde1_pre
 real(dp)  :: volconst,volconst_mult
 real(dp)  :: ii0, ii1, ii3
!arrays
 integer :: ind_k(4)
 real(dp), allocatable :: tweight_tmp(:,:,:)
 real(dp), allocatable :: dtweightde_tmp(:,:,:)
 real(dp)  :: eigen1_1tetra(4)
 real(dp)  :: eigen2_1tetra(4)

! *********************************************************************

 ierr = 0
 if (nene1 <= 1 .or. nene2 <= 1)  then
   !'get_dbl_tetra_weight: nene must be at least 2'
   ierr = 1
   return
 else
   deltaene1 = (enemax1-enemin1) / (nene1-1)
   deltaene2 = (enemax2-enemin2) / (nene2-1)
 end if

 TETRA_ALLOCATE(tweight_tmp, (4, nene2, nene1))
 TETRA_ALLOCATE(dtweightde_tmp, (4, nene2, nene1))

!print *, "warning: for the moment, heaviside weights are 0. The delta function / DOS weights are the only ones calculated "

 volconst = tetra%vv/4.d0

 ! for each tetrahedron
 do itetra=1,tetra%ntetra
   ! these are for 1 tetrahedron only.
   tweight_tmp = zero
   dtweightde_tmp = zero

   volconst_mult = max_occ*volconst*dble(tetra%tetra_mult(itetra))

   ! Here we need the original ordering to reference the correct irred kpoints
   ! ind_k refers to the index in the full k list of the summits of the present tetrahedra
   ! we can forget the order of the summits within the tetrahedron, because eigen1 fixes that
   ! order with its increasing value
   ind_k(1) = tetra%tetra_full(1,1,itetra)
   ind_k(2) = tetra%tetra_full(2,1,itetra)
   ind_k(3) = tetra%tetra_full(3,1,itetra)
   ind_k(4) = tetra%tetra_full(4,1,itetra)
   eigen1_1tetra(1) = eigen1_in(ind_k(1))
   eigen1_1tetra(2) = eigen1_in(ind_k(2))
   eigen1_1tetra(3) = eigen1_in(ind_k(3))
   eigen1_1tetra(4) = eigen1_in(ind_k(4))
   call sort_tetra(4,eigen1_1tetra,ind_k,tol14)

   ! re-sort eigen2 values according to order chosen for eigen1. Eigen2 are _not_ in increasing order!
   eigen2_1tetra(1) = eigen2_in(ind_k(1))
   eigen2_1tetra(2) = eigen2_in(ind_k(2))
   eigen2_1tetra(3) = eigen2_in(ind_k(3))
   eigen2_1tetra(4) = eigen2_in(ind_k(4))

   ! the epsilons are energy differences for the two eigenvalue sets
   epsilon1 = zero
   epsilon2 = zero
   do ieps1 = 1, 4
     do ieps2 = ieps1+1, 4
       epsilon1(ieps1,ieps2) = eigen1_1tetra(ieps1)-eigen1_1tetra(ieps2)
       epsilon1(ieps2,ieps1) = -epsilon1(ieps1,ieps2)
       epsilon2(ieps1,ieps2) = eigen2_1tetra(ieps1)-eigen2_1tetra(ieps2)
       epsilon2(ieps2,ieps1) = -epsilon2(ieps1,ieps2)
     end do
   end do

   ! we precalculate the inverses to avoid doing tons of divisions in the energy loops below
   ! Allen formulae only require the inverses of the differences of eigen1 + the a b c below
   inv_epsilon1 = zero
   do ieps1 = 1, 4
     do ieps2 = ieps1+1, 4
       if (abs(epsilon1(ieps1,ieps2)) > tol6) then
         inv_epsilon1(ieps1,ieps2) = 1.d0 / epsilon1(ieps1,ieps2)
         inv_epsilon1(ieps2,ieps1) = -inv_epsilon1(ieps1,ieps2)
       end if
     end do
   end do

   ! these bounds determine the intervals for Omega in Allen paper, and cases A, B, C
   nn1_1 = int((eigen1_1tetra(1)-enemin1)/deltaene1)+1
   nn1_2 = int((eigen1_1tetra(2)-enemin1)/deltaene1)+1
   nn1_3 = int((eigen1_1tetra(3)-enemin1)/deltaene1)+1
   nn1_4 = int((eigen1_1tetra(4)-enemin1)/deltaene1)+1

   nn1_1 = max(1,nn1_1)
   nn1_1 = min(nn1_1,nene1)
   nn1_2 = max(1,nn1_2)
   nn1_2 = min(nn1_2,nene1)
   nn1_3 = max(1,nn1_3)
   nn1_3 = min(nn1_3,nene1)
   nn1_4 = max(1,nn1_4)
   nn1_4 = min(nn1_4,nene1)

   ! calculate Allen a_i b_i and c_i parameters
   ! sort the a_i b_i c_i
   !
   ! NOTE: indices here go from 1 to 4 instead of 0 to 3 as in Allen...
   aa(1) = epsilon2(2,1) * inv_epsilon1(2,1)
   aa(2) = epsilon2(3,1) * inv_epsilon1(3,1)
   aa(3) = epsilon2(4,1) * inv_epsilon1(4,1)
   ind_a = (/2,3,4/)
   call sort_tetra(3,aa,ind_a,tol14)
   ! aa are now in order a_s a_m a_l !!! Preserve the hash function ind_a to order the positions of k below
   delaa(1) = aa(2)-aa(1)
   delaa(2) = aa(3)-aa(1)
   delaa(3) = aa(3)-aa(2)
   inv_delaa = zero
   if(delaa(1)> tol6) inv_delaa(1)= 1.0d0 / delaa(1)
   if(delaa(2)> tol6) inv_delaa(2)= 1.0d0 / delaa(2)
   if(delaa(3)> tol6) inv_delaa(3)= 1.0d0 / delaa(3)

   bb(1) = epsilon2(1,2) * inv_epsilon1(1,2)
   bb(2) = epsilon2(3,2) * inv_epsilon1(3,2)
   bb(3) = epsilon2(4,2) * inv_epsilon1(4,2)
   ind_b = (/1,3,4/)
   call sort_tetra(3,bb,ind_b,tol14)
   delbb(1) = bb(2)-bb(1)
   delbb(2) = bb(3)-bb(1)
   delbb(3) = bb(3)-bb(2)
   inv_delbb = zero
   if(delbb(1)> tol6) inv_delbb(1)= 1.0d0 / delbb(1)
   if(delbb(2)> tol6) inv_delbb(2)= 1.0d0 / delbb(2)
   if(delbb(3)> tol6) inv_delbb(3)= 1.0d0 / delbb(3)

   cc(1) = epsilon2(1,4) * inv_epsilon1(1,4)
   cc(2) = epsilon2(2,4) * inv_epsilon1(2,4)
   cc(3) = epsilon2(3,4) * inv_epsilon1(3,4)
   ind_c = (/1,2,3/)
   call sort_tetra(3,cc,ind_c,tol14)
   delcc(1) = cc(2)-cc(1)
   delcc(2) = cc(3)-cc(1)
   delcc(3) = cc(3)-cc(2)
   inv_delcc = zero
   if(delcc(1)> tol6) inv_delcc(1)= 1.0d0 / delcc(1)
   if(delcc(2)> tol6) inv_delcc(2)= 1.0d0 / delcc(2)
   if(delcc(3)> tol6) inv_delcc(3)= 1.0d0 / delcc(3)

   !----------------------------------------------------------------------
   ! start main loop A B C over eps1
   !----------------------------------------------------------------------

   !
   !  interval enemin1 < eps1 < e1 nothing to do
   !
   !
   !  interval e1 < eps1 < e3   CASE A in Allen + first term in B
   !
   ! NB: eps1 is not updated inside the loop, only between the loops
   eps1 = enemin1+nn1_1*deltaene1
   deleps1 = eps1-eigen1_1tetra(1) ! this is Omega - omega_0
   dccde1_pre = 6.d0*volconst_mult*inv_epsilon1(2,1)*inv_epsilon1(3,1)*inv_epsilon1(4,1)

   ! note we go to nn1_3
   do ieps1=nn1_1+1,nn1_3

     dccde1 = dccde1_pre * deleps1  ! this is f_0(Omega)*6*v

     ! at fixed ieps1 we can find the pivot indices for the ieps2 loop
     nn2_1 = int((eigen2_1tetra(1)+deleps1*aa(1) -enemin2)/deltaene2)+1
     nn2_2 = int((eigen2_1tetra(1)+deleps1*aa(2) -enemin2)/deltaene2)+1
     nn2_3 = int((eigen2_1tetra(1)+deleps1*aa(3) -enemin2)/deltaene2)+1

     nn2_1 = max(1,nn2_1)
     nn2_1 = min(nn2_1,nene2)
     nn2_2 = max(1,nn2_2)
     nn2_2 = min(nn2_2,nene2)
     nn2_3 = max(1,nn2_3)
     nn2_3 = min(nn2_3,nene2)

     inv_deleps1 = 1.0d0 / deleps1

     eps2 = enemin2+nn2_1*deltaene2 ! this is E
     deleps2 = eps2 - eigen2_1tetra(1) ! this is E-epsilon_0

     !-----------------------------------------------------------------------
     ! This is case AI
     !-----------------------------------------------------------------------
     do ieps2 = nn2_1+1, nn2_2
       ! calculate running value of del "a"  = a-a_s: first term should really mix eps1 and eps2
       delaa0 = deleps2*inv_deleps1 - aa(1) ! a - a_s

       ii0 = dccde1*delaa0*inv_delaa(1)*inv_delaa(2) ! this is I_0(Omega E)

       dtweightde_tmp(1,ieps2,ieps1) = dtweightde_tmp(1,ieps2,ieps1) + &
&         ii0*(1.d0 + 0.5d0*deleps1*inv_epsilon1(ind_a(1),1)* &
&                 (-2.0d0 + delaa0*inv_delaa(1)*epsilon1(ind_a(2),ind_a(1))*inv_epsilon1(ind_a(2),1) &
&                         + delaa0*inv_delaa(2)*epsilon1(ind_a(3),ind_a(1))*inv_epsilon1(ind_a(3),1)))
       dtweightde_tmp(ind_a(1),ieps2,ieps1) = dtweightde_tmp(ind_a(1),ieps2,ieps1) + &
&         ii0*0.5d0*deleps1*inv_epsilon1(ind_a(1),1)*(2.0d0 - delaa0*inv_delaa(1) - delaa0*inv_delaa(2))
       dtweightde_tmp(ind_a(2),ieps2,ieps1) = dtweightde_tmp(ind_a(2),ieps2,ieps1) + &
&         ii0*0.5d0*delaa0*inv_delaa(1)*deleps1*inv_epsilon1(ind_a(2),1)
       dtweightde_tmp(ind_a(3),ieps2,ieps1) = dtweightde_tmp(ind_a(3),ieps2,ieps1) + &
&         ii0*0.5d0*delaa0*inv_delaa(2)*deleps1*inv_epsilon1(ind_a(3),1)
       deleps2 = deleps2 + deltaene2
     end do


     eps2 = enemin2+nn2_2*deltaene2 ! this is E
     deleps2 = eps2 - eigen2_1tetra(1)  ! E-E_0

     !-----------------------------------------------------------------------
     ! This is case AII
     !-----------------------------------------------------------------------
     do ieps2 = nn2_2+1, nn2_3
       ! calculate running value of del "a"  = a_l-a: first term should really mix eps1 and eps2
       delaa0 = aa(3) - deleps2*inv_deleps1 ! a_l - a

       ii0 = dccde1*delaa0*inv_delaa(3)*inv_delaa(2) ! this is I_0(Omega E)

       dtweightde_tmp(1,ieps2,ieps1) = dtweightde_tmp(1,ieps2,ieps1) + &
&         ii0*(1.d0 + 0.5d0*deleps1*inv_epsilon1(ind_a(3),1)* &
&                 (-2.0d0 + delaa0*inv_delaa(3)*epsilon1(ind_a(2),ind_a(3))*inv_epsilon1(ind_a(2),1) &
&                         + delaa0*inv_delaa(2)*epsilon1(ind_a(1),ind_a(3))*inv_epsilon1(ind_a(1),1)))
       dtweightde_tmp(ind_a(3),ieps2,ieps1) = dtweightde_tmp(ind_a(3),ieps2,ieps1) + &
&         ii0*0.5d0*deleps1*inv_epsilon1(ind_a(3),1)*(2.0d0 - delaa0*inv_delaa(3) - delaa0*inv_delaa(2))
       dtweightde_tmp(ind_a(2),ieps2,ieps1) = dtweightde_tmp(ind_a(2),ieps2,ieps1) + &
&         ii0*0.5d0*delaa0*inv_delaa(3)*deleps1*inv_epsilon1(ind_a(2),1)
       dtweightde_tmp(ind_a(1),ieps2,ieps1) = dtweightde_tmp(ind_a(1),ieps2,ieps1) + &
&         ii0*0.5d0*delaa0*inv_delaa(2)*deleps1*inv_epsilon1(ind_a(1),1)

       deleps2 = deleps2 + deltaene2
     end do
     deleps1 = deleps1 + deltaene1
   end do
   !
   !  interval e2 < eps < e3
   !
   eps1 = eps1 + (nn1_2-nn1_1)*deltaene1

   deleps1 = eps1-eigen1_1tetra(2) ! Omega - omega_1

   dccde1_pre = 6.d0*volconst_mult*inv_epsilon1(2,1)*inv_epsilon1(3,2)*inv_epsilon1(4,2) ! f1 function
   do ieps1=nn1_2+1,nn1_3

     dccde1 = dccde1_pre * deleps1 ! f2(Omega) * 6 * v

     ! at fixed ieps1 we can find the pivot indices for the ieps2 loop
     nn2_1 = int((eigen2_1tetra(2)+deleps1*bb(1) -enemin2)/deltaene2)+1
     nn2_2 = int((eigen2_1tetra(2)+deleps1*bb(2) -enemin2)/deltaene2)+1
     nn2_3 = int((eigen2_1tetra(2)+deleps1*bb(3) -enemin2)/deltaene2)+1

     nn2_1 = max(1,nn2_1)
     nn2_1 = min(nn2_1,nene2)
     nn2_2 = max(1,nn2_2)
     nn2_2 = min(nn2_2,nene2)
     nn2_3 = max(1,nn2_3)
     nn2_3 = min(nn2_3,nene2)

     inv_deleps1 = 1.0d0 / deleps1

     eps2 = enemin2+nn2_1*deltaene2 ! starting value for E
     deleps2 = eps2 - eigen2_1tetra(2) ! E - epsilon_1

     !-----------------------------------------------------------------------
     ! This is case BI
     !-----------------------------------------------------------------------
     do ieps2 = nn2_1+1, nn2_2
       ! calculate running value of del "b"  = b-b_s: first term should really mix eps1 and eps2
       delbb0 = deleps2*inv_deleps1 - bb(1)

       ii1 = dccde1*delbb0*inv_delbb(1)*inv_delbb(2) ! this is I_1(Omega E)

       ! note negative sign here - we are correcting the I0 a0 term already calculated above
       dtweightde_tmp(2,ieps2,ieps1) = dtweightde_tmp(2,ieps2,ieps1) - &
&         ii1*(1.d0 + 0.5d0*deleps1*inv_epsilon1(ind_b(1),2)* &
&                 (-2.0d0 + delbb0*inv_delbb(1)*epsilon1(ind_b(2),ind_b(1))*inv_epsilon1(ind_b(2),2) &
&                         + delbb0*inv_delbb(2)*epsilon1(ind_b(3),ind_b(1))*inv_epsilon1(ind_b(3),2)))
       dtweightde_tmp(ind_b(1),ieps2,ieps1) = dtweightde_tmp(ind_b(1),ieps2,ieps1) - &
&         ii1*0.5d0*deleps1*inv_epsilon1(ind_b(1),2)*(2.0d0 - delbb0*inv_delbb(1) - delbb0*inv_delbb(2))
       dtweightde_tmp(ind_b(2),ieps2,ieps1) = dtweightde_tmp(ind_b(2),ieps2,ieps1) - &
&         ii1*0.5d0*delbb0*inv_delbb(1)*deleps1*inv_epsilon1(ind_b(2),2)
       dtweightde_tmp(ind_b(3),ieps2,ieps1) = dtweightde_tmp(ind_b(3),ieps2,ieps1) - &
&         ii1*0.5d0*delbb0*inv_delbb(2)*deleps1*inv_epsilon1(ind_b(3),2)
       deleps2 = deleps2 + deltaene2
     end do

     eps2 = enemin2+nn2_2*deltaene2
     deleps2 = eps2 - eigen2_1tetra(2)

     !-----------------------------------------------------------------------
     ! This is case BII
     !-----------------------------------------------------------------------
     do ieps2 = nn2_2+1, nn2_3
       ! calculate running value of del "b"  = b_l-b: first term should really mix eps1 and eps2
       delbb0 = bb(3) - deleps2*inv_deleps1

       ii1 = dccde1*delbb0*inv_delbb(3)*inv_delbb(2) ! this is I_1(Omega E)

       ! note negative sign here - we are correcting the I0 a0 term already calculated above
       dtweightde_tmp(2,ieps2,ieps1) = dtweightde_tmp(2,ieps2,ieps1) - &
&         ii1*(1.d0 + 0.5d0*deleps1*inv_epsilon1(ind_b(3),2)* &
&                 (-2.0d0 + delbb0*inv_delbb(3)*epsilon1(ind_b(2),ind_b(3))*inv_epsilon1(ind_b(2),2) &
&                         + delbb0*inv_delbb(2)*epsilon1(ind_b(1),ind_b(3))*inv_epsilon1(ind_b(1),2)))
       dtweightde_tmp(ind_b(3),ieps2,ieps1) = dtweightde_tmp(ind_b(3),ieps2,ieps1) - &
&         ii1*0.5d0*deleps1*inv_epsilon1(ind_b(3),2)*(2.0d0 - delbb0*inv_delbb(3) - delbb0*inv_delbb(2))
       dtweightde_tmp(ind_b(2),ieps2,ieps1) = dtweightde_tmp(ind_b(2),ieps2,ieps1) - &
&         ii1*0.5d0*delbb0*inv_delbb(3)*deleps1*inv_epsilon1(ind_b(2),2)
       dtweightde_tmp(ind_b(1),ieps2,ieps1) = dtweightde_tmp(ind_b(1),ieps2,ieps1) - &
&         ii1*0.5d0*delbb0*inv_delbb(2)*deleps1*inv_epsilon1(ind_b(1),2)

       deleps2 = deleps2 + deltaene2
     end do

     deleps1 = deleps1 + deltaene1
   end do

   !
   !  interval e3 < eps < e4
   !
   eps1 = eps1 + (nn1_3-nn1_2)*deltaene1
   deleps1 = eps1-eigen1_1tetra(4)
   dccde1_pre = 6.d0*volconst_mult*inv_epsilon1(4,1)*inv_epsilon1(4,2)*inv_epsilon1(4,3)
   do ieps1=nn1_3+1,nn1_4
     ! note - sign from definition of f3
     dccde1 = -dccde1_pre *  deleps1 ! f3(Omega) * 6 * v

     ! at fixed ieps1 we can find the pivot indices for the ieps2 loop
     ! NB: order is inverted for cc because deleps1 is defined negative (Omega is always less than omega_3)
     nn2_1 = int((eigen2_1tetra(4)+deleps1*cc(3) -enemin2)/deltaene2)+1
     nn2_2 = int((eigen2_1tetra(4)+deleps1*cc(2) -enemin2)/deltaene2)+1
     nn2_3 = int((eigen2_1tetra(4)+deleps1*cc(1) -enemin2)/deltaene2)+1

     nn2_1 = max(1,nn2_1)
     nn2_1 = min(nn2_1,nene2)
     nn2_2 = max(1,nn2_2)
     nn2_2 = min(nn2_2,nene2)
     nn2_3 = max(1,nn2_3)
     nn2_3 = min(nn2_3,nene2)
     inv_deleps1 = 1.0d0 / deleps1

     eps2 = enemin2+nn2_1*deltaene2 ! starting value for E
     deleps2 = eps2 - eigen2_1tetra(4) ! E - epsilon_3

     !-----------------------------------------------------------------------
     ! This is case CII
     !-----------------------------------------------------------------------
     do ieps2 = nn2_1+1, nn2_2
       ! calculate running value of del "c"  = c_l-c: first term should really mix eps1 and eps2
       delcc0 = cc(3) - deleps2*inv_deleps1

       ii3 = dccde1*delcc0*inv_delcc(3)*inv_delcc(2) ! this is I_3(Omega E)

       dtweightde_tmp(4,ieps2,ieps1) = dtweightde_tmp(4,ieps2,ieps1) + &
&         ii3*(1.d0 + 0.5d0*deleps1*inv_epsilon1(ind_c(3),4)* &
&                 (-2.0d0 + delcc0*inv_delcc(3)*epsilon1(ind_c(2),ind_c(3))*inv_epsilon1(ind_c(2),4) &
&                         + delcc0*inv_delcc(2)*epsilon1(ind_c(1),ind_c(3))*inv_epsilon1(ind_c(1),4)))
       dtweightde_tmp(ind_c(3),ieps2,ieps1) = dtweightde_tmp(ind_c(3),ieps2,ieps1) + &
&         ii3*0.5d0*deleps1*inv_epsilon1(ind_c(3),4)*(2.0d0 - delcc0*inv_delcc(3) - delcc0*inv_delcc(2))
       dtweightde_tmp(ind_c(2),ieps2,ieps1) = dtweightde_tmp(ind_c(2),ieps2,ieps1) + &
&         ii3*0.5d0*delcc0*inv_delcc(3)*deleps1*inv_epsilon1(ind_c(2),4)
       dtweightde_tmp(ind_c(1),ieps2,ieps1) = dtweightde_tmp(ind_c(1),ieps2,ieps1) + &
&         ii3*0.5d0*delcc0*inv_delcc(2)*deleps1*inv_epsilon1(ind_c(1),4)

       deleps2 = deleps2 + deltaene2
     end do


     eps2 = enemin2+nn2_2*deltaene2
     deleps2 = eps2 - eigen2_1tetra(4)

     !-----------------------------------------------------------------------
     ! This is case CI
     !-----------------------------------------------------------------------
     do ieps2 = nn2_2+1, nn2_3
       ! calculate running value of del "c"  = c-c_s: first term should really mix eps1 and eps2
       delcc0 = deleps2*inv_deleps1 - cc(1) ! c - c_s

       ii3 = dccde1*delcc0*inv_delcc(1)*inv_delcc(2) ! this is I_3(Omega E)

       dtweightde_tmp(4,ieps2,ieps1) = dtweightde_tmp(4,ieps2,ieps1) + &
&         ii3*(1.d0 + 0.5d0*deleps1*inv_epsilon1(ind_c(1),4)* &
&                 (-2.0d0 + delcc0*inv_delcc(1)*epsilon1(ind_c(2),ind_c(1))*inv_epsilon1(ind_c(2),4) &
&                         + delcc0*inv_delcc(2)*epsilon1(ind_c(3),ind_c(1))*inv_epsilon1(ind_c(3),4)))
       dtweightde_tmp(ind_c(1),ieps2,ieps1) = dtweightde_tmp(ind_c(1),ieps2,ieps1) + &
&         ii3*0.5d0*deleps1*inv_epsilon1(ind_c(1),4)*(2.0d0 - delcc0*inv_delcc(1) - delcc0*inv_delcc(2))
       dtweightde_tmp(ind_c(2),ieps2,ieps1) = dtweightde_tmp(ind_c(2),ieps2,ieps1) + &
&         ii3*0.5d0*delcc0*inv_delcc(1)*deleps1*inv_epsilon1(ind_c(2),4)
       dtweightde_tmp(ind_c(3),ieps2,ieps1) = dtweightde_tmp(ind_c(3),ieps2,ieps1) + &
&         ii3*0.5d0*delcc0*inv_delcc(2)*deleps1*inv_epsilon1(ind_c(3),4)
       deleps2 = deleps2 + deltaene2
     end do

     deleps1 = deleps1 + deltaene1
   end do

   eps1 = eps1 + (nn1_4-nn1_3)*deltaene1
   !
   !
   !  interval e4 < eps < enemax
   !
   do ieps1=nn1_4+1,nene1
     ! dtweightde unchanged by this tetrahedron
   end do

   ! if we have a fully degenerate tetrahedron,
   ! 1) the tweight is a Heaviside (step) function, which is correct above, but
   ! 2) the dtweightde should contain a Dirac function: add a Gaussian here

   ! TODO: add treatment in double tetra case
   !  end degenerate tetrahedron if

   ! the following blas calls are not working systematically, or do not give speed ups, strange...
   !call daxpy (nene, 1.d0, dtweightde_tmp(:,1), 1, dtweightde_t(:,ind_ibz(1)), 1)
   !call daxpy (nene, 1.d0, dtweightde_tmp(:,2), 1, dtweightde_t(:,ind_ibz(2)), 1)
   !call daxpy (nene, 1.d0, dtweightde_tmp(:,3), 1, dtweightde_t(:,ind_ibz(3)), 1)
   !call daxpy (nene, 1.d0, dtweightde_tmp(:,4), 1, dtweightde_t(:,ind_ibz(4)), 1)

   do ieps2 = 1, nene2
     dtweightde(ind_k(1),:,ieps2) = dtweightde(ind_k(1),:,ieps2) + dtweightde_tmp(1,ieps2,:)
     dtweightde(ind_k(2),:,ieps2) = dtweightde(ind_k(2),:,ieps2) + dtweightde_tmp(2,ieps2,:)
     dtweightde(ind_k(3),:,ieps2) = dtweightde(ind_k(3),:,ieps2) + dtweightde_tmp(3,ieps2,:)
     dtweightde(ind_k(4),:,ieps2) = dtweightde(ind_k(4),:,ieps2) + dtweightde_tmp(4,ieps2,:)
     !tweight(nkpt,nene1,nene2)
   end do

 end do ! itetra

 ! transpose: otherwise the data access is crap and the code slows by an order of magnitude
 TETRA_DEALLOCATE(tweight_tmp)
 TETRA_DEALLOCATE(dtweightde_tmp)

end subroutine get_dbl_tetra_weight
!!***

!!****f* m_tetrahedron/sort_tetra
!! NAME
!!  sort_tetra
!!
!! FUNCTION
!!  Sort double precision array list(n) into ascending numerical order using Heapsort
!!  algorithm, while making corresponding rearrangement of the integer
!!  array iperm. Consider that two double precision numbers
!!  within tolerance tol are equal.
!!
!! INPUTS
!!  n        intent(in)    dimension of the list
!!  tol      intent(in)    numbers within tolerance are equal
!!  list(n)  intent(inout) list of double precision numbers to be sorted
!!  iperm(n) intent(inout) iperm(i)=i (very important)
!!
!! OUTPUT
!!  list(n)  sorted list
!!  iperm(n) index of permutation given the right ascending order
!!
!! PARENTS
!!      m_tetrahedron
!!
!! CHILDREN
!!      get_onetetra_,sort_tetra
!!
!! SOURCE


subroutine sort_tetra(n,list,iperm,tol)

 integer, intent(in) :: n
 integer, intent(inout) :: iperm(n)
 real(dp), intent(inout) :: list(n)
 real(dp), intent(in) :: tol

 integer :: l,ir,iap,i,j
 real(dp) :: ap
 character(len=500) :: msg

 if (n==1) then
   ! Accomodate case of array of length 1: already sorted!
   return
 else if (n<1) then
  ! Should not call with n<1
  write(msg,1000) n
  1000  format(/,' sort_tetra has been called with array length n=',i12,/, &
&  ' having a value less than 1. This is not allowed.')
  TETRA_ERROR(msg)

 else ! n>1

  ! Conduct the usual sort
  l=n/2+1
  ir=n

  do   ! Infinite do-loop
   if (l>1) then
    l=l-1
    ap=list(l)
    iap=iperm(l)

   else ! l<=1
    ap=list(ir)
    iap=iperm(ir)
    list(ir)=list(1)
    iperm(ir)=iperm(1)
    ir=ir-1

    if (ir==1) then
     list(1)=ap
     iperm(1)=iap
     exit   ! This is the end of this algorithm
    end if
   end if ! l>1

   i=l
   j=l+l

   do while (j<=ir)
    if (j<ir) then
     if ( list(j)<list(j+1)-tol .or.  &
&        (list(j)<list(j+1)+tol.and.iperm(j)<iperm(j+1))) j=j+1
    endif
    if (ap<list(j)-tol .or. (ap<list(j)+tol.and.iap<iperm(j))) then
     list(i)=list(j)
     iperm(i)=iperm(j)
     i=j
     j=j+j
    else
     j=ir+1
    end if
   enddo

   list(i)=ap
   iperm(i)=iap

  enddo ! End infinite do-loop

 end if ! n>1

end subroutine sort_tetra
!!***

!----------------------------------------------------------------------

!!****f* m_tetrahedron/tetralib_has_mpi
!! NAME
!! tetralib_has_mpi
!!
!! FUNCTION
!! Return True if library has been compiled with MPI support
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

logical function tetralib_has_mpi() result(ans)

  ans = .False.
#ifdef HAVE_MPI
  ans = .True.
#endif

end function tetralib_has_mpi
!!***

!----------------------------------------------------------------------

!!****f* m_tetrahedron/split_work
!! NAME
!!  split_work
!!
!! FUNCTION
!!  Splits the number of tasks, ntasks, among nprocs processors. Used for the MPI parallelization of simple loops.
!!
!! INPUTS
!!  ntasks=number of tasks
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  nprocs=Number of MPI processes in the communicator.
!!  my_start,my_stop= indices defining the initial and final task for this processor
!!  ierr=Exit status.
!!
!! NOTES
!!  If nprocs>ntasks then :
!!    my_start=ntasks+1
!!    my_stop=ntask
!!
!!  In this particular case, loops of the form
!!
!!  do ii=my_start,my_stop
!!   ...
!!  end do
!!
!!  are not executed. Moreover allocation such as foo(my_start:my_stop) will generate a zero-sized array.
!!
!! PARENTS
!!      m_tetrahedron
!!
!! CHILDREN
!!      get_onetetra_,sort_tetra
!!
!! SOURCE

subroutine split_work(ntasks,comm,nprocs,my_start,my_stop,ierr)

!Arguments ------------------------------------
 integer,intent(in)  :: ntasks,comm
 integer,intent(out) :: nprocs,my_start,my_stop,ierr

!Local variables-------------------------------
 integer :: res,my_rank,block_p1,block,mpierr

! *************************************************************************

 nprocs = 1; my_start = 1; my_stop = ntasks; ierr = 1
#ifdef HAVE_MPI
 call MPI_COMM_SIZE(comm,nprocs,mpierr); if (mpierr /= MPI_SUCCESS) return
 call MPI_COMM_RANK(comm,my_rank,mpierr); if (mpierr /= MPI_SUCCESS) return

 block   = ntasks/nprocs
 res     = MOD(ntasks,nprocs)
 block_p1= block+1

 if (my_rank<res) then
   my_start =  my_rank   *block_p1+1
   my_stop  = (my_rank+1)*block_p1
 else
   my_start = res*block_p1 + (my_rank-res  )*block + 1
   my_stop  = res*block_p1 + (my_rank-res+1)*block
 end if
#endif
 ierr = 0

end subroutine split_work
!!***

!----------------------------------------------------------------------

!!****f* m_tetrahedron/get_onetetra_
!! NAME
!! get_onetetra_
!!
!! FUNCTION
!! Private function to calculate the contributions to the weights due to a single tetrahedron.
!! Extracted from get_tetra_weight
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure subroutine get_onetetra_(tetra,itetra,eigen_1tetra,enemin,enemax,max_occ,nene,bcorr, &
&  tweight_tmp,dtweightde_tmp)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nene,bcorr,itetra
 type(t_tetrahedron), intent(in) :: tetra
 real(dp) ,intent(in) :: enemax,enemin,max_occ
!arrays
 ! MGTODO: This layout is not optimal (lots of cache thrashing, I will optimize it later on)
 real(dp), intent(out) ::  tweight_tmp(nene, 4)
 real(dp), intent(out) :: dtweightde_tmp(nene, 4)
 real(dp),intent(in)  :: eigen_1tetra(4)

!Local variables-------------------------------
!  needed for gaussian replacement of Dirac functions
!  the three coefficients of the DOS as quadratic form,
!    in the interval [eig(ikpt-1), eig(ikpt)]
!    for ikpt = 1 we add a point below eigen(1) which doesnt
!    contribute to the DOS in any tetrahedron
!scalars
 integer :: ieps,nn1,nn2,nn3,nn4
 real(dp)  :: cc,cc1,cc2,cc3,dcc1de,dcc2de,dcc3de,dccde,deltaene,eps
 real(dp)  :: epsilon21,epsilon31,epsilon32,epsilon41,epsilon42,epsilon43
 real(dp)  :: gau_prefactor,gau_width,gau_width2,inv_epsilon21,inv_epsilon31,gval
 real(dp)  :: inv_epsilon32,inv_epsilon41,inv_epsilon42,inv_epsilon43
 real(dp)  :: deleps1, deleps2, deleps3, deleps4
 real(dp)  :: invepsum, cc_pre, dccde_pre
 real(dp)  :: cc1_pre, cc2_pre, cc3_pre
 real(dp)  :: cc_tmp, dccde_tmp
 real(dp)  :: dcc1de_pre, dcc2de_pre, dcc3de_pre
 real(dp)  :: tmp,volconst,volconst_mult

! *********************************************************************

 volconst = tetra%vv/4.d0

 deltaene = (enemax-enemin) / (nene-1)

 ! This is output
 tweight_tmp = zero; dtweightde_tmp = zero

 volconst_mult = max_occ*volconst*dble(tetra%tetra_mult(itetra))

 ! all notations are from Blochl PRB 49 16223 [[cite:Bloechl1994a]] Appendix B
 epsilon21 = eigen_1tetra(2)-eigen_1tetra(1)
 epsilon31 = eigen_1tetra(3)-eigen_1tetra(1)
 epsilon41 = eigen_1tetra(4)-eigen_1tetra(1)
 epsilon32 = eigen_1tetra(3)-eigen_1tetra(2)
 epsilon42 = eigen_1tetra(4)-eigen_1tetra(2)
 epsilon43 = eigen_1tetra(4)-eigen_1tetra(3)
 inv_epsilon21 = zero; if (epsilon21 > tol6) inv_epsilon21 = 1.d0 / epsilon21
 inv_epsilon31 = zero; if (epsilon31 > tol6) inv_epsilon31 = 1.d0 / epsilon31
 inv_epsilon41 = zero; if (epsilon41 > tol6) inv_epsilon41 = 1.d0 / epsilon41
 inv_epsilon32 = zero; if (epsilon32 > tol6) inv_epsilon32 = 1.d0 / epsilon32
 inv_epsilon42 = zero; if (epsilon42 > tol6) inv_epsilon42 = 1.d0 / epsilon42
 inv_epsilon43 = zero; if (epsilon43 > tol6) inv_epsilon43 = 1.d0 / epsilon43

 nn1 = int((eigen_1tetra(1)-enemin)/deltaene)+1
 nn2 = int((eigen_1tetra(2)-enemin)/deltaene)+1
 nn3 = int((eigen_1tetra(3)-enemin)/deltaene)+1
 nn4 = int((eigen_1tetra(4)-enemin)/deltaene)+1

 nn1 = max(1,nn1)
 nn1 = min(nn1,nene)
 nn2 = max(1,nn2)
 nn2 = min(nn2,nene)
 nn3 = max(1,nn3)
 nn3 = min(nn3,nene)
 nn4 = max(1,nn4)
 nn4 = min(nn4,nene)

 eps = enemin+nn1*deltaene
 !
 !interval enemin < eps < e1 nothing to do
 !
 !
 !interval e1 < eps < e2
 !
 deleps1 = eps-eigen_1tetra(1)
 cc_pre = volconst_mult*inv_epsilon21*inv_epsilon31*inv_epsilon41
 invepsum = inv_epsilon21+inv_epsilon31+inv_epsilon41
 dccde_pre = 3.d0*volconst_mult*inv_epsilon21*inv_epsilon31*inv_epsilon41
 do ieps=nn1+1,nn2
   cc = cc_pre * deleps1*deleps1*deleps1
   tweight_tmp(ieps,1) = tweight_tmp(ieps,1) + cc*(4.d0-deleps1*invepsum)
   tweight_tmp(ieps,2) = tweight_tmp(ieps,2) + cc*deleps1*inv_epsilon21
   tweight_tmp(ieps,3) = tweight_tmp(ieps,3) + cc*deleps1*inv_epsilon31
   tweight_tmp(ieps,4) = tweight_tmp(ieps,4) + cc*deleps1*inv_epsilon41

   dccde = dccde_pre * deleps1*deleps1
   dtweightde_tmp(ieps,1) = dtweightde_tmp(ieps,1) + dccde*(4.d0 - deleps1*invepsum) -cc*invepsum
   dtweightde_tmp(ieps,2) = dtweightde_tmp(ieps,2) + (dccde*deleps1 + cc) * inv_epsilon21
   dtweightde_tmp(ieps,3) = dtweightde_tmp(ieps,3) + (dccde*deleps1 + cc) * inv_epsilon31
   dtweightde_tmp(ieps,4) = dtweightde_tmp(ieps,4) + (dccde*deleps1 + cc) * inv_epsilon41

   if (bcorr == 1) then
     ! bxu, correction terms based on Bloechl's paper
     tweight_tmp(ieps,1) = tweight_tmp(ieps,1) + &
&     4.d0*dccde_pre*deleps1*deleps1*(epsilon21+epsilon31+epsilon41)/40.d0
     tweight_tmp(ieps,2) = tweight_tmp(ieps,2) + &
&     4.d0*dccde_pre*deleps1*deleps1*(-epsilon21+epsilon32+epsilon42)/40.d0
     tweight_tmp(ieps,3) = tweight_tmp(ieps,3) + &
&     4.d0*dccde_pre*deleps1*deleps1*(-epsilon31-epsilon32+epsilon43)/40.d0
     tweight_tmp(ieps,4) = tweight_tmp(ieps,4) + &
&     4.d0*dccde_pre*deleps1*deleps1*(-epsilon41-epsilon42-epsilon43)/40.d0

     dtweightde_tmp(ieps,1) = dtweightde_tmp(ieps,1) + &
&     8.d0*dccde_pre*deleps1*(epsilon21+epsilon31+epsilon41)/40.d0
     dtweightde_tmp(ieps,2) = dtweightde_tmp(ieps,2) + &
&     8.d0*dccde_pre*deleps1*(-epsilon21+epsilon32+epsilon42)/40.d0
     dtweightde_tmp(ieps,3) = dtweightde_tmp(ieps,3) + &
&     8.d0*dccde_pre*deleps1*(-epsilon31-epsilon32+epsilon43)/40.d0
     dtweightde_tmp(ieps,4) = dtweightde_tmp(ieps,4) + &
&     8.d0*dccde_pre*deleps1*(-epsilon41-epsilon42-epsilon43)/40.d0
   end if

   deleps1 = deleps1 + deltaene
 end do

 eps = eps + (nn2-nn1)*deltaene
 !
 !  interval e2 < eps < e3
 !
 deleps1 = eps-eigen_1tetra(1)
 deleps2 = eps-eigen_1tetra(2)
 deleps3 = eigen_1tetra(3)-eps
 deleps4 = eigen_1tetra(4)-eps

 cc1_pre = volconst_mult*inv_epsilon31*inv_epsilon41
 cc2_pre = volconst_mult*inv_epsilon41*inv_epsilon32*inv_epsilon31
 cc3_pre = volconst_mult*inv_epsilon42*inv_epsilon32*inv_epsilon41

 dcc1de_pre = 2.d0*cc1_pre
 dcc2de_pre = cc2_pre
 dcc3de_pre = cc3_pre
 do ieps=nn2+1,nn3
   cc1 = cc1_pre * deleps1*deleps1
   cc2 = cc2_pre * deleps1*deleps2*deleps3
   cc3 = cc3_pre * deleps2*deleps2*deleps4

   tweight_tmp(ieps,1) = tweight_tmp(ieps,1) + &
&   cc1 + (cc1+cc2)*deleps3*inv_epsilon31 + (cc1+cc2+cc3)*deleps4*inv_epsilon41
   tweight_tmp(ieps,2) = tweight_tmp(ieps,2) + &
&   cc1+cc2+cc3+(cc2+cc3)*deleps3*inv_epsilon32 + cc3*deleps4*inv_epsilon42
   tweight_tmp(ieps,3) = tweight_tmp(ieps,3) + &
&   (cc1+cc2)*deleps1*inv_epsilon31 + (cc2+cc3)*deleps2*inv_epsilon32
   tweight_tmp(ieps,4) = tweight_tmp(ieps,4) + &
&   (cc1+cc2+cc3)*deleps1*inv_epsilon41 + cc3*deleps2*inv_epsilon42


   dcc1de = dcc1de_pre * deleps1
   dcc2de = dcc2de_pre * (-deleps1*deleps2  +deleps1*deleps3  +deleps2*deleps3)
   dcc3de = dcc3de_pre * (2.d0*deleps2*deleps4  -deleps2*deleps2)

   dtweightde_tmp(ieps,1) = dtweightde_tmp(ieps,1) &
&   + dcc1de &
&   + ((dcc1de+dcc2de)*deleps3 -(cc1+cc2)) * inv_epsilon31 &
&   + ((dcc1de+dcc2de+dcc3de)*deleps4 -(cc1+cc2+cc3)) * inv_epsilon41
   dtweightde_tmp(ieps,2) = dtweightde_tmp(ieps,2) &
&   + dcc1de+dcc2de+dcc3de &
&   + ((dcc2de+dcc3de)*deleps3 -(cc2+cc3) ) * inv_epsilon32 &
&   + (dcc3de*deleps4  -cc3 ) * inv_epsilon42
   dtweightde_tmp(ieps,3) = dtweightde_tmp(ieps,3) &
&   + ((dcc1de+dcc2de)*deleps1 + (cc1+cc2) ) * inv_epsilon31 &
&   + ((dcc2de+dcc3de)*deleps2 + (cc2+cc3) ) * inv_epsilon32
   dtweightde_tmp(ieps,4) = dtweightde_tmp(ieps,4) &
&   + ((dcc1de+dcc2de+dcc3de)*deleps1 + (cc1+cc2+cc3) ) * inv_epsilon41 &
&   + (dcc3de*deleps2 + cc3) * inv_epsilon42

 if (bcorr == 1) then
   ! bxu, correction terms based on Bloechl's paper
   ! The correction terms may cause the dtweightde become negative
   tweight_tmp(ieps,1) = tweight_tmp(ieps,1) + &
&   4.d0*cc1_pre* &
&   (3.d0*epsilon21+6.d0*deleps2-3.d0*(epsilon31+epsilon42)*deleps2**2.d0*inv_epsilon32*inv_epsilon42)* &
&   (epsilon21+epsilon31+epsilon41)/40.d0
   tweight_tmp(ieps,2) = tweight_tmp(ieps,2) + &
&   4.d0*cc1_pre* &
&   (3.d0*epsilon21+6.d0*deleps2-3.d0*(epsilon31+epsilon42)*deleps2**2.d0*inv_epsilon32*inv_epsilon42)* &
&   (-epsilon21+epsilon32+epsilon42)/40.d0
   tweight_tmp(ieps,3) = tweight_tmp(ieps,3) + &
&   4.d0*cc1_pre* &
&   (3.d0*epsilon21+6.d0*deleps2-3.d0*(epsilon31+epsilon42)*deleps2**2.d0*inv_epsilon32*inv_epsilon42)* &
&   (-epsilon31-epsilon32+epsilon43)/40.d0
   tweight_tmp(ieps,4) = tweight_tmp(ieps,4) + &
&   4.d0*cc1_pre* &
&   (3.d0*epsilon21+6.d0*deleps2-3.d0*(epsilon31+epsilon42)*deleps2**2.d0*inv_epsilon32*inv_epsilon42)* &
&   (-epsilon41-epsilon42-epsilon43)/40.d0

   dtweightde_tmp(ieps,1) = dtweightde_tmp(ieps,1) + &
&   4.d0*cc1_pre* &
&   (6.d0-6.d0*(epsilon31+epsilon42)*deleps2*inv_epsilon32*inv_epsilon42)* &
&   (epsilon21+epsilon31+epsilon41)/40.d0
   dtweightde_tmp(ieps,2) = dtweightde_tmp(ieps,2) + &
&   4.d0*cc1_pre* &
&   (6.d0-6.d0*(epsilon31+epsilon42)*deleps2*inv_epsilon32*inv_epsilon42)* &
&   (-epsilon21+epsilon32+epsilon42)/40.d0
   dtweightde_tmp(ieps,3) = dtweightde_tmp(ieps,3) + &
&   4.d0*cc1_pre* &
&   (6.d0-6.d0*(epsilon31+epsilon42)*deleps2*inv_epsilon32*inv_epsilon42)* &
&   (-epsilon31-epsilon32+epsilon43)/40.d0
   dtweightde_tmp(ieps,4) = dtweightde_tmp(ieps,4) + &
&   4.d0*cc1_pre* &
&   (6.d0-6.d0*(epsilon31+epsilon42)*deleps2*inv_epsilon32*inv_epsilon42)* &
&   (-epsilon41-epsilon42-epsilon43)/40.d0
  end if

  deleps1 = deleps1 + deltaene
  deleps2 = deleps2 + deltaene
  deleps3 = deleps3 - deltaene
  deleps4 = deleps4 - deltaene
 end do

 eps = eps + (nn3-nn2)*deltaene
 !
 !  interval e3 < eps < e4
 !
 deleps4 = eigen_1tetra(4)-eps
 cc_pre = volconst_mult*inv_epsilon41*inv_epsilon42*inv_epsilon43
 invepsum = inv_epsilon41+inv_epsilon42+inv_epsilon43
 dccde_pre = -3.d0*cc_pre
 do ieps=nn3+1,nn4
   cc = cc_pre * deleps4*deleps4*deleps4
   cc_tmp = cc * deleps4
   tweight_tmp(ieps,1) = tweight_tmp(ieps,1) + volconst_mult - cc_tmp*inv_epsilon41
   tweight_tmp(ieps,2) = tweight_tmp(ieps,2) + volconst_mult - cc_tmp*inv_epsilon42
   tweight_tmp(ieps,3) = tweight_tmp(ieps,3) + volconst_mult - cc_tmp*inv_epsilon43
   tweight_tmp(ieps,4) = tweight_tmp(ieps,4) + volconst_mult - cc*4.d0 + cc_tmp*invepsum

   dccde = dccde_pre * deleps4*deleps4
   dccde_tmp = -dccde*deleps4 + cc
   dtweightde_tmp(ieps,1) = dtweightde_tmp(ieps,1) + dccde_tmp * inv_epsilon41
   dtweightde_tmp(ieps,2) = dtweightde_tmp(ieps,2) + dccde_tmp * inv_epsilon42
   dtweightde_tmp(ieps,3) = dtweightde_tmp(ieps,3) + dccde_tmp * inv_epsilon43
   dtweightde_tmp(ieps,4) = dtweightde_tmp(ieps,4) - dccde*4.d0 - dccde_tmp*invepsum

   if (bcorr == 1) then
     ! bxu, correction terms based on Bloechl's paper
     ! The correction terms may cause the dtweightde become negative
     tweight_tmp(ieps,1) = tweight_tmp(ieps,1) + &
&     12.d0*cc_pre*deleps4*deleps4*(epsilon21+epsilon31+epsilon41)/40.d0
     tweight_tmp(ieps,2) = tweight_tmp(ieps,2) + &
&     12.d0*cc_pre*deleps4*deleps4*(-epsilon21+epsilon32+epsilon42)/40.d0
     tweight_tmp(ieps,3) = tweight_tmp(ieps,3) + &
&     12.d0*cc_pre*deleps4*deleps4*(-epsilon31-epsilon32+epsilon43)/40.d0
     tweight_tmp(ieps,4) = tweight_tmp(ieps,4) + &
&     12.d0*cc_pre*deleps4*deleps4*(-epsilon41-epsilon42-epsilon43)/40.d0

     dtweightde_tmp(ieps,1) = dtweightde_tmp(ieps,1) - &
&     24.d0*cc_pre*deleps4*(epsilon21+epsilon31+epsilon41)/40.d0
     dtweightde_tmp(ieps,2) = dtweightde_tmp(ieps,2) - &
&     24.d0*cc_pre*deleps4*(-epsilon21+epsilon32+epsilon42)/40.d0
     dtweightde_tmp(ieps,3) = dtweightde_tmp(ieps,3) - &
&     24.d0*cc_pre*deleps4*(-epsilon31-epsilon32+epsilon43)/40.d0
     dtweightde_tmp(ieps,4) = dtweightde_tmp(ieps,4) - &
&     24.d0*cc_pre*deleps4*(-epsilon41-epsilon42-epsilon43)/40.d0
   end if

   deleps4 = deleps4 - deltaene
 end do
 eps = eps + (nn4-nn3)*deltaene
 !
 !
 !  interval e4 < eps < enemax
 !
 do ieps=nn4+1,nene
   tweight_tmp(ieps,1) = tweight_tmp(ieps,1) + volconst_mult
   tweight_tmp(ieps,2) = tweight_tmp(ieps,2) + volconst_mult
   tweight_tmp(ieps,3) = tweight_tmp(ieps,3) + volconst_mult
   tweight_tmp(ieps,4) = tweight_tmp(ieps,4) + volconst_mult
   ! dtweightde unchanged by this tetrahedron
 end do

 !
 !  if we have a fully degenerate tetrahedron,
 !  1) the tweight is a Heaviside (step) function, which is correct above, but
 !  2) the dtweightde should contain a Dirac function: add a Gaussian here
 !
 if (epsilon41 < tol6) then

   !  to ensure the gaussian will integrate properly:
   !  WARNING: this smearing could be problematic if too large
   !  and doesnt integrate well if its too small
   gau_width = 10.0d0*deltaene
   gau_width2 = 1.0 / gau_width / gau_width
   gau_prefactor = volconst_mult / gau_width / sqrtpi
   !
   ! average position since bracket for epsilon41 is relatively large
   cc = (eigen_1tetra(1)+eigen_1tetra(2)+eigen_1tetra(3)+eigen_1tetra(4))/4.d0
   eps = enemin
   do ieps=1,nene
     tmp = eps - cc
     gval = gau_prefactor*exp(-tmp*tmp*gau_width2)
     ! MG TODO: I think this is not correct, because we have divided by 4 so
     ! the other points should be accumulated as well.
     ! There are however changes in the unit tests if I activate these lines...
     !dtweightde_tmp(ieps,1) = dtweightde_tmp(ieps,1) + gval
     !dtweightde_tmp(ieps,2) = dtweightde_tmp(ieps,2) + gval
     !dtweightde_tmp(ieps,3) = dtweightde_tmp(ieps,3) + gval
     dtweightde_tmp(ieps,4) = dtweightde_tmp(ieps,4) + gval
     eps = eps + deltaene
   end do
 end if ! end degenerate tetrahedron if

end subroutine get_onetetra_
!!***

!----------------------------------------------------------------------

!!****f* m_tetrahedron/tetra_get_onewk
!! NAME
!! tetra_get_onewk
!!
!! FUNCTION
!! Calculate integration weights and their derivatives for a single k-point in the IBZ.
!!
!! INPUTS
!! tetra<t_tetrahedron>=Object with tables for tetrahedron method.
!! ik_ibz=Index of the k-point in the IBZ array
!! bcorr=1 to include Blochl correction else 0.
!! nene=number of energies for DOS
!! nibz=number of irreducible kpoints
!! eigen_ibz(nkibz)=eigenenergies for each k point
!! enemin=minimal energy for DOS
!! enemax=maximal energy for DOS
!! max_occ=maximal occupation number (2 for nsppol=1, 1 for nsppol=2)
!!
!! OUTPUT
!!  weights(nene,2) = integration weights for
!!    Dirac delta (derivative of theta wrt energy) and Theta (Heaviside function)
!!    for a given (band, k-point, spin).
!!
!! PARENTS
!!
!! CHILDREN
!!      get_onetetra_,sort_tetra
!!
!! SOURCE

subroutine tetra_get_onewk(tetra,ik_ibz,bcorr,nene,nkibz,eig_ibz,enemin,enemax,max_occ,weights)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,nene,nkibz,bcorr
 type(t_tetrahedron), intent(in) :: tetra
 real(dp) ,intent(in) :: enemin,enemax,max_occ
!arrays
 real(dp),intent(in) :: eig_ibz(nkibz)
 real(dp),intent(out) :: weights(nene,2)

!Local variables-------------------------------
!scalars
 integer :: itetra,ii
!arrays
 integer :: ind_ibz(4)
 real(dp) :: tweight_tmp(nene,4),dtweightde_tmp(nene,4),eigen_1tetra(4)

! *********************************************************************

 weights = zero

 ! For each tetrahedron
 do itetra=1,tetra%ntetra

   ! Here we need the original ordering to reference the correct irred kpoints
   ind_ibz(:) = tetra%tetra_full(:,1,itetra)
   ! Cycle if this tetra does not contribute to this k-point.
   if (all(ind_ibz /= ik_ibz)) cycle

   ! Sort energies before calling get_onetetra_
   eigen_1tetra(:) = eig_ibz(ind_ibz(:))
   call sort_tetra(4, eigen_1tetra, ind_ibz, tol14)

   call get_onetetra_(tetra, itetra, eigen_1tetra, enemin, enemax, max_occ, nene, bcorr, &
     tweight_tmp, dtweightde_tmp)

   ! Accumulate contributions to ik_ibz (there might be multiple vertexes that map onto ik_ibz)
   do ii=1,4
     if (ind_ibz(ii) == ik_ibz) then
       weights(:,1) = weights(:,1) + dtweightde_tmp(:,ii)
       weights(:,2) = weights(:,2) + tweight_tmp(:,ii)
     end if
   end do
 end do ! itetra

end subroutine tetra_get_onewk
!!***

!----------------------------------------------------------------------

!!****f* m_tetrahedron/tetra_get_onewk_wvals
!! NAME
!! tetra_get_onewk_wvals
!!
!! FUNCTION
!! Calculate integration weights and their derivatives for a single k-point in the IBZ.
!!
!! INPUTS
!! tetra<t_tetrahedron>=Object with tables for tetrahedron method.
!! ik_ibz=Index of the k-point in the IBZ array
!! bcorr=1 to include Blochl correction else 0.
!! nw=number of energies in wvals
!! nibz=number of irreducible kpoints
!! wvals(nw)=Frequency points.
!! eigen_ibz(nkibz)=eigenenergies for each k point
!! [wtol]: If present, frequency points that differ by less that wtol are treated as equivalent.
!!  and the tetrahedron integration is performed only once per frequency point.
!!
!! OUTPUT
!!  weights(nw,2) = integration weights for
!!    Dirac delta (derivative of theta wrt energy) and Theta (Heaviside function)
!!    for a given (band, k-point, spin).
!!
!! PARENTS
!!
!! CHILDREN
!!      get_onetetra_,sort_tetra
!!
!! SOURCE

subroutine tetra_get_onewk_wvals(tetra, ik_ibz, bcorr, nw, wvals, nkibz, eig_ibz, weights, wtol)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,nw,nkibz,bcorr
 real(dp), optional, intent(in) :: wtol
 type(t_tetrahedron), intent(in) :: tetra
!arrays
 real(dp),intent(in) :: wvals(nw)
 real(dp),intent(in) :: eig_ibz(nkibz)
 real(dp),intent(out) :: weights(nw, 2)

!Local variables-------------------------------
!scalars
 !integer,save :: done = 0
 integer,parameter :: nene=3
 integer :: itetra,ii,jj,iw,ie
 logical :: samew
 real(dp),parameter :: max_occ1 = one
 real(dp) :: enemin, enemax
!arrays
 integer :: ind_ibz(4)
 real(dp) :: theta_tmp(nene,4), delta_tmp(nene,4), eigen_1tetra(4)

! *********************************************************************

 weights = zero

 ! For each tetrahedron
 do jj=1,tetra%ibz_tetra_count(ik_ibz)
   itetra = tetra%ibz_tetra_mapping(ik_ibz,jj)

   ! Here we need the original ordering to reference the correct irred kpoints
   ind_ibz(:) = tetra%tetra_full(:,1,itetra)

   ! Sort energies before calling get_onetetra_
   eigen_1tetra(:) = eig_ibz(ind_ibz(:))
   call sort_tetra(4, eigen_1tetra, ind_ibz, tol14)

   do iw=1,nw
     samew = .False.
     if (present(wtol)) then
       if (iw > 1) samew = abs(wvals(iw) - wvals(iw - 1)) < wtol
     end if
     if (.not. samew) then
         enemin = wvals(iw) - 0.01; enemax = wvals(iw) + 0.01
         ie = nene / 2 + 1
         call get_onetetra_(tetra, itetra, eigen_1tetra, enemin, enemax, max_occ1, nene, bcorr, &
            theta_tmp, delta_tmp)
     end if

     ! Accumulate contributions to ik_ibz (there might be multiple vertexes that map onto ik_ibz)
     do ii=1,4
       if (ind_ibz(ii) == ik_ibz) then
         weights(iw, 1) = weights(iw, 1) + delta_tmp(ie, ii)
         weights(iw, 2) = weights(iw, 2) + theta_tmp(ie, ii)
       end if
     end do
   end do ! iw
 end do ! itetra

end subroutine tetra_get_onewk_wvals
!!***

!----------------------------------------------------------------------

!!****f* m_tetrahedron/tetra_get_onetetra_wvals
!! NAME
!! tetra_get_onetetra_wvals
!!
!! FUNCTION
!! Calculate integration weights and their derivatives for a single k-point in the IBZ.
!!
!! INPUTS
!! tetra<t_tetrahedron>=Object with tables for tetrahedron method.
!! ik_ibz=Index of the k-point in the IBZ array
!! bcorr=1 to include Blochl correction else 0.
!! nw=number of energies in wvals
!! nibz=number of irreducible kpoints
!! wvals(nw)=Frequency points.
!! eigen_ibz(nkibz)=eigenenergies for each k point
!! [wtol]: If present, frequency points that differ by less that wtol are treated as equivalent.
!!  and the tetrahedron integration is performed only once per frequency point.
!!
!! OUTPUT
!!  weights(nw,2) = integration weights for
!!    Dirac delta (derivative of theta wrt energy) and Theta (Heaviside function)
!!    for a given (band, k-point, spin).
!!
!! PARENTS
!!
!! CHILDREN
!!      get_onetetra_,sort_tetra
!!
!! SOURCE

subroutine tetra_get_onetetra_wvals(tetra, itetra, eigen_1tetra, bcorr, nw, wvals, weights, wtol)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nw,bcorr
 real(dp), optional, intent(in) :: wtol
 type(t_tetrahedron), intent(in) :: tetra
!arrays
 real(dp),intent(in) :: wvals(nw)
 real(dp),intent(out) :: weights(nw, 2, 4)

!Local variables-------------------------------
!scalars
 !integer,save :: done = 0
 integer,parameter :: nene3=3
 integer :: itetra,ii,idx,iw,ie
 integer :: ind(4)
 logical :: samew
 real(dp),parameter :: max_occ1 = one
 real(dp) :: enemin, enemax
!arrays
 real(dp) :: theta_tmp(nene3,4), delta_tmp(nene3,4), eigen_1tetra(4)

! *********************************************************************

 ind = [1,2,3,4]
 call sort_tetra(4, eigen_1tetra, ind, tol14)
 weights = 0

 !for all the frequencies
 do iw=1,nw
   samew = .False.
   if (present(wtol)) then
     if (iw > 1) samew = abs(wvals(iw) - wvals(iw - 1)) < wtol
   end if
   if (.not. samew) then
     enemin = wvals(iw) - 0.01
     enemax = wvals(iw) + 0.01
     ie = nene3 / 2 + 1
     call get_onetetra_(tetra, itetra, eigen_1tetra, enemin, enemax, max_occ1, nene3, bcorr, &
        theta_tmp, delta_tmp)
   end if

   ! Accumulate contributions to ik_ibz (there might be multiple vertexes that map onto ik_ibz)
   do ii=1,4
     idx = ind(ii)
     weights(iw, 1, idx) = weights(iw, 1, idx) + delta_tmp(ie, ii)
     weights(iw, 2, idx) = weights(iw, 2, idx) + theta_tmp(ie, ii)
   end do
 end do !iw

end subroutine tetra_get_onetetra_wvals
!!***

end module m_tetrahedron
!!***
