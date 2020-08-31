!!****m* ABINIT/m_supercell
!!
!! NAME
!! m_supercell
!!
!! FUNCTION
!! Module for using a supercell, in particular for phonon displacement freezing.
!! Container type is defined, and destruction, print subroutines
!! as well as the central init_supercell
!!
!! COPYRIGHT
!! Copyright (C) 2010-2020 ABINIT group (MJV, DJA)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_supercell

 use defs_basis
 use m_errors
 use m_abicore

 use m_symtk,    only : matr3inv
 use m_copy,     only : alloc_copy
 use m_io_tools, only : open_file
 use m_fstrings, only : int2char4, write_num

 implicit none

 private

!!***

!!****t* m_supercell/supercell_type
!! NAME
!! supercell_type
!!
!! FUNCTION
!! structure for a supercell constructed from a basic rprimd and xcart, with indexing to original atoms
!! The supercell may not be oriented the same way as the original cell, if you can reduce it by symmetry
!!
!! SOURCE

 type, public :: supercell_type
   integer :: natom_primcell                        ! number of atoms in primitive cell
   integer :: natom                                 ! number of atoms in supercell
   integer :: ntypat                                ! number of atom types
   integer :: ncells                                ! number of unit cells in supercell
   integer :: rlatt(3,3)                            ! matrix for multiplicity of supercell
   real(dp) :: rprimd(3,3)                          ! new lattice vectors for supercell
   real(dp) :: qphon(3)                             ! phonon q vector used to generate scell, if any
   real(dp), allocatable :: xcart(:,:)              ! (3, natom) positions of atoms
   real(dp), allocatable :: xcart_ref(:,:)          ! (3, natom) equilibrium positions of atoms
   integer, allocatable :: atom_indexing(:)         ! (natom) indexes original atom: 1..natom_primcell
   integer, allocatable :: uc_indexing(:,:)         ! (3, natom) indexes unit cell atom is in:
   integer, allocatable :: typat(:)                 ! (3, natom) positions of atoms
   real(dp), allocatable :: znucl(:)                ! (ntypat) nuclear charges of species
   integer, allocatable :: rvecs(:,:)               ! supercell vectors

 end type supercell_type

 public :: init_supercell_for_qpt
 public :: init_supercell
 public :: freeze_displ_supercell
 public :: prt_supercell_for_qpt
 public :: prt_supercell
 public :: copy_supercell
 public :: distance_supercell
 public :: findBound_supercell
 public :: getPBCIndexes_supercell
 public :: destroy_supercell

 public :: mksupercell  !  computes atomic positons, magnetic ordering of supercell
!!***

CONTAINS  !===========================================================================================

!!****f* m_supercell/init_supercell_for_qpt
!!
!! NAME
!! init_supercell_for_qpt
!!
!! FUNCTION
!! Initialize scell structure, from unit cell vectors, and atoms, based on qpoint chosen
!!
!! INPUTS
!! natom_primcell = number of atoms in primitive cell
!! qphon(3) = phonon wavevector
!!      find smallest supercell which will accomodate phonon qphon = (1/2,1/2,1/2)
!! rprimd_primcell(3,3) = real space lattice vectors (bohr)
!! typat_primcell = types of atoms
!! xcart_primcell(3,natom) = cartesian positions of atoms in primitive cell
!! znucl = nuclear charges for all species
!! ordering = if true,  typat will be 1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 ....
!!            if false, typat will be 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 ....
!!
!! OUTPUT
!! scell = supercell structure to be initialized
!!
!! PARENTS
!!      m_phonons
!!
!! CHILDREN
!!
!! SOURCE

subroutine init_supercell_for_qpt(natom_primcell, qphon, rprimd_primcell, &
&    typat_primcell, xcart_primcell, znucl, scell, ordering)

 implicit none
!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom_primcell
 type(supercell_type), intent(out) :: scell
 logical,optional,intent(in) :: ordering
!arrays
 integer , intent(in) :: typat_primcell(natom_primcell)
 real(dp), intent(in) :: qphon(3)
 real(dp), intent(in) :: znucl(:)
 real(dp), intent(in) :: rprimd_primcell(3,3)
 real(dp), intent(in) :: xcart_primcell(3,natom_primcell)

!Local variables-------------------------------
!scalar
 integer :: ii, maxsc, iscmult
 real(dp) :: qbymult
!arrays
 integer :: rlatt(3,3) ! number of primitive cells in each direction for the supercell
 character(len=500) :: msg
! *************************************************************************


! maximum number of unit cells in a given direction
 maxsc = 10

! find smallest supercell which will accomodate phonon.
! FIXME: for the moment, just get smallest multiple along each direction, with an upper bound
 rlatt = 0
 rlatt(1,1) = -1
 rlatt(2,2) = -1
 rlatt(3,3) = -1
 do ii=1,3
   do iscmult=1,maxsc
     qbymult = qphon(ii)*iscmult
     if (abs(qbymult - int(qbymult)) < tol10) then
       rlatt(ii,ii) = iscmult
       exit
     end if
   end do
   if (rlatt(ii,ii) == -1) then
     write(msg,'(a,I4,a,I7,2a,3E20.10)')' No supercell found with less than ', &
&           maxsc,' unit cells in direction ', &
&           ii, ch10, ' qphon = ', qphon
     MSG_ERROR(msg)
   end if
 end do

 if (present(ordering)) then
   call init_supercell(natom_primcell, rlatt, rprimd_primcell, typat_primcell, xcart_primcell, znucl, scell, ordering)
 else
   call init_supercell(natom_primcell, rlatt, rprimd_primcell, typat_primcell, xcart_primcell, znucl, scell)
 end if

 scell%qphon = qphon

end subroutine init_supercell_for_qpt
!!***

!!****f* m_supercell/init_supercell
!!
!! NAME
!! init_supercell
!!
!! FUNCTION
!! Initialize scell structure, from unit cell vectors, and atoms, based on rlatt multiplicity matrix
!!
!! INPUTS
!! natom_primcell = number of atoms in primitive cell
!! rlatt(3,3) = multiplicity of primtive unit cells in supercell
!! rprimd_primcell(3,3) = real space lattice vectors (bohr)
!! typat_primcell(natom) = types of all atoms in primitive cell
!! xcart_primcell(3,natom) = cartesian positions of atoms in primitive cell
!! znucl = nuclear charges for all species
!!
!! ordering = if true,  typat will be 1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 ....
!!            if false, typat will be 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 ....
!!
!! OUTPUT
!! scell = supercell structure to be initialized
!!
!! PARENTS
!!      m_effective_potential,m_fit_polynomial_coeff,m_phonons,m_supercell
!!
!! CHILDREN
!!
!! SOURCE

subroutine init_supercell(natom_primcell, rlatt, rprimd_primcell, typat_primcell, xcart_primcell, znucl, scell, ordering)

 implicit none
!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom_primcell
 type(supercell_type), intent(out) :: scell
 logical,optional,intent(in) :: ordering
!arrays
 integer , intent(in) :: rlatt(3,3)
 integer , intent(in) :: typat_primcell(natom_primcell)
 real(dp), intent(in) :: znucl(:)
 real(dp), intent(in) :: rprimd_primcell(3,3)
 real(dp), intent(in) :: xcart_primcell(3,natom_primcell)

!local
!scalars
 integer :: iatom_supercell, i1,i2,i3, iatom, icell

!arrays

 scell%natom_primcell = natom_primcell
 scell%rlatt = rlatt
 scell%ncells = rlatt(1,1)*rlatt(2,2)*rlatt(3,3)
 scell%rprimd(:,1) = rprimd_primcell(:,1) * rlatt(1,1)
 scell%rprimd(:,2) = rprimd_primcell(:,2) * rlatt(2,2)
 scell%rprimd(:,3) = rprimd_primcell(:,3) * rlatt(3,3)

 scell%ntypat = size(znucl)
 ABI_ALLOCATE(scell%znucl,(scell%ntypat))
 scell%znucl(:) = znucl(:)

!number of atoms in full supercell
 scell%natom= natom_primcell*scell%ncells
 ABI_ALLOCATE(scell%xcart,(3,scell%natom))
 ABI_ALLOCATE(scell%xcart_ref,(3,scell%natom))
 ABI_ALLOCATE(scell%typat,(scell%natom))
 ABI_ALLOCATE(scell%atom_indexing,(scell%natom))
 ABI_ALLOCATE(scell%uc_indexing,(3,scell%natom))
 ABI_ALLOCATE(scell%rvecs, (3, scell%ncells))

 iatom_supercell = 0
 icell =0
 do i1 = 1, rlatt(1,1)
   do i2 = 1, rlatt(2,2)
     do i3 = 1, rlatt(3,3)

       icell=icell+1
       scell%rvecs(:,icell)=(/i1-1, i2-1, i3-1/)


       do iatom = 1, natom_primcell
         iatom_supercell = iatom_supercell + 1
         scell%uc_indexing(:,iatom_supercell) = (/i1-1,i2-1,i3-1/)
         scell%xcart_ref(:,iatom_supercell) = xcart_primcell(:,iatom) &
&            + matmul(rprimd_primcell,scell%uc_indexing(:,iatom_supercell))
         scell%atom_indexing(iatom_supercell) = iatom
         scell%typat(iatom_supercell) = typat_primcell(iatom)
       end do
     end do
   end do
 end do

 ABI_CHECK(iatom_supercell == scell%natom, "iatom_supercell /= scell%natom")
 if(iatom_supercell /= scell%natom) then
    write(std_out,*)"iatom_supercell /= scell%natom"
 endif



 scell%xcart = scell%xcart_ref
 scell%qphon = zero

 if (present(ordering)) then
   if (ordering) then
     call order_supercell_typat(scell)
   end if
 end if

end subroutine init_supercell
!!***

!!****f* m_supercell/order_supercell_typat
!!
!! NAME
!! order_supercell_typat
!!
!! FUNCTION
!! Re-order atoms in place for types
!!
!! INPUTS
!! scell = supercell structure with reference atomic positions etc...
!!
!! OUTPUT
!! scell = supercell structure: typat, xcart and so on will be updated
!!
!! PARENTS
!!      m_supercell
!!
!! CHILDREN
!!
!! SOURCE

subroutine order_supercell_typat (scell)

 implicit none

!Arguments ------------------------------------
!scalars
 type(supercell_type), intent(inout) :: scell

! local tmp variables
 integer :: itypat, iatom_supercell, iatom
 type(supercell_type) :: scell_tmp

 call copy_supercell (scell,scell_tmp)

 iatom_supercell = 0
 do itypat = 1, scell%ntypat
   do iatom = 1, scell%natom
     if (scell_tmp%typat(iatom) /= itypat) cycle
     iatom_supercell = iatom_supercell + 1
     scell%xcart(:,iatom_supercell) = scell_tmp%xcart(:,iatom)
     scell%xcart_ref(:,iatom_supercell) = scell_tmp%xcart_ref(:,iatom)
     scell%atom_indexing(iatom_supercell) = scell_tmp%atom_indexing(iatom)
     scell%uc_indexing(:,iatom_supercell) = scell_tmp%uc_indexing(:,iatom)
     scell%typat(iatom_supercell) = scell_tmp%typat(iatom)
   end do
 end do

 call destroy_supercell(scell_tmp)

end subroutine order_supercell_typat
!!***


!!****f* m_supercell/freeze_displ_supercell
!!
!! NAME
!! freeze_displ_supercell
!!
!! FUNCTION
!! Freeze a specific displacement phonon field into the supercell scell
!!
!! INPUTS
!! displ = phonon displacement vectors for this mode
!! freeze_displ = desired amplitude for phonon displacement along displ.
!!    for thermal displacement use sqrt[ (1/2 + bose_einstein(freq,T)) / freq ]
!! scell = supercell structure with reference atomic positions etc...
!!
!! OUTPUT
!! scell = supercell structure: xcart will be updated with phonon displacement
!!
!! PARENTS
!!      m_phonons
!!
!! CHILDREN
!!
!! SOURCE

subroutine freeze_displ_supercell (displ,freeze_displ,scell)

 implicit none

!Arguments ------------------------------------
!scalars
 type(supercell_type), intent(inout) :: scell
 real(dp), intent(in) :: freeze_displ
!arrays
 real(dp), intent(in) :: displ(2,3*scell%natom_primcell)

!Local variables-------------------------------
!scalar
 integer :: iatom, ipratom
 complex(dpc) :: expqdotr, j=cmplx(zero,one)
 complex(dpc) :: phase
 complex(dpc) :: zdispl(3,scell%natom_primcell)
! *************************************************************************

 zdispl = (cmplx(reshape(displ(1,:), (/3,scell%natom_primcell/)),&
&                reshape(displ(2,:), (/3,scell%natom_primcell/))))

 ! fix gauge by imposing real displacement for first atom in first direction
 ! multiply by normalized complex conjugate of first element
 ! NB 6 March 2018: this may be imposing a positive (not just real) displacement for 1st atom along x!!!
 ! That might be problematic below, though for the thermal displacement method freeze_displ swaps sign for each new mode
 phase = cmplx(one,zero)
 if (abs(zdispl(1,1)) > tol10) then
   phase = conjg(zdispl(1,1)) / abs(zdispl(1,1))
 end if

 do iatom = 1, scell%natom
   expqdotr = exp(j*two_pi*(scell%qphon(1)*scell%uc_indexing(1,iatom) &
&                          +scell%qphon(2)*scell%uc_indexing(2,iatom) &
&                          +scell%qphon(3)*scell%uc_indexing(3,iatom)))

! this is offset in zdispl vector due to primitive cell atom position
   ipratom = scell%atom_indexing(iatom)

!add real part of displacement times Bloch phase
   scell%xcart(:,iatom) = scell%xcart(:,iatom) &
&        + freeze_displ * real(expqdotr * zdispl(:,ipratom) * phase)

!   scell%xcart(:,iatom) = scell%xcart(:,iatom) &
!&        + freeze_displ * cos(qdotr) * displ(1,ipratom+1:ipratom+3) &
!&        - freeze_displ * sin(qdotr) * displ(2,ipratom+1:ipratom+3)
 end do

end subroutine freeze_displ_supercell
!!***

!****f* m_supercell/prt_supercell_for_qpt
!!
!! NAME
!! prt_supercell_for_qpt
!!
!! FUNCTION
!! output atomic positions, supercell vectors, etc... to a file
!! single qpoint and mode
!!
!! INPUTS
!! freq = phonon frequency for mode jmode
!! jmode = mode which has been frozen into xcart contained in scell
!! outfile_radix = radix of file name to be written to
!! scell = supercell structure with data to be output
!!
!! OUTPUT
!! printing to file
!!
!! PARENTS
!!      m_phonons
!!
!! CHILDREN
!!
!! SOURCE

subroutine prt_supercell_for_qpt (freq, jmode, outfile_radix, scell)

  implicit none

!Arguments ------------------------------------
!scalars
  real(dp), intent(in) :: freq
  type(supercell_type), intent(in) :: scell
  integer, intent(in) :: jmode
  character(len=fnlen), intent(in) :: outfile_radix

!Local variables-------------------------------
!scalar
  character(len=fnlen) :: filename
  character(len=10) :: jmodestring
  character(len=80) :: title1, title2
  character(len=5) :: qphonstring1, qphonstring2, qphonstring3

! add suffix with mode and qpoint
  call int2char4(jmode, jmodestring)
  ABI_CHECK((jmodestring(1:1)/='#'),'Bug: string length too short!')
! qphonstring should be like 0.000_0.000_0.000
  call write_num(scell%qphon(1),qphonstring1,'(F5.3)')
  call write_num(scell%qphon(2),qphonstring2,'(F5.3)')
  call write_num(scell%qphon(3),qphonstring3,'(F5.3)')
  filename = trim(outfile_radix) // "_qpt_" // qphonstring1 // "_" // qphonstring2 // &
&              "_" // qphonstring3 // "_mode_" // trim(jmodestring)

  write (title1, '(a,3E20.10)') '# phonon q point : ', scell%qphon
  write (title2, '(a,I7,a,E20.10)') '# phonon mode number : ', jmode, ' frequency ', freq

  call prt_supercell(filename, scell, title1, title2)

end subroutine prt_supercell_for_qpt
!!***


!****f* m_supercell/prt_supercell
!!
!! NAME
!! prt_supercell
!!
!! FUNCTION
!! output atomic positions, supercell vectors, etc... to a file
!!
!! INPUTS
!! filename = filename
!! title1 = first line of description of contents
!! title2 = second line of description of contents
!! scell = supercell structure with data to be output
!!
!! OUTPUT
!! printing to file
!!
!! PARENTS
!!      m_phonons,m_supercell
!!
!! CHILDREN
!!
!! SOURCE

subroutine prt_supercell (filename, scell, title1, title2)

  implicit none

!Arguments ------------------------------------
!scalars
  type(supercell_type), intent(in) :: scell
  character(len=fnlen), intent(in) :: filename
  character(len=80), intent(in) :: title1
  character(len=80), intent(in) :: title2


!Local variables-------------------------------
!scalar
  integer :: scunit, iatom
  character(len=500) :: msg
  real(dp) :: xred(3), gprimd(3,3)

! *************************************************************************

  if (open_file(filename, msg, newunit=scunit) /= 0) then
    MSG_ERROR(msg)
  end if

! print header
  write (scunit, '(a)') '#'
  write (scunit, '(a)') '# anaddb file with frozen phonon mode in supercell'
  write (scunit, '(a)') '# !!!   Do not forget to adjust nband   !!! '
  write (scunit, '(a)') '#'
  write (scunit, '(a)') title1
  write (scunit, '(a)') title2
  write (scunit, '(a,3(3I7,2x))') '# supercell rlatt is ', scell%rlatt
  write (scunit, '(a,I7,a)') '# and has ', scell%ncells, ' primitive unit cells '
  write (scunit, '(a)') '#'
  write (scunit, '(a)') '# lattice vectors for supercell :'
  write (scunit, '(a,I7)') 'natom ', scell%natom
  write (scunit, *)
  write (scunit, '(a)') 'znucl '
  do iatom = 1, size(scell%znucl)
    write (scunit, '(I5)', ADVANCE="NO") int(scell%znucl(iatom))
    if (mod(iatom,6) == 0) write (scunit, *)
  end do
  write (scunit, *)
  write (scunit, *)
  write (scunit, '(a,I7)') 'ntypat', scell%ntypat
  write (scunit, '(a)') 'typat '
  do iatom = 1, scell%natom
    write (scunit, '(I5)', ADVANCE="NO") scell%typat(iatom)
    if (mod(iatom,6) == 0) write (scunit, *)
  end do
  write (scunit, *)
  write (scunit, '(a)') 'acell 1.0 1.0 1.0'
  write (scunit, '(a)') 'rprim'
  write (scunit, '(3E20.10)') scell%rprimd(:,1)
  write (scunit, '(3E20.10)') scell%rprimd(:,2)
  write (scunit, '(3E20.10)') scell%rprimd(:,3)
  write (scunit, *)
  write (scunit, '(a)') 'xcart'
  do iatom = 1, scell%natom
    write (scunit, '(3E20.10)') scell%xcart(:,iatom)
  end do
  ! for information, also print xred for atoms inside full supercell
  call matr3inv(scell%rprimd, gprimd)
  ! TODO: check this transpose is correct in some asymetric case
  gprimd = transpose(gprimd)
  write (scunit, '(a)') '# for information, add xred as well'
  write (scunit, '(a)') '# xred'
  do iatom = 1, scell%natom
    xred = matmul (gprimd, scell%xcart(:,iatom))
    write (scunit, '(a, 3E20.10)') '#  ', xred
  end do

! close file
  close(scunit)

end subroutine prt_supercell
!!***

!****f* m_supercell/copy_supercell
!!
!! NAME
!! copy_supercell
!!
!! FUNCTION
!! copy supercell structure
!!
!! INPUTS
!! scell_in = supercell structure with data to copy
!!
!! OUTPUT
!! scell = supercell structure with data to be output
!!
!! PARENTS
!!      m_effective_potential,m_supercell
!!
!! CHILDREN
!!
!! SOURCE

subroutine copy_supercell (scell_in,scell_copy)

  implicit none

!Arguments ------------------------------------
!scalars
  type(supercell_type), intent(in) :: scell_in
  type(supercell_type), intent(inout) :: scell_copy

! *************************************************************************

  call destroy_supercell(scell_copy)
  scell_copy%natom_primcell = scell_in%natom_primcell
  scell_copy%natom = scell_in%natom
  scell_copy%ntypat = scell_in%ntypat
  scell_copy%ncells = scell_in%ncells
  scell_copy%rlatt = scell_in%rlatt
  scell_copy%rprimd = scell_in%rprimd
  scell_copy%qphon = scell_in%qphon
  call alloc_copy(scell_in%xcart        , scell_copy%xcart)
  call alloc_copy(scell_in%xcart_ref    , scell_copy%xcart_ref)
  call alloc_copy(scell_in%atom_indexing, scell_copy%atom_indexing)
  call alloc_copy(scell_in%uc_indexing  , scell_copy%uc_indexing)
  call alloc_copy(scell_in%typat        , scell_copy%typat)
  call alloc_copy(scell_in%znucl        , scell_copy%znucl)

end subroutine copy_supercell
!!***

!!****f* m_effective_potential/getPBCIndexes_supercell
!! NAME
!!
!! FUNCTION
!! Get the index of the cell by using PBC
!!
!! INPUTS
!! index  = index of the cell into the supercell
!! ncell = number of total cell
!!
!! OUTPUT
!! index  = index of the cell into the supercell with PBC
!!
!! SOURCE

subroutine getPBCIndexes_supercell(index,ncell)

 implicit none

!Arguments ---------------------------------------------
  integer, intent(inout)  :: index(3)
  integer, intent(in) :: ncell(3)
!Local variables ---------------------------------------
  integer :: ii
! *********************************************************************

  do ii=1,3
    do while (index(ii) > ncell(ii))
      index(ii) = index(ii) - ncell(ii)
    end do
    do while (index(ii) <= 0)
      index(ii) = index(ii) + ncell(ii)
    end do
  end do

end subroutine getPBCIndexes_supercell
!!***

!****f* m_supercell/findBound_supercell
!! NAME
!!  findBound_supercell
!!
!! FUNCTION
!!  compute the bound of the supercell by considering the 0 0 0 (reference)
!!  in the center of the supercell.
!!  for example: (4 4 4) => min = -1 and max = 2
!!
!! INPUTS
!! ncell(3) = size of the supercell (for example 3 3 3)
!!
!! OUTPUT
!! min = minimun of the range
!! max = maximum of the range
!!
!! PARENTS
!!      m_effective_potential,m_polynomial_coeff
!!
!! CHILDREN
!!
!! SOURCE

subroutine findBound_supercell(min,max,ncell)

 implicit none

!Arguments ---------------------------------------------
  integer, intent(inout) :: min,max
  integer, intent(in) :: ncell
!Local variables ---------------------------------------
  if(abs(max)>abs(min)) then
    max=(ncell)/2; min=-max;  if(mod(ncell,2)==0) max = max -1
  else
    min=-(ncell)/2; max=-min; if(mod(ncell,2)==0)  min= min +1
  end if

! *********************************************************************
end subroutine findBound_supercell
!!***

!!****f* m_supercell/distance_supercell
!! NAME
!!
!! FUNCTION
!! compute the distance_supercell betwen 2 atoms in different cell
!!
!! INPUTS
!! xcart1(3) = cartesian coordinates of the first atom
!! xcart1(3) = cartesian coordinates of the second atom
!! rprimd(3,3) = primitive lattice vectors
!! cell1(3) = index of the cell of the first atom (for example -1 0 2)
!! cell2(3) = index of the cell of the second atom (for example  0 0 2)
!!
!! OUTPUT
!! distance_supercell = distance_supercell between the 2 atoms
!!
!! SOURCE
!!

function distance_supercell(xcart1,xcart2,rprimd,cell1,cell2) result(dist)

!Arguments ------------------------------------
!scalar
!array
  real(dp),intent(in):: rprimd(3,3)
  real(dp),intent(in):: xcart1(3),xcart2(3)
  integer,intent(in) :: cell1(3),cell2(3)
  real(dp) :: dist
!Local variables -------------------------------
  real(dp) :: rpt1(3),rpt2(3)
  integer  :: mu
!! *************************************************************************
  do mu=1,3
    rpt1(mu) = cell1(1)*rprimd(mu,1)+cell1(2)*rprimd(mu,2)+cell1(3)*rprimd(mu,3)
    rpt2(mu) = cell2(1)*rprimd(mu,1)+cell2(2)*rprimd(mu,2)+cell2(3)*rprimd(mu,3)
  end do

  dist = ((xcart2(1)+rpt2(1)-xcart1(1)-rpt1(1))**2+&
&         (xcart2(2)+rpt2(2)-xcart1(2)-rpt1(2))**2+&
&         (xcart2(3)+rpt2(3)-xcart1(3)-rpt1(3))**2)**0.5

end function distance_supercell
!!***

!****f* m_supercell/destroy_supercell
!!
!! NAME
!! destroy_supercell
!!
!! FUNCTION
!! deallocate all dynamic memory for this supercell structure
!!
!! INPUTS
!!
!! OUTPUT
!! scell = supercell structure with data to be output
!!
!! PARENTS
!!      m_effective_potential,m_fit_polynomial_coeff,m_phonons,m_supercell
!!
!! CHILDREN
!!
!! SOURCE

subroutine destroy_supercell (scell)

  implicit none

!Arguments ------------------------------------
!scalars
  type(supercell_type), intent(inout) :: scell

! *************************************************************************

  if(allocated(scell%xcart))  then
    ABI_DEALLOCATE(scell%xcart)
  end if
  if(allocated(scell%xcart_ref))  then
    ABI_DEALLOCATE(scell%xcart_ref)
  end if
  if(allocated(scell%typat))  then
    ABI_DEALLOCATE(scell%typat)
  end if
  if(allocated(scell%atom_indexing))  then
    ABI_DEALLOCATE(scell%atom_indexing)
  end if
  if(allocated(scell%uc_indexing))  then
    ABI_DEALLOCATE(scell%uc_indexing)
  end if
  if(allocated(scell%znucl))  then
    ABI_DEALLOCATE(scell%znucl)
  end if
   if(allocated(scell%rvecs))  then
    ABI_DEALLOCATE(scell%rvecs)
  end if

end subroutine destroy_supercell
!!***

!!****f* m_supercell/mksupercell
!! NAME
!!  mksupercell
!!
!! FUNCTION
!!  computes atomic positons, magnetic ordering of supercell
!!
!! INPUTS
!!  magv_org (optional) magnetic ordering of atoms in primitive cell,
!!   ordering of atoms given als 1 and -1, if not given fm is assumed
!!  xred_org relative position of atoms in primitive cell
!!  rprimd_org unit cell dimensions of primitive cell
!!  natom=number of atoms in unit cell
!!  option= 1 output ion-ion distances / 2 output ordering of ion-ion distances / 3 output variables in varlist
!!           according to ion-ion distances * magnetic ordering
!!
!! OUTPUT
!!  magv_sc magnetic ordering of atoms in supercell
!!  xred_sc relative position of atoms in supercell
!!  rprimd_sc unit cell dimensions of supercell
!!
!! PARENTS
!!      m_paw_uj
!!
!! CHILDREN
!!
!! SOURCE

subroutine mksupercell(xred_org,magv_org,rprimd_org,nat_org,nat_sc,xred_sc,magv_sc,rprimd_sc,ext,prtvol)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)              :: nat_org,nat_sc
 integer,intent(in),optional     :: prtvol
!arrays
 real(dp),intent(in)             :: rprimd_org(3,3)
 integer,intent(in)              :: ext(3)
 real(dp),intent(in)             :: xred_org(3,nat_org)
 real(dp),intent(out)            :: xred_sc(3,nat_sc)
 real(dp),intent(out)            :: magv_sc(nat_sc)
 real(dp),intent(out)            :: rprimd_sc(3,3)
 integer,intent(in),optional     :: magv_org(nat_org)

!Local variables-------------------------------
!scalars
 integer :: prtvoll,ix,iy,iz,nprcl,iprcl,jdim,iatom
!arrays
 real(dp) :: magvv_org(nat_org)
 real(dp),allocatable :: transv(:,:,:)

! *************************************************************************

 if (present(magv_org)) then
   magvv_org=magv_org
 else
   magvv_org=(/ (1, iatom=1,nat_org)  /)
 end if

 if (present(prtvol)) then
   prtvoll=prtvol
 else
   prtvoll=1
 end if

 rprimd_sc=reshape((/ (rprimd_org(ix,:)*ext(ix) ,ix=1,3) /),(/3,3 /))
 nprcl=product(ext)
 ABI_ALLOCATE(transv,(3,nat_org,nprcl))

 transv=reshape((/ (((((/ ix,iy,iz /),iatom=1,nat_org),ix=0,ext(1)-1),iy=0,ext(2)-1),iz=0,ext(3)-1) /), (/ 3, nat_org,nprcl/) )

!DEBUG
!write(std_out,*)'mksupercell: xred_org ' ,xred_org
!END DEBUG

 do iprcl=1,nprcl
   xred_sc(:,1+(iprcl-1)*nat_org:iprcl*nat_org)=xred_org+transv(:,:,iprcl)
   magv_sc(1+(iprcl-1)*nat_org:iprcl*nat_org)=magv_org
 end do


 do jdim=1,3
   xred_sc(jdim,:)=xred_sc(jdim,:)/ext(jdim)
 end do

!DEBUG
!write(std_out,*)'mksupercell: xred_sc ', xred_sc
!write(std_out,*)'mksupercell: magv_sc ', magv_sc
!END DEBUG

 ABI_DEALLOCATE(transv)

end subroutine mksupercell
!!***

end module m_supercell
!!***
