!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_phonon_supercell
!!
!! NAME
!! m_phonon_supercell
!!
!! FUNCTION
!! Module for using a phonon supercell
!! Container type is defined, and destruction, print subroutines 
!! as well as the central init_supercell
!!
!! COPYRIGHT
!! Copyright (C) 2010-2016 ABINIT group (MJV)
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

module m_phonon_supercell

 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_io_tools, only : open_file
 use m_fstrings, only : int2char4, write_num

 implicit none

 private
!!***

!!***

!!****t* defs_abitypes/supercell_type
!! NAME
!! supercell_type
!!
!! FUNCTION
!! structure for a supercell constructed from a basic rprimd and xcart, with indexing to original atoms
!! The supercell may not be oriented the same way as the original cell, if you can reduce it by symmetry
!!
!! SOURCE

 type, public :: supercell_type
   integer :: natom                                 ! number of atoms in primitive cell
   integer :: natom_supercell                       ! number of atoms in supercell
   real(dp) :: rprimd_supercell(3,3)                ! new lattice vectors for supercell
   real(dp) :: qphon(3)                             ! phonon q vector
   real(dp), allocatable :: xcart_supercell(:,:)        ! (3, natom_supercell) positions of atoms
   real(dp), allocatable :: xcart_supercell_ref(:,:)    ! (3, natom_supercell) equilibrium positions of atoms
   integer, allocatable :: atom_indexing_supercell(:)   ! (natom_supercell) indexes original atom: 1..natom 
   integer, allocatable :: uc_indexing_supercell(:,:)   ! (3, natom_supercell) indexes unit cell atom is in:
 end type supercell_type

 public :: init_supercell
 public :: freeze_displ_supercell
 public :: prt_supercell
 public :: destroy_supercell
!!***

CONTAINS  !===========================================================================================

!!****f* m_phonon_supercell/init_supercell
!!
!! NAME
!! init_supercell
!!
!! FUNCTION
!! Initialize scell structure, from unit cell vectors, qpoint chosen, and atoms
!!
!! INPUTS
!! natom = number of atoms in primitive cell
!! qphon(3) = phonon wavevector
!! rprimd(3,3) = real space lattice vectors (bohr)
!! xcart(3,natom) = cartesian positions of atoms in primitive cell
!!
!! OUTPUT
!! scell = supercell structure to be initialized
!!
!! PARENTS
!!      freeze_displ_allmodes
!!
!! CHILDREN
!!
!! SOURCE

subroutine init_supercell(natom, qphon, rprimd, xcart, scell)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_supercell'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom
 type(supercell_type), intent(out) :: scell
!arrays
 real(dp), intent(in) :: qphon(3)
 real(dp), intent(in) :: rprimd(3,3)
 real(dp), intent(in) :: xcart(3,natom)

!Local variables-------------------------------
!scalar
 integer :: ii, maxsc, i1,i2,i3, iscmult
 integer :: iatom, iatom_supercell
 real(dp) :: qbymult
!arrays
 integer :: r_cell(3)
 integer :: supercell(3) ! number of primitive cells in each direction for the supercell
 character(len=500) :: msg
! *************************************************************************

! maximum number of unit cells in a given direction
 maxsc = 10  

! find smallest supercell which will accomodate phonon.
! FIXME: for the moment, just get smallest multiple along each direction, with an upper bound
 supercell = -1
 do ii=1,3
   do iscmult=1,maxsc
     qbymult = qphon(ii)*iscmult
     if (abs(qbymult - int(qbymult)) < tol10) then
       supercell(ii) = iscmult
       exit
     end if
   end do
   if (supercell(ii) == -1) then
     write(msg,'(a,I4,a,I7,2a,3E20.10)')' No supercell found with less than ', &
&             maxsc,' unit cells in direction ', &
&             ii, ch10, ' qphon = ', qphon
     MSG_ERROR(msg)
   end if
 end do
 
 scell%natom = natom
 scell%qphon = qphon
 scell%rprimd_supercell(:,1) = rprimd(:,1) * supercell(1) 
 scell%rprimd_supercell(:,2) = rprimd(:,2) * supercell(2) 
 scell%rprimd_supercell(:,3) = rprimd(:,3) * supercell(3) 

!number of atoms in full supercell
 scell%natom_supercell = natom*supercell(1)*supercell(2)*supercell(3)
 ABI_ALLOCATE(scell%xcart_supercell,(3,scell%natom_supercell))
 ABI_ALLOCATE(scell%xcart_supercell_ref,(3,scell%natom_supercell))
 ABI_ALLOCATE(scell%atom_indexing_supercell,(scell%natom_supercell))
 ABI_ALLOCATE(scell%uc_indexing_supercell,(3,scell%natom_supercell))

 iatom_supercell = 0
 do i1 = 1, supercell(1) 
   do i2 = 1, supercell(2) 
     do i3 = 1, supercell(3) 
       do iatom = 1, natom
         iatom_supercell = iatom_supercell + 1
         r_cell = (/i1-1,i2-1,i3-1/)
         scell%xcart_supercell_ref(:,iatom_supercell) = xcart(:,iatom) + matmul(rprimd,r_cell)
         scell%atom_indexing_supercell(iatom_supercell) = iatom
         scell%uc_indexing_supercell(:,iatom_supercell) = r_cell
       end do
     end do
   end do
 end do

 ABI_CHECK(iatom_supercell == scell%natom_supercell, "iatom_supercell /= scell%natom_supercell")

 scell%xcart_supercell = scell%xcart_supercell_ref

end subroutine init_supercell 
!!***

!!****f* m_phonon_supercell/freeze_displ_supercell
!!
!! NAME
!! freeze_displ_supercell
!!
!! FUNCTION
!! Freeze a specific displacement phonon field into the supercell scell
!!
!! INPUTS
!! displ = phonon displacement vectors for all modes
!! freeze_displ = desired amplitude for phonon displacement along displ
!! jmode = phonon mode number to be used
!! scell = supercell structure with reference atomic positions etc...
!!
!! OUTPUT
!! scell = supercell structure: xcart will be updated with phonon displacement
!!
!! PARENTS
!!      freeze_displ_allmodes
!!
!! CHILDREN
!!
!! SOURCE

subroutine freeze_displ_supercell (displ,freeze_displ,jmode,scell)

  use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'freeze_displ_supercell'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(supercell_type), intent(inout) :: scell
 integer, intent(in) :: jmode
 real(dp), intent(in) :: freeze_displ
!arrays
 real(dp), intent(in) :: displ(2,3*scell%natom,3*scell%natom)

!Local variables-------------------------------
!scalar
 integer :: iatom, ipratom
 real(dp) :: qdotr
! *************************************************************************

 do iatom = 1, scell%natom_supercell
   qdotr = two_pi*(scell%qphon(1)*scell%uc_indexing_supercell(1,iatom) &
&                 +scell%qphon(2)*scell%uc_indexing_supercell(2,iatom) &
&                 +scell%qphon(3)*scell%uc_indexing_supercell(3,iatom))

!add real part of displacement times Bloch phase
   ipratom = (scell%atom_indexing_supercell(iatom)-1)*3
   scell%xcart_supercell(:,iatom) = scell%xcart_supercell_ref(:,iatom) &
&        + freeze_displ * cos(qdotr) * displ(1,ipratom+1:ipratom+3,jmode) &
&        - freeze_displ * sin(qdotr) * displ(2,ipratom+1:ipratom+3,jmode) 
 end do

end subroutine freeze_displ_supercell
!!***

!****f* m_phonon_supercell/prt_supercell
!!
!! NAME
!! prt_supercell
!!
!! FUNCTION
!! output atomic positions, supercell vectors, etc... to a file
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
!!      freeze_displ_allmodes
!!
!! CHILDREN
!!
!! SOURCE

subroutine prt_supercell (freq, jmode, outfile_radix, scell, typat, znucl)

  use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prt_supercell'
 use interfaces_32_util
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  real(dp), intent(in) :: freq
  type(supercell_type), intent(in) :: scell
  integer, intent(in) :: jmode
  character(len=fnlen), intent(in) :: outfile_radix
!arrays
  integer, intent(in) :: typat(scell%natom)
  real(dp), intent(in) :: znucl(:)

!Local variables-------------------------------
!scalar
  character(len=fnlen) :: filename
  integer :: scunit, iatom
  character(len=10) :: jmodestring
  character(len=5) :: qphonstring1, qphonstring2, qphonstring3
  character(len=500) :: msg
  real(dp) :: xred(3), gprimd(3,3)
! *************************************************************************

! add suffix with mode and qpoint
  call int2char4(jmode, jmodestring)
  ABI_CHECK((jmodestring(1:1)/='#'),'Bug: string length too short!')
! qphonstring should be like 0.000_0.000_0.000
  call write_num(scell%qphon(1),qphonstring1,'(F5.3)')
  call write_num(scell%qphon(2),qphonstring2,'(F5.3)')
  call write_num(scell%qphon(3),qphonstring3,'(F5.3)')
  filename = trim(outfile_radix) // "_qpt_" // qphonstring1 // "_" // qphonstring2 // &
&              "_" // qphonstring3 // "_mode_" // trim(jmodestring)
  if (open_file(filename, msg, newunit=scunit) /= 0) then
    MSG_ERROR(msg)
  end if

! print header
  write (scunit, '(a)') '#'
  write (scunit, '(a)') '# anaddb file with frozen phonon mode in supercell'
  write (scunit, '(a)') '#'
  write (scunit, '(a,3E20.10)') '# phonon q point : ', scell%qphon
  write (scunit, '(a,I7,a,E20.10)') '# phonon mode number : ', jmode, ' frequency ', freq
  write (scunit, '(a)') '#'
  write (scunit, '(a)') '# lattice vectors for supercell :'
  write (scunit, '(a,I7)') 'natom ', scell%natom_supercell
  write (scunit, *)
  write (scunit, '(a)') 'znucl '
  do iatom = 1, size(znucl)
    write (scunit, '(I5)', ADVANCE="NO") int(znucl(iatom))
    if (mod(iatom,6) == 0) write (scunit, *)
  end do
  write (scunit, *)
  write (scunit, *)
  write (scunit, '(a,I7)') 'ntypat', size(znucl)
  write (scunit, '(a)') 'typat '
  do iatom = 1, scell%natom_supercell
    write (scunit, '(I5)', ADVANCE="NO") typat(scell%atom_indexing_supercell(iatom))
    if (mod(iatom,6) == 0) write (scunit, *)
  end do
  write (scunit, *)
  write (scunit, '(a)') 'acell 1.0 1.0 1.0'
  write (scunit, '(a)') 'rprim'
  write (scunit, '(3E20.10)') scell%rprimd_supercell(:,1)
  write (scunit, '(3E20.10)') scell%rprimd_supercell(:,2)
  write (scunit, '(3E20.10)') scell%rprimd_supercell(:,3)
  write (scunit, *)
  write (scunit, '(a)') 'xcart'
  do iatom = 1, scell%natom_supercell
    write (scunit, '(3E20.10)') scell%xcart_supercell(:,iatom)
  end do
  ! for information, also print xred for atoms inside full supercell
  call matr3inv(scell%rprimd_supercell, gprimd)
  ! TODO: check this transpose is correct in some asymetric case
  gprimd = transpose(gprimd)
  write (scunit, '(a)') '# for information, add xred as well'
  write (scunit, '(a)') '# xred'
  do iatom = 1, scell%natom_supercell
    xred = matmul (gprimd, scell%xcart_supercell(:,iatom))
    write (scunit, '(a, 3E20.10)') '#  ', xred
  end do

! close file
  close(scunit)

end subroutine prt_supercell
!!***

!****f* m_phonon_supercell/destroy_supercell
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
!!      freeze_displ_allmodes
!!
!! CHILDREN
!!
!! SOURCE
 
subroutine destroy_supercell (scell)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_supercell'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  type(supercell_type), intent(inout) :: scell !vz_i

! *************************************************************************

  if(allocated(scell%xcart_supercell))  then
    ABI_DEALLOCATE(scell%xcart_supercell)
  end if
  if(allocated(scell%xcart_supercell_ref))  then
    ABI_DEALLOCATE(scell%xcart_supercell_ref)
  end if
  if(allocated(scell%atom_indexing_supercell))  then
    ABI_DEALLOCATE(scell%atom_indexing_supercell)
  end if
  if(allocated(scell%uc_indexing_supercell))  then
    ABI_DEALLOCATE(scell%uc_indexing_supercell)
  end if

end subroutine destroy_supercell
!!***

end module m_phonon_supercell
!!***
