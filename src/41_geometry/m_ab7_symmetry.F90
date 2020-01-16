!!****m* ABINIT/m_ab7_symmetry
!! NAME
!! m_ab7_symmetry
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group (DC)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_ab7_symmetry

  use defs_basis
  use m_abicore

  use m_symtk,     only : mati3inv, mati3det, symatm, symcharac
  use m_symfind,   only : symfind, symanal, symlatt
  use m_geometry,  only : metric
  use m_spgdata,   only : spgdata

  implicit none

  private

  integer, parameter, public :: AB7_MAX_SYMMETRIES = 384

  type, public :: symmetry_type
     ! The input characteristics
     real(dp) :: tolsym
     real(dp) :: rprimd(3,3), gprimd(3,3), rmet(3,3)
     integer :: nAtoms
     integer, pointer :: typeAt(:)
     real(dp), pointer :: xRed(:,:)

     logical :: withField
     real(dp) :: field(3)

     logical :: withJellium

     integer :: nzchempot

     integer :: withSpin
     real(dp), pointer :: spinAt(:,:)

     logical :: withSpinOrbit

     integer :: vacuum(3)

     ! The output characteristics
     ! The bravais parameters
     integer :: nBravSym
     integer :: bravais(11), bravSym(3, 3, AB7_MAX_SYMMETRIES)
     ! The symmetry matrices
     logical  :: auto
     integer  :: nSym
     integer, pointer  :: sym(:,:,:)
     real(dp), pointer :: transNon(:,:)
     integer, pointer  :: symAfm(:)
     ! Some additional information
     integer          :: multiplicity
     real(dp)         :: genAfm(3)
     integer          :: spaceGroup, pointGroupMagn
     integer, pointer :: indexingAtoms(:,:,:)
  end type symmetry_type

  ! We store here a list of symmetry objects to be able to
  ! call several symmetry operations on different objects.
  ! The simplest portable way to do it, is to create
  ! a list of Fortran structure and to use the list index
  ! as an identifier that can be given to the other languages.
  type, private :: symmetry_list
     integer                       :: id
     type(symmetry_list),  pointer :: next
     type(symmetry_type)           :: data
  end type symmetry_list
  type(symmetry_list), pointer :: my_symmetries
  integer :: n_symmetries = 0

  logical, private, parameter :: AB_DBG = .false.

  public :: symmetry_new
  public :: symmetry_free
  public :: symmetry_set_tolerance
  public :: symmetry_set_lattice
  public :: symmetry_set_structure
  public :: symmetry_set_collinear_spin
  public :: symmetry_set_spin
  public :: symmetry_set_spin_orbit
  public :: symmetry_set_field
  public :: symmetry_set_jellium
  public :: symmetry_set_periodicity
  public :: symmetry_set_n_sym

  public :: symmetry_get_from_id
  public :: symmetry_get_n_atoms
  public :: symmetry_get_n_sym
  public :: symmetry_get_multiplicity
  public :: symmetry_get_bravais
  public :: symmetry_get_matrices
  public :: symmetry_get_matrices_p
  public :: symmetry_get_group
  public :: symmetry_get_equivalent_atom
  public :: symmetry_get_type

contains

  subroutine new_item(token)

    type(symmetry_list), pointer :: token

    ! We allocate a new list token and prepend it.
    if (AB_DBG) write(std_err,*) "AB symmetry: create a new token."

    ! Init case, very first call.
    if (n_symmetries == 0) then
       nullify(my_symmetries)
    end if

    ! Normal treatment.
    n_symmetries = n_symmetries + 1

    ABI_DATATYPE_ALLOCATE(token,)
    token%id = n_symmetries
    call new_symmetry(token%data)
    token%next => my_symmetries

    my_symmetries => token
    if (AB_DBG) write(std_err,*) "AB symmetry: creation OK with id ", token%id
  end subroutine new_item

  subroutine free_item(token)

    type(symmetry_list), pointer :: token

    type(symmetry_list), pointer :: tmp

    if (.not. associated(token)) then
       return
    end if

    call free_symmetry(token%data)

    if (AB_DBG) write(std_err,*) "AB symmetry: free request on token ", token%id
    ! We remove token from the list.
    if (my_symmetries%id == token%id) then
       my_symmetries => token%next
    else
       tmp => my_symmetries
       do
          if (.not.associated(tmp)) then
             return
          end if
          if (associated(tmp%next) .and. tmp%next%id == token%id) then
             exit
          end if
          tmp => tmp%next
       end do
       tmp%next => token%next
    end if
    ABI_DATATYPE_DEALLOCATE(token)
    if (AB_DBG) write(std_err,*) "AB symmetry: free done"
  end subroutine free_item

  subroutine get_item(token, id)

    type(symmetry_list), pointer :: token
    integer, intent(in) :: id

    type(symmetry_list), pointer :: tmp

    if (AB_DBG) write(std_err,*) "AB symmetry: request list element ", id
    nullify(token)

    tmp => my_symmetries
    do
       if (.not. associated(tmp)) then
          exit
       end if
       if (tmp%id == id) then
          token => tmp
          return
       end if
       tmp => tmp%next
    end do
  end subroutine get_item

  subroutine symmetry_get_from_id(sym, id, errno)

    type(symmetry_type), pointer :: sym
    integer, intent(in) :: id
    integer, intent(out) :: errno

    type(symmetry_list), pointer :: token

    errno = AB7_NO_ERROR
    call get_item(token, id)
    if (associated(token)) then
       sym => token%data
       if (sym%nSym <= 0) then
          ! We do the computation of the matrix part.
          call compute_matrices(sym, errno)
       end if
    else
       errno = AB7_ERROR_OBJ
       nullify(sym)
    end if
  end subroutine symmetry_get_from_id

  subroutine new_symmetry(sym)

    type(symmetry_type), intent(out) :: sym

    if (AB_DBG) write(std_err,*) "AB symmetry: create a new symmetry object."
    nullify(sym%xRed)
    nullify(sym%spinAt)
    nullify(sym%typeAt)
    sym%tolsym   = tol8
    sym%auto     = .true.
    sym%nSym     = 0
    nullify(sym%sym)
    nullify(sym%symAfm)
    nullify(sym%transNon)
    sym%nBravSym = -1
    sym%withField   = .false.
    sym%withJellium = .false.
    sym%nzchempot = 0
    sym%withSpin = 1
    sym%withSpinOrbit = .false.
    sym%multiplicity = -1
    nullify(sym%indexingAtoms)
    sym%vacuum = 0
  end subroutine new_symmetry

  subroutine free_symmetry(sym)

    type(symmetry_type), intent(inout) :: sym

    if (AB_DBG) write(std_err,*) "AB symmetry: free a symmetry."

    if (associated(sym%xRed))  then
      ABI_DEALLOCATE(sym%xRed)
    end if
    if (associated(sym%spinAt))  then
      ABI_DEALLOCATE(sym%spinAt)
    end if
    if (associated(sym%typeAt))  then
      ABI_DEALLOCATE(sym%typeAt)
    end if
    if (associated(sym%indexingAtoms))  then
      ABI_DEALLOCATE(sym%indexingAtoms)
    end if
    if (associated(sym%sym))  then
      ABI_DEALLOCATE(sym%sym)
    end if
    if (associated(sym%symAfm))  then
      ABI_DEALLOCATE(sym%symAfm)
    end if
    if (associated(sym%transNon))  then
      ABI_DEALLOCATE(sym%transNon)
    end if
  end subroutine free_symmetry





  subroutine symmetry_new(id)

    integer, intent(out) :: id

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(std_err,*) "AB symmetry: call new symmetry."
    call new_item(token)
    id = token%id
  end subroutine symmetry_new

  subroutine symmetry_free(id)

    integer, intent(in) :: id

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(std_err,*) "AB symmetry: call free symmetry."

    call get_item(token, id)
    if (associated(token)) then
      call free_item(token)
    end if
  end subroutine symmetry_free

  subroutine symmetry_set_tolerance(id, tolsym, errno)

    integer, intent(in) :: id
    real(dp), intent(in) :: tolsym
    integer, intent(out) :: errno

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(std_err,*) "AB symmetry: call set tolerance."

    errno = AB7_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB7_ERROR_OBJ
       return
    end if

    token%data%tolsym = tolsym

    ! We unset all the computed symmetries
    token%data%nBravSym = -1
    if (token%data%auto) then
       token%data%nSym  = 0
    end if
  end subroutine symmetry_set_tolerance

  subroutine symmetry_set_lattice(id, rprimd, errno)

    integer, intent(in) :: id
    real(dp), intent(in) :: rprimd(3,3)
    integer, intent(out) :: errno

    type(symmetry_list), pointer :: token
    real(dp) :: ucvol
    real(dp) :: gmet(3,3)

    if (AB_DBG) write(std_err,*) "AB symmetry: call set lattice."
    if (AB_DBG) write(std_err, "(A,3F12.6,A)") "  (", rprimd(:,1), ")"
    if (AB_DBG) write(std_err, "(A,3F12.6,A)") "  (", rprimd(:,2), ")"
    if (AB_DBG) write(std_err, "(A,3F12.6,A)") "  (", rprimd(:,3), ")"

    errno = AB7_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB7_ERROR_OBJ
       return
    end if

    token%data%rprimd = rprimd
    call metric(gmet, token%data%gprimd, -1, token%data%rmet, rprimd, ucvol)

    ! We unset all the computed symmetries
    token%data%nBravSym = -1
    if (token%data%auto) then
       token%data%nSym  = 0
    end if
  end subroutine symmetry_set_lattice

  subroutine symmetry_set_structure(id, nAtoms, typeAt, xRed, errno)

    integer, intent(in) :: id
    integer, intent(in) :: nAtoms
    integer, intent(in) :: typeAt(nAtoms)
    real(dp), intent(in) :: xRed(3,nAtoms)
    integer, intent(out) :: errno

    type(symmetry_list), pointer :: token
    integer :: i

    if (AB_DBG) write(std_err,*) "AB symmetry: call set structure."
    if (AB_DBG) write(std_err, "(A,I3,A)") "  ", nAtoms, " atoms"
    if (AB_DBG) then
       do i = 1, nAtoms, 1
          if (AB_DBG) write(std_err,"(A,3F12.6,I3)") "  ", xRed(:, i), typeAt(i)
       end do
    end if

    errno = AB7_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB7_ERROR_OBJ
       return
    end if

    token%data%nAtoms =  nAtoms
    ABI_ALLOCATE(token%data%typeAt,(nAtoms))
    token%data%typeAt = typeAt
    ABI_ALLOCATE(token%data%xRed,(3, nAtoms))
    token%data%xRed   = xRed

    ! We unset only the symmetries
    if (token%data%auto) then
       token%data%nSym = 0
    end if
    if (associated(token%data%indexingAtoms)) then
      ABI_DEALLOCATE(token%data%indexingAtoms)
    end if
  end subroutine symmetry_set_structure

  subroutine symmetry_set_spin(id, nAtoms, spinAt, errno)

    integer, intent(in) :: id
    integer, intent(in) :: nAtoms
    real(dp), intent(in) :: spinAt(3,nAtoms)
    integer, intent(out) :: errno

    type(symmetry_list), pointer :: token
    integer :: i

    if (AB_DBG) write(std_err,*) "AB symmetry: call set spin."
    if (AB_DBG) then
       do i = 1, nAtoms, 1
          if (AB_DBG) write(std_err,"(A,3F12.6)") "  ", spinAt(:, i)
       end do
    end if

    errno = AB7_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB7_ERROR_OBJ
       return
    end if
    if (token%data%nAtoms /= nAtoms) then
       errno = AB7_ERROR_ARG
       return
    end if

    token%data%withSpin = 4
    ABI_ALLOCATE(token%data%spinAt,(3, nAtoms))
    token%data%spinAt = spinAt

    ! We unset only the symmetries
    if (token%data%auto) then
       token%data%nSym  = 0
    end if
  end subroutine symmetry_set_spin

  subroutine symmetry_set_collinear_spin(id, nAtoms, spinAt, errno)

    integer, intent(in) :: id
    integer, intent(in) :: nAtoms
    integer, intent(in) :: spinAt(nAtoms)
    integer, intent(out) :: errno

    type(symmetry_list), pointer :: token
    integer :: i

    if (AB_DBG) write(std_err,*) "AB symmetry: call set collinear spin."
    if (AB_DBG) then
       do i = 1, nAtoms, 1
          if (AB_DBG) write(std_err,"(A,I3)") "  ", spinAt(i)
       end do
    end if

    errno = AB7_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB7_ERROR_OBJ
       return
    end if
    if (token%data%nAtoms /= nAtoms) then
       errno = AB7_ERROR_ARG
       return
    end if

    token%data%withSpin = 2
    ABI_ALLOCATE(token%data%spinAt,(3,nAtoms))
    token%data%spinAt = 0._dp
    token%data%spinAt(3, :) = real(spinAt)

    ! We unset only the symmetries
    if (token%data%auto) then
       token%data%nSym  = 0
    end if
  end subroutine symmetry_set_collinear_spin

  subroutine symmetry_set_spin_orbit(id, withSpinOrbit, errno)

    integer, intent(in) :: id
    logical, intent(in) :: withSpinOrbit
    integer, intent(out) :: errno

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(std_err,*) "AB symmetry: call set spin orbit."

    errno = AB7_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB7_ERROR_OBJ
       return
    end if

    token%data%withSpinOrbit = withSpinOrbit

    ! We unset only the symmetries
    if (token%data%auto) then
       token%data%nSym  = 0
    end if
  end subroutine symmetry_set_spin_orbit

  subroutine symmetry_set_field(id, field, errno)

    integer, intent(in) :: id
    real(dp), intent(in) :: field(3)
    integer, intent(out) :: errno

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(std_err,*) "AB symmetry: call set field."

    errno = AB7_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB7_ERROR_OBJ
       return
    end if

    token%data%withField = .true.
    token%data%field = field

    ! We unset all the computed symmetries
    token%data%nBravSym = -1
    if (token%data%auto) then
       token%data%nSym  = 0
    end if
  end subroutine symmetry_set_field

  subroutine symmetry_set_jellium(id, jellium, errno)

    integer, intent(in) :: id
    logical, intent(in) :: jellium
    integer, intent(out) :: errno

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(std_err,*) "AB symmetry: call set jellium."

    errno = AB7_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB7_ERROR_OBJ
       return
    end if

    token%data%withJellium = jellium

    ! We unset only the symmetries
    if (token%data%auto) then
       token%data%nSym  = 0
    end if
  end subroutine symmetry_set_jellium


  subroutine symmetry_set_nzchempot(id, nzchempot, errno)

    integer, intent(in) :: id
    integer, intent(in) :: nzchempot
    integer, intent(out) :: errno

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(std_err,*) "AB symmetry: call set nzchempot."

    errno = AB7_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB7_ERROR_OBJ
       return
    end if

    token%data%nzchempot = nzchempot

    ! We unset only the symmetries
    if (token%data%auto) then
       token%data%nSym  = 0
    end if
  end subroutine symmetry_set_nzchempot

  subroutine symmetry_set_periodicity(id, periodic, errno)

    integer, intent(in) :: id
    logical, intent(in) :: periodic(3)
    integer, intent(out) :: errno

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(std_err,*) "AB symmetry: call set periodicity."
    if (AB_DBG) write(std_err, "(A,3L1,A)") "  (", periodic, ")"

    errno = AB7_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB7_ERROR_OBJ
       return
    end if

    token%data%vacuum = 0
    if (.not. periodic(1)) token%data%vacuum(1) = 1
    if (.not. periodic(2)) token%data%vacuum(2) = 1
    if (.not. periodic(3)) token%data%vacuum(3) = 1
  end subroutine symmetry_set_periodicity





  subroutine symmetry_get_n_atoms(id, nAtoms, errno)
    !scalars
    integer, intent(in) :: id
    integer, intent(out) :: errno
    integer, intent(out) :: nAtoms

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(std_err,*) "AB symmetry: call get nAtoms."

    errno = AB7_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB7_ERROR_OBJ
       return
    end if

    nAtoms = token%data%nAtoms
  end subroutine symmetry_get_n_atoms

  subroutine compute_bravais(sym)

    type(symmetry_type), intent(inout) :: sym

    integer :: berryopt

    ! We do the computation
    if (sym%withField) then
       berryopt = 4
    else
       berryopt = 0
    end if
    if (AB_DBG) write(std_err,*) "AB symmetry: call ABINIT symlatt."
    call symlatt(sym%bravais, AB7_MAX_SYMMETRIES, &
         & sym%nBravSym, sym%bravSym, sym%rprimd, sym%tolsym)
    if (AB_DBG) write(std_err,*) "AB symmetry: call ABINIT OK."
    if (AB_DBG) write(std_err, "(A,I3)") "  nSymBrav :", sym%nBravSym
    if (AB_DBG) write(std_err, "(A,I3)") "  holohedry:", sym%bravais(1)
    if (AB_DBG) write(std_err, "(A,I3)") "  center   :", sym%bravais(2)
  end subroutine compute_bravais

  subroutine symmetry_get_bravais(id, bravais, holohedry, center, &
       & nBravSym, bravSym, errno)
    !scalars
    integer, intent(in) :: id
    integer, intent(out) :: errno
    integer, intent(out) :: nBravSym, holohedry, center
    !arrays
    integer, intent(out) :: bravais(3,3), bravSym(3, 3, AB7_MAX_SYMMETRIES)

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(std_err,*) "AB symmetry: call get bravais."

    errno = AB7_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB7_ERROR_OBJ
       return
    end if

    if (token%data%nBravSym < 0) then
       ! We do the computation
       call compute_bravais(token%data)
    end if

    holohedry = token%data%bravais(1)
    center    = token%data%bravais(2)
    bravais   = reshape(token%data%bravais(3:11), (/ 3,3 /))
    nBravSym  = token%data%nBravSym
    bravSym(:, :, 1:nBravSym) = token%data%bravSym(:, :, 1:nBravSym)
  end subroutine symmetry_get_bravais

  subroutine compute_matrices(sym, errno)

    type(symmetry_type), intent(inout) :: sym
    integer, intent(out) :: errno

    integer :: berryopt, jellslab, noncol
    integer :: use_inversion
    real(dp), pointer :: spinAt_(:,:)
    integer  :: sym_(3, 3, AB7_MAX_SYMMETRIES)
    real(dp) :: transNon_(3, AB7_MAX_SYMMETRIES)
    integer  :: symAfm_(AB7_MAX_SYMMETRIES)

    errno = AB7_NO_ERROR

    if (sym%nBravSym < 0) then
       ! We do the computation of the Bravais part.
       call compute_bravais(sym)
    end if

    if (sym%withField) then
       berryopt = 4
    else
       berryopt = 0
    end if
    if (sym%withJellium) then
       jellslab = 1
    else
       jellslab = 0
    end if
    if (sym%withSpin == 4) then
       noncol = 1
       spinAt_ => sym%spinAt
    else if (sym%withSpin == 2) then
       noncol = 0
       spinAt_ => sym%spinAt
    else
       noncol = 0
       ABI_ALLOCATE(spinAt_,(3, sym%nAtoms))
       spinAt_ = 0
    end if
    if (sym%withSpinOrbit) then
       use_inversion = 0
    else
       use_inversion = 1
    end if

    if (sym%nsym == 0) then
       if (AB_DBG) write(std_err,*) "AB symmetry: call ABINIT symfind."
       call symfind(berryopt, sym%field, sym%gprimd, jellslab, AB7_MAX_SYMMETRIES, &
            & sym%nAtoms, noncol, sym%nBravSym, sym%nSym, sym%nzchempot, 0, sym%bravSym, spinAt_, &
            & symAfm_, sym_, transNon_, sym%tolsym, sym%typeAt, &
            & use_inversion, sym%xRed)
       if (AB_DBG) write(std_err,*) "AB symmetry: call ABINIT OK."
       if (AB_DBG) write(std_err, "(A,I3)") "  nSym:", sym%nSym
       if (associated(sym%sym))  then
          ABI_DEALLOCATE(sym%sym)
       end if
       if (associated(sym%symAfm))  then
          ABI_DEALLOCATE(sym%symAfm)
       end if
       if (associated(sym%transNon))  then
          ABI_DEALLOCATE(sym%transNon)
       end if
       ABI_ALLOCATE(sym%sym,(3, 3, sym%nSym))
       sym%sym(:,:,:) = sym_(:,:, 1:sym%nSym)
       ABI_ALLOCATE(sym%symAfm,(sym%nSym))
       sym%symAfm(:) = symAfm_(1:sym%nSym)
       ABI_ALLOCATE(sym%transNon,(3, sym%nSym))
       sym%transNon(:,:) = transNon_(:, 1:sym%nSym)
    else if (sym%nsym < 0) then
       sym%nsym = -sym%nsym
       sym_(:,:, 1:sym%nSym) = sym%sym(:,:,:)
       transNon_(:, 1:sym%nSym) = sym%transNon(:,:)
       symAfm_(1:sym%nSym) = sym%symAfm(:)
    end if

    if (sym%withSpin == 1) then
       ABI_DEALLOCATE(spinAt_)
    end if

    if (AB_DBG) write(std_err,*) "AB symmetry: call ABINIT symanal."
    call symanal(sym%bravais, 0, sym%genAfm, AB7_MAX_SYMMETRIES, sym%nSym, &
         & sym%pointGroupMagn, sym%rprimd, sym%spaceGroup, symAfm_, &
         & sym_, transNon_, sym%tolsym)
    if (AB_DBG) write(std_err,*) "AB symmetry: call ABINIT OK."
    sym%transNon(:,:) = transNon_(:, 1:sym%nSym)

    if (sym%bravais(1) < 0) then
       sym%multiplicity = 2
    else
       sym%multiplicity = 1
    end if
    if (AB_DBG) write(std_err, "(A,I3)") "  multi:", sym%multiplicity
    if (AB_DBG) write(std_err, "(A,I3)") "  space:", sym%spaceGroup
  end subroutine compute_matrices

  subroutine symmetry_get_n_sym(id, nSym, errno)
    !scalars
    integer, intent(in) :: id
    integer, intent(out) :: errno
    integer, intent(out) :: nSym

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(std_err,*) "AB symmetry: call get nSym."

    errno = AB7_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB7_ERROR_OBJ
       return
    end if

    if (token%data%nSym <= 0) then
       ! We do the computation of the matrix part.
       call compute_matrices(token%data, errno)
    end if

    nSym = token%data%nSym
  end subroutine symmetry_get_n_sym

  subroutine symmetry_set_n_sym(id, nSym, sym, transNon, symAfm, errno)
    !scalars
    integer, intent(in)  :: id
    integer, intent(in)  :: nSym
    integer, intent(in)  :: sym(3, 3, nSym)
    real(dp), intent(in) :: transNon(3, nSym)
    integer, intent(in)  :: symAfm(nSym)
    integer, intent(out) :: errno

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(std_err,*) "AB symmetry: call get nSym."

    errno = AB7_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB7_ERROR_OBJ
       return
    end if

    if (nSym <= 0) then
       errno = AB7_ERROR_ARG
       return
    else
       ABI_ALLOCATE(token%data%sym, (3, 3, nSym))
       token%data%sym(:,:,:) = sym(:,:,:)
       ABI_ALLOCATE(token%data%symAfm, (nSym))
       token%data%symAfm(:) = symAfm(:)
       ABI_ALLOCATE(token%data%transNon, (3, nSym))
       token%data%transNon(:,:) = transNon(:,:)

       token%data%auto = .false.
       token%data%nsym = -nSym
    end if

    ! We do the computation of the matrix part.
    call compute_matrices(token%data, errno)
  end subroutine symmetry_set_n_sym

  subroutine symmetry_get_matrices(id, nSym, sym, transNon, symAfm, errno)

    integer, intent(in) :: id
    integer, intent(out) :: errno
    integer, intent(out) :: nSym
    integer, intent(out)  :: sym(3, 3, AB7_MAX_SYMMETRIES)
    integer, intent(out)  :: symAfm(AB7_MAX_SYMMETRIES)
    real(dp), intent(out) :: transNon(3, AB7_MAX_SYMMETRIES)

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(std_err,*) "AB symmetry: call get matrices."

    errno = AB7_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB7_ERROR_OBJ
       return
    end if

    if (token%data%nSym <= 0) then
       ! We do the computation of the matrix part.
       call compute_matrices(token%data, errno)
    end if

    nSym                = token%data%nSym
    sym(:, :, 1:nSym)   = token%data%sym(:, :,:)
    symAfm(1:nSym)      = token%data%symAfm(:)
    transNon(:, 1:nSym) = token%data%transNon(:,:)
  end subroutine symmetry_get_matrices

  subroutine symmetry_get_matrices_p(id, nSym, sym, transNon, symAfm, errno)

    integer, intent(in) :: id
    integer, intent(out) :: errno
    integer, intent(out) :: nSym
    integer, pointer  :: sym(:,:,:)
    integer, pointer  :: symAfm(:)
    real(dp), pointer :: transNon(:,:)

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(std_err,*) "AB symmetry: call get matrices as pointers."

    errno = AB7_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB7_ERROR_OBJ
       return
    end if

    if (token%data%nSym <= 0) then
       ! We do the computation of the matrix part.
       call compute_matrices(token%data, errno)
    end if

    nSym     =  token%data%nSym
    sym      => token%data%sym
    symAfm   => token%data%symAfm
    transNon => token%data%transNon
  end subroutine symmetry_get_matrices_p

  subroutine symmetry_get_multiplicity(id, multiplicity, errno)

    integer, intent(in) :: id
    integer, intent(out) :: multiplicity, errno

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(std_err,*) "AB symmetry: call get multiplicity."

    errno = AB7_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB7_ERROR_OBJ
       return
    end if

    if (token%data%multiplicity < 0) then
       ! We do the computation of the matrix part.
       call compute_matrices(token%data, errno)
    end if
    multiplicity = token%data%multiplicity
  end subroutine symmetry_get_multiplicity

  subroutine symmetry_get_group(id, spaceGroup, spaceGroupId, &
       & pointGroupMagn, genAfm, errno)

    integer, intent(in)            :: id
    integer, intent(out)           :: errno
    real(dp), intent(out)          :: genAfm(3)
    character(len=15), intent(out) :: spaceGroup
    integer, intent(out)           :: spaceGroupId, pointGroupMagn

    type(symmetry_list), pointer  :: token
    integer :: sporder
    character(len=1)  :: brvLattice
    character(len=15) :: ptintsb,ptschsb,schsb,spgrp
    character(len=35) :: intsbl

    if (AB_DBG) write(std_err,*) "AB symmetry: call get group."

    errno = AB7_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB7_ERROR_OBJ
       return
    end if

    if (token%data%multiplicity < 0) then
       ! We do the computation of the matrix part.
       call compute_matrices(token%data, errno)
    end if

    if (token%data%multiplicity /= 1) then
       errno = AB7_ERROR_SYM_NOT_PRIMITIVE
       return
    end if

    call spgdata(brvLattice,spgrp,intsbl,ptintsb,ptschsb,&
         &  schsb,1,token%data%spaceGroup,sporder,1)

    write(spaceGroup, "(3A)") brvLattice, " ", trim(spgrp(1:13))
    pointGroupMagn = token%data%pointGroupMagn
    spaceGroupId   = token%data%spaceGroup
    genAfm         = token%data%genAfm
  end subroutine symmetry_get_group

  subroutine compute_equivalent_atoms(sym)

    type(symmetry_type), intent(inout) :: sym

    integer, allocatable :: symrec(:,:,:)
    integer :: isym

    if (.not. associated(sym%indexingAtoms)) then
      ABI_ALLOCATE(sym%indexingAtoms,(4, sym%nSym, sym%nAtoms))
    end if

    !Get the symmetry matrices in terms of reciprocal basis
    ABI_ALLOCATE(symrec,(3, 3, sym%nSym))
    do isym = 1, sym%nSym, 1
       call mati3inv(sym%sym(:,:,isym), symrec(:,:,isym))
    end do

    !Obtain a list of rotated atom labels:
    call symatm(sym%indexingAtoms, sym%nAtoms, sym%nSym, symrec, &
         & sym%transNon, sym%tolsym, sym%typeAt, sym%xRed)

    ABI_DEALLOCATE(symrec)
  end subroutine compute_equivalent_atoms

  subroutine symmetry_get_equivalent_atom(id, equiv, iAtom, errno)

    integer, intent(in)  :: id
    integer, intent(in)  :: iAtom
    integer, intent(out) :: equiv(4, AB7_MAX_SYMMETRIES)
    integer, intent(out) :: errno

    type(symmetry_list), pointer  :: token

    if (AB_DBG) write(std_err,*) "AB symmetry: call get equivalent."

    errno = AB7_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB7_ERROR_OBJ
       return
    end if

    if (iAtom < 1 .or. iAtom > token%data%nAtoms) then
       errno = AB7_ERROR_ARG
       return
    end if

    if (.not. associated(token%data%indexingAtoms)) then
       ! We do the computation of the matrix part.
       call compute_equivalent_atoms(token%data)
    end if

    equiv(:, 1:token%data%nSym) = token%data%indexingAtoms(:,:,iAtom)
  end subroutine symmetry_get_equivalent_atom

  subroutine symmetry_get_type(id, iSym, label, type, errno)

    integer, intent(in)  :: id, iSym
    character(len = 128), intent(out) :: label
    integer, intent(out) :: errno, type

    integer :: det
    type(symmetry_list), pointer :: token
    type(symmetry_type), pointer :: sym

    if (AB_DBG) write(std_err,*) "AB symmetry: call get type."

    errno = AB7_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB7_ERROR_OBJ
       return
    end if
    sym => token%data

    if (iSym < 1 .or. iSym > sym%nSym) then
       errno = AB7_ERROR_ARG
       return
    end if

    if (sym%multiplicity < 0) then
       ! We do the computation of the matrix part.
       call compute_matrices(sym, errno)
    end if

    ! Calculate determinant.
    call mati3det(sym%sym(:,:,iSym),det)
    call symcharac(sym%bravais(2), det, sym%bravais(1), iSym, label, &
         sym%sym(:,:,iSym), sym%transNon(:, iSym), type)
  end subroutine symmetry_get_type

end module m_ab7_symmetry
!!***
