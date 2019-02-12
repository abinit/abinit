module m_pair_list

  use iso_c_binding
  use m_type_pair_list
  implicit none

  integer,parameter :: TC_EMPTY=-2, TC_NOTFOUND=-1, TC_INT=0, TC_REAL=1, TC_STRING=2

  private
  public :: pair_list_set, pair_list_get, pair_list_free
  public :: pair_list_next, pair_list_look, pair_list_iter, pair_list_restart
  public :: pair_list
  public :: TC_EMPTY, TC_NOTFOUND, TC_INT, TC_REAL, TC_STRING

  type,public :: pair_list
    type(c_pair_list) :: plc
    contains
      procedure :: set => pair_list_set
      procedure :: get => pair_list_get
      procedure :: free => pair_list_free
      procedure :: next => pair_list_next
      procedure :: look => pair_list_look
      procedure :: iter => pair_list_iter
      procedure :: restart => pair_list_restart
      procedure :: length => pair_list_length
  end type

! -------------------------------------------------------------------------------
! -                                                                             -
! -                        Private C function binding                           -
! -                                                                             -
! -------------------------------------------------------------------------------
  interface

    subroutine pair_list_next_c(pl) bind(C, name="pair_list_next")
      use m_type_pair_list
      type(c_pair_list),intent(in) :: pl
    end subroutine pair_list_next_c

    subroutine pair_list_free_c(pl) bind(C, name="pair_list_free")
      use m_type_pair_list
      type(c_pair_list),intent(inout) :: pl
    end subroutine pair_list_free_c
    subroutine pair_list_seti(pl, key, i, len) bind(C, name="pair_list_seti")
      use m_type_pair_list
      type(c_pair_list) :: pl
      character(kind=c_char) :: key(*)
      integer(kind=c_int) :: i, len
    end subroutine pair_list_seti

    subroutine pair_list_setr(pl, key, r, len) bind(C, name="pair_list_setr")
      use m_type_pair_list
      type(c_pair_list) :: pl
      character(kind=c_char) :: key(*)
      integer(kind=c_int) :: len
      real(kind=c_double) :: r
    end subroutine pair_list_setr

    subroutine pair_list_sets(pl, key, s, len, len_s) bind(C, name="pair_list_sets")
      use m_type_pair_list
      type(c_pair_list) :: pl
      character(kind=c_char) :: key(*), s(*)
      integer(kind=c_int) :: len, len_s
      real(kind=c_double) :: r
    end subroutine pair_list_sets

    subroutine pair_list_get_c(pl, key, type_code, i, r, s, len, len_s) bind(C, name="pair_list_get_")
      use m_type_pair_list
      type(c_pair_list) :: pl
      character(kind=c_char) :: key(*), s(*)
      integer(kind=c_int) :: i, type_code, len, len_s
      real(kind=c_double) :: r
    end subroutine pair_list_get_c

    subroutine pair_list_look_c(pl, key, type_code, i, r, s, len, len_s) bind(C, name="pair_list_look_")
      use m_type_pair_list
      type(c_pair_list) :: pl
      integer(kind=c_int) :: type_code, i, len, len_s
      character(kind=c_char) :: key(len), s(len_s)
      real(kind=c_double) :: r
    end subroutine pair_list_look_c

  end interface

! -------------------------------------------------------------------------------
! -                                                                             -
! -                          Pure Fortran Wrapper                               -
! -                                                                             -
! -------------------------------------------------------------------------------
  contains

! pair_list_length
! arg
!   type(pair_list) :: pl
! get the number of pair stored in pl
  function pair_list_length(pl) result(length)
    class(pair_list),intent(in) :: pl
    integer :: length
    length = pl%plc%length
  end function

! pair_list_get
! arg
!   type(pair_list) :: pl
!   character(len=*) :: key
!   integer :: type_code, i
!   real :: r
! get the value associated with a key, only one of i and r is modified
! type_code is 0 if the value was an integer (and so that i is setted)
! type_code is 1 if the value was a real number (and so that r is setted)
! type_code is 2 if the value was a string (and so that s is setted)
! type_code is -1 if the key was not present (neither i nor r are setted)
! type_code is -2 if the list is empty (neither i nor r are setted)
  subroutine pair_list_get(pl, key, type_code, i, r, s)
    class(pair_list),intent(in) :: pl
    character(kind=c_char,len=*),intent(in) :: key, s
    integer(kind=c_int),intent(out) :: i, type_code
    real(kind=c_double),intent(out) :: r
    call pair_list_get_c(pl%plc, trim(key), type_code, i, r, s, len_trim(key), len(s))
  end subroutine pair_list_get

! pair_list_look
! arg
!   class(pair_list) :: pl
!   character(*) :: key
!   integer :: type_code
!   integer :: i
!   real :: r
! pair_list variables have a cursor wich point onto an arbitrary element
! of the list. pair_list_look allow to extract the key-value pair from
! that element
!
! If key is shorter than the actual key of the pair, only available space
! is used resulting in truncated key
! If key is longer than the actual key remaining space is filled with spaces
!
! type_code is 1 if the value was a real number (and so that r is setted)
! type_code is 0 if the value was an integer (and so that i is setted)
! type_code is -2 if the cursor is null (list is empty or end have been reached)
  subroutine pair_list_look(pl, key, type_code, i, r, s)
    use m_type_pair_list
    class(pair_list),intent(in) :: pl
    character(kind=c_char,len=*),intent(out) :: key, s
    integer(kind=c_int),intent(out) :: type_code, i
    real(kind=c_double),intent(out) :: r
    call pair_list_look_c(pl%plc, key, type_code, i, r, s, len(key), len(s))
  end subroutine pair_list_look

! pair_list_next
! arg
!   type(pair_list) :: pl
! have the cursor (cf: pair_list_look) moving forward of one element.
    subroutine pair_list_next(pl)
      class(pair_list),intent(in) :: pl
      call pair_list_next_c(pl%plc)
    end subroutine pair_list_next

! pair_list_free
! arg
!   type(pair_list) :: pl
! free memory occupied by the list (not the pair_list variable itself !)
! and reset the pair_list variable (it can be reused as an empty list) 
    subroutine pair_list_free(pl)
      class(pair_list),intent(inout) :: pl
      call pair_list_free_c(pl%plc)
    end subroutine pair_list_free

! pair_list_set
! arg:
!   class(pair_list): pl
!   character(*): key
! optional:
!   integer: i
!   real: r
! set a key-value par into the list. If the key is already presen, the
! corresponding pair is updated. If not the pair is created.
! Only one of i and r should be provided (i is the default if both are
! provided). Nothing happen if none of them are provided.
  subroutine pair_list_set(pl, key, i, r, s)
    class(pair_list),intent(in) :: pl
    character(len=*),intent(in) :: key
    integer,intent(in),optional :: i
    real(kind=c_double),intent(in),optional :: r
    character(len=*),intent(in),optional :: s
    if(present(i)) then
      call pair_list_seti(pl%plc, trim(key), i, len_trim(key))
    elseif(present(r)) then
        call pair_list_setr(pl%plc, trim(key), r, len_trim(key))
    elseif(present(s)) then
        call pair_list_sets(pl%plc, trim(key), s, len_trim(key), len_trim(s))
    end if
  end subroutine pair_list_set

! pair_list_restart
! arg
!   class(pair_list) :: pl
! have the cursor going back to the first element (cf: pair_list_next)
  subroutine pair_list_restart(pl)
    class(pair_list),intent(inout) :: pl
    pl%plc%cursor = pl%plc%first;
  end subroutine pair_list_restart

! pair_list_next
! arg
!   class(pair_list) :: pl
!   character(*) :: key
!   integer :: type_code
!   integer :: i
!   real :: r
! equivalent to pair_list_look followed by pair_list_next
  subroutine pair_list_iter(pl, key, type_code, i, r, s)
    class(pair_list),intent(in) :: pl
    character(len=*),intent(out) :: key
    integer,intent(out) :: type_code
    integer,intent(out) :: i
    real(kind=c_double),intent(out) :: r
    character(len=*),intent(out) :: s
    call pair_list_look(pl, key, type_code, i, r, s)
    if(type_code >= 0) then
      call pair_list_next_c(pl%plc);
    end if
  end subroutine pair_list_iter
end module m_pair_list
