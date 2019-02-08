module m_type_pair_list
  use iso_c_binding
! type pair_list
! Represent a list of key-value pairs
! value can be either integerof double precision real
! Never manipulate member directly or you will corrupt
! data. Use routnes of this module.
  type, bind(C) :: pair_list
    type(c_ptr) :: first = C_NULL_PTR
    type(c_ptr) :: cursor = C_NULL_PTR
    integer(c_int) :: length = 0
  end type pair_list
  contains
end module m_type_pair_list
