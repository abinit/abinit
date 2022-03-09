program have_fc_iso_c_binding

  use, intrinsic :: iso_c_binding

  implicit none

  integer(kind=c_int32_t) :: data

  data = 42

  write(*,*) "data=",data

end program have_fc_iso_c_binding
