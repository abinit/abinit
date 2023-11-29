program have_fc_iso_fortran_env

  use ISO_FORTRAN_ENV, only : int16,int32,int64

  implicit none

  integer(kind=int64) :: data

  data = 42

  write(*,*) "data=",data

end program have_fc_iso_fortran_env
