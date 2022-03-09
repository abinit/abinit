program have_netcdf_fortran

  use netcdf

  implicit none

  character(len=*), parameter :: path = "dummy"
  integer :: mode, ncerr, ncid
  ncerr = nf90_open(path,mode,ncid)

end program have_netcdf_fortran
