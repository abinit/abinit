program have_fc_command_argument

  implicit none

  integer :: ii
  character(len=500) :: arg

  call get_command(arg)

  do ii=1,command_argument_count()
    call get_command_argument(ii, arg)
  end do


end program have_fc_command_argument
