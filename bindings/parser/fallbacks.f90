subroutine wrtout(unit,message,paral)

  use defs_basis

  implicit none

  !Arguments ------------------------------------
  integer,intent(in) :: unit
  character(len=500),intent(inout) :: message
  character(len=4),intent(in) :: paral

  if (unit == ab_out) then
     write(unit, *) "[AB] ", trim(paral), " - ", trim(message)
  end if
end subroutine wrtout

subroutine leave_new()

  implicit none

  write(0, *) "Sorry, I'm dead."
  write(0, *) "R.I.P."
  stop
end subroutine leave_new

subroutine timab(nn,option,tottim)

  use defs_basis
  use defs_time

  implicit none

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: nn,option
  !arrays
  real(dp),intent(out) :: tottim(2)

  tottim = real((/nn, option /), dp)
end subroutine timab

