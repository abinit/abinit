program have_fc_int_quad

  implicit none

  character(len=500) :: errmsg, dummy

  open(unit=11, file="iomsg.demo", status="NEW", iomsg=errmsg)
  write(11,iomsg=errmsg) "first"
  read(11, iomsg=errmsg) dummy
  close(unit=11, iomsg=errmsg)

end program have_fc_int_quad
