program have_mkl_omatadd

  implicit none

  call mkl_somatadd
  call mkl_domatadd
  call mkl_comatadd
  call mkl_zomatadd

end program have_mkl_omatadd
