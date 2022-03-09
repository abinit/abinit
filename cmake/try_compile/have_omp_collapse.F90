program test_omp_collapse
  implicit none
  integer :: data(100)
  integer :: i,j

  !$omp parallel do collapse(2)
  do j=1,10
    do i=1,10
      data((i-1)+(j-1)*10) = i+j
    end do
  end do
  !$omp end parallel do

  !write(*,*) data(1), data(2), data(3), data(4), data(5)

end program test_omp_collapse
