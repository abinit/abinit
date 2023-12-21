program test_openmp_offload
  use omp_lib
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  integer :: alfa,i,levels,arr(7),device_id,rc
  logical :: is_bool
  type(c_ptr) :: ptr,ptr2

  levels = 7
  arr = 0
  is_bool = omp_is_initial_device()
  device_id = omp_get_default_device()
  rc = omp_target_is_present(ptr, device_id)
  ptr2 = omp_get_mapped_ptr(ptr, device_id)

  ! Not an actual GPU offload test, as users may not compile on host with GPU.
  ! Only compiler capabilities are checked.
  ! This program should compile with compilers supporting OpenMP v5 with usual
  ! OpenMP flags. Other compilers will treat "!$OMP TARGET" as syntax error.
  if (levels == 8) then
!$OMP TARGET ENTER DATA MAP(to:arr)
!$OMP TARGET PARALLEL DO COLLAPSE(2) PRIVATE(alfa,i) MAP(to:arr)
    do alfa=1,1
      do i=1,levels
        arr(i)=i
      end do
    end do
!$OMP END TARGET PARALLEL DO

!$OMP TARGET DATA USE_DEVICE_ADDR(arr)
    alfa=1
!$OMP END TARGET DATA
!$OMP TARGET EXIT DATA MAP(from:arr)
  end if

end program test_openmp_offload
