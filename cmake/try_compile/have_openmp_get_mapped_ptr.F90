program test_openmp_get_mapped_ptr
        use omp_lib
        use, intrinsic :: iso_c_binding, only : c_ptr
        integer :: device_id
        type(c_ptr) :: ptr,ptr2

        device_id = omp_get_default_device()

        ptr2 = omp_get_mapped_ptr(ptr, device_id)

end program test_openmp_get_mapped_ptr
