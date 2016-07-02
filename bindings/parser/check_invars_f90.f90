program toto

  use defs_basis
  use m_ab7_invars

  implicit none

  integer :: id, errno

  integer :: natom, ndtset, idtset, nsym, i
  double precision :: rprimd(3, 3)
  double precision, allocatable :: coord(:)
  integer :: dims(7), ndims
  character(len = 256) :: filename

  call getarg(1, filename)

  call ab7_invars_new_from_file(id, filename, len(filename))

  call ab7_invars_get_ndtset(id, ndtset, errno)
  if (errno /= AB7_NO_ERROR) goto 1000
  do idtset = 0, ndtset, 1
     call ab7_invars_get_integer(id, natom, ab7_invars_natom, idtset, errno)
     if (errno /= AB7_NO_ERROR) goto 1000

     call ab7_invars_get_shape(id, dims, ndims, ab7_invars_xred_orig, idtset, errno)
     if (errno /= AB7_NO_ERROR) goto 1000
     allocate(coord(product(dims(1:ndims))))
     call ab7_invars_get_real_array(id, coord, size(coord), ab7_invars_xred_orig, idtset, errno)
     if (errno /= AB7_NO_ERROR) goto 1000

     call ab7_invars_get_real_array(id, rprimd, 9, ab7_invars_rprimd_orig, idtset, errno)
     if (errno /= AB7_NO_ERROR) goto 1000

     call ab7_invars_get_integer(id, nsym, ab7_invars_nsym, idtset, errno)
     if (errno /= AB7_NO_ERROR) goto 1000

     write(*, "(A,I0,A,I0,A)") "### DATASET ", idtset, "/", ndtset, " ###"
     write(*, "(A,I0,A,I0)") "Number of atoms in dataset ", idtset, ": ", natom
     write(*, "(A,3F12.6,A)") "box definition: (", rprimd(:, 1), ")"
     write(*, "(A,3F12.6,A)") "                (", rprimd(:, 2), ")"
     write(*, "(A,3F12.6,A)") "                (", rprimd(:, 3), ")"
     write(*, "(A,I0,A,I0)") "Size of coordiantes array in dataset ", idtset, ": ", size(coord)
     write(*, "(A,I0,A)") "Coordinates in dataset ", idtset, ":"
     do i = 0, size(coord) / 3 - 1, 1
        write(*, "(3F12.6)") coord(i * 3 + 1), coord(i * 3 + 2), coord(i * 3 + 3)
     end do
     write(*, "(A,I0,A,I0)") "Number of symmetries in dataset ", idtset, ": ", nsym
     write(*,*)

     deallocate(coord)
  end do

  1000 continue
  call ab7_invars_free(id)
  if (allocated(coord)) deallocate(coord)
  if (errno /= AB7_NO_ERROR) then
     write(0, *) "Error!", errno
  end if

end program toto
