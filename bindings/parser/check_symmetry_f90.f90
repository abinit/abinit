program lalalalayou

  use defs_basis
  use m_ab7_symmetry
  use m_ab7_kpoints

  implicit none

  integer :: symObj, errno, i
  real(dp) :: rprimd(3, 3)
  real(dp) :: xRed(3, 2)
  real(dp) :: spinAt(3, 2)

  integer :: nSym, multiplicity
  integer :: sym(3, 3, AB7_MAX_SYMMETRIES)
  integer :: symAfm(AB7_MAX_SYMMETRIES)
  real(dp) :: transNon(3, AB7_MAX_SYMMETRIES)
  real(dp)         :: genAfm(3)
  character(len=15):: spaceGroup
  integer          :: spaceGroupId, pointGroupMagn
  integer :: equiv(4, AB7_MAX_SYMMETRIES)
  integer :: nkpt
  integer, parameter :: ngkpt(3) = (/ 2, 2, 2 /)
  real(dp), pointer :: kpt(:,:)
  real(dp), pointer :: wkpt(:)
  real(dp) :: shiftk(3, 1)

  rprimd = reshape((/ 0., 1., 1., 0.5, 0., 1., 0.5, 0.5, 0. /), (/ 3, 3 /))
  xRed = reshape((/ 0., 0., 0., 0.25, 0.25, 0.25 /), (/ 3, 2 /))
  spinAt = reshape((/ 0., 0., 1., 0., 0., -1. /), (/ 3, 2 /))

  call symmetry_new(symObj)
  call symmetry_set_lattice(symObj, rprimd, errno)
  call symmetry_set_structure(symObj, 2, (/ 1, 1 /), xRed, errno)
  call symmetry_set_spin(symObj, 2, spinAt, errno)

  call symmetry_get_matrices(symObj, nSym, sym, transNon, symAfm, errno)
  write(*,"(A,I3)") "nSym:", nSym
  do i = 1, nSym, 1
     write(*, "(A,I4.3,A)") "sym", i, ":"
     write(*, "(3I3,F12.6,I3)") sym(:, 1, i), transNon(1, i), symAfm(i)
     write(*, "(3I3,F12.6,I3)") sym(:, 2, i), transNon(2, i), symAfm(i)
     write(*, "(3I3,F12.6,I3)") sym(:, 3, i), transNon(3, i), symAfm(i)
  end do
  call symmetry_get_multiplicity(symObj, multiplicity, errno)
  write(*,"(A,I3)") "multiplicity:", multiplicity
  call symmetry_get_group(symObj, spaceGroup, &
       & spaceGroupId, pointGroupMagn, genAfm, errno)
  write(*, "(A,A,I6)") "space group:", trim(spaceGroup), spaceGroupId
  if (pointGroupMagn > 0) then
     write(*, "(3F12.6)") genAfm
  end if
  call symmetry_get_equivalent_atom(symObj, equiv, 1, errno)

  shiftk = reshape((/ 0.5, 0.5, 0.5 /), (/ 3, 1 /))
  call kpoints_get_mp_k_grid(symObj, nkpt, kpt, wkpt, ngkpt, 1, shiftk, errno)
  write(*,"(A,I3)") "k-points (MP 2x2x2):", nkpt
  do i = 1, nkpt, 1
     write(*, "(I3,A,3F10.6,F12.6)") i - 1, ":", kpt(:, i), wkpt(i)
  end do
  deallocate(kpt)
  deallocate(wkpt)

  call kpoints_get_auto_k_grid(symObj, nkpt, kpt, wkpt, 2._dp, errno)
  write(*,"(A,I3)") "k-points (kptrlen = 2):", nkpt
  do i = 1, nkpt, 1
     write(*, "(I3,A,3F10.6,F12.6)") i - 1, ":", kpt(:, i), wkpt(i)
  end do
  deallocate(kpt)
  deallocate(wkpt)

  call symmetry_free(symObj)
end program lalalalayou
