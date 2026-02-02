!!****p* ABINIT/testdiago
!! NAME
!! testdiago
!!
!! FUNCTION
!! Testing algorithms related to diagonalisation of Hermitian matrices.
!!
!! COPYRIGHT
!! Copyright (C) 2017-2026 ABINIT group (IML)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!!
!! SOURCE

program testdiago

    use iso_fortran_env, only: output_unit
    implicit none

    !Local variables-------------------------------
    !scalars
    real :: xTx
    !arrays
    real :: x(10)

!**********************************************************************

    call random_number(x) 
    xTx = dot_product(x, x)
    write(output_unit,*) 'testdiago OK', xTx

!**********************************************************************

 end program testdiago
!!***
