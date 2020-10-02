!!****m* ABINIT/m_libpaw_defs
!! NAME
!! m_libpaw_defs
!!
!! FUNCTION
!! Several definitions used in libPAW: named constants, physical constants, datatypes
!!
!! COPYRIGHT
!! Copyright (C) 2000-2020 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  This file comes directly from defs_basis.F90 module delivered with ABINIT.
!!
!!  FOR DEVELOPPERS: in order to preserve the portability of libPAW library,
!!  please consult ~abinit/src/??_libpaw/libpaw-coding-rules.txt
!!
!! SOURCE

module m_libpaw_defs

 implicit none

!Double precision real/complex subtypes
!-----------------------------------------------
 integer, parameter :: dp=kind(1.0d0)                   !Nb of bytes related to DP real numbers
 integer, parameter :: dpc=kind((1.0_dp,1.0_dp))        !Nb of bytes related to DP complex numbers

!Real constants
!-----------------------------------------------
 real(dp), parameter :: zero=0._dp
 real(dp), parameter :: one=1._dp
 real(dp), parameter :: two=2._dp
 real(dp), parameter :: three=3._dp
 real(dp), parameter :: four=4._dp
 real(dp), parameter :: half=0.50_dp
 real(dp), parameter :: third=one/three
 real(dp), parameter :: quarter=0.25_dp
 real(dp), parameter :: eighth=0.125_dp
 real(dp), parameter :: sqrt2=1.4142135623730950488016887242096939_dp
 real(dp), parameter :: sqrt3=1.7320508075688772935274463415058739_dp
 real(dp), parameter :: sqrthalf=0.70710678118654752440084436210484697_dp
 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
 real(dp), parameter :: two_pi=two*pi
 real(dp), parameter :: four_pi=four*pi
 real(dp), parameter :: tol3= 0.001_dp
 real(dp), parameter :: tol6= 0.000001_dp
 real(dp), parameter :: tol8= 0.00000001_dp
 real(dp), parameter :: tol9= 0.000000001_dp
 real(dp), parameter :: tol10=0.0000000001_dp
 real(dp), parameter :: tol12=0.000000000001_dp
 real(dp), parameter :: tol14=0.00000000000001_dp
 real(dp), parameter :: tol16=0.0000000000000001_dp

!Complex constants
!-----------------------------------------------
 complex(dpc), parameter :: czero=(0._dp,0._dp)         ! 0 (complex)
 complex(dpc), parameter :: cone =(1._dp,0._dp)         ! 1 (complex)

!Character constants
!-----------------------------------------------
 character(len=1), parameter :: ch10 = char(10)         ! carriage return
 integer, parameter :: fnlen=264                        ! maximum length of file name variables
 integer, parameter :: strlen=2000000                   ! maximum length of input string

!UNIX unit numbers
!-----------------------------------------------
 integer, save      :: ab_out= 7                        ! output file
 integer, save      :: std_out=6                        ! standard output
 integer, parameter :: std_err=0                        ! standard error
 integer, parameter :: tmp_unit=9,tmp_unit2=10          ! units for temporary files

!Real physical constants
!-----------------------------------------------
 real(dp), parameter :: Bohr_Ang=0.52917720859_dp       ! 1 Bohr, in Angstrom
 real(dp), parameter :: Ha_eV=27.21138386_dp            ! 1 Hartree, in eV
 real(dp), parameter :: InvFineStruct=137.035999679_dp  ! Inverse of fine structure constant
 real(dp), parameter :: FineStructureConstant2=0.000053251354478_dp ! Square of fine structure

!A collection of small datatypes for ragged arrays
!-----------------------------------------------
 type coeffi1_type                    !A small datatype for ragged integer 1D-arrays
  integer, allocatable :: value(:) 
 end type coeffi1_type
 type coeff1_type                     !A small datatype for ragged real 1D-arrays
  real(dp), allocatable :: value(:) 
 end type coeff1_type
 type coeff2_type                     !A small datatype for ragged real 2D-arrays
  real(dp), allocatable :: value(:,:)  
 end type coeff2_type
 type coeff3_type                     !A small datatype for ragged real 3D-arrays
  real(dp), allocatable :: value(:,:,:) 
 end type coeff3_type

!Small functions used in cpp macros
!-----------------------------------------------
 public :: to_array
 interface to_array
  module procedure to_array1
  module procedure to_array2
  module procedure to_array3
  module procedure to_array4
  module procedure to_array5
  module procedure to_array6
 end interface to_array

 contains

 function to_array1(i1) result(arr)

  integer :: i1,arr(1)
  arr=(/i1/)
 end function to_array1

 function to_array2(i1,i2) result(arr)

  integer :: i1,i2,arr(2)
  arr=(/i1,i2/)
 end function to_array2

 function to_array3(i1,i2,i3) result(arr)

  integer :: i1,i2,i3,arr(3)
  arr=(/i1,i2,i3/)
 end function to_array3

 function to_array4(i1,i2,i3,i4) result(arr)

  integer :: i1,i2,i3,i4,arr(4)
  arr=(/i1,i2,i3,i4/)
 end function to_array4

 function to_array5(i1,i2,i3,i4,i5) result(arr)

  integer :: i1,i2,i3,i4,i5,arr(5)
  arr=(/i1,i2,i3,i4,i5/)
 end function to_array5

 function to_array6(i1,i2,i3,i4,i5,i6) result(arr)

  integer :: i1,i2,i3,i4,i5,i6,arr(6)
  arr=(/i1,i2,i3,i4,i5,i6/)
 end function to_array6

end module m_libpaw_defs
!!***
