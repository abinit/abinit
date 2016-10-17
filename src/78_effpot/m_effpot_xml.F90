!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_effpot_xml
!! NAME
!!  m_effpot_xml
!!
!! FUNCTION
!!  Module containing interfaces to the Effpot_xml library
!!
!! COPYRIGHT
!! Copyright (C) 2015-2015 ABINIT group (MO, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!

!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_effpot_xml

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_effective_potential
 use iso_c_binding

 implicit none
 private

#if defined HAVE_LIBXML
 public :: effpot_xml_checkXML
 public :: effpot_xml_read
 public :: effpot_xml_getValue
 public :: effpot_xml_getAttribute
 public :: effpot_xml_getDims
#endif
 
 public :: char_f2c
 public :: char_c2f

#if defined HAVE_LIBXML
 interface
   subroutine effpot_xml_read(filename,natom,&
&     ntypat,nrpt,nqpt,amu,atmfrc,cell,dynmat,elastic_constants,energy,&
&     epsilon_inf,ewald_atmfrc,internal_strain,phfrq,rprimd,qph1l,&
&     short_atmfrc,typat,xcart,zeff)&
&                          bind(C,name="effpot_xml_read")
     use iso_c_binding, only : C_CHAR,C_DOUBLE,C_INT
     integer(C_INT) :: natom,ntypat,nrpt,nqpt
     integer(C_INT) :: typat(natom)
     integer(C_INT) :: cell(3,nrpt)
     real(C_DOUBLE) :: energy
     real(C_DOUBLE) :: dynmat(2,3,natom,3,natom,nqpt)
     real(C_DOUBLE) :: internal_strain(6,natom,3)
     real(C_DOUBLE) :: phfrq(3*natom,nqpt),qph1l(3,nqpt)
     real(C_DOUBLE) :: atmfrc(2,3,natom,3,natom,nrpt)
     real(C_DOUBLE) :: short_atmfrc(2,3,natom,3,natom,nrpt)
     real(C_DOUBLE) :: ewald_atmfrc(2,3,natom,3,natom,nrpt)
     real(C_DOUBLE) :: amu(ntypat),rprimd(3,3),epsilon_inf(3,3)
     real(C_DOUBLE) :: zeff(3,3,natom)
     real(C_DOUBLE) :: elastic_constants(6,6),xcart(3,natom)
     character(kind=C_CHAR) :: filename(*)
   end subroutine effpot_xml_read
 end interface

 interface
   subroutine effpot_xml_getDims(filename,natom,ntypat,nqpt,nrpt)&
&                          bind(C,name="effpot_xml_getDims")
     use iso_c_binding, only : C_CHAR,C_INT
     integer(C_INT) :: natom,ntypat,nqpt,nrpt
     character(kind=C_CHAR) :: filename(*)
   end subroutine effpot_xml_getDims
 end interface

 interface
   subroutine effpot_xml_checkXML(filename,name_root) &
&                          bind(C,name="effpot_xml_checkXML")
     use iso_c_binding, only : C_CHAR
     character(kind=C_CHAR) :: filename(*),name_root(*)
   end subroutine effpot_xml_checkXML
 end interface

 interface
   subroutine effpot_xml_getValue(filename,name_value,value_result) &
 &                          bind(C,name="effpot_xml_getValue")
      use iso_c_binding, only : C_CHAR
      implicit none
      character(kind=C_CHAR) :: filename(*),name_value(*)
      character(kind=C_CHAR) :: value_result
    end subroutine effpot_xml_getValue
  end interface
 
 interface
   subroutine effpot_xml_getAttribute(filename,name_value,name_attribute) &
&                          bind(C,name="effpot_xml_getAttribute")
     use iso_c_binding, only : C_CHAR
     character(kind=C_CHAR) :: filename(*),name_value(*),name_attribute(*)
   end subroutine effpot_xml_getAttribute
 end interface

#endif
!!***
 
CONTAINS  !===========================================================================================

!!****f* m_effpot_xml/char_f2c
!! NAME
!!  char_f_to_c
!!
!! FUNCTION
!! Helper function to convert a Fortran string to a C string
!! Based on a routine by Joseph M. Krahn
!!
!! INPUTS
!!  f_string=Fortran string
!!
!! OUTPUT
!!  c_string=C string
!!
!! SOURCE

function char_f2c(f_string) result(c_string)
!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'char_f2c'
!End of the abilint section

 character(len=*),intent(in) :: f_string
 character(kind=C_CHAR,len=1) :: c_string(len_trim(f_string)+1)
!Local variables -------------------------------
 integer :: ii,strlen
!! *************************************************************************
 strlen=len_trim(f_string)
 forall(ii=1:strlen)
   c_string(ii)=f_string(ii:ii)
 end forall
 c_string(strlen+1)=C_NULL_CHAR
end function char_f2c
!!***

!----------------------------------------------------------------------

!!****f* m_effpot_xml/char_c2f
!! NAME
!!  char_c_to_f
!!
!! FUNCTION
!! Helper function to convert a C string to a Fortran string
!! Based on a routine by Joseph M. Krahn
!!
!! INPUTS
!!  c_string=C string
!!
!! OUTPUT
!!  f_string=Fortran string
!!
!! PARENTS
!!      m_libpaw_libxc
!!
!! CHILDREN
!!
!! SOURCE

subroutine char_c2f(c_string,f_string)
!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'char_c2f'
!End of the abilint section

 character(kind=C_CHAR,len=1),intent(in) :: c_string(*)
 character(len=*),intent(out) :: f_string
!Local variables -------------------------------
 integer :: ii
!! *************************************************************************
 ii=1
 do while(c_string(ii)/=C_NULL_CHAR.and.ii<=len(f_string))
   f_string(ii:ii)=c_string(ii) ; ii=ii+1
 end do
 if (ii<len(f_string)) f_string(ii:)=' '
end subroutine char_c2f
!!***

end module m_effpot_xml
