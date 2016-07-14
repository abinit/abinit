!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_libxml
!! NAME
!!  m_libxml
!!
!! FUNCTION
!!  Module containing interfaces to the LibXML library
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

module m_libxml

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_effective_potential
 use iso_c_binding

 implicit none
 private

#if HAVE_LIBXML
 public :: libxml_checkXML
 public :: libxml_readEffPot
 public :: libxml_getValue
 public :: libxml_getAttribute
 public :: libxml_getDimsEffPot
#endif
 
 public :: char_f2c
 public :: char_c2f

#if HAVE_LIBXML
 interface
   subroutine libxml_readEffPot(filename,natom,&
&     ntypat,nrpt,nph1l,amu,atmfrc,cell,dynmat,elastic_constants,energy,&
&     epsilon_inf,ewald_atmfrc,internal_strain,phfrq,rprimd,qph1l,&
&     short_atmfrc,typat,xcart,zeff)&
&                          bind(C,name="libxml_readEffPot")
     use iso_c_binding, only : C_CHAR,C_DOUBLE,C_INT
     integer(C_INT) :: natom,ntypat,nrpt,nph1l
     integer(C_INT) :: typat(natom)
     integer(C_INT) :: cell(nrpt,3)
     real(C_DOUBLE) :: energy
     real(C_DOUBLE) :: dynmat(2,3,natom,3,natom,nph1l)
     real(C_DOUBLE) :: internal_strain(6,natom,3)
     real(C_DOUBLE) :: phfrq(3*natom,nph1l),qph1l(3,nph1l)
     real(C_DOUBLE) :: atmfrc(2,3,natom,3,natom,nrpt)
     real(C_DOUBLE) :: short_atmfrc(2,3,natom,3,natom,nrpt)
     real(C_DOUBLE) :: ewald_atmfrc(2,3,natom,3,natom,nrpt)
     real(C_DOUBLE) :: amu(ntypat),rprimd(3,3),epsilon_inf(3,3)
     real(C_DOUBLE) :: zeff(3,3,natom)
     real(C_DOUBLE) :: elastic_constants(6,6),xcart(3,natom)
     character(kind=C_CHAR) :: filename(*)
   end subroutine libxml_readEffPot
 end interface

 interface
   subroutine libxml_getDimsEffPot(filename,natom,ntypat,nph1l,nrpt)&
&                          bind(C,name="libxml_getDimsEffPot")
     use iso_c_binding, only : C_CHAR,C_INT
     integer(C_INT) :: natom,ntypat,nph1l,nrpt
     character(kind=C_CHAR) :: filename(*)
   end subroutine libxml_getDimsEffPot
 end interface

 interface
   subroutine libxml_checkXML(filename,name_root) &
&                          bind(C,name="libxml_checkXML")
     use iso_c_binding, only : C_CHAR
     character(kind=C_CHAR) :: filename(*),name_root(*)
   end subroutine libxml_checkXML
 end interface

 interface
   subroutine libxml_getValue(filename,name_value,value_result) &
 &                          bind(C,name="libxml_getValue")
      use iso_c_binding, only : C_CHAR
      implicit none
      character(kind=C_CHAR) :: filename(*),name_value(*)
      character(kind=C_CHAR) :: value_result
    end subroutine libxml_getValue
  end interface
 
 interface
   subroutine libxml_getAttribute(filename,name_value,name_attribute) &
&                          bind(C,name="libxml_getAttribute")
     use iso_c_binding, only : C_CHAR
     character(kind=C_CHAR) :: filename(*),name_value(*),name_attribute(*)
   end subroutine libxml_getAttribute
 end interface

#endif
!!***
 
CONTAINS  !===========================================================================================

!!****f* m_libxml/char_f2c
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

!!****f* m_libxml/char_c2f
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

end module m_libxml


