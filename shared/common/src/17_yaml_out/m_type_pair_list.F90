!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_type_pair_list
!! NAME
!!  m_pair_list
!!
!! FUNCTION
!!  This module define the interface to the C implementation
!!  of pair_list. See m_pair_list.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2019 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! PARENTS
!!   m_pair_list
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module m_type_pair_list
  use iso_c_binding
! type c_pair_list
! Represent a list of key-value pairs value can be either
! integer, double precision real or string.
! Never manipulate member directly or you will corrupt
! data. Use routines of m_pair_list instead
  type, bind(C) :: c_pair_list
    type(c_ptr) :: first = C_NULL_PTR
    type(c_ptr) :: cursor = C_NULL_PTR
    integer(c_int) :: length = 0
  end type c_pair_list
  contains
end module m_type_pair_list
