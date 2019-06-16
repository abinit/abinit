!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_lattice_harmonic_potential
!! NAME
!! m_lattice_harmonic_potential
!!
!! FUNCTION
!! This module contains an harmonic lattice potential. 
!!
!! Datatypes:
!!
!! Subroutines:
!! TODO: add this when F2003 doc style is determined.
!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2019 ABINIT group (hexu)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE


#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


module m_lattice_harmonic_potential
  use defs_basis
  use m_errors
  use m_abicore
  use m_abstract_potential, only: abstract_potential_t
  use m_spmat_coo, only: COO_mat_t
  use m_spmat_csr, only: CSR_mat_t
  implicit none

  private

  type, public, extends(abstract_potential_t) :: lattice_harmonic_potential_t
     integer :: natom
     real(dp), allocatable :: ref_xcart(:,:)

  end type lattice_harmonic_potential_t

end module m_lattice_harmonic_potential
