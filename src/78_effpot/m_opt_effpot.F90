!!****m* ABINIT/m_opt_effpot
!!
!! NAME
!! m_opt_effpot
!!
!! FUNCTION
!!      
!! 
!!
!! COPYRIGHT
!! Copyright (C) 2010-2018 ABINIT group (AM)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE



#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_opt_effpot

 use defs_basis
 use m_errors
 use m_abicore
 use m_polynomial_coeff
 use m_atomdata
 use m_xmpi
 use m_supercell

 use m_special_funcs,only : factorial
 use m_geometry,       only : xred2xcart
 use m_crystal,only : symbols_crystal
 use m_strain,only : strain_type,strain_get
 use m_effective_potential,only : effective_potential_type, effective_potential_evaluate
 use m_effective_potential,only : effective_potential_freeCoeffs,effective_potential_setCoeffs
 use m_effective_potential,only : effective_potential_getDisp, effective_potential_writeAnhHead
 use m_effective_potential_file, only : effective_potential_file_mapHistToRef
 use m_io_tools,   only : open_file,get_unit
 use m_abihist, only : abihist,abihist_free,abihist_init,abihist_copy,write_md_hist,var2hist
 use m_random_zbq
 use m_fit_data
 use m_geometry, only: metric 

 implicit none

public :: opt_effpot 

!!****
CONTAINS 
      
!!****f* m_opt_effpot/opt_effpot 
!!
!! NAME
!! opt_effpot
!!
!! FUNCTION
!! Optimize Effective Potential by fitting the value of certain
!! coefficients while keeping the values of the others     
!!
!! INPUTS
!! eff_pot<type(effective_potential)> = effective potential
!!
!! opt_coeff(opt_ncoeff) = list of terms whose coefficients are to be
!! optimized
!!
!! hist<type(abihist)> = Training set Data(or snapshot of DFT)
!! comm = MPI communicator
!!
!! OUTPUT
!! eff_pot<type(effective_potential)> = effective potential datatype with new fitted coefficients
!!
!! PARENTS
!! multibinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine opt_effpot(eff_pot,opt_coeff,hist,comm) 

 implicit none  

!Arguments ------------------------------------
!scalars

!arrays 

!Logicals 
!Strings 
!Local variables ------------------------------
!scalars

!arrays 

!Logicals 
!Strings 
! *************************************************************************



write(*,*) "I was here, nothing to optimize yet"
      
end subroutine 


end module m_opt_effpot

      
