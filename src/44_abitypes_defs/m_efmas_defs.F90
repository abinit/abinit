!!****m* ABINIT/m_efmas_defs
!! NAME
!! m_efmas_defs
!!
!! FUNCTION
!! This module contains datatypes for efmas functionalities.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2020 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_efmas_defs

 use defs_basis

 implicit none
 public
 !Put eventual module variables here...
!!***

!!****t* m_efmas_defs/efmasval_type
!! NAME
!! efmasval_type
!!
!! FUNCTION
!! The efmasval_type structured datatype, related to one band or one degenerated set of bands, for one k-point
!!
!! SOURCE

 type efmasval_type

   !For k-point
   complex(dpc),allocatable :: ch2c(:,:,:,:) ! ch2c(mdim,mdim,1:deg_dim,1:deg_dim) 
                                             ! where mdim=3 labels reciprocal space directions
                                             ! See Eq.(50) of Laflamme2016 : 2nd-order Hamiltonian contribution
                                             ! Two first indices are for number of directions
                                             ! Two last indices are for band indices within degenerate subspace
   complex(dpc),allocatable :: eig2_diag(:,:,:,:) ! eig2_diag(mdim,mdim,1:deg_dim,1:deg_dim) 
                                             ! where mdim=3 labels reciprocal space directions
                                             ! See Eq.(50) of Laflamme2016 : generalized second-order k-derivative

 end type efmasval_type
!!***

!!****t* m_efmas_defs/efmasdeg_type
!! NAME
!! efmasdeg_type
!!
!! FUNCTION
!! The efmasdeg_type structured datatype, related to one k-point
!!
!! SOURCE

 type efmasdeg_type

   !For k-point
   integer :: nband                           ! Number of bands (related to one specific k point)
   integer :: ndegs                           ! Number of (degenerate) sets of eigenvalues (related to one specific k point)
   integer, allocatable :: degs_bounds(:,:)   ! degs_bounds(2,ndegs) 
                                              ! Minimal and maximal band indices for each possibly degenerate set of eigenvalues
                                              ! actually the second dimension is declared as nband_k
   !For band
   integer :: deg_range(2)                    ! Indices of the sets that corresponds to the interval of bands for which
                                              ! the generalized second-order k-derivative eig2_diag is computed,
                                              ! possibly extended due to the degeneracies.
   integer,allocatable :: ideg(:)             ! ideg(nband_k)  index of the set to which a particular band belongs

 end type efmasdeg_type
!!***

end module m_efmas_defs
!!***
