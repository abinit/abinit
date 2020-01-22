!!****m* ABINIT/defs_wannier90
!! NAME
!! defs_wannier90
!!
!! FUNCTION
!! This module contains interfaces for wannier90 lib.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2020 ABINIT group (BA, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module defs_wannier90

 implicit none

 interface
#ifndef HAVE_WANNIER90_V1
!WANNIER90_V2
 subroutine wannier_setup(seed__name,mp_grid_loc,num_kpts_loc,&
     real_lattice_loc,recip_lattice_loc,kpt_latt_loc,num_bands_tot, &
     num_atoms_loc,atom_symbols_loc,atoms_cart_loc,gamma_only_loc,spinors_loc,&
     nntot_loc,nnlist_loc,nncell_loc,num_bands_loc,num_wann_loc, &
     proj_site_loc,proj_l_loc,proj_m_loc,proj_radial_loc,proj_z_loc, &
     proj_x_loc,proj_zona_loc,exclude_bands_loc, &
&    proj_s_loc,proj_s_qaxis_loc)
#else
!WANNIER90_V1
  subroutine wannier_setup(seed__name,mp_grid_loc,num_kpts_loc,&
     real_lattice_loc,recip_lattice_loc,kpt_latt_loc,num_bands_tot, &
     num_atoms_loc,atom_symbols_loc,atoms_cart_loc,gamma_only_loc,spinors_loc,&
     nntot_loc,nnlist_loc,nncell_loc,num_bands_loc,num_wann_loc, &
     proj_site_loc,proj_l_loc,proj_m_loc,proj_radial_loc,proj_z_loc, &
     proj_x_loc,proj_zona_loc,exclude_bands_loc)
#endif
  use defs_basis
  implicit none
  character(len=*), intent(in) :: seed__name
  integer, dimension(3), intent(in) :: mp_grid_loc
  integer, intent(in) :: num_kpts_loc
  real(dp), dimension(3,3), intent(in) :: real_lattice_loc
  real(dp), dimension(3,3), intent(in) :: recip_lattice_loc
  real(dp), dimension(3,num_kpts_loc), intent(in) :: kpt_latt_loc
  integer, intent(in) :: num_bands_tot
  integer, intent(in) :: num_atoms_loc
  character(len=*), dimension(num_atoms_loc), intent(in) :: atom_symbols_loc
  real(dp), dimension(3,num_atoms_loc), intent(in) :: atoms_cart_loc
  logical, intent(in) :: gamma_only_loc
  logical, intent(in) :: spinors_loc
  integer, intent(out) :: nntot_loc
  integer, dimension(num_kpts_loc,12), intent(out) :: nnlist_loc
  integer,dimension(3,num_kpts_loc,12), intent(out) :: nncell_loc
  integer, intent(out) :: num_bands_loc
  integer, intent(out) :: num_wann_loc
  real(dp), dimension(3,num_bands_tot), intent(out) :: proj_site_loc
  integer, dimension(num_bands_tot), intent(out) :: proj_l_loc
  integer, dimension(num_bands_tot), intent(out) :: proj_m_loc
  integer, dimension(num_bands_tot), intent(out) :: proj_radial_loc
  real(dp), dimension(3,num_bands_tot), intent(out) :: proj_z_loc
  real(dp), dimension(3,num_bands_tot), intent(out) :: proj_x_loc
  real(dp), dimension(num_bands_tot), intent(out) :: proj_zona_loc
  integer, dimension(num_bands_tot), intent(out) :: exclude_bands_loc
#ifndef HAVE_WANNIER90_V1
  !WANNIER90_V2
  integer, dimension(num_bands_tot), optional, intent(out) :: proj_s_loc  
  real(dp), dimension(3,num_bands_tot), optional, intent(out) :: proj_s_qaxis_loc
#endif
  end subroutine wannier_setup
 end interface

 interface
  subroutine wannier_run(seed__name,mp_grid_loc,num_kpts_loc, &
     real_lattice_loc,recip_lattice_loc,kpt_latt_loc,num_bands_loc, &
     num_wann_loc,nntot_loc,num_atoms_loc,atom_symbols_loc, &
     atoms_cart_loc,gamma_only_loc,M_matrix_loc,A_matrix_loc,eigenvalues_loc, &
     U_matrix_loc,U_matrix_opt_loc,lwindow_loc,wann_centres_loc, &
     wann_spreads_loc,spread_loc)
   use defs_basis
   character(len=*), intent(in) :: seed__name
   integer, dimension(3), intent(in) :: mp_grid_loc
   integer, intent(in) :: num_kpts_loc
   real(dp), dimension(3,3), intent(in) :: real_lattice_loc
   real(dp), dimension(3,3), intent(in) :: recip_lattice_loc
   real(dp), dimension(3,num_kpts_loc), intent(in) :: kpt_latt_loc
   integer, intent(in) :: num_bands_loc
   integer, intent(in) :: num_wann_loc
   integer, intent(in) :: nntot_loc
   integer, intent(in) :: num_atoms_loc
   character(len=*), dimension(num_atoms_loc), intent(in) :: atom_symbols_loc
   real(dp), dimension(3,num_atoms_loc), intent(in) :: atoms_cart_loc
   logical, intent(in)::gamma_only_loc
   complex(dpc), dimension(num_bands_loc,num_bands_loc,nntot_loc,num_kpts_loc), intent(in) :: M_matrix_loc
   complex(dpc), dimension(num_bands_loc,num_wann_loc,num_kpts_loc), intent(in) :: A_matrix_loc
   real(dp), dimension(num_bands_loc,num_kpts_loc), intent(in) :: eigenvalues_loc
   complex(dpc), dimension(num_wann_loc,num_wann_loc,num_kpts_loc), intent(out) :: U_matrix_loc
   complex(dpc), dimension(num_bands_loc,num_wann_loc,num_kpts_loc), optional, intent(out) :: U_matrix_opt_loc
   logical, dimension(num_bands_loc,num_kpts_loc), optional, intent(out) :: lwindow_loc
   real(dp), dimension(3,num_wann_loc), optional, intent(out) :: wann_centres_loc
   real(dp), dimension(num_wann_loc), optional, intent(out) :: wann_spreads_loc
   real(dp), dimension(3), optional, intent(out) :: spread_loc
  end subroutine wannier_run
 end interface

end module defs_wannier90
!!***
