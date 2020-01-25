!!****m* ABINIT/m_multibinit_io_xml
!! NAME
!!  m_multibinit_io_xml
!!
!!
!! FUNCTION
!!   This module contains interface to the C functions for the potential xml file parser.
!!
!! Subroutines:
!!
!! COPYRIGHT
!! Copyright (C) 2001-2020 ABINIT group (hexu)
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

module m_multibinit_io_xml
  use iso_c_binding
  use m_mathfuncs
  use defs_basis
  use m_abicore
  use m_errors
  implicit none

!!*** 

  interface
     !-------------------------------------------------------------------!
     ! Read the spin xml file and save to the variables.
     ! C function:
     ! void xml_read_spin(char *fname, double *ref_energy, double *unitcell[9],
     ! int *natoms, double *masses[], int *nspin,
     ! int *index_spin[], double *gyroratios[], double *damping_factors[],
     ! double *positions[], double *spinat[],
     ! // exchange
     ! int *exc_nnz, int *exc_ilist[],
     ! int *exc_jlist[], int *exc_Rlist[],
     ! double *exc_vallist[],
     ! //dmi
     ! int *dmi_nnz, int *dmi_ilist[],
     ! int *dmi_jlist[], int *dmi_Rlist[],
     ! double *dmi_vallist[],
     ! //uniaxial SIA
     ! int *uni_nnz, int *uni_ilist[],
     ! double *uni_amplitude_list[],
     ! double *uni_direction_list[],
     ! //bilinear
     ! int *bi_nnz, int *bi_ilist[],
     ! int *bi_jlist[], int *bi_Rlist[],
     ! double *bi_vallist[])
     subroutine xml_read_spin(xml_fname, ref_energy, unitcell,                 &
          natoms, masses, nspin, index_spin, gyroratios, damping_factors, positions, spinat, &
          exc_nnz, exc_ilist, exc_jlist, exc_Rlist, exc_vallist, &
          dmi_nnz, dmi_ilist, dmi_jlist, dmi_Rlist, dmi_vallist, &
          uni_nnz, uni_ilist, uni_amplitude_list, uni_direction_list, &
          bi_nnz, bi_ilist, bi_jilst, bi_Rlist, bi_vallist) bind(C, name="xml_read_spin")
       import
       character(c_char), intent(in) :: xml_fname(*)
       integer (c_int), intent(out):: natoms, nspin, exc_nnz, dmi_nnz, uni_nnz, bi_nnz
       real  (c_double), intent(out) :: ref_energy
       type(c_ptr)::  unitcell,  &
            masses,  index_spin, gyroratios, damping_factors, positions, spinat, &
            exc_ilist, exc_jlist, exc_Rlist, exc_vallist, &
            dmi_ilist, dmi_jlist, dmi_Rlist, dmi_vallist, &
            uni_ilist, uni_amplitude_list, uni_direction_list, &
            bi_ilist, bi_jilst, bi_Rlist, bi_vallist
     end subroutine xml_read_spin
 end interface

 interface
    !-------------------------------------------------------------------!
    ! Interface to the C code which free the memory.
    !-------------------------------------------------------------------!
     subroutine xml_free_spin(xml_fname, ref_energy, unitcell,                 &
          natoms, masses, nspin, index_spin, gyroratios, damping_factors, positions, spinat, &
          exc_nnz, exc_ilist, exc_jlist, exc_Rlist, exc_vallist, &
          dmi_nnz, dmi_ilist, dmi_jlist, dmi_Rlist, dmi_vallist, &
          uni_nnz, uni_ilist, uni_amplitude_list, uni_direction_list, &
          bi_nnz, bi_ilist, bi_jilst, bi_Rlist, bi_vallist) bind(C, name="xml_free_spin")
       import
       character(c_char), intent(in) :: xml_fname(*)
       integer (c_int), intent(out):: natoms, nspin, exc_nnz, dmi_nnz, uni_nnz, bi_nnz
       real  (c_double), intent(out) :: ref_energy
       type(c_ptr)::  unitcell,  &
            masses,  index_spin, gyroratios, damping_factors, positions, spinat, &
            exc_ilist, exc_jlist, exc_Rlist, exc_vallist, &
            dmi_ilist, dmi_jlist, dmi_Rlist, dmi_vallist, &
            uni_ilist, uni_amplitude_list, uni_direction_list, &
            bi_ilist, bi_jilst, bi_Rlist, bi_vallist
     end subroutine xml_free_spin
  end interface

end module m_multibinit_io_xml
