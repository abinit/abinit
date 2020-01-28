!!****m* ABINIT/m_pseudo_pwscf
!! NAME
!!  m_pseudo_pwscf
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (MVer)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module pseudo_pwscf

  implicit none
  !
  ! All variables to be read from the UPF file
  ! (UPF = unified pseudopotential format)
  !
  integer ,parameter :: npsx = 6
  ! npsx  : maximum number of different pseudopotentials
  integer, parameter :: lmaxx  = 3, nchix  = 6, ndm = 2000
  ! lmaxx : maximum non local angular momentum in PP      
  ! nchix : maximum number of atomic wavefunctions per PP
  ! ndm   : maximum number of points in the radial mesh
  integer, parameter :: nbrx = 8, lqmax = 5, nqfx = 8
  ! nbrx  : maximum number of beta functions         
  ! lqmax : maximum number of angular momentum of Q  
  ! nqfx  : maximum number of coefficients in Q smoothing
  !
  ! pp_header
  character (len=80):: generated, date_author, comment
  character (len=2) :: psd(npsx), pseudotype
  character (len=20):: dft(npsx)
  integer :: lmax(npsx), mesh(npsx), nbeta(npsx), ntwfc(npsx)
  logical :: nlcc(npsx), isus(npsx)
  real(8) :: zp(npsx), ecutrho, ecutwfc, etotps
  real(8) :: oc(nchix,npsx)
  character(len=2) :: els(nchix,npsx)
  integer :: lchi(nchix,npsx)
  !
  ! pp_mesh
  real(8) :: r(ndm,npsx), rab(ndm,npsx)
  !   pp_nlcc
  real(8) :: rho_atc(ndm,npsx)
  !
  ! pp_local
  real(8) ::  vloc0(ndm,npsx)
  !
  ! pp_nonlocal
  ! pp_beta
  real(8) :: betar(ndm, nbrx, npsx)
  integer :: lll(nbrx,npsx), ikk2(nbrx,npsx)  
  ! pp_dij
  real(8) :: dion(nbrx,nbrx,npsx)
  ! pp_qij
  integer ::  nqf(npsx), nqlc(npsx)
  real(8) :: rinner(lqmax,npsx), qqq(nbrx,nbrx,npsx), &
       qfunc(ndm,nbrx,nbrx,npsx)
  ! pp_qfcoef
  real(8) :: qfcoef(nqfx,lqmax,nbrx,nbrx,npsx)
  !
  ! pp_pswfc
  real(8) :: chi(ndm,nchix,npsx)
  !
  ! pp_rhoatom
  real(8) :: rho_at(ndm,npsx)

end module pseudo_pwscf
!!***
