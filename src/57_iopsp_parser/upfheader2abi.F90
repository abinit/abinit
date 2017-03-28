!{\src2tex{textfont=tt}}
!!****f* ABINIT/upfheader2abi
!! NAME
!! upfheader2abi
!!
!! FUNCTION
!!  This routine wraps a call to a PWSCF module, which reads in
!!  a UPF (PWSCF / Espresso) format pseudopotential, then transfers
!!  data for the HEADER of abinit psps only!
!!
!! COPYRIGHT
!! Copyright (C) 2009-2017 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  filpsp = name of file with UPF data
!!
!! OUTPUT
!!  pspxc = index of xc functional for this pseudo
!!  lmax_ = maximal angular momentum
!!  znucl = charge of species nucleus
!!  zion = valence charge
!!  n1xccc = default number of points. Set to 0 if no nlcc is present
!!  nproj_l= number of projectors for each channel
!!  nprojso_l= number of projectors for each channel for SO correction projectors
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      inpspheads
!!
!! CHILDREN
!!      set_dft_from_indices,set_dft_from_name
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine upfheader2abi (filpsp, znucl, zion, pspxc, lmax_, n1xccc, nproj_l, nprojso_l)

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_atomdata
 use pseudo_pwscf ! pwscf module with all data explicit!

 use m_io_tools,  only : open_file

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'upfheader2abi'
 use interfaces_11_qespresso_ext
 use interfaces_57_iopsp_parser, except_this_one => upfheader2abi
!End of the abilint section

  implicit none

!Arguments -------------------------------
  character(len=fnlen), intent(in) :: filpsp
  integer, intent(inout) :: n1xccc 
  integer, intent(out) :: pspxc, lmax_
  real(dp), intent(out) :: znucl, zion
  !arrays
  integer, intent(out) :: nproj_l(0:3)
  integer, intent(out) :: nprojso_l(1:3)

!Local variables -------------------------
  integer :: iproj, ll, iunit
  character(len=500) :: msg
  type(atomdata_t) :: atom

! *********************************************************************

!call pwscf routine for reading in UPF
 if (open_file(filpsp, msg, newunit=iunit, status='old',form='formatted') /= 0) then
   MSG_ERROR(msg)
 end if

!read in psp data to static data in pseudo module, for ipsx == 1
 call read_pseudo(1,iunit)
 close (iunit)

!copy over to abinit internal arrays and vars
 call upfxc2abi(dft(1), pspxc)
 lmax_ = lmax(1)
 call atomdata_from_symbol(atom,psd(1))
 znucl = atom%znucl
 zion = zp(1)

 nproj_l = 0
 do iproj = 1, nbeta(1)
   ll = lll(iproj,1)
   nproj_l(ll) = nproj_l(ll) + 1
 end do

 nprojso_l = 0 !FIXME deal with so
!do iproj = 1, nbeta(1)
!nprojso_l(ll+1) = nprojso_l(ll+1) + 1
!end do

 if (.not. nlcc(1)) n1xccc = 0

end subroutine upfheader2abi 
!!***

!!****f* ABINIT/upfxc2abi
!! NAME
!! upfxc2abi
!!
!! FUNCTION
!!  This routine wraps a call to an OCTOPUS module, which reformats
!!  a UPF (PWSCF / Espresso) string describing XC functionals,
!!  and returns the abinit internal code pspxc
!!
!! COPYRIGHT
!! Copyright (C) 2009-2017 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dft = string with x/c functionals from PWSCF format
!!
!! OUTPUT
!!  pspxc = index of xc functional for this pseudo
!!
!! SIDE EFFECTS
!!
!! NOTES
!!   FIXME: extend to more functionals with libxc
!!   Could be included in separate module, eg read_upf_pwscf or funct_pwscf
!!   Left without defs_basis or calls to abinit routines ON PURPOSE
!!
!! PARENTS
!!      upf2abinit,upfheader2abi
!!
!! CHILDREN
!!      set_dft_from_indices,set_dft_from_name
!!
!! SOURCE
subroutine upfxc2abi(dft, pspxc)

 use defs_basis
 use funct_pwscf ! pwscf module for naming xc functionals
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'upfxc2abi'
!End of the abilint section

  implicit none
 
!Arguments -------------------------------
  character(len=20), intent(in) :: dft
  integer, intent(out) :: pspxc

!Local variables -------------------------
  integer :: iexch,icorr,igcx,igcc
  integer :: totalindex, offset

! *********************************************************************
!extract from char*20 :: dft(:)
!###  The following has been copied from pwscf src/Modules/upf_to_internal.f90:
!workaround for rrkj format - it contains the indices, not the name
 if ( dft(1:6)=='INDEX:') then
   read( dft(7:10), '(4i1)') iexch,icorr,igcx,igcc
   call set_dft_from_indices(iexch,icorr,igcx,igcc)
 else
   call set_dft_from_name( dft )
   iexch = get_iexch()
   icorr = get_icorr()
   igcx = get_igcx()
   igcc = get_igcc()
 end if
!reset dft string to avoid stray spaces
 call set_dft_from_indices(iexch,icorr,igcx,igcc)
 write(std_out,'(a)') ' upf2abinit: XC string from pseudopotential is :'
 write(std_out,'(3a)') '>', dft, '<'

 offset = 100
 totalindex = offset*offset*offset*iexch + offset*offset*icorr + offset*igcx + igcc
 select case (totalindex)
 case (00000000)  !(" NOX  NOC NOGX NOGC") ! no xc
   pspxc = 0  
 case (01010000)  !(" SLA   PZ NOGX NOGC") ! slater exchange + Perdew Zunger
   pspxc = 2  
 case (01050000)  !(" SLA  WIG NOGX NOGC") ! slater exchange + Wigner corr
   pspxc = 4  
 case (01060000)  !(" SLA   HL NOGX NOGC") ! Hedin + Lundqvist
   pspxc = 5  
 case (02000000)  !(" SL1  NOC NOGX NOGC") ! full slater exchange
   pspxc = 6  
 case (01040000)  !(" SLA   PW NOGX NOGC") ! slater exchange + Perdew Wang
   pspxc = 7  
 case (01000000)  !(" SLA  NOC NOGX NOGC") ! Perdew Wang + no corr
   pspxc = 8  
 case (01040304)  !(" SLA   PW  PBX  PBC") ! LDA + PBE GGA
   pspxc = 11 ! PBE
 case (01000300)  !(" SLA  NOC  PBX NOGC") ! exchange part of PBE GGA
   pspxc = 12 
 case (01040404)  !(" SLA   PW  RPB  PBC") ! rev PBE
   pspxc = 14 
 case (00000505)  !(" NOX  NOC HTCH HTCH") ! HTCH 120
   pspxc = 17 
 case (01030103)  !(" SLA  LYP  B88 BLYP") ! BLYP 
   pspxc = -106131
 case (01040101)  !(" SLA   PW  B88  P86") ! BP86
   pspxc = -106132
 case (00030603)  !(" NOX  LYP OPTX BLYP") ! OLYP
   pspxc = -110131
!    FIXME: important cases left to be patched with libxc:
!    vosko wilkins nusair
!    ortiz ballone
!    pbe0
!    Gunnarson-Lunqvist
!    make general approach: check gradient parts first, then lda.
!    event. check if they are consistent.
 case default
   MSG_ERROR('upf2abinit: XC functional not recognized')
 end select

end subroutine upfxc2abi
!!***
