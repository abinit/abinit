!{\src2tex{textfont=tt}}
!!****f* ABINIT/dtfil_init_time
!!
!! NAME
!! dtfil_init_time
!!
!! FUNCTION
!! Inside the itimimage, iimage and itime loops (this is only needed for optdriver=0),
!! initialize the remaining parts of dtfil.
!!
!!
!! COPYRIGHT
!! Copyright (C) 2010-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! iapp=indicates the eventual suffix to be appended to the generic output root
!!         if 0 : no suffix to be appended (called directly from gstate)
!!         if positive : append "_TIM//iapp" (called from move or brdmin)
!!         if -1 : append "_TIM0" (called from brdmin)
!!         if -2, -3, -4, -5: append "_TIMA", ... ,"_TIMD", (called from move)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! dtfil=<type datafiles_type>infos about file names, file unit numbers
!!  (part of which were initialized previously)
!!
!! PARENTS
!!      gstate,mover
!!
!! CHILDREN
!!      fappnd
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dtfil_init_time(dtfil,iapp)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dtfil_init_time'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: iapp
 type(datafiles_type),intent(inout) :: dtfil

!Local variables-------------------------------
!scalars
 character(len=fnlen) :: filapp,filprot

!******************************************************************

 DBG_ENTER("COLL")

!--------------------------------------------------------
!Names based on dtfil%filnam_ds(4)+iapp

!Prepare the name of the auxiliary files DOS, EIG...
 call fappnd(filapp,dtfil%filnam_ds(4),iapp)
 dtfil%fnameabo_app=trim(filapp)
 dtfil%fnameabo_app_atmden_core=trim(filapp)//'_ATMDEN_CORE'
 dtfil%fnameabo_app_atmden_val=trim(filapp)//'_ATMDEN_VAL'
 dtfil%fnameabo_app_atmden_full=trim(filapp)//'_ATMDEN_FULL'
 dtfil%fnameabo_app_n_tilde=trim(filapp)//'_N_TILDE'
 dtfil%fnameabo_app_n_one=trim(filapp)//'_N_ONE'
 dtfil%fnameabo_app_nt_one=trim(filapp)//'_NT_ONE'
 dtfil%fnameabo_app_bxsf=trim(filapp)//'_BXSF'
 dtfil%fnameabo_app_cif=trim(filapp)//'.cif'
 dtfil%fnameabo_app_den=trim(filapp)//'_DEN'
 dtfil%fnameabo_app_dos=trim(filapp)//'_DOS'
 dtfil%fnameabo_app_eig=trim(filapp)//'_EIG'
 dtfil%fnameabo_app_elf=trim(filapp)//'_ELF'
 dtfil%fnameabo_app_elf_down=trim(filapp)//'_ELF_DOWN'
 dtfil%fnameabo_app_elf_up=trim(filapp)//'_ELF_UP'
 dtfil%fnameabo_app_fatbands=trim(filapp)//'_FATBANDS'
 dtfil%fnameabo_app_gden1=trim(filapp)//'_GDEN1'
 dtfil%fnameabo_app_gden2=trim(filapp)//'_GDEN2'
 dtfil%fnameabo_app_gden3=trim(filapp)//'_GDEN3'
 dtfil%fnameabo_app_geo=trim(filapp)//'_GEO'
 dtfil%fnameabo_app_kden=trim(filapp)//'_KDEN'
 dtfil%fnameabo_app_lden=trim(filapp)//'_LDEN'
 dtfil%fnameabo_app_nesting=trim(filapp)//'_NEST'
 dtfil%fnameabo_app_opt=trim(filapp)//'_OPT'
 dtfil%fnameabo_app_opt2=trim(filapp)//'_OPT2'
 dtfil%fnameabo_app_pawden=trim(filapp)//'_PAWDEN'
 dtfil%fnameabo_app_pot=trim(filapp)//'_POT'
 dtfil%fnameabo_app_stm=trim(filapp)//'_STM'
 dtfil%fnameabo_app_vclmb=trim(filapp)//'_VCLMB'
 dtfil%fnameabo_app_vha=trim(filapp)//'_VHA'
 dtfil%fnameabo_app_vhxc=trim(filapp)//'_VHXC'
 dtfil%fnameabo_app_vpsp=trim(filapp)//'_VPSP'
 dtfil%fnameabo_app_vxc=trim(filapp)//'_VXC'
 dtfil%fnameabo_app_wfk=trim(filapp)//'_WFK'
 dtfil%fnameabo_app_vha_1dm=trim(filapp)//'_VHA_1DM'
 dtfil%fnameabo_app_vclmb_1dm=trim(filapp)//'_VCLMB_1DM'
 dtfil%fnameabo_app_1dm=trim(filapp)//'_1DM'

!--------------------------------------------------------
!Names based on dtfil%filnam_ds(5)+iapp

!Prepare the name of the auxiliary files for protection
 call fappnd(filprot,dtfil%filnam_ds(5),iapp)
 dtfil%fnametmp_app_den=trim(filprot)//'_DEN'
 dtfil%fnametmp_app_kden=trim(filprot)//'_KDEN'

 DBG_EXIT("COLL")

end subroutine dtfil_init_time
!!***
