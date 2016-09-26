!{\src2tex{textfont=tt}}
!!****f* ABINIT/indefo1
!! NAME
!! indefo1
!!
!! FUNCTION
!! Initialisation phase : defaults values for a first batch of input variables
!! (especially dimensions, needed to allocate other parts of dtsets, as well
!!  as other input variables whose existence is needed for other initialisations to proceed).
!!
!! COPYRIGHT
!! Copyright (C) 2009-2016 ABINIT group (XG,MM,FF)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!  dtset=<type datafiles_type>contains all input variables for one dataset,
!!   some of which are given a default value here.
!!
!! PARENTS
!!      invars1m
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine indefo1(dtset)

 use defs_basis
 use m_profiling_abi
 use m_errors

 use defs_abitypes,  only : dataset_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'indefo1'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 type(dataset_type),intent(inout) :: dtset !vz_i

!Local variables -------------------------------
!scalars

!******************************************************************
!
!Set up default values. All variables to be output in outvars.f
!should have a default, even if a nonsensible one can be
!chosen to garantee print in that routine.

 DBG_ENTER("COLL")

!Use alphabetic order

!A
 dtset%acell_orig(:,:)=zero
 dtset%algalch(:)=1
 dtset%amu_orig(:,:)=-one
 dtset%autoparal=0
!B
 dtset%bandpp=1
 dtset%berryopt=0
 dtset%berrysav=0
 dtset%bfield(:)=zero
!C
 dtset%cd_customnimfrqs=0
!D
 dtset%densty(:,:)=zero
 dtset%dfield(:)=zero    !!HONG
 dtset%dynimage(:)=1
!E
 dtset%efield(:)=zero
 dtset%efmas_calc_dirs=0
 dtset%efmas_n_dirs=0
!F
!G
 dtset%ga_n_rules=1
 dtset%gw_customnfreqsp=0
 dtset%gw_nqlwl=0
 dtset%gwls_n_proj_freq=0
!I
 dtset%iatfix(:,:)=0
 dtset%icoulomb=0
 dtset%imgmov=0
!J
 dtset%jellslab=0
 dtset%jfielddir(:)=0
!K
 dtset%kptopt=0
!L
 dtset%lexexch(:)=-1
 dtset%lpawu(:)=-1
!M
 dtset%maxestep=0.005d0
 dtset%mixalch_orig(:,:,:)=zero
 dtset%mkmem=-1
 dtset%mkqmem=-1
 dtset%mk1mem=-1
!N
 dtset%nnos=0
 dtset%natpawu=0
 dtset%natsph=0
 dtset%natsph_extra=0
 dtset%natvshift=0
 dtset%nconeq=0
 dtset%ndynimage=1
 dtset%nkpt=-1
 dtset%nkptgw=0
 dtset%nkpthf=0
 dtset%npband=1
 dtset%npfft=1
 dtset%nphf=1
 dtset%npimage=1
 dtset%npkpt=1
 dtset%nppert=1
 dtset%npspalch=0
 dtset%npspinor=1
 dtset%nqptdm=0
 dtset%nspden=1
 dtset%nspinor=1
 dtset%nsppol=1
 dtset%nsym=0     ! Actually, this default value is not used : it is to be reimposed before each call to ingeo in invars1
 dtset%ntimimage=1
 dtset%ntypalch=0
 dtset%ntyppure=-1
 dtset%nucdipmom(:,:)=zero
 dtset%np_slk=1000000
!O
 dtset%optdriver=0
!P
 dtset%paral_rf=0
 dtset%pawspnorb=0  ! will be changed to 1 as soon as usepaw==1 and nspinor==2
 dtset%pimass(:)=-one
!Q
 dtset%qptn=zero
!R
 dtset%red_efield(:)=zero
 dtset%red_dfield(:)=zero
 dtset%red_efieldbar(:)=zero
 dtset%rprim_orig(:,:,:)=zero
 dtset%rprim_orig(1,1,:)=one
 dtset%rprim_orig(2,2,:)=one
 dtset%rprim_orig(3,3,:)=one
!S
 dtset%slabzbeg=zero
 dtset%slabzend=zero
 dtset%so_psp(:)=1
 dtset%spinat(:,:)=zero
 dtset%symmorphi=1
!T
 dtset%tfkinfunc=0
 dtset%typat(:)=0
!U
 dtset%usedmatpu=0
 dtset%usedmft=0
 dtset%useexexch=0
 dtset%usepawu=0
 dtset%usepotzero=0
 dtset%use_slk=0
!V
 dtset%vel_orig(:,:,:)=zero
 dtset%vel_cell_orig(:,:,:)=zero
!W
 dtset%wtq=0
 if (dtset%usepaw==0) dtset%wfoptalg=0
 if (dtset%usepaw/=0) dtset%wfoptalg=10
 if (dtset%optdriver==RUNL_GSTATE.and.dtset%paral_kgb>0) dtset%wfoptalg=14

!X
 dtset%xred_orig(:,:,:)=zero
!Y
!Z
 dtset%zeemanfield(:)=zero

 DBG_EXIT("COLL")

end subroutine indefo1
!!***
