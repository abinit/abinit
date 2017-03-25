!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_mkrho
!! NAME
!! wvl_mkrho
!!
!! FUNCTION
!! This method is just a wrapper around the BigDFT routine to compute the
!! density from the wavefunctions.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=input variables.
!!  mpi_enreg=informations about MPI parallelization
!!  occ(dtset%mband)=occupation numbers.
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  wvl_wfs <type(wvl_projector_type)>=wavefunctions informations for wavelets.
!!
!! OUTPUT
!!  rhor(dtset%nfft)=electron density in r space
!!
!! SIDE EFFECTS
!!  proj <type(wvl_projector_type)>=projectors informations for wavelets.
!!   | proj(OUT)=computed projectors.
!!
!! PARENTS
!!      afterscfloop,gstate,mkrho,mover,vtorho
!!
!! CHILDREN
!!      communicate_density,sumrho,wrtout,wvl_rho_abi2big
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine wvl_mkrho(dtset, irrzon, mpi_enreg, phnons, rhor, wvl_wfs, wvl_den)

 use defs_basis
 use defs_abitypes
 use defs_wvltypes
 use m_profiling_abi
 use m_errors
 use m_abi2big
 use m_xmpi
#if defined HAVE_BIGDFT
  use BigDFT_API, only : sumrho, symmetry_data, ELECTRONIC_DENSITY, communicate_density
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_mkrho'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(wvl_wf_type),intent(inout) :: wvl_wfs
 type(wvl_denspot_type), intent(inout) :: wvl_den
!arrays
 real(dp),intent(inout) :: rhor(dtset%nfft,dtset%nspden)
 integer, target, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,  &
&               (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 real(dp), target, intent(in) :: phnons(2,(dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3))**(1-1/dtset%nsym),  &
&                                 (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))

!Local variables-------------------------------
#if defined HAVE_BIGDFT
!scalars
 character(len=500) :: message
 integer :: comm,me,nproc
 type(symmetry_data) :: sym
 !for debugging:
 !integer::ifile,ierr
#endif

! *************************************************************************

 DBG_ENTER("COLL")

#if defined HAVE_BIGDFT
 comm=mpi_enreg%comm_wvl
 me=xmpi_comm_rank(comm)
 nproc=xmpi_comm_size(comm)

 sym%symObj = wvl_den%symObj
 sym%irrzon => irrzon
 sym%phnons => phnons

 call sumrho(wvl_den%denspot%dpbox,wvl_wfs%ks%orbs,wvl_wfs%ks%Lzd,&
& wvl_wfs%GPU,sym,wvl_den%denspot%rhod,wvl_den%denspot%xc,&
& wvl_wfs%ks%psi,wvl_den%denspot%rho_psi)

 call communicate_density(wvl_den%denspot%dpbox,wvl_wfs%ks%orbs%nspin,&
& wvl_den%denspot%rhod,wvl_den%denspot%rho_psi,wvl_den%denspot%rhov,.false.)

 wvl_den%denspot%rhov_is = ELECTRONIC_DENSITY
 write(message, '(a,a,a,a)' ) ch10, ' wvl_mkrho : but why are you copying me :..o('
 call wrtout(std_out,message,'COLL')

 call wvl_rho_abi2big(2,rhor,wvl_den)

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) mpi_enreg%me,dtset%nstep,wvl_wfs%ks,wvl_den%symObj,&
& rhor(1,1),irrzon(1,1,1),phnons(1,1,1)
#endif

 DBG_EXIT("COLL")

end subroutine wvl_mkrho
!!***
