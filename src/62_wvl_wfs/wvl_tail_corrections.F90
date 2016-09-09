!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_tail_corrections
!! NAME
!! wvl_tail_corrections
!!
!! FUNCTION
!! Perform a minimization on the wavefunctions (especially the treatment
!! of the kinetic operator) with exponentialy decreasing functions on
!! boundaries.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2016 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      afterscfloop
!!
!! CHILDREN
!!      calculatetailcorrection,dcopy,wrtout,xmpi_allgatherv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_tail_corrections(dtset, energies, etotal, mpi_enreg, psps, wvl, xcart)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_xmpi
 use m_errors
 use m_profiling_abi

 use m_energies, only : energies_type
#if defined HAVE_BIGDFT
  use BigDFT_API, only: CalculateTailCorrection
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_tail_corrections'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(out) :: etotal
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(energies_type),intent(inout) :: energies
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_data),intent(inout) :: wvl
!arrays
 real(dp),intent(in) :: xcart(3,dtset%natom)

!Local variables-------------------------------
#if defined HAVE_BIGDFT
!scalars
 integer :: ierr,me,nbuf,nproc,nsize,spaceComm
 real(dp) :: ekin_sum,epot_sum,eproj_sum
 logical :: parallel
 character(len=500) :: message
!arrays
 integer :: ntails(3)
 real(dp) :: atails(3)
#endif

! *************************************************************************

#if defined HAVE_BIGDFT

 spaceComm=mpi_enreg%comm_wvl
 me=xmpi_comm_rank(spaceComm)
 nproc=xmpi_comm_size(spaceComm)
 parallel = (nproc > 1)

!Write a message with the total energy before tail corrections.
 etotal = energies%e_kinetic + energies%e_hartree + energies%e_xc + &
&         energies%e_localpsp + energies%e_corepsp + energies%e_fock+&
&         energies%e_entropy + energies%e_elecfield + energies%e_magfield+&
&         energies%e_ewald + energies%e_vdw_dftd
 if (dtset%usepaw==0) etotal = etotal + energies%e_nonlocalpsp
 if (dtset%usepaw/=0) etotal = etotal + energies%e_paw

 write(message,'(a,2x,e19.12)') ' Total energy before tail correction', etotal
 call wrtout(std_out, message, 'COLL')

!Calculate kinetic energy correction due to boundary conditions
 nbuf = nint(dtset%tl_radius / dtset%wvl_hgrid)
 ntails = (/ wvl%descr%Glr%d%n1, wvl%descr%Glr%d%n2, wvl%descr%Glr%d%n3 /) + 2 * nbuf
 atails = real(ntails, dp) * dtset%wvl_hgrid
 write(message,'(a,a,i6,a,A,A,3F12.6,A,A,3I12,A)') ch10,&
& ' Tail requires ',nbuf,' additional grid points around cell.', ch10, &
& '  | new acell:', atails, ch10, &
& '  | new box size for wavelets:', ntails, ch10
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')


!Calculate energy correction due to finite size effects
!---reformat potential
 nsize = wvl%descr%Glr%d%n1i * wvl%descr%Glr%d%n2i
 ABI_ALLOCATE(wvl%den%denspot%pot_work, (nsize * wvl%descr%Glr%d%n3i * dtset%nsppol))

 if (parallel) then
   call xmpi_allgatherv(wvl%den%denspot%rhov, &
&   nsize * wvl%den%denspot%dpbox%nscatterarr(me, 2), &
&   wvl%den%denspot%pot_work, nsize * wvl%den%denspot%dpbox%ngatherarr(:,2), &
&   nsize * wvl%den%denspot%dpbox%ngatherarr(:,3),spaceComm,ierr)
 else
   call dcopy(wvl%descr%Glr%d%n1i * wvl%descr%Glr%d%n2i * &
&   wvl%descr%Glr%d%n3i * dtset%nsppol,wvl%den%denspot%rhov,1,wvl%den%denspot%pot_work,1)
 end if

 if(dtset%usepaw==1) then
   call CalculateTailCorrection(me, nproc, wvl%descr%atoms, dtset%tl_radius, &
&   wvl%wfs%ks%orbs, wvl%wfs%ks%lzd%Glr, wvl%projectors%nlpsp, dtset%tl_nprccg, &
&   wvl%den%denspot%pot_work, dtset%wvl_hgrid, xcart, psps%gth_params%radii_cf, &
&   dtset%wvl_crmult, dtset%wvl_frmult, dtset%nsppol, &
&   wvl%wfs%ks%psi, .false., ekin_sum, epot_sum, eproj_sum,&
&   wvl%projectors%G,wvl%descr%paw)
 else
   call CalculateTailCorrection(me, nproc, wvl%descr%atoms, dtset%tl_radius, &
&   wvl%wfs%ks%orbs, wvl%wfs%ks%lzd%Glr, wvl%projectors%nlpsp, dtset%tl_nprccg, &
&   wvl%den%denspot%pot_work, dtset%wvl_hgrid, xcart, psps%gth_params%radii_cf, &
&   dtset%wvl_crmult, dtset%wvl_frmult, dtset%nsppol, &
&   wvl%wfs%ks%psi, .false., ekin_sum, epot_sum, eproj_sum)
 end if

 ABI_DEALLOCATE(wvl%den%denspot%pot_work)

 energies%e_kinetic = ekin_sum
 energies%e_localpsp = epot_sum - two * energies%e_hartree
 energies%e_nonlocalpsp = eproj_sum
 energies%e_corepsp = zero
#if defined HAVE_BIGDFT
 energies%e_localpsp = energies%e_localpsp - wvl%e%energs%evxc
#endif

 write(message,'(a,3(1x,e18.11))') ' ekin_sum,epot_sum,eproj_sum',  &
 ekin_sum,epot_sum,eproj_sum
 call wrtout(std_out, message, 'COLL')
 write(message,'(a,2(1x,e18.11))') ' ehart,eexcu', &
& energies%e_hartree,energies%e_xc
 call wrtout(std_out, message, 'COLL')

 etotal = energies%e_kinetic + energies%e_hartree + energies%e_xc + &
&         energies%e_localpsp + energies%e_corepsp + energies%e_fock+&
&         energies%e_entropy + energies%e_elecfield + energies%e_magfield+&
&         energies%e_ewald + energies%e_vdw_dftd
 if (dtset%usepaw==0) etotal = etotal + energies%e_nonlocalpsp
 if (dtset%usepaw/=0) etotal = etotal + energies%e_paw

 write(message,'(a,2x,e19.12)') ' Total energy with tail correction', etotal
 call wrtout(std_out, message, 'COLL')

!--- End if of tail calculation

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) etotal,mpi_enreg%nproc,dtset%nstep,energies%e_ewald,psps%npsp,&
& wvl%wfs%ks,xcart(1,1)
#endif

end subroutine wvl_tail_corrections
!!***
