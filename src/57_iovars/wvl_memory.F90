!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_memory
!! NAME
!! wvl_memory
!!
!! FUNCTION
!! Estimation of the memory needed for waelet based computation job.
!! According to the value of the option variable,
!! might also try to allocate this amount of memory, and if it fails,
!! might estimate the available memory.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset=<type datafiles_type>contains all input variables.
!!  idtset=number of the current dataset
!!  mpi_enreg=informations about MPI parallelization
!!  npsp=number of pseudopotentials
!!  option : if 0 , no test of available memory
!!           if 1 , the routine tries to allocate the estimated memory, for testing
!!                    purposes, and if a failure occurs, the routine stops.
!!           if 2 , like 1, but before stopping, the routine will provide
!!                    an estimation of the available memory.
!!  pspheads(npsp)=<type pspheader_type>all the important information from the
!!   pseudopotential file header, as well as the psp file name
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!! The estimator is the one provided by BigDFT.
!!
!! PARENTS
!!      memory_eval
!!
!! CHILDREN
!!      atomic_info,createwavefunctionsdescriptors,deallocate_lr
!!      memoryestimator,mkradim,wrtout,wvl_descr_atoms_set,wvl_descr_free
!!      wvl_setboxgeometry,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_memory(dtset, idtset, mpi_enreg, npsp, option, pspheads)
  
 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_profiling_abi
 use m_errors
 use m_xmpi

#if defined HAVE_BIGDFT
 use BigDFT_API, only: MemoryEstimator, createWavefunctionsDescriptors, deallocate_lr, &
      & atomic_info, memory_estimation
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_memory'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
 use interfaces_43_wvl_wrappers
!End of the abilint section

  implicit none

!Arguments ------------------------------------
  !scalars
  integer,intent(in) :: idtset, npsp, option
  type(dataset_type),intent(in) :: dtset
  type(MPI_type),intent(in) :: mpi_enreg
  !arrays
  type(pspheader_type),intent(in) :: pspheads(npsp)

!Local variables-------------------------------
#if defined HAVE_BIGDFT
  !scalars
  integer :: ityp, i, mu, nstates, me, nproc, comm
  character(len=500) :: message
  real(dp) :: ehomo, radfine
  type(wvl_internal_type) :: wvl
  type(memory_estimation) :: peakmem
  !arrays
  real(dp) :: acell(3), rprimd(3,3), rprim(3,3)
  real(dp), allocatable :: radii_cf(:,:)
  real(dp), allocatable :: xred(:,:), xcart(:,:)
#endif

! **************************************************************************
 
#if defined HAVE_BIGDFT

 comm=mpi_enreg%comm_wvl
 me=xmpi_comm_rank(comm)
 nproc=xmpi_comm_size(comm)

 if(option<0 .or. option>2)then
   write(message, '(A,A,A,A,I0,A)') ch10,&
&   ' wvl_memory : BUG -',ch10,&
&   '  option=',option,' while the only allowed values are 0, 1, or 2.'
   call wrtout(std_out,message,'COLL')
 end if

 wvl%paw%usepaw=0 !no PAW here
 nullify(wvl%rholoc%d)
 nullify(wvl%rholoc%msz)
 nullify(wvl%rholoc%rad)
 nullify(wvl%rholoc%radius)
 nullify(wvl%paw%spsi)
 nullify(wvl%paw%indlmn)
 nullify(wvl%paw%spsi)
 nullify(wvl%paw%indlmn)

 write(message,*)' wvl_memory : analysis of memory needs '
 call wrtout(std_out,message,'COLL')

 if(idtset>=100)then
   write(message,'(80a,a,a,i5,a)')('=',mu=1,80),ch10,&
&   ' Values of the parameters that define the memory need for DATASET', idtset,&
&   ' (WVL).'
 else if(idtset/=0)then
   write(message,'(80a,a,a,i3,a)')('=',mu=1,80),ch10,&
&   ' Values of the parameters that define the memory need for DATASET', idtset,&
&   ' (WVL).'
 else
   write(message,'(80a,a,a,a)')('=',mu=1,80),ch10,&
&   ' Values of the parameters that define the memory need of the present run',&
&   ' (WVL).'
 end if
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 write(message,'( a,f7.3,a,i7,2(a,F7.3),a,a,f7.3,a,i7 )' ) &
& '  wvl_hgrid =', dtset%wvl_hgrid , '   nwfshist =', dtset%nwfshist, &
& ' wvl_crmult =', dtset%wvl_crmult, ' wvl_frmult =', dtset%wvl_frmult, ch10,&
& '  tl_radius =', dtset%tl_radius , '  tl_nprccg =', dtset%tl_nprccg
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 if (dtset%nsppol == 2) then
   nstates = dtset%nelect
 else
   nstates = dtset%mband
 end if
 write(message,'(4(a,i7))')&
& '      natom =', dtset%natom, '     ntypat =', dtset%ntypat, &
& '    nstates =', nstates,     '     nsppol =', dtset%nsppol
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 write(message,'(80a)') ('=',mu=1,80)
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

!First, use eleconf to get radii_cf().
 ABI_ALLOCATE(radii_cf,(npsp, 3))
 do ityp = 1, npsp, 1
   call atomic_info(int(pspheads(ityp)%znuclpsp), int(pspheads(ityp)%zionpsp), ehomo = ehomo)

!  new method for assigning the radii
   radii_cf(ityp, 1) = one / sqrt(abs(two * ehomo))
   radfine = 100.d0
   do i = 0, 4, 1
     if (pspheads(ityp)%GTHradii(i) /= zero) then
       radfine = min(radfine, pspheads(ityp)%GTHradii(i))
     end if
   end do
   radii_cf(ityp,2) = radfine
 end do

!Compute the shifted positions and acell
 acell = dtset%acell_orig(1:3,1)
 call wvl_descr_atoms_set(acell, dtset%icoulomb, dtset%natom, dtset%ntypat, dtset%typat, wvl)
 ABI_ALLOCATE(xred,(3, dtset%natom))
 xred = dtset%xred_orig(:,:,1)
 rprimd = dtset%rprimd_orig(1:3,1:3,1)
 wvl%h(:) = dtset%wvl_hgrid
 call wvl_setBoxGeometry(1, radii_cf, rprimd, xred, &
& wvl, dtset%wvl_crmult, dtset%wvl_frmult)
!Compute acell and rprim from rprimd
 call mkradim(acell,rprim,rprimd)
 ABI_ALLOCATE(xcart,(3, dtset%natom))
 call xred2xcart(dtset%natom, rprimd, xcart, xred)
 call createWavefunctionsDescriptors(me, wvl%h(1), wvl%h(2), wvl%h(3), &
& wvl%atoms, xcart, radii_cf, dtset%wvl_crmult, dtset%wvl_frmult, wvl%Glr)
 call MemoryEstimator(nproc, dtset%nwfshist, wvl%Glr, &
& dtset%mband, dtset%nspinor, dtset%nkpt, 0, dtset%nsppol, &
& 0, dtset%iscf, peakmem)

 call deallocate_lr(wvl%Glr)
 call wvl_descr_free(wvl)
 ABI_DEALLOCATE(radii_cf)
 ABI_DEALLOCATE(xred)
 ABI_DEALLOCATE(xcart)

 write(message,'(80a,a)') ('=',mu=1,80), ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) idtset,npsp,option,dtset%nstep,mpi_enreg%nproc,pspheads(1)%zionpsp
#endif

end subroutine wvl_memory
!!***
