!!****m* ABINIT/m_gwr_driver
!! NAME
!!  m_gwr_driver
!!
!! FUNCTION
!!   Driver for GWR calculations
!!
!! COPYRIGHT
!!  Copyright (C) 2021-2022 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_gwr_driver

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_xomp
 use m_hdr
 use m_crystal
 use m_ebands
 use m_dtset
 use m_dtfil
 use m_wfk
 use m_distribfft
!#ifdef HAVE_NETCDF
! use netcdf
!#endif
 use m_nctk

 use defs_datatypes,    only : pseudopotential_type, ebands_t
 use defs_abitypes,     only : MPI_type
 use m_io_tools,        only : file_exists, open_file
 use m_time,            only : cwtime, cwtime_report
 use m_fstrings,        only : strcat, sjoin, ftoa, itoa
 use m_fftcore,         only : print_ngfft
 use m_mpinfo,          only : destroy_mpi_enreg, initmpi_seq
 use m_pawang,          only : pawang_type
 use m_pawrad,          only : pawrad_type
 use m_pawtab,          only : pawtab_type
 use m_paw_an,          only : paw_an_type, paw_an_free !, paw_an_nullify, paw_an_init,
 use m_paw_ij,          only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify
 use m_pawfgrtab,       only : pawfgrtab_type, pawfgrtab_free, pawfgrtab_init
 use m_pawrhoij,        only : pawrhoij_type, pawrhoij_alloc, pawrhoij_copy, pawrhoij_free, pawrhoij_symrhoij
 use m_pawfgr,          only : pawfgr_type, pawfgr_init, pawfgr_destroy
 use m_pspini,          only : pspini
 use m_gwr,             only : gwr_new, gwr_t
 !use m_ephtk,          only : ephtk_update_ebands

 implicit none

 private
!!***

 public :: gwr_driver
!!***

contains
!!***

!!****f* m_gwr_driver/gwr_driver
!! NAME
!!  gwr_driver
!!
!! FUNCTION
!! Main routine for GWR calculations.
!!
!! INPUTS
!! acell(3)=Length scales of primitive translations (bohr)
!! codvsn=Code version
!! dtfil<datafiles_type>=Variables related to files.
!! dtset<dataset_type>=All input variables for this dataset.
!! pawang<pawang_type)>=PAW angular mesh and related data.
!! pawrad(ntypat*usepaw)<pawrad_type>=Paw radial mesh and related data.
!! pawtab(ntypat*usepaw)<pawtab_type>=Paw tabulated starting data.
!! psps<pseudopotential_type>=Variables related to pseudopotentials.
!!   Before entering the first time in the routine, a significant part of Psps has been initialized :
!!   the integers dimekb,lmnmax,lnmax,mpssang,mpssoang,mpsso,mgrid,ntypat,n1xccc,usepaw,useylm,
!!   and the arrays dimensioned to npsp. All the remaining components of Psps are to be initialized in
!!   the call to pspini. The next time the code enters bethe_salpeter, Psps might be identical to the
!!   one of the previous Dtset, in which case, no reinitialisation is scheduled in pspini.F90.
!! rprim(3,3)=Dimensionless real space primitive translations.
!! xred(3,natom)=Reduced atomic coordinates.
!!
!! PARENTS
!!      m_driver
!!
!! NOTES
!!
!! ON THE USE OF FFT GRIDS:
!! =================
!! In case of PAW:
!! ---------------
!!    Two FFT grids are used:
!!    - A "coarse" FFT grid (defined by ecut) for the application of the Hamiltonian on the plane waves basis.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      Hamiltonian, wave-functions, density related to WFs (rhor here), ... are expressed on this grid.
!!    - A "fine" FFT grid (defined) by ecutdg) for the computation of the density inside PAW spheres.
!!      It is defined by nfftf, ngfftf, mgfftf, ... Total density, potentials, ... are expressed on this grid.
!! In case of norm-conserving:
!! ---------------------------
!!    - Only the usual FFT grid (defined by ecut) is used. It is defined by nfft, ngfft, mgfft, ...
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf) are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! CHILDREN
!!
!! SOURCE

subroutine gwr_driver(acell, codvsn, dtfil, dtset, pawang, pawrad, pawtab, psps, rprim, xred)

!Arguments ------------------------------------
!scalars
 character(len=8),intent(in) :: codvsn
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 real(dp),intent(in) :: acell(3),rprim(3,3),xred(3,dtset%natom)
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: ii, comm, nprocs, my_rank, psp_gencond, mgfftf, nfftf
 integer :: omp_ncpus, work_size, nks_per_proc, ierr
 real(dp):: eff, mempercpu_mb, max_wfsmem_mb, nonscal_mem
!#ifdef HAVE_NETCDF
! integer :: ncid, ncerr
!#endif
 real(dp) :: ecore, ecut_eff, ecutdg_eff, gsqcutc_eff, gsqcutf_eff
 real(dp) :: cpu, wall, gflops
 logical :: use_wfk
 character(len=500) :: msg
 character(len=fnlen) :: wfk0_path, path
 type(hdr_type) :: wfk0_hdr
 type(crystal_t) :: cryst
 type(ebands_t) :: ebands
 type(pawfgr_type) :: pawfgr
 type(mpi_type) :: mpi_enreg
 type(gwr_t) :: gwr
!arrays
 integer :: ngfftc(18), ngfftf(18)
 real(dp),parameter :: k0(3) = zero
 !type(pawfgrtab_type),allocatable :: pawfgrtab(:)
 !type(paw_ij_type),allocatable :: paw_ij(:)
 !type(paw_an_type),allocatable :: paw_an(:)

!************************************************************************

 ! This part performs the initialization of the basic objects used to perform e-ph calculations:
 !
 !     1) Crystal structure `cryst`
 !     2) Ground state band energies: `ebands`
 !     5) Pseudos and PAW basic objects.
 !
 ! Once we have these objects, we can call specialized routines for e-ph calculations.
 ! Notes:
 !
 !   * Any modification to the basic objects mentioned above should be done here (e.g. change of efermi)
 !   * This routines shall not allocate big chunks of memory. The CPU-demanding sections should be
 !     performed in the subdriver that will employ different MPI distribution schemes optimized for that particular task.

 if (psps%usepaw == 1) then
   ABI_ERROR("PAW not implemented")
   ABI_UNUSED((/pawang%nsym, pawrad(1)%mesh_size/))
 end if

 ! abirules!
 if (.False.) write(std_out,*)acell,codvsn,rprim,xred

 comm = xmpi_world; nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 ! Initialize filenames
 wfk0_path = dtfil%fnamewffk
 use_wfk = .True.

 if (my_rank == master) then
   ! Accept WFK file in Fortran or netcdf format.
   if (use_wfk .and. nctk_try_fort_or_ncfile(wfk0_path, msg) /= 0) then
     ABI_ERROR(sjoin("Cannot find GS WFK file:", wfk0_path, ". Error:", msg))
   end if
 end if ! master

 ! Broadcast filenames (needed because they might have been changed if we are using netcdf files)
 if (use_wfk) then
   call xmpi_bcast(wfk0_path, master, comm, ierr)
   call wrtout(ab_out, sjoin("- Reading GS states from WFK file:", wfk0_path))
 end if
 call wrtout(ab_out, ch10//ch10)

 ! autoparal section
 ! TODO: This just to activate autoparal in AbiPy. Lot of things should be improved.
 if (dtset%max_ncpus /= 0) then
   write(ab_out,'(a)')"--- !Autoparal"
   write(ab_out,"(a)")"# Autoparal section for GWR runs"
   write(ab_out,"(a)")   "info:"
   write(ab_out,"(a,i0)")"    autoparal: ",dtset%autoparal
   write(ab_out,"(a,i0)")"    max_ncpus: ",dtset%max_ncpus
   write(ab_out,"(a,i0)")"    nkpt: ",dtset%nkpt
   write(ab_out,"(a,i0)")"    nsppol: ",dtset%nsppol
   write(ab_out,"(a,i0)")"    nspinor: ",dtset%nspinor
   write(ab_out,"(a,i0)")"    mband: ",dtset%mband
   write(ab_out,"(a,i0)")"    eph_task: ",dtset%eph_task

   work_size = dtset%nkpt * dtset%nsppol
   ! Non-scalable memory in Mb i.e. memory that is not distributed with MPI.
   nonscal_mem = zero
   max_wfsmem_mb = (two * dp * dtset%mpw * dtset%mband * dtset%nkpt * dtset%nsppol * dtset%nspinor * b2Mb) * 1.1_dp

   ! List of configurations.
   ! Assuming an OpenMP implementation with perfect speedup!
   write(ab_out,"(a)")"configurations:"

   do ii=1,dtset%max_ncpus
     nks_per_proc = work_size / ii
     nks_per_proc = nks_per_proc + mod(work_size, ii)
     eff = (one * work_size) / (ii * nks_per_proc)
     ! Add the non-scalable part and increase by 10% to account for other datastructures.
     mempercpu_mb = (max_wfsmem_mb + nonscal_mem) * 1.1_dp

     do omp_ncpus=1,1 !xomp_get_max_threads()
       write(ab_out,"(a,i0)")"    - tot_ncpus: ",ii * omp_ncpus
       write(ab_out,"(a,i0)")"      mpi_ncpus: ",ii
       write(ab_out,"(a,i0)")"      omp_ncpus: ",omp_ncpus
       write(ab_out,"(a,f12.9)")"      efficiency: ",eff
       write(ab_out,"(a,f12.2)")"      mem_per_cpu: ",mempercpu_mb
     end do
   end do
   write(ab_out,'(a)')"..."
   ABI_STOP("Stopping now!")
 end if

 call cwtime(cpu, wall, gflops, "start")
 call cwtime_report(" eph%init", cpu, wall, gflops)

 if (use_wfk) then
   ! Construct crystal and ebands from the GS WFK file.
   ebands = wfk_read_ebands(wfk0_path, comm, out_hdr=wfk0_hdr)
   call wfk0_hdr%vs_dtset(dtset)

   cryst = wfk0_hdr%get_crystal()
   call cryst%print(header="crystal structure from WFK file")

   ! Here we change the GS bands (Fermi level, scissors operator ...)
   ! All the modifications to ebands should be done here.
   !call ephtk_update_ebands(dtset, ebands, "Ground state energies")

   ! TODO: Make sure that ef is inside the gap if semiconductor.
 end if

 call pawfgr_init(pawfgr, dtset, mgfftf, nfftf, ecut_eff, ecutdg_eff, ngfftc, ngfftf, &
                  gsqcutc_eff=gsqcutc_eff, gsqcutf_eff=gsqcutf_eff, gmet=cryst%gmet, k0=k0)

 call print_ngfft(ngfftc, header='Coarse FFT mesh used for the wavefunctions')
 call print_ngfft(ngfftf, header='Dense FFT mesh used for densities and potentials')

 ! Fake MPI_type for the sequential part.
 call initmpi_seq(mpi_enreg)
 call init_distribfft_seq(mpi_enreg%distribfft, 'c', ngfftc(2), ngfftc(3), 'all')
 call init_distribfft_seq(mpi_enreg%distribfft, 'f', ngfftf(2), ngfftf(3), 'all')

 ! ===========================================
 ! === Open and read pseudopotential files ===
 ! ===========================================
 call pspini(dtset, dtfil, ecore, psp_gencond, gsqcutc_eff, gsqcutf_eff, pawrad, pawtab, psps, cryst%rprimd, comm_mpi=comm)

 gwr = gwr_new(dtset, cryst, psps, pawtab, ebands, mpi_enreg, comm)

 if (use_wfk) then
   call gwr%build_gtau_from_wfk(wfk0_path)
 else
 end if

 ! ====================================================
 ! === This is the real GWR stuff once all is ready ===
 ! ====================================================

 !select case (dtset%gwr_task)
 !case ("RPA_ENERGY")
 !  call gwr%rpa_energy()
 !case ("G0W0")
 !  call gwr%g0w0()
 !case default
 !  !ABI_ERROR(sjoin("Invalid value of gwr_task:", itoa(dtset%gwr_task)))
 !end select

 !=====================
 !==== Free memory ====
 !=====================
 call gwr%free()
 call cryst%free()
 call wfk0_hdr%free()
 call ebands_free(ebands)
 call pawfgr_destroy(pawfgr)
 call destroy_mpi_enreg(mpi_enreg)

 ! Deallocation for PAW.
 if (dtset%usepaw == 1) then
   !call pawrhoij_free(pawrhoij)
   !ABI_FREE(pawrhoij)
   !call pawfgrtab_free(pawfgrtab)
   !ABI_FREE(pawfgrtab)
   !call paw_ij_free(paw_ij)
   !ABI_FREE(paw_ij)
   !call paw_an_free(paw_an)
   !ABI_FREE(paw_an)
 end if

end subroutine gwr_driver
!!***

end module m_gwr_driver
!!***
