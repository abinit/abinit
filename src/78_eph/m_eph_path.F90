!!****m* ABINIT/m_eph_path
!! NAME
!!  m_eph_path
!!
!! FUNCTION
!!  Compute e-ph matrix elements along a path
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2024 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_eph_path

 use, intrinsic :: iso_c_binding
 use defs_basis
 use m_abicore
 use m_xmpi
 use m_mpinfo
 use m_errors
 use m_ifc
 use m_ddb
 use m_dvdb
! use m_fft
 use m_hamiltonian
! use m_pawcprj
! use m_wfd
! use m_hdr
! use m_sigtk
! use m_ephtk
! use netcdf
! use m_nctk
 use m_dtset
 use m_dtfil
! use m_clib
! use m_mkffnl
!
 use defs_abitypes,    only : mpi_type
 use defs_datatypes,   only : ebands_t, pseudopotential_type
 use m_time,           only : cwtime, cwtime_report, timab, sec2str
! use m_fstrings,       only : itoa, ftoa, sjoin, ktoa, ltoa, strcat
! use m_numeric_tools,  only : arth, c2r, get_diag, linfit, iseven, simpson_cplx, simpson, print_arr, inrange
! use m_io_tools,       only : iomode_from_fname, file_exists, is_open, open_file, flush_unit
! use m_fftcore,        only : ngfft_seq, sphereboundary, get_kg, kgindex
! use m_cgtools,        only : cg_zdotc, cg_real_zdotc, cg_zgemm, fxphas_seq
 use m_crystal,        only : crystal_t
! use m_getgh1c,        only : getgh1c, rf_transgrid_and_pack, getgh1c_setup
 use m_pawang,         only : pawang_type
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type
! use m_pawrhoij,       only : pawrhoij_type
 use m_pawfgr,         only : pawfgr_type
! use m_phonons,        only : phstore_t, phstore_new
 use m_pstat,          only : pstat_t
 use m_cgwf,           only : nscf_t

 implicit none

 private
!!***

 public :: eph_path_run


contains  !=====================================================
!!***

!----------------------------------------------------------------------

!!****f* m_eph_path/eph_path_run
!! NAME
!!  eph_path_run
!!
!! FUNCTION
!!  Compute e-ph matrix elements along a path.
!!
!! INPUTS
!! ngfft(18),ngfftf(18)=Coarse and Fine FFT meshes.
!! dtset<dataset_type>=All input variables for this dataset.
!! dvdb<dbdb_type>=Database with the DFPT SCF potentials.
!! ifc<ifc_type>=interatomic force constants and corresponding real space grid info.
!! pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!! pawang<pawang_type)>=PAW angular mesh and related data.
!! pawrad(ntypat*usepaw)<pawrad_type>=Paw radial mesh and related data.
!! pawtab(ntypat*usepaw)<pawtab_type>=Paw tabulated starting data.
!! psps<pseudopotential_type>=Variables related to pseudopotentials.
!! comm=MPI communicator.
!!
!! OUTPUT
!!
!! SOURCE

subroutine eph_path_run(dtfil, ngfft, ngfftf, dtset, cryst, dvdb, ifc, &
                        pawfgr, pawang, pawrad, pawtab, psps, mpi_enreg, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(dvdb_t),intent(inout) :: dvdb
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
 type(pawfgr_type),intent(in) :: pawfgr
 type(ifc_type),intent(in) :: ifc
 type(mpi_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18),ngfftf(18)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables ------------------------------
!scalars
  integer,parameter :: istwfk_1 = 1 ! ,tim_getgh1c1 = 1, berryopt0 = 0, ider0 = 0, idir0 = 0,
! integer,parameter :: useylmgr = 0, useylmgr1 =0, master = 0, ndat1 = 1, cplex1 = 1, pawread0 = 0
  integer :: spin, iqpt, nqpt, nband, ierr
  real(dp) :: cpu_all,wall_all,gflops_all
! complex(dp) :: cfact,dka,dkap,dkpa,dkpap, cnum, sig_cplx, cfact2
 type(gs_hamiltonian_type) :: gs_hamk, gs_hamkq
 type(nscf_t) :: nscf
 character(len=5000) :: msg
!!arrays
 integer :: units(2)
 real(dp) :: kpt(3), qpt(3), kq(3)
! !real(dp),allocatable :: phfreqs_qibz(:,:), pheigvec_qibz(:,:,:,:), eigvec_qpt(:,:,:)
! real(dp) :: ylmgr_dum(1,1,1)
!************************************************************************

 if (psps%usepaw == 1) then
   ABI_ERROR("PAW not implemented")
   ABI_UNUSED((/pawang%nsym, pawrad(1)%mesh_size/))
 end if

 !my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 call cwtime(cpu_all, wall_all, gflops_all, "start")

 units = [std_out, ab_out]

 !call pstat%from_pid(); call pstat%print([std_out], header="Memory at the beginning of eph_path")

 !! Copy important dimensions
 !natom = cryst%natom; natom3 = 3 * natom; nsppol = ebands%nsppol; nspinor = ebands%nspinor
 !nspden = dtset%nspden; nkpt = ebands%nkpt

 call nscf%init(dtset, dtfil, cryst, comm)
 nband = dtset%mband

 do spin=1,dtset%nsppol
   kpt = zero
   call nscf%solve(spin, kpt, istwfk_1, nband, cryst, dtset, dtfil, psps, pawtab, pawfgr, gs_hamk, msg, ierr)
   ABI_CHECK(ierr == 0, msg)
   nqpt = 1
   do iqpt=1,nqpt
     qpt = zero
     kq = kpt + qpt
     call nscf%solve(spin, kq, istwfk_1, nband, cryst, dtset, dtfil, psps, pawtab, pawfgr, gs_hamkq, msg, ierr)
     ABI_CHECK(ierr == 0, msg)
     call gs_hamkq%free()
   end do ! iqpt
   call gs_hamk%free()
 end do ! spin

 call nscf%free()

 call cwtime_report(" eph_path: MPI barrier before returning.", cpu_all, wall_all, gflops_all, end_str=ch10, comm=comm)

end subroutine eph_path_run
!!***

end module m_eph_path
!!***
