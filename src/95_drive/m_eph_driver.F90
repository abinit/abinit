!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_eph_driver
!! NAME
!!  m_eph_driver
!!
!! FUNCTION
!!   Driver for EPH calculations
!!
!! COPYRIGHT
!!  Copyright (C) 2009-2019 ABINIT group (MG, MVer,GA)
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

module m_eph_driver

 use defs_basis
 use m_errors
 use m_abicore
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_xomp
 use m_hdr
 use m_crystal
 use m_ebands
 use m_dtset
 use m_efmas_defs
 use m_ddk
 use m_ddb
 use m_dvdb
 use m_ifc
 use m_phonons
 use m_nctk
 use m_wfk
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_io_tools,        only : file_exists
 use m_time,            only : cwtime, cwtime_report
 use m_fstrings,        only : strcat, sjoin, ftoa, itoa
 use m_fftcore,         only : print_ngfft
 use m_frohlichmodel,   only : frohlichmodel
 !use m_special_funcs,   only : levi_civita_3
 use m_transport,       only : transport
 use m_mpinfo,          only : destroy_mpi_enreg, initmpi_seq
 use m_pawang,          only : pawang_type
 use m_pawrad,          only : pawrad_type
 use m_pawtab,          only : pawtab_type
 use m_paw_an,          only : paw_an_type, paw_an_free !, paw_an_nullify, paw_an_init,
 use m_paw_ij,          only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify
 use m_pawfgrtab,       only : pawfgrtab_type, pawfgrtab_free, pawfgrtab_init
 use m_pawrhoij,        only : pawrhoij_type, pawrhoij_alloc, pawrhoij_copy, pawrhoij_free, pawrhoij_symrhoij
 use m_pawfgr,          only : pawfgr_type, pawfgr_init, pawfgr_destroy
 use m_phgamma,         only : eph_phgamma
 use m_efmas,           only : efmasdeg_free_array, efmasval_free_array, efmas_ncread
 use m_gkk,             only : eph_gkk, ncwrite_v1qnu
 use m_phpi,            only : eph_phpi
 use m_sigmaph,         only : sigmaph
 use m_pspini,          only : pspini

 implicit none

 private
!!***

 public :: eph
!!***

contains
!!***

!!****f* m_eph_driver/eph
!! NAME
!!  eph
!!
!! FUNCTION
!! Main routine to compute electron phonon coupling matrix elements and
!! calculate related properties - superconductin Tc, phonon linewidths, electronic renormalization
!! due to phonons and temperature effects...
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
!!      driver
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
!!      crystal_free,crystal_from_hdr,crystal_print,cwtime,ddb_free
!!      ddb_from_file,ddk_free,ddk_init,destroy_mpi_enreg,dvdb_free,dvdb_init
!!      dvdb_interpolate_and_write,dvdb_list_perts,dvdb_print,ebands_free
!!      ebands_print,ebands_prtbltztrp,ebands_set_fermie,ebands_set_scheme
!!      ebands_update_occ,ebands_write,edos_free,edos_print,edos_write,eph_gkk
!!      eph_phgamma,eph_phpi,hdr_free,hdr_vs_dtset,ifc_free,ifc_init,ifc_mkphbs
!!      ifc_outphbtrap,ifc_print,ifc_printbxsf
!!      init_distribfft_seq,initmpi_seq,mkphdos,pawfgr_destroy,pawfgr_init
!!      phdos_free,phdos_ncwrite,phdos_print,phdos_print_thermo,print_ngfft
!!      pspini,sigmaph,wfk_read_eigenvalues,wrtout,xmpi_bcast,xmpi_end
!!
!! SOURCE

subroutine eph(acell,codvsn,dtfil,dtset,pawang,pawrad,pawtab,psps,rprim,xred)

!Arguments ------------------------------------
!scalars
 character(len=6),intent(in) :: codvsn
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(pawang_type),intent(inout) :: pawang
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 real(dp),intent(in) :: acell(3),rprim(3,3),xred(3,dtset%natom)
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: master = 0, natifc0 = 0, timrev2 = 2, selectz0 = 0, nsphere0 = 0, prtsrlr0 = 0
 integer :: ii,comm,nprocs,my_rank,psp_gencond,mgfftf,nfftf !,jj
 integer :: iblock_dielt_zeff, iblock_dielt, ddb_nqshift,ierr,brav1
 integer :: omp_ncpus, work_size, nks_per_proc, comm_rpt, method
 real(dp):: eff,mempercpu_mb,max_wfsmem_mb,nonscal_mem
#ifdef HAVE_NETCDF
 integer :: ncid,ncerr
#endif
 real(dp),parameter :: rifcsph0=zero
 real(dp) :: ecore,ecut_eff,ecutdg_eff,gsqcutc_eff,gsqcutf_eff
 real(dp) :: cpu,wall,gflops
 logical :: use_wfk,use_wfq,use_dvdb
 character(len=500) :: msg
 character(len=fnlen) :: wfk0_path, wfq_path, ddb_path, dvdb_path, path
 type(hdr_type) :: wfk0_hdr, wfq_hdr
 type(crystal_t) :: cryst,cryst_ddb
 type(ebands_t) :: ebands, ebands_kq
 type(ddb_type) :: ddb
 type(dvdb_t) :: dvdb
 type(ddk_t) :: ddk
 type(ifc_type) :: ifc, ifc_coarse
 type(pawfgr_type) :: pawfgr
 type(mpi_type) :: mpi_enreg
 type(phonon_dos_type) :: phdos
!arrays
 integer :: ngfftc(18), ngfftf(18), ngqpt_coarse(3), count_wminmax(2)
 integer,allocatable :: dummy_atifc(:)
 real(dp),parameter :: k0(3)=zero
 real(dp) :: wminmax(2), dielt(3,3), zeff(3,3,dtset%natom), zeff_raw(3,3,dtset%natom)
 real(dp),pointer :: gs_eigen(:,:,:)
 real(dp),allocatable :: ddb_qshifts(:,:), kpt_efmas(:,:)
 character(len=fnlen) :: ddk_path(3)
 type(efmasdeg_type),allocatable :: efmasdeg(:)
 type(efmasval_type),allocatable :: efmasval(:,:)
 !type(pawfgrtab_type),allocatable :: pawfgrtab(:)
 !type(paw_ij_type),allocatable :: paw_ij(:)
 !type(paw_an_type),allocatable :: paw_an(:)

!************************************************************************

 ! This part performs the initialization of basic objects used to perform e-ph calculations:
 !
 !     1) Crystal structure `cryst`
 !     2) Ground state band energies: `ebands`
 !     3) Interatomic force constants: `ifc`
 !     4) DVDB database with the dvscf potentials
 !     5) Pseudos and PAW basic objects.
 !
 ! Once we have these objects, we can call specialized routines for e-ph calculations.
 ! Notes:
 !
 !   * Any modification to the basic objects mentioned above should be done here (e.g. change of efermi)
 !   * This routines shall not allocate big chunks of memory. The CPU-demanding sections should be
 !     performed in the subdriver that will employ different MPI distribution schemes optimized for that particular task.

 DBG_ENTER('COLL')

 if (psps%usepaw == 1) then
   MSG_ERROR("PAW not implemented")
   ABI_UNUSED((/pawang%nsym, pawrad(1)%mesh_size/))
 end if

 ! abirules!
 if (.False.) write(std_out,*)acell,codvsn,rprim,xred

 comm = xmpi_world; nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

#ifndef HAVE_MPI_IBCAST
 do ii=1,20
   MSG_WARNING("Your MPI library does not provide MPI_IBCAST. Calculations will be slow")
 end do
#endif

 ! Initialize filenames
 wfk0_path = dtfil%fnamewffk
 wfq_path = dtfil%fnamewffq
 ddb_path = dtfil%filddbsin
 ! Use the ddb file as prefix if getdvdb or irddvb are not given in the input.
 dvdb_path = dtfil%fildvdbin
 if (dvdb_path == ABI_NOFILE) then
   dvdb_path = dtfil%filddbsin; ii=len_trim(dvdb_path); dvdb_path(ii-2:ii+1) = "DVDB"
 end if
 use_wfk = (dtset%eph_task /= 5)
 use_wfq = (dtset%irdwfq /= 0 .or. dtset%getwfq /= 0 .and. dtset%eph_frohlichm /= 1)
 ! If eph_task is needed and ird/get variables are not provided we assume WFQ == WFK
 if (any(dtset%eph_task == [2, -2, 3]) .and. .not. use_wfq) then
   wfq_path = wfk0_path
   use_wfq = .True.
   write(msg, "(4a)")&
       "eph_task requires WFQ but neither irdwfq nor getwfq are specified in the input.", ch10, &
       "Will read WFQ wavefunctions from WFK file:", trim(wfk0_path)
   MSG_COMMENT(msg)
 end if
 use_dvdb = (dtset%eph_task /= 0 .and. dtset%eph_frohlichm /= 1 .and. dtset%eph_task /= 7)

 ddk_path(1) = strcat(dtfil%fnamewffddk, itoa(3*dtset%natom+1))
 ddk_path(2) = strcat(dtfil%fnamewffddk, itoa(3*dtset%natom+2))
 ddk_path(3) = strcat(dtfil%fnamewffddk, itoa(3*dtset%natom+3))

 if (my_rank == master) then
   if (.not. file_exists(ddb_path)) MSG_ERROR(sjoin("Cannot find DDB file:", ddb_path))
   if (use_dvdb .and. .not. file_exists(dvdb_path)) MSG_ERROR(sjoin("Cannot find DVDB file:", dvdb_path))

   ! Accept WFK file in Fortran or netcdf format.
   if (use_wfk .and. nctk_try_fort_or_ncfile(wfk0_path, msg) /= 0) then
     MSG_ERROR(sjoin("Cannot find GS WFK file:", wfk0_path, msg))
   end if
   ! WFQ file
   if (use_wfq) then
     if (nctk_try_fort_or_ncfile(wfq_path, msg) /= 0) then
       MSG_ERROR(sjoin("Cannot find GS WFQ file:", wfq_path, msg))
     end if
   end if

   if (dtset%eph_transport > 0) then
     do ii=1,3
       if (nctk_try_fort_or_ncfile(ddk_path(ii), msg) /= 0) then
         MSG_ERROR(sjoin("Cannot find DDK file:", ddk_path(ii), msg))
       end if
     end do
   end if
 end if ! master

 ! Broadcast filenames (needed because they might have been changed if we are using netcdf files)
 if (use_wfk) then
   call xmpi_bcast(wfk0_path,master,comm,ierr)
   call wrtout(ab_out, sjoin("- Reading GS states from WFK file:", wfk0_path))
 end if
 if (use_wfq) then
   call xmpi_bcast(wfq_path,master,comm,ierr)
   call wrtout(ab_out, sjoin("- Reading GS states from WFQ file:", wfq_path) )
 end if
 call wrtout(ab_out, sjoin("- Reading DDB from file:", ddb_path))
 if (use_dvdb) call wrtout(ab_out, sjoin("- Reading DVDB from file:", dvdb_path))
 if (dtset%eph_transport > 0) then
   call xmpi_bcast(ddk_path,master,comm,ierr)
   call wrtout(ab_out, sjoin("- Reading DDK x from file:", ddk_path(1)))
   call wrtout(ab_out, sjoin("- Reading DDK y from file:", ddk_path(2)))
   call wrtout(ab_out, sjoin("- Reading DDK z from file:", ddk_path(3)))
   ! Read header in DDK files and init basic dimensions.
   ! subdrivers will use ddk to get the matrix elements from file.
   call ddk_init(ddk, ddk_path, comm)
   ! TODO: Should perform consistency check
   !call hdr_vs_dtset(ddk_hdr(ii), dtset)
 end if

 if (dtset%eph_frohlichm /= 1) call wrtout(ab_out, sjoin("- Reading EFMAS information from file:", dtfil%fnameabi_efmas))

 ! autoparal section
 ! TODO: This just to activate autoparal in AbiPy. Lot of things should be improved.
 if (dtset%max_ncpus /=0) then
   write(ab_out,'(a)')"--- !Autoparal"
   write(ab_out,"(a)")"# Autoparal section for EPH runs"
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
   MSG_ERROR_NODUMP("aborting now")
 end if

 call cwtime(cpu, wall, gflops, "start")

 ! Construct crystal and ebands from the GS WFK file.
 if (use_wfk) then
   call wfk_read_eigenvalues(wfk0_path, gs_eigen, wfk0_hdr, comm)
   call hdr_vs_dtset(wfk0_hdr, dtset)

   cryst = hdr_get_crystal(wfk0_hdr, timrev2)
   call crystal_print(cryst,header="crystal structure from WFK file")

   ebands = ebands_from_hdr(wfk0_hdr, maxval(wfk0_hdr%nband), gs_eigen)
   ABI_FREE(gs_eigen)
 end if

 ! Read WFQ and construct ebands on the shifted grid.
 if (use_wfq) then
   call wfk_read_eigenvalues(wfq_path, gs_eigen, wfq_hdr, comm)
   ! GKA TODO: Have to construct a header with the proper set of q-shifted k-points then compare against file.
   !call hdr_vs_dtset(wfq_hdr, dtset)
   ebands_kq = ebands_from_hdr(wfq_hdr, maxval(wfq_hdr%nband), gs_eigen)
   call hdr_free(wfq_hdr)
   ABI_FREE(gs_eigen)
 end if

 ! Here we change the GS bands (Fermi level, scissors operator ...)
 ! All the modifications to ebands should be done here.
 if (use_wfk) then
   if (dtset%occopt /= ebands%occopt .or. abs(dtset%tsmear - ebands%tsmear) > tol12) then
     write(msg,"(2a,2(a,i0,a,f14.6,a))")&
     " Changing occupation scheme as input occopt and tsmear differ from those read from WFK file.",ch10,&
     "   From WFK file: occopt = ",ebands%occopt,", tsmear = ",ebands%tsmear,ch10,&
     "   From input:    occopt = ",dtset%occopt,", tsmear = ",dtset%tsmear,ch10
     call wrtout(ab_out, msg)
     call ebands_set_scheme(ebands, dtset%occopt, dtset%tsmear, dtset%spinmagntarget, prtvol=dtset%prtvol)
     if (use_wfq) then
       call ebands_set_scheme(ebands_kq, dtset%occopt, dtset%tsmear, dtset%spinmagntarget, prtvol=dtset%prtvol)
     end if
   end if

   ! Default value of eph_fermie is zero hence no tolerance is used!
   if (dtset%eph_fermie /= zero) then
     ABI_CHECK(abs(dtset%eph_extrael) <= tol12, "eph_fermie and eph_extrael are mutually exclusive")
     call wrtout(ab_out, sjoin(" Fermi level set by the user at:", ftoa(dtset%eph_fermie)))
     call ebands_set_fermie(ebands, dtset%eph_fermie, msg)
     call wrtout(ab_out, msg)
     if (use_wfq) then
       call ebands_set_fermie(ebands_kq, dtset%eph_fermie, msg)
       call wrtout(ab_out, msg)
     end if

   else if (abs(dtset%eph_extrael) > tol12) then
     call ebands_set_scheme(ebands, dtset%occopt, dtset%tsmear, dtset%spinmagntarget, dtset%prtvol)
     call ebands_set_nelect(ebands, ebands%nelect + dtset%eph_extrael, dtset%spinmagntarget, msg)
     call wrtout(ab_out, msg)
     if (use_wfq) then
       call ebands_set_scheme(ebands_kq, dtset%occopt, dtset%tsmear, dtset%spinmagntarget, dtset%prtvol)
       call ebands_set_nelect(ebands_kq, ebands%nelect + dtset%eph_extrael, dtset%spinmagntarget, msg)
       call wrtout(ab_out, msg)
     end if
   end if

   ! Recompute occupations. This is needed if WFK files have been produced in a NSCF run
   ! since occ are set to zero, and fermie is taken from the previous density.
   if (dtset%kptopt > 0) then
     call ebands_update_occ(ebands, dtset%spinmagntarget, prtvol=dtset%prtvol)
     call ebands_print(ebands,header="Ground state energies", prtvol=dtset%prtvol)
     if (use_wfq) then
       call ebands_update_occ(ebands_kq, dtset%spinmagntarget, prtvol=dtset%prtvol)
       call ebands_print(ebands_kq,header="Ground state energies (K+Q)", prtvol=dtset%prtvol)
     end if
   end if
 end if ! use_wfk

 call cwtime_report(" eph%init", cpu, wall, gflops)

 ! =======================================
 ! Output useful info on electronic bands
 ! =======================================
 call cwtime(cpu, wall, gflops, "start")
 if (my_rank == master) then
   ! Fermi Surface
   if (dtset%prtfsurf /= 0) then
     path = strcat(dtfil%filnam_ds(4), "_BXSF")
     call wrtout(ab_out, sjoin("- Writing Fermi surface to file:", path))
     if (ebands_write_bxsf(ebands,cryst,path) /= 0) then
       msg = "Cannot produce file for Fermi surface, check log file for more info"
       MSG_WARNING(msg)
       call wrtout(ab_out,msg)
     end if
   end if

   ! Nesting factor (requires qpath)
   if (dtset%prtnest /= 0 .and. dtset%ph_nqpath > 0) then
     path = strcat(dtfil%filnam_ds(4), "_NEST")
     call wrtout(ab_out, sjoin("- Writing nesting factor to file:", path))
     if (ebands_write_nesting(ebands,cryst,path,dtset%prtnest,&
         dtset%tsmear,dtset%fermie_nest,dtset%ph_qpath(:,1:dtset%ph_nqpath),msg) /= 0) then
       MSG_WARNING(msg)
       call wrtout(ab_out,msg)
     end if
   end if

   if (use_wfk) call ebands_write(ebands, dtset%prtebands, dtfil%filnam_ds(4))
 end if

 call cwtime_report(" eph%ebands_postprocess:", cpu, wall, gflops)

 ! Read the DDB file.
 ABI_CALLOC(dummy_atifc, (dtset%natom))

 ! TODO
 ! WARNING. This choice is only to insure backwards compatibility with the tests,
 ! while eph is developed. Actually, should be switched to brav1=1 as soon as possible ...
 brav1 = 1; if (dtset%eph_transport > 0) brav1 = -1

 if (use_wfk) then
   call ddb_from_file(ddb,ddb_path,brav1,dtset%natom,natifc0,dummy_atifc,cryst_ddb,comm, prtvol=dtset%prtvol)
   call cryst_ddb%free()
 else
   call ddb_from_file(ddb,ddb_path,brav1,dtset%natom,natifc0,dummy_atifc,cryst,comm, prtvol=dtset%prtvol)
 end if
 ABI_FREE(dummy_atifc)

 ! Set the q-shift for the DDB (well we mainly use gamma-centered q-meshes)
 ddb_nqshift = 1
 ABI_CALLOC(ddb_qshifts, (3, ddb_nqshift))
 ddb_qshifts(:,1) = dtset%ddb_shiftq(:)

 ! Get Dielectric Tensor
 iblock_dielt = ddb_get_dielt(ddb, dtset%rfmeth, dielt)

 ! Get Dielectric Tensor and Effective Charges
 ! (initialized to one_3D and zero if the derivatives are not available in the DDB file)
 iblock_dielt_zeff = ddb_get_dielt_zeff(ddb, cryst, dtset%rfmeth, dtset%chneut, selectz0, dielt, zeff, zeff_raw=zeff_raw)
 if (my_rank == master) then
   if (iblock_dielt_zeff == 0) then
     call wrtout(ab_out, sjoin("- Cannot find dielectric tensor and Born effective charges in DDB file:", ddb_path))
     call wrtout(ab_out, "Values initialized with zeros")
   else
     call wrtout(ab_out, sjoin("- Found dielectric tensor and Born effective charges in DDB file:", ddb_path))
   end if
 end if
 !do ii=1,3; do jj=1,3; if (ii /= jj) dielt(ii, jj) = zero; enddo; enddo

 if (any(dtset%ddb_qrefine > 1)) then
   ! Gaal-Nagy's algorithm in PRB 73 014117 [[cite:GaalNagy2006]]
   ! Build the IFCs using the coarse q-mesh.
   ngqpt_coarse = dtset%ddb_ngqpt / dtset%ddb_qrefine
   call ifc_init(ifc_coarse, cryst, ddb, &
     brav1, dtset%asr, dtset%symdynmat, dtset%dipdip, dtset%rfmeth, ngqpt_coarse, ddb_nqshift, ddb_qshifts, dielt, zeff, &
     nsphere0, rifcsph0, prtsrlr0, dtset%enunit, comm)

   ! Now use the coarse q-mesh to fill the entries in dynmat(q)
   ! on the dense q-mesh that cannot be obtained from the DDB file.
   call ifc_init(ifc, cryst, ddb, &
     brav1, dtset%asr, dtset%symdynmat, dtset%dipdip, dtset%rfmeth, dtset%ddb_ngqpt, ddb_nqshift, ddb_qshifts, dielt, zeff, &
     nsphere0, rifcsph0, prtsrlr0, dtset%enunit, comm, ifc_coarse=ifc_coarse)
   call ifc_free(ifc_coarse)

 else
   call ifc_init(ifc, cryst, ddb, &
     brav1, dtset%asr, dtset%symdynmat, dtset%dipdip, dtset%rfmeth, dtset%ddb_ngqpt, ddb_nqshift, ddb_qshifts, dielt, zeff, &
     nsphere0, rifcsph0, prtsrlr0, dtset%enunit, comm)
 end if

 ABI_FREE(ddb_qshifts)
 call ifc_print(ifc, unit=std_out)

 ! Output phonon band structure (requires qpath)
 if (dtset%prtphbands /= 0) call ifc_mkphbs(ifc, cryst, dtset, dtfil%filnam_ds(4), comm)

 if (dtset%prtphdos == 1) then
   ! Compute Phonon Density of States.
   wminmax = zero
   do
     call mkphdos(phdos, cryst, ifc, dtset%ph_intmeth, dtset%ph_wstep, dtset%ph_smear, dtset%ph_ngqpt, &
       dtset%ph_nqshift, dtset%ph_qshift, dtfil%filnam_ds(4), wminmax, count_wminmax, comm)
     if (all(count_wminmax == 0)) exit
     wminmax(1) = wminmax(1) - abs(wminmax(1)) * 0.05
     wminmax(2) = wminmax(2) + abs(wminmax(2)) * 0.05
     call phdos_free(phdos)
     write(msg, "(a, 2f8.5)")"Initial frequency mesh not large enough. Recomputing PHDOS with wmin, wmax: ",wminmax
     call wrtout(std_out, msg)
   end do

   if (my_rank == master) then
     ! Disabled by default because it's slow and we use netcdf that is much better.
     if (dtset%prtvol > 0) then
       path = strcat(dtfil%filnam_ds(4), "_PHDOS")
       call wrtout(ab_out, sjoin("- Writing phonon DOS to file:", path))
       call phdos_print(phdos, path)
       !call phdos_print_debye(phdos, cryst%ucvol)
     end if

#ifdef HAVE_NETCDF
     path = strcat(dtfil%filnam_ds(4), "_PHDOS.nc")
     !call wrtout(ab_out, sjoin("- Writing phonon DOS to netcdf file:", path))
     ncerr = nctk_open_create(ncid, path, xmpi_comm_self)
     NCF_CHECK_MSG(ncerr, sjoin("Creating PHDOS.nc file:", path))
     NCF_CHECK(cryst%ncwrite(ncid))
     call phdos_ncwrite(phdos, ncid)
     NCF_CHECK(nf90_close(ncid))
#endif
   end if
   call phdos_free(phdos)
 end if ! prtphdos

 if (dtset%prtbltztrp == 1 .and. my_rank == master) then
   call ifc_outphbtrap(ifc,cryst,dtset%ph_ngqpt,dtset%ph_nqshift,dtset%ph_qshift,dtfil%filnam_ds(4))
   ! BoltzTraP output files in GENEric format
   call ebands_prtbltztrp(ebands, cryst, dtfil%filnam_ds(4))
 end if

 ! Output phonon isosurface in Xcrysden format.
 if (dtset%prtphsurf == 1) then
   path = strcat(dtfil%filnam_ds(4), "_PH.bxsf")
   call wrtout(ab_out, sjoin("- Writing phonon frequencies in Xcrysden format to file:", path))
   call ifc_printbxsf(ifc, cryst, dtset%ph_ngqpt, dtset%ph_nqshift, dtset%ph_qshift, path, comm)
 end if

 call cwtime_report(" eph%ifc:", cpu, wall, gflops)

 ! Initialize the object used to read DeltaVscf (required if eph_task /= 0)
 if (use_dvdb) then
   dvdb = dvdb_new(dvdb_path, comm)
   if (dtset%prtvol > 10) dvdb%debug = .True.
   ! This to symmetrize the DFPT potentials.
   dvdb%symv1 = dtset%symv1scf

   ! Set qdamp from frohl_params
   if (dtset%frohl_params(4)/=0) then
     dvdb%qdamp = dtset%frohl_params(4)
   end if

   !dvdb%diel = 13.103 * dvdb%dielt
   !do ii=1,dvdb%natom
   !  dvdb%qstar(:,:,:,ii) = (-1 ** (ii + 1)) * abs(levi_civita_3()) * 13.368_dp
   !end do
   !dvdb%has_quadrupoles = .True.

   ! Set dielectric tensor, BECS and associated flags.
   ! This activates automatically the treatment of the long-range term in the Fourier interpolation
   ! of the DFPT potentials except when dvdb_add_lr == 0
   dvdb%add_lr = dtset%dvdb_add_lr
   if (iblock_dielt /= 0) then
     dvdb%has_dielt = .True.; dvdb%dielt = dielt
   end if
   if (iblock_dielt_zeff /= 0) then
     dvdb%has_zeff = .True.; dvdb%zeff = zeff; dvdb%zeff_raw = zeff_raw
   end if
   if (.not. dvdb%has_dielt .or. .not. (dvdb%has_zeff .or. dvdb%has_quadrupoles)) then
     if (dvdb%add_lr /= 0) then
       dvdb%add_lr = 0
       !call wrtout([std_out, ab_out], &
       !  " WARNING: Setting dvdb_add_lr to 0. Long-range term won't be substracted in Fourier interpolation.")
     end if
   end if

   if (my_rank == master) then
     call dvdb%print()
     call dvdb%list_perts([-1, -1, -1], unit=ab_out)
   end if

 end if

 ! TODO Recheck getng, should use same trick as that used in screening and sigma.
 call pawfgr_init(pawfgr,dtset,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfftc,ngfftf,&
 gsqcutc_eff=gsqcutc_eff,gsqcutf_eff=gsqcutf_eff,gmet=cryst%gmet,k0=k0)

 call print_ngfft(ngfftc,header='Coarse FFT mesh used for the wavefunctions')
 call print_ngfft(ngfftf,header='Dense FFT mesh used for densities and potentials')

 ! Fake MPI_type for the sequential part.
 call initmpi_seq(mpi_enreg)
 call init_distribfft_seq(mpi_enreg%distribfft,'c',ngfftc(2),ngfftc(3),'all')
 call init_distribfft_seq(mpi_enreg%distribfft,'f',ngfftf(2),ngfftf(3),'all')

!I am not sure yet the EFMAS file will be needed as soon as eph_frohlichm/=0. To be decided later.
 if (dtset%eph_frohlichm /= 0) then
#ifdef HAVE_NETCDF
   NCF_CHECK(nctk_open_read(ncid, dtfil%fnameabi_efmas, xmpi_comm_self))
   call efmas_ncread(efmasdeg, efmasval, kpt_efmas, ncid)
   NCF_CHECK(nf90_close(ncid))
#else
   MSG_ERROR("netcdf support not enabled")
#endif
 endif

 ! ===========================================
 ! === Open and read pseudopotential files ===
 ! ===========================================
 call pspini(dtset, dtfil, ecore, psp_gencond, gsqcutc_eff, gsqcutf_eff, pawrad, pawtab, psps, cryst%rprimd, comm_mpi=comm)

 ! ====================================================
 ! === This is the real epc stuff once all is ready ===
 ! ====================================================

 ! TODO: Make sure that all subdrivers work with useylm == 1
 ABI_CHECK(dtset%useylm == 0, "useylm != 0 not implemented/tested")

 ! Relase nkpt-based arrays in dtset to decreased memory requirement if dense sampling.
 ! EPH routines should not access them after this point.
 if (dtset%eph_task /= 6) call dtset_free_nkpt_arrays(dtset)

 select case (dtset%eph_task)
 case (0)
   continue

 case (1)
   ! Compute phonon linewidths in metals.
   call eph_phgamma(wfk0_path,dtfil,ngfftc,ngfftf,dtset,cryst,ebands,dvdb,ddk,ifc,&
     pawfgr,pawang,pawrad,pawtab,psps,mpi_enreg,comm)

 case (2, -2)
   ! Compute electron-phonon matrix elements
   call eph_gkk(wfk0_path,wfq_path,dtfil,ngfftc,ngfftf,dtset,cryst,ebands,ebands_kq,dvdb,ifc,&
     pawfgr,pawang,pawrad,pawtab,psps,mpi_enreg,comm)

 case (3)
   ! Compute phonon self-energy
   call eph_phpi(wfk0_path,wfq_path,dtfil,ngfftc,ngfftf,dtset,cryst,ebands,ebands_kq,dvdb,ifc,&
     pawfgr,pawang,pawrad,pawtab,psps,mpi_enreg,comm)

 case (4, -4)
   ! Compute electron self-energy (phonon contribution)
   call sigmaph(wfk0_path, dtfil, ngfftc, ngfftf, dtset, cryst, ebands, dvdb, ifc, wfk0_hdr, &
     pawfgr, pawang, pawrad, pawtab, psps, mpi_enreg, comm)

   if (dtset%eph_task == -4) call transport(dtfil, dtset, ebands, cryst, comm)

 case (5, -5)
   ! Interpolate the phonon potential
   call dvdb%interpolate_and_write(dtset, dtfil%fnameabo_dvdb, ngfftc, ngfftf, cryst, &
     ifc%ngqpt, ifc%nqshft, ifc%qshft, comm)

 case (6)
   ! Compute ZPR and temperature-dependent electronic structure using the Frohlich model
   call frohlichmodel(cryst, dtset, efmasdeg, efmasval, ifc)

 case (7)
   ! Compute phonon limited transport from SIGEPH file
   call transport(dtfil, dtset, ebands, cryst, comm)

 case (15)
   comm_rpt = xmpi_comm_self
   method = dtset%userid
   path = strcat(dtfil%filnam_ds(4), "_WRMAX")
   call wrtout(std_out, sjoin("Saving W(r, R) to file:", path))
   call dvdb%open_read(ngfftf, xmpi_comm_self)
   ! This to symmetrize the DFPT potentials.
   dvdb%symv1 = dtset%symv1scf
   call dvdb%print()
   !call dvdb%list_perts([-1, -1, -1])
   call dvdb%ftinterp_setup(dtset%ddb_ngqpt, dtset%ddb_qrefine, 1, dtset%ddb_shiftq, nfftf, ngfftf, method, path, comm_rpt)

 case (16)
   ! Compute \delta V_{q,nu)(r) and dump results to netcdf file.
   if (my_rank == master) then
     call ncwrite_v1qnu(dvdb, cryst, ifc, dvdb%nqpt, dvdb%qpts, dtset%prtvol, strcat(dtfil%filnam_ds(4), "_V1QNU.nc"))
   end if

 case default
   MSG_ERROR(sjoin("Unsupported value of eph_task:", itoa(dtset%eph_task)))
 end select

 !=====================
 !==== Free memory ====
 !=====================
 call cryst%free()
 call dvdb%free()
 call ddb_free(ddb)
 call ddk_free(ddk)
 call ifc_free(ifc)
 call hdr_free(wfk0_hdr)
 if (use_wfk) call ebands_free(ebands)
 if (use_wfq) call ebands_free(ebands_kq)
 call pawfgr_destroy(pawfgr)
 call destroy_mpi_enreg(mpi_enreg)
 if (allocated(efmasdeg)) call efmasdeg_free_array(efmasdeg)
 if (allocated(efmasval)) call efmasval_free_array(efmasval)
 ABI_SFREE(kpt_efmas)

 ! XG20180810: please do not remove. Otherwise, I get an error on my Mac.
 !write(std_out,*)' eph : after free efmasval and kpt_efmas'

 ! Deallocation for PAW.
 if (dtset%usepaw==1) then
   !call pawrhoij_free(pawrhoij)
   !ABI_DT_FREE(pawrhoij)
   !call pawfgrtab_free(pawfgrtab)
   !ABI_DT_FREE(pawfgrtab)
   !call paw_ij_free(paw_ij)
   !ABI_DT_FREE(paw_ij)
   !call paw_an_free(paw_an)
   !ABI_DT_FREE(paw_an)
 end if

 DBG_EXIT('COLL')

end subroutine eph
!!***

end module m_eph_driver
!!***
