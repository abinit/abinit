!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_eph_driver
!! NAME
!!  m_eph_driver
!!
!! FUNCTION
!!   Driver for EPH calculations
!!
!! COPYRIGHT
!!  Copyright (C) 2009-2018 ABINIT group (MG, MVer,GA)
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
!!      ifc_outphbtrap,ifc_print,ifc_printbxsf,ifc_test_phinterp
!!      init_distribfft_seq,initmpi_seq,mkphdos,pawfgr_destroy,pawfgr_init
!!      phdos_free,phdos_ncwrite,phdos_print,phdos_print_thermo,print_ngfft
!!      pspini,sigmaph,wfk_read_eigenvalues,wrtout,xmpi_bcast,xmpi_end
!!
!! SOURCE

subroutine eph(acell,codvsn,dtfil,dtset,pawang,pawrad,pawtab,psps,rprim,xred)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_abicore
 use m_xmpi
 use m_xomp
 use m_errors
 use m_hdr
 use m_crystal
 use m_crystal_io
 use m_ebands
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
 use m_time,            only : cwtime
 use m_fstrings,        only : strcat, sjoin, ftoa, itoa
 use m_fftcore,         only : print_ngfft
 use m_frohlichmodel,   only : frohlichmodel
 use m_mpinfo,          only : destroy_mpi_enreg, initmpi_seq
 use m_pawang,          only : pawang_type
 use m_pawrad,          only : pawrad_type
 use m_pawtab,          only : pawtab_type, pawtab_print, pawtab_get_lsize
 use m_paw_an,          only : paw_an_type, paw_an_init, paw_an_free, paw_an_nullify
 use m_paw_ij,          only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify
 use m_pawfgrtab,       only : pawfgrtab_type, pawfgrtab_free, pawfgrtab_init
 use m_pawrhoij,        only : pawrhoij_type, pawrhoij_alloc, pawrhoij_copy, pawrhoij_free, symrhoij
 use m_pawfgr,          only : pawfgr_type, pawfgr_init, pawfgr_destroy
 use m_phgamma,         only : eph_phgamma
 use m_efmas,           only : efmasdeg_free_array, efmasval_free_array, efmas_ncread
 use m_gkk,             only : eph_gkk, ncwrite_v1qnu
 use m_phpi,            only : eph_phpi
 use m_sigmaph,         only : sigmaph

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eph'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=6),intent(in) :: codvsn
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(pawang_type),intent(inout) :: pawang
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 real(dp),intent(in) :: acell(3),rprim(3,3),xred(3,dtset%natom)
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: master=0,natifc0=0,timrev2=2,selectz0=0
 integer,parameter :: brav1=-1 ! WARNING. This choice is only to insure backwards compatibility with the tests,
!while eph is developed. Actually, should be switched to brav1=1 as soon as possible ...
 integer,parameter :: nsphere0=0,prtsrlr0=0
 integer :: ii,comm,nprocs,my_rank,psp_gencond,mgfftf,nfftf !,nfftf_tot
 integer :: iblock,ddb_nqshift,ierr
 integer :: omp_ncpus, work_size, nks_per_proc
 real(dp):: eff,mempercpu_mb,max_wfsmem_mb,nonscal_mem !,ug_mem,ur_mem,cprj_mem
#ifdef HAVE_NETCDF
 integer :: ncid,ncerr
#endif
 real(dp),parameter :: rifcsph0=zero
 real(dp) :: ecore,ecut_eff,ecutdg_eff,gsqcutc_eff,gsqcutf_eff
 real(dp) :: cpu,wall,gflops
 logical :: use_wfk,use_wfq,use_dvdb
 character(len=500) :: msg
 character(len=fnlen) :: wfk0_path,wfq_path,ddb_path,dvdb_path,efmas_path,path
 character(len=fnlen) :: ddk_path(3)
 type(hdr_type) :: wfk0_hdr, wfq_hdr
 type(crystal_t) :: cryst,cryst_ddb
 type(ebands_t) :: ebands, ebands_kq
 type(ddb_type) :: ddb
 type(dvdb_t) :: dvdb
 type(ddk_t) :: ddk
 type(ifc_type) :: ifc
 type(pawfgr_type) :: pawfgr
 type(mpi_type) :: mpi_enreg
 type(phonon_dos_type) :: phdos
!arrays
 integer :: ngfftc(18),ngfftf(18)
 integer,allocatable :: dummy_atifc(:)
 integer :: count_wminmax(2)
 real(dp),parameter :: k0(3)=zero
 real(dp) :: dielt(3,3),zeff(3,3,dtset%natom)
 real(dp) :: wminmax(2)
 real(dp),pointer :: gs_eigen(:,:,:) !,gs_occ(:,:,:)
 real(dp),allocatable :: ddb_qshifts(:,:)
 real(dp),allocatable :: kpt_efmas(:,:)
 type(efmasdeg_type),allocatable :: efmasdeg(:)
 type(efmasval_type),allocatable :: efmasval(:,:)
 !real(dp) :: tsec(2)
 !type(pawfgrtab_type),allocatable :: pawfgrtab(:)
 !type(paw_ij_type),allocatable :: paw_ij(:)
 !type(paw_an_type),allocatable :: paw_an(:)

!************************************************************************

 ! This part performs the initialization of basic objects used to perform e-ph calculations i.e:
 !
 ! 1) Crystal structure `cryst`
 ! 2) Ground state band energies: `ebands`
 ! 3) Interatomic force constants: `ifc`
 ! 4) DVDB database with the dvscf potentials
 ! 5) Pseudos and PAW basic objects.
 !
 ! Once we have these objects, we can call specialized routines for e-ph calculations.
 ! Notes:
 !   * Any modification to the basic objects mentioned above should be done here (e.g. change of efermi)
 !   * This routines shall not allocate big chunks of memory. The CPU-demanding sections should be
 !     performed in the specialized routines that will employ different MPI distribution schemes.

 DBG_ENTER('COLL')

 if (psps%usepaw == 1) then
   MSG_ERROR("PAW not implemented")
   ABI_UNUSED((/pawang%nsym, pawrad(1)%mesh_size/))
 end if

 ! abirules!
 if (.False.) write(std_out,*)acell,codvsn,rprim,xred

 comm = xmpi_world; nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 ! Initialize filenames
 wfk0_path = dtfil%fnamewffk
 wfq_path = dtfil%fnamewffq
 ddb_path = dtfil%filddbsin
 efmas_path = dtfil%fnameabi_efmas
 dvdb_path = dtfil%filddbsin; ii=len_trim(dvdb_path); dvdb_path(ii-2:ii+1) = "DVDB"
 use_wfk = (dtset%eph_task /= 5)
 use_wfq = (dtset%irdwfq/=0 .or. dtset%getwfq/=0 .and. dtset%eph_frohlichm/=1)
 use_dvdb = (dtset%eph_task /= 0  .and. dtset%eph_frohlichm/=1)

 if(dtset%eph_frohlichm/=1)then
   efmas_path = dtfil%fnameabi_efmas
 endif

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
 if (dtset%eph_frohlichm/=1) then
   call xmpi_bcast(efmas_path,master,comm,ierr)
   call wrtout(ab_out, sjoin("- Reading EFMAS information from file:", efmas_path) )
 end if

 ! autoparal section
 ! TODO: This just to activate autoparal in abipy. Lot of things should be improved.
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

 call cwtime(cpu,wall,gflops,"start")

 ! Construct crystal and ebands from the GS WFK file.
 if (use_wfk) then
   call wfk_read_eigenvalues(wfk0_path,gs_eigen,wfk0_hdr,comm) !,gs_occ)
   call hdr_vs_dtset(wfk0_hdr,dtset)

   call crystal_from_hdr(cryst,wfk0_hdr,timrev2)
   call crystal_print(cryst,header="crystal structure from WFK file")

   ebands = ebands_from_hdr(wfk0_hdr,maxval(wfk0_hdr%nband),gs_eigen)
   call hdr_free(wfk0_hdr)
   ABI_FREE(gs_eigen)
 end if

 ! Read WFQ and construct ebands on the shifted grid.
 if (use_wfq) then
   call wfk_read_eigenvalues(wfq_path,gs_eigen,wfq_hdr,comm) !,gs_occ)
   ! GKA TODO: Have to construct a header with the proper set of q-shifted k-points then compare against file.
   !call hdr_vs_dtset(wfq_hdr,dtset)
   ebands_kq = ebands_from_hdr(wfq_hdr,maxval(wfq_hdr%nband),gs_eigen)
   call hdr_free(wfq_hdr)
   ABI_FREE(gs_eigen)
 end if

 ! Here we change the GS bands (fermi level, scissors operator ...)
 ! All the modifications to ebands should be done here.
 if (use_wfk) then

   if (dtset%occopt /= ebands%occopt .or. abs(dtset%tsmear - ebands%tsmear) > tol12) then
     write(msg,"(2a,2(a,i0,a,f14.6,a))")&
     " Changing occupation scheme as input occopt and tsmear differ from those read from WFK file.",ch10,&
     "   From WFK file: occopt = ",ebands%occopt,", tsmear = ",ebands%tsmear,ch10,&
     "   From input:    occopt = ",dtset%occopt,", tsmear = ",dtset%tsmear,ch10
     call wrtout(ab_out,msg)
     call ebands_set_scheme(ebands,dtset%occopt,dtset%tsmear,dtset%spinmagntarget,dtset%prtvol)
     if (use_wfq) then
       call ebands_set_scheme(ebands_kq,dtset%occopt,dtset%tsmear,dtset%spinmagntarget,dtset%prtvol)
     end if
   end if

   if (dtset%eph_fermie /= zero) then ! default value of eph_fermie is zero hence no tolerance is used!
     ABI_CHECK(abs(dtset%eph_extrael) <= tol12, "eph_fermie and eph_extrael are mutually exclusive")
     call wrtout(ab_out, sjoin(" Fermi level set by the user at:",ftoa(dtset%eph_fermie)))
     call ebands_set_fermie(ebands, dtset%eph_fermie, msg)
     call wrtout(ab_out,msg)
     if (use_wfq) then
       call ebands_set_fermie(ebands_kq, dtset%eph_fermie, msg)
       call wrtout(ab_out,msg)
     end if

   else if (abs(dtset%eph_extrael) > tol12) then
     !NOT_IMPLEMENTED_ERROR()
     ! TODO: Be careful with the trick used in elphon for passing the concentration
     call ebands_set_scheme(ebands, dtset%occopt, dtset%tsmear, dtset%spinmagntarget, dtset%prtvol)
     call ebands_set_nelect(ebands, dtset%eph_extrael, dtset%spinmagntarget, msg)
     call wrtout(ab_out,msg)
     if (use_wfq) then
       call ebands_set_scheme(ebands_kq, dtset%occopt, dtset%tsmear, dtset%spinmagntarget, dtset%prtvol)
       call ebands_set_nelect(ebands_kq, dtset%eph_extrael, dtset%spinmagntarget, msg)
       call wrtout(ab_out,msg)
     end if
   end if

   ! Recompute occupations. This is needed if WFK files have been produced in a NSCF run
   ! since occ are set to zero, and fermie is taken from the previous density.
   if (dtset%kptopt > 0) then
     call ebands_update_occ(ebands, dtset%spinmagntarget, prtvol=dtset%prtvol)
     call ebands_print(ebands,header="Ground state energies",prtvol=dtset%prtvol)
     if (use_wfq) then
       call ebands_update_occ(ebands_kq, dtset%spinmagntarget, prtvol=dtset%prtvol)
       call ebands_print(ebands_kq,header="Ground state energies (K+Q)", prtvol=dtset%prtvol)
     end if
   end if

 end if

 call cwtime(cpu,wall,gflops,"stop")
 write(msg,'(2(a,f8.2))')"eph%init: cpu: ",cpu,", wall: ",wall
 call wrtout(std_out, msg, do_flush=.True.)
 call cwtime(cpu,wall,gflops,"start")

 ! =======================================
 ! Output useful info on electronic bands
 ! =======================================
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

   if (use_wfk) then
     call ebands_write(ebands, dtset%prtebands, dtfil%filnam_ds(4))
   end if
 end if

 !if (.False.) then
 !!if (.True.) then
 !  !call ebands_set_interpolator(ebands, cryst, bstart, bcount, mode, espline_ords, eskw_ratio, comm)
 !  call ebands_test_interpolator(ebands, dtset, cryst, dtfil%filnam_ds(4), comm)
 !  MSG_ERROR("interpolation done")
 !end if

 call cwtime(cpu,wall,gflops,"stop")
 write(msg,'(2(a,f8.2))')"eph%edos: cpu:",cpu,", wall: ",wall
 call wrtout(std_out, msg, do_flush=.True.)
 call cwtime(cpu,wall,gflops,"start")

 ! Read the DDB file.
 ABI_CALLOC(dummy_atifc, (dtset%natom))

 if (use_wfk) then
   call ddb_from_file(ddb,ddb_path,brav1,dtset%natom,natifc0,dummy_atifc,cryst_ddb,comm, prtvol=dtset%prtvol)
   call crystal_free(cryst_ddb)
 else
   call ddb_from_file(ddb,ddb_path,brav1,dtset%natom,natifc0,dummy_atifc,cryst,comm, prtvol=dtset%prtvol)
 end if
 ABI_FREE(dummy_atifc)

 ddb_nqshift = 1
 ABI_CALLOC(ddb_qshifts, (3,ddb_nqshift))

 ! Set the q-shift for the DDB
 ddb_qshifts(:,1) = dtset%ddb_shiftq(:)

 ! Get Dielectric Tensor and Effective Charges
 ! (initialized to one_3D and zero if the derivatives are not available in the DDB file)
 iblock = ddb_get_dielt_zeff(ddb,cryst,dtset%rfmeth,dtset%chneut,selectz0,dielt,zeff)
 if (my_rank == master) then
   if (iblock == 0) then
     call wrtout(ab_out, sjoin("- Cannot find dielectric tensor and Born effective charges in DDB file:", ddb_path))
     call wrtout(ab_out, "Values initialized with zeros")
   else
     call wrtout(ab_out, sjoin("- Found dielectric tensor and Born effective charges in DDB file:", ddb_path))
   end if
 end if

 ! Build the inter-atomic force constants.
 ! WARNING : brav1 has been set to -1 at the initialisation of eph.F90, see the message there. Should be turned to 1 as soon as possible.
 call ifc_init(ifc,cryst,ddb,&
 brav1,dtset%asr,dtset%symdynmat,dtset%dipdip,dtset%rfmeth,dtset%ddb_ngqpt,ddb_nqshift,ddb_qshifts,dielt,zeff,&
 nsphere0,rifcsph0,prtsrlr0,dtset%enunit,comm)
 ABI_FREE(ddb_qshifts)
 call ifc_print(ifc, unit=std_out)

 ! Test B-spline interpolation of phonons
 !if (.True.) then
 if (.False.) then
   call ifc_test_phinterp(ifc, cryst, [8,8,8], 1, [zero,zero,zero], [3,3,3], comm)
   !call ifc_set_interpolator(ifc, cryst, nustart, nucount, mode, phspline_ords, phskw_ratio, comm)
   !call ifc_test_intepolator(ifc, dtset, dtfil, comm)
   call xmpi_end()
 end if

 ! Output phonon band structure (requires qpath)
 if (dtset%prtphbands /= 0) call ifc_mkphbs(ifc, cryst, dtset, dtfil%filnam_ds(4), comm)

 if (dtset%prtphdos == 1) then

   ! Phonon Density of States.
   wminmax = zero
   do
     call mkphdos(phdos, cryst, ifc, dtset%ph_intmeth, dtset%ph_wstep, dtset%ph_smear, dtset%ph_ngqpt, &
       dtset%ph_nqshift, dtset%ph_qshift, wminmax, count_wminmax, comm)
     if (all(count_wminmax == 0)) exit
     wminmax(1) = wminmax(1) - abs(wminmax(1)) * 0.05
     wminmax(2) = wminmax(2) + abs(wminmax(2)) * 0.05
     call phdos_free(phdos)
   end do

   if (my_rank == master) then
     path = strcat(dtfil%filnam_ds(4), "_PHDOS")
     call wrtout(ab_out, sjoin("- Writing phonon DOS to file:", path))
     call phdos_print(phdos, path)

     !call phdos_print_debye(phdos, crystal%ucvol)

!TODO: do we want to pass the temper etc... from anaddb_dtset into the full dtset for abinit?
! Otherwise just leave these defaults.
!MG: 1) Disabled for the time being because of SIGFPE in v8[41]
!    2) I've added a new abinit variable (tmesh) to specifiy the list of temperatures.
     path = strcat(dtfil%filnam_ds(4), "_MSQD_T")
!MG: Disabled for the time being because of SIGFPE in v8[41]
     !call phdos_print_msqd(phdos, path, 1000, one, one)
     path = strcat(dtfil%filnam_ds(4), "_THERMO")
     call phdos_print_thermo(PHdos, path, 1000, zero, one)

#ifdef HAVE_NETCDF
     path = strcat(dtfil%filnam_ds(4), "_PHDOS.nc")
     ncerr = nctk_open_create(ncid, path, xmpi_comm_self)
     NCF_CHECK_MSG(ncerr, sjoin("Creating PHDOS.nc file:", path))
     NCF_CHECK(crystal_ncwrite(cryst, ncid))
     call phdos_ncwrite(phdos, ncid)
     NCF_CHECK(nf90_close(ncid))
#endif
   end if
   call phdos_free(phdos)
 end if

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

 call cwtime(cpu,wall,gflops,"stop")
 write(msg,'(2(a,f8.2))')"eph%ifc: cpu:",cpu,", wall: ",wall
 call wrtout(std_out, msg, do_flush=.True.)
 call cwtime(cpu,wall,gflops,"start")

 ! Initialize the object used to read DeltaVscf (required if eph_task /= 0)
 if (use_dvdb) then
   call dvdb_init(dvdb, dvdb_path, comm)
   if (my_rank == master) then
     call dvdb_print(dvdb)
     call dvdb_list_perts(dvdb, [-1,-1,-1], unit=ab_out)
   end if

   if (iblock /= 0) then
     dvdb%dielt = dielt
     dvdb%zeff = zeff
     dvdb%has_dielt_zeff = .True.
   end if

   ! Compute \delta V_{q,nu)(r) and dump results to netcdf file.
   if (.False. .and. my_rank == master) then
     call ncwrite_v1qnu(dvdb, cryst, ifc, dvdb%nqpt, dvdb%qpts, dtset%prtvol, strcat(dtfil%filnam_ds(4), "_V1QNU.nc"))
   end if
 end if

 if(dtset%eph_frohlichm/=1)then

   ! TODO Recheck getng, should use same trick as that used in screening and sigma.
   call pawfgr_init(pawfgr,dtset,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfftc,ngfftf,&
   gsqcutc_eff=gsqcutc_eff,gsqcutf_eff=gsqcutf_eff,gmet=cryst%gmet,k0=k0)

   call print_ngfft(ngfftc,header='Coarse FFT mesh used for the wavefunctions')
   call print_ngfft(ngfftf,header='Dense FFT mesh used for densities and potentials')

   ! Fake MPI_type for the sequential part.
   call initmpi_seq(mpi_enreg)
   call init_distribfft_seq(mpi_enreg%distribfft,'c',ngfftc(2),ngfftc(3),'all')
   call init_distribfft_seq(mpi_enreg%distribfft,'f',ngfftf(2),ngfftf(3),'all')

 endif

!I am not sure yet the EFMAS file will be needed as soon as eph_frohlichm/=0. To be decided later.
 if(dtset%eph_frohlichm/=0)then
   
#ifdef HAVE_NETCDF
   NCF_CHECK(nctk_open_read(ncid, efmas_path, xmpi_comm_self))
   call efmas_ncread(efmasdeg,efmasval,kpt_efmas,ncid)
   NCF_CHECK(nf90_close(ncid))
#else
   MSG_ERROR("netcdf support not enabled")
#endif

 endif

 ! ===========================================
 ! === Open and read pseudopotential files ===
 ! ===========================================

 if(dtset%eph_frohlichm/=1)then

   call pspini(dtset,dtfil,ecore,psp_gencond,gsqcutc_eff,gsqcutf_eff,&
&   pawrad,pawtab,psps,cryst%rprimd,comm_mpi=comm)

 endif

 ! ====================================================
 ! === This is the real epc stuff once all is ready ===
 ! ====================================================
! TODO: decide whether to make several driver functions.
!  before that, however, need to encapsulate most of the functionalities in eph_phgamma
!  otherwise there will be tons of duplicated code

 ! TODO: Make sure that all subdrivers work with useylm == 1
 ABI_CHECK(dtset%useylm == 0, "useylm != 0 not implemented/tested")
 select case (dtset%eph_task)
 case (0)
   continue

 case (1)
   ! Compute phonon linewidths in metals.
   call eph_phgamma(wfk0_path,dtfil,ngfftc,ngfftf,dtset,cryst,ebands,dvdb,ddk,ifc,&
   pawfgr,pawang,pawrad,pawtab,psps,mpi_enreg,comm)

 case (2)
   ! Compute electron-phonon matrix elements
   call eph_gkk(wfk0_path,wfq_path,dtfil,ngfftc,ngfftf,dtset,cryst,ebands,ebands_kq,dvdb,ifc,&
   pawfgr,pawang,pawrad,pawtab,psps,mpi_enreg,comm)

 case (3)
   ! Compute phonon self-energy
   call eph_phpi(wfk0_path,wfq_path,dtfil,ngfftc,ngfftf,dtset,cryst,ebands,ebands_kq,dvdb,ifc,&
   pawfgr,pawang,pawrad,pawtab,psps,mpi_enreg,comm)

 case (4)
   ! Compute electron self-energy (phonon contribution)
   call sigmaph(wfk0_path,dtfil,ngfftc,ngfftf,dtset,cryst,ebands,dvdb,ifc,&
   pawfgr,pawang,pawrad,pawtab,psps,mpi_enreg,comm)

 case (5)
   ! Interpolate the phonon potential
   call dvdb_interpolate_and_write(dtfil,ngfftc,ngfftf,cryst,dvdb,&
&   ifc%ngqpt,ifc%nqshft,ifc%qshft, &
&   dtset%eph_ngqpt_fine,dtset%qptopt,mpi_enreg,comm)

 case (6)
   ! Compute ZPR and temperature-dependent electronic structure using the Frohlich model

   call frohlichmodel(cryst,dtfil,dtset,ebands,efmasdeg,efmasval,ifc)

 case default
   MSG_ERROR(sjoin("Unsupported value of eph_task:", itoa(dtset%eph_task)))
 end select

 !=====================
 !==== Free memory ====
 !=====================

 call crystal_free(cryst)
 call dvdb_free(dvdb)
 call ddb_free(ddb)
 call ddk_free(ddk)
 call ifc_free(ifc)
 if (use_wfk) call ebands_free(ebands)
 if (use_wfq) call ebands_free(ebands_kq)
 if(dtset%eph_frohlichm/=1)then
   call pawfgr_destroy(pawfgr)
   call destroy_mpi_enreg(mpi_enreg)
 endif
 if(allocated(efmasdeg))then
   call efmasdeg_free_array(efmasdeg)
 endif
 if( allocated (efmasval))then
   call efmasval_free_array(efmasval)
 endif
 if(allocated(kpt_efmas))then
   ABI_DEALLOCATE(kpt_efmas)
 endif

!XG20180810: please do not remove. Otherwise, I get an error on my Mac.
 write(std_out,*)' eph : after free efmasval and kpt_efmas'

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
