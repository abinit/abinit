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
 use m_profiling_abi

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
 use m_profiling_abi
 use m_xmpi
 use m_xomp
 use m_errors
 use m_hdr
 use m_crystal
 use m_crystal_io
 use m_ebands
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
 use m_gkk,             only : eph_gkk, ncwrite_v1qnu
 use m_phpi,            only : eph_phpi
 use m_sigmaph,         only : sigmaph, eph_double_grid_t
 use m_ephwg,           only : ephwg_test
 use m_kpts,            only : listkk

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eph'
 use interfaces_14_hidewrite
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
 integer :: interp_kmult(3), interp_side(3), nkpt_coarse(3), nkpt_dense(3), band_block(2)
 integer :: intp_kptrlatt(3,3), intp_nshiftk
 integer :: ii,jj,kk,i1,i2,i3,i_coarse,i_dense,i_subdense,this_dense
 integer,parameter :: master=0,natifc0=0,timrev2=2,selectz0=0,sppoldbl1=1,timrev1=1
 integer,parameter :: nsphere0=0,prtsrlr0=0
 integer :: ii,comm,nprocs,my_rank,psp_gencond,mgfftf,nfftf !,nfftf_tot
 integer :: iblock,ddb_nqshift,ierr,brav1
 integer :: omp_ncpus, work_size, nks_per_proc
 real(dp):: eff,mempercpu_mb,max_wfsmem_mb,nonscal_mem !,ug_mem,ur_mem,cprj_mem
 real(dp):: params(3), intp_shiftk(3)
#ifdef HAVE_NETCDF
 integer :: ncid,ncerr
#endif
 real(dp),parameter :: rifcsph0=zero
 real(dp) :: ecore,ecut_eff,ecutdg_eff,gsqcutc_eff,gsqcutf_eff
 real(dp) :: cpu,wall,gflops
 logical :: use_wfk,use_wfq,use_dvdb,use_dg
 character(len=500) :: msg
 character(len=fnlen) :: wfk0_path,wfq_path,ddb_path,dvdb_path,path,wfk_fname_dense
 character(len=fnlen) :: ddk_path(3)
 type(hdr_type) :: wfk0_hdr, wfq_hdr
 type(crystal_t) :: cryst,cryst_ddb
 type(ebands_t) :: ebands, ebands_kq, ebands_dense
 type(ddb_type) :: ddb
 type(dvdb_t) :: dvdb
 type(ddk_t) :: ddk
 type(ifc_type) :: ifc
 type(pawfgr_type) :: pawfgr
 type(mpi_type) :: mpi_enreg
 type(phonon_dos_type) :: phdos
 type(hdr_type) :: hdr_wfk_dense
 type(eph_double_grid_t) :: eph_dg
!arrays
 integer :: ngfftc(18),ngfftf(18)
 integer,allocatable :: dummy_atifc(:)
 integer :: count_wminmax(2)
 real(dp) :: wminmax(2)
 integer,allocatable :: indqq(:,:)
 real(dp),parameter :: k0(3)=zero
 real(dp) :: dksqmax, dksqmin, dksqmean, maxfreq, error
 real(dp) :: dielt(3,3),zeff(3,3,dtset%natom), qpt(3)
 real(dp),pointer :: energies_dense(:,:,:), displ_cart(:,:,:,:), kq_kpts(:,:)
 real(dp),pointer :: gs_eigen(:,:,:) !,gs_occ(:,:,:)
 real(dp),allocatable :: ddb_qshifts(:,:), phfreq_bz(:), phfreq_ibz(:)
 !real(dp) :: tsec(2)
 !type(pawfgrtab_type),allocatable :: pawfgrtab(:)
 !type(paw_ij_type),allocatable :: paw_ij(:)
 !type(paw_an_type),allocatable :: paw_an(:)

!************************************************************************

 ! This part performs the initialization of basic objects used to perform e-ph calculations i.e:
 !
 ! 1) crystal structure `cryst`
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
 dvdb_path = dtfil%filddbsin; ii=len_trim(dvdb_path); dvdb_path(ii-2:ii+1) = "DVDB"
 use_wfk = (dtset%eph_task /= 5)
 use_wfq = (dtset%irdwfq/=0 .or. dtset%getwfq/=0)
 use_dvdb = (dtset%eph_task /= 0)
 use_dg = .false.

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
     call ebands_set_scheme(ebands, dtset%occopt, dtset%tsmear, dtset%spinmagntarget, prtvol=dtset%prtvol)
     if (use_wfq) then
       call ebands_set_scheme(ebands_kq, dtset%occopt, dtset%tsmear, dtset%spinmagntarget, prtvol=dtset%prtvol)
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

 call cwtime(cpu,wall,gflops,"stop")
 write(msg,'(2(a,f8.2))')"eph%edos: cpu:",cpu,", wall: ",wall
 call wrtout(std_out, msg, do_flush=.True.)
 call cwtime(cpu,wall,gflops,"start")

 ! Read the DDB file.
 ABI_CALLOC(dummy_atifc, (dtset%natom))

 ! WARNING. This choice is only to insure backwards compatibility with the tests,
 ! while eph is developed. Actually, should be switched to brav1=1 as soon as possible ...
 brav1 = 1; if (dtset%eph_transport > 0) brav1 = -1

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
   ! Compute Phonon Density of States.
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
!    2) I've addeded a new abinit variable (tmesh) to specifiy the list of temperatures.
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

 ! TODO Recheck getng, should use same trick as that used in screening and sigma.
 call pawfgr_init(pawfgr,dtset,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfftc,ngfftf,&
 gsqcutc_eff=gsqcutc_eff,gsqcutf_eff=gsqcutf_eff,gmet=cryst%gmet,k0=k0)

 call print_ngfft(ngfftc,header='Coarse FFT mesh used for the wavefunctions')
 call print_ngfft(ngfftf,header='Dense FFT mesh used for densities and potentials')

 ! Fake MPI_type for the sequential part.
 call initmpi_seq(mpi_enreg)
 call init_distribfft_seq(mpi_enreg%distribfft,'c',ngfftc(2),ngfftc(3),'all')
 call init_distribfft_seq(mpi_enreg%distribfft,'f',ngfftf(2),ngfftf(3),'all')

 ! ===========================================
 ! === Open and read pseudopotential files ===
 ! ===========================================
 call pspini(dtset,dtfil,ecore,psp_gencond,gsqcutc_eff,gsqcutf_eff,&
& pawrad,pawtab,psps,cryst%rprimd,comm_mpi=comm)

 ! =======================================
 ! === Prepare Double grid integration ===
 ! =======================================

 ! 1. Read kpoints and eig from a WFK/EIG/GSR file
 ! 2. Initialize a double grid object for the claculations
 ! of the integration weights on a coarse grid using the
 ! values of energy in a fine grid
 !
 ! -------------------
 ! |. . .|. . .|. . .| 
 ! |. x .|. x .|. x .| 
 ! |. . .|. . .|. . .|
 ! -------------------
 ! . = double grid
 ! x = coarse grid
 !
 ! The fine grid is used to evaluate the weights on the coarse grid
 ! The points of the fine grid are associated to the points of the
 ! coarse grid according to proximity
 ! The integration weights are returned on the coarse grid.
 !
 ! Steps of the implementation
 ! 1. Load the fine k grid from file
 ! 2. Find the matching between the k_coarse and k_dense using the double_grid object
 ! 3. Calculate the phonon frequencies on the dense mesh and store them on a array
 ! 4. Create an array to bring the points in the full brillouin zone to the irreducible brillouin zone
 ! 5. Create a scatter array between the points in the fine grid
 ! 

 !1.
 if ((dtset%getwfkfine /= 0 .and. dtset%irdwfkfine ==0) .or.&
     (dtset%getwfkfine == 0 .and. dtset%irdwfkfine /=0) )  then

   wfk_fname_dense = trim(dtfil%fnameabi_wfkfine)//'FINE'
   call wrtout(std_out,"EPH Interpolation: will read energies from: "//trim(wfk_fname_dense),"COLL")

   if (nctk_try_fort_or_ncfile(wfk_fname_dense, msg) /= 0) then
     MSG_ERROR(msg)
   end if

   call wfk_read_eigenvalues(wfk_fname_dense,energies_dense,hdr_wfk_dense,comm)
   call wfk_read_eigenvalues(wfk0_path,gs_eigen,wfk0_hdr,comm)
   ebands_dense = ebands_from_hdr(hdr_wfk_dense,maxval(hdr_wfk_dense%nband),energies_dense)
   eph_dg%ebands_dense = ebands_dense

   !TODO add a check for consistency
   ! number of bands and kpoints (comensurability)
   ABI_CHECK(hdr_wfk_dense%mband == wfk0_hdr%mband, 'Inconsistent number of bands for the fine and dense grid')

   nkpt_coarse = [wfk0_hdr%kptrlatt(1,1),&
                  wfk0_hdr%kptrlatt(2,2),&
                  wfk0_hdr%kptrlatt(3,3)]
   nkpt_dense  = [hdr_wfk_dense%kptrlatt(1,1),&
                  hdr_wfk_dense%kptrlatt(2,2),&
                  hdr_wfk_dense%kptrlatt(3,3)]
   interp_kmult = nkpt_dense/nkpt_coarse
   use_dg = .true.

 !read bs_interpmult
 else if (dtset%bs_interp_kmult(1) /= 0 .or.&
          dtset%bs_interp_kmult(2) /= 0 .or.&
          dtset%bs_interp_kmult(3) /= 0 ) then

   call wrtout(std_out,"EPH Interpolation: will use star functions interpolation","COLL")

   nkpt_coarse = [wfk0_hdr%kptrlatt(1,1),&
                  wfk0_hdr%kptrlatt(2,2),&
                  wfk0_hdr%kptrlatt(3,3)]
   interp_kmult = dtset%bs_interp_kmult
   nkpt_dense = nkpt_coarse*interp_kmult

   ! Interpolate band energies with star-functions
   params = 0; params(1) = 1; params(2) = 5
   !TODO: mband should be min of nband
   band_block = [1, ebands%mband]
   intp_kptrlatt = reshape([nkpt_dense(1), 0, 0, 0, nkpt_dense(2), 0, 0, 0, nkpt_dense(3)], [3, 3])
   intp_shiftk = zero
   intp_nshiftk = 1
   ebands_dense = ebands_interp_kmesh(ebands, cryst, params, intp_kptrlatt,&
                                      intp_nshiftk, intp_shiftk, band_block, comm)
   eph_dg%ebands_dense = ebands_dense
   use_dg = .true.
 end if

 if (use_dg) then
 !2.

   !debug
   write(std_out,*) 'dense ngkpt:',ebands_dense%nkpt
   write(std_out,*) 'dense mband:',ebands_dense%mband
   !end debug

   eph_dg%dense_nbz = nkpt_dense(1)*nkpt_dense(2)*nkpt_dense(3)
   eph_dg%coarse_nbz = nkpt_coarse(1)*nkpt_coarse(2)*nkpt_coarse(3)
   eph_dg%interp_kmult = interp_kmult
   eph_dg%nkpt_coarse = nkpt_coarse
   eph_dg%nkpt_dense = nkpt_dense
   
   ! microzone is the set of points in the fine grid belonging to a certain coarse point
   ! we have to consider a side of a certain size around the coarse point
   ! to make sure the microzone is centered around it point.
   ! The fine points shared by multiple microzones should have weights
   ! according to in how many microzones they appear
   !
   ! double grid:
   ! ----------------- interp_kmult 2
   ! |. .|.|.|.|. . .| side 1
   ! |. x|.|x|.|. x .| size 3 (2*side+1)
   ! |. .|.|.|.|. . .|
   ! -----------------
   !
   ! triple grid:
   ! ------------------- interp_kmult 3
   ! |. . .|. . .|. . .| side 1
   ! |. x .|. x .|. x .| size 3 (2*side+1)
   ! |. . .|. . .|. . .|
   ! -------------------
   !
   ! quadruple grid:
   ! --------------------------- interp_kmult 4
   ! |. . . .|.|. . .|.|. . . .| side 2
   ! |. . . .|.|. . .|.|. . . .| size 5  (2*side+1)
   ! |. . x .|.|. x .|.|. x . .|
   ! |. . . .|.|. . .|.|. . . .|
   ! |. . . .|.|. . .|.|. . . .|
   ! ---------------------------
   !
   ! and so on
   ! this is integer division
   interp_side = interp_kmult/2

   eph_dg%ndiv = (2*interp_side(1)+1)*&
                 (2*interp_side(2)+1)*&
                 (2*interp_side(3)+1)

   write(std_out,*) 'coarse:      ', nkpt_coarse
   write(std_out,*) 'dense:       ', nkpt_dense
   write(std_out,*) 'interp_kmult:', interp_kmult
   write(std_out,*) 'ndiv:        ', eph_dg%ndiv

   ABI_MALLOC(eph_dg%kpts_coarse,(3,eph_dg%coarse_nbz))
   ABI_MALLOC(eph_dg%kpts_dense,(3,eph_dg%dense_nbz))
   ABI_MALLOC(eph_dg%coarse_to_dense,(eph_dg%coarse_nbz,eph_dg%ndiv))

   ABI_MALLOC(eph_dg%dense_to_indexes,(3,eph_dg%dense_nbz))
   ABI_MALLOC(eph_dg%indexes_to_dense,(nkpt_dense(1),nkpt_dense(2),nkpt_dense(3)))

   ABI_MALLOC(eph_dg%coarse_to_indexes,(3,eph_dg%dense_nbz))
   ABI_MALLOC(eph_dg%indexes_to_coarse,(nkpt_coarse(1),nkpt_coarse(2),nkpt_coarse(3)))

   ABI_MALLOC(eph_dg%weights_dense,(eph_dg%dense_nbz))

   write(std_out,*) 'create dense to coarse mapping'
   ! generate mapping of points in dense bz to the dense bz
   ! coarse loop
   i_dense = 0
   i_coarse = 0
   do kk=1,nkpt_coarse(3)
     do jj=1,nkpt_coarse(2)
       do ii=1,nkpt_coarse(1)
         i_coarse = i_coarse + 1
         !calculate reduced coordinates of point in coarse mesh
         eph_dg%kpts_coarse(:,i_coarse) = [dble(ii-1)/nkpt_coarse(1),&
                                           dble(jj-1)/nkpt_coarse(2),&
                                           dble(kk-1)/nkpt_coarse(3)]
         !create the fine mesh
         do i3=1,interp_kmult(3)
           do i2=1,interp_kmult(2)
             do i1=1,interp_kmult(1)
               i_dense = i_dense + 1
               !calculate reduced coordinates of point in dense mesh
               eph_dg%kpts_dense(:,i_dense) =  &
                    [dble((ii-1)*interp_kmult(1)+i1-1)/(nkpt_coarse(1)*interp_kmult(1)),&
                     dble((jj-1)*interp_kmult(2)+i2-1)/(nkpt_coarse(2)*interp_kmult(2)),&
                     dble((kk-1)*interp_kmult(3)+i3-1)/(nkpt_coarse(3)*interp_kmult(3))]
               !integer indexes mapping
               eph_dg%indexes_to_dense((ii-1)*interp_kmult(1)+i1,&
                                       (jj-1)*interp_kmult(2)+i2,&
                                       (kk-1)*interp_kmult(3)+i3) = i_dense
               eph_dg%dense_to_indexes(:,i_dense) = [(ii-1)*interp_kmult(1)+i1,&
                                                     (jj-1)*interp_kmult(2)+i2,&
                                                     (kk-1)*interp_kmult(3)+i3]
             enddo
           enddo
         enddo
         eph_dg%indexes_to_coarse(ii,jj,kk) = i_coarse
         eph_dg%coarse_to_indexes(:,i_coarse) = [ii,jj,kk]
       enddo
     enddo
   enddo 

   ! here we need to iterate again because we can have points of the dense grid
   ! belonging to multiple coarse points
   i_coarse = 0
   do kk=1,nkpt_coarse(3)
     do jj=1,nkpt_coarse(2)
       do ii=1,nkpt_coarse(1)
         i_coarse = i_coarse + 1

         !create a mapping from coarse to dense
         i_subdense = 0
         do i3=-interp_side(3),interp_side(3)
           do i2=-interp_side(2),interp_side(2)
             do i1=-interp_side(1),interp_side(1)
               i_subdense = i_subdense + 1
               !integer indexes mapping
               this_dense = eph_dg%indexes_to_dense(&
                      mod((ii-1)*interp_kmult(1)+i1+nkpt_dense(1),nkpt_dense(1))+1,&
                      mod((jj-1)*interp_kmult(2)+i2+nkpt_dense(2),nkpt_dense(2))+1,&
                      mod((kk-1)*interp_kmult(3)+i3+nkpt_dense(3),nkpt_dense(3))+1)

               !array indexes mapping
               eph_dg%coarse_to_dense(i_coarse,i_subdense) = this_dense
             enddo
           enddo
         enddo
       enddo
     enddo
   enddo 

   ABI_CHECK(i_dense == eph_dg%dense_nbz, 'dense mesh mapping is incomplete') 

   !calculate the weights of each fine point
   !different methods to distribute the weights might lead to better convergence
   !loop over coarse points
   eph_dg%weights_dense = 0
   do ii=1,eph_dg%coarse_nbz
     !loop over points in the microzone
     do jj=1,eph_dg%ndiv
       i_dense = eph_dg%coarse_to_dense(ii,jj)
       eph_dg%weights_dense(i_dense) = eph_dg%weights_dense(i_dense) + 1
     end do
   end do
   !weights_dense is array, ndiv is scalar
   eph_dg%weights_dense = 1/eph_dg%weights_dense/(interp_kmult(1)*interp_kmult(2)*interp_kmult(3))

 !5.
#if 0
   write(std_out,*) 'calculate scatering'
   !from any two k and k' find q such that k - k' = q (in bz) k -> k'+q
   do ii=1,eph_dg%dense_nbz !k'
     do jj=1,eph_dg%dense_nbz !k
       !calculate indexes of k and k'
       indexes_jk = eph_dg%dense_to_indexes(:,jj) !k
       indexes_ik = eph_dg%dense_to_indexes(:,ii) !k'
       !calcualte indexes of q 
       indexes_qq = indexes_jk - indexes_ik !q = k - k'
       !bring to first bz
       indexes_qq(1) = mod(indexes_qq(1)+nkpt_dense(1),nkpt_dense(1))+1
       indexes_qq(2) = mod(indexes_qq(2)+nkpt_dense(2),nkpt_dense(2))+1
       indexes_qq(3) = mod(indexes_qq(3)+nkpt_dense(3),nkpt_dense(3))+1
       !calculate given two indexes of k and k' give index of q such that k -> k'+q
       eph_dg%scatter_dense(jj,ii) = eph_dg%indexes_to_dense(indexes_qq(1),&
                                                             indexes_qq(2),&
                                                             indexes_qq(3))
     enddo
   enddo
#endif
 !3.
   eph_dg%kpts_dense_ibz = ebands_dense%kptns
   eph_dg%dense_nibz = ebands_dense%nkpt


#if 0
   write(std_out,*) 'calculate phonon frequencies'
   !calculate the phonon frequencies at the q-points on the ibz of the dense q-grid
   ABI_MALLOC(displ_cart,(2,3,cryst%natom,3*cryst%natom))
#if 1
   !ibz version
   ! HM: I noticed that the fourier interpolation sometimes breaks the symmetries
   !     for low q-point sampling
   ABI_MALLOC(eph_dg%phfrq_dense,(3*cryst%natom,eph_dg%dense_nibz))
   do ii=1,eph_dg%dense_nibz
     qpt = eph_dg%kpts_dense_ibz(:,ii)
     ! Get phonon frequencies and displacements in reduced coordinates for this q-point
     call ifc_fourq(ifc, cryst, qpt, eph_dg%phfrq_dense(:,ii), displ_cart )
   enddo
#else
   !bz version
   ABI_MALLOC(eph_dg%phfrq_dense,(3*cryst%natom,eph_dg%dense_nbz))
   do ii=1,eph_dg%dense_nbz
     qpt = eph_dg%kpts_dense(:,ii)
     ! Get phonon frequencies and displacements in reduced coordinates for this q-point
     call ifc_fourq(ifc, cryst, qpt, eph_dg%phfrq_dense(:,ii), displ_cart )
   enddo
#endif
   ABI_FREE(displ_cart)
#endif
 
 !4. 
   write(std_out,*) 'map bz -> ibz'
   ABI_MALLOC(eph_dg%bz2ibz_dense,(eph_dg%dense_nbz))
   ABI_MALLOC(indqq,(eph_dg%dense_nbz,6))
   call listkk(dksqmax,cryst%gmet,    indqq,&
               eph_dg%kpts_dense_ibz, eph_dg%kpts_dense,&
               eph_dg%dense_nibz,     eph_dg%dense_nbz, &
               cryst%nsym,sppoldbl1,cryst%symafm,cryst%symrec,timrev1,use_symrec=.True.)
   ABI_CHECK(dksqmax < tol6, 'Problem creating a bz to ikbz kpoint mapping')
   eph_dg%bz2ibz_dense(:) = indqq(:,1)

#if 0
   !debug
   open (unit = 2, file = "ibz.dat")
   do ii=1,eph_dg%dense_nibz
     write(2,*) eph_dg%kpts_dense_ibz(:,ii)
   end do

   open (unit = 2, file = "bz.dat")
   do ii=1,eph_dg%dense_nbz
     write(2,*) eph_dg%kpts_dense(:,ii), eph_dg%bz2ibz_dense(ii)
   end do

   ABI_MALLOC(phfreq_bz,(cryst%natom*3))
   ABI_MALLOC(phfreq_ibz,(cryst%natom*3))
   ABI_MALLOC(displ_cart,(2,3,cryst%natom,3*cryst%natom))
   open (unit = 2, file = "phbz.dat")
   open (unit = 3, file = "phibz.dat")
   dksqmax = 0
   dksqmin = 1e8
   maxfreq = 0
   do ii=1,eph_dg%dense_nbz
     call ifc_fourq(ifc, cryst, eph_dg%kpts_dense(:,ii), phfreq_bz, displ_cart )
     call ifc_fourq(ifc, cryst, eph_dg%kpts_dense_ibz(:,eph_dg%bz2ibz_dense(ii)), phfreq_ibz, displ_cart )
     write(2,*) phfreq_bz
     write(3,*) phfreq_ibz
     do jj=1,cryst%natom*3
       error = abs(phfreq_bz(jj)-phfreq_ibz(jj))
       if (dksqmax < error) dksqmax = error
       if (dksqmin > error) dksqmin = error
       if (maxfreq < phfreq_bz(jj)) maxfreq = phfreq_bz(jj)
       dksqmean = dksqmean + error
     enddo
   end do
   write(*,*) 'bz2ibz phonon error min: ', dksqmin
   write(*,*) 'bz2ibz phonon error max: ', dksqmax, dksqmax/maxfreq
   write(*,*) 'bz2ibz phonon error mean:', dksqmean/eph_dg%dense_nbz, &
                                           dksqmean/eph_dg%dense_nbz/maxfreq
   ABI_FREE(displ_cart)

   open (unit = 2, file = "coarse2dense.dat")
   do ii=1,eph_dg%coarse_nbz
     write(2,*)
     write(2,*)
     write(2,*) eph_dg%kpts_coarse(:,ii)
     do jj=1,eph_dg%ndiv
       i_dense = eph_dg%coarse_to_dense(ii,jj)
       write(2,*) eph_dg%kpts_dense(:,i_dense), eph_dg%weights_dense(i_dense)
     end do
   end do

   open (unit = 2, file = "coarse.dat")
   do ii=1,eph_dg%coarse_nbz
     write(2,*) eph_dg%kpts_coarse(:,ii)
   end do
   close(2)

   open (unit = 2, file = "dense.dat")
   do ii=1,eph_dg%dense_nbz
     write(2,*) eph_dg%kpts_dense(:,ii)
   end do
   close(2)
   !end debug
#endif


 !6. we will call sigmaph here for testing purposes only
   if (dtset%eph_task == 4) then
     call sigmaph(wfk0_path,dtfil,ngfftc,ngfftf,dtset,cryst,ebands,dvdb,ifc,&
                  pawfgr,pawang,pawrad,pawtab,psps,mpi_enreg,comm,eph_dg)
   endif

   !ABI_FREE(eph_dg%phfrq_dense)
   ABI_FREE(eph_dg%weights_dense)
   ABI_FREE(eph_dg%bz2ibz_dense)
   ABI_FREE(eph_dg%kpts_coarse)
   ABI_FREE(eph_dg%kpts_dense)
   ABI_FREE(eph_dg%coarse_to_dense)
   ABI_FREE(eph_dg%dense_to_indexes)
   ABI_FREE(eph_dg%indexes_to_dense)
   ABI_FREE(eph_dg%coarse_to_indexes)
   ABI_FREE(eph_dg%indexes_to_coarse)

 end if

 ! end double grid stuff

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
   if ( .not. use_dg) then
     call sigmaph(wfk0_path,dtfil,ngfftc,ngfftf,dtset,cryst,ebands,dvdb,ifc,&
     pawfgr,pawang,pawrad,pawtab,psps,mpi_enreg,comm)
   endif

 case (5)
   ! Interpolate the phonon potential
   call dvdb_interpolate_and_write(dtfil,ngfftc,ngfftf,cryst,dvdb,&
&   ifc%ngqpt,ifc%nqshft,ifc%qshft, &
&   dtset%eph_ngqpt_fine,dtset%qptopt,mpi_enreg,comm)

 case (-4)
   call ephwg_test(dtset, cryst, ebands, ifc, dtfil%filnam_ds(4), comm)

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
 call pawfgr_destroy(pawfgr)
 call destroy_mpi_enreg(mpi_enreg)

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
