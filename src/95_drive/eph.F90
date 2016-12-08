!{\src2tex{textfont=tt}}
!!****f* ABINIT/eph
!! NAME
!!  eph
!!
!! FUNCTION
!! Main routine that computes the electron phonon coupling matrix elements and
!! calculates related properties - Tc, phonon linewidths...
!!
!! COPYRIGHT
!! Copyright (C) 2009-2016 ABINIT group (MG, MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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
!!      dvdb_list_perts,dvdb_print,ebands_free,ebands_print,ebands_set_fermie
!!      ebands_set_scheme,ebands_update_occ,edos_free,edos_init,edos_write
!!      eph_phgamma,hdr_free,hdr_vs_dtset,ifc_free,ifc_init,ifc_outphbtrap
!!      init_distribfft_seq,initmpi_seq,mkphdos,pawfgr_destroy,pawfgr_init
!!      phdos_free,phdos_print,print_ngfft,ebands_prtbltztrp,pspini
!!      wfk_read_eigenvalues,wrtout,xmpi_bcast
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine eph(acell,codvsn,dtfil,dtset,pawang,pawrad,pawtab,psps,rprim,xred)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
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
 use m_mpinfo,          only : destroy_mpi_enreg
 use m_pawang,          only : pawang_type
 use m_pawrad,          only : pawrad_type
 use m_pawtab,          only : pawtab_type, pawtab_print, pawtab_get_lsize
 use m_paw_an,          only : paw_an_type, paw_an_init, paw_an_free, paw_an_nullify
 use m_paw_ij,          only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify
 use m_pawfgrtab,       only : pawfgrtab_type, pawfgrtab_free, pawfgrtab_init
 use m_pawrhoij,        only : pawrhoij_type, pawrhoij_alloc, pawrhoij_copy, pawrhoij_free, symrhoij
 !use m_pawdij,          only : pawdij, symdij
 use m_pawfgr,          only : pawfgr_type, pawfgr_init, pawfgr_destroy
 use m_phgamma,         only : eph_phgamma
 use m_gkk,             only : eph_gkk
 use m_phpi,            only : eph_phpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eph'
 use interfaces_14_hidewrite
 use interfaces_51_manage_mpi
 use interfaces_56_io_mpi
 use interfaces_64_psp
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
 integer,parameter :: master=0,level40=40,natifc0=0,brav1=1,timrev2=2,selectz0=0
 integer,parameter :: nsphere0=0,prt_ifc0=0,ifcana0=0,ifcout0=0,prtsrlr0=0
 integer :: ii,comm,nprocs,my_rank,psp_gencond,mgfftf,nfftf !,nfftf_tot
 integer :: iblock,ddb_nqshift,ierr,edos_intmeth
#ifdef HAVE_NETCDF
 integer :: ncid,ncerr
#endif
 real(dp),parameter :: rifcsph0=zero,spinmagntarget=-99.99_dp
 real(dp) :: ecore,ecut_eff,ecutdg_eff,gsqcutc_eff,gsqcutf_eff
 real(dp) :: edos_step,edos_broad
 real(dp) :: cpu,wall,gflops
 character(len=500) :: msg
 character(len=fnlen) :: wfk0_path,wfq_path,ddb_path,dvdb_path,path
 character(len=fnlen) :: ddk_path(3)
 character(len=10) :: strddk
 type(hdr_type) :: wfk0_hdr, wfq_hdr
 type(crystal_t) :: cryst,cryst_ddb
 type(ebands_t) :: ebands, ebands_kq
 type(edos_t) :: edos
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
 real(dp),parameter :: k0(3)=zero
 real(dp) :: dielt(3,3),zeff(3,3,dtset%natom),n0(dtset%nsppol)
 real(dp),pointer :: gs_eigen(:,:,:) !,gs_occ(:,:,:)
 real(dp),allocatable :: ddb_qshifts(:,:)
 logical :: usewfq
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
 usewfq = (dtset%irdwfq/=0 .or. dtset%getwfq/=0)

 ddk_path(1) = strcat(dtfil%fnamewffddk, itoa(3*dtset%natom+1))
 ddk_path(2) = strcat(dtfil%fnamewffddk, itoa(3*dtset%natom+2))
 ddk_path(3) = strcat(dtfil%fnamewffddk, itoa(3*dtset%natom+3))

 if (my_rank == master) then
   if (.not. file_exists(ddb_path)) MSG_ERROR(sjoin("Cannot find DDB file:", ddb_path))
   if (.not. file_exists(dvdb_path)) MSG_ERROR(sjoin("Cannot find DVDB file:", dvdb_path))

   ! Accept WFK file in Fortran or netcdf format.
   if (nctk_try_fort_or_ncfile(wfk0_path, msg) /= 0) then
     MSG_ERROR(sjoin("Cannot find GS WFK file:", wfk0_path, msg))
   end if
   ! WFQ file
   if (usewfq) then
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
 call xmpi_bcast(wfk0_path,master,comm,ierr)
 call wrtout(ab_out, sjoin("- Reading GS states from WFK file:", wfk0_path))
 if (usewfq) then
   call xmpi_bcast(wfq_path,master,comm,ierr)
   call wrtout(ab_out, sjoin("- Reading GS states from WFQ file:", wfq_path) )
 end if
 call wrtout(ab_out, sjoin("- Reading DDB from file:", ddb_path))
 call wrtout(ab_out, sjoin("- Reading DVDB from file:", dvdb_path))
 if (dtset%eph_transport > 0) then
   call xmpi_bcast(ddk_path,master,comm,ierr)
   call wrtout(ab_out, sjoin("- Reading DDK x from file:", ddk_path(1)))
   call wrtout(ab_out, sjoin("- Reading DDK y from file:", ddk_path(2)))
   call wrtout(ab_out, sjoin("- Reading DDK z from file:", ddk_path(3)))
   ! Read header in DDK files and init basic dimensions.
   ! subdrivers will use ddk to get the matrix elements from file.
   ! TODO: Should perform consistency check
   call ddk_init(ddk, ddk_path, comm)
 end if

 call cwtime(cpu,wall,gflops,"start")

 ! Construct crystal and ebands from the GS WFK file.
 call wfk_read_eigenvalues(wfk0_path,gs_eigen,wfk0_hdr,comm) !,gs_occ)
 call hdr_vs_dtset(wfk0_hdr,dtset)

 call crystal_from_hdr(cryst,wfk0_hdr,timrev2)
 call crystal_print(cryst,header="crystal structure from WFK file")

 ebands = ebands_from_hdr(wfk0_hdr,maxval(wfk0_hdr%nband),gs_eigen)
 call hdr_free(wfk0_hdr)
 ABI_FREE(gs_eigen)

 ! TODO:
 ! Make sure everything is OK if WFK comes from a NSCF run since occ are set to zero
 ! fermie is set to 0 if nscf!

 ! Read WFQ and construct ebands on the shifted grid.
 if (usewfq) then
   call wfk_read_eigenvalues(wfq_path,gs_eigen,wfq_hdr,comm) !,gs_occ)
   !call hdr_vs_dtset(wfq_hdr,dtset)  ! GKA TODO: Have to construct a header with the proper set of q-shifted k-points
                                      !           then compare against file.
   ebands_kq = ebands_from_hdr(wfq_hdr,maxval(wfq_hdr%nband),gs_eigen)
   call hdr_free(wfq_hdr)
   ABI_FREE(gs_eigen)
 end if

 ! Here we change the GS bands (fermi level, scissors operator ...)
 ! All the modifications to ebands should be done here.

 if (dtset%occopt /= ebands%occopt .or. abs(dtset%tsmear - ebands%tsmear) > tol12) then
   write(msg,"(2a,2(a,i0,a,f14.6,a))")&
   " Changing occupation scheme as input occopt and tsmear differ from those read from WFK file.",ch10,&
   "   From WFK file: occopt = ",ebands%occopt,", tsmear = ",ebands%tsmear,ch10,&
   "   From input:    occopt = ",dtset%occopt,", tsmear = ",dtset%tsmear,ch10
   call wrtout(ab_out,msg)
   call ebands_set_scheme(ebands,dtset%occopt,dtset%tsmear,spinmagntarget,dtset%prtvol)
   if (usewfq) then
     call ebands_set_scheme(ebands_kq,dtset%occopt,dtset%tsmear,spinmagntarget,dtset%prtvol)
   end if
 end if

 if (dtset%eph_fermie /= zero) then ! default value of eph_fermie is zero hence no tolerance is used!
   ABI_CHECK(abs(dtset%eph_extrael) <= tol12, "eph_fermie and eph_extrael are mutually exclusive")
   call wrtout(ab_out, sjoin(" Fermi level set by the user at:",ftoa(dtset%eph_fermie)))
   call ebands_set_fermie(ebands, dtset%eph_fermie, msg)
   call wrtout(ab_out,msg)
   if (usewfq) then
     call ebands_set_fermie(ebands_kq, dtset%eph_fermie, msg)
     call wrtout(ab_out,msg)
   end if

 else if (abs(dtset%eph_extrael) > tol12) then
   NOT_IMPLEMENTED_ERROR()
   ! TODO: Be careful with the trick used in elphon for passing the concentration
   !call ebands_set_nelect(ebands, dtset%eph_extrael, spinmagntarget, msg)
   call wrtout(ab_out,msg)
 end if

 call ebands_update_occ(ebands, spinmagntarget)
 call ebands_print(ebands,header="Ground state energies",prtvol=dtset%prtvol)
 if (usewfq) then
   call ebands_update_occ(ebands_kq, spinmagntarget)
   call ebands_print(ebands_kq,header="Ground state energies (K+Q)",prtvol=dtset%prtvol)
 end if

 call cwtime(cpu,wall,gflops,"stop")
 write(msg,'(2(a,f8.2))')"eph%init: cpu: ",cpu,", wall: ",wall
 call wrtout(std_out,msg,"COLL",do_flush=.True.)
 call cwtime(cpu,wall,gflops,"start")

 ! Compute electron DOS.
 ! TODO: Optimize this part. Really slow if tetra and lots of points
 ! Could just do DOS around efermi
 edos_intmeth = 2; if (dtset%prtdos == 1) edos_intmeth = 1
 edos_step = dtset%dosdeltae; edos_broad = dtset%tsmear
 edos_step = 0.01 * eV_Ha; edos_broad = 0.3 * eV_Ha
 call edos_init(edos,ebands,cryst,edos_intmeth,edos_step,edos_broad,comm,ierr)
 ABI_CHECK(ierr==0, "Error in edos_init, see message above.")

 ! Store DOS per spin channels
 n0(:) = edos%gef(1:edos%nsppol)
 if (my_rank == master) then
   call edos_print(edos, unit=ab_out)
   path = strcat(dtfil%filnam_ds(4), "_EDOS")
   call wrtout(ab_out, sjoin("- Writing electron DOS to file:", path))
   call edos_write(edos, path)
 end if

 call edos_free(edos)

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
       call wrtout(ab_out,msg,'COLL')
     end if
   end if

   ! Nesting factor (requires qpath)
   if (dtset%prtnest /= 0 .and. dtset%ph_nqpath > 0) then
     path = strcat(dtfil%filnam_ds(4), "_NEST")
     call wrtout(ab_out, sjoin("- Writing nesting factor to file:", path))
     if (ebands_write_nesting(ebands,cryst,path,dtset%prtnest,&
     dtset%tsmear,dtset%fermie_nest,dtset%ph_qpath(:,1:dtset%ph_nqpath),msg) /= 0) then
       MSG_WARNING(msg)
       call wrtout(ab_out,msg,'COLL')
     end if
   end if
 end if ! master

 call cwtime(cpu,wall,gflops,"stop")
 write(msg,'(2(a,f8.2))')"eph%edos: cpu:",cpu,", wall: ",wall
 call wrtout(std_out,msg,"COLL",do_flush=.True.)
 call cwtime(cpu,wall,gflops,"start")

 ! Read the DDB file.
 ABI_CALLOC(dummy_atifc, (cryst%natom))

 call ddb_from_file(ddb,ddb_path,brav1,cryst%natom,natifc0,dummy_atifc,cryst_ddb,comm)
 call crystal_free(cryst_ddb)
 ABI_FREE(dummy_atifc)

 ddb_nqshift = 1
 ABI_CALLOC(ddb_qshifts, (3,ddb_nqshift))

 ! Set the q-shift for the DDB
 ddb_qshifts(:,1) = dtset%ddb_shiftq(:)

 ! Get Dielectric Tensor and Effective Charges
 ! (initialized to one_3D and zero if the derivatives are not available in the DDB file)
 iblock = ddb_get_dielt_zeff(ddb,cryst,dtset%rfmeth,dtset%chneut,selectz0,dielt,zeff)
 if (iblock == 0) then
   call wrtout(std_out,"DDB does not contain the dielectric tensor and the effective charges. Init with zeros", "COLL")
 end if

 ! Build the inter-atomic force constants.
 call ifc_init(ifc,cryst,ddb,&
 brav1,dtset%asr,dtset%symdynmat,dtset%dipdip,dtset%rfmeth,dtset%ddb_ngqpt,ddb_nqshift,ddb_qshifts,dielt,zeff,&
 nsphere0,rifcsph0,prtsrlr0,dtset%enunit)
 ABI_FREE(ddb_qshifts)

 if (dtset%prtphdos == 1) then
   ! Phonon Density of States.

   ! FIXME: mkphdos expects qshift(3) instead of qshift(3, nqshift)
   call mkphdos(phdos,cryst,ifc,dtset%ph_intmeth,dtset%ph_wstep,dtset%ph_smear,dtset%ph_ngqpt,dtset%ph_qshift)

   !call phdos_print_debye(phdos, cryst%ucvol)
   if (my_rank == master) then
     path = strcat(dtfil%filnam_ds(4), "_PHDOS")
     call wrtout(ab_out, sjoin("- Writing phonon dos to file:", path))
     call phdos_print(phdos, path)
     !call phdos_print_debye(phdos, crystal%ucvol)


!TODO: do we want to pass the temper etc... from anaddb_dtset into the full dtset for abinit? Otherwise just leave these defaults.
     path = strcat(dtfil%filnam_ds(4), "_MSQD_T")
     call phdos_print_msqd(phdos, path, 1000, one, one)

#ifdef HAVE_NETCDF
     path = strcat(dtfil%filnam_ds(4), "_PHDOS.nc")
     ncerr = nctk_open_create(ncid, path, xmpi_comm_self)
     NCF_CHECK_MSG(ncerr, "Creating PHDOS.nc file")
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

 call cwtime(cpu,wall,gflops,"stop")
 write(msg,'(2(a,f8.2))')"eph%ifc: cpu:",cpu,", wall: ",wall
 call wrtout(std_out,msg,"COLL",do_flush=.True.)
 call cwtime(cpu,wall,gflops,"start")

 ! Initialize the object used to read DeltaVscf
 call dvdb_init(dvdb, dvdb_path, comm)
 call dvdb_print(dvdb)
 if (dtset%prtvol > 0) call dvdb_list_perts(dvdb, [-1,-1,-1])

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
 call pspini(dtset,dtfil,ecore,psp_gencond,gsqcutc_eff,gsqcutf_eff,level40,&
&  pawrad,pawtab,psps,cryst%rprimd,comm_mpi=comm)

 ! ====================================================
 ! === This is the real epc stuff once all is ready ===
 ! ====================================================
 if (dtset%eph_task == 1) then
   ! Compute phonon linewidths in metals.
   call eph_phgamma(wfk0_path,dtfil,ngfftc,ngfftf,dtset,cryst,ebands,dvdb,ddk,ifc,&
&   pawfgr,pawang,pawrad,pawtab,psps,mpi_enreg,n0,comm)
 end if
! TODO: decide whether to make several driver functions.
!  before that, however, need to encapsulate most of the functionalities in eph_phgamma
!  otherwise there will be tons of duplicated code

 if (dtset%eph_task == 2) then
   ! Compute electron-phonon matrix elements
   call eph_gkk(wfk0_path,wfq_path,dtfil,ngfftc,ngfftf,dtset,cryst,ebands,ebands_kq,dvdb,ifc,&
   pawfgr,pawang,pawrad,pawtab,psps,mpi_enreg,n0,comm)
 end if

 if (dtset%eph_task == 3) then
   ! Compute phonon self-energy
   call eph_phpi(wfk0_path,wfq_path,dtfil,ngfftc,ngfftf,dtset,cryst,ebands,ebands_kq,dvdb,ifc,&
   pawfgr,pawang,pawrad,pawtab,psps,mpi_enreg,n0,comm)
 end if

 !=====================
 !==== Free memory ====
 !=====================
 call crystal_free(cryst)
 call dvdb_free(dvdb)
 call ddb_free(ddb)
 call ddk_free(ddk)
 call ifc_free(ifc)
 call ebands_free(ebands)
 if (usewfq) call ebands_free(ebands_kq)
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
