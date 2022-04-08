!!****m* ABINIT/m_rttddft_tdks
!! NAME
!!  m_rttddft_tdks
!!
!! FUNCTION
!!  Contains the main object (tdks) to propagate
!!  the time-dependent Kohn-Sham equations in RT-TDDFT
!!
!! COPYRIGHT
!!  Copyright (C) 2021-2022 ABINIT group (FB, MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! PARENTS
!!  m_rttddft_driver
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_rttddft_tdks

 use defs_basis
 use defs_abitypes,      only: MPI_type
 use defs_datatypes,     only: pseudopotential_type, ebands_t
 use defs_wvltypes,      only: wvl_data, nullify_wvl_data

 use libxc_functionals,  only: libxc_functionals_get_hybridparams
 use m_bandfft_kpt,      only: bandfft_kpt, bandfft_kpt_init1, &
                             & bandfft_kpt_destroy_array
 use m_cgprj,            only: ctocprj
 use m_common,           only: setup1
 use m_dtfil,            only: datafiles_type
 use m_dtset,            only: dataset_type
 use m_ebands,           only: ebands_from_dtset, ebands_free, unpack_eneocc
 use m_energies,         only: energies_type, energies_init
 use m_errors,           only: msg_hndl, assert
 use m_extfpmd,          only: extfpmd_type
 use m_fourier_interpol, only: transgrid
 use m_gemm_nonlop,      only: init_gemm_nonlop, destroy_gemm_nonlop
 use m_geometry,         only: fixsym
 use m_hdr,              only: hdr_type, hdr_init
 use m_initylmg,         only: initylmg
 use m_invovl,           only: init_invovl, destroy_invovl
 use m_io_tools,         only: open_file
 use m_inwffil,          only: inwffil
 use m_kg,               only: kpgio, getph, getcut
 use m_mkrho,            only: mkrho
 use m_mpinfo,           only: proc_distrb_cycle
 use m_occ,              only: newocc
 use m_paw_an,           only: paw_an_type, paw_an_init, paw_an_free, &
                             & paw_an_nullify
 use m_pawang,           only: pawang_type
 use m_pawcprj,          only: pawcprj_type,pawcprj_free,pawcprj_alloc, &
                             & pawcprj_getdim
 use m_paw_dmft,         only: init_sc_dmft,destroy_sc_dmft,paw_dmft_type
 use m_pawfgr,           only: pawfgr_type, pawfgr_init, pawfgr_destroy
 use m_pawfgrtab,        only: pawfgrtab_type, pawfgrtab_init, pawfgrtab_free
 use m_paw_init,         only: pawinit,paw_gencond
 use m_paw_ij,           only: paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify
 use m_paw_mkrho,        only: pawmkrho
 use m_paw_nhat,         only: nhatgrid
 use m_paw_occupancies,  only: initrhoij, pawmkrhoij
 use m_pawrad,           only: pawrad_type
 use m_pawrhoij,         only: pawrhoij_type, pawrhoij_copy, pawrhoij_free, &
                             & pawrhoij_alloc, pawrhoij_inquire_dim
 use m_paw_sphharm,      only: setsym_ylm
 use m_pawtab,           only: pawtab_type, pawtab_get_lsize
 use m_paw_tools,        only: chkpawovlp
 use m_pspini,           only: pspini
 use m_spacepar,         only: setsym
 use m_symtk,            only: symmetrize_xred
 use m_wffile,           only: wffile_type, WffClose
 use m_xmpi,             only: xmpi_bcast, xmpi_sum

 implicit none

 private
!!***

! NAME
! tdks_type: Time Dependent Kohn-Sham type
! Object containing the TD KS orbitals and all other 
! important variables required to run RT-TDDFT
!
! SOURCE
 type,public :: tdks_type

  !scalars
   integer                          :: bantot      !total number of bands
   integer                          :: first_step  !start propagation from first_step (for restart)
   integer                          :: mband_cprj  !nb of band per proc (for cprj)
   integer                          :: mcg         !nb of WFs (cg) coeffs
   integer                          :: mcprj       !nb of cprj (projectors applied to WF)
!  integer                          :: my_nspinor  !nb of spinors treated by proc              !FB: Needed?
   integer                          :: nfftf       !nb of FFT grid pts (fine grid)
   integer                          :: nfft        !nb of FFT grid pts (coarse grid)
   integer                          :: nhatgrdim   !dimension of nhatgr array
   integer                          :: ngrvdw      !dimension of grvdw array
   integer                          :: ntime       !max nb of time steps
   integer                          :: tdener_unit !unit nb of the energy file
   integer                          :: tdrestart_unit !unit nb of the restart file
   integer                          :: unpaw       !paw data tmp file unit
   integer                          :: usexcnhat   !use nhat in the computation of the XC term
   real(dp)                         :: boxcut      !boxcut ratio (gcut(box)/gcut(sphere))
   real(dp)                         :: dt          !propagation time step
   real(dp)                         :: ecore       !core energy
   real(dp)                         :: etot        !total energy
   real(dp)                         :: gsqcut      !cut-off on G^2
   real(dp)                         :: ucvol       !primitive cell volume
   real(dp)                         :: zion        !total ionic charge
   logical                          :: gemm_nonlop_use_gemm !use efficient BLAS call
                                                   !for computing  non local potential
   type(energies_type)              :: energies    !contains various energy values
   type(hdr_type)                   :: hdr         !header: contains various info
   type(paw_dmft_type)              :: paw_dmft    !paw_dmft object (unused but
                                                   !required by various routines)
   type(pawfgr_type)                :: pawfgr      !FFT fine grid in PAW sphere
   type(pawang_type),pointer        :: pawang => NULL() !angular grid in PAW sphere
   type(wvl_data)                   :: wvl         !wavelets ojects (unused but
                                                   !required by various routines)
   character(len=fnlen)             :: fname_tdener!Name of the TDENER file
   character(len=fnlen)             :: fname_wfk   !Name of the input WFK file to use
   !arrays
   integer,allocatable              :: atindx(:)   !index table of atom ordered by type
   integer,allocatable              :: atindx1(:)  !nb of the atom for each index in atindx
   integer,allocatable              :: dimcprj(:)  !Contains dimension for cprj array
   integer,allocatable              :: indsym(:,:,:) !atom indexing for symmetries
   integer,allocatable              :: irrzon(:,:,:) !irreducible Brillouin zone
   integer,allocatable              :: kg(:,:)     !red. coord. of G vecs
   integer,allocatable              :: nattyp(:)   !nb of atoms of different types
   integer,allocatable              :: npwarr(:)   !number of PW at each k-point
   integer,allocatable              :: symrec(:,:,:) !sym. operations in recip space
   real(dp)                         :: gprimd(3,3) !prim cell vectors in recip space
   real(dp)                         :: gmet(3,3)   !metric tensor in recip space
   real(dp)                         :: rprimd(3,3) !prim cell vectors in direct space
   real(dp)                         :: rmet(3,3)   !metric tensor in direct space
   real(dp),allocatable             :: cg(:,:)     !WF coefficients in PW basis <k+G|psi_nk>
   real(dp),allocatable             :: cg0(:,:)    !Initial WF coefficients in PW basis <k+G|psi_nk>
   real(dp),allocatable             :: eigen(:)    !eigen-energies
   real(dp),allocatable             :: eigen0(:)   !Initial eigen-energies (at t=0)
   real(dp),allocatable             :: grvdw(:,:)  !Gradient of the total energy coming
                                                   !from VDW dispersion correction             !FB: Needed?
   real(dp),allocatable             :: occ(:)      !occupation numbers
   real(dp),allocatable             :: occ0(:)     !Initial occupation numbers
   real(dp),allocatable             :: nhat(:,:)   !compensation charge density
   real(dp),allocatable             :: nhatgr(:,:,:) !gradient of nhat
   real(dp),allocatable             :: phnons(:,:,:) !For symmetries (nonsymmorphic translation phases)
   real(dp),allocatable             :: ph1d(:,:)   !Structure factor phase: exp(2Pi i G.xred)
                                                   !on coarse grid
   real(dp),allocatable             :: ph1df(:,:)  !Structure factor phase: exp(2Pi i G.xred) for G
                                                   !on fine grid
   real(dp),allocatable             :: rhog(:,:)   !charge density in recip space
   real(dp),allocatable             :: rhor(:,:)   !charge density in direct space
   real(dp),allocatable             :: taug(:,:)   !kin ener density in recip space            !FB: Needed?
   real(dp),allocatable             :: taur(:,:)   !kin ener density in direct space           !FB: Needed?
   real(dp),allocatable             :: vhartr(:)   !Hartree part of the potential
   real(dp),allocatable             :: vpsp(:)     !PSP part of the potential
   real(dp),allocatable             :: vtrial(:,:) !"Trial" potential
   real(dp),allocatable             :: vxc(:,:)    !XC part of the potential
   real(dp),allocatable             :: vxc_hybcomp(:,:) !Hybrid part of the xc potential       !FB:Needed?
   real(dp),allocatable             :: vxctau(:,:,:) !dV_{XC}/dtau (tau = kin. ener density)
                                                   !for mGGAs                                  !FB: Needed?
   real(dp),allocatable             :: xred(:,:,:) !red. coord. of atoms
   real(dp),allocatable             :: xccc3d(:)   !3D core electron density
                                                   !for XC core correction
   real(dp),allocatable             :: xcctau3d(:) !3D core electron kin ener density
                                                   !for XC core correction
   real(dp),allocatable             :: ylm(:,:)    !real spherical harmonics for each k+G
   real(dp),allocatable             :: ylmgr(:,:,:)!real spherical harmonics gradients         !FB: Needed?
   type(pawcprj_type),allocatable   :: cprj(:,:)   !projectors applied on WF <p_lmn|C_nk>
   type(paw_an_type),allocatable    :: paw_an(:)   !various arrays on angular mesh
   type(pawfgrtab_type),allocatable :: pawfgrtab(:) !PAW atomic data on fine grid
   type(paw_ij_type),allocatable    :: paw_ij(:)   !various arrays on partial waves (i,j channels)
   type(pawrad_type),pointer        :: pawrad(:) => NULL() !radial grid in PAW sphere
   type(pawrhoij_type),pointer      :: pawrhoij(:) !operator rho_ij= <psi|p_i><p_j|psi>
   type(pawtab_type),pointer        :: pawtab(:) => NULL() !tabulated PAW atomic data

    contains

    procedure :: init => tdks_init
    procedure :: free => tdks_free

 end type tdks_type
!***

contains

!!****f* m_rttddft_tdks/tdks_init
!!
!! NAME
!!  tdks_init
!!
!! FUNCTION
!!  Initialize the tdks object
!!
!! INPUTS
!!  codvsn = code version
!!  dtfil <type datafiles_type> = infos about file names, file unit numbers
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  pawang <type(pawang_type)> = paw angular mesh and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)> = paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)> = paw tabulated starting data
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!  tdks <class(tdks_type)> = the tdks object to initialize
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!  m_rttddft_driver/rttddft
!!
!! CHILDREN
!!
!! SOURCE
subroutine tdks_init(tdks ,codvsn, dtfil, dtset, mpi_enreg, pawang, pawrad, pawtab, psps)

 implicit none

 !Arguments ------------------------------------
 !scalars
 class(tdks_type),           intent(inout)        :: tdks
 character(len=8),           intent(in)           :: codvsn
 type(datafiles_type),       intent(in)           :: dtfil
 type(dataset_type),         intent(inout)        :: dtset
 type(MPI_type),             intent(inout)        :: mpi_enreg
 type(pawang_type),          intent(inout),target :: pawang
 type(pseudopotential_type), intent(inout)        :: psps
 !arrays
 type(pawrad_type),          intent(inout),target :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),          intent(inout),target :: pawtab(psps%ntypat*psps%usepaw)

 !Local variables-------------------------------
 !scalars
 integer                     :: ierr
 integer                     :: my_natom
 integer                     :: psp_gencond
 real(dp)                    :: ecut_eff
 character(len=500)          :: msg
 type(extfpmd_type),pointer  :: extfpmd => null()
 !arrays
 real(dp),allocatable        :: doccde(:)

! ***********************************************************************

 my_natom=mpi_enreg%my_natom

 !1) Various initializations & checks (MPI, PW, FFT, PSP, Symmetry ...)
 call first_setup(codvsn,dtfil,dtset,ecut_eff,mpi_enreg,pawrad, &
                & pawtab,psps,psp_gencond,tdks)

 !2) Deals with restart
 !FB: @MT Is this the proper way to read a file in abinit..?
 tdks%first_step = 1
 tdks%fname_tdener = dtfil%fnameabo_td_ener
 tdks%fname_wfk = dtfil%fnamewffk
 if (dtset%td_restart > 0) then
   if (mpi_enreg%me == 0) then
      if (open_file('TD_RESTART', msg, newunit=tdks%tdrestart_unit, status='old', form='formatted') /= 0) then
         write(msg,'(a,a,a)') 'Error while trying to open file TD_RESTART needed for restarting the calculation.'
         ABI_ERROR(msg)
      end if 
      read(tdks%tdrestart_unit,*) tdks%first_step
      tdks%first_step = tdks%first_step + 1
      read(tdks%tdrestart_unit,*) tdks%fname_tdener
      read(tdks%tdrestart_unit,*) tdks%fname_wfk
   end if
   !Send to all procs
   call xmpi_bcast(tdks%first_step,0,mpi_enreg%comm_world,ierr)
   call xmpi_bcast(tdks%fname_tdener,0,mpi_enreg%comm_world,ierr)
   call xmpi_bcast(tdks%fname_wfk,0,mpi_enreg%comm_world,ierr)
 else
   if (mpi_enreg%me == 0) then
      if (open_file('TD_RESTART', msg, newunit=tdks%tdrestart_unit, status='unknown', form='formatted') /= 0) then
         write(msg,'(a,a,a)') 'Error while trying to open file TD_RESTART.'
         ABI_ERROR(msg)
      end if 
   end if
 end if

 !3) Reads initial KS orbitals from file (calls inwffil)
 call read_wfk(dtfil,dtset,ecut_eff,mpi_enreg,tdks)
 
 !4) Init occupation numbers
 ABI_MALLOC(tdks%occ0,(dtset%mband*dtset%nkpt*dtset%nsppol))
 tdks%occ0(:)=dtset%occ_orig(:,1)
 !calc occupation number with metallic occupation using the previously read WF
 if (dtset%occopt>=3.and.dtset%occopt<=9) then  ! allowing for occopt 9
   ABI_MALLOC(doccde,(dtset%mband*dtset%nkpt*dtset%nsppol))
   call newocc(doccde,tdks%eigen,tdks%energies%entropy,tdks%energies%e_fermie, &
             & tdks%energies%e_fermih,dtset%ivalence,dtset%spinmagntarget,     &
             & dtset%mband,dtset%nband,dtset%nelect,dtset%ne_qFD,dtset%nh_qFD, &
             & dtset%nkpt,dtset%nspinor,dtset%nsppol,tdks%occ0,dtset%occopt,    &
             & dtset%prtvol,zero,dtset%tphysel,dtset%tsmear,dtset%wtk,extfpmd)
   ABI_FREE(doccde)
 end if

 !FB-TODO: only do this if needed ie. if printing of occupation is requested
 if (dtset%prtocc > 0) then
   ABI_MALLOC(tdks%cg0,(2,tdks%mcg))
   !FB Ouch! Could we avoid this...?
   tdks%cg0(:,:) = tdks%cg(:,:)
   ABI_MALLOC(tdks%occ,(dtset%mband*dtset%nkpt*dtset%nsppol))
   tdks%occ(:) = tdks%occ0(:)
 end if

 !5) Some further initialization (Mainly for PAW)
 call second_setup(dtset,mpi_enreg,pawang,pawrad,pawtab,psps,psp_gencond,tdks)

 !FB: That should be all for now but there were a few more initialization in
 !g_state.F90 in particular related to electric field, might want to check that out
 !once we reach the point of including external electric field

 !Keep some additional stuff in memory within the tdks object
 tdks%unpaw  = dtfil%unpaw
 tdks%dt     = dtset%dtele
 tdks%ntime  = dtset%ntime

 tdks%pawang => pawang
 tdks%pawrad => pawrad
 tdks%pawtab => pawtab

end subroutine tdks_init

!!****f* m_rttddft_tdks/tdks_free
!!
!! NAME
!!  tdks_free
!!
!! FUNCTION
!!  Free all the memory associated with the tdks object
!!
!! INPUTS
!!  tdks <class(tdks_type)> = the tdks object to free
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!  m_rttddft_driver/rttddft
!!
!! CHILDREN
!!
!! SOURCE
subroutine tdks_free(tdks,dtset,mpi_enreg,psps)

 implicit none

 !Arguments ------------------------------------
 !scalars
 class(tdks_type),           intent(inout) :: tdks
 type(dataset_type),         intent(inout) :: dtset
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(pseudopotential_type), intent(inout) :: psps

! ***********************************************************************

   !Destroy hidden save variables
   call bandfft_kpt_destroy_array(bandfft_kpt,mpi_enreg)
   if (psps%usepaw ==1) then 
      call destroy_invovl(dtset%nkpt)
   end if
   if(tdks%gemm_nonlop_use_gemm) then
      call destroy_gemm_nonlop(dtset%nkpt)
   end if

   !Call type destructors
   call destroy_sc_dmft(tdks%paw_dmft)
   call pawfgr_destroy(tdks%pawfgr)
   call tdks%hdr%free()

   !Nullify pointers
   if(associated(tdks%pawang)) tdks%pawang => null()
   if(associated(tdks%pawrad)) tdks%pawrad => null()
   if(associated(tdks%pawtab)) tdks%pawtab => null()
   if(associated(tdks%pawrhoij)) call pawrhoij_free(tdks%pawrhoij)

   !Deallocate allocatables
   if(allocated(tdks%atindx))      ABI_FREE(tdks%atindx)
   if(allocated(tdks%atindx1))     ABI_FREE(tdks%atindx1)
   if(allocated(tdks%cg))          ABI_FREE(tdks%cg)
   if(allocated(tdks%cg0))         ABI_FREE(tdks%cg0)
   if(allocated(tdks%dimcprj))     ABI_FREE(tdks%dimcprj)
   if(allocated(tdks%eigen))       ABI_FREE(tdks%eigen)
   if(allocated(tdks%eigen0))      ABI_FREE(tdks%eigen0)
   if(allocated(tdks%grvdw))       ABI_FREE(tdks%grvdw)
   if(allocated(tdks%indsym))      ABI_FREE(tdks%indsym)
   if(allocated(tdks%irrzon))      ABI_FREE(tdks%irrzon)
   if(allocated(tdks%kg))          ABI_FREE(tdks%kg)
   if(allocated(tdks%nattyp))      ABI_FREE(tdks%nattyp)
   if(allocated(tdks%nhat))        ABI_FREE(tdks%nhat)
   if(allocated(tdks%nhatgr))      ABI_FREE(tdks%nhatgr)
   if(allocated(tdks%npwarr))      ABI_FREE(tdks%npwarr)
   if(allocated(tdks%occ))         ABI_FREE(tdks%occ)
   if(allocated(tdks%occ0))        ABI_FREE(tdks%occ0)
   if(allocated(tdks%ph1d))        ABI_FREE(tdks%ph1d)
   if(allocated(tdks%ph1df))       ABI_FREE(tdks%ph1df)
   if(allocated(tdks%phnons))      ABI_FREE(tdks%phnons)
   if(allocated(tdks%rhog))        ABI_FREE(tdks%rhog)
   if(allocated(tdks%rhor))        ABI_FREE(tdks%rhor)
   if(allocated(tdks%symrec))      ABI_FREE(tdks%symrec)
   if(allocated(tdks%taug))        ABI_FREE(tdks%taug)
   if(allocated(tdks%taur))        ABI_FREE(tdks%taur)
   if(allocated(tdks%vhartr))      ABI_FREE(tdks%vhartr)
   if(allocated(tdks%vpsp))        ABI_FREE(tdks%vpsp)
   if(allocated(tdks%vtrial))      ABI_FREE(tdks%vtrial)
   if(allocated(tdks%vxc))         ABI_FREE(tdks%vxc)
   if(allocated(tdks%vxctau))      ABI_FREE(tdks%vxctau)
   if(allocated(tdks%vxc_hybcomp)) ABI_FREE(tdks%vxc_hybcomp)
   if(allocated(tdks%xred))        ABI_FREE(tdks%xred)
   if(allocated(tdks%xccc3d))      ABI_FREE(tdks%xccc3d)
   if(allocated(tdks%xcctau3d))    ABI_FREE(tdks%xcctau3d)
   if(allocated(tdks%ylm))         ABI_FREE(tdks%ylm)
   if(allocated(tdks%ylmgr))       ABI_FREE(tdks%ylmgr)

   if (allocated(tdks%cprj)) then       
      call pawcprj_free(tdks%cprj)
      ABI_FREE(tdks%cprj)
   end if
   if (allocated(tdks%paw_an)) then
      call paw_an_free(tdks%paw_an)
      ABI_FREE(tdks%paw_an)
   end if
   if(allocated(tdks%pawfgrtab)) then
      call pawfgrtab_free(tdks%pawfgrtab)
      ABI_FREE(tdks%pawfgrtab)
   end if
   if (allocated(tdks%paw_ij)) then
      call paw_ij_free(tdks%paw_ij)
      ABI_FREE(tdks%paw_ij)
   end if

end subroutine tdks_free

!!****f* m_rttddft_tdks/first_setup
!!
!! NAME
!!  first_setup
!!
!! FUNCTION
!!  Intialize many important quantities before running RT-TDDFT
!!  (PW, FFT, PSP, Symmetry etc.)
!!
!! INPUTS
!!  codvsn = code version
!!  dtfil <type datafiles_type> = infos about file names, file unit numbers
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  ecut_eff <real(dp)> = effective PW cutoff energy
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  pawrad(ntypat*usepaw) <type(pawrad_type)> = paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)> = paw tabulated starting data
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!  psp_gencond <integer> = store conditions for generating psp
!!  tdks <type(tdks_type)> = the tdks object to initialize
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!! It is worth to explain THE USE OF FFT GRIDS:
!! ============================================
!! In case of PAW:
!! ---------------
!!    Two FFT grids are used:
!!    - A "coarse" FFT grid (defined by ecut)
!!      for the application of the Hamiltonian on the plane waves basis.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      Hamiltonian, wave-functions, density related to WFs (rhor here), ...
!!      are expressed on this grid.
!!    - A "fine" FFT grid (defined) by ecutdg)
!!      for the computation of the density inside PAW spheres.
!!      It is defined by nfftf, ngfftf, mgfftf, ...
!!      Total density, potentials, ...
!!      are expressed on this grid.
!! In case of norm-conserving:
!! ---------------------------
!!    - Only the usual FFT grid (defined by ecut) is used.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf)
!!      are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! PARENTS
!!  m_rttddft_tdks/rttddft_init
!!
!! CHILDREN
!!
!! SOURCE
subroutine first_setup(codvsn,dtfil,dtset,ecut_eff,mpi_enreg,pawrad,pawtab,psps,psp_gencond,tdks)

 implicit none

 !Arguments ------------------------------------
 !scalars
 character(len=8),           intent(in)    :: codvsn
 integer,                    intent(out)   :: psp_gencond
 real(dp),                   intent(out)   :: ecut_eff
 type(datafiles_type),       intent(in)    :: dtfil
 type(dataset_type),         intent(inout) :: dtset
 type(pseudopotential_type), intent(inout) :: psps
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(tdks_type),            intent(inout) :: tdks
 !arrays
 type(pawrad_type),          intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),          intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

 !Local variables-------------------------------
 !scalars
 integer,parameter    :: response=0, cplex=1
 integer              :: comm_psp
 integer              :: gscase
 integer              :: iatom, ierr, itypat, indx
 integer              :: mgfftf, my_natom
 integer              :: npwmin, nfftot
 integer              :: ylm_option
 real(dp)             :: gsqcut_eff, gsqcutc_eff
 real(dp)             :: ecutdg_eff
 type(ebands_t)       :: bstruct
 !arrays
 character(len=500)   :: msg
 integer, allocatable :: npwarr_(:)
 integer              :: ngfft(18)
 integer              :: ngfftf(18)
 integer              :: npwtot(dtset%nkpt)

! ***********************************************************************

 my_natom=mpi_enreg%my_natom

 !** Init FFT grid(s) sizes (be careful !)
 !See NOTES in the comments at the beginning of this file.
 tdks%nfft = dtset%nfft
 call pawfgr_init(tdks%pawfgr,dtset,mgfftf,tdks%nfftf,ecut_eff,ecutdg_eff, &
                & ngfft,ngfftf)

 !** Init to zero different energies
 call energies_init(tdks%energies)
 tdks%ecore = zero
 tdks%etot = zero

 !** various additional setup mostly related to fft grids and the box (rprimd, metric..)
 call setup1(dtset%acell_orig,tdks%bantot,dtset,ecutdg_eff,ecut_eff,tdks%gmet, &
           & tdks%gprimd,gsqcut_eff,gsqcutc_eff,ngfftf,ngfft,dtset%nkpt,       &
           & dtset%nsppol,response,tdks%rmet,dtset%rprim_orig,tdks%rprimd,     &
           & tdks%ucvol,psps%usepaw)

!FB: @MT Needed?
!!In some cases (e.g. getcell/=0), the plane wave vectors have
!! to be generated from the original simulation cell
!rprimd_for_kg=rprimd
!if (dtset%getcell/=0.and.dtset%usewvl==0) rprimd_for_kg=args_gs%rprimd_orig
!call matr3inv(rprimd_for_kg,gprimd_for_kg)
!gmet_for_kg=matmul(transpose(gprimd_for_kg),gprimd_for_kg)

 !** Set up the basis sphere of planewaves
 ABI_MALLOC(tdks%npwarr,(dtset%nkpt))
 ABI_MALLOC(tdks%kg,(3,dtset%mpw*dtset%mkmem))
 call kpgio(ecut_eff,dtset%exchn2n3d,tdks%gmet,dtset%istwfk,tdks%kg,dtset%kptns, &
          & dtset%mkmem,dtset%nband,dtset%nkpt,'PERS',mpi_enreg,dtset%mpw,       &
          & tdks%npwarr,npwtot,dtset%nsppol)
 call bandfft_kpt_init1(bandfft_kpt,dtset%istwfk,tdks%kg,dtset%mgfft,dtset%mkmem, &
                      & mpi_enreg,dtset%mpw,dtset%nband,dtset%nkpt,tdks%npwarr,   &
                      & dtset%nsppol)

 !** Use efficient BLAS calls for computing the non local potential
 if(dtset%use_gemm_nonlop == 1 .and. dtset%use_gpu_cuda/=1) then
   ! set global variable
   tdks%gemm_nonlop_use_gemm = .true.
   call init_gemm_nonlop(dtset%nkpt)
 else
   tdks%gemm_nonlop_use_gemm = .false.
 end if

 !** Setup the Ylm for each k point
 if (psps%useylm==1) then
   ABI_MALLOC(tdks%ylm,(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm))
   ABI_MALLOC(tdks%ylmgr,(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm))
   ylm_option=0
   if ((dtset%prtstm==0.and.dtset%iscf>0.and.dtset%positron/=1) .or. &
   &   (dtset%berryopt==4 .and. dtset%optstress /= 0 .and. psps%usepaw==1) .or. &
   &   (dtset%orbmag<0 .and. psps%usepaw==1)) then
      ylm_option = 1 ! gradients of YLM
   end if
   call initylmg(tdks%gprimd,tdks%kg,dtset%kptns,dtset%mkmem,mpi_enreg,&
   & psps%mpsang,dtset%mpw,dtset%nband,dtset%nkpt,&
   & tdks%npwarr,dtset%nsppol,ylm_option,tdks%rprimd,tdks%ylm,tdks%ylmgr)
 else
   ABI_MALLOC(tdks%ylm,(0,0))
   ABI_MALLOC(tdks%ylmgr,(0,0,0))
 end if

 !** Open and read pseudopotential files
 comm_psp=mpi_enreg%comm_cell
 call pspini(dtset,dtfil,tdks%ecore,psp_gencond,gsqcutc_eff,gsqcut_eff,pawrad, &
           & pawtab,psps,tdks%rprimd,comm_mpi=comm_psp)

 !In case of isolated computations, ecore must be set to zero
 !because its contribution is counted in the ewald energy as the ion-ion interaction.
 if (dtset%icoulomb == 1) tdks%ecore = zero

 !Include core energy?
 select case(dtset%usepotzero)
 case(0,1)
   tdks%energies%e_corepsp   = tdks%ecore / tdks%ucvol
   tdks%energies%e_corepspdc = zero
 case(2)
   ! No need to include the PspCore energy since it is already included in the
   ! local pseudopotential  (vpsp)
   tdks%energies%e_corepsp   = zero
   tdks%energies%e_corepspdc = zero
 end select

 !** Initialize band structure datatype
 ABI_MALLOC(npwarr_,(dtset%nkpt))
 npwarr_(:)=tdks%npwarr(:)
 if (dtset%paral_kgb/=0) then
   call xmpi_sum(npwarr_,mpi_enreg%comm_bandfft,ierr)
 end if
 bstruct = ebands_from_dtset(dtset, npwarr_)
 ABI_FREE(npwarr_)
 call unpack_eneocc(dtset%nkpt,dtset%nsppol,bstruct%mband,bstruct%nband,dtset%occ_orig(:,1),bstruct%occ,val=zero)

 !** Initialize PAW atomic occupancies
 ABI_MALLOC(tdks%pawrhoij,(my_natom*psps%usepaw))
 if (psps%usepaw == 1) then
   call initrhoij(dtset%pawcpxocc,dtset%lexexch,dtset%lpawu,my_natom,dtset%natom, &
                & dtset%nspden,dtset%nspinor,dtset%nsppol,dtset%ntypat,           &
                & tdks%pawrhoij,dtset%pawspnorb,pawtab,cplex,dtset%spinat,        &
                & dtset%typat,comm_atom=mpi_enreg%comm_atom,                      &
                & mpi_atmtab=mpi_enreg%my_atmtab)
 end if

 !Nullify wvl_data. It is important to do so irregardless of the value of usewvl
 !Only needed here because hdr%init requires a wvl object for the wvl%descr input
 call nullify_wvl_data(tdks%wvl)

 !** Initialize header
 gscase=0
 call hdr_init(bstruct,codvsn,dtset,tdks%hdr,pawtab,gscase,psps,tdks%wvl%descr,&
             & comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

 !Clean band structure datatype
 call ebands_free(bstruct)

 !** PW basis set: test if the problem is ill-defined.
 npwmin=minval(tdks%hdr%npwarr(:))
 if (dtset%mband > npwmin) then
   ! No way we can solve the problem. Abort now!
   write(msg,"(2(a,i0),4a)") "Number of bands nband= ",dtset%mband, &
   & " > number of planewaves npw= ",npwmin,ch10,                   &
   & "The number of eigenvectors cannot be greater that the size of the Hamiltonian!",&
   & ch10, "Action: decrease nband or, alternatively, increase ecut"
   if (dtset%ionmov/=23) then
      ABI_ERROR(msg)
   else
      ABI_WARNING(msg)
   end if

 else if (dtset%mband >= 0.9 * npwmin) then
   ! Warn the user
   write(msg,"(a,i0,a,f6.1,4a)") "Number of bands nband= ",dtset%mband, &
   & " >= 0.9 * maximum number of planewaves= ",0.9*npwmin,ch10,&
   & "This could lead to some instabilities, you might want to decrease nband or increase ecut!", &
   & ch10,"Assume experienced user. Execution will continue."
   ABI_WARNING(msg)
 end if

 !** Initialize symmetry
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 ABI_MALLOC(tdks%irrzon,(nfftot**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
 ABI_MALLOC(tdks%phnons,(2,nfftot**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
 ABI_MALLOC(tdks%indsym,(4,dtset%nsym,dtset%natom))
 ABI_MALLOC(tdks%symrec,(3,3,dtset%nsym))
 tdks%irrzon(:,:,:)=0
 tdks%phnons(:,:,:)=zero
 tdks%indsym(:,:,:)=0
 tdks%symrec(:,:,:)=0

 !TODO FB: Should symmetry be used when ions are moving? Modify if Ehrenfest dynamics
 !Do symmetry stuff if nsym>1
 if (dtset%nsym>1) then
   call setsym(tdks%indsym,tdks%irrzon,dtset%iscf,dtset%natom, &
   & nfftot,ngfft,dtset%nspden,dtset%nsppol,dtset%nsym, &
   & tdks%phnons,dtset%symafm,tdks%symrec,dtset%symrel, &
   & dtset%tnons,dtset%typat,dtset%xred_orig)

   !Make sure dtset%iatfix does not break symmetry
   call fixsym(dtset%iatfix,tdks%indsym,dtset%natom,dtset%nsym)
 else
   !The symrec array is used by initberry even in case nsym = 1 - FB: @MT Needed ?
   tdks%symrec(:,:,1) = 0
   tdks%symrec(1,1,1) = 1 ; tdks%symrec(2,2,1) = 1 ; tdks%symrec(3,3,1) = 1
 end if

 !** Initialize and eventually symmetrize reduced atomic coordinates
 ABI_MALLOC(tdks%xred,(3,dtset%natom,dtset%nimage))
 tdks%xred = dtset%xred_orig
 !FB: Should we?
 !Eventually symmetrize atomic coordinates over space group elements
 call symmetrize_xred(dtset%natom,dtset%nsym,dtset%symrel,dtset%tnons,tdks%xred, &
                    & indsym=tdks%indsym)

 !** Create the atindx array
 !** index table of atoms, in order for them to be used type after type.
 ABI_MALLOC(tdks%atindx,(dtset%natom))
 ABI_MALLOC(tdks%atindx1,(dtset%natom))
 ABI_MALLOC(tdks%nattyp,(psps%ntypat))
 indx=1
 do itypat=1,psps%ntypat
   tdks%nattyp(itypat)=0
   do iatom=1,dtset%natom
      if(dtset%typat(iatom)==itypat)then
         tdks%atindx(iatom)=indx
         tdks%atindx1(indx)=iatom
         indx=indx+1
          tdks%nattyp(itypat)=tdks%nattyp(itypat)+1
      end if
   end do
 end do

 !** Calculate zion: the total positive charge acting on the valence electrons
 tdks%zion=zero
 do iatom=1,dtset%natom
   tdks%zion=tdks%zion+psps%ziontypat(dtset%typat(iatom))
 end do

 !FB: probably not needed
 !Further setup
 !call setup2(dtset,npwtot,start,tdks%wvl%wfs,tdks%xred)

end subroutine first_setup

!!****f* m_rttddft_tdks/second_setup
!!
!! NAME
!! second_setup
!!
!! FUNCTION
!! Further important initialization required after reading WFK and computing
!! occupation numbers in paticular related to PAW
!!
!! INPUTS
!! dtset <type(dataset_type)> = all input variables for this dataset
!! mpi_enreg <MPI_type> = MPI-parallelisation information
!! pawang <type(pawang_type)> = paw angular mesh and related data
!! pawrad(ntypat*usepaw) <type(pawrad_type)> = paw radial mesh and related data
!! pawtab(ntypat*usepaw) <type(pawtab_type)> = paw tabulated starting data
!! psps <type(pseudopotential_type)> = variables related to pseudopotentials
!! psp_gencond <integer> = store conditions for generating psp
!! tdks <type(tdks_type)> = the tdks object to initialize
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!  m_rttddft_tdks/rttddft_init
!!
!! CHILDREN
!!
!! SOURCE
subroutine second_setup(dtset, mpi_enreg, pawang, pawrad, pawtab, psps, psp_gencond, tdks)

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,                    intent(in)    :: psp_gencond
 type(pawang_type),          intent(inout) :: pawang
 type(dataset_type),         intent(inout) :: dtset
 type(pseudopotential_type), intent(inout) :: psps
 type(MPI_type),             intent(inout) :: mpi_enreg
 !arrays
 type(pawrad_type),          intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),          intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
 type(tdks_type),            intent(inout) :: tdks

 !Local variables-------------------------------
 !scalars
 logical             :: call_pawinit
 integer, parameter  :: cplex = 1
 integer             :: forces_needed
 integer             :: gnt_option
 integer             :: has_dijhat, has_vhartree, has_dijfock
 integer             :: has_dijnd, has_dijU, has_vxctau
 integer             :: my_natom, my_nspinor
 integer             :: ncpgr
 integer             :: optcut, optgr0, optgr1, optgr2, optrad
 integer             :: stress_needed
 integer             :: use_hybcomp
 real(dp)            :: gsqcut_shp
 real(dp)            :: hyb_range_fock
 real(dp),parameter  :: k0(3)=(/zero,zero,zero/)
 !arrays
 integer,allocatable :: l_size_atm(:)

! ***********************************************************************

 my_natom=mpi_enreg%my_natom
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

 !FB: @MT needed?
 if (psps%usepaw==1) then
    call pawrhoij_copy(tdks%hdr%pawrhoij,tdks%pawrhoij,comm_atom=mpi_enreg%comm_atom, &
                     & mpi_atmtab=mpi_enreg%my_atmtab)
 end if

 !FB: Needed because paw_dmft is required in mkrho
 !PAW related operations
 !Initialize paw_dmft, even if neither dmft not paw are used
 call init_sc_dmft(dtset%nbandkss,dtset%dmftbandi,dtset%dmftbandf,                &
                 & dtset%dmft_read_occnd,dtset%mband,dtset%nband,dtset%nkpt,      &
                 & dtset%nspden,dtset%nspinor,dtset%nsppol,tdks%occ0,dtset%usedmft,&
                 & tdks%paw_dmft,dtset%usedmft,dtset%dmft_solv,mpi_enreg)

 !*** Main PAW initialization ***
 tdks%mcprj=0;tdks%mband_cprj=0
 if(psps%usepaw==1) then
   gnt_option=1
   if (dtset%pawxcdev==2.or.(dtset%pawxcdev==1.and.dtset%positron/=0)) gnt_option=2

   !** Test if we have to call pawinit
   ! Some gen-cond have to be added...
   call paw_gencond(dtset,gnt_option,"test",call_pawinit)

   if (psp_gencond==1.or.call_pawinit) then
      gsqcut_shp=two*abs(dtset%diecut)*dtset%dilatmx**2/pi**2
      hyb_range_fock=zero
      if (dtset%ixc<0) then
         call libxc_functionals_get_hybridparams(hyb_range=hyb_range_fock)
      end if
      call pawinit(dtset%effmass_free,gnt_option,gsqcut_shp,hyb_range_fock,  &
                 & dtset%pawlcutd,dtset%pawlmix,psps%mpsang,dtset%pawnphi,   &
                 & dtset%nsym,dtset%pawntheta,pawang,pawrad,dtset%pawspnorb, &
                 & pawtab,dtset%pawxcdev,dtset%ixc,dtset%usepotzero)

      ! Update internal values
      call paw_gencond(dtset,gnt_option,"save",call_pawinit)
   end if
   psps%n1xccc=maxval(pawtab(1:psps%ntypat)%usetcore)
   call setsym_ylm(tdks%gprimd,pawang%l_max-1,dtset%nsym,dtset%pawprtvol, &
                 & tdks%rprimd,tdks%symrec,pawang%zarot)

   !** Initialisation of cprj
   tdks%mband_cprj=dtset%mband
   if (dtset%paral_kgb/=0) tdks%mband_cprj=tdks%mband_cprj/mpi_enreg%nproc_band
   tdks%mcprj=my_nspinor*tdks%mband_cprj*dtset%mkmem*dtset%nsppol
   ABI_MALLOC(tdks%cprj,(dtset%natom,tdks%mcprj))
   ncpgr=0
   !FB: @MT dimcprj_srt needed?
   ABI_MALLOC(tdks%dimcprj,(dtset%natom))
   !ABI_MALLOC(dimcprj_srt,(dtset%natom))
   call pawcprj_getdim(tdks%dimcprj,dtset%natom,tdks%nattyp,dtset%ntypat, &
                     & dtset%typat,pawtab,'R')
   !call pawcprj_getdim(dimcprj_srt,dtset%natom,tdks%nattyp,dtset%ntypat,  &
   !                  & dtset%typat,pawtab,'O')
   !call pawcprj_alloc(tdks%cprj,ncpgr,dimcprj_srt)
   call pawcprj_alloc(tdks%cprj,ncpgr,tdks%dimcprj)
   !ABI_FREE(dimcprj_srt)

   !** Variables/arrays related to the fine FFT grid
   ABI_MALLOC(tdks%pawfgrtab,(my_natom))
   if (my_natom>0) then
     call pawtab_get_lsize(pawtab,l_size_atm,my_natom,dtset%typat, &
                         & mpi_atmtab=mpi_enreg%my_atmtab)
     call pawfgrtab_init(tdks%pawfgrtab,cplex,l_size_atm,dtset%nspden,dtset%typat, &
                     & mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
     ABI_FREE(l_size_atm)
   end if
   tdks%usexcnhat=maxval(pawtab(:)%usexcnhat)

   !** Variables/arrays related to the PAW spheres
   ABI_MALLOC(tdks%paw_ij,(my_natom))
   ABI_MALLOC(tdks%paw_an,(my_natom))
   call paw_an_nullify(tdks%paw_an)
   call paw_ij_nullify(tdks%paw_ij)
   has_dijhat=0; if (dtset%iscf==22) has_dijhat=1
   has_vhartree=0; if (dtset%prtvha > 0 .or. dtset%prtvclmb > 0) has_vhartree=1
   has_dijnd=0;if(any(abs(dtset%nucdipmom)>tol8)) has_dijnd=1
   has_dijfock=0
   has_dijU=merge(0,1,dtset%usepawu>0) !Be careful on this!
   has_vxctau=dtset%usekden
   call paw_an_init(tdks%paw_an,dtset%natom,dtset%ntypat,0,0,dtset%nspden,        &
                  & cplex,dtset%pawxcdev,dtset%typat,pawang,pawtab,has_vxc=1,     &
                  & has_vxctau=has_vxctau,has_vxc_ex=1,has_vhartree=has_vhartree, &
                  & comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   call paw_ij_init(tdks%paw_ij,cplex,dtset%nspinor,dtset%nsppol,dtset%nspden,   &
                  & dtset%pawspnorb,dtset%natom,dtset%ntypat,dtset%typat,pawtab, &
                  & has_dij=1,has_dijfock=has_dijfock,has_dijhartree=1,          &
                  & has_dijnd=has_dijnd,has_dijso=1,has_dijhat=has_dijhat,       &
                  & has_dijU=has_dijU,has_pawu_occ=1,has_exexch_pot=1,           &
                  & nucdipmom=dtset%nucdipmom,comm_atom=mpi_enreg%comm_atom,     &
                  & mpi_atmtab=mpi_enreg%my_atmtab)

   !** Check for non-overlapping spheres
   call chkpawovlp(dtset%natom,psps%ntypat,dtset%pawovlp,pawtab,tdks%rmet, &
                 & dtset%typat,tdks%xred)

   !** Identify parts of the rectangular grid where the density has to be calculated
   optcut=0;optgr0=dtset%pawstgylm;optgr1=0;optgr2=0;optrad=1-dtset%pawstgylm
   forces_needed=0 !FB TODO Maybe needs to be changed if Ehrenfest?
   stress_needed=0
   if ((forces_needed==1) .or.                                                &
     & (dtset%xclevel==2 .and. dtset%pawnhatxc>0 .and. tdks%usexcnhat>0) .or. &
     & (dtset%positron/=0.and.forces_needed==2)) then
     optgr1=dtset%pawstgylm
     if (stress_needed==1) optrad=1; if (dtset%pawprtwf==1) optrad=1
   end if
   call nhatgrid(tdks%atindx1,tdks%gmet,my_natom,dtset%natom,                    &
               & tdks%nattyp,tdks%pawfgr%ngfft,psps%ntypat,optcut,optgr0,optgr1, &
               & optgr2,optrad,tdks%pawfgrtab,pawtab,tdks%rprimd,dtset%typat,    &
               & tdks%ucvol,tdks%xred,comm_atom=mpi_enreg%comm_atom,             &
               & mpi_atmtab=mpi_enreg%my_atmtab,comm_fft=mpi_enreg%comm_fft,     &
               & distribfft=mpi_enreg%distribfft)

   tdks%nhatgrdim=0;if (dtset%xclevel==2) tdks%nhatgrdim=tdks%usexcnhat*dtset%pawnhatxc
   if (tdks%nhatgrdim>0)   then
      ABI_MALLOC(tdks%nhatgr,(cplex*tdks%nfftf,dtset%nspden,3*tdks%nhatgrdim))
   else
      ABI_MALLOC(tdks%nhatgr,(0,0,0))
   end if

   ABI_MALLOC(tdks%nhat,(tdks%nfftf,dtset%nspden*psps%usepaw))

   !Required in the PAW case to compute the inverse of the overlap (invovl) operator
   call init_invovl(dtset%nkpt)
 else
   ABI_MALLOC(tdks%nhat,(0,0))
   ABI_MALLOC(tdks%nhatgr,(0,0,0))
   tdks%nhatgrdim=0
 end if

 !Allocate various required arrays for calculation of the Hamiltonian
 !Potentials
 ABI_MALLOC(tdks%vhartr,(tdks%nfftf))
 tdks%vhartr=zero
 ABI_MALLOC(tdks%vpsp,(tdks%nfftf))
 tdks%vpsp=zero
 ABI_MALLOC(tdks%vtrial,(tdks%nfftf,dtset%nspden))
 tdks%vtrial=zero
 ABI_MALLOC(tdks%vxc,(tdks%nfftf,dtset%nspden))
 tdks%vxc=zero
 if (psps%n1xccc/=0) then
    ABI_MALLOC(tdks%xccc3d,(tdks%nfftf))
 else
    ABI_MALLOC(tdks%xccc3d,(0))
 end if
 tdks%xccc3d=zero
 if (psps%usepaw==1) then
    ABI_MALLOC(tdks%xcctau3d,(tdks%nfftf*dtset%usekden))
    tdks%xcctau3d=zero
 endif
 !For mGGA
 ABI_MALLOC(tdks%vxctau,(tdks%nfftf,dtset%nspden,4*dtset%usekden))
 tdks%vxctau=zero
 !For hybrid functionals
 use_hybcomp=0
 if(mod(dtset%fockoptmix,100)==11) use_hybcomp=1
 ABI_MALLOC(tdks%vxc_hybcomp,(tdks%pawfgr%nfft,dtset%nspden*use_hybcomp))
 tdks%vxc_hybcomp=zero
 !For VDW corrected functionals
 tdks%ngrvdw=0
 if ((dtset%vdw_xc>=5.and.dtset%vdw_xc<=7)) then
   tdks%ngrvdw=dtset%natom
 end if
 ABI_MALLOC(tdks%grvdw,(3,tdks%ngrvdw))
 tdks%grvdw=zero

 !Compute large sphere G^2 cut-off (gsqcut) and box / sphere ratio
 if (psps%usepaw==1) then
   call getcut(tdks%boxcut,dtset%pawecutdg,tdks%gmet,tdks%gsqcut,dtset%iboxcut, &
             & std_out,k0,tdks%pawfgr%ngfft)
 else
   call getcut(tdks%boxcut,dtset%ecut,tdks%gmet,tdks%gsqcut,dtset%iboxcut, &
             & std_out,k0,tdks%pawfgr%ngfft)
 end if

 !Compute structure factor phases (exp(2Pi i G.xred)) on coarse and fine grid
 ABI_MALLOC(tdks%ph1d,(2,3*(2*tdks%pawfgr%mgfftc+1)*dtset%natom))
 ABI_MALLOC(tdks%ph1df,(2,3*(2*tdks%pawfgr%mgfft+1)*dtset%natom))
 call getph(tdks%atindx,dtset%natom,tdks%pawfgr%ngfftc(1),tdks%pawfgr%ngfftc(2), &
          & tdks%pawfgr%ngfftc(3),tdks%ph1d,tdks%xred)
 if (psps%usepaw==1.and.tdks%pawfgr%usefinegrid==1) then
   call getph(tdks%atindx,dtset%natom,tdks%pawfgr%ngfft(1),tdks%pawfgr%ngfft(2), &
            & tdks%pawfgr%ngfft(3),tdks%ph1df,tdks%xred)
 else
   tdks%ph1df(:,:)=tdks%ph1d(:,:)
 end if

!!FB: @MT Needed? If yes, then don't forget to put it back in the begining of
!! propagate_ele as well
!!if any nuclear dipoles are nonzero, compute the vector potential in real space (depends on
!!atomic position so should be done for nstep = 1 and for updated ion positions
!if ( any(abs(dtset%nucdipmom(:,:))>tol8) ) then
!   with_vectornd = 1
!else
!   with_vectornd = 0
!end if
!if(allocated(vectornd)) then
!   ABI_FREE(vectornd)
!end if
!ABI_MALLOC(vectornd,(with_vectornd*nfftf,3))
!vectornd=zero
!if(with_vectornd .EQ. 1) then
!   call make_vectornd(1,gsqcut,psps%usepaw,mpi_enreg,dtset%natom,nfftf,ngfftf,dtset%nucdipmom,&
!        & rprimd,vectornd,xred)
!endif

 !Allocate memory for density
 ABI_MALLOC(tdks%rhor,(tdks%nfftf,dtset%nspden))
 ABI_MALLOC(tdks%taur,(tdks%nfftf,dtset%nspden*dtset%usekden))
 ABI_MALLOC(tdks%rhog,(2,tdks%nfftf))
 ABI_MALLOC(tdks%taug,(2,tdks%nfftf*dtset%usekden))

end subroutine second_setup

!!****f* m_rttddft_tdks/read_wfk
!!
!! NAME
!! read_wfk
!!
!! FUNCTION
!! Reads initial wavefunctions (KS orbitals) in WFK file (call inwffil)
!!
!! INPUTS
!! dtfil <type datafiles_type> = infos about file names, file unit numbers
!! dtset <type(dataset_type)> = all input variables for this dataset
!! ecut_eff <real(dp)> = effective PW cutoff energy
!! mpi_enreg <MPI_type> = MPI-parallelisation information
!! tdks <type(tdks_type)> = the tdks object to initialize
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!  m_rttddft_tdks/rttddft_init
!!
!! CHILDREN
!!
!! SOURCE
subroutine read_wfk(dtfil, dtset, ecut_eff, mpi_enreg, tdks)

 implicit none

 !Arguments ------------------------------------
 !scalars
 real(dp),                   intent(in)    :: ecut_eff
 type(datafiles_type),       intent(in)    :: dtfil
 type(dataset_type),         intent(inout) :: dtset
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(tdks_type),            intent(inout) :: tdks

 !Local variables-------------------------------
 !scalars
 integer,parameter           :: formeig=0
 integer                     :: ask_accurate
 integer                     :: band
 integer                     :: cnt
 integer                     :: ierr, ikpt
 integer                     :: my_nspinor
 integer                     :: optorth
 integer                     :: spin
 type(wffile_type)           :: wff1, wffnow
 !arrays
 character(len=500)          :: msg

! ***********************************************************************

 !If paral_kgb == 0, it may happen that some processors are idle (no entry in proc_distrb)
 !but mkmem == nkpt and this can cause integer overflow in mcg or allocation error.
 !Here we count the number of states treated by the proc. if cnt == 0, mcg is then set to 0.
 cnt = 0
 do spin=1,dtset%nsppol
   do ikpt=1,dtset%nkpt
      do band=1,dtset%nband(ikpt + (spin-1) * dtset%nkpt)
         if (.not. proc_distrb_cycle(mpi_enreg%proc_distrb, ikpt, band, band, spin, mpi_enreg%me_kpt)) cnt = cnt + 1
      end do
   end do
 end do

 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 tdks%mcg=dtset%mpw*my_nspinor*dtset%mband*dtset%mkmem*dtset%nsppol
 if (cnt == 0) then
   tdks%mcg = 0
   write(msg,"(2(a,i0))")"rank: ",mpi_enreg%me, "does not have wavefunctions to treat. Setting mcg to: ",tdks%mcg
   ABI_WARNING(msg)
 end if

 if (dtset%usewvl == 0 .and. dtset%mpw > 0 .and. cnt /= 0)then
   if (my_nspinor*dtset%mband*dtset%mkmem*dtset%nsppol > floor(real(HUGE(0))/real(dtset%mpw) )) then
      ierr = 0
      write (msg,'(9a)')&
     & "Default integer is not wide enough to store the size of the wavefunction array (mcg).",ch10,&
     & "This usually happens when paral_kgb == 0 and there are not enough procs to distribute kpts and spins",ch10,&
     & "Action: if paral_kgb == 0, use nprocs = nkpt * nsppol to reduce the memory per node.",ch10,&
     & "If tdks does not solve the problem, use paral_kgb 1 with nprocs > nkpt * nsppol and use npfft/npband/npspinor",ch10,&
     & "to decrease the memory requirements. Consider also OpenMP threads."
      ABI_ERROR_NOSTOP(msg,ierr)
      write (msg,'(5(a,i0), 2a)')&
     & "my_nspinor: ",my_nspinor, ", mpw: ",dtset%mpw, ", mband: ",dtset%mband,&
     & ", mkmem: ",dtset%mkmem, ", nsppol: ",dtset%nsppol,ch10,&
     & 'Note: Compiling with large int (int64) requires a full software stack (MPI/FFTW/BLAS...) compiled in int64 mode'
      ABI_ERROR(msg)
   end if
 end if

 ! Alloc size for wfk and bands
 ABI_MALLOC_OR_DIE(tdks%cg,(2,tdks%mcg),ierr)
 ABI_MALLOC(tdks%eigen,(dtset%mband*dtset%nkpt*dtset%nsppol))
 ABI_MALLOC(tdks%eigen0,(dtset%mband*dtset%nkpt*dtset%nsppol))

 tdks%eigen(:) = zero
 ask_accurate=0

 !Actually read the intial KS orbitals here
 wff1%unwff=dtfil%unwff1
 optorth=0   !No need to orthogonalize the wfk
 tdks%hdr%rprimd=tdks%rprimd
 tdks%cg=0._dp
 call inwffil(ask_accurate,tdks%cg,dtset,dtset%ecut,ecut_eff,tdks%eigen,     &
            & dtset%exchn2n3d,formeig,tdks%hdr,1,dtset%istwfk,tdks%kg,       &
            & dtset%kptns,dtset%localrdwf,dtset%mband,tdks%mcg,dtset%mkmem,  &
            & mpi_enreg,dtset%mpw,dtset%nband,tdks%pawfgr%ngfft,dtset%nkpt,  &
            & tdks%npwarr,dtset%nsppol,dtset%nsym,dtset%occ_orig,optorth,    &
            & dtset%symafm,dtset%symrel,dtset%tnons,dtfil%unkg,wff1,wffnow,  &
            & dtfil%unwff1,tdks%fname_wfk,tdks%wvl)

 !Close wff1
 call WffClose(wff1,ierr)

 !Keep initial eigenvalues in memory
 tdks%eigen0(:) = tdks%eigen(:)

end subroutine read_wfk

end module m_rttddft_tdks
!!***
