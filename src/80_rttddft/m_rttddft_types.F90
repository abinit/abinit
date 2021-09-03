!!****m* ABINIT/m_rttddft_types
!! NAME
!!  m_rttddft_types
!!
!! FUNCTION
!!  Contains the main object (tdks) to propagate
!!  the time-dependent Kohn-Sham equations in RT-TDDFT
!!
!! COPYRIGHT
!!  Copyright (C) 2021 ABINIT group (FB, MT)
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

module m_rttddft_types

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
 use m_ebands,           only: ebands_from_dtset, ebands_free
 use m_energies,         only: energies_type, energies_init
 use m_errors,           only: msg_hndl, assert
 use m_extfpmd,          only: extfpmd_type
 use m_fock,             only: fock_type, fock_init, fock_destroy,    &
                             & fock_ACE_destroy, fock_common_destroy, &
                             & fock_BZ_destroy, fock_update_exc,      &
                             & fock_updatecwaveocc
 use m_fock_getghc,      only: fock2ACE
 use m_fourier_interpol, only: transgrid
 use m_gemm_nonlop,      only: init_gemm_nonlop, destroy_gemm_nonlop
 use m_geometry,         only: fixsym
 use m_hdr,              only: hdr_type, hdr_init
 use m_initylmg,         only: initylmg
 use m_invovl,           only: init_invovl, destroy_invovl
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
 use m_xmpi,             only: xmpi_sum

 implicit none

 private
!!***

! NAME
! tdks_type: Time Dependent Kohn-Sham type
! Object containing the TD KS orbitals and most of the important
! variables and subroutines required to run RT-TDDFT
!
! SOURCE
 type,public :: tdks_type

  !scalars
   integer                          :: bantot      !total number of bands
   integer                          :: mband_cprj  !nb of band per proc (for cprj?)
   integer                          :: mcg         !nb of WFs (cg) coeffs
   integer                          :: mcprj       !nb of cprj (projectors applied to WF)
   integer                          :: my_nspinor  !nb of spinors treated by proc
   integer                          :: nfftf       !nb of FFT grid pts (fine grid)
   integer                          :: nfft        !nb of FFT grid pts (coarse grid)
   integer                          :: nhatgrdim   !dimension of nhatgr array
   integer                          :: ngrvdw      !dimension of grvdw array
   integer                          :: ntime       !max nb of time steps
   integer                          :: usexcnhat   !max nb of time steps
   real(dp)                         :: boxcut      !boxcut ratio (gcut(box)/gcut(sphere))
   real(dp)                         :: dt          !propagation time step
   real(dp)                         :: ecore       !core energy
   real(dp)                         :: gsqcut      !cut-off on G^2
   real(dp)                         :: ucvol       !primitive cell volume
   real(dp)                         :: zion        !total ionic charge
   logical                          :: gemm_nonlop_use_gemm !use efficient BLAS call
                                                   !for computing  non local potential
   type(energies_type)              :: energies    !contains various energy values
   type(fock_type),pointer          :: fock => NULL() !Object for exact Fock exchange in Hybrid functionals  !FB: Needed?
   type(hdr_type)                   :: hdr         !header: contains various info
   type(paw_dmft_type)              :: paw_dmft    !paw_dmft object (unused but
                                                   !required by various routines)
   type(pawfgr_type)                :: pawfgr      !FFT fine grid in PAW sphere
   type(pawang_type),pointer        :: pawang => NULL() !angular grid in PAW sphere
   type(wvl_data)                   :: wvl         !wavelets ojects (unused but
                                                   !required by various routines)
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
   real(dp),allocatable             :: eigen(:)    !eigen-energies
   real(dp),allocatable             :: grvdw(:,:)  !Gradient of the total energy coming
                                                   !from VDW dispersion correction  !FB: Needed?
   real(dp),allocatable             :: occ(:)      !occupation numbers
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
   real(dp),allocatable             :: taur(:,:)   !kin ener density in direct space                !FB: Needed?
   real(dp),allocatable             :: vhartr(:)   !Hartree part of the potential
   real(dp),allocatable             :: vpsp(:)     !PSP part of the potential
   real(dp),allocatable             :: vtrial(:,:) !"Trial" potential
   real(dp),allocatable             :: vxc(:,:)    !XC part of the potential
   real(dp),allocatable             :: vxc_hybcomp(:,:) !Hybrid part of the xc potential
   real(dp),allocatable             :: vxctau(:,:,:) !dV_{XC}/dtau (tau = kin. ener density)
                                                   !for mGGAs    !FB: Needed?
   real(dp),allocatable             :: xred(:,:,:) !red. coord. of atoms
   real(dp),allocatable             :: xccc3d(:)   !3D core electron density
                                                   !for XC core correction
   real(dp),allocatable             :: xcctau3d(:) !3D core electron kin ener density
                                                   !for XC core correction
   real(dp),allocatable             :: ylm(:,:)    !real spherical harmonics for each k+G
   real(dp),allocatable             :: ylmgr(:,:,:)!real spherical harmonics gradients                    !FB: Needed?
   type(pawcprj_type),allocatable   :: cprj(:,:)   !projectors applied on WF <p_lmn|C_nk>
   type(paw_an_type),allocatable    :: paw_an(:)   !various arrays on angular mesh
   type(pawfgrtab_type),allocatable :: pawfgrtab(:) !PAW atomic data on fine grid
   type(paw_ij_type),allocatable    :: paw_ij(:)   !various arrays on partial waves (i,j channels)
   type(pawrad_type),pointer        :: pawrad(:) => NULL()   !radial grid in PAW sphere
   type(pawrhoij_type),pointer      :: pawrhoij(:) !operator rho_ij= <psi|p_i><p_j|psi>
   type(pawtab_type),pointer        :: pawtab(:) => NULL()   !tabulated PAW atomic data

    contains

    procedure :: init => tdks_init
    procedure :: free => tdks_free

 end type tdks_type
!***

contains

!!****f* m_rttddft_types/tdks_init
!!
!! NAME
!!  tdks_init
!!
!! FUNCTION
!!  Initialize the tdks object
!!
!! INPUTS
!!  tdks <class(tdks_type)> = the tdks object to initialize
!!  codvsn = code version
!!  dtfil <type datafiles_type> = infos about file names, file unit numbers
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  pawang <type(pawang_type)> = paw angular mesh and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)> = paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)> = paw tabulated starting data
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
subroutine tdks_init(tdks, codvsn, dtfil, dtset, mpi_enreg, pawang, pawrad, pawtab, psps)

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
 integer                     :: my_natom
 integer                     :: psp_gencond
 real(dp)                    :: entropy, ecut_eff
 type(extfpmd_type),pointer  :: extfpmd => null()
 !arrays
 real(dp),allocatable        :: doccde(:)

! ***********************************************************************

 my_natom=mpi_enreg%my_natom

 !1) Various initializations & checks (MPI, PW, FFT, PSP, Symmetry ...)
 call first_setup(tdks,codvsn,dtfil,dtset,ecut_eff,mpi_enreg,pawrad, &
                & pawtab,psps,psp_gencond)

 !2) Reads initial KS orbitals from file (calls inwffil)
 call read_wfk(tdks,dtfil,dtset,ecut_eff,mpi_enreg)

 !3) Compute occupation numbers from WF
 ABI_MALLOC(tdks%occ,(dtset%mband*dtset%nkpt*dtset%nsppol))
 ABI_MALLOC(doccde,(dtset%mband*dtset%nkpt*dtset%nsppol))
 call newocc(doccde,tdks%eigen,entropy,tdks%energies%e_fermie,               &
           & tdks%energies%e_fermih,dtset%ivalence,dtset%spinmagntarget,     &
           & dtset%mband,dtset%nband,dtset%nelect,dtset%ne_qFD,dtset%nh_qFD, &
           & dtset%nkpt,dtset%nspinor,dtset%nsppol,tdks%occ,dtset%occopt,    &
           & dtset%prtvol,zero,dtset%tphysel,dtset%tsmear,dtset%wtk,extfpmd)
 ABI_FREE(doccde)

 !4) Some further initialization (Mainly for PAW)
 call second_setup(tdks,dtset,dtfil,mpi_enreg,pawang,pawrad,pawtab,psps,psp_gencond)

 !5) Compute charge density from WFs
 call calc_density(tdks,dtfil,dtset,mpi_enreg,pawang,pawtab,psps)

 !TODO FB: That should be all for now but there were a few more initialization in
 !g_state.F90 in particular related to electric field, might want to check that out
 !once we reach the point of including external electric field

 !Keep some additional stuff in memory within the tdks object
 tdks%dt     = dtset%dtele
 tdks%ntime  = dtset%ntime

 tdks%pawang => pawang
 tdks%pawrad => pawrad
 tdks%pawtab => pawtab

end subroutine tdks_init

!!****f* m_rttddft_types/tdks_free
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
subroutine tdks_free(tdks,dtset,mpi_enreg)

 implicit none

 !Arguments ------------------------------------
 !scalars
 class(tdks_type),   intent(inout) :: tdks
 type(dataset_type), intent(inout) :: dtset
 type(MPI_type),     intent(inout) :: mpi_enreg

! ***********************************************************************

   !Destroy hidden save variables
   call bandfft_kpt_destroy_array(bandfft_kpt,mpi_enreg)
   call destroy_invovl(dtset%nkpt)
   if(tdks%gemm_nonlop_use_gemm) then
      call destroy_gemm_nonlop(dtset%nkpt)
   end if

   !Call type destructors
   call destroy_sc_dmft(tdks%paw_dmft)
   ! Deallocate exact exchange data at the end of the calculation
   if (dtset%usefock==1) then
      if (tdks%fock%fock_common%use_ACE/=0) call fock_ACE_destroy(tdks%fock%fockACE)
      call fock_common_destroy(tdks%fock%fock_common)
      call fock_BZ_destroy(tdks%fock%fock_BZ)
      call fock_destroy(tdks%fock)
      nullify(tdks%fock)
   end if
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
   if(allocated(tdks%dimcprj))     ABI_FREE(tdks%dimcprj)
   if(allocated(tdks%eigen))       ABI_FREE(tdks%eigen)
   if(allocated(tdks%grvdw))       ABI_FREE(tdks%grvdw)
   if(allocated(tdks%indsym))      ABI_FREE(tdks%indsym)
   if(allocated(tdks%irrzon))      ABI_FREE(tdks%irrzon)
   if(allocated(tdks%kg))          ABI_FREE(tdks%kg)
   if(allocated(tdks%nattyp))      ABI_FREE(tdks%nattyp)
   if(allocated(tdks%nhat))        ABI_FREE(tdks%nhat)
   if(allocated(tdks%nhatgr))      ABI_FREE(tdks%nhatgr)
   if(allocated(tdks%npwarr))      ABI_FREE(tdks%npwarr)
   if(allocated(tdks%occ))         ABI_FREE(tdks%occ)
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

!!****f* m_rttddft_types/first_setup
!!
!! NAME
!!  first_setup
!!
!! FUNCTION
!!  Intialize many important quantities before running RT-TDDFT
!!  in particular PW, FFT, PSP, Symmetry etc.
!!
!! INPUTS
!!  tdks <type(tdks_type)> = the tdks object to initialize
!!  codvsn = code version
!!  dtfil <type datafiles_type> = infos about file names, file unit numbers
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  ecut_eff <real(dp)> = effective PW cutoff energy
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  pawrad(ntypat*usepaw) <type(pawrad_type)> = paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)> = paw tabulated starting data
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!  psp_gencond <integer> = store conditions for generating psp
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!  m_rttddft_types/rttddft_init
!!
!! CHILDREN
!!
!! SOURCE
subroutine first_setup(tdks,codvsn,dtfil,dtset,ecut_eff,mpi_enreg,pawrad,pawtab,psps,psp_gencond)

 implicit none

 !Arguments ------------------------------------
 !scalars
 type(tdks_type),            intent(inout) :: tdks
 character(len=8),           intent(in)    :: codvsn
 integer,                    intent(out)   :: psp_gencond
 real(dp),                   intent(out)   :: ecut_eff
 type(datafiles_type),       intent(in)    :: dtfil
 type(dataset_type),         intent(inout) :: dtset
 type(pseudopotential_type), intent(inout) :: psps
 type(MPI_type),             intent(inout) :: mpi_enreg
 !arrays
 type(pawrad_type),          intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),          intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

 !Local variables-------------------------------
 !scalars
 integer,parameter   :: response=0, cplex=1
 integer             :: comm_psp
 integer             :: gscase
 integer             :: ierr, iatom, itypat, indx
 integer             :: mgfftf, my_natom
 integer             :: npwmin, nfftot
 integer             :: ylm_option
 real(dp)            :: gsqcut_eff, gsqcutc_eff
 real(dp)            :: ecutdg_eff
 type(ebands_t)      :: bstruct
 !arrays
 character(len=500)  :: msg
 integer             :: ngfft(18)
 integer             :: ngfftf(18)
 integer             :: npwtot(dtset%nkpt)

! ***********************************************************************

 !Init MPI data
 my_natom=mpi_enreg%my_natom

 !** Init FFT grid(s) sizes (be careful !)
 !See NOTES in the comments at the beginning of this file.
 tdks%nfft = dtset%nfft
 call pawfgr_init(tdks%pawfgr,dtset,mgfftf,tdks%nfftf,ecut_eff,ecutdg_eff, &
                & ngfft,ngfftf)

 !** Init to zero different energies
 call energies_init(tdks%energies)
 tdks%ecore = zero

 !** various additional setup
 call setup1(dtset%acell_orig,tdks%bantot,dtset,ecutdg_eff,ecut_eff,tdks%gmet, &
           & tdks%gprimd,gsqcut_eff,gsqcutc_eff,ngfftf,ngfft,dtset%nkpt,       &
           & dtset%nsppol,response,tdks%rmet,dtset%rprim_orig,tdks%rprimd,     &
           & tdks%ucvol,psps%usepaw)

!FB: Needed?
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

 !** Setup to use efficient BLAS calls for computing the non local potential
 if(dtset%use_gemm_nonlop == 1 .and. dtset%use_gpu_cuda/=1) then
   ! set global variable
   tdks%gemm_nonlop_use_gemm = .true.
   call init_gemm_nonlop(dtset%nkpt)
 else
   tdks%gemm_nonlop_use_gemm = .false.
 end if

 !** Set up the Ylm for each k point
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

 !In case of isolated computations, ecore must set to zero
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

 !FB: The npwarr_ pointer array doesn't seem necessary?
 !** Initialize band structure datatype
 !if (dtset%paral_kgb/=0) then     !  We decide to store total npw in bstruct,
 !  ABI_MALLOC(npwarr_,(dtset%nkpt))
 !  npwarr_(:)=tdks%npwarr(:)
 !  call xmpi_sum(npwarr_,mpi_enreg%comm_bandfft,ierr)
 !else
 !  npwarr_ => npwarr
 !end if
 !bstruct = ebands_from_dtset(dtset, npwarr_)
 !if (dtset%paral_kgb/=0)  then
 !  ABI_FREE(npwarr_)
 !end if
 !nullify(npwarr_)
 if (dtset%paral_kgb/=0) then     !  We decide to store total npw in bstruct,
   call xmpi_sum(tdks%npwarr,mpi_enreg%comm_bandfft,ierr)
 end if
 bstruct = ebands_from_dtset(dtset, tdks%npwarr)

 !** Initialize PAW atomic occupancies
 ABI_MALLOC(tdks%pawrhoij,(my_natom*psps%usepaw))
 if (psps%usepaw == 1) then
   call initrhoij(dtset%pawcpxocc,dtset%lexexch,dtset%lpawu,my_natom,dtset%natom, &
                & dtset%nspden,dtset%nspinor,dtset%nsppol,dtset%ntypat,           &
                & tdks%pawrhoij,dtset%pawspnorb,pawtab,cplex,dtset%spinat,       &
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

 !Clean band structure datatype (should use it more in the future !)
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

 !TODO FB: Symmetry should not be used when ions are moving. Modify for Ehrenfest
 !Do symmetry stuff if nsym>1
 if (dtset%nsym>1) then
   call setsym(tdks%indsym,tdks%irrzon,dtset%iscf,dtset%natom, &
   & nfftot,ngfft,dtset%nspden,dtset%nsppol,dtset%nsym, &
   & tdks%phnons,dtset%symafm,tdks%symrec,dtset%symrel, &
   & dtset%tnons,dtset%typat,dtset%xred_orig)

 !Make sure dtset%iatfix does not break symmetry
   call fixsym(dtset%iatfix,tdks%indsym,dtset%natom,dtset%nsym)
 else
   !The symrec array is used by initberry even in case nsym = 1 - FB: Needed ?
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

 !Update header, with evolving variables, when available
 !Here, rprimd, xred and occ are potentially available
 call tdks%hdr%update(tdks%bantot,tdks%hdr%etot,tdks%hdr%fermie,tdks%hdr%fermih,    &
                    & tdks%hdr%residm,tdks%rprimd,dtset%occ_orig,tdks%pawrhoij,     &
                    & dtset%xred_orig,dtset%amu_orig,comm_atom=mpi_enreg%comm_atom, &
                    & mpi_atmtab=mpi_enreg%my_atmtab)

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

!!****f* m_rttddft_types/second_setup
!!
!! NAME
!! second_setup
!!
!! FUNCTION
!! Further important initialization required after reading WFK and computing
!! occupation numbers in paticular related to PAW
!!
!! INPUTS
!! tdks <type(tdks_type)> = the tdks object to initialize
!! dtset <type(dataset_type)> = all input variables for this dataset
!! mpi_enreg <MPI_type> = MPI-parallelisation information
!! pawang <type(pawang_type)> = paw angular mesh and related data
!! pawrad(ntypat*usepaw) <type(pawrad_type)> = paw radial mesh and related data
!! pawtab(ntypat*usepaw) <type(pawtab_type)> = paw tabulated starting data
!! psps <type(pseudopotential_type)> = variables related to pseudopotentials
!! psp_gencond <integer> = store conditions for generating psp
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!  m_rttddft_types/rttddft_init
!!
!! CHILDREN
!!
!! SOURCE
subroutine second_setup(tdks, dtset, dtfil, mpi_enreg, pawang, pawrad, pawtab, psps, psp_gencond)

 implicit none

 !Arguments ------------------------------------
 !scalars
 type(tdks_type),           intent(inout) :: tdks
 integer,                    intent(in)    :: psp_gencond
 type(pawang_type),          intent(inout) :: pawang
 type(dataset_type),         intent(inout) :: dtset
 type(datafiles_type),       intent(in)    :: dtfil
 type(pseudopotential_type), intent(inout) :: psps
 type(MPI_type),             intent(inout) :: mpi_enreg
 !arrays
 type(pawrad_type),          intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),          intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

 !Local variables-------------------------------
 !scalars
 logical             :: call_pawinit
 integer, parameter  :: cplex = 1
 integer             :: cplex_hf
 integer             :: ctocprj_choice
 integer             :: forces_needed
 integer             :: gnt_option
 integer             :: has_dijhat, has_vhartree, has_dijfock
 integer             :: has_dijnd, has_dijU, has_vxctau
 integer             :: hyb_mixing, hyb_mixing_sr
 integer             :: iatom, idir
 integer             :: iorder_cprj
 integer             :: my_natom
 integer             :: ncpgr
 integer             :: optcut, optgr0, optgr1, optgr2, optrad
 integer             :: stress_needed
 integer             :: use_hybcomp, usecprj
 real(dp)            :: gsqcut_shp
 real(dp)            :: hyb_range_fock
 real(dp),parameter  :: k0(3)=(/zero,zero,zero/)
 !arrays
 integer,allocatable :: dimcprj_srt(:)
 integer,allocatable :: l_size_atm(:)

! ***********************************************************************

 my_natom=mpi_enreg%my_natom

 if (psps%usepaw==1) then
    call pawrhoij_copy(tdks%hdr%pawrhoij,tdks%pawrhoij,comm_atom=mpi_enreg%comm_atom, &
                     & mpi_atmtab=mpi_enreg%my_atmtab)
 end if

 !FB: Needed because paw_dmft is needed in mkrho
 !PAW related operations
 !Initialize paw_dmft, even if neither dmft not paw are used
 call init_sc_dmft(dtset%nbandkss,dtset%dmftbandi,dtset%dmftbandf,                &
                 & dtset%dmft_read_occnd,dtset%mband,dtset%nband,dtset%nkpt,      &
                 & dtset%nspden,dtset%nspinor,dtset%nsppol,tdks%occ,dtset%usedmft,&
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
   tdks%mcprj=tdks%my_nspinor*tdks%mband_cprj*dtset%mkmem*dtset%nsppol
   ABI_MALLOC(tdks%cprj,(dtset%natom,tdks%mcprj))
   ncpgr=0
   if (dtset%usefock==1) then
      if (dtset%optforces == 1) then
         ncpgr = 3
      end if
   end if
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
   !has_dijhat=0; if (dtset%iscf==22) has_dijhat=1
   has_dijhat=1
   has_vhartree=0; if (dtset%prtvha > 0 .or. dtset%prtvclmb > 0) has_vhartree=1
   has_dijfock=0; if (dtset%usefock==1) has_dijfock=1
   has_dijnd=0;if(any(abs(dtset%nucdipmom)>tol8)) has_dijnd=1
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
   forces_needed=0 !FB TODO needs to be changed if Ehrenfest?
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
 end if

 !Allocate various required arrays for calculation of the Hamiltonian
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
 if(mod(dtset%fockoptmix,100)==11)use_hybcomp=1
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

!!FB: Needed? If yes, then don't forget to put it back in the begining of
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

 !** Initialize data in the case of an Exact-exchange (Hartree-Fock) or hybrid XC calculation
 hyb_mixing=zero;hyb_mixing_sr=zero
 if (dtset%usefock==1) then
   ! Initialize data_type fock for the calculation
   cplex_hf=cplex
   if (psps%usepaw==1) cplex_hf=dtset%pawcpxocc
   call fock_init(tdks%atindx,cplex_hf,dtset,tdks%fock,tdks%gsqcut,tdks%kg,    &
                & mpi_enreg,tdks%nattyp,tdks%npwarr,pawang,tdks%pawfgr,pawtab, &
                & tdks%rprimd)
   if (tdks%fock%fock_common%usepaw==1) then
      optcut = 0 ! use rpaw to construct local_pawfgrtab
      optgr0 = 0; optgr1 = 0; optgr2 = 0 ! dont need gY terms locally
      optrad = 1 ! do store r-R
      call nhatgrid(tdks%atindx1,tdks%gmet,dtset%natom,dtset%natom,tdks%nattyp, &
                  & tdks%pawfgr%ngfft,psps%ntypat,optcut,optgr0,optgr1,optgr2,  &
                  & optrad,tdks%fock%fock_common%pawfgrtab,pawtab,tdks%rprimd,  &
                  & dtset%typat,tdks%ucvol,tdks%xred,typord=1)
      iatom=-1;idir=0
      if (dtset%optforces == 1) then
         ctocprj_choice = 2
      else
         ctocprj_choice = 1
      end if
      iorder_cprj = 0 ! cprj are sorted by atom type
      call ctocprj(tdks%atindx,tdks%cg,ctocprj_choice,tdks%cprj,tdks%gmet,  &
                 & tdks%gprimd,iatom,idir,iorder_cprj,dtset%istwfk,tdks%kg, &
                 & dtset%kptns,tdks%mcg,tdks%mcprj,dtset%mgfft,dtset%mkmem, &
                 & mpi_enreg,psps%mpsang,dtset%mpw,dtset%natom,tdks%nattyp, &
                 & dtset%nband,dtset%natom,tdks%pawfgr%ngfftc,dtset%nkpt,   &
                 & dtset%nloalg,tdks%npwarr,dtset%nspinor,dtset%nsppol,     &
                 & dtset%ntypat,dtset%paral_kgb,tdks%ph1d,psps,tdks%rmet,   &
                 & dtset%typat,tdks%ucvol,dtfil%unpaw,tdks%xred,tdks%ylm,   &
                 & tdks%ylmgr)
   end if

   !Fock energy & forces
   tdks%energies%e_exactX=zero
   if (tdks%fock%fock_common%optfor) tdks%fock%fock_common%forces=zero

   !Update data relative to the occupied states in fock
   call fock_updatecwaveocc(tdks%cg,tdks%cprj,dtset,tdks%fock,tdks%indsym,         &
                          & tdks%mcg,tdks%mcprj,mpi_enreg,tdks%nattyp,tdks%npwarr, &
                          & tdks%occ,tdks%ucvol)
   !Possibly (re)compute the ACE operator
   if(tdks%fock%fock_common%use_ACE/=0) then
      if (tdks%mcprj>0) then
         usecprj=1
      else
         usecprj=0
      end if
      call fock2ACE(tdks%cg,tdks%cprj,tdks%fock,dtset%istwfk,tdks%kg,dtset%kptns, &
                  & dtset%mband,tdks%mcg,tdks%mcprj,dtset%mgfft,dtset%mkmem,      &
                  & mpi_enreg,psps%mpsang,dtset%mpw,my_natom,dtset%natom,         &
                  & dtset%nband,dtset%nfft,tdks%pawfgr%ngfftc,dtset%nkpt,         &
                  & dtset%nloalg,tdks%npwarr,dtset%nspden,dtset%nspinor,          &
                  & dtset%nsppol,dtset%ntypat,tdks%occ,dtset%optforces,           &
                  & tdks%paw_ij,pawtab,tdks%ph1d,psps,tdks%rprimd,dtset%typat,    &
                  & usecprj,dtset%use_gpu_cuda,dtset%wtk,tdks%xred,tdks%ylm)
      tdks%energies%e_fock0=tdks%fock%fock_common%e_fock0
   end if

 end if

 !FB: Needed to compute the inverse of the overlap (invovl) operator
 !FB: This seems to init an involv_kpt object declared in the invovl module with the
 !save argument which does not sounds great..
 call init_invovl(dtset%nkpt)

end subroutine second_setup

!!****f* m_rttddft_types/read_wfk
!!
!! NAME
!! read_wfk
!!
!! FUNCTION
!! Reads initial wavefunctions (KS orbitals) in WFK file (call inwfill)
!!
!! INPUTS
!! tdks <type(tdks_type)> = the tdks object to initialize
!! dtfil <type datafiles_type> = infos about file names, file unit numbers
!! dtset <type(dataset_type)> = all input variables for this dataset
!! ecut_eff <real(dp)> = effective PW cutoff energy
!! mpi_enreg <MPI_type> = MPI-parallelisation information
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!  m_rttddft_types/rttddft_init
!!
!! CHILDREN
!!
!! SOURCE
subroutine read_wfk(tdks, dtfil, dtset, ecut_eff, mpi_enreg)

 implicit none

 !Arguments ------------------------------------
 !scalars
 type(tdks_type),           intent(inout) :: tdks
 real(dp),                   intent(in)    :: ecut_eff
 type(datafiles_type),       intent(in)    :: dtfil
 type(dataset_type),         intent(inout) :: dtset
 type(MPI_type),             intent(inout) :: mpi_enreg

 !Local variables-------------------------------
 !scalars
 integer,parameter           :: formeig=0
 integer                     :: ask_accurate
 integer                     :: band
 integer                     :: cnt
 integer                     :: ierr, ikpt
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

 tdks%my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 tdks%mcg=dtset%mpw*tdks%my_nspinor*dtset%mband*dtset%mkmem*dtset%nsppol
 if (cnt == 0) then
   tdks%mcg = 0
   write(msg,"(2(a,i0))")"rank: ",mpi_enreg%me, "does not have wavefunctions to treat. Setting mcg to: ",tdks%mcg
   ABI_WARNING(msg)
 end if

 if (dtset%usewvl == 0 .and. dtset%mpw > 0 .and. cnt /= 0)then
   if (tdks%my_nspinor*dtset%mband*dtset%mkmem*dtset%nsppol > floor(real(HUGE(0))/real(dtset%mpw) )) then
      ierr = 0
      write (msg,'(9a)')&
     & "Default integer is not wide enough to store the size of the wavefunction array (mcg).",ch10,&
     & "This usually happens when paral_kgb == 0 and there are not enough procs to distribute kpts and spins",ch10,&
     & "Action: if paral_kgb == 0, use nprocs = nkpt * nsppol to reduce the memory per node.",ch10,&
     & "If tdks does not solve the problem, use paral_kgb 1 with nprocs > nkpt * nsppol and use npfft/npband/npspinor",ch10,&
     & "to decrease the memory requirements. Consider also OpenMP threads."
      ABI_ERROR_NOSTOP(msg,ierr)
      write (msg,'(5(a,i0), 2a)')&
     & "my_nspinor: ",tdks%my_nspinor, ", mpw: ",dtset%mpw, ", mband: ",dtset%mband,&
     & ", mkmem: ",dtset%mkmem, ", nsppol: ",dtset%nsppol,ch10,&
     & 'Note: Compiling with large int (int64) requires a full software stack (MPI/FFTW/BLAS...) compiled in int64 mode'
      ABI_ERROR(msg)
   end if
 end if

 ! Alloc size for wfk and bands
 ABI_MALLOC_OR_DIE(tdks%cg,(2,tdks%mcg),ierr)
 ABI_MALLOC(tdks%eigen,(dtset%mband*dtset%nkpt*dtset%nsppol))

 tdks%eigen(:) = zero
 ask_accurate=0

 !Actually read the intial KS orbitals here
 wff1%unwff=dtfil%unwff1
 optorth=0   !No need to orthogonalize the wfk
 tdks%hdr%rprimd=tdks%rprimd
 call inwffil(ask_accurate,tdks%cg,dtset,dtset%ecut,ecut_eff,tdks%eigen,     &
            & dtset%exchn2n3d,formeig,tdks%hdr,dtfil%ireadwf,dtset%istwfk,   &
            & tdks%kg,dtset%kptns,dtset%localrdwf,dtset%mband,tdks%mcg,      &
            & dtset%mkmem,mpi_enreg,dtset%mpw,dtset%nband,tdks%pawfgr%ngfft, &
            & dtset%nkpt,tdks%npwarr,dtset%nsppol,dtset%nsym,dtset%occ_orig, &
            & optorth,dtset%symafm,dtset%symrel,dtset%tnons,dtfil%unkg,wff1, &
            & wffnow,dtfil%unwff1,dtfil%fnamewffk,tdks%wvl)

 !Close wff1
 call WffClose(wff1,ierr)

end subroutine read_wfk

!!****f* m_rttddft_types/calc_density
!!
!! NAME
!!  calc_density
!!
!! FUNCTION
!!  Compute electronic density (in 1/bohr^3) from the WF (cg coefficients)
!!
!! INPUTS
!!  tdks <type(tdks_type)> = the tdks object to initialize
!!  dtfil <type datafiles_type> = infos about file names, file unit numbers
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!  m_rttddft_types/rttddft_init
!!
!! CHILDREN
!!
!! SOURCE
subroutine calc_density(tdks, dtfil, dtset, mpi_enreg, pawang, pawtab, psps)

 implicit none

 !Arguments ------------------------------------
 !scalars
 type(tdks_type),            intent(inout) :: tdks
 type(datafiles_type),       intent(in)    :: dtfil
 type(dataset_type),         intent(inout) :: dtset
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(pawang_type),          intent(inout) :: pawang
 type(pseudopotential_type), intent(inout) :: psps
 !arrays
 type(pawtab_type),          intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

 !Local variables-------------------------------
 !scalars
 integer, parameter          :: cplex=1
 integer                     :: cplex_rhoij
 integer                     :: ipert, idir, ider, izero
 integer                     :: my_natom
 integer                     :: nspden_rhoij
 integer                     :: tim_mkrho
 real(dp)                    :: compch_fft
 !arrays
 real(dp)                    :: qpt(3)
 real(dp),allocatable        :: rhowfg(:,:), rhowfr(:,:)
 type(pawrhoij_type),pointer :: pawrhoij_unsym(:)

! ***********************************************************************

 my_natom=mpi_enreg%my_natom

 ABI_MALLOC(tdks%rhor,(tdks%nfftf,dtset%nspden))
 ABI_MALLOC(tdks%taur,(tdks%nfftf,dtset%nspden*dtset%usekden))
 ABI_MALLOC(tdks%rhog,(2,tdks%nfftf))
 ABI_MALLOC(tdks%taug,(2,tdks%nfftf*dtset%usekden))

 tim_mkrho=1

 if (psps%usepaw==1) then

   ABI_MALLOC(rhowfg,(2,dtset%nfft))
   ABI_MALLOC(rhowfr,(dtset%nfft,dtset%nspden))
   ABI_MALLOC(tdks%nhat,(tdks%nfftf,dtset%nspden*psps%usepaw))

   tdks%nhatgrdim=0;if (dtset%xclevel==2) tdks%nhatgrdim=tdks%usexcnhat*dtset%pawnhatxc
   ider=2*tdks%nhatgrdim;izero=0
   if (tdks%nhatgrdim>0)   then
      ABI_MALLOC(tdks%nhatgr,(cplex*tdks%nfftf,dtset%nspden,3*tdks%nhatgrdim))
   else
      ABI_MALLOC(tdks%nhatgr,(0,0,0))
   end if

   ! 1-Compute density from WFs (without compensation charge density nhat)
   call mkrho(tdks%cg,dtset,tdks%gprimd,tdks%irrzon,tdks%kg,tdks%mcg,mpi_enreg, &
            & tdks%npwarr,tdks%occ,tdks%paw_dmft,tdks%phnons,rhowfg,rhowfr,     &
            & tdks%rprimd,tim_mkrho,tdks%ucvol,tdks%wvl%den,tdks%wvl%wfs)
   ! transfer density from the coarse to the fine FFT grid
   call transgrid(1,mpi_enreg,dtset%nspden,+1,1,1,dtset%paral_kgb,tdks%pawfgr, &
                & rhowfg,tdks%rhog,rhowfr,tdks%rhor)

   write(99,*) tdks%rhor

   ! 2-Compute cprj = <\psi_{n,k}|p_{i,j}>
   call ctocprj(tdks%atindx,tdks%cg,1,tdks%cprj,tdks%gmet,tdks%gprimd,0,0,0,      &
              & dtset%istwfk,tdks%kg,dtset%kptns,tdks%mcg,tdks%mcprj,dtset%mgfft, &
              & dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,dtset%natom,          &
              & tdks%nattyp,dtset%nband,dtset%natom,dtset%ngfft,dtset%nkpt,       &
              & dtset%nloalg,tdks%npwarr,dtset%nspinor,dtset%nsppol,psps%ntypat,  &
              & dtset%paral_kgb,tdks%ph1d,psps,tdks%rmet,dtset%typat,tdks%ucvol,  &
              & dtfil%unpaw,tdks%xred,tdks%ylm,tdks%ylmgr)

   !paral atom
   if (my_natom/=dtset%natom) then
      ABI_MALLOC(pawrhoij_unsym,(dtset%natom))
      call pawrhoij_inquire_dim(cplex_rhoij=cplex_rhoij,nspden_rhoij=nspden_rhoij, &
                              & nspden=dtset%nspden,spnorb=dtset%pawspnorb,        &
                              & cpxocc=dtset%pawcpxocc)
      call pawrhoij_alloc(pawrhoij_unsym,cplex_rhoij,nspden_rhoij,dtset%nspinor, &
                        & dtset%nsppol,dtset%typat,pawtab=pawtab,use_rhoijp=0)
   else
      pawrhoij_unsym => tdks%pawrhoij
   end if

   ! 3-Compute pawrhoij = \rho_{i,j} = \sum_{n,k}f_{n,k} \tilde{c}^{i,*}_{n,k} \tilde{c}^{j}_{n,k}
   call pawmkrhoij(tdks%atindx,tdks%atindx1,tdks%cprj,tdks%dimcprj,dtset%istwfk,    &
                 & dtset%kptopt,dtset%mband,tdks%mband_cprj,tdks%mcprj,dtset%mkmem, &
                 & mpi_enreg,dtset%natom,dtset%nband,dtset%nkpt,dtset%nspinor,      &
                 & dtset%nsppol,tdks%occ,dtset%paral_kgb,tdks%paw_dmft,             &
                 & pawrhoij_unsym,dtfil%unpaw,dtset%usewvl,dtset%wtk)

   ! 4-Symetrize rhoij, compute nhat and add it to rhor
   ! Note pawrhoij_unsym and pawrhoij are the same, which means that pawrhoij
   ! cannot be distributed over different atomic sites.
   ipert=0; idir=0; qpt(:)=zero; compch_fft=-1e-5_dp
   tdks%nhat = zero
   call pawmkrho(1,compch_fft,cplex,tdks%gprimd,idir,tdks%indsym,ipert,mpi_enreg, &
               & my_natom,dtset%natom,dtset%nspden,dtset%nsym,dtset%ntypat,       &
               & dtset%paral_kgb,pawang,tdks%pawfgr,tdks%pawfgrtab,               &
               & dtset%pawprtvol,tdks%pawrhoij,pawrhoij_unsym,pawtab,qpt,         &
               & rhowfg,rhowfr,tdks%rhor,tdks%rprimd,dtset%symafm,tdks%symrec,    &
               & dtset%typat,tdks%ucvol,dtset%usewvl,tdks%xred,pawnhat=tdks%nhat, &
               & pawnhatgr=tdks%nhatgr,rhog=tdks%rhog)

   ! 5-Take care of kinetic energy density
   if(dtset%usekden==1)then
     call mkrho(tdks%cg,dtset,tdks%gprimd,tdks%irrzon,tdks%kg,tdks%mcg,mpi_enreg, &
              & tdks%npwarr,tdks%occ,tdks%paw_dmft,tdks%phnons,rhowfg,rhowfr,     &
              & tdks%rprimd,tim_mkrho,tdks%ucvol,tdks%wvl%den,tdks%wvl%wfs,option=1)
     call transgrid(1,mpi_enreg,dtset%nspden,+1,1,1,dtset%paral_kgb,tdks%pawfgr, &
                  & rhowfg,tdks%taug,rhowfr,tdks%taur)
   end if

   ABI_FREE(rhowfg)
   ABI_FREE(rhowfr)

   if (my_natom/=dtset%natom) then
      call pawrhoij_free(pawrhoij_unsym)
      ABI_FREE(pawrhoij_unsym)
   else
      pawrhoij_unsym => NULL()
   end if

 else

   ABI_MALLOC(tdks%nhat,(0,0))
   ABI_MALLOC(tdks%nhatgr,(0,0,0))
   tdks%nhatgrdim=0

   ! 1-Compute density from WFs
   call mkrho(tdks%cg,dtset,tdks%gprimd,tdks%irrzon,tdks%kg,tdks%mcg,mpi_enreg,   &
            & tdks%npwarr,tdks%occ,tdks%paw_dmft,tdks%phnons,tdks%rhog,tdks%rhor, &
            & tdks%rprimd,tim_mkrho,tdks%ucvol,tdks%wvl%den,tdks%wvl%wfs)
   ! 2-Take care of kinetic energy density
   if(dtset%usekden==1)then
     call mkrho(tdks%cg,dtset,tdks%gprimd,tdks%irrzon,tdks%kg,tdks%mcg,mpi_enreg,   &
              & tdks%npwarr,tdks%occ,tdks%paw_dmft,tdks%phnons,tdks%taug,tdks%taur, &
              & tdks%rprimd,tim_mkrho,tdks%ucvol,tdks%wvl%den,tdks%wvl%wfs,option=1)
   end if

 endif

 end subroutine calc_density

end module m_rttddft_types
!!***
