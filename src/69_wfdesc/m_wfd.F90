!****m* ABINIT/m_wfd
!! NAME
!!  m_wfd
!!
!! FUNCTION
!!  This module contains the declaration of the wfd_t object.
!!  The wfd_t is a container of Bloch states (wave_t).
!!  It provides a high-level API to perform FFT transforms G --> R, compute PAW projections
!!
!! COPYRIGHT
!! Copyright (C) 2008-2020 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_wfd

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_copy
 use m_errors
 use m_crystal
 use m_wfk
 use m_hdr
 use m_distribfft
 use iso_c_binding

 use defs_datatypes,   only : pseudopotential_type, ebands_t
 use defs_abitypes,    only : mpi_type
 use m_gwdefs,         only : one_gw
 use m_time,           only : cwtime, cwtime_report
 use m_fstrings,       only : toupper, firstchar, int2char10, sjoin, itoa, strcat, itoa
 use m_io_tools,       only : get_unit, iomode_from_fname, iomode2str, open_file
 use m_numeric_tools,  only : imin_loc, list2blocks, bool2index
 use m_hide_blas,      only : xcopy, xdotc
 use m_pptools,        only : printxsf
 use m_cgtools,        only : cg_zdotc
 use m_cgtk,           only : cgtk_change_gsphere, cgtk_rotate
 use m_fftcore,        only : print_ngfft, kgindex, sphereboundary, ngfft_seq
 use m_fft_mesh,       only : rotate_fft_mesh, calc_ceikr, check_rot_fft
 use m_fft,            only : fft_ug !, fft_ug_dpc, fft_ur_dpc
 use m_kg,             only : getph, ph1d3d, mkkpg
 use m_gsphere,        only : kg_map, make_istwfk_table
 use m_fftcore,        only : kpgsph, get_kg
 use m_mpinfo,         only : nullify_mpi_enreg, destroy_mpi_enreg, copy_mpi_enreg, initmpi_seq
 use m_bz_mesh,        only : kmesh_t, get_bz_item
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type, pawtab_get_lsize
 use m_pawfgrtab,      only : pawfgrtab_type, pawfgrtab_init, pawfgrtab_free, pawfgrtab_print
 use m_pawcprj,        only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_copy, paw_overlap
 use m_paw_pwaves_lmn, only : paw_pwaves_lmn_t, paw_pwaves_lmn_init, paw_pwaves_lmn_free
 use m_pawrhoij,       only : pawrhoij_type, pawrhoij_mpisum_unpacked, pawrhoij_print_rhoij
 use m_paw_nhat,       only : nhatgrid
 use m_paw_occupancies,only : pawaccrhoij
 use m_iterators,      only : iter2_t, iter_yield, iter_len, iter_free, iter_push, iter_alloc
 use m_spacepar,       only : symrhg, irrzg
 use m_initylmg,       only : initylmg
 use m_mkffnl,         only : mkffnl
 use m_cgprj,          only : getcprj

 implicit none

 private
!!***

 ! Flags giving the status of the local %ug, %ur %cprj buffers.
 ! Use 1-byte int to save memory as much as possible.
 integer(c_int8_t),public,parameter :: WFD_NOWAVE    = 0
 integer(c_int8_t),public,parameter :: WFD_ALLOCATED = 1
 integer(c_int8_t),public,parameter :: WFD_STORED    = 2

 integer(c_int8_t),public,parameter :: CPR_RANDOM    = 1
 integer(c_int8_t),public,parameter :: CPR_SORTED    = 2

!----------------------------------------------------------------------

!!****t* m_wfd/kdata_t
!! NAME
!! kdata_t
!!
!! FUNCTION
!! Datatype storing k-dependent quantities and tables needed
!! for performing the zero-padded FFT of wavefunctions.
!!
!! SOURCE

 type,public :: kdata_t

   integer :: istwfk
   ! Storage mode for this k point.

   integer :: npw
   ! Number of plane-waves for this k-point.

   integer :: useylm
   ! 1 if nonlocal part is applied using real spherical Harmonics. 0 for Legendre polynomial.

   integer :: has_ylm
   ! 0 if ylm is not used.
   ! 1 if ylm is allocated.
   ! 2 if ylm is already computed.

   integer,allocatable :: kg_k(:,:)
   ! kg_k(3,npw)
   ! G vector coordinates in reduced cordinates.

   integer,allocatable :: gbound(:,:)
   ! gbound(2*mgfft+8,2))
   ! The boundary of the basis sphere of G vectors at a given k point.
   ! for use in improved zero padding of FFTs in 3 dimensions.

   !% real(dp) :: kpoint(3)

   real(dp),allocatable :: ph3d(:,:,:)
   ! ph3d(2, npw, natom)
   ! 3-dim structure factors, for each atom and each plane wave.
   ! Available only for PAW.

   real(dp),allocatable :: phkxred(:,:)
   ! phkxred(2,natom))
   ! e^{ik.Ra} for each atom. Packed according to the atom type (atindx).

   real(dp),allocatable :: fnl_dir0der0(:,:,:,:)
   ! fnl_dir0der0(npw,1,lmnmax,ntypat)
   ! nonlocal form factors. Computed only if usepaw == 1.
   ! fnl(k+G).ylm(k+G) if PAW
   ! f_ln(k+G)/|k+G|^l if NC

   real(dp),allocatable :: ylm(:,:)
   ! ylm(npw,mpsang**2*useylm)
   ! Real spherical harmonics for each k+G

 end type kdata_t

 public :: kdata_init
 public :: kdata_free
 public :: kdata_copy

 interface kdata_free
   module procedure kdata_free_0D
   module procedure kdata_free_1D
 end interface kdata_free

 interface kdata_copy
   module procedure copy_kdata_0D
   module procedure copy_kdata_1D
 end interface kdata_copy
!!***

!----------------------------------------------------------------------

!!****t* m_wfd/wave_t
!! NAME
!! wave_t
!!
!! FUNCTION
!!  Structure used to store a single wavefunction in G-space and, optionally, its r-space representation.
!!
!! SOURCE

 type, public :: wave_t

  !! integer :: cplex
  ! 1 for real wavefunctions u(r)
  ! 2 for complex wavefunctions u(r).
  ! At gamma we always have real u(r) provided that time-reversal can be used.
  ! In systems with both time-reversal and spatial inversion, wavefunctions can be chosen to be real.
  ! One might use this to reduce memory in wave_t.

  integer(c_int8_t) :: has_ug = WFD_NOWAVE
  ! Flag giving the status of ug.

  integer(c_int8_t) :: has_ur = WFD_NOWAVE
  ! Flag giving the status of ur.

  integer(c_int8_t) :: has_cprj = WFD_NOWAVE
  ! Flag giving the status of cprj.

  integer(c_int8_t) :: cprj_order = CPR_RANDOM
  ! Flag defining whether cprj are sorted by atom type or ordered according
  ! to the typat variable used in the input file.

  complex(gwpc),allocatable :: ug(:)
  ! ug(npw_k*nspinor)
  ! The periodic part of the Bloch wavefunction in G-space.

  complex(gwpc),allocatable :: ur(:)
  ! ur(nfft*nspinor)
  ! The periodic part of the Bloch wavefunction in real space.

  type(pawcprj_type),allocatable :: Cprj(:,:)
  ! Cprj(natom,nspinor)
  ! PAW projected wave function <Proj_i|Cnk> with all NL projectors.

  contains

  procedure :: free => wave_free
  ! Free memory

  procedure :: copy => wave_copy
  ! Copy object.

 end type wave_t

 public :: wave_init
!!***

!----------------------------------------------------------------------

!!****t* m_wfd/kpt_store_t
!! NAME
!!  kpt_store_t
!!
!! FUNCTION
!!  Used to build ragged arrays of wave_t in compact form.
!!
!! SOURCE

 type :: kpt_store_t
   type(wave_t),allocatable :: b(:)
 end type kpt_store_t
!!***

!----------------------------------------------------------------------

!!****t* m_wfd/spin_store_t
!! NAME
!!  spin_store_t
!!
!! FUNCTION
!!  Used to build ragged arrays of wave_t in compact form.
!!
!! SOURCE

 type :: spin_store_t
   type(kpt_store_t),allocatable :: k(:)
 end type spin_store_t
!!***

!----------------------------------------------------------------------

!!****t* m_wfd/wfd_t
!! NAME
!! wfd_t
!!
!! FUNCTION
!! Container gathering information on the set of wavefunctions treated by
!! this node as well as their distribution inside the MPI communicator.
!!
!! SOURCE

 type,public :: wfd_t

  integer :: debug_level = 0    ! Internal flag defining the debug level.
  integer :: lmnmax
  integer :: mband              ! MAX(nband)
  integer :: mgfft              ! Maximum size of 1D FFTs i.e. MAXVAL(ngfft(1:3)), used to dimension some arrays.
  integer :: natom
  integer :: nfft               ! Number of FFT points treated by this processor
  integer :: nfftot             ! Total number of points in the FFT grid
  integer :: nkibz              ! Number of irreducible k-points
  integer :: nspden             ! Number of independent spin-density components
  integer :: nspinor            ! Number of spinor components
  integer :: nsppol             ! Number of independent spin polarizations
  integer :: ntypat
  integer :: paral_kgb          ! Option for kgb parallelism
  integer :: usepaw             ! 1 if PAW is used, 0 otherwise.
  integer :: prtvol             ! Verbosity level.
  integer :: pawprtvol          ! Verbosity level for PAW.
  integer :: usewvl             ! 1 if BigDFT is used, 0 otherwise.
  integer :: comm               ! The MPI communicator for this pool of processors.
  integer :: master             ! The rank of master node in comm.
  integer :: my_rank            ! The rank of my processor inside the MPI communicator comm.
  integer :: nproc              ! The number of processors in MPI comm.

  integer :: my_nspins          ! Number of spins treated by this MPI proc

  integer,allocatable :: my_nkspin(:)
  ! (%nsppol))
  ! Number of k-points treated by this MPI proc.

  logical :: rfft_is_symok      ! .TRUE. if the real space FFT mesh is compatible with the rotational
                                ! part of the space group.

  real(dp) :: dilatmx

  real(dp) :: ecut
   ! Cutoff for plane wave basis set.

  real(dp) :: ecutsm
   ! ecutsm=smearing energy for plane wave kinetic energy (Ha)
   ! Cutoff for plane wave basis set.

!arrays
  integer :: ngfft(18)
   ! Information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft

  integer :: nloalg(3)
   ! Governs the choice of the algorithm for nonlocal operator. See doc.

  integer,allocatable :: irottb(:,:)
   ! irottb(nfftot,nsym)
   ! Index of $R^{-1}(r-\tau)$ in the FFT box.

  integer,allocatable :: istwfk(:)
   ! istwfk(nkibz)
   ! Storage mode for this k-point.

  integer,allocatable :: nband(:,:)
   ! nband(nkibz,nsppol)
   ! Number of bands at each k-point and spin.

  integer,allocatable :: indlmn(:,:,:)
   ! indlmn(6,lmnmax,ntypat)
   ! array giving l,m,n,lm,ln,spin for i=ln  (if useylm=0)
   !                                or i=lmn (if useylm=1)

  integer,allocatable :: nlmn_atm(:)
   ! nlmn_atm(natom)
   ! Number of (n,l,m) channels for each atom. Only for PAW

  integer,allocatable :: nlmn_sort(:)
   ! nlmn_sort(natom)
   ! Number of (n,l,m) channels for each atom (sorted by atom type). Only for PAW

  integer,allocatable :: nlmn_type(:)
   ! nlmn_type(ntypat)
   ! Number of (n,l,m) channels for each type of atom. Only for PAW.

  integer,allocatable :: npwarr(:)
   ! npwarr(nkibz)
   ! Number of plane waves for this k-point.

  integer, allocatable :: bks2wfd(:,:,:,:)
   ! (3, %mband, %nkibz, %nsppol)
   ! Maps global (band, ik_ibz, spin) to index in the wave store.
   ! Set to 0 if the (b, k, s) state is not in the store.

  integer(c_int8_t), private, allocatable :: bks_tab(:,:,:,:)
   ! bks_tab(mband,nkibz,nsppol,0:nproc-1)
   ! Global table used to keep trace of the distribution of the (b,k,s) states on each node inside Wfd%comm.
   ! 1 if the node has this state. 0 otherwise.
   ! A node owns a wavefunction if the corresponding ug is allocated AND computed.
   ! If a node owns ur but not ug, or ug is just allocated then its entry in the table is zero.

  integer, allocatable :: bks_comm(:,:,:)
   ! TODO: To be removed
   ! spin_comm(0:mband,0:nkibz,0:nsppol)
   ! MPI communicators.
   ! bks_comm(0,0,spin) MPI communicator for spin
   ! bks_comm(0,ik_ibz,spin)  MPI communicator for k-points.

  real(dp),allocatable :: kibz(:,:)
   ! kibz(3,nkibz)
   ! Reduced coordinates of the k-points in the IBZ.

  real(dp),allocatable :: ph1d(:,:)
   ! ph1d(2,3*(2*mgfft+1)*natom)
   ! 1-dim structure factor phase information.

  logical,private, allocatable :: keep_ur(:,:,:)
   ! TODO: To be removed
   ! keep(mband,nkibz,nsppol)
   ! Storage strategy: keep or not keep calculated u(r) in memory.

  type(kdata_t),allocatable :: Kdata(:)
   ! Kdata(nkibz)
   ! datatype storing k-dependent quantities.

  type(spin_store_t),allocatable :: s(:)
   ! (%my_nsppol)
   ! wfd%s(is)%k(ik)%b(ib)

  type(MPI_type) :: MPI_enreg
   ! The MPI_type structured datatype gather different information about the MPI parallelisation :
   ! number of processors, the index of my processor, the different groups of processors, etc ...

  !type(pseudopotential_type), pointer :: psps
  !type(pawtab_type), pointer :: pawtab(:)

 contains

   procedure :: free => wfd_free
   ! Destructor.

   !procedure :: copy => wfd_copy
   ! Copy routine

   procedure :: norm2 => wfd_norm2
   ! Compute <u(g)|u(g)> for the same k-point and spin.

   procedure :: xdotc => wfd_xdotc
   ! Compute <u_{b1ks}|u_{b2ks}> in G-space

   procedure :: reset_ur_cprj => wfd_reset_ur_cprj
   ! Reinitialize memory storage of u(r) and <p_i|psi>

   procedure :: get_many_ur => wfd_get_many_ur
   ! Get many wavefunctions in real space from its (bands(:),k,s) indices.

   procedure :: copy_cg => wfd_copy_cg
   ! Return a copy of u(g) in a real(2,npw_k)) array (Abinit convention)

   procedure :: get_ur => wfd_get_ur
   ! Get one wavefunction in real space from its (b,k,s) indices.

   procedure :: get_cprj => wfd_get_cprj
   ! Get one PAW projection <Proj_i|Cnk> with all NL projectors from its (b,k,s) indices.

   procedure :: change_ngfft => wfd_change_ngfft
   ! Reinitialize internal FFT tables.

   procedure :: print => wfd_print
   ! Printout of basic info.

   procedure :: ug2cprj => wfd_ug2cprj
   ! Get PAW cprj from its (b,k,s) indices.

   procedure :: wave_free => wfd_wave_free
   ! Free internal buffers used to store the wavefunctions.

   procedure :: get_wave_ptr => wfd_get_wave_ptr
   ! Return pointer to wave_t from its (b,k,s) indices

   procedure :: push_ug => wfd_push_ug
   ! Modify the value of u(g)_ks stored in the object.

   procedure :: extract_cgblock => wfd_extract_cgblock
   ! Extract a block of wavefunctions for a given spin and k-points (uses the cg storage mode)

   procedure :: ihave_ug => wfd_ihave_ug
   ! True if the node has this ug with the specified status.

   procedure :: mybands => wfd_mybands
   ! Returns the list of band indices of the u(g) owned by this node at given (k,s).

   procedure :: show_bkstab => wfd_show_bkstab
   ! Print a table showing the distribution of the wavefunctions.

   procedure :: distribute_bands => wfd_distribute_bands
   ! Distribute a set of bands taking into account the distribution of the ug.

   procedure :: iterator_bks => wfd_iterator_bks
   ! Iterator used to loop over bands, k-points and spin indices

   procedure :: bks_distrb => wfd_bks_distrb
   ! Distribute bands, k-points and spins

   procedure :: update_bkstab => wfd_update_bkstab
   ! Update the internal table with info on the distribution of the ugs.

   procedure :: set_mpicomm => wfd_set_mpicomm

   procedure :: rotate => wfd_rotate
   ! Linear transformation of the wavefunctions stored in Wfd

   procedure :: sanity_check => wfd_sanity_check
   ! Debugging tool

   procedure :: distribute_bbp => wfd_distribute_bbp
   ! Distribute a set of (b,b') indices

   procedure :: distribute_kb_kpbp => wfd_distribute_kb_kpbp

   procedure :: test_ortho => wfd_test_ortho
   ! Test the orthonormalization of the wavefunctions.

   procedure :: sym_ur => wfd_sym_ur
   ! Symmetrize a wave function in real space
   ! This routine is deprecated, see wfd_sym_ug_kg for algo in G-space.

   procedure :: sym_ug_kg => wfd_sym_ug_kg
   ! Symmetrize a wave function in G-space
   ! TODO: Can be removed ?
   ! Used in phgamma only, see sigmaph for a more efficient version.

   procedure :: paw_get_aeur => wfd_paw_get_aeur
   ! Compute the AE PAW wavefunction in real space.

   procedure :: plot_ur => wfd_plot_ur
   ! Write u(r) to an external file in XSF format.

   procedure :: write_wfk => wfd_write_wfk
   ! Write u(g) to a WFK file.

   procedure :: read_wfk => wfd_read_wfk
   ! Read u(g) from the WFK file completing the initialization of the object.

   procedure :: mkrho => wfd_mkrho
   ! Calculate the charge density on the fine FFT grid in real space.

   procedure :: pawrhoij => wfd_pawrhoij

   procedure :: dump_errinfo => wfd_dump_errinfo

 end type wfd_t

 public :: wfd_init                ! Main creation method.
 public :: wfd_copy
 !public :: wfd_get_socpert
 public :: test_charge
!!***

CONTAINS  !==============================================================================

!!****f* m_wfd/kdata_init
!! NAME
!!  kdata_init
!!
!! FUNCTION
!!  Main creation method for the kdata_t datatype.
!!
!! PARENTS
!!      debug_tools,m_wfd
!!
!! CHILDREN
!!
!! SOURCE

subroutine kdata_init(Kdata,Cryst,Psps,kpoint,istwfk,ngfft,MPI_enreg,ecut,kg_k)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwfk
 real(dp),optional,intent(in) :: ecut
 type(crystal_t),intent(in) :: Cryst
 type(pseudopotential_type),intent(in) :: Psps
 type(kdata_t),intent(inout) :: Kdata
 type(MPI_type),intent(in) :: MPI_enreg
!arrays
 integer,optional,target,intent(in) :: kg_k(:,:)
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: kpoint(3)

!Local variables ------------------------------
!scalars
 integer,parameter :: ider0=0,idir0=0
 integer :: mpw_,npw_k,dimffnl,useylmgr,nkpg,iatom, mkmem_,nkpt_,optder,mgfft
 integer :: iatm,matblk
 real(dp) :: arg
!arrays
 integer :: nband_(1),npwarr_(1)
 real(dp),allocatable :: ylmgr_k(:,:,:),kpg_k(:,:),ph1d(:,:)

!************************************************************************

 !@kdata_t
 Kdata%istwfk = istwfk
 Kdata%useylm = Psps%useylm

 if (present(ecut)) then
  ! Calculate G-sphere from input ecut.
  ABI_CHECK(.not.allocated(Kdata%kg_k), "Kdata%kg_k is allocated!")
  call get_kg(kpoint,istwfk,ecut,Cryst%gmet,npw_k,Kdata%kg_k)

 else if (present(kg_k)) then
   ! Use input g-vectors.
   npw_k = SIZE(kg_k,DIM=2)
   ABI_MALLOC(Kdata%kg_k,(3,npw_k))
   Kdata%kg_k = kg_k
 else
   MSG_ERROR("Either ecut or kg_k must be present")
 end if
 Kdata%npw = npw_k

 mgfft = MAXVAL(ngfft(1:3))

 ! Finds the boundary of the basis sphere of G vectors (for this k point)
 ! for use in improved zero padding of ffts in 3 dimensions.
 ABI_MALLOC(Kdata%gbound,(2*mgfft+8, 2))
 call sphereboundary(Kdata%gbound, istwfk, Kdata%kg_k, mgfft, npw_k)

 ! Compute e^{ik.Ra} for each atom. Packed according to the atom type (atindx).
 ABI_MALLOC(Kdata%phkxred,(2, Cryst%natom))
 do iatom=1,Cryst%natom
   iatm=Cryst%atindx(iatom)
   arg=two_pi*(DOT_PRODUCT(kpoint,Cryst%xred(:,iatom)))
   Kdata%phkxred(1,iatm)=DCOS(arg)
   Kdata%phkxred(2,iatm)=DSIN(arg)
 end do

 ! TODO: Should avoid storing all this stuff in memory (risky if lots of k-points)
 ! Write method to prepare kdata inside loop

 ! Calculate 1-dim structure factor phase information.
 mgfft = MAXVAL(ngfft(1:3))
 ABI_MALLOC(ph1d,(2, 3*(2*mgfft+1)*Cryst%natom))
 call getph(Cryst%atindx,Cryst%natom,ngfft(1),ngfft(2),ngfft(3),ph1d,Cryst%xred)

 matblk = 0; if (psps%usepaw == 1) matblk = Cryst%natom
 ABI_MALLOC(Kdata%ph3d,(2, npw_k, matblk))
 if (psps%usepaw == 1) then
   call ph1d3d(1,Cryst%natom,Kdata%kg_k,matblk,Cryst%natom,npw_k,ngfft(1),ngfft(2),ngfft(3),Kdata%phkxred,ph1d,Kdata%ph3d)
 end if
 ABI_FREE(ph1d)

 ! Compute spherical harmonics if required.
 Kdata%has_ylm = 0
 ABI_MALLOC(Kdata%ylm, (npw_k, Psps%mpsang**2*Psps%useylm))
 useylmgr=0
 ABI_MALLOC(ylmgr_k,(npw_k, 3, Psps%mpsang**2*useylmgr))

 if (Kdata%useylm == 1) then
   mkmem_=1; mpw_=npw_k; nband_=0; nkpt_=1; npwarr_(1)=npw_k
   optder=0 ! only Ylm(K) are computed.

   call initylmg(Cryst%gprimd,Kdata%kg_k,kpoint,mkmem_,MPI_enreg,Psps%mpsang,mpw_,nband_,nkpt_,&
    npwarr_,1,optder,Cryst%rprimd,Kdata%ylm,ylmgr_k)

   Kdata%has_ylm = 2
 end if

 ! Compute (k+G) vectors.
 nkpg = 0
 ABI_MALLOC(kpg_k,(npw_k,nkpg))
 if (nkpg>0) call mkkpg(Kdata%kg_k,kpg_k,kpoint,nkpg,npw_k)

 ! Compute nonlocal form factors fnl_dir0der0 for all (k+G).
 dimffnl = 0
 if (psps%usepaw == 1) dimffnl = 1+3*ider0
 ABI_MALLOC(Kdata%fnl_dir0der0,(npw_k,dimffnl,Psps%lmnmax,Cryst%ntypat))

 if (dimffnl /= 0) then
   call mkffnl(Psps%dimekb,dimffnl,Psps%ekb,Kdata%fnl_dir0der0,Psps%ffspl,&
     Cryst%gmet,Cryst%gprimd,ider0,idir0,Psps%indlmn,Kdata%kg_k,kpg_k,kpoint,Psps%lmnmax,&
     Psps%lnmax,Psps%mpsang,Psps%mqgrid_ff,nkpg,npw_k,Cryst%ntypat,&
     Psps%pspso,Psps%qgrid_ff,Cryst%rmet,Psps%usepaw,Psps%useylm,Kdata%ylm,ylmgr_k)
 end if

 ABI_FREE(kpg_k)
 ABI_FREE(ylmgr_k)

end subroutine kdata_init
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/kdata_free_0D
!! NAME
!!  kdata_free_0D
!!
!! FUNCTION
!!  Deallocate memory
!!
!! PARENTS
!!      m_wfd
!!
!! CHILDREN
!!
!! SOURCE

subroutine kdata_free_0D(Kdata)

!Arguments ------------------------------------
!scalars
 type(kdata_t),intent(inout) :: Kdata

!************************************************************************

 !@kdata_t
 ABI_SFREE(Kdata%kg_k)
 ABI_SFREE(Kdata%gbound)

 ABI_SFREE(Kdata%ph3d)
 ABI_SFREE(Kdata%phkxred)
 ABI_SFREE(Kdata%fnl_dir0der0)
 ABI_SFREE(Kdata%ylm)

end subroutine kdata_free_0D
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/kdata_free_1D
!! NAME
!!  kdata_free_1D
!!
!! FUNCTION
!!   Deallocate memory.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine kdata_free_1D(Kdata)

!Arguments ------------------------------------
!scalars
 type(kdata_t),intent(inout) :: Kdata(:)

!Local variables ------------------------------
!scalars
 integer :: ik

!************************************************************************

 do ik=LBOUND(Kdata,DIM=1),UBOUND(Kdata,DIM=1)
   call kdata_free_0D(Kdata(ik))
 end do

end subroutine kdata_free_1D
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/copy_kdata_0D
!! NAME
!!  copy_kdata_0D
!!
!! FUNCTION
!!  Deallocate memory
!!
!! PARENTS
!!      m_wfd
!!
!! CHILDREN
!!
!! SOURCE

subroutine copy_kdata_0D(Kdata_in,Kdata_out)

!Arguments ------------------------------------
!scalars
 type(kdata_t),intent(in) :: Kdata_in
 type(kdata_t),intent(inout) :: Kdata_out

!************************************************************************

 !@kdata_t
 Kdata_out%istwfk  = Kdata_in%istwfk
 Kdata_out%npw     = Kdata_in%npw
 Kdata_out%useylm  = Kdata_in%useylm
 Kdata_out%has_ylm = Kdata_in%has_ylm

 call alloc_copy(Kdata_in%kg_k, Kdata_out%kg_k)
 call alloc_copy(Kdata_in%gbound, Kdata_out%gbound)

 call alloc_copy(Kdata_in%ph3d,Kdata_out%ph3d)
 call alloc_copy(Kdata_in%phkxred,Kdata_out%phkxred)
 call alloc_copy(Kdata_in%fnl_dir0der0,Kdata_out%fnl_dir0der0)
 call alloc_copy(Kdata_in%ylm,Kdata_out%ylm)

end subroutine copy_kdata_0D
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/copy_kdata_1D
!! NAME
!!  copy_kdata_1D
!!
!! FUNCTION
!!   Deallocate memory.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine copy_kdata_1D(Kdata_in, Kdata_out)

!Arguments ------------------------------------
!scalars
 type(kdata_t),intent(in) :: Kdata_in(:)
 type(kdata_t),intent(inout) :: Kdata_out(:)

!Local variables ------------------------------
!scalars
 integer :: ik

!************************************************************************

 if (size(Kdata_in,DIM=1) /= size(Kdata_out,DIM=1)) then
   MSG_ERROR("copy_kdata_1D: wrong sizes !")
 end if

 do ik=LBOUND(Kdata_in,DIM=1),UBOUND(Kdata_in,DIM=1)
   call copy_kdata_0d(Kdata_in(ik),Kdata_out(ik))
 end do

end subroutine copy_kdata_1D
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_init
!! NAME
!! wfd_init
!!
!! FUNCTION
!!  Initialize the object.
!!
!! INPUTS
!!  Cryst<crystal_t>=Object defining the unit cell and its symmetries.
!!  Pawtab(ntypat*usepaw)<type(pawtab_type)>=PAW tabulated starting data.
!!  Psps<Pseudopotential_type>=datatype storing data on the pseudopotentials.
!!  ngfft(18)=All needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nkibz=Number of irreducible k-points.
!!  nsppol=Number of independent spin polarizations.
!!  nspden=Number of density components.
!!  nspinor=Number of spinorial components.
!!  ecut=Cutoff energy in Hartree
!!  ecutsm=Smearing for kinetic energy
!!  dilatmx
!!  mband
!!  nband(nkibz,nsppol)
!!  keep_ur(mband,nkibz,nsppol)=Option for memory storage of u(r).
!!  istwfk(nkibz)=Storage mode.
!!  kibz(3,nkibz)=Reduced coordinates of the k-points.
!!  nloalg(3)=Governs the choice of the algorithm for nonlocal operator. See doc.
!!  prtvol=Verbosity level.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  Initialize the object with basic dimensions, allocate also memory for u(g) and u(r) according to keep_ur
!!    %ug in G-space are always allocated.
!!    %ur in r-space only if keep_ur.
!!
!! PARENTS
!!      bethe_salpeter,m_gkk,m_phgamma,m_phpi,m_sigmaph,m_wfd
!!      screening,sigma,wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_init(Wfd,Cryst,Pawtab,Psps,keep_ur,mband,nband,nkibz,nsppol,bks_mask,&
&  nspden,nspinor,ecut,ecutsm,dilatmx,istwfk,kibz,ngfft,nloalg,prtvol,pawprtvol,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,comm,prtvol,pawprtvol
 integer,intent(in) :: nkibz,nsppol,nspden,nspinor
 real(dp),intent(in) :: ecut,ecutsm,dilatmx
 type(crystal_t),intent(in) :: Cryst
 type(pseudopotential_type),intent(in) :: Psps
 type(wfd_t),intent(inout) :: Wfd
!array
 integer,intent(in) :: ngfft(18),istwfk(nkibz),nband(nkibz,nsppol)
 integer,intent(in) :: nloalg(3)
 real(dp),intent(in) :: kibz(3,nkibz)
 logical,intent(in) :: bks_mask(mband,nkibz,nsppol)
 logical,intent(in) :: keep_ur(mband,nkibz,nsppol)
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: nfft0=0,mpw0=0,ikg0=0
 integer :: ik_ibz,spin,band,mpw,exchn2n3d,istwf_k,npw_k,iatom,itypat,iat
 integer :: cnt_b, cnt_k, cnt_s, ierr
 real(dp) :: ug_size,ur_size,cprj_size, bks_size
 real(dp) :: cpu, wall, gflops
 logical :: iscompatibleFFT
 character(len=500) :: msg
!arrays
 integer :: dum_kg(3,0)
 real(dp) :: kpoint(3)
 !integer :: my_band_list(Wfd%mband)

!************************************************************************

 DBG_ENTER("COLL")

 call cwtime(cpu, wall, gflops, "start")

 ! MPI info
 Wfd%comm    = comm
 Wfd%my_rank = xmpi_comm_rank(Wfd%comm)
 Wfd%nproc   = xmpi_comm_size(Wfd%comm)
 Wfd%master  = 0

 ABI_MALLOC(Wfd%bks_comm,(0:mband,0:nkibz,0:nsppol))
 Wfd%bks_comm = xmpi_comm_null

 ! Sequential MPI datatype to be passed to abinit routines.
 call initmpi_seq(Wfd%MPI_enreg)
 call init_distribfft(Wfd%MPI_enreg%distribfft,'c',Wfd%MPI_enreg%nproc_fft,ngfft(2),ngfft(3))

 ! TODO: To simply high-level API.
 !wfd%cryst => cryst
 !wfd%psps => psps
 !wfd%pawtab => pawtab

 ! Basic dimensions
 Wfd%nkibz     = nkibz
 Wfd%nsppol    = nsppol
 Wfd%nspden    = nspden
 Wfd%nspinor   = nspinor
 Wfd%paral_kgb = 0
 Wfd%nloalg    = nloalg

 Wfd%usepaw = Psps%usepaw
 Wfd%usewvl = 0 ! wavelets are not supported.
 Wfd%natom  = Cryst%natom
 Wfd%ntypat = Cryst%ntypat
 Wfd%lmnmax = Psps%lmnmax
 Wfd%prtvol = prtvol
 Wfd%pawprtvol = pawprtvol

 Wfd%ecutsm  = ecutsm
 Wfd%dilatmx = dilatmx

 ABI_MALLOC(Wfd%indlmn,(6,Wfd%lmnmax,Wfd%ntypat))
 Wfd%indlmn = Psps%indlmn

 if (Wfd%usepaw==1) then
   ABI_MALLOC(Wfd%nlmn_atm,(Cryst%natom))
   ABI_MALLOC(Wfd%nlmn_type,(Cryst%ntypat))
   do iatom=1,Cryst%natom
     Wfd%nlmn_atm(iatom)=Pawtab(Cryst%typat(iatom))%lmn_size
   end do

   do itypat=1,Cryst%ntypat
     Wfd%nlmn_type(itypat)=Pawtab(itypat)%lmn_size
   end do

   ABI_MALLOC(Wfd%nlmn_sort,(Cryst%natom))
   iat=0 ! nlmn dims sorted by atom type.
   do itypat=1,Cryst%ntypat
     Wfd%nlmn_sort(iat+1:iat+Cryst%nattyp(itypat))=Pawtab(itypat)%lmn_size
     iat=iat+Cryst%nattyp(itypat)
   end do
 end if

 ABI_MALLOC(Wfd%keep_ur,(mband,nkibz,nsppol))
 Wfd%keep_ur=keep_ur

 ! Setup of the FFT mesh
 Wfd%ngfft  = ngfft
 Wfd%mgfft  = MAXVAL (Wfd%ngfft(1:3))
 Wfd%nfftot = PRODUCT(Wfd%ngfft(1:3))
 Wfd%nfft   = Wfd%nfftot ! At present no FFT parallelism.

 Wfd%ecut = ecut

 ! Precalculate the FFT index of $ R^{-1} (r-\tau) $ used to symmetrize u_Rk.
 ABI_MALLOC(Wfd%irottb,(Wfd%nfftot,Cryst%nsym))
 call rotate_FFT_mesh(Cryst%nsym,Cryst%symrel,Cryst%tnons,Wfd%ngfft,Wfd%irottb,iscompatibleFFT)

 if (.not.iscompatibleFFT) then
   msg = "FFT mesh is not compatible with symmetries. Wavefunction symmetrization might be affected by large errors!"
   MSG_WARNING(msg)
 end if

 ! Is the real space mesh compatible with the rotational part?
 Wfd%rfft_is_symok = check_rot_fft(Cryst%nsym,Cryst%symrel,Wfd%ngfft(1),Wfd%ngfft(2),Wfd%ngfft(3))

 ABI_MALLOC(Wfd%kibz, (3, Wfd%nkibz))
 Wfd%kibz = kibz
 ABI_MALLOC(Wfd%istwfk, (Wfd%nkibz))
 Wfd%istwfk = istwfk

 ! Get the number of planewaves npw_k
 ! TODO Here we should use ecut_eff instead of ecut
 ABI_ICALLOC(Wfd%npwarr, (Wfd%nkibz))
 exchn2n3d = 0
 do ik_ibz=1,Wfd%nkibz
   if (mod(ik_ibz, wfd%nproc) /= wfd%my_rank) cycle ! MPI parallelism.
   istwf_k = Wfd%istwfk(ik_ibz)
   kpoint  = Wfd%kibz(:,ik_ibz)
   call kpgsph(Wfd%ecut,exchn2n3d,Cryst%gmet,ikg0,ik_ibz,istwf_k,dum_kg,kpoint,0,Wfd%MPI_enreg,mpw0,npw_k)
   Wfd%npwarr(ik_ibz)= npw_k
 end do
 call xmpi_sum(wfd%npwarr, wfd%comm, ierr)

 mpw = MAXVAL(Wfd%npwarr)

 ABI_MALLOC(Wfd%nband, (nkibz,nsppol))
 Wfd%nband=nband; Wfd%mband = mband
 ABI_CHECK(MAXVAL(Wfd%nband)==mband,"wrong mband")

 ! Allocate u(g) and, if required, also u(r)
 ug_size = one*nspinor*mpw*COUNT(bks_mask)
 write(msg,'(a,f8.1,a)')' Memory needed for Fourier components u(G) : ',two*gwpc*ug_size*b2Mb,' [Mb] <<< MEM'
 call wrtout(std_out, msg)
#ifdef HAVE_GW_DPC
 call wrtout(std_out, ' Storing wavefunctions in double precision array as `enable_gw_dpc="no"`')
 call wrtout(std_out, ' Recompile the code with `enable_gw_dpc="no"` to halve the memory requirements for the WFs')
#else
 call wrtout(std_out, ' Storing wavefunctions in single precision array as `enable_gw_dpc="no"`')
#endif

 if (Wfd%usepaw==1) then
   cprj_size = one * nspinor*SUM(Wfd%nlmn_atm)*COUNT(bks_mask)
   write(msg,'(a,f8.1,a)')' Memory needed for PAW projections Cprj: ',dp*cprj_size*b2Mb,' [Mb] <<< MEM'
   call wrtout(std_out, msg)
 end if

 ur_size = one*nspinor*Wfd%nfft*COUNT(Wfd%keep_ur)
 write(msg,'(a,f8.1,a)')' Memory needed for real-space u(r): ',two*gwpc*ur_size*b2Mb,' [Mb] <<< MEM'
 call wrtout(std_out, msg)

 ! Count the number of spins treated by this proc.
 wfd%my_nspins = 0
 do spin=1,wfd%nsppol
   if (any(bks_mask(:,:,spin))) wfd%my_nspins = wfd%my_nspins + 1
 end do
 !write(std_out, *)"my_nspins", wfd%my_nspins
 ABI_MALLOC(wfd%s, (wfd%my_nspins))

 ! Count the number of kpts in the IBZ treated by this proc. may be spin-dependent.
 ABI_ICALLOC(wfd%my_nkspin, (wfd%nsppol))
 cnt_s = 0
 do spin=1,wfd%nsppol
   do ik_ibz=1,wfd%nkibz
     if (any(bks_mask(:,ik_ibz,spin))) wfd%my_nkspin(spin) = wfd%my_nkspin(spin) + 1
   end do
   if (wfd%my_nkspin(spin) > 0) then
     cnt_s = cnt_s + 1
     ABI_MALLOC(wfd%s(cnt_s)%k, (wfd%my_nkspin(spin)))
   end if
 end do

 ! Allocate bands in packed format and use bks2wfd to go from global (b,k,s) index to local index.
 ABI_ICALLOC(wfd%bks2wfd, (3, wfd%mband, wfd%nkibz, wfd%nsppol))
 cnt_s = 0
 do spin=1,wfd%nsppol
   if (wfd%my_nkspin(spin) == 0) cycle
   cnt_s = cnt_s + 1
   cnt_k = 0
   do ik_ibz=1,wfd%nkibz
     cnt_b = count(bks_mask(:, ik_ibz, spin))
     if (cnt_b == 0) cycle
     cnt_k = cnt_k + 1
     ABI_MALLOC(wfd%s(cnt_s)%k(cnt_k)%b, (cnt_b))
     cnt_b = 0
     npw_k = Wfd%npwarr(ik_ibz)
     do band=1,Wfd%nband(ik_ibz, spin)
       if (bks_mask(band, ik_ibz, spin)) then
         cnt_b = cnt_b + 1
         call wave_init(wfd%s(cnt_s)%k(cnt_k)%b(cnt_b), &
                        Wfd%usepaw,npw_k,nfft0,Wfd%nspinor,Wfd%natom,Wfd%nlmn_atm,CPR_RANDOM)
         wfd%bks2wfd(:, band, ik_ibz, spin) = [cnt_b, cnt_k, cnt_s]
       end if
     end do
   end do
 end do

 bks_size = one * wfd%mband * wfd%nkibz * wfd%nsppol * wfd%nproc
 write(msg,'(a,f8.1,a)')' Memory needed for bks_tab: ',one * bks_size * b2Mb,' [Mb] <<< MEM'
 call wrtout(std_out, msg)

 !ABI_MALLOC(wfd%bks_ranks, (wfd%mband, nkibz, nsppol))

 ! Allocate the global table used to keep trace of the distribution, including a possible duplication.
 ABI_MALLOC(Wfd%bks_tab, (Wfd%mband, nkibz, nsppol, 0:Wfd%nproc-1))
 Wfd%bks_tab = WFD_NOWAVE

 ! Update the kbs table storing the distribution of the ug.
 call wfd%update_bkstab(show=-std_out)
 !
 ! Initialize the MPI communicators.
 ! init MPI communicators:cannot be done here since waves are not stored yet.
 !call wfd%set_mpicomm()
 !
 ! ===================================================
 ! ==== Precalculate nonlocal form factors for PAW ====
 ! ===================================================
 !
 ! Calculate 1-dim structure factor phase information.
 ABI_MALLOC(Wfd%ph1d, (2,3*(2*Wfd%mgfft+1)*Wfd%natom))
 call getph(Cryst%atindx,Wfd%natom,Wfd%ngfft(1),Wfd%ngfft(2),Wfd%ngfft(3),Wfd%ph1d,Cryst%xred)

 ! TODO: This one will require some memory if nkibz is large.
 ABI_MALLOC(Wfd%Kdata, (Wfd%nkibz))

 do ik_ibz=1,Wfd%nkibz
   kpoint  = Wfd%kibz(:,ik_ibz)
   istwf_k = Wfd%istwfk(ik_ibz)
   npw_k   = Wfd%npwarr(ik_ibz)
   if (any(wfd%bks2wfd(1, :, ik_ibz, :) /= 0)) then
     call kdata_init(Wfd%Kdata(ik_ibz), Cryst, Psps, kpoint, istwf_k, ngfft, Wfd%MPI_enreg, ecut=Wfd%ecut)
   end if
 end do

 call cwtime_report(" wfd_init", cpu, wall, gflops)

 DBG_EXIT("COLL")

end subroutine wfd_init
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_free
!! NAME
!!  wfd_free
!!
!! FUNCTION
!!  Free the memory allocated in the wfd_t data type.
!!
!! PARENTS
!!      bethe_salpeter,m_gkk,m_phgamma,m_phpi,m_sigmaph,screening
!!      sigma,wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_free(Wfd)

!Arguments ------------------------------------
!scalars
 class(wfd_t),intent(inout) :: Wfd

!Local variables ------------------------------
!scalars
 integer :: ib, ik, is

!************************************************************************

 ! integer.
 ABI_SFREE(Wfd%irottb)
 ABI_SFREE(Wfd%istwfk)
 ABI_SFREE(Wfd%nband)
 ABI_SFREE(Wfd%indlmn)
 ABI_SFREE(Wfd%nlmn_atm)
 ABI_SFREE(Wfd%nlmn_sort)
 ABI_SFREE(Wfd%nlmn_type)
 ABI_SFREE(Wfd%npwarr)
 ABI_SFREE(Wfd%bks_tab)

 if (allocated(wfd%s)) then
   do is=1,size(wfd%s)
     do ik=1,size(wfd%s(is)%k)
       do ib=1,size(wfd%s(is)%k(ik)%b)
         call wfd%s(is)%k(ik)%b(ib)%free()
       end do
       ABI_FREE(wfd%s(is)%k(ik)%b)
     end do
     ABI_FREE(wfd%s(is)%k)
   end do
   ABI_FREE(wfd%s)
 end if

 ABI_SFREE(wfd%my_nkspin)
 ABI_SFREE(wfd%bks2wfd)

 ! Free the MPI communicators.
 if (allocated(Wfd%bks_comm)) then
   call xmpi_comm_free(Wfd%bks_comm)
   ABI_FREE(Wfd%bks_comm)
 end if

 ! real arrays.
 ABI_SFREE(Wfd%kibz)
 ABI_SFREE(Wfd%ph1d)

 ! logical arrays.
 ABI_SFREE(Wfd%keep_ur)

 ! datatypes.
 if (allocated(Wfd%Kdata)) then
   call kdata_free(Wfd%Kdata)
   ABI_FREE(Wfd%Kdata)
 end if

 call destroy_mpi_enreg(Wfd%MPI_enreg)

end subroutine wfd_free
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_copy
!! NAME
!!  wfd_copy
!!
!! FUNCTION
!!  Copy a wfd_t data type.
!!
!! PARENTS
!!      screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_copy(Wfd_in, Wfd_out)

!Arguments ------------------------------------
 type(wfd_t),intent(inout) :: Wfd_in,Wfd_out

!Local variables ------------------------------
!scalars
 integer :: band, ik_ibz, spin, cnt_s, cnt_b, ib, ik, is

!************************************************************************

 !@wfd_t
 call deep_copy(Wfd_in%debug_level    ,Wfd_out%debug_level)
 call deep_copy(Wfd_in%lmnmax         ,Wfd_out%lmnmax)
 call deep_copy(Wfd_in%mband          ,Wfd_out%mband)
 call deep_copy(Wfd_in%mgfft          ,Wfd_out%mgfft)
 call deep_copy(Wfd_in%natom          ,Wfd_out%natom)
 call deep_copy(Wfd_in%nfft           ,Wfd_out%nfft)
 call deep_copy(Wfd_in%nfftot         ,Wfd_out%nfftot)
 call deep_copy(Wfd_in%nkibz          ,Wfd_out%nkibz)
 call deep_copy(Wfd_in%nspden         ,Wfd_out%nspden)
 call deep_copy(Wfd_in%nspinor        ,Wfd_out%nspinor)
 call deep_copy(Wfd_in%nsppol         ,Wfd_out%nsppol)
 call deep_copy(Wfd_in%ntypat         ,Wfd_out%ntypat)
 call deep_copy(Wfd_in%paral_kgb      ,Wfd_out%paral_kgb)
 call deep_copy(Wfd_in%usepaw         ,Wfd_out%usepaw)
 call deep_copy(Wfd_in%prtvol         ,Wfd_out%prtvol)
 call deep_copy(Wfd_in%pawprtvol      ,Wfd_out%pawprtvol)
 call deep_copy(Wfd_in%usewvl         ,Wfd_out%usewvl)
 call deep_copy(Wfd_in%comm           ,Wfd_out%comm)
 call deep_copy(Wfd_in%master         ,Wfd_out%master)
 call deep_copy(Wfd_in%my_rank        ,Wfd_out%my_rank)
 call deep_copy(Wfd_in%nproc          ,Wfd_out%nproc)
 call deep_copy(Wfd_in%my_nspins      ,Wfd_out%my_nspins)
 call deep_copy(Wfd_in%rfft_is_symok  ,Wfd_out%rfft_is_symok)
 call deep_copy(Wfd_in%dilatmx        ,Wfd_out%dilatmx)
 call deep_copy(Wfd_in%ecut           ,Wfd_out%ecut)
 call deep_copy(Wfd_in%ecutsm         ,Wfd_out%ecutsm)

!arrays
 Wfd_out%ngfft =Wfd_in%ngfft
 Wfd_out%nloalg=Wfd_in%nloalg

 call alloc_copy(Wfd_in%my_nkspin     ,Wfd_out%my_nkspin)
 call alloc_copy(Wfd_in%irottb        ,Wfd_out%irottb)
 call alloc_copy(Wfd_in%istwfk        ,Wfd_out%istwfk)
 call alloc_copy(Wfd_in%nband         ,Wfd_out%nband)
 call alloc_copy(Wfd_in%indlmn        ,Wfd_out%indlmn)
 call alloc_copy(Wfd_in%nlmn_atm      ,Wfd_out%nlmn_atm)
 call alloc_copy(Wfd_in%nlmn_sort     ,Wfd_out%nlmn_sort)
 call alloc_copy(Wfd_in%nlmn_type     ,Wfd_out%nlmn_type)
 call alloc_copy(Wfd_in%npwarr        ,Wfd_out%npwarr)
 call alloc_copy(Wfd_in%kibz          ,Wfd_out%kibz)
 call alloc_copy(Wfd_in%bks2wfd       ,Wfd_out%bks2wfd)
 call alloc_copy(Wfd_in%bks_tab       ,Wfd_out%bks_tab)
 call alloc_copy(Wfd_in%bks_comm      ,Wfd_out%bks_comm)
 call alloc_copy(Wfd_in%ph1d          ,Wfd_out%ph1d)
 call alloc_copy(Wfd_in%keep_ur       ,Wfd_out%keep_ur)

 ! types
 if (size(Wfd_in%Kdata,DIM=1) .ne. size(Wfd_out%Kdata,DIM=1)) then
   ABI_REMALLOC(Wfd_out%Kdata, (Wfd_out%nkibz))
 end if

 call kdata_copy(Wfd_in%Kdata, Wfd_out%Kdata)

 ! Allocate ragged array.
 ABI_MALLOC(wfd_out%s, (wfd_out%my_nspins))
 cnt_s = 0
 do spin=1,wfd_out%nsppol
   if (wfd_out%my_nkspin(spin) > 0) then
     cnt_s = cnt_s + 1
     ABI_MALLOC(wfd_out%s(cnt_s)%k, (wfd_out%my_nkspin(spin)))
     do ik=1,wfd_out%my_nkspin(spin)
        cnt_b = size(wfd_out%s(cnt_s)%k(ik)%b)
        ABI_MALLOC(wfd_out%s(cnt_s)%k(ik)%b, (cnt_b))
     end do
   end if
 end do

 ! Copy waves
 do spin=1,wfd_in%nsppol
   do ik_ibz=1,wfd_in%nkibz
     do band=1,wfd_in%nband(ik_ibz, spin)
       ib = wfd_in%bks2wfd(1, band, ik_ibz, spin)
       ik = wfd_in%bks2wfd(2, band, ik_ibz, spin)
       is = wfd_in%bks2wfd(3, band, ik_ibz, spin)
       if (ib /= 0) wfd_out%s(is)%k(ik)%b(ib) = wfd_in%s(is)%k(ik)%b(ib)%copy()
     end do
   end do
 end do

 call copy_mpi_enreg(Wfd_in%MPI_enreg, Wfd_out%MPI_enreg)

end subroutine wfd_copy
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_norm2
!! NAME
!!  wfd_norm2
!!
!! FUNCTION
!!   Compute <u_{bks}|u_{bks}> in G-space
!!
!! INPUTS
!!  Wfd<wfd_t>=the wavefunction descriptor.
!!  Cryst<crystal_t>=Structure describing the crystal structure and its symmetries.
!!  Pawtab(ntypat*usepaw)<type(pawtab_type)>=PAW tabulated starting data.
!!  band=Band index.
!!  ik_bz=Index of the k-point in the BZ.
!!  spin=Spin index
!!
!! PARENTS
!!
!! SOURCE

function wfd_norm2(Wfd,Cryst,Pawtab,band,ik_ibz,spin) result(norm2)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin
 real(dp) :: norm2
 type(crystal_t),intent(in) :: Cryst
 class(wfd_t),target,intent(inout) :: Wfd
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)

!Local variables ------------------------------
!scalars
 integer :: npw_k,istwf_k
 complex(dpc) :: cdum
 type(wave_t),pointer :: wave
 character(len=500) :: msg
!arrays
 real(dp) :: pawovlp(2)
 complex(gwpc),ABI_CONTIGUOUS pointer :: ug1(:)
 type(pawcprj_type),allocatable :: Cp1(:,:)

!************************************************************************

 ! Planewave part.
 npw_k   = Wfd%npwarr(ik_ibz)
 istwf_k = Wfd%istwfk(ik_ibz)

 ABI_CHECK(wfd%get_wave_ptr(band, ik_ibz, spin, wave, msg) == 0, msg)

 ug1 => wave%ug
 cdum = xdotc(Wfd%nspinor*npw_k,ug1,1,ug1,1)

 if (istwf_k>1) then
   cdum=two*DBLE(cdum)
   if (istwf_k==2) cdum=cdum-CONJG(ug1(1))*ug1(1)
 end if

 ! Paw on-site term.
 if (Wfd%usepaw==1) then

   ! Avoid the computation if Cprj are already in memory with the correct order.
   if (wave%has_cprj == WFD_STORED .and. wave%cprj_order == CPR_RANDOM) then

       pawovlp = paw_overlap(wave%Cprj, wave%Cprj, Cryst%typat, Pawtab)
       cdum = cdum + CMPLX(pawovlp(1),pawovlp(2), kind=dpc)

   else
     ! Compute Cproj
     ABI_MALLOC(Cp1,(Wfd%natom,Wfd%nspinor))
     call pawcprj_alloc(Cp1,0,Wfd%nlmn_atm)

     call wfd%get_cprj(band,ik_ibz,spin,Cryst,Cp1,sorted=.FALSE.)
     pawovlp = paw_overlap(Cp1,Cp1,Cryst%typat,Pawtab)
     cdum = cdum + CMPLX(pawovlp(1),pawovlp(2), kind=dpc)

     call pawcprj_free(Cp1)
     ABI_FREE(Cp1)
   end if
 end if

 norm2 = DBLE(cdum)

end function wfd_norm2
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_xdotc
!! NAME
!!  wfd_xdotc
!!
!! FUNCTION
!!   Compute <u_{b1ks}|u_{b2ks}> in G-space
!!
!! INPUTS
!!  Wfd<wfd_t>=the wavefunction descriptor.
!!  Cryst<crystal_t>=Structure describing the crystal structure and its symmetries.
!!  Pawtab(ntypat*usepaw)<type(pawtab_type)>=PAW tabulated starting data.
!!  band1, band2=Band indices.
!!  ik_bz=Index of the k-point in the BZ.
!!  spin=Spin index
!!
!! PARENTS
!!
!! SOURCE

function wfd_xdotc(Wfd,Cryst,Pawtab,band1,band2,ik_ibz,spin)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band1,band2,ik_ibz,spin
 complex(gwpc) :: wfd_xdotc
 class(wfd_t),target,intent(inout) :: Wfd
 type(crystal_t),intent(in) :: Cryst
!arrays
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)

!Local variables ------------------------------
!scalars
 integer :: npw_k,istwf_k
 type(wave_t),pointer :: wave1, wave2
 character(len=500) :: msg
!arrays
 real(dp) :: pawovlp(2)
 complex(gwpc),ABI_CONTIGUOUS pointer :: ug1(:),ug2(:)
 type(pawcprj_type),allocatable :: Cp1(:,:),Cp2(:,:)

!************************************************************************

 ! Planewave part.
 npw_k   = Wfd%npwarr(ik_ibz)
 istwf_k = Wfd%istwfk(ik_ibz)

 ABI_CHECK(wfd%get_wave_ptr(band1, ik_ibz, spin, wave1, msg) == 0, msg)
 ABI_CHECK(wfd%get_wave_ptr(band2, ik_ibz, spin, wave2, msg) == 0, msg)
 ug1 => wave1%ug
 ug2 => wave2%ug

 wfd_xdotc = xdotc(npw_k*Wfd%nspinor,ug1,1,ug2,1)
 if (istwf_k>1) then
   wfd_xdotc=two*DBLE(wfd_xdotc)
   if (istwf_k==2) wfd_xdotc = wfd_xdotc-CONJG(ug1(1))*ug2(1)
 end if

 ! Paw on-site term.
 if (Wfd%usepaw==1) then
   ! Avoid the computation if Cprj are already in memory with the correct order.
   if (wave1%has_cprj == WFD_STORED .and. wave1%cprj_order == CPR_RANDOM .and. &
       wave2%has_cprj == WFD_STORED .and. wave2%cprj_order == CPR_RANDOM) then

       pawovlp = paw_overlap(wave1%Cprj, wave2%Cprj,&
                             Cryst%typat,Pawtab,spinor_comm=Wfd%MPI_enreg%comm_spinor)
       wfd_xdotc = wfd_xdotc + CMPLX(pawovlp(1),pawovlp(2), kind=gwpc)

   else
     ! Compute Cprj
     ABI_MALLOC(Cp1,(Wfd%natom,Wfd%nspinor))
     call pawcprj_alloc(Cp1,0,Wfd%nlmn_atm)
     ABI_MALLOC(Cp2,(Wfd%natom,Wfd%nspinor))
     call pawcprj_alloc(Cp2,0,Wfd%nlmn_atm)

     call wfd%get_cprj(band1,ik_ibz,spin,Cryst,Cp1,sorted=.FALSE.)
     call wfd%get_cprj(band2,ik_ibz,spin,Cryst,Cp2,sorted=.FALSE.)

     pawovlp = paw_overlap(Cp1,Cp2,Cryst%typat,Pawtab,spinor_comm=Wfd%MPI_enreg%comm_spinor)
     wfd_xdotc = wfd_xdotc + CMPLX(pawovlp(1),pawovlp(2), kind=gwpc)

     call pawcprj_free(Cp1)
     ABI_FREE(Cp1)
     call pawcprj_free(Cp2)
     ABI_FREE(Cp2)
   end if
 end if

end function wfd_xdotc
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_reset_ur_cprj
!! NAME
!!  wfd_reset_ur_cprj
!!
!! FUNCTION
!!  Reinitialize the storage mode of the ur treated by this node.
!!
!! PARENTS
!!      bethe_salpeter,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_reset_ur_cprj(Wfd)

!Arguments ------------------------------------
 class(wfd_t),intent(inout) :: Wfd

!Local variables ------------------------------
!scalars
 integer :: ib, ik, is
!************************************************************************

 do is=1,size(wfd%s)
   do ik=1,size(wfd%s(is)%k)
     do ib=1,size(wfd%s(is)%k(ik)%b)
       if (wfd%s(is)%k(ik)%b(ib)%has_ur == WFD_STORED) wfd%s(is)%k(ik)%b(ib)%has_ur = WFD_ALLOCATED
       if (wfd%usepaw == 1) then
         if (wfd%s(is)%k(ik)%b(ib)%has_cprj == WFD_STORED) wfd%s(is)%k(ik)%b(ib)%has_cprj = WFD_ALLOCATED
       end if
     end do
   end do
 end do

end subroutine wfd_reset_ur_cprj
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_get_many_ur
!! NAME
!!  wfd_get_many_ur
!!
!! FUNCTION
!!  Get many wave functions in real space, either by doing a G-->R FFT
!!  or by just retrieving the data already stored in Wfd.
!!
!! INPUTS
!!  Wfd<wfd_t>=the wavefunction descriptor.
!!  ndat=Number of wavefunctions required
!!  bands(:)=Band indices.
!!  ik_ibz=Index of the k-point in the IBZ.
!!  spin=Spin index
!!
!! OUTPUT
!!  ur(Wfd%nfft*Wfd%nspinor*SIZE(bands))=The wavefunction in real space.
!!
!! PARENTS
!!      calc_sigc_me,calc_sigx_me,cohsex_me
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_get_many_ur(Wfd, bands, ik_ibz, spin, ur)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,spin
 class(wfd_t),intent(inout) :: Wfd
!arrays
 integer,intent(in) :: bands(:)
 complex(gwpc),intent(out) :: ur(Wfd%nfft*Wfd%nspinor*SIZE(bands))

!Local variables ------------------------------
!scalars
 integer :: dat,ptr,band
!************************************************************************

 do dat=1,SIZE(bands)
   band = bands(dat)
   ptr = 1 + (dat-1)*Wfd%nfft*Wfd%nspinor
   call wfd%get_ur(band,ik_ibz,spin,ur(ptr))
 end do

end subroutine wfd_get_many_ur
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_copy_cg
!! NAME
!!  wfd_copy_cg
!!
!! FUNCTION
!!  Return a copy u(g) in a real(dp) array. Useful if we have to interface
!!  the wavefunction descriptor with Abinit code expecting cg(2,npw_k*nspinor) arrays
!!  The routine takes also into account the fact that the ug in wfs could be stored in single-precision.
!!
!! INPUTS
!!  wfd<wfd_t>=the wavefunction descriptor.
!!  band=Band index.
!!  ik_ibz=Index of the k-point in the IBZ.
!!  spin=Spin index
!!
!! OUTPUT
!!  cg(npw_k*nspinor)=The wavefunction in real space in the Abinit cg convention.
!!
!! PARENTS
!!      m_gkk,m_phgamma,m_phpi,m_sigmaph
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_copy_cg(wfd, band, ik_ibz, spin, cg)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin
 class(wfd_t),intent(in) :: wfd
!arrays
 real(dp),intent(out) :: cg(2,*) ! npw_k*wfd%nspinor)

!Local variables ------------------------------
!scalars
 integer :: siz
 type(wave_t),pointer :: wave
 character(len=500) :: msg
!************************************************************************

 ABI_CHECK(wfd%get_wave_ptr(band, ik_ibz, spin, wave, msg) == 0, msg)

 if (.not. wave%has_ug == WFD_STORED) then
   write(msg,'(a,3(i0,1x),a)')" ug for (band, ik_ibz, spin): ",band,ik_ibz,spin," is not stored in memory!"
   MSG_ERROR(msg)
 end if

 siz = wfd%npwarr(ik_ibz) * wfd%nspinor
#ifdef HAVE_GW_DPC
 call zcopy(siz, wave%ug, 1, cg, 1)
#else
 cg(1,1:siz) = dble(wave%ug)
 cg(2,1:siz) = aimag(wave%ug)
#endif

end subroutine wfd_copy_cg
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_get_ur
!! NAME
!!  wfd_get_ur
!!
!! FUNCTION
!!  Get a wave function in real space, either by doing a G-->R FFT
!!  or by just retrieving the data already stored in Wfd.
!!
!! INPUTS
!!  Wfd<wfd_t>=the wavefunction descriptor.
!!  band=Band index.
!!  ik_ibz=Index of the k-point in the IBZ.
!!  spin=Spin index
!!
!! OUTPUT
!!  ur(Wfd%nfft*Wfd%nspinor)=The wavefunction in real space.
!!
!! PARENTS
!!      calc_sigc_me,calc_sigx_me,calc_vhxc_me,cchi0,cchi0q0,cchi0q0_intraband
!!      classify_bands,cohsex_me,exc_build_block,exc_build_ham,exc_den
!!      m_wfd,prep_calc_ucrpa,wfd_mkrho
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_get_ur(Wfd, band, ik_ibz, spin, ur)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin
 class(wfd_t),target,intent(inout) :: Wfd
!arrays
 complex(gwpc),intent(out) :: ur(Wfd%nfft*Wfd%nspinor)

!Local variables ------------------------------
!scalars
 integer,parameter :: npw0=0,ndat1=1
 integer :: npw_k, nfft, nspinor
 character(len=500) :: msg
 type(wave_t),pointer :: wave
!arrays
 integer,ABI_CONTIGUOUS pointer :: kg_k(:,:),gbound(:,:)
 complex(gwpc),ABI_CONTIGUOUS pointer :: ug(:)
!************************************************************************

 npw_k  = Wfd%npwarr(ik_ibz)
 nfft   = Wfd%nfft
 nspinor= Wfd%nspinor

 ABI_CHECK(wfd%get_wave_ptr(band, ik_ibz, spin, wave, msg) == 0, msg)

 select case (wave%has_ur)

 case (WFD_NOWAVE, WFD_ALLOCATED)
   ! FFT is required.
   if (.not. wave%has_ug == WFD_STORED) then
     write(msg,'(a,3(i0,1x),a)')" ug for (band, ik_ibz, spin): ",band,ik_ibz,spin," is not stored in memory!"
     MSG_ERROR(msg)
   end if

   ug => wave%ug
   kg_k    => Wfd%Kdata(ik_ibz)%kg_k
   gbound  => Wfd%Kdata(ik_ibz)%gbound(:,:)

   call fft_ug(npw_k,nfft,nspinor,ndat1,Wfd%mgfft,Wfd%ngfft,Wfd%istwfk(ik_ibz),kg_k,gbound,ug,ur)

   if (Wfd%keep_ur(band,ik_ibz,spin)) then
     ! Store results
     if (wave%has_ur == WFD_NOWAVE) then
       ! Alloc buffer for ur.
       call wave_init(wave, Wfd%usepaw,npw0,nfft,nspinor,Wfd%natom,Wfd%nlmn_atm,CPR_RANDOM)
     end if
     call xcopy(nfft*nspinor, ur, 1, wave%ur, 1)
     wave%has_ur = WFD_STORED
   end if

 case (WFD_STORED)
   ! copy it back.
   call xcopy(nfft*nspinor, wave%ur, 1, ur, 1)

 case default
   MSG_BUG(sjoin("Wrong has_ur:", itoa(wave%has_ur)))
 end select

end subroutine wfd_get_ur
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_print
!! NAME
!! wfd_print
!!
!! FUNCTION
!!  Print the content of a wfd_t datatype
!!
!! INPUTS
!!  Wfd<wfd_t>=The datatype.
!!  [header]=String to be printed as header for additional info.
!!  [unit]=Unit number for output
!!  [prtvol]=Verbosity level
!!  [mode_paral]=Either "COLL" or "PERS". Defaults to "COLL".
!!
!! OUTPUT
!!  Only printing
!!
!! PARENTS
!!      bethe_salpeter,m_gkk,m_phgamma,m_phpi,m_sigmaph,screening
!!      sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_print(Wfd, header, unit, prtvol, mode_paral)

!Arguments ------------------------------------
 integer,optional,intent(in) :: unit,prtvol
 character(len=4),optional,intent(in) :: mode_paral
 character(len=*),optional,intent(in) :: header
 class(wfd_t),intent(in) :: Wfd

!Local variables-------------------------------
!scalars
 integer :: my_prtvol, my_unt, mpw, ib, ik, is, ug_cnt, ur_cnt, cprj_cnt, spin, ik_ibz, band
 real(dp) :: ug_size, ur_size, cprj_size !,kdata_bsize
 character(len=4) :: my_mode
 character(len=500) :: msg
! *************************************************************************

 my_unt   =std_out; if (present(unit      )) my_unt   =unit
 my_prtvol=0      ; if (present(prtvol    )) my_prtvol=prtvol
 my_mode  ='COLL' ; if (present(mode_paral)) my_mode  =mode_paral

 msg=' ==== Info on the Wfd% object ==== '
 if (present(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 write(msg,'(3(a,i0,a),a,i0,2a,f5.1)')&
   '  Number of irreducible k-points ........ ',Wfd%nkibz,ch10,&
   '  Number of spinorial components ........ ',Wfd%nspinor,ch10,&
   '  Number of spin-density components ..... ',Wfd%nspden,ch10,&
   '  Number of spin polarizations .......... ',Wfd%nsppol,ch10,&
   '  Plane wave cutoff energy .............. ',Wfd%ecut
 call wrtout(my_unt, msg, my_mode)

 mpw = maxval(Wfd%npwarr)
 write(msg,'(3(a,i0,a))')&
   '  Max number of G-vectors ............... ',mpw,ch10,&
   '  Total number of FFT points ............ ',Wfd%nfftot,ch10,&
   '  Number of FFT points treated by me .... ',Wfd%nfft,ch10
 call wrtout(my_unt, msg, my_mode)

 call print_ngfft(Wfd%ngfft, 'FFT mesh for wavefunctions', my_unt, my_mode, my_prtvol)

 ug_cnt = 0; ur_cnt = 0; cprj_cnt = 0
 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz
     do band=1,Wfd%nband(ik_ibz,spin)
       ib = wfd%bks2wfd(1, band, ik_ibz, spin)
       if (ib  == 0) cycle
       ik = wfd%bks2wfd(2, band, ik_ibz, spin)
       is = wfd%bks2wfd(3, band, ik_ibz, spin)
       if (wfd%s(is)%k(ik)%b(ib)%has_ug >= WFD_ALLOCATED) ug_cnt = ug_cnt + 1
       if (wfd%s(is)%k(ik)%b(ib)%has_ur >= WFD_ALLOCATED) ur_cnt = ur_cnt + 1
       if (wfd%s(is)%k(ik)%b(ib)%has_cprj >= WFD_ALLOCATED) cprj_cnt = cprj_cnt + 1
     end do
   end do
 end do

 ! Info on memory needed for u(g), u(r) and PAW cprj
 write(msg, '(a,i0)')' Number of Bloch states treated by this rank: ', ug_cnt
 call wrtout(std_out, msg, pre_newlines=1)

 ug_size = one * Wfd%nspinor * mpw * ug_cnt
 write(msg,'(a,f8.1,a)')' Memory allocated for Fourier components u(G) = ',two*gwpc*ug_size*b2Mb,' [Mb] <<< MEM'
 call wrtout(std_out, msg)

 ur_size = one * Wfd%nspinor * Wfd%nfft * ur_cnt
 write(msg,'(a,f8.1,a)')' Memory allocated for real-space u(r) = ',two*gwpc*ur_size*b2Mb,' [Mb] <<< MEM'
 call wrtout(std_out, msg)

 if (wfd%usepaw==1) then
   cprj_size = one * Wfd%nspinor * sum(Wfd%nlmn_atm) * cprj_cnt
   write(msg,'(a,f8.1,a)')' Memory allocated for PAW projections Cprj = ',dp*cprj_size*b2Mb,' [Mb] <<< MEM'
   call wrtout(std_out, msg)
 end if

 !TODO
 ! Add additionanl info
 !kdata_bsize = nkibz * (four * (3 * mpw) + dp * two * mpw * natom)
 !write(msg,'(a,f8.1,a)')' Memory allocated for Kdata = ',kdata_bsize * b2Mb,' [Mb] <<< MEM'

 write(msg,'(a,f8.1,a)')' Memory needed for wfd%s datastructure: ',ABI_MEM_MB(wfd%s),' [Mb] <<< MEM'
 call wrtout(std_out, msg)
 write(msg,'(a,f8.1,a)')' Memory needed for wfd%s(0)%k datastructure: ',ABI_MEM_MB(wfd%s(1)%k),' [Mb] <<< MEM'
 call wrtout(std_out, msg)
 write(msg,'(a,f8.1,a)')' Memory allocated for Kdata array: ',ABI_MEM_MB(wfd%kdata),' [Mb] <<< MEM'
 call wrtout(std_out, msg, newlines=1)

end subroutine wfd_print
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_ug2cprj
!! NAME
!! wfd_ug2cprj
!!
!! FUNCTION
!!  Calculate the projected wave function <Proj_i|Cnk> with all NL projectors for a single
!!  k-point, band and spin.
!!
!! INPUTS
!!  Wfd<wfd_t>=Structure containing the wave functions for the GW.
!!  ik_ibz=Index of the required k-point
!!  spin=Required spin index.
!!  choice=chooses possible output:
!!    In addition to projected wave function:
!!    choice=1 => nothing else
!!          =2 => 1st gradients with respect to atomic position(s)
!!          =3 => 1st gradients with respect to strain(s)
!!          =23=> 1st gradients with respect to atm. pos. and strain(s)
!!          =4 => 2nd derivatives with respect to atomic pos.
!!          =24=> 1st and 2nd derivatives with respect to atomic pos.
!!          =5 => 1st gradients with respect to k wavevector
!!          =6 => 2nd derivatives with respect to strain and atm. pos.
!!  idir=direction of the derivative, i.e. dir. of - atom to be moved  in the case choice=2
!!                                                 - strain component  in the case choice=3
!!                                                 - k point direction in the case choice=5
!!       Compatible only with choice=2,3,5; if idir=0, all derivatives are computed
!!  natom
!!  Cryst
!!  [sorted]=Logical flags defining if the output Cprj has to be sorted by atom type or not.
!!    By default, Cprj matrix elements are unsorted.
!!
!! OUTPUT
!!  cwaveprj
!!
!! PARENTS
!!      classify_bands,m_wfd
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_ug2cprj(Wfd,band,ik_ibz,spin,choice,idir,natom,Cryst,cwaveprj,sorted)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: choice,idir,natom,band,ik_ibz,spin
 logical,optional,intent(in) :: sorted
 class(wfd_t),target,intent(inout) :: Wfd
 type(crystal_t),intent(in) :: Cryst
!arrays
 type(pawcprj_type),intent(inout) :: cwaveprj(natom,Wfd%nspinor)

!Local variables-------------------------------
!scalars
 integer :: cpopt,istwf_k,npw_k,nkpg
 integer :: ia,iatm,dimffnl,itypat,iatom,isp
 character(len=500) :: msg
 type(wave_t),pointer :: wave
 logical :: want_sorted
!arrays
 integer,ABI_CONTIGUOUS pointer :: kg_k(:,:)
 integer,allocatable :: dimcprj_srt(:)
 real(dp) :: kpoint(3)
 real(dp),ABI_CONTIGUOUS pointer :: phkxred(:,:)
 real(dp),allocatable :: cwavef(:,:), kpg(:,:)
 !real(dp),allocatable :: ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),ABI_CONTIGUOUS pointer :: ph3d(:,:,:)    ! ph3d(2,npw_k,matblk)
 real(dp),ABI_CONTIGUOUS pointer :: ffnl(:,:,:,:)  ! ffnl(npw_k,dimffnl,lmnmax,ntypat)
 type(pawcprj_type),allocatable :: Cprj_srt(:,:)

! *********************************************************************

 ! Different form factors have to be calculated and stored in Kdata.
 ABI_CHECK(choice == 1, "choice/=1 not coded")

 dimffnl = 1
 npw_k   = Wfd%npwarr(ik_ibz)
 istwf_k = Wfd%istwfk(ik_ibz)
 kpoint  = Wfd%kibz(:,ik_ibz)

 kg_k    => Wfd%Kdata(ik_ibz)%kg_k
 ph3d    => Wfd%Kdata(ik_ibz)%ph3d
 ffnl    => Wfd%Kdata(ik_ibz)%fnl_dir0der0
 phkxred => Wfd%Kdata(ik_ibz)%phkxred

 ! Compute (k+G) vectors
 nkpg=0
 !% if (choice==3.or.choice==2.or.choice==23) nkpg=3*Wfd%nloalg(3)
 !% if (choice==4.or.choice==24) nkpg=9*Wfd%nloalg(3)
 ABI_MALLOC(kpg,(npw_k,nkpg))
 if (nkpg>0) call mkkpg(kg_k,kpg,kpoint,nkpg,npw_k)

 ! Copy wavefunction in G-space
 ABI_MALLOC(cwavef, (2,npw_k*Wfd%nspinor))
 ABI_CHECK(wfd%get_wave_ptr(band, ik_ibz, spin, wave, msg) == 0, msg)
 cwavef(1,:) = DBLE (wave%ug)
 cwavef(2,:) = AIMAG(wave%ug)

 cpopt = 0 ! Nothing is already calculated.

 want_sorted=.FALSE.; if (present(sorted)) want_sorted=sorted

 if (want_sorted) then
   ! Output cprj are sorted.
   call getcprj(choice,cpopt,cwavef,cwaveprj,ffnl,&
     idir,Wfd%indlmn,istwf_k,kg_k,kpg,kpoint,Wfd%lmnmax,Wfd%mgfft,Wfd%MPI_enreg,&
     Cryst%natom,Cryst%nattyp,Wfd%ngfft,Wfd%nloalg,npw_k,Wfd%nspinor,Cryst%ntypat,&
     phkxred,Wfd%ph1d,ph3d,Cryst%ucvol,1)

 else
   ! Output cprj are unsorted.
   ABI_MALLOC(dimcprj_srt,(Cryst%natom))
   ia=0
   do itypat=1,Cryst%ntypat
     dimcprj_srt(ia+1:ia+Cryst%nattyp(itypat))=Wfd%nlmn_type(itypat)
     ia=ia+Cryst%nattyp(itypat)
   end do

   ABI_MALLOC(Cprj_srt,(natom,Wfd%nspinor))
   call pawcprj_alloc(Cprj_srt,0,dimcprj_srt)
   ABI_FREE(dimcprj_srt)

   ! Calculate sorted cprj.
   call getcprj(choice,cpopt,cwavef,Cprj_srt,ffnl,&
    idir,Wfd%indlmn,istwf_k,kg_k,kpg,kpoint,Wfd%lmnmax,Wfd%mgfft,Wfd%MPI_enreg,&
    Cryst%natom,Cryst%nattyp,Wfd%ngfft,Wfd%nloalg,npw_k,Wfd%nspinor,Cryst%ntypat,&
    phkxred,Wfd%ph1d,ph3d,Cryst%ucvol,1)

   ! Reorder cprj (sorted --> unsorted)
   do iatom=1,Cryst%natom
     iatm=Cryst%atindx(iatom)
     do isp=1,Wfd%nspinor
       cwaveprj(iatom,isp)%cp=Cprj_srt(iatm,isp)%cp
     end do
   end do

   call pawcprj_free(Cprj_srt)
   ABI_FREE(Cprj_srt)
 end if

 ABI_FREE(cwavef)
 ABI_FREE(kpg)

end subroutine wfd_ug2cprj
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wave_init
!! NAME
!!  wave_init
!!
!! FUNCTION
!!   Main creation method for the wave_t data type
!!
!! INPUTS
!!  usepaw=1 if PAW is used.
!!  npw =Number of plane-waves for ug
!!  nfft=Number of FFT points for the real space wavefunction.
!!  nspinor=Number of spinor components.
!!  natom=Number of atoms in cprj matrix elements.
!!  nlmn_size(natom)=Number of (n,l,m) channel for each atom. Ordering of atoms depends on cprj_order
!!  cprj_order=Flag defining the ordering of the atoms in the cprj matrix elements (CPR_RANDOM|CPR_SORTED).
!!    Use to know if we have to reorder the matrix elements when wfd_get_cprj is called.
!!
!! OUTPUT
!!  Wave<wave_t>=The structure fully initialized.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wave_init(Wave, usepaw, npw, nfft, nspinor, natom, nlmn_size, cprj_order)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw,nfft,nspinor,usepaw,natom
 integer(c_int8_t),intent(in) :: cprj_order
 type(wave_t),intent(inout) :: Wave
!arrays
 integer,intent(in) :: nlmn_size(:)

!Local variables ------------------------------
!scalars
 integer,parameter :: ncpgr0=0  ! For the time being, no derivatives

!************************************************************************

 !@wave_t
 if (npw >0) then
   ABI_MALLOC(Wave%ug, (npw*nspinor))
   Wave%has_ug = WFD_ALLOCATED
   Wave%ug = huge(one_gw)
   if (usepaw == 1) then
     ABI_MALLOC(Wave%Cprj, (natom,nspinor))
     call pawcprj_alloc(Wave%Cprj,ncpgr0,nlmn_size)
     Wave%has_cprj = WFD_ALLOCATED
     Wave%cprj_order = cprj_order
   end if
 end if

 if (nfft > 0) then
   ABI_MALLOC(Wave%ur, (nfft*nspinor))
   Wave%ur = huge(one_gw)
   Wave%has_ur = WFD_ALLOCATED
 end if

end subroutine wave_init
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wave_free
!! NAME
!!  wave_free
!!
!! FUNCTION
!!  Main destruction method for the wave_t datatype.
!!
!! INPUTS
!!  [what]=String defining what has to be freed.
!!     "A" =Both ug and ur and Cprj. Default.
!!     "G" =Only ug.
!!     "R" =Only ur
!!     "C" =Only PAW Cprj.
!!
!! SIDE EFFECTS
!!  Memory in Wave is deallocated depending on what
!!
!! PARENTS
!!      m_wfd
!!
!! CHILDREN
!!
!! SOURCE

subroutine wave_free(Wave, what)

!Arguments ------------------------------------
!scalars
 class(wave_t),intent(inout) :: Wave
 character(len=*),optional,intent(in) :: what

!Local variables ------------------------------
!scalars
 character(len=10) :: my_what

!************************************************************************

 my_what="ALL"; if (present(what)) my_what=toupper(what)

 if (.not.firstchar(my_what, ["A", "G", "R", "C"] )) then
   MSG_ERROR(sjoin("Unknow what:", what))
 end if

 if (firstchar(my_what, ["A", "G"])) then
    ABI_SFREE(Wave%ug)
   Wave%has_ug = WFD_NOWAVE
 end if

 if (firstchar(my_what, ["A", "R"])) then
   ABI_SFREE(Wave%ur)
   Wave%has_ur = WFD_NOWAVE
 end if

 if (firstchar(my_what, ["A", "C"])) then
   if (allocated(Wave%Cprj)) then
     call pawcprj_free(Wave%Cprj)
     ABI_FREE(Wave%Cprj)
   end if
   Wave%has_cprj = WFD_NOWAVE
 end if

end subroutine wave_free
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wave_copy
!! NAME
!!  wave_copy
!!
!! FUNCTION
!!  Copy method for the wave_t datatype.
!!
!! PARENTS
!!      m_wfd
!!
!! CHILDREN
!!
!! SOURCE

type(wave_t) function wave_copy(Wave_in) result(Wave_out)

!Arguments ------------------------------------
!scalars
 class(wave_t),intent(in) :: Wave_in

!Local variables ------------------------------
 integer :: natom,nspinor,iatom,ispinor

!************************************************************************

 Wave_out%has_ug = Wave_in%has_ug
 Wave_out%has_ur = Wave_in%has_ur
 Wave_out%has_cprj = Wave_in%has_cprj
 Wave_out%cprj_order = Wave_in%cprj_order

 ABI_MALLOC(Wave_out%ug, (SIZE(Wave_in%ug)))
 Wave_out%ug = Wave_in%ug
 ABI_MALLOC(Wave_out%ur, (SIZE(Wave_in%ur)))
 Wave_out%ur = Wave_in%ur

 natom   = size(Wave_in%Cprj,dim=1)
 nspinor = size(Wave_in%Cprj,dim=2)
 if ((size(Wave_out%Cprj,dim=1) .ne. natom) .or. (size(Wave_out%Cprj,dim=2) .ne. nspinor)) then
   if (allocated(Wave_out%Cprj))  then
     ABI_FREE(Wave_out%Cprj)
   end if
   ABI_MALLOC(Wave_out%Cprj,(natom,nspinor))
 end if

 do ispinor=1,nspinor
   do iatom=1,natom
    Wave_out%Cprj(iatom,ispinor)%ncpgr=Wave_in%Cprj(iatom,ispinor)%ncpgr
    Wave_out%Cprj(iatom,ispinor)%nlmn=Wave_in%Cprj(iatom,ispinor)%nlmn
    call alloc_copy(Wave_in%Cprj(iatom,ispinor)%cp,Wave_out%Cprj(iatom,ispinor)%cp)
    call alloc_copy(Wave_in%Cprj(iatom,ispinor)%dcp,Wave_out%Cprj(iatom,ispinor)%dcp)
   end do
 end do

end function wave_copy
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_get_wave_prt
!! NAME
!!  wfd_get_wave_prt
!!
!! FUNCTION
!!  Return pointer to the wave object corresponding to the given (band, ik_ibz, spin) indices.
!!  If the state is not treated ...
!!
!! INPUTS
!!   band=Band index.
!!   ik_ibz=k-point index
!!   spin=Spin index.
!!
!! PARENTS
!!      m_wfd
!!
!! CHILDREN
!!
!! SOURCE

integer function wfd_get_wave_ptr(wfd, band, ik_ibz, spin, wave_ptr, msg) result(ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz, spin, band
 class(wfd_t),target,intent(in) :: wfd
 type(wave_t),pointer :: wave_ptr
 character(len=*),intent(out) :: msg

!Local variables ------------------------------
!scalars
 integer :: ib, ik, is

!************************************************************************

 ierr = 1
 if (any(wfd%bks2wfd(:, band, ik_ibz, spin) == 0)) then
   write(msg,'(a,i0,a,3(i0,1x))')" MPI rank ",Wfd%my_rank," doesn't have ug for (band, ik_ibz, spin): ",band,ik_ibz,spin
   wave_ptr => null(); return
 end if

 ib = wfd%bks2wfd(1, band, ik_ibz, spin)
 ik = wfd%bks2wfd(2, band, ik_ibz, spin)
 is = wfd%bks2wfd(3, band, ik_ibz, spin)
 wave_ptr => wfd%s(is)%k(ik)%b(ib)
 !if (wave_ptr%has_ug /= WFD_STORED)

 ierr = 0

end function wfd_get_wave_ptr
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_push_ug
!! NAME
!!  wfd_push_ug
!!
!! FUNCTION
!!  This routine changes the status of the object by saving the wavefunction in the correct
!!  slot inside Wfd%Wave. It also set the corresponding has_ug flag to WFD_STORED.
!!  If the status of the corresponding ur is (WFD_STORED|WFD_ALLOCATED), then an G->R FFT transform
!!  is done (see also update_ur)
!!
!! INPUTS
!!   band=Band index.
!!   ik_ibz=k-point index
!!   spin=Spin index.
!!   Cryst<crystal_t>=Object defining the unit cell and its symmetries.
!!   ug(npw_k*Wfd%nspinor)=The ug to be saved.
!!   [update_ur]=If .FALSE.: no G-->R transform is done even if ur is (WFD_STORED|WFD_ALLOCATED) so be careful.
!!               Defaults to .TRUE.
!!   [update_cprj]=If .FALSE.: <C|p_i> matrix elements are not recalculatedd even
!!     if cprj is (WFD_STORED|WFD_ALLOCATED) so be careful. Defaults to .TRUE.
!!
!! SIDE EFFECTS
!!   Wfd<wfd_t>=See above.
!!
!! PARENTS
!!      m_wfd
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_push_ug(Wfd, band, ik_ibz, spin, Cryst, ug, update_ur, update_cprj)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,spin,band
 logical,optional,intent(in) :: update_ur,update_cprj
 class(wfd_t),target,intent(inout) :: Wfd
 type(crystal_t),intent(in) :: Cryst
!arrays
 complex(gwpc),intent(inout) :: ug(:)

!Local variables ------------------------------
!scalars
 integer,parameter :: choice1=1,idir0=0,tim_fourdp=5,ndat1=1
 integer :: npw_k, ib, ik, is
 logical :: do_update_ur,do_update_cprj,want_sorted
 character(len=500) :: msg
 type(wave_t),pointer :: wave

!************************************************************************

 if (size(ug) /= Wfd%npwarr(ik_ibz) * Wfd%nspinor) then
   MSG_ERROR("Wrong size in assumed shape array")
 end if

 if (any(wfd%bks2wfd(:, band, ik_ibz, spin) == 0)) then
   write(msg,'(a,i0,a,3(i0,1x))')" MPI rank ",Wfd%my_rank," doesn't have ug for (band, ik_ibz, spin): ",band,ik_ibz,spin
   MSG_ERROR(msg)
 end if

 ib = wfd%bks2wfd(1, band, ik_ibz, spin)
 ik = wfd%bks2wfd(2, band, ik_ibz, spin)
 is = wfd%bks2wfd(3, band, ik_ibz, spin)

 wave => wfd%s(is)%k(ik)%b(ib)
 wave%ug = ug
 wave%has_ug = WFD_STORED

 if (Wfd%debug_level>0) then
   if (wave%has_ug == WFD_NOWAVE) then
     write(msg,'(a,i0,a,3(i0,1x))')" MPI rank ",Wfd%my_rank," doesn't have ug for (band, ik_ibz, spin): ",band,ik_ibz,spin
     MSG_ERROR(msg)
   end if
 end if

 if (Wfd%usepaw==1) then
   ! Update the corresponding cprj if required.
   do_update_cprj=.TRUE.; if (present(update_cprj)) do_update_cprj=update_cprj
   if (do_update_cprj) then
     want_sorted = (wave%cprj_order == CPR_SORTED)
     call wfd%ug2cprj(band, ik_ibz, spin, choice1, idir0, wfd%natom, cryst, wave%cprj, sorted=want_sorted)
     wave%has_cprj = WFD_STORED
   else
     wave%has_cprj = WFD_ALLOCATED
   end if
 end if

 if (any(wave%has_ur == [WFD_STORED, WFD_ALLOCATED])) then
   ! Update the corresponding ur if required.
   do_update_ur=.TRUE.; if (present(update_ur)) do_update_ur=update_ur

   if (do_update_ur) then
     npw_k = Wfd%npwarr(ik_ibz)
     call fft_ug(npw_k,Wfd%nfft,Wfd%nspinor,ndat1,Wfd%mgfft,Wfd%ngfft,Wfd%istwfk(ik_ibz),&
       Wfd%Kdata(ik_ibz)%kg_k,Wfd%Kdata(ik_ibz)%gbound,ug,wave%ur)
     wave%has_ur = WFD_STORED
   else
     wave%has_ur = WFD_ALLOCATED
   end if
 end if

end subroutine wfd_push_ug
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_extract_cgblock
!! NAME
!!  wfd_extract_cgblock
!!
!! FUNCTION
!!  This routine extract a block of wavefunctions for a given spin and k-points.
!!  The wavefunctions are stored in a real(dp) array with the same convention
!!  as the one used in the GS part of Abinit, i.e cg_block(2,nspinor*npw_k*num_bands)
!!
!! INPUTS
!!   Wfd<wfd_t>=Wavefunction descriptor.
!!   band_list(:)=List of bands to extract
!!   ik_ibz=k-point index
!!   spin=Spin index.
!!
!! OUTPUT
!!   cgblock(nspinor*npw_k*num_bands)=A contiguous block of memory with the set of u(g)
!!
!! PARENTS
!!      m_wfd
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_extract_cgblock(Wfd,band_list,ik_ibz,spin,cgblock)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,spin
 class(wfd_t),intent(in) :: Wfd
!arrays
 integer,intent(in) :: band_list(:)
 real(dp),intent(out) :: cgblock(:,:)

!Local variables ------------------------------
!scalars
 integer :: ii,band,start,istop,npw_k
 character(len=500) :: msg
 type(wave_t),pointer :: wave

!************************************************************************

 npw_k = Wfd%npwarr(ik_ibz)

 if (size(cgblock, dim=1) /= 2) then
   MSG_ERROR("Wrong size(1) in assumed shape array")
 end if

 if (size(cgblock, dim=2) /= Wfd%nspinor* npw_k * size(band_list)) then
   MSG_ERROR("Wrong size in assumed shape array")
 end if

 start = 1
 do ii=1,size(band_list)
   band = band_list(ii)
   ABI_CHECK(wfd%get_wave_ptr(band, ik_ibz, spin, wave, msg) == 0, msg)
   if (wave%has_ug /= WFD_STORED) then
     write(msg,"(3(a,i0),a)")"u(g) for band: ",band,", ik_ibz: ",ik_ibz,", spin: ",spin," is not stored!"
     MSG_ERROR(msg)
   end if
   istop = start + Wfd%nspinor*npw_k - 1
   cgblock(1,start:istop) = REAL(wave%ug)
   cgblock(2,start:istop) = AIMAG(wave%ug)
   start = start + Wfd%nspinor * npw_k
 end do

end subroutine wfd_extract_cgblock
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_rank_has_ug
!! NAME
!!  wfd_rank_has_ug
!!
!! FUNCTION
!!  This function is used to ask a particular processor whether it has a particular ug and with which status.
!!
!! INPUTS
!!   rank=The MPI rank of the processor.
!!   band=Band index.
!!   ik_ibz=k-point index
!!   spin=Spin index.
!!
!! NOTES
!!   A zero index can be used to inquire the status of a bunch of states.
!!   Thus (band,ik_ibz,spin) = (0,1,1) means: Do you have at least one band for the first k-point and the first spin.
!!
!! PARENTS
!!
!! SOURCE

function wfd_rank_has_ug(Wfd,rank,band,ik_ibz,spin)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin,rank
 logical :: wfd_rank_has_ug
 class(wfd_t),intent(in) :: Wfd

!Local variables ------------------------------
!scalars
 integer :: nzeros
 integer(c_int8_t) :: bks_flag
!arrays
 integer :: indices(3)

!************************************************************************

 indices = [band,ik_ibz,spin]
 bks_flag = WFD_STORED

 if (ALL(indices/= [0,0,0])) then
   wfd_rank_has_ug = (Wfd%bks_tab(band,ik_ibz,spin,rank) == bks_flag); RETURN
 else
   nzeros = COUNT(indices==0)
   if (nzeros==3) MSG_ERROR("All indices are zero!")

   if (band==0) then
     if (nzeros==1) wfd_rank_has_ug = ANY( Wfd%bks_tab(:,ik_ibz,spin,rank)==bks_flag)
     if (ik_ibz==0) wfd_rank_has_ug = ANY( Wfd%bks_tab(:,:,spin,rank)     ==bks_flag)
     if (spin  ==0) wfd_rank_has_ug = ANY( Wfd%bks_tab(:,ik_ibz,:,rank)   ==bks_flag)

   else if (ik_ibz==0) then
     if (nzeros==1) wfd_rank_has_ug = ANY( Wfd%bks_tab(band,:,spin,rank)==bks_flag)
     if (band  ==0) wfd_rank_has_ug = ANY( Wfd%bks_tab(:,:,spin,rank)   ==bks_flag)
     if (spin  ==0) wfd_rank_has_ug = ANY( Wfd%bks_tab(band,:,:,rank)   ==bks_flag)

   else
     if (nzeros==1) wfd_rank_has_ug = ANY( Wfd%bks_tab(band,ik_ibz,:,rank)==bks_flag)
     if (ik_ibz==0) wfd_rank_has_ug = ANY( Wfd%bks_tab(band,:,:,rank)     ==bks_flag)
     if (band  ==0) wfd_rank_has_ug = ANY( Wfd%bks_tab(:,ik_ibz,:,rank)   ==bks_flag)
   end if
 end if

end function wfd_rank_has_ug
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_ihave_ug
!! NAME
!!  wfd_ihave_ug
!!
!! FUNCTION
!!  This function is used to ask the processor whether it has a particular ug and with which status.
!!
!! INPUTS
!!   band=Band index.
!!   ik_ibz=k-point index
!!   spin=Spin index.
!!   [how]=string defining which status is checked.
!!     Possible mutually exclusive values: "Allocated", "Stored".
!!     Only the first character is checked (no case-sensitive)
!!     By default the function returns .TRUE. if the wave is either WFD_ALLOCATED or WFD_STORED.
!!
!! NOTES
!!   A zero index can be used to inquire the status of a bunch of states.
!!   Thus (band,ik_ibz,spin) = (0,1,1) means: Do you have at least one band for the first k-point and the first spin.
!!
!! PARENTS
!!
!! SOURCE

pure function wfd_ihave_ug(Wfd, band, ik_ibz, spin, how)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin
 logical :: wfd_ihave_ug
 character(len=*),optional,intent(in) :: how
 class(wfd_t),intent(in) :: Wfd

!Local variables ------------------------------
!scalars
 integer :: ib, ik, is
 integer(c_int8_t) :: check2(2)

!************************************************************************

 check2 = [WFD_ALLOCATED, WFD_STORED]
 if (present(how)) then
   if (firstchar(how, ["A", "a"])) check2 = [WFD_ALLOCATED, WFD_ALLOCATED]
   if (firstchar(how, ["S", "s"])) check2 = [WFD_STORED, WFD_STORED]
 end if
 ib = wfd%bks2wfd(1, band, ik_ibz, spin)
 ik = wfd%bks2wfd(2, band, ik_ibz, spin)
 is = wfd%bks2wfd(3, band, ik_ibz, spin)
 wfd_ihave_ug = .False.
 if (ib /= 0) wfd_ihave_ug = any(wfd%s(is)%k(ik)%b(ib)%has_ug == check2)

end function wfd_ihave_ug
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_mybands
!! NAME
!!  wfd_mybands
!!
!! FUNCTION
!!  Return the list of band indices of the ug owned by this node at given (k,s).
!!
!! INPUTS
!!  ik_ibz=Index of the k-point in the IBZ
!!  spin=spin index
!!  [how]=string defining which status is checked.
!!    Possible mutually exclusive values: "Allocated", "Stored".
!!    Only the first character is checked (no case-sensitive)
!!    By default the list of bands whose status is either WFD_ALLOCATED or WFD_STORED is returned.
!!
!! OUTPUT
!!  how_manyb=The number of bands owned by this node
!!  my_band_list(Wfd%mband)=The first how_manyb values are the bands treated by this node.
!!
!! PARENTS
!!      m_wfd
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_mybands(Wfd, ik_ibz, spin, how_manyb, my_band_list, how)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,spin
 integer,intent(out) :: how_manyb
 character(len=*),optional,intent(in) :: how
 class(wfd_t),intent(in) :: Wfd
!arrays
 integer,intent(out) :: my_band_list(Wfd%mband)

!Local variables ------------------------------
!scalars
 integer :: band
 logical :: do_have

!************************************************************************

 how_manyb=0; my_band_list=-1
 do band=1,Wfd%nband(ik_ibz,spin)
   if (present(how)) then
     do_have = wfd%ihave_ug(band, ik_ibz, spin, how=how)
   else
     do_have = wfd%ihave_ug(band, ik_ibz, spin)
   end if
   if (do_have) then
     how_manyb = how_manyb + 1
     my_band_list(how_manyb) = band
   end if
 end do

end subroutine wfd_mybands
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_show_bkstab
!! NAME
!!  wfd_show_bkstab
!!
!! FUNCTION
!!  Print a table showing the distribution of the wavefunctions.
!!
!! PARENTS
!!      m_wfd
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_show_bkstab(Wfd, unit)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: unit
 class(wfd_t),intent(in) :: Wfd

!Local variables ------------------------------
!scalars
 integer :: ik_ibz,spin,band,nband_k,width
 character(len=1) :: chlist(0:Wfd%nproc-1)
 character(len=500) :: fmt

!************************************************************************

 width = max(80, Wfd%nproc)

 write(fmt,"(a,i0,a)")"(i5,3x,",Wfd%nproc,"(a))"

 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz
     write(unit,"(a)")repeat("=",width)
     write(unit,"(2(a,i0))")"Spin: ",spin,", ik_ibz: ",ik_ibz
     write(unit,"(a)")"MPI rank ----> (A=allocated, S=Stored, N=NoWave)."
     nband_k = Wfd%nband(ik_ibz, spin)
     do band=1,nband_k
       where (Wfd%bks_tab(band, ik_ibz, spin,:) == WFD_NOWAVE)
         chlist = "N"
       elsewhere (Wfd%bks_tab(band, ik_ibz, spin,:) == WFD_ALLOCATED)
         chlist = "A"
       elsewhere (Wfd%bks_tab(band, ik_ibz, spin,:) == WFD_STORED)
         chlist = "S"
       end where
       write(unit,fmt)band,chlist(:)
     end do
     write(unit,"(a)")repeat("=",width)
   end do
 end do

end subroutine wfd_show_bkstab
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_bands_of_rank
!! NAME
!!  wfd_bands_of_rank
!!
!! FUNCTION
!!  Return the list of band index of the ug owned by a given processor at given (k,s).
!!
!! INPUTS
!!  Wfd
!!  rank=The MPI rank of the processor.
!!  ik_ibz=Index of the k-point in the IBZ
!!  spin=spin index
!!
!! OUTPUT
!!  how_manyb=The number of bands owned by this node
!!  rank_band_list(Wfd%mband)=The first how_manyb values are the bands treated by the node.
!!
!! PARENTS
!!      m_wfd
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_bands_of_rank(Wfd,rank,ik_ibz,spin,how_manyb,rank_band_list)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,spin,rank
 integer,intent(out) :: how_manyb
 class(wfd_t),intent(in) :: Wfd
!arrays
 integer,intent(out) :: rank_band_list(Wfd%mband)

!Local variables ------------------------------
!scalars
 integer :: band
 logical :: it_has

!************************************************************************

 how_manyb=0; rank_band_list=-1
 do band=1,Wfd%nband(ik_ibz,spin)
   it_has = wfd_rank_has_ug(Wfd, rank, band, ik_ibz, spin)
   if (it_has) then
     how_manyb = how_manyb +1
     rank_band_list(how_manyb)=band
   end if
 end do

end subroutine wfd_bands_of_rank
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_get_ug
!! NAME
!!  wfd_get_ug
!!
!! FUNCTION
!!  Get a **copy** of a wave function in G-space.
!!
!! INPUTS
!!  Wfd<wfd_t>=the data type
!!  band=the index of the band.
!!  ik_ibz=Index of the k-point in the IBZ
!!  spin=spin index
!!
!! OUTPUT
!!  ug(npw_k*Wfd%nspinor)=The required wavefunction in G-space
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_get_ug(Wfd, band, ik_ibz, spin, ug)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin
 class(wfd_t),intent(inout) :: Wfd
!arrays
 complex(gwpc),intent(out) :: ug(Wfd%npwarr(ik_ibz)*Wfd%nspinor)

!Local variables ------------------------------
!scalars
 integer :: npw_k
 character(len=500) :: msg
 type(wave_t),pointer :: wave
!************************************************************************

 ABI_CHECK(wfd%get_wave_ptr(band, ik_ibz, spin, wave, msg) == 0, msg)

 if (.not. wave%has_ug == WFD_STORED) then
   write(msg,'(a,i0,a,3i0)')" Node ",Wfd%my_rank," doesn't have (band,ik_ibz,spin): ",band,ik_ibz,spin
   MSG_BUG(msg)
 end if

 npw_k = Wfd%npwarr(ik_ibz)
 call xcopy(npw_k*Wfd%nspinor, wave%ug, 1, ug, 1)

end subroutine wfd_get_ug
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_wave_free
!! NAME
!!  wfd_wave_free
!!
!! FUNCTION
!!  Collection procedure that frees the set of waves specified by mask.
!!
!! INPUTS
!!  mask(Wfd%mband,Wfd%nkibz,Wfd%nsppol)=.TRUE. if the memory allocated for
!!    this state has to be freed
!!  [what]=String specifying which array have to be deallocated.
!!    Possible values (no case-sensitive).
!!      "All"= To free both ug and ur and PAW Cprj, if any. Default
!!      "G"  = Only ug
!!      "R"  = Only ur.
!!      "C"  = Only PAW Cprj.
!!
!! SIDE EFFECTS
!!  Wfd<wfd_t>=See above.
!!
!! PARENTS
!!      bethe_salpeter,m_haydock
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_wave_free(Wfd, what, bks_mask)

!Arguments ------------------------------------
!scalars
 class(wfd_t),target,intent(inout) :: Wfd
 character(len=*),optional,intent(in) :: what
!arrays
 logical,optional,intent(in) :: bks_mask(Wfd%mband,Wfd%nkibz,Wfd%nsppol)

!Local variables ------------------------------
!scalars
 integer :: ik_ibz, spin, band, ib, ik, is
 logical :: do_free
 type(wave_t),pointer :: wave
 !character(len=500) :: msg
 character(len=10) :: my_what
!************************************************************************

 my_what="ALL"; if (present(what)) my_what=toupper(what)

 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz
     do band=1,Wfd%nband(ik_ibz,spin)
        do_free=.TRUE.; if (present(bks_mask)) do_free=bks_mask(band,ik_ibz,spin)
        if (do_free) then
          ib = wfd%bks2wfd(1, band, ik_ibz, spin)
          ik = wfd%bks2wfd(2, band, ik_ibz, spin)
          is = wfd%bks2wfd(3, band, ik_ibz, spin)
          if (ib /= 0) then
            wave => wfd%s(is)%k(ik)%b(ib)
            call wave%free(what=my_what)
          end if
          ! Update the associated flags if we release the G-space.
          if ( firstchar(my_what, ["A", "G"])) Wfd%bks_tab(band, ik_ibz, spin, Wfd%my_rank) = WFD_NOWAVE
        end if
     end do
   end do
 end do

 ! Reinit the MPI communicators.
 call wfd%set_mpicomm()

end subroutine wfd_wave_free
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_who_has_ug
!! NAME
!!  wfd_who_has_ug
!!
!! FUNCTION
!!  Return the number of processors having a particular (b,k,s) state as well as their MPI rank.
!!  Warning: Wfd%bks_tab is supposed to be up-to-date (see wfd_update_bkstab).
!!
!! INPUTS
!!  band=the index of the band.
!!  ik_ibz=Index of the k-point in the IBZ
!!  spin=spin index
!!
!! OUTPUT
!!  how_many=The number of nodes owing this ug state.
!!  proc_ranks(1:how_many)=Gives the MPI rank of the nodes owing the state.
!!
!! PARENTS
!!      m_wfd
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_who_has_ug(Wfd,band,ik_ibz,spin,how_many,proc_ranks)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin
 integer,intent(out) :: how_many
 class(wfd_t),intent(in) :: Wfd
!arrays
 integer,intent(out) :: proc_ranks(Wfd%nproc)

!Local variables ------------------------------
!scalars
 integer :: irank
 logical :: bks_select,spin_select,kpt_select
 character(len=500) :: msg

!************************************************************************

 bks_select  = (band/=0.and.ik_ibz/=0.and.spin/=0)
 spin_select = (band==0.and.ik_ibz==0.and.spin/=0)
 kpt_select = (band==0.and.ik_ibz/=0.and.spin/=0)

 how_many=0; proc_ranks=-1

 if (bks_select) then
   ! List the proc owining this (b,k,s) state.
   do irank=0,Wfd%nproc-1
     if (Wfd%bks_tab(band, ik_ibz, spin, irank) == WFD_STORED) then
       how_many = how_many +1
       proc_ranks(how_many)=irank
     end if
   end do

 else if (spin_select) then
   ! List the proc owining at least one state with this spin.
   do irank=0,Wfd%nproc-1
     if ( ANY(Wfd%bks_tab(:,:,spin,irank)==WFD_STORED) ) then
       how_many = how_many +1
       proc_ranks(how_many)=irank
     end if
   end do

 else if (kpt_select) then
   ! List the proc owining at least one state with this (k-point, spin).
   do irank=0,Wfd%nproc-1
     if ( ANY(Wfd%bks_tab(:,ik_ibz,spin,irank)==WFD_STORED) ) then
       how_many = how_many +1
       proc_ranks(how_many)=irank
     end if
   end do

 else
   write(msg,'(a,3(i0,1x))')" Wrong value for (b,k,s): ",band,ik_ibz,spin
   MSG_ERROR(msg)
 end if

end subroutine wfd_who_has_ug
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_update_bkstab
!! NAME
!!  wfd_update_bkstab
!!
!! FUNCTION
!!  This routine should be called by all the nodes before any MPI operation involving the object.
!!  It updates the bks_tab storing information on the distribution of ug.
!!
!! INPUT
!!  [show]=If present and > 0, print tabs to unit show.
!!
!! SIDE EFFECTS
!!  Wfd%bks_tab
!!
!! PARENTS
!!      m_sigma,m_wfd,wfd_mkrho
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_update_bkstab(Wfd, show)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: show
 class(wfd_t),intent(inout) :: Wfd

!Local variables ------------------------------
!scalars
 integer :: ierr, nelem, spin, ik_ibz, band, is, ik, ib
 integer(c_int8_t),allocatable :: my_vtab(:),gather_vtabs(:)
 !logical,allocatable :: tab_ranks(:)

!************************************************************************

 ! Fill my slice of the global table.
 do spin=1,wfd%nsppol
   do ik_ibz=1,wfd%nkibz
     do band=1,Wfd%nband(ik_ibz, spin)
       ib = wfd%bks2wfd(1, band, ik_ibz, spin)
       ik = wfd%bks2wfd(2, band, ik_ibz, spin)
       is = wfd%bks2wfd(3, band, ik_ibz, spin)
       if (ib /= 0) then
         Wfd%bks_tab(band, ik_ibz, spin, Wfd%my_rank) = wfd%s(is)%k(ik)%b(ib)%has_ug
       else
         Wfd%bks_tab(band, ik_ibz, spin, Wfd%my_rank) = WFD_NOWAVE
       end if
     end do
   end do
 end do

 ! Gather flags on each node.
 nelem = Wfd%mband*Wfd%nkibz*Wfd%nsppol
 ABI_MALLOC(my_vtab, (nelem))
 my_vtab(:) = reshape(Wfd%bks_tab(:,:,:,Wfd%my_rank), [nelem])

 ABI_MALLOC(gather_vtabs, (nelem*Wfd%nproc))

 call xmpi_allgather(my_vtab,nelem,gather_vtabs,Wfd%comm,ierr)

 Wfd%bks_tab(:,:,:,:) = reshape(gather_vtabs, [Wfd%mband, Wfd%nkibz, Wfd%nsppol, Wfd%nproc])
 ABI_FREE(my_vtab)
 ABI_FREE(gather_vtabs)

#if 0
 ! This is gonna be slow but if lot of k-points as I cannot assume bands or k-points have been filtered
 ! Need to introduce global_filter_ikibz_spin in wfd_init ...
 ABI_MALLOC(tab_ranks, (wfd%nproc))
 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz
     !if wfd%global_filter_ikibz_spin(ik_ibz, spin) cycle
     do band=1,Wfd%nband(ik_ibz, spin)
       tab_ranks = .False.
       if (len(wfd%bks_ranks(band, ik_ibz, spin) > 0) then
         if (any(wfd%my_rank == wfd%bks_ranks(band, ik_ibz, spin)) tab_ranks(wfd%my_rank) = .True.
       end if
       call xmpi_lor(tab_ranks, wfd%comm)
       call bool2index(tab_ranks, wfd%bks_ranks(band, ik_ibz, spin))
     end do
   end do
 end do
 ABI_FREE(tab_ranks)
#endif

 if (present(show)) then
   if (show >= 0) call wfd%show_bkstab(unit=show)
 end if

end subroutine wfd_update_bkstab
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_set_mpicomm
!! NAME
!!  wfd_set_mpicomm
!!
!! FUNCTION
!!
!! PARENTS
!!      m_wfd
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_set_mpicomm(Wfd)

!Arguments ------------------------------------
!scalars
 class(wfd_t),intent(inout) :: Wfd

!Local variables ------------------------------
!scalars
 integer :: spin,ierr,how_many,spin_comm
 integer :: world_group,spin_group
 !character(len=500) :: msg
!arrays
 integer :: proc_ranks(Wfd%nproc)

!************************************************************************

 ! First free the old communicators.
 call xmpi_comm_free(Wfd%bks_comm)

 ! Update the bks_tab.
 call wfd%update_bkstab()

 call xmpi_comm_group(Wfd%comm,world_group,ierr)

 ! Init spin communicators.
 do spin=1,Wfd%nsppol
   ! The list of procs owining at least one state with this spin.
   call wfd_who_has_ug(Wfd,0,0,spin,how_many,proc_ranks)

   if (how_many>0) then
     call xmpi_group_incl(world_group,how_many,proc_ranks,spin_group,ierr)
     call xmpi_comm_create(Wfd%comm,spin_group,spin_comm,ierr)
     Wfd%bks_comm(0,0,spin) = spin_comm
     call xmpi_group_free(spin_group)
   else
     MSG_WARNING(sjoin("Nobody has spin:",itoa(spin)))
     Wfd%bks_comm(0,0,spin) = xmpi_comm_null
   end if

 end do

 call xmpi_group_free(world_group)

end subroutine wfd_set_mpicomm
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_distribute_bands
!! NAME
!!  wfd_distribute_bands
!!
!! FUNCTION
!!  Distribute a set of bands taking into account the distribution of the ug.
!!
!! INPUTS
!!  band=the index of the band.
!!  ik_ibz=Index of the k-point in the IBZ
!!  spin=spin index
!!  [got(Wfd%nproc)]=The number of tasks already assigned to the nodes.
!!  [bmask(Wfd%mband)]=The routine will raise an error if one band index
!!    is not treated by any processor. bmask can be used to select the subset of
!!    indices that are expected to be available.
!!
!! OUTPUT
!!   my_nband=The number of bands that will be treated by this node.
!!   my_band_list(1:my_nband)=The band indices for this node
!!
!! PARENTS
!!      cchi0q0_intraband,m_sigma,m_wfd,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_distribute_bands(Wfd,ik_ibz,spin,my_nband,my_band_list,got,bmask)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,spin
 integer,intent(out) :: my_nband
 class(wfd_t),intent(in) :: Wfd
!arrays
 integer,intent(out) :: my_band_list(Wfd%mband)
 integer,optional,intent(inout) :: got(Wfd%nproc)
 logical,optional,intent(in) :: bmask(Wfd%mband)

!Local variables ------------------------------
!scalars
 integer :: band,how_many,idle
 character(len=500) :: msg
!arrays
 integer :: proc_ranks(Wfd%nproc),get_more(Wfd%nproc)
 logical :: rank_mask(Wfd%nproc)

!************************************************************************

 my_nband=0; my_band_list=0
 get_more=0; if (present(got)) get_more = got

 do band=1,Wfd%nband(ik_ibz,spin)
   if (present(bmask)) then
     if (.not.bmask(band)) CYCLE
   end if

   call wfd_who_has_ug(Wfd, band, ik_ibz, spin, how_many, proc_ranks)

   if (how_many == 1) then
     ! I am the only one owing this band. Add it to list.
     if (proc_ranks(1) == Wfd%my_rank) then
       my_nband=my_nband + 1
       my_band_list(my_nband) = band
     end if
   else if (how_many > 1) then
     ! This band is duplicated. Assign it trying to obtain a good load distribution.
     rank_mask=.FALSE.; rank_mask(proc_ranks(1:how_many)+1)=.TRUE.
     idle = imin_loc(get_more,mask=rank_mask)
     get_more(idle) = get_more(idle) + 1
     if (Wfd%my_rank==idle-1) then
       my_nband=my_nband + 1
       my_band_list(my_nband) = band
     end if
   else
     write(msg,'(a,3(i0,1x))')" No processor has (band, ik_ibz, spin): ",band,ik_ibz,spin
     MSG_ERROR(msg)
   end if
 end do

 if (present(got)) got = get_more

end subroutine wfd_distribute_bands
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_rotate
!! NAME
!! wfd_rotate
!!
!! FUNCTION
!!  This routine performs a linear transformation of the wavefunctions stored in Wfd
!!  taking into account memory distribution. The transformation is done in G-space
!!  therefore all the ug should be available. Wavefunctions in real space are then
!!  obtained via FFT. The implementation assumes that the matrix associated to the
!!  linear transformation is sparse (No BLAS-3 calls here).
!!
!! INPUTS
!!  Cryst<crystal_t>=Object defining the unit cell and its symmetries.
!!  m_ks_to_qp(mband,mband,nkibz,nsppol)= expansion of the QP amplitudes in terms of KS wavefunctions.
!!  [bmask(mband,nkibz,nsppol)]=The routine will raise an error if one band index
!!    is not treated by any processor. bmask can be used to select the subset of
!!    indices that are expected to be available.
!!
!! SIDE EFFECTS
!!   Wfd<wfd_t>=See above.
!!
!! PARENTS
!!      bethe_salpeter,screening,sigma
!!
!! CHILDREN
!!
!! SOURCE
!!

subroutine wfd_rotate(Wfd, Cryst, m_ks_to_qp, bmask)

!Arguments ------------------------------------
!scalars
 class(wfd_t),intent(inout) :: Wfd
 type(crystal_t),intent(in) :: Cryst
!arrays
 complex(dpc),target,intent(in) :: m_ks_to_qp(Wfd%mband,Wfd%mband,Wfd%nkibz,Wfd%nsppol)
 logical,optional,intent(in) :: bmask(Wfd%mband,Wfd%nkibz,Wfd%nsppol)

!Local variables-------------------------------
!scalars
 integer :: band,ik_ibz,spin,ierr,icol,nnew,inew,my_nband,ib,npw_k,istwf_k
 character(len=500) :: msg
 type(wave_t),pointer :: wave
!arrays
 integer :: new_list(Wfd%mband),my_band_list(Wfd%mband)
 complex(dpc),ABI_CONTIGUOUS pointer :: umat_sk(:,:)
 complex(gwpc) :: mcol(Wfd%mband)
 complex(gwpc),allocatable :: new_ug(:,:) !, new_ur(:)

!************************************************************************

 DBG_ENTER("COLL")

 ! Update the distribution table, first.
 call wfd%update_bkstab()

 ! Calculate: $\Psi^{QP}_{r,b} = \sum_n \Psi^{KS}_{r,n} M_{n,b}$
 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz
     npw_k  = Wfd%npwarr(ik_ibz)
     istwf_k = Wfd%istwfk(ik_ibz)
     if (istwf_k /= 1) then
       MSG_WARNING("wfd_rotate with istwfk /= 1")
     end if
     umat_sk => m_ks_to_qp(:,:,ik_ibz,spin)

     ! Select only those states that are mixed by the (sparse) m_ks_to_qp.
     nnew=0; new_list=0
     do icol=1,Wfd%nband(ik_ibz,spin)
       mcol = m_ks_to_qp(:,icol,ik_ibz,spin)
       mcol(icol) = mcol(icol) - cone
       if (ANY(ABS(mcol)>tol12)) then  ! Avoid a simple copy.
         nnew=nnew+1
         new_list(nnew)=icol
       end if
     end do
     if (nnew==0) CYCLE ! Nothing to do.

     ! Retrieve the set of band indices that have to be treated by
     ! this node taking into account a possible duplication.
     if (present(bmask)) then
       call wfd%distribute_bands(ik_ibz,spin,my_nband,my_band_list,bmask=bmask(:,ik_ibz,spin))
     else
       call wfd%distribute_bands(ik_ibz,spin,my_nband,my_band_list)
     end if

     !if (my_nband>0) then
     !  write(std_out,*)" At (ik_ibz,spin) ",ik_ibz,spin,&
     !  & ", rank ",Wfd%my_rank," will sum ",my_nband," bands, my_band_list: ",my_band_list(1:my_nband)
     !end if
     ABI_MALLOC(new_ug,(npw_k*Wfd%nspinor,nnew))
     new_ug=czero
     do inew=1,nnew
       icol = new_list(inew)
       do ib=1,my_nband
         band = my_band_list(ib)
         if (ABS(umat_sk(band,icol))>tol12) then
            ABI_CHECK(wfd%get_wave_ptr(band, ik_ibz, spin, wave, msg) == 0, msg)
            new_ug(:,inew) = new_ug(:,inew) + umat_sk(band, icol) * wave%ug
         end if
       end do
     end do

     !if (istwf_k /= 1) then
     !  ABI_MALLOC(new_ur, (wfd%nfft * wfd%nspinor * nnew))
     !  call fft_ug_dpc(npw_k, wfd%nfft, wfd%nspinor, nnew, wfd%mgfft, wfd%ngfft, istwf_k, &
     !                  wfd%kdata(ik_ibz)%kg_k, wfd%kdata(ik_ibz)%gbound, new_ug, new_ur)
     !  new_ur = real(new_ur)
     !  call fft_ur_dpc(npw_k, wfd%nfft, wfd%nspinor, nnew, wfd%mgfft, wfd%ngfft, istwf_k, &
     !                  wfd%kdata(ik_ibz)%kg_k, wfd%kdata(ik_ibz)%gbound, new_ur, new_ug)
     !  ABI_FREE(new_ur)
     !end if

     call xmpi_sum(new_ug,Wfd%comm,ierr)

     ! Update the input wave functions
     do inew=1,nnew
       band = new_list(inew)
       if (wfd%ihave_ug(band, ik_ibz, spin)) call wfd%push_ug(band, ik_ibz, spin, Cryst, new_ug(:,inew))
     end do

     ABI_FREE(new_ug)
   end do !ik_ibz
 end do !spin

 ! Reinit the storage mode of Wfd as ug have been changed.
 ! This is needed only if FFTs are not done in wfd_push_ug. Do not know which one is faster.
 !call wfd%reset_ur_cprj()
 call xmpi_barrier(Wfd%comm)

 DBG_EXIT("COLL")

end subroutine wfd_rotate
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_iterator_bks
!! NAME
!!  wfd_iterator_bks
!!
!! FUNCTION
!!  Iterator used to loop over bands, k-points and spin indices
!!  taking into account the distribution of the ug.
!!
!! INPUTS
!!  Wfd<wfd_t>=
!!  bks_mask(Wfd%mband.Wfd%nkibz,Wfd%nsppol)= mask used to select the (b,k,s) indices.
!!
!! OUTPUT
!!  iter_bks<iter2_t>=Iterator over the bands treated by this node for each k-point and spin.
!!
!! PARENTS
!!
!! SOURCE

type(iter2_t) function wfd_iterator_bks(Wfd, bks_mask) result(iter_bks)

!Arguments ------------------------------------
!scalars
 class(wfd_t),intent(in) :: Wfd
!arrays
 logical,optional,intent(in) :: bks_mask(Wfd%mband,Wfd%nkibz,Wfd%nsppol)

!Local variables ------------------------------
!scalars
 integer :: ik_ibz,spin,my_nband
!arrays
 integer :: my_band_list(Wfd%mband)

!************************************************************************

 call iter_alloc(iter_bks,(/Wfd%nkibz,Wfd%nsppol/))

 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz
     if (present(bks_mask)) then
       call wfd%distribute_bands(ik_ibz,spin,my_nband,my_band_list,bmask=bks_mask(:,ik_ibz,spin))
     else
       call wfd%distribute_bands(ik_ibz,spin,my_nband,my_band_list)
     end if
     call iter_push(iter_bks,ik_ibz,spin,my_band_list(1:my_nband))
   end do
 end do

end function wfd_iterator_bks
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_bks_distrb
!! NAME
!!  wfd_bks_distrb
!!
!! FUNCTION
!!  Build a local logical table indexed by bands, k-points and spin that defines
!!  the distribution of the load inside the loops according to the availability of the ug.
!!
!! INPUTS
!!  Wfd<wfd_t>=
!!  [bks_mask(Wfd%mband,Wfd%nkibz,Wfd%nsppol)]=Mask used to skip selecter (b,k,s) entries.
!!  [got(Wfd%nproc)]=The number of tasks already assigned to the nodes.
!!
!! OUTPUT
!!  bks_distrbk(Wfd%mband,Wfd%nkibz,Wfd%nsppol)=Global table with the rank of the node treating (b,k,s)
!!
!! PARENTS
!!      wfd_pawrhoij
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_bks_distrb(Wfd, bks_distrb, got, bks_mask)

!Arguments ------------------------------------
!scalars
 class(wfd_t),intent(in) :: Wfd
!arrays
 integer,intent(out) :: bks_distrb(Wfd%mband,Wfd%nkibz,Wfd%nsppol)
 integer,optional,intent(inout) :: got(Wfd%nproc)
 logical,optional,intent(in) :: bks_mask(Wfd%mband,Wfd%nkibz,Wfd%nsppol)

!Local variables ------------------------------
!scalars
 integer :: ik_ibz,spin,band,how_many,idle
 character(len=500) :: msg
!arrays
 integer :: get_more(Wfd%nproc),proc_ranks(Wfd%nproc)
 logical :: rank_mask(Wfd%nproc)

!************************************************************************

 get_more=0; if (present(got)) get_more=got

 ! Initialize the table here to avoid problem with the cycle instruction below.
 bks_distrb = xmpi_undefined_rank

 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz
     do band=1,Wfd%nband(ik_ibz,spin)
       if (present(bks_mask)) then
         if (.not.bks_mask(band, ik_ibz, spin)) CYCLE
       end if

       call wfd_who_has_ug(Wfd, band, ik_ibz, spin, how_many, proc_ranks)

       if (how_many == 1) then
         ! I am the only one owing this band. Add it to list.
         bks_distrb(band, ik_ibz, spin) = proc_ranks(1)

       else if (how_many>1) then
         ! This band is duplicated. Assign it trying to obtain a good load distribution.
         rank_mask=.FALSE.; rank_mask(proc_ranks(1:how_many)+1)=.TRUE.
         idle = imin_loc(get_more,mask=rank_mask)
         get_more(idle) = get_more(idle) + 1
         bks_distrb(band,ik_ibz,spin) = proc_ranks(idle)

       else
         call wfd%dump_errinfo()
         write(msg,'(a,3(i0,1x))')" Nobody has (band, ik_ibz, spin): ",band,ik_ibz,spin
         MSG_ERROR(msg)
       end if

     end do
   end do
 end do

 if (present(got)) got=get_more

end subroutine wfd_bks_distrb
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_sanity_check
!! NAME
!!  wfd_sanity_check
!!
!! FUNCTION
!!  Debugging tool
!!
!! INPUTS
!!  Wfd<wfd_t>=
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_sanity_check(Wfd)

!Arguments ------------------------------------
!scalars
 class(wfd_t),intent(inout) :: Wfd

!Local variables ------------------------------
!scalars
 integer :: ik_ibz,spin,band,mpi_ierr,ierr,how_manyb,unt_dbg,irank
 character(len=500) :: msg
!arrays
 integer :: my_band_list(Wfd%mband)

!************************************************************************

 call wfd%update_bkstab()
 ierr=0

 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz
      do band=1,Wfd%nband(ik_ibz,spin)
        if (Wfd%bks_tab(band, ik_ibz, spin, Wfd%my_rank) == WFD_STORED .and. &
           .not. wfd%ihave_ug(band, ik_ibz, spin, how="Stored") ) then
          write(msg,'(a,3(i0,1x))')" Found inconsistency in bks_tab for (band, ik_ibz, spin): ",band,ik_ibz,spin
          call wrtout(std_out, msg)
          ierr = ierr+1
        end if
     end do
   end do
 end do

 call xmpi_sum(ierr,Wfd%comm,mpi_ierr)

 if (ierr/=0) then
   if (open_file("__WFD_DEBUG__",msg,newunit=unt_dbg,form="formatted") /= 0) then
     MSG_ERROR(msg)
   end if

   do irank=0,Wfd%nproc-1
     if (irank==Wfd%my_rank) then
       write(unt_dbg,*)" (k,b,s) states owned by rank: ",Wfd%my_rank

       do spin=1,Wfd%nsppol
         do ik_ibz=1,Wfd%nkibz
            write(unt_dbg,*)" (spin,ik_ibz) ",spin,ik_ibz
            call wfd%mybands(ik_ibz, spin, how_manyb, my_band_list, how="Stored")
            write(unt_dbg,*) (my_band_list(band),band=1,how_manyb)
          end do
       end do

     end if
   end do
   close(unt_dbg)
   call xmpi_barrier(Wfd%comm)
   MSG_ERROR("Sanity check failed. Check WFD_DEBUG")
 end if

end subroutine wfd_sanity_check
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_dump_errinfo
!! NAME
!!  wfd_dump_errinfo
!!
!! FUNCTION
!!
!! INPUTS
!!  Wfd<wfd_t>=
!!
!! OUTPUT
!!
!! PARENTS
!!      m_wfd
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_dump_errinfo(Wfd,onfile)

!Arguments ------------------------------------
!scalars
 logical,optional,intent(in) :: onfile
 class(wfd_t),intent(in) :: Wfd

!Local variables ------------------------------
!scalars
 integer :: ik_ibz,spin,band,how_manyb,unt_dbg
 character(len=10) :: strank
 character(len=500) :: msg
 character(len=fnlen) :: fname_dbg
!arrays
 integer :: my_band_list(Wfd%mband)

!************************************************************************

 unt_dbg=std_out

 if (present(onfile)) then
   if (onfile) then
     call int2char10(Wfd%my_rank,strank)
     fname_dbg = "WFD_DEBUG_RANK"//TRIM(strank)
     if (open_file(fname_dbg,msg,newunit=unt_dbg,form="formatted") /= 0) then
       MSG_ERROR(msg)
     end if
   end if
 end if

 write(unt_dbg,*)" (k,b,s) states owned by rank: ",Wfd%my_rank
 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz
      write(unt_dbg,*)" ug stored at (ik_ibz, spin) ",ik_ibz,spin
      call wfd%mybands(ik_ibz, spin, how_manyb, my_band_list, how="Stored")
      write(unt_dbg,*) (my_band_list(band),band=1,how_manyb)
    end do
 end do

end subroutine wfd_dump_errinfo
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_distribute_bbp
!! NAME
!!  wfd_distribute_bbp
!!
!! FUNCTION
!!  Distribute a set of (b,b') indices taking into account the MPI distribution of the ug.
!!  It is used to calculate matrix elements of the form <b,k,s|O|b',k,s>
!!
!! INPUTS
!!  Wfd<wfd_t>=
!!  ik_ibz=The index of the k-point in the IBZ.
!!  spin=Spin index.
!!  allup=String used to select or not the upper triangle. Possible values:
!!    "All"  =Entire (b,b') matrix will be distributed.
!!    "Upper"=Only the upper triangle is distributed.
!!  [got(%nproc)]=The number of tasks already assigned to the nodes. Used to optimize the work load.
!!    Be careful when this routine is called inside several loops since each node should call the routine
!!    at each iteration with the same (local) copy of got so that bbp_distrb will assume the same value on each node.
!!  [bbp_mask(%mband,%mband)]= mask used to select a subset of (b,b') indices.
!!
!! OUTPUT
!!  my_nbbp=The number of (b,b') indices treated by this node.
!!  bbp_distrb(%mband%mband)=The rank of the node that will treat (b,b').
!!
!! PARENTS
!!      calc_optical_mels,calc_vhxc_me,cchi0q0
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_distribute_bbp(Wfd,ik_ibz,spin,allup,my_nbbp,bbp_distrb,got,bbp_mask)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,spin
 integer,intent(out) :: my_nbbp
 class(wfd_t),intent(in) :: Wfd
 character(len=*),intent(in) :: allup
!arrays
 integer,intent(out) :: bbp_distrb(Wfd%mband,Wfd%mband)
 integer,optional,intent(inout) :: got(Wfd%nproc)
 logical,optional,intent(in) :: bbp_mask(Wfd%mband,Wfd%mband)

!Local variables ------------------------------
!arrays
 integer :: loc_got(Wfd%nproc)

!************************************************************************

 ! Just a wrapper around wfd_distribute_kb_kpbp.
 loc_got=0; if (present(got)) loc_got = got

 if (present(bbp_mask)) then
   call wfd%distribute_kb_kpbp(ik_ibz,ik_ibz,spin,allup,my_nbbp,bbp_distrb,loc_got,bbp_mask)
 else
   call wfd%distribute_kb_kpbp(ik_ibz,ik_ibz,spin,allup,my_nbbp,bbp_distrb,loc_got)
 end if

end subroutine wfd_distribute_bbp
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_distribute_kb_kpbp
!! NAME
!!  wfd_distribute_kb_kpbp
!!
!! FUNCTION
!!  This routines distributes as set of (b,b') indices taking into account the MPI distribution of the ug.
!!  It is used to calculate matrix elements of the form <b,k,s|O|b',k',s>
!!
!! INPUTS
!!  Wfd<wfd_t>=
!!  ik_ibz =The index of the k-point k  in the IBZ.
!!  ikp_ibz=The index of the k-point k' in the IBZ.
!!  spin=Spin index.
!!  allup=String used to select the upper triangle of the (b,b') matrix. Possible values:
!!    "All"  =Entire (b,b') matrix will be distributed.
!!    "Upper"=Only the upper triangle is distributed.
!!  [got(%nproc)]=The number of tasks already assigned to the nodes. Used to optimize the distribution of the tasks.
!!    Be careful when this routine is called inside several loops since each node should call the routine
!!    at each iteration with the same (local) copy of got so that bbp_distrb will assume the same value on each node.
!!  [bbp_mask(%mband,%mband)]= mask used to select a subset of (b,b') indices.
!!
!! OUTPUT
!!  my_nbbp=The number of (b,b') indices treated by this node.
!!  bbp_distrb(%mband%mband)=The rank of the node that will treat (b,b').
!!
!! PARENTS
!!      cchi0,m_wfd
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_distribute_kb_kpbp(Wfd,ik_ibz,ikp_ibz,spin,allup,my_nbbp,bbp_distrb,got,bbp_mask)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,ikp_ibz,spin
 integer,intent(out) :: my_nbbp
 class(wfd_t),intent(in) :: Wfd
 character(len=*),intent(in) :: allup
!arrays
 integer,intent(out) :: bbp_distrb(Wfd%mband,Wfd%mband)
 integer,optional,intent(inout) :: got(Wfd%nproc)
 logical,optional,intent(in) :: bbp_mask(Wfd%mband,Wfd%mband)

!Local variables ------------------------------
!scalars
 integer :: my_nband,ib1,ib2,pcb2,pcb1,howmany_b,howmany_bp,workload_min
 integer :: rank,ncpus,idle,b1_stop,ierr
 character(len=500) :: msg
!arrays
 integer :: rank_bandlist_k(Wfd%mband),rank_bandlist_kp(Wfd%mband)
 integer :: get_more(Wfd%nproc),my_band_list_k(Wfd%mband)
 integer,allocatable :: whocan_k(:,:),whocan_kp(:,:)
 logical :: b_mask(Wfd%mband)

!************************************************************************

 ABI_MALLOC_OR_DIE(whocan_k ,(Wfd%mband,Wfd%nproc), ierr)
 ABI_MALLOC_OR_DIE(whocan_kp,(Wfd%mband,Wfd%nproc), ierr)
 whocan_k =0 !  Will be set to 1 if this node can calculate something containing (k,b)
 whocan_kp=0 !  Will be set to 1 if this node can calculate something containing (kp,bp)

 do rank=0,Wfd%nproc-1

   call wfd_bands_of_rank(Wfd,rank,ik_ibz ,spin,howmany_b, rank_bandlist_k )
   do pcb1=1,howmany_b
     ib1 = rank_bandlist_k(pcb1)
     whocan_k(ib1,rank+1) = 1
   end do

   call wfd_bands_of_rank(Wfd,rank,ikp_ibz,spin,howmany_bp,rank_bandlist_kp)
   do pcb2=1,howmany_bp
     ib2 = rank_bandlist_kp(pcb2)
     whocan_kp(ib2,rank+1) = 1
   end do

 end do

 get_more=0; if (present(got)) get_more=got
 b1_stop=Wfd%nband(ik_ibz,spin)

 bbp_distrb = xmpi_undefined_rank

 do ib2=1,Wfd%nband(ikp_ibz,spin)
   b_mask = .TRUE.; if (present(bbp_mask)) b_mask = bbp_mask(:,ib2)
   if (ANY(b_mask)) then
     my_nband=0; my_band_list_k=0
     ! Only the upper triangle of the (b1,b2) matrix.
     if (firstchar(allup, ["U","u"])) b1_stop = MIN(ib2,Wfd%nband(ik_ibz,spin))

     do ib1=1,b1_stop
       if (b_mask(ib1)) then
         !
         ! find which CPUs can do the calculation (k,b)->(kp,bp)
         ! find the one which is less busy
         ncpus=0
         workload_min=HUGE(0)
         do rank=0,Wfd%nproc-1
           if( whocan_k(ib1,rank+1)==1 .AND.  whocan_kp(ib2,rank+1)==1 ) then
             ncpus=ncpus+1
             if( get_more(rank+1) < workload_min ) then
               idle=rank+1
               workload_min=get_more(idle)
             end if

           end if
         end do

         if(ncpus>0) then
           bbp_distrb(ib1,ib2)=idle-1
           get_more(idle) = get_more(idle) + 1

         else
           call wfd%dump_errinfo()
           write(msg,'(a,5(i0,1x))')" Nobody has (band1, ik_ibz) (band2, ikp_ibz) spin: ",ib1,ik_ibz,ib2,ikp_ibz,spin
           MSG_ERROR(msg)
         end if

       end if
     end do ! ib1
   end if
 end do ! ib2

 ABI_FREE(whocan_k)
 ABI_FREE(whocan_kp)

 my_nbbp = COUNT(bbp_distrb==Wfd%my_rank)
 if (present(got)) got=get_more

end subroutine wfd_distribute_kb_kpbp
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_get_cprj
!! NAME
!!  wfd_get_cprj
!!
!! FUNCTION
!!  Return a copy of Cprj either by calculating it on-the-fly or by just retrieving the data already stored in the data type.
!!
!! INPUTS
!!  Wfd<wfd_t>=the wavefunction descriptor.
!!  band=Band index.
!!  ik_ibz=Index of the k-point in the IBZ.
!!  spin=Spin index
!!  sorted=.TRUE. if the output cprj matrix elements have to be sorted by atom type.
!!
!! OUTPUT
!!  Cprj_out(Wfd%natom,Wfd%nspinor) <type(pawcprj_type)>=Unsorted matrix elements.
!!
!! PARENTS
!!      calc_optical_mels,calc_sigc_me,calc_sigx_me,calc_vhxc_me,cchi0,cchi0q0
!!      cchi0q0_intraband,cohsex_me,debug_tools,exc_build_block,exc_build_ham
!!      m_wfd,prep_calc_ucrpa,sigma,wfd_pawrhoij,wfd_vnlpsi
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_get_cprj(Wfd, band, ik_ibz, spin, Cryst, Cprj_out, sorted)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin
 logical,intent(in) :: sorted
 class(wfd_t),intent(inout) :: Wfd
 type(crystal_t),intent(in) :: Cryst
!arrays
 type(pawcprj_type),intent(inout) :: Cprj_out(Wfd%natom,Wfd%nspinor)

!Local variables ------------------------------
!scalars
 integer,parameter :: choice1=1,idir0=0
 integer :: want_order,iatom,sidx
 character(len=500) :: msg
 type(wave_t),pointer :: wave

!************************************************************************

 want_order=CPR_RANDOM; if (sorted) want_order=CPR_SORTED

 ABI_CHECK(wfd%get_wave_ptr(band, ik_ibz, spin, wave, msg) == 0, msg)

 select case (wave%has_cprj)

 case (WFD_NOWAVE, WFD_ALLOCATED)
   ! Have to calculate it!
   if (.not. wave%has_ug == WFD_STORED) then
     write(msg,'(a,3(i0,1x),a)')" ug for (band, ik_ibz, spin): ",band,ik_ibz,spin," is not stored in memory!"
     MSG_ERROR(msg)
   end if
   ! Get cprj.
   call wfd%ug2cprj(band,ik_ibz,spin,choice1,idir0,Wfd%natom,Cryst,Cprj_out,sorted=sorted)

   if (wave%has_cprj == WFD_ALLOCATED) then
     ! Store it.
     if (want_order == wave%cprj_order) then
       call pawcprj_copy(Cprj_out, wave%Cprj)
       wave%has_cprj = WFD_STORED

     else
       ! Have to reorder cprj_out
       select case (want_order)
       case (CPR_SORTED)
         do iatom=1,Cryst%natom
           sidx = Cryst%atindx(iatom) ! random --> sorted table.
           call pawcprj_copy(Cprj_out(sidx:sidx,:), wave%Cprj(iatom:iatom,:))
         end do
       case (CPR_RANDOM)
         do sidx=1,Cryst%natom
           iatom = Cryst%atindx1(sidx) ! sorted --> random table.
           call pawcprj_copy(Cprj_out(iatom:iatom,:), wave%Cprj(sidx:sidx,:))
         end do
       case default
         MSG_ERROR(sjoin("Wrong value for want_order:", itoa(want_order)))
       end select
     end if
   end if

 case (WFD_STORED)
   ! copy it back.
   if (want_order == wave%cprj_order) then
     call pawcprj_copy(wave%Cprj,Cprj_out)

   else
     select case (want_order)
     case (CPR_SORTED)
       do iatom=1,Cryst%natom
         sidx = Cryst%atindx(iatom) ! random --> sorted table.
         call pawcprj_copy(wave%Cprj(iatom:iatom,:),Cprj_out(sidx:sidx,:))
       end do
     case (CPR_RANDOM)
       do sidx=1,Cryst%natom
         iatom = Cryst%atindx1(sidx) ! sorted --> random table.
         call pawcprj_copy(wave%Cprj(sidx:sidx,:),Cprj_out(iatom:iatom,:))
       end do
     case default
       MSG_ERROR(sjoin("Wrong value for want_order:", itoa(want_order)))
     end select
   end if

 case default
   MSG_BUG(sjoin("Wrong has_cprj: ", itoa(wave%has_cprj)))
 end select

end subroutine wfd_get_cprj
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_change_ngfft
!! NAME
!!  wfd_change_ngfft
!!
!! FUNCTION
!!   Reallocate and reinitialize internal tables for performing FFTs of wavefunctions.
!!
!! INPUTS
!!  Cryst<crystal_t>=Info on unit cell.
!!  Psps<pseudopotential_type>=Pseudopotential info.
!!  new_ngfft(18)=FFT descriptor for the new FFT mesh.
!!
!!  SIDE EFFECTS
!!  Wfd<wfd_t>=Wavefunction descriptor with new internal tables for FFT defined by new_ngfft.
!!
!! PARENTS
!!      calc_sigc_me,calc_sigx_me,calc_vhxc_me,cchi0,cchi0q0,cchi0q0_intraband
!!      classify_bands,cohsex_me,exc_build_block,exc_build_ham,exc_plot
!!      m_bseinterp,m_wfd,prep_calc_ucrpa,screening,sigma,wfd_mkrho
!!      wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_change_ngfft(Wfd,Cryst,Psps,new_ngfft)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: new_ngfft(18)
 type(crystal_t),intent(in) :: Cryst
 type(pseudopotential_type),intent(in) :: Psps
 class(wfd_t),intent(inout) :: Wfd

!Local variables ------------------------------
!scalars
 integer,parameter :: npw0=0
 integer :: npw_k, ik_ibz, istwf_k, is, ik, ib
 logical :: iscompatibleFFT
 character(len=500) :: msg
!arrays
 integer,allocatable :: kg_k(:,:)

!************************************************************************

 DBG_ENTER("COLL")

 !@wfd_t
 if ( ALL(Wfd%ngfft(1:3) == new_ngfft(1:3)) ) RETURN ! Nothing to do.

 if (Wfd%prtvol > 0) then
   write(msg,"(a,3(i0,1x),a,3(i0,1x),a)")"Changing FFT mesh: [",Wfd%ngfft(1:3),"] ==> [",new_ngfft(1:3),"]"
   call wrtout(std_out, msg)
 end if

 ! Change FFT dimensions.
 Wfd%ngfft  = new_ngfft
 Wfd%mgfft  = MAXVAL(new_ngfft(1:3))
 Wfd%nfftot = PRODUCT(new_ngfft(1:3))
 Wfd%nfft   = Wfd%nfftot ! No FFT parallelism.

 ! Re-initialize fft distribution
 call destroy_distribfft(Wfd%MPI_enreg%distribfft)
 call init_distribfft(Wfd%MPI_enreg%distribfft,'c',Wfd%MPI_enreg%nproc_fft,new_ngfft(2),new_ngfft(3))

 ABI_REMALLOC(Wfd%ph1d,(2,3*(2*Wfd%mgfft+1)*Cryst%natom))
 call getph(Cryst%atindx,Cryst%natom,Wfd%ngfft(1),Wfd%ngfft(2),Wfd%ngfft(3),Wfd%ph1d,Cryst%xred)

 ! Recalculate FFT tables.
 ! Calculate the FFT index of $ R^{-1} (r-\tau) $ used to symmetrize u_Rk.
 ABI_REMALLOC(Wfd%irottb, (Wfd%nfftot,Cryst%nsym))
 call rotate_FFT_mesh(Cryst%nsym,Cryst%symrel,Cryst%tnons,Wfd%ngfft,Wfd%irottb,iscompatibleFFT)

 if (.not.iscompatibleFFT) then
   MSG_WARNING("FFT mesh not compatible with symmetries. Wavefunction symmetrization should not be done in r-space!")
 end if

 ! Is the new real space FFT mesh compatible with the rotational part?
 Wfd%rfft_is_symok = check_rot_fft(Cryst%nsym,Cryst%symrel,Wfd%ngfft(1),Wfd%ngfft(2),Wfd%ngfft(3))

 ! Reallocate ur buffers with correct dimensions.
 do is=1,size(wfd%s)
   do ik=1,size(wfd%s(is)%k)
     do ib=1,size(wfd%s(is)%k(ik)%b)
       call wfd%s(is)%k(ik)%b(ib)%free(what="R")
     end do
   end do
 end do

 ! Reinit Kdata_t
 do ik_ibz=1,Wfd%nkibz
   if (any(wfd%bks2wfd(1, :, ik_ibz, :) /= 0)) then
     istwf_k = Wfd%istwfk(ik_ibz)
     npw_k   = Wfd%Kdata(ik_ibz)%npw
     ABI_MALLOC(kg_k, (3,npw_k))
     kg_k = Wfd%Kdata(ik_ibz)%kg_k
     call kdata_free(Wfd%Kdata(ik_ibz))
     call kdata_init(Wfd%Kdata(ik_ibz),Cryst,Psps,Wfd%kibz(:,ik_ibz),istwf_k,new_ngfft,Wfd%MPI_enreg,kg_k=kg_k)
     ABI_FREE(kg_k)
   end if
 end do

 DBG_EXIT("COLL")

end subroutine wfd_change_ngfft
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_test_ortho
!! NAME
!! wfd_test_ortho
!!
!! FUNCTION
!!  Test the orthonormalization of the wavefunctions stored in Wfd.
!!
!! INPUTS
!!  Wfd<wfd_t>=wavefunction descriptor.
!!  Cryst<crystal_t>=Object defining the unit cell and its symmetries.
!!  Pawtab(ntypat*usepaw)<type(pawtab_type)>=PAW tabulated starting data.
!!
!! OUTPUT
!!   Only writing.
!!
!! PARENTS
!!      bethe_salpeter,m_gkk,m_phgamma,m_phpi,m_sigmaph,screening
!!      sigma,wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_test_ortho(Wfd,Cryst,Pawtab,unit,mode_paral)

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: unit
 character(len=4),optional,intent(in) :: mode_paral
 type(crystal_t),intent(in) :: Cryst
 class(wfd_t),target,intent(inout) :: Wfd
!array
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)

!Local variables ------------------------------
!scalars
 integer :: ik_ibz,spin,band,band1,band2,ib,ib1,ib2,ierr,how_manyb,my_unt,npw_k,istwf_k
 real(dp) :: glob_cinf,my_cinf,glob_csup,my_csup,glob_einf,min_norm2,glob_esup,max_norm2
 complex(dpc) :: cdum
 logical :: bands_are_spread
 character(len=4) :: my_mode
 character(len=500) :: msg
 type(wave_t),pointer :: wave1, wave2
!arrays
 integer :: my_bandlist(Wfd%mband)
 real(dp) :: pawovlp(2)
 complex(gwpc),ABI_CONTIGUOUS pointer :: ug1(:),ug2(:)
 !complex(gwpc) :: ur(Wfd%nfft*Wfd%nspinor)
 character(len=6) :: tag_spin(2)
 type(pawcprj_type),allocatable :: Cp1(:,:),Cp2(:,:)

!************************************************************************

 tag_spin(:)=(/'      ','      '/); if (Wfd%nsppol==2) tag_spin(:)=(/' UP   ',' DOWN '/)

 my_unt   =std_out; if (present(unit      )) my_unt   =unit
 my_mode  ='COLL' ; if (present(mode_paral)) my_mode  =mode_paral

 if (Wfd%usepaw==1) then
   ABI_MALLOC(Cp1,(Wfd%natom,Wfd%nspinor))
   call pawcprj_alloc(Cp1,0,Wfd%nlmn_atm)
   ABI_MALLOC(Cp2,(Wfd%natom,Wfd%nspinor))
   call pawcprj_alloc(Cp2,0,Wfd%nlmn_atm)
 end if

 bands_are_spread = .FALSE.

 do spin=1,Wfd%nsppol
   min_norm2=greatest_real; max_norm2=-greatest_real
   my_cinf=greatest_real;  my_csup=-greatest_real
   do ik_ibz=1,Wfd%nkibz
     istwf_k = Wfd%istwfk(ik_ibz)
     npw_k   = Wfd%npwarr(ik_ibz)

     ! Select my band indices.
     call wfd%mybands(ik_ibz,spin,how_manyb,my_bandlist, how="Stored")
     if (how_manyb/=Wfd%nband(ik_ibz,spin)) bands_are_spread = .TRUE.

     ! 1) Normalization.
     do ib=1,how_manyb
       band = my_bandlist(ib)
       ABI_CHECK(wfd%get_wave_ptr(band, ik_ibz, spin, wave1, msg) == 0, msg)
       ug1 => wave1%ug
       cdum = xdotc(npw_k*Wfd%nspinor,ug1,1,ug1,1)
       if (istwf_k>1) then
         cdum=two*DBLE(cdum)
         if (istwf_k==2) cdum=cdum-CONJG(ug1(1))*ug1(1)
       end if
       if (Wfd%usepaw==1) then
         call wfd%get_cprj(band,ik_ibz,spin,Cryst,Cp1,sorted=.FALSE.)
         pawovlp = paw_overlap(Cp1,Cp1,Cryst%typat,Pawtab,spinor_comm=Wfd%MPI_enreg%comm_spinor)
         cdum = cdum + CMPLX(pawovlp(1),pawovlp(2), kind=dpc)
       end if
       !write(std_out,*)"ik_ibz, band, spin, cdum: ",ik_ibz,band,spin,cdum
       if (REAL(cdum)<min_norm2) min_norm2=REAL(cdum)
       if (REAL(cdum)>max_norm2) max_norm2=REAL(cdum)
     end do

     ! TODO should use the communicator for this spin
     call xmpi_min(min_norm2,glob_einf,Wfd%comm,ierr)
     call xmpi_max(max_norm2,glob_esup,Wfd%comm,ierr)

     ! 2) Orthogonality of wavefunctions.
     do ib1=1,how_manyb
       band1 = my_bandlist(ib1)
       ABI_CHECK(wfd%get_wave_ptr(band1, ik_ibz, spin, wave1, msg) == 0, msg)
       ug1 => wave1%ug
       if (Wfd%usepaw==1) call wfd%get_cprj(band1,ik_ibz,spin,Cryst,Cp1,sorted=.FALSE.)

       do ib2=ib1+1,how_manyb
         band2 = my_bandlist(ib2)
         ABI_CHECK(wfd%get_wave_ptr(band2, ik_ibz, spin, wave2, msg) == 0, msg)
         ug2 => wave2%ug
         if (Wfd%usepaw==1) call wfd%get_cprj(band2,ik_ibz,spin,Cryst,Cp2,sorted=.FALSE.)

         cdum = xdotc(npw_k*Wfd%nspinor,ug1,1,ug2,1)
         if (istwf_k>1) then
           cdum=two*DBLE(cdum)
           if (istwf_k==2) cdum=cdum-CONJG(ug1(1))*ug2(1)
         end if
         if (Wfd%usepaw==1) then
           pawovlp = paw_overlap(Cp1,Cp2,Cryst%typat,Pawtab,spinor_comm=Wfd%MPI_enreg%comm_spinor)
           cdum = cdum + CMPLX(pawovlp(1),pawovlp(2), kind=dpc)
         end if

         if (ABS(cdum)<my_cinf) my_cinf=ABS(cdum)
         if (ABS(cdum)>my_csup) my_csup=ABS(cdum)
         !if (ABS(cdum) > 0.1) write(std_out,*)" ib1,ib2,ABS_dotprod: ",ib1,ib2,ABS(cdum)
       end do !ib2
     end do !ib

     ! TODO should use the communicator for this spin
     call xmpi_min(my_cinf,glob_cinf,Wfd%comm,ierr)
     call xmpi_max(my_csup,glob_csup,Wfd%comm,ierr)
   end do ! ik_ibz

   ! Output results for this spin
   write(msg,'(2a)')ch10,' test on the normalization of the wavefunctions'
   if (Wfd%nsppol==2) write(msg,'(3a)')ch10,' test on the normalization of the wavefunctions with spin ',tag_spin(spin)
   call wrtout(my_unt,msg,mode_paral)
   write(msg,'(a,f9.6,a,a,f9.6)')&
     ' min sum_G |a(n,k,G)| = ',glob_einf,ch10,&
     ' max sum_G |a(n,k,G)| = ',glob_esup
   call wrtout(my_unt,msg,mode_paral)

   write(msg,'(a)')' test on the orthogonalization of the wavefunctions (NB: this is not invariant for degenerate states)'
   if (Wfd%nsppol==2) write(msg,'(2a)')' test on the orthogonalization of the wavefunctions with spin ',tag_spin(spin)
   call wrtout(my_unt,msg,mode_paral)
   write(msg,'(a,f9.6,a,a,f9.6,a)')&
     '- min sum_G a(n,k,G)a(n",k,G) = ',glob_cinf,ch10,&
     '- max sum_G a(n,k,G)a(n",k,G) = ',glob_csup,ch10
   call wrtout(my_unt,msg,mode_paral)

 end do ! spin

 if (bands_are_spread) then
   write(msg,'(3a)')&
     'Note that the test on the orthogonalization is not complete ',ch10,&
     'since bands are spread among different processors'
   call wrtout(my_unt,msg,mode_paral)
 end if

 if (Wfd%usepaw==1) then
   call pawcprj_free(Cp1)
   ABI_FREE(Cp1)
   call pawcprj_free(Cp2)
   ABI_FREE(Cp2)
 end if

end subroutine wfd_test_ortho
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_sym_ur
!! NAME
!!  wfd_sym_ur
!!
!! FUNCTION
!!  Symmetrize a wave function in real space
!!
!! INPUTS
!!  Wfd<wfd_t>=the wavefunction descriptor.
!!  Cryst<crystal_t>=Structure describing the crystal structure and its symmetries.
!!  Kmesh<kmesh_t>=Structure describing the BZ sampling
!!  band=Band index.
!!  ik_bz=Index of the k-point in the BZ.
!!  spin=Spin index
!!  [trans] = "N" if only the symmetried wavefunction is needed, "C" if the complex conjugate is required.
!!            Default is "N"
!!  [with_umklp] = Optional flag. If .True. (Default) the umklapp G0 vector in the relation kbz = Sk + G0
!!                 is taken into account when constructing u_kbz.
!!
!! NOTES
!!  This method is deprecated. See wfd_sym_ug_kg for symmetrization in G-space
!!
!! OUTPUT
!!  ur_kbz(Wfd%nfft*Wfd%nspinor)=The symmetrized wavefunction in real space.
!!  [ur_kibz(Wfd%nfft*Wfd%nspinor)]= Optional output: u(r) in the IBZ.
!!
!! PARENTS
!!      debug_tools,exc_plot,m_bseinterp
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_sym_ur(Wfd,Cryst,Kmesh,band,ik_bz,spin,ur_kbz,trans,with_umklp,ur_kibz)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_bz,spin
 character(len=*),optional,intent(in) :: trans
 logical,optional,intent(in) :: with_umklp
 type(crystal_t),intent(in) :: Cryst
 type(kmesh_t),intent(in) :: Kmesh
 class(wfd_t),intent(inout) :: Wfd
!arrays
 complex(gwpc),intent(out) :: ur_kbz(Wfd%nfft*Wfd%nspinor)
 complex(gwpc),optional,intent(out) :: ur_kibz(Wfd%nfft*Wfd%nspinor)

!Local variables ------------------------------
!scalars
 integer :: ik_ibz,isym_k,itim_k,nr,ispinor,spad,ir,ir2
 integer :: fft_idx,ix,iy,iz,nx,ny,nz,irot
 real(dp) :: gdotr
 complex(dpc) :: ph_mkt,u2b,u2a
 complex(gwpc) :: gwpc_ph_mkt
 logical :: isirred,my_with_umklp
 character(len=1) :: my_trans
 !character(len=500) :: msg
!arrays
 integer :: umklp(3)
 real(dp) :: kbz(3),spinrot_k(4)
 complex(dpc) :: spinrot_mat(2,2)
 complex(gwpc),allocatable :: ur(:)

!************************************************************************

 my_trans = "N"; if (present(trans)) my_trans = toupper(trans(1:1))
 my_with_umklp = .TRUE.; if (present(with_umklp)) my_with_umklp = with_umklp

 ! k_bz =  S k - G0 ==> u_{k_bz} =  e^{iG0.r} u_{Sk}
 ! k_bz = -S k - G0 ==> u_{k_bz} =  e^{iG0.r} u_{Sk}^*

 ! u(r,b,kbz)=e^{-2i\pi kibz.(R^{-1}t} u (R{^-1}(r-t),b,kibz)
 !           =e^{+2i\pi kibz.(R^{-1}t} u*({R^-1}(r-t),b,kibz) for time-reversal
 !
 ! Get ik_ibz, non-symmorphic phase, ph_mkt, and symmetries from ik_bz.
 call get_BZ_item(Kmesh,ik_bz,kbz,ik_ibz,isym_k,itim_k,ph_mkt,umklp,isirred)
 gwpc_ph_mkt = ph_mkt

 if (isirred) then
   ! Avoid symmetrization if this point is irreducible.
   call wfd%get_ur(band,ik_ibz,spin,ur_kbz)
   if (present(ur_kibz)) call xcopy(Wfd%nfft*Wfd%nspinor,ur_kbz,1,ur_kibz,1)
   if (my_trans=="C") ur_kbz = GWPC_CONJG(ur_kbz)
   RETURN
 end if

 ! Reconstruct ur in the BZ from the corresponding wavefunction in IBZ.
 ABI_MALLOC(ur, (Wfd%nfft*Wfd%nspinor))

 call wfd%get_ur(band,ik_ibz,spin,ur)
 if (present(ur_kibz)) call xcopy(Wfd%nfft*Wfd%nspinor,ur,1,ur_kibz,1)

 ! Wfd%irottb(:,isym_k) is the table for rotated FFT points
 SELECT CASE (Wfd%nspinor)

 CASE (1)
   ! Rotation in real space
   do ir=1,Wfd%nfft
     irot = Wfd%irottb(ir,isym_k)
     ur_kbz(ir) = ur(irot) * gwpc_ph_mkt
   end do

   ! Apply time-reversal symmetry if needed.
   if (itim_k==2) ur_kbz = GWPC_CONJG(ur_kbz)

   ! Take into account a possible umklapp.
   if (ANY(umklp/=0).and. my_with_umklp) then
     ! Compute ur_kbz = ur_kbz*eig0r
     nx = Wfd%ngfft(1); ny = Wfd%ngfft(2); nz = Wfd%ngfft(3)
     fft_idx=0
     do iz=0,nz-1
       do iy=0,ny-1
         do ix=0,nx-1
           gdotr= two_pi*( umklp(1)*(ix/DBLE(nx)) &
                          +umklp(2)*(iy/DBLE(ny)) &
                          +umklp(3)*(iz/DBLE(nz)) )
           fft_idx = fft_idx+1
           ur_kbz(fft_idx) = ur_kbz(fft_idx) * DCMPLX(DCOS(gdotr),DSIN(gdotr))
         end do
       end do
     end do
   end if

   if (my_trans=="C") ur_kbz = GWPC_CONJG(ur_kbz)

 CASE (2)
   MSG_ERROR("Implementation has to be tested")

   nr = Wfd%nfft
   spinrot_k = Cryst%spinrot(:,isym_k)
   !
   ! ==== Apply Time-reversal if required ====
   ! \psi_{-k}^1 =  (\psi_k^2)^*
   ! \psi_{-k}^2 = -(\psi_k^1)^*
   if (itim_k==1) then
     ur_kbz = ur
   else if (itim_k==2) then
     ur_kbz(1:nr)     = GWPC_CONJG(ur(nr+1:2*nr))
     ur_kbz(nr+1:2*nr)=-GWPC_CONJG(ur(1:nr))
   else
     MSG_ERROR('Wrong i2 in spinor')
   end if
   !
   ! Rotate wavefunctions in real space.
   do ispinor=1,Wfd%nspinor
     spad=(ispinor-1)*nr
     do ir=1,nr
       ir2 = Wfd%irottb(ir,isym_k)
       ur_kbz(ir+spad) = ur_kbz(ir2+spad) * gwpc_ph_mkt
     end do
   end do
   !
   ! Rotation in spinor space.
   spinrot_mat(1,1)= spinrot_k(1) + j_dpc*spinrot_k(4)
   spinrot_mat(1,2)= spinrot_k(3) + j_dpc*spinrot_k(2)
   spinrot_mat(2,1)=-spinrot_k(3) + j_dpc*spinrot_k(2)
   spinrot_mat(2,2)= spinrot_k(1) - j_dpc*spinrot_k(4)

   do ir=1,nr
     u2a=ur_kbz(ir)
     u2b=ur_kbz(ir+nr)
     ur_kbz(ir)   =spinrot_mat(1,1)*u2a+spinrot_mat(1,2)*u2b
     ur_kbz(ir+nr)=spinrot_mat(2,1)*u2a+spinrot_mat(2,2)*u2b
   end do

   if (ANY(umklp /=0)) then
     !ur_kbz(1:Wfd%nfft)  = ur_kbz(1:Wfd%nfft) *eig0r
     !ur_kbz(Wfd%nfft+1:) = ur_kbz(Wfd%nfft+1:)*eig0r
   end if

 CASE DEFAULT
   MSG_ERROR(sjoin("Wrong value for nspinor: ", itoa(Wfd%nspinor)))
 END SELECT

 ABI_FREE(ur)

end subroutine wfd_sym_ur
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_sym_ug_kg
!! NAME
!!  wfd_sym_ug_kg
!!
!! FUNCTION
!!  Use crystalline symmetries and time reversal to reconstruct wavefunctions at kk_bz from the IBZ image kk_ibz.
!!  Return periodic part in G-space as well as list of G-vectors belonging to the G-sphere centered on kk_bz
!!
!! INPUTS
!!  ecut: Cutoff energy for planewave basis set.
!!  kk_bz: k-point in the BZ for output wavefunctions and G-vectors.
!!  kk_ibz: Symmetrical image of kk_bz in the IBZ.
!!  bstart: Initial band
!!  nband: Number of bands to symmetrize.
!!  spin: Spin index
!!  mpw: Maximum number of planewaves used to dimension arrays.
!!  indkk: Symmetry map kk_bz -> kk_ibz as computed by listkk.
!!  cryst: Crystalline structure and symmetries
!!  work_ngfft: Define the size of the workspace array work
!!  work: Workspace array used to symmetrize wavefunctions
!!
!! OUTPUT
!!  istwf_kbz: Time-reversal flag associated to output wavefunctions
!!  npw_kbz: Number of G-vectors in kk_bz G-sphere
!!  kg_kbz: G-vectors in reduced coordinates.
!!  cgs_kbz: Periodic part of wavefunctions at kk_bz
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_sym_ug_kg(self, ecut, kk_bz, kk_ibz, bstart, nband, spin, mpw, indkk, cryst, &
                         work_ngfft, work, istwf_kbz, npw_kbz, kg_kbz, cgs_kbz)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bstart, nband, spin, mpw
 type(crystal_t),intent(in) :: cryst
 class(wfd_t),intent(in) :: self
 integer,intent(out) :: istwf_kbz, npw_kbz
 real(dp),intent(in) :: ecut
!arrays
 integer :: work_ngfft(18)
 integer,intent(in) :: indkk(6)
 integer,intent(out) :: kg_kbz(3, mpw)
 real(dp),intent(in) :: kk_bz(3), kk_ibz(3)
 real(dp),intent(out) :: cgs_kbz(2, mpw*self%nspinor, nband)
 real(dp),intent(out) :: work(2, work_ngfft(4), work_ngfft(5), work_ngfft(6))

!Local variables ------------------------------
!scalars
 integer,parameter :: ndat1 = 1
 integer :: ik_ibz, isym_k, trev_k, ib, band, istwf_kirr, npw_kirr
 logical :: isirr_k
!arrays
 integer :: g0_k(3)
 integer,allocatable :: gtmp(:,:)
 real(dp),allocatable :: cg_kirr(:,:)

!************************************************************************

 ! As reported by listkk via symrel
 ik_ibz = indkk(1); isym_k = indkk(2); trev_k = indkk(6); g0_k = indkk(3:5)
 isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))

 ! Get npw_kbz, kg_kbz and symmetrize wavefunctions from IBZ (if needed).
 ! Be careful with time-reversal symmetry.
 if (isirr_k) then
   ! Copy u_k(G)
   istwf_kbz = self%istwfk(ik_ibz); npw_kbz = self%npwarr(ik_ibz)
   ABI_CHECK(mpw >= npw_kbz, "mpw < npw_kbz")
   kg_kbz(:,1:npw_kbz) = self%kdata(ik_ibz)%kg_k
   do ib=1,nband
     band = ib + bstart - 1
     call self%copy_cg(band, ik_ibz, spin, cgs_kbz(1,1,ib))
   end do
 else
   ! Reconstruct u_k(G) from the IBZ image.
   istwf_kbz = 1
   call get_kg(kk_bz, istwf_kbz, ecut, cryst%gmet, npw_kbz, gtmp)
   ABI_CHECK(mpw >= npw_kbz, "mpw < npw_kbz")
   kg_kbz(:,1:npw_kbz) = gtmp(:,:npw_kbz)
   ABI_FREE(gtmp)

   ! Use cg_kirr as workspace array, results stored in cgs_kbz.
   istwf_kirr = self%istwfk(ik_ibz); npw_kirr = self%npwarr(ik_ibz)
   ABI_MALLOC(cg_kirr, (2, npw_kirr*self%nspinor))
   do ib=1,nband
     band = ib + bstart - 1
     call self%copy_cg(band, ik_ibz, spin, cg_kirr)
     call cgtk_rotate(cryst, kk_ibz, isym_k, trev_k, g0_k, self%nspinor, ndat1, &
                      npw_kirr, self%kdata(ik_ibz)%kg_k, &
                      npw_kbz, kg_kbz, istwf_kirr, istwf_kbz, cg_kirr, cgs_kbz(:,:,ib), work_ngfft, work)
   end do
   ABI_FREE(cg_kirr)
 end if

end subroutine wfd_sym_ug_kg
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_write_wfk
!! NAME
!! wfd_write_wfk
!!
!! FUNCTION
!!  This routine writes the wavefunctions to the specified WFK file
!!  All the wavefunction are stored on each node, only the spin is distributed.
!!
!! INPUTS
!!  Wfd<wfd_t>=Initialized wavefunction descritptor.
!!  wfk_fname=Name of the WFK file.
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_write_wfk(Wfd,Hdr,Bands,wfk_fname)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: wfk_fname
 class(wfd_t),intent(in) :: Wfd
 type(Hdr_type),intent(in) :: Hdr
 type(ebands_t),intent(in) :: Bands

!Local variables ------------------------------
!scalars
 integer,parameter :: formeig0=0,master=0
 integer :: nprocs,my_rank,iomode,cgsize,npw_k,ik_ibz,spin,nband_k,band,ii
 integer :: blk,nblocks,how_many,ierr,how_manyb
 real(dp) :: cpu,wall,gflops
 logical :: iam_master
 character(len=500) :: msg
 type(wfk_t) :: Wfkfile
!arrays
 integer :: band_block(2),proc_ranks(Wfd%nproc),my_band_list(Wfd%mband)
 integer,allocatable :: blocks(:,:) !band_list(:),
 real(dp),allocatable :: cg_k(:,:)

!************************************************************************

 DBG_ENTER("COLL")

 nprocs = xmpi_comm_size(Wfd%comm); my_rank = xmpi_comm_rank(Wfd%comm)
 iam_master = (my_rank == master)

 ! Select the IO library from the file extension.
 iomode = iomode_from_fname(wfk_fname)
 call wrtout(std_out, sjoin('Writing GS WFK file: ',wfk_fname,", with iomode ",iomode2str(iomode)))

 if (nprocs > 1 .and. iomode /= IO_MODE_MPI) then
   MSG_ERROR("You need MPI-IO to write wavefunctions in parallel")
 end if
 !
 ! Check consistency between Wfd and Header!
 ! The ideal approach would be to generate the header from the Wfd but a lot of info are missing
 ABI_CHECK(Wfd%nkibz == Hdr%nkpt,"Different number of k-points")
 ABI_CHECK(Wfd%nsppol == Hdr%nsppol,"Different number of spins")
 ABI_CHECK(Wfd%nspinor == Hdr%nspinor,"Different number of spinors")

 if (any(Wfd%nband /= reshape(Hdr%nband, [Wfd%nkibz, Wfd%nsppol]))) then
   MSG_ERROR("Wfd%nband /= Hdr%nband")
 end if

 ! Use bks_tab to decide who will write the data. Remember
 ! integer,allocatable :: bks_tab(:,:,:,:)
 ! Wfd%bks_tab(mband,nkibz,nsppol,0:nproc-1)
 ! Global table used to keep trace of the distribution of the (b,k,s) states on each node inside Wfd%comm.
 ! 1 if the node has this state. 0 otherwise.
 ! A node owns a wavefunction if the corresponding ug is allocated AND computed.
 ! If a node owns ur but not ug, or ug is just allocated then its entry in the table is zero.
 ! The main difficulties here are:
 !
 ! 1) FFT parallelism (not coded, indeed)
 ! 2) Wavefunctions that are replicated, i.e. the same (b,k,s) is treated by more than one node.

 ierr = 0
 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz
     do band=1,Wfd%nband(ik_ibz,spin)
       call wfd_who_has_ug(Wfd,band,ik_ibz,spin,how_many,proc_ranks)
       if (how_many /= 1) then
         ierr = ierr + 1
         write(msg,'(a,3(i0,1x))')" Found replicated state (b,k,s) ",band,ik_ibz,spin
         MSG_WARNING(msg)
       end if
     end do
   end do
 end do

 if (ierr /= 0) then
   MSG_ERROR("Cannot write WFK file when wavefunctions are replicated")
 end if

 call cwtime(cpu,wall,gflops,"start")

 ! Master node opens the file and writes the Abinit header.
 if (iam_master) then
   call wfkfile%open_write(Hdr,wfk_fname,formeig0,iomode,get_unit(),xmpi_comm_self,write_hdr=.TRUE.,write_frm=.FALSE.)
 end if

 ! Other nodes wait here before opening the same file.
 call xmpi_barrier(Wfd%comm)
 if (.not.iam_master) then
   call wfkfile%open_write(Hdr,wfk_fname,formeig0,iomode,get_unit(),xmpi_comm_self,write_hdr=.FALSE.,write_frm=.FALSE.)
 end if

 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz
     if (.not. wfd%ihave_ug(band, ik_ibz, spin, how="Stored")) cycle
     !if (.not. wave%has_ug == WFD_STORED) then
     nband_k = Wfd%nband(ik_ibz,spin)
     npw_k   = Wfd%npwarr(ik_ibz)

     ! Compute my block of bands for this k-point and spin.
     call wfd%mybands(ik_ibz, spin, how_manyb, my_band_list, how="Stored")
     call list2blocks(my_band_list(1:how_manyb), nblocks, blocks)

     !if (proc_distrb_cycle(mpi_enreg%proc_distrb,ik_ibz,1,nband_k,spin,my_rank)) CYCLE
     !call mask2blocks(mpi_enreg%proc_distrb(ik_ibz,:,spin)==my_rank, nblocks,blocks)

     ABI_CHECK(nblocks==1,"nblocks !=1")
     write(msg,"(a,3(i0,2x))")"Will write (ik_ibz, spin, nblocks)",ik_ibz,spin,nblocks
     call wrtout(std_out, msg)

     ! Extract the block of wavefunctions from Wfd.
     ! Try to allocate all u(g) first,
     ! TODO If not enough memory fallback to a blocked algorithm.
     cgsize = Wfd%nspinor * npw_k * how_manyb
     ABI_MALLOC_OR_DIE(cg_k, (2,cgsize), ierr)

     ! Extract the set of u(g) for this (kpoint,spin)
     ! This works only if all the bands are on the same node.
     !band_block = [1, nband_k]
     !call wfd_extract_cgblock(Wfd,[(ii, ii=1,nband_k)],ik_ibz,spin,cg_k)

     do blk=1,nblocks
       band_block = blocks(:,blk)
       call wfd_extract_cgblock(Wfd,[(ii, ii=band_block(1),band_block(2))],ik_ibz,spin,cg_k)

       if (band_block(1)==1) then
         ! Write also kg_k, eig_k and occ_k
         call wfkfile%write_band_block(band_block,ik_ibz,spin,xmpio_single,&
            kg_k=Wfd%Kdata(ik_ibz)%kg_k,cg_k=cg_k, &
            eig_k=Bands%eig(:,ik_ibz,spin),occ_k=Bands%occ(:,ik_ibz,spin))
       else
         MSG_ERROR("This should not happen in the present version!")
         !call wfkfile%write_band_block(band_block,ik_ibz,spin,xmpio_single,cg_k=cg_k(:,1+icg:))
       end if
     end do

     ABI_FREE(cg_k)
     ABI_FREE(blocks)
   end do
 end do

 ! Close the file.
 call wfkfile%close()

 call cwtime_report(" write all cg" , cpu, wall, gflops)

 DBG_EXIT("COLL")

end subroutine wfd_write_wfk
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_read_wfk
!! NAME
!! wfd_read_wfk
!!
!! FUNCTION
!!  This routine reads the WFK file completing the initialization of the wavefunction descriptor
!!
!! INPUTS
!!  wfk_fname=Name of the WFK file.
!!  iomode=Option specifying the fileformat as well as the IO mode to be used.
!!
!! OUTPUT
!!  [out_hdr]=Header of the WFK file.
!!
!! SIDE EFFECTS
!!  Wfd<wfd_t>=All the states owned by this node whose status is (STORED|ALLOCATED) read.
!!
!! PARENTS
!!      bethe_salpeter,m_gkk,m_phgamma,m_phpi,m_sigmaph,m_wfd,screening,sigma
!!      wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_read_wfk(Wfd, wfk_fname, iomode, out_hdr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iomode
 character(len=*),intent(in) :: wfk_fname
 class(wfd_t),target,intent(inout) :: Wfd
 type(Hdr_type),optional,intent(out) :: out_hdr

!Local variables ------------------------------
!scalars
 integer,parameter :: formeig0=0, optkg1=1, method = 2
 integer :: wfk_unt,npw_disk,nmiss,ig,sc_mode,ii, enough
 integer :: comm,master,my_rank,spin,ik_ibz,fform,ierr ! ,igp
 integer :: mcg,nband_wfd,nband_disk,band,mband_disk,bcount,istwfk_disk
 integer :: spinor,cg_spad,gw_spad,icg,igw,cg_bpad,ib, ik, is
 logical :: change_gsphere
 real(dp) :: cpu, wall, gflops, cpu_ks, wall_ks, gflops_ks
 character(len=500) :: msg
 type(Wfk_t) :: Wfk
 type(Hdr_type) :: Hdr
 type(wave_t),pointer :: wave
!arrays
 integer,allocatable :: gf2wfd(:),kg_k(:,:)
 integer :: work_ngfft(18),gmax_wfd(3),gmax_disk(3),gmax(3), all_countks(wfd%nkibz, wfd%nsppol)
 real(dp),allocatable :: eig_k(:),cg_k(:,:) !occ_k(:),
 real(dp),allocatable :: out_cg(:,:), work(:,:,:,:)
 logical,allocatable :: my_readmask(:,:,:)
 character(len=6) :: tag_spin(2)

!************************************************************************

 DBG_ENTER("COLL")

 if (any(iomode == [IO_MODE_NETCDF, IO_MODE_FORTRAN_MASTER])) then
   MSG_ERROR(sjoin("Unsupported value for iomode: ",itoa(iomode)))
 end if

 sc_mode = xmpio_collective
 comm = Wfd%comm; my_rank = Wfd%my_rank; master = Wfd%master

 tag_spin(:)=(/'      ','      '/); if (Wfd%nsppol==2) tag_spin(:)=(/' UP   ',' DOWN '/)
 call wrtout(std_out, sjoin(" wfd_read_wfk: reading", wfk_fname, "with iomode:", iomode2str(iomode)), &
             do_flush=.True., pre_newlines=2)

 wfk_unt = get_unit()
 call wfk_open_read(Wfk,wfk_fname,formeig0,iomode,wfk_unt,Wfd%comm,Hdr_out=Hdr)
 if (present(out_hdr)) call hdr_copy(hdr, out_hdr)

 ! TODO: Perform consistency check btw Hdr and Wfd.
 ! Output the header of the GS wavefunction file.
 fform = 0
 if (wfd%prtvol /= 0 .and. wfd%my_rank == 0) call hdr%echo(fform, 4, unit=std_out)

 mband_disk = MAXVAL(Hdr%nband)
 ABI_CHECK_ILEQ(Wfd%mband, mband_disk, "Not enough bands stored on WFK file")

 ! Each node will read the waves whose status if (WFD_ALLOCATED|WFD_STORED).
 ! all_countks is a global array used to skip (ik_ibz, spin) if all MPI procs do not need bands for this (k, s)
 ABI_MALLOC(my_readmask,(mband_disk, Wfd%nkibz, Wfd%nsppol))
 my_readmask=.FALSE.
 all_countks = 0; enough = 0
 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz
     do band=1,Wfd%nband(ik_ibz,spin)
       if (wfd%ihave_ug(band, ik_ibz, spin)) then
         my_readmask(band,ik_ibz,spin) = .TRUE.
         all_countks(ik_ibz, spin) = 1
         if (wfd%ihave_ug(band, ik_ibz, spin, how="Stored")) then
           enough = enough + 1
           if (enough < 30) then
             MSG_WARNING("Wavefunction is already stored!")
           end if
         end if
       end if
     end do
   end do
 end do

 ! All procs must agree when skipping (k, s)
 call xmpi_sum(all_countks, wfd%comm, ierr)

 write(msg,'(a,i0,a)')" Reading ",COUNT(my_readmask)," (b,k,s) states ..."
 call wrtout(std_out, msg)
 if (wfd%prtvol > 0) call wrtout(std_out,' k       eigenvalues [eV]','COLL')
 call cwtime(cpu, wall, gflops, "start")

 if (method == 1) then
  do spin=1,Wfd%nsppol
    do ik_ibz=1,Wfd%nkibz
      if (all_countks(ik_ibz, spin) == 0) cycle
      npw_disk   = Hdr%npwarr(ik_ibz)
      nband_disk = Hdr%nband(ik_ibz+(spin-1)*Hdr%nkpt)
      istwfk_disk = hdr%istwfk(ik_ibz)
      change_gsphere = istwfk_disk /= wfd%istwfk(ik_ibz)
      ABI_CHECK(.not. change_gsphere, "different istwfk values are not coded")

      nband_wfd  = Wfd%nband(ik_ibz,spin)
      if (nband_wfd > nband_disk) then
        write(msg,'(a,2(i0,1x))')&
         " nband_wfd to be read cannot be greater than nband_disk while: ",nband_wfd,nband_disk
        MSG_ERROR(msg)
      end if

      mcg = npw_disk*Wfd%nspinor*nband_wfd

      ABI_MALLOC(eig_k,((2*Wfk%mband)**formeig0*Wfk%mband))

      ABI_MALLOC(kg_k,(3,optkg1*npw_disk))
      ABI_MALLOC_OR_DIE(cg_k,(2,mcg), ierr)

      call wfk%read_band_block([1,nband_wfd] , ik_ibz, spin, sc_mode, kg_k=kg_k, cg_k=cg_k, eig_k=eig_k)

      if (wfd%prtvol > 0 .and. Wfd%my_rank==Wfd%master) then
        if (Wfd%nsppol==2) then
          write(std_out,'(i3,a,10f7.2/50(10x,10f7.2/))') ik_ibz,tag_spin(spin),(eig_k(ib)*Ha_eV,ib=1,nband_wfd)
        else
          write(std_out,'(i3,7x,10f7.2/50(10x,10f7.2/))')ik_ibz,(eig_k(ib)*Ha_eV,ib=1,nband_wfd)
        end if
      end if

      ! Table with the correspondence btw the k-centered sphere of the WFK file
      ! and the one used in Wfd (possibly smaller due to ecutwfn).
      ABI_MALLOC(gf2wfd,(npw_disk))
      if (any(my_readmask(:,ik_ibz,spin))) then
        call kg_map(wfd%npwarr(ik_ibz), wfd%kdata(ik_ibz)%kg_k, npw_disk, kg_k, gf2wfd, nmiss)
      end if
      !if (nmiss/=0) then
      !  write(msg,'(a,2(1x,i0),a,i0)')" For (k,s) ",ik_ibz,spin," the number of missing G is ",nmiss
      !  MSG_WARNING(msg)
      !end if

      ! Conversion of the basis set.
      do band=1,Wfd%nband(ik_ibz,spin)
        if (my_readmask(band, ik_ibz, spin)) then

          ABI_CHECK(all(wfd%bks2wfd(:, band, ik_ibz, spin) /= 0), "state in not allocated")
          ib = wfd%bks2wfd(1, band, ik_ibz, spin)
          ik = wfd%bks2wfd(2, band, ik_ibz, spin)
          is = wfd%bks2wfd(3, band, ik_ibz, spin)
          wave => wfd%s(is)%k(ik)%b(ib)
          wave%ug = czero

          cg_bpad=npw_disk*Wfd%nspinor*(band-1)
          do spinor=1,Wfd%nspinor
            cg_spad=(spinor-1)*npw_disk
            gw_spad=(spinor-1)*Wfd%npwarr(ik_ibz)
            do ig=1,npw_disk
              icg = ig+cg_spad+cg_bpad
              igw = gf2wfd(ig)+gw_spad
              if (gf2wfd(ig) /= 0) then
                wave%ug(igw) = CMPLX(cg_k(1,icg), cg_k(2,icg), kind=gwpc)
              end if
            end do
          end do
          wave%has_ug = WFD_STORED

        end if
      end do

      ABI_FREE(eig_k)
      ABI_FREE(kg_k)
      ABI_FREE(cg_k)
      ABI_FREE(gf2wfd)
    end do !ik_ibz
  end do !spin

 else if (method==2) then
  ! DEFAULT ALGO: This seems to be the most efficient one.

  do spin=1,Wfd%nsppol
    do ik_ibz=1,Wfd%nkibz
      if (all_countks(ik_ibz, spin) == 0) cycle
      call cwtime(cpu_ks, wall_ks, gflops_ks, "start")
      !write(std_out,*)"about to read ik_ibz: ",ik_ibz,", spin: ",spin
      npw_disk   = Hdr%npwarr(ik_ibz)
      nband_disk = Hdr%nband(ik_ibz+(spin-1)*Hdr%nkpt)
      istwfk_disk = hdr%istwfk(ik_ibz)
      change_gsphere = istwfk_disk /= wfd%istwfk(ik_ibz)

      nband_wfd  = Wfd%nband(ik_ibz,spin)

      if (nband_wfd > nband_disk) then
        write(msg,'(a,2(i0,1x))')&
         "nband_wfd to be read cannot be greater than nband_disk while: ",nband_wfd,nband_disk
        MSG_ERROR(msg)
      end if

      ABI_MALLOC(eig_k,((2*nband_disk)**formeig0*nband_disk))
      ABI_MALLOC(kg_k,(3,optkg1*npw_disk))

      mcg = npw_disk*Wfd%nspinor*COUNT(my_readmask(:,ik_ibz,spin))
      ABI_MALLOC_OR_DIE(cg_k,(2,mcg), ierr)

      call wfk%read_bmask(my_readmask(:,ik_ibz,spin),ik_ibz,spin,sc_mode,kg_k=kg_k,cg_k=cg_k,eig_k=eig_k)

      if (Wfd%my_rank == Wfd%master .and. wfd%prtvol > 0) then
        if (Wfd%nsppol==2) then
          write(std_out,'(i3,a,10f7.2/50(10x,10f7.2/))') ik_ibz,tag_spin(spin),(eig_k(ib)*Ha_eV,ib=1,nband_wfd)
        else
          write(std_out,'(i3,7x,10f7.2/50(10x,10f7.2/))')ik_ibz,(eig_k(ib)*Ha_eV,ib=1,nband_wfd)
        end if
      end if

      ! Table with the correspondence btw the k-centered sphere of the WFK file
      ! and the one used in Wfd (possibly smaller due to ecutwfn).
      ! TODO: Here I should treat the case in which istwfk in wfd differs from the one on disk.
      ABI_MALLOC(gf2wfd,(npw_disk))
      if (any(my_readmask(:,ik_ibz,spin))) then
        call kg_map(wfd%npwarr(ik_ibz), wfd%kdata(ik_ibz)%kg_k, npw_disk, kg_k, gf2wfd, nmiss)
      end if
      !if (nmiss/=0) then
      !  write(msg,'(a,2(1x,i0),a,i0)')" For (k,s) ",ik_ibz,spin," the number of missing G is ",nmiss
      !  MSG_WARNING(msg)
      !end if

      if (change_gsphere .and. any(my_readmask(:,ik_ibz,spin))) then
        ! Prepare call to ctgk_change_sphere
        ! FFT box must enclose the two spheres (wfd(k), wfk(k))
        gmax_wfd = maxval(abs(wfd%kdata(ik_ibz)%kg_k), dim=2)
        gmax_disk = maxval(abs(kg_k), dim=2)
        do ii=1,3
          gmax(ii) = max(gmax_wfd(ii), gmax_disk(ii))
        end do
        gmax = 2*gmax + 1
        call ngfft_seq(work_ngfft, gmax)
        ABI_MALLOC(work, (2, work_ngfft(4),work_ngfft(5),work_ngfft(6)))
        ABI_MALLOC(out_cg, (2, wfd%npwarr(ik_ibz) * wfd%nspinor))
     end if

      ! Conversion of the basis set.
      bcount = 0
      do band=1,Wfd%nband(ik_ibz,spin)
        if (my_readmask(band,ik_ibz,spin)) then

          ib = wfd%bks2wfd(1, band, ik_ibz, spin)
          ik = wfd%bks2wfd(2, band, ik_ibz, spin)
          is = wfd%bks2wfd(3, band, ik_ibz, spin)
          ABI_CHECK(all(wfd%bks2wfd(:, band, ik_ibz, spin) /= 0), "state in not allocated")

          wave => wfd%s(is)%k(ik)%b(ib)
          wave%ug = czero

          bcount = bcount + 1
          cg_bpad=npw_disk*Wfd%nspinor*(bcount-1)

          if (change_gsphere) then
            ! Different istwfk storage.
            call cgtk_change_gsphere(wfd%nspinor, &
               npw_disk, istwfk_disk, kg_k, cg_k(:, cg_bpad+1:), &
               wfd%npwarr(ik_ibz), wfd%istwfk(ik_ibz), wfd%kdata(ik_ibz)%kg_k, out_cg, work_ngfft, work)

            wave%ug(:) = CMPLX(out_cg(1, :), out_cg(2, :), kind=gwpc)
            !call wfd%push_ug(band, ik_ibz, spin, cryst, out_cg)
          else
            do spinor=1,Wfd%nspinor
              cg_spad=(spinor-1)*npw_disk
              gw_spad=(spinor-1)*Wfd%npwarr(ik_ibz)
              do ig=1,npw_disk
                icg = ig+cg_spad+cg_bpad
                igw = gf2wfd(ig)+gw_spad
                if (gf2wfd(ig) /= 0) then
                  wave%ug(igw) = CMPLX(cg_k(1,icg),cg_k(2,icg), kind=gwpc)
                end if
              end do
              !call wfd%push_ug(band, ik_ibz, spin, cryst, out_cg)
            end do
          end if

          wave%has_ug = WFD_STORED
        end if
      end do

      ABI_FREE(eig_k)
      ABI_FREE(kg_k)
      ABI_FREE(cg_k)
      ABI_FREE(gf2wfd)
      ABI_SFREE(work)
      ABI_SFREE(out_cg)

      if (ik_ibz <= 10 .or. mod(ik_ibz, 200) == 0) then
        write(msg,'(4(a,i0),a)') " Reading k-point [", ik_ibz, "/", wfd%nkibz, "] spin [", spin, "/", wfd%nsppol, "]"
        call cwtime_report(msg, cpu_ks, wall_ks, gflops_ks)
      end if
    end do !ik_ibz
  end do !spin

 else
   MSG_ERROR(sjoin("Wrong method: ", itoa(method)))
 end if

 call wfk%close()
 call Hdr%free()

 ABI_FREE(my_readmask)

 ! Update the kbs table storing the distribution of the ug and set the MPI communicators.
 call wfd%set_mpicomm()
 !call wfd%update_bkstab()

 call cwtime_report(" WFK IO", cpu, wall, gflops, end_str=ch10)

 DBG_EXIT("COLL")

end subroutine wfd_read_wfk
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_paw_get_aeur
!! NAME
!! wfd_paw_get_aeur
!!
!! FUNCTION
!!   Compute the AE PAW wavefunction in real space.
!!
!! INPUTS
!!   band,ik_ibz,spin=indices specifying the band, the k-point and the spin.
!!   Psps<pseudopotential_type>=variables related to pseudopotentials
!!   Cryst<crystal_t>= data type gathering info on symmetries and unit cell.
!!   Wfd<wfd_t>=wavefunction descriptor.
!!   Pawtab(ntypat*usepaw)<type(pawtab_type)>=paw tabulated starting data.
!!   Pawfgrtab(natom)<pawfgrtab_type>= atomic data given on fine rectangular grid.
!!     NB: rpaw should be used in nhatgrid to initialize the datatype (optcut=1 option) instead of the radius for the
!!     shape functions (rpaw /= rshp).
!!   Paw_onsite(natom)<paw_pwaves_lmn_t>=3D PAW partial waves in real space for each FFT point in the PAW spheres.
!!
!! OUTPUT
!! ur_ae(Wfd%nfft*Wfd%nspinor)=AE PAW wavefunction in real space.
!! [ur_ae_onsite(Wfd%nfft*Wfd%nspinor)]
!! [ur_ps_onsite(Wfd%nfft*Wfd%nspinor)]
!!
!! NOTES
!!  (1) The true wavefunction integrates in real space to the unit cell volume.
!!      The definition of the cprj matrix elements includes the term 1/SQRT(ucvol) that comes
!!      from the use of a normalized planewave e^(iG.r)/SQRT(omega) in the FFT transform G-->R (see e.g. opernla_ylm)
!!      On the contrary, the convention for the G-->R transform employed in the FFT routines used in abinit is
!!      u(r) = sum_G u(G) e^(iG.r); u(G) = one/omega \int u(r) e^(-iG.r)dr.
!!      Hence we have to multiply the onsite part by SQRT(uvol) before adding the smooth FFT part in real space.
!!
!!  (2) Care has to be taken in the calculation of the onsite contribution when the FFT point belongs to the PAW
!!      sphere of a periodically repeated atom. In this case one evaluates the onsite term associated to the
!!      atom in the first unit cell then the contribution has to be multiplied by a k- dependent
!!      phase factor to account for the wrapping of the real-space point in the first unit cell.
!!
!! PARENTS
!!      calc_sigc_me,calc_sigx_me,cchi0,cchi0q0,classify_bands,m_wfd
!!      prep_calc_ucrpa,wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_paw_get_aeur(Wfd,band,ik_ibz,spin,Cryst,Paw_onsite,Psps,Pawtab,Pawfgrtab,ur_ae,ur_ae_onsite,ur_ps_onsite)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin
 type(pseudopotential_type),intent(in) :: Psps
 type(crystal_t),intent(in) :: Cryst
 class(wfd_t),intent(inout) :: Wfd
!arrays
 type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat)
 type(pawfgrtab_type),intent(in) :: Pawfgrtab(Cryst%natom)
 type(paw_pwaves_lmn_t),intent(in) :: Paw_onsite(Cryst%natom)
 complex(gwpc),intent(out) :: ur_ae(Wfd%nfft*Wfd%nspinor)
 complex(gwpc),optional,intent(out) :: ur_ae_onsite(Wfd%nfft*Wfd%nspinor)
 complex(gwpc),optional,intent(out) :: ur_ps_onsite(Wfd%nfft*Wfd%nspinor)

!Local variables-------------------------------
!scalars
 integer :: itypat,ln_size,lmn_size,iatom,spinor
 integer :: nfgd,ifgd,jlmn,jl,jm,ifftsph
 real(dp) :: phj,tphj,arg,re_cp,im_cp
 complex(dpc) :: cp,cnorm
!arrays
 real(dp) :: kpoint(3)
 complex(dpc),allocatable :: ceikr(:),phk_atm(:)
 type(pawcprj_type),allocatable :: Cp1(:,:)

! *************************************************************************

 ! TODO ngfft should be included in pawfgrtab_type
 !% if (ANY(Wfd%ngfft(1:3)/=Pawfgrtab%ngfft(1:3)) then
 !!  MSG_ERROR("Wfd%ngfft(1:3)/=Pawfgrtab%ngfft(1:3)")
 !% end if

 call wfd%get_ur(band,ik_ibz,spin,ur_ae)

 kpoint = Wfd%kibz(:,ik_ibz)

 ABI_MALLOC(ceikr,(Wfd%nfftot))

 call calc_ceikr(kpoint,Wfd%nfftot,Wfd%ngfft,ceikr)
 ur_ae = ur_ae * ceikr

 ABI_MALLOC(Cp1,(Wfd%natom,Wfd%nspinor))
 call pawcprj_alloc(Cp1,0,Wfd%nlmn_atm)

 call wfd%get_cprj(band,ik_ibz,spin,Cryst,Cp1,sorted=.FALSE.)

 ! Add onsite term on the augmented FFT mesh.
 if (present(ur_ae_onsite)) ur_ae_onsite = czero
 if (present(ur_ps_onsite)) ur_ps_onsite = czero

 ABI_CHECK(Wfd%nspinor==1,"nspinor==1 not coded")

 do iatom=1,Cryst%natom
   itypat  =Cryst%typat(iatom)
   lmn_size=Pawtab(itypat)%lmn_size
   ln_size =Pawtab(itypat)%basis_size   ! no. of nl elements in PAW basis.
   nfgd    =Pawfgrtab(iatom)%nfgd       ! no. of points in the fine grid for this PAW sphere.

   ABI_MALLOC(phk_atm,(nfgd))
   do ifgd=1,nfgd
     arg = -two_pi* DOT_PRODUCT(Paw_onsite(iatom)%r0shift(:,ifgd),kpoint)
     phk_atm(ifgd) = DCMPLX(COS(arg),SIN(arg))
   end do

   do spinor=1,Wfd%nspinor
     do jlmn=1,lmn_size
       jl=Psps%indlmn(1,jlmn,itypat)
       jm=Psps%indlmn(2,jlmn,itypat)
       re_cp = Cp1(iatom,spinor)%cp(1,jlmn)
       im_cp = Cp1(iatom,spinor)%cp(2,jlmn)
       cp = DCMPLX(re_cp, im_cp) * SQRT(Cryst%ucvol) ! Pay attention here. see (1).

       do ifgd=1,nfgd ! loop over fine grid points in current PAW sphere.
         ifftsph = Pawfgrtab(iatom)%ifftsph(ifgd) ! index of the point on the grid
         phj  = Paw_onsite(iatom)% phi(ifgd,jlmn)
         tphj = Paw_onsite(iatom)%tphi(ifgd,jlmn)
         ur_ae(ifftsph)           = ur_ae(ifftsph) + cp * (phj-tphj) * phk_atm(ifgd)
         if (present(ur_ae_onsite)) ur_ae_onsite(ifftsph) = ur_ae_onsite(ifftsph) + cp *  phj * phk_atm(ifgd)
         if (present(ur_ps_onsite)) ur_ps_onsite(ifftsph) = ur_ps_onsite(ifftsph) + cp * tphj * phk_atm(ifgd)
       end do
     end do !jlmn
   end do !spinor

   ABI_FREE(phk_atm)
 end do !iatom

 ! Remove the phase e^{ikr}, u(r) is returned.
 ur_ae = ur_ae * CONJG(ceikr)
 cnorm = xdotc(Wfd%nfft*Wfd%nspinor,ur_ae,1,ur_ae,1)/Wfd%nfft
 !write(std_out,*)" AE PAW norm: (b,k,s)",band,ik_ibz,spin,REAL(cnorm)

 if (present(ur_ae_onsite)) ur_ae_onsite = ur_ae_onsite * CONJG(ceikr)
 if (present(ur_ps_onsite)) ur_ps_onsite = ur_ps_onsite * CONJG(ceikr)

 call pawcprj_free(Cp1)
 ABI_FREE(Cp1)
 ABI_FREE(ceikr)

end subroutine wfd_paw_get_aeur
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_plot_ur
!! NAME
!! wfd_plot_ur
!!
!! FUNCTION
!!  This routine writes the squared modulus of the wavefunctions in real space
!!  to an external files, one for each (k,b,s). File are written in the XSF format (Xcrysden).
!!  A subset of (b,k,s) states can be specified via the bks_mask. The routine is MPI parallelized.
!!
!! INPUTS
!!  Wfd<wfd_t>=Wavefunction descriptor.
!!  Cryst<crystal_t>= Information on symmetries and unit cell.
!!  Psps<pseudopotential_type>=Pseudopotential info.
!!  Pawtab(ntypat*usepaw)<type(pawtab_type)>=PAW tabulated starting data.
!!  Pawrad(ntypat*usepaw)<type(pawrad_type)>=paw radial mesh and related data.
!!  ngfftf(18)=The FFT mesh used for plotting |wfr|**2, it can differ from the one internally used in Wfd.
!!    For example, PAW wavefunctions should be plotted on a much finer FFT mesh.
!!  bks_mask(mband,nkibz,nsppol)=logical mask used to select the states to be plotted.
!!
!! OUTPUT
!!  Output is written on file.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_plot_ur(Wfd,Cryst,Psps,Pawtab,Pawrad,ngfftf,bks_mask)

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: Cryst
 type(Pseudopotential_type),intent(in) :: Psps
 class(wfd_t),intent(inout) :: Wfd
!arrays
 integer,intent(in) :: ngfftf(18)
 logical,target,intent(in) :: bks_mask(Wfd%mband,Wfd%nkibz,Wfd%nsppol)
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
 type(Pawrad_type),intent(in) :: Pawrad(Cryst%ntypat*Wfd%usepaw)

!Local variables ------------------------------
!scalars
 integer :: spin,band,ik_ibz,optcut,optgr0,optgr1,optgr2,optrad
 integer :: n1,n2,n3,my_nplots,plot,funt,my_nband,cplex
 character(len=500) :: msg
 character(len=fnlen) :: xsf_fname
!arrays
 integer :: got(Wfd%nproc)
 integer,allocatable :: l_size_atm(:),my_plot_list(:,:)
 integer :: my_band_list(Wfd%mband)
 real(dp),allocatable :: data_plot(:)
 logical,ABI_CONTIGUOUS pointer :: bmask(:)
 complex(gwpc),allocatable :: ur_ae(:),nc_ur(:)
 type(Pawfgrtab_type),allocatable :: Pawfgrtab(:)
 type(paw_pwaves_lmn_t),allocatable :: Paw_onsite(:)

!************************************************************************

 if (ALL(.not.bks_mask)) RETURN

 DBG_ENTER("COLL")

 call wrtout(std_out," Plotting |wfs|^2 ...")
 !
 ! Change the FFT mesh if needed because we want u(r) on the ngfftf mesh (pawecutd for PAW).
 call wfd%change_ngfft(Cryst,Psps,ngfftf)
 n1 = ngfftf(1); n2 = ngfftf(2); n3 = ngfftf(3)

 ! Distribute the plots among the nodes taking into account the distribution of the waves.
 ! my_plot_list gives the list of (b,k,s) states plotted by this node.
 ABI_MALLOC(my_plot_list,(3,Wfd%mband*Wfd%nkibz*Wfd%nsppol))

 my_nplots=0; got=0
 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz
     bmask => bks_mask(:,ik_ibz,spin)
     call wfd%distribute_bands(ik_ibz,spin,my_nband,my_band_list,got,bmask)

     if (my_nband>0) then
       my_plot_list(1,my_nplots+1:my_nplots+my_nband) = my_band_list(1:my_nband)
       my_plot_list(2,my_nplots+1:my_nplots+my_nband) = ik_ibz
       my_plot_list(3,my_nplots+1:my_nplots+my_nband) = spin
       my_nplots = my_nplots + my_nband
     end if
   end do
 end do

 if (Wfd%usepaw==1) then
   MSG_WARNING("Testing the calculation of AE PAW wavefunctions.")
   ! Use a local pawfgrtab to make sure we use the correction in the paw spheres
   ! the usual pawfgrtab uses r_shape which may not be the same as r_paw.
   cplex=1
   call pawtab_get_lsize(Pawtab,l_size_atm,Cryst%natom,Cryst%typat)
   ABI_MALLOC(Pawfgrtab,(Cryst%natom))
   call pawfgrtab_init(Pawfgrtab,cplex,l_size_atm,Wfd%nspden,Cryst%typat)
   ABI_FREE(l_size_atm)

   optcut=1                     ! use rpaw to construct Pawfgrtab.
   optgr0=0; optgr1=0; optgr2=0 ! dont need gY terms locally.
   optrad=1                     ! do store r-R.

   call nhatgrid(Cryst%atindx1,Cryst%gmet,Cryst%natom,Cryst%natom,Cryst%nattyp,ngfftf,Cryst%ntypat,&
     optcut,optgr0,optgr1,optgr2,optrad,Pawfgrtab,Pawtab,Cryst%rprimd,Cryst%typat,Cryst%ucvol,Cryst%xred)

   !Pawfgrtab is ready to use

   if (Wfd%pawprtvol>0) then
     call pawfgrtab_print(Pawfgrtab,natom=Cryst%natom,unit=std_out,&
                          prtvol=Wfd%pawprtvol,mode_paral="COLL")
   end if

   ABI_MALLOC(Paw_onsite,(Cryst%natom))
   call paw_pwaves_lmn_init(Paw_onsite,Cryst%natom,Cryst%natom,Cryst%ntypat,&
                            Cryst%rprimd,Cryst%xcart,Pawtab,Pawrad,Pawfgrtab)

   ABI_MALLOC(ur_ae,(Wfd%nfft*Wfd%nspinor))
   ABI_MALLOC(data_plot,(Wfd%nfft))

   do plot=1,my_nplots
     band  =my_plot_list(1,plot)
     ik_ibz=my_plot_list(2,plot)
     spin  =my_plot_list(3,plot)

     call wfd%paw_get_aeur(band,ik_ibz,spin,Cryst,Paw_onsite,Psps,Pawtab,Pawfgrtab,ur_ae)

     data_plot = DBLE(ur_ae(1:Wfd%nfft)*CONJG(ur_ae(1:Wfd%nfft)))/Cryst%ucvol
     if (Wfd%nspinor==2) data_plot = data_plot + DBLE(ur_ae(Wfd%nfft+1:)*CONJG(ur_ae(Wfd%nfft+1:)))/Cryst%ucvol

     write(xsf_fname,'(3(a,i0),a)') 'PAW_AE_wfk2_sp',spin,'_kpt',ik_ibz,'_bd',band,'.xsf'
     if (open_file(xsf_fname,msg,newunit=funt,status='unknown',form='formatted') /= 0) then
       MSG_ERROR(msg)
     end if

     call printxsf(n1,n2,n3,data_plot,Cryst%rprimd,(/zero,zero,zero/),&
       Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xcart,Cryst%znucl,funt,0)

     close(funt)
   end do

   ABI_FREE(ur_ae)
   ABI_FREE(data_plot)

   call pawfgrtab_free(Pawfgrtab)
   ABI_FREE(Pawfgrtab)
   call paw_pwaves_lmn_free(Paw_onsite)
   ABI_FREE(Paw_onsite)

 else
   ! NC case. Just a simple FFT G-->R and then dump the results.
   ABI_MALLOC(nc_ur,(Wfd%nfft*Wfd%nspinor))
   ABI_MALLOC(data_plot,(Wfd%nfft))

   do plot=1,my_nplots
     band  =my_plot_list(1,plot)
     ik_ibz=my_plot_list(2,plot)
     spin  =my_plot_list(3,plot)

     call wfd%get_ur(band,ik_ibz,spin,nc_ur)

     data_plot = DBLE(nc_ur(1:Wfd%nfft)*CONJG(nc_ur(1:Wfd%nfft)))/Cryst%ucvol
     if (Wfd%nspinor==2) data_plot = data_plot + DBLE(nc_ur(Wfd%nfft+1:)*CONJG(nc_ur(Wfd%nfft+1:)))/Cryst%ucvol

     write(xsf_fname,'(3(a,i0),a)') 'NC_wfk2_sp',spin,'_kpt',ik_ibz,'_bd',band,'.xsf'
     if (open_file(xsf_fname,msg,newunit=funt,status='unknown',form='formatted') /= 0) then
       MSG_ERROR(msg)
     end if
     call printxsf(n1,n2,n3,data_plot,Cryst%rprimd,(/zero,zero,zero/),&
       Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xcart,Cryst%znucl,funt,0)

     close(funt)
   end do

   ABI_FREE(nc_ur)
   ABI_FREE(data_plot)
 end if

 ABI_FREE(my_plot_list)

 DBG_EXIT("COLL")

end subroutine wfd_plot_ur
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/wfd_get_socpert
!! NAME
!! wfd_get_socpert
!!
!! FUNCTION
!!
!! INPUTS
!! cryst<crystal_t>= data type gathering info on symmetries and unit cell
!! psps<pseudopotential_type>=variables related to pseudopotentials
!! pawtab(psps%ntypat) <type(pawtab_type)>=paw tabulated starting data
!! paw_ij(natom)<type(paw_ij_type)>=data structure containing PAW arrays given on (i,j) channels.
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

!!!  subroutine wfd_get_socpert(wfd, cryst, psps, pawtab, bks_mask, osoc_bks)
!!!
!!!   !use m_pawcprj
!!!   use m_hamiltonian,    only : destroy_hamiltonian, init_hamiltonian, &
!!!                                load_spin_hamiltonian,load_k_hamiltonian, gs_hamiltonian_type
!!!
!!!   implicit none
!!!
!!!  !Arguments ------------------------------------
!!!  !scalars
!!!   type(wfd_t),target,intent(inout) :: wfd
!!!   type(crystal_t),intent(in) :: cryst
!!!   type(pseudopotential_type),intent(in) :: psps
!!!  ! arrays
!!!   logical,intent(in) :: bks_mask(wfd%mband, wfd%nkibz, wfd%nsppol)
!!!   real(dp),allocatable,intent(out) :: osoc_bks(:, :, :)
!!!   type(Pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
!!!   !type(paw_ij_type),intent(in) :: paw_ij(cryst%natom*psps%usepaw)
!!!
!!!  !Local variables ------------------------------
!!!  !scalars
!!!   integer,parameter :: nspinor2=2,nspden4=4,nsppol1=1,spin1=1
!!!   integer,parameter :: ndat1=1,nnlout0=0,tim_nonlop0=0,idir0=0 !,ider0=0,
!!!   integer :: natom,band,spin,ik_ibz,npw_k,istwf_k,nkpg !,ig,optder,matblk,mkmem_,nkpg,dimffnl,nspinortot
!!!   integer :: choice,cpopt,cp_dim,paw_opt,signs,ierr
!!!   !character(len=500) :: msg
!!!   type(gs_hamiltonian_type) :: ham_k
!!!  !arrays
!!!   integer :: bks_distrb(wfd%mband, wfd%nkibz, wfd%nsppol)
!!!   integer, ABI_CONTIGUOUS pointer :: kg_k(:,:)
!!!   !real(dp) :: kptns_(3,1),ylmgr_dum(1,1,1),shifts(3)
!!!   !real(dp),allocatable :: ylm_k(:,:),dum_ylm_gr_k(:,:,:)
!!!   !real(dp),pointer :: ffnl_k(:,:,:,:)
!!!   real(dp) :: kpoint(3),dum_enlout(0),dummy_lambda(1),soc(2)
!!!   real(dp),allocatable :: kpg_k(:,:),vnl_psi(:,:),vectin(:,:) !,s_psi(:,:)
!!!   real(dp),allocatable :: opaw_psi(:,:) !2, npw_k*wfd%nspinor*wfd%usepaw) ! <G|1+S|Cnk>
!!!   real(dp),ABI_CONTIGUOUS pointer :: ffnl_k(:,:,:,:),ph3d_k(:,:,:)
!!!   type(pawcprj_type),allocatable :: cprj(:,:)
!!!
!!!  !************************************************************************
!!!
!!!   DBG_ENTER("COLL")
!!!   ABI_CHECK(wfd%paral_kgb == 0, "paral_kgb not coded")
!!!
!!!   natom = cryst%natom
!!!
!!!   signs  = 2  ! => apply the non-local operator to a function in G-space.
!!!   choice = 1  ! => <G|V_nonlocal|vectin>.
!!!   cpopt  =-1; paw_opt= 0
!!!   if (wfd%usepaw==1) then
!!!     paw_opt=4 ! both PAW nonlocal part of H (Dij) and overlap matrix (Sij)
!!!     cpopt=3   ! <p_lmn|in> are already in memory
!!!
!!!     cp_dim = ((cpopt+5) / 5)
!!!     ABI_MALLOC(cprj, (natom, nspinor2*cp_dim))
!!!     call pawcprj_alloc(cprj, 0, wfd%nlmn_sort)
!!!   end if
!!!
!!!   ! Initialize the Hamiltonian on the coarse FFT mesh.
!!!   call init_hamiltonian(ham_k, psps, pawtab, nspinor2, nsppol1, nspden4, natom, cryst%typat, cryst%xred, &
!!!      wfd%nfft, wfd%mgfft, wfd%ngfft, cryst%rprimd, wfd%nloalg)
!!!   !ham_k%ekb(:,:,1) = zero
!!!
!!!   ! Continue to prepare the GS Hamiltonian (note spin1)
!!!   call load_spin_hamiltonian(ham_k, spin1, with_nonlocal=.True.)
!!!
!!!   ! Distribute (b, k, s) states.
!!!   call wfd%bks_distrb(bks_distrb, bks_mask=bks_mask)
!!!
!!!   ABI_CALLOC(osoc_bks, (wfd%mband, wfd%nkibz, wfd%nsppol))
!!!   osoc_bks = zero
!!!
!!!   do spin=1,wfd%nsppol
!!!     do ik_ibz=1,wfd%nkibz
!!!       if (all(bks_distrb(:, ik_ibz, spin) /= wfd%my_rank)) cycle
!!!
!!!       kpoint = wfd%kibz(:, ik_ibz)
!!!       npw_k = wfd%Kdata(ik_ibz)%npw; istwf_k = wfd%istwfk(ik_ibz)
!!!       ABI_CHECK(istwf_k == 1, "istwf_k must be 1 if SOC term is computed with perturbation theory.")
!!!       kg_k => wfd%kdata(ik_ibz)%kg_k
!!!       ffnl_k => wfd%Kdata(ik_ibz)%fnl_dir0der0
!!!       ph3d_k => wfd%Kdata(ik_ibz)%ph3d
!!!
!!!       ABI_MALLOC(vectin, (2, npw_k * nspinor2))
!!!       ABI_MALLOC(vnl_psi, (2, npw_k * nspinor2))
!!!       !ABI_MALLOC(cvnl_psi, (npw_k * nspinor2))
!!!       !ABI_MALLOC(s_psi, (2, npw_k * nspinor2 * psps%usepaw))
!!!
!!!       ! Compute (k+G) vectors (only if psps%useylm=1)
!!!       nkpg = 3 * wfd%nloalg(3)
!!!       ABI_MALLOC(kpg_k, (npw_k, nkpg))
!!!       if (nkpg > 0) then
!!!         call mkkpg(kg_k, kpg_k, kpoint, nkpg, npw_k)
!!!       end if
!!!
!!!       ! Load k-dependent part in the Hamiltonian datastructure
!!!       !matblk = min(NLO_MINCAT, maxval(ham_k%nattyp)); if (wfd%nloalg(2) > 0) matblk = natom
!!!       !ABI_MALLOC(ph3d_k,(2, npw_k, matblk))
!!!       call load_k_hamiltonian(ham_k, kpt_k=kpoint, npw_k=npw_k, istwf_k=istwf_k, kg_k=kg_k, &
!!!                               kpg_k=kpg_k, ffnl_k=ffnl_k, ph3d_k=ph3d_k, compute_ph3d=(wfd%paral_kgb/=1))
!!!
!!!       ! THIS PART IS NEEDED FOR THE CALL TO opernl although some quantities won't be used.
!!!       ! Now I do things cleanly then we try to pass zero-sized arrays!
!!!       !ABI_MALLOC(ylm_k, (npw_k, psps%mpsang**2 * psps%useylm))
!!!       !if (psps%useylm == 1) then
!!!       !  kptns_(:,1) = k4intp; optder = 0; mkmem_ = 1
!!!       !  ABI_MALLOC(dum_ylm_gr_k,(npw_k,3+6*(optder/2),psps%mpsang**2))
!!!       !  ! Here mband is not used if paral_compil_kpt=0
!!!       !  call initylmg(cryst%gprimd, kg_k, kptns_, mkmem_, wfd%MPI_enreg, psps%mpsang, npw_k, [1], 1,&
!!!       !    [npw_k], 1, optder, cryst%rprimd, ylm_k, dum_ylm_gr_k)
!!!       !  ABI_FREE(dum_ylm_gr_k)
!!!       !end if
!!!
!!!       ! ========================================================
!!!       ! ==== Compute nonlocal form factors ffnl at all (k+G) ====
!!!       ! ========================================================
!!!       !dimffnl = 1 + ider0 ! Derivatives are not needed.
!!!       !ABI_MALLOC(ffnl_k, (npw_k, dimffnl, psps%lmnmax, psps%ntypat))
!!!       !! ffnl_k => Kdata%fnl_dir0der0
!!!       !call mkffnl(psps%dimekb, dimffnl, psps%ekb, ffnl_k, psps%ffspl, cryst%gmet, cryst%gprimd, ider0, idir0, psps%indlmn,&
!!!       !   kg_k, kpg_k, k4intp, psps%lmnmax, psps%lnmax, psps%mpsang, psps%mqgrid_ff, nkpg, npw_k, &
!!!       !   psps%ntypat, psps%pspso, psps%qgrid_ff, cryst%rmet, psps%usepaw, psps%useylm, ylm_k, ylmgr_dum)
!!!       !ABI_FREE(ylm_k)
!!!
!!!       ! Calculate <G|Vnl|psi> for this k-point
!!!       do band=1,wfd%nband(ik_ibz, spin)
!!!         if (bks_distrb(band, ik_ibz, spin) /= wfd%my_rank) cycle
!!!
!!!         ! Input wavefunction coefficients <G|Cnk>.
!!!         ! vectin, (2, npw_k * nspinor2))
!!!         if (spin == 1) then
!!!           vectin(1, 1:npw_k) = dble(wfd%wave(band, ik_ibz, spin)%ug)
!!!           vectin(2, 1:npw_k) = aimag(wfd%wave(band, ik_ibz, spin)%ug)
!!!           vectin(:, npw_k+1:) = zero
!!!         else
!!!           vectin(:, 1:npw_k) = zero
!!!           vectin(1, npw_k+1:) = dble(wfd%wave(band, ik_ibz, spin)%ug)
!!!           vectin(2, npw_k+1:) = aimag(wfd%wave(band, ik_ibz, spin)%ug)
!!!         end if
!!!
!!!         if (wfd%usepaw == 1) call wfd%get_cprj(band, ik_ibz, spin, cryst, cprj, sorted=.True.)
!!!
!!!         ! TODO: consistency check for only_SO
!!!         call nonlop(choice, cpopt, cprj, dum_enlout, ham_k, idir0, dummy_lambda, wfd%mpi_enreg, ndat1, nnlout0, &
!!!                     paw_opt, signs, opaw_psi, tim_nonlop0, vectin, vnl_psi, only_SO=1)
!!!
!!!         soc = cg_zdotc(npw_k * nspinor2, vectin, vnl_psi)
!!!         write(std_out,*)soc * Ha_eV, "for (b, k, s)",band, ik_ibz, spin
!!!         osoc_bks(band, ik_ibz, spin) = soc(1)
!!!       end do ! band
!!!
!!!       !ABI_FREE(ffnl_k)
!!!       !ABI_FREE(ph3d_k)
!!!       ABI_FREE(vectin)
!!!       ABI_FREE(vnl_psi)
!!!       ABI_FREE(kpg_k)
!!!       !ABI_FREE(cvnl_psi)
!!!       !ABI_FREE(s_psi)
!!!     end do ! ik_ibz
!!!   end do ! spin
!!!
!!!   call xmpi_sum(osoc_bks, wfd%comm, ierr)
!!!
!!!   call destroy_hamiltonian(ham_k)
!!!
!!!   if (wfd%usepaw == 1) then
!!!     call pawcprj_free(cprj)
!!!     ABI_FREE(cprj)
!!!   end if
!!!
!!!   DBG_EXIT("COLL")
!!!
!!!  end subroutine wfd_get_socpert
!!***

!!****f* m_wfd/wfd_mkrho
!! NAME
!! wfd_mkrho
!!
!! FUNCTION
!! Calculate the charge density on the fine FFT grid in real space.
!!
!! INPUTS
!!  ngfftf(18)=array containing all the information for the "fine" FFT.
!!  Cryst<crystal_t> Info on the crystalline structure
!!  optcalc=option for calculation. If =0 (default value) perform calculation
!!    of electronic density. If =1, perform calculation of kinetic energy density.
!!    In both cases, the result is returned in rhor.
!!  Psps<type(pseudopotential_type)>=variables related to pseudopotentials
!!  nfftf=Total number of points on the fine FFT grid (for this processor)
!!  Kmesh<kmesh_t>= Info on the k-sampling:
!!  Wfd<wfd_t)=datatype gathering info on the wavefunctions.
!! [optcalc]=Optional option used to calculate the kinetic energy density. Defaults to 0.
!!
!! OUTPUT
!!  rhor(nfftf,nspden)=The density in the real space on the fine FFT grid.
!!   If nsppol==2, total charge in first half, spin-up component in second half.
!!   (for non-collinear magnetism, first element: total density, 3 next ones: mx,my,mz in units of hbar/2)
!!   If optcalc==1 (optional argument, default value is 0), then rhor will actually
!!   contains kinetic energy density (taur) instead of electronic density.
!!
!! NOTES
!! In the case of PAW calculations:
!!    All computations are done on the fine FFT grid.
!!    All variables (nfftf,ngfftf,mgfftf) refer to this fine FFT grid.
!!    All arrays (densities/potentials...) are computed on this fine FFT grid.
!!    Developers have to be careful when introducing others arrays:
!!      they have to be stored on the fine FFT grid.
!! In the case of norm-conserving calculations:
!!    The mesh is the usual augmented FFT grid to treat correctly the convolution.
!!
!! PARENTS
!!      bethe_salpeter,screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_mkrho(Wfd,Cryst,Psps,Kmesh,Bands,ngfftf,nfftf,rhor,&
                     optcalc) ! optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftf
 integer,intent(in),optional :: optcalc
 type(ebands_t),intent(in) :: Bands
 type(kmesh_t),intent(in) :: Kmesh
 type(crystal_t),intent(in) :: Cryst
 type(Pseudopotential_type),intent(in) :: Psps
 class(wfd_t),intent(inout) :: Wfd
!arrays
 integer,intent(in) :: ngfftf(18)
 real(dp),intent(out) :: rhor(nfftf,Wfd%nspden)

!Local variables ------------------------------
!scalars
 integer,parameter :: ndat1=1
 integer :: cplex,ib,ib_iter,ierr,ik,ir,is,n1,n2,n3,nfftotf
 integer :: alpha,nalpha,ipw,myoptcalc
 real(dp) :: kpt_cart,kg_k_cart,gp2pi1,gp2pi2,gp2pi3,cwftmp,bks_weight
 character(len=500) :: msg
 type(wave_t),pointer :: wave
!arrays
 integer,allocatable :: irrzon(:,:,:)
 real(dp),allocatable :: phnons(:,:,:),rhog(:,:),rhor_down(:),rhor_mx(:),rhor_my(:),cwavef(:,:)
 complex(dpc),allocatable :: wfr_x(:),wfr_y(:)
 complex(gwpc),allocatable :: gradug(:),work(:)
 complex(gwpc),allocatable,target :: wfr(:)
 complex(gwpc), ABI_CONTIGUOUS pointer :: cwavef1(:),cwavef2(:)
 type(iter2_t) :: Iter_bks

!*************************************************************************

 DBG_ENTER("COLL")

 ! Consistency check.
 ABI_CHECK(Wfd%nsppol == Bands%nsppol, "Mismatch in nsppol")

 if (ANY(ngfftf(1:3) /= Wfd%ngfft(1:3))) call wfd%change_ngfft(Cryst,Psps,ngfftf)

 ! Calculate IBZ contribution to the charge density.
 ABI_MALLOC(wfr, (nfftf*Wfd%nspinor))

 if (wfd%nspden == 4) then
   ABI_MALLOC(wfr_x, (nfftf))
   ABI_MALLOC(wfr_y, (nfftf))
   ABI_MALLOC(rhor_down, (nfftf))
   ABI_MALLOC(rhor_mx, (nfftf))
   ABI_MALLOC(rhor_my, (nfftf))
   rhor_down = zero; rhor_mx = zero; rhor_my = zero
 end if

 ! Update the (b,k,s) distribution table.
 call wfd%update_bkstab()

 ! Calculate the unsymmetrized density.
 rhor=zero
 myoptcalc=0; if (present(optcalc)) myoptcalc=optcalc
 nalpha=1; if (myoptcalc==1) nalpha=3
 if (myoptcalc == 1 .and. wfd%nspinor == 2) then
   MSG_ERROR("kinetic energy density with nspinor == 2 not implemented")
 end if

 ! Build the iterator that will distribute the states in an automated way.
 Iter_bks = wfd_iterator_bks(Wfd,bks_mask=ABS(Bands%occ)>=tol8)

 do alpha=1,nalpha
   do is=1,Wfd%nsppol
     do ik=1,Wfd%nkibz
       do ib_iter=1,iter_len(Iter_bks,ik,is)
         ib = iter_yield(Iter_bks,ib_iter,ik,is)
         bks_weight = Bands%occ(ib,ik,is) * Kmesh%wt(ik) / Cryst%ucvol

         call wfd%get_ur(ib,ik,is,wfr)

         cwavef1 => wfr(1:nfftf)
         if (myoptcalc == 1) then
           ABI_MALLOC(gradug,(Wfd%Kdata(ik)%npw))
           ABI_MALLOC(cwavef,(2,Wfd%Kdata(ik)%npw))
           ABI_MALLOC(work,(nfftf))

           ABI_CHECK(wfd%get_wave_ptr(ib, ik, is, wave, msg) == 0, msg)
           cwavef(1,:)= REAL(wave%ug(:))
           cwavef(2,:)=AIMAG(wave%ug(:))

           ! Multiplication by 2pi i (k+G)_alpha
           gp2pi1=Cryst%gprimd(alpha,1)*two_pi
           gp2pi2=Cryst%gprimd(alpha,2)*two_pi
           gp2pi3=Cryst%gprimd(alpha,3)*two_pi
           kpt_cart=gp2pi1*Wfd%kibz(1,ik)+gp2pi2*Wfd%kibz(2,ik)+gp2pi3*Wfd%kibz(3,ik)
           do ipw=1,Wfd%Kdata(ik)%npw
             kg_k_cart= gp2pi1*Wfd%Kdata(ik)%kg_k(1,ipw) + &
                        gp2pi2*Wfd%Kdata(ik)%kg_k(2,ipw) + &
                        gp2pi3*Wfd%Kdata(ik)%kg_k(3,ipw)+kpt_cart
!             ipwsp=ipw!+(ispinor-1)*Wfd%Kdata(ik)%npw
             cwftmp=-cwavef(2,ipw)*kg_k_cart
             cwavef(2,ipw)=cwavef(1,ipw)*kg_k_cart
             cwavef(1,ipw)=cwftmp
           end do
           gradug(:)=CMPLX(cwavef(1,:),cwavef(2,:),gwpc)
           call fft_ug(Wfd%npwarr(ik),nfftf,Wfd%nspinor,ndat1,Wfd%mgfft,Wfd%ngfft,&
             Wfd%istwfk(ik),Wfd%Kdata(ik)%kg_k,Wfd%Kdata(ik)%gbound,gradug,work)
           cwavef1(:)=work(:)
           ABI_FREE(work)
           ABI_FREE(cwavef)
           ABI_FREE(gradug)
         end if

!$OMP PARALLEL DO
         do ir=1,nfftf
           rhor(ir,is) = rhor(ir,is) + CONJG(cwavef1(ir)) * cwavef1(ir) * bks_weight
         end do
         !call cplx_addtorho(n1,n2,n3,n4,n5,n6,ndat,weight_r,ur,rho)

         if (wfd%nspinor == 2 .and. wfd%nspden == 1) then
           cwavef2 => wfr(1+nfftf:2*nfftf)
           do ir=1,nfftf
             rhor(ir, 1) = rhor(ir, 1) + CONJG(cwavef2(ir)) * cwavef2(ir) * bks_weight
           end do
         end if

         if (wfd%nspinor == 2 .and. wfd%nspden == 4) then
           cwavef2 => wfr(1+nfftf:2*nfftf)
           wfr_x(:) = cwavef1(:) + cwavef2(:)       ! $(\Psi^{1}+\Psi^{2})$
           wfr_y(:) = cwavef1(:) -j_dpc*cwavef2(:)  ! $(\Psi^{1}-i\Psi^{2})$
!$OMP PARALLEL DO
           do ir=1,nfftf
             rhor_down(ir) = rhor_down(ir) + CONJG(cwavef2(ir)) * cwavef2(ir) * bks_weight
             rhor_mx(ir) = rhor_mx(ir) + CONJG(wfr_x(ir)) * wfr_x(ir) * bks_weight
             rhor_my(ir) = rhor_my(ir) + CONJG(wfr_y(ir)) * wfr_y(ir) * bks_weight
           end do
         end if

       end do
     end do
   end do

 end do ! enddo alpha

 call iter_free(Iter_bks)

 select case (myoptcalc)
 case (0)
   ! density
   if (wfd%nspden == 4) then
     rhor(:, 2) = rhor_mx
     rhor(:, 3) = rhor_my
     rhor(:, 4) = rhor_down
   end if
 case (1)
   ! convention for taur = 1/2 Sum_i |grad phi_i|^2
   rhor(:,:)=half*rhor(:,:)

 case default
   MSG_ERROR(sjoin("Wrong myoptcalc:", itoa(myoptcalc)))
 end select

 call xmpi_sum(rhor,Wfd%comm,ierr)

 ! Symmetrization in G-space implementing also the AFM case
 n1=ngfftf(1); n2=ngfftf(2); n3=ngfftf(3); nfftotf=n1*n2*n3

 ABI_MALLOC(irrzon,(nfftotf**(1-1/Cryst%nsym),2,(Wfd%nspden/Wfd%nsppol)-3*(Wfd%nspden/4)))
 ABI_MALLOC(phnons,(2,nfftotf,(Wfd%nspden/Wfd%nsppol)-3*(Wfd%nspden/4)))

 if (Cryst%nsym/=1) then
   call irrzg(irrzon,Wfd%nspden,Wfd%nsppol,Cryst%nsym,n1,n2,n3,phnons,Cryst%symafm,Cryst%symrel,Cryst%tnons)
 end if

 ! Symmetrize rho(r), and pack nspden components following abinit conventions.
 cplex=1
 ABI_MALLOC(rhog,(2,cplex*nfftf))

 call symrhg(cplex,Cryst%gprimd,irrzon,Wfd%MPI_enreg,nfftf,nfftotf,ngfftf,Wfd%nspden,Wfd%nsppol,&
             Cryst%nsym,phnons,rhog,rhor,Cryst%rprimd,Cryst%symafm,Cryst%symrel,Cryst%tnons)

 ABI_FREE(rhog)
 ABI_FREE(phnons)
 ABI_FREE(irrzon)

 ! Find and print minimum and maximum total electron density
 ! (or total kinetic energy density, or total element of kinetic energy density tensor) and locations
 !call wrtout(std_out,'mkrho: echo density (plane-wave part only)','COLL')
 !call prtrhomxmn(std_out,wfd%mpi_enreg,nfftf,ngfftf,wfd%nspden,1,rhor,optrhor=optcalc,ucvol=crystl%ucvol)

 write(msg,'(a,f9.4)')' planewave contribution to nelect: ',SUM(rhor(:,1))*Cryst%ucvol/nfftf
 call wrtout(std_out, msg)

 if (Wfd%nspden==4) then
   write(msg,'(a,3f9.4)')&
     ' mx, my, mz: ',SUM(rhor(:,2))*Cryst%ucvol/nfftf,SUM(rhor(:,3))*Cryst%ucvol/nfftf,SUM(rhor(:,4))*Cryst%ucvol/nfftf
   call wrtout(std_out, msg)
 end if

 ABI_FREE(wfr)

 if (Wfd%nspden == 4) then
   ABI_FREE(wfr_x)
   ABI_FREE(wfr_y)
   ABI_FREE(rhor_down)
   ABI_FREE(rhor_mx)
   ABI_FREE(rhor_my)
 end if

 DBG_EXIT("COLL")

end subroutine wfd_mkrho
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/test_charge
!! NAME
!! test_charge
!!
!! FUNCTION
!!  Reports info on the electronic charge as well as Drude plasma frequency.
!!  Mainly used in the GW part.
!!
!! INPUTS
!!  nelectron_exp=Expected total number of electrons (used to normalize the charge)
!!
!! OUTPUT
!!
!! PARENTS
!!      bethe_salpeter,mrgscr,screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine test_charge(nfftf,nelectron_exp,nspden,rhor,ucvol,&
& usepaw,usexcnhat,usefinegrid,compch_sph,compch_fft,omegaplasma)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftf,nspden,usefinegrid,usepaw,usexcnhat
 real(dp),intent(in) :: compch_fft,compch_sph,ucvol,nelectron_exp
 real(dp),intent(out) :: omegaplasma
!arrays
 real(dp),intent(inout) :: rhor(nfftf,nspden)

!Local variables ------------------------------
!scalars
 real(dp) :: nelectron_tot,nelectron_fft
 real(dp) :: nelectron_pw,nelectron_sph,rhoav,rs,nratio
 character(len=500) :: msg

!*************************************************************************

! ABI_UNUSED(usexcnhat)
if (usexcnhat==0)then
end if

 ! === For PAW output of compensation charges ===
 if (usepaw==1) then
!if (usepaw==1.and.usexcnhat>0) then ! TODO I still dont understand this if!
   write(msg,'(4a)')ch10,' PAW TEST:',ch10,' ==== Compensation charge inside spheres ============'
   if (compch_sph<greatest_real.and.compch_fft<greatest_real) &
&    write(msg,'(3a)')TRIM(msg),ch10,' The following values must be close...'
   if (compch_sph<greatest_real) &
&    write(msg,'(3a,f22.15)')TRIM(msg),ch10,' Compensation charge over spherical meshes = ',compch_sph
   if (compch_fft<greatest_real) then
     if (usefinegrid==1) then
       write(msg,'(3a,f22.15)')TRIM(msg),ch10,' Compensation charge over fine fft grid    = ',compch_fft
     else
       write(msg,'(3a,f22.15)')TRIM(msg),ch10,' Compensation charge over fft grid         = ',compch_fft
     end if
   end if
   call wrtout([std_out, ab_out], msg)
   write(msg,'(a)')ch10
   call wrtout([std_out, ab_out], msg)
 end if !PAW

 nelectron_pw =SUM(rhor(:,1))*ucvol/nfftf
 nelectron_tot=nelectron_pw
 nratio       =nelectron_exp/nelectron_tot

 if (usepaw==1) then
   nelectron_sph=nelectron_pw+compch_sph
   nelectron_fft=nelectron_pw+compch_fft
   nelectron_tot=nelectron_sph
   nratio=(nelectron_exp-nelectron_sph)/nelectron_pw
 end if

 rhoav=nelectron_tot/ucvol ; rs=(three/(four_pi*rhoav))**third
 if (usepaw==0) then
  write(msg,'(2(a,f9.4))')&
   ' Number of electrons calculated from density = ',nelectron_tot,'; Expected = ',nelectron_exp
 else
   write(msg,'(2(a,f9.4),a)')&
   ' Total number of electrons per unit cell = ',nelectron_sph,' (Spherical mesh), ',nelectron_fft,' (FFT mesh)'
 end if
 call wrtout([std_out, ab_out], msg)

!$write(msg,'(a,f9.4)')' Renormalizing smooth charge density using nratio = ',nratio
!! rhor(:,:)=nratio*rhor(:,:)

 write(msg,'(a,f9.6)')' average of density, n = ',rhoav
 call wrtout([std_out, ab_out], msg)
 write(msg,'(a,f9.4)')' r_s = ',rs
 call wrtout([std_out, ab_out], msg)
 omegaplasma=SQRT(four_pi*rhoav)
 write(msg,'(a,f9.4,2a)')' omega_plasma = ',omegaplasma*Ha_eV,' [eV]',ch10
 call wrtout([std_out, ab_out], msg)

end subroutine test_charge
!!***

!!****f* m_wfd/wfd_pawrhoij
!! NAME
!! wfd_pawrhoij
!!
!! FUNCTION
!! Calculate the PAW quantities rhoij (augmentation occupancies)
!! Remember:for each atom, rho_ij=Sum_{n,k} {occ(n,k)*<Cnk|p_i><p_j|Cnk>}
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  cprj(natom,nspinor*mband*mkmem*nsppol)= wave functions projected with non-local projectors:
!!                                   cprj_nk(i)=<p_i|Cnk> where p_i is a non-local projector.
!!  istwfk(nkpt)=parameter that describes the storage of wfs
!!  kptopt=option for the generation of k points
!!  mband=maximum number of bands
!!  natom=number of atoms in cell
!!  nkpt=number of k points
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ(mband*nkpt*nsppol)=occupation number for each band for each k
!!  pawprtvol=control print volume and debugging output for PAW
!!
!! SIDE EFFECTS
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  On input: arrays dimensions
!!  On output:
!!    pawrhoij(:)%rhoij_(lmn2_size,nspden)=
!!          Sum_{n,k} {occ(n,k)*conjugate[cprj_nk(ii)].cprj_nk(jj)} (non symetrized)
!!
!! PARENTS
!!      paw_qpscgw
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_pawrhoij(Wfd,Cryst,Bst,kptopt,pawrhoij,pawprtvol)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: kptopt,pawprtvol
 type(crystal_t),intent(in) :: Cryst
 class(wfd_t),intent(inout) :: Wfd
 type(ebands_t),intent(in) :: Bst
!arrays
 type(pawrhoij_type),intent(inout) :: pawrhoij(Wfd%natom)

!Local variables ---------------------------------------
!scalars
 integer :: cplex,cplex_rhoij,qphase,iatom,band,ik_ibz
 integer :: spin,natinc,nband_k,option,lmn2_size,nspden
 logical :: usetimerev
 real(dp) :: occup,wtk_k
 character(len=500) :: msg
!arrays
 !real(dp) :: tsec(2)
 character(len=8),parameter :: dspin(6)=(/"up      ","down    ","dens (n)","magn (x)","magn (y)","magn (z)"/)
 type(pawcprj_type),allocatable :: cwaveprj(:,:)
 integer :: bks_distrb(Wfd%mband,Wfd%nkibz,Wfd%nsppol)
 integer :: got(Wfd%nproc)
 logical :: bks_mask(Wfd%mband,Wfd%nkibz,Wfd%nsppol)

!************************************************************************

 DBG_ENTER("COLL")

 ! Allocate temporary cwaveprj storage (sorted by atom type)
 ABI_DATATYPE_ALLOCATE(cwaveprj,(Wfd%natom,Wfd%nspinor))
 call pawcprj_alloc(cwaveprj,0,Wfd%nlmn_sort)

 ! Initialize output quantities if not already done.
 do iatom=1,Wfd%natom
   if (pawrhoij(iatom)%use_rhoij_==0) then
     cplex_rhoij= pawrhoij(iatom)%cplex_rhoij
     qphase     = pawrhoij(iatom)%qphase
     lmn2_size  = pawrhoij(iatom)%lmn2_size
     nspden     = pawrhoij(iatom)%nspden
     ABI_ALLOCATE(pawrhoij(iatom)%rhoij_,(cplex_rhoij*qphase*lmn2_size,nspden))
     pawrhoij(iatom)%use_rhoij_=1
   end if
   pawrhoij(iatom)%rhoij_=zero
 end do

 option=1; usetimerev=(kptopt>0.and.kptopt<3)

 ! Distribute (b,k,s).
 where (ABS(Bst%occ)>tol8)
   bks_mask=.TRUE.
 else where
   bks_mask=.FALSE.
 end where
 got = 0

 call wfd%bks_distrb(bks_distrb,got,bks_mask)

 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz

     nband_k=Wfd%nband(ik_ibz,spin)
     wtk_k=Bst%wtk(ik_ibz)

     cplex=2; if (Wfd%istwfk(ik_ibz)>1) cplex=1

     do band=1,nband_k

       if (bks_distrb(band,ik_ibz,spin) == Wfd%my_rank) then
         !locc_test = (abs(Bst%occ(band,ik_ibz,spin))>tol8)
         occup = Bst%occ(band,ik_ibz,spin)

          ! Extract cprj for current band cwaveprj are sorted by atom type.
          call wfd%get_cprj(band,ik_ibz,spin,Cryst,cwaveprj,sorted=.TRUE.)

          ! Accumulate contribution from (occupied) current band
          !if (locc_test) then
           call pawaccrhoij(Cryst%atindx,cplex,cwaveprj,cwaveprj ,0,spin,Wfd%natom,Wfd%natom,&
&            Wfd%nspinor,occup,option,pawrhoij,usetimerev,wtk_k)
          !end if
       end if
     end do !band

   end do !ik_ibz
 end do !spin
 !
 ! Free temporary cwaveprj storage.
 call pawcprj_free(cwaveprj)
 ABI_FREE(cwaveprj)
 !
 !==========================================
 ! MPI: need to exchange arrays between procs
 ! TODO it should be tested.
 call pawrhoij_mpisum_unpacked(pawrhoij,Wfd%comm)

 ! Print info.
 if (abs(pawprtvol)>=1) then
   natinc=1; if(Wfd%natom>1.and.pawprtvol>=0) natinc=Wfd%natom-1
   write(msg, '(7a)') ch10," PAW TEST:",ch10,' ========= Values of RHOIJ in wfd_pawrhoij =========',ch10
   call wrtout(std_out,msg,'COLL')
   do iatom=1,Cryst%natom,natinc
     call pawrhoij_print_rhoij(pawrhoij(iatom)%rhoij_,pawrhoij(iatom)%cplex_rhoij,&
                  pawrhoij(iatom)%qphase,iatom,Cryst%natom,&
                  unit=std_out,opt_prtvol=pawprtvol)
  end do
 end if

 DBG_EXIT("COLL")

end subroutine wfd_pawrhoij
!!***

end module m_wfd
!!***
