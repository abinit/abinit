!!****m* ABINIT/defs_elphon
!!
!! NAME
!! defs_elphon
!!
!! FUNCTION
!! This module contains the datastructures for elphon
!!  the different (huge) matrices will either be allocated and
!!  used, or be written to disk. All combinations should be feasible.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2020 ABINIT group (MVer, MG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!!  Contains the following datastructures:
!!   1) elph_type contains data and dimensions for the kpoints near the
!!      fermi surface and the $g_{k k+q}$ matrix elements
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

module defs_elphon

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_krank
 use m_crystal

 implicit none

 private
!!***

 public :: gam_mult_displ
 public :: complete_gamma
 public :: complete_gamma_tr
 public :: mkqptequiv

!----------------------------------------------------------------------
!!****t* defs_elphon/elph_kgrid_type
!! NAME
!! elph_kgrid_type
!!
!! FUNCTION
!! elph_kgrid_type contains k-point grid data and dimensions
!!  this is a sub object of elph_type
!!
!! SOURCE

  type,public :: elph_kgrid_type

   integer :: nband                           ! number of bands for weights
   integer :: nsppol                          ! number of spin pol for weights
   integer :: nsym                            ! number of symmetry operations
   integer :: nkpt                            ! number of k-points in full grid
   integer :: nkptirr                         ! number of k-points in irreducible grid
   integer :: new_nkptirr                     ! number of k-points in irreducible grid
   integer :: my_nkpt                         ! number of k-points on present processor

   type(krank_t) :: krank            ! ranking of all kpoints on phonon calculation grid, and inverse rank

   integer, allocatable :: irr2full(:)            ! correspondence of irred kpoints to a full one
   integer, allocatable :: full2irr(:,:)          ! correspondence of full k to one irred kpoints through sym and timrev
   integer, allocatable :: full2full(:,:,:)       ! symmetry mapping of kpoints
   integer, allocatable :: my_kpt(:)              ! flag for k-points belonging to present proc (= me index of proc for each k-point)
   integer, allocatable :: my_ikpt(:)             ! flag for k-points belonging to present proc (= me index of proc for each k-point)
   integer, allocatable :: irredtoGS(:)           ! (nkptirr)
   integer, allocatable :: new_irredtoGS(:)       ! (new_nkptirr)

   real(dp),allocatable :: kpt(:,:)               ! coordinates of the full kpoints from phonon calculation
   real(dp),allocatable :: kptirr(:,:)            ! irreducible k-points, for preliminary set up
   real(dp),allocatable :: new_kptirr(:,:)        ! irreducible k-points, for preliminary set up
   real(dp),allocatable :: wtk(:,:,:)             ! integration weights (see also gkk_intweight)
   real(dp),allocatable :: wtq(:,:,:)             ! integration weights (see also gkk_intweight)
   real(dp),allocatable :: wtkirr(:)              ! weights for irreducible kpoints, to sum over _whole_ BZ (not just Fermi Surface)
   real(dp),allocatable :: new_wtkirr(:)          ! weights for irreducible kpoints, to sum over _whole_ BZ (not just Fermi Surface)
   real(dp),allocatable :: velocwtk(:,:,:,:)      ! (nFSband,nkpt_fine,3,nsppol), v(k)*wtk
   real(dp),allocatable :: vvelocwtk(:,:,:,:,:)   ! (nFSband,nkpt_fine,3,3,nsppol), v(k)*v(k)*wtk

  end type elph_kgrid_type

  public :: elph_k_copy
  public :: elph_k_procs
  public :: elph_k_destroy
!!***

!----------------------------------------------------------------------

!!****t* defs_elphon/elph_type
!! NAME
!! elph_type
!!
!! FUNCTION
!! elph_type contains data and dimensions for the kpoints near the
!! fermi surface and the $g_{k k+q}$ matrix elements
!!
!! SOURCE

  type,public :: elph_type

   type(elph_kgrid_type) :: k_phon               ! object for k-grid of phonon calculation
   type(elph_kgrid_type) :: k_fine               ! object for fine k-grid for FS integration

   integer :: natom,nbranch,nFSband,nband
   integer :: minFSband,maxFSband                !Index of lower and upper bands used for the FS integration

   integer :: ngkkband                           !Number of bands kept in final gkk matrix elements:
                                                 !either 1 if sum is performed immediately
                                                 !or = nband if all elements are kept based on flag ep_keepbands

   integer :: nenergy                            ! number of points in energy
                                                 !   space for the electron band energies
                                                 !   energies from
                                                 !   ef-nenergy*delta_e to
                                                 !   ef+nenergy*delta_e

   integer :: n_pair                             ! number of pairs considered

   integer :: nqpt_full                          !number of q in full BZ
   integer :: nqptirred                          !number of irred q-points


   integer :: unita2f,unit_gkk2,unit_gkk_rpt
   integer :: unitgkq                            !units for file output

   integer :: gkqwrite
   integer :: gkk2write
   integer :: gkk_rptwrite

   integer :: ep_scalprod                        !flag to perform the scalar product
   integer :: symgkq                             !flag to symmetrize gkq matrix elements
   integer :: ep_keepbands                       !flag to sum over bands or not
   integer :: ep_lova                            ! 1 for lova, and 0 for general
   integer :: tuniformgrid                       !flag to expect uniform grid of q or not
   integer :: prtbltztrp                         !flag to output BoltzTraP input files

   integer :: na2f                               !dimensions and increments for a2F function
   integer :: nsppol                             ! number of spin polarization channels
   integer :: nspinor                            ! number of spinorial components
   integer :: telphint                           ! flag for integration over the FS with 0=tetrahedra 1=gaussians
   integer :: ep_nspline                         ! scale factor for spline interpolation in RTA
   integer :: ntemper                            ! number of temperature points
   integer :: use_k_fine                         ! flag for using fine k-grids for eigenvalues and velocities. 0=no 1=yes
   integer :: ep_int_gkk                         ! flag for interpolate gkk(1) or gamma (0)
   integer :: ep_b_min                           ! first band taken into account in FS integration (if telphint==2)
   integer :: ep_b_max                           ! last band taken into account in FS integration (if telphint==2)
   integer :: kptrlatt(3,3)                      ! kpoint grid generating vectors, as in abinit
   integer :: kptrlatt_fine(3,3)                 ! kpoint grid generating vectors, for fine grid used in FS integration

   real(dp) :: delta_e                           ! step in electronic energies, around Fermi level
   real(dp) :: omega_min,omega_max
   real(dp) :: a2fsmear,domega
   real(dp) :: nelect                            ! number of electrons per unit cell, eventually with extra charges for carriers in semiconductors.
   real(dp) :: occ_factor                        ! normalization for integrals over FS, for num of spins, spinors, etc...

   real(dp) :: mustar                            ! mustar parameter
   real(dp) :: fermie                            ! Fermi energy (Ha), either comes from wfk file or from anaddb input file
   real(dp) :: elphsmear                         ! smearing width for gaussian integration or buffer in energy for
                                                 ! calculations with tetrahedra (telphint=0)
   real(dp) :: tempermin                         ! minimum temperature at which resistivity etc are calculated (in K)
   real(dp) :: temperinc                         ! interval temperature grid on which resistivity etc are calculated (in K)

   character(len=fnlen) :: elph_base_name        !base name for output files

   integer,allocatable :: qirredtofull(:)            !mapping between the qpoints found in the GGK file
                                                 !and the array of qpoints generated by the code

   real(dp),allocatable :: wtq(:)                    !weight for each qpoint in the full grid spqt
                                                 !if a point is not in the IBZ ==>  wtq=0
                                                 !MG we can also use indqpt

   real(dp),allocatable :: n0(:)                     !DOS at the Fermi level (states/Ha/spin)
   real(dp),allocatable :: qpt_full(:,:)             !special q points obtained by the Monkhorst & Pack method,
                                                 !in reduced coordinates


   real(dp),allocatable :: gkk_intweight(:,:,:)      ! (nFSband,nkpt_fine,nsppol)
                                                 !integration weights for gkk matrix elements on FS:
                                                 !if ep_keepbands == 0 all are 1
                                                 !if ep_keepbands == 1 then = to wtk_phon in elphon
                                                 !DOES NOT INCLUDE FACTOR OF 1/nkpt_phon

   real(dp),allocatable :: gkk_velocwtk(:,:,:)      ! (nFSband,nkpt_fine,nsppol)

   real(dp),allocatable :: gkk_vvelocwtk(:,:,:)      ! (nFSband,nkpt_fine,nsppol)

   real(dp),allocatable :: gkk_qpt(:,:,:,:,:,:)      ! (2, ngkkband*ngkkband, nbranch*nbranch, nkpt_phon, nsppol, nqptirred)
                                                 !Now gkq contains gkk2 matrices on basic qpts,
                                                 !summed over bands if ngkkband==1


   real(dp),allocatable :: gkk_rpt(:,:,:,:,:,:)      ! (2, ngkkband**2, nbranch**2, nkpt_phon, nsppol, nrpt)
                                                 !For the moment, gkk_rpt in memory is out of the question
   real(dp),allocatable :: gkk2(:,:,:,:,:,:)         ! (nbranch, ngkkband,ngkkband, nkpt_phon, nkpt_phon, nsppol)

   real(dp),allocatable :: gamma_qpt(:,:,:,:)        !gamma matrices integrated over kpoint coeff
                                                 !  and bands: still depends on qpt
                                                 ! dims= 2, nbranch**2, nsppol, nqpt
   real(dp),allocatable :: gamma_rpt(:,:,:,:)
                                                 ! dims= 2, nbranch**2, nsppol, nrpt
!NOTE: choice to put nsppol before or after nqpt is a bit arbitrary
!   abinit uses nband,nkpt,nsppol, but here for convenience nkpt_phon,nsppol,nqpt
!   as interpolation is on qpt

   real(dp),allocatable :: phfrq(:,:)                !phonon frequencies
   real(dp),allocatable :: a2f(:,:,:)                !a2f function

   real(dp),allocatable :: qgrid_data(:,:,:,:)       !e-ph values calculated over the irreducible part of the q-grid:
                                                 !first entry  =  index of the q-point,
                                                 !second index =  branch index
                                                 !the third slice contains the frequency, the linewidth and lambda(q,nu)
                                                 !for that particular phonon mode
                                                 ! dims= nqptirred,elph_ds%nbranch,nsppol,3

 end type elph_type

 public :: elph_ds_clean
!!***

!----------------------------------------------------------------------

!!****t* defs_elphon/elph_tr_type
!! NAME
!! elph_tr_type
!!
!! FUNCTION
!! elph_tr_ds contains the necessary data for the transport properties
!!
!! SOURCE

  type,public :: elph_tr_type

     integer :: ifltransport
     integer :: unitgkq_trin,unitgkq_trout
     integer :: gkqwrite,gkqexist
     integer :: onegkksize

     character(len=fnlen) :: ddkfilename

     real(dp),allocatable :: dos_n0(:,:)                  ! (nT,nsppol) DOS at the Fermi level (states/Ha/spin) at input temperatures
     real(dp),allocatable :: dos_n(:,:)                   ! (nE,nsppol) DOS at the selected energies (states/Ha/spin)
     real(dp),allocatable :: en_all(:,:)                  ! (nE,nsppol) selected energies
     real(dp),allocatable :: de_all(:,:)                  ! (nE,nsppol) differences between selected energies
     real(dp),allocatable :: veloc_sq0(:,:,:)             ! (3,nsppol,nT)
     real(dp),allocatable :: veloc_sq(:,:,:)              ! (3,nsppol,nE)

     real(dp),allocatable :: el_veloc(:,:,:,:)        ! nkpt nband 3 nsppol
! the 9 = 3x3 is for the full tensorial transport coefficients
     real(dp),allocatable :: gamma_qpt_tr(:,:,:,:,:)    ! 2 9 branches**2 nsppol qpt
     real(dp),allocatable :: gamma_qpt_trin(:,:,:,:,:)  !idem
     real(dp),allocatable :: gamma_qpt_trout(:,:,:,:,:) !idem

     real(dp),allocatable :: gamma_rpt_tr(:,:,:,:,:,:,:)    !idem
     real(dp),allocatable :: gamma_rpt_trin(:,:,:,:,:)  !idem
     real(dp),allocatable :: gamma_rpt_trout(:,:,:,:,:) !idem

     real(dp),allocatable :: a2f_1d_tr(:,:,:,:,:,:)           ! nfreq 9 nsppol 4 n_pair ntemp
     real(dp),allocatable :: a2f_1d_trin(:,:,:)
     real(dp),allocatable :: a2f_1d_trout(:,:,:)

     real(dp),allocatable :: FSelecveloc_sq(:,:)       ! 3 nsppol

     real(dp),allocatable :: tmp_gkk_intweight(:,:,:,:)
     real(dp),allocatable :: tmp_gkk_intweight1(:,:,:)
     real(dp),allocatable :: tmp_gkk_intweight2(:,:,:)

     real(dp),allocatable :: tmp_velocwtk(:,:,:,:,:)
     real(dp),allocatable :: tmp_velocwtk1(:,:,:,:)
     real(dp),allocatable :: tmp_velocwtk2(:,:,:,:)

     real(dp),allocatable :: tmp_vvelocwtk(:,:,:,:,:,:)
     real(dp),allocatable :: tmp_vvelocwtk1(:,:,:,:,:)
     real(dp),allocatable :: tmp_vvelocwtk2(:,:,:,:,:)

  end type elph_tr_type

 public :: elph_tr_ds_clean
!!***

!----------------------------------------------------------------------

CONTAINS
!!***

!----------------------------------------------------------------------

!!****f* defs_elphon/elph_ds_clean
!!
!! NAME
!!   elph_ds_clean
!!
!! FUNCTION
!!   deallocate remaining arrays in the elph_ds datastructure
!!
!! INPUTS
!!  elph_ds = elphon datastructure
!!
!! PARENTS
!!      m_elphon
!!
!! CHILDREN
!!      krank%free,wrtout
!!
!! NOTES
!!
!! SOURCE

subroutine elph_ds_clean(elph_ds)

!Arguments ------------------------------------
!scalars
 type(elph_type), intent(inout) :: elph_ds

! *************************************************************************

 !@elph_type
 ABI_SFREE(elph_ds%qirredtofull)
 ABI_SFREE(elph_ds%wtq)
 ABI_SFREE(elph_ds%n0)
 ABI_SFREE(elph_ds%qpt_full)
 ABI_SFREE(elph_ds%gkk_intweight)
 ABI_SFREE(elph_ds%gkk_qpt)
 ABI_SFREE(elph_ds%gkk_rpt)
 ABI_SFREE(elph_ds%gkk2)
 ABI_SFREE(elph_ds%gamma_qpt)
 ABI_SFREE(elph_ds%gamma_rpt)
 ABI_SFREE(elph_ds%phfrq)
 ABI_SFREE(elph_ds%a2f)
 ABI_SFREE(elph_ds%qgrid_data)

 call elph_k_destroy (elph_ds%k_phon)
 call elph_k_destroy (elph_ds%k_fine)

 call elph_ds%k_fine%krank%free()

end subroutine elph_ds_clean
!!***

!----------------------------------------------------------------------

!!****f* defs_elphon/elph_tr_ds_clean
!!
!! NAME
!!   elph_tr_ds_clean
!!
!! FUNCTION
!!   deallocate remaining arrays in the elph_tr_ds datastructure
!!
!! INPUTS
!!  elph_tr_ds = elphon transport datastructure
!!
!! PARENTS
!!      m_elphon
!!
!! CHILDREN
!!      krank%free,wrtout
!!
!! NOTES
!!
!! SOURCE

subroutine elph_tr_ds_clean(elph_tr_ds)

!Arguments ------------------------------------
!scalars
 type(elph_tr_type), intent(inout) :: elph_tr_ds

! *************************************************************************

 !@elph_tr_type
 ABI_SFREE(elph_tr_ds%el_veloc)
 ABI_SFREE(elph_tr_ds%FSelecveloc_sq)
 ABI_SFREE(elph_tr_ds%veloc_sq0)
 ABI_SFREE(elph_tr_ds%veloc_sq)
 ABI_SFREE(elph_tr_ds%dos_n0)
 ABI_SFREE(elph_tr_ds%dos_n)
 ABI_SFREE(elph_tr_ds%en_all)
 ABI_SFREE(elph_tr_ds%de_all)
 ABI_SFREE(elph_tr_ds%gamma_qpt_tr)
 ABI_SFREE(elph_tr_ds%gamma_qpt_trin)
 ABI_SFREE(elph_tr_ds%gamma_qpt_trout)
 ABI_SFREE(elph_tr_ds%gamma_rpt_tr)
 ABI_SFREE(elph_tr_ds%gamma_rpt_trin)
 ABI_SFREE(elph_tr_ds%gamma_rpt_trout)
 ABI_SFREE(elph_tr_ds%a2f_1d_tr)
 ABI_SFREE(elph_tr_ds%a2f_1d_trin)
 ABI_SFREE(elph_tr_ds%a2f_1d_trout)
 ABI_SFREE(elph_tr_ds%tmp_gkk_intweight)
 ABI_SFREE(elph_tr_ds%tmp_gkk_intweight1)
 ABI_SFREE(elph_tr_ds%tmp_gkk_intweight2)
 ABI_SFREE(elph_tr_ds%tmp_velocwtk)
 ABI_SFREE(elph_tr_ds%tmp_velocwtk1)
 ABI_SFREE(elph_tr_ds%tmp_velocwtk2)
 ABI_SFREE(elph_tr_ds%tmp_vvelocwtk)
 ABI_SFREE(elph_tr_ds%tmp_vvelocwtk1)
 ABI_SFREE(elph_tr_ds%tmp_vvelocwtk2)

end subroutine elph_tr_ds_clean
!!***

!----------------------------------------------------------------------

!!****f* defs_elphon/elph_k_copy
!!
!! NAME
!!   elph_k_copy
!!
!! FUNCTION
!!   allocate and copy arrays in the elph_k datastructure
!!
!! INPUTS
!!  elph_k = elphon k-points datastructure
!!
!! PARENTS
!!
!! CHILDREN
!!      krank%free,wrtout
!!
!! NOTES
!!
!! SOURCE

subroutine elph_k_copy(elph_k_in, elph_k_out)

!Arguments ------------------------------------
!scalars
 type(elph_kgrid_type), intent(in) :: elph_k_in
 type(elph_kgrid_type), intent(out) :: elph_k_out

! *************************************************************************

 !@elph_kgrid_type
 elph_k_out%nband = elph_k_in%nband
 elph_k_out%nsppol = elph_k_in%nsppol
 elph_k_out%nsym = elph_k_in%nsym

 elph_k_out%nkpt = elph_k_in%nkpt
 elph_k_out%nkptirr = elph_k_in%nkptirr

 elph_k_out%my_nkpt = elph_k_in%my_nkpt

 ABI_ALLOCATE(elph_k_out%my_kpt,(elph_k_out%nkpt))
 elph_k_out%my_kpt = elph_k_in%my_kpt

 ABI_ALLOCATE(elph_k_out%my_ikpt,(elph_k_out%my_nkpt))
 elph_k_out%my_ikpt = elph_k_in%my_ikpt

 ABI_ALLOCATE(elph_k_out%kptirr,(3,elph_k_out%nkptirr))
 elph_k_out%kptirr = elph_k_in%kptirr
 ABI_ALLOCATE(elph_k_out%wtkirr,(elph_k_out%nkptirr))
 elph_k_out%wtkirr = elph_k_in%wtkirr

 ABI_ALLOCATE(elph_k_out%wtk,(elph_k_out%nband,elph_k_out%nkpt,elph_k_out%nsppol))
 elph_k_out%wtk = elph_k_in%wtk
 ABI_ALLOCATE(elph_k_out%kpt,(3,elph_k_out%nkpt))
 elph_k_out%kpt = elph_k_in%kpt

 elph_k_out%krank = elph_k_in%krank%copy()

 ABI_ALLOCATE(elph_k_out%irr2full,(elph_k_out%nkptirr))
 elph_k_out%irr2full = elph_k_in%irr2full
 ABI_ALLOCATE(elph_k_out%full2irr,(3,elph_k_out%nkpt))
 elph_k_out%full2irr = elph_k_in%full2irr
 ABI_ALLOCATE(elph_k_out%full2full,(2,elph_k_out%nsym,elph_k_out%nkpt))
 elph_k_out%full2full = elph_k_in%full2full

 ABI_ALLOCATE(elph_k_out%irredtoGS,(elph_k_out%nkptirr))
 elph_k_out%irredtoGS = elph_k_in%irredtoGS

end subroutine elph_k_copy
!!***

!----------------------------------------------------------------------

!!****f* defs_elphon/elph_k_destroy
!!
!! NAME
!!   elph_k_destroy
!!
!! FUNCTION
!!   deallocate arrays in the elph_k datastructure
!!
!! INPUTS
!!  elph_k = elphon k-points datastructure
!!
!! PARENTS
!!      defs_elphon
!!
!! CHILDREN
!!      krank%free,wrtout
!!
!! NOTES
!!
!! SOURCE

subroutine elph_k_destroy(elph_k)

!Arguments ------------------------------------
!scalars
 type(elph_kgrid_type), intent(inout) :: elph_k

! *************************************************************************

 !@elph_kgrid_type
 ABI_SFREE(elph_k%irr2full)
 ABI_SFREE(elph_k%full2irr)
 ABI_SFREE(elph_k%full2full)
 ABI_SFREE(elph_k%irredtoGS)
 ABI_SFREE(elph_k%new_irredtoGS)
 ABI_SFREE(elph_k%kpt)
 ABI_SFREE(elph_k%kptirr)
 ABI_SFREE(elph_k%new_kptirr)
 ABI_SFREE(elph_k%my_kpt)
 ABI_SFREE(elph_k%my_ikpt)
 ABI_SFREE(elph_k%wtk)
 ABI_SFREE(elph_k%wtq)
 ABI_SFREE(elph_k%wtkirr)
 ABI_SFREE(elph_k%new_wtkirr)
 ABI_SFREE(elph_k%velocwtk)
 ABI_SFREE(elph_k%vvelocwtk)

 call elph_k%krank%free()

end subroutine elph_k_destroy
!!***

!----------------------------------------------------------------------

!!****f* defs_elphon/elph_k_procs
!!
!! NAME
!!   elph_k_procs
!!
!! FUNCTION
!!   allocate kpt to processors, in the elph_k datastructure
!!
!! INPUTS
!!  nproc = number of k-parallel processors
!!  elph_k = elphon k-points datastructure
!!
!! PARENTS
!!      m_elphon
!!
!! CHILDREN
!!      krank%free,wrtout
!!
!! NOTES
!!
!! SOURCE

subroutine elph_k_procs(nproc, elph_k)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: nproc
 type(elph_kgrid_type), intent(inout) :: elph_k

 integer :: ikpt, me, ik_this_proc
! *************************************************************************

 if (allocated(elph_k%my_kpt)) then
   ABI_DEALLOCATE (elph_k%my_kpt)
 end if
 ABI_ALLOCATE (elph_k%my_kpt, (elph_k%nkpt))

 elph_k%my_kpt = 0
 elph_k%my_nkpt = 0
 me = xmpi_comm_rank(xmpi_world)
 do ikpt = 1, elph_k%nkpt
   elph_k%my_kpt(ikpt) = MOD(ikpt-1, nproc)
   if (elph_k%my_kpt(ikpt) == me) elph_k%my_nkpt = elph_k%my_nkpt + 1
 end do

! create inverse mapping from ik_this_proc to ikpt
 if (allocated(elph_k%my_ikpt)) then
   ABI_DEALLOCATE (elph_k%my_ikpt)
 end if
 ABI_ALLOCATE (elph_k%my_ikpt, (elph_k%my_nkpt))
 elph_k%my_ikpt = 0

 ik_this_proc = 0
 do ikpt = 1, elph_k%nkpt
   if (elph_k%my_kpt(ikpt) == me) then
     ik_this_proc = ik_this_proc + 1
     elph_k%my_ikpt(ik_this_proc) = ikpt
   end if
 end do
 ABI_CHECK(ik_this_proc == elph_k%my_nkpt, 'found inconsistent k distribution in processors')

 write (std_out,*) 'elph_k_procs : nkpt, distrib = ', elph_k%my_nkpt
 write (std_out,*) elph_k%my_kpt
 write (std_out,*) elph_k%my_ikpt

end subroutine elph_k_procs
!!***

!----------------------------------------------------------------------

!!****f* defs_elphon/gam_mult_displ
!!
!! NAME
!! gam_mult_displ
!!
!! FUNCTION
!! This routine takes the bare gamma matrices and multiplies them
!!  by the displ_red matrices (related to the scalprod variable)
!!
!! INPUTS
!!   nbranch = number of phonon branches (3*natom)
!!   displ_red = phonon mode displacement vectors in reduced coordinates.
!!   gam_bare = bare gamma matrices before multiplication
!!
!! OUTPUT
!!   gam_now = output gamma matrices multiplied by displacement matrices
!!
!! PARENTS
!!      m_a2ftr,m_elphon,m_iogkk
!!
!! CHILDREN
!!      krank%free,wrtout
!!
!! SOURCE

subroutine gam_mult_displ(nbranch, displ_red, gam_bare, gam_now)

!Arguments -------------------------------
 integer, intent(in)  :: nbranch
 real(dp), intent(in)  :: displ_red(2,nbranch,nbranch)
 real(dp), intent(in)  :: gam_bare(2,nbranch,nbranch)
 real(dp), intent(out) :: gam_now(2,nbranch,nbranch)

!Local variables -------------------------
 real(dp) :: zgemm_tmp_mat(2,nbranch,nbranch)

! *********************************************************************

 gam_now = zero

 call zgemm('c','n',nbranch,nbranch,nbranch,cone,displ_red,nbranch,gam_bare,nbranch,czero,zgemm_tmp_mat,nbranch)

 call zgemm('n','n',nbranch,nbranch,nbranch,cone,zgemm_tmp_mat,nbranch,displ_red,nbranch,czero,gam_now,nbranch)

end subroutine gam_mult_displ
!!***

!!****f* ABINIT/complete_gamma
!!
!! NAME
!! complete_gamma
!!
!! FUNCTION
!! Use the set of special q points calculated by the Monkhorst & Pack Technique.
!! Check if all the informations for the q points are present in the input gamma matrices.
!! Generate the gamma matrices (already summed over the FS) of the set of q points which
!! samples homogeneously the entire Brillouin zone.
!!
!! INPUTS
!! qpttoqpt = qpoint index mapping under symops
!!
!! OUTPUT
!! gamma_qpt = in/out: set of gamma matrix elements completed and symmetrized
!!    gamma_qpt(2,nbranch**2,nsppol,nqpt_full)
!!
!! PARENTS
!!      m_elphon,m_phgamma
!!
!! CHILDREN
!!      krank%free,wrtout
!!
!! SOURCE

subroutine complete_gamma(Cryst,nbranch,nsppol,nqptirred,nqpt_full,ep_scalprod,qirredtofull,qpttoqpt,gamma_qpt)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsppol,nbranch,nqptirred,nqpt_full,ep_scalprod
 type(crystal_t),intent(in) :: Cryst
!arrays
 integer,intent(in) :: qirredtofull(nqptirred)
 integer,intent(in) :: qpttoqpt(2,Cryst%nsym,nqpt_full)
 real(dp), intent(inout) :: gamma_qpt(2,nbranch**2,nsppol,nqpt_full)

!Local variables-------------------------------
!scalars
 integer :: ibranch,ieqqpt,ii,natom,nsym,iqpt,isppol,isym
 integer :: itim,jbranch,jj,kk,ll,neqqpt,iatom,ancestor_iatom,iqpt_fullbz
!arrays
 integer :: symmetrized_qpt(nqpt_full)
 integer :: gkk_flag(nbranch,nbranch,nsppol,nqpt_full)
 real(dp) :: ss(3,3)
 real(dp) :: tmp_mat(2,nbranch,nbranch)
 real(dp) :: tmp_mat2(2,nbranch,nbranch)
 real(dp) :: ss_allatoms(2,nbranch,nbranch)
 complex(dpc) :: c_one, c_zero
 real(dp),allocatable :: gkk_qpt_new(:,:,:),gkk_qpt_tmp(:,:,:)

! *********************************************************************

 c_one = dcmplx(one,zero)
 c_zero = dcmplx(zero,zero)

 natom = Cryst%natom
 nsym  = Cryst%nsym

!Generation of the gkk matrices relative to the q points
!of the set which samples the entire Brillouin zone

!set up flags for gamma_qpt matrices we have
 gkk_flag = -1
 do iqpt=1,nqptirred
   iqpt_fullbz = qirredtofull(iqpt)
   gkk_flag(:,:,:,iqpt_fullbz) = 1
 end do

 symmetrized_qpt(:) = -1

 ABI_ALLOCATE(gkk_qpt_new,(2,nbranch**2,nsppol))
 ABI_ALLOCATE(gkk_qpt_tmp,(2,nbranch**2,nsppol))

 do iqpt=1,nqpt_full
!
!  Already symmetrized?
   if (symmetrized_qpt(iqpt) == 1) cycle

   gkk_qpt_new(:,:,:) = zero

!  loop over qpoints equivalent to iqpt
   neqqpt=0
!  do not use time reversal symmetry to complete the qpoints:
!  do not know what happens to the gamma matrices
!  11/2011: MJV: time reversal is needed here if inversion is absent
!  - used in read_gkk and all reductions of q-points by symmetry.

   do itim=1,2
     do isym=1,nsym
!      ieqqpt is sent onto iqpt by itim/isym
       ieqqpt = qpttoqpt(itim,isym,iqpt)


       if (gkk_flag(1,1,1,ieqqpt) == -1) cycle
!      if we have information on this qpt
!      iqpt is equivalent to ieqqpt: get it from file or memory
       gkk_qpt_tmp(:,:,:) = gamma_qpt(:,:,:,ieqqpt)

       neqqpt=neqqpt+1

!
!      MJV note 02/2010:
!      the correspondence of symrel and symrec in the different cases, symmetrizing there
!      and back, has been fixed in the cases with and without scalprod (ie cartesian
!      and reduced real space coordinates) with respect to a calculation with no symmetries
!      I believe everything is settled, but still do not know why the 2 versions of the ss
!      matrices here use different rel/rec, instead of just being multiplied by the rprim gprim...
!
       if (ep_scalprod==1) then
         do ii=1,3
           do jj=1,3
             ss(ii,jj)=zero
             do kk=1,3
               do ll=1,3
                 ss(ii,jj)=ss(ii,jj)+Cryst%rprimd(ii,kk)*Cryst%symrel(kk,ll,isym)*Cryst%gprimd(ll,jj)
               end do
             end do
           end do
         end do
       else
         do ii=1,3
           do jj=1,3
             ss(ii,jj) = Cryst%symrec(ii,jj,isym)
           end do
         end do
       end if

       ss_allatoms(:,:,:) = zero
       do iatom=1,natom
         ancestor_iatom = Cryst%indsym(4,isym,iatom)
         ss_allatoms(1, (ancestor_iatom-1)*3+1:(ancestor_iatom-1)*3+3, (iatom-1)*3+1:(iatom-1)*3+3) = ss(1:3,1:3)
       end do


!      NOTE   ssinv(ii,jj)=ssinv(ii,jj)+Cryst%gprimd(ii,kk)*rprimd(jj,ll)*Cryst%symrec(ll,kk,isym)

       do isppol=1,nsppol
!        multiply by the ss matrices
         tmp_mat2(:,:,:) = zero
         tmp_mat(:,:,:) = reshape(gkk_qpt_tmp(:,:,isppol),(/2,nbranch,nbranch/))

         call ZGEMM ('N','N',nbranch,nbranch,nbranch,&
&         c_one,ss_allatoms,nbranch,tmp_mat,nbranch,c_zero,tmp_mat2,nbranch)

         call ZGEMM ('N','T',nbranch,nbranch,nbranch,&
&         c_one,tmp_mat2,nbranch,ss_allatoms,nbranch,c_zero,tmp_mat,nbranch)

!        add to gkk_qpt_new
         do ibranch =1,nbranch
           do jbranch =1,nbranch
             gkk_qpt_new(:,(jbranch-1)*nbranch+ibranch,isppol) = &
&             gkk_qpt_new(:,(jbranch-1)*nbranch+ibranch,isppol) + tmp_mat(:,jbranch,ibranch)
           end do
         end do
       end do ! isppol
!
     end do ! isym
   end do ! itim
!
   ABI_CHECK(neqqpt>0,'no q-points found equivalent to iqpt ')
!  Divide by number of equivalent qpts found.
   gkk_qpt_new(:,:,:) = gkk_qpt_new(:,:,:)/neqqpt

!  copy the symmetrized version into all the equivalent qpoints, appropriately transformed
   do itim=1,2
     do isym=1,nsym
!      ieqqpt is sent onto iqpt by itim/isym
       ieqqpt = qpttoqpt(itim,isym,iqpt)

       if (symmetrized_qpt(ieqqpt) /= -1) cycle
       gkk_qpt_tmp(:,:,:) = zero

!      use symrec matrices to get inverse transform from isym^{-1}
       if (ep_scalprod==1) then
         do ii=1,3
           do jj=1,3
             ss(ii,jj)=zero
             do kk=1,3
               do ll=1,3
!                Use inverse of symop matrix here to get back to ieqqpt (inv+transpose is in symrec and in gprimd)
                 ss(ii,jj)=ss(ii,jj)+Cryst%rprimd(ii,kk)*Cryst%symrec(ll,kk,isym)*Cryst%gprimd(ll,jj)
               end do
             end do
           end do
         end do
       else
         do ii=1,3
           do jj=1,3
             ss(ii,jj) = Cryst%symrel(jj,ii,isym)
           end do
         end do
       end if

       ss_allatoms(:,:,:) = zero
       do iatom=1,natom
         ancestor_iatom = Cryst%indsym(4,isym,iatom)
         ss_allatoms(1, (ancestor_iatom-1)*3+1:(ancestor_iatom-1)*3+3, (iatom-1)*3+1:(iatom-1)*3+3) = ss(1:3,1:3)
       end do

!      ! Use inverse of symop matrix here to get back to ieqqpt
!      ssinv(ii,jj)=ssinv(ii,jj)+gprimd(ii,kk)*rprimd(jj,ll)*Cryst%symrel(kk,ll,isym)

       do isppol=1,nsppol
!        multiply by the ss^{-1} matrices
         tmp_mat2(:,:,:) = zero
         tmp_mat(:,:,:) = reshape(gkk_qpt_new(:,:,isppol),(/2,nbranch,nbranch/))

         call ZGEMM ('N','N',nbranch,nbranch,nbranch,&
&         c_one,ss_allatoms,nbranch,tmp_mat,nbranch,c_zero,tmp_mat2,nbranch)

         call ZGEMM ('N','T',nbranch,nbranch,nbranch,&
&         c_one,tmp_mat2,nbranch,ss_allatoms,nbranch,c_zero,tmp_mat,nbranch)

!        FIXME: the following could just be a reshape
         do ibranch =1,nbranch
           do jbranch =1,nbranch
             gkk_qpt_tmp(:,(jbranch-1)*nbranch+ibranch,isppol) =&
&             tmp_mat(:,jbranch,ibranch)
           end do
         end do
         if (gkk_flag (1,1,isppol,ieqqpt) == -1) gkk_flag (:,:,isppol,ieqqpt) = 0
       end do ! end isppol do

!      save symmetrized matrices for qpt ieqqpt
       gamma_qpt(:,:,:,ieqqpt) = gkk_qpt_tmp(:,:,:)

       symmetrized_qpt(ieqqpt) = 1

     end do !isym
   end do !itim
 end do !iqpt

 ABI_DEALLOCATE(gkk_qpt_new)
 ABI_DEALLOCATE(gkk_qpt_tmp)

end subroutine complete_gamma
!!***

!!****f* ABINIT/complete_gamma_tr
!!
!! NAME
!! complete_gamma_tr
!!
!! FUNCTION
!! Use the set of special q points calculated by the Monkhorst & Pack Technique.
!! Check if all the informations for the q points are present in
!! the input gamma transport matrices.
!! Generate the gamma transport matrices (already summed over the FS) of the set of q points which
!! samples homogeneously the entire Brillouin zone.
!!
!! INPUTS
!! crystal<crystal_t>=data type gathering info on the crystalline structure.
!! ep_scalprod= flag for scalar product of gkk with phonon displacement vectors
!! nbranch=number of phonon branches = 3*natom
!! nqptirred=nqpt irred BZ
!! nqpt_full=nqpt full BZ
!! nsppol=number of spins
!! qirredtofull= mapping irred to full qpoints
!! qpttoqpt = qpoint index mapping under symops
!!
!! OUTPUT
!! gamma_qpt_tr = in/out: set of gamma matrix elements completed and symmetrized
!!    gamma_qpt_tr(2,9,nbranch*nbranch,nsppol,nqpt_full)
!!
!! PARENTS
!!      m_elphon
!!
!! CHILDREN
!!      krank%free,wrtout
!!
!! SOURCE

subroutine complete_gamma_tr(crystal,ep_scalprod,nbranch,nqptirred,nqpt_full,nsppol,gamma_qpt_tr,qirredtofull,qpttoqpt)

 use m_linalg_interfaces

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: nbranch,nqptirred,nqpt_full,nsppol, ep_scalprod
 type(crystal_t),intent(in) :: crystal
!arrays
 integer,intent(in) :: qpttoqpt(2,crystal%nsym,nqpt_full)
 integer,intent(in) :: qirredtofull(nqptirred)
 real(dp), intent(inout) :: gamma_qpt_tr(2,9,nbranch*nbranch,nsppol,nqpt_full)

!Local variables-------------------------------
!scalars
 integer :: ieqqpt,ii,iqpt,isppol,isym
 integer :: itim,jj,kk,ll,neqqpt
 integer :: iatom,ancestor_iatom
 integer :: iqpt_fullbz,imode, itensor
 real(dp),parameter :: tol=2.d-8
!arrays
 integer :: symrel(3,3,crystal%nsym),symrec(3,3,crystal%nsym)
 integer :: symmetrized_qpt(nqpt_full)
 integer :: gkk_flag(nbranch,nbranch,nsppol,nqpt_full)
 real(dp) :: gprimd(3,3),rprimd(3,3)
 real(dp) :: ss(3,3), sscart(3,3)
 real(dp) :: tmp_mat(2,nbranch,nbranch)
 real(dp) :: tmp_mat2(2,nbranch,nbranch)
 real(dp) :: tmp_tensor(2,3,3)
 real(dp) :: tmp_tensor2(2,3,3)
 real(dp) :: ss_allatoms(nbranch,nbranch)
 real(dp),allocatable :: gkk_qpt_new(:,:,:,:),gkk_qpt_tmp(:,:,:,:)

! *********************************************************************

 gprimd = crystal%gprimd
 rprimd = crystal%rprimd

 symrec =  crystal%symrec
 symrel =  crystal%symrel

!Generation of the gkk matrices relative to the q points
!of the set which samples the entire Brillouin zone

!set up flags for gamma_qpt matrices we have
 gkk_flag = -1
 do iqpt=1,nqptirred
   iqpt_fullbz = qirredtofull(iqpt)
   gkk_flag(:,:,:,iqpt_fullbz) = 1
 end do

 symmetrized_qpt(:) = -1
! isppol=1

 ABI_ALLOCATE(gkk_qpt_new,(2,9,nbranch*nbranch, nsppol))
 ABI_ALLOCATE(gkk_qpt_tmp,(2,9,nbranch*nbranch, nsppol))

 do iqpt=1,nqpt_full

!  Already symmetrized?
   if (symmetrized_qpt(iqpt) == 1) cycle

   gkk_qpt_new(:,:,:,:) = zero

!  loop over qpoints equivalent to iqpt
   neqqpt=0
!  do not use time reversal symmetry to complete the qpoints:
!  do not know what happens to the gamma matrices

   do itim=1,2
     do isym=1,crystal%nsym
!      ieqqpt is sent onto iqpt by itim/isym
       ieqqpt = qpttoqpt(itim,isym,iqpt)

       if (gkk_flag(1,1,1,ieqqpt) == -1) cycle
!      if we have information on this qpt
!      iqpt is equivalent to ieqqpt: get it from file or memory
       gkk_qpt_tmp(:,:,:,:) = gamma_qpt_tr(:,:,:,:,ieqqpt)

       neqqpt=neqqpt+1

!
!      MJV note 02/2010:
!      the correspondence of symrel and symrec in the different cases, symmetrizing there
!      and back, has been fixed in the cases with and without scalprod (ie cartesian
!      and reduced real space coordinates) with respect to a calculation with no symmetries
!      I believe everything is settled, but still do not know why the 2 versions of the ss
!      matrices here use different rel/rec, instead of just being multiplied by the rprim gprim...
!
       do ii=1,3
         do jj=1,3
           sscart(ii,jj)=0.0_dp
           do kk=1,3
             do ll=1,3
               sscart(ii,jj)=sscart(ii,jj)+rprimd(ii,kk)*symrel(kk,ll,isym)*gprimd(ll,jj)
!              sscart(ii,jj)=sscart(ii,jj)+rprimd(ii,kk)*symrel(kk,ll,isym)*gprimd(ll,jj)
             end do
           end do
         end do
       end do
       if (ep_scalprod==1) then
         do ii=1,3
           do jj=1,3
             ss(ii,jj)=0.0_dp
             do kk=1,3
               do ll=1,3
                 ss(ii,jj)=ss(ii,jj)+rprimd(ii,kk)*symrel(kk,ll,isym)*gprimd(ll,jj)
               end do
             end do
           end do
         end do
       else
         do ii=1,3
           do jj=1,3
             ss(ii,jj) = symrec(ii,jj,isym)
           end do
         end do
       end if

       ss_allatoms(:,:) = zero
       do iatom=1,crystal%natom
         ancestor_iatom = crystal%indsym(4,isym,iatom)
         ss_allatoms((ancestor_iatom-1)*3+1:(ancestor_iatom-1)*3+3,&
&         (iatom-1)*3+1:         (iatom-1)*3+3) = ss(1:3,1:3)
       end do


!      NOTE   ssinv(ii,jj)=ssinv(ii,jj)+gprimd(ii,kk)*rprimd(jj,ll)*symrec(ll,kk,isym)

       do isppol=1,nsppol

!        for each tensor component, rotate the cartesian directions of phonon modes
         do itensor = 1, 9
!          multiply by the ss matrices
           tmp_mat2(:,:,:) = zero
           tmp_mat(:,:,:) = reshape(gkk_qpt_tmp(:,itensor,:,isppol),&
&           (/2,nbranch,nbranch/))
           call DGEMM ('N','N',nbranch,nbranch,nbranch,&
&           one,ss_allatoms,nbranch,tmp_mat(1,:,:),nbranch,zero,&
&           tmp_mat2(1,:,:),nbranch)
           call DGEMM ('N','N',nbranch,nbranch,nbranch,&
&           one,ss_allatoms,nbranch,tmp_mat(2,:,:),nbranch,zero,&
&           tmp_mat2(2,:,:),nbranch)

           call DGEMM ('N','T',nbranch,nbranch,nbranch,&
&           one,tmp_mat2(1,:,:),nbranch,ss_allatoms,nbranch,zero,&
&           tmp_mat(1,:,:),nbranch)
           call DGEMM ('N','T',nbranch,nbranch,nbranch,&
&           one,tmp_mat2(2,:,:),nbranch,ss_allatoms,nbranch,zero,&
&           tmp_mat(2,:,:),nbranch)

           gkk_qpt_tmp(:,itensor,:,isppol) = reshape (tmp_mat, (/2,nbranch*nbranch/))
         end do ! itensor

!        for each cartesian direction/phonon mode, rotate the tensor components
         do imode = 1, nbranch*nbranch
           tmp_tensor2(:,:,:) = zero
           tmp_tensor(:,:,:) = reshape(gkk_qpt_tmp(:,:,imode,isppol),&
&           (/2,3,3/))
           call DGEMM ('N','N',3,3,3,&
&           one,sscart,3,tmp_tensor(1,:,:),3,zero,&
&           tmp_tensor2(1,:,:),3)
           call DGEMM ('N','T',3,3,3,&
&           one,tmp_tensor2(1,:,:),3,sscart,3,zero,&
&           tmp_tensor(1,:,:),3)

           call DGEMM ('N','N',3,3,3,&
&           one,sscart,3,tmp_tensor(2,:,:),3,zero,&
&           tmp_tensor2(2,:,:),3)
           call DGEMM ('N','T',3,3,3,&
&           one,tmp_tensor2(2,:,:),3,sscart,3,zero,&
&           tmp_tensor(2,:,:),3)

           gkk_qpt_tmp(:,:,imode,isppol) = reshape (tmp_tensor, (/2,9/)) ! modified by BX
         end do ! imode

!        add to gkk_qpt_new
         gkk_qpt_new(:,:,:,isppol) = gkk_qpt_new(:,:,:,isppol) + gkk_qpt_tmp(:,:,:,isppol)

       end do ! end isppol do

     end do ! end isym do
   end do ! end itim do

   ABI_CHECK(neqqpt>0,'no q-points found equivalent to iqpt ')

!  divide by number of equivalent qpts found
   gkk_qpt_new = gkk_qpt_new/neqqpt


!  copy the symmetrized version into all the equivalent qpoints, appropriately transformed
   do itim=1,2
     do isym=1,crystal%nsym
!      ieqqpt is sent onto iqpt by itim/isym
       ieqqpt = qpttoqpt(itim,isym,iqpt)

       if (symmetrized_qpt(ieqqpt) /= -1) cycle
       gkk_qpt_tmp = zero

!      use symrec matrices to get inverse transform from isym^{-1}
       do ii=1,3
         do jj=1,3
           sscart(ii,jj)=0.0_dp
           do kk=1,3
             do ll=1,3
!              Use inverse of symop matrix here to get back to ieqqpt (inv+transpose is in symrec and in gprimd)
!              sscart(ii,jj)=sscart(ii,jj)+rprimd(ii,kk)*symrec(ll,kk,isym)*gprimd(ll,jj)
               sscart(ii,jj)=sscart(ii,jj)+rprimd(ii,kk)*symrec(ll,kk,isym)*gprimd(ll,jj)
             end do
           end do
         end do
       end do
       if (ep_scalprod==1) then
         do ii=1,3
           do jj=1,3
             ss(ii,jj)=0.0_dp
             do kk=1,3
               do ll=1,3
!                Use inverse of symop matrix here to get back to ieqqpt (inv+transpose is in symrec and in gprimd)
                 ss(ii,jj)=ss(ii,jj)+rprimd(ii,kk)*symrec(ll,kk,isym)*gprimd(ll,jj)
               end do
             end do
           end do
         end do
       else
         do ii=1,3
           do jj=1,3
             ss(ii,jj) = symrel(jj,ii,isym)
           end do
         end do
       end if

       ss_allatoms(:,:) = zero
       do iatom=1,crystal%natom
         ancestor_iatom = crystal%indsym(4,isym,iatom)
         ss_allatoms((ancestor_iatom-1)*3+1:(ancestor_iatom-1)*3+3,&
&         (iatom-1)*3+1:          (iatom-1)*3+3) = ss(1:3,1:3)
       end do

!      ! Use inverse of symop matrix here to get back to ieqqpt
!      ssinv(ii,jj)=ssinv(ii,jj)+gprimd(ii,kk)*rprimd(jj,ll)*symrel(kk,ll,isym)

       do isppol=1,nsppol
         do itensor = 1, 9
!          multiply by the ss^{-1} matrices
           tmp_mat2(:,:,:) = zero
           tmp_mat(:,:,:) = reshape(gkk_qpt_new(:,itensor,:,isppol),&
&           (/2,nbranch,nbranch/))


           call DGEMM ('N','N',nbranch,nbranch,nbranch,&
&           one,ss_allatoms,nbranch,tmp_mat(1,:,:),nbranch,zero,&
&           tmp_mat2(1,:,:),nbranch)
           call DGEMM ('N','N',nbranch,nbranch,nbranch,&
&           one,ss_allatoms,nbranch,tmp_mat(2,:,:),nbranch,zero,&
&           tmp_mat2(2,:,:),nbranch)

           call DGEMM ('N','T',nbranch,nbranch,nbranch,&
&           one,tmp_mat2(1,:,:),nbranch,ss_allatoms,nbranch,zero,&
&           tmp_mat(1,:,:),nbranch)
           call DGEMM ('N','T',nbranch,nbranch,nbranch,&
&           one,tmp_mat2(2,:,:),nbranch,ss_allatoms,nbranch,zero,&
&           tmp_mat(2,:,:),nbranch)


           gkk_qpt_tmp(:,itensor,:,isppol) = reshape (tmp_mat, (/2,nbranch*nbranch/))
         end do ! itensor

!        for each cartesian direction/phonon mode, rotate the tensor components
         do imode = 1, nbranch*nbranch
           tmp_tensor2(:,:,:) = zero
           tmp_tensor(:,:,:) = reshape(gkk_qpt_tmp(:,:,imode,isppol),&
&           (/2,3,3/))
           call DGEMM ('N','N',3,3,3,&
&           one,sscart,3,tmp_tensor(1,:,:),3,zero,&
&           tmp_tensor2(1,:,:),3)
           call DGEMM ('N','T',3,3,3,&
&           one,tmp_tensor2(1,:,:),3,sscart,3,zero,&
&           tmp_tensor(1,:,:),3)

           call DGEMM ('N','N',3,3,3,&
&           one,sscart,3,tmp_tensor(2,:,:),3,zero,&
&           tmp_tensor2(2,:,:),3)
           call DGEMM ('N','T',3,3,3,&
&           one,tmp_tensor2(2,:,:),3,sscart,3,zero,&
&           tmp_tensor(2,:,:),3)

!          gkk_qpt_new(:,:,imode,isppol) = reshape (tmp_tensor, (/2,9/)) ! Modified by BX
           gkk_qpt_tmp(:,:,imode,isppol) = reshape (tmp_tensor, (/2,9/)) ! Modified by BX
         end do ! imode

         if (gkk_flag (1,1,isppol,ieqqpt) == -1) then
           gkk_flag (:,:,isppol,ieqqpt) = 0
         end if

       end do ! end isppol do


!      save symmetrized matrices for qpt ieqqpt
       gamma_qpt_tr(:,:,:,:,ieqqpt) = gkk_qpt_tmp(:,:,:,:)

       symmetrized_qpt(ieqqpt) = 1

     end do ! end isym do
   end do ! end itim do

 end do
!end iqpt do

 ABI_DEALLOCATE(gkk_qpt_new)
 ABI_DEALLOCATE(gkk_qpt_tmp)

end subroutine complete_gamma_tr
!!***

!----------------------------------------------------------------------

!!****f* m_fstab/mkqptequiv
!! NAME
!! mkqptequiv
!!
!! FUNCTION
!! This routine determines the equivalence between
!!   1) qpoints and fermi surface kpoints
!!   2) qpoints under symmetry operations
!!
!! INPUTS
!!   Cryst<crystal_t>=Info on unit cell and symmetries.
!!   kpt_phon = fermi surface kpoints
!!   nkpt_phon = number of kpoints in the full FS set
!!   nqpt = number of qpoints
!!   qpt_full = qpoint coordinates
!!
!! OUTPUT
!!   FSfullpqtofull = mapping of k + q onto k' for k and k' in full BZ
!!   qpttoqpt(itim,isym,iqpt) = qpoint index which transforms to iqpt under isym and with time reversal itim.
!!
!! NOTES
!!   REMOVED 3/6/2008: much too large matrix, and not used at present
!!       FStoqpt = mapping of kpoint pairs (1 irreducible and 1 full) to qpoints
!!
!! PARENTS
!!      m_a2ftr,m_elphon
!!
!! CHILDREN
!!      krank%free,wrtout
!!
!! SOURCE

subroutine mkqptequiv(FSfullpqtofull,Cryst,kpt_phon,nkpt_phon,nqpt,qpttoqpt,qpt_full,mqtofull)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt_phon,nqpt
 type(crystal_t),intent(in) :: Cryst
!arrays
 integer,intent(out) :: FSfullpqtofull(nkpt_phon,nqpt),qpttoqpt(2,Cryst%nsym,nqpt)
 integer,intent(out),optional :: mqtofull(nqpt)
 real(dp),intent(in) :: kpt_phon(3,nkpt_phon),qpt_full(3,nqpt)

!Local variables-------------------------------
!scalars
 integer :: ikpt_phon,iFSqpt,iqpt,isym,symrankkpt_phon
 !character(len=500) :: message
 type(krank_t) :: krank
!arrays
 real(dp) :: tmpkpt(3),gamma_kpt(3)

! *************************************************************************

 call wrtout(std_out,' mkqptequiv : making rankkpt_phon and invrankkpt_phon',"COLL")

 krank = krank_new(nkpt_phon, kpt_phon)

 FSfullpqtofull = -999
 gamma_kpt(:) = zero

 do ikpt_phon=1,nkpt_phon
   do iqpt=1,nqpt
     ! tmpkpt = jkpt = ikpt + qpt
     tmpkpt(:) = kpt_phon(:,ikpt_phon) + qpt_full(:,iqpt)

     ! which kpt is it among the full FS kpts?
     symrankkpt_phon = krank%get_rank(tmpkpt)

     FSfullpqtofull(ikpt_phon,iqpt) = krank%invrank(symrankkpt_phon)
     if (FSfullpqtofull(ikpt_phon, iqpt) == -1) then
       MSG_ERROR("looks like no kpoint equiv to k+q !!!")
     end if

   end do
 end do

 if (present(mqtofull)) then
   do iqpt=1,nqpt
     tmpkpt(:) = gamma_kpt(:) - qpt_full(:,iqpt)

     ! which kpt is it among the full FS kpts?
     symrankkpt_phon = krank%get_rank(tmpkpt)

     mqtofull(iqpt) = krank%invrank(symrankkpt_phon)
     if (mqtofull(iqpt) == -1) then
       MSG_ERROR("looks like no kpoint equiv to -q !!!")
     end if
   end do
 end if

 call krank%free()

 ! start over with q grid
 call wrtout(std_out,' mkqptequiv : FSfullpqtofull made. Do qpttoqpt',"COLL")

 krank = krank_new(nqpt, qpt_full)

 qpttoqpt(:,:,:) = -1
 do iFSqpt=1,nqpt
   do isym=1,Cryst%nsym
     tmpkpt(:) =  Cryst%symrec(:,1,isym)*qpt_full(1,iFSqpt) &
                + Cryst%symrec(:,2,isym)*qpt_full(2,iFSqpt) &
                + Cryst%symrec(:,3,isym)*qpt_full(3,iFSqpt)

     symrankkpt_phon = krank%get_rank(tmpkpt)
     if (krank%invrank(symrankkpt_phon) == -1) then
       MSG_ERROR("looks like no kpoint equiv to q by symmetry without time reversal!!!")
     end if
     qpttoqpt(1,isym,krank%invrank(symrankkpt_phon)) = iFSqpt

     tmpkpt = -tmpkpt
     symrankkpt_phon = krank%get_rank(tmpkpt)
     if (krank%invrank(symrankkpt_phon) == -1) then
       MSG_ERROR('looks like no kpoint equiv to q by symmetry with time reversal!!!')
     end if
     qpttoqpt(2,isym,krank%invrank(symrankkpt_phon)) = iFSqpt
   end do
 end do

 call krank%free()

end subroutine mkqptequiv
!!***

end module defs_elphon
!!***
