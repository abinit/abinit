!{\src2tex{textfont=tt}}
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
!! Copyright (C) 2004-2017 ABINIT group (MVer, MG)
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
 use m_profiling_abi
 use m_errors

 use m_kptrank,  only : kptrank_type, destroy_kptrank, copy_kptrank

 implicit none

 private 
!!***

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

   type(kptrank_type) :: kptrank_t            ! ranking of all kpoints on phonon calculation grid, and inverse rank

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

CONTAINS  !=========================================================================================================================
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
!!      elphon
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

subroutine elph_ds_clean(elph_ds)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elph_ds_clean'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(elph_type), intent(inout) :: elph_ds

! *************************************************************************

 !@elph_type
 if (allocated(elph_ds%qirredtofull))  then
   ABI_DEALLOCATE(elph_ds%qirredtofull)
 end if
 if (allocated(elph_ds%wtq))  then
   ABI_DEALLOCATE(elph_ds%wtq)
 end if
 if (allocated(elph_ds%n0))  then
   ABI_DEALLOCATE(elph_ds%n0)
 end if
 if (allocated(elph_ds%qpt_full))  then
   ABI_DEALLOCATE(elph_ds%qpt_full)
 end if
 if (allocated(elph_ds%gkk_intweight))  then
   ABI_DEALLOCATE(elph_ds%gkk_intweight)
 end if
 if (allocated(elph_ds%gkk_qpt))  then
   ABI_DEALLOCATE(elph_ds%gkk_qpt)
 end if
 if (allocated(elph_ds%gkk_rpt))  then
   ABI_DEALLOCATE(elph_ds%gkk_rpt)
 end if
 if (allocated(elph_ds%gkk2))  then
   ABI_DEALLOCATE(elph_ds%gkk2)
 end if
 if (allocated(elph_ds%gamma_qpt))  then
   ABI_DEALLOCATE(elph_ds%gamma_qpt)
 end if
 if (allocated(elph_ds%gamma_rpt))  then
   ABI_DEALLOCATE(elph_ds%gamma_rpt)
 end if
 if (allocated(elph_ds%phfrq))  then
   ABI_DEALLOCATE(elph_ds%phfrq)
 end if
 if (allocated(elph_ds%a2f))  then
   ABI_DEALLOCATE(elph_ds%a2f)
 end if
 if (allocated(elph_ds%qgrid_data))  then
   ABI_DEALLOCATE(elph_ds%qgrid_data)
 end if

 call elph_k_destroy (elph_ds%k_phon)
 call elph_k_destroy (elph_ds%k_fine)

 call destroy_kptrank (elph_ds%k_fine%kptrank_t)

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
!!      elphon
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

subroutine elph_tr_ds_clean(elph_tr_ds)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elph_tr_ds_clean'
!End of the abilint section

 implicit none 

!Arguments ------------------------------------
!scalars
 type(elph_tr_type), intent(inout) :: elph_tr_ds
 
! *************************************************************************

 !@elph_tr_type
 if (allocated(elph_tr_ds%el_veloc))  then
   ABI_DEALLOCATE(elph_tr_ds%el_veloc)
 end if
 if (allocated(elph_tr_ds%FSelecveloc_sq))  then
   ABI_DEALLOCATE(elph_tr_ds%FSelecveloc_sq)
 end if
 if (allocated(elph_tr_ds%veloc_sq0))  then
   ABI_DEALLOCATE(elph_tr_ds%veloc_sq0)
 end if
 if (allocated(elph_tr_ds%veloc_sq))  then
   ABI_DEALLOCATE(elph_tr_ds%veloc_sq)
 end if
 if (allocated(elph_tr_ds%dos_n0))  then
   ABI_DEALLOCATE(elph_tr_ds%dos_n0)
 end if
 if (allocated(elph_tr_ds%dos_n))  then
   ABI_DEALLOCATE(elph_tr_ds%dos_n)
 end if
 if (allocated(elph_tr_ds%en_all))  then
   ABI_DEALLOCATE(elph_tr_ds%en_all)
 end if
 if (allocated(elph_tr_ds%de_all))  then
   ABI_DEALLOCATE(elph_tr_ds%de_all)
 end if
 if (allocated(elph_tr_ds%gamma_qpt_tr))  then
   ABI_DEALLOCATE(elph_tr_ds%gamma_qpt_tr)
 end if
 if (allocated(elph_tr_ds%gamma_qpt_trin))  then
   ABI_DEALLOCATE(elph_tr_ds%gamma_qpt_trin)
 end if
 if (allocated(elph_tr_ds%gamma_qpt_trout))  then
   ABI_DEALLOCATE(elph_tr_ds%gamma_qpt_trout)
 end if
 if (allocated(elph_tr_ds%gamma_rpt_tr))  then
   ABI_DEALLOCATE(elph_tr_ds%gamma_rpt_tr)
 end if
 if (allocated(elph_tr_ds%gamma_rpt_trin))  then
   ABI_DEALLOCATE(elph_tr_ds%gamma_rpt_trin)
 end if
 if (allocated(elph_tr_ds%gamma_rpt_trout))  then
   ABI_DEALLOCATE(elph_tr_ds%gamma_rpt_trout)
 end if
 if (allocated(elph_tr_ds%a2f_1d_tr))  then
   ABI_DEALLOCATE(elph_tr_ds%a2f_1d_tr)
 end if
 if (allocated(elph_tr_ds%a2f_1d_trin))  then
   ABI_DEALLOCATE(elph_tr_ds%a2f_1d_trin)
 end if
 if (allocated(elph_tr_ds%a2f_1d_trout))  then
   ABI_DEALLOCATE(elph_tr_ds%a2f_1d_trout)
 end if
 if (allocated(elph_tr_ds%tmp_gkk_intweight))  then
   ABI_DEALLOCATE(elph_tr_ds%tmp_gkk_intweight)
 end if
 if (allocated(elph_tr_ds%tmp_gkk_intweight1))  then
   ABI_DEALLOCATE(elph_tr_ds%tmp_gkk_intweight1)
 end if
 if (allocated(elph_tr_ds%tmp_gkk_intweight2))  then
   ABI_DEALLOCATE(elph_tr_ds%tmp_gkk_intweight2)
 end if
 if (allocated(elph_tr_ds%tmp_velocwtk))  then
   ABI_DEALLOCATE(elph_tr_ds%tmp_velocwtk)
 end if
 if (allocated(elph_tr_ds%tmp_velocwtk1))  then
   ABI_DEALLOCATE(elph_tr_ds%tmp_velocwtk1)
 end if
 if (allocated(elph_tr_ds%tmp_velocwtk2))  then
   ABI_DEALLOCATE(elph_tr_ds%tmp_velocwtk2)
 end if
 if (allocated(elph_tr_ds%tmp_vvelocwtk))  then
   ABI_DEALLOCATE(elph_tr_ds%tmp_vvelocwtk)
 end if
 if (allocated(elph_tr_ds%tmp_vvelocwtk1))  then
   ABI_DEALLOCATE(elph_tr_ds%tmp_vvelocwtk1)
 end if
 if (allocated(elph_tr_ds%tmp_vvelocwtk2))  then
   ABI_DEALLOCATE(elph_tr_ds%tmp_vvelocwtk2)
 end if

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
!!
!! NOTES
!!
!! SOURCE

subroutine elph_k_copy(elph_k_in, elph_k_out)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elph_k_copy'
!End of the abilint section

 implicit none

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

 call copy_kptrank(elph_k_in%kptrank_t, elph_k_out%kptrank_t)

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
!!
!! NOTES
!!
!! SOURCE

subroutine elph_k_destroy(elph_k)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elph_k_destroy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(elph_kgrid_type), intent(inout) :: elph_k

! *************************************************************************

 !@elph_kgrid_type
 if(allocated(elph_k%irr2full)) then
   ABI_DEALLOCATE (elph_k%irr2full)
 end if
 if(allocated(elph_k%full2irr)) then
   ABI_DEALLOCATE (elph_k%full2irr)
 end if
 if(allocated(elph_k%full2full)) then
   ABI_DEALLOCATE (elph_k%full2full)
 end if
 if(allocated(elph_k%irredtoGS)) then
   ABI_DEALLOCATE (elph_k%irredtoGS)
 end if
 if(allocated(elph_k%new_irredtoGS)) then
   ABI_DEALLOCATE (elph_k%new_irredtoGS)
 end if
 if(allocated(elph_k%kpt)) then
   ABI_DEALLOCATE (elph_k%kpt)
 end if
 if(allocated(elph_k%kptirr)) then
   ABI_DEALLOCATE (elph_k%kptirr)
 end if
 if(allocated(elph_k%new_kptirr)) then
   ABI_DEALLOCATE (elph_k%new_kptirr)
 end if
 if(allocated(elph_k%my_kpt)) then
   ABI_DEALLOCATE (elph_k%my_kpt)
 end if
 if(allocated(elph_k%my_ikpt)) then
   ABI_DEALLOCATE (elph_k%my_ikpt)
 end if
 if(allocated(elph_k%wtk)) then
   ABI_DEALLOCATE (elph_k%wtk)
 end if
 if(allocated(elph_k%wtq)) then
   ABI_DEALLOCATE (elph_k%wtq)
 end if
 if(allocated(elph_k%wtkirr)) then
   ABI_DEALLOCATE (elph_k%wtkirr)
 end if
 if(allocated(elph_k%new_wtkirr)) then
   ABI_DEALLOCATE (elph_k%new_wtkirr)
 end if
 if(allocated(elph_k%velocwtk)) then
   ABI_DEALLOCATE (elph_k%velocwtk)
 end if
 if(allocated(elph_k%vvelocwtk)) then
   ABI_DEALLOCATE (elph_k%vvelocwtk)
 end if

 call destroy_kptrank (elph_k%kptrank_t)

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
!!      elphon
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

subroutine elph_k_procs(nproc, elph_k)

 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elph_k_procs'
!End of the abilint section

 implicit none

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

end module defs_elphon
!!***

