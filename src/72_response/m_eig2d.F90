!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_eig2d
!! NAME
!!  m_eig2d
!!
!! FUNCTION
!!  This module contains utilities to analyze and retrieve information
!!  from the second order derivative of the eigen-energies wrt
!!  displacements.
!!
!! COPYRIGHT
!! Copyright (C) 2014-2019 ABINIT group (SP, PB, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!  dfpt_looppert
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_eig2d

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_abicore
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_xmpi
 use m_ebands
 use m_cgtools

 use m_time,       only : timab
 use m_fstrings,   only : strcat
 use m_crystal,    only : crystal_init,  crystal_t
 use m_pawtab,     only : pawtab_type
 use m_ddb,        only : DDB_VERSION
 use m_ddb_hdr,    only : ddb_hdr_type, ddb_hdr_init, ddb_hdr_free, ddb_hdr_open_write
 use m_double_grid,only : kptfine_av
 use m_mpinfo,     only : distrb2, proc_distrb_cycle

 implicit none

 private

 public :: eigr2d_init              ! Main creation method of EIG2D.nc files.
 public :: eigr2d_ncwrite           ! Dump the object into NETCDF file.
 public :: eigr2d_free              ! Destruction method.
 public :: fan_init                 ! Main creation method of Fan.nc files.
 public :: fan_ncwrite              ! Dump the object into NETCDF file.
 public :: fan_free                 ! Destruction method.
 public :: gkk_init                 ! Main creation method of GKK.nc files.
 public :: gkk_ncwrite              ! Dump the object into NETCDF file.
 public :: gkk_free                 ! Destruction method.

 public :: eig2tot                  ! This routine calculates the second-order eigenvalues.
 public :: outbsd                   ! output bsd file for one perturbation (used for elphon calculations in anaddb)
 public :: eig2stern
 public :: elph2_fanddw             ! Calculates the zero-point motion corrections

!!***

!!****t* m_eig2d/eigr2d_t
!! NAME
!! eig2d_t
!!
!! FUNCTION
!! It contains information about the second-order derivative of the
!! eigenenergies wrt atomic displacement
!!
!! SOURCE

 type,public :: eigr2d_t

! WARNING : if you modify this datatype, please check whether there might be
! creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your
! modification.

  integer :: mband                 ! Max number of bands i.e MAXVAL(nband) (to dimension arrays)
  integer :: nsppol                ! number of spin-polarization
  integer :: nkpt                  ! number of k points
  integer :: natom                 ! number of atoms

  real(dp),allocatable :: eigr2d(:,:,:,:,:,:,:)
  ! eigr2d(2,mband*nsppol,nkpt,3,natom,3,natom)
  ! Second-order derivative of eigenergies (real,im) at each
  ! spin,band,k-point,dir1,dir2,natom1,natom2 .


 end type eigr2d_t
!!***

!!****t* m_eig2d/fan_t
!! NAME
!! fan_t
!!
!! FUNCTION
!! It contains information about the second-order derivative of the
!! eigenenergies wrt atomic displacement
!!
!! SOURCE

 type,public :: fan_t

! WARNING : if you modify this datatype, please check whether there might be
! creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your
! modification.

  integer :: mband                 ! Max number of bands i.e MAXVAL(nband) (to dimension arrays)
  integer :: nsppol                ! number of spin-polarization
  integer :: nkpt                  ! number of k points
  integer :: natom                 ! number of atoms

  real(dp),allocatable :: fan2d(:,:,:,:,:,:,:)
  ! fan2d(2*mband*nsppol,nkpt,3,natom,3,natom,mband)
  ! Second-order derivative of the eigenergies (real,im) at each
  ! ispin,iband(real,im),k-point,dir1,dir2,natom1,natom2,jband

 end type fan_t
!!***

!!****t* m_eig2d/gkk_t
!! NAME
!! gkk_t
!!
!! FUNCTION
!! It contains information about the second-order derivative of the
!! eigenenergies wrt atomic displacement
!!
!! SOURCE

 type,public :: gkk_t

! WARNING : if you modify this datatype, please check whether there might be
! creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your
! modification.

  integer :: mband                 ! Max number of bands i.e MAXVAL(nband) (to dimension arrays)
  integer :: nsppol                ! number of spin-polarization
  integer :: nkpt                  ! number of k points
  integer :: natom                 ! number of atoms
  integer :: ncart                 ! number of cartesian directions

  real(dp),allocatable :: gkk2d(:,:,:,:,:)
  ! gkk2d(2*mband*nsppol,nkpt,ncart,natom,mband)
  ! Second-order derivative of the eigenergies (real,im) at each
  ! ispin,iband(real,im),k-point,dir1,natom1,jband

 end type gkk_t
!!***

CONTAINS
!!***

!----------------------------------------------------------------------

!!****f* m_eig2d/eigr2d_init
!! NAME
!! eigr2d_init
!!
!! FUNCTION
!! This subroutine initializes the eigr2d_t structured datatype
!!
!! INPUTS
!! mbands=maximum number of bands
!! nkpt=number of k points
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! natom=number of atoms
!! eig2nkq(2,mband*nsppol,nkpt,3,natom,3,natom)=second-order derivative of the
!!     eigenenergies wrt phononic displacements
!!
!! OUTPUT
!! eigr2d<eigr2d_t>=the eigr2d_t datatype
!!
!! PARENTS
!!      dfpt_looppert,eig2tot
!!
!! CHILDREN
!!
!! SOURCE

subroutine eigr2d_init(eig2nkq,eigr2d,mband,nsppol,nkpt,natom)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) ::mband,nsppol,nkpt,natom
 type(eigr2d_t),intent(out) :: eigr2d
!arrays
 real(dp), intent(in) :: eig2nkq(2,mband*nsppol,nkpt,3,natom,3,natom)

! *************************************************************************

 eigr2d%mband = mband
 eigr2d%nsppol = nsppol
 eigr2d%nkpt = nkpt
 eigr2d%natom = natom

 ABI_MALLOC(eigr2d%eigr2d  ,(2,mband*nsppol,nkpt,3,natom,3,natom))
 eigr2d%eigr2d=eig2nkq

end subroutine eigr2d_init
!!***

!----------------------------------------------------------------------

!!****f* m_eig2d/eigr2d_ncwrite
!! NAME
!! eigr2d_ncwrite
!!
!! FUNCTION
!!  Writes the content of a eigr2d_t object to a NETCDF file
!!  according to the ETSF-IO specifications.
!!
!! INPUTS
!!  ncid =NC file handle
!!
!! OUTPUT
!!
!! PARENTS
!!      dfpt_looppert,eig2tot
!!
!! CHILDREN
!!
!! SOURCE

subroutine eigr2d_ncwrite(eigr2d,iqpt,wtq,ncid)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) ::ncid
 real(dp),intent(in) :: iqpt(3),wtq
 type(eigr2d_t),intent(in) :: eigr2d

!Local variables-------------------------------
#ifdef HAVE_NETCDF
 integer :: ncerr
 integer :: cplex,cart_dir,one_dim
 character(len=200) :: temp

! *************************************************************************

 ! ==============================================
 ! === Write the dimensions specified by ETSF ===
 ! ==============================================
 one_dim=1; cplex=2; cart_dir=3

 ncerr = nctk_def_dims(ncid, [&
   nctkdim_t('max_number_of_states', eigr2d%mband),&
   nctkdim_t('number_of_spins', eigr2d%nsppol),&
   nctkdim_t('number_of_kpoints', eigr2d%nkpt),&
   nctkdim_t('number_of_atoms', eigr2d%natom),&
   nctkdim_t('number_of_cartesian_directions', cart_dir),&
   nctkdim_t('current_one_dim', one_dim),&
   nctkdim_t('cplex', cplex),&
   nctkdim_t('product_mband_nsppol', eigr2d%mband*eigr2d%nsppol)], defmode=.True.)
 NCF_CHECK(ncerr)

 temp='cplex,product_mband_nsppol,number_of_kpoints,number_of_cartesian_directions,number_of_atoms,' //&
      'number_of_cartesian_directions , number_of_atoms'
 ncerr = nctk_def_arrays(ncid, [&
   nctkarr_t('current_q_point', "dp", 'number_of_cartesian_directions'), &
   nctkarr_t('current_q_point_weight', "dp", 'current_one_dim'), &
   nctkarr_t('second_derivative_eigenenergies', "dp", temp )])
!   nctkarr_t('second_derivative_eigenenergies', "dp",&
!   &'cplex, product_mband_nsppol, number_of_kpoints, number_of_cartesian_directions, number_of_atoms,&
!   &number_of_cartesian_directions, number_of_atoms')])
 NCF_CHECK(ncerr)

! Write data
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, vid('current_q_point'), iqpt))
 NCF_CHECK(nf90_put_var(ncid, vid('current_q_point_weight'), wtq))
 NCF_CHECK(nf90_put_var(ncid, vid('second_derivative_eigenenergies'), eigr2d%eigr2d))

#else
 MSG_ERROR("ETSF-IO support is not activated. ")
#endif


contains
 integer function vid(vname)
   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
 end function vid

end subroutine eigr2d_ncwrite
!!***

!----------------------------------------------------------------------

!!****f* m_eig2d/eigr2d_free
!! NAME
!! eigr2d_free
!!
!! FUNCTION
!! Deallocates the components of the eigr2d_t structured datatype
!!
!! INPUTS
!!  eigr2d<eigr2d_t>=The data type to be deallocated.
!!
!! OUTPUT
!!  Deallocate the dynamic arrays in the ebands_t type.
!!  (only deallocate)
!!
!! PARENTS
!!      dfpt_looppert,eig2tot
!!
!! CHILDREN
!!
!! SOURCE

subroutine eigr2d_free(eigr2d)

 implicit none

!Arguments ------------------------------------
!scalars
 type(eigr2d_t),intent(inout) :: eigr2d
! *************************************************************************
DBG_ENTER("COLL")

!Deallocate all components of bstruct

 if (allocated(eigr2d%eigr2d)) then
   ABI_FREE(eigr2d%eigr2d)
 end if

 DBG_EXIT("COLL")

end subroutine eigr2d_free
!!***

!!****f* m_eig2d/fan_init
!! NAME
!! fan_init
!!
!! FUNCTION
!! This subroutine initializes the fan_t structured datatype
!!
!! INPUTS
!! mbands=maximum number of bands
!! nkpt=number of k points
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! natom=number of atoms
!! fan2d(2*mband*nsppol,nkpt,3,natom,3,natom,mband*nsppol)=second-order derivative of the
!!     eigenenergies wrt phononic displacements
!!
!! OUTPUT
!! fan2d<fan_t>=the fan_t datatype
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      eig2tot
!!
!! CHILDREN
!!
!! SOURCE

subroutine fan_init(fan,fan2d,mband,nsppol,nkpt,natom)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) ::mband,nsppol,nkpt,natom
 type(fan_t),intent(out) :: fan2d
!arrays
 real(dp), intent(in) :: fan(2*mband*nsppol,nkpt,3,natom,3,natom,mband)
! *************************************************************************

 fan2d%mband = mband
 fan2d%nsppol = nsppol
 fan2d%nkpt = nkpt
 fan2d%natom = natom

 ABI_MALLOC(fan2d%fan2d,(2*mband*nsppol,nkpt,3,natom,3,natom,mband))
 fan2d%fan2d=fan

end subroutine fan_init
!!***

!!****f* m_eig2d/gkk_init
!! NAME
!! gkk_init
!!
!! FUNCTION
!! This subroutine initializes the gkk_t structured datatype
!!
!! INPUTS
!! mbands=maximum number of bands
!! nkpt=number of k points
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! natom=number of atoms
!! gkk2d(2*mband*nsppol,nkpt,3,natom,mband*nsppol)=second-order derivative of the
!!     eigenenergies wrt phononic displacements
!!
!! OUTPUT
!! gkk2d<gkk_t>=the gkk_t datatype
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      dfpt_looppert,eig2tot,m_gkk
!!
!! CHILDREN
!!
!! SOURCE

subroutine gkk_init(gkk,gkk2d,mband,nsppol,nkpt,natom,ncart)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) ::mband,nsppol,nkpt,natom,ncart
 type(gkk_t),intent(out) :: gkk2d
!arrays
 real(dp), intent(in) :: gkk(2*mband*nsppol,nkpt,ncart,natom,mband)
! *************************************************************************

 gkk2d%mband = mband
 gkk2d%nsppol = nsppol
 gkk2d%nkpt = nkpt
 gkk2d%natom = natom
 gkk2d%ncart = ncart

 ABI_MALLOC(gkk2d%gkk2d,(2*mband*nsppol,nkpt,ncart,natom,mband))
 gkk2d%gkk2d=gkk

end subroutine gkk_init
!!***

!----------------------------------------------------------------------

!!****f* m_eig2d/fan_ncwrite
!! NAME
!! fan_ncwrite
!!
!! FUNCTION
!!  Writes the content of a fan_t object to a NETCDF file
!!  according to the ETSF-IO specifications.
!!
!! INPUTS
!!  ncid =NC file handle
!!
!! OUTPUT
!!
!! PARENTS
!!      eig2tot
!!
!! CHILDREN
!!
!! SOURCE

subroutine fan_ncwrite(fan2d,iqpt,wtq,ncid)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) ::ncid
 real(dp),intent(in) :: iqpt(3),wtq
 type(fan_t),intent(in) :: fan2d

!Local variables-------------------------------
#ifdef HAVE_NETCDF
 integer :: ncerr
 integer :: cplex,cart_dir,one_dim
 character(len=200) :: temp


! *************************************************************************

 ! ==============================================
 ! === Write the dimensions specified by ETSF ===
 ! ==============================================
 one_dim=1; cplex=2; cart_dir=3

 ncerr = nctk_def_dims(ncid, [&
   nctkdim_t('max_number_of_states',fan2d%mband),&
   nctkdim_t('number_of_spins',fan2d%nsppol),&
   nctkdim_t('number_of_kpoints',fan2d%nkpt),&
   nctkdim_t('number_of_atoms',fan2d%natom),&
   nctkdim_t('3_number_of_atoms',3*fan2d%natom),&     ! TODO: not sure that variables can start with digits
   nctkdim_t('number_of_cartesian_directions',cart_dir),&
   nctkdim_t('current_one_dim',one_dim),&
   nctkdim_t('cplex',cplex),&
   nctkdim_t('product_mband_nsppol',fan2d%mband*fan2d%nsppol),&
   nctkdim_t('product_mband_nsppol2',fan2d%mband*fan2d%nsppol*2) &
 ], defmode=.True.)
 NCF_CHECK(ncerr)

 temp= 'product_mband_nsppol2, number_of_kpoints, number_of_cartesian_directions,' //&
   'number_of_atoms, number_of_cartesian_directions, number_of_atoms, max_number_of_states'
 ncerr = nctk_def_arrays(ncid, [&
   nctkarr_t('current_q_point', "dp", 'number_of_cartesian_directions'),&
   nctkarr_t('current_q_point_weight', "dp", 'current_one_dim'),&
   nctkarr_t('second_derivative_eigenenergies_actif', "dp", temp )])
!   nctkarr_t('second_derivative_eigenenergies_actif', "dp",&
!   &'product_mband_nsppol2, number_of_kpoints, number_of_cartesian_directions,&
!   &number_of_atoms, number_of_cartesian_directions, number_of_atoms, max_number_of_states')])
 NCF_CHECK(ncerr)

! Write data
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, vid('current_q_point'), iqpt))
 NCF_CHECK(nf90_put_var(ncid, vid('current_q_point_weight'), wtq))
 NCF_CHECK(nf90_put_var(ncid, vid('second_derivative_eigenenergies_actif'), fan2d%fan2d))

#else
 MSG_ERROR("netcdf support is not activated. ")
#endif

contains
 integer function vid(vname)
   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
 end function vid

end subroutine fan_ncwrite
!!***

!----------------------------------------------------------------------

!!****f* m_eig2d/gkk_ncwrite
!! NAME
!! gkk_ncwrite
!!
!! FUNCTION
!!  Writes the content of a gkk_t object to a NETCDF file
!!  according to the ETSF-IO specifications.
!!
!! INPUTS
!!  ncid =NC file handle
!!
!! OUTPUT
!!
!! PARENTS
!!      dfpt_looppert,eig2tot,m_gkk
!!
!! CHILDREN
!!
!! SOURCE

subroutine gkk_ncwrite(gkk2d,iqpt,wtq,ncid)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) ::ncid
 real(dp), intent(in) :: iqpt(3),wtq
 type(gkk_t),intent(in) :: gkk2d

!Local variables-------------------------------
#ifdef HAVE_NETCDF
 integer :: cplex,one_dim,ncerr,vid_

! *************************************************************************

 ! ==============================================
 ! === Write the dimensions specified by ETSF ===
 ! ==============================================
 one_dim=1; cplex=2

 ncerr = nctk_def_dims(ncid, [ &
&   nctkdim_t('max_number_of_states', gkk2d%mband), &
&   nctkdim_t('number_of_spins', gkk2d%nsppol), &
&   nctkdim_t('number_of_kpoints', gkk2d%nkpt), &
&   nctkdim_t('number_of_atoms_for_gkk', gkk2d%natom), &
&   nctkdim_t('3_number_of_atoms', 3*gkk2d%natom), &
&   nctkdim_t('number_of_cartesian_directions_for_gkk', gkk2d%ncart), &
&   nctkdim_t('current_one_dim', one_dim), &
&   nctkdim_t('cplex', cplex), &
&   nctkdim_t('product_mband_nsppol', gkk2d%mband*gkk2d%nsppol), &
&   nctkdim_t('product_mband_nsppol2', gkk2d%mband*gkk2d%nsppol*2) &
& ], defmode=.True.)
 NCF_CHECK(ncerr)

!arrays
 ncerr = nctk_def_arrays(ncid, [&
&   nctkarr_t('current_q_point', "dp", "number_of_cartesian_directions"), &
&   nctkarr_t('current_q_point_weight', "dp", 'current_one_dim'), &
&   nctkarr_t('second_derivative_eigenenergies_actif', "dp", &
&     'product_mband_nsppol2, number_of_kpoints, number_of_cartesian_directions_for_gkk,'// &
&     'number_of_atoms_for_gkk, max_number_of_states') &
& ])
 NCF_CHECK(ncerr)

 NCF_CHECK(nctk_set_datamode(ncid))
 vid_=vid('current_q_point')
 NCF_CHECK(nf90_put_var(ncid, vid_, iqpt))
 vid_=vid('current_q_point_weight')
 NCF_CHECK(nf90_put_var(ncid, vid_, wtq))
 vid_=vid('second_derivative_eigenenergies_actif')
 NCF_CHECK(nf90_put_var(ncid, vid_, gkk2d%gkk2d))

#else
 MSG_ERROR("netcdf support is not activated. ")
#endif

contains
 integer function vid(vname)
   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
 end function vid

end subroutine gkk_ncwrite
!!***

!----------------------------------------------------------------------

!!****f* m_eig2d/fan_free
!! NAME
!! fan_free
!!
!! FUNCTION
!! Deallocates the components of the fan_t structured datatype
!!
!! INPUTS
!!  fan2d<fan_t>=The data type to be deallocated.
!!
!! OUTPUT
!!  Deallocate the dynamic arrays in the fan_t type.
!!  (only deallocate)
!!
!! PARENTS
!!      eig2tot
!!
!! CHILDREN
!!
!! SOURCE

subroutine fan_free(fan2d)

 implicit none

!Arguments ------------------------------------
!scalars
 type(fan_t),intent(inout) :: fan2d
! *************************************************************************
DBG_ENTER("COLL")

!Deallocate all components of bstruct

 if (allocated(fan2d%fan2d)) then
   ABI_FREE(fan2d%fan2d)
 end if

 DBG_EXIT("COLL")

end subroutine fan_free
!!***

!----------------------------------------------------------------------

!!****f* m_eig2d/gkk_free
!! NAME
!! gkk_free
!!
!! FUNCTION
!! Deallocates the components of the gkk_t structured datatype
!!
!! INPUTS
!!  gkk2d<gkk_t>=The data type to be deallocated.
!!
!! OUTPUT
!!  Deallocate the dynamic arrays in the gkk_t type.
!!  (only deallocate)
!!
!! PARENTS
!!      dfpt_looppert,eig2tot,m_gkk
!!
!! CHILDREN
!!
!! SOURCE

subroutine gkk_free(gkk2d)

 implicit none

!Arguments ------------------------------------
!scalars
 type(gkk_t),intent(inout) :: gkk2d
! *************************************************************************
DBG_ENTER("COLL")

!Deallocate all components of bstruct

 if (allocated(gkk2d%gkk2d)) then
   ABI_FREE(gkk2d%gkk2d)
 end if

 DBG_EXIT("COLL")

end subroutine gkk_free
!!***

!!****f* ABINIT/eig2stern
!! NAME
!! eig2stern
!!
!! FUNCTION
!! This routine calculates the second-order eigenvalues.
!! The output eig2nkq is this quantity for the input k points.
!!
!! INPUTS
!!  bdeigrf = number of bands for which to calculate the second-order eigenvalues.
!!  clflg(3,mpert)= array on calculated perturbations for eig2rf.
!!  dim_eig2nkq = 1 if eig2nkq is to be computed.
!!  cg1_pert(2,mpw1*nspinor*mband*mk1mem*nsppol,3,mpert) = first-order wf in G
!!            space for each perturbation. The wavefunction is orthogonal to the
!!            active space.
!!  gh0c1_pert(2,mpw1*nspinor*mband*mk1mem*nsppol,3,mpert) = matrix containing the
!!            vector:  <G|H(0)|psi(1)>, for each perturbation.
!!  gh1c_pert(2,mpw1*nspinor*mband*mk1mem*nsppol,3,mpert)) = matrix containing the
!!            vector:  <G|H(1)|n,k>, for each perturbation. The wavefunction is
!!            orthogonal to the active space.
!!  eigbrd(2,mband*nsppol,nkpt,3,natom,3,natom) = broadening factors for the
!!            electronic eigenvalues (optional).
!!  eigen0(nkpt_rbz*mband*nsppol) = 0-order eigenvalues at all K-points:
!!            <k,n'|H(0)|k,n'> (hartree).
!!  eigenq(nkpt_rbz*mband*nsppol) = 0-order eigenvalues at all shifted K-points:
!!            <k+Q,n'|H(0)|k+Q,n'> (hartree).
!!  eigen1(nkpt_rbz*2*nsppol*mband**2,3,mpert) = matrix of first-order:
!!            <k+Q,n'|H(1)|k,n> (hartree) (calculated in dfpt_cgwf).
!!  eig2nkq(2,mband*nsppol,nkpt,3,natom,3,natom*dim_eig2nkq) = second derivatives of
!!            the electronic eigenvalues.
!!  elph2_imagden = imaginary part of the denominator of the sum-over-state expression
!!            for the electronic eigenenergy shift due to second-order electron-phonon
!!            interation.
!!  ieig2rf = integer for calculation type.
!!  indsym(4,nsym,natom) = indirect indexing array for atom labels
!!            (not used yet, but will be used with symmetries).
!!  istwfk_pert(nkpt_rbz,3,mpert) = integer for choice of storage of wavefunction at
!!            each k point for each perturbation.
!!  mband = maximum number of bands.
!!  mk1mem = maximum number of k points which can fit in memory (RF data);
!!            0 if use disk.
!!  mpert = maximum number of perturbations.
!!  natom = number of atoms in the unit cell.
!!  npert = number of phonon perturbations, without taking into account directions:
!!            natom.
!!  nsym = number of symmetries (not used yet).
!!  mpi_enreg = information about MPI parallelization.
!!  mpw1 = maximum number of planewaves used to represent first-order wavefunctions.
!!  nkpt_rbz = number of k-points for each perturbation.
!!  npwar1(nkpt_rbz,mpert) = number of planewaves at k-point for first-order.
!!  nspinor = number of spinorial components of the wavefunctions.
!!  nsppol = 1 for unpolarized, 2 for spin-polarized.
!!  occ(mband*nkpt*nsppol)=occup number for each band (often 2) at each k point
!!  smdelta = integer controling the calculation of electron lifetimes.
!!  symq(4,2,nsym) = 1 if symmetry preserves present qpoint. From littlegroup_q (not used yet).
!!  symrec(3,3,nsym) = 3x3 matrices of the group symmetries (reciprocal space)
!!            (not used yet).
!!  symrel(3,3,nsym) = array containing the symmetries in real space (not used yet).
!!  timrev = 1 if time-reversal preserves the q wavevector; 0 otherwise
!!            (not in use yet).
!!  dtset = OPTIONAL, dataset structure containing the input variable of the
!!            calculation. This is required to use the k-interpolation routine.
!!  eigenq_fine(mband_fine,mkpt_fine,nsppol_fine) = OPTIONAL, 0-order eigenvalues
!!            at all shifted K-points: <k+Q,n'|H(0)|k+Q,n'> (hartree) of the
!!            fine grid. This information is read from the WF dense k-grid file.
!!  hdr_fine = OPTIONAL, header of the WF file of the fine k-point grid. This
!!            variable is required for the k-interpolation routine.
!!  hdr0     = OPTIONAL, header of the GS WF file of the corse k-point grid. This
!!            variable is required for the k-interpolation routine.
!!
!! OUTPUT
!!  eig2nkq(2,mband*nsppol,nkpt_rbz,3,npert,3,npert)= diagonal part of the
!!            second-order eigenvalues: E^{(2),diag}_{k,q,j}.
!!  eigbrd(2,mband*nsppol,nkpt_rbz,3,npert,3,npert)= OPTIONAL, array containing the
!!            electron lifetimes.
!!
!! PARENTS
!!      dfpt_looppert
!!
!! CHILDREN
!!      distrb2,dotprod_g,kptfine_av,smeared_delta,timab,wrtout,xmpi_sum
!!
!! SOURCE

subroutine eig2stern(occ,bdeigrf,clflg,cg1_pert,dim_eig2nkq,dim_eig2rf,eigen0,eigenq,&
&  eigen1,eig2nkq,elph2_imagden,esmear,gh0c1_pert,gh1c_pert,ieig2rf,istwfk_pert,&
&  mband,mk1mem,mpert,npert,mpi_enreg,mpw1,nkpt_rbz,npwar1,nspinor,nsppol,smdelta,&
&  dtset,eigbrd,eigenq_fine,hdr_fine,hdr0)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bdeigrf,dim_eig2nkq,dim_eig2rf,ieig2rf,mband,mk1mem,mpert,mpw1,nkpt_rbz
 integer,intent(in) :: npert,nspinor,nsppol,smdelta
 real(dp),intent(in) :: elph2_imagden,esmear
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: clflg(3,mpert)
 integer,intent(in) :: istwfk_pert(nkpt_rbz,3,mpert)
 integer,intent(in) :: npwar1(nkpt_rbz,mpert)
 real(dp),intent(in) :: cg1_pert(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf,3,mpert)
 real(dp),intent(in) :: gh0c1_pert(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf,3,mpert)
 real(dp),intent(in) :: gh1c_pert(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf,3,mpert)
 real(dp),intent(inout) :: eigen0(nkpt_rbz*mband*nsppol)
 real(dp),intent(in) :: eigen1(nkpt_rbz*2*nsppol*mband**2,3,mpert)
 real(dp),intent(inout) :: eigenq(nkpt_rbz*mband*nsppol)
 real(dp),intent(out) :: eig2nkq(2,mband*nsppol,nkpt_rbz,3,npert,3,npert*dim_eig2nkq)
 real(dp),intent(out),optional :: eigbrd(2,mband*nsppol,nkpt_rbz,3,npert,3,npert)
 real(dp),intent(in),pointer,optional :: eigenq_fine(:,:,:)
 real(dp), intent(in) :: occ(mband*nkpt_rbz*nsppol)
 type(dataset_type), intent(in) :: dtset
 type(hdr_type),intent(in),optional :: hdr_fine,hdr0

!Local variables-------------------------------
!tolerance for non degenerated levels
!scalars
 integer :: band2tot_index,band_index,bandtot_index,iband,icg2,idir1,idir2
 integer :: ikpt,ipert1,ipert2,isppol,istwf_k,jband,npw1_k,nkpt_sub,ikpt2
!integer :: ipw
 integer :: master,me,spaceworld,ierr
!real(dp),parameter :: etol=1.0d-3
 real(dp),parameter :: etol=1.0d-6
!real(dp),parameter :: etol=zero
 real(dp) :: ar,ai,deltae,den,dot2i,dot2r,dot3i,dot3r,doti,dotr,eig1_i1,eig1_i2
 real(dp) :: eig1_r1,eig1_r2,eig2_diai,den_av
 real(dp) :: wgt_int
 real(dp) :: eig2_diar,eigbrd_i,eigbrd_r
 character(len=500) :: message
 character(len=500) :: msg
!DBSP
! character(len=300000) :: message2
!END
 logical :: test_do_band
!arrays
 integer, allocatable :: nband_rbz(:),icg2_rbz(:,:)
 integer,pointer      :: kpt_fine_sub(:)
 real(dp)             :: tsec(2)
 real(dp),allocatable :: cwavef(:,:),cwavef2(:,:),center(:),eigen0tmp(:),eigenqtmp(:)
 real(dp) :: eigen(mband*nsppol),eigen_prime(mband*nsppol)
 real(dp),allocatable :: gh(:,:),gh1(:,:),ghc(:,:)
 real(dp),allocatable :: smdfun(:,:)
 real(dp),pointer     :: wgt_sub(:)

! *********************************************************************

!Init parallelism
 master =0
 spaceworld=mpi_enreg%comm_cell
 me=mpi_enreg%me_kpt
!DEBUG
!write(std_out,*)' eig2stern : enter '
!write(std_out,*)' mpw1=',mpw1
!write(std_out,*)' mband=',mband
!write(std_out,*)' nsppol=',nsppol
!write(std_out,*)' nkpt_rbz=',nkpt_rbz
!write(std_out,*)' npert=',npert
!ENDDEBUG

!Init interpolation method
 if(present(eigenq_fine))then
   ABI_ALLOCATE(center,(3))
 end if

 call timab(148,1,tsec)

 if(nsppol==2)then
   message = 'nsppol=2 is still under development. Be careful when using it ...'
   MSG_COMMENT(message)
 end if

 band2tot_index =0
 bandtot_index=0
 band_index=0

!Add scissor shift to eigenenergies
 if (dtset%dfpt_sciss > tol6 ) then
   write(msg,'(a,f7.3,2a)')&
&   ' A scissor operator of ',dtset%dfpt_sciss*Ha_eV,' [eV] has been applied to the eigenenergies',ch10
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   ABI_ALLOCATE(eigen0tmp,(nkpt_rbz*mband*nsppol))
   ABI_ALLOCATE(eigenqtmp,(nkpt_rbz*mband*nsppol))
   eigen0tmp =   eigen0(:)
   eigenqtmp =   eigenq(:)
   eigen0 = zero
   eigenq = zero
 end if

 if(ieig2rf > 0) then
   eig2nkq(:,:,:,:,:,:,:) = zero
 end if
 if(present(eigbrd))then
   eigbrd(:,:,:,:,:,:,:) = zero
 end if

 if(xmpi_paral==1) then
   ABI_ALLOCATE(mpi_enreg%proc_distrb,(nkpt_rbz,mband,nsppol))
   ABI_ALLOCATE(nband_rbz,(nkpt_rbz*nsppol))
   if (allocated(mpi_enreg%my_kpttab)) then
     ABI_DEALLOCATE(mpi_enreg%my_kpttab)
   end if
   ABI_ALLOCATE(mpi_enreg%my_kpttab,(nkpt_rbz))
!  Assume the number of bands is the same for all k points.
   nband_rbz(:)=mband
   call distrb2(mband,nband_rbz,nkpt_rbz,mpi_enreg%nproc_cell,nsppol,mpi_enreg)
 end if

 icg2=0
 ipert1=1 ! Suppose that the situation is the same for all perturbations
 ABI_ALLOCATE(icg2_rbz,(nkpt_rbz,nsppol))
 do isppol=1,nsppol
   do ikpt=1,nkpt_rbz
     icg2_rbz(ikpt,isppol)=icg2
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,mband,isppol,me)) cycle
     icg2 = icg2 + npwar1(ikpt,ipert1)*nspinor*mband
   end do
 end do

 do isppol=1,nsppol
   do ikpt =1,nkpt_rbz

     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,mband,isppol,me)) then
       band2tot_index = band2tot_index + 2*mband**2
       bandtot_index = bandtot_index + mband
       cycle
     end if

     if(present(eigenq_fine))then
       write(std_out,*) 'Start of the energy denominator interpolation method.'
       nkpt_sub = 0
!      center is the k+q point around which we will average the kpt_fine
       center = hdr0%kptns(:,ikpt)+ dtset%qptn(:)

       call kptfine_av(center,dtset%qptrlatt,hdr_fine%kptns,hdr_fine%nkpt,kpt_fine_sub,nkpt_sub,wgt_sub)
       write(std_out,'(a,3f8.4,a,i3)') 'Number of k-points of the fine grid &
&       around the k+Q point ',center,' is:',nkpt_sub
       write(std_out,'(a,f10.5)') 'The sum of the weights of the k-points is: ',SUM(wgt_sub)
     end if

!    Add scissor shift to eigenenergies
     if (dtset%dfpt_sciss > tol6 ) then
       do iband=1,mband
         if (occ(iband+bandtot_index) < tol6) then
           eigen0(iband+bandtot_index) = eigen0tmp(iband+bandtot_index) + dtset%dfpt_sciss
           eigenq(iband+bandtot_index) = eigenqtmp(iband+bandtot_index) + dtset%dfpt_sciss
         else
           eigen0(iband+bandtot_index) = eigen0tmp(iband+bandtot_index)
           eigenq(iband+bandtot_index) = eigenqtmp(iband+bandtot_index)
         end if
       end do
     end if


     if(smdelta >0) then   !broadening
       if(.not.allocated(smdfun))  then
         ABI_ALLOCATE(smdfun,(mband,mband))
       end if
       smdfun(:,:) = zero
       do iband=1,mband
         eigen(iband) = eigen0(iband+bandtot_index)
         eigen_prime(iband) =eigenq(iband+bandtot_index)
       end do
       if(esmear>tol6) then
         call smeared_delta(eigen,eigen_prime,esmear,mband,smdelta,smdfun)
       end if
     end if
     icg2=icg2_rbz(ikpt,isppol)

     ipert1=1 ! Suppose all perturbations lead to the same number of planewaves
     npw1_k = npwar1(ikpt,ipert1)
     ABI_ALLOCATE(cwavef,(2,npw1_k*nspinor))
     ABI_ALLOCATE(cwavef2,(2,npw1_k*nspinor))
     ABI_ALLOCATE(gh,(2,npw1_k*nspinor))
     ABI_ALLOCATE(gh1,(2,npw1_k*nspinor))
     ABI_ALLOCATE(ghc,(2,npw1_k*nspinor))

     do iband=1,bdeigrf

!      If the k point and band belong to me, compute the contribution
       test_do_band=.true.
       if(mpi_enreg%proc_distrb(ikpt,iband,isppol)/=me)test_do_band=.false.

       if(test_do_band)then

         do ipert1=1,npert

           do idir1=1,3
             if(clflg(idir1,ipert1)==0)cycle
             istwf_k = istwfk_pert(ikpt,idir1,ipert1)

             do ipert2=1,npert
               do idir2=1,3
                 if(clflg(idir2,ipert2)==0)cycle

                 eig2_diar = zero ; eig2_diai = zero ; eigbrd_r = zero ; eigbrd_i = zero

                 do jband=1,mband
                   eig1_r1 = eigen1(2*jband-1+(iband-1)*2*mband+band2tot_index,idir1,ipert1)
                   eig1_r2 = eigen1(2*jband-1+(iband-1)*2*mband+band2tot_index,idir2,ipert2)
                   eig1_i1 = eigen1(2*jband+(iband-1)*2*mband+band2tot_index,idir1,ipert1)
                   eig1_i2 = - eigen1(2*jband+(iband-1)*2*mband+band2tot_index,idir2,ipert2) !the negative sign is from the CC
!                  If no interpolation, fallback on to the previous
!                  implementation
                   if(.not. present(eigenq_fine))then
                     deltae=eigenq(jband+bandtot_index)-eigen0(iband+bandtot_index)
                   end if
                   ar=eig1_r1*eig1_r2-eig1_i1*eig1_i2
                   ai=eig1_r1*eig1_i2+eig1_i1*eig1_r2

!                  Sum over all active space to retrieve the diagonal gauge
                   if(ieig2rf == 1 .or. ieig2rf ==2 ) then
!                    if(abs(deltae)>etol) then ! This is commented because
!                    there is no problem with divergencies with elph2_imag != 0
                     if( present(eigenq_fine))then
                       den_av = zero
                       wgt_int = zero
                       do ikpt2=1,nkpt_sub
                         deltae=eigenq_fine(jband,kpt_fine_sub(ikpt2),1)&
&                         -eigen0(iband+bandtot_index)
                         den_av = den_av-(wgt_sub(ikpt2)*deltae)/(deltae**2+elph2_imagden**2)
                         wgt_int = wgt_int+wgt_sub(ikpt2)
                       end do
                       den = den_av/wgt_int
                     else
                       if(abs(elph2_imagden) < etol) then
                         if(abs(deltae)>etol) then
                           den=-one/(deltae**2+elph2_imagden**2)
                         else
                           den= zero
                         end if
                       else
                         den=-one/(deltae**2+elph2_imagden**2)
                       end if
                     end if

!                    The following should be the most general implementation of the presence of elph2_imagden
!                    eig2_diar=eig2_diar+(ar*deltae+ai*elph2_imagden)*den
!                    eig2_diai=eig2_diai+(ai*deltae-ar*elph2_imagden)*den
!                    This gives back the implementation without elph2_imagden
!                    eig2_diar=eig2_diar+ar*deltae*den
!                    eig2_diai=eig2_diai+ai*deltae*den
!                    This is what Samuel had implemented
!                    eig2_diar=eig2_diar+ar*deltae*den
!                    eig2_diai=eig2_diai+ai*elph2_imagden*den
!                    Other possibility : throw away the broadening part, that is actually treated separately.
                     if( present(eigenq_fine))then
                       eig2_diar=eig2_diar+ar*den
                       eig2_diai=eig2_diai+ai*den
                     else
                       eig2_diar=eig2_diar+ar*deltae*den
                       eig2_diai=eig2_diai+ai*deltae*den
!DBSP
!                       if (iband+band_index==2 .and. ikpt==1 .and. idir1==1 .and. ipert1==1 .and. idir2==1 .and. ipert2==1) then
!                         write(message2,*) 'eig2_diar1=',eig2_diar,' ar=',ar,' deltae=',deltae,' den=',den
!                         call wrtout(std_out,message2,'PERS')
!                       endif
!END

                     end if
                   end if ! ieig2rf==1 or 2

                   if(present(eigbrd))then
                     if(smdelta >0) then   !broadening
                       eigbrd_r = eigbrd_r + ar*smdfun(iband,jband)
                       eigbrd_i = eigbrd_i + ai*smdfun(iband,jband)
                     end if
                   end if

                 end do !jband

!                Add the contribution of non-active bands, if DFPT calculation (= Sternheimer)
                 if(ieig2rf == 1 .or. ieig2rf ==3 .or. ieig2rf ==4 .or. ieig2rf==5 ) then
!                  if(ieig2rf == 1   ) then

                   dotr=zero ; doti=zero
                   dot2r=zero ; dot2i=zero
                   dot3r=zero ; dot3i=zero


                   cwavef(:,:) = cg1_pert(:,1+(iband-1)*npw1_k*nspinor+icg2:iband*npw1_k*nspinor+icg2,idir2,ipert2)
                   cwavef2(:,:)= cg1_pert(:,1+(iband-1)*npw1_k*nspinor+icg2:iband*npw1_k*nspinor+icg2,idir1,ipert1)
                   gh1(:,:)    = gh1c_pert(:,1+(iband-1)*npw1_k*nspinor+icg2:iband*npw1_k*nspinor+icg2,idir1,ipert1)
                   gh(:,:)     = gh1c_pert(:,1+(iband-1)*npw1_k*nspinor+icg2:iband*npw1_k*nspinor+icg2,idir2,ipert2)
                   ghc(:,:)    = gh0c1_pert(:,1+(iband-1)*npw1_k*nspinor+icg2:iband*npw1_k*nspinor+icg2,idir1,ipert1)

!                  The first two dotprod corresponds to:  <Psi(1)nkq|H(1)k+q,k|Psi(0)nk> and <Psi(0)nk|H(1)k,k+q|Psi(1)nkq>
!                  They are calculated using wavefunctions <Psi(1)| that are orthogonal to the active space.
                   call dotprod_g(dotr,doti,istwf_k,npw1_k*nspinor,2,cwavef,gh1,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
                   call dotprod_g(dot2r,dot2i,istwf_k,npw1_k*nspinor,2,gh,cwavef2,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

!                  This dotprod corresponds to : <Psi(1)nkq|H(0)k+q- E(0)nk|Psi(1)nkq>
!                  It is calculated using wavefunctions that are orthogonal to the active space.
!                  Should work for metals. (But adiabatic approximation is bad in this case...)
                   call dotprod_g(dot3r,dot3i,istwf_k,npw1_k*nspinor,2,cwavef,ghc,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

                   eig2_diar= eig2_diar + dotr + dot2r + dot3r
                   eig2_diai= eig2_diai + doti + dot2i + dot3i

                 end if

!                Store the contribution
                 if(ieig2rf > 0) then
                   eig2nkq(1,iband+band_index,ikpt,idir1,ipert1,idir2,ipert2) = eig2_diar
                   eig2nkq(2,iband+band_index,ikpt,idir1,ipert1,idir2,ipert2) = eig2_diai
                 end if

                 if(present(eigbrd))then
                   if(smdelta >0) then   !broadening
                     eigbrd(1,iband+band_index,ikpt,idir1,ipert1,idir2,ipert2) = eigbrd_r
                     eigbrd(2,iband+band_index,ikpt,idir1,ipert1,idir2,ipert2) = eigbrd_i
                   end if
                 end if

               end do !idir2
             end do !ipert2
           end do  !idir1
         end do   !ipert1

       end if ! Selection of processor

     end do !iband

     ABI_DEALLOCATE(cwavef)
     ABI_DEALLOCATE(cwavef2)
     ABI_DEALLOCATE(gh)
     ABI_DEALLOCATE(gh1)
     ABI_DEALLOCATE(ghc)
     band2tot_index = band2tot_index + 2*mband**2
     bandtot_index = bandtot_index + mband

     if(present(eigenq_fine))then
       ABI_DEALLOCATE(kpt_fine_sub) ! Deallocate the variable
       ABI_DEALLOCATE(wgt_sub)
     end if

   end do    !ikpt
   band_index = band_index + mband
 end do !isppol

!Accumulate eig2nkq and/or eigbrd
 if(xmpi_paral==1) then
   if(ieig2rf == 1 .or. ieig2rf == 2) then
     call xmpi_sum(eig2nkq,spaceworld,ierr)
     if (dtset%dfpt_sciss > tol6 ) then
       call xmpi_sum(eigen0,spaceworld,ierr)
       call xmpi_sum(eigenq,spaceworld,ierr)
     end if
   end if
   if(present(eigbrd) .and. (ieig2rf == 1 .or. ieig2rf == 2))then
     if(smdelta >0) then
       call xmpi_sum(eigbrd,spaceworld,ierr)
     end if
   end if
   ABI_DEALLOCATE(nband_rbz)
   ABI_DEALLOCATE(mpi_enreg%proc_distrb)
   ABI_DEALLOCATE(mpi_enreg%my_kpttab)
 end if

 if(ieig2rf==1 .or. ieig2rf==2 ) then
   write(ab_out,'(a)')' Components of second-order derivatives of the electronic energy, EIGR2D.'
   write(ab_out,'(a)')' For automatic tests, printing the matrix for the first k-point, first band, first atom.'
   do idir1=1,3
     do idir2=1,3
       ar=eig2nkq(1,1,1,idir1,1,idir2,1) ; if(abs(ar)<tol10)ar=zero
       ai=eig2nkq(2,1,1,idir1,1,idir2,1) ; if(abs(ai)<tol10)ai=zero
       write (ab_out,'(4i4,2es20.10)') idir1,1,idir2,1,ar,ai
     end do ! idir2
   end do ! idir1
 end if

 if(present(eigbrd))then
   if(smdelta >0) then   !broadening
     write(ab_out,'(a)')' '
     write(ab_out,'(a)')' Components of second-order derivatives of the electronic energy, EIGI2D.'
     write(ab_out,'(a)')' For automatic tests, printing the matrix for the first k-point, first band, first atom.'
     do idir1=1,3
       do idir2=1,3
         ar=eigbrd(1,1,1,idir1,1,idir2,1) ; if(abs(ar)<tol10)ar=zero
         ai=eigbrd(2,1,1,idir1,1,idir2,1) ; if(abs(ai)<tol10)ai=zero
         write (ab_out,'(4i4,2es20.10)') idir1,1,idir2,1,ar,ai
       end do
     end do !nband
   end if
 end if

 if(allocated(smdfun))  then
   ABI_DEALLOCATE(smdfun)
 end if
 ABI_DEALLOCATE(icg2_rbz)
 if(present(eigenq_fine))then
   ABI_DEALLOCATE(center)
 end if
 if (dtset%dfpt_sciss > tol6 ) then
   ABI_DEALLOCATE(eigen0tmp)
   ABI_DEALLOCATE(eigenqtmp)
 end if

 call timab(148,2,tsec)

!DEBUG
!write(std_out,*)' eig2stern: exit'
!ENDDEBUG

end subroutine eig2stern
!!***

!!****f* m_eig2d/eig2tot
!! NAME
!! eig2tot
!!
!! FUNCTION
!! This routine calculates the second-order eigenvalues.
!! The output eig2nkq is this quantity for the input k points.
!!
!! INPUTS
!!  bdeigrf = number of bands for which to calculate the second-order eigenvalues.
!!  clflg(3,mpert)= array on calculated perturbations for eig2rf.
!!  dim_eig2nkq = 1 if eig2nkq is to be computed.
!!  eigbrd(2,mband*nsppol,nkpt,3,natom,3,natom) = broadening factors for the
!!            electronic eigenvalues (optional).
!!  eigen0(nkpt_rbz*mband*nsppol) = 0-order eigenvalues at all K-points:
!!            <k,n'|H(0)|k,n'> (hartree).
!!  eigenq(nkpt_rbz*mband*nsppol) = 0-order eigenvalues at all shifted K-points:
!!            <k+Q,n'|H(0)|k+Q,n'> (hartree).
!!  eigen1(nkpt_rbz*2*nsppol*mband**2,3,mpert) = matrix of first-order:
!!            <k+Q,n'|H(1)|k,n> (hartree) (calculated in dfpt_cgwf).
!!  eig2nkq(2,mband*nsppol,nkpt,3,natom,3,natom*dim_eig2nkq) = second derivatives of
!!            the electronic eigenvalues.
!!  elph2_imagden = imaginary part of the denominator of the sum-over-state expression
!!            for the electronic eigenenergy shift due to second-order electron-phonon
!!            interation.
!!  ieig2rf = integer for calculation type.
!!  indsym(4,nsym,natom) = indirect indexing array for atom labels
!!            (not used yet, but will be used with symmetries).
!!  mband = maximum number of bands.
!!  mpert = maximum number of perturbations.
!!  natom = number of atoms in the unit cell.
!!  npert = number of phonon perturbations, without taking into account directions:
!!            natom.
!!  nsym = number of symmetries (not used yet).
!!  mpi_enreg = information about MPI parallelization.
!!  nkpt_rbz = number of k-points for each perturbation.
!!  nsppol = 1 for unpolarized, 2 for spin-polarized.
!!  smdelta = integer controling the calculation of electron lifetimes.
!!  symq(4,2,nsym) = 1 if symmetry preserves present qpoint. From littlegroup_q (not used yet).
!!  symrec(3,3,nsym) = 3x3 matrices of the group symmetries (reciprocal space)
!!            (not used yet).
!!  symrel(3,3,nsym) = array containing the symmetries in real space (not used yet).
!!  timrev = 1 if time-reversal preserves the q wavevector; 0 otherwise
!!            (not in use yet).
!!  dtset = OPTIONAL, dataset structure containing the input variable of the
!!            calculation. This is required to use the k-interpolation routine.
!!  eigenq_fine(mband_fine,mkpt_fine,nsppol_fine) = OPTIONAL, 0-order eigenvalues
!!            at all shifted K-points: <k+Q,n'|H(0)|k+Q,n'> (hartree) of the
!!            fine grid. This information is read from the WF dense k-grid file.
!!  hdr_fine = OPTIONAL, header of the WF file of the fine k-point grid. This
!!            variable is required for the k-interpolation routine.
!!  hdr0     = header of the GS WF file of the corse k-point grid.
!!
!!
!! OUTPUT
!!  eig2nkq(2,mband*nsppol,nkpt_rbz,3,npert,3,npert)= diagonal part of the
!!            second-order eigenvalues: E^{(2),diag}_{k,q,j}.
!!  eigbrd(2,mband*nsppol,nkpt_rbz,3,npert,3,npert)= OPTIONAL, array containing the
!!            electron lifetimes.
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      crystal_free,crystal_init,ddb_hdr_free,ddb_hdr_init,ddb_hdr_open_write
!!      distrb2,ebands_free,ebands_init,eigr2d_free,eigr2d_init,eigr2d_ncwrite
!!      fan_free,fan_init,fan_ncwrite,gkk_free,gkk_init,gkk_ncwrite,kptfine_av
!!      outbsd,smeared_delta,timab,xmpi_sum
!!
!! SOURCE

subroutine eig2tot(dtfil,xred,psps,pawtab,natom,bdeigrf,clflg,dim_eig2nkq,eigen0,eigenq,eigen1,eig2nkq,&
&  elph2_imagden,esmear,ieig2rf,mband,mpert,npert,mpi_enreg,doccde,&
&  nkpt_rbz,nsppol,smdelta,rprimd,dtset,occ_rbz,hdr0,eigbrd,eigenq_fine,hdr_fine)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bdeigrf,dim_eig2nkq,ieig2rf,mband,mpert,natom,nkpt_rbz
 integer,intent(in) :: npert,nsppol,smdelta
 real(dp),intent(in) :: elph2_imagden,esmear
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type), intent(in) :: dtfil
 type(pseudopotential_type), intent(inout) :: psps
!arrays
 type(dataset_type), intent(in) :: dtset
 integer,intent(in) :: clflg(3,mpert)
 real(dp),intent(in) :: doccde(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: eigen0(nkpt_rbz*mband*nsppol)
 real(dp),intent(in) :: eigen1(nkpt_rbz*2*nsppol*mband**2,3,mpert)
 real(dp),intent(in) :: eigenq(nkpt_rbz*mband*nsppol)
 real(dp),intent(in) :: occ_rbz(mband*nkpt_rbz*nsppol)
 real(dp),intent(inout) :: eig2nkq(2,mband*nsppol,nkpt_rbz,3,npert,3,npert*dim_eig2nkq)
 real(dp),intent(in) :: rprimd(3,3),xred(3,natom)
 real(dp),intent(inout),optional :: eigbrd(2,mband*nsppol,nkpt_rbz,3,npert,3,npert)
 real(dp),intent(in),pointer,optional :: eigenq_fine(:,:,:)
 type(pawtab_type), intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
 type(hdr_type),intent(in) :: hdr0
 type(hdr_type),intent(in),optional :: hdr_fine

!Local variables-------------------------------
!tolerance for non degenerated levels
!scalars
 integer :: band2tot_index,band_index,bantot,bandtot_index,iband,idir1,idir2
 integer :: ikpt,ipert1,ipert2,isppol,jband,nkpt_sub,ikpt2,unitout,ncid
!integer :: ipw
 character(len=fnlen) :: dscrpt,fname
 integer :: master,me,spaceworld,ierr
! real(dp),parameter :: etol=1.0d-6
 real(dp),parameter :: etol=1.0d-7
!real(dp),parameter :: etol=zero
 real(dp) :: ar,ai,deltae,den,eig1_i1,eig1_i2,eigen_corr
 real(dp) :: eig1_r1,eig1_r2,eig2_diai,den_av
 real(dp) :: eig2_diar,eigbrd_i,eigbrd_r,wgt_int
 character(len=500) :: message
 logical :: remove_inv,test_do_band
 type(crystal_t) :: Crystal
 type(ebands_t)  :: Bands
 type(eigr2d_t)  :: eigr2d,eigi2d
 type(fan_t)     :: fan2d
 type(gkk_t)     :: gkk2d
 type(ddb_hdr_type) :: ddb_hdr
!arrays
 integer, allocatable :: nband_rbz(:)
 integer,pointer      :: kpt_fine_sub(:)
 real(dp)             :: tsec(2)
 real(dp),allocatable :: center(:)
 real(dp) :: eigen(mband*nsppol),eigen_prime(mband*nsppol)
 real(dp),allocatable :: fan(:,:,:,:,:,:,:)
 real(dp),allocatable :: gkk(:,:,:,:,:)
 real(dp),allocatable :: eig2nkq_tmp(:,:,:,:,:,:,:)
 real(dp),allocatable :: smdfun(:,:)
 real(dp),pointer     :: wgt_sub(:)

! *********************************************************************

!Init parallelism
 master =0
 spaceworld=mpi_enreg%comm_cell
 me=mpi_enreg%me_kpt
!DEBUG
!write(std_out,*)' eig2tot : enter '
!write(std_out,*)' mband=',mband
!write(std_out,*)' nsppol=',nsppol
!write(std_out,*)' nkpt_rbz=',nkpt_rbz
!write(std_out,*)' npert=',npert
!ENDDEBUG

!Init interpolation method
 if(present(eigenq_fine))then
   ABI_ALLOCATE(center,(3))
 end if

 call timab(148,1,tsec)

 if(nsppol==2)then
   message = 'nsppol=2 is still under development. Be careful when using it ...'
   MSG_COMMENT(message)
 end if

 band2tot_index =0
 bandtot_index=0
 band_index=0

 if(xmpi_paral==1) then
   ABI_ALLOCATE(mpi_enreg%proc_distrb,(nkpt_rbz,mband,nsppol))
   ABI_ALLOCATE(nband_rbz,(nkpt_rbz*nsppol))
   if (allocated(mpi_enreg%my_kpttab)) then
     ABI_DEALLOCATE(mpi_enreg%my_kpttab)
   end if
   ABI_ALLOCATE(mpi_enreg%my_kpttab,(nkpt_rbz))
!  Assume the number of bands is the same for all k points.
   nband_rbz(:)=mband
   call distrb2(mband,nband_rbz,nkpt_rbz,mpi_enreg%nproc_cell,nsppol,mpi_enreg)
 end if

 if(ieig2rf == 4 ) then
   ABI_MALLOC_OR_DIE(fan,(2*mband*nsppol,dtset%nkpt,3,natom,3,natom*dim_eig2nkq,mband), ierr)
   fan(:,:,:,:,:,:,:) = zero
   ABI_MALLOC_OR_DIE(eig2nkq_tmp,(2,mband*nsppol,dtset%nkpt,3,natom,3,natom*dim_eig2nkq), ierr)
   eig2nkq_tmp(:,:,:,:,:,:,:) = zero
!  This is not efficient because double the memory. Alternative: use buffer and
!  print part by part.
   eig2nkq_tmp = eig2nkq
   if(present(eigbrd))then
     eigbrd(:,:,:,:,:,:,:)=zero
   end if
   eigen_corr = 0
 end if

 if(ieig2rf == 5 ) then
   ABI_MALLOC_OR_DIE(gkk,(2*mband*nsppol,dtset%nkpt,3,natom,mband), ierr)
   gkk(:,:,:,:,:) = zero
   ABI_MALLOC_OR_DIE(eig2nkq_tmp,(2,mband*nsppol,dtset%nkpt,3,natom,3,natom*dim_eig2nkq), ierr)
   eig2nkq_tmp(:,:,:,:,:,:,:) = zero
!  This is not efficient because double the memory. Alternative: use buffer and
!  print part by part.
   eig2nkq_tmp = eig2nkq
   if(present(eigbrd))then
     eigbrd(:,:,:,:,:,:,:)=zero
   end if
   eigen_corr = 0
 end if

 do isppol=1,nsppol
   do ikpt =1,nkpt_rbz

     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,mband,isppol,me)) then
       band2tot_index = band2tot_index + 2*mband**2
       bandtot_index = bandtot_index + mband
       cycle
     end if

     if(present(eigenq_fine))then
       write(std_out,*) 'Start of the energy denominator interpolation method.'
       nkpt_sub = 0
!      center is the k+q point around which we will average the kpt_fine
       center = hdr0%kptns(:,ikpt)+ dtset%qptn(:)

       call kptfine_av(center,dtset%qptrlatt,hdr_fine%kptns,hdr_fine%nkpt,kpt_fine_sub,nkpt_sub,wgt_sub)
       write(std_out,'(a,3f8.4,a,i3)') 'Number of k-points of the fine grid &
&       around the k+Q point ',center,' is:',nkpt_sub
       write(std_out,'(a,f10.5)') 'The sum of the weights of the k-points is: ',SUM(wgt_sub)
     end if

     if(smdelta >0) then   !broadening
       if(.not.allocated(smdfun))  then
         ABI_ALLOCATE(smdfun,(mband,mband))
       end if
       smdfun(:,:) = zero
       do iband=1,mband
         eigen(iband) = eigen0(iband+bandtot_index)
         eigen_prime(iband) =eigenq(iband+bandtot_index)
       end do
       if(esmear>tol6) then
         call smeared_delta(eigen,eigen_prime,esmear,mband,smdelta,smdfun)
       end if
     end if

     ipert1=1 ! Suppose all perturbations lead to the same number of planewaves

     do iband=1,bdeigrf

!      If the k point and band belong to me, compute the contribution
       test_do_band=.true.
       if(mpi_enreg%proc_distrb(ikpt,iband,isppol)/=me)test_do_band=.false.

       if(test_do_band)then
!        ------------------------------------------------------------------------------------------------------!
!        ------- ieig2rf ==3 : Non dynamic traditional AHC theory with Sternheimer (computed in eig2stern.F90)-!
!        ------------------------------------------------------------------------------------------------------!
!        Note that ieig2rf==4 and ieig2rf==5 also goes into that part only for later printing of the ZPR in the ouput of abinit
!        later in the code
         if(ieig2rf==3 .or. ieig2rf==4 .or. ieig2rf==5) then
           do ipert1=1,npert
             do idir1=1,3
               if(clflg(idir1,ipert1)==0) cycle
               do ipert2=1,npert
                 do idir2=1,3
                   if(clflg(idir2,ipert2)==0)cycle
                   eig2_diar = zero ; eig2_diai = zero ; eigbrd_r = zero ; eigbrd_i = zero
                   do jband=1,mband
                     eig1_r1 = eigen1(2*jband-1+(iband-1)*2*mband+band2tot_index,idir1,ipert1)
                     eig1_r2 = eigen1(2*jband-1+(iband-1)*2*mband+band2tot_index,idir2,ipert2)
                     eig1_i1 = eigen1(2*jband+(iband-1)*2*mband+band2tot_index,idir1,ipert1)
                     eig1_i2 = - eigen1(2*jband+(iband-1)*2*mband+band2tot_index,idir2,ipert2) !the negative sign is from the CC
!                    If no interpolation, fallback on to the previous
!                    implementation
                     if(.not. present(eigenq_fine))then
                       deltae=eigenq(jband+bandtot_index)-eigen0(iband+bandtot_index)
                     end if
                     ar=eig1_r1*eig1_r2-eig1_i1*eig1_i2
                     ai=eig1_r1*eig1_i2+eig1_i1*eig1_r2

!                    Sum over all active space to retrieve the diagonal gauge
!                    if(abs(deltae)>etol) then ! This is commented because
!                    there is no problem with divergencies with elph2_imag != 0
                     if( present(eigenq_fine))then
                       den_av = zero
                       wgt_int = zero
                       do ikpt2=1,nkpt_sub
                         deltae=eigenq_fine(jband,kpt_fine_sub(ikpt2),1)&
&                         -eigen0(iband+bandtot_index)
                         den_av = den_av-(wgt_sub(ikpt2)*deltae)/(deltae**2+elph2_imagden**2)
                         wgt_int = wgt_int+wgt_sub(ikpt2)
                       end do
                       den = den_av/wgt_int
                     else
                       if(abs(elph2_imagden) < etol) then
                         if(abs(deltae)>etol) then
                           den=-one/(deltae**2+elph2_imagden**2)
                         else
                           den= zero
                         end if
                       else
                         den=-one/(deltae**2+elph2_imagden**2)
                       end if
                     end if

                     if( present(eigenq_fine))then
                       eig2_diar=eig2_diar+ar*den
                       eig2_diai=eig2_diai+ai*den
                     else
                       eig2_diar=eig2_diar+ar*deltae*den
                       eig2_diai=eig2_diai+ai*deltae*den
                     end if

                     if(present(eigbrd))then
                       if(smdelta >0) then   !broadening
                         eigbrd_r = eigbrd_r + ar*smdfun(iband,jband)
                         eigbrd_i = eigbrd_i + ai*smdfun(iband,jband)
                       end if
                     end if
                   end do !jband

!                  Store the contribution
                   eig2nkq(1,iband+band_index,ikpt,idir1,ipert1,idir2,ipert2) = &
&                   eig2nkq(1,iband+band_index,ikpt,idir1,ipert1,idir2,ipert2) + eig2_diar
                   eig2nkq(2,iband+band_index,ikpt,idir1,ipert1,idir2,ipert2) = &
&                   eig2nkq(2,iband+band_index,ikpt,idir1,ipert1,idir2,ipert2) + eig2_diai

                   if(present(eigbrd))then
                     if(smdelta >0) then   !broadening
                       eigbrd(1,iband+band_index,ikpt,idir1,ipert1,idir2,ipert2) = eigbrd_r
                       eigbrd(2,iband+band_index,ikpt,idir1,ipert1,idir2,ipert2) = eigbrd_i
                     end if
                   end if

                 end do !idir2
               end do !ipert2
             end do  !idir1
           end do   !ipert1
         end if !ieig2rf 3

!        -------------------------------------------------------------------------------------------!
!        ------- ieig2rf ==4  Dynamic AHC using second quantization and Sternheimer from eig2stern -!
!        -------------------------------------------------------------------------------------------!
         if(ieig2rf ==4 ) then
           do ipert1=1,npert
             do idir1=1,3
               if(clflg(idir1,ipert1)==0) cycle
               do ipert2=1,npert
                 do idir2=1,3
                   if(clflg(idir2,ipert2)==0)cycle
                   do jband=1,mband
                     eig1_r1 = eigen1(2*jband-1+(iband-1)*2*mband+band2tot_index,idir1,ipert1)
                     eig1_r2 = eigen1(2*jband-1+(iband-1)*2*mband+band2tot_index,idir2,ipert2)
                     eig1_i1 = eigen1(2*jband+(iband-1)*2*mband+band2tot_index,idir1,ipert1)
                     eig1_i2 = - eigen1(2*jband+(iband-1)*2*mband+band2tot_index,idir2,ipert2) !the negative sign is from the CC
                     ar=eig1_r1*eig1_r2-eig1_i1*eig1_i2
                     ai=eig1_r1*eig1_i2+eig1_i1*eig1_r2
!                  Store the contribution
                     fan(2*iband-1+2*band_index,ikpt,idir1,ipert1,idir2,ipert2,jband) = &
&                     fan(2*iband-1+2*band_index,ikpt,idir1,ipert1,idir2,ipert2,jband) + ar
                     fan(2*iband+2*band_index,ikpt,idir1,ipert1,idir2,ipert2,jband) = &
&                     fan(2*iband+2*band_index,ikpt,idir1,ipert1,idir2,ipert2,jband) + ai
                   end do !jband
                 end do !idir2
               end do !ipert2
             end do  !idir1
           end do   !ipert1
         end if !ieig2rf 4
!        --------------------------------------------------------------------------------!
!        ------- ieig2rf ==5  Dynamic AHC with Sternheimer from eig2stern but print GKK -!
!        --------------------------------------------------------------------------------!
         if(ieig2rf ==5 ) then
           do ipert1=1,npert
             do idir1=1,3
               if(clflg(idir1,ipert1)==0) cycle
               do jband=1,mband
                 eig1_r1 = eigen1(2*jband-1+(iband-1)*2*mband+band2tot_index,idir1,ipert1)
                 eig1_i1 = eigen1(2*jband+(iband-1)*2*mband+band2tot_index,idir1,ipert1)
!              Store the contribution
                 gkk(2*iband-1+2*band_index,ikpt,idir1,ipert1,jband) = &
&                 gkk(2*iband-1+2*band_index,ikpt,idir1,ipert1,jband) + eig1_r1
                 gkk(2*iband+2*band_index,ikpt,idir1,ipert1,jband) = &
&                 gkk(2*iband+2*band_index,ikpt,idir1,ipert1,jband) + eig1_i1
               end do !jband
             end do  !idir1
           end do   !ipert1
         end if !ieig2rf 5
       end if ! Selection of processor
     end do !iband

     band2tot_index = band2tot_index + 2*mband**2
     bandtot_index = bandtot_index + mband

     if(present(eigenq_fine))then
       ABI_DEALLOCATE(kpt_fine_sub) ! Deallocate the variable
       ABI_DEALLOCATE(wgt_sub)
     end if

   end do    !ikpt
   band_index = band_index + mband
 end do !isppol

!Accumulate eig2nkq and/or eigbrd
 if(xmpi_paral==1) then
   if(ieig2rf == 3) then
     call xmpi_sum(eig2nkq,spaceworld,ierr)
   end if
   if(ieig2rf == 4) then
     call xmpi_sum(eig2nkq,spaceworld,ierr)
     call xmpi_sum(eig2nkq_tmp,spaceworld,ierr)
     call xmpi_sum(fan,spaceworld,ierr)
   end if
   if(ieig2rf == 5) then
     call xmpi_sum(eig2nkq,spaceworld,ierr)
     call xmpi_sum(eig2nkq_tmp,spaceworld,ierr)
     call xmpi_sum(gkk,spaceworld,ierr)
   end if
   if(present(eigbrd) .and. (ieig2rf == 3 .or. ieig2rf == 4 .or. ieig2rf == 5))then
     if(smdelta >0) then
       call xmpi_sum(eigbrd,spaceworld,ierr)
     end if
   end if
   ABI_DEALLOCATE(nband_rbz)
   ABI_DEALLOCATE(mpi_enreg%proc_distrb)
   ABI_DEALLOCATE(mpi_enreg%my_kpttab)
 end if

 if(ieig2rf > 2) then
   write(ab_out,'(a)')' Components of second-order derivatives of the electronic energy, EIGR2D.'
   write(ab_out,'(a)')' For automatic tests, printing the matrix for the first k-point, first band, first atom.'
   band_index = 0
   do isppol=1,dtset%nsppol
     do idir1=1,3
       do idir2=1,3
         ar=eig2nkq(1,1+band_index,1,idir1,1,idir2,1) ; if(abs(ar)<tol10)ar=zero
         ai=eig2nkq(2,1+band_index,1,idir1,1,idir2,1) ; if(abs(ai)<tol10)ai=zero
         write (ab_out,'(4i4,2es20.10)') idir1,1,idir2,1,ar,ai
       end do ! idir2
     end do ! idir1
     band_index = band_index + mband
     write(ab_out,'(a)')' '
   end do
 end if

 if(present(eigbrd))then
   if(smdelta >0) then   !broadening
     write(ab_out,'(a)')' Components of second-order derivatives of the electronic energy, EIGI2D.'
     write(ab_out,'(a)')' For automatic tests, printing the matrix for the first k-point, first band, first atom.'
     band_index = 0
     do isppol=1,dtset%nsppol
       do idir1=1,3
         do idir2=1,3
           ar=eigbrd(1,1+band_index,1,idir1,1,idir2,1) ; if(abs(ar)<tol10)ar=zero
           ai=eigbrd(2,1+band_index,1,idir1,1,idir2,1) ; if(abs(ai)<tol10)ai=zero
           write (ab_out,'(4i4,2es20.10)') idir1,1,idir2,1,ar,ai
         end do
       end do
       band_index = band_index + mband
       write(ab_out,'(a)')' '
     end do
   end if
 end if

 if(allocated(smdfun))  then
   ABI_DEALLOCATE(smdfun)
 end if
 if(present(eigenq_fine))then
   ABI_DEALLOCATE(center)
 end if

 master=0
 if (me==master) then
!  print _EIGR2D file for this perturbation in the case of ieig2rf 3 or 4 or 5
   if (ieig2rf == 3 .or. ieig2rf == 4 .or. ieig2rf == 5) then

     dscrpt=' Note : temporary (transfer) database '
     unitout = dtfil%unddb

     call ddb_hdr_init(ddb_hdr,dtset,psps,pawtab,DDB_VERSION,dscrpt,&
&     1,xred=xred,occ=occ_rbz)

     call ddb_hdr_open_write(ddb_hdr, dtfil%fnameabo_eigr2d, unitout)

     call ddb_hdr_free(ddb_hdr)

   end if
   if(ieig2rf == 3 ) then
     call outbsd(bdeigrf,dtset,eig2nkq,dtset%natom,nkpt_rbz,unitout)
   end if
   if(ieig2rf == 4 .or. ieig2rf == 5 ) then
     call outbsd(bdeigrf,dtset,eig2nkq_tmp,dtset%natom,nkpt_rbz,unitout)
   end if
!  Output of the EIGR2D.nc file.
   fname = strcat(dtfil%filnam_ds(4),"_EIGR2D.nc")
!  Crystalline structure.
   remove_inv=.false.
   if(dtset%nspden==4 .and. dtset%usedmft==1) remove_inv=.true.
   call crystal_init(dtset%amu_orig(:,1),Crystal,dtset%spgroup,dtset%natom,dtset%npsp,psps%ntypat, &
&   dtset%nsym,rprimd,dtset%typat,xred,dtset%ziontypat,dtset%znucl,1,&
&   dtset%nspden==2.and.dtset%nsppol==1,remove_inv,hdr0%title,&
&   dtset%symrel,dtset%tnons,dtset%symafm)
!  Electronic band energies.
   bantot= dtset%mband*dtset%nkpt*dtset%nsppol
   call ebands_init(bantot,Bands,dtset%nelect,doccde,eigen0,hdr0%istwfk,hdr0%kptns,&
&   hdr0%nband, hdr0%nkpt,hdr0%npwarr,hdr0%nsppol,hdr0%nspinor,&
&   hdr0%tphysel,hdr0%tsmear,hdr0%occopt,hdr0%occ,hdr0%wtk,&
&   hdr0%charge, hdr0%kptopt, hdr0%kptrlatt_orig, hdr0%nshiftk_orig, hdr0%shiftk_orig, &
&   hdr0%kptrlatt, hdr0%nshiftk, hdr0%shiftk)

!  Second order derivative EIGR2D (real and Im)
   if(ieig2rf == 3 ) then
     call eigr2d_init(eig2nkq,eigr2d,dtset%mband,hdr0%nsppol,nkpt_rbz,dtset%natom)
   end if
   if(ieig2rf == 4 .or. ieig2rf == 5 ) then
     call eigr2d_init(eig2nkq_tmp,eigr2d,dtset%mband,hdr0%nsppol,nkpt_rbz,dtset%natom)
   end if
#ifdef HAVE_NETCDF
   NCF_CHECK_MSG(nctk_open_create(ncid, fname, xmpi_comm_self), "Creating EIGR2D file")
   NCF_CHECK(crystal%ncwrite(ncid))
   NCF_CHECK(ebands_ncwrite(Bands, ncid))
   call eigr2d_ncwrite(eigr2d,dtset%qptn(:),dtset%wtq,ncid)
   NCF_CHECK(nf90_close(ncid))
#else
   ABI_UNUSED(ncid)
#endif

!  print _FAN file for this perturbation. Note that the Fan file will only be produced if
!  abinit is compiled with netcdf.
   if(ieig2rf == 4 ) then
!    Output of the Fan.nc file.
#ifdef HAVE_NETCDF
     fname = strcat(dtfil%filnam_ds(4),"_FAN.nc")
     call fan_init(fan,fan2d,dtset%mband,hdr0%nsppol,nkpt_rbz,dtset%natom)
     NCF_CHECK_MSG(nctk_open_create(ncid, fname, xmpi_comm_self), "Creating FAN file")
     NCF_CHECK(crystal%ncwrite(ncid))
     NCF_CHECK(ebands_ncwrite(Bands, ncid))
     call fan_ncwrite(fan2d,dtset%qptn(:),dtset%wtq, ncid)
     NCF_CHECK(nf90_close(ncid))
#else
     MSG_ERROR("Dynamical calculation with ieig2rf 4 only work with NETCDF support.")
     ABI_UNUSED(ncid)
#endif
     ABI_DEALLOCATE(fan)
     ABI_DEALLOCATE(eig2nkq_tmp)
   end if
!  print _GKK.nc file for this perturbation. Note that the GKK file will only be produced if
!  abinit is compiled with netcdf.
   if(ieig2rf == 5 ) then
!    Output of the GKK.nc file.
#ifdef HAVE_NETCDF
     fname = strcat(dtfil%filnam_ds(4),"_GKK.nc")
     call gkk_init(gkk,gkk2d,dtset%mband,hdr0%nsppol,nkpt_rbz,dtset%natom,3)
     NCF_CHECK_MSG(nctk_open_create(ncid, fname, xmpi_comm_self), "Creating GKK file")
     NCF_CHECK(crystal%ncwrite(ncid))
     NCF_CHECK(ebands_ncwrite(Bands, ncid))
     call gkk_ncwrite(gkk2d,dtset%qptn(:),dtset%wtq, ncid)
     NCF_CHECK(nf90_close(ncid))
#else
     MSG_ERROR("Dynamical calculation with ieig2rf 5 only work with NETCDF support.")
     ABI_UNUSED(ncid)
#endif
     ABI_DEALLOCATE(gkk)
     ABI_DEALLOCATE(eig2nkq_tmp)
   end if
!  print _EIGI2D file for this perturbation
   if (ieig2rf /= 5 ) then
     if(smdelta>0) then
       unitout = dtfil%unddb
       dscrpt=' Note : temporary (transfer) database '

       call ddb_hdr_init(ddb_hdr,dtset,psps,pawtab,DDB_VERSION,dscrpt,&
&       1,xred=xred,occ=occ_rbz)

       call ddb_hdr_open_write(ddb_hdr, dtfil%fnameabo_eigi2d, unitout)

       call ddb_hdr_free(ddb_hdr)

       call outbsd(bdeigrf,dtset,eigbrd,dtset%natom,nkpt_rbz,unitout)

!      Output of the EIGI2D.nc file.
       fname = strcat(dtfil%filnam_ds(4),"_EIGI2D.nc")
!      Broadening EIGI2D (real and Im)
       call eigr2d_init(eigbrd,eigi2d,dtset%mband,hdr0%nsppol,nkpt_rbz,dtset%natom)
#ifdef HAVE_NETCDF
       NCF_CHECK_MSG(nctk_open_create(ncid, fname, xmpi_comm_self), "Creating EIGI2D file")
       NCF_CHECK(crystal%ncwrite(ncid))
       NCF_CHECK(ebands_ncwrite(Bands, ncid))
       call eigr2d_ncwrite(eigi2d,dtset%qptn(:),dtset%wtq,ncid)
       NCF_CHECK(nf90_close(ncid))
#else
       ABI_UNUSED(ncid)
#endif
     end if !smdelta
   end if
 end if

 if (allocated(fan)) then
   ABI_DEALLOCATE(fan)
 end if
 if (allocated(eig2nkq_tmp)) then
   ABI_DEALLOCATE(eig2nkq_tmp)
 end if
 if (allocated(gkk)) then
   ABI_DEALLOCATE(gkk)
 end if

 call crystal%free()
 call ebands_free(Bands)
 call eigr2d_free(eigr2d)
 call eigr2d_free(eigi2d)
 call fan_free(fan2d)
 call gkk_free(gkk2d)


 call timab(148,2,tsec)
!DEBUG
!write(std_out,*)' eig2tot: exit'
!ENDDEBUG

end subroutine eig2tot
!!***

!!****f* m_eig2d/outbsd
!! NAME
!! outbsd
!!
!! FUNCTION
!! output bsd file for one perturbation (used for elphon calculations in anaddb)
!!
!! INPUTS
!!  bdeigrf=number of bands for which the derivatives of the eigenvalues have been computed
!!  dtset = dataset variable for run flags
!!  eig2nkq= second ordre eigenvalue (or electron lifetime) that must be printed out
!!  mpert= maximum number of perturbations
!!  nkpt_rbz= number of k-points for perturbation
!!  unitout= writting unit of file
!!
!! OUTPUTS
!!  to file
!!
!! PARENTS
!!      dfpt_looppert,eig2tot
!!
!! CHILDREN
!!
!! SOURCE

subroutine outbsd(bdeigrf,dtset,eig2nkq,mpert,nkpt_rbz,unitout)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bdeigrf,mpert,nkpt_rbz,unitout
 type(dataset_type),intent(in) :: dtset
!arrays
 real(dp),intent(in) :: eig2nkq(2,dtset%mband*dtset%nsppol,nkpt_rbz,3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer :: bandtot_index,iband,idir1,idir2,ikpt,ipert1,ipert2,isppol

! *********************************************************************

!DEBUG
!write(std_out,*)' outbsd : enter'
!write(std_out,*)' eig2nkq(1,1,1,1,1,1)=',eig2nkq(1,1,1,1,1,1,1)
!ENDDEBUG

!output information in this file
 write(unitout,*)
 write(unitout,'(a,i8)') ' 2nd eigenvalue derivatives   - # elements :', 9*dtset%natom**2
 write(unitout,'(a,3es16.8,a)') ' qpt', dtset%qptn(:), ' 1.0'

!output RF eigenvalues

 do ikpt=1,nkpt_rbz
!  bandtot_index differs from zero only in the spin-polarized case
   bandtot_index=0
   write (unitout,'(a,3es16.8)') ' K-point:', dtset%kptns(:,ikpt)
   do isppol=1,dtset%nsppol
     do iband=1,bdeigrf
       write (unitout,'(a,i5)') ' Band:', iband+bandtot_index
!      write (unitout,*) 'ipert1     ','idir1     ','ipert2     ','idir2    ','Real    ','Im    '
       do ipert2=1,mpert
         do idir2=1,3
           do ipert1=1,mpert
             do idir1=1,3
               write (unitout,'(4i4,2d22.14)') idir1,ipert1,idir2,ipert2,&
&               eig2nkq(1,iband+bandtot_index,ikpt,idir1,ipert1,idir2,ipert2),&
&               eig2nkq(2,iband+bandtot_index,ikpt,idir1,ipert1,idir2,ipert2)
             end do !idir2
           end do !ipert2
         end do !idir1
       end do !ipert1
     end do !iband
     bandtot_index = bandtot_index + dtset%mband
   end do !isppol
 end do !ikpt

!close bsd file
 close (unitout)

end subroutine outbsd
!!***

!!****f* m-eig2d/smeared_delta
!! NAME
!! smeared_delta
!!
!! FUNCTION
!! This subroutine calculates the smeared delta that weights matrix elements:
!! \delta (\epsilon_{kn}-\epsilon_{k+Q,n'})
!!
!! INPUTS
!! eigen0(mband*nsppol) : Eigenvalues at point K
!! eigenq(mband*nsppol)  : Eigenvalues at point K+Q
!! mband : maximum number of bands
!! smdelta : Variable controlling the smearinf scheme
!!
!! OUTPUT
!! smdfunc(mband,mband) : Smeared delta function weight corresponding to \delta(\epsilon_{n,k} - \epsilon_{n',k+Q})
!!
!! PARENTS
!!      eig2stern,eig2tot
!!
!! CHILDREN
!!
!! SOURCE

subroutine smeared_delta(eigen0,eigenq,esmear,mband,smdelta,smdfunc)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,smdelta
!arrays
 real(dp),intent(in) :: eigen0(mband),eigenq(mband),esmear
 real(dp),intent(out) :: smdfunc(mband,mband)

!Local variables-------------------------------
!tolerance for non degenerated levels
!scalars
 integer :: ii,jj
 real(dp) :: aa,dsqrpi,gauss,xx
 character(len=500) :: message

! *********************************************************************


!---------------------------------------------------------
!Ordinary (unique) smearing function
!---------------------------------------------------------

 if(smdelta==1)then

!  Fermi-Dirac
   do ii=1,mband
     do jj= 1,mband
       xx= ( eigen0(ii) - eigenq(jj) )/esmear
       smdfunc(ii,jj)=0.25_dp/esmear/(cosh(xx/2.0_dp))**2
     end do
   end do

 else if(smdelta==2 .or. smdelta==3)then

!  Cold smearing of Marzari, two values of the "a" parameter being possible
!  first value gives minimization of the bump
   if(smdelta==2)aa=-.5634
!  second value gives monotonic occupation function
   if(smdelta==3)aa=-.8165

   dsqrpi=1.0_dp/sqrt(pi)
   do ii=1,mband
     do jj=1,mband
       xx= ( eigen0(ii) - eigenq(jj) ) / esmear
       gauss=dsqrpi*exp(-xx**2)/esmear
       smdfunc(ii,jj)=gauss*(1.5_dp+xx*(-aa*1.5_dp+xx*(-1.0_dp+aa*xx)))
     end do
   end do

 else if(smdelta==4)then

!  First order Hermite-Gaussian of Paxton and Methfessel
   dsqrpi=1.0_dp/sqrt(pi)
   do ii=1,mband
     do jj=1,mband
       xx= ( eigen0(ii) - eigenq (jj) ) / esmear
       smdfunc(ii,jj)=dsqrpi*(1.5_dp-xx**2)*exp(-xx**2)/esmear
     end do
   end do

 else if(smdelta==5)then

!  Gaussian smearing
   dsqrpi=1.0_dp/sqrt(pi)
   do ii=1,mband
     do jj=1,mband
       xx= ( eigen0(ii) - eigenq (jj) ) / esmear
       smdfunc(ii,jj)=dsqrpi*exp(-xx**2)/esmear
     end do
   end do

 else
   write(message, '(a,i0,a)' )'  Smdelta= ',smdelta,' is not allowed in smdfunc'
   MSG_BUG(message)
 end if

end subroutine smeared_delta
!!***

!!****f* m_eig2d/elph2_fanddw
!! NAME
!! elph2_fanddw
!!
!! FUNCTION
!! This routine calculates the zero-point motion corrections
!! due to the Fan term or to the DDW term..
!!
!! INPUTS
!!  dim_eig2nkq=1 if eig2nkq is to be computed
!!  displ(2*3*natom*3*natom)=the displacements of atoms in cartesian coordinates.
!!  eig2nkq(2,mband*nsppol,nkpt,3,natom,3,natom*dim_eig2nkq)=one half second derivatives of the electronic eigenvalues
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!  mband= maximum number of bands
!!  natom= number of atoms in the unit cell
!!  nkpt= number of k-points
!!  nsppol= 1 for unpolarized, 2 for spin-polarized
!!  option 1 for Fan term, 2 for DDW term
!!  phfrq(3*natom)=phonon frequencies
!!  (prtvol > 4) if the mode decomposition is to be printed
!!
!! OUTPUT
!!  eigen_corr(mband*nkpt*nsppol)= T=0 correction to the electronic eigenvalues, due to the Fan term.
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine elph2_fanddw(dim_eig2nkq,displ,eig2nkq,eigen_corr,gprimd,mband,natom,nkpt,nsppol,option,phfrq,prtvol)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dim_eig2nkq,mband,natom,nkpt,nsppol,option,prtvol

!arrays
 real(dp) :: gprimd(3,3)
 real(dp),intent(in) :: displ(2*3*natom*3*natom)
 real(dp),intent(in) :: eig2nkq(2,mband*nsppol,nkpt,3,natom,3,natom*dim_eig2nkq)
 real(dp),intent(in) :: phfrq(3*natom)
 real(dp),intent(out) :: eigen_corr(mband*nkpt*nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: neigs_per_line=6
 integer :: iatom1,iatom2,idir1,idir2,iband,ikpt,imode,index,isppol, imin, ii
 real(dp) :: d_at1_dir1_re,d_at1_dir1_im
 real(dp) :: d_at1_dir2_re,d_at1_dir2_im
 real(dp) :: d_at2_dir1_re,d_at2_dir1_im
 real(dp) :: d_at2_dir2_re,d_at2_dir2_im
 real(dp) :: e2_im,e2_re
 real(dp), allocatable :: eigen_corr_mode(:)
 character(len=500) :: message
 character(len=20) :: eig_format, line_format
!arrays
 real(dp) :: displ2cart(2,3,3),displ2red(2,3,3),tmp_displ2(2,3,3)

! *********************************************************************

!DEBUG
!write(std_out,*)' elph2_fanddw : enter '
!write(std_out,*)' option=',option
!ENDDEBUG

 if(option/=1 .and. option/=2)then
   write(message,'(a,i0)')' The argument option should be 1 or 2, while it is found that option=',option
   MSG_BUG(message)
 end if

 !printing options
 eig_format='f16.8'
 write(line_format,'(a,i1,a6,a)') '(',neigs_per_line,eig_format,')'

 if (prtvol > 4) then
   write(message,'(a,a)')ch10,' ================================================================================'
   call wrtout(ab_out_default,message,'COLL')
   if (option==1) then
     write(message,'(a)') ' ---- Begin Fan contributions to eigenvalues renormalization by mode ----'
     call wrtout(ab_out_default,message,'COLL')
   else if (option==2) then
     write(message,'(a)') ' ---- Begin DDW contributions to eigenvalues renormalization by mode ----'
     call wrtout(ab_out_default,message,'COLL')
   end if
 end if

 ABI_ALLOCATE(eigen_corr_mode,(mband*nkpt*nsppol))

 eigen_corr(:)=zero
 do imode=1,3*natom
   eigen_corr_mode(:)=zero

!  DEBUG
!  write(std_out,*)' Contribution of mode ',imode,' with frequency=',phfrq(imode),' and displacements :'
!  write(std_out,'(2f14.7)' ) displ(1+2*3*natom*(imode-1):2*3*natom*imode)
!  ENDDEBUG

   if (phfrq(imode)>tol6) then
     do iatom1=1,natom
       do iatom2=1,natom
!        DEBUG
!        write(std_out,*)' iatom1,iatom2=',iatom1,iatom2
!        ENDDEBUG

         do idir1=1,3
           do idir2=1,3
!            Compute the mean cartesian displacements
             d_at1_dir1_re=displ(1 + 2*(idir1-1 +3*(iatom1-1 +natom*(imode-1))))
             d_at1_dir1_im=displ(2 + 2*(idir1-1 +3*(iatom1-1 +natom*(imode-1))))
             d_at2_dir2_re=displ(1 + 2*(idir2-1 +3*(iatom2-1 +natom*(imode-1))))
             d_at2_dir2_im=displ(2 + 2*(idir2-1 +3*(iatom2-1 +natom*(imode-1))))

!            DEBUG
!            write(std_out,*)' idir1,idir2=',iatom1,iatom2,idir1,idir2
!            write(std_out,'(a,4f12.5)' )' d_at1_dir1 re,d_at2_dir2 re=',d_at1_dir1_re,d_at2_dir2_re
!            ENDDEBUG

             if(option==1)then
!              Compute the mean displacement correlation at T=0.
!              Consistent with Eqs.(7) and (8) of PRB51, 8610 (1995) [[cite:Lee1995]], specialized for the contribution of one q point.
!              but generalized to two different atoms. Note that the complex conjugate is taken on the second direction.
               displ2cart(1,idir1,idir2)=(d_at1_dir1_re*d_at2_dir2_re+ &
&               d_at1_dir1_im*d_at2_dir2_im )/(two*phfrq(imode))
               displ2cart(2,idir1,idir2)=(d_at1_dir1_im*d_at2_dir2_re- &
&               d_at1_dir1_re*d_at2_dir2_im )/(two*phfrq(imode))
             else if(option==2)then
!              Compute the mean square displacement correlation of each atom at T=0, and take mean over iatom1 and iatom2.
!              See Eqs.(7) and (8) of PRB51, 8610 (1995) [[cite:Lee1995]], specialized for the contribution of one q point.
!              Note that the complex conjugate is taken on the second direction.
!              Also, note the overall negative sign, to make it opposite to the Fan term.
               d_at1_dir2_re=displ(1 + 2*(idir2-1 +3*(iatom1-1 +natom*(imode-1))))
               d_at1_dir2_im=displ(2 + 2*(idir2-1 +3*(iatom1-1 +natom*(imode-1))))
               d_at2_dir1_re=displ(1 + 2*(idir1-1 +3*(iatom2-1 +natom*(imode-1))))
               d_at2_dir1_im=displ(2 + 2*(idir1-1 +3*(iatom2-1 +natom*(imode-1))))
               displ2cart(1,idir1,idir2)=-(d_at1_dir1_re*d_at1_dir2_re+ &
&               d_at1_dir1_im*d_at1_dir2_im+ &
&               d_at2_dir1_re*d_at2_dir2_re+ &
&               d_at2_dir1_im*d_at2_dir2_im )/(four*phfrq(imode))
               displ2cart(2,idir1,idir2)=-(d_at1_dir1_im*d_at1_dir2_re- &
&               d_at1_dir1_re*d_at1_dir2_im+ &
&               d_at2_dir1_im*d_at2_dir2_re- &
&               d_at2_dir1_re*d_at2_dir2_im )/(four*phfrq(imode))
             end if
           end do
         end do
!        Switch to reduced coordinates in two steps
         tmp_displ2(:,:,:)=zero
         do idir1=1,3
           do idir2=1,3
             tmp_displ2(:,:,idir1)=tmp_displ2(:,:,idir1)+displ2cart(:,:,idir2)*gprimd(idir2,idir1)
           end do
         end do
         displ2red(:,:,:)=zero
         do idir1=1,3
           do idir2=1,3
             displ2red(:,idir1,:)=displ2red(:,idir1,:)+tmp_displ2(:,idir2,:)*gprimd(idir2,idir1)
           end do
         end do
!        Compute the T=0 shift due to this q point
         do idir1=1,3
           do idir2=1,3
             do ikpt=1,nkpt
               do isppol=1,nsppol
                 do iband=1,mband
                   index=iband+mband*(isppol-1 + nsppol*(ikpt-1))
                   e2_re=eig2nkq(1,iband+mband*(isppol-1),ikpt,idir1,iatom1,idir2,iatom2)
                   e2_im=eig2nkq(2,iband+mband*(isppol-1),ikpt,idir1,iatom1,idir2,iatom2)
                   eigen_corr(index)=eigen_corr(index)+&
&                   e2_re*displ2red(1,idir1,idir2)-e2_im*displ2red(2,idir1,idir2)
                   eigen_corr_mode(index)=eigen_corr_mode(index)+&
&                   e2_re*displ2red(1,idir1,idir2)-e2_im*displ2red(2,idir1,idir2)
                 end do  ! band
               end do  ! spin
             end do  ! kpt
           end do  ! dir2
         end do  ! dir1
       end do  ! atom2
     end do  ! atom1
   end if

   if (prtvol > 4) then
     ! Print the corrections by mode
     write(message,'(a,i1)') ' imode= ',imode
     call wrtout(ab_out_default,message,'COLL')

     do ikpt=1,nkpt
       do isppol=1,nsppol
         write(message,'(a,i4,a,i1)')' ikpt= ',ikpt,' ispin= ',isppol
         call wrtout(ab_out_default,message,'COLL')

         imin = mband * (isppol-1 + nsppol*(ikpt-1))
         do ii=0, (mband-1)/neigs_per_line
           write(message, line_format) (eigen_corr_mode(iband+imin), &
&           iband = 1 + ii * neigs_per_line, min(mband, (ii+1)*neigs_per_line))
           call wrtout(ab_out_default,message,'COLL')
         end do
       end do
     end do
   end if

 end do  ! mode

 if (prtvol > 4) then
   if (option==1) then
     write(message,'(a)') ' ---- End Fan contribution to eigenvalues renormalization by mode ----'
     call wrtout(ab_out_default,message,'COLL')
   else if (option==2) then
     write(message,'(a)') ' ---- End DDW contribution to eigenvalues renormalization by mode ----'
     call wrtout(ab_out_default,message,'COLL')
   end if
   write(message,'(a,a)')' ================================================================================', ch10
   call wrtout(ab_out_default,message,'COLL')
 end if

 ABI_DEALLOCATE(eigen_corr_mode)

!DEBUG
!write(std_out,*)' elph2_fanddw : exit'
!write(std_out,*)' eigen_corr(1)=',eigen_corr(1)
!ENDDEBUG

end subroutine elph2_fanddw
!!***

END MODULE m_eig2d
!!***
