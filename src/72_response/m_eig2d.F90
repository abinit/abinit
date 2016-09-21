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
!! Copyright (C) 2014-2016 ABINIT group (SP)
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
 use m_errors
 use m_profiling_abi
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif

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

!!***

!!****t* m_eig2d/eigr2d_t
!! NAME
!! eig2d_t
!!
!! FUNCTION
!! It contains informations about the second-order derivative of the
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
!! It contains informations about the second-order derivative of the 
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
!! It contains informations about the second-order derivative of the 
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
!=====================================================================================
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eigr2d_init'
!End of the abilint section

 implicit none

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eigr2d_init'
!End of the abilint section

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eigr2d_ncwrite'
!End of the abilint section

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

 ncerr = nctk_def_arrays(ncid, [&
&  nctkarr_t('current_q_point', "dp", 'number_of_cartesian_directions'),&
&  nctkarr_t('current_q_point_weight', "dp", 'current_one_dim'),&
&  nctkarr_t('second_derivative_eigenenergies', "dp", &
&    'cplex, product_mband_nsppol, number_of_kpoints, number_of_cartesian_directions, number_of_atoms,&
& number_of_cartesian_directions, number_of_atoms')])
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vid'
!End of the abilint section

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eigr2d_free'
!End of the abilint section

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fan_init'
!End of the abilint section

 implicit none

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fan_init'
!End of the abilint section

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
!!      dfpt_looppert,eig2tot
!!
!! CHILDREN
!!
!! SOURCE

subroutine gkk_init(gkk,gkk2d,mband,nsppol,nkpt,natom,ncart)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gkk_init'
!End of the abilint section

 implicit none

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gkk_init'
!End of the abilint section

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fan_ncwrite'
!End of the abilint section

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

 ncerr = nctk_def_arrays(ncid, [&
   nctkarr_t('current_q_point', "dp", 'number_of_cartesian_directions'),&
   nctkarr_t('current_q_point_weight', "dp", 'current_one_dim'),&
   nctkarr_t('second_derivative_eigenenergies_actif', "dp", &
&'product_mband_nsppol2, number_of_kpoints, number_of_cartesian_directions, &
&number_of_atoms, number_of_cartesian_directions, number_of_atoms, max_number_of_states')])
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

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vid'
!End of the abilint section

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
!!      dfpt_looppert,eig2tot
!!
!! CHILDREN
!!
!! SOURCE

subroutine gkk_ncwrite(gkk2d,iqpt,wtq,ncid)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gkk_ncwrite'
!End of the abilint section

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vid'
!End of the abilint section

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fan_free'
!End of the abilint section

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
!!      dfpt_looppert,eig2tot
!!
!! CHILDREN
!!
!! SOURCE

subroutine gkk_free(gkk2d)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gkk_free'
!End of the abilint section

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


END MODULE m_eig2d
!!***

