!{\src2tex{textfont=tt}}
!!****f* ABINIT/outddbnc
!!
!! NAME
!! outddbnc
!!
!! FUNCTION
!! Write Derivative DataBase in NetCDF format.
!! See ~abinit/scripts/post_processing/merge_ddb_nc.py
!! for a merging utility.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group (XG,MT,GA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  filename=name of the DDB.nc file to be written.
!!  Crystal <type(crystal_t)>=info on the crystal
!!  natom = number of atoms in the unit cell.
!!  mpert = maximum number of perturbations.
!!  qpt(3)=curret q-point, in reduced coordinates.
!!  d2matr(2,3,mpert,3,mpert)=second-order derivative of the total energy
!!  blkflg(3,mpert,3,mpert)=mask telling whether each element is computed (1) or not (0).
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      ddb_interpolate,respfn
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine outddbnc (filename, mpert, d2matr, blkflg, qpt, Crystal)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_xmpi
 use m_profiling_abi
 use m_nctk
 use m_crystal
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_crystal_io, only : crystal_ncwrite

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'outddbnc'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 character(len=*),intent(in) :: filename
 integer,intent(in) :: mpert
 integer,intent(in) :: blkflg(3,mpert,3,mpert)
 real(dp),intent(in) :: d2matr(2,3,mpert,3,mpert)
 real(dp),intent(in) :: qpt(3)
 type(crystal_t), intent(in) :: Crystal

!Local variables -------------------------
!scalars
 integer :: natom
 integer :: ncid, ncerr
 integer :: cplex, cart_dir, one_dim
 integer :: ipert1, ipert2, idir1, idir2
 integer,allocatable :: dynmat_mask(:,:,:,:)
 integer,allocatable :: born_effective_charge_tensor_mask(:,:,:)
 real(dp),allocatable :: dynmat(:,:,:,:,:)
 real(dp),allocatable :: born_effective_charge_tensor(:,:,:)

! *********************************************************************

#ifdef HAVE_NETCDF

 natom = Crystal%natom

 ABI_ALLOCATE(dynmat, (2,3,natom,3,natom))
 ABI_ALLOCATE(dynmat_mask, (3,natom,3,natom))
 ABI_ALLOCATE(born_effective_charge_tensor, (3,natom,3))
 ABI_ALLOCATE(born_effective_charge_tensor_mask, (3,natom,3))

 ! Initialize NetCDF file.
 NCF_CHECK(nctk_open_create(ncid, filename, xmpi_comm_self))

 ! ------------------------------
 ! Construct local DDB
 ! ------------------------------

 ! Construct the dynamical matrix
 do ipert1=1,natom
   do idir1=1,3
     do ipert2=1,natom
       do idir2=1,3
         dynmat_mask(idir1,ipert1,idir2,ipert2) = blkflg(idir1,ipert1,idir2,ipert2)
         if (blkflg(idir1,ipert1,idir2,ipert2)==1) then
           dynmat(1,idir1,ipert1,idir2,ipert2) = d2matr(1,idir1,ipert1,idir2,ipert2)
           dynmat(2,idir1,ipert1,idir2,ipert2) = d2matr(2,idir1,ipert1,idir2,ipert2)
         else
           dynmat(1,idir1,ipert1,idir2,ipert2) = zero
           dynmat(2,idir1,ipert1,idir2,ipert2) = zero
         end if
       end do
     end do
   end do
 end do 

 ! Construct the Born effective charge tensor
 ipert2 = natom + 2
 do ipert1=1,natom
   do idir1=1,3
     do idir2=1,3
       born_effective_charge_tensor_mask(idir1,ipert1,idir2) = blkflg(idir1,ipert1,idir2,ipert2)
       if (blkflg(idir1,ipert1,idir2,ipert2)==1) then
            ! This is a real quantity
         born_effective_charge_tensor(idir1,ipert1,idir2) = d2matr(1,idir1,ipert1,idir2,ipert2)
       else
         born_effective_charge_tensor(idir1,ipert1,idir2) = zero
       end if
     end do
   end do
 end do 

 ! TODO: also store the dielectric matrix

 ! ------------------------------
 ! Write crystal info
 ! ------------------------------
 ncerr = crystal_ncwrite(Crystal, ncid)
 NCF_CHECK(ncerr)


 ! ------------------------------
 ! Write DDB
 ! ------------------------------

 ! Write the dimensions specified by ETSF
 one_dim = 1
 cplex = 2
 cart_dir = 3

 ncerr = nctk_def_dims(ncid, [&
& nctkdim_t('current_one_dim', one_dim), &
& nctkdim_t('number_of_atoms', natom), &
& nctkdim_t('number_of_cartesian_directions', cart_dir), &
& nctkdim_t('number_of_perturbations', mpert), &
& nctkdim_t('cplex',cplex)], defmode=.True.)
 NCF_CHECK(ncerr)

 ! Create the arrays
 ncerr = nctk_def_arrays(ncid, [&
 nctkarr_t('atomic_masses_amu', "dp", 'number_of_atom_species'),&
 nctkarr_t('q_point_reduced_coord', "dp", 'number_of_cartesian_directions'),&
 nctkarr_t('second_derivative_of_energy', "dp", 'cplex, &
& number_of_cartesian_directions, number_of_atoms, &
& number_of_cartesian_directions, number_of_atoms'), &
 nctkarr_t('second_derivative_of_energy_mask', "i", '&
& number_of_cartesian_directions, number_of_atoms, &
& number_of_cartesian_directions, number_of_atoms'), &
 nctkarr_t('born_effective_charge_tensor', "dp", '&
& number_of_cartesian_directions, number_of_atoms, &
& number_of_cartesian_directions'), &
 nctkarr_t('born_effective_charge_tensor_mask', "i", ' &
& number_of_cartesian_directions, number_of_atoms, &
& number_of_cartesian_directions')])
 NCF_CHECK(ncerr)

! Write data
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, vid('atomic_masses_amu'), Crystal%amu))
 NCF_CHECK(nf90_put_var(ncid, vid('q_point_reduced_coord'), qpt))
 NCF_CHECK(nf90_put_var(ncid, vid('second_derivative_of_energy'), dynmat))
 NCF_CHECK(nf90_put_var(ncid, vid('second_derivative_of_energy_mask'), dynmat_mask))
 NCF_CHECK(nf90_put_var(ncid, vid('born_effective_charge_tensor'), born_effective_charge_tensor))
 NCF_CHECK(nf90_put_var(ncid, vid('born_effective_charge_tensor_mask'), born_effective_charge_tensor_mask))

 ! Close file
 NCF_CHECK(nf90_close(ncid))

 ! Deallocate stuff

 ABI_FREE(dynmat)
 ABI_FREE(dynmat_mask)
 ABI_FREE(born_effective_charge_tensor)
 ABI_FREE(born_effective_charge_tensor_mask)

#else
 MSG_ERROR("NETCDF support required to write DDB.nc file.")
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

end subroutine outddbnc
!!***
