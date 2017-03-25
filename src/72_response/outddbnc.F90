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
!! Copyright (C) 1999-2017 ABINIT group (XG,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  hdr0 <type(hdr_type)>=contains header info
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  natom = number of atoms in the unit cell.
!!  mpert = maximum number of perturbations.
!!  rprimd(3,3)=primitive vectors
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  qpt(3)=curret q-point, in reduced coordinates.
!!  d2matr(2,3,mpert,3,mpert)=second-order derivative of the total energy
!!  blkflg(3,mpert,3,mpert)=mask telling whether each element is computed (1) or not (0).
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      crystal_free,crystal_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine outddbnc (dtfil, dtset, hdr0, psps, natom, mpert, rprimd, xred, qpt, d2matr, blkflg)

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

 use m_fstrings,   only : strcat
 use m_crystal_io, only : crystal_ncwrite

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'outddbnc'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert, natom
 real(dp),intent(in) :: rprimd(3,3), xred(3,natom)
 type(dataset_type),intent(in) :: dtset
 type(datafiles_type), intent(in) :: dtfil
 type(hdr_type),intent(in) :: hdr0
 type(pseudopotential_type),intent(in) :: psps
 integer,intent(in) :: blkflg(3,mpert,3,mpert)
 real(dp),intent(in) :: d2matr(2,3,mpert,3,mpert)
 real(dp),intent(in) :: qpt(3)

!Local variables -------------------------
!scalars
 character(len=fnlen) :: fname
 integer :: ntypat
 integer :: ncid, ncerr
 integer :: cplex, cart_dir, one_dim
 integer :: ipert1, ipert2, idir1, idir2, ii
 type(crystal_t) :: Crystal
 real(dp), allocatable :: amu(:)
 real(dp) :: dynmat(2,3,natom,3,natom)
 integer :: dynmat_mask(3,natom,3,natom)
 real(dp) :: born_effective_charge_tensor(3,natom,3)
 integer :: born_effective_charge_tensor_mask(3,natom,3)

! *********************************************************************

#ifdef HAVE_NETCDF

 ! Initialize NetCDF file.
 fname = strcat(dtfil%filnam_ds(4),"_DDB.nc")
 NCF_CHECK(nctk_open_create(ncid, fname, xmpi_comm_self))

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

 ! Retrieve atomic masses
 ntypat = psps%ntypat
 ABI_ALLOCATE(amu,(ntypat))
 do ii=1,ntypat
   amu(ii) = dtset%amu_orig(ii,1)
 end do


 ! ------------------------------
 ! Write crystal info
 ! ------------------------------
 call crystal_init(dtset%amu_orig(:,1),Crystal, dtset%spgroup, dtset%natom, dtset%npsp, psps%ntypat, &
& dtset%nsym, rprimd, dtset%typat, xred, dtset%ziontypat, dtset%znucl, 1, &
& dtset%nspden==2.and.dtset%nsppol==1, .false., hdr0%title, &
& dtset%symrel, dtset%tnons, dtset%symafm)
 ncerr = crystal_ncwrite(Crystal, ncid)
 NCF_CHECK(ncerr)

 call crystal_free(Crystal)


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
 NCF_CHECK(nf90_put_var(ncid, vid('atomic_masses_amu'), amu))
 NCF_CHECK(nf90_put_var(ncid, vid('q_point_reduced_coord'), qpt))
 NCF_CHECK(nf90_put_var(ncid, vid('second_derivative_of_energy'), dynmat))
 NCF_CHECK(nf90_put_var(ncid, vid('second_derivative_of_energy_mask'), dynmat_mask))
 NCF_CHECK(nf90_put_var(ncid, vid('born_effective_charge_tensor'), born_effective_charge_tensor))
 NCF_CHECK(nf90_put_var(ncid, vid('born_effective_charge_tensor_mask'), born_effective_charge_tensor_mask))

 ! Close file
 NCF_CHECK(nf90_close(ncid))

 ! Deallocate stuff
 ABI_DEALLOCATE(amu)

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
