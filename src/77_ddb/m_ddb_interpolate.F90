!!****m* ABINIT/m_ddb_interpolate
!! NAME
!!  m_ddb_interpolate
!!
!! FUNCTION
!! Interpolate the ddb onto a fine set of q-points using
!! the interatomic force constants and write the result in a DDB file.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (GA)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_ddb_interpolate

 use defs_basis
 use m_errors
 use m_xmpi
 use m_abicore
 use m_ddb
 use m_ddb_hdr
 use m_ifc
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_anaddb_dataset, only : anaddb_dataset_type
 use m_crystal,        only : crystal_t
 use m_io_tools,       only : get_unit
 use m_fstrings,       only : strcat
 use m_dynmat,         only : gtdyn9, d2cart_to_red

 implicit none

 private
!!***

 public :: ddb_interpolate
 public :: outddbnc
!!***

contains
!!***

!!****f* ABINIT/ddb_interpolate
!!
!! NAME
!! ddb_interpolate
!!
!! FUNCTION
!! Interpolate the ddb onto a fine set of q-points using
!! the interatomic force constants and write the result in a DDB file.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      d2cart_to_red,ddb_free,ddb_hdr_open_write,ddb_malloc,ddb_write_block
!!      gtblk9,gtdyn9,outddbnc,wrtout
!!
!! SOURCE

subroutine ddb_interpolate(ifc, crystal, inp, ddb, ddb_hdr, asrq0, prefix, comm)

!Arguments -------------------------------
!scalars
 type(ifc_type),intent(in) :: ifc
 type(anaddb_dataset_type),target,intent(inout) :: inp
 type(crystal_t),intent(in) :: crystal
 type(ddb_type),intent(inout) :: ddb
 type(ddb_hdr_type),intent(inout) :: ddb_hdr
 type(asrq0_t),intent(inout) :: asrq0
 integer,intent(in) :: comm
 character(len=*),intent(in) :: prefix
!arrays

!Local variables -------------------------
!scalars
 integer,parameter :: master=0
 integer :: unddb
 integer :: nsym,natom,ntypat,mband,nqpt_fine
 integer :: msize,nsize,mpert,nblok,mtyp
 integer :: rftyp,choice
 integer :: ii,iblok,jblok,iqpt,ipert1,ipert2,idir1,idir2
 integer :: nprocs,my_rank
 character(len=500) :: msg
 character(len=fnlen) :: ddb_out_filename, ddb_out_nc_filename
 type(ddb_type) :: ddb_new
!arrays
 integer :: rfphon(4),rfelfd(4),rfstrs(4)
 integer,allocatable :: blkflg(:,:,:,:)
 real(dp) :: qpt(3), qptnrm(3), qpt_padded(3,3), qred(3)
 real(dp),allocatable :: d2cart(:,:,:,:,:),d2red(:,:,:,:,:)
 real(dp),pointer :: qpt_fine(:,:)

! *********************************************************************


 ! Only master works for the time being
 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 if (my_rank /= master) return

 ! ===================================
 ! Copy dimensions and allocate arrays
 ! ===================================

 write(msg, '(a,(80a),a,a,a,a)' ) ch10,('=',ii=1,80),ch10,ch10,' Treat the first list of vectors ',ch10
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

 nullify(qpt_fine)
 nqpt_fine = inp%nph1l
 qpt_fine => inp%qph1l

 rftyp=inp%rfmeth

 nsym = Crystal%nsym
 natom = Crystal%natom
 ntypat = Crystal%ntypat

 mband = ddb_hdr%mband

 mtyp = max(ddb_hdr%mblktyp, 2)  ! Limited to 2nd derivatives of total energy
 ddb_hdr%mblktyp = mtyp

 mpert = natom + 6
 msize = 3 * mpert * 3 * mpert  !; if (mtyp==3) msize=msize*3*mpert
 nsize = 3 * mpert * 3 * mpert
 nblok = nqpt_fine

 ddb_new%nblok = nblok
 call ddb_new%malloc(msize,nblok,natom,ntypat)
 ddb_new%flg = 0
 ddb_new%amu = ddb%amu
 ddb_new%typ = 1
 ddb_new%qpt = zero
 ddb_new%nrm = one

 ABI_MALLOC(d2cart,(2,3,mpert,3,mpert))
 ABI_MALLOC(d2red,(2,3,mpert,3,mpert))
 ABI_MALLOC(blkflg,(3,mpert,3,mpert))

 blkflg = 1

 rfphon(1:2)=1; rfelfd(1:2)=0; rfstrs(1:2)=0
 qpt_padded = zero

 ddb_out_filename = strcat(prefix, "_DDB")

 ddb_hdr%dscrpt = 'Interpolated DDB using interatomic force constants'
 ddb_hdr%nblok = nblok

 ! ================================================
 ! Interpolate the dynamical matrix at each q-point
 ! ================================================

 qptnrm = one

 do iqpt=1,nqpt_fine

   ! Initialisation of the phonon wavevector
   qpt(:)=qpt_fine(:,iqpt)

   if (inp%nph1l /= 0) qptnrm(1) = inp%qnrml1(iqpt)

   ! Look for the information in the DDB
   qpt_padded(:,1) = qpt
   call ddb%get_block(iblok,qpt_padded,qptnrm,rfphon,rfelfd,rfstrs,rftyp)

   if (iblok /= 0) then
     ! DEBUG
     !write(*,*) 'DDB found in file for qpt=', qpt
     ! END DEBUG

     ! q-point is present in the ddb. No interpolation needed.

     !d2cart(:,1:msize) = ddb%val(:,:,iblok)
     d2cart(1,:,:,:,:) = reshape(ddb%val(1,:,iblok), shape = (/3,mpert,3,mpert/))
     d2cart(2,:,:,:,:) = reshape(ddb%val(2,:,iblok), shape = (/3,mpert,3,mpert/))
   else

     ! Get d2cart using the interatomic forces and the
     ! long-range coulomb interaction through Ewald summation
     call gtdyn9(ddb%acell,Ifc%atmfrc,Ifc%dielt,Ifc%dipdip,Ifc%dyewq0,d2cart, &
&     crystal%gmet,ddb%gprim,ddb%mpert,natom,Ifc%nrpt,qptnrm(1), &
&     qpt, crystal%rmet,ddb%rprim,Ifc%rpt,Ifc%trans,crystal%ucvol, &
&     Ifc%wghatm,crystal%xred,ifc%zeff, xmpi_comm_self)

   end if

   ! Eventually impose the acoustic sum rule based on previously calculated d2asr
   call asrq0%apply(natom, ddb%mpert, ddb%msize, crystal%xcart, d2cart)

   ! Transform d2cart into reduced coordinates.
   call d2cart_to_red(d2cart,d2red,crystal%gprimd,crystal%rprimd,mpert, &
&   natom,ntypat,crystal%typat,crystal%ucvol,crystal%zion)

   ! Store the dynamical matrix into a block of the new ddb
   jblok = iqpt
   ddb_new%val(1,1:nsize,jblok) = reshape(d2red(1,:,:,:,:), shape = (/nsize/))
   ddb_new%val(2,1:nsize,jblok) = reshape(d2red(2,:,:,:,:), shape = (/nsize/))

   ! Store the q-point
   ddb_new%qpt(1:3,jblok) = qpt
   ddb_new%nrm(1,jblok) = qptnrm(1)

   ! Set the flags
   ii=0
   do ipert2=1,mpert
     do idir2=1,3
       do ipert1=1,mpert
         do idir1=1,3
           ii=ii+1
           if (ipert1<=natom.and.ipert2<=natom) then
             ddb_new%flg(ii,jblok) = 1
           end if
         end do
       end do
     end do
   end do

 end do ! iqpt

 ! Copy the flags for Gamma
 qpt_padded(:,1) = zero
 qptnrm = one
 call ddb%get_block(iblok,qpt_padded,qptnrm,rfphon,rfelfd,rfstrs,rftyp)
 call ddb_new%get_block(jblok,qpt_padded,qptnrm,rfphon,rfelfd,rfstrs,rftyp)

 if (iblok /= 0 .and. jblok /= 0) then

   ii=0
   do ipert2=1,mpert
     do idir2=1,3
       do ipert1=1,mpert
         do idir1=1,3
           ii=ii+1
           ddb_new%flg(ii,jblok) = ddb%flg(ii,iblok)
         end do
       end do
     end do
   end do

 end if


 ! =======================
 ! Write the new DDB files
 ! =======================

 if (my_rank == master) then

   unddb = get_unit()

   ! Write the DDB header
   call ddb_hdr_open_write(ddb_hdr, ddb_out_filename, unddb)

   ! Write the whole database
   call wrtout(std_out,' write the DDB ','COLL')
   choice=2
   do iblok=1,nblok
     call ddb_new%write_block(iblok,choice,ddb_hdr%mband,mpert,msize,ddb_hdr%nkpt,unddb)
   end do

   ! Also write summary of bloks at the end
   write(unddb, '(/,a)' )' List of bloks and their characteristics '
   choice=3
   do iblok=1,nblok
     call ddb_new%write_block(iblok,choice,ddb_hdr%mband,mpert,msize,ddb_hdr%nkpt,unddb)
   end do

   close (unddb)

#ifdef HAVE_NETCDF
   ! Write the DDB.nc files (one for each q-point)
   do jblok=1,nblok

     write(ddb_out_nc_filename,'(2a,i5.5,a)') trim(prefix),'_qpt_',jblok,'_DDB.nc'

     d2red(1,:,:,:,:) = reshape(ddb_new%val(1,1:nsize,jblok), &
&     shape = (/3,mpert,3,mpert/))
     d2red(2,:,:,:,:) = reshape(ddb_new%val(2,1:nsize,jblok), &
&     shape = (/3,mpert,3,mpert/))
     blkflg(:,:,:,:)  = reshape(ddb_new%flg(1:nsize,jblok), &
&     shape = (/3,mpert,3,mpert/))

     do ii=1,3
       qred(ii) = ddb_new%qpt(ii,jblok) / ddb_new%nrm(1,jblok)
     end do

     call outddbnc(ddb_out_nc_filename,mpert,d2red,blkflg,qred,crystal)

   end do
#endif
 end if

 ! ===========
 ! Free memory
 ! ===========

 call ddb_new%free()
 ABI_FREE(d2cart)
 ABI_FREE(d2red)
 ABI_FREE(blkflg)

end subroutine ddb_interpolate
!!***


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

subroutine outddbnc (filename, mpert, d2matr, blkflg, qpt, Crystal)

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
 NCF_CHECK(crystal%ncwrite(ncid))

 ! ------------------------------
 ! Write DDB
 ! ------------------------------

 ! Write the dimensions specified by ETSF
 one_dim = 1
 cplex = 2
 cart_dir = 3

 ncerr = nctk_def_dims(ncid, [&
  nctkdim_t('current_one_dim', one_dim), &
  nctkdim_t('number_of_atoms', natom), &
  nctkdim_t('number_of_cartesian_directions', cart_dir), &
  nctkdim_t('number_of_perturbations', mpert), &
  nctkdim_t('cplex',cplex)], defmode=.True.)
 NCF_CHECK(ncerr)

! Create the arrays
 ncerr = nctk_def_arrays(ncid, [&
 &nctkarr_t('atomic_masses_amu', "dp", 'number_of_atom_species'),&
 nctkarr_t('q_point_reduced_coord', "dp", 'number_of_cartesian_directions'),&
 nctkarr_t('second_derivative_of_energy', "dp", &
 &'cplex, number_of_cartesian_directions, number_of_atoms, number_of_cartesian_directions, number_of_atoms'),&
 nctkarr_t('second_derivative_of_energy_mask', "i",&
 &'number_of_cartesian_directions, number_of_atoms, number_of_cartesian_directions, number_of_atoms'),&
 nctkarr_t('born_effective_charge_tensor', "dp",'number_of_cartesian_directions,number_of_atoms,number_of_cartesian_directions'),&
 nctkarr_t('born_effective_charge_tensor_mask', "i",&
 &'number_of_cartesian_directions, number_of_atoms, number_of_cartesian_directions')])
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
   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
 end function vid

end subroutine outddbnc
!!***

end module m_ddb_interpolate
!!***
