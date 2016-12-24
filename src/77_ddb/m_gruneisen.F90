!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_gruneisen
!! NAME
!!  m_gruneisen
!!
!! FUNCTION
!!  This module contains the declaration of data types and methods
!!  to compute Gruneisen parameters.
!!
!! COPYRIGHT
!! Copyright (C) 2011-2016 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_gruneisen

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_xmpi
 use m_crystal
 use m_tetrahedron
 use m_ddb
 use m_ifc
 use m_cgtools
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_io_tools,            only : get_unit
 use m_time,                only : cwtime
 use m_fstrings,            only : sjoin, itoa, ltoa, ftoa
 use m_numeric_tools,       only : central_finite_diff
 use m_kpts,                only : kpts_ibz_from_kptrlatt, tetra_from_kptrlatt
 use m_bz_mesh,             only : kpath_t, kpath_init, kpath_free
 use m_anaddb_dataset,      only : anaddb_dataset_type
 use m_dynmat,              only : massmult_and_breaksym

 implicit none

 private
!!***

!!****t* m_gruneisen/gruns_t
!! NAME
!! gruns_t
!!
!! FUNCTION
!!  Contains the necessary data to interpolate the
!!  phonon bandstructure and eigenvectors in reciprocal space (ie.
!!  interatomic force constants and corresponding real space grid info).
!!
!! SOURCE

 type,public :: gruns_t

   integer :: natom3
    ! 3 * natom

   integer :: nvols
    ! Number of volumes.

   integer :: iv0
    ! Index of the DDB file associated to the equilibrium volume

   real(dp) :: v0
    ! Equilibrium volume.

   real(dp) :: delta_vol
    ! Uniform grid spacing for finite difference

   type(crystal_t),allocatable :: cryst_vol(:)
    ! cryst_vol(nvols)
    ! crystalline structure for the different volumes

   type(ddb_type),allocatable :: ddb_vol(:)
    ! dbb_vol(nvols)
    ! DDB objects for the different volumes

   type(ifc_type),allocatable :: ifc_vol(:)
    ! ifc_vol(nvols)
    ! interatomic force constants for the different volumes

 end type gruns_t

 public :: gruns_new        ! Constructor
 public :: gruns_qpath      ! Compute Grunesein parameters on a q-path
 public :: gruns_qmesh      ! Compute Grunesein parameters on a q-mesh.
 public :: gruns_free       ! Release memory
 public :: gruns_anaddb
!!***

contains  !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_gruneisen/gruns_new
!! NAME
!!  gruns_new
!!
!! FUNCTION
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(gruns_t) function gruns_new(ddb_paths, inp, comm) result(new)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_free'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: comm
 type(anaddb_dataset_type),intent(in) :: inp
!arrays
 character(len=*),intent(in) :: ddb_paths(:)

!Local variables-------------------------------
 integer,parameter :: natifc0=0,master=0
 integer :: ivol,iblock,msym,dimekb,lmnmax,mband,mblktyp,natom,nblok,nkpt,ntypat,usepaw,ddbun
 integer :: nprocs,my_rank,ierr
 character(len=500) :: msg
!arrays
 integer,allocatable :: atifc0(:)
 real(dp) :: dielt(3,3)
 real(dp),allocatable :: zeff(:,:,:)

! ************************************************************************

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 new%nvols = size(ddb_paths)
 ABI_MALLOC(new%cryst_vol, (new%nvols))
 ABI_MALLOC(new%ddb_vol, (new%nvols))
 ABI_MALLOC(new%ifc_vol, (new%nvols))

 call wrtout(ab_out, "Computation of Gruneisen parameter with central finite difference:")

 ddbun = get_unit()
 do ivol=1,new%nvols
   call wrtout(ab_out, sjoin(" Reading DDB file:", ddb_paths(ivol)))
   call ddb_getdims(dimekb,ddb_paths(ivol),lmnmax,mband,mblktyp,msym,natom,nblok,nkpt,ntypat,ddbun,usepaw,DDB_VERSION,comm)
   ABI_MALLOC(atifc0, (natom))
   atifc0 = 0
   call ddb_from_file(new%ddb_vol(ivol), ddb_paths(ivol), inp%brav, natom, natifc0, atifc0, new%cryst_vol(ivol), comm)
   ABI_FREE(atifc0)
   if (my_rank == master) then
     call crystal_print(new%cryst_vol(ivol), header=sjoin("Structure for ivol:", itoa(ivol)), unit=ab_out, prtvol=-1)
   end if

   ! Get Dielectric Tensor and Effective Charges
   ! (initialized to one_3D and zero if the derivatives are not available in the DDB file)
   ABI_MALLOC(zeff, (3,3,natom))
   iblock = ddb_get_dielt_zeff(new%ddb_vol(ivol), new%cryst_vol(ivol), inp%rfmeth, inp%chneut, inp%selectz, dielt, zeff)
   if (iblock == 0) then
     call wrtout(ab_out, sjoin("- Cannot find dielectric tensor and Born effective charges in DDB file:", ddb_paths(ivol)))
     call wrtout(ab_out, "Values initialized with zeros")
   else
     call wrtout(ab_out, sjoin("- Found dielectric tensor and Born effective charges in DDB file:", ddb_paths(ivol)))
   end if

   call ifc_init(new%ifc_vol(ivol), new%cryst_vol(ivol), new%ddb_vol(ivol),&
     inp%brav,inp%asr,inp%symdynmat,inp%dipdip,inp%rfmeth,inp%ngqpt(1:3),inp%nqshft,inp%q1shft,dielt,zeff,&
     inp%nsphere,inp%rifcsph,inp%prtsrlr,inp%enunit)
   ABI_FREE(zeff)
 end do

 ! Consistency check
 ! TODO
 ABI_CHECK(any(new%nvols == [3, 5, 7, 9]), "Central finite difference requires [3,5,7,9] DDB files")

 new%natom3 = 3 * new%cryst_vol(1)%natom
 new%iv0 = 1 + new%nvols / 2
 new%v0 = new%cryst_vol(new%iv0)%ucvol
 new%delta_vol = new%cryst_vol(new%iv0+1)%ucvol - new%cryst_vol(new%iv0)%ucvol

 ierr = 0
 do ivol=1,new%nvols
   !write(std_out,*)"ucvol",new%cryst_vol(ivol)%ucvol
   if (abs(new%cryst_vol(1)%ucvol + new%delta_vol * (ivol-1) - new%cryst_vol(ivol)%ucvol) > tol4) then
      write(std_out,*)"ucvol, delta_vol, diff_vol", new%cryst_vol(ivol)%ucvol, new%delta_vol, &
        abs(new%cryst_vol(1)%ucvol + new%delta_vol * (ivol-1) - new%cryst_vol(ivol)%ucvol)
      ierr = ierr + 1
   end if
 end do
 if (ierr /= 0) then
   msg = ltoa([(new%cryst_vol(ivol)%ucvol, ivol=1,new%nvols)])
   MSG_ERROR(sjoin("Gruneisen calculations requires linear mesh of volumes but received:", msg))
 end if

end function gruns_new
!!***

!----------------------------------------------------------------------

!!****f* m_gruneisen/gruns_fourq
!! NAME
!!  gruns_fourq
!!
!! FUNCTION
!!
!! INPUTS
!!  qpt(3)=q-point in reduced coordinates.
!!
!! OUTPUT
!!  phfrq(3*natom) = Gruneisen parameters
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine gruns_fourq(gruns, qpt, wv0, gvals)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(gruns_t),intent(in) :: gruns
!arrays
 real(dp),intent(in) :: qpt(3)
 real(dp),intent(out) :: wv0(gruns%natom3),gvals(gruns%natom3)

!Local variables-------------------------------
!scalars
 integer :: ivol,ii,natom3,nu
 real(dp) :: fact
!arrays
 real(dp) :: wvols(gruns%natom3,gruns%nvols),dot(2)
 real(dp) :: eigvec(2,gruns%natom3,gruns%natom3,gruns%nvols),d2cart(2,gruns%natom3,gruns%natom3,gruns%nvols)
 real(dp) :: dddv(2,gruns%natom3,gruns%natom3),displ_cart(2,gruns%natom3,gruns%natom3)
 real(dp) :: omat(2,gruns%natom3, gruns%natom3)

! ************************************************************************

 natom3 = gruns%natom3
 do ivol=1,gruns%nvols
   call ifc_fourq(gruns%ifc_vol(ivol), gruns%cryst_vol(ivol), qpt, wvols(:,ivol), displ_cart, &
                  out_d2cart=d2cart(:,:,:,ivol), out_eigvec=eigvec(:,:,:,ivol))

   call massmult_and_breaksym(gruns%cryst_vol(ivol)%natom, gruns%cryst_vol(ivol)%ntypat, &
     gruns%cryst_vol(ivol)%typat, gruns%ifc_vol(ivol)%amu, d2cart(:,:,:,ivol))

   !call zgemm('N','N',natom3,natom3,natom3,cone,d2cart(:,:,:,ivol),natom3,eigvec(:,:,:,ivol),natom3,czero,omat,natom3)
   !do nu=1,natom3
   !  write(std_out,*)"H|psi> - w**2 |psi>",maxval(abs(omat(:,:,nu) - wvols(nu,ivol) ** 2 * eigvec(:,:,nu,ivol)))
   !end do
 end do
 wv0 = wvols(:, gruns%iv0)

 ! Compute dD(q)/dV with finite difference.
 dddv = zero
 do ivol=1,gruns%nvols
   fact = central_finite_diff(1, ivol, gruns%nvols)
   if (fact /= zero) dddv = dddv + fact * d2cart(:,:,:,ivol)
 end do
 dddv = dddv / gruns%delta_vol

 ! Compute -volume/(2w(q)**2) <u(q)|dD(q)/dq|u(q)>
 call zgemm('N','N',natom3,natom3,natom3,cone,dddv,natom3,eigvec(:,:,:,gruns%iv0),natom3,czero,omat,natom3)
 do nu=1,natom3
   if (abs(wv0(nu)) > tol12) then
     dot = cg_zdotc(natom3, eigvec(1,1,nu, gruns%iv0), omat(1,1,nu))
     gvals(nu) = -gruns%v0 * dot(1) / (two * wv0(nu)**2)
     write(std_out,*)gvals(nu)
   else
     gvals(nu) = zero
   end if
 end do

end subroutine gruns_fourq
!!***

!----------------------------------------------------------------------

!!****f* m_gruneisen/gruns_qpath
!! NAME
!!  gruns_qpath
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine gruns_qpath(gruns, qpath, comm)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 type(gruns_t),intent(in) :: gruns
 type(kpath_t),intent(in) :: qpath

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: nprocs,my_rank,iqpt,ierr
!arrays
 real(dp),allocatable :: gvals(:,:),wv0(:,:)

! ************************************************************************

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 ABI_CALLOC(wv0, (gruns%natom3, qpath%npts))
 ABI_CALLOC(gvals, (gruns%natom3, qpath%npts))

 do iqpt=1,qpath%npts
   if (mod(iqpt, nprocs) /= my_rank) cycle ! mpi-parallelism
   call gruns_fourq(gruns, qpath%points(:,iqpt), wv0(:,iqpt), gvals(:,iqpt))
 end do

 call xmpi_sum(wv0, comm, ierr)
 call xmpi_sum(gvals, comm, ierr)
 ABI_FREE(wv0)
 ABI_FREE(gvals)

end subroutine gruns_qpath
!!***

!----------------------------------------------------------------------

!!****f* m_gruneisen/gruns_qmesh
!! NAME
!!  gruns_qmesh
!!
!! FUNCTION
!!
!! INPUTS
!!  ngqpt(3)=q-mesh divisions
!!  nqshift=Number of shifts used to generated the ab-initio q-mesh.
!!  qshift(3,nqshift)=The shifts of the ab-initio q-mesh.
!!  comm=MPI communicator
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine gruns_qmesh(gruns, ngqpt, nqshift, qshift, comm)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_free'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nqshift,comm
 type(gruns_t),intent(in) :: gruns
!arrays
 integer,intent(in) :: ngqpt(3)
 real(dp),intent(in) :: qshift(3,nqshift)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0,qptopt1=1,bcorr0=0
 integer :: nprocs,my_rank,iqibz,nqbz,nqibz,ierr,ii,nu
 real(dp) :: gavg
 type(t_tetrahedron) :: tetra
 character(len=500) :: msg
!arrays
 integer :: qptrlatt(3,3)
 real(dp),allocatable :: gvals(:,:),wv0(:,:)
 real(dp),allocatable :: qibz(:,:),qbz(:,:),wtq(:)
 !real(dp),allocatable :: wdt(:,:),eig_dos(:,:)

! ************************************************************************

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 ! Generate the q-mesh by finding the IBZ and the corresponding weights.
 ABI_CHECK(all(ngqpt > 0), sjoin("invalid ngqpt:", ltoa(ngqpt)))
 qptrlatt = 0
 do ii=1,3
   qptrlatt(ii,ii) = ngqpt(ii)
 end do

 ! Get IBZ and BZ.
 call kpts_ibz_from_kptrlatt(gruns%cryst_vol(gruns%iv0), qptrlatt, qptopt1, nqshift, qshift, &
   nqibz, qibz, wtq, nqbz, qbz)

 ! Build tetrahedra
 tetra = tetra_from_kptrlatt(gruns%cryst_vol(gruns%iv0), qptopt1, qptrlatt, nqshift, qshift, nqibz, qibz, msg, ierr)
 if (ierr /= 0) MSG_ERROR(msg)

 ABI_CALLOC(wv0, (gruns%natom3, nqibz))
 ABI_CALLOC(gvals, (gruns%natom3, nqibz))

 gavg = zero
 do iqibz=1,nqibz
   if (mod(iqibz, nprocs) /= my_rank) cycle ! mpi-parallelism
   call gruns_fourq(gruns, qibz(:,iqibz), wv0(:,iqibz), gvals(:,iqibz))
   gavg = gavg + wtq(iqibz) * sum(gvals(:,iqibz))
 end do
 gavg = gavg / gruns%natom3

 call xmpi_sum(gavg, comm, ierr)
 call xmpi_sum(wv0, comm, ierr)
 call xmpi_sum(gvals, comm, ierr)

#if 0
 !ABI_MALLOC(wv0_ibz, (nqibz))
 !ABI_MALLOC(wdt, (nene, 2))
 !ABI_CALLOC(eig_dos, (nene, 2))

 ! Accumulate total DOS.
 cnt = 0
 do iqibz=1,nqibz
   do nu=1,gruns%natom3
     cnt = cnt + 1
     if (mod(cnt, nprocs) /= my_rank) cycle ! mpi-parallelism
     wv0_ibz = wv0(nu,:)
     !call tetra_get_onewk(tetra,iqibz,bcorr0,nene,nqibz,wv0_ibz,enemin,enemax,one,wdt)
     !eig_dos = eig_dos + wdt
     wv0_ibz = gvals(nu,:)
   end do
 end do
 !ABI_FREE(wv0_ibz)
 !ABI_FREE(wdt)
 !ABI_FREE(eig_dos)
#endif

 if (my_rank == master) then
   call wrtout(ab_out, sjoin(" Average Gruneisen parameter:", ftoa(gavg, fmt="f8.5")))
 end if

 ABI_FREE(qibz)
 ABI_FREE(wtq)
 ABI_FREE(qbz)
 ABI_FREE(wv0)
 ABI_FREE(gvals)

 call destroy_tetra(tetra)

end subroutine gruns_qmesh
!!***

!----------------------------------------------------------------------

!!****f* m_gruneisen/gruns_free
!! NAME
!!  gruns_free
!!
!! FUNCTION
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine gruns_free(gruns)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(gruns_t),intent(inout) :: gruns

!Local variables-------------------------------
!scalars
 integer :: ii

! ************************************************************************

 if (allocated(gruns%ifc_vol)) then
   do ii=1,size(gruns%cryst_vol)
     call crystal_free(gruns%cryst_vol(ii))
   end do
   ABI_FREE(gruns%cryst_vol)
 end if

 if (allocated(gruns%ddb_vol)) then
   do ii=1,size(gruns%ddb_vol)
     call ddb_free(gruns%ddb_vol(ii))
   end do
   ABI_FREE(gruns%ddb_vol)
 end if

 if (allocated(gruns%ifc_vol)) then
   do ii=1,size(gruns%ifc_vol)
     call ifc_free(gruns%ifc_vol(ii))
   end do
   ABI_FREE(gruns%ifc_vol)
 end if

end subroutine gruns_free
!!***

!----------------------------------------------------------------------

!!****f* m_gruneisen/gruns_anaddb
!! NAME
!!  gruns_anaddb
!!
!! FUNCTION
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine gruns_anaddb(inp, comm)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_free'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: comm
 type(anaddb_dataset_type) :: inp

!Local variables-------------------------------
!scalars
 integer :: ii
 real(dp) :: cpu,wall,gflops
 type(gruns_t) :: gruns
 type(kpath_t) :: qpath
 character(len=500) :: msg

! ************************************************************************

 call cwtime(cpu, wall, gflops, "start")

#ifdef HAVE_NETCDF
 !NCF_CHECK_MSG(nctk_open_create(ncid, strcat(prefix, "_GRUNS.nc"), xmpi_comm_self), "Creating GRUNS.nc")
 !NCF_CHECK(crystal_ncwrite(cryst, ncid))
 !call phonons_ncwrite(ncid,natom,nfineqpath,save_qpoints,weights,save_phfrq,save_phdispl_cart)
 !NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('atomic_mass_units', "dp", "number_of_atom_species")],defmode=.True.))
 !NCF_CHECK(nctk_set_datamode(ncid))
 !NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'atomic_mass_units'), ddb%amu))
 !NCF_CHECK(nf90_close(ncid))
#endif

 gruns = gruns_new(inp%gruns_ddbs, inp, comm)

 ! Compute grunesein on the q-mesh
 if (all(inp%ng2qpt /= 0)) then
   call gruns_qmesh(gruns, inp%ng2qpt, 1, inp%q2shft, comm)
 else
   MSG_WARNING("Cannot compute Grunesein parameters on q-mesh because ng2qpt == 0")
 end if

 ! Compute grunesein on the q-path
 if (inp%nqpath /= 0) then
   call kpath_init(qpath, inp%qpath, gruns%cryst_vol(gruns%iv0)%gprimd, inp%ndivsm)
   call gruns_qpath(gruns, qpath, comm)
   call kpath_free(qpath)
 else
   MSG_WARNING("Cannot compute Grunesein parameters on q-path because nqpath == 0")
 end if

 call gruns_free(gruns)

 call cwtime(cpu,wall,gflops,"stop")
 write(msg,'(2(a,f8.2))')"gruns_anaddb completed. cpu:",cpu,", wall:",wall
 call wrtout(std_out, msg, do_flush=.True.)

end subroutine gruns_anaddb
!!***

!----------------------------------------------------------------------

end module m_gruneisen
