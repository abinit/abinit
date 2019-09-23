!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_gruneisen
!! NAME
!!  m_gruneisen
!!
!! FUNCTION
!!  Objects and methods to compute Gruneisen parameters with central finite differences
!!  of dynamical matrices obtained with different unit cell volumes.
!!
!! COPYRIGHT
!! Copyright (C) 2011-2019 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_gruneisen

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_crystal
 use m_htetra
 use m_ddb
 use m_ddb_hdr
 use m_ifc
 use m_cgtools
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_io_tools,            only : get_unit, open_file
 use m_time,                only : cwtime, cwtime_report
 use m_fstrings,            only : sjoin, itoa, ltoa, ftoa, strcat
 use m_numeric_tools,       only : central_finite_diff, arth
 use m_kpts,                only : kpts_ibz_from_kptrlatt, tetra_from_kptrlatt
 use m_bz_mesh,             only : kpath_t, kpath_new
 use m_anaddb_dataset,      only : anaddb_dataset_type
 use m_dynmat,              only : massmult_and_breaksym, dfpt_phfrq, gtdyn9

 implicit none

 private
!!***

!!****t* m_gruneisen/gruns_t
!! NAME
!! gruns_t
!!
!! FUNCTION
!!  Contains the interatomic force constants for different volumes.
!!  Provides methods to compute phonon bandstructures and gruneisen parameters.
!!
!! SOURCE

 type,public :: gruns_t

   integer :: natom3
    ! 3 * natom

   integer :: nvols
    ! Number of volumes.

   integer :: iv0
    ! Index of the DDB file corresponding to the equilibrium volume V0.

   real(dp) :: v0
    ! Equilibrium volume.

   real(dp) :: delta_vol
    ! Uniform grid spacing for finite difference.

   type(crystal_t),allocatable :: cryst_vol(:)
    ! cryst_vol(nvols)
    ! crystalline structure for the different volumes.

   type(ddb_type),allocatable :: ddb_vol(:)
    ! dbb_vol(nvols)
    ! DDB objects for the different volumes.

   type(ifc_type),allocatable :: ifc_vol(:)
    ! ifc_vol(nvols)
    ! interatomic force constants for the different volumes.

 end type gruns_t

 public :: gruns_new        ! Constructor.
 public :: gruns_qpath      ! Compute Gruneisen parameters on a q-path.
 public :: gruns_qmesh      ! Compute Gruneisen parameters on a q-mesh.
 public :: gruns_free       ! Release memory.
 public :: gruns_anaddb     ! Driver routine called in anaddb.
!!***

contains  !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_gruneisen/gruns_new
!! NAME
!!  gruns_new
!!
!! FUNCTION
!!  Construct new object from a list of DDB files.
!!
!! INPUTS
!!  ddb_paths(:)=Paths of the DDB files (must be ordered by volume)
!!  inp<anaddb_dataset_type>=Anaddb dataset with input variables
!!  comm=MPI communicator
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(gruns_t) function gruns_new(ddb_paths, inp, comm) result(new)

!Arguments ------------------------------------
 integer,intent(in) :: comm
 type(anaddb_dataset_type),intent(in) :: inp
!arrays
 character(len=*),intent(in) :: ddb_paths(:)

!Local variables-------------------------------
 integer,parameter :: natifc0=0,master=0
 integer :: ivol,iblock,natom,ddbun
 integer :: nprocs,my_rank,ierr
 character(len=500) :: msg
 type(ddb_hdr_type) :: ddb_hdr
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
   call ddb_hdr_open_read(ddb_hdr,ddb_paths(ivol),ddbun,DDB_VERSION, dimonly=1)
   natom = ddb_hdr%natom

   call ddb_hdr_free(ddb_hdr)

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
   iblock = new%ddb_vol(ivol)%get_dielt_zeff(new%cryst_vol(ivol), inp%rfmeth, inp%chneut, inp%selectz, dielt, zeff)
   if (iblock == 0) then
     call wrtout(ab_out, sjoin("- Cannot find dielectric tensor and Born effective charges in DDB file:", ddb_paths(ivol)))
     call wrtout(ab_out, "Values initialized with zeros")
   else
     call wrtout(ab_out, sjoin("- Found dielectric tensor and Born effective charges in DDB file:", ddb_paths(ivol)))
   end if

   call ifc_init(new%ifc_vol(ivol), new%cryst_vol(ivol), new%ddb_vol(ivol),&
     inp%brav,inp%asr,inp%symdynmat,inp%dipdip,inp%rfmeth,inp%ngqpt(1:3),inp%nqshft,inp%q1shft,dielt,zeff,&
     inp%nsphere,inp%rifcsph,inp%prtsrlr,inp%enunit,comm)
   ABI_FREE(zeff)
 end do

 ! Consistency check
 ! TODO: Add more tests.
 ABI_CHECK(any(new%nvols == [3, 5, 7, 9]), "Central finite difference requires [3,5,7,9] DDB files")

 new%natom3 = 3 * new%cryst_vol(1)%natom
 new%iv0 = 1 + new%nvols / 2
 new%v0 = new%cryst_vol(new%iv0)%ucvol
 new%delta_vol = new%cryst_vol(new%iv0+1)%ucvol - new%cryst_vol(new%iv0)%ucvol

 ierr = 0
 do ivol=1,new%nvols
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
!!  Compute gruneisen parameters at an arbitrary q-point.
!!
!! INPUTS
!!  qpt(3)=q-point in reduced coordinates.
!!
!! OUTPUT
!!  wvols(3*natom, nvols) = Phonon frequencies for the different volumen.
!!  gvals(3*natom)=Gruneisen parameters evaluated at V0.
!!  dwdq(3,3*natom)=Group velocities at V0 in Cartesian coordinates.
!!  phdispl_cart(2, natom3, natom3, nvols)=Phonon displacement in Cartesian coordinates for the different volumes
!!
!! NOTES
!!
!!  The Gruneisen parameter is given by:
!!
!!     gamma(q,nu) = - (V / w(q,nu)) dw(q,nu)/dV
!!
!!  Using w*2 = <u|D|u> and the Hellmann-Feynmann theorem, one obtains:
!!
!!     gamma(q,nu) = - (V / 2 w(q,nu)**2) <u(q,nu)|dD(q)/dV|u(q,nu)>
!!
!!  The derivative dD/dV is computed via central finite difference.
!!
!! PARENTS
!!      m_gruneisen
!!
!! CHILDREN
!!
!! SOURCE

subroutine gruns_fourq(gruns, qpt, wvols, gvals, dwdq, phdispl_cart)

!Arguments ------------------------------------
!scalars
 type(gruns_t),intent(in) :: gruns
!arrays
 real(dp),intent(in) :: qpt(3)
 real(dp),intent(out) :: wvols(gruns%natom3,gruns%nvols),gvals(gruns%natom3),dwdq(3,gruns%natom3)
 real(dp),intent(out) :: phdispl_cart(2, gruns%natom3, gruns%natom3, gruns%nvols)

!Local variables-------------------------------
!scalars
 integer :: ivol,natom3,nu
 real(dp) :: fact
!arrays
 real(dp) :: dot(2)
 real(dp) :: eigvec(2,gruns%natom3,gruns%natom3,gruns%nvols),d2cart(2,gruns%natom3,gruns%natom3,gruns%nvols)
 real(dp) :: dddv(2,gruns%natom3,gruns%natom3)
 real(dp) :: omat(2,gruns%natom3, gruns%natom3)

! ************************************************************************

 natom3 = gruns%natom3
 do ivol=1,gruns%nvols
   if (ivol == gruns%iv0) then
     ! Compute group velocities for V=V0
     call gruns%ifc_vol(ivol)%fourq(gruns%cryst_vol(ivol), qpt, wvols(:,ivol), phdispl_cart(:,:,:,ivol), &
                    out_d2cart=d2cart(:,:,:,ivol), out_eigvec=eigvec(:,:,:,ivol), dwdq=dwdq)
   else
     call gruns%ifc_vol(ivol)%fourq(gruns%cryst_vol(ivol), qpt, wvols(:,ivol), phdispl_cart(:,:,:,ivol), &
                    out_d2cart=d2cart(:,:,:,ivol), out_eigvec=eigvec(:,:,:,ivol))
   end if

   call massmult_and_breaksym(gruns%cryst_vol(ivol)%natom, gruns%cryst_vol(ivol)%ntypat, &
     gruns%cryst_vol(ivol)%typat, gruns%ifc_vol(ivol)%amu, d2cart(:,:,:,ivol))

   !call zgemm('N','N',natom3,natom3,natom3,cone,d2cart(:,:,:,ivol),natom3,eigvec(:,:,:,ivol),natom3,czero,omat,natom3)
   !do nu=1,natom3
   !  write(std_out,*)"H|psi> - w**2 |psi>",maxval(abs(omat(:,:,nu) - wvols(nu,ivol) ** 2 * eigvec(:,:,nu,ivol)))
   !end do
 end do

 ! Compute dD(q)/dV with central finite difference.
 dddv = zero
 do ivol=1,gruns%nvols
   fact = central_finite_diff(1, ivol, gruns%nvols)
   if (fact /= zero) dddv = dddv + fact * d2cart(:,:,:,ivol)
 end do
 dddv = dddv / gruns%delta_vol

 ! Compute -V0/(2w(q)**2) <u(q)|dD(q)/dq|u(q)>
 call zgemm('N','N',natom3,natom3,natom3,cone,dddv,natom3,eigvec(:,:,:,gruns%iv0),natom3,czero,omat,natom3)
 do nu=1,natom3
   if (abs(wvols(nu, gruns%iv0)) > tol12) then
     dot = cg_zdotc(natom3, eigvec(1,1,nu, gruns%iv0), omat(1,1,nu))
     ! Must change sign if we have a purely imaginary solution i.e. z = iw
     gvals(nu) = -sign(one, wvols(nu, gruns%iv0)) * gruns%v0 * dot(1) / (two * wvols(nu, gruns%iv0)**2)
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
!!  Compute gruneisen parameters and group velocities on a q-path
!!  Write results to file.
!!
!! INPUTS
!!  prefix=Prefix for output files.
!!  qpath<kpath_t>=Object describing the q-path.
!!  ncid=netcdf file id.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      m_gruneisen
!!
!! CHILDREN
!!
!! SOURCE

subroutine gruns_qpath(gruns, prefix, qpath, ncid, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid,comm
 type(gruns_t),intent(in) :: gruns
 type(kpath_t),intent(in) :: qpath
 character(len=*),intent(in) :: prefix

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: nprocs,my_rank,iqpt,ierr,ncerr,unt,iv0,nu,ii
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: gvals_qpath(:,:),wvols_qpath(:,:,:),dwdq_qpath(:,:,:)
 real(dp),allocatable :: phdispl_cart_qpath(:,:,:,:,:)

! ************************************************************************

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 write(msg,'(a,(80a),4a)')ch10,('=',ii=1,80),ch10,ch10,' Calculation of Gruneisen parameters along q-path ',ch10
 call wrtout(std_out, msg)
 !call wrtout(ab_out, msg)

 iv0 = gruns%iv0
 ABI_CALLOC(wvols_qpath, (gruns%natom3, gruns%nvols, qpath%npts))
 ABI_CALLOC(gvals_qpath, (gruns%natom3, qpath%npts))
 ABI_CALLOC(dwdq_qpath, (3, gruns%natom3, qpath%npts))
 ABI_CALLOC(phdispl_cart_qpath, (2, gruns%natom3, gruns%natom3, gruns%nvols, qpath%npts))

 do iqpt=1,qpath%npts
   if (mod(iqpt, nprocs) /= my_rank) cycle ! mpi-parallelism
   call gruns_fourq(gruns, qpath%points(:,iqpt), wvols_qpath(:,:,iqpt), gvals_qpath(:,iqpt), &
                    dwdq_qpath(:,:,iqpt), phdispl_cart_qpath(:,:,:,:,iqpt))
 end do

 call xmpi_sum(wvols_qpath, comm, ierr)
 call xmpi_sum(gvals_qpath, comm, ierr)
 call xmpi_sum(dwdq_qpath, comm, ierr)
 call xmpi_sum(phdispl_cart_qpath, comm, ierr)

 ! Write text files with phonon frequencies and gruneisen on the path.
 if (my_rank == master) then
   if (open_file(strcat(prefix, "_GRUNS_QPATH"), msg, newunit=unt, form="formatted", action="write") /= 0) then
     MSG_ERROR(msg)
   end if
   write(unt,'(a)')'# Phonon band structure, Gruneisen parameters and group velocity'
   write(unt,'(a)')"# Energy in Hartree, DOS in states/Hartree"
   call qpath%print(unit=unt, pre="#")
   write(unt,'(5a)')&
     "# phfreq(mode=1) gruneisen(mode=1) velocity(mode=1)    phfreq(mode=2) gruneisen(mode=2) velocity(mode=2)   ..."
   do iqpt=1,qpath%npts
     do nu=1,gruns%natom3
       write(unt, "(3es17.8)", advance="no") &
         wvols_qpath(nu, iv0, iqpt), gvals_qpath(nu, iqpt), sum(dwdq_qpath(1:3, nu, iqpt) ** 2)
     end do
     write(unt, "(a)")" "
   end do
   close(unt)
 end if

#ifdef HAVE_NETCDF
 if (my_rank == master .and. ncid /= nctk_noid) then
   ncerr = nctk_def_dims(ncid, [nctkdim_t("gruns_nqpath", qpath%npts)], defmode=.True.)
   NCF_CHECK(ncerr)

   ncerr = nctk_def_arrays(ncid, [ &
     ! q-points of the path
     nctkarr_t("gruns_qpath", "dp", "three, gruns_nqpath"), &
     ! gruneisen parameters on the path
     nctkarr_t("gruns_gvals_qpath", "dp", "number_of_phonon_modes, gruns_nqpath"), &
     ! phonon frequencies at the different volumes
     nctkarr_t("gruns_wvols_qpath", "dp", "number_of_phonon_modes, gruns_nvols, gruns_nqpath"), &
     ! group velocities at V0 in Cartesian coordinates.
     nctkarr_t("gruns_dwdq_qpath", "dp", "three, number_of_phonon_modes, gruns_nqpath"), &
     ! displacements for the different volumes.
     nctkarr_t("gruns_phdispl_cart_qpath", "dp", &
         "two, number_of_phonon_modes, number_of_phonon_modes, gruns_nvols, gruns_nqpath") &
   ])
   NCF_CHECK(ncerr)

   ! Write data.
   NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gruns_qpath"), qpath%points))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gruns_gvals_qpath"), gvals_qpath))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gruns_wvols_qpath"), wvols_qpath))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gruns_dwdq_qpath"), dwdq_qpath))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gruns_phdispl_cart_qpath"), phdispl_cart_qpath))
 end if
#endif

 ABI_FREE(wvols_qpath)
 ABI_FREE(gvals_qpath)
 ABI_FREE(dwdq_qpath)
 ABI_FREE(phdispl_cart_qpath)

end subroutine gruns_qpath
!!***

!----------------------------------------------------------------------

!!****f* m_gruneisen/gruns_qmesh
!! NAME
!!  gruns_qmesh
!!
!! FUNCTION
!!  Compute gruneisen parameters and group velocities on a q-mesh.
!!  Save results to file.
!!
!! INPUTS
!!  prefix=Prefix for output files.
!!  dosdeltae=Step for the frequency mesh.
!!  ngqpt(3)=q-mesh divisions
!!  nshiftq=Number of shifts used to generated the ab-initio q-mesh.
!!  shiftq(3,nshiftq)=The shifts of the ab-initio q-mesh.
!!  ncid=netcdf file id.
!!  comm=MPI communicator
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gruneisen
!!
!! CHILDREN
!!
!! SOURCE

subroutine gruns_qmesh(gruns, prefix, dosdeltae, ngqpt, nshiftq, shiftq, ncid, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nshiftq,ncid,comm
 real(dp),intent(in) :: dosdeltae !,dossmear
 type(gruns_t),intent(in) :: gruns
 character(len=*),intent(in) :: prefix
!arrays
 integer,intent(in) :: ngqpt(3)
 real(dp),intent(in) :: shiftq(3,nshiftq)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0,qptopt1=1,bcorr0=0
 integer :: nprocs,my_rank,iqibz,nqbz,nqibz,ierr,ii,nu,ncerr,nomega,cnt,unt,io
 real(dp) :: gavg,omega_min,omega_max,v2
 type(htetra_t) :: tetra
 character(len=500) :: msg
!arrays
 integer :: qptrlatt(3,3)
 real(dp),allocatable :: gvals_qibz(:,:),wvols_qibz(:,:,:),dwdq_qibz(:,:,:)
 real(dp),allocatable :: qibz(:,:),qbz(:,:),wtq(:)
 real(dp),allocatable :: wdt(:,:),wdos(:,:),grdos(:,:),gr2dos(:,:),wibz(:),omega_mesh(:)
 real(dp),allocatable :: vdos(:,:),v2dos(:,:)
 real(dp),allocatable :: phdispl_cart_qibz(:,:,:,:,:)

! ************************************************************************

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 write(msg,'(a,(80a),4a)')ch10,('=',ii=1,80),ch10,ch10,' Calculation of Gruneisen DOSes ',ch10
 call wrtout(std_out, msg)

 ! Generate the q-mesh by finding the IBZ and the corresponding weights.
 ABI_CHECK(all(ngqpt > 0), sjoin("invalid ngqpt:", ltoa(ngqpt)))
 qptrlatt = 0
 do ii=1,3
   qptrlatt(ii,ii) = ngqpt(ii)
 end do

 ! Get IBZ and BZ.
 call kpts_ibz_from_kptrlatt(gruns%cryst_vol(gruns%iv0), qptrlatt, qptopt1, nshiftq, shiftq, &
   nqibz, qibz, wtq, nqbz, qbz)

 ! Build tetrahedra
 tetra = tetra_from_kptrlatt(gruns%cryst_vol(gruns%iv0), qptopt1, qptrlatt, nshiftq, shiftq, nqibz, qibz, comm, msg, ierr)
 if (ierr /= 0) MSG_ERROR(msg)

 ABI_CALLOC(wvols_qibz, (gruns%natom3, gruns%nvols, nqibz))
 ABI_CALLOC(gvals_qibz, (gruns%natom3, nqibz))
 ABI_CALLOC(dwdq_qibz, (3, gruns%natom3, nqibz))
 ABI_CALLOC(phdispl_cart_qibz, (2, gruns%natom3, gruns%natom3, gruns%nvols, nqibz))

 gavg = zero
 do iqibz=1,nqibz
   if (mod(iqibz, nprocs) /= my_rank) cycle ! mpi-parallelism
   call gruns_fourq(gruns, qibz(:,iqibz), wvols_qibz(:,:,iqibz), gvals_qibz(:,iqibz), &
                    dwdq_qibz(:,:,iqibz), phdispl_cart_qibz(:,:,:,:,iqibz))
   gavg = gavg + wtq(iqibz) * sum(gvals_qibz(:,iqibz))
 end do
 gavg = gavg / gruns%natom3

 call xmpi_sum(gavg, comm, ierr)
 call xmpi_sum(wvols_qibz, comm, ierr)
 call xmpi_sum(gvals_qibz, comm, ierr)
 call xmpi_sum(dwdq_qibz, comm, ierr)
 call xmpi_sum(phdispl_cart_qibz, comm, ierr)

 omega_min = gruns%ifc_vol(gruns%iv0)%omega_minmax(1)
 omega_max = gruns%ifc_vol(gruns%iv0)%omega_minmax(2)
 nomega = nint((omega_max - omega_min) / dosdeltae) + 1
 nomega = max(6, nomega) ! Ensure Simpson integration will be ok

 ABI_MALLOC(omega_mesh, (nomega))
 omega_mesh = arth(omega_min, dosdeltae, nomega)
 omega_max = omega_mesh(nomega)
 !write(std_out,*)"hello",omega_min,omega_max,dosdeltae,(omega_max-omega_min) / (nomega-1)
 ABI_MALLOC(wibz, (nqibz))
 ABI_MALLOC(wdt, (nomega, 2))
 ABI_CALLOC(wdos, (nomega, 2))
 ABI_CALLOC(grdos, (nomega, 2))
 ABI_CALLOC(gr2dos, (nomega, 2))
 ABI_CALLOC(vdos, (nomega, 2))
 ABI_CALLOC(v2dos, (nomega, 2))

 ! Compute DOSes.
 cnt = 0
 do iqibz=1,nqibz
   do nu=1,gruns%natom3
     cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle ! mpi-parallelism
     wibz = wvols_qibz(nu, gruns%iv0, :)
     call tetra%get_onewk(iqibz,bcorr0,nomega,nqibz,wibz,omega_min,omega_max,one,wdt)
     wdt = wdt*wtq(iqibz)
     wdos = wdos + wdt
     grdos = grdos + wdt * gvals_qibz(nu,iqibz)
     gr2dos = gr2dos + wdt * gvals_qibz(nu,iqibz) ** 2
     v2 = sum(dwdq_qibz(1:3,nu,iqibz) ** 2)
     vdos = vdos + wdt * sqrt(v2)
     v2dos = v2dos + wdt * v2
   end do
 end do

 call xmpi_sum(wdos, comm, ierr)
 call xmpi_sum(grdos, comm, ierr)
 call xmpi_sum(gr2dos, comm, ierr)
 call xmpi_sum(vdos, comm, ierr)
 call xmpi_sum(v2dos, comm, ierr)

 if (my_rank == master) then
   call wrtout(ab_out, sjoin(" Average Gruneisen parameter:", ftoa(gavg, fmt="f8.5")))

   ! Write text files with Gruneisen and DOSes.
   if (open_file(strcat(prefix, "_GRUNS_DOS"), msg, newunit=unt, form="formatted", action="write") /= 0) then
     MSG_ERROR(msg)
   end if
   write(unt,'(a)')'# Phonon density of states, Gruneisen DOS and phonon group velocity DOS'
   write(unt,'(a)')"# Energy in Hartree, DOS in states/Hartree"
   write(unt,'(a,i0)')'# Tetrahedron method with nqibz= ',nqibz
   write(unt,"(a,f8.5)")"# Average Gruneisen parameter:", gavg
   write(unt,'(5a)') &
     "# omega PH_DOS Gruns_DOS Gruns**2_DOS Vel_DOS  Vel**2_DOS  PH_IDOS Gruns_IDOS Gruns**2_IDOS Vel_IDOS Vel**2_IDOS"
   do io=1,nomega
     write(unt, "(11es17.8)")omega_mesh(io), &
       wdos(io,1), grdos(io,1), gr2dos(io,1), vdos(io,1), v2dos(io,1), &
       wdos(io,2), grdos(io,2), gr2dos(io,2), vdos(io,2), v2dos(io,2)
   end do
   close(unt)
 end if

#ifdef HAVE_NETCDF
 ! Write netcdf files.
 if (my_rank == master .and. ncid /= nctk_noid) then
   ncerr = nctk_def_dims(ncid, [ &
     nctkdim_t("gruns_nqibz", nqibz), nctkdim_t('gruns_nshiftq', nshiftq), &
     nctkdim_t('gruns_nomega', nomega)], defmode=.True.)
   NCF_CHECK(ncerr)

   ncerr = nctk_def_arrays(ncid, [ &
    ! q-point sampling in IBZ,
    nctkarr_t("gruns_qptrlatt", "int", "three, three"), &
    nctkarr_t("gruns_shiftq", "dp", "three, gruns_nshiftq"), &
    nctkarr_t("gruns_qibz", "dp", "three, gruns_nqibz"), &
    nctkarr_t("gruns_wtq", "dp", "gruns_nqibz"), &
    ! gruneisen parameters in IBZ
    ! phonon frequencies at the different volumes,
    ! group velocities at V0 in Cartesian coordinates.
    nctkarr_t("gruns_gvals_qibz", "dp", "number_of_phonon_modes, gruns_nqibz"), &
    nctkarr_t("gruns_wvols_qibz", "dp", "number_of_phonon_modes, gruns_nvols, gruns_nqibz"), &
    nctkarr_t("gruns_dwdq_qibz", "dp", "three, number_of_phonon_modes, gruns_nqibz"), &
    ! displacements for the different volumes.
    nctkarr_t("gruns_phdispl_cart_qibz", "dp", &
        "two, number_of_phonon_modes, number_of_phonon_modes, gruns_nvols, gruns_nqibz"), &
    ! DOSes and IDOSes
    nctkarr_t("gruns_omega_mesh", "dp", "gruns_nomega"), &
    nctkarr_t("gruns_wdos", "dp", "gruns_nomega, two"), &
    nctkarr_t("gruns_grdos", "dp", "gruns_nomega, two"), &
    nctkarr_t("gruns_gr2dos", "dp", "gruns_nomega, two"), &
    nctkarr_t("gruns_v2dos", "dp", "gruns_nomega, two"), &
    nctkarr_t("gruns_vdos", "dp", "gruns_nomega, two") &
   ])
   NCF_CHECK(ncerr)

   ! Write data.
   NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gruns_qptrlatt"), qptrlatt))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gruns_shiftq"), shiftq))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gruns_qibz"), qibz))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gruns_wtq"), wtq))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gruns_gvals_qibz"), gvals_qibz))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gruns_wvols_qibz"), wvols_qibz))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gruns_dwdq_qibz"), dwdq_qibz))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gruns_phdispl_cart_qibz"), phdispl_cart_qibz))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gruns_omega_mesh"), omega_mesh))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gruns_wdos"), wdos))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gruns_grdos"), grdos))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gruns_gr2dos"), gr2dos))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gruns_v2dos"), v2dos))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gruns_vdos"), vdos))
 end if
#endif

 ABI_FREE(qibz)
 ABI_FREE(wtq)
 ABI_FREE(qbz)
 ABI_FREE(wvols_qibz)
 ABI_FREE(gvals_qibz)
 ABI_FREE(dwdq_qibz)
 ABI_FREE(phdispl_cart_qibz)
 ABI_FREE(omega_mesh)
 ABI_FREE(wibz)
 ABI_FREE(wdt)
 ABI_FREE(wdos)
 ABI_FREE(grdos)
 ABI_FREE(gr2dos)
 ABI_FREE(v2dos)
 ABI_FREE(vdos)

 call tetra%free()

end subroutine gruns_qmesh
!!***

!----------------------------------------------------------------------

!!****f* m_gruneisen/gruns_free
!! NAME
!!  gruns_free
!!
!! FUNCTION
!!  Free dynamic memory.
!!
!! PARENTS
!!      m_gruneisen
!!
!! CHILDREN
!!
!! SOURCE

subroutine gruns_free(gruns)

!Arguments ------------------------------------
!array
 type(gruns_t),intent(inout) :: gruns

!Local variables-------------------------------
!scalars
 integer :: ii

! ************************************************************************

 if (allocated(gruns%ifc_vol)) then
   do ii=1,size(gruns%cryst_vol)
     call gruns%cryst_vol(ii)%free()
   end do
   ABI_FREE(gruns%cryst_vol)
 end if

 if (allocated(gruns%ddb_vol)) then
   do ii=1,size(gruns%ddb_vol)
     call gruns%ddb_vol(ii)%free()
   end do
   ABI_FREE(gruns%ddb_vol)
 end if

 if (allocated(gruns%ifc_vol)) then
   do ii=1,size(gruns%ifc_vol)
     call gruns%ifc_vol(ii)%free()
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
!!  Driver routine called in anaddb to compute Gruneisen parameters.
!!
!! INPUTS
!!  inp<anaddb_dataset_type: :: anaddb input variables.
!!  prefix=Prefix for output files
!!  comm=MPI communicator
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!
!! SOURCE

subroutine gruns_anaddb(inp, prefix, comm)

!Arguments ------------------------------------
 integer,intent(in) :: comm
 character(len=*),intent(in) :: prefix
 type(anaddb_dataset_type) :: inp

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: ii,nprocs,my_rank,ncid,iv0
#ifdef HAVE_NETCDF
 integer :: ncerr
#endif
 real(dp) :: cpu,wall,gflops
 type(gruns_t),target :: gruns
 type(kpath_t) :: qpath

! ************************************************************************

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 ABI_CHECK(inp%ifcflag == 1, "Gruneisen requires ifcflag == 1")

 call cwtime(cpu, wall, gflops, "start")

 gruns = gruns_new(inp%gruns_ddbs, inp, comm)
 iv0 = gruns%iv0

 ncid = nctk_noid
#ifdef HAVE_NETCDF
 if (my_rank == master) then
   NCF_CHECK_MSG(nctk_open_create(ncid, strcat(prefix, "_GRUNS.nc"), xmpi_comm_self), "Creating _GRUNS.nc")

   ! Write structure corresponding to iv0
   NCF_CHECK(gruns%cryst_vol(iv0)%ncwrite(ncid))

   ! Add important dimensions and additional metadata.
   ncerr = nctk_def_dims(ncid, [ &
     nctkdim_t("gruns_nvols", gruns%nvols), &
     nctkdim_t('number_of_phonon_modes', gruns%natom3)], defmode=.True.)
   NCF_CHECK(ncerr)
   ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "gruns_iv0"])
   NCF_CHECK(ncerr)

   ! Add lattice parameters and positions of the `gruns_nvols` structures so
   ! that we can easily reconstruct structure objects in AbiPy.
   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("gruns_rprimd", "dp", "three, three, gruns_nvols"), &
     nctkarr_t("gruns_xred", "dp", "three, number_of_atoms, gruns_nvols") &
   ])
   NCF_CHECK(ncerr)

   ncerr = nctk_write_iscalars(ncid, [character(len=nctk_slen) :: "gruns_iv0"], [iv0], datamode=.True.)
   NCF_CHECK(ncerr)

   do ii=1,gruns%nvols
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gruns_rprimd"), gruns%cryst_vol(ii)%rprimd, start=[1,1,ii]))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gruns_xred"), gruns%cryst_vol(ii)%xred, start=[1,1,ii]))
   end do

   !call phonons_ncwrite(ncid,natom,nfineqpath,save_qpoints,weights,save_phfrq,save_phdispl_cart)
   !NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'atomic_mass_units'), ddb%amu))
 end if
#endif

 ! Compute gruneisen parameters on the q-mesh.
 if (all(inp%ng2qpt /= 0)) then
   call gruns_qmesh(gruns, prefix, inp%dosdeltae, inp%ng2qpt, 1, inp%q2shft, ncid, comm)
 else
   MSG_WARNING("Cannot compute Gruneisen parameters on q-mesh because ng2qpt == 0")
 end if

 ! Compute gruneisen on the q-path.
 if (inp%nqpath /= 0) then
   qpath = kpath_new(inp%qpath, gruns%cryst_vol(iv0)%gprimd, inp%ndivsm)
   call gruns_qpath(gruns, prefix, qpath, ncid, comm)
   call qpath%free()
 else
   MSG_WARNING("Cannot compute Gruneisen parameters on q-path because nqpath == 0")
 end if

 ! Compute speed of sound for V0.
 if (inp%vs_qrad_tolkms(1) > zero) then
   call gruns%ifc_vol(iv0)%speedofsound(gruns%cryst_vol(iv0), inp%vs_qrad_tolkms, ncid, comm)
 end if

 ! Now treat the second list of vectors (only at the Gamma point, but can include non-analyticities)
 if (my_rank == master .and. inp%nph2l /= 0 .and. inp%ifcflag == 1) then
   call gruns%ifc_vol(iv0)%calcnwrite_nana_terms(gruns%cryst_vol(iv0), inp%nph2l, inp%qph2l, inp%qnrml2, ncid)
 end if

#ifdef HAVE_NETCDF
 if (my_rank == master) then
   NCF_CHECK(nf90_close(ncid))
 end if
#endif

 call gruns_free(gruns)
 call cwtime_report("gruns_anaddb", cpu, wall, gflops)

end subroutine gruns_anaddb
!!***

!----------------------------------------------------------------------

end module m_gruneisen
