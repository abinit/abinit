!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_fstab
!! NAME
!!
!! FUNCTION
!!  Tools for the management of a set of Fermi surface k-points
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group (MG, MVer)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_fstab

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_krank
 use m_htetra
 use m_ebands
 use m_crystal
 use m_dtset

 use m_time,           only : cwtime, cwtime_report
 use m_fstrings,       only : itoa, sjoin, ktoa
 use m_numeric_tools,  only : bisect, isdiagmat, get_diag
 use m_symtk,          only : matr3inv
 use defs_datatypes,   only : ebands_t
 use m_special_funcs,  only : gaussian
 use m_kpts,           only : kpts_timrev_from_kptopt, listkk, smpbz

 implicit none

 private
!!***

!!****t* m_fstab/fstab_t
!! NAME
!! fstab_t
!!
!! FUNCTION
!!  Tables with the correspondence between points of the Fermi surface (FS) and the k-points in the
!!  IBZ (k-points found in ebands_t). We use `nsppol` fstab_t objects to account for spin polarization.
!!
!! SOURCE

 type,public :: fstab_t

   integer :: spin
    ! Spin index

   integer :: nkfs
    ! Number of k-points on the Fermi-surface (FS-BZ).

   integer :: nktot
    ! Total number of k-points in the initial mesh.

   integer :: nkibz
    ! Number of points in the IBZ

   integer :: bmin, bmax
    ! Min and max band index included in the calculation.
    ! Note that these values are obtained by taking the union over the k-points in the FS
    ! For each k-point, we usually have a different number of states crossing eF given by bstart_cnt_ibz

   integer :: maxnb
   ! Max number of bands on the FS (bmax - bmin + 1)

   integer :: eph_intmeth
   ! Integration method. 1 for gaussian, 2 for tetrahedra

   integer :: nene
   ! Number of chemical potential values used for inelastic integration

   real(dp) :: eph_fsmear
   ! Gaussian broadening. Negative value activates adaptive gaussian broadening.
   ! https://journals.aps.org/prb/pdf/10.1103/PhysRevB.92.075405

   real(dp) :: min_smear = tol9
   ! Used for the adaptive gaussian broadening: use min_smear if the broadening computed from the group velocity
   ! is smaller than this value to avoid divergences in the gaussian

   real(dp) :: enemin
   ! Minimal chemical potential value used for inelastic integration

   real(dp) :: deltaene
   ! Chemical potential increment for inelastic integration

   type(krank_t) :: krank
   ! rank/inverse_rank pair for the k-points on the FS (kpts).

   integer,allocatable :: indkk_fs(:,:)
   ! (6, nkfs)
   ! Tables giving the correspondence between a point in the FS-BZ and the IBZ computed in listkk.
   !
   !   indkk_fs(1,:)      Mapping FS-BZ --> k-points in the IBZ (taken from ebands_t)
   !   indkk_fs(2,:)      The index of the symmetry S such that kfs = tim_sign * S(k_ibz) + G0
   !   indkk_fs(3:5,:)    The reduced components of G0.
   !   indkk_fs(6,:)      1 if time-reversal was used to generate the k-point, 0 otherwise

   integer,allocatable :: bstart_cnt_ibz(:,:)
    ! (2, nkibz)
    ! The indices of the bands within the energy window (depends on fsk)
    ! Note that we use the k-point index in the IBZ.
    !
    !   bstcnt(1, :) The index of the first band inside the energy window (start)
    !   bstcnt(2, :) Number of bands on the FS (count)

   real(dp) :: kmesh_cartvec(3,3)
    ! vectors defining the k-mesh (stored as column vector in Cartesian coords.
    ! Used to implement the adaptive gaussian broadening.
    ! NB: two_pi is included.

   real(dp),allocatable :: kpts(:,:)
   ! (3, nkfs)
   ! Reduced coordinates of the k-points on the Fermi surface.

   real(dp),allocatable :: vk(:,:), vkq(:,:)
   ! (3, mnb)
   ! Velocities in cartesian cordinates. Used to implement the adaptive gaussian broadening
   ! Values are filled by the called (phgamma) inside the loop over k-points.

   real(dp),allocatable :: tetra_wtk(:,:)
   ! (maxnb, nkibz)
   ! Weights for FS integration with tetrahedron method
   ! Note that the weights are dimensioned with nkibz
   ! (1, :) corresponds to %bmin

   real(dp),allocatable :: tetra_wtk_ene(:,:,:)
   ! (maxnb, nkibz, nene)
   ! Weights for FS integration with tetrahedron method for all chemical potentials
   ! Note that the weights are dimensioned with nkibz
   ! (1, :) corresponds to %bmin

   real(dp),allocatable :: dbldelta_tetra_weights_kfs(:,:,:)
   ! (maxnb, maxnb, nkfs)
   ! (1, 1, :) corresponds to %bmin

 contains

 procedure :: free => fstab_free
  ! Free memory.

 procedure :: findkg0 => fstab_findkg0
  ! Find the index of the k-point on the FS

 procedure :: setup_qpoiint => fstab_setup_qpoint

 procedure :: get_dbldelta_weights => fstab_get_dbldelta_weights

 end type fstab_t

 public :: fstab_init    ! Initialize the object.
 public :: fstab_print   ! Print the object
!!***

 !FIXME These routines are deprecated
 public :: mkqptequiv

!----------------------------------------------------------------------

contains  !============================================================
!!***

!!****f* m_fstab/fstab_free
!! NAME
!!  fstab_free
!!
!! FUNCTION
!!  Free memory
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_phgamma
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine fstab_free(fstab)

!Arguments ------------------------------------
 class(fstab_t),intent(inout) :: fstab

! ************************************************************************

 ! integer
 ABI_SFREE(fstab%indkk_fs)
 ABI_SFREE(fstab%bstart_cnt_ibz)

 ! real
 ABI_SFREE(fstab%kpts)
 ABI_SFREE(fstab%vk)
 ABI_SFREE(fstab%vkq)
 ABI_SFREE(fstab%tetra_wtk)
 ABI_SFREE(fstab%tetra_wtk_ene)
 ABI_SFREE(fstab%dbldelta_tetra_weights_kfs)

 ! types
 call fstab%krank%free()

end subroutine fstab_free
!!***

!----------------------------------------------------------------------

!!****f* m_fstab/fstab_init
!! NAME
!!  fstab_init
!!
!! FUNCTION
!!  Initialize the tables for the FS integration.
!!
!! INPUTS
!!  ebands<ebands_t>=The object describing the band structure.
!!  cryst<crystal_t>=Info on the crystalline structure.
!!  dtset:
!!    eph_fsewin=Energy window in Hartree. Only states in [efermi-fsewin, efermi+fsewin] are included.
!!    eph_intmeth=Flag selecting the integration method.
!!    kptrlatt(3,3)=k-point lattice specification
!!    nshiftk= number of shift vectors.
!!    shiftk(3,nshiftk)=shift vectors for k point generation
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  fstab(nsppol)=Tables with the correspondence between points of the Fermi surface (FS)
!!     and the k-points in ebands_t.
!!
!! TODO
!!  Use a different algorithm to select k-points if tetra. First compute tetra weights
!!  then k-points contributing to FS integral are selected according to some threshold.
!!
!! PARENTS
!!      m_phgamma
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine fstab_init(fstab, ebands, cryst, dtset, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 type(ebands_t),intent(in) :: ebands
 type(crystal_t),intent(in) :: cryst
 type(dataset_type),intent(in) :: dtset
!arrays
 type(fstab_t),target,intent(out) :: fstab(ebands%nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: option0 = 0, brav1 = 1, bcorr0 = 0
 integer :: nkfs,spin,band,nband_k,i1,i2,ib,blow,ik_bz,ik_ibz,nkibz,sppoldbl,timrev
 integer :: ik,mkpt,nkbz,ierr, bstart_k,bstop_k,nene,ifermi
 real(dp),parameter :: max_occ1=one
 real(dp) :: elow,ehigh,ebis,enemin,enemax,deltaene,dksqmax,cpu,wall,gflops
 logical :: inwin
 character(len=80) :: errstr
 character(len=500) :: msg
 type(fstab_t),pointer :: fs
 type(htetra_t) :: tetra
!arrays
 integer :: kptrlatt(3,3)
 integer,allocatable :: full2ebands(:,:),bz2ibz(:), fs2bz(:),indkk(:,:) !,fs2ibz(:)
 real(dp) :: rlatt(3,3), klatt(3,3)
 real(dp),allocatable :: kbz(:,:), tmp_eigen(:),bdelta(:,:),btheta(:,:)

! *************************************************************************

 !@fstab_t
 call cwtime(cpu, wall, gflops, "start")

 if (any(cryst%symrel(:,:,1) /= identity_3d) .and. any(abs(cryst%tnons(:,1)) > tol10) ) then
  MSG_ERROR('The first symmetry is not the identity operator!')
 end if

 nkibz = ebands%nkpt

 kptrlatt = dtset%kptrlatt
 !call kpts_ibz_from_kptrlatt(cryst, kptrlatt, ebands%kptopt, dtset%nshiftk, dtset%shiftk, &
 ! nkibz, kibz, wtk, nkbz, kbz, &
 ! new_kptrlatt, new_shiftk)  ! Optional

 ! Call smpbz to get the full grid of k-points `kbz`
 ! brav1=1 is able to treat all bravais lattices (same option used in getkgrid)
 mkpt= kptrlatt(1,1)*kptrlatt(2,2)*kptrlatt(3,3) &
   +kptrlatt(1,2)*kptrlatt(2,3)*kptrlatt(3,1) &
   +kptrlatt(1,3)*kptrlatt(2,1)*kptrlatt(3,2) &
   -kptrlatt(1,2)*kptrlatt(2,1)*kptrlatt(3,3) &
   -kptrlatt(1,3)*kptrlatt(2,2)*kptrlatt(3,1) &
   -kptrlatt(1,1)*kptrlatt(2,3)*kptrlatt(3,2)

 ABI_MALLOC_OR_DIE(kbz, (3, mkpt), ierr)

 call smpbz(brav1, std_out, kptrlatt, mkpt, nkbz, dtset%nshiftk, option0, dtset%shiftk, kbz)

 ! Find correspondence BZ --> IBZ
 ! Note that we use symrel so these tables can be used to symmetrize wavefunctions.
 timrev = kpts_timrev_from_kptopt(ebands%kptopt)
 sppoldbl = 1; if (any(cryst%symafm == -1) .and. ebands%nsppol == 1) sppoldbl=2
 ABI_MALLOC(indkk, (nkbz*sppoldbl, 6))

 call listkk(dksqmax, cryst%gmet, indkk, ebands%kptns, kbz, ebands%nkpt, nkbz, cryst%nsym,&
    sppoldbl, cryst%symafm, cryst%symrel, timrev, comm, exit_loop=.True., use_symrec=.False.)

 if (dksqmax > tol12) then
   write(msg, '(7a,es16.6,4a)' ) &
   'The WFK file cannot be used to start the present calculation ',ch10, &
   'It was asked that the wavefunctions be accurate, but',ch10, &
   'at least one of the k points could not be generated from a symmetrical one.',ch10, &
   'dksqmax= ',dksqmax,ch10, &
   'Action: check your WFK file and k-point input variables',ch10, &
   '        (e.g. kptopt or shiftk might be wrong in the present dataset or the preparatory one.'
   MSG_ERROR(msg)
 end if

 call cwtime_report(" fstab_init%listkk", cpu, wall, gflops)

 ABI_MALLOC(full2ebands, (6, nkbz))
 full2ebands = 0

 do ik_bz=1,nkbz
   full2ebands(1, ik_bz) = indkk(ik_bz, 1)     ! ik_ibz
   full2ebands(2, ik_bz) = indkk(ik_bz, 2)     ! isym
   full2ebands(3:5, ik_bz) = indkk(ik_bz, 3:5) ! g0
   full2ebands(6, ik_bz) = indkk(ik_bz, 6)     ! itimrev
 end do

 ! Select only the k-points in the BZ that are sufficiently close to the FS.
 ! FIXME: Do not know why but lamda depends on eph_fsewin if gaussian
 ABI_CHECK(dtset%eph_fsewin > tol12, "dtset%eph_fsewin < tol12")
 elow = ebands%fermie - dtset%eph_fsewin
 ehigh = ebands%fermie + dtset%eph_fsewin
 ebis = elow - abs(elow) * 0.001_dp

 ! Allocate workspace arrays.
 !ABI_MALLOC(fs2ibz, (nkbz))
 ABI_MALLOC(fs2bz, (nkbz))

 do spin=1,ebands%nsppol
   fs => fstab(spin)
   fs%spin = spin
   ABI_MALLOC(fs%bstart_cnt_ibz, (2, nkibz))
   fs%bstart_cnt_ibz = -1

   ! Find k-points on the FS associated to this spin.
   nkfs = 0
   do ik_bz=1,nkbz
     ik_ibz = full2ebands(1, ik_bz)
     nband_k = ebands%nband(ik_ibz + (spin-1)*nkibz)

     blow = bisect(ebands%eig(:nband_k, ik_ibz, spin), ebis)
     if (blow == 0) blow = 1
     !if (blow == nband_k .or. blow == 0) cycle ! out of range
     !write(std_out,*)"here with blow: ", blow,nband_k
     !write(std_out,*)"eig_blow, eig_max, elow, ehigh:", ebands%eig(blow, ik_ibz, spin), ebands%eig(nband_k, ik_ibz, spin), elow,ehigh

     inwin = .False.; i1 = huge(1); i2 = -1
     do band=blow,nband_k
        !if (ebands%eig(band, ik_ibz, spin) > ehigh) exit
        !write(std_out,*)band, ebands%eig(band, ik_ibz, spin) >= elow, ebands%eig(band, ik_ibz, spin) <= ehigh
        if (ebands%eig(band, ik_ibz, spin) >= elow .and. ebands%eig(band, ik_ibz, spin) <= ehigh) then
          inwin = .True.; i1 = min(i1, band); i2 = max(i2, band)
        end if
     end do

     if (inwin) then
       ! Add this k-point and the corresponding bands.
       !write(std_out,*)"in win"
       nkfs = nkfs + 1
       !fs2ibz(nkfs) = ik_ibz
       fs2bz(nkfs) = ik_bz
       if (any(fs%bstart_cnt_ibz(:, ik_ibz) /= [-1, -1])) then
         ABI_CHECK(all(fs%bstart_cnt_ibz(:, ik_ibz) == [i1, i2-i1+1]), "bstart_cnt_ibz!")
       end if
       fs%bstart_cnt_ibz(:, ik_ibz) = [i1, i2-i1+1]
     end if
   end do ! ik_bz

   ! @fstab_t
   ! Build fstab_t for this spin.
   fs%nkibz = nkibz; fs%nkfs = nkfs; fs%nktot = nkbz
   ABI_MALLOC(fs%kpts, (3, nkfs))
   ABI_MALLOC(fs%indkk_fs, (6, nkfs))
   do ik=1,nkfs
     !ik_ibz = fs2ibz(ik)
     ik_bz = fs2bz(ik)
     fs%kpts(:,ik) = kbz(:, ik_bz)
     fs%indkk_fs(:, ik) = full2ebands(:, ik_bz)
   end do

   ! Define band indices enclosing states on the FS.
   ! Note that we need all k-points for a given band when computing weights with tetrahedron.
   ! This means that we have to be carefull when selecting the weight associated to a given pair
   ! (band_kq, kq), (band_k, k).
   ! Then we have to rearrange the weights
   fs%bmin = huge(1); fs%bmax = -huge(1)
   do ik_ibz=1,nkibz
     if (fs%bstart_cnt_ibz(1, ik_ibz) /= -1) then
       fs%bmin = min(fs%bmin, fs%bstart_cnt_ibz(1,ik_ibz))
     end if
     if (fs%bstart_cnt_ibz(2, ik_ibz) /= -1) then
       fs%bmax = max(fs%bmax, fs%bstart_cnt_ibz(1,ik_ibz) + fs%bstart_cnt_ibz(2,ik_ibz) - 1)
     end if
   end do
   !write(std_out,*)"bmin, bmax for tetra: ",fs%bmin, fs%bmax
   ABI_CHECK(fs%bmin /= huge(1) .and. fs%bmax /= -huge(1), "No point on the Fermi surface!")
   fs%maxnb = fs%bmax - fs%bmin + 1

   ABI_CALLOC(fs%vk, (3, fs%maxnb))
   ABI_CALLOC(fs%vkq, (3, fs%maxnb))

   fs%krank = krank_new(nkfs, fs%kpts)
 end do ! spin

 call cwtime_report(" fstab_init%fs_build:", cpu, wall, gflops)

 ! fix window around fermie for tetrahedron or gaussian weight calculation
 ! this is spin independent
 nene = 100 ! TODO: make this variable and maybe temperature dependent???
 deltaene = 2 * dtset%eph_fsewin / dble(nene-1)
 ifermi = int(nene / 2)
 enemin = ebands%fermie - dble(ifermi-1)*deltaene
 enemax = enemin + dble(nene-1)*deltaene

 rlatt = kptrlatt
 call matr3inv(rlatt, klatt)

 ! Setup FS integration
 do spin=1,ebands%nsppol
   fs => fstab(spin)
   fs%nene = nene
   fs%enemin = enemin
   fs%deltaene = deltaene
   fs%eph_intmeth = dtset%eph_intmeth
   fs%eph_fsmear = dtset%eph_fsmear

   fs%kmesh_cartvec(:, 1) = cryst%gprimd(:,1)*klatt(1,1) + cryst%gprimd(:,2)*klatt(2,1) + cryst%gprimd(:,3)*klatt(3,1)
   fs%kmesh_cartvec(:, 2) = cryst%gprimd(:,1)*klatt(1,2) + cryst%gprimd(:,2)*klatt(2,2) + cryst%gprimd(:,3)*klatt(3,2)
   fs%kmesh_cartvec(:, 3) = cryst%gprimd(:,1)*klatt(1,3) + cryst%gprimd(:,2)*klatt(2,3) + cryst%gprimd(:,3)*klatt(3,3)
   fs%kmesh_cartvec = two_pi * fs%kmesh_cartvec
   !do i1=1,3
   !  write(std_out, *)"klatt:", klatt(:, i1)
   !  write(std_out, *)"gprimd:", cryst%gprimd(:, i1)
   !  write(std_out, *)"cartvec:", fs%kmesh_cartvec(:, i1)
   !end do
 end do

 ! TODO: compute weights on the fly to reduce memory? nene should be set to zero if not used!
 if (dtset%eph_intmeth == 2) then
   ABI_MALLOC(bz2ibz, (nkbz))
   bz2ibz = full2ebands(1, :)

   call htetra_init(tetra, bz2ibz, cryst%gprimd, klatt, kbz, nkbz, ebands%kptns, nkibz, ierr, errstr, comm)
   ABI_CHECK(ierr == 0, errstr)
   ABI_FREE(bz2ibz)

   ABI_MALLOC(tmp_eigen, (nkibz))
   ABI_MALLOC(btheta, (nene, nkibz))
   ABI_MALLOC(bdelta, (nene, nkibz))

   do spin=1,ebands%nsppol
     fs => fstab(spin)

     ! Allocate tables used to store tetrahedron weights.
     ABI_CALLOC(fs%dbldelta_tetra_weights_kfs, (fs%maxnb, fs%maxnb, fs%nkfs))

     ABI_CALLOC(fs%tetra_wtk, (fs%maxnb, nkibz))
     ABI_CALLOC(fs%tetra_wtk_ene, (fs%maxnb, nkibz, fs%nene))

     do band=fs%bmin,fs%bmax
       ! Get the contribution of this band
       tmp_eigen = ebands%eig(band, :nkibz, spin)

       ! Calculate general integration weights at each irred kpoint
       ! as in Blochl et al PRB 49 16223 [[cite:Bloechl1994a]]
       call tetra%blochl_weights(tmp_eigen, enemin, enemax, max_occ1, fs%nene, nkibz, &
         bcorr0, btheta, bdelta, xmpi_comm_self)

       ! Save weights in the correct position.
       ib = band - fs%bmin + 1
       do ik_ibz=1,nkibz
         fs%tetra_wtk(ib, ik_ibz) = bdelta(ifermi, ik_ibz) * nkibz
         fs%tetra_wtk_ene(ib, ik_ibz, 1:fs%nene) = bdelta(1:fs%nene, ik_ibz) * nkibz
       end do
     end do ! band
   end do ! spin

   ABI_FREE(tmp_eigen)
   ABI_FREE(btheta)
   ABI_FREE(bdelta)
   call tetra%free()
 end if

 !ABI_FREE(fs2ibz)
 ABI_FREE(fs2bz)
 ABI_FREE(kbz)
 ABI_FREE(full2ebands)
 ABI_FREE(indkk)

 call cwtime_report(" fstab_init%fs_weights:", cpu, wall, gflops)

end subroutine fstab_init
!!***

!----------------------------------------------------------------------

!!****f* m_fstab/fstab_findkg0
!! NAME
!!  fstab_findkg0
!!
!! FUNCTION
!!  Return the index `ikfs` of the k-point `kpt` in the FS-BZ. Return -1 if not found.
!!
!! INPUTS
!!  kpt(3)=K-point in reduced coordinates
!!
!! OUTPUT
!!   g0=Reciprocal lattice vector such that kpt = fstab%kpts(:, ikfs) + g0
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function fstab_findkg0(fstab, kpt, g0) result(ikfs)

!Arguments ------------------------------------
!scalars
 class(fstab_t),intent(in) :: fstab
!arrays
 integer,intent(out) :: g0(3)
 real(dp),intent(in) :: kpt(3)

! *************************************************************************

 ikfs = fstab%krank%get_index(kpt)
 if (ikfs /= -1) then
   g0 = nint(kpt - fstab%kpts(:, ikfs))
 else
   g0 = huge(1)
 end if

end function fstab_findkg0
!!***

!----------------------------------------------------------------------

!!****f* m_fstab/fstab_setup_qpoint
!! NAME
!! fstab_setup_qpoint
!!
!! FUNCTION
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine fstab_setup_qpoint(fs, cryst, ebands, spin, ltetra, qpt, comm)

 use libtetrabz_dbldelta_mod

!Arguments ------------------------------------
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 class(fstab_t),intent(inout) :: fs
 integer,intent(in) :: spin, ltetra, comm
 real(dp),intent(in) :: qpt(3)

!Local variables-------------------------------
!scalars
 integer,parameter :: enough = 50
 integer :: nkbz, ierr, ib, nb, ik_bz, ik_ibz, ik_fs, i1, i2, i3, band
 character(len=500) :: msg
 type(krank_t) :: krank
!arrays
 integer :: nge(3), ngw(3)
 !integer,allocatable :: symrecfm(:,:,:)
 real(dp) :: kpt(3)
 real(dp),allocatable :: eig_k(:,:), eig_kq(:,:), wght_bz(:,:,:)
 !real(dp),allocatable,intent(out) :: wght(:,:,:) !(nb*nb,PRODUCT(ngw(1:3)))

! *************************************************************************

 if (fs%eph_intmeth /= 2) return
 return

 fs%dbldelta_tetra_weights_kfs = zero

 ! The double delta is ill-defined for q == 0.
 if (fs%eph_intmeth == 2 .and. all(abs(qpt) < tol12)) then
   MSG_COMMENT("tetrahedron method with q = 0 is ill-defined. Returning")
   return
 end if

 ABI_CHECK(isdiagmat(ebands%kptrlatt), "kptrlatt must be diagonal when tetra is used.")
 ABI_CHECK(ebands%nshiftk == 1, "nshiftk must be 1 when tetra is used")
 nge = get_diag(ebands%kptrlatt)
 ngw = nge

 nb = fs%maxnb
 nkbz = product(ngw(1:3))

 ABI_MALLOC(eig_k, (nb, nkbz))
 ABI_MALLOC(eig_kq, (nb, nkbz))

 ! TODO: Handle symmetries in a cleaner way
 !timrev = 0; if (use_tr) timrev=1
 krank = krank_new(ebands%nkpt, ebands%kptns, nsym=cryst%nsym, symrec=cryst%symrec, time_reversal=.True.)

 ierr = 0
 ik_bz = 0
 do i3=0,nge(3) - 1
   do i2=0,nge(2) - 1
     do i1 =0,nge(1) - 1
       ik_bz = ik_bz + 1

       ! Find correspondence between the grid and the IBZ
       kpt = ([i1, i2, i3] + ebands%shiftk(:, 1)) / nge(:)
       !ik_ibz = krank%find(kpt, msg)
       ik_ibz = krank%invrank(krank%get_rank(kpt))

       if (ik_ibz < 1) then
         if (ierr <= enough) then
           write(msg,'(3a,i0,a)') &
            'kpt: ',trim(ktoa(kpt)),' with rank: ', krank%get_rank(kpt),' has no symmetric among the k-points'
           MSG_WARNING(msg)
         end if
         ierr = ierr + 1
         cycle
       end if

       eig_k(:, ik_bz) = ebands%eig(fs%bmin:fs%bmax, ik_ibz, spin)

       ! Find correspondence between the grid and the IBZ
       kpt = kpt + qpt
       ik_ibz = krank%invrank(krank%get_rank(kpt))

       if (ik_ibz < 1) then
         if (ierr <= enough) then
           write(msg,'(3a,i0,a)') &
            'kpt + qpt: ',trim(ktoa(kpt)),' with rank: ', krank%get_rank(kpt),' has no symmetric among the k-points'
           MSG_WARNING(msg)
         end if
         ierr = ierr + 1
         cycle
       end if

       eig_kq(:, ik_bz) = ebands%eig(fs%bmin:fs%bmax, ik_ibz, spin)
     end do
   end do
 end do

 ABI_CHECK(ierr == 0, "See above warnings")
 call krank%free()

 ! Compute weights for double delta. Note that libtetra assumes Ef set to zero.
 eig_k = eig_k - ebands%fermie
 eig_kq = eig_kq - ebands%fermie
 ABI_MALLOC(wght_bz, (nb, nb, nkbz))

 write(std_out,*)" Calling libtetrabz_dbldelta with ltetra:", ltetra
 write(std_out,*)" Q-point", qpt
 call libtetrabz_dbldelta(ltetra, cryst%gprimd, nb, nge, eig_k, eig_kq, ngw, wght_bz) !, comm=comm)

 ABI_FREE(eig_k)
 ABI_FREE(eig_kq)

 ! Convert from full BZ to fs% kpoints
 ! TODO: Average over degenerate states?
 krank = krank_new(fs%nkfs, fs%kpts)

 ik_bz = 0
 do i3=0,nge(3) - 1
    do i2=0,nge(2) - 1
       do i1 =0,nge(1) - 1
         ik_bz = ik_bz + 1
         kpt = ([i1, i2, i3] + ebands%shiftk(:, 1)) / nge(:)
         ! Then we have to rearrange the weights
         ik_fs = krank%invrank(krank%get_rank(kpt))
         if (ik_fs /= -1) then
           fs%dbldelta_tetra_weights_kfs(:,:,ik_fs) = wght_bz(:,:,ik_bz) !* fs%nktot
           !write(std_out,*)"dbldelta wts:", wght_bz(:,:,ik_bz)
         else
           !write(std_out,*)"should be zero :", wght_bz(:,:,ik_bz)
         end if
       end do
    end do
 end do

 call krank%free()
 ABI_FREE(wght_bz)

end subroutine fstab_setup_qpoint
!!***

!----------------------------------------------------------------------

!!****f* m_fstab/fstab_get_dbldelta_weights
!! NAME
!!  fstab_get_dbldelta_weights
!!
!! FUNCTION
!!  Return the weights for the integration on the Fermi-surface
!!
!! INPUTS
!!  ebands<ebands_type>=GS band structure.
!!  ik_ibz=Index of the k-point in the IBZ
!!  spin=Spin index
!!
!! OUTPUT
!!   wtk(fs%maxnb)=Weights for FS integration.
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine fstab_get_dbldelta_weights(fs, ebands, ikfs, ik_ibz, ikq_ibz, spin, wtk)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikfs, ik_ibz, ikq_ibz, spin
 class(fstab_t),intent(in) :: fs
 type(ebands_t),intent(in) :: ebands
!arrays
 real(dp),intent(out) :: wtk(fs%maxnb, fs%maxnb)

!Local variables-------------------------------
!scalars
 integer :: bstart_k, nband_k, bstart_kq, nband_kq, ib1, band1, ib2, band2, ii
 real(dp) :: g1, g2, sigma

! *************************************************************************

 bstart_k = fs%bstart_cnt_ibz(1, ik_ibz); nband_k = fs%bstart_cnt_ibz(2, ik_ibz)
 ABI_CHECK(nband_k >= 1 .and. nband_k <= fs%maxnb, "Wrong nband_k")

 bstart_kq = fs%bstart_cnt_ibz(1, ikq_ibz); nband_kq = fs%bstart_cnt_ibz(2, ikq_ibz)
 ABI_CHECK(nband_kq >= 1 .and. nband_kq <= fs%maxnb, "Wrong nband_kq")

 wtk = zero
 select case (fs%eph_intmeth)
 case (1)
   sigma = fs%eph_fsmear
   do ib2=1,nband_k
     band2 = ib2 + bstart_k - 1
     if (fs%eph_fsmear < zero) then
       sigma = max(maxval([(abs(dot_product(fs%vk(:, ib2), fs%kmesh_cartvec(:,ii))), ii=1,3)]), fs%min_smear)
     end if
     g2 = gaussian(ebands%eig(band2, ik_ibz, spin) - ebands%fermie, sigma)
     do ib1=1,nband_kq
       band1 = ib1 + bstart_kq - 1
       if (fs%eph_fsmear < zero) then
         sigma = max(maxval([(abs(dot_product(fs%vkq(:, ib1), fs%kmesh_cartvec(:,ii))), ii=1,3)]), fs%min_smear)
       end if
       g1 = gaussian(ebands%eig(band1, ikq_ibz, spin) - ebands%fermie, sigma)
       wtk(ib1, ib2) = (g1 * g2) / fs%nktot
     end do
   end do

 case (2)
   ! Copy weights in the correct position.
   do ib2=1,nband_k
     band2 = ib2 + bstart_k - fs%bmin
     do ib1=1,nband_kq
       band1 = ib1 + bstart_kq - fs%bmin
       ! This is the old version (WRONG)
       wtk(ib1, ib2) = fs%tetra_wtk(band1, ikq_ibz) * fs%tetra_wtk(band2, ik_ibz) / fs%nktot

       !write(std_out,*)wtk(ib1, ib2), fs%dbldelta_tetra_weights_kfs(band1, band2, ikfs), &
       !                abs(wtk(ib1, ib2) - fs%dbldelta_tetra_weights_kfs(band1, band2, ikfs))

       ! libtetrabz_dbldelta seems to report Weights in this order.
       !wtk(ib1, ib2) = fs%dbldelta_tetra_weights_kfs(band1, band2, ikfs)
     end do
   end do

 case default
   MSG_ERROR(sjoin("Wrong integration method:", itoa(fs%eph_intmeth)))
 end select

end subroutine fstab_get_dbldelta_weights
!!***

!----------------------------------------------------------------------

!!****f* m_fstab/fstab_print
!! NAME
!!  fstab_print
!!
!! FUNCTION
!!  Print info on the object.
!!
!! INPUTS
!! [unit]=the unit number for output
!! [prtvol]=verbosity level
!!
!! OUTPUT
!!  Only printing.
!!
!! PARENTS
!!      m_phgamma
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine fstab_print(fstab, header, unit, prtvol)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: prtvol,unit
 character(len=*),optional,intent(in) :: header
 class(fstab_t),target,intent(in) :: fstab(:)

!Local variables-------------------------------
!scalars
 integer :: my_unt,my_prtvol,spin
 type(fstab_t),pointer :: fs
 character(len=500) :: msg

! *************************************************************************

 my_unt = std_out; if (present(unit)) my_unt = unit
 my_prtvol = 0; if (present(prtvol)) my_prtvol = prtvol

 msg=' ==== Info on the fstab% object ==== '
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 write(my_unt, "(a)")trim(msg)

 if (fstab(1)%eph_intmeth == 1) then
   write(my_unt,"(a)")" FS integration done with gaussian method"
 else if (fstab(1)%eph_intmeth == 2) then
   write(my_unt,"(a)")" FS integration done with tetrahedron method"
 end if
 write(my_unt,"(a,i0)")" Total number of k-points in the full mesh: ",fstab(1)%nktot

 do spin=1,size(fstab)
   fs => fstab(spin)
   write(my_unt,"(a,i0)")" For spin: ",spin
   write(my_unt,"(a,i0,a,f5.1,a)") &
     "  Number of BZ k-points close to the Fermi surface: ",fs%nkfs," [", (100.0_dp * fs%nkfs) / fs%nktot, " %]"
   write(my_unt,"(a,i0)")"  Maximum number of bands crossing the Fermi level: ",fs%maxnb
   write(my_unt,"(2(a,i0))")"  min band: ",minval(fs%bstart_cnt_ibz(1,:), mask=fs%bstart_cnt_ibz(1,:) /= -1)
   write(my_unt,"(2(a,i0))")"  Max band: ",maxval(fs%bstart_cnt_ibz(1,:)+fs%bstart_cnt_ibz(2,:) - 1, &
                                                  mask=fs%bstart_cnt_ibz(1,:) /= -1)
 end do

end subroutine fstab_print
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
!!      elphon,get_tau_k
!!
!! CHILDREN
!!      destroy_kptrank,get_rank,mkkptrank,wrtout
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

end module m_fstab
!!***
