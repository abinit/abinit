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

 use m_time,           only : cwtime, cwtime_report
 use m_fstrings,       only : itoa, sjoin
 use m_numeric_tools,  only : bisect
 use m_symtk,          only : matr3inv
 use defs_datatypes,   only : ebands_t
 use m_crystal,        only : crystal_t
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

   integer :: nkfs
    ! Number of k-points on the Fermi-surface (FS-BZ).

   integer :: nktot
    ! Total number of k-points in the initial mesh.

   integer :: nkibz
    ! Number of points in the IBZ

   integer :: maxnb
   ! Max number of bands on the FS.
   ! TODO: Maybe maxnbfs

   ! real(dp) :: fermie
   ! Fermi energy

   integer :: integ_method
   ! Integration method. 1 for gaussian, 2 for tetrahedra

   integer :: nsig
   ! Number of smearing values used for Gaussian integration

   integer :: nene
   ! Number of chemical potential values used for inelastic integration

   real(dp) :: enemin
   ! Minimal chemical potential value used for inelastic integration

   real(dp) :: deltaene
   ! Chemical potential increment for inelastic integration

   type(krank_t) :: krank
   ! rank/inverse_rank pair for the k-points on the FS (kpts).

   integer,allocatable :: indkk_fs(:,:)
   ! indkk_fs(6, nkfs)
   ! Tables giving the correspondence between a point in the FS-BZ and the IBZ computed in listkk.
   !
   !   indkk_fs(1,:)      Mapping FS-BZ --> k-points in the IBZ (taken from ebands_t)
   !   indkk_fs(2,:)      The index of the symmetry S such that kfs = tim_sign * S(k_ibz) + G0
   !   indkk_fs(3:5,:)    The reduced components of G0.
   !   indkk_fs(6,:)      1 if time-reversal was used to generate the k-point, 0 otherwise

   integer,allocatable :: bstcnt_ibz(:,:)
    ! bstcnt_ibz(2, nkibz)
    ! The indices of the bands within the energy window (depends on fsk)
    ! Note that we use the k-point index in the IBZ.
    !
    !   bstcnt(1, :) The index of the first band inside the energy window (start)
    !   bstcnt(2, :) Number of bands on the FS (count)

   real(dp),allocatable :: kpts(:,:)
   ! kpts(3,nkfs)
   ! Reduced coordinates of the k-points on the Fermi surface.

   real(dp),allocatable :: tetra_wtk(:,:)
   ! tetra_wtk(maxnb, nkibz)
   ! Weights for FS integration with tetrahedron method
   ! Note that the weights are dimensioned with nkibz

   real(dp),allocatable :: tetra_wtk_ene(:,:,:)
   ! tetra_wtk_ene(maxnb, nkibz, nene)
   ! Weights for FS integration with tetrahedron method
   ! for all chemical potentials
   ! Note that the weights are dimensioned with nkibz

 contains

   procedure :: free => fstab_free
     ! Free memory.

   procedure :: findkg0 => fstab_findkg0
     ! Find the index of the k-point on the FS

   procedure :: get_weights_ibz => fstab_get_weights_ibz
     ! Compute weights for FS integration.

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

 !@fstab_t

 ! integer
 ABI_SFREE(fstab%indkk_fs)
 ABI_SFREE(fstab%bstcnt_ibz)

 ! real
 ABI_SFREE(fstab%kpts)
 ABI_SFREE(fstab%tetra_wtk)
 ABI_SFREE(fstab%tetra_wtk_ene)

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
!!  fsewin=Energy window in Hartree. Only states in [efermi-fsewin, efermi+fsewin] are included.
!!  integ_method=Flag selecting the integration method.
!!  kptrlatt(3,3)=k-point lattice specification
!!  nshiftk= number of shift vectors.
!!  shiftk(3,nshiftk)=shift vectors for k point generation
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

subroutine fstab_init(fstab, ebands, cryst, fsewin, integ_method, kptrlatt, nshiftk, shiftk, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nshiftk,integ_method,comm
 real(dp),intent(in) :: fsewin
 type(ebands_t),intent(in) :: ebands
 type(crystal_t),intent(in) :: cryst
!arrays
 integer,intent(in) :: kptrlatt(3,3)
 real(dp),intent(in) :: shiftk(3,nshiftk)
 type(fstab_t),target,intent(out) :: fstab(ebands%nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: option0 = 0, brav1 = 1, bcorr0 = 0
 integer :: nkfs,spin,band,nband_k,i1,i2,ib,blow,ik_bz,ik_ibz,nkibz,sppoldbl,timrev
 integer :: ik,mkpt,nkbz,ierr, bstart_k,bstop_k,nene,ifermi,bmin,bmax
 real(dp) :: elow,ehigh,ebis,enemin,enemax,deltaene,max_occ,dksqmax,cpu,wall,gflops
 logical :: inwin
 character(len=80) :: errstr
 character(len=500) :: msg
 type(fstab_t),pointer :: fs
 type(htetra_t) :: tetra
!arrays
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

 !call kpts_ibz_from_kptrlatt(cryst, kptrlatt, ebands%kptopt, nshiftk, shiftk, nkibz, kibz, wtk, nkbz, kbz, &
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

 call smpbz(brav1, std_out, kptrlatt, mkpt, nkbz, nshiftk, option0, shiftk, kbz)

 ! Find correspondence BZ --> IBZ
 ! Note that we use symrel so these tables can be used to symmetrize wavefunctions.
 timrev = kpts_timrev_from_kptopt(ebands%kptopt)
 sppoldbl = 1; if (any(cryst%symafm == -1) .and. ebands%nsppol == 1) sppoldbl=2
 ABI_MALLOC(indkk, (nkbz*sppoldbl, 6))

 call listkk(dksqmax, cryst%gmet, indkk, ebands%kptns, kbz, ebands%nkpt, nkbz, cryst%nsym,&
    sppoldbl, cryst%symafm, cryst%symrel, timrev, comm, exit_loop=.True., use_symrec=.False.)

 if (dksqmax > tol12) then
   write(msg, '(7a,es16.6,4a)' )&
   'The WFK file cannot be used to start thee present calculation ',ch10,&
   'It was asked that the wavefunctions be accurate, but',ch10,&
   'at least one of the k points could not be generated from a symmetrical one.',ch10,&
   'dksqmax= ',dksqmax,ch10,&
   'Action: check your WFK file and k-point input variables',ch10,&
   '        (e.g. kptopt or shiftk might be wrong in the present dataset or the preparatory one.'
   MSG_ERROR(msg)
 end if

 call cwtime_report("fstab_init%listkk", cpu, wall, gflops)

 ABI_MALLOC(full2ebands, (6, nkbz))
 full2ebands = 0

 do ik_bz=1,nkbz
   full2ebands(1, ik_bz) = indkk(ik_bz, 1)     ! ik_ibz
   full2ebands(2, ik_bz) = indkk(ik_bz, 2)     ! isym
   !full2ebands(3, ik_bz) = indkk(ik_bz, 6)     ! itimrev
   !full2ebands(4:6, ik_bz) = indkk(ik_bz, 3:5) ! g0
   full2ebands(3:5, ik_bz) = indkk(ik_bz, 3:5)  ! g0
   full2ebands(6, ik_bz) = indkk(ik_bz, 6)  ! itimrev
 end do

 ! Select only those k-points in the BZ close to the FS.
 ABI_CHECK(fsewin > tol12, "fsewin < tol12")
 elow = ebands%fermie - fsewin
 ehigh = ebands%fermie + fsewin
 ebis = elow - abs(elow) * 0.001_dp

 ! Allocate workspace arrays.
 !ABI_MALLOC(fs2ibz, (nkbz))
 ABI_MALLOC(fs2bz, (nkbz))

 do spin=1,ebands%nsppol
   fs => fstab(spin)
   ABI_MALLOC(fs%bstcnt_ibz, (2, nkibz))
   fs%bstcnt_ibz = -1

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
       if (any(fs%bstcnt_ibz(:, ik_ibz) /= [-1, -1])) then
         ABI_CHECK(all(fs%bstcnt_ibz(:, ik_ibz) == [i1, i2-i1+1]), "bstcnt_ibz!")
       end if
       fs%bstcnt_ibz(:, ik_ibz) = [i1, i2-i1+1]
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
   fs%maxnb = maxval(fs%bstcnt_ibz(2, :))
   fs%krank = krank_new(nkfs, fs%kpts)

 end do ! spin

 call cwtime_report("fstab_init%fs_build:", cpu, wall, gflops)

 ! fix window around fermie for tetrahedron or gaussian weight calculation
 ! this is spin independent
 nene = 100 ! TODO: make this variable and maybe temperature dependent???
 deltaene = 2*fsewin/dble(nene-1)
 ifermi = int(nene/2)
 enemin = ebands%fermie - dble(ifermi-1)*deltaene
 enemax = enemin + dble(nene-1)*deltaene

 ! Setup FS integration
 do spin=1,ebands%nsppol
   fs => fstab(spin)
   fs%nene = nene
   fs%enemin = enemin
   fs%deltaene = deltaene
   fs%nsig = 1
   fs%integ_method = integ_method
 end do

 ! TODO: compute weights on the fly to reduce memory? nene should be set to zero if not used!
 if (integ_method == 2) then
   rlatt = kptrlatt
   call matr3inv(rlatt,klatt)

   ABI_MALLOC(bz2ibz, (nkbz))
   bz2ibz = full2ebands(1, :)

   call htetra_init(tetra, bz2ibz, cryst%gprimd, klatt, kbz, nkbz, ebands%kptns, nkibz, ierr, errstr, comm)
   ABI_CHECK(ierr == 0, errstr)
   ABI_FREE(bz2ibz)

   ABI_MALLOC(tmp_eigen, (nkibz))
   ABI_MALLOC(btheta, (nene, nkibz))
   ABI_MALLOC(bdelta, (nene, nkibz))

   max_occ = one
   ! in spinor or spin polarized case, orbitals have occupation <= 1 instead of 2
   !if (ebands%nsppol > 1) max_occ = one
   ! this accounts for the doubling of the num of bands, even though spin channels are not well defined
   !if (ebands%nspinor == 2) max_occ = half

   do spin=1,ebands%nsppol
     fs => fstab(spin)
     ABI_CALLOC(fs%tetra_wtk, (fs%maxnb, nkibz))
     ABI_CALLOC(fs%tetra_wtk_ene, (fs%maxnb, nkibz, fs%nene))

     ! we have to pass full bands to get_tetra_weight
     ! Here I create a full set of bands enclosing the states that will be used in the FS integration.
     ! Then we have to rearrange the weights
     bmin = huge(1); bmax = -huge(1)
     do ik_ibz=1,nkibz
       if (fs%bstcnt_ibz(1,ik_ibz) /= -1) then
         bmin = min(bmin, fs%bstcnt_ibz(1,ik_ibz))
       end if
       if (fs%bstcnt_ibz(2,ik_ibz) /= -1) then
         bmax = max(bmax, fs%bstcnt_ibz(1,ik_ibz) + fs%bstcnt_ibz(2,ik_ibz) - 1)
       end if
     end do
     !write(std_out,*)"bmin, bmax for tetra: ",bmin,bmax
     ABI_CHECK(bmin /= huge(1) .and. bmax /= -huge(1), "No point on the Fermi surface!")

     !call libtetrabz_dbldelta(ltetra, bvec, nb, nge, eig1, eig2, ngw, wght_bz, comm)

     do band=bmin,bmax
       ! Get the contribution of this band
       tmp_eigen = ebands%eig(band, :nkibz, spin)

       ! Calculate general integration weights at each irred kpoint
       ! as in Blochl et al PRB 49 16223 [[cite:Bloechl1994a]]
       call tetra%blochl_weights(tmp_eigen, enemin, enemax, max_occ, fs%nene, nkibz, &
         bcorr0, btheta, bdelta, xmpi_comm_self)

       do ik_ibz=1,nkibz
         bstart_k = fs%bstcnt_ibz(1, ik_ibz); bstop_k = bstart_k + fs%bstcnt_ibz(2, ik_ibz) - 1
         if (band >= bstart_k .and. band <= bstop_k) then
           ! Save weights in the correct position.
           ib = band - bstart_k + 1
           fs%tetra_wtk(ib, ik_ibz) = bdelta(ifermi, ik_ibz) * nkibz
           fs%tetra_wtk_ene(ib, ik_ibz, 1:fs%nene) = bdelta(1:fs%nene, ik_ibz) * nkibz
         end if
       end do ! ik_ibz
     end do ! band
   end do ! sppol

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

 call cwtime_report("fstab_init%fs_weights:", cpu, wall, gflops)

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

!!****f* m_fstab/fstab_get_weights_ibz
!! NAME
!!  fstab_get_weights_ibz
!!
!! FUNCTION
!!  Return the weights for the integration on the Fermi-surface
!!
!! INPUTS
!!  ebands<ebands_type>=GS band structure.
!!  ik_ibz=Index of the k-point in the IBZ
!!  spin=Spin index
!!  sigmas
!!  [iene]
!!
!! OUTPUT
!!   wtk(fs%nsig, fs%maxnb)=Weights for FS integration.
!!
!! PARENTS
!!      m_ddk,m_phgamma
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine fstab_get_weights_ibz(fs, ebands, ik_ibz, spin, sigmas, wtk, iene)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,spin
 integer,intent(in),optional :: iene
 class(fstab_t),intent(in) :: fs
 type(ebands_t),intent(in) :: ebands
!arrays
 real(dp),intent(in) :: sigmas(:) !fs%nsig)
 real(dp),intent(out) :: wtk(fs%nsig,fs%maxnb)

!Local variables-------------------------------
!scalars
 integer :: ib, bstart_k, nband_k, band, isig
 real(dp) :: arg, chempot

! *************************************************************************

 bstart_k = fs%bstcnt_ibz(1, ik_ibz); nband_k = fs%bstcnt_ibz(2, ik_ibz)
 ABI_CHECK(nband_k >= 1 .and. nband_k <= fs%maxnb, "wrong nband_k")

 ! TODO: add iene looping for chemical potential in gaussian case too
 select case (fs%integ_method)
 case (1)
   chempot = ebands%fermie; if (present(iene)) chempot = fs%enemin + (iene-1)*fs%deltaene
   do ib=1,nband_k
     band = ib + bstart_k - 1
     arg = ebands%eig(band,ik_ibz,spin) - chempot
     do isig=1,fs%nsig
       wtk(isig,ib) = gaussian(arg, sigmas(isig))
     end do
   end do

 case (2)
   if (present(iene)) then
     wtk(1,1:nband_k) = fs%tetra_wtk_ene(1:nband_k, ik_ibz, iene)
   else
     wtk(1,1:nband_k) = fs%tetra_wtk(1:nband_k, ik_ibz)
   end if

 case default
   MSG_ERROR(sjoin("Wrong integration method:", itoa(fs%integ_method)))
 end select

end subroutine fstab_get_weights_ibz
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
!! [mode_paral]=either "COLL" or "PERS"
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

subroutine fstab_print(fstab, header, unit, prtvol, mode_paral)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: prtvol,unit
 character(len=4),optional,intent(in) :: mode_paral
 character(len=*),optional,intent(in) :: header
 class(fstab_t),target,intent(in) :: fstab(:)

!Local variables-------------------------------
!scalars
 integer :: my_unt,my_prtvol,spin
 type(fstab_t),pointer :: fs
 character(len=4) :: my_mode
 character(len=500) :: msg

! *************************************************************************

 my_unt =std_out; if (present(unit)) my_unt = unit
 my_prtvol=0    ; if (present(prtvol)) my_prtvol = prtvol
 my_mode='COLL' ; if (present(mode_paral)) my_mode = mode_paral

 msg=' ==== Info on the fstab% object ==== '
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 if (fstab(1)%integ_method == 1) then
   write(std_out,"(a,i0)")"FS integration done with gaussian method and nsig: ",fstab(1)%nsig
 else if (fstab(1)%integ_method == 2) then
   write(std_out,"(a)")"FS integration done with tetrahedron method"
 end if
 write(std_out,"(a,i0)")"Total number of points in the full mesh: ",fstab(1)%nktot

 do spin=1,size(fstab)
   fs => fstab(spin)
   write(std_out,"(a,i0)")"For spin: ",spin
   write(std_out,"(a,i0,a,f5.1,a)")&
     "  Number of BZ k-points close to the Fermi surface: ",fs%nkfs," [",(100.0_dp*fs%nkfs)/fs%nktot," %]"
   write(std_out,"(a,i0)")"  Maximum number of bands crossing the Fermi level: ",fs%maxnb
   write(std_out,"(2(a,i0))")"  min band: ",minval(fs%bstcnt_ibz(1,:), mask=fs%bstcnt_ibz(1,:)/=-1)
   write(std_out,"(2(a,i0))")"  Max band: ",maxval(fs%bstcnt_ibz(1,:)+fs%bstcnt_ibz(2,:)-1, mask=fs%bstcnt_ibz(1,:)/=-1)
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
