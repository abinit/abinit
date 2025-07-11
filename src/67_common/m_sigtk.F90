!!****m* ABINIT/m_sigtk
!! NAME
!!  m_sigtk
!!
!! FUNCTION
!!  Helper functions common to electron self-energy calculations. Provides tools to:
!!
!!      - Define list of k-points and bands in sel-energy matrix elements from input variables.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2025 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_sigtk

 use defs_basis
 use m_abicore
 use m_errors
 use m_ebands
 use m_crystal
 use m_xmpi
 use netcdf
 use m_nctk
 use m_hdr
 use m_dtset
 use m_krank

 use m_build_info,   only : abinit_version
 use m_fstrings,     only : sjoin, ltoa, strcat, itoa, ftoa
 use m_io_tools,     only : open_file
 use defs_datatypes, only : pseudopotential_type
 use defs_wvltypes,  only : wvl_internal_type
 use m_pawtab,       only : pawtab_type
 use m_kpts,         only : kpts_ibz_from_kptrlatt, kpts_timrev_from_kptopt, kpts_map

 implicit none

 private

 public :: sigtk_kcalc_from_nkptgw
 public :: sigtk_kcalc_from_qprange
 public :: sigtk_kcalc_from_gaps
 public :: sigtk_kcalc_from_erange
 public :: sigtk_kpts_in_erange
 public :: sigtk_sigma_tables
!!***


 ! Tables for degenerated KS states.
 type, public :: bids_t
   integer, allocatable :: vals(:)
 end type bids_t

 type, public :: degtab_t
   type(bids_t), allocatable :: bids(:)
 end type degtab_t

 public :: degtab_array_free   ! Free array of degtab_t objects.
!!***

contains  !=====================================================
!!***

!!****f* m_sigtk/sigtk_kcalc_from_nkptgw
!! NAME
!!  sigtk_kcalc_from_nkptgw
!!
!! FUNCTION
!!  Initialize list of k-points and bands for self-energy matrix elements from nkptgw.
!!
!! INPUT
!!  dtset<dataset_type>=All input variables for this dataset.
!!  mband: Max number of bands.
!!
!! OUTPUT
!!  nkcalc: Number of k-points in self-energy matrix elements.
!!  kcalc(3, nkcalc): List of k-points where the self-energy is computed.
!!  bstart_ks(nkcalc, nsppol): Initial KS band index included in self-energy matrix elements for each k-point in kcalc.
!!  nbcalc_ks(nkcalc, nsppol): Number of bands included in self-energy matrix elements for each k-point in kcalc.
!!
!! SOURCE

subroutine sigtk_kcalc_from_nkptgw(dtset, mband, nkcalc, kcalc, bstart_ks, nbcalc_ks)

!Arguments ------------------------------------
 type(dataset_type),intent(in) :: dtset
 integer,intent(in) :: mband
 integer,intent(out) :: nkcalc
!arrays
 real(dp),allocatable,intent(out) :: kcalc(:,:)
 integer,allocatable,intent(out) :: bstart_ks(:,:)
 integer,allocatable,intent(out) :: nbcalc_ks(:,:)

!Local variables ------------------------------
!scalars
 integer :: spin, ierr, ikcalc
 character(len=500) :: msg

! *************************************************************************

 call wrtout(std_out, " Generating list of k-points for self-energy from kptgw and bdgw.")

 nkcalc = dtset%nkptgw
 ABI_MALLOC(kcalc, (3, nkcalc))
 ABI_MALLOC(bstart_ks, (nkcalc, dtset%nsppol))
 ABI_MALLOC(nbcalc_ks, (nkcalc, dtset%nsppol))

 kcalc = dtset%kptgw(:,1:nkcalc)
 do spin=1,dtset%nsppol
   bstart_ks(:,spin) = dtset%bdgw(1,1:nkcalc,spin)
   nbcalc_ks(:,spin) = dtset%bdgw(2,1:nkcalc,spin) - dtset%bdgw(1,1:nkcalc,spin) + 1
 end do

 ! Consistency check on bdgw and mband
 ierr = 0
 do spin=1,dtset%nsppol
   do ikcalc=1,nkcalc
     if (dtset%bdgw(2,ikcalc,spin) > mband) then
       ierr = ierr + 1
       write(msg,'(a,2(i0,1x),2(a,i0))')&
        "For (k, s) ",ikcalc,spin," bdgw= ",dtset%bdgw(2,ikcalc,spin), " > mband = ",mband
       ABI_WARNING(msg)
     end if
   end do
 end do
 ABI_CHECK(ierr == 0, "Not enough bands in WFK file. See messages above. Aborting now.")

end subroutine sigtk_kcalc_from_nkptgw
!!***

!!****f* m_sigtk/sigtk_kcalc_from_qprange
!! NAME
!!  sigtk_kcalc_from_qprange
!!
!! FUNCTION
!! Use qprange to select the interesting k-points and the corresponding bands.
!!
!!    0 --> Compute the QP corrections only for the fundamental and the direct gap.
!! +num --> Compute the QP corrections for all the k-points in the irreducible zone and include `num`
!!           bands above and below the Fermi level.
!! -num --> Compute the QP corrections for all the k-points in the irreducible zone.
!!          Include all occupied states and `num` empty states.
!!
!! INPUT
!!  dtset<dataset_type>=All input variables for this dataset.
!!  cryst<crystal_t>=Crystalline structure
!!  ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!!  qprange: See above description
!!
!! OUTPUT
!!  nkcalc: Number of k-points in self-energy matrix elements.
!!  kcalc(3, nkcalc): List of k-points where the self-energy is computed.
!!  bstart_ks(nkcalc, nsppol): Initial KS band index included in self-energy matrix elements for each k-point in kcalc.
!!  nbcalc_ks(nkcalc, nsppol): Number of bands included in self-energy matrix elements for each k-point in kcalc.
!!
!! SOURCE

subroutine sigtk_kcalc_from_qprange(dtset, cryst, ebands, qprange, nkcalc, kcalc, bstart_ks, nbcalc_ks)

!Arguments ------------------------------------
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 integer,intent(in) :: qprange
 integer,intent(out) :: nkcalc
!arrays
 real(dp),allocatable,intent(out) :: kcalc(:,:)
 integer,allocatable,intent(out) :: bstart_ks(:,:), nbcalc_ks(:,:)

!Local variables ------------------------------
!scalars
 integer :: spin, ik, bstop, mband, sigma_nkbz
!arrays
 integer :: kptrlatt(3,3)
 integer :: val_indices(ebands%nkpt, ebands%nsppol)
 real(dp),allocatable :: sigma_wtk(:),sigma_kbz(:,:)

! *************************************************************************

 mband = ebands%mband

 val_indices = ebands%get_valence_idx()

 if (any(dtset%sigma_ngkpt /= 0)) then
    call wrtout(std_out, " Generating list of k-points for self-energy from sigma_ngkpt and qprange.")
    ABI_CHECK(qprange /= 0, "qprange must be != 0")
    ! Get %kcalc from sigma_ngkpt
    kptrlatt = 0
    kptrlatt(1,1) = dtset%sigma_ngkpt(1); kptrlatt(2,2) = dtset%sigma_ngkpt(2); kptrlatt(3,3) = dtset%sigma_ngkpt(3)
    call kpts_ibz_from_kptrlatt(cryst, kptrlatt, dtset%kptopt, dtset%sigma_nshiftk, dtset%sigma_shiftk, &
                                nkcalc, kcalc, sigma_wtk, sigma_nkbz, sigma_kbz)
    ABI_FREE(sigma_kbz)
    ABI_FREE(sigma_wtk)
 else
    ! Include all the k-points in the IBZ.
    ! Note that kcalc == ebands%kptns so we can use a single ik index in the loop over k-points.
    ! No need to map kcalc onto ebands%kptns.
    call wrtout(std_out, " nkptgw set to 0 ==> Include all k-points in the IBZ for Sigma_nk.")
    nkcalc = ebands%nkpt
    ABI_MALLOC(kcalc, (3, nkcalc))
    kcalc = ebands%kptns
 end if

 ABI_MALLOC(bstart_ks, (nkcalc, dtset%nsppol))
 ABI_MALLOC(nbcalc_ks, (nkcalc, dtset%nsppol))

 if (qprange > 0) then
   call wrtout(std_out, " Using buffer of bands above and below the Fermi level.")
   do spin=1,dtset%nsppol
     do ik=1,nkcalc
       bstart_ks(ik,spin) = max(val_indices(ik,spin) - qprange, 1)
       bstop = min(val_indices(ik,spin) + qprange, mband)
       nbcalc_ks(ik,spin) = bstop - bstart_ks(ik,spin) + 1
     end do
   end do

 else
   call wrtout(std_out, " Including all occupied states and -qprange empty states.")
   bstart_ks = 1
   do spin=1,dtset%nsppol
     do ik=1,nkcalc
       nbcalc_ks(ik,spin) = min(val_indices(ik,spin) - qprange, mband)
     end do
   end do
 end if

end subroutine sigtk_kcalc_from_qprange
!!***

!!****f* m_sigtk/sigtk_kcalc_from_gaps
!! NAME
!!  sigtk_kcalc_from_gaps
!!
!! FUNCTION
!!  Select list of k-points and bands for self-energy matrix elements so that fundamental and direct gaps are included.
!!
!! INPUT
!!  dtset<dataset_type>=All input variables for this dataset.
!!  ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!!  gaps<gaps_t>=Store location of gaps.
!!
!! OUTPUT
!!  nkcalc: Number of k-points in self-energy matrix elements.
!!  kcalc(3, nkcalc): List of k-points where the self-energy is computed.
!!  bstart_ks(nkcalc, nsppol): Initial KS band index included in self-energy matrix elements for each k-point in kcalc.
!!  nbcalc_ks(nkcalc, nsppol): Number of bands included in self-energy matrix elements for each k-point in kcalc.
!!
!! SOURCE

subroutine sigtk_kcalc_from_gaps(dtset, ebands, gaps, nkcalc, kcalc, bstart_ks, nbcalc_ks)

!Arguments ------------------------------------
 type(dataset_type),intent(in) :: dtset
 type(ebands_t),intent(in) :: ebands
 type(gaps_t) :: gaps
 integer,intent(out) :: nkcalc
!arrays
 real(dp),allocatable,intent(out) :: kcalc(:,:)
 integer,allocatable,intent(out) :: bstart_ks(:,:)
 integer,allocatable,intent(out) :: nbcalc_ks(:,:)

!Local variables ------------------------------
!scalars
 integer :: spin, nsppol, ii, ik, nk_found, ifo, jj
 logical :: found
!arrays
 integer :: val_indices(ebands%nkpt, ebands%nsppol)
 integer :: kpos(6)

! *************************************************************************

 ABI_UNUSED((/dtset%natom/))

 call wrtout(std_out, " Including direct and fundamental KS gap in Sigma_nk")
 ABI_CHECK(maxval(gaps%ierr) == 0, "qprange 0 cannot be used because I cannot find the gap (gap_err !=0)")

 nsppol = ebands%nsppol
 val_indices = ebands%get_valence_idx()

 ! Include the direct and the fundamental KS gap.
 ! The problem here is that kptgw and nkptgw do not depend on the spin and therefore
 ! we have compute the union of the k-points where the fundamental and the direct gaps are located.
 nk_found = 1; kpos(1) = gaps%fo_kpos(1,1)

 ! Find the list of `interesting` kpoints.
 do spin=1,nsppol
   do ifo=1,3
     ik = gaps%fo_kpos(ifo, spin)
     found = .False.; jj = 0
     do while (.not. found .and. jj < nk_found)
       jj = jj + 1; found = (kpos(jj) == ik)
     end do
     if (.not. found) then
       nk_found = nk_found + 1; kpos(nk_found) = ik
     end if
   end do
 end do

 ! Now we can define the list of k-points and the bands range.
 nkcalc = nk_found
 ABI_MALLOC(kcalc, (3, nkcalc))
 ABI_MALLOC(bstart_ks, (nkcalc, nsppol))
 ABI_MALLOC(nbcalc_ks, (nkcalc, nsppol))

 do ii=1,nkcalc
   ik = kpos(ii)
   kcalc(:,ii) = ebands%kptns(:,ik)
   do spin=1,nsppol
     bstart_ks(ii,spin) = val_indices(ik,spin)
     nbcalc_ks(ii,spin) = 2
   end do
 end do

end subroutine sigtk_kcalc_from_gaps
!!***

!!****f* m_sigtk/sigtk_kcalc_from_erange
!! NAME
!!  sigtk_kcalc_from_erange
!!
!! FUNCTION
!!  Select list of k-points and bands for self-energy matrix elements on the basis of their positions
!!  wrt to the (band edges|fermi level) and the value of sigma_erange.
!!  Useful when computing electron-lifetimes for transport calculations.
!!
!! INPUT
!!  dtset<dataset_type>=All input variables for this dataset.
!!  cryst<crystal_t>=Crystalline structure
!!  ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!!  gaps<gaps_t>=Store location of gaps.
!!  comm: MPI communicator.
!!
!! OUTPUT
!!  nkcalc: Number of k-points in self-energy matrix elements.
!!  kcalc(3, nkcalc): List of k-points where the self-energy is computed.
!!  bstart_ks(nkcalc, nsppol): Initial KS band index included in self-energy matrix elements for each k-point in kcalc.
!!  nbcalc_ks(nkcalc, nsppol): Number of bands included in self-energy matrix elements for each k-point in kcalc.
!!
!! SOURCE

subroutine sigtk_kcalc_from_erange(dtset, cryst, ebands, gaps, nkcalc, kcalc, bstart_ks, nbcalc_ks, comm)

!Arguments ------------------------------------
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 type(gaps_t),intent(in) :: gaps
 integer,intent(in) :: comm
 integer,intent(out) :: nkcalc
!arrays
 real(dp),allocatable,intent(out) :: kcalc(:,:)
 integer,allocatable,intent(out) :: bstart_ks(:,:)
 integer,allocatable,intent(out) :: nbcalc_ks(:,:)

!Local variables ------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: spin, ik, band, ii, ic, nsppol, tmp_nkpt, sigma_nkbz, my_rank
 logical :: found
 real(dp) :: cmin, vmax, ee
 logical :: assume_gap
 character(len=500) :: msg
 type(krank_t) :: krank
!arrays
 integer :: kptrlatt(3,3), units(1)
 integer,allocatable :: ib_work(:,:,:), sigmak2ebands(:), indkk(:,:)
 integer :: kpos(ebands%nkpt)
 real(dp),allocatable :: sigma_wtk(:),sigma_kbz(:,:),tmp_kcalc(:,:)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm) !; nprocs = xmpi_comm_size(comm)
 units = [std_out]
 assume_gap = .not. all(dtset%sigma_erange < zero)

 if (my_rank == master) then
   write(std_out, "(a)")" Selecting k-points and bands according to their position wrt the band edges (sigma_erange)."
   write(std_out, "(a, 2(f6.3, 1x), a)")" sigma_erange: ", dtset%sigma_erange(:) * Ha_eV, " (eV)"
   if (assume_gap) then
     call gaps%print([std_out])
     ABI_CHECK(maxval(gaps%ierr) == 0, "sigma_erange 0 cannot be used because I cannot find the gap (gap_err !=0)")
   end if
 end if

 if (any(dtset%sigma_ngkpt /= 0)) then
    call wrtout(std_out, sjoin(" Generating initial list of k-points from sigma_nkpt:", ltoa(dtset%sigma_ngkpt)))
    ! Get tentative tmp_nkpt and tmp_kcalc from sigma_ngkpt.
    kptrlatt = 0
    kptrlatt(1,1) = dtset%sigma_ngkpt(1); kptrlatt(2,2) = dtset%sigma_ngkpt(2); kptrlatt(3,3) = dtset%sigma_ngkpt(3)
    call kpts_ibz_from_kptrlatt(cryst, kptrlatt, dtset%kptopt, dtset%sigma_nshiftk, dtset%sigma_shiftk, &
                                tmp_nkpt, tmp_kcalc, sigma_wtk, sigma_nkbz, sigma_kbz)

    ABI_FREE(sigma_kbz)
    ABI_FREE(sigma_wtk)

    ! Map tmp_kcalc to ebands%kpts

    ABI_MALLOC(indkk, (6, tmp_nkpt))

    krank = krank_from_kptrlatt(ebands%nkpt, ebands%kptns, ebands%kptrlatt, compute_invrank=.False.)

    if (kpts_map("symrec", ebands%kptopt, cryst, krank, tmp_nkpt, tmp_kcalc, indkk) /= 0) then
      write(msg, '(3a)' )&
        "At least one of the k-points could not be generated from a symmetrical one in the WFK.",ch10,&
        'Action: check your WFK file and the value of sigma_nkpt, sigma_shiftk in the input file.'
      ABI_ERROR(msg)
    end if

    call krank%free()

    ABI_MALLOC(sigmak2ebands, (tmp_nkpt))
    sigmak2ebands = indkk(1, :)
    ABI_FREE(tmp_kcalc)
    ABI_FREE(indkk)

 else
   ! Include all the k-points in the IBZ in the initial list.
   call wrtout(std_out, " Generating initial list of k-points from input ebands%kptns.")
   tmp_nkpt = ebands%nkpt
   ! Trivial map
   ABI_MALLOC(sigmak2ebands, (tmp_nkpt))
   sigmak2ebands = [(ii, ii=1, ebands%nkpt)]
 end if

 nsppol = ebands%nsppol
 ABI_MALLOC(ib_work, (2, tmp_nkpt, nsppol))

 do spin=1,nsppol

   if (assume_gap) then
     ! Get CBM and VBM with some tolerance
     vmax = gaps%vb_max(spin) + tol2 * eV_Ha
     cmin = gaps%cb_min(spin) - tol2 * eV_Ha
   else
     vmax = ebands%fermie
     cmin = ebands%fermie
   end if

   do ii=1,tmp_nkpt
     ! Index of k-point in ebands.
     ik = sigmak2ebands(ii)
     ! Will use this initial values to understand if k-point is in energy window.
     ib_work(1, ii, spin) = huge(1)
     ib_work(2, ii, spin) = -huge(1)
     do band=1,ebands%nband(ik + (spin-1) * ebands%nkpt)
        ee = ebands%eig(band, ik, spin)
        if (abs(dtset%sigma_erange(1)) > zero) then
          if (ee <= vmax .and. vmax - ee <= abs(dtset%sigma_erange(1))) then
            ib_work(1, ii, spin) = min(ib_work(1, ii, spin), band)
            ib_work(2, ii, spin) = max(ib_work(2, ii, spin), band)
            !write(std_out, *), "Adding valence band", band, " with ee [eV]: ", ee * Ha_eV
          end if
        end if
        if (abs(dtset%sigma_erange(2)) > zero) then
          if (ee >= cmin .and. ee - cmin <= abs(dtset%sigma_erange(2))) then
            ib_work(1, ii, spin) = min(ib_work(1, ii, spin), band)
            ib_work(2, ii, spin) = max(ib_work(2, ii, spin), band)
            !write(std_out, *)"Adding conduction band", band, " with ee [eV]: ", ee * Ha_eV
          end if
        end if
     end do
   end do
 end do

 ! Now we can define the list of k-points and the bands range.
 ! The main problem here is that kptgw and nkptgw do not depend on the spin and therefore
 ! we have to compute the union of the k-points.
 nkcalc = 0
 do ii=1,tmp_nkpt
   found = .False.
   do spin=1,nsppol
      if (ib_work(1, ii, spin) <= ib_work(2, ii, spin)) then
        found = .True.; exit
      end if
   end do
   if (found) then
     nkcalc = nkcalc + 1
     kpos(nkcalc) = ii
   end if
 end do

 ABI_MALLOC(kcalc, (3, nkcalc))
 ABI_MALLOC(bstart_ks, (nkcalc, nsppol))
 ABI_MALLOC(nbcalc_ks, (nkcalc, nsppol))

 do ic=1,nkcalc
   ! Index in the ib_work array
   ii = kpos(ic)
   ! Index in ebands.
   ik = sigmak2ebands(ii)
   kcalc(:,ic) = ebands%kptns(:,ik)
   do spin=1,nsppol
     bstart_ks(ic, spin) = 0
     nbcalc_ks(ic, spin) = 0
     if (ib_work(1, ii, spin) <= ib_work(2, ii, spin)) then
       bstart_ks(ic, spin) = ib_work(1, ii, spin)
       nbcalc_ks(ic, spin) = ib_work(2, ii, spin) - ib_work(1, ii, spin) + 1
     end if
     if (nbcalc_ks(ic, spin) == 0) then
       ABI_WARNING("Spin-polarized case with nbcalc_ks == 0, don't know if code can handle it!")
     end if
   end do
 end do

 if (my_rank == master) then
   ! Write info about k-points used in the calculation.
   write(msg, "(a, i0, a, 2(f6.3, 1x), a)") &
     " Found ", nkcalc, " k-points within sigma_erange: ", dtset%sigma_erange(:) * Ha_eV, " (eV)"
   call wrtout(units, msg)
   if (any(dtset%sigma_ngkpt /= 0)) then
     call wrtout(units, sjoin(" These k-points belong to the sigma_ngkpt k-mesh:", ltoa(dtset%sigma_ngkpt)))
   end if
   write(msg, "(2(a, i0))")" min(nbcalc_ks): ", minval(nbcalc_ks), " Max(nbcalc_ks): ", maxval(nbcalc_ks)
   call wrtout(units, msg)
 end if

 ABI_FREE(ib_work)
 ABI_FREE(sigmak2ebands)

end subroutine sigtk_kcalc_from_erange
!!***

!!****f* m_sigmaph/sigtk_kpts_in_erange
!! NAME
!!  sigtk_kpts_in_erange
!!
!! FUNCTION
!!  Use star functions interpolation and [[einterp]] to interpolate KS energies onto dense k-mesh
!!  defined by [[sigma_ngkpt]] and [[sigma_shiftk]].
!!  find k-points inside (electron/hole) pockets according to the values specifed by [[sigma_erange]].
!!  write kerange.nc file with the tables required by abinit to automate nscf band structure calculations
!!  mainly used to prepare eph calculations in which only selected k-points are nededed (imaginary part of self-energies).
!!
!! INPUTS
!!  dtset <dataset_type>=all input variables for this dataset
!!  cryst<crystal_t>=Crystalline structure
!!  ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!!  psps <pseudopotential_type>=all the information about psps
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  prefix=Prefix for output file.
!!  comm: MPI communicator.
!!
!! SOURCE

subroutine sigtk_kpts_in_erange(dtset, cryst, ebands, psps, pawtab, prefix, comm)

!Arguments ------------------------------------
!scalars
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 type(pseudopotential_type),intent(in) :: psps
 character(len=*),intent(in) :: prefix
 integer,intent(in) :: comm
!arrays
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: master = 0, pertcase0 = 0, image1 = 1
 integer :: ii, my_rank, nprocs, spin, ikf_ibz, band, nkpt_inerange, gap_err, unt, ncid, cnt, ncerr
 logical :: assume_gap
 real(dp) :: ee, cmin, vmax
 character(len=500) :: msg
 character(len=fnlen) :: path
 type(ebands_t) :: fine_ebands
 type(gaps_t) :: gaps, fine_gaps
 type(wvl_internal_type) :: dummy_wvl
 type(hdr_type) :: fine_hdr
!arrays
 integer :: fine_kptrlatt(3,3), band_block(2), units(2)
 integer,allocatable :: kshe_mask(:,:,:), krange2ibz(:)
 real(dp) :: params(4)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 ! (-num, -num) activate treatment of metals with energy window around Efermi.
 assume_gap = .not. all(dtset%sigma_erange < zero)
 units = [std_out, ab_out]

 if (my_rank == master) then
   call wrtout(units, sjoin(ch10, repeat("=", 92)))
   call wrtout(units, " Using SKW interpolation to interpolate KS energies onto dense k-mesh.")
   call wrtout(units, sjoin(" defined by sigma_ngkpt:", trim(ltoa(dtset%sigma_ngkpt))))
   ABI_CHECK(allocated(dtset%sigma_shiftk), "sigma_nshiftk must be specified in input.")
   write(std_out, "(2a)") " and sigma_shiftk shifts:"
   do ii=1,dtset%nshiftk
     call wrtout(units, sjoin(itoa(ii), ltoa(dtset%sigma_shiftk(:, ii))))
   end do

   if (assume_gap) then
     call wrtout(units, " Finding k-points inside (electron/hole) pockets (assuming semiconductor).")
   else
     call wrtout(units, " Finding k-points inside energy window around Fermi level (assuming metal).")
   end if
   write(msg, "(a, 2(f6.3, 1x), a)")" Using sigma_erange: ", dtset%sigma_erange(:) * Ha_eV, " (eV)"
   call wrtout(units, msg)
   call wrtout(units, sjoin(" SKW parameters (einterp): ", ltoa(dtset%einterp)))
   call wrtout(units, sjoin(repeat("=", 92), ch10))
   !call ebands%print([std_out], header, prtvol=dtset%prtvol)

   ! Consistency check.
   if (all(dtset%sigma_erange == zero)) then
     ABI_ERROR("sigma_erange must be specified in input when calling sigtk_kpts_in_erange.")
   end if
   if (all(dtset%sigma_ngkpt == 0)) then
     ABI_ERROR("sigma_ngkpt must be specified in input when calling sigtk_kpts_in_erange.")
   end if
 end if

 if (assume_gap) then
   ! Compute gaps using input ebands.
   gaps = ebands%get_gaps(gap_err)
   if (gap_err /= 0) then
     ABI_ERROR("Cannot compute fundamental and direct gap (likely metal).")
   end if

   if (my_rank == master) call gaps%print(units, header="Gaps from input WFK")
   call gaps%free()
 else
   call wrtout(units, sjoin("Using Fermi level:", ftoa(ebands%fermie * Ha_eV, fmt="f6.2"), " (eV)"))
 end if

 ! Interpolate band energies with star functions.
 ! In the EPH code, we will need eigens in the IBZ to compute efermi not just energies inside pockets.
 fine_kptrlatt = 0
 do ii=1,3
   fine_kptrlatt(ii, ii) = dtset%sigma_ngkpt(ii)
 end do
 band_block = [1, ebands%mband]
 params = 0; params(1) = 1; params(2) = 5; if (nint(dtset%einterp(1)) == 1) params = dtset%einterp

 fine_ebands = ebands%interp_kmesh(cryst, params, fine_kptrlatt, &
                                   dtset%sigma_nshiftk, dtset%sigma_shiftk, band_block, comm)
 fine_ebands%istwfk = 1

 call fine_ebands%update_occ(dtset%spinmagntarget, prtvol=dtset%prtvol)
 call fine_ebands%print([std_out], header="FINE EBANDS", prtvol=dtset%prtvol)

 if (assume_gap) then
   ! Compute gaps using fine_ebands.
   fine_gaps = fine_ebands%get_gaps(gap_err)
   if (gap_err /= 0) then
     ABI_ERROR("Cannot compute fundamental and direct gap (likely metal).")
   end if

   if (my_rank == master) call fine_gaps%print(units, header="Gaps from SKW interpolated eigenvalues")
 end if

 ! Build new header with fine k-mesh (note kptrlatt_orig == kptrlatt)
 call fine_hdr%init_lowlvl(fine_ebands, psps, pawtab, dummy_wvl, abinit_version, pertcase0, &
   dtset%natom, dtset%nsym, dtset%nspden, dtset%ecut, dtset%pawecutdg, dtset%ecutsm, dtset%dilatmx, &
   dtset%intxc, dtset%ixc, dtset%stmbias, dtset%usewvl, dtset%pawcpxocc, dtset%pawspnorb, dtset%ngfft, dtset%ngfftdg, &
   dtset%so_psp, dtset%qptn, cryst%rprimd, cryst%xred, cryst%symrel, cryst%tnons, cryst%symafm, cryst%typat, &
   dtset%amu_orig(:, image1), dtset%icoulomb, &
   dtset%kptopt, dtset%nelect, dtset%ne_qFD, dtset%nh_qFD, dtset%ivalence, dtset%cellcharge(1), &
   fine_kptrlatt, fine_kptrlatt, dtset%sigma_nshiftk, dtset%sigma_nshiftk, dtset%sigma_shiftk, dtset%sigma_shiftk)

 ! Find k-points inside sigma_erange energy window.
 ! Set entry to the number of states inside the pocket at (ikpt, spin)
 ! (last index discerns between hole and electron pockets)
 ABI_ICALLOC(kshe_mask, (fine_ebands%nkpt, ebands%nsppol, 2))

 do spin=1,ebands%nsppol
   ! Get CBM and VBM with some tolerance.
   if (assume_gap) then
     vmax = fine_gaps%vb_max(spin) + tol2 * eV_Ha
     cmin = fine_gaps%cb_min(spin) - tol2 * eV_Ha
   else
     ! Note that we use the Fermi level from ebands instead of fine_ebands.
     vmax = ebands%fermie
     cmin = ebands%fermie
   end if

   do ikf_ibz=1,fine_ebands%nkpt
     do band=1,ebands%mband
       ee = fine_ebands%eig(band, ikf_ibz, spin)
       ! Check whether the interpolated eigenvalue is inside the sigma_erange window.
       if (abs(dtset%sigma_erange(1)) > zero) then
         if (ee <= vmax .and. vmax - ee <= abs(dtset%sigma_erange(1))) then
           kshe_mask(ikf_ibz, spin, 1) = kshe_mask(ikf_ibz, spin, 1) + 1; exit
         end if
       end if
       if (abs(dtset%sigma_erange(2)) > zero) then
         if (ee >= cmin .and. ee - cmin <= abs(dtset%sigma_erange(2))) then
           kshe_mask(ikf_ibz, spin, 2) = kshe_mask(ikf_ibz, spin, 2) + 1; exit
         end if
       end if
     end do
   end do
 end do

 ! Build list of k-points inside pockets. Use over dimensioned array.
 cnt = count(kshe_mask /= 0)
 ABI_MALLOC(krange2ibz, (cnt))
 cnt = 0
 do ikf_ibz=1,fine_ebands%nkpt
   if (any(kshe_mask(ikf_ibz,:,:) /= 0)) then
     cnt = cnt + 1
     krange2ibz(cnt) = ikf_ibz
   end if
 end do
 nkpt_inerange = cnt

 ! Possible extensions that may be implemented at this level:
 !     1. Find image points in the BZ?
 !     2. Compute tetra and q-points for EPH calculation or use +/- wmax window and heuristic approach in sigmaph at runtime?
 !     3. Compute SKW 1st and 2nd derivatives needed to treat Frohlich?

 ! Write output files with k-point list.
 if (my_rank == master .and. len_trim(prefix) /= 0) then
   write(std_out, "(a,i0,a,f5.1,a)")" Found: ",  nkpt_inerange, " kpoints in sigma_erange energy windows. (nkeff / nkibz): ", &
       (100.0_dp * nkpt_inerange) / fine_ebands%nkpt, " [%]"

   ! Write text file with Abinit input variables (mainly for testing purposes).
   path = strcat(prefix, "_KERANGE")
   if (open_file(path, msg, newunit=unt, form="formatted") /= 0) then
     ABI_ERROR(msg)
   end if
   write(unt, "(a)")"kptopt 0"
   write(unt, "(a, i0)")"nkpt ", nkpt_inerange
   write(unt, "(a)")"kpt"
   do ii=1,nkpt_inerange
     write(unt, "(3(es16.8,1x))") fine_ebands%kptns(:, krange2ibz(ii))
   end do
   write(unt, "(a, i0)")"wtk"
   do ii=1,nkpt_inerange
     write(unt, "(es16.8)") fine_ebands%wtk(krange2ibz(ii))
   end do
   close(unt)

   ! Write netcdf file used to perform NSCF run and EPH calculations with eph_task = -4.
   path = strcat(prefix, "_KERANGE.nc")
   NCF_CHECK(nctk_open_create(ncid, path, xmpi_comm_self))
   ! Write crystalline structure, fine_hdr and fine_ebands defined on the fine k-mesh.
   ! fine_ebands will be used to compare with the ab-initio NSCF eigenvalues.
   !
   ! TODO: The size of the KERANGE.nc quickly increases with the k-mesh.
   ! It is ~700 Mb for a ~ 300^3 grid due to occ and eigens
   ! But these quantities are now used in inkpts so it may be possible to avoid writing them to disk.
   !
   NCF_CHECK(fine_hdr%ncwrite(ncid, fform_from_ext("KERANGE.nc"), nc_define=.True.))
   NCF_CHECK(cryst%ncwrite(ncid))
   NCF_CHECK(fine_ebands%ncwrite(ncid))
   NCF_CHECK(nctk_def_dims(ncid, [nctkdim_t("nkpt_inerange", nkpt_inerange)], defmode=.True.))
   ! Define extra arrays.
   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("kshe_mask", "int", "number_of_kpoints, number_of_spins, two"), &
     nctkarr_t("krange2ibz", "int", "nkpt_inerange"), &
     nctkarr_t("sigma_erange", "dp", "two"), &
     nctkarr_t("einterp", "dp", "four") &
   ], defmode=.True.)
   NCF_CHECK(ncerr)
   ! Write extra arrays.
   NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kshe_mask"), kshe_mask))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "krange2ibz"), krange2ibz(1:nkpt_inerange)))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "sigma_erange"), dtset%sigma_erange))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "einterp"), params))
   NCF_CHECK(nf90_close(ncid))
 end if

 ABI_FREE(kshe_mask)
 ABI_FREE(krange2ibz)

 call fine_gaps%free()
 call fine_ebands%free()
 call fine_hdr%free()

end subroutine sigtk_kpts_in_erange
!!***

subroutine degtab_array_free(degtab)
 type(degtab_t),intent(inout) :: degtab(:,:)

 integer :: jj, ii, ideg

 do jj=1,size(degtab, dim=2)
   do ii=1,size(degtab, dim=1)
     if (.not. allocated(degtab(ii, jj)%bids)) cycle
     do ideg=1,size(degtab(ii, jj)%bids)
       ABI_SFREE(degtab(ii, jj)%bids(ideg)%vals)
     end do
     ABI_SFREE(degtab(ii, jj)%bids)
   end do
 end do

end subroutine degtab_array_free
!!***

!----------------------------------------------------------------------

!!****f* m_sigtk/sigtk_sigma_tables
!! NAME
!!  sigtk_sigma_tables
!!
!! FUNCTION
!!  Build tables with the band indices used to compute the matrix elements of sigma_x and sigma_c
!!  taking into account the kind of self-energies and symmetries from esymm.
!!
!! INPUTS
!!  nkcalc: Number of k-points to compute
!!  nkibz: Number of k-points in the IBZ
!!  nsppol: Number of spins
!!  bstart_ks, bstop_ks: First and last band for each (ikcalc, spin)
!!  kcalc2ibz: Mapping kcalc --> IBZ
!!  only_diago: True if only diagonal matrix elements are wanted
!!  sigc_is_herm: True is Sigma_c is Hermitian
!!  [esymm]: Band symmetries
!!
!! OUTPUT
!!  sigxij_tab, sigcij_tab
!!
!! SOURCE

subroutine sigtk_sigma_tables(nkcalc, nkibz, nsppol, bstart_ks, bstop_ks, kcalc2ibz, &
                              only_diago, sigc_is_herm, sigxij_tab, sigcij_tab, esymm)

 use m_gwdefs,        only : sigijtab_t, sigijtab_free
 use m_esymm,         only : esymm_t, esymm_failed

!Arguments ------------------------------------
 integer,intent(in) :: nkcalc, nkibz, nsppol
 integer,intent(in) :: bstart_ks(nkcalc, nsppol), bstop_ks(nkcalc, nsppol)
 logical,intent(in) :: only_diago, sigc_is_herm
 integer,intent(in) :: kcalc2ibz(nkcalc)
 type(sigijtab_t),allocatable,intent(inout) :: Sigxij_tab(:,:), Sigcij_tab(:,:)
 type(esymm_t),optional,intent(in) :: esymm(nkibz, nsppol)

!Local variables-------------------------------
!scalars
 integer :: spin,ikcalc,ik_ibz,bmin,bmax,bcol,brow
 integer :: ii,idx_x,idx_c,irr_idx1,irr_idx2
!arrays
 integer,allocatable :: sigc_bidx(:), sigx_bidx(:)
 logical :: use_sym_at(nkibz, nsppol)

! *************************************************************************

 if (allocated(Sigxij_tab)) then
   call sigijtab_free(Sigxij_tab)
   ABI_FREE(Sigxij_tab)
 end if
 if (allocated(Sigcij_tab)) then
   call sigijtab_free(Sigcij_tab)
   ABI_FREE(Sigcij_tab)
 end if

 ABI_MALLOC(Sigcij_tab, (nkcalc, nsppol))
 ABI_MALLOC(Sigxij_tab, (nkcalc, nsppol))

 use_sym_at = .FALSE.
 if (present(esymm)) then
   ! Create the Sig_ij tables taking advantage of the classification of the bands.
   do spin=1,nsppol
     do ikcalc=1,nkcalc
      ik_ibz = kcalc2ibz(ikcalc)
      use_sym_at(ik_ibz, spin) = .not. esymm_failed(esymm(ik_ibz, spin))
     end do
   end do
 end if

 do spin=1,nsppol
   do ikcalc=1,nkcalc
     ik_ibz = kcalc2ibz(ikcalc)

     if (use_sym_at(ik_ibz, spin)) then
       if (only_diago) then
         ABI_ERROR("You should not be here!")
       end if

       bmin = bstart_ks(ikcalc, spin); bmax = bstop_ks(ikcalc, spin)
       ABI_MALLOC(Sigxij_tab(ikcalc, spin)%col, (bmin:bmax))
       ABI_MALLOC(Sigcij_tab(ikcalc, spin)%col, (bmin:bmax))

       do bcol=bmin,bmax
         ABI_MALLOC(sigc_bidx, (bmax - bmin + 1))
         ABI_MALLOC(sigx_bidx, (bmax - bmin + 1))

         if (esymm(ik_ibz,spin)%err_status /= 0) then
           ! Band classification failed.
           sigc_bidx = [(ii, ii=bmin, bmax)]
           idx_c = bmax - bmin + 1
           sigx_bidx = [(ii,ii=bmin,bcol)] ! Hermitian
           idx_x = bcol - bmin + 1
         else
           irr_idx2 = esymm(ik_ibz,spin)%b2irrep(bcol)
           idx_c = 0
           do brow=bmin,bmax
             irr_idx1 = esymm(ik_ibz,spin)%b2irrep(brow)
             if (sigc_is_herm .and. bcol < brow) CYCLE  ! Only the upper triangle for HF, SEX, or COHSEX.
             if (irr_idx1 == irr_idx2) then ! same character, add this row to the list.
               idx_c = idx_c + 1
               sigc_bidx(idx_c) = brow
             end if
           end do
           idx_x = 0
           do brow=bmin,bcol
             irr_idx1 = esymm(ik_ibz,spin)%b2irrep(brow)
             if (bcol<brow) CYCLE  ! Sig_x is always Hermitian.
             if (irr_idx1 == irr_idx2) then ! same character, add this row to the list.
               idx_x = idx_x +1
               sigx_bidx(idx_x) = brow
             end if
           end do
         end if

         ! Table for Sigma_x matrix elements taking into account symmetries of the bands.
         ABI_MALLOC(Sigxij_tab(ikcalc, spin)%col(bcol)%bidx, (idx_x))

         Sigxij_tab(ikcalc, spin)%col(bcol)%size1 = idx_x
         Sigxij_tab(ikcalc, spin)%col(bcol)%bidx(:) = sigx_bidx(1:idx_x)
         !write(std_out,*)" Sigxij_tab: ikcalc, spin, bcol ",ikcalc,spin,bcol
         !write(std_out,*)" size: ",idx_x,(Sigxij_tab(ikcalc,spin)%col(bcol)%bidx(ii),ii=1,idx_x)
         !
         ! Table for Sigma_c matrix elements taking into account symmetries of the bands.
         ABI_MALLOC(Sigcij_tab(ikcalc, spin)%col(bcol)%bidx, (idx_c))

         Sigcij_tab(ikcalc, spin)%col(bcol)%size1= idx_c
         Sigcij_tab(ikcalc, spin)%col(bcol)%bidx(:) = sigc_bidx(1:idx_c)
         !write(std_out,*)" Sigcij_tab: ikcalc, spin, bcol ",ikcalc,spin,bcol
         !write(std_out,*)" size: ",idx_c,(Sigcij_tab(ikcalc,spin)%col(bcol)%bidx(ii), ii=1,idx_c)

         ABI_FREE(sigx_bidx)
         ABI_FREE(sigc_bidx)
       end do ! bcol

     else
       ! Symmetries cannot be used for this (k,s).
       bmin = bstart_ks(ikcalc, spin); bmax = bstop_ks(ikcalc, spin)
       ABI_MALLOC(Sigcij_tab (ikcalc, spin)%col, (bmin:bmax))
       ABI_MALLOC(Sigxij_tab (ikcalc, spin)%col, (bmin:bmax))

       if (only_diago) then
         ! QP wavefunctions == KS, therefore only diagonal elements are calculated.
         do bcol=bmin,bmax
           ABI_MALLOC(Sigcij_tab(ikcalc, spin)%col(bcol)%bidx, (1:1))
           Sigcij_tab(ikcalc, spin)%col(bcol)%size1= 1
           Sigcij_tab(ikcalc, spin)%col(bcol)%bidx(1) = bcol

           ABI_MALLOC(Sigxij_tab(ikcalc, spin)%col(bcol)%bidx, (1:1))
           Sigxij_tab(ikcalc, spin)%col(bcol)%size1 = 1
           Sigxij_tab(ikcalc, spin)%col(bcol)%bidx(1) = bcol
         end do
       else
         ! Use QP wavefunctions, Sigma_ij matrix is sparse but we have to classify the states in sigma.
         ! The only thing we can do here is filling the entire matrix taking advantage of Hermiticity (if any).
         do bcol=bmin,bmax
           ABI_MALLOC(Sigxij_tab(ikcalc, spin)%col(bcol)%bidx, (bcol-bmin+1))
           Sigxij_tab(ikcalc, spin)%col(bcol)%size1= bcol-bmin+1
           Sigxij_tab(ikcalc, spin)%col(bcol)%bidx(:) = [(ii, ii=bmin,bcol)] ! Sigma_x is Hermitian.
           !write(std_out,*)"Sigxij_tab: ikcalc, spin, bcol ",ikcalc,spin,bcol,Sigxij_tab(ikcalc,spin)%col(bcol)%bidx(:)

           ABI_MALLOC(sigc_bidx, (bmax-bmin+1))
           idx_c = 0
           do brow=bmin,bmax
             if (sigc_is_herm .and. bcol < brow) CYCLE  ! Only the upper triangle of Sigc_ij is needed (SEX, COHSEX).
             idx_c = idx_c +1
             sigc_bidx(idx_c) = brow
           end do
           ABI_MALLOC(Sigcij_tab(ikcalc, spin)%col(bcol)%bidx,(idx_c))
           Sigcij_tab(ikcalc, spin)%col(bcol)%size1= idx_c
           Sigcij_tab(ikcalc, spin)%col(bcol)%bidx(:) = sigc_bidx(1:idx_c)
           ABI_FREE(sigc_bidx)
           !write(std_out,*)"Sigcij_tab: ikcalc, spin, bcol ",ikcalc,spin,bcol,Sigcij_tab(ikcalc,spin)%col(bcol)%bidx(:)
         end do
       end if
     end if

   end do !ikcalc
 end do !spin

end subroutine sigtk_sigma_tables
!!***

end module m_sigtk
!!***
