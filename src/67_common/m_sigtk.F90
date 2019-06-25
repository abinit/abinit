!{\src2tex{textfont=tt}}
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
!!  Copyright (C) 2008-2018 ABINIT group (MG)
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

module m_sigtk

 use defs_basis
 use m_abicore
 use m_errors
 use m_ebands
 use m_crystal
 use m_xmpi
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_nctk
 use m_hdr

 use m_fstrings,     only : sjoin, ltoa, strcat
 use m_io_tools,     only : open_file
 use defs_datatypes, only : ebands_t, pseudopotential_type
 use defs_abitypes,  only : dataset_type, hdr_type
 use defs_wvltypes,  only : wvl_internal_type
 use m_pawtab,       only : pawtab_type
 use m_kpts,         only : kpts_ibz_from_kptrlatt, kpts_timrev_from_kptopt, listkk

 implicit none

 private

 public :: sigtk_kcalc_from_nkptgw
 public :: sigtk_kcalc_from_qprange
 public :: sigtk_kcalc_from_gaps
 public :: sigtk_kcalc_from_erange
 public :: sigtk_kpts_in_erange
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
!! PARENTS
!!
!! CHILDREN
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
       MSG_WARNING(msg)
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
!! PARENTS
!!
!! CHILDREN
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
 integer,allocatable,intent(out) :: bstart_ks(:,:)
 integer,allocatable,intent(out) :: nbcalc_ks(:,:)

!Local variables ------------------------------
!scalars
 integer :: spin, ik, bstop, mband, sigma_nkbz
!arrays
 integer :: kptrlatt(3,3)
 integer :: val_indeces(ebands%nkpt, ebands%nsppol)
 real(dp),allocatable :: sigma_wtk(:),sigma_kbz(:,:)

! *************************************************************************

 mband = ebands%mband

 val_indeces = get_valence_idx(ebands)

 if (any(dtset%sigma_ngkpt /= 0)) then
    call wrtout(std_out, " Generating list of k-points for self-energy from sigma_nkpt and qprange.")
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
       bstart_ks(ik,spin) = max(val_indeces(ik,spin) - qprange, 1)
       bstop = min(val_indeces(ik,spin) + qprange, mband)
       nbcalc_ks(ik,spin) = bstop - bstart_ks(ik,spin) + 1
     end do
   end do

 else
   call wrtout(std_out, " Including all occupied states and -qprange empty states.")
   bstart_ks = 1
   do spin=1,dtset%nsppol
     do ik=1,nkcalc
       nbcalc_ks(ik,spin) = min(val_indeces(ik,spin) - qprange, mband)
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
!! PARENTS
!!
!! CHILDREN
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
 integer :: val_indeces(ebands%nkpt, ebands%nsppol)
 integer :: kpos(6)

! *************************************************************************

 ABI_UNUSED((/dtset%natom/))

 call wrtout(std_out, " Including direct and fundamental KS gap in Sigma_nk")
 ABI_CHECK(maxval(gaps%ierr) == 0, "qprange 0 cannot be used because I cannot find the gap (gap_err !=0)")

 nsppol = ebands%nsppol
 val_indeces = get_valence_idx(ebands)

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
     bstart_ks(ii,spin) = val_indeces(ik,spin)
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
!!  wrt to the band edges and the value of sigma_erange. Useful when computing electron-lifetimes
!!  for transport calculations.
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
!! PARENTS
!!
!! CHILDREN
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
 integer :: spin, ik, band, ii, ic, nsppol, tmp_nkpt, timrev, sigma_nkbz, my_rank
 logical :: found
 real(dp) :: cmin, vmax, ee, dksqmax
 character(len=500) :: msg
!arrays
 integer :: kptrlatt(3,3)
 integer,allocatable :: ib_work(:,:,:), sigmak2ebands(:), indkk(:,:)
 integer :: kpos(ebands%nkpt)
 real(dp),allocatable :: sigma_wtk(:),sigma_kbz(:,:),tmp_kcalc(:,:)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm) !; nprocs = xmpi_comm_size(comm)
 if (my_rank == master) then
   write(std_out, "(a)")" Selecting k-points and bands according to their position wrt band edges (sigma_erange)."
   write(std_out, "(a, 2(f6.3, 1x), a)")" sigma_erange: ", dtset%sigma_erange(:) * Ha_eV, " (eV)"
   call gaps%print(unit=std_out)
 end if

 ABI_CHECK(maxval(gaps%ierr) == 0, "erange 0 cannot be used because I cannot find the gap (gap_err !=0)")

 if (any(dtset%sigma_ngkpt /= 0)) then
    call wrtout(std_out, sjoin(" Generating initial list of k-points from sigma_nkpt.", ltoa(dtset%sigma_ngkpt)))
    ! Get tentative tmp_nkpt and tmp_kcalc from sigma_ngkpt.
    kptrlatt = 0
    kptrlatt(1,1) = dtset%sigma_ngkpt(1); kptrlatt(2,2) = dtset%sigma_ngkpt(2); kptrlatt(3,3) = dtset%sigma_ngkpt(3)
    call kpts_ibz_from_kptrlatt(cryst, kptrlatt, dtset%kptopt, dtset%sigma_nshiftk, dtset%sigma_shiftk, &
      tmp_nkpt, tmp_kcalc, sigma_wtk, sigma_nkbz, sigma_kbz)

    ABI_FREE(sigma_kbz)
    ABI_FREE(sigma_wtk)

    ! Map tmp_kcalc to ebands%kpts
    timrev = kpts_timrev_from_kptopt(ebands%kptopt)
    ABI_MALLOC(indkk, (tmp_nkpt,  6))
    call listkk(dksqmax, cryst%gmet, indkk, ebands%kptns, tmp_kcalc, ebands%nkpt, tmp_nkpt, cryst%nsym, &
         1, cryst%symafm, cryst%symrec, timrev, comm, exit_loop=.True., use_symrec=.True.)
    if (dksqmax > tol12) then
      write(msg, '(a,es16.6,2a)' )&
        "At least one of the k-points could not be generated from a symmetrical one in the WFK. dksqmax: ",dksqmax, ch10,&
        'Action: check your WFK file and the value of sigma_nkpt, sigma_shiftk in the input file.'
      MSG_ERROR(msg)
    end if

    ABI_MALLOC(sigmak2ebands, (tmp_nkpt))
    sigmak2ebands = indkk(:, 1)
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
   ! Get cmb and vbm with some tolerance
   vmax = gaps%vb_max(spin) + tol2 * eV_Ha
   cmin = gaps%cb_min(spin) - tol2 * eV_Ha
   do ii=1,tmp_nkpt
     ! Index of k-point in ebands.
     ik = sigmak2ebands(ii)
     ib_work(1, ii, spin) = huge(1)
     ib_work(2, ii, spin) = -huge(1)
     do band=1,ebands%nband(ik + (spin-1) * ebands%nkpt)
        ee = ebands%eig(band, ik, spin)
        if (dtset%sigma_erange(1) > zero) then
          if (ee <= vmax .and. vmax - ee <= dtset%sigma_erange(1)) then
            ib_work(1, ii, spin) = min(ib_work(1, ii, spin), band)
            ib_work(2, ii, spin) = max(ib_work(2, ii, spin), band)
            !write(std_out, *), "Adding valence band", band, " with ee [eV]: ", ee * Ha_eV
          end if
        end if
        if (dtset%sigma_erange(2) > zero) then
          if (ee >= cmin .and. ee - cmin <= dtset%sigma_erange(2)) then
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
 ! we have compute the union of the k-points.
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
   ! Index in ib_work array
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
       MSG_WARNING("Spin-polarized case with nbcalc_ks == 0, don't know if code can handle it!")
     end if
   end do
 end do

 if (my_rank == master) then
   write(std_out, "(a, i0, a, 2(f6.3, 1x), a)") &
   " Found ", nkcalc, " k-points within erange: ", dtset%sigma_erange(:) * Ha_eV, " (eV)"
   write(std_out, "(2(a, i0))")" min(nbcalc_ks): ", minval(nbcalc_ks), " MAX(nbcalc_ks): ", maxval(nbcalc_ks)
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
!!  Use star functions and [[einterp]] to interpolate electron energies onto fine dense defined 
!!  by [[sigma_ngkpt]] and [[sigma_shiftk]].
!!  Find k-points inside (electron/hole) pockets according to the values specifed in [[sigma_erange]].
!!  Write KERANGE.nc file with the tables required by the code to automate NSCF band structure calculations
!!  and electron lifetime computation in the EPH code.
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
!! PARENTS
!!
!! CHILDREN
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
 real(dp) :: ee, cmin, vmax
 character(len=500) :: msg
 character(len=fnlen) :: path
 type(ebands_t) :: fine_ebands
 type(gaps_t) :: gaps, fine_gaps
 type(wvl_internal_type) :: dummy_wvl
 type(hdr_type) :: fine_hdr
!arrays
 integer :: fine_kptrlatt(3,3), band_block(2)
 integer,allocatable :: kshe_mask(:,:,:), krange2ibz(:)
 real(dp) :: params(4)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 if (my_rank == master) then
   write(std_out, "(2a)")ch10, repeat("=", 92)
   write(std_out, "(a)") " Finding k-points inside (electron/hole) pockets."
   write(std_out, "(a, 2(f6.3, 1x), a)")" sigma_erange: ", dtset%sigma_erange(:) * Ha_eV, " (eV)"
   write(std_out, "(2a)") "  Interpolating eigenvalues onto dense K-mesh using sigma_ngkpt: ", trim(ltoa(dtset%sigma_ngkpt))
   write(std_out, "(2a)") "  and sigma_shiftk shifts:"
   ABI_CHECK(allocated(dtset%sigma_shiftk), "sigma_shiftk and sigma_nshiftk must be specified in input.")
   do ii=1,dtset%nshiftk
     write(std_out, "(a, 3(f2.1, 1x))")"   sigma_shiftk:", dtset%sigma_shiftk(:, ii)
   end do
   write(std_out, "(2a)")repeat("=", 92),ch10
   !call ebands_print(ebands, header, unit=std_out, prtvol=dtset%prtvol)

   ! Consistency check.
   if (.not. any(dtset%sigma_erange > zero)) then
     MSG_ERROR("sigma_erange must be specified in input when calling sigtk_kpts_in_erange.")
   end if
   if (all(dtset%sigma_ngkpt == 0)) then
     MSG_ERROR("sigma_ngkpt must be specified in input when calling sigtk_kpts_in_erange.")
   end if
 end if

 ! Compute gaps using input ebands.
 gap_err = get_gaps(ebands, gaps)
 if (gap_err /= 0) then
   MSG_ERROR("Cannot compute fundamental and direct gap (likely metal).")
 end if
 call gaps%print(header="Gaps from input ebands", unit=std_out)
 call gaps%free()

 ! Interpolate band energies with star-functions.
 ! In the EPH code, we will need eigens in the IBZ to compute efermi not just energies inside pockets.
 fine_kptrlatt = 0
 do ii=1,3
   fine_kptrlatt(ii, ii) = dtset%sigma_ngkpt(ii)
 end do
 band_block = [1, ebands%mband]
 params = 0; params(1) = 1; params(2) = 5; if (nint(dtset%einterp(1)) == 1) params = dtset%einterp

 fine_ebands = ebands_interp_kmesh(ebands, cryst, params, fine_kptrlatt, &
                                   dtset%sigma_nshiftk, dtset%sigma_shiftk, band_block, comm)

 call ebands_update_occ(fine_ebands, dtset%spinmagntarget, prtvol=dtset%prtvol)
 call ebands_print(fine_ebands, header="FINE EBANDS", unit=std_out, prtvol=dtset%prtvol)

 ! Interpolate bands on k-path.
 !if (nint(dtset%einterp(1)) == 0)
 !call ebands_interpolate_kpath(ebands, dtset, cryst, band_block, prefix, comm)
 !end do

 ! Compute gaps using fine_ebands.
 gap_err = get_gaps(fine_ebands, fine_gaps)
 if (gap_err /= 0) then
   MSG_ERROR("Cannot compute fundamental and direct gap (likely metal).")
 end if
 call fine_gaps%print(header="Gaps from interpolated eigenvalues", unit=std_out)

 ! Build new header with fine k-mesh (note kptrlatt_orig == kptrlatt)
 call hdr_init_lowlvl(fine_hdr, fine_ebands, psps, pawtab, dummy_wvl, ABINIT_VERSION, pertcase0, &
   dtset%natom, dtset%nsym, dtset%nspden, dtset%ecut, dtset%pawecutdg, dtset%ecutsm, dtset%dilatmx, &
   dtset%intxc, dtset%ixc, dtset%stmbias, dtset%usewvl, dtset%pawcpxocc, dtset%pawspnorb, dtset%ngfft, dtset%ngfftdg, &
   dtset%so_psp, dtset%qptn, cryst%rprimd, cryst%xred, cryst%symrel, cryst%tnons, cryst%symafm, cryst%typat, &
   dtset%amu_orig(:, image1), dtset%icoulomb, &
   dtset%kptopt, dtset%nelect, dtset%charge, fine_kptrlatt, fine_kptrlatt, &
   dtset%sigma_nshiftk, dtset%sigma_nshiftk, dtset%sigma_shiftk, dtset%sigma_shiftk)

 ! Find k-points inside sigma_erange energy window.
 ! Set entry to 1 if (ikpt, spin) is inside the pocket (last index discerns between hole and electron pockets)
 ABI_ICALLOC(kshe_mask, (fine_ebands%nkpt, ebands%nsppol, 2))

 do spin=1,ebands%nsppol
   ! Get CBM and VBM with some tolerance.
   vmax = fine_gaps%vb_max(spin) + tol2 * eV_Ha
   cmin = fine_gaps%cb_min(spin) - tol2 * eV_Ha
   do ikf_ibz=1,fine_ebands%nkpt
     do band=1,ebands%mband
       ee = fine_ebands%eig(band, ikf_ibz, spin)
       ! Check whether the interpolated eigenvalue is inside the sigma_erange window.
       if (dtset%sigma_erange(1) > zero) then
         if (ee <= vmax .and. vmax - ee <= dtset%sigma_erange(1)) then
           kshe_mask(ikf_ibz, spin, 1) = kshe_mask(ikf_ibz, spin, 1)  + 1
           exit
         end if
       end if
       if (dtset%sigma_erange(2) > zero) then
         if (ee >= cmin .and. ee - cmin <= dtset%sigma_erange(2)) then
           kshe_mask(ikf_ibz, spin, 2) = kshe_mask(ikf_ibz, spin, 2)  + 1
           exit
         end if
       end if
     end do
   end do
 end do

 ! Build list of k-points inside pockets.
 nkpt_inerange = count(kshe_mask /= 0)
 ABI_MALLOC(krange2ibz, (nkpt_inerange))
 cnt = 0
 do ikf_ibz=1,fine_ebands%nkpt
   if (any(kshe_mask(ikf_ibz,:,:) /= 0)) then
     cnt = cnt + 1
     krange2ibz(cnt) = ikf_ibz
   end if
 end do

 ! Possible extensions that may be implemented at this level:
 !     Find image points in the BZ?
 !     Compute tetra and q-points for EPH calculation or use +/- wmax window and heuristic approach in sigmaph at runtime?
 !     Compute SKW 1st and 2nd derivatives needed to treat Frohlich?

 ! Write output files with k-point list.
 if (my_rank == master .and. len_trim(prefix) /= 0) then
   write(std_out, "(a,i0,a,f5.1,a)")" Found: ",  nkpt_inerange, " kpoints in sigma_erange energy windows. (nkeff / nkibz): ", &
       (100.0_dp * nkpt_inerange) / fine_ebands%nkpt, " [%]"
   ! Write text file with Abinit input variables (mainly for testing purposes).
   path = strcat(prefix, "_KERANGE")
   if (open_file(path, msg, newunit=unt, form="formatted") /= 0) then
     MSG_ERROR(msg)
   end if
   write(unt, "(a)")"kptopt 0"
   write(unt, "(a, i0)")"nkpt ", nkpt_inerange
   write(unt, "(a)")"kpt"
   do ii=1,nkpt_inerange
     write(unt, "(3(es16.8,1x))")fine_ebands%kptns(:, krange2ibz(ii))
   end do
   write(unt, "(a, i0)")"wtk"
   do ii=1,nkpt_inerange
     write(unt, "(es16.8)")fine_ebands%wtk(krange2ibz(ii))
   end do
   close(unt)

   ! Write netcdf file used to perform NSCF run and EPH calculations with eph_task = -4.
   path = strcat(prefix, "_KERANGE.nc")
#ifdef HAVE_NETCDF
   NCF_CHECK(nctk_open_create(ncid, path, xmpi_comm_self))
   ! Write crystalline structure, fine_hdr and fine_ebands defined on the fine k-mesh.
   ! fine_ebands will be used to compare with the ab-initio NSCF eigenvalues.
   NCF_CHECK(hdr_ncwrite(fine_hdr, ncid, fform_from_ext("KERANGE.nc"), nc_define=.True.))
   NCF_CHECK(cryst%ncwrite(ncid))
   NCF_CHECK(ebands_ncwrite(fine_ebands, ncid))
   ncerr = nctk_def_dims(ncid, [nctkdim_t("nkpt_inerange", nkpt_inerange)], defmode=.True.)
   NCF_CHECK(ncerr)
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
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "krange2ibz"), krange2ibz))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "sigma_erange"), dtset%sigma_erange))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "einterp"), params))
   NCF_CHECK(nf90_close(ncid))
#endif
 end if

 ABI_FREE(kshe_mask)
 ABI_FREE(krange2ibz)

 call fine_gaps%free()
 call ebands_free(fine_ebands)
 call hdr_free(fine_hdr)

end subroutine sigtk_kpts_in_erange
!!***

end module m_sigtk
!!***
