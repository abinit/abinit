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

 use m_fstrings,     only : sjoin, ltoa
 use defs_datatypes, only : ebands_t
 use defs_abitypes,  only : dataset_type
 use m_kpts,         only : kpts_ibz_from_kptrlatt, kpts_timrev_from_kptopt, listkk

 implicit none

 private

 public :: sigtk_kcalc_from_nkptgw
 public :: sigtk_kcalc_from_qprange
 public :: sigtk_kcalc_from_gaps
 public :: sigtk_kcalc_from_erange
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
       write(msg,'(a,2i0,2(a,i0))')&
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
 if (my_rank == 0) then
   write(std_out, "(a)")" Selecting k-points and bands according to their position wrt band edges (sigma_erange)."
   !call gaps%print()
   write(std_out, "(a, 2(f6.3, 1x), a)")" sigma_erange: ", dtset%sigma_erange(:) * Ha_eV, " (eV)"
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
         1, cryst%symafm, cryst%symrec, timrev, comm, use_symrec=.True.)
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
   ! index in ib_work array
   ii = kpos(ic)
   ! index in ebands.
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

 if (my_rank == 0) then
   write(std_out, "(a, i0, a, 2(f6.3, 1x), a)") &
   " Found ", nkcalc, " k-points within erange: ", dtset%sigma_erange(:) * Ha_eV, " (eV)"
   write(std_out, "(2(a, i0))")" min(nbcalc_ks): ", minval(nbcalc_ks), " MAX(nbcalc_ks): ", maxval(nbcalc_ks)
 end if

 ABI_FREE(ib_work)
 ABI_FREE(sigmak2ebands)

end subroutine sigtk_kcalc_from_erange
!!***

end module m_sigtk
!!***
