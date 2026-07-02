!!****m* ABINIT/m_gstore_converters
!! NAME
!! m_gstore_converters
!!
!! FUNCTION
!!  Convert data from gstore.nc to other formats.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2026 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_gstore_converters

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_clib

 use m_io_tools,       only : open_file
 use m_dtset,          only : dataset_type
 use m_dtfil,          only : datafiles_type
 use m_fstrings,       only : sjoin, itoa, strcat
 use m_crystal,        only : crystal_t
 use m_ebands,         only : ebands_t, gaps_t
 use m_ifc,            only : ifc_type
 use m_gstore,         only : gstore_t, GSTORE_GMODE_PHONON, gstore_read_gtype

 implicit none

 private

 public :: gstore_convert
 ! Convert data from gstore.nc to other formats

!!***

contains
!!***

!----------------------------------------------------------------------

!!****f* m_gstore_converters/gstore_convert
!! NAME
!! gstore_convert
!!
!! FUNCTION
!! Convert data from gstore.nc to other formats
!!
!! INPUTS
!! gstore_path=Filename of the output GSTORE.nc file
!!
!! SOURCE

subroutine gstore_convert(gstore_path, dtset, dtfil, cryst, ebands, ifc, comm)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: gstore_path
 type(dataset_type),target,intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil
 class(crystal_t),target,intent(in) :: cryst
 class(ebands_t),target,intent(in) :: ebands
 class(ifc_type),target,intent(in) :: ifc
 integer,intent(in) :: comm

!Local variables-------------------------------
!scalars
 integer :: nprocs, my_rank, nsppol, spin, nmodes, this_comm, unt, ib, nu, i, j, ierr
 integer :: with_cplex, ik_ibz, my_is, my_ik, my_iq, iq_glob, natom, itypat, lstr_j
 integer :: ibrav, idir, iat, ipert, unt_ascii, ik_glob, band_kq, band_k, mu
 logical :: with_g2dw, q_is_gamma, lborn, ascii_write
 real(dp),parameter :: Ha2Ry = two
 real(dp) :: weight_qq
 character(len=5000) :: msg
 character(len=abi_slen) :: with_gmode, gvals_name, gtype
 character(len=fnlen) :: fname, fname_ascii, elphmat_dir, prefix
 character(len=3) :: band_i
 type(gstore_t) :: gstore
!arrays
 integer :: units(2)
 real(dp) :: qpt(3), kk_bz(3), kk_ibz(3), celldm(6)
 real(dp) :: wfact
 real(dp),allocatable :: tau_cart(:,:)
 complex(dp),allocatable :: dyn_qe(:,:), g_cart(:,:,:)
 character(len=3) :: atm(cryst%ntypat)
!----------------------------------------------------------------------

 ! Only master works here as performance is not crucial.
 ! In principle one could activate the q-point/spin parallelism just to distribute
 ! the memory for the g's to avoid going OOM.
 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 if (my_rank /= 0) return

 units = [std_out, ab_out]
 natom = cryst%natom; nmodes = 3 * cryst%natom; nsppol = ebands%nsppol; this_comm = xmpi_comm_self

 ! QE Bravais lattice metadata.
 ! We do not classify the QE Bravais lattice from the ABINIT cell, hence ibrav = 0
 ! (free lattice) and celldm(1) = alat = |a1| in Bohr (the remaining celldm are unused
 ! for ibrav == 0). The full lattice is otherwise defined by the crystal structure.
 ibrav = 0
 celldm(:) = zero
 celldm(1) = sqrt(sum(cryst%rprimd(:, 1) ** 2))

 ! Work array for the phonon displacement matrix (per q-point).
 ABI_MALLOC(dyn_qe, (nmodes, nmodes))

 ! QE tau: ionic positions in Cartesian coordinates, in units of alat (= celldm(1)).
 ABI_MALLOC(tau_cart, (3, natom))
 do iat=1,natom
   tau_cart(:, iat) = matmul(cryst%rprimd, cryst%xred(:, iat)) / celldm(1)
 end do

 ! If ascii_write is .True. a human-readable copy of each binary elph.mat.q_* file is
 ! also written, with the same content (record by record) and a ".ascii" suffix.
 ! Toggle this flag (or wire it to an input variable) to disable the extra files.
 ascii_write = .True.

 ! Preliminary consistency check.
 call wrtout(units, sjoin(" Begin conversion GSTORE --> ", dtset%gstore_convert))
 ABI_CHECK(gstore_path /= ABI_NOFILE, sjoin("Invalid gstore_path:", gstore_path))
 ABI_CHECK(dtset%gstore_convert == "epiq", "only gstore_convert == 'epiq' is supported.")
 !ABI_CHECK(.not. cryst%isalchemical(), "'epiq' format does not support alchemical pseudos.")

 ! Read g(k,q) from GSTORE and store them in gqk%my_g(nu, im_kq, my_iq, in_k, ik)
 !
 ! Shape of array is:
 !
 ! (my_npert, nb_kq, my_nq, nb_k, my_nk)
 ! (       p, b1_kq,     q, b2_k, k)  -->  <k+q, b1| D_{q,p}H |k, b2>

 ! The g's are complex and in the phonon representation.
 ! Note <k+q| for the final state and |k> for the initial state. I guess epiq uses the same convention.
 ! All quantities are in atomic units (Hartree and Bohr).
 !
 ! The gstore file produced by the GWPT code has both GWPT and KS g.
 ! In this case, we convert the GWPT matrix elements but one can still select
 ! the KS e-ph vertex via gstore_gname.

 call gstore_read_gtype(gstore_path, gtype, this_comm)
 gvals_name = "gvals"
 if (gtype == "gwpt" .and. dtset%gstore_gname == "gvals_ks") gvals_name = "gvals_ks"

 with_cplex = 2; with_gmode = GSTORE_GMODE_PHONON; with_g2dw = .False.

 call gstore%from_ncpath(gstore_path, with_cplex, dtset, dtfil, cryst, ebands, ifc, &
                         with_gmode, gvals_name, with_g2dw, this_comm)

 ! Consistency check.
 ABI_CHECK(nsppol == 1, "Don't know how to convert spin-polarized g to epiq format!")

 ! For wannierization, we need the same number of bands for m and n.
 ! Also, k and q must be in the BZ without any filter.
 ! Once the symmetrization of the g's has been implemented, this routine
 ! will receive a gstore file in which all g(k,q) matrix elements in the BZ
 ! have been reconstructed using symmetry operations.
 ABI_CHECK(gstore%same_nbands(msg), msg)
 if (gstore%check_cplex_qkzone_gmode(2, "bz", "bz", "phonon", kfilter="none") /= 0) then
   ABI_ERROR("GSTORE.nc should have both k and q in the full BZ. See messages above.")
 end if

 ! Create directory to store output files.
 prefix = "epiq"
 elphmat_dir = strcat(dtfil%filnam_ds(4), "_", prefix)
 call wrtout(units, sjoin(" Output files written to directory:", elphmat_dir))
 call execute_command_line(sjoin("rm -rf", elphmat_dir), exitstat=ierr)
 call clib_mkdir_if_needed(elphmat_dir, ierr)
 ABI_CHECK(ierr == 0, "mkdir returned ierr /= 0")

 ! Write the EPIQ input namelist (&Diff_Start_Param + KPOINTS), the analog of QE's
 ! print_ph_input2epiq called inside ep_matrix_element_wannier.
 call write_epiq_input(dtset, ebands, strcat(elphmat_dir, "/scf_dfpt.2epiq.in"))
 call wrtout(units, sjoin(" EPIQ input namelist written to:", strcat(elphmat_dir, "/scf_dfpt.2epiq.in")))

 ! Write the dynq0 file: q-mesh, number of irreducible q-points and their positions
 ! in Cartesian coordinates (2pi/alat), the analog of QE's dynq0 output.
 call write_dynq0(gstore%ngqpt, gstore%qibz, cryst%gprimd, celldm(1), strcat(elphmat_dir, "/dynq0"))
 call wrtout(units, sjoin(" EPIQ dynq0 file written to:", strcat(elphmat_dir, "/dynq0")))

 ! NB: atm is character(len=3) while symbol_type returns character(len=2).
 do itypat=1, cryst%ntypat
   atm(itypat)(1:2) = cryst%symbol_type(itypat)
   atm(itypat)(3:3) = ""
 end do

 ! Write one dynamical-matrix file (dynq<iq>) per irreducible q-point, in QE format
 ! (same ordering as the dynq0 list).
 do iq_glob=1,gstore%nqibz
   call write_dynq(cryst, ifc, gstore%qibz(:,iq_glob), celldm(1), atm, &
                   strcat(elphmat_dir, "/dynq", itoa(iq_glob)))
 end do
 call wrtout(units, sjoin(" EPIQ dynq<iq> dynamical-matrix files written to directory:", elphmat_dir))

 ! TODO: Need helper function to get ibrav, celldm from crystal
 ! Some input from the EPIC developers would be greatly appreciated.
 !call cryst%get_ibrav_celldm(ibrav, celldm)

 ! Loop over collinear spins.
 do my_is=1,gstore%my_nspins
   spin = gstore%my_spins(my_is)
   associate (gqk => gstore%gqk(my_is))

   ! Buffer for the e-ph matrix elements rotated to the Cartesian atomic-displacement basis.
   ABI_MALLOC(g_cart, (nmodes, gqk%nb_kq, gqk%nb_k))

   !num_bands = gqk%nb_k

   ! Loop over q-points in the BZ.
   do my_iq=1, gqk%my_nq
     iq_glob = my_iq + gqk%my_qstart - 1

     call gqk%myqpt(my_iq, gstore, weight_qq, qpt); q_is_gamma = sum(qpt**2) < tol14
     call define_band_string(iq_glob, band_i, lstr_j)

     fname = trim(elphmat_dir)//"/"//trim(prefix)//'_elph.mat.q_'//band_i(1:lstr_j)
     if (open_file(fname, msg, newunit=unt, form="unformatted", status="unknown", action="write") /= 0) then
       ABI_ERROR(msg)
     end if

     ! Optional human-readable companion file with the same content.
     if (ascii_write) then
       fname_ascii = trim(fname)//".ascii"
       if (open_file(fname_ascii, msg, newunit=unt_ascii, form="formatted", status="unknown", action="write") /= 0) then
         ABI_ERROR(msg)
       end if
     end if

     ! EPIQ format. See https://gitlab.com/the-epiq-team/epiq/-/blob/develop/src/io_matelem.F90
     !read(unt) (xq_r(j),j=1,3)
     !if (.not.fet) read(unt) noncolin, nspin, lborn !REMOVE FOR VERSION 5.1 FET
     !read(unt) nel_aux
     !read(unt) nbnd_min, nbnd_max, nbnd_r !nbnd_r = total number of bands in pw
     !read(unt) nmodes, nk_r, nat, ntyp
     !read(unt) ibrav,(celldm(j), j=1,6)
     !read(unt) (atm(j),j=1,ntyp),(amass(j),j=1,ntyp), &
     !           (ityp(j),j=1,nat),((tau(j,i),j=1,3),i=1,nat)
     !read(unt) (w2 (nu,iqph) , nu=1,nmodes)
     !read(unt) ((zz(i,j,iqph), i=1,nmodes),j=1,nmodes)   !eigenvectors in the QE basis
     !read(unt) ((dyn(i,j,iqph),i=1,nmodes),j=1,nmodes)  ! eigenvectors divided by masses
     !do k=1,num_kpts
     !  read(unt) (xk_r(i,k),i=1,3)
     !  read(unt) (eig(i,k),i=1,num_bands)
     !  do nu=1,nmodes
     !    write(unt) ((g_matrix(j, i, nu, k, iq),j=1,num_bands),i=1,num_bands)
     !  end do
     !end do

     ! Conventions used below (matched to QE pw/PHonon elphsum_wannier and the EPIQ reader):
     ! - w2 is omega^2 (squared phonon frequency) in Ry^2, signed (< 0 for imaginary modes).
     ! - lborn = .False.: BECS / dynamical quadrupoles are not exported yet.
     ! - tau are atom positions in Cartesian coordinates, in units of alat (= xcart/alat).
     ! - nel_aux (ebands%nelect) is the (real) number of electrons including possible doping.
     ! - The e-ph matrix elements are written in the Cartesian atomic-displacement basis
     !   (as in QE ep_matrix_element_wannier), obtained by rotating the ABINIT mode-basis
     !   g with the phonon eigenvectors. Consequently zz (the QE pattern matrix) is the
     !   identity, while dyn keeps the Cartesian phonon displacements (= my_displ_cart).
     lborn = .False.

     ! Record 1: q-point in reduced (crystal) coordinates.
     write(unt) qpt
     if (ascii_write) then
       write(unt_ascii, '(a)') "# Record 1: q-point (reduced coordinates)"
       write(unt_ascii, '(3es24.15)') qpt
     end if

     ! Record 2: noncolin, nspin, lborn. BECS/quadrupoles not exported yet => lborn = .False.
     write(unt) ebands%nspinor == 2, nsppol, lborn
     if (ascii_write) then
       write(unt_ascii, '(a)') "# Record 2: noncolin, nspin, lborn"
       write(unt_ascii, *) ebands%nspinor == 2, nsppol, lborn
     end if

     ! Record 3: number of electrons (including possible doping).
     write(unt) ebands%nelect
     if (ascii_write) then
       write(unt_ascii, '(a)') "# Record 3: number of electrons"
       write(unt_ascii, '(es24.15)') ebands%nelect
     end if

     ! Record 4: first band, last band, total number of bands in the pw calculation.
     write(unt) gqk%bstart_k, gqk%bstop_k, ebands%mband
     if (ascii_write) then
       write(unt_ascii, '(a)') "# Record 4: nbnd_min, nbnd_max, nbnd_total"
       write(unt_ascii, '(3i8)') gqk%bstart_k, gqk%bstop_k, ebands%mband
     end if

     ! Record 5: number of modes, number of k-points, number of atoms, number of atom types.
     write(unt) nmodes, gqk%glob_nk, natom, cryst%ntypat
     if (ascii_write) then
       write(unt_ascii, '(a)') "# Record 5: nmodes, nkpt, natom, ntypat"
       write(unt_ascii, '(4i8)') nmodes, gqk%glob_nk, natom, cryst%ntypat
     end if

     ! Record 6: Bravais lattice index and cell dimensions (see ibrav/celldm comment above).
     write(unt) ibrav, (celldm(j), j=1,6)
     if (ascii_write) then
       write(unt_ascii, '(a)') "# Record 6: ibrav, celldm(1:6)"
       write(unt_ascii, '(i6,6es24.15)') ibrav, (celldm(j), j=1,6)
     end if

     ! Record 7: atom symbols, atomic masses (atomic mass unit), atom types and
     !           Cartesian atomic positions tau in units of alat.
     write(unt) (atm(j), j=1,cryst%ntypat), (cryst%amu(j), j=1,cryst%ntypat), &
                (cryst%typat(j), j=1,natom), ((tau_cart(j,i), j=1,3), i=1,natom)
     if (ascii_write) then
       write(unt_ascii, '(a)') "# Record 7: atom symbols, masses (amu), types, tau (Cartesian, alat units)"
       write(unt_ascii, '(*(a3,1x))') (atm(j), j=1,cryst%ntypat)
       write(unt_ascii, '(*(es24.15,1x))') (cryst%amu(j), j=1,cryst%ntypat)
       write(unt_ascii, '(*(i6,1x))') (cryst%typat(j), j=1,natom)
       do i=1,natom
         write(unt_ascii, '(3es24.15)') (tau_cart(j,i), j=1,3)
       end do
     end if

     ! Record 8: squared phonon frequencies omega^2 in Ry^2.
     ! my_wnuq is the signed frequency omega in Ha (negative for imaginary modes), so
     ! omega^2 = sign(omega) * (omega*Ha2Ry)^2 = omega*|omega|*Ha2Ry^2, matching QE's w2.
     write(unt) (gqk%my_wnuq(nu, my_iq) * abs(gqk%my_wnuq(nu, my_iq)) * Ha2Ry ** 2, nu=1,nmodes)
     if (ascii_write) then
       write(unt_ascii, '(a)') "# Record 8: phonon frequencies squared omega^2 (Ry^2)"
       write(unt_ascii, '(*(es24.15,1x))') &
         (gqk%my_wnuq(nu, my_iq) * abs(gqk%my_wnuq(nu, my_iq)) * Ha2Ry ** 2, nu=1,nmodes)
     end if

     ! Phonon displacements for this q-point (composite index ipert = idir + 3*(iat-1),
     ! Cartesian direction fast, atom slow, matching the QE mode/perturbation ordering):
     !   dyn_qe = phonon displacements (eigenvectors divided by sqrt(mass)) = my_displ_cart.
     ! These are written as 'dyn' (record 10) and used below to rotate g to the Cartesian basis.
     do nu=1,nmodes
       do iat=1,natom
         do idir=1,3
           ipert = idir + 3 * (iat - 1)
           dyn_qe(ipert, nu) = cmplx(gqk%my_displ_cart(1, idir, iat, nu, my_iq), &
                                     gqk%my_displ_cart(2, idir, iat, nu, my_iq), kind=dp)
         end do
       end do
     end do

     ! Record 9: eigenvectors in the QE basis (zz).
     ! Like QE's ep_matrix_element_wannier, the e-ph matrix elements are written in the
     ! Cartesian (atomic-displacement) basis (see the rotation in the k-loop below), so the
     ! QE pattern matrix u reduces to the identity. We therefore write the identity matrix.
     write(unt) ((cmplx(merge(one, zero, i == j), zero, kind=dp), i=1,nmodes), j=1,nmodes)
     if (ascii_write) then
       write(unt_ascii, '(a)') "# Record 9: eigenvectors in the QE basis zz = identity (Cartesian g)"
       do j=1,nmodes
         do i=1,nmodes
           write(unt_ascii, '(2i6,2es24.15)') i, j, cmplx(merge(one, zero, i == j), zero, kind=dp)
         end do
       end do
     end if

     ! Record 10: eigenvectors divided by masses, i.e. phonon displacements (dyn).
     write(unt) ((dyn_qe(i, j), i=1,nmodes), j=1,nmodes)
     if (ascii_write) then
       write(unt_ascii, '(a)') "# Record 10: displacements dyn(component, mode)"
       do j=1,nmodes
         do i=1,nmodes
           write(unt_ascii, '(2i6,2es24.15)') i, j, dyn_qe(i, j)
         end do
       end do
     end if

     ! Records 11+: per k-point, the k-point (reduced coords), band energies (Ry) and,
     !              for each Cartesian atomic perturbation, the e-ph matrix elements g (Ry).
     do my_ik=1,gqk%my_nk
       kk_bz = gqk%my_kpts(:, my_ik)
       ik_ibz = gqk%my_k2ibz(1, my_ik)
       kk_ibz = ebands%kptns(:,ik_ibz)
       ik_glob = my_ik + gqk%my_kstart - 1

       write(unt) kk_bz
       write(unt) (ebands%eig(ib, ik_ibz, spin) * Ha2Ry, ib=gqk%bstart_k, gqk%bstop_k)
       if (ascii_write) then
         write(unt_ascii, '(a,i0,a)') "# k-point ", ik_glob, " (reduced coordinates)"
         write(unt_ascii, '(3es24.15)') kk_bz
         write(unt_ascii, '(a)') "#   band energies (Ry)"
         write(unt_ascii, '(*(es24.15,1x))') (ebands%eig(ib, ik_ibz, spin) * Ha2Ry, ib=gqk%bstart_k, gqk%bstop_k)
       end if

       ! Build the Cartesian deformation potential expected by EPIQ:
       !   d_cart(mu) = M_kappa * sum_nu [ g_mode(nu) * sqrt(2 w_nu) ] * conjg(displ_cart(mu,nu))
       ! This inverts ABINIT's atom->mode transform g_mode(nu) = sum_mu d_cart(mu) e(mu,nu)/sqrt(2 w_nu)
       ! (see m_ephtk/ephtk_gkknu_from_atm), recovering the bare <k+q,m| dV/du^cart_mu |k,n>:
       ! the sqrt(2 w_nu) and mass factors exactly cancel ABINIT's normalization.
       ! mu = idir + 3*(iat-1); M_kappa = amu(iat) * amu_emass (electron-mass units).
       ! Modes with w_nu <= 0 (acoustic/imaginary) carry g_mode = 0 and are skipped.
       g_cart = (zero, zero)
       do nu=1,nmodes
         if (gqk%my_wnuq(nu, my_iq) <= tol6) cycle
         wfact = sqrt(two * gqk%my_wnuq(nu, my_iq))
         do mu=1,nmodes
           g_cart(mu, :, :) = g_cart(mu, :, :) &
             + (wfact * conjg(dyn_qe(mu, nu))) * gqk%my_g(nu, :, my_iq, :, my_ik)
         end do
       end do
       do mu=1,nmodes
         iat = (mu - 1) / 3 + 1
         g_cart(mu, :, :) = g_cart(mu, :, :) * (cryst%amu(cryst%typat(iat)) * amu_emass)
       end do

       do mu=1,nmodes
         ! d_matrix(m, n) = <k+q, m| dV/du^cart_mu |k, n>, with the k+q (bra) band index m
         ! running fastest, matching QE's el_ph_mat(jbnd, ibnd) write order.
         ! The 2nd dim (nb_kq) is the k+q band, the 3rd dim (nb_k) the k band.
         write(unt) ((g_cart(mu, j, i) * Ha2Ry, j=1,gqk%nb_kq), i=1,gqk%nb_k)
         if (ascii_write) then
           ! Columns: k index, Cartesian perturbation, band(k+q), band(k), Re(d), Im(d) in Ry.
           do i=1,gqk%nb_k
             band_k = i + gqk%bstart_k - 1
             do j=1,gqk%nb_kq
               band_kq = j + gqk%bstart_kq - 1
               write(unt_ascii, '(4i6,2es30.15)') &
                 ik_glob, mu, band_kq, band_k, g_cart(mu, j, i) * Ha2Ry
             end do
           end do
         end if
       end do
     end do ! my_ik

     ! Trailing block: symmetry operations and star of q (QE elphsum_wannier layout).
     call write_qe_symmetry(cryst, qpt, celldm(1), unt, ascii_write, unt_ascii)

     close(unt)
     if (ascii_write) close(unt_ascii)
   end do ! my_iq

   ABI_FREE(g_cart)
   end associate
 end do ! spin

 ! TODO:
 ! Output BECS, dynamical quadrupoles, dynamical matrix, group velocities
 ! See m_ifc
 !ifc%zeff
 !fc%qdrp_cart
 !ifc%eta
 !ifc%rpt
 !ifc%wghatm
 !ifc%dynmat
 !ifc%short_atmfrc

 ABI_FREE(dyn_qe)
 ABI_FREE(tau_cart)

 call gstore%free()

end subroutine gstore_convert
!!***

!----------------------------------------------------------------------

!!****f* m_gstore_converters/write_epiq_input
!! NAME
!! write_epiq_input
!!
!! FUNCTION
!!  Write the EPIQ input file with the &Diff_Start_Param namelist and the KPOINTS
!!  section, reproducing QE's print_ph_input2epiq (called from ep_matrix_element_wannier).
!!  The SCF/DFPT parameters are taken from the ABINIT dataset and band structure.
!!
!! INPUTS
!!  dtset<dataset_type>=input variables (occopt, tsmear, k-grid).
!!  ebands<ebands_t>=band structure (Fermi level, nelect, k-points, gaps).
!!  fname=name of the output file.
!!
!! SOURCE

subroutine write_epiq_input(dtset, ebands, fname)

!Arguments ------------------------------------
 type(dataset_type),intent(in) :: dtset
 class(ebands_t),intent(in) :: ebands
 character(len=*),intent(in) :: fname

!Local variables-------------------------------
!scalars
 integer :: unt, ik, ngauss, gap_err, nk1, nk2, nk3, k1, k2, k3
 logical :: is_metal, automatic
 real(dp) :: homo, lumo
 character(len=500) :: msg
 character(len=24) :: smear_label
 type(gaps_t) :: gaps
!----------------------------------------------------------------------

 if (open_file(fname, msg, newunit=unt, form="formatted", status="unknown", action="write") /= 0) then
   ABI_ERROR(msg)
 end if

 ! occopt >= 3 => metallic occupation with smearing; otherwise fixed occupations (insulator).
 is_metal = dtset%occopt >= 3

 write(unt, '(a)') "! parameter of the SCF DFPT calculation useful for EPIq"
 write(unt, '(a)') "&Diff_Start_Param"
 write(unt, '(3x,a,f12.6,a)') "efermi=", ebands%fermie * Ha_eV, ", ! in (eV)"
 write(unt, '(3x,a,f12.6,a)') "nel_r=", ebands%nelect, ","

 if (.not. is_metal) then
   ! Insulator: report the HOMO and LUMO levels (in eV).
   gaps = ebands%get_gaps(gap_err)
   if (gap_err == 0) then
     homo = gaps%vb_max(1); lumo = gaps%cb_min(1)
   else
     ! Could not determine a gap (semimetal?): fall back to the Fermi level.
     homo = ebands%fermie; lumo = ebands%fermie
   end if
   call gaps%free()
   write(unt, '(3x,a,f12.6,a)') "homo=", homo * Ha_eV, ", ! in (eV)"
   write(unt, '(3x,a,f12.6,a)') "lumo=", lumo * Ha_eV, ", ! in (eV)"
 else
   ! Metal: report the smearing width (Rydberg) and the QE ngauss code.
   ! Map ABINIT occopt onto QE ngauss (see Modules input conventions):
   !   3 -> -99 (Fermi-Dirac), 4/5 -> -1 (cold/Marzari), 6 -> 1 (Methfessel-Paxton), 7 -> 0 (Gaussian)
   select case (dtset%occopt)
   case (3);        ngauss = -99; smear_label = "fd"
   case (4, 5);     ngauss =  -1; smear_label = "cold"
   case (6);        ngauss =   1; smear_label = "mp"
   case (7);        ngauss =   0; smear_label = "gauss"
   case default;    ngauss = -66; smear_label = "unknown"
   end select
   write(unt, '(3x,a,f12.6,a)') "sigma_ph=", dtset%tsmear * two, ", ! in (Rydberg)"
   write(unt, '(3x,a,i3,a)') "ngauss_ph=", ngauss, ", ! "//trim(smear_label)
 end if
 write(unt, '(a)') "/"

 ! KPOINTS section. Use the "automatic" (Monkhorst-Pack) form when the k-mesh is a
 ! diagonal kptrlatt with a single shift; otherwise dump the explicit list.
 automatic = (dtset%nshiftk == 1) .and. &
   all([dtset%kptrlatt(1,2), dtset%kptrlatt(1,3), dtset%kptrlatt(2,1), &
        dtset%kptrlatt(2,3), dtset%kptrlatt(3,1), dtset%kptrlatt(3,2)] == 0) .and. &
   dtset%kptrlatt(1,1) > 0 .and. dtset%kptrlatt(2,2) > 0 .and. dtset%kptrlatt(3,3) > 0

 write(unt, '(/,a)') "KPOINTS"
 if (automatic) then
   nk1 = dtset%kptrlatt(1,1); nk2 = dtset%kptrlatt(2,2); nk3 = dtset%kptrlatt(3,3)
   ! ABINIT shiftk are fractional shifts; a half-grid shift (0.5) maps to the QE flag 1.
   k1 = nint(two * dtset%shiftk(1,1)); k2 = nint(two * dtset%shiftk(2,1)); k3 = nint(two * dtset%shiftk(3,1))
   write(unt, '(a)') "automatic"
   write(unt, '(3(2x,i6),3x,3(2x,i3))') nk1, nk2, nk3, k1, k2, k3
 else
   ! Explicit list in reduced (crystal) coordinates.
   write(unt, '(a)') "crystal"
   write(unt, '(6x,i9)') ebands%nkpt
   do ik=1,ebands%nkpt
     write(unt, '(3x,4(es20.10,2x))') ebands%kptns(:,ik), ebands%wtk(ik)
   end do
 end if

 close(unt)

end subroutine write_epiq_input
!!***

!----------------------------------------------------------------------

!!****f* m_gstore_converters/write_dynq0
!! NAME
!! write_dynq0
!!
!! FUNCTION
!!  Write the EPIQ "dynq0" file: the q-mesh, the number of irreducible q-points
!!  contained in the elph.mat files, and the list of those q-points in Cartesian
!!  coordinates (units of 2pi/alat, i.e. QE tpiba). Analog of QE's dynq0 output.
!!
!! INPUTS
!!  ngqpt(3)=dimensions of the q-mesh.
!!  qibz(3,nqibz)=irreducible q-points in reduced (crystal) coordinates.
!!  gprimd(3,3)=reciprocal lattice vectors (Bohr^-1), columns G_i/2pi (ABINIT convention).
!!  alat=lattice parameter in Bohr (celldm(1)).
!!  fname=name of the output file.
!!
!! SOURCE

subroutine write_dynq0(ngqpt, qibz, gprimd, alat, fname)

!Arguments ------------------------------------
 integer,intent(in) :: ngqpt(3)
 real(dp),intent(in) :: qibz(:,:), gprimd(3,3), alat
 character(len=*),intent(in) :: fname

!Local variables-------------------------------
 integer :: unt, iq, nqibz
 real(dp) :: qcart(3)
 character(len=500) :: msg
!----------------------------------------------------------------------

 nqibz = size(qibz, 2)

 if (open_file(fname, msg, newunit=unt, form="formatted", status="unknown", action="write") /= 0) then
   ABI_ERROR(msg)
 end if

 ! Line 1: q-mesh. Line 2: number of irreducible q-points.
 write(unt, '(3i4)') ngqpt(1), ngqpt(2), ngqpt(3)
 write(unt, '(i4)') nqibz

 ! One line per irreducible q-point in Cartesian coordinates (2pi/alat units):
 ! q_cart[tpiba] = alat * matmul(gprimd, q_red).
 do iq=1,nqibz
   qcart = alat * matmul(gprimd, qibz(:,iq))
   write(unt, '(3e24.15)') qcart(1), qcart(2), qcart(3)
 end do

 close(unt)

end subroutine write_dynq0
!!***

!----------------------------------------------------------------------

!!****f* m_gstore_converters/write_dynq
!! NAME
!! write_dynq
!!
!! FUNCTION
!!  Write a Quantum ESPRESSO dynamical-matrix file (dynq<iq>) for a single q-point,
!!  reproducing the layout of QE's write_dyn_on_file + dyndiag: header (cell, atoms),
!!  the dynamical matrix in Cartesian axes, and the diagonalization (frequencies and
!!  eigenvectors).
!!
!!  Frequencies, displacements and orthonormal eigenvectors are obtained from ifc%fourq
!!  (the gauge is irrelevant here since the dynamical matrix is gauge-invariant). The
!!  Cartesian dynamical matrix is rebuilt as
!!    phi(ka,k'b) = sqrt(M_k M_k') * sum_nu z(ka,nu) * w2(nu) * conjg(z(k'b,nu))
!!  with z the orthonormal eigenvectors, w2 = signed omega^2 in Ry^2 and M the QE
!!  Rydberg atomic masses (amu * amu_emass/2). This matches QE's convention.
!!
!! INPUTS
!!  cryst<crystal_t>=crystal structure.
!!  ifc<ifc_type>=interatomic force constants (for Fourier interpolation at q).
!!  qpt_red(3)=q-point in reduced (crystal) coordinates.
!!  alat=lattice parameter in Bohr (celldm(1)).
!!  atm(ntypat)=atomic symbols.
!!  fname=name of the output file.
!!
!! SOURCE

subroutine write_dynq(cryst, ifc, qpt_red, alat, atm, fname)

!Arguments ------------------------------------
 class(crystal_t),intent(in) :: cryst
 class(ifc_type),intent(in) :: ifc
 real(dp),intent(in) :: qpt_red(3), alat
 character(len=3),intent(in) :: atm(cryst%ntypat)
 character(len=*),intent(in) :: fname

!Local variables-------------------------------
!scalars
 integer :: natom, nmodes, ntypat, unt, na, nb, icar, jcar, nu, it
 real(dp),parameter :: Ha2Ry = two
 real(dp) :: znorm, freq_cm, freq_thz
 complex(dp) :: zi, zj, cs
 character(len=500) :: msg
!arrays
 real(dp) :: qcart(3), at(3,3), celldm(6)
 real(dp),allocatable :: phfrq(:), displ_cart(:,:,:,:), eigvec(:,:,:,:), w2(:), amass_qe(:), tau(:,:)
 complex(dp),allocatable :: phi(:,:,:,:)
!----------------------------------------------------------------------

 natom = cryst%natom; nmodes = 3 * natom; ntypat = cryst%ntypat

 ABI_MALLOC(phfrq, (nmodes))
 ABI_MALLOC(displ_cart, (2, 3, natom, nmodes))
 ABI_MALLOC(eigvec, (2, 3, natom, nmodes))
 ABI_MALLOC(w2, (nmodes))
 ABI_MALLOC(amass_qe, (ntypat))
 ABI_MALLOC(tau, (3, natom))
 ABI_MALLOC(phi, (3, 3, natom, natom))

 ! Phonon frequencies, displacements and orthonormal eigenvectors at this q.
 call ifc%fourq(cryst, qpt_red, phfrq, displ_cart, out_eigvec=eigvec)

 ! QE masses in Rydberg atomic units (amu * amu_ry, amu_ry = amu_emass/2).
 do it=1,ntypat
   amass_qe(it) = cryst%amu(it) * amu_emass * half
 end do

 ! Signed squared phonon frequencies in Ry^2.
 do nu=1,nmodes
   w2(nu) = phfrq(nu) * abs(phfrq(nu)) * Ha2Ry ** 2
 end do

 ! Cartesian quantities (alat / tpiba units).
 qcart = alat * matmul(cryst%gprimd, qpt_red)
 do na=1,natom
   tau(:,na) = matmul(cryst%rprimd, cryst%xred(:,na)) / alat
 end do
 do nu=1,3
   at(:,nu) = cryst%rprimd(:,nu) / alat
 end do
 celldm = zero; celldm(1) = alat

 ! Dynamical matrix in Cartesian axes (QE convention, see header).
 do nb=1,natom
   do na=1,natom
     do jcar=1,3
       do icar=1,3
         cs = czero
         do nu=1,nmodes
           zi = cmplx(eigvec(1,icar,na,nu), eigvec(2,icar,na,nu), kind=dp)
           zj = cmplx(eigvec(1,jcar,nb,nu), eigvec(2,jcar,nb,nu), kind=dp)
           cs = cs + zi * w2(nu) * conjg(zj)
         end do
         phi(icar,jcar,na,nb) = sqrt(amass_qe(cryst%typat(na)) * amass_qe(cryst%typat(nb))) * cs
       end do
     end do
   end do
 end do

 if (open_file(fname, msg, newunit=unt, form="formatted", status="unknown", action="write") /= 0) then
   ABI_ERROR(msg)
 end if

 ! ---- Header ----
 write(unt, '(a)') "Dynamical matrix file"
 write(unt, '(a)') "Converted from ABINIT GSTORE"
 ! ntyp, nat, ibrav=0 (free lattice) followed by celldm; with ibrav=0 the basis vectors follow.
 write(unt, '(i3,i5,i4,6f11.7)') ntypat, natom, 0, (celldm(it), it=1,6)
 write(unt, '(a)') "Basis vectors"
 do nu=1,3
   write(unt, '(2x,3f15.9)') at(1,nu), at(2,nu), at(3,nu)
 end do
 do it=1,ntypat
   write(unt, *) it, " '"//atm(it)//"' ", amass_qe(it)
 end do
 do na=1,natom
   write(unt, '(2i5,3f18.10)') na, cryst%typat(na), tau(1,na), tau(2,na), tau(3,na)
 end do

 ! ---- Dynamical matrix in Cartesian axes ----
 write(unt, '(/,5x,a)') "Dynamical  Matrix in cartesian axes"
 write(unt, '(/,5x,a,3f14.9,a,/)') "q = ( ", qcart(1), qcart(2), qcart(3), " ) "
 do na=1,natom
   do nb=1,natom
     write(unt, '(2i5)') na, nb
     do icar=1,3
       write(unt, '(3(2f12.8,2x))') (phi(icar,jcar,na,nb), jcar=1,3)
     end do
   end do
 end do

 ! ---- Diagonalization: frequencies and eigenvectors ----
 write(unt, '(/,5x,a)') "Diagonalizing the dynamical matrix"
 write(unt, '(/,5x,a,3f14.9,a,/)') "q = ( ", qcart(1), qcart(2), qcart(3), " ) "
 write(unt, '(1x,74("*"))')
 do nu=1,nmodes
   freq_cm = phfrq(nu) * Ha_cmm1
   freq_thz = phfrq(nu) * Ha_THz
   write(unt, '(5x,a,i5,a,f15.6,a,f15.6,a)') "freq (", nu, ") = ", freq_thz, " [THz] = ", freq_cm, " [cm-1]"
   znorm = zero
   do na=1,natom
     do icar=1,3
       znorm = znorm + eigvec(1,icar,na,nu)**2 + eigvec(2,icar,na,nu)**2
     end do
   end do
   znorm = sqrt(znorm); if (znorm < tol12) znorm = one
   do na=1,natom
     write(unt, '(1x,a,3(f10.6,1x,f10.6,3x),a)') "( ", &
       (eigvec(1,icar,na,nu)/znorm, eigvec(2,icar,na,nu)/znorm, icar=1,3), ")"
   end do
 end do
 write(unt, '(1x,74("*"))')

 close(unt)

 ABI_FREE(phfrq)
 ABI_FREE(displ_cart)
 ABI_FREE(eigvec)
 ABI_FREE(w2)
 ABI_FREE(amass_qe)
 ABI_FREE(tau)
 ABI_FREE(phi)

end subroutine write_dynq
!!***

!----------------------------------------------------------------------

!!****f* m_gstore_converters/write_qe_symmetry
!! NAME
!! write_qe_symmetry
!!
!! FUNCTION
!!  Append the symmetry-operations block and the star of q to an already open
!!  elph.mat.q_* file, reproducing the layout written by Quantum ESPRESSO's
!!  elphsum_wannier (PHonon/PH/elphon.f90).
!!
!!  The QE quantities are reconstructed from the ABINIT crystal object using the
!!  following correspondence (see PW/src/symm_base.f90, PHonon/PH/obsolete.f90 and
!!  LR_Modules/star_q.f90 in the bundled QE sources):
!!
!!    QE s(:,:,isym)   = cryst%symrec(:,:,isym)   (rotations on reduced reciprocal coords;
!!                                                 s^T = symrel^-1 acts on reduced positions)
!!    QE irt(isym,na)  = cryst%indsym(4,isym,na)  (atom na -> atom irt under symrel^-1)
!!    QE at(:,j)       = cryst%rprimd(:,j) / alat (direct lattice, alat units)
!!    QE bg(:,j)       = cryst%gprimd(:,j) * alat (reciprocal lattice, 2pi/alat units)
!!    QE tau / xau     = cryst%xred               (reduced atomic positions)
!!
!!  invs, rtau and the star (nq, sxq, isq, imq) are then obtained by porting QE's
!!  inverse_s, sgam_ph and star_q. The symmetry operations are written in the native
!!  ABINIT order (1..nsym); only the order differs from QE (which sorts the small
!!  group of q first), the set and the per-operation data are equivalent.
!!
!! INPUTS
!!  cryst<crystal_t>=crystal structure (with symmetries).
!!  qpt(3)=q-point in reduced (crystal) coordinates.
!!  alat=lattice parameter in Bohr (celldm(1)).
!!  unt=Fortran unit of the (open) unformatted file.
!!  ascii_write=if .True. also mirror the records in the formatted file unt_ascii.
!!  unt_ascii=Fortran unit of the (open) formatted file (used only if ascii_write).
!!
!! SOURCE

subroutine write_qe_symmetry(cryst, qpt, alat, unt, ascii_write, unt_ascii)

!Arguments ------------------------------------
!scalars
 class(crystal_t),intent(in) :: cryst
 integer,intent(in) :: unt, unt_ascii
 logical,intent(in) :: ascii_write
 real(dp),intent(in) :: alat
!arrays
 real(dp),intent(in) :: qpt(3)

!Local variables-------------------------------
!scalars
 integer :: nsym, natom, isym, jsym, ism1, i, j, k, na, nb, nq, imq, iq
 logical :: found
 real(dp),parameter :: accep = 1.0e-5_dp
!arrays
 integer :: ss(3,3)
 integer,parameter :: identity(3,3) = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])
 integer,allocatable :: s(:,:,:), invs(:), irt(:,:), isq(:), nsq(:)
 real(dp) :: at(3,3), bg(3,3), aq(3), raq(3), ft(3), dq(3)
 real(dp),allocatable :: rtau(:,:,:), sxq(:,:), saq(:,:)
!----------------------------------------------------------------------

 nsym = cryst%nsym; natom = cryst%natom

 ABI_MALLOC(s, (3, 3, nsym))
 ABI_MALLOC(invs, (nsym))
 ABI_MALLOC(irt, (nsym, natom))
 ABI_MALLOC(isq, (nsym))
 ABI_MALLOC(nsq, (nsym))
 ABI_MALLOC(rtau, (3, nsym, natom))
 ABI_MALLOC(sxq, (3, nsym))
 ABI_MALLOC(saq, (3, nsym))

 ! Lattice vectors in QE units (at . bg^T = identity).
 do j=1,3
   at(:,j) = cryst%rprimd(:,j) / alat
   bg(:,j) = cryst%gprimd(:,j) * alat
 end do

 ! Rotations in crystal axis and atom mapping.
 do isym=1,nsym
   s(:,:,isym) = cryst%symrec(:,:,isym)
   do na=1,natom
     irt(isym, na) = cryst%indsym(4, isym, na)
   end do
 end do

 ! invs(isym): index of the inverse operation (ported from QE inverse_s).
 do isym=1,nsym
   found = .False.
   do jsym=1,nsym
     ss = matmul(s(:,:,jsym), s(:,:,isym))
     if (all(ss == identity)) then
       invs(isym) = jsym; found = .True.; exit
     end if
   end do
   ABI_CHECK(found, "write_qe_symmetry: symmetry operations do not form a group.")
 end do

 ! rtau(:,isym,na) = S.tau_na - tau_nb in Cartesian coords (alat units), with nb = irt.
 ! Ported from QE sgam_ph using xau = reduced atomic coordinates = cryst%xred.
 rtau = zero
 do isym=1,nsym
   do na=1,natom
     nb = irt(isym, na)
     do i=1,3
       ft(i) = s(1,i,isym) * cryst%xred(1,na) + s(2,i,isym) * cryst%xred(2,na) &
             + s(3,i,isym) * cryst%xred(3,na) - cryst%xred(i,nb)
     end do
     do i=1,3
       rtau(i, isym, na) = at(i,1) * ft(1) + at(i,2) * ft(2) + at(i,3) * ft(3)
     end do
   end do
 end do

 ! Star of q (ported from QE star_q). aq is q in reduced (crystal) coordinates.
 aq(:) = qpt(:)
 nsq(:) = 0; isq(:) = 0; saq(:,:) = zero; sxq(:,:) = zero; nq = 0
 do isym=1,nsym
   ism1 = invs(isym)
   do i=1,3
     raq(i) = s(i,1,ism1) * aq(1) + s(i,2,ism1) * aq(2) + s(i,3,ism1) * aq(3)
   end do
   do iq=1,nq
     dq(:) = raq(:) - saq(:,iq)
     if (all(abs(dq - nint(dq)) < accep)) then
       isq(isym) = iq; nsq(iq) = nsq(iq) + 1
     end if
   end do
   if (isq(isym) == 0) then
     nq = nq + 1; nsq(nq) = 1; isq(isym) = nq; saq(:,nq) = raq(:)
     do i=1,3
       sxq(i,nq) = bg(i,1) * saq(1,nq) + bg(i,2) * saq(2,nq) + bg(i,3) * saq(3,nq)
     end do
   end if
 end do

 ! imq: index of -q in the star (0 if absent).
 imq = 0
 do iq=1,nq
   dq(:) = -aq(:) - saq(:,iq)
   if (all(abs(dq - nint(dq)) < accep)) imq = iq
 end do

 ! Sanity check on the star degeneracy (as in QE star_q).
 do iq=1,nq
   if (nsq(iq) * nq /= nsym) then
     ABI_WARNING(sjoin("write_qe_symmetry: unexpected star-of-q degeneracy for iq=", itoa(iq)))
   end if
 end do

 ! ---- Write the block in the QE elphsum_wannier order ----
 do j=1,3
   write(unt) (at(i,j), i=1,3)
 end do
 do j=1,3
   write(unt) (bg(i,j), i=1,3)
 end do
 write(unt) nsym, nq, imq
 do i=1,nsym
   write(unt) i, invs(i), isq(i)
   do j=1,3
     do k=1,3
       write(unt) k, j, s(k,j,i)
     end do
   end do
   do j=1,natom
     write(unt) j, irt(i,j)
   end do
   do j=1,3
     do k=1,natom
       write(unt) j, i, rtau(j,i,k)
     end do
   end do
   do j=1,3
     write(unt) j, sxq(j,i)
   end do
 end do

 ! ---- Optional human-readable mirror ----
 if (ascii_write) then
   write(unt_ascii, '(a)') "# Symmetry: direct lattice vectors at(:,j) (alat units)"
   do j=1,3
     write(unt_ascii, '(3es24.15)') (at(i,j), i=1,3)
   end do
   write(unt_ascii, '(a)') "# Symmetry: reciprocal lattice vectors bg(:,j) (2pi/alat units)"
   do j=1,3
     write(unt_ascii, '(3es24.15)') (bg(i,j), i=1,3)
   end do
   write(unt_ascii, '(a)') "# Symmetry: nsym, nq, imq"
   write(unt_ascii, '(3i8)') nsym, nq, imq
   do i=1,nsym
     write(unt_ascii, '(a,3i6)') "# isym, invs, isq: ", i, invs(i), isq(i)
     write(unt_ascii, '(a)') "#   rotation s(row k, col j)"
     do j=1,3
       do k=1,3
         write(unt_ascii, '(3i6)') k, j, s(k,j,i)
       end do
     end do
     write(unt_ascii, '(a)') "#   irt(atom)"
     do j=1,natom
       write(unt_ascii, '(2i6)') j, irt(i,j)
     end do
     write(unt_ascii, '(a)') "#   rtau(coord, isym, atom) (alat units)"
     do j=1,3
       do k=1,natom
         write(unt_ascii, '(2i6,es24.15)') j, i, rtau(j,i,k)
       end do
     end do
     write(unt_ascii, '(a)') "#   sxq(coord) (2pi/alat units)"
     do j=1,3
       write(unt_ascii, '(i6,es24.15)') j, sxq(j,i)
     end do
   end do
 end if

 ABI_FREE(s)
 ABI_FREE(invs)
 ABI_FREE(irt)
 ABI_FREE(isq)
 ABI_FREE(nsq)
 ABI_FREE(rtau)
 ABI_FREE(sxq)
 ABI_FREE(saq)

end subroutine write_qe_symmetry
!!***

!----------------------------------------------------------------------

! Helper function copied from epic/src/io_matelem.F90
subroutine define_band_string(index, string, lstr)
  integer,intent(in) :: index
  integer,intent(out) :: lstr
  character(len=3),intent(out) :: string

! here put a check on the string length

  string=' '
  if(index < 10) then
     WRITE( string(1:1), '(I1)' ) index
     lstr=1
  elseif(index < 100) then
     WRITE( string(1:2), '(I2)' ) index
     lstr=2
  elseif(index < 1000) then
     WRITE( string(1:3), '(I3)' ) index
     lstr=3
  endif

  string=trim(adjustl(string))
end subroutine define_band_string

end module m_gstore_converters
