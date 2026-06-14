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
 use m_ebands,         only : ebands_t
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
 logical :: with_g2dw, q_is_gamma, lborn
 real(dp),parameter :: Ha2Ry = two
 real(dp) :: weight_qq
 character(len=5000) :: msg
 character(len=abi_slen) :: with_gmode, gvals_name, gtype
 character(len=fnlen) :: fname, elphmat_dir, prefix
 character(len=3) :: band_i
 type(gstore_t) :: gstore
!arrays
 integer :: units(2)
 real(dp) :: qpt(3), kk_bz(3), kk_ibz(3)
 character(len=3) :: atm(cryst%ntypat)
!----------------------------------------------------------------------

 ! Only master works here as performance is not crucial.
 ! In principle one could activate the q-point/spin parallelism just to distribute
 ! the memory for the g's to avoid going OOM.
 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 if (my_rank /= 0) return

 units = [std_out, ab_out]
 natom = cryst%natom; nmodes = 3 * cryst%natom; nsppol = ebands%nsppol; this_comm = xmpi_comm_self

 ! Preliminary consistency check.
 call wrtout(units, sjoin(" Begin conversion GSTORE --> ", dtset%gstore_convert))
 ABI_CHECK(gstore_path /= ABI_NOFILE, sjoin("Invalid gstore_path:", gstore_path))
 ABI_CHECK(dtset%gstore_convert == "epiq", "only gstore_convert == 'epiq' is supported.")
 ABI_CHECK(.not. cryst%isalchemical(), "'epiq' format does not support alchemical pseudos.")

 ! Read g(k,q) from GSTORE and store them in gqk%my_g(nu, ib_kq, my_iq, in_k, ik)
 !
 ! (my_npert, nb_kq, my_nq, nb_k, my_nk)
 ! (       p, b1_kq,     q, b2_k, k)  -->  <k+q, b1| D_{q,p}H |k, b2>
 !gqk%my_g(nu, ib_kq, my_iq, in_k, ik)
 !
 ! The g are complex and in the phonon representation.
 ! All quantities are in atomic units (Hartree and Bohr).
 !
 ! The gstore file produced by GWPT has both GWPT and KS g.
 ! By default, we convert the GWPT matrix elements but one can still select
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
 ABI_CHECK(gstore%same_nbands(msg), msg)
 if (gstore%check_cplex_qkzone_gmode(2, "bz", "bz", "phonon", kfilter="none") /= 0) then
   ABI_ERROR("GSTORE.nc should have both k and q in the full BZ. See messages above.")
 end if

 ! Create directory to host output files.
 prefix = "epiq"
 elphmat_dir = strcat(dtfil%filnam_ds(4), "_", prefix)
 call wrtout(units, sjoin(" Output files written to directory:", elphmat_dir))
 call execute_command_line(sjoin("rm -rf", elphmat_dir), exitstat=ierr)
 call clib_mkdir_if_needed(elphmat_dir, ierr)
 ABI_CHECK(ierr == 0, "mkdir returned ierr /= 0")

 ! NB: atm is character(len=3) while symbol_type returns character(len=2).
 do itypat=1, cryst%ntypat
   atm(itypat)(1:2) = cryst%symbol_type(itypat)
   atm(itypat)(3:3) = ""
 end do

 ! TODO: Need helper function to get ibrav, celldm from crystal
 ! Some help from the EPIC developers would be great.
 !call cryst%get_ibrav_celldm(ibrav, celldm)

 ! Loop over collinear spins.
 do my_is=1,gstore%my_nspins
   spin = gstore%my_spins(my_is)
   associate (gqk => gstore%gqk(my_is))

   !num_bands = gqk%nb_k

   ! Loop over q-points: MG TODO: Here I assume the q-points are in the BZ, right?
   do my_iq=1, gqk%my_nq
     iq_glob = my_iq + gqk%my_qstart - 1

     call gqk%myqpt(my_iq, gstore, weight_qq, qpt); q_is_gamma = sum(qpt**2) < tol14
     call define_band_string(iq_glob, band_i, lstr_j)

     fname = trim(elphmat_dir)//"/"//trim(prefix)//'_elph.mat.q_'//band_i(1:lstr_j)
     if (open_file(fname, msg, newunit=unt, form="unformatted", status="unknown", action="write") /= 0) then
       ABI_ERROR(msg)
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

     ! MG TODO:
     ! - I assume w2 is omega and not omega^2, right?
     ! - Is lborn used, how do you read BECS, dynamical quadrupoles?
     ! - I assume tau are atom positions in reduced coords.
     ! - I assume nel_aux is the number of electrons including possible doping (real variable)
     ! - I don't know the conventions for zz and I assume dyn are in Cartesian coords in Bohr.
     lborn = .False.
     write(unt) qpt
     write(unt) ebands%nspinor == 2, nsppol, lborn
     write(unt) ebands%nelect
     write(unt) gqk%bstart_k, gqk%bstop_k, gqk%bstop_k - gqk%bstart_k + 1
     write(unt) nmodes, gqk%glob_nk, natom, cryst%ntypat
     !write(unt) ibrav, (celldm(j), j=1,6) ! FIXME
     ! amu are the mass of the atoms (atomic mass unit).
     write(unt) (atm(j), j=1,cryst%ntypat), (cryst%amu(j), j=1,cryst%ntypat), &
                (cryst%typat(j), j=1,natom), ((cryst%xred(j,i), j=1,3), i=1,natom)
     write(unt) (gqk%my_wnuq(nu, my_iq) * Ha2Ry, nu=1,nmodes)
     !read(unt) ((zz(i,j,iqph), i=1,nmodes),j=1,nmodes)   !eigenvectors in the QE basis ! FIXME
     write(unt) (gqk%my_displ_cart(:,:,:,:,my_iq))

     do my_ik=1,gqk%my_nk
       kk_bz = gqk%my_kpts(:, my_ik)
       ik_ibz = gqk%my_k2ibz(1, my_ik)
       kk_ibz = ebands%kptns(:,ik_ibz)

       write(unt) kk_bz
       write(unt) (ebands%eig(ib, ik_ibz, spin) * Ha2Ry, ib=gqk%bstart_k, gqk%bstop_k)
       do nu=1,nmodes
         ! MG FIXME: Here I need to know how bands are ordered. Is it (m, n) or (n, m)?
         ! The code below Assumes (m, n).
         !write(unt) ((g_matrix(j, i, nu, k, iq),j=1,num_bands),i=1,num_bands)
         write(unt) ((gqk%my_g(nu, j, my_iq, i, my_ik) * Ha2Ry, i=1,gqk%nb_k), j=1,gqk%nb_k)
       end do
     end do ! my_ik

     close(unt)
   end do ! my_iq
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

 call gstore%free()

end subroutine gstore_convert
!!***

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

