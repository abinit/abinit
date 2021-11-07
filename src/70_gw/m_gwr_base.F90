!!****m* ABINIT/m_gwr_base
!! NAME
!!  m_gwr_base
!!
!! FUNCTION
!!  Base objects and tools to implement GW in real-space.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2021 ABINIT group (MG)
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

module m_gwr_base

 use defs_basis
 !use m_abicore
 use m_errors
 use m_xmpi
 use m_slk
 !use m_xomp
 !use m_sort

 !use defs_datatypes, only : ebands_t
 use m_fstrings, only : sjoin, itoa
 use m_crystal,  only : crystal_t
 use m_dtset,    only : dataset_type
 !use m_fft_core, only : get_kg
 use m_fft,      only : fft_ug

 implicit none

 private

 !public ::
!!***

!----------------------------------------------------------------------

!!****t* m_gwr_base/grdesc_t
!! NAME
!! gedesc_t
!!
!! FUNCTION
!!
!! SOURCE

 type,public :: grdesc_t

   !integer :: istwfk
   ! Storage mode for this k point.

   integer :: npw = -1
   ! Number of plane-waves for this k-point.

   integer :: mgfft = -1
   ! maximum size of 1D FFTs

   integer :: nfft = -1
   ! Number of fft points

   integer :: ngfft(18) = -1
   ! FFT mesh used to transform G(g,g')

   integer,allocatable :: gbound(:,:)
   !  gbound(2*green_mgfft+8,2)
   !  sphere boundary info for reciprocal to real space (local)

   integer,allocatable :: gvec(:,:)
   ! gvec(3, gvec_size)
   ! For each k, the g-vectors are ordered according to |k+g|
   ! so that we can use the same array for Green, chi and Sigma_x
   ! Note that this array is global i.e. it is not MPI-distributed in G-space, only in k-space.

   integer,allocatable :: kg_k(:,:)
   ! kg_k(3,npw)
   ! G vector coordinates in reduced cordinates.

   !% real(dp) :: kpoint(3)

 contains

   !procedure :: set_ngfft => gr_desc_set_ngfft

   procedure :: free => gr_desc_free
   ! Free memory.

 end type grdesc_t
!!***

type, extends(matrix_scalapack) :: twopoint_func_t
   integer, pointer :: gvec(:,:) => null()

   type(grdesc_t), pointer :: desc => null()

   !type(matrix_scalapack),allocatable :: g_g(:)
   ! g_g(my_nab)

   !type(matrix_scalapack),allocatable :: g_r(:)
   ! g_g(my_nab)

   !type(matrix_scalapack),allocatable :: r_r(:)
   ! r_r(my_nab)

end type twopoint_func_t

!----------------------------------------------------------------------

!!****t* m_gwr_base/green_t
!! NAME
!! green_t
!!
!! FUNCTION
!!
!! SOURCE

 type,public :: green_t

   !integer :: glob_npw = - 1
   !integer :: my_npw = -1
   !integer :: max_my_npw = -1
   !integer :: glob_nfft = -1
   !integer :: my_nfft = -1
   !integer :: max_my_nfft = -1

   !integer :: ngfft(18)
   ! FFT mesh in the unit cell.

   !integer mgfft
   ! maximum size of 1D FFTs (local)

   ! gboundin(2*mgfft+8,2)
   !  sphere boundary info for reciprocal to real space (local)

   !integer :: cplex
   ! 1 if denpot is real, 2 if complex

   !integer,allocatable :: myab2glob(:)
   ! myab2glob(my_nab)

   type(matrix_scalapack),allocatable :: g_g(:)
   ! g_g(my_nab)

   type(matrix_scalapack),allocatable :: g_r(:)
   ! g_g(my_nab)

   type(matrix_scalapack),allocatable :: r_r(:)
   ! r_r(my_nab)

   !type(grdesc_t), pointer :: desc => null()

 contains

   !procedure :: fft_gg => green_fft_gg
   !procedure :: set_ngfft => green_set_ngfft

   procedure :: free => green_free
   ! Free memory.

 end type green_t
!!***


!----------------------------------------------------------------------

!!****t* m_gwr_base/gwr_t
!! NAME
!! gwr_t
!!
!! FUNCTION
!!
!! SOURCE

 type, public :: gwr_t

   integer :: nsppol
   ! Global Number of independent spin polarizations.

   integer :: nspinor

   integer :: my_nsppol
   ! Number of independent spin polarizations treated by this MPI proc

   integer :: nkbz
   ! Number of k-points in the full BZ

   integer :: my_nkbz
   ! Number of k-points in the full BZ treated by this MPI proc

   integer :: nkibz
   ! Number of k-points in the IBZ

   integer :: ntau, niw
   ! Total number of imaginary time/frequency points

   integer :: my_ntau, my_niw
   ! Number of imaginary time/frequency points treated by this MPI rank.

   logical :: use_supercell = .True.
   ! True if we are using the supercell formalism.
   ! False if we are using the mixed-space approach with convolutions in k-space.

   !integer :: ngkpt(3)

   !integer,allocatable :: my_spin_to_glob()
   ! (my_nspin)

   !integer :: myit_to_glob(:)
   ! (my_ntau)

   !real(dp),allocatable :: tau_mesh(:)
   ! tau_mesh(ntau)

   !real(dp),allocatable :: t2w_cos_wgs(:,:)
   ! cos_weights

   !type(grdesc_t),allocatable :: green_desc_mykbz(:)
   ! green_desck_mykbz(my_nbz)

   !integer,allocatable :: gvec(:,:)
   ! gvec(3, gvec_size, my_nkbz)
   ! For each k, the gvectors are ordered according to |k+g|
   ! so that we can use the same array for G, chi and Sigma_x
   ! Note that this array is global i.e. is not MPI-distributed in G-space, only in k-space.

   !integer, allocatable :: npw_green(:)
   !integer, allocatable :: npw_chi(:)
   !integer, allocatable :: npw_sigma(:)
   ! Arrays of shape (my_nkbz) giving the number of G-vectors in gvec.
   ! These quantities are computed from the values of ecut, ecuteps and ecutsigx
   ! specified in input.

   !integer :: coords_pqbks(5) = 0
   ! Cartesian coordinates of this processor in the Cartesian grid.

   integer :: comm
   ! MPI communicator with all procs involved in the calculation

   !type(xcomm_t) :: spin_comm
   ! MPI communicator for parallelism over spins (high-level)

   !type(xcomm_t) :: k_comm
   ! MPI communicator for k-point distribution

   !type(xcomm_t) :: gr_comm
   ! MPI communicator for gr distribution

   !type(xcomm_t) :: ab_comm
   ! MPI communicator for spinor-space distribution

   !type(xcomm_t) :: tau_comm
   ! MPI communicator for imag time distribution

   !type(xcomm_t) :: omega_comm
   ! MPI communicator for imag frequency distribution

   type(dataset_type), pointer :: dtset => null()

   !type(crystal_t), pointer :: crystal => null()

   ! Tow possible approaches
   type(green_t), allocatable :: go(:,:,:)
   type(green_t), allocatable :: ge(:,:,:)
   ! (my_nkbz, my_ntau, my_nsppol * my_nab)

   type(twopoint_func_t),allocatable :: go_gg(:,:,:)
   !type(twopoint_func_t),allocatable :: go_gr(:,:,:)
   !type(twopoint_func_t),allocatable :: go_rr(:,:,:)
   type(twopoint_func_t),allocatable :: ge_gg(:,:,:)
   !type(twopoint_func_t),allocatable :: ge_gr(:,:,:)
   !type(twopoint_func_t),allocatable :: ge_rr(:,:,:)
   ! (my_nkbz, my_ntau, my_nsppol * my_nab)

   !type(matrix_scalapack), allocatable :: chiw_gg(:,:,:)
   ! (my_nqbz, my_niw, my_nsppol * my_nab)

   !character(len=fnlen) :: wfk_filepath
   ! Path to the WFK file with the KS wavefunctions.

  !real(dp),allocatable :: qbz(:,:)
  ! ! qbz(3, nqbz)
  ! ! Reduced coordinates of the q-points in the full BZ.

  !real(dp),allocatable :: qibz(:,:)
  ! ! qibz(3, nqibz)
  ! ! Reduced coordinates of the q-points in the IBZ (full simmetry of the system).

  !real(dp),allocatable :: wtq(:)
   ! wtq(nqibz)
   ! Weights of the q-points in the IBZ (normalized to one).

   !integer,allocatable:: indkk_kq(:, :)
   ! (6, %nqibz_k))
   ! Mapping k+q --> initial IBZ. Depends on ikcalc.
   ! These table used the conventions for the symmetrization of the wavefunctions expected by cgtk_rotate.
   ! In this case listkk has been called with symrel and use_symrec=False

 contains

   !procedure :: new => gwr_new

   procedure :: ggp_to_gpr  => gwr_ggp_to_gpr

   procedure :: free => gwr_free
   ! Free memory.

 end type gwr_t
!!***

contains
!!***

!----------------------------------------------------------------------

!!****f* m_gwr_base/gr_desc_new
!! NAME
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

type(grdesc_t) function gr_desc_new(ecut, kpoint, gmet) result (new)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: ecut
!arrays
 real(dp),intent(in) :: kpoint(3), gmet(3,3)
 !integer,intent(out) :: npw_k
 !integer,allocatable,intent(out) :: kg_k(:,:)

!Local variables-------------------------------
!scalars
 integer, parameter :: istwf_1 = 1

! *************************************************************************

 ! Calculate G-sphere from input ecut.
 !call get_kg(kpoint, istwf_1, ecut, gmet, new%npw_k, new%kg_k)

 !if (present(ngfft))
 !  !mgfft = maxval(ngfft(1:3))
 !  !ABI_MALLOC(gbound, (2*mgfft+8, 2))
 !  !call sphereboundary(gbound, istwf_1, kg_k, mgfft, npw_k)
 !end if

end function gr_desc_new
!!***

!----------------------------------------------------------------------

!!****f* m_gwr_base/gr_desc_free
!! NAME
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

subroutine gr_desc_free(self)

!Arguments ------------------------------------
 class(grdesc_t) :: self

! *************************************************************************

 ABI_SFREE(self%gvec)

end subroutine gr_desc_free
!!***

!----------------------------------------------------------------------

!!****f* m_gwr_base/gwr_new
!! NAME
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

type(gwr_t) function gwr_new(dtset, comm) result (gwr)

!Arguments ------------------------------------
!scalars
 type(dataset_type),target,intent(in) :: dtset
 integer,intent(in) :: comm
!arrays

!Local variables-------------------------------
!scalars
 !integer :: nkbz
 integer :: my_spin, my_it, my_ikbz, ir, my_nr
 type(grdesc_t),pointer :: desc
 !type(processor_scalapack) :: slk_processor
!arrays

! *************************************************************************

 gwr%dtset => dtset
 gwr%comm = comm
 !gwr%crystal = dtset%get_crystal(dtset, 1)

 gwr%nspinor = dtset%nspinor
 gwr%nsppol = dtset%nsppol
 gwr%my_nsppol = 1
 gwr%nkbz = 1
 gwr%nkibz = 1
 gwr%my_nkbz = 1

 gwr%ntau = 1
 gwr%my_ntau = 1
 gwr%niw = 1
 gwr%my_niw = 1

 !ABI_MALLOC(gwr%green_desc, (gwr%my_nkbz))
 !do my_ikbz=1,gwr%my_nkbz
 !  ik_bz = gwr%my_kbz_to_bz(my_ikbz)
 !  ik_ibz = gwr%my_kbz_to_ibz(my_ikbz)
 !  gwr%green_desc(my_ikbz) = gr_desc_new(ecuteff, kpoint, cryst%gmet, ngfft=green_ngfft)
 !end do

 ! Allocate G_k(g, g') in tau space
 ABI_MALLOC(gwr%go, (gwr%my_nkbz, gwr%my_ntau, gwr%my_nsppol))
 ABI_MALLOC(gwr%ge, (gwr%my_nkbz, gwr%my_ntau, gwr%my_nsppol))

 !call init_scalapack(slk_processor, gwr%gr_comm%value) # grid=[1, nprocs]
 ABI_MALLOC(gwr%go_gg, (gwr%my_nkbz, gwr%my_ntau, gwr%my_nsppol))

 !do my_it=1,gwr%my_ntau
 !  do my_ikbz=1,gwr%my_nkbz
 !  do my_spin=1,gwr%my_nsppol
 !    do my_ab=1,gwr%my_nab
 !       call gwr%go_gg, (my_ikbz, my_ntau, my_ab))
 !       col_bsize = green_npw_k / gwr%gr_comm%nproc
 !       if (mod(green_npw_k, gwr%gr_comm%nproc) /= 0) col_bsize = col_bsize + 1
 !       call init_matrix_scalapack(gg,  &
 !          green_npw_k, green_npw_k, slk_processor, 1, tbloc=[greep_npw_k, col_bsize])
 !    end do
 !  end do
 !end do

 ! Allocate chi_k(g, g') in tau space.
 !ABI_MALLOC(gwr%chi, (gwr%my_nqbz, gwr%my_ntau, gwr%my_nsppol))

 ! Prepare FFT in the supercell from R' -> G'
 !integer :: sc_ngfft(18)
 !sc_ngfft(1:3) = gwr%ngkpt * ngfft(1:3)
 !sc_nr = product(sc_ngfft(1:3))
 !call fftbox_plan3(sc_plan, sc_ngfft, fftalg, isign)
 !cal fftbox_plan3_init(plan,ndat,dims,embed,fftalg,fftcache,isign)
 !ABI_MALLOC(sc_chi, (sc_nr))
 !plan%ldxyz*plan%ndat
 !ABI_FREE(sc_chi)

 do my_spin=1,gwr%my_nsppol
   !spin = gwr%my_spin_to_glob(my_spin)
   do my_it=1,gwr%my_ntau
     !it_glob = gwr%myit_to_glob(my_it)

     ! 1) FFT Transform the first index (local): G(g,g',it) --> G(r,g',it)
     ! 2) MPI transposition: G(r,g',it) --> G^*(g',r,it)
     do my_ikbz=1,gwr%my_nkbz
       !ik_bz = gwr%my_kbz_to_glob(my_ikbz)
       !desc => gwr%green_desc(my_ikbz)
       !call gwr%go(my_ikbz, my_it, my_spin)%fft_gg(desc)
       !call gwr%ge(my_ikbz, my_it, my_spin)%fft_gg(desc)
       call gwr%ggp_to_gpr("go", my_ikbz, my_it, my_spin, desc)
       call gwr%ggp_to_gpr("ge", my_ikbz, my_it, my_spin, desc)
     end do

     !my_nr = gwr%go(my_ikbz, my_it, my_spin)%g_r(1)%sizeb_local(2)
     !my_nr = gwr%go_gr(my_ikbz, my_it, my_spin)%sizeb_local(2)

     ! ============================
     ! GWr algorithm with supercell
     ! ============================
     do ir=1,my_nr
       ! Insert Green's functions in G'-space in the supercell FFT box.
       ! Note that we need to take the union of (k, g') for k in the BZ.
       ! This requires MPI communication inside the BZ communicator.

       do my_ikbz=1,gwr%my_nkbz
         !do ii=1,npw_k
         !    sc_gvec(:,ii) = nint(kpt_bz * ngkpt) + ngkpt * gvec(:,ii)
         !end do
         !gwr%go_gr(my_ikbz, my_it, my_spin)%buffer_cplx(:, ir)
         !gwr%ge_gr(my_ikbz, my_it, my_spin)%buffer_cplx(:, ir)
         !call sphere(cg,ndat,npw,cfft_o,n1,n2,n3,n4,n5,n6,kg_k,istwf_k,iflag,me_g0,shiftg,symm,xnorm)
         !call sphere(cg,ndat,npw,cfft_e,n1,n2,n3,n4,n5,n6,kg_k,istwf_k,iflag,me_g0,shiftg,symm,xnorm)
       end do

       !call fftbox_execute_ip_dpc(sc_plan, cfft_o)
       !call fftbox_execute_ip_dpc(sc_plan, cfft_e)

       ! Multiply by phase

       ! Now we have chi(R', r, it)
       !sc_chi = cfft_o * cfft_e

       ! Go back to G' = k + g' space immediately and redistribute (k, g') in k-communicator
       ! Also note that the G-sphere for chi is not necessarly equal to the one used for the Green's function.
       !call fftbox_execute_ip_dpc(sc_plan, sc_chi)

       !do my_ikbz=1,gwr%my_nkbz
         !gsuper = kk + gvec
         !call sphere(cg,ndat,npw,cfft_o,n1,n2,n3,n4,n5,n6,kg_k,istwf_k,iflag,me_g0,shiftg,symm,xnorm)
         !chi_gr(my_ikbz, my_spin)%buffer_cplx(:, ir) = ??
       !end do

     end do ! ir

     ! These arrays are not needed anymore. Make some space for chi
     !call gwr%go_gr(my_ikbz, my_it, my_spin)%free()
     !call gwr%ge_gr(my_ikbz, my_it, my_spin)%free()

     !chi_gr(my_ikbz, my_spin)%buffer_cplx(:, ir)
     !call transpose(chi_gr, chi_rg)
     !call chi_gr%free()
     ! chi_q(r, g') --> chi_q(g, g')
     !do ig2=1,chi_rg%sizeb_local(2)
     !   chi_rg%buffer_cplx(:, ig2)
     !   gwr%chi_gg(my_ikbz, my_it, my_spin)%buffer_cplx(:, ig2) = ??
     !end do
     !call chi_rg%free()

     ! ===================================================
     ! Mixed-space algorithm in unit cell and convolutions
     ! ===================================================
     ! I assume no k-point distribution
     ! Each core has all the G_k in the IBZ and we only need to rotate to get the BZ on the fly.

#if 0
     do iq_ibz=1,gwr%nq_ibz
       qq_ibz = ??
       do my_ikbz=1,gwr%my_nkbz
         kk = ??
         kq = kk + qq_ibz
         ik_ibz = ??
         ikq_ibz = ??
         ! Use symmetries to get Go_kq and Ge_k from the images in the IBZ.
         go_kq_ggp = symmetrize_from(ikq_ibz, symms)
         go_kq_gpr = fftg1_and_mpi_transpose(go_kq_ggp)
         call go_kq_ggp%free()
         ge_k_ggp = symmetrize_from(ik_ibz, symms)
         ge_k_gpr = fftg1_and_mpi_transpose(ge_kq_ggp)
         call ge_kq_ggp%free()
         do ir=1,my_nr
            !go_rpr = FFT_gr(go_kq_gpr%buffer_cplx(:, ir))
            !ge_rpr = FFT_gr(ge_k_gpr%buffer_cplx(:, ir))
            !chi_fft = go_rpr * ge_rpr
            ! from r' to g'
            !chi_gr(:, ir) = FFT_rg(chi_fft)
         end do
       end do
       !call go_kq_gpr%free()
       !call ge_k_ggp%free()

       !chi_rg = chi_gr%transpose()
       !chi_ggp(:, :, iq_ibz, my_it, my_spin) = FFT_rg
     end do ! iq_ibz
#endif

   end do ! it

   ! Transform irred chi from tau to omega
   !call gwr%cos_transform("t2w", my_ikbz, my_it, my_spin)

 end do ! spin

end function gwr_new
!!***

!----------------------------------------------------------------------

!!****f* m_gwr_base/gwr_free
!! NAME
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

subroutine gwr_free(self)

!Arguments ------------------------------------
 class(gwr_t), intent(inout) :: self

!Local variables-------------------------------
 integer :: ii, jj, kk

! *************************************************************************

 !call self%crystal%free()

 do kk=1,size(self%go, 3)
   do jj=1,size(self%go, 2)
     do ii=1,size(self%go, 1)
       call self%go(ii, jj, kk)%free()
       call self%ge(ii, jj, kk)%free()
     end do
   end do
 end do

 ABI_SFREE(self%go)
 ABI_SFREE(self%ge)

end subroutine gwr_free
!!***

!----------------------------------------------------------------------

!!****f* m_gwr_base/gwr_ggp_to_gpr
!! NAME
!!
!! FUNCTION
!!
!!   1) FFT Transform the first index (local): G(g,g',it) --> G(r,g',it)
!!   2) MPI transposition: G(r,g',it) --> G^*(g',r,it)
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

subroutine gwr_ggp_to_gpr(gwr, what, my_ikbz, my_it, my_spin, desc)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr
 character(len=*),intent(in) :: what
 integer,intent(in) :: my_ikbz, my_it, my_spin
 class(grdesc_t),intent(in) :: desc

!Local variables-------------------------------
!scalars
 !integer,parameter :: ndat1 = 1
 integer :: iab, ig2
 type(matrix_scalapack), pointer :: g_g, g_r
 type(matrix_scalapack) :: r_g
 !type(processor_scalapack) :: slk_processor

! *************************************************************************

 !desc%ngfft
 !npw = desc%npw
 !nfft = desc%nfft
 !nfftot = desc%nfftot

 !select case (what)
 !case ("go")
 !  g_g => gwr%go_g_g(my_ikbz, my_it, my_spin)
 !  g_r => gwr%go_g_r(my_ikbz, my_it, my_spin)
 !  desc => gwr%green_desc(my_ikbz)
 !case ("ge")
 !  g_g => gwr%ge_g_g(my_ikbz, my_it, my_spin)
 !  g_r => gwr%ge_g_r(my_ikbz, my_it, my_spin)
 !  desc => gwr%green_desc(my_ikbz)
 !case ("chi")
 !  g_g => gwr%chi_g_g(my_ikbz, my_it, my_spin)
 !  g_r => gwr%chi_g_r(my_ikbz, my_it, my_spin)
 !  desc => gwr%chi_desc(my_ikbz)
 !case default
 !  MSG_ERROR(sjoin("Invalid value for what:", what))
 !end select

 !call init_matrix_scalapack(r_g, nfftot * gwr%nspinor, g_g%sizeb_global(2), g_g%processor, 1, tbloc=tbloc)

 do ig2=1,g_g%sizeb_local(2)
   !call fft_ug(npw, nfft, gwr%nspinor, ndat1, mgfft, ngfft, istwf_k, kg_k, gbound_k, &
               !ug => g_g%buffer_cplx(:, ig2),
               !ur => r_g%buffer_cplx(:, ig2)
               !ug, ur)

   ! Multiply by k-dependent phase factor.
   !r_g%buffer_cplx(:, ig2) =
 end do

 ! Transpose and release r_g
 !call g_r%free()
 !call init_matrix_scalapack(g_r, g_g%sizeb_global(2), nfftot, g_g%processor, 1, tbloc=tbloc)

 !call r_g%transpose("C", g_r)
 !call r_g%free()

end subroutine gwr_ggp_to_gpr
!!***

!----------------------------------------------------------------------

!!****f* m_gwr_base/gwr_chiwt_transform
!! NAME
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

subroutine gwr_chiwt_transform(gwr, mode, my_ikbz, my_spin)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr
 character(len=*),intent(in) :: mode
 integer,intent(in) :: my_ikbz, my_spin

!Local variables-------------------------------
!scalars
 !integer,parameter :: ndat1 = 1
 integer :: ig1, ig2, it, iw
 !type(matrix_scalapack), pointer :: g_g !, g_r
 !type(matrix_scalapack) :: r_g
 !type(processor_scalapack) :: slk_processor
!arrays
 real(dp),allocatable :: my_weights(:,:)
 complex(dp),allocatable :: cwork_t(:), cwork_wglb(:)

! *************************************************************************

 !if what == "chi" then
 !  chit(:) => gwr%chi(gwr%my_nqbz, 1:gwr%my_ntau, my_spin))
 !  desc => gwr%chi_desc(my_ikbz)
 !else
 !  MSG_ERROR(sjoin("Invalid value for what:", what))
 !end if

 ABI_MALLOC(my_weights, (gwr%ntau, gwr%my_ntau))
 ABI_MALLOC(cwork_t, (gwr%my_ntau))
 ABI_MALLOC(cwork_wglb, (gwr%my_ntau))

 ! Extract my weights from global array according to mode.
 !do it=1,gwr%my_ntau
 !  it_glob = gwr%myit_to_glob(it)
 !  my_weights(:, it) = cos(omega * tau) * global_weighs(:, it_glob)
 !  select case (mode)
 !  case ('w2t")
 !  case ('t2w")
 !  case default
!     ABI_ERROR(sjoin("Wrong mode", mode))
 !  end select
 !end do

 !do ig2=1,chit(1)%sizeb_local(2)
 !  do ig1=1,chit(1)%sizeb_local(1)
 !    ! Extract matrix elements as a function of tau
 !    ! Here we can block over ig1 and call zgemm
 !    ! in order to reduced the number of MPI communications.
 !    do it=1,gwr%my_ntau
 !      cwork_t(it) = chit(it)%buffer_cplx(ig1, ig2)
 !    end do
 !    do iw=1,gwr%ntau
 !      cwork_w(iw) = dot_product(my_weights(iw, :), cwork_t)
 !    end do

 !    !call xmpi_sum(cwork_w, gwr%tau_comm%value, ierr)

 !    ! update values for this (g1, g2), now we have a slice of gg' in frequency space.
 !    !do it=1,gwr%my_ntau
 !    !  !it_glob = gwr%myit_to_glob(it)
 !    !  chit(it)%buffer_cplx(ig1, ig2) = cwork_t(it_glob)
 !    !end do

 !  end do ! ig1
 !end do ! ig2

 ABI_FREE(cwork_t)
 ABI_FREE(cwork_wglb)
 ABI_FREE(my_weights)

end subroutine gwr_chiwt_transform
!!***

!----------------------------------------------------------------------

!!****f* m_gwr_base/green_free
!! NAME
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

subroutine green_free(self)

!Arguments ------------------------------------
 class(green_t), intent(inout) :: self

!Local arguments
 integer :: iab

! *************************************************************************

 !call self%g_g(iab)%free()
 !call self%g_r(iab)%free()
 !call self%r_r(iab)%free()

 ABI_SFREE(self%g_g)
 ABI_SFREE(self%g_r)
 ABI_SFREE(self%r_r)

end subroutine green_free
!!***

!----------------------------------------------------------------------

!!****f* m_gwr_base/slk_alloc_transpose
!! NAME
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

type(matrix_scalapack) function slk_init_transpose(self) result (out_mat)

!Arguments ------------------------------------
 class(matrix_scalapack), intent(in) :: self
 !class(matrix_scalapack), intent(out) :: out_mat

!Local variables-------------------------------
 !type(processor_scalapack) :: slk_processor

! *************************************************************************

 !call init_matrix_scalapack(out_mat, self%sizeb_global(2), self%sizeb_global(1), slk_processor, 1, tbloc=tbloc)

end function slk_init_transpose
!!***

!----------------------------------------------------------------------

!!****f* m_gwr_base/slk_transpose
!! NAME
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

subroutine slk_transpose(self, trans, out_mat)

!Arguments ------------------------------------
 class(matrix_scalapack), intent(in) :: self
 character(len=1),intent(in) :: trans
 class(matrix_scalapack), intent(inout) :: out_mat

! *************************************************************************

 ! Transposes a complex distributed matrix, conjugated.
 ! sub(C):=beta*sub(C) + alpha*conjg(sub(A)'),
 !if (trans == "C") then
 !call pztranc(self%sizeb_global(2), self%sizeb_global(1), cone, self%buffer_cplx, 1, 1,
 !             self%descript%tab,, czero, out_mat%buffer_cplx, 1, 1, out_mat%descript%tab)

 ! sub(C):=beta*sub(C) + alpha*sub(A)',
  !call pztranu(m, n, alpha, a, ia, ja, desca, beta, c, ic, jc, descc)

  !call pdtran(m, n, alpha, a, ia, ja, desca, beta, c, ic, jc, descc)

end subroutine slk_transpose
!!***

end module m_gwr_base
!!***
