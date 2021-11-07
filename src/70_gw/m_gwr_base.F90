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

   integer :: nfft = -1
   ! Number of fft points

   integer :: mgfft = -1
   ! maximum size of 1D FFTs

   integer :: ngfft(18) = -1
   ! FFT mesh used to transform G(g,g')

   integer,allocatable :: gbound(:,:)
   !  gbound(2*mgfft+8,2)
   !  sphere boundary info for reciprocal to real space (local)

   integer,allocatable :: gvec(:,:)
   ! gvec(3, gvec_size)
   ! For each k, the g-vectors are ordered according to |k+g|
   ! so that we can use the same array for Green, chi and Sigma_x
   ! Note that this array is global i.e. it is not MPI-distributed in G-space, only in k-space.

   !integer,allocatable :: kg_k(:,:)
   ! kg_k(3,npw)
   ! G vector coordinates in reduced cordinates.

   !% real(dp) :: kpoint(3)

 contains

   !procedure :: set_ngfft => gr_desc_set_ngfft

   procedure :: free => gr_desc_free
   ! Free memory.

 end type grdesc_t
!!***

!type, extends(matrix_scalapack) :: twopoint_func_t
!   type(grdesc_t), pointer :: desc => null()
!
!   type(matrix_scalapack) :: g_g
!   type(matrix_scalapack) :: r_r
!   type(matrix_scalapack) :: g_r
!
!end type twopoint_func_t

!----------------------------------------------------------------------

!!****t* m_gwr_base/gwr_t
!! NAME
!! gwr_t
!!
!! FUNCTION
!!
!! SOURCE

 type, public :: gwr_t

   integer :: nspinor

   integer :: nsppol
   ! Global number of independent spin polarizations.

   integer :: my_nsppol
   ! Number of independent spin polarizations treated by this MPI proc

   integer :: nkbz
   ! Number of k-points in the full BZ

   integer :: nkibz
   ! Number of k-points in the IBZ

   integer :: nk
   ! nkbz if GWr, nkibz if GWrk

   integer :: nq

   integer :: nqbz
   ! Number of q-points in the full BZ

   integer :: nqibz
   ! Number of q-points in the IBZ

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

   integer :: gf_mpw = -1
   ! Max number of g-vectors for Green's function over k-points

   integer :: chi_mpw = -1
   ! Max number of g-vectors for Chi over q-points

   type(grdesc_t),allocatable :: gf_desc_k(:)
   ! gf_desc_k(nk)

   type(grdesc_t),allocatable :: chi_desc_q(:)
   ! chi_desc_q(nq)

   !integer,allocatable :: gvec(:,:)
   ! gvec(3, gvec_size, nk)
   ! For each k, the gvectors are ordered according to |k+g|
   ! so that we can use the same array for G, chi and Sigma_x
   ! Note that this array is global i.e. is not MPI-distributed in G-space, only in k-space.

   !integer, allocatable :: npw_green(:)
   !integer, allocatable :: npw_chi(:)
   !integer, allocatable :: npw_sigma(:)
   ! Arrays of shape (nk) giving the number of G-vectors in gvec.
   ! These quantities are computed from the values of ecut, ecuteps and ecutsigx
   ! specified in input.

   !integer :: coords_pqbks(5) = 0
   ! Cartesian coordinates of this processor in the Cartesian grid.

   integer :: comm
   ! MPI communicator with all procs involved in the calculation

   type(xcomm_t) :: spin_comm
   ! MPI communicator for parallelism over spins (high-level)

   !type(xcomm_t) :: k_comm
   ! MPI communicator for k-point distribution

   type(xcomm_t) :: gr_comm
   ! MPI communicator for gr distribution

   type(xcomm_t) :: tau_comm
   ! MPI communicator for imag time distribution

   !type(xcomm_t) :: omega_comm
   ! MPI communicator for imag frequency distribution

   type(dataset_type), pointer :: dtset => null()

   type(crystal_t) :: cryst

   type(matrix_scalapack),allocatable :: go_ggp(:,:,:)
   !type(matrix_scalapack),allocatable :: go_gr(:,:,:)
   !type(matrix_scalapack),allocatable :: go_rrp(:,:,:)
   type(matrix_scalapack),allocatable :: ge_ggp(:,:,:)
   !type(matrix_scalapack),allocatable :: ge_gr(:,:,:)
   !type(matrix_scalapack),allocatable :: ge_rrp(:,:,:)
   ! (nk, my_ntau, my_nsppol)

   type(matrix_scalapack),allocatable :: chi_ggp(:,:,:)
   ! (nq, my_ntau, my_nsppol)

   !character(len=fnlen) :: wfk_filepath
   ! Path to the WFK file with the KS wavefunctions.

  real(dp),allocatable :: kbz(:,:)
  ! kbz(3, nkbz)
  ! Reduced coordinates of the k-points in the full BZ.

  real(dp),allocatable :: kibz(:,:)
   ! kibz(3, nkibz)
   ! Reduced coordinates of the k-points in the IBZ (full simmetry of the system).

  real(dp),allocatable :: wt_kibz(:)
   ! wt_kibz(nqibz)
   ! Weights of the k-points in the IBZ (normalized to one).

  real(dp),allocatable :: qbz(:,:)
  ! qbz(3, nqbz)
  ! Reduced coordinates of the q-points in the full BZ.

  real(dp),allocatable :: qibz(:,:)
   ! qibz(3, nqibz)
   ! Reduced coordinates of the q-points in the IBZ (full simmetry of the system).

  real(dp),allocatable :: wt_qibz(:)
   ! wtq_qibz(nqibz)
   ! Weights of the q-points in the IBZ (normalized to one).

   !integer,allocatable:: indkk_kq(:, :)
   ! (6, %nqibz_k))
   ! Mapping k+q --> initial IBZ. Depends on ikcalc.
   ! These table used the conventions for the symmetrization of the wavefunctions expected by cgtk_rotate.
   ! In this case listkk has been called with symrel and use_symrec=False

 contains

   !procedure :: new => gwr_new

   procedure :: ggp_to_gpr  => gwr_ggp_to_gpr
   !procedure :: cos_transform_q  => gwr_cos_transform_q
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
 ABI_SFREE(self%gbound)

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

function gwr_new(dtset, comm) result (gwr)

!Arguments ------------------------------------
!scalars
 type(dataset_type),target,intent(in) :: dtset
 integer,intent(in) :: comm
 type(gwr_t),target :: gwr
!arrays

!Local variables-------------------------------
!scalars
 integer, parameter :: istwf_1 = 1
 integer :: my_spin, my_it, ik, iq, ir, ii, my_nr, npwsp, col_bsize
 integer :: iq_ibz, mgfft, nfft
 real(dp) :: ecuteff
 type(grdesc_t),pointer :: desc_k, desc_kq, desc_q
 type(processor_scalapack) :: slk_processor
!arrays
 integer,allocatable :: gfsc_gvec(:,:), chisc_gvec(:,:)
 integer :: tmp_ngfft(18), ngfft(18)
 real(dp) :: kk(3), qq(3), kq(3)
 type(matrix_scalapack),allocatable :: go_gpr(:), ge_gpr(:)

! *************************************************************************

 gwr%dtset => dtset
 gwr%comm = comm
 gwr%cryst = dtset%get_crystal(img=1)

 gwr%nspinor = dtset%nspinor
 gwr%nsppol = dtset%nsppol
 gwr%use_supercell = .True.
 gwr%my_nsppol = 1

 ! =======================
 ! Setup k-mesh and q-mesh
 ! =======================
 ! TODO: Should rearrange kbz so that the first nkibz items are equal to kbz
 ! and we can simply loop like
 ! do ik=1,gwr%nk
 !   kk = gwr%kbz(ik)
 ! end do
 ! use similar approach for the q-points.

 gwr%nkbz = 1
 gwr%nkibz = 1
 gwr%nqbz = 1
 gwr%nqibz = 1

 ABI_MALLOC(gwr%kbz, (3, gwr%nkbz))
 ABI_MALLOC(gwr%kibz, (3, gwr%nkibz))
 ABI_MALLOC(gwr%wt_kibz, (gwr%nkibz))

 ABI_MALLOC(gwr%qbz, (3, gwr%nqbz))
 ABI_MALLOC(gwr%qibz, (3, gwr%nqibz))
 ABI_MALLOC(gwr%wt_qibz, (gwr%nqibz))

 if (gwr%use_supercell) then
   gwr%nk = gwr%nkbz
   gwr%nq = gwr%nqbz
 else
   gwr%nk = gwr%nkibz
   gwr%nq = gwr%nqibz
 end if

 ! ====================
 ! Setup tau/omega mesh
 ! ====================
 gwr%ntau = 1
 gwr%my_ntau = 1
 gwr%niw = 1
 gwr%my_niw = 1

 ! TODO: Set MPI grid and communicators.
 call gwr%spin_comm%set_to_self()
 call gwr%gr_comm%set_to_self()
 call gwr%tau_comm%set_to_self()

 ! Build FFT descriptors for Green's functions
 ABI_MALLOC(gwr%gf_desc_k, (gwr%nk))
 do ik=1,gwr%nk
   kk = gwr%kbz(:, ik)
   desc_k => gwr%gf_desc_k(ik)
   !gwr%gf_desc_k(ik) = gr_desc_new(ecuteff, kpoint, gwr%cryst%gmet, ngfft=green_ngfft)
   ecuteff = dtset%ecut * dtset%dilatmx**2
   call get_kg(kk, istwf_1, ecuteff, gwr%cryst%gmet, desc_k%npw, desc_k%gvec)
   !call getng(boxcutmin,ecut,gmet,kpt,me_fft,mgfft,nfft,ngfft,nproc_fft,nsym,paral_fft,symrel,&
   !           ngfftc,use_gpu_cuda,unit) ! optional
   gwr%gf_mpw = max(gwr%gf_mpw, gwr%gf_desc_k(ik)%npw)
 end do

 ! Build FFT descriptors for chi
 ABI_MALLOC(gwr%chi_desc_q, (gwr%nq))
 do iq=1,gwr%nq
   qq = gwr%qbz(:, iq)
   desc_q => gwr%chi_desc_q(iq)
   !gwr%chi_desc_q(iq) = gr_desc_new(ecuteff, qpoint, gwr%cryst%gmet, ngfft=chi_ngfft)
   call get_kg(qq, istwf_1, dtset%ecuteps, gwr%cryst%gmet, desc_q%npw, desc_q%gvec)
   gwr%chi_mpw = max(gwr%chi_mpw, gwr%chi_desc_q(iq)%npw)
 end do

 ! Now we know the value of ngfft. Setup tables for zero-padded FFTs.
 !ngfft = ??
 nfft = product(ngfft(1:3)); mgfft = maxval(ngfft(1:3))

 do ik=1,gwr%nk
   desc_k => gwr%gf_desc_k(ik)
   desc_k%nfft = nfft; desc_k%mgfft = mgfft; desc_k%ngfft = ngfft
   ABI_MALLOC(desc_k%gbound, (2 * mgfft + 8, 2))
   call sphereboundary(desc_k%gbound, istwf_1, desc_k%gvec, mgfft, desc_k%npw)
 end do

 do iq=1,gwr%nq
   desc_q => gwr%chi_desc_q(ik)
   desc_q%nfft = nfft; desc_q%mgfft = mgfft; desc_q%ngfft = ngfft
   ABI_MALLOC(desc_q%gbound, (2 * mgfft + 8, 2))
   call sphereboundary(desc_q%gbound, istwf_1, desc_q%gvec, mgfft, desc_q%npw)
 end do

 ! G_k and chi_q
 ABI_MALLOC(gwr%go_ggp, (gwr%nk, gwr%my_ntau, gwr%my_nsppol))
 ABI_MALLOC(gwr%ge_ggp, (gwr%nk, gwr%my_ntau, gwr%my_nsppol))
 ABI_MALLOC(gwr%chi_ggp, (gwr%nq, gwr%my_ntau, gwr%my_nsppol))

 call init_scalapack(slk_processor, gwr%gr_comm%value) !# grid_shape=[1, gwr%gr_comm%value]

 ABI_MALLOC(gfsc_gvec, (3, gwr%gf_mpw))
 ABI_MALLOC(chisc_gvec, (3, gwr%chi_mpw))
 ABI_FREE(gfsc_gvec)
 ABI_FREE(chisc_gvec)

 do my_spin=1,gwr%my_nsppol
   do my_it=1,gwr%my_ntau

     ! Allocate G_k(g, g')
     do ik=1,gwr%nk
       npwsp = gwr%gf_desc_k(ik)%npw * gwr%nspinor
       col_bsize = npwsp / gwr%gr_comm%nproc
       if (mod(npwsp, gwr%gr_comm%nproc) /= 0) col_bsize = col_bsize + 1
       !col_bsize = 50
       call init_matrix_scalapack(gwr%go_ggp(ik, my_it, my_spin), &
          npwsp, npwsp, slk_processor, 1) ! , tbloc=[npwsp, col_bsize])
       call init_matrix_scalapack(gwr%ge_ggp(ik, my_it, my_spin), &
          npwsp, npwsp, slk_processor, 1) !, tbloc=[npwsp, col_bsize])
     end do

     ! Allocate chi_k(g, g')
    do iq=1,gwr%nq
      npwsp = gwr%chi_desc_q(iq)%npw * gwr%nspinor
      col_bsize = npwsp / gwr%gr_comm%nproc
      if (mod(npwsp, gwr%gr_comm%nproc) /= 0) col_bsize = col_bsize + 1
      !col_bsize = 50
      call init_matrix_scalapack(gwr%chi_ggp(iq, my_it, my_spin), &
         npwsp, npwsp, slk_processor, 1) ! , tbloc=[npwsp, col_bsize])
      call init_matrix_scalapack(gwr%chi_ggp(iq, my_it, my_spin), &
         npwsp, npwsp, slk_processor, 1) !, tbloc=[npwsp, col_bsize])
    end do

   end do
 end do

 ! Prepare FFT in the supercell from R' -> G'
 ! Be careful when using the FFT plan with ndat as ndata can change inside the loop.
 ! Perhaps the safest approach would be to generate the plan on the fly.
 !integer :: sc_ngfft(18)
 !sc_ngfft(1:3) = gwr%ngkpt * ngfft(1:3)
 !sc_nr = product(sc_ngfft(1:3))
 !call fftbox_plan3(sc_plan, sc_ngfft, fftalg, isign)
 !cal fftbox_plan3_init(plan,ndat,dims,embed,fftalg,fftcache,isign)
 !ABI_MALLOC(sc_chi, (sc_nr))
 !plan%ldxyz*plan%ndat
 !ABI_FREE(sc_chi)

 ABI_MALLOC(go_gpr, (gwr%nk))
 ABI_MALLOC(ge_gpr, (gwr%nk))

 do my_spin=1,gwr%my_nsppol
   !spin = gwr%my_spin_to_glob(my_spin)
   do my_it=1,gwr%my_ntau
     !it_glob = gwr%myit_to_glob(my_it)

     do ik=1,gwr%nk
       call gwr%ggp_to_gpr("go", ik, my_it, my_spin, go_gpr(ik))
       call gwr%ggp_to_gpr("ge", ik, my_it, my_spin, ge_gpr(ik))
     end do

     my_nr = go_gpr(1)%sizeb_local(2)

     ! ============================
     ! GWr algorithm with supercell
     ! ============================
     do ir=1,my_nr
       ! Insert Green's functions in G'-space in the supercell FFT box.
       ! Note that we need to take the union of (k, g') for k in the BZ.

       do ik=1,gwr%nk
         desc_k => gwr%gf_desc_k(ik)
         do ii=1,desc_k%npw
           !gfsc_gvec(:,ii) = nint(gwr%kbz(ik) * gwr%ngkpt) + gwr%ngkpt * desc_k%gvec(:,ii)
         end do

         !go_gpr(ik)%buffer_cplx(:, ir)
         !call sphere(cg,ndat,npw,cgo,n1,n2,n3,n4,n5,n6,kg_k,istwf_k,iflag,me_g0,shiftg,symm,xnorm)

         !ge_gpr(ik)%buffer_cplx(:, ir)
         !call sphere(cg,ndat,npw,cge,n1,n2,n3,n4,n5,n6,kg_k,istwf_k,iflag,me_g0,shiftg,symm,xnorm)
       end do

       !call fftbox_execute_ip_dpc(sc_plan, cgo)
       !call fftbox_execute_ip_dpc(sc_plan, cge)

       ! Multiply by the phase if needed.

       ! Now we have chi(R', r, it)
       !sc_chi = cgo * cge

       ! Go back to G' = k + g' space immediately and redistribute (k, g') in k-communicator
       ! Also note that the G-sphere for chi is not necessarly equal to the one used for the Green's function.
       !call fftbox_execute_ip_dpc(sc_plan, sc_chi)

       ! Note that here are using the g-sphere for chi.
       do iq=1,gwr%nq
         desc_q => gwr%chi_desc_q(iq)
         do ii=1,desc_q%npw
           !chisc_gvec(:,ii) = nint(gwr%qbz(iq) * gwr%ngkpt) + gwr%ngkpt * desc_q%gvec(:,ii)
         end do

         !call sphere(cg,ndat,npw,cfft_o,n1,n2,n3,n4,n5,n6,kg_k,istwf_k,iflag,me_g0,shiftg,symm,xnorm)
         !chi_gpr(iq, my_spin)%buffer_cplx(:, ir) = ??
       end do

     end do ! ir

     do iq=1,gwr%nq
       !chi_gpr(iq, my_spin)%buffer_cplx(:, ir)
       !if (iq == 1) then
       !  chi_rgp = chi_gpr(iq, my_spin)%alloc_transpose()
       !end if
       !call chi_gpr(iq, my_spin)%ptran("N", chi_rgp)
       !call chi_gpr(iq, my_spin)%free()
       ! FFT along the first dimesions: chi_q(r, g') --> chi_q(g, g')
       !do ig2=1,chi_rgp%sizeb_local(2)
       !   chi_rgp%buffer_cplx(:, ig2)
       !   gwr%chi_gg(iq, my_it, my_spin)%buffer_cplx(:, ig2) = ??
       !end do
     end do
     !call chi_rgp%free()

     ! ===================================================
     ! Mixed-space algorithm in unit cell and convolutions
     ! ===================================================
     ! I assume no k-point distribution
     ! Each core has all the G_k in the IBZ and we only need to rotate to get the BZ on the fly.

     do iq_ibz=1,gwr%nqibz
       !qq_ibz = gwr%qbz(:, iq_ibz)
       !desc_q => gwr%chi_desc_q(iq_ibz)
       do ik=1,gwr%nk
         ! Use symmetries to get Go_kq and Ge_k from the images in the IBZ.
         !kk = gwr%kbz(:, ik)
         !kq = kk + qq_ibz
         !ik_ibz = ??
         !ikq_ibz = ??
         !desc_k => gwr%gf_desc_k(ik)
         !desc_kq => gwr%gf_desc_k(ik)
         !gf_ibz => gwr%go_ggp(ikq_ibz, my_it, my_spin)
         !go_kq_ggp = symmetrize_from(ikq_ibz, symms)
         !go_kq_gpr = fftg1_and_mpi_transpose(go_kq_ggp)
         !call go_kq_ggp%free()
         !gf_ibz => gwr%go_ggp(ik_ibz, my_it, my_spin)
         !ge_k_ggp = symmetrize_from(ik_ibz, symms)
         !ge_k_gpr = fftg1_and_mpi_transpose(ge_kq_ggp)
         !call ge_kq_ggp%free()
         !my_nr = ?
         do ir=1,my_nr
            !go_rpr = FFT_gr(go_kq_gpr%buffer_cplx(:, ir))
            !ge_rpr = FFT_gr(ge_k_gpr%buffer_cplx(:, ir))
            !chi_fft = go_rpr * ge_rpr
            ! from r' to g'
            !chi_gpr(:, ir) = FFT_rg(chi_fft)
         end do
         !call go_kq_gpr%free()
         !call ge_k_ggp%free()

         !chi_rgp = chi_gpr%transpose()
         !gwr%chi_ggp(:, :, iq_ibz, my_it, my_spin) = FFT_rg
       end do ! ik
     end do ! iq_ibz

   end do ! it
 end do ! spin

 do ik=1,gwr%nk
   call go_gpr(ik)%free()
   call ge_gpr(ik)%free()
 end do

 ABI_FREE(go_gpr)
 ABI_FREE(ge_gpr)

 ! Transform irred chi from tau to omega
 ! and sum the collinear components to get total chi
 do iq=1,gwr%nq
   do my_spin=1,gwr%my_nsppol
     !call gwr%cos_transform_q("chi", "t2w", iq, my_spin)
   end do
   !if (gwr%nsppol == 2) then
   !  call xmpi_sum(cwork_wglb, gwr%spin_comm%value, ierr)
   !end if
 end do

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

 call self%cryst%free()

 ABI_SFREE(self%kbz)
 ABI_SFREE(self%kibz)
 ABI_SFREE(self%wt_kibz)

 ABI_SFREE(self%qbz)
 ABI_SFREE(self%qibz)
 ABI_SFREE(self%wt_qibz)

 do kk=1,self%nk
   call self%gf_desc_k(kk)%free()
 end do
 ABI_SFREE(self%gf_desc_k)

 do kk=1,self%nq
   call self%chi_desc_q(kk)%free()
 end do
 ABI_SFREE(self%chi_desc_q)

 do kk=1,size(self%go_ggp, 3)
   do jj=1,size(self%go_ggp, 2)
     do ii=1,size(self%go_ggp, 1)
       call self%go_ggp(ii, jj, kk)%free()
       call self%ge_ggp(ii, jj, kk)%free()
     end do
   end do
 end do
 ABI_SFREE(self%go_ggp)
 ABI_SFREE(self%ge_ggp)

 do kk=1,size(self%chi_ggp, 3)
   do jj=1,size(self%chi_ggp, 2)
     do ii=1,size(self%chi_ggp, 1)
       call self%chi_ggp(ii, jj, kk)%free()
     end do
   end do
 end do
 ABI_SFREE(self%chi_ggp)

 ! Free MPI communicators
 call self%spin_comm%free()
 call self%gr_comm%free()
 call self%tau_comm%free()

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

subroutine gwr_ggp_to_gpr(gwr, what, ik, my_it, my_spin, gpr)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr
 character(len=*),intent(in) :: what
 integer,intent(in) :: ik, my_it, my_spin
 type(matrix_scalapack),intent(out) :: gpr

!Local variables-------------------------------
!scalars
 integer,parameter :: ndat1 = 1, istwf_1 = 1
 integer :: iab, ig2
 type(matrix_scalapack), pointer :: ggp
 type(matrix_scalapack) :: rgp, g_r
 type(grdesc_t),pointer :: desc

! *************************************************************************

 select case (what)
 case ("go")
   ggp => gwr%go_ggp(ik, my_it, my_spin)
   desc => gwr%gf_desc_k(ik)
 case ("ge")
   ggp => gwr%ge_ggp(ik, my_it, my_spin)
   desc => gwr%gf_desc_k(ik)
 case ("chi")
   ggp => gwr%chi_ggp(ik, my_it, my_spin)
   desc => gwr%chi_desc_q(ik)
 case default
   ABI_ERROR(sjoin("Invalid value for what:", what))
 end select

 ! X(r, g')
 call init_matrix_scalapack(rgp, desc%nfft * gwr%nspinor, ggp%sizeb_global(2), ggp%processor, istwf_1) !, tbloc=tbloc)

 do ig2=1,ggp%sizeb_local(2)
   !call fft_ug(desc%npw, desc%nfft, gwr%nspinor, ndat1,
   !            desc%mgfft, desc%ngfft, istwf_1, desc%kg_k, desc%gbound_k, &
   !            ggp%buffer_cplx(:, ig2), rgp%buffer_cplx(:, ig2)

   ! Multiply by k-dependent phase factor.
   !rgp%buffer_cplx(:, ig2) =
 end do

 ! Transpose: X(r, g') -> Y(g', r)
 call gpr%free()
 call init_matrix_scalapack(gpr, ggp%sizeb_global(2), desc%nfft * gwr%nspinor, ggp%processor, istwf_1) !, tbloc=tbloc)

 !call rgp%ptran("C", gpr)
 call rgp%free()

end subroutine gwr_ggp_to_gpr
!!***

!----------------------------------------------------------------------

!!****f* m_gwr_base/gwr_cos_transform_q
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

subroutine gwr_cos_transform_q(gwr, what, mode, iq, my_spin)

!Arguments ------------------------------------
 class(gwr_t),target,intent(in) :: gwr
 character(len=*),intent(in) :: what, mode
 integer,intent(in) :: iq, my_spin

!Local variables-------------------------------
!scalars
 integer :: ig1, ig2, it, iw, ierr
 type(grdesc_t),pointer :: desc
 type(matrix_scalapack), pointer :: chit(:)
!arrays
 real(dp),allocatable :: my_weights(:,:)
 complex(dp),allocatable :: cwork_t(:), cwork_wglb(:)

! *************************************************************************

 select case (what)
 case ("chi")
   chit => gwr%chi_ggp(iq, 1:gwr%my_ntau, my_spin)
   desc => gwr%chi_desc_q(iq)
 case default
   ABI_ERROR(sjoin("Invalid value for what:", what))
 end select

 ABI_MALLOC(my_weights, (gwr%ntau, gwr%my_ntau))
 ABI_MALLOC(cwork_t, (gwr%my_ntau))
 ABI_MALLOC(cwork_wglb, (gwr%my_ntau))

 ! Extract my weights from global array according to mode.
 do it=1,gwr%my_ntau
   !it_glob = gwr%myit_to_glob(it)
   !my_weights(:, it) = cos(omega * tau) * global_weighs(:, it_glob)
   select case (mode)
   case ("w2t")
   case ("t2w")
   case default
     ABI_ERROR(sjoin("Wrong mode", mode))
   end select
 end do

 do ig2=1,chit(1)%sizeb_local(2)
   do ig1=1,chit(1)%sizeb_local(1)
     ! Extract matrix elements as a function of tau
     ! Here we can block over ig1 and call zgemm
     ! in order to reduced the number of MPI communications.
     do it=1,gwr%my_ntau
       cwork_t(it) = chit(it)%buffer_cplx(ig1, ig2)
     end do
     do iw=1,gwr%ntau
       cwork_wglb(iw) = dot_product(my_weights(iw, :), cwork_t)
     end do

     call xmpi_sum(cwork_wglb, gwr%tau_comm%value, ierr)

     ! update values for this (g1, g2), now we have a slice of gg' in frequency space.
     do it=1,gwr%my_ntau
       !it_glob = gwr%myit_to_glob(it)
       !chit(it)%buffer_cplx(ig1, ig2) = cwork_wglb(it_glob)
     end do

   end do ! ig1
 end do ! ig2

 ABI_FREE(cwork_t)
 ABI_FREE(cwork_wglb)
 ABI_FREE(my_weights)

end subroutine gwr_cos_transform_q
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

!!****f* m_gwr_base/slk_ptran
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

subroutine slk_ptran(self, trans, out_mat)

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

end subroutine slk_ptran
!!***

end module m_gwr_base
!!***
