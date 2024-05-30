!!****m* ABINIT/m_varpeq
!! NAME
!!  m_varpeq
!!
!! FUNCTION
!!  Variational Polaron Equations
!!  Datatype with methods providing the variatonal treatment of the adiabatic
!! polaron formation, i.e. electron/hole self-localization. At the current
!! state of implementation, it requires preliminary comptation of gstore, a
!! datatype containing MPI-distributed electron-phobon matrix elements and
!! other useful quantities related to the electron-phonon calculations.
!! Ssee, for example, src/78_eph/m_gstore.F90 for the details of implementation
!! and src/78_eph/m_berry_curvature.F90 for the usage in different context.
!!
!!  The parallelism in this module is strogly related to the distribition of
!! gstore components in the memory. In general, they can be distributed over
!! several dimensions: k-points, phonon wavectors q, perturbations, number of
!! electronic bands.
!!  Currently ONLY the parallelism over q-points is assumed and hard-coded.
!!
!!  Additionally, one can not fully rely on the symmetries provided by the
!! crystal's point group in order to reduce the number of k/q-points required
!! for the BZ integration. Indeed, the polaron formation may break the crystal
!! symmetry, and in the variational approach one can safely assume only the
!! inversion symmetry, i.e.
!!
!!     A_nk = A_n{-k} for the electronic variational coefficients,
!!     B_q\nu = B_-q\nu for the phonon vibrational coefficients,
!!
!! provided the crystal itself posseses spatial inversion.
!!
!!  Moreover, felectron-phonon matrix elements can be reduced by symmetry only
!! in k or q-dimension, but NOT in both, as for any symmmetry operation S
!!
!!     g_mn\nu (Sk, q) = g_mn\nu (k, S^{-1}q).
!!
!!  Hence, one has to be cautions when it comes to symmetries used to simplify
!! the problem. In this implementation, the following is assumed:
!!     * groundstate properties (electronic energies, phonons) keep their
!!       initial symmetries and are defined in the irreducible BZ;
!!     * electronic coeffeiceints posses the inversion (time-reversal) symmetry
!!       A_nk = A_n{-k};
!!     * vibrational coefficients are defined in the full BZ in q-space, which
!!       is disribuied over several MPI processes by a gstore object.
!!
!! COPYRIGHT
!!  Copyright (C) 2023-2024 ABINIT group (VV)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_varpeq

 use defs_basis
 use m_abicore
 use m_dtset
 use m_crystal
 use m_ebands
 use m_errors
 use m_ifc
 use m_krank
 use m_xmpi

 use defs_datatypes,    only : ebands_t
 use m_fstrings,        only : sjoin, itoa, ktoa
 use m_gstore,          only : gstore_t, gqk_t
 use m_kpts,            only : kpts_ibz_from_kptrlatt, kpts_map, kpts_timrev_from_kptopt
 use m_symkpt,          only : symkpt
 use m_symtk,           only : mati3inv

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_varpeq/varpeq_t
!! NAME
!!  varpeq_t
!!
!! FUNCTION
!!  Description
!!
!! SOURCE

 type, public :: varpeq_t


 ! Electroncic subspace -----------------------

  integer :: nsppol = -1
   ! Spin polarization:
   ! 1 -> unpolarized; 2 -> polarized calculations

  integer :: nb = -1
   ! Number of bands participating in the polaron formation

  integer :: nkpt = -1
   ! Total number of k-points defining the electronic coeffeicients A_nk

  integer :: nkpt_g = -1

  real(dp), allocatable :: kpts(:, :)
   ! k-points defining the electronic coefficients A_nk, possibly reduced by
   ! inversion symmetry: A_nk = A_n{-k} provided the crystal's space group has
   ! spatial inversion. No symmetries are assumed otherwise.
   ! (3, nkpt)

  real(dp), allocatable :: kpts_g(:, :)

  real(dp), allocatable :: wtk(:)
   ! Weights of the k-points (normalized to one)
   ! (nkpt)

  complex(dp), allocatable :: a(:, :, :)
   ! Variational coefficients in the electornic subspace defining the charge
   ! localization process
   ! (nbands, nkpt, nsppol)

  complex(dp), allocatable :: grad_a(:, :, :)
   ! Gradient of the polaron energy functionl wrt electronic coefficients
   ! (nbands, nkpt, nsppol)


 ! Vibrational (phonon) subspace --------------

  integer :: nqpt = -1
   ! Total number of q-points defining the vibrational coeffeicients B_q\nu

  real(dp), allocatable :: qpts(:, :)
   ! q-points defining the vibrational coefficients B_q\nu
   ! (3, nqpt)

  complex(dp), allocatable :: b(:, :)
   ! Variational coeffeicients in the vibrational subspace defining the
   ! deformation potential; treated by this MPI rank
   ! (my_npert, my_nq)

  complex(dp), allocatable :: bq_glob(:, :)

 ! Electron-phonon subspace -------------------

  class(gstore_t), pointer :: gstore => null()
   ! gstore datatype storing MPI-distributed electron-phonon matrix elements
   ! and other related quantities


 ! Symmetry tables ----------------------------

  integer :: nkbz = -1
   ! Size of the BvK supercell (total number of k-points in full BZ)

  integer :: nsym = -1
   ! Number of symmetries defining effective IBZ

  integer :: timrev = -1
   ! Time-reversal symmetry

  type(crystal_t) :: cryst
   ! structure defining geometry w/o symmetries except identity and (possibly)
   ! inversion

  type(krank_t) :: krank_kpts
   ! krank datatype used to obtain k+q->k' & k-q->k'' mappings

  type(krank_t) :: krank_kpts_g

  type(krank_t) :: krank_qpts
   ! krank datatype used to obtain q-indices

  integer, allocatable :: map_k_kibz(:, :)
  integer, allocatable :: map_invk_kibz(:, :)
   ! mappings: k-points from the electronic subspace of the problem -> k-points
   ! from IBZ; second array maps an inverse set of k-points
   ! (6, nkpt)

  integer, allocatable :: map_q_qibz(:, :)

  integer :: nkpt_half
  real(dp),allocatable :: kpts_half(:,:)
  real(dp),allocatable :: wtk_half(:)
  integer,allocatable :: map_k_khalf(:,:)
  type(krank_t) :: krank_kpts_half

  integer :: nqpt_half
  real(dp),allocatable :: qpts_half(:,:)
  real(dp),allocatable :: wtq_half(:)
  integer,allocatable :: map_q_qhalf(:,:)
  type(krank_t) :: krank_qpts_half

 ! Optimization variables ---------------------

  real(dp) :: theta
   ! Line minimization parameter


 ! Polaronic properties -----------------------

  real(dp) :: eps
   ! Lagrange multiplier of the optimization problem, equal to the energy
   ! of a localized polaronic state

  real(dp) :: enterms(3)
   ! Polaron formation energy terms
   ! enpol(1) -> electronic, enpol(2) -> vibrational, enpol(3) -> el-phonon

  contains

    procedure :: init => varpeq_init
     ! Initialiation of basic dimensions, parameters and symmetry relations

    procedure :: free => varpeq_free
     ! Free allocated memory

    procedure :: get_mintheta => varpeq_get_mintheta
     ! Calculatae line minimization parameter theta

    procedure :: orth_elgrad => varpeq_orth_elgrad
     ! Orthogonalization of the electronic gradient wrt to the electronic
     ! vector A_nk

    procedure :: get_elgrad => varpeq_get_elgrad
     ! Calculate the gradient of polaron energy functional wrt electronic
     ! subspace coefficients

    procedure :: get_enterms => varpeq_get_enterms
     ! Calculate polaron formation energy terms

    procedure :: b_from_a => varpeq_b_from_a
     ! Calculate deformation potential, adjusted to the charge localization
     ! denstiy

    procedure :: seed_a => varpeq_seed_a
     ! Seed initial charge localization (vector of electronic coefficients)

    procedure :: test => varpeq_test
     ! Test subroutine used solely for the debug purpose

    procedure :: get_invsym_k => varpeq_get_invsym_k
    procedure :: get_invsym_khalf => varpeq_get_invsym_khalf
    procedure :: get_invsym_q => varpeq_get_invsym_q
    procedure :: get_invsym_qhalf => varpeq_get_invsym_qhalf

    procedure :: get_mapping => varpeq_get_mapping

    procedure :: get_gavg => varpeq_get_gavg

 end type varpeq_t
!!***


contains !=====================================================================
!!***


!!****f* m_varpeq/varpeq_test
!! NAME
!!  varpeq_test
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_test(self)

!Arguments ------------------------------------
!scalars
 class(varpeq_t), intent(inout) :: self

!Local variables-------------------------------
!scalars
 integer :: comm, nproc, my_rank, ierr
 integer :: ik, ik_ibz, ik_ibz_g, tsign_k
 integer :: trev_k, isym_k
 integer :: my_iq, iq_glob, my_is
 integer :: ikq, ikq_ibz, iq_ibz
 integer :: iSq_glob, iSk, imq_glob
 integer :: my_pert
 integer :: ii, jj, pp
 integer :: is, ib, jb
 integer :: bstart
 integer :: ib_ndeg, jb_ndeg, pert_ndeg
 integer :: iq_ibz_g, tsign_q
 integer :: ikq_ibz_g, tsign_kq
 integer :: isk_ibz_g, tsign_sk
 integer :: isq_ibz_g, tsign_sq
 integer :: iskq_ibz_g, tsign_skq
 integer :: iskq
 real(dp) :: eig
 real(dp) :: gradres
 real(dp) :: costh, sinth
 real(dp) :: bden
 real(dp) :: g2, gsym2
 complex(dp) :: g, gsym
 complex(dp) :: bnum, beta
 class(ebands_t), pointer :: ebands
 class(gqk_t), pointer :: gqk
 logical :: isirr_k, isirr_q, isirr_kq
 logical :: isirr_sk, isirr_sq, isirr_skq
!arrays
 integer :: sym_k(3,3), invsym_k(3,3), g0_k(3)
 integer :: sym_q(3,3), invsym_q(3,3), g0_q(3)
 integer :: sym_kq(3,3), invsym_kq(3,3), g0_kq(3)
 !DEBUG
 integer :: sym_sk(3,3), invsym_sk(3,3), g0_sk(3)
 integer :: sym_sq(3,3), invsym_sq(3,3), g0_sq(3)
 integer :: sym_skq(3,3), invsym_skq(3,3), g0_skq(3)
 integer, allocatable :: qpk_map(:,:), kpq_map(:,:)
 integer, allocatable :: ib_deginds(:), jb_deginds(:), pert_deginds(:)
 real(dp), allocatable :: ik_eig(:), ikq_eig(:), pert_freq(:)
 real(dp) :: kpt(3), qpt(3), Sqpt(3), Skpt(3), mqpt(3), kpqpt(3)
 real(dp) :: kpt_ibz(3), qpt_ibz(3)
 real(dp), allocatable :: pc(:, :, :)
 complex(dp), allocatable :: prev_a(:, :, :), prev_b(:, :), prev_grad(:, :, :)
 complex(dp), allocatable :: gk_gathered(:,:,:,:)
 complex(dp), allocatable :: gk_symgathered(:,:,:,:)

! *************************************************************************

 ebands => self%gstore%ebands
 gqk => self%gstore%gqk(1)
 comm = self%gstore%comm
 nproc = xmpi_comm_size(comm)
 my_rank = xmpi_comm_rank(comm)

 ABI_MALLOC(prev_a, (self%nb, self%nkpt, self%nsppol))
 ABI_MALLOC(prev_grad, (self%nb, self%nkpt, self%nsppol))
 ABI_MALLOC(pc, (self%nb, self%nkpt, self%nsppol))
 ABI_MALLOC(prev_b, (gqk%my_npert, gqk%my_nq))

 write(ab_out,'(a4,a13,2a12,a13,a13,a14/)') 'Step', 'E_pol', 'E_el', &
   'E_ph', 'E_elph', 'epsilon', '|gradient|^2'

 call self%seed_a()
 do ii=1,150
   call self%b_from_a()
   prev_a(:, :, :) = self%a(:, :, :)
   prev_b(:, :) = self%b(:, :)
   call self%get_enterms()
   call self%get_elgrad()
   gradres = sum(abs(self%grad_a(:, :, :))**2)


   if (mod(ii, 50) == 1) then
     do is=1,self%nsppol
       bstart = self%gstore%brange_spin(1, is)
       do ik=1,self%nkpt
         ik_ibz = self%map_k_kibz(1, ik)
         do ib=1,self%nb
           pc(ib, ik, is) = one/(-ebands%eig(bstart+ib-1, ik_ibz, is) - half*self%eps)
         enddo
       enddo
     enddo
   endif

   write(ab_out,'(i4,es13.5,2es12.4,es13.4,es13.4,es14.4)') ii, &
     sum(self%enterms(:)), self%enterms(1), self%enterms(2), self%enterms(3), &
     self%eps, gradres

   self%grad_a(:, :, :) = self%grad_a(:, :, :)*pc(:, :, :)
   call self%orth_elgrad()

   if (ii > 1) then
     bnum = sum(self%grad_a(:, :, :)*conjg(self%grad_a(:, :, :) - prev_grad(:, :, :)))
     bden = sum(abs(prev_grad(:, :, :))**2)
     beta = bnum/bden
     self%grad_a(:, :, :) = self%grad_a(:, :, :) + beta*prev_grad(:, :, :)
     call self%orth_elgrad()
   endif
   prev_grad(:, :, :) = self%grad_a(:, :, :)

   call self%get_mintheta()
   costh = cos(self%theta); sinth = sin(self%theta)
   self%a(:, :, :) = costh*prev_a(:, :, :) + sinth*self%grad_a(:, :, :)

 enddo
 ABI_FREE(prev_a)
 ABI_FREE(prev_grad)
 ABI_FREE(pc)
 ABI_FREE(prev_b)



! testing the g(Sk,q) = g(k,S^{-1}q) symmetry
 ABI_MALLOC(qpk_map, (6, self%nkpt))
 do is=1,self%nsppol
   my_is = self%gstore%spin2my_is(is); if (my_is == 0) cycle
   gqk => self%gstore%gqk(my_is)
   bstart = self%gstore%brange_spin(1, is)

   ABI_MALLOC(ik_eig, (self%nb))
   ABI_MALLOC(ikq_eig, (self%nb))
   ABI_MALLOC(pert_freq, (gqk%my_npert))

   do ik=1,self%nkpt
     kpt = self%kpts(:, ik)
     ik_ibz = self%map_k_kibz(1, ik)

     ! eigenergies at ik
     ik_eig(:) = -ebands%eig(bstart:bstart+self%nb-1, ik_ibz, is)

     ! for this k' q+k'->k scattering map
     call self%get_mapping("q+k'->k", kpt, qpk_map)

     ! for this k' get S^{-1} such as k'=S*k_ibz
     call self%get_invsym_khalf(ik, ik_ibz_g, tsign_k, invsym_k, g0_k, isirr_k)
     kpt_ibz(:) = self%kpts(:, ik_ibz_g)

     call self%get_invsym_khalf(ik_ibz_g, isk_ibz_g, tsign_sk, invsym_sk, g0_sk, isirr_sk)

     ! for this k': gather all matrix elements distributed over q
     call gqk%gather("k", ik, gk_gathered)
     call gqk%gather("k", ik_ibz_g, gk_symgathered)

     do my_iq=1,gqk%my_nq
       iq_glob = my_iq + gqk%my_qstart - 1
       qpt(:) = self%qpts(:, iq_glob)

       ! eigenergies at ik+q
       ikq = qpk_map(1, iq_glob)
       ikq_ibz = self%map_k_kibz(1, ikq)
       ikq_eig(:) = -ebands%eig(bstart:bstart+self%nb-1, ikq_ibz, is)

       call self%get_invsym_khalf(ikq, ikq_ibz_g, tsign_kq, invsym_kq, g0_kq, isirr_kq)
       call self%get_invsym_qhalf(iq_glob, iq_ibz_g, tsign_q, invsym_q, g0_q, isirr_q)

       ! phonon frequencies at iq
       pert_freq(:) = gqk%my_wnuq(:, my_iq)

       ! S^{-1}q index in global array
       Sqpt(:) = tsign_k*matmul(invsym_k, qpt - g0_k)
       iSq_glob = self%krank_qpts%get_index(Sqpt)

       kpqpt(:) = kpt_ibz(:) + Sqpt(:)
       iskq = self%krank_kpts_g%get_index(kpqpt)

       call self%get_invsym_khalf(iskq, iskq_ibz_g, tsign_skq, invsym_skq, g0_skq, isirr_skq)
       call self%get_invsym_qhalf(iSq_glob, isq_ibz_g, tsign_sq, invsym_sq, g0_sq, isirr_sq)

       do ib=1,self%nb
         ! degenerate states indices at ib,k
         call get_deginds(ib, ik_eig, self%nb, ib_deginds, ib_ndeg)

         do jb=1,self%nb
           ! degenerate states indices at jb,ik+q
           call get_deginds(jb, ikq_eig, self%nb, jb_deginds, jb_ndeg)

           do my_pert=1,gqk%my_npert
             ! degenerate phonons at pert,iq
             call get_deginds(my_pert, pert_freq, gqk%my_npert, &
               pert_deginds, pert_ndeg)


             g = zero; gsym = zero
             g2 = zero; gsym2 = zero
             do ii=1,ib_ndeg
               do jj=1,jb_ndeg
                 do pp=1,pert_ndeg
                   g = g + gk_gathered(pert_deginds(pp), jb_deginds(jj), &
                     ib_deginds(ii), iq_glob)
                   gsym = gsym + gk_symgathered(pert_deginds(pp), &
                     jb_deginds(jj), ib_deginds(ii), iSq_glob)

                   g2 = g2 + abs(gk_gathered(pert_deginds(pp), jb_deginds(jj), &
                     ib_deginds(ii), iq_glob))**2
                   gsym2 = gsym2 + abs(gk_symgathered(pert_deginds(pp), &
                     jb_deginds(jj), ib_deginds(ii), iSq_glob))**2
                 enddo
               enddo
             enddo

             g = g/(ib_ndeg*jb_ndeg*pert_ndeg)
             gsym = gsym/(ib_ndeg*jb_ndeg*pert_ndeg)

             g2 = g2/(ib_ndeg*jb_ndeg*pert_ndeg)
             gsym2 = gsym2/(ib_ndeg*jb_ndeg*pert_ndeg)

             if (abs(g2 - gsym2) > tol8) then
               write(ab_out, '(3i3,6f7.3,2es13.4)') my_pert, jb, ib, kpt(:), qpt(:), &
                 g2, gsym2
             endif

             g = gk_gathered(my_pert, jb, ib, iq_glob)
             gsym = gk_symgathered(my_pert, jb, ib, iSq_glob)
             !if (tsign_k == -1 .and. tsign_sk == 1) then
             !  if ((tsign_q == -1) .and. (tsign_kq == 1) .and. (tsign_sq == 1) .and. (tsign_skq == 1)) then
             !     gsym = -gk_symgathered(my_pert, jb, ib, iSq_glob)
             !  else if ((tsign_q == 1) .and. (tsign_kq == -1) .and. (tsign_sq == 1) .and. (tsign_skq == 1)) then
             !     gsym = -gk_symgathered(my_pert, jb, ib, isq_glob)
             !  else if ((tsign_q == 1) .and. (tsign_kq == 1) .and. (tsign_sq == -1) .and. (tsign_skq == 1)) then
             !     gsym = -gk_symgathered(my_pert, jb, ib, iSq_glob)
             !  else if ((tsign_q == 1) .and. (tsign_kq == 1) .and. (tsign_sq == 1) .and. (tsign_skq == -1)) then
             !     gsym = -gk_symgathered(my_pert, jb, ib, iSq_glob)
             !  else
             !     gsym = gk_symgathered(my_pert, jb, ib, iSq_glob)
             !  endif
             !else
             !  gsym = gk_symgathered(my_pert, jb, ib, iSq_glob)
             !endif

             !if (ib_ndeg == 1 .and. jb_ndeg == 1 .and. pert_ndeg == 1) then
             !if (ib_ndeg == 1 .and. pert_ndeg == 1) then
             !if (ib_ndeg == 1 .and. jb_ndeg == 1) then
             !if (jb_ndeg == 1 .and. pert_ndeg == 1) then
               !if (abs(abs(g)**2 - abs(gsym)**2) > tol8) then
               !if (abs(g - gsym)**2 > tol12) then
               !if (abs(g - gsym)**2 > tol12) then
               !  write(ab_out, '(4es13.4,a,3i3,a,3i3,a)') g, gsym, &
               !    '(', tsign_k, tsign_q, tsign_kq, ') (', tsign_sk, tsign_sq, tsign_skq, ')'

               !  !write(ab_out, *) '(', isirr_k, isirr_q, isirr_kq, ') (', isirr_sk, isirr_sq, isirr_skq, ')'

               !endif
             !endif

             ABI_FREE(pert_deginds)
           enddo
           ABI_FREE(jb_deginds)
         enddo
         ABI_FREE(ib_deginds)
       enddo
     enddo
     ABI_FREE(gk_gathered)
     ABI_FREE(gk_symgathered)
   enddo
   ABI_FREE(ik_eig)
   ABI_FREE(ikq_eig)
   ABI_FREE(pert_freq)
 enddo
 ABI_FREE(qpk_map)

end subroutine varpeq_test
!!***


!!****f* m_varpeq/varpeq_orth_elgrad
!! NAME
!!  varpeq_orth_elgrad
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_orth_elgrad(self)

!Arguments ------------------------------------
!scalars
 class(varpeq_t), intent(inout) :: self
!arrays

!Local variables-------------------------------
 real(dp) :: normgrad2
 real(dp) :: nkbzinv
 real(dp) :: factor
 complex(dp) :: inprod

! *************************************************************************

 inprod = zero
 normgrad2 = zero
 nkbzinv = one/self%nkbz

 inprod = sum(self%a(:,:,:)*conjg(self%grad_a(:,:,:)))
 factor = inprod*nkbzinv
 self%grad_a(:, :, :) = self%grad_a(:, :, :) - factor*self%a(:, :, :)

 call norm_a(self%grad_a, self%nsppol, self%nkpt, self%nb, self%nkbz, self%wtk)

end subroutine varpeq_orth_elgrad
!!***


!!****f* m_varpeq/varpeq_get_mintheta
!! NAME
!!  varpeq_get_mintheta
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_get_mintheta(self)

!Arguments ------------------------------------
!scalars
 class(varpeq_t), intent(inout) :: self

!Local variables-------------------------------
!scalars
 integer :: comm, nproc, my_rank, ierr
 integer :: bstart
 integer :: is, ik, ib, jb
 integer :: ik_ibz, ik_ibz_g, tsign_k
 integer :: ik_forw, ikq_ibz
 integer :: my_iq, my_pert, iq_glob, my_is
 integer :: iSq_glob
 complex(dp) :: a_from, a_forw
 complex(dp) :: d_from, d_forw
 real(dp) :: kweight
 real(dp) :: eig
 complex(dp) :: bqnu, g_forw
 complex(dp) :: term3_tmp, term4_tmp
 class(ebands_t), pointer :: ebands
 class(gqk_t), pointer :: gqk
 logical :: isirr
!arrays
 integer, allocatable :: qpk_map(:, :)
 integer :: invsym_k(3,3), g0_k(3)
 real(dp) :: kpt(3), qpt(3), Sqpt(3)
 real(dp) :: terms(4)
 complex(dp), allocatable :: gk_gathered(:,:,:,:)

! *************************************************************************

 comm = self%gstore%comm
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 gqk => self%gstore%gqk(1)
 ebands => self%gstore%ebands

 terms(:) = zero

 ! non-scattering dependent terms
 do is=1,self%nsppol
   bstart = self%gstore%brange_spin(1, is)
   do ik=1,self%nkpt
     ik_ibz = self%map_k_kibz(1, ik)
     kweight = self%wtk(ik)
     do ib=1,self%nb
       eig = -ebands%eig(bstart+ib-1, ik_ibz, is)
       a_from = self%a(ib, ik, is)
       d_from = self%grad_a(ib, ik, is)
       ! sin^2(theta) term
       terms(1) = terms(1) + kweight*eig*abs(d_from)**2
       ! sin(theta)cos(theta) term
       terms(2) = terms(2) + &
         kweight*eig*real(d_from*conjg(a_from) + conjg(d_from)*a_from)
     enddo
   enddo
 enddo

 ! scattering dependent terms
 term3_tmp = zero; term4_tmp = zero
 ABI_MALLOC(qpk_map, (6, self%nkpt))

 do is=1,self%nsppol
   my_is = self%gstore%spin2my_is(is); if (my_is == 0) cycle
   gqk => self%gstore%gqk(my_is)

   ! forward scattering part
   do ik=1,self%nkpt
     kpt = self%kpts(:, ik)
     kweight = self%wtk(ik)
     ik_ibz = self%map_k_kibz(1, ik)

     ! for this k' q+k'->k scattering map
     call self%get_mapping("q+k'->k", kpt, qpk_map)

     ! for this k' get S^{-1} such as k'=S*k_ibz
     !call self%get_invsym_k(ik, ik_ibz_g, tsign_k, invsym_k, g0_k, isirr)
     call self%get_invsym_khalf(ik, ik_ibz_g, tsign_k, invsym_k, g0_k, isirr)

     ! for this k': gather all matrix elements distributed over q
     !call gqk%gather("k", ik, gk_gathered)
     call gqk%gather("k", ik_ibz_g, gk_gathered)

     do my_iq=1,gqk%my_nq
       iq_glob = my_iq + gqk%my_qstart - 1
       qpt(:) = self%qpts(:, iq_glob)

       ! S^{-1}q index in global array
       Sqpt(:) = tsign_k*matmul(invsym_k, qpt - g0_k)
       iSq_glob = self%krank_qpts%get_index(Sqpt)

       ! k+q index
       ik_forw = qpk_map(1, iq_glob)
       ikq_ibz = self%map_k_kibz(1, ik_forw)

       do ib=1,self%nb
         a_from = self%a(ib, ik, is)
         d_from = self%grad_a(ib, ik, is)

         do jb=1,self%nb
           a_forw = self%a(jb, ik_forw, is)
           d_forw = self%grad_a(jb, ik_forw, is)

           do my_pert=1,gqk%my_npert
             bqnu = self%bq_glob(my_pert, iq_glob)
             !g_forw = gk_gathered(my_pert, jb, ib, iq_glob)
             g_forw = gk_gathered(my_pert, jb, ib, iSq_glob)
             !g_forw = self%get_gavg(gk_gathered, gqk, my_pert, jb, ib, ik_ibz, &
             !  my_iq, iq_glob, ikq_ibz)

             ! sin^2(theta) term
             term3_tmp = term3_tmp + &
               kweight*d_from*conjg(bqnu)*g_forw*conjg(d_forw)
             ! sin(theta)cos(theta) term
             term4_tmp = term4_tmp + &
               kweight*conjg(bqnu)*g_forw*(a_from*conjg(d_forw) + &
               d_from*conjg(a_forw))

           enddo
         enddo
       enddo
     enddo

     ABI_FREE(gk_gathered)
   enddo
 enddo

 ABI_FREE(qpk_map)

 term3_tmp = term3_tmp + conjg(term3_tmp)
 term4_tmp = term4_tmp + conjg(term4_tmp)
 call xmpi_sum(term3_tmp, comm, ierr)
 call xmpi_sum(term4_tmp, comm, ierr)
 terms(3) = -real(term3_tmp) / (one*self%nqpt)
 terms(4) = -real(term4_tmp) / (one*self%nqpt)

 self%theta = half*atan2(-(terms(2) + terms(4)), (terms(1) + terms(3) - self%eps))

end subroutine varpeq_get_mintheta
!!***


!!****f* m_varpeq/varpeq_get_elgrad
!! NAME
!!  varpeq_get_elgrad
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_get_elgrad(self)

!Arguments ------------------------------------
!scalars
 class(varpeq_t), intent(inout) :: self

!Local variables-------------------------------
!scalars
 integer :: comm, nproc, my_rank, ierr
 integer :: bstart
 integer :: is, ik, jk, ib, jb
 integer :: ik_ibz, ik_ibz_g, tsign_k
 integer :: ik_forw, ik_back, ikq_ibz
 integer :: my_iq, my_pert, iq_glob, my_is
 integer :: iSq_glob
 real(dp) :: eig
 real(dp) :: nkbzinv, nkbz2inv
 complex(dp) grad_tmp
 complex(dp) :: a_forw, a_back
 complex(dp) :: bqnu, g_forw, g_back
 class(ebands_t), pointer :: ebands
 class(gqk_t), pointer :: gqk
 logical :: isirr
!arrays
 integer, allocatable :: qpk_map(:, :), mqpk_map(:,:), kmk_map(:, :)
 integer :: invsym_k(3,3), g0_k(3)
 real(dp) :: kpt(3), qpt(3), Sqpt(3)
 complex(dp), allocatable :: gk_gathered(:,:,:,:)

! *************************************************************************

 comm = self%gstore%comm
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 gqk => self%gstore%gqk(1)
 ebands => self%gstore%ebands

 ABI_MALLOC(qpk_map, (6, self%nqpt))
 ABI_MALLOC(mqpk_map, (6, self%nqpt))
 ABI_MALLOC(kmk_map, (6, self%nkpt))

 nkbzinv = one/(one*self%nkbz)
 nkbz2inv = nkbzinv**2

 self%grad_a(:, :, :) = zero

 do is=1,self%nsppol
   my_is = self%gstore%spin2my_is(is); if (my_is == 0) cycle
   gqk => self%gstore%gqk(my_is)

   ! forward scattering part
   do ik=1,self%nkpt
     kpt = self%kpts(:, ik)
     ik_ibz = self%map_k_kibz(1, ik)

     ! for this k' q+k'->k scattering map
     call self%get_mapping("q+k'->k", kpt, qpk_map)

     ! for this k' q+k'->k scattering map
     call self%get_mapping("-q+k'->k", kpt, mqpk_map)

     ! for this k' get S^{-1} such as k'=S*k_ibz
     !call self%get_invsym_k(ik, ik_ibz_g, tsign_k, invsym_k, g0_k, isirr)
     call self%get_invsym_khalf(ik, ik_ibz_g, tsign_k, invsym_k, g0_k, isirr)

     ! for this k': gather all matrix elements distributed over q
     !call gqk%gather("k", ik, gk_gathered)
     call gqk%gather("k", ik_ibz_g, gk_gathered)

     do my_iq=1,gqk%my_nq
       iq_glob = my_iq + gqk%my_qstart - 1
       qpt(:) = self%qpts(:, iq_glob)

       ! S^{-1}q index in global array
       Sqpt(:) = tsign_k*matmul(invsym_k, qpt - g0_k)
       iSq_glob = self%krank_qpts%get_index(Sqpt)

       ! k+q index
       ik_forw = qpk_map(1, iq_glob)
       ikq_ibz = self%map_k_kibz(1, ik_forw)

       ! k-q index
       ik_back = mqpk_map(1, iq_glob)

       do ib=1,self%nb
         do jb=1,self%nb
           a_forw = self%a(jb, ik_forw, is)
           a_back = self%a(jb, ik_back, is)

           do my_pert=1,gqk%my_npert
             bqnu = self%bq_glob(my_pert, iq_glob)
             !g_forw = gk_gathered(my_pert, jb, ib, iq_glob)
             g_forw = gk_gathered(my_pert, jb, ib, iSq_glob)
             !g_forw = self%get_gavg(gk_gathered, gqk, my_pert, jb, ib, ik_ibz, &
             !  my_iq, iq_glob, ikq_ibz)
             !g_back = gqk%my_g(my_pert, ib, my_iq, jb, ik)


             self%grad_a(ib, ik, is) = self%grad_a(ib, ik, is) - &
               two*nkbz2inv*a_forw*bqnu*conjg(g_forw)
           enddo
         enddo
       enddo
     enddo

     ABI_FREE(gk_gathered)
   enddo
 enddo
 call xmpi_sum(self%grad_a, comm, ierr)


 do is=1,self%nsppol
   my_is = self%gstore%spin2my_is(is); if (my_is == 0) cycle
   gqk => self%gstore%gqk(my_is)

   ! backwards scattering part
   do ik=1,self%nkpt
     kpt = self%kpts(:, ik)
     ik_ibz = self%map_k_kibz(1, ik)

     ! for this k' k-k'->q scattering map
     call self%get_mapping("k+k'->q", -kpt, kmk_map)

     ! for this k' get S^{-1} such as k'=S*k_ibz
     !call self%get_invsym_k(ik, ik_ibz_g, tsign_k, invsym_k, g0_k, isirr)
     call self%get_invsym_khalf(ik, ik_ibz_g, tsign_k, invsym_k, g0_k, isirr)

     ! for this k': gather all matrix elements distributed over q
     !call gqk%gather("k", ik, gk_gathered)
     call gqk%gather("k", ik_ibz_g, gk_gathered)

     do jk=1,self%nkpt
       iq_glob = kmk_map(1, jk)
       qpt(:) = self%qpts(:, iq_glob)
       ikq_ibz = self%map_k_kibz(1, jk)

       ! S^{-1}q index in global array
       Sqpt(:) = tsign_k*matmul(invsym_k, qpt - g0_k)
       iSq_glob = self%krank_qpts%get_index(Sqpt)

       do ib=1,self%nb
         a_back = self%a(ib, ik, is)
         do jb=1,self%nb
           do my_pert=1,gqk%my_npert
             bqnu = self%bq_glob(my_pert, iq_glob)
             !g_back = gk_gathered(my_pert, jb, ib, iq_glob)
             g_back = gk_gathered(my_pert, jb, ib, iSq_glob)
             !g_back = self%get_gavg(gk_gathered, gqk, my_pert, jb, ib, ik_ibz, &
             !  my_iq, iq_glob, ikq_ibz)

             self%grad_a(jb, jk, is) = self%grad_a(jb, jk, is) - &
               two*nkbz2inv*a_back*conjg(bqnu)*g_back
           enddo
         enddo
       enddo
     enddo

     ABI_FREE(gk_gathered)
   enddo
 enddo


 ! Scattering-independent part
 do is=1,self%nsppol
   bstart = self%gstore%brange_spin(1, is)
   do ik=1,self%nkpt
     ik_ibz = self%map_k_kibz(1, ik)
     do ib=1,self%nb
       eig = -ebands%eig(bstart+ib-1, ik_ibz, is)
       self%grad_a(ib, ik, is) = self%grad_a(ib, ik, is) + &
         two*nkbzinv*(eig - self%eps)*self%a(ib, ik, is)
     enddo
   enddo
 enddo

 ABI_FREE(qpk_map)
 ABI_FREE(mqpk_map)
 ABI_FREE(kmk_map)

end subroutine varpeq_get_elgrad
!!***

!!****f* m_varpeq/varpeq_get_enterms
!! NAME
!!  varpeq_get_enterms
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_get_enterms(self)

!Arguments ------------------------------------
!scalars
 class(varpeq_t), intent(inout) :: self

!Local variables-------------------------------
!scalars
 integer :: comm, nproc, my_rank, ierr
 integer :: bstart
 integer :: is, ik, ib, jb
 integer :: ik_ibz, ik_ibz_g, tsign_k
 integer :: ik_forw, ikq_ibz
 integer :: my_iq, my_pert, iq_glob, my_is
 integer :: iSq_glob
 complex(dp) :: a_from, a_forw
 real(dp) :: enel, envib, enelph
 real(dp) :: kweight, phfreq
 real(dp) :: eig
 complex(dp) :: enelph_tmp
 complex(dp) :: bqnu, g_forw
 class(ebands_t), pointer :: ebands
 class(gqk_t), pointer :: gqk
 logical :: isirr
!arrays
 integer, allocatable :: qpk_map(:, :)
 integer :: invsym_k(3,3), g0_k(3)
 real(dp) :: kpt(3), qpt(3), Sqpt(3)
 complex(dp), allocatable :: gk_gathered(:,:,:,:)

! *************************************************************************

 comm = self%gstore%comm
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 gqk => self%gstore%gqk(1)
 ebands => self%gstore%ebands

 ! Electronic energy
 enel = zero
 do is=1,self%nsppol
   bstart = self%gstore%brange_spin(1, is)
   do ik=1,self%nkpt
     ik_ibz = self%map_k_kibz(1, ik)
     kweight = self%wtk(ik)
     do ib=1,self%nb
       eig = -ebands%eig(bstart+ib-1, ik_ibz, is)
       enel = enel + kweight*eig*abs(self%a(ib, ik, is))**2
     enddo
   enddo
 enddo

 ! Vibrational energy
 envib = zero
 do my_iq=1,gqk%my_nq
   iq_glob = my_iq + gqk%my_qstart - 1
   do my_pert=1,gqk%my_npert
     phfreq = gqk%my_wnuq(my_pert, my_iq)
     envib = envib + phfreq*abs(self%bq_glob(my_pert, iq_glob))**2
   enddo
 enddo
 call xmpi_sum(envib, comm, ierr)
 envib = envib/(one*self%nqpt)

 ! Electron-phonon energy
 enelph_tmp = (zero, zero)
 ABI_MALLOC(qpk_map, (6, self%nqpt))
 do is=1,self%nsppol
   my_is = self%gstore%spin2my_is(is); if (my_is == 0) cycle
   gqk => self%gstore%gqk(my_is)

   do ik=1,self%nkpt
     kweight = self%wtk(ik)
     kpt = self%kpts(:, ik)
     ik_ibz = self%map_k_kibz(1, ik)

     ! for this k': q+k'-> k scattering map
     call self%get_mapping("q+k'->k", kpt, qpk_map)

     ! for this k' get S^{-1} such as k'=S*k_ibz
     !call self%get_invsym_k(ik, ik_ibz_g, tsign_k, invsym_k, g0_k, isirr)
     call self%get_invsym_khalf(ik, ik_ibz_g, tsign_k, invsym_k, g0_k, isirr)

     ! for this k': gather all matrix elements distributed over q
     !call gqk%gather("k", ik, gk_gathered)
     call gqk%gather("k", ik_ibz_g, gk_gathered)

     do my_iq=1,gqk%my_nq
       iq_glob = my_iq + gqk%my_qstart - 1
       qpt(:) = self%qpts(:, iq_glob)

       ! S^{-1}q index in global array
       Sqpt(:) = tsign_k*matmul(invsym_k, qpt - g0_k)
       iSq_glob = self%krank_qpts%get_index(Sqpt)

       ! k+q index
       ik_forw = qpk_map(1, iq_glob)
       ikq_ibz = self%map_k_kibz(1, ik_forw)

       do ib=1,self%nb
         a_from = self%a(ib, ik, is)
         do jb=1,self%nb
           a_forw = self%a(jb, ik_forw, is)
           do my_pert=1,gqk%my_npert
             bqnu = self%bq_glob(my_pert, iq_glob)
             !g_forw = gk_gathered(my_pert, jb, ib, iq_glob)
             g_forw = gk_gathered(my_pert, jb, ib, iSq_glob)
             !g_forw = self%get_gavg(gk_gathered, gqk, my_pert, jb, ib, ik_ibz, &
             !  my_iq, iq_glob, ikq_ibz)

             enelph_tmp = enelph_tmp + &
               kweight*a_from*conjg(bqnu)*g_forw*conjg(a_forw)

           enddo
         enddo
       enddo

     enddo
     ABI_FREE(gk_gathered)
   enddo
 enddo
 ABI_FREE(qpk_map)

 enelph_tmp = enelph_tmp + conjg(enelph_tmp)
 call xmpi_sum(enelph_tmp, comm, ierr)
 enelph = -real(enelph_tmp) / (one*self%nqpt)

 ! Setting the corresponding state variables of the varpeq type
 self%enterms(1) = enel; self%enterms(2) = envib; self%enterms(3) = enelph
 self%eps = enel + enelph

end subroutine varpeq_get_enterms
!!***


!!****f* m_varpeq/varpeq_b_from_a
!! NAME
!!  varpeq_b_from_a
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_b_from_a(self)

!Arguments ------------------------------------
!scalars
 class(varpeq_t), intent(inout) :: self

!Local variables-------------------------------
!scalars
 integer :: ierr
 integer :: is, ik, ib, jb
 integer :: my_iq, my_pert, my_is
 integer :: iq_glob, iSq_glob
 integer :: ik_ibz_g, tsign_k
 integer :: ik_forw, ik_ibz, ikq_ibz
 real(dp) :: kweight
 complex(dp) :: a_from, a_forw
 complex(dp) :: g_forw
 class(gqk_t), pointer :: gqk
 logical :: isirr
!arrays
 integer, allocatable :: qpk_map(:,:)
 integer :: invsym_k(3,3), g0_k(3)
 real(dp) :: kpt(3), qpt(3), Sqpt(3)
 complex(dp), allocatable :: gk_gathered(:,:,:,:)

! *************************************************************************

 ABI_MALLOC(qpk_map, (6, self%nqpt))

 self%bq_glob(:,:) = zero
 do is=1,self%nsppol
   my_is = self%gstore%spin2my_is(is); if (my_is == 0) cycle
   gqk => self%gstore%gqk(my_is)

   do ik=1,self%nkpt
     kweight = self%wtk(ik)
     kpt = self%kpts(:, ik)
     ik_ibz = self%map_k_kibz(1, ik)

     ! for this k': q+k'-> k scattering map
     call self%get_mapping("q+k'->k", kpt, qpk_map)

     ! for this k' get S^{-1} such as k'=S*k_ibz
     !call self%get_invsym_k(ik, ik_ibz_g, tsign_k, invsym_k, g0_k, isirr)
     call self%get_invsym_khalf(ik, ik_ibz_g, tsign_k, invsym_k, g0_k, isirr)

     ! for this k': gather all matrix elements distributed over q
     !call gqk%gather("k", ik, gk_gathered)
     call gqk%gather("k", ik_ibz_g, gk_gathered)

     do my_iq=1,gqk%my_nq
       iq_glob = my_iq + gqk%my_qstart - 1
       qpt(:) = self%qpts(:, iq_glob)

       ! S^{-1}q index in global array
       Sqpt(:) = tsign_k*matmul(invsym_k, qpt - g0_k)
       iSq_glob = self%krank_qpts%get_index(Sqpt)

       ! k+q index
       ik_forw = qpk_map(1, iq_glob)
       ikq_ibz = self%map_k_kibz(1, ik_forw)

       do ib=1,self%nb
         a_from = self%a(ib, ik, is)
         do jb=1,self%nb
           a_forw = self%a(jb, ik_forw, is)
           do my_pert=1,gqk%my_npert
             !g_forw = gk_gathered(my_pert, jb, ib, iq_glob)
             g_forw = gk_gathered(my_pert, jb, ib, iSq_glob)
             !g_forw = self%get_gavg(gk_gathered, gqk, my_pert, jb, ib, ik_ibz, &
             !  my_iq, iq_glob, ikq_ibz)

             if (gqk%my_wnuq(my_pert, my_iq) == 0) then
               self%bq_glob(my_pert, iq_glob) = zero
             else
               self%bq_glob(my_pert, iq_glob) = self%bq_glob(my_pert, iq_glob) + &
                 conjg(a_forw)*g_forw*a_from*kweight/gqk%my_wnuq(my_pert, my_iq)
             endif

           enddo
         enddo
       enddo
     enddo
     ABI_FREE(gk_gathered)
   enddo
   call xmpi_sum(self%bq_glob, gqk%qpt_comm%value, ierr)
 enddo
 ABI_FREE(qpk_map)

end subroutine varpeq_b_from_a
!!***


!!****f* m_varpeq/varpeq_seed_a
!! NAME
!!  varpeq_seed_a
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_seed_a(self)

!Arguments ------------------------------------
!scalars
 class(varpeq_t), intent(inout) :: self

!Local variables-------------------------------
!scalars
 integer :: is, ik, ib
 integer :: ik_ibz
 integer :: bstart
 complex(dp) :: imeig
 class(ebands_t), pointer :: ebands

! *************************************************************************

 ebands => self%gstore%ebands

 do is=1,self%nsppol
   bstart = self%gstore%brange_spin(1, is)
   do ik=1,self%nkpt
     ik_ibz = self%map_k_kibz(1, ik)
     do ib=1,self%nb
       imeig = -one*(ebands%eig(bstart+ib-1, ik_ibz, is))
       self%a(ib, ik, is) = exp(-imeig)
     enddo
   enddo
 enddo

 call norm_a(self%a, self%nsppol, self%nkpt, self%nb, self%nkbz, self%wtk)

end subroutine varpeq_seed_a
!!***


!!****f* m_varpeq/varpeq_init
!! NAME
!!  varpeq_init
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_init(self, gstore, dtset)

!Arguments ------------------------------------
!scalars
 class(varpeq_t), intent(inout) :: self
 type(dataset_type), intent(in) :: dtset
 type(gstore_t), target, intent(in) :: gstore
!arrays

!Local variables-------------------------------
!scalars
 integer :: comm, nproc, my_rank, ierr
 integer :: kptopt, nkbz
 integer :: my_iq, iq_glob
 integer :: my_ik, ik_glob
 integer :: my_pert
 real(dp) :: wtq, my_freq
 type(crystal_t) :: cryst_fake
 class(crystal_t), pointer :: cryst
 class(ebands_t), pointer :: ebands
 class(gqk_t), pointer :: gqk
 type(krank_t) :: krank_ibz
!arrays
 real(dp) :: my_qpt(3), my_kpt(3)
 real(dp), allocatable :: kbz(:, :)

! *************************************************************************

 comm = gstore%comm
 nproc = xmpi_comm_size(comm)
 my_rank = xmpi_comm_rank(comm)

 ! Consistency check
 ierr = 0
 ABI_CHECK_NOSTOP(gstore%kzone == "bz", "kzone = 'bz' is required", ierr)
 ABI_CHECK_NOSTOP(gstore%qzone == "bz", "qzone = 'bz' is required", ierr)
 ABI_CHECK_NOSTOP(gstore%gqk(1)%kpt_comm%nproc == 1, "varpeq is not &
   compatible with parallelism over k-points", ierr)
 ABI_CHECK_NOSTOP(gstore%gqk(1)%cplex == 2, "cplex = 2 is required", ierr)
 ABI_CHECK(ierr == 0, "The gstore object is incosistent with varpeq. &
   See messages above.")
 self%gstore => gstore

 ! Initalization of the electronic subspace & related symmetry tables
 cryst => gstore%cryst
 ebands => gstore%ebands
 self%nsppol = gstore%nsppol
 self%nb = gstore%gqk(1)%nb

 ! Do not use symmetries for the electronic and phonon coefficients
 ! TODO: add support for A_{-k}^* = A_k & B_{-q}^* = B_q symmetries

 kptopt = 3
 self%cryst = cryst%new_without_symmetries()
 self%timrev = kpts_timrev_from_kptopt(kptopt)

 ! Generate k-points in the effective BZ and allocate electronic vector
 call kpts_ibz_from_kptrlatt(cryst, ebands%kptrlatt, kptopt, ebands%nshiftk, &
   ebands%shiftk, self%nkpt, self%kpts, self%wtk, nkbz, kbz)
 self%nkbz = nkbz

 ABI_MALLOC(self%a, (self%nb, self%nkpt, ebands%nsppol))
 ABI_MALLOC(self%grad_a, (self%nb, self%nkpt, ebands%nsppol))

 ! krank for kpts
 self%krank_kpts = krank_from_kptrlatt(self%nkpt, self%kpts, ebands%kptrlatt, &
   compute_invrank = .True.)

 ! Mappings for k -> k_IBZ and -k -> k_IBZ wrt full set of symmetries
 ABI_MALLOC(self%map_k_kibz, (6, self%nkpt))
 ABI_MALLOC(self%map_invk_kibz, (6, self%nkpt))

 krank_ibz = krank_from_kptrlatt(ebands%nkpt, ebands%kptns, ebands%kptrlatt, &
   compute_invrank=.True.)

 if (kpts_map("symrel", 1, cryst, krank_ibz, self%nkpt, self%kpts, &
   self%map_k_kibz) /= 0) then
   ABI_ERROR("Varpeq initialization: cannot map k-space to IBZ")
 end if
 if (kpts_map("symrel", 1, cryst, krank_ibz, self%nkpt, -self%kpts, &
   self%map_invk_kibz) /=0) then
   ABI_ERROR("Varpeq initialization: cannot map inverse k-space to IBZ")
 end if

 ! DEBUG: k-points in half BZ
 kptopt = 2
 cryst_fake = cryst%new_inversion_only()
 call kpts_ibz_from_kptrlatt(cryst_fake, ebands%kptrlatt, kptopt, ebands%nshiftk, &
   ebands%shiftk, self%nkpt_half, self%kpts_half, self%wtk_half, nkbz, kbz)

 self%krank_kpts_half = krank_from_kptrlatt(self%nkpt_half, self%kpts_half, &
   ebands%kptrlatt, compute_invrank=.True.)

 ABI_MALLOC(self%map_k_khalf, (6, self%nkpt))
 if (kpts_map("symrel", 1, cryst_fake, self%krank_kpts_half, self%nkpt, self%kpts, &
   self%map_k_khalf) /= 0) then
   ABI_ERROR("Varpeq initialization: cannot map k-space to half BZ")
 end if


 ! Initialization of the vibrational subspace dimensions
 ! generate q-points in the full BZ
 gqk => self%gstore%gqk(1)
 self%nqpt = gstore%nqbz

 ABI_MALLOC(self%b, (gqk%my_npert, gqk%my_nq))
 ABI_MALLOC(self%bq_glob, (gqk%my_npert, self%nqpt))
 ABI_MALLOC(self%qpts, (3, self%nqpt))
 ABI_MALLOC(self%map_q_qibz, (6, self%nqpt))

 self%qpts(:, :) = zero
 self%map_q_qibz(:, :) = zero
 do my_iq=1,gqk%my_nq
   iq_glob = my_iq + gqk%my_qstart - 1
   call gqk%myqpt(my_iq, gstore, wtq, my_qpt)
   self%qpts(:, iq_glob) = my_qpt(:)
   self%map_q_qibz(:, iq_glob) = gqk%my_q2ibz(:, my_iq)
 enddo
 call xmpi_sum(self%qpts, gqk%qpt_comm%value, ierr)
 call xmpi_sum(self%map_q_qibz, gqk%qpt_comm%value, ierr)

 ! krank for qpts
 self%krank_qpts = krank_from_kptrlatt(self%nqpt, self%qpts, ebands%kptrlatt, &
   compute_invrank = .True.)

 ! setting k-zone from ggk
 self%nkpt_g = gqk%glob_nk
 ABI_MALLOC(self%kpts_g, (3, gqk%glob_nk))
 do my_ik=1,gqk%my_nk
   ik_glob = my_ik + gqk%my_kstart - 1
   my_kpt(:) = gqk%my_kpts(:, my_ik)
   self%kpts_g(:, ik_glob) = my_kpt(:)
 enddo
 call xmpi_sum(self%kpts_g, gqk%kpt_comm%value, ierr)
 self%krank_kpts_g = krank_from_kptrlatt(self%nkpt_g, self%kpts_g, ebands%kptrlatt, &
   compute_invrank = .True.)

! q-points in half BZ
 call kpts_ibz_from_kptrlatt(cryst_fake, ebands%kptrlatt, 2, ebands%nshiftk, &
   ebands%shiftk, self%nqpt_half, self%qpts_half, self%wtq_half, nkbz, kbz)

 self%krank_qpts_half = krank_from_kptrlatt(self%nqpt_half, self%qpts_half, &
   ebands%kptrlatt, compute_invrank=.True.)

 ABI_MALLOC(self%map_q_qhalf, (6, self%nqpt))
 if (kpts_map("symrec", 1, cryst_fake, self%krank_qpts_half, self%nqpt, self%qpts, &
   self%map_q_qhalf) /= 0) then
   ABI_ERROR("Varpeq initialization: cannot map q-space to half BZ")
 end if

 ABI_FREE(kbz)
 call krank_ibz%free()

end subroutine varpeq_init
!!***


!!****f* m_varpeq/varpeq_free
!! NAME
!!  varpeq_free
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_free(self)

!Arguments ------------------------------------
!scalars
 class(varpeq_t), intent(inout) :: self

!Local variables-------------------------------

! *************************************************************************

 ! Electronic subspace
 ABI_SFREE(self%kpts)
 ABI_SFREE(self%wtk)
 ABI_SFREE(self%a)
 ABI_SFREE(self%grad_a)

 ! Vibrational subspace
 ABI_SFREE(self%qpts)
 ABI_SFREE(self%b)
 ABI_SFREE(self%bq_glob)

 ! Symmetry tables
 ABI_SFREE(self%map_k_kibz)
 ABI_SFREE(self%map_invk_kibz)
 call self%krank_kpts%free()
 call self%krank_qpts%free()
 call self%cryst%free()

 ABI_SFREE(self%map_k_khalf)
 ABI_SFREE(self%kpts_half)
 ABI_SFREE(self%wtk_half)
 call self%krank_kpts_half%free()

 ABI_SFREE(self%map_q_qhalf)
 ABI_SFREE(self%qpts_half)
 ABI_SFREE(self%wtq_half)
 call self%krank_qpts_half%free()

end subroutine varpeq_free
!!***


!!****f* m_varpeq/varpeq_get_mapping
!! NAME
!!  vapreq_get_mapping
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_get_mapping(self, mode, pt, map)

!Arguments ------------------------------------
!scalars
 class(varpeq_t), target, intent(inout) :: self
 character(len=*),intent(in) :: mode
!arrays
 integer, intent(out) :: map(:,:)
 real(dp), intent(in) :: pt(3)

!Local variables-------------------------------
!scalars
 integer :: sigma
 character(len=6) :: symmode
 integer :: ierr
 integer :: nscatter
 type(krank_t), pointer :: scatter_to
!arrays
 real(dp), pointer :: scatter_from(:,:)

! *************************************************************************

 select case(mode)
 case("q+k'->k")
   scatter_to => self%krank_kpts
   scatter_from => self%qpts
   nscatter = self%nqpt
   symmode = "symrel"
   sigma = one
 case("-q+k'->k")
   scatter_to => self%krank_kpts
   scatter_from => self%qpts
   nscatter = self%nqpt
   symmode = "symrel"
   sigma = -one
 case("k+q'->k")
   scatter_to => self%krank_kpts
   scatter_from => self%kpts
   nscatter = self%nkpt
   symmode = "symrel"
   sigma = one
 case("k+k'->q")
   scatter_to => self%krank_qpts
   scatter_from => self%kpts
   nscatter = self%nkpt
   symmode = "symrec"
   sigma = one
 case default
   ABI_ERROR(sjoin("VarPEq, scattering mapping. Unsupported mode: ", mode))
 end select

 ierr = kpts_map(symmode, self%timrev, self%cryst, scatter_to, nscatter, &
   sigma*scatter_from, map, qpt=pt)

 if (ierr /= 0) then
   ABI_ERROR(sjoin("VarPEq: cannot map ", mode, "with point ", ktoa(pt)))
 endif

end subroutine varpeq_get_mapping
!!***


!!****f* m_varpeq/varpeq_get_invsym_khalf
!! NAME
!!  vapreq_get_invsym_khalf
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_get_invsym_khalf(self, ik, ik_ibz_g, tsign_k, invsym_k, g0_k, isirr)

!Arguments ------------------------------------
!scalars
 class(varpeq_t), intent(inout) :: self
 integer, intent(in) :: ik
 integer, intent(out) :: ik_ibz_g, tsign_k
 logical, intent(out) :: isirr
!arrays
 integer, intent(out) :: invsym_k(3,3)
 integer, intent(out) :: g0_k(3)

!Local variables-------------------------------
!scalars
 type(crystal_t) :: cryst_fake
 integer :: ik_ibz, isym_k, trev_k
 logical :: isirr_k
!arrays
 integer :: sym_k(3,3)
 real(dp) :: kpt_ibz(3)

! *************************************************************************

 cryst_fake = self%gstore%cryst%new_inversion_only()

 ik_ibz = self%map_k_khalf(1, ik); isym_k = self%map_k_khalf(2, ik)
 trev_k = self%map_k_khalf(6, ik); g0_k = self%map_k_khalf(3:5, ik)
 isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
 tsign_k = 1; if (trev_k == 1) tsign_k = -1

 kpt_ibz(:) = self%kpts_half(:, ik_ibz)
 ik_ibz_g = self%krank_kpts_g%get_index(kpt_ibz)

 sym_k(:,:) = transpose(cryst_fake%symrel(:,:,isym_k))
 call mati3inv(sym_k, invsym_k)
 invsym_k(:,:) = transpose(invsym_k)

 isirr = isirr_k
end subroutine varpeq_get_invsym_khalf
!!***


!!****f* m_varpeq/varpeq_get_invsym_k
!! NAME
!!  vapreq_get_invsym_k
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_get_invsym_k(self, ik, ik_ibz_g, tsign_k, invsym_k, g0_k, isirr)

!Arguments ------------------------------------
!scalars
 class(varpeq_t), intent(inout) :: self
 integer, intent(in) :: ik
 integer, intent(out) :: ik_ibz_g, tsign_k
 logical, intent(out) :: isirr
!arrays
 integer, intent(out) :: invsym_k(3,3)
 integer, intent(out) :: g0_k(3)

!Local variables-------------------------------
!scalars
 integer :: ik_ibz, isym_k, trev_k
 logical :: isirr_k
!arrays
 integer :: sym_k(3,3)
 real(dp) :: kpt_ibz(3)

! *************************************************************************

 ik_ibz = self%map_k_kibz(1, ik); isym_k = self%map_k_kibz(2, ik)
 trev_k = self%map_k_kibz(6, ik); g0_k = self%map_k_kibz(3:5, ik)
 isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
 tsign_k = 1; if (trev_k == 1) tsign_k = -1

 kpt_ibz(:) = self%gstore%ebands%kptns(:, ik_ibz)
 ik_ibz_g = self%krank_kpts_g%get_index(kpt_ibz)

 sym_k(:,:) = transpose(self%gstore%cryst%symrel(:,:,isym_k))
 call mati3inv(sym_k, invsym_k)
 invsym_k(:,:) = transpose(invsym_k)

 isirr = isirr_k
end subroutine varpeq_get_invsym_k
!!***


!!****f* m_varpeq/varpeq_get_invsym_q
!! NAME
!!  vapreq_get_invsym_q
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_get_invsym_q(self, iq, iq_ibz_g, tsign_q, invsym_q, g0_q, isirr)

!Arguments ------------------------------------
!scalars
 class(varpeq_t), intent(inout) :: self
 integer, intent(in) :: iq
 integer, intent(out) :: iq_ibz_g, tsign_q
 logical, intent(out) :: isirr
!arrays
 integer, intent(out) :: invsym_q(3,3)
 integer, intent(out) :: g0_q(3)

!Local variables-------------------------------
!scalars
 integer :: iq_ibz, isym_q, trev_q
 logical :: isirr_q
!arrays
 integer :: sym_q(3,3)
 real(dp) :: qpt_ibz(3)

! *************************************************************************

 iq_ibz = self%map_q_qibz(1, iq); isym_q = self%map_q_qibz(2, iq)
 trev_q = self%map_q_qibz(6, iq); g0_q = self%map_q_qibz(3:5, iq)
 isirr_q = (isym_q == 1 .and. trev_q == 0 .and. all(g0_q == 0))
 tsign_q = 1; if (trev_q == 1) tsign_q = -1

 qpt_ibz(:) = self%gstore%qibz(:, iq_ibz)
 iq_ibz_g = self%krank_qpts%get_index(qpt_ibz)

 sym_q(:,:) = self%gstore%cryst%symrec(:,:,isym_q)
 call mati3inv(sym_q, invsym_q)
 invsym_q(:,:) = transpose(invsym_q)

 isirr = isirr_q
end subroutine varpeq_get_invsym_q
!!***


!!****f* m_varpeq/varpeq_get_invsym_qhalf
!! NAME
!!  vapreq_get_invsym_qhalf
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_get_invsym_qhalf(self, iq, iq_ibz_g, tsign_q, invsym_q, g0_q, isirr)

!Arguments ------------------------------------
!scalars
 class(varpeq_t), intent(inout) :: self
 integer, intent(in) :: iq
 integer, intent(out) :: iq_ibz_g, tsign_q
 logical, intent(out) :: isirr
!arrays
 integer, intent(out) :: invsym_q(3,3)
 integer, intent(out) :: g0_q(3)

!Local variables-------------------------------
!scalars
 type(crystal_t) :: cryst_fake
 integer :: iq_ibz, isym_q, trev_q
 logical :: isirr_q
!arrays
 integer :: sym_q(3,3)
 real(dp) :: qpt_ibz(3)

! *************************************************************************

 cryst_fake = self%gstore%cryst%new_inversion_only()

 iq_ibz = self%map_q_qhalf(1, iq); isym_q = self%map_q_qhalf(2, iq)
 trev_q = self%map_q_qhalf(6, iq); g0_q = self%map_q_qhalf(3:5, iq)
 isirr_q = (isym_q == 1 .and. trev_q == 0 .and. all(g0_q == 0))
 tsign_q = 1; if (trev_q == 1) tsign_q = -1

 qpt_ibz(:) = self%qpts_half(:, iq_ibz)
 iq_ibz_g = self%krank_qpts%get_index(qpt_ibz)

 sym_q(:,:) = cryst_fake%symrec(:,:,isym_q)
 call mati3inv(sym_q, invsym_q)
 invsym_q(:,:) = transpose(invsym_q)

 isirr = isirr_q
end subroutine varpeq_get_invsym_qhalf
!!***


!!****f* m_varpeq/norm_a
!! NAME
!!  norm_a
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine norm_a(a, nsppol, nkpt, nb, nkbz, wtk)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: nsppol
 integer, intent(in) :: nkpt
 integer, intent(in) :: nb
 integer, intent(in) :: nkbz
!arrays
 real(dp), intent(in) :: wtk(:)
 complex(dp), intent(out) :: a(:, :, :)

!Local variables-------------------------------
 integer :: ik, ib, is
 real(dp) :: anorm2

! *************************************************************************

 anorm2 = zero
 do is=1,nsppol
   do ik=1,nkpt
     do ib=1,nb
       anorm2 = anorm2 + wtk(ik)*abs(a(ib, ik, is))**2
     enddo
   enddo
 enddo
 a(:, :, :) = a(:, :, :) * sqrt(one/anorm2)

end subroutine norm_a
!!***


!!****f* m_varpeq/get_deginds
!! NAME
!!  get_deginds
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine get_deginds(ib, eig, nb, deginds, ndeg)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ib
 integer, intent(in) :: nb
 integer, intent(out) :: ndeg
!arrays
 integer, allocatable :: deginds(:)
 real(dp), intent(in) :: eig(:)

!Local variables-------------------------------
!scalars
 integer :: ii, jj
 real(dp) :: eval

! *************************************************************************

 ndeg = 0
 eval = eig(ib)
 do ii=1,nb
   if (abs(eig(ii) - eval) < tol8) ndeg = ndeg + 1
 enddo

 ABI_MALLOC(deginds, (ndeg))
 jj = 1
 do ii=1,nb
   if (abs(eig(ii) - eval) < tol6) then
     deginds(jj) = ii
     jj = jj + 1
   endif
 enddo

end subroutine get_deginds
!!***

!!****f* m_varpeq/varpeq_get_gavg
!! NAME
!! vapreq_get_gavg
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

complex(dp) function varpeq_get_gavg(self, g_gathered, gqk, &
    pert, jb, ib, ik_ibz, iq, iq_glob, ikq_ibz) result(gavg)

!Arguments ------------------------------------
!scalars
 class(varpeq_t),intent(in) :: self
 class(gqk_t),intent(in) :: gqk
 integer,intent(in) :: pert
 integer,intent(in) :: jb
 integer,intent(in) :: ib
 integer,intent(in) :: ik_ibz
 integer,intent(in) :: iq
 integer,intent(in) :: iq_glob
 integer,intent(in) :: ikq_ibz
!arrays
 complex(dp),intent(in) :: g_gathered(:,:,:,:)

!Local variables-------------------------------
!scalars
 class(ebands_t), pointer :: ebands
 integer :: ii, jj, pp
 integer :: bstart
 integer :: ib_ndeg, jb_ndeg, pert_ndeg
 complex(dp) :: gavg_tmp
!arrays
 integer, allocatable :: ib_deginds(:), jb_deginds(:), pert_deginds(:)
 real(dp), allocatable :: ik_eig(:)
 real(dp), allocatable :: ikq_eig(:)
 real(dp), allocatable :: pert_freq(:)

!----------------------------------------------------------------------

 ebands => self%gstore%ebands
 bstart = self%gstore%brange_spin(1, 1)

 ABI_MALLOC(ik_eig, (self%nb))
 ABI_MALLOC(ikq_eig, (self%nb))
 ABI_MALLOC(pert_freq, (gqk%my_npert))

 ik_eig(:) = -ebands%eig(bstart:bstart+self%nb-1, ik_ibz, 1)
 ikq_eig(:) = -ebands%eig(bstart:bstart+self%nb-1, ikq_ibz, 1)
 pert_freq(:) = gqk%my_wnuq(:, iq)

 call get_deginds(ib, ik_eig, self%nb, ib_deginds, ib_ndeg)
 call get_deginds(jb, ikq_eig, self%nb, jb_deginds, jb_ndeg)
 call get_deginds(pert, pert_freq, gqk%my_npert, pert_deginds, pert_ndeg)

 gavg_tmp = zero
 do ii=1,ib_ndeg
   do jj=1,jb_ndeg
     do pp=1,pert_ndeg
       gavg_tmp = gavg_tmp + g_gathered(pert_deginds(pp), jb_deginds(jj), &
         ib_deginds(ii), iq_glob)
     enddo
   enddo
 enddo
 gavg_tmp = gavg_tmp / (ib_ndeg*jb_ndeg*pert_ndeg)
 gavg = gavg_tmp

 ABI_FREE(ib_deginds)
 ABI_FREE(jb_deginds)
 ABI_FREE(pert_deginds)

 ABI_FREE(ik_eig)
 ABI_FREE(ikq_eig)
 ABI_FREE(pert_freq)

end function varpeq_get_gavg
!!***

end module m_varpeq
!!***
