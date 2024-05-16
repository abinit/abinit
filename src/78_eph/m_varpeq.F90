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

  real(dp), allocatable :: kpts(:, :)
   ! k-points defining the electronic coefficients A_nk, possibly reduced by
   ! inversion symmetry: A_nk = A_n{-k} provided the crystal's space group has
   ! spatial inversion. No symmetries are assumed otherwise.
   ! (3, nkpt)

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

  complex(dp), allocatable :: bfull(:, :)


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

  integer, allocatable :: symafm(:)
   ! (Anti)Ferromagnetic symmetries +1/-1, helper array required for mappings
   ! (nsym)

  integer, allocatable :: symrec(:, :, :)
   ! Effective symmetry operations in reciprocal space (reduced coordinates)
   ! (3, 3, nsym)

  type(krank_t) :: krank_kpts
   ! krank datatype used to obtain k+q->k' & k-q->k'' mappings

  type(krank_t) :: krank_qpts
   ! krank datatype used to obtain q-indices

  integer, allocatable :: map_k_kibz(:, :)
  integer, allocatable :: map_invk_kibz(:, :)
   ! mappings: k-points from the electronic subspace of the problem -> k-points
   ! from IBZ; second array maps an inverse set of k-points
   ! (6, nkpt)


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
 integer :: ik, ik_ibz, imk_ibz, ik_forw, ik_back
 integer :: isym_k, trev_k, tsign_k, g0_k(3)
 integer :: my_iq, iq_glob, my_is
 integer :: my_nu
 integer :: ii
 integer :: is, ib, jb
 integer :: bstart
 complex(dp) :: a_from, a_forw, a_back
 real(dp) :: wtq, dksqmax
 real(dp) :: anorm2
 real(dp) :: enel, eig
 real(dp) :: phfreq, envib
 real(dp) :: kweight
 real(dp) :: enelph
 real(dp) :: gradres
 real(dp) :: costh, sinth
 real(dp) :: bden
 complex(dp) :: bnum, beta
 complex(dp) :: g_forw, g_back
 complex(dp) :: tmp
 complex(dp) :: bqnu
 complex(dp) :: enelph_tmp
 complex(dp) :: inprod
 logical :: isirr_k
 class(crystal_t), pointer :: cryst
 class(ebands_t), pointer :: ebands
 class(gqk_t), pointer :: gqk
!arrays
 integer, allocatable :: k_plus_q_map(:, :), k_minus_q_map(:, :)
 real(dp) :: qpt(3), my_qpt(3)
 real(dp) :: kpt_from_varpeq(3), kpt_from_gstore(3)
 real(dp) :: qpt_from_varpeq(3), qpt_from_gstore(3)
 real(dp), allocatable :: pc(:, :, :)
 complex(dp), allocatable :: prev_a(:, :, :), prev_b(:, :), prev_grad(:, :, :)


 integer :: nkbz
 integer :: ekptopt
 integer :: enkpt
 integer :: etimrev
 integer :: ik_ebz, ik_sym
 type(krank_t) :: krank_kebz
 integer :: nsym
 integer :: symrec(3,3,2)
 integer :: symafm(2)

 integer,allocatable :: k2kebz(:,:)
 real(dp) :: kpt(3), kpt_ebz(3)
 real(dp),allocatable :: ekpts(:,:), kbz(:,:)
 real(dp),allocatable :: ewtk(:)
 complex(dp),allocatable :: aeff(:,:,:)


! *************************************************************************

 cryst => self%gstore%cryst
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
 integer :: ik_ibz, imk_ibz
 integer :: ik_forw, ik_back
 integer :: my_iq, my_nu, iq_glob, my_is
 complex(dp) :: a_from, a_forw, a_back
 complex(dp) :: d_from, d_forw, d_back
 real(dp) :: kweight, phfreq
 real(dp) :: dksqmax
 real(dp) :: eig
 complex(dp) :: bqnu, g_forw, g_back
 complex(dp) :: term3_tmp, term4_tmp
 class(crystal_t), pointer :: cryst
 class(ebands_t), pointer :: ebands
 class(gqk_t), pointer :: gqk
!arrays
 integer, allocatable :: k_plus_q_map(:, :), k_minus_q_map(:, :)
 real(dp) :: my_qpt(3)
 real(dp) :: terms(4)

! *************************************************************************

 comm = self%gstore%comm
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 gqk => self%gstore%gqk(1)
 cryst => self%gstore%cryst
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
 ABI_MALLOC(k_plus_q_map, (6, self%nkpt))
 ABI_MALLOC(k_minus_q_map, (6, self%nkpt))

 do is=1,self%nsppol
   my_is = self%gstore%spin2my_is(is); if (my_is == 0) cycle
   gqk => self%gstore%gqk(my_is)

   ! Loop over q-points treated by this processor
   do my_iq=1,gqk%my_nq
     iq_glob = my_iq + gqk%my_qstart - 1
     my_qpt(:) = self%qpts(:, iq_glob)

     ! Forward (k+q) band backward (k-q) scattering mappings
     call self%krank_kpts%get_mapping(self%nkpt, self%kpts, dksqmax, &
       cryst%gmet, k_plus_q_map, self%nsym, self%symafm, self%symrec, &
       self%timrev, use_symrec=.True., qpt=my_qpt)

     ierr = merge(1, 0, dksqmax > tol12)
     if (ierr /= 0) then
       ABI_ERROR(sjoin("VarPEq, setting B-coefficients: cannot map k+q or k-q &
         to the effective IBZ with qpt ", ktoa(my_qpt)))
     endif

     ! Loop over k-points in the effective BZ
     do ik=1,self%nkpt
       kweight = self%wtk(ik)
       ik_forw = k_plus_q_map(1, ik)

       ! Loop over bands
       do ib=1,self%nb
         a_from = self%a(ib, ik, is)
         d_from = self%grad_a(ib, ik, is)

         do jb=1,self%nb
           a_forw = self%a(jb, ik_forw, is)
           d_forw = self%grad_a(jb, ik_forw, is)

           ! Loop over perturbations
           do my_nu=1,gqk%my_npert
             bqnu = self%b(my_nu, my_iq)
             !bqnu = self%bfull(my_nu, iq_glob)
             g_forw = gqk%my_g(my_nu, jb, my_iq, ib, ik)
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
   enddo
 enddo

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
 integer :: is, ik, ib, jb
 integer :: ik_ibz, ikmq_ibz
 integer :: ik_forw, ik_back
 integer :: my_iq, my_nu, iq_glob, my_is
 real(dp) :: eig
 real(dp) :: dksqmax
 real(dp) :: nkbzinv, nkbz2inv
 complex(dp) grad_tmp
 complex(dp) :: a_forw, a_back
 complex(dp) :: bqnu, g_forw, g_back
 class(crystal_t), pointer :: cryst
 class(ebands_t), pointer :: ebands
 class(gqk_t), pointer :: gqk
!arrays
 integer, allocatable :: k_plus_q_map(:, :), k_minus_q_map(:, :)
 real(dp) :: my_qpt(3)

! *************************************************************************

 comm = self%gstore%comm
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 gqk => self%gstore%gqk(1)
 cryst => self%gstore%cryst
 ebands => self%gstore%ebands

 ABI_MALLOC(k_plus_q_map, (6, self%nkpt))
 ABI_MALLOC(k_minus_q_map, (6, self%nkpt))

 nkbzinv = one/(one*self%nkbz)
 nkbz2inv = nkbzinv**2

 self%grad_a(:, :, :) = zero

 do is=1,self%nsppol
   my_is = self%gstore%spin2my_is(is); if (my_is == 0) cycle
   gqk => self%gstore%gqk(my_is)

   do my_iq=1,gqk%my_nq
     iq_glob = my_iq + gqk%my_qstart - 1
     my_qpt(:) = self%qpts(:, iq_glob)

     ! Forward (k+q) band backward (k-q) scattering mappings
     call self%krank_kpts%get_mapping(self%nkpt, self%kpts, dksqmax, &
       cryst%gmet, k_plus_q_map, self%nsym, self%symafm, self%symrec, &
       self%timrev, use_symrec=.True., qpt=my_qpt)

     call self%krank_kpts%get_mapping(self%nkpt, self%kpts, dksqmax, &
       cryst%gmet, k_minus_q_map, self%nsym, self%symafm, self%symrec, &
       self%timrev, use_symrec=.True., qpt=-my_qpt)

     ierr = merge(1, 0, dksqmax > tol12)
     if (ierr /= 0) then
       ABI_ERROR(sjoin("VarPEq, setting gradient: cannot map k+q or k-q &
         to the effective IBZ with qpt ", ktoa(my_qpt)))
     endif

     do ik=1,self%nkpt
       ik_forw = k_plus_q_map(1, ik)
       ik_back = k_minus_q_map(1, ik)

       do ib=1,self%nb
         do jb=1,self%nb
           a_forw = self%a(jb, ik_forw, is)
           a_back = self%a(jb, ik_back, is)

           do my_nu=1,gqk%my_npert
             bqnu = self%b(my_nu, my_iq)
             !bqnu = self%bfull(my_nu, iq_glob)
             g_forw = gqk%my_g(my_nu, jb, my_iq, ib, ik)
             g_back = gqk%my_g(my_nu, ib, my_iq, jb, ik_back)

             self%grad_a(ib, ik, is) = self%grad_a(ib, ik, is) - &
               two*nkbz2inv*(a_back*conjg(bqnu)*g_back + a_forw*bqnu*conjg(g_forw))
           enddo
         enddo
       enddo
     enddo
   enddo
 enddo
 call xmpi_sum(self%grad_a, comm, ierr)

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

 ABI_FREE(k_plus_q_map)
 ABI_FREE(k_minus_q_map)

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
 integer :: ik_ibz, imk_ibz
 integer :: ik_forw, ik_back
 integer :: my_iq, my_nu, iq_glob, my_is
 complex(dp) :: a_from, a_forw, a_back
 real(dp) :: enel, envib, enelph
 real(dp) :: kweight, phfreq
 real(dp) :: dksqmax
 real(dp) :: eig
 complex(dp) :: enelph_tmp
 complex(dp) :: bqnu, g_forw, g_back
 class(crystal_t), pointer :: cryst
 class(ebands_t), pointer :: ebands
 class(gqk_t), pointer :: gqk
!arrays
 integer, allocatable :: k_plus_q_map(:, :), k_minus_q_map(:, :)
 real(dp) :: my_qpt(3)

! *************************************************************************

 comm = self%gstore%comm
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 gqk => self%gstore%gqk(1)
 cryst => self%gstore%cryst
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
   do my_nu=1,gqk%my_npert
     phfreq = gqk%my_wnuq(my_nu, my_iq)
     envib = envib + phfreq*abs(self%b(my_nu, my_iq))**2
     !envib = envib + phfreq*abs(self%bfull(my_nu, iq_glob))**2
   enddo
 enddo
 call xmpi_sum(envib, comm, ierr)
 envib = envib/(one*self%nqpt)

 ! Electron-phonon energy
 enelph_tmp = (zero, zero)
 ABI_MALLOC(k_plus_q_map, (6, self%nkpt))
 ABI_MALLOC(k_minus_q_map, (6, self%nkpt))
 do is=1,self%nsppol
   my_is = self%gstore%spin2my_is(is); if (my_is == 0) cycle
   gqk => self%gstore%gqk(my_is)

   ! Loop over q-points treated by this processor
   do my_iq=1,gqk%my_nq
     iq_glob = my_iq + gqk%my_qstart - 1
     my_qpt(:) = self%qpts(:, iq_glob)

     ! Forward (k+q) band backward (k-q) scattering mappings
     call self%krank_kpts%get_mapping(self%nkpt, self%kpts, dksqmax, &
       cryst%gmet, k_plus_q_map, self%nsym, self%symafm, self%symrec, &
       self%timrev, use_symrec=.True., qpt=my_qpt)

     ierr = merge(1, 0, dksqmax > tol12)
     if (ierr /= 0) then
       ABI_ERROR(sjoin("VarPEq, setting B-coefficients: cannot map k+q or k-q &
         to the effective IBZ with qpt ", ktoa(my_qpt)))
     endif

     ! Loop over k-points in the effective BZ
     do ik=1,self%nkpt
       kweight = self%wtk(ik)
       ik_forw = k_plus_q_map(1, ik)

       ! Loop over bands
       do ib=1,self%nb
         a_from = self%a(ib, ik, is)

         do jb=1,self%nb
           a_forw = self%a(jb, ik_forw, is)

           ! Loop over perturbations
           do my_nu=1,gqk%my_npert
             bqnu = self%b(my_nu, my_iq)
             !bqnu = self%bfull(my_nu, iq_glob)
             g_forw = gqk%my_g(my_nu, jb, my_iq, ib, ik)
             enelph_tmp = enelph_tmp + &
               kweight*a_from*conjg(bqnu)*g_forw*conjg(a_forw)
           enddo

         enddo
       enddo
     enddo
   enddo
 enddo
 enelph_tmp = enelph_tmp + conjg(enelph_tmp)
 call xmpi_sum(enelph_tmp, comm, ierr)
 enelph = -real(enelph_tmp) / (one*self%nqpt)

 ! Setting the corresponding state variables of the varpeq type
 self%enterms(1) = enel; self%enterms(2) = envib; self%enterms(3) = enelph
 self%eps = enel + enelph

 ABI_FREE(k_plus_q_map)
 ABI_FREE(k_minus_q_map)

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
 integer :: comm, nproc, my_rank, ierr
 integer :: is, ik, ib, jb
 integer :: my_iq, my_nu, my_is
 integer :: my_nq, my_npert
 integer :: iq_glob
 integer :: ik_forw, ik_back, ik_ibz, imk_ibz, ik_ibz_glob
 integer :: isym_k, trev_k, g0_k(3), tsign_k
 integer :: iq_new
 logical :: isirr_k
 real(dp) :: my_phfreqinv, kweight
 real(dp) :: dksqmax
 complex(dp) :: b_tmp
 complex(dp) :: a_from, a_forw, a_back
 complex(dp) :: g_forw, g_back
 class(crystal_t), pointer :: cryst
 class(ebands_t), pointer :: ebands
 class(gqk_t), pointer :: gqk

!arrays
 integer, allocatable :: k_plus_q_map(:, :), k_minus_q_map(:, :)
 integer, allocatable :: q_plus_k_map(:, :)
 integer :: sym(3,3), syminv(3,3)
 real(dp) :: my_qpt(3), my_kpt(3), my_kpt_ibz(3), new_qpt(3)
 complex(dp), allocatable :: gq_gathered(:, :, :, :)
 complex(dp), allocatable :: gq_symgathered(:, :, :, :)

! *************************************************************************

 comm = self%gstore%comm
 nproc = xmpi_comm_size(comm)
 my_rank = xmpi_comm_rank(comm)

 cryst => self%gstore%cryst
 ebands => self%gstore%ebands

 ABI_MALLOC(k_plus_q_map, (6, self%nkpt))
 ABI_MALLOC(k_minus_q_map, (6, self%nkpt))
 ABI_MALLOC(q_plus_k_map, (6, self%nqpt))

 self%b(:, :) = zero
 self%bfull(:, :) = zero

 gqk => self%gstore%gqk(1)
 my_nq = gqk%my_nq; my_npert = gqk%my_npert

 ABI_MALLOC(gq_gathered, (my_npert, self%nb, self%nb, self%nkpt))
 ABI_MALLOC(gq_symgathered, (my_npert, self%nb, self%nb, self%nkpt))

 ! Loop over q-points
 do my_iq=1,my_nq
   iq_glob = my_iq + gqk%my_qstart - 1
   my_qpt = self%qpts(:, iq_glob)

   ! Forward (k+q) scattering mapping
   call self%krank_kpts%get_mapping(self%nkpt, self%kpts, dksqmax, &
     cryst%gmet, k_plus_q_map, self%nsym, self%symafm, self%symrec, &
     self%timrev, use_symrec=.True., qpt=my_qpt)

   ierr = merge(1, 0, dksqmax > tol12)
   if (ierr /= 0) then
     ABI_ERROR(sjoin("VarPEq, setting B-coeffeicient: cannot map k+q or k-q &
       to the effective IBZ with qpt ", ktoa(my_qpt)))
   endif

   do my_nu=1,my_npert
     b_tmp = (zero, zero)

     ! for a (nu,q) collect the mnk\sigma sum
     do is=1,self%nsppol
       my_is = self%gstore%spin2my_is(is); if (my_is == 0) cycle
       gqk => self%gstore%gqk(my_is)

       do ik=1,self%nkpt
         kweight = self%wtk(ik)
         ik_forw = k_plus_q_map(1, ik)

         do ib=1,self%nb
           a_from = self%a(ib, ik, is)

           do jb=1,self%nb
             a_forw = self%a(jb, ik_forw, is)
             g_forw = gqk%my_g(my_nu, jb, my_iq, ib, ik)

             b_tmp = b_tmp + conjg(a_forw)*g_forw*a_from*kweight
           enddo

         enddo
       enddo
     enddo

     if (gqk%my_wnuq(my_nu, my_iq) == 0) then
       self%b(my_nu, my_iq) = zero
       !self%bfull(my_nu, iq_glob) = zero
     else
       self%b(my_nu, my_iq) = b_tmp/gqk%my_wnuq(my_nu, my_iq)
       !self%bfull(my_nu, iq_glob) = b_tmp/gqk%my_wnuq(my_nu, my_iq)
     endif

   enddo
 enddo

 ABI_FREE(k_plus_q_map)
 ABI_FREE(k_minus_q_map)
 ABI_FREE(q_plus_k_map)
 ABI_FREE(gq_gathered)
 ABI_FREE(gq_symgathered)

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
       !imeig = -j_dpc*ebands%eig(bstart+ib-1, ik_ibz, is)
       !imeig = j_dpc*ebands%eig(bstart+ib-1, ik_ibz, is)
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
 integer :: iinv
 real(dp) :: wtq
 class(crystal_t), pointer :: cryst
 class(ebands_t), pointer :: ebands
 class(gqk_t), pointer :: gqk
 type(krank_t) :: krank_ibz
!arrays
 real(dp) :: my_qpt(3)
 real(dp), allocatable :: kbz(:, :)

! *************************************************************************

 comm = gstore%comm
 nproc = xmpi_comm_size(comm)
 my_rank = xmpi_comm_rank(comm)

 ! Consistency check
 ierr = 0
 ABI_CHECK_NOSTOP(gstore%kzone == "bz", "kzone = 'bz' is required", ierr)
 ABI_CHECK_NOSTOP(gstore%qzone == "bz", "qzone = 'bz' is required", ierr)
 ABI_CHECK_NOSTOP(gstore%gqk(1)%cplex == 2, "cplex = 2 is required", ierr)
 ABI_CHECK_NOSTOP(gstore%gqk(1)%kpt_comm%nproc == 1, "varpeq is not &
   compatible with parallelism over k-points", ierr)
 ABI_CHECK(ierr == 0, "The gstore object is incosistent with varpeq. &
   See messages above.")
 self%gstore => gstore

 ! Initalization of the electronic subspace dimensions
 ! and related symmetry tables
 cryst => gstore%cryst
 ebands => gstore%ebands
 self%nsppol = gstore%nsppol
 self%nb = gstore%gqk(1)%nb

 ! Time-reversal symmetry -> generate k-points using only spatial inversion
 ! no time-reveral symmetry -> k-points in the full BZ
 !! self%timrev = kpts_timrev_from_kptopt(ebands%kptopt)
 !! iinv = cryst%idx_spatial_inversion()
 ! hack to artificially disable any symmetry except for identity
 self%timrev = 0
 if (self%timrev /= 0 .and. iinv /= 0)  then
   kptopt = 2; self%nsym = 2
 else
   kptopt = 3; self%nsym = 1
 endif

 ABI_MALLOC(self%symrec, (3, 3, self%nsym))
 ABI_MALLOC(self%symafm, (self%nsym))
 self%symafm(:) = 1
 self%symrec(:, :, :) = zero
 ! Identity
 self%symrec(1, 1, 1) = 1; self%symrec(2, 2, 1) = 1; self%symrec(3, 3, 1) = 1
 ! Inversion symmetry if present
 if (self%nsym == 2) then
   self%symafm(2) = cryst%symafm(iinv)
   self%symrec(:, :, 2) = cryst%symrec(:, :, iinv)
 end if

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

 ! Initialization of the vibrational subspace dimensions
 ! generate q-points in the full BZ
 gqk => self%gstore%gqk(1)
 self%nqpt = gstore%nqbz
 ABI_MALLOC(self%b, (gstore%gqk(1)%my_npert, gstore%gqk(1)%my_nq))
 ABI_MALLOC(self%qpts, (3, self%nqpt))
 ABI_MALLOC(self%bfull, (gstore%gqk(1)%my_npert, self%nqpt))
 self%qpts(:, :) = zero

 do my_iq=1,gqk%my_nq
   iq_glob = my_iq + gqk%my_qstart - 1
   call gqk%myqpt(my_iq, gstore, wtq, my_qpt)
   self%qpts(:, iq_glob) = my_qpt(:)
 enddo
 call xmpi_sum(self%qpts, comm, ierr)

 ! krank for qpts
 self%krank_qpts = krank_from_kptrlatt(self%nqpt, self%qpts, ebands%kptrlatt, &
   compute_invrank = .True.)

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

 ! Symmetry tables
 ABI_SFREE(self%symafm)
 ABI_SFREE(self%symrec)
 ABI_SFREE(self%map_k_kibz)
 ABI_SFREE(self%map_invk_kibz)
 call self%krank_kpts%free()

end subroutine varpeq_free
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

end module m_varpeq
!!***
