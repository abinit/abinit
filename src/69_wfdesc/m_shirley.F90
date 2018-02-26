!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_shirley
!! NAME
!!  m_shirley
!!
!! FUNCTION
!!  Procedures and objects for the interpolation of
!!  KS eigenvalues and eigenvectors with Shirley's method.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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

MODULE m_shirley

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_xmpi
 use m_profiling_abi
 use m_fft
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_fstrings,       only : firstchar
 use m_time,           only : cwtime
 use m_io_tools,       only : flush_unit
 use m_numeric_tools,  only : print_arr, hermitianize, imax_loc, stats_t, stats_eval, wrap2_pmhalf
 use m_blas,           only : xdotc, xgemm, blas_cholesky_ortho
 use m_abilasi,        only : xheev, xhegv, xheevx, xhegvx
 use m_fft_mesh,       only : fft_check_rotrans, setmesh
 use m_geometry,       only : normv, vdotw
 use m_fftcore,        only : get_kg
 use m_crystal,        only : crystal_t
 use m_crystal_io,     only : crystal_ncwrite
 use m_bz_mesh,        only : kmesh_t, get_BZ_item, make_mesh, kmesh_free, kpath_t
 use m_ebands,         only : pack_eneocc, ebands_init, ebands_ncwrite, ebands_free
 use m_hamiltonian,    only : ddiago_ctl_type, init_ddiago_ctl, destroy_hamiltonian, init_hamiltonian, &
&                             load_spin_hamiltonian,load_k_hamiltonian, gs_hamiltonian_type
 use m_wfd,            only : wfd_get_ur, wfd_t, wfd_barrier, wfd_get_cprj, wfd_barrier, wfd_change_ngfft, wfd_print,&
&                             wfd_sym_ur, wfd_init, wfd_push_ug, wfd_reset_ur_cprj, wfd_ihave_ur, wfd_free,&
&                             wfd_test_ortho, kdata_t, kdata_init, kdata_free, wfd_set_mpicomm, wfd_itreat_spin, wfd_xdotc
 use m_pawang,         only : pawang_type
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type
 use m_paw_ij,         only : paw_ij_type
 use m_pawrhoij,       only : pawrhoij_type
 use m_pawcprj,        only : pawcprj_type, pawcprj_alloc, pawcprj_free, paw_overlap
 use m_pawpwij,       only : pawpwff_t, pawpwff_init, pawpwff_free, pawpwij_t, pawpwij_init,&
&                             pawpwij_free, paw_rho_tw_g
 use m_pawfgr,         only : pawfgr_type

 implicit none

 private

 public :: wfd_bloch_to_shirley ! Produce new wavefunction descriptor with the optical basis set from a set of KS states,
 public :: shirley_interp       ! Interpolate the KS Hamiltonian in the Shirley basis set.
 public :: shirley_hks          ! Compute <i|H_ks|j> in the Shirley basis set for a single k-point and spin
 public :: shirley_window       ! Find an optimal basis set for the interpolation of the KS states in a energy window,
 public :: shirley_bands        ! Compute interpolated bands (no eigenvectors)
 !public :: wfd_shirley_to_eh   ! Under development
!!***

!----------------------------------------------------------------------

!!****t* m_shirley/shparams_t
!! NAME
!!  shparams_t
!!
!! FUNCTION
!!  Store the parameters controlling the algorith for the
!!  Shirley interpolation of the KS band energies and wavefunctions.
!!
!! SOURCE

 type,public :: shparams_t

   integer :: ndivsm=20
     ! Number of division for the smallest segment of the k-path.
     ! Used for the interpolation of band structures

     ! Options used for the generazion of the K-mesh if sampling == "mesh"
     ! See abinit documentation
   integer :: kptopt = 3
   integer :: kptrlatt(3,3) = -1
   integer,allocatable :: shiftk(:,:)

   real(dp) :: ewin(2) = [smallest_real, greatest_real]
     ! Energy window (default: unused)

   real(dp) :: dos_tol = tol6
     ! Tolerance on the converged of the DOS (used if method=="iterative")

   character(len=500) :: method="one-shot"
     ! Method used for selecting an optimal basis set inside the energy window
     ! "one-shot"
     ! "iterative"

   character(len=500) :: sampling="path"
     ! Type of sampling
     ! "path"
     ! "mesh"

   character(len=500) :: jobz = "Normal"
     ! jobz="V" if eigevectors are wanted. "N" if only eigenvalues are required.

   real(dp),allocatable :: kptbounds(:,:)
     ! List of vertex for band structure interpolation.

 end type shparams_t

 !public :: shparams_init
!!***

!----------------------------------------------------------------------

! Flags used to define the content of ovlp%mat.
 integer,private,parameter :: TYPE_NONE    = 0
 integer,private,parameter :: TYPE_OVERLAP = 1
 integer,private,parameter :: TYPE_EIGVEC  = 2

 logical,private,parameter :: DEBUGME=.False.

!----------------------------------------------------------------------

!!****t* m_shirley/shres_t
!! NAME
!!  shres_t
!!
!! FUNCTION
!! The results of the Shirley interpolation for a given k-point and spin.
!! TODO fix issue with the indexing as we might want to read a subset of bands
!! for ouri nterpolation moreover one might used vl and vu in the direct diagonalization
!! to select the energy window we are interested in.
!!
!! SOURCE

 type,public :: shres_t

   integer :: sh_size
   ! Number of Shirley optimal basis set elements.

   integer :: nband_k
   ! Number of interpolated Kohn-Sham bands.

   real(dp),allocatable :: ene(:)
   ! ene(nband_k)
   ! Interpolated energies

   complex(dpc),allocatable :: obloch(:,:)
   ! obloch(sh_size,nband_k)
   ! Matrix with the coeffients <O_i|Bloch_b>

 end type shres_t

 public :: shres_init
 public :: shres_free

 interface shres_free
   module procedure shres_free_0D
   module procedure shres_free_2D
 end interface shres_free
!!***

!----------------------------------------------------------------------

!!****t* m_shirley/ovlp_t
!! NAME
!!  ovlp_t
!!
!! FUNCTION
!!  Store the overlap matrix <u_bk|u_{b'k'}> for a given spin.
!!
!! SOURCE

 type,public :: ovlp_t

  integer :: size=0
  ! The size of the overlap matrix.

  integer :: mband=0
  ! Maximum number of bands.

  integer :: nkpt=0
  ! Number of k-points.

  integer :: mat_type = TYPE_NONE
  ! Matrix type (overlap or eigenvectors)

  real(dp) :: min_ene = -HUGE(one)
  ! Min energy included.

  real(dp) :: max_ene = +HUGE(one)
  ! Max energy included

  integer(i8b),allocatable :: bk2idx(:,:)
  ! bk2idx(mband,nkpt)
  ! Mapping (b,k) --> i
  ! 0 if the state has been excluded from the energy window.

  integer,allocatable :: idx2bk(:,:)
  ! idx2bk(2,size))
  ! Mapping i --> (b,k) for i=1,size. k is always the in the BZ.

  complex(dpc),allocatable :: mat(:,:)
  ! mat(size,size)
  ! Stores O_{ij} = <u_i|u_j>
  ! The matrix is Hermitian hence one could use a packed matrix to save memory,
  ! but the Lapack routines for the diagonalization are usually slower.

  real(dp),allocatable :: eigene(:)
  ! eigene(size)
  ! The eigenvalues of the overlap operator.

 end type ovlp_t

 public :: ovlp_init               ! Creation method
 public :: ovlp_free               ! Deallocate memory
 public :: ovlp_diago_and_prune    ! Diagonalize the overlap matrix and select the important subspace.
!!***

CONTAINS  !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_shirley/shres_init
!! NAME
!!  shres_init
!!
!! FUNCTION
!!
!! PARENTS
!!
!! CHILDREN
!!      blas_cholesky_ortho,cwtime,fft_ur,fftbox_execute,fftbox_plan3_many
!!      flush_unit,get_kg,kgindex,ovlp_diago_and_prune,ovlp_free
!!      wfd_change_ngfft,wfd_get_ur,wfd_init,wfd_print,wfd_push_ug
!!      wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine shres_init(Shres,sh_size,nband_k,intp_ene,obloch)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'shres_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: sh_size,nband_k
 type(shres_t),intent(inout) :: Shres
!arrays
 real(dp),intent(in) :: intp_ene(nband_k)
 complex(dpc),intent(in) :: obloch(sh_size,nband_k)

! *************************************************************************

 !@shres_t
 Shres%sh_size = sh_size
 Shres%nband_k = nband_k

 ABI_MALLOC(Shres%ene,(nband_k))
 Shres%ene = intp_ene

 ABI_MALLOC(Shres%obloch,(sh_size,nband_k))
 Shres%obloch = obloch

end subroutine shres_init
!!***

!----------------------------------------------------------------------

!!****f* m_shirley/shres_free_0D
!! NAME
!!  shres_free_0D
!!
!! FUNCTION
!!  Free memory.
!!
!! PARENTS
!!      m_shirley
!!
!! CHILDREN
!!      blas_cholesky_ortho,cwtime,fft_ur,fftbox_execute,fftbox_plan3_many
!!      flush_unit,get_kg,kgindex,ovlp_diago_and_prune,ovlp_free
!!      wfd_change_ngfft,wfd_get_ur,wfd_init,wfd_print,wfd_push_ug
!!      wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine shres_free_0D(Shres)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'shres_free_0D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(shres_t),intent(inout) :: Shres

! *************************************************************************

 !@shres_t
 ! real
 if (allocated(Shres%ene)) then
   ABI_FREE(Shres%ene)
 end if
 !
 ! complex
 if (allocated(Shres%obloch)) then
   ABI_FREE(Shres%obloch)
 end if

end subroutine shres_free_0D
!!***

!----------------------------------------------------------------------

!!****f* m_shirley/shres_free_2D
!! NAME
!!  shres_free_2D
!!
!! FUNCTION
!!  Free memory
!!
!! PARENTS
!!
!! CHILDREN
!!      blas_cholesky_ortho,cwtime,fft_ur,fftbox_execute,fftbox_plan3_many
!!      flush_unit,get_kg,kgindex,ovlp_diago_and_prune,ovlp_free
!!      wfd_change_ngfft,wfd_get_ur,wfd_init,wfd_print,wfd_push_ug
!!      wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine shres_free_2D(Shres)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'shres_free_2D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(shres_t),intent(inout) :: Shres(:,:)

!Local variables
 integer :: ii,jj

! *************************************************************************

 do jj=1,size(Shres,dim=2)
   do ii=1,size(Shres,dim=1)
     call shres_free_0d(Shres(ii,jj))
   end do
 end do

end subroutine shres_free_2D
!!***

!----------------------------------------------------------------------

!!****f* m_shirley/ovlp_free
!! NAME
!!  ovlp_free
!!
!! FUNCTION
!!  Free the memory.
!!
!! PARENTS
!!      m_shirley
!!
!! CHILDREN
!!      blas_cholesky_ortho,cwtime,fft_ur,fftbox_execute,fftbox_plan3_many
!!      flush_unit,get_kg,kgindex,ovlp_diago_and_prune,ovlp_free
!!      wfd_change_ngfft,wfd_get_ur,wfd_init,wfd_print,wfd_push_ug
!!      wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine ovlp_free(Ovlp)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ovlp_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ovlp_t),intent(inout) :: Ovlp

! *************************************************************************

 DBG_ENTER("COLL")

 !@ovlp_t
 ! integer
 if (allocated(Ovlp%bk2idx)) then
   ABI_FREE(Ovlp%bk2idx)
 end if
 if (allocated(Ovlp%idx2bk)) then
   ABI_FREE(Ovlp%idx2bk)
 end if
 !
 ! real
 if (allocated(Ovlp%eigene)) then
   ABI_FREE(Ovlp%eigene)
 end if
 !
 ! complex
 if (allocated(Ovlp%mat)) then
   ABI_FREE(Ovlp%mat)
 end if

 DBG_EXIT("COLL")

end subroutine ovlp_free
!!***

!----------------------------------------------------------------------

!!****f* m_shirley/ovlp_init
!! NAME
!! ovlp_init
!!
!! FUNCTION
!!  Calculates the upper triangle of the overlap matrix <u_i|u_j> for a given spin and for
!!  all the possible combinations (k,b|k',b') with k and k' in the full Brillouin zone.
!!  The u functions are the periodic part of the Bloch wavefunctions hence there is
!!  no selection rule in k-space for the matrix elements when k != k'
!!  The diagonal matrix elements equals one provided that the input wavefunctions are correctly normalized.
!!
!! INPUTS
!! ewin(2,Wfd%nsppol)=Window energies for the two spins.
!! ov_ngfft(18)=FFT mesh used for the computation of the overlap in real space.
!! use_sym=logical flag defining whether symmetries have to be used for reconstructing the overlap matrix.
!! Wfd<wfd_t>=Datatype gathering info on the wavefunctions.
!! Cryst<crystal_t>= data type gathering info on symmetries and unit cell
!! Kmesh<kmesh_t>=Datatype describing the k-point sampling used for the wavefunctions.
!! Bands<ebands_t>=Band structure energies.
!! Psps <type(pseudopotential_type)>=variables related to pseudopotentials
!! Pawtab(Psps%ntypat) <type(pawtab_type)>=paw tabulated starting data
!! Pawang<pawang_type> angular mesh discretization and related data.
!! Pawrad(ntypat*usepaw)<type(pawrad_type)>=paw radial mesh and related data.
!!
!! OUTPUT
!! O(Wfd%nsppol) <type(ovlp_t)>=Object with the overlap matrix elements.
!!
!! PARENTS
!!      m_shirley
!!
!! CHILDREN
!!      blas_cholesky_ortho,cwtime,fft_ur,fftbox_execute,fftbox_plan3_many
!!      flush_unit,get_kg,kgindex,ovlp_diago_and_prune,ovlp_free
!!      wfd_change_ngfft,wfd_get_ur,wfd_init,wfd_print,wfd_push_ug
!!      wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine ovlp_init(O,ewin,use_sym,ov_ngfft,Wfd,Cryst,Kmesh,Bands,Psps,Pawtab,Pawang,Pawrad)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ovlp_init'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
 use interfaces_65_paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 logical,intent(in) :: use_sym
 type(crystal_t),intent(in) :: Cryst
 type(Pawang_type),intent(in) :: Pawang
 type(Pseudopotential_type),intent(in) :: Psps
 type(kmesh_t),intent(in) :: Kmesh
 type(ebands_t),intent(in) :: Bands
 type(wfd_t),target,intent(inout) :: Wfd
 type(ovlp_t),intent(out) :: O(Wfd%nsppol)
!arrays
 integer,intent(in) :: ov_ngfft(18)
 real(dp),intent(in) :: ewin(2,Wfd%nsppol)
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
 type(Pawrad_type),intent(in) :: Pawrad(Cryst%ntypat*Wfd%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: dim_rtwg1=1,option_lob0=0
 integer :: ik1_bz,ik2_bz,ik1_ibz,ik2_ibz,ierr,rhoxsp_method,nqpt,iqpt,rk
 integer :: band1,band2,nband_k,nband_k2,nband_k1,npw_k1,npw_k2,spin
 integer :: with_sym,without_sym
 integer(i8b) :: row,col,ovlp_size,os_size
!MPI
 integer(i8b) :: tot_nels,oidx,my_hsize,bsize_my_block
 integer(i8b) :: my_cols(2),my_rows(2)
 integer(i8b),allocatable :: t_start(:),t_stop(:),hsize_of(:)
!ENDMPI
 integer :: k1_sym,k2_sym,k1_tim,k2_tim,inv_k2_sym,inv_k1_sym
 integer :: band1_stop,nspinor,lk,rak,iq_found,comm_spin,nproc_spin,rank_spin
 integer :: nsppol,mband,nkpt !,npw_k,istwf_k,!,ii,itypat,klmn,klmn_size
 real(dp),parameter :: tnons_tol=tol8,tol_kdiff=tol6
 real(dp) :: fft_fact,ksq,qsq_max,ene_bks,cpu,wall,gflops
 complex(dpc) :: blk_ovlp,covlp,k1_eimkt,k2_eimkt,tnons_fact
 complex(gwpc) :: paw_ovlp(1)
 logical :: k1_isirred,k2_isirred,can_use_sym,take_cnj
 character(len=500) :: msg
!arrays
 integer,parameter :: g_gamma(3)=(/0,0,0/)
 integer :: k1_umklp(3),k2_umklp(3)
 integer,allocatable :: tmp_idx2bk(:,:)
 integer,allocatable :: toinv(:,:),multable(:,:,:),nq_spl(:),k1mk2q(:,:)
 real(dp) :: kpoint(3),fft_err(3,Cryst%nsym) !,onsite(2)
 real(dp) :: kk1(3),kk2(3),k1mk2(3) !,ovlp_paw(2) !,r1_tau3(3)
 real(dp),allocatable :: qpts(:,:),qmax(:)
 complex(gwpc),allocatable :: ur1(:),ur2(:)
 complex(gwpc),pointer :: ug2(:) !ug1(:)
 complex(dpc),allocatable :: ovlp_ikfk(:,:,:,:)
 logical :: k_needs_tr(2)
 type(pawcprj_type),pointer :: Cp_k1(:,:),Cp_k2(:,:)
 type(pawpwij_t),allocatable :: Pwij(:,:)
 type(pawpwff_t),allocatable :: Paw_pwff(:)

!************************************************************************

 DBG_ENTER("COLL")

 ABI_CHECK(Wfd%nspinor==1,"nspinor==2 not coded")

 ABI_UNUSED(Pawrad(1)%mesh_size)

 call cwtime(cpu,wall,gflops,"start")

 ! Change FFT mesh.
 call wfd_change_ngfft(Wfd,Cryst,Psps,ov_ngfft)

 ! Check if the mesh preserves the symmetries
 if (.not.fft_check_rotrans(Cryst%nsym,Cryst%symrel,Cryst%tnons,Wfd%ngfft,fft_err)) then
   write(msg,'(a,3(i0,1x),a)')" Real space FFT mesh ",Wfd%ngfft(1:3)," is not symmetric. Cannot symmetrize in real space"
   MSG_ERROR(msg)
 end if

 nsppol = Wfd%nsppol; nspinor = Wfd%nspinor
 mband = Wfd%mband

 ! Use the total number of k-points in the BZ.
 nkpt = Kmesh%nbz
 !
 ! Calculate the overlap matrix
 ! ======================================================
 ! ====  <u_{ik b1}| u_{fk b2}>, ik in IBZ, fk in BZ ====
 ! ======================================================

 ABI_MALLOC(multable,(4,Cryst%nsym,Cryst%nsym))
 ABI_MALLOC(toinv,(4,Cryst%nsym))
 call sg_multable(Cryst%nsym,Cryst%symafm,Cryst%symrel,Cryst%tnons,tnons_tol,ierr,multable=multable,toinv=toinv)
 ABI_CHECK(ierr==0,"Group error, cannot continue")
 !
 ! For PAW we have to recalculate the projections in the IBZ setting k=0 in the exponential.
 ! TODO: For the time being, these terms are calculated in the BZ, symmetrization will be added afterwards.
 if (Wfd%usepaw==1) then
   !
   ! <u_k2|u_k1> = <phi_k2|e^{-i(k1-k2).r}|psi_k1>
   ! Find the list of q-points: q=k1-k2 and create table (k1, k2) --> q = k1-k2
   ABI_MALLOC(qpts,(3,nkpt**2))
   ABI_MALLOC(k1mk2q,(nkpt,nkpt))

   nqpt=0; qsq_max=zero
   do ik1_bz=1,nkpt
     kk1 = Kmesh%bz(:,ik1_bz)
     do ik2_bz=1,nkpt
       kk2 = Kmesh%bz(:,ik2_bz)

       k1mk2 = kk1 - kk2
       !k1mk2 = kk2 - kk1
       ksq = normv(k1mk2,Cryst%gmet,"G")
       qsq_max = MAX(ksq,qsq_max)

       iq_found=0
       do iqpt=1,nqpt
         if (ALL( ABS(k1mk2 - qpts(:,iqpt)) < tol_kdiff) ) then
           iq_found=iqpt
           EXIT
         end if
       end do

       if (iq_found==0) then
         ! Add it to the list.
         nqpt = nqpt+1
         !ABI_CHECK(nqpt <= size(qpts,dim=2),"too many qpts!")
         qpts(:,nqpt) = k1mk2
         iq_found = nqpt
       end if

       k1mk2q(ik1_bz,ik2_bz) = iq_found
     end do
   end do

   ! TODO: This one should be correct but it's too small, why?
   !qsq_max = qsq_max/(two_pi**2)
   qsq_max = qsq_max/(two*pi**2)
   !
   ! Set up q-grid for form factors, make qmax 20% larger than the largest expected.
   ABI_MALLOC(nq_spl,(Cryst%ntypat))
   ABI_MALLOC(qmax,(Cryst%ntypat))
   nq_spl = Psps%mqgrid_ff
   qmax = SQRT(qsq_max)*1.2d0
   write(std_out,*)"Using nqpt:",nqpt," nq_spl ",nq_spl," qmax_type=",qmax

   ! TODO: pass input variable to control the computation.
   rhoxsp_method=1 ! Arnaud-Alouani
   !rhoxsp_method=2 ! Shiskin-Kresse

   ABI_DT_MALLOC(Paw_pwff,(Psps%ntypat))
   call pawpwff_init(Paw_pwff,rhoxsp_method,nq_spl,qmax,Cryst%gmet,Pawrad,Pawtab,Psps)

   ABI_FREE(nq_spl)
   ABI_FREE(qmax)

   ABI_DT_MALLOC(Pwij,(Psps%ntypat,nqpt))
   do iqpt=1,nqpt
     call pawpwij_init(Pwij(:,iqpt),1,qpts(:,iqpt),g_gamma,Cryst%rprimd,Psps,Pawtab,Paw_pwff)

     ! DEBUG
     !if (iqpt==1) then
     !  write(std_out,*)qpts(:,iqpt)
     !  do itypat=1,Cryst%ntypat
     !    klmn_size = Pwij(itypat,iqpt)%lmn_size*(Pwij(itypat,iqpt)%lmn_size+1)/2
     !    do klmn=1,klmn_size
     !      !write(std_out,*)" mqpgij: ",Pwij(itypat,iqpt)%mqpgij(:,1,klmn),pawtab(itypat)%sij(klmn)
     !      write(std_out,*)" diff mqpg0ij: ",Pwij(itypat,iqpt)%mqpgij(1,1,klmn)-pawtab(itypat)%sij(klmn)
     !      Pwij(itypat,iqpt)%mqpgij(1,1,klmn) = pawtab(itypat)%sij(klmn)
     !      Pwij(itypat,iqpt)%mqpgij(2,1,klmn) = zero
     !    end do
     !  end do
     !end if
   end do

   call pawpwff_free(Paw_pwff)
   ABI_DT_FREE(Paw_pwff)
   ABI_FREE(qpts)

   ! Allocate cprj datastructures.
   ABI_DT_MALLOC(Cp_k1 ,(Cryst%natom,nspinor))
   call pawcprj_alloc(Cp_k1,0,Wfd%nlmn_atm)

   ABI_DT_MALLOC(Cp_k2 ,(Cryst%natom,nspinor))
   call pawcprj_alloc(Cp_k2,0,Wfd%nlmn_atm)
 end if ! PAW
 !
 ! ==============================================================
 ! ==== Reconstruct full <u_{kb}| u_{k'b'}> matrix in the BZ ====
 ! ==============================================================
 ! 1) Symmetrization is done in real space. Easier, especially when k-centered G-sphere are used.
 !     u(r,b,kbz)=e^{-2i\pi kibz.(R^{-1}t} u (R{^-1}(r-t),b,kibz)
 !               =e^{+2i\pi kibz.(R^{-1}t} u*({R^-1}(r-t),b,kibz) for time-reversal
 ! 2) Matrix is Hermitian.
 ! 3) <u_{Sk b}| u_{Sk b'}> are obtained from the previously calculated <u_{kb}| u_{kb'}> table.
 !
 ! {A,a} {B,b} = {AB, a+Ab}
 !
 fft_fact = one/Wfd%nfft
 ABI_MALLOC(ur1,(Wfd%nfft*nspinor))
 ABI_MALLOC(ur2,(Wfd%nfft*nspinor))

 do spin=1,nsppol
   if (.not. wfd_itreat_spin(Wfd,spin,comm_spin,rank_spin,nproc_spin)) cycle

   ! Compute table (b,k) --> indx and indx --> (b,k) where (b,k) is
   ! the set of bands and kpoints included in the calculation of the overlap.
   O(spin)%mband   = mband
   O(spin)%nkpt    = nkpt
   O(spin)%min_ene = ewin(1,spin)
   O(spin)%max_ene = ewin(2,spin)
   !
   ! The size of the overlap matrix and useful tables.
   ABI_CALLOC(O(spin)%bk2idx,(mband,nkpt))
   ABI_CALLOC(tmp_idx2bk,(2,mband*nkpt))

   ovlp_size=0
   do ik2_bz=1,nkpt
     ik2_ibz  = Kmesh%tab(ik2_bz)
     nband_k2 = Wfd%nband(ik2_ibz,spin)
     do band2=1,nband_k2
       ene_bks = Bands%eig(band2,ik2_ibz,spin)
       ! Is it inside the window?
       if (ene_bks >= O(spin)%min_ene .and. ene_bks <= O(spin)%max_ene) then
         ovlp_size = ovlp_size + 1
         O(spin)%bk2idx(band2,ik2_bz) = ovlp_size
         tmp_idx2bk(1,ovlp_size) = band2
         tmp_idx2bk(2,ovlp_size) = ik2_bz
       end if
     end do
   end do

   O(spin)%size = ovlp_size
   if (O(spin)%size /= int(ovlp_size,kind=i4b)) then
     ! Stop here since we need Lapack i64!
     MSG_ERROR("Size of overlap matrix too large for a default integer!")
   end if

   ABI_MALLOC(O(spin)%idx2bk,(2,ovlp_size))
   O(spin)%idx2bk = tmp_idx2bk(:,1:ovlp_size)
   ABI_FREE(tmp_idx2bk)
   !
   ! Allocate overlap matrix. Could use packed matrix to save memory, but Lapack call is slower.
   write(msg,'(a,i0,a,f12.1,a,i0)')" Memory for the overlap matrix(spin=",spin,"): ",&
&    two*dpc*ovlp_size**2*b2Mb," Mb; Matrix size = ",ovlp_size
   call wrtout(std_out,msg,"COLL")

   ABI_MALLOC(O(spin)%mat, (ovlp_size,ovlp_size))
   O(spin)%mat = czero
   O(spin)%mat_type = TYPE_OVERLAP

   ! Distribute the calculation of the matrix elements among the nodes.
   ! * tstart and t_stop give the initial and final transition index treated by each node.
   ! * my_hsize is the number of transitions treated by this processor
   ! * my_cols(1:2) gives the initial and final column treated by this node.
   !
   ! Distribution of the matrix elements among the nodes.
   os_size = O(spin)%size; tot_nels =  os_size * (os_size + 1_i8b) / 2
   ABI_MALLOC(t_start,(0:nproc_spin-1))
   ABI_MALLOC(t_stop,(0:nproc_spin-1))
   !t_start = 1; t_stop = tot_nels
   !my_rows = [1_i8b,os_size]; my_cols = [1_i8b,os_size]

   call xmpi_split_work2_i8b(tot_nels,nproc_spin,t_start,t_stop,msg,ierr)
   if (ierr==2) MSG_WARNING(msg)

   ABI_MALLOC(hsize_of,(0:nproc_spin-1))
   hsize_of=0
   do rak=0,nproc_spin-1
     if (t_stop(rak)>=t_start(rak)) hsize_of(rak) = t_stop(rak)-t_start(rak)+1
     !write(std_out,*)"tot_nels: ",tot_nels,hsize_of(rak)
   end do

   my_hsize = hsize_of(rank_spin)
   if (my_hsize<=0) then
     write(msg,'(a,i0)')"Wrong number of transitions: my_hsize= ",my_hsize
     MSG_ERROR(msg)
   end if
   if (my_hsize /= int(my_hsize,kind=i4b)) then
     write(msg,'(a,i0)')"Size of local block too large for a default integer, Increase the number of CPUs: my_hsize= ",my_hsize
     MSG_ERROR(msg)
   end if
   ABI_FREE(hsize_of)

   my_cols=0
   do col=1,os_size
     do row=1,col
       oidx = row + col*(col-1_i8b)/2
       if (oidx==t_start(rank_spin)) then
         my_rows(1) = row
         my_cols(1) = col
       end if
       if (oidx==t_stop(rank_spin)) then
         my_rows(2) = row
         my_cols(2) = col
       end if
     end do
   end do
   !
   ! * Announce the treatment of submatrix treated by each node.
   bsize_my_block = 2*dpc*my_hsize
   write(msg,'(4(a,i0))')&
&     ' Treating ',my_hsize,'/',tot_nels,' matrix elements, from column ',my_cols(1),' up to column ',my_cols(2)
   call wrtout(std_out,msg,'PERS')

   ! Computation of the overlap matrix.
   SELECT CASE (use_sym)

   CASE (.FALSE.)
     call wrtout(std_out," Version without symmetries","COLL")

     do ik2_bz=1,nkpt
       call get_BZ_item(Kmesh,ik2_bz,kk2,ik2_ibz,k2_sym,k2_tim,k2_eimkt,k2_umklp,k2_isirred)
       nband_k2 = Wfd%nband(ik2_ibz,spin)
       npw_k2   = Wfd%npwarr(ik2_ibz)

       do band2=1,nband_k2
         col = O(spin)%bk2idx(band2,ik2_bz)

         ! Check is this col is treated by me.
         if (my_cols(2)<col .or. my_cols(1)>col) cycle

         ug2 => Wfd%Wave(band2,ik2_ibz,spin)%ug
         ! Get the symmetrized wavefunction in real space.
         call wfd_sym_ur(Wfd,Cryst,Kmesh,band2,ik2_bz,spin,ur2)

         if (Wfd%usepaw==1) then
           call wfd_get_cprj(Wfd,band2,ik2_ibz,spin,Cryst,Cp_k2,sorted=.FALSE.)
           if (.not.k2_isirred) then
             call paw_symcprj(ik2_bz,nspinor,1,Cryst,Kmesh,Pawtab,Pawang,Cp_k2)
           endif
         end if

         do ik1_bz=1,ik2_bz
           call get_BZ_item(Kmesh,ik1_bz,kk1,ik1_ibz,k1_sym,k1_tim,k1_eimkt,k1_umklp,k1_isirred)
           npw_k1     = Wfd%npwarr(ik1_ibz)
           nband_k1   = Wfd%nband(ik1_ibz,spin)
           band1_stop = nband_k1; if (ik1_bz==ik2_bz) band1_stop = band2
           if (Wfd%usepaw==1) iqpt = k1mk2q(ik1_bz,ik2_bz)

           do band1=1,band1_stop
             row = O(spin)%bk2idx(band1,ik1_bz)

             ! Check if this element is treated by me.
             oidx = row + col*(col-1)/2
             if (oidx<t_start(rank_spin) .or. oidx>t_stop(rank_spin)) cycle

             if (ik1_bz==ik2_bz) then
               ! Save overlap here and skip the computation. Assume input u(g) are normalized!
               blk_ovlp = czero; if (band1==band2) blk_ovlp = cone
               O(spin)%mat(row,col) = blk_ovlp
               !blk_ovlp = wfd_norm2(Wfd,Cryst,Pawtab,band,ik1_ibz,spin) result(norm2)
               cycle
             end if

             if (ik2_bz==ik1_bz .and. k1_isirred) then
               ! Save the overlap here and skip the symmetrization.
               blk_ovlp = wfd_xdotc(Wfd,Cryst,Pawtab,band1,band2,ik1_ibz,spin)
               cycle
             end if

             ! Get the symmetrized wavefunction in real space.
             ! Here we have the SIGFAULT
             call wfd_sym_ur(Wfd,Cryst,Kmesh,band1,ik1_bz,spin,ur1)

             ! Scalar product
             blk_ovlp = xdotc(Wfd%nfft*nspinor,ur1,1,ur2,1)/Wfd%nfft

             if (Wfd%usepaw==1) then
               ! Add onsite contribution.
               call wfd_get_cprj(Wfd,band1,ik1_ibz,spin,Cryst,Cp_k1,sorted=.FALSE.)
               if (.not.k1_isirred) then
                 call paw_symcprj(ik1_bz,nspinor,1,Cryst,Kmesh,Pawtab,Pawang,Cp_k1)
               endif

               !if (iqpt==1) then
               !  onsite = paw_overlap(Cp_k1,Cp_k2,Cryst%typat,Pawtab)
               !  blk_ovlp = blk_ovlp + DCMPLX(onsite(1),onsite(2))
               !else
               paw_ovlp = czero
               call paw_rho_tw_g(1,dim_rtwg1,nspinor,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xred,g_gamma,&
&                Cp_k1,Cp_k2,Pwij(:,iqpt),paw_ovlp)

               ! FIXME not addinf the onsite term works just fine!
               !paw_ovlp = czero
               blk_ovlp = blk_ovlp + paw_ovlp(1)
               !end if
             end if
             !end if
             !
             ! Save overlap.
             O(spin)%mat(row,col) = blk_ovlp
           end do
         end do
       end do
     end do

   CASE (.TRUE.)
     ! FIXME does not work yet. Problems with umklapp symmorphic operations and symmetry tables somewhere.

     call wrtout(std_out," Version with symmetries","COLL")
     MSG_WARNING(" Using version with symmetries (still under development")

     ABI_CHECK(Wfd%usepaw==0,"PAW not coded yet")
     ABI_MALLOC(ovlp_ikfk,(mband,Wfd%nkibz,mband,nkpt))
     ovlp_ikfk=+HUGE(one)

     do ik2_bz=1,nkpt
       call get_BZ_item(Kmesh,ik2_bz,kk2,ik2_ibz,k2_sym,k2_tim,k2_eimkt,k2_umklp)
       nband_k2 = Wfd%nband(ik2_ibz,spin)

       do band2=1,nband_k2
         call wfd_sym_ur(Wfd,Cryst,Kmesh,band2,ik2_bz,spin,ur2)

         do ik1_ibz=1,Wfd%nkibz
           nband_k = Wfd%nband(ik1_ibz,spin)
           do band1=1,nband_k

             call wfd_get_ur(Wfd,band1,ik1_ibz,spin,ur1)
             covlp = xdotc(Wfd%nfft*nspinor,ur1,1,ur2,1) * fft_fact

             if (Wfd%usepaw==1) then ! TODO Add onsite term after IBZ-->BZ symmetrization
             end if

             ovlp_ikfk(band1,ik1_ibz,band2,ik2_bz) = covlp
           end do
         end do
         !
       end do
     end do
     !
     ! FIXME: Equations are not completed. non-symmorphic phases are missing!
     ! Let <r|k> indicate the periodic part of the Bloch wavefunction that transforms according to:
     !   u_{Sk} = e^{-iSk.t} u_k(S^{-1}(r-t))
     !   S3 = S1^{-1} S2
     !
     ! 1) <S1 k1 | S2 k2> = <k1| S1^{-1}   S2 k2>  e^{i (S1 k1 - S2 k2).t1}
     !
     ! 2) <T S1 k1 | T S2 k2> = <S1 k1| S2 k2>*

     ! <S1 k1   | T S2 k2> = <k1| S1^{-1} T S2 k2>
     ! <T S1 k1 |   S2 k2> = <k2| S2^{-1} T S1 k1>^*
     ! <T S1 k1 | T S2 k2> = <k2| S2^{-1}   S1 k1>     ! Problematic if the BZ mesh is not invariant under inversion.
     !                                                   e.g. randomly shifted k-meshes. In this case one should use kptopt=3
     call wrtout(std_out,"Using version with symmetries","COLL")
     with_sym=0; without_sym=0
     do ik2_bz=1,nkpt

       call get_BZ_item(Kmesh,ik2_bz,kk2,ik2_ibz,k2_sym,k2_tim,k2_eimkt,k2_umklp)

       !$ik2_ibz = Kmesh%tab(ik2_bz)
       nband_k2 = Wfd%nband(ik2_ibz,spin)

       inv_k2_sym = toinv(1,k2_sym)
       k_needs_tr(2) = (k2_tim==2)

       do ik1_bz=1,ik2_bz !nkpt
         call get_BZ_item(Kmesh,ik1_bz,kk1,ik1_ibz,k1_sym,k1_tim,k1_eimkt,k1_umklp)

         nband_k1 = Wfd%nband(ik1_ibz,spin)

         inv_k1_sym = toinv(1,k1_sym)

         k_needs_tr(1) = (k1_tim==2)

         !sym3 = multable(1,inv_k1_sym,k2_sym)

         rk= Kmesh%rottb(ik2_bz,1,inv_k1_sym)

         can_use_sym = .TRUE.
         !can_use_sym = ( ALL(ABS(Cryst%tnons(:,k2_sym))<tol6) .and. ALL(ABS(Cryst%tnons(:,k1_sym))<tol6) )
         !can_use_sym = can_use_sym .and. ( ALL(k2_umklp == g_gamma) .and. ALL(k1_umklp == g_gamma))

         can_use_sym = can_use_sym .and. ( ALL(k2_umklp == g_gamma) .and. ALL(k1_umklp == g_gamma) &
&          .and. ALL(multable(2:4,inv_k1_sym,k2_sym)==g_gamma) )

         !can_use_sym = ( can_use_sym .and. &
         !& ALL( ABS ( -Kmesh%bz(:,rk) + MATMUL(Cryst%symrec(:,:,inv_k1_sym),Kmesh%bz(:,ik2_bz)) ) < tol6 ) )

         !can_use_sym = ( can_use_sym .and. ALL(ABS(Cryst%tnons(:,k1_sym)) <tol6)  .and. ALL(ABS(Cryst%tnons(:,k2_sym)) <tol6) )

         kpoint = kk1 - kk2
         !do ii=1,3 ! Wrap in the first BZ thus enforcing traslational invariance.
         !  call wrap2_pmhalf(kk1(ii)-kk2(ii),kpoint(ii),shift)  ! TODO overloaded interface.
         !end do

         if (ANY (ABS(Cryst%tnons(:,k1_sym)) > tol6) ) then
           !tnons_fact = EXP(j_dpc*two_pi*DOT_PRODUCT(kk1-kk2,Cryst%tnons(:,k1_sym)))
           tnons_fact = EXP(j_dpc*two_pi*DOT_PRODUCT(kpoint,Cryst%tnons(:,k1_sym)))
         else
           tnons_fact = cone
         end if

         take_cnj=ALL(k_needs_tr)

         if (ALL(k_needs_tr) .or. ALL(.not.k_needs_tr) ) then
           lk = ik1_ibz
           rk= Kmesh%rottb(ik2_bz,1,inv_k1_sym)
           !rk= Kmesh%rottbm1(ik2_bz,1,k1_sym)
         else
           MSG_ERROR("Need TR")
           take_cnj=.FALSE.
           !if (k_needs_tr(2)) then
           !lk = ik1_ibz
           !rk=Kmesh%rottb(ik2_bz,1,inv_k1_sym)
           !tnons_fact = CONJG(k1_eimkt) * CONJG(k2_eimkt) * EXP(-j_dpc*two_pi*DOT_PRODUCT(kk2,r1_tau3))
           !end if
         end if

         if (can_use_sym) then

           with_sym=with_sym+1
           do band2=1,nband_k2
             col=O(spin)%bk2idx(band2,ik2_bz)

             band1_stop = nband_k1; if (ik1_bz==ik2_bz) band1_stop = band2
             do band1=1,band1_stop
               blk_ovlp = tnons_fact * ovlp_ikfk(band1,lk,band2,rk)
               if (take_cnj) blk_ovlp = DCONJG(blk_ovlp)
               row=O(spin)%bk2idx(band1,ik1_bz)
               if (col>=row) O(spin)%mat(row,col) = blk_ovlp
             end do
           end do

         else
           without_sym=without_sym+1
           do band2=1,nband_k2
             col=O(spin)%bk2idx(band2,ik2_bz)

             call wfd_sym_ur(Wfd,Cryst,Kmesh,band2,ik2_bz,spin,ur2)

             band1_stop = nband_k1; if (ik1_bz==ik2_bz) band1_stop = band2

             do band1=1,band1_stop
               call wfd_sym_ur(Wfd,Cryst,Kmesh,band1,ik1_bz,spin,ur1)

               blk_ovlp = xdotc(Wfd%nfft,ur1,1,ur2,1) * fft_fact
               row=O(spin)%bk2idx(band1,ik1_bz)
               if (col>=row) O(spin)%mat(row,col) = blk_ovlp
             end do
           end do
         end if

       end do
     end do

     write(std_out,*)"with_sym",with_sym," without_sym",without_sym
     ABI_FREE(ovlp_ikfk)
   END SELECT

   ABI_FREE(t_start)
   ABI_FREE(t_stop)
   !
   ! Collect ovlp_mat on each node inside comm_spin.
   call xmpi_sum(O(spin)%mat,comm_spin,ierr)
 end do ! spin

 ABI_FREE(multable)
 ABI_FREE(toinv)
 ABI_FREE(ur1)
 ABI_FREE(ur2)

 ! Deallocate PAW stuff.
 if (Wfd%usepaw==1) then
   ABI_FREE(k1mk2q)
   call pawpwij_free(Pwij)
   ABI_DT_FREE(Pwij)
   call pawcprj_free(Cp_k1)
   ABI_DT_FREE(Cp_k1)
   call pawcprj_free(Cp_k2)
   ABI_DT_FREE(Cp_k2)
 end if

 call cwtime(cpu,wall,gflops,"stop")
 write(std_out,*)"SHTIME: Ovlp_init, cpu_time: ",cpu,", wall_time: ",wall

 DBG_EXIT("COLL")

end subroutine ovlp_init
!!***

!----------------------------------------------------------------------

!!****f* m_shirley/ovlp_diago_and_prune
!! NAME
!! ovlp_diago_and_prune
!!
!! FUNCTION
!!  Diagonalize the overlap matrix for a single spin, and selected the optimal
!!  basis set according the value of sh_coverage.
!!
!! INPUTS
!!  sh_coverage=Coverage factor.
!!
!! OUTPUT
!!  sh_size=Size of the optimal basis set
!!  base=Index of the first eigenvector that is excluded from the optimal basis set.
!!
!! SIDE EFFECTS
!!  O <type(ovlp_t)>= Compute the eigenvalues and store the results in O%eigene
!!
!! PARENTS
!!      m_shirley
!!
!! CHILDREN
!!      blas_cholesky_ortho,cwtime,fft_ur,fftbox_execute,fftbox_plan3_many
!!      flush_unit,get_kg,kgindex,ovlp_diago_and_prune,ovlp_free
!!      wfd_change_ngfft,wfd_get_ur,wfd_init,wfd_print,wfd_push_ug
!!      wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine ovlp_diago_and_prune(O,sh_coverage,sh_size,base)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ovlp_diago_and_prune'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: sh_size,base
 real(dp),intent(in) :: sh_coverage
 type(ovlp_t),intent(inout) :: O

!Local variables ------------------------------
!scalars
 integer :: ii
 real(dp) :: sum_eige,trace,cpu,wall,gflops
 character(len=500) :: msg

!************************************************************************

 DBG_ENTER("COLL")

 ! @ovlp_t
 call cwtime(cpu,wall,gflops,"start")
 !
 ! Diagonalize the overlap matrix.
 ABI_MALLOC(O%eigene,(O%size))

 call xheev("Vectors","Upper",O%size,O%mat,O%eigene)
 O%mat_type = TYPE_EIGVEC

 trace = SUM(O%eigene)
 write(msg,'(3(a,f8.2))')" Trace: ",trace," Min eig: ",MINVAL(O%eigene)," Max eig: ",MAXVAL(O%eigene)
 call wrtout(std_out,msg,"COLL")

 ! Select the optimal subspace.
 sum_eige=zero; ii=O%size
 do while (sum_eige < trace * sh_coverage .and. ii/=1)
   sum_eige = sum_eige + O%eigene(ii)
   ii = ii-1
 end do
 base=ii; sh_size=O%size - base

 ! Include degenerate states (if any) in order to have a symmetric basis set.
 write(std_out,'(a, i0)')"first base: ",base
 do ii=base,1,-1
   if (abs(O%eigene(ii)- O%eigene(base+1)) > tol8) exit
 end do
 base = ii; sh_size=O%size - base
 write(std_out,'(a, i0)')"new base: ",base

 call cwtime(cpu,wall,gflops,"stop")
 write(std_out,*)"SHTIME: Ovlp_diago, cpu_time: ",cpu,", wall_time: ",wall

 write(msg,'(2(a,i0),a)')" The Shirley optimal basis set contains ",sh_size,"/",O%size," elements."
 call wrtout(std_out,msg,"COLL")

 DBG_EXIT("COLL")

end subroutine ovlp_diago_and_prune
!!***

!----------------------------------------------------------------------

!!****f* m_shirley/ovlp_diff
!! NAME
!! ovlp_diff
!!
!! FUNCTION
!!  Debugging tool used to compare the overlap matrices calculated with and without symmetries.
!!
!! INPUTS
!!
!! PARENTS
!!      m_shirley
!!
!! CHILDREN
!!      blas_cholesky_ortho,cwtime,fft_ur,fftbox_execute,fftbox_plan3_many
!!      flush_unit,get_kg,kgindex,ovlp_diago_and_prune,ovlp_free
!!      wfd_change_ngfft,wfd_get_ur,wfd_init,wfd_print,wfd_push_ug
!!      wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine ovlp_diff(O1,O2,Kmesh,tol,header,unit,mode_paral,prtvol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ovlp_diff'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: unit,prtvol
 real(dp),intent(in) :: tol
 character(len=4),optional,intent(in) :: mode_paral
 character(len=*),optional,intent(in) :: header
 type(ovlp_t),intent(in) :: O1,O2
 type(kmesh_t),intent(in) :: Kmesh

!Local variables ------------------------------
!scalars
 integer :: row,col,my_unt,my_prtvol
 integer :: ik1_bz,ik2_bz,ik1_ibz,ik2_ibz,band1
 integer :: k1_sym,k2_sym,k1_tim,k2_tim
 complex(dpc) :: k1_eimkt,k2_eimkt
 character(len=4) :: my_mode
 character(len=500) :: msg
!arrays
 integer :: k1_umklp(3),k2_umklp(3)
 real(dp) :: kk1(3),kk2(3)

!************************************************************************

 MSG_ERROR("Not tested")

 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 msg=' ==== Diff O1-O2  ==== '
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 ABI_CHECK(O1%size==O1%size,"Different sizes")
 !
 ! Compare the overlap matrices
 do col=1,O2%size
   do row=1,O1%size
     !
     if (ABS(O1%mat(row,col)-O2%mat(row,col))>tol) then
       ik1_bz = O1%idx2bk(1,row)
       band1  = O1%idx2bk(2,row)
       call get_BZ_item(Kmesh,ik1_bz,kk1,ik1_ibz,k1_sym,k1_tim,k1_eimkt,k1_umklp)
       !
       ik2_bz = O1%idx2bk(1,col)
       band1  = O1%idx2bk(2,col)
       call get_BZ_item(Kmesh,ik2_bz,kk2,ik2_ibz,k2_sym,k2_tim,k2_eimkt,k2_umklp)

       write(my_unt,'(2i3,2(2x,a,3f7.3),4f8.4)')&
&        row,col," bz1 ",Kmesh%bz(:,ik1_bz)," bz2 ",Kmesh%bz(:,ik2_bz),O1%mat(row,col),O2%mat(row,col)
       write(my_unt,'(a,i3,3i3,2f4.1,3i3)')"k2 ",ik2_bz,ik2_ibz,k2_sym,k2_tim,k2_eimkt,k2_umklp
       write(my_unt,'(a,i3,3i3,2f4.1,3i3)')"k1 ",ik1_bz,ik1_ibz,k1_sym,k1_tim,k1_eimkt,k1_umklp

     end if
   end do
 end do

end subroutine ovlp_diff
!!***

!----------------------------------------------------------------------

!!****f* m_shirley/wfd_bloch_to_shirley
!! NAME
!! wfd_bloch_to_shirley
!!
!! FUNCTION
!!  Build a new wavefunction descriptor containing the  Shirley basis set.
!!
!! INPUT
!! Wfd<wfd_t>=Input wavefunctions.
!! Cryst<crystal_t>= data type gathering info on symmetries and unit cell.
!! Kmesh<kmesh_t>=K-points of the wavefunctions stored in Wfd.
!! Bands<ebands_t>=Input energies.
!! Psps<pseudopotential_type)>=variables related to pseudopotentials.
!! Pawtab(Psps%ntypat)<pawtab_type>=paw tabulated starting data
!! Pawang<pawang_type>=angular mesh discretization and related data:
!! Pawrad<Pawrad_type>=For PAW, RADial mesh discretization and related data
!! min_bsize=Minimum number of states in the optimal basis set. < 0 if no constraint should be enforced.
!! sh_coverage=Parameter defining the coverage of the optimal basis set:
!! ewin(2,nsppol)=Energy window for the spin channels. Use ewin(1,:)=smallest_real and ewin(2,:)=greatest_real
!! to disable the use of the energy window
!!
!! OUTPUT
!!  oWsh<wfd_t>=New wavefunction descriptor with the optimal basis set.
!!    Note that Wfs contains a single k-point, i.e. gamma
!!    oWsh contains N vectors where N is the first elements for which we have:
!!
!!    \sum_i^N e_i >= sh_coverage * M
!!
!!    with e_i being the eigenvalues of the overlap matrix ordered in descending order.
!!    M indicates the dimension of the overlap (related to the number of k-points in the BZ and
!!    the number of bands in the input Wfd descriptor.
!!
!! TODO
!!   Add periodic replica.
!!
!! PARENTS
!!      m_shirley
!!
!! CHILDREN
!!      blas_cholesky_ortho,cwtime,fft_ur,fftbox_execute,fftbox_plan3_many
!!      flush_unit,get_kg,kgindex,ovlp_diago_and_prune,ovlp_free
!!      wfd_change_ngfft,wfd_get_ur,wfd_init,wfd_print,wfd_push_ug
!!      wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine wfd_bloch_to_shirley(Wfd,Cryst,Kmesh,Bands,Psps,Pawtab,Pawang,Pawrad,min_bsize,sh_coverage,ewin,oWsh)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_bloch_to_shirley'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: min_bsize
 real(dp),intent(in) :: sh_coverage
 type(crystal_t),intent(in) :: Cryst
 type(Pawang_type),intent(in) :: Pawang
 type(kmesh_t),intent(in) :: Kmesh
 type(Pseudopotential_type),intent(in) :: Psps
 type(wfd_t),intent(inout) :: Wfd
 type(ebands_t),intent(in) :: Bands
 type(wfd_t),intent(out) :: oWsh
!arrays
 real(dp) :: ewin(2,Wfd%nsppol)
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
 type(Pawrad_type),intent(in) :: Pawrad(Cryst%ntypat*Wfd%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: istwf1=1,k1=1,ndat1=1,method1=1,enforce_sym1=1
 integer :: ik1_bz,band1,natom,sh,midx,my_start,my_stop
 integer :: npw_gamma,spin !,row,col,useylm_,!fft_idx,ig
 integer :: npw_k,nspinor,nsppol,sh_nkibz,sh_mband,ierr
 integer :: ii,ik_ibz,istwf_k,tmp_nfft,comm_spin,nproc_spin,rank_spin
 real(dp) :: fft_fact !,norm1 !sqrt_norm1,sqrt_norm2
 real(dp) :: cpu,wall,gflops
 logical :: use_sym
 character(len=500) :: msg
!arrays
 integer :: ov_ngfft(18),trial_ngfft(18)
 integer :: sh_istwfk(1),sh_size(Wfd%nsppol),base_spin(Wfd%nsppol)
 integer,allocatable :: kg_gamma(:,:)  !kg_k(:,:),
 integer,allocatable :: sh_nband(:,:)
 real(dp) :: sh_kibz(3,1),gamma_point(3)=(/zero,zero,zero/)
 real(dp) :: pawovlp(2),fft_err(3,Cryst%nsym)
 complex(dpc) :: cdum
 !complex(dpc),allocatable :: dpc_tmp(:)
 complex(gwpc),allocatable :: ur1(:) !,ur2(:)
 complex(gwpc),allocatable :: sh_ug(:),sh_ur(:,:)
 logical,allocatable :: sh_keep_ur(:,:,:),sh_bks_mask(:,:,:)
 type(pawcprj_type),allocatable :: Cp_sh1(:,:)
 type(ovlp_t) :: O(Wfd%nsppol),Osym(Wfd%nsppol)

!************************************************************************

 DBG_ENTER("COLL")

 write(msg,'(a,f9.6)')" Starting transformation Bloch ==> Shirley with coverage: ",sh_coverage
 call wrtout(std_out,msg,"COLL")

 ! Consistency check.
 ABI_CHECK(Wfd%nspinor==1,"nspinor==2 not coded.")
 ABI_CHECK(wfd%paral_kgb==0,"paral_kgb/=0 not coded.")
 ABI_CHECK(Wfd%rfft_is_symok,"Real space FFT is not symmetric.")
 ABI_CHECK(Wfd%nsppol==1,"nsppol==2 not coded")
 ABI_CHECK(all(ewin(2,:) > ewin(1,:)), "Wrong energy window")

 if (Wfd%usepaw==1) MSG_WARNING("Shirley with PAW is still under testing")

 ABI_UNUSED(Pawrad(1)%mesh_size)

 nsppol  = Wfd%nsppol; nspinor = Wfd%nspinor
 natom   = Wfd%natom
 !
 ! 1) Get the overlap matrix <u_i|u_j> for this spin.
 ! TODO
 ! *) The Overlap can be calculated using a coarse real space FFT mesh provided that it is compatible
 !    with the space group symmetries  (enforce_sym1=1)
 !    We only have to make sure that ngfft encloses the k-centered G-spheres.
 ! *) Use method1 to get a coarse FFT mesh (no problem due the convolution here)
 ov_ngfft(1:3) = 0; ov_ngfft(7:) = Wfd%ngfft(7:)  ! Have to preserve the previous fftalg options

 do ik_ibz=1,Wfd%nkibz
   istwf_k = Wfd%istwfk(ik_ibz)
   ABI_CHECK(istwf_k==1,"istwf_k/=1 not coded")
   npw_k = Wfd%npwarr(ik_ibz)
   trial_ngfft(7:) = Wfd%ngfft(7:)

   !kg_k => Wfd%Kdata(ik_ibz)%kg_k
   call setmesh(Cryst%gmet,Wfd%Kdata(ik_ibz)%kg_k ,trial_ngfft,npw_k,1,npw_k,tmp_nfft,&
&    method1,[0,0,0],Cryst,enforce_sym1,unit=dev_null)
   do ii=1,3
     ov_ngfft(ii) = MAX(ov_ngfft(ii),trial_ngfft(ii))
   end do
 end do
 !ov_ngfft(4)=2*(ov_ngfft(1)/2)+1
 !ov_ngfft(5)=2*(ov_ngfft(2)/2)+1
 !ov_ngfft(6)=   ov_ngfft(3)

!THIS WAS THE CAUSE OF THE BUG in shirley_window!!!!!!!!!!!!!!!
 ov_ngfft(:) = Wfd%ngfft(:)

 ! TODO
 ! 1) Be careful here in parallel,
 ! 2) One should always include all the degenerated states in the calculation of the
 !    overlap so that the optimal basis set will preserve the degeneracies of the
 !    interpolated eigenvalues.
 use_sym =.FALSE.
 call ovlp_init(O,ewin,use_sym,ov_ngfft,Wfd,Cryst,Kmesh,Bands,Psps,Pawtab,Pawang,Pawrad)

 if (.FALSE.) then
   ! TODO Work in progress
   use_sym =.TRUE.
   call ovlp_init(Osym,ewin,use_sym,ov_ngfft,Wfd,Cryst,Kmesh,Bands,Psps,Pawtab,Pawang,Pawrad)

   do spin=1,nsppol
     call ovlp_diff(O(spin),Osym(spin),Kmesh,tol6,"Error in matrix elements",unit=std_out)
     call ovlp_free(Osym(spin))
   end do
   MSG_ERROR("Check done")
 end if
 !
 ! 2) Diagonalize the overlap matrix selecting the optimal subspace: [base_spin(spin):ovlp_size]
 ! In exit O(spin)%mat stores the eigenvectors (only on those processors treating this spin).
 sh_size = -1; base_spin=HUGE(0)
 do spin=1,nsppol

   if (.not. wfd_itreat_spin(Wfd,spin,comm_spin,rank_spin,nproc_spin)) cycle
   call ovlp_diago_and_prune(O(spin),sh_coverage,sh_size(spin),base_spin(spin))
   !
   ! Make sure we have enough states.
   if (sh_size(spin)<min_bsize) then
     if (O(spin)%size<min_bsize) then
       write(msg,'(2(a,i0),2a)')&
&        "Overlap size is ",O(spin)%size," whereas min_bsize is ",min_bsize,ch10,&
&        "Decrease the number of bands to be interpolated or increase the number of ab-initio input states."
       MSG_ERROR(msg)
     end if
     sh_size(spin) = min_bsize
     write(msg,'(a,2i0)')" Had to enlarge Shirley subspace since input sh_size < min_bsize: ",sh_size(spin),min_bsize
     MSG_COMMENT(msg)
   end if
 end do
 !
 ! 3) Init a new wavefunction descriptor to store the Shirley basis set.
 !
 !    *) oWsh must be allocated here since sh_size is known only after the diagonalization of the overlap.
 !    *) Keep the optimal wavefunctions on each node (if possible) to facilitate the interpolation over the fine k-mesh.
 !    *) Use Gamma-centered basis set to facilitate the operations in reciprocal space.
 !    *) The new basis is orthogonal, but not normalized since <U_i|U_j> = delta_ij e_i.
 !
 ! The optimal basis set is centered on gamma and uses istwfk==1.
 call cwtime(cpu,wall,gflops,"start")

 call get_kg(gamma_point,istwf1,Wfd%ecut,Cryst%gmet,npw_gamma,kg_gamma)
 sh_istwfk = istwf1; sh_nkibz = 1; sh_kibz(:,1) = gamma_point

 ! Compute number of elements in the Shirley basis set and initialize the descriptor oWsh
 ! TODO: BE careful in parallel when nsppol==2. I should recreate the communicators.
 ABI_MALLOC(sh_nband,(sh_nkibz,nsppol))
 do spin=1,nsppol
   sh_nband(:,spin)=sh_size(spin)
 end do
 sh_mband=MAXVAL(sh_nband)

 ABI_MALLOC(sh_bks_mask,(sh_mband,sh_nkibz,nsppol))
 ABI_MALLOC(sh_keep_ur ,(sh_mband,sh_nkibz,nsppol))
 ! Set sh_keep_ur to False. as sh_u(r) are only needed for the computation of vloc so that we can save memory.
 !sh_bks_mask=.TRUE.; sh_keep_ur =.False.
 sh_bks_mask=.TRUE.; sh_keep_ur =.TRUE.

 call wfd_init(oWsh,Cryst,Pawtab,Psps,sh_keep_ur,Wfd%paral_kgb,npw_gamma,sh_mband,sh_nband,sh_nkibz,nsppol,&
&  sh_bks_mask,Wfd%nspden,nspinor,Wfd%ecutsm,Wfd%dilatmx,sh_istwfk,sh_kibz,Wfd%ngfft,kg_gamma,Wfd%nloalg,&
&  Wfd%prtvol,Wfd%pawprtvol,Wfd%comm)

 if (oWsh%prtvol > 0) then
   call wfd_print(oWsh,header="Shirley wavefunction descriptor")
 end if

 ABI_FREE(sh_keep_ur)
 ABI_FREE(sh_bks_mask)
 ABI_FREE(sh_nband)

!DEBUG
! call wfd_change_ngfft(Wfd,Cryst,Psps,ov_ngfft)
! call wfd_change_ngfft(oWsh,Cryst,Psps,ov_ngfft)
!DEBUG
 !
 ! =====================================================================
 ! ==== Rotate the input wavefunctions to get the optimal basis set ====
 ! =====================================================================
 ABI_CHECK(Wfd%nfft == oWsh%nfft, "Different nfft for Wfd and Wsh!")
 fft_fact = one/Wfd%nfft
 ABI_MALLOC(ur1,(Wfd%nfft*nspinor))

 if (.not.fft_check_rotrans(Cryst%nsym,Cryst%symrel,Cryst%tnons,Wfd%ngfft,fft_err)) then
   write(msg,'(a,3(i0,1x),a)')" Real space FFT mesh ",Wfd%ngfft(1:3)," is not symmetric. Cannot symmetrize in real space"
   MSG_ERROR(msg)
 end if

 do spin=1,nsppol
   if (.not. wfd_itreat_spin(Wfd,spin,comm_spin,rank_spin,nproc_spin)) cycle
   write(msg,'(a,f12.1,a)')' Memory needed for Shirley u(r) = ',two*gwpc*Wfd%nfft*nspinor*sh_size(spin)*b2Mb,' [Mb]'
   call wrtout(std_out,msg,'PERS')

   ! Allocate space for Shirley wavefunctions in real space.
   ABI_MALLOC(sh_ur, (Wfd%nfft*nspinor,sh_size(spin)))
   sh_ur = czero

   ! MPI loop over the index of the overlap operator.
   call xmpi_split_work(O(spin)%size,comm_spin,my_start,my_stop,msg,ierr)

   do midx=my_start,my_stop
     ! Retrieve the k-point index in the BZ and the band
     band1  = O(spin)%idx2bk(1,midx)
     ik1_bz = O(spin)%idx2bk(2,midx)

     ! Symmetrize the wavefunction in real space.
     call wfd_sym_ur(Wfd,Cryst,Kmesh,band1,ik1_bz,spin,ur1)

     ! Construct the new optimal basis set.
     do sh=1,sh_size(spin)
       sh_ur(:,sh) = sh_ur(:,sh) + O(spin)%mat(midx,base_spin(spin)+sh) * ur1
     end do
   end do

   ! Gather results
   call xmpi_sum(sh_ur,comm_spin,ierr)

   ! NC pseudos: Normalize the basis set using the eigenvalues of the overlap matrix.
   do sh=1,sh_size(spin)
     sh_ur(:,sh) = sh_ur(:,sh)/SQRT(O(spin)%eigene(base_spin(spin)+sh))
     !norm1 = xdotc(Wfd%nfft*nspinor,sh_ur(:,sh),1,sh_ur(:,sh),1) * fft_fact
     !sh_ur(:,sh) = sh_ur(:,sh)/SQRT(norm1)
     !write(std_out,*)" sh_ur integrates to: ", xdotc(Wfd%nfft*nspinor,sh_ur(:,sh),1,sh_ur(:,sh),1) * fft_fact
   end do
   !
   ! Construct new optimal basis set in G-space and save results in oWsh.
   ABI_MALLOC(sh_ug,(npw_gamma*nspinor))

   do sh=1,sh_size(spin)
     ! Transform sh_u(r) from the FFT mesh to the G-sphere.
     ! Don't try to update u(r) and PAW cpcrj since they are not used here.

     ! gbound_k => oWsh%Kdata(k1)%gbound
     call fft_ur(npw_gamma,oWsh%nfft,nspinor,ndat1,oWsh%mgfft,oWsh%ngfft,istwf1,&
&      kg_gamma,oWsh%Kdata(k1)%gbound,sh_ur(:,sh),sh_ug)
     ! NC: Normalize the basis set using the eigenvalues of the overlap matrix.
     !if (Wfd%usepaw==0) sh_ug = sh_ug/SQRT(O(spin)%eigene(base_spin(spin)+sh))

     ! Store Shirley u(G).
     call wfd_push_ug(oWsh,sh,k1,spin,Cryst,sh_ug,update_ur=.FALSE.,update_cprj=.FALSE.)
   end do
   !
   ! PAW: Normalize the basis set using the eigenvalues of the overlap matrix.
   !if (.FALSE. .and. Wfd%usepaw==1) then
   if (Wfd%usepaw==1) then
     call wfd_reset_ur_cprj(oWsh)
     ABI_DT_MALLOC(Cp_sh1,(natom,nspinor))
     call pawcprj_alloc(Cp_sh1,0,Wfd%nlmn_atm)

     do sh=1,sh_size(spin)
       ! TODO: Write function to compute the scalar product at the same k, taking into account istwfk.
       !ug1 => oWsh%Wave(sh,k1,spin)%ug
       !cdum = xdotc(npw_gamma*nspinor,ug1,1,ug1,1) * Cryst%ucvol
       cdum = xdotc(npw_gamma*nspinor,oWsh%Wave(sh,k1,spin)%ug,1,oWsh%Wave(sh,k1,spin)%ug,1)
       !if (istwf_k>1) then
       !  cdum=two*DBLE(cdum)
       !  if (istwf_k==2) cdum=cdum-CONJG(ug1(1))*ug1(1)
       !end if

       call wfd_get_cprj(oWsh,sh,k1,spin,Cryst,Cp_sh1,sorted=.FALSE.)

       pawovlp = paw_overlap(Cp_sh1,Cp_sh1,Cryst%typat,Pawtab)
       !pawovlp = pawovlp * Cryst%ucvol
       cdum = cdum + CMPLX(pawovlp(1),pawovlp(2))

       !norm1 = xdotc(Wfd%nfft*nspinor,sh_ur(:,sh),1,sh_ur(:,sh),1) * fft_fact * Cryst%ucvol
       !sh_ur(:,sh) = sh_ur(:,sh)/SQRT(O(spin)%eigene(base_spin(spin)+sh))
       !write(std_out,*)" sh_ur integrates to: ", xdotc(Wfd%nfft*nspinor,sh_ur(:,sh),1,sh_ur(:,sh),1) * fft_fact
       !write(std_out,*)" PAW test: ",DBLE(cdum),O(spin)%eigene(base_spin(spin)+sh)

       sh_ug = oWsh%Wave(sh,k1,spin)%ug/SQRT(DBLE(cdum))
       call wfd_push_ug(oWsh,sh,k1,spin,Cryst,sh_ug,update_ur=.FALSE.,update_cprj=.FALSE.)
     end do

     call pawcprj_free(Cp_sh1)
     ABI_DT_FREE(Cp_sh1)
   end if

   ABI_FREE(sh_ur)
   ABI_FREE(sh_ug)
   !
   ! Free the overlap eigenvectors.
   call ovlp_free(O(spin))
 end do ! spin

 do spin=1,nsppol ! Just to be on the safe side.
   call ovlp_free(O(spin))
 end do

 ABI_FREE(kg_gamma)
 ABI_FREE(ur1)

 ! Set the MPI communicators.
 call wfd_set_mpicomm(oWsh)

 call cwtime(cpu,wall,gflops,"stop")
 write(std_out,*)"SHTIME: Rotation, cpu_time: ",cpu,", wall_time: ",wall

 DBG_EXIT("COLL")

end subroutine wfd_bloch_to_shirley
!!***

!----------------------------------------------------------------------

!!****f* m_shirley/shirley_hks
!! NAME
!! shirley_hks
!!
!! FUNCTION
!!  Compute <i|H_ks|j> in the Shirley basis set for a single k-point and spin
!!  This routine is called by a single MPI node who has the entire set
!!  of Shirley wavefunctions in memory. MPI parallelization is done in the caller
!!  at the level of k-points and spins.
!!
!! INPUTS
!! Wsh=Description with Shirley savefunction
!! kpt(3)=K-point in reduced coordinates. Note: The k-point will be wrapped to [-1/2,1/2[
!! spin=Spin index.
!! Ham_k=Ground state Hamiltonian for this k-point and spin.
!! Cryst<crystal_t>= data type gathering info on symmetries and unit cell
!! Psps <type(pseudopotential_type)>=variables related to pseudopotentials
!! Pawtab(Psps%ntypat) <type(pawtab_type)>=paw tabulated starting data
!! Pawang<pawang_type> angular mesh discretization and related data:
!! Paw_ij(natom)<type(paw_ij_type)>=data structure containing PAW arrays given on (i,j) channels.
!! sh_size=Size of the Shirley basis set.
!! vloc_ij(sh_size,sh_size)=Matrix element of the local part in the Shirley basis set.
!!
!! OUTPUT
!!  hk_ij(sh_size,sh_size)= <i|H_ks|j>
!!  sk_ij(sh_size,sh_size*Psps%usepaw)=PAW matrix.
!!
!! PARENTS
!!      m_shirley
!!
!! CHILDREN
!!      blas_cholesky_ortho,cwtime,fft_ur,fftbox_execute,fftbox_plan3_many
!!      flush_unit,get_kg,kgindex,ovlp_diago_and_prune,ovlp_free
!!      wfd_change_ngfft,wfd_get_ur,wfd_init,wfd_print,wfd_push_ug
!!      wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine shirley_hks(Wsh,kpt,spin,Ham_k,Cryst,Psps,Pawtab,Pawang,Paw_ij,sh_size,vloc_ij,hk_ij,sk_ij)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'shirley_hks'
 use interfaces_56_recipspace
 use interfaces_66_nonlocal
 use interfaces_69_wfdesc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spin,sh_size
 type(crystal_t),intent(in) :: Cryst
 type(Pawang_type),intent(in) :: Pawang
 type(Pseudopotential_type),intent(in) :: Psps
 type(wfd_t),target,intent(inout) :: Wsh
 type(gs_hamiltonian_type),intent(inout) :: Ham_k
!arrays
 real(dp),intent(in) :: kpt(3)
 complex(dpc),intent(in) :: vloc_ij(sh_size,sh_size)
 complex(dpc),intent(out) ::  hk_ij(sh_size,sh_size),sk_ij(sh_size,sh_size*Psps%usepaw)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
 type(paw_ij_type),intent(in) :: Paw_ij(Cryst%natom*Psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: idir0=0,ider0=0,nnlout0=0,tim_nonlop0=0,k1=1
 integer :: ig,natom,npw_k,optder,matblk,mkmem_,nkpg,dimffnl,nspinortot
 integer :: blc1,blc2,istwf_k,useylm_,iat,nspinor
 complex(dpc) :: kin_ij,vnl_ij
 !character(len=500) :: msg
 type(kdata_t) :: Kdata
!arrays
 integer :: nloalg(3)
 integer,pointer :: kg_k(:,:)
 integer,allocatable :: nlmn_sort(:)
 real(dp) :: k4intp(3),kptns_(3,1),ylmgr_dum(1,1,1),shifts(3)
 real(dp),allocatable :: kdotg(:),half_gsq(:),ylm_k(:,:),dum_ylm_gr_k(:,:,:)
 real(dp),pointer :: ffnl(:,:,:,:)
 real(dp),allocatable :: kpg_k(:,:),ph3d(:,:,:),vnl_psi(:,:),vectin(:,:),s_psi(:,:)
 complex(gwpc),pointer :: ug1(:)
 complex(gwpc),allocatable :: cvnl_psi(:),cs_psi(:),wsg(:)
 type(pawcprj_type),allocatable :: Cp_blc2(:,:)

!************************************************************************

 DBG_ENTER("COLL")

 ABI_UNUSED(Paw_ij%cplex)

 ABI_CHECK(Wsh%nspinor==1,"nspinor==2 not coded")
 ABI_CHECK(Wsh%nsppol==1,"Wsh%nsppol must be 1")
 ABI_CHECK(Wsh%paral_kgb/=1,"paral_kgb not coded")

 natom   = Cryst%natom
 useylm_ = Psps%useylm
 nloalg  = Wsh%nloalg
 nspinor = Wsh%nspinor
 nspinortot = nspinor

 ! Gamma-centered basis set.
 npw_k   = Wsh%Kdata(k1)%npw
 istwf_k = Wsh%istwfk(k1)
 ABI_CHECK(istwf_k==1,"istwfk/=0 not coded")
 kg_k => Wsh%Kdata(k1)%kg_k

 ! First wrap in the first BZ thus enforcing traslational invariance.
 call wrap2_pmhalf(kpt(:),k4intp(:),shifts(:))

 ! Continue to prepare the GS Hamiltonian.
 call load_spin_hamiltonian(Ham_k,spin,with_nonlocal=.true.)

 call kdata_init(Kdata,Cryst,Psps,k4intp,istwf_k,Wsh%ngfft,Wsh%MPI_enreg,kg_k=kg_k)

 ABI_MALLOC(nlmn_sort,(Cryst%natom))
 iat=0 ! nlmn dims sorted by atom type.

 if (Psps%usepaw==1) then
   nlmn_sort = Wsh%nlmn_sort
 else  ! FIXME here lmnmax == lnmax if useylm_==0
   nlmn_sort = 9
   if (Psps%useylm==1) then
     nlmn_sort=Psps%lmnmax
   else
     MSG_ERROR("useylm==0 not coded")
   end if
   !do itypat=1,Cryst%ntypat
   !  nlmn_sort(iat+1:iat+Cryst%nattyp(itypat))=Pawtab(itypat)%lmn_size
   !  iat=iat+Cryst%nattyp(itypat)
   !end do
   !write(std_out,*)" hacking nlmn_sort",nlmn_sort
   !write(std_out,*)" Psps%lmnmax is ",Psps%lmnmax
 end if
 !
 ! ============================================================================
 ! ==== Evaluate <p_lmn|e^(ikr)U_i> for each k on the k-grid  and each U_i ====
 ! ============================================================================
 !
 ! Here I assume that the G-sphere is gamma-centered.
 ! Real Spherical Harmonics are always used to apply the non-local part even for NC pseudos.
 ! I did non find any easy way to extract only <p_nl|psi> from nonlop_pl.
 !useylm_=1

 ABI_DT_MALLOC(Cp_blc2 ,(natom,nspinor))
 call pawcprj_alloc(Cp_blc2, 0,nlmn_sort)

 ! =======================
 ! ==== Interpolation ====
 ! =======================
 !
 ! Prepare the kinetic term.
 ABI_MALLOC(kdotg,(npw_k))
 ABI_MALLOC(half_gsq,(npw_k))
 ABI_MALLOC(wsg,(npw_k))

 ! TODO Add new overloaded interface. effmass_free option!
 do ig=1,npw_k
   kdotg(ig)    = two_pi**2 * DOT_PRODUCT(k4intp,MATMUL(Cryst%gmet,kg_k(:,ig)))
   half_gsq(ig) = half * vdotw(one*kg_k(:,ig),one*kg_k(:,ig),Cryst%gmet,"G")
 end do

 ABI_MALLOC(vectin,(2, npw_k*nspinor))
 ABI_MALLOC(vnl_psi,(2,npw_k*nspinor))
 ABI_MALLOC(cvnl_psi,(npw_k*nspinor))

 ABI_MALLOC(s_psi,(2,npw_k*nspinor*Psps%usepaw))
 ABI_MALLOC(cs_psi,(npw_k*nspinor*Psps%usepaw))

 ! THIS PART IS NEEDED FOR THE CALL TO opernl although some quantities won't be used.
 ! Now I do things cleanly then we try to pass zero-sized arrays!
 ABI_MALLOC(ylm_k,(npw_k,Psps%mpsang**2*useylm_))

 if (useylm_==1) then
   kptns_(:,1)=k4intp; optder=0; mkmem_=1
   ABI_MALLOC(dum_ylm_gr_k,(npw_k,3+6*(optder/2),Psps%mpsang**2))

   !  Here mband is not used if paral_compil_kpt=0
   call initylmg(Cryst%gprimd,kg_k,kptns_,mkmem_,Wsh%MPI_enreg,Psps%mpsang,npw_k,(/1/),1,&
&    (/npw_k/),1,optder,Cryst%rprimd,ylm_k,dum_ylm_gr_k)

   ABI_FREE(dum_ylm_gr_k)
 end if
 !
 ! Compute (k+G) vectors (only if useylm_=1)
 nkpg=3*nloalg(3)
 ABI_MALLOC(kpg_k,(npw_k,nkpg))
 if (nkpg>0) then
   call mkkpg(kg_k,kpg_k,k4intp,nkpg,npw_k)
 endif
 !
 ! ========================================================
 ! ==== Compute nonlocal form factors ffnl at all (k+G) ====
 ! ========================================================
 dimffnl=1+ider0 ! Derivatives are not needed.
 ABI_MALLOC(ffnl,(npw_k,dimffnl,Psps%lmnmax,Psps%ntypat))
 !ffnl => Kdata%fnl_dir0der0
 call mkffnl(Psps%dimekb,dimffnl,Psps%ekb,ffnl,Psps%ffspl,Cryst%gmet,Cryst%gprimd,ider0,idir0,Psps%indlmn,&
&  kg_k,kpg_k,k4intp,Psps%lmnmax,Psps%lnmax,Psps%mpsang,Psps%mqgrid_ff,nkpg,npw_k,&
&  Psps%ntypat,Psps%pspso,Psps%qgrid_ff,Cryst%rmet,Psps%usepaw,useylm_,ylm_k,ylmgr_dum)
 ABI_FREE(ylm_k)
 !
 ! Load k-dependent part in the Hamiltonian datastructure
 matblk=min(NLO_MINCAT,maxval(Ham_k%nattyp)) ; if (nloalg(2)>0) matblk=natom
 ABI_MALLOC(ph3d,(2,npw_k,matblk))
 call load_k_hamiltonian(Ham_k,kpt_k=k4intp,npw_k=npw_k,istwf_k=istwf_k,kg_k=kg_k,&
&                        kpg_k=kpg_k,ffnl_k=ffnl,ph3d_k=ph3d,compute_ph3d=(Wsh%paral_kgb/=1))

!END BORING CODE NEEDED TO CALL opernl

 !
 ! Calculate the upper triangle of <blc2| H_k |blc1>.
 do blc2=1,sh_size

   ! Calculate <G|Vnl|psi> for this k-point
   call wfd_vnlpsi(Wsh,blc2,k1,spin,npw_k,Cryst,Ham_k,vnl_psi,s_psi,Kext=Kdata)
   cvnl_psi = DCMPLX(vnl_psi(1,:),vnl_psi(2,:))
   if (Psps%usepaw==1) cs_psi = DCMPLX(s_psi(1,:),s_psi(2,:))

   ! Upper triangle of the hk_ij matrix
   do blc1=1,blc2
     ug1 => Wsh%Wave(blc1,k1,1)%ug

     ! Treat the matrix elements of the Vnl operator.
     vnl_ij = xdotc(npw_k*nspinor,ug1,1,cvnl_psi,1)
     !
     ! ===================================================
     ! ==== Assemble final Hamiltonian T + Vloc + Vnl ====
     ! ===================================================
     !
     ! Kinetic energy.
     wsg = kdotg * Wsh%Wave(blc2,k1,1)%ug
     kin_ij = xdotc(npw_k,ug1,1,wsg,1)

     wsg = half_gsq * Wsh%Wave(blc2,k1,1)%ug
     kin_ij = kin_ij + xdotc(npw_k,ug1,1,wsg,1)
     if (blc1==blc2) kin_ij = kin_ij + half * vdotw(k4intp,k4intp,Cryst%gmet,"G")
     !
     ! Total Hamiltonian.
     hk_ij(blc1,blc2) = kin_ij + vloc_ij(blc1,blc2) + vnl_ij
     !
     ! PAW Overlap operator.
     if (Psps%usepaw==1) sk_ij(blc1,blc2) = xdotc(npw_k*nspinor,ug1,1,cs_psi,1)
   end do ! blc1
 end do ! blc2

!DEALLOCATE BORING STUFF
 ABI_FREE(ffnl)
 ABI_FREE(kpg_k)
 ABI_FREE(ph3d)
!END BORING STUFF

 ABI_FREE(kdotg)
 ABI_FREE(half_gsq)
 ABI_FREE(wsg)
 ABI_FREE(vectin)
 ABI_FREE(vnl_psi)
 ABI_FREE(cvnl_psi)
 ABI_FREE(s_psi)
 ABI_FREE(cs_psi)

 call pawcprj_free(Cp_blc2 )
 ABI_DT_FREE(Cp_blc2)
 call kdata_free(Kdata)

 ABI_FREE(nlmn_sort)

 DBG_EXIT("COLL")

 RETURN

 ABI_UNUSED(Pawang%ngnt)
 ABI_UNUSED(Pawtab(1)%basis_size)

end subroutine shirley_hks
!!***

!----------------------------------------------------------------------

!!****f* m_shirley/shirley_interp
!! NAME
!! shirley_interp
!!
!! FUNCTION
!!  Compute interpolated KS eigenvalues and optionally KS eigenvectors by direct
!!  diagonalization of the KS Hamiltonian in the Shirley basis set.
!!  Main entry point for client code.
!!
!! INPUTS
!! jobz="V" if eigevectors are wanted. "N" if only eigenvalues are required.
!! Dtset=Input variables for this run (it will be removed later on)
!! Cryst<crystal_t>= data type gathering info on symmetries and unit cell
!! Psps <pseudopotential_type>=variables related to pseudopotentials
!! Pawtab(Psps%ntypat) <type(pawtab_type)>=paw tabulated starting data
!! Pawfgr<Pawfgr_type>
!! Pawang<pawang_type> angular mesh discretization and related data:
!! Pawrad<Pawrad_type>
!! Pawrhoij
!! Paw_ij(natom)<paw_ij_type>=data structure containing PAW arrays given on (i,j) channels.
!! ngfftc(18)=Information about the coarse 3D FFT.
!! ngfftf(18)=Information about the dense 3D FFT used for vtrial.
!! nfftf=Number of points in the FFT grid in vtrial. Might differ from the FFT mesh used for the wavefunctions.
!! vtrial(nfftf,nspden)= Trial potential (Hartree)
!! intp_nband(intp_nk,Wsh%nsppol)=Number of interpolated bands
!! intp_mband=Max number of interpolated bands (used for dimensioning arrays)
!! intp_nk=Number of interpolated k-points.
!! intp_kpt(3,intp_nk)=Reduced coordinates of the interpolated k-points.
!! sh_fname=String with the filename of the netcdf file to be produced. No file is created if empty string.
!! [kweights(intp_nk)]=Weights of the k-points used in oBands, Assume 0 if not present
!!
!! OUTPUT
!!  oBands=Interpolated band energies. Note: weights are initialized to zero.
!!  [oWfd]
!!
!! PARENTS
!!      m_shirley,wfk_analyze
!!
!! CHILDREN
!!      blas_cholesky_ortho,cwtime,fft_ur,fftbox_execute,fftbox_plan3_many
!!      flush_unit,get_kg,kgindex,ovlp_diago_and_prune,ovlp_free
!!      wfd_change_ngfft,wfd_get_ur,wfd_init,wfd_print,wfd_push_ug
!!      wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine shirley_interp(Wsh,jobz,Dtset,Cryst,Psps,Pawtab,Pawfgr,Pawang,Pawrad,&
&  Pawrhoij,Paw_ij,ngfftc,ngfftf,nfftf,vtrial,&
&  intp_nband,intp_mband,intp_nk,intp_kpt,sh_fname,oBands,&
&  kweights,oWfd) ! optional


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'shirley_interp'
 use interfaces_14_hidewrite
 use interfaces_67_common
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftf,intp_nk,intp_mband
 character(len=*),intent(in) :: jobz,sh_fname
 type(wfd_t),intent(inout) :: Wsh
 type(crystal_t),intent(in) :: Cryst
 type(Pawang_type),intent(in) :: Pawang
 type(Pseudopotential_type),intent(in) :: Psps
 type(dataset_type),intent(in) :: Dtset
 type(Pawfgr_type),intent(in) :: Pawfgr
 type(ebands_t),intent(out) :: oBands
!arrays
 integer,intent(in) :: ngfftf(18),ngfftc(18)
 integer,intent(in) :: intp_nband(intp_nk,Wsh%nsppol)
 real(dp),intent(in) :: vtrial(nfftf,Wsh%nspden),intp_kpt(3,intp_nk)
 real(dp),optional,intent(in) :: kweights(intp_nk)
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wsh%usepaw)
 type(paw_ij_type),intent(in) :: Paw_ij(Cryst%natom*Wsh%usepaw)
 type(Pawrad_type),intent(in) :: Pawrad(Cryst%ntypat*Wsh%usepaw)
 type(pawrhoij_type),intent(in) :: Pawrhoij(Cryst%natom*Wsh%usepaw)
 type(wfd_t),optional,intent(out) :: oWfd

!Local variables ------------------------------
!scalars
 integer,parameter :: istwf1=1,k1=1
#ifdef HAVE_NETCDF
 integer,parameter :: dummy_nshiftk1=1,dummy_nshiftk_orig1=1
#endif
 integer :: ii,ib,jj,ierr,nband_k,ikpt,natom,nefound,band
 integer :: sh1,sh2,ldz,prtvol,pawprtvol,istwf_k
 integer :: mgfftc,nfftc,onband_diago,usepaw,spin
 integer :: comm_spin,nproc_spin,rank_spin,root
 integer :: npw_k,nspinor,nsppol,nspden,sh_size,intp_bantot
 integer :: my_shstart,my_shstop,my_kstart,my_kstop
 real(dp) :: fft_fact,ene_fact,cpu,wall,gflops
 !logical,parameter :: debug_with_diago=.TRUE.
 logical,parameter :: debug_with_diago=.FALSE.
 logical :: want_eigenvectors,do_full_diago
 character(len=500) :: msg,frmt1,frmt2
 type(ddiago_ctl_type) :: Diago_ctl
 type(gs_hamiltonian_type) :: Ham_k
 type(stats_t) :: Stats
!arrays
#ifdef HAVE_NETCDF
 integer,parameter :: dummy_ngkpt(3)=0,dummy_kptrlatt(3,3)=0
 real(dp),parameter :: dummy_shiftk(3,dummy_nshiftk1)=zero,dummy_shiftk_orig(3,dummy_nshiftk_orig1)=zero
#endif
 integer :: nloalg(3),intp_npwarr(intp_nk),intp_istwfk(intp_nk)
 integer,allocatable :: nlmn_sort(:)
 real(dp) :: intp_ene(intp_mband,intp_nk,Wsh%nsppol)
 real(dp) :: intp_wtk(intp_nk),kpoint(3)
 real(dp),pointer :: diag_ene(:),diag_vec(:,:,:)
 real(dp),allocatable :: enek_ij(:),vnl_psi(:,:),opaw_psi(:,:)
 real(dp),allocatable :: copy_vtrial(:,:)
 real(dp),allocatable :: intp_doccde(:),intp_occ(:),ugly_ene(:)
 complex(gwpc),allocatable :: ur1(:),ur2(:),vloc_psi(:)
 complex(dpc),allocatable :: hk_ij(:,:),sk_ij(:,:),eig_vec(:,:),vloc_ij(:,:)
 character(len=10) :: spin_name(2)
 type(pawcprj_type),pointer :: diag_Cprj(:,:)
!BEGIN For the output wfd.
 integer :: sh_npw
 logical,allocatable :: keep_ur(:,:,:),bks_mask(:,:,:)
 complex(gwpc),allocatable :: ug(:)
!END For the output wfd.
#ifdef HAVE_NETCDF
 integer :: ncid
#endif

!************************************************************************

 DBG_ENTER("COLL")

 !call wrtout(std_out,"Starting Shirley interpolation on: "asctime(),"COLL")

 ABI_UNUSED(Pawrhoij(1)%cplex)
 ABI_UNUSED(Pawrad(1)%mesh_size)

 ABI_CHECK(Wsh%nspinor==1,"nspinor==2 not coded")
 ABI_CHECK(Wsh%paral_kgb==0,"paral_kgb/=0 not coded")
 ABI_CHECK(Wsh%rfft_is_symok,"FFT not symmetric in real space")

 if (maxval(intp_nband) > maxval(Wsh%npwarr)) then
    write(msg,'(a,i0,2a,i0)')&
&    "The number of bands to be interpolated: ",maxval(intp_nband),ch10,&
&    "cannot be greater than the size of the Shirley basis: ",maxval(Wsh%npwarr)
    MSG_ERROR(msg)
 end if

 call cwtime(cpu,wall,gflops,"start")

 ! Copy basic dimensions and parameters.
 nspinor   = Wsh%nspinor
 nsppol    = Wsh%nsppol
 nspden    = Wsh%nspden
 nloalg    = Wsh%nloalg
 usepaw    = Wsh%usepaw
 natom     = Wsh%natom
 prtvol    = Wsh%prtvol
 pawprtvol = Wsh%pawprtvol

 ! The coarse mesh used for the Hamiltonian.
 nfftc  = PRODUCT(ngfftc(1:3))
 mgfftc = MAXVAL(ngfftc(1:3))

 if (nsppol==1) spin_name=(/'          ','          '/)
 if (nsppol==2) spin_name=(/'SPIN_UP   ','SPIN_DOWN '/)

 want_eigenvectors = firstchar(jobz,"V",csens=.FALSE.)
 if (want_eigenvectors) then
   call wrtout(std_out," Eigenvectors will be calculated." ,"COLL")
 else
   call wrtout(std_out," Interpolating energies only. Eigenvectors won't be calculated.","COLL")
 end if

!BEGIN DEBUG: Test on the orthogonalization of the input wavefunctions.
! call wfd_reset_ur_cprj(Wsh)
! call wfd_test_ortho(Wsh,Cryst,Pawtab,unit=std_out,mode_paral="COLL")
!END DEBUG

 ! vtrial might be given on a FFT mesh that is denser than the FFT used for Wsh.
 ! If the two FFTs differ, change the mesh for the wavefunctions.
 ! Another possibility would be moving vtrial from the dense to the coarse mesh.
 call wfd_change_ngfft(Wsh,Cryst,Psps,ngfftf)

 ! Init a new wavefunction descriptor to store the interpolated KS states
 ! The optimal basis set is given on the gamma centered basis set with istwfk==1.
 ! Note that the basis set in oWfd is not k-centered since its u(g) are given
 ! in terms of a linear combination of Shirley u(g) that are Gamma-centered
 if (present(oWfd)) then
   ABI_CHECK(want_eigenvectors,"When oWfd is present, want_eigenvectors must be true")

   ! TODO: BE careful in parallel when nsppol==2. I should recreate the communicators.
   ABI_MALLOC(bks_mask,(intp_mband,intp_nk,nsppol))
   ABI_MALLOC(keep_ur ,(intp_mband,intp_nk,nsppol))
   bks_mask=.TRUE.; keep_ur =.TRUE.

   ! Build new wavefunction descriptor (Gamma centered)
   !sh_gvec => Wsh%Kdata(1)%kg_k
   intp_istwfk = istwf1; sh_npw = Wsh%npwarr(1)

   call wfd_init(oWfd,Cryst,Pawtab,Psps,keep_ur,Wsh%paral_kgb,sh_npw,intp_mband,intp_nband,intp_nk,Wsh%nsppol,bks_mask,&
&    Wsh%nspden,Wsh%nspinor,Wsh%ecutsm,Wsh%dilatmx,intp_istwfk,intp_kpt,Wsh%ngfft,&
&    Wsh%Kdata(k1)%kg_k,Wsh%nloalg,Wsh%prtvol,Wsh%pawprtvol,Wsh%comm)

   if (oWfd%prtvol > 0) then
     call wfd_print(oWfd,header="New wavefunction descriptor with interpolated states")
   end if

   ABI_FREE(bks_mask)
   ABI_FREE(keep_ur)
 end if
 !
 ! =======================
 ! ==== Interpolation ====
 ! =======================
 fft_fact = one/Wsh%nfft
 ABI_MALLOC(ur1,(Wsh%nfft*nspinor))
 ABI_MALLOC(ur2,(Wsh%nfft*nspinor))

 intp_ene=zero

 ! Precompute <sh1| vloc_spin |sh2> as it is not k-dependent.
 do spin=1,nsppol
   if (.not. wfd_itreat_spin(Wsh,spin,comm_spin,rank_spin,nproc_spin)) cycle

   ! Split the k-points inside comm_spin
   !my_kstart=1; my_kstop=intp_nk
   call xmpi_split_work(intp_nk,comm_spin,my_kstart,my_kstop,msg,ierr)
   if (ierr==2) MSG_WARNING(msg)
   write(std_out,*)"Will treat from my_kstart ",my_kstart,"to my_kstop ",my_kstop
   !
   ! Compute the upper triangle of the <sh2|vloc|sh1> matrix.
   ! All MPI(spin) processors participate
   ! Be careful here as |\tpsi> is not normalized.
   ABI_MALLOC(vloc_psi,(Wsh%nfft*nspinor))

   sh_size = Wsh%nband(k1,spin)
   ABI_MALLOC(vloc_ij,(sh_size,sh_size))
   vloc_ij = czero

   !my_shstart=1; my_shstop=sh_size
   call xmpi_split_work(sh_size,comm_spin,my_shstart,my_shstop,msg,ierr)
   if (ierr==2) MSG_WARNING(msg)

   do sh2=my_shstart,my_shstop
     call wfd_get_ur(Wsh,sh2,k1,spin,ur2)
     if (nspinor==1) then
       vloc_psi = ur2*vtrial(:,spin)
     else
       MSG_ERROR("vloc_psi doesn't support nspinor==2")
     end if
     !
     ! Diagonal element.
     vloc_ij(sh2,sh2) = xdotc(Wsh%nfft,ur2,1,vloc_psi,1)*fft_fact
     !
     ! Upper triangle.
     do sh1=1,sh2-1
       call wfd_get_ur(Wsh,sh1,k1,spin,ur1)
       vloc_ij(sh1,sh2) = xdotc(Wsh%nfft,ur1,1,vloc_psi,1)*fft_fact
     end do
   end do !sh2
   ABI_FREE(vloc_psi)

   ! Gather results.
   call xmpi_sum(vloc_ij,comm_spin,ierr)
   !
   ! =============================================
   ! ==== Loop over the interpolated k-points ====
   ! =============================================

   call init_hamiltonian(Ham_k,Psps,Pawtab,nspinor,nsppol,nspden,natom,&
&    Cryst%typat,Cryst%xred,Wsh%nfft,Wsh%mgfft,Wsh%ngfft,Cryst%rprimd,nloalg)

   ABI_MALLOC(hk_ij, (sh_size,sh_size))
   ABI_MALLOC(sk_ij,(sh_size,sh_size*usepaw))

   npw_k = Wsh%Kdata(k1)%npw

   ABI_MALLOC(vnl_psi,(2,npw_k*nspinor))
   ABI_MALLOC(opaw_psi,(2,npw_k*nspinor*usepaw))

   ! Loop over the K-points for the interpolation (MPI loop)
   ! (NB: k will be wrapped in ]-1/2,1/2]).
   do ikpt=my_kstart,my_kstop
     istwf_k = 1 ! FIXME no time-reversal tricks!
     kpoint  = intp_kpt(:,ikpt)
     nband_k = intp_nband(ikpt,spin)

     ! Compute <i|H_k|j> in the Shirley basis set.
     call shirley_hks(Wsh,kpoint,spin,Ham_k,Cryst,Psps,Pawtab,Pawang,Paw_ij,sh_size,vloc_ij,hk_ij,sk_ij)
     !
     ! Diagonalize Hk_ij in the optimal Bloch supspace.
     ABI_MALLOC(enek_ij,(sh_size))

     do_full_diago = (nband_k == sh_size)
     if (.not.do_full_diago) then
       ldz=1; if (want_eigenvectors) ldz=sh_size
       ABI_MALLOC(eig_vec,(ldz,nband_k))
     end if
     !do_full_diago = .True.
     !write(std_out,*)"full diago",do_full_diago

     if (usepaw==0) then
       ! Solve H*v = e*v
       if (do_full_diago) then
         call xheev(jobz,"Upper",sh_size,hk_ij,enek_ij)
       else
         call xheevx(jobz,"Irange","Upper",sh_size,hk_ij,zero,zero,1,nband_k,-tol8,nefound,enek_ij,eig_vec,ldz)
         if (want_eigenvectors) hk_ij(:,1:ldz) = eig_vec
       end if
     else
       ! Solve H*v = e*S*v
       if (do_full_diago) then
         call xhegv(1,jobz,"Upper",sh_size,hk_ij,sk_ij,enek_ij)
       else
         call xhegvx(1,jobz,"Irange","Upper",sh_size,hk_ij,sk_ij,zero,zero,1,nband_k,-tol8,nefound,enek_ij,eig_vec,ldz)
         if (want_eigenvectors) hk_ij(:,1:ldz) = eig_vec
       end if
     end if
     !
     ! Store the interpolated eigenvalues and eigenstates in the optimal basis set.
     if (want_eigenvectors) then
       ! Compute interpolated eigenvectors in G-space.
       ABI_CHECK(oWfd%npwarr(ikpt) == Wsh%Kdata(k1)%npw, "Mismatch oWfd, Wsh")
       ABI_MALLOC(ug, (nspinor*oWfd%npwarr(ikpt)))
       do band=1,nband_k
         ug = czero
         do ii=1,sh_size
            ug(:) = ug(:) + hk_ij(ii,band) * Wsh%Wave(ii,k1,spin)%ug
         end do
         call wfd_push_ug(oWfd,band,ikpt,spin,Cryst,ug,update_ur=.FALSE.,update_cprj=.FALSE.)
       end do
       ABI_FREE(ug)
     end if

     if (allocated(eig_vec)) then
       ABI_FREE(eig_vec)
     end if

     if (prtvol>0) then
       ! Write interpolated energies.
       ene_fact=Ha_eV; frmt1='(i4,4x,9(1x,f7.4))'; frmt2='(8x,9(1x,f7.4))'
       write(msg,'(a,3es16.8,2a)')' Eigenvalues in eV for kpt= ( ',kpoint," ), spin ",spin_name(spin)
       call wrtout(std_out,msg,'COLL')

       write(msg,frmt1)ikpt,(enek_ij(ib)*ene_fact,ib=1,MIN(9,nband_k))
       call wrtout(std_out,msg,'COLL')

       if (nband_k>9) then
         do jj=10,nband_k,9
           write(msg,frmt2) (enek_ij(ib)*ene_fact,ib=jj,MIN(jj+8,nband_k))
           call wrtout(std_out,msg,'COLL')
         end do
       end if
       call flush_unit(std_out)
     end if
     !
     ! Save the interpolated energies.
     intp_ene(1:nband_k,ikpt,spin) = enek_ij(1:nband_k)

     if (debug_with_diago) then ! Compare interpolated energies wrt direct diago results.

       ! FIXME here there is a problem with Wd%ecut and Dtset%ecut.
       !if (ikpt==1) write(std_out,*)" CHECK: Dtset%ecut=",Dtset%ecut," Wsh%ecut= ",Wsh%ecut
       !call init_ddiago_ctl(Diago_ctl,"No Vectors",spin,nspinor,Wsh%ecut,kpoint,nloalg,Cryst%gmet,&
       !&   nband_k=nband_k,effmass_free=Dtset%effmass_free,istwf_k=istwf1,prtvol=prtvol)

       call init_ddiago_ctl(Diago_ctl,"No Vectors",spin,nspinor,Dtset%ecut,kpoint,nloalg,Cryst%gmet,&
&        nband_k=nband_k,effmass_free=Dtset%effmass_free,istwf_k=istwf1,prtvol=prtvol)

       nullify(diag_ene)
       nullify(diag_vec)
       nullify(diag_Cprj)

       ABI_MALLOC(copy_vtrial,(nfftf,nspden))
       copy_vtrial = vtrial ! To avoid intent inout

       call ks_ddiago(Diago_ctl,nband_k,nfftc,mgfftc,ngfftc,natom,&
&        Cryst%typat,nfftf,nspinor,nspden,nsppol,Pawtab,Pawfgr,Paw_ij,&
&        Psps,Cryst%rprimd,copy_vtrial,Cryst%xred,onband_diago,diag_ene,&
&        diag_vec,diag_Cprj,xmpi_comm_self,ierr)
       ABI_CHECK(ierr==0,"Fatal error. Cannot continue!")

       ABI_FREE(copy_vtrial)

       if (.TRUE..or.prtvol>0) then
         ene_fact=Ha_eV
         write(msg,'(a,3es16.8,2a)')' Diago Eigenvalues in eV for kpt= ( ',kpoint," ), spin ",spin_name(spin)
         call wrtout(std_out,msg,'COLL')

         write(msg,'(i4,4x,9(1x,f7.4))')ikpt,(ene_fact*diag_ene(ii),ii=1,MIN(9,nband_k))
         call wrtout(std_out,msg,'COLL')

         if (nband_k>9) then
           do jj=10,nband_k,9
             write(msg,'(8x,9(1x,f7.4))')(ene_fact*diag_ene(ii),ii=jj,MIN(jj+8,nband_k))
             call wrtout(std_out,msg,'COLL')
           end do
         end if
       end if

       ! Compute statistical parameters
       Stats = stats_eval(ABS(diag_ene(1:nband_k-2) - enek_ij(1:nband_k-2))*Ha_eV*1000)

       write(std_out,'(a,f7.2,a,i0,a,i0)')&
&        " MAX abs error diag_ene - intp_ene ",MAXVAL(ABS(diag_ene(1:nband_k) - enek_ij(1:nband_k)) )*Ha_eV*1000,&
&        " [meV] for band ",imax_loc(ABS(diag_ene(1:nband_k) - enek_ij(1:nband_k)) ),"/",nband_k
       write(std_out,'(4(a,es9.1),a)')&
&        " min: ",Stats%min,", Max: ",Stats%max,", mean: ",Stats%mean,", stdev: ",Stats%stdev," [meV]"
       call flush_unit(std_out)

       if (associated(diag_ene)) then
         ABI_FREE(diag_ene)
       end if
       if (associated(diag_vec)) then
         ABI_FREE(diag_vec)
       end if
       if (associated(diag_Cprj)) then
       end if
     end if

     ABI_FREE(enek_ij)
   end do ! intp_kpt
   !
   ! Interpolation completed, deallocate memory.
   ABI_FREE(hk_ij)
   ABI_FREE(vloc_ij)
   ABI_FREE(sk_ij)
   ABI_FREE(vnl_psi)

   ABI_FREE(opaw_psi)
   call destroy_hamiltonian(Ham_k)
 end do ! spin
 !
 ! Collect energies on each node.
 call xmpi_sum(intp_ene,Wsh%comm,ierr)
 !
 ! Free memory
 ABI_FREE(ur1)
 ABI_FREE(ur2)
 if (allocated(nlmn_sort)) then
   ABI_FREE(nlmn_sort)
 end if

 ! TODO
 ! If vectors we have to write the unitary transformation with MPI-IO

 ! Create ebands_t for storing the interpolated energies.
 intp_bantot=SUM(intp_nband)
 ABI_CALLOC(intp_doccde,(intp_bantot))
 ABI_CALLOC(intp_occ,(intp_bantot))
 intp_istwfk=1; intp_npwarr=Wsh%npwarr(k1) ! Useless but oh well

 ! Have to reshape intp_ene --> ugly_ene due to the ugly convention used in abinit to store energies.
 ABI_MALLOC(ugly_ene,(intp_bantot))
 ugly_ene(:)=HUGE(zero)

 call pack_eneocc(intp_nk,nsppol,intp_mband,intp_nband,intp_bantot,intp_ene,ugly_ene)

 ! Warning: weigths are set to zero if kweights are not specified.
 intp_wtk = zero; if (present(kweights)) intp_wtk = kweights

 MSG_ERROR("Fix ebands_init")

 !call ebands_init(intp_bantot,oBands,Dtset%nelect,intp_doccde,ugly_ene,intp_istwfk,intp_kpt,&
 ! intp_nband,intp_nk,intp_npwarr,nsppol,nspinor,Dtset%tphysel,Dtset%tsmear,Dtset%occopt,intp_occ,intp_wtk,&
 ! charge, kptopt, kptrlatt_orig, nshiftk_orig, shiftk_orig, kptrlatt, nshiftk, shiftk)

 ABI_FREE(intp_doccde)
 ABI_FREE(ugly_ene)
 ABI_FREE(intp_occ)

 if (present(oWfd)) then
   ABI_MALLOC(bks_mask, (oWfd%mband,oWfd%nkibz,oWfd%nsppol))
   ! TODO
   ! Collect all the wavefunctions on each MPI node.
   do spin=1,oWfd%nsppol
     if (.not. wfd_itreat_spin(oWfd,spin,comm_spin,rank_spin,nproc_spin)) cycle
     do root=0,nproc_spin-1
       bks_mask=.False. !do ik_ibz=my_kstart,my_kstop
       bks_mask(:,my_kstart:my_kstop,spin) = .True.
       !call wfd_bcast_waves(Wfd,what,bks_mask,root,spin_comm,reinit_comms=.False.)
     end do
   end do
   ! Set the MPI communicators in the output wavefunction descriptor.
   ! Master nodes in the spin communicator must communicate the bks_table
   ! and broadcast the updated values to the nodes in the spin communicator.
   if (oWfd%nsppol == 2) MSG_ERROR("Not implemented error")
   call wfd_set_mpicomm(oWfd)
   ABI_FREE(bks_mask)
 end if

 call cwtime(cpu,wall,gflops,"stop")
 write(std_out,*)"SHTIME: Interpolation, cpu_time: ",cpu,", wall_time: ",wall
 call flush_unit(std_out)

 ! Write netcdf file with results.
 ! FIXME: k-point info are wrong. This trick is needed so that
 ! abipy will detect a path instead of a BZ mesh.
#ifdef HAVE_NETCDF
 if (Wsh%my_rank == 0 .and. len_trim(sh_fname)>0) then
   NCF_CHECK_MSG(nctk_open_create(ncid, sh_fname, xmpi_comm_self), "Creating Shirley file")
   NCF_CHECK(crystal_ncwrite(Cryst, ncid))
   NCF_CHECK(ebands_ncwrite(oBands, ncid))
   NCF_CHECK(nf90_close(ncid))
 end if
#endif

 DBG_EXIT("COLL")

end subroutine shirley_interp
!!***

!----------------------------------------------------------------------

!!****f* m_shirley/shirley_bands
!! NAME
!!  shirley_bands
!!
!! FUNCTION
!!  Helper function for the interpolation of the eigenvalues with the Shirley method
!!
!! INPUT
!! iWfd<wfd_t>=Input wavefunctions.
!! Dtset=Input variables for this run (it will be removed later on)
!! Cryst<crystal_t>= data type gathering info on symmetries and unit cell.
!! iKmesh<kmesh_t>=K-points of the wavefunctions stored in Wfd.
!! iBands<ebands_t>=Input energies.
!! Psps<pseudopotential_type)>=variables related to pseudopotentials.
!! Pawtab(Psps%ntypat)<pawtab_type>=paw tabulated starting data
!! Pawang<pawang_type>=angular mesh discretization and related data:
!! Pawrad<Pawrad_type>=For PAW, RADial mesh discretization and related data
!! Pawfgr<Pawfgr_type>
!! Pawang<pawang_type> angular mesh discretization and related data:
!! Pawrad<Pawrad_type>
!! Pawrhoij
!! Paw_ij(natom)<paw_ij_type>=data structure containing PAW arrays given on (i,j) channels.
!! ngfftc(18)=Information about the coarse 3D FFT.
!! ngfftf(18)=Information about the dense 3D FFT used for vtrial.
!! nfftf=Number of points in the FFT grid in vtrial. Might differ from the FFT mesh used for the wavefunctions.
!! vtrial(nfftf,nspden)= Trial potential (Hartree)
!! intp_mband=Number of bands to intepolate.
!! sh_coverage=Parameter defining the coverage of the optimal basis set:
!!   Wsh contains N vectors where N is the first elements for which we have:
!!
!!   \sum_i^N e_i >= sh_coverage * M
!!
!!   with e_i being the eigenvalues of the overlap matrix ordered in descending order.
!!   M indicates the dimension of the overlap (related to the number of k-points in the BZ and
!!   the number of bands in the input iWfd descriptor.
!! ewin(2,nsppol)=Energy window for the spin channels. Use ewin(1,:)=smallest_real and ewin(2,:)=greatest_real
!! to disable the use of the energy window
!! min_bsize=Minimum number of states in the optimal basis set.
!! sh_fname=String with the filename of the netcdf file to be produced. No file is created if empty string.
!! kpts(3,nk)=List of k-points for the interpolation.
!! kweights(nk)=K-point weights (stored in oBands)
!!
!! OUTPUT
!!  oWsh<wfd_t>=New wavefunction descriptor with the optimal basis set.
!!  Note that Wfs contains a single (fake) k-point.
!!
!! PARENTS
!!      wfk_analyze
!!
!! CHILDREN
!!      blas_cholesky_ortho,cwtime,fft_ur,fftbox_execute,fftbox_plan3_many
!!      flush_unit,get_kg,kgindex,ovlp_diago_and_prune,ovlp_free
!!      wfd_change_ngfft,wfd_get_ur,wfd_init,wfd_print,wfd_push_ug
!!      wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine shirley_bands(iWfd,Dtset,Cryst,iKmesh,iBands,Psps,Pawtab,Pawang,Pawrad,Pawfgr,Pawrhoij,Paw_ij,&
&  sh_coverage,ewin,ngfftc,ngfftf,nfftf,vtrial,intp_mband,kpts,kweights,sh_fname,oBands,oWsh)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'shirley_bands'
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftf,intp_mband
 real(dp),intent(in) :: sh_coverage
 character(len=fnlen),intent(in) :: sh_fname
 type(crystal_t),intent(in) :: Cryst
 type(Pawang_type),intent(in) :: Pawang
 type(kmesh_t),intent(in) :: iKmesh
 type(Pseudopotential_type),intent(in) :: Psps
 type(wfd_t),intent(inout) :: iWfd
 type(ebands_t),intent(in) :: iBands
 type(dataset_type),intent(in) :: Dtset
 type(Pawfgr_type),intent(in) :: Pawfgr
 type(wfd_t),intent(out) :: oWsh
 type(ebands_t),intent(out) :: oBands
!arrays
 integer,intent(in) :: ngfftf(18),ngfftc(18)
 real(dp),intent(in) :: vtrial(nfftf,iWfd%nspden),kpts(:,:),kweights(:)
 real(dp),intent(in) :: ewin(2,iWfd%nsppol)
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*iWfd%usepaw)
 type(Pawrad_type),intent(in) :: Pawrad(Cryst%ntypat*iWfd%usepaw)
 type(paw_ij_type),intent(in) :: Paw_ij(Cryst%natom*iWfd%usepaw)
 type(pawrhoij_type),intent(in) :: Pawrhoij(Cryst%natom*iWfd%usepaw)

!Local variables ------------------------------
!scalars
 integer :: min_bsize,npts
!arrays
 integer,allocatable :: intp_nband(:,:)

!************************************************************************

 DBG_ENTER("COLL")

 ! Compute shirley basis set
 min_bsize=intp_mband
 call wfd_bloch_to_shirley(iWfd,Cryst,iKmesh,iBands,Psps,Pawtab,Pawang,Pawrad,min_bsize,sh_coverage,ewin,oWsh)

 npts = size(kpts,dim=2)
 ABI_MALLOC(intp_nband, (npts,iWfd%nsppol))
 intp_nband=intp_mband

! Energies only
 call shirley_interp(oWsh,"N",Dtset,Cryst,Psps,Pawtab,Pawfgr,Pawang,Pawrad,&
&  Pawrhoij,Paw_ij,ngfftc,ngfftf,nfftf,vtrial,&
&  intp_nband,intp_mband,npts,kpts,sh_fname,oBands,kweights)

 ! Free memory
 ABI_FREE(intp_nband)

 DBG_EXIT("COLL")

end subroutine shirley_bands
!!***

!----------------------------------------------------------------------

!!****f* m_shirley/shirley_window
!! NAME
!!  shirley_window
!!
!! FUNCTION
!!  Find an optimal basis set for the interpolation of the KS states in a energy window, e.g states around Ef.
!!  The algorithm starts with a relatively coarse k-mesh and an initial energy window (sufficiently large so that
!!  the initial basis set has enough flexibility). Constructs the initial Shirley set and uses it to interpolate
!!  the KS states on a denser mesh. The energy window is progressively reduced so that the states far from Ef
!!  are progressively filtered out (this trick allows one to reduced the size of the overlap matrix and therefore the CPU time
!!  needed for its diagonalization. The iterative algorithm stops when the DOS is converged within a given tolerance.
!!  Returns a new wavefunction descriptor with the Shirley basis set.
!!
!! INPUT
!! iWfd<wfd_t>=Input wavefunctions.
!! Dtset=Input variables for this run (it will be removed later on)
!! Cryst<crystal_t>= data type gathering info on symmetries and unit cell.
!! iKmesh<kmesh_t>=K-points of the wavefunctions stored in Wfd.
!! iBands<ebands_t>=Input energies.
!! Psps<pseudopotential_type)>=variables related to pseudopotentials.
!! Pawtab(Psps%ntypat)<pawtab_type>=paw tabulated starting data
!! Pawang<pawang_type>=angular mesh discretization and related data:
!! Pawrad<Pawrad_type>=For PAW, RADial mesh discretization and related data
!! Pawfgr<Pawfgr_type>
!! Pawang<pawang_type> angular mesh discretization and related data:
!! Pawrad<Pawrad_type>
!! Pawrhoij
!! Paw_ij(natom)<paw_ij_type>=data structure containing PAW arrays given on (i,j) channels.
!! ngfftc(18)=Information about the coarse 3D FFT.
!! ngfftf(18)=Information about the dense 3D FFT used for vtrial.
!! nfftf=Number of points in the FFT grid in vtrial. Might differ from the FFT mesh used for the wavefunctions.
!! vtrial(nfftf,nspden)= Trial potential (Hartree)
!! sh_coverage=Parameter defining the coverage of the optimal basis set:
!!   Wsh contains N vectors where N is the first elements for which we have:
!!
!!   \sum_i^N e_i >= sh_coverage * M
!!
!!   with e_i being the eigenvalues of the overlap matrix ordered in descending order.
!!   M indicates the dimension of the overlap (related to the number of k-points in the BZ and
!!   the number of bands in the input iWfd descriptor.
!! ewin(2,nsppol)=Energy window for the spin channels. Use ewin(1,:)=smallest_real and ewin(2,:)=greatest_real
!! to disable the use of the energy window
!! min_bsize=Minimum number of states in the optimal basis set.
!!
!! SIDE EFFECTS
!!  ewin(2,nsppol)
!!
!! OUTPUT
!!  oWsh<wfd_t>=New wavefunction descriptor with the optimal basis set.
!!  Note that the descriptor contains a single (fake) k-point.
!!
!! PARENTS
!!      wfk_analyze
!!
!! CHILDREN
!!      blas_cholesky_ortho,cwtime,fft_ur,fftbox_execute,fftbox_plan3_many
!!      flush_unit,get_kg,kgindex,ovlp_diago_and_prune,ovlp_free
!!      wfd_change_ngfft,wfd_get_ur,wfd_init,wfd_print,wfd_push_ug
!!      wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine shirley_window(iWfd,Dtset,Cryst,iKmesh,iBands,Psps,Pawtab,Pawang,Pawrad,Pawfgr,Pawrhoij,Paw_ij,&
&  sh_coverage,ewin,ngfftc,ngfftf,nfftf,vtrial,oWsh)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'shirley_window'
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftf
 real(dp),intent(in) :: sh_coverage
 type(crystal_t),intent(in) :: Cryst
 type(Pawang_type),intent(in) :: Pawang
 type(kmesh_t),intent(in) :: iKmesh
 type(Pseudopotential_type),intent(in) :: Psps
 type(wfd_t),intent(inout) :: iWfd
 type(ebands_t),intent(in) :: iBands
 type(dataset_type),intent(in) :: Dtset
 type(Pawfgr_type),intent(in) :: Pawfgr
 type(wfd_t),intent(out) :: oWsh
!arrays
 integer,intent(in) :: ngfftf(18),ngfftc(18)
 real(dp),intent(in) :: vtrial(nfftf,iWfd%nspden)
 real(dp),intent(inout) :: ewin(2,iWfd%nsppol)
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*iWfd%usepaw)
 type(Pawrad_type),intent(in) :: Pawrad(Cryst%ntypat*iWfd%usepaw)
 type(paw_ij_type),intent(in) :: Paw_ij(Cryst%natom*iWfd%usepaw)
 type(pawrhoij_type),intent(in) :: Pawrhoij(Cryst%natom*iWfd%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: max_niter=2
 integer :: min_bsize,intp_mband,intp_nk,iter,nsppol !,ierr
 logical :: converged
 character(len=500) :: sh_fname
 type(wfd_t) :: tmpWfd
 !type(edos_t),pointer :: newEdos,oldEdos
 !type(edos_t),target :: Doses(2)
 type(kmesh_t) :: oKmesh ! <<<<
 type(ebands_t) :: oBands ! <<<<
!arrays
 integer :: kptrlatt(3,3)
 integer,allocatable :: intp_nband(:,:)

!************************************************************************

 DBG_ENTER("COLL")

 ! Compute the DOS from the input bands and k-points.
 !call edos_init(Doses(1),iBands,iKmesh,method,step,broad)
 !call edos_free(Doses(1))
 !call edos_free(Doses(2))
 !oldEdos => Doses(1)
 !newEdos => Doses(2)

 nsppol=iWfd%nsppol
 min_bsize=-1

 ! Compute Shirley basis set with the initial Kmesh and the initial window energy.
 call wfd_bloch_to_shirley(iWfd,Cryst,iKmesh,iBands,Psps,Pawtab,Pawang,Pawrad,min_bsize,sh_coverage,ewin,oWsh)

 do iter=1,max_niter
   ! Densify the k-mesh
   !if (iter==1) kptrlatt = reshape([5,0,0,0,5,0,0,0,5], [3,3])
   !if (iter==2) kptrlatt = reshape([6,0,0,0,6,0,0,0,6], [3,3])
   if (iter==1) kptrlatt = reshape([7,0,0,0,7,0,0,0,7], [3,3])
   if (iter==2) kptrlatt = reshape([8,0,0,0,8,0,0,0,8], [3,3])

   !if (iter==1) kptrlatt = reshape([2,0,0,0,2,0,0,0,2], [3,3])
   !if (iter==2) kptrlatt = reshape([3,0,0,0,3,0,0,0,3], [3,3])
   !if (iter==3) kptrlatt = reshape([4,0,0,0,4,0,0,0,4], [3,3])
   !if (iter==4) kptrlatt = reshape([5,0,0,0,5,0,0,0,5], [3,3])

   call make_mesh(oKmesh,Cryst,iKmesh%kptopt,kptrlatt,Dtset%nshiftk,Dtset%shiftk)
   write(std_out,*)"Moving to denser Kmesh with kptrlatt ",kptrlatt

   ! Interpolate energies and vectors on the denser mesh.
   intp_nk=oKmesh%nibz

   ABI_MALLOC(intp_nband, (intp_nk,nsppol))
   ! this has a huge impact on the CPU time (full versus partial diago)
   ! Cannot have more eigenvectors than the number of basis elements, doh!
   intp_nband=10
   if (maxval(intp_nband) > oWsh%nband(1,1)) intp_nband = oWsh%nband(1,1)
   intp_mband=maxval(intp_nband)

   ! Set the weights in oBands
   write(sh_fname,"(a,i0,a)")"SHWITER",iter,"_GSR.nc"
   call shirley_interp(oWsh,"V",Dtset,Cryst,Psps,Pawtab,Pawfgr,Pawang,Pawrad,&
&    Pawrhoij,Paw_ij,ngfftc,ngfftf,nfftf,vtrial,&
&    intp_nband,intp_mband,oKmesh%nibz,oKmesh%ibz,sh_fname,oBands,kweights=oKmesh%wt,oWfd=tmpWfd)

   if (tmpWfd%prtvol > 0) then
     call wfd_print(tmpWfd,header="Shirley wavefunction descriptor")
   end if
   if (DEBUGME) then
     call wfd_test_ortho(tmpWfd,Cryst,Pawtab,unit=std_out,mode_paral="COLL")
   end if

   ! Free memory
   ABI_FREE(intp_nband)

   converged = .False.
   ! Compute DOS and DOS(e0) and compare wrt previous one. Use same mesh as oldEdos
   ! Decrease the energy window if not converged and iterate
   !call edos_init(newEdos,oBands,oKmesh,method,step,broad,ierr,mesh=oldEdos%mesh)
   !ws = bisect(newEdos%mesh, ene_min)
   !we = bisect(newEdos%mesh, ene_max)
   !if (iter>1) then
   !converged = (l2norm(newEdos%dos(ws:we) - oldEdos%dos(ws:we)) < dos_atol .and. &
   !&           abs(newEdos%gef(0) - oldEdos%gef(0)) < dos_atol)
   !call edos_free(newEdos)
   if (converged) exit

   ewin = ewin
   !ewin = 0.8*ewin
   write(std_out,*)"Finding new Shirley basis set with energy window. Iteration: ",iter
   write(std_out,*)"Energy window: ",ewin*Ha_eV," [eV]"

   call wfd_free(oWsh)
   call wfd_bloch_to_shirley(tmpWfd,Cryst,oKmesh,oBands,Psps,Pawtab,Pawang,Pawrad,min_bsize,sh_coverage,ewin,oWsh)

   call wfd_free(tmpWfd)

   !if (iter==3) exit
   !if (.not. converged) then
   call kmesh_free(oKmesh)
   call ebands_free(oBands)
   !end if
 end do

 DBG_EXIT("COLL")

end subroutine shirley_window
!!***

!----------------------------------------------------------------------

!!****f* m_shirley/wfd_shirley_to_eh
!! NAME
!! wfd_shirley_to_eh
!!
!! FUNCTION
!!  Return a new wavefunction descriptor containing the basis set for the e-h manifold.
!!
!! INPUT
!! Wfd<wfd_t>
!! Cryst<crystal_t>= data type gathering info on symmetries and unit cell
!! Psps<pseudopotential_type)>=variables related to pseudopotentials
!! Pawtab(Psps%ntypat)<pawtab_type>=paw tabulated starting data
!! Pawang <pawang_type>=angular mesh discretization and related data:
!! eh_coverage
!!
!! OUTPUT
!!  Weh<wfd_t>
!!
!! PARENTS
!!
!! CHILDREN
!!      blas_cholesky_ortho,cwtime,fft_ur,fftbox_execute,fftbox_plan3_many
!!      flush_unit,get_kg,kgindex,ovlp_diago_and_prune,ovlp_free
!!      wfd_change_ngfft,wfd_get_ur,wfd_init,wfd_print,wfd_push_ug
!!      wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine wfd_shirley_to_eh(Wsh,Cryst,Psps,Pawtab,Pawang,Pawrad,min_bsize,eh_coverage,Weh,sh2eh)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_shirley_to_eh'
 use interfaces_14_hidewrite
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: min_bsize
 real(dp),intent(in) :: eh_coverage
 type(crystal_t),intent(in) :: Cryst
 type(Pawang_type),intent(in) :: Pawang
 type(Pseudopotential_type),intent(in) :: Psps
 type(wfd_t),target,intent(inout) :: Wsh
 type(wfd_t),target,intent(out) :: Weh
!arrays
 !integer,intent(in) :: ov_ngfft(18)
 complex(gwpc),pointer :: sh2eh(:,:)
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wsh%usepaw)
 type(Pawrad_type),intent(in) :: Pawrad(Cryst%ntypat*Wsh%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: istwf1=1,k1=1,s1=1,nkpt1=1,ndat1=1
 integer :: ij,band1,band2,natom,eh,midx,ig,mband,nstates,ovlp_size
 integer :: npw_gamma,usepaw !,spin !,row,col,useylm_,!,ierr
 integer :: fft_idx,nspinor,eh_nkibz,eh_mband,npw,col,band3,band4,nsppol
 real(dp) :: fft_fact,norm1 !sqrt_norm1,sqrt_norm2
 real(dp) :: cpu,wall,gflops
 character(len=500) :: msg
 type(fftbox_plan3_t) :: plan
 type(ovlp_t) :: Oeh
!arrays
 integer :: ov_ngfft(18)
 integer :: eh_istwfk(1),eh_size,base
 integer,allocatable :: kg_gamma(:,:) !,gbound_k(:,:),kg_k(:,:),
 integer,allocatable :: eh_nband(:,:)
 integer,pointer :: igfft0(:)
 real(dp) :: eh_kibz(3,1),gamma_point(3)=(/zero,zero,zero/)
 !real(dp) :: pawovlp(2)
 complex(dpc) :: cfft_fact
 complex(dpc),allocatable :: dpc_tmp(:)
 complex(gwpc),allocatable :: ur1(:),ur2(:),ur12(:),cf_ovlp(:,:)
 complex(gwpc),target,allocatable :: ur3(:),ur4(:) !ur1(:),ur2(:),
 !complex(gwpc),pointer :: ug1(:),pt_ur1(:),pt_ur2(:)
 complex(gwpc),pointer :: pt_ur3(:),pt_ur4(:)
 complex(dpc),allocatable :: ur34_big(:,:)
 complex(gwpc),allocatable :: ehg(:),eh_ur(:,:),eh_ug(:,:),gwpc_tmp(:)
 logical,allocatable :: eh_keep_ur(:,:,:),eh_bks_mask(:,:,:),kg_mask(:)

!************************************************************************

 DBG_ENTER("COLL")

 ABI_CHECK(Wsh%nspinor==1,"nspinor==2 not coded.")
 ABI_CHECK(Wsh%paral_kgb==0,"paral_kgb/=0 not coded.")
 ABI_CHECK(Wsh%rfft_is_symok,"Real space FFT is not symmetric.")
 ABI_CHECK(Wsh%nsppol==1,"nsppol must be 1!")

 ABI_UNUSED(Pawang%l_max)
 ABI_UNUSED(Pawrad(1)%mesh_size)

 nspinor   = Wsh%nspinor
 usepaw    = Wsh%usepaw
 natom     = Wsh%natom
 nsppol    = Wsh%nsppol
 mband     = Wsh%mband
 nstates   = Wsh%nband(k1,s1)
 ABI_CHECK(xmpi_comm_size(Wsh%comm)==1,"ovlp2_init is not parallelized")

 write(msg,'(a,f9.6)')" Transformation Shirley --> e-h basis set with eh_coverage : ",eh_coverage
 call wrtout(std_out,msg,"COLL")

 ! 1) Get the overlap matrix <S_i^* S_j|S_k^* S_l> for this spin.
 ! The Overlap is calculated using the coarse real space FFT ov_ngfft
 ov_ngfft = Wsh%ngfft
 call wfd_change_ngfft(Wsh,Cryst,Psps,ov_ngfft)

 ovlp_size   = nstates**2
 Oeh%mband   = mband
 Oeh%nkpt    = nkpt1
 Oeh%min_ene = smallest_real
 Oeh%max_ene = greatest_real
 !
 ! The size of the overlap matrix and useful tables.
 ABI_MALLOC(Oeh%bk2idx,(mband,nkpt1))
 Oeh%bk2idx=0

 Oeh%size = ovlp_size
 ABI_MALLOC(Oeh%idx2bk,(2,ovlp_size))
 Oeh%idx2bk = 0
 !
 ! Allocate overlap matrix. Could use packed matrix to save memory, but Lapack call is slower.
 write(msg,'(a,f8.2,a)')" out of memory in Oeh%mat, requiring :",two*dpc*ovlp_size**2*b2Gb," Gb"
 ABI_MALLOC(Oeh%mat, (ovlp_size,ovlp_size))

 !Oeh%mat = -HUGE(one)
 write(msg,'(a,f12.1,a,i0)')" Memory required for the overlap matrix: ",two*dpc*ovlp_size**2*b2Mb," Mb; Matrix size= ",ovlp_size
 call wrtout(std_out,msg,"COLL")
 call flush_unit(std_out)
 !
 ! Calculate the overlap matrix  --------------------------------------------------------------------------
 ! 1) Symmetrization in k-space is not needed here
 ! 2) Matrix is Hermitian.
 !
 fft_fact = one/Wsh%nfft
 ABI_MALLOC(ur3,(Wsh%nfft*nspinor))
 ABI_MALLOC(ur4,(Wsh%nfft*nspinor))

 ! TODO Temporary implementation used to speed up this part.
 ABI_MALLOC(ur34_big,(Wsh%nfft,nstates**2))

 Oeh%mat_type = TYPE_OVERLAP
 npw = Wsh%npwarr(k1)

 call cwtime(cpu,wall,gflops,"start")
 col = 0
 do band4=1,nstates
   !
   if (wfd_ihave_ur(Wsh,band4,k1,s1,how="Stored")) then
     pt_ur4 =>  Wsh%Wave(band4,k1,s1)%ur
   else
     call wfd_get_ur(Wsh,band4,k1,s1,ur4)
     pt_ur4 =>  ur4
   end if
   !
   do band3=1,nstates
     col = col+1
     if (wfd_ihave_ur(Wsh,band3,k1,s1,how="Stored")) then
       pt_ur3 =>  Wsh%Wave(band3,k1,s1)%ur
     else
       call wfd_get_ur(Wsh,band3,k1,s1,ur3)
       pt_ur3 =>  ur3
     end if
     !if (Wsh%usepaw==1) then
     !call wfd_get_cprj(Wsh,band3,k1,s1,Cryst,Cp_k3,sorted=.FALSE.)
     !end if
     ur34_big(:,col) = CONJG(pt_ur3) * pt_ur4
     Oeh%idx2bk(1,col) = band3
     Oeh%idx2bk(2,col) = band4
   end do
 end do

 cfft_fact = cone/Wsh%nfft
 !Oeh%mat = cfft_fact * MATMUL( CONJG(TRANSPOSE(ur34_big)), ur34_big)
 call xgemm("C","N",nstates**2,nstates**2,Wsh%nfft,cfft_fact,ur34_big,Wsh%nfft,ur34_big,Wsh%nfft,czero,Oeh%mat,nstates**2)

 !call print_arr(Oeh%mat,max_r=10,max_c=12,unit=std_out,mode_paral="COLL")

 ABI_FREE(ur3)
 ABI_FREE(ur4)

 call cwtime(cpu,wall,gflops,"stop")
 write(std_out,*)"SHTIME: Ovlp2_build, cpu_time: ",cpu,", wall_time: ",wall
 call flush_unit(std_out)

 ! 2) Diagonalize the overlap matrix selecting the optimal subspace: [ base(spin):ovlp_size ]
 call ovlp_diago_and_prune(Oeh,eh_coverage,eh_size,base) ! In exit Oeh%mat stores the eigenvectors.
 !
 ! Make sure we have enough states.
 if (eh_size < min_bsize) then
   if (Oeh%size<min_bsize) then
     write(msg,'(2(a,i0),2a)')&
&      " Overlap size is ",Oeh%size," whereas min_bsize is ",min_bsize,ch10,&
&      " Decrease the number of bands to be interpolated or increase the number of ab-initio input states."
     MSG_ERROR(msg)
   end if
   eh_size = min_bsize
   write(msg,'(a,2i0)')" Had to enlarge Shirley subspace since input eh_size < min_bsize: ",eh_size,min_bsize
   MSG_COMMENT(msg)
 end if

 call cwtime(cpu,wall,gflops,"start")
 !
 ! 3) Init a new wavefunction descriptor to store the optimal basis set.
 !    *) Weh must be allocated here since eh_size is know only after the diagonalization of the overlap.
 !    *) Keep the optimal wavefunctions on each node (if possible) to facilitate the interpolation over the fine k-mesh.
 !    *) Use Gamma-centered basis set to facilitate the operations in reciprocal space.
 !    *) The new basis is orthogonal, but not normalized since <U_i|U_j> = delta_ij e_i.
 !
 ! The optimal basis set is given on the gamma centered basis set with istwfk==1.
 ! FIXME temporary hacking. There is a bug somewhere in kpgsph
 call get_kg(gamma_point,istwf1,Wsh%ecut,Cryst%gmet,npw_gamma,kg_gamma)
 !call get_kg(gamma_point,istwf1,14.0_dp,Cryst%gmet,npw_gamma,kg_gamma)
 !
 ! * Index of the G-sphere in the FFT box.
 ABI_MALLOC(igfft0,(npw_gamma))
 ABI_MALLOC(kg_mask,(npw_gamma))
 call kgindex(igfft0,kg_gamma,kg_mask,Wsh%MPI_enreg,Wsh%ngfft,npw_gamma)

 ABI_CHECK(ALL(kg_mask),"FFT para not yet implemented")
 ABI_FREE(kg_mask)

 eh_istwfk=istwf1; eh_nkibz= 1
 eh_kibz(:,1) = gamma_point

 ! TODO: BE careful in parallel when nsppol==2. I should recreate the communicators.
 ABI_MALLOC(eh_nband,(eh_nkibz,nsppol))
 eh_nband = eh_size
 !call wfd_change_ngfft(Wsh,Cryst,Psps,ov_ngfft)

 eh_mband=MAXVAL(eh_nband)
 ABI_MALLOC(eh_bks_mask,(eh_mband,eh_nkibz,nsppol))
 ABI_MALLOC(eh_keep_ur ,(eh_mband,eh_nkibz,nsppol))
 eh_bks_mask=.TRUE.; eh_keep_ur =.TRUE.

 call wfd_init(Weh,Cryst,Pawtab,Psps,eh_keep_ur,Wsh%paral_kgb,npw_gamma,eh_mband,eh_nband,eh_nkibz,nsppol,&
&  eh_bks_mask,Wsh%nspden,nspinor,Wsh%ecutsm,Wsh%dilatmx,eh_istwfk,eh_kibz,Wsh%ngfft,kg_gamma,Wsh%nloalg,&
&  Wsh%prtvol,Wsh%pawprtvol,Wsh%comm)

 if (Weh%prtvol > 0) then
   call wfd_print(Weh,header="Shirley wavefunction descriptor")
 end if

 ABI_FREE(eh_keep_ur)
 ABI_FREE(eh_bks_mask)
 ABI_FREE(eh_nband)
 !
 ! =====================================================================
 ! ==== Rotate the input wavefunctions to get the optimal basis set ====
 ! =====================================================================

 fft_fact = one/Wsh%nfft
 ABI_MALLOC(ur1,(Wsh%nfft*nspinor))
 ABI_MALLOC(ur2,(Wsh%nfft*nspinor))
 ABI_MALLOC(ur12,(Wsh%nfft*nspinor))
 !
 write(msg,'(a,f12.1,a)')' Memory needed for storing eh_ur= ',two*gwpc*Wsh%nfft*nspinor*eh_size*b2Mb,' [Mb]'
 call wrtout(std_out,msg,'PERS')

 ABI_MALLOC(eh_ur,(Wsh%nfft*nspinor,eh_size))
 eh_ur = czero

 ABI_MALLOC(eh_ug,(npw_gamma*nspinor,eh_size))
 eh_ug = czero

 do midx=1,ovlp_size ! Loop over the single particle orbitals.
   band1  = Oeh%idx2bk(1,midx) ! TODO to be removed.
   band2  = Oeh%idx2bk(2,midx)

   call wfd_get_ur(Wsh,band1,k1,s1,ur1)
   call wfd_get_ur(Wsh,band2,k1,s1,ur2)
   ur12 = CONJG(ur1) * ur2   ! Use same convention as the one used in ovlp2_init.
   !
   ! Construct the new optimal basis set.
   do eh=1,eh_size
     eh_ur(:,eh) = eh_ur(:,eh) + Oeh%mat(midx,base+eh) * ur12
   end do
   !
 end do
 !
 ! NC: Normalize the basis set.
 do eh=1,eh_size
   norm1 = xdotc(Weh%nfft*nspinor,eh_ur(:,eh),1,eh_ur(:,eh),1) * fft_fact
   eh_ur(:,eh) = eh_ur(:,eh)/SQRT(norm1)
   !write(std_out,*)" eh_ur integrates to: ", xdotc(Weh%nfft*nspinor,eh_ur(:,eh),1,eh_ur(:,eh),1) * fft_fact
 end do
 !
 ! From the FFT mesh to the G-sphere.
 ABI_MALLOC(ehg,(npw_gamma*nspinor))
 ABI_MALLOC(gwpc_tmp,(Wsh%nfft*nspinor))
 !
 ! ============================================================================
 ! ==== Construct new optimal basis set in G-space and save results in Weh ====
 ! ============================================================================
 do eh=1,eh_size
   gwpc_tmp = eh_ur(:,eh)

#if 1
   ABI_MALLOC(dpc_tmp,(Wsh%nfft*nspinor))
   dpc_tmp  = eh_ur(:,eh)

   call fftbox_plan3_many(plan,Wsh%nspinor*ndat1,Wsh%ngfft(1:3),Wsh%ngfft(1:3),Wsh%ngfft(7),-1)
   call fftbox_execute(plan,dpc_tmp)
   !
   do ig=1,npw_gamma
     fft_idx = igfft0(ig)
     if (fft_idx/=0) then ! G-G0 belong to the FFT mesh.
       if (fft_idx>Wsh%nfft .or.fft_idx<0) then
         MSG_ERROR("fft_idx bug")
       end if
       ehg(ig)=dpc_tmp(fft_idx)
     else                 ! Set this component to zero.
       MSG_ERROR("fft_idx bug")
       ehg(ig)=czero
     end if
   end do
   ABI_FREE(dpc_tmp)
#else
   ! FIXME does not work anymore.
     ! TODO Check with the new version
   gbound_k => Weh%Kdata(k1)%gbound
   call fft_ur(npw_gamma,Wsh%nfft,nspinor,ndat1,Wsh%mgfft,Wsh%ngfft,istwf1,kg_gamma,gbound_k,gwpc_tmp,ehg)
#endif
   ! NC: Normalize the basis set using the eigenvalues of the overlap matrix.
   !if (usepaw==0) then
   !  ehg = ehg/SQRT(Oeh%eigene(base+eh))
   !end if
   !%call wfd_push_ug(Weh,eh,k1,s1,Cryst,ehg,update_ur=.FALSE.,update_cprj=.FALSE.)
   eh_ug(:,eh) = ehg
 end do
 !
 ! ======================================
 ! ==== Orthonormalize the basis set ====
 ! ======================================
 ABI_MALLOC(cf_ovlp,(eh_size,eh_size))

 call blas_cholesky_ortho(npw_gamma,eh_size,eh_ug,cf_ovlp)

 ABI_FREE(cf_ovlp)
 !
 ! Push data.
 do eh=1,eh_size
   call wfd_push_ug(Weh,eh,k1,s1,Cryst,eh_ug(:,eh),update_ur=.FALSE.,update_cprj=.FALSE.)
 end do

 ABI_FREE(gwpc_tmp)
 ABI_FREE(ehg)
 ABI_FREE(eh_ur)
 ABI_FREE(eh_ug)
 ABI_FREE(kg_gamma)
 ABI_FREE(igfft0)
 ABI_FREE(ur1)
 ABI_FREE(ur2)
 ABI_FREE(ur12)

 call cwtime(cpu,wall,gflops,"stop")
 write(std_out,*)"SHTIME: EH_Rotation, cpu_time: ",cpu,", wall_time:",wall
 !
 ! ==================================
 ! ==== sh2eh = <S_i* S_j| EH_a> ====
 ! ==================================
 call cwtime(cpu,wall,gflops,"start")
 ABI_MALLOC(sh2eh,(ovlp_size,eh_size))
 sh2eh=czero

 call wfd_change_ngfft(Weh,Cryst,Psps,Wsh%ngfft) ! Make sure the two set of wave are on the same mesh.
 ABI_MALLOC(ur1,(Weh%nfft*nspinor))

 do eh=1,eh_size
   call wfd_get_ur(Weh,eh,k1,s1,ur1)
   do ij=1,ovlp_size
     sh2eh(ij,eh) = DOT_PRODUCT(ur34_big(:,ij),ur1)
   end do
 end do
 ABI_FREE(ur1)

 ABI_FREE(ur34_big)

 call cwtime(cpu,wall,gflops,"stop")
 write(std_out,*)"SHTIME: EH_Projection, cpu_time: ",cpu,",wall_time:",wall
 !
 if (DEBUGME) then
   call wfd_test_ortho(Weh,Cryst,Pawtab,unit=std_out,mode_paral="COLL")
 end if
 !
 ! Deallocate memory.
 call ovlp_free(Oeh)

 DBG_EXIT("COLL")

end subroutine wfd_shirley_to_eh
!!***

END MODULE m_shirley
!!***
