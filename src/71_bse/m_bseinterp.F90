!!****m* ABINIT/m_bseinterp
!! NAME
!! m_bseinterp
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2020 ABINIT group (M.Giantomassi, Y. Gillet)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_bseinterp

 use defs_basis
 use m_abicore
 use m_bs_defs
 use m_xmpi
 use m_errors
 use m_nctk
 use m_haydock_io
 use m_linalg_interfaces
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_fstrings,          only : indent, strcat, sjoin, itoa
 use defs_datatypes,      only : pseudopotential_type
 use m_hide_blas,         only : xdotc
 use m_fft_mesh,          only : calc_ceigr
 use m_crystal,           only : crystal_t
 use m_bz_mesh,           only : kmesh_t
 use m_double_grid,       only : double_grid_t, get_kpt_from_indices_coarse
 use m_wfd,               only : wfd_t
 use m_pawtab,            only : pawtab_type

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_haydock/interpolator_t
!! NAME
!! interpolator_t
!!
!! FUNCTION
!!  Store the overlap matrix elements needed for the interpolation of the BSE Hamiltonian
!!
!! TODO
!!  Decide if we want to make the number of bands k-dependent.
!!
!! SOURCE

 type,public :: interpolator_t

    integer :: nvert=8
    ! Number of vertices for interpolation

    integer :: method
    ! Interpolation method (YG or Rohlfing & Louie or ...)

    integer :: mband_dense, mband_coarse
    ! Max number of bands dense and coarse

    integer :: nsppol
    ! Number of spin channels

    integer, allocatable :: corresp(:,:,:)
    ! corresp(max_nreh,nvert,spin)
    ! it_coarse, idiv -> it_coarse (idiv-th neighbour)

    real(dp),allocatable :: interp_factors(:,:)
    ! interp_factors(nvert,ndiv)
    ! index_in_fine_box -> k-point in Trans_interp

    complex(gwpc),allocatable :: overlaps(:,:,:,:,:)
    ! Overlaps between dense and coarse mesh
    ! overlaps(mband_coarse,mband_dense,ivertex_coarse,double_grid%nkpt_dense,spin)

    complex(dpc),allocatable :: btemp(:), ctemp(:)
    ! Temporary arrays for work

    ! Pointers to datatypes that are already in memory
    type(double_grid_t),pointer :: double_grid => null()
    ! Mapping between coarse and dense mesh

 end type interpolator_t

 public :: interpolator_init    ! Construct the object
 public :: interpolator_free    ! Free memory
 public :: interpolator_normalize ! Normalize the overlaps
 public :: int_alloc_work       ! Alloc temp memory
 public :: int_free_work        ! Free temp memory

!!***

!----------------------------------------------------------------------

CONTAINS  !=======================================================================
!!***

!!****f* m_bseinterp/interpolator_init
!! NAME
!! interpolator_init
!!
!! FUNCTION
!! Construct the interpolator object
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_hexc
!!
!! CHILDREN
!!
!! SOURCE

subroutine interpolator_init(interpolator, double_grid, Wfd_dense, Wfd_coarse, &
&    Kmesh_dense, Kmesh_coarse, BSp, Cryst, Psps, Pawtab, method)

!Arguments ---------------------------
!scalars
 integer,intent(in) :: method
 type(interpolator_t),intent(inout) :: interpolator
 type(double_grid_t),intent(in),target :: double_grid
 type(wfd_t),intent(inout) :: Wfd_dense, Wfd_coarse
 type(kmesh_t),intent(in) :: Kmesh_dense, Kmesh_coarse
 type(excparam),intent(in) :: BSp
 type(crystal_t),intent(in) :: Cryst
 type(pseudopotential_type),intent(in) :: Psps
!arrays
 type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd_coarse%usepaw)

!Local variables ---------------------
!scalars
 integer :: nsppol, nvert
 integer :: maxnreh, nreh1, nreh2
 integer :: mbandc, mbandd, nbzd
 real(dp),parameter :: threshold = 0.1_dp
!arrays
 character(len=500) :: msg

!*****************************************************************************

 ABI_CHECK(Wfd_coarse%usepaw==0, "PAW not yet supported")
 ABI_CHECK(BSp%nsppol==1, "nsppol != 1 not yet implemented")
 ABI_CHECK(Wfd_coarse%nspinor==1, "nspinor != 1 not supported")

 ABI_UNUSED(Pawtab(1)%basis_size)
 !paw_overlap(cprj1,cprj2,typat,pawtab,spinor_comm) result(onsite)

 interpolator%double_grid => double_grid
 interpolator%mband_dense = Wfd_dense%mband
 interpolator%mband_coarse = Wfd_coarse%mband
 interpolator%method = method
 interpolator%nsppol = BSp%nsppol

 SELECT CASE(method)
 CASE (BSE_INTERP_YG)
   nvert = 8
 CASE (BSE_INTERP_RL2)
   nvert = 2
 CASE (BSE_INTERP_RL)
   nvert = 1
 CASE DEFAULT
   write(msg,'(a,i0)') "Wrong interpolation method: ",method
   MSG_ERROR(msg)
 END SELECT

 interpolator%nvert = nvert

 mbandc = interpolator%mband_coarse
 mbandd = interpolator%mband_dense
 nbzd = double_grid%nbz_dense
 nsppol = interpolator%nsppol
 ABI_MALLOC(interpolator%overlaps,(mbandc,mbandd,nvert,nbzd,nsppol))

 call int_compute_overlaps(interpolator,double_grid, Wfd_dense, Wfd_coarse, Kmesh_dense, &
&   Kmesh_coarse, BSp, Cryst, Psps, Pawtab)

 ABI_MALLOC(interpolator%interp_factors,(nvert,double_grid%ndiv))

 call int_preprocess_tables(interpolator,double_grid)

 nreh1 = BSp%nreh(1)
 nreh2 = nreh1; if(BSp%nsppol == 2) nreh2 = BSp%nreh(2)
 maxnreh = MAX(nreh1,nreh2)

 ABI_MALLOC(interpolator%corresp,(maxnreh,interpolator%nvert,interpolator%nsppol))

 call int_compute_corresp(interpolator,BSp,double_grid)

end subroutine interpolator_init
!!***

!-------------------------------------------------------------------

!!****f* m_bseinterp/int_alloc_work
!! NAME
!! int_alloc_work
!!
!! FUNCTION
!! Allocate temporary arrays
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_hexc
!!
!! CHILDREN
!!
!! SOURCE

subroutine int_alloc_work(interpolator, work_size)

!Arguments ---------------------------
!scalars
 integer,intent(in) :: work_size
 type(interpolator_t),intent(inout) :: interpolator

!*****************************************************************************

 ABI_MALLOC(interpolator%btemp,(work_size))
 ABI_MALLOC(interpolator%ctemp,(work_size))

end subroutine int_alloc_work
!!***

!-------------------------------------------------------------------

!!****f* m_bseinterp/int_free_work
!! NAME
!! int_free_work
!!
!! FUNCTION
!! Deallocate temporary arrays
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_hexc
!!
!! CHILDREN
!!
!! SOURCE

subroutine int_free_work(interpolator)

!Arguments ---------------------------
!scalars
 type(interpolator_t),intent(inout) :: interpolator

!*****************************************************************************

 if( allocated(interpolator%btemp)) then
   ABI_FREE(interpolator%btemp)
 end if
 if( allocated(interpolator%ctemp)) then
   ABI_FREE(interpolator%ctemp)
 end if

end subroutine int_free_work
!!***

!-------------------------------------------------------------------

!!****f* m_bseinterp/int_compute_overlaps
!! NAME
!! int_compute_overlaps
!!
!! FUNCTION
!! Compute the overlaps prefactors
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_bseinterp
!!
!! CHILDREN
!!
!! SOURCE

subroutine int_compute_overlaps(interpolator, double_grid, Wfd_dense, Wfd_coarse, &
&   Kmesh_dense, Kmesh_coarse, BSp, Cryst, Psps, Pawtab)

!Arguments ---------------------------
!scalars
 type(interpolator_t),intent(inout) :: interpolator
 type(double_grid_t),intent(in),target :: double_grid
 type(wfd_t),intent(inout) :: Wfd_dense, Wfd_coarse
 type(kmesh_t),intent(in) :: Kmesh_dense, Kmesh_coarse
 type(excparam),intent(in) :: BSp
 type(crystal_t),intent(in) :: Cryst
 type(pseudopotential_type),intent(in) :: Psps
!arrays
 type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd_coarse%usepaw)

!Local variables ---------------------
!scalars
 integer :: nprocs, my_rank
 integer :: ierr
 integer :: nfft, nspinor, nsppol, nvert
 integer :: ib_coarse, ib_dense
 integer :: ik_coarse, ik_dense
 integer :: spin
 integer :: iorder
 integer :: ivertex, ix, iy, iz
 integer :: bstart, bstop
 real(dp),parameter :: threshold = 0.1_dp
 complex(gwpc) :: ovlp
!arrays
 integer :: curindices_dense(6), curindices_coarse(3)
 integer :: neighbour(3)
 integer :: g0(3),g01(3),diffg0(3)
 complex(gwpc),allocatable :: ur_coarse(:),ur_dense(:)
 complex(gwpc),allocatable :: ceigr(:)
!arrays

!*****************************************************************************

 nprocs = xmpi_comm_size(Wfd_coarse%comm)
 my_rank = xmpi_comm_rank(Wfd_coarse%comm)

 ABI_UNUSED(Pawtab(1)%basis_size)

 ! Ensure Wfd and Wfd_coarse use the same FFT mesh.
 call wfd_dense%change_ngfft(Cryst,Psps,Wfd_coarse%ngfft)
 nfft = Wfd_coarse%nfft
 nspinor = Wfd_coarse%nspinor
 nsppol = Bsp%nsppol
 nvert = interpolator%nvert

 ! Allocate workspace for wavefunctions in real space.
 ABI_MALLOC(ur_coarse,(nfft*nspinor))
 ABI_MALLOC(ur_dense,(nfft*nspinor))
 ABI_MALLOC(ceigr,(nfft*nspinor))

 interpolator%overlaps = czero

 ! TODO
 ! 1) Choose whether we want to compute only dvv, dcc or all dbb
 ! 2) Check the ordering of the loops
 ! 3) Improve vertex -> neighbour (in double_grid ?)
 do spin = 1,nsppol
   do ik_dense = 1,double_grid%nbz_dense

     ! MPI parallelization
     ! We assume that each node owns in memory the full set of wavefunctions
     ! both coarse and dense k-mesh and both spins.
     if (mod(ik_dense, nprocs) /= my_rank) cycle

     ! From ik_dense -> indices_dense
     iorder = double_grid%iktoint_dense(ik_dense)
     g01 = double_grid%g0_dense(:,iorder)
     curindices_dense = double_grid%indices_dense(:,iorder)

     do ivertex = 1,nvert

       ! From vertex to neighbour
       ! TODO improve this part + permit to choose other neighbour (e.g. nearest neighbour for RL)
       if(nvert > 1) then
         ix = (ivertex-1)/4
         iy = (ivertex-ix*4-1)/2
         iz = (ivertex-ix*4-iy*2-1)
       else
         ix = (BSp%rl_nb-1)/4
         iy = (BSp%rl_nb-ix*4-1)/2
         iz = (BSp%rl_nb-ix*4-iy*2-1)
       end if

       neighbour = [ix,iy,iz]

       ! From indices_dense -> indices_coarse
       curindices_coarse = curindices_dense(1:3) + neighbour(:)

       ! From indices_coarse -> ik_ibz in the coarse mesh
       call get_kpt_from_indices_coarse(curindices_coarse,double_grid%maxcomp_coarse,&
&        double_grid%inttoik_coarse,double_grid%g0_coarse,double_grid%nbz_closedcoarse,ik_coarse,g0)

       ! Take into account a possible umklapp between k_dense and k_coarse
       diffg0 = g0 - g01

       if (ANY(diffg0/=0)) then
         ! WARNING works only with nspinor = 1 !!!
         call calc_ceigr(diffg0,nfft,nspinor,Wfd_coarse%ngfft,ceigr)
       end if

       do ib_dense = BSp%lomo_spin(spin), BSp%humo_spin(spin)
         ! ur(ib_dense, ik_dense)
         call wfd_dense%sym_ur(Cryst,Kmesh_dense,ib_dense,ik_dense,spin,ur_dense)

         if (ANY(diffg0/=0)) then
           !ur_kbz = ur_kbz*e(ig0r)
           ur_dense(:) = ur_dense(:)*ceigr(:)
         end if

         ! Uncomment for the complete overlap
         !bstart = BSp%lomo_spin(spin); bstop = BSp%humo_spin(spin)

         ! Compute only dvv or dcc
         if (ib_dense <= BSp%homo_spin(spin)) then
           ! if ib_dense is a valence band => loop on valence bands
           bstart = BSp%lomo_spin(spin); bstop = BSp%homo_spin(spin)
         else
           ! if ib_dense is a conduction band => loop on conduction bands
           bstart = BSp%lumo_spin(spin); bstop = BSp%humo_spin(spin)
         end if

         do ib_coarse = bstart, bstop
           ! ur(ib_coarse, ik_coarse)
           call wfd_coarse%sym_ur(Cryst,Kmesh_coarse,ib_coarse,ik_coarse,spin,ur_coarse)

           ! ovlp = < u_{ib_coarse,ik_coarse} | u_{ib_dense,ik_dense} >
           ovlp =  xdotc(nfft,ur_coarse,1,ur_dense,1)/nfft

           ! Filter too low values
           if (ABS(ovlp) < threshold) ovlp = czero

           interpolator%overlaps(ib_coarse,ib_dense,ivertex,ik_dense,spin) = ovlp
         end do ! ib_coarse

         !DBYG
         ! write(std_out,*) "nb = ",neighbour
         ! write(std_out,*) "(i1,i2,i3,j1,j2,j3) = ",curindices_dense
         ! write(std_out,*) "ib = ",ib_dense
         ! write(std_out,*) "Sum of dbb = ",REAL(SUM(GWPC_CONJG(interpolator%overlaps(:,ib_dense,ivertex,ik_dense,spin))*interpolator%overlaps(:,ib_dense,ivertex,ik_dense,spin))); call flush(std_out)
         !ENDDBYG

       end do ! ib_dense
     end do ! ivertex
   end do ! ik_dense
 end do ! spin

 ABI_FREE(ur_coarse)
 ABI_FREE(ur_dense)
 ABI_FREE(ceigr)

 ! Gather results on each node.
 call xmpi_sum(interpolator%overlaps,Wfd_coarse%comm,ierr)

end subroutine int_compute_overlaps
!!***

!----------------------------------------------------------------------

!!****f* m_bseinterp/int_preprocess_tables
!! NAME
!! int_preprocess_tables
!!
!! FUNCTION
!! Pre-process tables to improve interpolation technique
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_bseinterp
!!
!! CHILDREN
!!
!! SOURCE

subroutine int_preprocess_tables(interpolator,double_grid)

!Argument ------------------------------------
!scalars
 type(interpolator_t),intent(inout) :: interpolator
 type(double_grid_t),intent(in) :: double_grid
!arrays

!Local variables -----------------------------
!scalars
 integer :: iorder,ik_dense,ik_coarse
 integer :: ix,iy,iz,ineighbour,curdim, curj
 real(dp) :: interp_factor
!arrays
 integer :: allxyz(3),curindices_dense(6)
 integer,allocatable :: curindex(:)

!*********************************************
 ABI_MALLOC(curindex,(double_grid%nbz_coarse))
 curindex = 1

 interpolator%interp_factors = zero

 do ik_dense = 1,double_grid%nbz_dense

   ! From ik_ibz in the dense mesh -> indices_dense
   iorder = double_grid%iktoint_dense(ik_dense)
   !g01 = double_grid%g0_dense(:,iorder)

   ! From indices_dense -> indices_coarse
   curindices_dense = double_grid%indices_dense(:,iorder)

   ik_coarse = double_grid%dense_to_coarse(ik_dense)

   ! Compute multi-linear interpolation factors
   ! Loop over the neighbours
   do ineighbour = 1,interpolator%nvert
     !TODO helper function from [ix,iy,iz] -> ineighbour and vice versa
     ix = (ineighbour-1)/4
     iy = (ineighbour-ix*4-1)/2
     iz = (ineighbour-ix*4-iy*2-1)
     allxyz = [ix,iy,iz]
     interp_factor = one
     do curdim = 1,3
       if (interpolator%method == BSE_INTERP_RL) then
         cycle
       else if(interpolator%method == BSE_INTERP_RL2) then
         if (curdim /= 3) cycle
       end if
       curj = curindices_dense(3+curdim)
       interp_factor = interp_factor*((allxyz(curdim)*(curj*1.0/double_grid%kmult(curdim)))&
&                               +((1-allxyz(curdim))*(1-(curj*1.0/double_grid%kmult(curdim)))))
     end do
     interpolator%interp_factors(ineighbour,curindex(ik_coarse)) = interp_factor
   end do

   curindex(ik_coarse) = curindex(ik_coarse) + 1
 end do

 ABI_FREE(curindex)

end subroutine int_preprocess_tables
!!***

!-------------------------------------------------------------------

!!****f* m_haydock/int_compute_corresp
!! NAME
!! int_compute_corresp
!!
!! FUNCTION
!!
!! INPUTS
!! BSp<type(excparam)=The parameter for the Bethe-Salpeter run.
!! grid <double_grid_t>=Correspondence between coarse and fine k-grid
!! spin=Spin index.
!!
!! OUTPUT
!! corresp(Bsp%nreh(spin),8)= Correspondence between a transition on the
!!   coarse mesh and its i-th neighbour for i in [1,2,..,8].
!!
!! TODO:
!!  Some operations are faster if we allocate with shape (8,nreh(spin))
!!
!! PARENTS
!!      m_bseinterp
!!
!! CHILDREN
!!
!! SOURCE

subroutine int_compute_corresp(interpolator,BSp,double_grid)

!Arguments ------------------------------------
!scalars
 type(excparam),intent(in) :: BSp
 type(double_grid_t),intent(in) :: double_grid
 type(interpolator_t),intent(inout) :: interpolator

!Local variables ------------------------------
!scalars
 integer :: spin
 integer :: itt,ik_dense,ik_coarse,iorder,it_coarse
 integer :: ic,iv,ik_coarse0,it_coarse0,iovlp,ix,iy,iz
!arrays
 integer :: curindices_dense(6),curindices_coarse(3),g0(3),g01(3),neighbour(3)

!************************************************************************

 do spin=1,interpolator%nsppol
   do itt=1,BSp%nreh_interp(spin)
     ! From dense itt -> ik_dense, ic, iv
     ik_dense = BSp%Trans_interp(itt,spin)%k
     ic = BSp%Trans_interp(itt,spin)%c
     iv = BSp%Trans_interp(itt,spin)%v

     ! From ik_dense -> indices_dense
     iorder = double_grid%iktoint_dense(ik_dense)
     g01 = double_grid%g0_dense(:,iorder)

     ! Index of the k-point in the coarse mesh.
     ik_coarse0 = double_grid%dense_to_coarse(ik_dense)
     it_coarse0 = BSp%vcks2t(iv,ic,ik_coarse0,spin)

     ! From indices_dense -> indices_coarse
     curindices_dense = double_grid%indices_dense(:,iorder)

     ! Loop over the 8 neighbors.
     do iovlp = 1,interpolator%nvert

       !TODO : helper function from [ix,iy,iz] -> iovlp and vice versa
       if(interpolator%nvert > 1) then
         ix = (iovlp-1)/4
         iy = (iovlp-ix*4-1)/2
         iz = (iovlp-ix*4-iy*2-1)
       else
         ix = (BSp%rl_nb-1)/4
         iy = (BSp%rl_nb-ix*4-1)/2
         iz = (BSp%rl_nb-ix*4-iy*2-1)
       end if
       neighbour = [ix,iy,iz]

       curindices_coarse = curindices_dense(1:3) + neighbour(:)

       ! From indices_coarse -> ik_ibz in the coarse mesh
       call get_kpt_from_indices_coarse(curindices_coarse,double_grid%maxcomp_coarse,&
&        double_grid%inttoik_coarse,double_grid%g0_coarse,double_grid%nbz_closedcoarse,ik_coarse,g0)

       ! From ik_coarse, ic, iv to it_coarse
       it_coarse = BSp%vcks2t(iv,ic,ik_coarse,spin)

       interpolator%corresp(it_coarse0,iovlp,spin) = it_coarse
     end do
   end do ! itt
 end do

end subroutine int_compute_corresp
!!***

!----------------------------------------------------------------------

!!****f* m_bseinterp/interpolator_normalize
!! NAME
!! interpolator_normalize
!!
!! FUNCTION
!! Normalize the overlaps so that \sum_{ib} | d_{kk'}^{b,ib} | ^2 = 1
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_hexc
!!
!! CHILDREN
!!
!! SOURCE

subroutine interpolator_normalize(interpolator)

!Arguments ---------------------------
!scalars
 type(interpolator_t),intent(inout) :: interpolator

!Local variables ---------------------
!scalars
 integer :: spin, ivertex, ib_dense, ik_dense
 complex(gwpc) :: sum_ovlp
!arrays
 complex(gwpc),allocatable :: overlaps(:)

!*****************************************************************************

 ABI_MALLOC(overlaps,(interpolator%mband_coarse))
 do spin = 1, interpolator%nsppol
   do ivertex = 1, interpolator%nvert
     do ib_dense = 1, interpolator%mband_dense
       do ik_dense = 1, interpolator%double_grid%nbz_dense
         overlaps(:) = interpolator%overlaps(:,ib_dense,ivertex,ik_dense,spin)
         sum_ovlp = SQRT(REAL(SUM(GWPC_CONJG(overlaps(:))*overlaps(:))))
         if (ABS(sum_ovlp) > tol6) then
           overlaps(:) = overlaps(:)/sum_ovlp
         else
           overlaps(:) = czero
         end if
         interpolator%overlaps(:,ib_dense,ivertex,ik_dense,spin) = overlaps(:)
       end do
     end do
   end do
 end do
 ABI_FREE(overlaps)

end subroutine interpolator_normalize
!!***

!-------------------------------------------------------------------

!!****f* m_bseinterp/interpolator_free
!! NAME
!! interpolator_free
!!
!! FUNCTION
!! Destroy the interpolator object in memory
!!
!! INPUTS
!!
!!
!! OUTPUT
!!
!!
!! PARENTS
!!      m_hexc
!!
!! CHILDREN
!!
!! SOURCE

subroutine interpolator_free(interpolator)

!Arguments ---------------------------
 type(interpolator_t),intent(inout) :: interpolator

!*****************************************************************************

 if( allocated(interpolator%overlaps) ) then
   ABI_FREE(interpolator%overlaps)
 end if

 if (allocated(interpolator%corresp)) then
   ABI_FREE(interpolator%corresp)
 end if

 if (allocated(interpolator%interp_factors)) then
   ABI_FREE(interpolator%interp_factors)
 end if

 if( associated(interpolator%double_grid) ) then
   nullify(interpolator%double_grid)
 end if

end subroutine interpolator_free
!!***

!-------------------------------------------------------------------

end module m_bseinterp
!!***
