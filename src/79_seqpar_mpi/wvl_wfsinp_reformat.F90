!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_wfsinp_reformat
!! NAME
!! wvl_wfsinp_reformat
!!
!! FUNCTION
!! This method allocates and initialises wavefunctions with values from disk.
!! See wvl_wfsinp_scratch() or wvl_wfsinp_reformat() from other initialisation
!! routines.
!! 
!! When initialised from scratch or from disk, wvl%wfs%ks%[h]psi comes unallocated
!! and will be allocated inside this routine.
!! When initialised from memory (reformating), wvl%wfs%ks%[h]psi will be reallocated.
!! The projectors are also recomputed.
!!
!! The scalar arrays should be reallocated using dtset%nfft after a call to
!! this routine.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2016 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      copy_old_wavefunctions,deallocate_wfd,first_orthon
!!      local_potential_dimensions,nullify_gaussian_basis,psolver_kernel
!!      reformatmywaves,wrtout,wvl_denspot_free,wvl_denspot_set
!!      wvl_projectors_free,wvl_projectors_set,wvl_setboxgeometry,wvl_setngfft
!!      wvl_wfs_free,wvl_wfs_lr_copy,wvl_wfs_set,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine wvl_wfsinp_reformat(dtset, mpi_enreg, psps, rprimd, wvl, xred, xred_old)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_profiling_abi
 use m_errors
 use m_xmpi

#if defined HAVE_DFT_BIGDFT
 use BigDFT_API, only : copy_old_wavefunctions, reformatmywaves, first_orthon, &
& deallocate_wfd, wavefunctions_descriptors, deallocate_lr, &
& local_potential_dimensions, copy_coulomb_operator, &
& deallocate_coulomb_operator, nullify_gaussian_basis
 use dynamic_memory 
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_wfsinp_reformat'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
 use interfaces_43_wvl_wrappers
 use interfaces_62_poisson
!End of the abilint section

  implicit none

!Arguments ------------------------------------
  type(dataset_type), intent(inout)      :: dtset
  type(MPI_type), intent(inout)          :: mpi_enreg
  type(pseudopotential_type), intent(in) :: psps
  type(wvl_data), intent(inout)          :: wvl
  real(dp), intent(inout)                :: rprimd(3,3)
  real(dp), intent(inout)                :: xred_old(3, dtset%natom)
  real(dp), intent(inout)                :: xred(3, dtset%natom)

!Local variables-------------------------------
#if defined HAVE_DFT_BIGDFT
  integer                  :: itypat
  integer                  :: nSize_old(3)
  real(dp)                 :: hgrid_old(3)
  real(dp), allocatable    :: xcart(:,:), xcart_old(:,:)
  real(dp), pointer        :: psi_old(:), eigen_old(:)
  integer :: comm,me,nproc,icoulomb
  type(coulomb_operator)::kernel
  type(wavefunctions_descriptors) :: keys_old
  character(len=500)       :: message
#if defined DEBUG_MODE
 !integer, parameter :: ndebug = 5  !5 will not work for wavelets compiling with debug=naughty
 integer,parameter :: ndebug = 0
#else
 integer, parameter :: ndebug = 0
#endif
#endif

! *********************************************************************

#if defined HAVE_DFT_BIGDFT

 write(message, '(a,a)' ) ch10,&
& ' wvl_wfsinp_reformat: reformat the wavefunctions.'
 call wrtout(std_out, message, 'COLL')

 comm=mpi_enreg%comm_wvl
 me=xmpi_comm_rank(comm)
 nproc=xmpi_comm_size(comm)
!Convert input xred_old (reduced coordinates) to xcart_old (cartesian)
 ABI_ALLOCATE(xcart_old,(3, dtset%natom))
 call xred2xcart(dtset%natom, rprimd, xcart_old, xred_old)

!Copy current to old.
 ABI_ALLOCATE(eigen_old,(wvl%wfs%ks%orbs%norb))
 eigen_old = wvl%wfs%ks%orbs%eval
 hgrid_old = wvl%descr%h
 call copy_old_wavefunctions(nproc, wvl%wfs%ks%orbs, &
& wvl%descr%Glr%d%n1, wvl%descr%Glr%d%n2, wvl%descr%Glr%d%n3, &
& wvl%wfs%ks%lzd%Glr%wfd, wvl%wfs%ks%psi, nSize_old(1), nSize_old(2), nSize_old(3), &
& keys_old, psi_old)
!Patch because copy_old_wavefunctions() free wvl%wfs%ks%lzd%Glr%wfd but don't nullify it.
 nullify(wvl%wfs%ks%lzd%glr%wfd%keyglob)
 nullify(wvl%wfs%ks%lzd%glr%wfd%keygloc)
 nullify(wvl%wfs%ks%lzd%glr%wfd%keyvloc)
 nullify(wvl%wfs%ks%lzd%glr%wfd%keyvglob)

!We deallocate the previous projectors.
 call wvl_projectors_free(wvl%projectors)

!Deallocate old wavefunctions
 call wvl_wfs_free(wvl%wfs)

!Deallocate old denspot
 call wvl_denspot_free(wvl%den)

!We change the box geometry.
 call wvl_setBoxGeometry(dtset%prtvol, psps%gth_params%radii_cf, rprimd, xred, &
& wvl%descr, dtset%wvl_crmult, dtset%wvl_frmult)
 call wvl_denspot_set(wvl%den, psps%gth_params, dtset%ixc, dtset%natom, dtset%nsppol, &
& rprimd, wvl%descr, dtset%wvl_crmult, dtset%wvl_frmult, comm, xred)
 if (wvl%descr%atoms%astruct%geocode == "F") then
   icoulomb = 1
 else if (wvl%descr%atoms%astruct%geocode == "S") then
   icoulomb = 2
 else
   icoulomb = 0
 end if
!calculation of the Poisson kernel anticipated to reduce memory peak for small systems
 call psolver_kernel( wvl%den%denspot%dpbox%hgrids, 2, icoulomb, me, kernel, &
& comm, wvl%den%denspot%dpbox%ndims, nproc, dtset%nscforder)
 ! Shallow copy of the kernel (still owned by ABINIT).
 wvl%den%denspot%pkernel = kernel
 wvl%den%denspot%pkernelseq = kernel
!Associate the denspot distribution into mpi_enreg.
 mpi_enreg%nscatterarr => wvl%den%denspot%dpbox%nscatterarr
 mpi_enreg%ngatherarr => wvl%den%denspot%dpbox%ngatherarr
 mpi_enreg%ngfft3_ionic = wvl%den%denspot%dpbox%n3pi
 call wvl_setngfft(me, dtset%mgfft, dtset%nfft, &
& dtset%ngfft, nproc, wvl%den%denspot%dpbox%ndims(1), &
& wvl%den%denspot%dpbox%ndims(2), &
& wvl%den%denspot%dpbox%ndims(3),wvl%den%denspot%dpbox%n3d)

!We copy the geometry structure.
 call wvl_wfs_lr_copy(wvl%wfs, wvl%descr)
!Reallocate them with new size.
 call wvl_wfs_set(dtset%strprecon,dtset%spinmagntarget, dtset%kpt, me, dtset%natom, sum(dtset%nband), dtset%nkpt, &
& nproc, dtset%nspinor, dtset%nsppol, dtset%nwfshist, dtset%occ_orig, psps, rprimd, &
& wvl%wfs, dtset%wtk, wvl%descr, dtset%wvl_crmult, dtset%wvl_frmult, xred)

!Recopy old eval for precond.
 wvl%wfs%ks%orbs%eval = eigen_old
 ABI_DEALLOCATE(eigen_old)

!We allocate psi.
!ABI_ALLOCATE(wvl%wfs%ks%psi,( max(wvl%wfs%ks%orbs%npsidim_comp,wvl%wfs%ks%orbs%npsidim_orbs)+ndebug) )
 wvl%wfs%ks%psi=f_malloc_ptr(max(wvl%wfs%ks%orbs%npsidim_comp,wvl%wfs%ks%orbs%npsidim_orbs)+ndebug,id='psi')
 write(message, '(a,a,a,a,I0)' ) ch10, &
& ' wvl_wfsinp_reformat: allocate wavefunctions,', ch10, &
& '  size of the compressed array per proc: ', &
& product(shape(wvl%wfs%ks%psi))
 call wrtout(std_out,message,'COLL')

!Convert input xred (reduced coordinates) to xcart (cartesian)
 ABI_ALLOCATE(xcart,(3, dtset%natom))
 call xred2xcart(dtset%natom, rprimd, xcart, xred)

!We transfer the old wavefunctions to the new ones.
 call reformatmywaves(me, wvl%wfs%ks%orbs, wvl%descr%atoms, &
& hgrid_old(1), hgrid_old(2), hgrid_old(3), nSize_old(1), nSize_old(2), &
& nSize_old(3), xcart_old, keys_old, psi_old, wvl%descr%h(1), wvl%descr%h(2), &
& wvl%descr%h(3), wvl%descr%Glr%d%n1, wvl%descr%Glr%d%n2, wvl%descr%Glr%d%n3, xcart, &
& wvl%wfs%ks%lzd%Glr%wfd, wvl%wfs%ks%psi)
 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(xcart_old)

!We free the old descriptors and arrays.
 ABI_DEALLOCATE(psi_old)
 call deallocate_wfd(keys_old)

 call local_potential_dimensions(me,wvl%wfs%ks%lzd,wvl%wfs%ks%orbs,wvl%den%denspot%xc,&
& wvl%den%denspot%dpbox%ngatherarr(0,1))

!it seems that the table "wvl%projectors%G" is no more used
!but it's not allocated -> fortran runtime error 
#if defined HAVE_DFT_BIGDFT
 ABI_DATATYPE_ALLOCATE(wvl%projectors%G,(dtset%ntypat))
 do itypat=1,dtset%ntypat
   call nullify_gaussian_basis(wvl%projectors%G(itypat))
 end do
#endif

!Reallocate projectors for the new positions.
 call wvl_projectors_set(me, dtset%natom, wvl%projectors, psps, rprimd, &
& wvl%wfs, wvl%descr, dtset%wvl_frmult, xred)

!Orthogonilise new wavefunctions.
 call first_orthon(me, nproc, wvl%wfs%ks%orbs, wvl%wfs%ks%lzd, wvl%wfs%ks%comms, &
& wvl%wfs%ks%psi, wvl%wfs%ks%hpsi, wvl%wfs%ks%psit, wvl%wfs%ks%orthpar,wvl%descr%paw)

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) dtset%nstep,mpi_enreg%me,psps%npsp,wvl%wfs%ks,rprimd(1,1),&
& xred_old(1,1),xred(1,1)
#endif

end subroutine wvl_wfsinp_reformat
!!***
