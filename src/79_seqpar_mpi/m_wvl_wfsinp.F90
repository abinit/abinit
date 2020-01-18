!!****m* ABINIT/m_wvl_wfsinp
!! NAME
!!  m_wvl_wfsinp
!!
!! FUNCTION
!!  Routines to initialize (wavelet) wavefunctions
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (DC)
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

module m_wvl_wfsinp

 use defs_basis
 use defs_wvltypes
 use m_wffile
 use m_abicore
 use m_errors
 use m_xmpi
 use m_hdr
 use m_dtset

 use defs_datatypes, only : pseudopotential_type
 use defs_abitypes,  only : MPI_type
 use m_geometry,  only : xred2xcart
 use m_abi2big,   only : wvl_occ_abi2big, wvl_occopt_abi2big, wvl_setngfft, wvl_setboxgeometry
 use m_psolver,   only : psolver_kernel
 use m_wvl_rwwf,  only : wvl_read
 use m_mklocl_realspace, only : mklocl_wavelets
 use m_wvl_wfs,          only : wvl_wfs_set, wvl_wfs_free, wvl_wfs_lr_copy
 use m_wvl_denspot,      only : wvl_denspot_set, wvl_denspot_free
 use m_wvl_projectors,   only : wvl_projectors_set, wvl_projectors_free

 implicit none

 private
!!***

 public :: wvl_wfsinp_disk
 public :: wvl_wfsinp_reformat
 public :: wvl_wfsinp_scratch
!!***

contains
!!***

!!****f* ABINIT/wvl_wfsinp_disk
!! NAME
!! wvl_wfsinp_disk
!!
!! FUNCTION
!! This method allocates and initialises wavefunctions with values from disk.
!! See wvl_wfsinp_scratch() or wvl_wfsinp_reformat() from other initialisation
!! routines.
!!
!! When initialised from scratch or from disk, wvl%wfs%[h]psi comes unallocated
!! and will be allocated inside this routine.
!! When initialised from memory (reformating), wvl%wfs%[h]psi will be reallocated.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=input variables.
!!  hdr0 <type(hdr_type)>=the header of wf, den and pot files (read from restart)
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  mpi_enreg=informations about MPI parallelization
!!  option=1 for reading a file following ABINIT format, -1 for a BigDFT format.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  wff <type(wffile_type)>= structure with informations on wf file.
!!  xred(3,natom)=reduced dimensionless atomic coordinates.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  wfs <type(wvl_projector_type)>=wavefunctions informations for wavelets.
!!
!! PARENTS
!!      inwffil
!!
!! CHILDREN
!!      first_orthon,wrtout,wvl_occ_abi2big,wvl_read
!!
!! SOURCE

subroutine wvl_wfsinp_disk(dtset, hdr0, hdr, mpi_enreg, occ, option, rprimd, wff, wfs, wvl, xred)

#if defined HAVE_BIGDFT
 use BigDFT_API, only : first_orthon,sumrho,communicate_density,plot_density
 use dynamic_memory
#endif

!Arguments -------------------------------
  !scalars
  integer, intent(in)                       :: option
  type(dataset_type), intent(in)            :: dtset
  type(hdr_type), intent(in)                :: hdr0
  type(hdr_type), intent(in)                :: hdr
  type(MPI_type), intent(in)                :: mpi_enreg
  type(wffile_type), intent(in)             :: wff
  type(wvl_wf_type), intent(inout)          :: wfs
  type(wvl_internal_type), intent(inout)       :: wvl
  !type(wvl_denspot_type), intent(inout)       :: wvl_den
  !arrays
  real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), intent(in)                      :: rprimd(3, 3)
  real(dp), intent(in)                      :: xred(3, dtset%natom)

!Local variables-------------------------------
#if defined HAVE_BIGDFT
 logical :: wvlbigdft=.false.
 integer :: comm,me,nproc
 character(len = 500)  :: message
 !real(dp), allocatable :: xcart(:,:)
#if defined DEBUG_MODE
 !integer, parameter :: ndebug = 5  !5 will not work for wavelets compiling with debug=naughty
 integer,parameter :: ndebug = 0
#else
 integer, parameter :: ndebug = 0
#endif
#endif

! *********************************************************************

#if defined HAVE_BIGDFT

 write(message, '(a,a)' ) ch10,' wvl_wfsinp_disk: wavefunction initialisation.'
 call wrtout(std_out,message,'COLL')

!If usewvl: wvlbigdft indicates that the BigDFT workflow will be followed
 wvlbigdft=(dtset%usewvl==1.and.dtset%wvl_bigdft_comp==1)

 comm=mpi_enreg%comm_wvl
 me=xmpi_comm_rank(comm)
 nproc=xmpi_comm_size(comm)
!We allocate psi.
!ABI_ALLOCATE(wfs%ks%psi,( max(wfs%ks%orbs%npsidim_comp,wfs%ks%orbs%npsidim_orbs)+ndebug) )
 wfs%ks%psi=f_malloc_ptr(max(wfs%ks%orbs%npsidim_comp,wfs%ks%orbs%npsidim_orbs)+ndebug,id='psi')

 write(message, '(a,a,a,a,I0)' ) ch10, &
& ' wvl_wfsinp_disk: allocate wavefunctions,', ch10, &
& '  size of the compressed array per proc: ', &
& product(shape(wfs%ks%psi))
 call wrtout(std_out,message,'COLL')

 call wvl_read(dtset, hdr0, hdr, mpi_enreg, option, rprimd, wff, wfs, wvl, xred)

!We orthogonalise,only for NC.
 if(wvl%paw%usepaw==0 .and. wvlbigdft) then
   call first_orthon(me, nproc, wfs%ks%orbs, wfs%ks%lzd, wfs%ks%comms, &
&   wfs%ks%psi, wfs%ks%hpsi, wfs%ks%psit, wfs%ks%orthpar,wvl%paw)
 else
!  ABI_ALLOCATE(wfs%ks%hpsi,(max(wfs%ks%orbs%npsidim_orbs,wfs%ks%orbs%npsidim_comp)))
   wfs%ks%hpsi=f_malloc_ptr(max(wfs%ks%orbs%npsidim_orbs,wfs%ks%orbs%npsidim_comp),id='hpsi')
   if(wvl%paw%usepaw==1) then
     ABI_ALLOCATE(wvl%paw%spsi,(max(wfs%ks%orbs%npsidim_orbs,wfs%ks%orbs%npsidim_comp)))
   end if

!  Set orbs%eval=-0.5.
!  This will be done in LDiagHam
!  For the moment we skip this, since hpsi is not yet calculated
!  and it an input argument in LDiagHam.
   wfs%ks%orbs%eval(:)=-0.5d0

!  Copy occupations from BigDFT objects to ABINIT
   call wvl_occ_abi2big(dtset%mband,dtset%nkpt,dtset%nsppol,occ,2,wfs)


 end if

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) option,dtset%nstep,hdr0%ecut,hdr%ecut,mpi_enreg%me,wff%me,&
& wfs%ks,wvl%h(1),occ(1),rprimd(1,1),xred(1,1)
#endif

!for testing
!!Plot the density
!call sumrho(wvl_den%denspot%dpbox,wfs%ks%orbs,&
!& wfs%GPU,& wvl%atoms%sym,&
!& wvl_den%denspot%rhod,wfs%psi,wvl_den%denspot%rho_psi)
!call communicate_density(wvl_den%denspot%dpbox,wfs%ks%orbs%nspin,&
!& wvl_den%denspot%rhod,wvl_den%denspot%dpcom%nscatterarr,&
!& wvl_den%denspot%rho_psi,wvl_den%denspot%rhov)
!call plot_density('electronic_density',&
!& me,nproc,wfs%Lzd%Glr%d%n1,wfs%Lzd%Glr%d%n2,wfs%Lzd%Glr%d%n3,&
!& wfs%Lzd%Glr%d%n1i,wfs%Lzd%Glr%d%n2i,wfs%Lzd%Glr%d%n3i,&
!& wvl_den%denspot%dpcom%nscatterarr(me,2),  &
!& wfs%orbs%nspin,&
!& wvl_den%denspot%hgrids(1),wvl_den%denspot%hgrids(2),wvl_den%denspot%hgrids(3),&
!& wvl%atoms,xcart,wvl_den%denspot%dpcom%ngatherarr,&
!& wvl_den%denspot%rhov(1+wvl_den%denspot%dpcom%nscatterarr(me,4)*wfs%Lzd%Glr%d%n1i*wfs%Lzd%Glr%d%n2i))
!ABI_DEALLOCATE(xcart)
!end of debug


end subroutine wvl_wfsinp_disk
!!***

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


subroutine wvl_wfsinp_reformat(dtset, mpi_enreg, psps, rprimd, wvl, xred, xred_old)

#if defined HAVE_BIGDFT
 use BigDFT_API, only : copy_old_wavefunctions, reformatmywaves, first_orthon, &
& deallocate_wfd, wavefunctions_descriptors, deallocate_lr, &
& local_potential_dimensions, copy_coulomb_operator, &
& deallocate_coulomb_operator, nullify_gaussian_basis
 use dynamic_memory
#endif

!Arguments ------------------------------------
  type(dataset_type), intent(inout)      :: dtset
  type(MPI_type), intent(inout)          :: mpi_enreg
  type(pseudopotential_type), intent(in) :: psps
  type(wvl_data), intent(inout)          :: wvl
  real(dp), intent(inout)                :: rprimd(3,3)
  real(dp), intent(inout)                :: xred_old(3, dtset%natom)
  real(dp), intent(inout)                :: xred(3, dtset%natom)

!Local variables-------------------------------
#if defined HAVE_BIGDFT
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

#if defined HAVE_BIGDFT

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
& nproc, dtset%nspinor, dtset%nsppol, dtset%nwfshist, dtset%occ_orig(:,1), psps, rprimd, &
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
#if defined HAVE_BIGDFT
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

!!****f* ABINIT/wvl_wfsinp_scratch
!! NAME
!! wvl_wfsinp_scratch
!!
!! FUNCTION
!! This method allocates and initialises wavefunctions with values from input guess.
!! See wvl_wfsinp_disk() or wvl_wfsinp_reformat() from other initialisation
!! routines.
!!
!! When initialised from scratch or from disk, wvl%wfs%[h]psi comes unallocated
!! and will be allocated inside this routine.
!! When initialised from memory (reformating), wvl%wfs%[h]psi will be reallocated.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=input variables.
!!  hdr0 <type(hdr_type)>=the header of wf, den and pot files (read from restart)
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  ireadwf=1 for reading from file, 0 otherwise.
!!  mpi_enreg=informations about MPI parallelization
!!  option=1 for reading a file following ABINIT format, -1 for a BigDFT format.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  wff <type(wffile_type)>= structure with informations on wf file.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  wvl <type(wvl_data)>=wavefunctions & projectors informations for wavelets.
!!
!! PARENTS
!!      inwffil
!!
!! CHILDREN
!!      input_wf_diag,mklocl_wavelets,wrtout,wvl_occ_abi2big,wvl_occopt_abi2big
!!      xred2xcart
!!
!! SOURCE

subroutine wvl_wfsinp_scratch(dtset, mpi_enreg, occ, rprimd, wvl, xred)

#if defined HAVE_BIGDFT
 use BigDFT_API, only : createIonicPotential, input_wf_diag, gaussian_basis, &
      & input_variables, calculate_rhocore, deallocate_Lzd_except_Glr, INPUT_IG_OFF,&
      & SMEARING_DIST_ERF, PSPCODE_PAW
#endif

!Arguments -------------------------------
  !scalars
  type(dataset_type), intent(in)        :: dtset
  type(MPI_type), intent(inout)         :: mpi_enreg
  type(wvl_data), intent(inout)         :: wvl
  !arrays
  real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), intent(in)                  :: rprimd(3, 3)
  real(dp), intent(in)                  :: xred(3, dtset%natom)

!Local variables-------------------------------
#if defined HAVE_BIGDFT
  character(len = 500)  :: message
  integer               :: comm,me,nproc
  integer               :: iscf_local
  integer               :: nvirt
  integer               :: ii,shift_vpsp,size_vpsp
  logical               :: onlywf=.false. ! find the wavefunctions and return
  logical               :: wvlbigdft=.false.
  real(dp), allocatable :: xcart(:,:)
  real(dp), allocatable :: rhor(:,:)
  real(dp), pointer     :: vpsp(:)
  real(dp):: elecfield(3)
  type(gaussian_basis) :: Gvirt
  type(input_variables) :: in  ! To be removed, waiting for BigDFT upgrade
#endif

! *********************************************************************

#if defined HAVE_BIGDFT

 elecfield=zero

!If usewvl: wvlbigdft indicates that the BigDFT workflow will be followed
 wvlbigdft=(dtset%usewvl==1.and.dtset%wvl_bigdft_comp==1)

 write(message, '(a,a)' ) ch10,&
& ' wvl_wfsinp_scratch: wavefunction initialisation.'
 call wrtout(std_out,message,'COLL')

 comm=mpi_enreg%comm_wvl
 me=xmpi_comm_rank(comm)
 nproc=xmpi_comm_size(comm)
!Store xcart for each atom
 ABI_ALLOCATE(xcart,(3, dtset%natom))
 call xred2xcart(dtset%natom, rprimd, xcart, xred)

!We allocate temporary arrays for rho and vpsp.
!allocate ionic potential
 if (wvl%den%denspot%dpbox%n3pi > 0) then
   size_vpsp=wvl%descr%Glr%d%n1i*wvl%descr%Glr%d%n2i*wvl%den%denspot%dpbox%n3pi
   shift_vpsp=wvl%den%denspot%dpbox%ndims(1)*wvl%den%denspot%dpbox%ndims(2) &
&   *wvl%den%denspot%dpbox%nscatterarr(me,4)
   ABI_ALLOCATE(vpsp,(size_vpsp+shift_vpsp))
 else
   ABI_ALLOCATE(vpsp,(1))
 end if

 if(.not. onlywf) then
   ABI_ALLOCATE(rhor,(dtset%nfft,dtset%nspden))
   call mklocl_wavelets(elecfield, xcart, mpi_enreg, dtset%natom, &
&   dtset%nfft, dtset%nspden, 1, rprimd, vpsp, &
&   wvl%den, wvl%descr, xcart)
   ABI_DEALLOCATE(rhor)
 end if

! IMPORTANT: onlywf=.true. does not work yet, do not change this:
! if(.not. wvlbigdft) onlywf=.true. !do not apply Hamiltonian inside input_wf_diag.
 if(dtset%usepaw==1) wvl%descr%atoms%npspcode(:)=wvl%descr%npspcode_paw_init_guess(:)
 iscf_local=dtset%iscf
 if(.not. wvlbigdft) iscf_local=0  !important to have good occ values

!This routine allocates psi, hpsi and psit.
 nvirt = 0
 in%linear = INPUT_IG_OFF
 in%nspin = dtset%nsppol
 in%exctxpar = wvl%descr%exctxpar
 in%itrpmax = dtset%nnsclo
 in%iscf = iscf_local !dtset%iscf
 in%Tel = dtset%tsmear
! if (dtset%iscf == 0) in%Tel = zero
 if (iscf_local == 0) in%Tel = zero
 in%SIC = wvl%wfs%ks%SIC
 in%orthpar = wvl%wfs%ks%orthpar

 in%occopt=dtset%occopt
 if(dtset%occopt>2 .and. dtset%occopt<7)then
   call wvl_occopt_abi2big(in%occopt,in%occopt,1)
 else
!  This will be used only for the initial wavefunctions:
   in%occopt=SMEARING_DIST_ERF
 end if

!Note: check if all required in "in" is passed.
!remove argument "in" from input_wf_diag
 call input_wf_diag(me, nproc, &
& wvl%descr%atoms, wvl%den%denspot,&
& wvl%wfs%ks%orbs, nvirt, wvl%wfs%ks%comms, &
& wvl%wfs%ks%lzd, wvl%e%energs, xcart, &
& wvl%projectors%nlpsp, &
& dtset%ixc, wvl%wfs%ks%psi, wvl%wfs%ks%hpsi, wvl%wfs%ks%psit, &
& Gvirt, in%nspin, wvl%wfs%GPU, in, onlywf, wvl%projectors%G, wvl%descr%paw)

!This provisory: wvl%descr%paw could be passed as optional
!to input_wf_diag to allocate spsi inside this routine
 if(dtset%usepaw==1) then
   wvl%descr%atoms%npspcode(:)=PSPCODE_PAW
   ABI_ALLOCATE(wvl%descr%paw%spsi,(max(wvl%wfs%ks%orbs%npsidim_orbs,wvl%wfs%ks%orbs%npsidim_comp)))
   do ii=1,size(wvl%wfs%ks%psi)
     wvl%descr%paw%spsi(ii)=wvl%wfs%ks%psi(ii)
     wvl%wfs%ks%hpsi(ii)=wvl%wfs%ks%psi(ii)
   end do
 end if

 if(wvlbigdft ) then
!  Copy occupations from BigDFT objects to ABINIT
   call wvl_occ_abi2big(dtset%mband,dtset%nkpt,dtset%nsppol,occ,2,wvl%wfs)
 else
!  Copy occupations from ABINIT to BigDFT objects
   call wvl_occ_abi2big(dtset%mband,dtset%nkpt,dtset%nsppol,occ,1,wvl%wfs)
 end if

 write(message, '(a)' ) &
& '  | wavefunctions have been calculated.'
 call wrtout(std_out,message,'COLL')

 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(vpsp)

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) dtset%nstep,mpi_enreg%me,wvl%wfs%ks,occ(1),rprimd(1,1),xred(1,1)
#endif

end subroutine wvl_wfsinp_scratch
!!***

end module m_wvl_wfsinp
!!***
