!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_hpsitopsi
!! NAME
!! wvl_hpsitopsi
!!
!! FUNCTION
!! Heart of the wavelet resolution, compute new wavefunctions mixed witf previous
!! by computing the gradient of the wavefunctions knowing the external potential.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=input variables.
!!  istep=id of the current iteration (first is 1).
!!  mpi_enreg=informations about MPI parallelization
!!  proj <type(wvl_projector_type)>=projectors informations for wavelets.
!!  vtrial(dtset%nfft)=external potential.
!!  xcart(3,natom)=cartesian atomic coordinates
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  energies <type(energies_type)>=storage for energies computed here :
!!   | e_kinetic(OUT)=kinetic energy part of total energy
!!   | e_localpsp(OUT)=local pseudopotential part of total energy
!!   | e_nonlocalpsp(OUT)=nonlocal pseudopotential part of total energy
!!  residm=max value for gradient in the minimisation process.
!!  rhor(dtset%nfft)=electron density in r space
!!  wfs <type(wvl_projector_type)>=wavefunctions informations for wavelets.
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!      calculate_energy_and_gradient,hpsitopsi,pawcprj_alloc,wrtout
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine wvl_hpsitopsi(cprj,dtset, energies, istep, mcprj,mpi_enreg, &
&  residm, wvl,xcart)

 use m_profiling_abi

  use defs_basis
  use defs_abitypes
  use defs_wvltypes
  use m_errors
  use m_xmpi

  use m_energies, only : energies_type
  use m_pawcprj, only : pawcprj_type, pawcprj_alloc
#if defined HAVE_BIGDFT
  use BigDFT_API, only : hpsitopsi, calculate_energy_and_gradient
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_hpsitopsi'
 use interfaces_14_hidewrite
!End of the abilint section

  implicit none

!Arguments -------------------------------
  type(dataset_type), intent(in)         :: dtset
  type(energies_type), intent(inout)     :: energies
  integer, intent(in)                    :: istep,mcprj
  type(MPI_type), intent(in)             :: mpi_enreg
  real(dp), intent(inout)                :: residm
  type(wvl_data), intent(inout)          :: wvl
  !arrays
  real(dp), intent(in) :: xcart(3, dtset%natom)
  type(pawcprj_type),dimension(dtset%natom,mcprj),intent(inout)::cprj

!Local variables-------------------------------
#if defined HAVE_BIGDFT
  integer               :: iatom,icprj
  character(len = 500)  :: message
  logical               :: wvlbigdft=.false.
  real(dp), save        :: etotal_local
  integer, save         :: ids
  real(dp)              :: gnrm_zero
  integer               :: comm,me,nproc
  real(dp)              :: eproj_sum
  integer :: nlmn(dtset%natom)
 !debug
 ! integer::ierr,ii,jj,ifile
 !end debug
#endif


! *********************************************************************

 DBG_ENTER("COLL")

#if defined HAVE_BIGDFT

!If usewvl: wvlbigdft indicates that the BigDFT workflow will be followed
 if(dtset%wvl_bigdft_comp==1) wvlbigdft=.true.

 if(wvl%wfs%ks%orthpar%methOrtho .ne. 0) then
   write(message,'(2a)') ch10,&
&   'wvl_hpsitopsi: the only orthogonalization method supported for PAW+WVL is Cholesky'
   MSG_ERROR(message)
 end if


!Variables now grouped in BigDFT:
!energies%e_exactX        =eexctX
!energies%e_xc            =eexcu
!energies%e_hartree       =ehart
!energies%e_kinetic       =ekin_sum
!energies%e_localpsp      =epot_sum
!energies%e_sicdc         =eSIC_DC
!energies%e_xcdc,          =vexcu
!energies%e_localpsp = energies%e_localpsp - real(2, dp) * energies%e_hartree
!if(NC) : energies%e_nonlocalpsp   =eproj_sum 
!if(PAW): energies%e_paw           =eproj_sum 
 if(dtset%usepaw==1) eproj_sum=energies%e_paw
 if(dtset%usepaw==0) eproj_sum=energies%e_nonlocalpsp

 write(message, '(a,a)' ) ch10,&
& ' wvl_hpsitopsi: compute the new wavefunction from the trial potential.'
 call wrtout(std_out,message,'COLL')

 comm=mpi_enreg%comm_wvl
 me=xmpi_comm_rank(comm)
 nproc=xmpi_comm_size(comm)

!write(*,*)'wvl_hpsitopsi l138,erase me'
!ifile=me+40
!jj=size(wvl%wfs%ks%hpsi)
!ii=size(wvl%wfs%ks%psi)
!write(ifile,*)"# ",jj,ii
!do ii=1,jj
!write(ifile,'(2f20.8)')wvl%wfs%ks%hpsi(ii),wvl%wfs%ks%psi(ii)
!end do
!call mpi_barrier(comm,ierr)
!
!write(*,*)'erase me wvl_hpsitopsi, stop'
!stop


!Initialisation of mixing parameter
 if (istep == 1) then
   etotal_local = real(1.d100, dp)
   ids          = dtset%nwfshist
 end if

!WARNING! e_hartree is taken from the previous iteration as e_xc
!Update physical values
 if(wvlbigdft) energies%e_corepsp = zero

!Precondition, minimise (DIIS or steepest descent) and ortho.
!Compute also the norm of the gradient.
 if(dtset%usepaw==1) then
   call calculate_energy_and_gradient(istep, me, nproc, wvl%wfs%GPU, dtset%wvl_nprccg, &
&   dtset%iscf, wvl%e%energs, wvl%wfs%ks, residm, gnrm_zero,wvl%descr%paw)
 else
   call calculate_energy_and_gradient(istep, me, nproc, wvl%wfs%GPU, dtset%wvl_nprccg, &
&   dtset%iscf, wvl%e%energs, wvl%wfs%ks, residm, gnrm_zero)
 end if
 etotal_local = wvl%wfs%ks%diis%energy

 if(dtset%usepaw==1) then
   call hpsitopsi(me, nproc, istep, ids, wvl%wfs%ks,&
&   wvl%descr%atoms,wvl%projectors%nlpsp,&
&   wvl%descr%paw,xcart,eproj_sum,wvl%projectors%G)
 else
   call hpsitopsi(me, nproc, istep, ids, wvl%wfs%ks,&
&   wvl%descr%atoms,wvl%projectors%nlpsp)
 end if

 if(dtset%usepaw==1) energies%e_paw=eproj_sum
 if(dtset%usepaw==0) energies%e_nonlocalpsp=eproj_sum

 if(dtset%usepaw==1) then
!  PENDING : cprj should not be copied'

   ! Cannot use pawcprj_copy because cprj and paw%cprj are not the same objects
   ! Get nlmn from bigdft cprj, and allocate our copy
   do iatom=1,dtset%natom
     nlmn(iatom) = wvl%descr%paw%cprj(iatom,1)%nlmn
   end do
   call pawcprj_alloc(cprj,mcprj,nlmn)

   do iatom=1,dtset%natom
     do icprj=1,mcprj
       cprj(iatom,icprj)%cp(:,:)= wvl%descr%paw%cprj(iatom,icprj)%cp(:,:)
     end do
   end do


 end if

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) dtset%nstep,energies%e_ewald,istep,mcprj,mpi_enreg%nproc,residm,&
& wvl%wfs%ks,xcart(1,1),cprj(1,1)%nlmn
#endif

 DBG_EXIT("COLL")

end subroutine wvl_hpsitopsi
!!***
