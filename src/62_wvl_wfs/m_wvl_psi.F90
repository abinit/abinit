!!****m* ABINIT/m_wvl_psi
!! NAME
!!  m_wvl_psi
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (DC, MT)
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

module m_wvl_psi

  use defs_basis
  use defs_wvltypes
  use m_errors
  use m_xmpi
  use m_abicore
  use m_dtset

  use defs_datatypes, only : pseudopotential_type
  use defs_abitypes,  only : MPI_type
  use m_energies, only : energies_type
  use m_pawcprj,  only : pawcprj_type, pawcprj_alloc
  use m_abi2big,  only : wvl_vxc_abi2big, wvl_vtrial_abi2big

 implicit none

 private
!!***

 public :: wvl_hpsitopsi
 public :: wvl_psitohpsi
 public :: wvl_nl_gradient
 public :: wvl_tail_corrections
!!***

contains
!!***

!!****f* ABINIT/wvl_hpsitopsi
!! NAME
!! wvl_hpsitopsi
!!
!! FUNCTION
!! Heart of the wavelet resolution, compute new wavefunctions mixed witf previous
!! by computing the gradient of the wavefunctions knowing the external potential.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=input variables.
!!  istep=id of the current iteration (first is 1).
!!  mpi_enreg=information about MPI parallelization
!!  proj <type(wvl_projector_type)>=projectors information for wavelets.
!!  vtrial(dtset%nfft)=external potential.
!!  xcart(3,natom)=cartesian atomic coordinates
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  energies <type(energies_type)>=storage for energies computed here :
!!   | e_kinetic(OUT)=kinetic energy part of total energy
!!   | e_localpsp(OUT)=local pseudopotential part of total energy
!!   | e_nlpsp_vfock(OUT)=nonlocal psp + potential Fock ACE part of total energy
!!  residm=max value for gradient in the minimisation process.
!!  rhor(dtset%nfft)=electron density in r space
!!  wfs <type(wvl_projector_type)>=wavefunctions information for wavelets.
!!
!! PARENTS
!!      m_vtorho
!!
!! CHILDREN
!!      calculatetailcorrection,dcopy,wrtout,xmpi_allgatherv
!!
!! SOURCE


subroutine wvl_hpsitopsi(cprj,dtset,energies,istep,mcprj,mpi_enreg,residm,wvl,xcart)

#if defined HAVE_BIGDFT
  use BigDFT_API, only : hpsitopsi, calculate_energy_and_gradient
#endif
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
  real(dp), save        :: etotal_local
  integer, save         :: ids
  real(dp)              :: gnrm_zero
  integer               :: comm,me,nproc
  integer :: nlmn(dtset%natom)
#endif


! *********************************************************************

 DBG_ENTER("COLL")

#if defined HAVE_BIGDFT

 if(wvl%wfs%ks%orthpar%methOrtho .ne. 0) then
   write(message,'(2a)') ch10,&
&   'wvl_hpsitopsi: the only orthogonalization method supported for PAW+WVL is Cholesky'
   MSG_ERROR(message)
 end if

 write(message, '(a,a)' ) ch10,&
& ' wvl_hpsitopsi: compute the new wavefunction from the trial potential.'
 call wrtout(std_out,message,'COLL')

 comm=mpi_enreg%comm_wvl
 me=xmpi_comm_rank(comm)
 nproc=xmpi_comm_size(comm)

!Initialisation of mixing parameter
 if (istep == 1) then
   etotal_local = real(1.d100, dp)
   ids          = dtset%nwfshist
 end if

!WARNING! e_hartree is taken from the previous iteration as e_xc
!Update physical values

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
&   wvl%descr%paw,xcart,energies%e_nlpsp_vfock,wvl%projectors%G)
 else
   call hpsitopsi(me, nproc, istep, ids, wvl%wfs%ks,&
&   wvl%descr%atoms,wvl%projectors%nlpsp)
 end if

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

!!****f* ABINIT/wvl_psitohpsi
!! NAME
!! wvl_psitohpsi
!!
!! FUNCTION
!! Compute new trial potential and calculate the hamiltionian application into hpsi.
!!
!! INPUTS
!!  mpi_enreg=information about MPI parallelization
!!
!! OUTPUT
!!  vxc(nfft,nspden)=exchange-correlation potential (hartree)
!!  vtrial(nfft,nspden)=new potential
!!
!! NOTES
!!
!! PARENTS
!!      m_afterscfloop,m_rhotov,m_setvtr,m_vtorho
!!
!! CHILDREN
!!      calculatetailcorrection,dcopy,wrtout,xmpi_allgatherv
!!
!! SOURCE


subroutine wvl_psitohpsi(alphamix,eexctX, eexcu, ehart, ekin_sum, epot_sum, eproj_sum, eSIC_DC, &
     & itrp, iter, iscf, me, natom, nfft, nproc, nspden, rpnrm, scf, &
     & vexcu, wvl, wvlbigdft, xcart, xcstr,vtrial,vxc)

#if defined HAVE_BIGDFT
 use BigDFT_API, only: psitohpsi, KS_POTENTIAL, total_energies
#endif
 implicit none

!Arguments-------------------------------
!scalars
 integer, intent(in) :: me, nproc, itrp, iter, iscf, natom, nfft, nspden
 real(dp), intent(in) :: alphamix
 real(dp), intent(out) :: rpnrm
 logical, intent(in) :: scf
 logical, intent(in) :: wvlbigdft
 type(wvl_data), intent(inout) :: wvl
 real(dp), intent(inout) :: eexctX,eSIC_DC,ehart,eexcu,vexcu, ekin_sum, epot_sum, eproj_sum
 real(dp), dimension(6), intent(out) :: xcstr
 real(dp), intent(inout) :: xcart(3, natom)
!arrays
 real(dp),intent(out), optional :: vxc(nfft,nspden)
 real(dp),intent(out), optional :: vtrial(nfft,nspden)

!Local variables-------------------------------
!scalars
#if defined HAVE_BIGDFT
 character(len=500) :: message
 integer :: linflag = 0
 character(len=3), parameter :: unblock_comms = "OFF"
#endif

! *************************************************************************

 DBG_ENTER("COLL")

#if defined HAVE_BIGDFT

 if(wvl%descr%atoms%npspcode(1)==7) then
   call psitohpsi(me,nproc,wvl%descr%atoms,scf,wvl%den%denspot, &
&   itrp, iter, iscf, alphamix,&
&   wvl%projectors%nlpsp,xcart,linflag,unblock_comms, &
&   wvl%wfs%GPU,wvl%wfs%ks,wvl%e%energs,rpnrm,xcstr,&
&   wvl%projectors%G,wvl%descr%paw)
 else
   call psitohpsi(me,nproc,wvl%descr%atoms,scf,wvl%den%denspot, &
&   itrp, iter, iscf, alphamix,&
&   wvl%projectors%nlpsp,xcart,linflag,unblock_comms, &
&   wvl%wfs%GPU,wvl%wfs%ks,wvl%e%energs,rpnrm,xcstr)
 end if

 if(scf) then
   ehart     = wvl%e%energs%eh
   eexcu     = wvl%e%energs%exc
   vexcu     = wvl%e%energs%evxc
 end if
 eexctX    = wvl%e%energs%eexctX
 eSIC_DC   = wvl%e%energs%evsic
 ekin_sum  = wvl%e%energs%ekin
 eproj_sum = wvl%e%energs%eproj
 epot_sum  = wvl%e%energs%epot

!Correct local potential, since in BigDFT
!this variable contains more terms
!Do the following only if sumpion==.true. in psolver_rhohxc.
!For the moment it is set to false.

 epot_sum=epot_sum-real(2,dp)*wvl%e%energs%eh
 epot_sum=epot_sum-wvl%e%energs%evxc

 if(wvlbigdft) then
   call total_energies(wvl%e%energs, iter, me)
 end if

!Note: if evxc is not rested here,
!we have to rest this from etotal in prtene, afterscfcv and etotfor.
!check ABINIT-6.15.1.

 if(scf) then
   if (present(vxc)) then
     write(message, '(a,a,a,a)' ) ch10, ' wvl_psitohpsi : but why are you copying vxc :..o('
     call wrtout(std_out,message,'COLL')
     call wvl_vxc_abi2big(2,vxc,wvl%den)
   end if
   if (wvl%den%denspot%rhov_is == KS_POTENTIAL .and. present(vtrial)) then
     write(message, '(a,a,a,a)' ) ch10, ' wvl_psitohpsi : but why are you copying vtrial :..o('
     call wrtout(std_out,message,'COLL')
     call wvl_vtrial_abi2big(2,vtrial,wvl%den)
   end if
 end if

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) me,nproc,itrp,iter,iscf,natom,nfft,nspden,alphamix,rpnrm,scf,&
& wvlbigdft,wvl%wfs%ks,eexctX,eSIC_DC,ehart,eexcu,vexcu,ekin_sum,&
& epot_sum,eproj_sum,xcstr(1),xcart(1,1),vxc(1,1),vtrial(1,1)
#endif

 DBG_EXIT("COLL")

end subroutine wvl_psitohpsi
!!***

!!****f* ABINIT/wvl_nl_gradient
!! NAME
!! wvl_nl_gradient
!!
!! FUNCTION
!! Compute the non local part of the wavefunction gradient.
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
!!      m_forstr,m_vtorho
!!
!! CHILDREN
!!      calculatetailcorrection,dcopy,wrtout,xmpi_allgatherv
!!
!! SOURCE

subroutine wvl_nl_gradient(grnl, mpi_enreg, natom, rprimd, wvl, xcart)

#if defined HAVE_BIGDFT
 use BigDFT_API, only: nonlocal_forces
#endif
 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom
 type(MPI_type),intent(in) :: mpi_enreg
 type(wvl_data),intent(inout) :: wvl
!arrays
 real(dp),intent(in) :: xcart(3,natom),rprimd(3,3)
 real(dp),intent(inout) :: grnl(3,natom)

!Local variables-------------------------------
#if defined HAVE_BIGDFT
!scalars
 integer :: ia,ierr,igeo,me,nproc,spaceComm
 character(len=500) :: message
!arrays
 real(dp),allocatable :: gxyz(:,:)
 real(dp)::strtens(6,4)
#endif

! *************************************************************************

#if defined HAVE_BIGDFT

!Compute forces
 write(message, '(a,a)' ) ' wvl_nl_gradient(): compute non-local part to gradient.'
 call wrtout(std_out,message,'COLL')

!Nullify output arrays.
 grnl(:, :) = zero
 strtens(:,:)=zero

 ABI_ALLOCATE(gxyz,(3, natom))
 gxyz(:,:) = zero

!Add the nonlocal part of the forces to grtn (BigDFT routine)
 spaceComm=mpi_enreg%comm_wvl
 me=xmpi_comm_rank(spaceComm)
 nproc=xmpi_comm_size(spaceComm)
 call nonlocal_forces(wvl%descr%Glr, &
& wvl%descr%h(1), wvl%descr%h(2), wvl%descr%h(3), wvl%descr%atoms, &
& xcart, wvl%wfs%ks%orbs, wvl%projectors%nlpsp, wvl%wfs%ks%Lzd%Glr%wfd, &
& wvl%wfs%ks%psi, gxyz, .true.,strtens(1,2), &
& proj_G=wvl%projectors%G,paw=wvl%descr%paw)

 if (nproc > 1) then
   call xmpi_sum(gxyz, spaceComm, ierr)
 end if

!Forces should be in reduced coordinates.
 do ia = 1, natom, 1
   do igeo = 1, 3, 1
     grnl(igeo, ia) = - rprimd(1, igeo) * gxyz(1, ia) - &
&     rprimd(2, igeo) * gxyz(2, ia) - &
&     rprimd(3, igeo) * gxyz(3, ia)
   end do
 end do
 ABI_DEALLOCATE(gxyz)

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) natom,mpi_enreg%nproc,wvl%wfs%ks,xcart(1,1),rprimd(1,1),grnl(1,1)
#endif

end subroutine wvl_nl_gradient
!!***

!!****f* ABINIT/wvl_tail_corrections
!! NAME
!! wvl_tail_corrections
!!
!! FUNCTION
!! Perform a minimization on the wavefunctions (especially the treatment
!! of the kinetic operator) with exponentialy decreasing functions on
!! boundaries.
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
!!      m_afterscfloop
!!
!! CHILDREN
!!      calculatetailcorrection,dcopy,wrtout,xmpi_allgatherv
!!
!! SOURCE

subroutine wvl_tail_corrections(dtset, energies, etotal, mpi_enreg, psps, wvl, xcart)

#if defined HAVE_BIGDFT
  use BigDFT_API, only: CalculateTailCorrection
#endif
 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(out) :: etotal
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(energies_type),intent(inout) :: energies
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_data),intent(inout) :: wvl
!arrays
 real(dp),intent(in) :: xcart(3,dtset%natom)

!Local variables-------------------------------
#if defined HAVE_BIGDFT
!scalars
 integer :: ierr,me,nbuf,nproc,nsize,spaceComm
 real(dp) :: ekin_sum,epot_sum,eproj_sum
 logical :: parallel
 character(len=500) :: message
!arrays
 integer :: ntails(3)
 real(dp) :: atails(3)
#endif

! *************************************************************************

#if defined HAVE_BIGDFT

 spaceComm=mpi_enreg%comm_wvl
 me=xmpi_comm_rank(spaceComm)
 nproc=xmpi_comm_size(spaceComm)
 parallel = (nproc > 1)

!Write a message with the total energy before tail corrections.
 etotal = energies%e_kinetic + energies%e_hartree + energies%e_xc + &
& energies%e_localpsp + energies%e_corepsp + energies%e_fock+&
& energies%e_entropy + energies%e_elecfield + energies%e_magfield+&
& energies%e_ewald + energies%e_chempot + energies%e_vdw_dftd
 if (dtset%usepaw==0) etotal = etotal + energies%e_nlpsp_vfock
 if (dtset%usepaw/=0) etotal = etotal + energies%e_paw
 write(message,'(a,2x,e19.12)') ' Total energy before tail correction', etotal
 call wrtout(std_out, message, 'COLL')

!Calculate kinetic energy correction due to boundary conditions
 nbuf = nint(dtset%tl_radius / dtset%wvl_hgrid)
 ntails = (/ wvl%descr%Glr%d%n1, wvl%descr%Glr%d%n2, wvl%descr%Glr%d%n3 /) + 2 * nbuf
 atails = real(ntails, dp) * dtset%wvl_hgrid
 write(message,'(a,a,i6,a,A,A,3F12.6,A,A,3I12,A)') ch10,&
& ' Tail requires ',nbuf,' additional grid points around cell.', ch10, &
& '  | new acell:', atails, ch10, &
& '  | new box size for wavelets:', ntails, ch10
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')


!Calculate energy correction due to finite size effects
!---reformat potential
 nsize = wvl%descr%Glr%d%n1i * wvl%descr%Glr%d%n2i
 ABI_ALLOCATE(wvl%den%denspot%pot_work, (nsize * wvl%descr%Glr%d%n3i * dtset%nsppol))

 if (parallel) then
   call xmpi_allgatherv(wvl%den%denspot%rhov, &
&   nsize * wvl%den%denspot%dpbox%nscatterarr(me, 2), &
&   wvl%den%denspot%pot_work, nsize * wvl%den%denspot%dpbox%ngatherarr(:,2), &
&   nsize * wvl%den%denspot%dpbox%ngatherarr(:,3),spaceComm,ierr)
 else
   call dcopy(wvl%descr%Glr%d%n1i * wvl%descr%Glr%d%n2i * &
&   wvl%descr%Glr%d%n3i * dtset%nsppol,wvl%den%denspot%rhov,1,wvl%den%denspot%pot_work,1)
 end if

 if(dtset%usepaw==1) then
   call CalculateTailCorrection(me, nproc, wvl%descr%atoms, dtset%tl_radius, &
&   wvl%wfs%ks%orbs, wvl%wfs%ks%lzd%Glr, wvl%projectors%nlpsp, dtset%tl_nprccg, &
&   wvl%den%denspot%pot_work, dtset%wvl_hgrid, xcart, psps%gth_params%radii_cf, &
&   dtset%wvl_crmult, dtset%wvl_frmult, dtset%nsppol, &
&   wvl%wfs%ks%psi, .false., ekin_sum, epot_sum, eproj_sum,&
&   wvl%projectors%G,wvl%descr%paw)
 else
   call CalculateTailCorrection(me, nproc, wvl%descr%atoms, dtset%tl_radius, &
&   wvl%wfs%ks%orbs, wvl%wfs%ks%lzd%Glr, wvl%projectors%nlpsp, dtset%tl_nprccg, &
&   wvl%den%denspot%pot_work, dtset%wvl_hgrid, xcart, psps%gth_params%radii_cf, &
&   dtset%wvl_crmult, dtset%wvl_frmult, dtset%nsppol, &
&   wvl%wfs%ks%psi, .false., ekin_sum, epot_sum, eproj_sum)
 end if

 ABI_DEALLOCATE(wvl%den%denspot%pot_work)

 energies%e_kinetic = ekin_sum
 energies%e_localpsp = epot_sum - two * energies%e_hartree
 energies%e_nlpsp_vfock = eproj_sum
 energies%e_corepsp = zero
 energies%e_chempot = zero
#if defined HAVE_BIGDFT
 energies%e_localpsp = energies%e_localpsp - wvl%e%energs%evxc
#endif

 write(message,'(a,3(1x,e18.11))') ' ekin_sum,epot_sum,eproj_sum',  &
 ekin_sum,epot_sum,eproj_sum
 call wrtout(std_out, message, 'COLL')
 write(message,'(a,2(1x,e18.11))') ' ehart,eexcu', &
& energies%e_hartree,energies%e_xc
 call wrtout(std_out, message, 'COLL')

 etotal = energies%e_kinetic + energies%e_hartree + energies%e_xc + &
& energies%e_localpsp + energies%e_corepsp + energies%e_fock+&
& energies%e_entropy + energies%e_elecfield + energies%e_magfield+&
& energies%e_ewald + energies%e_vdw_dftd
 if (dtset%usepaw==0) etotal = etotal + energies%e_nlpsp_vfock
 if (dtset%usepaw/=0) etotal = etotal + energies%e_paw

 write(message,'(a,2x,e19.12)') ' Total energy with tail correction', etotal
 call wrtout(std_out, message, 'COLL')

!--- End if of tail calculation

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) etotal,mpi_enreg%nproc,dtset%nstep,energies%e_ewald,psps%npsp,&
& wvl%wfs%ks,xcart(1,1)
#endif

end subroutine wvl_tail_corrections
!!***

end module m_wvl_psi
!!***
