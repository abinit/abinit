!{\src2tex{textfont=tt}}
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
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine wvl_wfsinp_scratch(dtset, mpi_enreg, occ, rprimd, wvl, xred)

 use defs_basis
 use defs_abitypes
 use defs_wvltypes
 use m_wffile
 use m_profiling_abi
 use m_errors
 use m_ab7_kpoints
 use m_ab7_symmetry
 use m_xmpi

 use m_abi2big,          only : wvl_occ_abi2big,wvl_occopt_abi2big

#if defined HAVE_DFT_BIGDFT
 use BigDFT_API, only : createIonicPotential, input_wf_diag, gaussian_basis, &
      & input_variables, calculate_rhocore, deallocate_Lzd_except_Glr, INPUT_IG_OFF,&
      & SMEARING_DIST_ERF, PSPCODE_PAW
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_wfsinp_scratch'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
 use interfaces_67_common
!End of the abilint section

  implicit none

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
#if defined HAVE_DFT_BIGDFT
  character(len = 500)  :: message
  integer               :: comm,me,nproc
  integer               :: iscf_local
  integer               :: nvirt
  integer               :: shift_vpsp,size_vpsp
  logical               :: onlywf=.false. ! find the wavefuncitons and return
  logical               :: wvlbigdft=.false.
  real(dp), allocatable :: xcart(:,:)
  real(dp), allocatable :: rhor(:,:)
  real(dp), pointer     :: vpsp(:)
  real(dp):: elecfield(3)
  type(gaussian_basis) :: Gvirt  
  type(input_variables) :: in  ! To be removed, waiting for BigDFT upgrade
#endif

! *********************************************************************

#if defined HAVE_DFT_BIGDFT

 elecfield=zero

!If usewvl: wvlbigdft indicates that the BigDFT workflow will be followed
 if(dtset%usewvl==1 .and. dtset%wvl_bigdft_comp==1) wvlbigdft=.true.

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
   ABI_ALLOCATE(wvl%descr%paw%spsi,(max(wvl%wfs%ks%orbs%npsidim_orbs,wvl%wfs%ks%orbs%npsidim_comp)))
   wvl%descr%atoms%npspcode(:)=PSPCODE_PAW
!if( onlywf)
!     wvl%wfs%ks%orbs%eval=-0.5_dp
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
