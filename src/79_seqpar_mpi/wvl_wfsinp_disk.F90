!{\src2tex{textfont=tt}}
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
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine wvl_wfsinp_disk(dtset, hdr0, hdr, mpi_enreg, occ, option, rprimd, wff, wfs, wvl, xred)

 use defs_basis
 use defs_abitypes
 use defs_wvltypes
 use m_wffile
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use m_abi2big,          only : wvl_occ_abi2big

#if defined HAVE_BIGDFT
 use BigDFT_API, only : first_orthon,sumrho,communicate_density,plot_density
 use dynamic_memory
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_wfsinp_disk'
 use interfaces_14_hidewrite
 use interfaces_62_wvl_wfs
!End of the abilint section

  implicit none

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
 if(dtset%usewvl==1 .and. dtset%wvl_bigdft_comp==1) wvlbigdft=.true.

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
