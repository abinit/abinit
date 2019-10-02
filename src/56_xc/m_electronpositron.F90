!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_electronpositron
!! NAME
!!  m_electronpositron
!!
!! FUNCTION
!!  This module provides the definition of the electronpositron_type used
!!  used to store data for the electron-positron two-component DFT
!!  as methods to operate on it.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2019 ABINIT group (MT, GJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_electronpositron

 use defs_basis
 use m_abicore
 use m_errors
 use m_energies
 use m_xmpi
 use m_cgtools
 use m_dtset

 use defs_abitypes, only : MPI_type
 use m_pawtab,   only : pawtab_type
 use m_paw_an,   only : paw_an_type
 use m_pawrhoij, only : pawrhoij_type, pawrhoij_alloc, pawrhoij_free, pawrhoij_copy
 use m_pawcprj,  only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_copy
 use m_mpinfo,   only : proc_distrb_cycle
 use m_xcpositron, only : xcpositron
 use m_drivexc,    only : mkdenpos
 use m_xctk,       only : xcden
 use m_fft,        only : fourdp

 implicit none

 private

! public constants
 integer,public,parameter :: EP_NOTHING  =-1
 integer,public,parameter :: EP_ELECTRON = 0
 integer,public,parameter :: EP_POSITRON = 1
!!***

!!****t* m_electronpositron/electronpositron_type
!! NAME
!!
!! FUNCTION
!!
!! NOTES
!!
!! SOURCE

 type, public :: electronpositron_type

! Integer scalars
  integer :: calctype        ! type of electron-positron calculation:
                             !   0: no calculation
                             !   1: positron in the electrons potential
                             !   2: electrons in the positron potential
  integer :: particle        ! current particle stored in electronpositron%xxx_ep arrays
                             !                 -1: no particle, 0: electron, 1: positron
  integer :: dimcg           ! Dimension of cg array dimcg=dtset%mpw*dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol
  integer :: dimcprj         ! Dimension of cprj array dimcprj=dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol*usecprj
  integer :: dimeigen        ! Dimension of eigen array dimeigen=dtset%mband*dtset%nkpt*dtset%nsppol
  integer :: dimocc          ! Dimension of occ array dimocc=dtset%mband*dtset%nkpt*dtset%nsppol
  integer :: has_pawrhoij_ep ! flag for pawrhoij_ep (0: not allocated, 1: allocated, 2: computed)
  integer :: has_pos_ham     ! flag: 1 if current Hamiltonian in memory (vtrial, vpsp, vhartr, vxc, paw_ij%dij)
!                                    is the positronic hamiltonian, 0 is it is the electronic one
  integer :: ixcpositron     ! XC type for electron-positron correlation
  integer :: istep           ! Current index of TC-DFT SCF step
  integer :: istep_scf       ! Current index of DFT SCF step  in current electron/positron minimization
  integer :: lmmax           ! Max. number of (l,m) moments over all types of atom
  integer :: natom           ! Number of atoms
  integer :: nfft            ! Number of points in FFT grid
  integer :: nspden          ! Number of spin density components
  integer :: nstep           ! Max. number of steps for the TC-DFT SCF cycle

! Logical scalars
  logical :: posdensity0_limit ! True if we are in the zero positron density limit
  logical :: scf_converged     ! True if the SCF cycle is converged for a positronic/electronic GS calculation

! Real(dp) scalars
  real(dp) :: e_hartree      !  Hartree electron-positron interaction energy
  real(dp) :: e_xc           !  XC electron-positron interaction energy
  real(dp) :: e_xcdc         !  Double-counting XC electron-positron interaction energy
  real(dp) :: e_paw          !  PAW electron-positron interaction energy
  real(dp) :: e_pawdc        !  Double-counting PAW electron-positron interaction energy
  real(dp) :: e0             !  Energy only due to particle(s) currently evolving
                                  !   calctype=1, energy due to positron  only
                                  !   calctype=2, energy due to electrons only
  real(dp) :: etotal_prev    !  Total energy of the previous GS calculation
  real(dp) :: lambda         ! Electron-positron annihilation rate
  real(dp) :: lifetime       ! Positron lifetime
  real(dp) :: maxfor_prev    ! Max. force of the previous GS calculation
  real(dp) :: posocc         ! Occupation number for the positron
  real(dp) :: postoldfe      ! Tolerance on total energy for the TC-DFT SCF cycle
  real(dp) :: postoldff      ! Tolerance on max. force for the TC-DFT SCF cycle

! Other scalars
  type(energies_type) :: energies_ep  !  Energies of the previous electronic/positronic SCF step

! Logical pointers
  logical, allocatable :: lmselect_ep(:,:)
!  lmselect_ep(lmmax,my_natom)
!  flags selecting the non-zero LM-moments of on-site densities

! Real(dp) pointers
  real(dp), allocatable :: cg_ep(:,:)
!  cg_ep(2,dimcg)
!  if typecalc=1: electronic wavefunctions
!  if typecalc=2: positronic wavefunctions

  real(dp), allocatable :: eigen_ep(:)
!  eigen(dimeigen)
!  if typecalc=1: electronic eigen energies
!  if typecalc=2: positronic eigen energies

  real(dp), allocatable :: fred_ep(:,:)
!  fred_ep(3,natom)
!  if typecalc=1: forces only due to electrons
!  if typecalc=2: forces only due to positron

  real(dp), allocatable :: nhat_ep(:,:)
!  nhat_ep(nfft,nspden)
!  if typecalc=1: electronic compensation charge density in real space
!  if typecalc=2: positronic compensation charge density in real space

  real(dp), allocatable :: occ_ep(:)
!  occ(dimocc)
!  if typecalc=1: electronic occupations
!  if typecalc=2: positronic occupations

  real(dp), allocatable :: rhor_ep(:,:)
!  rhor_ep(nfft,nspden)
!  if typecalc=1: electronic density in real space
!  if typecalc=2: positronic density in real space

  real(dp), allocatable :: stress_ep(:)
!  stress_ep(6)
!  if typecalc=1: stresses only due to electrons
!  if typecalc=2: stresses only due to positron

  real(dp), allocatable :: vha_ep(:)
!  vha_ep(nfft)
!  if typecalc=1: electronic Hartree potential
!  if typecalc=2: positronic Hartree potential

! Other pointers
  type(pawrhoij_type), allocatable :: pawrhoij_ep(:)
!  pawrhoij_ep(natom)
!  Relevant only if PAW
!  if typecalc=1: electronic PAW occupation matrix associated with rhor_ep
!  if typecalc=2: positronic PAW occupation matrix associated with rhor_ep

  type(pawcprj_type), allocatable :: cprj_ep(:,:)
!  cprj_ep(natom,dimcprj)
!  Relevant only if PAW
!  if typecalc=1: electronic WF projected on nl projectors <p_i|Cnk>
!  if typecalc=2: positronic WF projected on nl projectors <p_i|Cnk>

 end type electronpositron_type

! public procedures
 public :: init_electronpositron
 public :: destroy_electronpositron
 public :: exchange_electronpositron
 public :: electronpositron_calctype
 public :: rhohxcpositron

CONTAINS

!===========================================================
!!***

!!****f* m_electronpositron/init_electronpositron
!! NAME
!!  init_electronpositron
!!
!! FUNCTION
!!  Init all scalars and pointers in the structure.
!!
!! INPUTS
!!  ireadwf=if 1, read the wavefunction
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  mpi_enreg=information about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  pawrhoij(natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic occupancies
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!
!! SIDE EFFECTS
!!  electronpositron=<type(electronpositron_type)>=electronpositron datastructure
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      energies_copy,fourdp,pawcprj_alloc,pawcprj_copy,pawcprj_free
!!      pawrhoij_alloc,pawrhoij_copy,pawrhoij_free
!!
!! SOURCE

subroutine init_electronpositron(ireadwf,dtset,electronpositron,mpi_enreg,nfft,pawrhoij,pawtab)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ireadwf,nfft
 type(dataset_type),intent(in) :: dtset
 type(electronpositron_type),pointer :: electronpositron
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 type(pawrhoij_type), intent(in) :: pawrhoij(mpi_enreg%my_natom*dtset%usepaw)
 type(pawtab_type),intent(in)  :: pawtab(dtset%ntypat*dtset%usepaw)

!Local variables-------------------------------
!scalars
 integer :: ii,my_nspinor,ncpgr,optfor,optstr
 logical,parameter :: include_nhat_in_gamma=.false.
!arrays
 integer,allocatable :: nlmn(:)

!************************************************************************

 !@electronpositron_type

 if (dtset%positron/=0) then

  ABI_DATATYPE_ALLOCATE(electronpositron,)

  electronpositron%calctype=0
  electronpositron%particle=-1

  electronpositron%ixcpositron=dtset%ixcpositron
  electronpositron%natom=dtset%natom
  electronpositron%nfft=nfft
  electronpositron%nspden=dtset%nspden
  electronpositron%istep=0
  electronpositron%istep_scf=0

  electronpositron%posocc=dtset%posocc
  electronpositron%nstep=dtset%posnstep
  electronpositron%postoldfe=dtset%postoldfe
  electronpositron%postoldff=dtset%postoldff
  electronpositron%posdensity0_limit=(dtset%ixcpositron/=2)
  electronpositron%scf_converged=.false.
  electronpositron%has_pos_ham=0

  call energies_init(electronpositron%energies_ep)

  electronpositron%e_hartree  =zero
  electronpositron%e_xc       =zero
  electronpositron%e_xcdc     =zero
  electronpositron%e_paw      =zero
  electronpositron%e_pawdc    =zero
  electronpositron%e0         =zero
  electronpositron%etotal_prev=zero
  electronpositron%maxfor_prev=zero

  electronpositron%lambda=zero
  electronpositron%lifetime=zero

  ABI_ALLOCATE(electronpositron%rhor_ep,(nfft,dtset%nspden))
  ABI_ALLOCATE(electronpositron%vha_ep,(nfft))
  ABI_DATATYPE_ALLOCATE(electronpositron%pawrhoij_ep,(mpi_enreg%my_natom*dtset%usepaw))

  if (dtset%usepaw==1) then
   electronpositron%has_pawrhoij_ep=1
   if (mpi_enreg%my_natom>0) then
    call pawrhoij_alloc(electronpositron%pawrhoij_ep,pawrhoij(1)%cplex_rhoij,pawrhoij(1)%nspden,&
&                    pawrhoij(1)%nspinor,pawrhoij(1)%nsppol,dtset%typat,&
&                    mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom,&
&                    pawtab=pawtab,ngrhoij=pawrhoij(1)%ngrhoij,nlmnmix=pawrhoij(1)%lmnmix_sz,&
&                    qphase=pawrhoij(1)%qphase,use_rhoij_=pawrhoij(1)%use_rhoij_,&
&                    use_rhoijres=pawrhoij(1)%use_rhoijres)
   end if
   electronpositron%lmmax=0
   do ii=1,dtset%ntypat
    electronpositron%lmmax=max(electronpositron%lmmax,pawtab(ii)%lcut_size**2)
   end do
   ABI_ALLOCATE(electronpositron%lmselect_ep,(electronpositron%lmmax,mpi_enreg%my_natom))
   if (maxval(pawtab(1:dtset%ntypat)%usexcnhat)==0.or.(.not.include_nhat_in_gamma)) then
     ABI_ALLOCATE(electronpositron%nhat_ep,(nfft,dtset%nspden))
   end if
  else
   electronpositron%has_pawrhoij_ep=0
   electronpositron%lmmax=0
  end if

  if (dtset%positron<=-10.or.dtset%posdoppler>0) then
   my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
   electronpositron%dimcg=dtset%mpw*my_nspinor*dtset%mband*dtset%mkmem*dtset%nsppol
   electronpositron%dimocc=dtset%mband*dtset%nkpt*dtset%nsppol
   electronpositron%dimeigen=dtset%mband*dtset%nkpt*dtset%nsppol
   ABI_ALLOCATE(electronpositron%cg_ep,(2,electronpositron%dimcg))
   ABI_ALLOCATE(electronpositron%eigen_ep,(electronpositron%dimeigen))
   ABI_ALLOCATE(electronpositron%occ_ep,(electronpositron%dimocc))
   electronpositron%dimcprj=0
!  if (.false.) then !TEMPORARY: will be activated later
   if (dtset%usepaw==1.and.dtset%pawusecp>0.and.dtset%posdoppler>0) then
    electronpositron%dimcprj=my_nspinor*dtset%mband*dtset%mkmem*dtset%nsppol
    if (mpi_enreg%paral_kgb/=0) electronpositron%dimcprj=electronpositron%dimcprj/mpi_enreg%nproc_band
    ABI_DATATYPE_ALLOCATE(electronpositron%cprj_ep,(dtset%natom,electronpositron%dimcprj))
    ABI_ALLOCATE(nlmn,(dtset%natom))
    ncpgr=0
    do ii=1,dtset%natom;nlmn(ii)=pawtab(dtset%typat(ii))%lmn_size;end do
    call pawcprj_alloc(electronpositron%cprj_ep,ncpgr,nlmn)
    ABI_DEALLOCATE(nlmn)
   else
    ABI_DATATYPE_ALLOCATE(electronpositron%cprj_ep,(dtset%natom,electronpositron%dimcprj))
   end if
  else
   electronpositron%dimcg   =0
   electronpositron%dimcprj =0
   electronpositron%dimeigen=0
   electronpositron%dimocc  =0
  end if

  optfor=0;optstr=0
  if ((dtset%optforces>0.or.dtset%ionmov/=0.or.abs(dtset%toldff)>tiny(0._dp))) optfor=1
  if (dtset%optstress>0.and.dtset%iscf>0.and.(dtset%nstep>0.or.ireadwf==1)) optstr=1

  if (optfor>0) then
   ABI_ALLOCATE(electronpositron%fred_ep,(3,dtset%natom))
   electronpositron%fred_ep(:,:)=zero
  end if

  if (optstr>0) then
   ABI_ALLOCATE(electronpositron%stress_ep,(6))
   electronpositron%stress_ep(:)=zero
  end if

 else !dtset%positron==0
  nullify(electronpositron)
 end if

end subroutine init_electronpositron
!!***

!----------------------------------------------------------------------

!!****f* m_electronpositron/destroy_electronpositron
!! NAME
!!  destroy_electronpositron
!!
!! FUNCTION
!!  Clean and destroy electronpositron datastructure
!!
!! SIDE EFFECTS
!!  electronpositron=<type(electronpositron_type)>=electronpositron datastructure
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      energies_copy,fourdp,pawcprj_alloc,pawcprj_copy,pawcprj_free
!!      pawrhoij_alloc,pawrhoij_copy,pawrhoij_free
!!
!! SOURCE

subroutine destroy_electronpositron(electronpositron)

!Arguments ------------------------------------
!scalars
 type(electronpositron_type),pointer :: electronpositron

!************************************************************************

 !@electronpositron_type

 if (associated(electronpositron)) then

  if (allocated(electronpositron%cg_ep))        then
    ABI_DEALLOCATE(electronpositron%cg_ep)
  end if
  if (allocated(electronpositron%eigen_ep))     then
    ABI_DEALLOCATE(electronpositron%eigen_ep)
  end if
  if (allocated(electronpositron%occ_ep))       then
    ABI_DEALLOCATE(electronpositron%occ_ep)
  end if
  if (allocated(electronpositron%rhor_ep))      then
    ABI_DEALLOCATE(electronpositron%rhor_ep)
  end if
  if (allocated(electronpositron%nhat_ep))      then
    ABI_DEALLOCATE(electronpositron%nhat_ep)
  end if
  if (allocated(electronpositron%vha_ep))       then
    ABI_DEALLOCATE(electronpositron%vha_ep)
  end if
  if (allocated(electronpositron%lmselect_ep))  then
    ABI_DEALLOCATE(electronpositron%lmselect_ep)
  end if
  if (allocated(electronpositron%fred_ep))      then
    ABI_DEALLOCATE(electronpositron%fred_ep)
  end if
  if (allocated(electronpositron%stress_ep))    then
    ABI_DEALLOCATE(electronpositron%stress_ep)
  end if

  if (electronpositron%has_pawrhoij_ep/=0) then
   call pawrhoij_free(electronpositron%pawrhoij_ep)
  end if
  if (allocated(electronpositron%pawrhoij_ep))  then
    ABI_DATATYPE_DEALLOCATE(electronpositron%pawrhoij_ep)
  end if

  if (electronpositron%dimcprj/=0) then
   call pawcprj_free(electronpositron%cprj_ep)
  end if
  if (allocated(electronpositron%cprj_ep))  then
    ABI_DATATYPE_DEALLOCATE(electronpositron%cprj_ep)
  end if

  electronpositron%calctype       =0
  electronpositron%particle       =-1
  electronpositron%dimcg          =0
  electronpositron%dimcprj        =0
  electronpositron%dimeigen       =0
  electronpositron%dimocc         =0
  electronpositron%has_pawrhoij_ep=0
  electronpositron%has_pos_ham    =0
  electronpositron%istep          =0
  electronpositron%istep_scf      =0

  electronpositron%posdensity0_limit=.false.
  electronpositron%scf_converged=.false.

  ABI_DATATYPE_DEALLOCATE(electronpositron)

 end if

end subroutine destroy_electronpositron
!!***

!----------------------------------------------------------------------

!!****f* m_electronpositron/exchange_electronpositron
!! NAME
!!  exchange_electronpositron
!!
!! FUNCTION
!!  Invert electron and positron quantities between an electronpositron datastructure
!!  and current evoving variables
!!  Example: exchange electronpositron%rhor_ep and rhor
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  mpi_enreg=information about MPI parallelization
!!  my_natom=number of atoms treated by current proc
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  usecprj= 1 if cprj array is stored in memory
!!
!! SIDE EFFECTS
!!  cg(2,mcg)=wavefunctions
!!  cprj(natom,mcprj*usecprj)= wave functions projected with non-local projectors:
!!                             cprj(n,k,i)=<p_i|Cnk> where p_i is a non-local projector.
!!  electronpositron=<type(electronpositron_type)>=electronpositron datastructure
!!  energies <type(energies_type)>=all part of total energy.
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  fred(3,natom)=forces in reduced coordinates
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  occ(mband*nkpt*nsppol)=occupation number for each band at each k point
!!  paw_an(natom) <type(paw_an_type)>=paw arrays given on angular mesh
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  rhog(2,nfft)=Fourier transform of total electron/positron density
!!  rhor(nfft,nspden)=total electron/positron density (el/bohr**3)
!!  stress(6)=components of the stress tensor (hartree/bohr^3) for the
!!  vhartr(nfftf)=array for holding Hartree potential
!!
!! PARENTS
!!      afterscfloop
!!
!! CHILDREN
!!      energies_copy,fourdp,pawcprj_alloc,pawcprj_copy,pawcprj_free
!!      pawrhoij_alloc,pawrhoij_copy,pawrhoij_free
!!
!! SOURCE

subroutine exchange_electronpositron(cg,cprj,dtset,eigen,electronpositron,energies,fred,mcg,mcprj,&
&                                    mpi_enreg,my_natom,nfft,ngfft,nhat,npwarr,occ,paw_an,pawrhoij,&
&                                    rhog,rhor,stress,usecprj,vhartr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mcg,mcprj,my_natom,nfft,usecprj
 type(dataset_type),intent(in) :: dtset
 type(electronpositron_type),pointer :: electronpositron
 type(energies_type),intent(inout) :: energies
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18),npwarr(dtset%nkpt)
 real(dp),intent(inout) :: cg(2,mcg)
 real(dp),intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(inout) :: fred(3,dtset%natom),nhat(nfft,dtset%nspden)
 real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(inout) :: rhog(2,nfft),rhor(nfft,dtset%nspden)
 real(dp),intent(inout) :: stress(6),vhartr(nfft)
 type(pawcprj_type) :: cprj(dtset%natom,mcprj*usecprj)
 type(paw_an_type),intent(inout) :: paw_an(my_natom*dtset%usepaw)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*dtset%usepaw)

!Local variables-------------------------------
!scalars
 integer :: comm,iatom,ib,ibsp,icg,icgb,ifft,ii,ilm,ikpt
 integer :: ispden,isppol,ispinor,me,my_nspinor,nband_k,npw_k,sz1,sz2,sz3
 logical :: ltmp
 real(dp) :: rtmp
 type(energies_type) :: energies_tmp
!arrays
 integer,allocatable :: nlmn(:),typ(:)
 real(dp) :: ctmp(2)
 type(pawcprj_type),allocatable :: cprj_tmp(:,:)
 type(pawrhoij_type),allocatable :: pawrhoij_tmp(:)

!*********************************************************************

 if (associated(electronpositron)) then
  if (electronpositron%particle/=EP_NOTHING) then

!  Type of particle stored
   if (electronpositron%particle==EP_ELECTRON) then
     electronpositron%particle=EP_POSITRON
   else if (electronpositron%particle==EP_POSITRON) then
     electronpositron%particle=EP_ELECTRON
   end if

!  Energies
   ctmp(1)=energies%e_electronpositron
!  ctmp(2)=energies%edc_electronpositron
   call energies_copy(electronpositron%energies_ep,energies_tmp)
   call energies_copy(energies,electronpositron%energies_ep)
   call energies_copy(energies_tmp,energies)
   energies%e_electronpositron=ctmp(1)
!  energies%edc_electronpositron=ctmp(2)
   energies%e0_electronpositron=electronpositron%e0
   electronpositron%e0=electronpositron%energies_ep%e0_electronpositron

!  Density and PAW occupation matrix
   do ispden=1,dtset%nspden
     do ifft=1,nfft
       rtmp=rhor(ifft,ispden)
       rhor(ifft,ispden)=electronpositron%rhor_ep(ifft,ispden)
       electronpositron%rhor_ep(ifft,ispden)=rtmp
     end do
     if (allocated(electronpositron%nhat_ep).and.size(nhat,2)>0) then
       do ifft=1,nfft
         rtmp=nhat(ifft,ispden)
         nhat(ifft,ispden)=electronpositron%nhat_ep(ifft,ispden)
         electronpositron%nhat_ep(ifft,ispden)=rtmp
       end do
     end if
   end do
   call fourdp(1,rhog,rhor,-1,mpi_enreg,nfft,1,ngfft,0)
   if (dtset%usepaw==1.and.my_natom>0) then
    if (electronpositron%has_pawrhoij_ep==1) then
      ABI_DATATYPE_ALLOCATE(pawrhoij_tmp,(my_natom))
      ABI_ALLOCATE(typ,(my_natom))
      ABI_ALLOCATE(nlmn,(my_natom))
      do iatom=1,my_natom
        typ(iatom)=iatom
        nlmn(iatom)=pawrhoij(iatom)%lmn_size
      end do
!     Be careful: parallelism over atoms is ignored...
      call pawrhoij_alloc(pawrhoij_tmp,pawrhoij(1)%cplex_rhoij,pawrhoij(1)%nspden,&
&                      pawrhoij(1)%nspinor,pawrhoij(1)%nsppol,typ, &
&                      lmnsize=nlmn,ngrhoij=pawrhoij(1)%ngrhoij,nlmnmix=pawrhoij(1)%lmnmix_sz,&
&                      qphase=pawrhoij(1)%qphase,use_rhoij_=pawrhoij(1)%use_rhoij_,&
&                      use_rhoijres=pawrhoij(1)%use_rhoijres)
      ABI_DEALLOCATE(typ)
      ABI_DEALLOCATE(nlmn)
      call pawrhoij_copy(pawrhoij,pawrhoij_tmp)
      call pawrhoij_copy(electronpositron%pawrhoij_ep,pawrhoij)
      call pawrhoij_copy(pawrhoij_tmp,electronpositron%pawrhoij_ep)
      if (pawrhoij_tmp(1)%ngrhoij>0.and.pawrhoij(1)%ngrhoij==0) then
        do iatom=1,my_natom
          sz1=pawrhoij_tmp(iatom)%ngrhoij
          sz2=pawrhoij_tmp(iatom)%cplex_rhoij*pawrhoij_tmp(iatom)%qphase*pawrhoij_tmp(iatom)%lmn2_size
          sz3=pawrhoij_tmp(iatom)%nspden
          ABI_ALLOCATE(pawrhoij(iatom)%grhoij,(sz1,sz2,sz3))
          pawrhoij(iatom)%grhoij(:,:,:)=pawrhoij_tmp(iatom)%grhoij(:,:,:)
        end do
      end if
      if (pawrhoij_tmp(1)%use_rhoijres>0.and.pawrhoij(1)%use_rhoijres==0) then
        do iatom=1,my_natom
          sz1=pawrhoij_tmp(iatom)%cplex_rhoij*pawrhoij_tmp(iatom)%qphase*pawrhoij_tmp(iatom)%lmn2_size
          sz2=pawrhoij_tmp(iatom)%nspden
          ABI_ALLOCATE(pawrhoij(iatom)%rhoijres,(sz1,sz2))
          pawrhoij(iatom)%rhoijres(:,:)=pawrhoij_tmp(iatom)%rhoijres(:,:)
        end do
      end if
      if (pawrhoij_tmp(1)%use_rhoij_>0.and.pawrhoij(1)%use_rhoij_==0) then
        do iatom=1,my_natom
          sz1=pawrhoij_tmp(iatom)%cplex_rhoij*pawrhoij_tmp(iatom)%qphase*pawrhoij_tmp(iatom)%lmn2_size
          sz2=pawrhoij_tmp(iatom)%nspden
          ABI_ALLOCATE(pawrhoij(iatom)%rhoij_,(sz1,sz2))
          pawrhoij(iatom)%rhoij_(:,:)=pawrhoij_tmp(iatom)%rhoij_(:,:)
        end do
      end if
      if (pawrhoij_tmp(1)%lmnmix_sz>0.and.pawrhoij(1)%lmnmix_sz==0) then
        do iatom=1,my_natom
          ABI_ALLOCATE(pawrhoij(iatom)%kpawmix,(pawrhoij_tmp(iatom)%lmnmix_sz))
          pawrhoij(iatom)%kpawmix(:)=pawrhoij_tmp(iatom)%kpawmix(:)
        end do
      end if
      call pawrhoij_free(pawrhoij_tmp)
      ABI_DATATYPE_DEALLOCATE(pawrhoij_tmp)
    else
      do iatom=1,my_natom
        pawrhoij(iatom)%rhoijp=zero
      end do
    end if
   end if

!  Hartree potential
   do ifft=1,nfft
    rtmp=vhartr(ifft)
    vhartr(ifft)=electronpositron%vha_ep(ifft)
    electronpositron%vha_ep(ifft)=rtmp
   end do

!  PAW LM-moment selection flags
   if (dtset%usepaw==1.and.my_natom>0) then
    if (electronpositron%lmmax>0) then
     do iatom=1,my_natom
      do ilm=1,paw_an(iatom)%lm_size
       ltmp=electronpositron%lmselect_ep(ilm,iatom)
       electronpositron%lmselect_ep(ilm,iatom)=paw_an(iatom)%lmselect(ilm)
       paw_an(iatom)%lmselect(ilm)=ltmp
      end do
     end do
    else
     do iatom=1,my_natom
      paw_an(iatom)%lmselect(:)=.true.
     end do
    end if
   end if

!  Wave-functions
   if (electronpositron%dimcg>0) then
    do ii=1,electronpositron%dimcg
     ctmp(1:2)=electronpositron%cg_ep(1:2,ii)
     electronpositron%cg_ep(1:2,ii)=cg(1:2,ii)
     cg(1:2,ii)=ctmp(1:2)
    end do
   else
    icg=0
    my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
    comm=mpi_enreg%comm_cell
    me=xmpi_comm_rank(comm)
    do isppol=1,dtset%nsppol
     do ikpt=1,dtset%nkpt
      npw_k=npwarr(ikpt);nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
      if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) cycle
      icgb=icg;ibsp=0
      do ib=1,nband_k
       cg(:,icgb+1:icgb+my_nspinor*npw_k)=zero
       do ispinor=1,my_nspinor
        ibsp=ibsp+1;if (ibsp<my_nspinor*npw_k) cg(1,icgb+ibsp)=one
       end do
       icgb=icgb+my_nspinor*npw_k
      end do
      if (dtset%mkmem/=0) icg=icg+my_nspinor*npw_k*nband_k
     end do
    end do
   end if
   if (dtset%usepaw==1) then
    if(electronpositron%dimcprj>0) then
     ABI_ALLOCATE(nlmn,(dtset%natom))
     ABI_DATATYPE_ALLOCATE(cprj_tmp,(dtset%natom,electronpositron%dimcprj))
     do iatom=1,dtset%natom;nlmn(iatom)=cprj(iatom,1)%nlmn;end do
     call pawcprj_alloc(cprj_tmp,cprj(1,1)%ncpgr,nlmn)
     ABI_DEALLOCATE(nlmn)
     call pawcprj_copy(electronpositron%cprj_ep,cprj_tmp)
     call pawcprj_copy(cprj,electronpositron%cprj_ep)
     call pawcprj_copy(cprj_tmp,cprj)
     call pawcprj_free(cprj_tmp)
     ABI_DATATYPE_DEALLOCATE(cprj_tmp)
    else
!TO BE ACTIVATED WHEN cprj IS PRESENT
!    call pawcprj_set_zero(cprj)
    end if
   end if

!  Eigenvalues
   if (electronpositron%dimeigen>0) then
    do ii=1,electronpositron%dimeigen
     rtmp=eigen(ii)
     eigen(ii)=electronpositron%eigen_ep(ii)
     electronpositron%eigen_ep(ii)=rtmp
    end do
   else
    eigen(:)=9.99999_dp
   end if

!  Occupations
   if (electronpositron%dimocc>0) then
    do ii=1,electronpositron%dimocc
     rtmp=occ(ii)
     occ(ii)=electronpositron%occ_ep(ii)
     electronpositron%occ_ep(ii)=rtmp
    end do
   else
    occ(:)=9.99999_dp
   end if

!  Forces
   if (allocated(electronpositron%fred_ep)) then
    do iatom=1,dtset%natom
     electronpositron%fred_ep(1:3,iatom)=fred(1:3,iatom)-electronpositron%fred_ep(1:3,iatom)
    end do
   end if

!  Stresses
   if (allocated(electronpositron%stress_ep)) then
    electronpositron%stress_ep(1:6)=stress(1:6)-electronpositron%stress_ep(1:6)
   end if

  end if
 end if

end subroutine exchange_electronpositron
!!***

!----------------------------------------------------------------------

!!****f* m_electronpositron/electronpositron_calctype
!! NAME
!!  electronpositron_calctype
!!
!! FUNCTION
!!  Returns the value of the calculation type from an electronpositron
!!  structure (can be eventually unassociated)
!!
!! INPUTS
!!  electronpositron=<type(electronpositron_type)>=electronpositron datastructure
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function electronpositron_calctype(electronpositron)

!Arguments ------------------------------------
!scalars
 type(electronpositron_type),pointer :: electronpositron

!************************************************************************

 if (associated(electronpositron)) then
  electronpositron_calctype=electronpositron%calctype
 else
  electronpositron_calctype=0
 end if


end function electronpositron_calctype
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/rhohxcpositron
!! NAME
!! rhohxcpositron
!!
!! FUNCTION
!! Calculate the electrons/positron correlation term for the positron
!!
!! INPUTS
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  mpi_enreg=information about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nhat(nfft,nspden*usepaw)= -PAW only- compensation density
!!  nkxc=second dimension of the array kxc, see rhotoxc.f for a description
!!  nspden=number of spin density components
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  paral_kgb=flag for (k,band,FFT) parallelism
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3.
!!  ucvol = unit cell volume (Bohr**3)
!!  usexcnhat= -PAW only- flag controling use of compensation density in Vxc
!!  usepaw=flag for PAW
!!  xccc3d(n3xccc)=3D core electron density for XC core correction (bohr^-3)
!!  xc_denpos= lowest allowed density (usually for the computation of the XC functionals)
!!
!! OUTPUT
!!  electronpositron%e_xc=electron-positron XC energy
!!  electronpositron%e_xcdc=Double-counting electron-positron XC energy
!!  strsxc(6)= contribution of xc to stress tensor (hartree/bohr^3),
!!  vhartr(nfft)=Hartree potential (returned if option/=0 and option/=10)
!!  vxcapn=XC electron-positron XC potential for the positron
!!  vxcavg=unit cell average of Vxc = (1/ucvol) Int [Vxc(r) d^3 r].
!!  kxcapn(nfft,nkxc)=electron-positron XC kernel (returned only if nkxc/=0)
!!
!! SIDE EFFECTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!
!! PARENTS
!!      energy,rhotov,setvtr
!!
!! CHILDREN
!!      mean_fftr,mkdenpos,xcden,xcpositron,xmpi_sum
!!
!! SOURCE

subroutine rhohxcpositron(electronpositron,gprimd,kxcapn,mpi_enreg,nfft,ngfft,nhat,nkxc,nspden,n3xccc,&
&                         paral_kgb,rhor,strsxc,ucvol,usexcnhat,usepaw,vhartr,vxcapn,vxcavg,xccc3d,xc_denpos)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,nkxc,nspden,n3xccc,paral_kgb,usexcnhat,usepaw
 real(dp),intent(in) :: ucvol,xc_denpos
 real(dp),intent(out) :: vxcavg
 type(electronpositron_type),pointer :: electronpositron
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(in) :: nhat(nfft,nspden*usepaw),rhor(nfft,nspden),xccc3d(n3xccc)
 real(dp),intent(out) :: kxcapn(nfft,nkxc),strsxc(6),vhartr(nfft),vxcapn(nfft,nspden)
 type(MPI_type),intent(in) :: mpi_enreg

!Local variables-------------------------------
!scalars
 integer :: cplex,ierr,ifft,ishift,iwarn,iwarnp,nfftot,ngr,ngrad,nspden_ep
 real(dp) :: exc,excdc,strdiag
 character(len=500) :: message
!arrays
 real(dp),parameter :: qphon(3)=(/0._dp,0._dp,0._dp/)
 real(dp) :: vxcavg_tmp(1)
 real(dp),allocatable :: fxcapn(:),grho2apn(:),rhoe(:,:,:),rhop(:,:),rhotote(:),vxc_ep(:),vxcgr_ep(:)

! *************************************************************************

 if (electronpositron_calctype(electronpositron)/=1) then
   message = 'Only electronpositron%calctype=1 allowed !'
   MSG_BUG(message)
 end if

 if (nkxc>3) then
   message = 'nkxc>3 (Kxc for GGA) not yet implemented !'
   MSG_ERROR(message)
 end if

!Hartree potential of the positron is zero
 vhartr=zero

!Some allocations/inits
 ngrad=1;if (electronpositron%ixcpositron==3.or.electronpositron%ixcpositron==31) ngrad=2
 ngr=0;if (ngrad==2) ngr=nfft
 ABI_ALLOCATE(fxcapn,(nfft))
 ABI_ALLOCATE(grho2apn,(ngr))
 nspden_ep=1;cplex=1;ishift=0
 iwarn=0;iwarnp=1

!Compute total electronic density
 ABI_ALLOCATE(rhotote,(nfft))
 rhotote(:)=electronpositron%rhor_ep(:,1)
 if (n3xccc>0) rhotote(:)=rhotote(:)+xccc3d(:)
 if (usepaw==1.and.usexcnhat==0) rhotote(:)=rhotote(:)-electronpositron%nhat_ep(:,1)

!Extra total electron/positron densities; compute gradients for GGA
 ABI_ALLOCATE(rhoe,(nfft,nspden_ep,ngrad**2))
 ABI_ALLOCATE(rhop,(nfft,nspden_ep))
 call xcden(cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden_ep,qphon,rhotote,rhoe)
 if (ngrad==2) grho2apn(:)=rhoe(:,1,2)**2+rhoe(:,1,3)**2+rhoe(:,1,4)**2
 rhop(:,1)=rhor(:,1);if (usepaw==1.and.usexcnhat==0) rhop(:,1)=rhop(:,1)-nhat(:,1)
 ABI_DEALLOCATE(rhotote)

!Make the densities positive
 call mkdenpos(iwarn ,nfft,nspden_ep,1,rhoe(:,1,1),xc_denpos)
 if (.not.electronpositron%posdensity0_limit) then
   call mkdenpos(iwarnp,nfft,nspden_ep,1,rhop,xc_denpos)
 end if

!Compute electron-positron Vxc_pos, Vxc_el, Fxc, Kxc, ...
 ABI_ALLOCATE(vxc_ep,(nfft))
 ABI_ALLOCATE(vxcgr_ep,(ngr))
 if (nkxc==0) then
   call xcpositron(fxcapn,grho2apn,electronpositron%ixcpositron,ngr,nfft,electronpositron%posdensity0_limit,&
&   rhoe(:,1,1),rhop(:,1),vxc_ep,vxcgr_ep,vxcapn)
 else
   call xcpositron(fxcapn,grho2apn,electronpositron%ixcpositron,ngr,nfft,electronpositron%posdensity0_limit,&
&   rhoe(:,1,1),rhop(:,1),vxc_ep,vxcgr_ep,vxcapn,dvxce=kxcapn)
 end if
 ABI_DEALLOCATE(rhoe)
 ABI_DEALLOCATE(vxc_ep)
 ABI_DEALLOCATE(vxcgr_ep)
 ABI_DEALLOCATE(grho2apn)

!Store Vxc and Kxc according to spin components
 if (nspden>=2) vxcapn(:,2)=vxcapn(:,1)
 if (nspden==4) vxcapn(:,3:4)=zero
 if (nkxc==3) then
   kxcapn(:,1)=two*kxcapn(:,1)
   kxcapn(:,2)=kxcapn(:,1)
   kxcapn(:,3)=kxcapn(:,1)
 end if

!Compute XC energies and contribution to stress tensor
 electronpositron%e_xc  =zero
 electronpositron%e_xcdc=zero
 strdiag=zero
 nfftot=PRODUCT(ngfft(1:3))
 do ifft=1,nfft
   electronpositron%e_xc  =electronpositron%e_xc  +fxcapn(ifft)
   electronpositron%e_xcdc=electronpositron%e_xcdc+vxcapn(ifft,1)*rhor(ifft,1)
!  strdiag=strdiag+fxcapn(ifft)   ! Already stored in rhotoxc !
   strdiag=strdiag-vxcapn(ifft,1)*rhop(ifft,1)
 end do
 if (usepaw==1.and.usexcnhat==0) then
   do ifft=1,nfft
     electronpositron%e_xcdc=electronpositron%e_xcdc-vxcapn(ifft,1)*nhat(ifft,1)
   end do
 end if
 electronpositron%e_xc  =electronpositron%e_xc  *ucvol/dble(nfftot)
 electronpositron%e_xcdc=electronpositron%e_xcdc*ucvol/dble(nfftot)
 strdiag=strdiag/dble(nfftot)
 ABI_DEALLOCATE(fxcapn)
 ABI_DEALLOCATE(rhop)

!Reduction in case of parallelism
 if(mpi_enreg%paral_kgb==1)then
   if(paral_kgb/=0)then
     exc=electronpositron%e_xc;excdc=electronpositron%e_xcdc
     call xmpi_sum(exc  ,mpi_enreg%comm_fft,ierr)
     call xmpi_sum(excdc,mpi_enreg%comm_fft,ierr)
     electronpositron%e_xc=exc;electronpositron%e_xcdc=excdc
     call xmpi_sum(strsxc,mpi_enreg%comm_fft,ierr)
   end if
 end if

!Store stress tensor
 strsxc(1:3)=strdiag
 strsxc(4:6)=zero

!Compute vxcavg
 call mean_fftr(vxcapn(:,1),vxcavg_tmp,nfft,nfftot,1)
 vxcavg=vxcavg_tmp(1)

end subroutine rhohxcpositron
!!***

END MODULE m_electronpositron
!!***
