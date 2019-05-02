!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_paw_ij
!! NAME
!!  m_paw_ij
!!
!! FUNCTION
!!  This module contains the definition of the paw_ij_type structured datatype,
!!  as well as related functions and methods.
!!  paw_ij_type variables contain various arrays given on (i,j) (partial waves) channels
!!  for a given atom.
!!
!! COPYRIGHT
!! Copyright (C) 2013-2019 ABINIT group (MT, FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  FOR DEVELOPPERS: in order to preserve the portability of libPAW library,
!!  please consult ~abinit/src/??_libpaw/libpaw-coding-rules.txt
!!
!! SOURCE

#include "libpaw.h"

MODULE m_paw_ij

 USE_DEFS
 USE_MSG_HANDLING
 USE_MPI_WRAPPERS
 USE_MEMORY_PROFILING

 use m_paral_atom, only : get_my_atmtab, free_my_atmtab, get_my_natom
 use m_pawtab,     only : pawtab_type
 use m_paw_io,     only : pawio_print_ij

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_paw_ij/paw_ij_type
!! NAME
!! paw_ij_type
!!
!! FUNCTION
!! For PAW, various arrays given on (i,j) (partial waves) channels
!!
!! SOURCE

 type,public :: paw_ij_type

!Integer scalars

  integer :: cplex_dij
   ! cplex_dij=1 if dij are real
   ! cplex_dij=2 if dij are complex (spin-orbit, non-collinear magnetism, magnetic field, ...)

  integer :: has_dij=0
   ! 1 if dij is allocated
   ! 2 if dij is already computed

  integer :: has_dij0=0
   ! 1 if dij0 is allocated
   ! 2 if dij0 is already computed

  integer :: has_dijexxc=0
   ! 1 if dijexxc is associated and used, 0 otherwise
   ! 2 if dijexxc is already computed

  integer :: has_dijfock=0
   ! 1 if dijfock is allocated
   ! 2 if dijfock is already computed

  integer :: has_dijfr=0
   ! 1 if dijfr is allocated
   ! 2 if dijfr is already computed

  integer :: has_dijhartree=0
   ! 1 if dijhartree is allocated
   ! 2 if dijhartree is already computed

  integer :: has_dijhat=0
   ! 1 if dijhat is allocated
   ! 2 if dijhat is already computed

  integer :: has_dijnd=0
   ! on site term due to nuclear dipole moment
   ! 1 if dijnd is associated and used, 0 otherwise
   ! 2 if dijnd is already computed

  integer :: has_dijso=0
   ! 1 if dijso is associated and used, 0 otherwise
   ! 2 if dijso is already computed

  integer :: has_dijU=0
   ! 1 if dijU is associated and used, 0 otherwise
   ! 2 if dijU is already computed

  integer :: has_dijxc=0
   ! 1 if dijxc is associated and used, 0 otherwise
   ! 2 if dijxc is already computed

  integer :: has_dijxc_hat=0
   ! 1 if dijxc_hat is associated and used, 0 otherwise
   ! 2 if dijxc_hat is already computed

  integer :: has_dijxc_val=0
   ! 1 if dijxc_val is associated and used, 0 otherwise
   ! 2 if dijxc_val is already computed

  integer :: has_exexch_pot=0
   ! 1 if PAW+(local exact exchange) potential is allocated

  integer :: has_pawu_occ=0
   ! 1 if PAW+U occupations are allocated

  integer :: itypat
   ! itypat=type of the atom

  integer :: lmn_size
   ! Number of (l,m,n) elements for the paw basis

  integer :: lmn2_size
   ! lmn2_size=lmn_size*(lmn_size+1)/2
   ! where lmn_size is the number of (l,m,n) elements for the paw basis

  integer :: ndij
   ! Number of components of dij
   ! Usually ndij=nspden, except for nspinor==2 (where ndij=nspinor**2)

  integer :: nspden
   ! Number of spin-density components (may be different from dtset%nspden if spin-orbit)

  integer :: nsppol
   ! Number of independant spin-components

  integer :: qphase
   ! qphase=2 if dij contain a exp(-i.q.r) phase (as in the q<>0 RF case), 1 if not
   ! (this may change the ij symmetry)

!Real (real(dp)) arrays

  real(dp), allocatable :: dij(:,:)
   ! dij(cplex_dij*qphase*lmn2_size,ndij)
   ! Dij term (non-local operator)
   ! May be complex if cplex_dij=2 or qphase=2
   ! ==== Storage for the 1st dimension ====
   ! For each klmn=ij:
   !   When Dij is complex (cplex_dij=2):
   !     dij(2*ij-1,:) contains the real part
   !     dij(2*ij  ,:) contains the imaginary part
   !   When a exp(-i.q.r) phase is included (qphase=2):
   !     dij(1:cplex_dij*lmn2_size,:)
   !         contains the real part of the phase, i.e. D_ij*cos(q.r)
   !     dij(cplex_dij*lmn2_size+1:2*cplex_dij*lmn2_size,:)
   !         contains the imaginary part of the phase, i.e. D_ij*sin(q.r)
   ! ==== Storage for the 2nd dimension ====
   !   dij(:,1) contains Dij^up-up
   !   dij(:,2) contains Dij^dn-dn
   !   dij(:,3) contains Dij^up-dn (only if nspinor=2)
   !   dij(:,4) contains Dij^dn-up (only if nspinor=2)

  real(dp), allocatable :: dij0(:)
   ! dij0(lmn2_size)
   ! Atomic part of Dij (read from PAW dataset)
   ! Same storage as Dij (see above); always real, spin independent

  real(dp), allocatable :: dijexxc(:,:)
   ! dijexxc(cplex_dij*lmn2_size,ndij)
   ! On-site matrix elements of the Fock operator (Local Exact exchange implementation)
   ! Same storage as Dij (see above); not available for RF (i.e. qphase=2)

  real(dp), allocatable :: dijfock(:,:)
   ! dijfock(cplex_dij*lmn2_size,ndij)
   ! Dij_fock term
   ! Contains all contributions to Dij from Fock exchange
   ! Same storage as Dij (see above); not available for RF (i.e. qphase=2)

  real(dp), allocatable :: dijfr(:,:)
   ! dijhat(cplex_dij*qphase*lmn2_size,ndij)
   ! For response function calculation only
   ! RF Frozen part of Dij (depends on q vector but not on 1st-order wave function)
   ! Same storage as Dij (see above)

  real(dp), allocatable :: dijhartree(:)
   ! dijhartree(qphase*lmn2_size)
   ! Dij_hartree term; contains all contributions to Dij from hartree
   ! Warning: dimensioned only by qphase (exp(-iqr)), not cplex_dij
   ! Same storage as Dij (see above); spin independent

  real(dp), allocatable :: dijhat(:,:)
   ! dijhat(cplex_dij*qphase*lmn2_size,ndij)
   ! Dij_hat term (non-local operator) i.e \sum_LM \int_FFT Q_{ij}^{LM} vtrial
   ! Same storage as Dij (see above)
   ! Same storage as Dij (see above)

  real(dp), allocatable :: dijnd(:,:)
   ! dijnd(cplex_dij*lmn2_size,ndij)
   ! On-site matrix elements of -\frac{1}{c}\mu\cdot L/r^3
   ! Same storage as Dij (see above); not available for RF (i.e. qphase=2)

  real(dp), allocatable :: dijso(:,:)
   ! dijso(cplex_dij*qphase*lmn2_size,ndij)
   ! On-site matrix elements of L.S i.e <phi_i|L.S|phi_j>
   ! Same storage as Dij (see above)
   ! Same storage as Dij (see above); not available for RF (i.e. qphase=2)

  real(dp), allocatable :: dijU(:,:)
   ! dijU(cplex_dij*qphase*lmn2_size,ndij)
   ! On-site matrix elements of the U part of the PAW Hamiltonian.
   ! Same storage as Dij (see above); not available for RF (i.e. qphase=2)

  real(dp), allocatable :: dijxc(:,:)
   ! dijxc(cplex_dij*qphase*lmn2_size,ndij)
   ! On-site matrix elements of vxc i.e
   !   <phi_i|vxc[n1+nc]|phi_j> - <tphi_i|vxc(tn1+nhat+tnc]|tphi_j>
   ! Same storage as Dij (see above)

  real(dp), allocatable :: dijxc_hat(:,:)
   ! dijxc_hat(cplex_dij*lmn2_size,ndij)
   ! Dij_hat term i.e \sum_LM \int_FFT Q_{ij}^{LM} Vxc
   ! Same storage as Dij (see above); not available for RF (i.e. qphase=2)

  real(dp), allocatable :: dijxc_val(:,:)
   ! dijxc_val(cplex_dij*lmn2_size,ndij)
   ! Onsite matrix elements of valence-only vxc i.e
   ! <phi_i|vxc[n1]|phi_j> - <tphi_i|vxc(tn1+nhat]|tphi_j>
   ! Same storage as Dij (see above); not available for RF (i.e. qphase=2)

  real(dp), allocatable :: noccmmp(:,:,:,:)
   ! noccmmp(cplex_dij,2*lpawu+1,2*lpawu+1,nocc_nspden)
   ! cplex_dij=1 if collinear
   ! cplex_dij=2 if spin orbit is used
   ! cplex_dij=2 is used if non-collinear (for coherence, it is not necessary in this case, however)
   ! gives occupation matrix for lda+u (computed in setnoccmmp)
   ! Stored as: noccmmp(:,:,1)=   n^{up,up}_{m,mp}
   !            noccmmp(:,:,2)=   n^{dn,dn}_{m,mp}
   !            noccmmp(:,:,3)=   n^{up,dn}_{m,mp}
   !            noccmmp(:,:,4)=   n^{dn,up}_{m,mp}
   ! noccmmp(m,mp,:) is computed from rhoij(klmn) with  m=klmntomn(2)>mp=klmntomn(1)

  real(dp), allocatable :: nocctot(:)
   ! nocctot(ndij)
   ! gives trace of occupation matrix for lda+u (computed in pawdenpot)
   ! for each value of ispden (1 or 2)

  real(dp), allocatable :: vpawx(:,:,:)
   ! vpawx(1,2*lexexch+1,nspden)
   ! exact exchange potential

 end type paw_ij_type

!public procedures.
 public :: paw_ij_init           ! Creation method
 public :: paw_ij_free           ! Free memory
 public :: paw_ij_nullify
 public :: paw_ij_copy           ! Copy object
 public :: paw_ij_print          ! Printout of the object
 public :: paw_ij_gather         ! MPI gather
 public :: paw_ij_redistribute   ! MPI redistribute
 public :: paw_ij_reset_flags    ! Resent the internal flags.

!private procedures.
 private :: paw_ij_isendreceive_getbuffer
 private :: paw_ij_isendreceive_fillbuffer
!!***

CONTAINS

!===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_paw_ij/paw_ij_init
!! NAME
!! paw_ij_init
!!
!! FUNCTION
!!  Initialize a Paw_ij data type.
!!
!! INPUTS
!!  cplex=1 if no phase is applied (GS), 2 if a exp(-iqr) phase is applied (Response Function at q<>0)
!!  natom=Number of atoms.
!!  ntypat=Number of types of atoms in cell.
!!  nspinor=number of spinor components
!!  nsppol=Number of independent spin polarizations.
!!  nspden=Number of spin-density components
!!  pawspnorb=1 if spin-orbit coupling is activated
!!  typat(natom)=Type of each atom
!!  Pawtab(ntypat)<type(pawtab_type)>=PAW tabulated starting data
!!
!! OPTIONAL INPUTS
!!  has_dij=1 to allocate Paw_ij%dij, 0 otherwise (default)
!!  has_dij0=1 to allocate Paw_ij%dij0, 0 otherwise (default)
!!  has_dijfr=1 to allocate Paw_ij%dijfr, 0 otherwise (default)
!!  has_dijhat=1 to allocate Paw_ij%dijhat, 0 otherwise (default)
!!  has_dijxc=1 to allocate Paw_ij%dijxc, 0 otherwise (default)
!!  has_dijxc_hat=1 to allocate Paw_ij%dijxc_hat, 0 otherwise (default)
!!  has_dijxc_val=1 to allocate Paw_ij%dijxc_val, 0 otherwise (default)
!!  has_dijhartree=1 to allocate Paw_ij%dijhartree, 0 otherwise (default)
!!  has_dijfock=1 to allocate Paw_ij%dijfock, 0 otherwise (default)
!!  has_dijnd=1 to allocate Paw_ij%dijnd, used only if some nucdipmom /= 0; otherwise 0 (default)
!!  has_dijso=1 to allocate Paw_ij%dijso, used only if pawspnorb>0. 0 otherwise (default)
!!  has_dijU=1 to allocate Paw_ij%dijU, used only if Pawtab(itypat)%usepawu/=0. 0 otherwise (default).
!!  has_dijexxc=to allocate Paw_ij%dijxx, 0 otherwise (default)
!!  has_exexch_pot=1 to allocate potential used in PAW+(local exact exchange) formalism, 0 otherwise (default)
!!  has_pawu_occ=1 to allocate occupations used in PAW+U formalism, 0 otherwise (default)
!!  nucdipmom(3,natom)= (optional) array of nuclear dipole moments at atomic sites
!!  mpi_atmtab(:)=indexes of the atoms treated by current proc
!!  comm_atom=MPI communicator over atoms
!!
!! OUTPUT
!!  Paw_ij(natom)<type(paw_ij_type)>=data structure containing PAW arrays given on (i,j) channels.
!!   In output all the basic dimensions are defined and the arrays are allocated
!!   according to the input variables.
!!
!! PARENTS
!!      bethe_salpeter,d2frnl,dfpt_nstpaw,dfpt_rhofermi,dfpt_scfcv,ldau_self
!!      m_energy,paw_qpscgw,respfn,scfcv,screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine paw_ij_init(Paw_ij,cplex,nspinor,nsppol,nspden,pawspnorb,natom,ntypat,typat,Pawtab,&
&                      has_dij,has_dij0,has_dijfock,has_dijfr,has_dijhartree,has_dijhat,& ! Optional
&                      has_dijxc,has_dijxc_hat,has_dijxc_val,has_dijnd,has_dijso,has_dijU,has_dijexxc,&  ! Optional
&                      has_exexch_pot,has_pawu_occ,nucdipmom,& ! Optional
&                      mpi_atmtab,comm_atom) ! optional arguments (parallelism)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nspinor,nspden,nsppol,natom,ntypat,pawspnorb
 integer,optional,intent(in) :: has_dij,has_dij0,has_dijfr,has_dijhat,has_dijxc,has_dijxc_hat,has_dijxc_val
 integer,optional,intent(in) :: has_dijnd,has_dijso,has_dijhartree,has_dijfock,has_dijU,has_dijexxc
 integer,optional,intent(in) :: has_exexch_pot,has_pawu_occ
 integer,optional,intent(in) :: comm_atom

!arrays
 integer,intent(in) :: typat(natom)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),optional,intent(in) :: nucdipmom(3,natom)
 type(Paw_ij_type),intent(inout) :: Paw_ij(:)
 type(Pawtab_type),intent(in) :: Pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: cplex_dij,iat,iat_tot,itypat,lmn2_size,my_comm_atom,my_natom,ndij,qphase
 logical :: my_atmtab_allocated,paral_atom,with_nucdipmom
!arrays
 integer,pointer :: my_atmtab(:)

! *************************************************************************

!@Paw_ij_type

 with_nucdipmom=.false.;if (present(nucdipmom)) with_nucdipmom=any(abs(nucdipmom)>tol8)

!Set up parallelism over atoms
 my_natom=size(paw_ij);if (my_natom==0) return
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

 do iat=1,my_natom
  iat_tot=iat;if (paral_atom) iat_tot=my_atmtab(iat)
  itypat=typat(iat_tot)

  cplex_dij=1;if ((nspinor==2).or.(with_nucdipmom)) cplex_dij=2
  qphase=cplex

  lmn2_size              =Pawtab(itypat)%lmn2_size
  Paw_ij(iat)%qphase     =qphase
  Paw_ij(iat)%cplex_dij  =cplex_dij
  Paw_ij(iat)%itypat     =itypat
  Paw_ij(iat)%nspden     =nspden
  Paw_ij(iat)%nsppol     =nsppol
  Paw_ij(iat)%lmn_size   =Pawtab(itypat)%lmn_size
  Paw_ij(iat)%lmn2_size  =lmn2_size
  Paw_ij(iat)%ndij       =MAX(nspinor**2,nspden)

  ndij=Paw_ij(iat)%ndij

  ! ==================================
  ! === Allocations (all optional) ===
  ! ==================================

  ! === Allocation for total Dij ===
  Paw_ij(iat)%has_dij=0
  if (PRESENT(has_dij)) then
    if (has_dij/=0) then
      Paw_ij(iat)%has_dij=1
      LIBPAW_ALLOCATE(Paw_ij(iat)%dij,(cplex_dij*qphase*lmn2_size,ndij))
      Paw_ij(iat)%dij(:,:)=zero
    end if
  end if

  ! === Allocation for atomic Dij ===
  Paw_ij(iat)%has_dij0=0
  if (PRESENT(has_dij0)) then
    if (has_dij0/=0) then
      Paw_ij(iat)%has_dij0=1
      LIBPAW_ALLOCATE(Paw_ij(iat)%dij0,(lmn2_size))
      Paw_ij(iat)%dij0(:)=zero
    end if
  end if

  ! === Allocation for Dij local exact exchange ===
  Paw_ij(iat)%has_dijexxc=0
  if (PRESENT(has_dijexxc)) then
    if (has_dijexxc/=0.and.Pawtab(itypat)%useexexch/=0) then
      Paw_ij(iat)%has_dijexxc=1
      LIBPAW_ALLOCATE(Paw_ij(iat)%dijexxc,(cplex_dij*lmn2_size,ndij))
      Paw_ij(iat)%dijexxc(:,:)=zero
    end if
  end if

 ! === Allocation for Dij_Fock ===
  Paw_ij(iat)%has_dijfock=0
  if (PRESENT(has_dijfock)) then
    if (has_dijfock/=0) then
      Paw_ij(iat)%has_dijfock=1
      LIBPAW_ALLOCATE(Paw_ij(iat)%dijfock,(cplex_dij*lmn2_size,ndij))
      Paw_ij(iat)%dijfock(:,:)=zero
    end if
  end if

  ! === Allocation for frozen part of 1st-order Dij ===
  Paw_ij(iat)%has_dijfr=0
  if (PRESENT(has_dijfr)) then
    if (has_dijfr/=0) then
      Paw_ij(iat)%has_dijfr=1
      LIBPAW_ALLOCATE(Paw_ij(iat)%dijfr,(cplex_dij*qphase*lmn2_size,ndij))
      Paw_ij(iat)%dijfr(:,:)=zero
    end if
  end if

  ! === Allocation for Dij_Hartree ===
  Paw_ij(iat)%has_dijhartree=0
  if (PRESENT(has_dijhartree)) then
    if (has_dijhartree/=0) then
      Paw_ij(iat)%has_dijhartree=1
      LIBPAW_ALLOCATE(Paw_ij(iat)%dijhartree,(qphase*lmn2_size))
      Paw_ij(iat)%dijhartree(:)=zero
    end if
  end if

  ! === Allocation for Dij_hat ===
  Paw_ij(iat)%has_dijhat=0
  if (PRESENT(has_dijhat)) then
    if (has_dijhat/=0) then
      Paw_ij(iat)%has_dijhat=1
      LIBPAW_ALLOCATE(Paw_ij(iat)%dijhat,(cplex_dij*qphase*lmn2_size,ndij))
      Paw_ij(iat)%dijhat(:,:)=zero
    end if
  end if

  ! === Allocation for Dij nuclear dipole moment ===
  Paw_ij(iat)%has_dijnd=0
  if (PRESENT(has_dijnd)) then
    if (has_dijnd/=0.and.with_nucdipmom) then
      Paw_ij(iat)%has_dijnd=1
      LIBPAW_ALLOCATE(Paw_ij(iat)%dijnd,(cplex_dij*lmn2_size,ndij))
      Paw_ij(iat)%dijnd(:,:)=zero
     end if
  end if

  ! === Allocation for Dij_SO ===
  Paw_ij(iat)%has_dijso=0
  if (PRESENT(has_dijso)) then
    if (has_dijso/=0.and.pawspnorb>0) then
      Paw_ij(iat)%has_dijso=1
      LIBPAW_ALLOCATE(Paw_ij(iat)%dijso,(cplex_dij*qphase*lmn2_size,ndij))
      Paw_ij(iat)%dijso(:,:)=zero
     end if
  end if

  ! === Allocation for Dij_U_val ===
  Paw_ij(iat)%has_dijU=0
  if (PRESENT(has_dijU)) then
    if (has_dijU/=0) then
      Paw_ij(iat)%has_dijU=1
      LIBPAW_ALLOCATE(Paw_ij(iat)%dijU,(cplex_dij*qphase*lmn2_size,ndij))
       Paw_ij(iat)%dijU(:,:)=zero
    end if
  end if

  ! === Allocation for total Dij_XC ===
  Paw_ij(iat)%has_dijxc=0
  if (PRESENT(has_dijxc)) then
    if (has_dijxc/=0) then
      Paw_ij(iat)%has_dijxc=1
      LIBPAW_ALLOCATE(Paw_ij(iat)%dijxc,(cplex_dij*qphase*lmn2_size,ndij))
      Paw_ij(iat)%dijxc(:,:)=zero
    end if
  end if

  ! === Allocation for total Dij_XC_hat ===
  Paw_ij(iat)%has_dijxc_hat=0
  if (PRESENT(has_dijxc_hat)) then
    if (has_dijxc_hat/=0) then
      Paw_ij(iat)%has_dijxc_hat=1
      LIBPAW_ALLOCATE(Paw_ij(iat)%dijxc_hat,(cplex_dij*lmn2_size,ndij))
     Paw_ij(iat)%dijxc_hat(:,:)=zero
    end if
  end if

  ! === Allocation for total Dij_XC_val ===
  Paw_ij(iat)%has_dijxc_val=0
  if (PRESENT(has_dijxc_val)) then
    if (has_dijxc_val/=0) then
      Paw_ij(iat)%has_dijxc_val=1
      LIBPAW_ALLOCATE(Paw_ij(iat)%dijxc_val,(cplex_dij*lmn2_size,ndij))
      Paw_ij(iat)%dijxc_val(:,:)=zero
    end if
  end if

  ! === Allocation for PAW+U occupations ===
  Paw_ij(iat)%has_pawu_occ=0
  if (PRESENT(has_pawu_occ)) then
    if (has_pawu_occ/=0.and.Pawtab(itypat)%usepawu/=0) then
      Paw_ij(iat)%has_pawu_occ=1
      LIBPAW_ALLOCATE(Paw_ij(iat)%noccmmp,(cplex_dij,2*Pawtab(itypat)%lpawu+1,2*Pawtab(itypat)%lpawu+1,ndij))
      LIBPAW_ALLOCATE(Paw_ij(iat)%nocctot,(ndij))
     Paw_ij(iat)%noccmmp(:,:,:,:)=zero
     Paw_ij(iat)%nocctot(:)=zero
    end if
  end if

  ! === Allocation for PAW+LEXX potential ===
  Paw_ij(iat)%has_exexch_pot=0
  if (PRESENT(has_exexch_pot)) then
    if (has_exexch_pot/=0.and.Pawtab(itypat)%useexexch/=0) then
      Paw_ij(iat)%has_exexch_pot=1
    ! TODO solve issue with first dimension
      LIBPAW_ALLOCATE(Paw_ij(iat)%vpawx,(1,lmn2_size,nspden))
      Paw_ij(iat)%vpawx(:,:,:)=zero
     end if
  end if

 end do

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

end subroutine paw_ij_init
!!***

!----------------------------------------------------------------------

!!****f* m_paw_ij/paw_ij_free
!! NAME
!!  paw_ij_free
!!
!! FUNCTION
!!  Deallocate pointers and nullify flags in a paw_ij structure
!!
!! SIDE EFFECTS
!!  paw_ij(:)<type(paw_ij_type)>=paw arrays given on (i,j) channels
!!
!! PARENTS
!!      bethe_salpeter,d2frnl,dfpt_nstpaw,dfpt_rhofermi,dfpt_scfcv,ldau_self
!!      m_energy,m_paral_pert,m_paw_ij,pawprt,respfn,scfcv,screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine paw_ij_free(Paw_ij)

!Arguments ------------------------------------
!arrays
 type(Paw_ij_type),intent(inout) :: Paw_ij(:)

!Local variables-------------------------------
 integer :: iat,natom

! *************************************************************************

!@Paw_ij_type

 natom=SIZE(Paw_ij);if (natom==0) return

 do iat=1,natom
  if (allocated(Paw_ij(iat)%dij       ))  then
    LIBPAW_DEALLOCATE(Paw_ij(iat)%dij)
  end if
  if (allocated(Paw_ij(iat)%dij0      ))  then
    LIBPAW_DEALLOCATE(Paw_ij(iat)%dij0)
  end if
  if (allocated(Paw_ij(iat)%dijexxc   ))  then
    LIBPAW_DEALLOCATE(Paw_ij(iat)%dijexxc)
  end if
  if (allocated(Paw_ij(iat)%dijfock   ))  then
    LIBPAW_DEALLOCATE(Paw_ij(iat)%dijfock)
  end if
  if (allocated(Paw_ij(iat)%dijfr     ))  then
    LIBPAW_DEALLOCATE(Paw_ij(iat)%dijfr)
  end if
  if (allocated(Paw_ij(iat)%dijhartree))  then
    LIBPAW_DEALLOCATE(Paw_ij(iat)%dijhartree)
  end if
  if (allocated(Paw_ij(iat)%dijhat    ))  then
    LIBPAW_DEALLOCATE(Paw_ij(iat)%dijhat)
  end if
  if (allocated(Paw_ij(iat)%dijnd     ))  then
    LIBPAW_DEALLOCATE(Paw_ij(iat)%dijnd)
  end if
  if (allocated(Paw_ij(iat)%dijU      ))  then
    LIBPAW_DEALLOCATE(Paw_ij(iat)%dijU)
  end if
  if (allocated(Paw_ij(iat)%dijso     ))  then
    LIBPAW_DEALLOCATE(Paw_ij(iat)%dijso)
  end if
  if (allocated(Paw_ij(iat)%dijxc     ))  then
    LIBPAW_DEALLOCATE(Paw_ij(iat)%dijxc)
  end if
  if (allocated(Paw_ij(iat)%dijxc_hat ))  then
    LIBPAW_DEALLOCATE(Paw_ij(iat)%dijxc_hat)
  end if
  if (allocated(Paw_ij(iat)%dijxc_val ))  then
    LIBPAW_DEALLOCATE(Paw_ij(iat)%dijxc_val)
  end if
  if (allocated(Paw_ij(iat)%noccmmp   ))  then
    LIBPAW_DEALLOCATE(Paw_ij(iat)%noccmmp)
  end if
  if (allocated(Paw_ij(iat)%nocctot   ))  then
    LIBPAW_DEALLOCATE(Paw_ij(iat)%nocctot)
  end if
  if (allocated(Paw_ij(iat)%vpawx     ))  then
    LIBPAW_DEALLOCATE(Paw_ij(iat)%vpawx)
  end if

  ! === Reset all has_* flags ===
  Paw_ij(iat)%has_dij       =0
  Paw_ij(iat)%has_dij0      =0
  Paw_ij(iat)%has_dijexxc   =0
  Paw_ij(iat)%has_dijfock   =0
  Paw_ij(iat)%has_dijfr     =0
  Paw_ij(iat)%has_dijhartree=0
  Paw_ij(iat)%has_dijhat    =0
  Paw_ij(iat)%has_dijnd     =0
  Paw_ij(iat)%has_dijso     =0
  Paw_ij(iat)%has_dijU      =0
  Paw_ij(iat)%has_dijxc     =0
  Paw_ij(iat)%has_dijxc_hat =0
  Paw_ij(iat)%has_dijxc_val =0
  Paw_ij(iat)%has_exexch_pot=0
  Paw_ij(iat)%has_pawu_occ  =0
 end do

 !call paw_ij_nullify(Paw_ij)

end subroutine paw_ij_free
!!***

!----------------------------------------------------------------------

!!****f* m_paw_ij/paw_ij_nullify
!! NAME
!!  paw_ij_nullify
!!
!! FUNCTION
!!  Reset all flags in a paw_ij structure
!!
!! SIDE EFFECTS
!!  Paw_ij(:)<type(paw_ij_type)>=PAW arrays given on (i,j) channels.
!!
!! PARENTS
!!      bethe_salpeter,d2frnl,dfpt_nstpaw,dfpt_rhofermi,dfpt_scfcv,ldau_self
!!      m_energy,m_paw_ij,paw_qpscgw,pawprt,respfn,scfcv,screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine paw_ij_nullify(Paw_ij)

!Arguments ------------------------------------
!arrays
 type(Paw_ij_type),intent(inout) :: Paw_ij(:)

!Local variables-------------------------------
 integer :: iat,natom

! *************************************************************************

 !@Paw_ij_type

 ! MGPAW: This one could be removed/renamed,
 ! variables can be initialized in the datatype declaration
 ! Do we need to expose this in the public API?

 natom=SIZE(Paw_ij(:));if (natom==0) return

 ! Set all has_* flags to zero.
 do iat=1,natom
   ! === Set all has_* flags to zero ===
   Paw_ij(iat)%has_dij       =0
   Paw_ij(iat)%has_dij0      =0
   Paw_ij(iat)%has_dijexxc   =0
   Paw_ij(iat)%has_dijfock   =0
   Paw_ij(iat)%has_dijfr     =0
   Paw_ij(iat)%has_dijhartree=0
   Paw_ij(iat)%has_dijhat    =0
   Paw_ij(iat)%has_dijnd     =0
   Paw_ij(iat)%has_dijso     =0
   Paw_ij(iat)%has_dijU      =0
   Paw_ij(iat)%has_dijxc     =0
   Paw_ij(iat)%has_dijxc_hat =0
   Paw_ij(iat)%has_dijxc_val =0
   Paw_ij(iat)%has_exexch_pot=0
   Paw_ij(iat)%has_pawu_occ  =0
 end do !iat

end subroutine paw_ij_nullify
!!***

!----------------------------------------------------------------------

!!****f* m_paw_ij/paw_ij_copy
!! NAME
!!  paw_ij_copy
!!
!! FUNCTION
!!  Copy one paw_ij datastructure into another
!!  Can take into accound changes of dimensions
!!  Can copy a shared paw_ij into distributed ones (when parallelism is activated)
!!
!! INPUTS
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  paw_ij_in(:)<type(paw_ij_type)>= input paw_ij datastructure
!!
!! SIDE EFFECTS
!!  paw_ij_cpy(:)<type(paw_ij_type)>= output paw_ij datastructure
!!
!! NOTES
!!  paw_ij_cpy must have been allocated in the calling function.
!!
!! PARENTS
!!      m_paw_ij
!!
!! CHILDREN
!!
!! SOURCE

subroutine paw_ij_copy(paw_ij_in,paw_ij_cpy, &
&                      mpi_atmtab,comm_atom) ! optional arguments (parallelism)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: comm_atom
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 type(paw_ij_type),intent(in) :: paw_ij_in(:)
 type(paw_ij_type),intent(inout),target :: paw_ij_cpy(:)

!Local variables-------------------------------
!scalars
integer :: ij,ij1,my_comm_atom,my_paw_ij,npaw_ij_in,npaw_ij_max,npaw_ij_out,paral_case,sz1,sz2,sz3,sz4
logical :: my_atmtab_allocated,paral_atom
character(len=500) :: msg
!arrays
 integer,pointer :: my_atmtab(:)
 type(paw_ij_type),pointer :: paw_ij_out(:)

! *************************************************************************

!@Paw_ij_type

!Retrieve sizes
 npaw_ij_in=size(paw_ij_in);npaw_ij_out=size(paw_ij_cpy)

!Set up parallelism over atoms
 paral_atom=(present(comm_atom));if (paral_atom) paral_atom=(xmpi_comm_size(comm_atom)>1)
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 my_atmtab_allocated=.false.

!Determine in which case we are (parallelism, ...)
!No parallelism: a single copy operation
 paral_case=0;npaw_ij_max=npaw_ij_in
 paw_ij_out => paw_ij_cpy
 if (paral_atom) then
   if (npaw_ij_out<npaw_ij_in) then ! Parallelism: the copy operation is a scatter
     call get_my_natom(my_comm_atom,my_paw_ij,npaw_ij_in)
     if (my_paw_ij==npaw_ij_out) then
       call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,npaw_ij_in)
       paral_case=1;npaw_ij_max=npaw_ij_out
       paw_ij_out => paw_ij_cpy
     else
       msg=' npaw_ij_out should be equal to my_paw_ij !'
       MSG_BUG(msg)
     end if
   else                            ! Parallelism: the copy operation is a gather
     call get_my_natom(my_comm_atom,my_paw_ij,npaw_ij_out)
     if (my_paw_ij==npaw_ij_in) then
       paral_case=2;npaw_ij_max=npaw_ij_in
     else
       msg=' npaw_ij_in should be equal to my_paw_ij !'
       MSG_BUG(msg)
     end if
   end if
 end if

!First case: a simple copy or a scatter
 if (npaw_ij_max>0.and.((paral_case==0).or.(paral_case==1))) then
   call paw_ij_nullify(paw_ij_out)

!  Loop on paw_ij components
   do ij1=1,npaw_ij_max
     ij=ij1; if (paral_case==1) ij=my_atmtab(ij1)

     paw_ij_out(ij1)%qphase=paw_ij_in(ij)%qphase
     paw_ij_out(ij1)%cplex_dij=paw_ij_in(ij)%cplex_dij
     paw_ij_out(ij1)%has_dij=paw_ij_in(ij)%has_dij
     paw_ij_out(ij1)%has_dij0=paw_ij_in(ij)%has_dij0
     paw_ij_out(ij1)%has_dijexxc=paw_ij_in(ij)%has_dijexxc
     paw_ij_out(ij1)%has_dijfock=paw_ij_in(ij)%has_dijfock
     paw_ij_out(ij1)%has_dijfr=paw_ij_in(ij)%has_dijfr
     paw_ij_out(ij1)%has_dijhartree=paw_ij_in(ij)%has_dijhartree
     paw_ij_out(ij1)%has_dijhat=paw_ij_in(ij)%has_dijhat
     paw_ij_out(ij1)%has_dijnd=paw_ij_in(ij)%has_dijnd
     paw_ij_out(ij1)%has_dijso=paw_ij_in(ij)%has_dijso
     paw_ij_out(ij1)%has_dijU=paw_ij_in(ij)%has_dijU
     paw_ij_out(ij1)%has_dijxc=paw_ij_in(ij)%has_dijxc
     paw_ij_out(ij1)%has_dijxc_hat=paw_ij_in(ij)%has_dijxc_hat
     paw_ij_out(ij1)%has_dijxc_val=paw_ij_in(ij)%has_dijxc_val
     paw_ij_out(ij1)%has_exexch_pot=paw_ij_in(ij)%has_exexch_pot
     paw_ij_out(ij1)%has_pawu_occ=paw_ij_in(ij)%has_pawu_occ
     paw_ij_out(ij1)%itypat=paw_ij_in(ij)%itypat
     paw_ij_out(ij1)%lmn_size=paw_ij_in(ij)%lmn_size
     paw_ij_out(ij1)%lmn2_size=paw_ij_in(ij)%lmn2_size
     paw_ij_out(ij1)%ndij=paw_ij_in(ij)%ndij
     paw_ij_out(ij1)%nspden=paw_ij_in(ij)%nspden
     paw_ij_out(ij1)%nsppol=paw_ij_in(ij)%nsppol
     if (paw_ij_in(ij)%has_dij>=1) then
       sz1=size(paw_ij_in(ij)%dij,1);sz2=size(paw_ij_in(ij)%dij,2)
       LIBPAW_ALLOCATE(paw_ij_out(ij1)%dij,(sz1,sz2))
       if (paw_ij_in(ij)%has_dij==2) &
&          paw_ij_out(ij1)%dij(:,:)=paw_ij_in(ij)%dij(:,:)
     end if
     if (paw_ij_in(ij)%has_dij0>=1) then
       sz1=size(paw_ij_in(ij)%dij0,1)
       LIBPAW_ALLOCATE(paw_ij_out(ij1)%dij0,(sz1))
      if (paw_ij_in(ij)%has_dij0 ==2) &
&         paw_ij_out(ij1)%dij0(:)=paw_ij_in(ij)%dij0(:)
     end if
     if (paw_ij_in(ij)%has_dijexxc>=1) then
       sz1=size(paw_ij_in(ij)%dijexxc,1);sz2=size(paw_ij_in(ij)%dijexxc,2)
       LIBPAW_ALLOCATE(paw_ij_out(ij1)%dijexxc,(sz1,sz2))
       if (paw_ij_in(ij)%has_dijexxc==2) &
&          paw_ij_out(ij1)%dijexxc(:,:)=paw_ij_in(ij)%dijexxc(:,:)
     end if
     if (paw_ij_in(ij)%has_dijfock>=1) then
       sz1=size(paw_ij_in(ij)%dijfock,1);sz2=size(paw_ij_in(ij)%dijfock,2)
       LIBPAW_ALLOCATE(paw_ij_out(ij1)%dijfock,(sz1,sz2))
       if (paw_ij_in(ij)%has_dijfock==2) &
&          paw_ij_out(ij1)%dijfock(:,:)=paw_ij_in(ij)%dijfock(:,:)
     end if
     if (paw_ij_in(ij)%has_dijfr>=1) then
       sz1=size(paw_ij_in(ij)%dijfr,1);sz2=size(paw_ij_in(ij)%dijfr,2)
       LIBPAW_ALLOCATE(paw_ij_out(ij1)%dijfr,(sz1,sz2))
       if (paw_ij_in(ij)%has_dijfr==2) &
&          paw_ij_out(ij1)%dijfr(:,:)=paw_ij_in(ij)%dijfr(:,:)
     end if
     if (paw_ij_in(ij)%has_dijhartree>=1) then
       sz1=size(paw_ij_in(ij)%dijhartree,1)
       LIBPAW_ALLOCATE(paw_ij_out(ij1)%dijhartree,(sz1))
       if (paw_ij_in(ij)%has_dijhartree==2) &
&          paw_ij_out(ij1)%dijhartree(:)=paw_ij_in(ij)%dijhartree(:)
     end if
     if (paw_ij_in(ij)%has_dijhat>=1) then
       sz1=size(paw_ij_in(ij)%dijhat,1);sz2=size(paw_ij_in(ij)%dijhat,2)
       LIBPAW_ALLOCATE(paw_ij_out(ij1)%dijhat,(sz1,sz2))
       if (paw_ij_in(ij)%has_dijhat==2) &
 &         paw_ij_out(ij1)%dijhat(:,:)=paw_ij_in(ij)%dijhat(:,:)
     end if
     if (paw_ij_in(ij)%has_dijnd>=1) then
       sz1=size(paw_ij_in(ij)%dijnd,1);sz2=size(paw_ij_in(ij)%dijnd,2)
       LIBPAW_ALLOCATE(paw_ij_out(ij1)%dijnd,(sz1,sz2))
       if (paw_ij_in(ij)%has_dijnd==2) &
&          paw_ij_out(ij1)%dijnd(:,:)=paw_ij_in(ij)%dijnd(:,:)
     end if
     if (paw_ij_in(ij)%has_dijU>=1) then
       sz1=size(paw_ij_in(ij)%dijU,1);sz2=size(paw_ij_in(ij)%dijU,2)
       LIBPAW_ALLOCATE(paw_ij_out(ij1)%dijU,(sz1,sz2))
       if (paw_ij_in(ij)%has_dijU==2) &
&          paw_ij_out(ij1)%dijU(:,:)=paw_ij_in(ij)%dijU(:,:)
     end if
     if (paw_ij_in(ij)%has_dijso>=1) then
       sz1=size(paw_ij_in(ij)%dijso,1);sz2=size(paw_ij_in(ij)%dijso,2)
       LIBPAW_ALLOCATE(paw_ij_out(ij1)%dijso,(sz1,sz2))
       if (paw_ij_in(ij)%has_dijso==2) &
&          paw_ij_out(ij1)%dijso(:,:)=paw_ij_in(ij)%dijso(:,:)
     end if
     if (paw_ij_in(ij)%has_dijxc>=1) then
       sz1=size(paw_ij_in(ij)%dijxc,1);sz2=size(paw_ij_in(ij)%dijxc,2)
       LIBPAW_ALLOCATE(paw_ij_out(ij1)%dijxc,(sz1,sz2))
       if (paw_ij_in(ij)%has_dijxc==2) &
&          paw_ij_out(ij1)%dijxc(:,:)=paw_ij_in(ij)%dijxc(:,:)
     end if
     if (paw_ij_in(ij)%has_dijxc_hat>=1) then
       sz1=size(paw_ij_in(ij)%dijxc_hat,1);sz2=size(paw_ij_in(ij)%dijxc_hat,2)
       LIBPAW_ALLOCATE(paw_ij_out(ij1)%dijxc_hat,(sz1,sz2))
       if (paw_ij_in(ij)%has_dijxc_hat==2) &
&          paw_ij_out(ij1)%dijxc_hat(:,:)=paw_ij_in(ij)%dijxc_hat(:,:)
     end if
     if (paw_ij_in(ij)%has_dijxc_val>=1) then
       sz1=size(paw_ij_in(ij)%dijxc_val,1);sz2=size(paw_ij_in(ij)%dijxc_val,2)
       LIBPAW_ALLOCATE(paw_ij_out(ij1)%dijxc_val,(sz1,sz2))
       if (paw_ij_in(ij)%has_dijxc_val==2) &
&          paw_ij_out(ij1)%dijxc_val(:,:)=paw_ij_in(ij)%dijxc_val(:,:)
     end if
     if (paw_ij_in(ij)%has_pawu_occ>=1) then
       sz1=size(paw_ij_in(ij)%noccmmp,1);sz2=size(paw_ij_in(ij)%noccmmp,2)
       sz3=size(paw_ij_in(ij)%noccmmp,3);sz4=size(paw_ij_in(ij)%noccmmp,4)
       LIBPAW_ALLOCATE(paw_ij_out(ij1)%noccmmp,(sz1,sz2,sz3,sz4))
       sz1=size(paw_ij_in(ij)%nocctot,1)
       LIBPAW_ALLOCATE(paw_ij_out(ij1)%nocctot,(sz1))
       if (paw_ij_in(ij)%has_pawu_occ==2) then
         paw_ij_out(ij1)%noccmmp(:,:,:,:)=paw_ij_in(ij)%noccmmp(:,:,:,:)
         paw_ij_out(ij1)%nocctot(:)=paw_ij_in(ij)%nocctot(:)
        end if
     end if
     if (paw_ij_in(ij)%has_exexch_pot >= 1) then
       sz1=size(paw_ij_in(ij)%vpawx,1);sz2=size(paw_ij_in(ij)%vpawx,2)
       sz3=size(paw_ij_in(ij)%vpawx,3)
       LIBPAW_ALLOCATE(paw_ij_out(ij1)%vpawx,(sz1,sz2,sz3))
       if (paw_ij_in(ij)%has_exexch_pot==2) &
&          paw_ij_out(ij1)%vpawx(:,:,:)=paw_ij_in(ij)%vpawx(:,:,:)
     end if

   end do
 end if

!Second case: a gather
 if (paral_case==2) then
   call paw_ij_gather(paw_ij_in,paw_ij_cpy,-1,my_comm_atom)
 end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

end subroutine paw_ij_copy
!!***

!----------------------------------------------------------------------

!!****f* m_paw_ij/paw_ij_print
!! NAME
!! paw_ij_print
!!
!! FUNCTION
!!  Print out the content of a paw_ij datastructure (Dij only)
!!
!! INPUTS
!!  [enunit]=governs the units to be used for the output of total Dij (0:Ha, 1:Ha+eV)
!!  [ipert]=only for DFPT: index of the perturbation (0 for ground state)
!!  [unit]=the unit number for output
!!  [pawprtvol]=verbosity level
!!  [pawspnorb]=1 if spin-orbit coupling is activated
!!  [mode_paral]=either "COLL" or "PERS"
!!  [mpi_atmtab(:)]=indexes of the atoms treated by current proc (can be computed here)
!!  [comm_atom]=MPI communicator over atoms (needed if parallelism over atoms is activated)
!!  [natom]=total number of atom (needed if parallelism over atoms is activated)
!!          if Paw_ij is distributed, natom is different from size(Paw_ij).
!!
!! OUTPUT
!! (Only writing)
!!
!! NOTES
!!
!! PARENTS
!!      m_pawdij,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine paw_ij_print(Paw_ij,unit,pawprtvol,pawspnorb,mode_paral,enunit,ipert, &
&                       mpi_atmtab,comm_atom,natom)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: enunit,ipert
 integer,optional,intent(in) :: comm_atom,natom
 integer,optional,intent(in) :: pawprtvol,pawspnorb
 integer,optional,intent(in) :: unit
 character(len=4),optional,intent(in) :: mode_paral
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 type(Paw_ij_type),target,intent(in) :: Paw_ij(:)

!Local variables-------------------------------
 character(len=7),parameter :: dspin(6)=(/"up     ","down   ","up-up  ","dwn-dwn","up-dwn ","dwn-up "/)
!scalars
 integer :: cplex_dij,iatom,iatom_tot,idij,idij_sym,lmn2_size,lmn_size,my_comm_atom,my_natom,nspden !klmn,
 integer :: nsploop,nsppol,my_unt,ndij,qphase,tmp_cplex_dij,my_ipert,my_enunit,my_prtvol,size_paw_ij
 logical :: my_atmtab_allocated,paral_atom
 character(len=4) :: my_mode
 character(len=2000) :: msg
!arrays
 integer :: idum(0)
 integer,pointer :: my_atmtab(:)
 real(dp),allocatable,target :: dij(:),dijs(:),dijh(:,:)
 real(dp),pointer :: dij2p(:),dij2p_(:)

! *************************************************************************

 if (.False.) write(std_out,*)"pawspnorb:",pawspnorb

!@Paw_ij_type
 size_paw_ij=SIZE(Paw_ij);if (size_paw_ij==0) return

 my_unt   =std_out   ; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0         ; if (PRESENT(pawprtvol )) my_prtvol=pawprtvol
 my_mode  ='COLL'    ; if (PRESENT(mode_paral)) my_mode  =mode_paral
 my_ipert =0         ; if (PRESENT(ipert))      my_ipert =ipert
 my_enunit=0         ; if (PRESENT(enunit))     my_enunit=enunit
 my_natom=size_paw_ij; if (PRESENT(natom))      my_natom=natom

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.my_natom/=size_paw_ij)
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,my_natom,my_natom_ref=size_paw_ij)

 if (abs(my_prtvol)>=1) then
   if (my_ipert==0) then
     write(msg,'(4a)')ch10,' ==== Values of psp strength Dij (Hartree) ============'
   else
     write(msg,'(4a)')ch10,' ==== Values of psp strength Dij(1) (Hartree) ========='
   end if
   call wrtout(my_unt,msg,my_mode)
 end if

 nsppol = Paw_ij(1)%nsppol
 nspden = Paw_ij(1)%nspden
 nsploop= nsppol; if (Paw_ij(1)%ndij==4) nsploop=4

 do iatom=1,size_paw_ij

  iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)

  lmn_size  = Paw_ij(iatom)%lmn_size
  lmn2_size = Paw_ij(iatom)%lmn2_size
  cplex_dij = Paw_ij(iatom)%cplex_dij
  qphase    = Paw_ij(iatom)%qphase
  ndij      = Paw_ij(iatom)%ndij

  ! ====================================
  ! === Loop over density components ===
  ! ====================================
  do idij=1,nsploop

   idij_sym=idij;if (ndij==4.and.idij>2) idij_sym=7-idij
   if (qphase==2) then
     LIBPAW_ALLOCATE(dij,(2*lmn2_size))
     LIBPAW_ALLOCATE(dijs,(2*lmn2_size))
   end if

!  =================== Detailed output =====================================
   if (ABS(my_prtvol)>=1.and.(iatom_tot==1.or.iatom_tot==my_natom.or.my_prtvol<0)) then

     !Title
     if (nspden==2.and.nsppol==1) then
       write(msg,'(2a,i3,3a)')ch10,&
&       ' >>>>>>>>>> Atom ',iatom_tot,':',ch10,&
&       ' (antiferromagnetism case: only one spin component)'
     else if (paw_ij(iatom)%ndij==1) then
       write(msg, '(2a,i3,a)') ch10,&
&      ' >>>>>>>>>> Atom ',iatom_tot,':'
     else
       write(msg,'(2a,i3,3a)') ch10,&
&       ' >>>>>>>>>> Atom ',iatom_tot,' (component ',TRIM(dspin(idij+2*(nsploop/4))),'):'
     end if
     call wrtout(my_unt,msg,my_mode)

     !Dij atomic
     if (Paw_ij(iatom)%has_dij0/=0.and.idij<=2.and.my_ipert<=0) then
       write(msg,'(a)') '   ************ Dij atomic (Dij0) ***********'
       call wrtout(my_unt,msg,my_mode)
       call pawio_print_ij(my_unt,Paw_ij(iatom)%dij0,lmn2_size,1,lmn_size,-1,idum,0,my_prtvol,idum,-1.d0,1,&
&                   opt_sym=2,mode_paral=my_mode)
     end if

     !Dij Local Exact Exchange
     if (Paw_ij(iatom)%has_dijexxc/=0.and.(idij<=2.or.nspden==4).and.my_ipert<=0) then
       write(msg,'(a)') '   ************* Dij_Local Exact exchange **********'
       call wrtout(my_unt,msg,my_mode)
       call get_dij_parts(cplex_dij,1,Paw_ij(iatom)%dijexxc)
       call pawio_print_ij(my_unt,dij2p,lmn2_size,tmp_cplex_dij,lmn_size,-1,idum,0,&
&                   my_prtvol,idum,-1.d0,1,opt_sym=2,asym_ij=dij2p_,mode_paral=my_mode)
     end if

     !Dij Fock
     if (Paw_ij(iatom)%has_dijfock/=0.and.(idij<=2.or.nspden==4).and.my_ipert<=0) then
       write(msg,'(a)') '   ************* Dij_Fock **********'
       call wrtout(my_unt,msg,my_mode)
       call get_dij_parts(cplex_dij,1,Paw_ij(iatom)%dijfock)
       call pawio_print_ij(my_unt,dij2p,lmn2_size,tmp_cplex_dij,lmn_size,-1,idum,0,&
&                   my_prtvol,idum,-1.d0,1,opt_sym=2,asym_ij=dij2p_,mode_paral=my_mode)
     end if

     !Dij Frozen (RF)
     if (Paw_ij(iatom)%has_dijfr/=0.and.(idij<=2.or.nspden==4).and.my_ipert>0) then
       write(msg,'(a)') '   ************** Dij(1) Frozen **************'
       call wrtout(my_unt,msg,my_mode)
       call get_dij_parts(cplex_dij,qphase,Paw_ij(iatom)%dijfr)
       call pawio_print_ij(my_unt,dij2p,lmn2_size,tmp_cplex_dij,lmn_size,-1,idum,0,&
&                   my_prtvol,idum,-1.d0,1,opt_sym=2,asym_ij=dij2p_,mode_paral=my_mode)
     end if

     !Dij Hartree
     if (Paw_ij(iatom)%has_dijhartree/=0.and.idij<=2) then
       if (my_ipert==0) then
         write(msg,'(a)') '   ************** Dij Hartree ***************'
       else
         write(msg,'(a)') '   ************* Dij(1) Hartree *************'
       end if
       call wrtout(my_unt,msg,my_mode)
       LIBPAW_ALLOCATE(dijh,(qphase*lmn2_size,1))
       dijh(:,1)=Paw_ij(iatom)%dijhartree(:)
       call get_dij_parts(1,qphase,dijh)
       call pawio_print_ij(my_unt,dij2p,lmn2_size,tmp_cplex_dij,lmn_size,-1,idum,0, &
&                   my_prtvol,idum,-1.d0,1,opt_sym=2,asym_ij=dij2p_,mode_paral=my_mode)
       LIBPAW_DEALLOCATE(dijh)
     end if

     !Dij Hat
     if (Paw_ij(iatom)%has_dijhat/=0.and.(idij<=2.or.nspden==4)) then
       if (my_ipert==0) then
         write(msg,'(a)') '   **************** Dij_hat *****************'
       else
         write(msg,'(a)') '   ***** Dij_hat(1) (incl. frozen Dij) ******'
       end if
       call wrtout(my_unt,msg,my_mode)
       call get_dij_parts(cplex_dij,qphase,Paw_ij(iatom)%dijhat)
       call pawio_print_ij(my_unt,dij2p,lmn2_size,tmp_cplex_dij,lmn_size,-1,idum,0,&
&                   my_prtvol,idum,-1.d0,1,opt_sym=2,asym_ij=dij2p_,mode_paral=my_mode)
     end if

     !Dij nuclear dipole
     if (Paw_ij(iatom)%has_dijnd/=0.and.my_ipert<=0) then
       write(msg,'(a)') '   *********** Dij Nuclear Dipole **********'
       call wrtout(my_unt,msg,my_mode)
       call get_dij_parts(cplex_dij,1,Paw_ij(iatom)%dijnd,always_img=.true.)
       call pawio_print_ij(my_unt,dij2p,lmn2_size,tmp_cplex_dij,lmn_size,-1,idum,0,&
&                   my_prtvol,idum,-1.d0,1,opt_sym=2,asym_ij=dij2p_,mode_paral=my_mode)
     end if

     !Dij spin-orbit
     if (Paw_ij(iatom)%has_dijso/=0.and.my_ipert<=0) then
       write(msg,'(a)') '   ************** Dij SpinOrbit ************'
       call wrtout(my_unt,msg,my_mode)
       call get_dij_parts(cplex_dij,qphase,Paw_ij(iatom)%dijso,always_img=.true.)
       call pawio_print_ij(my_unt,dij2p,lmn2_size,tmp_cplex_dij,lmn_size,-1,idum,0,&
&                   my_prtvol,idum,-1.d0,1,opt_sym=2,asym_ij=dij2p_,mode_paral=my_mode)
     end if

     !Dij LDA+U
     if (Paw_ij(iatom)%has_dijU/=0.and.(idij<=2.or.nspden==4).and.my_ipert<=0) then
       write(msg,'(a)') '   ************* Dij_LDA+U (dijpawu) **********'
       call wrtout(my_unt,msg,my_mode)
       call get_dij_parts(cplex_dij,qphase,Paw_ij(iatom)%diju)
       call pawio_print_ij(my_unt,dij2p,lmn2_size,tmp_cplex_dij,lmn_size,-1,idum,0,&
&                   my_prtvol,idum,-1.d0,1,opt_sym=2,asym_ij=dij2p_,mode_paral=my_mode)
     end if

     !Dij XC
     if (Paw_ij(iatom)%has_dijxc/=0.and.(idij<=2.or.nspden==4)) then
       if (my_ipert<=0) then
         write(msg,'(a)') '   ***************** Dij_xc *****************'
       else
         write(msg,'(a)') '   **************** Dij(1)_xc ***************'
       end if
       call wrtout(my_unt,msg,my_mode)
       call get_dij_parts(cplex_dij,qphase,Paw_ij(iatom)%dijxc)
       call pawio_print_ij(my_unt,dij2p,lmn2_size,tmp_cplex_dij,lmn_size,-1,idum,0,&
&                   my_prtvol,idum,-1.d0,1,opt_sym=2,asym_ij=dij2p_,mode_paral=my_mode)
     end if

     !Dij hat XC
     if (Paw_ij(iatom)%has_dijxc_hat/=0.and.(idij<=2.or.nspden==4).and.my_ipert<=0) then
       if (my_ipert<=0) then
         write(msg,'(a)') '   *************** Dijhat_xc ****************'
       else
         write(msg,'(a)') '   ************** Dij(1)hat_xc **************'
       end if
       call wrtout(my_unt,msg,my_mode)
       call get_dij_parts(cplex_dij,1,Paw_ij(iatom)%dijxc_hat)
       call pawio_print_ij(my_unt,dij2p,lmn2_size,tmp_cplex_dij,lmn_size,-1,idum,0,&
&                   my_prtvol,idum,-1.d0,1,opt_sym=2,asym_ij=dij2p_,mode_paral=my_mode)
     end if

     !Dij XC val
     if (Paw_ij(iatom)%has_dijxc_val/=0.and.(idij<=2.or.nspden==4).and.my_ipert<=0) then
       write(msg, '(a)') '   *************** Dij_xc_val ***************'
       call wrtout(my_unt,msg,my_mode)
       call get_dij_parts(cplex_dij,1,Paw_ij(iatom)%dijxc_val)
       call pawio_print_ij(my_unt,dij2p,lmn2_size,tmp_cplex_dij,lmn_size,-1,idum,0,&
&                   my_prtvol,idum,-1.d0,1,opt_sym=2,asym_ij=dij2p_,mode_paral=my_mode)
     end if

     !Dij TOTAL
     if (Paw_ij(iatom)%has_dij/=0) then
       if (my_ipert<=0) then
         write(msg,'(a)') '   **********    TOTAL Dij in Ha   **********'
       else
         write(msg,'(a)') '   **********  TOTAL Dij(1) in Ha  **********'
       end if
       call wrtout(my_unt,msg,my_mode)
       call get_dij_parts(cplex_dij,qphase,Paw_ij(iatom)%dij,always_img=.true.)
       call pawio_print_ij(my_unt,dij2p,lmn2_size,tmp_cplex_dij,lmn_size,-1,idum,0,&
&           my_prtvol,idum,50.d0*dble(3-2*idij),1,opt_sym=2,asym_ij=dij2p_,mode_paral=my_mode)
       if (my_enunit>0) then
         if (my_ipert<=0) then
          write(msg,'(a)') '   **********    TOTAL Dij in eV   **********'
         else
          write(msg,'(a)') '   **********  TOTAL Dij(1) in eV  **********'
         end if
         call pawio_print_ij(my_unt,dij2p,lmn2_size,tmp_cplex_dij,lmn_size,-1,idum,0,&
&           my_prtvol,idum,-1._dp,2,opt_sym=2,asym_ij=dij2p_,mode_paral=my_mode)
       end if
     end if

     if (allocated(dij)) then
       LIBPAW_DEALLOCATE(dij)
     end if
     if (allocated(dijs)) then
       LIBPAW_DEALLOCATE(dijs)
     end if

   end if   !(ABS(my_prtvol)>=1.and.(iatom_tot==1.or.iatom_tot==my_natom.or.my_prtvol<0)

!  =================== Standard output =====================================
   if ((abs(my_prtvol)==0).and.(iatom_tot==1.or.iatom_tot==my_natom)) then

     !Title
     if (idij==1) then
       if (my_ipert<=0) then
         write(msg, '(2a,i6,a)') ch10,' ****** Psp strength Dij in Ha (atom ',iatom_tot,') *****'
       else
         write(msg, '(2a,i6,a)') ch10,' **** Psp strength  Dij(1) in Ha (atom ',iatom_tot,') *****'
       end if
       if (nspden==2.and.nsppol==1) then
         write(msg,'(4a)') trim(msg),') *****',ch10,' (antiferromagnetism case: only one spin component)'
       end if
       call wrtout(my_unt,msg,my_mode)
     end if
     if (paw_ij(iatom)%ndij/=1) then
       write(msg,'(3a)') ' Component ',trim(dspin(idij+2*(nsploop/4))),':'
       call wrtout(my_unt,msg,my_mode)
     end if

     !Dij TOTAL
     if (Paw_ij(iatom)%has_dij/=0) then
       call get_dij_parts(cplex_dij,qphase,Paw_ij(iatom)%dij,always_img=.true.)
       call pawio_print_ij(my_unt,dij2p,lmn2_size,tmp_cplex_dij,lmn_size,-1,idum,0,&
&           my_prtvol,idum,50.d0*dble(3-2*idij),1,opt_sym=2,asym_ij=dij2p_,mode_paral=my_mode)
     end if

   end if

!  =================== End main loops =====================================

   if (allocated(dij)) then
     LIBPAW_DEALLOCATE(dij)
   end if
   if (allocated(dijs)) then
     LIBPAW_DEALLOCATE(dijs)
   end if

  end do !idij
 end do !iat

 call wrtout(my_unt,' ',my_mode)

!Small helper function
 contains

!Real and imaginary parts of phase.
   subroutine get_dij_parts(my_cplex_dij,my_qphase,my_dij,always_img)

     integer,intent(in) :: my_cplex_dij,my_qphase
     logical,intent(in),optional :: always_img
     real(dp),intent(in),target :: my_dij(:,:)
     integer :: my_idij,my_idij_sym,kk
     logical :: always_img_
     always_img_=.false.;if(present(always_img)) always_img_=always_img
     my_idij=min(size(my_dij,2),idij)
     my_idij_sym=min(size(my_dij,2),idij_sym)
     if (my_qphase==1) then
       if ((idij<=nsppol.or.idij==2).and.(.not.always_img_))then
         tmp_cplex_dij=1
         dij2p  => my_dij(1:my_cplex_dij*lmn2_size:my_cplex_dij,my_idij)
         dij2p_ => dij2p
       else
         tmp_cplex_dij=my_cplex_dij
         dij2p  => my_dij(1:my_cplex_dij*lmn2_size:1,my_idij)
         dij2p_ => my_dij(1:my_cplex_dij*lmn2_size:1,my_idij_sym)
       end if
     else
       tmp_cplex_dij=2
       if (my_cplex_dij==1) then
         do kk=1,lmn2_size
           dij(2*kk-1)= my_dij(kk,my_idij)
           dij(2*kk  )= my_dij(kk+lmn2_size,my_idij)
           dijs(2*kk-1)= my_dij(kk,my_idij)
           dijs(2*kk  )=-my_dij(kk+lmn2_size,my_idij)
         end do
       else
         do kk=1,lmn2_size
           dij(2*kk-1)= my_dij(2*kk-1,idij)-my_dij(2*kk  +2*lmn2_size,my_idij)
           dij(2*kk  )= my_dij(2*kk  ,idij)+my_dij(2*kk-1+2*lmn2_size,my_idij)
           dijs(2*kk-1)= my_dij(2*kk-1,idij_sym)+my_dij(2*kk  +2*lmn2_size,my_idij_sym)
           dijs(2*kk  )= my_dij(2*kk  ,idij_sym)-my_dij(2*kk-1+2*lmn2_size,my_idij_sym)
         end do
       end if
       dij2p => dij ; dij2p_ => dijs
     end if
 end subroutine get_dij_parts

end subroutine paw_ij_print
!!***

!----------------------------------------------------------------------

!!****f* m_paw_ij/paw_ij_gather
!! NAME
!! paw_ij_gather
!!
!! FUNCTION
!!   (All)Gather paw_ij datastructures
!!
!! INPUTS
!!  master=master receiving data ; if -1 do a ALLGATHER
!!  comm_atom= communicator
!!  paw_ij_in(:)<type(paw_ij_type)>= input paw_ij datastructures on every process
!!
!! OUTPUT
!!  paw_ij_gathered(:)<type(paw_ij_type)>= output paw_oij datastructure
!!
!! PARENTS
!!      m_paw_ij,outkss,pawprt
!!
!! CHILDREN
!!
!! SOURCE

subroutine paw_ij_gather(paw_ij_in,paw_ij_gathered,master,comm_atom)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: master,comm_atom
!arrays
 type(paw_ij_type),intent(in) :: paw_ij_in(:)
 type(paw_ij_type),intent(inout) :: paw_ij_gathered(:)

!Local variables-------------------------------
!scalars
 integer :: buf_dp_size,buf_dp_size_all,buf_int_size,buf_int_size_all
 integer :: cplxdij_lmn2_size,cplxdijq_lmn2_size,cplxq_lmn2_size
 integer :: iat,ii,ierr,ij,indx_dp,indx_int,lmn2_size
 integer :: me_atom,ndij,nocc,nocc1,nocc2,nocc3,nocc4,npaw_ij_in,npaw_ij_in_sum
 integer :: npaw_ij_out,nproc_atom,nspden,sz1,sz2,sz3,sz4
 logical :: my_atmtab_allocated,paral_atom
 character(len=500) :: msg
!arrays
 integer :: bufsz(2)
 integer,allocatable :: buf_int(:),buf_int_all(:)
 integer,allocatable :: count_dp(:),count_int(:),count_tot(:),displ_dp(:),displ_int(:)
 integer,pointer :: my_atmtab(:)
 real(dp),allocatable :: buf_dp(:),buf_dp_all(:)

! *************************************************************************

!@Paw_ij_type

 npaw_ij_in=size(paw_ij_in);npaw_ij_out=size(paw_ij_gathered)

 nproc_atom=xmpi_comm_size(comm_atom)
 me_atom=xmpi_comm_rank(comm_atom)

!Special case 1 process
 if (nproc_atom==1) then
   if (master==-1.or.me_atom==master) then
     call paw_ij_free(paw_ij_gathered)
     call paw_ij_nullify(paw_ij_gathered)
     do iat=1,npaw_ij_in
       paw_ij_gathered(iat)%cplex_dij  =paw_ij_in(iat)%cplex_dij
       paw_ij_gathered(iat)%qphase     =paw_ij_in(iat)%qphase
       Paw_ij_gathered(iat)%has_dij    =paw_ij_in(iat)%has_dij
       Paw_ij_gathered(iat)%has_dij0   =paw_ij_in(iat)%has_dij0
       Paw_ij_gathered(iat)%has_dijexxc =paw_ij_in(iat)%has_dijexxc
       Paw_ij_gathered(iat)%has_dijfock =paw_ij_in(iat)%has_dijfock
       Paw_ij_gathered(iat)%has_dijfr   =paw_ij_in(iat)%has_dijfr
       Paw_ij_gathered(iat)%has_dijhartree=paw_ij_in(iat)%has_dijhartree
       Paw_ij_gathered(iat)%has_dijhat =paw_ij_in(iat)%has_dijhat
       Paw_ij_gathered(iat)%has_dijnd  =paw_ij_in(iat)%has_dijnd
       Paw_ij_gathered(iat)%has_dijso  =paw_ij_in(iat)%has_dijso
       Paw_ij_gathered(iat)%has_dijU   =paw_ij_in(iat)%has_dijU
       Paw_ij_gathered(iat)%has_dijxc  =paw_ij_in(iat)%has_dijxc
       Paw_ij_gathered(iat)%has_dijxc_hat =paw_ij_in(iat)%has_dijxc_hat
       Paw_ij_gathered(iat)%has_dijxc_val =paw_ij_in(iat)%has_dijxc_val
       Paw_ij_gathered(iat)%has_exexch_pot=paw_ij_in(iat)%has_exexch_pot
       Paw_ij_gathered(iat)%has_pawu_occ  =paw_ij_in(iat)%has_pawu_occ
       paw_ij_gathered(iat)%itypat     =paw_ij_in(iat)%itypat
       paw_ij_gathered(iat)%lmn_size   =paw_ij_in(iat)%lmn_size
       paw_ij_gathered(iat)%lmn2_size  =paw_ij_in(iat)%lmn2_size
       paw_ij_gathered(iat)%ndij       =paw_ij_in(iat)%ndij
       paw_ij_gathered(iat)%nspden     =paw_ij_in(iat)%nspden
       paw_ij_gathered(iat)%nsppol     =paw_ij_in(iat)%nsppol
       lmn2_size=paw_ij_gathered(iat)%lmn2_size
       cplxdij_lmn2_size=paw_ij_gathered(iat)%cplex_dij*lmn2_size
       cplxq_lmn2_size=paw_ij_gathered(iat)%qphase*lmn2_size
       cplxdijq_lmn2_size=cplxdij_lmn2_size*paw_ij_gathered(iat)%qphase
       ndij=paw_ij_gathered(iat)%ndij
       nspden=paw_ij_gathered(iat)%nspden
       if (paw_ij_gathered(iat)%has_dij>=1) then
         LIBPAW_ALLOCATE(paw_ij_gathered(iat)%dij,(cplxdijq_lmn2_size,ndij))
         if (paw_ij_in(iat)%has_dij==2) then
           paw_ij_gathered(iat)%dij=paw_ij_in(iat)%dij
         end if
       end if
       if (paw_ij_gathered(iat)%has_dij0 >=1) then
         LIBPAW_ALLOCATE(paw_ij_gathered(iat)%dij0,(lmn2_size))
         if (paw_ij_in(iat)%has_dij0==2) then
           paw_ij_gathered(iat)%dij0=paw_ij_in(iat)%dij0
         end if
       end if
       if (paw_ij_gathered(iat)%has_dijexxc >=1) then
         LIBPAW_ALLOCATE(paw_ij_gathered(iat)%dijexxc,(cplxdij_lmn2_size,ndij))
         if (paw_ij_in(iat)%has_dijexxc==2) then
           paw_ij_gathered(iat)%dijexxc(:,:)=paw_ij_in(iat)%dijexxc(:,:)
         end if
       end if
       if (paw_ij_gathered(iat)%has_dijfock >=1) then
         LIBPAW_ALLOCATE(paw_ij_gathered(iat)%dijfock,(cplxdij_lmn2_size,ndij))
         if (paw_ij_in(iat)%has_dijfock==2) then
           paw_ij_gathered(iat)%dijfock(:,:)=paw_ij_in(iat)%dijfock(:,:)
         end if
       end if
       if (paw_ij_gathered(iat)%has_dijfr >=1) then
         LIBPAW_ALLOCATE(paw_ij_gathered(iat)%dijfr,(cplxdijq_lmn2_size,ndij))
         if (paw_ij_in(iat)%has_dijfr==2) then
           paw_ij_gathered(iat)%dijfr(:,:)=paw_ij_in(iat)%dijfr(:,:)
         end if
       end if
         if (paw_ij_gathered(iat)%has_dijhartree >=1) then
         LIBPAW_ALLOCATE(paw_ij_gathered(iat)%dijhartree,(cplxq_lmn2_size))
         if (paw_ij_in(iat)%has_dijhartree==2) then
           paw_ij_gathered(iat)%dijhartree(:)=paw_ij_in(iat)%dijhartree(:)
         end if
       end if
       if (paw_ij_gathered(iat)%has_dijhat >=1) then
         LIBPAW_ALLOCATE(paw_ij_gathered(iat)%dijhat,(cplxdijq_lmn2_size,ndij))
         if (paw_ij_in(iat)%has_dijhat==2) then
           paw_ij_gathered(iat)%dijhat(:,:)=paw_ij_in(iat)%dijhat(:,:)
         end if
       end if
       if (paw_ij_gathered(iat)%has_dijnd >=1) then
         LIBPAW_ALLOCATE(paw_ij_gathered(iat)%dijnd,(cplxdij_lmn2_size,ndij))
         if (paw_ij_in(iat)%has_dijnd==2) then
           paw_ij_gathered(iat)%dijnd(:,:)=paw_ij_in(iat)%dijnd(:,:)
         end if
       end if
       if (paw_ij_gathered(iat)%has_dijU >=1) then
         LIBPAW_ALLOCATE(paw_ij_gathered(iat)%dijU,(cplxdijq_lmn2_size,ndij))
         if (paw_ij_in(iat)%has_dijU==2) then
         paw_ij_gathered(iat)%dijU(:,:)=paw_ij_in(iat)%dijU(:,:)
         end if
       end if
       if (paw_ij_gathered(iat)%has_dijso >=1) then
         LIBPAW_ALLOCATE(paw_ij_gathered(iat)%dijso,(cplxdijq_lmn2_size,ndij))
         if (paw_ij_in(iat)%has_dijso==2) then
           paw_ij_gathered(iat)%dijso(:,:)=paw_ij_in(iat)%dijso(:,:)
         end if
       end if
       if (paw_ij_gathered(iat)%has_dijxc >=1) then
         LIBPAW_ALLOCATE(paw_ij_gathered(iat)%dijxc,(cplxdijq_lmn2_size,ndij))
         if (paw_ij_in(iat)%has_dijxc==2) then
           paw_ij_gathered(iat)%dijxc(:,:)=paw_ij_in(iat)%dijxc(:,:)
         end if
       end if
       if (paw_ij_gathered(iat)%has_dijxc_hat >=1) then
         LIBPAW_ALLOCATE(paw_ij_gathered(iat)%dijxc_hat,(cplxdij_lmn2_size,ndij))
         if (paw_ij_in(iat)%has_dijxc_hat==2) then
           paw_ij_gathered(iat)%dijxc_hat(:,:)=paw_ij_in(iat)%dijxc_hat(:,:)
         end if
       end if
       if (paw_ij_gathered(iat)%has_dijxc_val >=1) then
         LIBPAW_ALLOCATE(paw_ij_gathered(iat)%dijxc_val,(cplxdij_lmn2_size,ndij))
         if (paw_ij_in(iat)%has_dijxc_val==2) then
           paw_ij_gathered(iat)%dijxc_val(:,:)=paw_ij_in(iat)%dijxc_val(:,:)
         end if
       end if
       if (paw_ij_gathered(iat)%has_pawu_occ >=1) then
         LIBPAW_ALLOCATE(paw_ij_gathered(iat)%nocctot,(ndij))
         if (paw_ij_gathered(iat)%has_pawu_occ==2) then
           paw_ij_gathered(iat)%nocctot(:)=paw_ij_in(iat)%nocctot(:)
           if (allocated(paw_ij_in(iat)%noccmmp)) then
             sz1=size(paw_ij_in(iat)%noccmmp,1);sz2=size(paw_ij_in(iat)%noccmmp,2)
             sz3=size(paw_ij_in(iat)%noccmmp,3);sz4=size(paw_ij_in(iat)%noccmmp,4)
             LIBPAW_ALLOCATE(paw_ij_gathered(iat)%noccmmp,(sz1,sz2,sz3,sz4))
             paw_ij_gathered(iat)%noccmmp(:,:,:,:)=paw_ij_in(iat)%noccmmp(:,:,:,:)
           end if
         end if
       end if
       if (paw_ij_in(iat)%has_exexch_pot >=1) then
         sz1=size(paw_ij_in(iat)%vpawx,1);sz2=size(paw_ij_in(iat)%vpawx,2)
         sz3=size(paw_ij_in(iat)%vpawx,3)
         LIBPAW_ALLOCATE(paw_ij_gathered(iat)%vpawx,(sz1,sz2,sz3))
         if (paw_ij_in(iat)%has_exexch_pot==2) then
           paw_ij_gathered(iat)%vpawx(:,:,:)=paw_ij_in(iat)%vpawx(:,:,:)
         end if
       end if
     end do !  iat
   end if
   return
 end if !nproc_atom =1

!Test on sizes
 npaw_ij_in_sum=npaw_ij_in

 call xmpi_sum(npaw_ij_in_sum,comm_atom,ierr)
  if (master==-1) then
   if (npaw_ij_out/=npaw_ij_in_sum) then
     msg='Wrong sizes sum[npaw_ij_ij]/=npaw_ij_out !'
     MSG_BUG(msg)
   end if
 else
   if (me_atom==master.and.npaw_ij_out/=npaw_ij_in_sum) then
     msg='(2) paw_ij_gathered wrongly allocated !'
     MSG_BUG(msg)
   end if
 end if

!Retrieve table of atoms
 paral_atom=.true.;nullify(my_atmtab)
 call get_my_atmtab(comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,npaw_ij_in_sum)

!Compute sizes of buffers
 buf_int_size=0;buf_dp_size=0
 do ij=1,npaw_ij_in
   lmn2_size=paw_ij_in(ij)%lmn2_size
   cplxdij_lmn2_size=paw_ij_in(ij)%cplex_dij*lmn2_size
   cplxq_lmn2_size=paw_ij_in(ij)%qphase*lmn2_size
   cplxdijq_lmn2_size=cplxdij_lmn2_size*paw_ij_in(ij)%qphase
   ndij=paw_ij_in(ij)%ndij
   buf_int_size=buf_int_size+24
   if (paw_ij_in(ij)%has_dij==2) then
     buf_dp_size=buf_dp_size +cplxdijq_lmn2_size*ndij
   end if
   if (paw_ij_in(ij)%has_dij0==2) then
     buf_dp_size=buf_dp_size +lmn2_size
   end if
   if (paw_ij_in(ij)%has_dijexxc==2) then
     buf_dp_size=buf_dp_size +cplxdij_lmn2_size*ndij
   end if
   if (paw_ij_in(ij)%has_dijfock==2) then
     buf_dp_size=buf_dp_size +cplxdij_lmn2_size*ndij
   end if
   if (paw_ij_in(ij)%has_dijfr==2) then
     buf_dp_size=buf_dp_size +cplxdijq_lmn2_size*ndij
   end if
   if (paw_ij_in(ij)%has_dijhartree==2) then
     buf_dp_size=buf_dp_size +cplxq_lmn2_size
   end if
   if (paw_ij_in(ij)%has_dijhat==2) then
     buf_dp_size=buf_dp_size +cplxdijq_lmn2_size*ndij
   end if
   if (paw_ij_in(ij)%has_dijnd==2) then
     buf_dp_size=buf_dp_size +cplxdij_lmn2_size*ndij
   end if
   if (paw_ij_in(ij)%has_dijso==2) then
     buf_dp_size=buf_dp_size +cplxdijq_lmn2_size*ndij
   end if
   if (paw_ij_in(ij)%has_dijU==2) then
     buf_dp_size=buf_dp_size +cplxdijq_lmn2_size*ndij
   end if
   if (paw_ij_in(ij)%has_dijxc==2) then
     buf_dp_size=buf_dp_size +cplxdijq_lmn2_size*ndij
   end if
   if (paw_ij_in(ij)%has_dijxc_hat==2) then
     buf_dp_size=buf_dp_size +cplxdij_lmn2_size*ndij
   end if
   if (paw_ij_in(ij)%has_dijxc_val==2) then
     buf_dp_size=buf_dp_size +cplxdij_lmn2_size*ndij
   end if
   if (paw_ij_in(ij)%has_pawu_occ>=1) then
     buf_int_size=buf_int_size+4
     buf_dp_size=buf_dp_size &
&        +size(paw_ij_in(ij)%nocctot) &
&        +size(paw_ij_in(ij)%noccmmp)
   end if
   if (paw_ij_in(ij)%has_exexch_pot>=1) then
     buf_int_size=buf_int_size+3
     buf_dp_size=buf_dp_size +size(paw_ij_in(ij)%vpawx)
   end if
 end do

!Fill input buffers
 LIBPAW_ALLOCATE(buf_int,(buf_int_size))
 LIBPAW_ALLOCATE(buf_dp ,(buf_dp_size))
 indx_int=1;indx_dp=1
 do ij=1,npaw_ij_in
   nspden=paw_ij_in(ij)%nspden
   buf_int(indx_int)=my_atmtab(ij) ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij_in(ij)%cplex_dij ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij_in(ij)%qphase ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij_in(ij)%itypat ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij_in(ij)%nspden ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij_in(ij)%nsppol ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij_in(ij)%lmn_size ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij_in(ij)%lmn2_size ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij_in(ij)%ndij ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij_in(ij)%has_dij ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij_in(ij)%has_dij0 ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij_in(ij)%has_dijexxc ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij_in(ij)%has_dijfock ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij_in(ij)%has_dijfr ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij_in(ij)%has_dijhartree ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij_in(ij)%has_dijhat ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij_in(ij)%has_dijnd ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij_in(ij)%has_dijso ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij_in(ij)%has_dijU ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij_in(ij)%has_dijxc ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij_in(ij)%has_dijxc_hat ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij_in(ij)%has_dijxc_val ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij_in(ij)%has_exexch_pot ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij_in(ij)%has_pawu_occ ;indx_int=indx_int+1
   lmn2_size=paw_ij_in(ij)%lmn2_size
   cplxdij_lmn2_size=paw_ij_in(ij)%cplex_dij*lmn2_size
   cplxq_lmn2_size=paw_ij_in(ij)%qphase*lmn2_size
   cplxdijq_lmn2_size=cplxdij_lmn2_size*paw_ij_in(ij)%qphase
   ndij=paw_ij_in(ij)%ndij
   if (paw_ij_in(ij)%has_dij==2) then
     ii=cplxdijq_lmn2_size*ndij
     buf_dp(indx_dp:indx_dp+ii-1)=reshape(paw_ij_in(ij)%dij,(/ii/))
     indx_dp=indx_dp+ii
   end if
   if (paw_ij_in(ij)%has_dij0==2) then
     ii=lmn2_size
     buf_dp(indx_dp:indx_dp+ii-1)=paw_ij_in(ij)%dij0(1:ii)
     indx_dp=indx_dp+ii
   end if
   if (paw_ij_in(ij)%has_dijexxc==2) then
     ii=cplxdij_lmn2_size*ndij
     buf_dp(indx_dp:indx_dp+ii-1)=reshape(paw_ij_in(ij)%dijexxc,(/ii/))
     indx_dp=indx_dp+ii
   end if
   if (paw_ij_in(ij)%has_dijfock==2) then
     ii=cplxdij_lmn2_size*ndij
     buf_dp(indx_dp:indx_dp+ii-1)=reshape(paw_ij_in(ij)%dijfock,(/ii/))
     indx_dp=indx_dp+ii
   end if
   if (paw_ij_in(ij)%has_dijfr==2) then
     ii=cplxdijq_lmn2_size*ndij
     buf_dp(indx_dp:indx_dp+ii-1)=reshape(paw_ij_in(ij)%dijfr,(/ii/))
     indx_dp=indx_dp+ii
   end if
   if (paw_ij_in(ij)%has_dijhartree==2) then
     ii=cplxq_lmn2_size
     buf_dp(indx_dp:indx_dp+ii-1)=paw_ij_in(ij)%dijhartree(1:ii)
     indx_dp=indx_dp+ii
   end if
   if (paw_ij_in(ij)%has_dijhat==2) then
     ii=cplxdijq_lmn2_size*ndij
     buf_dp(indx_dp:indx_dp+ii-1)=reshape(paw_ij_in(ij)%dijhat,(/ii/))
     indx_dp=indx_dp+ii
   end if
   if (paw_ij_in(ij)%has_dijnd==2) then
     ii=cplxdij_lmn2_size*ndij
     buf_dp(indx_dp:indx_dp+ii-1)=reshape(paw_ij_in(ij)%dijnd,(/ii/))
     indx_dp=indx_dp+ii
   end if
   if (paw_ij_in(ij)%has_dijso==2) then
     ii=cplxdijq_lmn2_size*ndij
     buf_dp(indx_dp:indx_dp+ii-1)=reshape(paw_ij_in(ij)%dijso,(/ii/))
     indx_dp=indx_dp+ii
   end if
   if (paw_ij_in(ij)%has_dijU==2) then
     ii=cplxdijq_lmn2_size*ndij
     buf_dp(indx_dp:indx_dp+ii-1)=reshape(paw_ij_in(ij)%dijU,(/ii/))
     indx_dp=indx_dp+ii
   end if
   if (paw_ij_in(ij)%has_dijxc==2) then
     ii=cplxdijq_lmn2_size*ndij
     buf_dp(indx_dp:indx_dp+ii-1)=reshape(paw_ij_in(ij)%dijxc,(/ii/))
     indx_dp=indx_dp+ii
   end if
   if (paw_ij_in(ij)%has_dijxc_hat==2) then
     ii=cplxdij_lmn2_size*ndij
     buf_dp(indx_dp:indx_dp+ii-1)=reshape(paw_ij_in(ij)%dijxc_hat,(/ii/))
     indx_dp=indx_dp+ii
   end if
   if (paw_ij_in(ij)%has_dijxc_val==2) then
     ii=cplxdij_lmn2_size*ndij
     buf_dp(indx_dp:indx_dp+ii-1)=reshape(paw_ij_in(ij)%dijxc_val,(/ii/))
     indx_dp=indx_dp+ii
   end if
   if (paw_ij_in(ij)%has_pawu_occ>=1)then
     buf_int(indx_int)=size(paw_ij_in(ij)%noccmmp,1) ;indx_int=indx_int+1
     buf_int(indx_int)=size(paw_ij_in(ij)%noccmmp,2) ;indx_int=indx_int+1
     buf_int(indx_int)=size(paw_ij_in(ij)%noccmmp,3) ;indx_int=indx_int+1
     buf_int(indx_int)=size(paw_ij_in(ij)%noccmmp,4) ;indx_int=indx_int+1
     nocc=paw_ij_in(ij)%ndij
     buf_dp(indx_dp:indx_dp+nocc-1)=paw_ij_in(ij)%nocctot(1:nocc)
     indx_dp=indx_dp+nocc
     nocc=size(paw_ij_in(ij)%noccmmp)
     buf_dp(indx_dp:indx_dp+nocc-1)=reshape(paw_ij_in(ij)%noccmmp,(/nocc/))
     indx_dp=indx_dp+nocc
   end if
   if (paw_ij_in(ij)%has_exexch_pot>=1) then
     sz1=size(paw_ij_in(ij)%vpawx,1);sz2=size(paw_ij_in(ij)%vpawx,2)
     sz3=size(paw_ij_in(ij)%vpawx,3)
     buf_int(indx_int)=sz1; indx_int=indx_int+1
     buf_int(indx_int)=sz2; indx_int=indx_int+1
     buf_int(indx_int)=sz3; indx_int=indx_int+1
     ii=sz1*sz2*sz3
     buf_dp(indx_dp:indx_dp+ii-1)=reshape(paw_ij_in(ij)%vpawx,(/(ii)/))
     indx_dp=indx_dp+ii
   end if
 end do

!Check
 indx_int=indx_int-1;indx_dp=indx_dp-1
 if ((indx_int/=buf_int_size).or.(indx_dp/=buf_dp_size)) then
   write(msg,*) 'Wrong buffer sizes: buf_int_size=',buf_int_size, &
     ' indx_int_size=',indx_int,' buf_dp_size=',buf_dp_size , ' indx_dp=', indx_dp
   MSG_BUG(msg)
 end if

!Communicate (1 gather for integers, 1 gather for reals)
 LIBPAW_ALLOCATE(count_int,(nproc_atom))
 LIBPAW_ALLOCATE(displ_int,(nproc_atom))
 LIBPAW_ALLOCATE(count_dp ,(nproc_atom))
 LIBPAW_ALLOCATE(displ_dp ,(nproc_atom))
 LIBPAW_ALLOCATE(count_tot,(2*nproc_atom))
 bufsz(1)=buf_int_size; bufsz(2)=buf_dp_size
 call xmpi_allgather(bufsz,2,count_tot,comm_atom,ierr)
 do ij=1,nproc_atom
   count_int(ij)=count_tot(2*ij-1)
   count_dp (ij)=count_tot(2*ij)
 end do
 displ_int(1)=0;displ_dp(1)=0
 do ij=2,nproc_atom
   displ_int(ij)=displ_int(ij-1)+count_int(ij-1)
   displ_dp (ij)=displ_dp (ij-1)+count_dp (ij-1)
 end do
 buf_int_size_all=sum(count_int)
 buf_dp_size_all =sum(count_dp)
 LIBPAW_DEALLOCATE(count_tot)
 if (master==-1.or.me_atom==master) then
   LIBPAW_ALLOCATE(buf_int_all,(buf_int_size_all))
   LIBPAW_ALLOCATE(buf_dp_all ,(buf_dp_size_all))
 else
   LIBPAW_ALLOCATE(buf_int_all,(0))
   LIBPAW_ALLOCATE(buf_dp_all ,(0))
 end if
 if (master==-1) then
   call xmpi_allgatherv(buf_int,buf_int_size,buf_int_all,count_int,displ_int,comm_atom,ierr)
   call xmpi_allgatherv(buf_dp ,buf_dp_size ,buf_dp_all ,count_dp ,displ_dp ,comm_atom,ierr)
 else
   call xmpi_gatherv(buf_int,buf_int_size,buf_int_all,count_int,displ_int,master,comm_atom,ierr)
   call xmpi_gatherv(buf_dp ,buf_dp_size ,buf_dp_all ,count_dp ,displ_dp ,master,comm_atom,ierr)
 end if
 LIBPAW_DEALLOCATE(count_int)
 LIBPAW_DEALLOCATE(displ_int)
 LIBPAW_DEALLOCATE(count_dp)
 LIBPAW_DEALLOCATE(displ_dp)

!Retrieve data from output buffer
 if (master==-1.or.me_atom==master) then
   indx_int=1;indx_dp=1
   call paw_ij_free(paw_ij_gathered)
   call paw_ij_nullify(paw_ij_gathered)
   do ij=1,npaw_ij_out
     iat=buf_int_all(indx_int) ;indx_int=indx_int+1
     paw_ij_gathered(iat)%cplex_dij=buf_int_all(indx_int) ;indx_int=indx_int+1
     paw_ij_gathered(iat)%qphase=buf_int_all(indx_int) ;indx_int=indx_int+1
     paw_ij_gathered(iat)%itypat=buf_int_all(indx_int) ;indx_int=indx_int+1
     paw_ij_gathered(iat)%nspden=buf_int_all(indx_int) ;indx_int=indx_int+1
     paw_ij_gathered(iat)%nsppol=buf_int_all(indx_int) ;indx_int=indx_int+1
     paw_ij_gathered(iat)%lmn_size=buf_int_all(indx_int) ;indx_int=indx_int+1
     paw_ij_gathered(iat)%lmn2_size=buf_int_all(indx_int) ;indx_int=indx_int+1
     paw_ij_gathered(iat)%ndij=buf_int_all(indx_int) ;indx_int=indx_int+1
     paw_ij_gathered(iat)%has_dij=buf_int_all(indx_int) ;indx_int=indx_int+1
     paw_ij_gathered(iat)%has_dij0=buf_int_all(indx_int) ;indx_int=indx_int+1
     paw_ij_gathered(iat)%has_dijexxc=buf_int_all(indx_int) ;indx_int=indx_int+1
     paw_ij_gathered(iat)%has_dijfock=buf_int_all(indx_int) ;indx_int=indx_int+1
     paw_ij_gathered(iat)%has_dijfr=buf_int_all(indx_int) ;indx_int=indx_int+1
     paw_ij_gathered(iat)%has_dijhartree=buf_int_all(indx_int) ;indx_int=indx_int+1
     paw_ij_gathered(iat)%has_dijhat=buf_int_all(indx_int) ;indx_int=indx_int+1
     paw_ij_gathered(iat)%has_dijnd=buf_int_all(indx_int) ;indx_int=indx_int+1
     paw_ij_gathered(iat)%has_dijso=buf_int_all(indx_int) ;indx_int=indx_int+1
     paw_ij_gathered(iat)%has_dijU=buf_int_all(indx_int) ;indx_int=indx_int+1
     paw_ij_gathered(iat)%has_dijxc=buf_int_all(indx_int) ;indx_int=indx_int+1
     paw_ij_gathered(iat)%has_dijxc_hat=buf_int_all(indx_int) ;indx_int=indx_int+1
     paw_ij_gathered(iat)%has_dijxc_val=buf_int_all(indx_int) ;indx_int=indx_int+1
     paw_ij_gathered(iat)%has_exexch_pot=buf_int_all(indx_int) ;indx_int=indx_int+1
     paw_ij_gathered(iat)%has_pawu_occ=buf_int_all(indx_int) ;indx_int=indx_int+1
     if (paw_ij_gathered(iat)%has_pawu_occ>=1) then
       nocc1=buf_int_all(indx_int) ;indx_int=indx_int+1
       nocc2=buf_int_all(indx_int) ;indx_int=indx_int+1
       nocc3=buf_int_all(indx_int) ;indx_int=indx_int+1
       nocc4=buf_int_all(indx_int) ;indx_int=indx_int+1
     else
       nocc1=0;nocc2=0;nocc3=0;nocc4=0
     end if
     lmn2_size=paw_ij_gathered(iat)%lmn2_size
     cplxdij_lmn2_size=paw_ij_gathered(iat)%cplex_dij*lmn2_size
     cplxq_lmn2_size=paw_ij_gathered(iat)%qphase*lmn2_size
     cplxdijq_lmn2_size=cplxdij_lmn2_size*paw_ij_gathered(iat)%qphase
     ndij=paw_ij_gathered(iat)%ndij

     if (paw_ij_gathered(iat)%has_dij>=1) then
       ii=cplxdijq_lmn2_size
       LIBPAW_ALLOCATE(paw_ij_gathered(iat)%dij,(ii,ndij))
       if (paw_ij_gathered(iat)%has_dij==2) then
         paw_ij_gathered(iat)%dij(:,:)= &
 &         reshape(buf_dp_all(indx_dp:indx_dp+ii*ndij-1),(/ii,ndij/))
         indx_dp=indx_dp+ii*ndij
       end if
     end if
     if (paw_ij_gathered(iat)%has_dij0 >=1) then
       ii=lmn2_size
       LIBPAW_ALLOCATE(paw_ij_gathered(iat)%dij0,(ii))
       if (paw_ij_gathered(iat)%has_dij0==2) then
         paw_ij_gathered(iat)%dij0(:)=buf_dp_all(indx_dp:indx_dp+ii-1)
         indx_dp=indx_dp+ii
       end if
     end if
     if (paw_ij_gathered(iat)%has_dijexxc >=1) then
       ii=cplxdij_lmn2_size
       LIBPAW_ALLOCATE(paw_ij_gathered(iat)%dijexxc,(ii,ndij))
       if (paw_ij_gathered(iat)%has_dijexxc==2) then
         paw_ij_gathered(iat)%dijexxc(:,:)= &
&          reshape(buf_dp_all(indx_dp:indx_dp+ii*ndij-1),(/ii,ndij/))
         indx_dp=indx_dp+ii*ndij
       end if
     end if
     if (paw_ij_gathered(iat)%has_dijfock >=1) then
       ii=cplxdij_lmn2_size
       LIBPAW_ALLOCATE(paw_ij_gathered(iat)%dijfock,(ii,ndij))
       if (paw_ij_gathered(iat)%has_dijfock==2) then
         paw_ij_gathered(iat)%dijfock(:,:)= &
&          reshape(buf_dp_all(indx_dp:indx_dp+ii*ndij-1),(/ii,ndij/))
         indx_dp=indx_dp+ii*ndij
       end if
     end if
     if (paw_ij_gathered(iat)%has_dijfr >=1) then
       ii=cplxdijq_lmn2_size
       LIBPAW_ALLOCATE(paw_ij_gathered(iat)%dijfr,(ii,ndij))
       if (paw_ij_gathered(iat)%has_dijfr==2) then
         paw_ij_gathered(iat)%dijfr(:,:)= &
&          reshape(buf_dp_all(indx_dp:indx_dp+ii*ndij-1),(/ii,ndij/))
         indx_dp=indx_dp+ii*ndij
       end if
     end if
     if (paw_ij_gathered(iat)%has_dijhartree >=1) then
       ii=cplxq_lmn2_size
       LIBPAW_ALLOCATE(paw_ij_gathered(iat)%dijhartree,(ii))
       if (paw_ij_gathered(iat)%has_dijhartree==2) then
         paw_ij_gathered(iat)%dijhartree(:)=buf_dp_all(indx_dp:indx_dp+ii-1)
         indx_dp=indx_dp+ii
       end if
     end if
     if (paw_ij_gathered(iat)%has_dijhat >=1) then
       ii=cplxdijq_lmn2_size
       LIBPAW_ALLOCATE(paw_ij_gathered(iat)%dijhat,(ii,ndij))
       if (paw_ij_gathered(iat)%has_dijhat==2) then
         paw_ij_gathered(iat)%dijhat(:,:)= &
&          reshape(buf_dp_all(indx_dp:indx_dp+ii*ndij-1),(/ii,ndij/))
         indx_dp=indx_dp+ii*ndij
       end if
     endif
     if (paw_ij_gathered(iat)%has_dijnd >=1) then
       ii=cplxdij_lmn2_size
       LIBPAW_ALLOCATE(paw_ij_gathered(iat)%dijnd,(ii,ndij))
       if (paw_ij_gathered(iat)%has_dijnd==2) then
         paw_ij_gathered(iat)%dijnd(:,:)= &
&          reshape(buf_dp_all(indx_dp:indx_dp+ii*ndij-1),(/ii,ndij/))
         indx_dp=indx_dp+ii*ndij
       end if
     end if
     if (paw_ij_gathered(iat)%has_dijso >=1) then
       ii=cplxdijq_lmn2_size
       LIBPAW_ALLOCATE(paw_ij_gathered(iat)%dijso,(ii,ndij))
       if (paw_ij_gathered(iat)%has_dijso==2) then
         paw_ij_gathered(iat)%dijso(:,:)= &
&          reshape(buf_dp_all(indx_dp:indx_dp+ii*ndij-1),(/ii,ndij/))
         indx_dp=indx_dp+ii*ndij
       end if
     end if
     if (paw_ij_gathered(iat)%has_dijU >=1) then
       ii=cplxdijq_lmn2_size
       LIBPAW_ALLOCATE(paw_ij_gathered(iat)%dijU,(ii,ndij))
       if (paw_ij_gathered(iat)%has_dijU==2) then
         paw_ij_gathered(iat)%dijU(:,:)= &
&          reshape(buf_dp_all(indx_dp:indx_dp+ii*ndij-1),(/ii,ndij/))
         indx_dp=indx_dp+ii*ndij
       end if
     end if
     if (paw_ij_gathered(iat)%has_dijxc >=1) then
       ii=cplxdijq_lmn2_size
       LIBPAW_ALLOCATE(paw_ij_gathered(iat)%dijxc,(ii,ndij))
       if (paw_ij_gathered(iat)%has_dijxc==2) then
         paw_ij_gathered(iat)%dijxc(:,:)= &
&          reshape(buf_dp_all(indx_dp:indx_dp+ii*ndij-1),(/ii,ndij/))
         indx_dp=indx_dp+ii*ndij
       end if
     end if
     if (paw_ij_gathered(iat)%has_dijxc_hat >=1) then
       ii=cplxdij_lmn2_size
       LIBPAW_ALLOCATE(paw_ij_gathered(iat)%dijxc_hat,(ii,ndij))
       if (paw_ij_gathered(iat)%has_dijxc_hat==2) then
         paw_ij_gathered(iat)%dijxc_hat(:,:)= &
&          reshape(buf_dp_all(indx_dp:indx_dp+ii*ndij-1),(/ii,ndij/))
         indx_dp=indx_dp+ii*ndij
       end if
     end if
     if (paw_ij_gathered(iat)%has_dijxc_val >=1) then
       ii=cplxdij_lmn2_size
       LIBPAW_ALLOCATE(paw_ij_gathered(iat)%dijxc_val,(ii,ndij))
       if (paw_ij_gathered(iat)%has_dijxc_val==2) then
         paw_ij_gathered(iat)%dijxc_val(:,:)= &
&          reshape(buf_dp_all(indx_dp:indx_dp+ii*ndij-1),(/ii,ndij/))
         indx_dp=indx_dp+ii*ndij
       end if
     end if
     if (paw_ij_gathered(iat)%has_pawu_occ >=1) then
       nocc=paw_ij_gathered(iat)%ndij
       LIBPAW_ALLOCATE(paw_ij_gathered(iat)%nocctot,(nocc))
       paw_ij_gathered(iat)%nocctot(1:nocc)=buf_dp_all(indx_dp:indx_dp+nocc-1)
       indx_dp=indx_dp+nocc
       nocc=nocc1*nocc2*nocc3*nocc4
       LIBPAW_ALLOCATE(paw_ij_gathered(iat)%noccmmp,(nocc1,nocc2,nocc3,nocc4))
       paw_ij_gathered(iat)%noccmmp(1:nocc1,1:nocc2,1:nocc3,1:nocc4)= &
&        reshape(buf_dp_all(indx_dp:indx_dp+nocc-1),(/nocc1,nocc2,nocc3,nocc4/))
       indx_dp=indx_dp+nocc
     end if
     if (paw_ij_gathered(iat)%has_exexch_pot >=1) then
       sz1=buf_int_all(indx_int);indx_int=indx_int+1
       sz2=buf_int_all(indx_int);indx_int=indx_int+1
       sz3=buf_int_all(indx_int);indx_int=indx_int+1
       LIBPAW_ALLOCATE(paw_ij_gathered(iat)%vpawx,(sz1,sz2,sz3))
       if (paw_ij_gathered(iat)%has_exexch_pot == 2) then
           paw_ij_gathered(iat)%vpawx(:,:,:)=&
&          reshape(buf_dp_all(indx_dp:indx_dp+sz1*sz2*sz3-1),(/sz1,sz2,sz3/))
         indx_dp=indx_dp+sz1*sz2*sz3
       end if
     end if
   end do
   indx_int=indx_int-1;indx_dp=indx_dp-1
   if ((indx_int/=buf_int_size_all).or.(indx_dp/=buf_dp_size_all)) then
     write(msg,*) 'Wrong buffer sizes: buf_int_size_all=',buf_int_size_all, &
&      ' indx_int=',indx_int, ' buf_dp_size_all=',buf_dp_size_all,' indx_dp=',indx_dp
     MSG_BUG(msg)
   end if
 end if

!Free memory
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)
 LIBPAW_DEALLOCATE(buf_int)
 LIBPAW_DEALLOCATE(buf_dp)
 LIBPAW_DEALLOCATE(buf_int_all)
 LIBPAW_DEALLOCATE(buf_dp_all)

end subroutine paw_ij_gather
!!***

!----------------------------------------------------------------------

!!****f* m_paw_ij/paw_ij_redistribute
!! NAME
!! paw_ij_redistribute
!!
!! FUNCTION
!!   Redistribute an array of paw_ij datastructures
!!   Input paw_ij is given on a MPI communicator
!!   Output paw_ij is redistributed on another MPI communicator
!!
!! INPUTS
!!  mpi_comm_in= input MPI (atom) communicator
!!  mpi_comm_out= output MPI (atom) communicator
!!  mpi_atmtab_in= --optional-- indexes of the input paw_ij treated by current proc
!!                 if not present, will be calculated in the present routine
!!  mpi_atmtab_out= --optional-- indexes of the output paw_ij treated by current proc
!!                  if not present, will be calculated in the present routine
!!  natom= --optional-- total number of atoms
!!  ----- Optional arguments used only for asynchronous communications -----
!!    RecvAtomProc(:)= rank of processor from which I expect atom (in mpi_comm_in)
!!    RecvAtomList(:)= indexes of atoms to be received by me
!!      RecvAtomList(irecv) are the atoms I expect from RecvAtomProc(irecv)
!!    SendAtomProc(:)= ranks of process destination of atom (in mpi_comm_in)
!!    SendAtomList(:)= indexes of atoms to be sent by me
!!      SendAtomList(isend) are the atoms sent to SendAtomProc(isend)
!!
!! OUTPUT
!!  [paw_ij_out(:)]<type(paw_ij_type)>= --optional--
!!                    if present, the redistributed datastructure does not replace
!!                    the input one but is delivered in paw_ij_out
!!                    if not present, input and output datastructure are the same.
!!
!! SIDE EFFECTS
!!  paw_ij(:)<type(paw_ij_type)>= input (and eventually output) paw_ij datastructures
!!
!! PARENTS
!!      m_paral_pert
!!
!! CHILDREN
!!
!! SOURCE

subroutine paw_ij_redistribute(paw_ij,mpi_comm_in,mpi_comm_out,&
&                 natom,mpi_atmtab_in,mpi_atmtab_out,paw_ij_out,&
&                 SendAtomProc,SendAtomList,RecvAtomProc,RecvAtomList)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpi_comm_in,mpi_comm_out
 integer,optional,intent(in) :: natom
!arrays
 integer,intent(in),optional,target :: mpi_atmtab_in(:),mpi_atmtab_out(:)
 type(paw_ij_type),allocatable,intent(inout) :: paw_ij(:)
 type(paw_ij_type),pointer,optional :: paw_ij_out(:)      !vz_i
 integer,intent(in),optional :: SendAtomProc(:),SendAtomList(:),RecvAtomProc(:),RecvAtomList(:)

!Local variables-------------------------------
!scalars
 integer :: algo_option,i1,iatom,iat_in,iat_out,ierr,iircv,iisend,imsg,imsg_current
 integer :: imsg1,iproc_rcv,iproc_send,ireq,me_exch,mpi_comm_exch,my_natom_in,my_natom_out,my_tag,natom_tot
 integer :: nb_dp,nb_int,nb_msg,nbmsg_incoming,nbrecv,nbrecvmsg,nbsendreq,nbsent,nbsend,next,npaw_ij_sent
 integer :: nproc_in,nproc_out
 logical :: flag,in_place,message_yet_prepared,my_atmtab_in_allocated,my_atmtab_out_allocated,paral_atom
!arrays
 integer :: buf_size(3),request1(3)
 integer,pointer :: my_atmtab_in(:),my_atmtab_out(:)
 integer,allocatable :: atmtab_send(:),atm_indx_in(:),atm_indx_out(:),From(:),buf_int1(:),request(:)
 integer,allocatable,target:: buf_int(:)
 integer,pointer :: buf_ints(:)
 logical,allocatable :: msg_pick(:)
 real(dp),allocatable :: buf_dp1(:)
 real(dp),allocatable,target :: buf_dp(:)
 real(dp),pointer :: buf_dps(:)
 type(coeffi1_type),target,allocatable :: tab_buf_int(:),tab_buf_atom(:)
 type(coeff1_type),target,allocatable :: tab_buf_dp(:)
 type(paw_ij_type),allocatable :: paw_ij_all(:)
 type(paw_ij_type),pointer :: paw_ij_out1(:)

! *************************************************************************

!@Paw_ij_type

 in_place=(.not.present(paw_ij_out))
 my_natom_in=size(paw_ij)

!If not "in_place", destroy the output datastructure
 if (.not.in_place) then
   if (associated(paw_ij_out)) then
     call paw_ij_free(paw_ij_out)
     LIBPAW_DATATYPE_DEALLOCATE(paw_ij_out)
   end if
 end if

!Special sequential case
 if (mpi_comm_in==xmpi_comm_self.and.mpi_comm_out==xmpi_comm_self) then
   if ((.not.in_place).and.(my_natom_in>0)) then
     LIBPAW_DATATYPE_ALLOCATE(paw_ij_out,(my_natom_in))
     call paw_ij_nullify(paw_ij_out)
     call paw_ij_copy(paw_ij,paw_ij_out)
   end if
   return
 end if

!Get total natom
 if (present(natom)) then
   natom_tot=natom
 else
   natom_tot=my_natom_in
   call xmpi_sum(natom_tot,mpi_comm_in,ierr)
 end if

!Select input distribution
 if (present(mpi_atmtab_in)) then
   my_atmtab_in => mpi_atmtab_in
   my_atmtab_in_allocated=.false.
 else
   call get_my_atmtab(mpi_comm_in,my_atmtab_in,my_atmtab_in_allocated,&
&                     paral_atom,natom_tot,my_natom_in)
 end if

!Select output distribution
 if (present(mpi_atmtab_out)) then
   my_natom_out=size(mpi_atmtab_out)
   my_atmtab_out => mpi_atmtab_out
   my_atmtab_out_allocated=.false.
 else
   call get_my_natom(mpi_comm_out,my_natom_out,natom_tot)
   call get_my_atmtab(mpi_comm_out,my_atmtab_out,my_atmtab_out_allocated,&
&                     paral_atom,natom_tot)
 end if

!Select algo according to optional input arguments
 algo_option=1
 if (present(SendAtomProc).and.present(SendAtomList).and.&
&    present(RecvAtomProc).and.present(RecvAtomList)) algo_option=2


!Brute force algorithm (allgather + scatter)
!---------------------------------------------------------
 if (algo_option==1) then

   LIBPAW_DATATYPE_ALLOCATE(paw_ij_all,(natom_tot))
   call paw_ij_nullify(paw_ij_all)
   call paw_ij_copy(paw_ij,paw_ij_all,comm_atom=mpi_comm_in,mpi_atmtab=my_atmtab_in)
   if (in_place) then
    call paw_ij_free(paw_ij)
    LIBPAW_DATATYPE_DEALLOCATE(paw_ij)
    LIBPAW_DATATYPE_ALLOCATE(paw_ij,(my_natom_out))
    call paw_ij_nullify(paw_ij)
    call paw_ij_copy(paw_ij_all,paw_ij,comm_atom=mpi_comm_out,mpi_atmtab=my_atmtab_out)
   else
     LIBPAW_DATATYPE_ALLOCATE(paw_ij_out,(my_natom_out))
     call paw_ij_nullify(paw_ij_out)
     call paw_ij_copy(paw_ij_all,paw_ij_out,comm_atom=mpi_comm_out,mpi_atmtab=my_atmtab_out)
   end if
   call paw_ij_free(paw_ij_all)
   LIBPAW_DATATYPE_DEALLOCATE(paw_ij_all)


!Asynchronous algorithm (asynchronous communications)
!---------------------------------------------------------
 else if (algo_option==2) then

   nbsend=size(SendAtomProc) ; nbrecv=size(RecvAtomProc)

   if (in_place) then
     if (my_natom_out > 0) then
       LIBPAW_DATATYPE_ALLOCATE(paw_ij_out1,(my_natom_out))
       call paw_ij_nullify(paw_ij_out1)
     else
       LIBPAW_DATATYPE_ALLOCATE(paw_ij_out1,(0))
     end if
   else
     LIBPAW_DATATYPE_ALLOCATE(paw_ij_out,(my_natom_out))
     call paw_ij_nullify(paw_ij_out)
      paw_ij_out1=>paw_ij_out
   end if

   nproc_in=xmpi_comm_size(mpi_comm_in)
   nproc_out=xmpi_comm_size(mpi_comm_out)
   if (nproc_in<=nproc_out) mpi_comm_exch=mpi_comm_out
   if (nproc_in>nproc_out) mpi_comm_exch=mpi_comm_in
   me_exch=xmpi_comm_rank(mpi_comm_exch)


!  Dimension put to the maximum to send
   LIBPAW_ALLOCATE(atmtab_send,(nbsend))
   LIBPAW_ALLOCATE(atm_indx_in,(natom_tot))
   atm_indx_in=-1
   do iatom=1,my_natom_in
     atm_indx_in(my_atmtab_in(iatom))=iatom
   end do
   LIBPAW_ALLOCATE(atm_indx_out,(natom_tot))
   atm_indx_out=-1
   do iatom=1,my_natom_out
     atm_indx_out(my_atmtab_out(iatom))=iatom
   end do

   LIBPAW_DATATYPE_ALLOCATE(tab_buf_int,(nbsend))
   LIBPAW_DATATYPE_ALLOCATE(tab_buf_dp,(nbsend))
   LIBPAW_DATATYPE_ALLOCATE(tab_buf_atom,(nbsend))
   LIBPAW_ALLOCATE(request,(3*nbsend))

!  A send buffer in an asynchrone communication couldn't be deallocate before it has been receive
   nbsent=0 ; ireq=0 ; iisend=0 ; nbsendreq=0 ; nb_msg=0
   do iisend=1,nbsend
     iproc_rcv=SendAtomProc(iisend)
     next=-1
     if (iisend < nbsend) next=SendAtomProc(iisend+1)
     if (iproc_rcv /= me_exch) then
       nbsent=nbsent+1
       atmtab_send(nbsent)=SendAtomList(iisend) ! we groups the atoms sends to the same process
       if (iproc_rcv /= next) then
         if (nbsent > 0 ) then
!          Check if message has been yet prepared
           message_yet_prepared=.false.
           do imsg=1,nb_msg
             if (size(tab_buf_atom(imsg)%value) /= nbsent) then
               cycle
             else
               do imsg1=1,nbsent
                 if (tab_buf_atom(imsg)%value(imsg1)/=atmtab_send(imsg1)) exit
                 message_yet_prepared=.true.
                 imsg_current=imsg
               end do
             end if
           enddo
!          Create the message
           if (.not.message_yet_prepared) then
             nb_msg=nb_msg+1
             call paw_ij_isendreceive_fillbuffer( &
&                   paw_ij,atmtab_send,atm_indx_in,nbsent,buf_int,nb_int,buf_dp,nb_dp)
             LIBPAW_ALLOCATE(tab_buf_int(nb_msg)%value,(nb_int))
             LIBPAW_ALLOCATE(tab_buf_dp(nb_msg)%value,(nb_dp))
             tab_buf_int(nb_msg)%value(1:nb_int)=buf_int(1:nb_int)
             tab_buf_dp(nb_msg)%value(1:nb_dp)=buf_dp(1:nb_dp)
             LIBPAW_DEALLOCATE(buf_int)
             LIBPAW_DEALLOCATE(buf_dp)
             LIBPAW_ALLOCATE(tab_buf_atom(nb_msg)%value, (nbsent))
             tab_buf_atom(nb_msg)%value(1:nbsent)=atmtab_send(1:nbsent)
             imsg_current=nb_msg
           end if
!          Communicate
           buf_size(1)=size(tab_buf_int(imsg_current)%value)
           buf_size(2)=size(tab_buf_dp(imsg_current)%value)
           buf_size(3)=nbsent
           buf_ints=>tab_buf_int(imsg_current)%value
           buf_dps=>tab_buf_dp(imsg_current)%value
           my_tag=200
           ireq=ireq+1
           call xmpi_isend(buf_size,iproc_rcv,my_tag,mpi_comm_exch,request(ireq),ierr)
           my_tag=201
           ireq=ireq+1
           call xmpi_isend(buf_ints,iproc_rcv,my_tag,mpi_comm_exch,request(ireq),ierr)
           my_tag=202
           ireq=ireq+1
           call xmpi_isend(buf_dps,iproc_rcv,my_tag,mpi_comm_exch,request(ireq),ierr)
           nbsendreq=ireq
           nbsent=0
         end if
       end if
     else ! Just a renumbering, not a sending
       iat_in=atm_indx_in(SendAtomList(iisend))
       iat_out=atm_indx_out(my_atmtab_in(iat_in))
       call paw_ij_copy(paw_ij(iat_in:iat_in),paw_ij_out1(iat_out:iat_out))
       nbsent=0
     end if
   end do

   LIBPAW_ALLOCATE(From,(nbrecv))
   From(:)=-1 ; nbrecvmsg=0
   do iircv=1,nbrecv
     iproc_send=RecvAtomProc(iircv) !receive from (RcvAtomProc is sorted by growing process)
     next=-1
     if (iircv < nbrecv) next=RecvAtomProc(iircv+1)
     if (iproc_send /= me_exch .and. iproc_send/=next) then
       nbrecvmsg=nbrecvmsg+1
       From(nbrecvmsg)=iproc_send
     end if
   end do

   LIBPAW_ALLOCATE(msg_pick,(nbrecvmsg))
   msg_pick=.false.
   nbmsg_incoming=nbrecvmsg
   do while (nbmsg_incoming > 0)
     do i1=1,nbrecvmsg
       if (.not.msg_pick(i1)) then
         iproc_send=From(i1)
         flag=.false.
         my_tag=200
         call xmpi_iprobe(iproc_send,my_tag,mpi_comm_exch,flag,ierr)
         if (flag) then
           msg_pick(i1)=.true.
           call xmpi_irecv(buf_size,iproc_send,my_tag,mpi_comm_exch,request1(1),ierr)
           call xmpi_wait(request1(1),ierr)
           nb_int=buf_size(1)
           nb_dp=buf_size(2)
           npaw_ij_sent=buf_size(3)
           LIBPAW_ALLOCATE(buf_int1,(nb_int))
           LIBPAW_ALLOCATE(buf_dp1,(nb_dp))
           my_tag=201
           call xmpi_irecv(buf_int1,iproc_send,my_tag,mpi_comm_exch,request1(2),ierr)
           my_tag=202
           call xmpi_irecv(buf_dp1,iproc_send,my_tag,mpi_comm_exch,request1(3),ierr)
           call xmpi_waitall(request1(2:3),ierr)
           call paw_ij_isendreceive_getbuffer(paw_ij_out1,npaw_ij_sent,atm_indx_out,buf_int1,buf_dp1)
           nbmsg_incoming=nbmsg_incoming-1
           LIBPAW_DEALLOCATE(buf_int1)
           LIBPAW_DEALLOCATE(buf_dp1)
         end if
       end if
     end do
   end do
   LIBPAW_DEALLOCATE(msg_pick)

   if (in_place) then
     call paw_ij_free(paw_ij)
     LIBPAW_DATATYPE_DEALLOCATE(paw_ij)
     LIBPAW_DATATYPE_ALLOCATE(paw_ij,(my_natom_out))
     call paw_ij_copy(paw_ij_out1,paw_ij)
     call paw_ij_free(paw_ij_out1)
    LIBPAW_DATATYPE_DEALLOCATE(paw_ij_out1)
  end if

!  Wait for deallocating arrays that all sending operations has been realized
   if (nbsendreq > 0) then
     call xmpi_waitall(request(1:nbsendreq),ierr)
   end if

!  Deallocate buffers
   do i1=1,nb_msg
     LIBPAW_DEALLOCATE(tab_buf_int(i1)%value)
     LIBPAW_DEALLOCATE(tab_buf_dp(i1)%value)
     LIBPAW_DEALLOCATE(tab_buf_atom(i1)%value)
   end do
   LIBPAW_DATATYPE_DEALLOCATE(tab_buf_int)
   LIBPAW_DATATYPE_DEALLOCATE(tab_buf_dp)
   LIBPAW_DATATYPE_DEALLOCATE(tab_buf_atom)
   LIBPAW_DEALLOCATE(From)
   LIBPAW_DEALLOCATE(request)
   LIBPAW_DEALLOCATE(atmtab_send)
   LIBPAW_DEALLOCATE(atm_indx_in)
   LIBPAW_DEALLOCATE(atm_indx_out)

 end if !algo_option

!Eventually release temporary pointers
 call free_my_atmtab(my_atmtab_in,my_atmtab_in_allocated)
 call free_my_atmtab(my_atmtab_out,my_atmtab_out_allocated)

end subroutine paw_ij_redistribute
!!***

!----------------------------------------------------------------------

!!****f* m_paw_ij/paw_ij_reset_flags
!! NAME
!! paw_ij_reset_flags
!!
!! FUNCTION
!!  Set all paw_ij flags to 1 (force the recomputation of all arrays)
!!
!! SIDE EFFECTS
!!  Paw_ij<type(paw_ij_type)>=paw_ij structure
!!
!! PARENTS
!!      d2frnl,dfpt_nstpaw,dfpt_scfcv,scfcv
!!
!! CHILDREN
!!
!! SOURCE

subroutine paw_ij_reset_flags(Paw_ij,all,dijhartree,self_consistent)

!Arguments ------------------------------------
!scalars
 logical,optional,intent(in) :: all,dijhartree,self_consistent
!arrays
 type(Paw_ij_type),intent(inout) :: Paw_ij(:)

!Local variables-------------------------------
 integer :: iat,natom
 logical :: all_,dijhartree_,self_consistent_

! *************************************************************************

!@Paw_ij_type

 natom=SIZE(Paw_ij);if (natom==0) return
 all_=.true.;if (present(all)) all_=all
 dijhartree_=.false.;if (present(dijhartree)) dijhartree_=dijhartree
 self_consistent_=.false.;if (present(self_consistent)) self_consistent_=self_consistent

 if (dijhartree_) then
   do iat=1,natom
     if (Paw_ij(iat)%has_dijhartree>0) Paw_ij(iat)%has_dijhartree=1
   end do

 else if (self_consistent_) then
   do iat=1,natom
     if (Paw_ij(iat)%has_dij       >0) Paw_ij(iat)%has_dij       =1
     if (Paw_ij(iat)%has_dijexxc   >0) Paw_ij(iat)%has_dijexxc   =1
     if (Paw_ij(iat)%has_dijfock   >0) Paw_ij(iat)%has_dijfock   =1
     if (Paw_ij(iat)%has_dijhartree>0) Paw_ij(iat)%has_dijhartree=1
     if (Paw_ij(iat)%has_dijhat    >0) Paw_ij(iat)%has_dijhat    =1
     if (Paw_ij(iat)%has_dijnd     >0) Paw_ij(iat)%has_dijnd     =1
     if (Paw_ij(iat)%has_dijso     >0) Paw_ij(iat)%has_dijso     =1
     if (Paw_ij(iat)%has_dijU      >0) Paw_ij(iat)%has_dijU      =1
     if (Paw_ij(iat)%has_dijxc     >0) Paw_ij(iat)%has_dijxc     =1
     if (Paw_ij(iat)%has_dijxc_hat >0) Paw_ij(iat)%has_dijxc_hat =1
     if (Paw_ij(iat)%has_dijxc_val >0) Paw_ij(iat)%has_dijxc_val =1
     if (Paw_ij(iat)%has_exexch_pot>0) Paw_ij(iat)%has_exexch_pot=1
     if (Paw_ij(iat)%has_pawu_occ  >0) Paw_ij(iat)%has_pawu_occ  =1
   end do

 else if (all_) then
   do iat=1,natom
     if (Paw_ij(iat)%has_dij       >0) Paw_ij(iat)%has_dij       =1
     if (Paw_ij(iat)%has_dij0      >0) Paw_ij(iat)%has_dij0      =1
     if (Paw_ij(iat)%has_dijexxc   >0) Paw_ij(iat)%has_dijexxc   =1
     if (Paw_ij(iat)%has_dijfock   >0) Paw_ij(iat)%has_dijfock   =1
     if (Paw_ij(iat)%has_dijfr     >0) Paw_ij(iat)%has_dijfr     =1
     if (Paw_ij(iat)%has_dijhartree>0) Paw_ij(iat)%has_dijhartree=1
     if (Paw_ij(iat)%has_dijhat    >0) Paw_ij(iat)%has_dijhat    =1
     if (Paw_ij(iat)%has_dijnd     >0) Paw_ij(iat)%has_dijnd     =1
     if (Paw_ij(iat)%has_dijso     >0) Paw_ij(iat)%has_dijso     =1
     if (Paw_ij(iat)%has_dijU      >0) Paw_ij(iat)%has_dijU      =1
     if (Paw_ij(iat)%has_dijxc     >0) Paw_ij(iat)%has_dijxc     =1
     if (Paw_ij(iat)%has_dijxc_hat >0) Paw_ij(iat)%has_dijxc_hat =1
     if (Paw_ij(iat)%has_dijxc_val >0) Paw_ij(iat)%has_dijxc_val =1
     if (Paw_ij(iat)%has_exexch_pot>0) Paw_ij(iat)%has_exexch_pot=1
     if (Paw_ij(iat)%has_pawu_occ  >0) Paw_ij(iat)%has_pawu_occ  =1
   end do
 end if

end subroutine paw_ij_reset_flags
!!***

!----------------------------------------------------------------------

!!****f* m_paw_ij/paw_ij_isendreceive_getbuffer
!! NAME
!!  paw_ij_isendreceive_getbuffer
!!
!! FUNCTION
!!  Fill a paw_ij structure with the buffers received in a receive operation
!!  This buffer should have been first extracted by a call to paw_ij_isendreceive_fillbuffer
!!
!! INPUTS
!!  atm_indx_recv(1:total number of atoms)= array for receive operation
!!                 Given an index of atom in global numbering, give its index
!!                 in the table of atoms treated by current processor
!!                 or -1 if the atoms is not treated by current processor
!!  buf_int= buffer of receive integers
!!  buf_dp= buffer of receive double precision numbers
!!  npaw_ij_send= number of sent atoms
!!
!! OUTPUT
!!  paw_ij= output datastructure filled with buffers receive in a receive operation
!!
!! PARENTS
!!      m_paw_ij
!!
!! CHILDREN
!!
!! SOURCE

subroutine paw_ij_isendreceive_getbuffer(paw_ij,npaw_ij_send,atm_indx_recv,buf_int,buf_dp)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npaw_ij_send
!arrays
 integer,intent(in):: atm_indx_recv(:),buf_int(:)
 real(dp),intent(in):: buf_dp(:)
 type(paw_ij_type),target,intent(inout) :: paw_ij(:)

!Local variables-------------------------------
!scalars
 integer :: buf_dp_size,buf_int_size
 integer :: cplxq_lmn2_size,cplxdij_lmn2_size,cplxdijq_lmn2_size
 integer :: iat,iatom_tot,ii,ij,indx_dp,indx_int,ndij,nocc,nocc1,nocc2,nocc3,nocc4
 integer :: lmn2_size,sz1,sz2,sz3
 character(len=500) :: msg
 type(Paw_ij_type),pointer :: paw_ij1
!arrays

! *********************************************************************

 buf_int_size=size(buf_int)
 buf_dp_size=size(buf_dp)
 indx_int=1;indx_dp=1

 do ij=1,npaw_ij_send
   iatom_tot=buf_int(indx_int) ;indx_int=indx_int+1
   iat= atm_indx_recv(iatom_tot)
   paw_ij1=>paw_ij(iat)
   paw_ij1%cplex_dij=buf_int(indx_int) ;indx_int=indx_int+1
   paw_ij1%qphase=buf_int(indx_int) ;indx_int=indx_int+1
   paw_ij1%itypat=buf_int(indx_int) ;indx_int=indx_int+1
   paw_ij1%nspden=buf_int(indx_int) ;indx_int=indx_int+1
   paw_ij1%nsppol=buf_int(indx_int) ;indx_int=indx_int+1
   paw_ij1%lmn_size=buf_int(indx_int) ;indx_int=indx_int+1
   paw_ij1%lmn2_size=buf_int(indx_int) ;indx_int=indx_int+1
   paw_ij1%ndij=buf_int(indx_int) ;indx_int=indx_int+1
   paw_ij1%has_dij=buf_int(indx_int) ;indx_int=indx_int+1
   paw_ij1%has_dij0=buf_int(indx_int) ;indx_int=indx_int+1
   paw_ij1%has_dijexxc=buf_int(indx_int) ;indx_int=indx_int+1
   paw_ij1%has_dijfock=buf_int(indx_int) ;indx_int=indx_int+1
   paw_ij1%has_dijfr=buf_int(indx_int) ;indx_int=indx_int+1
   paw_ij1%has_dijhartree=buf_int(indx_int) ;indx_int=indx_int+1
   paw_ij1%has_dijhat=buf_int(indx_int) ;indx_int=indx_int+1
   paw_ij1%has_dijnd=buf_int(indx_int) ;indx_int=indx_int+1
   paw_ij1%has_dijso=buf_int(indx_int) ;indx_int=indx_int+1
   paw_ij1%has_dijU=buf_int(indx_int) ;indx_int=indx_int+1
   paw_ij1%has_dijxc=buf_int(indx_int) ;indx_int=indx_int+1
   paw_ij1%has_dijxc_hat=buf_int(indx_int) ;indx_int=indx_int+1
   paw_ij1%has_dijxc_val=buf_int(indx_int) ;indx_int=indx_int+1
   paw_ij1%has_exexch_pot=buf_int(indx_int) ;indx_int=indx_int+1
   paw_ij1%has_pawu_occ=buf_int(indx_int) ;indx_int=indx_int+1
   if (paw_ij1%has_pawu_occ>=1) then
     nocc1=buf_int(indx_int) ;indx_int=indx_int+1
     nocc2=buf_int(indx_int) ;indx_int=indx_int+1
     nocc3=buf_int(indx_int) ;indx_int=indx_int+1
     nocc4=buf_int(indx_int) ;indx_int=indx_int+1
   else
     nocc1=0;nocc2=0;nocc3=0;nocc4=0
   end if
   lmn2_size=paw_ij1%lmn2_size
   cplxdij_lmn2_size=paw_ij1%cplex_dij*lmn2_size
   cplxq_lmn2_size=paw_ij1%qphase*lmn2_size
   cplxdijq_lmn2_size=cplxq_lmn2_size*paw_ij1%qphase
   ndij=paw_ij1%ndij

   if (paw_ij1%has_dij>=1) then
     ii=cplxdijq_lmn2_size
     LIBPAW_ALLOCATE(paw_ij1%dij,(ii,ndij))
     if (paw_ij1%has_dij==2) then
       paw_ij1%dij(:,:)= &
 &       reshape(buf_dp(indx_dp:indx_dp+ii*ndij-1),(/ii,ndij/))
       indx_dp=indx_dp+ii*ndij
     end if
   end if
   if (paw_ij1%has_dij0 >=1) then
     ii=lmn2_size
     LIBPAW_ALLOCATE(paw_ij1%dij0,(ii))
     if (paw_ij1%has_dij0==2) then
       paw_ij1%dij0(:)=buf_dp(indx_dp:indx_dp+ii-1)
       indx_dp=indx_dp+ii
     end if
   end if
   if (paw_ij1%has_dijexxc >=1) then
     ii=cplxdij_lmn2_size
     LIBPAW_ALLOCATE(paw_ij1%dijexxc,(ii,ndij))
     if (paw_ij1%has_dijexxc==2) then
       paw_ij1%dijexxc(:,:)= &
&        reshape(buf_dp(indx_dp:indx_dp+ii*ndij-1),(/ii,ndij/))
       indx_dp=indx_dp+ii*ndij
     end if
   end if
   if (paw_ij1%has_dijfock >=1) then
     ii=cplxdij_lmn2_size
     LIBPAW_ALLOCATE(paw_ij1%dijfock,(ii,ndij))
     if (paw_ij1%has_dijfock==2) then
       paw_ij1%dijfock(:,:)= &
&        reshape(buf_dp(indx_dp:indx_dp+ii*ndij-1),(/ii,ndij/))
       indx_dp=indx_dp+ii*ndij
     end if
   end if
   if (paw_ij1%has_dijfr >=1) then
     ii=cplxdijq_lmn2_size
     LIBPAW_ALLOCATE(paw_ij1%dijfr,(ii,ndij))
     if (paw_ij1%has_dijfr==2) then
       paw_ij1%dijfr(:,:)= &
&        reshape(buf_dp(indx_dp:indx_dp+ii*ndij-1),(/ii,ndij/))
       indx_dp=indx_dp+ii*ndij
     end if
   end if
   if (paw_ij1%has_dijhartree >=1) then
     ii=cplxq_lmn2_size
     LIBPAW_ALLOCATE(paw_ij1%dijhartree,(ii))
     if (paw_ij1%has_dijhartree==2) then
       paw_ij1%dijhartree(:)=buf_dp(indx_dp:indx_dp+ii-1)
       indx_dp=indx_dp+ii
     end if
   end if
   if (paw_ij1%has_dijhat >=1) then
     ii=cplxdijq_lmn2_size
     LIBPAW_ALLOCATE(paw_ij1%dijhat,(ii,ndij))
     if (paw_ij1%has_dijhat==2) then
       paw_ij1%dijhat(:,:)= &
&        reshape(buf_dp(indx_dp:indx_dp+ii*ndij-1),(/ii,ndij/))
       indx_dp=indx_dp+ii*ndij
     end if
   end if
   if (paw_ij1%has_dijnd >=1) then
     ii=cplxdij_lmn2_size
     LIBPAW_ALLOCATE(paw_ij1%dijnd,(ii,ndij))
     if (paw_ij1%has_dijnd==2) then
       paw_ij1%dijnd(:,:)= &
&        reshape(buf_dp(indx_dp:indx_dp+ii*ndij-1),(/ii,ndij/))
       indx_dp=indx_dp+ii*ndij
     end if
   end if
   if (paw_ij1%has_dijso >=1) then
     ii=cplxdijq_lmn2_size
     LIBPAW_ALLOCATE(paw_ij1%dijso,(ii,ndij))
     if (paw_ij1%has_dijso==2) then
       paw_ij1%dijso(:,:)= &
&        reshape(buf_dp(indx_dp:indx_dp+ii*ndij-1),(/ii,ndij/))
       indx_dp=indx_dp+ii*ndij
     end if
   end if
   if (paw_ij1%has_dijU >=1) then
     ii=cplxdijq_lmn2_size
     LIBPAW_ALLOCATE(paw_ij1%dijU,(ii,ndij))
     if (paw_ij1%has_dijU==2) then
       paw_ij1%dijU(:,:)= &
&        reshape(buf_dp(indx_dp:indx_dp+ii*ndij-1),(/ii,ndij/))
       indx_dp=indx_dp+ii*ndij
     end if
   end if
   if (paw_ij1%has_dijxc >=1) then
     ii=cplxdijq_lmn2_size
     LIBPAW_ALLOCATE(paw_ij1%dijxc,(ii,ndij))
     if (paw_ij1%has_dijxc==2) then
       paw_ij1%dijxc(:,:)= &
&        reshape(buf_dp(indx_dp:indx_dp+ii*ndij-1),(/ii,ndij/))
       indx_dp=indx_dp+ii*ndij
     end if
   end if
   if (paw_ij1%has_dijxc_hat >=1) then
     ii=cplxdij_lmn2_size
     LIBPAW_ALLOCATE(paw_ij1%dijxc_hat,(ii,ndij))
       if (paw_ij1%has_dijxc_hat==2) then
       paw_ij1%dijxc_hat(:,:)= &
&        reshape(buf_dp(indx_dp:indx_dp+ii*ndij-1),(/ii,ndij/))
         indx_dp=indx_dp+ii*ndij
     end if
   end if
   if (paw_ij1%has_dijxc_val >=1) then
     ii=cplxdij_lmn2_size
     LIBPAW_ALLOCATE(paw_ij1%dijxc_val,(ii,ndij))
     if (paw_ij1%has_dijxc_val==2) then
       paw_ij1%dijxc_val(:,:)= &
&        reshape(buf_dp(indx_dp:indx_dp+ii*ndij-1),(/ii,ndij/))
       indx_dp=indx_dp+ii*ndij
     end if
   end if
   if (paw_ij1%has_pawu_occ >=1) then
     nocc=paw_ij1%ndij
     LIBPAW_ALLOCATE(paw_ij1%nocctot,(nocc))
     paw_ij1%nocctot(1:nocc)=buf_dp(indx_dp:indx_dp+nocc-1)
     indx_dp=indx_dp+nocc
     nocc=nocc1*nocc2*nocc3*nocc4
     LIBPAW_ALLOCATE(paw_ij1%noccmmp,(nocc1,nocc2,nocc3,nocc4))
     paw_ij1%noccmmp(1:nocc1,1:nocc2,1:nocc3,1:nocc4)= &
&      reshape(buf_dp(indx_dp:indx_dp+nocc-1),(/nocc1,nocc2,nocc3,nocc4/))
     indx_dp=indx_dp+nocc
   end if
   if (paw_ij1%has_exexch_pot >=1) then
     sz1=buf_int(indx_int);indx_int=indx_int+1
     sz2=buf_int(indx_int);indx_int=indx_int+1
     sz3=buf_int(indx_int);indx_int=indx_int+1
     LIBPAW_ALLOCATE(paw_ij1%vpawx,(sz1,sz2,sz3))
     if (paw_ij1%has_exexch_pot == 2) then
       paw_ij1%vpawx(:,:,:)=&
&        reshape(buf_dp(indx_dp:indx_dp+sz1*sz2*sz3-1),(/sz1,sz2,sz3/))
       indx_dp=indx_dp+sz1*sz2*sz3
     end if
   end if
 end do

 if ((indx_int/=1+buf_int_size).or.(indx_dp /=1+buf_dp_size)) then
   write(msg,'(a,i10,a,i10)') 'Wrong buffer sizes: buf_int_size=',buf_int_size,' buf_dp_size=',buf_dp_size
   MSG_BUG(msg)
 end if

end subroutine paw_ij_isendreceive_getbuffer
!!***

!----------------------------------------------------------------------

!!****f* m_paw_an/paw_ij_isendreceive_fillbuffer
!! NAME
!!  paw_ij_isendreceive_fillbuffer
!!
!! FUNCTION
!!  Extract from paw_ij and from the global index of atoms
!!  the buffers to send in a sending operation
!!  This function has to be coupled with a call to paw_ij_isendreceive_getbuffer
!!
!! INPUTS
!!  atm_indx_send(1:total number of atoms)= array for send operation,
!!                 Given an index of atom in global numbering, give its index
!!                 in the table of atoms treated by current processor
!!                 or -1 if the atoms is not treated by current processor
!!  npaw_ij_send= number of sent atoms
!!  paw_ij= data structure from which are extract buffer int and buffer dp
!!
!! OUTPUT
!!  buf_int= buffer of integers to be sent
!!  buf_int_size= size of buffer of integers
!!  buf_dp= buffer of double precision numbers to be sent
!!  buf_dp_size= size of buffer of double precision numbers
!!
!! PARENTS
!!      m_paw_ij
!!
!! CHILDREN
!!
!! SOURCE

subroutine paw_ij_isendreceive_fillbuffer(paw_ij,atmtab_send,atm_indx_send,npaw_ij_send,&
&                                         buf_int,buf_int_size,buf_dp,buf_dp_size)

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: buf_int_size,buf_dp_size
 integer,intent(in) :: npaw_ij_send
!arrays
 integer,intent(in) :: atmtab_send(:),atm_indx_send(:)
 integer,allocatable,intent(out) :: buf_int(:)
 real(dp),allocatable,intent(out):: buf_dp(:)
 type(paw_ij_type),target,intent(in) :: paw_ij(:)

!Local variables-------------------------------
!scalars
 integer :: cplxdij_lmn2_size,cplxdijq_lmn2_size,cplxq_lmn2_size
 integer :: iatom_tot,ii,ij,indx_dp,indx_int,ipaw_ij_send
 integer :: lmn2_size,ndij,nocc,nspden,sz1,sz2,sz3
 character(len=500) :: msg
 type(Paw_ij_type),pointer :: paw_ij1
!arrays

! *********************************************************************

!Compute sizes of buffers
 buf_int_size=0;buf_dp_size=0
 do ipaw_ij_send=1,npaw_ij_send
   iatom_tot=atmtab_send(ipaw_ij_send)
   ij = atm_indx_send(iatom_tot)
   paw_ij1=>paw_ij(ij)
   lmn2_size=paw_ij1%lmn2_size
   cplxdij_lmn2_size=paw_ij1%cplex_dij*lmn2_size
   cplxq_lmn2_size=paw_ij1%qphase*lmn2_size
   cplxdijq_lmn2_size=cplxdij_lmn2_size*paw_ij1%qphase
   ndij=paw_ij1%ndij
   buf_int_size=buf_int_size+24
   if (paw_ij1%has_dij==2) then
     buf_dp_size=buf_dp_size +cplxdijq_lmn2_size*ndij
   end if
   if (paw_ij1%has_dij0==2) then
     buf_dp_size=buf_dp_size +lmn2_size
   end if
   if (paw_ij1%has_dijexxc==2) then
     buf_dp_size=buf_dp_size +cplxdij_lmn2_size*ndij
   end if
   if (paw_ij1%has_dijfock==2) then
     buf_dp_size=buf_dp_size +cplxdij_lmn2_size*ndij
   end if
   if (paw_ij1%has_dijfr==2) then
     buf_dp_size=buf_dp_size +cplxdijq_lmn2_size*ndij
   end if
   if (paw_ij1%has_dijhartree==2) then
     buf_dp_size=buf_dp_size +cplxq_lmn2_size
   end if
   if (paw_ij1%has_dijhat==2) then
     buf_dp_size=buf_dp_size +cplxdijq_lmn2_size*ndij
   end if
   if (paw_ij1%has_dijnd==2) then
     buf_dp_size=buf_dp_size +cplxdij_lmn2_size*ndij
   end if
   if (paw_ij1%has_dijso==2) then
     buf_dp_size=buf_dp_size +cplxdijq_lmn2_size*ndij
   end if
   if (paw_ij1%has_dijU==2) then
     buf_dp_size=buf_dp_size +cplxdijq_lmn2_size*ndij
   end if
   if (paw_ij1%has_dijxc==2) then
     buf_dp_size=buf_dp_size +cplxdijq_lmn2_size*ndij
   end if
   if (paw_ij1%has_dijxc_hat==2) then
     buf_dp_size=buf_dp_size +cplxdij_lmn2_size*ndij
   end if
   if (paw_ij1%has_dijxc_val==2) then
     buf_dp_size=buf_dp_size +cplxdij_lmn2_size*ndij
   end if
   if (paw_ij1%has_pawu_occ>=1) then
     buf_int_size=buf_int_size+4
     buf_dp_size=buf_dp_size &
&         +size(paw_ij1%nocctot) &
&         +size(paw_ij1%noccmmp)
   end if
   if (paw_ij1%has_exexch_pot>=1) then
     buf_int_size=buf_int_size+3
     buf_dp_size=buf_dp_size +size(paw_ij1%vpawx)
   end if
 end do

!Fill input buffers
 LIBPAW_ALLOCATE(buf_int,(buf_int_size))
 LIBPAW_ALLOCATE(buf_dp,(buf_dp_size))
 indx_int=1;indx_dp=1
 do ipaw_ij_send=1,npaw_ij_send
   iatom_tot=atmtab_send(ipaw_ij_send)
   ij = atm_indx_send(iatom_tot)
   paw_ij1=>paw_ij(ij)
   nspden=paw_ij1%nspden
   buf_int(indx_int)=iatom_tot ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij1%cplex_dij ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij1%qphase ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij1%itypat ;indx_int=indx_int+1
   buf_int(indx_int)=nspden ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij1%nsppol ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij1%lmn_size ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij1%lmn2_size ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij1%ndij ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij1%has_dij ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij1%has_dij0 ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij1%has_dijexxc ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij1%has_dijfock ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij1%has_dijfr ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij1%has_dijhartree ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij1%has_dijhat ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij1%has_dijnd ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij1%has_dijso ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij1%has_dijU ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij1%has_dijxc ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij1%has_dijxc_hat ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij1%has_dijxc_val ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij1%has_exexch_pot ;indx_int=indx_int+1
   buf_int(indx_int)=paw_ij1%has_pawu_occ ;indx_int=indx_int+1
   lmn2_size=paw_ij1%lmn2_size
   cplxdij_lmn2_size=paw_ij1%cplex_dij*lmn2_size
   cplxq_lmn2_size=paw_ij1%qphase*lmn2_size
   cplxdijq_lmn2_size=cplxdij_lmn2_size*paw_ij1%qphase
   ndij=paw_ij1%ndij
   if (paw_ij1%has_dij==2) then
     ii=cplxdijq_lmn2_size*ndij
     buf_dp(indx_dp:indx_dp+ii-1)=reshape(paw_ij1%dij,(/ii/))
     indx_dp=indx_dp+ii
   end if
   if (paw_ij1%has_dij0==2) then
     ii=lmn2_size
     buf_dp(indx_dp:indx_dp+lmn2_size-1)=paw_ij1%dij0(:)
     indx_dp=indx_dp+lmn2_size
   end if
   if (paw_ij1%has_dijexxc==2) then
     ii=cplxdij_lmn2_size*ndij
     buf_dp(indx_dp:indx_dp+ii-1)=reshape(paw_ij1%dijexxc,(/ii/))
     indx_dp=indx_dp+ii
   end if
   if (paw_ij1%has_dijfock==2) then
     ii=cplxdij_lmn2_size*ndij
     buf_dp(indx_dp:indx_dp+ii-1)=reshape(paw_ij1%dijfock,(/ii/))
     indx_dp=indx_dp+ii
   end if
   if (paw_ij1%has_dijfr==2) then
     ii=cplxdijq_lmn2_size*ndij
     buf_dp(indx_dp:indx_dp+ii-1)=reshape(paw_ij1%dijfr,(/ii/))
     indx_dp=indx_dp+ii
   end if
   if (paw_ij1%has_dijhartree==2) then
     ii=cplxq_lmn2_size
     buf_dp(indx_dp:indx_dp+ii-1)=paw_ij1%dijhartree(:)
     indx_dp=indx_dp+ii
   end if
   if (paw_ij1%has_dijhat==2) then
     ii=cplxdijq_lmn2_size*ndij
     buf_dp(indx_dp:indx_dp+ii-1)=reshape(paw_ij1%dijhat,(/ii/))
     indx_dp=indx_dp+ii
   end if
   if (paw_ij1%has_dijnd==2) then
     ii=cplxdij_lmn2_size*ndij
     buf_dp(indx_dp:indx_dp+ii-1)=reshape(paw_ij1%dijnd,(/ii/))
     indx_dp=indx_dp+ii
   end if
   if (paw_ij1%has_dijso==2) then
     ii=cplxdijq_lmn2_size*ndij
     buf_dp(indx_dp:indx_dp+ii-1)=reshape(paw_ij1%dijso,(/ii/))
     indx_dp=indx_dp+ii
   end if
   if (paw_ij1%has_dijU==2) then
     ii=cplxdijq_lmn2_size*ndij
     buf_dp(indx_dp:indx_dp+ii-1)=reshape(paw_ij1%dijU,(/ii/))
     indx_dp=indx_dp+ii
   end if
   if (paw_ij1%has_dijxc==2) then
     ii=cplxdijq_lmn2_size*ndij
     buf_dp(indx_dp:indx_dp+ii-1)=reshape(paw_ij1%dijxc,(/ii/))
     indx_dp=indx_dp+ii
   end if
   if (paw_ij1%has_dijxc_hat==2) then
     ii=cplxdij_lmn2_size*ndij
     buf_dp(indx_dp:indx_dp+ii-1)=reshape(paw_ij1%dijxc_hat,(/ii/))
     indx_dp=indx_dp+ii
   end if
   if (paw_ij1%has_dijxc_val==2) then
     ii=cplxdij_lmn2_size*ndij
     buf_dp(indx_dp:indx_dp+ii-1)=reshape(paw_ij1%dijxc_val,(/ii/))
     indx_dp=indx_dp+ii
   end if
   if (paw_ij1%has_pawu_occ>=1)then
     buf_int(indx_int)=size(paw_ij1%noccmmp,1) ;indx_int=indx_int+1
     buf_int(indx_int)=size(paw_ij1%noccmmp,2) ;indx_int=indx_int+1
     buf_int(indx_int)=size(paw_ij1%noccmmp,3) ;indx_int=indx_int+1
     buf_int(indx_int)=size(paw_ij1%noccmmp,4) ;indx_int=indx_int+1
     nocc=paw_ij1%ndij
     buf_dp(indx_dp:indx_dp+nocc-1)=paw_ij1%nocctot(1:nocc)
     indx_dp=indx_dp+nocc
     nocc=size(paw_ij1%noccmmp)
     buf_dp(indx_dp:indx_dp+nocc-1)=reshape(paw_ij1%noccmmp,(/nocc/))
     indx_dp=indx_dp+nocc
   end if
   if (paw_ij1%has_exexch_pot>=1) then
     sz1=size(paw_ij1%vpawx,1);sz2=size(paw_ij1%vpawx,2)
     sz3=size(paw_ij1%vpawx,3)
     buf_int(indx_int)=sz1; indx_int=indx_int+1
     buf_int(indx_int)=sz2; indx_int=indx_int+1
     buf_int(indx_int)=sz3; indx_int=indx_int+1
     ii=sz1*sz2*sz3
     buf_dp(indx_dp:indx_dp+ii-1)=reshape(paw_ij1%vpawx,(/(ii)/))
     indx_dp=indx_dp+ii
   end if
 end do
 indx_int=indx_int-1;indx_dp=indx_dp-1
 if ((indx_int/=buf_int_size).or.(indx_dp/=buf_dp_size)) then
   write(msg,'(a,i10,a,i10)') 'Wrong buffer sizes: buf_int_size=',buf_int_size,' buf_dp_size=',buf_dp_size
   MSG_BUG(msg)
 end if

end subroutine paw_ij_isendreceive_fillbuffer
!!***

!----------------------------------------------------------------------

END MODULE m_paw_ij
!!***
