!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_pawrhoij
!! NAME
!!  m_pawrhoij
!!
!! FUNCTION
!!  This module contains the definition of the pawrhoij_type structured datatype,
!!  as well as related functions and methods.
!!  pawrhoij_type variables define rhoij occupancies matrixes used within PAW formalism.
!!
!! COPYRIGHT
!! Copyright (C) 2012-2018 ABINIT group (MT, FJ)
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

MODULE m_pawrhoij

 USE_DEFS
 USE_MSG_HANDLING
 USE_MPI_WRAPPERS
 USE_MEMORY_PROFILING
#ifdef LIBPAW_HAVE_NETCDF
 use netcdf
#endif

 use m_libpaw_tools, only : libpaw_flush, libpaw_to_upper

 use m_paw_io,     only : pawio_print_ij
 use m_pawang,     only : pawang_type
 use m_pawtab,     only : pawtab_type
 use m_paral_atom, only : get_my_atmtab, free_my_atmtab, get_my_natom

 implicit none

 private

!public procedures.
 public :: pawrhoij_alloc
 public :: pawrhoij_free
 public :: pawrhoij_nullify
 public :: pawrhoij_copy
 public :: pawrhoij_gather
 public :: pawrhoij_bcast
 public :: pawrhoij_redistribute
 public :: pawrhoij_io
 public :: pawrhoij_unpack
 public :: pawrhoij_init_unpacked
 public :: pawrhoij_free_unpacked
 public :: pawrhoij_get_nspden
 public :: symrhoij

 public :: pawrhoij_mpisum_unpacked

 interface pawrhoij_mpisum_unpacked
   module procedure pawrhoij_mpisum_unpacked_1D
   module procedure pawrhoij_mpisum_unpacked_2D
 end interface pawrhoij_mpisum_unpacked

!private procedures.
 private :: pawrhoij_isendreceive_getbuffer
 private :: pawrhoij_isendreceive_fillbuffer
!!***

!!****t* m_pawrhoij/pawrhoij_type
!! NAME
!! pawrhoij_type
!!
!! FUNCTION
!! This structured datatype contains rhoij quantities (occucpancies)
!! and related data, used in PAW calculations.
!!
!! SOURCE

 type,public :: pawrhoij_type

!Integer scalars

  integer :: cplex
   ! cplex=1 if rhoij are real, 2 if rhoij are complex
   ! For GS calculations: rhoij are real provided that time-reversal can be used.

  integer :: itypat
   ! itypat=type of the atom

  integer :: lmn_size
   ! Number of (l,m,n) elements for the paw basis

  integer :: lmn2_size
   ! lmn2_size=lmn_size*(lmn_size+1)/2
   ! where lmn_size is the number of (l,m,n) elements for the paw basis

  integer :: lmnmix_sz=0
   ! lmnmix_sz=number of (lmn,lmn_prime) verifying l<=lmix and l_prime<=lmix
   !           i.e. number of rhoij elements being mixed during SCF cycle
   ! lmnmix_sz=0 if mixing data are not used

  integer :: ngrhoij=0
   ! First dimension of array grhoij

  integer :: nrhoijsel=0
   ! nrhoijsel
   ! Number of non-zero value of rhoij
   ! This is the size of rhoijp(:,:) (see below in this datastructure)

  integer :: nspden
   ! Number of spin-density components for rhoij (may be different from nspden for density)

  integer :: nspinor
   ! Number of spinorial components

  integer :: nsppol
   ! Number of independent spin-components

  integer :: use_rhoij_=0
   ! 1 if pawrhoij%rhoij_ is allocated

  integer :: use_rhoijp=0
   ! 1 if pawrhoij%rhoijp and pawrhoij%rhoijselect are allocated

  integer :: use_rhoijres=0
   ! 1 if pawrhoij%rhoijres is allocated

  integer :: use_rhoijim=0
   ! 1 if pawrhoij%rhoijim is allocated

!Integer arrays

  integer, allocatable :: kpawmix(:)
   ! kpawmix(lmnmix_sz)
   ! Indirect array selecting the elements of rhoij
   ! being mixed during SCF cycle

  integer, allocatable :: rhoijselect(:)
   ! rhoijselect(lmn2_size)
   ! Indirect array selecting the non-zero elements of rhoij:
   ! rhoijselect(isel,ispden)=klmn if rhoij(klmn,ispden) is non-zero

!Real (real(dp)) arrays

  real(dp), allocatable :: grhoij (:,:,:)
   ! grhoij(ngrhoij,cplex*lmn2_size,nspden)
   ! Gradients of Rho_ij wrt xred, strains, ... (non-packed storage)

  real(dp), allocatable :: rhoij_ (:,:)
   ! rhoij_(cplex*lmn2_size,nspden)
   ! Array used to (temporary) store Rho_ij in a non-packed storage mode

  real(dp), allocatable :: rhoijp (:,:)
   ! rhoijp(cplex*lmn2_size,nspden)
   ! Augmentation waves occupancies Rho_ij in PACKED STORAGE (only non-zero elements are stored)

  real(dp), allocatable :: rhoijres (:,:)
   ! rhoijres(cplex*lmn2_size,nspden)
   ! Rho_ij residuals during SCF cycle (non-packed storage)

  real(dp), allocatable :: rhoijim (:,:)
   ! rhoijim(lmn2_size,nspden)
   ! Missing term of the imaginary part of Rho_ij (non-packed storage)

 end type pawrhoij_type
!!***

CONTAINS

!===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/pawrhoij_alloc
!! NAME
!! pawrhoij_alloc
!!
!! FUNCTION
!! Initialize and allocate a pawrhoij datastructure
!!
!! INPUTS
!! [comm_atom] = communicator over atoms  (OPTIONAL)
!! cplex=1 if rhoij is REAL,2 if COMPLEX
!! [my_atmtab(:)] = Index of atoms treated by current proc (OPTIONAL)
!! nspden=number of spin-components for rhoij
!! nsppol=number of spinorial components for rhoij
!! nsppol=number of independant spin-components for rhoij
!! typat(:)=types of atoms
!! [lmnsize(:)]=array of (l,m,n) sizes for rhoij for each type of atom (OPTIONAL)
!!              must be present if [pawtab] argument is not passed
!! [ngrhoij]=number of gradients to be allocated (OPTIONAL, default=0)
!! [nlmnmix]=number of rhoij elements to be mixed during SCF cycle (OPTIONAL, default=0)
!! [pawtab(:)] <type(pawtab_type)>=paw tabulated starting data (OPTIONAL)
!!             must be present if [lmnsize(:)] argument is not passed
!! [use_rhoij_]=1 if pawrhoij(:)%rhoij_ has to be allocated (OPTIONAL, default=0)
!! [use_rhoijp]=1 if pawrhoij(:)%rhoijp has to be allocated (OPTIONAL, default=1)
!!              (in that case, pawrhoij%rhoijselect is also allocated)
!! [use_rhoijres]=1 if pawrhoij(:)%rhoijres has to be allocated (OPTIONAL, default=0)
!! [use_rhoijim]=1 if pawrhoij(:)%rhoijim has to be allocated (OPTIONAL, default=0)
!!
!! SIDE EFFECTS
!! pawrhoij(:)<type(pawrhoij_type)>= rhoij datastructure
!!
!! NOTES
!! One of the two optional arguments lmnsize(:) or pawtab(:) must be present !
!! If both are present, only pawtab(:) is used.
!!
!! PARENTS
!!      bethe_salpeter,dfpt_looppert,dfpt_nstpaw,dfpt_rhofermi,dfpt_scfcv
!!      dfpt_vtorho,energy,extraprho,initrhoij,m_electronpositron,m_hdr,m_ioarr
!!      m_pawrhoij,m_qparticles,paw_qpscgw,posdoppler,respfn,screening
!!      setup_bse,setup_positron,setup_screening,setup_sigma,sigma,vtorho
!!      wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawrhoij_alloc(pawrhoij,cplex,nspden,nspinor,nsppol,typat,&
&  lmnsize,ngrhoij,nlmnmix,pawtab,use_rhoij_,use_rhoijp,& ! Optional
&  use_rhoijres,use_rhoijim,comm_atom,mpi_atmtab)         ! Optional


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawrhoij_alloc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nspden,nspinor,nsppol
 integer,optional,intent(in):: comm_atom,ngrhoij,nlmnmix,use_rhoij_,use_rhoijp,use_rhoijres,use_rhoijim
 integer,optional,target,intent(in)  :: mpi_atmtab(:)
!arrays
 integer,intent(in) :: typat(:)
 integer,optional,target,intent(in) :: lmnsize(:)
 type(pawrhoij_type),intent(inout) :: pawrhoij(:)
 type(pawtab_type),optional,intent(in)  :: pawtab(:)

!Local variables-------------------------------
!scalars
 integer :: irhoij,irhoij_,itypat,lmn2_size,my_comm_atom,nn1,natom,nrhoij
 logical :: has_rhoijp,my_atmtab_allocated,paral_atom
 character(len=500) :: msg
!array
 integer,pointer :: lmn_size(:),my_atmtab(:)

! *************************************************************************

 nrhoij=size(pawrhoij);natom=size(typat)
 if (nrhoij>natom) then
   msg=' wrong sizes (1) !'
   MSG_BUG(msg)
 end if

!Select lmn_size for each atom type
 if (present(pawtab)) then
   nn1=size(pawtab)
   if (maxval(typat)>nn1) then
     msg=' wrong sizes (2) !'
     MSG_BUG(msg)
   end if
   LIBPAW_POINTER_ALLOCATE(lmn_size,(nn1))
   do itypat=1,nn1
     lmn_size(itypat)=pawtab(itypat)%lmn_size
   end do
 else if (present(lmnsize)) then
   nn1=size(lmnsize)
   if (maxval(typat)>nn1) then
     msg=' wrong sizes (3) !'
     MSG_BUG(msg)
   end if
   lmn_size => lmnsize
 else
   msg=' one of the 2 arguments pawtab or lmnsize must be present !'
   MSG_BUG(msg)
 end if

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(nrhoij/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom)

 if (nrhoij>0) then
   do irhoij=1,nrhoij
     irhoij_=irhoij;if (paral_atom) irhoij_=my_atmtab(irhoij)
     itypat=typat(irhoij_)

     lmn2_size=lmn_size(itypat)*(lmn_size(itypat)+1)/2

!    Scalars initializations
     pawrhoij(irhoij)%cplex=cplex
     pawrhoij(irhoij)%itypat=itypat
     pawrhoij(irhoij)%lmn_size=lmn_size(itypat)
     pawrhoij(irhoij)%lmn2_size=lmn2_size
     pawrhoij(irhoij)%nspden=nspden
     pawrhoij(irhoij)%nspinor=nspinor
     pawrhoij(irhoij)%nsppol=nsppol
     pawrhoij(irhoij)%nrhoijsel=0
     pawrhoij(irhoij)%lmnmix_sz=0
     pawrhoij(irhoij)%ngrhoij=0
     pawrhoij(irhoij)%use_rhoij_=0
     pawrhoij(irhoij)%use_rhoijres=0
     pawrhoij(irhoij)%use_rhoijim=0

!    Arrays allocations
     has_rhoijp=.true.; if (present(use_rhoijp)) has_rhoijp=(use_rhoijp>0)
     if (has_rhoijp) then
       pawrhoij(irhoij)%use_rhoijp=1
       LIBPAW_ALLOCATE(pawrhoij(irhoij)%rhoijselect,(lmn2_size))
       LIBPAW_ALLOCATE(pawrhoij(irhoij)%rhoijp,(cplex*lmn2_size,nspden))
       pawrhoij(irhoij)%rhoijselect(:)=0
       pawrhoij(irhoij)%rhoijp(:,:)=zero
     end if

     if (present(ngrhoij)) then
       if (ngrhoij>0) then
         pawrhoij(irhoij)%ngrhoij=ngrhoij
         LIBPAW_ALLOCATE(pawrhoij(irhoij)%grhoij,(ngrhoij,cplex*lmn2_size,nspden))
         pawrhoij(irhoij)%grhoij=zero
       end if
     end if
     if (present(nlmnmix)) then
       if (nlmnmix>0) then
         pawrhoij(irhoij)%lmnmix_sz=nlmnmix
         LIBPAW_ALLOCATE(pawrhoij(irhoij)%kpawmix,(nlmnmix))
         pawrhoij(irhoij)%kpawmix=0
       end if
     end if
     if (present(use_rhoij_)) then
       if (use_rhoij_>0) then
         pawrhoij(irhoij)%use_rhoij_=use_rhoij_
         LIBPAW_ALLOCATE(pawrhoij(irhoij)%rhoij_,(cplex*lmn2_size,nspden))
         pawrhoij(irhoij)%rhoij_=zero
       end if
     end if
     if (present(use_rhoijres)) then
       if (use_rhoijres>0) then
         pawrhoij(irhoij)%use_rhoijres=use_rhoijres
         LIBPAW_ALLOCATE(pawrhoij(irhoij)%rhoijres,(cplex*lmn2_size,nspden))
         pawrhoij(irhoij)%rhoijres=zero
       end if
     end if
     if (present(use_rhoijim)) then
       if (use_rhoijim>0) then
         pawrhoij(irhoij)%use_rhoijim=use_rhoijim
         LIBPAW_ALLOCATE(pawrhoij(irhoij)%rhoijim,(lmn2_size,nspden))
         pawrhoij(irhoij)%rhoijim=zero
       end if
     end if

   end do
 end if

 if (present(pawtab)) then
   LIBPAW_POINTER_DEALLOCATE(lmn_size)
 end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

end subroutine pawrhoij_alloc
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/pawrhoij_free
!! NAME
!! pawrhoij_free
!!
!! FUNCTION
!! Destroy a pawrhoij datastructure
!!
!! SIDE EFFECTS
!! pawrhoij(:)<type(pawrhoij_type)>= rhoij datastructure
!!
!! PARENTS
!!      bethe_salpeter,d2frnl,dfpt_looppert,dfpt_nstpaw,dfpt_rhofermi
!!      dfpt_scfcv,dfpt_vtorho,energy,gstate,m_electronpositron,m_hdr,m_ioarr
!!      m_paral_pert,m_pawrhoij,m_scf_history,mrgscr,pawgrnl,pawmkrho
!!      pawmkrhoij,pawprt,posdoppler,respfn,screening,setup_bse,setup_positron
!!      setup_screening,setup_sigma,sigma,vtorho,wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawrhoij_free(pawrhoij)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawrhoij_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(pawrhoij_type),intent(inout) :: pawrhoij(:)

!Local variables-------------------------------
!scalars
 integer :: irhoij,nrhoij

! *************************************************************************

 nrhoij=size(pawrhoij)

 if (nrhoij>0) then
   do irhoij=1,nrhoij
     pawrhoij(irhoij)%nrhoijsel=0
     pawrhoij(irhoij)%ngrhoij=0
     pawrhoij(irhoij)%lmnmix_sz=0
     pawrhoij(irhoij)%use_rhoij_=0
     pawrhoij(irhoij)%use_rhoijp=0
     pawrhoij(irhoij)%use_rhoijres=0
     pawrhoij(irhoij)%use_rhoijim=0
     if (allocated(pawrhoij(irhoij)%rhoijp))       then
       LIBPAW_DEALLOCATE(pawrhoij(irhoij)%rhoijp)
     end if
     if (allocated(pawrhoij(irhoij)%rhoijselect))  then
       LIBPAW_DEALLOCATE(pawrhoij(irhoij)%rhoijselect)
     end if
     if (allocated(pawrhoij(irhoij)%grhoij))       then
       LIBPAW_DEALLOCATE(pawrhoij(irhoij)%grhoij)
     end if
     if (allocated(pawrhoij(irhoij)%kpawmix))      then
       LIBPAW_DEALLOCATE(pawrhoij(irhoij)%kpawmix)
     end if
     if (allocated(pawrhoij(irhoij)%rhoij_))       then
       LIBPAW_DEALLOCATE(pawrhoij(irhoij)%rhoij_)
     end if
     if (allocated(pawrhoij(irhoij)%rhoijres))     then
       LIBPAW_DEALLOCATE(pawrhoij(irhoij)%rhoijres)
     end if
     if (allocated(pawrhoij(irhoij)%rhoijim))     then
       LIBPAW_DEALLOCATE(pawrhoij(irhoij)%rhoijim)
     end if
   end do
 end if

end subroutine pawrhoij_free
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/pawrhoij_nullify
!! NAME
!! pawrhoij_nullify
!!
!! FUNCTION
!! Nullify (initialize to null) a pawrhoij datastructure
!!
!! SIDE EFFECTS
!! pawrhoij(:)<type(pawrhoij_type)>= rhoij datastructure
!!
!! PARENTS
!!      d2frnl,dfpt_looppert,m_ioarr,m_pawrhoij,m_scf_history,outscfcv,pawgrnl
!!      pawmkrho,pawprt,posdoppler,respfn
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawrhoij_nullify(pawrhoij)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawrhoij_nullify'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(pawrhoij_type),intent(inout) :: pawrhoij(:)

!Local variables-------------------------------
!scalars
 integer :: irhoij,nrhoij

! *************************************************************************

 ! MGPAW: This one could be removed/renamed,
 ! variables can be initialized in the datatype declaration
 ! Do we need to expose this in the public API?

 nrhoij=size(pawrhoij)

 if (nrhoij>0) then
   do irhoij=1,nrhoij
     pawrhoij(irhoij)%nrhoijsel=0
     pawrhoij(irhoij)%ngrhoij=0
     pawrhoij(irhoij)%lmnmix_sz=0
     pawrhoij(irhoij)%use_rhoij_=0
     pawrhoij(irhoij)%use_rhoijp=0
     pawrhoij(irhoij)%use_rhoijres=0
     pawrhoij(irhoij)%use_rhoijim=0
   end do
 end if

end subroutine pawrhoij_nullify
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/pawrhoij_copy
!! NAME
!! pawrhoij_copy
!!
!! FUNCTION
!! Copy one pawrhoij datastructure into another
!! Can take into accound changes of dimensions
!! Can copy a shared pawrhoij into distributed ones (when parallelism is activated)
!!
!! INPUTS
!!  keep_cplex= optional argument (logical, default=.TRUE.)
!!              if .TRUE. pawrhoij_out(:)%cplex is NOT MODIFIED,
!!              even if different from pawrhoij_in(:)%cplex
!!  keep_itypat= optional argument (logical, default=.FALSE.)
!!               if .TRUE. pawrhoij_out(:)%ityp is NOT MODIFIED,
!!               even if different from pawrhoij_in(:)%ityp
!!  keep_nspden= optional argument (logical, default=.TRUE.)
!!               if .TRUE. pawrhoij_out(:)%nspden is NOT MODIFIED,
!!               even if different from pawrhoij_in(:)%nspden
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  pawrhoij_in(:)<type(pawrhoij_type)>= input rhoij datastructure
!!
!! SIDE EFFECTS
!!  pawrhoij_out(:)<type(pawrhoij_type)>= output rhoij datastructure
!!
!! NOTES
!!  In case of a single copy operation pawrhoij_out must have been allocated.
!!
!! PARENTS
!!      bethe_salpeter,dfpt_looppert,gstate,inwffil,m_electronpositron,m_hdr
!!      m_ioarr,m_pawrhoij,m_wfk,outscfcv,pawmkrho,respfn,screening,setup_bse
!!      setup_positron,setup_screening,setup_sigma,sigma,wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawrhoij_copy(pawrhoij_in,pawrhoij_cpy, &
&          keep_cplex,keep_itypat,keep_nspden,& ! optional arguments
&          mpi_atmtab,comm_atom)            ! optional arguments (parallelism)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawrhoij_copy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: comm_atom
 logical,intent(in),optional :: keep_cplex,keep_itypat,keep_nspden
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 type(pawrhoij_type),intent(in) :: pawrhoij_in(:)
 type(pawrhoij_type),intent(inout),target :: pawrhoij_cpy(:)

!Local variables-------------------------------
!scalars
 integer :: cplex_in,cplex_out,dplex,dplex_in,dplex_out,i_in,i_out,ilmn
 integer :: irhoij,ispden,jrhoij,lmn2_size_out,lmnmix,my_comm_atom,my_nrhoij
 integer :: ngrhoij,nrhoij_in,nrhoij_max,nrhoij_out,nselect,nselect_out
 integer :: nspden_in,nspden_out,paral_case,use_rhoij_,use_rhoijp,use_rhoijres,use_rhoijim
 logical :: change_dim,keep_cplex_,keep_itypat_,keep_nspden_,my_atmtab_allocated,paral_atom
 character(len=500) :: msg
!arrays
 integer,pointer :: my_atmtab(:)
 integer,allocatable :: nlmn(:),typat(:)
 type(pawrhoij_type),pointer :: pawrhoij_out(:)

! *************************************************************************

!Retrieve sizes
 nrhoij_in=size(pawrhoij_in);nrhoij_out=size(pawrhoij_cpy)

!Init flags
 keep_cplex_=.true.
 if (present(keep_cplex)) keep_cplex_=keep_cplex
 keep_itypat_=.false.
 if (present(keep_itypat)) keep_itypat_=keep_itypat
 keep_nspden_=.true.
 if (present(keep_nspden)) keep_nspden_=keep_nspden

!Set up parallelism over atoms
 paral_atom=(present(comm_atom));if (paral_atom) paral_atom=(xmpi_comm_size(comm_atom)>1)
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 my_atmtab_allocated=.false.

!Determine in which case we are (parallelism, ...)
!No parallelism: a single copy operation
 paral_case=0;nrhoij_max=nrhoij_in
 pawrhoij_out => pawrhoij_cpy
 if (paral_atom) then
   if (nrhoij_out<nrhoij_in) then ! Parallelism: the copy operation is a scatter
     call get_my_natom(my_comm_atom,my_nrhoij,nrhoij_in)
     if (my_nrhoij==nrhoij_out) then
       call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,nrhoij_in)
       paral_case=1;nrhoij_max=nrhoij_out
       pawrhoij_out => pawrhoij_cpy
     else
       msg=' nrhoij_out should be equal to my_natom !'
       MSG_BUG(msg)
     end if
   else                            ! Parallelism: the copy operation is a gather
     call get_my_natom(my_comm_atom,my_nrhoij,nrhoij_out)
     if (my_nrhoij==nrhoij_in) then
       paral_case=2;nrhoij_max=nrhoij_in
       LIBPAW_DATATYPE_ALLOCATE(pawrhoij_out,(nrhoij_in))
       call pawrhoij_nullify(pawrhoij_out)
       if (nrhoij_in>0) then
         LIBPAW_ALLOCATE(typat,(nrhoij_in))
         LIBPAW_ALLOCATE(nlmn,(nrhoij_in))
         do irhoij=1,nrhoij_in
           typat(irhoij)=irhoij;nlmn(irhoij)=pawrhoij_in(irhoij)%lmn_size
         end do
         call pawrhoij_alloc(pawrhoij_out,pawrhoij_cpy(1)%cplex,pawrhoij_cpy(1)%nspden,&
&         pawrhoij_cpy(1)%nspinor,pawrhoij_cpy(1)%nsppol,typat,&
&         lmnsize=nlmn,ngrhoij=pawrhoij_cpy(1)%ngrhoij,nlmnmix=pawrhoij_cpy(1)%lmnmix_sz,&
&         use_rhoij_=pawrhoij_cpy(1)%use_rhoij_,&
&         use_rhoijp=pawrhoij_cpy(1)%use_rhoijp,&
&         use_rhoijres=pawrhoij_cpy(1)%use_rhoijres,&
&         use_rhoijim=pawrhoij_cpy(1)%use_rhoijim)
         LIBPAW_DEALLOCATE(typat)
         LIBPAW_DEALLOCATE(nlmn)
       end if
     else
       msg=' nrhoij_in should be equal to my_natom !'
       MSG_BUG(msg)
     end if
   end if
 end if

!Loop on rhoij components
 if (nrhoij_max>0) then
   do irhoij=1,nrhoij_max
     jrhoij=irhoij;if (paral_case==1) jrhoij=my_atmtab(irhoij)

     lmn2_size_out=pawrhoij_in(jrhoij)%lmn2_size
     cplex_in=pawrhoij_in(jrhoij)%cplex
     cplex_out=cplex_in;if(keep_cplex_)cplex_out=pawrhoij_out(irhoij)%cplex
     nspden_in=pawrhoij_in(jrhoij)%nspden
     nselect=pawrhoij_in(jrhoij)%nrhoijsel
     nselect_out=pawrhoij_out(irhoij)%nrhoijsel
     nspden_out=nspden_in;if(keep_nspden_)nspden_out=pawrhoij_out(irhoij)%nspden

     change_dim=(pawrhoij_out(irhoij)%cplex/=cplex_out.or. &
&     pawrhoij_out(irhoij)%lmn2_size/=lmn2_size_out.or. &
&     pawrhoij_out(irhoij)%nspden/=nspden_out.or. &
&     nselect/=nselect_out)
     dplex_in=cplex_in-1;dplex_out=cplex_out-1
     dplex=min(dplex_in,dplex_out)

!    Scalars
     pawrhoij_out(irhoij)%cplex=cplex_out+0
     pawrhoij_out(irhoij)%nspden=nspden_out+0
     pawrhoij_out(irhoij)%lmn2_size=lmn2_size_out+0
     pawrhoij_out(irhoij)%lmn_size=pawrhoij_in(jrhoij)%lmn_size+0
     if(.not.keep_itypat_) pawrhoij_out(irhoij)%itypat =pawrhoij_in(jrhoij)%itypat+0
     if(.not.keep_nspden_) pawrhoij_out(irhoij)%nsppol =pawrhoij_in(jrhoij)%nsppol+0
     if(.not.keep_nspden_) pawrhoij_out(irhoij)%nspinor=pawrhoij_in(jrhoij)%nspinor+0
     pawrhoij_out(irhoij)%nrhoijsel=nselect+0
!    if (pawrhoij_out(irhoij)%itypat/=pawrhoij_in(jrhoij)%itypat) then
!    write(unit=msg,fmt='(a,i3,a)') 'Type of atom ',jrhoij,' is different (dont copy it) !'
!    MSG_COMMENT(msg)
!    end if

!    Optional pointer: non-zero elements of rhoij
     use_rhoijp=pawrhoij_in(jrhoij)%use_rhoijp
     if (pawrhoij_out(irhoij)%use_rhoijp/=use_rhoijp) then
       if (pawrhoij_out(irhoij)%use_rhoijp>0)  then
         LIBPAW_DEALLOCATE(pawrhoij_out(irhoij)%rhoijp)
         LIBPAW_DEALLOCATE(pawrhoij_out(irhoij)%rhoijselect)
       end if
       if (use_rhoijp>0)  then
         LIBPAW_ALLOCATE(pawrhoij_out(irhoij)%rhoijp,(cplex_out*lmn2_size_out,nspden_out))
         LIBPAW_ALLOCATE(pawrhoij_out(irhoij)%rhoijselect,(lmn2_size_out))
       end if
       pawrhoij_out(irhoij)%use_rhoijp=use_rhoijp
     end if
     if (use_rhoijp>0) then
       if (change_dim) then
         if(allocated(pawrhoij_out(irhoij)%rhoijp)) then
           LIBPAW_DEALLOCATE(pawrhoij_out(irhoij)%rhoijp)
         end if
         LIBPAW_ALLOCATE(pawrhoij_out(irhoij)%rhoijp,(cplex_out*lmn2_size_out,nspden_out))
       end if
       if (cplex_out==cplex_in.and.nspden_out==nspden_in) then
         do ispden=1,nspden_out
           pawrhoij_out(irhoij)%rhoijp(1:cplex_out*nselect,ispden)=pawrhoij_in(jrhoij)%rhoijp(1:cplex_out*nselect,ispden)+zero
           if ( (nselect<lmn2_size_out) .and. &
&               (size(pawrhoij_out(irhoij)%rhoijp(:,ispden))==cplex_out*lmn2_size_out) ) then
             pawrhoij_out(irhoij)%rhoijp(cplex_out*nselect+1:cplex_out*lmn2_size_out,ispden)=zero
           end if
         end do
       else
         pawrhoij_out(irhoij)%rhoijp(:,:)=zero
         if (nspden_out==1) then
           if (nspden_in==2) then
             do ilmn=1,nselect
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijp(i_out:i_out+dplex,1)=pawrhoij_in(jrhoij)%rhoijp(i_in:i_in+dplex,1) &
&                                                              +pawrhoij_in(jrhoij)%rhoijp(i_in:i_in+dplex,2)+zero
             end do
           else ! nspden_in==1 or nspden_in=4
             do ilmn=1,nselect
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijp(i_out:i_out+dplex,1)=pawrhoij_in(jrhoij)%rhoijp(i_in:i_in+dplex,1)+zero
             end do
           end if
         else if (nspden_out==2) then
           if (nspden_in==1) then
             do ilmn=1,nselect
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijp(i_out:i_out+dplex,1)=half*pawrhoij_in(jrhoij)%rhoijp(i_in:i_in+dplex,1)+zero
               pawrhoij_out(irhoij)%rhoijp(i_out:i_out+dplex,2)=half*pawrhoij_in(jrhoij)%rhoijp(i_in:i_in+dplex,1)+zero
             end do
           else if (nspden_in==2) then
             do ilmn=1,nselect
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijp(i_out:i_out+dplex,1:2)=pawrhoij_in(jrhoij)%rhoijp(i_in:i_in+dplex,1:2)+zero
             end do
           else ! nspden_in==4
             do ilmn=1,nselect
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijp(i_out:i_out+dplex,1)=half*(pawrhoij_in(jrhoij)%rhoijp(i_in:i_in+dplex,1) &
&                                                                    +pawrhoij_in(jrhoij)%rhoijp(i_in:i_in+dplex,4))+zero
               pawrhoij_out(irhoij)%rhoijp(i_out:i_out+dplex,2)=half*(pawrhoij_in(jrhoij)%rhoijp(i_in:i_in+dplex,1) &
&                                                                    -pawrhoij_in(jrhoij)%rhoijp(i_in:i_in+dplex,4))+zero
             end do
           end if
         else if (nspden_out==4) then
           if (nspden_in==1) then
             do ilmn=1,nselect
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijp(i_out:i_out+dplex,1)=pawrhoij_in(jrhoij)%rhoijp(i_in:i_in+dplex,1)+zero
               pawrhoij_out(irhoij)%rhoijp(i_out:i_out+dplex,2:4)=zero
             end do
           else if (nspden_in==2) then
             do ilmn=1,nselect
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijp(i_out:i_out+dplex,1)=pawrhoij_in(jrhoij)%rhoijp(i_in:i_in+dplex,1) &
&                                                              +pawrhoij_in(jrhoij)%rhoijp(i_in:i_in+dplex,2)+zero
               pawrhoij_out(irhoij)%rhoijp(i_out:i_out+dplex,4)=pawrhoij_in(jrhoij)%rhoijp(i_in:i_in+dplex,1) &
&                                                              -pawrhoij_in(jrhoij)%rhoijp(i_in:i_in+dplex,2)+zero
               pawrhoij_out(irhoij)%rhoijp(i_out:i_out+dplex,2:3)=zero
             end do
           else ! nspden_in==4
             do ilmn=1,nselect
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijp(i_out:i_out+dplex,1:4)=pawrhoij_in(jrhoij)%rhoijp(i_in:i_in+dplex,1:4)+zero
             end do
           end if
         end if
       end if
     end if

!    Optional pointer: indexes for non-zero elements selection
     if (use_rhoijp>0) then
       if (change_dim) then
         if(allocated(pawrhoij_out(irhoij)%rhoijselect)) then
           LIBPAW_DEALLOCATE(pawrhoij_out(irhoij)%rhoijselect)
         end if
         LIBPAW_ALLOCATE(pawrhoij_out(irhoij)%rhoijselect,(lmn2_size_out))
       end if
       pawrhoij_out(irhoij)%rhoijselect(1:nselect)=pawrhoij_in(jrhoij)%rhoijselect(1:nselect)+0
       if ( (nselect<lmn2_size_out) .and. &
&           (size(pawrhoij_out(irhoij)%rhoijselect(:))==lmn2_size_out) ) then
         pawrhoij_out(irhoij)%rhoijselect(nselect+1:lmn2_size_out)=0
       end if
     end if

!    Optional pointer: indexes of rhoij to be mixed
     lmnmix=pawrhoij_in(jrhoij)%lmnmix_sz
     if (pawrhoij_out(irhoij)%lmnmix_sz/=lmnmix) then
       if (pawrhoij_out(irhoij)%lmnmix_sz>0)  then
         LIBPAW_DEALLOCATE(pawrhoij_out(irhoij)%kpawmix)
       end if
       if (lmnmix>0)  then
         LIBPAW_ALLOCATE(pawrhoij_out(irhoij)%kpawmix,(lmnmix))
       end if
       pawrhoij_out(irhoij)%lmnmix_sz=lmnmix
     end if
     if (lmnmix>0) pawrhoij_out(irhoij)%kpawmix(1:lmnmix)=pawrhoij_in(jrhoij)%kpawmix(1:lmnmix)

!    Optional pointer: gradients of rhoij
     ngrhoij=pawrhoij_in(jrhoij)%ngrhoij
     if (pawrhoij_out(irhoij)%ngrhoij/=ngrhoij) then
       if (pawrhoij_out(irhoij)%ngrhoij>0)  then
         LIBPAW_DEALLOCATE(pawrhoij_out(irhoij)%grhoij)
       end if
       if (ngrhoij>0)  then
         LIBPAW_ALLOCATE(pawrhoij_out(irhoij)%grhoij,(ngrhoij,cplex_out*lmn2_size_out,nspden_out))
       end if
       pawrhoij_out(irhoij)%ngrhoij=ngrhoij
     end if
     if (ngrhoij>0) then
       if (change_dim) then
         LIBPAW_DEALLOCATE(pawrhoij_out(irhoij)%grhoij)
         LIBPAW_ALLOCATE(pawrhoij_out(irhoij)%grhoij,(ngrhoij,cplex_out*lmn2_size_out,nspden_out))
       end if
       if (cplex_out==cplex_in.and.nspden_out==nspden_in) then
         do ispden=1,nspden_out
           do ilmn=1,cplex_out*lmn2_size_out
             pawrhoij_out(irhoij)%grhoij(1:ngrhoij,ilmn,ispden)= &
&               pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,ilmn,ispden)+zero
           end do
         end do
       else
         pawrhoij_out(irhoij)%grhoij(:,:,:)=zero
         if (nspden_out==1) then
           if (nspden_in==2) then
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out:i_out+dplex,1)= &
&                 pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in:i_in+dplex,1) &
&                +pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in:i_in+dplex,2)+zero
             end do
           else ! nspden_in==4
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out:i_out+dplex,1)= &
&                 pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in:i_in+dplex,1)+zero
             end do
           end if
         else if (nspden_out==2) then
           if (nspden_in==1) then
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out:i_out+dplex,1)= &
&                 half*pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in:i_in+dplex,1)+zero
               pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out:i_out+dplex,2)= &
&                 half*pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in:i_in+dplex,1)+zero
             end do
           else if (nspden_in==2) then
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out:i_out+dplex,1:2)= &
&                 pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in:i_in+dplex,1:2)+zero
             end do
           else ! nspden_in==4
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out:i_out+dplex,1)= &
&                 half*(pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in:i_in+dplex,1) &
&                      +pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in:i_in+dplex,4))+zero
               pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out:i_out+dplex,2)= &
&                 half*(pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in:i_in+dplex,1) &
&                      -pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in:i_in+dplex,4))+zero
             end do
           end if
         else if (nspden_out==4) then
           if (nspden_in==1) then
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out:i_out+dplex,1)= &
&                 pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in:i_in+dplex,1)+zero
               pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out:i_out+dplex,2:4)=zero
             end do
           else if (nspden_in==2) then
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out:i_out+dplex,1)= &
&                 pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in:i_in+dplex,1) &
&                +pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in:i_in+dplex,2)+zero
               pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out:i_out+dplex,4)= &
&                 pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in:i_in+dplex,1) &
&                -pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in:i_in+dplex,2)+zero
               pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out:i_out+dplex,2:3)=zero
             end do
           else ! nspden_in==4
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out:i_out+dplex,1:4)= &
&                 pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in:i_in+dplex,1:4)+zero
             end do
           end if
         end if
       end if
     end if

!    Optional pointer: residuals of rhoij
     use_rhoijres=pawrhoij_in(jrhoij)%use_rhoijres
     if (pawrhoij_out(irhoij)%use_rhoijres/=use_rhoijres) then
       if (pawrhoij_out(irhoij)%use_rhoijres>0)  then
         LIBPAW_DEALLOCATE(pawrhoij_out(irhoij)%rhoijres)
       end if
       if (use_rhoijres>0)  then
         LIBPAW_ALLOCATE(pawrhoij_out(irhoij)%rhoijres,(cplex_out*lmn2_size_out,nspden_out))
       end if
       pawrhoij_out(irhoij)%use_rhoijres=use_rhoijres
     end if
     if (use_rhoijres>0) then
       if (change_dim) then
         LIBPAW_DEALLOCATE(pawrhoij_out(irhoij)%rhoijres)
         LIBPAW_ALLOCATE(pawrhoij_out(irhoij)%rhoijres,(cplex_out*lmn2_size_out,nspden_out))
       end if
       if (cplex_out==cplex_in.and.nspden_out==nspden_in) then
         do ispden=1,nspden_out
           do ilmn=1,cplex_out*lmn2_size_out
             pawrhoij_out(irhoij)%rhoijres(ilmn,ispden)=pawrhoij_in(jrhoij)%rhoijres(ilmn,ispden)+zero
           end do
         end do
       else
         pawrhoij_out(irhoij)%rhoijres(:,:)=zero
         if (nspden_out==1) then
           if (nspden_in==2) then
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijres(i_out:i_out+dplex,1)=pawrhoij_in(jrhoij)%rhoijres(i_in:i_in+dplex,1) &
&                                                                +pawrhoij_in(jrhoij)%rhoijres(i_in:i_in+dplex,2)+zero
             end do
           else ! nspden_in==4
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijres(i_out:i_out+dplex,1)=pawrhoij_in(jrhoij)%rhoijres(i_in:i_in+dplex,1)+zero
             end do
           end if
         else if (nspden_out==2) then
           if (nspden_in==1) then
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijres(i_out:i_out+dplex,1)=half*pawrhoij_in(jrhoij)%rhoijres(i_in:i_in+dplex,1)+zero
               pawrhoij_out(irhoij)%rhoijres(i_out:i_out+dplex,2)=half*pawrhoij_in(jrhoij)%rhoijres(i_in:i_in+dplex,1)+zero
             end do
           else if (nspden_in==2) then
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijres(i_out:i_out+dplex,1:2)=pawrhoij_in(jrhoij)%rhoijres(i_in:i_in+dplex,1:2)+zero
             end do
           else ! nspden_in==4
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijres(i_out:i_out+dplex,1)=half*(pawrhoij_in(jrhoij)%rhoijres(i_in:i_in+dplex,1) &
&                                                                      +pawrhoij_in(jrhoij)%rhoijres(i_in:i_in+dplex,4))+zero
               pawrhoij_out(irhoij)%rhoijres(i_out:i_out+dplex,2)=half*(pawrhoij_in(jrhoij)%rhoijres(i_in:i_in+dplex,1) &
&                                                                      -pawrhoij_in(jrhoij)%rhoijres(i_in:i_in+dplex,4))+zero
             end do
           end if
         else if (nspden_out==4) then
           if (nspden_in==1) then
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijres(i_out:i_out+dplex,1)=pawrhoij_in(jrhoij)%rhoijres(i_in:i_in+dplex,1)+zero
               pawrhoij_out(irhoij)%rhoijres(i_out:i_out+dplex,2:4)=zero
             end do
           else if (nspden_in==2) then
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijres(i_out:i_out+dplex,1)=pawrhoij_in(jrhoij)%rhoijres(i_in:i_in+dplex,1) &
&                                                                +pawrhoij_in(jrhoij)%rhoijres(i_in:i_in+dplex,2)+zero
               pawrhoij_out(irhoij)%rhoijres(i_out:i_out+dplex,4)=pawrhoij_in(jrhoij)%rhoijres(i_in:i_in+dplex,1) &
&                                                                -pawrhoij_in(jrhoij)%rhoijres(i_in:i_in+dplex,2)+zero
               pawrhoij_out(irhoij)%rhoijres(i_out:i_out+dplex,2:3)=zero
             end do
           else ! nspden_in==4
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijres(i_out:i_out+dplex,1:4)=pawrhoij_in(jrhoij)%rhoijres(i_in:i_in+dplex,1:4)+zero
             end do
           end if
         end if
       end if
     end if

!    Optional pointer: imaginary part of rhoij
     use_rhoijim=pawrhoij_in(jrhoij)%use_rhoijim
     if (pawrhoij_out(irhoij)%use_rhoijim/=use_rhoijim) then
       if (pawrhoij_out(irhoij)%use_rhoijim>0)  then
         LIBPAW_DEALLOCATE(pawrhoij_out(irhoij)%rhoijim)
       end if
       if (use_rhoijim>0)  then
         LIBPAW_ALLOCATE(pawrhoij_out(irhoij)%rhoijim,(lmn2_size_out,nspden_out))
       end if
       pawrhoij_out(irhoij)%use_rhoijim=use_rhoijim
     end if
     if (use_rhoijim>0) then
       if (change_dim) then
         LIBPAW_DEALLOCATE(pawrhoij_out(irhoij)%rhoijim)
         LIBPAW_ALLOCATE(pawrhoij_out(irhoij)%rhoijim,(lmn2_size_out,nspden_out))
       end if
       if (nspden_out==nspden_in) then
         do ispden=1,nspden_out
           do ilmn=1,lmn2_size_out
             pawrhoij_out(irhoij)%rhoijim(ilmn,ispden)=pawrhoij_in(jrhoij)%rhoijim(ilmn,ispden)+zero
           end do
         end do
       else
         pawrhoij_out(irhoij)%rhoijim(:,:)=zero
         if (nspden_out==1) then
           if (nspden_in==2) then
             do ilmn=1,lmn2_size_out
               i_in=ilmn;i_out=ilmn
               pawrhoij_out(irhoij)%rhoijim(i_out,1)=pawrhoij_in(jrhoij)%rhoijim(i_in,1) &
&                                                                +pawrhoij_in(jrhoij)%rhoijim(i_in,2)+zero
             end do
           else ! nspden_in==4
             do ilmn=1,lmn2_size_out
               i_in=ilmn;i_out=ilmn
               pawrhoij_out(irhoij)%rhoijim(i_out,1)=pawrhoij_in(jrhoij)%rhoijim(i_in,1)+zero
             end do
           end if
         else if (nspden_out==2) then
           if (nspden_in==1) then
             do ilmn=1,lmn2_size_out
               i_in=ilmn;i_out=ilmn
               pawrhoij_out(irhoij)%rhoijim(i_out,1)=half*pawrhoij_in(jrhoij)%rhoijim(i_in,1)+zero
               pawrhoij_out(irhoij)%rhoijim(i_out,2)=half*pawrhoij_in(jrhoij)%rhoijim(i_in,1)+zero
             end do
           else if (nspden_in==2) then
             do ilmn=1,lmn2_size_out
               i_in=ilmn;i_out=ilmn
               pawrhoij_out(irhoij)%rhoijim(i_out,1:2)=pawrhoij_in(jrhoij)%rhoijim(i_in,1:2)+zero
             end do
           else ! nspden_in==4
             do ilmn=1,lmn2_size_out
               i_in=ilmn;i_out=ilmn
               pawrhoij_out(irhoij)%rhoijim(i_out,1)=half*(pawrhoij_in(jrhoij)%rhoijim(i_in,1) &
&                                                                      +pawrhoij_in(jrhoij)%rhoijim(i_in,4))+zero
               pawrhoij_out(irhoij)%rhoijim(i_out,2)=half*(pawrhoij_in(jrhoij)%rhoijim(i_in,1) &
&                                                                      -pawrhoij_in(jrhoij)%rhoijim(i_in,4))+zero
             end do
           end if
         else if (nspden_out==4) then
           if (nspden_in==1) then
             do ilmn=1,lmn2_size_out
               i_in=ilmn;i_out=ilmn
               pawrhoij_out(irhoij)%rhoijim(i_out,1)=pawrhoij_in(jrhoij)%rhoijim(i_in,1)+zero
               pawrhoij_out(irhoij)%rhoijim(i_out,2:4)=zero
             end do
           else if (nspden_in==2) then
             do ilmn=1,lmn2_size_out
               i_in=ilmn;i_out=ilmn
               pawrhoij_out(irhoij)%rhoijim(i_out,1)=pawrhoij_in(jrhoij)%rhoijim(i_in,1) &
&                                                                +pawrhoij_in(jrhoij)%rhoijim(i_in,2)+zero
               pawrhoij_out(irhoij)%rhoijim(i_out,4)=pawrhoij_in(jrhoij)%rhoijim(i_in,1) &
&                                                                -pawrhoij_in(jrhoij)%rhoijim(i_in,2)+zero
               pawrhoij_out(irhoij)%rhoijim(i_out,2:3)=zero
             end do
           else ! nspden_in==4
             do ilmn=1,lmn2_size_out
               i_in=ilmn;i_out=ilmn
               pawrhoij_out(irhoij)%rhoijim(i_out,1:4)=pawrhoij_in(jrhoij)%rhoijim(i_in,1:4)+zero
             end do
           end if
         end if
       end if
     end if

!    Optional pointer: non-symmetrized rhoij
     use_rhoij_=pawrhoij_in(jrhoij)%use_rhoij_
     if (pawrhoij_out(irhoij)%use_rhoij_/=use_rhoij_) then
       if (pawrhoij_out(irhoij)%use_rhoij_>0)  then
         LIBPAW_DEALLOCATE(pawrhoij_out(irhoij)%rhoij_)
       end if
       if (use_rhoij_>0)  then
         LIBPAW_ALLOCATE(pawrhoij_out(irhoij)%rhoij_,(cplex_out*lmn2_size_out,nspden_out))
       end if
       pawrhoij_out(irhoij)%use_rhoij_=use_rhoij_
     end if
     if (use_rhoij_>0) then
       if (change_dim) then
         if(allocated(pawrhoij_out(irhoij)%rhoij_)) then
           LIBPAW_DEALLOCATE(pawrhoij_out(irhoij)%rhoij_)
         end if
         LIBPAW_ALLOCATE(pawrhoij_out(irhoij)%rhoij_,(cplex_out*lmn2_size_out,nspden_out))
       end if
       if (cplex_out==cplex_in.and.nspden_out==nspden_in) then
         do ispden=1,nspden_out
           do ilmn=1,cplex_out*lmn2_size_out
             pawrhoij_out(irhoij)%rhoij_(ilmn,ispden)=pawrhoij_in(jrhoij)%rhoij_(ilmn,ispden)+zero
           end do
         end do
       else
         pawrhoij_out(irhoij)%rhoij_(:,:)=zero
         if (nspden_out==1) then
           if (nspden_in==2) then
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoij_(i_out:i_out+dplex,1)=pawrhoij_in(jrhoij)%rhoij_(i_in:i_in+dplex,1) &
&                                                              +pawrhoij_in(jrhoij)%rhoij_(i_in:i_in+dplex,2)+zero
             end do
           else ! nspden_in==1 or nspden_in==4
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoij_(i_out:i_out+dplex,1)=pawrhoij_in(jrhoij)%rhoij_(i_in:i_in+dplex,1)+zero
             end do
           end if
         else if (nspden_out==2) then
           if (nspden_in==1) then
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoij_(i_out:i_out+dplex,1)=half*pawrhoij_in(jrhoij)%rhoij_(i_in:i_in+dplex,1)+zero
               pawrhoij_out(irhoij)%rhoij_(i_out:i_out+dplex,2)=half*pawrhoij_in(irhoij)%rhoij_(i_in:i_in+dplex,1)+zero
             end do
           else if (nspden_in==2) then
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoij_(i_out:i_out+dplex,1:2)=pawrhoij_in(jrhoij)%rhoij_(i_in:i_in+dplex,1:2)+zero
             end do
           else ! nspden_in==4
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoij_(i_out:i_out+dplex,1)=half*(pawrhoij_in(jrhoij)%rhoij_(i_in:i_in+dplex,1) &
&                                                                    +pawrhoij_in(jrhoij)%rhoij_(i_in:i_in+dplex,4))+zero
               pawrhoij_out(irhoij)%rhoij_(i_out:i_out+dplex,2)=half*(pawrhoij_in(jrhoij)%rhoij_(i_in:i_in+dplex,1) &
&                                                                    -pawrhoij_in(jrhoij)%rhoij_(i_in:i_in+dplex,4))+zero
             end do
           end if
         else if (nspden_out==4) then
           if (nspden_in==1) then
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoij_(i_out:i_out+dplex,1)=pawrhoij_in(jrhoij)%rhoij_(i_in:i_in+dplex,1)+zero
               pawrhoij_out(irhoij)%rhoij_(i_out:i_out+dplex,2:4)=zero
             end do
           else if (nspden_in==2) then
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoij_(i_out:i_out+dplex,1)=pawrhoij_in(jrhoij)%rhoij_(i_in:i_in+dplex,1) &
&                                                              +pawrhoij_in(jrhoij)%rhoij_(i_in:i_in+dplex,2)+zero
               pawrhoij_out(irhoij)%rhoij_(i_out:i_out+dplex,4)=pawrhoij_in(jrhoij)%rhoij_(i_in:i_in+dplex,1) &
&                                                              -pawrhoij_in(jrhoij)%rhoij_(i_in:i_in+dplex,2)+zero
               pawrhoij_out(irhoij)%rhoij_(i_out:i_out+dplex,2:3)=zero
             end do
           else ! nspden_in==4
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoij_(i_out:i_out+dplex,1:4)=pawrhoij_in(jrhoij)%rhoij_(i_in:i_in+dplex,1:4)+zero
             end do
           end if
         end if
       end if
     end if

   end do ! irhoij
 end if

!Parallel case: do a gather if needed
 if (paral_case==2) then
   call pawrhoij_free(pawrhoij_cpy)
   call pawrhoij_gather(pawrhoij_out,pawrhoij_cpy,-1,my_comm_atom)
   call pawrhoij_free(pawrhoij_out)
   LIBPAW_DATATYPE_DEALLOCATE(pawrhoij_out)

!  Sequential case: fill missing elements
 else if (paral_case==0) then
   if (nrhoij_in<nrhoij_out) then
     do irhoij=nrhoij_in+1,nrhoij_out
       pawrhoij_cpy(irhoij)%nrhoijsel=0
       if (pawrhoij_cpy(irhoij)%use_rhoijp>0) then
         pawrhoij_cpy(irhoij)%rhoijselect=0
         pawrhoij_cpy(irhoij)%rhoijp=zero
       end if
       if (pawrhoij_cpy(irhoij)%lmnmix_sz>0) pawrhoij_cpy(irhoij)%kpawmix=0
       if (pawrhoij_cpy(irhoij)%ngrhoij>0) pawrhoij_cpy(irhoij)%grhoij=zero
       if (pawrhoij_cpy(irhoij)%use_rhoij_>0) pawrhoij_cpy(irhoij)%rhoij_=zero
       if (pawrhoij_cpy(irhoij)%use_rhoijres>0) pawrhoij_cpy(irhoij)%rhoijres=zero
       if (pawrhoij_cpy(irhoij)%use_rhoijim>0) pawrhoij_cpy(irhoij)%rhoijim=zero
     end do
   end if
 end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

end subroutine pawrhoij_copy
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/pawrhoij_gather
!! NAME
!! pawrhoij_gather
!!
!! FUNCTION
!!   (All)Gather pawrhoij datastructures
!!
!! INPUTS
!!  master=master communicator receiving data ; if -1 do a ALLGATHER
!!  comm_atom= communicator
!!  pawrhoij_in(:)<type(pawrhoij_type)>= input rhoij datastructures on every process
!!  with_grhoij   : optional argument (logical, default=.TRUE.)
!!                  TRUE if pawrhoij%grhoij field is included in the gather operation
!!  with_lmnmix   : optional argument (logical, default=.TRUE.)
!!                  TRUE if pawrhoij%lmnmix field is included in the gather operation
!!  with_rhoijp   : optional argument (logical, default=.TRUE.)
!!                  TRUE if pawrhoij%rhoijp and pawrhoij%rhoijselect fields
!!                       are included in the gather operation
!!  with_rhoijres : optional argument (logical, default=.TRUE.)
!!                 TRUE if pawrhoij%rhoijres field is included in the gather operation
!!  with_rhoijim : optional argument (logical, default=.TRUE.)
!!                 TRUE if pawrhoij%rhoijim field is included in the gather operation
!!  with_rhoij_   : optional argument (logical, default=.TRUE.)
!!                  TRUE if pawrhoij%rhoij_ field is included in the gather operation
!!
!! OUTPUT
!!  pawrhoij_gathered(:)<type(pawrhoij_type)>= output rhoij datastructure
!!
!! NOTES
!!  The gathered structure are ordered like in sequential mode.
!!
!! PARENTS
!!      d2frnl,m_pawrhoij,pawgrnl,pawprt,posdoppler
!!
!! CHILDREN
!!
!! SOURCE


 subroutine pawrhoij_gather(pawrhoij_in,pawrhoij_gathered,master,comm_atom, &
&    with_grhoij,with_lmnmix,with_rhoijp,with_rhoijres,with_rhoijim,with_rhoij_) ! optional arguments


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawrhoij_gather'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: master,comm_atom
 logical,intent(in),optional :: with_grhoij,with_lmnmix,with_rhoijp,with_rhoijres,with_rhoijim,with_rhoij_
!arrays
 type(pawrhoij_type),intent(in) :: pawrhoij_in(:)
 type(pawrhoij_type),intent(inout) :: pawrhoij_gathered(:)
!Local variables-------------------------------
!scalars
 integer :: buf_dp_size,buf_dp_size_all,buf_int_size,buf_int_size_all
 integer :: cplex,ierr,ii,indx_dp,indx_int,irhoij,isp,jrhoij,lmn2_size,lmnmix,me_atom
 integer :: ngrhoij,nproc_atom,nrhoij_in,nrhoij_in_sum,nrhoij_out,nselect,nspden
 integer :: rhoij_size2,use_rhoijp,use_rhoijres,use_rhoijim,use_rhoij_
 logical :: my_atmtab_allocated,paral_atom
 logical :: with_grhoij_,with_lmnmix_,with_rhoijp_,with_rhoijres_,with_rhoijim_,with_rhoij__
 character(len=500) :: msg
!arrays
 integer :: bufsz(2)
 integer,allocatable :: buf_int(:),buf_int_all(:)
 integer,allocatable :: count_dp(:),count_int(:),count_tot(:),displ_dp(:),displ_int(:)
 integer, pointer :: my_atmtab(:)
 real(dp),allocatable :: buf_dp(:),buf_dp_all(:)

! *************************************************************************

 nrhoij_in=size(pawrhoij_in);nrhoij_out=size(pawrhoij_gathered)

 nproc_atom=xmpi_comm_size(comm_atom)
 me_atom=xmpi_comm_rank(comm_atom)

 if (nproc_atom==1) then
   if (master==-1.or.me_atom==master) then
     call pawrhoij_copy(pawrhoij_in,pawrhoij_gathered,.false.,.false.,.false.)
   end if
   return
 end if

!Test on sizes
 nrhoij_in_sum=nrhoij_in
 call xmpi_sum(nrhoij_in_sum,comm_atom,ierr)
 if (master==-1) then
   if (nrhoij_out/=nrhoij_in_sum) then
     msg='Wrong sizes sum[nrhoij_ij]/=nrhoij_out !'
     MSG_BUG(msg)
   end if
 else
   if (me_atom==master.and.nrhoij_out/=nrhoij_in_sum) then
     msg='(2) pawrhoij_gathered wrongly allocated !'
     MSG_BUG(msg)
   end if
 end if

!Optional arguments
 with_grhoij_  =.true.;if (present(with_grhoij))  with_grhoij_  =with_grhoij
 with_lmnmix_  =.true.;if (present(with_lmnmix))  with_lmnmix_  =with_lmnmix
 with_rhoijp_  =.true.;if (present(with_rhoijp))  with_rhoijp_  =with_rhoijp
 with_rhoijres_=.true.;if (present(with_rhoijres))with_rhoijres_=with_rhoijres
 with_rhoijim_ =.true.;if (present(with_rhoijim)) with_rhoijim_ =with_rhoijim
 with_rhoij__  =.true.;if (present(with_rhoij_))  with_rhoij__  =with_rhoij_

!Retrieve table of atoms
 paral_atom=.true.;nullify(my_atmtab)
 call get_my_atmtab(comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,nrhoij_in_sum, &
& my_natom_ref=nrhoij_in)

!Compute sizes of buffers
 buf_int_size=0;buf_dp_size=0
 nselect=0;lmnmix=0;ngrhoij=0;rhoij_size2=0
 use_rhoijp=0;use_rhoijres=0;use_rhoijim=0;use_rhoij_=0
 do irhoij=1,nrhoij_in
   cplex    =pawrhoij_in(irhoij)%cplex
   lmn2_size=pawrhoij_in(irhoij)%lmn2_size
   nspden   =pawrhoij_in(irhoij)%nspden
   if (with_lmnmix_) lmnmix=pawrhoij_in(irhoij)%lmnmix_sz
   if (with_grhoij_) ngrhoij=pawrhoij_in(irhoij)%ngrhoij
   if (with_rhoijp_) use_rhoijp=pawrhoij_in(irhoij)%use_rhoijp
   if (with_rhoijres_)use_rhoijres=pawrhoij_in(irhoij)%use_rhoijres
   if (with_rhoijim_)use_rhoijim=pawrhoij_in(irhoij)%use_rhoijim
   if (with_rhoij__) use_rhoij_=pawrhoij_in(irhoij)%use_rhoij_
   buf_int_size=buf_int_size+16
   if (use_rhoijp>0) then
     nselect=pawrhoij_in(irhoij)%nrhoijsel
     buf_int_size=buf_int_size+nselect
     buf_dp_size=buf_dp_size + cplex*nselect*nspden
   end if
   if (lmnmix>0)       buf_int_size=buf_int_size+lmnmix
   if (ngrhoij>0)      buf_dp_size=buf_dp_size + cplex*lmn2_size*nspden*ngrhoij
   if (use_rhoijres>0) buf_dp_size=buf_dp_size + cplex*lmn2_size*nspden
   if (use_rhoijim>0)  buf_dp_size=buf_dp_size + lmn2_size*nspden
   if (use_rhoij_>0) then
     rhoij_size2=size(pawrhoij_in(irhoij)%rhoij_,dim=2)
     buf_dp_size=buf_dp_size + cplex*lmn2_size*rhoij_size2
   end if
 end do

!Fill input buffers
 LIBPAW_ALLOCATE(buf_int,(buf_int_size))
 LIBPAW_ALLOCATE(buf_dp ,(buf_dp_size))
 indx_int=1;indx_dp =1
 lmnmix=0;ngrhoij=0;nselect=0;rhoij_size2=0
 use_rhoijp=0;use_rhoijres=0;use_rhoijim=0;use_rhoij_=0
 do irhoij=1,nrhoij_in
   cplex    =pawrhoij_in(irhoij)%cplex
   lmn2_size=pawrhoij_in(irhoij)%lmn2_size
   nspden   =pawrhoij_in(irhoij)%nspden
   if (with_lmnmix_) lmnmix=pawrhoij_in(irhoij)%lmnmix_sz
   if (with_grhoij_) ngrhoij=pawrhoij_in(irhoij)%ngrhoij
   if (with_rhoijp_) use_rhoijp=pawrhoij_in(irhoij)%use_rhoijp
   if (use_rhoijp > 0) nselect=pawrhoij_in(irhoij)%nrhoijsel
   if (with_rhoijres_)use_rhoijres=pawrhoij_in(irhoij)%use_rhoijres
   if (with_rhoijim_) use_rhoijim =pawrhoij_in(irhoij)%use_rhoijim
   if (with_rhoij__) use_rhoij_  =pawrhoij_in(irhoij)%use_rhoij_
   if (use_rhoij_> 0) rhoij_size2 =size(pawrhoij_in(irhoij)%rhoij_,dim=2)
   buf_int(indx_int)=my_atmtab(irhoij)               ;indx_int=indx_int+1
   buf_int(indx_int)=cplex                           ;indx_int=indx_int+1
   buf_int(indx_int)=lmn2_size                       ;indx_int=indx_int+1
   buf_int(indx_int)=nspden                          ;indx_int=indx_int+1
   buf_int(indx_int)=nselect                         ;indx_int=indx_int+1
   buf_int(indx_int)=lmnmix                          ;indx_int=indx_int+1
   buf_int(indx_int)=ngrhoij                         ;indx_int=indx_int+1
   buf_int(indx_int)=use_rhoijp                      ;indx_int=indx_int+1
   buf_int(indx_int)=use_rhoijres                    ;indx_int=indx_int+1
   buf_int(indx_int)=use_rhoijim                     ;indx_int=indx_int+1
   buf_int(indx_int)=use_rhoij_                      ;indx_int=indx_int+1
   buf_int(indx_int)=rhoij_size2                     ;indx_int=indx_int+1
   buf_int(indx_int)=pawrhoij_in(irhoij)%itypat      ;indx_int=indx_int+1
   buf_int(indx_int)=pawrhoij_in(irhoij)%lmn_size    ;indx_int=indx_int+1
   buf_int(indx_int)=pawrhoij_in(irhoij)%nsppol      ;indx_int=indx_int+1
   buf_int(indx_int)=pawrhoij_in(irhoij)%nspinor     ;indx_int=indx_int+1
   if (use_rhoijp>0) then
     buf_int(indx_int:indx_int+nselect-1)=pawrhoij_in(irhoij)%rhoijselect(1:nselect)
     indx_int=indx_int+nselect
     do isp=1,nspden
       buf_dp(indx_dp:indx_dp+cplex*nselect-1)=pawrhoij_in(irhoij)%rhoijp(1:cplex*nselect,isp)
       indx_dp=indx_dp+cplex*nselect
     end do
   end if
   if (lmnmix>0) then
     buf_int(indx_int:indx_int+lmnmix-1)=pawrhoij_in(irhoij)%kpawmix(1:lmnmix)
     indx_int=indx_int+lmnmix
   end if
   if (ngrhoij>0) then
     do isp=1,nspden
       do ii=1,cplex*lmn2_size
         buf_dp(indx_dp:indx_dp+ngrhoij-1)=pawrhoij_in(irhoij)%grhoij(1:ngrhoij,ii,isp)
         indx_dp=indx_dp+ngrhoij
       end do
     end do
   end if
   if (use_rhoijres>0) then
     do isp=1,nspden
       buf_dp(indx_dp:indx_dp+cplex*lmn2_size-1)=pawrhoij_in(irhoij)%rhoijres(1:cplex*lmn2_size,isp)
       indx_dp=indx_dp+cplex*lmn2_size
     end do
   end if
   if (use_rhoijim>0) then
     do isp=1,nspden
       buf_dp(indx_dp:indx_dp+lmn2_size-1)=pawrhoij_in(irhoij)%rhoijim(1:lmn2_size,isp)
       indx_dp=indx_dp+lmn2_size
     end do
   end if
   if (use_rhoij_>0) then
     do isp=1,rhoij_size2
       buf_dp(indx_dp:indx_dp+cplex*lmn2_size-1)=pawrhoij_in(irhoij)%rhoij_(1:cplex*lmn2_size,isp)
       indx_dp=indx_dp+cplex*lmn2_size
     end do
   end if
 end do

!Check
 if ((indx_int-1/=buf_int_size).or.(indx_dp-1/=buf_dp_size)) then
   write(msg,*) 'Wrong buffer sizes: buf_int_size=',buf_int_size,' buf_dp_size=',buf_dp_size
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
 do ii=1,nproc_atom
   count_int(ii)=count_tot(2*ii-1)
   count_dp (ii)=count_tot(2*ii)
 end do
 displ_int(1)=0;displ_dp(1)=0
 do ii=2,nproc_atom
   displ_int(ii)=displ_int(ii-1)+count_int(ii-1)
   displ_dp (ii)=displ_dp (ii-1)+count_dp (ii-1)
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
   call pawrhoij_free(pawrhoij_gathered)
   do irhoij=1,nrhoij_out
     jrhoij      =buf_int_all(indx_int)    ;indx_int=indx_int+1
     cplex       =buf_int_all(indx_int)    ;indx_int=indx_int+1
     lmn2_size   =buf_int_all(indx_int)    ;indx_int=indx_int+1
     nspden      =buf_int_all(indx_int)    ;indx_int=indx_int+1
     nselect     =buf_int_all(indx_int)    ;indx_int=indx_int+1
     lmnmix      =buf_int_all(indx_int)    ;indx_int=indx_int+1
     ngrhoij     =buf_int_all(indx_int)    ;indx_int=indx_int+1
     use_rhoijp  =buf_int_all(indx_int)    ;indx_int=indx_int+1
     use_rhoijres=buf_int_all(indx_int)    ;indx_int=indx_int+1
     use_rhoijim =buf_int_all(indx_int)    ;indx_int=indx_int+1
     use_rhoij_  =buf_int_all(indx_int)    ;indx_int=indx_int+1
     rhoij_size2 =buf_int_all(indx_int)    ;indx_int=indx_int+1
     pawrhoij_gathered(jrhoij)%itypat=buf_int_all(indx_int)   ;indx_int=indx_int+1
     pawrhoij_gathered(jrhoij)%lmn_size=buf_int_all(indx_int) ;indx_int=indx_int+1
     pawrhoij_gathered(jrhoij)%nsppol=buf_int_all(indx_int)   ;indx_int=indx_int+1
     pawrhoij_gathered(jrhoij)%nspinor=buf_int_all(indx_int)  ;indx_int=indx_int+1
     pawrhoij_gathered(jrhoij)%cplex=cplex
     pawrhoij_gathered(jrhoij)%lmn2_size=lmn2_size
     pawrhoij_gathered(jrhoij)%nspden=nspden
     pawrhoij_gathered(jrhoij)%nrhoijsel=nselect
     pawrhoij_gathered(jrhoij)%lmnmix_sz=lmnmix
     pawrhoij_gathered(jrhoij)%ngrhoij=ngrhoij
     pawrhoij_gathered(jrhoij)%use_rhoijp=use_rhoijp
     pawrhoij_gathered(jrhoij)%use_rhoijres=use_rhoijres
     pawrhoij_gathered(jrhoij)%use_rhoijim =use_rhoijim
     pawrhoij_gathered(jrhoij)%use_rhoij_=use_rhoij_
     if (use_rhoijp>0) then
       LIBPAW_ALLOCATE(pawrhoij_gathered(jrhoij)%rhoijselect,(lmn2_size))
       pawrhoij_gathered(jrhoij)%rhoijselect(1:nselect)=buf_int_all(indx_int:indx_int+nselect-1)
       if (nselect < lmn2_size )pawrhoij_gathered(jrhoij)%rhoijselect(nselect+1:lmn2_size)=zero
       indx_int=indx_int+nselect
       LIBPAW_ALLOCATE(pawrhoij_gathered(jrhoij)%rhoijp,(cplex*lmn2_size,nspden))
       do isp=1,nspden
         pawrhoij_gathered(jrhoij)%rhoijp(1:cplex*nselect,isp)=buf_dp_all(indx_dp:indx_dp+cplex*nselect-1)
         if (nselect < lmn2_size )pawrhoij_gathered(jrhoij)%rhoijp(cplex*nselect+1:cplex*lmn2_size,isp)=zero
         indx_dp=indx_dp+cplex*nselect
       end do
     end if
     if (lmnmix>0) then
       LIBPAW_ALLOCATE(pawrhoij_gathered(jrhoij)%kpawmix,(lmnmix))
       pawrhoij_gathered(jrhoij)%kpawmix(1:lmnmix)=buf_int_all(indx_int:indx_int+lmnmix-1)
       indx_int=indx_int+lmnmix
     end if
     if (ngrhoij>0) then
       LIBPAW_ALLOCATE(pawrhoij_gathered(jrhoij)%grhoij,(ngrhoij,cplex*lmn2_size,nspden))
       do isp=1,nspden
         do ii=1,cplex*lmn2_size
           pawrhoij_gathered(jrhoij)%grhoij(1:ngrhoij,ii,isp)=buf_dp_all(indx_dp:indx_dp+ngrhoij-1)
           indx_dp=indx_dp+ngrhoij
         end do
       end do
     end if
     if (use_rhoijres>0) then
       LIBPAW_ALLOCATE(pawrhoij_gathered(jrhoij)%rhoijres,(cplex*lmn2_size,nspden))
       do isp=1,nspden
         pawrhoij_gathered(jrhoij)%rhoijres(1:cplex*lmn2_size,isp)=buf_dp_all(indx_dp:indx_dp+cplex*lmn2_size-1)
         indx_dp=indx_dp+cplex*lmn2_size
       end do
     end if
     if (use_rhoijim>0) then
       LIBPAW_ALLOCATE(pawrhoij_gathered(jrhoij)%rhoijim,(lmn2_size,nspden))
       do isp=1,nspden
         pawrhoij_gathered(jrhoij)%rhoijim(1:lmn2_size,isp)=buf_dp_all(indx_dp:indx_dp+lmn2_size-1)
         indx_dp=indx_dp+lmn2_size
       end do
     end if
     if (use_rhoij_>0) then
       LIBPAW_ALLOCATE(pawrhoij_gathered(jrhoij)%rhoij_,(cplex*lmn2_size,rhoij_size2))
       do isp=1,rhoij_size2
         pawrhoij_gathered(jrhoij)%rhoij_(1:cplex*lmn2_size,isp)=buf_dp_all(indx_dp:indx_dp+cplex*lmn2_size-1)
         indx_dp=indx_dp+cplex*lmn2_size
       end do
     end if
   end do
   if ((indx_int/=1+buf_int_size_all).or.(indx_dp/=1+buf_dp_size_all)) then
     write(msg,*) 'Wrong buffer sizes: buf_int_size_all=',buf_int_size_all,' buf_dp_size_all=',buf_dp_size_all
     MSG_BUG(msg)
   end if
 end if

!Free memory
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)
 LIBPAW_DEALLOCATE(buf_int)
 LIBPAW_DEALLOCATE(buf_dp)
 LIBPAW_DEALLOCATE(buf_int_all)
 LIBPAW_DEALLOCATE(buf_dp_all)

end subroutine pawrhoij_gather
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/pawrhoij_bcast
!! NAME
!! pawrhoij_bcast
!!
!! FUNCTION
!!   Broadcast pawrhoij datastructures
!!   Can take into account a distribution of data over a "atom" communicator
!!
!! INPUTS
!!  master=master communicator receiving data
!!  mpicomm= MPI communicator
!!  comm_atom= --optional-- MPI communicator over atoms
!!  pawrhoij_in(:)<type(pawrhoij_type)>= input rhoij datastructures on master process
!!
!! OUTPUT
!!  pawrhoij_out(:)<type(pawrhoij_type)>= output rhoij datastructure on every process
!!    Eventually distributed according to comm_atom communicator
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!
!! SOURCE

 subroutine pawrhoij_bcast(pawrhoij_in,pawrhoij_out,master,mpicomm,comm_atom)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawrhoij_bcast'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: master,mpicomm
 integer,intent(in),optional :: comm_atom
!arrays
 type(pawrhoij_type),intent(in) :: pawrhoij_in(:)
 type(pawrhoij_type),intent(inout) :: pawrhoij_out(:)
!Local variables-------------------------------
!scalars
 integer :: buf_dp_size,buf_dp_size_all,buf_int_size,buf_int_size_all
 integer :: cplex,ierr,ii,indx_dp,indx_int,iproc,irhoij,isp,jrhoij,lmn2_size,lmnmix,me,me_atom
 integer :: my_comm_atom,ngrhoij,nproc,nproc_atom,nrhoij_in,nrhoij_out,nrhoij_out_all
 integer :: nselect,nspden,rhoij_size2,use_rhoijp,use_rhoijres,use_rhoijim,use_rhoij_
 logical :: my_atmtab_allocated,paral_atom
 character(len=500) :: msg
!arrays
 integer :: buf_size(2)
 integer,allocatable :: atmtab(:),buf_int_size_i(:),buf_dp_size_i(:)
 integer,allocatable :: count_dp(:),count_int(:),disp_dp(:),disp_int(:),nrhoij_out_i(:)
 integer,allocatable,target :: buf_int(:)
 integer,pointer :: buf_int_all(:),my_atmtab(:)
 real(dp),allocatable,target :: buf_dp(:)
 real(dp),pointer :: buf_dp_all(:)

! *************************************************************************

!Load MPI "atom" distribution data
 my_comm_atom=xmpi_comm_self;nproc_atom=1;me_atom=0
 if (present(comm_atom)) then
   my_comm_atom=comm_atom
   me_atom=xmpi_comm_rank(my_comm_atom)
   nproc_atom=xmpi_comm_size(my_comm_atom)
   paral_atom=(nproc_atom>1)
   if (my_comm_atom/=mpicomm.and.nproc_atom/=1) then
     msg='wrong comm_atom communicator !'
     MSG_BUG(msg)
   end if
 end if

!Load global MPI data
 me=xmpi_comm_rank(mpicomm)
 nproc=xmpi_comm_size(mpicomm)

!Just copy in case of a sequential run
 if (nproc==1.and.nproc_atom==1) then
   call pawrhoij_copy(pawrhoij_in,pawrhoij_out,.false.,.false.,.false.)
   return
 end if

!Retrieve and test pawrhoij sizes
 nrhoij_in=0;if (me==master) nrhoij_in=size(pawrhoij_in)
 nrhoij_out=size(pawrhoij_out);nrhoij_out_all=nrhoij_out
 if (paral_atom) then
   LIBPAW_ALLOCATE(nrhoij_out_i,(nproc_atom))
   buf_size(1)=nrhoij_out
   call xmpi_allgather(buf_size,1,nrhoij_out_i,my_comm_atom,ierr)
   nrhoij_out_all=sum(nrhoij_out_i)
 end if
 if (me==master.and.nrhoij_in/=nrhoij_out_all) then
   msg='pawrhoij_in or pawrhoij_out wrongly allocated!'
   MSG_BUG(msg)
 end if

!Retrieve table(s) of atoms (if necessary)
 if (paral_atom) then
   nullify(my_atmtab)
   call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom, &
&                     nrhoij_out_all,my_natom_ref=nrhoij_out)
   LIBPAW_ALLOCATE(disp_int,(nproc_atom))
   disp_int(1)=0
   do iproc=2,nproc_atom
     disp_int(iproc)=disp_int(iproc-1)+nrhoij_out_i(iproc-1)
   end do
   LIBPAW_ALLOCATE(atmtab,(nrhoij_in))
   call xmpi_gatherv(my_atmtab,nrhoij_out,atmtab,nrhoij_out_i,disp_int,&
&                    master,my_comm_atom,ierr)
   LIBPAW_DEALLOCATE(disp_int)
 end if

!Compute sizes of input buffers and broadcast them
 LIBPAW_ALLOCATE(buf_int_size_i,(nrhoij_out_all))
 LIBPAW_ALLOCATE(buf_dp_size_i ,(nrhoij_out_all))
 if (me==master) then
   buf_int_size_i(:)=0;buf_dp_size_i(:)=0
   do irhoij=1,nrhoij_in
     jrhoij=irhoij;if (paral_atom) jrhoij=atmtab(irhoij)
     cplex       =pawrhoij_in(jrhoij)%cplex
     lmn2_size   =pawrhoij_in(jrhoij)%lmn2_size
     nspden      =pawrhoij_in(jrhoij)%nspden
     lmnmix      =pawrhoij_in(jrhoij)%lmnmix_sz
     ngrhoij     =pawrhoij_in(jrhoij)%ngrhoij
     use_rhoijp  =pawrhoij_in(jrhoij)%use_rhoijp
     use_rhoijres=pawrhoij_in(jrhoij)%use_rhoijres
     use_rhoijim =pawrhoij_in(jrhoij)%use_rhoijim
     use_rhoij_  =pawrhoij_in(jrhoij)%use_rhoij_
     buf_int_size_i(irhoij)=buf_int_size_i(irhoij)+16
     if (ngrhoij>0) buf_dp_size_i(irhoij)=buf_dp_size_i(irhoij)+cplex*lmn2_size*nspden*ngrhoij
     if (use_rhoijres>0) buf_dp_size_i(irhoij)=buf_dp_size_i(irhoij)+cplex*lmn2_size*nspden
     if (use_rhoijim>0)  buf_dp_size_i(irhoij)=buf_dp_size_i(irhoij)+lmn2_size*nspden
     if (use_rhoijp>0) then
       nselect=pawrhoij_in(jrhoij)%nrhoijsel
       buf_int_size_i(irhoij)=buf_int_size_i(irhoij)+nselect
       buf_dp_size_i(irhoij)=buf_dp_size_i(irhoij)+cplex*nselect*nspden
     end if
     if (use_rhoij_>0) then
       rhoij_size2=size(pawrhoij_in(jrhoij)%rhoij_,dim=2)
       buf_dp_size_i(irhoij)=buf_dp_size_i(irhoij)+cplex*lmn2_size*rhoij_size2
     end if
   end do
 end if
 call xmpi_bcast(buf_int_size_i,master,mpicomm,ierr)
 call xmpi_bcast(buf_dp_size_i,master,mpicomm,ierr)
 buf_int_size_all=sum(buf_int_size_i) ; buf_dp_size_all=sum(buf_dp_size_i)

!Prepare buffers/tabs for communication
 if (paral_atom) then
   LIBPAW_ALLOCATE(count_int,(nproc_atom))
   LIBPAW_ALLOCATE(count_dp,(nproc_atom))
   LIBPAW_ALLOCATE(disp_int,(nproc_atom))
   LIBPAW_ALLOCATE(disp_dp,(nproc_atom))
   indx_int=0
   do iproc=1,nproc_atom
     ii=nrhoij_out_i(iproc)
     count_int(iproc)=sum(buf_int_size_i(indx_int+1:indx_int+ii))
     count_dp (iproc)=sum(buf_dp_size_i (indx_int+1:indx_int+ii))
     indx_int=indx_int+ii
   end do
   disp_int(1)=0;disp_dp(1)=0
   do iproc=2,nproc_atom
     disp_int(iproc)=disp_int(iproc-1)+count_int(iproc-1)
     disp_dp (iproc)=disp_dp (iproc-1)+count_dp (iproc-1)
   end do
   if (buf_int_size_all/=sum(count_int).or.buf_dp_size_all/=sum(count_dp)) then
     msg='(1) Wrong buffer sizes !'
     MSG_BUG(msg)
   end if
   buf_int_size=count_int(me_atom+1)
   buf_dp_size =count_dp(me_atom+1)
   LIBPAW_ALLOCATE(buf_int,(buf_int_size))
   LIBPAW_ALLOCATE(buf_dp ,(buf_dp_size))
   if (me==master) then
     LIBPAW_POINTER_ALLOCATE(buf_int_all,(buf_int_size_all))
     LIBPAW_POINTER_ALLOCATE(buf_dp_all ,(buf_dp_size_all))
   else
     LIBPAW_POINTER_ALLOCATE(buf_int_all,(1))
     LIBPAW_POINTER_ALLOCATE(buf_dp_all ,(1))
   end if
   LIBPAW_DEALLOCATE(nrhoij_out_i)
 else
   buf_int_size=buf_int_size_all
   buf_dp_size =buf_dp_size_all
   LIBPAW_ALLOCATE(buf_int,(buf_int_size))
   LIBPAW_ALLOCATE(buf_dp ,(buf_dp_size))
   buf_int_all => buf_int
   buf_dp_all  => buf_dp
 end if
 LIBPAW_DEALLOCATE(buf_int_size_i)
 LIBPAW_DEALLOCATE(buf_dp_size_i)

!Fill input buffers
 if (me==master) then
   indx_int=1;indx_dp =1
   do irhoij=1,nrhoij_in
     jrhoij=irhoij;if (paral_atom) jrhoij=atmtab(irhoij)
     cplex       =pawrhoij_in(jrhoij)%cplex
     lmn2_size   =pawrhoij_in(jrhoij)%lmn2_size
     nspden      =pawrhoij_in(jrhoij)%nspden
     lmnmix      =pawrhoij_in(jrhoij)%lmnmix_sz
     ngrhoij     =pawrhoij_in(jrhoij)%ngrhoij
     use_rhoijp  =pawrhoij_in(jrhoij)%use_rhoijp
     nselect     =pawrhoij_in(jrhoij)%nrhoijsel
     use_rhoijres=pawrhoij_in(jrhoij)%use_rhoijres
     use_rhoijim =pawrhoij_in(jrhoij)%use_rhoijim
     use_rhoij_  =pawrhoij_in(jrhoij)%use_rhoij_
     rhoij_size2 =size(pawrhoij_in(jrhoij)%rhoij_,dim=2)
     buf_int_all(indx_int)=jrhoij                      ;indx_int=indx_int+1 ! Not used !
     buf_int_all(indx_int)=cplex                       ;indx_int=indx_int+1
     buf_int_all(indx_int)=lmn2_size                   ;indx_int=indx_int+1
     buf_int_all(indx_int)=nspden                      ;indx_int=indx_int+1
     buf_int_all(indx_int)=nselect                     ;indx_int=indx_int+1
     buf_int_all(indx_int)=lmnmix                      ;indx_int=indx_int+1
     buf_int_all(indx_int)=ngrhoij                     ;indx_int=indx_int+1
     buf_int_all(indx_int)=use_rhoijp                  ;indx_int=indx_int+1
     buf_int_all(indx_int)=use_rhoijres                ;indx_int=indx_int+1
     buf_int_all(indx_int)=use_rhoijim                 ;indx_int=indx_int+1
     buf_int_all(indx_int)=use_rhoij_                  ;indx_int=indx_int+1
     buf_int_all(indx_int)=rhoij_size2                 ;indx_int=indx_int+1
     buf_int_all(indx_int)=pawrhoij_in(jrhoij)%itypat  ;indx_int=indx_int+1
     buf_int_all(indx_int)=pawrhoij_in(jrhoij)%lmn_size;indx_int=indx_int+1
     buf_int_all(indx_int)=pawrhoij_in(jrhoij)%nsppol  ;indx_int=indx_int+1
     buf_int_all(indx_int)=pawrhoij_in(jrhoij)%nspinor ;indx_int=indx_int+1
     if (use_rhoijp>0) then
       buf_int_all(indx_int:indx_int+nselect-1)=pawrhoij_in(jrhoij)%rhoijselect(1:nselect)
       indx_int=indx_int+nselect
       do isp=1,nspden
         buf_dp_all(indx_dp:indx_dp+cplex*nselect-1)=pawrhoij_in(jrhoij)%rhoijp(1:cplex*nselect,isp)
         indx_dp=indx_dp+cplex*nselect
       end do
     end if
     if (lmnmix>0) then
       buf_int_all(indx_int:indx_int+lmnmix-1)=pawrhoij_in(jrhoij)%kpawmix(1:lmnmix)
       indx_int=indx_int+lmnmix
     end if
     if (ngrhoij>0) then
       do isp=1,nspden
         do ii=1,cplex*lmn2_size
           buf_dp_all(indx_dp:indx_dp+ngrhoij-1)=pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,ii,isp)
           indx_dp=indx_dp+ngrhoij
         end do
       end do
     end if
     if (use_rhoijres>0) then
       do isp=1,nspden
         buf_dp_all(indx_dp:indx_dp+cplex*lmn2_size-1)=pawrhoij_in(jrhoij)%rhoijres(1:cplex*lmn2_size,isp)
         indx_dp=indx_dp+cplex*lmn2_size
       end do
     end if
     if (use_rhoijim>0) then
       do isp=1,nspden
         buf_dp_all(indx_dp:indx_dp+lmn2_size-1)=pawrhoij_in(jrhoij)%rhoijim(1:lmn2_size,isp)
         indx_dp=indx_dp+lmn2_size
       end do
     end if
     if (use_rhoij_>0) then
       do isp=1,rhoij_size2
         buf_dp_all(indx_dp:indx_dp+cplex*lmn2_size-1)=pawrhoij_in(jrhoij)%rhoij_(1:cplex*lmn2_size,isp)
         indx_dp=indx_dp+cplex*lmn2_size
       end do
     end if
   end do
! Check
  if ((indx_int-1/=buf_int_size_all).or.(indx_dp-1/=buf_dp_size_all)) then
    msg='(2) Wrong buffer sizes !'
    MSG_BUG(msg)
  end if
 end if ! me=master

!Communicate
 if (paral_atom) then
   call xmpi_scatterv(buf_int_all,count_int,disp_int,buf_int,buf_int_size,master,mpicomm,ierr)
   call xmpi_scatterv(buf_dp_all ,count_dp ,disp_dp ,buf_dp ,buf_dp_size ,master,mpicomm,ierr)
 else
   call xmpi_bcast(buf_int,master,mpicomm,ierr)
   call xmpi_bcast(buf_dp ,master,mpicomm,ierr)
 end if

!Retrieve data from output buffer
 indx_int=1;indx_dp=1
 call pawrhoij_free(pawrhoij_out)
 do irhoij=1,nrhoij_out
   jrhoij      =buf_int(indx_int);indx_int=indx_int+1 ! Not used !
   cplex       =buf_int(indx_int);indx_int=indx_int+1
   lmn2_size   =buf_int(indx_int);indx_int=indx_int+1
   nspden      =buf_int(indx_int);indx_int=indx_int+1
   nselect     =buf_int(indx_int);indx_int=indx_int+1
   lmnmix      =buf_int(indx_int);indx_int=indx_int+1
   ngrhoij     =buf_int(indx_int);indx_int=indx_int+1
   use_rhoijp  =buf_int(indx_int);indx_int=indx_int+1
   use_rhoijres=buf_int(indx_int);indx_int=indx_int+1
   use_rhoijim =buf_int(indx_int);indx_int=indx_int+1
   use_rhoij_  =buf_int(indx_int);indx_int=indx_int+1
   rhoij_size2 =buf_int(indx_int);indx_int=indx_int+1
   pawrhoij_out(irhoij)%itypat=buf_int(indx_int)  ;indx_int=indx_int+1
   pawrhoij_out(irhoij)%lmn_size=buf_int(indx_int);indx_int=indx_int+1
   pawrhoij_out(irhoij)%nsppol=buf_int(indx_int)  ;indx_int=indx_int+1
   pawrhoij_out(irhoij)%nspinor=buf_int(indx_int) ;indx_int=indx_int+1
   pawrhoij_out(irhoij)%cplex=cplex
   pawrhoij_out(irhoij)%lmn2_size=lmn2_size
   pawrhoij_out(irhoij)%nspden=nspden
   pawrhoij_out(irhoij)%nrhoijsel=nselect
   pawrhoij_out(irhoij)%lmnmix_sz=lmnmix
   pawrhoij_out(irhoij)%ngrhoij=ngrhoij
   pawrhoij_out(irhoij)%use_rhoijp=use_rhoijp
   pawrhoij_out(irhoij)%use_rhoijres=use_rhoijres
   pawrhoij_out(irhoij)%use_rhoijim =use_rhoijim
   pawrhoij_out(irhoij)%use_rhoij_=use_rhoij_
   if (use_rhoijp>0) then
     LIBPAW_ALLOCATE(pawrhoij_out(irhoij)%rhoijselect,(nselect))
     pawrhoij_out(irhoij)%rhoijselect(1:nselect)=buf_int(indx_int:indx_int+nselect-1)
     indx_int=indx_int+nselect
     LIBPAW_ALLOCATE(pawrhoij_out(irhoij)%rhoijp,(cplex*nselect,nspden))
     do isp=1,nspden
       pawrhoij_out(irhoij)%rhoijp(1:cplex*nselect,isp)=buf_dp(indx_dp:indx_dp+cplex*nselect-1)
       indx_dp=indx_dp+cplex*nselect
     end do
   end if
   if (lmnmix>0) then
     LIBPAW_ALLOCATE(pawrhoij_out(irhoij)%kpawmix,(lmnmix))
     pawrhoij_out(irhoij)%kpawmix(1:lmnmix)=buf_int(indx_int:indx_int+lmnmix-1)
     indx_int=indx_int+lmnmix
   end if
   if (ngrhoij>0) then
     LIBPAW_ALLOCATE(pawrhoij_out(irhoij)%grhoij,(ngrhoij,cplex*lmn2_size,nspden))
     do isp=1,nspden
       do ii=1,cplex*lmn2_size
         pawrhoij_out(irhoij)%grhoij(1:ngrhoij,ii,isp)=buf_dp(indx_dp:indx_dp+ngrhoij-1)
         indx_dp=indx_dp+ngrhoij
       end do
     end do
   end if
   if (use_rhoijres>0) then
     LIBPAW_ALLOCATE(pawrhoij_out(irhoij)%rhoijres,(cplex*lmn2_size,nspden))
     do isp=1,nspden
       pawrhoij_out(irhoij)%rhoijres(1:cplex*lmn2_size,isp)=buf_dp(indx_dp:indx_dp+cplex*lmn2_size-1)
       indx_dp=indx_dp+cplex*lmn2_size
     end do
   end if
   if (use_rhoijim>0) then
     LIBPAW_ALLOCATE(pawrhoij_out(irhoij)%rhoijim,(lmn2_size,nspden))
     do isp=1,nspden
       pawrhoij_out(irhoij)%rhoijim(1:cplex*lmn2_size,isp)=buf_dp(indx_dp:indx_dp+lmn2_size-1)
       indx_dp=indx_dp+lmn2_size
     end do
   end if
   if (use_rhoij_>0) then
     LIBPAW_ALLOCATE(pawrhoij_out(irhoij)%rhoij_,(cplex*lmn2_size,rhoij_size2))
     do isp=1,rhoij_size2
       pawrhoij_out(irhoij)%rhoij_(1:cplex*lmn2_size,isp)=buf_dp(indx_dp:indx_dp+cplex*lmn2_size-1)
       indx_dp=indx_dp+cplex*lmn2_size
     end do
   end if
 end do
!Check
 if ((indx_int/=1+buf_int_size).or.(indx_dp/=1+buf_dp_size)) then
   msg='(3) Wrong buffer sizes !'
   MSG_BUG(msg)
 end if

!Free memory
 LIBPAW_DEALLOCATE(buf_int)
 LIBPAW_DEALLOCATE(buf_dp)
 if (paral_atom) then
   LIBPAW_POINTER_DEALLOCATE(buf_int_all)
   LIBPAW_POINTER_DEALLOCATE(buf_dp_all)
   LIBPAW_DEALLOCATE(count_int)
   LIBPAW_DEALLOCATE(count_dp)
   LIBPAW_DEALLOCATE(disp_int)
   LIBPAW_DEALLOCATE(disp_dp)
   LIBPAW_DEALLOCATE(atmtab)
   call free_my_atmtab(my_atmtab,my_atmtab_allocated)
 end if

end subroutine pawrhoij_bcast
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/pawrhoij_redistribute
!! NAME
!! pawrhoij_redistribute
!!
!! FUNCTION
!!   Redistribute an array of pawrhoij datastructures
!!   Input pawrhoij is given on a MPI communicator
!!   Output pawrhoij is redistributed on another MPI communicator
!!
!! INPUTS
!!  mpi_comm_in= input MPI (atom) communicator
!!  mpi_comm_out= output MPI (atom) communicator
!!  mpi_atmtab_in= --optional-- indexes of the input pawrhoij treated by current proc
!!                 if not present, will be calculated in the present routine
!!  mpi_atmtab_out= --optional-- indexes of the output pawrhoij treated by current proc
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
!!  [pawrhoij_out(:)]<type(pawrhoij_type)>= --optional--
!!                    if present, the redistributed datastructure does not replace
!!                    the input one but is delivered in pawrhoij_out
!!                    if not present, input and output datastructure are the same.
!!
!! SIDE EFFECTS
!!  pawrhoij(:)<type(pawrhoij_type)>= input (and eventually output) pawrhoij datastructures
!!
!! PARENTS
!!      m_paral_pert
!!
!! CHILDREN
!!
!! SOURCE

 subroutine pawrhoij_redistribute(pawrhoij,mpi_comm_in,mpi_comm_out,&
&           natom,mpi_atmtab_in,mpi_atmtab_out,pawrhoij_out,&
&           SendAtomProc,SendAtomList,RecvAtomProc,RecvAtomList)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawrhoij_redistribute'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpi_comm_in,mpi_comm_out
 integer,optional,intent(in) :: natom
!arrays
 integer,intent(in),optional,target :: mpi_atmtab_in(:),mpi_atmtab_out(:)
 integer,intent(in),optional :: SendAtomProc(:),SendAtomList(:),RecvAtomProc(:),RecvAtomList(:)
 type(pawrhoij_type),allocatable,intent(inout) :: pawrhoij(:)
 type(pawrhoij_type),pointer,intent(out),optional :: pawrhoij_out(:)

!Local variables-------------------------------
!scalars
 integer :: algo_option,i1,iat_in,iat_out,iatom,ierr,ireq,iircv,iisend,imsg,imsg_current,imsg1
 integer :: iproc_rcv,iproc_send,me_exch,mpi_comm_exch,my_natom_in,my_natom_out,my_tag,natom_tot,nb_dp,nb_int
 integer :: nb_msg,nbmsg_incoming,nbrecv,nbrecvmsg,nbsend,nbsendreq,nbsent,next,npawrhoij_sent
 integer :: nproc_in,nproc_out
 logical :: flag,in_place,message_yet_prepared,my_atmtab_in_allocated,my_atmtab_out_allocated,paral_atom
!arrays
 integer :: buf_size(3),request1(3)
 integer,pointer :: my_atmtab_in(:),my_atmtab_out(:)
 integer,allocatable :: atmtab_send(:),atm_indx_in(:),atm_indx_out(:),buf_int1(:),From(:),request(:)
 integer,allocatable,target :: buf_int(:)
 integer,pointer:: buf_ints(:)
 logical,allocatable :: msg_pick(:)
 real(dp),allocatable :: buf_dp1(:)
 real(dp),allocatable,target :: buf_dp(:)
 real(dp),pointer :: buf_dps(:)
 type(coeffi1_type),target,allocatable :: tab_buf_int(:),tab_buf_atom(:)
 type(coeff1_type),target,allocatable :: tab_buf_dp(:)
 type(pawrhoij_type),pointer :: pawrhoij_out1(:)
 type(pawrhoij_type),allocatable :: pawrhoij_all(:)

! *************************************************************************

!@pawrhoij_type

 in_place=(.not.present(pawrhoij_out))
 my_natom_in=size(pawrhoij)

!If not "in_place", destroy the output datastructure
 if (.not.in_place) then
   if (associated(pawrhoij_out)) then
     call pawrhoij_free(pawrhoij_out)
     LIBPAW_DATATYPE_DEALLOCATE(pawrhoij_out)
   end if
 end if

!Special sequential case
 if (mpi_comm_in==xmpi_comm_self.and.mpi_comm_out==xmpi_comm_self) then
   if ((.not.in_place).and.(my_natom_in>0)) then
     LIBPAW_DATATYPE_ALLOCATE(pawrhoij_out,(my_natom_in))
     call pawrhoij_nullify(pawrhoij_out)
     call pawrhoij_copy(pawrhoij,pawrhoij_out,&
&                    keep_cplex=.false.,keep_itypat=.false.,keep_nspden=.false.)
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

   LIBPAW_DATATYPE_ALLOCATE(pawrhoij_all,(natom_tot))
   call pawrhoij_nullify(pawrhoij_all)
   call pawrhoij_copy(pawrhoij,pawrhoij_all,comm_atom=mpi_comm_in,mpi_atmtab=my_atmtab_in,&
&                  keep_cplex=.false.,keep_itypat=.false.,keep_nspden=.false.)
   if (in_place) then
     call pawrhoij_free(pawrhoij)
     LIBPAW_DATATYPE_DEALLOCATE(pawrhoij)
     LIBPAW_DATATYPE_ALLOCATE(pawrhoij,(my_natom_out))
     call pawrhoij_nullify(pawrhoij)
     call pawrhoij_copy(pawrhoij_all,pawrhoij,comm_atom=mpi_comm_out,mpi_atmtab=my_atmtab_out,&
&                       keep_cplex=.false.,keep_itypat=.false.,keep_nspden=.false.)
   else
     LIBPAW_DATATYPE_ALLOCATE(pawrhoij_out,(my_natom_out))
     call pawrhoij_nullify(pawrhoij_out)
     call pawrhoij_copy(pawrhoij_all,pawrhoij_out,comm_atom=mpi_comm_out,mpi_atmtab=my_atmtab_out,&
&                       keep_cplex=.false.,keep_itypat=.false.,keep_nspden=.false.)
   end if
   call pawrhoij_free(pawrhoij_all)
   LIBPAW_DATATYPE_DEALLOCATE(pawrhoij_all)


!Asynchronous algorithm (asynchronous communications)
!---------------------------------------------------------
 else if (algo_option==2) then

   nbsend=size(SendAtomProc) ; nbrecv=size(RecvAtomProc)

   if (in_place) then
     if ( my_natom_out > 0 ) then
       LIBPAW_DATATYPE_ALLOCATE(pawrhoij_out1,(my_natom_out))
       call pawrhoij_nullify(pawrhoij_out1)
     else
       LIBPAW_DATATYPE_ALLOCATE(pawrhoij_out1,(0))
     end if
   else
     LIBPAW_DATATYPE_ALLOCATE(pawrhoij_out,(my_natom_out))
     call pawrhoij_nullify(pawrhoij_out)
     pawrhoij_out1=>pawrhoij_out
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
         if (nbsent > 0) then
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
           end do
!          Create the message
           if (.not.message_yet_prepared) then
             nb_msg=nb_msg+1
             call pawrhoij_isendreceive_fillbuffer( &
&                   pawrhoij,atmtab_send,atm_indx_in,nbsent,buf_int,nb_int,buf_dp,nb_dp)
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
           my_tag=100
           ireq=ireq+1
           call xmpi_isend(buf_size,iproc_rcv,my_tag,mpi_comm_exch,request(ireq),ierr)
           my_tag=101
           ireq=ireq+1
           call xmpi_isend(buf_ints,iproc_rcv,my_tag,mpi_comm_exch,request(ireq),ierr)
           my_tag=102
           ireq=ireq+1
           call xmpi_isend(buf_dps,iproc_rcv,my_tag,mpi_comm_exch,request(ireq),ierr)
           nbsendreq=ireq
           nbsent=0
         end if
       end if
     else ! Just a renumbering, not a sending
       iat_in=atm_indx_in(SendAtomList(iisend))
       iat_out=atm_indx_out(my_atmtab_in(iat_in))
       call pawrhoij_copy(pawrhoij(iat_in:iat_in),pawrhoij_out1(iat_out:iat_out), &
&                         keep_cplex=.false.,keep_itypat=.false.,keep_nspden=.false.)
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
         my_tag=100
         call xmpi_iprobe(iproc_send,my_tag,mpi_comm_exch,flag,ierr)
         if (flag) then
           msg_pick(i1)=.true.
           call xmpi_irecv(buf_size,iproc_send,my_tag,mpi_comm_exch,request1(1),ierr)
           call xmpi_wait(request1(1),ierr)
           nb_int=buf_size(1)
           nb_dp=buf_size(2)
           npawrhoij_sent=buf_size(3)
           LIBPAW_ALLOCATE(buf_int1,(nb_int))
           LIBPAW_ALLOCATE(buf_dp1,(nb_dp))
           my_tag=101
           call xmpi_irecv(buf_int1,iproc_send,my_tag,mpi_comm_exch,request1(2),ierr)
           my_tag=102
           call xmpi_irecv(buf_dp1,iproc_send,my_tag,mpi_comm_exch,request1(3),ierr)
           call xmpi_waitall(request1(2:3),ierr)
           call pawrhoij_isendreceive_getbuffer(pawrhoij_out1,npawrhoij_sent,atm_indx_out,buf_int1,buf_dp1)
           nbmsg_incoming=nbmsg_incoming-1
           LIBPAW_DEALLOCATE(buf_int1)
           LIBPAW_DEALLOCATE(buf_dp1)
         end if
       end if
     end do
   end do
   LIBPAW_DEALLOCATE(msg_pick)

   if (in_place) then
     call pawrhoij_free(pawrhoij)
     LIBPAW_DATATYPE_DEALLOCATE(pawrhoij)
     LIBPAW_DATATYPE_ALLOCATE(pawrhoij,(my_natom_out))
     call pawrhoij_nullify(pawrhoij)
     call pawrhoij_copy(pawrhoij_out1,pawrhoij, &
&         keep_cplex=.false.,keep_itypat=.false.,keep_nspden=.false.)
     call pawrhoij_free(pawrhoij_out1)
     LIBPAW_DATATYPE_DEALLOCATE(pawrhoij_out1)
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

end subroutine pawrhoij_redistribute
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/pawrhoij_io
!! NAME
!! pawrhoij_io
!!
!! FUNCTION
!! IO method for pawrhoij datastructures.
!!
!! INPUTS
!!  unitfi=Unit number for IO file or netcdf file handler (already opened in the caller).
!!  nsppol_in=Number of independent spin polarizations. Only used for reading.
!!  nspinor_in=Number of spinorial components. Only used for reading.
!!  nspden_in=Number of spin-density components. only used for reading.
!!  nlmn_type(ntypat)= Number of (l,m,n) elements for the paw basis for each type of atom. Only used for reading.
!!  typat(natom) =Type of each atom.
!!  headform=Format of the abinit header (only used for reading as we need to know how to read
!!    the data. Writing is always done using the latest headform.
!!  rdwr_mode(len=*)=String defining the IO mode. Possible values (not case sensitive):
!!    "W"= For writing to unitfi
!!    "R"= For reading from unitfi
!!    "E"= For echoing.
!!    "D"= for debug
!!  [form(len=*)]= String defining the file format. Defaults to Fortran binary mode i.e., "unformatted"
!!  Other possible values are (case insensitive):
!!    "formatted"=For IO on a file open in formatted mode.
!!    "netcdf"=For IO on a netcdf file.
!!  [natinc]=Defines the increment in the loop over natom used for echoing the pawrhoij(natom) datastructures.
!!    If not specified, only the first and the last atom will be printed.
!!
!! SIDE EFFECTS
!!  pawrhoij(:)<type(pawrhoij_type)>= rhoij datastructure.
!!   if rdwr_mode="W", it will be written on unit unitfi using the file format given by form.
!!   if rdwr_mode="R", pawrhoij will be read and initialized from unit unitfi that has been
!!      opened with form=form.
!!   if rdwr_mode="E", the routines only echoes the content of the structure.
!!
!! PARENTS
!!      m_hdr,m_qparticles
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawrhoij_io(pawrhoij,unitfi,nsppol_in,nspinor_in,nspden_in,nlmn_type,typat,&
&                   headform,rdwr_mode,form,natinc,mpi_atmtab)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawrhoij_io'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: unitfi,headform,nspden_in,nspinor_in,nsppol_in
 integer,optional,intent(in) :: natinc
 character(len=*),intent(in) :: rdwr_mode
 character(len=*),optional,intent(in) :: form
 integer, intent(in), optional, pointer :: mpi_atmtab(:)
!arrays
 integer,intent(in) :: typat(:),nlmn_type(:)
 type(pawrhoij_type),intent(inout) :: pawrhoij(:)

!Local variables-------------------------------
!scalars
 integer,parameter :: fort_formatted=1,fort_binary=2,netcdf_io=3
 integer :: cplex,i1,i2,iatom,iatom_tot,natom,ispden,bsize,ii,jj,lmn2_size
 integer :: nselect,my_cplex,my_natinc,my_natom,my_nspden,ngrhoijmx,size_rhoij2
 integer :: iomode,ncid,natom_id,cplex_id,nspden_id,nsel56_id,buffer_id,ibuffer_id,ncerr
 integer :: bsize_id,bufsize_id
 logical :: paral_atom
 character(len=500) :: msg
!arrays
 integer,allocatable :: ibuffer(:),nsel44(:,:),nsel56(:)
 integer,pointer :: my_atmtab(:)
 real(dp), allocatable :: buffer(:)

! *************************************************************************

 my_natom=SIZE(pawrhoij);if (my_natom==0) return
 my_nspden=nspden_in
 natom=size(typat)
 paral_atom=(my_natom/=natom)
 if (present(mpi_atmtab)) then
   if (.not.associated(mpi_atmtab)) then
     msg='mpi_atmtab not associated (pawrhoij_io)'
     MSG_BUG(msg)
   end if
   my_atmtab=>mpi_atmtab
 else if (my_natom/=natom) then
   msg='my_natom /=natom, mpi_atmtab should be in argument (pawrhoij_io)'
   MSG_BUG(msg)
 end if

 iomode = fort_binary
 if (PRESENT(form)) then
   select case (libpaw_to_upper(form))
   case ("FORMATTED")
     iomode = fort_formatted
   case ("NETCDF")
     iomode = netcdf_io
   case default
     MSG_ERROR("Wrong form: "//trim(form))
   end select
 end if

#ifndef LIBPAW_HAVE_NETCDF
 if (iomode == netcdf_io) then
   MSG_ERROR("iomode == netcdf_io but netcdf library is missing.")
 end if
#endif
 ncid = unitfi

 select case (rdwr_mode(1:1))

   case ("R","r") ! Reading the Rhoij tab

     if ((headform>=44).and.(headform<56)) then
       LIBPAW_ALLOCATE(nsel44,(nspden_in,natom))
       if (iomode == fort_binary) then
         read(unitfi  ) ((nsel44(ispden,iatom),ispden=1,nspden_in),iatom=1,natom)
       else if (iomode == fort_formatted) then
         read(unitfi,*) ((nsel44(ispden,iatom),ispden=1,nspden_in),iatom=1,natom)
#ifdef LIBPAW_HAVE_NETCDF
       else if (iomode == netcdf_io) then
         MSG_ERROR("header in 44-56 not compatible with Netcdf")
#endif
       end if
       call pawrhoij_alloc(pawrhoij,1,nspden_in,nspinor_in,nsppol_in,typat,lmnsize=nlmn_type)
       do iatom=1,natom
         pawrhoij(iatom)%nrhoijsel=nsel44(1,iatom)
       end do
       bsize=sum(nsel44)
       LIBPAW_ALLOCATE(ibuffer,(bsize))
       LIBPAW_ALLOCATE(buffer,(bsize))
       if (iomode == fort_binary) then
         read(unitfi  ) ibuffer(:),buffer(:)
       else if (iomode == fort_formatted) then
         read(unitfi,*) ibuffer(:),buffer(:)
       end if
       ii=0
       do iatom=1,natom
         nselect=nsel44(1,iatom)
         pawrhoij(iatom)%rhoijselect(1:nselect)=ibuffer(ii+1:ii+nselect)
         do ispden=1,nspden_in
           pawrhoij(iatom)%rhoijp(1:nselect,ispden)=buffer(ii+1:ii+nselect)
           ii=ii+nselect
         end do
       end do
       LIBPAW_DEALLOCATE(ibuffer)
       LIBPAW_DEALLOCATE(buffer)
       LIBPAW_DEALLOCATE(nsel44)
     else if (headform>=56) then
       LIBPAW_ALLOCATE(nsel56,(natom))
       if (headform==56) then
         if (iomode == fort_binary) then
           read(unitfi  ) (nsel56(iatom),iatom=1,natom),my_cplex
         else if (iomode == fort_formatted) then
           read(unitfi,*) (nsel56(iatom),iatom=1,natom),my_cplex
#ifdef LIBPAW_HAVE_NETCDF
       else if (iomode == netcdf_io) then
          NCF_CHECK(nf90_inq_dimid(ncid, "pawrhoij_cplex", cplex_id))
          NCF_CHECK(nf90_inquire_dimension(ncid, cplex_id, len=my_cplex))

          NCF_CHECK(nf90_inq_varid(ncid, "rhoijsel_atoms", nsel56_id))
          NCF_CHECK(nf90_get_var(ncid, nsel56_id, nsel56))
#endif
         end if
       else
         if (iomode == fort_binary) then
           read(unitfi  ) (nsel56(iatom),iatom=1,natom),my_cplex,my_nspden
         else if (iomode == fort_formatted) then
           read(unitfi,*) (nsel56(iatom),iatom=1,natom),my_cplex,my_nspden
#ifdef LIBPAW_HAVE_NETCDF
       else if (iomode == netcdf_io) then
           NCF_CHECK(nf90_inq_dimid(ncid, "pawrhoij_cplex", cplex_id))
           NCF_CHECK(nf90_inquire_dimension(ncid, cplex_id, len=my_cplex))
           NCF_CHECK(nf90_inq_dimid(ncid, "pawrhoij_nspden", nspden_id))
           NCF_CHECK(nf90_inquire_dimension(ncid, nspden_id, len=my_nspden))
           NCF_CHECK(nf90_inq_varid(ncid, "nrhoijsel_atoms", nsel56_id))
           NCF_CHECK(nf90_get_var(ncid, nsel56_id, nsel56))
#endif
         end if
       end if
       call pawrhoij_alloc(pawrhoij,my_cplex,my_nspden,nspinor_in,nsppol_in,typat,lmnsize=nlmn_type)
       do iatom=1,natom
         pawrhoij(iatom)%nrhoijsel=nsel56(iatom)
       end do
       bsize=sum(nsel56)
       LIBPAW_ALLOCATE(ibuffer,(bsize))
       LIBPAW_ALLOCATE(buffer,(bsize*my_nspden*my_cplex))
       if (iomode == fort_binary) then
         read(unitfi  ) ibuffer(:),buffer(:)
       else if (iomode == fort_formatted) then
         read(unitfi,*) ibuffer(:),buffer(:)
#ifdef LIBPAW_HAVE_NETCDF
       else if (iomode == netcdf_io) then
         if (bsize > 0) then
           NCF_CHECK(nf90_inq_varid(ncid, "rhoijselect_atoms", ibuffer_id))
           NCF_CHECK(nf90_get_var(ncid, ibuffer_id, ibuffer))
           NCF_CHECK(nf90_inq_varid(ncid, "rhoijp_atoms", buffer_id))
           NCF_CHECK(nf90_get_var(ncid, buffer_id, buffer))
         end if
#endif
       end if
       ii=0;jj=0
       do iatom=1,natom
         nselect=nsel56(iatom)
         pawrhoij(iatom)%rhoijselect(1:nselect)=ibuffer(ii+1:ii+nselect)
         ii=ii+nselect
         do ispden=1,my_nspden
           pawrhoij(iatom)%rhoijp(1:my_cplex*nselect,ispden)=buffer(jj+1:jj+my_cplex*nselect)
           jj=jj+my_cplex*nselect
         end do
       end do
       LIBPAW_DEALLOCATE(ibuffer)
       LIBPAW_DEALLOCATE(buffer)
       LIBPAW_DEALLOCATE(nsel56)
     end if

   case ("W","w") ! Writing the Rhoij tab (latest format is used)

     LIBPAW_ALLOCATE(nsel56,(natom))
     my_cplex =pawrhoij(1)%cplex
     my_nspden=pawrhoij(1)%nspden
     do iatom=1,natom
       nsel56(iatom)=pawrhoij(iatom)%nrhoijsel
     end do
     bsize=sum(nsel56)

     if (iomode == fort_binary) then
       write(unitfi  ) (nsel56(iatom),iatom=1,natom),my_cplex,my_nspden
     else if (iomode == fort_formatted) then
       write(unitfi,*) (nsel56(iatom),iatom=1,natom),my_cplex,my_nspden
#ifdef LIBPAW_HAVE_NETCDF
     else if (iomode == netcdf_io) then
       ncerr = nf90_redef(ncid)

       ! Define dimensions.
       ncerr = nf90_inq_dimid(ncid, "number_of_atoms", natom_id)
       if (ncerr /= nf90_noerr) then
         NCF_CHECK(nf90_def_dim(ncid, "number_of_atoms", natom, natom_id))
       end if
       NCF_CHECK(nf90_def_var(ncid, "nrhoijsel_atoms", NF90_INT, natom_id, nsel56_id))

       NCF_CHECK(nf90_def_dim(ncid, "pawrhoij_cplex", my_cplex, cplex_id))
       NCF_CHECK(nf90_def_dim(ncid, "pawrhoij_nspden", my_nspden, nspden_id))
       !write(std_out,*)"bsize,my_nspden,my_cplex",bsize, my_nspden, my_cplex
       if (bsize > 0) then
         NCF_CHECK(nf90_def_dim(ncid, "rhoijselect_atoms_dim", bsize, bsize_id))
         NCF_CHECK(nf90_def_dim(ncid, "rhoijp_atoms_dim", bsize*my_nspden*my_cplex, bufsize_id))
         ! Define variables.
         NCF_CHECK(nf90_def_var(ncid, "rhoijselect_atoms", NF90_INT, bsize_id, ibuffer_id))
         NCF_CHECK(nf90_def_var(ncid, "rhoijp_atoms", NF90_DOUBLE, bufsize_id, buffer_id))
       else
         ! This happens in v5[40] and bsize == 0 corresponds to NC_UNLIMITED
         MSG_COMMENT("All rhoij entries are zero. No netcdf entry produced")
       end if

       ! Write nsel56
       NCF_CHECK(nf90_enddef(ncid))
       NCF_CHECK(nf90_put_var(ncid, nsel56_id, nsel56))
#endif
     end if

     LIBPAW_ALLOCATE(ibuffer,(bsize))
     LIBPAW_ALLOCATE(buffer,(bsize*my_nspden*my_cplex))
     ii=0;jj=0
     do iatom=1,natom
       nselect=nsel56(iatom)
       ibuffer(ii+1:ii+nselect)=pawrhoij(iatom)%rhoijselect(1:nselect)
       ii=ii+nselect
       do ispden=1,my_nspden
         buffer(jj+1:jj+my_cplex*nselect)=pawrhoij(iatom)%rhoijp(1:my_cplex*nselect,ispden)
         jj=jj+my_cplex*nselect
       end do
     end do
     if (iomode == fort_binary) then
       write(unitfi  ) ibuffer(:),buffer(:)
     else if (iomode == fort_formatted) then
       write(unitfi,*) ibuffer(:),buffer(:)
#ifdef LIBPAW_HAVE_NETCDF
     else if (iomode == netcdf_io) then
       if (bsize > 0) then
         NCF_CHECK(nf90_put_var(ncid, ibuffer_id, ibuffer))
         NCF_CHECK(nf90_put_var(ncid, buffer_id,  buffer))
       end if
#endif
     end if
     LIBPAW_DEALLOCATE(ibuffer)
     LIBPAW_DEALLOCATE(buffer)
     LIBPAW_DEALLOCATE(nsel56)

   case ("E","e") ! Echoing the Rhoij tab

     my_natinc=1; if(natom>1) my_natinc=natom-1
     if (PRESENT(natinc)) my_natinc = natinc ! user-defined increment.
     LIBPAW_ALLOCATE(ibuffer,(0))
     do iatom=1,my_natom,my_natinc
       iatom_tot=iatom;if(paral_atom)iatom_tot=my_atmtab(iatom)
       do ispden=1,pawrhoij(iatom)%nspden
         write(unitfi, '(a,i4,a,i1,a)' ) ' rhoij(',iatom_tot,',',ispden,')=  (max 12 non-zero components will be written)'
         call pawio_print_ij(unitfi,pawrhoij(iatom)%rhoijp(:,ispden),&
&         pawrhoij(iatom)%nrhoijsel,&
&         pawrhoij(iatom)%cplex,&
&         pawrhoij(iatom)%lmn_size,-1,ibuffer,1,0,&
&         pawrhoij(iatom)%rhoijselect(:),-1.d0,1,mode_paral='PERS')
       end do
     end do
     LIBPAW_DEALLOCATE(ibuffer)

  case ("D","d") ! Debug

     write(unitfi,'(a,i4)' ) 'size pawmkrhoij , natom = ' , my_natom
     my_natinc=1;  if(natom>1) my_natinc=natom-1
     if (PRESENT(natinc)) my_natinc = natinc ! user-defined increment.
     LIBPAW_ALLOCATE(ibuffer,(0))
     do iatom=1,my_natom,my_natinc
       iatom_tot=iatom;if(paral_atom) iatom_tot=my_atmtab(iatom)
       if (iatom_tot/=1) cycle
       write(unitfi,'(a,i4,a)' ) ' *******  rhoij (Atom # ',iatom_tot,' ********)'
       write(unitfi,'(a,i4,a,i4)' ) 'cplx=',pawrhoij(iatom)%cplex, ' nselect=', pawrhoij(iatom)%nrhoijsel
       write(unitfi,'(a,i4,a,i4)' ) 'nspden=',pawrhoij(iatom)%nspden, ' lmn2size=',pawrhoij(iatom)%lmn2_size
       write(unitfi,'(a,i4,a,i4)' ) 'lmnmix=',pawrhoij(iatom)%lmnmix_sz, ' ngrhoij=',pawrhoij(iatom)%ngrhoij
       write(unitfi,'(a,i4,a,i4)' ) 'use_rhoijres=',pawrhoij(iatom)%use_rhoijres, &
&                                   'use_rhoij_=',pawrhoij(iatom)%use_rhoij_
       write(unitfi,'(a,i4)' )      'use_rhoijim=',pawrhoij(iatom)%use_rhoijim
       write(unitfi,'(a,i4,a,i4)' ) 'itypat=',pawrhoij(iatom)%itypat, ' lmn_size=',pawrhoij(iatom)%lmn_size
       write(unitfi,'(a,i4,a,i4)' ) 'nsppol=',pawrhoij(iatom)%nsppol, ' nspinor=',pawrhoij(iatom)%nspinor
       cplex=pawrhoij(iatom)%cplex
       lmn2_size=pawrhoij(iatom)%lmn2_size
       do i2=1,pawrhoij(iatom)%nrhoijsel
         write(unitfi,'(a,i4,a,i4,a,i4,a,f9.5)') 'rhoijselect(,',i2,')=',&
&          pawrhoij(iatom)%rhoijselect(i2)
       end do
       if (pawrhoij(iatom)%ngrhoij>0) then
         ngrhoijmx=2; if (pawrhoij(iatom)%ngrhoij<ngrhoijmx) ngrhoijmx=pawrhoij(iatom)%ngrhoij
         do ispden=1,pawrhoij(iatom)%nspden
           do i1=ngrhoijmx,ngrhoijmx
             do i2=cplex*lmn2_size,cplex*lmn2_size
               write(unitfi,'(a,i4,a,i4,a,i4,a,f9.5)') ' grhoij(,',i1,',',i2,',',ispden,')=',&
&                pawrhoij(iatom)%grhoij(i1,i2,ispden)
             end do
           end do
         end do
         call libpaw_flush(unitfi)
       end if
       if (pawrhoij(iatom)%use_rhoijres>0) then
         do ispden=1,pawrhoij(iatom)%nspden
           do i2=cplex*lmn2_size,cplex*lmn2_size
             write(unitfi,'(a,i4,a,i4,a,f9.5)') ' rhoijres(,',i2,',ispden=',ispden,')=',&
&              pawrhoij(iatom)%rhoijres(i2,ispden)
           end do
         end do
         call libpaw_flush(unitfi)
       end if
       if (pawrhoij(iatom)%use_rhoijim>0) then
         do ispden=1,pawrhoij(iatom)%nspden
           do i2=lmn2_size,lmn2_size
             write(unitfi,'(a,i4,a,i4,a,f9.5)') ' rhoijim(,',i2,',ispden=',ispden,')=',&
&              pawrhoij(iatom)%rhoijim(i2,ispden)
           end do
         end do
         call libpaw_flush(unitfi)
       end if
       if (pawrhoij(iatom)%nrhoijsel>0) then
         do ispden=1,pawrhoij(iatom)%nspden
           do i2=cplex*pawrhoij(iatom)%nrhoijsel ,cplex*pawrhoij(iatom)%nrhoijsel
             write(unitfi,'(a,i4,a,i4,a,f9.5)') ' rhoijp!(nrhoijselec=,',i2,',ispden=',ispden,')=',&
&              pawrhoij(iatom)%rhoijp(i2,ispden)
           end do
         end do
         call libpaw_flush(unitfi)
       end if
       if (pawrhoij(iatom)%use_rhoij_>0) then
         size_rhoij2=size(pawrhoij(iatom)%rhoij_,2)
         do ispden=1,size_rhoij2
           do i2=cplex*lmn2_size,cplex*lmn2_size
             write(unitfi,'(a,i4,a,i4,a,f9.5)') ' rhoij_(,',i2,',ispden=',ispden,')=',&
&              pawrhoij(iatom)%rhoij_(i2,ispden)
           end do
         end do
       end if
       call libpaw_flush(unitfi)
       if (pawrhoij(iatom)%lmnmix_sz>0) then
         write(unitfi,'(a)') 'kpawmix '
         write(unitfi,'(i4,i4,i4,i4)') pawrhoij(iatom)%kpawmix(:)
       end if
       call libpaw_flush(unitfi)
     end do

   case default
     msg='Wrong rdwr_mode'//TRIM(rdwr_mode)
     MSG_ERROR(msg)

 end select

end subroutine pawrhoij_io
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/pawrhoij_unpack
!! NAME
!! pawrhoij_unpack
!!
!! FUNCTION
!!  Unpack the values store in rhoijp copying them to the rhoij_ array.
!!
!! SIDE EFFECTS
!!  rhoij(:) <pawrhoij_type)>= input/output datastructure
!!   * In output the rhoij_ array is filled with the values stored in the packed array rhoijp.
!!   * If use_rhoij_/=1, rhoij_ is allocated and the corresponding flag is set to 1.
!!
!! PARENTS
!!      paw_qpscgw
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawrhoij_unpack(rhoij)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawrhoij_unpack'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 type(pawrhoij_type),intent(inout) :: rhoij(:)

!Local variables-------------------------------
 integer :: natom,iat,lmn2_size,isel,klmn,nspden,cplex

! *************************************************************************

 natom  = SIZE(rhoij) ; if (natom==0) return
 nspden = rhoij(1)%nspden    ! MT jan-2010: this should not be nspden but nsppol or 4 if nspden=4
 cplex  = rhoij(1)%cplex

 do iat=1,natom

   lmn2_size =rhoij(iat)%lmn2_size

   if (rhoij(iat)%use_rhoij_/=1) then ! Have to allocate rhoij
     LIBPAW_ALLOCATE(rhoij(iat)%rhoij_,(cplex*lmn2_size,nspden))
     rhoij(iat)%use_rhoij_=1
   end if
   rhoij(iat)%rhoij_ = zero

   do isel=1,rhoij(iat)%nrhoijsel ! Looping over non-zero ij elements.
     klmn = rhoij(iat)%rhoijselect(isel)
     rhoij(iat)%rhoij_(klmn,:) = rhoij(iat)%rhoijp(isel,:)
   end do

 end do ! natom

end subroutine pawrhoij_unpack
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/pawrhoij_init_unpacked
!! NAME
!! pawrhoij_init_unpacked
!!
!! FUNCTION
!!  Initialize field of rhoij datastructure for unpacked values (pawrhoij%rhoij_ array)
!!
!! SIDE EFFECTS
!!  rhoij(:) <pawrhoij_type)>= input/output datastructure
!!   * In output the rhoij_ array is allocated
!!
!! PARENTS
!!      dfpt_nstpaw,dfpt_rhofermi,dfpt_vtorho,energy,pawmkrhoij
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawrhoij_init_unpacked(rhoij)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawrhoij_init_unpacked'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 type(pawrhoij_type),intent(inout) :: rhoij(:)

!Local variables-------------------------------
 integer :: iat,nrhoij,nsp2

! *************************************************************************

 nrhoij=SIZE(rhoij);if (nrhoij==0) return
 nsp2=rhoij(1)%nsppol;if (rhoij(1)%nspden==4) nsp2=4

 do iat=1,nrhoij

   if (allocated(rhoij(iat)%rhoij_))  then
     LIBPAW_DEALLOCATE(rhoij(iat)%rhoij_)
   end if
   LIBPAW_ALLOCATE(rhoij(iat)%rhoij_,(rhoij(iat)%cplex*rhoij(iat)%lmn2_size,nsp2))
   rhoij(iat)%use_rhoij_=1
   rhoij(iat)%rhoij_=zero

 end do

end subroutine pawrhoij_init_unpacked
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/pawrhoij_free_unpacked
!! NAME
!! pawrhoij_free_unpacked
!!
!! FUNCTION
!!  Destroy field of rhoij datastructure for unpacked values (pawrhoij%rhoij_ array)
!!
!! SIDE EFFECTS
!!  rhoij(:) <pawrhoij_type)>= input/output datastructure
!!   * In output the rhoij_ array is deallocated
!!
!! PARENTS
!!      dfpt_rhofermi,energy,pawmkrho
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawrhoij_free_unpacked(rhoij)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawrhoij_free_unpacked'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 type(pawrhoij_type),intent(inout) :: rhoij(:)

!Local variables-------------------------------
 integer :: iat,nrhoij

! *************************************************************************

 nrhoij=SIZE(rhoij);if (nrhoij==0) return

 do iat=1,nrhoij

   if (allocated(rhoij(iat)%rhoij_))  then
     LIBPAW_DEALLOCATE(rhoij(iat)%rhoij_)
   end if
   rhoij(iat)%use_rhoij_=0

 end do

end subroutine pawrhoij_free_unpacked
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/pawrhoij_mpisum_unpacked_1D
!! NAME
!! pawrhoij_mpisum_unpacked_1D
!!
!! FUNCTION
!! Build the MPI sum of the unsymmetrized PAW rhoij_ (augmentation occupancies)
!! Remember:for each atom, rho_ij=Sum_{n,k} {occ(n,k)*<Cnk|p_i><p_j|Cnk>}
!! Target: 1D array of pawrhoij datastructures
!!
!! INPUTS
!!  comm1=MPI communicator. Data will be MPI summed inside comm1
!!  [comm2]=second MPI communicator. If present, rhoij_ will be
!!          MPI summed inside comm2 after the collective sum in comm1.
!!
!! SIDE EFFECTS
!!  pawrhoij(:) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  Input: the data calculateed by this processor.
!!  Output: the final MPI sum over comm1 and comm2.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawrhoij_mpisum_unpacked_1D(pawrhoij,comm1,comm2)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawrhoij_mpisum_unpacked_1D'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: comm1
 integer,optional,intent(in) :: comm2
!arrays
 type(pawrhoij_type),intent(inout) :: pawrhoij(:)

!Local variables ---------------------------------------
!scalars
 integer :: bufdim,iatom,ierr,isppol,jdim,nsp2,natom
 integer :: nproc1,nproc2
 !character(len=500) :: msg
!arrays
 integer,allocatable :: dimlmn(:)
 real(dp),allocatable :: buffer(:)

!************************************************************************

 natom=SIZE(pawrhoij);if (natom==0) return

 nproc1 = xmpi_comm_size(comm1)
 nproc2=1; if (PRESENT(comm2)) nproc2 = xmpi_comm_size(comm2)
 if (nproc1==1.and.nproc2==1) RETURN

!Fill the MPI buffer from the local rhoij_
 LIBPAW_ALLOCATE(dimlmn,(natom))
 dimlmn(1:natom)=pawrhoij(1:natom)%cplex*pawrhoij(1:natom)%lmn2_size
 nsp2=pawrhoij(1)%nsppol; if (pawrhoij(1)%nspden==4) nsp2=4
 bufdim=sum(dimlmn)*nsp2
 LIBPAW_ALLOCATE(buffer,(bufdim))
 jdim=0
 do iatom=1,natom
   do isppol=1,nsp2
     buffer(jdim+1:jdim+dimlmn(iatom))=pawrhoij(iatom)%rhoij_(:,isppol)
     jdim=jdim+dimlmn(iatom)
   end do
 end do

!Build sum of pawrhoij%rhoij_
 call xmpi_sum(buffer,comm1,ierr)   ! Sum over the first communicator.
 if (PRESENT(comm2)) then
   call xmpi_sum(buffer,comm2,ierr) ! Sum over the second communicator.
 end if

!Unpack the MPI packet filling rhoij_
 jdim=0
 do iatom=1,natom
   do isppol=1,nsp2
     pawrhoij(iatom)%rhoij_(:,isppol)=buffer(jdim+1:jdim+dimlmn(iatom))
     jdim=jdim+dimlmn(iatom)
   end do
 end do

 LIBPAW_DEALLOCATE(buffer)
 LIBPAW_DEALLOCATE(dimlmn)

end subroutine pawrhoij_mpisum_unpacked_1D
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/pawrhoij_mpisum_unpacked_2D
!! NAME
!! pawrhoij_mpisum_unpacked_2D
!!
!! FUNCTION
!! Build the MPI sum of the unsymmetrized PAW rhoij_ (augmentation occupancies)
!! Remember:for each atom, rho_ij=Sum_{n,k} {occ(n,k)*<Cnk|p_i><p_j|Cnk>}
!! Target: 2D array of pawrhoij datastructures
!!
!! INPUTS
!!  comm1=MPI communicator. Data will be MPI summed inside comm1
!!  [comm2]=second MPI communicator. If present, rhoij_ will be
!!          MPI summed inside comm2 after the collective sum in comm1.
!!
!! SIDE EFFECTS
!!  pawrhoij(:,:) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  Input: the data calculateed by this processor.
!!  Output: the final MPI sum over comm1 and comm2.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawrhoij_mpisum_unpacked_2D(pawrhoij,comm1,comm2)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawrhoij_mpisum_unpacked_2D'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: comm1
 integer,optional,intent(in) :: comm2
!arrays
 type(pawrhoij_type),intent(inout) :: pawrhoij(:,:)

!Local variables ---------------------------------------
!scalars
 integer :: bufdim,iatom,ierr,irhoij,isppol,jdim,nsp2,natom,nrhoij
 integer :: nproc1,nproc2
 !character(len=500) :: msg
!arrays
 integer,allocatable :: dimlmn(:,:)
 real(dp),allocatable :: buffer(:)

!************************************************************************

 natom =SIZE(pawrhoij,1);if (natom ==0) return
 nrhoij=SIZE(pawrhoij,2);if (nrhoij==0) return

 nproc1 = xmpi_comm_size(comm1)
 nproc2=1; if (PRESENT(comm2)) nproc2 = xmpi_comm_size(comm2)
 if (nproc1==1.and.nproc2==1) RETURN

!Fill the MPI buffer from the local rhoij_
 LIBPAW_ALLOCATE(dimlmn,(natom,nrhoij))
 dimlmn(1:natom,1:nrhoij)=pawrhoij(1:natom,1:nrhoij)%cplex*pawrhoij(1:natom,1:nrhoij)%lmn2_size
 nsp2=pawrhoij(1,1)%nsppol; if (pawrhoij(1,1)%nspden==4) nsp2=4
 bufdim=sum(dimlmn)*nsp2
 LIBPAW_ALLOCATE(buffer,(bufdim))
 jdim=0
 do irhoij=1,nrhoij
   do iatom=1,natom
     do isppol=1,nsp2
       buffer(jdim+1:jdim+dimlmn(iatom,irhoij))=pawrhoij(iatom,irhoij)%rhoij_(:,isppol)
       jdim=jdim+dimlmn(iatom,irhoij)
     end do
   end do
 end do

!Build sum of pawrhoij%rhoij_
 call xmpi_sum(buffer,comm1,ierr)   ! Sum over the first communicator.
 if (PRESENT(comm2)) then
   call xmpi_sum(buffer,comm2,ierr) ! Sum over the second communicator.
 end if

!Unpack the MPI packet filling rhoij_
 jdim=0
 do irhoij=1,nrhoij
   do iatom=1,natom
     do isppol=1,nsp2
       pawrhoij(iatom,irhoij)%rhoij_(:,isppol)=buffer(jdim+1:jdim+dimlmn(iatom,irhoij))
       jdim=jdim+dimlmn(iatom,irhoij)
     end do
   end do
 end do

 LIBPAW_DEALLOCATE(buffer)
 LIBPAW_DEALLOCATE(dimlmn)

end subroutine pawrhoij_mpisum_unpacked_2D
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/pawrhoij_get_nspden
!! NAME
!! pawrhoij_get_nspden
!!
!! FUNCTION
!! Compute the value of nspden (number of spin component) for a pawrhoij datastructure.
!! This one depends on the number of spinorial components, the number of components
!! of the density and use of the spin-orbit coupling.
!!
!! INPUTS
!!  nspden= number of spin-density components
!!  nspinor= number of spinorial components
!!  spnorb= flag: 1 if spin-orbit coupling is activated
!!
!! OUTPUT
!!  pawrhoij_get_nspden= value of nspden associated to pawrhoij
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function pawrhoij_get_nspden(nspden,nspinor,spnorb)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawrhoij_get_nspden'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: nspden,nspinor,spnorb
 integer :: pawrhoij_get_nspden
!arrays

!Local variables ---------------------------------------

!************************************************************************

 pawrhoij_get_nspden=nspden
 if (nspinor==2.and.spnorb>0) pawrhoij_get_nspden=4

end function pawrhoij_get_nspden
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/symrhoij
!! NAME
!! symrhoij
!!
!! FUNCTION
!! Symmetrize rhoij quantities (augmentation occupancies) and/or gradients
!! Compute also rhoij residuals (new-old values of rhoij and gradients)
!!
!! INPUTS
!!  choice=select then type of rhoij gradients to symmetrize.
!!         choice=1 => no gradient
!!         choice=2 => gradient with respect to atomic position(s)
!!               =3 => a gradient with respect to strain(s)
!!               =4 => 2nd gradient with respect to atomic position(s)
!!               =23=> a gradient with respect to atm. pos. and strain(s)
!!               =24=> 1st and 2nd gradient with respect to atomic position(s)
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space(bohr^-1).
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  ipert=index of perturbation if pawrhoij is a pertubed rhoij
!!        no meaning for ground-state calculations (should be 0)
!!  [mpi_atmtab(:)]=--optional-- indexes of the atoms treated by current proc
!!  [comm_atom]=--optional-- MPI communicator over atoms
!!  natom=number of atoms in cell
!!  nsym=number of symmetry elements in space group
!!  ntypat=number of types of atoms in unit cell.
!!  optrhoij= 1 if rhoij quantities have to be symmetrized
!!  pawrhoij_unsym(:) <type(pawrhoij_type)>= datastructure containing PAW rhoij occupancies
!!    (and related data) non symmetrized in an unpacked storage (pawrhoij_unsym%rhoij_)
!!    Contains eventually unsymmetrized rhoij gradients (grhoij)
!!  pawang <type(pawang_type)>=angular mesh discretization and related data
!!  pawprtvol=control print volume and debugging output for PAW
!!            Note: if pawprtvol=-10001, nothing is printed out
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  [qphon(3)]=--optional-- (RF calculations only) - wavevector of the phonon
!!  rprimd(3,3)=real space primitive translations.
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrec(3,3,nsym)=symmetries of group in terms of operations on
!!                   reciprocal space primitive translations
!!  typat(natom)=type for each atom
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  pawrhoij(natom) <type(pawrhoij_type)>= datastructure containing PAW rhoij occupancies
!!    (and related data) SYMMETRIZED in a packed storage (pawrhoij%rhoijp).
!!    if (optrhoij==1)
!!      pawrhoij(natom)%nrhoijsel=number of non-zero values of rhoij
!!      pawrhoij(natom)%rhoijp(cplex*lmn2_size,nspden)=symetrized paw rhoij quantities in PACKED STORAGE (only non-zero values)
!!      pawrhoij(natom)%rhoijres(cplex*lmn2_size,nspden)=paw rhoij quantities residuals (new values - old values)
!!      pawrhoij(natom)%rhoijselect(lmn2_size)=select the non-zero values of rhoij!!
!!    if (pawrhoij(:)%ngrhoij>0) (equivalent to choice>1)
!!      pawrhoij(natom)%grhoij(ngrhoij,cplex*lmn2_size,nspden)=symetrized gradients of rhoij
!!
!! NOTES
!!  pawrhoij and pawrhoij_unsym can be identical (refer to the same pawrhoij datastructure).
!!  They should be different only if pawrhoij is distributed over atomic sites
!!  (in that case pawrhoij_unsym should not be distributed over atomic sites).
!!
!! PARENTS
!!      d2frnl,energy,paw_qpscgw,pawmkrho,posdoppler
!!
!! CHILDREN
!!
!! SOURCE

subroutine symrhoij(pawrhoij,pawrhoij_unsym,choice,gprimd,indsym,ipert,natom,nsym,&
&                   ntypat,optrhoij,pawang,pawprtvol,pawtab,rprimd,symafm,symrec,typat, &
&                   mpi_atmtab,comm_atom,qphon) ! optional arguments (parallelism)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symrhoij'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: choice,ipert,natom,nsym,ntypat,optrhoij,pawprtvol
 integer,optional,intent(in) :: comm_atom
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: indsym(4,nsym,natom)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 integer,intent(in) :: symafm(nsym),symrec(3,3,nsym),typat(natom)
 real(dp),intent(in) :: gprimd(3,3),rprimd(3,3)
 real(dp),intent(in),optional :: qphon(3)
 type(pawrhoij_type),intent(inout) :: pawrhoij(:)
 type(pawrhoij_type),target,intent(inout) :: pawrhoij_unsym(:)
 type(pawtab_type),target,intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
 character(len=8),parameter :: dspin(6)=(/"up      ","down    ","dens (n)","magn (x)","magn (y)","magn (z)"/)
!scalars
 integer :: at_indx,cplex,cplex_eff,iafm,iatm,iatom,idum1,il,il0,ilmn,iln,iln0,ilpm,indexi
 integer :: indexii,indexj,indexjj,indexjj0,indexk,indexk1,iplex,irhoij,irot,ishift2,ishift3
 integer :: ishift4,ispden,itypat,j0lmn,jj,jl,jl0,jlmn,jln,jln0,jlpm,jrhoij,jspden,klmn,klmn1,kspden
 integer :: lmn_size,mi,mj,my_comm_atom,mu,mua,mub,mushift,natinc,ngrhoij,nrhoij,nrhoij1,nrhoij_unsym
 integer :: nselect,nselect1,nspinor,nu,nushift,sz1,sz2
 logical,parameter :: afm_noncoll=.true.  ! TRUE if antiferro symmetries are used with non-collinear magnetism
 real(dp) :: arg,factafm,syma,zarot2
 logical :: antiferro,have_phase,my_atmtab_allocated,noncoll,paral_atom,paral_atom_unsym,use_afm,use_res
 character(len=8) :: pertstrg,wrt_mode
 character(len=500) :: msg
!arrays
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 integer :: nsym_used(2)
 integer, pointer :: indlmn(:,:)
 integer,pointer :: my_atmtab(:)
 integer :: idum(0)
 real(dp) :: factsym(2),phase(2),rhoijc(2),ro(2),sumrho(2,2),sum1(2),rotrho(2,2),xsym(3)
 real(dp),allocatable :: rotgr(:,:,:),rotmag(:,:),rotmaggr(:,:,:)
 real(dp),allocatable :: sumgr(:,:),summag(:,:),summaggr(:,:,:),symrec_cart(:,:,:),work1(:,:,:)
 real(dp),pointer :: grad(:,:,:)
 type(coeff3_type),target,allocatable :: tmp_grhoij(:)
 type(pawrhoij_type),pointer :: pawrhoij_unsym_all(:)

! *********************************************************************

!Sizes of pawrhoij datastructures
 nrhoij=size(pawrhoij)
 nrhoij_unsym=size(pawrhoij_unsym)

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(nrhoij/=natom))
 paral_atom_unsym=(present(comm_atom).and.(nrhoij_unsym/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom)

!Symetrization occurs only when nsym>1
 if (nsym>1) then

!  Test: consistency between choice and ngrhoij
   ngrhoij=0
   if (nrhoij>0) then
     ngrhoij=pawrhoij(1)%ngrhoij
     if(choice==2) ngrhoij=min(3,pawrhoij(1)%ngrhoij)
     if(choice==3.or.choice==4)   ngrhoij=min(6,pawrhoij(1)%ngrhoij)
     if(choice==23.or.choice==24) ngrhoij=min(9,pawrhoij(1)%ngrhoij)
     if ((choice==1.and.ngrhoij/=0).or.(choice==2.and.ngrhoij/=3).or. &
&     (choice==3.and.ngrhoij/=6).or.(choice==23.and.ngrhoij/=9).or. &
&     (choice==4.and.ngrhoij/=6).or.(choice==24.and.ngrhoij/=9) ) then
       msg='Inconsistency between variables choice and ngrhoij !'
       MSG_BUG(msg)
     end if
   end if

!  Symetrization of gradients not compatible with nspden=4
   if (nrhoij>0) then
     if (choice>2.and.pawrhoij(1)%nspden==4) then
       msg='For the time being, choice>2 is not compatible with nspden=4 !'
       MSG_BUG(msg)
     end if
   end if

!  Symetry matrixes must be in memory
   if (pawang%nsym==0) then
     msg='pawang%zarot must be allocated !'
     MSG_BUG(msg)
   end if

!  Antiferro case ?
   antiferro=.false.;if (nrhoij>0) antiferro=(pawrhoij(1)%nspden==2.and.pawrhoij(1)%nsppol==1)
!  Non-collinear case
   noncoll=.false.;if (nrhoij>0) noncoll=(pawrhoij(1)%nspden==4)
!  Do we use antiferro symmetries ?
   use_afm=((antiferro).or.(noncoll.and.afm_noncoll))

   cplex_eff=1
   if (nrhoij>0.and.(ipert>0.or.antiferro.or.noncoll)) cplex_eff=pawrhoij(1)%cplex ! Does not symmetrize imaginary part for GS calculations

!  Do we have a phase due to q-vector (phonons only) ?
   have_phase=.false.
   if (ipert>0.and.present(qphon).and.nrhoij>0) then
     have_phase=(abs(qphon(1))>tol8.or.abs(qphon(2))>tol8.or.abs(qphon(3))>tol8)
     if (choice>1) then
       msg='choice>1 not compatible with q-phase !'
       MSG_BUG(msg)
     end if
     if (have_phase.and.cplex_eff==1) then
       msg='Should have cplex_dij=2 for a non-zero q!'
       MSG_BUG(msg)
     end if
   end if

!  Several inits/allocations
   if (noncoll.and.optrhoij==1)  then
     LIBPAW_ALLOCATE(summag,(cplex_eff,3))
     LIBPAW_ALLOCATE(rotmag,(cplex_eff,3))
   end if
   if (noncoll) then
     LIBPAW_ALLOCATE(symrec_cart,(3,3,nsym))
     do irot=1,nsym
       symrec_cart(:,:,irot)=symrhoij_symcart(gprimd,rprimd,symrec(:,:,irot))
     end do
   end if
   ishift2=0;ishift3=0;ishift4=0
   if (choice>1) then
     LIBPAW_ALLOCATE(sumgr,(cplex_eff,ngrhoij))
     if (choice>2)  then
       LIBPAW_ALLOCATE(work1,(cplex_eff,3,3))
     end if
     if (antiferro) then
       LIBPAW_ALLOCATE(rotgr,(cplex_eff,ngrhoij,2))
     else
       LIBPAW_ALLOCATE(rotgr,(cplex_eff,ngrhoij,1))
     end if
     if (noncoll) then
       LIBPAW_ALLOCATE(summaggr,(cplex_eff,ngrhoij,3))
       LIBPAW_ALLOCATE(rotmaggr,(cplex_eff,ngrhoij,3))
     end if
     if (choice==23) ishift2=6
     if (choice==24) ishift4=3
     if (.not.paral_atom_unsym) then
!      Have to make a temporary copy of grhoij
       LIBPAW_DATATYPE_ALLOCATE(tmp_grhoij,(nrhoij))
       do iatm=1,nrhoij
         sz1=pawrhoij_unsym(iatm)%cplex*pawrhoij_unsym(iatm)%lmn2_size
         sz2=pawrhoij_unsym(iatm)%nspden
         LIBPAW_ALLOCATE(tmp_grhoij(iatm)%value,(ngrhoij,sz1,sz2))
         tmp_grhoij(iatm)%value(1:ngrhoij,1:sz1,1:sz2)=pawrhoij_unsym(iatm)%grhoij(1:ngrhoij,1:sz1,1:sz2)
       end do
     end if
   end if

!  In case of paw_rhoij_unsym distributed over atomic sites, gather it
   if (paral_atom_unsym) then
     LIBPAW_DATATYPE_ALLOCATE(pawrhoij_unsym_all,(natom))
     call pawrhoij_nullify(pawrhoij_unsym_all)
     call pawrhoij_gather(pawrhoij_unsym,pawrhoij_unsym_all,-1,my_comm_atom,&
&     with_lmnmix=.false.,with_rhoijp=.false.,&
&     with_rhoijres=.false.,with_rhoijim=.false.,with_grhoij=(choice>1))
     nrhoij1=natom
   else
     pawrhoij_unsym_all=>pawrhoij_unsym
     nrhoij1=nrhoij_unsym
   end if

!  Loops over atoms and spin components
!  ------------------------------------
   do iatm=1,nrhoij
     iatom=iatm;if (paral_atom) iatom=my_atmtab(iatm)
     if (nrhoij==1.and.ipert>0.and.ipert<=natom.and.(.not.paral_atom)) iatom=ipert
     itypat=typat(iatom)
     lmn_size=pawrhoij(iatm)%lmn_size
     nspinor=pawrhoij(iatm)%nspinor
     cplex=pawrhoij(iatm)%cplex
     cplex_eff=1;if (ipert>0.or.antiferro.or.noncoll) cplex_eff=cplex
     use_res=(pawrhoij(iatm)%use_rhoijres>0)
     indlmn => pawtab(itypat)%indlmn

     nselect=0
     do ispden=1,pawrhoij(iatm)%nsppol
       jspden=min(3-ispden,pawrhoij(iatm)%nsppol)

!      Store old -rhoij in residual
       if (optrhoij==1.and.use_res) then
         pawrhoij(iatm)%rhoijres(:,ispden)=zero
         if (cplex==1) then
           do irhoij=1,pawrhoij(iatm)%nrhoijsel
             klmn=pawrhoij(iatm)%rhoijselect(irhoij)
             pawrhoij(iatm)%rhoijres(klmn,ispden)=-pawrhoij(iatm)%rhoijp(irhoij,ispden)
           end do
         else
           do irhoij=1,pawrhoij(iatm)%nrhoijsel
             klmn1=2*pawrhoij(iatm)%rhoijselect(irhoij);jrhoij=2*irhoij
             pawrhoij(iatm)%rhoijres(klmn1-1:klmn1,ispden)=-pawrhoij(iatm)%rhoijp(jrhoij-1:jrhoij,ispden)
           end do
         end if
         if (noncoll) then
           pawrhoij(iatm)%rhoijres(:,2:4)=zero
           if (cplex==1) then
             do mu=2,4
               do irhoij=1,pawrhoij(iatm)%nrhoijsel
                 klmn=pawrhoij(iatm)%rhoijselect(irhoij)
                 pawrhoij(iatm)%rhoijres(klmn,mu)=-pawrhoij(iatm)%rhoijp(irhoij,mu)
               end do
             end do
           else
             do mu=2,4
               do irhoij=1,pawrhoij(iatm)%nrhoijsel
                 klmn1=2*pawrhoij(iatm)%rhoijselect(irhoij);jrhoij=2*irhoij
                 pawrhoij(iatm)%rhoijres(klmn1-1:klmn1,mu)=-pawrhoij(iatm)%rhoijp(jrhoij-1:jrhoij,mu)
               end do
             end do
           end if
         end if
         if (antiferro) then
           pawrhoij(iatm)%rhoijres(:,2)=zero
           if (cplex==1) then
             do irhoij=1,pawrhoij(iatm)%nrhoijsel
               klmn=pawrhoij(iatm)%rhoijselect(irhoij)
               pawrhoij(iatm)%rhoijres(klmn,2)=-pawrhoij(iatm)%rhoijp(irhoij,2)
             end do
           else
             do irhoij=1,pawrhoij(iatm)%nrhoijsel
               klmn1=2*pawrhoij(iatm)%rhoijselect(irhoij);jrhoij=2*irhoij
               pawrhoij(iatm)%rhoijres(klmn1-1:klmn1,2)=-pawrhoij(iatm)%rhoijp(jrhoij-1:jrhoij,2)
             end do
           end if
         end if
       end if

!      Loops over (il,im) and (jl,jm)
!      ------------------------------
       jl0=-1;jln0=-1;indexj=1
       do jlmn=1,lmn_size
         jl=indlmn(1,jlmn)
         jlpm=1+jl+indlmn(2,jlmn)
         jln=indlmn(5,jlmn)
         if (jln/=jln0) indexj=indexj+2*jl0+1
         j0lmn=jlmn*(jlmn-1)/2
         il0=-1;iln0=-1;indexi=1
         do ilmn=1,jlmn
           il=indlmn(1,ilmn)
           ilpm=1+il+indlmn(2,ilmn)
           iln=indlmn(5,ilmn)
           if (iln/=iln0) indexi=indexi+2*il0+1
           klmn=j0lmn+ilmn;klmn1=cplex*klmn

           nsym_used(:)=0
           if (optrhoij==1) rotrho(:,:)=zero
           if (optrhoij==1.and.noncoll) rotmag(:,:)=zero
           if (choice>1) rotgr(:,:,:)=zero
           if (choice>1.and.noncoll) rotmaggr(:,:,:)=zero

!          Loop over symmetries
!          --------------------
           do irot=1,nsym

             if ((symafm(irot)/=1).and.(.not.use_afm)) cycle
             kspden=ispden;if (symafm(irot)==-1) kspden=jspden
             iafm=1;if ((antiferro).and.(symafm(irot)==-1)) iafm=2
             factafm=dble(symafm(irot))

             nsym_used(iafm)=nsym_used(iafm)+1
             at_indx=min(indsym(4,irot,iatom),nrhoij1)

             if (have_phase) then
               arg=two_pi*(qphon(1)*indsym(1,irot,iatom)+qphon(2)*indsym(2,irot,iatom) &
&                         +qphon(3)*indsym(3,irot,iatom))
               phase(1)=cos(arg);phase(2)=sin(arg)
             end if

!            Accumulate values over (mi,mj)
!            ------------------------------
             if (optrhoij==1) sumrho(:,:)=zero
             if (optrhoij==1.and.noncoll) summag(:,:)=zero
             if (choice>1) sumgr(:,:)=zero
             if (choice>1.and.noncoll) summaggr(:,:,:)=zero
             do mj=1,2*jl+1
               indexjj=indexj+mj;indexjj0=indexjj*(indexjj-1)/2
               do mi=1,2*il+1
                 factsym(:)=one
                 indexii=indexi+mi
                 if (indexii<=indexjj) then
                   indexk=indexjj0+indexii
                   if(cplex_eff==2.and.nspinor==2) factsym(2)=one
                 else
                   indexk=indexii*(indexii-1)/2+indexjj
                   if(cplex_eff==2.and.nspinor==2) factsym(2)=-one
                 end if
!                Be careful: use here R_rel^-1 in term of spherical harmonics
!                which is tR_rec in term of spherical harmonics
!                so, use transpose[zarot]
                 zarot2=pawang%zarot(mi,ilpm,il+1,irot)*pawang%zarot(mj,jlpm,jl+1,irot)
!                zarot2=pawang%zarot(ilpm,mi,il+1,irot)*pawang%zarot(jlpm,mj,jl+1,irot)

!                Rotate rhoij
                 if (optrhoij==1) then
                   if (cplex==1) then
                     sumrho(1,iafm)=sumrho(1,iafm)+zarot2*pawrhoij_unsym_all(at_indx)%rhoij_(indexk,kspden)
                   else
                     indexk1=2*(indexk-1)
                     sumrho(1,iafm)=sumrho(1,iafm)+factsym(1)*zarot2*pawrhoij_unsym_all(at_indx)%rhoij_(indexk1+1,kspden)
                     if(cplex_eff==2) sumrho(2,iafm)= &
&                           sumrho(2,iafm)+factsym(2)*factafm*zarot2*pawrhoij_unsym_all(at_indx)%rhoij_(indexk1+2,kspden)
                   end if
                 end if

!                If non-collinear case, rotate rhoij magnetization
                 if (optrhoij==1.and.noncoll) then
                   if (cplex==1) then
                     do mu=1,3
                       summag(1,mu)=summag(1,mu)+zarot2*factafm*pawrhoij_unsym_all(at_indx)%rhoij_(indexk,1+mu)
                     end do
                   else
                     indexk1=2*(indexk-1)
                     do mu=1,3
                       summag(1,mu)=summag(1,mu) &
&                       +zarot2*factsym(1)*factafm*pawrhoij_unsym_all(at_indx)%rhoij_(indexk1+1,1+mu)
                       if(cplex_eff==2) summag(2,mu)=summag(2,mu)&
&                       +zarot2*factsym(2)*pawrhoij_unsym_all(at_indx)%rhoij_(indexk1+2,1+mu)
                     end do
                   end if
                 end if

!                Rotate gradients of rhoij
                 if (choice>1) then
                   if (paral_atom_unsym) then
                     grad => pawrhoij_unsym_all(at_indx)%grhoij
                   else
                     grad => tmp_grhoij(at_indx)%value
                   end if
                   if (cplex==1) then
                     do mu=1,ngrhoij
                       sumgr(1,mu)=sumgr(1,mu)+zarot2*grad(mu,indexk,kspden)
                     end do
                     if (noncoll) then
                       do mu=1,3
                         do nu=1,ngrhoij
                           summaggr(1,nu,mu)=summaggr(1,nu,mu)+zarot2*factafm*grad(nu,indexk,1+mu)
                         end do
                       end do
                     end if
                   else
                     indexk1=2*(indexk-1)
                     do mu=1,ngrhoij
                       sumgr(1:cplex_eff,mu)=sumgr(1:cplex_eff,mu) &
&                       +zarot2*factsym(1:cplex_eff)*grad(mu,indexk1+1:indexk1+cplex_eff,kspden)
                     end do
                     if (noncoll) then
                       do mu=1,3
                         do nu=1,ngrhoij
                           summaggr(1:cplex_eff,nu,mu)=summaggr(1:cplex_eff,nu,mu) &
&                           +zarot2*factsym(1:cplex_eff)*factafm*grad(nu,indexk1+1:indexk1+cplex_eff,1+mu)
                         end do
                       end do
                     end if
                   end if
                 end if

               end do
             end do

!            Apply phase for phonons
             if (have_phase) then
               if(optrhoij==1) then
                 rhoijc(1:2)=sumrho(1:2,iafm)
                 sumrho(1,iafm)=phase(1)*rhoijc(1)-phase(2)*rhoijc(2)
                 sumrho(2,iafm)=phase(1)*rhoijc(2)+phase(2)*rhoijc(1)
               end if
               if (noncoll) then
                 do mu=1,3
                   rhoijc(1:2)=summag(1:2,mu)
                   summag(1,mu)=phase(1)*rhoijc(1)-phase(2)*rhoijc(2)
                   summag(2,mu)=phase(1)*rhoijc(2)+phase(2)*rhoijc(1)
                 end do
               end if
               if (choice>1) then
                 do mu=1,ngrhoij
                   rhoijc(1:2)=sumgr(1:2,mu)
                   sumgr(1,mu)=phase(1)*rhoijc(1)-phase(2)*rhoijc(2)
                   sumgr(2,mu)=phase(1)*rhoijc(2)+phase(2)*rhoijc(1)
                 end do
                 if (noncoll) then
                   do mu=1,3
                     do nu=1,ngrhoij
                       rhoijc(1:2)=summaggr(1:2,nu,mu)
                       summaggr(1,nu,mu)=phase(1)*rhoijc(1)-phase(2)*rhoijc(2)
                       summaggr(2,nu,mu)=phase(1)*rhoijc(2)+phase(2)*rhoijc(1)
                     end do
                   end do
                 end if
               end if
             end if

!            Add contribution of this rotation
             if (optrhoij==1) then
               rotrho(1:cplex_eff,iafm)=rotrho(1:cplex_eff,iafm)+sumrho(1:cplex_eff,iafm)
             end if

!            Rotate vector fields in real space (forces, magnetization, etc...)
!            Should use symrel^-1 but use transpose[symrec] instead
!            ---------------------------------
!            ===== Rhoij magnetization ====
             if (noncoll.and.optrhoij==1) then
               do nu=1,3
                 do mu=1,3
                   rotmag(1:cplex_eff,mu)=rotmag(1:cplex_eff,mu)+symrec_cart(mu,nu,irot)*summag(1:cplex_eff,nu)
                 end do
               end do
             end if
!            ===== Derivatives vs atomic positions ====
             if (choice==2.or.choice==23.or.choice==24) then
               do nu=1,3
                 nushift=nu+ishift2
                 do mu=1,3
                   mushift=mu+ishift2
                   rotgr(1:cplex_eff,mushift,iafm)=&
&                   rotgr(1:cplex_eff,mushift,iafm)+dble(symrec(mu,nu,irot))*sumgr(1:cplex_eff,nushift)
                 end do
               end do
               if (noncoll) then
                 do mub=1,3 ! Loop on magnetization components
                   do mua=1,3 ! Loop on gradients
                     mushift=mua+ishift2
                     sum1(:)=zero;xsym(1:3)=dble(symrec(mua,1:3,irot))
                     do nu=1,3
                       syma=symrec_cart(mub,nu,irot)
                       sum1(1:cplex_eff)=sum1(1:cplex_eff)+syma*(summaggr(1:cplex_eff,ishift2+1,nu)*xsym(1) &
&                       +summaggr(1:cplex_eff,ishift2+2,nu)*xsym(2) &
&                       +summaggr(1:cplex_eff,ishift2+3,nu)*xsym(3))
                     end do
                     rotmaggr(1:cplex_eff,mushift,mub)=rotmaggr(1:cplex_eff,mushift,mub)+sum1(1:cplex_eff)
                   end do
                 end do
               end if
             end if
!            ===== Derivatives vs strain ====
             if (choice==3.or.choice==23) then
               work1(1:cplex_eff,1,1)=sumgr(1:cplex_eff,1+ishift3);work1(1:cplex_eff,2,2)=sumgr(1:cplex_eff,2+ishift3)
               work1(1:cplex_eff,3,3)=sumgr(1:cplex_eff,3+ishift3);work1(1:cplex_eff,2,3)=sumgr(1:cplex_eff,4+ishift3)
               work1(1:cplex_eff,1,3)=sumgr(1:cplex_eff,5+ishift3);work1(1:cplex_eff,1,2)=sumgr(1:cplex_eff,6+ishift3)
               work1(1:cplex_eff,3,1)=work1(1:cplex_eff,1,3);work1(1:cplex_eff,3,2)=work1(1:cplex_eff,2,3)
               work1(1:cplex_eff,2,1)=work1(1:cplex_eff,1,2)
               do mu=1,6
                 mushift=mu+ishift3
                 mua=alpha(mu);mub=beta(mu)
                 sum1(:)=zero;xsym(1:3)=dble(symrec(mub,1:3,irot))
                 do nu=1,3
                   syma=dble(symrec(mua,nu,irot))
                   sum1(1:cplex_eff)=sum1(1:cplex_eff)+syma*(work1(1:cplex_eff,nu,1)*xsym(1) &
&                   +work1(1:cplex_eff,nu,2)*xsym(2) &
&                   +work1(1:cplex_eff,nu,3)*xsym(3))
                 end do
                 rotgr(1:cplex_eff,mushift,iafm)=rotgr(1:cplex_eff,mushift,iafm)+sum1(1:cplex_eff)
               end do
             end if
!            ===== Second derivatives vs atomic positions ====
             if (choice==4.or.choice==24) then
               work1(1:cplex_eff,1,1)=sumgr(1:cplex_eff,1+ishift4);work1(1:cplex_eff,2,2)=sumgr(1:cplex_eff,2+ishift4)
               work1(1:cplex_eff,3,3)=sumgr(1:cplex_eff,3+ishift4);work1(1:cplex_eff,2,3)=sumgr(1:cplex_eff,4+ishift4)
               work1(1:cplex_eff,1,3)=sumgr(1:cplex_eff,5+ishift4);work1(1:cplex_eff,1,2)=sumgr(1:cplex_eff,6+ishift4)
               work1(1:cplex_eff,3,1)=work1(1:cplex_eff,1,3);work1(1:cplex_eff,3,2)=work1(1:cplex_eff,2,3)
               work1(1:cplex_eff,2,1)=work1(1:cplex_eff,1,2)
               do mu=1,6
                 mushift=mu+ishift4
                 mua=alpha(mu);mub=beta(mu)
                 sum1(:)=zero
                 xsym(1:3)=dble(symrec(mub,1:3,irot))
                 do nu=1,3
                   syma=dble(symrec(mua,nu,irot))
                   sum1(1:cplex_eff)=sum1(1:cplex_eff)+syma*(work1(1:cplex_eff,nu,1)*xsym(1) &
&                   +work1(1:cplex_eff,nu,2)*xsym(2) &
&                   +work1(1:cplex_eff,nu,3)*xsym(3))
                 end do
                 rotgr(1:cplex_eff,mushift,iafm)=rotgr(1:cplex_eff,mushift,iafm)+sum1(1:cplex_eff)
               end do
             end if

           end do ! End loop over symmetries

!          Store average result (over symmetries)
!          --------------------------------------
           if (optrhoij==1) then

!            Mean value for rhoij
             if (cplex==1) then
               ro(1)=rotrho(1,1)/nsym_used(1)
               if (abs(ro(1))>tol10) then
                 pawrhoij(iatm)%rhoijp(klmn,ispden)=ro(1)
                 if (use_res) pawrhoij(iatm)%rhoijres(klmn,ispden)=pawrhoij(iatm)%rhoijres(klmn,ispden)+ro(1)
               else
                 pawrhoij(iatm)%rhoijp(klmn,ispden)=zero
               end if
             else
               ro(1)=rotrho(1,1)/nsym_used(1)
               if (cplex_eff==2) then
                 ro(2)=rotrho(2,1)/nsym_used(1)
               else
                 ro(2)=pawrhoij_unsym_all(iatom)%rhoij_(klmn1,ispden)
               end if
               if (any(abs(ro(1:2))>tol10)) then
                 pawrhoij(iatm)%rhoijp(klmn1-1,ispden)=ro(1)
                 pawrhoij(iatm)%rhoijp(klmn1  ,ispden)=ro(2)
                 if (use_res) then
                   pawrhoij(iatm)%rhoijres(klmn1-1,ispden)=pawrhoij(iatm)%rhoijres(klmn1-1,ispden)+ro(1)
                   pawrhoij(iatm)%rhoijres(klmn1  ,ispden)=pawrhoij(iatm)%rhoijres(klmn1  ,ispden)+ro(2)
                 end if
               else
                 pawrhoij(iatm)%rhoijp(klmn1-1,ispden)=zero
                 pawrhoij(iatm)%rhoijp(klmn1  ,ispden)=zero
               end if
             end if

!            Non-collinear case: mean value for rhoij magnetization
             if (noncoll) then
!              Select on-zero elements
               if (cplex==1) then
                 do mu=2,4
                   ro(1)=rotmag(1,mu-1)/nsym_used(1)
                   if (abs(ro(1))>tol10) then
                     pawrhoij(iatm)%rhoijp(klmn,mu)=ro(1)
                     if (use_res) pawrhoij(iatm)%rhoijres(klmn,mu)=pawrhoij(iatm)%rhoijres(klmn,mu)+ro(1)
                   else
                     pawrhoij(iatm)%rhoijp(klmn,mu)=zero
                   end if
                 end do
               else
                 do mu=2,4
                   ro(1)=rotmag(1,mu-1)/nsym_used(1)
                   if (cplex_eff==2) then
                     ro(2)=rotmag(2,mu-1)/nsym_used(1)
                   else
                     ro(2)=pawrhoij_unsym_all(iatom)%rhoij_(klmn1,mu)
                   end if
                   if (any(abs(ro(1:2))>tol10)) then
                     pawrhoij(iatm)%rhoijp(klmn1-1,mu)=ro(1)
                     pawrhoij(iatm)%rhoijp(klmn1  ,mu)=ro(2)
                     if (use_res) then
                       pawrhoij(iatm)%rhoijres(klmn1-1,mu)=pawrhoij(iatm)%rhoijres(klmn1-1,mu)+ro(1)
                       pawrhoij(iatm)%rhoijres(klmn1  ,mu)=pawrhoij(iatm)%rhoijres(klmn1  ,mu)+ro(2)
                     end if
                   else
                     pawrhoij(iatm)%rhoijp(klmn1-1,mu)=zero
                     pawrhoij(iatm)%rhoijp(klmn1  ,mu)=zero
                   end if
                 end do
               end if
             end if

!            Antiferro case: mean value for down component
             if (antiferro.and.nsym_used(2)>0) then
               if (cplex==1) then
                 ro(1)=rotrho(1,2)/nsym_used(2)
                 if (abs(ro(1))>tol10) then
                   pawrhoij(iatm)%rhoijp(klmn,2)=ro(1)
                   if (use_res) pawrhoij(iatm)%rhoijres(klmn,2)=pawrhoij(iatm)%rhoijres(klmn,2)+ro(1)
                 else
                   pawrhoij(iatm)%rhoijp(klmn,2)=zero
                 end if
               else
                 ro(1:cplex_eff)=rotrho(1:cplex_eff,2)/nsym_used(2)
                 if (any(abs(ro(1:2))>tol10)) then
                   pawrhoij(iatm)%rhoijp(klmn1-1,2)=ro(1)
                   pawrhoij(iatm)%rhoijp(klmn1  ,2)=ro(2)
                   if (use_res) then
                     pawrhoij(iatm)%rhoijres(klmn1-1,2)=pawrhoij(iatm)%rhoijres(klmn1-1,2)+ro(1)
                     pawrhoij(iatm)%rhoijres(klmn1  ,2)=pawrhoij(iatm)%rhoijres(klmn1  ,2)+ro(2)
                   end if
                 else
                   pawrhoij(iatm)%rhoijp(klmn1-1,2)=zero
                   pawrhoij(iatm)%rhoijp(klmn1  ,2)=zero
                 end if
               end if
             end if

!            Select non-zero elements of rhoij
             if (ispden==pawrhoij(iatm)%nsppol) then
               if (cplex==1) then
                 if (any(abs(pawrhoij(iatm)%rhoijp(klmn,:))>tol10)) then
                   nselect=nselect+1
                   pawrhoij(iatm)%rhoijselect(nselect)=klmn
                   do jj=1,pawrhoij(iatm)%nspden
                     pawrhoij(iatm)%rhoijp(nselect,jj)=pawrhoij(iatm)%rhoijp(klmn,jj)
                   end do
                 end if
               else
                 if (any(abs(pawrhoij(iatm)%rhoijp(klmn1-1:klmn1,:))>tol10)) then
                   nselect=nselect+1;nselect1=2*nselect
                   pawrhoij(iatm)%rhoijselect(nselect)=klmn
                   do jj=1,pawrhoij(iatm)%nspden
                     pawrhoij(iatm)%rhoijp(nselect1-1,jj)=pawrhoij(iatm)%rhoijp(klmn1-1,jj)
                     pawrhoij(iatm)%rhoijp(nselect1  ,jj)=pawrhoij(iatm)%rhoijp(klmn1  ,jj)
                   end do
                 else if (pawrhoij(iatm)%use_rhoijim/=0) then ! for saving values of klmn which have non-zero rhoijim
                   if (any(abs(pawrhoij(iatm)%rhoijim(klmn1,:))>tol10)) then
                     nselect=nselect+1;nselect1=2*nselect
                     pawrhoij(iatm)%rhoijselect(nselect)=klmn
                     do jj=1,pawrhoij(iatm)%nspden
                       pawrhoij(iatm)%rhoijp(nselect1-1,jj)=pawrhoij(iatm)%rhoijp(klmn1-1,jj)
                       pawrhoij(iatm)%rhoijp(nselect1  ,jj)=pawrhoij(iatm)%rhoijp(klmn1  ,jj)
                     end do
                   end if
                 end if
               end if
             end if

           end if ! optrhoij==1

!          Store average result (over symmetries) for gradients
           if (choice>1) then
             do iplex=1,cplex_eff
               do mu=1,ngrhoij
                 pawrhoij(iatm)%grhoij(mu,klmn1+iplex-cplex,ispden)=rotgr(iplex,mu,1)/nsym_used(1)
               end do
             end do
             if (antiferro.and.nsym_used(2)>0) then
               do iplex=1,cplex_eff
                 do mu=1,ngrhoij
                   pawrhoij(iatm)%grhoij(mu,klmn1+iplex-cplex,2)=rotgr(iplex,mu,2)/nsym_used(2)
                 end do
               end do
             end if
             if (noncoll) then
               do nu=1,3
                 do iplex=1,cplex_eff
                   do mu=1,ngrhoij
                     pawrhoij(iatm)%grhoij(mu,klmn1+iplex-cplex,1+nu)=rotmaggr(iplex,mu,nu)/nsym_used(1)
                   end do
                 end do
               end do
             end if
           end if

           il0=il;iln0=iln  ! End loops over (il,im) and (jl,jm)
         end do
         jl0=jl;jln0=jln
       end do

     end do  ! End loop over ispden

!    Store number of non-zero values of rhoij
     if (optrhoij==1) pawrhoij(iatm)%nrhoijsel=nselect

   end do ! End loop over iatm

   if (noncoll.and.optrhoij==1)  then
     LIBPAW_DEALLOCATE(summag)
     LIBPAW_DEALLOCATE(rotmag)
   end if
   if (noncoll)  then
     LIBPAW_DEALLOCATE(symrec_cart)
   end if
   if (choice>1) then
     if (.not.paral_atom_unsym) then
       do iatm=1,nrhoij
         LIBPAW_DEALLOCATE(tmp_grhoij(iatm)%value)
       end do
       LIBPAW_DATATYPE_DEALLOCATE(tmp_grhoij)
     end if
     LIBPAW_DEALLOCATE(sumgr)
     LIBPAW_DEALLOCATE(rotgr)
     if (choice>2)  then
       LIBPAW_DEALLOCATE(work1)
     end if
     if (noncoll)  then
       LIBPAW_DEALLOCATE(summaggr)
       LIBPAW_DEALLOCATE(rotmaggr)
     end if
   end if
   if(paral_atom_unsym) then
     call pawrhoij_free(pawrhoij_unsym_all)
     LIBPAW_DATATYPE_DEALLOCATE(pawrhoij_unsym_all)
   end if


 else  ! nsym>1

!  *********************************************************************
!  If nsym==1, only copy rhoij_ into rhoij
!  also has to fill rhoijselect array

   if (nrhoij>0) then
     if(pawrhoij(1)%nspden==2.and.pawrhoij(1)%nsppol==1) then
       msg=' In the antiferromagnetic case, nsym cannot be 1'
       MSG_BUG(msg)
     end if
   end if
   if (optrhoij==1) then
     do iatm=1,nrhoij
       iatom=iatm;if ((paral_atom).and.(.not.paral_atom_unsym)) iatom=my_atmtab(iatm)
       cplex=pawrhoij(iatm)%cplex
       use_res=(pawrhoij(iatm)%use_rhoijres>0)
       if (use_res) then
         pawrhoij(iatm)%rhoijres(:,:)=zero
         if (cplex==1) then
           do ispden=1,pawrhoij(iatm)%nspden
             do irhoij=1,pawrhoij(iatm)%nrhoijsel
               klmn=pawrhoij(iatm)%rhoijselect(irhoij)
               pawrhoij(iatm)%rhoijres(klmn,ispden)=-pawrhoij(iatm)%rhoijp(irhoij,ispden)
             end do
           end do
         else
           do ispden=1,pawrhoij(iatm)%nspden
             do irhoij=1,pawrhoij(iatm)%nrhoijsel
               klmn1=2*pawrhoij(iatm)%rhoijselect(irhoij);jrhoij=2*irhoij
               pawrhoij(iatm)%rhoijres(klmn1-1,ispden)=-pawrhoij(iatm)%rhoijp(jrhoij-1,ispden)
               pawrhoij(iatm)%rhoijres(klmn1  ,ispden)=-pawrhoij(iatm)%rhoijp(jrhoij  ,ispden)
             end do
           end do
         end if
       end if
       nselect=0
       if (cplex==1) then
         do klmn=1,pawrhoij(iatm)%lmn2_size
           if (any(abs(pawrhoij_unsym(iatom)%rhoij_(klmn,:))>tol10)) then
             nselect=nselect+1
             pawrhoij(iatm)%rhoijselect(nselect)=klmn
             do jj=1,pawrhoij(iatm)%nspden
               ro(1)=pawrhoij_unsym(iatom)%rhoij_(klmn,jj)
               pawrhoij(iatm)%rhoijp(nselect,jj)=ro(1)
               if (use_res) pawrhoij(iatm)%rhoijres(klmn,jj)=pawrhoij(iatm)%rhoijres(klmn,jj)+ro(1)
             end do
           end if
         end do
       else
         do klmn=1,pawrhoij(iatm)%lmn2_size
           klmn1=2*klmn
           if (any(abs(pawrhoij_unsym(iatom)%rhoij_(klmn1-1:klmn1,:))>tol10)) then
             nselect=nselect+1;nselect1=2*nselect
             pawrhoij(iatm)%rhoijselect(nselect)=klmn
             do jj=1,pawrhoij(iatm)%nspden
               ro(1)=pawrhoij_unsym(iatom)%rhoij_(klmn1-1,jj)
               ro(2)=pawrhoij_unsym(iatom)%rhoij_(klmn1  ,jj)
               pawrhoij(iatm)%rhoijp(nselect1-1,jj)=ro(1)
               pawrhoij(iatm)%rhoijp(nselect1  ,jj)=ro(2)
               if (use_res) then
                 pawrhoij(iatm)%rhoijres(klmn1-1,jj)=pawrhoij(iatm)%rhoijres(klmn1-1,jj)+ro(1)
                 pawrhoij(iatm)%rhoijres(klmn1  ,jj)=pawrhoij(iatm)%rhoijres(klmn1  ,jj)+ro(2)
               end if
             end do
           end if
         end do
       end if
       pawrhoij(iatm)%nrhoijsel=nselect
     end do
   end if

 end if

!*********************************************************************
!Printing of Rhoij
 if (optrhoij==1.and.pawprtvol/=-10001) then
   wrt_mode='COLL';if (paral_atom) wrt_mode='PERS'
   pertstrg="RHOIJ";if (ipert>0) pertstrg="RHOIJ(1)"
   natinc=1;if(nrhoij>1.and.pawprtvol>=0) natinc=nrhoij-1
   do iatm=1,nrhoij,natinc
     iatom=iatm;if (paral_atom) iatom=my_atmtab(iatm)
     if (nrhoij==1.and.ipert>0.and.ipert<=natom) iatom=ipert
     idum1=2;if (pawrhoij(iatm)%cplex==2.and.pawrhoij(iatm)%nspinor==1) idum1=1
     if (abs(pawprtvol)>=1) then
       write(msg, '(6a,i3,a)') ch10," PAW TEST:",ch10,&
&       ' ====== Values of ',trim(pertstrg),' in symrhoij (iatom=',iatom,') ======'
       call wrtout(std_out,msg,wrt_mode)
     end if
     do ispden=1,pawrhoij(iatm)%nspden
       if (abs(pawprtvol)>=1.and.pawrhoij(iatm)%nspden/=1) then
         write(msg, '(3a)') '   Component ',trim(dspin(ispden+2*(pawrhoij(iatm)%nspden/4))),':'
       else if (pawrhoij(iatm)%nspden/=1) then
         if (pawrhoij(iatm)%nspden/=4) write(msg, '(4a,i3,a,i1,a)') ch10,&
&         ' *********** ',trim(pertstrg),' (atom ',iatom,', ispden=',ispden,') **********'
         if (pawrhoij(iatm)%nspden==4) write(msg, '(4a,i3,3a)') ch10,&
&         ' *********** ',trim(pertstrg),' (atom ',iatom,' - ',&
&         trim(dspin(ispden+2*(pawrhoij(iatm)%nspden/4))),') **********'
       else
         write(msg, '(4a,i3,a)') ch10,&
&         ' *********** ',trim(pertstrg),' (atom ',iatom,') **********'
       end if
       call wrtout(std_out,msg,wrt_mode)
       call pawio_print_ij(std_out,pawrhoij(iatm)%rhoijp(:,ispden),&
&       pawrhoij(iatm)%nrhoijsel,&
&       pawrhoij(iatm)%cplex,&
&       pawrhoij(iatm)%lmn_size,-1,idum,1,pawprtvol,&
&       pawrhoij(iatm)%rhoijselect(:),&
&       10.d0*dble(3-2*(ispden+ipert)),1,&
&       opt_sym=idum1,mode_paral=wrt_mode)
     end do
   end do
   msg=''
   call wrtout(std_out,msg,wrt_mode)
 end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

!*********************************************************************
!Small function: convert a symmetry operation
!from reduced coordinates (integers) to cartesian coordinates (reals)
 contains
   function symrhoij_symcart(aprim,bprim,symred)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symrhoij_symcart'
!End of the abilint section

   implicit none
   real(dp) :: symrhoij_symcart(3,3)
   integer,intent(in) :: symred(3,3)
   real(dp),intent(in) :: aprim(3,3),bprim(3,3)
   integer :: ii,jj,kk
   real(dp) :: tmp(3,3)
   symrhoij_symcart=zero;tmp=zero
   do kk=1,3
     do jj=1,3
       do ii=1,3
         tmp(ii,jj)=tmp(ii,jj)+bprim(ii,kk)*dble(symred(jj,kk))
       end do
     end do
   end do
   do kk=1,3
     do jj=1,3
       do ii=1,3
         symrhoij_symcart(ii,jj)=symrhoij_symcart(ii,jj)+aprim(ii,kk)*tmp(jj,kk)
       end do
     end do
   end do
   end function symrhoij_symcart

end subroutine symrhoij
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/pawrhoij_isendreceive_getbuffer
!! NAME
!! pawrhoij_isendreceive_getbuffer
!!
!! FUNCTION
!!  Fill a pawrhoij structure with the buffers received in a receive operation
!!  This buffer should have been first extracted by a call to pawrhoij_isendreceive_fillbuffer
!!
!! INPUTS
!!  atm_indx_recv(1:total number of atoms)= array for receive operation
!!                 Given an index of atom in global numbering, give its index
!!                 in the table of atoms treated by current processor
!!                 or -1 if the atoms is not treated by current processor
!!  buf_int= buffer of receive integers
!!  buf_dp= buffer of receive double precision numbers
!!  nrhoij_send= number of sent atoms
!!
!! OUTPUT
!!  pawrhoij= output datastructure filled with buffers receive in a receive operation
!!
!! PARENTS
!!      m_pawrhoij
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawrhoij_isendreceive_getbuffer(pawrhoij,nrhoij_send,atm_indx_recv,buf_int,buf_dp)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawrhoij_isendreceive_getbuffer'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nrhoij_send
!arrays
 integer,intent(in) ::atm_indx_recv(:),buf_int(:)
 real(dp),intent(in) :: buf_dp(:)
 type(pawrhoij_type),target,intent(inout) :: pawrhoij(:)

!Local variables-------------------------------
!scalars
 integer :: buf_dp_size,buf_int_size,cplex,ii,indx_int,indx_dp,iatom_tot,irhoij_send
 integer :: isp,jrhoij,lmn2_size,lmnmix,ngrhoij,nselect,nspden,rhoij_size2,use_rhoijp
 integer :: use_rhoijres,use_rhoijim,use_rhoij_
 character(len=500) :: msg
 type(pawrhoij_type),pointer :: pawrhoij1
!arrays

! *********************************************************************

 buf_int_size=size(buf_int)
 buf_dp_size=size(buf_dp)
 indx_int=1;indx_dp=1

 do irhoij_send=1,nrhoij_send

   iatom_tot=buf_int(indx_int) ; indx_int=indx_int+1
   jrhoij= atm_indx_recv(iatom_tot)
   if (jrhoij==-1)  then
     msg="Error in pawrhoij_isendreceive_getbuffer atom not found"
     MSG_BUG(msg)
   end if
   pawrhoij1=>pawrhoij(jrhoij)

   cplex       =buf_int(indx_int)    ;indx_int=indx_int+1
   lmn2_size   =buf_int(indx_int)    ;indx_int=indx_int+1
   nspden      =buf_int(indx_int)    ;indx_int=indx_int+1
   nselect     =buf_int(indx_int)    ;indx_int=indx_int+1
   lmnmix      =buf_int(indx_int)    ;indx_int=indx_int+1
   ngrhoij     =buf_int(indx_int)    ;indx_int=indx_int+1
   use_rhoijp  =buf_int(indx_int)    ;indx_int=indx_int+1
   use_rhoijres=buf_int(indx_int)    ;indx_int=indx_int+1
   use_rhoijim =buf_int(indx_int)    ;indx_int=indx_int+1
   use_rhoij_  =buf_int(indx_int)    ;indx_int=indx_int+1
   rhoij_size2 =buf_int(indx_int)    ;indx_int=indx_int+1
   pawrhoij1%itypat=buf_int(indx_int)   ;indx_int=indx_int+1
   pawrhoij1%lmn_size=buf_int(indx_int) ;indx_int=indx_int+1
   pawrhoij1%nsppol=buf_int(indx_int)   ;indx_int=indx_int+1
   pawrhoij1%nspinor=buf_int(indx_int)  ;indx_int=indx_int+1
   pawrhoij1%cplex=cplex
   pawrhoij1%lmn2_size=lmn2_size
   pawrhoij1%nspden=nspden
   pawrhoij1%nrhoijsel=nselect
   pawrhoij1%lmnmix_sz=lmnmix
   pawrhoij1%ngrhoij=ngrhoij
   pawrhoij1%use_rhoijp=use_rhoijp
   pawrhoij1%use_rhoijres=use_rhoijres
   pawrhoij1%use_rhoijim =use_rhoijim
   pawrhoij1%use_rhoij_=use_rhoij_
   if (use_rhoijp>0) then
     LIBPAW_ALLOCATE(pawrhoij1%rhoijselect,(lmn2_size))
     pawrhoij1%rhoijselect(1:nselect)=buf_int(indx_int:indx_int+nselect-1)
     if (nselect < lmn2_size )pawrhoij1%rhoijselect(nselect+1:lmn2_size)=zero
     indx_int=indx_int+nselect
     LIBPAW_ALLOCATE(pawrhoij1%rhoijp,(cplex*lmn2_size,nspden))
     do isp=1,nspden
       pawrhoij1%rhoijp(1:cplex*nselect,isp)=buf_dp(indx_dp:indx_dp+cplex*nselect-1)
       if (nselect < lmn2_size )pawrhoij1%rhoijp(cplex*nselect+1:cplex*lmn2_size,isp)=zero
       indx_dp=indx_dp+cplex*nselect
     end do
   end if
   if (lmnmix>0) then
     LIBPAW_ALLOCATE(pawrhoij1%kpawmix,(lmnmix))
     pawrhoij1%kpawmix(1:lmnmix)=buf_int(indx_int:indx_int+lmnmix-1)
     indx_int=indx_int+lmnmix
   end if
   if (ngrhoij>0) then
     LIBPAW_ALLOCATE(pawrhoij1%grhoij,(ngrhoij,cplex*lmn2_size,nspden))
     do isp=1,nspden
       do ii=1,cplex*lmn2_size
         pawrhoij1%grhoij(1:ngrhoij,ii,isp)=buf_dp(indx_dp:indx_dp+ngrhoij-1)
         indx_dp=indx_dp+ngrhoij
       end do
     end do
   end if
   if (use_rhoijres>0) then
     LIBPAW_ALLOCATE(pawrhoij1%rhoijres,(cplex*lmn2_size,nspden))
     do isp=1,nspden
       pawrhoij1%rhoijres(1:cplex*lmn2_size,isp)=buf_dp(indx_dp:indx_dp+cplex*lmn2_size-1)
       indx_dp=indx_dp+cplex*lmn2_size
     end do
   end if
   if (use_rhoijim>0) then
     LIBPAW_ALLOCATE(pawrhoij1%rhoijim,(lmn2_size,nspden))
     do isp=1,nspden
       pawrhoij1%rhoijim(1:lmn2_size,isp)=buf_dp(indx_dp:indx_dp+lmn2_size-1)
       indx_dp=indx_dp+lmn2_size
     end do
   end if
   if (use_rhoij_>0) then
     LIBPAW_ALLOCATE(pawrhoij1%rhoij_,(cplex*lmn2_size,rhoij_size2))
     do isp=1,rhoij_size2
       pawrhoij1%rhoij_(1:cplex*lmn2_size,isp)=buf_dp(indx_dp:indx_dp+cplex*lmn2_size-1)
       indx_dp=indx_dp+cplex*lmn2_size
     end do
   end if
 end do !irhoij_send
 if ((indx_int/=1+buf_int_size).or.(indx_dp/=1+buf_dp_size)) then
   write(msg,'(a,i10,a,i10)') 'Wrong buffer sizes: buf_int_size=',buf_int_size,' buf_dp_size=',buf_dp_size
   MSG_BUG(msg)
 end if

end subroutine pawrhoij_isendreceive_getbuffer
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/pawrhoij_isendreceive_fillbuffer
!! NAME
!!  pawrhoij_isendreceive_fillbuffer
!!
!! FUNCTION
!!  Extract from pawrhoij and from the global index of atoms
!!  the buffers to send in a sending operation
!!  This function has to be coupled with a call to pawrhoij_isendreceive_getbuffer
!!
!! INPUTS
!!  atm_indx_send(1:total number of atoms)= array for send operation,
!!                 Given an index of atom in global numbering, give its index
!!                 in the table of atoms treated by current processor
!!                 or -1 if the atoms is not treated by current processor
!!  nrhoij_send= number of sent atoms
!!  pawrhoij= data structure from which are extract buffer int and buffer dp
!!
!! OUTPUT
!!  buf_int= buffer of integers to be sent
!!  buf_int_size= size of buffer of integers
!!  buf_dp= buffer of double precision numbers to be sent
!!  buf_dp_size= size of buffer of double precision numbers
!!
!! PARENTS
!!      m_pawrhoij
!!
!! CHILDREN
!!
!! SOURCE
!!
subroutine pawrhoij_isendreceive_fillbuffer(pawrhoij,atmtab_send, atm_indx_send,nrhoij_send,&
&                                           buf_int,buf_int_size,buf_dp,buf_dp_size)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawrhoij_isendreceive_fillbuffer'
!End of the abilint section

implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) ::  buf_dp_size,buf_int_size
 integer,intent(in) :: nrhoij_send
!arrays
 integer,intent(in) :: atmtab_send(:),atm_indx_send(:)
 integer,intent(out),allocatable :: buf_int(:)
 real(dp),intent(out),allocatable :: buf_dp(:)
 type(pawrhoij_type),target,intent(in) :: pawrhoij(:)

!Local variables-------------------------------
!scalars
 integer :: cplex,ii,indx_int,indx_dp, iatom_tot,irhoij,irhoij_send,isp,lmn2_size,lmnmix
 integer :: ngrhoij,nselect,nspden,rhoij_size2,use_rhoijp,use_rhoijres,use_rhoijim,use_rhoij_
 character(len=500) :: msg
 type(pawrhoij_type),pointer :: pawrhoij1
!arrays

! *********************************************************************

!Compute sizes of buffers
 buf_int_size=0;buf_dp_size=0
 nselect=0;lmnmix=0;ngrhoij=0;rhoij_size2=0
 use_rhoijp=0;use_rhoijres=0;use_rhoijim=0;use_rhoij_=0
 do irhoij_send=1,nrhoij_send
   iatom_tot=atmtab_send(irhoij_send)
   irhoij=atm_indx_send(iatom_tot)
   if (irhoij == -1) then
     msg="Error in pawrhoij_isendreceive_fillbuffer atom not found"
     MSG_BUG(msg)
   end if
   pawrhoij1=>pawrhoij(irhoij)
   cplex    =pawrhoij1%cplex
   lmn2_size=pawrhoij1%lmn2_size
   nspden   =pawrhoij1%nspden
   lmnmix=pawrhoij1%lmnmix_sz
   ngrhoij=pawrhoij1%ngrhoij
   use_rhoijp=pawrhoij1%use_rhoijp
   use_rhoijres=pawrhoij1%use_rhoijres
   use_rhoijim=pawrhoij1%use_rhoijim
   use_rhoij_=pawrhoij1%use_rhoij_
   buf_int_size=buf_int_size+16
   if (use_rhoijp>0) then
     nselect=pawrhoij1%nrhoijsel
     buf_int_size=buf_int_size+nselect
     buf_dp_size=buf_dp_size + cplex*nselect*nspden
   end if
   if (lmnmix>0)       buf_int_size=buf_int_size+lmnmix
   if (ngrhoij>0)      buf_dp_size=buf_dp_size + cplex*lmn2_size*nspden*ngrhoij
   if (use_rhoijres>0) buf_dp_size=buf_dp_size + cplex*lmn2_size*nspden
   if (use_rhoijim>0)  buf_dp_size=buf_dp_size + lmn2_size*nspden
   if (use_rhoij_>0) then
     rhoij_size2=size(pawrhoij1%rhoij_,dim=2)
     buf_dp_size=buf_dp_size + cplex*lmn2_size*rhoij_size2
   end if
 end do

!Fill input buffers
 LIBPAW_ALLOCATE(buf_int,(buf_int_size))
 LIBPAW_ALLOCATE(buf_dp,(buf_dp_size))
 indx_int=1;indx_dp =1
 lmnmix=0;ngrhoij=0;nselect=0;rhoij_size2=0
 use_rhoijp=0;use_rhoijres=0;use_rhoijim=0;use_rhoij_=0
 do irhoij_send=1,nrhoij_send
   iatom_tot=atmtab_send(irhoij_send)
   irhoij=atm_indx_send(iatom_tot)
   pawrhoij1=>pawrhoij(irhoij)
   cplex    =pawrhoij1%cplex
   lmn2_size=pawrhoij1%lmn2_size
   nspden   =pawrhoij1%nspden
   lmnmix=pawrhoij1%lmnmix_sz
   ngrhoij=pawrhoij1%ngrhoij
   use_rhoijp=pawrhoij1%use_rhoijp
   nselect=pawrhoij1%nrhoijsel
   use_rhoijres=pawrhoij1%use_rhoijres
   use_rhoijim =pawrhoij1%use_rhoijim
   use_rhoij_  =pawrhoij1%use_rhoij_
   rhoij_size2 =size(pawrhoij1%rhoij_,dim=2)
   buf_int(indx_int)=atmtab_send(irhoij_send)        ;indx_int=indx_int+1
   buf_int(indx_int)=cplex                           ;indx_int=indx_int+1
   buf_int(indx_int)=lmn2_size                       ;indx_int=indx_int+1
   buf_int(indx_int)=nspden                          ;indx_int=indx_int+1
   buf_int(indx_int)=nselect                         ;indx_int=indx_int+1
   buf_int(indx_int)=lmnmix                          ;indx_int=indx_int+1
   buf_int(indx_int)=ngrhoij                         ;indx_int=indx_int+1
   buf_int(indx_int)=use_rhoijp                      ;indx_int=indx_int+1
   buf_int(indx_int)=use_rhoijres                    ;indx_int=indx_int+1
   buf_int(indx_int)=use_rhoijim                     ;indx_int=indx_int+1
   buf_int(indx_int)=use_rhoij_                      ;indx_int=indx_int+1
   buf_int(indx_int)=rhoij_size2                     ;indx_int=indx_int+1
   buf_int(indx_int)=pawrhoij1%itypat      ;indx_int=indx_int+1
   buf_int(indx_int)=pawrhoij1%lmn_size    ;indx_int=indx_int+1
   buf_int(indx_int)=pawrhoij1%nsppol      ;indx_int=indx_int+1
   buf_int(indx_int)=pawrhoij1%nspinor     ;indx_int=indx_int+1
   if (use_rhoijp>0) then
     buf_int(indx_int:indx_int+nselect-1)=pawrhoij1%rhoijselect(1:nselect)
     indx_int=indx_int+nselect
     do isp=1,nspden
       buf_dp(indx_dp:indx_dp+cplex*nselect-1)=pawrhoij1%rhoijp(1:cplex*nselect,isp)
       indx_dp=indx_dp+cplex*nselect
     end do
   end if
   if (lmnmix>0) then
     buf_int(indx_int:indx_int+lmnmix-1)=pawrhoij1%kpawmix(1:lmnmix)
     indx_int=indx_int+lmnmix
   end if
   if (ngrhoij>0) then
     do isp=1,nspden
       do ii=1,cplex*lmn2_size
         buf_dp(indx_dp:indx_dp+ngrhoij-1)=pawrhoij1%grhoij(1:ngrhoij,ii,isp)
         indx_dp=indx_dp+ngrhoij
       end do
     end do
   end if
   if (use_rhoijres>0) then
     do isp=1,nspden
       buf_dp(indx_dp:indx_dp+cplex*lmn2_size-1)=pawrhoij1%rhoijres(1:cplex*lmn2_size,isp)
       indx_dp=indx_dp+cplex*lmn2_size
     end do
   end if
   if (use_rhoijim>0) then
     do isp=1,nspden
       buf_dp(indx_dp:indx_dp+lmn2_size-1)=pawrhoij1%rhoijim(1:lmn2_size,isp)
       indx_dp=indx_dp+lmn2_size
     end do
   end if
   if (use_rhoij_>0) then
     do isp=1,rhoij_size2
       buf_dp(indx_dp:indx_dp+cplex*lmn2_size-1)=pawrhoij1%rhoij_(1:cplex*lmn2_size,isp)
       indx_dp=indx_dp+cplex*lmn2_size
     end do
   end if
 end do !irhoij_send

!Check
 if ((indx_int-1/=buf_int_size).or.(indx_dp-1/=buf_dp_size)) then
   write(msg,'(a,i10,a,i10)') 'Wrong buffer sizes: buf_int_size=',buf_int_size,' buf_dp_size=',buf_dp_size
   MSG_BUG(msg)
 end if

end subroutine pawrhoij_isendreceive_fillbuffer
!!***

!----------------------------------------------------------------------

END MODULE m_pawrhoij
!!***
