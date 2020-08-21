!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_pawfgrtab
!! NAME
!!  m_pawfgrtab
!!
!! FUNCTION
!!  This module contains the definition of the pawfgrtab_type structured datatype,
!!  as well as related functions and methods.
!!  pawfgrtab_type variables contain various data expressed on the fine grid
!!  for a given atom.
!!
!! COPYRIGHT
!! Copyright (C) 2013-2020 ABINIT group (MT, FJ)
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

MODULE m_pawfgrtab

 USE_DEFS
 USE_MSG_HANDLING
 USE_MPI_WRAPPERS
 USE_MEMORY_PROFILING

 use m_paral_atom, only : get_my_atmtab, free_my_atmtab, get_my_natom

 implicit none

 private
!!***

!----------------------------------------------------------------------


!!****t* m_pawfgrtab/pawfgrtab_type
!! NAME
!! pawfgrtab_type
!!
!! FUNCTION
!! For PAW, various arrays giving data related to fine grid for a given atom
!!
!! SOURCE

 type,public :: pawfgrtab_type

!Integer scalars

  integer :: cplex
   ! cplex=1 if potentials/densities are real, 2 if they are complex

  integer :: expiqr_allocated
   ! 1 if expiqr() is allocated (and computed)

  integer :: itypat
   ! itypat=type of the atom

  integer :: l_size
   ! 1+maximum value of l leading to non zero Gaunt coeffs
   ! for the considered atom type

  integer :: gylm_allocated
   ! 1 if gylm() is allocated (and computed)

  integer :: gylmgr_allocated
   ! 1 if gylmgr() is allocated (and computed)

  integer :: gylmgr2_allocated
   ! 1 if gylmgr2() is allocated (and computed)

  integer :: nfgd
   ! Number of Fine rectangular GriD points
   ! in the paw sphere around considered atom

  integer :: nhatfr_allocated
   ! 1 if nhatfr() is allocated (and computed)

  integer :: nhatfrgr_allocated
   ! 1 if nhatfrgr() is allocated (and computed)

  integer :: nspden
   ! Number of spin-density components

  integer :: rfgd_allocated
   ! 1 if rfgd() is allocated (and computed)

!Integer arrays
  !integer :: ngfftf(18)
  ! Info on the dense FFT mesh.

  integer, allocatable :: ifftsph(:)
   ! ifftsph(nfgd)
   ! Array giving the FFT index (fine grid) of a point in the paw
   ! sphere around considered atom (ifftsph=ix+n1*(iy-1+n2*(iz-1))

!Real (real(dp)) arrays

  real(dp), allocatable :: expiqr(:,:)
   ! expiqr(2,nfgd)
   ! Gives exp(-i.q.r) on the fine rectangular grid
   ! where q is the phonons wavevector

  real(dp), allocatable :: gylm(:,:)
   ! gylm(nfgd,l_size*l_size)
   ! Gives g_l(r)*Y_lm(r) on the fine rectangular grid
   ! around considered atom

  real(dp), allocatable :: gylmgr(:,:,:)
   ! gylmgr(3,nfgd,l_size*l_size)
   ! Gives the gradient of g_l(r)*Y_lm(r) wrt cart. coordinates
   ! on the fine rectangular grid around considered atom

  real(dp), allocatable :: gylmgr2(:,:,:)
   ! gylmgr(6,nfgd,l_size*l_size)
   ! Gives the second gradient of g_l(r)*Y_lm(r) wrt cart. coordinates
   ! on the fine rectangular grid around considered atom

  real(dp), allocatable :: nhatfr(:,:)
   ! nhatfr(nfgd,nspden)
   ! Gives the "frozen part" of 1st-order compensation charge density
   ! on the fine rectangular grid around considered atom
   ! Only used in response function calculations

  real(dp), allocatable :: nhatfrgr(:,:,:)
   ! nhatfrgr(3,nfgd,nspden)
   ! Gives the gradients (wrt cart. coordinates)
   ! of "frozen part" of 1st-order compensation charge density
   ! on the fine rectangular grid around considered atom
   ! Only used in response function calculations

  real(dp), allocatable :: rfgd(:,:)
   ! r(3,nfgd)
   ! Gives all R vectors (r-r_atom) on the Fine rectangular GriD
   ! around considered atom in Cartesian coordinates.

 end type pawfgrtab_type

!public procedures.
 public :: pawfgrtab_init           ! Constructor
 public :: pawfgrtab_free           ! Free memory
 public :: pawfgrtab_nullify
 public :: pawfgrtab_copy           ! Copy the object
 public :: pawfgrtab_print          ! Print the content of the object.
 public :: pawfgrtab_gather         ! MPI gather
 public :: pawfgrtab_redistribute   ! MPI redistribution

!private procedures.
 private :: pawfgrtab_isendreceive_getbuffer
 private :: pawfgrtab_isendreceive_fillbuffer
!!***

CONTAINS

!===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_pawfgrtab/pawfgrtab_init
!! NAME
!! pawfgrtab_init
!!
!! FUNCTION
!!  Initialize a pawfgrtab datastructure
!!
!! OUTPUT
!!  Pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!
!! PARENTS
!!      m_bethe_salpeter,m_classify_bands,m_d2frnl,m_exc_analyze,m_fock
!!      m_nonlinear,m_paw_mkaewf,m_paw_mkrho,m_respfn_driver,m_scfcv_core
!!      m_screening_driver,m_sigma_driver,m_wfd,m_wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawfgrtab_init(Pawfgrtab,cplex,l_size_atm,nspden,typat,&
&                         mpi_atmtab,comm_atom) ! optional arguments (parallelism)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nspden
!arrays
 integer,intent(in) :: l_size_atm(:),typat(:)
 integer,optional,intent(in) :: comm_atom
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 type(Pawfgrtab_type),intent(inout) :: Pawfgrtab(:)

!Local variables-------------------------------
!scalars
 integer :: iat,iatom_tot,my_comm_atom,my_natom,natom
 logical :: my_atmtab_allocated,paral_atom
 character(len=500) :: msg
!arrays
 integer,pointer :: my_atmtab(:)

! *************************************************************************

!@Pawfgrtab_type

 my_natom=SIZE(Pawfgrtab);if (my_natom==0) return
 natom=SIZE(typat)
 if (my_natom/=SIZE(l_size_atm)) then
  msg='Sizes of assumed shape arrays do not match'
  MSG_BUG(msg)
 end if

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

 call pawfgrtab_nullify(Pawfgrtab)

 do iat=1,my_natom
  iatom_tot=iat; if (paral_atom) iatom_tot=my_atmtab(iat)

  Pawfgrtab(iat)%cplex             = cplex
  Pawfgrtab(iat)%itypat            = typat(iatom_tot)
  Pawfgrtab(iat)%nspden            = nspden
  Pawfgrtab(iat)%l_size            = l_size_atm(iat)
  Pawfgrtab(iat)%nfgd              = 0
  LIBPAW_ALLOCATE(Pawfgrtab(iat)%ifftsph,(0))
  Pawfgrtab(iat)%gylm_allocated    = 0
  LIBPAW_ALLOCATE(Pawfgrtab(iat)%gylm,(0,0))
  Pawfgrtab(iat)%gylmgr_allocated  = 0
  LIBPAW_ALLOCATE(Pawfgrtab(iat)%gylmgr,(0,0,0))
  Pawfgrtab(iat)%gylmgr2_allocated = 0
  LIBPAW_ALLOCATE(Pawfgrtab(iat)%gylmgr2,(0,0,0))
  Pawfgrtab(iat)%nhatfr_allocated  = 0
  LIBPAW_ALLOCATE(Pawfgrtab(iat)%nhatfr,(0,0))
  Pawfgrtab(iat)%nhatfrgr_allocated  = 0
  LIBPAW_ALLOCATE(Pawfgrtab(iat)%nhatfrgr,(0,0,0))
  Pawfgrtab(iat)%rfgd_allocated    = 0
  LIBPAW_ALLOCATE(Pawfgrtab(iat)%rfgd,(0,0))
  Pawfgrtab(iat)%expiqr_allocated  = 0
  LIBPAW_ALLOCATE(Pawfgrtab(iat)%expiqr,(0,0))
 end do

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

end subroutine pawfgrtab_init
!!***

!----------------------------------------------------------------------

!!****f* m_pawfgrtab/pawfgrtab_free
!! NAME
!! pawfgrtab_free
!!
!! FUNCTION
!!  Free all dynamic memory stored in a pawfgrtab datastructure
!!
!! SIDE EFFECTS
!!  Pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!
!! PARENTS
!!      m_bethe_salpeter,m_classify_bands,m_d2frnl,m_exc_analyze,m_fock
!!      m_nonlinear,m_paral_pert,m_paw_dfpt,m_paw_mkaewf,m_paw_mkrho
!!      m_pawfgrtab,m_respfn_driver,m_scfcv_core,m_screening_driver
!!      m_sigma_driver,m_wfd,m_wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawfgrtab_free(Pawfgrtab)

!Arguments ------------------------------------
!arrays
 type(Pawfgrtab_type),intent(inout) :: Pawfgrtab(:)

!Local variables-------------------------------
!scalars
 integer :: iat,natom

! *************************************************************************

!@Pawfgrtab_type

 natom=SIZE(Pawfgrtab);if (natom==0) return
 do iat=1,natom
  if (allocated(Pawfgrtab(iat)%ifftsph))  then
    LIBPAW_DEALLOCATE(Pawfgrtab(iat)%ifftsph)
  end if
  if (allocated(Pawfgrtab(iat)%gylm   ))  then
    LIBPAW_DEALLOCATE(Pawfgrtab(iat)%gylm)
  end if
  if (allocated(Pawfgrtab(iat)%gylmgr ))  then
    LIBPAW_DEALLOCATE(Pawfgrtab(iat)%gylmgr)
  end if
  if (allocated(Pawfgrtab(iat)%gylmgr2))  then
    LIBPAW_DEALLOCATE(Pawfgrtab(iat)%gylmgr2)
  end if
  if (allocated(Pawfgrtab(iat)%nhatfr ))  then
    LIBPAW_DEALLOCATE(Pawfgrtab(iat)%nhatfr)
  end if
  if (allocated(Pawfgrtab(iat)%nhatfrgr))  then
    LIBPAW_DEALLOCATE(Pawfgrtab(iat)%nhatfrgr)
  end if
  if (allocated(Pawfgrtab(iat)%rfgd   ))  then
    LIBPAW_DEALLOCATE(Pawfgrtab(iat)%rfgd)
  end if
  if (allocated(Pawfgrtab(iat)%expiqr))   then
    LIBPAW_DEALLOCATE(Pawfgrtab(iat)%expiqr)
  end if
  Pawfgrtab(iat)%gylm_allocated=0
  Pawfgrtab(iat)%gylmgr_allocated=0
  Pawfgrtab(iat)%gylmgr2_allocated=0
  Pawfgrtab(iat)%nhatfr_allocated=0
  Pawfgrtab(iat)%nhatfrgr_allocated=0
  Pawfgrtab(iat)%rfgd_allocated=0
  Pawfgrtab(iat)%expiqr_allocated=0
 end do

end subroutine pawfgrtab_free
!!***

!----------------------------------------------------------------------

!!****f* m_pawfgrtab/pawfgrtab_nullify
!! NAME
!! pawfgrtab_nullify
!!
!! FUNCTION
!!  Nullify the pointers in a pawfgrtab datastructure
!!
!! SIDE EFFECTS
!!  Pawfgrtab(:) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!
!! PARENTS
!!      m_fock,m_paw_dfpt,m_pawfgrtab
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawfgrtab_nullify(Pawfgrtab)

!Arguments ------------------------------------
!arrays
 type(Pawfgrtab_type),intent(inout) :: Pawfgrtab(:)

!Local variables-------------------------------
!scalars
 integer :: iat,natom

! *************************************************************************

!@Pawfgrtab_type

! MGPAW: This one could be removed/renamed,
! variables can be initialized in the datatype declaration
! Do we need to expose this in the public API?

 natom=SIZE(Pawfgrtab);if (natom==0) return
 do iat=1,natom
  Pawfgrtab(iat)%nfgd              = 0
  Pawfgrtab(iat)%gylm_allocated    = 0
  Pawfgrtab(iat)%gylmgr_allocated  = 0
  Pawfgrtab(iat)%gylmgr2_allocated = 0
  Pawfgrtab(iat)%nhatfr_allocated  = 0
  Pawfgrtab(iat)%nhatfrgr_allocated= 0
  Pawfgrtab(iat)%rfgd_allocated    = 0
  Pawfgrtab(iat)%expiqr_allocated  = 0
 end do

end subroutine pawfgrtab_nullify
!!***

!----------------------------------------------------------------------

!!****f* m_pawfgrtab/pawfgrtab_copy
!! NAME
!!  pawfgrtab_copy
!!
!! FUNCTION
!!  Copy one pawfgrtab datastructure into another
!!  Can take into accound changes of dimensions
!!  Can copy a shared pawfgrtab into distributed ones (when parallelism is activated)
!!
!! INPUTS
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  pawfgrtab_in(:)<type(pawfgrtab_type)>= input paw_an datastructure
!!
!! SIDE EFFECTS
!!  pawfgrtab_copy(:)<type(pawfgrtab_type)>= output pawfgrtab datastructure
!!
!! NOTES
!!  paw_fgrtab_copy must have been allocated in the calling function.
!!
!! PARENTS
!!      m_pawfgrtab
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawfgrtab_copy(pawfgrtab_in,pawfgrtab_cp, &
&                         mpi_atmtab,comm_atom) ! optional arguments

!Arguments ------------------------------------
!scalars
integer,optional,intent(in) :: comm_atom
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 type(Pawfgrtab_type),intent(in) :: pawfgrtab_in(:)
 type(Pawfgrtab_type),intent(inout),target :: pawfgrtab_cp(:)

!Local variables-------------------------------
!scalars
 integer :: ij,ij1,istat,l_size,l_size2,my_comm_atom,my_natom,nfgd,nspden,paral_case
 integer :: siz_in, siz_max, siz_out
 logical ::  my_atmtab_allocated,paral_atom
 character(len=500) :: msg
!arrays
 integer,pointer :: my_atmtab(:)
 type(Pawfgrtab_type), pointer :: pawfgrtab_out(:)

! *************************************************************************

!@Pawfgrtab_type

!Retrieve sizes
 siz_in=size(pawfgrtab_in);siz_out=size(pawfgrtab_cp)

!Set up parallelism over atoms
 paral_atom=(present(comm_atom));if (paral_atom) paral_atom=(xmpi_comm_size(comm_atom)>1)
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 my_atmtab_allocated=.false.

!Determine in which case we are (parallelism, ...)
!No parallelism: a single copy operation
 paral_case=0;siz_max=siz_in
 pawfgrtab_out => pawfgrtab_cp
 if (paral_atom) then
   if (siz_out<siz_in) then ! Parallelism: the copy operation is a scatter
     call get_my_natom(my_comm_atom,my_natom,siz_in)
     if (my_natom==siz_out) then
       call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,siz_in)
       paral_case=1;siz_max=siz_out
       pawfgrtab_out => pawfgrtab_cp
     else
       msg=' siz_out should be equal to my_natom !'
       MSG_BUG(msg)
     end if
   else                            ! Parallelism: the copy operation is a gather
     call get_my_natom(my_comm_atom,my_natom,siz_out)
     if (my_natom==siz_in) then
       paral_case=2;siz_max=siz_in
     else
       msg=' siz_in should be equal to my_natom !'
       MSG_BUG(msg)
     end if
   end if
 end if

!First case: a simple copy or a scatter
 if (siz_max > 0 .and. ((paral_case== 0).or.(paral_case==1))) then
   call pawfgrtab_nullify(pawfgrtab_out)

!  Loop on pawfgrtab components
   do ij1=1,siz_max
     ij=ij1; if (paral_case==1)ij=my_atmtab(ij1)

     pawfgrtab_out(ij1)%cplex=pawfgrtab_in(ij)%cplex
     pawfgrtab_out(ij1)%expiqr_allocated=pawfgrtab_in(ij)%expiqr_allocated
     pawfgrtab_out(ij1)%itypat=pawfgrtab_in(ij)%itypat
     l_size=pawfgrtab_in(ij)%l_size
     l_size2=l_size*l_size
     pawfgrtab_out(ij1)%l_size=l_size
     pawfgrtab_out(ij1)%gylm_allocated=pawfgrtab_in(ij)%gylm_allocated
     pawfgrtab_out(ij1)%gylmgr_allocated=pawfgrtab_in(ij)%gylmgr_allocated
     pawfgrtab_out(ij1)%gylmgr2_allocated=pawfgrtab_in(ij)%gylmgr2_allocated
     pawfgrtab_out(ij1)%nfgd=pawfgrtab_in(ij)%nfgd
     pawfgrtab_out(ij1)%nhatfr_allocated=pawfgrtab_in(ij)%nhatfr_allocated
     pawfgrtab_out(ij1)%nhatfrgr_allocated=pawfgrtab_in(ij)%nhatfrgr_allocated
     nspden=pawfgrtab_in(ij)%nspden
     pawfgrtab_out(ij1)%nspden=nspden
     pawfgrtab_out(ij1)%rfgd_allocated=pawfgrtab_in(ij)%rfgd_allocated
     nfgd=pawfgrtab_in(ij)%nfgd
     LIBPAW_ALLOCATE(pawfgrtab_out(ij1)%ifftsph,(nfgd))
     pawfgrtab_out(ij1)%ifftsph(:)=pawfgrtab_in(ij)%ifftsph(:)
     if (pawfgrtab_out(ij1)%expiqr_allocated>0) then
       LIBPAW_ALLOCATE(pawfgrtab_out(ij1)%expiqr,(2,nfgd))
       pawfgrtab_out(ij1)%expiqr(:,:)=pawfgrtab_in(ij)%expiqr(:,:)
     end if
     if (pawfgrtab_out(ij1)%gylm_allocated>0) then
       LIBPAW_ALLOCATE(pawfgrtab_out(ij1)%gylm,(nfgd,l_size2))
       pawfgrtab_out(ij1)%gylm(:,:)=pawfgrtab_in(ij)%gylm(:,:)
     end if
     if (pawfgrtab_out(ij1)%gylmgr_allocated>0) then
       LIBPAW_ALLOCATE(pawfgrtab_out(ij1)%gylmgr,(3,nfgd,l_size2))
       pawfgrtab_out(ij1)%gylmgr(:,:,:)=pawfgrtab_in(ij)%gylmgr(:,:,:)
     end if
     if (pawfgrtab_out(ij1)%gylmgr2_allocated>0) then
       LIBPAW_ALLOCATE(pawfgrtab_out(ij1)%gylmgr2,(6,nfgd,l_size2))
       pawfgrtab_out(ij1)%gylmgr2(:,:,:)=pawfgrtab_in(ij)%gylmgr2(:,:,:)
     end if
     if (pawfgrtab_out(ij1)%nhatfrgr_allocated>0) then
       LIBPAW_ALLOCATE(pawfgrtab_out(ij1)%nhatfrgr,(3,nfgd,nspden))
       pawfgrtab_out(ij1)%nhatfrgr(:,:,:)=pawfgrtab_in(ij)%nhatfrgr(:,:,:)
      end if
      if (pawfgrtab_out(ij1)%rfgd_allocated>0) then
        LIBPAW_ALLOCATE(pawfgrtab_out(ij1)%rfgd,(3,nfgd))
        pawfgrtab_out(ij1)%rfgd(:,:)=pawfgrtab_in(ij)%rfgd(:,:)
      end if

    end do
 end if

!Second case: a gather
 if (paral_case==2) then
   call pawfgrtab_gather(pawfgrtab_in,pawfgrtab_cp,my_comm_atom,istat,my_atmtab)
 end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

end subroutine pawfgrtab_copy
!!***

!----------------------------------------------------------------------

!!****f* m_pawfgrtab/pawfgrtab_print
!! NAME
!! pawfgrtab_print
!!
!! FUNCTION
!!  Reports basic info on the pawfgrtab datatype.
!!
!! INPUTS
!!  Pawfgrtab<pawfgrtab_type>=The datatype to be printed
!!  [mode_paral]=either "COLL" or "PERS", "COLL" if None.
!!  [unit]=Unit number for output, std_out if None.
!!  [prtvol]=Verbosity level, lowest if None.
!!  [mpi_atmtab(:)]=indexes of the atoms treated by current proc
!!  [comm_atom]=MPI communicator over atoms
!!  [natom]=total number of atom (needed if parallelism over atoms is activated)
!!          if Pawfgrtab is distributed, natom is different from size(Pawfgrtab).
!!
!! OUTPUT
!! (only writing)
!!
!! PARENTS
!!      m_exc_analyze,m_paw_mkaewf,m_sigma_driver,m_wfd,m_wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawfgrtab_print(Pawfgrtab,natom,unit,prtvol,mode_paral,mpi_atmtab,comm_atom)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: comm_atom,natom,prtvol,unit
 character(len=4),optional,intent(in) :: mode_paral
!arrays
 type(Pawfgrtab_type),intent(inout) :: Pawfgrtab(:)
 integer,optional,target,intent(in)  :: mpi_atmtab(:)

!Local variables-------------------------------
!scalars
 integer :: iat,iatom_tot,my_comm_atom,my_natom,my_unt,my_prtvol,size_Pawfgrtab
 logical :: my_atmtab_allocated,paral_atom
 character(len=4) :: my_mode
 character(len=500) :: msg
!arrays
 integer,pointer :: my_atmtab(:)

! *************************************************************************

!@Pawfgrtab_type

 size_Pawfgrtab=SIZE(Pawfgrtab);if (size_Pawfgrtab==0) return

 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol
 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_mode  ='PERS' ; if (PRESENT(mode_paral)) my_mode  =mode_paral
 my_natom=size_Pawfgrtab; if (PRESENT(natom))  my_natom=natom

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.my_natom/=size_Pawfgrtab)
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,&
&                   my_natom,my_natom_ref=size_Pawfgrtab)

 write(msg,'(3a)')ch10,' === Content of the pawfgrtab datatype === ',ch10
 call wrtout(my_unt,msg,my_mode)
 do iat=1,my_natom
   iatom_tot=iat;if(paral_atom) iatom_tot=my_atmtab(iat)

   write(msg,'(6(2a,i0))')ch10,&
&    ' > For atom number : ',iatom_tot,ch10,&
&    '    cplex= ',Pawfgrtab(iat)%cplex,ch10,&
&    '    # of spins components for density= ',Pawfgrtab(iat)%nspden,ch10,&
&    '    Type of atom= ',Pawfgrtab(iat)%itypat,ch10,&
&    '    1+ Max l in Gaunt coefficients= ',Pawfgrtab(iat)%l_size,ch10,&
&    '    Number of fine FFT points in PAW sphere= ',Pawfgrtab(iat)%nfgd
   call wrtout(my_unt,msg,my_mode)

   if (my_prtvol>=3) then
     write(msg,'(a,7(a,i2,a))')ch10,&
&      '    rfgd_allocated    : ',Pawfgrtab(iat)%rfgd_allocated,ch10,&
&      '    gylm_allocated    : ',Pawfgrtab(iat)%gylm_allocated,ch10,&
&      '    gylmgr_allocated  : ',Pawfgrtab(iat)%gylmgr_allocated,ch10,&
&      '    gylmgr2_allocated : ',Pawfgrtab(iat)%gylmgr2_allocated,ch10,&
&      '    nhatgr_allocated  : ',Pawfgrtab(iat)%nhatfr_allocated,ch10,&
&      '    nhatfrgr_allocated: ',Pawfgrtab(iat)%nhatfrgr_allocated,ch10,&
&      '    expiqr_allocated  : ',Pawfgrtab(iat)%expiqr_allocated,ch10
     call wrtout(my_unt,msg,my_mode)

   end if

!  These huge arrays are not printed out!
!  Pawfgrtab(iat)%ifftsph
!  Pawfgrtab(iat)%rfgd
!  Pawfgrtab(iat)%gylm
!  Pawfgrtab(iat)%gylmgr
!  Pawfgrtab(iat)%gylmgr2
!  Pawfgrtab(ia)%nhatfr
!  Pawfgrtab(ia)%nhatfrgr
!  Pawfgrtab(ia)%expiqr
 end do

!Destroy atom table
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

end subroutine pawfgrtab_print
!!***

!----------------------------------------------------------------------

!!****f* m_pawfgrtab/pawfgrtab_gather
!! NAME
!! pawfgrtab_gather
!!
!! FUNCTION
!! gather a pawfgrtab
!! istat =1 if pawfgrtab is yet gathered
!!
!! INPUTS
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  mpicomm=MPI communicator
!!  pawfgrtab(:) =  pawfgrtab to be gathered
!!
!! OUTPUT
!!  istat : 1 if yet gathered in that case, pawfgrtab_gathered contains nothing
!!  pawfgrtab_gathered : pawfgrtab gathered between comm_atom
!!
!! PARENTS
!!      m_paw_dfpt,m_pawfgrtab
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawfgrtab_gather(pawfgrtab,pawfgrtab_gathered,comm_atom,istat, &
&                           mpi_atmtab) ! optional argument

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm_atom
 integer, intent(out) :: istat
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 type(pawfgrtab_type),intent(in) :: pawfgrtab(:)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab_gathered(:)

!Local variables-------------------------------
!scalars
 integer :: buf_int_size,buf_int_size_all,buf_dp_size,buf_dp_size_all,i1,i2,iat,iatot
 integer :: ierr,ij,indx_int,indx_dp,l_size,l_size2,my_natom,natom,nfgd,nproc_atom,nspden
 logical :: my_atmtab_allocated,paral_atom
 character(len=500) :: msg
!arrays
 integer :: bufsz(2)
 integer,allocatable :: buf_int(:),buf_int_all(:)
 integer,allocatable :: count_dp(:),count_int(:),count_tot(:),displ_dp(:),displ_int(:)
 integer,pointer :: my_atmtab(:)
 real(dp),allocatable :: buf_dp(:),buf_dp_all(:)

! *************************************************************************

!@Pawfgrtab_type

 istat=0
 my_natom=size(pawfgrtab);natom=size(pawfgrtab_gathered)

!Set up parallelism over atoms
 paral_atom=(my_natom/=natom)
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_atmtab_allocated = .False.
 if (paral_atom) then
   call get_my_atmtab(comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)
 end if
 nproc_atom=xmpi_comm_size(comm_atom)

!Without parallelism, just copy input to output
 if (.not.paral_atom) then
   istat=1;return
 end if

!Compute size of buffers
 buf_int_size=0
 do iat=1,my_natom
   buf_int_size=buf_int_size+pawfgrtab(iat)%nfgd+13
 end do
 buf_dp_size=0
 do iat=1,my_natom
   if (pawfgrtab(iat)%expiqr_allocated>0) then
     buf_dp_size=buf_dp_size + size(pawfgrtab(iat)%expiqr(:,:))
   end if
   if (pawfgrtab(iat)%gylm_allocated>0) then
     buf_dp_size=buf_dp_size+size(pawfgrtab(iat)%gylm(:,:))
   end if
   if (pawfgrtab(iat)%gylmgr_allocated>0) then
     buf_dp_size=buf_dp_size+size(pawfgrtab(iat)%gylmgr(:,:,:))
   end if
   if (pawfgrtab(iat)%gylmgr2_allocated>0) then
     buf_dp_size=buf_dp_size+size(pawfgrtab(iat)%gylmgr2(:,:,:))
   end if
   if (pawfgrtab(iat)%nhatfr_allocated>0) then
     buf_dp_size=buf_dp_size+size(pawfgrtab(iat)%nhatfr(:,:))
   end if
   if (pawfgrtab(iat)%nhatfrgr_allocated>0) then
     buf_dp_size=buf_dp_size+size(pawfgrtab(iat)%nhatfrgr(:,:,:))
   end if
   if (pawfgrtab(iat)%rfgd_allocated>0) then
     buf_dp_size=buf_dp_size+size(pawfgrtab(iat)%rfgd(:,:))
   end if
 end do

!Fill in input buffers
 LIBPAW_ALLOCATE(buf_int,(buf_int_size))
 LIBPAW_ALLOCATE(buf_dp,(buf_dp_size))
 indx_int=1;indx_dp=1
 do iat=1,my_natom
   l_size=pawfgrtab(iat)%l_size
   nfgd=pawfgrtab(iat)%nfgd
   nspden=pawfgrtab(iat)%nspden
   buf_int(indx_int)=my_atmtab(iat); indx_int=indx_int+1;
   buf_int(indx_int)=pawfgrtab(iat)%cplex; indx_int=indx_int+1;
   buf_int(indx_int)=pawfgrtab(iat)%itypat; indx_int=indx_int+1;
   buf_int(indx_int)=pawfgrtab(iat)%expiqr_allocated; indx_int=indx_int+1;
   buf_int(indx_int)=pawfgrtab(iat)%l_size; indx_int=indx_int+1;
   buf_int(indx_int)=pawfgrtab(iat)%gylm_allocated; indx_int=indx_int+1;
   buf_int(indx_int)=pawfgrtab(iat)%gylmgr_allocated; indx_int=indx_int+1;
   buf_int(indx_int)=pawfgrtab(iat)%gylmgr2_allocated; indx_int=indx_int+1;
   buf_int(indx_int)=pawfgrtab(iat)%nfgd; indx_int=indx_int+1;
   buf_int(indx_int)=pawfgrtab(iat)%nhatfr_allocated; indx_int=indx_int+1;
   buf_int(indx_int)=pawfgrtab(iat)%nhatfrgr_allocated; indx_int=indx_int+1;
   buf_int(indx_int)=pawfgrtab(iat)%nspden; indx_int=indx_int+1;
   buf_int(indx_int)=pawfgrtab(iat)%rfgd_allocated; indx_int=indx_int+1;
   if (nfgd>0) then
     buf_int(indx_int:indx_int+nfgd-1)=pawfgrtab(iat)%ifftsph(1:nfgd)
     indx_int=indx_int+nfgd;
   end if
   if (pawfgrtab(iat)%expiqr_allocated>0) then
     do i1=1,nfgd
       buf_dp(indx_dp:indx_dp+1)=pawfgrtab(iat)%expiqr(1:2,i1)
       indx_dp=indx_dp+2
     end do
   end if
   l_size2=l_size*l_size
   if (pawfgrtab(iat)%gylm_allocated>0) then
     do i1=1,l_size2
       buf_dp(indx_dp:indx_dp+nfgd-1)=pawfgrtab(iat)%gylm(1:nfgd,i1)
      indx_dp=indx_dp+nfgd
     end do
   end if
   if (pawfgrtab(iat)%gylmgr_allocated>0) then
     do i1=1,l_size2
       do i2=1,nfgd
         buf_dp(indx_dp:indx_dp+2)=pawfgrtab(iat)%gylmgr(1:3,i2,i1)
         indx_dp=indx_dp+3
       end do
     end do
   end if
   if (pawfgrtab(iat)%gylmgr2_allocated>0) then
     do i1=1,l_size2
       do i2=1,nfgd
         buf_dp(indx_dp:indx_dp+5)=pawfgrtab(iat)%gylmgr2(1:6,i2,i1)
         indx_dp=indx_dp+6
       end do
     end do
   end if
   if (pawfgrtab(iat)%nhatfr_allocated>0) then
     do i1=1,nspden
       buf_dp(indx_dp:indx_dp+nfgd-1)=pawfgrtab(iat)%nhatfr(1:nfgd,i1)
       indx_dp=indx_dp+nfgd
     end do
   end if
   if (pawfgrtab(iat)%nhatfrgr_allocated>0) then
     do i1=1,nspden
       do i2=1,nfgd
         buf_dp(indx_dp:indx_dp+2)=pawfgrtab(iat)%nhatfrgr(1:3,i2,i1)
         indx_dp=indx_dp+3
       end do
     end do
   end if
   if (pawfgrtab(iat)%rfgd_allocated>0) then
     do i1=1,nfgd
       buf_dp(indx_dp:indx_dp+2)=pawfgrtab(iat)%rfgd(1:3,i1)
       indx_dp=indx_dp+3
     end do
   end if
 end do
 if (indx_int/=1+buf_int_size) then
   msg='Error (1) in pawfgrtab_gather: wrong buffer sizes !'
   MSG_BUG(msg)
 end if
 if (indx_dp/=1+buf_dp_size) then
   msg='Error (2) in pawfgrtab_gather: wrong buffer sizes !'
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
 LIBPAW_ALLOCATE(buf_int_all,(buf_int_size_all))
 LIBPAW_ALLOCATE(buf_dp_all ,(buf_dp_size_all))
 call xmpi_allgatherv(buf_int,buf_int_size,buf_int_all,count_int,displ_int,comm_atom,ierr)
 call xmpi_allgatherv(buf_dp ,buf_dp_size ,buf_dp_all ,count_dp ,displ_dp ,comm_atom,ierr)
 LIBPAW_DEALLOCATE(count_int)
 LIBPAW_DEALLOCATE(displ_int)
 LIBPAW_DEALLOCATE(count_dp)
 LIBPAW_DEALLOCATE(displ_dp)

!Retrieve output datastructures
 indx_int=1; indx_dp=1
 call pawfgrtab_free(pawfgrtab_gathered)
 call pawfgrtab_nullify(pawfgrtab_gathered)
 do iat=1,natom
   iatot=buf_int_all(indx_int); indx_int=indx_int+1
   pawfgrtab_gathered(iatot)%cplex=buf_int_all(indx_int); indx_int=indx_int+1;
   pawfgrtab_gathered(iatot)%itypat=buf_int_all(indx_int); indx_int=indx_int+1;
   pawfgrtab_gathered(iatot)%expiqr_allocated=buf_int_all(indx_int); indx_int=indx_int+1;
   pawfgrtab_gathered(iatot)%l_size=buf_int_all(indx_int); indx_int=indx_int+1;
   pawfgrtab_gathered(iatot)%gylm_allocated=buf_int_all(indx_int); indx_int=indx_int+1;
   pawfgrtab_gathered(iatot)%gylmgr_allocated=buf_int_all(indx_int); indx_int=indx_int+1;
   pawfgrtab_gathered(iatot)%gylmgr2_allocated=buf_int_all(indx_int); indx_int=indx_int+1;
   pawfgrtab_gathered(iatot)%nfgd=buf_int_all(indx_int); indx_int=indx_int+1;
   pawfgrtab_gathered(iatot)%nhatfr_allocated=buf_int_all(indx_int); indx_int=indx_int+1;
   pawfgrtab_gathered(iatot)%nhatfrgr_allocated=buf_int_all(indx_int); indx_int=indx_int+1;
   pawfgrtab_gathered(iatot)%nspden=buf_int_all(indx_int); indx_int=indx_int+1;
   pawfgrtab_gathered(iatot)%rfgd_allocated=buf_int_all(indx_int); indx_int=indx_int+1;
   nfgd=pawfgrtab_gathered(iatot)%nfgd
   l_size=pawfgrtab_gathered(iatot)%l_size
   nspden= pawfgrtab_gathered(iatot)%nspden
   l_size2=l_size*l_size
   if (nfgd>0) then
     LIBPAW_ALLOCATE(pawfgrtab_gathered(iatot)%ifftsph,(nfgd))
     pawfgrtab_gathered(iatot)%ifftsph(1:nfgd)=buf_int_all(indx_int:indx_int+nfgd-1)
     indx_int=indx_int+nfgd
   end if
   if (pawfgrtab_gathered(iatot)%expiqr_allocated>0) then
     LIBPAW_ALLOCATE(pawfgrtab_gathered(iatot)%expiqr,(2,nfgd))
     do i1=1,nfgd
       pawfgrtab_gathered(iatot)%expiqr(1:2,i1)= buf_dp_all(indx_dp:indx_dp+1)
       indx_dp=indx_dp+2
     end do
   end if
   if (pawfgrtab_gathered(iatot)%gylm_allocated>0) then
     LIBPAW_ALLOCATE(pawfgrtab_gathered(iatot)%gylm,(nfgd,l_size2))
     do i1=1,l_size2
       pawfgrtab_gathered(iatot)%gylm(1:nfgd,i1)=buf_dp_all(indx_dp:indx_dp+nfgd-1)
       indx_dp=indx_dp+nfgd
     end do
   end if
   if (pawfgrtab_gathered(iatot)%gylmgr_allocated>0) then
     LIBPAW_ALLOCATE(pawfgrtab_gathered(iatot)%gylmgr,(3,nfgd,l_size2))
     do i1=1,l_size2
       do i2=1,nfgd
         pawfgrtab_gathered(iatot)%gylmgr(1:3,i2,i1)=buf_dp_all(indx_dp:indx_dp+2)
         indx_dp=indx_dp+3
       end do
     end do
   end if
   if (pawfgrtab_gathered(iatot)%gylmgr2_allocated>0) then
     LIBPAW_ALLOCATE(pawfgrtab_gathered(iatot)%gylmgr2,(6,nfgd,l_size2))
     do i1=1,l_size2
       do i2=1,nfgd
         pawfgrtab_gathered(iatot)%gylmgr2(1:6,i2,i1)=buf_dp_all(indx_dp:indx_dp+5)
         indx_dp=indx_dp+6
       end do
     end do
   end if
   if (pawfgrtab_gathered(iatot)%nhatfr_allocated>0) then
     LIBPAW_ALLOCATE(pawfgrtab_gathered(iatot)%nhatfr,(nfgd,nspden))
     do i1=1,nspden
       pawfgrtab_gathered(iatot)%nhatfr(1:nfgd,i1)=buf_dp_all(indx_dp:indx_dp+nfgd-1)
       indx_dp=indx_dp+nfgd
     end do
   end if
   if (pawfgrtab_gathered(iatot)%nhatfrgr_allocated>0) then
     LIBPAW_ALLOCATE(pawfgrtab_gathered(iatot)%nhatfrgr,(3,nfgd,nspden))
     do i1=1,nspden
       do i2=1,nfgd
         pawfgrtab_gathered(iatot)%nhatfrgr(1:3,i2,i1)=buf_dp_all(indx_dp:indx_dp+2)
         indx_dp=indx_dp+3
       end do
     end do
   end if
   if (pawfgrtab_gathered(iatot)%rfgd_allocated>0) then
     LIBPAW_ALLOCATE(pawfgrtab_gathered(iatot)%rfgd,(3,nfgd))
     do i1=1,nfgd
       pawfgrtab_gathered(iatot)%rfgd(1:3,i1)=buf_dp_all(indx_dp:indx_dp+2)
       indx_dp=indx_dp+3
     end do
   end if
 end do

!Free memory
 LIBPAW_DEALLOCATE(buf_int)
 LIBPAW_DEALLOCATE(buf_int_all)
 LIBPAW_DEALLOCATE(buf_dp)
 LIBPAW_DEALLOCATE(buf_dp_all)

!Destroy atom table
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

end subroutine pawfgrtab_gather
!!***

!----------------------------------------------------------------------

!!****f* m_pawfgrtab/pawfgrtab_redistribute
!! NAME
!! pawfgrtab_redistribute
!!
!! FUNCTION
!!   Redistribute an array of pawfgrtab datastructures
!!   Input pawfgrtab is given on a MPI communicator
!!   Output pawfgrtab is redistributed on another MPI communicator
!!
!! INPUTS
!!  mpi_comm_in= input MPI (atom) communicator
!!  mpi_comm_out= output MPI (atom) communicator
!!  mpi_atmtab_in= --optional-- indexes of the input pawfgrtab treated by current proc
!!                 if not present, will be calculated in the present routine
!!  mpi_atmtab_out= --optional-- indexes of the output pawfgrtab treated by current proc
!!                  if not present, will be calculated in the present routine
!!  natom= --optional-- total number of atoms
!!  ----- Optional arguments used only for asynchronous communications -----
!!    RecvAtomProc(:)= rank of processor from which I expect atom (in mpi_comm_in)
!!    RecvAtomList(:)= indexes of atoms to be received by me
!!      RecvAtomList(irecv) are the atoms I expect from RecvAtomProc(irecv)
!!    SendAtomProc(:)= ranks of process destination of atom (in mpi_comm_in)
!!    SendAtomList(:)= indexes of atoms to be sent by me
!!      SendAtomList(isend) are the atoms sent to SendAtomProc(isend)

!! OUTPUT
!!  [pawfgrtab_out(:)]<type(pawfgrtab_type)>= --optional--
!!                    if present, the redistributed datastructure does not replace
!!                    the input one but is delivered in pawfgrtab_out
!!                    if not present, input and output datastructure are the same.
!!
!! SIDE EFFECTS
!!  pawfgrtab(:)<type(pawfgrtab_type)>= input (and eventually output) pawfgrtab datastructures
!!
!! PARENTS
!!      m_paral_pert
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawfgrtab_redistribute(pawfgrtab,mpi_comm_in,mpi_comm_out,&
&                    natom,mpi_atmtab_in,mpi_atmtab_out,pawfgrtab_out,&
&                    SendAtomProc,SendAtomList,RecvAtomProc,RecvAtomList)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpi_comm_in,mpi_comm_out
 integer,optional,intent(in) :: natom
!arrays
 integer,intent(in),optional,target :: mpi_atmtab_in(:),mpi_atmtab_out(:)
 integer,intent(in),optional :: SendAtomProc(:),SendAtomList(:),RecvAtomProc(:),RecvAtomList(:)
 type(pawfgrtab_type),allocatable,intent(inout) :: pawfgrtab(:)
 type(pawfgrtab_type),pointer,optional :: pawfgrtab_out(:)

!Local variables-------------------------------
!scalars
 integer :: algo_option,iat_in,iat_out,iatom,i1,ierr,iircv,iisend,imsg,imsg1,imsg_current
 integer :: iproc_rcv,iproc_send,ireq,me_exch,mpi_comm_exch,my_natom_in,my_natom_out,my_tag,natom_tot,nb_dp
 integer :: nb_int,nb_msg,nbmsg_incoming,nbrecvmsg,nbsendreq,nbsend,nbsent,nbrecv,next,npawfgrtab_sent
 integer :: nproc_in,nproc_out
 logical :: flag,in_place,message_yet_prepared,my_atmtab_in_allocated,my_atmtab_out_allocated,paral_atom
!arrays
 integer :: buf_size(3),request1(3)
 integer,pointer :: my_atmtab_in(:),my_atmtab_out(:)
 integer,allocatable :: atmtab_send(:),atm_indx_in(:),atm_indx_out(:),buf_int1(:),From(:),request(:)
 integer,allocatable,target:: buf_int(:)
 integer,pointer :: buf_ints(:)
 logical,allocatable :: msg_pick(:)
 real(dp),allocatable :: buf_dp1(:)
 real(dp),allocatable,target :: buf_dp(:)
 real(dp),pointer :: buf_dps(:)
 type(coeffi1_type),target,allocatable :: tab_buf_int(:),tab_buf_atom(:)
 type(coeff1_type),target,allocatable :: tab_buf_dp(:)
 type(pawfgrtab_type),allocatable :: pawfgrtab_all(:)
 type(pawfgrtab_type),pointer :: pawfgrtab_out1(:)

! *************************************************************************

!@pawfgrtab_type

 in_place=(.not.present(pawfgrtab_out))
 my_natom_in=size(pawfgrtab)

!If not "in_place", destroy ) then
 if (.not.in_place) then
   if (associated(pawfgrtab_out)) then
     call pawfgrtab_free(pawfgrtab_out)
     LIBPAW_DATATYPE_DEALLOCATE(pawfgrtab_out)
   end if
 end if

!Special sequential case
 if (mpi_comm_in==xmpi_comm_self.and.mpi_comm_out==xmpi_comm_self) then
   if ((.not.in_place).and.(my_natom_in>0)) then
     LIBPAW_DATATYPE_ALLOCATE(pawfgrtab_out,(my_natom_in))
     call pawfgrtab_nullify(pawfgrtab_out)
     call pawfgrtab_copy(pawfgrtab,pawfgrtab_out)
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

   LIBPAW_DATATYPE_ALLOCATE(pawfgrtab_all,(natom_tot))
   call pawfgrtab_nullify(pawfgrtab_all)
   call pawfgrtab_copy(pawfgrtab,pawfgrtab_all,comm_atom=mpi_comm_in,mpi_atmtab=my_atmtab_in)
   if (in_place) then
    call pawfgrtab_free(pawfgrtab)
    LIBPAW_DATATYPE_DEALLOCATE(pawfgrtab)
    LIBPAW_DATATYPE_ALLOCATE(pawfgrtab,(my_natom_out))
    call pawfgrtab_nullify(pawfgrtab)
    call pawfgrtab_copy(pawfgrtab_all,pawfgrtab,comm_atom=mpi_comm_out,mpi_atmtab=my_atmtab_out)
   else
     LIBPAW_DATATYPE_ALLOCATE(pawfgrtab_out,(my_natom_out))
     call pawfgrtab_nullify(pawfgrtab_out)
     call pawfgrtab_copy(pawfgrtab_all,pawfgrtab_out,comm_atom=mpi_comm_out,mpi_atmtab=my_atmtab_out)
   end if
   call pawfgrtab_free(pawfgrtab_all)
   LIBPAW_DATATYPE_DEALLOCATE(pawfgrtab_all)


!Asynchronous algorithm (asynchronous communications)
!---------------------------------------------------------
 else if (algo_option==2) then

   nbsend=size(SendAtomProc) ; nbrecv=size(RecvAtomProc)

   if (in_place) then
     if (my_natom_out > 0) then
       LIBPAW_DATATYPE_ALLOCATE(pawfgrtab_out1,(my_natom_out))
       call pawfgrtab_nullify(pawfgrtab_out1)
     else
       LIBPAW_DATATYPE_ALLOCATE(pawfgrtab_out1,(0))
     end if
   else
     LIBPAW_DATATYPE_ALLOCATE(pawfgrtab_out,(my_natom_out))
     call pawfgrtab_nullify(pawfgrtab_out)
     pawfgrtab_out1=>pawfgrtab_out
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
           enddo
!          Create the message
           if (.not.message_yet_prepared) then
             nb_msg=nb_msg+1
             call pawfgrtab_isendreceive_fillbuffer( &
&                   pawfgrtab,atmtab_send,atm_indx_in,nbsent,buf_int,nb_int,buf_dp,nb_dp)
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
           my_tag=400
           ireq=ireq+1
           call xmpi_isend(buf_size,iproc_rcv,my_tag,mpi_comm_exch,request(ireq),ierr)
           my_tag=401
           ireq=ireq+1
           call xmpi_isend(buf_ints,iproc_rcv,my_tag,mpi_comm_exch,request(ireq),ierr)
           my_tag=402
           ireq=ireq+1
           call xmpi_isend(buf_dps,iproc_rcv,my_tag,mpi_comm_exch,request(ireq),ierr)
           nbsendreq=ireq
           nbsent=0
         end if
       end if
     else ! Just a renumbering, not a sending
       iat_in=atm_indx_in(SendAtomList(iisend))
       iat_out=atm_indx_out(my_atmtab_in(iat_in))
       call pawfgrtab_copy(pawfgrtab(iat_in:iat_in),pawfgrtab_out1(iat_out:iat_out))
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
         my_tag=400
         call xmpi_iprobe(iproc_send,my_tag,mpi_comm_exch,flag,ierr)
         if (flag) then
           msg_pick(i1)=.true.
           call xmpi_irecv(buf_size,iproc_send,my_tag,mpi_comm_exch,request1(1),ierr)
           call xmpi_wait(request1(1),ierr)
           nb_int=buf_size(1)
           nb_dp=buf_size(2)
           npawfgrtab_sent=buf_size(3)
           LIBPAW_ALLOCATE(buf_int1,(nb_int))
           LIBPAW_ALLOCATE(buf_dp1,(nb_dp))
           my_tag=401
           call xmpi_irecv(buf_int1,iproc_send,my_tag,mpi_comm_exch,request1(2),ierr)
           my_tag=402
           call xmpi_irecv(buf_dp1,iproc_send,my_tag,mpi_comm_exch,request1(3),ierr)
           call xmpi_waitall(request1(2:3),ierr)
           call pawfgrtab_isendreceive_getbuffer(pawfgrtab_out1,npawfgrtab_sent,atm_indx_out, buf_int1,buf_dp1)
           nbmsg_incoming=nbmsg_incoming-1
           LIBPAW_DEALLOCATE(buf_int1)
           LIBPAW_DEALLOCATE(buf_dp1)
         end if
       end if
     end do
   end do
   LIBPAW_DEALLOCATE(msg_pick)

   if (in_place) then
     call pawfgrtab_free(pawfgrtab)
     LIBPAW_DATATYPE_DEALLOCATE(pawfgrtab)
     LIBPAW_DATATYPE_ALLOCATE(pawfgrtab,(my_natom_out))
     call pawfgrtab_nullify(pawfgrtab)
     call pawfgrtab_copy(pawfgrtab_out1,pawfgrtab)
     call pawfgrtab_free(pawfgrtab_out1)
     LIBPAW_DATATYPE_DEALLOCATE(pawfgrtab_out1)
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

end subroutine pawfgrtab_redistribute
!!***

!----------------------------------------------------------------------

!!****f* m_pawfgrtab/pawfgrtab_isendreceive_getbuffer
!! NAME
!!  pawfgrtab_isendreceive_getbuffer
!!
!! FUNCTION
!!  Fill a pawfgrtab structure with the buffers received in a receive operation
!!  This buffer should have been first extracted by a call to pawfgrtab_isendreceive_fillbuffer
!!
!! INPUTS
!!  atm_indx_recv(1:total number of atoms)= array for receive operation
!!                 Given an index of atom in global numbering, give its index
!!                 in the table of atoms treated by current processor
!!                 or -1 if the atoms is not treated by current processor
!!  buf_int= buffer of receive integers
!!  buf_dp= buffer of receive double precision numbers
!!  npawfgrtab_send= number of sent atoms
!!
!! OUTPUT
!!  pawfgrtab= output datastructure filled with buffers receive in a receive operation
!!
!! PARENTS
!!      m_pawfgrtab
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawfgrtab_isendreceive_getbuffer(pawfgrtab,npawfgrtab_send,atm_indx_recv,buf_int,buf_dp)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npawfgrtab_send
!arrays
 integer,intent(in) ::atm_indx_recv(:),buf_int(:)
 real(dp),intent(in) :: buf_dp(:)
 type(pawfgrtab_type),target,intent(inout) :: pawfgrtab(:)

!Local variables-------------------------------
!scalars
 integer :: buf_dp_size,buf_int_size,i1,i2,iat,iatom_tot,ij,indx_int,indx_dp
 integer :: l_size,l_size2,nfgd,nspden
 character(len=500) :: msg
 type(pawfgrtab_type),pointer :: pawfgrtab1

! *********************************************************************

 buf_int_size=size(buf_int)
 buf_dp_size=size(buf_dp)
 indx_int=1; indx_dp=1

 do ij=1,npawfgrtab_send
   iatom_tot=buf_int(indx_int); indx_int=indx_int+1
   iat= atm_indx_recv(iatom_tot)
   pawfgrtab1=>pawfgrtab(iat)
   pawfgrtab1%cplex=buf_int(indx_int); indx_int=indx_int+1;
   pawfgrtab1%itypat=buf_int(indx_int); indx_int=indx_int+1;
   pawfgrtab1%expiqr_allocated=buf_int(indx_int); indx_int=indx_int+1;
   pawfgrtab1%l_size=buf_int(indx_int); indx_int=indx_int+1;
   pawfgrtab1%gylm_allocated=buf_int(indx_int); indx_int=indx_int+1;
   pawfgrtab1%gylmgr_allocated=buf_int(indx_int); indx_int=indx_int+1;
   pawfgrtab1%gylmgr2_allocated=buf_int(indx_int); indx_int=indx_int+1;
   pawfgrtab1%nfgd=buf_int(indx_int); indx_int=indx_int+1;
   pawfgrtab1%nhatfr_allocated=buf_int(indx_int); indx_int=indx_int+1;
   pawfgrtab1%nhatfrgr_allocated=buf_int(indx_int); indx_int=indx_int+1;
   pawfgrtab1%nspden=buf_int(indx_int); indx_int=indx_int+1;
   pawfgrtab1%rfgd_allocated=buf_int(indx_int); indx_int=indx_int+1;
   nfgd=pawfgrtab1%nfgd
   l_size=pawfgrtab1%l_size
   nspden= pawfgrtab1%nspden
   l_size2=l_size*l_size
   if (nfgd>0) then
     LIBPAW_ALLOCATE(pawfgrtab1%ifftsph,(nfgd))
     pawfgrtab1%ifftsph(1:nfgd)=buf_int(indx_int:indx_int+nfgd-1)
     indx_int=indx_int+nfgd
   end if
   if (pawfgrtab1%expiqr_allocated>0) then
     LIBPAW_ALLOCATE(pawfgrtab1%expiqr,(2,nfgd))
     do i1=1,nfgd
       pawfgrtab1%expiqr(1:2,i1)= buf_dp(indx_dp:indx_dp+1)
       indx_dp=indx_dp+2
     end do
   end if
   if (pawfgrtab1%gylm_allocated>0) then
     LIBPAW_ALLOCATE(pawfgrtab1%gylm,(nfgd,l_size2))
     do i1=1,l_size2
       pawfgrtab1%gylm(1:nfgd,i1)=buf_dp(indx_dp:indx_dp+nfgd-1)
       indx_dp=indx_dp+nfgd
     end do
   end if
   if (pawfgrtab1%gylmgr_allocated>0) then
     LIBPAW_ALLOCATE(pawfgrtab1%gylmgr,(3,nfgd,l_size2))
     do i1=1,l_size2
       do i2=1,nfgd
         pawfgrtab1%gylmgr(1:3,i2,i1)=buf_dp(indx_dp:indx_dp+2)
         indx_dp=indx_dp+3
       end do
     end do
   end if
   if (pawfgrtab1%gylmgr2_allocated>0) then
     LIBPAW_ALLOCATE(pawfgrtab1%gylmgr2,(6,nfgd,l_size2))
     do i1=1,l_size2
       do i2=1,nfgd
         pawfgrtab1%gylmgr2(1:6,i2,i1)=buf_dp(indx_dp:indx_dp+5)
         indx_dp=indx_dp+6
       end do
     end do
   end if
   if (pawfgrtab1%nhatfr_allocated>0) then
     LIBPAW_ALLOCATE(pawfgrtab1%nhatfr,(nfgd,nspden))
     do i1=1,nspden
       pawfgrtab1%nhatfr(1:nfgd,i1)=buf_dp(indx_dp:indx_dp+nfgd-1)
       indx_dp=indx_dp+nfgd
     end do
   end if
   if (pawfgrtab1%nhatfrgr_allocated>0) then
     LIBPAW_ALLOCATE(pawfgrtab1%nhatfrgr,(3,nfgd,nspden))
     do i1=1,nspden
       do i2=1,nfgd
         pawfgrtab1%nhatfrgr(1:3,i2,i1)=buf_dp(indx_dp:indx_dp+2)
         indx_dp=indx_dp+3
       end do
     end do
   end if
   if (pawfgrtab1%rfgd_allocated>0) then
     LIBPAW_ALLOCATE(pawfgrtab1%rfgd,(3,nfgd))
     do i1=1,nfgd
       pawfgrtab1%rfgd(1:3,i1)=buf_dp(indx_dp:indx_dp+2)
       indx_dp=indx_dp+3
     end do
   end if
 end do
 if ((indx_int/=1+buf_int_size).or.(indx_dp/=1+buf_dp_size)) then
   write(msg,'(a,i10,a,i10)') 'Wrong buffer sizes: buf_int_size=',buf_int_size,' buf_dp_size=',buf_dp_size
   MSG_BUG(msg)
 end if

end subroutine pawfgrtab_isendreceive_getbuffer
!!***

!----------------------------------------------------------------------
!!****f* m_pawfgrtab/pawfgrtab_isendreceive_fillbuffer
!! NAME
!!  pawfgrtab_isendreceive_fillbuffer
!!
!! FUNCTION
!!  Extract from pawfgrtab and from the global index of atoms
!!  the buffers to send in a sending operation
!!  This function has to be coupled with a call to pawfgrtab_isendreceive_getbuffer
!!
!! INPUTS
!!  atm_indx_send(1:total number of atoms)= array for send operation,
!!                 Given an index of atom in global numbering, give its index
!!                 in the table of atoms treated by current processor
!!                 or -1 if the atoms is not treated by current processor
!!  npawfgrtab_send= number of sent atoms
!!  pawfgrtab= data structure from which are extract buffer int and buffer dp
!!
!! OUTPUT
!!  buf_int= buffer of integers to be sent
!!  buf_int_size= size of buffer of integers
!!  buf_dp= buffer of double precision numbers to be sent
!!  buf_dp_size= size of buffer of double precision numbers
!!
!! PARENTS
!!      m_pawfgrtab
!!
!! CHILDREN
!!
!! SOURCE
!!
subroutine pawfgrtab_isendreceive_fillbuffer(pawfgrtab, atmtab_send,atm_indx_send,npawfgrtab_send,&
&                                            buf_int,buf_int_size,buf_dp,buf_dp_size)

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: buf_int_size,buf_dp_size
 integer,intent(in) :: npawfgrtab_send
!arrays
 integer,intent(in) :: atmtab_send(:),atm_indx_send(:)
 integer,intent(out),allocatable :: buf_int(:)
 real(dp),intent(out),allocatable :: buf_dp(:)
 type(pawfgrtab_type),target,intent(inout) :: pawfgrtab(:)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,iat,iat1,iatom_tot,indx_int,indx_dp,l_size,l_size2,nfgd,nspden
 character(len=500) :: msg
 type(pawfgrtab_type),pointer :: pawfgrtab1

! *********************************************************************

!Compute size of buffers
 buf_int_size=0;buf_dp_size=0
 do iat1=1,npawfgrtab_send
   iatom_tot=atmtab_send(iat1)
   iat=atm_indx_send(iatom_tot)
   pawfgrtab1=>pawfgrtab(iat)
   buf_int_size=buf_int_size+pawfgrtab1%nfgd+13
   if (pawfgrtab1%expiqr_allocated>0) then
     buf_dp_size=buf_dp_size + size(pawfgrtab1%expiqr(:,:))
   end if
   if (pawfgrtab1%gylm_allocated>0) then
     buf_dp_size=buf_dp_size+size(pawfgrtab1%gylm(:,:))
   end if
   if (pawfgrtab1%gylmgr_allocated>0) then
     buf_dp_size=buf_dp_size+size(pawfgrtab1%gylmgr(:,:,:))
   end if
   if (pawfgrtab1%gylmgr2_allocated>0 ) then
     buf_dp_size=buf_dp_size+size(pawfgrtab1%gylmgr2(:,:,:))
   end if
   if (pawfgrtab1%nhatfr_allocated>0) then
     buf_dp_size=buf_dp_size+size(pawfgrtab1%nhatfr(:,:))
   end if
   if (pawfgrtab1%nhatfrgr_allocated>0) then
     buf_dp_size=buf_dp_size+size(pawfgrtab1%nhatfrgr(:,:,:))
   end if
   if (pawfgrtab1%rfgd_allocated>0) then
     buf_dp_size=buf_dp_size+size(pawfgrtab1%rfgd(:,:))
   end if
 end do

!Fill in input buffers
 LIBPAW_ALLOCATE(buf_int,(buf_int_size))
 LIBPAW_ALLOCATE(buf_dp,(buf_dp_size))
 indx_int=1;indx_dp=1
 do iat1=1,npawfgrtab_send
    iatom_tot=atmtab_send(iat1)
   iat=atm_indx_send(iatom_tot)
   pawfgrtab1=>pawfgrtab(iat)
   l_size=pawfgrtab1%l_size
   nfgd=pawfgrtab1%nfgd
   nspden=pawfgrtab1%nspden
   buf_int(indx_int)=iatom_tot; indx_int=indx_int+1;
   buf_int(indx_int)=pawfgrtab1%cplex; indx_int=indx_int+1;
   buf_int(indx_int)=pawfgrtab1%itypat; indx_int=indx_int+1;
   buf_int(indx_int)=pawfgrtab1%expiqr_allocated; indx_int=indx_int+1;
   buf_int(indx_int)=pawfgrtab1%l_size; indx_int=indx_int+1;
   buf_int(indx_int)=pawfgrtab1%gylm_allocated; indx_int=indx_int+1;
   buf_int(indx_int)=pawfgrtab1%gylmgr_allocated; indx_int=indx_int+1;
   buf_int(indx_int)=pawfgrtab1%gylmgr2_allocated; indx_int=indx_int+1;
   buf_int(indx_int)=pawfgrtab1%nfgd; indx_int=indx_int+1;
   buf_int(indx_int)=pawfgrtab1%nhatfr_allocated; indx_int=indx_int+1;
   buf_int(indx_int)=pawfgrtab1%nhatfrgr_allocated; indx_int=indx_int+1;
   buf_int(indx_int)=pawfgrtab1%nspden; indx_int=indx_int+1;
   buf_int(indx_int)=pawfgrtab1%rfgd_allocated; indx_int=indx_int+1;
   if (nfgd>0) then
     buf_int(indx_int:indx_int+nfgd-1)=pawfgrtab1%ifftsph(1:nfgd)
     indx_int=indx_int+nfgd;
   end if
   if (pawfgrtab1%expiqr_allocated>0) then
     do i1=1,nfgd
       buf_dp(indx_dp:indx_dp+1)=pawfgrtab1%expiqr(1:2,i1)
       indx_dp=indx_dp+2
     end do
   end if
   l_size2=l_size*l_size
   if (pawfgrtab1%gylm_allocated>0) then
     do i1=1,l_size2
       buf_dp(indx_dp:indx_dp+nfgd-1)=pawfgrtab1%gylm(1:nfgd,i1)
      indx_dp=indx_dp+nfgd
     end do
   end if
   if (pawfgrtab1%gylmgr_allocated>0) then
     do i1=1,l_size2
       do i2=1,nfgd
         buf_dp(indx_dp:indx_dp+2)=pawfgrtab1%gylmgr(1:3,i2,i1)
         indx_dp=indx_dp+3
       end do
     end do
   end if
   if (pawfgrtab1%gylmgr2_allocated>0) then
     do i1=1,l_size2
       do i2=1,nfgd
         buf_dp(indx_dp:indx_dp+5)=pawfgrtab1%gylmgr2(1:6,i2,i1)
         indx_dp=indx_dp+6
       end do
     end do
   end if
   if (pawfgrtab1%nhatfr_allocated>0) then
     do i1=1,nspden
       buf_dp(indx_dp:indx_dp+nfgd-1)=pawfgrtab1%nhatfr(1:nfgd,i1)
       indx_dp=indx_dp+nfgd
     end do
   end if
   if (pawfgrtab1%nhatfrgr_allocated>0) then
     do i1=1,nspden
       do i2=1,nfgd
         buf_dp(indx_dp:indx_dp+2)=pawfgrtab1%nhatfrgr(1:3,i2,i1)
         indx_dp=indx_dp+3
       end do
     end do
   end if
   if (pawfgrtab1%rfgd_allocated>0) then
     do i1=1,nfgd
       buf_dp(indx_dp:indx_dp+2)=pawfgrtab1%rfgd(1:3,i1)
       indx_dp=indx_dp+3
     end do
   end if
 end do
 if ((indx_int-1/=buf_int_size).or.(indx_dp-1/=buf_dp_size)) then
   write(msg,'(a,i10,a,i10)') 'Wrong buffer sizes: buf_int_size=',buf_int_size,' buf_dp_size=',buf_dp_size
   MSG_BUG(msg)
 end if

end subroutine pawfgrtab_isendreceive_fillbuffer
!!***

!----------------------------------------------------------------------

END MODULE m_pawfgrtab
!!***
