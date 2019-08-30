!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_entropyDMFT
!! NAME
!!  m_entropyDMFT
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_entropyDMFT

  use defs_basis
  use m_errors
  use m_abicore
  use m_xmpi
  use m_dtset

  use m_energies, only : energies_type, energies_eval_eint
  use m_splines, only : spline_integrate, spline, splint
  use m_pawang, only : pawang_type
  use m_pawrad, only : pawrad_type, simp_gen, poisson
  use m_pawtab, only : pawtab_type
  use m_paw_correlations,only : pawpuxinit
  use m_io_tools, only : get_unit
  use m_data4entropyDMFT

  implicit none

  private

  public :: entropyDMFT_init
  public :: entropyDMFT_destroy
  public :: entropyDMFT_nextLambda
  public :: entropyDMFT_addIntegrand
  public :: entropyDMFT_computeEntropy

  integer, parameter :: E_U0       = 1
  integer, parameter :: E_UU       = 2
  integer, parameter :: E_DIRECT   = 1
  integer, parameter :: E_DC       = 2
  integer, parameter :: AC_NOTHING = 0
  integer, parameter :: AC_ETOT    = 1
  character(len=21), parameter :: HDR_NAME = "DATA FOR ETOT DMFT v="

!!***

!!****t* m_entropyDMFT/entropyDMFT
!! NAME
!!  entropyDMFT
!!
!! FUNCTION
!!  This structured datatype contains the necessary data
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

  type, public :: entropyDMFT_t
    logical               :: isset=.FALSE.  ! flag to be sure we are initialize
    integer               :: spacecomm      ! MPI comm
    integer               :: rank           ! rank in the comm
    integer               :: comm_size      ! Number of cpus in the comm
    integer               :: action         ! what to do in gstate
    integer               :: mylambda       ! Current lambda
    integer               :: natom          ! number of atoms
    integer               :: ncatom         ! number of correlated atoms
    integer               :: ntypat         ! number of type of atoms
    integer               :: nctypat        ! number of type of correlated atoms
    integer               :: nlambda        ! number of integration points
    integer               :: ofile          ! unit for file output data
    character(len=fnlen)  :: filename       ! name fo the file user readable
    character(len=fnlen)  :: filedata       ! name fo the file for restart purposes
    character(len=fnlen)  :: ofilename      ! ofstream prefix
    character(len=fnlen)  :: ifilename      ! ifstream prefix
    real(dp)              :: temp           ! temperature
    real(dp)              :: entropy0       ! entropy for lambda=0
    real(dp)              :: energies(2,2)  ! internal energy for lambda=0 and 1
    integer , allocatable :: index_atom(:)  ! index for correlated atoms
    integer , allocatable :: index_typat(:) ! index for correlated types
    integer , allocatable :: lpawu(:)       ! orbital moment to treat (ntypat,2)
    integer , allocatable :: typat(:)       ! type of each correlated atom
    real(dp), allocatable :: U_input(:)     ! U from input file (ntypat)
    real(dp), allocatable :: J_input(:)     ! J from input file (ntypat)
    real(dp), allocatable :: lambda(:)      ! ilamda
    real(dp), allocatable :: docc(:,:,:)    ! n_In_j, natom, ilamda
    real(dp), allocatable :: e_dc(:,:)      ! natom, ilamda
    real(dp), allocatable :: uij(:,:)       ! uij to compute <uij n_i n_j>, ntypat
  end type entropyDMFT_t
!!***


contains
!!***

!!****f* ABINIT/m_entropyDMFT/entropyDMFT_init
!! NAME
!!  entropyDMFT_init
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_scfcv
!!
!! CHILDREN
!!
!! SOURCE

subroutine entropyDMFT_init(e_t,dt,pawtab,spacecomm,ifilename,ofilename)

!Arguments ------------------------------------
    type(entropyDMFT_t) , intent(inout) :: e_t
    type(dataset_type)  , intent(in   ) :: dt
    type(pawtab_type)   , intent(in   ) :: pawtab(:)
    integer             , intent(in   ) :: spacecomm
    character(len=fnlen), intent(in   ) :: ifilename
    character(len=fnlen), intent(in   ) :: ofilename
!Local variables ------------------------------
    logical :: doRestart
    integer :: natom, ncatom
    integer :: iatom, icatom
    integer :: nctypat,itypat, ictypat
    integer :: ilambda
    integer, allocatable :: maptypat(:)
    character(len=500) :: message

    e_t%action = dt%dmft_entropy
    e_t%mylambda = e_t%action-1   ! should be usefull to start form a given value of lambda
                                  ! -1 added to start from 1 at the first call
                                  ! of nextLambda

    if ( e_t%action == AC_NOTHING ) then
      e_t%isset = .TRUE.
      return
    endif

    e_t%action = AC_ETOT

    natom = dt%natom
    ncatom = dt%natpawu

    icatom = 0
    do iatom=1,natom
      if ( pawtab(dt%typat(iatom))%lpawu /= -1 ) icatom = icatom + 1
    end do

    if ( icatom /= ncatom ) MSG_ERROR("Inconsistent number of correlated atoms")

    nctypat = 0
    do itypat=1,dt%ntypat
      if ( pawtab(itypat)%lpawu /= -1 ) nctypat = nctypat + 1
    end do

    e_t%natom  = natom
    e_t%ncatom = ncatom
    e_t%ntypat = dt%ntypat
    e_t%nctypat = nctypat

    if ( dt%dmft_nlambda < 3 ) then
      write(message,'(2a,i4,2a)') "DMFT must have dmft_nlamda >= 3 to compute entropy", &
        "whereas its value is dmft_nlambda = ",dt%dmft_nlambda,ch10,&
        "Action : check you input variable dmft_nlambda"
      MSG_ERROR(message)
    end if

    e_t%nlambda = dt%dmft_nlambda

    doRestart = .false.
    if ( e_t%nlambda < (e_t%mylambda+1) ) then
      MSG_ERROR("Restart calculation of DMFT entropy with a value of dmft_entropy greater than dmft_nlambda")
    else if ( e_t%mylambda > 0 ) then
      doRestart = .true.
    end if

    call entropyDMFT_allocateAll(e_t)

    e_t%entropy0      = zero
    e_t%energies(:,:) = zero

    e_t%lambda(:)     = (/ (DBLE(ilambda-1)/DBLE(e_t%nlambda-1), ilambda=1,e_t%nlambda) /)

    e_t%temp          = dt%tsmear

    ! Save each value of U and J for each correlated type
    ictypat = 1
    ABI_ALLOCATE(maptypat,(1:dt%ntypat))
    do itypat=1,dt%ntypat
      if ( pawtab(itypat)%lpawu /= -1 ) then
        e_t%lpawu(ictypat) = pawtab(itypat)%lpawu
        if(dt%dmft_t2g==1.and.e_t%lpawu(ictypat)==2) then
          e_t%lpawu(ictypat)=1
        end if
        e_t%U_input(ictypat) = pawtab(itypat)%upawu
        e_t%J_input(ictypat) = pawtab(itypat)%jpawu
        e_t%index_typat(ictypat) = itypat
        maptypat(itypat) = ictypat
        ictypat = ictypat + 1
      end if
    end do

    ! Save type local and global of each correlated atom
    ! Get the correct value of lpawu for each type
    icatom = 1
    do iatom=1,e_t%ncatom
      if ( pawtab(dt%typat(iatom))%lpawu /= -1 ) then
        e_t%typat(icatom) = maptypat(dt%typat(iatom))
        e_t%index_atom(icatom) = iatom
        icatom = icatom + 1
      end if
    end do
    ABI_DEALLOCATE(maptypat)

    write(message,'(a,1x,78a)') ch10,"+",(/ ("-",ilambda=1,76) /), "+"
    call wrtout(std_out,message,"COLL")
    call wrtout(ab_out,message,"COLL")
    write(message,'(1x,a)') "|             Calculation of entropy within  the DMFT Framework              |"
    call wrtout(std_out,message,"COLL")
    call wrtout(ab_out,message,"COLL")
    write(message,'(1x,40a)') "+",(/ ("- ",ilambda=1,38) /), "+"
    call wrtout(std_out,message,"COLL")
    do ilambda = 1, e_t%nlambda
      write(message,'(1x,a,i4,a11,f6.4,55x,a)') &
        "|", ilambda, ") lambda = ", e_t%lambda(ilambda), "|"
      call wrtout(std_out,message,"COLL")
      do ictypat=1,e_t%nctypat
        write(message,'(1x,a,6x,a12,i4,4x,a3,5x,2(3x,a4,f6.4,1x,a2),10x,a)') "|", &
          "- Atom type ", e_t%index_typat(ictypat), "->", &
          "U = ",e_t%U_input(ictypat)*e_t%lambda(ilambda), "Ha", &
          "J = ",e_t%J_input(ictypat)*e_t%lambda(ilambda), "Ha","|"
        call wrtout(std_out,message,"COLL")
      end do
    end do
    write(message,'(1x,78a)') "+",(/ ("-",ilambda=1,76) /), "+"
    call wrtout(std_out,message,"COLL")
    call wrtout(ab_out,message,"COLL")

    ! Set up MPI
    e_t%spacecomm = spacecomm
    e_t%rank      = xmpi_comm_rank(spacecomm)
    e_t%comm_size = xmpi_comm_size(spacecomm)

    e_t%ofile = get_unit()
    e_t%ofilename = ofilename
    e_t%ifilename = ifilename
    e_t%filename = TRIM(e_t%ofilename)//"_EntropyDMFT"
    e_t%filedata = TRIM(e_t%ofilename)//"_data4EntropyDMFT"

    if ( doRestart .eqv. .true. ) then
      ! If restart fails, then nothing changes and the full calculation is
      ! perform. Otherwise, we complete as much a possible the structure.
      call entropyDMFT_restart(e_t)
    end if

    ! Rewrite the files with the previous data
    call entropyDMFT_dump(e_t)

    e_t%isset = .TRUE.


  end subroutine entropyDMFT_init
!!***

!!****f* ABINIT/m_entropyDMFT/entropyDMFT_allocateAll
!! NAME
!!  entropyDMFT_allocateAll
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_entropyDMFT
!!
!! CHILDREN
!!
!! SOURCE

  subroutine entropyDMFT_allocateAll(e_t)

!Arguments ------------------------------------
    type(entropyDMFT_t), intent(inout) :: e_t

    ABI_ALLOCATE(e_t%index_atom, (1:e_t%natom))
    e_t%index_atom = 0
    ABI_ALLOCATE(e_t%index_typat,(1:e_t%nctypat))
    e_t%index_atom = 0
    ABI_ALLOCATE(e_t%typat,      (1:e_t%ncatom))
    e_t%typat      = 0
    ABI_ALLOCATE(e_t%lpawu,      (1:e_t%nctypat))
    e_t%lpawu      = 0
    ABI_ALLOCATE(e_t%U_input,    (1:e_t%nctypat))
    e_t%U_input    = zero
    ABI_ALLOCATE(e_t%J_input,    (1:e_t%nctypat))
    e_t%J_input    = zero
    ABI_ALLOCATE(e_t%lambda,     (1:e_t%nlambda))
    e_t%lambda     = 0
    ABI_ALLOCATE(e_t%docc,       (1:(14*13)/2,1:e_t%ncatom,1:e_t%nlambda))
    e_t%docc       = zero
    ABI_ALLOCATE(e_t%e_dc,       (1:e_t%ncatom,1:e_t%nlambda))
    e_t%e_dc       = zero
    ABI_ALLOCATE(e_t%uij,        (1:(14*13)/2,1:e_t%nctypat))
    e_t%uij        = zero
  end subroutine entropyDMFT_allocateAll
!!***

!!****f* ABINIT/m_entropyDMFT/entropyDMFT_restart
!! NAME
!!  entropyDMFT_restart
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_entropyDMFT
!!
!! CHILDREN
!!
!! SOURCE

  subroutine entropyDMFT_restart(e_t)

!Arguments ------------------------------------
    type(entropyDMFT_t), intent(inout) :: e_t
!Local variables-------------------------------
    logical            :: doBcast
    character(len=200) :: msg
    character(len=21 ) :: hdr
    integer            :: iostate
    integer            :: ilambda
    integer            :: tlambda
    integer            :: ictypat
    integer            :: ndim
    integer            :: icouple
    integer            :: iflavor1
    integer            :: iflavor2
    integer            :: icatom
    integer            :: iatom
    real(dp)           :: lambda

    write(msg,'(a,1x,2a,1x,a,f6.4,a)') ch10,"EtotDMFT will try to restart calculation from a previous run", ch10, &
    "and this calculation should start at lambda = ",e_t%lambda(e_t%mylambda+1), ch10
    call wrtout(std_out,msg,"COLL")
    call wrtout(ab_out,msg,"COLL")

    doBcast = .true.

    if ( e_t%rank == 0 ) then
      tlambda = 1
      inquire(file=e_t%filedata,iostat=iostate)

      if ( iostate /= 0 ) then
        write(msg,'(5a)') "File ", trim(e_t%filedata), " does not exist or is not accessible", ch10, &
        "-> No restart performed but full calculation."
        MSG_WARNING(msg)
        e_t%mylambda = 0
      end if

      open(unit=e_t%ofile,file=e_t%filedata,action="read",form="unformatted")
      read(e_t%ofile,end=42) hdr, iostate

      if ( hdr /= HDR_NAME .or. iostate /= 1 ) then
        write(msg,'(5a)') "File ", trim(e_t%filedata), " does not contain the proper header", ch10, &
        "-> No restart performed but full calculation."
        MSG_WARNING(msg)
        e_t%mylambda = 0
      end if

      do ilambda = 1, e_t%mylambda
        tlambda = ilambda
        read(e_t%ofile,end=42) lambda
        if ( ABS(lambda - e_t%lambda(ilambda)) >= tol9 ) then
          write(msg,'(5a,f6.4,a,f6.4)') "File ", trim(e_t%filedata), " is wrong:", ch10, &
          "Lambda values are differente: in file ", lambda, " instead of ", e_t%lambda(ilambda)
          MSG_WARNING(msg)
          goto 42
        end if

        if ( ilambda == 1 ) then
          read(e_t%ofile,end=42) lambda
          if ( lambda /= e_t%temp ) then
            write(msg,'(7a,f6.4)') "File ", trim(e_t%filedata), " is wrong:", ch10, &
            "Temperature is different than the value of tsmear", ch10, &
            "-> No restart performed but full calculation."
            MSG_WARNING(msg)
            goto 42
          end if
          read(e_t%ofile,end=42) e_t%entropy0
          read(e_t%ofile,end=42) e_t%energies(E_DC,E_U0)
        else if ( ilambda == e_t%nlambda ) then ! should never happend ?!
          read(e_t%ofile,end=42) e_t%energies(E_DC,E_UU)
          do ictypat = 1, e_t%nctypat
            ndim = 2*(2*e_t%lpawu(ictypat)+1)
            icouple = 0
            do iflavor1 = 1, ndim
              do iflavor2 = iflavor1+1, ndim
                icouple = icouple + 1
                read(e_t%ofile,end=42) e_t%uij(icouple,ictypat)
              end do
            end do
          end do
        end if

        do icatom = 1, e_t%ncatom
          iatom = e_t%index_atom(icatom)
          ndim = 2*(2*e_t%lpawu(e_t%typat(icatom))+1)
          icouple = 0
          read(e_t%ofile,end=42) e_t%e_dc(icatom,ilambda)
          do iflavor1 = 1, ndim
            do iflavor2 = iflavor1+1, ndim
              icouple = icouple + 1
              read(e_t%ofile,end=42) e_t%docc(icouple,icatom,ilambda)
            end do
          end do
        end do
      end do
      close(e_t%ofile)
      goto 43
42    write(msg,'(5a,f6.4)') "File ", trim(e_t%filedata), " is wrong or incomplete", ch10, &
          "-> Restart calculation will restart at lambda = ",e_t%lambda(tlambda)
      MSG_WARNING(msg)
      close(e_t%ofile)
      e_t%mylambda = tlambda-1 ! -1 to go to previous lambda
    end if
    ! MPI BDCAST
43  call xmpi_bcast(e_t%mylambda,0, e_t%spacecomm, ictypat)
    call xmpi_bcast(e_t%entropy0,0, e_t%spacecomm, ictypat)
    call xmpi_bcast(e_t%energies,0, e_t%spacecomm, ictypat)
    !call xmpi_bcast(e_t%uij,0, e_t%spacecomm, ictypat) ! No need since it
    ! restart always perform the last lambda and this is calculation at the end
    ! of last lambda
    call xmpi_bcast(e_t%e_dc,0, e_t%spacecomm, ictypat)
    call xmpi_bcast(e_t%docc,0, e_t%spacecomm, ictypat)

  end subroutine entropyDMFT_restart
!!***

!!****f* ABINIT/m_entropyDMFT/entropyDMFT_dump
!! NAME
!!  entropyDMFT_dump
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_entropyDMFT
!!
!! CHILDREN
!!
!! SOURCE

  subroutine entropyDMFT_dump(e_t)

!Arguments ------------------------------------
    type(entropyDMFT_t), intent(inout) :: e_t
!Local variables-------------------------------
    integer            :: ilambda
    integer            :: ictypat
    integer            :: ndim
    integer            :: icouple
    integer            :: iflavor1
    integer            :: iflavor2
    integer            :: icatom
    integer            :: iatom

    if ( e_t%rank /= 0 ) return

    ! Dump _data4EtotDMFT
    open(unit=e_t%ofile,file=e_t%filedata,action="write",form="unformatted")
    write(e_t%ofile) HDR_NAME, 1

    do ilambda = 1, e_t%mylambda
      write(e_t%ofile) e_t%lambda(ilambda)

      if ( ilambda == 1 ) then
        write(e_t%ofile) e_t%temp
        write(e_t%ofile) e_t%entropy0
        write(e_t%ofile) e_t%energies(E_DC,E_U0)
      else if ( ilambda == e_t%nlambda ) then ! should never happend ?!
        write(e_t%ofile) e_t%energies(E_DC,E_UU)
        do ictypat = 1, e_t%nctypat
          ndim = 2*(2*e_t%lpawu(ictypat)+1)
          icouple = 0
          do iflavor1 = 1, ndim
            do iflavor2 = iflavor1+1, ndim
              icouple = icouple + 1
              write(e_t%ofile) e_t%uij(icouple,ictypat)
            end do
          end do
        end do
      end if

      do icatom = 1, e_t%ncatom
        iatom = e_t%index_atom(icatom)
        ndim = 2*(2*e_t%lpawu(e_t%typat(icatom))+1)
        icouple = 0
        write(e_t%ofile) e_t%e_dc(icatom,ilambda)
        do iflavor1 = 1, ndim
          do iflavor2 = iflavor1+1, ndim
            icouple = icouple + 1
            write(e_t%ofile) e_t%docc(icouple,icatom,ilambda)
          end do
        end do
      end do
    end do
    close(e_t%ofile)

    ! Dump _EtotDMFT
    open(unit=e_t%ofile,file=e_t%filename)
    write(e_t%ofile,'(2a)') "# Data for entropy calculation in DMFT",ch10
    do ilambda = 1, e_t%mylambda
      if ( ilambda == 1 ) then
        write(e_t%ofile,'(a)') "# Temperature [Ha]:"
        write(e_t%ofile,'(es22.14)') e_t%temp
        write(e_t%ofile,'(a)') "# Entropy for lambda=0 [kb]:"
        write(e_t%ofile,'(es22.14)') e_t%entropy0
        write(e_t%ofile,'(a)') "# Internal energy for lambda=0 [Ha]:"
        write(e_t%ofile,'(es22.14)') e_t%energies(E_DC,E_U0)
      else if ( e_t%mylambda == e_t%nlambda ) then
        write(e_t%ofile,'(a)') "# Internal energy for lambda=1 [Ha]:"
        write(e_t%ofile,'(es22.14)') e_t%energies(E_DC,E_UU)
        do ictypat = 1, e_t%nctypat
          ndim = 2*(2*e_t%lpawu(ictypat)+1)
          icouple = 0
          write(e_t%ofile,'(a,f7.5,1x,a,i4)') "# Interaction Matrix normalized by U=",e_t%U_input(ictypat) , &
          "[Ha] for atom type", e_t%index_typat(ictypat)
          do iflavor1 = 1, ndim
            write(e_t%ofile,'(14(a21,2x))',advance="no") (/ ( "...", iflavor2=1,iflavor1 ) /)
            do iflavor2 = iflavor1+1, ndim
              icouple = icouple + 1
              write(e_t%ofile,'(14(es21.14,2x))',advance="no") e_t%uij(icouple,ictypat)
            end do
            write(e_t%ofile,*)
          end do
        end do
      end if
    end do
    close(e_t%ofile)
  end subroutine entropyDMFT_dump
!!***

!!****f* ABINIT/m_entropyDMFT/entropyDMFT_nextLambda
!! NAME
!!  entropyDMFT_nextLambda
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

  function entropyDMFT_nextLambda(e_t,dt,pawtab,pawang,pawrad) result(nextstep)

!Arguments ------------------------------------
    type(entropyDMFT_t) , intent(inout) :: e_t
    type(dataset_type)  , intent(in) :: dt
    type(pawtab_type)   , intent(inout) :: pawtab(:)
    type(pawang_type)   , intent(in   ) :: pawang
    type(pawrad_type)   , intent(inout) :: pawrad(:)
!Local variables ------------------------------
    logical :: nextstep
    integer :: itypat
    integer :: mylambda
    logical :: is_dfpt=.false.
    real(dp),allocatable :: upawu(:),jpawu(:)
    character(len=100) :: message

    if ( e_t%isset .eqv. .FALSE. ) &
      MSG_ERROR("entropyDMFT is not initialized")

    ! go to next lambda
    mylambda = e_t%mylambda + 1
    !if ( present(ilambda) ) mylambda = ilambda
    e_t%mylambda = mylambda

    if ( e_t%action == AC_NOTHING .and. mylambda >= 1 ) then
      nextstep = .FALSE.
      ! We do nothing and return false, one scfvc has already performed
    else if ( e_t%action == AC_NOTHING .and. mylambda < 1 ) then
      nextstep = .TRUE.
      ! We do nothing but return true to perform scfcv a least once
    else if ( e_t%action == AC_ETOT .and. mylambda <= e_t%nlambda ) then ! we iterate over lambda
      nextstep = .TRUE.
      ABI_ALLOCATE(upawu,(dt%ntypat))
      ABI_ALLOCATE(jpawu,(dt%ntypat))
      do itypat = 1, e_t%nctypat
        upawu(e_t%index_typat(itypat)) = e_t%lambda(mylambda) * e_t%U_input(itypat)
        jpawu(e_t%index_typat(itypat)) = e_t%lambda(mylambda) * e_t%J_input(itypat)
      end do
    else ! we did all lambda values
      nextstep = .FALSE.
    endif

    if ( e_t%action == AC_ETOT .and. nextstep .eqv. .true. ) then
      write(message,'(a,1x,78a)') ch10,"+",(/ ("-",itypat=1,76) /), "+"
      call wrtout(std_out,message,"COLL")
      call wrtout(ab_out,message,"COLL")
      write(message,'(1x,a,i4,a11,f6.4,55x,a)') &
        "|", mylambda, ") lambda = ", e_t%lambda(mylambda), "|"
      call wrtout(std_out,message,"COLL")
      call wrtout(ab_out,message,"COLL")
      do itypat=1,e_t%nctypat
        write(message,'(1x,a,6x,a12,i4,4x,a3,5x,2(3x,a4,f6.4,1x,a2),10x,a)') "|",&
          "- Atom type ", e_t%index_typat(itypat), "->", &
          "U = ",e_t%U_input(itypat)*e_t%lambda(mylambda), "Ha", &
          "J = ",e_t%J_input(itypat)*e_t%lambda(mylambda), "Ha","|"
          call wrtout(std_out,message,"COLL")
          call wrtout(ab_out,message,"COLL")
      end do
      write(message,'(1x,78a)') "+",(/ ("-",itypat=1,76) /), "+"
      call wrtout(std_out,message,"COLL")
      call wrtout(ab_out,message,"COLL")
      call pawpuxinit(dt%dmatpuopt,dt%exchmix,dt%f4of2_sla,dt%f6of2_sla,&
&        is_dfpt,jpawu,dt%lexexch,dt%lpawu,dt%ntypat,pawang,dt%pawprtvol,&
&        pawrad,pawtab,upawu,dt%usedmft,dt%useexexch,dt%usepawu)
      ABI_DEALLOCATE(upawu)
      ABI_DEALLOCATE(jpawu)
    end if

  end function entropyDMFT_nextLambda
!!***

!!****f* ABINIT/m_entropyDMFT/entropyDMFT_addIntegrand
!! NAME
!!  entropyDMFT_addIntegrand
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_scfcv
!!
!! CHILDREN
!!
!! SOURCE

  subroutine entropyDMFT_addIntegrand(e_t,dt,energies,data4etot)

!Arguments ------------------------------------
    type(entropyDMFT_t)   , intent(inout) :: e_t
    type(dataset_type) , intent(in   ) :: dt
    type(energies_type), intent(in   ) :: energies
    type(data4entropyDMFT_t)  , intent(in   ) :: data4etot
!Local variables ------------------------------
    integer :: optdc
    integer :: ictypat
    integer :: iatom, icatom
    integer :: iflavor1, iflavor2, ndim
    integer :: icouple

    if ( e_t%action == AC_NOTHING ) return

    ! Write lambda for restart
    if ( e_t%rank == 0 ) then
      open(unit=e_t%ofile,file=e_t%filedata,position="append",form="unformatted")
      write(e_t%ofile) e_t%lambda(e_t%mylambda)
      close(e_t%ofile)
    end if

    if ( e_t%mylambda == 1 ) then
      ! Save entropy and internal energy for U=0
      e_t%entropy0 = energies%entropy
      ! 1 is for usepaw that is 1 in DMFT, optdc is to know if the DC scheme is
      ! calculated.
      call energies_eval_eint(energies,dt,1,optdc,e_t%energies(E_DIRECT,E_U0),e_t%energies(E_DC,E_U0))
      if ( e_t%rank == 0 ) then
        open(unit=e_t%ofile,file=e_t%filename,position="append")
        write(e_t%ofile,'(a)') "# Temperature [Ha]:"
        write(e_t%ofile,'(es22.14)') e_t%temp
        write(e_t%ofile,'(a)') "# Entropy for lambda=0 [kb]:"
        write(e_t%ofile,'(es22.14)') e_t%entropy0
        write(e_t%ofile,'(a)') "# Internal energy for lambda=0 [Ha]:"
        write(e_t%ofile,'(es22.14)') e_t%energies(E_DC,E_U0)
        close(e_t%ofile)
        open(unit=e_t%ofile,file=e_t%filedata,position="append",form="unformatted")
        write(e_t%ofile) e_t%temp
        write(e_t%ofile) e_t%entropy0
        write(e_t%ofile) e_t%energies(E_DC,E_U0)
        close(e_t%ofile)
      end if
    else if ( e_t%mylambda == e_t%nlambda ) then
      ! Save internal energy for U=Umax
      call energies_eval_eint(energies,dt,1,optdc,e_t%energies(E_DIRECT,E_UU),e_t%energies(E_DC,E_UU))
      if ( e_t%rank == 0 ) then
        open(unit=e_t%ofile,file=e_t%filename,position="append")
        write(e_t%ofile,'(a)') "# Internal energy for lambda=1 [Ha]:"
        write(e_t%ofile,'(es22.14)') e_t%energies(E_DC,E_UU)
      end if
      ! Compute U_max,  and uij
      do ictypat = 1, e_t%nctypat
        ndim = 2*(2*e_t%lpawu(ictypat)+1)
        icouple = 0
        !write(*,*) matU
        if ( e_t%rank == 0 ) &
           write(e_t%ofile,'(a,f7.5,1x,a,i4)') "# Interaction Matrix normalized by U=",e_t%U_input(ictypat) , &
           "[Ha] for atom type", e_t%index_typat(ictypat)
        do iflavor1 = 1, ndim
          if ( e_t%rank == 0 ) &
            write(e_t%ofile,'(14(a21,2x))',advance="no") (/ ( "...", iflavor2=1,iflavor1 ) /)
          do iflavor2 = iflavor1+1, ndim
            icouple = icouple + 1
            e_t%uij(icouple,ictypat) = data4etot%hu_dens(iflavor1,iflavor2,e_t%index_typat(ictypat)) &
                                      /e_t%U_input(ictypat)   ! to modify, this is just the idea
            if ( e_t%rank == 0 ) &
              write(e_t%ofile,'(14(es21.14,2x))',advance="no") e_t%uij(icouple,ictypat)
          end do
          if ( e_t%rank == 0 ) write(e_t%ofile,*)
        end do
      end do

      if ( e_t%rank == 0 ) then
        close(e_t%ofile)
        open(unit=e_t%ofile,file=e_t%filedata,position="append",form="unformatted")
        write(e_t%ofile) e_t%energies(E_DC,E_UU)
        do ictypat = 1, e_t%nctypat
          ndim = 2*(2*e_t%lpawu(ictypat)+1)
          icouple = 0
          do iflavor1 = 1, ndim
            do iflavor2 = iflavor1+1, ndim
              icouple = icouple + 1
              write(e_t%ofile) e_t%uij(icouple,ictypat)
            end do
          end do
        end do
        close(e_t%ofile)
      end if
    endif

    ! For all lambda
    ! Save Docc, Nup and Ndwn
    if ( e_t%rank == 0 ) then
      open(unit=e_t%ofile,file=e_t%filedata,position="append",form="unformatted")
    end if
    do icatom = 1, e_t%ncatom
      iatom = e_t%index_atom(icatom)
      ndim = 2*(2*e_t%lpawu(e_t%typat(icatom))+1)
      icouple = 0
      e_t%e_dc(icatom,e_t%mylambda) = data4etot%e_dc(iatom)
      if ( e_t%rank == 0 ) then
        write(e_t%ofile) e_t%e_dc(icatom,e_t%mylambda)
      end if
      do iflavor1 = 1, ndim
        do iflavor2 = iflavor1+1, ndim
          icouple = icouple + 1
          e_t%docc(icouple,icatom,e_t%mylambda) = data4etot%Docc(iflavor1,iflavor2,iatom)
          if ( e_t%rank == 0 ) then
            write(e_t%ofile) e_t%docc(icouple,icatom,e_t%mylambda)
          end if
        end do
      end do
    end do

    if ( e_t%rank == 0 ) then
      close(e_t%ofile)
    end if

    !call entropyDMFT_computeIntegrand(e_t)

  end subroutine entropyDMFT_addIntegrand
!!***

!!****f* ABINIT/m_entropyDMFT/entropyDMFT_computeEntropy
!! NAME
!!  entropyDMFT_computeEntropy
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_scfcv
!!
!! CHILDREN
!!
!! SOURCE

  subroutine entropyDMFT_computeEntropy(e_t,entropy)

!Arguments ------------------------------------
    type(entropyDMFT_t), intent(inout) :: e_t
    real(dp)        , intent(inout) :: entropy !inout in case we do nothing, avoid to change entropy in NaN
!Local variables ------------------------------
    integer :: ilambda
    integer :: nlambda
    integer :: icatom
    integer :: ncatom
    integer :: icouple
    integer :: ndim
    real(dp), allocatable :: integrand(:,:)
    real(dp) :: docc
    real(dp) :: integral
    real(dp) :: entropyDirect
    character(len=500) :: string
    character(len=50) :: i2str

    if ( e_t%action == AC_NOTHING ) return

    nlambda = e_t%nlambda
    ncatom  = e_t%ncatom

    ABI_ALLOCATE(integrand,(1:nlambda,1:ncatom))
    integrand(1:nlambda,1:ncatom) = zero

    if ( e_t%rank == 0 ) then
      open(unit=e_t%ofile,file=e_t%filename,position="append")
      write(e_t%ofile,'(2a)') ch10,"#Decomposition for each lambda per atom"
    end if
    do ilambda = 1, nlambda
      if ( e_t%rank == 0 ) then
        write(e_t%ofile,'(a,i2,a)') "# ------| lambda ",ilambda," |------ #"
        write(e_t%ofile,'(a5,a23,a23)') "#Atom","DC","<uij*n_i*n_j>"
      end if
      do icatom = 1, ncatom
        integrand(ilambda,icatom) = -e_t%e_dc(icatom,ilambda)
        ndim = 2*(2*e_t%lpawu(e_t%typat(e_t%index_atom(icatom)))+1) ! nflavors
        ndim = ndim*(ndim-1)/2 ! number of couples
        docc = zero
        do icouple = 1, ndim
          docc = docc + e_t%uij(icouple,e_t%typat(icatom)) * e_t%docc(icouple,icatom,ilambda)
        end do
        integrand(ilambda,icatom) = integrand(ilambda,icatom) + docc
        if ( e_t%rank == 0 ) &
          write(e_t%ofile,'(i5,2es23.14)') e_t%index_atom(icatom),e_t%e_dc(icatom,ilambda),docc
      end do
    end do

    ! print for test purpose
    if ( e_t%rank == 0 ) then
      write(i2str,'(i6)') e_t%ncatom
      string='(a,a7,'//TRIM(ADJUSTL(i2str))//'(12x,"Atom",i6))'
      !open(unit=e_t%ofile,file=e_t%filename,position="append")
      write(e_t%ofile,string) ch10,"#lambda", (/ (e_t%index_atom(icatom),icatom=1,ncatom) /)
      string='(1x,f6.4,'//TRIM(ADJUSTL(i2str))//'(1x,e21.14))'
      do ilambda = 1, nlambda
        write(e_t%ofile,string) e_t%lambda(ilambda), (/ (integrand(ilambda,icatom),icatom=1,ncatom) /)
      end do
      close(e_t%ofile)
    end if

    ! Integrate the integrand for all correlated atoms
    call entropyDMFT_integrate(e_t,integrand,integral)
    ABI_DEALLOCATE(integrand)

    write(string,'(a,1x,78a)') ch10,"+",(/ ("-",ilambda=1,76) /), "+"
    call wrtout(std_out,string,"COLL")
    call wrtout(ab_out,string,"COLL")
    write(string,'(1x,a)') "|             Calculation of entropy within  the DMFT Framework              |"
    call wrtout(std_out,string,"COLL")
    call wrtout(ab_out,string,"COLL")
    write(string,'(1x,40a)') "+",(/ ("- ",ilambda=1,38) /), "+"
    call wrtout(std_out,string,"COLL")

    entropyDirect = e_t%entropy0 + ( e_t%energies(E_DIRECT,E_UU) - e_t%energies(E_DIRECT,E_U0) - integral ) /e_t%temp
    entropy = e_t%entropy0 + ( e_t%energies(E_DC,E_UU) - e_t%energies(E_DC,E_U0) - integral ) /e_t%temp

    write(i2str,'(a)') '(1x,a,19x,a11,es21.14,1x,a4,20x,a)'
    write(string,i2str) "|","Integral = ", integral, "[Ha]", "|"
    call wrtout(std_out,string,"COLL")
    write(string,i2str) "|","E(0)     = ", e_t%energies(E_DC,E_U0), "[Ha]", "|"
    call wrtout(std_out,string,"COLL")
    write(string,i2str) "|","E(U)     = ", e_t%energies(E_DC,E_UU), "[Ha]", "|"
    call wrtout(std_out,string,"COLL")
    write(string,i2str) "|","S(0)     = ", e_t%entropy0, "[kb]", "|"
    call wrtout(std_out,string,"COLL")
    write(string,i2str) "|","S(U)     = ", entropy, "[kb]", "|"
    call wrtout(std_out,string,"COLL")
    write(string,'(1x,40a)') "+",(/ ("- ",ilambda=1,38) /), "+"
    call wrtout(std_out,string,"COLL")

    write(string,'(1x,a,16x,a22,es21.14,1x,a4,12x,a)') "|","-kT*entropy is set to ", -e_t%temp*entropy, "[Ha]", "|"
    call wrtout(std_out,string,"COLL")
    call wrtout(ab_out,string,"COLL")

    write(string,'(1x,78a)') "+",(/ ("-",ilambda=1,76) /), "+"
    call wrtout(std_out,string,"COLL")
    call wrtout(ab_out,string,"COLL")

    if ( entropy < zero ) then
      write(string,'(3a)') "Entropy is negative !!!!",ch10,&
      "It does not make any sense"
      MSG_WARNING(string)
    end if

    if ( abs(entropy-entropyDirect) >= tol3 ) then
      write(string,'(1x,a,1x,f8.3,1x,a,1x,f8.3,2a)') "Difference between Direct and DC entropies is", abs(entropy-entropyDirect), &
        "which is greater than", tol3,ch10,"Action : converge better the DMFT and/or DFT loops"
      MSG_WARNING(string)
    end if



  end subroutine entropyDMFT_computeEntropy
!!***

!!****f* ABINIT/m_entropyDMFT/entropyDMFT_integrate
!! NAME
!!  entropyDMFT_integrate
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_entropyDMFT
!!
!! CHILDREN
!!
!! SOURCE

  subroutine entropyDMFT_integrate(e_t,integrand,integral)

!Arguments ------------------------------------
    type(entropyDMFT_t), intent(inout) :: e_t
    real(dp)        , intent(in   ) :: integrand(:,:)
    real(dp)        , intent(  out) :: integral
!Local variables ------------------------------
    real(dp), allocatable :: ypp(:)
    real(dp), allocatable :: fitx(:)
    real(dp), allocatable :: fity(:)
    real(dp) :: integral1
    !real(dp) :: integral2
    !real(dp) :: dx
    integer :: unit
    integer :: Nfit
    integer :: icatom
    integer :: i

    Nfit = 1000

    ABI_ALLOCATE(ypp,(1:e_t%nlambda))
    ABI_ALLOCATE(fitx,(1:Nfit))
    ABI_ALLOCATE(fity,(1:Nfit))

    integral = zero

    unit = get_unit()
    if ( e_t%rank == 0 ) open(unit=unit,file="fit.log")
    do icatom = 1, e_t%ncatom
      integral1 = zero
      !integral2 = zero
      !dx = e_t%U_input(icatom)/dble(e_t%nlambda-1)

      call spline(e_t%lambda,integrand(:,icatom),e_t%nlambda,&
        (integrand(2,icatom)-integrand(1,icatom))/e_t%lambda(2),&   ! first derivative left
        (integrand(e_t%nlambda,icatom)-integrand(e_t%nlambda-1,icatom))/e_t%lambda(2),& ! first derivative right
        ypp)
      fitx(1:Nfit)=(/ (dble(i-1)/dble(Nfit-1),i=1,Nfit) /)
      call splint(e_t%nlambda,e_t%lambda,integrand(:,icatom),ypp,Nfit,fitx,fity)
      integral1 = (sum(fity(1:Nfit))-(fity(1)+fity(Nfit))*half)*e_t%U_input(e_t%typat(icatom))/dble(Nfit-1)
      if ( e_t%rank == 0 ) then
        do i=1,Nfit
          write(unit,'(2ES21.14)') fitx(i),fity(i)
        end do
      end if
      !Unfortunately I don't trust this function and I don't have time to understand it.
      !call spline_integrate(integral2,Nfit,dx,integrand(:,icatom))
      !write(*,*) integral1, integral2

      !if ( abs(integral2-integral1) >= tol6 ) then
      !  write(msg,'(1x,a,1x,f8.6,1x,a,1x,f8.6,1x,a,i4)') "Difference between two different ways of integration is", abs(integral2-integral1), &
      !    "which is greater than", tol6, "for correlated atom",e_t%index_atom(icatom)
      !  MSG_WARNING(msg)
      !end if

      integral = integral+integral1
    end do
    if ( e_t%rank == 0 ) close(unit)

    ABI_DEALLOCATE(ypp)
    ABI_DEALLOCATE(fitx)
    ABI_DEALLOCATE(fity)

  end subroutine entropyDMFT_integrate
!!***

!!****f* ABINIT/m_entropyDMFT/entropyDMFT_destroy
!! NAME
!!  entropyDMFT_destroy
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_scfcv
!!
!! CHILDREN
!!
!! SOURCE

  subroutine entropyDMFT_destroy(e_t)

!Arguments ------------------------------------
    type(entropyDMFT_t), intent(inout) :: e_t

    if ( allocated(e_t%index_atom) ) then
      ABI_DEALLOCATE(e_t%index_atom)
    endif
    if ( allocated(e_t%index_typat)) then
      ABI_DEALLOCATE(e_t%index_typat)
    endif
    if ( allocated(e_t%typat     ) ) then
      ABI_DEALLOCATE(e_t%typat     )
    endif
    if ( allocated(e_t%lpawu     ) ) then
      ABI_DEALLOCATE(e_t%lpawu)
    endif
    if ( allocated(e_t%U_input   ) ) then
      ABI_DEALLOCATE(e_t%U_input)
    endif
    if ( allocated(e_t%J_input   ) ) then
      ABI_DEALLOCATE(e_t%J_input)
    endif
    if ( allocated(e_t%lambda    ) ) then
      ABI_DEALLOCATE(e_t%lambda)
    endif
    if ( allocated(e_t%docc      ) ) then
      ABI_DEALLOCATE(e_t%docc)
    endif
    if ( allocated(e_t%e_dc      ) ) then
      ABI_DEALLOCATE(e_t%e_dc)
    endif
    if ( allocated(e_t%uij       ) ) then
      ABI_DEALLOCATE(e_t%uij)
    endif
    e_t%isset = .FALSE.
  end subroutine entropyDMFT_destroy
!!***
end module m_entropyDMFT
!!***
