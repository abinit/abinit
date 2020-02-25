!!****m* ABINIT/m_data4entropyDMFT
!! NAME
!!  m_data4entropyDMFT
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2020 ABINIT group (J. Bieder)
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

module m_data4entropyDMFT

  use defs_basis
  use m_errors
  use m_abicore

  implicit none

  private

  public :: data4entropyDMFT_init
  public :: data4entropyDMFT_destroy
  public :: data4entropyDMFT_setDocc            ! Must be call for each lambda
  public :: data4entropyDMFT_setHu              ! Hu density
  public :: data4entropyDMFT_setDc

!!***

!!****t* m_data4entropyDMFT/data4entropyDMFT
!! NAME
!!  data4entropyDMFT
!!
!! FUNCTION
!!  This structured datatype contains the necessary data
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

  type, public :: data4entropyDMFT_t
    logical               :: isset = .false.! Are we initialized ?
    integer               :: maxlpawu       ! maximal value for lpawu
    integer               :: natom          ! number of atoms
    integer               :: ntypat         ! number of types of atoms
    real(dp), allocatable :: docc(:,:,:)    ! double occupation for each atom
    real(dp), allocatable :: J_over_U(:)    ! calculate J/U for each atom
    real(dp), allocatable :: e_dc(:)        ! double counting energy calculated for u=1 and j=u/j
    real(dp), allocatable :: hu_dens(:,:,:) ! interaction matrice in density representation
  end type data4entropyDMFT_t
!!***


contains
!!***

!!****f* ABINIT/m_data4entropyDMFT/data4entropyDMFT_init
!! NAME
!!  data4entropyDMFT_init
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2020 ABINIT group (J. Bieder)
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
!!      gstate
!!
!! CHILDREN
!!
!! SOURCE

subroutine data4entropyDMFT_init(this,natom,typat,lpawu,uset2g,upawu,jpawu)

!Arguments ------------------------------------
  type(data4entropyDMFT_t) , intent(inout) :: this
  integer               , intent(in   ) :: natom
  integer , dimension(:), intent(in   ) :: typat
  integer , dimension(:), intent(in   ) :: lpawu
  logical               , intent(in   ) :: uset2g
  real(dp), dimension(:), intent(in   ) :: upawu
  real(dp), dimension(:), intent(in   ) :: jpawu
!Local variables ------------------------------
  integer :: maxlpawu
  integer :: iatom
  integer :: ilpawu
  integer :: nlpawu
  integer :: ityp
  character(len=500) :: message

  this%natom = natom

  if ( size(typat) .ne. natom ) then
    write(message,'(a,i5,a,a,i5,a)') "Disagreement between number of atoms (",natom,")", &
     " and the number of atom types (",size(typat),")."
    MSG_ERROR(message)
  end if

  this%ntypat = maxval(typat) !!! Carefull This should always work but can we have
  ! one type that is not use (ntypat = 5; typat = 1 2 3 4)?
  if ( this%ntypat .ne. size(upawu) .or. this%ntypat .ne. size(jpawu) ) then
    write(message,'(a)') "Disagreement between size of ntypat,upawu and jpawu"
    MSG_ERROR(message)
  end if

  maxlpawu = -1
  nlpawu = size(lpawu)
  do iatom=1,natom
    ityp=typat(iatom)
    if (ityp.le.0 .or. ityp.gt.nlpawu) then
      write(message,'(a)') "Try to access the lpawu value of an atom type that has not a lpawu value."
      MSG_ERROR(message)
    end if
    ilpawu=lpawu(ityp)
    if(uset2g.and.ilpawu==2) ilpawu=1
    if ( ilpawu > maxlpawu ) maxlpawu = ilpawu
  enddo
  this%maxlpawu = maxlpawu

  ABI_ALLOCATE(this%docc,(1:2*(2*maxlpawu+1),1:2*(2*maxlpawu+1),1:natom))
  this%docc(:,:,:) = zero

  ABI_ALLOCATE(this%hu_dens,(1:2*(2*maxlpawu+1),1:2*(2*maxlpawu+1),1:this%ntypat))
  this%hu_dens(:,:,:) = zero

  ABI_ALLOCATE(this%e_dc,(1:natom))
  this%e_dc(:) = zero

  ABI_ALLOCATE(this%J_over_U,(1:natom))
  this%J_over_U(:) = zero
  do iatom=1,natom
    ityp=typat(iatom) ! no need to check since already done once before
    if ( lpawu(ityp) /= -1 .and. upawu(ityp) /= zero) then
      this%J_over_U(iatom) = jpawu(ityp) / upawu(ityp)
    end if
  enddo

  this%isset = .true.
end subroutine data4entropyDMFT_init
!!***

!!****f* ABINIT/m_data4entropyDMFT/data4entropyDMFT_setDocc
!! NAME
!!  data4entropyDMFT_setDocc
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2020 ABINIT group (J. Bieder)
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
!!      qmc_prep_ctqmc
!!
!! CHILDREN
!!
!! SOURCE

subroutine data4entropyDMFT_setDocc(this,iatom,Docc,Nocc)

!Arguments ------------------------------------
    type(data4entropyDMFT_t), intent(inout) :: this
    integer           , intent(in   ) :: iatom
    real(dp), optional, intent(in   ) :: Docc(:,:) !iflavor,iflavor
    real(dp), optional, intent(in   ) :: Nocc(:)   !iflavor
!Local variables ------------------------------
    integer            :: maxnflavor
    integer            :: iflavor1
    integer            :: iflavor2
    character(len=500) :: message

    if ( .not. this%isset ) then
      MSG_ERROR("data4entropyDMFT type not initialized")
    end if

    if ( iatom .gt. this%natom ) then
      write(message,'(a,i4,a,i4,a)') "Value of iatom (",iatom, &
        ") is greater than the number of atom natom(",this%natom,")."
      MSG_ERROR(message)
    end if

    if ( .not. present(Docc) .and. .not. present(Nocc) ) then
      write(message,'(2a)') "Neither Docc nor Nocc is present to set double", &
      "occupancy. Should have one and only one of those."
      MSG_ERROR(message)
    end if

    if ( present(Docc) .and. present(Nocc) ) then
      write(message,'(2a)') "Both Docc and Nocc are present to set double", &
      "occupancy. Should have one and only one of those."
      MSG_ERROR(message)
    end if

    maxnflavor=2*(2*this%maxlpawu+1)
    if ( present(Docc) ) then
      if ( size(Docc,1) .gt. maxnflavor .or. size(Docc,2) .gt. maxnflavor &
          .or. size(Docc,1) .ne. size(Docc,2) ) then
        write(message,'(a,i2,a,i2,a,i2)') "Problem with Docc shape/size : dim1=",size(Docc,1), &
                              " dim2=",size(Docc,2), " max=", maxnflavor
        MSG_ERROR(message)
      end if
      this%docc(1:size(Docc,1),1:size(Docc,1),iatom) = Docc(:,:)
    else if ( present(Nocc) ) then ! Need to compute n_i*n_j (only used for DFT+U)
      if ( size(Nocc,1) .gt. maxnflavor) then
        write(message,'(a,i2,a,i2)') "Problem with Nocc size : dim1=",size(Nocc,1), &
                              " maxnflavor=", maxnflavor
        MSG_ERROR(message)
      end if

      do iflavor1 = 1, (2*size(Nocc,1)+1)
        do iflavor2 = 1, (2*size(Nocc,1)+1)
          this%docc(iflavor2,iflavor1,iatom) = Nocc(iflavor1)*Nocc(iflavor2)
        end do
      end do
    end if

end subroutine data4entropyDMFT_setDocc
!!***

!!****f* ABINIT/m_data4entropyDMFT/data4entropyDMFT_setHu
!! NAME
!!  data4entropyDMFT_setHu
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2020 ABINIT group (J. Bieder)
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
!!      m_dmft
!!
!! CHILDREN
!!
!! SOURCE

subroutine data4entropyDMFT_setHu(this,itypat,hu)

!Arguments ------------------------------------
    type(data4entropyDMFT_t), intent(inout) :: this
    integer           , intent(in   ) :: itypat
    real(dp)          , intent(in   ) :: hu(:,:)   !iflavor
!Local variables ------------------------------
    integer            :: maxnflavor
    character(len=500) :: message

    if ( .not. this%isset ) then
      MSG_ERROR("data4entropyDMFT type not initialized")
    end if

    if ( itypat .gt. this%ntypat ) then
      write(message,'(a,i4,a,i4,a)') "Value of itypat (",itypat, &
        ") is greater than the number of types of atoms (",this%ntypat,")."
      MSG_ERROR(message)
    end if

    maxnflavor=2*(2*this%maxlpawu+1)
    if ( size(hu,1) .gt. maxnflavor .or. size(hu,1) .ne. size(hu,2) ) then
      write(message,'(a,i2,a,i2,a,i2,a,i2)') "Problem with hu size : dim1=",size(hu,1), &
                            " dim2=", size(hu,2), " max=", maxnflavor
      MSG_ERROR(message)
    end if
    this%hu_dens(1:size(hu,1),1:size(hu,1),itypat) = hu(:,:)

end subroutine data4entropyDMFT_setHu
!!***

!!****f* ABINIT/m_data4entropyDMFT/data4entropyDMFT_setDc
!! NAME
!!  data4entropyDMFT_setHu
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2020 ABINIT group (J. Bieder)
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

subroutine data4entropyDMFT_setDc(this,dc)

!Arguments ------------------------------------
    type(data4entropyDMFT_t) , intent(inout) :: this
    real(dp), dimension(:), intent(in   ) :: dc
!Local variables ------------------------------
    character(len=500) :: message

    if ( .not. this%isset ) then
      MSG_ERROR("data4entropyDMFT type not initialized")
    end if

    if ( size(dc,1) .gt. this%natom ) then
      write(message,'(a,i4,a,i4,a)') "Size of dc (",size(dc,1), &
        ") is greater than the number of atom natom(",this%natom,")."
      MSG_ERROR(message)
    end if

    this%e_dc(:) = dc(:)

end subroutine data4entropyDMFT_setDc
!!***

!!****f* ABINIT/m_data4entropyDMFT/data4etotdmf_destroy
!! NAME
!!  data4entropyDMFT_destroy
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2020 ABINIT group (J. Bieder)
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
!!      gstate
!!
!! CHILDREN
!!
!! SOURCE

subroutine data4entropyDMFT_destroy(this)

  !Arguments ------------------------------------
  type(data4entropyDMFT_t), intent(inout) :: this

  if ( .not. this%isset ) return
  if (allocated(this%docc))  then
    ABI_DEALLOCATE(this%docc)
  endif
  if (allocated(this%J_over_U))  then
    ABI_DEALLOCATE(this%J_over_U)
  endif
  if (allocated(this%e_dc))  then
    ABI_DEALLOCATE(this%e_dc)
  endif
  if (allocated(this%hu_dens))  then
    ABI_DEALLOCATE(this%hu_dens)
  endif
  this%maxlpawu = 0
  this%natom = 0
  this%ntypat = 0
  this%isset = .FALSE.
end subroutine data4entropyDMFT_destroy
!!***
end module m_data4entropyDMFT
!!***
