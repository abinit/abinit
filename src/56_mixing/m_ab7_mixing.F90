!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ab7_mixing
!! NAME
!! m_ab7_mixing
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2017 ABINIT group (XG, DC)
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


module m_ab7_mixing

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_io_tools,  only : open_file

 implicit none

 private
!!***

 integer, parameter, public :: AB7_MIXING_NONE        = 0
 integer, parameter, public :: AB7_MIXING_EIG         = 1
 integer, parameter, public :: AB7_MIXING_SIMPLE      = 2
 integer, parameter, public :: AB7_MIXING_ANDERSON    = 3
 integer, parameter, public :: AB7_MIXING_ANDERSON_2  = 4
 integer, parameter, public :: AB7_MIXING_CG_ENERGY   = 5
 integer, parameter, public :: AB7_MIXING_CG_ENERGY_2 = 6
 integer, parameter, public :: AB7_MIXING_PULAY       = 7

 integer, parameter, public :: AB7_MIXING_POTENTIAL  = 0
 integer, parameter, public :: AB7_MIXING_DENSITY    = 1

 integer, parameter, public :: AB7_MIXING_REAL_SPACE     = 1
 integer, parameter, public :: AB7_MIXING_FOURRIER_SPACE = 2

 type, public :: ab7_mixing_object
    integer :: iscf
    integer :: nfft, nspden, kind, space

    logical :: useprec
    integer :: mffmem
    character(len = fnlen) :: diskCache
    integer :: n_index, n_fftgr, n_pulayit, n_pawmix

    integer, dimension(:), pointer :: i_rhor, i_vtrial, i_vresid, i_vrespc
    real(dp), dimension(:,:,:), pointer :: f_fftgr, f_atm
    real(dp), dimension(:,:), pointer :: f_paw

    ! Private
    integer :: n_atom
    real(dp), pointer :: xred(:,:), dtn_pc(:,:)
 end type ab7_mixing_object

 public :: ab7_mixing_new
 public :: ab7_mixing_deallocate

 public :: ab7_mixing_use_disk_cache
 public :: ab7_mixing_use_moving_atoms
 public :: ab7_mixing_copy_current_step

 public :: ab7_mixing_eval_allocate
 public :: ab7_mixing_eval
 public :: ab7_mixing_eval_deallocate
!!***

contains
!!***


!!****f* m_ab7_mixing/init_
!! NAME
!!  init_
!!
!! FUNCTION
!!  Initialize the object
!!
!! PARENTS
!!      m_ab7_mixing
!!
!! CHILDREN
!!      nullify_
!!
!! SOURCE

subroutine init_(mix)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ab7_mixing_object), intent(out) :: mix
! *************************************************************************

 ! Default values.
 mix%iscf      = AB7_MIXING_NONE
 mix%mffmem    = 1
 mix%n_index   = 0
 mix%n_fftgr   = 0
 mix%n_pulayit = 7
 mix%n_pawmix  = 0
 mix%n_atom    = 0
 mix%useprec   = .true.

 call nullify_(mix)

end subroutine init_
!!***

!!****f* m_ab7_mixing/nullify
!! NAME
!!  nullify_
!!
!! FUNCTION
!!  Nullify the pointers
!!
!! PARENTS
!!      m_ab7_mixing
!!
!! CHILDREN
!!      nullify_
!!
!! SOURCE

subroutine nullify_(mix)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ab7_mixing_object), intent(inout) :: mix
! *************************************************************************

 ! Nullify internal pointers.
 nullify(mix%i_rhor)
 nullify(mix%i_vtrial)
 nullify(mix%i_vresid)
 nullify(mix%i_vrespc)
 nullify(mix%f_fftgr)
 nullify(mix%f_atm)
 nullify(mix%f_paw)

end subroutine nullify_
!!***

!!****f* m_ab7_mixing/ab7_mixing_new
!! NAME
!!  ab7_mixing_new
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      dfpt_scfcv,scfcv
!!
!! CHILDREN
!!      nullify_
!!
!! SOURCE

subroutine ab7_mixing_new(mix, iscf, kind, space, nfft, nspden, &
&  npawmix, errid, errmess, npulayit, useprec)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ab7_mixing_new'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ab7_mixing_object), intent(out) :: mix
 integer, intent(in) :: iscf, kind, space, nfft, nspden, npawmix
 integer, intent(out) :: errid
 character(len = 500), intent(out) :: errmess
 integer, intent(in), optional :: npulayit
 logical, intent(in), optional :: useprec

!Local variables-------------------------------
!scalars
 integer :: ii !, i_stat
 character(len = *), parameter :: subname = "ab7_mixing_new"
! *************************************************************************

 ! Set default values.
 call init_(mix)

 ! Argument checkings.
 if (kind /= AB7_MIXING_POTENTIAL .and. kind /= AB7_MIXING_DENSITY) then
    errid = AB7_ERROR_MIXING_ARG
    write(errmess, '(a,a,a,a)' )ch10,&
         & ' ab7_mixing_set_arrays: ERROR -',ch10,&
         & '  Mixing must be done on density or potential only.'
    return
 end if
 if (space /= AB7_MIXING_REAL_SPACE .and. space /= AB7_MIXING_FOURRIER_SPACE) then
    errid = AB7_ERROR_MIXING_ARG
    write(errmess, '(a,a,a,a)' )ch10,&
         & ' ab7_mixing_set_arrays: ERROR -',ch10,&
         & '  Mixing must be done in real or Fourrier space only.'
    return
 end if
 if (iscf /= AB7_MIXING_EIG .and. iscf /= AB7_MIXING_SIMPLE .and. &
      & iscf /= AB7_MIXING_ANDERSON .and. &
      & iscf /= AB7_MIXING_ANDERSON_2 .and. &
      & iscf /= AB7_MIXING_CG_ENERGY .and. &
      & iscf /= AB7_MIXING_PULAY .and. &
      & iscf /= AB7_MIXING_CG_ENERGY_2) then
    errid = AB7_ERROR_MIXING_ARG
    write(errmess, "(A,I0,A)") "Unknown mixing scheme (", iscf, ")."
    return
 end if
 errid = AB7_NO_ERROR

 ! Mandatory arguments.
 mix%iscf     = iscf
 mix%kind     = kind
 mix%space    = space
 mix%nfft     = nfft
 mix%nspden   = nspden
 mix%n_pawmix = npawmix

 ! Optional arguments.
 if (present(useprec)) mix%useprec = useprec

 ! Set-up internal dimensions.
 !These arrays are needed only in the self-consistent case
 if (iscf == AB7_MIXING_EIG) then
    !    For iscf==1, five additional vectors are needed
    !    The index 1 is attributed to the old trial potential,
    !    The new residual potential, and the new
    !    preconditioned residual potential receive now a temporary index
    !    The indices number 4 and 5 are attributed to work vectors.
    mix%n_fftgr=5 ; mix%n_index=1
 else if(iscf == AB7_MIXING_SIMPLE) then
    !    For iscf==2, three additional vectors are needed.
    !    The index number 1 is attributed to the old trial vector
    !    The new residual potential, and the new preconditioned
    !    residual potential, receive now a temporary index.
    mix%n_fftgr=3 ; mix%n_index=1
    if (.not. mix%useprec) mix%n_fftgr = 2
 else if(iscf == AB7_MIXING_ANDERSON) then
    !    For iscf==3 , four additional vectors are needed.
    !    The index number 1 is attributed to the old trial vector
    !    The new residual potential, and the new and old preconditioned
    !    residual potential, receive now a temporary index.
    mix%n_fftgr=4 ; mix%n_index=2
    if (.not. mix%useprec) mix%n_fftgr = 3
 else if (iscf == AB7_MIXING_ANDERSON_2) then
    !    For iscf==4 , six additional vectors are needed.
    !    The indices number 1 and 2 are attributed to two old trial vectors
    !    The new residual potential, and the new and two old preconditioned
    !    residual potentials, receive now a temporary index.
    mix%n_fftgr=6 ; mix%n_index=3
    if (.not. mix%useprec) mix%n_fftgr = 5
 else if(iscf == AB7_MIXING_CG_ENERGY .or. iscf == AB7_MIXING_CG_ENERGY_2) then
    !    For iscf==5 or 6, ten additional vectors are needed
    !    The index number 1 is attributed to the old trial vector
    !    The index number 6 is attributed to the search vector
    !    Other indices are attributed now. Altogether ten vectors
    mix%n_fftgr=10 ; mix%n_index=3
 else if(iscf == AB7_MIXING_PULAY) then
    !    For iscf==7, lot of additional vectors are needed
    !    The index number 1 is attributed to the old trial vector
    !    The index number 2 is attributed to the old residual
    !    The indices number 2 and 3 are attributed to two old precond. residuals
    !    Other indices are attributed now.
    if (present(npulayit)) mix%n_pulayit = npulayit
    mix%n_fftgr=2+2*mix%n_pulayit ; mix%n_index=1+mix%n_pulayit
    if (.not. mix%useprec) mix%n_fftgr = 1+2*mix%n_pulayit
 end if ! iscf cases

 ! Allocate new arrays.
 !allocate(mix%i_rhor(mix%n_index), stat = i_stat)
 !call memocc_abi(i_stat, mix%i_rhor, 'mix%i_rhor', subname)
 ABI_DATATYPE_ALLOCATE(mix%i_rhor,(mix%n_index))
 mix%i_rhor(:)=0
 !allocate(mix%i_vtrial(mix%n_index), stat = i_stat)
 !call memocc_abi(i_stat, mix%i_vtrial, 'mix%i_vtrial', subname)
 ABI_DATATYPE_ALLOCATE(mix%i_vtrial,(mix%n_index))
 mix%i_vtrial(:)=0
 !allocate(mix%i_vresid(mix%n_index), stat = i_stat)
 !call memocc_abi(i_stat, mix%i_vresid, 'mix%i_vresid', subname)
 ABI_DATATYPE_ALLOCATE(mix%i_vresid,(mix%n_index))
 mix%i_vresid(:)=0
 !allocate(mix%i_vrespc(mix%n_index), stat = i_stat)
 !call memocc_abi(i_stat, mix%i_vrespc, 'mix%i_vrespc', subname)
 ABI_DATATYPE_ALLOCATE(mix%i_vrespc,(mix%n_index))
 mix%i_vrespc(:)=0

 ! Setup initial values.
 if (iscf == AB7_MIXING_EIG) then
    mix%i_vtrial(1)=1 ; mix%i_vresid(1)=2 ; mix%i_vrespc(1)=3
 else if(iscf == AB7_MIXING_SIMPLE) then
    mix%i_vtrial(1)=1 ; mix%i_vresid(1)=2 ; mix%i_vrespc(1)=3
    if (.not. mix%useprec) mix%i_vrespc(1)=2
 else if(iscf == AB7_MIXING_ANDERSON) then
    mix%i_vtrial(1)=1 ; mix%i_vresid(1)=2
    if (mix%useprec) then
       mix%i_vrespc(1)=3 ; mix%i_vrespc(2)=4
    else
       mix%i_vrespc(1)=2 ; mix%i_vrespc(2)=3
    end if
 else if (iscf == AB7_MIXING_ANDERSON_2) then
    mix%i_vtrial(1)=1 ; mix%i_vtrial(2)=2
    mix%i_vresid(1)=3
    if (mix%useprec) then
       mix%i_vrespc(1)=4 ; mix%i_vrespc(2)=5 ; mix%i_vrespc(3)=6
    else
       mix%i_vrespc(1)=3 ; mix%i_vrespc(2)=4 ; mix%i_vrespc(3)=5
    end if
 else if(iscf == AB7_MIXING_CG_ENERGY .or. iscf == AB7_MIXING_CG_ENERGY_2) then
    mix%n_fftgr=10 ; mix%n_index=3
    mix%i_vtrial(1)=1
    mix%i_vresid(1)=2 ; mix%i_vresid(2)=4 ; mix%i_vresid(3)=7
    mix%i_vrespc(1)=3 ; mix%i_vrespc(2)=5 ; mix%i_vrespc(3)=8
    mix%i_rhor(2)=9 ; mix%i_rhor(3)=10
 else if(iscf == AB7_MIXING_PULAY) then
    do ii=1,mix%n_pulayit
       mix%i_vtrial(ii)=2*ii-1 ; mix%i_vrespc(ii)=2*ii
    end do
    mix%i_vrespc(mix%n_pulayit+1)=2*mix%n_pulayit+1
    mix%i_vresid(1)=2*mix%n_pulayit+2
    if (.not. mix%useprec) mix%i_vresid(1)=2
 end if ! iscf cases

end subroutine ab7_mixing_new
!!***

!!****f* m_ab7_mixing/ab7_mixing_use_disk_cache
!! NAME
!!  ab7_mixing_use_disk_cache
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!  Obsolete?
!!
!! PARENTS
!!      dfpt_scfcv,scfcv
!!
!! CHILDREN
!!      nullify_
!!
!! SOURCE

subroutine ab7_mixing_use_disk_cache(mix, fnametmp_fft)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ab7_mixing_use_disk_cache'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ab7_mixing_object), intent(inout) :: mix
 character(len = *), intent(in) :: fnametmp_fft
! *************************************************************************

 if (len(trim(fnametmp_fft)) > 0) then
    mix%mffmem = 0
    write(mix%diskCache, "(A)") fnametmp_fft
 else
    mix%mffmem = 1
 end if

end subroutine ab7_mixing_use_disk_cache
!!***


!!****f* m_ab7_mixing/ab7_mixing_use_moving_atoms
!! NAME
!!  ab7_mixing_use_moving_atoms
!!
!! FUNCTION
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
!!      newrho,newvtr
!!
!! CHILDREN
!!      nullify_
!!
!! SOURCE

subroutine ab7_mixing_use_moving_atoms(mix, natom, xred, dtn_pc)

!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ab7_mixing_use_moving_atoms'
!End of the abilint section

 type(ab7_mixing_object), intent(inout) :: mix
 integer, intent(in) :: natom
 real(dp), intent(in), target :: dtn_pc(3, natom)
 real(dp), intent(in), target :: xred(3, natom)

! *************************************************************************

 mix%n_atom = natom
 mix%dtn_pc => dtn_pc
 mix%xred => xred

end subroutine ab7_mixing_use_moving_atoms
!!***


!!****f* m_ab7_mixing/ab7_mixing_copy_current_step
!! NAME
!!  ab7_mixing_copy_current_step
!!
!! FUNCTION
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
!!      dfpt_newvtr,newrho,newvtr
!!
!! CHILDREN
!!      nullify_
!!
!! SOURCE
subroutine ab7_mixing_copy_current_step(mix, arr_resid, errid, errmess, &
&  arr_respc, arr_paw_resid, arr_paw_respc, arr_atm)

!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ab7_mixing_copy_current_step'
!End of the abilint section

 type(ab7_mixing_object), intent(inout) :: mix
 real(dp), intent(in) :: arr_resid(mix%space * mix%nfft, mix%nspden)
 integer, intent(out) :: errid
 character(len = 500), intent(out) :: errmess
 real(dp), intent(in), optional :: arr_respc(mix%space * mix%nfft, mix%nspden)
 real(dp), intent(in), optional :: arr_paw_resid(mix%n_pawmix), arr_paw_respc(mix%n_pawmix)
 real(dp), intent(in), optional :: arr_atm(3, mix%n_atom)
! *************************************************************************

 if (.not. associated(mix%f_fftgr)) then
    errid = AB7_ERROR_MIXING_ARG
    write(errmess, '(a,a,a,a)' )ch10,&
         & ' ab7_mixing_set_arr_current_step: ERROR -',ch10,&
         & '  Working arrays not yet allocated.'
    return
 end if
 errid = AB7_NO_ERROR

 mix%f_fftgr(:,:,mix%i_vresid(1)) = arr_resid(:,:)
 if (present(arr_respc)) mix%f_fftgr(:,:,mix%i_vrespc(1)) = arr_respc(:,:)
 if (present(arr_paw_resid)) mix%f_paw(:, mix%i_vresid(1)) = arr_paw_resid(:)
 if (present(arr_paw_respc)) mix%f_paw(:, mix%i_vrespc(1)) = arr_paw_respc(:)
 if (present(arr_atm)) mix%f_atm(:,:, mix%i_vresid(1)) = arr_atm(:,:)

end subroutine ab7_mixing_copy_current_step
!!***


!!****f* m_ab7_mixing/ab7_mixing_eval_allocate
!! NAME
!!  ab7_mixing_eval_allocate
!!
!! FUNCTION
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
!!      dfpt_newvtr,newrho,newvtr
!!
!! CHILDREN
!!      nullify_
!!
!! SOURCE

subroutine ab7_mixing_eval_allocate(mix, istep)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ab7_mixing_eval_allocate'
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ab7_mixing_object), intent(inout) :: mix
 integer, intent(in), optional :: istep

!Local variables-------------------------------
!scalars
 integer :: istep_,temp_unit !, i_stat
 real(dp) :: tsec(2)
 character(len = *), parameter :: subname = "ab7_mixing_eval_allocate"
 character(len=500) :: msg

! *************************************************************************

 istep_ = 1
 if (present(istep)) istep_ = istep

 ! Allocate work array.
 if (.not. associated(mix%f_fftgr)) then
   !allocate(mix%f_fftgr(mix%space * mix%nfft,mix%nspden,mix%n_fftgr), stat = i_stat)
   !call memocc_abi(i_stat, mix%f_fftgr, 'mix%f_fftgr', subname)
   ABI_ALLOCATE(mix%f_fftgr,(mix%space * mix%nfft,mix%nspden,mix%n_fftgr))
   mix%f_fftgr(:,:,:)=zero
   if (mix%mffmem == 0 .and. istep_ > 1) then
     call timab(83,1,tsec)
     if (open_file(mix%diskCache,msg,newunit=temp_unit,form='unformatted',status='old') /= 0) then
       MSG_ERROR(msg)
     end if
     rewind(temp_unit)
     read(temp_unit) mix%f_fftgr
     if (mix%n_pawmix == 0) close(unit=temp_unit)
     call timab(83,2,tsec)
   end if
 end if
 ! Allocate PAW work array.
 if (.not. associated(mix%f_paw)) then
    !allocate(mix%f_paw(mix%n_pawmix,mix%n_fftgr), stat = i_stat)
    !call memocc_abi(i_stat, mix%f_paw, 'mix%f_paw', subname)
    ABI_ALLOCATE(mix%f_paw,(mix%n_pawmix,mix%n_fftgr))
    if (mix%n_pawmix > 0) then
      mix%f_paw(:,:)=zero
      if (mix%mffmem == 0 .and. istep_ > 1) then
        read(temp_unit) mix%f_paw
        close(unit=temp_unit)
        call timab(83,2,tsec)
      end if
    end if
 end if
 ! Allocate atom work array.
 if (.not. associated(mix%f_atm)) then
    !allocate(mix%f_atm(3,mix%n_atom,mix%n_fftgr), stat = i_stat)
    !call memocc_abi(i_stat, mix%f_atm, 'mix%f_atm', subname)
    ABI_ALLOCATE(mix%f_atm,(3,mix%n_atom,mix%n_fftgr))
 end if

 end subroutine ab7_mixing_eval_allocate
!!***


!!****f* m_ab7_mixing/ab7_mixing_eval_deallocate
!! NAME
!!  ab7_mixing_eval_deallocate
!!
!! FUNCTION
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
!!      dfpt_newvtr,newrho,newvtr
!!
!! CHILDREN
!!      nullify_
!!
!! SOURCE

 subroutine ab7_mixing_eval_deallocate(mix)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ab7_mixing_eval_deallocate'
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ab7_mixing_object), intent(inout) :: mix

!Local variables-------------------------------
!scalars
 integer :: temp_unit !i_all, i_stat
 real(dp) :: tsec(2)
 character(len = *), parameter :: subname = "ab7_mixing_eval_deallocate"
 character(len=500) :: msg

! *************************************************************************

 ! Save on disk and deallocate work array in case on disk cache only.
 if (mix%mffmem == 0) then
    call timab(83,1,tsec)
    if (open_file(mix%diskCache,msg,newunit=temp_unit,form='unformatted',status='unknown') /= 0) then
      MSG_ERROR(msg)
    end if
    rewind(temp_unit)
    ! VALGRIND complains not all of f_fftgr_disk is initialized
     write(temp_unit) mix%f_fftgr
    if (mix%n_pawmix > 0) then
      write(temp_unit) mix%f_paw
    end if
    close(unit=temp_unit)
    call timab(83,2,tsec)
    ABI_DEALLOCATE(mix%f_fftgr)
    nullify(mix%f_fftgr)
    if (associated(mix%f_paw)) then
       ABI_DEALLOCATE(mix%f_paw)
       nullify(mix%f_paw)
    end if
 end if

end subroutine ab7_mixing_eval_deallocate
!!***


!!****f* m_ab7_mixing/ab7_mixing_eval
!! NAME
!!  ab7_mixing_eval
!!
!! FUNCTION
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
!!      dfpt_newvtr,newrho,newvtr
!!
!! CHILDREN
!!      nullify_
!!
!! SOURCE

 subroutine ab7_mixing_eval(mix, arr, istep, nfftot, ucvol, &
& mpi_comm, mpi_summarize, errid, errmess, &
& reset, isecur, pawarr, pawopt, response, etotal, potden, &
& resnrm, comm_atom)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ab7_mixing_eval'
 use interfaces_56_mixing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ab7_mixing_object), intent(inout) :: mix
 integer, intent(in) :: istep, nfftot, mpi_comm
 real(dp), intent(in) :: ucvol
 real(dp), intent(inout) :: arr(mix%space * mix%nfft,mix%nspden)
 logical, intent(in) :: mpi_summarize
 integer, intent(out) :: errid
 character(len = 500), intent(out) :: errmess
 logical, intent(in), optional :: reset
 integer, intent(in), optional :: isecur, comm_atom, pawopt, response
 real(dp), intent(inout), optional :: pawarr(mix%n_pawmix)
 real(dp), intent(in), optional :: etotal
 real(dp), intent(in), optional :: potden(mix%space * mix%nfft,mix%nspden)
 real(dp), intent(out), optional :: resnrm

!Local variables-------------------------------
!scalars
 integer :: moveAtm, dbl_nnsclo, initialized, isecur_
 integer :: usepaw, pawoptmix_, response_
 real(dp) :: resnrm_
! *************************************************************************

 ! Argument checkings.
 if (mix%iscf == AB7_MIXING_NONE) then
    errid = AB7_ERROR_MIXING_ARG
    write(errmess, '(a,a,a,a)' )ch10,&
         & ' ab7_mixing_eval: ERROR -',ch10,&
         & '  No method has been chosen.'
    return
 end if
 if (mix%n_pawmix > 0 .and. .not. present(pawarr)) then
    errid = AB7_ERROR_MIXING_ARG
    write(errmess, '(a,a,a,a)' )ch10,&
         & ' ab7_mixing_eval: ERROR -',ch10,&
         & '  PAW is used, but no pawarr argument provided.'
    return
 end if
 if (mix%n_atom > 0 .and. (.not. associated(mix%dtn_pc) .or. .not. associated(mix%xred))) then
    errid = AB7_ERROR_MIXING_ARG
    write(errmess, '(a,a,a,a)' )ch10,&
         & ' ab7_mixing_eval: ERROR -',ch10,&
         & '  Moving atoms is used, but no xred or dtn_pc attributes provided.'
    return
 end if
 errid = AB7_NO_ERROR

 ! Miscellaneous
 moveAtm = 0
 if (mix%n_atom > 0) moveAtm = 1
 initialized = 1
 if (present(reset)) then
    if (reset) initialized = 0
 end if
 isecur_ = 0
 if (present(isecur)) isecur_ = isecur
 usepaw = 0
 if (mix%n_pawmix > 0) usepaw = 1
 pawoptmix_ = 0
 if (present(pawopt)) pawoptmix_ = pawopt
 response_ = 0
 if (present(response)) response_ = response

 ! Do the mixing.
 resnrm_ = 0.d0
 if (mix%iscf == AB7_MIXING_EIG) then
    !  This routine compute the eigenvalues of the SCF operator
    call scfeig(istep, mix%space * mix%nfft, mix%nspden, &
         & mix%f_fftgr(:,:,mix%i_vrespc(1)), arr, &
         & mix%f_fftgr(:,:,1), mix%f_fftgr(:,:,4:5), errid, errmess)
 else if (mix%iscf == AB7_MIXING_SIMPLE .or. &
      & mix%iscf == AB7_MIXING_ANDERSON .or. &
      & mix%iscf == AB7_MIXING_ANDERSON_2 .or. &
      & mix%iscf == AB7_MIXING_PULAY) then
    if (present(comm_atom)) then
      call scfopt(mix%space, mix%f_fftgr,mix%f_paw,mix%iscf,istep,&
         & mix%i_vrespc,mix%i_vtrial, &
         & mpi_comm,mpi_summarize,mix%nfft,mix%n_pawmix,mix%nspden, &
         & mix%n_fftgr,mix%n_index,mix%kind,pawoptmix_,usepaw,pawarr, &
         & resnrm_, arr, errid, errmess, comm_atom=comm_atom)
    else
      call scfopt(mix%space, mix%f_fftgr,mix%f_paw,mix%iscf,istep,&
         & mix%i_vrespc,mix%i_vtrial, &
         & mpi_comm,mpi_summarize,mix%nfft,mix%n_pawmix,mix%nspden, &
         & mix%n_fftgr,mix%n_index,mix%kind,pawoptmix_,usepaw,pawarr, &
         & resnrm_, arr, errid, errmess)
    end if
    !  Change atomic positions
    if((istep==1 .or. mix%iscf==AB7_MIXING_SIMPLE) .and. mix%n_atom > 0)then
       !    GAF: 2009-06-03
       !    Apparently there are not reason
       !    to restrict iscf=2 for ionmov=5
       mix%xred(:,:) = mix%xred(:,:) + mix%dtn_pc(:,:)
    end if
 else if (mix%iscf == AB7_MIXING_CG_ENERGY .or.  mix%iscf == AB7_MIXING_CG_ENERGY_2) then
    !  Optimize next vtrial using an algorithm based
    !  on the conjugate gradient minimization of etotal
    if (.not. present(etotal) .or. .not. present(potden)) then
       errid = AB7_ERROR_MIXING_ARG
       write(errmess, '(a,a,a,a)' )ch10,&
            & ' ab7_mixing_eval: ERROR -',ch10,&
            & '  Arguments etotal or potden are missing for CG on energy methods.'
       return
    end if
    if (mix%n_atom == 0) then
       ABI_ALLOCATE(mix%xred,(3,0))
       ABI_ALLOCATE(mix%dtn_pc,(3,0))
    end if
    call scfcge(mix%space,dbl_nnsclo,mix%dtn_pc,etotal,mix%f_atm,&
         & mix%f_fftgr,initialized,mix%iscf,isecur_,istep,&
         & mix%i_rhor,mix%i_vresid,mix%i_vrespc,moveAtm,&
         & mpi_comm,mpi_summarize,mix%n_atom,mix%nfft,nfftot,&
         & mix%nspden,mix%n_fftgr,mix%n_index,mix%kind,&
         & response_,potden,ucvol,arr,mix%xred, errid, errmess)
    if (mix%n_atom == 0) then
       ABI_DEALLOCATE(mix%xred)
       ABI_DEALLOCATE(mix%dtn_pc)
    end if
    if (dbl_nnsclo == 1) errid = AB7_ERROR_MIXING_INC_NNSLOOP
 end if

 if (present(resnrm)) resnrm = resnrm_

end subroutine ab7_mixing_eval
!!***


!!****f* m_ab7_mixing/ab7_mixing_deallocate
!! NAME
!!  ab7_mixing_deallocate
!!
!! FUNCTION
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
!!      dfpt_scfcv,scfcv
!!
!! CHILDREN
!!      nullify_
!!
!! SOURCE

subroutine ab7_mixing_deallocate(mix)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ab7_mixing_deallocate'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ab7_mixing_object), intent(inout) :: mix

!Local variables-------------------------------
!scalars
 character(len = *), parameter :: subname = "ab7_mixing_deallocate"
! *************************************************************************

 if (associated(mix%i_rhor)) then
    ABI_DATATYPE_DEALLOCATE(mix%i_rhor)
 end if
 if (associated(mix%i_vtrial)) then
    ABI_DATATYPE_DEALLOCATE(mix%i_vtrial)
 end if
 if (associated(mix%i_vresid)) then
    ABI_DATATYPE_DEALLOCATE(mix%i_vresid)
 end if
 if (associated(mix%i_vrespc)) then
    ABI_DATATYPE_DEALLOCATE(mix%i_vrespc)
 end if
 if (associated(mix%f_fftgr)) then
    ABI_DEALLOCATE(mix%f_fftgr)
 end if
 if (associated(mix%f_paw)) then
    ABI_DEALLOCATE(mix%f_paw)
 end if
 if (associated(mix%f_atm)) then
    ABI_DEALLOCATE(mix%f_atm)
 end if

 call nullify_(mix)

end subroutine ab7_mixing_deallocate
!!***

end module m_ab7_mixing
!!***
