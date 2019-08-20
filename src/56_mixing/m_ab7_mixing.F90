!!****m* ABINIT/m_ab7_mixing
!! NAME
!! m_ab7_mixing
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group (XG, DC, GMR)
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
 use m_abicore
 use m_errors
 use m_linalg_interfaces
 use m_xmpi

 use m_time,      only : timab
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

!!****f* m_ab7_mixing/scfcge
!!
!! NAME
!! scfcge
!!
!! FUNCTION
!! Compute the next vtrial of the SCF cycle.
!! Uses a conjugate gradient minimization of the total energy
!! Can move only the trial potential (if moved_atm_inside==0), or
!! move the trial atomic positions as well (if moved_atm_inside==1).
!!
!! INPUTS
!!  cplex= if 1, real space functions on FFT grid are REAL, if 2, COMPLEX
!!  dtn_pc(3,natom)=preconditioned change of atomic position, in reduced
!!    coordinates. Will be quickly transferred to f_atm(:,:,i_vrespc(1))
!!  etotal=the actual total energy
!!  initialized= if 0, the initialization of the gstate run is not yet finished
!!  iscf =5 => SCF cycle, CG based on estimation of energy gradient
!!       =6 => SCF cycle, CG based on true minimization of the energy
!!  isecur=level of security of the computation
!!  istep= number of the step in the SCF cycle
!!  moved_atm_inside: if==1, the atoms are allowed to move.
!!  mpicomm=the mpi communicator used for the summation
!!  mpi_summarize=set it to .true. if parallelisation is done over FFT
!!  natom=number of atoms
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nfftot=total number of FFT grid points
!!  nspden=number of spin-density components
!!  n_fftgr=third dimension of the array f_fftgr
!!  n_index=dimension for indices of potential/density (see i_vresid, ivrespc, i_rhor...)
!!  opt_denpot= 0 vtrial (and also f_fftgr) really contains the trial potential
!!              1 vtrial (and also f_fftgr) actually contains the trial density
!!  response= if 0, GS calculation, if 1, RF calculation, intrinsically harmonic !
!!  rhor(cplex*nfft,nspden)=actual density
!!  ucvol=unit cell volume in bohr**3
!!
!! OUTPUT
!! dbl_nnsclo=1 if nnsclo has to be doubled to secure the convergence.
!!
!! SIDE EFFECTS
!! Input/Output:
!!  vtrial(cplex*nfft,nspden)= at input, it is the trial potential that gave
!!       the input residual of the potential and Hellman-Feynman forces
!!                       at output, it is the new trial potential .
!!  xred(3,natom)=(needed if moved_atm_inside==1)
!!      reduced dimensionless atomic coordinates
!!      at input, those that generated the input residual of the potential
!!      and Hellman-Feynman forces, at output, these are the new ones.
!!  f_fftgr(cplex*nfft,nspden,n_fftgr)=different functions defined on the fft grid :
!!   The input vtrial is transferred, at output, in f_fftgr(:,:,1).
!!   The input f_fftgr(:,:,i_vresid(1)) contains the last residual.
!!     the value of i_vresid(1) is transferred to i_vresid(2) at output.
!!   The input f_fftgr(:,:,i_vresid(2)) contains the old residual.
!!     the value of i_vresid(2) is transferred to i_vresid(3) at output.
!!   The input f_fftgr(:,:,i_vresid(3)) contains the previous last residual.
!!   For the preconditioned potential residual, the same logic as for the
!!     the potential residual is used, with i_vrespc replacing i_vresid.
!!   The input rhor is transferred, at output, in f_fft(:,:,i_rhor(2)).
!!   The old density is input in f_fft(:,:,i_rhor(2)), and the value of
!!      i_rhor(2) is transferred to i_rhor(3) before the end of the routine.
!!   The input/output search vector is stored in f_fftgr(:,:,6)
!!  f_atm(3,natom,n_fftgr)=different functions defined for each atom :
!!   The input xred is transferred, at output, in f_atm(:,:,1).
!!   The input f_atm(:,:,i_vresid(1)) contains minus the HF forces.
!!     the value of i_vresid(1) is transferred to i_vresid(2) at output.
!!   The input f_atm(:,:,i_vresid(2)) contains minus the old HF forces.
!!     the value of i_vresid(2) is transferred to i_vresid(3) at output.
!!   The input f_atm(:,:,i_vresid(3)) contains minus the previous old HF forces.
!!   For the preconditioned change of atomic positions, the same logic as for the
!!     the potential residual is used, with i_vrespc replacing i_vresid.
!!   The input/output search vector is stored in f_atm(:,:,6)
!!  i_rhor(2:3)=index of the density (past and previous past) in the array f_fftgr
!!  i_vresid(3)=index of the residual potentials (present, past and previous
!!   past) in the array f_fftgr; also similar index for minus Hellman-Feynman
!!   forces in the array f_atm .
!!  i_vrespc(3)=index of the preconditioned residual potentials
!!                  (present, past and previous past) in the array f_fftgr ;
!!   also similar index for the preconditioned change of atomic position (dtn_pc).
!!
!! TODO
!! This routine is much too difficult to read ! Should be rewritten ...
!! Maybe make separate subroutines for line search and CG step ?!
!!
!! PARENTS
!!      m_ab7_mixing
!!
!! CHILDREN
!!      aprxdr,findminscf,sqnormm_v,wrtout
!!
!! SOURCE

subroutine scfcge(cplex,dbl_nnsclo,dtn_pc,etotal,f_atm,&
& f_fftgr,initialized,iscf,isecur,istep,&
& i_rhor,i_vresid,i_vrespc,moved_atm_inside,mpicomm,mpi_summarize,&
& natom,nfft,nfftot,nspden,n_fftgr,n_index,opt_denpot,response,rhor,ucvol,vtrial,xred,errid,errmess)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,initialized,iscf,isecur,istep,moved_atm_inside,mpicomm
 integer,intent(in) :: n_fftgr,n_index,natom,nfft,nfftot,nspden,opt_denpot,response
 integer,intent(out) :: dbl_nnsclo, errid
 character(len = 500), intent(out) :: errmess
 logical, intent(in) :: mpi_summarize
 real(dp),intent(in) :: etotal,ucvol
!arrays
 integer,intent(inout) :: i_rhor(n_index),i_vresid(n_index),i_vrespc(n_index)
 real(dp),intent(in) :: dtn_pc(3,natom),rhor(cplex*nfft,nspden)
 real(dp),intent(inout) :: f_atm(3,natom,n_fftgr)
 real(dp),intent(inout) :: f_fftgr(cplex*nfft,nspden,n_fftgr)
 real(dp),intent(inout) :: vtrial(cplex*nfft,nspden),xred(3,natom)

!Local variables-------------------------------
!mlinmin gives the maximum number of steps in the line minimization
!   after which the algorithm is restarted (with a decrease of the
!   adaptative trial step length). This number should not be large,
!   since if the potential landscape is harmonic, the number of
!   search steps should be small. If it is large, we are not in the
!   harmonic region, and the CG algorithm will not be really useful,
!   so one can just restart the algorithm ...
!scalars
 integer,parameter :: mlinmin=5
 integer,save :: end_linmin,iline_cge,ilinear,ilinmin,isecur_eff,nlinear
 integer,save :: number_of_restart,status
 integer :: choice,iatom,idir,ifft,iline_cge_input,ilinmin_input,isp
 integer :: testcg,tmp,errid_
 real(dp),save :: d2edv2_old2,d_lambda_old2,dedv_old2,etotal_old
 real(dp),save :: etotal_previous,lambda_adapt,lambda_new,lambda_old,resid_old
 real(dp) :: d2e11,d2e12,d2e22,d2edv2_new,d2edv2_old
 real(dp) :: d2edv2_predict,d_lambda,de1,de2,dedv_mix
 real(dp) :: dedv_new,dedv_old,dedv_predict,determ,etotal_input
 real(dp) :: etotal_predict,gamma,lambda_input,lambda_predict2
 real(dp) :: lambda_predict=1.0_dp,ratio,reduction
 real(dp) :: resid_input,temp
 character(len=500) :: message
!arrays
 real(dp) :: resid_new(1)
 real(dp), allocatable :: tmp_fft1(:,:)

! *************************************************************************

 errid = AB7_NO_ERROR
 dbl_nnsclo = 0

!reduction gives the level of reduction of the error in
!the line minimization to be reached for the minimization to be
!considered successfull
 reduction=0.1_dp

!nlinear increases with the number of times the 2D minimization succeded
!to reach the true minimum directly. It is a measure of the
!degree of parabolicity of the problem, and is used to
!skip some steps by performing extrapolation.
 if(istep==1)then

!  Skipping some steps is sometimes unsecure, so it is possible
!  to make nlinear start at a negative value - if isecur is positive
   isecur_eff=isecur
   nlinear=min(-isecur_eff,0)
   ilinear=0

!  Response function calculation are intrinsically harmonic, so one
!  can shift isecur (by -2), and start with a positive nlinear
   if(response==1)then
     isecur_eff=isecur-2
     nlinear=-isecur_eff
     ilinear=nlinear
   end if

   iline_cge=0
   ilinmin=0
 end if

!Compute actual residual resid_new (residual of f_fftgr(:,:,i_vrespc(1))
 call sqnormm_v(cplex,i_vrespc(1),mpicomm,mpi_summarize,1,nfft,resid_new,n_fftgr,nspden,opt_denpot,f_fftgr)

!Save input residual and ilinmin for final printing
 resid_input=resid_new(1)
 etotal_input=etotal
 ilinmin_input=ilinmin
 iline_cge_input=iline_cge
!Transfer dtn_pc in f_atm
 if(moved_atm_inside==1)then
   f_atm(:,:,i_vrespc(1))=dtn_pc(:,:)
 end if

!=======================================================================
!Now the routine is decomposed in three mutually exclusive parts :
!if(istep==1)then initialize the algorithm
!else if(ilinmin>0)then perform the line minimisation
!else if(ilinmin==0)then determine the new search direction (CG step)
!=======================================================================


!--------------------------------------
!Here initialize the algorithm
 if(istep==1)then

!  At the beginning of each gstate run, lambda_adapt is forced to have the
!  same value, that is 1.0_dp. In the other cases when istep=1 (at different
!  broyden steps, for example), the previously obtained
!  adaptive value is kept.
   if(initialized==0)lambda_adapt=1.0_dp
   lambda_old=0.0_dp
   lambda_input=0.0_dp
   number_of_restart=0
   lambda_new=lambda_adapt

   f_fftgr(:,:,1)=vtrial(:,:)
   f_fftgr(:,:,i_rhor(2))=rhor(:,:)

!  This copy must be written in F77, because of stack problems on the DECs
   do isp=1,nspden
     do ifft=1,cplex*nfft
       f_fftgr(ifft,isp,6)=f_fftgr(ifft,isp,i_vrespc(1))
     end do
   end do
   vtrial(:,:)=f_fftgr(:,:,1)+(lambda_new-lambda_old)*f_fftgr(:,:,6)
   if(moved_atm_inside==1)then
     f_atm(:,:,1)=xred(:,:)
     f_atm(:,:,i_rhor(2))=xred(:,:)
!    There shouldn t be problems with the stack size for this small array.
     f_atm(:,:,6)=f_atm(:,:,i_vrespc(1))
     xred(:,:)=f_atm(:,:,1)+(lambda_new-lambda_old)*f_atm(:,:,6)
   end if
   tmp=i_vrespc(2) ; i_vrespc(2)=i_vrespc(1) ; i_vrespc(1)=tmp
   tmp=i_vresid(2) ; i_vresid(2)=i_vresid(1) ; i_vresid(1)=tmp
   ilinmin=1
   resid_old=resid_new(1)
   etotal_old=etotal

   status=0

!  --------------------------------------

!  Here performs the line minimisation
 else if(ilinmin>0)then

   lambda_input=lambda_new

!  The choice with the Brent algorithm has been abandoned in version 1.6.m

!  Compute the approximate energy derivatives dedv_new and dedv_old,
!  from vresid and vresid_old
   choice=2
   call aprxdr(cplex,choice,dedv_mix,dedv_new,dedv_old,&
&   f_atm,f_fftgr,i_rhor(2),i_vresid,moved_atm_inside,mpicomm,mpi_summarize,&
&   natom,nfft,nfftot,nspden,n_fftgr,rhor,ucvol,xred)
   d_lambda=lambda_new-lambda_old
   dedv_old=dedv_old/d_lambda
   dedv_new=dedv_new/d_lambda

!  DEBUG
!  write(std_out,'(a,4es12.4,i3)' )' scfcge:lold,lnew,dold,dnew,status',  &
!  &  lambda_old,lambda_new,dedv_old,dedv_new,status
!  ENDDEBUG

   if(status==0 .or. status==3)then
!
!    Then, compute a predicted point along the line
!    The value of choice determines the minimization algorithm
!    choice=1 uses the two values of the derivative of the energy
!    choice=2 uses the two values of the energy, and and estimate of the
!    second derivative at the mid-point.

     choice=1
     if(iscf==6)choice=2
     call findminscf(choice,dedv_new,dedv_old,dedv_predict,&
&     d2edv2_new,d2edv2_old,d2edv2_predict,&
&     etotal,etotal_old,etotal_predict,&
&     lambda_new,lambda_old,lambda_predict,errid_,message)
     if (errid_ /= AB7_NO_ERROR) then
       call wrtout(std_out,message,'COLL')
     end if

!    Suppress the next line for debugging  (there is another such line)
     status=0

!    DEBUG
!    Keep this debugging feature : it gives access to the investigation of lines
!    in a different approach
!    if(response==1 .and. istep>8)then
!    lambda_predict=1.2d-2
!    if(istep>=15)lambda_predict=lambda_predict-0.002
!    if(istep>=14)stop
!    status=3
!    end if
!    ENDDEBUG

   else
     if(status/=-1)then
       status=-1
       lambda_predict=-2.5_dp
     else
       lambda_predict=lambda_predict+0.1_dp
     end if
   end if

!  If the predicted point is very close to the most recent
!  computed point, while this is the first trial on this line,
!  then we are in the linear regime :
!  nlinear is increased by one unit. For the time being, do this even when
!  moved_atm_inside==1 (the code still works when it is done, but it
!  seems to be a bit unstable). The maximal value of nlinear is 1, except
!  when isecur_eff is a negative number, less than -1.
   if( abs(lambda_predict-lambda_new)/&
&   (abs(lambda_predict)+abs(lambda_new)) < 0.01 .and. ilinmin==1  ) then
!    if(moved_atm_inside==0 .and. nlinear<max(1,-isecur_eff) )nlinear=nlinear+1
     if(nlinear<max(1,-isecur_eff))nlinear=nlinear+1
     ilinear=nlinear
   end if

!  If the predicted point is close to the most recent computed point,
!  or the previous one, set on the flag of end of line minization
   end_linmin=0
   if(abs(lambda_new-lambda_predict)*2.0_dp&
&   /(abs(lambda_predict)+abs(lambda_new)) <reduction) end_linmin=1
   if(abs(lambda_old-lambda_predict)*2.0_dp&
&   /(abs(lambda_predict)+abs(lambda_new)) <reduction) end_linmin=1

   if(status/=0)end_linmin=0

!  Save the closest old lambda, if needed,
!  also examine the reduction of the interval, and eventual stop
!  the present line minimisation, because of convergence (end_linmin=1)
!  Also treat the case in which the predicted value of lambda is negative,
!  or definitely too small in which case the algorithm has to be restarted
!  (not a very good solution, though ...)
!  Finally also treat the case where insufficiently converged
!  density at lambda=0.0_dp happens, which screws up the line minimisation.

!  Here restart the algorithm with the best vtrial.
!  Also make reduction in lambda_adapt
!  DEBUG
!  write(std_out,*)' scfcge : status=',status
!  ENDDEBUG
   if( end_linmin==0 .and. status==0 .and.                               &
&   (  (lambda_predict<0.005_dp*lambda_adapt .and. iscf==5)     .or.  &
&   (abs(lambda_predict)<0.005_dp*lambda_adapt .and. iscf==6).or.  &
&   ilinmin==mlinmin                                      )     )then
     if(number_of_restart>12)then
       errid = AB7_ERROR_MIXING_CONVERGENCE
       write(errmess,'(a,a,i0,a,a,a,a,a)')&
&       'Potential-based CG line minimization not',' converged after ',number_of_restart,' restarts. ',ch10,&
&       'Action : read the eventual warnings about lack of convergence.',ch10,&
&       'Some might be relevant. Otherwise, raise nband. Returning'
       MSG_WARNING(errmess)
       return
     end if
!    Make reduction in lambda_adapt (kind of steepest descent...)
     write(message,'(a,a,a)')&
&     'Potential-based CG line minimization has trouble to converge.',ch10,&
&     'The algorithm is restarted with more secure parameters.'
     MSG_WARNING(message)
     number_of_restart=number_of_restart+1
!    At the second restart, double the number of non-self consistent loops.
     if(number_of_restart>=2)dbl_nnsclo=1
     lambda_adapt=lambda_adapt*0.7_dp
     lambda_new=lambda_adapt
!    If the last energy is better than the old one, transfer the data.
!    Otherwise, no transfer must occur (very simple to code...)
     if(etotal<etotal_old .or. abs(lambda_old)<1.0d-8)then
       f_fftgr(:,:,1)=vtrial(:,:)
       f_fftgr(:,:,i_rhor(2))=rhor(:,:)
       do isp=1,nspden
         do ifft=1,cplex*nfft
           f_fftgr(ifft,isp,6)=f_fftgr(ifft,isp,i_vrespc(1))
         end do
       end do
       if(moved_atm_inside==1)then
         f_atm(:,:,1)=xred(:,:)
         f_atm(:,:,i_rhor(2))=xred(:,:)
         f_atm(:,:,6)=f_atm(:,:,i_vrespc(1))
       end if
       tmp=i_vrespc(2) ; i_vrespc(2)=i_vrespc(1) ; i_vrespc(1)=tmp
       tmp=i_vresid(2) ; i_vresid(2)=i_vresid(1) ; i_vresid(1)=tmp
       resid_old=resid_new(1)
       etotal_old=etotal
     end if
     lambda_old=0.0_dp
     ilinmin=1
!    Putting the flag to -1 avoids the usual actions taken with end_linmin=1
     end_linmin=-1
!    Also put ilinear and nlinear to 0
     ilinear=0
     nlinear=0

!    Here lambda_new is the closest to lambda_predict,
!    or lambda_old is still 0.0_dp, while the energy shows that the minimum
!    is away from 0.0_dp (insufficiently converged density at lambda=0.0_dp).
   else if( abs(lambda_new-lambda_predict)<abs(lambda_old-lambda_predict) &
&     .or.                                                           &
&     ( abs(lambda_old)<1.0d-6 .and.                               &
&     ilinmin>1              .and.                               &
&     etotal>etotal_previous         )                           &
&     )then
     f_fftgr(:,:,1)=vtrial(:,:)
     tmp=i_rhor(3) ; i_rhor(3)=i_rhor(2) ; i_rhor(2)=tmp
     f_fftgr(:,:,i_rhor(2))=rhor(:,:)
     tmp=i_vrespc(3) ; i_vrespc(3)=i_vrespc(2)
     i_vrespc(2)=i_vrespc(1); i_vrespc(1)=tmp;
     tmp=i_vresid(3); i_vresid(3)=i_vresid(2)
     i_vresid(2)=i_vresid(1) ; i_vresid(1)=tmp
     if(moved_atm_inside==1)then
       f_atm(:,:,1)=xred(:,:)
       f_atm(:,:,i_rhor(2))=xred(:,:)
     end if
     d_lambda_old2=lambda_old-lambda_new
     lambda_old=lambda_new
     etotal_old=etotal
     resid_old=resid_new(1)
     d2edv2_old2=d2edv2_new
     dedv_old=dedv_new
     dedv_old2=dedv_new
!    if(abs(lambda_new-lambda_predict)*2.0_dp&
!    &    /abs(lambda_new+lambda_predict)        <reduction) end_linmin=1

!    Here lambda_old is the closest to lambda_predict (except for avoiding
!    lambda_old==0.0_dp)
   else
     tmp=i_vresid(3) ; i_vresid(3)=i_vresid(1) ; i_vresid(1)=tmp
     f_fftgr(:,:,i_rhor(3))=rhor(:,:)
     if(moved_atm_inside==1) f_atm(:,:,i_rhor(3))=xred(:,:)
     tmp=i_vrespc(3) ; i_vrespc(3)=i_vrespc(1) ; i_vrespc(1)=tmp
     d_lambda_old2=lambda_new-lambda_old
     etotal_previous=etotal
     d2edv2_old2=d2edv2_old
     dedv_old2=dedv_old
!    if(abs(lambda_old-lambda_predict)*2.0_dp&
!    &    /abs(lambda_old+lambda_predict)        <reduction) end_linmin=1
   end if

!  If the interval has not yet been sufficiently reduced,
!  continue the search
   if(end_linmin==0)then
     lambda_new=lambda_predict

!    DEBUG
!    write(std_out,'(a,2es16.6)' )&
!    &   ' scfcge : continue search, lambda_old,lambda_new=',lambda_old,lambda_new
!    write(std_out,'(a,2es16.6)' )&
!    &   ' scfcge : f_fftgr(3:4,1,1)=',f_fftgr(3:4,1,1)
!    write(std_out,'(a,2es16.6)' )&
!    &   ' scfcge : f_fftgr(3:4,1,6)=',f_fftgr(3:4,1,6)
!    ENDDEBUG

     vtrial(:,:)=f_fftgr(:,:,1)+(lambda_new-lambda_old)*f_fftgr(:,:,6)
     if(moved_atm_inside==1)then
       xred(:,:)=f_atm(:,:,1)+(lambda_new-lambda_old)*f_atm(:,:,6)
     end if

     ilinmin=ilinmin+1
!
!    Here generates a starting point for next line search
   else
     iline_cge=iline_cge+1
     if(end_linmin==1)ilinmin=0
     lambda_old=0.0_dp

!    In order to generate the new step, take into account previous
!    optimal lambdas (including those of previous ion moves),
!    and the selected new one, if it is positive.
!    However, wait iline_cge>1 to select new ones.
!    lambda_adapt has been initialized at 1.0_dp
     if(iline_cge>1 .and. lambda_new>0.0_dp )then
!      Actually compute a geometric mean
       lambda_adapt= ( lambda_adapt**(dble(iline_cge-1)) * abs(lambda_new)) &
&       **(1.0_dp/dble(iline_cge))
!      In order to recover the previous algorithm, it is enough
!      to decomment the next line
!      lambda_adapt=1.0_dp
     end if
     lambda_new=lambda_adapt

     vtrial(:,:)=f_fftgr(:,:,1)+lambda_new*f_fftgr(:,:,i_vrespc(2))
     if(moved_atm_inside==1)then
       xred(:,:)=f_atm(:,:,1)+lambda_new*f_atm(:,:,i_vrespc(2))
     end if

!    End choice between continue line minim and determine new direction
   end if

!
!  -------------------------------

!  Here perform the CG step

 else if(ilinmin==0)then

!  Compute the approximate energy derivatives dedv_mix,dedv_new,dedv_old
   choice=3
   call aprxdr(cplex,choice,dedv_mix,dedv_new,dedv_old,&
&   f_atm,f_fftgr,i_rhor(2),i_vresid,moved_atm_inside,mpicomm,mpi_summarize,&
&   natom,nfft,nfftot,nspden,n_fftgr,rhor,ucvol,xred)

   dedv_mix=dedv_mix/lambda_new
   dedv_new=dedv_new/lambda_new
   dedv_old=dedv_old/lambda_new

!  DEBUG
!  write(message, '(a,3es12.4)' )' scfcge: lambda_adapt',&
!  &     lambda_adapt
!  call wrtout(std_out,message,'COLL')

!  write(message, '(a,3es12.4)' )' scfcge: dedv_old,dedv_new,dedv_mix',&
!  &     dedv_old,dedv_new,dedv_mix
!  call wrtout(std_out,message,'COLL')
!  ENDDEBUG

!  Then, compute a predicted point, either along the line,
!  or in a 2D plane
   testcg=1
   if(testcg==0)then
!    This part corresponds to steepest descent,
!    in which the line minimisation can be done
!    using different algorithms, varying with the value of choice
     choice=1
     if(iscf==6)choice=2
     call findminscf(choice,dedv_new,dedv_old,dedv_predict,&
&     d2edv2_new,d2edv2_old,d2edv2_predict,&
&     etotal,etotal_old,etotal_predict,&
&     lambda_new,lambda_old,lambda_predict,errid_,message)
     if (errid_ /= AB7_NO_ERROR) then
       call wrtout(std_out,message,'COLL')
     end if
     lambda_predict2=0.0_dp
!    Suppress the next line for debugging (there is another such line)
     status=0
   else
!    This part corresponds to conjugate gradient
!    A 2D minimisation is performed
!    oldest direction is labelled 2
!    newest direction is labelled 1
     de1=dedv_old ;  de2=dedv_old2
     d2e11=(dedv_new-dedv_old)/lambda_new
     d2e22=d2edv2_old2
     d2e12=(dedv_mix-dedv_old)/d_lambda_old2
!    The system to be solved is
!    0 = de1 + lambda1 d2e11 + lambda2 d2d12
!    0 = de2 + lambda1 d2e12 + lambda2 d2d22
     determ=d2e11*d2e22-d2e12*d2e12
     lambda_predict=-(de1*d2e22-de2*d2e12)/determ
     lambda_predict2=(de1*d2e12-de2*d2e11)/determ
     d2edv2_new=d2e11 ;  d2edv2_old=d2e11
   end if

!  DEBUG
!  write(message, '(a,5es11.3)' )' scfcge: de1,de2,d2e11,d2e22,d2e12',&
!  &               de1,de2,d2e11,d2e22,d2e12
!  call wrtout(std_out,message,'COLL')
!  write(std_out,'(a,2es12.4)' )' scfcge: la_predict,la_predict2',&
!  &               lambda_predict,lambda_predict2
!  -----
!  write(std_out,*)'residues ',
!  !$       de1+lambda_predict*d2e11+lambda_predict2*d2e12,
!  !$       de2+lambda_predict*d2e12+lambda_predict2*d2e22
!  if(.true.)stop
!  ENDDEBUG
!

!  Determine the region of the 2D search space
!  in which the predicted point is located,
!  or use linear indicator to decide interpolation
!  and advance to next 2D search.
   end_linmin=0
   write(message, '(a,2i3)' )' nlinear, ilinear',nlinear,ilinear
   call wrtout(std_out,message,'COLL')
   if(lambda_predict<0.0_dp)then
!    Something is going wrong. Just take a reasonable step
!    along the steepest descent direction (Region III).
!    Actually, Region I and region III are treated in the same way later.
!    In effect, this corresponds to restart the algorithm
     end_linmin=3
!    Also put ilinear and nlinear to 0
     ilinear=0
     nlinear=0
!    Decrease the adaptive step to predict next direction
     lambda_adapt=lambda_adapt*0.7_dp
   else if(ilinear>=1) then
!    Region IV : will do an interpolation
     end_linmin=4
     ilinear=ilinear-1
   else if(abs(lambda_predict2)>reduction          .or.&
&     lambda_predict<0.5_dp                .or.&
&     lambda_predict>2.5_dp                .or.&
&     lambda_predict-abs(lambda_predict2)/reduction <0.0_dp  ) then
!    Region II : lambda_predict is not too good, and not too bad.
     end_linmin=2
   else if (abs(1.0_dp-lambda_predict)<reduction)then
!    Region I, the out-of-line point is OK.
     end_linmin=1
   else
!    If everything fails, then region II.
     end_linmin=2
   end if

!  DEBUG
!  write(message, '(a,2es12.4,i2)' )&
!  &     ' scfcge : la_predict, la_predict2, region',&
!  &       lambda_predict,lambda_predict2,end_linmin
!  call wrtout(std_out,message,'COLL')
!  ENDDEBUG

!  Treat region I, in the same way as region III
   if(end_linmin==1 .or. end_linmin==3)then

!    In region I, the line search is
!    along vtrial-vtrial_old.
!    The closest point is the new point
!    thus to be transfered in the "old" locations

     do isp=1,nspden
       do ifft=1,cplex*nfft
         f_fftgr(ifft,isp,6)=(vtrial(ifft,isp)-f_fftgr(ifft,isp,1))/lambda_new
       end do
     end do
     f_fftgr(:,:,1)=vtrial(:,:)
     f_fftgr(:,:,i_rhor(2))=rhor(:,:)
     if(moved_atm_inside==1)then
       f_atm(:,:,6)=(xred(:,:)-f_atm(:,:,1))/lambda_new
       f_atm(:,:,1)=xred(:,:)
       f_atm(:,:,i_rhor(2))=xred(:,:)
     end if
     tmp=i_vrespc(2) ; i_vrespc(2)=i_vrespc(1) ; i_vrespc(1)=tmp
     tmp=i_vresid(3) ; i_vresid(3)=i_vresid(2)
     i_vresid(2)=i_vresid(1) ; i_vresid(1)=tmp
     d_lambda_old2=-lambda_new
     lambda_old=lambda_new
     etotal_old=etotal
     resid_old=resid_new(1)
     d2edv2_old=d2edv2_new
     dedv_old=dedv_new

!    Region I or III : one is close of the 2D minimum,
!    or lambda_predict was negative (indicate a problem of convergence)
!    Compute next trial potential along the
!    PC residual and not along this search direction.
     ilinmin=0
!    Question : isn t it here that one should prevent region I to called
!    itself more than 1 time ???
!    Here the small difference between region I and region III
     if(end_linmin==3)ilinmin=1
     lambda_old=0.0_dp
     lambda_new=lambda_adapt

     vtrial(:,:)=f_fftgr(:,:,1)+lambda_new*f_fftgr(:,:,i_vrespc(2))
     if(moved_atm_inside==1)then
       xred(:,:)=f_atm(:,:,1)+lambda_new*f_atm(:,:,i_vrespc(2))
     end if
!    The new vtrial has been generated

   else

!    Here region II or IV
     ilinmin=1
     if (lambda_predict==0._dp) then
       gamma=zero
     else
       gamma=lambda_predict2/lambda_predict
     end if
!    Compute new search direction and trial potential
     write(message,*)' compute new search direction '
     call wrtout(std_out,message,'COLL')
     do isp=1,nspden
       do ifft=1,cplex*nfft
         f_fftgr(ifft,isp,6)=(vtrial(ifft,isp)-f_fftgr(ifft,isp,1))/lambda_new+ &
&         gamma*f_fftgr(ifft,isp,6)
       end do
     end do
     vtrial(:,:)=f_fftgr(:,:,1)+ lambda_predict*f_fftgr(:,:,6)
     if(moved_atm_inside==1)then
       f_atm(:,:,6)=(xred(:,:)-f_atm(:,:,1))/lambda_new+ gamma*f_atm(:,:,6)
       xred(:,:)=f_atm(:,:,1)+ lambda_predict*f_atm(:,:,6)
     end if

!    If end_linmin==2, then this vtrial is the good one

     if(end_linmin==2)then

       lambda_old=0.0_dp
       lambda_new=lambda_predict

     else if(end_linmin==4)then

!      predict the result of the computation at the trial potential
!      defined in the end_linmin==2 case
       gamma=lambda_predict2/d_lambda_old2
       ratio=lambda_predict/lambda_new

!      Take care of vtrial
       f_fftgr(:,:,1)=vtrial(:,:)

       ABI_ALLOCATE(tmp_fft1,(cplex*nfft,nspden))
!      Take care of vresid
       tmp_fft1(:,:)=f_fftgr(:,:,i_vresid(2))
       f_fftgr(:,:,i_vresid(2))=tmp_fft1(:,:)&
&       +ratio*(f_fftgr(:,:,i_vresid(1))-tmp_fft1(:,:))&
&       +gamma*(f_fftgr(:,:,i_vresid(3))-tmp_fft1(:,:))
       f_fftgr(:,:,i_vresid(3))=tmp_fft1(:,:)

!      Take care of rhor
       tmp_fft1(:,:)=f_fftgr(:,:,i_rhor(2))
       f_fftgr(:,:,i_rhor(2))=tmp_fft1(:,:)&
&       +ratio*(rhor(:,:)-tmp_fft1(:,:))&
&       +gamma*(f_fftgr(:,:,i_rhor(3))-tmp_fft1(:,:))
       f_fftgr(:,:,i_rhor(3))=tmp_fft1(:,:)

!      Take care of vrespc
       tmp_fft1(:,:)=f_fftgr(:,:,i_vrespc(2))
       f_fftgr(:,:,i_vrespc(2))=tmp_fft1(:,:)&
&       +ratio*(f_fftgr(:,:,i_vrespc(1))-tmp_fft1(:,:))&
&       +gamma*(f_fftgr(:,:,i_vrespc(3))-tmp_fft1(:,:))
       f_fftgr(:,:,i_vrespc(3))=tmp_fft1(:,:)
       ABI_DEALLOCATE(tmp_fft1)

       if(moved_atm_inside==1)then
         do idir=1,3
           do iatom=1,natom

!            Take care of xred
             f_atm(idir,iatom,1)=xred(idir,iatom)

!            Take care of -HF forces
             temp=f_atm(idir,iatom,i_vresid(2))
             f_atm(idir,iatom,i_vresid(2))=f_atm(idir,iatom,i_vresid(2))&
&             +ratio*(f_atm(idir,iatom,i_vresid(1))-f_atm(idir,iatom,i_vresid(2)))&
&             +gamma*(f_atm(idir,iatom,i_vresid(3))-f_atm(idir,iatom,i_vresid(2)))
             f_atm(idir,iatom,i_vresid(3))=temp

!            Take care of old xreds
             temp=f_atm(idir,iatom,i_rhor(2))
             f_atm(idir,iatom,i_rhor(2))=f_atm(idir,iatom,i_rhor(2))&
&             +ratio*(   xred(idir,iatom)          -f_atm(idir,iatom,i_rhor(2)))&
&             +gamma*(f_atm(idir,iatom,i_rhor(3))-f_atm(idir,iatom,i_rhor(2)))
             f_atm(idir,iatom,i_rhor(3))=temp

!            Take care of preconditioned changes of atomic positions
             temp=f_atm(idir,iatom,i_vrespc(2))
             f_atm(idir,iatom,i_vrespc(2))=f_atm(idir,iatom,i_vrespc(2))&
&             +ratio*(f_atm(idir,iatom,i_vrespc(1))-f_atm(idir,iatom,i_vrespc(2)))&
&             +gamma*(f_atm(idir,iatom,i_vrespc(3))-f_atm(idir,iatom,i_vrespc(2)))
             f_atm(idir,iatom,i_vrespc(3))=temp

           end do
         end do
       end if

!      Since we are at the 2D minimum, the derivative is supposed
!      to vanish. Note that dedv_old should not change, by contrast.
       dedv_old2=0.0_dp
       d_lambda_old2=-lambda_predict
       d2edv2_old2=-dedv_old/lambda_predict
       lambda_old=lambda_predict
       ilinmin=0

!      So, jump to the next line
       iline_cge=iline_cge+1
       write(message,*)' energy CG update : after 2D interpolation,'
       call wrtout(std_out,message,'COLL')
       write(message,*)'    computation in the next plane '
       call wrtout(std_out,message,'COLL')
       write(message,*)
       call wrtout(std_out,message,'COLL')
       lambda_old=0.0_dp
       lambda_new=lambda_adapt

       vtrial(:,:)=f_fftgr(:,:,1)+lambda_new*f_fftgr(:,:,i_vrespc(2))
       if(moved_atm_inside==1)then
         xred(:,:)=f_atm(:,:,1)+lambda_new*f_atm(:,:,i_vrespc(2))
       end if

!      The new trial potential is now generated

!      End the specific treatment of region IV
     end if
!
!    End the choice between treatment of region I, II, or IV
   end if

!  End of choice between initialisation or more developed parts of the CG algorithm
 else
   errid = AB7_ERROR_MIXING_ARG
   errmess = 'scfcge : BUG You should not be here ! '
   return
 end if

!--------------------------------------

!Write information : it will be easy to read by typing  grep scfcge logfile

 if(istep==1)then
   write(message,'(a,a,a)') ' scfcge:',ch10,' scfcge:istep-iline_cge-ilinmin lambda      etot             resid '
   call wrtout(std_out,message,'COLL')
 end if

 if(ilinmin_input/=0 .or. istep==1)then
!  Usual line minimisation step

   if(iline_cge_input<10)then
     write(message, '(a,i4,a,i1,a,i1,es13.4,es20.12,es12.4)' )&
&     ' scfcge: actual  ',istep,'-',iline_cge_input,'-',ilinmin_input,lambda_input,etotal_input,resid_input
   else
     write(message, '(a,i3,a,i2,a,i1,es13.4,es20.12,es12.4)' )&
&     ' scfcge: actual  ',istep,'-',iline_cge_input,'-',ilinmin_input,lambda_input,etotal_input,resid_input
   end if
   call wrtout(std_out,message,'COLL')

   if( (end_linmin==1.or.end_linmin==-1) .and. istep/=1 )then

     if(end_linmin==1)then
       write(message, '(a,es13.4,a,i2,a,a)' )&
&       ' scfcge: predict         ',lambda_predict,&
&       ' suff. close => next line, ilinear=',ilinear,ch10,&
&       ' scfcge:'
     else if(end_linmin==-1)then
       write(message, '(a,es13.4,a,a,a)' )&
&       ' scfcge: predict         ',lambda_predict,&
&       ' restart the algorithm ',ch10,&
&       ' scfcge:'
     end if
     call wrtout(std_out,message,'COLL')

     if(iline_cge_input<9)then
       write(message, '(a,i4,a,i1,a,i1,es13.4,es20.12,es12.4)' ) &
&       ' scfcge: start   ',istep,'-',iline_cge,'-',0,0.0,etotal_old,resid_old
     else
       write(message, '(a,i3,a,i2,a,i1,es13.4,es20.12,es12.4)' ) &
&       ' scfcge: start   ',istep,'-',iline_cge,'-',0,0.0,etotal_old,resid_old
     end if
     call wrtout(std_out,message,'COLL')

   else if(istep/=1) then
     write(message, '(a,es13.4,a)' )&
&     ' scfcge: predict         ',lambda_predict,&
&     ' not close enough => continue minim.'
     call wrtout(std_out,message,'COLL')
   end if

 else
!  CG prediction
   if(iline_cge_input<10)then
     write(message, '(a,i4,a,i1,a,es11.4,es20.12,es12.4,a,i1)' )&
&     ' scfcge: actual  ',istep,'-',iline_cge_input,'-off',&
&     lambda_adapt,etotal_input,resid_input,', end=',end_linmin
   else
     write(message, '(a,i3,a,i2,a,es11.4,es20.12,es12.4,a,i1)' )&
&     ' scfcge: actual  ',istep,'-',iline_cge_input,'-off',&
&     lambda_adapt,etotal_input,resid_input,', end=',end_linmin
   end if
   call wrtout(std_out,message,'COLL')

   if(end_linmin==4)then
     write(message, '(a)' ) ' scfcge:'
     call wrtout(std_out,message,'COLL')
   end if

 end if

end subroutine scfcge
!!***

!!****f* ABINIT/scfeig
!! NAME
!! scfeig
!!
!! FUNCTION
!! Compute the largest eigenvalue and eigenvector of the SCF cycle.
!! A brute force algorithm is presently used.
!!
!! INPUTS
!!  istep= number of the step in the SCF cycle
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nspden=number of spin-density components
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  vtrial0(nfft,nspden)= contains vtrial at istep == 1
!!  vtrial(nfft,nspden)= at input, it is the trial potential that gave vresid .
!!       at output, it is an updated trial potential
!!  vrespc(nfft,nspden)=the input preconditioned residual potential
!!  work(nfft,nspden,2)=work space
!!
!! PARENTS
!!      m_ab7_mixing
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine scfeig(istep,nfft,nspden,vrespc,vtrial,vtrial0,work,errid,errmess)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istep,nfft,nspden
 integer,intent(out) :: errid
 character(len = 500), intent(out) :: errmess
!arrays
 real(dp),intent(inout) :: vtrial0(nfft,nspden),work(nfft,nspden,2)
 real(dp),intent(inout) :: vrespc(nfft,nspden)
 real(dp), intent(inout) :: vtrial(nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: ifft,isp
 real(dp) :: eigen_scf,factor,fix_resid,resid_new,resid_old
 character(len=500) :: message

! *************************************************************************

 errid = AB7_NO_ERROR

 if(nspden==4)then
   errid = AB7_ERROR_MIXING_ARG
   write(errmess, *) ' scfeig: does not work yet for nspden=4'
   return
 end if

!Set a fixed residual square for normalization of eigenvectors
 fix_resid=1.0d-4

!A few initialisations for the first istep
 if(istep==1)then

   write(message, '(a,es12.4,a,a,a,a,a,a,a)' )&
&   ' scfeig: fixed PC_residual square =',fix_resid,ch10,&
&   '    Note that fixed resid should always be much larger',ch10,&
&   '    than initial PC resid square, still sufficiently',ch10,&
&   '    small to reduce anharmonic effects ',ch10
   call wrtout(std_out,message,'COLL')

!  Compute the preconditioned residual
   resid_old=0.0_dp
   do isp=1,nspden
     do ifft=1,nfft
       resid_old=resid_old+vrespc(ifft,isp)**2
     end do
   end do
   write(message, '(a,es12.4)' )' scfeig: initial PC_residual square =',resid_old
   call wrtout(std_out,message,'COLL')
   if(resid_old>1.0d-8)then
     errid = AB7_ERROR_MIXING_ARG
     write(errmess,'(a,a,a,a,a,a,a,a,a,a)') ch10,&
&     ' scfeig : ERROR -',ch10,&
&     '  This value is not good enough to allow',ch10,&
&     '  the computation of the eigenvectors of the SCF cycle.',ch10,&
&     '  It should be better than 1.0d-8 .',ch10,&
&     '  Action : improve the accuracy of your starting wavefunctions.'
     return
   end if

!  Also transfer vtrial in vtrial_old
   vtrial0(:,:)=vtrial(:,:)

!  In order to start the search for eigenvectors,
!  use the tiny residual vector, renormalized
   factor=sqrt(fix_resid/resid_old)
   work(:,:,1)=vrespc(:,:)*factor
   vtrial(:,:)=vtrial0(:,:)+work(:,:,1)

!  If istep is not equal to 1
 else if(istep>=2)then
!
!  Compute the corresponding operator expectation value
!  And put the residual vector minus the difference
!  between vtrial and vtrial_old
!  (this is actually the action of the operator !) in vect(*,2)
   eigen_scf=0.0_dp
   do isp=1,nspden
     do ifft=1,nfft
       eigen_scf=eigen_scf+&
&       work(ifft,isp,1) * vrespc(ifft,isp)
     end do
   end do

   do isp=1,nspden
     do ifft=1,nfft
       vrespc(ifft,isp)=vrespc(ifft,isp)&
&       +vtrial(ifft,isp)-vtrial0(ifft,isp)
       work(ifft,isp,2)=vrespc(ifft,isp)
     end do
   end do
   eigen_scf=eigen_scf/fix_resid
   write(message, '(a,es12.4,a)' ) &
&   ' scfeig : Operator expectation value ',eigen_scf,' (extremal eigenvalue * diemix)'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
!
!  Compute residual of vect(*,2)
   resid_new=zero
   do isp=1,min(nspden,2)
     do ifft=1,nfft
       resid_new=resid_new+ work(ifft,isp,2) ** 2
     end do
   end do
   if (nspden==4) then
     do ifft=1,nfft
       resid_new=resid_new+two*(work(ifft,3,2)**2+work(ifft,4,2)**2)
     end do
   end if
   factor=sqrt(fix_resid/resid_new)
   if(eigen_scf<zero) then
     factor=-factor ! the new vector MAY be oposite to the old one
!    if(factor<-one) factor=-factor ! the new vector is not opposed to the old one
   end if
   write(message, '(a,es12.4)' ) &
&   ' scfeig : Inverse of renormalization factor ',one/factor
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
   write(message, '(a,es12.4)' ) &
&   ' scfeig : Convergence criterion value (->0 at convergency) ',one/factor-eigen_scf-one
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   work(:,:,1)=work(:,:,2)*factor
   vtrial(:,:)=vtrial0(:,:)+work(:,:,1)
!  End the different istep cases
 end if

end subroutine scfeig
!!***

!!****f* m_ab7_mixing/scfopt
!!
!! NAME
!! scfopt
!!
!! FUNCTION
!! Compute the next vtrial of the SCF cycle.
!! Possible algorithms are : simple mixing, Anderson (order 1 or 2), Pulay
!!
!! INPUTS
!!  cplex= if 1, real space functions on FFT grid are REAL, if 2, COMPLEX
!!  iscf= 2 => simple mixing
!!      = 3,4 => Anderson mixing
!!      = 7 => Pulay mixing
!!  istep= number of the step in the SCF cycle
!!  mpicomm=the mpi communicator used for the summation
!!  comm_atom=the mpi communicator over atoms ; PAW only (optional argument)
!!  mpi_summarize=set it to .true. if parallelisation is done over FFT
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  npawmix=-PAW only- number of spherical part elements to be mixed
!!  nspden=number of spin-density components
!!  n_fftgr=third dimension of the array f_fftgr
!!  n_index=dimension for indices of potential/density (see ivrespc, i_vtrial...)
!!  opt_denpot= 0 vtrial (and also f_fftgr) really contains the trial potential
!!              1 vtrial (and also f_fftgr) actually contains the trial density
!!  pawoptmix= - PAW only - 1 if the computed residuals include the PAW (rhoij) part
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  vtrial(cplex*nfft,nspden)= at input, it is the trial potential that gave
!!     the input preconditioned residual potential
!!     at output, it is the new trial potential .
!!  f_fftgr(cplex*nfft,nspden,n_fftgr)=different functions defined on the fft grid :
!!   The input vtrial is transferred, at output,in f_fftgr(:,:,i_vtrial(1)).
!!   The old vtrial is transferred, at output,in f_fftgr(:,:,i_vtrial(2)).
!!   The input preconditioned residual potential is in f_fftgr(:,:,i_vrespc(1))
!!   Two input old preconditioned residual potentials in f_fftgr(:,:,i_vrespc(2)) and f_fftgr(:,:,i_vrespc(3))
!!    Before output a permutation of i_vrespc(1), i_vrespc(2) and i_vrespc(3) occurs, without
!!    actually copying all the data (change of pointer).
!!  i_vrespc(n_index)=index of the preconditioned residual potentials (present and past) in the array f_fftgr
!!  i_vtrial(n_index)  =indices of the potential (present and past) in the array f_fftgr
!!  ==== if usepaw==1
!!    f_paw(npawmix,n_fftgr*mffmem*usepaw)=different functions used for PAW
!!                                           (same as f_fftgr but for spherical part)
!!    vpaw(npawmix*usepaw)=at input, the aug. occupancies (rhoij) that gave
!!                               the input preconditioned residual potential
!!                           at output, it is the new aug. occupancies.
!!
!! PARENTS
!!      m_ab7_mixing
!!
!! CHILDREN
!!      dgetrf,dgetri,dotprodm_v,sqnormm_v,wrtout,xmpi_sum
!!
!! SOURCE

subroutine scfopt(cplex,f_fftgr,f_paw,iscf,istep,i_vrespc,i_vtrial,&
& mpicomm,mpi_summarize,nfft,npawmix,nspden,n_fftgr,&
& n_index,opt_denpot,pawoptmix,usepaw,vpaw,vresid,vtrial,errid,errmess, &
& comm_atom) ! optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,iscf,istep,n_fftgr,n_index,nfft
 integer,intent(in) :: npawmix,nspden,opt_denpot,pawoptmix,usepaw,mpicomm
 integer, intent(in),optional :: comm_atom
 integer,intent(out) :: errid
 character(len = 500), intent(out) :: errmess
 logical, intent(in) :: mpi_summarize
 real(dp), intent(out) :: vresid
!arrays
 integer,intent(inout) :: i_vrespc(n_index),i_vtrial(n_index)
 real(dp),intent(inout) :: f_fftgr(cplex*nfft,nspden,n_fftgr)
 real(dp),intent(inout) :: f_paw(npawmix,n_fftgr*usepaw),vpaw(npawmix*usepaw)
 real(dp),intent(inout) :: vtrial(cplex*nfft,nspden)
!Local variables-------------------------------
!scalars
 integer,parameter :: npulaymax=50
 integer :: i_vstore,ierr,ifft,ii,index,isp,jj,comm_atom_,niter,npulay,tmp
 real(dp),save :: prod_resid_old,resid_old,resid_old2
 real(dp) :: aa1,aa2,bb,cc1,cc2,current,det,lambda,lambda2,resid_best
 character(len=500) :: message
!arrays
 integer,allocatable :: ipiv(:)
 real(dp),save :: amat(npulaymax+1,npulaymax+1)
 real(dp) :: mpibuff(2),prod_resid(1),prod_resid2(1),resid_new(1)
 real(dp),allocatable :: alpha(:),amatinv(:,:),amat_paw(:),rwork(:)

! *************************************************************************

!DEBUG
!write(std_out,*)' scfopt : enter ; istep,iscf ',istep,iscf
!ENDDEBUG

 errid = AB7_NO_ERROR

 comm_atom_=xmpi_comm_self; if(present(comm_atom)) comm_atom_=comm_atom

 i_vstore=i_vtrial(1)
 if (iscf==4) i_vstore=i_vtrial(2)
 if (iscf==7) then
   if (modulo(n_fftgr, 2) == 0 ) then
     npulay=(n_fftgr-2)/2
   else
     npulay=(n_fftgr-1)/2
   end if
   i_vstore=i_vtrial(npulay)
 else
   npulay=0
 end if

!Compute the new residual resid_new, from f_fftgr/f_paw(:,:,i_vrespc(1))
 call sqnormm_v(cplex,i_vrespc(1),mpicomm,mpi_summarize,1,nfft,resid_new,n_fftgr,nspden,opt_denpot,f_fftgr)
 if (usepaw==1.and.pawoptmix==1) then
   do index=1,npawmix
     resid_new(1)=resid_new(1)+f_paw(index,i_vrespc(1))**2
   end do
   call xmpi_sum(resid_new(1),comm_atom_,ierr)
 end if
 vresid = resid_new(1)

!_______________________________________________________________
!Here use only the preconditioning, or initialize the other algorithms

 if (istep==1 .or. iscf==2) then
   write(message,'(2a)') ch10,' Simple mixing update:'
   call wrtout(std_out,message,'COLL')

   write(message,*)' residual square of the potential: ',resid_new(1)
   call wrtout(std_out,message,'COLL')

   ! Store information for later use
   if (iscf==3.or.iscf==4) resid_old=resid_new(1)
   if (iscf==7) then
     amat(:,:)=zero
     amat(1,1)=resid_new(1)
   end if

   ! Compute new vtrial (and new rhoij if PAW)
   if (iscf/=2) f_fftgr(:,:,i_vstore)=vtrial(:,:)
   vtrial(:,:)=vtrial(:,:)+f_fftgr(:,:,i_vrespc(1))
   if (usepaw==1) then
     if (iscf/=2) f_paw(:,i_vstore)=vpaw(:)
     vpaw(:)=vpaw(:)+f_paw(:,i_vrespc(1))
   end if

!  _______________________________________________________________
!  Here Anderson algorithm using one previous iteration
 else if((istep==2 .or. iscf==3).and.iscf/=7)then

   write(message,'(2a)') ch10,' Anderson update:'
   call wrtout(std_out,message,'COLL')

   write(message,*)' residual square of the potential: ',resid_new(1)
   call wrtout(std_out,message,'COLL')

!  Compute prod_resid from f_fftgr/f_paw(:,:,i_vrespc(1)) and f_fftgr/f_paw(:,:,i_vrespc(2))
   call dotprodm_v(cplex,1,prod_resid,i_vrespc(1),i_vrespc(2),mpicomm,mpi_summarize,1,1,&
&   nfft,n_fftgr,n_fftgr,nspden,opt_denpot,f_fftgr,f_fftgr)
   if (usepaw==1.and.pawoptmix==1) then
     do index=1,npawmix
       prod_resid(1)=prod_resid(1)+f_paw(index,i_vrespc(1))*f_paw(index,i_vrespc(2))
     end do
     call xmpi_sum(prod_resid(1),comm_atom_,ierr)
   end if

!  Compute mixing factor
   lambda=(resid_new(1)-prod_resid(1))/(resid_new(1)+resid_old-2*prod_resid(1))
   write(message,*)' mixing of old trial potential: ',lambda
   call wrtout(std_out,message,'COLL')

!  Evaluate best residual square on the line
   resid_best=(1.0_dp-lambda)*(1.0_dp-lambda)*resid_new(1)&
&   +(1.0_dp-lambda)*lambda        *2*prod_resid(1)&
&   +lambda        *lambda        *resid_old
   write(message,*)' predicted best residual square on the line: ',resid_best
   call wrtout(std_out,message,'COLL')

!  Store information for later use
   if (iscf==4) then
     prod_resid_old=prod_resid(1)
     resid_old2=resid_old
   end if
   resid_old=resid_new(1)

!  Save latest trial potential and compute new trial potential
   do isp=1,nspden
     do ifft=1,cplex*nfft
       current=vtrial(ifft,isp)
       vtrial(ifft,isp)=(one-lambda)*(current                      +f_fftgr(ifft,isp,i_vrespc(1)))&
&       +lambda      *(f_fftgr(ifft,isp,i_vtrial(1))+f_fftgr(ifft,isp,i_vrespc(2)))
       f_fftgr(ifft,isp,i_vstore)=current
     end do
   end do

!  PAW: save latest rhoij and compute new rhoij
   do index=1,npawmix
     current=vpaw(index)
     vpaw(index)=(one-lambda)*(current                 +f_paw(index,i_vrespc(1)))&
&     +lambda      *(f_paw(index,i_vtrial(1))+f_paw(index,i_vrespc(2)))
     f_paw(index,i_vstore)=current
   end do

!  _______________________________________________________________
!  Here Anderson algorithm using two previous iterations
 else if(iscf==4.and.iscf/=7)then

   write(message,'(2a)') ch10,' Anderson (order 2) update:'
   call wrtout(std_out,message,'COLL')

   write(message,*)' residual square of the potential: ',resid_new(1)
   call wrtout(std_out,message,'COLL')

!  Compute prod_resid from f_fftgr/f_paw(:,:,i_vrespc(1)) and f_fftgr/f_paw(:,:,i_vrespc(2))
   call dotprodm_v(cplex,1,prod_resid,i_vrespc(1),i_vrespc(2),mpicomm,mpi_summarize,1,1,&
&   nfft,n_fftgr,n_fftgr,nspden,opt_denpot,f_fftgr,f_fftgr)
   if (usepaw==1.and.pawoptmix==1) then
     do index=1,npawmix
       prod_resid(1)=prod_resid(1)+f_paw(index,i_vrespc(1))*f_paw(index,i_vrespc(2))
     end do
   end if

!  Compute prod_resid2 from f_fftgr/f_paw(:,:,i_vrespc(1)) and f_fftgr/f_paw(:,:,i_vrespc(3))
   call dotprodm_v(cplex,1,prod_resid2,i_vrespc(1),i_vrespc(3),mpicomm,mpi_summarize,1,1,&
&   nfft,n_fftgr,n_fftgr,nspden,opt_denpot,f_fftgr,f_fftgr)
   if (usepaw==1.and.pawoptmix==1) then
     do index=1,npawmix
       prod_resid2(1)=prod_resid2(1)+f_paw(index,i_vrespc(1))*f_paw(index,i_vrespc(3))
     end do
!    MPI reduction
     mpibuff(1)=prod_resid(1);mpibuff(2)=prod_resid2(1)
     call xmpi_sum(mpibuff,comm_atom_,ierr)
     prod_resid(1)=mpibuff(1);prod_resid2(1)=mpibuff(2)
   end if

!  Compute mixing factors
   aa1=resid_new(1)+resid_old -two*prod_resid (1)
   aa2=resid_new(1)+resid_old2-two*prod_resid2(1)
   bb =resid_new(1)+prod_resid_old-prod_resid(1)-prod_resid2(1)
   cc1=resid_new(1)-prod_resid (1)
   cc2=resid_new(1)-prod_resid2(1)
   det=aa1*aa2-bb*bb
   lambda =(aa2*cc1-bb*cc2)/det
   lambda2=(aa1*cc2-bb*cc1)/det
   write(message,*)' mixing of old trial potentials: ',lambda,lambda2
   call wrtout(std_out,message,'COLL')

!  Store information for later use
   prod_resid_old=prod_resid(1)
   resid_old2=resid_old
   resid_old=resid_new(1)

!  Save latest trial potential and compute new trial potential
   do isp=1,nspden
     do ifft=1,cplex*nfft
       current=vtrial(ifft,isp)
       vtrial(ifft,isp)=&
&       (one-lambda-lambda2)*(current                      +f_fftgr(ifft,isp,i_vrespc(1)))&
&       +lambda             *(f_fftgr(ifft,isp,i_vtrial(1))+f_fftgr(ifft,isp,i_vrespc(2)))&
&       +lambda2            *(f_fftgr(ifft,isp,i_vtrial(2))+f_fftgr(ifft,isp,i_vrespc(3)))
       f_fftgr(ifft,isp,i_vstore)=current
     end do
   end do

!  PAW: save latest rhoij and compute new rhoij
   do index=1,npawmix
     current=vpaw(index)
     vpaw(index)=&
&     (one-lambda-lambda2)*(current                 +f_paw(index,i_vrespc(1)))&
&     +lambda             *(f_paw(index,i_vtrial(1))+f_paw(index,i_vrespc(2)))&
&     +lambda2            *(f_paw(index,i_vtrial(2))+f_paw(index,i_vrespc(3)))
     f_paw(index,i_vstore)=current
   end do

!  _______________________________________________________________
!  Here Pulay algorithm
 else if(iscf==7)then

   niter=min(istep,npulay+1)

   write(message,'(2a,i2,a)') ch10,' Pulay update with ',niter-1,' previous iterations:'
   call wrtout(std_out,message,'COLL')

   if (npulay>npulaymax) then
     errid = AB7_ERROR_MIXING_CONVERGENCE
     write(errmess, '(4a)' ) ch10,&
&     ' scfopt: ERROR - ',ch10,&
&     '  Too many iterations required for Pulay algorithm (<50) !'
     return
   end if

!  Compute "A" matrix
   if (istep>npulay+1) then
     do jj=1,niter-1
       do ii=1,niter-1
         amat(ii,jj)=amat(ii+1,jj+1)
       end do
     end do
   end if
   if (usepaw==1.and.pawoptmix==1) then
     ABI_ALLOCATE(amat_paw,(niter))
     amat_paw(:)=zero
     do ii=1,niter
       do index=1,npawmix
         amat_paw(ii)=amat_paw(ii)+f_paw(index,i_vrespc(1))*f_paw(index,i_vrespc(1+niter-ii))
       end do
     end do
     call xmpi_sum(amat_paw,comm_atom_,ierr)
   end if
   do ii=1,niter
     call dotprodm_v(cplex,1,amat(ii,niter),i_vrespc(1),i_vrespc(1+niter-ii),mpicomm,mpi_summarize,1,1,&
&     nfft,n_fftgr,n_fftgr,nspden,opt_denpot,f_fftgr,f_fftgr)
     if (usepaw==1.and.pawoptmix==1) amat(ii,niter)=amat(ii,niter)+amat_paw(ii)
     if (ii<niter) amat(niter,ii)=amat(ii,niter)
   end do
   if (usepaw==1.and.pawoptmix==1)then
     ABI_DEALLOCATE(amat_paw)
   end if

!  Invert "A" matrix
   ABI_ALLOCATE(amatinv,(niter,niter))
   amatinv(1:niter,1:niter)=amat(1:niter,1:niter)
   ABI_ALLOCATE(ipiv,(niter))
   ABI_ALLOCATE(rwork,(niter))
   call dgetrf(niter,niter,amatinv,niter,ipiv,ierr)
   call dgetri(niter,amatinv,niter,ipiv,rwork,niter,ierr)
   ABI_DEALLOCATE(ipiv)
   ABI_DEALLOCATE(rwork)

!  Compute "alpha" factors
   ABI_ALLOCATE(alpha,(niter))
   alpha=zero
   det=zero
   do ii=1,niter
     do jj=1,niter
       alpha(ii)=alpha(ii)+amatinv(jj,ii)
       det=det+amatinv(jj,ii)
     end do
   end do
   alpha(:)=alpha(:)/det
   ABI_DEALLOCATE(amatinv)
   write(message,'(a,5(1x,g10.3))')' mixing of old trial potential: alpha(m:m-4)=',(alpha(ii),ii=niter,max(1,niter-4),-1)
   call wrtout(std_out,message,'COLL')

!  Save latest trial potential and compute new trial potential
   do isp=1,nspden
     do ifft=1,cplex*nfft
       current=vtrial(ifft,isp)
       vtrial(ifft,isp)=alpha(niter)*(current+f_fftgr(ifft,isp,i_vrespc(1)))
       do ii=niter-1,1,-1
         vtrial(ifft,isp)=vtrial(ifft,isp)+alpha(ii) &
&         *(f_fftgr(ifft,isp,i_vtrial(niter-ii))+f_fftgr(ifft,isp,i_vrespc(1+niter-ii)))
       end do
       f_fftgr(ifft,isp,i_vstore)=current
     end do
   end do

!  PAW: save latest rhoij and compute new rhoij
   do index=1,npawmix
     current=vpaw(index)
     vpaw(index)=alpha(niter)*(current+f_paw(index,i_vrespc(1)))
     do ii=niter-1,1,-1
       vpaw(index)=vpaw(index)+alpha(ii) &
&       *(f_paw(index,i_vtrial(niter-ii))+f_paw(index,i_vrespc(1+niter-ii)))
     end do
     f_paw(index,i_vstore)=current
   end do

   ABI_DEALLOCATE(alpha)
!  _______________________________________________________________
!  End of choice of optimization method
 end if

!Permute potential indices
 if (iscf==3) then
   tmp=i_vrespc(2) ; i_vrespc(2)=i_vrespc(1) ; i_vrespc(1)=tmp
 else if (iscf==4) then
   tmp=i_vrespc(3) ; i_vrespc(3)=i_vrespc(2) ; i_vrespc(2)=i_vrespc(1) ; i_vrespc(1)=tmp
   tmp=i_vtrial(2) ; i_vtrial(2)=i_vtrial(1) ; i_vtrial(1)=tmp
 else if (iscf==7) then
   tmp=i_vtrial(  npulay)
   do ii=  npulay,2,-1
     i_vtrial(ii)=i_vtrial(ii-1)
   end do
   i_vtrial(1)=tmp
   tmp=i_vrespc(1+npulay)
   do ii=1+npulay,2,-1
     i_vrespc(ii)=i_vrespc(ii-1)
   end do
   i_vrespc(1)=tmp
 end if

end subroutine scfopt
!!***

!!****f* ABINIT/findminscf
!! NAME
!! findminscf
!!
!! FUNCTION
!! Compute the minimum of a function whose value
!! and derivative are known at two points, using different algorithms.
!! Also deduce different quantities at this predicted
!! point, and at the two other points
!!
!! INPUTS
!! choice=1,uses a linear interpolation of the derivatives
!!       =2,uses a quadratic interpolation based on the
!!        values of the function, and the second derivative at mid-point
!! etotal_1=first value of the function
!! etotal_2=second value of the function
!! dedv_1=first value of the derivative
!! dedv_2=second value of the derivative
!! lambda_1=first value of the argument
!! lambda_2=second value of the argument
!!
!! OUTPUT
!! dedv_predict=predicted value of the derivative (usually zero,
!!  except if choice=4, if it happens that a minimum cannot be located,
!!  and a trial step is taken)
!! d2edv2_predict=predicted value of the second derivative (not if choice=4)
!! d2edv2_1=first value of the second derivative (not if choice=4)
!! d2edv2_2=second value of the second derivative (not if choice=4)
!! etotal_predict=predicted value of the function
!! lambda_predict=predicted value of the argument
!! status= 0 if everything went normally ;
!!         1 if negative second derivative
!!         2 if some other problem
!!
!! PARENTS
!!      scfcge
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine findminscf(choice,dedv_1,dedv_2,dedv_predict,&
& d2edv2_1,d2edv2_2,d2edv2_predict,&
& etotal_1,etotal_2,etotal_predict,&
& lambda_1,lambda_2,lambda_predict,errid,errmess)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice
 integer,intent(out) :: errid
 character(len=500), intent(out) :: errmess
 real(dp),intent(in) :: dedv_1,dedv_2,etotal_1,etotal_2,lambda_1,lambda_2
 real(dp),intent(out) :: d2edv2_1,d2edv2_2,d2edv2_predict,dedv_predict
 real(dp),intent(out) :: etotal_predict,lambda_predict

!Local variables-------------------------------
!scalars
 real(dp) :: cc,d2edv2_mid,d_lambda,dedv_2bis
 real(dp) :: dedv_mid2,etotal_2bis
 character(len=500) :: message

! *************************************************************************

!DEBUG
!write(std_out,*)' findmin : enter'
!write(std_out,*)' choice,lambda_1,lambda_2=',choice,lambda_1,lambda_2
!ENDDEBUG

 errid = AB7_NO_ERROR
 d_lambda=lambda_1-lambda_2

 if(choice==1) then

!  Use the derivative information to predict lambda
   d2edv2_mid=(dedv_1-dedv_2)/d_lambda
   lambda_predict=lambda_2-dedv_2/d2edv2_mid
   dedv_predict=dedv_2+(lambda_predict-lambda_2)*d2edv2_mid
   d2edv2_1=d2edv2_mid
   d2edv2_2=d2edv2_mid
   d2edv2_predict=d2edv2_mid
!  also use the first energy to predict new energy
   etotal_predict=etotal_1+dedv_1*(lambda_predict-lambda_1)&
&   +0.5_dp*d2edv2_1*(lambda_predict-lambda_1)**2
   etotal_2bis=etotal_1+dedv_1*(lambda_2-lambda_1)&
&   +0.5_dp*d2edv2_1*(lambda_2-lambda_1)**2

   if(d2edv2_mid<0.0_dp)then
     errid = AB7_ERROR_MIXING_INTERNAL
     write(errmess,'(a,es18.10,a)')'The second derivative is negative, equal to ',d2edv2_mid,'.'
     MSG_WARNING(errmess)
   end if

 else if(choice==2) then

!  Use energies and first derivative information
!  etotal = aa + bb * lambda + cc * lambda**2
   dedv_mid2=(etotal_1-etotal_2)/d_lambda
   cc=(dedv_1-dedv_mid2)/d_lambda
   lambda_predict=lambda_1-0.5_dp*dedv_1/cc
   d2edv2_1=2*cc
   d2edv2_2=d2edv2_1
   d2edv2_predict=d2edv2_1
   if(d2edv2_predict<0.0_dp)then
     errid = AB7_ERROR_MIXING_INTERNAL
     write(errmess, '(a,es18.10,a,a,a)' )&
&     'The second derivative is negative, equal to',d2edv2_predict,'.',ch10,&
&     '=> Pivoting                     '
     MSG_WARNING(errmess)
     if(etotal_2 < etotal_1)then
       lambda_predict=lambda_2-0.5_dp*(lambda_1-lambda_2)
     else
       lambda_predict=lambda_1-0.5_dp*(lambda_2-lambda_1)
     end if
   end if
   dedv_predict=dedv_1+(lambda_predict-lambda_1)*d2edv2_1
   dedv_2bis=dedv_1+(lambda_2-lambda_1)*d2edv2_1
   etotal_predict=etotal_1+dedv_1*(lambda_predict-lambda_1)&
&   +0.5_dp*d2edv2_1*(lambda_predict-lambda_1)**2

 end if

 write(message, '(a,es12.4,a,es18.10)' ) &
& ' findmin : lambda_predict ',lambda_predict,' etotal_predict ',etotal_predict
 call wrtout(std_out,message,'COLL')

end subroutine findminscf
!!***

!!****f* ABINIT/dotprodm_v
!! NAME
!! dotprodm_v
!!
!! FUNCTION
!! For two sets of potentials,
!! compute dot product of each pair of two potentials (integral over FFT grid), to obtain
!! a series of square residual-like quantity (so the sum of product of values
!! is NOT divided by the number of FFT points, and NOT multiplied by the primitive cell volume).
!! Take into account the spin components of the potentials (nspden),
!! and sum over them.
!! Need the index of the first pair of potentials to be treated, in each array
!! of potentials, and the number of potentials to be treated.
!! Might be used to compute just one square of norm, in
!! a big array, such as to avoid copying a potential from a big array
!! to a temporary place.
!!
!! INPUTS
!!  cplex=if 1, real space functions on FFT grid are REAL, if 2, COMPLEX
!!  cpldot=if 1, the dot array is real, if 2, the dot array is complex
!!  index1=index of the first potential to be treated in the potarr1 array
!!  index2=index of the first potential to be treated in the potarr2 array
!!  mpicomm=the mpi communicator used for the summation
!!  mpi_summarize=set it to .true. if parallelisation is done over FFT
!!  mult1=number of potentials to be treated in the first set
!!  mult2=number of potentials to be treated in the second set
!!  nfft= (effective) number of FFT grid points (for this processor)
!!  npot1= third dimension of the potarr1 array
!!  npot2= third dimension of the potarr2 array
!!  nspden=number of spin-density components
!!  opt_storage: 0, if potentials are stored as V^up-up, V^dn-dn, Re[V^up-dn], Im[V^up-dn]
!!               1, if potentials are stored as V, B_x, B_y, Bz  (B=magn. field)
!!  potarr1(cplex*nfft,nspden,npot)=first array of real space potentials on FFT grid
!!    (if cplex=2 and cpldot=2, potarr1 is the array that will be complex conjugated)
!!  potarr2(cplex*nfft,nspden,npot)=second array of real space potentials on FFT grid
!!
!! OUTPUT
!!  dot(cpldot,mult1,mult2)= series of values of the dot product
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  Concerning storage when nspden=4:
!!   cplex=1:
!!     opt_storage=0: V are stored as : V^11, V^22, Re[V^12], Im[V^12] (complex, hermitian)
!!     opt_storage=1: V are stored as : V, B_x, B_y, B_z               (real)
!!   cplex=2:
!!     opt_storage=0: V are stored as : V^11, V^22, V^12, i.V^21 (complex)
!!     opt_storage=1: V are stored as : V, B_x, B_y, B_z         (complex)
!!
!! PARENTS
!!      scfopt
!!
!! CHILDREN
!!      timab,xmpi_sum
!!
!! SOURCE

subroutine dotprodm_v(cplex,cpldot,dot,index1,index2,mpicomm,mpi_summarize,&
&   mult1,mult2,nfft,npot1,npot2,nspden,opt_storage,potarr1,potarr2)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cpldot,cplex,index1,index2,mult1,mult2,nfft,npot1,npot2
 integer,intent(in) :: nspden,opt_storage,mpicomm
 logical, intent(in) :: mpi_summarize
!arrays
 real(dp),intent(in) :: potarr1(cplex*nfft,nspden,npot1)
 real(dp),intent(in) :: potarr2(cplex*nfft,nspden,npot2)
 real(dp),intent(out) :: dot(cpldot,mult1,mult2)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,ierr,ifft,ispden
 real(dp) :: ai,ar
!arrays
 real(dp) :: tsec(2)

! *************************************************************************

!Real or complex inputs are coded
 DBG_CHECK(ANY(cplex==(/1,2/)),"Wrong cplex")

!Real or complex outputs are coded
 DBG_CHECK(ANY(cpldot==(/1,2/)),"Wrong cpldot")
 DBG_CHECK(ANY(nspden==(/1,2,4/)),"Wrong nspden")
 DBG_CHECK( npot1-index1-mult1 >= -1,"npot1-index1-mult1")
 DBG_CHECK( npot2-index2-mult2 >= -1,"npot2-index2-mult2")

 if(cplex==1 .or. cpldot==1)then

   do i1=1,mult1
     do i2=1,mult2
       ar=zero
       do ispden=1,min(nspden,2)
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(cplex,i1,i2,index1,index2,ispden,nfft,potarr1,potarr2) REDUCTION(+:ar)
         do ifft=1,cplex*nfft
           ar=ar + potarr1(ifft,ispden,index1+i1-1)*potarr2(ifft,ispden,index2+i2-1)
         end do
       end do
       dot(1,i1,i2)=ar
       if (nspden==4) then
         ar=zero
         do ispden=3,4
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(cplex,i1,i2,index1,index2,ispden,nfft,potarr1,potarr2) REDUCTION(+:ar)
           do ifft=1,cplex*nfft
             ar=ar + potarr1(ifft,ispden,index1+i1-1)*potarr2(ifft,ispden,index2+i2-1)
           end do
         end do
         if (opt_storage==0) then
           if (cplex==1) then
             dot(1,i1,i2)=dot(1,i1,i2)+two*ar
           else
             dot(1,i1,i2)=dot(1,i1,i2)+ar
           end if
         else
           dot(1,i1,i2)=half*(dot(1,i1,i2)+ar)
         end if
       end if
     end do
   end do

 else ! if (cplex==2 .and. cpldot==2)

   do i1=1,mult1
     do i2=1,mult2
       ar=zero ; ai=zero
       do ispden=1,min(nspden,2)
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(cplex,i1,i2,index1,index2,ispden,nfft,potarr1,potarr2) REDUCTION(+:ar,ai)
         do ifft=1,nfft
           ar=ar + potarr1(2*ifft-1,ispden,index1+i1-1)*potarr2(2*ifft-1,ispden,index2+i2-1) &
&           + potarr1(2*ifft  ,ispden,index1+i1-1)*potarr2(2*ifft  ,ispden,index2+i2-1)
           ai=ai + potarr1(2*ifft-1,ispden,index1+i1-1)*potarr2(2*ifft  ,ispden,index2+i2-1) &
&           - potarr1(2*ifft  ,ispden,index1+i1-1)*potarr2(2*ifft-1,ispden,index2+i2-1)
         end do
       end do
       dot(1,i1,i2)=ar ; dot(2,i1,i2)=ai
       if (nspden==4) then
         ar=zero
         do ispden=3,4
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(cplex,i1,i2,index1,index2,ispden,nfft,potarr1,potarr2) REDUCTION(+:ar,ai)
           do ifft=1,nfft
             ar=ar + potarr1(2*ifft-1,ispden,index1+i1-1)*potarr2(2*ifft-1,ispden,index2+i2-1) &
&             + potarr1(2*ifft  ,ispden,index1+i1-1)*potarr2(2*ifft  ,ispden,index2+i2-1)
             ai=ai + potarr1(2*ifft-1,ispden,index1+i1-1)*potarr2(2*ifft  ,ispden,index2+i2-1) &
&             - potarr1(2*ifft  ,ispden,index1+i1-1)*potarr2(2*ifft-1,ispden,index2+i2-1)
           end do
         end do
         if (opt_storage==0) then
           dot(1,i1,i2)=dot(1,i1,i2)+ar
           dot(2,i1,i2)=dot(2,i1,i2)+ai
         else
           dot(1,i1,i2)=half*(dot(1,i1,i2)+ar)
           dot(2,i1,i2)=half*(dot(2,i1,i2)+ai)
         end if
       end if
     end do
   end do
 end if

!XG030513 : MPIWF reduction (addition) on dot is needed here
 if (mpi_summarize) then
   call timab(48,1,tsec)
   call xmpi_sum(dot,mpicomm ,ierr)
   call timab(48,2,tsec)
 end if

 if(cpldot==2 .and. cplex==1)dot(2,:,:)=zero

end subroutine dotprodm_v
!!***

!!****f* ABINIT/dotprodm_vn
!! NAME
!! dotprodm_vn
!!
!! FUNCTION
!! For a set of densities and a set of potentials,
!! compute the dot product (integral over FFT grid) of each pair, to obtain
!! a series of energy-like quantity (so the usual dotproduct is divided
!! by the number of FFT points, and multiplied by the primitive cell volume).
!! Take into account the spin components of the density and potentials (nspden),
!! and sum correctly over them. Note that the storage of densities and
!! potentials is different : for potential, one stores the matrix components,
!! while for the density, one stores the trace, and then, either the
!! spin-polarisation (if nspden=2), or the magnetization vector (if nspden=4).
!! Need the index of the first density/potential pair to be treated, in each array,
!! and the number of pairs to be treated.
!! Might be used to compute just one dot product, in
!! a big array, such as to avoid copying the density and potential from a big array
!! to a temporary place.
!!
!! INPUTS
!!  cplex=if 1, real space functions on FFT grid are REAL, if 2, COMPLEX
!!  cpldot=if 1, the dot array is real, if 2, the dot array is complex (not coded yet for nspden=4)
!!  denarr(cplex*nfft,nspden,nden)=real space density on FFT grid
!!  id=index of the first density to be treated in the denarr array
!!  ip=index of the first potential to be treated in the potarr array
!!  mpicomm=the mpi communicator used for the summation
!!  mpi_summarize=set it to .true. if parallelisation is done over FFT
!!  multd=number of densities to be treated
!!  multp=number of potentials to be treated
!!  nden=third dimension of the denarr array
!!  nfft= (effective) number of FFT grid points (for this processor)
!!  nfftot= total number of FFT grid points
!!  npot=third dimension of the potarr array
!!  nspden=number of spin-density components
!!  potarr(cplex*nfft,nspden,npot)=real space potential on FFT grid
!!                 (will be complex conjugated if cplex=2 and cpldot=2)
!!  ucvol=unit cell volume (Bohr**3)
!!
!! OUTPUT
!!  dot(cpldot,multp,multd)= series of values of the dot product potential/density
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  Concerning storage when nspden=4:
!!   cplex=1:
!!     V are stored as : V^11, V^22, Re[V^12], Im[V^12] (complex, hermitian)
!!     N are stored as : n, m_x, m_y, m_z               (real)
!!   cplex=2:
!!     V are stored as : V^11, V^22, V^12, i.V^21 (complex)
!!     N are stored as : n, m_x, m_y, mZ          (complex)
!!
!! PARENTS
!!      aprxdr
!!
!! CHILDREN
!!      timab,xmpi_sum
!!
!! SOURCE

subroutine dotprodm_vn(cplex,cpldot,denarr,dot,id,ip,mpicomm, mpi_summarize,multd,multp,&
& nden,nfft,nfftot,npot,nspden,potarr,ucvol)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cpldot,cplex,id,ip,multd,multp,nden,nfft,nfftot,npot
 integer,intent(in) :: nspden,mpicomm
 logical, intent(in) :: mpi_summarize
 real(dp),intent(in) :: ucvol
!arrays
 real(dp),intent(in) :: denarr(cplex*nfft,nspden,nden)
 real(dp),intent(in) :: potarr(cplex*nfft,nspden,npot)
 real(dp),intent(out) :: dot(cpldot,multp,multd)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,ierr,ir,jr
 real(dp) :: ai,ar,dim11,dim12,dim21,dim22,dim_dn,dim_up,dre11,dre12,dre21
 real(dp) :: dre22,dre_dn,dre_up,factor,pim11,pim12,pim21,pim22,pim_dn,pim_up
 real(dp) :: pre11,pre12,pre21,pre22,pre_dn,pre_up
!arrays
 real(dp) :: tsec(2)

! *************************************************************************

!Real or complex inputs are coded
 DBG_CHECK(ANY(cplex==(/1,2/)),"Wrong cplex")

!Real or complex outputs are coded
 DBG_CHECK(ANY(cpldot==(/1,2/)),"Wrong cpldot")
 DBG_CHECK(ANY(nspden==(/1,2,4/)),"Wrong nspden")

 DBG_CHECK(id >= 1,'Wrong id')
 DBG_CHECK(ip >= 1,'Wrong id')

 DBG_CHECK(multd >= 1,"wrong multd")
 DBG_CHECK(multp >= 1,"wrong multp")

 DBG_CHECK(nden-id-multd >=-1,'nden-id-multd')
 DBG_CHECK(npot-ip-multp >=-1,'npot-ip-multp')

 if(nspden==1)then

   if(cpldot==1 .or. cplex==1 )then

     do i2=1,multd
       do i1=1,multp
         ar=zero
!$OMP PARALLEL DO PRIVATE(ir) SHARED(id,i1,i2,ip,cplex,nfft,denarr,potarr) REDUCTION(+:ar)
         do ir=1,cplex*nfft
           ar=ar + potarr(ir,1,ip+i1-1)*denarr(ir,1,id+i2-1)
         end do
         dot(1,i1,i2)=ar
       end do ! i1
     end do ! i2

   else  ! cpldot==2 and cplex==2 : one builds the imaginary part, from complex den/pot

     do i2=1,multd
       do i1=1,multp
         ar=zero ; ai=zero
!$OMP PARALLEL DO PRIVATE(ir,jr) SHARED(id,i1,i2,ip,nfft,denarr,potarr) REDUCTION(+:ar,ai)
         do ir=1,nfft
           jr=2*ir
           ar=ar + potarr(jr-1,1,ip+i1-1)*denarr(jr-1,1,id+i2-1) &
&           + potarr(jr  ,1,ip+i1-1)*denarr(jr  ,1,id+i2-1)
           ai=ai + potarr(jr-1,1,ip+i1-1)*denarr(jr  ,1,id+i2-1) &
&           - potarr(jr  ,1,ip+i1-1)*denarr(jr-1,1,id+i2-1)
         end do
         dot(1,i1,i2)=ar ; dot(2,i1,i2)=ai
       end do ! i1
     end do ! i2

   end if

 else if(nspden==2)then

   if(cpldot==1 .or. cplex==1 )then

     do i2=1,multd
       do i1=1,multp
         ar=zero
!$OMP PARALLEL DO PRIVATE(ir) SHARED(id,i1,i2,ip,cplex,nfft,denarr,potarr) REDUCTION(+:ar)
         do ir=1,cplex*nfft
           ar=ar + potarr(ir,1,ip+i1-1)* denarr(ir,2,id+i2-1)               &       ! This is the spin up contribution
&          + potarr(ir,2,ip+i1-1)*(denarr(ir,1,id+i2-1)-denarr(ir,2,id+i2-1)) ! This is the spin down contribution
         end do
         dot(1,i1,i2)=ar
       end do ! i1
     end do ! i2

   else ! cpldot==2 and cplex==2 : one builds the imaginary part, from complex den/pot

     do i2=1,multd
       do i1=1,multp
         ar=zero ; ai=zero
!$OMP PARALLEL DO PRIVATE(ir,jr,dre_up,dim_up,dre_dn,dim_dn,pre_up,pim_up,pre_dn,pim_dn) &
!$OMP&SHARED(id,i1,i2,ip,nfft,denarr,potarr) REDUCTION(+:ar,ai)
         do ir=1,nfft
           jr=2*ir

           dre_up=denarr(jr-1,2,id+i2-1)
           dim_up=denarr(jr  ,2,id+i2-1)
           dre_dn=denarr(jr-1,1,id+i2-1)-dre_up
           dim_dn=denarr(jr  ,1,id+i2-1)-dim_up

           pre_up=potarr(jr-1,1,ip+i1-1)
           pim_up=potarr(jr  ,1,ip+i1-1)
           pre_dn=potarr(jr-1,2,ip+i1-1)
           pim_dn=potarr(jr  ,2,ip+i1-1)

           ar=ar + pre_up * dre_up &
&           + pim_up * dim_up &
&           + pre_dn * dre_dn &
&           + pim_dn * dim_dn
           ai=ai + pre_up * dim_up &
&           - pim_up * dre_up &
&           + pre_dn * dim_dn &
&           - pim_dn * dre_dn

         end do
         dot(1,i1,i2)=ar ; dot(2,i1,i2)=ai
       end do ! i1
     end do ! i2

   end if

 else if(nspden==4)then
!  \rho{\alpha,\beta} V^{\alpha,\beta} =
!  rho*(V^{11}+V^{22})/2$
!  + m_x Re(V^{12})- m_y Im{V^{12}}+ m_z(V^{11}-V^{22})/2
   if (cplex==1) then
     do i2=1,multd
       do i1=1,multp
         ar=zero
!$OMP PARALLEL DO PRIVATE(ir) SHARED(id,i1,i2,ip,cplex,nfft,denarr,potarr) REDUCTION(+:ar)
         do ir=1,cplex*nfft
           ar=ar+(potarr(ir,1,ip+i1-1)+potarr(ir,2,ip+i1-1))*half*denarr(ir,1,id+i2-1)& ! This is the density contrib
&          + potarr(ir,3,ip+i1-1)                                *denarr(ir,2,id+i2-1)& ! This is the m_x contrib
&          - potarr(ir,4,ip+i1-1)                                *denarr(ir,3,id+i2-1)& ! This is the m_y contrib
&          +(potarr(ir,1,ip+i1-1)-potarr(ir,2,ip+i1-1))*half*denarr(ir,4,id+i2-1)       ! This is the m_z contrib
         end do
         dot(1,i1,i2)=ar
       end do ! i1
     end do ! i2
   else ! cplex=2
!    Note concerning storage when cplex=2:
!    V are stored as : v^11, v^22, V^12, i.V^21 (each are complex)
!    N are stored as : n, m_x, m_y, mZ          (each are complex)
     if (cpldot==1) then
       do i2=1,multd
         do i1=1,multp
           ar=zero ; ai=zero
!$OMP PARALLEL DO PRIVATE(ir,jr,dre11,dim11,dre22,dim22,dre12,dim12,pre11,pim11,pre22,pim22,pre12,pim12) &
!$OMP&SHARED(id,i1,i2,ip,nfft,denarr,potarr) REDUCTION(+:ar)
           do ir=1,nfft
             jr=2*ir
             dre11=half*(denarr(jr-1,1,id+i2)+denarr(jr-1,4,id+i2))
             dim11=half*(denarr(jr  ,1,id+i2)+denarr(jr-1,4,id+i2))
             dre22=half*(denarr(jr-1,1,id+i2)-denarr(jr-1,4,id+i2))
             dim22=half*(denarr(jr  ,1,id+i2)-denarr(jr-1,4,id+i2))
             dre12=half*(denarr(jr-1,2,id+i2)+denarr(jr  ,3,id+i2))
             dim12=half*(denarr(jr  ,2,id+i2)-denarr(jr-1,3,id+i2))
             dre21=half*(denarr(jr-1,2,id+i2)-denarr(jr  ,3,id+i2))
             dim21=half*(denarr(jr  ,2,id+i2)+denarr(jr-1,3,id+i2))
             pre11= potarr(jr-1,1,ip+i1)
             pim11= potarr(jr  ,1,ip+i1)
             pre22= potarr(jr-1,2,ip+i1)
             pim22= potarr(jr  ,2,ip+i1)
             pre12= potarr(jr-1,3,ip+i1)
             pim12= potarr(jr  ,3,ip+i1)
             pre21= potarr(jr  ,4,ip+i1)
             pim21=-potarr(jr-1,4,ip+i1)
             ar=ar + pre11 * dre11 &
&             + pim11 * dim11 &
&             + pre22 * dre22 &
&             + pim22 * dim22 &
&             + pre12 * dre12 &
&             + pim12 * dim12 &
&             + pre21 * dre21 &
&             + pim21 * dim21
           end do
           dot(1,i1,i2)=ar
         end do ! i1
       end do ! i2
     else !cpldot=2
       do i2=1,multd
         do i1=1,multp
           ar=zero ; ai=zero
!$OMP PARALLEL DO PRIVATE(ir,jr,dre11,dim11,dre22,dim22,dre12,dim12,pre11,pim11,pre12,pim12,pre22,pim22) &
!$OMP&SHARED(id,i1,i2,ip,nfft,denarr,potarr) REDUCTION(+:ar,ai)
           do ir=1,nfft
             jr=2*ir
             dre11=half*(denarr(jr-1,1,id+i2)+denarr(jr-1,4,id+i2))
             dim11=half*(denarr(jr  ,1,id+i2)+denarr(jr-1,4,id+i2))
             dre22=half*(denarr(jr-1,1,id+i2)-denarr(jr-1,4,id+i2))
             dim22=half*(denarr(jr  ,1,id+i2)-denarr(jr-1,4,id+i2))
             dre12=half*(denarr(jr-1,2,id+i2)+denarr(jr  ,3,id+i2))
             dim12=half*(denarr(jr  ,2,id+i2)-denarr(jr-1,3,id+i2))
             dre21=half*(denarr(jr-1,2,id+i2)-denarr(jr  ,3,id+i2))
             dim21=half*(denarr(jr  ,2,id+i2)+denarr(jr-1,3,id+i2))
             pre11= potarr(jr-1,1,ip+i1)
             pim11= potarr(jr  ,1,ip+i1)
             pre22= potarr(jr-1,2,ip+i1)
             pim22= potarr(jr  ,2,ip+i1)
             pre12= potarr(jr-1,3,ip+i1)
             pim12= potarr(jr  ,3,ip+i1)
             pre21= potarr(jr  ,4,ip+i1)
             pim21=-potarr(jr-1,4,ip+i1)
             ar=ar + pre11 * dre11 &
&             + pim11 * dim11 &
&             + pre22 * dre22 &
&             + pim22 * dim22 &
&             + pre12 * dre12 &
&             + pim12 * dim12 &
&             + pre21 * dre21 &
&             + pim21 * dim21
             ai=ai + pre11 * dim11 &
&             - pim11 * dre11 &
&             + pre22 * dim22 &
&             - pim22 * dre22 &
&             + pre12 * dim12 &
&             - pim12 * dre12 &
&             + pre21 * dim21 &
&             - pim21 * dre21
           end do
           dot(1,i1,i2)=ar
           dot(2,i1,i2)=ai
         end do ! i1
       end do ! i2
     end if ! cpldot
   end if ! cplex
 end if ! nspden

 factor=ucvol/dble(nfftot)
 dot(:,:,:)=factor*dot(:,:,:)

!XG030513 : MPIWF reduction (addition) on dot is needed here
 if (mpi_summarize) then
   call timab(48,1,tsec)
   call xmpi_sum(dot,mpicomm ,ierr)
   call timab(48,2,tsec)
 end if

 if(cpldot==2 .and. cplex==1)dot(2,:,:)=zero

end subroutine dotprodm_vn
!!***

!!****f* ABINIT/sqnormm_v
!! NAME
!! sqnormm_v
!!
!! FUNCTION
!! For a series of potentials,
!! compute square of the norm (integral over FFT grid), to obtain
!! a square residual-like quantity (so the sum of product of values
!! is NOT divided by the number of FFT points, and NOT multiplied by the primitive cell volume).
!! Take into account the spin components of the density and potentials (nspden), and sum over them.
!! Need the index of the first potential to be treated, in the provided array
!! of potentials, and the number of potentials to be treated.
!! Might be used to compute just one square of norm, in a big array, such as to avoid
!! copying a potential from a big array to a temporary place.
!!
!! INPUTS
!!  cplex=if 1, real space function on FFT grid is REAL, if 2, COMPLEX
!!  index=index of the first potential to be treated
!!  mpicomm=the mpi communicator used for the summation
!!  mpi_summarize=set it to .true. if parallelisation is done over FFT
!!  mult=number of potentials to be treated
!!  nfft= (effective) number of FFT grid points (for this processor)
!!  npot= third dimension of the potarr array
!!  nspden=number of spin-density components
!!  opt_storage: 0, if potential is stored as V^up-up, V^dn-dn, Re[V^up-dn], Im[V^up-dn]
!!               1, if potential is stored as V, B_x, B_y, Bz  (B=magn. field)
!!  potarr(cplex*nfft,nspden,npot)=array of real space potentials on FFT grid
!!
!! OUTPUT
!!  norm2(mult)= value of the square of the norm of the different potentials
!!
!! PARENTS
!!      scfcge,scfopt
!!
!! CHILDREN
!!      timab,xmpi_sum
!!
!! SOURCE

subroutine sqnormm_v(cplex,index,mpicomm, mpi_summarize,mult,nfft,norm2,npot,nspden,opt_storage,potarr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,index,mult,nfft,npot,nspden,opt_storage,mpicomm
 logical, intent(in) :: mpi_summarize
!arrays
 real(dp),intent(in) :: potarr(cplex*nfft,nspden,npot)
 real(dp),intent(out) :: norm2(mult)

!Local variables-------------------------------
!scalars
 integer :: ierr,ifft,ii,ispden
 real(dp) :: ar
!arrays
 real(dp) :: tsec(2)

! *************************************************************************

!Real or complex inputs are coded
 DBG_CHECK(ANY(cplex==(/1,2/)),"Wrong cplex")
 DBG_CHECK(ANY(nspden==(/1,2,4/)),"Wrong nspden")

 DBG_CHECK(index>=1,"wrong index")
 DBG_CHECK(mult>=1,"wrong mult")
 DBG_CHECK(npot>=1,"wrong npot")

 DBG_CHECK(npot-index-mult>=-1,'npot-index-mult')

 do ii=1,mult
   ar=zero
   do ispden=1,min(nspden,2)
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(cplex,ii,index,ispden,nfft,potarr) REDUCTION(+:ar)
     do ifft=1,cplex*nfft
       ar=ar + potarr(ifft,ispden,index+ii-1)**2
     end do
   end do
   norm2(ii)=ar
   if (nspden==4) then
     ar=zero
     do ispden=3,4
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(cplex,ii,index,ispden,nfft,potarr) REDUCTION(+:ar)
       do ifft=1,cplex*nfft
         ar=ar + potarr(ifft,ispden,index+ii-1)**2
       end do
     end do
     if (opt_storage==0) then
       if (cplex==1) then
         norm2(ii)=norm2(ii)+two*ar
       else
         norm2(ii)=norm2(ii)+ar
       end if
     else
       norm2(ii)=half*(norm2(ii)+ar)
     end if
   end if
 end do

!XG030513 : MPIWF reduction (addition) on norm2 is needed here
 if (mpi_summarize) then
   call timab(48,1,tsec)
   call xmpi_sum(norm2,mpicomm ,ierr)
   call timab(48,2,tsec)
 end if

end subroutine sqnormm_v
!!***

!!****f* ABINIT/aprxdr
!! NAME
!! aprxdr
!!
!! FUNCTION
!! Compute the approximative derivatives of the energy at different
!! points along the line search, thanks to a finite-difference formula.
!! This formula is the projection along the line search of the
!! Eq.(11) in PRB54, 4383 (1996) [[cite:Gonze1996]].
!!
!! INPUTS
!! cplex: if 1, real space functions on FFT grid are REAL, if 2, COMPLEX
!! choice= if==3, compute dedv_new, dedv_old, and dedv_mix,
!! if/=3, compute only dedv_new and dedv_old.
!! i_vresid and i_rhor, see the next lines.
!! f_fftgr(nfft,nspden,n_fftgr)=different functions defined on the fft grid :
!! The last residual potential is in f_fftgr(:,:,i_vresid(1)).
!! The old  residual potential is in f_fftgr(:,:,i_vresid(2)).
!! The previous old residual potential is in f_fftgr(:,:,i_vresid(3)).
!! (needed only when choice==3)
!! The old  density is in f_fftgr(:,:,i_rhor2).
!! f_atm(3,natom,n_fftgr)=different functions defined for each atom :
!! The last HF force is in f_atm(:,:,i_vresid(1)).
!! The old  HF force is in f_fftgr(:,:,i_vresid(2)).
!! The previous old HF force is in f_fftgr(:,:,i_vresid(3)).
!! (needed only when choice==3)
!! The old atomic positions are in f_atm(:,:,i_rhor2)
!! moved_atm_inside: if==1, the atoms are allowed to move.
!! mpicomm=the mpi communicator used for the summation
!! mpi_summarize=set it to .true. if parallelisation is done over FFT
!! natom=number of atoms in unit cell
!! nfft=(effective) number of FFT grid points (for this processor)
!! nfftot=total number of FFT grid points
!! nspden=number of spin-density components
!! rhor(nfft,nspden)=actual density
!! xred(3,natom)=reduced atomic coordinates
!!
!! OUTPUT
!! dedv_mix=approximate derivative from previous old residual
!! dedv_new=approximate derivative from new residual
!! dedv_old=approximate derivative from old residual (output only when choice==3)
!!
!! NOTES
!! Should be OpenMP parallelized
!!
!! PARENTS
!!      scfcge
!!
!! CHILDREN
!!      dotprodm_vn
!!
!! SOURCE

subroutine aprxdr(cplex,choice,dedv_mix,dedv_new,dedv_old,&
&  f_atm,f_fftgr,i_rhor2,i_vresid,moved_atm_inside,&
&  mpicomm,mpi_summarize,natom,nfft,nfftot,nspden,n_fftgr,rhor,ucvol,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,cplex,i_rhor2,moved_atm_inside,n_fftgr,natom,nfft
 integer,intent(in) :: mpicomm,nfftot,nspden
 logical, intent(in) :: mpi_summarize
 real(dp),intent(in) :: ucvol
 real(dp),intent(out) :: dedv_mix,dedv_new,dedv_old
!arrays
 integer,intent(in) :: i_vresid(3)
 real(dp),intent(in) :: f_atm(3,natom,n_fftgr)
 real(dp),intent(in) :: f_fftgr(cplex*nfft,nspden,n_fftgr)
 real(dp),intent(in) :: rhor(cplex*nfft,nspden),xred(3,natom)

!Local variables-------------------------------
!scalars
 integer :: iatom,idir
!arrays
 real(dp) :: dedv_temp(1)
 real(dp),allocatable :: ddens(:,:,:)

! *************************************************************************

 ABI_ALLOCATE(ddens,(cplex*nfft,nspden,1))

!Compute approximative derivative of the energy
!with respect to change of potential

 ddens(:,:,1)=rhor(:,:)-f_fftgr(:,:,i_rhor2)

!call dotprod_vn(cplex,1,ddens,dedv_old,nfft,nfftot,nspden,1,vresid,ucvol)
!Dot product ddens(:,:,1) f_fftgr(:,:,i_vresid(2))
 call dotprodm_vn(cplex,1,ddens,dedv_temp,1,i_vresid(2),mpicomm,mpi_summarize,1,1,1,&
& nfft,nfftot,n_fftgr,nspden,f_fftgr,ucvol)
 dedv_old = dedv_temp(1)

!Dot product ddens(:,:,1) f_fftgr(:,:,i_vresid(1))
 call dotprodm_vn(cplex,1,ddens,dedv_temp,1,i_vresid(1),mpicomm,mpi_summarize,1,1,1,&
& nfft,nfftot,n_fftgr,nspden,f_fftgr,ucvol)
 dedv_new= dedv_temp(1)

 if(choice==3)then
!  Dot product ddens(:,:,1) f_fftgr(:,:,i_vresid(3))
   call dotprodm_vn(cplex,1,ddens,dedv_temp,1,i_vresid(3),mpicomm,mpi_summarize,1,1,1,&
&   nfft,nfftot,n_fftgr,nspden,f_fftgr,ucvol)
   dedv_mix = dedv_temp(1)
 end if

 ABI_DEALLOCATE(ddens)

!-------------------------------------------------------

!Now, take care of eventual atomic displacements

 if(moved_atm_inside==1)then
   do idir=1,3
     do iatom=1,natom
       dedv_new=dedv_new+&
&       f_atm(idir,iatom,i_vresid(1))*(xred(idir,iatom)-f_atm(idir,iatom,i_rhor2))
       dedv_old=dedv_old+&
&       f_atm(idir,iatom,i_vresid(2))*(xred(idir,iatom)-f_atm(idir,iatom,i_rhor2))
       if(choice==3) dedv_mix=dedv_mix+&
&       f_atm(idir,iatom,i_vresid(3))*(xred(idir,iatom)-f_atm(idir,iatom,i_rhor2))
     end do
   end do
 end if

end subroutine aprxdr
!!***

end module m_ab7_mixing
!!***
