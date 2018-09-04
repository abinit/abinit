!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ab7_invars_f90
!! NAME
!! m_ab7_invars_f90
!!
!! FUNCTION
!! driver for the parser
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_ab7_invars

  use defs_basis
  use defs_datatypes
  use defs_abitypes
  use m_errors
  use m_abicore
  use m_xmpi

  use m_time,     only : timab
  use m_fstrings, only : inupper
  use m_parser,   only : intagm, importxyz, parsefile
  use m_dtset,    only : dtset_free, macroin, macroin2
  use m_dtfil,    only : status
  use m_pspheads, only : inpspheads, pspheads_comm
  use m_invars1,  only : invars1m, invars0, indefo
  use m_invars2,  only : invars2, invars2m

  implicit none

  private

  ! We store here a list of dtset arrays to be able to
  ! parse several ABINIT files without freeing it.
  ! The simplest portable way to do it, is to create
  ! a list of dtsets arrays and to use the list index
  ! as an identifier that can be given to the other languages.
  type, private :: dtsets_list
     integer                     :: id
     type(dtsets_list),  pointer :: next => null()
     type(dtsets_list),  pointer :: prev => null()
     type(dataset_type), pointer :: dtsets(:)
     type(pspheader_type), pointer :: pspheads(:)=>null()   !vz_z
     integer :: mxga_n_rules, mxgw_nqlwl, mxlpawu, &
         & mxmband, mxmband_upper, &
         & mxnatom, mxnatpawu, mxnatsph, mxnatsph_extra, mxnatvshift, &
         & mxnbandhf, mxnconeq, mxnfreqsp, &
         & mxnimage, mxnimfrqs, mxnkpt, mxnkptgw, mxnkpthf, &
         & mxnnos, mxnqptdm, mxnspinor, mxnsppol, mxnsym, &
         & mxntypat, mxnzchempot, &
         & mxn_efmas_dirs, mxn_projection_frequencies
     integer :: istatr, istatshft, dmatpuflag, papiopt, timopt
  end type dtsets_list

  type(dtsets_list), save, pointer :: my_dtsets => null()
  integer, save :: nb_dtsets = 0

  ! These flags should be .true. inside ABINIT.
  ! Use ab7_invars_set_flags() to change them.
  logical, save,private :: call_status = .false.
  character(len = fnlen), save,private :: opt_status_file
  logical, save, private :: call_timab = .false.
  real(dp), save, pointer, private :: opt_timab_tsec(:)

  ! These pointers are used only inside ABINIT and refers to their
  ! equivalent. The C binding don't support them.
  ! Use ab7_invars_get_abinit_vars() to get them.

  logical, private, parameter :: AB_DBG = .false.

#include "ab7_invars_f90.inc"

  ! The following group is used for Fortran bindings only,
  ! and specifically its usage inside ABINIT. They have no C or Python equivalent.
  public :: ab7_invars_set_flags
!  public :: ab7_invars_set_mpi
  public :: ab7_invars_get_abinit_vars
  public :: ab7_invars_load

  ! The following routines are the main creation routines, having also an equivalent in C or Python.
  public :: ab7_invars_new_from_file
  public :: ab7_invars_new_from_string
  public :: ab7_invars_free

  ! The following routines are the main getter functions, also available in C or Python.
  public :: ab7_invars_get_ndtset
  public :: ab7_invars_get_integer
  public :: ab7_invars_get_real
  public :: ab7_invars_get_shape
  public :: ab7_invars_get_integer_array
  public :: ab7_invars_get_real_array

contains
!!***

!----------------------------------------------------------------------

!!****f* m_ab7_invars/ab7_invars_set_flags
!! NAME
!!   ab7_invars_set_flags
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      abinit,multibinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine ab7_invars_set_flags(status, timab, status_file, timab_tsec)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ab7_invars_set_flags'
!End of the abilint section

 logical, intent(in) :: status, timab
 character(len = fnlen), intent(in), optional :: status_file
 real(dp), intent(in), target, optional :: timab_tsec(:)

 call_status = status
 if (present(status_file)) then
    write(opt_status_file, "(A)") status_file
 else
    write(opt_status_file, "(A)") "status"
 end if
 call_timab  = timab
 if (present(timab_tsec)) then
    opt_timab_tsec => timab_tsec
 else
    ABI_ALLOCATE(opt_timab_tsec,(2))
 end if

end subroutine ab7_invars_set_flags
!!***

!----------------------------------------------------------------------

!!****f* m_ab7_invars/ab7_invars_get_abinit_vars
!! NAME
!!   ab7_invars_get_abinit_vars
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!
!! SOURCE

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ab7_invars_get_abinit_vars'
!End of the abilint section

 subroutine ab7_invars_get_abinit_vars(dtsetsId, dtsets, pspheads, mxvals, papiopt, timopt, dmatpuflag)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ab7_invars_get_abinit_vars'
!End of the abilint section

 integer, intent(in) :: dtsetsId
 type(dataset_type), pointer :: dtsets(:)
 type(pspheader_type), pointer :: pspheads(:)
 type(ab_dimensions), intent(out) :: mxvals

 integer, intent(out) :: papiopt, timopt, dmatpuflag

 type(dtsets_list), pointer :: token

 call get_token(token, dtsetsId)

 if (associated(token)) then
    dtsets        => token%dtsets
    pspheads      => token%pspheads

    mxvals%ga_n_rules  = token%mxga_n_rules
    mxvals%gw_nqlwl    = token%mxgw_nqlwl
    mxvals%lpawu       = token%mxlpawu
    mxvals%mband       = token%mxmband
    mxvals%mband_upper = token%mxmband_upper
    mxvals%natom       = token%mxnatom
    mxvals%natpawu     = token%mxnatpawu
    mxvals%natsph      = token%mxnatsph
    mxvals%natsph_extra= token%mxnatsph_extra
    mxvals%natvshift   = token%mxnatvshift
    mxvals%nbandhf     = token%mxnbandhf
    mxvals%nconeq      = token%mxnconeq
    mxvals%n_efmas_dirs= token%mxn_efmas_dirs
    mxvals%nimage      = token%mxnimage
    mxvals%nimfrqs     = token%mxnimfrqs
    mxvals%nfreqsp     = token%mxnfreqsp
    mxvals%n_projection_frequencies = token%mxn_projection_frequencies
    mxvals%nkpt        = token%mxnkpt
    mxvals%nkpthf      = token%mxnkpthf
    mxvals%nkptgw      = token%mxnkptgw
    mxvals%nnos        = token%mxnnos
    mxvals%nqptdm      = token%mxnqptdm
    mxvals%nspinor     = token%mxnspinor
    mxvals%nsppol      = token%mxnsppol
    mxvals%nsym        = token%mxnsym
    mxvals%ntypat      = token%mxntypat
    mxvals%nzchempot   = token%mxnzchempot

    mxvals%nberry      = 20   ! This is presently a fixed value. Should be changed.

    papiopt       = token%papiopt
    timopt        = token%timopt
    dmatpuflag    = token%dmatpuflag
 else
    nullify(dtsets)
    nullify(pspheads)
 end if

end subroutine ab7_invars_get_abinit_vars
!!***

!----------------------------------------------------------------------

!!****f* m_ab7_invars/new_token
!! NAME
!!   new_token
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_ab7_invars_f90
!!
!! CHILDREN
!!
!! SOURCE

  subroutine new_token(token)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'new_token'
!End of the abilint section

type(dtsets_list), pointer :: token

 ! We allocate a new list token and prepend it.
 if (AB_DBG) write(std_err,*) "AB module: create a new token."
 nb_dtsets = nb_dtsets + 1

 ABI_DATATYPE_ALLOCATE(token,)
 token%id = nb_dtsets
 nullify(token%dtsets)
 nullify(token%pspheads)
 token%next => my_dtsets
 nullify(token%prev)

 my_dtsets => token
 if (AB_DBG) write(std_err,*) "AB module: creation OK with id ", token%id

end subroutine new_token
!!***

!----------------------------------------------------------------------

!!****f* m_ab7_invars/free_token
!! NAME
!!   free_token
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_ab7_invars_f90
!!
!! CHILDREN
!!
!! SOURCE

  subroutine free_token(token)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'free_token'
!End of the abilint section

type(dtsets_list), pointer :: token

integer :: idtset

 if (.not. associated(token)) then
    write (std_out,*) 'in m_ab7_invars:free_token : token not associated. Nothing doing.'
    return
 end if

 ! We free a token list.
 if (AB_DBG) write(std_err,*) "AB module: free request on dataset array ", token%id
 if (associated(token%dtsets)) then
    if (AB_DBG) write(std_err,*) " | ", size(token%dtsets), "dtsets found."
    do idtset = 0, size(token%dtsets) - 1, 1
       if (AB_DBG) write(std_err,*) " | free dtset ", idtset
       call dtset_free(token%dtsets(idtset))
       if (AB_DBG) write(std_err,*) " | free OK"
    end do
    ABI_DATATYPE_DEALLOCATE(token%dtsets)
    nullify(token%dtsets)
    ABI_DATATYPE_DEALLOCATE(token%pspheads)
    nullify(token%pspheads)
    if (AB_DBG) write(std_err,*) " | general free OK"

    ! We remove token from the list.
    if (associated(token%prev)) then
       token%prev%next => token%next
    else
       my_dtsets => token%next
    end if
    if (associated(token%next)) then
       token%next%prev => token%prev
    end if
    ABI_DATATYPE_DEALLOCATE(token)
    if (AB_DBG) write(std_err,*) " | token free OK"

 end if
 if (AB_DBG) write(std_err,*) "AB module: free done"

end subroutine free_token
!!***

!----------------------------------------------------------------------

!!****f* m_ab7_invars/get_token
!! NAME
!!  get_token
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_ab7_invars_f90
!!
!! CHILDREN
!!
!! SOURCE

 subroutine get_token(token, id)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_token'
!End of the abilint section

 type(dtsets_list), pointer :: token
 integer, intent(in) :: id

 type(dtsets_list), pointer :: tmpLst

 if (AB_DBG) write(std_err,*) "AB module: request list element ", id
 nullify(token)
 ! List element are prepended so element id is at (nb - id) position.
 tmpLst => my_dtsets
 do
    if (.not. associated(tmpLst)) then
       exit
    end if
    if (tmpLst%id == id .and. associated(tmpLst%dtsets)) then
       token => tmpLst
       return
    end if
    tmpLst => tmpLst%next
 end do

end subroutine get_token
!!***

!----------------------------------------------------------------------

!!****f* m_ab7_invars/ab7_invars_new_from_string
!! NAME
!!
!! FUNCTION
!!   Create a datasets object from a string containing the Abinit input.
!!
!! INPUTS
!!  instr=Input string
!!  len=len of the string
!!
!! OUTPUT
!!   dtsetsId=Datasets identifier
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

  subroutine ab7_invars_new_from_string(dtsetsId, instr, len)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ab7_invars_new_from_string'
!End of the abilint section

 integer, intent(out) :: dtsetsId
 integer, intent(in) :: len
 character(len = len), intent(in) :: instr

 character(len = strlen) :: string
 integer :: lenstr, ndtset
 integer :: marr, tread
 character(len = 30) :: token
 integer :: intarr(1)
 real(dp) :: dprarr(1)
 character(len=500) :: message

 dtsetsId = 0

 if (len > strlen) then
    return
 end if

 write(string,*) instr

 !To make case-insensitive, map characters of string to upper case:
 call inupper(string(1:len))

 !Might import data from xyz file(s) into string
 !Need string_raw to deal properly with xyz filenames
 lenstr = len
 call importxyz(lenstr, instr, string, len)

 !6) Take ndtset from the input string, then allocate
 !the arrays whose dimensions depends only on ndtset and msym.

 ndtset=0 ; marr=1
 token = 'ndtset'
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) ndtset=intarr(1)
 !Check that ndtset is not negative
 if (ndtset<0 .or. ndtset>99) then
    write(message, '(a,i12,a,a,a,a)' )&
&     'Input ndtset must be non-negative and < 100, but was ',ndtset,ch10,&
&     'This is not allowed.  ',ch10,&
&     'Action : modify ndtset in the input file.'
    MSG_ERROR(message)
 end if

 call ab7_invars_load(dtsetsId, string, lenstr, ndtset, .false., .false.)

end subroutine ab7_invars_new_from_string
!!***

!----------------------------------------------------------------------

!!****f* m_ab7_invars/ab7_invars_new_from_file
!! NAME
!!  ab7_invars_new_from_file
!!
!! FUNCTION
!!   Create a datasets object from a file with the Abinit input.
!!
!! INPUTS
!!  filename=String with the path to the file.
!!  n=len of filename.
!!  pspfiles=List of pseudopotential filenames.
!!  nsps=Number of pseudos
!!  comm=MPI communicator
!!
!! OUTPUT
!!   dtsetsId=Datasets identifier
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ab7_invars_new_from_file(dtsetsId, filename, n, pspfiles, npsp, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ab7_invars_new_from_file'
!End of the abilint section

 integer, intent(out) :: dtsetsId
 integer, intent(in) :: n, npsp
 integer,optional,intent(in) :: comm
 character(len=n), intent(in) :: filename
 character(len=fnlen), intent(in) :: pspfiles(npsp)

 character(len = strlen) :: string
 integer :: lenstr, ndtset,my_comm

 my_comm = xmpi_world; if (present(comm)) my_comm = comm

 dtsetsId = 0

 if (AB_DBG) write(std_err,*) "AB module: read '", trim(filename), "' to string."
 call parsefile(filename, lenstr, ndtset, string, my_comm)
 if (AB_DBG) write(std_err,*) "AB module: read OK, string length ", lenstr

 if (npsp == 0) then
    call ab7_invars_load(dtsetsId, string, lenstr, ndtset, .false., .false.)
 else
    call ab7_invars_load(dtsetsId, string, lenstr, ndtset, .true., .false., pspfiles)
 end if

end subroutine ab7_invars_new_from_file
!!***

!----------------------------------------------------------------------

!!****f* m_ab7_invars/ab7_invars_load
!! NAME
!!  ab7_invars_load
!!
!! FUNCTION
!!   Parse the input string and create a datasets object
!!
!! INPUTS
!!   string=String with the input.
!!   lenstr=Length of string.
!!   ndtset=Number of datasets.
!!   with_psp=True if pseudos must be read.
!!   with_mem
!!   [pspfilnam]=List of pseudopotential files
!!
!! OUTPUT
!!   dtsetsId=Datasets identifier
!!
!! PARENTS
!!      abinit,m_ab7_invars_f90
!!
!! CHILDREN
!!
!! SOURCE

 subroutine ab7_invars_load(dtsetsId, string, lenstr, ndtset, with_psp, with_mem, pspfilnam)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ab7_invars_load'
!End of the abilint section

 integer, intent(out) :: dtsetsId
 character(len = strlen), intent(inout) :: string
 integer, intent(in) :: lenstr, ndtset
 logical, intent(in) :: with_psp, with_mem
 character(len = fnlen), intent(in), optional :: pspfilnam(:)

 type(dtsets_list), pointer :: token
 integer, parameter :: level=3
 integer :: jdtset,nprocs,ipsp
 integer :: me,ndtset_alloc
 integer :: npsp, ii, idtset, msym, usepaw
 integer,allocatable :: mband_upper_(:)
 real(dp),allocatable :: zionpsp(:)
 real(dp) :: ecut_tmp(3,2,10)
 character(len = fnlen), allocatable :: pspfilnam_(:)
 !character(len=500) :: message

 ! We allocate a new list token and prepend it.
 if (AB_DBG) write(std_err,*) "AB module: allocate a new object."
 call new_token(token)
 dtsetsId = token%id

 ndtset_alloc=ndtset ; if(ndtset==0)ndtset_alloc=1
 ABI_DATATYPE_ALLOCATE(token%dtsets,(0:ndtset_alloc))
 if (AB_DBG) write(std_err,*) "AB module: allocation OK at ", dtsetsId

 if (AB_DBG) write(std_err,*) "AB module: call invars0()."
 token%timopt = 1
 if(xmpi_paral==1) token%timopt = 0

 !7) Continue to analyze the input string, get upper dimensions,
 !and allocate the remaining arrays.
 call invars0(token%dtsets,token%istatr,token%istatshft,lenstr,&
      & msym,token%mxnatom,token%mxnimage,token%mxntypat,ndtset,ndtset_alloc,npsp,&
      & token%papiopt, token%timopt, string)
 token%dtsets(:)%timopt=token%timopt
 token%dtsets(0)%timopt = 1
 if(xmpi_paral==1) token%dtsets(0)%timopt = 0

 !Be careful : at these fourth and fifth calls of status, istatr and istatshft taken
 !from the input variables will be saved definitively.
 if (call_status) then
   call status(0,opt_status_file,token%istatr,level,'init istatr   ')
   call status(0,opt_status_file,token%istatshft,level,'init istatshft')
 else
!  Fake call to disable status file
   call status(0,' ',0,3,' ')
 end if

 if (call_timab) then
   call timab(41,2,opt_timab_tsec)
   call timab(token%timopt,5,opt_timab_tsec)
 end if

 !8) Finish to read the "file" file completely, as npsp is known,
 !and also initialize pspheads, that contains the important information
 !from the pseudopotential headers, as well as the psp filename

 if (call_timab) then
   call timab(42,1,opt_timab_tsec)
 end if

 if (call_status) then
   call status(0,opt_status_file,99,level,'call iofn2    ')
 end if
 usepaw=0
 ABI_DATATYPE_ALLOCATE(token%pspheads,(npsp))
 if (npsp>10) then
    MSG_BUG('ecut_tmp is not well defined.')
 end if
 ecut_tmp=-one

 token%pspheads(:)%usewvl=token%dtsets(1)%usewvl
 if (with_psp) then
    if (AB_DBG) write(std_err,*) "AB module: call iofn2()."
    me=xmpi_comm_rank(xmpi_world)
    if (me == 0) then
       if (.not. present(pspfilnam)) then
          ABI_ALLOCATE(pspfilnam_,(npsp))
          call iofn2(npsp, pspfilnam_)
          call inpspheads(pspfilnam_,npsp,token%pspheads,ecut_tmp)
!      write(std_out,*)' ab7_invars_f90 : token%pspheads(1)%nproj(0:3)=',token%pspheads(1)%nproj(0:3)
          ABI_DEALLOCATE(pspfilnam_)
       else
          call inpspheads(pspfilnam,npsp,token%pspheads,ecut_tmp)
       end if
       if(minval(abs(token%pspheads(1:npsp)%pspcod-7))==0) usepaw=1
       if(minval(abs(token%pspheads(1:npsp)%pspcod-17))==0) usepaw=1
    end if
    !Communicate pspheads to all processors
    call pspheads_comm(npsp,token%pspheads,usepaw)
 else
    ! No psp files are given, we put default values into pspheads.
    token%pspheads(:)%zionpsp = 1
    token%pspheads(:)%pspxc   = token%dtsets(1)%ixc
    token%pspheads(:)%pspso   = 0
    token%pspheads(:)%xccc    = 0
 end if

 !If (all) pspcod are 7 then this is a PAW calculation. Initialize (default) the value of ratsph
 do idtset=0,ndtset_alloc
    token%dtsets(idtset)%usepaw=usepaw
    if(usepaw==0)then
      token%dtsets(idtset)%ratsph(:)=two
    else
!     Note that the following coding assumes that npsp=ntypati for PAW, which is true as of now (XG20101024).
      !token%dtsets(idtset)%ratsph(1:npsp)=token%pspheads(1:npsp)%pawheader%rpaw
      do ipsp=1,npsp
        token%dtsets(idtset)%ratsph(ipsp)=token%pspheads(ipsp)%pawheader%rpaw
      end do
    endif
 end do

 !Take care of other dimensions, and part of the content of dtsets
 !that is or might be needed early.
 !zion_max=maxval(pspheads(1:npsp)%zionpsp) ! This might not work properly with HP compiler

! zion_max=token%pspheads(1)%zionpsp
! do ii=1,npsp
!    zion_max=max(token%pspheads(ii)%zionpsp,zion_max)
! end do
 ABI_ALLOCATE(zionpsp,(npsp))
 do ii=1,npsp
  zionpsp(ii)=token%pspheads(ii)%zionpsp
 end do

 if (AB_DBG) write(std_err,*) "AB module: OK."

 ABI_ALLOCATE(mband_upper_ ,(  0:ndtset_alloc))

 if (AB_DBG) write(std_err,*) "AB module: call invars1m()."

! write(std_out,*)' ab7_invars_f90 , before invars1m : token%pspheads(1)%nproj(0:3)=',token%pspheads(1)%nproj(0:3)

 call invars1m(token%dmatpuflag,token%dtsets,ab_out,lenstr,mband_upper_,&
   & msym,token%mxga_n_rules,token%mxgw_nqlwl,token%mxlpawu,&
   & token%mxmband_upper,&
   & token%mxnatom,token%mxnatpawu,token%mxnatsph,token%mxnatsph_extra,&
   & token%mxnatvshift,&
   & token%mxnconeq,token%mxnimage,token%mxn_efmas_dirs,token%mxnkpt,token%mxnkptgw,token%mxnkpthf,token%mxnnos,&
   & token%mxnqptdm,&
   & token%mxnspinor,token%mxnsppol,token%mxnsym,token%mxntypat,token%mxnimfrqs,&
   & token%mxnfreqsp,token%mxnzchempot,token%mxn_projection_frequencies,ndtset,&
   & ndtset_alloc,string,npsp,zionpsp)

 ABI_DEALLOCATE(zionpsp)
 if (call_timab) then
   call timab(42,2,opt_timab_tsec)
 end if

 if (call_timab) then
   call timab(43,3,opt_timab_tsec)
 end if

 if (AB_DBG) write(std_err,*) "AB module: OK."

 !9) Provide defaults for the variables that have not yet been initialized.
 if (AB_DBG) write(std_err,*) "AB module: call indefo()."
 if (call_status) then
   call status(0,opt_status_file,99,level,'call indefo   ')
 end if

 nprocs = xmpi_comm_size(xmpi_world)
 call indefo(token%dtsets,ndtset_alloc,nprocs)
 if (AB_DBG) write(std_err,*) "AB module: OK."

 if (call_status) then
   call status(0,opt_status_file,99,level,'call macroin  ')
 end if

 call macroin(token%dtsets,ecut_tmp,lenstr,ndtset_alloc,string)

 !10) Perform some global initialization, depending on the value of
 ! pseudopotentials, parallelism variables, or macro input variables

 !If all the pseudopotentials have the same pspxc, override the default
 !value for dtsets 1 to ndtset
 if(with_psp .and. minval(abs((token%pspheads(1:npsp)%pspxc-token%pspheads(1)%pspxc)))==0)then
    token%dtsets(1:ndtset_alloc)%ixc=token%pspheads(1)%pspxc
 end if

 !11) Call the main input routine.
 if (AB_DBG) write(std_err,*) "AB module: call invars2()."

 if (call_status) then
   call status(0,opt_status_file,99,level,'call invars2m ')
 end if

 if (with_mem) then
   !write(std_out,*)' ab7_invars_f90 : token%pspheads(1)%nproj(0:3)=',token%pspheads(1)%nproj(0:3)
   call invars2m(token%dtsets,ab_out,lenstr,mband_upper_,msym,ndtset,ndtset_alloc,npsp,token%pspheads,string)
 else
   do idtset = 1, ndtset_alloc, 1
      jdtset=token%dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=0
      call invars2(token%dtsets(idtset)%bravais,token%dtsets(idtset),ab_out,jdtset,lenstr,&
         & mband_upper_(idtset),msym,npsp,string,usepaw,&
         & token%pspheads(1:npsp)%zionpsp)
   end do
 end if

 if (AB_DBG) write(std_err,*) "AB module: OK."

 if (call_status) then
   call status(0,opt_status_file,99,level,'call macroin2  ')
 end if

 call macroin2(token%dtsets,ndtset_alloc)

 !mxmband=maxval(dtsets(1:ndtset_alloc)%mband) ! This might not work with the HP compiler
 token%mxmband=token%dtsets(1)%mband
 do ii=1,ndtset_alloc
    token%mxmband=max(token%dtsets(ii)%mband,token%mxmband)
 end do

 if (call_timab) then
   call timab(43,2,opt_timab_tsec)
 end if

 ABI_DEALLOCATE(mband_upper_)

 if (call_status) then
   call status(0,opt_status_file,99,level,'exit')
 end if

end subroutine ab7_invars_load
!!***

!----------------------------------------------------------------------

!!****f* m_ab7_invars/ab7_invars_free
!! NAME
!!  ab7_invars_free
!!
!! FUNCTION
!!  Free the memory allocated for a dataset.
!!
!! INPUTS
!!   dtsetsId=Dataset identifier
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!
!! SOURCE

  subroutine ab7_invars_free(dtsetsId)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ab7_invars_free'
!End of the abilint section

 integer, intent(in) :: dtsetsId

 type(dtsets_list), pointer :: token

 nullify(token)
 call get_token(token, dtsetsId)
 call free_token(token)

end subroutine ab7_invars_free
!!***

!----------------------------------------------------------------------

!!****f* m_ab7_invars/ab7_invars_get_ndtset
!! NAME
!!  ab7_invars_get_ndtset
!!
!! FUNCTION
!!   Returns the number of datasets.
!!
!! INPUTS
!!   dtsetsId=Dataset identifier
!!
!! OUTPUT
!!  value=Number of datasets
!!  errno=Error code
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ab7_invars_get_ndtset(dtsetsId, value, errno)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ab7_invars_get_ndtset'
!End of the abilint section

 integer, intent(in) :: dtsetsId
 integer, intent(out) :: value
 integer, intent(out) :: errno

 type(dtsets_list), pointer :: token

 call get_token(token, dtsetsId)
 if (associated(token)) then
    value = size(token%dtsets) - 1
    errno = AB7_NO_ERROR
 else
    errno = AB7_ERROR_OBJ
 end if

end subroutine ab7_invars_get_ndtset
!!***

!!****f* m_ab7_invars/iofn2
!! NAME
!! iofn2
!!
!! FUNCTION
!! First, read and echo pseudopotential filenames from standard input unit
!! Store them in an array.
!!
!! INPUTS
!!  npsp=number of pseudopotentials
!!
!! OUTPUT
!!  pspheads(npsp)=<type pspheader_type>=all the important information from the
!!   pseudopotential file headers, as well as the psp file names
!!
!! PARENTS
!!      m_ab7_invars_f90
!!
!! CHILDREN
!!
!! SOURCE

subroutine iofn2(npsp,filnam)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'iofn2'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npsp
!arrays
 character(len=fnlen), intent(out) :: filnam(npsp)

!Local variables-------------------------------
!scalars
 integer :: ios,ipsp
 character(len=500) :: message
 character(len=fnlen) :: filpsp

!*************************************************************************

 do ipsp=1,npsp
   ! Read the name of the psp file
   write(std_out,'(/,a)' )' Please give name of formatted atomic psp file'
   read (std_in, '(a)' , iostat=ios ) filpsp
   filnam(ipsp)=trim(filpsp)

   !  It might be that a file name is missing
   if (ios/=0) then
     write(message, '(a,a,a,a,a,a,a)' )&
&     'There are not enough names of pseudopotentials',ch10,&
&     'provided in the files file.',ch10,&
&     'Action: check first the variable ntypat (and/or npsp) in the input file;',ch10,&
&     'if they are correct, complete your files file.'
     MSG_ERROR(message)
   end if

   write(std_out,'(a,i0,2a)' )' iofn2 : for atom type ',ipsp,', psp file is ',trim(filpsp)
 end do ! ipsp=1,npsp

end subroutine iofn2
!!***

#include "ab7_invars_f90_get.f90"

end module m_ab7_invars
!!***
