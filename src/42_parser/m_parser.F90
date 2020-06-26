!!****m* ABINIT/m_parser
!! NAME
!! m_parser
!!
!! FUNCTION
!! This module contains (low-level) procedures to parse and validate input files.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2020 ABINIT group (XG, MJV, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_parser

 use defs_basis
 use m_abicore
 use m_errors
 use m_atomdata
 use m_xmpi
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_nctk
 !use m_nctk,      only : write_var_netcdf    ! FIXME Deprecated

 use m_io_tools,  only : open_file
 use m_fstrings,  only : sjoin, strcat, itoa, inupper, ftoa, tolower, toupper, next_token, &
                         endswith, char_count, find_digit !, startswith,
 use m_geometry,  only : xcart2xred, det3r, mkrdim

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* defs_abitypes/ab_dimensions
!! NAME
!! ab_dimensions
!!
!! FUNCTION
!! One record for each dimension of arrays used in ABINIT.
!! Will be used to e.g.:
!! - contain the maximum size attained over all datasets (mxvals)
!! - indicate whether this dimension is the same for all datasets or not (multivals).
!! Used for example inside outvars
!!
!! SOURCE

 type,public :: ab_dimensions

    integer :: ga_n_rules   ! maximal value of input ga_n_rules for all the datasets
    integer :: gw_nqlwl     ! maximal value of input gw_nqlwl for all the datasets
    integer :: lpawu        ! maximal value of input lpawu for all the datasets
    integer :: mband
    integer :: mband_upper ! maximal value of input nband for all the datasets
                           ! Maybe this one could be removed
    integer :: natom
    integer :: natpawu     ! maximal value of number of atoms on which +U is applied for all the datasets
    integer :: natsph      ! maximal value of input natsph for all the datasets
    integer :: natsph_extra  ! maximal value of input natsph_extra for all the datasets
    integer :: natvshift   ! maximal value of input natvshift for all the datasets
    integer :: nberry = 20 ! This is presently a fixed value. Should be changed.
    integer :: nbandhf
    integer :: nconeq      ! maximal value of input nconeq for all the datasets
    integer :: n_efmas_dirs
    integer :: nfreqsp
    integer :: n_projection_frequencies
    integer :: nimage
    integer :: nimfrqs
    integer :: nkpt       ! maximal value of input nkpt for all the datasets
    integer :: nkptgw     ! maximal value of input nkptgw for all the datasets
    integer :: nkpthf     ! maximal value of input nkpthf for all the datasets
    integer :: nnos       ! maximal value of input nnos for all the datasets
    integer :: nqptdm     ! maximal value of input nqptdm for all the datasets
    integer :: nshiftk
    integer :: nsp
    integer :: nspinor    ! maximal value of input nspinor for all the datasets
    integer :: nsppol     ! maximal value of input nsppol for all the datasets
    integer :: nsym       ! maximum number of symmetries
    integer :: ntypalch
    integer :: ntypat     ! maximum number of types of atoms
    integer :: nzchempot  ! maximal value of input nzchempot for all the datasets

 end type ab_dimensions
!!***

 public :: parsefile
 public :: inread
 public :: instrng
 public :: incomprs
 public :: intagm
 public :: importxyz

 public :: chkdpr         ! Checks the value of an input real(dp) variable.
 public :: chkint         ! Checks the value of an input integer variable.
 public :: chkint_eq      ! Checks the value of an input integer variable against a list.
 public :: chkint_ge      ! Checks the value of an input integer variable, expected to be greater than some value.
 public :: chkint_le      ! Checks the value of an input integer variable, expected to be lower than some value.
 public :: chkint_ne      ! Checks the value of an input integer variable against a list.
 !public :: chkint_prt

 public :: prttagm             ! Print the content of intarr or dprarr.
 public :: prttagm_images      ! Extension to prttagm to include the printing of images information.
 public :: chkvars_in_string   ! Analyze variable names in string. Abort if name is not recognized.
 public :: get_acell_rprim     ! Get acell and rprim from string


!----------------------------------------------------------------------

!!****t* m_parser/geo_t
!! NAME
!! geo_t
!!
!! FUNCTION
!!  Small object describing the crystalline structure read from an external file
!!  or a string given in the input file.
!!
!! SOURCE

 type,public :: geo_t

  integer :: natom = 0
  ! Number of atoms

  integer :: ntypat = 0
  ! Number of type of atoms

  character(len=500) :: title = ""
  ! Optional title read for external file e.g. POSCAR

  character(len=500) :: fileformat = ""
  ! (poscar, netcdf)

  integer,allocatable :: typat(:)
  ! typat(natom)
  ! Type of each natom.

  real(dp) :: rprimd(3,3)

  real(dp),allocatable :: xred(:,:)
  ! xred(3,natom)
  ! Reduced coordinates.

  real(dp),allocatable :: znucl(:)
  ! znucl(ntypat)
  ! Nuclear charge for each type of pseudopotential
  ! Note that ntypat must be equal to npsp --> no alchemical mixing

 contains

   procedure :: free => geo_free
   ! Free memory.

   procedure :: malloc => geo_malloc
   ! Allocate memory

   procedure :: bcast => geo_bcast
   ! Brodcast object

   procedure :: print_abivars => geo_print_abivars
   !  Print Abinit variables corresponding to POSCAR

 end type geo_t

 public :: geo_from_abivar_string   ! Build object form abinit variable
 public :: geo_from_poscar_path     ! Build object from POSCAR filepath.

 public :: intagm_img   !  Read input file variables according to images path definition (1D array)

 interface intagm_img
   module procedure intagm_img_1D
   module procedure intagm_img_2D
 end interface intagm_img


CONTAINS  !===========================================================
!!***

!!****f* m_parser/parsefile
!! NAME
!! parsefile
!!
!! FUNCTION
!!  Glue function, to read the given file, put it into a string,
!!  change everything to uppercase, remove carriage returns and
!!  non significant blank characters. May also read a XYZ input
!!  file if specified. Finally read ndtset input variable.
!!
!! INPUTS
!!  filnamin= the file to read
!!  comm=MPI communicator
!!
!! OUTPUT
!!  lenstr= the length of the resulting string.
!!  ndtset= the number of declared datasets.
!!  string= contains on output the content of the file, ready for parsing.
!!
!! PARENTS
!!      abinit,m_ab7_invars_f90,ujdet
!!
!! CHILDREN
!!      importxyz,instrng,intagm,inupper,xmpi_bcast
!!
!! SOURCE

subroutine parsefile(filnamin, lenstr, ndtset, string, comm)

!Arguments ------------------------------------
 character(len=*),intent(in) :: filnamin
 integer,intent(in) :: comm
 integer,intent(out) :: ndtset,lenstr
 character(len=strlen),intent(out) :: string

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0, option1= 1
 integer :: marr,tread,lenstr_noxyz,ierr
 character(len=strlen) :: string_raw
 character(len=500) :: msg
!arrays
 integer :: intarr(1)
 real(dp) :: dprarr(1)

! *************************************************************************

 ! Read the input file, and store the information in a long string of characters
 ! Note: this is done only by me=0, and then string and other output vars are BCASTED

 if (xmpi_comm_rank(comm) == master) then

   ! strlen from defs_basis module
   call instrng(filnamin, lenstr, option1, strlen, string)

   ! Copy original file, without change of case
   string_raw=string

   ! To make case-insensitive, map characters of string to upper case.
   call inupper(string(1:lenstr))

   ! Might import data from xyz file(s) into string
   ! Need string_raw to deal properly with xyz filenames
   lenstr_noxyz = lenstr
   call importxyz(lenstr,string_raw,string,strlen)

   ! Make sure we don't have unmatched quotation marks
   if (mod(char_count(string, '"'), 2) /= 0) then
     MSG_ERROR('Your input file contains unmatched quotation marks `"`. This confuses the parser. Check your input.')
   end if

   ! Take ndtset from the input string
   ndtset=0; marr=1
   call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),"ndtset",tread,'INT')
   if (tread==1) ndtset=intarr(1)
   ! Check that ndtset is within bounds
   if (ndtset<0 .or. ndtset>9999) then
     write(msg, '(a,i0,4a)' )&
     'Input ndtset must be non-negative and < 10000, but was ',ndtset,ch10,&
     'This is not allowed.',ch10,'Action: modify ndtset in the input file.'
     MSG_ERROR(msg)
   end if
 end if ! master

 if (xmpi_comm_size(comm) > 1) then
   ! Broadcast data.
   call xmpi_bcast(lenstr, master, comm, ierr)
   call xmpi_bcast(ndtset, master, comm, ierr)
   call xmpi_bcast(string, master, comm, ierr)
   call xmpi_bcast(string_raw, master, comm, ierr)
 end if

 ! Save input string in global variable so that we can access it in ntck_open_create
 INPUT_STRING = string_raw

end subroutine parsefile
!!***

!!****f* m_parser/inread
!! NAME
!! inread
!!
!! FUNCTION
!! Carry out internal read from input character string, starting
!! at first character in string, reading ndig digits (including possible
!! sign, decimal, and exponent) by computing the appropriate format and
!! performing a formatted read (list-directed read would be perfect for
!! this application but is inconsistent with internal read according to
!! Fortran90 standard).
!! In case of a real number, this routine
!! is also able to read SQRT(number): return the square root of the number.
!!
!! INPUTS
!!  string=character string.
!!  ndig=length of field to be read (including signs, decimals, and exponents).
!!  typevarphys=variable type (might indicate the physical meaning for dimensionality purposes)
!!   'INT'=>integer
!!   'DPR','LEN','ENE'=>real(dp) (no special treatment)
!!   'LOG'=>integer, but read logical variable T,F,.true., or .false.
!!   'KEY'=>character, returned in token
!!
!! OUTPUT
!!  outi or outr (integer or real respectively)
!!  errcod, =0 for success, 1,2 for ini, inr failure resp.
!!
!! PARENTS
!!      adini,inarray
!!
!! CHILDREN
!!
!! SOURCE

subroutine inread(string,ndig,typevarphys,outi,outr,errcod)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndig
 integer,intent(out) :: errcod,outi
 real(dp),intent(out) :: outr
 character(len=*),intent(in) :: string
 character(len=*),intent(in) :: typevarphys

!Local variables-------------------------------
!scalars
 integer :: done,idig,index_slash,sign
 real(dp) :: den,num
 logical :: logi
 character(len=500) :: msg,iomsg

! *************************************************************************

 !write(std_out,*)'inread: enter with string(1:ndig): ',string(1:ndig)
 !write(std_out,*)'typevarphys: ',typevarphys

 if (typevarphys=='INT') then

   ! integer input section
   read(unit=string(1:ndig), fmt=*, iostat=errcod, iomsg=iomsg) outi

   if(errcod/=0)then
     ! integer reading error
     write(msg,'(a,i0,7a)' ) &
       "Attempted to read ndig: ",ndig," integer digits", ch10, &
       "from string(1:ndig)= `",string(1:ndig),"` to initialize an integer variable",ch10,&
       "iomsg: ", trim(iomsg)
     MSG_WARNING(msg)
     errcod=1
   end if

 else if (typevarphys=='DPR' .or. typevarphys=='LEN' .or. typevarphys=='ENE' &
         .or. typevarphys=='BFI' .or. typevarphys=='TIM') then

   ! real(dp) input section
   ! Special treatment of SQRT(xxx) or -SQRT(xxx) chains of characters, where xxx can be a fraction
   done=0
   if (ndig>5) then
     if(string(1:5)=='SQRT(' .and. string(ndig:ndig)==')')then
       done=1 ; sign=1
     else if(string(1:6)=='-SQRT(' .and. string(ndig:ndig)==')')then
       done=1 ; sign=2
     end if

     if(done==1)then
       index_slash=index(string(5+sign:ndig-1),'/')
       if(index_slash==0)then
         read (unit=string(5+sign:ndig-1),fmt=*,iostat=errcod, iomsg=iomsg) outr
       else if(index_slash/=0)then
         read (unit=string(5+sign:5+sign+index_slash-2),fmt=*,iostat=errcod, iomsg=iomsg) num
         if(errcod==0)then
           read (unit=string(5+sign+index_slash:ndig-1),fmt=*,iostat=errcod, iomsg=iomsg) den
           if(errcod==0)then
             if(abs(den)<tol12)then
               errcod=1
             else
               outr=num/den
             end if
           end if
         end if
       end if
       if(outr<-tol12)then
         errcod=1
       else
         outr=sqrt(outr)
         if(sign==2)outr=-outr
       end if
     end if
   end if

   ! Special treatment of fractions
   if(done==0)then
     index_slash=index(string(1:ndig),'/')
     if(index_slash/=0)then
       done=1
       read (unit=string(1:index_slash-1), fmt=*, iostat=errcod, iomsg=iomsg) num
       if(errcod==0)then
         read (unit=string(index_slash+1:ndig), fmt=*, iostat=errcod, iomsg=iomsg) den
         if(errcod==0)then
           if(abs(den)<tol12)then
             errcod=1
           else
             outr=num/den
           end if
         end if
       end if
     end if
   end if

   ! Normal treatment of floats
   if(done==0) read (unit=string(1:ndig), fmt=*, iostat=errcod, iomsg=iomsg) outr

   ! Treatment of errors
   if(errcod/=0)then
     ! real(dp) data reading error
     write(msg,'(a,i0,8a)' ) &
        'Attempted to read ndig: ',ndig,' floating point digits,',ch10, &
        'from string(1:ndig): `',string(1:ndig),'` to initialize a floating variable.',ch10, &
        "iomsg: ", trim(iomsg)
     MSG_WARNING(msg)
     errcod=2
   end if

 else if (typevarphys=='LOG') then

   read (unit=string(1:ndig), fmt=*, iostat=errcod, iomsg=iomsg) logi

   if(errcod/=0)then
     ! integer reading error
     write(msg,'(a,i0,8a)' ) &
       "Attempted to read ndig: ",ndig," integer digits", ch10, &
       "from string(1:ndig): `",string(1:ndig),"` to initialize a logical variable.",ch10,&
       "iomsg: ", trim(iomsg)
     MSG_WARNING(msg)
     errcod=3
   end if

   if(logi)outi=1
   if(.not.logi)outi=0

 else
   write(msg,'(4a)' ) &
   'Argument typevarphys must be INT, DPR, LEN, ENE, BFI, TIM or LOG ',ch10,&
   'but input value was: ',trim(typevarphys)
   MSG_ERROR(msg)
 end if

 if (errcod /= 0)then
   do idig=1,ndig
     if( string(idig:idig) == 'O' )then
       write(msg,'(3a)' ) &
       'Note that this string contains the letter O. ',ch10,&
       'It is likely that this letter should be replaced by the number 0.'
       MSG_WARNING(msg)
       exit
     end if
   end do
 end if

end subroutine inread
!!***

!!****f* m_parser/instrng
!! NAME
!! instrng
!!
!! FUNCTION
!! Read the input file, and product a string of character,
!! with all data, to be analyzed in later routines. The length
!! of this string is lenstr. This number is checked to be smaller
!! than the dimension of the string of character, namely strln.
!!
!! INPUTS
!!  filnam=name of the input file, to be read
!!  option= if 0, simple storing of the character string,
!!             no special treatment for ABINIT (comment delimiters, checks, include ...)
!!          if 1, suppresses text after an ABINIT comment delimiter (! or #),
!!             checks that a minus sign is followed by a number ...
!!                check for INCLUDE statement:
!!                if present, add string from included file
!!  strln=maximal number of character of string, as declared in the calling routine
!!
!! OUTPUT
!!  lenstr=actual number of character in string
!!  string*(strln)=string of character
!!
!! PARENTS
!!      anaddb,importcml,localorb_S,lwf,parsefile
!!
!! CHILDREN
!!      incomprs,instrng,wrtout
!!
!! SOURCE

recursive subroutine instrng(filnam, lenstr, option, strln, string)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: option,strln
 integer,intent(out) :: lenstr
 character(len=*),intent(in) :: filnam
 character(len=*),intent(out) :: string

!Local variables-------------------------------
 character :: blank=' '
!scalars
 integer,save :: include_level=-1
 integer :: ii,ii1,ii2,ij,iline,ios,iost,lenc,lenstr_inc,mline,nline1,input_unit
 logical :: include_found, ex
 character(len=1) :: string1
 character(len=3) :: string3
 character(len=500) :: filnam_inc,msg
 character(len=fnlen+20) :: line
 character(len=strlen),pointer :: string_inc

!************************************************************************

 DBG_ENTER("COLL")

!%%%%%%%%%%%%%%%%%%%%%%%%
!read in string from file
!%%%%%%%%%%%%%%%%%%%%%%%%

 ! The file can be included in another (prevent too many include levels)
 include_level=include_level+1
 if (include_level>2) then
   write(msg, '(3a)' ) &
   'At least 4 levels of included files are present in input file !',ch10,&
   'This is not allowed. Action: change your input file.'
   MSG_ERROR(msg)
 end if

 ! Open data file and read one line at a time, compressing data
 ! and concatenating into single string:
 if (open_file(filnam,msg,newunit=input_unit,form="formatted",status="old",action="read") /= 0) then
   MSG_ERROR(msg)
 end if
 rewind (unit=input_unit)

 ! Initialize string to blanks
 string=blank
 lenstr=1

 ! Set maximum number lines to be read to some large number
 mline=500000
 do iline=1,mline

   ! Keeps reading lines until end of input file
   read (unit=input_unit,fmt= '(a)' ,iostat=ios) line(1:fnlen+20)
   !  Hello ! This is a commentary. Please, do not remove me.
   !  In fact, this commentary protect tests_v4 t47 for miscopying
   !  the input file into the output string. It _is_ strange.
   !  The number of lines in the commentary is also resulting from
   !  a long tuning..

   !  DEBUG
   !  write(std_out,*)' instrng, iline=',iline,' ios=',ios,' echo :',trim(line(1:fnlen+20))
   !  ENDDEBUG

   ! Exit the reading loop when arrived at the end
   if(ios/=0)then
     backspace(input_unit)
     read (unit=input_unit,fmt= '(a1)' ,iostat=ios) string1
     if(ios/=0)exit
     backspace(input_unit)
     read (unit=input_unit,fmt= '(a3)' ,iostat=ios) string3
     if(string3=='end')exit
     write(msg, '(3a,i0,11a)' ) &
      'It is observed in the input file: ',TRIM(filnam),', line number ',iline,',',ch10,&
      'that there is a non-zero IO signal.',ch10,&
      'This is normal when the file is completely read.',ch10,&
      'However, it seems that the error appears while your file has not been completely read.',ch10,&
      'Action: correct your file. If your file seems correct, then,',ch10,&
      'add the keyword ''end'' at the very beginning of the last line of your input file.'
     MSG_ERROR(msg)
   end if

   ! TODO: Ignore sections inside TEST_INFO markers so that we don't need to prepend comment markers.
   !in_testinfo = 0
   !if startswith(line, "#%%<BEGIN TEST_INFO") in_testinfo = 1
   !if (in_testinfo /= 0) cycle
   !if startswith(line, "#%%<END TEST_INFO> ") then
   !  in_testinfo = 0; cycle
   !end if

   ! Find length of input line ignoring delimiter characters (# or !)
   ! and any characters beyond it (allows for comments beyond # or !)
   ii1=index(line(1:fnlen+20),'#')
   ii2=index(line(1:fnlen+20),'!')
   if ( (ii1==0 .and. ii2==0) .or. option==0 ) then
     ! delimiter character was not found on line so use full line
     ii=fnlen+20
   else if(ii1==0)then
     ! ii will represent length of line up to but not including !
     ii=ii2-1
   else if(ii2==0)then
     ! ii will represent length of line up to but not including #
     ii=ii1-1
   else
     ii=min(ii1,ii2)-1
   end if

   ! Checks that nothing is left beyond fnlen
   if(ii>fnlen)then
     do ij=fnlen+1,ii
       if(line(ij:ij)/=' ')then
         write(msg,'(3a,i0,3a,i0,3a)' ) &
          'It is observed in the input file: ',TRIM(filnam),' line number ',iline,',',ch10,&
          'that more than ',fnlen,' columns are used.',ch10,&
          'This is not allowed. Change this line of your input file.'
         MSG_ERROR(msg)
       end if
     end do
   end if

   if (ii>0) then
     ! Check for the occurence of a minus sign followed by a blank
     ij=index(line(1:ii),'- ')
     if (ij>0 .and. option==1) then
       write(msg, '(3a,i0,11a)' ) &
       'It is observed in the input file:, ',TRIM(filnam),' line number ',iline,',',ch10,&
       'the occurence of a minus sign followed',ch10,&
       'by a blank. This is forbidden.',ch10,&
       'If the minus sign is meaningful, do not leave a blank',ch10,&
       'between it and the number to which it applies.',ch10,&
       'Otherwise, remove it.'
       MSG_ERROR(msg)
     end if
     ! Check for the occurence of a tab
     ij=index(line(1:ii),char(9))
     if (ij>0 .and. option==1 ) then
       write(msg, '(3a,i0,3a)' ) &
        'The occurence of a tab, in the input file: ',TRIM(filnam),' line number ',iline,',',ch10,&
        'is observed. This sign is confusing, and has been forbidden.'
       MSG_ERROR(msg)
     end if

     ! Check for the occurence of a include statement
     include_found=.false.
     if (option==1) then
       ! Look for include statement
       ii1=index(line(1:ii),"include");ii2=index(line(1:ii),"INCLUDE")
       include_found=(ii1>0.or.ii2>0)
       if (include_found) then
         ij=max(ii1,ii2);ii1=0;ii2=0
         ! Look for quotes (ascii 34)
         ii1=index(line(ij+7:ii),char(34))
         if (ii1>1) ii2=index(line(ij+7+ii1:ii),char(34))
         ! Look for quotes (ascii 39)
         if (ii1==0.and.ii2==0) then
           ii1=index(line(ij+7:ii),char(39))
           if (ii1>1) ii2=index(line(ij+7+ii1:ii),char(39))
         end if
         ! Check if quotes are correctly set
         ex=(ii1<=1.or.ii2<=1)
         if (.not.ex) then
           msg=line(ij+7:ij+5+ii1)
           call incomprs(msg(1:ii1-1),lenc)
           ex=(len(trim(msg))/=0)
         end if
         if (ex) then
           write(msg, '(6a)' ) &
            'A "include" statement has been found in input file: ',TRIM(filnam),ch10,&
            'but there must be a problem with the quotes.',ch10,&
            'Action: change your input file.'
           MSG_ERROR(msg)
         end if
         ! Store included file name
         filnam_inc=line(ij+7+ii1:ij+5+ii1+ii2)
         ! Extract include statement from line
         lenc=ii1+ii2+7
         msg(1:ii-lenc)=line(1:ij-1)//line(ij+lenc:ii)
         ii=ii-lenc;line(1:ii)=msg(1:ii)
       end if
     end if

     ! Compress: remove repeated blanks, make all ASCII characters
     ! less than a blank (and '=') to become a blank.
     call incomprs(line(1:ii),lenc)

   else
     ! ii=0 means line starts with #, is entirely a comment line
     lenc=0;include_found=.false.
   end if

   ! Check resulting total string length
   if (lenstr+lenc>strln) then
     write(msg, '(8a)' ) &
      'The size of your input file: ',TRIM(filnam),' is such that the internal',ch10,&
      'character string that should contain it is too small.',ch10,&
      'Action: decrease the size of your input file,',ch10,&
      'or contact the ABINIT group.'
     MSG_ERROR(msg)
   end if

   if (lenc>0) then
     ! Concatenate new compressed characters
     ! with previous part of compressed string (unless all blank)
     string(lenstr+1:lenstr+lenc)=line(1:lenc)
   end if
   ! Keep track of total string length
   lenstr=lenstr+lenc

   ! Eventually (recursively) read included file
   if (include_found) then
     ! Check file existence
     inquire(file=filnam_inc ,iostat=iost,exist=ex)
     if (.not. ex .or. iost /= 0) then
       write(msg, '(5a)' ) &
        'Input file: ',TRIM(filnam),' reading: the included file ',trim(filnam_inc),' cannot be found !'
       MSG_ERROR(msg)
     end if
     ! Read included file (warning: recursive call !)
     ABI_ALLOCATE(string_inc,)
     call instrng(trim(filnam_inc),lenstr_inc,option,strln-lenstr,string_inc)
     ! Check resulting total string length
     if (lenstr+lenstr_inc>strln) then
       write(msg, '(6a)' ) &
        'The size of your input file: ',TRIM(filnam),' (including included files) is such that',ch10,&
        'the internal character string that should contain it is too small !',ch10,&
        'Action: decrease the size of your input file.'
       MSG_ERROR(msg)
     end if
     ! Concatenate total string
     string(lenstr+1:lenstr+lenstr_inc)=string_inc(1:lenstr_inc)
     lenstr=lenstr+lenstr_inc
     ABI_FREE(string_inc)
   end if

   ! If mline is reached, something is wrong
   if (iline>=mline) then
     write(msg, '(a,i0,2a,i0,4a)' ) &
     'The number of lines already read from input file: ',iline,ch10,&
     'is equal or greater than maximum allowed mline: ',mline,ch10,&
     'Action: you could decrease the length of the input file, or',ch10,&
     'increase mline in this routine.'
     MSG_ERROR(msg)
   end if

 end do !  End loop on iline. Note that there is an "exit" instruction in the loop

 nline1=iline-1
 close (unit=input_unit)

 ! Make sure we don't have unmatched quotation marks
 if (mod(char_count(string, '"'), 2) /= 0) then
   MSG_ERROR('Your input file contains unmatched quotation marks `"`. This confuses the parser. Check your input.')
 end if

 include_level = include_level - 1

 write(msg,'(a,i0,3a)')'-instrng: ',nline1,' lines of input have been read from file ',trim(filnam),ch10
 call wrtout(std_out,msg)
 !write(std_out, "(3a)")"string after instrng:", ch10, trim(string)

 DBG_EXIT("COLL")

end subroutine instrng
!!***

!!****f* m_parser/inreplsp
!! NAME
!! inreplsp
!!
!! FUNCTION
!! Replace all occurrences of characters lexically less than SP (blank)
!! by SP in the input string, returning modified string of same length.
!! Also replace a '=' by a SP.
!!
!! INPUTS
!!  string=character string to be modified
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  string=same character string with ASCII (decimal) 0-31 replaced by 32.
!!
!! PARENTS
!!      incomprs
!!
!! CHILDREN
!!
!! SOURCE

subroutine inreplsp(string)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(inout) :: string

!Local variables-------------------------------
!scalars
 integer :: ilenth,length

! *************************************************************************

 ! Get length of string. Proceed only if string has nonzero length
 length=len(string); if (length == 0) return

 !  Do replacement by going through input character string one character at a time
 do ilenth=1,length
   if (llt(string(ilenth:ilenth),' ')) string(ilenth:ilenth)=' '
   if (string(ilenth:ilenth)=='=') string(ilenth:ilenth)=' '
 end do

end subroutine inreplsp
!!***

!!****f* m_parser/incomprs
!! NAME
!! incomprs
!!
!! FUNCTION
!! Compresses input character string into the following form:
!! (1) Replaces tabs and all other characters lexically less than
!! SP (blank) with SP (blank), where lexically less than refers to
!! the ASCII collating sequence (SP is hex 20, dec 32).
!! The use of llt is needed e.g. on the IBM 9000 because it does not
!! handle tab characters sensibly in its AIX fortran.
!! Also replace occurences of '=' by a SP.
!! (2) Removes all repeated blanks, ignoring trailing blanks
!! after first (returns nontrailing final length in arg 'length').
!! (3) Makes first character in string NONBLANK.  This is done
!! to prevent double blanks from occurring when compressed string
!! is concatenated with other compressed strings.
!! (4) Makes last character (string(length:length)) a blank.
!! If input string is entirely blank or tabs, simply returns with length=0.
!!
!! INPUTS
!!  (see side effects)
!!
!! OUTPUT
!!  length=nonblank, nontab length of string as defined above
!!
!! SIDE EFFECT
!!  string=at input:  character string
!!         at output: repeated blanks and tabs have been removed and
!!                    remaining tabs have been replaced by blanks
!!
!! PARENTS
!!      importxyz,instrng
!!
!! CHILDREN
!!      inreplsp
!!
!! SOURCE

subroutine incomprs(string,length)

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: length
 character(len=*),intent(inout) :: string

!Local variables-------------------------------
 character(len=1) :: blank=' '
!scalars
 integer :: bb,f1,ii,jj,kk,l1,lbef,lcut,lold,stringlen
 character(len=500) :: msg

! *************************************************************************

 ! String length determined by calling program declaration of "string"
 stringlen=len(string)
 length=stringlen

 ! Only proceed if string has nonzero length
 if (length>0) then
   ! Find last nonblank character (i.e. nonblank and nontab length)
   length=len_trim(string)
   if (length==0) then
     ! Line is all blanks or tabs so do not proceed
     ! write(std_out,*)' incomprs: blank line encountered'
   else

     ! Replace all characters lexically less than SP, and '=', by SP (blank)
     call inreplsp(string(1:length))

     ! Continue with parsing
     ! l1 is set to last nonblank, nontab character position
     l1=length
     do ii=1,l1
       if (string(ii:ii)/=blank) exit
     end do

     ! f1 is set to first nonblank, nontab character position
     f1=ii
     ! lbef is number of characters in string starting at
     ! first nonblank, nontab and going to last
     lbef=l1-f1+1

     ! Process characters one at a time from right to left:
     bb=0
     lcut=lbef
     do ii=1,lbef
       jj=lbef+f1-ii
       ! set bb=position of next blank coming in from right
       if (string(jj:jj)==blank) then
         if (bb==0) bb=jj
       else
         if (bb/=0) then
           ! if several blanks in a row were found, cut from string
           if (jj<bb-1) then
             ! lold becomes string length before cutting blanks
             lold=lcut
             ! lcut will be new string length
             lcut=lcut-(bb-1-jj)
             ! redefine string with repeated blanks gone
             do kk=1,f1+lcut-1-jj
               string(jj+kk:jj+kk)=string(kk+bb-1:kk+bb-1)
             end do
           end if
           bb=0
         end if
       end if
     end do

     ! Remove initial blanks in string if any
     if (f1>1) string(1:lcut)=string(f1:f1+lcut-1)

     ! Add blank on end unless string had no extra space
     if (lcut==stringlen) then
       write(msg,'(a,i7,a,a,a,a,a,a,a,a)')&
       'For input file, with data forming a string of',stringlen,' characters,',ch10,&
       'no double blanks or tabs were found.',ch10,&
       'This is unusual for an input file (or any file),',ch10,&
       'and may cause parsing trouble.  Is this a binary file?',ch10
       MSG_WARNING(msg)
     else
       length=lcut+1
       string(length:length)=blank
     end if
   end if
 end if

end subroutine incomprs
!!***

!!****f* m_parser/intagm
!! NAME
!! intagm
!!
!! FUNCTION
!! Search input 'string' for specific 'token'. Search depends on
!! input dataset through 'jdtset'. Then, return the information mentioned after 'token'.
!! See the "notes" section
!!
!! INPUTS
!!  jdtset=see the notes section
!!  marr=dimension of the intarr and dprarr arrays, as declared in the calling subroutine.
!!  narr=actual size of array to be read in.
!!  string=character string containing 'tags' and data.
!!  token=character string for 'tag'.
!!  typevarphys= variable type (might indicate the physical meaning of for dimensionality purposes)
!!   'INT'=>integer
!!   'DPR'=>real(dp) (no special treatment)
!!   'LEN'=>real(dp) (expect a "length", identify bohr, au, nm or angstrom,
!!       and return in au -atomic units=bohr- )
!!   'ENE'=>real(dp) (expect a "energy", identify Ha, hartree, eV, Ry, Rydberg)
!!   'LOG'=>integer, but read logical variable T,F,.true., or .false.
!!   'KEY'=>character, returned in key_value
!!
!! OUTPUT
!!  intarr(1:narr), dprarr(1:narr)
!!   integer or real(dp) arrays, respectively (see typevarphys),
!!   into which data is read if typevarphys/='KEY'. Use these arrays even for scalars.
!!  tread is an integer: tread = 0 => no data was read
!!                       tread = 1 => data was read
!!  ds_input is an optional integer flag:
!!           ds_input = 0 => value was found which is not specific to jdtset
!!           ds_input > 0 => value was found which is specific to jdtset
!!   one could add more information, eg whether a ? or a : was used, etc...
!!   [key_value]=Stores the value of key if typevarphys=="KEY".
!!      The string must be large enough to contain the output. fnlen is OK in many cases
!!      except when reading a list of files. The routine aborts is key_value cannot store the output.
!!      Output string is left justified.
!!
!! NOTES
!!
!! If jdtset==0:
!!
!!  Search compressed 'string' for blank//'token'//blank and
!!  read input data beside 'token', to be read into appropriate variable.
!!  For this routine to find a given token, the token has to be preceded
!!  and followed by blanks--i.e. the first token should not start out as
!!  the first character in the input file.  This is checked in the calling
!!  subroutine 'input'. Calls inread which performs internal read from
!!  specified string.  Also calls upper which maps characters to all upper case.
!!  Also checks whether there is an occurence of blank//'token'//digit,
!!  in which case the input file might be erroneous, so stops.
!!
!! If jdtset is a positive number:
!!
!!  (1) First search for modified string, blank//'token'//jdtset//blank
!!
!!  (2a) if the occurence of (1) is not found,
!!       look for other modified strings,
!!       blank//'token'//'?'//unities//blank
!!       or
!!       blank//'token'//dozens//'?'//blank
!!       (issue an error message if more than one occurs)
!!       where jdtset=dozens*10+unities (decimal decomposition of jdtset)
!!       if one of them exists, just take the value
!!       Note that unities is a one-digit number, while dozens might be bigger than 9.
!!
!!  (2b-2c) search for a series, with the following tokens :
!!       (issue an error message if more than one occurs, or
!!       goto (3) if none exist)
!!
!!      blank//'token'//':'//blank
!!      if it exists, then a series might have been defined in the input file
!!      must thus find either the increment, blank//'token'//'+'//blank,
!!      or the multiplicative factor, blank//'token'//'*'//blank
!!
!!      blank//'token'//'?'//':'//blank
!!      if it exists, then a series for the inner loop
!!      might have been defined in the input file
!!      must thus find either the increment, blank//'token'//'?'//'+'//blank,
!!      or the multiplicative factor, blank//'token'//'?'//'*'//blank
!!
!!      blank//'token'//':'//'?'//blank
!!      if it exists, then a series for the outer loop
!!      might have been defined in the input file
!!      must thus find either the increment, blank//'token'//'+'//'?'//blank,
!!      or the multiplicative factor, blank//'token'//'*'//'?'//blank
!!
!!  (3) if neither (1) nor (2) are found, search for the 'normal'
!!       string, blank//'token'//blank
!!
!!
!! PARENTS
!!      ingeo,ingeobld,inkpts,inqpt,invacuum,invars0,invars1,invars2
!!      m_ab7_invars_f90,m_anaddb_dataset,m_band2eps_dataset,m_parser
!!      m_multibinit_dataset,m_scup_dataset,macroin,mpi_setup,parsefile,ujdet
!!
!! CHILDREN
!!      appdig,inarray,inupper,wrtout
!!
!! SOURCE

subroutine intagm(dprarr,intarr,jdtset,marr,narr,string,token,tread,typevarphys,ds_input,key_value)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: jdtset,marr,narr
 integer,intent(out) :: tread
 integer,intent(out),optional :: ds_input
 character(len=*),intent(in) :: string
 character(len=*),intent(in) :: token
 character(len=*),intent(in) :: typevarphys
 character(len=*),optional,intent(out) :: key_value
!arrays
 integer,intent(inout) :: intarr(marr)
 real(dp),intent(inout) :: dprarr(marr)

!Local variables-------------------------------
 character(len=1), parameter :: blank=' '
!scalars
 integer :: b1,b2,b3,cs1len,cslen,dozens,ier,itoken,itoken1,itoken2,itoken2_1colon
 integer :: itoken2_1plus,itoken2_1times,itoken2_2colon,itoken2_2plus
 integer :: itoken2_2times,itoken2_colon,itoken2_plus,itoken2_times
 integer :: itoken_1colon,itoken_1plus,itoken_1times,itoken_2colon,itoken_2plus
 integer :: itoken_2times,itoken_colon,itoken_plus,itoken_times,number,opttoken
 integer :: sum_token,toklen,trial_cslen,trial_jdtset,unities
 integer :: ds_input_
 character(len=4) :: appen
 character(len=3) :: typevar
 character(len=500) :: msg
 character(len=fnlen) :: cs,cs1,cs1colon,cs1plus,cs1times,cs2colon,cs2plus
 character(len=fnlen) :: cs2times,cscolon,csplus,cstimes,trial_cs
!arrays
 integer,allocatable :: int1(:),int2(:)
 real(dp),allocatable :: dpr1(:),dpr2(:)

! *************************************************************************

 ABI_CHECK(marr >= narr, sjoin("marr", itoa(marr)," < narr ", itoa(narr), "for token:", token))

 ds_input_ = -1
 dozens=jdtset/10
 unities=jdtset-10*dozens

 if(jdtset<0)then
   write(msg,'(a,i0,a)')' jdtset: ',jdtset,', while it should be non-negative.'
   MSG_ERROR(msg)
 end if

 if(jdtset > 9999)then
   write(msg,'(a,i0,a)')' jdtset: ',jdtset,', while it must be lower than 10000.'
   MSG_ERROR(msg)
 end if

 ! Default values: nothing has been read
 itoken=0
 opttoken=0
 ! Initialise flags in case of opttoken >= 2 later.
 itoken_times=0
 itoken_plus=0
 itoken_colon=0
 cslen=1

 if (narr/=0) then

   toklen=len_trim(token)

   ! --------------------------------------------------------------------------
   ! (1) try to find the token with dataset number appended
   if (jdtset > 0) then

     call appdig(jdtset,'',appen)
     cs=blank//token(1:toklen)//trim(appen)//blank
     if(jdtset<10) then
       cslen=toklen+3
     else if(jdtset<100) then
       cslen=toklen+4
     else if(jdtset<1000) then
       cslen=toklen+5
     else if(jdtset<10000)then
       cslen=toklen+6
     end if
     ! Map token to all upper case (make case-insensitive):
     call inupper(cs)
     ! Absolute index of blank//token//blank in string:
     itoken=index(string,cs(1:cslen))
     ! Look for another occurence of the same token in string, if so, leaves:
     itoken2=index(string,cs(1:cslen), BACK=.true. )
     if(itoken/=itoken2)then
       write(msg, '(7a)' )&
       'There are two occurences of the keyword "',cs(1:cslen),'" in the input file.',ch10,&
       'This is confusing, so it has been forbidden.',ch10,&
       'Action: remove one of the two occurences.'
       MSG_ERROR(msg)
     end if

     if(itoken/=0) then
       opttoken=1
       ds_input_=jdtset
     end if
   end if

   ! --------------------------------------------------------------------------
   ! (2a) try to find the token appended with a string that contains the metacharacter "?".
   if (jdtset>0 .and. opttoken==0)then

     ! Use the metacharacter for the dozens, and save in cs and itoken
     write(appen,'(i1)')unities
     cs=blank//token(1:toklen)//'?'//trim(appen)//blank
     cslen=toklen+4
     ! Map token to all upper case (make case-insensitive):
     call inupper(cs)
     ! Absolute index of blank//token//blank in string:
     itoken=index(string,cs(1:cslen))
     ! Look for another occurence of the same token in string, if so, leaves:
     itoken2=index(string,cs(1:cslen), BACK=.true. )
     if(itoken/=itoken2)then
       write(msg, '(7a)' )&
        'There are two occurences of the keyword: "',cs(1:cslen),'" in the input file.',ch10,&
        'This is confusing, so it has been forbidden.',ch10,&
        'Action: remove one of the two occurences.'
       MSG_ERROR(msg)
     end if
     if(itoken/=0) then
       opttoken=1
       ds_input_=jdtset
     end if

     ! Use the metacharacter for the units, and save in cs1 and itoken1
     write(appen,'(i1)')dozens
     cs1=blank//token(1:toklen)//trim(appen)//'?'//blank
     ! Map token to all upper case (make case-insensitive):
     call inupper(cs1)
     ! Absolute index of blank//token//blank in string:
     itoken1=index(string,cs1(1:cslen))
     ! Look for another occurence of the same token in string, if so, leaves:
     itoken2=index(string,cs1(1:cslen), BACK=.true. )
     if(itoken1/=itoken2)then
       write(msg, '(7a)' )&
       'There are two occurences of the keyword "',cs1(1:cslen),'" in the input file.',ch10,&
       'This is confusing, so it has been forbidden.',ch10,&
       'Action: remove one of the two occurences.'
       MSG_ERROR(msg)
     end if

     if(itoken/=0 .and. itoken1/=0)then
       write(msg, '(9a)' )&
       'The keywords: "',cs(1:cslen),'" and: "',cs1(1:cslen),'"',ch10,&
       'cannot be used together in the input file.',ch10,&
       'Action: remove one of the two keywords.'
       MSG_ERROR(msg)
     end if

     if(itoken1/=0)then
       opttoken=1
       itoken=itoken1
       cs=cs1
       ds_input_=jdtset
     end if

   end if

   ! --------------------------------------------------------------------------
   ! (2b) try to find the tokens defining a series
   if (opttoken==0) then

     cs=token(1:toklen)

     cslen=toklen+3
     cs1len=toklen+4

     cscolon=blank//token(1:toklen)//':'//blank
     csplus=blank//token(1:toklen)//'+'//blank
     cstimes=blank//token(1:toklen)//'*'//blank

     cs1colon=blank//token(1:toklen)//'?'//':'//blank
     cs1plus=blank//token(1:toklen)//'?'//'+'//blank
     cs1times=blank//token(1:toklen)//'?'//'*'//blank

     cs2colon=blank//token(1:toklen)//':'//'?'//blank
     cs2plus=blank//token(1:toklen)//'+'//'?'//blank
     cs2times=blank//token(1:toklen)//'*'//'?'//blank

     ! Map token to all upper case (make case-insensitive):
     call inupper(cscolon)
     call inupper(csplus)
     call inupper(cstimes)
     call inupper(cs1colon)
     call inupper(cs1plus)
     call inupper(cs1times)
     call inupper(cs2colon)
     call inupper(cs2plus)
     call inupper(cs2times)

     ! Absolute index of tokens in string:
     itoken_colon=index(string,cscolon(1:cslen))
     itoken_plus=index(string,csplus(1:cslen))
     itoken_times=index(string,cstimes(1:cslen))
     itoken_1colon=index(string,cs1colon(1:cs1len))
     itoken_1plus=index(string,cs1plus(1:cs1len))
     itoken_1times=index(string,cs1times(1:cs1len))
     itoken_2colon=index(string,cs2colon(1:cs1len))
     itoken_2plus=index(string,cs2plus(1:cs1len))
     itoken_2times=index(string,cs2times(1:cs1len))

     ! Look for another occurence of the same tokens in string
     itoken2_colon=index(string,cscolon(1:cslen), BACK=.true. )
     itoken2_plus=index(string,csplus(1:cslen), BACK=.true. )
     itoken2_times=index(string,cstimes(1:cslen), BACK=.true. )
     itoken2_1colon=index(string,cs1colon(1:cs1len), BACK=.true. )
     itoken2_1plus=index(string,cs1plus(1:cs1len), BACK=.true. )
     itoken2_1times=index(string,cs1times(1:cs1len), BACK=.true. )
     itoken2_2colon=index(string,cs2colon(1:cs1len), BACK=.true. )
     itoken2_2plus=index(string,cs2plus(1:cs1len), BACK=.true. )
     itoken2_2times=index(string,cs2times(1:cs1len), BACK=.true. )

     if(jdtset==0)then

       ! If the multi-dataset mode is not used, no token should have been found
       if(itoken_colon+itoken_plus+itoken_times+ itoken_2colon+itoken_2plus+itoken_2times > 0 ) then
         write(msg,'(a,a,a,a,a,a,a,a,a,a,a,a, a)' )&
         'Although the multi-dataset mode is not activated,',ch10,&
         'the keyword "',trim(cs),'" has been found',ch10,&
         'appended with  + * or :  .',ch10,&
         'This is not allowed.',ch10,&
         'Action: remove the appended keyword, or',ch10,&
         'use the multi-dataset mode (ndtset/=0).'
         MSG_ERROR(msg)
       end if
       if(itoken_1colon+itoken_1plus+itoken_1times > 0 ) then
         write(msg, '(a,a,a,a,a,a,a,a,a,a,a,a,a)' )&
         'Although the multi-dataset mode is not activated,',ch10,&
         'the keyword "',trim(cs),'" has been found',ch10,&
         'appended with ? , then + * or :  .',ch10,&
         'This is not allowed.',ch10,&
         'Action: remove the appended keyword, or',ch10,&
         'use the multi-dataset mode (ndtset/=0).'
         MSG_ERROR(msg)
       end if

     else

       ! If the multi-dataset mode is used, exactly zero or two token must be found
       sum_token=0
       if(itoken_colon/=0)sum_token=sum_token+1
       if(itoken_plus /=0)sum_token=sum_token+1
       if(itoken_times/=0)sum_token=sum_token+1
       if(itoken_1colon/=0)sum_token=sum_token+1
       if(itoken_1plus /=0)sum_token=sum_token+1
       if(itoken_1times/=0)sum_token=sum_token+1
       if(itoken_2colon/=0)sum_token=sum_token+1
       if(itoken_2plus /=0)sum_token=sum_token+1
       if(itoken_2times/=0)sum_token=sum_token+1

       if(sum_token/=0 .and. sum_token/=2) then
         write(msg, '(a,a,a,a,a,i0,a,a,a,a,a,a,a)' )&
         'The keyword "',trim(cs),'" has been found to take part',ch10,&
         'to series definition in the multi-dataset mode  ',sum_token,' times.',ch10,&
         'This is not allowed, since it should be used once with ":",',ch10,&
         'and once with "+" or "*".',ch10,&
         'Action: change the number of occurences of this keyword.'
         MSG_ERROR(msg)
       end if

       ! If the multi-dataset mode is used, make sure that no twice the same combined keyword happens
       ier=0
       if(itoken_colon/=itoken2_colon)then
         ier=1 ; cs=cscolon
       end if
       if(itoken_plus/=itoken2_plus)then
         ier=1 ; cs=csplus
       end if
       if(itoken_times/=itoken2_times)then
         ier=1 ; cs=cstimes
       end if
       if(itoken_1colon/=itoken2_1colon)then
         ier=1 ; cs=cs1colon
       end if
       if(itoken_1plus/=itoken2_1plus)then
         ier=1 ; cs=cs1plus
       end if
       if(itoken_1times/=itoken2_1times)then
         ier=1 ; cs=cs1times
       end if
       if(itoken_2colon/=itoken2_2colon)then
         ier=1 ; cs=cs2colon
       end if
       if(itoken_2plus/=itoken2_2plus)then
         ier=1 ; cs=cs2plus
       end if
       if(itoken_2times/=itoken2_2times)then
         ier=1 ; cs=cs2times
       end if
       if(ier==1)then
         write(msg, '(a,a,a,a,a,a,a)' )&
         'There are two occurences of the keyword "',cs(1:cslen),'" in the input file.',ch10,&
         'This is confusing, so it has been forbidden.',ch10,&
         'Action: remove one of the two occurences.'
         MSG_ERROR(msg)
       end if

       ! Select the series according to the presence of a colon flag
       if(itoken_colon>0)then
         opttoken=2
         ds_input_=jdtset
       else if(itoken_1colon>0)then
         opttoken=3
         cscolon=cs1colon ; csplus=cs1plus ; cstimes=cs1times
         itoken_colon=itoken_1colon
         itoken_plus=itoken_1plus ; itoken_times=itoken_1times
         cslen=cs1len
         ds_input_=jdtset
       else if(itoken_2colon>0)then
         opttoken=4
         cscolon=cs2colon ; csplus=cs2plus ; cstimes=cs2times
         itoken_colon=itoken_2colon
         itoken_plus=itoken_2plus ; itoken_times=itoken_2times
         cslen=cs1len
         ds_input_=jdtset
       end if

       ! Make sure that the proper combination of : + and * is found .
       if(itoken_colon > 0 .and. (itoken_plus==0 .and. itoken_times==0) )then
         write(msg, '(13a)' )&
         'The keyword "',cscolon(1:cslen),'" initiate a series,',ch10,&
         'but there is no occurence of "',csplus(1:cslen),'" or "',cstimes(1:cslen),'".',ch10,&
         'Action: either suppress the series, or make the increment',ch10,&
         'or the factor available.'
         MSG_ERROR(msg)
       end if
       if(itoken_plus/=0 .and. itoken_times/=0)then
         write(msg, '(a,a, a,a,a,a,a)' )&
         'The combined occurence of keywords "',csplus(1:cslen),'" and "',cstimes(1:cslen),'" is not allowed.',ch10,&
         'Action: suppress one of them in your input file.'
         MSG_ERROR(msg)
       end if
       if(itoken_colon==0 .and. (itoken_plus/=0 .or. itoken_times/=0) ) then
         cs=csplus
         if(itoken_times/=0)cs=cstimes
         write(msg, '(a,a,a,a,a,a,a,a,a,a,a)' )&
         'The keyword "',cscolon(1:cslen),'" does not appear in the input file.',ch10,&
         'However, the keyword "',cs(1:cslen),'" appears.',ch10,&
         'This is forbidden.',ch10,&
         'Action: make the first appear, or suppress the second.'
         MSG_ERROR(msg)
       end if

       ! At this stage, either
       !    - itoken_colon vanish as well as itoken_plus and itoken_times
       !    - itoken_colon does not vanish,
       ! as well as one of itoken_plus or itoken_times

     end if ! End the condition of multi-dataset mode
   end if ! End the check on existence of a series

   ! --------------------------------------------------------------------------
   ! (3) if not found, try to find the token with non-modified string
   if (opttoken==0) then

     cs=blank//token(1:toklen)//blank
     cslen=toklen+2

     ! Map token to all upper case (make case-insensitive):
     call inupper(cs)

     ! Absolute index of blank//token//blank in string:
     itoken=index(string,cs(1:cslen))

     ! Look for another occurence of the same token in string, if so, leaves:
     itoken2=index(string,cs(1:cslen), BACK=.true. )
     if (itoken/=itoken2) then
       write(msg, '(a,a,a,a,a,a,a)' )&
       'There are two occurences of the keyword "',cs(1:cslen),'" in the input file.',ch10,&
       'This is confusing, so it has been forbidden.',ch10,&
       'Action: remove one of the two occurences.'
       MSG_ERROR(msg)
     end if

     if(itoken/=0) then
       opttoken=1
       ds_input_=0
     end if

   end if

   ! --------------------------------------------------------------------------
   ! If jdtset==0, means that the multi-dataset mode is not used, so
   ! checks whether the input file contains a multi-dataset keyword,
   ! and if this occurs, stop. Check also the forbidden occurence of
   ! use of 0 as a multi-dataset index.
   ! Note that the occurence of series initiators has already been checked.

   do trial_jdtset=0,9
     if(jdtset==0 .or. trial_jdtset==0)then
       write(appen,'(i1)')trial_jdtset
       trial_cs=blank//token(1:toklen)//trim(appen)
       trial_cslen=toklen+2
       ! Map token to all upper case (make case-insensitive):
       call inupper(trial_cs)
       ! Look for an occurence of this token in string, if so, leaves:
       itoken2=index(string,trial_cs(1:trial_cslen))
       if(itoken2/=0)then
         if(trial_jdtset==0)then
           write(msg, '(a,a,a,a,a,a,a)' )&
           'There is an occurence of the keyword "',trim(token),'" appended with 0 in the input file.',ch10,&
           'This is forbidden.',ch10,&
           'Action: remove this occurence.'
         else
           write(msg, '(a,a,a,a,a,i1,a,a,a,a,a)' )&
           'In the input file, there is an occurence of the ',ch10,&
           'keyword "',trim(token),'", appended with the digit "',trial_jdtset,'".',ch10,&
           'This is forbidden when ndtset==0 .',ch10,&
           'Action: remove this occurence, or change ndtset.'
         end if
         MSG_ERROR(msg)
       end if
     end if
   end do

 end if

 !===========================================================================
 ! At this stage, the location of the keyword string is known, as well
 ! as its length. So, can read the data.
 ! Usual reading if opttoken==1 (need itoken).
 ! If opttoken>=2, the characteristics of a series must be read
 ! (need itoken_colon and either itoken_plus or itoken_times)

 tread = 0
 typevar='INT'

 if(typevarphys=='LOG')typevar='INT'
 if(typevarphys=='DPR' .or. typevarphys=='LEN' .or. typevarphys=='ENE' .or. &
    typevarphys=='BFI' .or. typevarphys=='TIM') typevar='DPR'

 if (typevarphys=='KEY') then
   ! Consistency check for keyword (no multidataset, no series)
   if (opttoken>=2) then
     write(msg, '(9a)' )&
     'For the keyword "',cs(1:cslen),'", of KEY type,',ch10,&
     'a series has been defined in the input file.',ch10,&
     'This is forbidden.',ch10,'Action: check your input file.'
     MSG_ERROR(msg)
   end if
   if (narr>=2) then
     write(msg, '(9a)' )&
     'For the keyword "',cs(1:cslen),'", of KEY type,',ch10,&
     'the number of data requested is larger than 1.',ch10,&
     'This is forbidden.',ch10,'Action: check your input file.'
     MSG_ERROR(msg)
   end if
 end if

 ! There is something to be read if opttoken>=1
 if (opttoken==1) then

   ! write(std_out,*)' intagm : opttoken==1 , token has been found, will read '
   ! Absolute location in string of blank which follows token:
   b1 = itoken + cslen - 1

   if (typevarphys == 'KEY') then
     ! In case of typevarphys='KEY', the chain of character will be returned in cs.
     ABI_CHECK(present(key_value), "typevarphys == KEY requires optional argument key_value")
     b2 = index(string(b1+1:), '"')
     ABI_CHECK(b2 /= 0, sjoin('Cannot find first " defining string for token:', token))
     b2 = b1 + b2 + 1
     b3 = index(string(b2:), '"')
     ABI_CHECK(b3 /= 0, sjoin('Cannot find second " defining string for token:', token))
     b3 = b3 + b2 - 2
     if ((b3 - b2 + 1) > len(key_value)) then
       MSG_ERROR("Len of key_value too small to contain value parsed from file")
     end if
     key_value = adjustl(string(b2:b3))

   else
     ! Read the array (or eventual scalar) that follows the blank
     call inarray(b1,cs,dprarr,intarr,marr,narr,string,typevarphys)
   end if

   ! if this point is reached then data has been read in successfully
   tread = 1

 else if(opttoken>=2) then

   ! write(std_out,*)' intagm : opttoken>=2 , token has been found, will read '
   ABI_ALLOCATE(dpr1,(narr))
   ABI_ALLOCATE(dpr2,(narr))
   ABI_ALLOCATE(int1,(narr))
   ABI_ALLOCATE(int2,(narr))

   ! Absolute location in string of blank which follows token//':':
   b1=itoken_colon+cslen-1
   call inarray(b1,cscolon,dpr1,int1,narr,narr,string,typevarphys)

   ! Initialise number even if the if series treat all cases.
   number=1
   ! Define the number of the term in the series
   if(opttoken==2)number=jdtset-1
   if(opttoken==3)number=unities-1
   if(opttoken==4)number=dozens-1

   ! Distinguish additive and multiplicative series
   if(itoken_plus/=0)then

     b1=itoken_plus+cslen-1
     call inarray(b1,csplus,dpr2,int2,narr,narr,string,typevarphys)

     if(typevar=='INT')then
       intarr(1:narr)=int1(:)+int2(:)*number
     else if(typevar=='DPR')then
       dprarr(1:narr)=dpr1(:)+dpr2(:)*number
     end if

   else if(itoken_times/=0)then

     b1=itoken_times+cslen-1
     call inarray(b1,cstimes,dpr2,int2,narr,narr,string,typevarphys)
     if(typevar=='INT')then
       intarr(1:narr)=int1(:)*int2(:)**number
     else if(typevar=='DPR')then
       dprarr(1:narr)=dpr1(:)*dpr2(:)**number
     end if

   end if

   tread = 1

   ABI_FREE(dpr1)
   ABI_FREE(dpr2)
   ABI_FREE(int1)
   ABI_FREE(int2)
 end if

 if(present(ds_input)) ds_input = ds_input_

 !write(std_out,*) ' intagm : exit value tread=',tread
 !write(std_out,*) ' intarr =',intarr(1:narr)
 !write(std_out,*) ' dprarr =',dprarr(1:narr)

end subroutine intagm
!!***

!----------------------------------------------------------------------

!!****f* m_parser/ingeo_img_1D
!! NAME
!!  intagm_img_1D
!!
!! FUNCTION
!!  Read input file variables according to images path definition (1D array)
!!
!!  This function is exposed through generic interface that allows to
!!  initialize some of the geometry variables in the case of "images".
!!  Set up: acell, scalecart, rprim, angdeg, xred, xcart, vel
!!  These variables can be defined for a set of images of the cell.
!!  They also can be be defined along a path (in the configuration space).
!!  The path must be defined with its first and last points, but also
!!  with intermediate points.
!!
!! INPUTS
!!  iimage=index of the current image
!!  jdtset=number of the dataset looked for
!!  lenstr=actual length of the input string
!!  nimage=number of images
!!  size1,size2, ...: size of array to be read (dp_data)
!!  string=character string containing 'tags' and data.
!!  token=character string for tagging the data to be read in input string
!!  typevarphys= variable type (for dimensionality purposes)
!!
!! SIDE EFFECTS
!!  dp_data(size1,size2,...)=data to be read (double precision)
!!  tread_ok=flag to be set to 1 if the data have been found in input string
!!
!! NOTES
!! The routine is a generic interface calling subroutine according to the
!! number of arguments of the variable to be read
!!
!! PARENTS
!!
!! CHILDREN
!!      intagm
!!
!! SOURCE

subroutine intagm_img_1D(dp_data,iimage,jdtset,lenstr,nimage,size1,string,token,tread_ok,typevarphys)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iimage,jdtset,lenstr,nimage,size1
 integer,intent(inout) :: tread_ok
 real(dp),intent(inout) :: dp_data(size1)
 character(len=*),intent(in) :: typevarphys
 character(len=*),intent(in) :: token
 character(len=*),intent(in) :: string
!arrays

!Local variables-------------------------------
!scalars
 integer :: iimage_after,iimage_before,marr,tread_after,tread_before,tread_current
 real(dp) :: alpha
 character(len=10) :: stringimage
 character(len=3*len(token)+10) :: token_img
!arrays
 integer, allocatable :: intarr(:)
 real(dp),allocatable :: dprarr(:),dp_data_after(:),dp_data_before(:)

! *************************************************************************

!Nothing to do in case of a single image
 if (nimage<=1) return

 marr=size1
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

!First, try to read data for current image
 tread_current=0
 write(stringimage,'(i10)') iimage
 token_img=trim(token)//'_'//trim(adjustl(stringimage))//'img'
 call intagm(dprarr,intarr,jdtset,marr,size1,string(1:lenstr),&
&            token_img,tread_current,typevarphys)
 if (tread_current==1)then
   dp_data(1:size1)=dprarr(1:size1)
   tread_ok=1
 end if
 if (tread_current==0.and.iimage==nimage) then
!  If the image is the last one, try to read data for last image (_lastimg)
   token_img=trim(token)//'_lastimg'
   call intagm(dprarr,intarr,jdtset,marr,size1,string(1:lenstr),&
&              token_img,tread_current,typevarphys)
   if (tread_current==1)then
     dp_data(1:size1)=dprarr(1:size1)
     tread_ok=1
   end if
 end if

 if (tread_current==0) then

!  The current image is not directly defined in the input string
   ABI_ALLOCATE(dp_data_before,(size1))
   ABI_ALLOCATE(dp_data_after,(size1))

!  Find the nearest previous defined image
   tread_before=0;iimage_before=iimage
   do while (iimage_before>1.and.tread_before/=1)
     iimage_before=iimage_before-1
     write(stringimage,'(i10)') iimage_before
     token_img=trim(token)//'_'//trim(adjustl(stringimage))//'img'
     call intagm(dprarr,intarr,jdtset,marr,size1,string(1:lenstr),&
&                token_img,tread_before,typevarphys)
     if (tread_before==1) dp_data_before(1:size1)=dprarr(1:size1)
   end do
   if (tread_before==0) then
     iimage_before=1
     dp_data_before(1:size1)=dp_data(1:size1)
   end if

!  Find the nearest following defined image
   tread_after=0;iimage_after=iimage
   do while (iimage_after<nimage.and.tread_after/=1)
     iimage_after=iimage_after+1
     write(stringimage,'(i10)') iimage_after
     token_img=trim(token)//'_'//trim(adjustl(stringimage))//'img'
     call intagm(dprarr,intarr,jdtset,marr,size1,string(1:lenstr),&
&                token_img,tread_after,typevarphys)
     if (tread_after==1) dp_data_after(1:size1)=dprarr(1:size1)
     if (tread_after==0.and.iimage_after==nimage) then
       token_img=trim(token)//'_lastimg'
       call intagm(dprarr,intarr,jdtset,marr,size1,string(1:lenstr),&
&                  token_img,tread_after,typevarphys)
       if (tread_after==1) dp_data_after(1:size1)=dprarr(1:size1)
     end if
   end do
   if (tread_after==0) then
     iimage_after=nimage
     dp_data_after(1:size1)=dp_data(1:size1)
   end if

!  Interpolate image data
   if (tread_before==1.or.tread_after==1) then
     alpha=real(iimage-iimage_before,dp)/real(iimage_after-iimage_before,dp)
     dp_data(1:size1)=dp_data_before(1:size1) &
&                    +alpha*(dp_data_after(1:size1)-dp_data_before(1:size1))
     tread_ok=1
   end if

   ABI_DEALLOCATE(dp_data_before)
   ABI_DEALLOCATE(dp_data_after)

 end if

 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(dprarr)

end subroutine intagm_img_1D
!!***

!----------------------------------------------------------------------

!!****f* m_parser/ingeo_img_2D
!! NAME
!!  intagm_img_2D
!!
!! FUNCTION
!!  Read input file variables according to images path definition (2D array)
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!      intagm
!!
!! SOURCE

subroutine intagm_img_2D(dp_data,iimage,jdtset,lenstr,nimage,size1,size2,string,token,tread_ok,typevarphys)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iimage,jdtset,lenstr,nimage,size1,size2
 integer,intent(inout) :: tread_ok
 real(dp),intent(inout) :: dp_data(size1,size2)
 character(len=*),intent(in) :: typevarphys
 character(len=*),intent(in) :: token
 character(len=*),intent(in) :: string
!arrays

!Local variables-------------------------------
!scalars
 integer :: iimage_after,iimage_before,marr,tread_after,tread_before,tread_current
 real(dp) :: alpha
 character(len=10) :: stringimage
 character(len=3*len(token)+10) :: token_img
!arrays
 integer, allocatable :: intarr(:)
 real(dp),allocatable :: dprarr(:),dp_data_after(:,:),dp_data_before(:,:)

! *************************************************************************

!Nothing to do in case of a single image
 if (nimage<=1) return

 marr=size1*size2
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

!First, try to read data for current image
 tread_current=0
 write(stringimage,'(i10)') iimage
 token_img=trim(token)//'_'//trim(adjustl(stringimage))//'img'
 call intagm(dprarr,intarr,jdtset,marr,size1*size2,string(1:lenstr),&
&            token_img,tread_current,typevarphys)
 if (tread_current==1)then
   dp_data(1:size1,1:size2)=reshape( dprarr(1:size1*size2),(/size1,size2/) )
   tread_ok=1
 end if
 if (tread_current==0.and.iimage==nimage) then
!  In the image is the last one, try to read data for last image (_lastimg)
   token_img=trim(token)//'_lastimg'
   call intagm(dprarr,intarr,jdtset,marr,size1*size2,string(1:lenstr),&
&              token_img,tread_current,typevarphys)
   if (tread_current==1)then
     dp_data(1:size1,1:size2)=reshape( dprarr(1:size1*size2),(/size1,size2/) )
     tread_ok=1
   end if
 end if

 if (tread_current==0) then

!  The current image is not directly defined in the input string
   ABI_ALLOCATE(dp_data_before,(size1,size2))
   ABI_ALLOCATE(dp_data_after,(size1,size2))

!  Find the nearest previous defined image
   tread_before=0;iimage_before=iimage
   do while (iimage_before>1.and.tread_before/=1)
     iimage_before=iimage_before-1
     write(stringimage,'(i10)') iimage_before
     token_img=trim(token)//'_'//trim(adjustl(stringimage))//'img'
     call intagm(dprarr,intarr,jdtset,marr,size1*size2,string(1:lenstr),&
&                token_img,tread_before,typevarphys)
     if (tread_before==1) &
&      dp_data_before(1:size1,1:size2)=reshape( dprarr(1:size1*size2),(/size1,size2/) )
   end do
   if (tread_before==0) then
     iimage_before=1
     dp_data_before(1:size1,1:size2)=dp_data(1:size1,1:size2)
   end if

!  Find the nearest following defined image
   tread_after=0;iimage_after=iimage
   do while (iimage_after<nimage.and.tread_after/=1)
     iimage_after=iimage_after+1
     write(stringimage,'(i10)') iimage_after
     token_img=trim(token)//'_'//trim(adjustl(stringimage))//'img'
     call intagm(dprarr,intarr,jdtset,marr,size1*size2,string(1:lenstr),&
&                token_img,tread_after,typevarphys)
     if (tread_after==1) &
&      dp_data_after(1:size1,1:size2)=reshape( dprarr(1:size1*size2),(/size1,size2/) )
     if (tread_after==0.and.iimage_after==nimage) then
       token_img=trim(token)//'_lastimg'
       call intagm(dprarr,intarr,jdtset,marr,size1*size2,string(1:lenstr),&
&                  token_img,tread_after,typevarphys)
       if (tread_after==1) &
&        dp_data_after(1:size1,1:size2)=reshape( dprarr(1:size1*size2),(/size1,size2/) )
     end if
   end do
   if (tread_after==0) then
     iimage_after=nimage
     dp_data_after(1:size1,1:size2)=dp_data(1:size1,1:size2)
   end if

!  Interpolate image data
   if (tread_before==1.or.tread_after==1) then
     alpha=real(iimage-iimage_before,dp)/real(iimage_after-iimage_before,dp)
     dp_data(1:size1,1:size2)=dp_data_before(1:size1,1:size2) &
&       +alpha*(dp_data_after(1:size1,1:size2)-dp_data_before(1:size1,1:size2))
     tread_ok=1
   end if

   ABI_DEALLOCATE(dp_data_before)
   ABI_DEALLOCATE(dp_data_after)

 end if

 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(dprarr)

end subroutine intagm_img_2D
!!***

!!****f* m_parser/inarray
!! NAME
!! inarray
!!
!! FUNCTION
!! Read the array of narr numbers located immediately after a specified blank in a string of character.
!! Might read instead one word, after the specified blank. Takes care of multipliers.
!!
!! INPUTS
!!  cs=character token
!!  marr=dimension of the intarr and dprarr arrays, as declared in the
!!   calling subroutine.
!!  narr=actual size of array to be read in  (if typevarphys='KEY', only narr=1 is allowed)
!!  string=character string containing the data.
!!  typevarphys=variable type (might indicate the physical meaning of
!!   for dimensionality purposes)
!!   'INT' => integer
!!   'DPR' => real(dp) (no special treatment)
!!   'LEN' => real(dp) (expect a "length", identify bohr, au, nm or angstrom,
!!            and return in au -atomic units=bohr- )
!!   'ENE' => real(dp) (expect a "energy", identify Ha, hartree, eV, Ry, Rydberg)
!!   'BFI' => real(dp) (expect a "magnetic field", identify T, Tesla)
!!   'TIM' => real(dp) (expect a "time", identify S, Second)
!!   'LOG' => integer, but read logical variable T,F,.true., or .false.
!!
!! OUTPUT
!!  intarr(1:narr), dprarr(1:narr)
!!   integer or real(dp) arrays, respectively into which data is read. Use these arrays even for scalars.
!!  errcod: if /= 0, then something went wrong in subroutine "inread"
!!
!! SIDE EFFECT
!!   b1=absolute location in string of blank which follows the token (will be modified in the execution)
!!
!! PARENTS
!!      intagm
!!
!! CHILDREN
!!      inread,wrtout
!!
!! SOURCE

subroutine inarray(b1,cs,dprarr,intarr,marr,narr,string,typevarphys)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: marr,narr
 integer,intent(inout) :: b1
 character(len=*),intent(in) :: string
 character(len=*),intent(in) :: typevarphys
 character(len=*),intent(in) :: cs
!arrays
 integer,intent(inout) :: intarr(marr) !vz_i
 real(dp),intent(out) :: dprarr(marr)

!Local variables-------------------------------
 character(len=1), parameter :: blank=' '
!scalars
 integer :: asciichar,b2,errcod,ii,integ,istar,nrep,strln
 real(dp) :: factor,real8
 character(len=3) :: typevar
 character(len=500*4) :: msg

! *************************************************************************

 !write(std_out,'(2a)' )' inarray: token: ',trim(cs)
 !write(std_out,'(2a)' )'          string: ',trim(string(b1:))
 !write(std_out,'(a,i0)' )'        narr: ',narr
 !write(std_out,'(2a)' )'          typevarphys: ',typevarphys

 ii = 0
 typevar='INT'
 if(typevarphys=='LOG') typevar='INT'
 if(typevarphys=='DPR' .or. typevarphys=='LEN' .or. typevarphys=='ENE' .or. &
    typevarphys=='BFI' .or. typevarphys=='TIM') typevar='DPR'

 strln=len_trim(string)

 do while (ii < narr)

   ! Relative location of next blank after data
   ! b1 is the last character of the string
   if (b1>=strln) exit

   b2 = index(string(b1+1:),blank)
   ! If no second blank is found put the second blank just beyond strln
   if(b2==0) b2=strln-b1+1

   ! nrep tells how many times to repeat input in array:
   nrep=1

   ! Check for *, meaning repeated input (as in list-directed input):
   istar=index(string(b1+1:b1+b2-1),'*')
   if (istar/=0) then
     if (istar==1) then ! Simply fills the array with the data, repeated as many times as needed
       nrep=narr-ii
       errcod=0
     else
       call inread(string(b1+1:b1+istar-1),istar-1,'INT',nrep,real8,errcod)
     end if
     if (errcod/=0) exit
     ! Shift starting position of input field:
     b1=b1+istar
     b2=b2-istar
   end if

   ! Read data internally by calling inread at entry ini:
   call inread(string(b1+1:b1+b2-1),b2-1,typevarphys,integ,real8,errcod)
   if (errcod/=0) exit

   ! Allow for list-directed input with repeat number nrep:
   if(typevar=='INT')then
     intarr(1+ii:min(nrep+ii,narr))=integ
   else if(typevar=='DPR')then
     dprarr(1+ii:min(nrep+ii,narr))=real8
   else
     MSG_BUG('Disallowed typevar: '//typevar)
   end if
   ii=min(ii+nrep,narr)

   !  Find new absolute location of next element of array:
   b1=b1+b2

 end do ! while (ii<narr). Note "exit" instructions within loop.

 if (errcod /= 0) then
   write(msg, '(5a,i0,12a)' ) &
   'An error occurred reading data for keyword `',trim(cs),'`,',ch10,&
   'looking for ',narr,' elements.', ch10, &
   'There is a problem with the input string:',ch10,trim(string(b1:)), ch10, &
   'Maybe a disagreement between the declared dimension of the array,',ch10,&
   'and the number of items provided. ',ch10,&
   'Action: correct your input file and especially the keyword: ', trim(cs)
   MSG_ERROR(msg)
 end if

 ! In case of 'LEN', 'ENE', 'BFI', or 'TIM', try to identify the unit
 if (typevarphys=='LEN' .or. typevarphys=='ENE' .or. typevarphys=='BFI' .or. typevarphys=='TIM') then
   do
     ! Relative location of next blank after data
     if(b1>=strln)exit   ! b1 is the last character of the string
     b2=index(string(b1+1:),blank)
     ! If no second blank is found put the second blank just beyond strln
     if(b2==0) b2=strln-b1+1

     ! write(std_out,*)' inarray : strln=',strln
     ! write(std_out,*)' inarray : b1=',b1, b2=',b2
     ! write(std_out,*)' inarray : string(b1+1:)=',string(b1+1:)
     ! write(std_out,*)' typevarphys==',typevarphys

     ! Identify the presence of a non-digit character
     asciichar=iachar(string(b1+1:b1+1))
     if(asciichar<48 .or. asciichar>57)then
       factor=one
       if(typevarphys=='LEN' .and. b2>=3)then
         if(string(b1+1:b1+6)=='ANGSTR')then
           factor=one/Bohr_Ang
         else if(string(b1+1:b1+3)=='NM ')then
           factor=ten/Bohr_Ang
         end if
       else if(typevarphys=='ENE' .and. b2>=3)then
         if(string(b1+1:b1+3)=='RY ')then
           factor=half
         else if(string(b1+1:b1+3)=='EV ')then
           factor=one/Ha_eV
         end if
       else if(typevarphys=='ENE' .and. b2>=2)then
         if(string(b1+1:b1+2)=='K ') factor=kb_HaK
       else if(typevarphys=='BFI' .and. b2>=2)then
         if(string(b1+1:b1+2)=='T ' .or. string(b1+1:b1+2)=='TE') factor=BField_Tesla
       else if (typevarphys=='TIM' .and. b2>=2) then
         if( string(b1+1:b1+2)=='SE' .or. string(b1+1:b1+2)=='S ') factor=one/Time_Sec
       endif

       dprarr(1:narr)=dprarr(1:narr)*factor
       exit
     else
       ! A digit has been observed, go to the next sequence
       b1=b2
       cycle
     end if

   end do
 end if

!write(std_out,*)' dprarr(1:narr)==',dprarr(1:narr)
!write(std_out,*)' inarray : exit '

end subroutine inarray
!!***

!!****f* m_parser/importxyz
!! NAME
!! importxyz
!!
!! FUNCTION
!! Examine the input string, to see whether data from xyz
!! file(s) has to be incorporated.
!! For each such xyz file, translate the relevant
!! information into intermediate input variables compatible
!! with the usual ABINIT formatting, then append it
!! to the input string.
!!
!! INPUTS
!!  string_raw*(strln)=raw string of character from input file (with original case)
!!  strln=maximal number of character of string, as declared in the calling routine
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  lenstr=actual number of character in string
!!  string_upper*(strln)=string of character
!!   the string (with upper case) from the input file, to which the xyz data are appended to it
!!
!! PARENTS
!!      m_ab7_invars_f90,parsefile
!!
!! CHILDREN
!!      append_xyz,incomprs,wrtout
!!
!! SOURCE

subroutine importxyz(lenstr,string_raw,string_upper,strln)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: strln
 integer,intent(inout) :: lenstr
 character(len=*),intent(in) :: string_raw
 character(len=*),intent(inout) :: string_upper

!Local variables-------------------------------
 character :: blank=' '
!scalars
 integer :: dtset_len,ixyz,ii,index_already_done,index_xyz_fname
 integer :: index_xyz_fname_end,index_xyz_token,kk
 character(len=2) :: dtset_char
 character(len=500) :: msg
 character(len=fnlen) :: xyz_fname

!************************************************************************

 index_already_done=1
 ixyz=0

 do
   ! Infinite do-loop, to identify the presence of the xyzFILE token
   index_xyz_token=index(string_upper(index_already_done:lenstr),"XYZFILE")
   if(index_xyz_token==0)exit

   ixyz=ixyz+1
   if(ixyz==1)then
     write(msg,'(80a)')('=',ii=1,80)
     call wrtout(ab_out,msg)
   end if

   ! The xyzFILE token has been identified
   index_xyz_token=index_already_done+index_xyz_token-1

   ! Find the related dataset tag, and length
   dtset_char=string_upper(index_xyz_token+7:index_xyz_token+8)
   if(dtset_char(1:1)==blank)dtset_char(2:2)=blank
   dtset_len=len_trim(dtset_char)

   ! Find the name of the xyz file
   index_xyz_fname=index_xyz_token+8+dtset_len
   index_xyz_fname_end=index(string_upper(index_xyz_fname:lenstr),blank)

   if(index_xyz_fname_end ==0 )then
     write(msg, '(5a,i4,2a)' )&
     'Could not find the name of the xyz file.',ch10,&
     'index_xyz_fname_end should be non-zero, while it is :',ch10,&
     'index_xyz_fname_end=',index_xyz_fname_end,ch10,&
     'Action: check the filename that was provided after the XYZFILE input variable keyword.'
     MSG_ERROR(msg)
   end if

   index_xyz_fname_end=index_xyz_fname_end+index_xyz_fname-1

   index_already_done=index_xyz_fname_end

   ! Initialize xyz_fname to a blank line
   xyz_fname=repeat(blank,fnlen)
   xyz_fname=string_raw(index_xyz_fname:index_xyz_fname_end-1)

   write(msg, '(3a)') ch10, ' importxyz : Identified token XYZFILE, referring to file ',trim(xyz_fname)
   call wrtout([std_out, ab_out],msg)

   ! Append the data from the xyz file to the string, and update the length of the string
   call append_xyz(dtset_char,lenstr,string_upper,xyz_fname,strln)

   ! erase the file name from string_upper
   string_upper(index_xyz_fname:index_xyz_fname_end-1) = blank
 end do

 if (index_already_done > 1) then
   ! Initialize xyz_fname to a blank line
   xyz_fname=repeat(blank,fnlen)
   call append_xyz("-1",lenstr,string_upper,xyz_fname,strln)
 end if

 if(ixyz/=0)then
   call incomprs(string_upper,lenstr)
   ! A blank is needed at the beginning of the string
   do kk=lenstr,1,-1
     string_upper(kk+1:kk+1)=string_upper(kk:kk)
   end do
   string_upper(1:1)=blank
   lenstr=lenstr+1
   write(msg,'(a,80a,a)')ch10,('=',ii=1,80),ch10
   call wrtout(ab_out,msg)
 end if

end subroutine importxyz
!!***

!!****f* m_parser/append_xyz
!! NAME
!! append_xyz
!!
!! FUNCTION
!! Translate the data from a xyz file (xyz_fname),
!! and add it at the end of the usual ABINIT input data string (string),
!! taking into account the dtset (dtset_char)
!!
!! INPUTS
!!  dtset_char*2=possible dtset label
!!  xyz_fname = name of the xyz file
!!  strln=maximal number of characters of string, as declared in the calling routine
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  lenstr=actual number of characters in string
!!  string*(strln)=string of characters  (upper case) to which the xyz data are appended
!!
!! PARENTS
!!      importxyz
!!
!! CHILDREN
!!      atomdata_from_symbol,wrtout
!!
!! SOURCE

subroutine append_xyz(dtset_char,lenstr,string,xyz_fname,strln)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: strln
 integer,intent(inout) :: lenstr
 character(len=2),intent(in) :: dtset_char
 character(len=fnlen),intent(in) :: xyz_fname
 character(len=strln),intent(inout) :: string

!Local variables-------------------------------
 character :: blank=' '
!scalars
 integer :: unitxyz, iatom, natom, mu
 integer :: lenstr_new
 integer :: lenstr_old
 integer :: ntypat
 real(dp) :: znucl
 character(len=5) :: string5
 character(len=20) :: string20
 character(len=500) :: msg
 type(atomdata_t) :: atom
!arrays
 real(dp),allocatable :: xcart(:,:)
 integer, save :: atomspecies(200) = 0
 character(len=500), save :: znuclstring = ""
 character(len=2),allocatable :: elementtype(:)

!************************************************************************

 lenstr_new=lenstr

 if (dtset_char == "-1") then
   ! write znucl
   lenstr_old=lenstr_new
   lenstr_new=lenstr_new+7+len_trim(znuclstring)+1
   string(lenstr_old+1:lenstr_new)=" ZNUCL"//blank//trim(znuclstring)//blank

   ! write ntypat
   ntypat = sum(atomspecies)
   write(string20,'(i10)') ntypat
   lenstr_old=lenstr_new
   lenstr_new=lenstr_new+8+len_trim(string20)+1
   string(lenstr_old+1:lenstr_new)=" NTYPAT"//blank//trim(string20)//blank

   return
 end if

 ! open file with xyz data
 if (open_file(xyz_fname, msg, newunit=unitxyz, status="unknown") /= 0) then
   MSG_ERROR(msg)
 end if
 write(msg, '(3a)')' importxyz : Opened file ',trim(xyz_fname),'; content stored in string_xyz'
 call wrtout(std_out,msg)

 ! check number of atoms is correct
 read(unitxyz,*) natom

 write(string5,'(i5)')natom
 lenstr_old=lenstr_new
 lenstr_new=lenstr_new+7+len_trim(dtset_char)+1+5
 string(lenstr_old+1:lenstr_new)=" _NATOM"//trim(dtset_char)//blank//string5

 ABI_ALLOCATE(xcart,(3,natom))
 ABI_ALLOCATE(elementtype,(natom))

 ! read dummy line
 read(unitxyz,*)

 ! read atomic types and positions
 do iatom = 1, natom
   read(unitxyz,*) elementtype(iatom), xcart(:,iatom)
   xcart(:,iatom)=xcart(:,iatom)/Bohr_Ang
   ! extract znucl for each atom type
   call atomdata_from_symbol(atom,elementtype(iatom))
   znucl = atom%znucl
   if (znucl > 200) then
     write (msg,'(5a)')&
     'found element beyond Z=200 ', ch10,&
     'Solution: increase size of atomspecies in append_xyz', ch10
     MSG_ERROR(msg)
   end if
   ! found a new atom type
   if (atomspecies(int(znucl)) == 0) then
     write(string20,'(f10.2)') znucl
     znuclstring = trim(znuclstring) // " " // trim(string20) // " "
   end if
   atomspecies(int(znucl)) = 1
 end do
 close (unitxyz)


 !Write the element types
 lenstr_old=lenstr_new
 lenstr_new=lenstr_new+7+len_trim(dtset_char)+1
 string(lenstr_old+1:lenstr_new)=" _TYPAX"//trim(dtset_char)//blank
 do iatom=1,natom
   lenstr_old=lenstr_new
   lenstr_new=lenstr_new+3
   string(lenstr_old+1:lenstr_new)=elementtype(iatom)//blank
 end do
 lenstr_old=lenstr_new
 lenstr_new=lenstr_new+3
 string(lenstr_old+1:lenstr_new)="XX " ! end card for TYPAX

 !Write the coordinates
 lenstr_old=lenstr_new
 lenstr_new=lenstr_new+8+len_trim(dtset_char)+1
 string(lenstr_old+1:lenstr_new)=" _XCART"//trim(dtset_char)//blank

 do iatom=1,natom
   do mu=1,3
     write(string20,'(f20.12)')xcart(mu,iatom)
     lenstr_old=lenstr_new
     lenstr_new=lenstr_new+20
     string(lenstr_old+1:lenstr_new)=string20
   end do
 end do

 ABI_FREE(elementtype)
 ABI_FREE(xcart)

 !Check the length of the string
 if(lenstr_new>strln)then
   write(msg,'(3a)')&
   'The maximal size of the input variable string has been exceeded.',ch10,&
   'The use of a xyz file is more character-consuming than the usual input file. Sorry.'
   MSG_BUG(msg)
 end if

 !Update the length of the string
 lenstr=lenstr_new

end subroutine append_xyz
!!***

!!****f* m_parser/chkdpr
!! NAME
!! chkdpr
!!
!! FUNCTION
!! Checks the value of an input real(dp) variable, and
!! write a sophisticated error message when it is erroneous.
!! A few conditions might have been checked before calling chkdpr,
!! and these are mentioned in the error message.
!!
!! INPUTS
!! advice_change_cond= if 1, and if an error is detected, will
!!  advice to change the value of the conditions.
!! cond_number= number of conditions checked before calling chkdpr.
!! cond_string(cond_number)= name of the variables associated to the conditions.
!! cond_values(cond_number)= value of the variables associated to the conditions.
!! input_name=name of the input variable to be checked
!! input_value=value of the input variable to be checked
!! minimal_flag=if 0, the reference_value must be matched within 1.0d-10
!!              if 1, admit values larger or equal to reference_value
!!              if -1, admit values smaller or equal to reference_value
!! reference_value=see the description of minimal_flag
!! unit=unit number for clean output file
!!
!! OUTPUT
!!  (only side effect)
!!
!! SIDE EFFECTS
!! ierr= switch it to 1 if an error was detected. No action otherwise.
!!
!! NOTES
!! cond_values(cond_number)
!! must be between -99 and 999 to be printed correctly.
!! for the time being, at most 3 conditions are allowed.
!!
!! PARENTS
!!      chkinp
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine chkdpr(advice_change_cond,cond_number,cond_string,cond_values,&
&  ierr,input_name,input_value,minimal_flag,reference_value,unit)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: advice_change_cond,cond_number,minimal_flag,unit
 integer,intent(inout) :: ierr
 real(dp),intent(in) :: input_value,reference_value
 character(len=*),intent(in) :: input_name
!arrays
 integer,intent(in) :: cond_values(4)
 character(len=*),intent(in) :: cond_string(4)

!Local variables-------------------------------
!scalars
 integer :: icond,ok
 character(len=500) :: msg

!******************************************************************

 if(cond_number<0 .or. cond_number>4)then
   write(msg,'(a,i0,a)' )'The value of cond_number is ',cond_number,'but it should be positive and < 5.'
   MSG_BUG(msg)
 end if

!Checks the allowed values
 ok=0
 if(minimal_flag==1 .and. input_value>=reference_value-tol10)      ok=1
 if(minimal_flag==-1 .and. input_value<=reference_value+tol10)     ok=1
 if(minimal_flag==0 .and. abs(input_value-reference_value)<=tol10) ok=1

 ! If there is something wrong, compose the message, and print it
 if(ok==0)then
   ierr=1
   write(msg, '(a,a)' ) ch10,' chkdpr: ERROR -'
   if(cond_number/=0)then
     do icond=1,cond_number
       ! The following format restricts cond_values(icond) to be between -99 and 999
       write(msg, '(2a,a,a,a,i4,a)' ) trim(msg),ch10,&
       '  Context : the value of the variable ',trim(cond_string(icond)),' is',cond_values(icond),'.'
     end do
   end if
   write(msg, '(2a,a,a,a,es20.12,a)' ) trim(msg),ch10,&
    '  The value of the input variable ',trim(input_name),' is',input_value,','
   if(minimal_flag==0)then
     write(msg, '(2a,a,es20.12,a)' ) trim(msg),ch10,'  while it must be equal to ',reference_value,'.'
   else if(minimal_flag==1)then
     write(msg, '(2a,a,es20.12,a)' ) trim(msg),ch10,'  while it must be larger or equal to',reference_value,'.'
   else if(minimal_flag==-1)then
     write(msg, '(2a,a,es20.12,a)' ) trim(msg),ch10,'  while it must be smaller or equal to',reference_value,'.'
   end if

   if(cond_number==0 .or. advice_change_cond==0)then
     write(msg, '(2a,a,a,a)' ) trim(msg),ch10,&
     '  Action: you should change the input variable ',trim(input_name),'.'
   else if(cond_number==1)then
     write(msg, '(2a,a,a,a,a,a)' ) trim(msg),ch10,&
     '  Action: you should change the input variables ',trim(input_name),' or ',trim(cond_string(1)),'.'
   else if(cond_number==2)then
     write(msg, '(2a,a,a,a,a,a,a,a,a,a)' ) trim(msg),ch10,&
     '  Action: you should change one of the input variables ',trim(input_name),',',ch10,&
     '   ',trim(cond_string(1)),' or ',trim(cond_string(2)),'.'
   else if(cond_number==3)then
     write(msg, '(2a,a,a,a,a,a,a,a,a,a,a,a)' ) trim(msg),ch10,&
     '  Action: you should change one of the input variables ',trim(input_name),',',ch10,&
     '   ',trim(cond_string(1)),', ',trim(cond_string(2)),' or ',trim(cond_string(3)),'.'
   end if

   call wrtout(unit,msg)
   MSG_WARNING(msg)
 end if

end subroutine chkdpr
!!***

!!****f* m_parser/chkint
!! NAME
!! chkint
!!
!! FUNCTION
!! Checks the value of an input integer variable, and
!! write a sophisticated error message when it is erroneous.
!! A few conditions might have been checked before calling chkint,
!! and these are mentioned in the error message.
!! See the examples in the NOTES
!!
!! INPUTS
!! advice_change_cond= if 1, and if an error is detected, will
!!  advice to change the value of the conditions.
!! cond_number= number of conditions checked before calling chkint.
!! cond_string(cond_number)= name of the variables associated to the conditions.
!! cond_values(cond_number)= value of the variables associated to the conditions.
!! input_name=name of the input variable to be checked
!! input_value=value of the input variable to be checked
!! list_number=number of allowed values (maximum 40).
!! list_values=list of allowed values
!! minmax_flag=if 0, only values in the list are allowed
!!              if 1, admit values larger or equal to minmax_value
!!              if -1, admit values smaller or equal to minmax_value
!! minmax_value=see the description of minmax_flag
!! unit=unit number for clean output file
!!
!! OUTPUT
!!  (only side effect)
!!
!! SIDE EFFECT
!! ierr= switch it to 1 if an error was detected. No action otherwise.
!!
!! NOTES
!! cond_values(cond_number) or list_values(list_number)
!! must be between -99 and 999 to be printed correctly.
!!
!! for the time being, at most 3 conditions are allowed.
!!
!! in order to ask only for a minimal value, set list_number
!! as well as minmax_flag to 1, and put the minimal value in both
!! list_values and minmax_value.
!!
!! Examples :
!!  List of values - ionmov must be equal to 0, 1, 3, 8, or 9
!!   call chkint(0,0,cond_string,cond_values,ierr,'ionmov',ionmov,5,(/0,1,3,8,9/),0,0,iout)
!!
!!  Larger or equal to a given value - nberry >= limit
!!   call chkint(0,0,cond_string,cond_values,ierr,'nberry',nberry,1,(/limit/),1,limit,iout)
!!
!!  Smaller or equal to a given value - nberry <= limit
!!   call chkint(0,0,cond_string,cond_values,ierr,'nberry',nberry,1,(/limit/),-1,limit,iout)
!!
!!  Conditional cases (examples to be provided - see chkinp.f for the time being)
!!
!! PARENTS
!!      chkinp
!!
!! CHILDREN
!!      chkint_prt
!!
!! SOURCE

subroutine chkint(advice_change_cond,cond_number,cond_string,cond_values,&
                  ierr,input_name,input_value,list_number,list_values,minmax_flag,minmax_value,unit)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: advice_change_cond,cond_number,input_value,list_number
 integer,intent(in) :: minmax_flag,minmax_value,unit
 integer,intent(inout) :: ierr
 character(len=*),intent(in) :: input_name
!arrays
 integer,intent(in) :: cond_values(4),list_values(list_number)
 character(len=*),intent(inout) :: cond_string(4)

!Local variables-------------------------------
!scalars
 integer :: ilist,ok

!******************************************************************

 ! Checks the allowed values
 ok=0
 if(list_number>0)then
   do ilist=1,list_number
     if(input_value == list_values(ilist))ok=1
   end do
 end if
 if(minmax_flag==1 .and. input_value>=minmax_value)ok=1
 if(minmax_flag==-1 .and. input_value<=minmax_value)ok=1

 ! If there is something wrong, compose the message, and print it
 if(ok==0)then
   call chkint_prt(advice_change_cond,cond_number,cond_string,cond_values,&
    ierr,input_name,input_value,&
    list_number,list_values,minmax_flag,minmax_value,unit)
 end if

 ! reset all cond_strings
 cond_string(:)='#####'

end subroutine chkint
!!***

!!****f* m_parser/chkint_eq
!! NAME
!! chkint_eq
!!
!! FUNCTION
!! Checks the value of an input integer variable against a list, and
!! write a sophisticated error message when the value does not appear
!! A few conditions might have been checked before calling chkint,
!! and these are mentioned in the error message.
!!
!! See the examples in the NOTES
!!
!! INPUTS
!! advice_change_cond= if 1, and if an error is detected, will
!!  advice to change the value of the conditions.
!! cond_number= number of conditions checked before calling chkint.
!! cond_string(cond_number)= name of the variables associated to the conditions.
!! cond_values(cond_number)= value of the variables associated to the conditions.
!! input_name=name of the input variable to be checked
!! input_value=value of the input variable to be checked
!! list_number=number of allowed values (maximum 40).
!! list_values=list of allowed values
!! unit=unit number for clean output file
!!
!! OUTPUT
!!  (only side effect)
!!
!! SIDE EFFECT
!! ierr= switch it to 1 if an error was detected. No action otherwise.
!!
!! NOTES
!! cond_values(cond_number) or list_values(list_number)
!! must be between -99 and 999 to be printed correctly.
!!
!! for the time being, at most 3 conditions are allowed.
!!
!! PARENTS
!!      chkinp,m_psps
!!
!! CHILDREN
!!      chkint_prt
!!
!! SOURCE

subroutine chkint_eq(advice_change_cond,cond_number,cond_string,cond_values,&
                     ierr,input_name,input_value,list_number,list_values,unit)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: advice_change_cond,cond_number,input_value,list_number
 integer,intent(in) :: unit
 integer,intent(inout) :: ierr
 character(len=*),intent(in) :: input_name
!arrays
 integer,intent(in) :: cond_values(4),list_values(list_number)
 character(len=*),intent(inout) :: cond_string(4)

!Local variables-------------------------------
!scalars
 integer :: ilist,minmax_flag,minmax_value,ok

!******************************************************************

 !Checks the allowed values
 ok=0
 if(list_number>0)then
   do ilist=1,list_number
     if(input_value == list_values(ilist))ok=1
   end do
 end if
 minmax_flag=0
 minmax_value=0

 !If there is something wrong, compose the message, and print it
 if(ok==0)then
   call chkint_prt(advice_change_cond,cond_number,cond_string,cond_values,&
     ierr,input_name,input_value,&
     list_number,list_values,minmax_flag,minmax_value,unit)
 end if

! reset all cond_strings
 cond_string(:)='#####'

end subroutine chkint_eq
!!***

!!****f* m_parser/chkint_ge
!! NAME
!! chkint_ge
!!
!! FUNCTION
!! Checks the value of an input integer variable, expected to be greater than some value, and
!! write a sophisticated error message when it is erroneous.
!! A few conditions might have been checked before calling chkint_ge,
!! and these are mentioned in the error message.
!!
!! See the examples in the NOTES
!!
!! INPUTS
!! advice_change_cond= if 1, and if an error is detected, will
!!  advice to change the value of the conditions.
!! cond_number= number of conditions checked before calling chkint_ge.
!! cond_string(cond_number)= name of the variables associated to the conditions.
!! cond_values(cond_number)= value of the variables associated to the conditions.
!! input_name=name of the input variable to be checked
!! input_value=value of the input variable to be checked
!! minmax_value=see the description of minmax_flag
!! unit=unit number for clean output file
!!
!! OUTPUT
!!  (only side effect)
!!
!! SIDE EFFECT
!! ierr= switch it to 1 if an error was detected. No action otherwise.
!!
!! NOTES
!! cond_values(cond_number) or list_values(list_number)
!! must be between -99 and 999 to be printed correctly.
!!
!! for the time being, at most 3 conditions are allowed.
!!
!! PARENTS
!!      chkinp,invars1
!!
!! CHILDREN
!!      chkint_prt
!!
!! SOURCE

subroutine chkint_ge(advice_change_cond,cond_number,cond_string,cond_values,&
                     ierr,input_name,input_value,minmax_value,unit)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: advice_change_cond,cond_number,input_value
 integer,intent(in) :: minmax_value,unit
 integer,intent(inout) :: ierr
 character(len=*),intent(in) :: input_name
!arrays
 integer,intent(in) :: cond_values(4)
 character(len=*),intent(inout) :: cond_string(4)

!Local variables-------------------------------
!scalars
 integer :: list_number,minmax_flag,ok
 integer, allocatable :: list_values(:)

!******************************************************************

 !Checks the allowed values
 ok=0
 minmax_flag=1
 if(input_value>=minmax_value)ok=1
 list_number=1
 ABI_ALLOCATE(list_values,(1))
 list_values=minmax_value

 !If there is something wrong, compose the message, and print it
 if(ok==0)then
   call chkint_prt(advice_change_cond,cond_number,cond_string,cond_values,&
     ierr,input_name,input_value,&
     list_number,list_values,minmax_flag,minmax_value,unit)
 end if

 ABI_FREE(list_values)

 ! reset all cond_strings
 cond_string(:)='#####'

end subroutine chkint_ge
!!***

!!****f* m_parser/chkint_le
!! NAME
!! chkint_le
!!
!! FUNCTION
!! Checks the value of an input integer variable, expected to be lower than some value, and
!! write a sophisticated error message when it is erroneous.
!! A few conditions might have been checked before calling chkint_le,
!! and these are mentioned in the error message.
!!
!! See the examples in the NOTES
!!
!! INPUTS
!! advice_change_cond= if 1, and if an error is detected, will
!!  advice to change the value of the conditions.
!! cond_number= number of conditions checked before calling chkint_le.
!! cond_string(cond_number)= name of the variables associated to the conditions.
!! cond_values(cond_number)= value of the variables associated to the conditions.
!! input_name=name of the input variable to be checked
!! input_value=value of the input variable to be checked
!! minmax_value=see the description of minmax_flag
!! unit=unit number for clean output file
!!
!! OUTPUT
!!  (only side effect)
!!
!! SIDE EFFECT
!! ierr= switch it to 1 if an error was detected. No action otherwise.
!!
!! NOTES
!! cond_values(cond_number) or list_values(list_number)
!! must be between -99 and 999 to be printed correctly.
!!
!! for the time being, at most 3 conditions are allowed.
!!
!! PARENTS
!!      chkinp
!!
!! CHILDREN
!!      chkint_prt
!!
!! SOURCE

subroutine chkint_le(advice_change_cond,cond_number,cond_string,cond_values,&
                     ierr,input_name,input_value,minmax_value,unit)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: advice_change_cond,cond_number,input_value
 integer,intent(in) :: minmax_value,unit
 integer,intent(inout) :: ierr
 character(len=*),intent(in) :: input_name
!arrays
 integer,intent(in) :: cond_values(4)
 character(len=*),intent(inout) :: cond_string(4)

!Local variables-------------------------------
!scalars
 integer :: list_number,minmax_flag,ok
 integer, allocatable :: list_values(:)

!******************************************************************

 !Checks the allowed values
 ok=0
 minmax_flag=-1
 if(input_value<=minmax_value)ok=1
 !write(std_out,*)' chkint_le : input_value,minmax_value=',input_value,minmax_value

 list_number=1
 ABI_ALLOCATE(list_values,(1))
 list_values=minmax_value

 !If there is something wrong, compose the message, and print it
 if(ok==0)then
   call chkint_prt(advice_change_cond,cond_number,cond_string,cond_values,&
     ierr,input_name,input_value,list_number,list_values,minmax_flag,minmax_value,unit)
 end if

 ABI_FREE(list_values)

 ! reset all cond_strings
 cond_string(:)='#####'

end subroutine chkint_le
!!***

!!****f* m_parser/chkint_ne
!! NAME
!! chkint_ne
!!
!! FUNCTION
!! Checks the value of an input integer variable against a list, and
!! write a sophisticated error message when the value appears in the list.
!! A few conditions might have been checked before calling chkint,
!! and these are mentioned in the error message.
!!
!! See the examples in the NOTES
!!
!! INPUTS
!! advice_change_cond= if 1, and if an error is detected, will
!!  advice to change the value of the conditions.
!! cond_number= number of conditions checked before calling chkint.
!! cond_string(cond_number)= name of the variables associated to the conditions.
!! cond_values(cond_number)= value of the variables associated to the conditions.
!! input_name=name of the input variable to be checked
!! input_value=value of the input variable to be checked
!! list_number=number of NOT allowed values (maximum 40).
!! list_values=list of allowed values
!! unit=unit number for clean output file
!!
!! OUTPUT
!!  (only side effect)
!!
!! SIDE EFFECT
!! ierr= switch it to 1 if an error was detected. No action otherwise.
!!
!! NOTES
!! cond_values(cond_number) or list_values(list_number)
!! must be between -99 and 999 to be printed correctly.
!!
!! for the time being, at most 3 conditions are allowed.
!!
!! PARENTS
!!      chkinp
!!
!! CHILDREN
!!      chkint_prt
!!
!! SOURCE

subroutine chkint_ne(advice_change_cond,cond_number,cond_string,cond_values,&
                     ierr,input_name,input_value, list_number,list_values,unit)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: advice_change_cond,cond_number,input_value,list_number
 integer,intent(in) :: unit
 integer,intent(inout) :: ierr
 character(len=*),intent(in) :: input_name
!arrays
 integer,intent(in) :: cond_values(4),list_values(list_number)
 character(len=*),intent(inout) :: cond_string(4)

!Local variables-------------------------------
!scalars
 integer :: ilist,minmax_flag,minmax_value,ok

!******************************************************************

 !Checks the allowed values
 ok=1
 if(list_number>0)then
   do ilist=1,list_number
     if(input_value == list_values(ilist))ok=0
   end do
 end if
 minmax_flag=2
 minmax_value=0

 !If there is something wrong, compose the message, and print it
 if(ok==0)then
   call chkint_prt(advice_change_cond,cond_number,cond_string,cond_values,&
     ierr,input_name,input_value,&
     list_number,list_values,minmax_flag,minmax_value,unit)
 end if

 ! reset all cond_strings
 cond_string(:)='#####'

end subroutine chkint_ne
!!***

!!****f* m_parser/chkint_prt
!! NAME
!! chkint_prt
!!
!! FUNCTION
!! During the checking of the value of a variable,
!! write a sophisticated error message when it is erroneous.
!! A few conditions might have been checked before calling chkval,
!! and these are mentioned in the error message.
!!
!! See the examples in the NOTES
!!
!! INPUTS
!! advice_change_cond= if 1, and if an error is detected, will
!!  advice to change the value of the conditions.
!! cond_number= number of conditions checked before calling chkint.
!! cond_string(cond_number)= name of the variables associated to the conditions.
!! cond_values(cond_number)= value of the variables associated to the conditions.
!! input_name=name of the input variable to be checked
!! input_value=value of the input variable to be checked
!! list_number=number of allowed values (maximum 40).
!! list_values=list of allowed values
!! minmax_flag=if 0, only values in the list are allowed
!!              if 1, admit values larger or equal to minmax_value
!!              if -1, admit values smaller or equal to minmax_value
!!              if 2, values in the list are not allowed
!! minmax_value=see the description of minmax_flag
!! unit=unit number for clean output file
!!
!! OUTPUT
!!  (only side effect)
!!
!! SIDE EFFECT
!! ierr= switch it to 1 if an error was detected. No action otherwise.
!!
!! NOTES
!! cond_values(cond_number) or list_values(list_number)
!! must be between -99 and 999 to be printed correctly.
!!
!! for the time being, at most 3 conditions are allowed.
!! In order to ask only for a minimal value, set list_number
!! as well as minmax_flag to 1, and put the minimal value in both
!! list_values and minmax_value.
!!
!! Examples:
!!  List of values - ionmov must be equal to 0, 1, 3, 8, or 9
!!   call chkint_prt(0,0,cond_string,cond_values,ierr,'ionmov',ionmov,5,(/0,1,3,8,9/),0,0,iout)
!!
!!  Larger or equal to a given value - nberry >= limit
!!   call chkint_prt(0,0,cond_string,cond_values,ierr,'nberry',nberry,1,(/limit/),1,limit,iout)
!!
!!  Smaller or equal to a given value - nberry <= limit
!!   call chkint_prt(0,0,cond_string,cond_values,ierr,'nberry',nberry,1,(/limit/),-1,limit,iout)
!!
!!  Conditional cases (examples to be provided - see chkinp.f for the time being)
!!
!! PARENTS
!!      chkint,chkint_eq,chkint_ge,chkint_le,chkint_ne
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine chkint_prt(advice_change_cond,cond_number,cond_string,cond_values,&
                      ierr,input_name,input_value,list_number,list_values,minmax_flag,minmax_value,unit)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: advice_change_cond,cond_number,input_value,list_number
 integer,intent(in) :: minmax_flag,minmax_value,unit
 integer,intent(inout) :: ierr
 character(len=*),intent(in) :: input_name
!arrays
 integer,intent(in) :: cond_values(4),list_values(list_number)
 character(len=*),intent(in) :: cond_string(4)

!Local variables-------------------------------
!scalars
 integer :: icond
 character(len=500) :: msg

!******************************************************************

 if(cond_number<0 .or. cond_number>4)then
   write(msg,'(a,i0,a)' )'The value of cond_number is ',cond_number,' but it should be positive and < 5.'
   MSG_BUG(msg)
 end if

 if(list_number<0 .or. list_number>40)then
   write(msg,'(a,i0,a)' )'The value of list_number is',list_number,' but it should be between 0 and 40.'
   MSG_BUG(msg)
 end if

 !Compose the message, and print it
 ierr=1
 write(msg, '(2a)' ) ch10,' chkint_prt: ERROR -'
 if(cond_number/=0)then
   do icond=1,cond_number
     ! The following format restricts cond_values(icond) to be between -99 and 999
     write(msg, '(5a,i0,a)' ) trim(msg),ch10,&
      ' Context: the value of the variable ',trim(cond_string(icond)),' is ',cond_values(icond),'.'
   end do
 end if
 write(msg, '(5a,i0,a)' ) trim(msg),ch10,&
  '  The value of the input variable ',trim(input_name),' is ',input_value,', while it must be'
 if(minmax_flag==2)then
   write(msg, '(3a,20(i0,1x))' ) trim(msg),ch10,&
   '  different from one of the following: ',list_values(1:list_number)
 else if(list_number>1 .or. minmax_flag==0 .or. list_values(1)/=minmax_value )then
   ! The following format restricts list_values to be between -99 and 999
   if(list_number/=1)then
     write(msg, '(3a,40(i0,1x))' ) trim(msg),ch10,&
     '  equal to one of the following: ',list_values(1:list_number)
   else
     write(msg, '(3a,40(i0,1x))' ) trim(msg),ch10,'  equal to ',list_values(1)
   end if
   if(minmax_flag==1)then
     ! The following format restricts minmax_value to be between -99 and 999
     write(msg, '(3a,i0,a)' ) trim(msg),ch10,'  or it must be larger or equal to ',minmax_value,'.'
   else if(minmax_flag==-1)then
     write(msg, '(3a,i0,a)' ) trim(msg),ch10,'  or it must be smaller or equal to ',minmax_value,'.'
   end if
 else if(minmax_flag==1)then
   ! The following format restricts minmax_value to be between -99 and 999
   write(msg, '(3a,i0,a)' ) trim(msg),ch10,'  larger or equal to ',minmax_value,'.'
 else if(minmax_flag==-1)then
   ! The following format restricts minmax_value to be between -99 and 999
   write(msg, '(3a,i0,a)' ) trim(msg),ch10,'  smaller or equal to ',minmax_value,'.'
 end if
 if(cond_number==0 .or. advice_change_cond==0)then
   write(msg, '(5a)' ) trim(msg),ch10,'  Action: you should change the input variable ',trim(input_name),'.'
 else if(cond_number==1)then
   write(msg, '(7a)' ) trim(msg),ch10,&
    '  Action: you should change the input variables ',trim(input_name),' or ',trim(cond_string(1)),'.'
 else if(cond_number==2)then
   write(msg, '(11a)' ) trim(msg),ch10,&
    '  Action: you should change one of the input variables ',trim(input_name),',',ch10,&
    '   ',trim(cond_string(1)),' or ',trim(cond_string(2)),'.'
 else if(cond_number==3)then
   write(msg, '(13a)' ) trim(msg),ch10,&
    '  Action: you should change one of the input variables ',trim(input_name),',',ch10,&
    '   ',trim(cond_string(1)),', ',trim(cond_string(2)),' or ',trim(cond_string(3)),'.'
 end if
 call wrtout([unit, std_out], msg)

end subroutine chkint_prt
!!***

!!****f* m_parser/prttagm
!!
!! NAME
!! prttagm
!!
!! FUNCTION
!! Eventually print the content of dprarr (if typevarphys='DPR','LEN', 'ENE', 'TIM' and 'BFI'),
!! or intarr (if typevarphys='INT'), arrays of effective dimensions narr and 0:ndtset_alloc
!! For the second dimension, the 0 index relates to a default.
!! Print the array only if the content for at least one value of the second
!! index is different from the default.
!! Print a generic value if the non-default values are all equal.
!! Print the detail of all values otherwise.
!! The input variable 'length' controls the print format, and, in the case
!! of the real(dp) variable, the way two numbers are determined to be
!! different or not.
!!
!! INPUTS
!!  intarr(1:marr,0:ndtset_alloc), dprarr(1:marr,0:ndtset_alloc)
!!   integer or real(dp) arrays, respectively,
!!   containing the data to be printed. Use these arrays even for scalars.
!!   For the first index, only the range 1:narr is relevant.
!!  iout=unit number for echoed output
!!  jdtset_(0:ndtset_alloc)=list of dataset indices.
!!  length= if 1, short format for printing, if 2, long format for printing
!!     special formats: if 3, INT : for symrel or kptrlatt
!!                      if 4, INT : for type
!!                      if 5, INT : for mkmem, mkqmem, mk1mem
!!                      if 6, INT : for kptrlatt
!!                      if 3, DPR : for tnons
!!                      if 4, DPR : for wtk and znucl
!!                      if 5, DPR : for atvshift
!!                      if 6, DPR : very short format for printing
!!     If the typevarphys is 'DPR', a negative value of 'length' will request that
!!        the equality of real(dp) numbers is determined by an ABSOLUTE
!!        difference criterion only. The absolute value of length is used
!!        to determine the format, as above.
!!
!!  marr=first dimension of the intarr and dprarr arrays, as declared in the
!!   calling subroutine.
!!  narr=actual first dimension of intarr and dprarr.
!!  narrm=used when the effective first dimension of intarr is variable
!!        in this case narrm(0:ndtset_alloc)
!!  ncid= NETCDF id
!!  ndtset_alloc=govern second dimension of intarr and dprarr
!!  token=character string for 'tag'.  Assumed no longer than 9 characters
!!  typevarphys=physical variable type (might indicate the physical meaning of
!!   for dimensionality purposes)
!!   'INT'=>integer
!!   'DPR'=>real(dp) (no special treatment)
!!   'LEN'=>real(dp) (output in bohr and angstrom)
!!   'ENE'=>real(dp) (output in hartree and eV)
!!   'BFI'=>real(dp) (output in Tesla)
!!   'TIM'=>real(dp) (output in second)
!!  use_narrm= if 0, use of scalar 'narr' instead of array 'narrm'
!!  [firstchar]= (optional) first character of the line (default=' ')
!!  [forceprint]= (optional) control if output is forced even if a variable is equal to its default value:
!!                0: not printed out if equal to default value
!!                1: output forced even if equal to default value in both TEXT and NETCDF file
!!                2: output forced even if equal to default value in NETCDF file only
!!                3: output forced even if equal to default value in TEXT file only
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      outvar_a_h,outvar_i_n,outvar_o_z,pawuj_det,prttagm_images
!!
!! CHILDREN
!!      appdig,write_var_netcdf
!!
!! SOURCE

subroutine prttagm(dprarr,intarr,iout,jdtset_,length,&
                    marr,narr,narrm,ncid,ndtset_alloc,token,typevarphys,use_narrm,&
                    firstchar,forceprint)  ! optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,length,marr,narr,ndtset_alloc,ncid,use_narrm
 integer,intent(in),optional :: forceprint
 character(len=*),intent(in) :: token
 character(len=3),intent(in) :: typevarphys
 character(len=1),intent(in),optional :: firstchar
!arrays
 integer,intent(in) :: intarr(marr,0:ndtset_alloc)
 integer,intent(in) :: jdtset_(0:ndtset_alloc)
 integer,intent(in) :: narrm(0:ndtset_alloc)
 real(dp),intent(in) :: dprarr(marr,0:ndtset_alloc)

!Local variables-------------------------------
!character(len=*), parameter :: long_beg     ='(a,a16,a,1x,(t22,'
 character(len=*), parameter :: format_1     ='",a16,a,t22,'
 character(len=*), parameter :: format_2     ='",t22,'
 character(len=*), parameter :: short_int    ='10i5)'
 character(len=*), parameter :: long_int     ='8i8)'
 character(len=*), parameter :: veryshort_dpr='f11.5)'
 character(len=*), parameter :: short_dpr    ='es16.8)'
 character(len=*), parameter :: long_dpr     ='es18.10)'
 character(len=*), parameter :: veryshort_dim='f11.5),a'
 character(len=*), parameter :: short_dim    ='es16.8),a'
 character(len=*), parameter :: long_dim     ='es18.10),a'
 character(len=*), parameter :: f_symrel     ='3(3i3,1x),4x,3(3i3,1x))'
 character(len=*), parameter :: f_type       ='20i3)'
 character(len=*), parameter :: f_mem        ='8i8)'
 character(len=*), parameter :: f_tnons      ='3f11.7,3x,3f11.7)'
 character(len=*), parameter :: f_wtk        ='6f11.5)'
 character(len=*), parameter :: f_atvshift   ='5f11.5)'
 character(len=*), parameter :: f_kptrlatt   ='3(3i5,2x))'
!scalars
 integer :: iarr,idtset,jdtset,multi,ndtset_eff,narr_eff
 logical :: print_netcdf,print_out
 real(dp),parameter :: tol21=1.0d-21
 real(dp) :: diff,scale_factor,sumtol
 character(len=4) :: digit
 character(len=1) :: first_column
 character(len=4) :: appen
 character(len=8) :: out_unit
 character(len=50) :: format_dp,format_int,full_format
 character(len=500) :: msg

! *************************************************************************

!###########################################################
!### 01. Check consistency of input

 if(len_trim(token)>16)then
   write(msg, '(3a,i0,2a)' )&
   'The length of the name of the input variable ',trim(token),' is ',len_trim(token),ch10,&
   'This exceeds 16 characters, the present maximum in routine prttagm.'
   MSG_ERROR(msg)
 end if

 if(ndtset_alloc<1)then
   write(msg, '(a,i0,a,a,a,a,a)' )&
   'ndtset_alloc=',ndtset_alloc,', while it should be >= 1.',ch10,&
   'This happened for token=',token,'.'
   MSG_BUG(msg)
 end if

 if(ndtset_alloc>9999)then
   write(msg, '(a,i0,a,a,a,a,a)' )&
   'ndtset_alloc=',ndtset_alloc,', while it must be lower than 10000.',ch10,&
   'This happened for token=',token,'.'
   MSG_BUG(msg)
 end if

 if(narr>99 .and. (typevarphys=='ENE'.or.typevarphys=='LEN'))then
   write(msg, '(3a,i0,a)' )' typevarphys=',typevarphys,' with narr=',narr,'  is not allowed.'
   MSG_BUG(msg)
 end if

 if ((narr>0).or.(use_narrm/=0)) then

   print_out=.true.;print_netcdf=.true.
   multi=0

!  ###########################################################
!  ### 02. Treatment of integer 'INT'

   if(typevarphys=='INT')then

!    Determine whether the different non-default occurences are all equal

     if (use_narrm==0) then ! use of scalar 'narr' instead of array 'narrm'
       if(ndtset_alloc>1)then
         do idtset=1,ndtset_alloc
           do iarr=1,narr
             if(intarr(iarr,1)/=intarr(iarr,idtset))multi=1
           end do
         end do
       end if
     else
!      If the sizes of the arrays are different we can not compare them
!      So we have to assume they are different
       multi=1
     end if

!    If they are all equal, then determine whether they are equal to the default
     if(multi==0)then
       print_out=.false.
       do iarr=1,narr
         if (intarr(iarr,1)/=intarr(iarr,0)) print_out=.true.
       end do
       print_netcdf=print_out
     end if

     if (present(forceprint)) then
       if (forceprint==1.or.forceprint==3) print_out=.true.
       if (forceprint==1.or.forceprint==2) print_netcdf=.true.
     end if

!    Print only if the values differ from the default
     if (print_out.or.print_netcdf.or.(ncid<0))then
       ndtset_eff=ndtset_alloc
       if((multi==0).or.(ncid<0)) ndtset_eff=1
       do idtset=1,ndtset_eff

!        Initialize the character in the first column
         first_column=' ';if (present(firstchar)) first_column=firstchar
         if(abs(length)==5)first_column='P'
!        Initialize the format
         if(abs(length)==1)format_int=trim(short_int)
         if(abs(length)==2)format_int=trim(long_int)
         if(abs(length)==3)format_int=trim(f_symrel)
         if(abs(length)==4)format_int=trim(f_type)
         if(abs(length)==5)format_int=trim(f_mem)
         if(abs(length)==6)format_int=trim(f_kptrlatt)
!        Initialize the dataset number string, and print
         if((multi==0).or.(ncid<0))then
           appen=' '
         else
           jdtset=jdtset_(idtset)
           call appdig(jdtset,'',appen)
         end if
!        full_format=trim(long_beg)//trim(format_int)
         full_format='("'//first_column//trim(format_1)//'("'// first_column//trim(format_2)//trim(format_int)//")"

!        narr_eff could be narr or narrm(idtset)
!        It depends if the size is variable for different datasets
         if (use_narrm==0)then
           narr_eff=narr
         else
           narr_eff=narrm(idtset)
         end if

         if (narr_eff/=0) then

           if (print_out) write(iout,full_format) token,trim(appen),intarr(1:narr_eff,idtset)
#ifdef HAVE_NETCDF
           if (print_netcdf) then
             call write_var_netcdf(intarr(1:narr_eff,idtset),&
&             dprarr(1:narr_eff,idtset),marr,narr_eff,abs(ncid),typevarphys,token//appen)
           end if
#endif
         end if

       end do
     end if !(print==1)

!    ###########################################################
!    ### 03. Treatment of real 'DPR', 'LEN', 'ENE', 'BFI', 'TIM'

   else if (typevarphys=='DPR' .or. typevarphys=='LEN' .or. typevarphys=='ENE' .or. typevarphys=='BFI' .or. typevarphys=='TIM') then

     if((ndtset_alloc>1).and.(use_narrm==0))then
       do idtset=1,ndtset_alloc
         do iarr=1,narr
!          The determination of effective equality is more difficult than in the
!          integer case :
!          - if length > 0, ask for a relative accuracy, and also include
!          the case of zero values, thanks to tol21.
!          - if length < 0, ask for absolute accuracy.
           diff=abs( dprarr(iarr,1)-dprarr(iarr,idtset) )
           if(length>0)then
             sumtol=abs(dprarr(iarr,1))+abs(dprarr(iarr,idtset))+10*tol21
             if(diff>sumtol*tol11)multi=1
           else
             if(diff>tol14)multi=1
           end if
         end do
       end do
     elseif (use_narrm/=0) then
       multi=1 ! Assume that values could not be compared between different datasets.
!      Nevertheless, checks whether not all dataset might be equal to the default, despite varying dimensions (e.g. all zeroes)
       print_out=.false.
       do idtset=1,ndtset_alloc
         if(narrm(idtset)>narrm(0))then
           print_out=.true.
         else
           do iarr=1,narrm(idtset)
             diff=abs( dprarr(iarr,idtset)-dprarr(iarr,0) )
             if(length>0)then
               sumtol=abs(dprarr(iarr,idtset))+abs(dprarr(iarr,0))+10*tol21
               if(diff>sumtol*tol11)print_out=.true.
             else
               if(diff>tol14)print_out=.true.
             end if
           end do
         end if
       end do
       print_netcdf=print_out
     end if

     if(multi==0)then
       print_out=.false.
       do iarr=1,narr
         diff=abs( dprarr(iarr,1)-dprarr(iarr,0) )
         if(length>0)then
           sumtol=abs(dprarr(iarr,1))+abs(dprarr(iarr,0))+10*tol21
           if(diff>sumtol*tol11)print_out=.true.
         else
           if(diff>tol14)print_out=.true.
         end if
       end do
       print_netcdf=print_out
     end if

     if (present(forceprint)) then
       if (forceprint==1.or.forceprint==3) print_out=.true.
       if (forceprint==1.or.forceprint==2) print_netcdf=.true.
     end if

     if(print_out.or.print_netcdf.or.(ncid<0))then
!      Select the proper format
       ndtset_eff=ndtset_alloc
       if((multi==0).or.(ncid<0))ndtset_eff=1
       narr_eff=narr
       if(use_narrm/=0)then
         narr_eff=maxval(narrm(1:ndtset_eff))
       end if
       if(abs(length)==1 .or. abs(length)==2 .or. abs(length)==6)then
         if(typevarphys=='DPR')then
           digit='3'
           if(abs(length)==1)format_dp=digit//short_dpr
           if(abs(length)==2)format_dp=digit//long_dpr
           if(abs(length)==6)format_dp=digit//veryshort_dpr
   else if(typevarphys=='ENE' .or. typevarphys=='LEN' .or. typevarphys=='BFI' .or. typevarphys=='TIM')then
           if (narr<10) write(digit,'(i1)')narr_eff
           if (narr> 9) write(digit,'(i2)')narr_eff
           if(abs(length)==1)format_dp=digit//short_dim
           if(abs(length)==2)format_dp=digit//long_dim
           if(abs(length)==6)format_dp=digit//veryshort_dim
         end if
       else
         if(abs(length)==3)format_dp=f_tnons
         if(abs(length)==4)format_dp=f_wtk
         if(abs(length)==5)format_dp=f_atvshift
       end if
       do idtset=1,ndtset_eff

!        narr_eff could be narr or narrm(idtset)
!        It depends if the size is variable for different datasets
         if (use_narrm==0)then
           narr_eff=narr
         else
           narr_eff=narrm(idtset)
         end if

         if (narr_eff/=0) then

!          Initialize the character in the first column
           first_column=' ';if (present(firstchar)) first_column=firstchar
!          Define scale_factor
           scale_factor=one !EB to what this is still usefull ???
!          EB remove           if(typevarphys=='BFI')scale_factor=one/BField_Tesla
!          Define out_unit
           if(typevarphys=='ENE')out_unit=' Hartree'
           if(typevarphys=='LEN')out_unit=' Bohr   '
           if(typevarphys=='BFI')out_unit='   ' !EB remove Tesla unit
           if(typevarphys=='TIM')out_unit=' Second'
!          Format, according to the length of the dataset string
           if((multi==0).or.(ncid<0))then
             appen=' '
           else
             jdtset=jdtset_(idtset)
             call appdig(jdtset,'',appen)
           end if
           ! full_format=trim(long_beg)//trim(format_dp)
           full_format='("'//first_column//trim(format_1)//'("'// first_column//trim(format_2)//trim(format_dp)//")"
           ! write(ab_out,*)' trim(long_beg)=',trim(long_beg)
           ! write(ab_out,*)' trim(format_dp)=',trim(format_dp)
           ! write(ab_out,*)' trim(full_format)=',trim(full_format)
           if(typevarphys=='DPR')then
             if (print_out) write(iout,full_format) token,trim(appen),dprarr(1:narr_eff,idtset)*scale_factor
           else
             if (print_out) write(iout,full_format) token,trim(appen),dprarr(1:narr_eff,idtset)*scale_factor,trim(out_unit)
           end if
#ifdef HAVE_NETCDF
           if (print_netcdf) then
             call write_var_netcdf(intarr(1:narr_eff,idtset),dprarr(1:narr_eff,idtset),&
               marr,narr_eff,abs(ncid),'DPR',token//trim(appen))
           end if
#endif

         end if

       end do
     end if

!    ###########################################################
!    ### 04. The type is neither 'INT' nor 'DPR','ENE','LEN','BFI','TIM'
   else
     MSG_BUG('Disallowed typevarphys = '//TRIM(typevarphys))
   end if

 end if ! End condition of narr>0

end subroutine prttagm
!!***

!!****f* m_parser/prttagm_images
!!
!! NAME
!! prttagm_images
!!
!! FUNCTION
!! Extension to prttagm to include the printing of
!! images information, in those cases the same variable
!! is printed several times for each dataset
!!
!! Cases where images information are relevant includes xcart, xred, acell, fcart.
!!
!! INPUT
!! (see prttagm.F90)
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      outvar_a_h,outvar_i_n,outvar_o_z
!!
!! CHILDREN
!!      appdig,prttagm,write_var_netcdf
!!
!! SOURCE

subroutine prttagm_images(dprarr_images,iout,jdtset_,length,&
& marr,narrm,ncid,ndtset_alloc,token,typevarphys,&
& mxnimage,nimagem,ndtset,prtimg,strimg,firstchar,forceprint)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,length,marr,ndtset_alloc,ncid
 integer,intent(in) :: mxnimage,ndtset
 integer,intent(in),optional :: forceprint
 character(len=*),intent(in) :: token
 character(len=3),intent(in) :: typevarphys
 character(len=1),intent(in),optional :: firstchar
!arrays
 integer,intent(in) :: prtimg(mxnimage,0:ndtset_alloc)
 integer,intent(in) :: jdtset_(0:ndtset_alloc)
 integer,intent(in) :: nimagem(0:ndtset_alloc)
 character(len=8),intent(in) :: strimg(mxnimage)
 integer,intent(in) :: narrm(0:ndtset_alloc)
 real(dp),intent(in) :: dprarr_images(marr,mxnimage,0:ndtset_alloc)

!Local variables-------------------------------
 integer :: iarr,idtset,iimage,jdtset,multi_narr,narr
 integer :: intarr_images(marr,mxnimage,0:ndtset_alloc)
 integer,allocatable :: intarr(:,:)
 real(dp), allocatable :: dprarr(:,:)
 logical :: print_out,print_netcdf,test_multiimages
 character(len=1) :: first_column
 character(len=4) :: appen
 character(len=16) :: keywd
 character(len=50) :: full_format
 character(len=*), parameter :: format_1  ='",a16,t22,'
 character(len=*), parameter :: format_1a ='",a16,a,t22,'
 character(len=*), parameter :: format_2  ='",t22,'
 character(len=*), parameter :: long_dpr  ='3es18.10)'

! *************************************************************************

!Test whether for this variable, the content of different images differ.
!test_multiimages=.false. if, for all datasets, the content is identical.
 test_multiimages=.false.
 do idtset=1,ndtset_alloc
   if(nimagem(idtset)>1)then
     do iarr=1,narrm(idtset)
       if(sum(abs( dprarr_images(iarr,2:nimagem(idtset),idtset)- &
&       dprarr_images(iarr,1              ,idtset)))>tol12)then
         test_multiimages=.true.
       end if
     end do
   end if
 end do

 if(nimagem(0)==0)test_multiimages=.true.

!If there is no differences between images, one is back to the usual prttagm routine.
!Note the treatment of firstchar and forceprint has to be transmitted to prttagm.
 if(.not.test_multiimages)then

   narr=narrm(1)
   ABI_ALLOCATE(intarr,(marr,0:ndtset_alloc))
   ABI_ALLOCATE(dprarr,(marr,0:ndtset_alloc))
   do idtset=0,ndtset_alloc
     dprarr(1:narrm(idtset),idtset)=dprarr_images(1:narrm(idtset),1,idtset)
   end do
   multi_narr=0
   if(ndtset_alloc>1)then
     do idtset=1,ndtset_alloc
       if(narrm(1)/=narrm(idtset))multi_narr=1
     end do
   end if
   if (present(firstchar).and.present(forceprint)) then
     call prttagm(dprarr,intarr,iout,jdtset_,length,marr,narr,&
       narrm,ncid,ndtset_alloc,token,typevarphys,multi_narr,&
       firstchar=firstchar,forceprint=forceprint)
   else if (present(firstchar)) then
     call prttagm(dprarr,intarr,iout,jdtset_,length,marr,narr,&
       narrm,ncid,ndtset_alloc,token,typevarphys,multi_narr,&
       firstchar=firstchar)
   else if (present(forceprint)) then
     call prttagm(dprarr,intarr,iout,jdtset_,length,marr,narr,&
       narrm,ncid,ndtset_alloc,token,typevarphys,multi_narr,&
       forceprint=forceprint)
   else
     call prttagm(dprarr,intarr,iout,jdtset_,length,marr,narr,&
       narrm,ncid,ndtset_alloc,token,typevarphys,multi_narr)
   end if
   ABI_FREE(intarr)
   ABI_FREE(dprarr)

 else

   first_column=' ';if (present(firstchar)) first_column=firstchar

   do idtset=1,ndtset_alloc

     if (narrm(idtset)>0)then
       do iimage=1,nimagem(idtset)

         print_out=.true.
         if (prtimg(iimage,idtset)==0) print_out=.false.
         if (nimagem(0)>=nimagem(idtset)) then
           if (sum(abs(dprarr_images(1:narrm(idtset),iimage,idtset) &
&           -dprarr_images(1:narrm(idtset),iimage,0)))<tol12) print_out=.false.
         end if
         print_netcdf=print_out

         if (present(forceprint)) then
           if (forceprint==1.or.forceprint==3) print_out=.true.
           if (forceprint==1.or.forceprint==2) print_netcdf=.true.
         end if

         if (print_out.or.print_netcdf.or.(ncid<0))then
           keywd=token//trim(strimg(iimage))

           if(ndtset>0)then
             jdtset=jdtset_(idtset)
             call appdig(jdtset,'',appen)
             if (print_out) then
               full_format='("'//first_column//trim(format_1a)//'("'// &
&               first_column//trim(format_2)//trim(long_dpr)//")"
               write(iout,full_format) &
&               trim(keywd),appen,dprarr_images(1:narrm(idtset),iimage,idtset)
             end if
#ifdef HAVE_NETCDF
             if (print_netcdf) then
               call write_var_netcdf(intarr_images(1:narrm(idtset),iimage,idtset),&
&               dprarr_images(1:narrm(idtset),iimage,idtset),&
&               marr,narrm(idtset),ncid,'DPR',trim(keywd)//appen)
             end if
#endif
           else

             if (print_out) then
               full_format='("'//first_column//trim(format_1)//'("'// &
&               first_column//trim(format_2)//trim(long_dpr)//")"
               write(iout,full_format) &
&               trim(keywd),dprarr_images(1:narrm(idtset),iimage,idtset)
             end if
#ifdef HAVE_NETCDF
             if (print_netcdf) then
               call write_var_netcdf(intarr_images(1:narrm(idtset),iimage,idtset),&
&               dprarr_images(1:narrm(idtset),iimage,idtset),&
&               marr,narrm(idtset),abs(ncid),'DPR',trim(keywd))
             end if
#endif

           end if
         end if
       end do
     end if
   end do

 end if

end subroutine prttagm_images
!!***

!!****f* m_parser/chkvars_in_string
!! NAME
!!  chkvars_in_string
!!
!! FUNCTION
!!  Analyze variable names in string. Ignore tokens withing double quotation marks.
!!  Abort if name is not recognized.
!!
!! INPUTS
!!  protocol=
!!    0 if parser does not accept multiple datasets and +* syntax (e.g. anaddb)
!!    1 if parser accepts multiple datasets and +* syntax (e.g. abinit)
!!
!!  list_vars(len=*)=string with the (upper case) names of the variables (excluding logicals and chars).
!!  list_logicals(len=*)=string with the (upper case) names of the logical variables.
!!  list_strings(len=*)=string with the (upper case) names of the character variables.
!!  string(len=*)=string (with upper case) from the input file.
!!
!! OUTPUT
!!  Abort if variable name is not recognized.
!!
!! PARENTS
!!      chkvars,m_anaddb_dataset
!!
!! CHILDREN
!!
!! SOURCE

subroutine chkvars_in_string(protocol, list_vars, list_logicals, list_strings, string)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: protocol
 character(len=*),intent(in) :: string
 character(len=*),intent(in) :: list_logicals,list_strings,list_vars

!Local variables-------------------------------
 character,parameter :: blank=' '
!scalars
 integer :: index_blank,index_current,index_endword,index_endwordnow,index_list_vars
 character(len=500) :: msg

!************************************************************************

 !write(std_out,"(3a)")"Checking vars in string:", ch10, trim(string)

 index_current=1
 do
   ! Infinite do-loop, to identify the presence of each potential variable names

   if(len_trim(string)<=index_current)exit
   index_blank=index(string(index_current:),blank)+index_current-1

   if(index('ABCDEFGHIJKLMNOPQRSTUVWXYZ',string(index_current:index_current))/=0)then

     index_endword = index_blank -1
     if (protocol == 1) then
       ! Skip characters like : + or the digits at the end of the word
       ! Start from the blank that follows the end of the word
       do index_endword=index_blank-1,index_current,-1
         if(index('ABCDEFGHIJKLMNOPQRSTUVWXYZ',string(index_endword:index_endword))/=0)exit
       end do
     end if
     !write(std_out,*)"Will analyze:", string(index_current:index_endword)

     ! Find the index of the potential variable name in the list of variables
     index_list_vars=index(list_vars,blank//string(index_current:index_endword)//blank)

     ! Treat the complications due to the possibility of images
     if (index_list_vars==0 .and. protocol==1) then

       ! Treat possible LASTIMG appendix
       if(index_endword-6>=1)then
         if(string(index_endword-6:index_endword)=='LASTIMG')index_endword=index_endword-7
       end if

       ! Treat possible IMG appendix
       if(index_endword-2>=1)then
         if(string(index_endword-2:index_endword)=='IMG')index_endword=index_endword-3
       end if

       index_endwordnow=index_endword

       ! Again skip characters like : + or the digits before IMG
       ! Start from the blank that follows the end of the word
       do index_endword=index_endwordnow,index_current,-1
         if(index('ABCDEFGHIJKLMNOPQRSTUVWXYZ',string(index_endword:index_endword))/=0)exit
       end do

       ! Find the index of the potential variable name in the list of variables
       index_list_vars=index(list_vars,blank//string(index_current:index_endword)//blank)
     end if

     if(index_list_vars==0)then

       ! Treat possible logical input variables
       if(index(list_logicals,blank//string(index_current:index_endword)//blank)/=0)then
         !write(std_out,*)"Found logical variable: ",string(index_current:index_endword)
         index_blank=index(string(index_current:),blank)+index_current-1
         if(index(' F T ',string(index_blank:index_blank+2))==0)then
           write(msg, '(8a)' )&
            'Found token `',string(index_current:index_endword),'` in the input file.',ch10,&
            'This variable should be given a logical value (T or F), but the following string was found:',&
            string(index_blank:index_blank+2),ch10,&
            'Action: check your input file. You likely misused the input variable.'
            MSG_ERROR(msg)
         else
           index_blank=index_blank+2
         end if

       else if(index(list_strings,blank//string(index_current:index_endword)//blank)/=0)then
         ! Treat possible string input variables
         ! Every following string is accepted
         !write(std_out,*)"Found string variable: ",string(index_current:index_endword)
         !write(std_out,*)"in string: ",trim(string(index_current:))
         index_current=index(string(index_current:),blank)+index_current
         index_blank=index(string(index_current:),blank)+index_current-1
         !write(std_out,*)"next:: ",string(index_current:index_endword)

       else
         ! If still not admitted, then there is a problem
         write(msg, '(7a)' )&
         'Found token: `',string(index_current:index_endword),'` in the input file.',ch10,&
         'This name is not one of the registered input variable names (see https://docs.abinit.org/).',ch10,&
         'Action: check your input file. You likely mistyped the input variable.'
         MSG_ERROR(msg)
       end if
     end if
   end if

   index_current=index_blank+1

   if (string(index_current:index_current) == '"') then
     do
       index_current = index_current + 1
       if (string(index_current:index_current) == '"') exit
       if (index_current > len_trim(string)) then
         MSG_ERROR('Cannot find closing quotation mark " in string. You likely forgot to close a string')
       end if
     end do

   end if

 end do

end subroutine chkvars_in_string
!!***

!!****f* m_parser/geo_from_abivar_string
!! NAME
!!  geo_from_abivars_string
!!
!! FUNCTION
!!  Build object form abinit `structure` variable
!!
!! INPUTS
!!  comm=MPI communicator. Used for performing IO.
!!
!! SOURCE

type(geo_t) function geo_from_abivar_string(string, comm) result(new)
!type(geo_t) function geo_from_structure_string(string, comm) result(new)

!Arguments ------------------------------------
 character(len=*),intent(in) :: string
 integer,intent(in) :: comm

!Local variables-------------------------------
 integer :: ii
 character(len=len(string)) :: prefix

!************************************************************************

 !print *, "in geo_from_abivar_string: `", trim(string), "`"

 ii = index(string, ":")
 ABI_CHECK(ii > 0, sjoin("Expecting string of the form `type:content`, got:", string))
 prefix = adjustl(string(1:ii-1))

 select case (prefix)

 case ("poscar")
   ! Build geo ifrom POSCAR from file.
   new = geo_from_poscar_path(trim(string(ii+1:)), comm)

 case ("abivars")
   ! Build geo from from file with Abinit variables.
   new = geo_from_abivars_path(trim(string(ii+1:)), comm)

 case ("abifile")
   if (endswith(string(ii+1:), ".nc")) then
     ! Build geo from netcdf file.
     new = geo_from_netcdf_path(trim(string(ii+1:)), comm)
   else
     ! Assume Fortran file with Abinit header.
     MSG_ERROR("structure variable with Fortran file is not yet implemented.")
     !new = geo_from_fortran_file_with_hdr(string(ii+1:), comm)
     !cryst = crystal_from_file(string(ii+1:), comm)
     !if (cryst%isalchemical()) then
     !  MSG_ERROR("Alchemical mixing is not compatibile with `structure` input variable!")
     !end if
     !new%natom = cryst%natom
     !new%ntypat = cryst%ntypat
     !new%rprimd = cryst%rprimd
     !call alloc_copy(cryst%typat, new%typat)
     !call alloc_copy(cryst%xred, new%xred)
     !call alloc_copy(cryst%znucl, new%znucl)
     !call cryst%free()
   end if

 case default
   MSG_ERROR(sjoin("Invalid prefix: `", prefix, "`"))
 end select

end function geo_from_abivar_string
!!***

!!****f* m_parser/geo_from_abivars_path
!! NAME
!!  geo_from_abivars_path
!!
!! FUNCTION
!!
!! SOURCE

type(geo_t) function geo_from_abivars_path(path, comm) result(new)

!Arguments ------------------------------------
 character(len=*),intent(in) :: path
 integer,intent(in) :: comm

!Local variables-------------------------------
 integer,parameter :: master = 0, option1 = 1
 integer :: jdtset, iimage, nimage, iatom, itypat
 integer :: my_rank, lenstr, ierr, ii, start, tread, marr
 !character(len=500) :: msg
 character(len=strlen) :: string
!arrays
 integer,allocatable :: intarr(:)
 real(dp) :: acell(3), rprim(3,3)
 real(dp),allocatable :: dprarr(:)
 character(len=5),allocatable :: symbols(:)

!************************************************************************

 ! Master node reads string and broadcasts
 my_rank = xmpi_comm_rank(comm)

 if (my_rank == master) then
   ! Below part copied from `parsefile`. strlen from defs_basis module
   call instrng(path, lenstr, option1, strlen, string)
   ! To make case-insensitive, map characters of string to upper case.
   call inupper(string(1:lenstr))
   !call chkvars_in_string(protocol1, list_vars, list_logicals, list_strings, string)
 end if

 if (xmpi_comm_size(comm) > 1) then
   call xmpi_bcast(string, master, comm, ierr)
   call xmpi_bcast(lenstr, master, comm, ierr)
 end if

 ! ==============================
 ! Now all procs parse the string
 ! ==============================

 jdtset = 0; iimage = 0; nimage = 0

 ! Get the number of atom in the unit cell. Read natom from string
 marr = 1
 ABI_MALLOC(intarr, (marr))
 ABI_MALLOC(dprarr, (marr))

 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'natom', tread, 'INT')
 ABI_CHECK(tread /= 0, sjoin("natom is required in file:", path))
 new%natom = intarr(1)

 marr = max(12, 3*new%natom)
 ABI_REMALLOC(intarr, (marr))
 ABI_REMALLOC(dprarr, (marr))

 ! Set up unit cell from acell, rprim, angdeg
 call get_acell_rprim(lenstr, string, jdtset, iimage, nimage, marr, acell, rprim)

 ! Compute different matrices in real and reciprocal space, also checks whether ucvol is positive.
 call mkrdim(acell, rprim, new%rprimd)

 ! Parse atomic positions.
 ! Only xcart is supported here because it makes life easier and we don't need to handle symbols + Units
 ii = index(string(1:lenstr), "XRED_SYMBOLS")
 ABI_CHECK(ii /= 0, "In structure mode only `xred_symbols` with coords followed by element symbol are supported")

 new%fileformat = "abivars"
 ABI_MALLOC(new%xred, (3, new%natom))

 ABI_MALLOC(symbols, (new%natom))
 start = ii + len("XRED_SYMBOLS")
 do iatom=1,new%natom
   call inarray(start, "xred_symbols", dprarr, intarr, marr, 3, string, "DPR")
   new%xred(:, iatom) = dprarr(1:3)
   ABI_CHECK(next_token(string, start, symbols(iatom)) == 0, "Error while reading element symbol.")
   symbols(iatom) = tolower(symbols(iatom))
   symbols(iatom)(1:1) = toupper(symbols(iatom)(1:1))
   !write(std_out, *)"xred", new%xred(:, iatom), "symbol:", trim(symbols(iatom))
 end do

 call typat_from_symbols(symbols, new%ntypat, new%typat)

 ! Note that the first letter should be capitalized, rest must be lower case
 ABI_MALLOC(new%znucl, (new%ntypat))
 do iatom=1,new%natom
   itypat = new%typat(iatom)
   new%znucl(itypat) = symbol2znucl(symbols(iatom))
 end do

 ABI_FREE(symbols)
 ABI_FREE(intarr)
 ABI_FREE(dprarr)

 !call new%print_abivars(std_out)

contains

subroutine typat_from_symbols(symbols, ntypat, typat)

!Arguments ------------------------------------
 character(len=*),intent(in) :: symbols(:)
 integer,intent(out) :: ntypat
 integer,allocatable,intent(out) :: typat(:)

!Local variables-------------------------------
 integer :: ii, jj, nstr, found

!************************************************************************

 nstr = size(symbols)
 ABI_ICALLOC(typat, (nstr))

 typat(1) = 1
 ntypat = 1
 do ii=2, nstr
   found = 0
   do jj=1, ntypat
     if (symbols(ii) == symbols(typat(jj))) then
       found = jj; exit
     end if
   end do
   if (found == 0) then
     ntypat = ntypat + 1
     typat(ii) = ntypat
   else
     typat(ii) = found
   end if
 end do

end subroutine typat_from_symbols

end function geo_from_abivars_path
!!***

!!****f* m_parser/geo_from_poscar_path
!! NAME
!!  geo_from_poscar_path
!!
!! FUNCTION
!!
!! SOURCE

type(geo_t) function geo_from_poscar_path(path, comm) result(new)

!Arguments ------------------------------------
 character(len=*),intent(in) :: path
 integer,intent(in) :: comm

!Local variables-------------------------------
 integer,parameter :: master = 0
 integer :: unt, my_rank
 character(len=500) :: msg

!************************************************************************

 my_rank = xmpi_comm_rank(comm)

 if (my_rank == master) then
   if (open_file(path, msg, newunit=unt, form='formatted', status='old', action="read") /= 0) then
     MSG_ERROR(msg)
   end if
   new = geo_from_poscar_unit(unt)
   close(unt)
 end if

 if (xmpi_comm_size(comm) > 1) call new%bcast(master, comm)

end function geo_from_poscar_path
!!***

!!****f* m_parser/geo_from_poscar_unit
!! NAME
!!  geo_from_poscar_unit
!!
!! FUNCTION
!!  Build object from string with seperator `sep`. Usually sep = newline = ch10
!!
!! SOURCE

type(geo_t) function geo_from_poscar_unit(unit) result(new)

!Arguments ------------------------------------
 integer,intent(in) :: unit

!Local variables-------------------------------
 !integer,parameter :: marr = 3
 integer :: beg, iatom, itypat, ierr, ii, cnt
 real(dp) :: scaling_constant
 character(len=500) :: line, system, iomsg
 character(len=5) :: symbol
!arrays
 integer,allocatable :: nattyp(:)
 logical,allocatable :: duplicated(:)
 character(len=5),allocatable :: symbols(:), dupe_symbols(:)
 real(dp),allocatable :: xcart(:,:)

!************************************************************************

 ! Example of POSCAR (with 6 figures --> space group won't be recognized by Abinit
 ! See also https://github.com/ExpHP/vasp-poscar/blob/master/doc/format.md

 ! Mg1 B2
 ! 1.0
 ! 2.672554 1.543000 0.000000
 ! -2.672554 1.543000 0.000000
 ! 0.000000 0.000000 3.523000
 ! Mg B
 ! 1 2
 ! direct
 ! 0.000000 0.000000 0.000000 Mg
 ! 0.333333 0.666667 0.500000 B
 ! 0.666667 0.333333 0.500000 B

 new%fileformat = "poscar"
 read(unit, "(a)", err=10, iomsg=iomsg) new%title
 read(unit, *, err=10, iomsg=iomsg) scaling_constant
 do ii=1,3
   read(unit, *, err=10, iomsg=iomsg) new%rprimd(:, ii)
 end do

 ! Read line with the names of the atoms.
 read(unit, "(a)", err=10, iomsg=iomsg) line
 !print *, "line:", trim(line)

 new%ntypat = 0
 do ii=1,2
   if (ii == 2) then
     ABI_MALLOC(symbols, (new%ntypat))
   end if
   itypat = 0; beg = 1
   do
     ierr = next_token(line, beg, symbol)
     !print *, "ierr:", ierr, "beg:", beg, "symbol:", trim(symbol)
     if (ierr /= 0) exit
     if (ii == 1) new%ntypat = new%ntypat + 1
     if (ii == 2) then
       itypat = itypat + 1
       symbols(itypat) = trim(symbol)
     end if
   end do
 end do
 !write(std_out, *)"ntypat: ", new%ntypat, "symbols: ", symbols

 ! TODO: Handle case in which not all atoms are not grouped by type
 ABI_MALLOC(duplicated, (new%ntypat))
 duplicated = .False.
 do itypat=1,new%ntypat-1
   do ii=itypat+1, new%ntypat
     if (symbols(itypat) == symbols(ii)) duplicated(ii) = .True.
   end do
 end do

 ! number of atoms of each type.
 ! NOTE: Assuming ntypat == npsp thus alchemical mixing is not supported.
 ! There's a check in the main parser though.
 ABI_MALLOC(nattyp, (new%ntypat))
 read(unit, *, err=10, iomsg=iomsg) nattyp
 new%natom = sum(nattyp)
 ABI_FREE(nattyp)

 if (any(duplicated)) then
   ! Need to recompute ntypat and symbols taking into account duplication.
   MSG_WARNING("Found POSCAR with duplicated symbols")
   ABI_MOVE_ALLOC(symbols, dupe_symbols)
   new%ntypat = count(.not. duplicated)
   ABI_MALLOC(symbols, (new%ntypat))
   cnt = 0
   do ii=1,size(duplicated)
     if (.not. duplicated(ii)) then
       cnt = cnt + 1; symbols(cnt) = dupe_symbols(ii)
     end if
   end do
   ABI_FREE(dupe_symbols)
 end if

 ! At this point, we can allocate Abinit arrays.
 call new%malloc()

 ! Note that first letter should be capitalized, rest must be lower case
 do itypat=1,new%ntypat
   new%znucl(itypat) = symbol2znucl(symbols(itypat))
 end do

 read(unit, *, err=10, iomsg=iomsg) system
 system = tolower(system)
 if (system /= "cartesian" .and. system /= "direct") then
   MSG_ERROR(sjoin("Expecting `cartesian` or `direct` for the coordinate system but got:", system))
 end if

 ! Parse atomic positions.
 do iatom=1,new%natom

   ! This should implement the POSCAR format.
   read(unit, *, err=10, iomsg=iomsg) new%xred(:, iatom), symbol
   if (len_trim(symbol) == 0) then
     if (new%ntypat == 1) then
       MSG_COMMENT("POTCAR without element symbol after coords but this is not critical because ntypat == 1")
       symbol = symbols(1)
     else
       MSG_ERROR("POTCAR positions should be followed by element symbol.")
     end if
   end if

   ! This to handle symbol + oxidation state e.g. Li1+
   !print *, symbol
   ii = find_digit(symbol)
   if (ii /= 0) symbol = symbol(:ii-1)

   do itypat=1, new%ntypat
     if (symbols(itypat) == symbol) then
       new%typat(iatom) = itypat; exit
     end if
   end do
   if (itypat == new%ntypat + 1) then
     MSG_ERROR(sjoin("Cannot find symbol:`", symbol, " `in initial symbol list. Typo or POSCAR without symbols?."))
   end if
 end do

 ! Convert ang -> bohr
 if (scaling_constant > zero) then
   new%rprimd = scaling_constant * new%rprimd * Ang_Bohr
 else if (scaling_constant < zero) then
   ! A negative scale factor is treated as a volume. translate scaling_constant to a lattice vector scaling.
   new%rprimd = Ang_Bohr * new%rprimd * (-scaling_constant / abs(det3r(new%rprimd))) ** (one / three)
 else
   ABI_CHECK(scaling_constant > zero, sjoin("scaling constant must be /= 0 but found:", ftoa(scaling_constant)))
 end if

 if (system == "cartesian") then
   ! Go from cartesian to reduced.
   ABI_MALLOC(xcart, (3, new%natom))
   xcart = new%xred
   call xcart2xred(new%natom, new%rprimd, xcart, new%xred)
   ABI_FREE(xcart)
 end if

 ABI_FREE(symbols)
 ABI_FREE(duplicated)
 return

 10 MSG_ERROR(sjoin("Error while parsing POSCAR file,", ch10, "iomsg:", trim(iomsg)))

end function geo_from_poscar_unit
!!***

!!****f* m_parser/geo_print_abivars
!! NAME
!!  geo_print_abivars
!!
!! FUNCTION
!!  Print Abinit variables corresponding to POSCAR
!!
!! SOURCE

subroutine geo_print_abivars(self, unit)

!Arguments ------------------------------------
 class(geo_t),intent(in) :: self
 integer,intent(in) :: unit

!Local variables-------------------------------
 integer :: ii, iatom, itypat

!************************************************************************

 if (unit == dev_null) return

 write(unit, "(2a)")"# fileformat: ", trim(self%fileformat)
 if (len_trim(self%title) > 0) write(unit, "(2a)")"# ",trim(self%title)
 write(unit, "(a, i0)")" natom ", self%natom
 write(unit, "(a, i0)")" ntypat ", self%ntypat
 write(unit, sjoin("(a, ", itoa(self%natom), "(i0,1x))")) " typat ", self%typat
 write(unit, sjoin("(a, ", itoa(self%ntypat), "(f5.1,1x))")) " znucl ", self%znucl
 write(unit, "(a)")" acell 1 1 1 Bohr"
 write(unit, "(a)")" rprim "
 do ii=1,3
   write(unit, "(2x, 3(f11.7,1x))") self%rprimd(:, ii)
 end do
 write(unit, "(a)")" xred"
 do iatom=1,self%natom
   itypat = self%typat(iatom)
   write(unit, "(2x, 3(f11.7,1x),3x,2a)") self%xred(:, iatom) , " # ", trim(znucl2symbol(self%znucl(itypat)))
 end do

end subroutine geo_print_abivars
!!***

!!****f* m_parser/geo_from_netdf_path
!! NAME
!!  geo_from_netdf_path
!!
!! FUNCTION
!!
!! SOURCE

type(geo_t) function geo_from_netcdf_path(path, comm) result(new)

!Arguments ------------------------------------
 character(len=*),intent(in) :: path
 integer,intent(in) :: comm

!Local variables-------------------------------
 integer, parameter :: master = 0
 integer :: ncid, npsp, dimid, itime
 logical :: has_nimage

!************************************************************************

 new%fileformat = "netcdf"

#ifdef HAVE_NETCDF
 if (xmpi_comm_rank(comm) == master) then
   NCF_CHECK(nctk_open_read(ncid, path, xmpi_comm_self))

   if (endswith(path, "_HIST.nc")) then
     ! See def_file_hist.
     !MSG_ERROR("Cannot yet read structure from HIST.nc file")
     NCF_CHECK(nctk_get_dim(ncid, "natom", new%natom))
     NCF_CHECK(nctk_get_dim(ncid, "ntypat", new%ntypat))

     NCF_CHECK(nctk_get_dim(ncid, "npsp", npsp))
     ABI_CHECK(npsp == new%ntypat, 'Geo from HIST file with alchemical mixing!')
     has_nimage = nf90_inq_dimid(ncid, "nimage", dimid) == nf90_noerr
     ABI_CHECK(.not. has_nimage, "Cannot initialize structure from HIST.nc when file contains images.")

     call new%malloc()

     NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid, "typat"), new%typat))
     NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid, "znucl"), new%znucl))

     ! time is NF90_UNLIMITED
     NCF_CHECK(nctk_get_dim(ncid, "time", itime))

     ! dim3 = [xyz_id, xyz_id, time_id]
     NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid, "rprimd"), new%rprimd, start=[1,1,itime]))

     ! dim3 = [xyz_id, natom_id, time_id]
     NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid, "xred"), new%xred, start=[1,1,itime]))

   else
     ! Assume netcdf file produced by calling crystal%ncwrite
     NCF_CHECK(nctk_get_dim(ncid, "number_of_atoms", new%natom))
     NCF_CHECK(nctk_get_dim(ncid, "number_of_atom_species", new%ntypat))

     ! Test if alchemical. NB: nsps added in crystal_ncwrite in v9.
     if (nf90_inq_dimid(ncid, "number_of_pseudopotentials", dimid) == nf90_noerr) then
       NCF_CHECK(nf90_inquire_dimension(ncid, dimid, len=npsp))
       ABI_CHECK(npsp == new%ntypat, 'Geo from HIST file with alchemical mixing!')
     end if

     call new%malloc()

     NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid, "primitive_vectors"), new%rprimd))
     NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid, "atom_species"), new%typat))
     NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid, "atomic_numbers"), new%znucl))
     NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid, "reduced_atom_positions"), new%xred))
   end if

   NCF_CHECK(nf90_close(ncid))
 end if
#endif

 call new%bcast(master, comm)
 !call new%print_abivars(std_out)

end function geo_from_netcdf_path
!!***

!!****f* m_parser/geo_bcast
!! NAME
!!  geo_bcast
!!
!! FUNCTION
!!  Brodcast object
!!
!! SOURCE

subroutine geo_bcast(self, master, comm)

!Arguments ------------------------------------
 class(geo_t),intent(inout) :: self
 integer,intent(in) :: master, comm

!Local variables-------------------------------
 integer :: ierr, my_rank, list_int(2)

!************************************************************************

 if (xmpi_comm_size(comm) == 1) return
 my_rank = xmpi_comm_rank(comm)

 if (my_rank == master) list_int = [self%natom, self%ntypat]
 call xmpi_bcast(list_int, master, comm, ierr)

 if (my_rank /= master) then
   self%natom = list_int(1); self%ntypat = list_int(2)
   call self%malloc()
 end if

 call xmpi_bcast(self%rprimd, master, comm, ierr)
 call xmpi_bcast(self%xred, master, comm, ierr)
 call xmpi_bcast(self%typat, master, comm, ierr)
 call xmpi_bcast(self%znucl, master, comm, ierr)
 call xmpi_bcast(self%title, master, comm, ierr)
 call xmpi_bcast(self%fileformat, master, comm, ierr)

end subroutine geo_bcast
!!***

!!****f* m_parser/geo_malloc
!! NAME
!!  geo_malloc
!!
!! FUNCTION
!!  Allocate memory once %natom and %ntypat are know
!!
!! SOURCE

subroutine geo_malloc(self)

!Arguments ------------------------------------
 class(geo_t),intent(inout) :: self

!************************************************************************

 ABI_MALLOC(self%typat, (self%natom))
 ABI_MALLOC(self%xred, (3, self%natom))
 ABI_MALLOC(self%znucl, (self%ntypat))

end subroutine geo_malloc
!!***

!!****f* m_parser/geo_free
!! NAME
!!  geo_free
!!
!! FUNCTION
!!  Free memory.
!!
!! SOURCE

subroutine geo_free(self)

!Arguments ------------------------------------
 class(geo_t),intent(inout) :: self

!************************************************************************

 ABI_SFREE(self%typat)
 ABI_SFREE(self%xred)
 ABI_SFREE(self%znucl)

end subroutine geo_free
!!***

!!****f* m_parser/get_acell_rprim
!! NAME
!!  get_acell_rprim
!!
!! FUNCTION
!!  Get acell and rprim from string
!!
!! INPUTS
!! string*(*)=character string containing all the input data. Initialized previously in instrng.
!! jdtset=number of the dataset looked for
!! iimage= index of the current image
!! nimage=Number of images.
!! marr=dimension of the intarr and dprarr arrays, as declared in the calling subroutine.
!!
!! OUTPUT
!! acell(3)=length of primitive vectors
!! rprim(3,3)=dimensionless real space primitive translations
!!
!! FUNCTION
!!
!! SOURCE

subroutine get_acell_rprim(lenstr, string, jdtset, iimage, nimage, marr, acell, rprim)

!Arguments ------------------------------------
 integer,intent(in) :: lenstr, jdtset, iimage, nimage, marr
 character(len=*),intent(in) :: string
 real(dp),intent(out) :: acell(3)
 real(dp),intent(out) :: rprim(3,3)

!Local variables-------------------------------
 integer :: tacell, tangdeg, tread, trprim, mu
 real(dp) :: a2, aa, cc, cosang
 character(len=500) :: msg
!arrays
 integer,allocatable :: intarr(:)
 real(dp) :: angdeg(3)
 real(dp),allocatable :: dprarr(:)

!************************************************************************

 ABI_MALLOC(intarr, (marr))
 ABI_MALLOC(dprarr, (marr))

 acell(1:3) = one
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'acell',tacell,'LEN')
 if(tacell==1) acell(1:3)=dprarr(1:3)
 call intagm_img(acell,iimage,jdtset,lenstr,nimage,3,string,"acell",tacell,'LEN')

 ! Check that input length scales acell(3) are > 0
 do mu=1,3
   if(acell(mu) <= zero) then
     write(msg, '(a,i0,a, 1p,e14.6,4a)' )&
      'Length scale ',mu,' is input as acell: ',acell(mu),ch10,&
      'However, length scales must be > 0 ==> stop',ch10,&
      'Action: correct acell in input file.'
     MSG_ERROR(msg)
   end if
 end do

 ! Initialize rprim, or read the angles
 tread=0
 call intagm(dprarr,intarr,jdtset,marr,9,string(1:lenstr),'rprim',trprim,'DPR')
 if (trprim==1) rprim(:,:) = reshape( dprarr(1:9), [3, 3])
 call intagm_img(rprim,iimage,jdtset,lenstr,nimage,3,3,string,"rprim",trprim,'DPR')

 if(trprim==0)then
   ! If none of the rprim were read ...
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'angdeg',tangdeg,'DPR')
   angdeg(:)=dprarr(1:3)
   call intagm_img(angdeg,iimage,jdtset,lenstr,nimage,3,string,"angdeg",tangdeg,'DPR')

   if(tangdeg==1)then
     !call wrtout(std_out,' ingeo: use angdeg to generate rprim.')

     ! Check that input angles are positive
     do mu=1,3
       if(angdeg(mu)<=0.0_dp) then
         write(msg, '(a,i0,a,1p,e14.6,a,a,a,a)' )&
          'Angle number ',mu,' is input as angdeg: ',angdeg(mu),ch10,&
          'However, angles must be > 0 ==> stop',ch10,&
          'Action: correct angdeg in the input file.'
         MSG_ERROR(msg)
       end if
     end do

     ! Check that the sum of angles is smaller than 360 degrees
     if(angdeg(1)+angdeg(2)+angdeg(3)>=360.0_dp) then
       write(msg, '(a,a,a,es14.4,a,a,a)' )&
        'The sum of input angles (angdeg(1:3)) must be lower than 360 degrees',ch10,&
        'while it is: ',angdeg(1)+angdeg(2)+angdeg(3),'.',ch10,&
        'Action: correct angdeg in the input file.'
       MSG_ERROR(msg)
     end if

     if( abs(angdeg(1)-angdeg(2))<tol12 .and. &
         abs(angdeg(2)-angdeg(3))<tol12 .and. &
         abs(angdeg(1)-90._dp)+abs(angdeg(2)-90._dp)+abs(angdeg(3)-90._dp)>tol12 )then
       ! Treat the case of equal angles (except all right angles):
       ! generates trigonal symmetry wrt third axis
       cosang=cos(pi*angdeg(1)/180.0_dp)
       a2=2.0_dp/3.0_dp*(1.0_dp-cosang)
       aa=sqrt(a2)
       cc=sqrt(1.0_dp-a2)
       rprim(1,1)=aa        ; rprim(2,1)=0.0_dp                 ; rprim(3,1)=cc
       rprim(1,2)=-0.5_dp*aa ; rprim(2,2)= sqrt(3.0_dp)*0.5_dp*aa ; rprim(3,2)=cc
       rprim(1,3)=-0.5_dp*aa ; rprim(2,3)=-sqrt(3.0_dp)*0.5_dp*aa ; rprim(3,3)=cc
       ! write(std_out,*)' ingeo: angdeg=',angdeg(1:3), aa,cc=',aa,cc
     else
       ! Treat all the other cases
       rprim(:,:)=0.0_dp
       rprim(1,1)=1.0_dp
       rprim(1,2)=cos(pi*angdeg(3)/180.0_dp)
       rprim(2,2)=sin(pi*angdeg(3)/180.0_dp)
       rprim(1,3)=cos(pi*angdeg(2)/180.0_dp)
       rprim(2,3)=(cos(pi*angdeg(1)/180.0_dp)-rprim(1,2)*rprim(1,3))/rprim(2,2)
       rprim(3,3)=sqrt(1.0_dp-rprim(1,3)**2-rprim(2,3)**2)
     end if

   end if
 end if ! No problem if neither rprim nor angdeg are defined: use default rprim

 ABI_FREE(intarr)
 ABI_FREE(dprarr)

end subroutine get_acell_rprim
!!***

end module m_parser
!!***
