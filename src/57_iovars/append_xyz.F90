!{\src2tex{textfont=tt}}
!!****f* ABINIT/append_xyz
!! NAME
!! append_xyz
!!
!! FUNCTION
!! Translate the data from a xyz file (xyz_fname),
!! and add it at the end of the usual ABINIT input data string (string),
!! taking into account the dtset (dtset_char)
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group (MJV).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine append_xyz(dtset_char,lenstr,string,xyz_fname,strln)

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_atomdata

 use m_io_tools,       only : open_file

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'append_xyz'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

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
 character(len=500) :: message
 type(atomdata_t) :: atom
!arrays
 real(dp),allocatable :: xangst(:,:)
 integer, save :: atomspecies(200) = 0
 character(len=500), save :: znuclstring = ""
 character(len=2),allocatable :: elementtype(:)

!************************************************************************

 lenstr_new=lenstr

 if (dtset_char == "-1") then
!  write znucl
   lenstr_old=lenstr_new
   lenstr_new=lenstr_new+7+len_trim(znuclstring)+1
   string(lenstr_old+1:lenstr_new)=" ZNUCL"//blank//trim(znuclstring)//blank

!  write ntypat
   ntypat = sum(atomspecies)
   write(string20,'(i10)') ntypat
   lenstr_old=lenstr_new
   lenstr_new=lenstr_new+8+len_trim(string20)+1
   string(lenstr_old+1:lenstr_new)=" NTYPAT"//blank//trim(string20)//blank

   return
 end if

!open file with xyz data
 if (open_file(xyz_fname, message, newunit=unitxyz, status="unknown") /= 0) then
   MSG_ERROR(message)
 end if
 write(message, '(3a)')' importxyz : Opened file ',trim(xyz_fname),'; content stored in string_xyz'
 call wrtout(std_out,message,'COLL')

!check number of atoms is correct
 read(unitxyz,*) natom
 
 write(string5,'(i5)')natom
 lenstr_old=lenstr_new
 lenstr_new=lenstr_new+7+len_trim(dtset_char)+1+5
 string(lenstr_old+1:lenstr_new)=" _NATOM"//trim(dtset_char)//blank//string5

 ABI_ALLOCATE(xangst,(3,natom))
 ABI_ALLOCATE(elementtype,(natom))

!read dummy line
 read(unitxyz,*)

!read atomic types and positions
 do iatom = 1, natom
   read(unitxyz,*) elementtype(iatom), xangst(:,iatom)
!  extract znucl for each atom type
   call atomdata_from_symbol(atom,elementtype(iatom))
   znucl = atom%znucl
   if (znucl > 200) then
     write (message,'(5a)')&
&     'Error: found element beyond Z=200 ', ch10,&
&     'Solution: increase size of atomspecies in append_xyz', ch10
     MSG_ERROR(message)
   end if
!  found a new atom type 
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
 string(lenstr_old+1:lenstr_new)=" _XANGST"//trim(dtset_char)//blank

 do iatom=1,natom
   do mu=1,3
     write(string20,'(f20.12)')xangst(mu,iatom)
     lenstr_old=lenstr_new
     lenstr_new=lenstr_new+20
     string(lenstr_old+1:lenstr_new)=string20
   end do
 end do

 ABI_DEALLOCATE(elementtype)
 ABI_DEALLOCATE(xangst)

!Check the length of the string
 if(lenstr_new>strln)then
   write(message,'(3a)')&
&   'The maximal size of the input variable string has been exceeded.',ch10,&
&   'The use of a xyz file is more character-consuming than the usual input file. Sorry.'
   MSG_BUG(message)
 end if

!Update the length of the string
 lenstr=lenstr_new

end subroutine append_xyz
!!***

