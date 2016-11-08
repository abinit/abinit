!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_effective_potential_file
!! NAME
!! m_effective_potential_file
!!
!! FUNCTION
!! This  module contains all routine to get the effective potential from files
!! (XML or DDB)
!!
!! COPYRIGHT
!! Copyright (C) 2000-2015 ABINIT group (AM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_effective_potential_file
 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_xmpi

 use iso_c_binding
 use m_harmonics_terms
 use m_anharmonics_terms
 use m_effective_potential
 use m_crystal,        only : crystal_t, crystal_init, crystal_free
 use m_ifc
 use m_io_tools, only : open_file

 implicit none

 public :: effective_potential_file_getDim
 public :: effective_potential_file_getType
 public :: effective_potential_file_read
 public :: effective_potential_file_readDisplacement
 public :: effective_potential_file_readStrain

 private :: char_f2c
 private :: char_c2f
 private :: xml_getdims
 private :: xml2effpot
 private :: ddb2effpot

#ifndef HAVE_LIBXML
 private :: rdfromline
 private :: rmtabfromline
 private :: rdfromline_value
#endif

#if defined HAVE_LIBXML
 public :: effpot_xml_checkXML
 public :: effpot_xml_read
 public :: effpot_xml_getValue
 public :: effpot_xml_getAttribute
 public :: effpot_xml_getDims

 interface
   subroutine effpot_xml_read(filename,natom,&
&     ntypat,nrpt,nqpt,amu,atmfrc,cell,dynmat,elastic_constants,energy,&
&     epsilon_inf,ewald_atmfrc,internal_strain,phfrq,rprimd,qph1l,&
&     short_atmfrc,typat,xcart,zeff)&
&                          bind(C,name="effpot_xml_read")
     use iso_c_binding, only : C_CHAR,C_DOUBLE,C_INT
     integer(C_INT) :: natom,ntypat,nrpt,nqpt
     integer(C_INT) :: typat(natom)
     integer(C_INT) :: cell(3,nrpt)
     real(C_DOUBLE) :: energy
     real(C_DOUBLE) :: dynmat(2,3,natom,3,natom,nqpt)
     real(C_DOUBLE) :: internal_strain(6,natom,3)
     real(C_DOUBLE) :: phfrq(3*natom,nqpt),qph1l(3,nqpt)
     real(C_DOUBLE) :: atmfrc(2,3,natom,3,natom,nrpt)
     real(C_DOUBLE) :: short_atmfrc(2,3,natom,3,natom,nrpt)
     real(C_DOUBLE) :: ewald_atmfrc(2,3,natom,3,natom,nrpt)
     real(C_DOUBLE) :: amu(ntypat),rprimd(3,3),epsilon_inf(3,3)
     real(C_DOUBLE) :: zeff(3,3,natom)
     real(C_DOUBLE) :: elastic_constants(6,6),xcart(3,natom)
     character(kind=C_CHAR) :: filename(*)
   end subroutine effpot_xml_read
 end interface

 interface
   subroutine effpot_xml_getDims(filename,natom,ntypat,nqpt,nrpt)&
&                          bind(C,name="effpot_xml_getDims")
     use iso_c_binding, only : C_CHAR,C_INT
     integer(C_INT) :: natom,ntypat,nqpt,nrpt
     character(kind=C_CHAR) :: filename(*)
   end subroutine effpot_xml_getDims
 end interface

 interface
   subroutine effpot_xml_checkXML(filename,name_root) &
&                          bind(C,name="effpot_xml_checkXML")
     use iso_c_binding, only : C_CHAR
     character(kind=C_CHAR) :: filename(*),name_root(*)
   end subroutine effpot_xml_checkXML
 end interface

 interface
   subroutine effpot_xml_getValue(filename,name_value,value_result) &
 &                          bind(C,name="effpot_xml_getValue")
      use iso_c_binding, only : C_CHAR
      implicit none
      character(kind=C_CHAR) :: filename(*),name_value(*)
      character(kind=C_CHAR) :: value_result
    end subroutine effpot_xml_getValue
  end interface
 
 interface
   subroutine effpot_xml_getAttribute(filename,name_value,name_attribute) &
&                          bind(C,name="effpot_xml_getAttribute")
     use iso_c_binding, only : C_CHAR
     character(kind=C_CHAR) :: filename(*),name_value(*),name_attribute(*)
   end subroutine effpot_xml_getAttribute
 end interface
#endif

 integer,parameter :: XML_RECL = 50000
!!***

CONTAINS  !===========================================================================================


!!****f* m_effective_potential_file/rdfromline
!! NAME
!! rdfromline
!!
!! FUNCTION
!! Read the value of a keyword from a XML line
!! Same function than m_pawxmlps/paw_rdfromline.F90
!!
!! INPUTS
!!  keyword= keyword which value has to be read
!!  line= string from which the data are read (line from a XML)
!!
!! OUTPUT
!!  output= (string) value of the keyword
!!
!! PARENTS
!!      xml2effpot
!!
!! CHILDREN
!!      paw_rdfromline
!!
!! SOURCE

 subroutine rdfromline(keyword,line,output)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rdfromline'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
  character(len=*), intent(in) :: keyword,line
  character(len=*), intent(out) :: output
!Local variables ---------------------------------------
  character(len=len(line)) :: temp
  integer :: pos,pos2

! *********************************************************************

 output=""
 pos=index(line,trim(keyword))
 if (pos>0) then
   temp=line(pos+len_trim(keyword):len_trim(line))
   pos=index(temp,char(34))
   if (pos>0) then
     pos2=index(temp(pos+1:len_trim(temp)),char(34))
     if (pos2>0) then
       output=temp(pos+1:pos+pos2-1)
     end if
   end if
 end if

 end subroutine rdfromline
!!***


!!****f* m_effective_potential_file/rmtabfromline
!! NAME
!! rmtabfromline
!!
!! FUNCTION
!! Read remove tab from the begining of line
!!
!! INPUTS
!!  line= string from which the data are read (line from a XML)
!!
!! OUTPUT
!!  output= line without tab
!!
!! PARENTS
!!    xml2effpot
!!
!! CHILDREN
!!    rmtabfromline
!!
!! SOURCE

recursive subroutine rmtabfromline(line)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rmtabfromline'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
  character(len=*), intent(inout) :: line
!Local variables ---------------------------------------
  integer :: pos

! *********************************************************************

 pos=index(line,char(9))
 if (pos==1) then
   line = line(2:len_trim(line))//" "
   call rmtabfromline(line)
 end if

 end subroutine rmtabfromline
!!***



!!****f* m_effective_potential_file/rdfromline_value
!! NAME
!! rdfromline
!!
!! FUNCTION
!! Read the value of a keyword from a XML line
!!
!! INPUTS
!!  keyword= keyword which value has to be read
!!  line= string from which the data are read (line from a XML)
!!
!! OUTPUT
!!  output= (string) value of the keyword
!!
!! PARENTS
!!      xml2effpot
!!
!! CHILDREN
!!      paw_rdfromline
!!
!! SOURCE

 subroutine rdfromline_value(keyword,line,output)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rdfromline_value'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
  character(len=*), intent(in) :: keyword,line
  character(len=*), intent(out) :: output
!Local variables ---------------------------------------
  character(len=len(line)) :: temp
  integer :: pos,pos2

! *********************************************************************

 output=""
 pos=index(line,trim(keyword))
 if (pos==2) then
   pos=pos+len_trim(keyword)
   pos=pos+index(line(pos:len_trim(line)),char(62))
   temp=line(pos:len_trim(line))
   pos2=index(temp,char(60))
   if (pos2>0) then
     output=line(pos:pos+pos2-2)
   else
     output=line(pos:len_trim(line))
   end if
 else
   if(pos>2)then
     output=line(1:pos-3)
   end if
 end if
end subroutine rdfromline_value
!!***


!!****f* m_effpot_xml/char_f2c
!! NAME
!!  char_f_to_c
!!
!! FUNCTION
!! Helper function to convert a Fortran string to a C string
!! Based on a routine by Joseph M. Krahn
!!
!! INPUTS
!!  f_string=Fortran string
!!
!! OUTPUT
!!  c_string=C string
!!
!! SOURCE

function char_f2c(f_string) result(c_string)
!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'char_f2c'
!End of the abilint section

 character(len=*),intent(in) :: f_string
 character(kind=C_CHAR,len=1) :: c_string(len_trim(f_string)+1)
!Local variables -------------------------------
 integer :: ii,strlen
!! *************************************************************************
 strlen=len_trim(f_string)
 forall(ii=1:strlen)
   c_string(ii)=f_string(ii:ii)
 end forall
 c_string(strlen+1)=C_NULL_CHAR
end function char_f2c
!!***

!----------------------------------------------------------------------

!!****f* m_effpot_xml/char_c2f
!! NAME
!!  char_c_to_f
!!
!! FUNCTION
!! Helper function to convert a C string to a Fortran string
!! Based on a routine by Joseph M. Krahn
!!
!! INPUTS
!!  c_string=C string
!!
!! OUTPUT
!!  f_string=Fortran string
!!
!! PARENTS
!!      m_libpaw_libxc
!!
!! CHILDREN
!!
!! SOURCE

subroutine char_c2f(c_string,f_string)
!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'char_c2f'
!End of the abilint section

 character(kind=C_CHAR,len=1),intent(in) :: c_string(*)
 character(len=*),intent(out) :: f_string
!Local variables -------------------------------
 integer :: ii
!! *************************************************************************
 ii=1
 do while(c_string(ii)/=C_NULL_CHAR.and.ii<=len(f_string))
   f_string(ii:ii)=c_string(ii) ; ii=ii+1
 end do
 if (ii<len(f_string)) f_string(ii:)=' '
end subroutine char_c2f
!!***


!!****f* m_effective_potential_file/xml_getdims
!! NAME
!! xml_getdims
!!
!! FUNCTION
!! Open xml file of effective potentiel, then reads the variables that
!! must be known in order to dimension the arrays before complete reading
!!
!! INPUTS
!! character(len=*) filnam: name of input or output file
!!
!! OUTPUT
!! natom=number of atoms
!! ntypat=number of atom types
!! nrpt  =number of real space points used to integrate IFC 
!  nph1l =number of wavevectors for phonon 
!!
!! PARENTS
!!      multibinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine xml_getdims(filename,natom,ntypat,nph1l,nrpt)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xml_getdims'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

 !Arguments ------------------------------------
 !scalars
  character(len=fnlen),intent(in) :: filename
  integer, intent(out) :: natom,ntypat,nph1l,nrpt
 !arrays
 !Local variables-------------------------------
 !scalar
  integer :: nrpt2
  real :: itypat
  character(len=500) :: message
#ifndef HAVE_LIBXML
  integer :: funit = 1
  integer :: iatom
  logical :: endfile,found
  character (len=XML_RECL) :: line,readline
  character (len=XML_RECL) :: strg,strg1
#endif
  !arrays
#ifndef HAVE_LIBXML
  real,allocatable :: typat(:)
#endif

 ! *************************************************************************

!Open the atomicdata XML file for reading
 write(message,'(a,a,a,a,a)') ' xml_getdims :',&
&    ' Opening the file ',trim(filename),' to read dimensions',&
&    ' (before initialisation)'

 call wrtout(std_out,message,'COLL')

 natom = zero
 ntypat= zero
 nph1l = zero
 nrpt  = zero
 nrpt2 = zero
 itypat= zero

!Open the atomicdata XML file for reading

#if defined HAVE_LIBXML
!Read with libxml
 call effpot_xml_getDims(char_f2c(trim(filename)),natom,ntypat,nph1l,nrpt)
#else
!Read by hand

!Start a reading loop
 endfile=.false.
 found=.false.

 if (open_file(filename,message,unit=funit,form="formatted",status="old",&
&              action="read") /= 0) then
   MSG_ERROR(message)
 end if

!First parse to know the number of atoms 
 do while ((.not.endfile).and.(.not.found))
   read(funit,'(a)',err=10,end=10) readline
   call rmtabfromline(readline)
   line=adjustl(readline);goto 20
   10 line="";endfile=.true.
   20 continue

!Need test with char(9) because the old version of XML file
!from old script includes tarbulation at the begining of each line
   if (line(1:5)==char(60)//'atom') then
     natom=natom+1
     cycle
   end if

   if (line(1:21)==char(60)//'total_force_constant') then
     nrpt = nrpt+1
     cycle
   end if

   if (line(1:21)==char(60)//'local_force_constant') then
     nrpt2 = nrpt2+1
     cycle
   end if

   if (line(1:7)==char(60)//'qpoint') then
     nph1l =  nph1l+1
     cycle
   end if   
 end do
 
 if (nrpt2/=nrpt) then
   write(message, '(2a,I5,3a,I5,5a)' )ch10,&
&   ' WARNING: the number of total IFC  (',nrpt,') is not equal to the  ',ch10,&
&   '          the number of short range IFC (',nrpt2,') in ',filename,ch10,&
&   '          the missing ifc will be set to zero',ch10
   call wrtout(std_out,message,'COLL')
   nrpt = max(nrpt2,nrpt)
 end if

 
!second parse to get the number of typat
 ABI_ALLOCATE(typat,(natom))
 typat = zero
 iatom = 0

 rewind(funit)
!Start a reading loop
 endfile=.false.
 found=.false.

 do while ((.not.endfile).and.(.not.found))
   read(funit,'(a)',err=40,end=40) readline
   call rmtabfromline(readline)
   line=adjustl(readline);goto 50
   40 line="";endfile=.true.
   50 continue

   if (line(1:5)==char(60)//'atom') then
     iatom = iatom + 1 
     call rdfromline("mass",line,strg)
     strg1=trim(strg)
     read(unit=strg1,fmt=*) itypat
     if (.not.any(typat==itypat)) then
       ntypat= ntypat+1
     end if
     typat(iatom) = itypat
   end if
   
   if (line(1:6)==char(60)//'local') then
     found=.true.
   end if
   
 end do

 close(funit)
 ABI_DEALLOCATE(typat)

#endif

end subroutine xml_getdims
!!***

!!****f* m_effective_potential_file/xml2effpot
!! NAME
!! xml2effpot
!!
!! FUNCTION
!! Open xml file of effective potentiel, then reads the variables 
!! and store them in effective potentential type
!!
!! INPUTS
!! eff_pot = effective potential type
!! comm=MPI communicator
!! character(len=*) filnam: name of input or output file
!!
!! OUTPUT
!! eff_pot = effective potential type 
!!
!! PARENTS
!!      multibinit
!!
!! CHILDREN
!!
!! SOURCE

 subroutine xml2effpot(eff_pot,filename,comm)

 use m_atomdata
 use m_effective_potential, only : effective_potential_type
 use m_multibinit_dataset, only : multibinit_dataset_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xml2effpot'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
!End of the abilint section

 implicit none

 !Arguments ------------------------------------
 !scalars
 character(len=*),intent(in) :: filename
 integer, intent(in) :: comm
 !arrays
 type(effective_potential_type), intent(inout) :: eff_pot

 !Local variables-------------------------------
 !scalar
 integer :: ierr,ii,itypat,my_rank,msym,natom,nrpt,ntypat
 integer :: nph1l,npsp,nproc,nsym,space_group,timrev
 real(dp):: energy,ucvol
 character(len=500) :: message
 integer,parameter :: master=0
 logical :: has_anharmonics = .FALSE.
 logical :: iam_master
#ifndef HAVE_LIBXML
 integer :: funit = 1
 integer :: iatom,iamu,iph1l,irpt,mu,nu,voigt
 real(dp):: amu
 logical :: endfile,found,short_range,total_range
 character (len=XML_RECL) :: line,readline
 character (len=XML_RECL) :: strg,strg1
 logical :: has_straincoupling= .FALSE.
#endif
 !arrays
 integer,allocatable :: typat(:)
 integer,allocatable  :: symrel(:,:,:),symafm(:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3)
 real(dp) :: elastic_constants(6,6),epsilon_inf(3,3)
 real(dp),allocatable :: all_amu(:),dynmat(:,:,:,:,:,:)
 real(dp),allocatable :: internal_strain(:,:,:),phfrq(:,:),qph1l(:,:),tnons(:,:)
 real(dp),allocatable :: xcart(:,:),xred(:,:),zeff(:,:,:),znucl(:),zion(:)
 character(len=132),allocatable :: title(:)  
 type(ifc_type) :: ifcs
 type(ifc_type),dimension(:),allocatable :: phonon_strain
 type(crystal_t)  :: crystal
 type(atomdata_t) :: atom
#ifndef HAVE_LIBXML
 real(dp),allocatable :: work2(:,:)
#endif

! *************************************************************************


 !Open the atomicdata XML file for reading
 write(message,'(a,a)')' Opening the file ',filename

 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 iam_master = (my_rank == master)

!Get Dimention of system and allocation/initialisation of array
 call effective_potential_file_getDim(filename,natom,ntypat,nph1l,nrpt,0)
 gmet= zero; gprimd = zero; rmet = zero; rprimd = zero
 elastic_constants = zero; epsilon_inf = zero
 ABI_ALLOCATE(all_amu,(ntypat))
 ABI_ALLOCATE(ifcs%atmfrc,(2,3,natom,3,natom,nrpt))
 ABI_ALLOCATE(ifcs%cell,(3,nrpt))
 ABI_ALLOCATE(ifcs%short_atmfrc,(2,3,natom,3,natom,nrpt))
 ABI_ALLOCATE(ifcs%ewald_atmfrc,(2,3,natom,3,natom,nrpt))
 ABI_ALLOCATE(internal_strain,(6,natom,3))
 ABI_ALLOCATE(dynmat,(2,3,natom,3,natom,nph1l))
 ABI_ALLOCATE(typat,(natom))
 ABI_ALLOCATE(phfrq,(3*natom,nph1l))
 ABI_ALLOCATE(qph1l,(3,nph1l))
 ABI_ALLOCATE(xcart,(3,natom))
 ABI_ALLOCATE(xred,(3,natom))
 ABI_ALLOCATE(zeff,(3,3,natom))
 ABI_ALLOCATE(zion,(ntypat))
 ABI_ALLOCATE(znucl,(ntypat))

 ABI_DATATYPE_ALLOCATE(phonon_strain,(6))
 do ii = 1,6
   ABI_ALLOCATE(phonon_strain(ii)%atmfrc,(2,3,natom,3,natom,nrpt))
   ABI_ALLOCATE(phonon_strain(ii)%cell,(3,nrpt))
   phonon_strain(ii)%nrpt   = nrpt
   phonon_strain(ii)%atmfrc = zero
   phonon_strain(ii)%cell   = zero
 end do

 all_amu(:) = zero
 dynmat(:,:,:,:,:,:)  = zero
 ifcs%nrpt = nrpt
 ifcs%atmfrc(:,:,:,:,:,:)  = zero
 ifcs%cell(:,:)  = zero
 ifcs%ewald_atmfrc(:,:,:,:,:,:) = zero
 ifcs%short_atmfrc(:,:,:,:,:,:) = zero
 internal_strain(:,:,:) = zero
 phfrq = zero
 qph1l = zero
 xcart = zero
 zeff  = zero
 znucl = zero

 if(iam_master)then
!Open the atomicdata XML file for reading
#if defined HAVE_LIBXML

   write(message,'(a,a,a,a)')'-Reading the file ',trim(filename),&
&   ' with LibXML library'
  
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
  !Read with libxml librarie
   call effpot_xml_read(char_f2c(trim(filename)),natom,ntypat,nrpt,nph1l,all_amu,&
&                       ifcs%atmfrc,ifcs%cell,&
&                       dynmat,elastic_constants,energy,&
&                       epsilon_inf,ifcs%ewald_atmfrc,&
&                       internal_strain,phfrq,rprimd,&
&                       qph1l,ifcs%short_atmfrc,typat,&
&                       xcart,zeff)

!  convert atomic mass unit to znucl
   do itypat=1,ntypat
     do ii=1,103
       call atomdata_from_znucl(atom,real(ii,dp))
       if (abs((real(atom%amu,sp)-real(all_amu(itypat),sp))&
&         /real(all_amu(itypat),sp)*100)<0.1) then
         znucl(itypat) = atom%znucl
         exit
       end if
     end do
   end do

#else

! Read by hand
   write(message,'(a,a,a,a)')'-Reading the file ',trim(filename),&
& ' with Fortran'

   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   if (open_file(filename,message,unit=funit,form="formatted",&
&               status="old",action="read") /= 0) then
     MSG_ERROR(message)
   end if

!Start a reading loop in fortran
   rewind(unit=funit)
 
   endfile=.false.
   found=.false.

   iatom  = 1 
   iamu   = 1
   itypat = 1
   irpt   = 1
   iph1l  = 1
   amu    = zero
   short_range  = .false.
   total_range  = .false.

   do while ((.not.endfile).and.(.not.found))
     read(funit,'(a)',err=10,end=10) readline
     call rmtabfromline(readline)
     line=adjustl(readline);goto 20
10   line="";endfile=.true.
20   continue

     if (.not.has_straincoupling) then
       
       if ((line(1:7)=='<energy')) then
         call rdfromline_value('energy',line,strg)
         if (strg/="") then 
           strg1=trim(strg)
           read(strg1,*) energy
         else
           read(funit,'(a)',err=10,end=10) readline
           call rmtabfromline(readline)
           line=adjustl(readline)
           call rdfromline_value('energy',line,strg)
           if (strg/="") then 
             strg1=trim(strg)
             read(strg1,*) energy
           else
             strg1=trim(line)
             read(strg1,*) energy
           end if
         end if
         cycle
       end if
       
       if ((line(1:10)=='<unit_cell')) then
         call rdfromline_value('unit_cell',line,strg)
         if (strg/="") then 
           strg1=trim(strg)
           read(strg1,*) (rprimd(mu,1),mu=1,3)
           read(funit,*) (rprimd(mu,2),mu=1,3)
         else
           do nu=1,2
             read(funit,*) (rprimd(mu,nu),mu=1,3)
           end do
         end if
         read(funit,'(a)',err=10,end=10) readline
         call rmtabfromline(readline)
         line=adjustl(readline)
         call rdfromline_value('unit_cell',line,strg)
         if (strg/="") then 
           strg1=trim(strg)
           read(strg1,*) (rprimd(mu,3),mu=1,3)
         else
           strg1=trim(line)
           read(strg1,*) (rprimd(mu,3),mu=1,3)
         end if
         cycle
       end if
     
       if ((line(1:12)=='<epsilon_inf')) then
         call rdfromline_value('epsilon_inf',line,strg)
         if (strg/="") then 
           strg1=trim(strg)
           read(strg1,*) (epsilon_inf(mu,1),mu=1,3)
           read(funit,*) (epsilon_inf(mu,2),mu=1,3)
         else
           do nu=1,2
             read(funit,*) (epsilon_inf(mu,nu),mu=1,3)
           end do
         end if
         read(funit,'(a)',err=10,end=10) readline
         call rmtabfromline(readline)
         line=adjustl(readline)
         call rdfromline_value('epsilon_inf',line,strg)
         if (strg/="") then 
           strg1=trim(strg)
           read(strg1,*) (epsilon_inf(mu,3),mu=1,3)
         else
           strg1=trim(line)
           read(strg1,*) (epsilon_inf(mu,3),mu=1,3)
         end if
         cycle
       end  if
       
       if ((line(1:8)=='<elastic')) then
         call rdfromline_value('elastic',line,strg)
         if (strg/="") then 
           strg1=trim(strg)
           read(strg1,*) (elastic_constants(mu,1),mu=1,6)
           do nu=2,5
             read(funit,*) (elastic_constants(mu,nu),mu=1,6)
           end do
         else
           do nu=1,5
             read(funit,*) (elastic_constants(mu,nu),mu=1,6)
           end do
         end if
         read(funit,'(a)',err=10,end=10) readline
         call rmtabfromline(readline)
         line=adjustl(readline)
         call rdfromline_value('elastic',line,strg)
         if (strg/="") then 
           strg1=trim(strg)
           read(strg1,*) (elastic_constants(mu,6),mu=1,6)
         else
           strg1=trim(line)
           read(strg1,*) (elastic_constants(mu,6),mu=1,6)
         end if
         cycle
       end if
       
       if ((line(1:5)=='<atom')) then
         call rdfromline("mass",line,strg)
         strg1=trim(strg)
         read(strg1,*) amu
         if (.not.any(all_amu==amu)) then
           all_amu(iamu) = amu
           typat(iatom) = amu
           !convert atomic mass unit to znucl
           do ii=1,103
             call atomdata_from_znucl(atom,real(ii,dp))
             if (abs((real(atom%amu,sp)-real(amu,sp))&
               &           /real(amu,sp)*100)<0.1) then
               znucl(iamu) = atom%znucl
               exit
             end if
           end do
           iamu = iamu +1
         end if
         do itypat=1,ntypat
           if(amu==all_amu(itypat)) then
             typat(iatom) = itypat
           end if
         end do
         cycle
       end if
       
       if ((line(1:9)=='<position')) then
         call rdfromline_value('position',line,strg)
         if (strg/="") then 
           strg1=trim(strg)
           read(strg1,*)(xcart(mu,iatom),mu=1,3)
         else
           read(funit,'(a)',err=10,end=10) readline
           call rmtabfromline(readline)
           line=adjustl(readline)
           call rdfromline_value('position',line,strg)
           if (strg/="") then 
             strg1=trim(strg)
             read(strg1,*)(xcart(mu,iatom),mu=1,3)
           else
             strg1=trim(line)
             read(strg1,*)(xcart(mu,iatom),mu=1,3)
           end if
         end if
         cycle
       end if
       
       if ((line(1:11)=='<borncharge')) then
         call rdfromline_value('borncharge',line,strg)
         if (strg/="") then 
           strg1=trim(strg)
           read(strg1,*) (zeff(mu,1,iatom),mu=1,3)
           read(funit,*) (zeff(mu,2,iatom),mu=1,3)
         else
           do nu=1,2
             read(funit,*) (zeff(mu,nu,iatom),mu=1,3)
           end do
         end if
         read(funit,'(a)',err=10,end=10) readline
         line=adjustl(readline)
         call rdfromline_value('borncharge',line,strg)
         if (strg/="") then 
           strg1=trim(strg)
           read(strg1,*) (zeff(mu,3,iatom),mu=1,3)
         else
           strg1=trim(line)
           read(strg1,*) (zeff(mu,3,iatom),mu=1,3)
         end if
         cycle
       end if
       
       if ((line(1:7)==char(60)//char(47)//'atom'//char(62))) then
         iatom=iatom+1
         cycle
       end if
       
       if ((line(1:12)=='<local_force')) then
         read(funit,'(a)',err=10,end=10) readline
         call rmtabfromline(readline)
         line=adjustl(readline)
         call rdfromline_value('data',line,strg)
         if (strg/="") then 
           ABI_ALLOCATE(work2,(3*natom,3*natom))
           strg1=trim(strg)
           read(strg1,*) (work2(1,nu),nu=1,3*natom)
           do mu=2,3*natom-1
             read(funit,*)(work2(mu,nu),nu=1,3*natom)
           end do
           read(funit,'(a)',err=10,end=10) readline
           call rmtabfromline(readline)
           line=adjustl(readline)
           call rdfromline_value('data',line,strg)
           if (strg/="") then 
             strg1=trim(strg)
             read(strg1,*) (work2(3*natom,nu),nu=1,3*natom)
           else
             strg1=trim(line)
             read(strg1,*) (work2(3*natom,nu),nu=1,3*natom)
           end if
           ifcs%short_atmfrc(1,:,:,:,:,irpt) =&
&                           reshape(work2,(/3,natom,3,natom/))
           ABI_DEALLOCATE(work2)
         else
           ABI_ALLOCATE(work2,(3*natom,3*natom))
           do mu=1,3*natom
             read(funit,*)(work2(mu,nu),nu=1,3*natom)
           end do
           ifcs%short_atmfrc(1,:,:,:,:,irpt) = &
&                         reshape(work2,(/3,natom,3,natom/))
           ABI_DEALLOCATE(work2)
         end if
         short_range = .true.
       end if
       
       if ((line(1:12)=='<total_force')) then
         read(funit,'(a)',err=10,end=10) readline
         call rmtabfromline(readline) 
         line=adjustl(readline)
         call rdfromline_value('data',line,strg)
         if (strg/="") then 
           ABI_ALLOCATE(work2,(3*natom,3*natom))
           strg1=trim(strg)
           read(strg1,*) (work2(1,nu),nu=1,3*natom)
           do mu=2,3*natom-1
             read(funit,*)(work2(mu,nu),nu=1,3*natom)
           end do
           read(funit,'(a)',err=10,end=10) readline
           call rmtabfromline(readline)
           line=adjustl(readline)
           call rdfromline_value('data',line,strg)
           if (strg/="") then 
             strg1=trim(strg)
             read(strg1,*) (work2(3*natom,nu),nu=1,3*natom)
           else
             strg1=trim(line)
             read(strg1,*) (work2(3*natom,nu),nu=1,3*natom)
           end if
           ifcs%atmfrc(1,:,:,:,:,irpt) = reshape(work2,(/3,natom,3,natom/))
           ABI_DEALLOCATE(work2)
         else
           ABI_ALLOCATE(work2,(3*natom,3*natom))
           do mu=1,3*natom
             read(funit,*)(work2(mu,nu),nu=1,3*natom)
           end do
           ifcs%atmfrc(1,:,:,:,:,irpt) = reshape(work2,(/3,natom,3,natom/))
           ABI_DEALLOCATE(work2)
         end if
         total_range = .true.
       end if
       
       if ((line(1:5)=='<cell')) then
         call rdfromline_value('cell',line,strg)
         if (strg/="") then 
           strg1=trim(strg)
           read(strg1,*)(ifcs%cell(mu,irpt),mu=1,3)
         else
           read(funit,*)(ifcs%cell(mu,irpt),mu=1,3)
         end if
         if (total_range)then
           if(short_range)then
!            retrive short range part
             ifcs%ewald_atmfrc(1,:,:,:,:,irpt) = ifcs%atmfrc(1,:,:,:,:,irpt) &
&                                              - ifcs%short_atmfrc(1,:,:,:,:,irpt)  
           else
             ifcs%ewald_atmfrc(1,:,:,:,:,irpt) = ifcs%atmfrc(1,:,:,:,:,irpt)
             ifcs%short_atmfrc(1,:,:,:,:,irpt) = zero
           end if
           irpt = irpt + 1
           total_range = .false.
           short_range  = .false.
         end if
       end if
       
       if ((line(1:7)=='<qpoint')) then
         call rdfromline_value('qpoint',line,strg)
         if (strg/="") then 
           strg1=trim(strg)
           read(strg1,*)(qph1l(mu,iph1l),mu=1,3)
         else        
           read(funit,*) (qph1l(mu,iph1l),mu=1,3)
         end if
       end if
       
       if ((line(1:12)=='<frequencies')) then
         call rdfromline_value('frequencies',line,strg)
         if (strg/="") then 
           strg1=trim(strg)
           read(strg1,*)(phfrq(mu,iph1l),mu=1,3*natom)
         else
           do nu=1,natom
             read(funit,*) (phfrq(((nu-1)*3)+mu,iph1l),mu=1,3)
           end do
         end if
       end if
       
       if ((line(1:17)=='<dynamical_matrix')) then
         call rdfromline_value('dynamical_matrix',line,strg)
         if (strg/="") then 
           ABI_ALLOCATE(work2,(3*natom,3*natom))
           strg1=trim(strg)
           read(strg1,*) (work2(nu,1),nu=1,3*natom)
           do mu=2,3*natom-1
             read(funit,*)(work2(nu,mu),nu=1,3*natom)
           end do
           read(funit,'(a)',err=10,end=10) readline
           call rmtabfromline(readline)
           line=adjustl(readline)
           call rdfromline_value('dynamical_matrix',line,strg)
           if (strg/="") then 
             strg1=trim(strg)
             read(strg1,*) (work2(3*natom,nu),nu=1,3*natom)
           else
             strg1=trim(line)
             read(strg1,*) (work2(3*natom,nu),nu=1,3*natom)
           end if
           dynmat(1,:,:,:,:,iph1l) = reshape(work2,(/3,natom,3,natom/))
           ABI_DEALLOCATE(work2)
         else
           ABI_ALLOCATE(work2,(3*natom,3*natom))
           do mu=1,3*natom
             read(funit,*)(work2(nu,mu),nu=1,3*natom)
           end do
           dynmat(1,:,:,:,:,iph1l) = reshape(work2,(/3,natom,3,natom/))
           ABI_DEALLOCATE(work2)
         end if
       end if
       
       if ((line(1:8)==char(60)//char(47)//'phonon')) then
         iph1l = iph1l +1
       end if
       
       if ((line(1:16)=='<strain_coupling')) then
         read(funit,'(a)',err=10,end=10) readline
         call rdfromline("voigt",line,strg)
         strg1=trim(strg)
         read(strg1,*) voigt 
         has_straincoupling = .true.
         irpt = 1 
       end if
       
     else
!    Now treat the strain phonon coupling part 
       if ((line(1:16)=='<strain_coupling')) then
         read(funit,'(a)',err=10,end=10) readline
         call rdfromline("voigt",line,strg)
         strg1=trim(strg)
         read(strg1,*) voigt
         irpt = 1 
         cycle
       end if
       
       if ((line(1:22)=='<correction_force unit')) then
         call rdfromline_value('correction_force',line,strg)
         if (strg/="") then 
           ABI_ALLOCATE(work2,(natom,3))
           strg1=trim(strg)
           read(strg1,*) (work2(1,nu),nu=1,3)
           do mu=2,natom-1
             read(funit,*)(work2(mu,nu),nu=1,3)
           end do
           read(funit,'(a)',err=10,end=10) readline
           call rmtabfromline(readline)
           line=adjustl(readline)
           call rdfromline_value('correction_force',line,strg)
           if (strg/="") then 
             strg1=trim(strg)
             read(strg1,*) (work2(natom,nu),nu=1,3)
           else
             strg1=trim(line)
             read(strg1,*) (work2(natom,nu),nu=1,3)
           end if
           internal_strain(voigt+1,:,:) = work2(:,:)
           ABI_DEALLOCATE(work2)
         else
           ABI_ALLOCATE(work2,(natom,3))
           do mu=1,natom
             read(funit,*)(work2(mu,nu),nu=1,3)
           end do
           internal_strain(voigt+1,:,:) = work2(:,:)
           ABI_DEALLOCATE(work2)
         end if
       end if
       
       if ((line(1:26)=='<correction_force_constant')) then
         read(funit,'(a)',err=10,end=10) readline
         call rmtabfromline(readline)
         line=adjustl(readline)
         call rdfromline_value('data',line,strg)
         if (strg/="") then 
           ABI_ALLOCATE(work2,(3*natom,3*natom))
           strg1=trim(strg)
           read(strg1,*) (work2(1,nu),nu=1,3*natom)
           do mu=2,3*natom-1
             read(funit,*)(work2(mu,nu),nu=1,3*natom)
           end do
           read(funit,'(a)',err=10,end=10) readline
           call rmtabfromline(readline)
           line=adjustl(readline)
           call rdfromline_value('data',line,strg)
           if (strg/="") then 
             strg1=trim(strg)
             read(strg1,*) (work2(3*natom,nu),nu=1,3*natom)
           else
             strg1=trim(line)
             read(strg1,*) (work2(3*natom,nu),nu=1,3*natom)
           end if
           phonon_strain(voigt+1)%atmfrc(1,:,:,:,:,irpt) = &
&                          reshape(work2,(/3,natom,3,natom/))
           ABI_DEALLOCATE(work2)
         else
           ABI_ALLOCATE(work2,(3*natom,3*natom))
           do mu=1,3*natom
             read(funit,*)(work2(mu,nu),nu=1,3*natom)
           end do
           phonon_strain(voigt+1)%atmfrc(1,:,:,:,:,irpt) =&
&              reshape(work2,(/3,natom,3,natom/))
           ABI_DEALLOCATE(work2)
         end if
         has_anharmonics = .true.
         cycle
       end if
       
       if ((line(1:5)=='<cell')) then
         call rdfromline_value('cell',line,strg)
         if (strg/="") then 
           strg1=trim(strg)
           read(strg1,*)(phonon_strain(voigt+1)%cell(mu,irpt),mu=1,3)
         else
           read(funit,*)(phonon_strain(voigt+1)%cell(mu,irpt),mu=1,3)
         end if
         irpt = irpt + 1
         cycle
       end if
       
       if ((line(1:17)==char(60)//char(47)//'strain_coupling')) then
!      set nrpt for the previous value of strain
         phonon_strain(voigt+1)%nrpt = irpt - 1
!      restart the calculation of nrpt
       end if

     end if
   end do
   
!  Do some checks
   if (any(typat==zero)) then
     write(message, '(a,a,a)' )&
&       ' Unable to read the type of atoms ',trim(filename),ch10
     MSG_ERROR(message)
   end if
 
   if (any(znucl==zero)) then
     write(message, '(a,a,a)' )&
&       ' Unable to read the atomic number ',trim(filename),ch10
     MSG_ERROR(message)
   end if
 
   if (any(all_amu==zero)) then
     write(message, '(a,a,a)' )&
&     ' Unable to read the atomic mass ',trim(filename),ch10
     MSG_ERROR(message)
   end if
  
   close(unit=funit)

#endif

 end if !End if master

!MPI BROADCAST
 call xmpi_bcast(all_amu,master, comm, ierr)
 call xmpi_bcast(dynmat,master, comm, ierr)
 call xmpi_bcast(elastic_constants,master, comm, ierr)
 call xmpi_bcast(epsilon_inf,master, comm, ierr)
 call xmpi_bcast(ifcs%nrpt,master, comm, ierr)
 call xmpi_bcast(ifcs%atmfrc,master, comm, ierr)
 call xmpi_bcast(ifcs%cell,master, comm, ierr)
 call xmpi_bcast(ifcs%ewald_atmfrc,master, comm, ierr)
 call xmpi_bcast(ifcs%short_atmfrc,master, comm, ierr)
 call xmpi_bcast(internal_strain,master, comm, ierr)
 call xmpi_bcast(phfrq,master, comm, ierr)
 call xmpi_bcast(qph1l,master, comm, ierr)
 call xmpi_bcast(typat,master, comm, ierr)
 call xmpi_bcast(rprimd,master, comm, ierr)
 call xmpi_bcast(xcart,master, comm, ierr)
 call xmpi_bcast(zeff,master, comm, ierr)
 call xmpi_bcast(znucl,master, comm, ierr)

!Fill somes others variables
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 call xcart2xred(natom,rprimd,xcart,xred)

!Initialisation of crystal 
 msym = 1; npsp = ntypat; space_group = 0; timrev = 2
 nsym = 1
 ABI_ALLOCATE(symrel,(3,3,msym))
 ABI_ALLOCATE(symafm,(msym))
 ABI_ALLOCATE(tnons,(3,msym))
 symrel = zero; 
 symafm = zero;
 tnons = zero;
 do ii=1,3
   symrel(ii,ii,1)=one
 end do
 ABI_ALLOCATE(title, (ntypat))
 do ii=1,ntypat
   write(title(ii),'(a,i0)')"No title for typat ",ii
 end do

!Warning znucl is dimension with ntypat = nspsp hence alchemy is not supported here
 call crystal_init(all_amu,Crystal,space_group,natom,npsp,ntypat,nsym,rprimd,typat,xred,&
&  zion,znucl,timrev,.FALSE.,.FALSE.,title,&
&  symrel=symrel,tnons=tnons,symafm=symafm) 

!amu is not fill in crystal_init...
 Crystal%amu(:) = all_amu(:)

 ABI_DEALLOCATE(symrel)
 ABI_DEALLOCATE(symafm)
 ABI_DEALLOCATE(tnons)

!Initialisation of eff_pot
 call effective_potential_init(crystal,dynmat,energy,eff_pot,epsilon_inf,&
&                              elastic_constants,has_anharmonics,ifcs,internal_strain,phfrq,qph1l,&
&                              nph1l,zeff,comm,phonon_strain=phonon_strain)

!DEALLOCATION OF ARRAYS
 ABI_DEALLOCATE(all_amu)
 ABI_DEALLOCATE(ifcs%atmfrc)
 ABI_DEALLOCATE(ifcs%cell)
 ABI_DEALLOCATE(ifcs%short_atmfrc)
 ABI_DEALLOCATE(ifcs%ewald_atmfrc)
 ABI_DEALLOCATE(dynmat)
 ABI_DEALLOCATE(internal_strain)
 ABI_DEALLOCATE(phfrq)
 ABI_DEALLOCATE(qph1l)
 ABI_DEALLOCATE(title)
 ABI_DEALLOCATE(typat)
 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(xred)
 ABI_DEALLOCATE(zeff)
 ABI_DEALLOCATE(zion)
 ABI_DEALLOCATE(znucl)
 do ii = 1,6
   phonon_strain(ii)%nrpt   = nrpt
   phonon_strain(ii)%atmfrc = zero
   phonon_strain(ii)%cell   = zero
   ABI_DEALLOCATE(phonon_strain(ii)%atmfrc)
   ABI_DEALLOCATE(phonon_strain(ii)%cell)
 end do
 ABI_FREE(phonon_strain)

!DEALLOCATION OF TYPES
 call ifc_free(ifcs)
 call crystal_free(crystal)

end subroutine xml2effpot
!!***


!!****f* m_phonon_effective_potential/ddb2effpot
!!
!! NAME
!! ddb2effpot
!!
!! FUNCTION
!!  Transfert ddb into effective potential structure.
!!  Also calculate the IFC
!!
!! INPUTS
!! crystal  = number of atoms in primitive cell
!! ddb  = number of type of atoms
!! inp  = input of multibinit
!! comm=MPI communicator
!! OUTPUT
!! effective_potantial = effective_potential structure to be initialized
!!
!! PARENTS
!!    multibinit
!!
!! CHILDREN
!!    asria_calc,gtdyn9,ifc_init,ddb_internalstr,gtblk9,dfpt_phfrq,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine ddb2effpot(crystal,ddb, effective_potential,inp,comm)

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_dynmat
 use m_xmpi

 use m_ddb
 use m_ifc
 use m_copy,            only : alloc_copy
 use m_crystal,         only : crystal_t,crystal_print
 use m_dynmat,          only : asrif9,gtdyn9,canat9
 use m_multibinit_dataset, only : multibinit_dataset_type
 use m_effective_potential, only : effective_potential_type, effective_potential_free

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb2effpot'
 use interfaces_14_hidewrite
 use interfaces_72_response
 use interfaces_77_ddb
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
!arrays
 type(ddb_type),intent(inout) :: ddb
 type(effective_potential_type), intent(inout) :: effective_potential
 type(crystal_t),intent(in) :: crystal
 type(multibinit_dataset_type),intent(in) :: inp

!Local variables-------------------------------
!scalar
 real(dp):: icount1,icount2
 integer :: chneut,i1,i2,i3,ia,ib,iblok,idir1,idir2,ierr,ii,ipert1,iphl1
 integer :: ipert2,irpt,irpt2,ivarA,ivarB,max1,max2,max3,min1,min2,min3
 integer :: msize,mpert,natom,nblok,nrpt_new,rftyp,selectz
 integer :: my_rank,nproc
 logical :: iam_master=.FALSE.
 integer,parameter :: master=0
!arrays
 integer :: cell_number(3),cell2(3),rfelfd(4),rfphon(4),rfstrs(4)
 integer :: shift(3)
 real(dp):: dielt(3,3),elast_clamped(6,6),fact
 real(dp):: red(3,3),qphnrm(3),qphon(3,3)
 real(dp),allocatable :: blkval(:,:,:,:,:,:),d2asr(:,:,:,:,:)
 real(dp),allocatable :: instrain(:,:),zeff(:,:,:)
 real(dp),pointer :: atmfrc_red(:,:,:,:,:,:),cell_red(:,:),wghatm_red(:,:,:)
 character(len=500) :: message
 type(ifc_type) :: ifc
 real(dp),allocatable :: d2cart(:,:,:,:,:),displ(:)
 real(dp),allocatable :: eigval(:,:),eigvec(:,:,:,:,:),phfrq(:)
! *************************************************************************

!0 MPI variables
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 iam_master = (my_rank == master)

!Free the eff_pot before filling
 call effective_potential_free(effective_potential)

!Initialisation of usefull values  
  natom = ddb%natom
  nblok = ddb%nblok
  mpert=natom+6
  msize=3*mpert*3*mpert;

!Tranfert the ddb into usable array (ipert and idir format like in abinit)
  ABI_ALLOCATE(blkval,(2,3,mpert,3,mpert,nblok))
  blkval = zero
  blkval = reshape(ddb%val,(/2,3,mpert,3,mpert,nblok/))

!**********************************************************************
! Transfert basics values 
!**********************************************************************

  effective_potential%crystal%natom  = crystal%natom
  effective_potential%crystal%ntypat = crystal%ntypat
  effective_potential%crystal%rprimd = crystal%rprimd
  effective_potential%crystal%ucvol  = crystal%ucvol

  ABI_ALLOCATE(effective_potential%crystal%amu,(crystal%ntypat))
  effective_potential%crystal%amu(:)    = ddb%amu(:)

  ABI_ALLOCATE(effective_potential%crystal%typat,(crystal%natom))
  effective_potential%crystal%typat(:)  = crystal%typat(:)

  ABI_ALLOCATE(effective_potential%crystal%xcart,(3,crystal%natom))
  effective_potential%crystal%xcart(:,:)  = crystal%xcart(:,:)

  ABI_ALLOCATE(effective_potential%crystal%znucl,(crystal%ntypat))
  effective_potential%crystal%znucl(:)  = crystal%znucl(:)

!**********************************************************************
! Transfert energy from input file
!**********************************************************************
  write(message, '(2a,(80a),6a)') ch10,('=',ii=1,80),ch10,ch10,&
&     ' Extraction of the energy of the structure (unit: Hartree)',ch10
  call wrtout(std_out,message,'COLL')
  call wrtout(ab_out,message,'COLL')
  if(inp%energy_reference==zero)then
   write(message,'(6a)')ch10,&
&    ' Warning : Energy of the reference structure is not specify in',&
&    ' the input file.',ch10,' Energy will set to zero):',ch10
    call wrtout(std_out,message,'COLL')
    call wrtout(ab_out,message,'COLL')

  else
    effective_potential%energy = inp%energy_reference
    write(message,'(a,es25.12)') ' Energy = ',&
&                  effective_potential%energy
    call wrtout(std_out,message,'COLL')
    call wrtout(ab_out,message,'COLL')
  end if

!**********************************************************************
! Dielectric Tensor and Effective Charges
!**********************************************************************
  ABI_ALLOCATE(zeff,(3,3,natom))
  ABI_ALLOCATE(effective_potential%harmonics_terms%zeff,(3,3,natom))
 
  rftyp   = 1 ! Blocks obtained by a non-stationary formulation.
  chneut  = 1 ! The ASR for effective charges is imposed
  selectz = 0 ! No selection of some parts of the effective charge tensor
  iblok = ddb_get_dielt_zeff(ddb,crystal,rftyp,chneut,selectz,dielt,zeff)
  if (iblok /=0) then
    effective_potential%harmonics_terms%epsilon_inf = dielt
    effective_potential%harmonics_terms%zeff = zeff
  else
    effective_potential%harmonics_terms%epsilon_inf(1,1) = one 
    effective_potential%harmonics_terms%epsilon_inf(2,2) = one 
    effective_potential%harmonics_terms%epsilon_inf(3,3) = one 
    effective_potential%harmonics_terms%zeff = zero
  end if

!**********************************************************************
! Look after the blok no. that contains the stress tensor
!**********************************************************************
  write(message, '(a,a,(80a),a,a,a,a,a)') ch10,('=',ii=1,80),ch10,ch10,&
&   ' Extraction of the stress tensor (unit: GPa) and forces (unit: Ha/bohr)'
  call wrtout(std_out,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  ABI_ALLOCATE(effective_potential%forces,(3,natom))
  effective_potential%forces = zero
  effective_potential%internal_stress = zero

  qphon(:,1)=zero
  qphnrm(1)=zero
  rfphon(1:2)=0
  rfelfd(1:2)=0
  rfstrs(1:2)=0
  rftyp=4

  call gtblk9(ddb,iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)

  if (iblok /=0) then
!  firts give the corect stress values store in hartree
!  diagonal parts
   effective_potential%internal_stress(1)=blkval(1,1,natom+3,1,1,iblok)
   effective_potential%internal_stress(2)=blkval(1,2,natom+3,1,1,iblok)
   effective_potential%internal_stress(3)=blkval(1,3,natom+3,1,1,iblok)
!  the shear parts
   effective_potential%internal_stress(4)=blkval(1,1,natom+4,1,1,iblok)
   effective_potential%internal_stress(5)=blkval(1,2,natom+4,1,1,iblok)
   effective_potential%internal_stress(6)=blkval(1,3,natom+4,1,1,iblok)

!  Get forces
   effective_potential%forces(:,1:natom) = blkval(1,:,1:natom,1,1,iblok)

   write(message, '(3a)' )ch10,&
&   ' Cartesian components of forces (hartree/bohr)',ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   do ii = 1, natom
     write(message, '(I4,a,3(e16.8))' ) &
&     ii,'   ',effective_potential%forces(:,ii)

     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
   end do

   write(message, '(a,a)' )ch10,&
&   ' Cartesian components of stress tensor (hartree/bohr^3)'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '  sigma(1 1)=',effective_potential%internal_stress(1),&
&   '  sigma(3 2)=',effective_potential%internal_stress(4)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '  sigma(2 2)=',effective_potential%internal_stress(2),&
&   '  sigma(3 1)=',effective_potential%internal_stress(5)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '  sigma(3 3)=',effective_potential%internal_stress(3),&
&   '  sigma(2 1)=',effective_potential%internal_stress(6)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a)' ) ' '
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

 else
   
    write(message,'(2a)')ch10,&
&    ' Warning : Stress Tensor and forces are set to zero (not available in the DDB)'
    call wrtout(std_out,message,'COLL')
    call wrtout(ab_out,message,'COLL')
   
  end if

!**********************************************************************
! Elastic tensors at Gamma Point
!**********************************************************************
  write(message, '(a,a,(80a),a,a,a,a,a,a)') ch10,('=',ii=1,80),ch10,ch10,&
&   ' Extraction of the clamped elastic tensor (unit:10^2GPa)',ch10
  call wrtout(std_out,message,'COLL')
  call wrtout(ab_out,message,'COLL')

! look after the blok no.iblok that contains the elastic tensor
  qphon(:,1)=zero
  qphnrm(1)=zero
  rfphon(1:2)=0
  rfelfd(1:2)=0
  rfstrs(1:2)=3 ! Need uniaxial  both stresses and  shear stresses
  rftyp=1 ! Blocks obtained by a non-stationary formulation.
! for both diagonal and shear parts
  call gtblk9(ddb,iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)

  if (iblok /=0) then
!   extraction of the elastic constants from the blkvals (GPa)
    do ivarA=1,6
      do ivarB=1,6
!       because the elastic constant is 6*6,
!       so we should judge if the idir is larger than 3
!       or not
        if(ivarA>3) then
          idir1=ivarA-3
          ipert1=natom+4  !for the shear modulus
        else if(ivarA<=3) then
          idir1=ivarA
          ipert1=natom+3  !for the diagonal part
        end if
        if(ivarB>3) then
          idir2=ivarB-3
          ipert2=natom+4  !for the shear modulus
        else if(ivarB<=3) then
          idir2=ivarB
          ipert2=natom+3  !for the diagonal part
        end if
        elast_clamped(ivarA,ivarB) = blkval(1,idir1,ipert1,idir2,ipert2,iblok)
      end do
    end do
    fact=HaBohr3_GPa / crystal%ucvol
    do ivarA=1,6
      write(message,'(6f12.7)')elast_clamped(ivarA,1)*fact/100.00_dp,&
&                              elast_clamped(ivarA,2)*fact/100.00_dp,&
&                              elast_clamped(ivarA,3)*fact/100.00_dp,&
&                              elast_clamped(ivarA,4)*fact/100.00_dp,&
&                              elast_clamped(ivarA,5)*fact/100.00_dp,&
&                              elast_clamped(ivarA,6)*fact/100.00_dp
    call wrtout(std_out,message,'COLL')
    call wrtout(ab_out,message,'COLL')
    end do
    
!   Set the clamped tensor into the effective potentiel
    effective_potential%harmonics_terms%elastic_constants = elast_clamped

  else
    
    write(message,'(3a)')ch10,&
&    ' Warning : Elastic Tensor is set to zero (not available in the DDB)'
    call wrtout(std_out,message,'COLL')
    call wrtout(ab_out,message,'COLL')

!   Set the clamped tensor to zero into the effective potentiel (not available in the DDB)
    effective_potential%harmonics_terms%elastic_constants = zero
  end if

!**********************************************************************
!   Acoustic Sum Rule
!***************************************************************************
! ASR-correction (d2asr) has to be determined here from the Dynamical matrix at Gamma.
  ABI_ALLOCATE(d2asr,(2,3,natom,3,natom))

  write(message, '(a,a,(80a),a,a,a,a,a,a)') ch10,('=',ii=1,80),ch10,ch10,&
&   ' Calculation of acoustic sum rule',ch10
  call wrtout(std_out,message,'COLL')
  call wrtout(ab_out,message,'COLL')

! Find the Gamma block in the DDB (no need for E-field entries)
  qphon(:,1)=zero
  qphnrm(1)=zero
  rfphon(1:2)=1
  rfelfd(:)=0
  rfstrs(:)=0
  rftyp=inp%rfmeth

  call gtblk9(ddb,iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)
  
  d2asr = zero
  if (iblok /=0) then
    call asria_calc(inp%asr,d2asr,ddb%val(:,:,iblok),ddb%mpert,ddb%natom)
  end if


!**********************************************************************
! Interatomic Forces Calculation
!**********************************************************************
! ifc to be calculated for interpolation
  write(message, '(a,a,(80a),a,a,a,a)' ) ch10,('=',ii=1,80),ch10,ch10,&
&   ' Calculation of the interatomic forces ',ch10
  call wrtout(std_out,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  call ifc_init(ifc,crystal,ddb,inp%brav,inp%asr,inp%symdynmat,inp%dipdip,inp%rfmeth,&
&   inp%ngqpt(1:3),inp%nqshft,inp%q1shft,dielt,effective_potential%harmonics_terms%zeff,&
&   inp%nsphere,inp%rifcsph,inp%prtsrlr,inp%enunit,prtfreq=.True.)

!***************************************************************************
! Dynamical matrix calculation for each qpoint give in input (default:gamma)
!***************************************************************************

  ABI_ALLOCATE(d2cart,(2,3,mpert,3,mpert))
  ABI_ALLOCATE(displ,(2*3*natom*3*natom))
  ABI_ALLOCATE(eigval,(3,natom))
  ABI_ALLOCATE(eigvec,(2,3,natom,3,natom))
  ABI_ALLOCATE(phfrq,(3*natom))

  ABI_ALLOCATE(effective_potential%harmonics_terms%dynmat,(2,3,natom,3,natom,inp%nph1l))
  ABI_ALLOCATE(effective_potential%harmonics_terms%phfrq,(3*natom,inp%nph1l))
  ABI_ALLOCATE(effective_potential%harmonics_terms%qpoints,(3,inp%nph1l))

  write(message,'(a,(80a),3a)')ch10,('=',ii=1,80),ch10,ch10,&
&     ' Calculation of dynamical matrix for each ph1l points '
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')
  
!Transfer value in effective_potential structure
  effective_potential%harmonics_terms%nqpt      = inp%nph1l
  effective_potential%harmonics_terms%qpoints(:,:) = inp%qph1l(:,:)
  
  do iphl1=1,inp%nph1l

   ! Initialisation of the phonon wavevector
    qphon(:,1)=inp%qph1l(:,iphl1)
    if (inp%nph1l /= 0) qphnrm(1) = inp%qnrml1(iphl1)

    ! Get d2cart using the interatomic forces and the
    ! long-range coulomb interaction through Ewald summation
    call gtdyn9(ddb%acell,ifc%atmfrc,ifc%dielt,ifc%dipdip,ifc%dyewq0,d2cart,crystal%gmet,&
&     ddb%gprim,mpert,natom,ifc%nrpt,qphnrm(1),qphon(:,1),crystal%rmet,ddb%rprim,ifc%rpt,&
&     ifc%trans,crystal%ucvol,ifc%wghatm,crystal%xred,zeff)

    ! Calculation of the eigenvectors and eigenvalues of the dynamical matrix
    call dfpt_phfrq(ddb%amu,displ,d2cart,eigval,eigvec,crystal%indsym,&
&     mpert,crystal%nsym,natom,crystal%nsym,crystal%ntypat,phfrq,qphnrm(1),qphon,&
&     crystal%rprimd,inp%symdynmat,crystal%symrel,crystal%symafm,crystal%typat,crystal%ucvol)

    ! Write the phonon frequencies
    call dfpt_prtph(displ,inp%eivec,inp%enunit,ab_out,natom,phfrq,qphnrm(1),qphon)
    
    effective_potential%harmonics_terms%dynmat(:,:,:,:,:,iphl1) = d2cart(:,:,:natom,:,:natom)
    effective_potential%harmonics_terms%phfrq(:,iphl1) = phfrq(:) * Ha_cmm1
    
  end do
  
  ABI_DEALLOCATE(d2cart)
  ABI_DEALLOCATE(displ)
  ABI_DEALLOCATE(eigval)
  ABI_DEALLOCATE(eigvec)
  ABI_DEALLOCATE(phfrq)

!**********************************************************************
! Transfert inter-atomic forces constants in reduced coordinates
!**********************************************************************

!Reorder cell from canonical coordinates to reduced coordinates (for multibinit)
!store the number of ifc before rearrangement  

!Set the maximum and the miminum for the bound of the cell
  max1 = maxval(ifc%cell(1,:));  min1 = minval(ifc%cell(1,:))
  max2 = maxval(ifc%cell(2,:));  min2 = minval(ifc%cell(2,:))
  max3 = maxval(ifc%cell(3,:));  min3 = minval(ifc%cell(3,:))
  cell_number(1) = max1 - min1 +1
  cell_number(2) = max2 - min2 +1
  cell_number(3) = max3 - min3 +1
! set the new number of cell, sometimes, in canonical coordinates,
! somme cell are delete but they exist in reduced coordinates.
  nrpt_new = product(cell_number(:))

! Allocate temporary array
  ABI_ALLOCATE(atmfrc_red,(2,3,natom,3,natom,nrpt_new))
  ABI_ALLOCATE(wghatm_red,(natom,natom,nrpt_new))
  ABI_ALLOCATE(cell_red,(3,nrpt_new))
  
  wghatm_red(:,:,:) = zero
  atmfrc_red(:,:,:,:,:,:) =  zero
  
  if(iam_master)then
    icount1 = 0
    do irpt=1,ifc%nrpt
      icount1 = icount1 + sum(ifc%wghatm(:,:,irpt))
    end do

    do ia=1,natom
      do ib=1,natom

!      Simple Lattice
        if (inp%brav==1) then
!      In this case, it is better to work in reduced coordinates
!      As rcan is in canonical coordinates, => multiplication by gprim
          do ii=1,3
            red(1,ii)=  ifc%rcan(1,ia)*ddb%gprim(1,ii) + &
&                       ifc%rcan(2,ia)*ddb%gprim(2,ii) + &
&                       ifc%rcan(3,ia)*ddb%gprim(3,ii)
            red(2,ii)=  ifc%rcan(1,ib)*ddb%gprim(1,ii) + &
                        ifc%rcan(2,ib)*ddb%gprim(2,ii) + &
&                       ifc%rcan(3,ib)*ddb%gprim(3,ii)
          end do
        end if

!       Get the shift of cell
        shift(:) = anint(red(2,:) - crystal%xred(:,ib)) - anint(red(1,:) - crystal%xred(:,ia))

        do irpt=1,ifc%nrpt
       
          cell2(:)= int(ifc%cell(:,irpt) + shift(:))
       
!         Use boundary condition to get the right cell
          if (cell2(1) < min1 .and. cell2(1) < max1) then
            cell2(1) = cell2(1) + cell_number(1)
          else if (cell2(1) > min1 .and. cell2(1) > max1) then
            cell2(1) = cell2(1) - cell_number(1)
          end if

          if (cell2(2) < min2 .and. cell2(2) < max2) then
            cell2(2) = cell2(2) + cell_number(2)
          else if (cell2(2) > min2 .and. cell2(2) > max2) then
            cell2(2) = cell2(2) - cell_number(2)
          end if

          if (cell2(3) < min3 .and. cell2(3) < max3) then
            cell2(3) = cell2(3) + cell_number(3)
          else if (cell2(3) > min3 .and. cell2(3) > max3) then
            cell2(3) = cell2(3) - cell_number(3)
          end if

          irpt2=1
          do i1=min1,max1
            do i2=min2,max2
              do i3=min3,max3
                if (i1  ==  cell2(1)  .and.&
                    i2  ==  cell2(2)  .and.&
                    i3  ==  cell2(3)) then
                  wghatm_red(ia,ib,irpt2) =  ifc%wghatm(ia,ib,irpt)
                  atmfrc_red(:,:,ia,:,ib,irpt2) = ifc%atmfrc(:,:,ia,:,ib,irpt)
                  cell_red(1,irpt2) = i1
                  cell_red(2,irpt2) = i2
                  cell_red(3,irpt2) = i3
                end if
                irpt2 = irpt2 + 1
              end do
            end do
          end do
        end do
      end do
    end do

    irpt2 = 0
    icount2 = 0
    do irpt = 1, nrpt_new
      icount2 = icount2 + sum(wghatm_red(:,:,irpt))
      if (sum(wghatm_red(:,:,irpt)) /= 0) then
        irpt2 = irpt2 + 1
      end if
    end do

    if (icount1 /= icount2) then
      write(message,'(2a,I7,a,I7,a)')'The total wghatm is no more the same',ch10,&
&                        icount1,' before and ', icount2, ' now.'
      MSG_BUG(message)
    end if
 
!  Apply weight on each R point
    do irpt=1,nrpt_new
      do ia=1,effective_potential%crystal%natom 
        do ib=1,effective_potential%crystal%natom 
          atmfrc_red(:,:,ia,:,ib,irpt) = atmfrc_red(:,:,ia,:,ib,irpt)*wghatm_red(ia,ib,irpt) 
        end do
      end do
    end do
  end if

  call xmpi_bcast(atmfrc_red,master, comm, ierr)
  call xmpi_bcast(wghatm_red,master, comm, ierr)
  call xmpi_bcast(cell_red,master, comm, ierr)
  
!  Copy ifc into effective potential
! !!Warning eff_pot%ifcs only contains atmfrc,short_atmfrc,ewald_atmfrc,,nrpt and cell!!
! rcan,ifc%rpt,wghatm and other quantities 
! are not needed for effective potential!!!
  call ifc_free(ifc)
  call ifc_free(effective_potential%harmonics_terms%ifcs)

  effective_potential%harmonics_terms%ifcs%nrpt = nrpt_new

  ABI_ALLOCATE(effective_potential%harmonics_terms%ifcs%atmfrc,(2,3,natom,3,natom,nrpt_new))
  ABI_ALLOCATE(effective_potential%harmonics_terms%ifcs%short_atmfrc,(2,3,natom,3,natom,nrpt_new))
  ABI_ALLOCATE(effective_potential%harmonics_terms%ifcs%ewald_atmfrc,(2,3,natom,3,natom,nrpt_new))
  ABI_ALLOCATE(effective_potential%harmonics_terms%ifcs%cell,(3,nrpt_new))
  
  effective_potential%harmonics_terms%ifcs%nrpt = nrpt_new
  effective_potential%harmonics_terms%ifcs%cell(:,:) = cell_red(:,:)
  effective_potential%harmonics_terms%ifcs%atmfrc(:,:,:,:,:,:) = atmfrc_red(:,:,:,:,:,:)
  if (inp%dipdip == 1) then
    effective_potential%harmonics_terms%ifcs%short_atmfrc(:,:,:,:,:,:) = atmfrc_red(:,:,:,:,:,:)
  else
    effective_potential%harmonics_terms%ifcs%short_atmfrc(:,:,:,:,:,:) = zero
  end if
  effective_potential%harmonics_terms%ifcs%short_atmfrc(:,:,:,:,:,:) = atmfrc_red(:,:,:,:,:,:)
  effective_potential%harmonics_terms%ifcs%ewald_atmfrc(:,:,:,:,:,:) = zero
 
  
  ABI_DEALLOCATE(atmfrc_red)
  ABI_DEALLOCATE(wghatm_red)
  ABI_DEALLOCATE(cell_red)
  

!**********************************************************************
! Internal strain tensors at Gamma point
!**********************************************************************
  write(message, '(a,a,(80a),a,a,a)') ch10,('=',ii=1,80),ch10,ch10,&
&   ' Calculation of the internal-strain  tensor'
  call wrtout(std_out,message,'COLL')
  call wrtout(ab_out,message,'COLL')
  ABI_ALLOCATE(instrain,(3*natom,6))
! looking after the no. of blok that contains the internal strain tensor
  qphon(:,1)=zero
  qphnrm(1)=zero
  rfphon(1:2)=0
  rfelfd(1:2)=0
  rfstrs(1:2)=3
  rftyp=1
  call gtblk9(ddb,iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)
  
  ABI_ALLOCATE(effective_potential%harmonics_terms%internal_strain,(6,natom,3))
  effective_potential%harmonics_terms%internal_strain = zero

  if (iblok /=0) then

!    then print the internal stain tensor
    call ddb_internalstr(inp%asr,ddb%val,d2asr,iblok,instrain,ab_out,mpert,natom,nblok)

    do ipert1=1,6
      do ipert2=1,natom
        do idir2=1,3
          ii=3*(ipert2-1)+idir2
          effective_potential%harmonics_terms%internal_strain(ipert1,ipert2,idir2)=instrain(ii,ipert1)
        end do
      end do
    end do
  else
    write(message,'(3a)')ch10,&
&    ' Warning : Internal strain is set to zero (not available in the DDB)'
    call wrtout(std_out,message,'COLL')
    call wrtout(ab_out,message,'COLL')
  end if
!-------------------------------------------------------------------------------------
! DEALLOCATION OF ARRAYS
  ABI_DEALLOCATE(blkval)
  ABI_DEALLOCATE(zeff)
  ABI_DEALLOCATE(instrain)
  ABI_DEALLOCATE(d2asr)

  write(message,'(a)')ch10
  call wrtout(std_out,message,'COLL')
  call wrtout(ab_out,message,'COLL')

end subroutine ddb2effpot
!!***

!!****f* m_phonon_effective_potential/effective_potential_file_getType
!!
!! NAME
!! effective_potential_file_getType
!!
!! FUNCTION
!! This routine test the xml or ddb file
!!
!! INPUTS
!! filename = names of the files
!!
!! OUTPUT
!! type_file  = 0 no type found
!!              1 DDB file
!!              2 XML file
!!
!! PARENTS
!!   multibinit
!!
!! CHILDREN
!!   
!!
!! SOURCE

subroutine effective_potential_file_getType(filename,filetype)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_file_getType'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=fnlen),intent(in) :: filename
 integer, intent(out) :: filetype
!arrays

!Local variables-------------------------------
!scalar
 integer :: ddbun = 666
 character(len=500) :: message
 character (len=1000) :: line,readline
!arrays
! *************************************************************************

 filetype = 0

 if (open_file(filename,message,unit=ddbun,form="formatted",status="old",action="read") /= 0) then
   MSG_ERROR(message)
 end if

!Check if the file is a XML file or a DDB and in this case, store the DDB code.
 read(ddbun,'(a)') readline
 line=adjustl(readline)
 if(line(3:13)=="xml version") then
   filetype = 2
 else
   read(ddbun,'(a)') readline
   line=adjustl(readline)
   if(line(6:24)=="DERIVATIVE DATABASE") then
     filetype = 1
   end if
 end if
 
 close(ddbun)

end subroutine effective_potential_file_getType
!!***

!!****f* m_phonon_effective_potential/get_natom_from_file
!!
!! NAME
!! get_natom_from_file
!!
!! FUNCTION
!! This routine test the xml or ddb file
!! Return the number of atoms/ntypat in the unit cell from ddb and xml
!! Return nqpt ans nrpt if the file is XML file 
!! In case of DDB file, you have to run bigbx9 to get nrpt
!!
!! INPUTS
!! filename = names of the files
!! comm=MPI communicator
!!
!! OUTPUT
!! natom = number of atoms
!!
!! PARENTS
!!   multibinit
!!
!! CHILDREN
!!   ddb_getdims,xml_getdims,effective_potential_file_getType
!!
!! SOURCE

subroutine effective_potential_file_getDim(filename,natom,ntypat,nqpt,nrpt,comm)

 use m_ddb

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_file_getDim'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=fnlen),intent(in) :: filename
 integer,intent(in) :: comm
 integer,intent(out) :: natom,ntypat,nqpt,nrpt
!arrays

!Local variables-------------------------------
 !scalar
 integer,parameter :: vrsio8=100401,vrsio8_old=010929,vrsio8_old_old=990527
 integer :: dimekb,filetype,lmnmax,mband,mtyp,msym,nblok,nkpt,usepaw
 integer :: ddbun = 666
 character(len=500) :: message
!arrays
! *************************************************************************

 natom = 0
 ntypat= 0
 nqpt = 0
 nrpt  = 0

 call effective_potential_file_getType(filename,filetype)

 if(filetype==1)then
   write(message, '(6a)' )ch10,' The file ',trim(filename),ch10,&
&                  ' is DDB file (extraction of the number of atoms)',ch10
   call wrtout(std_out,message,'COLL')

   write(message, '(8a)' )&
&   ' The file ',trim(filename),ch10,' is ddb file only the number of atoms is read,',&
&    'if you want to predic the number of cell (nrpt)',ch10,' use bigbx9 routines',ch10
   call wrtout(std_out,message,'COLL')

   call ddb_getdims(dimekb,filename,lmnmax,mband,mtyp,msym,natom,nblok,&
&                  nkpt,ntypat,ddbun,usepaw,DDB_VERSION,comm)

   write(message, '(a,a,a,a)' )&
&   ' WARNING: Unable to read the number of cell (nrpt) in ddb file, nrpt is set to 0',ch10
   call wrtout(std_out,message,'COLL')

   write(message, '(a,a,a,a)' )&
&   ' WARNING: Unable to read the number of qpoint (nqpt) in ddb file (not implemented)',ch10
   call wrtout(std_out,message,'COLL')

!  Must read some value to initialze  array (nprt for ifc)
!   call bigbx9(inp%brav,dummy_cell,0,1,inp%ngqpt,inp%nqshft,nrpt,ddb%rprim,dummy_rpt)

 else if (filetype==2) then
   write(message, '(a,a,a,a,a)' )ch10,' The file ',trim(filename),&
&                ' is XML file (extraction of all informations)',ch10
   call wrtout(std_out,message,'COLL')
   
   call xml_getdims(filename,natom,ntypat,nqpt,nrpt)

 else
   write(message, '(a,a,a,a)' )&
&   ' The file ',trim(filename),' is not compatible with multibinit',ch10
   MSG_ERROR(message)
 end if

! Do some checks
 if (natom < 1) then
   write(message, '(a,a,a,a,a)' )&
&   ' Unable to read the number of atom from ',trim(filename),ch10,&
&   'This file  is not compatible with multibinit',ch10
   MSG_ERROR(message)
 end if
 
 if (filetype==2) then 

   if (natom < 1) then
     write(message, '(a,a,a)' )&
&     ' Unable to read the number of atom from ',trim(filename),ch10
     MSG_ERROR(message)
   end if

   if (nrpt < 1) then
     write(message, '(a,a,a)' )&
&     ' Unable to read the number of rpt points ',trim(filename),ch10
     MSG_ERROR(message)
   end if

   if (ntypat < 1) then
     write(message, '(a,a,a)' )&
&     ' Unable to read the number of type of atoms ',trim(filename),ch10
     MSG_ERROR(message)
   end if

 end if

end subroutine effective_potential_file_getDim
!!***

!****f* m_effective_potential/effective_potential_file_read
!!
!! NAME
!! effective_potential_file_read
!!
!! FUNCTION
!! tranfert file (XML or DDB) in effective potential type
!!
!! INPUTS
!! comm=MPI communicator
!! filename = name of the file
!!
!! OUTPUT
!! eff_pot = supercell structure with data to be output
!!
!! PARENTS
!!   multibinit
!!
!! CHILDREN
!!  bigbx9,ddb_from_file,effective_potential_file_getType,effective_potential_file_getDim,
!!  effective_potential_free, effective_potential_init,
!!  effective_potential_print,xml2effpot,wrtout
!! SOURCE
 
subroutine effective_potential_file_read(filename,eff_pot,inp,comm)

  use m_effective_potential
  use m_multibinit_dataset
  use m_ddb, only : ddb_from_file,ddb_free
  use m_strain
  use m_crystal, only : crystal_t, crystal_free
  use m_dynmat, only : bigbx9

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_file_read'
 use interfaces_14_hidewrite
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer,intent(in) :: comm
  character(len=fnlen),intent(in) :: filename
!array
  type(effective_potential_type), intent(inout)  :: eff_pot
  type(multibinit_dataset_type),optional,intent(in) :: inp
  type(ddb_type) :: ddb
  type(crystal_t) :: Crystal
!Local variables------------------------------
!scalars
  integer :: filetype,natom,ntypat,nqpt,nrpt
  character(500) :: message
!array
  integer,allocatable :: atifc(:)
  integer :: dummy_cell(3)
  real(dp) :: dummy_rpt(3)

! *************************************************************************

  call effective_potential_file_getType(filename,filetype)

  if (filetype/=0) then 

    if(filetype ==1) then
      if (.not.(present(inp))) then
        write(message, '(4a)' )&
&        ' effective_potential_file_read: you need to give input file to compute ',&
&        'the response fonction from DDB file ',ch10
        MSG_ERROR(message)
      end if

!     Read the DDB information, also perform some checks, and symmetrize partially the DDB
      write(message, '(3a)' )' Read the DDB information of the reference',&
 &      ' system and perform some checks',ch10
      call wrtout(std_out,message,'COLL')
      call wrtout(ab_out,message,'COLL')

      call effective_potential_file_getDim(filename,natom,ntypat,nqpt,nrpt,comm)!

!     In anaddb, inp%atifc is set to 1 2 3 ... natom (see anaddb help).
!     Then in the next routine inp%atifc is convert with 0 or 1 (inp%atifc is now 1 1 1 1 0).
!     In multibinit the conversion is done directly in m_multibinit_dataset.
!     So in the next routine, we set natifc to 0 to ignore the conversion.
!     To keep the intent(in) of the inp parameters, we need to use local variables:
      ABI_ALLOCATE(atifc,(inp%natom))
      atifc = inp%atifc

      call ddb_from_file(ddb,filename,inp%brav,natom,0,atifc,Crystal,comm)

!     And finaly, we can check if the value of atifc is not change...
      if (.not.all(atifc.EQ.inp%atifc)) then
        write(message, '(3a)' )&
&        ' effective_potential_file_read: problem with atifc input variables ',&
&        'in ddb_from_file',ch10
        MSG_BUG(message)
      end if
      ABI_DEALLOCATE(atifc)
!     Must read some value to initialze  array (nprt for ifc)
      call bigbx9(inp%brav,dummy_cell,0,1,inp%ngqpt,inp%nqshft,nrpt,ddb%rprim,dummy_rpt)

!     Transfert the ddb to the effective potential
      call ddb2effpot(Crystal,ddb, eff_pot,inp,comm)

!     If needed, print the effective potential into XML file
      if (inp%prt_effpot>=3.or.inp%prt_effpot==-1) then
        call effective_potential_print(eff_pot,-1)
      end if

    else if (filetype==2) then

!     Free the effective potential before 
      call effective_potential_free(eff_pot)
      
      call xml2effpot(eff_pot,filename,comm)

!     If needed, print the effective potential
      call effective_potential_print(eff_pot,inp%prt_effpot)

    end if

!   Common initilisations
!   Generate long rage interation for the effective potential for both type and generate suppercell
    call effective_potential_generateDipDip(eff_pot,inp%n_cell,inp%dipdip,inp%asr,comm)

  else
    write(message, '(a,a,a,a)' )&
&      ' The file ',trim(filename),' is not readable'
    MSG_BUG(message)
  end if

! Deallocation of array
  call crystal_free(Crystal)
  call ddb_free(ddb)

end subroutine effective_potential_file_read
!!***

!****f* m_effective_potential/effective_potential_file_readDisplacement
!!
!! NAME
!! effective_potential_file_readDisplacement
!!
!! FUNCTION
!! Read a strain file
!!
!! INPUTS
!! filename = name of the file
!! natom    = number of atoms
!! ntime    = number of time
!!
!! OUTPUT
!! disp = array with all the strain
!!
!! PARENTS
!!   multibinit
!!
!! CHILDREN
!!
!!
!! SOURCE
 
subroutine effective_potential_file_readDisplacement(filename,disp,nstep,natom)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_file_readDisplacement'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom,nstep
 character(len=fnlen),intent(in) :: filename
!array
 real(dp),intent(out) :: disp(nstep,3,natom)
!Local variables------------------------------
!scalars
 integer :: ia,istep,mu
 character(500) :: message
 character (len=500000) :: line,readline
 integer :: funit = 666
!array

! *************************************************************************

 if (open_file(filename,message,unit=funit,form="formatted",&
   status="old",action="read") /= 0) then
   MSG_ERROR(message)
 end if

 write(message, '(2a)' ) " Read displacements from ", trim(filename)

 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')
 do istep=1,nstep
   do ia=1,natom
     read(funit,'(a)',err=10,end=10) readline
     line=adjustl(readline)
     read(unit=line,fmt=*) (disp(istep,mu,ia),mu=1,3)
     write(111,'(3(es21.12))') disp(istep,:,ia)
   end do
 end do

10   continue
 close(funit)

end subroutine effective_potential_file_readDisplacement
!!***

!****f* m_effective_potential/effective_potential_file_readStrain
!!
!! NAME
!! effective_potential_file_readStrain
!!
!! FUNCTION
!! Read a strain file
!!
!! INPUTS
!! filename = name of the file
!! natom    = number of atoms
!! ntime    = number of time
!!
!! OUTPUT
!! disp = array with all the strain
!!
!! PARENTS
!!   multibinit
!!
!! CHILDREN
!!
!!
!! SOURCE
 
 subroutine effective_potential_file_readStrain(filename,disp,ntime,natom)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_file_readStrain'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer, intent(in) :: natom,ntime
  character(len=fnlen),intent(in) :: filename
!array
  real(dp),intent(out) :: disp(ntime,3,natom)
!Local variables------------------------------
!scalars
  integer :: funit = 1
  character(500) :: message
!array

! *************************************************************************

   if (open_file(filename,message,unit=funit,form="formatted",&
     status="old",action="read") /= 0) then
     MSG_ERROR(message)
   end if

end subroutine effective_potential_file_readStrain
!!***

end module m_effective_potential_file
