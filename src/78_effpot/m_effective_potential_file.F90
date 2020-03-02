!!****m* ABINIT/m_effective_potential_file
!! NAME
!! m_effective_potential_file
!!
!! FUNCTION
!! This  module contains all routine to read the effective potential from files
!! Can also read coefficients from XML
!! (XML or DDB)
!!
!! COPYRIGHT
!! Copyright (C) 2000-2020 ABINIT group (AM)
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
 use m_abicore
 use m_xmpi
 use m_harmonics_terms
 use m_anharmonics_terms
 use m_effective_potential
 use m_ifc
#if defined HAVE_NETCDF
 use netcdf
#endif

 use m_io_tools,   only : open_file
 use m_geometry,   only : xcart2xred, xred2xcart, metric
 use m_symfind,    only : symfind, symlatt
 use m_crystal,    only : crystal_t, crystal_init
 use m_dynmat,     only : dfpt_prtph
 use m_abihist,    only : abihist,abihist_init,abihist_free,abihist_copy,read_md_hist
 use m_ddb_internalstr, only : ddb_internalstr

 implicit none

 public :: effective_potential_file_getDimCoeff
 public :: effective_potential_file_getDimMD
 public :: effective_potential_file_getDimSystem
 public :: effective_potential_file_getDimStrainCoupling
 public :: effective_potential_file_getType
 public :: effective_potential_file_mapHistToRef
 public :: effective_potential_file_read
 public :: effective_potential_file_readDisplacement
 public :: effective_potential_file_readMDfile
 private :: coeffs_xml2effpot
 private :: system_getDimFromXML
 private :: system_xml2effpot
 private :: system_ddb2effpot

#ifndef HAVE_XML
 private :: rdfromline
 private :: rmtabfromline
 private :: rdfromline_value
 private :: elementfromline
#endif

#if defined HAVE_XML
 public :: effpot_xml_checkXML
 public :: effpot_xml_getDimCoeff
 public :: effpot_xml_readSystem
 public :: effpot_xml_getValue
 public :: effpot_xml_getAttribute
 public :: effpot_xml_getDimSystem

 interface
   subroutine effpot_xml_readSystem(filename,natom,&
&     ntypat,nrpt,nqpt,amu,atmfrc,cell,dynmat,elastic_constants,&
&     energy,epsilon_inf,ewald_atmfrc,&
&     phfrq,rprimd,qph1l,short_atmfrc,typat,xcart,zeff)&
&                          bind(C,name="effpot_xml_readSystem")
     use iso_c_binding, only : C_CHAR,C_DOUBLE,C_INT
     integer(C_INT) :: natom,ntypat,nrpt,nqpt
     integer(C_INT) :: typat(natom)
     integer(C_INT) :: cell(3,nrpt)
     real(C_DOUBLE) :: energy
     real(C_DOUBLE) :: dynmat(2,3,natom,3,natom,nqpt)
     real(C_DOUBLE) :: phfrq(3*natom,nqpt),qph1l(3,nqpt)
     real(C_DOUBLE) :: atmfrc(3,natom,3,natom,nrpt)
     real(C_DOUBLE) :: short_atmfrc(3,natom,3,natom,nrpt)
     real(C_DOUBLE) :: ewald_atmfrc(3,natom,3,natom,nrpt)
     real(C_DOUBLE) :: amu(ntypat),rprimd(3,3),epsilon_inf(3,3)
     real(C_DOUBLE) :: zeff(3,3,natom)
     real(C_DOUBLE) :: elastic_constants(6,6),xcart(3,natom)
     character(kind=C_CHAR) :: filename(*)
   end subroutine effpot_xml_readSystem
 end interface

 interface
   subroutine effpot_xml_readStrainCoupling(filename,natom,&
&     nrpt,voigt,elastic3rd,elastic_displacement,&
&     strain_coupling,phonon_strain_atmfrc,phonon_straincell)&
&                          bind(C,name="effpot_xml_readStrainCoupling")
     use iso_c_binding, only : C_CHAR,C_DOUBLE,C_INT
     integer(C_INT) :: natom
     integer(C_INT) :: nrpt,voigt
     integer(c_INT) :: phonon_straincell(3,nrpt)
     real(C_DOUBLE) :: elastic3rd(6,6),elastic_displacement(6,3,natom)
     real(C_DOUBLE) :: strain_coupling(3,natom)
     real(C_DOUBLE) :: phonon_strain_atmfrc(3,natom,3,natom,nrpt)
     character(kind=C_CHAR) :: filename(*)
   end subroutine effpot_xml_readStrainCoupling
 end interface

 interface
   subroutine effpot_xml_readCoeff(filename,ncoeff,ndisp,nterm,&
&                                 coefficient,atindx,cell,direction,power_disp,&
&                                 power_strain,strain,weight)&
&                          bind(C,name="effpot_xml_readCoeff")
     use iso_c_binding, only : C_CHAR,C_DOUBLE,C_INT
     character(kind=C_CHAR) :: filename(*)
     integer(C_INT) :: ncoeff,ndisp,nterm
     integer(C_INT) :: atindx(ncoeff,nterm,2,ndisp)
     integer(C_INT) :: cell(ncoeff,nterm,3,2,ndisp)
     integer(C_INT) :: direction(ncoeff,nterm,ndisp)
     integer(C_INT) :: strain(ncoeff,nterm,ndisp)
     integer(C_INT) :: power_disp(ncoeff,nterm,ndisp)
     integer(C_INT) :: power_strain(ncoeff,nterm,ndisp)
     real(C_DOUBLE) :: coefficient(ncoeff)
     real(C_DOUBLE) :: weight(ncoeff,nterm)
   end subroutine effpot_xml_readCoeff
 end interface

 interface
   subroutine effpot_xml_getDimSystem(filename,natom,ntypat,nqpt,nrpt1,nrpt2)&
&                          bind(C,name="effpot_xml_getDimSystem")
     use iso_c_binding, only : C_CHAR,C_INT
     integer(C_INT) :: natom,ntypat,nqpt,nrpt1,nrpt2
     character(kind=C_CHAR) :: filename(*)
   end subroutine effpot_xml_getDimSystem
 end interface

 interface
   subroutine effpot_xml_getDimStrainCoupling(filename,nrpt,voigt)&
&                          bind(C,name="effpot_xml_getDimStrainCoupling")
     use iso_c_binding, only : C_CHAR,C_INT
     integer(C_INT) :: voigt
     integer(C_INT) :: nrpt
     character(kind=C_CHAR) :: filename(*)
   end subroutine effpot_xml_getDimStrainCoupling
 end interface

 interface
   subroutine effpot_xml_getDimCoeff(filename,ncoeff,nterm_max,ndisp_max)&
&                          bind(C,name="effpot_xml_getDimCoeff")
     use iso_c_binding, only : C_CHAR,C_DOUBLE,C_INT,C_PTR
     character(kind=C_CHAR) :: filename(*)
!     character(kind=C_CHAR) :: name(*)
     type(C_PTR) :: name
     integer(C_INT) :: ncoeff,ndisp_max,nterm_max
   end subroutine effpot_xml_getDimCoeff
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

 interface
   subroutine effpot_xml_getNumberKey(filename,name_value,number) &
&                          bind(C,name="effpot_xml_getNumberKey")
     use iso_c_binding, only : C_CHAR,C_INT
     character(kind=C_CHAR) :: filename(*),name_value(*)
     integer(C_INT) :: number
   end subroutine effpot_xml_getNumberKey
 end interface

#endif

 integer,parameter :: XML_RECL = 50000
!!***

CONTAINS  !===========================================================================================


!****f* m_effective_potential_file/effective_potential_file_read
!!
!! NAME
!! effective_potential_file_read
!!
!! FUNCTION
!! tranfert file (XML or DDB) in effective potential type
!! Also transfert coefficient from xml file for ahnarmonic part
!!
!! INPUTS
!! filename = path of the file
!! hist<type(abihist)> = optional,The history of the MD (or snapshot of DFT)
!! inp<type(multibinit_dtset_type)> = optional,datatype with all the input variables (mantadory to
!!                                      read DDB file)
!! comm=MPI communicator
!!
!! OUTPUT
!! eff_pot<type(effective_potential_type)> = datatype with all the informations for effective potential
!!
!! PARENTS
!!      compute_anharmonics,multibinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine effective_potential_file_read(filename,eff_pot,inp,comm,hist)

  use m_effective_potential
  use m_multibinit_dataset
  use m_ddb, only : ddb_from_file
  use m_strain
  use m_crystal, only : crystal_t
  use m_dynmat, only : bigbx9

!Arguments ------------------------------------
!scalars
  integer,intent(in) :: comm
  character(len=fnlen),intent(in) :: filename
!array
  type(effective_potential_type), intent(inout)  :: eff_pot
  type(multibinit_dtset_type),optional,intent(in) :: inp
  type(ddb_type) :: ddb
  type(crystal_t) :: Crystal
  type(abihist),optional :: hist
!Local variables------------------------------
!scalars
  integer :: ii,filetype,natom,ntypat,nqpt,nrpt
  character(500) :: message
!array
  integer,allocatable :: atifc(:)

! *************************************************************************

  call effective_potential_file_getType(filename,filetype)

  if (filetype/=0) then

    if (.not.(present(inp))) then
      write(message, '(4a)' )&
&        ' effective_potential_file_read: you need to give input file to compute ',&
&        'the response fonction from DDB file ',ch10
      MSG_ERROR(message)
    end if

    if(filetype ==1) then
!     Read the DDB information, also perform some checks, and symmetrize partially the DDB
      write(message, '(3a)' )' Read the DDB information of the reference',&
 &      ' system and perform some checks',ch10
      call wrtout(std_out,message,'COLL')
      call wrtout(ab_out,message,'COLL')

      call effective_potential_file_getDimSystem(filename,natom,ntypat,nqpt,nrpt)!

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

!     Transfert the ddb to the effective potential
      call system_ddb2effpot(Crystal,ddb, eff_pot,inp,comm)

!     Generate long rage interation for the effective potential for both type and generate supercell
      call effective_potential_generateDipDip(eff_pot,inp%dipdip_range,inp%dipdip,inp%asr,comm)

!     If needed, print the effective potential into the output
      if (inp%prt_model>=3.or.inp%prt_model==-1) then
        call effective_potential_print(eff_pot,-1)
      end if
    end if
    if (filetype==2 .or.filetype==23) then

!     Free the effective potential before
      call effective_potential_free(eff_pot)

      call system_xml2effpot(eff_pot,filename,comm,strcpling=inp%strcpling)


!     Assign the energy of the reference from input
      if(abs(inp%energy_reference)>tol16)then
        write(message,'(11a)') ch10,&
&      ' --- !WARNING',ch10,&
&      '     Energy of the reference structure is specify in ',ch10,&
&      '     the input file. The energy is set with',ch10,&
&      '     this value.',ch10,&
&      ' ---',ch10
        call wrtout(std_out,message,'COLL')
        eff_pot%energy = inp%energy_reference
      end if


!     Generate long rage interation for the effective potential for both type and generate supercell
      call effective_potential_generateDipDip(eff_pot,inp%dipdip_range,inp%dipdip,inp%asr,comm)

!     If needed, print the effective potential
      call effective_potential_print(eff_pot,inp%prt_model)
    end if
    if (filetype==3 .or. filetype==23) then
!     Read the  coefficient of the fit for the anharmonic part
      write(message, '(4a)' )ch10,' Read the coefficients of the polynomial fit from XML',&
 &      ' and perform some checks',ch10
      call wrtout(std_out,message,'COLL')
      call wrtout(ab_out,message,'COLL')

      if(eff_pot%anharmonics_terms%ncoeff/=0)then
        write(message,'(9a)') ch10,&
&      ' --- !WARNING',ch10,&
&      '     There is already fitted polynome set in the model',ch10,&
&      '     The previous coefficients will be remove',ch10,&
&      ' ---',ch10
        call wrtout(std_out,message,'COLL')
      end if

      call coeffs_xml2effpot(eff_pot,filename,comm)

!     Assign the coeff number from input
      if(inp%ncoeff==0)then
        write(message,'(12a)') ch10,&
&      ' --- !WARNING',ch10,&
&      '     The values of the coefficients are set to 0',&
&      ' in the input file.',ch10,&
&      '     The values of the coefficients will be read in the XML',ch10,&
&      '     or might be fitted',ch10,&
&      ' ---',ch10
        call wrtout(std_out,message,'COLL')
        if(inp%fit_coeff <= 0 .and. &
&          all(abs(eff_pot%anharmonics_terms%coefficients(:)%coefficient) <tol16)) then

          write(message,'(12a)') ch10,&
&          ' --- !WARNING',ch10,&
&          '     The input for the fit process is set to 0 or -1',&
&          ' in the input file.',ch10,&
&          '     However, the values of the coefficients in the XMF files are zero,',ch10,&
&          '     So the coefficients can not be used',ch10,&
&          ' ---',ch10
          call wrtout(std_out,message,'COLL')

        end if
      else
        if (eff_pot%anharmonics_terms%ncoeff /= inp%ncoeff)then
          write(message, '(5a)' )&
&            ' The number of coefficients in the XML file is superior to the ',ch10,&
&            'number of coefficients in the input ',ch10,&
&            'Action: correct your input file or change the file'
          MSG_ERROR(message)
        end if
        do ii = 1,eff_pot%anharmonics_terms%ncoeff
          call polynomial_coeff_setCoefficient(inp%coefficients(ii),&
&                                              eff_pot%anharmonics_terms%coefficients(ii))
        end do
      end if

    else if(filetype==4) then
      if(present(hist))then
        write(message,'(5a)')ch10,&
&         '-Reading the file ',trim(filename),ch10,&
&         ' with NetCDF in order to fit the polynomial coefficients'
        call wrtout(std_out,message,'COLL')
        call wrtout(ab_out,message,'COLL')
        call effective_potential_file_readMDfile(filename,hist,option=inp%ts_option)
      else
       write(message, '(3a)' )&
&         'There is no hist argument ',ch10,&
&         'Action: add hist argument'
       MSG_ERROR(message)
     end if
   end if
 else
   write(message, '(5a)' )&
&      ' The file ',trim(filename),' is not readable with Multibinit',ch10,&
&      ' Action: Change the file.'
   MSG_BUG(message)
 end if

! Deallocation of array
  call crystal%free()
  call ddb%free()

end subroutine effective_potential_file_read
!!***

!!****f* m_effective_potential_file/effective_potential_file_getType
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
!!              2 XML file the system definition and harmonic part
!!              3 XML file with polynomial coefficients
!!             23 XML file with both system definition and polynomial coefficients
!!             40 NetCDF file with history of MD or snapshot
!!             41 ASCII file with history of MD or snapshot
!!
!! PARENTS
!!      m_effective_potential_file,multibinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine effective_potential_file_getType(filename,filetype)

!Arguments ------------------------------------
!scalars
 character(len=fnlen),intent(in) :: filename
 integer, intent(out) :: filetype
!arrays
!Local variables-------------------------------
!scalar
 integer :: natom,nstep
 integer :: ddbun = 666,ios=0
 character(len=500) :: message
 character (len=1000) :: line,readline
#if defined HAVE_NETCDF
 integer :: natom_id,time_id,xyz_id,six_id
 integer :: ncid,ncerr
 logical :: md_file
#endif

!arrays
! *************************************************************************

 filetype = 0

 if (open_file(filename,message,unit=ddbun,form="formatted",status="old",action="read") /= 0) then
   MSG_ERROR(message)
 end if

!Check if the file is a XML file or a DDB and in this case, store the DDB code.
 ios = 0
 do while ((ios==0))
   read(ddbun,'(a)',iostat=ios) readline
   call rmtabfromline(readline)
   line=adjustl(readline)
   if(line(3:13)=="xml version") then
     do while ((ios==0))
       read(ddbun,'(a)',iostat=ios) readline
       call rmtabfromline(readline)
       line=adjustl(readline)
       if(line(1:16)==char(60)//"Heff_definition".or.&
&         line(1:17)==char(60)//"Terms_definition")then
         filetype = 3
         ios = -1
       end if
       if(line(1:18)==char(60)//"System_definition") then
         filetype = 2
         do while ((ios==0))
           read(ddbun,'(a)',iostat=ios) readline
           call rmtabfromline(readline)
           line=adjustl(readline)
           if(line(1:16)==char(60)//"Heff_definition".or.&
&             line(1:17)==char(60)//"Terms_definition")then
             filetype = 23
             ios = -1
           end if
         end do
       end if
     end do
   else  if(line(6:24)=="DERIVATIVE DATABASE") then
     filetype = 1
     ios = -1
   end if
 end do
 close(ddbun)

 if(filetype/=0) return

!try to read netcdf HIST file
#if defined HAVE_NETCDF
 ncerr=nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=ncid)
 if(ncerr == NF90_NOERR) then
   md_file = .TRUE.
   ncerr = nf90_inq_dimid(ncid,"natom",natom_id)
   if(ncerr /= NF90_NOERR)  md_file = .FALSE.
   ncerr = nf90_inq_dimid(ncid,"xyz",xyz_id)
   if(ncerr /= NF90_NOERR)  md_file = .FALSE.
   ncerr = nf90_inq_dimid(ncid,"time",time_id)
   if(ncerr /= NF90_NOERR)  md_file = .FALSE.
   ncerr = nf90_inq_dimid(ncid,"six",six_id)
   if(ncerr /= NF90_NOERR)  md_file = .FALSE.
   if (md_file) then
     filetype = 40
     return
   end if
 end if
 ncerr = nf90_close(ncid)
#endif

 if(filetype/=0) return

!Try to get the dim of MD ASCII file
 call effective_potential_file_getDimMD(filename,natom,nstep)
 if(natom /= 0 .and. nstep/=0) filetype = 41

end subroutine effective_potential_file_getType
!!***

!!****f* m_effective_potential_file/effective_potential_file_getDimSystem
!!
!! NAME
!! effective_potential_file_getDimSystem
!!
!! FUNCTION
!! This routine test the xml or ddb file
!! Return the number of atoms/ntypat in the unit cell from ddb and xml
!! Return natom/ntypat/nqpt and nrpt if the file is XML file
!! In case of DDB file, you have to run bigbx9 to get nrpt
!!
!! INPUTS
!! filename = names of the files
!!
!! OUTPUT
!! natom = number of atoms
!! ntypat= number of type of atoms
!! nqpt  = number of q points
!! nrpt  = number of rpt points
!!
!! PARENTS
!!      m_effective_potential_file,multibinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine effective_potential_file_getDimSystem(filename,natom,ntypat,nqpt,nrpt)

 use m_ddb
 use m_ddb_hdr

!Arguments ------------------------------------
!scalars
 character(len=fnlen),intent(in) :: filename
 integer,intent(out) :: natom,ntypat,nqpt,nrpt
!arrays

!Local variables-------------------------------
 !scalar
 integer :: filetype
! integer :: dimekb,lmnmax,mband,mtyp,msym,nblok,nkpt,usepaw
 integer :: ddbun = 666
 character(len=500) :: message
 type(ddb_hdr_type) :: ddb_hdr
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

   call ddb_hdr_open_read(ddb_hdr,filename,ddbun,DDB_VERSION,&
&                         dimonly=1)
   natom = ddb_hdr%natom
   ntypat = ddb_hdr%ntypat

   call ddb_hdr_free(ddb_hdr)

!  Must read some value to initialze  array (nprt for ifc)
!   call bigbx9(inp%brav,dummy_cell,0,1,inp%ngqpt,inp%nqshft,nrpt,ddb%rprim,dummy_rpt)

 else if (filetype==2 .or. filetype==23) then
   write(message, '(5a)' )ch10,' The file ',trim(filename),&
&                ' is XML file (extraction of all informations)'
   call wrtout(std_out,message,'COLL')

   call system_getDimFromXML(filename,natom,ntypat,nqpt,nrpt)

 else
   write(message, '(a,a,a,a)' )&
&   ' The file ',trim(filename),' is not compatible with multibinit',ch10
   MSG_ERROR(message)
 end if

! TODO hexu: temporarily disabled. Discuss with alex how to do this properly.
! Do some checks
! if (natom < 1) then
!   write(message, '(a,a,a,a,a)' )&
!&   ' Unable to read the number of atom from ',trim(filename),ch10,&
!&   'This file  is not compatible with multibinit',ch10
!   MSG_ERROR(message)
! end if
!
! if (filetype==2 .or. filetype==23) then
!
!   if (natom < 1) then
!     write(message, '(a,a,a)' )&
!&     ' Unable to read the number of atom from ',trim(filename),ch10
!     MSG_ERROR(message)
!   end if
!
!   if (nrpt < 1) then
!     write(message, '(a,a,a)' )&
!&     ' Unable to read the number of rpt points ',trim(filename),ch10
!     MSG_ERROR(message)
!   end if
!
!   if (ntypat < 1) then
!     write(message, '(a,a,a)' )&
!&     ' Unable to read the number of type of atoms ',trim(filename),ch10
!     MSG_ERROR(message)
!   end if
!
! end if

end subroutine effective_potential_file_getDimSystem
!!***

!!****f* m_effective_potential_file/effective_potential_file_getDimCoeff
!!
!! NAME
!! effective_potential_file_getDimCoeff
!!
!! FUNCTION
!! This routine test the xml with polynomial coefficients
!! Return the number of coefficients and the maximum number of displacement/strain
!!
!! INPUTS
!! filename = names of the files
!!
!! OUTPUT
!! ncoeff = number of coefficient for the polynome
!! nterm(ncoeff) = number terms per coefficient
!! ndisp(nterm,ncoeff) = number displacement per term
!!
!! PARENTS
!!      m_effective_potential_file
!!
!! CHILDREN
!!
!! SOURCE

subroutine effective_potential_file_getDimCoeff(filename,ncoeff,ndisp_max,nterm_max)

!Arguments ------------------------------------
!scalars
 character(len=fnlen),intent(in) :: filename
 integer,intent(out) :: ncoeff,ndisp_max,nterm_max
!Local variables-------------------------------
 !scalar
 integer ::  filetype
#ifndef HAVE_XML
 integer ::  count,count2
 integer :: funit = 1,ios=0
 logical :: found,found2
#endif
!arrays
#ifndef HAVE_XML
 character (len=XML_RECL) :: line,readline
#endif
 character(len=500) :: message

! *************************************************************************

 call effective_potential_file_getType(filename,filetype)

 if (filetype==3 .or. filetype==23) then
   write(message, '(2a)' )' Extraction of the number of coefficient in the XML ',&
&                         trim(filename)
   call wrtout(std_out,message,'COLL')

   ncoeff = 0
   nterm_max = 0
   ndisp_max = 0

#if defined HAVE_XML
!  Read with libxml the number of coefficient
   call effpot_xml_getDimCoeff(char_f2c(trim(filename)),ncoeff,nterm_max,ndisp_max)
#else
!  Read by hand
!  Start a reading loop
   found=.false.
   ncoeff = 0

   if (open_file(filename,message,unit=funit,form="formatted",status="old",&
&                action="read") /= 0) then
     MSG_ERROR(message)
   end if

!  First parse to know the number of coefficients
   ios = 0
   do while (ios == 0)
     read(funit,'(a)',iostat=ios) readline
     if(ios == 0)then
       call rmtabfromline(readline)
       line=adjustl(readline)
!      Need test with char(9) because the old version of XML file
!      from old script includes tarbulation at the begining of each line
       if (line(1:12)==char(60)//'coefficient') then
         ncoeff=ncoeff+1
         count = 0
         found = .false.
         do while(.not.found)
           read(funit,'(a)',iostat=ios) readline
           call rmtabfromline(readline)
           line=adjustl(readline)
           if (line(1:5)==char(60)//'term') then
             count = count +1
             found2 = .false.
             count2 = 0
             do while(.not.found2)
               read(funit,'(a)',iostat=ios) readline
               call rmtabfromline(readline)
               line=adjustl(readline)
               if (line(1:13)==char(60)//'displacement') then
                 count2 = count2 + 1
               else if (line(1:7)==char(60)//'strain') then
                 count2 = count2 + 1
               else if (line(1:6)==char(60)//'/term') then
                 if (count2 > ndisp_max) ndisp_max = count2
                 found2 = .true.
               else
                 cycle
               end if
             end do
           else  if (line(1:13)==char(60)//'/coefficient') then
             if (count > nterm_max) nterm_max = count
             found = .true.
           else
             cycle
           end if
         end do
         cycle
       end if
     end if
   end do

   close(funit)
#endif

 else
!  Maybe one day add an other type of file...
   write(message, '(a,a,a,a)' )&
&   ' The file ',trim(filename),' is not compatible with multibinit',ch10
   MSG_ERROR(message)
 end if

! Do some checks
 if (ncoeff < 1) then
   write(message, '(5a)' )&
&   ' Unable to read the number of coeff from ',trim(filename),ch10,&
&   ' This file is not compatible with multibinit',ch10
   MSG_ERROR(message)
 end if

end subroutine effective_potential_file_getDimCoeff
!!***


!!****f* m_effective_potential_file/effective_potential_file_getDimStrainCoupling
!!
!! NAME
!! effective_potential_file_getDimStrainCoupling
!!
!! FUNCTION
!! Return the number of nrpt for specific strain coupling from xml system file
!!
!! INPUTS
!! filename = names of the files
!! voigt    = voigt notation of the strain
!!
!! OUTPUT
!! nrpt  = number of rpt points
!!
!! PARENTS
!!      m_effective_potential_file
!!
!! CHILDREN
!!
!! SOURCE

subroutine effective_potential_file_getDimStrainCoupling(filename,nrpt,voigt)

!Arguments ------------------------------------
!scalars
 character(len=fnlen),intent(in) :: filename
 integer,intent(in) :: voigt
 integer,intent(out) :: nrpt
!Local variables-------------------------------
 !scalar
#ifndef HAVE_XML
 integer :: irpt,ivoigt
 integer :: funit = 1,ios=0
 logical :: found
#endif
!arrays
#ifndef HAVE_XML
 character (len=XML_RECL) :: line,readline,strg,strg1
 character(len=500) :: message
#endif

! *************************************************************************

   nrpt = 0

#if defined HAVE_XML
!  Read with libxml the number of coefficient
   call effpot_xml_getDimStrainCoupling(char_f2c(trim(filename)),nrpt,voigt)
#else
!  Read by hand
!  Start a reading loop
   found=.false.

   if (open_file(filename,message,unit=funit,form="formatted",status="old",&
&                action="read") /= 0) then
     MSG_ERROR(message)
   end if

!  First parse to know the number of atoms
   do while (ios == 0.and.(.not.found))
     read(funit,'(a)',iostat=ios) readline
     if(ios == 0)then
       call rmtabfromline(readline)
       line=adjustl(readline)
       if ((line(1:16)=='<strain_coupling')) then
         read(funit,'(a)',iostat=ios) readline
         call rdfromline("voigt",line,strg)
         strg1=trim(strg)
         read(strg1,*) ivoigt
         if (ivoigt == voigt)then
           irpt = 0
           do while (.not.found)
             read(funit,'(a)',iostat=ios) readline
             call rmtabfromline(readline)
             line=adjustl(readline)
             if ((line(1:26)=='<correction_force_constant')) then
               irpt = irpt + 1
               cycle
             end if
             if ((line(1:17)=='</strain_coupling')) then
               found = .TRUE.
               nrpt = irpt
               cycle
             end if
           end do
         else
           cycle
         end if
       end if
     end if
   end do

   close(funit)
#endif

end subroutine effective_potential_file_getDimStrainCoupling
!!***

!!****f* m_effective_potential_file/effective_potential_file_getDimMD
!!
!! NAME
!! effective_potential_file_getDimMD
!!
!! FUNCTION
!! Read MD FILE (HIST or ASCII) and return the dimensions
!! (natom and nstep)
!!
!! INPUTS
!! filename = path of the file
!!
!! OUTPUT
!! natom = number of atoms
!! nstep = number of MD steps
!!
!! PARENTS
!!      m_effective_potential_file
!!
!! CHILDREN
!!
!! SOURCE

subroutine effective_potential_file_getDimMD(filename,natom,nstep)

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: natom,nstep
!arrays
 character(len=fnlen),intent(in) :: filename
!Local variables-------------------------------
!scalar
 integer :: ia,natm_old,natm_new
 integer :: nenergy,nrprimd
 integer :: ios=0,ios2=0,ios3=0
 integer :: unit_md=24
 logical :: compatible,netcdf
#if defined HAVE_NETCDF
 integer :: natom_id,time_id,xyz_id,six_id
 integer :: ncid,ncerr
 character(len=5) :: char_tmp
#endif
!arrays
 character (len=10000) :: readline,line
 character(len=500) :: msg

! *************************************************************************

 natom = 0
 nstep = 0
!try to read netcdf
 netcdf = .false.
#if defined HAVE_NETCDF
 ncerr=nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=ncid)
 if(ncerr == NF90_NOERR) then
   netcdf = .TRUE.
   ncerr = nf90_inq_dimid(ncid,"natom",natom_id)
   if(ncerr /= NF90_NOERR)  netcdf = .FALSE.
   ncerr = nf90_inq_dimid(ncid,"xyz",xyz_id)
   if(ncerr /= NF90_NOERR)  netcdf = .FALSE.
   ncerr = nf90_inq_dimid(ncid,"time",time_id)
   if(ncerr /= NF90_NOERR)  netcdf = .FALSE.
   ncerr = nf90_inq_dimid(ncid,"six",six_id)
   if(ncerr /= NF90_NOERR)  netcdf = .FALSE.
   if(netcdf)then
     ncerr = nf90_inquire_dimension(ncid,natom_id,char_tmp,natom)
     NCF_CHECK_MSG(ncerr," inquire dimension ID for natom")
     ncerr = nf90_inquire_dimension(ncid,time_id,char_tmp,nstep)
     NCF_CHECK_MSG(ncerr," inquire dimension ID for time")
   end if
 end if
#endif

 if(.not.netcdf) then
!  try to read ASCII file...
   if (open_file(filename,msg,unit=unit_md,form="formatted",&
&       status="old",action="read") /= 0) then
     MSG_ERROR(msg)
   end if

!  Start a reading loop in fortran to get the dimension of the file
   rewind(unit=unit_md)
   ios = 0
   nstep   = 0
   nrprimd = 0
   natm_old= 0
   natm_new= 0
   nenergy = 0
   compatible = .TRUE.

   do while ((ios==0))
!    special treatment of the first step
     if(nstep==0)then
       ios2 = 0
       do while ((ios2==0))
         read(unit_md,'(a)',iostat=ios) readline
         line=adjustl(readline)
         call elementfromline(line,ia)
         if (ia==1)then
           nstep = nstep + 1
           ios2 = 1
         end if
       end do
     end if
     read(unit_md,'(a)',iostat=ios) readline
     if(ios == 0)then
       line=adjustl(readline)
       call elementfromline(line,ia)
       if (ia==1)then
         nenergy = nenergy + 1
         nrprimd = 0
       else if(ia==3)then
         nrprimd = nrprimd + 1
       end if
       if(nrprimd == 3)then
         ios3 = 0
         natm_new = 0
         do while ((ios3==0))
           read(unit_md,'(a)',iostat=ios3) readline
           if(ios3==0)then
             line=adjustl(readline)
             call elementfromline(line,ia)
             if(ia==1)then
               if(nstep==1) then
                 natm_old = natm_new
               else
                 if(natm_old /= natm_new) compatible = .FALSE.
               end if
               ios3 = 1
               ios2 = 1
               nstep = nstep + 1
             end if
             if(ia==6)then
               natm_new = natm_new + 1
             end if
           end if!end if ios3
         end do
       end if ! end if nrprimd
     end if! end if os1
   end do

   natom = natm_new - 1
   if(nstep /= nenergy) compatible = .FALSE.
   if(natom <= 0) compatible = .FALSE.
   if(nstep <= 0) compatible = .FALSE.

   if(.not.compatible)then
     natom = 0
     nstep = 0
   end if
 end if! end if not netcdf

end subroutine effective_potential_file_getDimMD
!!***

!!****f* m_effective_potential_file/system_getDimFromXML
!! NAME
!! system_getDimFromXML
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
!!      m_effective_potential_file
!!
!! CHILDREN
!!
!! SOURCE

subroutine system_getDimFromXML(filename,natom,ntypat,nph1l,nrpt)

 !Arguments ------------------------------------
 !scalars
  character(len=fnlen),intent(in) :: filename
  integer, intent(out) :: natom,ntypat,nph1l,nrpt
 !arrays
 !Local variables-------------------------------
 !scalar
  integer :: nrpt1,nrpt2
  real :: itypat
  character(len=500) :: message
#ifndef HAVE_XML
  integer :: funit = 1,ios = 0
  integer :: iatom
  logical :: found
  character (len=XML_RECL) :: line,readline
  character (len=XML_RECL) :: strg,strg1
#endif
  !arrays
#ifndef HAVE_XML
  integer,allocatable :: typat(:)
#endif

 ! *************************************************************************

!Open the atomicdata XML file for reading
 write(message,'(5a)') ' system_getDimFromXML :',&
&    '-Opening the file ',trim(filename),' to read dimensions',&
&    ' (before initialisation)'

 call wrtout(std_out,message,'COLL')

 natom = 0
 ntypat= 0
 nph1l = 0
 nrpt  = 0
 nrpt1 = 0
 nrpt2 = 0
 itypat= 0

!Open the atomicdata XML file for reading

#if defined HAVE_XML
!Read with libxml
 call effpot_xml_getDimSystem(char_f2c(trim(filename)),natom,ntypat,nph1l,nrpt1,nrpt2)
#else
!Read by hand

!Start a reading loop
 found=.false.

 if (open_file(filename,message,unit=funit,form="formatted",status="old",&
&              action="read") /= 0) then
   MSG_ERROR(message)
 end if

!First parse to know the number of atoms
 do while ((ios==0).and.(.not.found))
   read(funit,'(a)',iostat=ios) readline
   if(ios ==0)then
     call rmtabfromline(readline)
     line=adjustl(readline)

!Need test with char(9) because the old version of XML file
!from old script includes tarbulation at the begining of each line
     if (line(1:5)==char(60)//'atom') then
       natom=natom+1
       cycle
     end if

     if (line(1:21)==char(60)//'local_force_constant') then
       nrpt1 = nrpt1+1
       cycle
     end if

     if (line(1:21)==char(60)//'total_force_constant') then
       nrpt2 = nrpt2+1
       cycle
     end if

     if (line(1:7)==char(60)//'qpoint') then
       nph1l =  nph1l+1
       cycle
     end if
   end if
 end do

!second parse to get the number of typat
 ABI_ALLOCATE(typat,(natom))
 typat = 0
 iatom = 0

 rewind(funit)
!Start a reading loop
 ios   = 0
 found = .false.

 do while ((ios==0).and.(.not.found))
   read(funit,'(a)',iostat=ios) readline
   if(ios == 0)then
     call rmtabfromline(readline)
     line=adjustl(readline)

     if (line(1:5)==char(60)//'atom') then
       iatom = iatom + 1
       call rdfromline("mass",line,strg)
       strg1=trim(strg)
       read(unit=strg1,fmt=*) itypat
       if (.not.any(typat==int(itypat))) then
         ntypat= ntypat+1
       end if
       typat(iatom) = int(itypat)
     end if

     if (line(1:6)==char(60)//'local') then
       found=.true.
     end if
   end if
 end do

 close(funit)
 ABI_DEALLOCATE(typat)

#endif

!Check the RPT
 if (nrpt2/=nrpt1) then
   if(nrpt1> 0 .and. nrpt2== 0) then
     continue;
   else if (nrpt1==0.and.nrpt2>=0) then
     write(message, '(5a)' )ch10,&
&   ' WARNING: the number of local IFC is set to 0  ',ch10,&
&   '          Dipdip must be set to zero',ch10
     call wrtout(std_out,message,'COLL')
   else if (nrpt2 > nrpt1) then
     write(message, '(2a,I0,3a,I0,5a)' )ch10,&
&   ' WARNING: the number of total IFC  (',nrpt2,') is not equal to the  ',ch10,&
&   '          the number of short range IFC (',nrpt1,') in ',filename,ch10,&
&   '          the missing ifc will be set to zero',ch10
     call wrtout(std_out,message,'COLL')
   else if(nrpt1>nrpt2)then
     write(message, '(2a,I0,3a,I0,5a)' )ch10,&
&   ' The number of total IFC  (',nrpt2,') is inferior to  ',ch10,&
&   ' the number of short range IFC (',nrpt1,') in ',filename,ch10,&
&   ' This is not possible',ch10
     MSG_BUG(message)
   end if
 end if

!nrpt is the max between local and total:
 nrpt = max(nrpt1,nrpt2)


end subroutine system_getDimFromXML
!!***

!!****f* m_effective_potential_file/system_xml2effpot
!! NAME
!! system_xml2effpot
!!
!! FUNCTION
!! Open xml file of effective potentiel, then reads the variables
!! and store them in effective potentential type
!!
!! INPUTS
!! eff_pot<type(effective_potential_type)> = datatype with all the informations for effective potential
!! comm=MPI communicator
!! character(len=*) filnam: name of input or output file
!! strcpling = optional,logical to disable the strcpling
!!
!! OUTPUT
!! eff_pot<type(effective_potential_type)> = datatype with all the informations for effective potential
!!
!! PARENTS
!!      m_effective_potential_file
!!
!! CHILDREN
!!
!! SOURCE

 subroutine system_xml2effpot(eff_pot,filename,comm,strcpling)

 use m_atomdata
 use m_effective_potential, only : effective_potential_type
 use m_multibinit_dataset, only : multibinit_dtset_type
 use m_ab7_symmetry

 !Arguments ------------------------------------
 !scalars
 character(len=*),intent(in) :: filename
 integer, intent(in) :: comm
 integer, optional,intent(in) :: strcpling
 !arrays
 type(effective_potential_type), intent(inout) :: eff_pot

 !Local variables-------------------------------
 !scalar
 integer :: ierr,ii,itypat,my_rank,msym,natom,ncoeff,nrpt,nrpt_scoupling
 integer :: ntypat,nph1l,nptsym,npsp,nproc,nsym,space_group,timrev,use_inversion,voigt
 real(dp):: energy,tolsym,ucvol
 character(len=500) :: message
 integer,parameter :: master=0
 logical :: has_anharmonics = .FALSE.
 logical :: iam_master
#ifndef HAVE_XML
 integer :: funit = 1,ios=0
 integer :: iatom,iamu,iph1l,irpt,irpt1,irpt2,irpt3,jj,mu,nu
 real(dp):: amu
 logical :: found,found2,short_range,total_range
 character (len=XML_RECL) :: line,readline
 character (len=XML_RECL) :: strg,strg1
 logical :: has_straincoupling = .FALSE.
#endif
 !arrays
 integer :: bravais(11)
 integer,allocatable :: typat(:)
 integer,allocatable  :: symrel(:,:,:),symafm(:),ptsymrel(:,:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3)
 real(dp) :: elastic_constants(6,6),elastic3rd(6,6,6),epsilon_inf(3,3)
 real(dp),allocatable :: all_amu(:),cell_local(:,:),cell_total(:,:)
 real(dp),allocatable :: elastic_displacement(:,:,:,:),dynmat(:,:,:,:,:,:)
 real(dp),allocatable :: local_atmfrc(:,:,:,:,:),total_atmfrc(:,:,:,:,:)
 real(dp),allocatable :: spinat(:,:),strain_coupling(:,:,:),phfrq(:,:),qph1l(:,:),tnons(:,:)
 real(dp),allocatable :: xcart(:,:),xred(:,:),zeff(:,:,:),znucl(:),zion(:)
 character(len=132),allocatable :: title(:)
 type(ifc_type) :: ifcs
 type(ifc_type),dimension(:),allocatable :: phonon_strain
 type(crystal_t)  :: crystal
 type(atomdata_t) :: atom
#ifdef HAVE_XML
 real(dp),allocatable :: phonon_strain_atmfrc(:,:,:,:,:)
 integer,allocatable  :: phonon_straincell(:,:)
#endif
#ifndef HAVE_XML
 real(dp),allocatable :: work2(:,:)
#endif

! *************************************************************************


 !Open the atomicdata XML file for reading
 write(message,'(a,a)')'-Opening the file ',filename

 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 iam_master = (my_rank == master)

!Get Dimention of system and allocation/initialisation of array
 call effective_potential_file_getDimSystem(filename,natom,ntypat,nph1l,nrpt)
 gmet= zero; gprimd = zero; rmet = zero; rprimd = zero
 elastic_constants = zero; epsilon_inf = zero; ncoeff = 0
 ABI_ALLOCATE(all_amu,(ntypat))
 ABI_ALLOCATE(cell_local,(3,nrpt))
 ABI_ALLOCATE(cell_total,(3,nrpt))
 ABI_ALLOCATE(elastic_displacement,(6,6,3,natom))
 ABI_ALLOCATE(ifcs%atmfrc,(3,natom,3,natom,nrpt))
 ABI_ALLOCATE(ifcs%cell,(3,nrpt))
 ABI_ALLOCATE(ifcs%short_atmfrc,(3,natom,3,natom,nrpt))
 ABI_ALLOCATE(ifcs%ewald_atmfrc,(3,natom,3,natom,nrpt))
 ABI_ALLOCATE(strain_coupling,(6,3,natom))
 ABI_ALLOCATE(total_atmfrc,(3,natom,3,natom,nrpt))
 ABI_ALLOCATE(local_atmfrc,(3,natom,3,natom,nrpt))
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
 nrpt_scoupling = 0
 do ii = 1,6
!  Get The size of the strainPhonon-coupling
   call effective_potential_file_getDimStrainCoupling(filename,nrpt_scoupling,ii-1)
   ABI_ALLOCATE(phonon_strain(ii)%atmfrc,(3,natom,3,natom,nrpt_scoupling))
   ABI_ALLOCATE(phonon_strain(ii)%cell,(3,nrpt_scoupling))
   phonon_strain(ii)%nrpt   = nrpt_scoupling
   phonon_strain(ii)%atmfrc = zero
   phonon_strain(ii)%cell   = 0
 end do

 all_amu(:) = zero
 dynmat(:,:,:,:,:,:)  = zero
 cell_local(:,:) = 99D99
 cell_total(:,:) = 99D99
 elastic3rd(:,:,:) = zero
 elastic_displacement(:,:,:,:) = zero
 ifcs%nrpt = nrpt
 ifcs%atmfrc(:,:,:,:,:)  = zero
 ifcs%cell(:,:)  = 0
 ifcs%ewald_atmfrc(:,:,:,:,:) = zero
 ifcs%short_atmfrc(:,:,:,:,:) = zero
 strain_coupling(:,:,:) = zero
 phfrq = zero
 qph1l = 0
 xcart = zero
 zeff  = zero
 znucl = zero

 if(iam_master)then
!Open the atomicdata XML file for reading
#if defined HAVE_XML

   write(message,'(a,a,a,a)')'-Reading the file ',trim(filename),&
&   ' with LibXML library'

   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

!  Read with libxml library
   call effpot_xml_readSystem(char_f2c(trim(filename)),natom,ntypat,nrpt,nph1l,all_amu,&
&                       ifcs%atmfrc,ifcs%cell,dynmat,elastic_constants,energy,&
&                       epsilon_inf,ifcs%ewald_atmfrc,phfrq,rprimd,qph1l,&
&                       ifcs%short_atmfrc,typat,xcart,zeff)

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

!  Get the Phonon Strain coupling
   do voigt = 1,6
     nrpt_scoupling = phonon_strain(voigt)%nrpt
     ABI_ALLOCATE(phonon_straincell,(3,nrpt_scoupling))
     ABI_ALLOCATE(phonon_strain_atmfrc,(3,natom,3,natom,nrpt_scoupling))

!      Get The value
       call effpot_xml_readStrainCoupling(char_f2c(trim(filename)),natom,nrpt_scoupling,(voigt-1),&
&                                         elastic3rd(voigt,:,:),elastic_displacement(voigt,:,:,:),&
&                                         strain_coupling(voigt,:,:),&
&                                         phonon_strain_atmfrc,phonon_straincell)

!      Check if the 3rd order strain_coupling is present
       if(any(elastic3rd>tol10).or.any(elastic_displacement>tol10)) has_anharmonics = .TRUE.
       phonon_strain(voigt)%atmfrc(:,:,:,:,:) = phonon_strain_atmfrc(:,:,:,:,:)
       phonon_strain(voigt)%cell(:,:)   = phonon_straincell(:,:)
       if(any(phonon_strain(voigt)%atmfrc > tol10)) has_anharmonics = .TRUE.

       ABI_DEALLOCATE(phonon_straincell)
       ABI_DEALLOCATE(phonon_strain_atmfrc)
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
   found=.false.

   iatom  = 1
   iamu   = 1
   itypat = 1
   irpt   = 1
   irpt1  = 0
   irpt2  = 0
   iph1l  = 1
   amu    = zero
   short_range  = .false.
   total_range  = .false.

   do while ((ios==0).and.(.not.found))
     read(funit,'(a)',iostat=ios) readline
     if(ios == 0)then
       call rmtabfromline(readline)
       line=adjustl(readline)
       if (.not.has_straincoupling) then

         if ((line(1:7)=='<energy')) then
           call rdfromline_value('energy',line,strg)
           if (strg/="") then
             strg1=trim(strg)
             read(strg1,*) energy
           else
             read(funit,'(a)',iostat=ios) readline
             call rmtabfromline(readline)
             line=adjustl(readline)
             call rdfromline_value('energy',line,strg)
             if (strg/="") then
               strg1=trim(strg)
             else
               strg1=trim(line)
             end if
             read(strg1,*) energy
           end if
           cycle
         end if

         if ((line(1:10)=='<unit_cell')) then
           call rdfromline_value('unit_cell',line,strg)
           if (strg/="") then
             strg1=trim(strg)
             read(strg1,*) (rprimd(1,mu),mu=1,3)
             read(funit,*) (rprimd(2,mu),mu=1,3)
           else
             do nu=1,2
               read(funit,*) (rprimd(nu,mu),mu=1,3)
             end do
           end if
           read(funit,'(a)',iostat=ios) readline
           call rmtabfromline(readline)
           line=adjustl(readline)
           call rdfromline_value('unit_cell',line,strg)
           if (strg/="") then
             strg1=trim(strg)
           else
             strg1=trim(line)
           end if
           read(strg1,*) (rprimd(3,mu),mu=1,3)
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
           read(funit,'(a)',iostat=ios) readline
           call rmtabfromline(readline)
           line=adjustl(readline)
           call rdfromline_value('epsilon_inf',line,strg)
           if (strg/="") then
             strg1=trim(strg)
           else
             strg1=trim(line)
           end if
           read(strg1,*) (epsilon_inf(mu,3),mu=1,3)
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
           read(funit,'(a)',iostat=ios) readline
           call rmtabfromline(readline)
           line=adjustl(readline)
           call rdfromline_value('elastic',line,strg)
           if (strg/="") then
             strg1=trim(strg)
           else
             strg1=trim(line)
           end if
           read(strg1,*) (elastic_constants(mu,6),mu=1,6)
           cycle
         end if

         if ((line(1:5)=='<atom')) then
           call rdfromline("mass",line,strg)
           strg1=trim(strg)
           read(strg1,*) amu
           if (.not.any(abs(all_amu-amu)<tol16)) then
             all_amu(iamu) = amu
             typat(iatom) = int(amu)
             !convert atomic mass unit to znucl
             do ii=1,103
               call atomdata_from_znucl(atom,real(ii,dp))
               if (abs((real(atom%amu,sp)-real(amu,sp))&
&                 /real(amu,sp)*100)<0.1) then
                 znucl(iamu) = atom%znucl
                 exit
               end if
             end do
             iamu = iamu +1
           end if
           do itypat=1,ntypat
             if(abs(amu-all_amu(itypat))<tol16) then
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
             read(funit,'(a)',iostat=ios) readline
             call rmtabfromline(readline)
             line=adjustl(readline)
             call rdfromline_value('position',line,strg)
             if (strg/="") then
               strg1=trim(strg)
             else
               strg1=trim(line)
             end if
             read(strg1,*)(xcart(mu,iatom),mu=1,3)
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
           read(funit,'(a)',iostat=ios) readline
           line=adjustl(readline)
           call rdfromline_value('borncharge',line,strg)
           if (strg/="") then
             strg1=trim(strg)
           else
             strg1=trim(line)
           end if
             read(strg1,*) (zeff(mu,3,iatom),mu=1,3)
           cycle
         end if

         if ((line(1:7)==char(60)//char(47)//'atom'//char(62))) then
           iatom=iatom+1
           cycle
         end if

         if ((line(1:12)=='<local_force')) then
           found2 = .False.
           irpt1 = irpt1 + 1
           do while (.not.found2)
             read(funit,'(a)',iostat=ios) readline
             call rmtabfromline(readline)
             line=adjustl(readline)
             if ((line(1:5)=='<data')) then
               call rdfromline_value('data',line,strg)
               if (strg/="") then
                 ABI_ALLOCATE(work2,(3*natom,3*natom))
                 strg1=trim(strg)
                 read(strg1,*) (work2(1,nu),nu=1,3*natom)
                 do mu=2,3*natom-1
                   read(funit,*)(work2(mu,nu),nu=1,3*natom)
                 end do
                 read(funit,'(a)',iostat=ios) readline
                 call rmtabfromline(readline)
                 line=adjustl(readline)
                 call rdfromline_value('data',line,strg)
                 if (strg/="") then
                   strg1=trim(strg)
                 else
                   strg1=trim(line)
                 end if
                 read(strg1,*) (work2(3*natom,nu),nu=1,3*natom)
                 local_atmfrc(:,:,:,:,irpt1) = reshape(work2,(/3,natom,3,natom/))
                 ABI_DEALLOCATE(work2)
               else
                 ABI_ALLOCATE(work2,(3*natom,3*natom))
                 do mu=1,3*natom
                   read(funit,*)(work2(mu,nu),nu=1,3*natom)
                 end do
                 local_atmfrc(:,:,:,:,irpt1) =  reshape(work2,(/3,natom,3,natom/))
                 ABI_DEALLOCATE(work2)
               end if
             end if
             if ((line(1:5)=='<cell')) then
               call rdfromline_value('cell',line,strg)
               if (strg/="") then
                 strg1=trim(strg)
                 read(strg1,*)(cell_local(mu,irpt1),mu=1,3)
               else
                 read(funit,*)(cell_local(mu,irpt1),mu=1,3)
               end if
               found2 = .TRUE.
               cycle
             end if
           end do
         end if

         if ((line(1:12)=='<total_force')) then
           irpt2 = irpt2 + 1
           found2 = .False.
           do while (.not.found2)
             read(funit,'(a)',iostat=ios) readline
             call rmtabfromline(readline)
             line=adjustl(readline)
             if ((line(1:5)=='<data')) then
               call rdfromline_value('data',line,strg)
               if (strg/="") then
                 ABI_ALLOCATE(work2,(3*natom,3*natom))
                 strg1=trim(strg)
                 read(strg1,*) (work2(1,nu),nu=1,3*natom)
                 do mu=2,3*natom-1
                   read(funit,*)(work2(mu,nu),nu=1,3*natom)
                 end do
                 read(funit,'(a)',iostat=ios) readline
                 call rmtabfromline(readline)
                 line=adjustl(readline)
                 call rdfromline_value('data',line,strg)
                 if (strg/="") then
                   strg1=trim(strg)
                 else
                   strg1=trim(line)
                 end if
                 read(strg1,*) (work2(3*natom,nu),nu=1,3*natom)
                 total_atmfrc(:,:,:,:,irpt2) = reshape(work2,(/3,natom,3,natom/))
                 ABI_DEALLOCATE(work2)
               else
                 ABI_ALLOCATE(work2,(3*natom,3*natom))
                 do mu=1,3*natom
                   read(funit,*)(work2(mu,nu),nu=1,3*natom)
                 end do
                 total_atmfrc(:,:,:,:,irpt2) = reshape(work2,(/3,natom,3,natom/))
                 ABI_DEALLOCATE(work2)
               end if
             end if
             if ((line(1:5)=='<cell')) then
               call rdfromline_value('cell',line,strg)
               if (strg/="") then
                 strg1=trim(strg)
                 read(strg1,*)(cell_total(mu,irpt2),mu=1,3)
               else
                 read(funit,*)(cell_total(mu,irpt2),mu=1,3)
               end if
               found2 = .TRUE.
               cycle
             end if
           end do
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
             read(funit,'(a)',iostat=ios) readline
             call rmtabfromline(readline)
             line=adjustl(readline)
             call rdfromline_value('dynamical_matrix',line,strg)
             if (strg/="") then
               strg1=trim(strg)
             else
               strg1=trim(line)
             end if
             read(strg1,*) (work2(nu,3*natom),nu=1,3*natom)
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
           read(funit,'(a)',iostat=ios) readline
           call rdfromline("voigt",line,strg)
           strg1=trim(strg)
           read(strg1,*) voigt
           voigt = voigt + 1 ! 0 to 5 in the xml
           has_straincoupling = .true.
           irpt = 1
         end if

       else
!        Now treat the strain phonon coupling part
         if ((line(1:16)=='<strain_coupling')) then
           read(funit,'(a)',iostat=ios) readline
           call rdfromline("voigt",line,strg)
           strg1=trim(strg)
           read(strg1,*) voigt
           voigt = voigt + 1 ! 0 to 5 in the xml
           irpt = 1
           cycle
         end if

         if(voigt>6)then
           write(message, '(4a)' )ch10,&
&               ' WARNING: the number of strain phonon coupling is superior to 6 in ',filename,ch10
           call wrtout(std_out,message,'COLL')
           exit
         end if

         if ((line(1:22)=='<correction_force unit')) then
           call rdfromline_value('correction_force',line,strg)
           if (strg/="") then
             ABI_ALLOCATE(work2,(3,natom))
             strg1=trim(strg)
             read(strg1,*) (work2(nu,1),nu=1,3)
             do mu=2,natom-1
               read(funit,*)(work2(nu,mu),nu=1,3)
             end do
             read(funit,'(a)',iostat=ios) readline
             call rmtabfromline(readline)
             line=adjustl(readline)
             call rdfromline_value('correction_force',line,strg)
             if (strg/="") then
               strg1=trim(strg)
               read(strg1,*) (work2(nu,natom),nu=1,3)
             else
               strg1=trim(line)
               read(strg1,*) (work2(nu,natom),nu=1,3)
             end if
             strain_coupling(voigt,:,:) = work2(:,:)
             ABI_DEALLOCATE(work2)
           else
             ABI_ALLOCATE(work2,(3,natom))
             do mu=1,natom
               read(funit,*)(work2(nu,mu),nu=1,3)
             end do
             strain_coupling(voigt,:,:) = work2(:,:)
             ABI_DEALLOCATE(work2)
           end if
         end if

         if ((line(1:11)=='<elastic3rd')) then
           call rdfromline_value('elastic3rd',line,strg)
           if (strg/="") then
             strg1=trim(strg)
             read(strg1,*) (elastic3rd(voigt,mu,1),mu=1,6)
             do nu=2,5
               read(funit,*) (elastic3rd(voigt,mu,nu),mu=1,6)
             end do
           else
             do nu=1,5
               read(funit,*) (elastic3rd(voigt,mu,nu),mu=1,6)
             end do
           end if
           read(funit,'(a)',iostat=ios) readline
           call rmtabfromline(readline)
           line=adjustl(readline)
           call rdfromline_value('elastic3rd',line,strg)
           if (strg/="") then
             strg1=trim(strg)
             read(strg1,*) (elastic3rd(voigt,mu,6),mu=1,6)
           else
             strg1=trim(line)
             read(strg1,*) (elastic3rd(voigt,mu,6),mu=1,6)
           end if
           has_anharmonics = .true.
           cycle
         end if

         if ((line(1:29)=='<correction_strain_force unit')) then
           call rdfromline_value('correction_strain_force',line,strg)
           if (strg/="") then
             ABI_ALLOCATE(work2,(3*6,natom))
             strg1=trim(strg)
             read(strg1,*) (work2(nu,1),nu=1,3*6)
             do mu=2,natom-1
               read(funit,*)(work2(nu,mu),nu=1,3*6)
             end do
             read(funit,'(a)',iostat=ios) readline
             call rmtabfromline(readline)
             line=adjustl(readline)
             call rdfromline_value('correction_strain_force',line,strg)
             if (strg/="") then
               strg1=trim(strg)
               read(strg1,*) (work2(nu,natom),nu=1,3*6)
             else
               strg1=trim(line)
               read(strg1,*) (work2(nu,natom),nu=1,3*6)
             end if
             elastic_displacement(voigt,:,:,:) = reshape(work2(:,:),(/6,3,natom/))
             ABI_DEALLOCATE(work2)
           else
             ABI_ALLOCATE(work2,(3*6,natom))
             do mu=1,natom
               read(funit,*)(work2(nu,mu),nu=1,3*6)
             end do
             elastic_displacement(voigt,:,:,:) = reshape(work2(:,:),(/6,3,natom/))
             ABI_DEALLOCATE(work2)
           end if
         end if

         if ((line(1:26)=='<correction_force_constant')) then
           found2=.false.
           do while (.not.found2)
             read(funit,'(a)',iostat=ios) readline
             call rmtabfromline(readline)
             line=adjustl(readline)
             if ((line(1:5)=='<data')) then
               call rdfromline_value('data',line,strg)
               if (strg/="") then
                 ABI_ALLOCATE(work2,(3*natom,3*natom))
                 strg1=trim(strg)
                 read(strg1,*) (work2(1,nu),nu=1,3*natom)
                 do mu=2,3*natom-1
                   read(funit,*)(work2(mu,nu),nu=1,3*natom)
                 end do
                 read(funit,'(a)',iostat=ios) readline
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
                 phonon_strain(voigt)%atmfrc(:,:,:,:,irpt) = &
&                           reshape(work2,(/3,natom,3,natom/))
                 ABI_DEALLOCATE(work2)
               else
                 ABI_ALLOCATE(work2,(3*natom,3*natom))
                 do mu=1,3*natom
                   read(funit,*)(work2(mu,nu),nu=1,3*natom)
                 end do
                 phonon_strain(voigt)%atmfrc(:,:,:,:,irpt) =&
&              reshape(work2,(/3,natom,3,natom/))
                 ABI_DEALLOCATE(work2)
               end if
               has_anharmonics = .true.
             end if
             if ((line(1:5)=='<cell')) then
               call rdfromline_value('cell',line,strg)
               if (strg/="") then
                 strg1=trim(strg)
                 read(strg1,*)(phonon_strain(voigt)%cell(mu,irpt),mu=1,3)
               else
                 read(funit,*)(phonon_strain(voigt)%cell(mu,irpt),mu=1,3)
               end if
               irpt = irpt + 1
               found2=.true.
               cycle
             end if
           end do
         end if

         if ((line(1:17)==char(60)//char(47)//'strain_coupling')) then
!          set nrpt for the previous value of strain
           phonon_strain(voigt)%nrpt = irpt - 1
!         restart the calculation of nrpt
         end if
       end if
     end if
   end do


! Reorder the ATMFRC
! Case 1: only local in the xml
   if (irpt1>0 .and. irpt2==0) then
     ifcs%cell(:,:) = int(cell_local(:,:))
     ifcs%atmfrc(:,:,:,:,:)  = local_atmfrc(:,:,:,:,:)
     ifcs%short_atmfrc(:,:,:,:,:) = local_atmfrc(:,:,:,:,:)
     ifcs%ewald_atmfrc(:,:,:,:,:) = zero

! Case 2: only total in the xml
   else if(irpt1==0 .and. irpt2>0)then
     ifcs%cell(:,:) = int(cell_total(:,:))
     ifcs%atmfrc(:,:,:,:,:)  = total_atmfrc(:,:,:,:,:)
     ifcs%short_atmfrc(:,:,:,:,:) = zero
     ifcs%ewald_atmfrc(:,:,:,:,:) = total_atmfrc(:,:,:,:,:)

! Case 3: local + total in the xml
   else if (irpt1>0 .and. irpt2>0)then
     if(irpt1 <= irpt2)then
       irpt3 = 0
       do ii=1,irpt2
         ifcs%cell(:,ii) = int(cell_total(:,ii))
         ifcs%atmfrc(:,:,:,:,ii)  = total_atmfrc(:,:,:,:,ii)
         do jj=1,irpt1
           if (all(abs(int(cell_local(:,jj))-ifcs%cell(:,ii))<tol16)) then
             ifcs%short_atmfrc(:,:,:,:,ii) = local_atmfrc(:,:,:,:,jj)
             irpt3 = irpt3 + 1
           end if
         end do
       end do
       if(irpt3 /= irpt1)then
         write(message, '(4a)' )ch10,&
&         ' There is several similar short IFC in ',filename,ch10
         MSG_BUG(message)
       end if
     else
       write(message, '(2a,I5,3a,I5,5a)' )ch10,&
&     ' The number of total IFC  (',irpt2,') is inferior to  ',ch10,&
&     ' the number of short range IFC (',irpt1,') in ',filename,ch10,&
&     ' This is not possible',ch10
       MSG_BUG(message)
     end if
   end if

!  Do some checks
   if (any(typat==0)) then
     write(message, '(a,a,a)' )&
&      ' Unable to read the type of atoms ',trim(filename),ch10
     MSG_ERROR(message)
   end if

   if (any(abs(znucl)<tol16)) then
     write(message, '(a,a,a)' )&
&      ' Unable to read the atomic number ',trim(filename),ch10
     MSG_ERROR(message)
   end if

   if (any(abs(all_amu)<tol16)) then
     write(message, '(a,a,a)' )&
&     ' Unable to read the atomic mass ',trim(filename),ch10
     MSG_ERROR(message)
   end if

   close(unit=funit)

#endif

 end if !End if master

!MPI BROADCAST
 call xmpi_bcast(energy,master, comm, ierr)
 call xmpi_bcast(all_amu,master, comm, ierr)
 call xmpi_bcast(dynmat,master, comm, ierr)
 call xmpi_bcast(elastic_constants,master, comm, ierr)
 call xmpi_bcast(epsilon_inf,master, comm, ierr)
 call xmpi_bcast(ifcs%nrpt,master, comm, ierr)
 call xmpi_bcast(ifcs%atmfrc,master, comm, ierr)
 call xmpi_bcast(ifcs%cell,master, comm, ierr)
 call xmpi_bcast(ifcs%ewald_atmfrc,master, comm, ierr)
 call xmpi_bcast(ifcs%short_atmfrc,master, comm, ierr)
 call xmpi_bcast(strain_coupling,master, comm, ierr)
 call xmpi_bcast(phfrq,master, comm, ierr)
 call xmpi_bcast(qph1l,master, comm, ierr)
 call xmpi_bcast(typat,master, comm, ierr)
 call xmpi_bcast(rprimd,master, comm, ierr)
 call xmpi_bcast(xcart,master, comm, ierr)
 call xmpi_bcast(zeff,master, comm, ierr)
 call xmpi_bcast(znucl,master, comm, ierr)
 do ii = 1,6
   call xmpi_bcast(phonon_strain(ii)%nrpt   ,master, comm, ierr)
   call xmpi_bcast(phonon_strain(ii)%atmfrc ,master, comm, ierr)
   call xmpi_bcast(phonon_strain(ii)%cell   ,master, comm, ierr)
 end do
 call xmpi_bcast(elastic3rd   ,master, comm, ierr)
 call xmpi_bcast(elastic_displacement ,master, comm, ierr)
 call xmpi_bcast(has_anharmonics ,master, comm, ierr)

!Fill somes others variables
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 call xcart2xred(natom,rprimd,xcart,xred)

!Re-generate symmetry operations from the lattice and atomic coordinates
 tolsym=tol8
 msym = 384
 ABI_ALLOCATE(spinat,(3,natom))
 ABI_ALLOCATE(ptsymrel,(3,3,msym))
 ABI_ALLOCATE(symafm,(msym))
 ABI_ALLOCATE(symrel,(3,3,msym))
 ABI_ALLOCATE(tnons,(3,msym))
 use_inversion=1
 spinat = 0;
 symrel = 0;
 symafm = 0;
 tnons = 0 ;
 space_group = 0;
 call symlatt(bravais,msym,nptsym,ptsymrel,rprimd,tolsym)
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 call symfind(0,(/zero,zero,zero/),gprimd,0,msym,natom,0,nptsym,nsym,&
&  0,0,ptsymrel,spinat,symafm,symrel,tnons,tolsym,typat,use_inversion,xred)

!Initialisation of crystal
 npsp = ntypat; timrev = 1
 ABI_ALLOCATE(title, (ntypat))
 do ii=1,ntypat
   write(title(ii),'(a,i0)')"No title for typat ",ii
 end do

!Warning znucl is dimension with ntypat = nspsp hence alchemy is not supported here
 call crystal_init(all_amu,Crystal,space_group,natom,npsp,ntypat,nsym,rprimd,typat,xred,&
&  zion,znucl,timrev,.FALSE.,.FALSE.,title,&
&  symrel=symrel(:,:,1:nsym),tnons=tnons(:,1:nsym),symafm=symafm(1:nsym))

!amu is not fill in crystal_init...
 Crystal%amu(:) = all_amu(:)

 ABI_DEALLOCATE(symrel)
 ABI_DEALLOCATE(symafm)
 ABI_DEALLOCATE(tnons)
 ABI_DEALLOCATE(spinat)
 ABI_DEALLOCATE(ptsymrel)

!if strcpling is set to 0 by the user, need to set the flag to false for
!the initialisation of the effective potential
 if (present(strcpling))then
   if(strcpling == 0 )then
     has_anharmonics = .FALSE.
   end if
 end if

!Initialisation of eff_pot
 call effective_potential_init(crystal,eff_pot,energy,ifcs,ncoeff,nph1l,comm,&
&                              dynmat=dynmat,elastic_constants=elastic_constants,&
&                              elastic3rd=elastic3rd,elastic_displacement=elastic_displacement,&
&                              epsilon_inf=epsilon_inf,strain_coupling=strain_coupling,&
&                              phonon_strain=phonon_strain,phfrq=phfrq,qpoints=qph1l,&
&                              has_anharmonicsTerms=has_anharmonics,zeff=zeff)

!DEALLOCATION OF ARRAYS
 ABI_DEALLOCATE(all_amu)
 ABI_DEALLOCATE(cell_local)
 ABI_DEALLOCATE(cell_total)
 ABI_DEALLOCATE(total_atmfrc)
 ABI_DEALLOCATE(local_atmfrc)
 ABI_DEALLOCATE(ifcs%atmfrc)
 ABI_DEALLOCATE(ifcs%cell)
 ABI_DEALLOCATE(ifcs%short_atmfrc)
 ABI_DEALLOCATE(ifcs%ewald_atmfrc)
 ABI_DEALLOCATE(dynmat)
 ABI_DEALLOCATE(strain_coupling)
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
   phonon_strain(ii)%cell   = 0
   ABI_DEALLOCATE(phonon_strain(ii)%atmfrc)
   ABI_DEALLOCATE(phonon_strain(ii)%cell)
 end do
 ABI_FREE(phonon_strain)
 ABI_DEALLOCATE(elastic_displacement)

!DEALLOCATION OF TYPES
 call ifcs%free()
 call crystal%free()

end subroutine system_xml2effpot
!!***

!!****f* m_effective_potential_file/system_ddb2effpot
!!
!! NAME
!! system_ddb2effpot
!!
!! FUNCTION
!!  Transfert ddb into effective potential structure.
!!  Also calculate the IFC
!!
!! INPUTS
!! crytal<type(crystal_t)> = datatype with all the information for the crystal
!! ddb<type(ddb_type)> = datatype with the ddb
!! inp<type(multibinit_dtset_type)> = datatype with the input variables of multibinit
!! comm = MPI communicator
!!
!! OUTPUT
!! effective_potantial<type(effective_potential_type)> = effective_potential datatype to be initialized
!!
!! PARENTS
!!      m_effective_potential_file
!!
!! CHILDREN
!!
!! SOURCE

subroutine system_ddb2effpot(crystal,ddb, effective_potential,inp,comm)

 use defs_basis
 use m_errors
 use m_abicore
 use m_dynmat
 use m_xmpi

 use m_ddb
 use m_ifc
 use m_copy,            only : alloc_copy
 use m_crystal,         only : crystal_t
 use m_multibinit_dataset, only : multibinit_dtset_type
 use m_effective_potential, only : effective_potential_type, effective_potential_free

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
!arrays
 type(ddb_type),intent(inout) :: ddb
 type(effective_potential_type), intent(inout) :: effective_potential
 type(crystal_t),intent(in) :: crystal
 type(multibinit_dtset_type),intent(in) :: inp

!Local variables-------------------------------
!scalar
 real(dp):: wcount1,wcount2
 integer :: chneut,i1,i2,i3,ia,ib,iblok,idir1,idir2,ierr,ii,ipert1,iphl1
 integer :: ipert2,irpt,irpt2,ivarA,ivarB,max1,max2,max3,min1,min2,min3
 integer :: msize,mpert,natom,nblok,nrpt_new,nrpt_new2,rftyp,selectz
 integer :: my_rank,nproc,prt_internalstr
 logical :: iam_master=.FALSE.
 integer,parameter :: master=0
 integer :: nptsym,nsym
 integer :: msym = 384,  use_inversion = 1, space_group
 real(dp):: max_phfq,tolsym = tol8
!arrays
 integer :: bravais(11),cell_number(3),cell2(3)
 integer :: shift(3),rfelfd(4),rfphon(4),rfstrs(4)
 integer,allocatable :: cell_red(:,:)
 real(dp):: dielt(3,3),elast_clamped(6,6),fact
 real(dp):: red(3,3),qphnrm(3),qphon(3,3)
 real(dp),allocatable :: blkval(:,:,:,:,:,:),d2asr(:,:,:,:,:)
 real(dp),allocatable :: instrain(:,:),zeff(:,:,:)
 real(dp),pointer :: atmfrc_red(:,:,:,:,:),wghatm_red(:,:,:)
 character(len=500) :: message
 type(asrq0_t) :: asrq0
 type(ifc_type) :: ifc
 real(dp),allocatable :: d2cart(:,:,:,:,:),displ(:)
 real(dp),allocatable :: eigval(:,:),eigvec(:,:,:,:,:),phfrq(:)
 real(dp),allocatable :: spinat(:,:),tnons(:,:)
 integer,allocatable  :: symrel(:,:,:),symafm(:),ptsymrel(:,:,:)

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
  blkval = 0
  blkval = reshape(ddb%val,(/2,3,mpert,3,mpert,nblok/))

!**********************************************************************
! Transfert crystal values
!**********************************************************************
! Re-generate symmetry operations from the lattice and atomic coordinates
  ABI_ALLOCATE(spinat,(3,natom))
  ABI_ALLOCATE(ptsymrel,(3,3,msym))
  ABI_ALLOCATE(symafm,(msym))
  ABI_ALLOCATE(symrel,(3,3,msym))
  ABI_ALLOCATE(tnons,(3,msym))
  spinat = zero;  symrel = 0;  symafm = 0;  tnons = zero ; space_group = 0;
  call symlatt(bravais,msym,nptsym,ptsymrel,crystal%rprimd,tolsym)
  call symfind(0,(/zero,zero,zero/),crystal%gprimd,0,msym,crystal%natom,0,nptsym,nsym,&
&              0,0,ptsymrel,spinat,symafm,symrel,tnons,tolsym,&
&              crystal%typat,use_inversion,crystal%xred)
  if(crystal%nsym/=nsym)then
    write(message,'(4a,I0,3a,I0,3a)') ch10,&
&          ' --- !WARNING:',ch10,&
&          '     There is ',nsym,' found for the crystal',ch10,&
&          '     but ',crystal%nsym,' found in the DDB',ch10,&
&          ' ---'
      call wrtout(std_out,message,'COLL')
  end if
  call crystal_init(ddb%amu,effective_potential%crystal,&
&                   space_group,crystal%natom,crystal%npsp,&
&                   crystal%ntypat,nsym,crystal%rprimd,&
&                   crystal%typat,crystal%xred,crystal%zion,&
&                   crystal%znucl,crystal%timrev,crystal%use_antiferro,&
&                   .FALSE.,crystal%title,&
&                   symrel=symrel,tnons=tnons,&
&                   symafm=symafm)

  ABI_DEALLOCATE(spinat)
  ABI_DEALLOCATE(ptsymrel)
  ABI_DEALLOCATE(symafm)
  ABI_DEALLOCATE(symrel)
  ABI_DEALLOCATE(tnons)

!**********************************************************************
! Transfert energy from input file
!**********************************************************************
  write(message, '(2a,(80a),6a)') ch10,('=',ii=1,80),ch10,ch10,&
&     ' Extraction of the energy of the structure (unit: Hartree)',ch10
  call wrtout(std_out,message,'COLL')
  call wrtout(ab_out,message,'COLL')
  if (ddb%get_etotal(effective_potential%energy) == 0) then
    if(abs(inp%energy_reference) < tol16)then
      write(message,'(5a)')&
&      ' Warning : Energy of the reference structure is not specify in',&
&      ' the input file.',ch10,' Energy will set to zero',ch10
      call wrtout(std_out,message,'COLL')
      effective_potential%energy = zero
    else
      effective_potential%energy = inp%energy_reference
    end if
  else
    if(abs(inp%energy_reference) > tol16)then
      write(message,'(6a)')&
&      ' Warning : Energy of the reference structure is specify in',&
&      ' the input file.',ch10,' and in the DDB.',&
&      ' The value of the energy is set with the value from the input file',ch10
      call wrtout(std_out,message,'COLL')
      effective_potential%energy = inp%energy_reference
    end if
  end if
  write(message,'(a,es25.12)') ' Energy  = ',&
&                    effective_potential%energy
  call wrtout(std_out,message,'COLL')
  call wrtout(ab_out,message,'COLL')

!**********************************************************************
! Dielectric Tensor and Effective Charges
!**********************************************************************
  ABI_ALLOCATE(zeff,(3,3,natom))
  ABI_ALLOCATE(effective_potential%harmonics_terms%zeff,(3,3,natom))

  rftyp   = 1 ! Blocks obtained by a non-stationary formulation.
  chneut  = 1 ! The ASR for effective charges is imposed
  selectz = 0 ! No selection of some parts of the effective charge tensor
  iblok = ddb%get_dielt_zeff(crystal,rftyp,chneut,selectz,dielt,zeff)
  if (iblok /=0 .and. maxval(abs(dielt)) < 10000) then
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

  ABI_ALLOCATE(effective_potential%fcart,(3,natom))
  effective_potential%fcart = zero
  effective_potential%strten = zero

  qphon(:,1)=zero
  qphnrm(1)=zero
  rfphon(1:2)=0
  rfelfd(1:2)=0
  rfstrs(1:2)=0
  rftyp=4

  call ddb%get_block(iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)

  if (iblok /=0) then
   if(any(abs(inp%strten_reference)>tol16))then
     write(message,'(10a)') ch10,&
&          ' --- !WARNING:',ch10,&
&          '     The stress tensor of the reference structure is specify in the',ch10,&
&          '     input file and in the DDB. The value of the stress tensor is set',ch10,&
&          '     with the value from the input file',ch10,&
&          ' ---'
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')
     effective_potential%strten = inp%strten_reference
   else
!    firts give the corect stress values store in hartree
!    diagonal parts
     effective_potential%strten(1)=blkval(1,1,natom+3,1,1,iblok) *  crystal%ucvol
     effective_potential%strten(2)=blkval(1,2,natom+3,1,1,iblok) *  crystal%ucvol
     effective_potential%strten(3)=blkval(1,3,natom+3,1,1,iblok) *  crystal%ucvol
!    the shear parts
     effective_potential%strten(4)=blkval(1,1,natom+4,1,1,iblok) *  crystal%ucvol
     effective_potential%strten(5)=blkval(1,2,natom+4,1,1,iblok) *  crystal%ucvol
     effective_potential%strten(6)=blkval(1,3,natom+4,1,1,iblok) *  crystal%ucvol
   end if
!  Get forces
   effective_potential%fcart(:,1:natom) = blkval(1,:,1:natom,1,1,iblok)
 else
   if(all(abs(inp%strten_reference(:))<tol16))then
     write(message,'(8a)') ch10,&
&          ' --- !WARNING:',ch10,&
&          '     The stress tensor of the reference structure is not specify',ch10,&
&          '     The stress tensor will be set to zero',ch10,&
&          ' ---'
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')
     effective_potential%strten = zero
   else
     effective_potential%strten = inp%strten_reference
   end if
 end if

 if(any(abs(effective_potential%strten(:)) >tol16))then
   write(message, '(3a)' )ch10,&
&   ' Cartesian components of forces (hartree/bohr)',ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   do ii = 1, natom
     write(message, '(I4,a,3(e16.8))' ) &
&     ii,'   ',effective_potential%fcart(:,ii)

     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
   end do

   write(message, '(a,a)' )ch10,&
&   ' Cartesian components of stress tensor (hartree/bohr^3)'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '  sigma(1 1)=',effective_potential%strten(1) / crystal%ucvol,&
&   '  sigma(3 2)=',effective_potential%strten(4) / crystal%ucvol
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '  sigma(2 2)=',effective_potential%strten(2) / crystal%ucvol,&
&   '  sigma(3 1)=',effective_potential%strten(5) / crystal%ucvol
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '  sigma(3 3)=',effective_potential%strten(3) / crystal%ucvol,&
&   '  sigma(2 1)=',effective_potential%strten(6) / crystal%ucvol
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a)' ) ' '
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
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
  call ddb%get_block(iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)

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

  call ddb%get_block(iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)

  d2asr = zero
  if (iblok /=0) then
    call asria_calc(inp%asr,d2asr,ddb%val(:,:,iblok),ddb%mpert,ddb%natom)
  end if

  ! Acoustic Sum Rule
  ! In case the interatomic forces are not calculated, the
  ! ASR-correction (asrq0%d2asr) has to be determined here from the Dynamical matrix at Gamma.
  asrq0 = ddb%get_asrq0(inp%asr, inp%rfmeth, crystal%xcart)

!**********************************************************************
! Interatomic Forces Calculation
!**********************************************************************
! ifc to be calculated for interpolation
  write(message, '(a,a,(80a),a,a,a,a)' ) ch10,('=',ii=1,80),ch10,ch10,&
&   ' Calculation of the interatomic forces from DDB',ch10
  call wrtout(std_out,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  call ifc_init(ifc,crystal,ddb,inp%brav,inp%asr,inp%symdynmat,inp%dipdip,inp%rfmeth,&
&   inp%ngqpt(1:3),inp%nqshft,inp%q1shft,dielt,effective_potential%harmonics_terms%zeff,&
&   inp%nsphere,inp%rifcsph,inp%prtsrlr,inp%enunit,comm)

!***************************************************************************
! Interpolation of the dynamical matrix for each qpoint from ifc
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
  effective_potential%harmonics_terms%nqpt         = inp%nph1l
  effective_potential%harmonics_terms%qpoints(:,:) = inp%qph1l(:,:)

! Store the highest frequency
  max_phfq = zero

  do iphl1=1,inp%nph1l

   ! Initialisation of the phonon wavevector
    qphon(:,1)=inp%qph1l(:,iphl1)
    if (inp%nph1l /= 0) qphnrm(1) = inp%qnrml1(iphl1)

    ! Get d2cart using the interatomic forces and the
    ! long-range coulomb interaction through Ewald summation
    call gtdyn9(ddb%acell,ifc%atmfrc,ifc%dielt,ifc%dipdip,ifc%dyewq0,d2cart,crystal%gmet,&
&     ddb%gprim,mpert,natom,ifc%nrpt,qphnrm(1),qphon(:,1),crystal%rmet,ddb%rprim,ifc%rpt,&
&     ifc%trans,crystal%ucvol,ifc%wghatm,crystal%xred,zeff,xmpi_comm_self)

    ! Calculation of the eigenvectors and eigenvalues of the dynamical matrix
    call dfpt_phfrq(ddb%amu,displ,d2cart,eigval,eigvec,crystal%indsym,&
&     mpert,crystal%nsym,natom,crystal%nsym,crystal%ntypat,phfrq,qphnrm(1),qphon,&
&     crystal%rprimd,inp%symdynmat,crystal%symrel,crystal%symafm,crystal%typat,crystal%ucvol)

    ! Write the phonon frequencies
    call dfpt_prtph(displ,inp%eivec,inp%enunit,ab_out,natom,phfrq,qphnrm(1),qphon)

!   Store the highest frequency in cmm-1
    max_phfq = max(maxval(phfrq*Ha_cmm1),max_phfq)

    effective_potential%harmonics_terms%dynmat(:,:,:,:,:,iphl1) = d2cart(:,:,:natom,:,:natom)
    effective_potential%harmonics_terms%phfrq(:,iphl1) = phfrq(:) * Ha_cmm1

  end do

  write(message, '(2a,f15.7,a)' ) ch10,&
&   ' The highest frequency found is ',max_phfq,' cm-1'
  call wrtout(std_out,message,'COLL')

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

! Store the sum of the weight of IFC for the final check
  wcount1 = 0
  do irpt=1,ifc%nrpt
    wcount1 = wcount1 + sum(ifc%wghatm(:,:,irpt))
  end do

!Set the maximum and the miminum for the bound of the cell
  max1 = maxval(ifc%cell(1,:));  min1 = minval(ifc%cell(1,:))
  max2 = maxval(ifc%cell(2,:));  min2 = minval(ifc%cell(2,:))
  max3 = maxval(ifc%cell(3,:));  min3 = minval(ifc%cell(3,:))
  cell_number(1) = max1 - min1 + 1
  cell_number(2) = max2 - min2 + 1
  cell_number(3) = max3 - min3 + 1

! set the new number of cell, sometimes, in canonical coordinates,
! some cell are delete but they exist in reduced coordinates.
  nrpt_new = product(cell_number(:))

! Allocate temporary array
  ABI_ALLOCATE(atmfrc_red,(3,natom,3,natom,nrpt_new))
  ABI_ALLOCATE(wghatm_red,(natom,natom,nrpt_new))
  ABI_ALLOCATE(cell_red,(3,nrpt_new))

  wghatm_red(:,:,:) = zero

  if(iam_master)then
    do ia=1,natom
      do ib=1,natom

!       Simple Lattice
        if (inp%brav==1) then
!          In this case, it is better to work in reduced coordinates
!          As rcan is in canonical coordinates, => multiplication by gprim
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
        shift(:) = int(anint(red(2,:) - crystal%xred(:,ib)) - anint(red(1,:) - crystal%xred(:,ia)))

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
                   atmfrc_red(:,ia,:,ib,irpt2) = ifc%atmfrc(:,ia,:,ib,irpt)
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
  end if

  call xmpi_bcast(atmfrc_red,master, comm, ierr)
  call xmpi_bcast(wghatm_red,master, comm, ierr)
  call xmpi_bcast(cell_red,master, comm, ierr)

!  Copy ifc into effective potential
! !!Warning eff_pot%ifcs only contains atmfrc,short_atmfrc,ewald_atmfrc,,nrpt and cell!!
! rcan,ifc%rpt,wghatm and other quantities
! are not needed for effective potential!!!
  call ifc%free()
  call effective_potential%harmonics_terms%ifcs%free()

! Only conserve the necessary points in rpt
  nrpt_new2 = 0
  do irpt = 1, nrpt_new
    if (abs(sum(wghatm_red(:,:,irpt))) >tol16) then
      nrpt_new2 = nrpt_new2 + 1
    end if
  end do

! Set the new number of rpt
  effective_potential%harmonics_terms%ifcs%nrpt = nrpt_new2

! Allocation of the final arrays
  ABI_ALLOCATE(effective_potential%harmonics_terms%ifcs%atmfrc,(3,natom,3,natom,nrpt_new2))
  ABI_ALLOCATE(effective_potential%harmonics_terms%ifcs%short_atmfrc,(3,natom,3,natom,nrpt_new2))
  ABI_ALLOCATE(effective_potential%harmonics_terms%ifcs%ewald_atmfrc,(3,natom,3,natom,nrpt_new2))
  ABI_ALLOCATE(effective_potential%harmonics_terms%ifcs%cell,(3,nrpt_new2))
  ABI_ALLOCATE(effective_potential%harmonics_terms%ifcs%wghatm,(natom,natom,nrpt_new2))

  irpt2 = 0
  do irpt = 1,nrpt_new
    if (abs(sum(wghatm_red(:,:,irpt))) > tol16) then
      irpt2 = irpt2 + 1
!     Apply weight on each R point
      do ia=1,effective_potential%crystal%natom
        do ib=1,effective_potential%crystal%natom
          atmfrc_red(:,ia,:,ib,irpt) = atmfrc_red(:,ia,:,ib,irpt)*wghatm_red(ia,ib,irpt)
        end do
      end do
      effective_potential%harmonics_terms%ifcs%cell(:,irpt2) = cell_red(:,irpt)
      effective_potential%harmonics_terms%ifcs%atmfrc(:,:,:,:,irpt2) = atmfrc_red(:,:,:,:,irpt)
      if (inp%dipdip == 1) then
        effective_potential%harmonics_terms%ifcs%short_atmfrc(:,:,:,:,irpt2)=&
&                                                                     atmfrc_red(:,:,:,:,irpt)
      else
        effective_potential%harmonics_terms%ifcs%short_atmfrc(:,:,:,:,irpt2) = zero
      end if
      effective_potential%harmonics_terms%ifcs%short_atmfrc(:,:,:,:,irpt2)=atmfrc_red(:,:,:,:,irpt)
      effective_potential%harmonics_terms%ifcs%ewald_atmfrc(:,:,:,:,irpt2) = zero
      effective_potential%harmonics_terms%ifcs%wghatm(:,:,irpt2) =  wghatm_red(:,:,irpt)
    end if
  end do


  ABI_DEALLOCATE(atmfrc_red)
  ABI_DEALLOCATE(wghatm_red)
  ABI_DEALLOCATE(cell_red)

! Final check
  wcount2 = 0
  do irpt = 1, effective_potential%harmonics_terms%ifcs%nrpt
    wcount2 = wcount2 + sum(effective_potential%harmonics_terms%ifcs%wghatm(:,:,irpt))
  end do

  if (abs(wcount1-wcount2)/(wcount1+wcount2)>tol8) then
    write(message,'(2a,es15.4,a,es15.4,a,es15.4)')'The total wghatm has changed',ch10,&
&    wcount1,' before and ', wcount2, ' now, difference being ',wcount1-wcount2
    MSG_BUG(message)
  end if


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
  call ddb%get_block(iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)

  ABI_ALLOCATE(effective_potential%harmonics_terms%strain_coupling,(6,3,natom))
  effective_potential%harmonics_terms%strain_coupling = zero

  if (iblok /=0) then

!   then print the internal strain tensor (only the force one)
    prt_internalstr=1
    call ddb_internalstr(inp%asr,ddb%val,d2asr,iblok,instrain,&
&                        ab_out,mpert,natom,nblok,prt_internalstr)

    do ipert1=1,6
      do ipert2=1,natom
        do idir2=1,3
          ii=3*(ipert2-1)+idir2
            effective_potential%harmonics_terms%strain_coupling(ipert1,idir2,ipert2)=&
&                                                            (-1.0_dp)*instrain(ii,ipert1)
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
  call asrq0%free()

  write(message,'(a)')ch10
  call wrtout(std_out,message,'COLL')
  call wrtout(ab_out,message,'COLL')

end subroutine system_ddb2effpot
!!***

!!****f* m_effective_potential_file/coeffs_xml2effpot
!! NAME
!! coeffs_xml2effpot
!!
!! FUNCTION
!! Open xml file of effective potentiel, then reads the variables
!! and store them in effective potentential type
!!
!! INPUTS
!! filename = path of input or output file
!! comm=MPI communicator
!!
!! OUTPUT
!! eff_pot<type(effective_potential_type)> = effective_potential datatype
!!
!! PARENTS
!!      m_effective_potential_file
!!
!! CHILDREN
!!
!! SOURCE

subroutine coeffs_xml2effpot(eff_pot,filename,comm)

 use m_atomdata
 use m_effective_potential, only : effective_potential_type
 use m_polynomial_coeff
 use m_polynomial_term
 use m_crystal, only : symbols_crystal
#if defined HAVE_XML
 use iso_c_binding, only : C_CHAR,C_PTR,c_f_pointer
#endif

 !Arguments ------------------------------------
 !scalars
 character(len=*),intent(in) :: filename
 integer, intent(in) :: comm
 !arrays
 type(effective_potential_type), intent(inout) :: eff_pot

 !Local variables-------------------------------
 !scalar
 integer :: ii,jj,my_rank,ndisp,ncoeff,nterm_max,nstrain,ndisp_max,nproc,nterm
! character(len=200),allocatable :: name(:)
 character(len=200) :: name
#ifdef HAVE_XML
 integer :: icoeff,iterm
#endif

#ifndef HAVE_XML
 integer :: funit = 1,ios = 0
 integer :: icoeff,idisp,istrain,iterm,mu
 logical :: found,found2,displacement
 character (len=XML_RECL) :: line,readline
 character (len=XML_RECL) :: strg,strg1
#endif
 character(len=500) :: message
 character(len=5),allocatable :: symbols(:)
 integer,parameter :: master=0
 logical :: iam_master
 logical :: debug =.FALSE.
 !arrays
 real(dp),allocatable :: coefficient(:),weight(:,:)
 integer,allocatable :: atindx(:,:,:,:), cell(:,:,:,:,:),direction(:,:,:),power_disp(:,:,:)
 integer,allocatable :: strain(:,:,:),power_strain(:,:,:)
 type(polynomial_coeff_type),dimension(:),allocatable :: coeffs
 type(polynomial_term_type),dimension(:,:),allocatable :: terms

! *************************************************************************


 !Open the atomicdata XML file for reading
 write(message,'(a,a)')'-Opening the file ',filename

 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 iam_master = (my_rank == master)

!Get Dimention of system and allocation/initialisation of array
 ncoeff  = 0
 nterm   = 0
 ndisp   = 0
 nstrain = 0
 call effective_potential_file_getDimCoeff(filename,ncoeff,ndisp_max,nterm_max)

!  Do some checks
 if (nterm_max<=0) then
   write(message, '(a,a,a)' )&
&     ' Unable to read the number of terms in ',trim(filename),ch10
   MSG_ERROR(message)
 end if

  if (ndisp_max<=0) then
    write(message, '(a,a,a)' )&
&    ' Unable to read the number of displacement in ',trim(filename),ch10
    MSG_ERROR(message)
  end if

!Allocation ov the polynomial coeff type
 ABI_DATATYPE_ALLOCATE(coeffs,(ncoeff))

 if(iam_master)then

#if defined HAVE_XML
   write(message,'(3a)')'-Reading the file ',trim(filename),&
&   ' with LibXML library'
#else
   write(message,'(3a)')'-Reading the file ',trim(filename),&
&   ' with Fortran'
#endif
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')


   ABI_ALLOCATE(symbols,(eff_pot%crystal%natom))
!  Get the symbols arrays
   call symbols_crystal(eff_pot%crystal%natom,eff_pot%crystal%ntypat,&
&                       eff_pot%crystal%npsp,symbols,eff_pot%crystal%typat,eff_pot%crystal%znucl)


 !Read with libxml librarie
#if defined HAVE_XML

   ABI_DATATYPE_ALLOCATE(terms,(ncoeff,nterm_max))
   ABI_ALLOCATE(atindx,(ncoeff,nterm_max,2,ndisp_max))
   ABI_ALLOCATE(coefficient,(ncoeff))
   ABI_ALLOCATE(cell,(ncoeff,nterm_max,3,2,ndisp_max))
   ABI_ALLOCATE(direction,(ncoeff,nterm_max,ndisp_max))
   ABI_ALLOCATE(strain,(ncoeff,nterm_max,ndisp_max))
   ABI_ALLOCATE(power_disp,(ncoeff,nterm_max,ndisp_max))
   ABI_ALLOCATE(power_strain,(ncoeff,nterm_max,ndisp_max))
   ABI_ALLOCATE(weight,(ncoeff,nterm_max))

!  Read the values of this term with libxml
   call effpot_xml_readCoeff(char_f2c(trim(filename)),ncoeff,ndisp_max,nterm_max,&
&                            coefficient,atindx,cell,direction,power_disp,power_strain,&
&                            strain,weight)
!  In the XML the atom index begin to zero
!  Need to shift for fortran array
   atindx(:,:,:,:) = atindx(:,:,:,:) + 1

   do icoeff=1,ncoeff
     do iterm=1,nterm_max
!      Initialisation of the polynomial_term structure with the values from the
       call polynomial_term_init(atindx(icoeff,iterm,:,:),cell(icoeff,iterm,:,:,:),&
&                                direction(icoeff,iterm,:),ndisp_max,ndisp_max,terms(icoeff,iterm),&
&                                power_disp(icoeff,iterm,:),power_strain(icoeff,iterm,:),&
&                                strain(icoeff,iterm,:),weight(icoeff,iterm),check=.true.)
     end do
!    Initialisation of the polynomial_coefficent structure with the values
     call polynomial_coeff_init(coefficient(icoeff),nterm_max,coeffs(icoeff),&
&                               terms(icoeff,:),check=.true.)

!    Get the name of this coefficient  and set it
!    Try to find the index of the term corresponding to the interation in the
!    reference cell (000) in order to compute the name correctly...
!    If this coeff is not in the ref cell, take by default the first term
     if(coeffs(icoeff)%nterm > 0)then
       call polynomial_coeff_getName(name,coeffs(icoeff),symbols,recompute=.true.)
       call polynomial_coeff_setName(name,coeffs(icoeff))
     end if

!    Free them all
     do iterm=1,nterm_max
       call polynomial_term_free(terms(icoeff,iterm))
     end do
   end do

#else
   ABI_DATATYPE_ALLOCATE(terms,(1,nterm_max))
   ABI_ALLOCATE(atindx,(1,1,2,ndisp_max))
   ABI_ALLOCATE(coefficient,(1))
   ABI_ALLOCATE(cell,(1,1,3,2,ndisp_max))
   ABI_ALLOCATE(direction,(1,1,ndisp_max))
   ABI_ALLOCATE(strain,(1,1,ndisp_max))
   ABI_ALLOCATE(power_disp,(1,1,ndisp_max))
   ABI_ALLOCATE(power_strain,(1,1,ndisp_max))
   ABI_ALLOCATE(weight,(1,1))
!  Loop over the file
!  Read the values of all the terms with fortran
   if (open_file(filename,message,unit=funit,form="formatted",&
&              status="old",action="read") /= 0) then
     MSG_ERROR(message)
   end if

!    Start a reading loop in fortran
     rewind(unit=funit)
     ios  = 0
     found=.false.

!    Initialisation of counter
     icoeff  = 0

!    Parser
     do while (ios==0)
       read(funit,'(a)',iostat=ios) readline
       if (ios == 0) then
         call rmtabfromline(readline)
         line=adjustl(readline)
         if ((line(1:12)==char(60)//'coefficient')) then
!          Read headers of coefficient
           call rdfromline('text',line,strg)
           if (strg/="") then
             name=trim(strg)
           end if
           call rdfromline('value',line,strg)
           if (strg/="") then
             strg1=trim(strg)
             read(strg1,*) coefficient(1)
           else
             coefficient(1) = zero
           end if
!          End read headers of coefficient
!          Reset counter
           found  = .false.
           atindx = 0;  cell   = 0 ;  direction = 0
           strain = 0; power_strain = 0;  power_disp  = 0
           iterm   = 0
           idisp   = 0
           istrain = 0
           nterm   = 0
           do while (.not.found)
             read(funit,'(a)',iostat=ios) readline
             call rmtabfromline(readline)
             line=adjustl(readline)
             if ((line(1:13)==char(60)//'/coefficient')) then
               found= .true.
               cycle
             end if
             if ((line(1:5)==char(60)//'term')) then
               nterm = nterm + 1
               ndisp = 0
               nstrain = 0
               idisp = 0
               istrain = 0
               displacement = .true.
               call rdfromline('weight',line,strg)
               if (strg/="") then
                 strg1=trim(strg)
                 read(strg1,*) weight
               end if
               do while(displacement)
                 read(funit,'(a)',iostat=ios) readline
                 call rmtabfromline(readline)
                 line=adjustl(readline)
                 if ((line(1:6)==char(60)//'/term')) then
                   displacement = .false.
                 end if
                 if ((line(1:7)==char(60)//'strain')) then
                   nstrain = nstrain + 1
                   istrain = istrain + 1
                   call rdfromline('power',line,strg)
                   if (strg/="") then
                     strg1=trim(strg)
                     read(strg1,*) power_strain(1,1,istrain)
                   end if
                   call rdfromline('voigt',line,strg)
                   if (strg/="") then
                     strg1=trim(strg)
                     read(strg1,*) strain(1,1,istrain)
                   end if
                 end if
                 if ((line(1:18)==char(60)//'displacement_diff')) then
                   ndisp = ndisp + 1
                   idisp = idisp + 1
                   found2=.true.
                   call rdfromline('atom_a',line,strg)
                   if (strg/="") then
                     strg1=trim(strg)
                     read(strg1,*) atindx(1,1,1,idisp)
                   end if
                   call rdfromline('atom_b',line,strg)
                   if (strg/="") then
                     strg1=trim(strg)
                     read(strg1,*) atindx(1,1,2,idisp)
                   end if
                   call rdfromline('direction',line,strg)
                   if (strg/="") then
                     strg1=trim(strg)
                     if (trim(strg1).eq."x") direction(1,1,idisp) = 1
                     if (trim(strg1).eq."y") direction(1,1,idisp) = 2
                     if (trim(strg1).eq."z") direction(1,1,idisp) = 3
                   end if
                   call rdfromline('power',line,strg)
                   if (strg/="") then
                     strg1=trim(strg)
                     read(strg1,*) power_disp(1,1,idisp)
                   end if
                   do while(found2)
                     read(funit,'(a)',iostat=ios) readline
                     call rmtabfromline(readline)
                     line=adjustl(readline)
                     if ((line(1:7)==char(60)//'cell_a')) then
                       call rdfromline_value('cell_a',line,strg)
                       if (strg/="") then
                         strg1=trim(strg)
                         read(strg1,*) (cell(1,1,mu,1,idisp),mu=1,3)
                       else
                         read(funit,'(a)',iostat=ios) readline
                         call rmtabfromline(readline)
                         line=adjustl(readline)
                         call rdfromline_value('cell_a',line,strg)
                         if (strg/="") then
                           strg1=trim(strg)
                           read(strg1,*)(cell(1,1,mu,1,idisp),mu=1,3)
                         else
                           strg1=trim(line)
                           read(strg1,*)(cell(1,1,mu,1,idisp),mu=1,3)
                         end if
                       end  if
                     end if
                     if ((line(1:7)==char(60)//'cell_b')) then
                       call rdfromline_value('cell_b',line,strg)
                       if (strg/="") then
                         strg1=trim(strg)
                         read(strg1,*) (cell(1,1,mu,2,idisp),mu=1,3)
                       else
                         read(funit,'(a)',iostat=ios) readline
                         call rmtabfromline(readline)
                         line=adjustl(readline)
                         call rdfromline_value('cell_b',line,strg)
                         if (strg/="") then
                           strg1=trim(strg)
                           read(strg1,*)(cell(1,1,mu,2,idisp),mu=1,3)
                         else
                           strg1=trim(line)
                           read(strg1,*)(cell(1,1,mu,2,idisp),mu=1,3)
                         end if
                       end  if
                     end if
                     if ((line(1:19)==char(60)//'/displacement_diff')) then
                       found2=.false.
                     end if
                   end do
                 end if
               end do!end do while displacement
!              In the XML the atom index begin to zero
!              Need to shift for fortran array
               atindx(1,1,:,:) = atindx(1,1,:,:) + 1
!              Initialisation of the polynomial_term structure with the values from the
!              previous step
               iterm = iterm + 1
               call polynomial_term_init(atindx(1,1,:,:),cell(1,1,:,:,:),&
&                                        direction(1,1,:),ndisp,nstrain,terms(1,iterm),&
&                                        power_disp(1,1,:),power_strain(1,1,:),&
&                                        strain(1,1,:),weight(1,1),check=.true.)
             end if!end if term
           end do!end do while found (coeff)

!          Initialisation of the polynomial_coefficent structure with the values from the
!          previous step
           icoeff = icoeff + 1
           call polynomial_coeff_init(coefficient(1),nterm,coeffs(icoeff),terms(1,:))
           call polynomial_coeff_getName(name,coeffs(icoeff),symbols,recompute=.true.)
           call polynomial_coeff_setName(name,coeffs(icoeff))
!          Deallocation of the terms array for this coefficient
           do jj=1,nterm_max
             call polynomial_term_free(terms(1,jj))
           end do
         end if!end if line = coefficient
       end if!end if ios==0
     end do!end do while on file

     close(unit=funit)

#endif
     ABI_DATATYPE_DEALLOCATE(terms)
     ABI_DEALLOCATE(atindx)
     ABI_DEALLOCATE(coefficient)
     ABI_DEALLOCATE(cell)
     ABI_DEALLOCATE(direction)
     ABI_DEALLOCATE(strain)
     ABI_DEALLOCATE(power_disp)
     ABI_DEALLOCATE(power_strain)
     ABI_DEALLOCATE(weight)
     ABI_DEALLOCATE(symbols)
   end if !End if master

!9-MPI BROADCAST
 do ii=1,ncoeff
   call polynomial_coeff_broadcast(coeffs(ii),master, comm)
 end do

!10-checks

!11-debug print
 if(debug)then
   do ii=1,ncoeff
     do jj=1,coeffs(ii)%nterm
#if defined HAVE_XML
       write(200+my_rank,*)"ii,jj,ndisp,nterm",ii,jj,coeffs(ii)%nterm,coeffs(ii)%terms(jj)%ndisp
       write(200+my_rank,*)"atindx",coeffs(ii)%terms(jj)%atindx
       write(200+my_rank,*)"cell1",coeffs(ii)%terms(jj)%cell(:,1,:)
       write(200+my_rank,*)"cell2",coeffs(ii)%terms(jj)%cell(:,2,:)
       write(200+my_rank,*)"direction",coeffs(ii)%terms(jj)%direction
       write(200+my_rank,*)"power_disp",coeffs(ii)%terms(jj)%power_disp
       write(200+my_rank,*)"weight",coeffs(ii)%terms(jj)%weight
#else
       write(300+my_rank,*)"ii,jj,ndisp,nterm",ii,jj,coeffs(ii)%nterm,coeffs(ii)%terms(jj)%ndisp
       write(300+my_rank,*)"atindx",coeffs(ii)%terms(jj)%atindx
       write(300+my_rank,*)"cell1",coeffs(ii)%terms(jj)%cell(:,1,:)
       write(300+my_rank,*)"cell2",coeffs(ii)%terms(jj)%cell(:,2,:)
       write(300+my_rank,*)"direction",coeffs(ii)%terms(jj)%direction
       write(300+my_rank,*)"power_disp",coeffs(ii)%terms(jj)%power_disp
       write(300+my_rank,*)"weight",coeffs(ii)%terms(jj)%weight
#endif
     end do
   end do
#if defined HAVE_XML
   close(200+my_rank)
#else
   close(300+my_rank)
#endif
 end if

!12-Initialisation of eff_pot
 call effective_potential_setCoeffs(coeffs,eff_pot,ncoeff)

!13-Deallocation of type
 do ii=1,ncoeff
   call polynomial_coeff_free(coeffs(ii))
 end do
 ABI_DATATYPE_DEALLOCATE(coeffs)

end subroutine coeffs_xml2effpot
!!***

!!****f* m_effective_potential_file/effective_potential_file_readMDfile
!!
!! NAME
!! effective_potential_file_readMDfile
!!
!! FUNCTION
!! Read MD FILE (HIST or ASCII)
!!
!! INPUTS
!! filename = path of the file
!! option,optional   = 0 (default), the stress is printed in the MD File
!!                     1, the force on the cell is printed in the MD File (-1 * stress),
!!                        in this case, we multiply the stress by -1 in order to get the stresse
!!
!! OUTPUT
!! hist<type(abihist)> = datatype with the  history of the MD
!!
!! PARENTS
!!      m_effective_potential_file,multibinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine effective_potential_file_readMDfile(filename,hist,option)

!Arguments ------------------------------------
!scalars
 integer,optional :: option
!arrays
 type(abihist),intent(inout) :: hist
 character(len=fnlen),intent(in) :: filename
!Local variables-------------------------------
!scalar
 integer :: ia,ii,mu,nu,natom,nstep,type,option_in
 integer :: ios=0, unit_md=24
!arrays
 character (len=10000) :: readline,line
 real(dp) :: tmp(6)
 real(dp),allocatable :: xcart(:,:)

! *************************************************************************

 call effective_potential_file_getType(filename,type)

 option_in = 0
 if(present(option))then
   option_in = option
 end if
 if(type==40)then
!  Netcdf type
   call read_md_hist(filename,hist,.FALSE.,.FALSE.,.FALSE.)

 else if(type==41)then

!  ASCII file
   call  effective_potential_file_getDimMD(filename,natom,nstep)

   ii  = 1
   ios = 0

   ABI_ALLOCATE(xcart,(3,natom))
   call abihist_free(hist)
   call abihist_init(hist,natom,nstep,.FALSE.,.FALSE.)

!  Start a reading loop in fortran
   rewind(unit=unit_md)
   do while ((ios==0).and.ii<=nstep)
     read(unit_md,'(a)',iostat=ios) readline
     read(unit_md,'(a)',iostat=ios) readline
     line=adjustl(readline)
     read(line,*) hist%etot(ii)
     hist%etot(ii) = hist%etot(ii)
     do mu=1,3
       read(unit_md,'(a)',iostat=ios) readline
       line=adjustl(readline)
       read(line,*) (hist%rprimd(nu,mu,ii),nu=1,3)
     end do
     do ia=1,natom
       read(unit_md,'(a)',iostat=ios) readline
       line=adjustl(readline)
       read(line,*) (tmp(mu),mu=1,6)
       xcart(:,ia) = tmp(1:3)
       hist%fcart(:,ia,ii) = tmp(4:6)
     end do
     call xcart2xred(natom,hist%rprimd(:,:,ii),xcart(:,:),hist%xred(:,:,ii))
     read(unit_md,'(a)',iostat=ios) readline
     line=adjustl(readline)
     read(line,*) (hist%strten(mu,ii),mu=1,6)
     ii = ii + 1
   end do
   do ii=1,nstep
     do mu=1,3
       hist%acell(mu,:) = hist%rprimd(mu,mu,ii)
     end do
   end do
   close(unit_md)
   ABI_DEALLOCATE(xcart)

 end if!end if type

   if((type==40 .or. type==41).and.option == 1)then
!    multiply by -1 if the current strten -1*stress, we need only stress...
     hist%strten(:,:) = -1 * hist%strten(:,:)
   end if



end subroutine effective_potential_file_readMDfile
!!***

!!****f* m_effective_potential_file/effective_potential_file_mapHistToRef
!!
!! NAME
!! effective_potential_file_mapHistToRef
!!
!! FUNCTION
!! Generate the supercell in the effective potential according to the size of the
!! supercell in the hist file
!! Check if the hist file match to reference supercell in the effective potential
!! If not, the hist file is reordering
!!
!! INPUTS
!! eff_pot<type(effective_potential)> = effective potential
!! hist<type(abihist)> = The history of the MD
!! comm = MPI communicator
!!
!! OUTPUT
!! hist<type(abihist)> = The history of the MD
!!
!! PARENTS
!!      m_fit_polynomial_coeff,multibinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine effective_potential_file_mapHistToRef(eff_pot,hist,comm,verbose)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 logical,optional,intent(in) :: verbose
!arrays
 type(effective_potential_type),intent(inout) :: eff_pot
 type(abihist),intent(inout) :: hist
!Local variables-------------------------------
!scalar
 integer :: factE_hist,ia,ib,ii,jj,natom_hist,ncells,nstep_hist
 real(dp):: factor
 logical :: revelant_factor,need_map,need_verbose
!arrays
 real(dp) :: rprimd_hist(3,3),rprimd_ref(3,3)
 integer :: ncell(3),scale_cell(3)
 integer,allocatable  :: shift(:,:)
 integer,allocatable  :: list_map(:) !blkval(:),
 real(dp),allocatable :: xred_ref(:,:) ! xred_hist(:,:),
 real(dp),allocatable :: list_dist(:),list_reddist(:,:),list_absdist(:,:)
 character(len=500) :: msg
 type(abihist) :: hist_tmp
! *************************************************************************

!Set optional values
 need_verbose = .false.
 if (present(verbose)) need_verbose = verbose

 natom_hist = size(hist%xred,2)
 nstep_hist = size(hist%xred,3)

! Try to set the supercell according to the hist file
 rprimd_ref(:,:)  = eff_pot%crystal%rprimd
 rprimd_hist(:,:) = hist%rprimd(:,:,1)

 do ia=1,3
   scale_cell(:) = 0
   do ii=1,3
     if(abs(rprimd_ref(ii,ia)) > tol10)then
       scale_cell(ii) = nint(rprimd_hist(ii,ia) / rprimd_ref(ii,ia))
     end if
   end do
!  Check if the factor for the supercell is revelant
   revelant_factor = .TRUE.
   do ii=1,3
     if(abs(scale_cell(ii)) < tol10) cycle
     factor = abs(scale_cell(ii))
     do jj=ii,3
       if(abs(scale_cell(jj)) < tol10) cycle
       if(abs(abs(scale_cell(ii))-abs(scale_cell(jj))) > tol10) revelant_factor = .FALSE.
     end do
   end do
   if(.not.revelant_factor)then
     write(msg, '(3a)' )&
&         'unable to map the hist file ',ch10,&
&         'Action: check/change your MD file'
     MSG_ERROR(msg)
   else
     ncell(ia) = int(factor)
   end if
 end do

 ncells = product(ncell)

!Check if the energy store in the hist is revelant, sometimes some MD files gives
!the energy of the unit cell... This is not suppose to happen... But just in case...
 do ii=1,nstep_hist
   factE_hist = int(anint(hist%etot(ii) / eff_pot%energy))
   if(factE_hist == 1) then
!    In this case we mutiply the energy of the hist by the number of cell
     hist%etot(ii) = hist%etot(ii)  * ncells
   end if
   if(factE_hist /=1 .and. factE_hist /= ncells)then
     write(msg, '(4a,I0,a,I0,2a,I0,3a,I0,3a)' )ch10,&
&          ' --- !WARNING',ch10,&
&          '     The energy of the step ',ii,' seems to be with multiplicity of ',factE_hist,ch10,&
&          '     However, the multiplicity of the cell is ',ncells,'.',ch10,&
&          '     Please check the energy of the step ',ii,ch10,&
&          ' ---',ch10
     if(need_verbose) call wrtout(std_out,msg,'COLL')
   end if
 end do


!Set the new supercell datatype into the effective potential reference
 call effective_potential_setSupercell(eff_pot,comm,ncell)

!allocation
 ABI_ALLOCATE(shift,(3,natom_hist))
 ABI_ALLOCATE(list_map,(natom_hist))
 ABI_ALLOCATE(list_reddist,(3,natom_hist))
 ABI_ALLOCATE(list_absdist,(3,natom_hist))
 ABI_ALLOCATE(list_dist,(natom_hist))
 ABI_ALLOCATE(xred_ref,(3,natom_hist))

 !Putting maping list to zero
 list_map = 0

 !Fill xcart_ref/hist and xred_ref/hist

 call xcart2xred(eff_pot%supercell%natom,eff_pot%supercell%rprimd,&
&                eff_pot%supercell%xcart,xred_ref)                   ! Get xred_ref


do ia=1,natom_hist !Loop over all reference atoms
   ! Put temporary lists to zero
   list_reddist = 0
   list_absdist = 0
   list_dist = 0
   shift = 0
   do ib=1,natom_hist !Loop over all atoms of distorted structure
      !Calculate list of reduced distance between reference atom ia and all others
      list_reddist(:,ib) = hist%xred(:,ib,1) - xred_ref(:,ia)
      !If the distorted atom is further away than half the unit cell shift it.
      if(list_reddist(1,ib) > 0.5)then
         list_reddist(1,ib) = 1 -  list_reddist(1,ib)
         shift(1,ib) = -1
      end if
      if(list_reddist(2,ib) > 0.5)then
         list_reddist(2,ib) = 1 -  list_reddist(2,ib)
         shift(2,ib) = -1
      end if
      if(list_reddist(3,ib) > 0.5)then
         list_reddist(3,ib) = 1 -  list_reddist(3,ib)
         shift(3,ib) = -1
      end if
      if(list_reddist(1,ib) < -0.5)then
         list_reddist(1,ib) = -1 -  list_reddist(1,ib)
         shift(1,ib) = 1
      end if
      if(list_reddist(2,ib) < -0.5)then
         list_reddist(2,ib) = -1 -  list_reddist(2,ib)
         shift(2,ib) = 1
      end if
      if(list_reddist(3,ib) < -0.5)then
         list_reddist(3,ib) = -1 -  list_reddist(3,ib)
         shift(3,ib) = 1
      end if
      list_absdist(1,ib) = (rprimd_hist(1,1)+rprimd_hist(2,1)+rprimd_hist(3,1))*list_reddist(1,ib)
      list_absdist(2,ib) = (rprimd_hist(1,2)+rprimd_hist(2,2)+rprimd_hist(3,2))*list_reddist(2,ib)
      list_absdist(3,ib) = (rprimd_hist(1,3)+rprimd_hist(2,3)+rprimd_hist(3,3))*list_reddist(3,ib)
      list_dist(ib) = sqrt(abs(list_absdist(1,ib))**2 + abs(list_absdist(2,ib))**2 + abs(list_absdist(3,ib))**2 )
   end do !ib
   !find the closest atom ib
   list_map(ia) = minloc(list_dist,DIM=1)
   !If the closest atom ib was shifted, apply and store the shift
   if(any(shift(:,list_map(ia)) /= 0))then
      hist%xred(1,list_map(ia),:)= hist%xred(1,list_map(ia),:) + 1*shift(1,list_map(ia))
      hist%xred(2,list_map(ia),:)= hist%xred(2,list_map(ia),:) + 1*shift(2,list_map(ia))
      hist%xred(3,list_map(ia),:)= hist%xred(3,list_map(ia),:) + 1*shift(3,list_map(ia))
   end if
   !TEST MS
   !write(*,*) 'Atom', ia,' of reference is matche with', list_map(ia)
   !write(*,*) 'xred_ref(',xred_ref(:,ia),'), xred_hist(',hist%(:,list_map(ia)),')'
end do  ! ia

 if(need_verbose) then
   write(msg,'(2a,I3,a,I3,a,I3)') ch10,&
&       ' The size of the supercell for the fit is ',ncell(1),' ',ncell(2),' ',ncell(3)
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 end if


   if(any(list_map(:)==0))then
       write(msg, '(5a)' )&
&         'Unable to map the molecular dynamic file  ',ch10,&
&         'on the reference supercell structure',ch10,&
&         'Action: change the MD file'
       MSG_ERROR(msg)
   end if

 need_map = .FALSE.
 do ia=1,natom_hist
   if(list_map(ia) /= ia) need_map = .TRUE.
 end do
 if(need_map)then
   if(need_verbose) then
     write(msg, '(11a)' )ch10,&
&      ' --- !WARNING',ch10,&
&      '     The ordering of the atoms in the _HIST.nc file is different,',ch10,&
&      '     of the one built by multibinit. The _HIST.nc file will be mapped,',ch10,&
&      '     to the ordering of multibinit.',ch10,&
&      ' ---',ch10
     call wrtout(ab_out,msg,'COLL')
     call wrtout(std_out,msg,'COLL')
   end if

! Allocate hist datatype
   call abihist_init(hist_tmp,natom_hist,nstep_hist,.false.,.false.)
! copy all the information
   do ia=1,nstep_hist
     hist%ihist = ia
     hist_tmp%ihist = ia
     call abihist_copy(hist,hist_tmp)
   end do
   hist_tmp%mxhist = nstep_hist

! reoder array
   do ia=1,natom_hist
     hist_tmp%xred(:,ia,:)  = hist%xred(: ,list_map(ia),:)
     hist_tmp%fcart(:,ia,:) = hist%fcart(:,list_map(ia),:)
     hist_tmp%vel(:,ia,:)   = hist%vel(:,list_map(ia),:)
   end do

! free the old hist and reinit
   call abihist_free(hist)
   call abihist_init(hist,natom_hist,nstep_hist,.false.,.false.)
! copy the temporary hist into output
   do ia=1,nstep_hist
     hist%ihist = ia
     hist_tmp%ihist = ia
     call abihist_copy(hist_tmp,hist)
   end do
   hist_tmp%mxhist = nstep_hist
   call abihist_free(hist_tmp)
 end if !need map

!deallocation
 ABI_DEALLOCATE(shift)
 ABI_DEALLOCATE(list_map)
 ABI_DEALLOCATE(list_dist)
 ABI_DEALLOCATE(list_reddist)
 ABI_DEALLOCATE(list_absdist)
 ABI_DEALLOCATE(xred_ref)
end subroutine effective_potential_file_mapHistToRef
!!***


!****f* m_effective_potential_file/effective_potential_file_readDisplacement
!!
!! NAME
!! effective_potential_file_readDisplacement
!!
!! FUNCTION
!! Read a displacement ASCII file
!!
!! INPUTS
!! filename = path of the file
!! natom = number of atoms in the cell
!! nstep = number of time step
!!
!! OUTPUT
!! disp(3,natom_sc) = atomics displacement between configuration and the reference
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine effective_potential_file_readDisplacement(filename,disp,nstep,natom)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom,nstep
 character(len=fnlen),intent(in) :: filename
!array
 real(dp),intent(out) :: disp(nstep,3,natom)
!Local variables------------------------------
!scalars
 integer :: ios = 0
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
     read(funit,'(a)',iostat=ios) readline
     line=adjustl(readline)
     read(unit=line,fmt=*) (disp(istep,mu,ia),mu=1,3)
     write(111,'(3(es21.12))') disp(istep,:,ia)
   end do
 end do

 close(funit)

end subroutine effective_potential_file_readDisplacement
!!***

!!****f* m_effective_potential_file/elementfromline
!! NAME
!! elementfromline
!!
!! FUNCTION
!! Read the number of element of a line
!!
!! INPUTS
!!  line= string from which the data are read
!!
!! OUTPUT
!!  nelement = number of element in the line
!!
!! PARENTS
!!      m_effective_potential_file
!!
!! CHILDREN
!!
!! SOURCE

subroutine elementfromline(line,nelement)

!Arguments ---------------------------------------------
 character(len=*), intent(in) :: line
 integer, intent(out) :: nelement
!Local variables ---------------------------------------
 integer :: ii,n
 logical :: element
! *********************************************************************

!Set the output
 nelement = 0
 n = len_trim(line)
 element = .false.
 do ii=1,n
   if(.not.element.and.line(ii:ii) /="")  then
     element=.true.
   else
     if((element.and.line(ii:ii) =="")) then
       element=.false.
       nelement = nelement + 1
     end if
   end if
   if((element.and.ii==n)) nelement = nelement + 1
 end do

 end subroutine elementfromline
!!***

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
!!      m_effective_potential_file
!!
!! CHILDREN
!!
!! SOURCE

 subroutine rdfromline(keyword,line,output)

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
!!    system_xml2effpot
!!
!! CHILDREN
!!    rmtabfromline
!!
!! SOURCE

recursive subroutine rmtabfromline(line)

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
!!      system_xml2effpot
!!
!! CHILDREN
!!      paw_rdfromline
!!
!! SOURCE

 subroutine rdfromline_value(keyword,line,output)

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

#if defined HAVE_XML

function char_f2c(f_string) result(c_string)

 use iso_c_binding, only : C_CHAR,C_NULL_CHAR
!Arguments ------------------------------------
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

 use iso_c_binding, only : C_CHAR,C_NULL_CHAR
!Arguments ------------------------------------
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
#endif

end module m_effective_potential_file
!!***
