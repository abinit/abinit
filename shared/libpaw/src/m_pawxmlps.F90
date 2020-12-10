!!****m* ABINIT/m_pawxmlps
!! NAME
!! m_pawxmlps
!!
!! FUNCTION
!! This module reads a PAW pseudopotential file written in XML.
!! Can use either FoX or pure Fortran routines.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2020 ABINIT group (MT, FJ)
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

module m_pawxmlps

 USE_DEFS
 USE_MSG_HANDLING
 USE_MEMORY_PROFILING

#if defined LIBPAW_HAVE_FOX
 use fox_sax
#endif

 use m_pawrad     , only : pawrad_type, pawrad_init, pawrad_free, pawrad_ifromr, bound_deriv
 use m_paw_numeric, only : paw_spline, paw_splint

 implicit none

 private

!Procedures used for the Fortran reader
 public  ::  rdpawpsxml
 public  ::  rdpawpsxml_header
 public  ::  rdpawpsxml_core

!Procedures used for the FoX reader (called from xml_parser in response to events)
#if defined LIBPAW_HAVE_FOX
 public :: paw_begin_element1
 public :: paw_end_element1
 public :: pawdata_chunk
#endif

! private procedures and global variables.
 private ::  paw_rdfromline

! This type of real is used by the datatypes below
 integer,parameter,private :: dpxml = selected_real_kind(14)

! The maximum length of a record in a file connected for sequential access.
 integer,parameter,private :: XML_RECL = 50000
 integer,parameter,private :: NGAUSSIAN_MAX = 100
!!***

!!****t* m_pawxmlps/radial_grid_t
!! NAME
!! radial_grid_t
!!
!! FUNCTION
!! Radial grid type for the FoX XML reader
!!
!! SOURCE

type, public :: radial_grid_t
  logical           :: tread=.false.
  character(len=20) :: eq
  character(len=4)  :: id
  real(dpxml)       :: aa
  real(dpxml)       :: bb
  real(dpxml)       :: dd
  integer           :: nn
  integer           :: istart
  integer           :: iend
end type radial_grid_t
!!***

!-------------------------------------------------------------------------

!!****t* m_pawxmlps/radial_func_t
!! NAME
!! radial_func_t
!!
!! FUNCTION
!! Radial function type for the FoX XML reader
!!
!! SOURCE

type, public              :: radialfunc_t
  logical             :: tread=.false.
  character(len=6)    :: grid=' '  !vz_z
  character(len=6)    :: state=' ' !vz_z
  real(dpxml),allocatable :: data(:)
end type radialfunc_t
!!***

!-------------------------------------------------------------------------

!!****t* m_pawxmlps/gaussian_expansion_t
!! NAME
!! radial_func_t
!!
!! FUNCTION
!! Radial function type for the FoX XML reader
!!
!! SOURCE

type, public              :: gaussian_expansion_t
  logical             :: tread=.false.
  integer             :: ngauss=0
  character(len=6)    :: state=' '
  real(dpxml), dimension(2, NGAUSSIAN_MAX) :: factors
  real(dpxml), dimension(2, NGAUSSIAN_MAX) :: expos
end type gaussian_expansion_t
!!***

!-------------------------------------------------------------------------

!!****t* m_pawxmlps/shape_function_t
!! NAME
!! shape_function_t
!!
!! FUNCTION
!! Shape function type for the FoX XML reader
!!
!! SOURCE

type, public :: shape_function_t
  logical              :: tread=.false.
  character(len=20)    :: gtype
  real(dpxml)          :: rc=0
  character(len=6)     :: grid
  integer              :: lamb
  real(dpxml), allocatable :: data(:,:)
end type shape_function_t
!!***

!-------------------------------------------------------------------------

!!****t* m_pawxmlps/state_t
!! NAME
!! state_t
!!
!! FUNCTION
!! State type for the FoX XML reader
!!
!! SOURCE

type, public :: state_t
  logical          :: tread=.false.
  character(len=6) :: id
  real(dpxml)      :: ff
  real(dpxml)      :: rc
  real(dpxml)      :: ee
  integer          :: nn
  integer          :: ll
end type state_t
!!***

!-------------------------------------------------------------------------

!!****t* m_pawxmlps/valence_states_t
!! NAME
!! valence_states_t
!!
!! FUNCTION
!! Valence state type for the FoX XML reader
!!
!! SOURCE

type, public :: valence_states_t
  logical               :: tread=.false.
  integer               :: nval
  type(state_t),allocatable :: state(:)
end type valence_states_t
!!***

!-------------------------------------------------------------------------

!!****t* m_pawxmlps/generator_t
!! NAME
!! generator_t
!!
!! FUNCTION
!! Generator type for the FoX XML reader
!!
!! SOURCE

type, public :: generator_t
  logical           :: tread=.false.
  character(len=20) :: gen
  character(len=20) :: name
end type generator_t
!!***

!-------------------------------------------------------------------------

!!****t* m_pawxmlps/xc_functional_t
!! NAME
!! xc_functional_t function_t
!!
!! FUNCTION
!! XC functional type for the FoX XML reader
!!
!! SOURCE

type, public :: xc_functional_t
  logical           :: tread=.false.
  character(len=12) :: functionaltype
  character(len=100) :: name
end type xc_functional_t
!!***

!-------------------------------------------------------------------------

!!****t* m_pawxmlps/atom_t
!! NAME
!! atom_t
!!
!! FUNCTION
!! Atom type for the FoX XML reader
!!
!! SOURCE

type, public :: atom_t
  logical           :: tread=.false.
  character(len=2)  :: symbol
  real(dpxml)       :: znucl
  real(dpxml)       :: zion
  real(dpxml)       :: zval
end type atom_t
!!***

!-------------------------------------------------------------------------

!!****t* m_pawxmlps/paw_setup_t
!! NAME
!! paw_setup_t
!!
!! FUNCTION
!! PAW setup type (contain all the data for a PAW setup)
!!
!! SOURCE

type, public :: paw_setup_t
  character(len=3)             :: version
  logical                      :: tread=.false.
  integer                      :: ngrid
  real(dpxml)                  :: rpaw
  real(dpxml)                  :: ex_cc
  character(len=4)             :: idgrid
  character(len=12)            :: optortho
  type(atom_t)                 :: atom
  type(xc_functional_t)        :: xc_functional
  type(generator_t)            :: generator
  type(valence_states_t)       :: valence_states
  type(radial_grid_t), allocatable :: radial_grid(:)
  type(shape_function_t)       :: shape_function
  type(radialfunc_t)           :: ae_core_density
  type(radialfunc_t)           :: pseudo_core_density
  type(radialfunc_t)           :: pseudo_valence_density
  type(radialfunc_t)           :: zero_potential
  type(radialfunc_t)           :: LDA_minus_half_potential
  type(radialfunc_t)           :: ae_core_kinetic_energy_density
  type(radialfunc_t)           :: pseudo_core_kinetic_energy_density
  type(radialfunc_t),allocatable :: ae_partial_wave(:)
  type(radialfunc_t),allocatable :: pseudo_partial_wave(:)
  type(radialfunc_t),allocatable :: projector_function(:)
  type(gaussian_expansion_t),allocatable :: projector_fit(:)
  type(radialfunc_t)           :: kresse_joubert_local_ionic_potential
  type(radialfunc_t)           :: blochl_local_ionic_potential
  type(radialfunc_t)           :: kinetic_energy_differences
  type(radialfunc_t)           :: exact_exchange_matrix
end type paw_setup_t

 public :: paw_setup_free  ! Free memory
 public :: paw_setup_copy     ! Copy object
!!***


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!------------- PUBLIC AND PRIVATE VARIABLES ------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!Public variables (common to both readers)
integer,save,public,allocatable :: ipsp2xml(:)
integer,save,public :: npsp_pawxml
type(paw_setup_t),public,target,allocatable,save :: paw_setup(:)
type(paw_setup_t),public,target,save :: paw_setuploc

!Private variables (for the FoX reader)
#if defined LIBPAW_HAVE_FOX
logical,private,save  :: in_valenceStates = .false.,in_data=.false.
logical,private,save  :: in_generator =.false.
integer,private,save :: ndata
integer,private,save  :: ii,ival,igrid,ishpf,lmax,mesh_size
!Pointers to make it easier to manage the data
type(radialfunc_t),private,save,pointer  :: rp
type(state_t),private,save,pointer   :: valstate (:)
type(radial_grid_t),private,save,pointer   :: grids (:)
type(radialfunc_t),private,save,pointer :: shpf(:)
#endif
!!***


CONTAINS
!===========================================================
!!***

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!------------- ROUTINES AND FUNCTIONS FOR THE FOX READER -----------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
#if defined LIBPAW_HAVE_FOX

!!****f* m_pawxmlps/paw_begin_element1
!! NAME
!! begin_element
!!
!! FUNCTION
!!  Read an XML tag with a given name.
!!  Fills the present module private data.
!!
!! INPUTS
!!  namespaceURI = universal resource indicator for XML namespace??? Not used.
!!  localName = local equivalent of tag name?? Not used.
!!  name = name of XML tag which has been read in
!!  attributes = attributes of XML tag
!!
!! OUTPUT
!!  Fills private data in present module.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine paw_begin_element1(namespaceURI,localName,name,attributes)

character(len=*),intent(in)   :: namespaceURI,localName,name
type(dictionary_t),intent(in) :: attributes

character(len=100)  :: msg,value
integer ::iaewf=0,iproj=0,ipswf=0,iprojfit=0,igauss=0
!Just to fool abirules
 value=localName
 value=namespaceURI

select case(name)

      case ("paw_setup")
        paw_setuploc%tread=.true.
        igrid=0;ishpf=0
        paw_setuploc%rpaw=-1.d0
        LIBPAW_DATATYPE_ALLOCATE(grids,(10))
        LIBPAW_DATATYPE_ALLOCATE(shpf,(7))
        value = getValue(attributes,"version")
        write(std_out,'(3a)') "Processing a PSEUDO version ",trim(value)," XML file"
        paw_setuploc%version=trim(value)


      case ("atom")
         paw_setuploc%atom%tread=.true.
         value = getValue(attributes,"symbol")
         if (value == "" ) then
           msg="Cannot determine atomic symbol"
           LIBPAW_ERROR(msg)
         end if
         paw_setuploc%atom%symbol = trim(value)

         value = getValue(attributes,"Z") 
         if (value == "" ) then
           msg="Cannot determine znucl"
           LIBPAW_ERROR(msg)
         end if
         read(unit=value,fmt=*) paw_setuploc%atom%znucl

         value = getValue(attributes,"core")
         if (value == "" ) then
           msg="Cannot determine zion"
           LIBPAW_ERROR(msg)
         end if
         read(unit=value,fmt=*) paw_setuploc%atom%zion

         value = getValue(attributes,"valence")
         if (value == "" ) then
           msg="Cannot determine zval"
           LIBPAW_ERROR(msg)
         end if
         read(unit=value,fmt=*) paw_setuploc%atom%zval


      case ("xc_functional")
         paw_setuploc%xc_functional%tread=.true.
         value = getValue(attributes,"type")
         if (value == "" ) then
           msg="Cannot determine xc-functional-type"
           LIBPAW_ERROR(msg)
         end if
         paw_setuploc%xc_functional%functionaltype = trim(value)

         value = getValue(attributes,"name")
         if (value == "" ) then
           msg="Cannot determine xc-functional-name "
           LIBPAW_ERROR(msg)
         end if
         paw_setuploc%xc_functional%name= trim(value)

      case ("generator")
         paw_setuploc%generator%tread=.true.
         in_generator =.true.
         value = getValue(attributes,"type")
         if (value == "" ) value = "unknown"
         paw_setuploc%generator%gen = trim(value)

         value = getValue(attributes,"name")
         if (value == "" ) value = "unknown"
         paw_setuploc%generator%name = trim(value)

      case ("PAW_radius")
         value = getValue(attributes,"rpaw")
         if (value == "" ) then
           msg="Cannot determine rpaw"
           LIBPAW_ERROR(msg)
         end if
         read(unit=value,fmt=*) paw_setuploc%rpaw

      case ("valence_states")
         paw_setuploc%valence_states%tread=.true.
         in_valenceStates=.true.
         ival=0
         lmax=0
         LIBPAW_DATATYPE_ALLOCATE(valstate,(50))

      case ("state")
         ival=ival+1

         value = getValue(attributes,"n")
         if (value == "" ) then 
           valstate(ival)%nn=-1    
         else
           read(unit=value,fmt=*) valstate(ival)%nn
         end if
 
         value = getValue(attributes,"l")
         if (value == "" ) then
           msg="Cannot determine l"
           LIBPAW_ERROR(msg)
         end if
         read(unit=value,fmt=*) valstate(ival)%ll
         if(valstate(ival)%ll>lmax) lmax=valstate(ival)%ll

         value = getValue(attributes,"f")
         if (value == "" ) then 
           valstate(ival)%ff=-1.d0
         else
           read(unit=value,fmt=*) valstate(ival)%ff
         end if

         value = getValue(attributes,"rc")
         if (value == "" ) then
           msg="Cannot determine rc"
           LIBPAW_ERROR(msg)
         end if
         read(unit=value,fmt=*) valstate(ival)%rc

         value = getValue(attributes,"e")
         if (value == "" ) then
           msg="Cannot determine e"
           LIBPAW_ERROR(msg)
         end if
         read(unit=value,fmt=*) valstate(ival)%ee

         value = getValue(attributes,"id")
         if (value == "" ) value = "unknown"
         valstate(ival)%id = trim(value)

      case ("radial_grid")
         igrid=igrid+1
         value = getValue(attributes,"eq")
         if (value == "" ) value = "unknown"
         grids(igrid)%eq = trim(value)

         value = getValue(attributes,"a")
         if (value == "" ) then
           grids(igrid)%aa=0.d0
         else
           read(unit=value,fmt=*) grids(igrid)%aa
         end if

         value = getValue(attributes,"n")
         if (value == "" ) then
           grids(igrid)%nn=0
         else
           read(unit=value,fmt=*) grids(igrid)%nn
         end if

         value = getValue(attributes,"d")
         if (value == "" ) then
           grids(igrid)%dd=0.d0
         else
           read(unit=value,fmt=*) grids(igrid)%dd
         end if

         value = getValue(attributes,"b")
         if (value == "" ) then
           grids(igrid)%bb=0.d0
         else
           read(unit=value,fmt=*) grids(igrid)%bb
         end if

         value = getValue(attributes,"istart")
         if (value == "" ) then
           msg="Cannot determine istart"
           LIBPAW_ERROR(msg)
         end if
         read(unit=value,fmt=*) grids(igrid)%istart

         value = getValue(attributes,"iend")
         if (value == "" ) then
           msg="Cannot determine iend"
           LIBPAW_ERROR(msg)
         end if
         read(unit=value,fmt=*) grids(igrid)%iend

         value = getValue(attributes,"id")
         if (value == "" ) value = "unknown"
         grids(igrid)%id = trim(value)

end select

select case(name)
      case ("shape_function")
         paw_setuploc%shape_function%tread=.true.
         value = getValue(attributes,"type")
         if (value == "" ) value = "unknown"
         paw_setuploc%shape_function%gtype = trim(value)

         value = getValue(attributes,"grid")
         paw_setuploc%shape_function%grid=trim(value)
         if (value /= "" ) then
           paw_setuploc%shape_function%gtype ="num" 
           do ii=1,igrid
             if(trim(paw_setuploc%shape_function%grid)==trim(grids(ii)%id)) then
               mesh_size=grids(ii)%iend-grids(ii)%istart+1
             end if
           end do
           ishpf=ishpf+1
           LIBPAW_ALLOCATE(shpf(ishpf)%data,(mesh_size))
           rp=>shpf(ishpf)
           in_data=.true.
           ndata = 0
         end if

         value = getValue(attributes,"rc")
         if (value == "" ) then
           if(paw_setuploc%shape_function%gtype /="num") then
              msg="Cannot determine rc"
              LIBPAW_ERROR(msg)
           end if
         else
           read(unit=value,fmt=*) paw_setuploc%shape_function%rc
         end if

         value = getValue(attributes,"lamb")
         if (value == "" ) then
           paw_setuploc%shape_function%lamb=0
         else
           read(unit=value,fmt=*) paw_setuploc%shape_function%lamb
         end if

      case ("pseudo_partial_wave")
         ipswf=ipswf+1
         paw_setuploc%pseudo_partial_wave(ipswf)%tread=.true.
         value = getValue(attributes,"grid")
         if (value == "" ) value = "unknown"
         paw_setuploc%idgrid = trim(value)
         paw_setuploc%pseudo_partial_wave(ipswf)%grid=trim(value)

         value = getValue(attributes,"state")
         if (value == "" ) then
           msg="Cannot determine pseudo_partial_wave state"
           LIBPAW_ERROR(msg)
         end if
         paw_setuploc%pseudo_partial_wave(ipswf)%state=trim(value)

         do ii=1,igrid
           if(trim(paw_setuploc%pseudo_partial_wave(ipswf)%grid)==trim(grids(ii)%id)) then
             mesh_size=grids(ii)%iend-grids(ii)%istart+1
           end if
         end do

         LIBPAW_ALLOCATE(paw_setuploc%pseudo_partial_wave(ipswf)%data,(mesh_size))
         rp=>paw_setuploc%pseudo_partial_wave(ipswf)
         if(ipswf==paw_setuploc%valence_states%nval) ipswf=0
         in_data=.true.
         ndata = 0

      case ("ae_partial_wave")
         iaewf=iaewf+1
         paw_setuploc%ae_partial_wave(iaewf)%tread=.true.
         value = getValue(attributes,"grid")
         if (value == "" ) value = "unknown"
         paw_setuploc%ae_partial_wave(iaewf)%grid=trim(value)

         value = getValue(attributes,"state")
         if (value == "" ) then
           LIBPAW_ERROR("Cannot determine ae_partial_wave state")
         end if
         paw_setuploc%ae_partial_wave(iaewf)%state=trim(value)

         do ii=1,igrid
           if(trim(paw_setuploc%ae_partial_wave(iaewf)%grid)==trim(grids(ii)%id)) then
             mesh_size=grids(ii)%iend-grids(ii)%istart+1
           end if
         end do

         LIBPAW_ALLOCATE(paw_setuploc%ae_partial_wave(iaewf)%data,(mesh_size))
         rp=>paw_setuploc%ae_partial_wave(iaewf)
         if(iaewf==paw_setuploc%valence_states%nval) iaewf=0
         in_data=.true.
         ndata = 0

      case ("projector_function")
         iproj=iproj+1
         paw_setuploc%projector_function(iproj)%tread=.true.
         value = getValue(attributes,"grid")
         if (value == "" ) value = "unknown"
         paw_setuploc%projector_function(iproj)%grid=trim(value)

         value = getValue(attributes,"state")
         if (value == "" ) then
           msg="Cannot determine projector_function state"
           LIBPAW_ERROR(msg)
         end if
         paw_setuploc%projector_function(iproj)%state=trim(value)

         do ii=1,igrid
           if(trim(paw_setuploc%projector_function(iproj)%grid)==trim(grids(ii)%id)) then
             mesh_size=grids(ii)%iend-grids(ii)%istart+1
           end if
         end do

         LIBPAW_ALLOCATE(paw_setuploc%projector_function(iproj)%data,(mesh_size))
         rp=>paw_setuploc%projector_function(iproj)
         if(iproj==paw_setuploc%valence_states%nval) iproj=0
         in_data=.true.
         ndata = 0

      case ("projector_fit")
         if(.not.allocated(paw_setuploc%projector_fit)) then
            LIBPAW_DATATYPE_ALLOCATE(paw_setuploc%projector_fit,(paw_setuploc%valence_states%nval))
         end if

         iprojfit=iprojfit+1
         paw_setuploc%projector_fit(iprojfit)%tread=.true.
         value = getValue(attributes,"state")
         if (value == "" ) then
           msg="Cannot determine projector_fit state"
           LIBPAW_ERROR(msg)
         end if
         paw_setuploc%projector_fit(iprojfit)%state=trim(value)

         if(iprojfit==paw_setuploc%valence_states%nval) iprojfit=0
         igauss = 0

      case ("gaussian")
         igauss = igauss + 1
         value = getValue(attributes,"factor")
         if (value == "" ) then
           msg="Cannot determine gaussian factor"
           LIBPAW_ERROR(msg)
         end if
         read(value(2:100), *) paw_setuploc%projector_fit(iprojfit)%factors(1, igauss)
         read(value(index(value, ',') + 1:100), *) paw_setuploc%projector_fit(iprojfit)%factors(2, igauss)
         value = getValue(attributes,"exponent")
         if (value == "" ) then
           msg="Cannot determine gaussian exponent"
           LIBPAW_ERROR(msg)
         end if
         read(value(2:100), *) paw_setuploc%projector_fit(iprojfit)%expos(1, igauss)
         read(value(index(value, ',') + 1:100), *) paw_setuploc%projector_fit(iprojfit)%expos(2, igauss)

     case ("ae_core_density")
         paw_setuploc%ae_core_density%tread=.true.
          value = getValue(attributes,"grid")
         if (value == "" ) value = "unknown"
         paw_setuploc%ae_core_density%grid=trim(value)

         value = getValue(attributes,"state")
         if (value == "" ) value = "unknown"
         paw_setuploc%ae_core_density%state=trim(value)

         do ii=1,igrid
           if(trim(paw_setuploc%ae_core_density%grid)==trim(grids(ii)%id)) then
             mesh_size=grids(ii)%iend-grids(ii)%istart+1
           end if
         end do

         LIBPAW_ALLOCATE(paw_setuploc%ae_core_density%data,(mesh_size))
         rp=>paw_setuploc%ae_core_density
         in_data=.true.
         ndata = 0

     case ("pseudo_core_density")
         paw_setuploc%pseudo_core_density%tread=.true.
          value = getValue(attributes,"grid")
         if (value == "" ) value = "unknown"
         paw_setuploc%pseudo_core_density%grid=trim(value)

         value = getValue(attributes,"state")
         if (value == "" ) value = "unknown"
         paw_setuploc%pseudo_core_density%state=trim(value)

         do ii=1,igrid
           if(trim(paw_setuploc%pseudo_core_density%grid)==trim(grids(ii)%id)) then
             mesh_size=grids(ii)%iend-grids(ii)%istart+1
           end if
         end do

         LIBPAW_ALLOCATE(paw_setuploc%pseudo_core_density%data,(mesh_size))
         rp=>paw_setuploc%pseudo_core_density
         in_data=.true.
         ndata = 0

     case ("pseudo_valence_density")
         paw_setuploc%pseudo_valence_density%tread=.true.
          value = getValue(attributes,"grid")
         if (value == "" ) value = "unknown"
         paw_setuploc%pseudo_valence_density%grid=trim(value)

         value = getValue(attributes,"state")
         if (value == "" ) value = "unknown"
         paw_setuploc%pseudo_valence_density%state=trim(value)

         do ii=1,igrid
           if(trim(paw_setuploc%pseudo_valence_density%grid)==trim(grids(ii)%id)) then
             mesh_size=grids(ii)%iend-grids(ii)%istart+1
           end if
         end do

         LIBPAW_ALLOCATE(paw_setuploc%pseudo_valence_density%data,(mesh_size))
         rp=>paw_setuploc%pseudo_valence_density
         in_data=.true.
         ndata = 0

     case ("zero_potential")
         paw_setuploc%zero_potential%tread=.true.
          value = getValue(attributes,"grid")
         if (value == "" ) value = "unknown"
         paw_setuploc%zero_potential%grid=trim(value)

         value = getValue(attributes,"state")
         if (value == "" ) value = "unknown"
         paw_setuploc%zero_potential%state=trim(value)

         do ii=1,igrid
           if(trim(paw_setuploc%zero_potential%grid)==trim(grids(ii)%id)) then
             mesh_size=grids(ii)%iend-grids(ii)%istart+1
           end if
         end do

         LIBPAW_ALLOCATE(paw_setuploc%zero_potential%data,(mesh_size))
         rp=>paw_setuploc%zero_potential
         in_data=.true.
         ndata = 0

     case ("LDA_minus_half_potential")
         paw_setuploc%LDA_minus_half_potential%tread=.true.
          value = getValue(attributes,"grid")
         if (value == "" ) value = "unknown"
         paw_setuploc%LDA_minus_half_potential%grid=trim(value)

         value = getValue(attributes,"state")
         if (value == "" ) value = "unknown"
         paw_setuploc%LDA_minus_half_potential%state=trim(value)

         do ii=1,igrid
           if(trim(paw_setuploc%LDA_minus_half_potential%grid)==trim(grids(ii)%id)) then
             mesh_size=grids(ii)%iend-grids(ii)%istart+1
           end if
         end do

         LIBPAW_ALLOCATE(paw_setuploc%LDA_minus_half_potential%data,(mesh_size))
         rp=>paw_setuploc%LDA_minus_half_potential
         in_data=.true.
         ndata = 0

     case ("ae_core_kinetic_energy_density")
         paw_setuploc%ae_core_kinetic_energy_density%tread=.true.
          value = getValue(attributes,"grid")
         if (value == "" ) value = "unknown"
         paw_setuploc%ae_core_kinetic_energy_density%grid=trim(value)

         value = getValue(attributes,"state")
         if (value == "" ) value = "unknown"
         paw_setuploc%ae_core_kinetic_energy_density%state=trim(value)

         do ii=1,igrid
           if(trim(paw_setuploc%ae_core_kinetic_energy_density%grid)==trim(grids(ii)%id)) then
             mesh_size=grids(ii)%iend-grids(ii)%istart+1
           end if
         end do

         LIBPAW_ALLOCATE(paw_setuploc%ae_core_kinetic_energy_density%data,(mesh_size))
         rp=>paw_setuploc%ae_core_kinetic_energy_density
         in_data=.true.
         ndata = 0

     case ("pseudo_core_kinetic_energy_density")
         paw_setuploc%pseudo_core_kinetic_energy_density%tread=.true.
          value = getValue(attributes,"grid")
         if (value == "" ) value = "unknown"
         paw_setuploc%pseudo_core_kinetic_energy_density%grid=trim(value)

         value = getValue(attributes,"state")
         if (value == "" ) value = "unknown"
         paw_setuploc%pseudo_core_kinetic_energy_density%state=trim(value)

         do ii=1,igrid
           if(trim(paw_setuploc%pseudo_core_kinetic_energy_density%grid)==trim(grids(ii)%id)) then
             mesh_size=grids(ii)%iend-grids(ii)%istart+1
           end if
         end do

         LIBPAW_ALLOCATE(paw_setuploc%pseudo_core_kinetic_energy_density%data,(mesh_size))
         rp=>paw_setuploc%pseudo_core_kinetic_energy_density
         in_data=.true.
         ndata = 0

     case ("kresse_joubert_local_ionic_potential")
         paw_setuploc%kresse_joubert_local_ionic_potential%tread=.true.
          value = getValue(attributes,"grid")
         if (value == "" ) value = "unknown"
         paw_setuploc%kresse_joubert_local_ionic_potential%grid=trim(value)

         value = getValue(attributes,"state")
         if (value == "" ) value = "unknown"
         paw_setuploc%kresse_joubert_local_ionic_potential%state=trim(value)

         do ii=1,igrid
           if(trim(paw_setuploc%kresse_joubert_local_ionic_potential%grid)==trim(grids(ii)%id)) then
             mesh_size=grids(ii)%iend-grids(ii)%istart+1
           end if
         end do

         LIBPAW_ALLOCATE(paw_setuploc%kresse_joubert_local_ionic_potential%data,(mesh_size))
         rp=>paw_setuploc%kresse_joubert_local_ionic_potential
         in_data=.true.
         ndata = 0

     case ("blochl_local_ionic_potential")
         paw_setuploc%blochl_local_ionic_potential%tread=.true.
          value = getValue(attributes,"grid")
         if (value == "" ) value = "unknown"
         paw_setuploc%blochl_local_ionic_potential%grid=trim(value)

         value = getValue(attributes,"state")
         if (value == "" ) value = "unknown"
         paw_setuploc%blochl_local_ionic_potential%state=trim(value)

         do ii=1,igrid
           if(trim(paw_setuploc%blochl_local_ionic_potential%grid)==trim(grids(ii)%id)) then
             mesh_size=grids(ii)%iend-grids(ii)%istart+1
           end if
         end do

         LIBPAW_ALLOCATE(paw_setuploc%blochl_local_ionic_potential%data,(mesh_size))
         rp=>paw_setuploc%blochl_local_ionic_potential
         in_data=.true.
         ndata = 0

    case ("kinetic_energy_differences")
         paw_setuploc%kinetic_energy_differences%tread=.true.
         mesh_size=paw_setuploc%valence_states%nval*paw_setuploc%valence_states%nval
         LIBPAW_ALLOCATE(paw_setuploc%kinetic_energy_differences%data,(mesh_size))
         rp=>paw_setuploc%kinetic_energy_differences
         in_data=.true.
         ndata = 0

end select

end subroutine paw_begin_element1
!!***

!-------------------------------------------------------------------------

!!****f* m_pawxmlps/paw_end_element1
!! NAME
!! end_element
!!
!! FUNCTION
!!  End XML tag effect: switches flags in private data of this module
!!
!! INPUTS
!!  namespaceURI = universal resource indicator for XML namespace??? Not used.
!!  localName = local equivalent of tag name?? Not used.
!!  name = name of XML tag which has been read in
!!
!! OUTPUT
!!  side effect: private data flags in present module are turned to .false.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine paw_end_element1(namespaceURI,localName,name)

character(len=*),intent(in) :: namespaceURI,localName,name
character(len=100) :: msg,value

!Just to fool abirules 
 value=localName
 value=namespaceURI
 
select case(name)

      case ("generator")
         in_generator = .false.

      case ("valence_states")
        in_valenceStates = .false.
        if(ival>50) then
          msg="ival>50"
          LIBPAW_ERROR(msg)
        end if
        if(ival>0)then
          LIBPAW_DATATYPE_ALLOCATE(paw_setuploc%valence_states%state,(ival))
          paw_setuploc%valence_states%state(ival)%tread=.true.
          paw_setuploc%valence_states%nval=ival
          do ii=1,ival
            paw_setuploc%valence_states%state(ii)=valstate(ii)
          end do
        end if
        LIBPAW_DATATYPE_DEALLOCATE(valstate)
        if(.not.allocated(paw_setuploc%ae_partial_wave)) then
          LIBPAW_DATATYPE_ALLOCATE(paw_setuploc%ae_partial_wave,(paw_setuploc%valence_states%nval))
        end if
        if(.not.allocated(paw_setuploc%pseudo_partial_wave)) then
          LIBPAW_DATATYPE_ALLOCATE(paw_setuploc%pseudo_partial_wave,(paw_setuploc%valence_states%nval))
        end if
        if(.not.allocated(paw_setuploc%projector_function)) then
          LIBPAW_DATATYPE_ALLOCATE(paw_setuploc%projector_function,(paw_setuploc%valence_states%nval))
        end if

      case ("paw_setup")
        if(igrid>10) then
          msg="igrid>10"
          LIBPAW_ERROR(msg)
        end if
        LIBPAW_DATATYPE_ALLOCATE(paw_setuploc%radial_grid,(igrid))
        paw_setuploc%radial_grid(igrid)%tread=.true.
        paw_setuploc%ngrid=igrid
        do ii=1,igrid
          paw_setuploc%radial_grid(ii)=grids(ii)
        end do
        LIBPAW_DATATYPE_DEALLOCATE(grids)
        do ii=1,igrid
          if(trim(paw_setuploc%shape_function%grid)==trim(paw_setuploc%radial_grid(ii)%id)) then
            mesh_size=paw_setuploc%radial_grid(ii)%iend-paw_setuploc%radial_grid(ii)%istart+1
          end if
        end do
        if(ishpf>10) then
          msg="ishpf>7"
          LIBPAW_ERROR(msg)
        end if
        LIBPAW_ALLOCATE(paw_setuploc%shape_function%data,(mesh_size,ishpf))
        do ii=1,ishpf
          paw_setuploc%shape_function%data(:,ii)=shpf(ii)%data(:)
          LIBPAW_DEALLOCATE(shpf(ii)%data)
        end do
        LIBPAW_DATATYPE_DEALLOCATE(shpf)

      case ("shape_function")
        in_data=.false.

      case ("pseudo_partial_wave")
        in_data=.false.

      case ("ae_partial_wave")
        in_data=.false.

      case ("projector_function")
        in_data=.false.

      case ("ae_core_density")
        in_data=.false.

      case ("pseudo_core_density")
        in_data=.false.

      case ("pseudo_valence_density")
        in_data=.false.

      case ("zero_potential")
        in_data=.false.

      case ("LDA_minus_half_potential")
        in_data=.false.

      case ("ae_core_kinetic_energy_density")
        in_data=.false.

      case ("pseudo_core_kinetic_energy_density")
        in_data=.false.

      case ("kresse_joubert_local_ionic_potential")
        in_data=.false.

      case ("blochl_local_ionic_potential")
        in_data=.false.

      case ("kinetic_energy_differences")
        in_data=.false.

end select

end subroutine paw_end_element1
!!***

!-------------------------------------------------------------------------

!!****f* m_pawxmlps/pawdata_chunk
!! NAME
!! pawdata_chunk
!!
!! FUNCTION
!!   Take a string and turn it into useful data structure (reals)
!!
!! INPUTS
!!   chunk=raw data for chunk of XML data
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!   Copied and translated into module data (side effect)
!!
!! PARENTS
!!
!! CHILDREN
!!      bound_deriv,paw_rdfromline,paw_spline,paw_splint,pawrad_init
!!
!! SOURCE
subroutine pawdata_chunk(chunk)

character(len=*),intent(in) :: chunk

integer :: ii,ntokens,status,last_pos
logical :: in_token
character(len=len(chunk))  :: str
character(len=50)  :: msg
character(len=1)  :: cc
real(dpxml),pointer :: x(:)


if (len_trim(chunk) == 0) RETURN     ! skip empty chunk

if (in_data) then
  str = chunk ; x => rp%data

! Check the contents of the string and find the number of tokens it contains
! The standard separator is generalized whitespace (space, tab, CR, or LF)
  in_token=.false.;ntokens=0;last_pos=0
  do ii=1,len_trim(str)
    cc=str(ii:ii)
    if (in_token) then
      if (cc==char(9).or.cc==char(10).or.cc==char(13).or.cc==char(32)) then
        in_token = .false.
        if (cc==char(10).or.cc==char(13)) str(ii:ii) = " "
      else
        last_pos=ii
      end if
    else
      if (cc==char(9).or.cc==char(10).or.cc==char(13).or.cc==char(32)) then
        if (cc==char(10).or.cc==char(13)) str(ii:ii) = " "
      else
        in_token=.true.
        last_pos=ii
        ntokens=ntokens + 1
      end if
    end if
  end do

  if ((ndata+ntokens)>size(x)) then 
    msg="data array full"
    LIBPAW_ERROR(msg)
  end if

! Take the string and turn it into useful reals
  read(unit=str(1:last_pos),fmt=*,iostat=status) x(ndata+1:ndata+ntokens)
  if (status/=0) then
    msg="real conversion error"
    LIBPAW_ERROR(msg)
  end if
  ndata=ndata+ntokens

end if

end subroutine pawdata_chunk
!!***

#endif

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!------------- ROUTINES AND FUNCTIONS FOR THE FORTRAN READER -------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!!****f* m_pawxmlps/paw_setup_free
!! NAME
!! paw_setup_free
!!
!! FUNCTION
!!  Destroy a paw_setup datastructure
!!
!! SIDE EFFECTS
!!  paw_setup<paw_setup_type>=Datatype gathering information on XML paw setup.
!!
!! PARENTS
!!      m_pspheads,m_pspini
!!
!! CHILDREN
!!      bound_deriv,paw_rdfromline,paw_spline,paw_splint,pawrad_init
!!
!! SOURCE

subroutine paw_setup_free(paw_setupin)

!Arguments ------------------------------------
!scalars
 type(paw_setup_t),intent(inout) :: paw_setupin

!Local variables-------------------------------
 integer :: ii

! *********************************************************************
 
 paw_setupin%tread=.false.
 paw_setupin%atom%tread=.false.
 paw_setupin%xc_functional%tread=.false.
 paw_setupin%generator%tread=.false.
 paw_setupin%valence_states%tread=.false.
 paw_setupin%shape_function%tread=.false.
 paw_setupin%ae_core_density%tread=.false.
 paw_setupin%pseudo_core_density%tread=.false.
 paw_setupin%pseudo_valence_density%tread=.false.
 paw_setupin%zero_potential%tread=.false.
 paw_setupin%LDA_minus_half_potential%tread=.false.
 paw_setupin%ae_core_kinetic_energy_density%tread=.false.
 paw_setupin%pseudo_core_kinetic_energy_density%tread=.false.
 paw_setupin%kresse_joubert_local_ionic_potential%tread=.false.
 paw_setupin%blochl_local_ionic_potential%tread=.false.
 paw_setupin%kinetic_energy_differences%tread=.false.
 paw_setupin%exact_exchange_matrix%tread=.false.

 if(allocated( paw_setupin%shape_function%data)) then
   LIBPAW_DEALLOCATE(paw_setupin%shape_function%data)
 end if
 if(allocated( paw_setupin%ae_core_density%data)) then
   LIBPAW_DEALLOCATE(paw_setupin%ae_core_density%data)
 end if
 if(allocated( paw_setupin%pseudo_core_density%data)) then
   LIBPAW_DEALLOCATE(paw_setupin%pseudo_core_density%data)
 end if
 if(allocated( paw_setupin%pseudo_valence_density%data)) then
   LIBPAW_DEALLOCATE(paw_setupin%pseudo_valence_density%data)
 end if
 if(allocated( paw_setupin%zero_potential%data)) then
   LIBPAW_DEALLOCATE(paw_setupin%zero_potential%data)
 end if
 if(allocated( paw_setupin%LDA_minus_half_potential%data)) then
   LIBPAW_DEALLOCATE(paw_setupin%LDA_minus_half_potential%data)
 end if
 if(allocated( paw_setupin%ae_core_kinetic_energy_density%data)) then
   LIBPAW_DEALLOCATE(paw_setupin%ae_core_kinetic_energy_density%data)
 end if
 if(allocated( paw_setupin%pseudo_core_kinetic_energy_density%data)) then
   LIBPAW_DEALLOCATE(paw_setupin%pseudo_core_kinetic_energy_density%data)
 end if
 if(allocated( paw_setupin%kresse_joubert_local_ionic_potential%data)) then
   LIBPAW_DEALLOCATE(paw_setupin%kresse_joubert_local_ionic_potential%data)
 end if
 if(allocated( paw_setupin%blochl_local_ionic_potential%data)) then
   LIBPAW_DEALLOCATE(paw_setupin%blochl_local_ionic_potential%data)
 end if
 if(allocated( paw_setupin%kinetic_energy_differences%data)) then
   LIBPAW_DEALLOCATE(paw_setupin%kinetic_energy_differences%data)
 end if
 if(allocated( paw_setupin%exact_exchange_matrix%data)) then
   LIBPAW_DEALLOCATE(paw_setupin%exact_exchange_matrix%data)
 end if
 if (allocated( paw_setupin%ae_partial_wave)) then
   do ii=1,paw_setupin%valence_states%nval
     if(allocated( paw_setupin%ae_partial_wave(ii)%data)) then
       LIBPAW_DEALLOCATE(paw_setupin%ae_partial_wave(ii)%data)
     end if
   end do
   LIBPAW_DATATYPE_DEALLOCATE(paw_setupin%ae_partial_wave)
 end if
 if (allocated( paw_setupin%pseudo_partial_wave)) then
   do ii=1,paw_setupin%valence_states%nval
     if(allocated( paw_setupin%pseudo_partial_wave(ii)%data)) then
       LIBPAW_DEALLOCATE(paw_setupin%pseudo_partial_wave(ii)%data)
     end if
   end do
   LIBPAW_DATATYPE_DEALLOCATE(paw_setupin%pseudo_partial_wave)
 end if
 if (allocated( paw_setupin%projector_function)) then
   do ii=1,paw_setupin%valence_states%nval
     if(allocated( paw_setupin%projector_function(ii)%data)) then
       LIBPAW_DEALLOCATE(paw_setupin%projector_function(ii)%data)
     end if
   end do
   LIBPAW_DATATYPE_DEALLOCATE(paw_setupin%projector_function)
 end if
 if (allocated( paw_setupin%projector_fit)) then
   LIBPAW_DATATYPE_DEALLOCATE(paw_setupin%projector_fit)
 end if
 if(allocated(paw_setupin%valence_states%state)) then
   LIBPAW_DATATYPE_DEALLOCATE(paw_setupin%valence_states%state)
 end if
 if(allocated( paw_setupin%radial_grid)) then
   LIBPAW_DATATYPE_DEALLOCATE(paw_setupin%radial_grid)
 end if

end subroutine paw_setup_free
!!***

!-------------------------------------------------------------------------

!!****f* m_pawxmlps/paw_setup_copy
!! NAME
!! paw_setup_copy
!!
!! FUNCTION
!!  Copy a paw_setup datastructure into another
!!
!! INPUTS
!!  
!!  paw_setupin<paw_setup_type>=input paw_setup datastructure
!!
!! OUTPUT
!!  paw_setupout<paw_setup_type>=output paw_setup datastructure
!!
!! PARENTS
!!
!! CHILDREN
!!      bound_deriv,paw_rdfromline,paw_spline,paw_splint,pawrad_init
!!
!! SOURCE

subroutine paw_setup_copy(paw_setupin,paw_setupout)

!Arguments ------------------------------------
!scalars
 type(paw_setup_t),intent(in) :: paw_setupin
 type(paw_setup_t),intent(out) :: paw_setupout

!Local variables-------------------------------
!scalars
 integer :: ii,sz1,sz2

! *********************************************************************

!scalars
 paw_setupout%version=paw_setupin%version
 paw_setupout%tread=paw_setupin%tread
 paw_setupout%ngrid=paw_setupin%ngrid
 paw_setupout%idgrid=paw_setupin%idgrid
 paw_setupout%optortho=paw_setupin%optortho
 paw_setupout%rpaw=paw_setupin%rpaw
 paw_setupout%ex_cc=paw_setupin%ex_cc
 paw_setupout%atom%tread=paw_setupin%atom%tread
 paw_setupout%atom%symbol=paw_setupin%atom%symbol
 paw_setupout%atom%znucl=paw_setupin%atom%znucl
 paw_setupout%atom%zion=paw_setupin%atom%zion
 paw_setupout%atom%zval=paw_setupin%atom%zval
 paw_setupout%xc_functional%tread=paw_setupin%xc_functional%tread
 paw_setupout%xc_functional%functionaltype=paw_setupin%xc_functional%functionaltype
 paw_setupout%xc_functional%name=paw_setupin%xc_functional%name
 paw_setupout%generator%tread=paw_setupin%generator%tread
 paw_setupout%generator%gen=paw_setupin%generator%gen
 paw_setupout%generator%name=paw_setupin%generator%name
 paw_setupout%valence_states%tread=paw_setupin%valence_states%tread
 paw_setupout%valence_states%nval=paw_setupin%valence_states%nval
 paw_setupout%shape_function%tread=paw_setupin%shape_function%tread
 paw_setupout%shape_function%gtype=paw_setupin%shape_function%gtype
 paw_setupout%shape_function%grid=paw_setupin%shape_function%grid
 paw_setupout%shape_function%rc=paw_setupin%shape_function%rc
 paw_setupout%shape_function%lamb=paw_setupin%shape_function%lamb
 paw_setupout%ae_core_density%tread=paw_setupin%ae_core_density%tread
 paw_setupout%ae_core_density%grid=paw_setupin%ae_core_density%grid
 paw_setupout%ae_core_density%state=paw_setupin%ae_core_density%state
 paw_setupout%pseudo_core_density%tread=paw_setupin%pseudo_core_density%tread
 paw_setupout%pseudo_core_density%grid=paw_setupin%pseudo_core_density%grid
 paw_setupout%pseudo_core_density%state=paw_setupin%pseudo_core_density%state
 paw_setupout%pseudo_valence_density%tread=paw_setupin%pseudo_valence_density%tread
 paw_setupout%pseudo_valence_density%grid=paw_setupin%pseudo_valence_density%grid
 paw_setupout%pseudo_valence_density%state=paw_setupin%pseudo_valence_density%state
 paw_setupout%zero_potential%tread=paw_setupin%zero_potential%tread
 paw_setupout%zero_potential%grid=paw_setupin%zero_potential%grid
 paw_setupout%zero_potential%state=paw_setupin%zero_potential%state
 paw_setupout%LDA_minus_half_potential%tread=paw_setupin%LDA_minus_half_potential%tread
 paw_setupout%LDA_minus_half_potential%grid=paw_setupin%LDA_minus_half_potential%grid
 paw_setupout%LDA_minus_half_potential%state=paw_setupin%LDA_minus_half_potential%state
 paw_setupout%ae_core_kinetic_energy_density%tread=&
&     paw_setupin%ae_core_kinetic_energy_density%tread
 paw_setupout%ae_core_kinetic_energy_density%grid=&
&     paw_setupin%ae_core_kinetic_energy_density%grid
 paw_setupout%ae_core_kinetic_energy_density%state=&
&     paw_setupin%ae_core_kinetic_energy_density%state
 paw_setupout%pseudo_core_kinetic_energy_density%tread=&
&     paw_setupin%pseudo_core_kinetic_energy_density%tread
 paw_setupout%pseudo_core_kinetic_energy_density%grid=&
&     paw_setupin%pseudo_core_kinetic_energy_density%grid
 paw_setupout%pseudo_core_kinetic_energy_density%state=&
&     paw_setupin%pseudo_core_kinetic_energy_density%state
 paw_setupout%kresse_joubert_local_ionic_potential%tread=&
&    paw_setupin%kresse_joubert_local_ionic_potential%tread
 paw_setupout%kresse_joubert_local_ionic_potential%grid=&
&    paw_setupin%kresse_joubert_local_ionic_potential%grid
 paw_setupout%kresse_joubert_local_ionic_potential%state=&
&    paw_setupin%kresse_joubert_local_ionic_potential%state
 paw_setupout%blochl_local_ionic_potential%tread=&
&    paw_setupin%blochl_local_ionic_potential%tread
 paw_setupout%blochl_local_ionic_potential%grid=&
&    paw_setupin%blochl_local_ionic_potential%grid
 paw_setupout%blochl_local_ionic_potential%state=&
&    paw_setupin%blochl_local_ionic_potential%state
 paw_setupout%kinetic_energy_differences%tread=paw_setupin%kinetic_energy_differences%tread
 paw_setupout%kinetic_energy_differences%grid=paw_setupin%kinetic_energy_differences%grid
 paw_setupout%kinetic_energy_differences%state=paw_setupin%kinetic_energy_differences%state
 paw_setupout%exact_exchange_matrix%tread=paw_setupin%exact_exchange_matrix%tread
! allocatable arrays
 if (allocated(paw_setupin%shape_function%data)) then
   sz1=size(paw_setupin%shape_function%data,1)
   sz2=size(paw_setupin%shape_function%data,2)
   LIBPAW_ALLOCATE(paw_setupout%shape_function%data,(sz1,sz2))
   paw_setupout%shape_function%data=paw_setupin%shape_function%data
 end if
 if (allocated(paw_setupin%ae_core_density%data)) then
   sz1=size(paw_setupin%ae_core_density%data,1)
   LIBPAW_ALLOCATE(paw_setupout%ae_core_density%data,(sz1))
   paw_setupout%ae_core_density%data=paw_setupin%ae_core_density%data
 end if
 if (allocated(paw_setupin%pseudo_core_density%data)) then
   sz1=size(paw_setupin%pseudo_core_density%data,1)
   LIBPAW_ALLOCATE(paw_setupout%pseudo_core_density%data,(sz1))
   paw_setupout%pseudo_core_density%data=paw_setupin%pseudo_core_density%data
 end if
 if (allocated(paw_setupin%pseudo_valence_density%data)) then
   sz1=size(paw_setupin%pseudo_valence_density%data,1)
   LIBPAW_ALLOCATE(paw_setupout%pseudo_valence_density%data,(sz1))
   paw_setupout%pseudo_valence_density%data=paw_setupin%pseudo_valence_density%data
 end if
 if (allocated(paw_setupin%zero_potential%data)) then
   sz1=size(paw_setupin%zero_potential%data,1)
   LIBPAW_ALLOCATE(paw_setupout%zero_potential%data,(sz1))
   paw_setupout%zero_potential%data=paw_setupin%zero_potential%data
 end if
 if (allocated(paw_setupin%LDA_minus_half_potential%data)) then
   sz1=size(paw_setupin%LDA_minus_half_potential%data,1)
   LIBPAW_ALLOCATE(paw_setupout%LDA_minus_half_potential%data,(sz1))
   paw_setupout%LDA_minus_half_potential%data=paw_setupin%LDA_minus_half_potential%data
 end if
 if (allocated(paw_setupin%ae_core_kinetic_energy_density%data)) then
   sz1=size(paw_setupin%ae_core_kinetic_energy_density%data,1)
   LIBPAW_ALLOCATE(paw_setupout%ae_core_kinetic_energy_density%data,(sz1))
   paw_setupout%ae_core_kinetic_energy_density%data=paw_setupin%ae_core_kinetic_energy_density%data
 end if
 if (allocated(paw_setupin%pseudo_core_kinetic_energy_density%data)) then
   sz1=size(paw_setupin%pseudo_core_kinetic_energy_density%data,1)
   LIBPAW_ALLOCATE(paw_setupout%pseudo_core_kinetic_energy_density%data,(sz1))
   paw_setupout%pseudo_core_kinetic_energy_density%data=paw_setupin%pseudo_core_kinetic_energy_density%data
 end if
 if (allocated(paw_setupin%kresse_joubert_local_ionic_potential%data)) then
   sz1=size(paw_setupin%kresse_joubert_local_ionic_potential%data,1)
   LIBPAW_ALLOCATE(paw_setupout%kresse_joubert_local_ionic_potential%data,(sz1))
   paw_setupout%kresse_joubert_local_ionic_potential%data=paw_setupin%kresse_joubert_local_ionic_potential%data
 end if
 if (allocated(paw_setupin%blochl_local_ionic_potential%data)) then
   sz1=size(paw_setupin%blochl_local_ionic_potential%data,1)
   LIBPAW_ALLOCATE(paw_setupout%blochl_local_ionic_potential%data,(sz1))
   paw_setupout%blochl_local_ionic_potential%data=paw_setupin%blochl_local_ionic_potential%data
 end if
 if (allocated(paw_setupin%exact_exchange_matrix%data)) then
   sz1=size(paw_setupin%exact_exchange_matrix%data,1)
   LIBPAW_ALLOCATE(paw_setupout%exact_exchange_matrix%data,(sz1))
   paw_setupout%exact_exchange_matrix%data=paw_setupin%exact_exchange_matrix%data
 end if
 if (allocated(paw_setupin%kinetic_energy_differences%data)) then
   sz1=size(paw_setupin%kinetic_energy_differences%data,1)
   LIBPAW_ALLOCATE(paw_setupout%kinetic_energy_differences%data,(sz1))
   paw_setupout%kinetic_energy_differences%data=paw_setupin%kinetic_energy_differences%data
 end if
 if(allocated( paw_setupin%radial_grid)) then
   sz1=size(paw_setupin%radial_grid,1)
   LIBPAW_DATATYPE_ALLOCATE(paw_setupout%radial_grid,(sz1))
   paw_setupout%radial_grid=paw_setupin%radial_grid
 end if
 if(allocated(paw_setupin%valence_states%state)) then
   sz1=size(paw_setupin%valence_states%state,1)
   LIBPAW_DATATYPE_ALLOCATE(paw_setupout%valence_states%state,(sz1))
   paw_setupout%valence_states%state=paw_setupin%valence_states%state
 end if

 if (allocated( paw_setupin%ae_partial_wave)) then
   sz1=size(paw_setupin%ae_partial_wave,1)
   LIBPAW_DATATYPE_ALLOCATE(paw_setupout%ae_partial_wave,(sz1))
   do ii=1,paw_setupin%valence_states%nval
     paw_setupout%ae_partial_wave(ii)%tread=paw_setupin%ae_partial_wave(ii)%tread
     paw_setupout%ae_partial_wave(ii)%grid=paw_setupin%ae_partial_wave(ii)%grid
     paw_setupout%ae_partial_wave(ii)%state=paw_setupin%ae_partial_wave(ii)%state
     if(allocated( paw_setupin%ae_partial_wave(ii)%data)) then
       sz1=size(paw_setupin%ae_partial_wave(ii)%data,1)
       LIBPAW_ALLOCATE(paw_setupout%ae_partial_wave(ii)%data,(sz1))
       paw_setupout%ae_partial_wave(ii)%data=paw_setupin%ae_partial_wave(ii)%data
     end if
   end do
 end if 
 if (allocated( paw_setupin%pseudo_partial_wave)) then
   sz1=size(paw_setupin%pseudo_partial_wave,1)
   LIBPAW_DATATYPE_ALLOCATE(paw_setupout%pseudo_partial_wave,(sz1))
   do ii=1,paw_setupin%valence_states%nval
     paw_setupout%pseudo_partial_wave(ii)%tread=paw_setupin%pseudo_partial_wave(ii)%tread
     paw_setupout%pseudo_partial_wave(ii)%grid=paw_setupin%pseudo_partial_wave(ii)%grid
     paw_setupout%pseudo_partial_wave(ii)%state=paw_setupin%pseudo_partial_wave(ii)%state
     if(allocated( paw_setupin%pseudo_partial_wave(ii)%data)) then
       sz1=size(paw_setupin%pseudo_partial_wave(ii)%data,1)
       LIBPAW_ALLOCATE(paw_setupout%pseudo_partial_wave(ii)%data,(sz1))
       paw_setupout%pseudo_partial_wave(ii)%data=paw_setupin%pseudo_partial_wave(ii)%data
     end if
   end do
 end if 
  if (allocated( paw_setupin%projector_function)) then
   sz1=size(paw_setupin%projector_function,1)
   LIBPAW_DATATYPE_ALLOCATE(paw_setupout%projector_function,(sz1))
   do ii=1,paw_setupin%valence_states%nval
     paw_setupout%projector_function(ii)%tread=paw_setupin%projector_function(ii)%tread
     paw_setupout%projector_function(ii)%grid=paw_setupin%projector_function(ii)%grid
     paw_setupout%projector_function(ii)%state=paw_setupin%projector_function(ii)%state
     if(allocated( paw_setupin%projector_function(ii)%data)) then
       sz1=size(paw_setupin%projector_function(ii)%data,1)
       LIBPAW_ALLOCATE(paw_setupout%projector_function(ii)%data,(sz1))
       paw_setupout%projector_function(ii)%data=paw_setupin%projector_function(ii)%data
     end if
   end do
 end if 
  if (allocated( paw_setupin%projector_fit)) then
   sz1=size(paw_setupin%projector_fit,1)
   LIBPAW_DATATYPE_ALLOCATE(paw_setupout%projector_fit,(sz1))
   do ii=1,paw_setupin%valence_states%nval
     paw_setupout%projector_fit(ii)%tread=paw_setupin%projector_fit(ii)%tread
     paw_setupout%projector_fit(ii)%ngauss=paw_setupin%projector_fit(ii)%ngauss
     paw_setupout%projector_fit(ii)%state=paw_setupin%projector_fit(ii)%state
     paw_setupout%projector_fit(ii)%factors=paw_setupin%projector_fit(ii)%factors
     paw_setupout%projector_fit(ii)%expos=paw_setupin%projector_fit(ii)%expos
   end do
 end if 

end subroutine paw_setup_copy
!!***

!-------------------------------------------------------------------------

!!****f* m_pawxmlps/paw_rdfromline
!! NAME
!! paw_rdfromline
!!
!! FUNCTION
!! Read the value of a keyword from a XML line
!!
!! INPUTS
!!  keyword= keyword which value has to be read
!!  line= string from which the data are read (line from a XML)
!!
!! OUTPUT
!!  ierr= error code
!!  output= (string) value of the keyword
!!
!! PARENTS
!!      m_pawxmlps
!!
!! CHILDREN
!!      bound_deriv,paw_rdfromline,paw_spline,paw_splint,pawrad_init
!!
!! SOURCE

 subroutine paw_rdfromline(keyword,line,output,ierr)

!Arguments ---------------------------------------------
  character(len=*), intent(in) :: keyword,line
  character(len=*), intent(out) :: output
  integer, intent(out) :: ierr
!Local variables ---------------------------------------
  character(len=len(line)) :: temp
  integer :: pos,pos2

! *********************************************************************

 ierr=1;output=""
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

 end subroutine paw_rdfromline
!!***

!-------------------------------------------------------------------------

!!****f* m_pawxmlps/rdpawpsxml_header
!! NAME
!! rdpawpsxml_header
!!
!! FUNCTION
!! Read the header of a PAW pseudopotential XML file generated by AtomPAW
!!
!! INPUTS
!!  filename= input file name (atomicdata XML)
!!
!! OUTPUT
!!  paw_setup=pseudopotential data structure
!!
!! PARENTS
!!      m_pspheads
!!
!! CHILDREN
!!      bound_deriv,paw_rdfromline,paw_spline,paw_splint,pawrad_init
!!
!! SOURCE

 subroutine rdpawpsxml_header(ecut_tmp,filename,paw_setup)

!Arguments ---------------------------------------------
 
 character (len=fnlen),intent(in) :: filename
 real(dp), intent(inout) :: ecut_tmp(3,2)
 type(paw_setup_t),intent(inout) :: paw_setup
!Local variables ---------------------------------------
 integer :: funit,ii,ir,igrid,ival,ierr,ishpf,lmax,mesh_size
 logical :: endfile,found
 character(len=100) :: msg
 character (len=XML_RECL) :: line,readline
 character (len=XML_RECL) :: strg
 character (len=30) :: strg1
 real(dp) :: rc(6)
 real(dp), allocatable :: shpf(:,:)
 type(state_t), pointer :: valstate (:)
 type(radial_grid_t), pointer :: grids (:)

! *************************************************************************

!Open the atomicdata XML file for reading
 open(newunit=funit,file=filename,form='formatted',status='old', recl=XML_RECL)

!Start a reading loop
 endfile=.false.
 found=.false.
 paw_setup%rpaw=-1.d0
 rc=-1.d0

 do while ((.not.endfile).and.(.not.found))
   read(funit,'(a)',err=10,end=10) readline
   line=adjustl(readline);goto 20
   10 line="";endfile=.true.
   20 continue

!  --Read VERSION
   if ((line(1:10)=='<paw_setup').or.(line(1:12)=='<paw_dataset')) then
     paw_setup%tread=.true.
     igrid=0;ishpf=0
     LIBPAW_DATATYPE_ALLOCATE(grids,(10))

     call paw_rdfromline(" version",line,strg,ierr)
     paw_setup%version=trim(strg)
     cycle
   end if

!  --Read TITLE, ATOMIC CHARGE AND CORE CHARGE
   if (line(1:6)=='<atom ') then
     paw_setup%atom%tread=.true.
     call paw_rdfromline(" symbol",line,strg,ierr)
     paw_setup%atom%symbol=trim(strg)
     call paw_rdfromline(" Z",line,strg,ierr)
     if (len(trim(strg))<=30) then
       strg1=trim(strg)
       read(unit=strg1,fmt=*) paw_setup%atom%znucl
     else
       read(unit=strg,fmt=*) paw_setup%atom%znucl
     end if
     call paw_rdfromline(" core",line,strg,ierr)
     if (len(trim(strg))<=30) then
       strg1=trim(strg)
       read(unit=strg1,fmt=*) paw_setup%atom%zion
     else
       read(unit=strg,fmt=*) paw_setup%atom%zion
     end if
     call paw_rdfromline(" valence",line,strg,ierr)
     if (len(trim(strg))<=30) then
       strg1=trim(strg)
       read(unit=strg1,fmt=*) paw_setup%atom%zval
     else
       read(unit=strg,fmt=*) paw_setup%atom%zval
     end if
     cycle
   end if

!  --Read Ecut and Ecutdg
   if (line(1:8)=='<pw_ecut') then
     call paw_rdfromline(" low",line,strg,ierr)
     if (len(trim(strg))<=30) then
       strg1=trim(strg)
       read(unit=strg1,fmt=*) ecut_tmp(1,1)
     else
       read(unit=strg,fmt=*) ecut_tmp(1,1)
     end if
     call paw_rdfromline(" medium",line,strg,ierr)
     if (len(trim(strg))<=30) then
       strg1=trim(strg)
       read(unit=strg1,fmt=*) ecut_tmp(2,1)
     else
       read(unit=strg,fmt=*) ecut_tmp(2,1)
     end if
     call paw_rdfromline(" high",line,strg,ierr)
     if (len(trim(strg))<=30) then
       strg1=trim(strg)
       read(unit=strg1,fmt=*) ecut_tmp(3,1)
     else
       read(unit=strg,fmt=*) ecut_tmp(3,1)
     end if
     cycle
   end if

!  --Read EXCHANGE-CORRELATION TYPE
   if (line(1:14)=='<xc_functional') then
     paw_setup%xc_functional%tread=.true.
     call paw_rdfromline(" type",line,strg,ierr)
     paw_setup%xc_functional%functionaltype = trim(strg)
     call paw_rdfromline(" name",line,strg,ierr)
     paw_setup%xc_functional%name= trim(strg)
     cycle
   end if

!  --Read GENERATOR
   if (line(1:10)=='<generator') then
     paw_setup%generator%tread=.true.
     call paw_rdfromline(" type",line,strg,ierr)
     paw_setup%generator%gen  = trim(strg)
     call paw_rdfromline(" name",line,strg,ierr)
     paw_setup%generator%name= trim(strg)
     cycle
   end if

!  --Read PAW RADIUS
   if (line(1:11)=='<PAW_radius') then
     call paw_rdfromline(" rpaw",line,strg,ierr)
     if (len(trim(strg))<=30) then
       strg1=trim(strg)
       read(unit=strg1,fmt=*) paw_setup%rpaw
     else
       read(unit=strg,fmt=*) paw_setup%rpaw
     end if
     cycle
   end if
   if (line(1:11)=='<paw_radius') then
     call paw_rdfromline(" rc",line,strg,ierr)
     if (len(trim(strg))<=30) then
       strg1=trim(strg)
       read(unit=strg1,fmt=*) paw_setup%rpaw
     else
       read(unit=strg,fmt=*) paw_setup%rpaw
     end if
     cycle
   end if

!  --Read BASIS SIZE, ORBITALS, RC AND OCCUPATIONS/STATE IDs
   if (line(1:16)=='<valence_states>') then
     paw_setup%valence_states%tread=.true.
     LIBPAW_DATATYPE_ALLOCATE(valstate,(50))
     ival=0
     lmax=0
     do while (line(1:17)/='</valence_states>')
       read(funit,'(a)') readline;line=adjustl(readline)
       if (line(1:6)=='<state') then
         ival=ival+1
         if (ival>50) then
           close(funit)
           msg="Error in rdpawps1xml: basis size too large (>50)!"
           LIBPAW_ERROR(msg)
         end if
         call paw_rdfromline(" n",line,strg,ierr)
         if (strg == "" ) then 
           valstate(ival)%nn=-1    
         else
           if (len(trim(strg))<=30) then
             strg1=trim(strg)
             read(unit=strg1,fmt=*) valstate(ival)%nn
           else
             read(unit=strg,fmt=*) valstate(ival)%nn
           end if
         end if
         call paw_rdfromline(" l",line,strg,ierr)
         if (len(trim(strg))<=30) then
           strg1=trim(strg)
           read(unit=strg1,fmt=*) valstate(ival)%ll
         else
           read(unit=strg,fmt=*) valstate(ival)%ll
         end if
         if(valstate(ival)%ll>lmax) lmax=valstate(ival)%ll
         call paw_rdfromline(" f",line,strg,ierr)
         if (strg == "" ) then 
           valstate(ival)%ff=-1.d0
         else
           if (len(trim(strg))<=30) then
             strg1=trim(strg)
             read(unit=strg1,fmt=*) valstate(ival)%ff
           else
             read(unit=strg,fmt=*) valstate(ival)%ff
           end if
         end if
         call paw_rdfromline(" rc",line,strg,ierr)
         if (len(trim(strg))<=30) then
           strg1=trim(strg)
           read(unit=strg1,fmt=*) valstate(ival)%rc
         else
           read(unit=strg,fmt=*) valstate(ival)%rc
         end if
         call paw_rdfromline(" e",line,strg,ierr)
         if (len(trim(strg))<=30) then
           strg1=trim(strg)
           read(unit=strg1,fmt=*) valstate(ival)%ee
         else
           read(unit=strg,fmt=*) valstate(ival)%ee
         end if
         call paw_rdfromline(" id",line,strg,ierr)
         valstate(ival)%id = trim(strg)
       end if
     end do
     cycle
   end if

!  --Read MESH_STEP AND NUMBER OF POINTS
   if (line(1:12)=='<radial_grid')then
     igrid=igrid+1
     call paw_rdfromline(" eq",line,strg,ierr)
     grids(igrid)%eq = trim(strg)
     call paw_rdfromline(" a",line,strg,ierr)
     if (strg == "" ) then
       grids(igrid)%aa=0.d0
     else
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) grids(igrid)%aa
       else
         read(unit=strg,fmt=*) grids(igrid)%aa
       end if
     end if
     call paw_rdfromline(" n",line,strg,ierr)
     if (strg == "" ) then
       grids(igrid)%nn=0
     else
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) grids(igrid)%nn
       else
         read(unit=strg,fmt=*) grids(igrid)%nn
       end if
     end if
     call paw_rdfromline(" d",line,strg,ierr)
     if (strg == "" ) then
       grids(igrid)%dd=0.d0
     else
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) grids(igrid)%dd
       else
         read(unit=strg,fmt=*) grids(igrid)%dd
       end if
     end if
     call paw_rdfromline(" b",line,strg,ierr)
     if (strg == "" ) then
       grids(igrid)%bb=0.d0
     else
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) grids(igrid)%bb
       else
         read(unit=strg,fmt=*) grids(igrid)%bb
       end if
     end if
     call paw_rdfromline("istart",line,strg,ierr)
     if (len(trim(strg))<=30) then
       strg1=trim(strg)
       read(unit=strg1,fmt=*) grids(igrid)%istart
     else
       read(unit=strg,fmt=*) grids(igrid)%istart
     end if
     call paw_rdfromline("iend",line,strg,ierr)
     if (len(trim(strg))<=30) then
       strg1=trim(strg)
       read(unit=strg1,fmt=*) grids(igrid)%iend
     else
       read(unit=strg,fmt=*) grids(igrid)%iend
     end if
     call paw_rdfromline(" id",line,strg,ierr)
     grids(igrid)%id = trim(strg)
     if(igrid>10)then
       close(funit)
       msg="igrid>10"
       LIBPAW_ERROR(msg)
     end if
     cycle
   end if

!  --Read SHAPE TYPE
   if (line(1:15)=='<shape_function') then
     paw_setup%shape_function%tread=.true.
     call paw_rdfromline(" type",line,strg,ierr)
     paw_setup%shape_function%gtype = trim(strg)
     call paw_rdfromline(" rc",line,strg,ierr)
     if (strg /= "" ) then
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) paw_setup%shape_function%rc
       else
         read(unit=strg,fmt=*) paw_setup%shape_function%rc
       end if
     end if
     call paw_rdfromline(" lamb",line,strg,ierr)
     if (strg == "" ) then
       paw_setup%shape_function%lamb=0
     else
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) paw_setup%shape_function%lamb
       else
         read(unit=strg,fmt=*) paw_setup%shape_function%lamb
       end if
     end if
     found=paw_setup%shape_function%tread
     call paw_rdfromline("grid",line,strg,ierr)
     paw_setup%shape_function%grid=trim(strg)
     if (strg /= "" ) then
       paw_setup%shape_function%gtype ="num" 
       do ii=1,igrid
         if(trim(paw_setup%shape_function%grid)==trim(grids(ii)%id)) then
           mesh_size=grids(ii)%iend-grids(ii)%istart+1
           exit
         end if
       end do
       if(.not.allocated(shpf)) then
         LIBPAW_ALLOCATE(shpf,(mesh_size,7))
       end if
       ishpf=ishpf+1
       read(funit,*) (shpf(ir,ishpf),ir=1,mesh_size)
       call paw_rdfromline(" l",line,strg,ierr)
       if (strg /= "" ) then
         found=.false.
         if(paw_setup%valence_states%tread) then
           if(ishpf==2*lmax+1) found=.true.
         else
           write(msg,'(a,a,a)')"the grids and the states must be read before the shapefunction",ch10,&
&           "Action: Modify your XML PAW data file"
           LIBPAW_ERROR(msg)
         end if
       end if
     end if
     cycle
   end if

!  End of reading loop
 end do
 if(ival>0)then
   LIBPAW_DATATYPE_ALLOCATE(paw_setup%valence_states%state,(ival))
   paw_setup%valence_states%state(ival)%tread=.true.
   paw_setup%valence_states%nval=ival
   do ii=1,ival
     paw_setup%valence_states%state(ii)=valstate(ii)
   end do
 end if
 LIBPAW_DATATYPE_DEALLOCATE(valstate)
 LIBPAW_DATATYPE_ALLOCATE(paw_setup%radial_grid,(igrid))
 paw_setup%radial_grid(igrid)%tread=.true.
 paw_setup%ngrid=igrid
 do ii=1,igrid
   paw_setup%radial_grid(ii)=grids(ii)
 end do
 LIBPAW_DATATYPE_DEALLOCATE(grids)
 if(allocated(shpf)) then
   LIBPAW_DEALLOCATE(shpf)
 end if
 close(funit)

 end subroutine rdpawpsxml_header
!!***

!-------------------------------------------------------------------------

!!****f* m_pawxmlps/rdpawpsxml
!! NAME
!! rdpawpsxml
!!
!! FUNCTION
!! Read the PAW pseudopotential XML file generated by AtomPAW
!!
!! INPUTS
!!  filename= input file name (atomicdata XML)
!!
!! OUTPUT
!!  paw_setup=pseudopotential data structure
!!
!! PARENTS
!!      m_pspheads
!!
!! CHILDREN
!!      bound_deriv,paw_rdfromline,paw_spline,paw_splint,pawrad_init
!!
!! SOURCE

 subroutine rdpawpsxml(filename,paw_setup)

!Arguments ---------------------------------------------
 character (len=fnlen),intent(in) :: filename
 type(paw_setup_t),intent(inout) :: paw_setup
!Local variables ---------------------------------------
 integer :: funit, iaewf,ii,ipswf,iproj,ir,igrid,ival,ierr,ishpf,lmax,mesh_size,igauss,iprojfit
 logical :: endfile,found,endgauss
 character(len=100) :: msg
 character (len=XML_RECL) :: line,readline
 character (len=XML_RECL) :: strg
 character (len=30) :: strg1
 real(dp) :: rc(6)
 real(dp), allocatable :: shpf(:,:)
 type(state_t), pointer :: valstate (:)
 type(radial_grid_t), pointer :: grids (:)

! *************************************************************************

!Open the atomicdata XML file for reading
 open(newunit=funit,file=filename,form='formatted',status='old', recl=XML_RECL)

!Start a reading loop
 endfile=.false.
 found=.false.
 paw_setup%rpaw=-1.d0
 rc=-1.d0

 do while ((.not.endfile).and.(.not.found))
   read(funit,'(a)',err=10,end=10) readline
   line=adjustl(readline);goto 20
   10 line="";endfile=.true.
   20 continue

!  --Read VERSION
   if ((line(1:10)=='<paw_setup').or.(line(1:12)=='<paw_dataset')) then
     paw_setup%tread=.true.
     igrid=0;ishpf=0
     LIBPAW_DATATYPE_ALLOCATE(grids,(10))

     call paw_rdfromline(" version",line,strg,ierr)
     paw_setup%version=trim(strg)
     cycle
   end if

!  --Read TITLE, ATOMIC CHARGE AND CORE CHARGE
   if (line(1:6)=='<atom ') then
     paw_setup%atom%tread=.true.
     call paw_rdfromline(" symbol",line,strg,ierr)
     paw_setup%atom%symbol=trim(strg)
     call paw_rdfromline(" Z",line,strg,ierr)
     if (len(trim(strg))<=30) then
       strg1=trim(strg)
       read(unit=strg1,fmt=*) paw_setup%atom%znucl
     else
       read(unit=strg,fmt=*) paw_setup%atom%znucl
     end if
     call paw_rdfromline(" core",line,strg,ierr)
     if (len(trim(strg))<=30) then
       strg1=trim(strg)
       read(unit=strg1,fmt=*) paw_setup%atom%zion
     else
       read(unit=strg,fmt=*) paw_setup%atom%zion
     end if
     call paw_rdfromline(" valence",line,strg,ierr)
     if (len(trim(strg))<=30) then
       strg1=trim(strg)
       read(unit=strg1,fmt=*) paw_setup%atom%zval
     else
       read(unit=strg,fmt=*) paw_setup%atom%zval
     end if
     cycle
   end if

!  --Read EXCHANGE-CORRELATION TYPE
   if (line(1:14)=='<xc_functional') then
     paw_setup%xc_functional%tread=.true.
     call paw_rdfromline(" type",line,strg,ierr)
     paw_setup%xc_functional%functionaltype = trim(strg)
     call paw_rdfromline(" name",line,strg,ierr)
     paw_setup%xc_functional%name= trim(strg)
     cycle
   end if

!  --Read GENERATOR
   if (line(1:10)=='<generator') then
     paw_setup%generator%tread=.true.
     call paw_rdfromline(" type",line,strg,ierr)
     paw_setup%generator%gen  = trim(strg)
     call paw_rdfromline(" name",line,strg,ierr)
     paw_setup%generator%name= trim(strg)
     cycle
   end if

!  --Read PAW RADIUS
   if (line(1:11)=='<PAW_radius') then
     call paw_rdfromline(" rpaw",line,strg,ierr)
     if (len(trim(strg))<=30) then
       strg1=trim(strg)
       read(unit=strg1,fmt=*) paw_setup%rpaw
     else
       read(unit=strg,fmt=*) paw_setup%rpaw
     end if
     cycle
   end if
   if (line(1:11)=='<paw_radius') then
     call paw_rdfromline(" rc",line,strg,ierr)
     if (len(trim(strg))<=30) then
       strg1=trim(strg)
       read(unit=strg1,fmt=*) paw_setup%rpaw
     else
       read(unit=strg,fmt=*) paw_setup%rpaw
     end if
     cycle
   end if

!  --Read BASIS SIZE, ORBITALS, RC AND OCCUPATIONS/STATE IDs
   if (line(1:16)=='<valence_states>') then
     paw_setup%valence_states%tread=.true.
     LIBPAW_DATATYPE_ALLOCATE(valstate,(50))
     ival=0
     lmax=0
     do while (line(1:17)/='</valence_states>')
       read(funit,'(a)') readline;line=adjustl(readline)
       if (line(1:6)=='<state') then
         ival=ival+1
         if (ival>50) then
           close(funit)
           msg="Error in rdpawps1xml: basis size too large (>50)!"
           LIBPAW_ERROR(msg)
         end if
         call paw_rdfromline(" n",line,strg,ierr)
         if (strg == "" ) then 
           valstate(ival)%nn=-1    
         else
           if (len(trim(strg))<=30) then
             strg1=trim(strg)
             read(unit=strg1,fmt=*) valstate(ival)%nn
           else
             read(unit=strg,fmt=*) valstate(ival)%nn
           end if
         end if
         call paw_rdfromline(" l",line,strg,ierr)
         if (len(trim(strg))<=30) then
           strg1=trim(strg)
           read(unit=strg1,fmt=*) valstate(ival)%ll
         else
           read(unit=strg,fmt=*) valstate(ival)%ll
         end if
         if(valstate(ival)%ll>lmax) lmax=valstate(ival)%ll
         call paw_rdfromline(" f",line,strg,ierr)
         if (strg == "" ) then 
           valstate(ival)%ff=-1.d0
         else
           if (len(trim(strg))<=30) then
             strg1=trim(strg)
             read(unit=strg1,fmt=*) valstate(ival)%ff
           else
             read(unit=strg,fmt=*) valstate(ival)%ff
           end if
         end if
         call paw_rdfromline(" rc",line,strg,ierr)
         if (len(trim(strg))<=30) then
           strg1=trim(strg)
           read(unit=strg1,fmt=*) valstate(ival)%rc
         else
           read(unit=strg,fmt=*) valstate(ival)%rc
         end if
         call paw_rdfromline(" e",line,strg,ierr)
         if (len(trim(strg))<=30) then
           strg1=trim(strg)
           read(unit=strg1,fmt=*) valstate(ival)%ee
         else
           read(unit=strg,fmt=*) valstate(ival)%ee
         end if
         call paw_rdfromline(" id",line,strg,ierr)
         valstate(ival)%id = trim(strg)
       end if
     end do
     cycle
   end if

!  --Read MESH_STEP AND NUMBER OF POINTS
   if (line(1:12)=='<radial_grid')then
     igrid=igrid+1
     call paw_rdfromline(" eq",line,strg,ierr)
     grids(igrid)%eq = trim(strg)
     call paw_rdfromline(" a",line,strg,ierr)
     if (strg == "" ) then
       grids(igrid)%aa=0.d0
     else
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) grids(igrid)%aa
       else
         read(unit=strg,fmt=*) grids(igrid)%aa
       end if
     end if
     call paw_rdfromline(" n",line,strg,ierr)
     if (strg == "" ) then
       grids(igrid)%nn=0
     else
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) grids(igrid)%nn
       else
         read(unit=strg,fmt=*) grids(igrid)%nn
       end if
     end if
     call paw_rdfromline(" d",line,strg,ierr)
     if (strg == "" ) then
       grids(igrid)%dd=0.d0
     else
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) grids(igrid)%dd
       else
         read(unit=strg,fmt=*) grids(igrid)%dd
       end if
     end if
     call paw_rdfromline(" b",line,strg,ierr)
     if (strg == "" ) then
       grids(igrid)%bb=0.d0
     else
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) grids(igrid)%bb
       else
         read(unit=strg,fmt=*) grids(igrid)%bb
       end if
     end if
     call paw_rdfromline("istart",line,strg,ierr)
     if (len(trim(strg))<=30) then
       strg1=trim(strg)
       read(unit=strg1,fmt=*) grids(igrid)%istart
     else
       read(unit=strg,fmt=*) grids(igrid)%istart
     end if
     call paw_rdfromline("iend",line,strg,ierr)
     if (len(trim(strg))<=30) then
       strg1=trim(strg)
       read(unit=strg1,fmt=*) grids(igrid)%iend
     else
       read(unit=strg,fmt=*) grids(igrid)%iend
     end if
     call paw_rdfromline(" id",line,strg,ierr)
     grids(igrid)%id = trim(strg)
     if(igrid>10)then
       close(funit)
       msg="igrid>10"
       LIBPAW_ERROR(msg)
     end if
     cycle
   end if

!  --Read SHAPE TYPE
   if (line(1:15)=='<shape_function') then
     paw_setup%shape_function%tread=.true.
     call paw_rdfromline(" type",line,strg,ierr)
     paw_setup%shape_function%gtype = trim(strg)
     call paw_rdfromline(" rc",line,strg,ierr)
     if (strg /= "" ) then
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) paw_setup%shape_function%rc
       else
         read(unit=strg,fmt=*) paw_setup%shape_function%rc
       end if
     end if
     call paw_rdfromline(" lamb",line,strg,ierr)
     if (strg == "" ) then
       paw_setup%shape_function%lamb=0
     else
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) paw_setup%shape_function%lamb
       else
         read(unit=strg,fmt=*) paw_setup%shape_function%lamb
       end if
     end if
     found=paw_setup%shape_function%tread
     call paw_rdfromline("grid",line,strg,ierr)
     paw_setup%shape_function%grid=trim(strg)
     if (strg /= "" ) then
       paw_setup%shape_function%gtype ="num" 
       do ii=1,igrid
         if(trim(paw_setup%shape_function%grid)==trim(grids(ii)%id)) then
           mesh_size=grids(ii)%iend-grids(ii)%istart+1
           exit
         end if
       end do
       if(.not.allocated(shpf)) then
         LIBPAW_ALLOCATE(shpf,(mesh_size,7))
       end if
       ishpf=ishpf+1
       read(funit,*) (shpf(ir,ishpf),ir=1,mesh_size)
       call paw_rdfromline(" l",line,strg,ierr)
       if (strg /= "" ) then
         found=.false.
         if(paw_setup%valence_states%tread) then
           if(ishpf==2*lmax+1) found=.true.
         else
           write(msg,'(a,a,a)')"the grids and the states must be read before the shapefunction",ch10,&
&           "Action: Modify your XML PAW data file"
           LIBPAW_ERROR(msg)
         end if
       end if
     end if
     cycle
   end if

!  End of reading loop
 end do

 if(igrid==0.or.ival==0) then
   write(msg,'(a,a,a)')"the grids and the states must be read before the shapefunction",ch10,&
&   "Action: Modify your XML PAW data file"
   LIBPAW_ERROR(msg)
 end if
 if(ishpf>0)then
   LIBPAW_ALLOCATE(paw_setup%shape_function%data,(mesh_size,ishpf))
   do ii=1,ishpf
     paw_setup%shape_function%data(:,ii)=shpf(:,ii)
   end do
   LIBPAW_DEALLOCATE(shpf)
 end if

 if(ival>0)then
   LIBPAW_DATATYPE_ALLOCATE(paw_setup%valence_states%state,(ival))
   paw_setup%valence_states%state(ival)%tread=.true.
   paw_setup%valence_states%nval=ival
   do ii=1,ival
     paw_setup%valence_states%state(ii)=valstate(ii)
   end do
 end if
 LIBPAW_DATATYPE_DEALLOCATE(valstate)
 if(.not.allocated(paw_setup%ae_partial_wave)) then
   LIBPAW_DATATYPE_ALLOCATE(paw_setup%ae_partial_wave,(paw_setup%valence_states%nval))
 end if
 if(.not.allocated(paw_setup%pseudo_partial_wave)) then
   LIBPAW_DATATYPE_ALLOCATE(paw_setup%pseudo_partial_wave,(paw_setup%valence_states%nval))
 end if
 if(.not.allocated(paw_setup%projector_function)) then
   LIBPAW_DATATYPE_ALLOCATE(paw_setup%projector_function,(paw_setup%valence_states%nval))
 end if

 LIBPAW_DATATYPE_ALLOCATE(paw_setup%radial_grid,(igrid))
 paw_setup%radial_grid(igrid)%tread=.true.
 paw_setup%ngrid=igrid
 do ii=1,igrid
   paw_setup%radial_grid(ii)=grids(ii)
 end do
 LIBPAW_DATATYPE_DEALLOCATE(grids)

!Start a reading loop
 ipswf=0;iaewf=0;iproj=0;iprojfit=0
 endfile=.false.
 do while (.not.endfile)
   read(funit,'(a)',err=11,end=11) readline
   line=adjustl(readline);goto 21
   11 line="";endfile=.true.
   21 continue

!  --Read core density CORE_DENSITY
   if (line(1:16)=='<ae_core_density') then
     paw_setup%ae_core_density%tread=.true.
     call paw_rdfromline(" grid",line,strg,ierr)
     if (strg == "" ) strg = "unknown"
     paw_setup%ae_core_density%grid=trim(strg)
     do ii=1,paw_setup%ngrid
       if(trim(paw_setup%ae_core_density%grid)==trim(paw_setup%radial_grid(ii)%id)) then
         mesh_size=paw_setup%radial_grid(ii)%iend-paw_setup%radial_grid(ii)%istart+1
         exit
       end if
     end do
     call paw_rdfromline(" rc",line,strg,ierr)
     if (strg /= "" ) then
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) rc(1)
       else
         read(unit=strg,fmt=*) rc(1)
       end if
     end if
     LIBPAW_ALLOCATE(paw_setup%ae_core_density%data,(mesh_size))
     !MGNAG v7[62]
     ! Runtime Error: m_pawxmlps_cpp.f90, line 1657: 
     ! Record too long for input bufferProgram terminated by I/O error on unit 9 
     ! (File="/home/buildbot/ABINIT_OD/petrus_nag/gmatteo_7.7.1-training/tests/Psps_for_tests/Al.LDA",Formatted,Sequential)
     read(funit,*) (paw_setup%ae_core_density%data(ir),ir=1,mesh_size)
     cycle
   end if

!  --Read pseudized core density CORETAIL_DENSITY
   if (line(1:20)=='<pseudo_core_density') then
     paw_setup%pseudo_core_density%tread=.true.
     call paw_rdfromline(" grid",line,strg,ierr)
     if (strg == "" ) strg = "unknown"
     paw_setup%pseudo_core_density%grid=trim(strg)
     do ii=1,paw_setup%ngrid
       if(trim(paw_setup%pseudo_core_density%grid)==trim(paw_setup%radial_grid(ii)%id)) then
         mesh_size=paw_setup%radial_grid(ii)%iend-paw_setup%radial_grid(ii)%istart+1
         exit
       end if
     end do
     call paw_rdfromline(" rc",line,strg,ierr)
     if (strg /= "" ) then
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) rc(2)
       else
         read(unit=strg,fmt=*) rc(2)
       end if
     end if
     LIBPAW_ALLOCATE(paw_setup%pseudo_core_density%data,(mesh_size))
     read(funit,*) (paw_setup%pseudo_core_density%data(ir),ir=1,mesh_size)
     cycle
   end if

!  --Read core density CORE_DENSITY
   if (line(1:31)=='<ae_core_kinetic_energy_density') then
     paw_setup%ae_core_kinetic_energy_density%tread=.true.
     call paw_rdfromline(" grid",line,strg,ierr)
     if (strg == "" ) strg = "unknown"
     paw_setup%ae_core_kinetic_energy_density%grid=trim(strg)
     do ii=1,paw_setup%ngrid
       if(trim(paw_setup%ae_core_kinetic_energy_density%grid)==trim(paw_setup%radial_grid(ii)%id)) then
         mesh_size=paw_setup%radial_grid(ii)%iend-paw_setup%radial_grid(ii)%istart+1
         exit
       end if
     end do
     call paw_rdfromline(" rc",line,strg,ierr)
     if (strg /= "" ) then
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) rc(1)
       else
         read(unit=strg,fmt=*) rc(1)
       end if
     end if
     LIBPAW_ALLOCATE(paw_setup%ae_core_kinetic_energy_density%data,(mesh_size))
     !MGNAG v7[62]
     ! Runtime Error: m_pawxmlps_cpp.f90, line 1657: 
     ! Record too long for input bufferProgram terminated by I/O error on unit 9 
     ! (File="/home/buildbot/ABINIT_OD/petrus_nag/gmatteo_7.7.1-training/tests/Psps_for_tests/Al.LDA",Formatted,Sequential)
     read(funit,*) (paw_setup%ae_core_kinetic_energy_density%data(ir),ir=1,mesh_size)
     cycle
   end if

!  --Read pseudized core density CORETAIL_DENSITY
   if (line(1:35)=='<pseudo_core_kinetic_energy_density') then
     paw_setup%pseudo_core_kinetic_energy_density%tread=.true.
     call paw_rdfromline(" grid",line,strg,ierr)
     if (strg == "" ) strg = "unknown"
     paw_setup%pseudo_core_kinetic_energy_density%grid=trim(strg)
     do ii=1,paw_setup%ngrid
       if(trim(paw_setup%pseudo_core_kinetic_energy_density%grid)==trim(paw_setup%radial_grid(ii)%id)) then
         mesh_size=paw_setup%radial_grid(ii)%iend-paw_setup%radial_grid(ii)%istart+1
         exit
       end if
     end do
     call paw_rdfromline(" rc",line,strg,ierr)
     if (strg /= "" ) then
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) rc(2)
       else
         read(unit=strg,fmt=*) rc(2)
       end if
     end if
     LIBPAW_ALLOCATE(paw_setup%pseudo_core_kinetic_energy_density%data,(mesh_size))
     read(funit,*) (paw_setup%pseudo_core_kinetic_energy_density%data(ir),ir=1,mesh_size)
     cycle
   end if

!  --Read pseudized valence density PSEUDO_VALENCE_DENSITY
   if (line(1:23)=='<pseudo_valence_density') then
     paw_setup%pseudo_valence_density%tread=.true.
     call paw_rdfromline(" grid",line,strg,ierr)
     if (strg == "" ) strg = "unknown"
     paw_setup%pseudo_valence_density%grid=trim(strg)
     do ii=1,paw_setup%ngrid
       if(trim(paw_setup%pseudo_valence_density%grid)==trim(paw_setup%radial_grid(ii)%id)) then
         mesh_size=paw_setup%radial_grid(ii)%iend-paw_setup%radial_grid(ii)%istart+1
         exit
       end if
     end do
     call paw_rdfromline(" rc",line,strg,ierr)
     if (strg /= "" ) then
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) rc(3)
       else
         read(unit=strg,fmt=*) rc(3)
       end if
     end if
     LIBPAW_ALLOCATE(paw_setup%pseudo_valence_density%data,(mesh_size))
     read(funit,*) (paw_setup%pseudo_valence_density%data(ir),ir=1,mesh_size)
     cycle
   end if

!  --Read Vbare potential VLOCFUN
   if (line(1:15)=='<zero_potential') then
     paw_setup%zero_potential%tread=.true.
     call paw_rdfromline(" grid",line,strg,ierr)
     if (strg == "" ) strg = "unknown"
     paw_setup%zero_potential%grid=trim(strg)
     do ii=1,paw_setup%ngrid
       if(trim(paw_setup%zero_potential%grid)==trim(paw_setup%radial_grid(ii)%id)) then
         mesh_size=paw_setup%radial_grid(ii)%iend-paw_setup%radial_grid(ii)%istart+1
         exit
       end if
     end do
     call paw_rdfromline(" rc",line,strg,ierr)
     if (strg /= "" ) then
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) rc(4)
       else
         read(unit=strg,fmt=*) rc(4)
       end if
     end if
     LIBPAW_ALLOCATE(paw_setup%zero_potential%data,(mesh_size))
     read(funit,*) (paw_setup%zero_potential%data(ir),ir=1,mesh_size)
     cycle
   end if

!  --Read external potential
   if (line(1:25)=='<LDA_minus_half_potential') then
     paw_setup%LDA_minus_half_potential%tread=.true.
     call paw_rdfromline(" grid",line,strg,ierr)
     if (strg == "" ) strg = "unknown"
     paw_setup%LDA_minus_half_potential%grid=trim(strg)
     do ii=1,paw_setup%ngrid
       if(trim(paw_setup%LDA_minus_half_potential%grid)==trim(paw_setup%radial_grid(ii)%id)) then
         mesh_size=paw_setup%radial_grid(ii)%iend-paw_setup%radial_grid(ii)%istart+1
         exit
       end if
     end do
     call paw_rdfromline(" rc",line,strg,ierr)
     if (strg /= "" ) then
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) rc(4)
       else
         read(unit=strg,fmt=*) rc(4)
       end if
     end if
     LIBPAW_ALLOCATE(paw_setup%LDA_minus_half_potential%data,(mesh_size))
     read(funit,*) (paw_setup%LDA_minus_half_potential%data(ir),ir=1,mesh_size)
     cycle
   end if

!  --Read Vloc for Abinit potential VLOC_ION
   if (line(1:37)=='<kresse_joubert_local_ionic_potential') then
     paw_setup%kresse_joubert_local_ionic_potential%tread=.true.
     call paw_rdfromline(" grid",line,strg,ierr)
     if (strg == "" ) strg = "unknown"
     paw_setup%kresse_joubert_local_ionic_potential%grid=trim(strg)
     do ii=1,paw_setup%ngrid
       if(trim(paw_setup%kresse_joubert_local_ionic_potential%grid)==trim(paw_setup%radial_grid(ii)%id)) then
         mesh_size=paw_setup%radial_grid(ii)%iend-paw_setup%radial_grid(ii)%istart+1
         exit
       end if
     end do
     call paw_rdfromline(" rc",line,strg,ierr)
     if (strg /= "" ) then
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) rc(5)
       else
         read(unit=strg,fmt=*) rc(5)
       end if
     end if
     LIBPAW_ALLOCATE(paw_setup%kresse_joubert_local_ionic_potential%data,(mesh_size))
     read(funit,*) (paw_setup%kresse_joubert_local_ionic_potential%data(ir),ir=1,mesh_size)
     cycle
   end if
   if (line(1:29)=='<blochl_local_ionic_potential') then
     paw_setup%blochl_local_ionic_potential%tread=.true.
     call paw_rdfromline(" grid",line,strg,ierr)
     if (strg == "" ) strg = "unknown"
     paw_setup%blochl_local_ionic_potential%grid=trim(strg)
     do ii=1,paw_setup%ngrid
       if(trim(paw_setup%blochl_local_ionic_potential%grid)==trim(paw_setup%radial_grid(ii)%id)) then
         mesh_size=paw_setup%radial_grid(ii)%iend-paw_setup%radial_grid(ii)%istart+1
         exit
       end if
     end do
     call paw_rdfromline(" rc",line,strg,ierr)
     if (strg /= "" ) then
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) rc(6)
       else
         read(unit=strg,fmt=*) rc(6)
       end if
     end if
     LIBPAW_ALLOCATE(paw_setup%blochl_local_ionic_potential%data,(mesh_size))
     read(funit,*) (paw_setup%blochl_local_ionic_potential%data(ir),ir=1,mesh_size)
     cycle
   end if

!  --Read WAVE FUNCTIONS PHI
   if (line(1:16)=='<ae_partial_wave') then
     iaewf=iaewf+1
     paw_setup%ae_partial_wave(iaewf)%tread=.true.
     call paw_rdfromline(" grid",line,strg,ierr)
     if (strg == "" ) strg = "unknown"
     paw_setup%ae_partial_wave(iaewf)%grid=trim(strg)
     call paw_rdfromline(" state",line,strg,ierr)
     paw_setup%ae_partial_wave(iaewf)%state=trim(strg)
     do ii=1,paw_setup%ngrid
       if(trim(paw_setup%ae_partial_wave(iaewf)%grid)==trim(paw_setup%radial_grid(ii)%id)) then
         mesh_size=paw_setup%radial_grid(ii)%iend-paw_setup%radial_grid(ii)%istart+1
         exit
       end if
     end do
     LIBPAW_ALLOCATE(paw_setup%ae_partial_wave(iaewf)%data,(mesh_size))
     read(funit,*) (paw_setup%ae_partial_wave(iaewf)%data(ir),ir=1,mesh_size)
     cycle
   end if

!  --Read PSEUDO WAVE FUNCTIONS TPHI
   if (line(1:20)=='<pseudo_partial_wave') then
     ipswf=ipswf+1
     paw_setup%pseudo_partial_wave(ipswf)%tread=.true.
     call paw_rdfromline(" grid",line,strg,ierr)
     if (strg == "" ) strg = "unknown"
     paw_setup%idgrid = trim(strg)
     paw_setup%pseudo_partial_wave(ipswf)%grid=trim(strg)
     call paw_rdfromline(" state",line,strg,ierr)
     paw_setup%pseudo_partial_wave(ipswf)%state=trim(strg)
     do ii=1,paw_setup%ngrid
       if(trim(paw_setup%pseudo_partial_wave(ipswf)%grid)==trim(paw_setup%radial_grid(ii)%id)) then
         mesh_size=paw_setup%radial_grid(ii)%iend-paw_setup%radial_grid(ii)%istart+1
         exit
       end if
     end do
     LIBPAW_ALLOCATE(paw_setup%pseudo_partial_wave(ipswf)%data,(mesh_size))
     read(funit,*) (paw_setup%pseudo_partial_wave(ipswf)%data(ir),ir=1,mesh_size)
     cycle
   end if

!  --Read PROJECTORS TPROJ
   if (line(1:19)=='<projector_function') then
     iproj=iproj+1
     paw_setup%projector_function(iproj)%tread=.true.
     call paw_rdfromline(" grid",line,strg,ierr)
     if (strg == "" ) strg = "unknown"
     paw_setup%projector_function(iproj)%grid=trim(strg)
     call paw_rdfromline(" state",line,strg,ierr)
     paw_setup%projector_function(iproj)%state=trim(strg)
     do ii=1,paw_setup%ngrid
       if(trim(paw_setup%projector_function(iproj)%grid)==trim(paw_setup%radial_grid(ii)%id)) then
         mesh_size=paw_setup%radial_grid(ii)%iend-paw_setup%radial_grid(ii)%istart+1
         exit
       end if
     end do
     LIBPAW_ALLOCATE(paw_setup%projector_function(iproj)%data,(mesh_size))
     read(funit,*) (paw_setup%projector_function(iproj)%data(ir),ir=1,mesh_size)
     cycle
   end if

!  --Read PROJECTORS TPROJ as gaussian representations
   if (line(1:14)=='<projector_fit') then
     if(.not.allocated(paw_setup%projector_fit)) then
        LIBPAW_DATATYPE_ALLOCATE(paw_setup%projector_fit,(paw_setup%valence_states%nval))
     end if
     iprojfit=iprojfit+1
     paw_setup%projector_fit(iprojfit)%tread=.true.
     call paw_rdfromline(" state",line,strg,ierr)
     paw_setup%projector_fit(iprojfit)%state=trim(strg)
     igauss = 0
     endgauss = .false.
     do while(.not. endgauss)
        read(funit,'(a)',err=12,end=12) readline
        line=adjustl(readline);goto 22
12      line="";endgauss=.true.
22      continue
        endgauss = (line(1:15)=='</projector_fit')
        if (line(1:9)=='<gaussian') then
           igauss = igauss + 1
           call paw_rdfromline(" factor",line,strg,ierr)
           read(strg(index(strg, '{') + 1:index(strg, ',') - 1), *) &
                & paw_setup%projector_fit(iprojfit)%factors(1, igauss)
           read(strg(index(strg, ',') + 1:index(strg, '}') - 1), *) &
                & paw_setup%projector_fit(iprojfit)%factors(2, igauss)
           call paw_rdfromline(" exponent",line,strg,ierr)
           read(strg(index(strg, '{') + 1:index(strg, ',') - 1), *) &
                & paw_setup%projector_fit(iprojfit)%expos(1, igauss)
           read(strg(index(strg, ',') + 1:index(strg, '}') - 1), *) &
                & paw_setup%projector_fit(iprojfit)%expos(2, igauss)
        end if
     end do
     paw_setup%projector_fit(iprojfit)%ngauss = igauss
     cycle
   end if

!  --Read Kinetic term KINETIC_ENERGY_MATRIX
   if (line(1:28)=='<kinetic_energy_differences>') then
     paw_setup%kinetic_energy_differences%tread=.true.
     mesh_size=paw_setup%valence_states%nval*paw_setup%valence_states%nval
     LIBPAW_ALLOCATE(paw_setup%kinetic_energy_differences%data,(mesh_size))
     read(funit,*) (paw_setup%kinetic_energy_differences%data(ir),ir=1,mesh_size)
     cycle
   end if

!  --Read Exact exchange term EXACT_EXCHANGE_X_MATRIX
   if (line(1:25)=='<exact_exchange_X_matrix>') then
     paw_setup%exact_exchange_matrix%tread=.true.
     mesh_size=paw_setup%valence_states%nval*paw_setup%valence_states%nval
     LIBPAW_ALLOCATE(paw_setup%exact_exchange_matrix%data,(mesh_size))
     read(funit,*) (paw_setup%exact_exchange_matrix%data(ir),ir=1,mesh_size)
     cycle
   end if

!  --Read Exact exchange core-core energy
   if (line(1:25)=='<exact_exchange core-core') then
     call paw_rdfromline(" core-core",line,strg,ierr)
     if (len(trim(strg))<=30) then
       strg1=trim(strg)
       read(unit=strg1,fmt=*) paw_setup%ex_cc
     else
       read(unit=strg,fmt=*) paw_setup%ex_cc
     end if
     cycle
   end if

!  --Read orthogonalisation scheme
   if (line(1:18)=='<orthogonalisation') then
     call paw_rdfromline(" scheme",line,strg,ierr)
     if (len(trim(strg))<=30) then
       strg1=trim(strg)
       read(unit=strg1,fmt=*) paw_setup%optortho
     else
       read(unit=strg,fmt=*) paw_setup%optortho
     end if
     cycle
   end if

!  --Read the Atompaw input file
   ir=0
   if ((line(1:13)=='<!-- Program:').and.(ir==1)) then
     msg=" "
     do while ((msg(1:9)/=' Program:').and.(msg(1:8)/='Program:'))
       read(funit,'(a)') msg
       write(ab_out,'(a)') trim(msg)
     end do   
     cycle
   end if 

!  End of reading loop
 end do
 if(paw_setup%rpaw<0.d0) paw_setup%rpaw=maxval(rc)
!Close the XML atomicdata file
 close(funit)

!Test flags: is anything OK ?
 found=paw_setup%atom%tread.and.paw_setup%valence_states%tread.and.&
& paw_setup%xc_functional%tread.and.paw_setup%shape_function%tread

 if (.not.paw_setup%atom%tread) then
   msg="ATOM SYMBOL not found !"
   LIBPAW_WARNING(msg)
 end if
 if (.not.paw_setup%valence_states%tread) then
   msg="VALENCE STATES not found!"
   LIBPAW_WARNING(msg)
 end if
 if (.not.paw_setup%xc_functional%tread) then
   msg="EXCHANGE/CORRELATION not found !"
   LIBPAW_WARNING(msg)
 end if
 if (.not.paw_setup%shape_function%tread) then
   msg="SHAPE FUNCTION TYPE not found !"
   LIBPAW_WARNING(msg)
 end if

 if (.not.found) then
   msg="Aborting now"
   LIBPAW_ERROR(msg)
 end if
 
 end subroutine rdpawpsxml
!!***

!-------------------------------------------------------------------------

!!****f* m_pawxmlps/rdpawpsxml_core
!! NAME
!! rdpawpsxml_core
!!
!! FUNCTION
!! Read the core wavefunctions in the XML file generated by AtomPAW
!!
!! INPUTS
!!  filename= input file name (atomicdata XML)
!!  funit= input unit number
!!
!! OUTPUT
!!  paw_setup=pseudopotential data structure
!!
!! PARENTS
!!      m_pawpsp
!!
!! CHILDREN
!!      bound_deriv,paw_rdfromline,paw_spline,paw_splint,pawrad_init
!!
!! SOURCE

 subroutine rdpawpsxml_core(energy_cor,filename,lcor,ncor,nphicor,pawrad,phi_cor)

!Arguments ---------------------------------------------
 character (len=fnlen),intent(in) :: filename
 integer,intent(out) :: nphicor
!arrays
 integer,allocatable,intent(inout) :: lcor(:),ncor(:)
 real(dp),allocatable,intent(inout) :: phi_cor(:,:),energy_cor(:)
 type(pawrad_type),intent(in) :: pawrad


!Local variables ---------------------------------------
 integer :: funit,iaewf,ii,imeshae,imsh,ir,igrid,icor,ierr,maxmeshz,mesh_size,nmesh,shft
 logical :: endfile,found,tread
 real(dp) :: yp1,ypn
 character(len=100) :: msg,version
 character (len=XML_RECL) :: line,readline
 character (len=XML_RECL) :: strg
 character (len=30) :: strg1
 integer,allocatable :: mesh_shift(:)
 real(dp),allocatable :: work(:),phitmp(:,:)
 character (len=20), allocatable :: gridwf(:),statewf(:)
 type(state_t),allocatable   :: corestate (:)
 type(radial_grid_t),allocatable   :: grids (:)
 type(pawrad_type),allocatable :: radmesh(:)

! *************************************************************************

!Open the atomicdata XML file for reading
 funit=100
 open(unit=funit,file=filename,form='formatted',status='old', recl=XML_RECL)

!Start a reading loop
 endfile=.false.
 found=.false.


 do while ((.not.endfile).and.(.not.found))
   read(funit,'(a)',err=10,end=10) readline
   line=adjustl(readline);goto 20
   10 line="";endfile=.true.
   20 continue

!  --Read VERSION
   if (line(1:10)=='<paw_setup') then
     tread=.true.
     igrid=0
     LIBPAW_DATATYPE_ALLOCATE(grids,(10))

     call paw_rdfromline(" version",line,strg,ierr)
     version=trim(strg)
     cycle
   end if

!  --Read BASIS SIZE, ORBITALS, RC AND OCCUPATIONS/STATE IDs
   if (line(1:13)=='<core_states>') then
     tread=.true.
     LIBPAW_DATATYPE_ALLOCATE(corestate,(50))
     icor=0
     do while (line(1:14)/='</core_states>')
       read(funit,'(a)') readline;line=adjustl(readline)
       if (line(1:6)=='<state') then
         icor=icor+1
         if (icor>50) then
           close(funit)
           msg="basis size too large (>50)!"
           LIBPAW_ERROR(msg)
         end if
         call paw_rdfromline(" n",line,strg,ierr)
         if (strg == "" ) then 
           corestate(icor)%nn=-1    
         else
           if (len(trim(strg))<=30) then
             strg1=trim(strg)
             read(unit=strg1,fmt=*) corestate(icor)%nn
           else
             read(unit=strg,fmt=*) corestate(icor)%nn
           end if
         end if
         call paw_rdfromline(" l",line,strg,ierr)
         if (len(trim(strg))<=30) then
           strg1=trim(strg)
           read(unit=strg1,fmt=*) corestate(icor)%ll
         else
           read(unit=strg,fmt=*) corestate(icor)%ll
         end if
         call paw_rdfromline(" f",line,strg,ierr)
         if (strg == "" ) then 
           corestate(icor)%ff=-1.d0
         else
           if (len(trim(strg))<=30) then
             strg1=trim(strg)
             read(unit=strg1,fmt=*) corestate(icor)%ff
           else
             read(unit=strg,fmt=*) corestate(icor)%ff
           end if
         end if
         call paw_rdfromline(" rc",line,strg,ierr)
         if (strg == "" ) then 
           corestate(icor)%rc=-1    
         else
           if (len(trim(strg))<=30) then
             strg1=trim(strg)
             read(unit=strg1,fmt=*) corestate(icor)%rc
           else
             read(unit=strg,fmt=*) corestate(icor)%rc
           end if
         end if
         call paw_rdfromline(" e",line,strg,ierr)
         if (len(trim(strg))<=30) then
           strg1=trim(strg)
           read(unit=strg1,fmt=*) corestate(icor)%ee
         else
           read(unit=strg,fmt=*) corestate(icor)%ee
         end if
         call paw_rdfromline(" id",line,strg,ierr)
         corestate(icor)%id = trim(strg)
       end if
     end do
     cycle
   end if

!  --Read MESH_STEP AND NUMBER OF POINTS
   if (line(1:12)=='<radial_grid')then
     igrid=igrid+1
     call paw_rdfromline(" eq",line,strg,ierr)
     grids(igrid)%eq = trim(strg)
     call paw_rdfromline(" a",line,strg,ierr)
     if (strg == "" ) then
       grids(igrid)%aa=0.d0
     else
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) grids(igrid)%aa
       else
         read(unit=strg,fmt=*) grids(igrid)%aa
       end if
     end if
     call paw_rdfromline(" n",line,strg,ierr)
     if (strg == "" ) then
       grids(igrid)%nn=0
     else
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) grids(igrid)%nn
       else
         read(unit=strg,fmt=*) grids(igrid)%nn
       end if
     end if
     call paw_rdfromline(" d",line,strg,ierr)
     if (strg == "" ) then
       grids(igrid)%dd=0.d0
     else
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) grids(igrid)%dd
       else
         read(unit=strg,fmt=*) grids(igrid)%dd
       end if
     end if
     call paw_rdfromline(" b",line,strg,ierr)
     if (strg == "" ) then
       grids(igrid)%bb=0.d0
     else
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) grids(igrid)%bb
       else
         read(unit=strg,fmt=*) grids(igrid)%bb
       end if
     end if
     call paw_rdfromline("istart",line,strg,ierr)
     if (len(trim(strg))<=30) then
       strg1=trim(strg)
       read(unit=strg1,fmt=*) grids(igrid)%istart
     else
       read(unit=strg,fmt=*) grids(igrid)%istart
     end if
     call paw_rdfromline("iend",line,strg,ierr)
     if (len(trim(strg))<=30) then
       strg1=trim(strg)
       read(unit=strg1,fmt=*) grids(igrid)%iend
     else
       read(unit=strg,fmt=*) grids(igrid)%iend
     end if
     call paw_rdfromline(" id",line,strg,ierr)
     grids(igrid)%id = trim(strg)
     if(igrid>10)then
       close(funit) 
       msg="igrid>10"
       LIBPAW_ERROR(msg)
     end if
     found=.true.
     cycle
   end if

!  End of reading loop
 end do

 nphicor=icor
 nmesh=igrid
 if(nmesh>0)then
   LIBPAW_DATATYPE_ALLOCATE(radmesh,(nmesh))
   LIBPAW_ALLOCATE(mesh_shift,(nmesh))
   do imsh=1,nmesh
     radmesh(imsh)%mesh_type=-1
     radmesh(imsh)%rstep=zero
     radmesh(imsh)%lstep=zero
     mesh_shift(imsh)=0
     select case(trim(grids(imsh)%eq))
       case("r=a*exp(d*i)")
         mesh_shift(imsh)=1
         radmesh(imsh)%mesh_type=3
         radmesh(imsh)%mesh_size=grids(imsh)%iend-grids(imsh)%istart+1+mesh_shift(imsh)
         radmesh(imsh)%rstep=grids(imsh)%aa
         radmesh(imsh)%lstep=grids(imsh)%dd
       case("r=a*i/(1-b*i)")
         write(msg, '(3a)' )&
&         '  the grid r=a*i/(1-b*i) is not implemented in ABINIT !',ch10,&
&         '  Action: check your psp file.'
         LIBPAW_ERROR(msg)
       case("r=a*i/(n-i)")
         mesh_shift(imsh)=0
         radmesh(imsh)%mesh_type=5
         radmesh(imsh)%mesh_size=grids(imsh)%iend-grids(imsh)%istart+1+mesh_shift(imsh)
         radmesh(imsh)%rstep=grids(imsh)%aa
         radmesh(imsh)%lstep=dble(grids(imsh)%nn)
       case("r=a*(exp(d*i)-1)")
         mesh_shift(imsh)=0
         radmesh(imsh)%mesh_type=2
         radmesh(imsh)%mesh_size=grids(imsh)%iend-grids(imsh)%istart+1+mesh_shift(imsh)
         if(grids(imsh)%istart==1)radmesh(imsh)%mesh_size=radmesh(imsh)%mesh_size+1
         radmesh(imsh)%rstep=grids(imsh)%aa
         radmesh(imsh)%lstep=grids(imsh)%dd
       case("r=d*i")
         mesh_shift(imsh)=0
         radmesh(imsh)%mesh_type=1
         radmesh(imsh)%mesh_size=grids(imsh)%iend-grids(imsh)%istart+1+mesh_shift(imsh)
         if(grids(imsh)%istart==1)radmesh(imsh)%mesh_size=radmesh(imsh)%mesh_size+1
         radmesh(imsh)%rstep=grids(imsh)%dd
       case("r=(i/n+a)^5/a-a^4")
         write(msg, '(3a)' )&
&       '  the grid r=(i/n+a)^5/a-a^4 is not implemented in ABINIT !',ch10,&
&       '  Action: check your psp file.'
         LIBPAW_ERROR(msg)
     end select
   end do
 end if

!Initialize radial meshes
 do imsh=1,nmesh
   call pawrad_init(radmesh(imsh))
 end do

 maxmeshz=maxval(radmesh(:)%mesh_size)
 LIBPAW_DATATYPE_ALLOCATE(gridwf,(nphicor))
 LIBPAW_DATATYPE_ALLOCATE(statewf,(nphicor))
 LIBPAW_ALLOCATE(phitmp,(maxmeshz,nphicor))
 phitmp(:,:)=zero

!Start of reading loop
 iaewf=0 ; endfile=.false.
 do while (.not.endfile)
   read(funit,'(a)',err=11,end=11) readline
   line=adjustl(readline);goto 21
   11 line="";endfile=.true.
   21 continue

!  --Read CORE WAVE FUNCTIONS PHI
   if (line(1:21)=='<ae_core_wavefunction') then
     iaewf=iaewf+1
     tread=.true.
     call paw_rdfromline(" grid",line,strg,ierr)
     if (strg == "" ) strg = "unknown"
     gridwf(iaewf)=trim(strg)
     call paw_rdfromline(" state",line,strg,ierr)
     statewf(iaewf)=trim(strg)
     do ii=1,nmesh
       if(trim(gridwf(iaewf))==trim(grids(ii)%id)) then
         mesh_size=grids(ii)%iend-grids(ii)%istart+1
         exit
       end if
     end do
     read(funit,*) (phitmp(ir,iaewf),ir=1,mesh_size)
     cycle
   end if
!  End of reading loop
 end do

 if(nphicor>0)then
   LIBPAW_ALLOCATE(ncor,(nphicor))
   LIBPAW_ALLOCATE(lcor,(nphicor))
   LIBPAW_ALLOCATE(energy_cor,(nphicor))
   LIBPAW_ALLOCATE(phi_cor,(pawrad%mesh_size,nphicor))
   phi_cor(:,:)=zero
   do ii=1,nphicor
     ncor(ii)=corestate(ii)%nn
     lcor(ii)=corestate(ii)%ll
     energy_cor(ii)=corestate(ii)%ee
     do imsh=1,nmesh
       if(trim(gridwf(ii))==trim(grids(imsh)%id)) imeshae=imsh
     end do
     if ((pawrad%mesh_type/=radmesh(imeshae)%mesh_type) &
&    .or.(pawrad%rstep/=radmesh(imeshae)%rstep) &
&    .or.(pawrad%lstep/=radmesh(imeshae)%lstep)) then
       mesh_size=radmesh(imeshae)%mesh_size
       LIBPAW_ALLOCATE(work,(mesh_size))
       call bound_deriv(phitmp(1:mesh_size,ii),radmesh(imeshae),mesh_size,yp1,ypn)
       call paw_spline(radmesh(imeshae)%rad(1:mesh_size),phitmp(1:mesh_size,ii),mesh_size,yp1,ypn,work(1:mesh_size))
       ir=pawrad%mesh_size
       if (radmesh(imeshae)%rmax<pawrad%rmax+tol8) ir=pawrad_ifromr(pawrad,radmesh(imeshae)%rmax)-1
       call paw_splint(mesh_size,radmesh(imeshae)%rad(1:mesh_size),phitmp(1:mesh_size,ii),work(1:mesh_size),&
&                      ir,pawrad%rad(1:ir),phi_cor(1:ir,ii))
       phi_cor(1:ir,ii)=phi_cor(1:ir,ii)*pawrad%rad(1:ir)
       LIBPAW_DEALLOCATE(work)
     else
       shft=mesh_shift(imeshae)
       mesh_size=min(radmesh(imeshae)%mesh_size,pawrad%mesh_size)
       phi_cor(1+shft:mesh_size,ii)=phitmp(1:mesh_size-shft,ii)*radmesh(imeshae)%rad(1:mesh_size-shft)
       if (shft==1) phi_cor(1,ii)=zero
     end if
   end do
 end if

 LIBPAW_DATATYPE_DEALLOCATE(radmesh)
 LIBPAW_DATATYPE_DEALLOCATE(mesh_shift)

 LIBPAW_DATATYPE_DEALLOCATE(grids)
 LIBPAW_DATATYPE_DEALLOCATE(corestate)

 LIBPAW_DATATYPE_DEALLOCATE(gridwf)
 LIBPAW_DATATYPE_DEALLOCATE(statewf)
 LIBPAW_DEALLOCATE(phitmp)

!Close the XML atomicdata file
 close(funit)


 end subroutine rdpawpsxml_core
!!***

!-------------------------------------------------------------------------

end module m_pawxmlps
!!***
