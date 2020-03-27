!!****m* ABINIT/m_outxml
!! NAME
!! m_outxml
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt.
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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

module m_outxml

 use defs_basis
 use m_abicore
 use m_errors
 use m_dtset

 use m_io_tools,    only : open_file
 use m_geometry,    only : xcart2xred, xred2xcart
 use m_results_gs , only : results_gs_type

 implicit none

 private

 public :: outxml_open
 public :: outxml_finalise
 public :: out_resultsgs_XML
 public :: out_geometry_XML
!!***

contains

!!****f* m_outxml/outxml_open
!! NAME
!! outxml_open
!!
!! FUNCTION
!! Open the XML log file, and write the header inside.
!! (see extras/post_processing/abinitRun.dtd)
!! Warning: this method is not thread safe and should be called only by one thread.
!!
!! INPUTS
!!  filename=the name of the file to write to.
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      xred2xcart
!!
!! SOURCE

subroutine outxml_open(filename)

!Arguments -------------------------------
  character(len = *), intent(in) :: filename
!Local variables -------------------------
 character(len=500) :: msg

! *********************************************************************

!ABINIT has been compiled with XML output support, then we open the
!channel of the XML output file.
 if (open_file(trim(filename)//"_LOG.xml", msg, unit=ab_xml_out, form="formatted", action="write") /= 0) then
   MSG_ERROR(msg)
 end if

 write(ab_xml_out, "(A)") '<?xml version="1.0" encoding="utf-8"?>'
!DTD declaration : root element is "abinitRun", and DTD file is not public
!but given in the distribution by the file abinitRun.dtd.
 write(ab_xml_out, "(A)") '<!DOCTYPE abinitRun SYSTEM "extras/post_processing/abinitRun.dtd">'
!Creating the root markup.
 write(ab_xml_out, "(A)") '<abinitRun>'

end subroutine outxml_open
!!***

!!****f* m_outxml/outxml_finalise
!!
!! NAME
!! outxml_finalise
!!
!! FUNCTION
!! Close the XML log file, and write the timing information.
!! (see extras/post_processing/abinitRun.dtd)
!! Warning : this method is not thread safe and should be called
!! only by one thread.
!!
!! INPUTS
!!  tsec=the cpu time and the wall time in seconds.
!!  values=the date values returned by date_and_time() intrinsic Fortran routine.
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      xred2xcart
!!
!! SOURCE

subroutine outxml_finalise(tsec, values)

!Arguments -------------------------------
  integer, intent(in) :: values(8)
  real(dp), intent(in) :: tsec(2)
!Local variables -------------------------
  character(len=500) :: message

! *********************************************************************

 write(ab_xml_out, "(A)") '  <!-- timeInfo markup : cpu and wall attributes are given in seconds. -->'
 write(ab_xml_out, "(A,I0,A,I0,A,I0,A)", advance = "NO") &
& '  <timeInfo day="', values(3) , &
& '" month="', values(2) ,'" year="', values(1) ,'" '
 write(message, *) tsec(1)
 write(ab_xml_out, "(A,A,A)", advance = "NO") 'cpu="', trim(message) ,'"'
 write(message, *) tsec(2)
 write(ab_xml_out, "(A,A,A)") ' wall="', trim(message) ,'" />'
 write(ab_xml_out, "(A)") "</abinitRun>"
!ABINIT has been compiled with XML output support, then we close the channel of the XML output file.
 close(unit = ab_xml_out)

end subroutine outxml_finalise
!!***

!!****f* m_outxml/out_resultsgs_xml
!!
!! NAME
!! out_resultsgs_xml
!!
!! FUNCTION
!! Output in the XML file, the decomposition of the energy and
!! the forces after a scfcv loop.
!! (see extras/post_processing/abinitRun.dtd)
!! Warning : this method is not thread safe and should be called only by one thread.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  level=use for indentation of the XML, two spaces are added by level.
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!   forces and its components, the stress tensor) of a ground-state computation.
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      xred2xcart
!!
!! SOURCE

subroutine out_resultsgs_XML(dtset, level, results_gs, usepaw)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: level,usepaw
 type(dataset_type),intent(in) :: dtset
 type(results_gs_type),intent(inout) :: results_gs

!Local variables -------------------------
  character(len = 1), parameter :: axes_names(3) = (/ "x", "y", "z" /)
!scalars
 integer :: i,j
 character(len=128) :: spacer,value

! *************************************************************************

!Compute the spacer to put before each markup
 write(spacer, "(A,I0)") "A", 2 * level

!Begin with the energy part
 write(ab_xml_out, "("//trim(spacer)//",A)", advance = "NO") " ", '<energy'
 if (dtset%iscf < 10) then
   write(ab_xml_out, "(A)", advance = "NO") ' type="direct"'
   write(value, "(es20.8)") results_gs%energies%e_kinetic
   write(ab_xml_out, "(A,A,A)", advance = "NO") ' kinetic="', trim(value) ,'"'
   write(value, "(es20.8)") results_gs%energies%e_localpsp
   write(ab_xml_out, "(A,A,A)", advance = "NO") ' local="', trim(value) ,'"'
   write(value, "(es20.8)") results_gs%energies%e_nlpsp_vfock
   write(ab_xml_out, "(A,A,A)", advance = "NO") ' non-local="', trim(value) ,'"'
   if (usepaw == 1) then
     write(value, "(es20.8)") results_gs%energies%e_paw
     write(ab_xml_out, "(A,A,A)", advance = "NO") ' paw="', trim(value) ,'"'
   end if
 else
   write(ab_xml_out, "(A)", advance = "NO") ' type="double-counting"'
   write(value, "(es20.8)") results_gs%energies%e_eigenvalues
   write(ab_xml_out, "(A,A,A)", advance = "NO") ' eigen-values="', trim(value) ,'"'
   write(value, "(es20.8)") results_gs%energies%e_xcdc
   write(ab_xml_out, "(A,A,A)", advance = "NO") ' xcdc="', trim(value) ,'"'
   if (usepaw == 1) then
     write(value, "(es20.8)") results_gs%energies%e_pawdc
     write(ab_xml_out, "(A,A,A)", advance = "NO") ' pawdc="', trim(value) ,'"'
   end if
 end if
 if (dtset%berryopt == 4 .or. dtset%berryopt == 6 .or. dtset%berryopt == 7 .or. &
& dtset%berryopt ==14 .or. dtset%berryopt ==16 .or. dtset%berryopt ==17) then
   write(value, "(es20.8)") results_gs%energies%e_elecfield
   write(ab_xml_out, "(A,A,A)", advance = "NO") ' electric-field="', trim(value) ,'"'
 end if
 if(dtset%occopt >= 3 .and. dtset%occopt <= 8) then
   write(value, "(es20.8)") results_gs%energies%e_entropy
   write(ab_xml_out, "(A,A,A)", advance = "NO") ' entropy="', trim(value) ,'"'
 end if
 write(value, "(es20.8)") results_gs%energies%e_ewald
 if (dtset%icoulomb == 0) then
   write(ab_xml_out, "(A,A,A)", advance = "NO") ' ewald="', trim(value) ,'"'
 else
   write(ab_xml_out, "(A,A,A)", advance = "NO") ' ion-ion="', trim(value) ,'"'
 end if
 write(value, "(es20.8)") results_gs%energies%e_chempot
 write(ab_xml_out, "(A,A,A)", advance = "NO") ' chempot="', trim(value) ,'"'
 write(value, "(es20.8)") results_gs%energies%e_hartree
 write(ab_xml_out, "(A,A,A)", advance = "NO") ' hartree="', trim(value) ,'"'
 write(value, "(es20.8)") results_gs%energies%e_corepsp
 write(ab_xml_out, "(A,A,A)", advance = "NO") ' core="', trim(value) ,'"'
 write(value, "(es20.8)") results_gs%energies%e_xc
 write(ab_xml_out, "(A,A,A)", advance = "NO") ' xc="', trim(value) ,'"'
 write(value, "(es20.8)") results_gs%etotal
 write(ab_xml_out, "(A,A,A)", advance = "NO") ' total="', trim(value) ,'"'
 write(ab_xml_out, "(A)") ' />'


!finish with the forces part
 if (dtset%optforces == 1) then
   write(ab_xml_out, "("//trim(spacer)//",A)") " ", '<forces>'
   do i = 1, dtset%natom, 1
     write(ab_xml_out, "("//trim(spacer)//",A)", advance = "NO") " ", '  <force'
     write(ab_xml_out, "(A,I0,A,I0,A)", advance = "NO") ' atom="a_', dtset%jdtset ,'_', i ,'"'
     do j = 1, 3, 1
       write(value, "(es20.8)") results_gs%fcart(j, i)
       write(ab_xml_out, "(A,A,A,A,A)", advance = "NO") ' ', axes_names(j), '="', trim(value) ,'"'
     end do
     write(ab_xml_out, "(A)") ' />'
   end do
   write(ab_xml_out, "("//trim(spacer)//",A)") " ", '</forces>'
 end if

end subroutine out_resultsgs_XML
!!***

!!****f* m_outxml/out_geometry_xml
!!
!! NAME
!! out_geometry_xml
!!
!! FUNCTION
!! Output in the XML file, the box size and the atomic position.
!! (see extras/post_processing/abinitRun.dtd)
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  level=use for indentation of the XML, two spaces are added by level.
!!  natom=number of atoms.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      xred2xcart
!!
!! SOURCE

subroutine out_geometry_XML(dtset, level, natom, rprimd, xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: level,natom
 type(dataset_type),intent(in) :: dtset
!arrays
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: xred(3,natom)

!Local variables -------------------------
  character(len = 1), parameter :: rprimd_names(3) = (/ "x", "y", "z" /)
!scalars
 integer :: i,j
 character(len=128) :: spacer,value
!arrays
 real(dp),allocatable :: xcart(:,:)

! *************************************************************************

!Compute the spacer to put before each markup
 write(spacer, "(A,I0)") "A", 2 * level
 write(ab_xml_out, "("//trim(spacer)//",A)") " ", '<geometry>'
!Compute the cartesian coordinates of atoms
 ABI_ALLOCATE(xcart,(3, natom))
 call xred2xcart(natom, rprimd, xcart, xred)
!Ouput the rprimd matrix
 write(ab_xml_out, "("//trim(spacer)//",A)", advance = "NO") " ", '  <rprimd'
 do i = 1, 3, 1
   do j = 1, 3, 1
     write(value, "(es20.8)") rprimd(i, j)
     write(ab_xml_out, "(A,A,I0,A,A,A)", advance = "NO") ' ', rprimd_names(i), j, '="', trim(value) ,'"'
   end do
 end do
 write(ab_xml_out, "(A)") ' />'
!Output the atom position
 do i = 1, natom, 1
   write(ab_xml_out, "("//trim(spacer)//",A)", advance = "NO") " ", '  <position'
   write(ab_xml_out, "(A,I0,A,I0,A)", advance = "NO") ' atom="a_', dtset%jdtset ,'_', i ,'"'
   do j = 1, 3, 1
     write(value, "(es20.8)") xcart(j, i)
     write(ab_xml_out, "(A,A,A,A,A)", advance = "NO") ' ', rprimd_names(j), '="', trim(value) ,'"'
   end do
   write(ab_xml_out, "(A)") ' />'
 end do
 ABI_DEALLOCATE(xcart)
 write(ab_xml_out, "("//trim(spacer)//",A)") " ", '</geometry>'

end subroutine out_geometry_XML
!!***

end module m_outxml
!!***
