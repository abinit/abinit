!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_phonons
!! NAME
!! m_phonons
!!
!! FUNCTION
!! Module for the phonon density of states.
!! Container type is defined, and destruction, print subroutines
!! as well as the central mkphdos
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group (XG,MG,MJV)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_phonons

 use defs_basis
 use m_errors
 use m_xmpi
 use m_profiling_abi
 use m_tetrahedron
 use m_nctk
 use iso_c_binding
 use m_crystal_io
 use m_atprj
 use m_sortph
 use m_ddb
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_fstrings,        only : itoa, ftoa, sjoin, ktoa, strcat, basename, replace
 use m_numeric_tools,   only : simpson_int, wrap2_pmhalf
 use m_time,            only : cwtime
 use m_io_tools,        only : open_file
 use defs_abitypes,     only : dataset_type
 use m_dynmat,          only : gtdyn9, dfpt_phfrq
 use m_crystal,         only : crystal_t
 use m_bz_mesh,         only : isamek, make_path, kpath_t, kpath_new, kpath_free
 use m_ifc,             only : ifc_type, ifc_fourq, ifc_calcnwrite_nana_terms
 use m_anaddb_dataset,  only : anaddb_dataset_type
 use m_kpts,            only : kpts_ibz_from_kptrlatt
 use m_special_funcs,   only : bose_einstein
 use m_sort,            only : sort_dp
 use m_supercell

 implicit none

 private

! TODO Write object to store the bands
 public :: mkphbs                        ! Compute phonon band structure
 public :: phonons_write_xmgrace         ! Write phonons bands in Xmgrace format.
 public :: phonons_write_gnuplot         ! Write phonons bands in gnuplot format.
 public :: ifc_mkphbs                    ! Compute the phonon band structure from the IFC and write data to file(s)

 public :: thermal_supercell_make
 public :: thermal_supercell_free
 public :: thermal_supercell_print
!!***

!!****t* m_phonons/phonon_dos_type
!! NAME
!! phonon_dos_type
!!
!! FUNCTION
!! Container for phonon DOS and atom projected contributions
!!
!! SOURCE

 type,public :: phonon_dos_type

! Integer
  integer :: ntypat
  ! Number of type of atoms.

  integer :: natom
  ! Number of atoms is the unit cell.

  integer :: prtdos
  ! Option of DOS calculation (1 for Gaussian, 2 for tetrahedrons).

  integer :: nomega
  ! Number of frequency points in DOS mesh.

  integer :: nqibz
  ! Number of q-points in the IBZ.

! Reals
  real(dp) :: omega_min
  ! Min frequency for DOS calculation.

  real(dp) :: omega_max
  ! Max frequency for DOS calculation.

  real(dp) :: omega_step
  ! Frequency step.

  real(dp) :: dossmear
  ! Gaussian broadening.

! Real pointers
  real(dp),allocatable :: atom_mass(:)
   ! atom_mass(natom)

  real(dp),allocatable :: omega(:)
   ! omega(nomega)
   ! Frequency grid.

  real(dp),allocatable :: phdos(:)
   ! phdos(nomega)
   ! phonon DOS.

  real(dp),allocatable :: phdos_int(:)
   ! phdos_int(nomega)
   ! integrated phonon DOS

  real(dp),allocatable :: pjdos(:,:,:)
   ! pjdos(nomega,3,natom)
   ! projected DOS (over atoms and cartesian directions)

  real(dp),allocatable :: pjdos_int(:,:,:)
   ! pjdos_int(nomega,3,natom)
   ! Integrated atomic PJDOS along the three cartesian directions.

  real(dp),allocatable :: pjdos_type(:,:)
   ! pjdos_type(nomega,ntypat)
   ! phonon DOS contribution arising from a particular atom-type.

  real(dp),allocatable :: pjdos_type_int(:,:)
   ! pjdos_type_int(nomega,ntypat)
   ! Integrate phonon DOS contribution arising from a particular atom-type.

  real(dp),allocatable :: pjdos_rc_type(:,:,:)
   ! phdos(nomega,3,ntypat)
   ! phonon DOS contribution arising from a particular atom-type
   ! decomposed along the three cartesian directions.

  real(dp),allocatable :: msqd_dos_atom(:,:,:,:)
   ! msqd_dos_atom(nomega,3,3,natom)
   ! mean square displacement matrix, frequency dependent like a DOS, tensor in cartesian coords.
   ! allows to calculate Debye Waller factors by integration with 1/omega
   ! and the Bose Einstein factor

 end type phonon_dos_type

 public :: mkphdos
 public :: phdos_print
 public :: phdos_print_debye
 public :: phdos_print_msqd
 public :: phdos_print_thermo
 public :: phdos_free
 public :: phdos_ncwrite
!!**

CONTAINS  !===============================================================================
!!***

!!****f* m_phonons/phdos_print
!!
!! NAME
!! phdos_print
!!
!! FUNCTION
!! Print out phonon DOS (and partial DOS etc) in meV units
!!
!! INPUTS
!! PHdos= container object for phonon DOS
!! fname=File name for output
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      anaddb,eph,m_tdep_phdos
!!
!! CHILDREN
!!      ifc_fourq,kpath_free,phonons_ncwrite,phonons_write_gnuplot
!!      phonons_write_phfrq,phonons_write_xmgrace,xmpi_sum_master
!!
!! SOURCE

subroutine phdos_print(PHdos,fname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phdos_print'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=*),intent(in) :: fname
 type(phonon_dos_type),intent(in) :: PHdos

!Local variables-------------------------------
 integer :: io,itype,unt,unt_by_atom,unt_msqd,iatom
 real(dp) :: tens(3,3)
 character(len=500) :: msg
 character(len=500) :: msg_method
 character(len=fnlen) :: fname_by_atom
 character(len=fnlen) :: fname_msqd
 character(len=3) :: unitname

! *************************************************************************

! Use Ha units everywhere
 unitname='Ha'

 select case (PHdos%prtdos)
 case (1)
   write(msg_method,'(a,es16.8,2a,i0)')&
&   '# Gaussian method with smearing = ',PHdos%dossmear,unitname,', nqibz =',PHdos%nqibz
 case (2)
   write(msg_method,'(a,i0)')'# Tetrahedron method, nqibz= ',PHdos%nqibz
 case default
   MSG_ERROR(sjoin(" Wrong prtdos: ",itoa(PHdos%prtdos)))
 end select

! === Open external file and write results ===
! TODO Here I have to rationalize how to write all this stuff!!
!
 if (open_file(fname,msg,newunit=unt,form="formatted",action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 write(msg,'(3a)')'# ',ch10,'# Phonon density of states and atom type projected DOS'
 call wrtout(unt,msg,'COLL')
 write(msg,'(6a)')'# ',ch10,'# Energy in ',unitname,', DOS in states/',unitname
 call wrtout(unt,msg,'COLL')
 call wrtout(unt,msg_method,'COLL')
 write(msg,'(5a)')'# ',ch10,'# omega     PHDOS    INT_PHDOS   PJDOS[atom_type=1]  INT_PJDOS[atom_type=1] ...  ',ch10,'# '
 call wrtout(unt,msg,'COLL')
 do io=1,PHdos%nomega
   write(unt,'(3es17.8)',advance='NO')PHdos%omega(io),PHdos%phdos(io),PHdos%phdos_int(io)
   do itype=1,PHdos%ntypat
     write(unt,'(2es17.8,2x)',advance='NO')PHdos%pjdos_type(io,itype),PHdos%pjdos_type_int(io,itype)
   end do
   write(unt,*)
 end do
 close(unt)



 fname_by_atom = trim(fname) // "_by_atom"
 if (open_file(fname_by_atom,msg,newunit=unt_by_atom,form="formatted",action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 write(msg,'(3a)')'# ',ch10,'# Phonon density of states and atom projected DOS'
 call wrtout(unt_by_atom,msg,'COLL')
 write(msg,'(6a)')'# ',ch10,'# Energy in ',unitname,', DOS in states/',unitname
 call wrtout(unt_by_atom,msg,'COLL')
 call wrtout(unt_by_atom,msg_method,'COLL')
 write(msg,'(5a)')'# ',ch10,'# omega     PHDOS    PJDOS[atom=1]  PJDOS[atom=2] ...  ',ch10,'# '
 call wrtout(unt_by_atom,msg,'COLL')
 do io=1,PHdos%nomega
   write(unt_by_atom,'(2es17.8)',advance='NO')PHdos%omega(io),PHdos%phdos(io)
   do iatom=1,PHdos%natom
     write(unt_by_atom,'(1es17.8,2x)',advance='NO') sum(PHdos%pjdos(io,1:3,iatom))
   end do
   write(unt_by_atom,*)
 end do
 close(unt_by_atom)



 fname_msqd = trim(fname) // "_msqd"
 if (open_file(fname_msqd,msg,newunit=unt_msqd,form="formatted",action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 write(msg,'(3a)')'# ',ch10,'# Phonon density of states weighted msq displacement matrix (set to zero below 1e-12)'
 call wrtout(unt_msqd,msg,'COLL')
 write(msg,'(6a)')'# ',ch10,'# Energy in ',unitname,', DOS in bohr^2 states/',unitname
 call wrtout(unt_msqd,msg,'COLL')
 call wrtout(unt_msqd,msg_method,'COLL')
 write(msg,'(5a)')'# ',ch10,'# omega     MSQDisp[atom=1, xx, yy, zz, yz, xz, xy]  MSQDisp[atom=2, xx, yy,...] ...  ',ch10,'# '
 call wrtout(unt_msqd,msg,'COLL')
 do io=1,PHdos%nomega
   write(unt_msqd,'(2es17.8)',advance='NO')PHdos%omega(io)
   do iatom=1,PHdos%natom
     tens = PHdos%msqd_dos_atom(io,:,:,iatom)
     where (abs(tens) < tol12)
        tens = zero
     end where
     write(unt_msqd,'(6es17.8,2x)',advance='NO') &
&        tens(1,1), &
&        tens(2,2), &
&        tens(3,3), &
&        tens(2,3), &
&        tens(1,3), &
&        tens(1,2)
   end do
   write(unt_msqd,*)
 end do
 close(unt_msqd)


end subroutine phdos_print
!!***

!----------------------------------------------------------------------

!****f* m_phonons/phdos_print_debye
!!
!! NAME
!! phdos_print_debye
!!
!! FUNCTION
!! Print out global Debye temperature, force constant, etc... from phonon DOS
!!
!! INPUTS
!! phonon_dos= container object for phonon DOS
!! ucvol = unit cell volume
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      ifc_fourq,kpath_free,phonons_ncwrite,phonons_write_gnuplot
!!      phonons_write_phfrq,phonons_write_xmgrace,xmpi_sum_master
!!
!! SOURCE

subroutine phdos_print_debye(PHdos, ucvol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phdos_print_debye'
 use interfaces_14_hidewrite
!End of the abilint section

implicit none

!Arguments ------------------------------------
 real(dp), intent(in) :: ucvol
 type(phonon_dos_type),intent(in) :: PHdos

!Local variables-------------------------------
 integer :: io, iomax, iomin
 real(dp) :: avgom2dos, avgspeedofsound
 real(dp) :: debyefreq, meanfreq, meanfreq2
 character(len=500) :: msg
!arrays
 real(dp), allocatable :: om2dos(:), om1dos(:), intdos(:)

! *************************************************************************

! average speed of sound: coefficient of omega^2 in the DOS is = Volume / 2 pi^2 hbar^3 v_s^3
! first find how far out we can fit with a parabola
 ABI_MALLOC(om2dos,(PHdos%nomega))
 ABI_MALLOC(om1dos,(PHdos%nomega))
 ABI_MALLOC(intdos,(PHdos%nomega))
 avgom2dos = zero
 om1dos = zero
 om2dos = zero
 do io=1,PHdos%nomega
   if (abs(PHdos%omega(io)) > 1.e-8) then
     om1dos(io) = PHdos%phdos(io) / PHdos%omega(io)
     om2dos(io) = PHdos%phdos(io) / PHdos%omega(io)**2
   end if
 end do

! integrate dos / omega
 intdos = zero
 call simpson_int(PHdos%nomega,PHdos%omega_step,om1dos,intdos)
 meanfreq = intdos(PHdos%nomega)

! integrate dos / omega^2
 intdos = zero
 call simpson_int(PHdos%nomega,PHdos%omega_step,om2dos,intdos)
 meanfreq2 = intdos(PHdos%nomega)

 iomin = 1; iomax = PHdos%nomega
 do io = 1, PHdos%nomega
   ! skip eventual negative frequency modes
   if (PHdos%omega(io) <= tol10) then
     iomin = io
     cycle
   end if

   ! accumulate dos * om^2 to make an average
   avgom2dos = avgom2dos + om2dos(io)
   ! first deviation from initial value of more than 10 percent
   if (abs(one-om2dos(iomin)/om2dos(io)) > 0.1_dp) then
     iomax = io
     exit
   end if
 end do

 avgom2dos = avgom2dos / (iomax-iomin+1)
! this value is also useful for partial atomic DOS, related to kinetic energy and Force constant in Moessbauer

 avgspeedofsound = (ucvol / 2 / pi**2 / avgom2dos)**third
 write (msg,'(a,E20.10,3a,F16.4,2a)') ' Average speed of sound: ', avgspeedofsound, ' (at units) ',ch10,&
&             '-                      = ', avgspeedofsound * Bohr_Ang * 1.d-13 / Time_Sec, ' [km/s]',ch10
 call wrtout (ab_out,msg,"COLL")
 call wrtout (std_out,msg,"COLL")

! Debye frequency = vs * (6 pi^2 natom / ucvol)**1/3
 debyefreq = avgspeedofsound * (six*pi**2/ucvol)**(1./3.)
 write (msg,'(a,E20.10,3a,E20.10,a)') ' Debye frequency from DOS: ', debyefreq, ' (Ha) ',ch10,&
&                                    '-                        = ', debyefreq*Ha_THz, ' (THz)'
 call wrtout (ab_out,msg,"COLL")
 call wrtout (std_out,msg,"COLL")

! Debye temperature = hbar * Debye frequency / kb
 write (msg,'(a,E20.10,2a)') '-Debye temperature from DOS: ', debyefreq*Ha_K, ' (K)', ch10
 call wrtout (ab_out,msg,"COLL")
 call wrtout (std_out,msg,"COLL")

! average force constant
 ABI_FREE(om2dos)
 ABI_FREE(om1dos)
 ABI_FREE(intdos)

end subroutine phdos_print_debye
!!***

!----------------------------------------------------------------------

!****f* m_phonons/phdos_print_thermo
!!
!! NAME
!! phdos_print_thermo
!!
!! FUNCTION
!! Print out global thermodynamic quantities based on DOS
!!
!! INPUTS
!! phonon_dos= container object for phonon DOS
!! ucvol = unit cell volume
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      anaddb,eph
!!
!! CHILDREN
!!      ifc_fourq,kpath_free,phonons_ncwrite,phonons_write_gnuplot
!!      phonons_write_phfrq,phonons_write_xmgrace,xmpi_sum_master
!!
!! SOURCE

subroutine phdos_print_thermo(PHdos, fname, ntemper, tempermin, temperinc)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phdos_print_thermo'
 use interfaces_14_hidewrite
!End of the abilint section

implicit none

!Arguments ------------------------------------
 integer, intent(in) :: ntemper
 real(dp), intent(in) :: tempermin, temperinc
 type(phonon_dos_type),intent(in) :: PHdos
 character(len=*),intent(in) :: fname

!Local variables-------------------------------
 integer :: iomega, itemper, thermal_unit
 character(len=500) :: msg
 real(dp) :: wover2t, ln2shx, cothx, invsinh2
 real(dp) :: tmp, domega
!arrays
 real(dp), allocatable :: free(:), energy(:), entropy(:), spheat(:),wme(:)

 ABI_ALLOCATE(free,    (ntemper))
 ABI_ALLOCATE(energy,  (ntemper))
 ABI_ALLOCATE(entropy, (ntemper))
 ABI_ALLOCATE(spheat,  (ntemper))
 ABI_ALLOCATE(wme,     (ntemper))

!Put zeroes for F, E, S, Cv
 free(:)=zero
 energy(:)=zero
 entropy(:)=zero
 spheat(:)=zero
 wme(:)=zero

!open THERMO file
 if (open_file(fname,msg,newunit=thermal_unit,form="formatted",action="write") /= 0) then
   MSG_ERROR(msg)
 end if

! print header
 write(msg,'(a,a)') ch10,&
&  ' # At  T     F(J/mol-c)     E(J/mol-c)     S(J/(mol-c.K)) C(J/(mol-c.K)) Omega_mean(cm-1) from prtdos DOS'
 call wrtout(thermal_unit,msg,'COLL')
 msg = ' # (A mol-c is the abbreviation of a mole-cell, that is, the'
 call wrtout(thermal_unit,msg,'COLL')
 msg = ' #  number of Avogadro times the atoms in a unit cell)'
 call wrtout(thermal_unit,msg,'COLL')

 write(msg, '(a,a,a)' )&
&  ' phdos_print_thermo : thermodynamic functions calculated from prtdos DOS (not histogram)',ch10,&
&  '     see THERMO output file ...'
 call wrtout(std_out,msg,'COLL')

 domega = (PHdos%omega(2)-PHdos%omega(1))

 do itemper=1,ntemper

!  The temperature (tmp) is given in Ha
   tmp=(tempermin+temperinc*dble(itemper-1))*kb_HaK

   do iomega=1,PHdos%nomega
     if (abs(PHdos%phdos(iomega)) < 1.e-200_dp) cycle

!    wover2t= hbar*w / 2kT dimensionless
     wover2t = zero;     if(tmp > tol14) wover2t=PHdos%omega(iomega)*half/tmp
     ! should not be much of a problem for the log, but still put a check.
     ln2shx=zero;        if (wover2t > tol16 .and. wover2t < 100.0_dp) ln2shx=log(two * sinh(wover2t))
     cothx=zero;         if (wover2t > tol16) cothx=one/tanh(wover2t)
     invsinh2=zero;      if (wover2t > tol16 .and. wover2t < 100.0_dp) invsinh2=one/sinh(wover2t)**2

!    This matches the equations published in Lee & Gonze, PRB 51, 8610 (1995)
     free(itemper)   = free(itemper)    + PHdos%phdos(iomega)*tmp*ln2shx
     energy(itemper) = energy(itemper)  + PHdos%phdos(iomega)*half*PHdos%omega(iomega)*cothx
     spheat(itemper) = spheat(itemper)  + PHdos%phdos(iomega)*wover2t**2 * invsinh2
     entropy(itemper)= entropy(itemper) + PHdos%phdos(iomega)*(wover2t*cothx - ln2shx)
     wme(itemper)    = wme(itemper)     + PHdos%phdos(iomega)*PHdos%omega(iomega)*wover2t**2 * invsinh2

   end do ! iomega

   ! suppose homogeneous omega grid and multiply by domega
   free(itemper)   = free(itemper)    * domega
   energy(itemper) = energy(itemper)  * domega
   entropy(itemper)= entropy(itemper) * domega
   spheat(itemper) = spheat(itemper)  * domega
   wme(itemper)    = wme(itemper)     * domega

   if (abs(spheat(itemper))>tol8) wme(itemper)=wme(itemper)/spheat(itemper)

! do the printing to file
   write(msg,'(es11.3,5es15.7)') tmp/kb_HaK,&
&    Ha_J*Avogadro*free(itemper),&
&    Ha_J*Avogadro*energy(itemper),&
&    Ha_J*Avogadro*kb_HaK*entropy(itemper),&
&    Ha_J*Avogadro*kb_HaK*spheat(itemper),&
&    wme(itemper)*Ha_cmm1
   call wrtout(thermal_unit,msg,'COLL')
 end do ! itemper

 ABI_DEALLOCATE(free)
 ABI_DEALLOCATE(energy)
 ABI_DEALLOCATE(entropy)
 ABI_DEALLOCATE(spheat)
 ABI_DEALLOCATE(wme)

 close(thermal_unit)

end subroutine phdos_print_thermo
!!***
!----------------------------------------------------------------------

!!****f* m_phonons/phdos_free
!!
!! NAME
!! phdos_free
!!
!! FUNCTION
!! destructor function for phonon DOS object
!!
!! INPUTS
!! PHdos= container object for phonon DOS
!!
!! OUTPUT
!!
!! PARENTS
!!      anaddb,eph
!!
!! CHILDREN
!!      ifc_fourq,kpath_free,phonons_ncwrite,phonons_write_gnuplot
!!      phonons_write_phfrq,phonons_write_xmgrace,xmpi_sum_master
!!
!! SOURCE

subroutine phdos_free(PHdos)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phdos_free'
!End of the abilint section

 implicit none

!Arguments -------------------------------
 type(phonon_dos_type),intent(inout) ::PHdos

! *************************************************************************

 !@phonon_dos_type
 if (allocated(PHdos%atom_mass)) then
   ABI_FREE(PHdos%atom_mass)
 end if
 if (allocated(PHdos%omega)) then
   ABI_FREE(PHdos%omega)
 end if
 if (allocated(PHdos%phdos)) then
   ABI_FREE(PHdos%phdos)
 end if
 if (allocated(PHdos%phdos_int)) then
   ABI_FREE(PHdos%phdos_int)
 end if
 if (allocated(PHdos%pjdos)) then
   ABI_FREE(PHdos%pjdos)
 end if
 if (allocated(PHdos%pjdos_int)) then
   ABI_FREE(PHdos%pjdos_int)
 end if
 if (allocated(PHdos%pjdos_type)) then
   ABI_FREE(PHdos%pjdos_type)
 end if
 if (allocated(PHdos%pjdos_type_int)) then
   ABI_FREE(PHdos%pjdos_type_int)
 end if
 if (allocated(PHdos%pjdos_rc_type)) then
   ABI_FREE(PHdos%pjdos_rc_type)
 end if
 if (allocated(PHdos%msqd_dos_atom)) then
   ABI_FREE(PHdos%msqd_dos_atom)
 end if

end subroutine phdos_free
!!***

!----------------------------------------------------------------------

!!****f* m_phonons/mkphdos
!!
!! NAME
!! mkphdos
!!
!! FUNCTION
!! Function to calculate the phonon density of states as well as
!! the contributions associated to the different types of atoms in the unit cell.
!! Two methods are implemented: gaussian method and linear interpolation based on
!! tetrahedrons.
!!
!! INPUTS
!! Ifc<ifc_type>=Interatomic force constants
!! Crystal<crystal_t>=Info on the crystalline Structure.
!! prtdos=1 for Gaussian method, 2 for tetrahedra.
!! dosdeltae=Step for the frequency mesh.
!! dossmear=Gaussian broadening used if prtdos==1.
!! dos_ngqpt(3)=Divisions of the q-mesh used for computing the DOS
!! dos_qshift(3)=Shift of the q-mesh.
!! comm=MPI communicator
!!
!! OUTPUT
!! PHdos<phonon_dos_type>=Container with phonon DOS, IDOS and atom-projected DOS.
!!
!! NOTES
!! On the use of the q-grids:
!! Two different q-meshes are used in this subroutine. The first one is the coarse
!! mesh where the interatomic forces have been calculated during the DFPT run.
!! This q-grid is used to obtain an initial guess for the max and min frequency
!! value of the phonon spectrum. These values are, indeed, required to dimension
!! the array containing the PHDOS. The second (dense) grid is used to perform the
!! PHDOS calculation. If the Fourier interpolation on the second dense q-grid
!! generates a phonon frequency outside the initially calculated frequency mesh,
!! the mesh is enlarged and the calculation is restarted.
!!
!! PARENTS
!!      anaddb,eph,m_tdep_phdos
!!
!! CHILDREN
!!      ifc_fourq,kpath_free,phonons_ncwrite,phonons_write_gnuplot
!!      phonons_write_phfrq,phonons_write_xmgrace,xmpi_sum_master
!!
!! SOURCE

subroutine mkphdos(PHdos,Crystal,Ifc,prtdos,dosdeltae,dossmear,dos_ngqpt,&
&   dos_qshift, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkphdos'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: prtdos
 integer,intent(in) :: comm
 real(dp),intent(in) :: dosdeltae,dossmear
 type(crystal_t),intent(in) :: Crystal
 type(ifc_type),intent(in) :: Ifc
 type(phonon_dos_type),intent(inout) :: PHdos
!arrays
 integer,intent(in) :: dos_ngqpt(3)
 real(dp),intent(in) :: dos_qshift(3)

!Local variables -------------------------
!scalars
 integer,parameter :: brav1=1,chksymbreak0=0,bcorr0=0,qptopt1=1,master=0
 integer :: iat,jat,idir,imesh,imode,io,iq_ibz,itype,nkpt_fullbz
 integer :: nmesh,nqbz,nqshft,ierr,natom,nomega
 integer :: jdir, isym
 integer :: nprocs, my_rank
 integer :: ncid
 real(dp) :: nsmallq
 real(dp) :: dum,gaussfactor,gaussprefactor,gaussval,low_bound,max_occ !,pnorm
 real(dp) :: upr_bound,xx,gaussmaxarg
 real(dp) :: max_smallq = 0.0625_dp
 real(dp) :: normq
 real(dp) :: cpu, wall, gflops
 logical :: out_of_bounds
 character(len=500) :: msg
 character(len=80) :: errstr
 character(len=80) ::  prefix = "freq_displ"
 type(t_tetrahedron) :: tetraq
!arrays
 integer :: qptrlatt(3,3)
 integer,allocatable :: bz2ibz(:),ngqpt(:,:)
 real(dp) :: speedofsound(3)
 real(dp) :: speedofsound_(3)
 real(dp) :: debyefreq
 real(dp) :: displ(2*3*Crystal%natom*3*Crystal%natom)
 real(dp) :: eigvec(2,3,Crystal%natom,3*Crystal%natom),phfrq(3*Crystal%natom)
 real(dp) :: qlatt(3,3),rlatt(3,3)
 real(dp) :: msqd_atom_tmp(3,3),temp_33(3,3)
 real(dp) :: symcart(3,3,crystal%nsym)
 real(dp),allocatable :: dtweightde(:,:),full_eigvec(:,:,:,:,:),full_phfrq(:,:)
 real(dp),allocatable :: kpt_fullbz(:,:),qbz(:,:),qibz(:,:),qshft(:,:),tmp_phfrq(:),tweight(:,:)
 real(dp),allocatable :: wtqibz(:)
 real(dp), allocatable :: pjdos_tmp(:,:,:)
 real(dp),allocatable :: xvals(:), gvals_wtq(:), wdt(:,:)
 real(dp) :: syme2_xyza(3, crystal%natom)

! *********************************************************************

 DBG_ENTER("COLL")

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 ! Consistency check.
 if (ALL(prtdos /= [1,2])) then
   MSG_BUG(sjoin('prtdos should be 1 or 2, but received', itoa(prtdos)))
 end if
 if (dosdeltae<=zero) then
   MSG_BUG(sjoin('dosdeltae should be positive, but received', ftoa(dosdeltae)))
 end if
 if (prtdos==1.and.dossmear<=zero) then
   MSG_BUG(sjoin('dossmear should be positive but received', ftoa(dossmear)))
 end if

 call cwtime(cpu, wall, gflops, "start")

 natom = Crystal%natom
 gaussmaxarg = sqrt(-log(1.d-90))

 do isym = 1, Crystal%nsym
   call symredcart(Crystal%rprimd,Crystal%gprimd,symcart(:,:,isym),Crystal%symrel(:,:,isym))
 end do

 nomega = 1
 PHdos%ntypat     = crystal%ntypat
 PHdos%natom      = natom
 PHdos%prtdos     = prtdos
 PHdos%nomega     = 1
 PHdos%nqibz      = 1
 PHdos%omega_max  = smallest_real
 PHdos%omega_min  = greatest_real
 PHdos%omega_step = dosdeltae
 PHdos%dossmear   = dossmear

 ABI_MALLOC(PHdos%atom_mass, (natom))
 PHdos%atom_mass(:) = Crystal%amu(Crystal%typat(:))*amu_emass

 ! These arrays are independent of nqibz
 ABI_MALLOC(PHdos%omega, (nomega))
 ABI_MALLOC(PHdos%phdos, (nomega))
 ABI_MALLOC(PHdos%phdos_int, (nomega))
 ABI_MALLOC(PHdos%pjdos, (nomega,3,natom))
 ABI_MALLOC(PHdos%pjdos_int, (nomega,3,natom))
 ABI_MALLOC(PHdos%pjdos_type, (nomega,crystal%ntypat))
 ABI_MALLOC(PHdos%pjdos_type_int, (nomega,crystal%ntypat))
 ABI_MALLOC(PHdos%pjdos_rc_type, (nomega,3,crystal%ntypat))
 ABI_MALLOC(PHdos%msqd_dos_atom, (nomega,3,3,natom))
 ABI_CALLOC(pjdos_tmp,(2,3,natom))

 ! Parameters defining the gaussian approximant.
 if (prtdos==1) then
   ! TODO: use dirac_delta and update reference files.
   gaussprefactor=one/(dossmear*sqrt(two_pi))
   gaussfactor=one/(sqrt2*dossmear)
   write(msg,'(4a,f8.5,2a,f8.5)')ch10,&
&   ' mkphdos: calculating phonon DOS using gaussian method :',ch10,&
&   '    gaussian smearing [meV] = ',dossmear*Ha_meV,ch10,&
&   '    frequency step    [meV] = ',PHdos%omega_step*Ha_meV
 else if (prtdos==2) then
   write(msg,'(2a)')ch10,' mkphdos: calculating phonon DOS using tetrahedron method'
 end if
 call wrtout(std_out,msg,'COLL')

 ! Initial lower and upper bound of the phonon spectrum.
 low_bound=dble(huge(0))*half*PHdos%omega_step
 upr_bound=-low_bound

 ! TODO
 ! 1) fix bug in tetra if degenerate and update ref files
 ! 2) Get rid of double mesh algo and use min/max freq already computed in ifc
 ! 3) Add MPI Parallel version (easy after change in 2)

 nmesh=2
 ABI_MALLOC(ngqpt,(3,nmesh))
 do imesh=1,nmesh

   write(msg,'(a,I6)') ' Mesh number ', imesh
   call wrtout(std_out,msg,'COLL')

   if (imesh==1) then
     ! Coarse q-mesh used in RF calculation.
     ngqpt(:,imesh)=ifc%ngqpt(1:3)
     nqshft=ifc%nqshft
     ABI_MALLOC(qshft,(3,nqshft))
     ! TODO this has to be fixed  there is a small inconsistency in the dimension of q1shft
     qshft(:,1:nqshft)=ifc%qshft(:,1:nqshft)
   else
     ! Dense q-mesh used for the Fourier interpolation.
     !ngqpt(1:3,imesh)=inp%ng2qpt(1:3)
     ngqpt(1:3,imesh) = dos_ngqpt
     nqshft=1 !always 1
     ABI_MALLOC(qshft,(3,nqshft))
     !qshft(:,1)=inp%q2shft(:)  ! FIXME small inconsistency in the dimension of q1shft
     qshft(:,1)=dos_qshift(:)
   end if

   qptrlatt=0
   qptrlatt(1,1)=ngqpt(1,imesh); qptrlatt(2,2)=ngqpt(2,imesh); qptrlatt(3,3)=ngqpt(3,imesh)

   if(allocated(wtqibz)) then
     ABI_FREE(wtqibz)
   end if
   if(allocated(qibz)) then
     ABI_FREE(qibz)
   end if

   ! This call will set %nqibz, IBZ and BZ arrays
   call kpts_ibz_from_kptrlatt(crystal, qptrlatt, qptopt1, nqshft, qshft, &
     phdos%nqibz, qibz, wtqibz, nqbz, qbz) ! new_kptrlatt, new_shiftk)  ! Optional
   ABI_FREE(qshft)

   if (prtdos==2.and.imesh==2) then
     ! Second mesh with tetrahedron method
     ! convert kptrlatt to double and invert, qlatt here refer to the shortest qpt vectors
     rlatt = qptrlatt; call matr3inv(rlatt,qlatt)

     ABI_MALLOC(qshft,(3,nqshft))
     !qshft(:,1)=inp%q2shft(:)  ! FIXME small inconsistency in the dimension of q1shft
     qshft(:,1)= dos_qshift(:)
     nkpt_fullbz=nqbz
     ABI_MALLOC(bz2ibz,(nkpt_fullbz))
     ABI_MALLOC(kpt_fullbz,(3,nkpt_fullbz))

     ! Make full kpoint grid and get equivalence to irred kpoints.
     ! This routines scales **very badly** wrt nkpt_fullbz, should introduce check on the norm.
     call get_full_kgrid(bz2ibz,qibz,kpt_fullbz,qptrlatt,PHdos%nqibz,&
&      nkpt_fullbz,nqshft,Crystal%nsym,qshft,Crystal%symrel)

     ! Get tetrahedra, ie indexes of the full q-points at their summits
     call init_tetra(bz2ibz, crystal%gprimd, qlatt, kpt_fullbz, nqbz, tetraq, ierr, errstr)
     ABI_CHECK(ierr==0,errstr)

     ABI_FREE(qshft)
     ABI_FREE(bz2ibz)
     ABI_FREE(kpt_fullbz)

     ! Allocate arrays used to store the entire spectrum, Required to calculate tetra weights.
     ! this may change in the future if Matteo refactorizes the tetra weights as sums over k instead of sums over bands
     ABI_MALLOC(full_phfrq,(3*natom,PHdos%nqibz))
     ABI_STAT_MALLOC(full_eigvec,(2,3,natom,3*natom,PHdos%nqibz), ierr)
     ABI_CHECK(ierr==0, 'out-of-memory in full_eigvec')
   end if  ! prtdos==2.and.imesh==2

   ! This infinite loop is used to be sure that the frequency mesh is large enough to contain
   ! the entire phonon spectrum. The mesh is enlarged if, during the Fourier interpolation,
   ! a phonon frequency turns out to be outside the interval [omega_min:omega_max]
   do
     out_of_bounds=.FALSE.
     if (allocated(PHdos%omega)) then
       ABI_FREE(PHdos%omega)
     end if
     if (allocated(PHdos%phdos)) then
       ABI_FREE(PHdos%phdos)
     end if
     if (allocated(PHdos%pjdos)) then
       ABI_FREE(PHdos%pjdos)
     end if
     if (allocated(PHdos%msqd_dos_atom)) then
       ABI_FREE(PHdos%msqd_dos_atom)
     end if

     ! Frequency mesh.
     PHdos%omega_min=low_bound; if (ABS(PHdos%omega_min)<tol5) PHdos%omega_min=-tol5
     PHdos%omega_max=upr_bound
     PHdos%nomega=NINT((PHdos%omega_max-PHdos%omega_min)/PHdos%omega_step)+1
     PHdos%nomega=MAX(6,PHdos%nomega) ! Ensure Simpson integration will be ok

     ABI_MALLOC(PHdos%omega,(PHdos%nomega))
     do io=1,PHdos%nomega
       PHdos%omega(io)=PHdos%omega_min+PHdos%omega_step*(io-1)
     end do

     ABI_CALLOC(PHdos%phdos,(PHdos%nomega))
     ABI_CALLOC(PHdos%pjdos,(PHdos%nomega,3,natom))
     ABI_CALLOC(PHdos%msqd_dos_atom,(PHdos%nomega,3,3,natom))

     if (allocated(gvals_wtq)) then
       ABI_FREE(gvals_wtq)
     end if
     if (allocated(xvals)) then
       ABI_FREE(xvals)
     end if
     ABI_MALLOC(gvals_wtq, (phdos%nomega))
     ABI_MALLOC(xvals, (phdos%nomega))

     ! Sum over irreducible q-points
     nsmallq = zero
     speedofsound = zero
     do iq_ibz=1,PHdos%nqibz
       !if (mod(iq_ibz, nprocs) /= my_rank) cycle ! mpi-parallelism

       !TODO: this is where the mpi distribution should happen.
       ! then sync the following integrals:
       !   speedofsound
       !   nsmallq
       !   PHdos%phdos
       !   PHdos%msqd_dos_atom
       !   full_phfrq
       !   full_eigvec

       ! Fourier interpolation.
       call ifc_fourq(Ifc,Crystal,qibz(:,iq_ibz),phfrq,displ,out_eigvec=eigvec)

       dum=MINVAL(phfrq); PHdos%omega_min=MIN(PHdos%omega_min,dum)
       dum=MAXVAL(phfrq); PHdos%omega_max=MAX(PHdos%omega_max,dum)
       out_of_bounds = (PHdos%omega_min<low_bound .or. PHdos%omega_max>upr_bound)

       normq = sum(qibz(:,iq_ibz)**2)
       if (normq < max_smallq .and. normq > tol6) then
          call phdos_calc_vsound(eigvec, Crystal%gmet, natom, phfrq, qibz(:,iq_ibz), speedofsound_)
         speedofsound = speedofsound + speedofsound_*wtqibz(iq_ibz)
         nsmallq = nsmallq + wtqibz(iq_ibz)
       end if

       if (imesh>1.and..not.out_of_bounds) then
         select case (prtdos)
         case (1)
           do imode=1,3*natom
             ! Precompute \delta(w - w_{qnu}) * weight(q)
             xvals = (PHdos%omega(:) - phfrq(imode)) * gaussfactor
             where (abs(xvals) < gaussmaxarg)
               gvals_wtq = gaussprefactor * exp(-xvals*xvals) * wtqibz(iq_ibz)
             elsewhere
               gvals_wtq = zero
             end where

             ! Accumulate PHDOS
             PHdos%phdos(:) = PHdos%phdos(:) + gvals_wtq

             ! Rotate e(q) to get e(Sq) to account for symmetrical q-points in BZ.
             ! eigenvectors indeed are not invariant under rotation. See e.g. Eq 39-40 of PhysRevB.76.165108
             ! In principle there's a phase due to nonsymmorphic translations
             ! but we here need |e(Sq)_iatom|**2
             syme2_xyza = zero
             do iat=1,natom
               do isym=1, Crystal%nsym
                 jat = Crystal%indsym(4,isym,iat)
                 syme2_xyza(:,jat) = syme2_xyza(:,jat) + &
                   matmul(symcart(:,:,isym), eigvec(1,:,iat,imode)) ** 2 + &
                   matmul(symcart(:,:,isym), eigvec(2,:,iat,imode)) ** 2
               end do
             end do
             !syme2_xyza = syme2_xyza / crystal%nsym

             ! Accumulate PJDOS
             do iat=1,natom
               do idir=1,3
                 PHdos%pjdos(:,idir,iat) = PHdos%pjdos(:,idir,iat) + syme2_xyza(idir,iat) * gvals_wtq
               end do ! idir
             end do

             ! Accumulate outer product of displacement vectors
             ! NB: only accumulate real part. e(-q) = e(q)* thue full sum over the BZ guarantees Im=0
             ! this sum only does irreducible points: the matrix is symmetrized below
             ! msqd_atom_tmp has units of bohr^2 / Ha as gaussval ~ 1/smear ~ 1/Ha
             do iat=1,natom
               msqd_atom_tmp = zero
               do idir=1,3
                 do jdir=1,3
                   msqd_atom_tmp(jdir,idir) = msqd_atom_tmp(jdir,idir) + ( &
                         eigvec(1,idir,iat,imode)* eigvec(1,jdir,iat,imode) &
                      +  eigvec(2,idir,iat,imode)* eigvec(2,jdir,iat,imode) )
                 end do
               end do

               ! Symmetrize matrices to get full sum of tensor over all BZ, not just IBZ.
               ! the atom is not necessarily invariant under symops, so these contributions should be added to each iat separately
               ! normalization by nsym is done at the end outside the iqpt loop and after the tetrahedron clause
               ! NB: looks consistent with the sym in harmonic thermo, just used in opposite
               ! direction for symops: symrel here instead of symrec and the inverse of
               ! indsym in harmonic_thermo
               do isym=1, Crystal%nsym
                 temp_33 = matmul( (symcart(:,:,isym)), matmul(msqd_atom_tmp, transpose(symcart(:,:,isym))) )
                 jat = Crystal%indsym(4,isym,iat)
                 do idir=1,3
                   do jdir=1,3
                     PHdos%msqd_dos_atom(:,idir,jdir,jat) = PHdos%msqd_dos_atom(:,idir,jdir,jat) + &
                       temp_33(idir, jdir) * gvals_wtq
                   end do
                 end do
               end do

             end do ! iat
           end do ! imode

         case (2)
           ! Tetrahedrons
           ! * Save phonon frequencies and eigenvectors.
           ! * Sum is done after the loops over the two meshes.
           full_phfrq(:,iq_ibz)=phfrq(:)
           full_eigvec(:,:,:,:,iq_ibz)=eigvec

         case default
           MSG_ERROR(sjoin("Wrong value for prtdos:", itoa(prtdos)))
         end select

       end if !Second mesh and not out of boundaries
     end do !irred q-points

     if (out_of_bounds) then
       upr_bound=PHdos%omega_max+ABS(PHdos%omega_max/ten)
       low_bound=PHdos%omega_min-ABS(PHdos%omega_min/ten)
       write(msg,'(3a)')&
&       ' At least one phonon frequency falls outside the frequency mesh chosen',ch10,&
&       ' restarting the calculation with a larger frequency mesh '
       if (imesh>1) then
         MSG_COMMENT(msg)
       end if
     else
       EXIT !infinite loop
     end if
   end do !infinite loop

   if (nsmallq > tol10) then
     speedofsound = speedofsound/nsmallq

     write (msg,'(a,E20.10,3a,F16.4,2a)') ' Average speed of sound partial sums: ', third*sum(speedofsound), ' (at units)',ch10,&
&               '-                                   = ', third*sum(speedofsound) * Bohr_Ang * 1.d-13 / Time_Sec, ' [km/s]',ch10
     call wrtout (ab_out,msg,"COLL")
     call wrtout (std_out,msg,"COLL")

     ! Debye frequency = vs * (6 pi^2 natom / ucvol)**1/3
     debyefreq = third*sum(speedofsound) * (six*pi**2/Crystal%ucvol)**(1./3.)
     write (msg,'(a,E20.10,3a,E20.10,a)') ' Debye frequency from partial sums: ', debyefreq, ' (Ha)',ch10,&
&                                         '-                                 = ', debyefreq*Ha_THz, ' (THz)'
     call wrtout (ab_out,msg,"COLL")
     call wrtout (std_out,msg,"COLL")

     ! Debye temperature = hbar * Debye frequency / kb
     write (msg,'(a,E20.10,2a)') '-Debye temperature from partial sums: ', debyefreq*Ha_K, ' (K)', ch10
     call wrtout (ab_out,msg,"COLL")
     call wrtout (std_out,msg,"COLL")
   end if

   ABI_FREE(qbz)
 end do !imesh

 ABI_FREE(ngqpt)
 if (allocated(gvals_wtq)) then
   ABI_FREE(gvals_wtq)
 end if
 if (allocated(xvals)) then
   ABI_FREE(xvals)
 end if
 if (allocated(PHdos%phdos_int)) then
   ABI_FREE(PHdos%phdos_int)
 end if
 if (allocated(PHdos%pjdos_int)) then
   ABI_FREE(PHdos%pjdos_int)
 end if

 ABI_CALLOC(PHdos%phdos_int,(PHdos%nomega))

 if (prtdos==2) then
   ! Finalize integration with tetrahedra
   ! * All the data are contained in full_phfrq and full_eigvec.
   ! * low_bound and upr_bound contain the entire spectrum calculated on the dense mesh.
   ABI_MALLOC(tmp_phfrq,(PHdos%nqibz))
   ABI_MALLOC(PHdos%pjdos_int,(PHdos%nomega,3,natom))
   PHdos%phdos=zero; PHdos%pjdos=zero; PHdos%pjdos_int=zero
   max_occ=one

   ABI_MALLOC(wdt, (phdos%nomega, 2))
   ABI_MALLOC(tweight, (PHdos%nomega, PHdos%nqibz))
   ABI_MALLOC(dtweightde, (PHdos%nomega, PHdos%nqibz))

   !cnt = 0
   do imode=1,3*natom
     tmp_phfrq(:)= full_phfrq(imode,:)
     call tetra_blochl_weights(tetraq,tmp_phfrq,phdos%omega_min,phdos%omega_max,max_occ,phdos%nomega,&
        phdos%nqibz,bcorr0,tweight,dtweightde,comm)

     do iq_ibz=1,phdos%nqibz
       !cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle ! mpi-parallelism
       wdt(:, 2) = tweight(:, iq_ibz)
       wdt(:, 1) = dtweightde(:, iq_ibz)

       !call tetra_get_onewk(tetraq, iq_ibz, bcorr0, phdos%nomega, phdos%nqibz, &
       !  tmp_phfrq, phdos%omega_min, phdos%omega_max, max_occ, wdt)

       ! Accumulate DOS/IDOS
       PHdos%phdos(:) = PHdos%phdos(:) + wdt(:,1)
       PHdos%phdos_int(:) = PHdos%phdos_int(:) + wdt(:,2)

       ! Rotate e(q) to get e(Sq) to account for other q-points in BZ. See notes in gaussian branch
       syme2_xyza = zero
       do iat=1,natom
         do isym=1, Crystal%nsym
           jat = Crystal%indsym(4,isym,iat)
           syme2_xyza(:,jat) = syme2_xyza(:,jat) + &
             matmul(symcart(:,:,isym), full_eigvec(1,:,iat,imode,iq_ibz)) ** 2 + &
             matmul(symcart(:,:,isym), full_eigvec(2,:,iat,imode,iq_ibz)) ** 2
         end do
       end do
       !syme2_xyza = syme2_xyza / crystal%nsym

       do iat=1,natom
         do idir=1,3
           PHdos%pjdos(:,idir,iat) = PHdos%pjdos(:,idir,iat) + syme2_xyza(idir,iat) * wdt(:,1)
           PHdos%pjdos_int(:,idir,iat) = PHdos%pjdos_int(:,idir,iat) + syme2_xyza(idir,iat) * wdt(:,2)
         end do ! idir
       end do

       do iat=1,natom
         ! Accumulate outer product of displacement vectors
         msqd_atom_tmp = zero
         do idir=1,3
           do jdir=1,3
             msqd_atom_tmp(jdir,idir) = msqd_atom_tmp(jdir,idir) + ( &
                   full_eigvec(1,idir,iat,imode,iq_ibz)* full_eigvec(1,jdir,iat,imode,iq_ibz) &
                +  full_eigvec(2,idir,iat,imode,iq_ibz)* full_eigvec(2,jdir,iat,imode,iq_ibz) ) !* gvals_wtq
           end do ! jdie
         end do

         ! Symmetrize matrices to get full sum of tensor over all BZ, not just IBZ.
         ! the atom is not necessarily invariant under symops, so these contributions should be added to each iat separately
         ! normalization by nsym is done at the end outside the iqpt loop and after the tetrahedron clause
         ! from loops above only the eigvec are kept and not the displ, so we still have to divide by the masses
         ! TODO:  need to check the direction of the symcart vs transpose or inverse, given that jat is the pre-image of iat...
         do isym=1, Crystal%nsym
           temp_33 = matmul( (symcart(:,:,isym)), matmul(msqd_atom_tmp, transpose(symcart(:,:,isym))) )
           jat = Crystal%indsym(4,isym,iat)
           do idir=1,3
             do jdir=1,3
               PHdos%msqd_dos_atom(:,idir,jdir,jat) = PHdos%msqd_dos_atom(:,idir,jdir,jat) + &
                 temp_33(idir, jdir) * wdt(:,1)
             end do
           end do
         end do
       end do ! iat

     end do ! imode
   end do ! iq_ibz

   ABI_FREE(wdt)

   ! Make eigvec into displacements
   do iat = 1, natom
     full_eigvec(:,:,iat,:,:) = full_eigvec(:,:,iat,:,:)/sqrt(PHdos%atom_mass(iat))
   end do

   if (my_rank == master) then
#ifdef HAVE_NETCDF
     ! TODO: should pass prefix as arg.
     NCF_CHECK_MSG(nctk_open_create(ncid, strcat(prefix, "_PHIBZ.nc"), xmpi_comm_self), "Creating PHIBZ")
     NCF_CHECK(crystal_ncwrite(Crystal, ncid))
     call phonons_ncwrite(ncid,natom,PHdos%nqibz, qibz, wtqibz, full_phfrq, full_eigvec)
     NCF_CHECK(nf90_close(ncid))
#endif
   end if

! immediately free this - it contains displ and not eigvec at this stage
   ABI_FREE(full_eigvec)
   ABI_FREE(full_phfrq)
   ABI_FREE(tmp_phfrq)
   ABI_FREE(tweight)
   ABI_FREE(dtweightde)
   call destroy_tetra(tetraq)
 else
#ifdef HAVE_NETCDF
   MSG_WARNING('The netcdf PHIBZ file is only output for tetrahedron integration and DOS calculations')
#endif
 end if ! prtdos 2 = tetrahedra

! normalize by nsym : symmetrization is used in all prtdos cases
 PHdos%msqd_dos_atom = PHdos%msqd_dos_atom / Crystal%nsym
 PHdos%pjdos = PHdos%pjdos / Crystal%nsym
 if (prtdos == 2) PHdos%pjdos_int = PHdos%pjdos_int / Crystal%nsym

 ABI_FREE(pjdos_tmp)
 ABI_FREE(qibz)
 ABI_FREE(wtqibz)

! normalize by mass and factor of 2 ! now added in the printout to agree with harmonic_thermo
! do iat=1, natom
!   PHdos%msqd_dos_atom(:,:,:,iat) = PHdos%msqd_dos_atom(:,:,:,iat) * invmass(iat) * half
! end do ! iat

 ! =======================
 ! === calculate IPDOS ===
 ! =======================
 if (allocated(PHdos%pjdos_rc_type)) then
   ABI_FREE(PHdos%pjdos_rc_type)
 end if
 if (allocated(PHdos%pjdos_type)) then
   ABI_FREE(PHdos%pjdos_type)
 end if
 if (allocated(PHdos%pjdos_type_int)) then
   ABI_FREE(PHdos%pjdos_type_int)
 end if

 ABI_CALLOC(PHdos%pjdos_rc_type, (PHdos%nomega,3,Crystal%ntypat))
 ABI_CALLOC(PHdos%pjdos_type, (PHdos%nomega,Crystal%ntypat))
 ABI_CALLOC(PHdos%pjdos_type_int, (PHdos%nomega,Crystal%ntypat))

 do iat=1,natom
   itype=Crystal%typat(iat)
   do io=1,PHdos%nomega
     PHdos%pjdos_rc_type(io,:,itype)=PHdos%pjdos_rc_type(io,:,itype)+PHdos%pjdos(io,:,iat)
     PHdos%pjdos_type(io,itype)=PHdos%pjdos_type(io,itype)+sum(PHdos%pjdos(io,:,iat))
   end do
   if (prtdos==2) then
     do io=1,PHdos%nomega
       PHdos%pjdos_type_int(io,itype)=PHdos%pjdos_type_int(io,itype)+SUM(PHdos%pjdos_int(io,:,iat))
     end do
   end if
 end do

 ! Evaluate IDOS using simple simpson integration
 ! TODO should avoid the simpson rule using derf.F90, just to be consistent
 if (prtdos==1) then
   call simpson_int(PHdos%nomega,PHdos%omega_step,PHdos%phdos,PHdos%phdos_int)
   !do iat=1,natom
   !  do idir=1,3
   !    call simpson_int(PHdos%nomega,PHdos%omega_step,PHdos%pjdos(:,idir,iat),PHdos%pjdos_int(:,idir,iat))
   !  end do
   !end do
   do itype=1,Crystal%ntypat
     call simpson_int(PHdos%nomega,PHdos%omega_step,PHdos%pjdos_type(:,itype),PHdos%pjdos_type_int(:,itype))
   end do
 end if

 call cwtime(cpu,wall,gflops,"stop")
 write(msg,'(2(a,f8.2))')" mkphdos completed. cpu:",cpu,", wall:",wall
 call wrtout(std_out, msg, do_flush=.True.)

 DBG_EXIT("COLL")

end subroutine mkphdos
!!***

!----------------------------------------------------------------------

!!****f* m_phonons/thermal_supercell_make
!! NAME
!! thermal_supercell_make
!!
!! FUNCTION
!!  Construct an optimally thermalized supercell following Zacharias and Giustino
!! PRB 94 075125 (2016)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      ifc_fourq,kpath_free,phonons_ncwrite,phonons_write_gnuplot
!!      phonons_write_phfrq,phonons_write_xmgrace,xmpi_sum_master
!!
!! SOURCE

subroutine thermal_supercell_make(Crystal, Ifc, ntemper, &
&    rlatt, tempermin, temperinc, thm_scells)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'thermal_supercell_make'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ntemper
 integer, intent(in) :: rlatt(3,3)
 real(dp), intent(in) :: tempermin, temperinc
 type(crystal_t),intent(in) :: Crystal
 type(ifc_type),intent(in) :: Ifc
 type(supercell_type), intent(out) :: thm_scells(ntemper)

!Local variables-------------------------------
!scalars
 integer :: iq, nqibz, nqbz, qptopt1, imode, itemper, ierr, jmode
 real(dp) :: temperature_K, temperature, modesign, sigma, freeze_displ

!arrays
 integer, allocatable :: modeindex(:)
 real(dp), allocatable :: qshft(:,:) ! dummy with 2 dimensions for call to kpts_ibz_from_kptrlatt
 real(dp), allocatable :: qbz(:,:), qibz(:,:), wtqibz(:)
 real(dp), allocatable :: phfrq_allq(:), phdispl_allq(:,:,:,:,:)
 real(dp), allocatable :: phfrq(:), phdispl(:,:,:,:),pheigvec(:,:,:,:)
 real(dp), allocatable :: phdispl1(:,:,:)
 character (len=500) :: msg

! *************************************************************************

 ! check inputs
! TODO: add check that all rlatt are the same on input

 if (rlatt(1,2)/=0 .or.  rlatt(1,3)/=0 .or.  rlatt(2,3)/=0 .or. &
&    rlatt(2,1)/=0 .or.  rlatt(3,1)/=0 .or.  rlatt(3,2)/=0) then
   write (msg, '(4a, 9I6, a)') ' for the moment I have not implemented ', &
&    ' non diagonal supercells.',ch10,' rlatt for temp 1 = ', rlatt, ' Returning '
   MSG_WARNING(msg)
   return
 end if

 ! build qpoint grid used for the Fourier interpolation.
 !(use no syms for the moment!)
 qptopt1 = 3

 ! for the moment do not allow shifted q grids.
 ! We are interpolating anyway, so it will always work
 ABI_MALLOC(qshft,(3,1))
 qshft(:,1)=zero

 ! This call will set nqibz, IBZ and BZ arrays
 call kpts_ibz_from_kptrlatt(crystal, rlatt, qptopt1, 1, qshft, &
&   nqibz, qibz, wtqibz, nqbz, qbz) ! new_kptrlatt, new_shiftk)  ! Optional
 ABI_FREE(qshft)

 ! allocate arrays with all of the q, omega, and displacement vectors
 ABI_STAT_MALLOC(phfrq_allq, (3*Crystal%natom*nqibz), ierr)
 ABI_CHECK(ierr==0, 'out-of-memory in phfrq_allq')
 ABI_STAT_MALLOC(phdispl_allq, (2, 3, Crystal%natom, 3*Crystal%natom, nqibz), ierr)
 ABI_CHECK(ierr==0, 'out-of-memory in phdispl_allq')

 ABI_STAT_MALLOC(phfrq, (3*Crystal%natom), ierr)
 ABI_CHECK(ierr==0, 'out-of-memory in phfrq_allq')
 ABI_STAT_MALLOC(phdispl, (2, 3, Crystal%natom, 3*Crystal%natom), ierr)
 ABI_CHECK(ierr==0, 'out-of-memory in phdispl_allq')
 ABI_STAT_MALLOC(pheigvec, (2, 3, Crystal%natom, 3*Crystal%natom), ierr)
 ABI_CHECK(ierr==0, 'out-of-memory in phdispl_allq')

 ! loop over q to get all frequencies and displacement vectors
 ABI_ALLOCATE(modeindex, (nqibz*3*Crystal%natom))
 imode = 0
 do iq = 1, nqibz
   ! Fourier interpolation.
   call ifc_fourq(Ifc, Crystal, qibz(:,iq), phfrq, phdispl, out_eigvec=pheigvec)
   phfrq_allq((iq-1)*3*Crystal%natom+1 : iq*3*Crystal%natom) = phfrq
   phdispl_allq(1:2, 1:3, 1:Crystal%natom, 1:3*Crystal%natom, iq) = phdispl
   do jmode = 1, 3*Crystal%natom
     imode = imode + 1
     modeindex(imode) = imode
   end do
 end do
 ABI_FREE(phfrq)
 ABI_FREE(pheigvec)
 ABI_FREE(phdispl)

 ! sort modes in whole list: get indirect indexing for qbz and displ
 call sort_dp(nqibz*3*Crystal%natom, phfrq_allq, modeindex, tol10)
 ! NB: phfrq is sorted now, but displ and qibz will have to be indexed indirectly with modeindex

 ! only diagonal supercell case for the moment
 do itemper = 1, ntemper
   call init_supercell(Crystal%natom, rlatt, Crystal%rprimd, Crystal%typat, Crystal%xcart, Crystal%znucl, thm_scells(itemper))
 end do

 ! precalculate phase factors???

 ABI_STAT_MALLOC(phdispl1, (2, 3, Crystal%natom), ierr)
 ! for all modes at all q in whole list, sorted
 modesign=one
 do imode = 1, 3*Crystal%natom*nqibz
   ! skip modes with too low or negative frequency -> Bose factor explodes (eg acoustic at Gamma)
   if (phfrq_allq(imode) < tol10) cycle

   iq = ceiling(dble(modeindex(imode))/dble(3*Crystal%natom))
   jmode = modeindex(imode) - (iq-1)*3*Crystal%natom
   phdispl1 = phdispl_allq(:,:,:,jmode,iq)

   ! loop over temperatures
   do itemper = 1, ntemper
     temperature_K = tempermin + dble(itemper-1)*temperinc  ! this is in Kelvin
     temperature = temperature_K / Ha_K !=315774.65_dp

     ! trick supercell object into using present q point
     thm_scells(itemper)%qphon(:) = qibz(:,iq)

     ! find thermal displacement amplitude eq 4 of Zacharias
     !   combined with l_nu,q expression in paragraph before
     sigma = sqrt( (bose_einstein(phfrq_allq(imode), temperature) + half)/phfrq_allq(imode) )

     ! add displacement for this mode to supercell positions eq 5 of Zacharias
     freeze_displ = modesign * sigma
     call freeze_displ_supercell (phdispl1(:,:,:), freeze_displ, thm_scells(itemper))
   end do !itemper

   ! this is the prescription: flip sign for each successive mode in full
   ! spectrum, to cancel electron phonon coupling to 1st order
   ! (hopeflly 3rd order as well)
   modesign=-modesign

 end do !imode

 ABI_FREE(modeindex)
 ABI_FREE(phfrq_allq)
 ABI_FREE(phdispl_allq)
 ABI_FREE(phdispl1)
 ABI_FREE(qibz)
 ABI_FREE(qbz)
 ABI_FREE(wtqibz)

end subroutine thermal_supercell_make
!!***

!----------------------------------------------------------------------

!!****f* m_phonons/thermal_supercell_free
!! NAME
!! thermal_supercell_free
!!
!! FUNCTION
!!  deallocate thermal array of supercells
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      ifc_fourq,kpath_free,phonons_ncwrite,phonons_write_gnuplot
!!      phonons_write_phfrq,phonons_write_xmgrace,xmpi_sum_master
!!
!! SOURCE

subroutine thermal_supercell_free(ntemper, thm_scells)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'thermal_supercell_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ntemper
 type(supercell_type), allocatable, intent(inout) :: thm_scells(:)

! local
 integer :: itemp

 if(allocated(thm_scells)) then
   do itemp = 1, ntemper
     call destroy_supercell(thm_scells(itemp))
   end do
   ABI_FREE(thm_scells)
 end if
end subroutine thermal_supercell_free
!!***

!----------------------------------------------------------------------

!!****f* m_phonons/thermal_supercell_print
!! NAME
!! thermal_supercell_print
!!
!! FUNCTION
!!  print files with thermal array of supercells
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      ifc_fourq,kpath_free,phonons_ncwrite,phonons_write_gnuplot
!!      phonons_write_phfrq,phonons_write_xmgrace,xmpi_sum_master
!!
!! SOURCE

subroutine thermal_supercell_print(fname, ntemper, tempermin, temperinc, thm_scells)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'thermal_supercell_print'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ntemper
 real(dp), intent(in) :: tempermin
 real(dp), intent(in) :: temperinc
 type(supercell_type), intent(in) :: thm_scells(ntemper)
 character(len=fnlen), intent(in) :: fname

! local
 integer :: itemp
 character(len=80) :: title1, title2
 character(len=fnlen) :: filename
 real(dp) :: temper
 character(len=10) :: temper_str

 do itemp = 1, ntemper
   temper = dble(itemp-1)*temperinc+tempermin
   write (temper_str,'(I8)') int(temper)
   write (filename, '(3a)') trim(fname), "_T_", trim(adjustl(temper_str))
   write (title1, '(3a)') "#  thermalized supercell at temperature T= ", trim(temper_str), " Kelvin"
   title2 = "#  generated with alternating thermal displacements of all phonons"
   call prt_supercell (filename, thm_scells(itemp), title1, title2)
 end do

end subroutine thermal_supercell_print
!!***

!----------------------------------------------------------------------

!!****f* m_phonons/phdos_ncwrite
!! NAME
!! phdos_ncwrite
!!
!! FUNCTION
!!  Save the content of the object in a netcdf file.
!!
!! INPUTS
!!  ncid=NC file handle (open in the caller)
!!  phdos<phonon_dos_type>=Container object
!!
!! OUTPUT
!!  Only writing
!!
!! NOTES
!!  Frequencies are in eV, DOS are in states/eV.
!!
!! PARENTS
!!      anaddb,eph
!!
!! CHILDREN
!!      ifc_fourq,kpath_free,phonons_ncwrite,phonons_write_gnuplot
!!      phonons_write_phfrq,phonons_write_xmgrace,xmpi_sum_master
!!
!! SOURCE

subroutine phdos_ncwrite(phdos,ncid)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phdos_ncwrite'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(phonon_dos_type),intent(in) :: phdos
 integer,intent(in) :: ncid

!Local variables-------------------------------
!scalars
#ifdef HAVE_NETCDF
 integer :: ncerr

! *************************************************************************

! Define dimensions
 NCF_CHECK(nctk_def_basedims(ncid, defmode=.True.))

 ncerr = nctk_def_dims(ncid, [nctkdim_t("three", 3), nctkdim_t("number_of_atoms", phdos%natom),&
   nctkdim_t("number_of_atom_species", phdos%ntypat), nctkdim_t("number_of_frequencies", phdos%nomega)])
 NCF_CHECK(ncerr)

!scalars
 NCF_CHECK(nctk_def_iscalars(ncid, ["prtdos"]))
 NCF_CHECK(nctk_def_dpscalars(ncid, ["dossmear"]))

!arrays
 ncerr = nctk_def_arrays(ncid, [&
   nctkarr_t('wmesh', "dp", 'number_of_frequencies'),&
   nctkarr_t('phdos', "dp", 'number_of_frequencies'),&
   nctkarr_t('pjdos', "dp", 'number_of_frequencies, three, number_of_atoms'),&
   nctkarr_t('pjdos_type', "dp", 'number_of_frequencies, number_of_atom_species'),&
   nctkarr_t('pjdos_rc_type', "dp", 'number_of_frequencies, three, number_of_atom_species'), &
   nctkarr_t('msqd_dos_atom', "dp", 'number_of_frequencies, three, three, number_of_atoms') &
 ])
 NCF_CHECK(ncerr)

 ! Write variables. Note unit conversion.
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, vid("prtdos"), phdos%prtdos))
 NCF_CHECK(nf90_put_var(ncid, vid('dossmear'), phdos%dossmear*Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('wmesh'), phdos%omega*Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('phdos'), phdos%phdos/Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('pjdos'), phdos%pjdos/Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('pjdos_type'), phdos%pjdos_type/Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('pjdos_rc_type'), phdos%pjdos_rc_type/Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('msqd_dos_atom'), phdos%msqd_dos_atom/Ha_eV))

#else
 MSG_ERROR("netcdf support not enabled")
 ABI_UNUSED((/ncid, phdos%nomega/))
#endif

contains
 integer function vid(vname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vid'
!End of the abilint section

   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
 end function vid

end subroutine phdos_ncwrite
!!***

!----------------------------------------------------------------------

!!****f* m_phonons/mkphbs
!! NAME
!! mkphbs
!!
!! FUNCTION
!! Function to calculate the phonon band structure, from the IFC
!!
!! INPUTS
!! Ifc<ifc_type>=Interatomic force constants
!! crystal<type(crystal_t)> = Info on the crystalline structure.
!! inp= (derived datatype) contains all the input variables
!! ddb<type(ddb_type)>=Object storing the DDB results.
!! asrq0<asrq0_t>=Object for the treatment of the ASR based on the q=0 block found in the DDB file.
!! prefix=Prefix for output files.
!! dielt(3,3)=dielectric tensor
!! comm=MPI communicator
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      ifc_fourq,kpath_free,phonons_ncwrite,phonons_write_gnuplot
!!      phonons_write_phfrq,phonons_write_xmgrace,xmpi_sum_master
!!
!! SOURCE

subroutine mkphbs(Ifc,Crystal,inp,ddb,asrq0,prefix,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkphbs'
 use interfaces_14_hidewrite
 use interfaces_72_response
 use interfaces_77_ddb
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: comm
 character(len=*),intent(in) :: prefix
 type(ifc_type),intent(in) :: Ifc
 type(crystal_t),intent(in) :: Crystal
 type(anaddb_dataset_type),target,intent(in) :: inp
 type(ddb_type),intent(in) :: ddb
 type(asrq0_t),intent(inout) :: asrq0

!Local variables -------------------------
!scalars
 integer,parameter :: master=0
 integer :: unt
 integer :: iphl1,iblok,rftyp, ii,nfineqpath,nsym,natom,ncid,nprocs,my_rank
 integer :: natprj_bs,eivec,enunit,ifcflag
 real(dp) :: freeze_displ
 real(dp) :: cfact
 character(500) :: msg
 character(len=8) :: unitname

!arrays
 integer :: rfphon(4),rfelfd(4),rfstrs(4)
 integer :: nomega, imode, iomega
 integer,allocatable :: ndiv(:)
 real(dp) :: omega, omega_min, gaussmaxarg, gaussfactor, gaussprefactor, xx
 real(dp) :: speedofsound(3)
 real(dp) :: qphnrm(3), qphon(3), qphon_padded(3,3),res(3)
 real(dp) :: d2cart(2,ddb%msize),real_qphon(3)
 real(dp) :: displ(2*3*crystal%natom*3*crystal%natom),eigval(3,crystal%natom)
 real(dp),allocatable :: phfrq(:),eigvec(:,:,:,:,:)
 real(dp),allocatable :: save_phfrq(:,:),save_phdispl_cart(:,:,:,:),save_qpoints(:,:)
 real(dp),allocatable :: weights(:)
 real(dp),allocatable :: dos4bs(:)
 real(dp),allocatable,target :: alloc_path(:,:)
 real(dp),pointer :: fineqpath(:,:)
 type(atprj_type) :: atprj

! *********************************************************************

 ! Only master works for the time being
 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 if (my_rank /= master) return

 nsym = Crystal%nsym; natom = Crystal%natom

 ! Copy parameters from inp (then I will try to remove inp from the API so that I can call mkphbs in eph)
 ifcflag = inp%ifcflag
 natprj_bs = inp%natprj_bs
 freeze_displ = inp%freeze_displ
 eivec = inp%eivec; enunit = inp%enunit
 rftyp=inp%rfmeth

 nullify(fineqpath)
 nfineqpath = inp%nph1l
 fineqpath => inp%qph1l

 if(inp%nph1l==0) then
   if (inp%nqpath==0) then
     return ! if there is nothing to do, return
   else
     ! allow override of nph1l with nqpath if the former is not set
     ! allocate and compute path here and make fineqpath points to it
     ABI_MALLOC(ndiv,(inp%nqpath-1))
     call make_path(inp%nqpath,inp%qpath,Crystal%gmet,'G',inp%ndivsm,ndiv,nfineqpath,alloc_path,std_out)
     ABI_FREE(ndiv)
     fineqpath => alloc_path
   end if
 end if

 write(msg, '(a,(80a),a,a,a,a)' ) ch10,('=',ii=1,80),ch10,ch10,' Treat the first list of vectors ',ch10
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

 if (natprj_bs > 0) call atprj_init(atprj, natom, natprj_bs, inp%iatprj_bs, prefix)

 ABI_MALLOC(phfrq,(3*natom))
 ABI_MALLOC(eigvec,(2,3,natom,3,natom))
 ABI_MALLOC(save_qpoints,(3,nfineqpath))
 ABI_MALLOC(save_phfrq,(3*natom,nfineqpath))
 ABI_MALLOC(save_phdispl_cart,(2,3*natom,3*natom,nfineqpath))
 qphnrm = one

 do iphl1=1,nfineqpath

   ! Initialisation of the phonon wavevector
   qphon(:)=fineqpath(:,iphl1)

   if (inp%nph1l /= 0) qphnrm(1) = inp%qnrml1(iphl1)

   save_qpoints(:,iphl1) = qphon / qphnrm(1)

   ! Generation of the dynamical matrix in cartesian coordinates
   if (ifcflag == 1) then

     ! Get phonon frequencies and displacements in reduced coordinates for this q-point
     !call ifc_fourq(ifc, cryst, save_qpoints(:,iphl1), phfrq, displ, out_eigvec=eigvec)

     ! Get d2cart using the interatomic forces and the
     ! long-range coulomb interaction through Ewald summation
     call gtdyn9(ddb%acell,Ifc%atmfrc,Ifc%dielt,Ifc%dipdip,Ifc%dyewq0,d2cart,Crystal%gmet,ddb%gprim,ddb%mpert,natom,&
&     Ifc%nrpt,qphnrm(1),qphon,Crystal%rmet,ddb%rprim,Ifc%rpt,Ifc%trans,Crystal%ucvol,Ifc%wghatm,Crystal%xred,ifc%zeff)

   else if (ifcflag == 0) then

     !call ddb_diagoq(ddb, crystal, save_qpoints(:,iphl1), asrq0, ifc%symdynmat, rftyp, phfrq, displ, &
     !                out_eigvec=eigvec)

     ! Look for the information in the DDB (no interpolation here!)
     rfphon(1:2)=1; rfelfd(1:2)=0; rfstrs(1:2)=0
     qphon_padded = zero; qphon_padded(:,1) = qphon

     call gtblk9(ddb,iblok,qphon_padded,qphnrm,rfphon,rfelfd,rfstrs,rftyp)

     ! Copy the dynamical matrix in d2cart
     d2cart(:,1:ddb%msize)=ddb%val(:,:,iblok)

     ! Eventually impose the acoustic sum rule based on previously calculated d2asr
     call asrq0_apply(asrq0, natom, ddb%mpert, ddb%msize, crystal%xcart, d2cart)
   end if

   ! Use inp%symdynmat instead of ifc because of ifcflag
   ! Calculation of the eigenvectors and eigenvalues of the dynamical matrix
   call dfpt_phfrq(ddb%amu,displ,d2cart,eigval,eigvec,Crystal%indsym,&
&   ddb%mpert,Crystal%nsym,natom,nsym,Crystal%ntypat,phfrq,qphnrm(1),qphon,&
&   crystal%rprimd,inp%symdynmat,Crystal%symrel,Crystal%symafm,Crystal%typat,Crystal%ucvol)

   if (abs(freeze_displ) > tol10) then
     real_qphon = zero
     if (abs(qphnrm(1)) > tol8) real_qphon = qphon / qphnrm(1)
     call freeze_displ_allmodes(displ, freeze_displ, natom, prefix, phfrq, &
&     real_qphon, crystal%rprimd, Crystal%typat, crystal%xcart, crystal%znucl)
   end if

   ! If requested, output projection of each mode on given atoms
   if (natprj_bs > 0) call atprj_print(atprj, iphl1, phfrq, eigvec)

   ! In case eivec == 4, write output files for band2eps (visualization of phonon band structures)
   if (eivec == 4) then
     call sortph(eigvec,displ,strcat(prefix, "_B2EPS"),natom,phfrq)
   end if

   ! Write the phonon frequencies
   call dfpt_prtph(displ,eivec,enunit,ab_out,natom,phfrq,qphnrm(1),qphon)

   save_phfrq(:,iphl1) = phfrq
   save_phdispl_cart(:,:,:,iphl1) = RESHAPE(displ, [2, 3*natom, 3*natom])

   ! Determine the symmetries of the phonon mode at Gamma
   ! TODO: generalize for other q-point little groups.
   if (sum(abs(qphon)) < DDB_QTOL) then
     call dfpt_symph(ab_out,ddb%acell,eigvec,Crystal%indsym,natom,nsym,phfrq,ddb%rprim,Crystal%symrel)
   end if

   ! if we have an acoustic mode (small q and acoustic type displacements)
   ! extrapolate speed of sound in this direction, and Debye frequency
   call wrap2_pmhalf(qphon, real_qphon, res)
   if (sqrt(real_qphon(1)**2+real_qphon(2)**2+real_qphon(3)**2) < quarter .and. &
&   sqrt(real_qphon(1)**2+real_qphon(2)**2+real_qphon(3)**2) > tol6) then
     call phdos_calc_vsound(eigvec, Crystal%gmet, natom, phfrq, real_qphon, speedofsound)
     if (my_rank == master) call phdos_print_vsound(ab_out, Crystal%ucvol, speedofsound)
   end if

 end do ! iphl1

! calculate dos for the specific q points along the BS calculated
! only Gaussians are possible - no interpolation
 omega_min = minval(save_phfrq(:,:))
 nomega=NINT( (maxval(save_phfrq(:,:))-omega_min) / inp%dosdeltae ) + 1
 nomega=MAX(6,nomega) ! Ensure Simpson integration will be ok

 ABI_MALLOC(dos4bs,(nomega))
 dos4bs = zero
 gaussmaxarg = sqrt(-log(1.d-90))
 gaussprefactor = one/(inp%dossmear*sqrt(two_pi))
 gaussfactor    = one/(sqrt2*inp%dossmear)
 do iphl1=1,nfineqpath
   do imode=1,3*natom
     do iomega=1, nomega
       omega = omega_min + (iomega-1) * inp%dosdeltae
       xx = (omega - save_phfrq(imode,iphl1)) * gaussfactor
       if(abs(xx) < gaussmaxarg) then
         dos4bs(iomega) = dos4bs(iomega) + gaussprefactor*exp(-xx*xx)
       end if
     end do
   end do
 end do


!deallocate sortph array
 call end_sortph()

 if (natprj_bs > 0) call atprj_destroy(atprj)


! WRITE OUT FILES
 if (my_rank == master) then
   ABI_MALLOC(weights, (nfineqpath))
   weights = one

#ifdef HAVE_NETCDF
   NCF_CHECK_MSG(nctk_open_create(ncid, strcat(prefix, "_PHBST.nc"), xmpi_comm_self), "Creating PHBST")
   NCF_CHECK(crystal_ncwrite(Crystal, ncid))
   call phonons_ncwrite(ncid,natom,nfineqpath,save_qpoints,weights,save_phfrq,save_phdispl_cart)

   ! Now treat the second list of vectors (only at the Gamma point, but can include non-analyticities)
   if (inp%nph2l /= 0 .and. inp%ifcflag == 1) then
     call ifc_calcnwrite_nana_terms(ifc, crystal, inp%nph2l, inp%qph2l, inp%qnrml2, ncid)
   end if

   NCF_CHECK(nf90_close(ncid))
#endif

   call phonons_write_phfrq(strcat(prefix, "_PHFRQ"), natom,nfineqpath,save_qpoints,weights,save_phfrq,save_phdispl_cart)

   select case (inp%prtphbands)
   case (0)
     continue

   case (1)
     if (inp%nph1l == 0) then
       call phonons_write_xmgrace(strcat(prefix, "_PHBANDS.agr"), natom, nfineqpath, save_qpoints, save_phfrq, &
          qptbounds=inp%qpath)
     else
       call phonons_write_xmgrace(strcat(prefix, "_PHBANDS.agr"), natom, nfineqpath, save_qpoints, save_phfrq)
     end if

   case (2)
     if (inp%nph1l == 0) then
       call phonons_write_gnuplot(prefix, natom, nfineqpath, save_qpoints, save_phfrq, qptbounds=inp%qpath)
     else
       call phonons_write_gnuplot(prefix, natom, nfineqpath, save_qpoints, save_phfrq)
     end if

   !case (3)
     !call phonons_writeEPS(natom,nfineqpath,Crystal%ntypat,save_qpoints,Crystal%typat, &
     !  weights,save_phfrq,save_phdispl_cart)

   case default
     MSG_WARNING(sjoin("Don't know how to handle prtphbands:", itoa(inp%prtphbands)))
   end select

   ! write out DOS file for q along this path
   cfact=one
   unitname = 'Ha'
   if (open_file('PHBST_partial_DOS',msg,newunit=unt,form="formatted",action="write") /= 0) then
     MSG_ERROR(msg)
   end if
   write(msg,'(3a)')'# ',ch10,'# Partial phonon density of states for q along a band structure path'
   call wrtout(unt,msg,'COLL')
   write(msg,'(6a)')'# ',ch10,'# Energy in ',unitname,', DOS in states/',unitname
   call wrtout(unt,msg,'COLL')
   write(msg,'(a,E20.10,2a,i8)') '# Gaussian method with smearing = ',inp%dossmear*cfact,unitname, &
&        ', nq =', nfineqpath
   call wrtout(unt,msg,'COLL')
   write(msg,'(5a)')'# ',ch10,'# omega     PHDOS ',ch10,'# '
   call wrtout(unt,msg,'COLL')
   do iomega=1,nomega
     omega = omega_min + (iomega-1) * inp%dosdeltae
     write(unt,'(2es17.8)',advance='NO')omega*cfact,dos4bs(iomega)/cfact
     write(unt,*)
   end do
   close(unt)


   ABI_FREE(weights)
 end if

 ABI_FREE(save_qpoints)
 ABI_FREE(save_phfrq)
 ABI_FREE(save_phdispl_cart)
 ABI_FREE(phfrq)
 ABI_FREE(eigvec)
 ABI_FREE(dos4bs)

 if (allocated(alloc_path)) then
   ABI_FREE(alloc_path)
 end if

end subroutine mkphbs
!!***

!!****f* m_phonons/phdos_calc_vsound
!!
!! NAME
!! phdos_calc_vsound
!!
!! FUNCTION
!!  From the frequencies for acoustic modes at small q, estimate speed of sound (which also gives Debye temperature)
!!
!! INPUTS
!! eigvec(2,3*natom,3*natom) = phonon eigenvectors at present q-point
!! gmet(3,3) = metric tensor in reciprocal space.
!! natom = number of atoms in the unit cell
!! phfrq(3*natom) = phonon frequencies at present q-point
!! qphon(3) = phonon q-point
!! ucvol = unit cell volume
!!
!! OUTPUT
!!
!! PARENTS
!!      m_phonons
!!
!! CHILDREN
!!      ifc_fourq,kpath_free,phonons_ncwrite,phonons_write_gnuplot
!!      phonons_write_phfrq,phonons_write_xmgrace,xmpi_sum_master
!!
!! SOURCE

subroutine phdos_calc_vsound(eigvec,gmet,natom,phfrq,qphon, &
&   speedofsound)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phdos_calc_vsound'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalras
 integer, intent(in) :: natom
!arrays
 real(dp), intent(in) :: gmet(3,3),qphon(3)
 real(dp), intent(in) :: phfrq(3*natom),eigvec(2,3*natom,3*natom)

 real(dp), intent(out) :: speedofsound(3)

!Local variables -------------------------
 integer :: iatref,imode, iatom, isacoustic, imode_acoustic
! character(len=500) :: msg
 real(dp) :: qnormcart
 real(dp) :: qtmp(3)

! *********************************************************************

 imode_acoustic = 0
 do imode = 1, 3*natom

!  Check if this mode is acoustic like: scalar product of all displacement vectors are collinear
   isacoustic = 1
!  Find reference atom with non-zero displacement
   do iatom=1,natom
     if(sum(eigvec(:,(iatom-1)*3+1:(iatom-1)*3+3,imode)**2) >tol16)iatref=iatom
   enddo
!  Now compute scalar product, and check they are all positive
   do iatom = 1, natom
     if (sum(eigvec(:,(iatom-1)*3+1:(iatom-1)*3+3, imode)&
&           *eigvec(:,(iatref-1)*3+1:(iatref-1)*3+3, imode)) < tol16 ) isacoustic = 0
   end do
   if (isacoustic == 0) cycle
   imode_acoustic = imode_acoustic + 1

!   write (msg, '(a,I6,a,3F12.4)') ' Found acoustic mode ', imode, ' for |q| in red coord < 0.25 ; q = ', qphon
!   call wrtout(std_out,msg,'COLL')

   qtmp = matmul(gmet, qphon)
   qnormcart = two * pi * sqrt(sum(qphon*qtmp))
   speedofsound(imode_acoustic) = phfrq(imode) / qnormcart
 end do

end subroutine phdos_calc_vsound
!!***

!!****f* m_phonons/phdos_print_vsound
!!
!! NAME
!! phdos_print_vsound
!!
!! FUNCTION
!!  Print out estimate speed of sound and Debye temperature at this (small) q
!!  should only be called by master proc for the hard unit number
!!
!! INPUTS
!! unit=Fortran unit number
!! speedofsound(3)
!!
!! OUTPUT
!!
!! PARENTS
!!      m_phonons
!!
!! CHILDREN
!!      ifc_fourq,kpath_free,phonons_ncwrite,phonons_write_gnuplot
!!      phonons_write_phfrq,phonons_write_xmgrace,xmpi_sum_master
!!
!! SOURCE

subroutine phdos_print_vsound(iunit,ucvol,speedofsound)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phdos_print_vsound'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalras
 integer, intent(in) :: iunit
 real(dp), intent(in) :: ucvol
!arrays
 real(dp), intent(in) :: speedofsound(3)

!Local variables -------------------------
 integer :: imode_acoustic
 character(len=500) :: msg
 real(dp) :: tdebye

! *********************************************************************

 do imode_acoustic = 1, 3
!  from phonon frequency, estimate speed of sound by linear interpolation from Gamma
   write (msg, '(2a,a,E20.10,a,a,F20.5)') &
&   ' Speed of sound for this q and mode:',ch10,&
&   '   in atomic units: ', speedofsound(imode_acoustic), ch10,&
&   '   in units km/s: ', speedofsound(imode_acoustic) * Bohr_Ang * 1.d-13 / Time_Sec
   call wrtout(iunit,msg,'COLL')
   call wrtout(std_out,msg,'COLL')

!  also estimate partial Debye temperature, = energy if this band went to zone edge
   tdebye = speedofsound(imode_acoustic) * pi * (six / pi / ucvol)**(third)
   write (msg, '(2a,a,E20.10,a,a,F20.5)') &
&   ' Partial Debye temperature for this q and mode:',ch10,&
&   '   in atomic units: ', tdebye, ch10,&
&   '   in SI units K  : ', tdebye * Ha_K
   call wrtout(iunit,msg,'COLL')
   call wrtout(iunit,"",'COLL')
   call wrtout(std_out,msg,'COLL')
   call wrtout(std_out,"",'COLL')
 end do

end subroutine phdos_print_vsound
!!***




!----------------------------------------------------------------------

!!****f* m_phonons/phdos_print_msqd
!!
!! NAME
!! phdos_print_msqd
!!
!! FUNCTION
!!  Print out mean square displacement and velocity for each atom (trace and full matrix) as a function of T
!!  see for example https://atztogo.github.io/phonopy/thermal-displacement.html#thermal-displacement
!!
!! INPUTS
!!   PHdos structure
!!
!! OUTPUT
!!   to file only
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      ifc_fourq,kpath_free,phonons_ncwrite,phonons_write_gnuplot
!!      phonons_write_phfrq,phonons_write_xmgrace,xmpi_sum_master
!!
!! SOURCE

subroutine phdos_print_msqd(PHdos, fname, ntemper, tempermin, temperinc)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phdos_print_msqd'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer, intent(in) :: ntemper
 type(phonon_dos_type),intent(in) :: PHdos
 character(len=*),intent(in) :: fname
 real(dp), intent(in) :: tempermin, temperinc

!Local variables -------------------------
 integer :: io, iomin, itemp, iunit, junit, iatom
 real(dp) :: temper
 character(len=500) :: msg
 character(len=fnlen) :: fname_msqd
 character(len=fnlen) :: fname_veloc
!arrays
 real(dp), allocatable :: bose_msqd(:,:), tmp_msqd(:,:), integ_msqd(:,:)
 real(dp), allocatable :: bose_msqv(:,:), tmp_msqv(:,:), integ_msqv(:,:)

! *********************************************************************

 fname_msqd = trim(fname) //"_MSQD_T"
 if (open_file(fname_msqd, msg, newunit=iunit, form="formatted", status="unknown", action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 fname_veloc = trim(fname) // "_MSQV_T"
 if (open_file(fname_veloc, msg, newunit=junit, form="formatted", status="unknown", action="write") /= 0) then
   MSG_ERROR(msg)
 end if

! write a header
   write (msg, '(2a)') '# mean square displacement for each atom as a function of T (bohr^2)'

! NB: this call to wrtout does not seem to work from the eph executable, even in sequential, and whether within or outside a clause for me==master.
!  Do not change this to wrtout without checking extensively.
   !call wrtout(iunit, msg, 'COLL')
   write (iunit, '(a)') trim(msg)

   write (msg, '(2a)') '# mean square velocity for each atom as a function of T (bohr^2/atomic time unit^2)'
   write (junit, '(a)') trim(msg)

   write (msg, '(a,F18.10,a,F18.10,a)') '#  T in Kelvin, from ', tempermin, ' to ', tempermin+(ntemper-1)*temperinc
   !call wrtout(iunit, msg, 'COLL')
   write (iunit, '(a)') trim(msg)
   write (junit, '(a)') trim(msg)

   write (msg, '(2a)') '#    T             |u^2|                u_xx                u_yy                u_zz',&
&                                              '                u_yz                u_xz                u_xy in bohr^2'
   !call wrtout(iunit, msg, 'COLL')
   write (iunit, '(a)') trim(msg)
   write (msg, '(3a)') '#    T             |v^2|                v_xx                v_yy                v_zz',&
&                                              '                v_yz                v_xz                v_xy',&
&                                              ' in bohr^2/atomic time unit^2'
   !call wrtout(iunit, msg, 'COLL')
   write (junit, '(a)') trim(msg)

 ABI_ALLOCATE (tmp_msqd, (PHdos%nomega,9))
 ABI_ALLOCATE (tmp_msqv, (PHdos%nomega,9))
 ABI_ALLOCATE (integ_msqd, (9,ntemper))
 ABI_ALLOCATE (integ_msqv, (9,ntemper))
 ABI_ALLOCATE (bose_msqd, (PHdos%nomega, ntemper))
 ABI_ALLOCATE (bose_msqv, (PHdos%nomega, ntemper))

 do io=1, PHdos%nomega
     if (PHdos%omega(io) >= 1.d-10) exit
 end do
 iomin = io

 ! calculate bose only once for each atom (instead of for each atom)
 bose_msqd = zero
 bose_msqv = zero
 do itemp = 1, ntemper
   temper = tempermin + (itemp-1) * temperinc
   if (abs(temper) < 1.e-10) cycle
   do io = iomin, PHdos%nomega
! NB: factors follow convention in phonopy documentation
!   the 1/sqrt(omega) factor in phonopy is contained in the displacement vector definition
!   bose() is dimensionless
     !bose_msqd(io, itemp) =  (half + one  / ( exp(PHdos%omega(io)/(kb_HaK*temper)) - one )) / PHdos%omega(io)
     !bose_msqv(io, itemp) =  (half + one  / ( exp(PHdos%omega(io)/(kb_HaK*temper)) - one )) * PHdos%omega(io)
     bose_msqd(io, itemp) =  (half + bose_einstein(PHdos%omega(io),kb_HaK*temper)) / PHdos%omega(io)
     bose_msqv(io, itemp) =  (half + bose_einstein(PHdos%omega(io),kb_HaK*temper)) * PHdos%omega(io)
   end do
 end do

 do iatom=1,PHdos%natom
   write (msg, '(a,I8)') '# atom number ', iatom
   !call wrtout(iunit, msg, 'COLL')
   write (iunit, '(a)') trim(msg)
   write (junit, '(a)') trim(msg)

! for each T and each atom, integrate msqd matrix with Bose Einstein factor and output
   integ_msqd = zero
   tmp_msqd = reshape(PHdos%msqd_dos_atom(:,:,:,iatom), (/PHdos%nomega, 9/))
! perform all integrations as matrix multiplication: integ_msqd (idir, itemp) = [tmp_msqd(io,idir)]^T  * bose_msqd(io,itemp)
   call DGEMM('T','N', 9, ntemper, PHdos%nomega, one, tmp_msqd,PHdos%nomega,&
&      bose_msqd, PHdos%nomega, zero, integ_msqd, 9)
! NB: this presumes an equidistant omega grid
   integ_msqd = integ_msqd * (PHdos%omega(2)-PHdos%omega(1)) / PHdos%atom_mass(iatom)

   integ_msqv = zero
   tmp_msqv = reshape(PHdos%msqd_dos_atom(:,:,:,iatom), (/PHdos%nomega, 9/))
! perform all integrations as matrix multiplication: integ_msqv (idir, itemp) = [tmp_msqv(io,idir)]^T  * bose_msqv(io,itemp)
   call DGEMM('T','N', 9, ntemper, PHdos%nomega, one, tmp_msqv,PHdos%nomega,&
&      bose_msqv, PHdos%nomega, zero, integ_msqv, 9)
! NB: this presumes an equidistant omega grid
   integ_msqv = integ_msqv * (PHdos%omega(2)-PHdos%omega(1)) / PHdos%atom_mass(iatom)


! print out stuff
   do itemp = 1, ntemper
     temper = tempermin + (itemp-1) * temperinc
     write (msg, '(F10.2,4x,E22.10,2x,6E22.10)') temper, third*(integ_msqd(1,itemp)+integ_msqd(5,itemp)+integ_msqd(9,itemp)), &
&                      integ_msqd(1,itemp),integ_msqd(5,itemp),integ_msqd(9,itemp), &
&                      integ_msqd(6,itemp),integ_msqd(3,itemp),integ_msqd(2,itemp)
     !call wrtout(iunit, msg, 'COLL')
     write (iunit, '(a)') trim(msg)
     write (msg, '(F10.2,4x,E22.10,2x,6E22.10)') temper, third*(integ_msqv(1,itemp)+integ_msqv(5,itemp)+integ_msqv(9,itemp)), &
&                      integ_msqv(1,itemp),integ_msqv(5,itemp),integ_msqv(9,itemp), &
&                      integ_msqv(6,itemp),integ_msqv(3,itemp),integ_msqv(2,itemp)
     write (junit, '(a)') trim(msg)
   end do ! itemp

   !call wrtout(iunit, msg, 'COLL')
   write (iunit, '(a)') ''
   write (junit, '(a)') ''
 enddo ! iatom

 ABI_DEALLOCATE (tmp_msqd)
 ABI_DEALLOCATE (tmp_msqv)
 ABI_DEALLOCATE (bose_msqd)
 ABI_DEALLOCATE (bose_msqv)
 ABI_DEALLOCATE (integ_msqd)
 ABI_DEALLOCATE (integ_msqv)
 close(iunit)
 close(junit)

end subroutine phdos_print_msqd
!!***

!----------------------------------------------------------------------

!!****f* m_phonons/phonons_ncwrite
!! NAME
!! phonons_ncwrite
!!
!! FUNCTION
!!  Write phonon bandstructure in a netcdf file.
!!
!! INPUTS
!!  ncid =NC file handle
!!  natom=Number of atoms
!!  nqpts=Number of q-points.
!!  qpoints=List of q-points in reduced coordinates
!!  weights(nqpts)= q-point weights
!!  phfreq=Phonon frequencies
!!  phdispl_cart=Phonon displacementent in Cartesian coordinates.
!!
!! NOTES
!!  Input data is in a.u, whereas the netcdf files saves data in eV for frequencies
!!  and Angstrom for the displacements
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      m_phonons
!!
!! CHILDREN
!!      ifc_fourq,kpath_free,phonons_ncwrite,phonons_write_gnuplot
!!      phonons_write_phfrq,phonons_write_xmgrace,xmpi_sum_master
!!
!! SOURCE

subroutine phonons_ncwrite(ncid,natom,nqpts,qpoints,weights,phfreq,phdispl_cart)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phonons_ncwrite'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid,natom,nqpts
!arrays
 real(dp),intent(in) :: qpoints(3,nqpts),weights(nqpts)
 real(dp),intent(in) :: phfreq(3*natom,nqpts),phdispl_cart(2,3*natom,3*natom,nqpts)

!Local variables-------------------------------
!scalars
#ifdef HAVE_NETCDF
 integer :: nphmodes,ncerr

! *************************************************************************

 nphmodes = 3*natom

 NCF_CHECK(nctk_def_basedims(ncid, defmode=.True.))

 ncerr = nctk_def_dims(ncid, [&
   nctkdim_t("number_of_qpoints", nqpts), nctkdim_t('number_of_phonon_modes', nphmodes)])
 NCF_CHECK(ncerr)

! define arrays
 ncerr = nctk_def_arrays(ncid, [&
   nctkarr_t('qpoints', "dp" , 'number_of_reduced_dimensions, number_of_qpoints'),&
   nctkarr_t('qweights',"dp", 'number_of_qpoints'),&
   nctkarr_t('phfreqs',"dp", 'number_of_phonon_modes, number_of_qpoints'),&
   nctkarr_t('phdispl_cart',"dp", 'complex, number_of_phonon_modes, number_of_phonon_modes, number_of_qpoints')])
 NCF_CHECK(ncerr)

!Write variables.
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, vid('qpoints'), qpoints))
 NCF_CHECK(nf90_put_var(ncid, vid('qweights'), weights))
 NCF_CHECK(nf90_put_var(ncid, vid('phfreqs'), phfreq*Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, vid('phdispl_cart'), phdispl_cart*Bohr_Ang))
#endif

contains
 integer function vid(vname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vid'
!End of the abilint section

   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
 end function vid

end subroutine phonons_ncwrite
!!***

!----------------------------------------------------------------------

!!****f* m_phonons/phonons_write_phfrq
!! NAME
!! phonons_write_phfrq
!!
!! FUNCTION
!!  Write phonon bandstructure in a text file. Fixed file name for the moment
!!
!! INPUTS
!!  natom=Number of atoms
!!  nqpts=Number of q-points.
!!  qpoints=List of q-points in reduced coordinates
!!  weights(nqpts)= q-point weights
!!  phfreq=Phonon frequencies
!!  phdispl_cart=Phonon displacementent in Cartesian coordinates.
!!
!! NOTES
!!  Input data is in a.u, output too
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      m_phonons
!!
!! CHILDREN
!!      ifc_fourq,kpath_free,phonons_ncwrite,phonons_write_gnuplot
!!      phonons_write_phfrq,phonons_write_xmgrace,xmpi_sum_master
!!
!! SOURCE

 subroutine phonons_write_phfrq(path,natom,nqpts,qpoints,weights,phfreq,phdispl_cart)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phonons_write_phfrq'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nqpts
 character(len=*),intent(in) :: path
!arrays
 real(dp),intent(in) :: qpoints(3,nqpts),weights(nqpts)
 real(dp),intent(in) :: phfreq(3*natom,nqpts)
 real(dp),intent(in) :: phdispl_cart(2,3*natom,3*natom,nqpts)

!Local variables-------------------------------
!scalars
 integer :: nphmodes, iq, iunit, imod, icomp
 real(dp) :: dummy
 character(len=300) :: formt
 character(len=500) :: msg

! *************************************************************************

 nphmodes = 3*natom

 dummy = qpoints(1,1); dummy = weights(1)

 if (open_file(path, msg, newunit=iunit, form="formatted", status="unknown", action="write") /= 0) then
   MSG_ERROR(msg)
 end if

 write (iunit, '(a)')  '# ABINIT generated phonon band structure file. All in Ha atomic units'
 write (iunit, '(a)')  '# '
 write (iunit, '(a,i0)')  '# number_of_qpoints ', nqpts
 write (iunit, '(a,i0)')  '# number_of_phonon_modes ', nphmodes
 write (iunit, '(a)')  '# '

 write (formt,'(a,i0,a)') "(I5, ", nphmodes, "E20.10)"
 do iq= 1, nqpts
   write (iunit, formt)  iq, phfreq(:,iq)
 end do

 close(iunit)

 if (.False.) then
   if (open_file(strcat(path, "_PHDISPL"), msg, unit=iunit, form="formatted", status="unknown", action="write") /= 0) then
     MSG_ERROR(msg)
   end if

   write (iunit, '(a)')     '# ABINIT generated phonon displacements, along points in PHFRQ file. All in Ha atomic units'
   write (iunit, '(a)')     '# '
   write (iunit, '(a)')     '# displacements in cartesian coordinates, Re and Im parts '
   write (iunit, '(a,i0)')  '# number_of_qpoints ', nqpts
   write (iunit, '(a,i0)')  '# number_of_phonon_modes ', nphmodes
   write (iunit, '(a)')     '# '

   !write (formt,'(a,I3,a)') "( ", nphmodes, "(2E20.10,2x))"
   formt = "(2E20.10,2x)"

   do iq = 1, nqpts
     write (iunit, '(a, i0)') '# iq ', iq
     do imod = 1, nphmodes
       write (iunit, '(a, i0)') '# imode ', imod
       do icomp = 1, nphmodes
         write (iunit, formt, ADVANCE='NO') phdispl_cart(:,icomp,imod,iq)
       end do
       write (iunit, '(a)') ' '
     end do
   end do

   close(iunit)
 end if

end subroutine phonons_write_phfrq
!!***

!----------------------------------------------------------------------

!!****f* m_phonons/phonons_writeEPS
!! NAME
!! phonons_writeEPS
!!
!! FUNCTION
!!  Write phonons bands in EPS format. This routine should be called by a single processor.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      ifc_fourq,kpath_free,phonons_ncwrite,phonons_write_gnuplot
!!      phonons_write_phfrq,phonons_write_xmgrace,xmpi_sum_master
!!
!! SOURCE

subroutine phonons_writeEPS(natom,nqpts,ntypat,qpoints,typat,weights,phfreq,phdispl_cart)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phonons_writeEPS'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nqpts,ntypat
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: qpoints(3,nqpts),weights(nqpts)
 real(dp),intent(in) :: phfreq(3*natom,nqpts)
 real(dp),intent(in) :: phdispl_cart(2,3*natom,3*natom,nqpts)

!Local variables-------------------------------
!scalars
 integer :: cunits,EmaxN,EminN,gradRes,kmaxN,kminN,lastPos,pos,posk
 integer :: iatom,ii,imode,iqpt,jj,nqpt
 integer :: option,unt
 real(dp) :: E,Emax,Emin,deltaE
 real(dp) :: facUnit,norm,renorm
 character(len=500) :: msg
 logical :: set_color = .true.
!array
 complex(dpc) :: displcpx(3*natom,3*natom,nqpts)
 integer,allocatable :: nqptl(:)
 real(dp),allocatable :: phfrq(:),phfrqqm1(:),scale(:)
 real(dp),allocatable :: colorAtom(:,:),color(:,:)
 real(dp),allocatable :: displ(:,:)
 character(len=6),allocatable :: qname(:)

! *********************************************************************


 if (open_file("PHFRQ.eps", msg, unit=unt, form="formatted", status="unknown", action="write") /= 0) then
   MSG_ERROR(msg)
 end if

!Multiplication factor for units (from Hartree to cm-1 or THz)
 if(cunits==1) then
   facUnit=Ha_cmm1
 elseif(cunits==2) then
   facUnit=Ha_THz
 else
 end if

!Boundings of the plot (only the plot and not what is around)
 EminN=6900
 EmaxN=2400
 kminN=2400
 kmaxN=9600

!convert phdispl_cart in cpx array
 displcpx = dcmplx(phdispl_cart(1,:,:,:),phdispl_cart(2,:,:,:))

!Read the input file, and store the information in a long string of characters
!strlen from defs_basis module
 option = 1

!Allocate dynamique variables
 ABI_ALLOCATE(phfrqqm1,(3*natom))
 ABI_ALLOCATE(phfrq,(3*natom))
 ABI_ALLOCATE(color,(3,3*natom))
 ABI_ALLOCATE(qname,(nqpts+1))
 ABI_ALLOCATE(scale,(nqpts))
 ABI_ALLOCATE(nqptl,(nqpts))
 ABI_ALLOCATE(colorAtom,(3,natom))
!colorAtom(1,1:5) : atoms contributing to red (ex : [1 0 0 0 0])
!colorAtom(2,1:5) : atoms contributing to green (ex : [0 1 0 0 0])
!colorAtom(3,1:5) : atoms contributing to blue (ex : [0 0 1 1 1])
 ABI_ALLOCATE(displ,(natom,3*natom))



!TEST_AM TO DO
!Set Values
 if(ntypat /= 3) then
   set_color = .false.
 else
   color = zero
   do ii=1,natom
     if(typat(ii)==1) colorAtom(1,ii) = one
     if(typat(ii)==2) colorAtom(2,ii) = one
     if(typat(ii)==3) colorAtom(3,ii) = one
   end do
 end if

 Emin = -300.0
 Emax =   800.0
 gradRes = 8
 cunits = 1
 qname(:) = "T"
!Read end of input file
 ! read(21,*)
 ! read(21,*) (qname(ii),ii=1,nqpts+1)
 ! read(21,*)
 ! read(21,*) (nqptl(ii),ii=1,nqpts)
 ! read(21,*)
 ! read(21,*) (scale(ii),ii=1,nqpts)
 ! read(21,*)
 ! read(21,*)
 ! read(21,*)
 ! read(21,*) (colorAtom(1,ii),ii=1,natom)
 ! read(21,*)
 ! read(21,*) (colorAtom(2,ii),ii=1,natom)
 ! read(21,*)
 ! read(21,*) (colorAtom(3,ii),ii=1,natom)
!calculate nqpt
 nqpt=0
 do ii=1,nqpts
   nqpt=nqpt+nqptl(ii)
 end do
!compute normalisation factor
 renorm=0
 do ii=1,nqpts
   renorm=renorm+nqptl(ii)*scale(ii)
 end do
 renorm=renorm/nqpt
!Calculate Emin and Emax
 Emin=Emin/FacUnit
 Emax=Emax/FacUnit

!*******************************************************
!Begin to write some comments in the eps file
!This is based to 'xfig'

 write(unt,'(a)') '% !PS-Adobe-2.0 EPSF-2.0'
 write(unt,'(a)') '%%Title: band.ps'
 write(unt,'(a)') '%%BoundingBox: 0 0 581 310'
 write(unt,'(a)') '%%Magnification: 1.0000'

 write(unt,'(a)') '/$F2psDict 200 dict def'
 write(unt,'(a)') '$F2psDict begin'
 write(unt,'(a)') '$F2psDict /mtrx matrix put'
 write(unt,'(a)') '/col-1 {0 setgray} bind def'
 write(unt,'(a)') '/col0 {0.000 0.000 0.000 srgb} bind def'
 write(unt,'(a)') 'end'
 write(unt,'(a)') 'save'
 write(unt,'(a)') 'newpath 0 310 moveto 0 0 lineto 581 0 lineto 581 310 lineto closepath clip newpath'
 write(unt,'(a)') '-36.0 446.0 translate'
 write(unt,'(a)') '1 -1 scale'

 write(unt,'(a)') '/cp {closepath} bind def'
 write(unt,'(a)') '/ef {eofill} bind def'
 write(unt,'(a)') '/gr {grestore} bind def'
 write(unt,'(a)') '/gs {gsave} bind def'
 write(unt,'(a)') '/sa {save} bind def'
 write(unt,'(a)') '/rs {restore} bind def'
 write(unt,'(a)') '/l {lineto} bind def'
 write(unt,'(a)') '/m {moveto} bind def'
 write(unt,'(a)') '/rm {rmoveto} bind def'
 write(unt,'(a)') '/n {newpath} bind def'
 write(unt,'(a)') '/s {stroke} bind def'
 write(unt,'(a)') '/sh {show} bind def'
 write(unt,'(a)') '/slc {setlinecap} bind def'
 write(unt,'(a)') '/slj {setlinejoin} bind def'
 write(unt,'(a)') '/slw {setlinewidth} bind def'
 write(unt,'(a)') '/srgb {setrgbcolor} bind def'
 write(unt,'(a)') '/rot {rotate} bind def'
 write(unt,'(a)') '/sc {scale} bind def'
 write(unt,'(a)') '/sd {setdash} bind def'
 write(unt,'(a)') '/ff {findfont} bind def'
 write(unt,'(a)') '/sf {setfont} bind def'
 write(unt,'(a)') '/scf {scalefont} bind def'
 write(unt,'(a)') '/sw {stringwidth} bind def'
 write(unt,'(a)') '/tr {translate} bind def'
 write(unt,'(a)') '/tnt {dup dup currentrgbcolor'

 write(unt,'(a)') '4 -2 roll dup 1 exch sub 3 -1 roll mul add'
 write(unt,'(a)') '4 -2 roll dup 1 exch sub 3 -1 roll mul add'
 write(unt,'(a)') '4 -2 roll dup 1 exch sub 3 -1 roll mul add srgb}'
 write(unt,'(a)') 'bind def'
 write(unt,'(a)') '/shd {dup dup currentrgbcolor 4 -2 roll mul 4 -2 roll mul'
 write(unt,'(a)') ' 4 -2 roll mul srgb} bind def'
 write(unt,'(a)') '/$F2psBegin {$F2psDict begin /$F2psEnteredState save def} def'
 write(unt,'(a)') '/$F2psEnd {$F2psEnteredState restore end} def'
 write(unt,'(a)') '$F2psBegin'
 write(unt,'(a)') '%%Page: 1 1'
 write(unt,'(a)') '10 setmiterlimit'
 write(unt,'(a)') '0.06000 0.06000 sc'

!****************************************************************
!Begin of the intelligible part of the postcript document

 write(unt,'(a)') '%**************************************'
!****************************************************************
!Draw the box containing the plot
 write(unt,'(a)') '%****Big Box****'
 write(unt,'(a)') '16 slw'
 write(unt,'(a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a)') 'n ', kminN,' ', EmaxN,&
& ' m ', kmaxN,' ', EmaxN, ' l ', &
& kmaxN,' ', EminN, ' l ', kminN,' ', EminN, ' l'
 write(unt,'(a)') 'cp gs col0 s gr'

!****************************************************************
!Write unit on the middle left of the vertical axe
 write(unt,'(a)') '%****Units****'
 if(cunits==1) then
!  1/lambda
   write(unt,'(a)') '/Times-Roman ff 270.00 scf sf'
   write(unt,'(a)') '1425 5650 m'
   write(unt,'(3a)') 'gs 1 -1 sc  90.0 rot (Frequency ',achar(92),'(cm) col0 sh gr'
!  cm-1
   write(unt,'(a)') '/Times-Roman ff 200.00 scf sf'
   write(unt,'(a)') '1325 4030 m'
   write(unt,'(a)') 'gs 1 -1 sc 90.0 rot  (-1) col0 sh gr'
   write(unt,'(a)') '/Times-Roman ff 270.00 scf sf'
   write(unt,'(a)') '1425 3850 m'
   write(unt,'(3a)') 'gs 1 -1 sc  90.0 rot (',achar(92),')) col0 sh gr'
 else
!  Freq
   write(unt,'(a)') '/Times-Roman ff 270.00 scf sf'
   write(unt,'(a)') '825 4850 m'
   write(unt,'(a)') 'gs 1 -1 sc  90.0 rot (Freq) col0 sh gr'
!  THz
   write(unt,'(a)') '/Times-Roman ff 270.00 scf sf'
   write(unt,'(a)') '825 4350 m'
   write(unt,'(a)') 'gs 1 -1 sc 90.0 rot  (THz) col0 sh gr'
 end if
!*****************************************************************
!Write graduation on the vertical axe
 write(unt,'(a)') '%****Vertical graduation****'
 deltaE=(Emax-Emin)/gradRes

!Replacing do loop with real variables with standard g95 do loop
 E=Emin
 do
!  do E=Emin,(Emax-deltaE/2),deltaE
   if (E >= (Emax-deltaE/2)-tol6) exit
   pos=int(((EminN-EmaxN)*E &
&   +EmaxN*Emin -EminN*Emax)/(Emin-Emax))

!  write the value of energy(or frequence)
   write(unt,'(a)') '/Times-Roman ff 270.00 scf sf'
   write(unt,'(i4,a,i4,a)') kminN-800,' ',pos+60,' m'        !-1300 must be CHANGED
!  as a function of the width of E
   write(unt,'(a,i6,a)') 'gs 1 -1 sc (', nint(E*facUnit),') col0 sh gr'

!  write a little bar
   write(unt,'(a,i4,a,i4,a,i4,a,i4,a)') 'n ', kminN,' ',pos ,' m ', kminN+100,' ', pos, ' l'
   write(unt,'(a)') 'gs col0 s gr '

   E = E+deltaE
 end do

!do the same thing for E=Emax (floating point error)
 write(unt,'(a)') '/Times-Roman ff 270.00 scf sf'
 write(unt,'(i4,a,i4,a)') kminN-800,' ',EmaxN+60,' m'        !-1300 must be changed as E
 write(unt,'(a,i6,a)') 'gs 1 -1 sc (', nint(Emax*facUnit),') col0 sh gr'


!draw zero line
 E=0
 pos=int(((EminN-EmaxN)*E &
& +EmaxN*Emin -EminN*Emax)/(Emin-Emax))
 write(unt,'(a,i4,a,i4,a,i4,a,i4,a)') 'n ', kminN,' ',pos ,' m ', kmaxN,' ', pos, ' l'
 write(unt,'(a)') 'gs col0 s gr '


!******************************************************
!draw legend of horizontal axe
!+vertical line

 write(unt,'(a)') '%****Horizontal graduation****'

 lastPos=kminN

 do ii=0,nqpts

   if(ii/=0) then
     posk=int(((kminN-kmaxN)*(nqptl(ii))) &
&     *scale(ii)/renorm/(-nqpt))
   else
     posk=0
   end if

   posk=posk+lastPos
   lastPos=posk

   if(qname(ii+1)=='gamma') then             !GAMMA
     write(unt,'(a)') '/Symbol ff 270.00 scf sf'
     write(unt,'(i4,a,i4,a)') posk-100,' ', 7150, ' m'
     write(unt,'(a)') 'gs 1 -1 sc (G) col0 sh gr'
   elseif(qname(ii+1)=='lambda') then              !LAMBDA
     write(unt,'(a)') '/Symbol ff 270.00 scf sf'
     write(unt,'(i4,a,i4,a)') posk-100,' ', 7150, ' m'
     write(unt,'(a)') 'gs 1 -1 sc (L) col0 sh gr'
   else                                     !autre
     write(unt,'(a)') '/Times-Roman ff 270.00 scf sf'
     write(unt,'(i4,a,i4,a)') posk-100,' ', 7150, ' m'
     write(unt,'(a,a1,a)') 'gs 1 -1 sc (',qname(ii+1),') col0 sh gr'
   end if


!  draw vertical line
   write(unt,'(a,i4,a,i4,a,i4,a,i4,a)') 'n ', posk,' ',EminN ,' m ', posk,' ', EmaxN, ' l'
   write(unt,'(a)') 'gs col0 s gr '


 end do




!***********************************************************
!Write the bands (the most important part actually)

 write(unt,'(a)') '%****Write Bands****'


! read(19,*) (phfrqqm1(ii),ii=1,3*natom)
 jj = 1
 lastPos=kminN
 do iqpt=1,nqpts
!  Copy frequency of the qpoint
   phfrqqm1(:) = phfreq(:,iqpt)
!  Set displacement
   do ii=1,3*natom
     do iatom=1,natom
       displ(iatom,ii) =  real(sqrt(displcpx(3*(iatom-1)+1,ii,iqpt)*   &
           conjg(displcpx(3*(iatom-1)+1,ii,iqpt)) + &
&                displcpx(3*(iatom-1)+2,ii,iqpt)*   &
&          conjg(displcpx(3*(iatom-1)+2,ii,iqpt)) + &
&               displcpx(3*(iatom-1)+3,ii,iqpt)*   &
&          conjg(displcpx(3*(iatom-1)+3,ii,iqpt)) ))
     end do
   end do


   do imode=1,3*natom
!    normalize displ
     norm=0
     do iatom=1,natom
       norm=norm+displ(iatom,imode)
     end do

     do iatom=1,natom
       displ(iatom,imode)=displ(iatom,imode)/norm
     end do

!    Treat color
     color(:,imode)=0
     if(set_color)then
       do ii=1,natom
!        Red
         color(1,imode)=color(1,imode)+displ(ii,imode)*colorAtom(1,ii)
!        Green
         color(2,imode)=color(2,imode)+displ(ii,imode)*colorAtom(2,ii)
!        Blue
         color(3,imode)=color(3,imode)+displ(ii,imode)*colorAtom(3,ii)
       end do
     end if

     pos=int(((EminN-EmaxN)*phfrqqm1(imode) &
&     +EmaxN*Emin -EminN*Emax)/(Emin-Emax))

     posk=int(((kminN-kmaxN)*(iqpt-1) &
&        *scale(jj)/renorm/(-nqpts)))
     posk=posk+lastPos
     write(unt,'(a,i4,a,i4,a)') 'n ',posk,' ',pos,' m'
     pos=int(((EminN-EmaxN)*phfrq(imode) &
&       +EmaxN*Emin -EminN*Emax)/(Emin-Emax))
     posk=int(((kminN-kmaxN)*(iqpt) &
&       *scale(jj)/renorm/(-nqpts)))
     posk=posk+lastPos
     write(unt,'(i4,a,i4,a)') posk,' ',pos,' l gs'

     if(set_color) then     !(in color)
       write(unt,'(f6.3,a,f6.3,a,f6.3,a)') color(1,imode),' ', &
&        color(2,imode),' ',color(3,imode), ' srgb s gr'
     else
       write(unt,'(f6.3,a,f6.3,a,f6.3,a)') 0.0,' ', &
&        0.0,' ',0.0, ' srgb s gr'
     end if
   end do
   lastPos=posk
 end do


!**********************************************************
!Ending the poscript document
 write(unt,'(a)') '$F2psEnd'
 write(unt,'(a)') 'rs'

! *************************************************************************
 close(unt)

end subroutine phonons_writeEPS
!!***

!----------------------------------------------------------------------

!!****f* m_phonons/phonons_write_xmgrace
!! NAME
!! phonons_write_xmgrace
!!
!! FUNCTION
!!  Write phonons bands in Xmgrace format. This routine should be called by a single processor.
!!
!! INPUTS
!!  filename=Filename
!!  natom=Number of atoms
!!  nqpts=Number of q-points
!!  qpts(3,nqpts)=Q-points
!!  phfreqs(3*natom,nqpts)=Phonon frequencies.
!!  [qptbounds(:,:)]=Optional argument giving the extrema of the q-path.
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      m_phonons
!!
!! CHILDREN
!!      ifc_fourq,kpath_free,phonons_ncwrite,phonons_write_gnuplot
!!      phonons_write_phfrq,phonons_write_xmgrace,xmpi_sum_master
!!
!! SOURCE

subroutine phonons_write_xmgrace(filename, natom, nqpts, qpts, phfreqs, qptbounds)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phonons_write_xmgrace'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nqpts
 real(dp),intent(in) :: qpts(3,nqpts),phfreqs(3*natom,nqpts)
 character(len=*),intent(in) :: filename
!arrays
 real(dp),optional,intent(in) :: qptbounds(:,:)

!Local variables-------------------------------
!scalars
 integer :: unt,iq,nu,ii,start,nqbounds
 character(len=500) :: msg
!arrays
 integer :: g0(3)
 integer,allocatable :: bounds2qpt(:)

! *********************************************************************

 nqbounds = 0
 if (present(qptbounds)) then
   if (product(shape(qptbounds)) > 0 ) then
     ! Find correspondence between qptbounds and k-points in ebands.
     nqbounds = size(qptbounds, dim=2)
     ABI_MALLOC(bounds2qpt, (nqbounds))
     bounds2qpt = 1; start = 1
     do ii=1,nqbounds
        do iq=start,nqpts
          if (isamek(qpts(:, iq), qptbounds(:, ii), g0)) then
            bounds2qpt(ii) = iq; start = iq + 1; exit
          end if
        end do
     end do
   end if
 end if

 if (open_file(filename, msg, newunit=unt, form="formatted", action="write") /= 0) then
   MSG_ERROR(msg)
 end if

 write(unt,'(a)') "# Grace project file"
 write(unt,'(a)') "# Generated by Abinit"
 write(unt,'(2(a,i0))') "# natom: ",natom,", nqpt: ",nqpts
 write(unt,'(a)') "# Frequencies are in meV"
 write(unt,'(a)')"# List of q-points and their index (C notation i.e. count from 0)"
 do iq=1,nqpts
   write(unt, "(a)")sjoin("#", itoa(iq-1), ktoa(qpts(:,iq)))
 end do

 write(unt,'(a)') "@page size 792, 612"
 write(unt,'(a)') "@page scroll 5%"
 write(unt,'(a)') "@page inout 5%"
 write(unt,'(a)') "@link page off"
 write(unt,'(a)') "@with g0"
 write(unt,'(a)') "@world xmin 0.00"
 write(unt,'(a,i0)') '@world xmax ',nqpts
 write(unt,'(a,e16.8)') '@world ymin ',minval(phfreqs * Ha_meV)
 write(unt,'(a,e16.8)') '@world ymax ',maxval(phfreqs * Ha_meV)
 write(unt,'(a)') '@default linewidth 1.5'
 write(unt,'(a)') '@xaxis  tick on'
 write(unt,'(a)') '@xaxis  tick major 1'
 write(unt,'(a)') '@xaxis  tick major color 1'
 write(unt,'(a)') '@xaxis  tick major linestyle 3'
 write(unt,'(a)') '@xaxis  tick major grid on'
 write(unt,'(a)') '@xaxis  tick spec type both'
 write(unt,'(a)') '@xaxis  tick major 0, 0'
 if (nqbounds /= 0) then
   write(unt,'(a,i0)') '@xaxis  tick spec ',nqbounds
   do iq=1,nqbounds
     !write(unt,'(a,i0,a,a)') '@xaxis  ticklabel ',iq-1,',', "foo"
     write(unt,'(a,i0,a,i0)') '@xaxis  tick major ',iq-1,' , ',bounds2qpt(iq) - 1
   end do
 end if
 write(unt,'(a)') '@xaxis  ticklabel char size 1.500000'
 write(unt,'(a)') '@yaxis  tick major 10'
 write(unt,'(a)') '@yaxis  label "Phonon Energy [meV]"'
 write(unt,'(a)') '@yaxis  label char size 1.500000'
 write(unt,'(a)') '@yaxis  ticklabel char size 1.500000'
 do nu=1,3*natom
   write(unt,'(a,i0,a)') '@    s',nu-1,' line color 1'
 end do
 do nu=1,3*natom
   write(unt,'(a,i0)') '@target G0.S',nu-1
   write(unt,'(a)') '@type xy'
   do iq=1,nqpts
      write(unt,'(i0,1x,e16.8)') iq-1, phfreqs(nu, iq) * Ha_meV
   end do
   write(unt,'(a)') '&'
 end do

 close(unt)

 if (allocated(bounds2qpt)) then
   ABI_FREE(bounds2qpt)
 end if

end subroutine phonons_write_xmgrace
!!***

!----------------------------------------------------------------------

!!****f* m_phonons/phonons_write_gnuplot
!! NAME
!! phonons_write_gnuplot
!!
!! FUNCTION
!!  Write phonons bands in gnuplot format. This routine should be called by a single processor.
!!
!! INPUTS
!!  prefix=prefix for files (.data, .gnuplot)
!!  natom=Number of atoms
!!  nqpts=Number of q-points
!!  qpts(3,nqpts)=Q-points
!!  phfreqs(3*natom,nqpts)=Phonon frequencies.
!!  [qptbounds(:,:)]=Optional argument giving the extrema of the q-path.
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      m_phonons
!!
!! CHILDREN
!!      ifc_fourq,kpath_free,phonons_ncwrite,phonons_write_gnuplot
!!      phonons_write_phfrq,phonons_write_xmgrace,xmpi_sum_master
!!
!! SOURCE

subroutine phonons_write_gnuplot(prefix, natom, nqpts, qpts, phfreqs, qptbounds)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phonons_write_gnuplot'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nqpts
 real(dp),intent(in) :: qpts(3,nqpts),phfreqs(3*natom,nqpts)
 character(len=*),intent(in) :: prefix
!arrays
 real(dp),optional,intent(in) :: qptbounds(:,:)

!Local variables-------------------------------
!scalars
 integer :: unt,iq,ii,start,nqbounds,gpl_unt
 character(len=500) :: msg,fmt
 character(len=fnlen) :: datafile,basefile
!arrays
 integer :: g0(3)
 integer,allocatable :: bounds2qpt(:)

! *********************************************************************

 nqbounds = 0
 if (present(qptbounds)) then
   if (product(shape(qptbounds)) > 0 ) then
     ! Find correspondence between qptbounds and k-points in ebands.
     nqbounds = size(qptbounds, dim=2)
     ABI_MALLOC(bounds2qpt, (nqbounds))
     bounds2qpt = 1; start = 1
     do ii=1,nqbounds
        do iq=start,nqpts
          if (isamek(qpts(:, iq), qptbounds(:, ii), g0)) then
            bounds2qpt(ii) = iq; start = iq + 1; exit
          end if
        end do
     end do
   end if
 end if

 datafile = strcat(prefix, "_PHBANDS.data")
 if (open_file(datafile, msg, newunit=unt, form="formatted", action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 if (open_file(strcat(prefix, "_PHBANDS.gnuplot"), msg, newunit=gpl_unt, form="formatted", action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 basefile = basename(datafile)

 write(unt,'(a)') "# Phonon band structure data file"
 write(unt,'(a)') "# Generated by Abinit"
 write(unt,'(2(a,i0))') "# natom: ",natom,", nqpt: ",nqpts
 write(unt,'(a)') "# Frequencies are in meV"
 write(unt,'(a)')"# List of q-points and their index (C notation i.e. count from 0)"
 do iq=1,nqpts
   write(unt, "(a)")sjoin("#", itoa(iq-1), ktoa(qpts(:,iq)))
 end do

 fmt = sjoin("(i0,1x,", itoa(3*natom), "(es16.8,1x))")
 write(unt,'(a)')"# [kpt-index, mode_1, mode_2 ...]"
 do iq=1,nqpts
   write(unt, fmt) iq-1, phfreqs(:, iq) * Ha_meV
 end do

 ! gnuplot script file
  write(gpl_unt,'(a)') '# File to plot electron bandstructure with gnuplot'
  !write(gpl_unt,'(a)') "#set terminal postscript eps enhanced color font 'Times-Roman,26' lw 2"
  write(gpl_unt,'(a)') '#use the next lines to make a nice figure for a paper'
  write(gpl_unt,'(a)') '#set term postscript enhanced eps color lw 0.5 dl 0.5'
  write(gpl_unt,'(a)') '#set pointsize 0.275'
  write(gpl_unt,'(a)') 'set palette defined ( 0 "blue", 3 "green", 6 "yellow", 10 "red" )'
  write(gpl_unt,'(a)') 'unset key'
  write(gpl_unt,'(a)') '# can make pointsize smaller (~0.5). Too small and nothing is printed'
  write(gpl_unt,'(a)') 'set pointsize 0.8'
  write(gpl_unt,'(a)') 'set view 0,0'
  write(gpl_unt,'(a,i0,a)') 'set xrange [0:',nqpts-1,']'
  write(gpl_unt,'(2(a,es16.8),a)')&
    'set yrange [',minval(phfreqs * Ha_meV),':',maxval(phfreqs * Ha_meV),']'
  write(gpl_unt,'(a)') 'set xlabel "Momentum"'
  write(gpl_unt,'(a)') 'set ylabel "Energy [meV]"'
  write(gpl_unt,'(a)') strcat('set title "', replace(basefile, "_", "\\_"),'"')
  if (nqbounds == 0) then
     write(gpl_unt,'(a)') 'set grid xtics'
  else
    write(gpl_unt,"(a)")"# Add vertical lines in correspondence of high-symmetry points."
    write(gpl_unt,'(a)') 'unset xtics'
    do ii=1,nqbounds
      write(gpl_unt,"(a,2(i0,a))") &
        "set arrow from ",bounds2qpt(ii)-1,",graph(0,0) to ",bounds2qpt(ii)-1,",graph(1,1) nohead"
      !write(gpl_unt,"(a)")sjoin("set xtics add ('kname'", itoa(bounds2kpt(ii)-1), ")")
    end do
  end if
  write(gpl_unt,"(a)")sjoin("nbranch =", itoa(3*natom))
  write(gpl_unt,"(a)")strcat('plot for [i=2:nbranch] "', basefile, '" u 1:i every :1 with lines linetype -1')
  write(gpl_unt,"(a)")"pause -1"

 close(unt)
 close(gpl_unt)

 if (allocated(bounds2qpt)) then
   ABI_FREE(bounds2qpt)
 end if

end subroutine phonons_write_gnuplot
!!***

!!****f* m_phonons/ifc_mkphbs
!! NAME
!! ifc_mkphbs
!!
!! FUNCTION
!! Compute the phonon band structure from the IFC and write data to file(s)
!!
!! INPUTS
!! ifc<ifc_type>=Interatomic force constants
!! cryst<crystal_t> = Info on the crystalline structure.
!! dtset=<datasets_type>: input: all input variables initialized from the input file.
!! prefix=Prefix for output files.
!! comm=MPI communicator
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      eph
!!
!! CHILDREN
!!      ifc_fourq,kpath_free,phonons_ncwrite,phonons_write_gnuplot
!!      phonons_write_phfrq,phonons_write_xmgrace,xmpi_sum_master
!!
!! SOURCE

subroutine ifc_mkphbs(ifc, cryst, dtset, prefix, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_mkphbs'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: comm
 character(len=*),intent(in) :: prefix
 type(ifc_type),intent(in) :: ifc
 type(crystal_t),intent(in) :: cryst
 type(dataset_type),intent(in) :: dtset

!Local variables -------------------------
!scalars
 integer,parameter :: master=0
 integer :: iqpt,nqpts,natom,ncid,nprocs,my_rank,ierr
 !character(500) :: msg
 type(kpath_t) :: qpath
!arrays
 real(dp),allocatable :: eigvec(:,:,:,:,:),phfrqs(:,:),phdispl_cart(:,:,:,:),weights(:)

! *********************************************************************

 if (dtset%prtphbands == 0) return

 if (dtset%ph_nqpath <= 0 .or. dtset%ph_ndivsm <= 0) then
   MSG_WARNING("ph_nqpath <=0 or ph_ndivsm <= 0, phonon bands won't be produced. returning")
   return
 end if

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 natom = cryst%natom

 qpath = kpath_new(dtset%ph_qpath(:,1:dtset%ph_nqpath), cryst%gprimd, dtset%ph_ndivsm)
 nqpts = qpath%npts

 ABI_CALLOC(phfrqs, (3*natom,nqpts))
 ABI_CALLOC(phdispl_cart, (2,3*natom,3*natom,nqpts))
 ABI_CALLOC(eigvec, (2,3,natom,3,natom))

 do iqpt=1,nqpts
   if (mod(iqpt, nprocs) /= my_rank) cycle ! mpi-parallelism
   ! Get phonon frequencies and displacements in cartesian coordinates for this q-point
   call ifc_fourq(ifc, cryst, qpath%points(:,iqpt), phfrqs(:,iqpt), phdispl_cart(:,:,:,iqpt), out_eigvec=eigvec)
 end do

 call xmpi_sum_master(phfrqs, master, comm, ierr)
 call xmpi_sum_master(phdispl_cart, master, comm, ierr)

 if (my_rank == master) then
   ABI_MALLOC(weights, (nqpts))
   weights = one

#ifdef HAVE_NETCDF
   NCF_CHECK_MSG(nctk_open_create(ncid, strcat(prefix, "_PHBST.nc"), xmpi_comm_self), "Creating PHBST")
   NCF_CHECK(crystal_ncwrite(cryst, ncid))
   call phonons_ncwrite(ncid,natom,nqpts, qpath%points, weights, phfrqs, phdispl_cart)
   NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('atomic_mass_units', "dp", "number_of_atom_species")],defmode=.True.))
   NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'atomic_mass_units'), ifc%amu))
   NCF_CHECK(nf90_close(ncid))
#endif

   select case (dtset%prtphbands)
   case (1)
     call phonons_write_xmgrace(strcat(prefix, "_PHBANDS.agr"), natom, nqpts, qpath%points, phfrqs, qptbounds=qpath%bounds)
   case (2)
     call phonons_write_gnuplot(prefix, natom, nqpts, qpath%points, phfrqs, qptbounds=qpath%bounds)
   case (3)
     call phonons_write_phfrq(strcat(prefix, "_PHFRQ"), natom, nqpts, qpath%points, weights, phfrqs, phdispl_cart)
   case default
     MSG_WARNING(sjoin("Unsupported value for prtphbands:", itoa(dtset%prtphbands)))
   end select

   ABI_FREE(weights)
 end if ! master

 ABI_FREE(phfrqs)
 ABI_FREE(phdispl_cart)
 ABI_FREE(eigvec)

 call kpath_free(qpath)

end subroutine ifc_mkphbs
!!***

end module m_phonons
!!***
