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
!! Copyright (C) 1999-2019 ABINIT group (XG, MG, MJV, GMR)
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
 use m_abicore
 use m_htetra
 use m_numeric_tools
 use m_crystal
 use m_nctk
 use iso_c_binding
 use m_atprj
 use m_sortph
 use m_ddb
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_supercell
 use m_dtset

 use m_fstrings,        only : itoa, ftoa, sjoin, ktoa, strcat, basename, replace
 use m_symtk,           only : matr3inv
 use m_time,            only : cwtime, cwtime_report
 use m_io_tools,        only : open_file
 use m_geometry,        only : mkrdim, symredcart, normv
 use m_dynmat,          only : gtdyn9, dfpt_phfrq, dfpt_prtph
 use m_bz_mesh,         only : isamek, make_path, kpath_t, kpath_new
 use m_ifc,             only : ifc_type
 use m_anaddb_dataset,  only : anaddb_dataset_type
 use m_kpts,            only : kpts_ibz_from_kptrlatt, get_full_kgrid
 use m_special_funcs,   only : bose_einstein
 use m_sort,            only : sort_dp

 implicit none

 private

 public :: mkphbs                        ! Compute phonon band structure
 public :: phonons_write_xmgrace         ! Write phonons bands in Xmgrace format.
 public :: phonons_write_gnuplot         ! Write phonons bands in gnuplot format.
 public :: ifc_mkphbs                    ! Compute the phonon band structure from the IFC and write data to file(s)
 public :: dfpt_symph                    ! Determine the symmetry character of the different phonon modes at Gamma

 public :: zacharias_supercell_make
 public :: zacharias_supercell_print
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
   ! allows one to calculate Debye Waller factors by integration with 1/omega
   ! and the Bose Einstein factor

 contains

   procedure :: print => phdos_print
   procedure :: print_debye => phdos_print_debye
   procedure :: print_msqd => phdos_print_msqd
   procedure :: print_thermo => phdos_print_thermo
   procedure :: free => phdos_free
   procedure :: ncwrite => phdos_ncwrite

 end type phonon_dos_type

 public :: mkphdos
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

!Arguments ------------------------------------
 character(len=*),intent(in) :: fname
 class(phonon_dos_type),intent(in) :: PHdos

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

 ! Open external file and write results
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

!Arguments ------------------------------------
 real(dp), intent(in) :: ucvol
 class(phonon_dos_type),intent(in) :: PHdos

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
!! NAME
!! phdos_print_thermo
!!
!! FUNCTION
!! Print out global thermodynamic quantities based on DOS
!! Only master node should call this routine.
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
!!
!! SOURCE

subroutine phdos_print_thermo(PHdos, fname, ntemper, tempermin, temperinc)

!Arguments ------------------------------------
 integer, intent(in) :: ntemper
 real(dp), intent(in) :: tempermin, temperinc
 class(phonon_dos_type),intent(in) :: PHdos
 character(len=*),intent(in) :: fname

!Local variables-------------------------------
 integer :: iomega, itemper, tunt
 character(len=500) :: msg
 real(dp) :: wover2t, ln2shx, cothx, invsinh2
 real(dp) :: tmp, domega
!arrays
 real(dp), allocatable :: free(:), energy(:), entropy(:), spheat(:),wme(:)

! *********************************************************************

 ! Allocate and put zeroes for F, E, S, Cv
 ABI_CALLOC(free,    (ntemper))
 ABI_CALLOC(energy,  (ntemper))
 ABI_CALLOC(entropy, (ntemper))
 ABI_CALLOC(spheat,  (ntemper))
 ABI_CALLOC(wme,     (ntemper))

 ! open THERMO file
 if (open_file(fname, msg, newunit=tunt, form="formatted", action="write") /= 0) then
   MSG_ERROR(msg)
 end if

 write(msg, '(a,a,a)' )&
&  ' phdos_print_thermo: thermodynamic functions calculated from prtdos DOS (not histogram)',ch10,&
&  '     see THERMO output file ...'
 call wrtout(std_out,msg,'COLL')

 ! print header
 write(tunt,'(a,a)') ch10,&
&  ' # At  T     F(J/mol-c)     E(J/mol-c)     S(J/(mol-c.K)) C(J/(mol-c.K)) Omega_mean(cm-1) from prtdos DOS'
 write(tunt, "(a)")' # (A mol-c is the abbreviation of a mole-cell, that is, the'
 write(tunt, "(a)")' #  number of Avogadro times the atoms in a unit cell)'

 domega = phdos%omega_step

 do itemper=1,ntemper
   ! The temperature (tmp) is given in Ha
   tmp=(tempermin+temperinc*dble(itemper-1))*kb_HaK

   do iomega=1,PHdos%nomega
     if (abs(PHdos%phdos(iomega)) < 1.e-200_dp) cycle

     ! wover2t= hbar*w / 2kT dimensionless
     wover2t = zero;     if(tmp > tol14) wover2t=PHdos%omega(iomega)*half/tmp
     ! should not be much of a problem for the log, but still put a check.
     ln2shx=zero;        if (wover2t > tol16 .and. wover2t < 100.0_dp) ln2shx=log(two * sinh(wover2t))
     cothx=zero;         if (wover2t > tol16) cothx=one/tanh(wover2t)
     invsinh2=zero;      if (wover2t > tol16 .and. wover2t < 100.0_dp) invsinh2=one/sinh(wover2t)**2

     ! This matches the equations published in Lee & Gonze, PRB 51, 8610 (1995) [[cite:Lee1995]]
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
   write(tunt,'(es11.3,5es15.7)') tmp/kb_HaK,&
&    Ha_J*Avogadro*free(itemper),&
&    Ha_J*Avogadro*energy(itemper),&
&    Ha_J*Avogadro*kb_HaK*entropy(itemper),&
&    Ha_J*Avogadro*kb_HaK*spheat(itemper),&
&    wme(itemper)*Ha_cmm1
 end do ! itemper

 close(tunt)

 ABI_FREE(free)
 ABI_FREE(energy)
 ABI_FREE(entropy)
 ABI_FREE(spheat)
 ABI_FREE(wme)

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
!!
!! SOURCE

subroutine phdos_free(PHdos)

!Arguments -------------------------------
 class(phonon_dos_type),intent(inout) ::PHdos

! *************************************************************************

 !@phonon_dos_type
 ABI_SFREE(PHdos%atom_mass)
 ABI_SFREE(PHdos%omega)
 ABI_SFREE(PHdos%phdos)
 ABI_SFREE(PHdos%phdos_int)
 ABI_SFREE(PHdos%pjdos)
 ABI_SFREE(PHdos%pjdos_int)
 ABI_SFREE(PHdos%pjdos_type)
 ABI_SFREE(PHdos%pjdos_type_int)
 ABI_SFREE(PHdos%pjdos_rc_type)
 ABI_SFREE(PHdos%msqd_dos_atom)

end subroutine phdos_free
!!***

!--------------------------------------------------------------------------

!!****f* m_phonons/phdos_init
!!
!! NAME
!! phdos_init
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine phdos_init(phdos, crystal, ifc, dosdeltae, dossmear, wminmax, prtdos)

 ! Arguments ------------------------------------------------------
 type(crystal_t),intent(in) :: crystal
 type(phonon_dos_type),intent(out) :: phdos
 type(ifc_type),intent(in) :: ifc
 integer,intent(in) :: prtdos
 real(dp),intent(in) :: dosdeltae,dossmear
 real(dp),intent(in) :: wminmax(2)

 !Local variables -------------------------------------------------
 integer :: io

 phdos%ntypat     = crystal%ntypat
 phdos%natom      = crystal%natom
 phdos%prtdos     = prtdos
 phdos%dossmear   = dossmear
 phdos%omega_step = dosdeltae
 ! Use values stored in ifc (obtained with ab-initio q-mesh + pad)
 if (wminmax(2) > wminmax(1)) then
   phdos%omega_min = wminmax(1)
   phdos%omega_max = wminmax(2)
 else
   phdos%omega_min = ifc%omega_minmax(1)
   phdos%omega_max = ifc%omega_minmax(2)
 end if
 ! Must be consistent with mesh computed in tetra routines!
 phdos%nomega = nint((phdos%omega_max - phdos%omega_min) / phdos%omega_step) + 1
 ! Ensure Simpson integration will be ok
 phdos%nomega = max(6, phdos%nomega)

 ! Build frequency mesh.
 ABI_MALLOC(phdos%omega, (phdos%nomega))
 do io=1,phdos%nomega
   phdos%omega(io) = phdos%omega_min + phdos%omega_step * (io - 1)
 end do
 phdos%omega_min = phdos%omega(1)
 phdos%omega_max = phdos%omega(phdos%nomega)

 ! Allocate arrays that depend on nomega and set them to zero.
 ABI_CALLOC(phdos%phdos, (phdos%nomega))
 ABI_CALLOC(phdos%phdos_int, (phdos%nomega))
 ABI_CALLOC(phdos%pjdos, (phdos%nomega, 3, crystal%natom))
 ABI_CALLOC(phdos%pjdos_int, (phdos%nomega, 3, crystal%natom))
 ABI_CALLOC(phdos%msqd_dos_atom, (phdos%nomega, 3, 3, crystal%natom))
 ABI_MALLOC(phdos%atom_mass, (crystal%natom))
 phdos%atom_mass = crystal%amu(crystal%typat(:)) * amu_emass

end subroutine phdos_init
!!***

!---------------------------------------------------------------

!!****f* m_phonons/mkphdos
!!
!! NAME
!! mkphdos
!!
!! FUNCTION
!! Calculate the phonon density of states as well as
!! the contributions associated to the different types of atoms in the unit cell.
!! Two methods are implemented: gaussian method and linear interpolation based on tetrahedrons.
!!
!! INPUTS
!! ifc<ifc_type>=Interatomic force constants
!! crystal<crystal_t>=Info on the crystalline structure.
!! prtdos=1 for gaussian method, 2 for tetrahedra.
!! dosdeltae=Step of frequency mesh.
!! dossmear=Gaussian broadening, used if prtdos==1.
!! dos_ngqpt(3)=Divisions of the q-mesh used for computing the DOS
!! nqshift=Number of shifts in Q-mesh
!! dos_qshift(3, nqshift)=Shift of the q-mesh.
!! prefix=Prefix for output files.
!! comm=MPI communicator.
!!
!! OUTPUT
!! phdos<phonon_dos_type>=Container with phonon DOS, IDOS and atom-projected DOS.
!! count_wminmax(2)=Number of (interpolated) phonon frequencies that are outside
!!   input range (see wminmax). Client code can use count_wminmax and wminmax to
!!   enlarge the mesh and call the routine again to recompute the DOS
!!   if all frequencies should be included.
!!
!! SIDE EFFECTS
!! wminmax(2)=
!!   In input: min and max value of frequency mesh. Used only if minmax(2) > minmax(1)
!!    else values are taken from ifc%omega_minmax (computed from ab-initio mesh + pad)
!!   In output: min and max frequency obtained after interpolating the IFCs on the dense q-mesh dos_ngqpt
!!
!! PARENTS
!!      anaddb,eph,m_tdep_phdos
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkphdos(phdos, crystal, ifc, prtdos, dosdeltae_in, dossmear, dos_ngqpt, nqshft, dos_qshift, prefix, &
                   wminmax, count_wminmax, comm)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: prtdos,nqshft,comm
 real(dp),intent(in) :: dosdeltae_in,dossmear
 character(len=*),intent(in) ::  prefix
 type(crystal_t),intent(in) :: crystal
 type(ifc_type),intent(in) :: ifc
 type(phonon_dos_type),intent(out) :: phdos
!arrays
 integer,intent(in) :: dos_ngqpt(3)
 integer,intent(out) :: count_wminmax(2)
 real(dp),intent(in) :: dos_qshift(3,nqshft)
 real(dp),intent(inout) :: wminmax(2)

!Local variables -------------------------
!scalars
 integer,parameter :: brav1=1,chksymbreak0=0,bcorr0=0,qptopt1=1,master=0
 integer :: iat,jat,idir,imode,io,iq_ibz,itype,nkpt_fullbz
 integer :: nqbz,ierr,natom,nomega,jdir, isym, nprocs, my_rank, ncid
 logical :: refine_dosdeltae
 real(dp),parameter :: max_occ1=one, gaussmaxarg = sqrt(-log(1.d-90)), max_smallq = 0.0625_dp
 real(dp) :: nsmallq,gaussfactor,gaussprefactor,normq,debyefreq,rtmp
 real(dp) :: cpu, wall, gflops
 real(dp) :: dosdeltae, phdos_int
 character(len=500) :: msg
 character(len=80) :: errstr
 type(htetra_t) :: htetraq
!arrays
 integer :: in_qptrlatt(3,3),new_qptrlatt(3,3)
 integer,allocatable :: bz2ibz_smap(:,:), bz2ibz(:)
 real(dp) :: speedofsound(3),speedofsound_(3)
 real(dp) :: displ(2*3*Crystal%natom*3*Crystal%natom)
 real(dp) :: eigvec(2,3,Crystal%natom,3*Crystal%natom),phfrq(3*Crystal%natom)
 real(dp) :: qlatt(3,3),rlatt(3,3)
 real(dp) :: msqd_atom_tmp(3,3),temp_33(3,3)
 real(dp) :: symcart(3,3,crystal%nsym)
 real(dp) :: syme2_xyza(3, crystal%natom)
 real(dp),allocatable :: full_eigvec(:,:,:,:,:),full_phfrq(:,:),new_shiftq(:,:)
 real(dp),allocatable :: qbz(:,:),qibz(:,:),tmp_phfrq(:)
 real(dp),allocatable :: wtq_ibz(:),xvals(:), gvals_wtq(:), wdt(:,:)
 real(dp),allocatable :: energies(:)

! *********************************************************************

 DBG_ENTER("COLL")

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 ! Consistency check.
 if (all(prtdos /= [1, 2])) then
   MSG_BUG(sjoin('prtdos should be 1 or 2, but received', itoa(prtdos)))
 end if
 dosdeltae = dosdeltae_in
 refine_dosdeltae = .false.
 if (dosdeltae <= zero) then
   dosdeltae = -dosdeltae
   refine_dosdeltae = .true.
 end if
 if (prtdos == 1 .and. dossmear <= zero) then
   MSG_BUG(sjoin('dossmear should be positive but received', ftoa(dossmear)))
 end if

 call cwtime(cpu, wall, gflops, "start")

 ! Get symmetries in cartesian coordinates
 do isym = 1, crystal%nsym
   call symredcart(crystal%rprimd,crystal%gprimd,symcart(:,:,isym),crystal%symrel(:,:,isym))
 end do

 natom = crystal%natom
 call phdos_init(phdos, crystal, ifc, dosdeltae, dossmear, wminmax, prtdos)
 nomega = phdos%nomega

 ABI_MALLOC(gvals_wtq, (nomega))
 ABI_MALLOC(xvals, (nomega))

 ! Parameters defining the gaussian approximant.
 if (prtdos == 1) then
   ! TODO: use gaussian and update reference files.
   gaussprefactor = one / (dossmear * sqrt(two_pi))
   gaussfactor = one / (sqrt2 * dossmear)
   write(msg, '(4a,f8.5,2a,f8.5)') ch10, &
    ' mkphdos: calculating phonon DOS using gaussian method:',ch10, &
    '    gaussian smearing [meV] = ',dossmear*Ha_meV,ch10, &
    '    frequency step    [meV] = ',phdos%omega_step*Ha_meV
 else if (prtdos == 2) then
   write(msg,'(2a)')ch10,' mkphdos: calculating phonon DOS using tetrahedron method'
 end if
 call wrtout(std_out, msg, 'COLL')

 ! TODO
 ! 1) fix bug in tetra if degenerate and update ref files
 ! 2) nshift > 1?

 ! This call will set %nqibz and IBZ and BZ arrays
 in_qptrlatt = 0; in_qptrlatt(1, 1) = dos_ngqpt(1); in_qptrlatt(2, 2) = dos_ngqpt(2); in_qptrlatt(3, 3) = dos_ngqpt(3)

 call kpts_ibz_from_kptrlatt(crystal, in_qptrlatt, qptopt1, nqshft, dos_qshift, &
   phdos%nqibz, qibz, wtq_ibz, nqbz, qbz, new_kptrlatt=new_qptrlatt, new_shiftk=new_shiftq, bz2ibz=bz2ibz_smap)
 call cwtime_report(" kpts_ibz_from_kptrlatt", cpu, wall, gflops)

 if (prtdos == 2) then
   ! Prepare tetrahedron method including workspace arrays.
   ! Convert kptrlatt to double and invert, qlatt here refer to the shortest qpt vectors
   rlatt = new_qptrlatt; call matr3inv(rlatt, qlatt)

   nkpt_fullbz = nqbz
   ABI_MALLOC(bz2ibz, (nkpt_fullbz))
   bz2ibz = bz2ibz_smap(1,:)

   call htetra_init(htetraq, bz2ibz, crystal%gprimd, qlatt, qbz, nqbz, qibz, phdos%nqibz, ierr, errstr, comm)
   call cwtime_report(" init_tetra", cpu, wall, gflops)
   ABI_CHECK(ierr == 0, errstr)

   !ABI_FREE(kpt_fullbz)

   ! Allocate arrays used to store the entire spectrum, Required to calculate tetra weights.
   ! this may change in the future if Matteo refactorizes the tetra weights as sums over k instead of sums over bands
   ABI_CALLOC(full_phfrq, (3*natom, phdos%nqibz))
   ABI_MALLOC_OR_DIE(full_eigvec, (2, 3, natom, 3*natom, phdos%nqibz), ierr)
   full_eigvec = zero
 end if ! tetra
 ABI_SFREE(bz2ibz)
 ABI_SFREE(bz2ibz_smap)
 ABI_FREE(new_shiftq)

 ! MPI Sum over irreducible q-points then sync the following integrals:
 !   speedofsound, nsmallq
 !   wminmax and count_wminmax
 !   if gauss: %phdos, %msqd_dos_atom
 !   if tetra: full_phfrq, full_eigvec, %phdos_int

 nsmallq = zero; speedofsound = zero
 wminmax = [huge(one), -huge(one)]; count_wminmax = 0
 do iq_ibz=1,phdos%nqibz
   if (mod(iq_ibz, nprocs) /= my_rank) cycle ! mpi-parallelism

   ! Fourier interpolation (keep track of min/max to decide if initial mesh was large enough)
   call ifc%fourq(crystal, qibz(:,iq_ibz), phfrq, displ, out_eigvec=eigvec)
   wminmax(1) = min(wminmax(1), minval(phfrq))
   if (wminmax(1) < phdos%omega(1)) count_wminmax(1) = count_wminmax(1) + 1
   wminmax(2) = max(wminmax(2), maxval(phfrq))
   if (wminmax(2) > phdos%omega(nomega)) count_wminmax(2) = count_wminmax(2) + 1

   normq = sum(qibz(:,iq_ibz) ** 2)
   if (normq < max_smallq .and. normq > tol6) then
     call phdos_calc_vsound(eigvec, crystal%gmet, natom, phfrq, qibz(:,iq_ibz), speedofsound_)
     speedofsound = speedofsound + speedofsound_ * wtq_ibz(iq_ibz)
     nsmallq = nsmallq + wtq_ibz(iq_ibz)
   end if

   select case (prtdos)
   case (1)
     do imode=1,3*natom
       ! Precompute \delta(w - w_{qnu}) * weight(q)
       xvals = (phdos%omega(:) - phfrq(imode)) * gaussfactor
       where (abs(xvals) < gaussmaxarg)
         gvals_wtq = gaussprefactor * exp(-xvals*xvals) * wtq_ibz(iq_ibz)
       elsewhere
         gvals_wtq = zero
       end where

       ! Accumulate PHDOS
       phdos%phdos(:) = phdos%phdos(:) + gvals_wtq

       ! Rotate e(q) to get e(Sq) to account for symmetrical q-points in BZ.
       ! eigenvectors indeed are not invariant under rotation. See e.g. Eq 39-40 of PhysRevB.76.165108 [[cite:Giustino2007]].
       ! In principle there's a phase due to nonsymmorphic translations
       ! but we here need |e(Sq)_iatom|**2
       syme2_xyza = zero
       do iat=1,natom
         do isym=1, crystal%nsym
           jat = crystal%indsym(4,isym,iat)
           syme2_xyza(:,jat) = syme2_xyza(:,jat) + &
             matmul(symcart(:,:,isym), eigvec(1,:,iat,imode)) ** 2 + &
             matmul(symcart(:,:,isym), eigvec(2,:,iat,imode)) ** 2
         end do
       end do
       !syme2_xyza = syme2_xyza / crystal%nsym

       ! Accumulate PJDOS
       do iat=1,natom
         do idir=1,3
           phdos%pjdos(:,idir,iat) = phdos%pjdos(:,idir,iat) + syme2_xyza(idir,iat) * gvals_wtq
         end do
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
         do isym=1, crystal%nsym
           temp_33 = matmul( (symcart(:,:,isym)), matmul(msqd_atom_tmp, transpose(symcart(:,:,isym))) )
           jat = crystal%indsym(4,isym,iat)
           do idir=1,3
             do jdir=1,3
               phdos%msqd_dos_atom(:,idir,jdir,jat) = phdos%msqd_dos_atom(:,idir,jdir,jat) + &
                 temp_33(idir, jdir) * gvals_wtq
             end do
           end do
         end do

       end do ! iat
     end do ! imode

   case (2)
     ! Tetrahedra
     ! Save phonon frequencies and eigenvectors.
     ! Sum is done after the loops over the two meshes.
     full_phfrq(:,iq_ibz) = phfrq(:)
     full_eigvec(:,:,:,:,iq_ibz) = eigvec

   case default
     MSG_ERROR(sjoin("Wrong value for prtdos:", itoa(prtdos)))
   end select
 end do ! iq_ibz

 ABI_FREE(qbz)
 ABI_FREE(gvals_wtq)
 ABI_FREE(xvals)

 call xmpi_sum_master(nsmallq, master, comm, ierr)
 call xmpi_sum_master(speedofsound, master, comm, ierr)

 call cwtime_report(" phdos", cpu, wall, gflops)

 if (my_rank == master .and. nsmallq > tol10) then
   ! Write info about speed of sound
   speedofsound = speedofsound / nsmallq
   write (msg,'(a,E20.10,3a,F16.4,2a)') &
       ' Average speed of sound partial sums: ', third*sum(speedofsound), ' (at units)',ch10, &
       '-                                   = ', third*sum(speedofsound) * Bohr_Ang * 1.d-13 / Time_Sec, ' [km/s]',ch10
   call wrtout (ab_out,msg,"COLL")
   call wrtout (std_out,msg,"COLL")

   ! Debye frequency = vs * (6 pi^2 natom / ucvol)**1/3
   debyefreq = third*sum(speedofsound) * (six*pi**2/crystal%ucvol)**(1./3.)
   write (msg,'(a,E20.10,3a,E20.10,a)') &
      ' Debye frequency from partial sums: ', debyefreq, ' (Ha)',ch10, &
      '-                                 = ', debyefreq*Ha_THz, ' (THz)'
   call wrtout (ab_out,msg,"COLL")
   call wrtout (std_out,msg,"COLL")

   ! Debye temperature = hbar * Debye frequency / kb
   write (msg,'(a,E20.10,2a)') '-Debye temperature from partial sums: ', debyefreq*Ha_K, ' (K)', ch10
   call wrtout (ab_out,msg,"COLL")
   call wrtout (std_out,msg,"COLL")
 end if

 if (prtdos == 2) then
   ! Finalize integration with tetrahedra
   ! All the data are contained in full_phfrq and full_eigvec.
   call xmpi_sum(full_phfrq, comm, ierr)
   call xmpi_sum(full_eigvec, comm, ierr)

   ABI_MALLOC(tmp_phfrq, (phdos%nqibz))

   do
     ABI_MALLOC(wdt, (phdos%nomega, 2))
     ABI_MALLOC(energies, (phdos%nomega))
     energies = linspace(phdos%omega_min,phdos%omega_max,phdos%nomega)

     do iq_ibz=1,phdos%nqibz
       if (mod(iq_ibz, nprocs) /= my_rank) cycle ! mpi-parallelism

       ! Compute the weights for this q-point using tetrahedron
       do imode=1,3*natom
         tmp_phfrq(:) = full_phfrq(imode,:)
         call htetraq%get_onewk_wvals(iq_ibz,bcorr0,phdos%nomega,energies,max_occ1,phdos%nqibz,tmp_phfrq,wdt)
         wdt = wdt * wtq_ibz(iq_ibz)

         ! Accumulate DOS/IDOS
         phdos%phdos(:)     = phdos%phdos(:)     + wdt(:, 1)
         phdos%phdos_int(:) = phdos%phdos_int(:) + wdt(:, 2)

         ! Rotate e(q) to get e(Sq) to account for other q-points in BZ. See notes in gaussian branch
         syme2_xyza = zero
         do iat=1,natom
           do isym=1, crystal%nsym
             jat = crystal%indsym(4,isym,iat)
             syme2_xyza(:,jat) = syme2_xyza(:,jat) + &
               matmul(symcart(:,:,isym), full_eigvec(1,:,iat,imode,iq_ibz)) ** 2 + &
               matmul(symcart(:,:,isym), full_eigvec(2,:,iat,imode,iq_ibz)) ** 2
           end do
         end do
         !syme2_xyza = syme2_xyza / crystal%nsym

         do iat=1,natom
           do idir=1,3
             phdos%pjdos(:,idir,iat) = phdos%pjdos(:,idir,iat) + syme2_xyza(idir,iat) * wdt(:,1)
             phdos%pjdos_int(:,idir,iat) = phdos%pjdos_int(:,idir,iat) + syme2_xyza(idir,iat) * wdt(:,2)
           end do
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
           ! TODO: need to check the direction of the symcart vs transpose or inverse, given that jat is the pre-image of iat...
           do isym=1, crystal%nsym
             temp_33 = matmul( (symcart(:,:,isym)), matmul(msqd_atom_tmp, transpose(symcart(:,:,isym))) )
             jat = crystal%indsym(4,isym,iat)
             do idir=1,3
               do jdir=1,3
                 phdos%msqd_dos_atom(:,idir,jdir,jat) = phdos%msqd_dos_atom(:,idir,jdir,jat) + &
                   temp_33(idir, jdir) * wdt(:,1)
               end do
             end do
           end do
         end do ! iat

       end do ! imode
     end do ! iq_ibz

     if (refine_dosdeltae) then
       ! HM: Check if the integration of the DOS is correct, otherwise half dos%deltae and re-run
       call ctrap(phdos%nomega, phdos%phdos, phdos%omega_step, phdos_int)
       if (abs(phdos_int-crystal%natom*3)>tol2) then
         write(msg,'(a,f6.2,a,i4,2a,e10.3,a,e10.3)') "The value of the integral is", phdos_int, &
                      " but it should be", crystal%natom*3, ch10,&
                      "I will decrease dosdeltae from", dosdeltae, " to", dosdeltae/two
         MSG_WARNING(msg)
         dosdeltae = dosdeltae / two
         call phdos_free(phdos)
         call phdos_init(phdos, crystal, ifc, dosdeltae, dossmear, wminmax, prtdos)
         nomega = phdos%nomega
         ABI_FREE(wdt)
         ABI_FREE(energies)
         cycle
       endif
     endif
     exit
   end do
   call cwtime_report(" accumulate", cpu, wall, gflops)
   ABI_FREE(energies)
   ABI_FREE(wdt)

   ! Make eigvec into phonon displacements.
   do iat = 1, natom
     full_eigvec(:,:,iat,:,:) = full_eigvec(:,:,iat,:,:) / sqrt(phdos%atom_mass(iat))
   end do

   if (my_rank == master) then
#ifdef HAVE_NETCDF
     ! TODO: make it optional?
     NCF_CHECK_MSG(nctk_open_create(ncid, strcat(prefix, "_PHIBZ.nc"), xmpi_comm_self), "Creating PHIBZ")
     NCF_CHECK(crystal%ncwrite(ncid))
     call phonons_ncwrite(ncid, natom, phdos%nqibz, qibz, wtq_ibz, full_phfrq, full_eigvec)
     NCF_CHECK(nf90_close(ncid))
#endif
   end if

   ! immediately free this - it contains displ and not eigvec at this stage
   ABI_FREE(full_eigvec)
   ABI_FREE(full_phfrq)
   ABI_FREE(tmp_phfrq)
   call htetraq%free()
 else
#ifdef HAVE_NETCDF
   MSG_WARNING('The netcdf PHIBZ file is only output for tetrahedron integration and DOS calculations')
#endif
 end if ! tetrahedra

 ! Test if the initial mesh was large enough
 call xmpi_sum(count_wminmax, comm, ierr)
 call xmpi_min(wminmax(1), rtmp, comm, ierr); wminmax(1) = rtmp
 call xmpi_max(wminmax(2), rtmp, comm, ierr); wminmax(2) = rtmp

 call xmpi_sum(phdos%phdos, comm, ierr)
 call xmpi_sum(phdos%msqd_dos_atom, comm, ierr)
 call xmpi_sum(phdos%pjdos, comm, ierr)
 if (prtdos == 2) then
   ! Reduce integrals if tetra, gauss will compute idos with simpson.
   call xmpi_sum(phdos%phdos_int, comm, ierr)
   call xmpi_sum(phdos%pjdos_int, comm, ierr)
 end if

 ! normalize by nsym: symmetrization is used in all prtdos cases
 phdos%msqd_dos_atom = phdos%msqd_dos_atom / crystal%nsym
 phdos%pjdos = phdos%pjdos / crystal%nsym
 if (prtdos == 2) phdos%pjdos_int = phdos%pjdos_int / crystal%nsym
 ! normalize by mass and factor of 2, now added in the printout to agree with harmonic_thermo
 ! do iat=1, natom
 !   phdos%msqd_dos_atom(:,:,:,iat) = phdos%msqd_dos_atom(:,:,:,iat) * invmass(iat) * half
 ! end do ! iat

 ! ===============================
 ! === Compute Integrated PDOS ===
 ! ===============================
 ABI_CALLOC(phdos%pjdos_rc_type, (nomega, 3, crystal%ntypat))
 ABI_CALLOC(phdos%pjdos_type, (nomega, crystal%ntypat))
 ABI_CALLOC(phdos%pjdos_type_int, (nomega, crystal%ntypat))

 do iat=1,natom
   itype = crystal%typat(iat)
   do io=1,phdos%nomega
     phdos%pjdos_rc_type(io,:,itype) = phdos%pjdos_rc_type(io,:,itype) + phdos%pjdos(io,:,iat)
     phdos%pjdos_type(io,itype) = phdos%pjdos_type(io,itype) + sum(phdos%pjdos(io,:,iat))
   end do
   if (prtdos == 2) then
     do io=1,phdos%nomega
       phdos%pjdos_type_int(io,itype) = phdos%pjdos_type_int(io,itype) + sum(phdos%pjdos_int(io,:,iat))
     end do
   end if
 end do

 ! Evaluate IDOS using simple simpson integration
 ! In principle one could use derf.F90, just to be consistent ...
 if (prtdos == 1) then
   call simpson_int(phdos%nomega,phdos%omega_step, phdos%phdos, phdos%phdos_int)
   do iat=1,natom
     do idir=1,3
       call simpson_int(phdos%nomega, phdos%omega_step, phdos%pjdos(:,idir,iat), phdos%pjdos_int(:,idir,iat))
     end do
   end do
   do itype=1,crystal%ntypat
     call simpson_int(phdos%nomega, phdos%omega_step, phdos%pjdos_type(:,itype), phdos%pjdos_type_int(:,itype))
   end do
 end if

 ABI_FREE(qibz)
 ABI_FREE(wtq_ibz)

 call cwtime(cpu, wall, gflops, "stop")
 write(msg,'(2(a,f8.2))')" mkphdos completed. cpu:", cpu, ", wall:", wall
 call wrtout(std_out, msg, do_flush=.True.)

 DBG_EXIT("COLL")

end subroutine mkphdos
!!***

!----------------------------------------------------------------------

!!****f* m_phonons/zacharias_supercell_make
!! NAME
!! zacharias_supercell_make
!!
!! FUNCTION
!!  Construct an optimally thermalized supercell following Zacharias and Giustino
!!  PRB 94 075125 (2016) [[cite:Zacharias2016]]
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
!!
!! SOURCE

subroutine zacharias_supercell_make(Crystal, Ifc, ntemper, rlatt, tempermin, temperinc, thm_scells)

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
 real(dp), allocatable :: qbz(:,:), qibz(:,:), wtq_ibz(:)
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
&   nqibz, qibz, wtq_ibz, nqbz, qbz) ! new_kptrlatt, new_shiftk)  ! Optional
 ABI_FREE(qshft)

 ! allocate arrays with all of the q, omega, and displacement vectors
 ABI_MALLOC_OR_DIE(phfrq_allq, (3*Crystal%natom*nqibz), ierr)
 ABI_MALLOC_OR_DIE(phdispl_allq, (2, 3, Crystal%natom, 3*Crystal%natom, nqibz), ierr)

 ABI_MALLOC_OR_DIE(phfrq, (3*Crystal%natom), ierr)
 ABI_MALLOC_OR_DIE(phdispl, (2, 3, Crystal%natom, 3*Crystal%natom), ierr)
 ABI_MALLOC_OR_DIE(pheigvec, (2, 3, Crystal%natom, 3*Crystal%natom), ierr)

 ! loop over q to get all frequencies and displacement vectors
 ABI_ALLOCATE(modeindex, (nqibz*3*Crystal%natom))
 imode = 0
 do iq = 1, nqibz
   ! Fourier interpolation.
   call ifc%fourq(Crystal, qibz(:,iq), phfrq, phdispl, out_eigvec=pheigvec)
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

 ABI_MALLOC(phdispl1, (2, 3, Crystal%natom))
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
 ABI_FREE(wtq_ibz)

end subroutine zacharias_supercell_make
!!***

!----------------------------------------------------------------------

!!****f* m_phonons/thermal_supercell_make
!! NAME
!! thermal_supercell_make
!!
!! FUNCTION
!!  Construct an random thermalized supercell configuration, as in TDEP
!!  main function is for training set generation in multibinit
!!
!! INPUTS
!!   Crystal = crystal object with rprim etc...
!!   Ifc = interatomic force constants object from anaddb
!!   option = option to deal with negative frequency -> Bose factor explodes (eg acoustic at Gamma)
!!      several philosophies to be implemented for the unstable modes:
!!      option == 1 =>  ignore
!!      option == 2 =>  populate them according to a default amplitude
!!      option == 3 =>  populate according to their modulus squared
!!      option == 4 =>  USER defined value(s), require namplitude and amplitude
!!   nconfig = numer of requested configurations
!!   rlatt = matrix of conversion for supercell (3 0 0   0 3 0   0 0 3 for example)
!!   temperature_K =  temperature in Kelvin
!!   nqpt = number of q-point
!!   namplitude = number of amplitude provided by the user
!!   amplitudes(namplitude) = list of the amplitudes of the unstable phonons
!!                            amplitudes(1:3,iamplitude) = qpt
!!                            amplitudes(4,iamplitude)   = mode
!!                            amplitudes(5,iamplitude)   = amplitude
!!
!! OUTPUT
!!   thm_scells = array of configurations with thermalized supercells
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine thermal_supercell_make(amplitudes,Crystal, Ifc,namplitude, nconfig,option,&
&                                 rlatt, temperature_K, thm_scells)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: option,nconfig
 integer, intent(in) :: rlatt(3,3)
 real(dp), intent(in) :: temperature_K
 type(crystal_t),intent(in) :: Crystal
 type(ifc_type),intent(in) :: Ifc
 type(supercell_type), intent(out) :: thm_scells(nconfig)
 integer,intent(in) :: namplitude
!Local variables-------------------------------
!scalars
 integer :: iq, nqibz, nqbz, qptopt1, iampl ,imode, ierr, iconfig
 real(dp) :: temperature, sigma, freeze_displ
 real(dp) :: rand !, rand1, rand2
 real(dp),intent(in):: amplitudes(5,namplitude)
 !arrays
 real(dp), allocatable :: qshft(:,:) ! dummy with 2 dimensions for call to kpts_ibz_from_kptrlatt
 real(dp), allocatable :: qbz(:,:), qibz(:,:), wtqibz(:)
 real(dp), allocatable :: phfrq_allq(:,:), phdispl_allq(:,:,:,:,:)
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

 temperature = temperature_K /  Ha_K

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

 ! allocate arrays wzith all of the q, omega, and displacement vectors
 ABI_MALLOC_OR_DIE(phfrq_allq, (3*Crystal%natom, nqibz), ierr)
 ABI_MALLOC_OR_DIE(phdispl_allq, (2, 3, Crystal%natom, 3*Crystal%natom, nqibz), ierr)

 ABI_MALLOC_OR_DIE(phfrq, (3*Crystal%natom), ierr)
 ABI_MALLOC_OR_DIE(phdispl, (2, 3, Crystal%natom, 3*Crystal%natom), ierr)
 ABI_MALLOC_OR_DIE(pheigvec, (2, 3, Crystal%natom, 3*Crystal%natom), ierr)

 ! loop over q to get all frequencies and displacement vectors
 imode = 0
 do iq = 1, nqibz
   ! Fourier interpolation.
   call ifc%fourq(Crystal, qibz(:,iq), phfrq, phdispl, out_eigvec=pheigvec)
   phfrq_allq(1:3*Crystal%natom, iq) = phfrq
   phdispl_allq(1:2, 1:3, 1:Crystal%natom, 1:3*Crystal%natom, iq) = phdispl
 end do
 ABI_FREE(phfrq)
 ABI_FREE(pheigvec)
 ABI_FREE(phdispl)

 ! only diagonal supercell case for the moment
 do iconfig = 1, nconfig
   call init_supercell(Crystal%natom, rlatt, Crystal%rprimd, Crystal%typat, Crystal%xcart, Crystal%znucl, thm_scells(iconfig))
 end do

 ! precalculate phase factors???

 ABI_MALLOC_OR_DIE(phdispl1, (2, 3, Crystal%natom), ierr)

 ! for all modes at all q in whole list, sorted
 do iq = 1, nqibz
   do imode = 1, 3*Crystal%natom

     ! skip modes with too low or negative frequency -> Bose factor explodes (eg acoustic at Gamma)
     ! TODO: check the convergence wrt the tolerance
     ! several philosophies to be implemented for the unstable modes:
     ! 1) ignore
     ! 2) populate them according to a default amplitude
     ! 3) populate according to their modulus squared
     if (abs(phfrq_allq(imode, iq))<tol6) cycle

     phdispl1 = phdispl_allq(:,:,:,imode,iq)

     ! loop over configurations
     do iconfig = 1, nconfig

       ! trick supercell object into using present q point
       thm_scells(iconfig)%qphon(:) = qibz(:,iq)

       ! find thermal displacement amplitude eq 4 of Zacharias
       !   combined with l_nu,q expression in paragraph before
       if (phfrq_allq(imode, iq) > tol6) then
         sigma = sqrt( (bose_einstein(phfrq_allq(imode,iq), temperature) + half)/phfrq_allq(imode,iq))
       else
         !Treat negative frequencies
         select case (option)
         case(1)
           !Do not populate
           sigma = 0._dp
         case(2)
           !Default amplitude for all the frequencies
           sigma = 100._dp
         case(3)
           !Absolute value of the frequencies
           sigma=sqrt((bose_einstein(abs(phfrq_allq(imode,iq)),temperature)+half)/&
&                abs(phfrq_allq(imode,iq)))
         case(4)
           sigma = 0._dp
           !Search if the amplitude of this unstable phonon is in the input argument amplitudes
           do iampl=1,namplitude
             if(abs(thm_scells(iconfig)%qphon(1) - amplitudes(1,iampl)) < tol8.and.&
&               abs(thm_scells(iconfig)%qphon(2) - amplitudes(2,iampl)) < tol8.and.&
&               abs(thm_scells(iconfig)%qphon(3) - amplitudes(3,iampl)) < tol8.and.&
&               abs(imode - amplitudes(4,iampl)) < tol8) then
               sigma = amplitudes(5,iampl)
             end if
           end do
           !If not, the amplitude is zero
           if(abs(sigma) < tol8)then
             write (msg, '(a,I0,a,3es12.5,2a,I0)') ' The amplitude of the unstable mode ',&
&                int(imode),' of the qpt ',thm_scells(iconfig)%qphon(:), ch10,&
&                'is set to zero for the configuration ',iconfig
             MSG_WARNING(msg)
           end if
         end select
       end if

       ! add displacement for this mode to supercell positions eq 5 of Zacharias
       call RANDOM_NUMBER(rand)
       rand = two * rand - one

       ! from TDEP documentation for gaussian distribution of displacements
       !call RANDOM_NUMBER(rand1)
       !call RANDOM_NUMBER(rand2)
       ! rand = sqrt(-two*rand1) * sin(twopi*rand2)

       ! if (rand > half) then
       !   rand = one
       ! else
       !   rand = -one
       ! end if

       freeze_displ =  rand * sigma

       call freeze_displ_supercell (phdispl1(:,:,:), freeze_displ, thm_scells(iconfig))
     end do !iconfig
   end do !imode
 end do !iq

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
!!
!! SOURCE

subroutine thermal_supercell_free(nscells, thm_scells)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: nscells
 type(supercell_type), allocatable, intent(inout) :: thm_scells(:)

! local
 integer :: icell

 if(allocated(thm_scells)) then
   do icell = 1, nscells
     call destroy_supercell(thm_scells(icell))
   end do
 end if
end subroutine thermal_supercell_free
!!***

!----------------------------------------------------------------------

!!****f* m_phonons/zacharias_supercell_print
!! NAME
!! zacharias_supercell_print
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
!!
!! SOURCE

subroutine zacharias_supercell_print(fname, ntemper, tempermin, temperinc, thm_scells)

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
   write (title1, '(3a)') "#  Zacharias thermalized supercell at temperature T= ", trim(temper_str), " Kelvin"
   title2 = "#  generated with alternating thermal displacements of all phonons"
   call prt_supercell (filename, thm_scells(itemp), title1, title2)
 end do

end subroutine zacharias_supercell_print
!!***

!!****f* m_phonons/thermal_supercell_print
!! NAME
!! thermal_supercell_print
!!
!! FUNCTION
!!  print files with thermalized array of random supercell configurations
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
!!
!! SOURCE

subroutine thermal_supercell_print(fname, nconfig, temperature_K, thm_scells)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: nconfig
 type(supercell_type), intent(in) :: thm_scells(nconfig)
 character(len=fnlen), intent(in) :: fname
 real(dp), intent(in) :: temperature_K

! local
 integer :: iconfig,itemp
 character(len=80) :: title1, title2
 character(len=fnlen) :: filename
 character(len=10) :: config_str

 do iconfig = 1, nconfig
   write (config_str,'(I8)') iconfig
   write (filename, '(3a)') trim(fname), "_cf_", trim(adjustl(config_str))
   write (title1, '(a,I6,a)') "#  thermalized supercell at temperature T= ", temperature_K, " Kelvin"
   title2 = "#  generated with random thermal displacements of all phonons"
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
!!
!! SOURCE

subroutine phdos_ncwrite(phdos,ncid)

!Arguments ------------------------------------
!scalars
 class(phonon_dos_type),intent(in) :: phdos
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
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "prtdos"), phdos%prtdos))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'dossmear'), phdos%dossmear*Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'wmesh'), phdos%omega*Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'phdos'), phdos%phdos/Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'pjdos'), phdos%pjdos/Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'pjdos_type'), phdos%pjdos_type/Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'pjdos_rc_type'), phdos%pjdos_rc_type/Ha_eV))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'msqd_dos_atom'), phdos%msqd_dos_atom/Ha_eV))

#else
 MSG_ERROR("netcdf support not enabled")
 ABI_UNUSED((/ncid, phdos%nomega/))
#endif

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
!!
!! SOURCE

subroutine mkphbs(Ifc,Crystal,inp,ddb,asrq0,prefix,comm)

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
     !call ifc%fourq(cryst, save_qpoints(:,iphl1), phfrq, displ, out_eigvec=eigvec)

     ! Get d2cart using the interatomic forces and the
     ! long-range coulomb interaction through Ewald summation
     call gtdyn9(ddb%acell,Ifc%atmfrc,Ifc%dielt,Ifc%dipdip,Ifc%dyewq0,d2cart,Crystal%gmet,ddb%gprim,ddb%mpert,natom, &
      Ifc%nrpt,qphnrm(1),qphon,Crystal%rmet,ddb%rprim,Ifc%rpt,Ifc%trans,Crystal%ucvol,Ifc%wghatm,Crystal%xred,ifc%zeff,ifc%qdrp_cart,&
      xmpi_comm_self)

   else if (ifcflag == 0) then

     !call ddb_diagoq(ddb, crystal, save_qpoints(:,iphl1), asrq0, ifc%symdynmat, rftyp, phfrq, displ, &
     !                out_eigvec=eigvec)

     ! Look for the information in the DDB (no interpolation here!)
     rfphon(1:2)=1; rfelfd(1:2)=0; rfstrs(1:2)=0
     qphon_padded = zero; qphon_padded(:,1) = qphon

     call ddb%get_block(iblok,qphon_padded,qphnrm,rfphon,rfelfd,rfstrs,rftyp)

     ! Copy the dynamical matrix in d2cart
     d2cart(:,1:ddb%msize)=ddb%val(:,:,iblok)

     ! Eventually impose the acoustic sum rule based on previously calculated d2asr
     call asrq0%apply(natom, ddb%mpert, ddb%msize, crystal%xcart, d2cart)
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
   NCF_CHECK(crystal%ncwrite(ncid))
   call phonons_ncwrite(ncid,natom,nfineqpath,save_qpoints,weights,save_phfrq,save_phdispl_cart)

   ! Now treat the second list of vectors (only at the Gamma point, but can include non-analyticities)
   if (inp%nph2l /= 0 .and. inp%ifcflag == 1) then
     call ifc%calcnwrite_nana_terms(crystal, inp%nph2l, inp%qph2l, inp%qnrml2, ncid)
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
!!
!! SOURCE

subroutine phdos_calc_vsound(eigvec,gmet,natom,phfrq,qphon,speedofsound)

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
   imode_acoustic = min(imode_acoustic + 1, 3)

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
!!
!! SOURCE

subroutine phdos_print_vsound(iunit,ucvol,speedofsound)

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
!!  Only master node should call this routine.
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
!!
!! SOURCE

subroutine phdos_print_msqd(PHdos, fname, ntemper, tempermin, temperinc)

!Arguments -------------------------------
!scalars
 integer, intent(in) :: ntemper
 class(phonon_dos_type),intent(in) :: PHdos
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

 ! write the header
 write (iunit, '(a)') '# mean square displacement for each atom as a function of T (bohr^2)'
 write (junit, '(a)') "# mean square velocity for each atom as a function of T (bohr^2/atomic time unit^2)"

 write (msg, '(a,F18.10,a,F18.10,a)') '#  T in Kelvin, from ', tempermin, ' to ', tempermin+(ntemper-1)*temperinc
 write (iunit, '(a)') trim(msg)
 write (junit, '(a)') trim(msg)

 write (msg, '(2a)') '#    T             |u^2|                u_xx                u_yy                u_zz',&
&                                            '                u_yz                u_xz                u_xy in bohr^2'
 write (iunit, '(a)') trim(msg)
 write (msg, '(3a)') '#    T             |v^2|                v_xx                v_yy                v_zz',&
&                                            '                v_yz                v_xz                v_xy',&
&                                            ' in bohr^2/atomic time unit^2'
 write (junit, '(a)') trim(msg)

 ABI_ALLOCATE (tmp_msqd, (PHdos%nomega,9))
 ABI_ALLOCATE (tmp_msqv, (PHdos%nomega,9))
 ABI_ALLOCATE (integ_msqd, (9,ntemper))
 ABI_ALLOCATE (integ_msqv, (9,ntemper))
 ABI_ALLOCATE (bose_msqd, (PHdos%nomega, ntemper))
 ABI_ALLOCATE (bose_msqv, (PHdos%nomega, ntemper))

 do io=1, PHdos%nomega
   if ( PHdos%omega(io) >= 2._dp * 4.56d-6 ) exit ! 2 cm-1 TODO: make this an input parameter
 end do
 iomin = io

 ! calculate bose only once for each atom (instead of for each atom)
 bose_msqd = zero
 bose_msqv = zero
 do itemp = 1, ntemper
   temper = tempermin + (itemp-1) * temperinc
   if (temper < 1.e-3) cycle ! millikelvin at least to avoid exploding Bose factor(TM)
   do io = iomin, PHdos%nomega
     ! NB: factors follow convention in phonopy documentation
     ! the 1/sqrt(omega) factor in phonopy is contained in the displacement vector definition
     ! bose() is dimensionless
     !bose_msqd(io, itemp) =  (half + one  / ( exp(PHdos%omega(io)/(kb_HaK*temper)) - one )) / PHdos%omega(io)
     !bose_msqv(io, itemp) =  (half + one  / ( exp(PHdos%omega(io)/(kb_HaK*temper)) - one )) * PHdos%omega(io)
     bose_msqd(io, itemp) =  (half + bose_einstein(PHdos%omega(io),kb_HaK*temper)) / PHdos%omega(io)
     bose_msqv(io, itemp) =  (half + bose_einstein(PHdos%omega(io),kb_HaK*temper)) * PHdos%omega(io)
   end do
 end do

 do iatom=1,PHdos%natom
   write (msg, '(a,I8)') '# atom number ', iatom
   write (iunit, '(a)') trim(msg)
   write (junit, '(a)') trim(msg)

   ! for each T and each atom, integrate msqd matrix with Bose Einstein factor and output
   integ_msqd = zero
   tmp_msqd = reshape(PHdos%msqd_dos_atom(:,:,:,iatom), (/PHdos%nomega, 9/))

   ! perform all integrations as matrix multiplication: integ_msqd (idir, itemp) = [tmp_msqd(io,idir)]^T  * bose_msqd(io,itemp)
   call DGEMM('T','N', 9, ntemper, PHdos%nomega, one, tmp_msqd,PHdos%nomega, bose_msqd, PHdos%nomega, zero, integ_msqd, 9)
   ! NB: this presumes an equidistant omega grid
   integ_msqd = integ_msqd * (PHdos%omega(2)-PHdos%omega(1)) / PHdos%atom_mass(iatom)

   integ_msqv = zero
   tmp_msqv = reshape(PHdos%msqd_dos_atom(:,:,:,iatom), (/PHdos%nomega, 9/))

   ! perform all integrations as matrix multiplication: integ_msqv (idir, itemp) = [tmp_msqv(io,idir)]^T  * bose_msqv(io,itemp)
   call DGEMM('T','N', 9, ntemper, PHdos%nomega, one, tmp_msqv,PHdos%nomega, bose_msqv, PHdos%nomega, zero, integ_msqv, 9)
   ! NB: this presumes an equidistant omega grid
   integ_msqv = integ_msqv * (PHdos%omega(2)-PHdos%omega(1)) / PHdos%atom_mass(iatom)

   ! print out stuff
   do itemp = 1, ntemper
     temper = tempermin + (itemp-1) * temperinc
     write (msg, '(F10.2,4x,E22.10,2x,6E22.10)') temper, third*(integ_msqd(1,itemp)+integ_msqd(5,itemp)+integ_msqd(9,itemp)), &
&                      integ_msqd(1,itemp),integ_msqd(5,itemp),integ_msqd(9,itemp), &
&                      integ_msqd(6,itemp),integ_msqd(3,itemp),integ_msqd(2,itemp)
     write (iunit, '(a)') trim(msg)
     write (msg, '(F10.2,4x,E22.10,2x,6E22.10)') temper, third*(integ_msqv(1,itemp)+integ_msqv(5,itemp)+integ_msqv(9,itemp)), &
&                      integ_msqv(1,itemp),integ_msqv(5,itemp),integ_msqv(9,itemp), &
&                      integ_msqv(6,itemp),integ_msqv(3,itemp),integ_msqv(2,itemp)
     write (junit, '(a)') trim(msg)
   end do ! itemp

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
!!
!! SOURCE

subroutine phonons_ncwrite(ncid,natom,nqpts,qpoints,weights,phfreq,phdispl_cart)

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
!!
!! SOURCE

 subroutine phonons_write_phfrq(path,natom,nqpts,qpoints,weights,phfreq,phdispl_cart)

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
!!
!! SOURCE

subroutine phonons_write_xmgrace(filename, natom, nqpts, qpts, phfreqs, qptbounds)

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
!!
!! SOURCE

subroutine phonons_write_gnuplot(prefix, natom, nqpts, qpts, phfreqs, qptbounds)

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
!!
!! SOURCE

subroutine ifc_mkphbs(ifc, cryst, dtset, prefix, comm)

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
 integer :: iqpt,nqpts,natom,ncid,nprocs,my_rank,ierr,nph2l
 type(kpath_t) :: qpath
!arrays
 real(dp),allocatable :: qph2l(:,:), qnrml2(:)
 real(dp),allocatable :: eigvec(:,:,:,:,:),phfrqs(:,:),phdispl_cart(:,:,:,:),weights(:)

! *********************************************************************

 if (dtset%prtphbands == 0) return

 if (dtset%ph_nqpath <= 0 .or. dtset%ph_ndivsm <= 0) then
   MSG_COMMENT("ph_nqpath <= 0 or ph_ndivsm <= 0. Phonon bands won't be produced. Returning")
   return
 end if
 call wrtout(std_out, " Writing phonon bands, use prtphbands 0 to disable this part")

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
   call ifc%fourq(cryst, qpath%points(:,iqpt), phfrqs(:,iqpt), phdispl_cart(:,:,:,iqpt), out_eigvec=eigvec)
 end do

 call xmpi_sum_master(phfrqs, master, comm, ierr)
 call xmpi_sum_master(phdispl_cart, master, comm, ierr)

 if (my_rank == master) then
   ABI_MALLOC(weights, (nqpts))
   weights = one

   ! Compute directions for non-analytical behaviour.
   ! TODO: The same approach should be used in anaddb.
   ABI_MALLOC(qph2l, (3, 2*dtset%ph_nqpath))
   ABI_MALLOC(qnrml2, (2*dtset%ph_nqpath))
   nph2l = 0
   if (any(ifc%zeff /= zero)) then
     do iqpt=1,dtset%ph_nqpath
       if (sum(dtset%ph_qpath(:, iqpt)**2) < tol14) then
         nph2l = nph2l + 1
         if (iqpt == 1) then
           qph2l(:, nph2l) = dtset%ph_qpath(:, 2) - dtset%ph_qpath(:, 1)
         else if (iqpt == dtset%ph_nqpath) then
           qph2l(:, nph2l) = dtset%ph_qpath(:, dtset%ph_nqpath - 1) - dtset%ph_qpath(:, dtset%ph_nqpath)
         else
           qph2l(:, nph2l) = dtset%ph_qpath(:, iqpt - 1) - dtset%ph_qpath(:, iqpt)
           nph2l = nph2l + 1
           qph2l(:, nph2l) = dtset%ph_qpath(:, iqpt + 1) - dtset%ph_qpath(:, iqpt)
         end if
       end if
     end do
     ! Convert to Cartesian coordinates.
     do iqpt=1,nph2l
       qph2l(:, iqpt) = matmul(cryst%gprimd, qph2l(:, iqpt))
     end do
     qnrml2 = zero
   end if

#ifdef HAVE_NETCDF
   ! TODO: A similar piece of code is used in anaddb (mkpbs + ifc_calcnwrite_nana_terms).
   ! Should centralize everything in a single routine
   NCF_CHECK_MSG(nctk_open_create(ncid, strcat(prefix, "_PHBST.nc"), xmpi_comm_self), "Creating PHBST")
   NCF_CHECK(cryst%ncwrite(ncid))
   call phonons_ncwrite(ncid, natom, nqpts, qpath%points, weights, phfrqs, phdispl_cart)
   NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('atomic_mass_units', "dp", "number_of_atom_species")], defmode=.True.))
   NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'atomic_mass_units'), ifc%amu))
   if (nph2l /= 0) call ifc%calcnwrite_nana_terms(cryst, nph2l, qph2l, qnrml2, ncid=ncid)
   NCF_CHECK(nf90_close(ncid))
#endif

   ABI_FREE(qph2l)
   ABI_FREE(qnrml2)

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

 call qpath%free()

end subroutine ifc_mkphbs
!!***

!!****f* m_phonons/dfpt_symph
!! NAME
!! dfpt_symph
!!
!! FUNCTION
!! Determine the symmetry character of the different phonon modes.
!!
!! INPUTS
!! acell(3)=length scales of primitive translations (bohr)
!! eigvec(2*3*natom*3*natom)=eigenvectors of the dynamical matrix
!! indsym(4,nsym,natom)=indirect indexing array : for each
!!   isym,iatom, fourth element is label of atom into which iatom is sent by
!!   INVERSE of symmetry operation isym; first three elements are the primitive
!!   translations which must be subtracted after the transformation to get back
!!   to the original unit cell.
!! iout=unit number to which output is written
!! natom=number of atoms in unit cell
!! nsym=number of space group symmetries
!! phfrq(3*natom)=phonon frequencies (Hartree units)
!! rprim(3,3)=dimensionless primitive translations in real space
!! symrel(3,3,nsym)=matrices of the group symmetries (real space)
!!
!! OUTPUT
!!
!! PARENTS
!!      anaddb,m_phonons
!!
!! CHILDREN
!!      matr3inv,mkrdim,wrtout
!!
!! SOURCE

subroutine dfpt_symph(iout,acell,eigvec,indsym,natom,nsym,phfrq,rprim,symrel)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,natom,nsym
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symrel(3,3,nsym)
 real(dp),intent(in) :: acell(3),eigvec(2*3*natom*3*natom),phfrq(3*natom)
 real(dp),intent(in) :: rprim(3,3)

!Local variables -------------------------
!scalars
 integer :: iad1,iad2,iad3,iatom,idir,ii1,ii2,ii3,imode,isym,itol,jad,jatom,jj
 integer :: jmode,kk,ntol
 character(len=500) :: message
!arrays
 integer,allocatable :: degeneracy(:),integer_characters(:),symind(:,:)
 real(dp) :: gprimd(3,3),rprimd(3,3)
 real(dp),allocatable :: eigvtr(:),redvec(:),redvtr(:),symph(:,:)

!******************************************************************

!Compute dimensional primitive translations rprimd and its inverse gprimd
 call mkrdim(acell,rprim,rprimd)
 call matr3inv(rprimd,gprimd)

!Build the symmetry index (inverse of indsym(4,:,:))
 ABI_ALLOCATE(symind,(nsym,natom))
 do isym=1,nsym
   do iatom=1,natom
     symind(isym,indsym(4,isym,iatom))=iatom
   end do
 end do

 ABI_ALLOCATE(symph,(nsym,3*natom))
 ABI_ALLOCATE(redvec,(2*3*natom))
 ABI_ALLOCATE(redvtr,(2*3*natom))
 ABI_ALLOCATE(eigvtr,(2*3*natom))

!Loop over the vibration modes
 do imode=1,3*natom

!  Compute eigvec for this mode in reduced coordinates redvec
   do iatom=1,natom
     iad1=3*(iatom-1)+1
     ii1=2*3*natom*(imode-1)+2*(iad1-1)+1
     iad2=3*(iatom-1)+2
     ii2=2*3*natom*(imode-1)+2*(iad2-1)+1
     iad3=3*(iatom-1)+3
     ii3=2*3*natom*(imode-1)+2*(iad3-1)+1
     do idir=1,3
       jad=3*(iatom-1)+idir
       jj=2*(jad-1)+1
       redvec(jj)=gprimd(1,idir)*eigvec(ii1)+&
&       gprimd(2,idir)*eigvec(ii2)+&
&       gprimd(3,idir)*eigvec(ii3)
       redvec(jj+1)=gprimd(1,idir)*eigvec(ii1+1)+&
&       gprimd(2,idir)*eigvec(ii2+1)+&
&       gprimd(3,idir)*eigvec(ii3+1)
     end do !idir
   end do !iatom

!  Apply each transformation to redvec and store at the correct location in redvtr (iatom -> jatom)
   do isym=1,nsym
     do iatom=1,natom
       jatom=symind(isym,iatom)
       iad1=3*(iatom-1)+1
       ii1=2*(iad1-1)+1
       iad2=3*(iatom-1)+2
       ii2=2*(iad2-1)+1
       iad3=3*(iatom-1)+3
       ii3=2*(iad3-1)+1
       do idir=1,3
         jad=3*(jatom-1)+idir
         jj=2*(jad-1)+1
         redvtr(jj)=dble(symrel(idir,1,isym))*redvec(ii1)+&
&         dble(symrel(idir,2,isym))*redvec(ii2)+&
&         dble(symrel(idir,3,isym))*redvec(ii3)
         redvtr(jj+1)=dble(symrel(idir,1,isym))*redvec(ii1+1)+&
&         dble(symrel(idir,2,isym))*redvec(ii2+1)+&
&         dble(symrel(idir,3,isym))*redvec(ii3+1)

       end do !idir
     end do !iatom

!    Compute redvtr in cartesian coordinates eigvtr
     do iatom=1,natom
       iad1=3*(iatom-1)+1
       ii1=2*(iad1-1)+1
       iad2=3*(iatom-1)+2
       ii2=2*(iad2-1)+1
       iad3=3*(iatom-1)+3
       ii3=2*(iad3-1)+1
       do idir=1,3
         jad=3*(iatom-1)+idir
         jj=2*(jad-1)+1
         eigvtr(jj)=rprimd(idir,1)*redvtr(ii1)+&
&         rprimd(idir,2)*redvtr(ii2)+&
&         rprimd(idir,3)*redvtr(ii3)
         eigvtr(jj+1)=rprimd(idir,1)*redvtr(ii1+1)+&
&         rprimd(idir,2)*redvtr(ii2+1)+&
&         rprimd(idir,3)*redvtr(ii3+1)
       end do !idir
     end do !iatom

!    Compute scalar product...
     symph(isym,imode)=0.0_dp
     do jad=1,3*natom
       jj=2*(jad-1)+1
       kk=2*3*natom*(imode-1)+2*(jad-1)+1
       symph(isym,imode)=symph(isym,imode)+eigvtr(jj)*eigvec(kk)+eigvtr(jj+1)*eigvec(kk+1)
     end do

   end do !isym
 end do !imode

!Treat degeneracies (different tolerances will be tried)
!Compute the order of the degeneracy, and
!attribute it to the lowest of the degenerate modes
!Also attribute the characters to the lowest mode
!When all the characters are integers, consider that the
!mode is non-degenerate. The maximum difference in frequency
!that is tolerated is on the order of 4cm-1 (which is large...)
 ABI_ALLOCATE(degeneracy,(3*natom))
 ABI_ALLOCATE(integer_characters,(3*natom))
 degeneracy(:)=1
 integer_characters(:)=0
 do itol=1,20
   ntol=itol
   do imode=3*natom,2,-1
     if(integer_characters(imode)==0)then
       do jmode=imode-1,1,-1
         if(integer_characters(jmode)==0)then
           if(abs(phfrq(imode)-phfrq(jmode))<itol*tol6)then
             degeneracy(jmode)=degeneracy(jmode)+degeneracy(imode)
             degeneracy(imode)=0
             symph(:,jmode)=symph(:,jmode)+symph(:,imode)
             symph(:,imode)=0.0_dp
           end if
         end if !integer_characters(jmode)==0
       end do !jmode
     end if !integer_characters(imode)==0
   end do !imode
   do imode=1,3*natom
     if(maxval(abs( symph(:,imode)-nint(symph(:,imode)) ))<0.05_dp)then
       integer_characters(imode)=1
     end if
   end do
   if(sum(integer_characters(:))==3*natom)exit
 end do !itol

!write(std_out,*)' dfpt_symph : degeneracy=',degeneracy(:)

 write(message,'(a,a,es8.2,a)')ch10,&
& ' Analysis of degeneracies and characters (maximum tolerance=',ntol*tol6,' a.u.)'
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')

 do imode=1,3*natom
   if(degeneracy(imode)/=0)then
     write(message,'(a,i4)') ' Symmetry characters of vibration mode #',imode
     call wrtout(iout,message,'COLL')
     call wrtout(std_out,message,'COLL')
     if(degeneracy(imode)>=2)then
       if(degeneracy(imode)==2) write(message,'(a,i4)') &
'        degenerate with vibration mode #',imode+1
       if(degeneracy(imode)>=3) write(message,'(a,i4,a,i4)') &
&       '       degenerate with vibration modes #',imode+1,' to ',imode+degeneracy(imode)-1
       call wrtout(iout,message,'COLL')
       call wrtout(std_out,message,'COLL')
     end if
     do jj=1,(nsym-1)/16+1
       write(message,'(16f5.1)') (symph(isym,imode),isym=(jj-1)*16+1,min(nsym,jj*16))
       call wrtout(iout,message,'COLL')
       call wrtout(std_out,message,'COLL')
     end do
   end if
 end do !imode

 ABI_DEALLOCATE(degeneracy)
 ABI_DEALLOCATE(integer_characters)
 ABI_DEALLOCATE(eigvtr)
 ABI_DEALLOCATE(redvtr)
 ABI_DEALLOCATE(redvec)
 ABI_DEALLOCATE(symph)
 ABI_DEALLOCATE(symind)

end subroutine dfpt_symph
!!***

!!****f* m_phonons/freeze_displ_allmodes
!!
!! NAME
!! freeze_displ_allmodes
!!
!! FUNCTION
!!  From a given set of phonon modes, generate and output supercells and
!!  displaced configurations of atoms.
!!  Typically useful to follow soft modes and see distorsions of crystal structures
!!
!! INPUTS
!! amu(ntypat) = mass of the atoms (atomic mass unit)
!! displ(2,3*natom,3*natom) = phonon mode displacements (complex)
!! freeze_displ = amplitude of the displacement to freeze into the supercell
!! natom = number of atoms in the unit cell
!! ntypat = number of atom types
!! phfrq(3*natom) = phonon frequencies
!! qphnrm = norm of phonon q vector (should be 1 or 0)
!! qphon = phonon wavevector
!! rprimd(3,3) = dimensionfull primitive translations in real space
!! typat(natom) = integer label of each type of atom (1,2,...)
!! xcart(3,natom) = cartesian coords of atoms in unit cell (bohr)
!!
!! OUTPUT
!! for the moment only prints to file, but could also return pointer to supercell object, with
!! rprimd and atomic positions, for further use
!!
!! NOTES
!! freeze_displ could be determined automatically from a temperature and the phonon frequency,
!! as the average displacement of the mode with a Bose distribution.
!!
!! PARENTS
!!      m_phonons
!!
!! CHILDREN
!!      destroy_supercell,freeze_displ_supercell,init_supercell_for_qpt
!!      prt_supercell_for_qpt
!!
!! SOURCE
!!

subroutine freeze_displ_allmodes(displ, freeze_displ, natom, outfile_radix, phfreq,  &
&         qphon, rprimd, typat, xcart, znucl)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom
 character(len=*),intent(in) :: outfile_radix
 real(dp), intent(in) :: freeze_displ
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: displ(2,3*natom,3*natom)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: phfreq(3*natom)
 real(dp),intent(in) :: qphon(3)
 real(dp),intent(in) :: xcart(3,natom)
 real(dp),intent(in) :: znucl(:)

! local vars
 integer :: jmode
 type(supercell_type) :: scell

! *************************************************************************

!determine supercell needed to freeze phonon
 call init_supercell_for_qpt(natom, qphon, rprimd, typat, xcart, znucl, scell)

 do jmode = 1, 3*natom
! reset positions
   scell%xcart = scell%xcart_ref

!  displace atoms according to phonon jmode
   call freeze_displ_supercell(displ(:,:,jmode), freeze_displ, scell)

!  print out everything for this wavevector and mode
   call prt_supercell_for_qpt (phfreq(jmode), jmode, outfile_radix, scell)
 end do

 call destroy_supercell (scell)

end subroutine freeze_displ_allmodes
!!***

end module m_phonons
!!***
