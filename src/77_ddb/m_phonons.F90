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
!! Copyright (C) 1999-2016 ABINIT group (XG,MG,MJV)
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
 use m_bz_mesh
 use m_atprj
 use m_sortph
 use m_ddb
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_fstrings,        only : itoa, ftoa, sjoin
 use m_numeric_tools,   only : simpson_int, wrap2_pmhalf
 use m_io_tools,        only : open_file
 use m_dynmat,          only : asria_corr,asrprs, gtdyn9
 use m_crystal,         only : crystal_t
 use m_ifc,             only : ifc_type, ifc_fourq
 use m_anaddb_dataset,  only : anaddb_dataset_type

 implicit none

 private

 public :: mkphbs     ! Compute phonon band structure
! TODO Write object to store the bands
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
   ! projected DOS (over atoms)

  real(dp),allocatable :: pjdos_int(:,:,:)
   ! pjdos_int(nomega,3,natom)
   ! Integrated atomic PJDOS

  real(dp),allocatable :: pjdos_type(:,:)
   ! pjdos_type(nomega,ntypat)
   ! phonon DOS contribution arising from a particular atom-type.

  real(dp),allocatable :: pjdos_type_int(:,:)
   ! pjdos_type_int(nomega,ntypat)
   ! Integrate phonon DOS contribution arising from a particular atom-type.

  real(dp),allocatable :: pjdos_rc_type(:,:,:)
   ! phdos(nomega,3,ntypat)
   ! phonon DOS contribution arising from a particular atom-type 
   ! decomposed along the three reduced directions.

 end type phonon_dos_type

 public :: mkphdos
 public :: phdos_print
 public :: phdos_print_debye
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
!!      anaddb,eph
!!
!! CHILDREN
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
 integer :: io,itype,unt,unt_by_atom,iatom
 real(dp) :: cfact
 character(len=500) :: msg
 character(len=fnlen) :: fname_by_atom
 character(len=3) :: unitname

! *************************************************************************

! === Convert everything into meV ===
 cfact=Ha_eV*1000 ; unitname='meV'
! === Leave everything in Ha      ===
! this should be the abinit default!
 cfact=one        ; unitname='Ha'

! === Open external file and write results ===
! TODO Here I have to rationalize how to write all this stuff!!
! 
 if (open_file(fname,msg,newunit=unt,form="formatted",action="write") /= 0) then
   MSG_ERROR(msg)
 end if

 fname_by_atom = trim(fname) // "_by_atom"
 if (open_file(fname_by_atom,msg,newunit=unt_by_atom,form="formatted",action="write") /= 0) then
   MSG_ERROR(msg)
 end if

 write(msg,'(3a)')'# ',ch10,'# Phonon density of states and atom type projected DOS'
 call wrtout(unt,msg,'COLL')
 write(msg,'(6a)')'# ',ch10,'# Energy in ',unitname,', DOS in states/',unitname
 call wrtout(unt,msg,'COLL')

 write(msg,'(3a)')'# ',ch10,'# Phonon density of states and atom projected DOS'
 call wrtout(unt_by_atom,msg,'COLL')
 write(msg,'(6a)')'# ',ch10,'# Energy in ',unitname,', DOS in states/',unitname
 call wrtout(unt_by_atom,msg,'COLL')

 select case (PHdos%prtdos)
 case (1)
   write(msg,'(a,es16.8,2a,i0)')&
&   '# Gaussian method with smearing = ',PHdos%dossmear*cfact,unitname,', nqibz =',PHdos%nqibz
 case (2) 
   write(msg,'(a,i0)')'# Tetrahedron method, nqibz= ',PHdos%nqibz
 case default
   write(msg,'(a,i0)')" Wrong prtdos = ",PHdos%prtdos
   MSG_ERROR(msg)
 end select
 call wrtout(unt,msg,'COLL')
 call wrtout(unt_by_atom,msg,'COLL')

 write(msg,'(5a)')'# ',ch10,'# omega     PHDOS    INT_PHDOS   PJDOS[atom_type=1]  INT_PJDOS[atom_type1] ...  ',ch10,'# '
 call wrtout(unt,msg,'COLL')
 write(msg,'(5a)')'# ',ch10,'# omega     PHDOS    PJDOS[atom=1]  PJDOS[atom2] ...  ',ch10,'# '
 call wrtout(unt_by_atom,msg,'COLL')

 do io=1,PHdos%nomega
   write(unt,'(3es17.8)',advance='NO')PHdos%omega(io)*cfact,PHdos%phdos(io)/cfact,PHdos%phdos_int(io)/cfact 
   write(unt_by_atom,'(2es17.8)',advance='NO')PHdos%omega(io)*cfact,PHdos%phdos(io)/cfact
   do itype=1,PHdos%ntypat
     write(unt,'(2es17.8,2x)',advance='NO')PHdos%pjdos_type(io,itype)/cfact,PHdos%pjdos_type_int(io,itype)/cfact
   end do 
   do iatom=1,PHdos%natom
     write(unt_by_atom,'(1es17.8,2x)',advance='NO') sum(PHdos%pjdos(io,1:3,iatom))/cfact
   end do 
   write(unt,*)
   write(unt_by_atom,*)
 end do

 close(unt)
 close(unt_by_atom)

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
 do io=1,PHdos%nomega
   om1dos(io) = PHdos%omega(io)    * PHdos%phdos(io)
   om2dos(io) = PHdos%omega(io)**2 * PHdos%phdos(io)
 end do

! integrate dos * omega
 intdos = zero
 call simpson_int(PHdos%nomega,PHdos%omega_step,om1dos,intdos)
 meanfreq = intdos(PHdos%nomega)

! integrate dos * omega^2
 intdos = zero
 call simpson_int(PHdos%nomega,PHdos%omega_step,om2dos,intdos)
 meanfreq2 = intdos(PHdos%nomega)

! Debye frequency = sqrt(<omega^2>* 3/2)
 debyefreq = sqrt (meanfreq2 * three * half)
 write (msg,'(a,E20.10,a,E20.10,a)') ' Debye frequency from DOS: ', debyefreq, ' (Ha) = ', debyefreq*Ha_THz, ' (THz)'
 call wrtout (ab_out,msg,"COLL")
 call wrtout (std_out,msg,"COLL")

! Debye temperature = hbar * Debye frequency / kb
 write (msg,'(a,E20.10,2a)') ' Debye temperature from DOS: ', debyefreq*Ha_K, ' (K)', ch10
 call wrtout (ab_out,msg,"COLL")
 call wrtout (std_out,msg,"COLL")

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
   if (abs(one-om2dos(1)/om2dos(io)) > 0.1_dp) then
     iomax = io
     exit
   end if
 end do

 avgom2dos = avgom2dos / (iomax-iomin+1)
! this value is also useful for partial atomic DOS, related to kinetic energy and Force constant in Moessbauer

 avgspeedofsound = (ucvol / 2 / pi**2 / avgom2dos)**third
 write (msg,'(a,E20.10,a,F16.4,2a)') '- Average speed of sound: ', avgspeedofsound, ' (at units) = ', &
&    avgspeedofsound * Bohr_Ang * 1.d-10 / Time_Sec, ' (m/s)',ch10
 call wrtout (ab_out,msg,"COLL")
 call wrtout (std_out,msg,"COLL")

! average force constant
 ABI_FREE(om2dos)
 ABI_FREE(om1dos)
 ABI_FREE(intdos)

end subroutine phdos_print_debye
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
!!
!! OUTPUT
!! PHdos<phonon_dos_type>=Container with phonon DOS, IDOS and atom-projected DOS.
!!
!! NOTES
!! On the use of the q-grids : 
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
!!      anaddb,eph
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkphdos(PHdos,Crystal,Ifc,prtdos,dosdeltae,dossmear,dos_ngqpt,dos_qshift)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkphdos'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_56_recipspace
 use interfaces_61_occeig
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: prtdos
 real(dp),intent(in) :: dosdeltae,dossmear
 type(crystal_t),intent(in) :: Crystal
 type(ifc_type),intent(in) :: Ifc
 type(phonon_dos_type),intent(inout) :: PHdos
!arrays
 integer,intent(in) :: dos_ngqpt(3)
 real(dp),intent(in) :: dos_qshift(3)

!Local variables -------------------------
!scalars
 integer,parameter :: brav1=1,chksymbreak0=0,bcorr0=0
 integer :: facbrv,iat,idir,imesh,imode,io,iq_ibz,itype,nkpt_fullbz
 integer :: nmesh,nqbz,nqpt_max,nqshft,option,timrev,ierr,natom,nomega
 real(dp) :: dum,gaussfactor,gaussprefactor,gaussval,low_bound,max_occ,pnorm 
 real(dp) :: qphnrm,upr_bound,xx,gaussmaxarg
 logical :: out_of_bounds
 character(len=500) :: msg
 character(len=80) :: errstr
 type(t_tetrahedron) :: tetraq
!arrays
 integer :: qptrlatt(3,3)
 integer,allocatable :: bz2ibz(:),ibz2bz(:),ngqpt(:,:)
 real(dp) :: displ(2*3*Crystal%natom*3*Crystal%natom)
 real(dp) :: eigvec(2,3,Crystal%natom,3*Crystal%natom),phfrq(3*Crystal%natom)
 real(dp) :: qlatt(3,3),qphon(3),rlatt(3,3)
 real(dp),allocatable :: dtweightde(:,:),full_eigvec(:,:,:,:,:),full_phfrq(:,:),Prf3D(:,:,:)
 real(dp),allocatable :: kpt_fullbz(:,:),qbz(:,:),qibz(:,:),qshft(:,:),tmp_phfrq(:),tweight(:,:)
 real(dp),allocatable :: qibz2(:,:),qshft2(:,:),wtq(:),wtq_folded(:),wtqibz(:)

! *********************************************************************

 DBG_ENTER("COLL")

 ! Consistency check.
 if (ALL(prtdos /= [1,2])) then 
   MSG_BUG(sjoin('prtdos should be 1 or 2, but received',itoa(prtdos)))
 end if 

 if (dosdeltae<=zero) then 
   MSG_BUG(sjoin('dosdeltae should be positive, but received',ftoa(dosdeltae)))
 end if 

 if (prtdos==1.and.dossmear<=zero) then 
   MSG_BUG(sjoin('dossmear should be positive but received',ftoa(dossmear)))
 end if 

 natom = Crystal%natom
 gaussmaxarg = sqrt(-log(1.d-90))

 ! Initialize container type, but with minimal values
 !call init_phondos(PHdos,Crystal%ntypat,natom,prtdos,1,1,1,smallest_real,greatest_real,dosdeltae,dossmear)
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

 ABI_MALLOC(PHdos%omega,(nomega))
 ABI_MALLOC(PHdos%phdos,(nomega))
 ABI_MALLOC(PHdos%phdos_int,(nomega))
 ABI_MALLOC(PHdos%pjdos,(nomega,3,natom))
 ABI_MALLOC(PHdos%pjdos_int,(nomega,3,natom))
 ABI_MALLOC(PHdos%pjdos_type,(nomega,crystal%ntypat))
 ABI_MALLOC(PHdos%pjdos_type_int,(nomega,crystal%ntypat))
 ABI_MALLOC(PHdos%pjdos_rc_type,(nomega,3,crystal%ntypat))
 !
 ! === Parameters defining the gaussian approximant ===
 if (prtdos==1) then 
   ! TODO: use dirac_delta and update reference files.
   gaussprefactor=one/(dossmear*sqrt(two_pi))
   gaussfactor=one/(sqrt2*dossmear)
   write(msg,'(4a,f8.5,2a,f8.5)')ch10,&
&   'mkphdos: calculating phonon DOS using gaussian method :',ch10,&
&   '   gaussian smearing [meV] = ',dossmear*Ha_meV,ch10,&
&   '   frequency step    [meV] = ',PHdos%omega_step*Ha_meV
 else if (prtdos==2) then 
   write(msg,'(2a)')ch10,'mkphdos: calculating phonon DOS using tetrahedron method'
 end if 
 call wrtout(std_out,msg,'COLL')
 !
 ! Initial lower and upper bound of the phonon spectrum.
 low_bound=real(huge(0))*half*PHdos%omega_step
 upr_bound=-low_bound
 !
 ! Save memory during the generation of the q-mesh in the full BZ  
 ! Take into account the type of Bravais lattice
 facbrv=1
 if (brav1==2) facbrv=2
 if (brav1==3) facbrv=4

 nmesh=2
 ABI_MALLOC(ngqpt,(3,nmesh))
 do imesh=1,nmesh

   if (imesh==1) then  
     ! Coarse q-mesh used during RF calculation.
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
     if (prtdos == 2) then
       ABI_MALLOC(qshft2,(3,nqshft))
       qshft2(:,:)=qshft(:,:)
     end if
   end if 

   nqpt_max=(ngqpt(1,imesh)*ngqpt(2,imesh)*ngqpt(3,imesh)*nqshft)/facbrv
   ABI_MALLOC(qibz,(3,nqpt_max))
   ABI_MALLOC(qbz,(3,nqpt_max))
   if (prtdos == 2 .and. imesh == 2) then
     ABI_MALLOC(qibz2,(3,nqpt_max))
   endif

   qptrlatt(:,:)=0
   qptrlatt(1,1)=ngqpt(1,imesh)
   qptrlatt(2,2)=ngqpt(2,imesh)
   qptrlatt(3,3)=ngqpt(3,imesh)
   option=1 
   !
   ! here I noticed a problem in the declaration of q1shft in the anaddb datatype 
   ! FIXME we write on unit std_out just to avoid problem with automatic tests
   call smpbz(brav1,std_out,qptrlatt,nqpt_max,nqbz,nqshft,option,qshft,qbz)
   !  
   !  Reduce the number of such points by symmetrization.
   ABI_MALLOC(ibz2bz,(nqbz))
   ABI_MALLOC(wtq,(nqbz))
   ABI_MALLOC(wtq_folded,(nqbz))
   wtq(:)=one/nqbz         ! Weights sum up to one
   timrev=1; option=1     ! TODO timrev should be input 
   !
   ! This call will set PHdos%nqibz
   call symkpt(chksymbreak0,Crystal%gmet,ibz2bz,std_out,qbz,nqbz,&
&    PHdos%nqibz,Crystal%nsym,Crystal%symrec,timrev,wtq,wtq_folded)
   write(std_out,*) 'PHdos%nqibz = ', PHdos%nqibz

   ABI_MALLOC(wtqibz,(PHdos%nqibz))
   do iq_ibz=1,PHdos%nqibz
     wtqibz(iq_ibz)=wtq_folded(ibz2bz(iq_ibz))
     qibz(:,iq_ibz)=qbz(:,ibz2bz(iq_ibz))
     if (prtdos==2 .and. imesh==2)   qibz2(:,:)=qibz(:,:)
   end do
   ABI_FREE(wtq_folded)
   ABI_FREE(qshft)

   if (prtdos==2.and.imesh==2) then
     !    
     ! Second mesh with tetrahedron method
     ! convert kptrlatt to double and invert, qlatt here refer to the shortest qpt vectors
     rlatt(:,:)=qptrlatt(:,:)
     call matr3inv(rlatt,qlatt)

     ABI_MALLOC(qshft,(3,nqshft))
     !qshft(:,1)=inp%q2shft(:)  ! FIXME small inconsistency in the dimension of q1shft
     qshft(:,1)= dos_qshift(:)
     ABI_FREE(qshft2)
     ABI_MALLOC(qshft2,(3,nqshft))
     qshft2(:,:)=qshft(:,:)
     nkpt_fullbz=nqbz 
     ABI_MALLOC(bz2ibz,(nkpt_fullbz))
     ABI_MALLOC(kpt_fullbz,(3,nkpt_fullbz))
     !    
     ! Make full kpoint grid and get equivalence to irred kpoints.
     ! This routines scales **very badly** wrt nkpt_fullbz, should introduce check on the norm.
     call get_full_kgrid(bz2ibz,qibz,kpt_fullbz,qptrlatt,PHdos%nqibz,&
&      nkpt_fullbz,nqshft,Crystal%nsym,qshft,Crystal%symrel)
     !    
     ! Get tetrahedra, ie indexes of the full q-points at their summits
     call init_tetra(bz2ibz, crystal%gprimd, qlatt, kpt_fullbz, nqbz, tetraq, ierr, errstr)
     ABI_CHECK(ierr==0,errstr)

     ABI_FREE(qshft)
     ABI_FREE(bz2ibz)
     ABI_FREE(kpt_fullbz)
     !    
     ! Allocate arrays used to store the entire spectrum, Required to calculate tetra weights.
     ABI_MALLOC(full_phfrq,(3*natom,PHdos%nqibz))
     ABI_STAT_MALLOC(full_eigvec,(2,3,natom,3*natom,PHdos%nqibz), ierr)
     ABI_CHECK(ierr==0, 'out-of-memory in full_eigvec')

     !ABI_MALLOC(Prf3D,(3*natom,PHdos%nqibz,1))
   end if  ! prtdos==2.and.imesh==2
   !    
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
     !
     ! Frequency mesh.
     PHdos%omega_min=low_bound; if (ABS(PHdos%omega_min)<tol5) PHdos%omega_min=-tol5
     PHdos%omega_max=upr_bound 
     PHdos%nomega=NINT((PHdos%omega_max-PHdos%omega_min)/PHdos%omega_step)+1
     PHdos%nomega=MAX(6,PHdos%nomega) ! Ensure Simpson integration will be ok

     ABI_MALLOC(PHdos%omega,(PHdos%nomega))
     do io=1,PHdos%nomega
       PHdos%omega(io)=PHdos%omega_min+PHdos%omega_step*(io-1)
     end do

     if (imesh/=1) then 
       write(std_out,*)&
        'nomega = ',PHdos%nomega,' omega_min [cm-1] =',PHdos%omega_min*Ha_cmm1,' omega_max [cm-1] =',PHdos%omega_max*Ha_cmm1
     end if 

     ABI_CALLOC(PHdos%phdos,(PHdos%nomega))
     ABI_CALLOC(PHdos%pjdos,(PHdos%nomega,3,natom))
     !    
     ! === Sum over irreducible q-points ===
     do iq_ibz=1,PHdos%nqibz
       qphon(:)=qibz(:,iq_ibz); qphnrm=one

       ! Fourier interpolation.
       call ifc_fourq(Ifc,Crystal,qphon,phfrq,displ,out_eigvec=eigvec)

       !if (prtdos==2.and.imesh==2) Prf3D(:,iq_ibz,1)=phfrq(:)
       
       dum=MINVAL(phfrq); PHdos%omega_min=MIN(PHdos%omega_min,dum)
       dum=MAXVAL(phfrq); PHdos%omega_max=MAX(PHdos%omega_max,dum)
       out_of_bounds = (PHdos%omega_min<low_bound .or. PHdos%omega_max>upr_bound) 

       if (imesh>1.and..not.out_of_bounds) then
         select case (prtdos)
         case (1) 
           !
           ! Accumulate PHDOS and PJDOS with gaussian method.
           do imode=1,3*natom 
             do io=1,PHdos%nomega
               xx=(PHdos%omega(io)-phfrq(imode))*gaussfactor
               gaussval = zero
               if(abs(xx) < gaussmaxarg) gaussval=gaussprefactor*exp(-xx*xx)
               PHdos%phdos(io)=PHdos%phdos(io) + wtqibz(iq_ibz)*gaussval
               do iat=1,natom
                 do idir=1,3
                   pnorm=eigvec(1,idir,iat,imode)**2+eigvec(2,idir,iat,imode)**2
                   PHdos%pjdos(io,idir,iat)=PHdos%pjdos(io,idir,iat)+ pnorm*wtqibz(iq_ibz)*gaussval
                 end do
               end do
             end do 
           end do 

         case (2) 
           ! === Tetrahedrons ===
           !  * Save phonon frequencies and eigenvectors. 
           !  * Sum is done after the loops over the two meshes.
           full_phfrq(:,iq_ibz)=phfrq(:)
           full_eigvec(:,:,:,:,iq_ibz)=eigvec
           !if (prtdos==2.and.imesh==2) Prf3D(:,iq_ibz,1)=phfrq(:)
         case default
           write(msg,'(a,i0)')" Wrong value for prtdos= ",prtdos
           MSG_ERROR(msg)
         end select
       end if !Second mesh and not out of boundaries
       !
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

   ABI_FREE(ibz2bz)
   ABI_FREE(qibz)
   ABI_FREE(qbz)
   ABI_FREE(wtq)
   ABI_FREE(wtqibz)
 end do !imesh
 ABI_FREE(ngqpt)

 if (allocated(PHdos%phdos_int)) then
   ABI_FREE(PHdos%phdos_int)
 end if
 if (allocated(PHdos%pjdos_int)) then
   ABI_FREE(PHdos%pjdos_int)
 end if

 ABI_CALLOC(PHdos%phdos_int,(PHdos%nomega))

 if (prtdos==2) then 
   ! === Integrate using tetrahedrons ===
   !  * All the data are contained in full_phfrq and full_eigvec. 
   !  * low_bound and upr_bound contain the entire spectrum calculated on the dense mesh. 
   ABI_MALLOC(tmp_phfrq,(PHdos%nqibz))
   ABI_MALLOC(tweight,(PHdos%nqibz,PHdos%nomega))
   ABI_MALLOC(dtweightde,(PHdos%nqibz,PHdos%nomega))
   ABI_MALLOC(PHdos%pjdos_int,(PHdos%nomega,3,natom))
   PHdos%phdos=zero; PHdos%pjdos=zero; PHdos%pjdos_int=zero
   max_occ=one 

   do imode=1,3*natom 
     tmp_phfrq(:)=full_phfrq(imode,:)
     !    
     ! Calculate general integration weights at each irred kpoint as in Blochl et al PRB 49 16223.
     call get_tetra_weight(tmp_phfrq, low_bound, upr_bound,max_occ, PHdos%nomega, PHdos%nqibz, tetraq, bcorr0,&
&      tweight,dtweightde,xmpi_comm_self)

     do io=1,PHdos%nomega
       do iq_ibz=1,PHdos%nqibz
         PHdos%phdos(io)=PHdos%phdos(io)+dtweightde(iq_ibz,io)
         PHdos%phdos_int(io)=PHdos%phdos_int(io)+tweight(iq_ibz,io)
         do iat=1,natom
           do idir=1,3
             pnorm=full_eigvec(1,idir,iat,imode,iq_ibz)**2 + full_eigvec(2,idir,iat,imode,iq_ibz)**2
             PHdos%pjdos(io,idir,iat)=PHdos%pjdos(io,idir,iat) + pnorm*dtweightde(iq_ibz,io)
             PHdos%pjdos_int(io,idir,iat)=PHdos%pjdos_int(io,idir,iat) + pnorm*tweight(iq_ibz,io)         
           end do
         end do
       end do
     end do

   end do 
   ABI_FREE(tmp_phfrq)
   ABI_FREE(tweight)
   ABI_FREE(dtweightde)
 end if 
 !
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

 ABI_CALLOC(PHdos%pjdos_rc_type,(PHdos%nomega,3,Crystal%ntypat))
 ABI_CALLOC(PHdos%pjdos_type,(PHdos%nomega,Crystal%ntypat))
 ABI_CALLOC(PHdos%pjdos_type_int,(PHdos%nomega,Crystal%ntypat))

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
 !
 ! Evaluate IDOS using simple simpson integration
 ! TODO should avoid the simpson rule using derf.F90, just to be consistent
 if (prtdos==1) then 
   call simpson_int(PHdos%nomega,PHdos%omega_step,PHdos%phdos,PHdos%phdos_int)
   do itype=1,Crystal%ntypat
     call simpson_int(PHdos%nomega,PHdos%omega_step,PHdos%pjdos_type(:,itype),PHdos%pjdos_type_int(:,itype))
   end do
 end if 

!output phonon isosurface
! MG Commented this call on  Wed Jul 16 2014
! TODO: The calculation of the isosurface should be done in a specialized routine
 if (.False. .and. prtdos==2) then
   call printbxsf(Prf3D,zero,zero,crystal%gprimd,qptrlatt,3*natom,&
     PHdos%nqibz,qibz2,Crystal%nsym,.FALSE.,Crystal%symrec,Crystal%symafm,.TRUE.,1,qshft2,nqshft,'Phfrq3D',ierr)
   ABI_FREE(Prf3D)
 end if

 if (allocated(qibz2)) then
   ABI_FREE(qibz2)
 end if
 if (allocated(qshft2)) then
   ABI_FREE(qshft2)
 end if

 if (prtdos==2) then
   ABI_FREE(full_phfrq)
   ABI_FREE(full_eigvec)
   call destroy_tetra(tetraq)
 end if

 DBG_EXIT("COLL")

end subroutine mkphdos
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
!!      anaddb
!!
!! CHILDREN
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
 type(phonon_dos_type),target,intent(in) :: phdos
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
   nctkarr_t('pjdos_rc_type', "dp", 'number_of_frequencies, three, number_of_atom_species')])
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
!! dielt(3,3)=dielectric tensor
!! tcpui=initial cpu time
!! twalli=initial wall clock time
!! zeff(3,3,natom)=effective charge on each atom, versus electric field and atomic displacement
!! comm=MPI communicator
!!
!! OUTPUT
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkphbs(Ifc,Crystal,inp,ddb,d2asr,outfile_radix,singular,tcpui,twalli,uinvers,vtinvers,zeff,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkphbs'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_72_response
 use interfaces_77_ddb
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,optional,intent(in) :: comm
 real(dp),intent(in) :: tcpui,twalli
 character(len=*),intent(in) :: outfile_radix
 type(ifc_type),intent(in) :: Ifc
 type(crystal_t),intent(in) :: Crystal
 type(anaddb_dataset_type),target,intent(in) :: inp
 type(ddb_type),intent(in) :: ddb
!arrays
 real(dp),intent(in) :: zeff(3,3,ddb%natom)
 real(dp),intent(inout) :: singular(1:3*ddb%natom*(3*ddb%natom-1)/2)
 real(dp),intent(inout) :: uinvers(1:3*ddb%natom*(3*ddb%natom-1)/2,1:3*ddb%natom*(3*ddb%natom-1)/2)
 real(dp),intent(inout) :: vtinvers(1:3*ddb%natom*(3*ddb%natom-1)/2,1:3*ddb%natom*(3*ddb%natom-1)/2)
 real(dp),intent(inout) :: d2asr(2,3,ddb%natom,3,ddb%natom)

!Local variables -------------------------
!scalars
 integer,parameter :: master=0
 integer :: iphl1,iblok,rftyp, ii,nfineqpath,nsym,mpert,natom,ncid
 real(dp) :: tcpu,twall,res
 character(len=fnlen) :: tmpfilename
 character(500) :: msg
!arrays
 integer :: rfphon(4),rfelfd(4),rfstrs(4)
 integer,allocatable :: ndiv(:)
 real(dp) :: qphnrm(3), qphon(3), qphon_padded(3,3)
 real(dp) :: d2cart(2,ddb%msize),real_qphon(3) 
 real(dp) :: displ(2*3*ddb%natom*3*ddb%natom),eigval(3,ddb%natom)
 real(dp),allocatable :: phfrq(:),eigvec(:,:,:,:,:),save_phfrq(:,:),save_phdispl_cart(:,:,:,:),save_qpoints(:,:)
 real(dp),allocatable,target :: alloc_path(:,:)
 real(dp),pointer :: fineqpath(:,:)
 type(atprj_type) :: atprj

! *********************************************************************

 ! Only master works for the time being
 if(present(comm)) then
   if (xmpi_comm_rank(comm) /= master) return
 end if

 nsym = Crystal%nsym; natom = Crystal%natom
 mpert = ddb%mpert

 nullify(fineqpath)
 nfineqpath = inp%nph1l
 fineqpath => inp%qph1l

 if(inp%nph1l==0) then
   if (inp%nqpath==0) then
     return ! if there is nothing to do, return
   else
     ! allow override of nph1l with nqpath if the former is not set
     ! allocatae and compute path here and make fineqpath points to it
     ABI_MALLOC(ndiv,(inp%nqpath-1))
     call make_path(inp%nqpath,inp%qpath,Crystal%gmet,'G',inp%ndivsm,ndiv,nfineqpath,alloc_path,std_out)
     ABI_FREE(ndiv)
     fineqpath => alloc_path
   end if
 end if

 write(msg, '(a,(80a),a,a,a,a)' ) ch10,('=',ii=1,80),ch10,ch10,' Treat the first list of vectors ',ch10
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

 call timein(tcpu,twall)
 write(msg, '(a,f11.3,a,f11.3,a)' )'-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

 if (inp%natprj_bs > 0) then
   call atprj_init(atprj, natom, inp%natprj_bs, inp%iatprj_bs, outfile_radix)
 end if

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
   if(inp%ifcflag==1)then

     ! Get d2cart using the interatomic forces and the
     ! long-range coulomb interaction through Ewald summation
     call gtdyn9(ddb%acell,Ifc%atmfrc,Ifc%dielt,Ifc%dipdip,Ifc%dyewq0,d2cart,Crystal%gmet,ddb%gprim,mpert,natom,&
&     Ifc%nrpt,qphnrm(1),qphon,Crystal%rmet,ddb%rprim,Ifc%rpt,Ifc%trans,Crystal%ucvol,Ifc%wghatm,Crystal%xred,zeff)

   else if(inp%ifcflag==0)then

     ! Look for the information in the DDB (no interpolation here!)
     rfphon(1:2)=1
     rfelfd(1:2)=0
     rfstrs(1:2)=0
     rftyp=inp%rfmeth
     qphon_padded = zero
     qphon_padded(:,1) = qphon

     call gtblk9(ddb,iblok,qphon_padded,qphnrm,rfphon,rfelfd,rfstrs,rftyp)

     ! Copy the dynamical matrix in d2cart
     d2cart(:,1:ddb%msize)=ddb%val(:,:,iblok)

     ! Eventually impose the acoustic sum rule based on previously calculated d2asr
     select case (inp%asr)
     case (0)
       continue
     case (1,2,5)
       call asria_corr(inp%asr,d2asr,d2cart,mpert,natom)
     case (3,4)
       ! Impose acoustic sum rule plus rotational symmetry for 0D and 1D systems
       call asrprs(inp%asr,2,3,uinvers,vtinvers,singular,d2cart,mpert,natom,crystal%xcart)
     case default
       write(msg,'(a,i0)')"Wrong value for asr: ",inp%asr
       MSG_ERROR(msg)
     end select
   end if

   !  Calculation of the eigenvectors and eigenvalues of the dynamical matrix
   call dfpt_phfrq(ddb%amu,displ,d2cart,eigval,eigvec,Crystal%indsym,&
&   mpert,Crystal%nsym,natom,nsym,Crystal%ntypat,phfrq,qphnrm(1),qphon,&
&   crystal%rprimd,inp%symdynmat,Crystal%symrel,Crystal%symafm,Crystal%typat,Crystal%ucvol)


   !call ifc_fourq(Ifc,Crystal,qphon,phfrq,displ,out_eigvec=eigvec)

   if (abs(inp%freeze_displ) > tol10) then
     real_qphon = zero
     if (abs(qphnrm(1)) > tol8) then
       real_qphon = qphon / qphnrm(1)
     end if
     call freeze_displ_allmodes(displ, inp%freeze_displ, natom, outfile_radix, phfrq, &
&     real_qphon, crystal%rprimd, Crystal%typat, crystal%xcart, crystal%znucl)
   end if

   ! If requested, output projection of each mode on given atoms
   if (inp%natprj_bs > 0) then
     call atprj_print(atprj, iphl1, phfrq, eigvec)
   end if

   ! In case eivec == 4, write output files for band2eps (visualization of phonon band structures)
   if (inp%eivec == 4) then
     tmpfilename = trim(outfile_radix)//"_B2EPS"
     call sortph(eigvec,displ,tmpfilename,natom,phfrq)
   end if

   ! Write the phonon frequencies
   call dfpt_prtph(displ,inp%eivec,inp%enunit,ab_out,natom,phfrq,qphnrm(1),qphon)

   save_phfrq(:,iphl1) = phfrq
   save_phdispl_cart(:,:,:,iphl1) = RESHAPE(displ,(/2, 3*natom, 3*natom/))

   ! Determine the symmetries of the phonon mode at Gamma
   ! TODO: generalize for other q-point little groups.
   if(sum(abs(qphon(:)))<DDB_QTOL)then
     call dfpt_symph(ab_out,ddb%acell,eigvec,Crystal%indsym,natom,nsym,phfrq,ddb%rprim,Crystal%symrel)
   end if

   ! if we have an acoustic mode (small q and acoustic type displacements)
   ! extrapolate speed of sound in this direction, and Debye frequency
   call wrap2_pmhalf(qphon(1),real_qphon(1),res)
   call wrap2_pmhalf(qphon(2),real_qphon(2),res)
   call wrap2_pmhalf(qphon(3),real_qphon(3),res)

   if (sqrt(real_qphon(1)**2+real_qphon(2)**2+real_qphon(3)**2) < quarter .and. &
&   sqrt(real_qphon(1)**2+real_qphon(2)**2+real_qphon(3)**2) > tol6) then
     call prtvsound(ab_out,eigvec, Crystal%gmet, natom, phfrq, real_qphon, Crystal%ucvol)
   end if

 end do ! iphl1

!deallocate sortph array
 call end_sortph()

 if (inp%natprj_bs > 0) then
   call atprj_destroy(atprj)
 end if

#ifdef HAVE_NETCDF
 tmpfilename = trim(outfile_radix)//"_PHBST.nc"
 NCF_CHECK_MSG(nctk_open_create(ncid, tmpfilename, xmpi_comm_self), "Creating PHBST")
 NCF_CHECK(crystal_ncwrite(Crystal, ncid))
 call phonons_ncwrite(ncid,natom,nfineqpath,save_qpoints,(/(one, iphl1=1,nfineqpath)/),save_phfrq,save_phdispl_cart)
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('atomic_mass_units', "dp", "number_of_atom_species")],defmode=.True.))
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'atomic_mass_units'), ddb%amu))
 NCF_CHECK(nf90_close(ncid))
#endif

 call phonons_write(natom,nfineqpath,save_qpoints,(/(one, iphl1=1,nfineqpath)/),save_phfrq,save_phdispl_cart)

! call phonons_writeEPS(natom,nfineqpath,Crystal%ntypat,save_qpoints,Crystal%typat,(/(one, iphl1=1,nfineqpath)/),save_phfrq,save_phdispl_cart)

 ABI_FREE(save_qpoints)
 ABI_FREE(save_phfrq)
 ABI_FREE(save_phdispl_cart)
 ABI_FREE(phfrq)
 ABI_FREE(eigvec)

 if (allocated(alloc_path)) then 
   ABI_FREE(alloc_path)
 end if

contains 
!!***

!!****f* m_phonons/prtvsound
!!
!! NAME
!! prtvsound
!!
!! FUNCTION
!!  From the frequencies for acoustic modes at small q, estimate speed of sound and Debye temperature
!!
!! INPUTS
!! unit=Fortran unit number
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

subroutine prtvsound(unit,eigvec,gmet,natom,phfrq,qphon,ucvol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prtvsound'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalras
 integer, intent(in) :: natom,unit
 real(dp), intent(in) :: ucvol
!arrays
 real(dp), intent(in) :: gmet(3,3),qphon(3)
 real(dp), intent(in) :: phfrq(3*natom),eigvec(2,3*natom,3*natom)

!Local variables -------------------------
 integer :: iatref,imode, iatom, isacoustic
 character(len=500) :: msg
 real(dp) :: qnormcart, speedofsound, tdebye
 real(dp) :: qtmp(3)

! *********************************************************************

 do imode = 1, 3*natom
   
!  Check if this mode is acoustic like: scalar product of all displacement vectors are collinear
   isacoustic = 1
!  Find reference atom with non-zero displacement
   do iatom=1,natom
     if(sum(eigvec(:,(iatom-1)*3+1:(iatom-1)*3+3,imode)**2) >tol16)iatref=iatom
   enddo
!  Now compute scalar product, and check they are all positive
   do iatom = 1, natom
     if (sum(eigvec(:,(iatom-1)*3+1:(iatom-1)*3+3, imode)*eigvec(:,(iatref-1)*3+1:(iatref-1)*3+3, imode)) < tol16 ) isacoustic = 0
   end do
   if (isacoustic == 0) cycle

   write (msg, '(a,I6,a,3F12.4)') ' Found acoustic mode ', imode, ' for |q| in red coord < 0.25 ; q = ', qphon
   call wrtout(unit,msg,'COLL')
   call wrtout(std_out,msg,'COLL')

   qtmp = matmul(gmet, qphon)
   qnormcart = two * pi * sqrt(sum(qphon*qtmp))
   speedofsound = phfrq(imode) / qnormcart

!  from phonon frequency, estimate speed of sound by linear interpolation from Gamma
   write (msg, '(2a,a,E20.10,a,a,F20.5)') &
&   ' Speed of sound for this q and mode:',ch10,&
&   '   in atomic units: ', speedofsound, ch10,&
&   '   in SI units m/s: ', speedofsound * Bohr_Ang * 1.d-10 / Time_Sec
   call wrtout(unit,msg,'COLL')

!  also estimate partial Debye temperature, = energy if this band went to zone edge
   tdebye = speedofsound * pi * (six / pi / ucvol)**(third)
   write (msg, '(2a,a,E20.10,a,a,F20.5)') &
&   ' Partial Debye temperature for this q and mode:',ch10,&
&   '   in atomic units: ', tdebye, ch10,&
&   '   in SI units K  : ', tdebye * Ha_K
   call wrtout(unit,msg,'COLL')

   call wrtout(unit,"",'COLL')
 end do

end subroutine prtvsound
!!***

end subroutine mkphbs
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
 real(dp),target,intent(in) :: qpoints(3,nqpts),weights(nqpts)
 real(dp),target,intent(in) :: phfreq(3*natom,nqpts),phdispl_cart(2,3*natom,3*natom,nqpts)

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

!!****f* m_phonons/phonons_write
!! NAME
!! phonons_write
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

 subroutine phonons_write(natom,nqpts,qpoints,weights,phfreq,phdispl_cart)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phonons_write'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nqpts
!arrays
 real(dp),intent(in) :: qpoints(3,nqpts),weights(nqpts)
 real(dp),intent(in) :: phfreq(3*natom,nqpts)
 real(dp),intent(in) :: phdispl_cart(2,3*natom,3*natom,nqpts)

!Local variables-------------------------------
!scalars
 integer :: nphmodes, iq, iunit
! integer :: imod, icomp
 real(dp) :: dummy
 character(len=300) :: formt
 character(len=500) :: msg

! *************************************************************************

 nphmodes = 3*natom

 dummy = qpoints(1,1)
 dummy = weights(1)
 dummy = phdispl_cart(1,1,1,1)

 if (open_file("PHFRQ", msg, newunit=iunit, form="formatted", status="unknown", action="write") /= 0) then
   MSG_ERROR(msg)
 end if

 write (iunit, '(a)')  '# ANADDB generated phonon band structure file. All in Ha atomic units'
 write (iunit, '(a)')  '# '
 write (iunit, '(a,i0)')  '# number_of_qpoints ', nqpts
 write (iunit, '(a,i0)')  '# number_of_phonon_modes ', nphmodes
 write (iunit, '(a)')  '# '
 
 write (formt,'(a,i0,a)') "(I5, ", nphmodes, "E20.10)"

 do iq = 1, nqpts
   write (iunit, formt)  iq, phfreq(:,iq)
 end do

 close(iunit)

! if (open_file("PHDISPL", msg, unit=iunit, form="formatted", status="unknown", action="write") /= 0) then
!   MSG_ERROR(msg)
! end if
!
! write (iunit, '(a)')     '# ANADDB generated phonon displacements, along points in PHFRQ file. All in Ha atomic units'
! write (iunit, '(a)')     '# '
! write (iunit, '(a)')     '# displacements in cartesian coordinates, Re and Im parts '
! write (iunit, '(a,I5)')  '# number_of_qpoints ', nqpts
! write (iunit, '(a,I5)')  '# number_of_phonon_modes ', nphmodes
! write (iunit, '(a)')     '# '
! 
! !write (formt,'(a,I3,a)') "( ", nphmodes, "(2E20.10,2x))"
! formt = "(2E20.10,2x)"
!
! do iq = 1, nqpts
!   write (iunit, '(a, I5)') '# iq ', iq
!   do imod = 1, nphmodes
!     write (iunit, '(a, I5)') '# imode ', imod
!     do icomp = 1, nphmodes
!       write (iunit, formt, ADVANCE='NO') phdispl_cart(:,icomp,imod,iq)
!     end do
!     write (iunit, '(a)') ' '
!   end do
! end do
!
! close(iunit)


end subroutine phonons_write
!!***

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
!no_abirules
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
 displcpx = cmplx(phdispl_cart(1,:,:,:),phdispl_cart(2,:,:,:))

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

end module m_phonons
!!***
