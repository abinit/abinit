!!****m* ABINIT/m_generate_training_set
!! NAME
!!  m_generate_training_set
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group ()
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_generate_training_set

 implicit none

 private
!!***

 public :: generate_training_set
!!***

contains
!!***

!!****f* ABINIT/generate_training_set
!!
!! NAME
!! generate_training_set
!!
!! FUNCTION
!!
!! INPUTS
!! acell(3)=length scales by which rprim is to be multiplied
!! amplitudes(namplitude) = list of the amplitudes of the unstable phonons
!!                          amplitudes(1:3,iamplitude) = qpt
!!                          amplitudes(4,iamplitude)   = mode
!!                          amplitudes(5,iamplitude)   = amplitude
!! natom= Number of atoms in the supercell
!! nconfig = Number of configuration to generate
!! hist<type(abihist)> = The history of the MD
!! namplitude = number of amplitude provided by the user
!! ngqpt(3)= Grid of qpoint in the DDB
!! nqshft=number of shift vectors in the repeated cell
!! option = option to deal with negative frequency -> Bose factor explodes (eg acoustic at Gamma)
!!          several philosophies to be implemented for the unstable modes:
!!          option == 1 =>  ignore
!!          option == 2 =>  populate them according to a default amplitude
!!          option == 3 =>  populate according to their modulus squared
!!          option == 4 =>  USER defined value(s), require namplitude and amplitude
!! qshft(3,nqshft)=vectors that will be used to determine
!! rlatt(3,3)= size of the supercell
!! rprimd(3,3)=dimensional primitive translations (bohr)
!! temperature_K =  temperature in Kelvin
!! xred(3,natom) = reduced dimensionless atomic coordinates in the supercell
!! comm=MPI communicator
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! hist<type abihist>=Historical record of positions
!!
!!
!! PARENTS
!!      multibinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine generate_training_set(acell,add_strain,amplitudes,filename,hist,natom,namplitude,nconfig,&
&                                ngqpt,nqshift,option,qshift,rlatt,rprimd,temperature_k,xred,comm,DEBUG)

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_strain
 use m_abihist, only : abihist,var2hist,abihist_findIndex
 use m_ifc, only : ifc_type,ifc_init_fromFile
 use m_crystal,     only : crystal_t
 use m_supercell, only : supercell_type
 use m_geometry, only : xcart2xred
 use m_phonons ,only :thermal_supercell_make,thermal_supercell_free
 use m_xmpi

!Arguments ------------------------------------
  !scalars
  integer,intent(in) :: natom,nconfig,nqshift,option,comm
  logical,intent(in) :: DEBUG
  real(dp),intent(in):: temperature_k
  integer,intent(in) :: namplitude
  logical,intent(in) :: add_strain
  !arrays
  integer,intent(in) :: ngqpt(3)
  integer,intent(in) :: rlatt(3,3)
  real(dp),intent(in) :: amplitudes(5,namplitude)
  real(dp),intent(in) :: qshift(3,nqshift)
  real(dp),intent(in) :: acell(3)
  real(dp), intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,natom)
  character(len=*),intent(in) :: filename
  type(abihist),intent(inout) :: hist
!Local variables-------------------------------
!scalar
  real(dp):: delta,rand
  integer :: ii,iconfig,natom_uc,direction
  character(len=500) :: message
  INTEGER                            :: n
  INTEGER                            :: i, iseed
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed

!arrays
  real(dp) :: dielt(3,3)
  real(dp),allocatable :: zeff(:,:,:)
  real(dp) :: acell_next(3),xred_next(3,natom),rprimd_next(3,3)
  type(ifc_type) :: ifc
  type(crystal_t) :: crystal
  type(supercell_type),allocatable :: thm_scells(:)
  type(strain_type) :: strain

! *************************************************************************

 ABI_UNUSED((/rprimd(1,1), xred(1,1)/))

 if ( .not. DEBUG ) then
    CALL RANDOM_SEED(size = n)
    ABI_ALLOCATE(seed,(n))
    seed =  iseed + (/ (i - 1, i = 1, n) /)

    CALL RANDOM_SEED(PUT = seed+xmpi_comm_rank(xmpi_world))
  end if

  write(message,'(a,(80a),a)') ch10,('=',ii=1,80),ch10
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')

  write(message, '(a,a,a)' )' Generation of all the configuration for the training set',ch10
  call wrtout(std_out,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  call ifc_init_fromFile(dielt,trim(filename),ifc,natom_uc,ngqpt,nqshift,qshift,crystal,zeff,comm)

  write(message, '(a,I0,a,f10.2,02a)' )' Generation of ',nconfig,' at the temperature ',&
&                            temperature_K,' K from the phonons',ch10

  call wrtout(std_out,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  ABI_DATATYPE_ALLOCATE(thm_scells,(nconfig))

  call thermal_supercell_make(amplitudes,crystal, Ifc,namplitude, nconfig, option,int(rlatt),&
&                             temperature_K, thm_scells)

  do iconfig = 1,nconfig
!   Get the atomic position for the new configuration
    call xcart2xred(thm_scells(iconfig)%natom,thm_scells(iconfig)%rprimd,&
&                   thm_scells(iconfig)%xcart,xred_next)
!   Just fill acell with the reference value, we apply strain on rprimd
    acell_next(:) = acell(:)

!   Get the rprim for the new configuration
    if(.not.add_strain)then
      rprimd_next(:,:) = thm_scells(iconfig)%rprimd(:,:)
    else
      call random_number(rand)
      direction = int(rand*6+1)
      call random_number(rand)
      delta = rand/500
      call random_number(rand)
      delta = delta*(two*rand-1)

      call strain_init(strain,delta=delta,direction=direction)
      call strain_apply(thm_scells(iconfig)%rprimd,rprimd_next,strain)
    end if

!   Fill history with the values of xred, acell and rprimd
    call var2hist(acell_next,hist,natom,rprimd_next,xred_next,DEBUG)
    hist%ihist = abihist_findIndex(hist,+1)

  end do

! Restart ihist before to leave
  hist%ihist = 1

  call ifc%free()
  call crystal%free()
  call thermal_supercell_free(nconfig,thm_scells)
  ABI_DATATYPE_DEALLOCATE(thm_scells)
  ABI_DEALLOCATE(zeff)

end subroutine generate_training_set
!!***

end module m_generate_training_set
!!***
