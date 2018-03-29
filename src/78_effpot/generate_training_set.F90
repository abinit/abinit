!{\src2tex{textfont=tt}}
!!****f* ABINIT/generate_training_set
!!
!! NAME
!! generate_training_set
!!
!! FUNCTION
!!
!! INPUTS
!! need to be update
!! OUTPUT
!! hist<type(abihist)> = The history of the MD
!!
!! PARENTS
!!      multibinit
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine generate_training_set(acell,filename,hist,natom,nconfig,ngqpt,nqshift,qshift,rlatt,rprimd,&
&                                temperature_k,xred,DEBUG)

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_xmpi
 use m_abihist, only : abihist,var2hist,abihist_findIndex
 use m_ifc, only : ifc_type,ifc_init_fromFile,ifc_free
 use m_crystal,     only : crystal_t,crystal_free
 use m_supercell, only : supercell_type
 use m_phonons ,only :thermal_supercell_make,thermal_supercell_free

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'generate_training_set'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
!End of the abilint section

  implicit none

!Arguments ------------------------------------
  !scalars
  integer,intent(in) :: natom,nconfig,nqshift
  logical,intent(in) :: DEBUG
  real(dp),intent(in):: temperature_k
  !arrays
  integer,intent(in) :: ngqpt(3)
  integer,intent(in) :: rlatt(3,3)
  real(dp),intent(in) :: qshift(3,nqshift)
  real(dp),intent(in) :: acell(3)
  real(dp), intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,natom)
  character(len=*),intent(in) :: filename
  type(abihist),intent(inout) :: hist
!Local variables-------------------------------
!scalar
  integer :: ii,ia,iconfig,mu,nu,natom_uc
  real(dp):: rand
  character(len=500) :: message
!arrays
  real(dp) :: dielt(3,3)
  real(dp),allocatable :: zeff(:,:,:)
  real(dp) :: acell_next(3),xred_next(3,natom),rprimd_next(3,3)
  type(ifc_type) :: ifc
  type(crystal_t) :: crystal
  type(supercell_type),allocatable :: thm_scells(:)

! *************************************************************************

  write(message,'(a,(80a),a)') ch10,('=',ii=1,80),ch10
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')

  write(message, '(a,a,a)' )' Generation of all the configuration for the training set',ch10
  call wrtout(std_out,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  call ifc_init_fromFile(dielt,trim(filename),ifc,natom_uc,ngqpt,nqshift,qshift,crystal,zeff)

  write(message, '(a,I0,a,f10.2,02a)' )' Generation of ',nconfig,' at the temperature ',&
&                            temperature_K,' K from the phonons',ch10

  call wrtout(std_out,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  ABI_DATATYPE_ALLOCATE(thm_scells,(nconfig))
  call thermal_supercell_make(crystal, Ifc, nconfig,int(rlatt), temperature_K, thm_scells)
  
  do iconfig = 1,nconfig
!   Get the atomic position for the new configuration
    call xcart2xred(thm_scells(iconfig)%natom,thm_scells(iconfig)%rprimd,&
&                   thm_scells(iconfig)%xcart,xred_next)
    
!   Just fill acell with the reference value, we apply strain on rprimd    
    acell_next(:) = acell(:)
    
!   Get the rprim for the new configuration
    rprimd_next(:,:) = thm_scells(iconfig)%rprimd(:,:)

    
!   Fill history with the values of xred, acell and rprimd
    call var2hist(acell_next,hist,natom,rprimd_next,xred_next,DEBUG)
    hist%ihist = abihist_findIndex(hist,+1)
    
  end do

! Restart ihist before to leave
  hist%ihist = 1


  call ifc_free(ifc)
  call crystal_free(crystal)
  call thermal_supercell_free(nconfig,thm_scells)
  ABI_DATATYPE_DEALLOCATE(thm_scells)
  ABI_DEALLOCATE(zeff)
  
end subroutine generate_training_set
!!***
