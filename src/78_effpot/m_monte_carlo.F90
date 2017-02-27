!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_monte_carlo
!!
!! NAME
!! m_monte_carlo
!!
!! FUNCTION
!! Module for using a effective potential
!! Container type is defined, and destruction, print subroutines 
!! 
!! COPYRIGHT
!! Copyright (C) 2010-2016 ABINIT group (AM)
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

module m_monte_carlo

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_xmpi
 use m_strain
 use m_OurRng
 use m_ifc
 use m_effective_potential
 use m_multibinit_dataset
#ifdef HAVE_NETCDF
 use netcdf
#endif

 implicit none

 private

 public :: monte_carlo_init
 public :: monte_carlo_run
 public :: monte_carlo_step
!!*** 

!!****t* m_monte_carlo/monte_carlo_type
!! NAME
!! monte_carlo_type
!!
!! FUNCTION
!! structure for a effective potential constructed.
!!
!! SOURCE

 type, public :: monte_carlo_type

   integer :: natom                                 
!     Number of atoms in primitive cell

   integer :: natom_supercell          
!     Number of atoms in supercell


   integer, allocatable :: cell(:,:) 
!    cell(nrpt, 3) 
!    Give the relation between the number of the cell 
!    and the indeces of cell 

   real(dp), allocatable :: xcart(:,:)
!    xcart(3, natom) 
!    positions of atoms in cartesian coordinates

   type(ifc_type) :: ifcs
!   type with ifcs constants (short + ewald) 
!   also contains the number of cell and the indexes

   type(ifc_type),dimension(:),allocatable :: phonon_strain_coupling
!   Array of ifc with phonon_strain_coupling coupling for each strain 

 end type monte_carlo_type
!!***


CONTAINS  !===========================================================================================

!!****f* m_monte_carlo/monte_carlo_init
!! NAME
!! monte_carlo_init
!!
!! FUNCTION
!!  initilisation of monte carlo
!!
!! PARENTS
!!
!! CHILDREN
!!      effective_potential_getdeltaenergy,hist2var,metric,random_number
!!      var2hist,xred2xcart
!!
!! SOURCE

subroutine monte_carlo_init(eff_pot)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'monte_carlo_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(effective_potential_type), intent(out) :: eff_pot
!arrays
!Local variables-------------------------------
!scalar
!arrays

! *************************************************************************
end subroutine monte_carlo_init
!!***


!!****f* m_monte_carlo/monte_carlo_run
!! NAME
!! monte_carlo_init
!!
!! FUNCTION
!!  initilisation of monte carlo
!!
!! PARENTS
!!      
!!
!! CHILDREN
!!
!! SOURCE
subroutine monte_carlo_run(reference_effective_potential)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'monte_carlo_run'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(effective_potential_type), intent(inout) :: reference_effective_potential
!arrays
!Local variables-------------------------------
!scalar
 integer :: ii
 integer(8) :: seed
 real(dp) :: energy
!arrays
 real(dp) :: strain(6),jj
 real(dp) :: disp(3)
 character(len=500) :: message

! *************************************************************************
 write(message, '(a,a,(80a),a,a)' ) ch10,('=',ii=1,80),ch10,' Begin Monte Carlo',ch10
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')
 
 do ii=1,1
   CALL SYSTEM_CLOCK(count=seed)
   strain(:) = zero
   call OurRng(seed,jj)
   call OurRng(seed,strain(1))
   strain(1:3) = (strain(1) + one)*(-1)**ANINT(jj*1000)
   disp(1) = one
   
!   call effective_potential_getEnergy(reference_effective_potential,energy,&
!&    strain1=strain)

   write(message,*) 'Iteration:',ii,' Total Energy:',energy
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
   
   if(energy < reference_effective_potential%energy) then
     reference_effective_potential%energy = energy
     reference_effective_potential%internal_stress= strain
     write(message,'(a,a)') ' Configuration accepted',ch10
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')
   else
     write(message,'(a,a)') ' Configuration refused',ch10
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')
   end if

 end do
 
 write(message, '(a,a,(80a),a,a)' ) ch10,' End Monte Carlo',ch10,('=',ii=1,80),ch10
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

end subroutine monte_carlo_run
!!***  


!!****f* m_monte_carlo/monte_carlo_step
!! NAME
!! monte_carlo_step
!!
!! INPUTS
!! ab_mover <type(abimover)> : Datatype with all the information
!!                                needed by the preditor
!! itime  : Index of the present iteration
!! ntime  : Maximal number of iterations
!! zDEBUG : if true print some debugging information
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! hist <type(abihist)> : History of positions,forces
!!                               acell, rprimd, stresses
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      effective_potential_getdeltaenergy,hist2var,metric,random_number
!!      var2hist,xred2xcart
!!
!! SOURCE

subroutine monte_carlo_step(ab_mover,eff_pot,hist,itime,ntime,zDEBUG,iexit)

 use m_abimover
 use m_abihist
 use m_effective_potential

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'monte_carlo_step'
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itime
 integer,intent(in) :: ntime
 integer,intent(in) :: iexit
 logical,intent(in) :: zDEBUG
 type(abimover),intent(in)       :: ab_mover
 type(abihist),intent(inout) :: hist
 type(effective_potential_type),intent(in):: eff_pot
!Local variables-------------------------------
!scalars
 integer  ::  kk,ia,mu
 real(dp) ::  acc,delta,ucvol
 real(dp),parameter :: v2tol=tol8
 real(dp) :: de,etotal
!arrays
 real(dp),allocatable,save :: displacement(:,:)

 real(dp) :: acell(3),rprim(3,3),rprimd(3,3)
 real(dp) :: gprimd(3,3),gmet(3,3),rmet(3,3),fcart(3,ab_mover%natom)
 real(dp) :: fred(3,ab_mover%natom)
 real(dp) :: xcart(3,ab_mover%natom)
 real(dp) :: xred(3,ab_mover%natom)
 real(dp) :: vel(3,ab_mover%natom)
 real(dp) :: strten(6)

!***************************************************************************
!Beginning of executable session
!***************************************************************************

 if(iexit/=0)then
   if (allocated(displacement))       then
     ABI_DEALLOCATE(displacement)
   end if
   return
 end if

!write(std_out,*) 'monte carlo 01'
!##########################################################
!### 01. Debugging and Verbose

 if(zDEBUG)then
   write(std_out,'(a,3a,40a,37a)') ch10,('-',kk=1,3),&
&   'Debugging and Verbose for monte_carlo_step',('-',kk=1,37)
   write(std_out,*) 'ionmov: ',12
   write(std_out,*) 'itime:  ',itime
   write(std_out,*) 'ntime:  ',ntime
 end if

!write(std_out,*) 'monte carlo 02'
!##########################################################
!### 02. Allocate the vectors vin, vout and hessian matrix
!###     These arrays could be allocated from a previus
!###     dataset that exit before itime==ntime

 if(itime==1)then
   if (allocated(displacement))       then
     ABI_DEALLOCATE(displacement)
   end if
 end if
 if (.not.allocated(displacement))       then
   ABI_ALLOCATE(displacement,(3,ab_mover%natom))
 end if

!write(std_out,*) 'monte carlo 03'
!##########################################################
!### 03. Obtain the present values from the history

 call hist2var(acell,hist,ab_mover%natom,rprim,rprimd,xcart,xred,zDEBUG)

 fcart(:,:) =hist%histXF(:,:,3,hist%ihist)
 fred(:,:)  =hist%histXF(:,:,4,hist%ihist)
 vel(:,:)   =hist%histV(:,:,hist%ihist)
 strten(:)  =hist%histS(:,hist%ihist)
 etotal     =hist%histE(hist%ihist)

 if(zDEBUG)then
   write (std_out,*) 'fcart:'
   do kk=1,ab_mover%natom
     write (std_out,*) fcart(:,kk)
   end do
   write (std_out,*) 'fred:'
   do kk=1,ab_mover%natom
     write (std_out,*) fred(:,kk)
   end do
   write (std_out,*) 'vel:'
   do kk=1,ab_mover%natom
     write (std_out,*) vel(:,kk)
   end do
   write (std_out,*) 'strten:'
   write (std_out,*) strten(1:3),ch10,strten(4:6)
   write (std_out,*) 'etotal:'
   write (std_out,*) etotal
 end if

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!write(std_out,*) 'monte carlo 04'
!##########################################################
!### 04. Compute monte carlo

 if(itime>1) then

   if (zDEBUG)then
     write(std_out,*) 'Computation of the monte carlo'
   end if

   do ia=1,ab_mover%natom
     do mu=1,3
       call random_number(delta)
       call random_number(acc)
       delta = (2*delta-1)/1000
       call effective_potential_getDeltaEnergy(eff_pot,de,ia,mu,ab_mover%natom,rprimd,displacement)
       displacement(mu,ia) = displacement(mu,ia) + delta 
     end do
   end do

   xred =xred + displacement


 end if ! (if itime>1)

!write(std_out,*) 'monte carlo 05'
!##########################################################
!### 06. Update the history with the prediction

!Increase indexes
 hist%ihist=hist%ihist+1

!Compute xcart from xred, and rprimd
 call xred2xcart(ab_mover%natom,rprimd,xcart,xred)

!Fill the history with the variables
!xcart, xred, acell, rprimd
 call var2hist(acell,hist,ab_mover%natom,rprim,rprimd,xcart,xred,zDEBUG)

 if(zDEBUG)then
   write (std_out,*) 'fcart:'
   do kk=1,ab_mover%natom
     write (std_out,*) fcart(:,kk)
   end do
   write (std_out,*) 'fred:'
   do kk=1,ab_mover%natom
     write (std_out,*) fred(:,kk)
   end do
   write (std_out,*) 'vel:'
   do kk=1,ab_mover%natom
     write (std_out,*) vel(:,kk)
   end do
   write (std_out,*) 'strten:'
   write (std_out,*) strten(1:3),ch10,strten(4:6)
   write (std_out,*) 'etotal:'
   write (std_out,*) etotal
 end if

 hist%histV(:,:,hist%ihist)=vel(:,:)

end subroutine monte_carlo_step
!!***


end module m_monte_carlo
