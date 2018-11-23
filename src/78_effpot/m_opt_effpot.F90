!!****m* ABINIT/m_opt_effpot
!!
!! NAME
!! m_opt_effpot
!!
!! FUNCTION
!!      
!! 
!!
!! COPYRIGHT
!! Copyright (C) 2010-2018 ABINIT group (AM)
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

module m_opt_effpot

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_abicore
 use m_xmpi
 use m_effective_potential 
 use m_effective_potential_file, only : effective_potential_file_mapHistToRef
 use m_fit_data
 use m_fit_polynomial_coeff
 use m_polynomial_coeff

 implicit none

 public :: opt_effpot 

!!****
CONTAINS 
      
!!****f* m_opt_effpot/opt_effpot 
!!
!! NAME
!! opt_effpot
!!
!! FUNCTION
!! Optimize Effective Potential by fitting the value of certain
!! coefficients while keeping the values of the others     
!!
!! INPUTS
!! eff_pot<type(effective_potential)> = effective potential
!!
!! opt_coeff(opt_ncoeff) = list of terms whose coefficients are to be
!! optimized
!!
!! hist<type(abihist)> = Training set Data(or snapshot of DFT)
!! comm = MPI communicator
!!
!! OUTPUT
!! eff_pot<type(effective_potential)> = effective potential datatype with new fitted coefficients
!!
!! PARENTS
!! multibinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine opt_effpot(eff_pot,opt_ncoeff,opt_coeff,hist,comm) 

 implicit none  

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm,opt_ncoeff
 type(effective_potential_type),target,intent(inout) :: eff_pot
 type(abihist),intent(inout) :: hist
!arrays 
 integer,intent(in) :: opt_coeff(opt_ncoeff)
!Logicals 
!Strings 
!Local variables ------------------------------
!scalars
 integer :: ii, info,natom_sc,ntime
 real(dp) :: factor,mse,msef,mses
 real(dp),parameter :: HaBohr_meVAng = 27.21138386 / 0.529177249
!arrays 
 integer :: sc_size(3)
 integer :: test(2)
 type(fit_data_type) :: fit_data
 type(polynomial_coeff_type),allocatable :: my_coeffs(:)
 real(dp) :: coeff_values(opt_ncoeff)
 real(dp), allocatable :: energy_coeffs(:,:),fcart_coeffs(:,:,:,:)
 real(dp), allocatable :: strten_coeffs(:,:,:)
!Logicals 
!Strings 
 character(len=1000) :: message
! *************************************************************************

 !Setting/Initializing Variables
  ntime = hist%mxhist
  natom_sc = size(hist%xred,2)
  factor   = 1._dp/natom_sc
  ABI_ALLOCATE(my_coeffs,(eff_pot%anharmonics_terms_type%ncoeff))

 !if the number of atoms in reference supercell into effpot is not correct,
 !wrt to the number of atom in the hist, we set map the hist and set the good supercell
  if (natom_sc /= eff_pot%supercell%natom) then
    call effective_potential_file_mapHistToRef(eff_pot,hist,comm,verbose=.TRUE.)
  end if

 !we get the size of the supercell in the hist file
  do ii=1,3
    sc_size(ii) = int(anint(sqrt(eff_pot%supercell%rprimd(ii,1)**2+&
 &                               eff_pot%supercell%rprimd(ii,2)**2+&
 &                               eff_pot%supercell%rprimd(ii,3)**2) / &
 &                          sqrt(eff_pot%crystal%rprimd(ii,1)**2+&
 &                               eff_pot%crystal%rprimd(ii,2)**2+&
 &                               eff_pot%crystal%rprimd(ii,3)**2)))
  end do


 !Before the fit, compute constants with fit_data_compute.
 !Conpute the strain of each configuration.
 !Compute the displacmeent of each configuration.
 !Compute the variation of the displacement due to strain of each configuration.
 !Compute fixed forces and stresse and get the standard deviation.
 !Compute Sheppard and al Factors  \Omega^{2} see J.Chem Phys 136, 074103 (2012) [[cite:Sheppard2012]].
  call fit_data_compute(fit_data,eff_pot,hist,comm,verbose=.TRUE.)
 
 
 !Before deleting coefficients calculate MSD of initial model  
  call fit_polynomial_coeff_computeMSD(eff_pot,hist,mse,msef,mses,&
 &                                     natom_sc,ntime,fit_data%training_set%sqomega,&
 &                                      compute_anharmonic=.TRUE.)


 !  Print the standard deviation after the fit
      write(message,'(6a,ES24.16,6a,ES24.16,2a,ES24.16,2a,ES24.16,a)' )ch10,&
 &                    ' Mean Standard Deviation values of the effective-potential',ch10,&
 &                    ' with respect to the training-set before optimization (meV/atm):',&
 &               ch10,'   Energy          : ',&
 &               mse*Ha_EV*1000*factor ,ch10,&
 &                    ' Goal function values of the effective.potential',ch10,& 
 &                    ' with respect to the test-set (eV^2/A^2):',ch10,&
 &                    '   Forces+Stresses : ',&
 &               (msef+mses)*(HaBohr_meVAng)**2,ch10,&
 &                    '   Forces          : ',&
 &               msef*(HaBohr_meVAng)**2,ch10,&
 &                    '   Stresses        : ',&
 &               mses*(HaBohr_meVAng)**2,ch10
      call wrtout(ab_out,message,'COLL')
      call wrtout(std_out,message,'COLL')


 ! Delete coeffs from eff_pot to calculate fit_data
 do ii=1,opt_ncoeff
   call polynomial_coeff_free(eff_pot%anharmonics_terms%coefficients(opt_coeff(ii)))
 end do  

 !Before the fit, compute constants with fit_data_compute.
 !Conpute the strain of each configuration.
 !Compute the displacmeent of each configuration.
 !Compute the variation of the displacement due to strain of each configuration.
 !Compute fixed forces and stresse and get the standard deviation.
 !Compute Sheppard and al Factors  \Omega^{2} see J.Chem Phys 136, 074103 (2012) [[cite:Sheppard2012]].
  call fit_data_compute(fit_data,eff_pot,hist,comm,verbose=.TRUE.)

 !After deleting coefficients calculate MSD  
  call fit_polynomial_coeff_computeMSD(eff_pot,hist,mse,msef,mses,&
 &                                     natom_sc,ntime,fit_data%training_set%sqomega,&
 &                                      compute_anharmonic=.TRUE.)


 !  Print the standard deviation after deleting
      write(message,'(6a,ES24.16,6a,ES24.16,2a,ES24.16,2a,ES24.16,a)' )ch10,&
 &                    ' Mean Standard Deviation values of the effective-potential',ch10,&
 &                    ' with respect to the training-set after deleting selected terms (meV/atm):',&
 &               ch10,'   Energy          : ',&
 &               mse*Ha_EV*1000*factor ,ch10,&
 &                    ' Goal function values of the effective.potential',ch10,& 
 &                    ' with respect to the test-set (eV^2/A^2):',ch10,&
 &                    '   Forces+Stresses : ',&
 &               (msef+mses)*(HaBohr_meVAng)**2,ch10,&
 &                    '   Forces          : ',&
 &               msef*(HaBohr_meVAng)**2,ch10,&
 &                    '   Stresses        : ',&
 &               mses*(HaBohr_meVAng)**2,ch10
      call wrtout(ab_out,message,'COLL')
      call wrtout(std_out,message,'COLL')


      
   test = (1,2)
   ABI_ALLOCATE(energy_coeffs,(opt_ncoeff,ntime))
   ABI_ALLOCATE(fcart_coeffs,(3,natom_sc,opt_ncoeff,ntime))
   ABI_ALLOCATE(strten_coeffs,(6,ntime,opt_ncoeff))
   call fit_polynomial_coeff_getFS(my_coeffs,fit_data%training_set%du_delta,&
&                                 fit_data%training_set%displacement,&
&                                 energy_coeffs,fcart_coeffs,natom_sc,eff_pot%crystal%natom,&
&                                 opt_ncoeff,ntime,sc_size,fit_data%training_set%strain,&
&                                 strten_coeffs,fit_data%training_set%ucvol,test,opt_ncoeff)


!  call the fit process routine
!  This routine solves the linear system proposed 
!  by C.Escorihuela-Sayalero see PRB95,094115(2017) [[cite:Escorihuela-Sayalero2017]]
   call fit_polynomial_coeff_solve(coeff_values(1:opt_ncoeff),fcart_coeffs,fit_data%fcart_diff,&
&                                  energy_coeffs,fit_data%energy_diff,info,&
&                                  test,natom_sc,opt_ncoeff,opt_ncoeff,ntime,&
&                                  strten_coeffs,fit_data%strten_diff,&
&                                  fit_data%training_set%sqomega)

 
   
  !Transfer new fitted values to coefficients 
  do ii=1,opt_ncoeff 
     my_coeffs(ii)%coefficient = coeff_values(ii) 
  end do 

  !Transfer coefficients back to effective potential 

  call effective_potential_setCoeffs(my_coeffs,eff_pot,opt_ncoeff)

 !Deallocation of fitting variables
  do ii=1,opt_ncoeff 
     call polynomial_coeff_free(my_coeffs(ii))
  end do 
 ABI_DEALLOCATE(energy_coeffs) 
 ABI_DEALLOCATE(fcart_coeffs)
 ABI_DEALLOCATE(strten_coeffs)

 !Recalculate MSD of Final Model TODO apply printing of fit_diff files 

 !Before the fit, compute constants with fit_data_compute.
 !Conpute the strain of each configuration.
 !Compute the displacmeent of each configuration.
 !Compute the variation of the displacement due to strain of each configuration.
 !Compute fixed forces and stresse and get the standard deviation.
 !Compute Sheppard and al Factors  \Omega^{2} see J.Chem Phys 136, 074103 (2012) [[cite:Sheppard2012]].
  call fit_data_compute(fit_data,eff_pot,hist,comm,verbose=.TRUE.)


 !After deleting coefficients calculate MSD  
  call fit_polynomial_coeff_computeMSD(eff_pot,hist,mse,msef,mses,&
 &                                     natom_sc,ntime,fit_data%training_set%sqomega,&
 &                                      compute_anharmonic=.TRUE.)


 !  Print the standard deviation after optimization
      write(message,'(6a,ES24.16,6a,ES24.16,2a,ES24.16,2a,ES24.16,a)' )ch10,&
 &                    ' Mean Standard Deviation values of the effective-potential',ch10,&
 &                    ' with respect to the training-set after optimizing selected terms (meV/atm):',&
 &               ch10,'   Energy          : ',&
 &               mse*Ha_EV*1000*factor ,ch10,&
 &                    ' Goal function values of the effective.potential',ch10,& 
 &                    ' with respect to the test-set (eV^2/A^2):',ch10,&
 &                    '   Forces+Stresses : ',&
 &               (msef+mses)*(HaBohr_meVAng)**2,ch10,&
 &                    '   Forces          : ',&
 &               msef*(HaBohr_meVAng)**2,ch10,&
 &                    '   Stresses        : ',&
 &               mses*(HaBohr_meVAng)**2,ch10
      call wrtout(ab_out,message,'COLL')
      call wrtout(std_out,message,'COLL')

 call fit_data_free(fit_data)
end subroutine opt_effpot


end module m_opt_effpot

      
